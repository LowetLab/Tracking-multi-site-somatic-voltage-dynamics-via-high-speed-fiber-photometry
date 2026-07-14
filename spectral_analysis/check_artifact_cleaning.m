%% ============================================================================
%  ARTIFACT CLEANING DIAGNOSTIC SCRIPT - SESSION LEVEL
%  ============================================================================
%  This script validates artifact cleaning for an entire session by:
%    1. Loading all trials in a session
%    2. Cleaning each trial BEFORE concatenation (matching pipeline logic)
%    3. Showing original LFP + motion with artifact masks (top row)
%    4. Showing cleaned LFP (bottom row)
%
%  USAGE:
%    1. Modify the QUICK CONFIG section below
%    2. Run: check_artifact_cleaning
%
%  NOTE: This script uses the EXACT same logic as the main pipeline
%        (spectral_analysis.m) to ensure diagnostic accuracy.
%
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  QUICK CONFIG - MODIFY THESE
%  ============================================================================

% Animal and session to check
ANIMAL_ID = 'Animal01';
SESSION_ID = '01_01_25-R1';  % Change to any session
FIBER_INDEX = 1;  % Which fiber to load (1 = first fiber)

% Artifact cleaning parameters (should match run_spectral_pipeline.m)
% NOTE: If artifact blocks aren't fully removed, try increasing these values:
%   - Increase SMOOTH_WINDOW_SEC to merge nearby artifacts better (e.g., 0.100-0.200)
%   - Increase PRE_PAD_SEC/POST_PAD_SEC to remove more surrounding data (e.g., 0.150-0.200)
PRE_PAD_SEC = 0.100;      % Pre-padding around artifacts (seconds)
POST_PAD_SEC = 0.100;     % Post-padding around artifacts (seconds)
SMOOTH_WINDOW_SEC = 0.050; % Smoothing window to merge nearby artifacts (seconds)

%% ============================================================================
%  LOAD ANIMAL DATABASE
%  ============================================================================

addpath(fullfile(fileparts(mfilename('fullpath')), 'config'));
addpath(fullfile(fileparts(mfilename('fullpath')), 'core'));  % For helper functions if needed
animals = animal_session_database();

% Find the animal
animal_idx = find(strcmp({animals.mouse_id}, ANIMAL_ID));
if isempty(animal_idx)
    error('Animal "%s" not found in database', ANIMAL_ID);
end
animal = animals(animal_idx);

% Find the session
session_idx = find(strcmp({animal.sessions.session_id}, SESSION_ID));
if isempty(session_idx)
    error('Session "%s" not found for animal "%s"', SESSION_ID, ANIMAL_ID);
end
session = animal.sessions(session_idx);

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  ARTIFACT CLEANING DIAGNOSTIC - SESSION LEVEL\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  Animal: %s\n', ANIMAL_ID);
fprintf('  Session: %s\n', SESSION_ID);
fprintf('  Number of trials: %d\n', session.num_trials);
fprintf('════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
%  LOAD AND CLEAN EACH TRIAL (MATCHING PIPELINE LOGIC)
%  ============================================================================
%  CRITICAL: Clean each trial BEFORE concatenation, not after!
%  This matches the exact logic in spectral_analysis.m (lines 233-310)

fprintf('Loading and cleaning trials (matching pipeline logic)...\n');

% Initialize concatenated data (for visualization)
concat_lfp_original = [];
concat_lfp_cleaned = [];
concat_fiber_original = [];
concat_fiber_cleaned = [];
concat_motion_original = [];
concat_motion_cleaned = [];
concat_time_original = [];
concat_time_cleaned = [];
concat_artifact_mask_original = [];  % For visualization only
concat_artifact_mask_processed = [];  % For visualization only
trial_boundaries = [0];  % Sample indices where trials start
fs = 0;
fiber_index = 1;  % Default to fiber 1 (can be configured)

% Store trial info
trial_info = struct([]);

% Create config struct for cleaning
cfg_cleaning = struct();
cfg_cleaning.artifact = struct();
cfg_cleaning.artifact.pre_pad_sec = PRE_PAD_SEC;
cfg_cleaning.artifact.post_pad_sec = POST_PAD_SEC;
cfg_cleaning.artifact.smooth_window_sec = SMOOTH_WINDOW_SEC;

for loop_idx = 1:session.num_trials
    trial_path = session.trial_paths{loop_idx};
    
    if ~exist(trial_path, 'file')
        warning('Trial %d file not found: %s', loop_idx, trial_path);
        continue;
    end
    
    fprintf('  Trial %d/%d: ', loop_idx, session.num_trials);
    
    try
        % =====================================================================
        % STEP 1: Extract actual trial number from filename (CRITICAL!)
        % =====================================================================
        % The loop index might not match the trial number in artifact info
        % (especially for combined sessions). Extract from filename.
        actual_trial_num = extract_trial_number_from_path(trial_path);
        if isempty(actual_trial_num)
            actual_trial_num = loop_idx;  % Fallback to loop index
            warning('Could not extract trial number from path: %s', trial_path);
        end
        
        fprintf('Trial%d (from filename)...', actual_trial_num);
        
        % =====================================================================
        % STEP 2: Load trial data (matching load_trial_data logic)
        % =====================================================================
        trial_data = load_trial_data_for_diagnostic(trial_path, FIBER_INDEX);
        
        if isempty(trial_data)
            fprintf('LOAD FAILED\n');
            continue;
        end
        
        if fs == 0
            fs = trial_data.fs;
        end
        
        % =====================================================================
        % STEP 3: Apply artifact cleaning to THIS TRIAL (before concatenation)
        % =====================================================================
        [cleaned_trial_data, cleaning_info] = apply_artifact_cleaning_diagnostic(...
            trial_data, trial_path, actual_trial_num, cfg_cleaning);
        
        % =====================================================================
        % STEP 4: Store original data for visualization
        % =====================================================================
        % Calculate time offset for concatenation
        if isempty(concat_time_original)
            time_offset = 0;
        else
            time_offset = concat_time_original(end) + 1/fs;
        end
        
        % Store trial info
        trial_info(loop_idx).loop_idx = loop_idx;
        trial_info(loop_idx).actual_trial_num = actual_trial_num;
        trial_info(loop_idx).concat_start_idx_original = length(concat_lfp_original) + 1;
        trial_info(loop_idx).concat_end_idx_original = length(concat_lfp_original) + length(trial_data.lfp);
        trial_info(loop_idx).concat_start_idx_cleaned = length(concat_lfp_cleaned) + 1;
        trial_info(loop_idx).concat_end_idx_cleaned = length(concat_lfp_cleaned) + length(cleaned_trial_data.lfp);
        trial_info(loop_idx).duration_sec = trial_data.duration;
        trial_info(loop_idx).artifact_pct_original = cleaning_info.removed_pct_original;
        trial_info(loop_idx).artifact_pct_processed = cleaning_info.removed_pct;
        trial_info(loop_idx).cleaning_applied = cleaning_info.applied;
        
        % Concatenate ORIGINAL data (for visualization)
        concat_lfp_original = [concat_lfp_original; trial_data.lfp(:)]; %#ok<AGROW>
        concat_fiber_original = [concat_fiber_original; trial_data.gevi(:)]; %#ok<AGROW>
        concat_motion_original = [concat_motion_original; trial_data.speed(:)]; %#ok<AGROW>
        concat_time_original = [concat_time_original; trial_data.time(:) + time_offset]; %#ok<AGROW>
        
        % Concatenate artifact masks (for visualization)
        if ~isempty(cleaning_info.artifact_mask_original)
            concat_artifact_mask_original = [concat_artifact_mask_original; cleaning_info.artifact_mask_original(:)]; %#ok<AGROW>
        else
            % No mask found - pad with false
            concat_artifact_mask_original = [concat_artifact_mask_original; false(length(trial_data.lfp), 1)]; %#ok<AGROW>
        end
        
        if ~isempty(cleaning_info.artifact_mask_processed)
            concat_artifact_mask_processed = [concat_artifact_mask_processed; cleaning_info.artifact_mask_processed(:)]; %#ok<AGROW>
        else
            % No mask found - pad with false
            concat_artifact_mask_processed = [concat_artifact_mask_processed; false(length(trial_data.lfp), 1)]; %#ok<AGROW>
        end
        
        % Concatenate CLEANED data
        if cleaning_info.applied && ~isempty(cleaned_trial_data.lfp)
            % Calculate time offset for cleaned data (may have gaps)
            if isempty(concat_time_cleaned)
                cleaned_time_offset = 0;
            else
                cleaned_time_offset = concat_time_cleaned(end) + 1/fs;
            end
            concat_lfp_cleaned = [concat_lfp_cleaned; cleaned_trial_data.lfp(:)]; %#ok<AGROW>
            concat_fiber_cleaned = [concat_fiber_cleaned; cleaned_trial_data.gevi(:)]; %#ok<AGROW>
            concat_motion_cleaned = [concat_motion_cleaned; cleaned_trial_data.speed(:)]; %#ok<AGROW>
            concat_time_cleaned = [concat_time_cleaned; cleaned_trial_data.time(:) + cleaned_time_offset]; %#ok<AGROW>
        else
            % No cleaning applied or no cleaned data - use original
            if isempty(concat_time_cleaned)
                cleaned_time_offset = 0;
            else
                cleaned_time_offset = concat_time_cleaned(end) + 1/fs;
            end
            concat_lfp_cleaned = [concat_lfp_cleaned; trial_data.lfp(:)]; %#ok<AGROW>
            concat_fiber_cleaned = [concat_fiber_cleaned; trial_data.gevi(:)]; %#ok<AGROW>
            concat_motion_cleaned = [concat_motion_cleaned; trial_data.speed(:)]; %#ok<AGROW>
            concat_time_cleaned = [concat_time_cleaned; trial_data.time(:) + cleaned_time_offset]; %#ok<AGROW>
        end
        
        trial_boundaries = [trial_boundaries; length(concat_lfp_original)]; %#ok<AGROW>
        
        % Print status
        if cleaning_info.applied
            fprintf('OK [cleaned: %.1f%% removed (orig: %.1f%%)]\n', ...
                cleaning_info.removed_pct, cleaning_info.removed_pct_original);
        else
            fprintf('OK [no artifact info found]\n');
        end
        
    catch ME
        warning('Failed to process trial %d: %s', loop_idx, ME.message);
        continue;
    end
end

if isempty(concat_lfp_original)
    error('No valid trials loaded');
end

fprintf('\n  Session total: %.1f sec (%d samples original, %d samples cleaned)\n', ...
    concat_time_original(end) - concat_time_original(1), length(concat_lfp_original), length(concat_lfp_cleaned));
fprintf('  Total artifacts (original): %.1f%%\n', ...
    100 * sum(concat_artifact_mask_original) / length(concat_artifact_mask_original));

% ============================================================================
% VERIFY MASK ALIGNMENT AND CLEANING EFFECTIVENESS
% ============================================================================
fprintf('\n  VERIFICATION:\n');
fprintf('    Original mask length: %d\n', length(concat_artifact_mask_original));
fprintf('    Processed mask length: %d\n', length(concat_artifact_mask_processed));
fprintf('    Original data length: %d\n', length(concat_lfp_original));
fprintf('    Cleaned data length: %d\n', length(concat_lfp_cleaned));

% Check if mask lengths match
if length(concat_artifact_mask_processed) ~= length(concat_lfp_original)
    warning('MASK LENGTH MISMATCH: Processed mask (%d) != Original data (%d)', ...
        length(concat_artifact_mask_processed), length(concat_lfp_original));
    % Try to fix by truncating or padding
    if length(concat_artifact_mask_processed) > length(concat_lfp_original)
        concat_artifact_mask_processed = concat_artifact_mask_processed(1:length(concat_lfp_original));
        fprintf('    Fixed: Truncated processed mask to match data length\n');
    else
        concat_artifact_mask_processed = [concat_artifact_mask_processed; false(length(concat_lfp_original) - length(concat_artifact_mask_processed), 1)];
        fprintf('    Fixed: Padded processed mask to match data length\n');
    end
end

% Verify that what should be removed actually was removed
expected_removed = sum(concat_artifact_mask_processed);
actual_removed_lfp = length(concat_lfp_original) - length(concat_lfp_cleaned);
actual_removed_fiber = length(concat_fiber_original) - length(concat_fiber_cleaned);
actual_removed_motion = length(concat_motion_original) - length(concat_motion_cleaned);
fprintf('    Expected removed (from processed mask): %d samples (%.1f%%)\n', ...
    expected_removed, 100 * expected_removed / length(concat_lfp_original));
fprintf('    Actual removed - LFP: %d samples (%.1f%%)\n', ...
    actual_removed_lfp, 100 * actual_removed_lfp / length(concat_lfp_original));
fprintf('    Actual removed - Fiber: %d samples (%.1f%%)\n', ...
    actual_removed_fiber, 100 * actual_removed_fiber / length(concat_fiber_original));
fprintf('    Actual removed - Motion: %d samples (%.1f%%)\n', ...
    actual_removed_motion, 100 * actual_removed_motion / length(concat_motion_original));

% Verify all three signals were cleaned consistently
if abs(actual_removed_lfp - actual_removed_fiber) > 1 || abs(actual_removed_lfp - actual_removed_motion) > 1
    warning('INCONSISTENT CLEANING: LFP removed %d, Fiber removed %d, Motion removed %d', ...
        actual_removed_lfp, actual_removed_fiber, actual_removed_motion);
else
    fprintf('    ✓ All three signals (LFP, Fiber, Motion) cleaned consistently\n');
end

% Verify cleaned data lengths match
if length(concat_lfp_cleaned) ~= length(concat_fiber_cleaned) || ...
   length(concat_lfp_cleaned) ~= length(concat_motion_cleaned) || ...
   length(concat_lfp_cleaned) ~= length(concat_time_cleaned)
    warning('CLEANED DATA LENGTH MISMATCH: LFP=%d, Fiber=%d, Motion=%d, Time=%d', ...
        length(concat_lfp_cleaned), length(concat_fiber_cleaned), ...
        length(concat_motion_cleaned), length(concat_time_cleaned));
else
    fprintf('    ✓ All cleaned signals have matching lengths (%d samples)\n', length(concat_lfp_cleaned));
end

% Show smoothing/padding expansion
original_artifacts = sum(concat_artifact_mask_original);
processed_artifacts = sum(concat_artifact_mask_processed);
expansion_pct = 100 * (processed_artifacts - original_artifacts) / max(original_artifacts, 1);
fprintf('    Original artifacts: %d samples (%.1f%%)\n', ...
    original_artifacts, 100 * original_artifacts / length(concat_lfp_original));
fprintf('    Processed artifacts (after smoothing/padding): %d samples (%.1f%%)\n', ...
    processed_artifacts, 100 * processed_artifacts / length(concat_lfp_original));
fprintf('    Expansion from smoothing/padding: +%d samples (%.1f%% increase)\n', ...
    processed_artifacts - original_artifacts, expansion_pct);

if abs(expected_removed - actual_removed_lfp) > 10  % Allow small tolerance for rounding
    warning('CLEANING MISMATCH: Expected to remove %d samples, but actually removed %d (diff: %d)', ...
        expected_removed, actual_removed_lfp, abs(expected_removed - actual_removed_lfp));
    fprintf('    This suggests the cleaning may not be working correctly!\n');
else
    fprintf('    ✓ Cleaning verification passed (difference: %d samples)\n', abs(expected_removed - actual_removed_lfp));
end

% Check if smoothing/padding is sufficient
% If original artifacts are much smaller than processed, smoothing/padding is working
% But if they're similar, maybe parameters need to be increased
if original_artifacts > 0 && expansion_pct < 20
    fprintf('    ⚠ WARNING: Smoothing/padding only expanded artifacts by %.1f%%\n', expansion_pct);
    fprintf('      This might be insufficient. Consider increasing:\n');
    fprintf('        - SMOOTH_WINDOW_SEC (current: %.3f s)\n', SMOOTH_WINDOW_SEC);
    fprintf('        - PRE_PAD_SEC (current: %.3f s)\n', PRE_PAD_SEC);
    fprintf('        - POST_PAD_SEC (current: %.3f s)\n', POST_PAD_SEC);
end

%% ============================================================================
%  CREATE VISUALIZATION
%  ============================================================================

fprintf('\nGenerating diagnostic plot...\n');

% Publication-quality colors (matching artifact_removal_lfp_multisession.m)
color_lfp = [0.35, 0.25, 0.45];        % Dark purple for LFP
color_motion = [0.993, 0.7, 0.4];       % Warm orange for motion
artifact_color = [0.9, 0.15, 0.15];    % Bright red for artifact shading
artifact_alpha = 0.35;                 % Semi-transparent
trial_boundary_color = [0.2, 0.2, 0.2]; % Dark gray for trial boundaries
trial_boundary_style = '--';           % Dashed line
trial_boundary_width = 1.5;

fig = figure('Name', sprintf('Artifact Cleaning Diagnostic: %s %s', ...
    ANIMAL_ID, SESSION_ID), ...
    'Position', [100, 100, 1900, 1000], ...
    'Color', 'w');

% Layout: Two rows
% Top row: LFP (70%) + Motion (30%) with artifacts
% Bottom row: Cleaned LFP (full height)

margin_left = 0.10;
margin_right = 0.03;
margin_bottom = 0.08;
margin_top = 0.06;
row_spacing = 0.02;

% Row heights
top_row_height = 0.45;  % 45% for top row (LFP + Fiber + Motion)
bottom_row_height = 0.45;  % 45% for bottom row (Cleaned LFP + Fiber + Motion)

% Within each row: LFP 40%, Fiber 40%, Motion 20%
lfp_ratio = 0.40;
fiber_ratio = 0.40;
motion_ratio = 0.20;
subrow_spacing = 0.005;

lfp_height = top_row_height * lfp_ratio - subrow_spacing/2;
fiber_height = top_row_height * fiber_ratio - subrow_spacing/2;
motion_height = top_row_height * motion_ratio - subrow_spacing/2;

ax_width = 1 - margin_left - margin_right;

% ============================================================================
% TOP ROW: Original LFP + Fiber + Motion with Artifacts
% ============================================================================

% LFP subplot (top of top row)
lfp_bottom = margin_bottom + bottom_row_height + row_spacing + fiber_height + motion_height + 2*subrow_spacing;
ax_lfp = axes('Parent', fig, 'Position', [margin_left, lfp_bottom, ax_width, lfp_height]);
hold(ax_lfp, 'on');

% Z-score LFP for consistent scaling
lfp_z_original = zscore(concat_lfp_original);
lfp_ylim = [-5, 5];

% Plot artifact regions (behind trace)
% First, plot ORIGINAL mask (lighter, to show what was initially detected)
segments_original = find_artifact_segments(concat_artifact_mask_original);
for s = 1:size(segments_original, 1)
    seg_start = segments_original(s, 1);
    seg_end = segments_original(s, 2);
    if seg_start <= length(concat_time_original) && seg_end <= length(concat_time_original)
        fill(ax_lfp, [concat_time_original(seg_start), concat_time_original(seg_end), concat_time_original(seg_end), concat_time_original(seg_start)], ...
            [lfp_ylim(1), lfp_ylim(1), lfp_ylim(2), lfp_ylim(2)], ...
            [0.9, 0.7, 0.7], 'FaceAlpha', artifact_alpha * 0.5, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Then, plot PROCESSED mask (darker, to show what will actually be removed after smoothing/padding)
segments = find_artifact_segments(concat_artifact_mask_processed);
for s = 1:size(segments, 1)
    seg_start = segments(s, 1);
    seg_end = segments(s, 2);
    if seg_start <= length(concat_time_original) && seg_end <= length(concat_time_original)
        fill(ax_lfp, [concat_time_original(seg_start), concat_time_original(seg_end), concat_time_original(seg_end), concat_time_original(seg_start)], ...
            [lfp_ylim(1), lfp_ylim(1), lfp_ylim(2), lfp_ylim(2)], ...
            artifact_color, 'FaceAlpha', artifact_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Plot trial boundaries
for b = 2:length(trial_boundaries)-1
    if trial_boundaries(b) <= length(concat_time_original)
        boundary_time = concat_time_original(trial_boundaries(b));
        line(ax_lfp, [boundary_time, boundary_time], lfp_ylim, ...
            'Color', trial_boundary_color, 'LineStyle', trial_boundary_style, ...
            'LineWidth', trial_boundary_width, 'HandleVisibility', 'off');
    end
end

% Plot original LFP trace
plot(ax_lfp, concat_time_original, lfp_z_original, 'Color', color_lfp, 'LineWidth', 0.5);

% Configure LFP axis
xlim(ax_lfp, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_lfp, lfp_ylim);
ax_lfp.YColor = color_lfp;
ax_lfp.YTick = [-4, 0, 4];
ax_lfp.YTickLabel = {'-4σ', '0', '4σ'};
ax_lfp.FontSize = 9;
ax_lfp.XTickLabel = [];  % No x-labels on LFP row
ax_lfp.Box = 'on';
ax_lfp.XColor = [0.7, 0.7, 0.7];

% Session label
pct_original_artifacts = 100 * sum(concat_artifact_mask_original) / length(concat_artifact_mask_original);
ylabel(ax_lfp, sprintf('%s\n%.1f%% art.', ...
    strrep(SESSION_ID, '_', '-'), pct_original_artifacts), ...
    'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k', 'Rotation', 0, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Add trial labels at top
for t = 1:length(trial_info)
    if isfield(trial_info(t), 'concat_start_idx_original') && trial_info(t).concat_start_idx_original > 0
        start_idx = trial_info(t).concat_start_idx_original;
        end_idx = trial_info(t).concat_end_idx_original;
        if start_idx <= length(concat_time_original) && end_idx <= length(concat_time_original)
            trial_start_time = concat_time_original(start_idx);
            trial_end_time = concat_time_original(min(end_idx, length(concat_time_original)));
            mid_time = (trial_start_time + trial_end_time) / 2;
            
            label_str = sprintf('T%d: %.0f%%', trial_info(t).actual_trial_num, trial_info(t).artifact_pct_original);
            text(ax_lfp, mid_time, lfp_ylim(2) * 0.92, label_str, ...
                'Color', [0.15, 0.5, 0.15], 'FontSize', 8, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', [0.92, 1, 0.92], 'EdgeColor', [0.15, 0.5, 0.15], ...
                'Margin', 1, 'Clipping', 'on');
        end
    end
end

hold(ax_lfp, 'off');

% Fiber subplot (middle of top row)
fiber_bottom = margin_bottom + bottom_row_height + row_spacing + motion_height + subrow_spacing;
ax_fiber = axes('Parent', fig, 'Position', [margin_left, fiber_bottom, ax_width, fiber_height]);
hold(ax_fiber, 'on');

% Z-score fiber for consistent scaling
fiber_z_original = zscore(concat_fiber_original);
fiber_ylim = [-5, 5];

% Plot artifact regions (same as LFP)
for s = 1:size(segments_original, 1)
    seg_start = segments_original(s, 1);
    seg_end = segments_original(s, 2);
    if seg_start <= length(concat_time_original) && seg_end <= length(concat_time_original)
        fill(ax_fiber, [concat_time_original(seg_start), concat_time_original(seg_end), concat_time_original(seg_end), concat_time_original(seg_start)], ...
            [fiber_ylim(1), fiber_ylim(1), fiber_ylim(2), fiber_ylim(2)], ...
            [0.9, 0.7, 0.7], 'FaceAlpha', artifact_alpha * 0.5, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

for s = 1:size(segments, 1)
    seg_start = segments(s, 1);
    seg_end = segments(s, 2);
    if seg_start <= length(concat_time_original) && seg_end <= length(concat_time_original)
        fill(ax_fiber, [concat_time_original(seg_start), concat_time_original(seg_end), concat_time_original(seg_end), concat_time_original(seg_start)], ...
            [fiber_ylim(1), fiber_ylim(1), fiber_ylim(2), fiber_ylim(2)], ...
            artifact_color, 'FaceAlpha', artifact_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Plot trial boundaries
for b = 2:length(trial_boundaries)-1
    if trial_boundaries(b) <= length(concat_time_original)
        boundary_time = concat_time_original(trial_boundaries(b));
        line(ax_fiber, [boundary_time, boundary_time], fiber_ylim, ...
            'Color', trial_boundary_color, 'LineStyle', trial_boundary_style, ...
            'LineWidth', trial_boundary_width, 'HandleVisibility', 'off');
    end
end

% Plot original fiber trace
color_fiber = [0.8, 0.2, 0.2];  % Reddish for fiber
plot(ax_fiber, concat_time_original, fiber_z_original, 'Color', color_fiber, 'LineWidth', 0.5);

% Configure fiber axis
xlim(ax_fiber, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_fiber, fiber_ylim);
ax_fiber.YColor = color_fiber;
ax_fiber.YTick = [-4, 0, 4];
ax_fiber.YTickLabel = {'-4σ', '0', '4σ'};
ax_fiber.FontSize = 9;
ax_fiber.XTickLabel = [];  % No x-labels on fiber row
ax_fiber.Box = 'on';
ax_fiber.XColor = [0.7, 0.7, 0.7];

% Fiber label
ylabel(ax_fiber, 'Fiber (z-scored)', 'FontSize', 10, 'FontWeight', 'bold', ...
    'Color', color_fiber, 'Rotation', 0, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

hold(ax_fiber, 'off');

% Motion subplot (bottom of top row)
motion_bottom = margin_bottom + bottom_row_height + row_spacing;
ax_motion = axes('Parent', fig, 'Position', [margin_left, motion_bottom, ax_width, motion_height]);
hold(ax_motion, 'on');

% Motion y-limits
motion_min = 0;
motion_max = max(prctile(concat_motion_original, 99.5), 0.5);
motion_ylim = [motion_min, motion_max * 1.15];

% Plot artifact regions on motion too (for visual alignment)
for s = 1:size(segments, 1)
    seg_start = segments(s, 1);
    seg_end = segments(s, 2);
    if seg_start <= length(concat_time_original) && seg_end <= length(concat_time_original)
        fill(ax_motion, [concat_time_original(seg_start), concat_time_original(seg_end), concat_time_original(seg_end), concat_time_original(seg_start)], ...
            [motion_ylim(1), motion_ylim(1), motion_ylim(2), motion_ylim(2)], ...
            artifact_color, 'FaceAlpha', artifact_alpha * 0.5, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
end

% Plot trial boundaries on motion
for b = 2:length(trial_boundaries)-1
    if trial_boundaries(b) <= length(concat_time_original)
        boundary_time = concat_time_original(trial_boundaries(b));
        line(ax_motion, [boundary_time, boundary_time], motion_ylim, ...
            'Color', trial_boundary_color, 'LineStyle', trial_boundary_style, ...
            'LineWidth', trial_boundary_width * 0.7, 'HandleVisibility', 'off');
    end
end

% Plot motion trace
plot(ax_motion, concat_time_original, concat_motion_original, 'Color', color_motion, 'LineWidth', 0.6);

% Configure motion axis
xlim(ax_motion, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_motion, motion_ylim);
ax_motion.YColor = color_motion;
ax_motion.FontSize = 9;
ax_motion.Box = 'on';
ax_motion.XColor = [0.5, 0.5, 0.5];
xlabel(ax_motion, 'Time (s)', 'FontSize', 10);

% Motion label
ylabel(ax_motion, 'Speed', 'FontSize', 9, 'Color', color_motion, ...
    'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

hold(ax_motion, 'off');

% ============================================================================
% BOTTOM ROW: Cleaned LFP + Fiber + Motion
% ============================================================================

% Helper function to plot cleaned signal with gap handling
function plot_cleaned_with_gaps(ax, time_vec, signal_vec, color, ylims, fs_local, gap_color)
    if isempty(signal_vec) || length(time_vec) ~= length(signal_vec)
        return;
    end
    
    if length(time_vec) > 1
        time_diffs = diff(time_vec);
        expected_interval = 1 / fs_local;
        gap_threshold = 1.5 * expected_interval;
        gap_indices = find(time_diffs > gap_threshold);
        
        if isempty(gap_indices)
            plot(ax, time_vec, signal_vec, 'Color', color, 'LineWidth', 0.5);
        else
            segment_start = 1;
            for g = 1:length(gap_indices)
                segment_end = gap_indices(g);
                if segment_end >= segment_start
                    plot(ax, time_vec(segment_start:segment_end), ...
                        signal_vec(segment_start:segment_end), ...
                        'Color', color, 'LineWidth', 0.5);
                end
                gap_time = time_vec(segment_end);
                line(ax, [gap_time, gap_time], ylims, ...
                    'Color', gap_color, 'LineStyle', '--', 'LineWidth', 1.5);
                segment_start = gap_indices(g) + 1;
            end
            if segment_start <= length(time_vec)
                plot(ax, time_vec(segment_start:end), signal_vec(segment_start:end), ...
                    'Color', color, 'LineWidth', 0.5);
            end
        end
    end
end

% Cleaned LFP subplot (top of bottom row)
cleaned_lfp_bottom = margin_bottom + fiber_height + motion_height + 2*subrow_spacing;
ax_cleaned_lfp = axes('Parent', fig, 'Position', [margin_left, cleaned_lfp_bottom, ax_width, lfp_height]);
hold(ax_cleaned_lfp, 'on');

if ~isempty(concat_lfp_cleaned)
    lfp_z_cleaned = zscore(concat_lfp_cleaned);
    plot_cleaned_with_gaps(ax_cleaned_lfp, concat_time_cleaned, lfp_z_cleaned, color_lfp, lfp_ylim, fs, artifact_color);
end

xlim(ax_cleaned_lfp, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_cleaned_lfp, lfp_ylim);
ax_cleaned_lfp.YColor = color_lfp;
ax_cleaned_lfp.YTick = [-4, 0, 4];
ax_cleaned_lfp.YTickLabel = {'-4σ', '0', '4σ'};
ax_cleaned_lfp.FontSize = 9;
ax_cleaned_lfp.XTickLabel = [];
ax_cleaned_lfp.Box = 'on';
ax_cleaned_lfp.XColor = [0.7, 0.7, 0.7];
ylabel(ax_cleaned_lfp, 'LFP (z-scored)', 'FontSize', 10, 'FontWeight', 'bold', ...
    'Color', color_lfp, 'Rotation', 0, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
hold(ax_cleaned_lfp, 'off');

% Cleaned Fiber subplot (middle of bottom row)
cleaned_fiber_bottom = margin_bottom + motion_height + subrow_spacing;
ax_cleaned_fiber = axes('Parent', fig, 'Position', [margin_left, cleaned_fiber_bottom, ax_width, fiber_height]);
hold(ax_cleaned_fiber, 'on');

if ~isempty(concat_fiber_cleaned)
    fiber_z_cleaned = zscore(concat_fiber_cleaned);
    plot_cleaned_with_gaps(ax_cleaned_fiber, concat_time_cleaned, fiber_z_cleaned, color_fiber, fiber_ylim, fs, artifact_color);
end

xlim(ax_cleaned_fiber, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_cleaned_fiber, fiber_ylim);
ax_cleaned_fiber.YColor = color_fiber;
ax_cleaned_fiber.YTick = [-4, 0, 4];
ax_cleaned_fiber.YTickLabel = {'-4σ', '0', '4σ'};
ax_cleaned_fiber.FontSize = 9;
ax_cleaned_fiber.XTickLabel = [];
ax_cleaned_fiber.Box = 'on';
ax_cleaned_fiber.XColor = [0.7, 0.7, 0.7];
ylabel(ax_cleaned_fiber, 'Fiber (z-scored)', 'FontSize', 10, 'FontWeight', 'bold', ...
    'Color', color_fiber, 'Rotation', 0, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
hold(ax_cleaned_fiber, 'off');

% Cleaned Motion subplot (bottom of bottom row)
ax_cleaned_motion = axes('Parent', fig, 'Position', [margin_left, margin_bottom, ax_width, motion_height]);
hold(ax_cleaned_motion, 'on');

if ~isempty(concat_motion_cleaned) && length(concat_motion_cleaned) == length(concat_time_cleaned)
    plot_cleaned_with_gaps(ax_cleaned_motion, concat_time_cleaned, concat_motion_cleaned, color_motion, motion_ylim, fs, artifact_color);
end

xlim(ax_cleaned_motion, [concat_time_original(1), concat_time_original(end)]);
ylim(ax_cleaned_motion, motion_ylim);
ax_cleaned_motion.YColor = color_motion;
ax_cleaned_motion.FontSize = 9;
ax_cleaned_motion.Box = 'on';
ax_cleaned_motion.XColor = [0.5, 0.5, 0.5];
xlabel(ax_cleaned_motion, 'Time (s)', 'FontSize', 10);
ylabel(ax_cleaned_motion, 'Speed', 'FontSize', 9, 'Color', color_motion, ...
    'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');

% Calculate removal statistics
n_original = length(concat_lfp_original);
n_cleaned = length(concat_lfp_cleaned);
n_removed = n_original - n_cleaned;
pct_removed = 100 * n_removed / n_original;
num_segments = size(find_artifact_segments(concat_artifact_mask_processed), 1);

title(ax_cleaned_motion, sprintf('CLEANED DATA (%.1f%% removed, %d segments)', ...
    pct_removed, num_segments), ...
    'FontSize', 14, 'FontWeight', 'bold');
hold(ax_cleaned_motion, 'off');

% ============================================================================
% MAIN TITLE
% ============================================================================
sgtitle(sprintf('Artifact Cleaning Diagnostic: %s %s (%.1f%% artifacts, %.1f%% removed)', ...
    ANIMAL_ID, SESSION_ID, pct_original_artifacts, pct_removed), ...
    'FontSize', 16, 'FontWeight', 'bold');

fprintf('  Plot generated.\n');

%% ============================================================================
%  SUMMARY
%  ============================================================================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  SUMMARY\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  Original data: %d samples (%.1f seconds)\n', n_original, concat_time_original(end) - concat_time_original(1));
if ~isempty(concat_time_cleaned)
    fprintf('  Cleaned data: %d samples (%.1f seconds)\n', n_cleaned, concat_time_cleaned(end) - concat_time_cleaned(1));
else
    fprintf('  Cleaned data: 0 samples\n');
end
fprintf('  Removed: %d samples (%.1f%%)\n', n_removed, pct_removed);
fprintf('  Artifact segments: %d\n', num_segments);
fprintf('\n  Validation:\n');
if pct_removed > 0
    fprintf('    ✓ Artifacts were removed\n');
    if ~isempty(concat_lfp_cleaned)
        fprintf('    ✓ Cleaned LFP has data\n');
    else
        fprintf('    ✗ WARNING: Cleaned LFP is empty!\n');
    end
else
    fprintf('    ✗ WARNING: No artifacts were removed (check artifact mask)\n');
end
fprintf('════════════════════════════════════════════════════════════════════════\n');

%% ============================================================================
%  HELPER FUNCTIONS (Matching spectral_analysis.m exactly)
%  ============================================================================

function trial_num = extract_trial_number_from_path(trial_path)
%EXTRACT_TRIAL_NUMBER_FROM_PATH Extract trial number from filename
%  EXACT COPY from spectral_analysis.m (lines 2382-2435)

trial_num = [];

% Get just the filename
[~, filename, ~] = fileparts(trial_path);

% Try to match pattern "_Trial{N}_" in filename
tokens = regexp(filename, '_Trial(\d+)_', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Try pattern "Trial{N}_" at start of filename
tokens = regexp(filename, '^Trial(\d+)_', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Try to extract from folder path
tokens = regexp(trial_path, '[/\\]Trial(\d+)_', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Fallback: try any pattern with Trial followed by a number
tokens = regexp(filename, 'Trial(\d+)', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
end

end

function data = load_trial_data_for_diagnostic(trial_path, fiber_index)
%LOAD_TRIAL_DATA_FOR_DIAGNOSTIC Load trial data (matching pipeline logic)
%  Loads LFP, fiber (GEVI), and motion traces

data = [];

% Motion conversion constants (matching spectral_analysis.m)
WHEEL_DIAMETER_CM = 19.0;
WHEEL_CIRCUMFERENCE_CM = pi * WHEEL_DIAMETER_CM;
ENCODER_COUNTS_PER_REV = 1024;
EPHYS_SAMPLING_RATE = 30000;
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV;
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM;
MOTION_SMOOTH_SAMPLES = 10;

try
    loaded = load(trial_path);
    
    if ~isfield(loaded, 'FiberPhotometryAnalysis')
        warning('FiberPhotometryAnalysis structure not found');
        return;
    end
    
    FPA = loaded.FiberPhotometryAnalysis;
    
    % Extract time vector
    if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
        t = FPA.time.time_vector_seconds(:);
    else
        warning('Time vector not found');
        return;
    end
    
    % Extract sampling rate
    if isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
        fs = FPA.parameters.sampling_rate;
    elseif isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
        fs = FPA.time.sampling_rate;
    else
        fs = 1 / median(diff(t));
    end
    
    % Extract LFP trace
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_HP')
        lfp = FPA.ephys.lfp_raw_aligned_HP(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_mPFC')
        lfp = FPA.ephys.lfp_raw_aligned_mPFC(:);
    else
        warning('LFP trace not found');
        return;
    end
    
    % Extract fiber trace (GEVI) - matching spectral_analysis.m logic
    if isfield(FPA, 'signals') && isfield(FPA.signals, 'final_processed_traces')
        fiber_all = FPA.signals.final_processed_traces;
        if isstruct(fiber_all)
            warning('final_processed_traces is a struct');
            return;
        end
        if size(fiber_all, 2) >= fiber_index
            gevi = fiber_all(:, fiber_index);
        elseif size(fiber_all, 1) >= fiber_index
            gevi = fiber_all(fiber_index, :)';
        else
            warning('Fiber index %d out of range', fiber_index);
            return;
        end
        gevi = gevi(:);
    elseif isfield(FPA, 'signals') && isfield(FPA.signals, 'deltaF_F_traces')
        fiber_all = FPA.signals.deltaF_F_traces;
        if isstruct(fiber_all)
            warning('deltaF_F_traces is a struct');
            return;
        end
        gevi = fiber_all(:, min(fiber_index, size(fiber_all, 2)));
        gevi = gevi(:);
    else
        warning('Fiber trace not found');
        return;
    end
    
    % Extract motion trace
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
        motion_raw = FPA.ephys.running_velocity_smooth(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
        motion_raw = FPA.ephys.running_velocity(:);
    else
        motion_raw = zeros(size(lfp));
    end
    
    % Align lengths
    n = min([length(t), length(lfp), length(gevi), length(motion_raw)]);
    t = t(1:n);
    lfp = lfp(1:n);
    gevi = gevi(1:n);
    motion_raw = motion_raw(1:n);
    
    % Convert motion to cm/s
    speed_cm_s = motion_raw * MOTION_TO_CM_PER_S;
    
    % Smooth motion
    if MOTION_SMOOTH_SAMPLES > 1
        kernel = ones(MOTION_SMOOTH_SAMPLES, 1) / MOTION_SMOOTH_SAMPLES;
        speed_cm_s = conv(speed_cm_s, kernel, 'same');
    end
    
    % Store
    data = struct();
    data.time = t;
    data.fs = fs;
    data.lfp = lfp;
    data.gevi = gevi;
    data.speed = speed_cm_s;
    data.duration = t(end) - t(1);
    
catch ME
    warning('check_artifact_cleaning:loadFailed', 'Failed to load trial data: %s', ME.message);
    data = [];
end

end

function [cleaned_data, cleaning_info] = apply_artifact_cleaning_diagnostic(data, trial_path, trial_idx, cfg)
%APPLY_ARTIFACT_CLEANING_DIAGNOSTIC Apply artifact cleaning (EXACT COPY from spectral_analysis.m)
%  This is an exact copy of apply_artifact_cleaning (lines 1994-2222)

cleaning_info = struct();
cleaning_info.original_samples = length(data.time);
cleaning_info.cleaned_samples = cleaning_info.original_samples;
cleaning_info.removed_samples = 0;
cleaning_info.removed_pct = 0;
cleaning_info.removed_pct_original = 0;
cleaning_info.num_segments = 0;
cleaning_info.num_segments_original = 0;
cleaning_info.num_segments_processed = 0;
cleaning_info.artifact_mask_original = [];
cleaning_info.artifact_mask_processed = [];
cleaning_info.applied = false;

cleaned_data = data;

% Load artifact mask
[artifact_mask, mask_info] = load_artifact_mask_diagnostic(trial_path, trial_idx);

if isempty(artifact_mask)
    warning('artifact_cleaning:noMask', ...
        'No artifact mask found for trial %d. Data returned unchanged.', trial_idx);
    return;
end

n_artifacts_original = sum(artifact_mask);
if n_artifacts_original == 0
    warning('artifact_cleaning:noArtifacts', ...
        'Artifact mask for trial %d contains no artifacts. Data returned unchanged.', trial_idx);
    cleaning_info.artifact_mask_original = artifact_mask;
    cleaning_info.num_segments_original = 0;
    cleaning_info.removed_pct_original = 0;
    cleaning_info.applied = false;
    return;
end

cleaning_info.artifact_mask_original = artifact_mask;
cleaning_info.num_segments_original = mask_info.num_segments;

% Ensure mask length matches data length
n_data = length(data.time);
n_mask = length(artifact_mask);

if n_mask ~= n_data
    if n_mask > n_data
        artifact_mask = artifact_mask(1:n_data);
        warning('artifact_cleaning:lengthMismatch', ...
            'Artifact mask longer than data (%d vs %d). Truncating mask.', n_mask, n_data);
    else
        artifact_mask = [artifact_mask; false(n_data - n_mask, 1)];
        warning('artifact_cleaning:lengthMismatch', ...
            'Artifact mask shorter than data (%d vs %d). Padding with clean.', n_mask, n_data);
    end
end

% Get cleaning parameters
pre_pad_sec = 0.100;
post_pad_sec = 0.100;
smooth_window_sec = 0.050;

if isfield(cfg, 'artifact')
    if isfield(cfg.artifact, 'pre_pad_sec')
        pre_pad_sec = cfg.artifact.pre_pad_sec;
    end
    if isfield(cfg.artifact, 'post_pad_sec')
        post_pad_sec = cfg.artifact.post_pad_sec;
    end
    if isfield(cfg.artifact, 'smooth_window_sec')
        smooth_window_sec = cfg.artifact.smooth_window_sec;
    end
end

fs = data.fs;

% ============================================================================
% INTELLIGENT ARTIFACT MASK PROCESSING
% ============================================================================
% Apply intelligent processing to the artifact mask:
%   1. Optional smoothing to merge nearby artifacts
%   2. Ensure contiguous blocks are merged
%   3. Add padding around artifacts to remove potentially contaminated data
% ============================================================================

% Step 1: Optional smoothing to merge nearby artifacts
% This helps merge artifacts that are very close together (within smooth_window)
% Uses morphological dilation: if any sample in a window around a point is artifact,
% mark that point as artifact. This effectively extends each artifact by half the
% window size on each side, merging nearby artifacts.
if smooth_window_sec > 0
    smooth_samples = max(1, round(smooth_window_sec * fs));
    if smooth_samples > 1
        % Efficient dilation using convolution with a ones kernel
        % This is equivalent to: for each point, check if any sample in window is artifact
        half_win = floor(smooth_samples / 2);
        artifact_mask_dilated = false(size(artifact_mask));
        
        % Use a more efficient approach: find all artifact positions and dilate them
        artifact_indices = find(artifact_mask);
        for idx = artifact_indices(:)'
            start_idx = max(1, idx - half_win);
            end_idx = min(length(artifact_mask), idx + half_win);
            artifact_mask_dilated(start_idx:end_idx) = true;
        end
        
        artifact_mask = artifact_mask_dilated;
    end
end

% Step 2: Ensure contiguous blocks are merged (should already be done, but ensure)
% Find all artifact segments and merge any that are adjacent
artifact_segments = find_artifact_segments(artifact_mask);
if ~isempty(artifact_segments)
    % Reconstruct mask from segments (this ensures contiguous blocks)
    artifact_mask_merged = false(size(artifact_mask));
    for s = 1:size(artifact_segments, 1)
        start_idx = artifact_segments(s, 1);
        end_idx = artifact_segments(s, 2);
        artifact_mask_merged(start_idx:end_idx) = true;
    end
    artifact_mask = artifact_mask_merged;
end

% Step 3: Add padding around artifacts
% Expand each artifact segment by pre_pad before and post_pad after
pre_pad_samples = round(pre_pad_sec * fs);
post_pad_samples = round(post_pad_sec * fs);

if pre_pad_samples > 0 || post_pad_samples > 0
    artifact_mask_expanded = artifact_mask;
    
    % Find artifact segments
    segments = find_artifact_segments(artifact_mask);
    
    for s = 1:size(segments, 1)
        seg_start = segments(s, 1);
        seg_end = segments(s, 2);
        
        % Expand segment with padding
        expanded_start = max(1, seg_start - pre_pad_samples);
        expanded_end = min(length(artifact_mask), seg_end + post_pad_samples);
        
        % Mark expanded region as artifact
        artifact_mask_expanded(expanded_start:expanded_end) = true;
    end
    
    artifact_mask = artifact_mask_expanded;
end

cleaning_info.artifact_mask_processed = artifact_mask;
cleaning_info.num_segments_processed = size(find_artifact_segments(artifact_mask), 1);

clean_mask = ~artifact_mask;

n_clean = sum(clean_mask);
n_artifact = sum(artifact_mask);
n_original_artifact = sum(cleaning_info.artifact_mask_original);

cleaning_info.cleaned_samples = n_clean;
cleaning_info.removed_samples = n_artifact;
cleaning_info.removed_pct = 100 * n_artifact / n_data;
cleaning_info.removed_pct_original = 100 * n_original_artifact / n_data;  % Before processing
cleaning_info.num_segments = cleaning_info.num_segments_processed;
cleaning_info.applied = true;

% Diagnostic warning if no artifacts remain after processing (shouldn't happen, but check)
if n_artifact == 0 && n_original_artifact > 0
    warning('artifact_cleaning:processingRemovedAll', ...
        'WARNING: All artifacts were removed during processing for trial %d! Original: %d artifacts, Processed: 0. This may indicate a bug.', ...
        trial_idx, n_original_artifact);
end

% Diagnostic warning if very little data remains
if cleaning_info.removed_pct > 50
    warning('artifact_cleaning:excessiveRemoval', ...
        'WARNING: Trial %d has %.1f%% of data removed by artifact cleaning. This may affect analysis quality.', ...
        trial_idx, cleaning_info.removed_pct);
end

% Apply cleaning - HARD REMOVAL (samples dropped entirely, not NaN-masked)
% IMPORTANT: This creates discontinuities in time, but spectral analysis
% methods (PSD via pwelch, coherence via mscohere) work on the concatenated
% clean segments. These functions handle gaps by processing segments independently
% and then averaging, so hard removal is appropriate and more efficient than NaN masking.
%
% The cleaned data will have:
%   - Shorter time vector (only clean samples)
%   - Discontinuous time (gaps where artifacts were removed)
%   - All signals (LFP, GEVI, speed) aligned to the same indices
cleaned_data.time = data.time(clean_mask);
cleaned_data.lfp = data.lfp(clean_mask);
cleaned_data.gevi = data.gevi(clean_mask);
cleaned_data.speed = data.speed(clean_mask);
cleaned_data.fs = data.fs;  % Sampling rate unchanged
cleaned_data.duration = length(cleaned_data.time) / data.fs;  % Updated duration (actual clean data duration)

end

function [artifact_mask, mask_info] = load_artifact_mask_diagnostic(trial_path, trial_idx)
%LOAD_ARTIFACT_MASK_DIAGNOSTIC Load artifact mask (EXACT COPY from spectral_analysis.m)
%  This is an exact copy of load_artifact_mask (lines 2224-2343)

artifact_mask = [];
mask_info = struct('num_segments', 0, 'artifact_pct', 0);

trial_dir = fileparts(trial_path);

% First try: session folder (parent of trial folder)
session_dir = fileparts(trial_dir);
artifact_files = dir(fullfile(session_dir, '*_artifact_removal.mat'));

% Second try: directly in trial folder
if isempty(artifact_files)
    artifact_files = dir(fullfile(trial_dir, '*_artifact_removal.mat'));
end

% Third try: parent's parent
if isempty(artifact_files)
    parent_dir = fileparts(session_dir);
    artifact_files = dir(fullfile(parent_dir, '*_artifact_removal.mat'));
end

if isempty(artifact_files)
    return;
end

try
    artifact_path = fullfile(artifact_files(1).folder, artifact_files(1).name);
    loaded = load(artifact_path, 'ArtifactInfo');
    
    if ~isfield(loaded, 'ArtifactInfo')
        warning('artifact_cleaning:invalidFile', ...
            'ArtifactInfo struct not found in %s', artifact_path);
        return;
    end
    
    info = loaded.ArtifactInfo;
    
    if ~isfield(info, 'trials')
        warning('artifact_cleaning:noTrialsField', ...
            'ArtifactInfo has no "trials" field in %s', artifact_path);
        return;
    end
    
    if length(info.trials) < trial_idx
        warning('artifact_cleaning:trialNotFound', ...
            'Trial %d not found in artifact info (has %d trials, file: %s).', ...
            trial_idx, length(info.trials), artifact_path);
        return;
    end
    
    trial_info = info.trials(trial_idx);
    
    if isfield(trial_info, 'artifact_mask')
        artifact_mask = trial_info.artifact_mask(:);
        
        if isempty(artifact_mask)
            warning('artifact_cleaning:emptyMask', ...
                'Artifact mask for trial %d is empty in %s', trial_idx, artifact_path);
            return;
        end
        
        if ~islogical(artifact_mask)
            artifact_mask = logical(artifact_mask);
        end
        
        mask_info.artifact_pct = trial_info.artifact_pct;
        
        if isfield(trial_info, 'num_artifacts')
            mask_info.num_segments = trial_info.num_artifacts;
        elseif isfield(trial_info, 'artifact_segments')
            mask_info.num_segments = size(trial_info.artifact_segments, 1);
        else
            segments = find_artifact_segments(artifact_mask);
            mask_info.num_segments = size(segments, 1);
        end
    else
        warning('artifact_cleaning:noMask', ...
            'No artifact_mask field in trial %d of %s', trial_idx, artifact_path);
    end
    
catch ME
    warning('artifact_cleaning:loadFailed', ...
        'Failed to load artifact mask: %s', ME.message);
end

end

function segments = find_artifact_segments(artifact_mask)
%FIND_ARTIFACT_SEGMENTS Find contiguous artifact segments (EXACT COPY from spectral_analysis.m)

segments = [];

if isempty(artifact_mask) || ~any(artifact_mask)
    return;
end

mask_diff = diff([false; artifact_mask(:); false]);
starts = find(mask_diff == 1);
ends = find(mask_diff == -1) - 1;

if ~isempty(starts) && length(starts) == length(ends)
    segments = [starts, ends];
end

end
