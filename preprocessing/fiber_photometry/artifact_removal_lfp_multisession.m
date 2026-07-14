%% ============================================================================
%  LFP ARTIFACT REMOVAL TOOL - MULTI-ANIMAL, MULTI-SESSION VIEW
%  ============================================================================
%  Identifies and removes artifacts in LFP data across multiple sessions for
%  ALL animals. Generates one figure per animal, with all sessions as rows
%  and concatenated trials within each session.
%
%  FEATURES:
%    - Multi-animal processing (one figure per animal)
%    - Multi-session visualization (sessions as rows)
%    - Concatenated trials within each session
%    - Dotted vertical lines for trial boundaries
%    - Per-trial contamination statistics
%    - Interactive user approval system
%    - Non-destructive output (saves artifact info file)
%
%  ANIMALS PROCESSED:
%    Defined below in ALL_ANIMALS -- replace with your own cohort.
%
%  CREATED: 2025
%  ============================================================================

close all; clear; clc;

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'config'));
setup_lab_paths();

%% ============================================================================
%  ANIMAL DATABASE  -- EDIT THIS to your own cohort
%  ============================================================================
% One entry per animal. `sessions` is an Nx2 cell: {RECORDING_DATE, RECORDING_ID}.
% A session can also be a *combined* recording split across multiple IDs (see
% Animal01's second row) -- pass a 1x1 cell containing a cell array of IDs,
% e.g. {{'R1', 'R2'}}, when a session's trials were split across two runs.

BASE_PATH = fullfile(lab_paths().data_root, 'FiberVoltageImaging');

ALL_ANIMALS = struct();

% -------------------------------------------------------------------------
% ANIMAL 1 -- example with a combined (split) session
% -------------------------------------------------------------------------
ALL_ANIMALS(1).mouse_name = 'Animal01';
ALL_ANIMALS(1).base_path = BASE_PATH;
ALL_ANIMALS(1).sessions = {
    '01_01_25', 'R1';           % e.g. 6 trials
    '02_01_25', 'R1';           % e.g. 6 trials
    '03_01_25', {{'R1', 'R2'}}; % COMBINED: R1 + R2 recorded as one session
};

% -------------------------------------------------------------------------
% ANIMAL 2
% -------------------------------------------------------------------------
ALL_ANIMALS(2).mouse_name = 'Animal02';
ALL_ANIMALS(2).base_path = BASE_PATH;
ALL_ANIMALS(2).sessions = {
    '10_01_25', 'R1';
    '11_01_25', 'R1';
};

%% USER SELECTION - Which animals to process
% Set to empty [] to process ALL animals, or specify indices, e.g. [1, 2].
ANIMALS_TO_PROCESS = [];

% If empty, process all
if isempty(ANIMALS_TO_PROCESS)
    ANIMALS_TO_PROCESS = 1:length(ALL_ANIMALS);
end

%% ARTIFACT DETECTION PARAMETERS
% -------------------------------------------------------------------------
% LITERATURE-BASED RECOMMENDATIONS:
%
% 1. MAD (Median Absolute Deviation) - RECOMMENDED PRIMARY METHOD
%    - Most robust to outliers (unlike standard deviation)
%    - Used by Buzsáki lab, Harris lab for LFP artifact rejection
%    - Reference: Hampel (1974), Rousseeuw & Croux (1993)
%    - Threshold 4-6 MAD is typical for neural data
%
% 2. Rolling Variance - OPTIONAL (for sustained artifacts)
%    - Detects prolonged high-amplitude periods
%    - Useful for movement artifacts that span multiple samples
%
% 3. Z-score - NOT RECOMMENDED as sole method
%    - Sensitive to the very outliers you're trying to detect
%    - MAD is preferred; Z-score included for compatibility
%
% RECOMMENDED SETTINGS:
%   - For most data: MAD only (cfg.use_mad_detection = true, others false)
%   - For noisy data: MAD + Variance (both true)
% -------------------------------------------------------------------------

% MAD-based detection (PRIMARY - most robust)
cfg.mad_threshold = 6.0;              % Samples > threshold × MAD are flagged
                                       % Lower = more sensitive (try 4.0-5.0 if missing artifacts)
                                       % Higher = less sensitive (try 7.0-8.0 if over-detecting)
                                       % Current: 6.0 = moderate sensitivity (good default)

% Rolling variance detection (for sustained artifacts)
cfg.variance_window_sec = 0.5;         % Window size (seconds)
cfg.variance_threshold_factor = 4.0;   % Flag if variance > factor × median variance

% Z-score detection (less robust, included for compatibility)
cfg.zscore_threshold = 4.0;            % Flag samples with |z-score| > threshold

% Detection method selection - RECOMMENDED: MAD only
cfg.use_mad_detection = true;          % PRIMARY method - always use
cfg.use_variance_detection = false;    % Enable for sustained/movement artifacts
cfg.use_zscore_detection = false;      % Not recommended as primary method
cfg.combine_method = 'union';          % 'union' (any flags) or 'intersection' (all must flag)

% Temporal processing
cfg.extend_artifact_sec = 0.05;        % Extend artifacts by this amount on each side (50ms)
                                        % Set to 0 for no extension (less strict)
cfg.min_artifact_duration_sec = 0.05;  % Minimum artifact duration (50ms)

%% TRIAL EXCLUSION CRITERIA
% Literature: Buzsáki & Mizuseki (2014), typical range 20-50%
cfg.trial_exclusion_threshold = 0.30;  % Exclude trials with >30% contamination

%% OUTPUT OPTIONS
cfg.save_artifact_info = true;
cfg.output_suffix = '_artifact_removal';

%% VISUALIZATION OPTIONS
cfg.figure_width = 1900;
cfg.figure_height = 1000;

% Publication-quality colors (matching plot_figure1_gevi_lfp.py)
cfg.color_lfp = [0.35, 0.25, 0.45];        % Dark purple for LFP
cfg.color_motion = [0.993, 0.7, 0.4];      % Warm orange for motion (same as Python)
cfg.artifact_color = [0.9, 0.15, 0.15];    % Bright red for artifact shading
cfg.artifact_alpha = 0.35;                  % Semi-transparent
cfg.trial_boundary_color = [0.2, 0.2, 0.2]; % Dark gray for trial boundaries
cfg.trial_boundary_style = '--';           % Dashed line (more visible than dotted)
cfg.trial_boundary_width = 1.5;

%% ============================================================================
%  MAIN PROCESSING LOOP - ITERATE OVER ALL ANIMALS
%  ============================================================================

fprintf('\n');
fprintf('############################################################\n');
fprintf('  LFP ARTIFACT REMOVAL - MULTI-ANIMAL PROCESSING\n');
fprintf('############################################################\n');
fprintf('  Animals to process: %d\n', length(ANIMALS_TO_PROCESS));
for a_idx = ANIMALS_TO_PROCESS
    fprintf('    %d. %s\n', a_idx, ALL_ANIMALS(a_idx).mouse_name);
end
fprintf('############################################################\n\n');

% Store all animal data for final approval
all_animal_data = struct([]);

for animal_iter = 1:length(ANIMALS_TO_PROCESS)
    animal_idx = ANIMALS_TO_PROCESS(animal_iter);
    
    % Get current animal configuration
    MOUSE_NAME = ALL_ANIMALS(animal_idx).mouse_name;
    DATA_BASE_PATH = ALL_ANIMALS(animal_idx).base_path;
    SESSIONS = ALL_ANIMALS(animal_idx).sessions;
    
    fprintf('\n============================================================\n');
    fprintf('  ANIMAL %d/%d: %s\n', animal_iter, length(ANIMALS_TO_PROCESS), MOUSE_NAME);
    fprintf('============================================================\n');
    fprintf('  Base Path: %s\n', DATA_BASE_PATH);
    fprintf('  Sessions: %d\n', size(SESSIONS, 1));
    for s = 1:size(SESSIONS, 1)
        rec_id = SESSIONS{s, 2};
        if iscell(rec_id)
            % Combined session
            fprintf('    - %s-combined (%s)\n', SESSIONS{s, 1}, strjoin(rec_id{1}, '+'));
        else
            fprintf('    - %s-%s\n', SESSIONS{s, 1}, rec_id);
        end
    end
    fprintf('============================================================\n\n');
    
    %% ========================================================================
    %  LOAD ALL SESSION DATA FOR THIS ANIMAL
    %  ========================================================================
    
    fprintf('Loading session data for %s...\n\n', MOUSE_NAME);
    
    num_sessions = size(SESSIONS, 1);
    session_data = struct([]);
    
    for sess_idx = 1:num_sessions
        recording_date = SESSIONS{sess_idx, 1};
        recording_id_raw = SESSIONS{sess_idx, 2};
        
        % Handle combined sessions (multiple recording IDs)
        if iscell(recording_id_raw)
            % Combined session: e.g., {{'R1', 'R2'}}
            recording_ids = recording_id_raw{1};  % Extract the cell array of IDs
            session_id = sprintf('%s-combined', recording_date);
            is_combined = true;
        else
            % Single recording session
            recording_ids = {recording_id_raw};
            session_id = sprintf('%s-%s', recording_date, recording_id_raw);
            is_combined = false;
        end
        
        fprintf('Session %d/%d: %s', sess_idx, num_sessions, session_id);
        if is_combined
            fprintf(' (combined: %s)\n', strjoin(recording_ids, ' + '));
        else
            fprintf('\n');
        end
        
        % Collect trial files from all recording IDs in this session
        trial_files = [];
        data_paths = {};  % Store paths for each recording
        
        for r = 1:length(recording_ids)
            rec_id = recording_ids{r};
            rec_session_id = sprintf('%s-%s', recording_date, rec_id);
            data_path = fullfile(DATA_BASE_PATH, MOUSE_NAME, 'Fiber_Voltage_Processed', rec_session_id);
            data_paths{r} = data_path;
            
            if ~exist(data_path, 'dir')
                warning('  Path not found: %s', data_path);
                continue;
            end
            
            % Find trial files for this recording
            % Pattern 1: Files in Trial subfolders (standard structure from preprocessing)
            %   e.g., Trial1_fov1_baselineRecording_60sec_1/<MOUSE_NAME>-<DATE>-<ID>_Trial1_FiberPhotometry_Analysis.mat
            rec_trial_files = dir(fullfile(data_path, 'Trial*', '*_FiberPhotometry_Analysis.mat'));
            
            if isempty(rec_trial_files)
                % Pattern 2: Files directly in session folder (alternate structure)
                rec_trial_files = dir(fullfile(data_path, '*_FiberPhotometry_Analysis.mat'));
            end
            
            if isempty(rec_trial_files)
                % Pattern 3: Single file with session name only (no Trial subfolder)
                rec_trial_files = dir(fullfile(data_path, sprintf('%s-%s-%s_FiberPhotometry_Analysis.mat', ...
                    MOUSE_NAME, recording_date, rec_id)));
            end
            
            if ~isempty(rec_trial_files)
                trial_files = [trial_files; rec_trial_files]; %#ok<AGROW>
                fprintf('    %s: %d trial(s)\n', rec_id, length(rec_trial_files));
            else
                fprintf('    %s: No trial files found\n', rec_id);
            end
        end
        
        if isempty(trial_files)
            warning('  No trial files found for session');
            continue;
        end
        
        % Sort trials by name
        [~, sort_idx] = sort({trial_files.name});
        trial_files = trial_files(sort_idx);
        
        % Use first valid data path for saving artifact info
        data_path = data_paths{1};
        
        num_trials_in_session = length(trial_files);
        fprintf('  Total: %d trial(s)\n', num_trials_in_session);
        
        % Initialize session
        session_data(sess_idx).session_id = session_id;
        session_data(sess_idx).recording_date = recording_date;
        session_data(sess_idx).recording_ids = recording_ids;
        session_data(sess_idx).is_combined = is_combined;
        session_data(sess_idx).data_paths = data_paths;
        session_data(sess_idx).data_path = data_path;  % Primary path for saving
        session_data(sess_idx).num_trials = num_trials_in_session;
        session_data(sess_idx).trials = struct([]);
        
        % Concatenated data for this session
        concat_lfp = [];
        concat_motion = [];
        concat_time = [];
        trial_boundaries = [0]; % Sample indices where trials start
        fs = 0;  % Will be set from first valid trial
        
        for t = 1:num_trials_in_session
            trial_path = fullfile(trial_files(t).folder, trial_files(t).name);
            fprintf('    Trial %d: %s...', t, trial_files(t).name);
            
            try
                loaded = load(trial_path);
                FPA = loaded.FiberPhotometryAnalysis;
                
                % Extract data
                if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
                    time_vec = FPA.time.time_vector_seconds(:);
                else
                    error('Time vector not found');
                end
                
                if isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
                    fs = FPA.parameters.sampling_rate;
                elseif isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
                    fs = FPA.time.sampling_rate;
                else
                    fs = 1 / median(diff(time_vec));
                end
                
                if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_HP')
                    lfp = FPA.ephys.lfp_raw_aligned_HP(:);
                elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_mPFC')
                    lfp = FPA.ephys.lfp_raw_aligned_mPFC(:);
                else
                    error('LFP not found');
                end
                
                if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
                    motion = FPA.ephys.running_velocity_smooth(:);
                elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
                    motion = FPA.ephys.running_velocity(:);
                else
                    motion = zeros(size(lfp));
                end
                
                % Align lengths
                n = min([length(time_vec), length(lfp), length(motion)]);
                time_vec = time_vec(1:n);
                lfp = lfp(1:n);
                motion = motion(1:n);
                
                % Store trial info
                session_data(sess_idx).trials(t).filename = trial_files(t).name;
                session_data(sess_idx).trials(t).filepath = trial_path;
                session_data(sess_idx).trials(t).num_samples = n;
                session_data(sess_idx).trials(t).duration_sec = time_vec(end) - time_vec(1);
                session_data(sess_idx).trials(t).fs = fs;
                session_data(sess_idx).trials(t).local_indices = (1:n)';  % Indices within trial
                
                % Calculate time offset for concatenation
                if isempty(concat_time)
                    time_offset = 0;
                else
                    time_offset = concat_time(end) + 1/fs;
                end
                
                % Store concatenated indices
                session_data(sess_idx).trials(t).concat_start_idx = length(concat_lfp) + 1;
                session_data(sess_idx).trials(t).concat_end_idx = length(concat_lfp) + n;
                
                % Concatenate
                concat_lfp = [concat_lfp; lfp]; %#ok<AGROW>
                concat_motion = [concat_motion; motion]; %#ok<AGROW>
                concat_time = [concat_time; time_vec + time_offset]; %#ok<AGROW>
                trial_boundaries = [trial_boundaries; length(concat_lfp)]; %#ok<AGROW>
                
                fprintf(' OK (%.1f sec)\n', session_data(sess_idx).trials(t).duration_sec);
                
            catch ME
                warning('artifact_removal:loadFailed', 'Failed: %s', ME.message);
                continue;
            end
        end
        
        % Store session-level concatenated data
        session_data(sess_idx).concat_lfp = concat_lfp;
        session_data(sess_idx).concat_motion = concat_motion;
        session_data(sess_idx).concat_time = concat_time;
        session_data(sess_idx).trial_boundaries = trial_boundaries;
        session_data(sess_idx).fs = fs;
        session_data(sess_idx).total_samples = length(concat_lfp);
        if fs > 0 && ~isempty(concat_lfp)
            session_data(sess_idx).total_duration_sec = length(concat_lfp) / fs;
        else
            session_data(sess_idx).total_duration_sec = 0;
        end
        
        fprintf('  Session total: %.1f sec (%d samples)\n\n', ...
            session_data(sess_idx).total_duration_sec, session_data(sess_idx).total_samples);
    end
    
    % Remove empty sessions
    valid_sessions = arrayfun(@(x) ~isempty(x.concat_lfp), session_data);
    session_data = session_data(valid_sessions);
    num_sessions = length(session_data);
    
    if num_sessions == 0
        warning('No valid sessions loaded for animal %s. Skipping...', MOUSE_NAME);
        continue;
    end

    %% ========================================================================
    %  ARTIFACT DETECTION FOR THIS ANIMAL
    %  ========================================================================
    
    fprintf('------------------------------------------------------------\n');
    fprintf('  DETECTING ARTIFACTS for %s\n', MOUSE_NAME);
    fprintf('------------------------------------------------------------\n');
    fprintf('Methods: MAD(%.1f) + Variance(%.1f) + Z-score(%.1f)\n', ...
        cfg.mad_threshold, cfg.variance_threshold_factor, cfg.zscore_threshold);
    fprintf('Combine: %s\n\n', cfg.combine_method);
    
    for sess_idx = 1:num_sessions
    fprintf('Session %s:\n', session_data(sess_idx).session_id);
    
    lfp = session_data(sess_idx).concat_lfp;
    fs = session_data(sess_idx).fs;
    n = length(lfp);
    
    % Initialize masks
    artifact_mad = false(n, 1);
    artifact_var = false(n, 1);
    artifact_zscore = false(n, 1);
    
    % MAD-based detection
    if cfg.use_mad_detection
        med = median(lfp);
        mad_val = median(abs(lfp - med));
        if mad_val > 0
            artifact_mad = abs(lfp - med) > cfg.mad_threshold * mad_val;
        end
    end
    
    % Rolling variance detection
    if cfg.use_variance_detection
        win_samples = max(3, round(cfg.variance_window_sec * fs));
        lfp_sq = lfp.^2;
        kernel = ones(win_samples, 1) / win_samples;
        mean_sq = conv(lfp_sq, kernel, 'same');
        mean_val = conv(lfp, kernel, 'same');
        rolling_var = max(0, mean_sq - mean_val.^2);
        med_var = median(rolling_var);
        if med_var > 0
            artifact_var = rolling_var > cfg.variance_threshold_factor * med_var;
        end
    end
    
    % Z-score detection
    if cfg.use_zscore_detection
        lfp_z = zscore(lfp);
        artifact_zscore = abs(lfp_z) > cfg.zscore_threshold;
    end
    
    % Combine
    if strcmpi(cfg.combine_method, 'union')
        artifact_mask = artifact_mad | artifact_var | artifact_zscore;
    else
        artifact_mask = artifact_mad & artifact_var & artifact_zscore;
    end
    
    % Extend artifacts
    extend_samples = round(cfg.extend_artifact_sec * fs);
    if extend_samples > 0
        for i = 1:extend_samples
            artifact_mask = artifact_mask | [false; artifact_mask(1:end-1)] | [artifact_mask(2:end); false];
        end
    end
    
    % Remove short artifacts
    min_samples = round(cfg.min_artifact_duration_sec * fs);
    artifact_mask = remove_short_segments(artifact_mask, min_samples);
    
    % Store session-level results
    session_data(sess_idx).artifact_mask = artifact_mask;
    session_data(sess_idx).artifact_pct = 100 * sum(artifact_mask) / n;
    session_data(sess_idx).artifact_segments = find_segments(artifact_mask);
    
    fprintf('  Overall: %.1f%% contaminated\n', session_data(sess_idx).artifact_pct);
    
    % Calculate per-trial contamination
    for t = 1:session_data(sess_idx).num_trials
        start_idx = session_data(sess_idx).trials(t).concat_start_idx;
        end_idx = session_data(sess_idx).trials(t).concat_end_idx;
        
        trial_mask = artifact_mask(start_idx:end_idx);
        trial_n = end_idx - start_idx + 1;
        
        session_data(sess_idx).trials(t).artifact_mask = trial_mask;
        session_data(sess_idx).trials(t).artifact_pct = 100 * sum(trial_mask) / trial_n;
        session_data(sess_idx).trials(t).artifact_segments = find_segments(trial_mask);
        session_data(sess_idx).trials(t).num_artifacts = size(session_data(sess_idx).trials(t).artifact_segments, 1);
        
        % Recommendation
        if session_data(sess_idx).trials(t).artifact_pct > cfg.trial_exclusion_threshold * 100
            session_data(sess_idx).trials(t).recommendation = 'EXCLUDE';
        else
            session_data(sess_idx).trials(t).recommendation = 'CLEAN';
        end
        
        fprintf('    Trial %d: %.1f%% -> %s\n', t, ...
            session_data(sess_idx).trials(t).artifact_pct, ...
            session_data(sess_idx).trials(t).recommendation);
    end
    fprintf('\n');
end

    %% ========================================================================
    %  MULTI-SESSION VISUALIZATION FOR THIS ANIMAL (Publication-Quality)
    %  ========================================================================
    %  Each session has TWO separate rows:
    %    - Top row (purple): Z-scored LFP with artifact shading
    %    - Bottom row (orange): Motion trace
    %  This provides clean, non-overlapping visualization.
    
    fprintf('------------------------------------------------------------\n');
    fprintf('  GENERATING VISUALIZATION for %s\n', MOUSE_NAME);
    fprintf('------------------------------------------------------------\n');
    
    fig = figure('Name', sprintf('[%d/%d] LFP Artifact Review - %s', animal_iter, length(ANIMALS_TO_PROCESS), MOUSE_NAME), ...
    'NumberTitle', 'off', ...
    'Position', [30, 30, cfg.figure_width, cfg.figure_height], ...
    'Color', 'w');

% Layout parameters
margin_left = 0.10;
margin_right = 0.03;
margin_bottom = 0.08;
margin_top = 0.06;
session_spacing = 0.025;      % Space between sessions
subrow_spacing = 0.005;       % Space between LFP and motion within a session

% Each session has 2 sub-rows (LFP + motion)
% LFP row is 70% of session height, motion is 30%
lfp_ratio = 0.7;
motion_ratio = 0.3;

% Calculate session block height
total_session_spacing = (num_sessions - 1) * session_spacing;
available_height = 1 - margin_top - margin_bottom - total_session_spacing;
session_block_height = available_height / num_sessions;

% Sub-row heights within each session
lfp_height = session_block_height * lfp_ratio - subrow_spacing/2;
motion_height = session_block_height * motion_ratio - subrow_spacing/2;

ax_width = 1 - margin_left - margin_right;

% Store handles
ax_lfp_handles = gobjects(num_sessions, 1);
ax_motion_handles = gobjects(num_sessions, 1);

for sess_idx = 1:num_sessions
    % Calculate position for this session block (top to bottom)
    session_block_bottom = margin_bottom + (num_sessions - sess_idx) * (session_block_height + session_spacing);
    
    % Get session data
    time = session_data(sess_idx).concat_time;
    lfp = session_data(sess_idx).concat_lfp;
    motion = session_data(sess_idx).concat_motion;
    artifact_mask = session_data(sess_idx).artifact_mask;
    
    % Z-score LFP for consistent scaling across sessions
    lfp_zscore = zscore(lfp);
    
    % =====================================================================
    % LFP ROW (TOP of session block)
    % =====================================================================
    lfp_bottom = session_block_bottom + motion_height + subrow_spacing;
    ax_lfp = axes('Parent', fig, 'Position', [margin_left, lfp_bottom, ax_width, lfp_height]);
    ax_lfp_handles(sess_idx) = ax_lfp;
    hold(ax_lfp, 'on');
    
    % Y-limits for z-scored LFP
    lfp_ylim = [-5, 5];
    
    % Plot ARTIFACT regions (behind traces)
    segments = session_data(sess_idx).artifact_segments;
    for s = 1:size(segments, 1)
        seg_start = time(segments(s, 1));
        seg_end = time(segments(s, 2));
        fill(ax_lfp, [seg_start, seg_end, seg_end, seg_start], ...
            [lfp_ylim(1), lfp_ylim(1), lfp_ylim(2), lfp_ylim(2)], ...
            cfg.artifact_color, 'FaceAlpha', cfg.artifact_alpha, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    % Plot TRIAL BOUNDARIES
    trial_boundaries = session_data(sess_idx).trial_boundaries;
    for b = 2:length(trial_boundaries)-1
        if trial_boundaries(b) <= length(time)
            boundary_time = time(trial_boundaries(b));
            line(ax_lfp, [boundary_time, boundary_time], lfp_ylim, ...
                'Color', cfg.trial_boundary_color, 'LineStyle', cfg.trial_boundary_style, ...
                'LineWidth', cfg.trial_boundary_width, 'HandleVisibility', 'off');
        end
    end
    
    % Plot LFP trace
    plot(ax_lfp, time, lfp_zscore, 'Color', cfg.color_lfp, 'LineWidth', 0.5);
    
    % Configure LFP axis
    xlim(ax_lfp, [time(1), time(end)]);
    ylim(ax_lfp, lfp_ylim);
    ax_lfp.YColor = cfg.color_lfp;
    ax_lfp.YTick = [-4, 0, 4];
    ax_lfp.YTickLabel = {'-4σ', '0', '4σ'};
    ax_lfp.FontSize = 8;
    ax_lfp.XTickLabel = [];  % No x-labels on LFP row
    ax_lfp.Box = 'on';
    ax_lfp.XColor = [0.7, 0.7, 0.7];
    
    % Session label
    ylabel(ax_lfp, sprintf('%s\n%.1f%% art.', ...
        strrep(session_data(sess_idx).session_id, '_', '-'), ...
        session_data(sess_idx).artifact_pct), ...
        'FontSize', 9, 'FontWeight', 'bold', 'Color', 'k', 'Rotation', 0, ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    
    % Add TRIAL LABELS at top
    for t = 1:session_data(sess_idx).num_trials
        start_idx = session_data(sess_idx).trials(t).concat_start_idx;
        end_idx = session_data(sess_idx).trials(t).concat_end_idx;
        trial_start_time = time(start_idx);
        trial_end_time = time(min(end_idx, length(time)));
        mid_time = (trial_start_time + trial_end_time) / 2;
        
        if strcmp(session_data(sess_idx).trials(t).recommendation, 'EXCLUDE')
            label_color = [0.85, 0.15, 0.15];
            label_bg = [1, 0.92, 0.92];
        else
            label_color = [0.15, 0.5, 0.15];
            label_bg = [0.92, 1, 0.92];
        end
        
        label_str = sprintf('T%d: %.0f%%', t, session_data(sess_idx).trials(t).artifact_pct);
        text(ax_lfp, mid_time, lfp_ylim(2) * 0.92, label_str, ...
            'Color', label_color, 'FontSize', 7, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
            'BackgroundColor', label_bg, 'EdgeColor', label_color, ...
            'Margin', 1, 'Clipping', 'on');
    end
    
    hold(ax_lfp, 'off');
    
    % =====================================================================
    % MOTION ROW (BOTTOM of session block)
    % =====================================================================
    ax_motion = axes('Parent', fig, 'Position', [margin_left, session_block_bottom, ax_width, motion_height]);
    ax_motion_handles(sess_idx) = ax_motion;
    hold(ax_motion, 'on');
    
    % Motion y-limits
    motion_min = 0;
    motion_max = max(prctile(motion, 99.5), 0.5);
    motion_ylim = [motion_min, motion_max * 1.15];
    
    % Plot artifact regions on motion too (for visual alignment)
    for s = 1:size(segments, 1)
        seg_start = time(segments(s, 1));
        seg_end = time(segments(s, 2));
        fill(ax_motion, [seg_start, seg_end, seg_end, seg_start], ...
            [motion_ylim(1), motion_ylim(1), motion_ylim(2), motion_ylim(2)], ...
            cfg.artifact_color, 'FaceAlpha', cfg.artifact_alpha * 0.5, ...
            'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    
    % Plot trial boundaries on motion
    for b = 2:length(trial_boundaries)-1
        if trial_boundaries(b) <= length(time)
            boundary_time = time(trial_boundaries(b));
            line(ax_motion, [boundary_time, boundary_time], motion_ylim, ...
                'Color', cfg.trial_boundary_color, 'LineStyle', cfg.trial_boundary_style, ...
                'LineWidth', cfg.trial_boundary_width * 0.7, 'HandleVisibility', 'off');
        end
    end
    
    % Plot motion trace
    plot(ax_motion, time, motion, 'Color', cfg.color_motion, 'LineWidth', 0.6);
    
    % Configure motion axis
    xlim(ax_motion, [time(1), time(end)]);
    ylim(ax_motion, motion_ylim);
    ax_motion.YColor = cfg.color_motion;
    ax_motion.FontSize = 8;
    ax_motion.Box = 'on';
    ax_motion.XColor = [0.5, 0.5, 0.5];
    
    % Only show x-labels on bottom session
    if sess_idx == num_sessions
        xlabel(ax_motion, 'Time (s)', 'FontSize', 10);
    else
        ax_motion.XTickLabel = [];
    end
    
    % Motion label (compact)
    ylabel(ax_motion, 'Speed', 'FontSize', 8, 'Color', cfg.color_motion, ...
        'Rotation', 0, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle');
    
    hold(ax_motion, 'off');
end

% =========================================================================
% TITLE
% =========================================================================
annotation(fig, 'textbox', [0, 0.95, 1, 0.04], ...
    'String', sprintf('%s - Multi-Session LFP Artifact Review | MAD Threshold: %.1f | Exclusion: %.0f%%', ...
        MOUSE_NAME, cfg.mad_threshold, cfg.trial_exclusion_threshold * 100), ...
    'FontSize', 13, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'EdgeColor', 'none');

% =========================================================================
% LEGEND (bottom of figure)
% =========================================================================
legend_ax = axes('Parent', fig, 'Position', [0.25, 0.005, 0.5, 0.04], 'Visible', 'off');
hold(legend_ax, 'on');
h1 = plot(legend_ax, NaN, NaN, 'Color', cfg.color_lfp, 'LineWidth', 2);
h2 = plot(legend_ax, NaN, NaN, 'Color', cfg.color_motion, 'LineWidth', 2);
h3 = fill(legend_ax, NaN, NaN, cfg.artifact_color, 'FaceAlpha', cfg.artifact_alpha, 'EdgeColor', 'none');
h4 = line(legend_ax, NaN, NaN, 'Color', cfg.trial_boundary_color, 'LineStyle', cfg.trial_boundary_style, 'LineWidth', cfg.trial_boundary_width);
legend(legend_ax, [h1, h2, h3, h4], {'LFP (z-scored)', 'Motion/Speed', 'Artifact', 'Trial Boundary'}, ...
    'Orientation', 'horizontal', 'Location', 'north', 'FontSize', 9, 'Box', 'off');
hold(legend_ax, 'off');

drawnow;
fprintf('  Visualization complete for %s.\n\n', MOUSE_NAME);

% Store animal data and figure handle for later approval
all_animal_data(animal_iter).mouse_name = MOUSE_NAME;
all_animal_data(animal_iter).base_path = DATA_BASE_PATH;
all_animal_data(animal_iter).session_data = session_data;
all_animal_data(animal_iter).num_sessions = num_sessions;
all_animal_data(animal_iter).fig_handle = fig;

end  % END OF ANIMAL LOOP

%% ============================================================================
%  SUMMARY REPORT - ALL ANIMALS
%  ============================================================================

fprintf('\n############################################################\n');
fprintf('  SUMMARY REPORT - ALL ANIMALS\n');
fprintf('############################################################\n\n');

for animal_iter = 1:length(all_animal_data)
    fprintf('============================================================\n');
    fprintf('ANIMAL: %s\n', all_animal_data(animal_iter).mouse_name);
    fprintf('============================================================\n');
    
    session_data = all_animal_data(animal_iter).session_data;
    num_sessions = all_animal_data(animal_iter).num_sessions;
    
    total_exclude_count = 0;
    total_trial_count = 0;
    
    for sess_idx = 1:num_sessions
        fprintf('  Session: %s (%.1f%% artifacts)\n', ...
            session_data(sess_idx).session_id, session_data(sess_idx).artifact_pct);
        
        exclude_trials = [];
        clean_trials = [];
        for t = 1:session_data(sess_idx).num_trials
            if strcmp(session_data(sess_idx).trials(t).recommendation, 'EXCLUDE')
                exclude_trials = [exclude_trials, t]; %#ok<AGROW>
            else
                clean_trials = [clean_trials, t]; %#ok<AGROW>
            end
        end
        
        total_exclude_count = total_exclude_count + length(exclude_trials);
        total_trial_count = total_trial_count + session_data(sess_idx).num_trials;
        
        if ~isempty(exclude_trials)
            fprintf('    -> EXCLUDE trials: %s\n', ...
                strjoin(arrayfun(@num2str, exclude_trials, 'UniformOutput', false), ', '));
        else
            fprintf('    -> All trials CLEAN\n');
        end
    end
    
    fprintf('  SUMMARY: %d/%d trials recommended for exclusion\n\n', total_exclude_count, total_trial_count);
end

%% ============================================================================
%  USER APPROVAL - ALL ANIMALS
%  ============================================================================

fprintf('############################################################\n');
fprintf('  USER APPROVAL - ALL ANIMALS\n');
fprintf('############################################################\n');
fprintf('Options:\n');
fprintf('  [1] ACCEPT all algorithm recommendations for ALL animals\n');
fprintf('  [2] MANUAL selection (per animal, per session)\n');
fprintf('  [3] CANCEL - exit without saving\n\n');

user_choice = input('Enter choice (1/2/3): ', 's');

switch user_choice
    case '1'
        fprintf('\nAccepted all algorithm recommendations for all animals.\n');
        user_decision = 'algorithm';
        
        % Apply algorithm recommendations to all animals
        for animal_iter = 1:length(all_animal_data)
            session_data = all_animal_data(animal_iter).session_data;
            for sess_idx = 1:all_animal_data(animal_iter).num_sessions
                for t = 1:session_data(sess_idx).num_trials
                    session_data(sess_idx).trials(t).final_decision = session_data(sess_idx).trials(t).recommendation;
                end
            end
            all_animal_data(animal_iter).session_data = session_data;
        end
        
    case '2'
        user_decision = 'manual';
        
        for animal_iter = 1:length(all_animal_data)
            fprintf('\n------------------------------------------------------------\n');
            fprintf('ANIMAL: %s\n', all_animal_data(animal_iter).mouse_name);
            fprintf('------------------------------------------------------------\n');
            
            session_data = all_animal_data(animal_iter).session_data;
            
            for sess_idx = 1:all_animal_data(animal_iter).num_sessions
                fprintf('\n  Session %s:\n', session_data(sess_idx).session_id);
                fprintf('    Algorithm recommendations: ');
                exclude_list = [];
                for t = 1:session_data(sess_idx).num_trials
                    if strcmp(session_data(sess_idx).trials(t).recommendation, 'EXCLUDE')
                        exclude_list = [exclude_list, t]; %#ok<AGROW>
                    end
                end
                if isempty(exclude_list)
                    fprintf('None\n');
                else
                    fprintf('%s\n', strjoin(arrayfun(@num2str, exclude_list, 'UniformOutput', false), ', '));
                end
                
                manual_input = input('    Trials to exclude (comma-separated, or Enter for default): ', 's');
                
                if isempty(manual_input)
                    % Use default
                    for t = 1:session_data(sess_idx).num_trials
                        session_data(sess_idx).trials(t).final_decision = session_data(sess_idx).trials(t).recommendation;
                    end
                else
                    excluded = str2num(manual_input); %#ok<ST2NM>
                    for t = 1:session_data(sess_idx).num_trials
                        if ismember(t, excluded)
                            session_data(sess_idx).trials(t).final_decision = 'EXCLUDE';
                        else
                            session_data(sess_idx).trials(t).final_decision = 'CLEAN';
                        end
                    end
                end
            end
            
            all_animal_data(animal_iter).session_data = session_data;
        end
        
    case '3'
        fprintf('\nCancelled. No changes saved.\n');
        return;
        
    otherwise
        fprintf('\nInvalid choice. Exiting.\n');
        return;
end

%% ============================================================================
%  SAVE ARTIFACT INFORMATION - ALL ANIMALS
%  ============================================================================

if cfg.save_artifact_info
    fprintf('\n############################################################\n');
    fprintf('  SAVING ARTIFACT INFORMATION - ALL ANIMALS\n');
    fprintf('############################################################\n');
    
    for animal_iter = 1:length(all_animal_data)
        MOUSE_NAME = all_animal_data(animal_iter).mouse_name;
        DATA_BASE_PATH = all_animal_data(animal_iter).base_path;
        session_data = all_animal_data(animal_iter).session_data;
        num_sessions = all_animal_data(animal_iter).num_sessions;
        fig = all_animal_data(animal_iter).fig_handle;
        
        fprintf('\n  ANIMAL: %s\n', MOUSE_NAME);
        
        for sess_idx = 1:num_sessions
            % Create artifact info structure
            ArtifactInfo = struct();
            ArtifactInfo.mouse_name = MOUSE_NAME;
            ArtifactInfo.session_id = session_data(sess_idx).session_id;
            ArtifactInfo.recording_date = session_data(sess_idx).recording_date;
            ArtifactInfo.recording_ids = session_data(sess_idx).recording_ids;
            ArtifactInfo.analysis_date = datestr(now);
            ArtifactInfo.config = cfg;
            ArtifactInfo.user_decision = user_decision;
            
            % Summary
            ArtifactInfo.summary.total_trials = session_data(sess_idx).num_trials;
            ArtifactInfo.summary.total_samples = session_data(sess_idx).total_samples;
            ArtifactInfo.summary.total_duration_sec = session_data(sess_idx).total_duration_sec;
            ArtifactInfo.summary.total_artifact_pct = session_data(sess_idx).artifact_pct;
            ArtifactInfo.summary.is_combined = session_data(sess_idx).is_combined;
            ArtifactInfo.summary.recording_ids = session_data(sess_idx).recording_ids;
            
            excluded_trials = [];
            for t = 1:session_data(sess_idx).num_trials
                if strcmp(session_data(sess_idx).trials(t).final_decision, 'EXCLUDE')
                    excluded_trials = [excluded_trials, t]; %#ok<AGROW>
                end
            end
            ArtifactInfo.summary.excluded_trials = excluded_trials;
            
            % Per-trial information
            for t = 1:session_data(sess_idx).num_trials
                ArtifactInfo.trials(t).trial_number = t;
                ArtifactInfo.trials(t).filename = session_data(sess_idx).trials(t).filename;
                ArtifactInfo.trials(t).filepath = session_data(sess_idx).trials(t).filepath;
                ArtifactInfo.trials(t).num_samples = session_data(sess_idx).trials(t).num_samples;
                ArtifactInfo.trials(t).duration_sec = session_data(sess_idx).trials(t).duration_sec;
                ArtifactInfo.trials(t).artifact_pct = session_data(sess_idx).trials(t).artifact_pct;
                ArtifactInfo.trials(t).num_artifacts = session_data(sess_idx).trials(t).num_artifacts;
                ArtifactInfo.trials(t).artifact_segments = session_data(sess_idx).trials(t).artifact_segments;
                ArtifactInfo.trials(t).artifact_mask = session_data(sess_idx).trials(t).artifact_mask;
                ArtifactInfo.trials(t).algorithm_recommendation = session_data(sess_idx).trials(t).recommendation;
                ArtifactInfo.trials(t).final_decision = session_data(sess_idx).trials(t).final_decision;
            end
            
            % Save file(s) - for combined sessions, save to ALL recording paths
            output_filename = sprintf('%s-%s%s.mat', ...
                MOUSE_NAME, session_data(sess_idx).session_id, cfg.output_suffix);
            
            if session_data(sess_idx).is_combined
                % Save to each recording's data path so spectral analysis can find it
                for dp = 1:length(session_data(sess_idx).data_paths)
                    if exist(session_data(sess_idx).data_paths{dp}, 'dir')
                        output_path = fullfile(session_data(sess_idx).data_paths{dp}, output_filename);
                        save(output_path, 'ArtifactInfo', '-v7.3');
                        fprintf('    Saved: %s (in %s)\n', output_filename, session_data(sess_idx).recording_ids{dp});
                    end
                end
            else
                output_path = fullfile(session_data(sess_idx).data_path, output_filename);
                save(output_path, 'ArtifactInfo', '-v7.3');
                fprintf('    Saved: %s\n', output_filename);
            end
        end
        
        % Save figure for this animal
        if isvalid(fig)
            fig_filename = sprintf('%s_all_sessions%s.png', MOUSE_NAME, cfg.output_suffix);
            fig_path = fullfile(DATA_BASE_PATH, MOUSE_NAME, 'Fiber_Voltage_Processed', fig_filename);
            saveas(fig, fig_path);
            fprintf('    Figure saved: %s\n', fig_filename);
        end
    end
end

fprintf('\n############################################################\n');
fprintf('  ARTIFACT REMOVAL COMPLETE - ALL %d ANIMALS PROCESSED\n', length(all_animal_data));
fprintf('############################################################\n\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================

function segments = find_segments(mask)
segments = [];
in_segment = false;
start_idx = 0;
for i = 1:length(mask)
    if mask(i) && ~in_segment
        in_segment = true;
        start_idx = i;
    elseif ~mask(i) && in_segment
        in_segment = false;
        segments = [segments; start_idx, i-1]; %#ok<AGROW>
    end
end
if in_segment
    segments = [segments; start_idx, length(mask)];
end
end

function mask = remove_short_segments(mask, min_samples)
segments = find_segments(mask);
for s = 1:size(segments, 1)
    seg_length = segments(s, 2) - segments(s, 1) + 1;
    if seg_length < min_samples
        mask(segments(s, 1):segments(s, 2)) = false;
    end
end
end
