%% ============================================================================
%  LFP ARTIFACT REMOVAL TOOL
%  ============================================================================
%  Identifies and removes artifacts in LFP data from fiber photometry
%  preprocessing pipeline. Uses multiple detection methods for robust
%  artifact identification.
%
%  FEATURES:
%    - Multi-method artifact detection (MAD, rolling variance, z-score)
%    - Interactive visualization with trial boundaries
%    - User approval system (accept algorithm / manual exclusion)
%    - Saves artifact info separately (non-destructive approach)
%    - Optional: Create cleaned data struct
%
%  DETECTION METHODS:
%    1. MAD (Median Absolute Deviation) - robust to outliers
%    2. Rolling variance - detects sustained high-amplitude periods
%    3. Z-score - flags extreme deviations
%
%  LITERATURE-BASED EXCLUSION:
%    - Default: Exclude trials with >30% contamination
%    - Based on: Buzsáki & Mizuseki (2014), Harris et al. practices
%    - User configurable threshold
%
%  CREATED: 2025
%  ============================================================================

close all; clear; clc;

addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'config'));
setup_lab_paths();

%% ============================================================================
%  USER CONFIGURATION - MODIFY THESE PARAMETERS
%  ============================================================================

%% ANIMAL AND SESSION CONFIGURATION  -- EDIT THESE FOR YOUR RECORDING
MOUSE_NAME = 'Animal01';                    % Your animal/mouse identifier
RECORDING_DATE = '01_01_25';                % Format: DD_MM_YY
RECORDING_ID = 'R1';                        % Recording ID (R1, R2, etc.)

%% DATA PATHS
% If you keep different projects/cohorts under different base folders, you can
% branch on MOUSE_NAME here (or, simpler, just point DATA_BASE_PATH at a single
% root below and skip the switch). Example pattern for multiple projects:
%
% switch MOUSE_NAME
%     case {'Animal01', 'Animal02'}
%         DATA_BASE_PATH = fullfile(lab_paths().data_root, 'ProjectA');
%     case {'Animal03'}
%         DATA_BASE_PATH = fullfile(lab_paths().data_root, 'ProjectB');
%     otherwise
%         DATA_BASE_PATH = lab_paths().data_root;
% end

DATA_BASE_PATH = fullfile(lab_paths().data_root, 'FiberVoltageImaging');

% Construct path automatically
DATA_PATH = fullfile(DATA_BASE_PATH, MOUSE_NAME, 'Fiber_Voltage_Processed', ...
    [RECORDING_DATE '-' RECORDING_ID]);

%% ARTIFACT DETECTION PARAMETERS
% All thresholds are configurable for fine-tuning

% MAD-based detection (primary method)
cfg.mad_threshold = 5.0;              % Samples > MAD_THRESHOLD × MAD are flagged
                                       % Literature: 4-6 MAD is typical for artifact rejection

% Rolling variance detection
cfg.variance_window_sec = 0.5;         % Window size for variance calculation (seconds)
cfg.variance_threshold_factor = 3.0;   % Flag if variance > factor × median variance

% Z-score detection
cfg.zscore_threshold = 4.0;            % Flag samples with |z-score| > threshold

% Detection method selection (combine for robustness)
cfg.use_mad_detection = true;          % Enable MAD-based detection
cfg.use_variance_detection = true;     % Enable rolling variance detection
cfg.use_zscore_detection = true;       % Enable z-score detection
cfg.combine_method = 'union';          % 'union' (any method flags) or 'intersection' (all must flag)

% Temporal extension (artifacts often have buildup/decay)
cfg.extend_artifact_sec = 0.1;         % Extend artifact regions by this amount on each side

% Minimum artifact duration (ignore very brief spikes)
cfg.min_artifact_duration_sec = 0.05;  % Minimum 50ms to be considered artifact

%% TRIAL EXCLUSION CRITERIA
% Based on literature: Buzsáki & Mizuseki (2014), common practices in Harris lab
cfg.trial_exclusion_threshold = 0.30;  % Exclude trials with >30% contamination
                                        % Range in literature: 20-50% depending on study

%% OUTPUT OPTIONS
cfg.save_artifact_info = true;         % Save artifact information file
cfg.save_cleaned_data = false;         % Save cleaned data struct (with artifacts removed)
cfg.output_suffix = '_artifact_removal'; % Suffix for output files

%% VISUALIZATION OPTIONS
cfg.show_interactive_plot = true;      % Show interactive artifact review plot
cfg.figure_width = 1800;               % Figure width in pixels
cfg.figure_height = 900;               % Figure height in pixels
cfg.artifact_color = [1, 0.2, 0.2, 0.3]; % RGBA for artifact shading
cfg.trial_boundary_color = [0.5, 0.5, 0.5]; % Color for trial boundary lines

%% ============================================================================
%  INITIALIZATION
%  ============================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  LFP ARTIFACT REMOVAL TOOL\n');
fprintf('============================================================\n');
fprintf('  Mouse: %s\n', MOUSE_NAME);
fprintf('  Date: %s\n', RECORDING_DATE);
fprintf('  Recording: %s\n', RECORDING_ID);
fprintf('  Data path: %s\n', DATA_PATH);
fprintf('============================================================\n\n');

% Validate path exists
if ~exist(DATA_PATH, 'dir')
    error('Data path does not exist: %s', DATA_PATH);
end

%% ============================================================================
%  LOAD TRIAL DATA
%  ============================================================================

fprintf('Loading trial data...\n');

% Find all trial files in the directory
% Pattern 1: Files in Trial subfolders (standard structure from preprocessing)
%   e.g., Trial1_fov1_baselineRecording_60sec_1/<MOUSE_NAME>-<DATE>-<ID>_Trial1_FiberPhotometry_Analysis.mat
trial_files = dir(fullfile(DATA_PATH, 'Trial*', '*_FiberPhotometry_Analysis.mat'));

if isempty(trial_files)
    % Pattern 2: Files directly in session folder (alternate structure)
    trial_files = dir(fullfile(DATA_PATH, '*_FiberPhotometry_Analysis.mat'));
end

if isempty(trial_files)
    % Pattern 3: Single file with session name only (no Trial subfolder)
    trial_files = dir(fullfile(DATA_PATH, sprintf('%s-%s-%s_FiberPhotometry_Analysis.mat', ...
        MOUSE_NAME, RECORDING_DATE, RECORDING_ID)));
end

if isempty(trial_files)
    error('No trial files found in: %s\nTried patterns:\n  1. Trial*/\n  2. Direct files\n  3. Single session file', DATA_PATH);
end

% Sort by name to ensure correct trial order
[~, sort_idx] = sort({trial_files.name});
trial_files = trial_files(sort_idx);

num_trials = length(trial_files);
fprintf('  Found %d trial(s)\n\n', num_trials);

% Initialize storage
trial_data = struct([]);
total_samples = 0;

for t = 1:num_trials
    trial_path = fullfile(trial_files(t).folder, trial_files(t).name);
    fprintf('  Loading Trial %d: %s...', t, trial_files(t).name);
    
    try
        loaded = load(trial_path);
        
        if ~isfield(loaded, 'FiberPhotometryAnalysis')
            warning('FiberPhotometryAnalysis not found in %s', trial_path);
            continue;
        end
        
        FPA = loaded.FiberPhotometryAnalysis;
        
        % Extract required fields
        trial_data(t).filename = trial_files(t).name;
        trial_data(t).filepath = trial_path;
        
        % Time and sampling rate
        if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
            trial_data(t).time = FPA.time.time_vector_seconds(:);
        else
            error('Time vector not found');
        end
        
        if isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
            trial_data(t).fs = FPA.parameters.sampling_rate;
        elseif isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
            trial_data(t).fs = FPA.time.sampling_rate;
        else
            trial_data(t).fs = 1 / median(diff(trial_data(t).time));
        end
        
        % LFP data
        if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_HP')
            trial_data(t).lfp = FPA.ephys.lfp_raw_aligned_HP(:);
        elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_mPFC')
            trial_data(t).lfp = FPA.ephys.lfp_raw_aligned_mPFC(:);
        else
            error('LFP data not found');
        end
        
        % GEVI/Fiber data (for later trimming if needed)
        if isfield(FPA, 'signals') && isfield(FPA.signals, 'final_processed_traces')
            trial_data(t).gevi = FPA.signals.final_processed_traces;
        elseif isfield(FPA, 'signals') && isfield(FPA.signals, 'deltaF_F_traces')
            trial_data(t).gevi = FPA.signals.deltaF_F_traces;
        end
        
        % Motion/velocity data
        if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
            trial_data(t).motion = FPA.ephys.running_velocity_smooth(:);
        elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
            trial_data(t).motion = FPA.ephys.running_velocity(:);
        else
            trial_data(t).motion = zeros(size(trial_data(t).lfp));
        end
        
        % Align lengths
        n = min([length(trial_data(t).time), length(trial_data(t).lfp), ...
                 length(trial_data(t).motion)]);
        trial_data(t).time = trial_data(t).time(1:n);
        trial_data(t).lfp = trial_data(t).lfp(1:n);
        trial_data(t).motion = trial_data(t).motion(1:n);
        if ~isempty(trial_data(t).gevi) && size(trial_data(t).gevi, 1) > n
            trial_data(t).gevi = trial_data(t).gevi(1:n, :);
        end
        
        trial_data(t).num_samples = n;
        trial_data(t).duration_sec = trial_data(t).time(end) - trial_data(t).time(1);
        
        total_samples = total_samples + n;
        
        fprintf(' OK (%.1f sec, %d samples)\n', trial_data(t).duration_sec, n);
        
    catch ME
        warning('Failed to load trial %d: %s', t, ME.message);
        continue;
    end
end

% Remove empty trials
valid_trials = ~arrayfun(@isempty, {trial_data.lfp});
trial_data = trial_data(valid_trials);
num_trials = length(trial_data);

if num_trials == 0
    error('No valid trials loaded');
end

fprintf('\n  Total: %d valid trials, %d samples (%.1f sec)\n\n', ...
    num_trials, total_samples, total_samples / trial_data(1).fs);

%% ============================================================================
%  ARTIFACT DETECTION
%  ============================================================================

fprintf('Detecting artifacts...\n');
fprintf('  Methods: ');
methods_used = {};
if cfg.use_mad_detection; methods_used{end+1} = sprintf('MAD(%.1f)', cfg.mad_threshold); end
if cfg.use_variance_detection; methods_used{end+1} = sprintf('Variance(%.1f)', cfg.variance_threshold_factor); end
if cfg.use_zscore_detection; methods_used{end+1} = sprintf('Z-score(%.1f)', cfg.zscore_threshold); end
fprintf('%s\n', strjoin(methods_used, ' + '));
fprintf('  Combine: %s\n\n', cfg.combine_method);

% Process each trial
for t = 1:num_trials
    fprintf('  Trial %d: ', t);
    
    lfp = trial_data(t).lfp;
    fs = trial_data(t).fs;
    n = length(lfp);
    
    % Initialize artifact masks for each method
    artifact_mad = false(n, 1);
    artifact_var = false(n, 1);
    artifact_zscore = false(n, 1);
    
    % -------------------------------------------------------------------------
    % Method 1: MAD-based detection (robust to outliers)
    % -------------------------------------------------------------------------
    if cfg.use_mad_detection
        med = median(lfp);
        mad_val = median(abs(lfp - med));
        if mad_val > 0
            threshold = cfg.mad_threshold * mad_val;
            artifact_mad = abs(lfp - med) > threshold;
        end
    end
    
    % -------------------------------------------------------------------------
    % Method 2: Rolling variance detection
    % -------------------------------------------------------------------------
    if cfg.use_variance_detection
        win_samples = round(cfg.variance_window_sec * fs);
        if win_samples < 3
            win_samples = 3;
        end
        
        % Compute rolling variance using convolution
        lfp_sq = lfp.^2;
        kernel = ones(win_samples, 1) / win_samples;
        mean_sq = conv(lfp_sq, kernel, 'same');
        mean_val = conv(lfp, kernel, 'same');
        rolling_var = mean_sq - mean_val.^2;
        rolling_var(rolling_var < 0) = 0; % Numerical stability
        
        med_var = median(rolling_var);
        if med_var > 0
            artifact_var = rolling_var > cfg.variance_threshold_factor * med_var;
        end
    end
    
    % -------------------------------------------------------------------------
    % Method 3: Z-score detection
    % -------------------------------------------------------------------------
    if cfg.use_zscore_detection
        lfp_z = zscore(lfp);
        artifact_zscore = abs(lfp_z) > cfg.zscore_threshold;
    end
    
    % -------------------------------------------------------------------------
    % Combine methods
    % -------------------------------------------------------------------------
    switch lower(cfg.combine_method)
        case 'union'
            artifact_mask = artifact_mad | artifact_var | artifact_zscore;
        case 'intersection'
            artifact_mask = artifact_mad & artifact_var & artifact_zscore;
        otherwise
            artifact_mask = artifact_mad | artifact_var | artifact_zscore;
    end
    
    % -------------------------------------------------------------------------
    % Post-processing: Extend and filter artifact regions
    % -------------------------------------------------------------------------
    
    % Extend artifacts by specified amount
    extend_samples = round(cfg.extend_artifact_sec * fs);
    if extend_samples > 0
        artifact_extended = artifact_mask;
        for i = 1:extend_samples
            artifact_extended = artifact_extended | [false; artifact_mask(1:end-1)] | [artifact_mask(2:end); false];
        end
        artifact_mask = artifact_extended;
    end
    
    % Remove very short artifacts (likely noise spikes, not true artifacts)
    min_samples = round(cfg.min_artifact_duration_sec * fs);
    artifact_mask = remove_short_segments(artifact_mask, min_samples);
    
    % -------------------------------------------------------------------------
    % Store results
    % -------------------------------------------------------------------------
    trial_data(t).artifact_mask = artifact_mask;
    trial_data(t).artifact_pct = 100 * sum(artifact_mask) / n;
    trial_data(t).artifact_segments = find_segments(artifact_mask);
    trial_data(t).num_artifacts = size(trial_data(t).artifact_segments, 1);
    
    % Individual method contributions
    trial_data(t).artifact_mad = artifact_mad;
    trial_data(t).artifact_var = artifact_var;
    trial_data(t).artifact_zscore = artifact_zscore;
    
    % Determine recommendation
    if trial_data(t).artifact_pct > cfg.trial_exclusion_threshold * 100
        trial_data(t).recommendation = 'EXCLUDE';
    else
        trial_data(t).recommendation = 'CLEAN';
    end
    
    fprintf('%.1f%% contaminated (%d segments) -> %s\n', ...
        trial_data(t).artifact_pct, trial_data(t).num_artifacts, ...
        trial_data(t).recommendation);
end

%% ============================================================================
%  SUMMARY STATISTICS
%  ============================================================================

fprintf('\n');
fprintf('============================================================\n');
fprintf('  ARTIFACT DETECTION SUMMARY\n');
fprintf('============================================================\n');

total_artifact_samples = sum([trial_data.artifact_pct] .* [trial_data.num_samples] / 100);
total_samples = sum([trial_data.num_samples]);

fprintf('  Total contamination: %.1f%% (%d / %d samples)\n', ...
    100 * total_artifact_samples / total_samples, ...
    round(total_artifact_samples), total_samples);
fprintf('  Exclusion threshold: %.0f%%\n', cfg.trial_exclusion_threshold * 100);
fprintf('\n');

exclude_trials = find(strcmp({trial_data.recommendation}, 'EXCLUDE'));
clean_trials = find(strcmp({trial_data.recommendation}, 'CLEAN'));

fprintf('  Trials to EXCLUDE: %d / %d\n', length(exclude_trials), num_trials);
if ~isempty(exclude_trials)
    for i = 1:length(exclude_trials)
        t = exclude_trials(i);
        fprintf('    - Trial %d: %.1f%% contaminated\n', t, trial_data(t).artifact_pct);
    end
end

fprintf('  Trials to CLEAN: %d / %d\n', length(clean_trials), num_trials);
if ~isempty(clean_trials)
    for i = 1:length(clean_trials)
        t = clean_trials(i);
        fprintf('    - Trial %d: %.1f%% contaminated (%d segments to remove)\n', ...
            t, trial_data(t).artifact_pct, trial_data(t).num_artifacts);
    end
end

fprintf('============================================================\n\n');

%% ============================================================================
%  INTERACTIVE VISUALIZATION
%  ============================================================================

if cfg.show_interactive_plot
    fprintf('Generating interactive artifact review plot...\n');
    
    % Create figure
    fig = figure('Name', sprintf('LFP Artifact Review - %s %s-%s', MOUSE_NAME, RECORDING_DATE, RECORDING_ID), ...
        'NumberTitle', 'off', ...
        'Position', [50, 50, cfg.figure_width, cfg.figure_height], ...
        'Color', 'w');
    
    % Calculate layout
    num_rows = num_trials;
    row_height = 0.85 / num_rows;
    
    % Determine global y-axis limits for consistent scaling
    all_lfp = vertcat(trial_data.lfp);
    y_lim = [prctile(all_lfp, 1), prctile(all_lfp, 99)] * 1.2;
    
    % Also get motion scaling
    all_motion = vertcat(trial_data.motion);
    motion_scale = max(abs(all_motion)) / (y_lim(2) - y_lim(1)) * 0.3;
    if motion_scale == 0
        motion_scale = 1;
    end
    
    ax_handles = gobjects(num_trials, 1);
    
    for t = 1:num_trials
        % Create subplot for this trial
        ax = axes('Parent', fig, ...
            'Position', [0.06, 0.95 - t * row_height, 0.88, row_height * 0.85]);
        ax_handles(t) = ax;
        hold(ax, 'on');
        
        time = trial_data(t).time;
        lfp = trial_data(t).lfp;
        motion = trial_data(t).motion;
        artifact_mask = trial_data(t).artifact_mask;
        
        % Plot artifact regions as red shaded areas
        segments = trial_data(t).artifact_segments;
        for s = 1:size(segments, 1)
            seg_start = time(segments(s, 1));
            seg_end = time(segments(s, 2));
            fill(ax, [seg_start, seg_end, seg_end, seg_start], ...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
                cfg.artifact_color(1:3), 'FaceAlpha', cfg.artifact_color(4), ...
                'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        
        % Plot motion (scaled and shifted)
        motion_scaled = motion / motion_scale + y_lim(1) + (y_lim(2) - y_lim(1)) * 0.2;
        plot(ax, time, motion_scaled, 'Color', [1, 0.5, 0, 0.7], 'LineWidth', 0.5);
        
        % Plot LFP
        plot(ax, time, lfp, 'b', 'LineWidth', 0.5);
        
        % Configure axes
        xlim(ax, [time(1), time(end)]);
        ylim(ax, y_lim);
        
        % Labels
        if t == num_trials
            xlabel(ax, 'Time (s)', 'FontSize', 10);
        else
            set(ax, 'XTickLabel', []);
        end
        
        % Trial info on y-axis label
        ylabel_str = sprintf('Trial %d\n%.1f%% artifact', t, trial_data(t).artifact_pct);
        ylabel(ax, ylabel_str, 'FontSize', 9, 'Rotation', 0, 'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle');
        
        % Add recommendation badge
        if strcmp(trial_data(t).recommendation, 'EXCLUDE')
            text(ax, time(end), y_lim(2), ' EXCLUDE ', ...
                'Color', 'w', 'BackgroundColor', [0.8, 0.2, 0.2], ...
                'FontSize', 9, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        else
            text(ax, time(end), y_lim(2), ' CLEAN ', ...
                'Color', 'w', 'BackgroundColor', [0.2, 0.6, 0.2], ...
                'FontSize', 9, 'FontWeight', 'bold', ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top');
        end
        
        hold(ax, 'off');
        box(ax, 'on');
    end
    
    % Add title
    sgtitle(sprintf('LFP Artifact Review: %s | %s-%s | Total: %.1f%% contaminated', ...
        MOUSE_NAME, RECORDING_DATE, RECORDING_ID, ...
        100 * total_artifact_samples / total_samples), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Add legend
    legend_ax = axes('Parent', fig, 'Position', [0.06, 0.01, 0.88, 0.03], 'Visible', 'off');
    hold(legend_ax, 'on');
    h1 = plot(legend_ax, NaN, NaN, 'b', 'LineWidth', 2);
    h2 = plot(legend_ax, NaN, NaN, 'Color', [1, 0.5, 0], 'LineWidth', 2);
    h3 = fill(legend_ax, NaN, NaN, cfg.artifact_color(1:3), 'FaceAlpha', cfg.artifact_color(4), 'EdgeColor', 'none');
    legend(legend_ax, [h1, h2, h3], {'LFP', 'Motion', 'Artifact'}, ...
        'Orientation', 'horizontal', 'Location', 'south', 'FontSize', 10);
    hold(legend_ax, 'off');
    
    drawnow;
    fprintf('  Plot generated.\n\n');
end

%% ============================================================================
%  USER APPROVAL
%  ============================================================================

fprintf('============================================================\n');
fprintf('  USER APPROVAL\n');
fprintf('============================================================\n');
fprintf('  Review the plot and choose an action:\n');
fprintf('    [1] ACCEPT algorithm recommendations\n');
fprintf('    [2] MANUAL selection - exclude specific trials\n');
fprintf('    [3] CANCEL - exit without saving\n');
fprintf('\n');

user_choice = input('  Enter choice (1/2/3): ', 's');

switch user_choice
    case '1'
        % Accept algorithm recommendations
        fprintf('\n  Accepted algorithm recommendations.\n');
        user_decision = 'algorithm';
        excluded_trials = exclude_trials;
        
    case '2'
        % Manual selection
        fprintf('\n  Enter trial numbers to EXCLUDE (comma-separated, e.g., "1,3,5"):\n');
        fprintf('  Current recommendations: ');
        if isempty(exclude_trials)
            fprintf('None\n');
        else
            fprintf('%s\n', strjoin(arrayfun(@num2str, exclude_trials, 'UniformOutput', false), ', '));
        end
        
        manual_input = input('  Trials to exclude: ', 's');
        if isempty(manual_input)
            excluded_trials = [];
        else
            excluded_trials = str2num(manual_input); %#ok<ST2NM>
            excluded_trials = excluded_trials(excluded_trials >= 1 & excluded_trials <= num_trials);
        end
        
        fprintf('  Manual selection: excluding trials [%s]\n', ...
            strjoin(arrayfun(@num2str, excluded_trials, 'UniformOutput', false), ', '));
        user_decision = 'manual';
        
    case '3'
        fprintf('\n  Cancelled. No changes saved.\n');
        return;
        
    otherwise
        fprintf('\n  Invalid choice. Exiting.\n');
        return;
end

% Update recommendations based on user decision
for t = 1:num_trials
    if ismember(t, excluded_trials)
        trial_data(t).final_decision = 'EXCLUDE';
    else
        trial_data(t).final_decision = 'CLEAN';
    end
end

%% ============================================================================
%  SAVE ARTIFACT INFORMATION
%  ============================================================================

if cfg.save_artifact_info
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('  SAVING ARTIFACT INFORMATION\n');
    fprintf('============================================================\n');
    
    % Create artifact info structure
    ArtifactInfo = struct();
    ArtifactInfo.mouse_name = MOUSE_NAME;
    ArtifactInfo.recording_date = RECORDING_DATE;
    ArtifactInfo.recording_id = RECORDING_ID;
    ArtifactInfo.analysis_date = datestr(now);
    ArtifactInfo.config = cfg;
    
    % Summary
    ArtifactInfo.summary.total_trials = num_trials;
    ArtifactInfo.summary.total_samples = total_samples;
    ArtifactInfo.summary.total_artifact_pct = 100 * total_artifact_samples / total_samples;
    ArtifactInfo.summary.excluded_trials = excluded_trials;
    ArtifactInfo.summary.user_decision = user_decision;
    
    % Per-trial information
    for t = 1:num_trials
        ArtifactInfo.trials(t).trial_number = t;
        ArtifactInfo.trials(t).filename = trial_data(t).filename;
        ArtifactInfo.trials(t).num_samples = trial_data(t).num_samples;
        ArtifactInfo.trials(t).duration_sec = trial_data(t).duration_sec;
        ArtifactInfo.trials(t).artifact_pct = trial_data(t).artifact_pct;
        ArtifactInfo.trials(t).num_artifacts = trial_data(t).num_artifacts;
        ArtifactInfo.trials(t).artifact_segments = trial_data(t).artifact_segments;
        ArtifactInfo.trials(t).artifact_mask = trial_data(t).artifact_mask;
        ArtifactInfo.trials(t).algorithm_recommendation = trial_data(t).recommendation;
        ArtifactInfo.trials(t).final_decision = trial_data(t).final_decision;
    end
    
    % Generate output filename
    output_filename = sprintf('%s-%s-%s%s.mat', ...
        MOUSE_NAME, RECORDING_DATE, RECORDING_ID, cfg.output_suffix);
    output_path = fullfile(DATA_PATH, output_filename);
    
    save(output_path, 'ArtifactInfo', '-v7.3');
    fprintf('  Artifact info saved to:\n  %s\n', output_path);
    
    % Also save figure
    if cfg.show_interactive_plot && isvalid(fig)
        fig_filename = sprintf('%s-%s-%s%s.png', ...
            MOUSE_NAME, RECORDING_DATE, RECORDING_ID, cfg.output_suffix);
        fig_path = fullfile(DATA_PATH, fig_filename);
        saveas(fig, fig_path);
        fprintf('  Figure saved to:\n  %s\n', fig_path);
    end
end

%% ============================================================================
%  OPTIONAL: SAVE CLEANED DATA STRUCTS
%  ============================================================================

if cfg.save_cleaned_data
    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('  SAVING CLEANED DATA STRUCTS\n');
    fprintf('============================================================\n');
    
    for t = 1:num_trials
        if strcmp(trial_data(t).final_decision, 'EXCLUDE')
            fprintf('  Trial %d: EXCLUDED (no cleaned file saved)\n', t);
            continue;
        end
        
        fprintf('  Trial %d: Cleaning and saving...', t);
        
        % Load original file
        loaded = load(trial_data(t).filepath);
        FPA = loaded.FiberPhotometryAnalysis;
        
        % Get clean indices (non-artifact samples)
        clean_mask = ~trial_data(t).artifact_mask;
        clean_indices = find(clean_mask);
        
        if isempty(clean_indices)
            fprintf(' No clean samples, skipping.\n');
            continue;
        end
        
        % Trim all relevant fields
        % Time vector
        FPA.time.time_vector_seconds = FPA.time.time_vector_seconds(clean_indices);
        
        % Signals
        if isfield(FPA, 'signals')
            signal_fields = fieldnames(FPA.signals);
            for f = 1:length(signal_fields)
                field = signal_fields{f};
                if isnumeric(FPA.signals.(field)) && size(FPA.signals.(field), 1) == trial_data(t).num_samples
                    FPA.signals.(field) = FPA.signals.(field)(clean_indices, :);
                end
            end
        end
        
        % Ephys
        if isfield(FPA, 'ephys')
            ephys_fields = {'lfp_raw_aligned_HP', 'lfp_raw_aligned_mPFC', 'lfp_z_HP', ...
                           'running_velocity', 'running_velocity_smooth'};
            for f = 1:length(ephys_fields)
                field = ephys_fields{f};
                if isfield(FPA.ephys, field) && isnumeric(FPA.ephys.(field)) && ...
                   length(FPA.ephys.(field)) == trial_data(t).num_samples
                    FPA.ephys.(field) = FPA.ephys.(field)(clean_indices);
                end
            end
        end
        
        % Update parameters
        FPA.parameters.num_frames = length(clean_indices);
        FPA.parameters.recording_duration_sec = length(clean_indices) / trial_data(t).fs;
        
        % Add artifact removal metadata
        FPA.artifact_removal.applied = true;
        FPA.artifact_removal.date = datestr(now);
        FPA.artifact_removal.original_samples = trial_data(t).num_samples;
        FPA.artifact_removal.clean_samples = length(clean_indices);
        FPA.artifact_removal.removed_pct = trial_data(t).artifact_pct;
        FPA.artifact_removal.config = cfg;
        
        % Save with _cleaned suffix
        FiberPhotometryAnalysis = FPA;
        [~, base_name, ~] = fileparts(trial_data(t).filename);
        cleaned_filename = [base_name '_cleaned.mat'];
        cleaned_path = fullfile(DATA_PATH, cleaned_filename);
        
        save(cleaned_path, 'FiberPhotometryAnalysis', '-v7.3');
        fprintf(' Saved (%d -> %d samples)\n', trial_data(t).num_samples, length(clean_indices));
    end
end

fprintf('\n');
fprintf('============================================================\n');
fprintf('  ARTIFACT REMOVAL COMPLETE\n');
fprintf('============================================================\n\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================

function segments = find_segments(mask)
%FIND_SEGMENTS Find contiguous true segments in a logical mask
%  Returns Nx2 matrix where each row is [start_idx, end_idx]

segments = [];
in_segment = false;
start_idx = 0;

for i = 1:length(mask)
    if mask(i) && ~in_segment
        % Start of new segment
        in_segment = true;
        start_idx = i;
    elseif ~mask(i) && in_segment
        % End of segment
        in_segment = false;
        segments = [segments; start_idx, i-1]; %#ok<AGROW>
    end
end

% Handle segment that extends to end
if in_segment
    segments = [segments; start_idx, length(mask)];
end
end

function mask = remove_short_segments(mask, min_samples)
%REMOVE_SHORT_SEGMENTS Remove segments shorter than min_samples

segments = find_segments(mask);

for s = 1:size(segments, 1)
    seg_length = segments(s, 2) - segments(s, 1) + 1;
    if seg_length < min_samples
        mask(segments(s, 1):segments(s, 2)) = false;
    end
end
end
