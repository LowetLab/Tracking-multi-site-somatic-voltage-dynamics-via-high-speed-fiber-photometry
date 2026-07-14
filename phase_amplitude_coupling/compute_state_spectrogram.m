%% ============================================================================
%  RUNNING STATE PHASE-ALIGNED SPECTROGRAM ANALYSIS
%  ============================================================================
%  README - 配置参数与数据结构说明
%  ============================================================================
%
%  DESCRIPTION:
%  This script computes phase-aligned spectrograms separately for running and
%  non-running epochs. It reads pre-computed epoch data from all trials of a
%  single animal and averages epochs based on locomotion velocity thresholds.
%
%  USAGE:
%  1. Run the script
%  2. Select an animal folder (e.g., Process/phase_aligned_spectrogram/.../Animal01/)
%  3. Script will scan all trials, categorize epochs, and compute averages
%  4. Results saved to animal_folder/run_nonrun/LFP|Fiber1|Fiber2/running|non_running/
%
%  CONFIGURATION PARAMETERS:
%  -------------------------------------------------------------------------
%  Parameter                Default     Description
%  ----------------------   ---------   ------------------------------------
%  RUNNING_THRESHOLD        2.0         Speed threshold (cm/s) for running classification
%                                       (epoch is running if ALL points > 2 cm/s)
%  NON_RUNNING_THRESHOLD    0.1         Speed threshold (cm/s) for non-running classification
%                                       (epoch is rest if ALL points < 0.1 cm/s)
%  MIN_TIME_FRACTION        1.0         Fraction of epoch that must satisfy speed criterion
%                                       (1.0 = ALL timepoints must be continuously above/below threshold)
%  (Missing or empty/all-NaN velocity → error in classification; no mean_velocity fallback.)
%
%  OUTPUT FILES (per signal type and running state):
%  -------------------------------------------------------------------------
%  - PhaseAlignedSpectrogram_[running|non_running].mat (Python-friendly)
%  - PhaseAlignedSpectrogram_[running|non_running]_full.mat (complete results)
%  - PhaseAlignedSpectrogram_[running|non_running].png
%  - PhaseAlignedSpectrogram_[running|non_running].fig
%
%  PYTHON DATA STRUCTURE (PhaseAlignedSpectrogram_[state].mat):
%  -------------------------------------------------------------------------
%  phase_bins              [1 x 36]        Phase bin centers (radians)
%  phase_bins_deg          [1 x 36]        Phase bin centers (degrees)
%  freq_vector             [1 x N_freq]    Frequency vector (Hz)
%  mean_spectrogram        [N_freq x 36]   Averaged 1-cycle spectrogram
%  cycle1_mean_spectrogram [N_freq x 36]   Averaged cycle 1 spectrogram
%  cycle2_mean_spectrogram [N_freq x 36]   Averaged cycle 2 spectrogram
%  two_cycle_spectrogram   [N_freq x 72]   Concatenated 2-cycle view
%  n_epochs                scalar          Number of epochs used
%  mean_velocity           scalar          Mean velocity across all used epochs
%  running_threshold       scalar          Threshold used for running
%  non_running_threshold   scalar          Threshold used for non-running
%  n_trials                scalar          Number of trials contributing data
%  animal_name             string          Animal identifier
%  signal_type             string          'LFP' or 'Fiber1', 'Fiber2', etc.
%  running_state           string          'running' or 'non_running'
%
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  CONFIGURATION PARAMETERS
%  ============================================================================

% Velocity thresholds for epoch classification
% Running: ALL points must be > 2 cm/s
% Rest: ALL points must be < 0.1 cm/s
RUNNING_THRESHOLD = 2.0;     % cm/s - speed above this = running
NON_RUNNING_THRESHOLD = 0.1; % cm/s - speed below this = non-running (rest)

% Time fraction requirement (ensures consistent behavioral state within epoch)
% An epoch is classified as "running" only if ≥ this fraction of timepoints
% have speed > RUNNING_THRESHOLD (same logic for non-running)
MIN_TIME_FRACTION = 1.0;     % 100% = ALL timepoints must satisfy the speed criterion (absolute continuous)

% Artifact rejection
% Set EXCLUDE_ARTIFACT_EPOCHS = true to skip any epoch that contains at least
% one flagged sample (has_artifact = true in the results file).
% Epochs from trials without an artifact mask are always kept.
EXCLUDE_ARTIFACT_EPOCHS = true;

% Figure settings
SAVE_FIGURES = true;

%% ============================================================================
%  SELECT ANIMAL FOLDER
%  ============================================================================

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);

% Default starting path
default_path = fullfile(project_root, 'Process', 'phase_aligned_spectrogram');
if ~exist(default_path, 'dir')
    default_path = project_root;
end

% Select animal folder
animal_folder = uigetdir(default_path, 'Select Animal Folder (e.g., Animal01)');

if isequal(animal_folder, 0)
    disp('Folder selection cancelled.');
    return;
end

% Extract animal name from folder path
[~, animal_name] = fileparts(animal_folder);
fprintf('\n=== Processing Animal: %s ===\n', animal_name);
fprintf('Folder: %s\n', animal_folder);

lfp_ephys_field = '';

%% ============================================================================
%  FIND ALL RESULT FILES
%  ============================================================================

fprintf('\nScanning for PhaseAlignedSpectrogram result files...\n');

% Find all result .mat files (not ForPython versions)
result_files = dir(fullfile(animal_folder, '**', 'PhaseAlignedSpectrogram_Results_*.mat'));

if isempty(result_files)
    error('No PhaseAlignedSpectrogram_Results_*.mat files found in: %s', animal_folder);
end

fprintf('Found %d trial result files.\n', length(result_files));

%% ============================================================================
%  COLLECT ALL EPOCHS FROM ALL TRIALS
%  ============================================================================

fprintf('\nCollecting epochs from all trials...\n');

% Initialize storage structures
LFP_epochs = struct('spectrograms', [], 'velocities', [], 'phase_bins', [], ...
                    'freq_vector', [], 'cycle1_specs', [], 'cycle2_specs', []);
Fiber_epochs = {};  % Cell array for multiple fibers

n_trials_with_lfp = 0;
n_trials_with_fiber = 0;
max_fiber_idx = 0;

for file_idx = 1:length(result_files)
    file_path = fullfile(result_files(file_idx).folder, result_files(file_idx).name);
    fprintf('  Loading [%d/%d]: %s\n', file_idx, length(result_files), result_files(file_idx).name);
    
    try
        data = load(file_path);
        results = data.results;
        
        % Get metadata (use correct field names)
        if isempty(LFP_epochs.phase_bins) && isfield(results, 'metadata')
            if isfield(results.metadata, 'phase_bin_centers')
                LFP_epochs.phase_bins = results.metadata.phase_bin_centers;
            elseif isfield(results.metadata, 'phase_bins')
                LFP_epochs.phase_bins = results.metadata.phase_bins;
            end
            LFP_epochs.freq_vector = results.metadata.freq_vector;
            if isempty(lfp_ephys_field) && isfield(results.metadata, 'lfp_ephys_field')
                lfp_ephys_field = char(string(results.metadata.lfp_ephys_field));
            end
        end
        
        % Process LFP epochs
        if isfield(results, 'LFP') && isfield(results.LFP, 'all_epochs_raw')
            n_trials_with_lfp = n_trials_with_lfp + 1;
            lfp_data = results.LFP;
            n_epochs = length(lfp_data.all_epochs_raw);
            
            % Get phase_bins and freq_vector from LFP data if not yet set
            if isempty(LFP_epochs.phase_bins) && isfield(lfp_data, 'phase_bins')
                LFP_epochs.phase_bins = lfp_data.phase_bins;
            end
            if isempty(LFP_epochs.freq_vector) && isfield(lfp_data, 'freq_vector')
                LFP_epochs.freq_vector = lfp_data.freq_vector;
            end
            
            for ei = 1:n_epochs
                epoch = lfp_data.all_epochs_raw{ei};

                % Artifact rejection
                if EXCLUDE_ARTIFACT_EPOCHS && isfield(epoch, 'has_artifact') && epoch.has_artifact
                    continue;
                end

                velocity = epoch.mean_velocity;
                
                % Store epoch data
                LFP_epochs.velocities(end+1) = velocity;
                
                % Store raw epoch data for later processing
                if isempty(LFP_epochs.spectrograms)
                    LFP_epochs.spectrograms = {epoch};
                else
                    LFP_epochs.spectrograms{end+1} = epoch;
                end
            end
        end
        
        % Process Fiber epochs
        fiber_fields = fieldnames(results);
        fiber_fields = fiber_fields(startsWith(fiber_fields, 'Fiber'));
        
        for fi = 1:length(fiber_fields)
            field_name = fiber_fields{fi};
            fiber_idx = str2double(field_name(6:end));  % Extract number from 'FiberN'
            max_fiber_idx = max(max_fiber_idx, fiber_idx);
            
            % Initialize fiber storage if needed
            while length(Fiber_epochs) < fiber_idx
                Fiber_epochs{end+1} = struct('spectrograms', [], 'velocities', [], ...
                                             'phase_bins', [], 'freq_vector', [], ...
                                             'roi_info', '');
            end
            
            fiber_data = results.(field_name);
            
            if isfield(fiber_data, 'all_epochs_raw')
                if fi == 1
                    n_trials_with_fiber = n_trials_with_fiber + 1;
                end
                
                % Store metadata
                if isempty(Fiber_epochs{fiber_idx}.phase_bins)
                    Fiber_epochs{fiber_idx}.phase_bins = fiber_data.phase_bins;
                    Fiber_epochs{fiber_idx}.freq_vector = fiber_data.freq_vector;
                    if isfield(fiber_data, 'roi_info')
                        Fiber_epochs{fiber_idx}.roi_info = fiber_data.roi_info;
                    end
                end
                
                n_epochs = length(fiber_data.all_epochs_raw);
                for ei = 1:n_epochs
                    epoch = fiber_data.all_epochs_raw{ei};

                    % Artifact rejection
                    if EXCLUDE_ARTIFACT_EPOCHS && isfield(epoch, 'has_artifact') && epoch.has_artifact
                        continue;
                    end

                    velocity = epoch.mean_velocity;
                    
                    Fiber_epochs{fiber_idx}.velocities(end+1) = velocity;
                    
                    if isempty(Fiber_epochs{fiber_idx}.spectrograms)
                        Fiber_epochs{fiber_idx}.spectrograms = {epoch};
                    else
                        Fiber_epochs{fiber_idx}.spectrograms{end+1} = epoch;
                    end
                end
            end
        end
        
    catch ME
        warning('Error loading file %s: %s', result_files(file_idx).name, ME.message);
    end
end

fprintf('\nData collection complete:\n');
fprintf('  Trials with LFP: %d\n', n_trials_with_lfp);
fprintf('  Trials with Fiber: %d\n', n_trials_with_fiber);
fprintf('  Artifact epoch exclusion: %s\n', mat2str(EXCLUDE_ARTIFACT_EPOCHS));
fprintf('  Total LFP epochs kept: %d\n', length(LFP_epochs.velocities));
for fi = 1:length(Fiber_epochs)
    fprintf('  Total Fiber%d epochs kept: %d\n', fi, length(Fiber_epochs{fi}.velocities));
end

%% ============================================================================
%  CREATE OUTPUT DIRECTORY
%  ============================================================================

output_base = fullfile(animal_folder, 'run_nonrun');
if ~exist(output_base, 'dir')
    mkdir(output_base);
end
fprintf('\nOutput directory: %s\n', output_base);

%% ============================================================================
%  HELPER FUNCTION: CLASSIFY EPOCHS BY TIME FRACTION
%  ============================================================================

    function [is_running, is_non_running, running_fraction, non_running_fraction] = ...
            classify_epochs_by_time_fraction(epochs, run_thresh, non_run_thresh, min_fraction)
        % Classify epochs based on fraction of time spent above/below speed thresholds
        %
        % Inputs:
        %   epochs: cell array of epoch structures (each with .velocity field)
        %   run_thresh: speed threshold for running (cm/s)
        %   non_run_thresh: speed threshold for non-running (cm/s)
        %   min_fraction: minimum fraction of time that must satisfy criterion (e.g., 0.8)
        %
        % Outputs:
        %   is_running: logical array, true if epoch qualifies as running
        %   is_non_running: logical array, true if epoch qualifies as non-running
        %   running_fraction: fraction of time with speed > run_thresh for each epoch
        %   non_running_fraction: fraction of time with speed < non_run_thresh for each epoch
        
        n_epochs = length(epochs);
        is_running = false(1, n_epochs);
        is_non_running = false(1, n_epochs);
        running_fraction = zeros(1, n_epochs);
        non_running_fraction = zeros(1, n_epochs);
        
        for ei = 1:n_epochs
            epoch = epochs{ei};
            
            % Per-sample velocity required — error if missing or all NaN (no mean_velocity fallback)
            if ~isfield(epoch, 'velocity') || isempty(epoch.velocity)
                error('compute_state_spectrogram:MissingVelocity', ...
                    'Epoch %d: missing velocity vector (required for running/rest classification).', ei);
            end
            vel = epoch.velocity(:);
            vel = vel(~isnan(vel));
            n_points = length(vel);
            if n_points == 0
                error('compute_state_spectrogram:EmptyVelocity', ...
                    'Epoch %d: velocity is empty or all NaN (required for classification).', ei);
            end
            
            running_fraction(ei) = sum(vel > run_thresh) / n_points;
            non_running_fraction(ei) = sum(vel < non_run_thresh) / n_points;
            is_running(ei) = running_fraction(ei) >= min_fraction;
            is_non_running(ei) = non_running_fraction(ei) >= min_fraction;
        end
    end

%% ============================================================================
%  HELPER FUNCTION: COMPUTE AVERAGE SPECTROGRAM FROM SELECTED EPOCHS
%  ============================================================================

    function [mean_spec, cycle1_spec, cycle2_spec, n_used] = compute_average_from_selected_epochs(epochs, selected_idx, phase_bins, n_bins)
        % Compute average phase-aligned spectrogram from selected epochs
        %
        % Inputs:
        %   epochs: cell array of epoch structures
        %   selected_idx: logical or index array of epochs to include
        %   phase_bins: phase bin edges
        %   n_bins: number of phase bins
        %
        % Outputs:
        %   mean_spec: averaged spectrogram across all cycles
        %   cycle1_spec: averaged spectrogram for cycle 1
        %   cycle2_spec: averaged spectrogram for cycle 2
        %   n_used: number of epochs used
        
        % Convert logical to indices if needed
        if islogical(selected_idx)
            valid_idx = find(selected_idx);
        else
            valid_idx = selected_idx;
        end
        n_used = length(valid_idx);
        
        if n_used == 0
            mean_spec = [];
            cycle1_spec = [];
            cycle2_spec = [];
            return;
        end
        
        % Get dimensions from first valid epoch
        first_epoch = epochs{valid_idx(1)};
        n_freqs = size(first_epoch.spectrogram, 1);
        
        % Initialize accumulators
        bin_power_sum = zeros(n_freqs, n_bins);
        bin_count = zeros(1, n_bins);
        cycle1_power_sum = zeros(n_freqs, n_bins);
        cycle1_count = zeros(1, n_bins);
        cycle2_power_sum = zeros(n_freqs, n_bins);
        cycle2_count = zeros(1, n_bins);
        
        % Process each valid epoch
        for idx = 1:n_used
            ei = valid_idx(idx);
            epoch = epochs{ei};
            spec = epoch.spectrogram;
            phase = epoch.phase;
            
            n_timepoints = length(phase);
            
            % Find cycle boundary (where phase wraps from pi to -pi)
            phase_diff = diff(phase);
            cycle_boundary_idx = find(phase_diff < -pi, 1, 'first');
            
            if isempty(cycle_boundary_idx)
                cycle_boundary_idx = floor(n_timepoints / 2);
            end
            
            % Process all timepoints
            for ti = 1:n_timepoints
                ph = phase(ti);
                
                % Find which bin this phase belongs to
                bin_idx = find(ph >= phase_bins(1:end-1) & ph < phase_bins(2:end), 1);
                if isempty(bin_idx)
                    if ph >= phase_bins(end)
                        bin_idx = n_bins;
                    else
                        bin_idx = 1;
                    end
                end
                
                % Accumulate power
                bin_power_sum(:, bin_idx) = bin_power_sum(:, bin_idx) + spec(:, ti);
                bin_count(bin_idx) = bin_count(bin_idx) + 1;
                
                % Separate cycle 1 and cycle 2
                if ti <= cycle_boundary_idx
                    cycle1_power_sum(:, bin_idx) = cycle1_power_sum(:, bin_idx) + spec(:, ti);
                    cycle1_count(bin_idx) = cycle1_count(bin_idx) + 1;
                else
                    cycle2_power_sum(:, bin_idx) = cycle2_power_sum(:, bin_idx) + spec(:, ti);
                    cycle2_count(bin_idx) = cycle2_count(bin_idx) + 1;
                end
            end
        end
        
        % Compute means
        mean_spec = zeros(n_freqs, n_bins);
        cycle1_spec = zeros(n_freqs, n_bins);
        cycle2_spec = zeros(n_freqs, n_bins);
        
        for bi = 1:n_bins
            if bin_count(bi) > 0
                mean_spec(:, bi) = bin_power_sum(:, bi) / bin_count(bi);
            end
            if cycle1_count(bi) > 0
                cycle1_spec(:, bi) = cycle1_power_sum(:, bi) / cycle1_count(bi);
            end
            if cycle2_count(bi) > 0
                cycle2_spec(:, bi) = cycle2_power_sum(:, bi) / cycle2_count(bi);
            end
        end
    end

%% ============================================================================
%  HELPER FUNCTION: SAVE RESULTS AND GENERATE FIGURES
%  ============================================================================

    function save_results_and_figure(output_dir, signal_type, running_state, ...
            mean_spec, cycle1_spec, cycle2_spec, phase_bins, freq_vector, ...
            n_epochs, mean_velocity, n_trials, animal_name, ...
            running_thresh, non_running_thresh, lfp_ephys_field_arg)
        
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        % Phase bin centers
        phase_bin_centers = phase_bins(1:end-1) + diff(phase_bins)/2;
        
        % Create Python-friendly results
        python_results = struct();
        python_results.phase_bins = phase_bin_centers;
        python_results.phase_bins_deg = rad2deg(phase_bin_centers);
        python_results.freq_vector = freq_vector;
        python_results.mean_spectrogram = mean_spec;
        python_results.cycle1_mean_spectrogram = cycle1_spec;
        python_results.cycle2_mean_spectrogram = cycle2_spec;
        python_results.two_cycle_spectrogram = [cycle1_spec, cycle2_spec];
        python_results.n_epochs = n_epochs;
        python_results.mean_velocity = mean_velocity;
        python_results.running_threshold = running_thresh;
        python_results.non_running_threshold = non_running_thresh;
        python_results.n_trials = n_trials;
        python_results.animal_name = animal_name;
        python_results.signal_type = signal_type;
        python_results.running_state = running_state;
        if strcmp(signal_type, 'LFP')
            python_results.lfp_ephys_field = lfp_ephys_field_arg;
        end
        
        % Save Python-friendly version
        python_filename = sprintf('PhaseAlignedSpectrogram_%s.mat', running_state);
        save(fullfile(output_dir, python_filename), '-struct', 'python_results', '-v7');
        fprintf('    Saved: %s\n', python_filename);
        
        % Save full results
        full_results = python_results;
        full_results.phase_bin_edges = phase_bins;
        full_results.timestamp = datetime('now');
        
        full_filename = sprintf('PhaseAlignedSpectrogram_%s_full.mat', running_state);
        save(fullfile(output_dir, full_filename), '-struct', 'full_results', '-v7.3');
        fprintf('    Saved: %s\n', full_filename);
        
        % Generate figure
        fig = figure('Position', [100 100 1200 500], 'Color', 'w', 'Visible', 'off');
        
        % 1-cycle view
        subplot(1, 2, 1);
        imagesc(rad2deg(phase_bin_centers), freq_vector, mean_spec);
        axis xy;
        colorbar;
        colormap(gca, 'jet');
        xlabel('Phase (degrees)');
        ylabel('Frequency (Hz)');
        title(sprintf('1-Cycle Phase-Aligned Spectrogram\n%s - %s (%s)\nEpochs: %d | Mean Velocity: %.2f', ...
            animal_name, signal_type, strrep(running_state, '_', '-'), n_epochs, mean_velocity));
        
        % 2-cycle view
        subplot(1, 2, 2);
        two_cycle_spec = [cycle1_spec, cycle2_spec];
        two_cycle_phase = [rad2deg(phase_bin_centers), rad2deg(phase_bin_centers) + 360];
        imagesc(two_cycle_phase, freq_vector, two_cycle_spec);
        axis xy;
        colorbar;
        colormap(gca, 'jet');
        xlabel('Phase (degrees)');
        ylabel('Frequency (Hz)');
        title(sprintf('2-Cycle Phase-Aligned Spectrogram\nCycle 1 (-180° to 180°) | Cycle 2 (180° to 540°)'));
        
        % Add vertical line at cycle boundary
        hold on;
        xline(180, 'w--', 'LineWidth', 1.5);
        hold off;
        
        % Save figure
        fig_filename = sprintf('PhaseAlignedSpectrogram_%s', running_state);
        saveas(fig, fullfile(output_dir, [fig_filename, '.png']));
        saveas(fig, fullfile(output_dir, [fig_filename, '.fig']));
        fprintf('    Saved: %s.png/fig\n', fig_filename);
        
        close(fig);
    end

%% ============================================================================
%  PROCESS LFP
%  ============================================================================

% Variables to store classification results for summary
LFP_classification = struct('n_total', 0, 'n_running', 0, 'n_non_running', 0, 'n_intermediate', 0);
Fiber_classification = {};

if ~isempty(LFP_epochs.spectrograms)
    fprintf('\n=== Processing LFP ===\n');
    
    % Create phase bin edges
    n_bins = length(LFP_epochs.phase_bins);
    phase_bin_edges = linspace(-pi, pi, n_bins + 1);
    
    % Classify epochs using time fraction criterion
    fprintf('  Classifying epochs (%.0f%% time fraction criterion)...\n', MIN_TIME_FRACTION * 100);
    [is_running, is_non_running, run_frac, ~] = classify_epochs_by_time_fraction(...
        LFP_epochs.spectrograms, RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, MIN_TIME_FRACTION);
    
    n_running = sum(is_running);
    n_non_running = sum(is_non_running);
    n_total = length(is_running);
    n_intermediate = n_total - n_running - n_non_running;
    
    % Save for summary
    LFP_classification.n_total = n_total;
    LFP_classification.n_running = n_running;
    LFP_classification.n_non_running = n_non_running;
    LFP_classification.n_intermediate = n_intermediate;
    
    fprintf('    Running epochs: %d (%.1f%% of total)\n', n_running, 100*n_running/n_total);
    fprintf('    Non-running epochs: %d (%.1f%% of total)\n', n_non_running, 100*n_non_running/n_total);
    fprintf('    Intermediate (excluded): %d (%.1f%% of total)\n', n_intermediate, 100*n_intermediate/n_total);
    
    % Process RUNNING epochs
    fprintf('  Computing running spectrogram...\n');
    if n_running > 0
        [mean_spec_run, cycle1_run, cycle2_run, ~] = compute_average_from_selected_epochs(...
            LFP_epochs.spectrograms, is_running, phase_bin_edges, n_bins);
        
        % Calculate mean velocity for running epochs
        mean_vel_running = mean(LFP_epochs.velocities(is_running));
        
        output_dir = fullfile(output_base, 'LFP', 'running');
        save_results_and_figure(output_dir, 'LFP', 'running', ...
            mean_spec_run, cycle1_run, cycle2_run, ...
            phase_bin_edges, LFP_epochs.freq_vector, ...
            n_running, mean_vel_running, n_trials_with_lfp, animal_name, ...
            RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, lfp_ephys_field);
    else
        fprintf('    No running epochs found.\n');
    end
    
    % Process NON-RUNNING epochs
    fprintf('  Computing non-running spectrogram...\n');
    if n_non_running > 0
        [mean_spec_nonrun, cycle1_nonrun, cycle2_nonrun, ~] = compute_average_from_selected_epochs(...
            LFP_epochs.spectrograms, is_non_running, phase_bin_edges, n_bins);
        
        % Calculate mean velocity for non-running epochs
        mean_vel_non_running = mean(LFP_epochs.velocities(is_non_running));
        
        output_dir = fullfile(output_base, 'LFP', 'non_running');
        save_results_and_figure(output_dir, 'LFP', 'non_running', ...
            mean_spec_nonrun, cycle1_nonrun, cycle2_nonrun, ...
            phase_bin_edges, LFP_epochs.freq_vector, ...
            n_non_running, mean_vel_non_running, n_trials_with_lfp, animal_name, ...
            RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, lfp_ephys_field);
    else
        fprintf('    No non-running epochs found.\n');
    end
    
    % Process ALL epochs (no classification)
    fprintf('  Computing all-epochs spectrogram...\n');
    all_idx = true(1, n_total);
    [mean_spec_all, cycle1_all, cycle2_all, ~] = compute_average_from_selected_epochs(...
        LFP_epochs.spectrograms, all_idx, phase_bin_edges, n_bins);
    
    mean_vel_all = mean(LFP_epochs.velocities);
    
    output_dir = fullfile(output_base, 'LFP', 'all_epochs');
    save_results_and_figure(output_dir, 'LFP', 'all_epochs', ...
        mean_spec_all, cycle1_all, cycle2_all, ...
        phase_bin_edges, LFP_epochs.freq_vector, ...
        n_total, mean_vel_all, n_trials_with_lfp, animal_name, ...
        RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, lfp_ephys_field);
else
    fprintf('\nNo LFP data found.\n');
end

%% ============================================================================
%  PROCESS FIBERS
%  ============================================================================

for fiber_idx = 1:length(Fiber_epochs)
    fiber_data = Fiber_epochs{fiber_idx};
    
    if isempty(fiber_data.spectrograms)
        continue;
    end
    
    signal_type = sprintf('Fiber%d', fiber_idx);
    if ~isempty(fiber_data.roi_info)
        signal_type = sprintf('Fiber%d_%s', fiber_idx, fiber_data.roi_info);
    end
    
    fprintf('\n=== Processing %s ===\n', signal_type);
    
    % Create phase bin edges
    n_bins = length(fiber_data.phase_bins);
    phase_bin_edges = linspace(-pi, pi, n_bins + 1);
    
    % Classify epochs using time fraction criterion
    fprintf('  Classifying epochs (%.0f%% time fraction criterion)...\n', MIN_TIME_FRACTION * 100);
    [is_running, is_non_running, ~, ~] = classify_epochs_by_time_fraction(...
        fiber_data.spectrograms, RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, MIN_TIME_FRACTION);
    
    n_running = sum(is_running);
    n_non_running = sum(is_non_running);
    n_total = length(is_running);
    n_intermediate = n_total - n_running - n_non_running;
    
    % Save for summary
    Fiber_classification{fiber_idx} = struct('n_total', n_total, 'n_running', n_running, ...
        'n_non_running', n_non_running, 'n_intermediate', n_intermediate);
    
    fprintf('    Running epochs: %d, Non-running epochs: %d, Intermediate: %d\n', n_running, n_non_running, n_intermediate);
    
    % Process RUNNING epochs
    fprintf('  Computing running spectrogram...\n');
    if n_running > 0
        [mean_spec_run, cycle1_run, cycle2_run, ~] = compute_average_from_selected_epochs(...
            fiber_data.spectrograms, is_running, phase_bin_edges, n_bins);
        
        mean_vel_running = mean(fiber_data.velocities(is_running));
        
        % Use simple folder name
        folder_name = sprintf('Fiber%d', fiber_idx);
        output_dir = fullfile(output_base, folder_name, 'running');
        save_results_and_figure(output_dir, signal_type, 'running', ...
            mean_spec_run, cycle1_run, cycle2_run, ...
            phase_bin_edges, fiber_data.freq_vector, ...
            n_running, mean_vel_running, n_trials_with_fiber, animal_name, ...
            RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, '');
    else
        fprintf('    No running epochs found.\n');
    end
    
    % Process NON-RUNNING epochs
    fprintf('  Computing non-running spectrogram...\n');
    if n_non_running > 0
        [mean_spec_nonrun, cycle1_nonrun, cycle2_nonrun, ~] = compute_average_from_selected_epochs(...
            fiber_data.spectrograms, is_non_running, phase_bin_edges, n_bins);
        
        mean_vel_non_running = mean(fiber_data.velocities(is_non_running));
        
        folder_name = sprintf('Fiber%d', fiber_idx);
        output_dir = fullfile(output_base, folder_name, 'non_running');
        save_results_and_figure(output_dir, signal_type, 'non_running', ...
            mean_spec_nonrun, cycle1_nonrun, cycle2_nonrun, ...
            phase_bin_edges, fiber_data.freq_vector, ...
            n_non_running, mean_vel_non_running, n_trials_with_fiber, animal_name, ...
            RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, '');
    else
        fprintf('    No non-running epochs found.\n');
    end
    
    % Process ALL epochs (no classification)
    fprintf('  Computing all-epochs spectrogram...\n');
    all_idx = true(1, n_total);
    [mean_spec_all, cycle1_all, cycle2_all, ~] = compute_average_from_selected_epochs(...
        fiber_data.spectrograms, all_idx, phase_bin_edges, n_bins);
    
    mean_vel_all = mean(fiber_data.velocities);
    
    folder_name = sprintf('Fiber%d', fiber_idx);
    output_dir = fullfile(output_base, folder_name, 'all_epochs');
    save_results_and_figure(output_dir, signal_type, 'all_epochs', ...
        mean_spec_all, cycle1_all, cycle2_all, ...
        phase_bin_edges, fiber_data.freq_vector, ...
        n_total, mean_vel_all, n_trials_with_fiber, animal_name, ...
        RUNNING_THRESHOLD, NON_RUNNING_THRESHOLD, '');
end

%% ============================================================================
%  SUMMARY
%  ============================================================================

fprintf('\n========================================\n');
fprintf('       ANALYSIS COMPLETE               \n');
fprintf('========================================\n');
fprintf('Animal: %s\n', animal_name);
fprintf('Output folder: %s\n', output_base);
fprintf('\nVelocity Thresholds:\n');
fprintf('  Running: > %d\n', RUNNING_THRESHOLD);
fprintf('  Non-running: < %d\n', NON_RUNNING_THRESHOLD);
fprintf('\nEpoch Summary (using %.0f%% time fraction criterion):\n', MIN_TIME_FRACTION * 100);

if LFP_classification.n_total > 0
    fprintf('  LFP: %d total, %d running, %d non-running, %d intermediate\n', ...
        LFP_classification.n_total, LFP_classification.n_running, ...
        LFP_classification.n_non_running, LFP_classification.n_intermediate);
end

for fi = 1:length(Fiber_classification)
    if ~isempty(Fiber_classification{fi})
        fc = Fiber_classification{fi};
        fprintf('  Fiber%d: %d total, %d running, %d non-running, %d intermediate\n', ...
            fi, fc.n_total, fc.n_running, fc.n_non_running, fc.n_intermediate);
    end
end

fprintf('========================================\n');
fprintf('\nDone!\n');

