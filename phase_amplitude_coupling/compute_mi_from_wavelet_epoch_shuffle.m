%% ============================================================================
%  COMPUTE MI FROM WAVELET RESULTS (EPOCH SHUFFLING)
%  ============================================================================
%
%  DESCRIPTION:
%  This script computes Modulation Index (MI) using pre-computed wavelet
%  spectrogram results. It uses the 5-9 Hz phase (from LFP) and extracts
%  amplitude from the wavelet spectrogram for a specified frequency band.
%
%  SURROGATE METHOD: EPOCH SHUFFLING
%  - Randomly pairs phase from one epoch with amplitude from another epoch
%  - This preserves the statistical properties of both signals
%  - Only breaks the phase-amplitude coupling relationship
%  - This is one of the methods used in TortLab
%
%  The script:
%  1. Loads PhaseAlignedSpectrogram results
%  2. Filters running epochs based on velocity threshold
%  3. Extracts phase (5-9 Hz from LFP) and amplitude (target band from Fiber)
%  4. Pools all running epochs and computes a single MI value
%  5. Computes significance using EPOCH SHUFFLING surrogates
%  6. Generates surrogate distribution plot
%  ============================================================================

clear; clc;

% Shared MI/histogram/plotting helpers (were byte-identical local copies in
% each compute_mi_from_wavelet_*.m script; now live once in utils/).
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

%% ============================================================================
%  CONFIGURATION
%  ============================================================================

% Target frequency band for amplitude (coupling analysis)
AMP_FREQ_BAND = [30 60];         % Hz - amplitude frequency band (PAC / surrogates)

% Phase bins for MI calculation
N_PHASE_BINS = 18;

% Surrogate parameters
N_SURROGATES = 1000;             % Number of epoch-shuffled surrogates

% Running threshold (must match compute_state_spectrogram.m)
RUNNING_THRESHOLD = 2.0;         % cm/s - speed above this = running
MIN_TIME_FRACTION = 1.0;         % 100% of epoch must be above threshold

% Fiber selection
ONLY_FIBER1 = true;

% Plot settings
BIN_WIDTH_SCI = 0.2;
TARGET_BINS = 30;
MIN_BINS = 15;
N_TICKS = 5;

%% ============================================================================
%  SETUP PATHS
%  ============================================================================

script_dir = fileparts(mfilename('fullpath'));
project_dir = fileparts(script_dir);
process_dir = fullfile(project_dir, 'Process', 'phase_aligned_spectrogram');
output_dir = fullfile(project_dir, 'Process', 'mi_from_wavelet');

%% ============================================================================
%  SELECT ANIMAL FOLDER
%  ============================================================================

fprintf('Please select an animal folder in phase_aligned_spectrogram...\n');
animal_folder = uigetdir(process_dir, 'Select Animal Folder');

if isequal(animal_folder, 0)
    error('No folder selected. Exiting.');
end

% Extract animal and condition names
path_parts = strsplit(animal_folder, filesep);
animal_name = path_parts{end};
condition_name = path_parts{end-1};

fprintf('\n');
fprintf('============================================================\n');
fprintf('MI from Wavelet Results (EPOCH SHUFFLING)\n');
fprintf('============================================================\n');
fprintf('Animal: %s\n', animal_name);
fprintf('Condition: %s\n', condition_name);
fprintf('Amplitude band: %d-%d Hz\n', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
fprintf('Surrogates: %d (Epoch Shuffling)\n', N_SURROGATES);
fprintf('============================================================\n\n');

%% ============================================================================
%  FIND AND LOAD ALL TRIAL RESULTS
%  ============================================================================

% Find all PhaseAlignedSpectrogram_Results_*.mat files (recursive search)
result_files = dir(fullfile(animal_folder, '**', 'PhaseAlignedSpectrogram_Results_*.mat'));

% Exclude files in run_nonrun folder (those are aggregated results, not per-trial)
keep_idx = true(length(result_files), 1);
for i = 1:length(result_files)
    if contains(result_files(i).folder, 'run_nonrun')
        keep_idx(i) = false;
    end
end
result_files = result_files(keep_idx);

if isempty(result_files)
    error('No PhaseAlignedSpectrogram_Results_*.mat files found in: %s', animal_folder);
end

fprintf('Found %d trial result files\n', length(result_files));

% Storage for all running epochs (keep each epoch separate for shuffling)
all_epochs = struct();
all_epochs.phase = {};          % cell array, each cell is one epoch's phase
all_epochs.amplitude = {};      % cell array, each cell is one epoch's amplitude
all_epochs.lengths = [];        % track original lengths
all_epochs.trial_names = {};

% Process each trial
for fi = 1:length(result_files)
    result_path = fullfile(result_files(fi).folder, result_files(fi).name);
    fprintf('\nLoading trial %d/%d: %s\n', fi, length(result_files), result_files(fi).name);
    
    data = load(result_path);
    
    % Check required fields
    if ~isfield(data, 'results')
        fprintf('  Warning: No results field found, skipping...\n');
        continue;
    end
    
    results = data.results;
    
    % Check for LFP and Fiber data
    if ~isfield(results, 'LFP') || ~isfield(results.LFP, 'all_epochs_raw')
        fprintf('  Warning: No LFP all_epochs_raw found, skipping...\n');
        continue;
    end
    
    fiber_field = 'Fiber1';
    if ~isfield(results, fiber_field) || ~isfield(results.(fiber_field), 'all_epochs_raw')
        fprintf('  Warning: No %s all_epochs_raw found, skipping...\n', fiber_field);
        continue;
    end
    
    % Get frequency vector
    freq_vector = results.metadata.freq_vector;
    
    % Find amplitude band indices
    amp_band_idx = (freq_vector >= AMP_FREQ_BAND(1)) & (freq_vector <= AMP_FREQ_BAND(2));
    if sum(amp_band_idx) == 0
        fprintf('  Warning: Amplitude band %d-%d Hz not covered by freq_vector, skipping...\n', ...
            AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
        continue;
    end
    
    n_epochs = length(results.LFP.all_epochs_raw);
    fprintf('  Found %d epochs\n', n_epochs);
    
    % Process each epoch
    for ei = 1:n_epochs
        lfp_epoch = results.LFP.all_epochs_raw{ei};
        fiber_epoch = results.(fiber_field).all_epochs_raw{ei};
        
        % Check if epoch is running
        if isfield(lfp_epoch, 'velocity') && ~isempty(lfp_epoch.velocity)
            vel = lfp_epoch.velocity(:);
            running_fraction = sum(vel > RUNNING_THRESHOLD) / length(vel);
            is_running = running_fraction >= MIN_TIME_FRACTION;
        elseif isfield(lfp_epoch, 'mean_velocity')
            is_running = lfp_epoch.mean_velocity > RUNNING_THRESHOLD;
        else
            is_running = false;
        end
        
        if ~is_running
            continue;
        end
        
        % Extract phase from LFP epoch
        lfp_phase = lfp_epoch.phase(:)';
        
        % Extract amplitude from Fiber epoch spectrogram
        fiber_spec = fiber_epoch.spectrogram;  % [n_freqs x n_timepoints]
        fiber_amp = mean(fiber_spec(amp_band_idx, :), 1);  % average over frequency band
        
        % Ensure same length within this epoch
        min_len = min(length(lfp_phase), length(fiber_amp));
        lfp_phase = lfp_phase(1:min_len);
        fiber_amp = fiber_amp(1:min_len);
        
        % Store each epoch separately (for shuffling)
        all_epochs.phase{end+1} = lfp_phase;
        all_epochs.amplitude{end+1} = fiber_amp;
        all_epochs.lengths(end+1) = min_len;
        all_epochs.trial_names{end+1} = result_files(fi).name;
    end
end

n_total_epochs = length(all_epochs.phase);

fprintf('\n============================================================\n');
fprintf('Total running epochs collected: %d\n', n_total_epochs);
fprintf('============================================================\n');

if n_total_epochs < 2
    error('Need at least 2 epochs for epoch shuffling. Found: %d', n_total_epochs);
end

%% ============================================================================
%  COMPUTE OBSERVED MI (FULL DATA - for reporting)
%  ============================================================================
% This is the "true" observed MI using ALL available data

fprintf('\nComputing observed MI (full data)...\n');

% Concatenate all epochs with original lengths
phase_concat_full = [all_epochs.phase{:}];
amp_concat_full = [all_epochs.amplitude{:}];
n_datapoints_full = length(phase_concat_full);

observed_mi_full = compute_tort_mi(phase_concat_full, amp_concat_full, N_PHASE_BINS);
fprintf('  Observed MI (full) = %.8f  [%d datapoints]\n', observed_mi_full, n_datapoints_full);

%% ============================================================================
%  UNIFY EPOCH LENGTHS (for fair surrogate comparison)
%  ============================================================================
% To ensure surrogate MI uses the same amount of data as the matched observed MI,
% truncate all epochs to the minimum length across all epochs.

unified_length = min(all_epochs.lengths);
data_loss_pct = (1 - unified_length * n_total_epochs / n_datapoints_full) * 100;

fprintf('\nUnifying epoch lengths for surrogate comparison:\n');
fprintf('  min = %d, max = %d -> unified = %d\n', ...
    min(all_epochs.lengths), max(all_epochs.lengths), unified_length);
fprintf('  Data loss: %.2f%%\n', data_loss_pct);

% Create matched-length copies for surrogate comparison
all_epochs_matched = struct();
all_epochs_matched.phase = cell(1, n_total_epochs);
all_epochs_matched.amplitude = cell(1, n_total_epochs);

for ei = 1:n_total_epochs
    all_epochs_matched.phase{ei} = all_epochs.phase{ei}(1:unified_length);
    all_epochs_matched.amplitude{ei} = all_epochs.amplitude{ei}(1:unified_length);
end

%% ============================================================================
%  COMPUTE OBSERVED MI (MATCHED DATA - for significance testing)
%  ============================================================================

fprintf('\nComputing observed MI (matched length)...\n');

phase_concat_matched = [all_epochs_matched.phase{:}];
amp_concat_matched = [all_epochs_matched.amplitude{:}];
n_datapoints_matched = length(phase_concat_matched);

observed_mi_matched = compute_tort_mi(phase_concat_matched, amp_concat_matched, N_PHASE_BINS);
fprintf('  Observed MI (matched) = %.8f  [%d datapoints]\n', observed_mi_matched, n_datapoints_matched);

%% ============================================================================
%  COMPUTE SURROGATE DISTRIBUTION (EPOCH SHUFFLING)
%  ============================================================================

fprintf('\nComputing %d epoch-shuffled surrogates (using matched-length data)...\n', N_SURROGATES);

surrogate_mi = zeros(1, N_SURROGATES);
n_epochs = n_total_epochs;

for si = 1:N_SURROGATES
    if mod(si, 200) == 0
        fprintf('  Surrogate %d/%d\n', si, N_SURROGATES);
    end
    
    % Shuffle epoch indices for amplitude (keep phase indices fixed)
    shuffled_amp_idx = randperm(n_epochs);
    
    % Ensure no epoch is paired with itself
    while any(shuffled_amp_idx == 1:n_epochs)
        shuffled_amp_idx = randperm(n_epochs);
    end
    
    % Concatenate shuffled pairs
    % (All matched epochs have the same unified length)
    phase_shuffled = [];
    amp_shuffled = [];
    
    for ei = 1:n_epochs
        phase_shuffled = [phase_shuffled, all_epochs_matched.phase{ei}];
        amp_shuffled = [amp_shuffled, all_epochs_matched.amplitude{shuffled_amp_idx(ei)}];
    end
    
    % Compute MI with shuffled pairing
    surrogate_mi(si) = compute_tort_mi(phase_shuffled, amp_shuffled, N_PHASE_BINS);
end

fprintf('Done.\n');

%% ============================================================================
%  COMPUTE SIGNIFICANCE (using matched MI for fair comparison)
%  ============================================================================

% One-sided p-value: proportion of surrogates >= observed (matched)
p_value = (sum(surrogate_mi >= observed_mi_matched) + 1) / (N_SURROGATES + 1);

% Z-score (using matched MI)
z_score = (observed_mi_matched - mean(surrogate_mi)) / std(surrogate_mi);

fprintf('\n============================================================\n');
fprintf('SIGNIFICANCE RESULTS (Epoch Shuffling)\n');
fprintf('============================================================\n');
fprintf('Observed MI (FULL, for reporting):    %.8f  [%d pts]\n', observed_mi_full, n_datapoints_full);
fprintf('Observed MI (MATCHED, for testing):   %.8f  [%d pts]\n', observed_mi_matched, n_datapoints_matched);
fprintf('Surrogate mean:                       %.8f\n', mean(surrogate_mi));
fprintf('Surrogate std:                        %.8f\n', std(surrogate_mi));
fprintf('P-value (based on matched MI):        %.6f\n', p_value);
fprintf('Z-score (based on matched MI):        %.3f\n', z_score);
fprintf('============================================================\n');

%% ============================================================================
%  SAVE RESULTS
%  ============================================================================

output_path = fullfile(output_dir, condition_name, animal_name);
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

Results = struct();
Results.animal_name = animal_name;
Results.condition = condition_name;
Results.timestamp = datetime('now');
Results.surrogate_method = 'epoch_shuffling';

Results.parameters.amp_freq_band = AMP_FREQ_BAND;
Results.parameters.n_phase_bins = N_PHASE_BINS;
Results.parameters.n_surrogates = N_SURROGATES;
Results.parameters.running_threshold = RUNNING_THRESHOLD;
Results.parameters.min_time_fraction = MIN_TIME_FRACTION;

% Store both MI values
Results.observed_mi_full = observed_mi_full;           % True MI using all data (for reporting)
Results.observed_mi_matched = observed_mi_matched;     % MI using matched-length data (for testing)
Results.observed_mi = observed_mi_full;                % Default to full MI for backward compatibility

Results.surrogate_mi = surrogate_mi;
Results.p_value = p_value;
Results.z_score = z_score;
Results.n_epochs = n_total_epochs;
Results.n_datapoints_full = n_datapoints_full;
Results.n_datapoints_matched = n_datapoints_matched;
Results.unified_epoch_length = unified_length;
Results.data_loss_pct = data_loss_pct;

save(fullfile(output_path, 'MI_EpochShuffle_Results.mat'), 'Results', '-v7');
fprintf('\nSaved: %s\n', fullfile(output_path, 'MI_EpochShuffle_Results.mat'));

%% ============================================================================
%  PLOT SURROGATE DISTRIBUTION
%  ============================================================================

fprintf('\nGenerating surrogate distribution plot...\n');

fig = figure('Position', [200, 200, 600, 450], 'Visible', 'off');

% Use FULL observed MI for the plot (this is the true value we report)
% But p-value is calculated using MATCHED MI (for fair comparison)
observed_mi_for_plot = observed_mi_full;

% Compute unified bin width
all_values = [surrogate_mi(:); observed_mi_for_plot];
unified_bin_width = compute_unified_bin_width(all_values, BIN_WIDTH_SCI, TARGET_BINS, MIN_BINS);

% Build bin edges (starting from 0)
edges = build_bin_edges_with_obs(surrogate_mi, observed_mi_for_plot, unified_bin_width);

% Plot histogram
histogram(surrogate_mi, edges, 'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none');
hold on;

% Add observed MI line (FULL MI - the true reported value)
xline(observed_mi_for_plot, 'r-', 'LineWidth', 2);

% Apply simple ticks
max_val = max(all_values);
if max_val > 0
    unified_exp = floor(log10(max_val));
    max_scaled = max_val / 10^unified_exp;
    if max_scaled < 2
        unified_exp = unified_exp - 1;
    end
else
    unified_exp = -4;
end
apply_simple_ticks(gca, edges, N_TICKS, unified_exp);

xlabel(sprintf('Surrogate MI (\\times10^{%d})', unified_exp));
ylabel('Count');
title_str = sprintf('LFP (5-9 Hz) → Fiber (%d-%d Hz) - Running [Epoch Shuffle]', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
title(title_str, 'Interpreter', 'none');

% Add significance stars (top-left)
text(0.02, 0.95, p_to_stars(p_value), 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14);

% Add MI value annotation (top-right)
text(0.98, 0.95, sprintf('MI = %.6f', observed_mi_full), 'Units', 'normalized', ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', 10);

hold off;

% Save figure
saveas(fig, fullfile(output_path, sprintf('SurrogateDist_%d-%dHz_running_EpochShuffle.png', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2))));
saveas(fig, fullfile(output_path, sprintf('SurrogateDist_%d-%dHz_running_EpochShuffle.fig', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2))));
close(fig);

fprintf('Saved surrogate distribution plot.\n');
fprintf('\nDone!\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================
%  compute_tort_mi, compute_unified_bin_width, build_bin_edges_with_obs,
%  apply_simple_ticks, p_to_stars moved to utils/ (shared with the other
%  compute_mi_from_wavelet_*.m scripts; see addpath near the top of this file).
