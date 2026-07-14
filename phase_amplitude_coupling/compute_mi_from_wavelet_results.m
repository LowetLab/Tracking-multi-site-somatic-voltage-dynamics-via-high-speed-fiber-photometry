%% ============================================================================
%  COMPUTE MI FROM WAVELET RESULTS
%  ============================================================================
%
%  DESCRIPTION:
%  This script computes Modulation Index (MI) using pre-computed wavelet
%  spectrogram results. It uses the 5-9 Hz phase (from LFP) and extracts
%  amplitude from the wavelet spectrogram for a specified frequency band.
%
%  The script:
%  1. Loads PhaseAlignedSpectrogram results
%  2. Filters running epochs based on velocity threshold
%  3. Extracts phase (5-9 Hz from LFP) and amplitude (target band from Fiber)
%  4. Pools all running epochs and computes a single MI value
%  5. Computes significance using IAAFT surrogates
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
N_SURROGATES = 1000;
IAAFT_ITER = 50;

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
fprintf('MI from Wavelet Results\n');
fprintf('============================================================\n');
fprintf('Animal: %s\n', animal_name);
fprintf('Condition: %s\n', condition_name);
fprintf('Amplitude band: %d-%d Hz\n', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
fprintf('Surrogates: %d (IAAFT iters: %d)\n', N_SURROGATES, IAAFT_ITER);
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

% Storage for all running epochs
all_running_epochs = struct();
all_running_epochs.lfp_phase = [];      % concatenated phase
all_running_epochs.fiber_amp = [];      % concatenated amplitude
all_running_epochs.n_epochs = 0;
all_running_epochs.trial_names = {};

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
        
        % Ensure same length
        min_len = min(length(lfp_phase), length(fiber_amp));
        lfp_phase = lfp_phase(1:min_len);
        fiber_amp = fiber_amp(1:min_len);
        
        % Append to storage
        all_running_epochs.lfp_phase = [all_running_epochs.lfp_phase, lfp_phase];
        all_running_epochs.fiber_amp = [all_running_epochs.fiber_amp, fiber_amp];
        all_running_epochs.n_epochs = all_running_epochs.n_epochs + 1;
        all_running_epochs.trial_names{end+1} = result_files(fi).name;
    end
end

fprintf('\n============================================================\n');
fprintf('Total running epochs collected: %d\n', all_running_epochs.n_epochs);
fprintf('Total data points: %d\n', length(all_running_epochs.lfp_phase));
fprintf('============================================================\n');

if all_running_epochs.n_epochs == 0
    error('No running epochs found. Check velocity threshold or data.');
end

%% ============================================================================
%  COMPUTE OBSERVED MI
%  ============================================================================

fprintf('\nComputing observed MI...\n');

phase_data = all_running_epochs.lfp_phase;
amp_data = all_running_epochs.fiber_amp;

observed_mi = compute_tort_mi(phase_data, amp_data, N_PHASE_BINS);
fprintf('  Observed MI = %.6e\n', observed_mi);

%% ============================================================================
%  COMPUTE SURROGATE DISTRIBUTION (IAAFT)
%  ============================================================================

fprintf('\nComputing %d IAAFT surrogates...\n', N_SURROGATES);

surrogate_mi = zeros(1, N_SURROGATES);

for si = 1:N_SURROGATES
    if mod(si, 100) == 0
        fprintf('  Surrogate %d/%d\n', si, N_SURROGATES);
    end
    
    % Generate IAAFT surrogate for amplitude data
    amp_surrogate = iaaft_surrogate(amp_data, IAAFT_ITER);
    
    % Compute MI with surrogate amplitude
    surrogate_mi(si) = compute_tort_mi(phase_data, amp_surrogate, N_PHASE_BINS);
end

fprintf('Done.\n');

%% ============================================================================
%  COMPUTE SIGNIFICANCE
%  ============================================================================

% One-sided p-value: proportion of surrogates >= observed
p_value = (sum(surrogate_mi >= observed_mi) + 1) / (N_SURROGATES + 1);

% Z-score
z_score = (observed_mi - mean(surrogate_mi)) / std(surrogate_mi);

fprintf('\n============================================================\n');
fprintf('SIGNIFICANCE RESULTS\n');
fprintf('============================================================\n');
fprintf('Observed MI: %.6e\n', observed_mi);
fprintf('Surrogate mean: %.6e\n', mean(surrogate_mi));
fprintf('Surrogate std: %.6e\n', std(surrogate_mi));
fprintf('P-value: %.6f\n', p_value);
fprintf('Z-score: %.3f\n', z_score);
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

Results.parameters.amp_freq_band = AMP_FREQ_BAND;
Results.parameters.n_phase_bins = N_PHASE_BINS;
Results.parameters.n_surrogates = N_SURROGATES;
Results.parameters.iaaft_iter = IAAFT_ITER;
Results.parameters.running_threshold = RUNNING_THRESHOLD;
Results.parameters.min_time_fraction = MIN_TIME_FRACTION;

Results.observed_mi = observed_mi;
Results.surrogate_mi = surrogate_mi;
Results.p_value = p_value;
Results.z_score = z_score;
Results.n_epochs = all_running_epochs.n_epochs;
Results.n_datapoints = length(phase_data);

save(fullfile(output_path, 'MI_FromWavelet_Results.mat'), 'Results', '-v7');
fprintf('\nSaved: %s\n', fullfile(output_path, 'MI_FromWavelet_Results.mat'));

%% ============================================================================
%  PLOT SURROGATE DISTRIBUTION
%  ============================================================================

fprintf('\nGenerating surrogate distribution plot...\n');

fig = figure('Position', [200, 200, 600, 450], 'Visible', 'off');

% Compute unified bin width
all_values = [surrogate_mi(:); observed_mi];
unified_bin_width = compute_unified_bin_width(all_values, BIN_WIDTH_SCI, TARGET_BINS, MIN_BINS);

% Build bin edges (starting from 0)
edges = build_bin_edges_with_obs(surrogate_mi, observed_mi, unified_bin_width);

% Plot histogram
histogram(surrogate_mi, edges, 'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none');
hold on;

% Add observed MI line
xline(observed_mi, 'r-', 'LineWidth', 2);

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
title_str = sprintf('LFP (5-9 Hz) → Fiber (%d-%d Hz) - Running', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
title(title_str, 'Interpreter', 'none');

% Add significance stars
text(0.02, 0.95, p_to_stars(p_value), 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'FontSize', 14);

hold off;

% Save figure
saveas(fig, fullfile(output_path, sprintf('SurrogateDist_%d-%dHz_running.png', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2))));
saveas(fig, fullfile(output_path, sprintf('SurrogateDist_%d-%dHz_running.fig', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2))));
close(fig);

fprintf('Saved surrogate distribution plot.\n');
fprintf('\nDone!\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================
%  compute_tort_mi, compute_unified_bin_width, build_bin_edges_with_obs,
%  apply_simple_ticks, p_to_stars moved to utils/ (shared with the other
%  compute_mi_from_wavelet_*.m scripts; see addpath near the top of this file).

function surrogate = iaaft_surrogate(signal, n_iter)
    % Generate IAAFT surrogate
    % Preserves amplitude distribution and power spectrum
    
    signal = signal(:);
    n = length(signal);
    
    % Store original amplitude spectrum
    fft_orig = fft(signal);
    amp_spectrum = abs(fft_orig);
    
    % Store sorted original values
    sorted_signal = sort(signal);
    
    % Initialize with random shuffle
    surrogate = signal(randperm(n));
    
    for iter = 1:n_iter
        % Impose amplitude spectrum
        fft_surr = fft(surrogate);
        phase_surr = angle(fft_surr);
        fft_surr = amp_spectrum .* exp(1i * phase_surr);
        surrogate = real(ifft(fft_surr));
        
        % Impose amplitude distribution (rank ordering)
        [~, idx] = sort(surrogate);
        surrogate(idx) = sorted_signal;
    end
end
