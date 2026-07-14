%% ============================================================================
%  COMPUTE MI FROM WAVELET RESULTS (IAAFT PER EPOCH)
%  ============================================================================
%
%  DESCRIPTION:
%  This script computes Modulation Index (MI) using pre-computed wavelet
%  spectrogram results. It uses the 5-9 Hz phase (from LFP) and extracts
%  amplitude from the wavelet spectrogram for a specified frequency band.
%
%  SURROGATE METHOD: IAAFT PER EPOCH
%  - For each surrogate iteration:
%    - Apply IAAFT independently to each epoch's amplitude signal
%    - This preserves power spectrum and amplitude distribution within each epoch
%    - Randomizes phase relationships within each epoch
%  - Then concatenate all IAAFT-processed epochs to compute surrogate MI
%  - This method respects epoch boundaries and applies IAAFT to continuous signals
%
%  The script:
%  1. Loads PhaseAlignedSpectrogram results
%  2. Filters running epochs based on velocity threshold
%  3. Extracts phase (5-9 Hz from LFP) and amplitude (target band from Fiber)
%  4. Pools all running epochs and computes a single MI value
%  5. Computes significance using IAAFT surrogates (per epoch)
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
N_SURROGATES = 1000;             % Number of IAAFT surrogates
IAAFT_ITER = 50;                 % IAAFT iterations per surrogate

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
fprintf('MI from Wavelet Results (IAAFT PER EPOCH)\n');
fprintf('============================================================\n');
fprintf('Animal: %s\n', animal_name);
fprintf('Condition: %s\n', condition_name);
fprintf('Amplitude band: %d-%d Hz\n', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2));
fprintf('Surrogates: %d (IAAFT per epoch, %d iterations)\n', N_SURROGATES, IAAFT_ITER);
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

% Storage for all running epochs (keep each epoch separate)
all_epochs = struct();
all_epochs.phase = {};              % cell array, each cell is one epoch's LFP phase
all_epochs.fiber_amplitude = {};    % cell array, each cell is one epoch's Fiber amplitude
all_epochs.lfp_amplitude = {};      % cell array, each cell is one epoch's LFP amplitude
all_epochs.lengths = [];            % track original lengths
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
        
        % Extract amplitude from LFP epoch spectrogram
        lfp_spec = lfp_epoch.spectrogram;  % [n_freqs x n_timepoints]
        lfp_amp = mean(lfp_spec(amp_band_idx, :), 1);  % average over frequency band
        
        % Ensure same length within this epoch
        min_len = min([length(lfp_phase), length(fiber_amp), length(lfp_amp)]);
        lfp_phase = lfp_phase(1:min_len);
        fiber_amp = fiber_amp(1:min_len);
        lfp_amp = lfp_amp(1:min_len);
        
        % Store each epoch separately
        all_epochs.phase{end+1} = lfp_phase;
        all_epochs.fiber_amplitude{end+1} = fiber_amp;
        all_epochs.lfp_amplitude{end+1} = lfp_amp;
        all_epochs.lengths(end+1) = min_len;
        all_epochs.trial_names{end+1} = result_files(fi).name;
    end
end

n_total_epochs = length(all_epochs.phase);

fprintf('\n============================================================\n');
fprintf('Total running epochs collected: %d\n', n_total_epochs);
fprintf('============================================================\n');

if n_total_epochs < 1
    error('Need at least 1 epoch. Found: %d', n_total_epochs);
end

%% ============================================================================
%  COMPUTE OBSERVED MI (both LFP→Fiber and LFP→LFP)
%  ============================================================================
% Note: Unlike Epoch Shuffling, IAAFT per epoch does NOT require length unification
% because each epoch's phase and amplitude always stay paired within the same epoch.

fprintf('\nComputing observed MI...\n');

% Concatenate all epochs with original lengths
phase_concat = [all_epochs.phase{:}];
fiber_amp_concat = [all_epochs.fiber_amplitude{:}];
lfp_amp_concat = [all_epochs.lfp_amplitude{:}];
n_datapoints = length(phase_concat);

% LFP phase → Fiber amplitude
observed_mi_fiber = compute_tort_mi(phase_concat, fiber_amp_concat, N_PHASE_BINS);
fprintf('  Observed MI (LFP→Fiber) = %.8f  [%d datapoints]\n', observed_mi_fiber, n_datapoints);

% LFP phase → LFP amplitude
observed_mi_lfp = compute_tort_mi(phase_concat, lfp_amp_concat, N_PHASE_BINS);
fprintf('  Observed MI (LFP→LFP)   = %.8f  [%d datapoints]\n', observed_mi_lfp, n_datapoints);

fprintf('\nEpoch length statistics:\n');
fprintf('  min = %d, max = %d, mean = %.1f\n', ...
    min(all_epochs.lengths), max(all_epochs.lengths), mean(all_epochs.lengths));

%% ============================================================================
%  COMPUTE SURROGATE DISTRIBUTION (IAAFT PER EPOCH) - PARALLEL
%  ============================================================================

fprintf('\nComputing %d IAAFT surrogates (per epoch, %d iterations each)...\n', N_SURROGATES, IAAFT_ITER);

% Set and record random seed for reproducibility
rng_seed = randi(2^31 - 1);  % Generate a random seed
rng(rng_seed);               % Set the seed
fprintf('Random seed: %d\n', rng_seed);

% Check and start parallel pool
pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting parallel pool...\n');
    pool = parpool('local');
end
fprintf('Using %d workers for parallel computation.\n', pool.NumWorkers);

surrogate_mi_fiber = zeros(1, N_SURROGATES);
surrogate_mi_lfp = zeros(1, N_SURROGATES);
n_epochs = n_total_epochs;

% Pre-concatenate phase (it never changes across surrogates)
phase_all = phase_concat;  % Already computed earlier

% Extract amplitude data to simple cell arrays for parfor
fiber_amp_data = all_epochs.fiber_amplitude;
lfp_amp_data = all_epochs.lfp_amplitude;

% Local copies of parameters for parfor
iaaft_iter_local = IAAFT_ITER;
n_phase_bins_local = N_PHASE_BINS;

fprintf('Computing surrogates for both LFP→Fiber and LFP→LFP...\n');
tic;

parfor si = 1:N_SURROGATES
    % For each epoch, apply IAAFT to both amplitudes independently
    fiber_surr_parts = cell(1, n_epochs);
    lfp_surr_parts = cell(1, n_epochs);
    
    for ei = 1:n_epochs
        % IAAFT for Fiber amplitude
        fiber_surr_parts{ei} = iaaft_surrogate_parfor(fiber_amp_data{ei}, iaaft_iter_local);
        % IAAFT for LFP amplitude
        lfp_surr_parts{ei} = iaaft_surrogate_parfor(lfp_amp_data{ei}, iaaft_iter_local);
    end
    
    % Concatenate all parts
    fiber_amp_surr = [fiber_surr_parts{:}];
    lfp_amp_surr = [lfp_surr_parts{:}];
    
    % Compute MI with IAAFT-processed amplitudes
    surrogate_mi_fiber(si) = compute_tort_mi_parfor(phase_all, fiber_amp_surr, n_phase_bins_local);
    surrogate_mi_lfp(si) = compute_tort_mi_parfor(phase_all, lfp_amp_surr, n_phase_bins_local);
end

total_time = toc;
fprintf('Done. Total time: %.1f min (%.2f sec per surrogate)\n', total_time/60, total_time/N_SURROGATES);

%% ============================================================================
%  COMPUTE SIGNIFICANCE (for both LFP→Fiber and LFP→LFP)
%  ============================================================================

% LFP→Fiber significance
p_value_fiber = (sum(surrogate_mi_fiber >= observed_mi_fiber) + 1) / (N_SURROGATES + 1);
z_score_fiber = (observed_mi_fiber - mean(surrogate_mi_fiber)) / std(surrogate_mi_fiber);

% LFP→LFP significance
p_value_lfp = (sum(surrogate_mi_lfp >= observed_mi_lfp) + 1) / (N_SURROGATES + 1);
z_score_lfp = (observed_mi_lfp - mean(surrogate_mi_lfp)) / std(surrogate_mi_lfp);

fprintf('\n============================================================\n');
fprintf('SIGNIFICANCE RESULTS (IAAFT Per Epoch)\n');
fprintf('============================================================\n');
fprintf('LFP → Fiber:\n');
fprintf('  Observed MI:      %.8f  [%d pts]\n', observed_mi_fiber, n_datapoints);
fprintf('  Surrogate mean:   %.8f\n', mean(surrogate_mi_fiber));
fprintf('  Surrogate std:    %.8f\n', std(surrogate_mi_fiber));
fprintf('  P-value:          %.6f\n', p_value_fiber);
fprintf('  Z-score:          %.3f\n', z_score_fiber);
fprintf('------------------------------------------------------------\n');
fprintf('LFP → LFP:\n');
fprintf('  Observed MI:      %.8f  [%d pts]\n', observed_mi_lfp, n_datapoints);
fprintf('  Surrogate mean:   %.8f\n', mean(surrogate_mi_lfp));
fprintf('  Surrogate std:    %.8f\n', std(surrogate_mi_lfp));
fprintf('  P-value:          %.6f\n', p_value_lfp);
fprintf('  Z-score:          %.3f\n', z_score_lfp);
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
Results.surrogate_method = 'iaaft_per_epoch';

Results.parameters.amp_freq_band = AMP_FREQ_BAND;
Results.parameters.n_phase_bins = N_PHASE_BINS;
Results.parameters.n_surrogates = N_SURROGATES;
Results.parameters.iaaft_iterations = IAAFT_ITER;
Results.parameters.running_threshold = RUNNING_THRESHOLD;
Results.parameters.min_time_fraction = MIN_TIME_FRACTION;

% LFP→Fiber results
Results.LFP_to_Fiber.observed_mi = observed_mi_fiber;
Results.LFP_to_Fiber.surrogate_mi = surrogate_mi_fiber;
Results.LFP_to_Fiber.p_value = p_value_fiber;
Results.LFP_to_Fiber.z_score = z_score_fiber;

% LFP→LFP results
Results.LFP_to_LFP.observed_mi = observed_mi_lfp;
Results.LFP_to_LFP.surrogate_mi = surrogate_mi_lfp;
Results.LFP_to_LFP.p_value = p_value_lfp;
Results.LFP_to_LFP.z_score = z_score_lfp;

Results.n_epochs = n_total_epochs;
Results.n_datapoints = n_datapoints;
Results.epoch_lengths = all_epochs.lengths;
Results.computation_time_sec = total_time;
Results.rng_seed = rng_seed;  % Random seed for reproducibility

save(fullfile(output_path, 'MI_IAAFT_PerEpoch_Results.mat'), 'Results', '-v7');
fprintf('\nSaved: %s\n', fullfile(output_path, 'MI_IAAFT_PerEpoch_Results.mat'));

%% ============================================================================
%  PLOT SURROGATE DISTRIBUTIONS (LOG10 SCALE)
%  ============================================================================

fprintf('\nGenerating surrogate distribution plots (log10 scale)...\n');

% Helper function for plotting
plot_surrogate_dist = @(surr_mi, obs_mi, p_val, title_str, filename) plot_log_surrogate(...
    surr_mi, obs_mi, p_val, title_str, filename, output_path);

% Plot 1: LFP → Fiber
plot_surrogate_dist(surrogate_mi_fiber, observed_mi_fiber, p_value_fiber, ...
    sprintf('LFP (5-9 Hz) → Fiber (%d-%d Hz) - Running [IAAFT/Epoch]', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2)), ...
    sprintf('SurrogateDist_%d-%dHz_LFP-Fiber_IAAFT_PerEpoch', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2)));

% Plot 2: LFP → LFP
plot_surrogate_dist(surrogate_mi_lfp, observed_mi_lfp, p_value_lfp, ...
    sprintf('LFP (5-9 Hz) → LFP (%d-%d Hz) - Running [IAAFT/Epoch]', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2)), ...
    sprintf('SurrogateDist_%d-%dHz_LFP-LFP_IAAFT_PerEpoch', AMP_FREQ_BAND(1), AMP_FREQ_BAND(2)));

fprintf('Saved surrogate distribution plots.\n');
fprintf('\nDone!\n');

%% ============================================================================
%  PLOTTING FUNCTION
%  ============================================================================

function plot_log_surrogate(surrogate_mi, observed_mi, p_value, title_str, filename, output_path)
    fig = figure('Position', [200, 200, 600, 450], 'Visible', 'off');
    
    % Convert to log10 scale (filter out zeros/negatives if any)
    surrogate_mi_valid = surrogate_mi(surrogate_mi > 0);
    log_surrogate_mi = log10(surrogate_mi_valid);
    log_observed_mi = log10(observed_mi);
    
    % Compute bin edges for log scale
    all_log_values = [log_surrogate_mi(:); log_observed_mi];
    log_min = floor(min(all_log_values) * 10) / 10;  % Round down to 0.1
    log_max = ceil(max(all_log_values) * 10) / 10;   % Round up to 0.1
    bin_width = 0.1;  % 0.1 in log10 scale
    edges = log_min:bin_width:log_max;
    
    % Ensure we have enough bins
    if length(edges) < 15
        bin_width = (log_max - log_min) / 20;
        edges = log_min:bin_width:log_max;
    end
    
    % Plot histogram
    histogram(log_surrogate_mi, edges, 'FaceColor', [0.4 0.4 0.8], 'EdgeColor', 'none');
    hold on;
    
    % Add observed MI line
    xline(log_observed_mi, 'r-', 'LineWidth', 2);
    
    % Set x-axis limits
    xlim([log_min, log_max]);
    
    % Create nice tick labels (show as 10^x, integers only)
    tick_values = ceil(log_min):1:floor(log_max);
    % Ensure at least 2 ticks
    if length(tick_values) < 2
        tick_values = [floor(log_min), ceil(log_max)];
    end
    set(gca, 'XTick', tick_values);
    tick_labels_arr = cell(size(tick_values));
    for ti = 1:length(tick_values)
        tick_labels_arr{ti} = sprintf('10^{%d}', tick_values(ti));  % Integer only
    end
    set(gca, 'XTickLabel', tick_labels_arr);
    
    xlabel('log_{10}(MI)');
    ylabel('Count');
    title(title_str, 'Interpreter', 'none');
    
    % Add significance stars and observed MI legend (top-right, inside plot area)
    % Format MI as "X.XXX×10^Y"
    mi_str = format_scientific(observed_mi);
    
    % Create two-line legend box in top-right (inside plot area)
    % Line 1: significance stars, Line 2: Observed MI = value
    annotation_str = sprintf('%s\nObserved MI = %s', p_to_stars_local(p_value), mi_str);
    text(0.78, 0.92, annotation_str, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 10, 'BackgroundColor', 'white', 'EdgeColor', 'black', ...
        'LineWidth', 0.5, 'Margin', 4);
    
    function str = format_scientific(val)
        if val == 0
            str = '0';
        else
            exp_val = floor(log10(abs(val)));
            mantissa = val / 10^exp_val;
            str = sprintf('%.3f\\times10^{%d}', mantissa, exp_val);
        end
    end
    
    hold off;
    
    % Save figure
    saveas(fig, fullfile(output_path, [filename '.png']));
    saveas(fig, fullfile(output_path, [filename '.fig']));
    close(fig);
    
    function stars = p_to_stars_local(pv)
        if ~isfinite(pv)
            stars = 'n.s.';
        elseif pv < 0.001
            stars = '***';
        elseif pv < 0.01
            stars = '**';
        elseif pv < 0.05
            stars = '*';
        else
            stars = 'n.s.';
        end
    end
end

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================
%  compute_tort_mi, compute_unified_bin_width, build_bin_edges_with_obs,
%  apply_simple_ticks, p_to_stars moved to utils/ (shared with the other
%  compute_mi_from_wavelet_*.m scripts; see addpath near the top of this file).
%  compute_tort_mi_parfor below is a separate, parfor-safe duplicate -- see
%  the PARFOR-COMPATIBLE FUNCTIONS section.

function surrogate = iaaft_surrogate(data, n_iter)
    % IAAFT (Iterative Amplitude Adjusted Fourier Transform)
    % Generates a surrogate that preserves power spectrum and amplitude distribution
    %
    % data: input time series (row vector)
    % n_iter: number of iterations
    
    data = data(:)';
    n = length(data);
    
    % Store original sorted values and power spectrum
    sorted_data = sort(data);
    fft_orig = fft(data);
    amp_orig = abs(fft_orig);
    
    % Step 1: Start with random shuffle of original data
    surrogate = data(randperm(n));
    
    % Iterative refinement
    for iter = 1:n_iter
        % Step 2: FFT constraint - impose original power spectrum
        fft_surr = fft(surrogate);
        phase_surr = angle(fft_surr);
        fft_new = amp_orig .* exp(1i * phase_surr);
        surrogate = real(ifft(fft_new));
        
        % Step 3: Amplitude constraint - impose original amplitude distribution
        [~, sort_idx] = sort(surrogate);
        surrogate(sort_idx) = sorted_data;
    end
end

%% ============================================================================
%  PARFOR-COMPATIBLE FUNCTIONS
%  ============================================================================
% These are duplicates of the above functions, designed to work inside parfor
% (parfor requires functions to be defined at end of file, not nested)

function mi = compute_tort_mi_parfor(phase, amp, nbin)
    phase = phase(:)';
    amp = amp(:)';
    
    bin_edges = linspace(-pi, pi, nbin + 1);
    
    mean_amp = zeros(1, nbin);
    for bi = 1:nbin
        if bi < nbin
            idx = (phase >= bin_edges(bi)) & (phase < bin_edges(bi+1));
        else
            idx = (phase >= bin_edges(bi)) & (phase <= bin_edges(bi+1));
        end
        if sum(idx) > 0
            mean_amp(bi) = mean(amp(idx));
        end
    end
    
    if sum(mean_amp) > 0
        p = mean_amp / sum(mean_amp);
    else
        mi = 0;
        return;
    end
    
    q = ones(1, nbin) / nbin;
    
    nonzero_idx = (p > 0);
    kl_div = sum(p(nonzero_idx) .* log(p(nonzero_idx) ./ q(nonzero_idx)));
    
    mi = kl_div / log(nbin);
end

function surrogate = iaaft_surrogate_parfor(data, n_iter)
    data = data(:)';
    n = length(data);
    
    sorted_data = sort(data);
    fft_orig = fft(data);
    amp_orig = abs(fft_orig);
    
    surrogate = data(randperm(n));
    
    for iter = 1:n_iter
        fft_surr = fft(surrogate);
        phase_surr = angle(fft_surr);
        fft_new = amp_orig .* exp(1i * phase_surr);
        surrogate = real(ifft(fft_new));
        
        [~, sort_idx] = sort(surrogate);
        surrogate(sort_idx) = sorted_data;
    end
end
