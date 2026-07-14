%% EXTRACT_SPIKE_ML_DATA
%  Extracts a concise spike-centric data structure from a CellularAnalysis
%  output file, intended for sharing with ML collaborators for spike
%  detection model training.
%
%  INPUT:  CellularAnalysis .mat file from cellular_processing_multitrial_DBS.m
%  OUTPUT: SpikeMLData .mat file with spike snippets, full traces, and
%          per-neuron + pooled spike matrices.
%
%  USAGE:
%    1. Set INPUT_FILE and OUTPUT_FILE below
%    2. Optionally adjust SNIPPET_WINDOW_MS
%    3. Run the script
%
clear; clc;

%% ========================================================================
%  CONFIGURATION  -- EDIT THESE FOR YOUR SESSION
%  ========================================================================

INPUT_FILE = 'C:\PATH\TO\YOUR\CellularDataProcessed\01-01-25-R1\Animal01_01-01-25-R1_CellularAnalysis.mat';
OUTPUT_FILE = 'C:\PATH\TO\YOUR\CellularDataProcessed\01-01-25-R1\Animal01_01-01-25-R1_SpikeMLData.mat';

SNIPPET_WINDOW_MS = [-50, 100];  % [pre_peak, post_peak] in milliseconds

%% ========================================================================
%  LOAD SOURCE DATA
%  ========================================================================

fprintf('Loading CellularAnalysis from:\n  %s\n', INPUT_FILE);
S = load(INPUT_FILE, 'CellularAnalysis');
CA = S.CellularAnalysis;

num_trials = CA.metadata.num_trials;
num_neurons = CA.metadata.num_neurons;

fprintf('  Neurons: %d, Trials: %d\n', num_neurons, num_trials);

%% ========================================================================
%  BUILD SpikeMLData
%  ========================================================================

SpikeMLData = struct();

% --- INFO (one-time metadata) -------------------------------------------
info = struct();
info.mouse_name          = CA.metadata.mouse_name;
info.recording_date      = CA.metadata.recording_date;
info.recording_id        = CA.metadata.recording_id;
info.num_neurons         = num_neurons;
info.num_trials          = num_trials;
info.snippet_window_ms   = SNIPPET_WINDOW_MS;
info.source_file         = INPUT_FILE;
info.creation_date       = datestr(now, 'yyyy-mm-dd HH:MM:SS');

if isfield(CA.metadata, 'indicator_promoter')
    info.indicator          = CA.metadata.indicator_promoter;
    info.indicator_polarity = CA.metadata.indicator_polarity;
end

% Spike detection parameters from first non-empty trial
first_trial = CA.trials{1};
info.sampling_rate_hz             = first_trial.parameters.imaging_fs;
info.recording_duration_sec       = first_trial.parameters.recording_duration_sec;
info.spike_detection_method       = first_trial.spikes.detection_type;
info.spike_detection_thresholds   = struct( ...
    'up',   first_trial.parameters.spike_detection_params.up_threshold, ...
    'down', first_trial.parameters.spike_detection_params.down_threshold);

fs = info.sampling_rate_hz;

% Convert snippet window from ms to samples
pre_samples  = round(abs(SNIPPET_WINDOW_MS(1)) * fs / 1000);
post_samples = round(SNIPPET_WINDOW_MS(2) * fs / 1000);
snippet_length = pre_samples + post_samples + 1;  % +1 for the peak sample

info.snippet_window_samples = [-pre_samples, post_samples];
info.snippet_length         = snippet_length;

fprintf('  Sampling rate: %.2f Hz\n', fs);
fprintf('  Snippet window: [%d, +%d] samples (%d total = %.1f ms)\n', ...
    -pre_samples, post_samples, snippet_length, 1000*snippet_length/fs);

% --- Per-trial extraction -----------------------------------------------
all_pooled_snippets   = [];
all_pooled_neuron_ids = [];
all_pooled_times      = [];
all_pooled_snr        = [];
all_pooled_trial_ids  = [];

trials_out = cell(num_trials, 1);
total_spikes_all = 0;

for t = 1:num_trials
    trial = CA.trials{t};
    if isempty(trial), continue; end

    trial_fs = trial.parameters.imaging_fs;
    num_frames = trial.parameters.num_frames;
    corrected_traces = trial.signals.fluorescence_corrected;   % (frames x neurons)
    denoised_traces  = trial.signals.denoised_traces;          % (frames x neurons)
    time_vec         = trial.time.time_vector;

    fprintf('\n--- Trial %d ---\n', t);
    fprintf('  Frames: %d, Duration: %.2f s\n', num_frames, num_frames/trial_fs);

    % -- Traces (keep full traces for this trial) -------------------------
    trial_out = struct();
    trial_out.traces.corrected   = corrected_traces;
    trial_out.traces.denoised    = denoised_traces;
    trial_out.traces.time_vector = time_vec;

    % -- Per-neuron spike extraction --------------------------------------
    neurons_out = struct([]);

    for n = 1:num_neurons
        ndata = struct();
        ndata.neuron_id     = n;
        ndata.firing_rate_hz = trial.spikes.firing_rates_hz(n);
        ndata.noise_level    = trial.spikes.neuron_noise_levels(n);

        spike_idx = trial.spikes.spike_indices{n};
        if isempty(spike_idx)
            spike_idx = [];
        end

        % Discard spikes too close to edges for a full snippet
        valid = spike_idx > pre_samples & spike_idx <= (num_frames - post_samples);
        spike_idx = spike_idx(valid);

        ndata.num_spikes     = length(spike_idx);
        ndata.spike_indices  = spike_idx(:)';
        ndata.spike_times_sec = time_vec(spike_idx)';

        % SNR values (keep only those matching valid spikes)
        raw_snr = trial.spikes.spike_snr_values{n};
        if ~isempty(raw_snr) && length(raw_snr) >= length(spike_idx)
            ndata.spike_snr = raw_snr(valid)';
        else
            ndata.spike_snr = NaN(1, ndata.num_spikes);
        end

        % Extract wide snippets from corrected trace
        trace = corrected_traces(:, n);
        snippets = zeros(ndata.num_spikes, snippet_length);
        for k = 1:ndata.num_spikes
            peak = spike_idx(k);
            snippets(k, :) = trace(peak - pre_samples : peak + post_samples);
        end
        ndata.spike_snippets = snippets;

        % Mean waveform from the wide window
        if ndata.num_spikes > 0
            ndata.mean_waveform = mean(snippets, 1);
        else
            ndata.mean_waveform = NaN(1, snippet_length);
        end

        if isempty(neurons_out)
            neurons_out = ndata;
        else
            neurons_out(end+1) = ndata; %#ok<SAGROW>
        end

        % Accumulate into pooled arrays
        if ndata.num_spikes > 0
            all_pooled_snippets   = [all_pooled_snippets;   snippets]; %#ok<AGROW>
            all_pooled_neuron_ids = [all_pooled_neuron_ids;  repmat(n, ndata.num_spikes, 1)]; %#ok<AGROW>
            all_pooled_times      = [all_pooled_times;       ndata.spike_times_sec(:)]; %#ok<AGROW>
            all_pooled_snr        = [all_pooled_snr;         ndata.spike_snr(:)]; %#ok<AGROW>
            all_pooled_trial_ids  = [all_pooled_trial_ids;   repmat(t, ndata.num_spikes, 1)]; %#ok<AGROW>
        end

        total_spikes_all = total_spikes_all + ndata.num_spikes;

        fprintf('  Neuron %2d: %3d spikes (%.1f Hz)\n', n, ndata.num_spikes, ndata.firing_rate_hz);
    end

    trial_out.neurons = neurons_out;
    trials_out{t} = trial_out;
end

% --- Pooled snippets (ML-ready) -----------------------------------------
snippets_pooled = struct();
snippets_pooled.data         = all_pooled_snippets;       % (total_spikes x snippet_length)
snippets_pooled.neuron_id    = all_pooled_neuron_ids;     % (total_spikes x 1)
snippets_pooled.spike_time_sec = all_pooled_times;        % (total_spikes x 1)
snippets_pooled.snr          = all_pooled_snr;            % (total_spikes x 1)
snippets_pooled.trial_id     = all_pooled_trial_ids;      % (total_spikes x 1)

info.total_spikes = total_spikes_all;

% --- Snippet time axis (relative to peak at 0) --------------------------
snippet_time_ms = 1000 * (-pre_samples : post_samples) / fs;
info.snippet_time_axis_ms = snippet_time_ms;

% --- Assemble final struct -----------------------------------------------
SpikeMLData.info            = info;
SpikeMLData.trials          = trials_out;
SpikeMLData.snippets_pooled = snippets_pooled;

%% ========================================================================
%  SAVE
%  ========================================================================

% Ensure output directory exists
[out_dir, ~, ~] = fileparts(OUTPUT_FILE);
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

fprintf('\nSaving SpikeMLData to:\n  %s\n', OUTPUT_FILE);
save(OUTPUT_FILE, 'SpikeMLData', '-v7.3');

% Report file size
f_info = dir(OUTPUT_FILE);
fprintf('  File size: %.2f MB\n', f_info.bytes / 1e6);

%% ========================================================================
%  SUMMARY
%  ========================================================================

fprintf('\n=== SpikeMLData Summary ===\n');
fprintf('  Mouse: %s  |  Date: %s  |  Rec: %s\n', info.mouse_name, info.recording_date, info.recording_id);
fprintf('  Sampling rate: %.2f Hz\n', info.sampling_rate_hz);
fprintf('  Snippet window: [%.0f, +%.0f] ms  (%d samples)\n', ...
    SNIPPET_WINDOW_MS(1), SNIPPET_WINDOW_MS(2), snippet_length);
fprintf('  Neurons: %d\n', info.num_neurons);
fprintf('  Total spikes: %d\n', info.total_spikes);
fprintf('  Pooled matrix size: [%d x %d]\n', size(snippets_pooled.data));
fprintf('  Output: %s\n', OUTPUT_FILE);
fprintf('===========================\n');
