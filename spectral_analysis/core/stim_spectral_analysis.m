function results = stim_spectral_analysis(animal, session, trial_num, method, cfg)
%% ============================================================================
%  STIMULATION TRIAL SPECTRAL ANALYSIS
%  ============================================================================
%  Compute spectral metrics for a single stimulation trial.
%  This function analyzes ONE trial and returns comprehensive spectral results
%  for Python plotting.
%
%  INPUTS:
%    animal    - Animal struct from stim_database
%    session   - Session struct with stimulation parameters
%    trial_num - Trial number (1-based)
%    method    - 'mscohere' or 'fieldtrip'
%    cfg       - Configuration from stim_analysis_config()
%
%  OUTPUTS:
%    results   - Struct containing all spectral metrics
%  ============================================================================

results = struct();
results.success = false;

%% ============================================================================
%  BUILD FILE PATH AND LOAD DATA
%  ============================================================================

% Construct trial folder path
session_path = fullfile(animal.data_root, session.session_id);

% Find trial folder (handle pattern matching)
% First try the pattern from config
if isfield(session, 'trial_folder_pattern') && ~isempty(session.trial_folder_pattern)
    try
        trial_folder_name = sprintf(session.trial_folder_pattern, trial_num);
        trial_folder_path = fullfile(session_path, trial_folder_name);
    catch
        trial_folder_path = '';
    end
else
    trial_folder_path = '';
end

% If exact pattern doesn't exist, try glob search
if isempty(trial_folder_path) || ~exist(trial_folder_path, 'dir')
    % Try searching for Trial{N}_* pattern
    search_pattern = fullfile(session_path, sprintf('Trial%d_*', trial_num));
    found_dirs = dir(search_pattern);
    if ~isempty(found_dirs) && found_dirs(1).isdir
        trial_folder_path = fullfile(session_path, found_dirs(1).name);
        fprintf('    Found trial folder: %s\n', found_dirs(1).name);
    else
        warning('Trial folder not found for trial %d in: %s', trial_num, session_path);
        fprintf('    Searched for: %s\n', search_pattern);
        return;
    end
end

% Construct MAT file path
if isfield(session, 'mat_file_pattern') && ~isempty(session.mat_file_pattern)
    try
        mat_filename = sprintf(session.mat_file_pattern, session.session_id, trial_num);
    catch
        % Fallback pattern
        mat_filename = sprintf('%s-%s_Trial%d_FiberPhotometry_Analysis.mat', ...
            animal.mouse_id, session.session_id, trial_num);
    end
else
    % Default pattern
    mat_filename = sprintf('%s-%s_Trial%d_FiberPhotometry_Analysis.mat', ...
        animal.mouse_id, session.session_id, trial_num);
end
mat_file_path = fullfile(trial_folder_path, mat_filename);

if ~exist(mat_file_path, 'file')
    warning('MAT file not found: %s', mat_file_path);
    fprintf('    Expected file: %s\n', mat_filename);
    fprintf('    In folder: %s\n', trial_folder_path);
    return;
end

fprintf('    Loading: %s\n', mat_file_path);

% Load trial data
data = load_stim_trial_data(mat_file_path, cfg);
if isempty(data)
    warning('Failed to load data from: %s', mat_file_path);
    return;
end

%% ============================================================================
%  EXTRACT TIMING INFORMATION
%  ============================================================================

% Determine stimulation onset
% First try to get from data struct, else use pre_stim_duration as assumption
if isfield(data, 'stim_onset_sec') && ~isempty(data.stim_onset_sec)
    stim_onset_sec = data.stim_onset_sec;
else
    % Assume stim starts at pre_stim_duration from recording start
    stim_onset_sec = cfg.stim_timing.pre_stim_duration_sec;
end

% Calculate period boundaries (relative to recording start)
% For coherence analysis: use shorter, balanced windows for fair comparison
% Pre-stim: use last N seconds before stim onset (default: 1s)
% Post-stim: use first N seconds after stim offset (default: 1s)
if isfield(cfg.stim_timing, 'coherence_prestim_window_sec')
    pre_stim_window = cfg.stim_timing.coherence_prestim_window_sec;
else
    pre_stim_window = 1.0;  % Default: 1 second
end
if isfield(cfg.stim_timing, 'coherence_poststim_window_sec')
    post_stim_window = cfg.stim_timing.coherence_poststim_window_sec;
else
    post_stim_window = 1.0;  % Default: 1 second
end

% Pre-stim: last N seconds before stim onset
pre_stim_start = max(0, stim_onset_sec - pre_stim_window);
pre_stim_end = stim_onset_sec;
stim_start = stim_onset_sec;
stim_end = stim_onset_sec + cfg.stim_timing.stim_duration_sec;
% Post-stim: first N seconds after stim offset
post_stim_start = stim_end;
post_stim_end = stim_end + post_stim_window;

% Sub-periods within stimulation
transient_start = stim_start;
transient_end = stim_start + cfg.stim_timing.transient_end_sec;
sustained_start = stim_start + cfg.stim_timing.sustained_start_sec;
sustained_end = stim_end;

%% ============================================================================
%  STORE METADATA
%  ============================================================================

results.metadata = struct();
results.metadata.animal = animal.mouse_id;
results.metadata.project = animal.project;
results.metadata.session = session.session_id;
results.metadata.trial_num = trial_num;
results.metadata.condition = session.condition_name;
results.metadata.method = method;
results.metadata.mat_file = mat_file_path;
results.metadata.fs = data.fs;

% Stimulation parameters
results.stim_params = session.stim_params;

% Timing information (in seconds, relative to recording start)
results.timing = struct();
results.timing.stim_onset_sec = stim_onset_sec;
results.timing.pre_stim = [pre_stim_start, pre_stim_end];
results.timing.stim_full = [stim_start, stim_end];
results.timing.stim_transient = [transient_start, transient_end];
results.timing.stim_sustained = [sustained_start, sustained_end];
results.timing.post_stim = [post_stim_start, post_stim_end];

% Also store duration parameters for Python (matching config)
results.timing.pre_stim_duration_sec = cfg.stim_timing.pre_stim_duration_sec;
results.timing.stim_duration_sec = cfg.stim_timing.stim_duration_sec;
results.timing.post_stim_duration_sec = cfg.stim_timing.post_stim_duration_sec;
results.timing.transient_end_sec = cfg.stim_timing.transient_end_sec;
results.timing.sustained_start_sec = cfg.stim_timing.sustained_start_sec;
% Coherence window parameters (for balanced period comparison)
results.timing.coherence_prestim_window_sec = pre_stim_window;
results.timing.coherence_poststim_window_sec = post_stim_window;

% Recording info
results.timing.recording_duration_sec = data.time(end) - data.time(1);

%% ============================================================================
%  SEGMENT DATA INTO PERIODS
%  ============================================================================

fs = data.fs;
t = data.time;

% Create logical masks for each period
mask_pre = (t >= pre_stim_start) & (t < pre_stim_end);
mask_stim_full = (t >= stim_start) & (t < stim_end);
mask_transient = (t >= transient_start) & (t < transient_end);
mask_sustained = (t >= sustained_start) & (t < sustained_end);
mask_post = (t >= post_stim_start) & (t <= min(post_stim_end, t(end)));

% Store sample counts for each period (initialize structure first)
results.periods = struct();
results.periods.pre_stim = struct();
results.periods.stim_full = struct();
results.periods.transient = struct();
results.periods.sustained = struct();
results.periods.post_stim = struct();

results.periods.pre_stim.n_samples = sum(mask_pre);
results.periods.stim_full.n_samples = sum(mask_stim_full);
results.periods.transient.n_samples = sum(mask_transient);
results.periods.sustained.n_samples = sum(mask_sustained);
results.periods.post_stim.n_samples = sum(mask_post);

% Print period info for debugging
fprintf('      Period samples: pre=%d, trans=%d, sust=%d, post=%d\n', ...
    results.periods.pre_stim.n_samples, results.periods.transient.n_samples, ...
    results.periods.sustained.n_samples, results.periods.post_stim.n_samples);

%% ============================================================================
%  COMPUTE SPECTRAL METRICS - DISPATCH BY METHOD
%  ============================================================================

fprintf('      Computing spectral metrics (%s)...\n', method);
try
    switch lower(method)
        case 'mscohere'
            spectra = compute_stim_spectra_mscohere(data, mask_pre, mask_transient, ...
                mask_sustained, mask_stim_full, mask_post, cfg);
        case 'fieldtrip'
            spectra = compute_stim_spectra_fieldtrip(data, mask_pre, mask_transient, ...
                mask_sustained, mask_stim_full, mask_post, cfg);
        otherwise
            error('Unknown method: %s', method);
    end
    
    % Validate spectra structure
    if ~isfield(spectra, 'coherence') || ~isfield(spectra, 'psd_lfp') || ~isfield(spectra, 'psd_fiber')
        error('Spectral computation returned incomplete structure');
    end
    
    % Store spectral results
    results.coherence = spectra.coherence;
    results.psd_lfp = spectra.psd_lfp;
    results.psd_fiber = spectra.psd_fiber;
    results.freq = spectra.freq;
catch ME
    warning('Spectral computation failed: %s', ME.message);
    fprintf('      Error details: %s\n', getReport(ME));
    return;
end

%% ============================================================================
%  COMPUTE TIME-RESOLVED COHERENCE (for heatmaps)
%  ============================================================================

fprintf('      Computing time-resolved coherence...\n');
try
    time_resolved = compute_time_resolved_coherence(data, method, cfg);
    if isempty(time_resolved) || (~isfield(time_resolved, 'coherence') || isempty(time_resolved.coherence))
        warning('Time-resolved coherence computation returned empty result');
        time_resolved = struct();
        time_resolved.coherence = [];
        time_resolved.time = [];
        time_resolved.freq = [];
    end
    results.time_resolved = time_resolved;
catch ME
    warning('Time-resolved coherence computation failed: %s', ME.message);
    results.time_resolved = struct();
    results.time_resolved.coherence = [];
    results.time_resolved.time = [];
    results.time_resolved.freq = [];
end

%% ============================================================================
%  COMPUTE SPECTROGRAMS (for heatmaps)
%  ============================================================================

fprintf('      Computing spectrograms...\n');
try
    [spec_lfp, spec_fiber] = compute_spectrograms(data, cfg);
    if isempty(spec_lfp) || ~isfield(spec_lfp, 'power')
        warning('LFP spectrogram computation returned empty result');
        spec_lfp = struct();
        spec_lfp.power = [];
        spec_lfp.freq = [];
        spec_lfp.time = [];
    end
    if isempty(spec_fiber) || ~isfield(spec_fiber, 'power')
        warning('Fiber spectrogram computation returned empty result');
        spec_fiber = struct();
        spec_fiber.power = [];
        spec_fiber.freq = [];
        spec_fiber.time = [];
    end
    results.spectrogram_lfp = spec_lfp;
    results.spectrogram_fiber = spec_fiber;
catch ME
    warning('Spectrogram computation failed: %s', ME.message);
    results.spectrogram_lfp = struct();
    results.spectrogram_lfp.power = [];
    results.spectrogram_lfp.freq = [];
    results.spectrogram_lfp.time = [];
    results.spectrogram_fiber = struct();
    results.spectrogram_fiber.power = [];
    results.spectrogram_fiber.freq = [];
    results.spectrogram_fiber.time = [];
end

%% ============================================================================
%  STORE RAW TRACES (for Python plotting)
%  ============================================================================

% Always initialize traces structure (Python expects it)
results.traces = struct();
if cfg.output.save_raw_traces
    % Convert time to be relative to stim onset for easier plotting
    results.traces.time_sec = data.time - stim_onset_sec;  % 0 = stim onset
    results.traces.time_abs_sec = data.time;               % Original time
    results.traces.lfp = data.lfp;
    results.traces.fiber = data.gevi;
    results.traces.motion = data.speed;
    results.traces.fs = data.fs;
else
    % Still provide empty structure for consistency
    results.traces.time_sec = [];
    results.traces.time_abs_sec = [];
    results.traces.lfp = [];
    results.traces.fiber = [];
    results.traces.motion = [];
    results.traces.fs = data.fs;
end

results.success = true;
fprintf('      Trial %d completed successfully.\n', trial_num);

end

%% ============================================================================
%  HELPER: Load Stimulation Trial Data
%  ============================================================================
function data = load_stim_trial_data(mat_file_path, cfg)
%LOAD_STIM_TRIAL_DATA Load data from preprocessed MAT file

data = [];

try
    mat_data = load(mat_file_path);
    
    % Handle different struct names
    if isfield(mat_data, 'FiberPhotometryAnalysis')
        FPA = mat_data.FiberPhotometryAnalysis;
    elseif isfield(mat_data, 'FPA')
        FPA = mat_data.FPA;
    else
        fn = fieldnames(mat_data);
        if ~isempty(fn)
            FPA = mat_data.(fn{1});
        else
            warning('No data struct found in %s', mat_file_path);
            return;
        end
    end
    
    % Extract time vector
    if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
        t = FPA.time.time_vector_seconds(:);
    else
        warning('Time vector not found');
        return;
    end
    
    % Extract sampling rate
    if isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
        fs = FPA.time.sampling_rate;
    else
        fs = 1 / median(diff(t));
    end
    
    % Extract fiber trace
    fiber_index = cfg.fiber_index;
    if isfield(FPA, 'signals') && isfield(FPA.signals, 'final_processed_traces')
        fiber_all = FPA.signals.final_processed_traces;
        if size(fiber_all, 2) >= fiber_index
            fiber_trace = fiber_all(:, fiber_index);
        else
            fiber_trace = fiber_all(:, 1);
        end
    else
        warning('Fiber trace not found');
        return;
    end
    
    % Extract LFP trace
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_HP')
        lfp_trace = FPA.ephys.lfp_raw_aligned_HP(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_mPFC')
        lfp_trace = FPA.ephys.lfp_raw_aligned_mPFC(:);
    else
        warning('LFP trace not found');
        return;
    end
    
    % Extract motion trace
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
        motion_raw = FPA.ephys.running_velocity_smooth(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
        motion_raw = FPA.ephys.running_velocity(:);
    else
        warning('Motion trace not found');
        motion_raw = zeros(size(t));
    end
    
    % Try to get stimulation onset from data
    stim_onset_sec = [];
    if isfield(FPA, 'stimulation') && isfield(FPA.stimulation, 'onset_time_sec')
        stim_onset_sec = FPA.stimulation.onset_time_sec;
    elseif isfield(FPA, 'STIM_PERIOD') && ~isempty(FPA.STIM_PERIOD)
        stim_onset_sec = FPA.STIM_PERIOD(1);
    end
    
    % Align lengths
    n = min([length(t), length(fiber_trace), length(lfp_trace), length(motion_raw)]);
    t = t(1:n);
    fiber_trace = fiber_trace(1:n);
    lfp_trace = lfp_trace(1:n);
    motion_raw = motion_raw(1:n);
    
    % Convert motion to cm/s
    wheel_circumference = pi * cfg.motion.wheel_diameter_cm;
    distance_per_edge = wheel_circumference / cfg.motion.encoder_counts_per_rev;
    motion_to_cms = (cfg.motion.ephys_sampling_rate / 1000) * distance_per_edge;
    speed_cm_s = motion_raw * motion_to_cms;
    
    % Optional smoothing
    if cfg.motion.smooth_samples > 1
        kernel = ones(cfg.motion.smooth_samples, 1) / cfg.motion.smooth_samples;
        speed_cm_s = conv(speed_cm_s, kernel, 'same');
    end
    
    % Store
    data = struct();
    data.time = t;
    data.fs = fs;
    data.gevi = fiber_trace;
    data.lfp = lfp_trace;
    data.speed = speed_cm_s;
    data.stim_onset_sec = stim_onset_sec;
    
catch ME
    warning('Failed to load %s: %s', mat_file_path, ME.message);
    data = [];
end

end

%% ============================================================================
%  HELPER: Compute Spectra (mscohere method)
%  ============================================================================
function spectra = compute_stim_spectra_mscohere(data, mask_pre, mask_trans, mask_sust, mask_stim_full, mask_post, cfg)
%COMPUTE_STIM_SPECTRA_MSCOHERE Compute coherence and PSD using mscohere

fs = data.fs;
lfp = data.lfp;
gevi = data.gevi;

% Parameters for OVERALL coherence (full trial - uses standard segments)
segment_samples = round(cfg.coherence.mscohere.segment_sec * fs);
overlap_samples = round(cfg.coherence.mscohere.overlap_frac * segment_samples);
nfft = cfg.coherence.mscohere.nfft_factor * segment_samples;
freq_min = cfg.coherence.mscohere.freq_min;
freq_max = cfg.coherence.mscohere.freq_max;
window = hanning(segment_samples);

% Parameters for PERIOD coherence (shorter segments for 1s windows)
if isfield(cfg.coherence, 'period_segment_sec')
    period_segment_samples = round(cfg.coherence.period_segment_sec * fs);
    period_overlap_frac = cfg.coherence.period_overlap_frac;
else
    % Default: use 0.5s segments for period coherence
    period_segment_samples = round(0.5 * fs);
    period_overlap_frac = 0.5;
end
period_overlap_samples = round(period_overlap_frac * period_segment_samples);
period_nfft = cfg.coherence.mscohere.nfft_factor * period_segment_samples;
period_window = hanning(period_segment_samples);

spectra = struct();

% -------------------------------------------------------------------------
% OVERALL (full trial)
% -------------------------------------------------------------------------
[Cxy, F] = mscohere(lfp, gevi, window, overlap_samples, nfft, fs);
freq_mask = (F >= freq_min) & (F <= freq_max);
spectra.freq = F(freq_mask);
spectra.coherence.overall = Cxy(freq_mask);

% PSD
[Plfp, ~] = pwelch(lfp, window, overlap_samples, nfft, fs);
[Pgevi, ~] = pwelch(gevi, window, overlap_samples, nfft, fs);
spectra.psd_lfp.overall = 10*log10(Plfp(freq_mask) + 1e-10);
spectra.psd_fiber.overall = 10*log10(Pgevi(freq_mask) + 1e-10);

% -------------------------------------------------------------------------
% PRE-STIM (uses period-specific shorter segments for 1s window compatibility)
% -------------------------------------------------------------------------
if sum(mask_pre) >= period_segment_samples * 2
    [Cxy_pre, F_period] = mscohere(lfp(mask_pre), gevi(mask_pre), period_window, period_overlap_samples, period_nfft, fs);
    freq_mask_period = (F_period >= freq_min) & (F_period <= freq_max);
    % Interpolate to match overall frequency grid
    spectra.coherence.pre_stim = interp1(F_period(freq_mask_period), Cxy_pre(freq_mask_period), spectra.freq, 'linear', 'extrap');
    spectra.coherence.pre_stim = max(0, min(1, spectra.coherence.pre_stim));  % Clamp to [0,1]
    
    [Plfp_pre, ~] = pwelch(lfp(mask_pre), period_window, period_overlap_samples, period_nfft, fs);
    [Pgevi_pre, ~] = pwelch(gevi(mask_pre), period_window, period_overlap_samples, period_nfft, fs);
    spectra.psd_lfp.pre_stim = interp1(F_period(freq_mask_period), 10*log10(Plfp_pre(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
    spectra.psd_fiber.pre_stim = interp1(F_period(freq_mask_period), 10*log10(Pgevi_pre(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
else
    spectra.coherence.pre_stim = nan(size(spectra.freq));
    spectra.psd_lfp.pre_stim = nan(size(spectra.freq));
    spectra.psd_fiber.pre_stim = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% TRANSIENT (0-0.15s) - usually too short for reliable coherence
% -------------------------------------------------------------------------
min_samples_transient = max(32, round(fs * 0.1));  % At least 100ms
if sum(mask_trans) >= min_samples_transient
    % Use shorter window for transient period
    trans_window_samples = min(segment_samples, sum(mask_trans));
    trans_window = hanning(trans_window_samples);
    trans_overlap = round(0.5 * trans_window_samples);
    trans_nfft = 2^nextpow2(trans_window_samples * 2);
    
    try
        [Cxy_trans, F_trans] = mscohere(lfp(mask_trans), gevi(mask_trans), ...
            trans_window, trans_overlap, trans_nfft, fs);
        freq_mask_trans = (F_trans >= freq_min) & (F_trans <= freq_max);
        
        % Interpolate to common frequency grid
        spectra.coherence.transient = interp1(F_trans(freq_mask_trans), ...
            Cxy_trans(freq_mask_trans), spectra.freq, 'linear', 'extrap');
        
        [Plfp_trans, ~] = pwelch(lfp(mask_trans), trans_window, trans_overlap, trans_nfft, fs);
        [Pgevi_trans, ~] = pwelch(gevi(mask_trans), trans_window, trans_overlap, trans_nfft, fs);
        spectra.psd_lfp.transient = interp1(F_trans(freq_mask_trans), ...
            10*log10(Plfp_trans(freq_mask_trans) + 1e-10), spectra.freq, 'linear', 'extrap');
        spectra.psd_fiber.transient = interp1(F_trans(freq_mask_trans), ...
            10*log10(Pgevi_trans(freq_mask_trans) + 1e-10), spectra.freq, 'linear', 'extrap');
    catch
        spectra.coherence.transient = nan(size(spectra.freq));
        spectra.psd_lfp.transient = nan(size(spectra.freq));
        spectra.psd_fiber.transient = nan(size(spectra.freq));
    end
else
    spectra.coherence.transient = nan(size(spectra.freq));
    spectra.psd_lfp.transient = nan(size(spectra.freq));
    spectra.psd_fiber.transient = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% SUSTAINED (0.15-1.0s) - uses period-specific shorter segments
% NOTE: For period comparison, we use FULL STIM (1.0s) to match pre/post window sizes
% -------------------------------------------------------------------------
% Use full stim period (mask_stim_full) for consistent 1s comparison
if sum(mask_stim_full) >= period_segment_samples * 2
    [Cxy_sust, F_period] = mscohere(lfp(mask_stim_full), gevi(mask_stim_full), period_window, period_overlap_samples, period_nfft, fs);
    freq_mask_period = (F_period >= freq_min) & (F_period <= freq_max);
    spectra.coherence.sustained = interp1(F_period(freq_mask_period), Cxy_sust(freq_mask_period), spectra.freq, 'linear', 'extrap');
    spectra.coherence.sustained = max(0, min(1, spectra.coherence.sustained));
    
    [Plfp_sust, ~] = pwelch(lfp(mask_stim_full), period_window, period_overlap_samples, period_nfft, fs);
    [Pgevi_sust, ~] = pwelch(gevi(mask_stim_full), period_window, period_overlap_samples, period_nfft, fs);
    spectra.psd_lfp.sustained = interp1(F_period(freq_mask_period), 10*log10(Plfp_sust(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
    spectra.psd_fiber.sustained = interp1(F_period(freq_mask_period), 10*log10(Pgevi_sust(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
else
    spectra.coherence.sustained = nan(size(spectra.freq));
    spectra.psd_lfp.sustained = nan(size(spectra.freq));
    spectra.psd_fiber.sustained = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% POST-STIM (uses period-specific shorter segments for 1s window compatibility)
% -------------------------------------------------------------------------
if sum(mask_post) >= period_segment_samples * 2
    [Cxy_post, F_period] = mscohere(lfp(mask_post), gevi(mask_post), period_window, period_overlap_samples, period_nfft, fs);
    freq_mask_period = (F_period >= freq_min) & (F_period <= freq_max);
    spectra.coherence.post_stim = interp1(F_period(freq_mask_period), Cxy_post(freq_mask_period), spectra.freq, 'linear', 'extrap');
    spectra.coherence.post_stim = max(0, min(1, spectra.coherence.post_stim));
    
    [Plfp_post, ~] = pwelch(lfp(mask_post), period_window, period_overlap_samples, period_nfft, fs);
    [Pgevi_post, ~] = pwelch(gevi(mask_post), period_window, period_overlap_samples, period_nfft, fs);
    spectra.psd_lfp.post_stim = interp1(F_period(freq_mask_period), 10*log10(Plfp_post(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
    spectra.psd_fiber.post_stim = interp1(F_period(freq_mask_period), 10*log10(Pgevi_post(freq_mask_period) + 1e-10), spectra.freq, 'linear', 'extrap');
else
    spectra.coherence.post_stim = nan(size(spectra.freq));
    spectra.psd_lfp.post_stim = nan(size(spectra.freq));
    spectra.psd_fiber.post_stim = nan(size(spectra.freq));
end

end

%% ============================================================================
%  HELPER: Compute Spectra (FieldTrip method)
%  ============================================================================
function spectra = compute_stim_spectra_fieldtrip(data, mask_pre, mask_trans, mask_sust, mask_stim_full, mask_post, cfg)
%COMPUTE_STIM_SPECTRA_FIELDTRIP Compute coherence and PSD using FieldTrip

% Check FieldTrip availability
if ~exist('ft_freqanalysis', 'file')
    warning('FieldTrip not available, falling back to mscohere');
    spectra = compute_stim_spectra_mscohere(data, mask_pre, mask_trans, mask_sust, mask_post, cfg);
    return;
end

fs = data.fs;
lfp = data.lfp;
gevi = data.gevi;

% Parameters for OVERALL coherence (standard pseudotrials)
EPOCH_LENGTH_SEC = cfg.coherence.fieldtrip.pseudotrial_length_sec;
EPOCH_OVERLAP_SEC = cfg.coherence.fieldtrip.pseudotrial_overlap_sec;

% Parameters for PERIOD coherence (shorter pseudotrials for 1s windows)
if isfield(cfg.coherence.fieldtrip, 'period_pseudotrial_length_sec')
    PERIOD_EPOCH_LENGTH_SEC = cfg.coherence.fieldtrip.period_pseudotrial_length_sec;
    PERIOD_EPOCH_OVERLAP_SEC = cfg.coherence.fieldtrip.period_pseudotrial_overlap_sec;
else
    % Default: 0.5s pseudotrials with 50% overlap
    PERIOD_EPOCH_LENGTH_SEC = 0.5;
    PERIOD_EPOCH_OVERLAP_SEC = 0.25;
end

spectra = struct();
spectra.freq = (cfg.coherence.fieldtrip.foi_min:cfg.coherence.fieldtrip.foi_step:cfg.coherence.fieldtrip.foi_max)';

% FieldTrip configs
ft_cfg = [];
ft_cfg.method = cfg.coherence.fieldtrip.method;
ft_cfg.taper = cfg.coherence.fieldtrip.taper;
ft_cfg.foi = spectra.freq';
ft_cfg.keeptrials = 'yes';
ft_cfg.output = 'fourier';
ft_cfg.pad = 'nextpow2';
if strcmpi(cfg.coherence.fieldtrip.taper, 'dpss')
    ft_cfg.tapsmofrq = cfg.coherence.fieldtrip.tapsmofrq;
end

coh_cfg = [];
coh_cfg.method = 'coh';
coh_cfg.complex = 'abs';

ft_cfg_pow = ft_cfg;
ft_cfg_pow.output = 'pow';
ft_cfg_pow.keeptrials = 'no';

% -------------------------------------------------------------------------
% HELPER: Create pseudo-trials (supports different epoch lengths)
% -------------------------------------------------------------------------
    function ft_data = create_pseudotrials(lfp_seg, gevi_seg, fs_local, epoch_len_sec, epoch_overlap_sec)
        % Use default parameters if not provided
        if nargin < 4
            epoch_len_sec = EPOCH_LENGTH_SEC;
            epoch_overlap_sec = EPOCH_OVERLAP_SEC;
        end
        
        n_samples = length(lfp_seg);
        epoch_samples = round(epoch_len_sec * fs_local);
        epoch_step = round((epoch_len_sec - epoch_overlap_sec) * fs_local);
        n_epochs = floor((n_samples - epoch_samples) / epoch_step) + 1;
        
        if n_epochs < 2
            ft_data = [];
            return;
        end
        
        ft_data = struct();
        ft_data.label = {'LFP'; 'GEVI'};
        ft_data.fsample = fs_local;
        ft_data.trial = cell(1, n_epochs);
        ft_data.time = cell(1, n_epochs);
        
        for ep = 1:n_epochs
            start_idx = (ep - 1) * epoch_step + 1;
            end_idx = start_idx + epoch_samples - 1;
            if end_idx > n_samples
                break;
            end
            ft_data.trial{ep} = [lfp_seg(start_idx:end_idx)'; gevi_seg(start_idx:end_idx)'];
            ft_data.time{ep} = (0:epoch_samples-1) / fs_local;
        end
    end

% -------------------------------------------------------------------------
% OVERALL
% -------------------------------------------------------------------------
ft_data_all = create_pseudotrials(lfp, gevi, fs);
if ~isempty(ft_data_all)
    try
        freq_all = ft_freqanalysis(ft_cfg, ft_data_all);
        coh_all = ft_connectivityanalysis(coh_cfg, freq_all);
        
        % Extract coherence - handle both 2D and 3D cohspctrm
        coh_vals = squeeze(coh_all.cohspctrm(1,2,:));
        if isempty(coh_vals)
            % Try alternative dimension order
            coh_vals = squeeze(coh_all.cohspctrm);
            if size(coh_vals, 1) == 1 && size(coh_vals, 2) > 1
                coh_vals = coh_vals(:);
            end
        end
        
        % Interpolate to match spectra.freq if needed
        if length(coh_vals) == length(coh_all.freq) && length(coh_vals) ~= length(spectra.freq)
            % Interpolate to match the expected frequency axis
            spectra.coherence.overall = interp1(coh_all.freq(:), coh_vals(:), spectra.freq(:), 'linear', 'extrap');
        else
            spectra.coherence.overall = coh_vals(:);
        end
        
        % Ensure correct size
        if length(spectra.coherence.overall) ~= length(spectra.freq)
            warning('FieldTrip coherence size mismatch: coh=%d, freq=%d. Using NaN.', ...
                length(spectra.coherence.overall), length(spectra.freq));
            spectra.coherence.overall = nan(size(spectra.freq));
        end
        
        freq_pow = ft_freqanalysis(ft_cfg_pow, ft_data_all);
        spectra.psd_lfp.overall = 10*log10(squeeze(freq_pow.powspctrm(1,:))' + 1e-10);
        spectra.psd_fiber.overall = 10*log10(squeeze(freq_pow.powspctrm(2,:))' + 1e-10);
    catch ME
        warning('FieldTrip overall failed: %s', ME.message);
        spectra.coherence.overall = nan(size(spectra.freq));
        spectra.psd_lfp.overall = nan(size(spectra.freq));
        spectra.psd_fiber.overall = nan(size(spectra.freq));
    end
else
    spectra.coherence.overall = nan(size(spectra.freq));
    spectra.psd_lfp.overall = nan(size(spectra.freq));
    spectra.psd_fiber.overall = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% PRE-STIM (uses period-specific shorter pseudotrials for 1s window)
% -------------------------------------------------------------------------
if sum(mask_pre) >= round(PERIOD_EPOCH_LENGTH_SEC * fs * 2)
    ft_data_pre = create_pseudotrials(lfp(mask_pre), gevi(mask_pre), fs, PERIOD_EPOCH_LENGTH_SEC, PERIOD_EPOCH_OVERLAP_SEC);
    if ~isempty(ft_data_pre)
        try
            freq_pre = ft_freqanalysis(ft_cfg, ft_data_pre);
            coh_pre = ft_connectivityanalysis(coh_cfg, freq_pre);
            
            % Extract coherence with interpolation if needed
            coh_vals = squeeze(coh_pre.cohspctrm(1,2,:));
            if length(coh_vals) == length(coh_pre.freq) && length(coh_vals) ~= length(spectra.freq)
                spectra.coherence.pre_stim = interp1(coh_pre.freq(:), coh_vals(:), spectra.freq(:), 'linear', 'extrap');
            else
                spectra.coherence.pre_stim = coh_vals(:);
            end
            if length(spectra.coherence.pre_stim) ~= length(spectra.freq)
                spectra.coherence.pre_stim = nan(size(spectra.freq));
            end
            
            freq_pow_pre = ft_freqanalysis(ft_cfg_pow, ft_data_pre);
            spectra.psd_lfp.pre_stim = 10*log10(squeeze(freq_pow_pre.powspctrm(1,:))' + 1e-10);
            spectra.psd_fiber.pre_stim = 10*log10(squeeze(freq_pow_pre.powspctrm(2,:))' + 1e-10);
        catch
            spectra.coherence.pre_stim = nan(size(spectra.freq));
            spectra.psd_lfp.pre_stim = nan(size(spectra.freq));
            spectra.psd_fiber.pre_stim = nan(size(spectra.freq));
        end
    else
        spectra.coherence.pre_stim = nan(size(spectra.freq));
        spectra.psd_lfp.pre_stim = nan(size(spectra.freq));
        spectra.psd_fiber.pre_stim = nan(size(spectra.freq));
    end
else
    spectra.coherence.pre_stim = nan(size(spectra.freq));
    spectra.psd_lfp.pre_stim = nan(size(spectra.freq));
    spectra.psd_fiber.pre_stim = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% TRANSIENT - use mscohere fallback (too short for FieldTrip)
% -------------------------------------------------------------------------
spectra.coherence.transient = nan(size(spectra.freq));
spectra.psd_lfp.transient = nan(size(spectra.freq));
spectra.psd_fiber.transient = nan(size(spectra.freq));

if sum(mask_trans) >= 32
    try
        lfp_trans = lfp(mask_trans);
        gevi_trans = gevi(mask_trans);
        n_trans = length(lfp_trans);
        
        % For transient period (0.15s = 75 samples at 500Hz), we need a window
        % that's smaller than the data to get multiple segments
        % Use window that's at most 50% of data length to ensure at least 2 segments
        trans_seg = round(min(n_trans * 0.5, fs * 0.1));  % Max 100ms or 50% of data
        % Ensure minimum window size for frequency resolution
        if trans_seg < 32
            trans_seg = 32;  % Minimum 32 samples for basic frequency resolution
        end
        % Ensure window is smaller than data to get multiple segments
        if trans_seg >= n_trans
            trans_seg = max(16, floor(n_trans * 0.4));  % Use 40% of data, min 16 samples
        end
        
        trans_ovlp = round(0.5 * trans_seg);
        % Ensure overlap is less than segment length
        if trans_ovlp >= trans_seg
            trans_ovlp = max(1, round(0.3 * trans_seg));  % Use 30% overlap if 50% too large
        end
        
        % Check we can get at least 2 segments
        n_segments = floor((n_trans - trans_ovlp) / (trans_seg - trans_ovlp));
        if n_segments < 2 && n_trans > trans_seg
            % Reduce window size to get more segments
            trans_seg = floor(n_trans / 3);  % Aim for 3 segments
            trans_ovlp = round(0.3 * trans_seg);
        end
        
        trans_nfft = 2^nextpow2(trans_seg * 2);
        trans_win = hanning(trans_seg);
        
        % Check if we have enough data and can get multiple segments
        if length(lfp_trans) >= trans_seg && length(gevi_trans) >= trans_seg && n_segments >= 1
            [Cxy, F] = mscohere(lfp_trans, gevi_trans, trans_win, trans_ovlp, trans_nfft, fs);
            
            % Validate coherence values (should be between 0 and 1)
            % Check if all values are exactly 1.0 (indicates unreliable estimate - window too large)
            if all(Cxy == 1.0) || (all(isfinite(Cxy)) && all(abs(Cxy - 1.0) < 1e-10))
                warning('Transient coherence: All values are 1.0 - likely unreliable (window=%d, data=%d samples, segments=%d)', ...
                    trans_seg, n_trans, n_segments);
                % Use NaN instead of perfect coherence
                spectra.coherence.transient = nan(size(spectra.freq));
            else
                Cxy = max(0, min(1, Cxy));  % Clamp to [0, 1]
                
                % Interpolate only if we have valid data
                if ~isempty(Cxy) && ~isempty(F) && any(isfinite(Cxy)) && any(isfinite(F))
                    spectra.coherence.transient = interp1(F, Cxy, spectra.freq, 'linear', 'extrap');
                    % Clamp interpolated values to [0, 1]
                    spectra.coherence.transient = max(0, min(1, spectra.coherence.transient));
                else
                    spectra.coherence.transient = nan(size(spectra.freq));
                end
            end
            
            [Plfp, ~] = pwelch(lfp_trans, trans_win, trans_ovlp, trans_nfft, fs);
            [Pgevi, ~] = pwelch(gevi_trans, trans_win, trans_ovlp, trans_nfft, fs);
            
            if ~isempty(Plfp) && ~isempty(F) && any(isfinite(Plfp))
                spectra.psd_lfp.transient = interp1(F, 10*log10(Plfp + 1e-10), spectra.freq, 'linear', 'extrap');
            else
                spectra.psd_lfp.transient = nan(size(spectra.freq));
            end
            
            if ~isempty(Pgevi) && ~isempty(F) && any(isfinite(Pgevi))
                spectra.psd_fiber.transient = interp1(F, 10*log10(Pgevi + 1e-10), spectra.freq, 'linear', 'extrap');
            else
                spectra.psd_fiber.transient = nan(size(spectra.freq));
            end
        else
            % Not enough data
            spectra.coherence.transient = nan(size(spectra.freq));
            spectra.psd_lfp.transient = nan(size(spectra.freq));
            spectra.psd_fiber.transient = nan(size(spectra.freq));
        end
    catch ME
        warning('Transient period computation failed: %s', ME.message);
        % Keep NaN
    end
end

% -------------------------------------------------------------------------
% SUSTAINED - Uses FULL STIM period (0-1.0s) for fair comparison with pre/post
% -------------------------------------------------------------------------
% Use mscohere with period-specific parameters (0.5s segments) to ensure
% compatibility with 1s pre-stim and post-stim windows
spectra.coherence.sustained = nan(size(spectra.freq));
spectra.psd_lfp.sustained = nan(size(spectra.freq));
spectra.psd_fiber.sustained = nan(size(spectra.freq));

% Use mask_stim_full (full 1s stim period) for fair comparison
n_samples_stim_full = sum(mask_stim_full);
fprintf('      Stim period (full): %d samples available (%.2fs at %.1f Hz)\n', ...
    n_samples_stim_full, n_samples_stim_full/fs, fs);

% Use period-specific parameters (0.5s segments with 50% overlap)
period_seg = round(PERIOD_EPOCH_LENGTH_SEC * fs);
period_ovlp = round(PERIOD_EPOCH_OVERLAP_SEC * fs);

% Minimum data: 2 * period_seg = 1.0s with 0.5s segments
if n_samples_stim_full >= period_seg * 2
    try
        period_nfft = 2^nextpow2(period_seg * 2);
        period_win = hanning(period_seg);
        
        lfp_stim = lfp(mask_stim_full);
        gevi_stim = gevi(mask_stim_full);
        
        % Check if we have enough data
        if length(lfp_stim) >= period_seg && length(gevi_stim) >= period_seg
            % Get frequency range from config
            freq_min = cfg.coherence.fieldtrip.foi_min;
            freq_max = cfg.coherence.fieldtrip.foi_max;
            
            [Cxy, F] = mscohere(lfp_stim, gevi_stim, period_win, period_ovlp, period_nfft, fs);
            
            % Validate coherence values (should be between 0 and 1)
            Cxy = max(0, min(1, Cxy));  % Clamp to [0, 1]
            
            % Interpolate only if we have valid data
            freq_mask = (F >= freq_min) & (F <= freq_max);
            if ~isempty(Cxy) && ~isempty(F) && any(freq_mask) && any(isfinite(Cxy(freq_mask)))
                spectra.coherence.sustained = interp1(F(freq_mask), Cxy(freq_mask), spectra.freq, 'linear', 'extrap');
                % Clamp interpolated values to [0, 1]
                spectra.coherence.sustained = max(0, min(1, spectra.coherence.sustained));
            else
                spectra.coherence.sustained = nan(size(spectra.freq));
            end
            
            [Plfp, F_psd] = pwelch(lfp_stim, period_win, period_ovlp, period_nfft, fs);
            [Pgevi, ~] = pwelch(gevi_stim, period_win, period_ovlp, period_nfft, fs);
            
            freq_mask_psd = (F_psd >= freq_min) & (F_psd <= freq_max);
            if ~isempty(Plfp) && ~isempty(F_psd) && any(freq_mask_psd) && any(isfinite(Plfp(freq_mask_psd)))
                spectra.psd_lfp.sustained = interp1(F_psd(freq_mask_psd), ...
                    10*log10(Plfp(freq_mask_psd) + 1e-10), spectra.freq, 'linear', 'extrap');
            else
                spectra.psd_lfp.sustained = nan(size(spectra.freq));
            end
            
            if ~isempty(Pgevi) && ~isempty(F_psd) && any(freq_mask_psd) && any(isfinite(Pgevi(freq_mask_psd)))
                spectra.psd_fiber.sustained = interp1(F_psd(freq_mask_psd), ...
                    10*log10(Pgevi(freq_mask_psd) + 1e-10), spectra.freq, 'linear', 'extrap');
            else
                spectra.psd_fiber.sustained = nan(size(spectra.freq));
            end
            
            fprintf('      Stim (full): mscohere succeeded - coh size=%d, lfp_psd size=%d, fiber_psd size=%d\n', ...
                length(spectra.coherence.sustained), length(spectra.psd_lfp.sustained), ...
                length(spectra.psd_fiber.sustained));
        else
            % Not enough data
            spectra.coherence.sustained = nan(size(spectra.freq));
            spectra.psd_lfp.sustained = nan(size(spectra.freq));
            spectra.psd_fiber.sustained = nan(size(spectra.freq));
        end
    catch ME
        warning('Stim period computation failed: %s', ME.message);
        spectra.coherence.sustained = nan(size(spectra.freq));
        spectra.psd_lfp.sustained = nan(size(spectra.freq));
        spectra.psd_fiber.sustained = nan(size(spectra.freq));
    end
else
    % Not enough samples for stim period
    fprintf('      Stim (full): Not enough samples (%d < %d)\n', n_samples_stim_full, period_seg * 2);
    spectra.coherence.sustained = nan(size(spectra.freq));
    spectra.psd_lfp.sustained = nan(size(spectra.freq));
    spectra.psd_fiber.sustained = nan(size(spectra.freq));
end

% -------------------------------------------------------------------------
% POST-STIM (uses period-specific shorter pseudotrials for 1s window)
% -------------------------------------------------------------------------
if sum(mask_post) >= round(PERIOD_EPOCH_LENGTH_SEC * fs * 2)
    ft_data_post = create_pseudotrials(lfp(mask_post), gevi(mask_post), fs, PERIOD_EPOCH_LENGTH_SEC, PERIOD_EPOCH_OVERLAP_SEC);
    if ~isempty(ft_data_post)
        try
            freq_post = ft_freqanalysis(ft_cfg, ft_data_post);
            coh_post = ft_connectivityanalysis(coh_cfg, freq_post);
            
            % Extract coherence with interpolation if needed
            coh_vals = squeeze(coh_post.cohspctrm(1,2,:));
            if length(coh_vals) == length(coh_post.freq) && length(coh_vals) ~= length(spectra.freq)
                spectra.coherence.post_stim = interp1(coh_post.freq(:), coh_vals(:), spectra.freq(:), 'linear', 'extrap');
            else
                spectra.coherence.post_stim = coh_vals(:);
            end
            if length(spectra.coherence.post_stim) ~= length(spectra.freq)
                spectra.coherence.post_stim = nan(size(spectra.freq));
            end
            
            freq_pow_post = ft_freqanalysis(ft_cfg_pow, ft_data_post);
            spectra.psd_lfp.post_stim = 10*log10(squeeze(freq_pow_post.powspctrm(1,:))' + 1e-10);
            spectra.psd_fiber.post_stim = 10*log10(squeeze(freq_pow_post.powspctrm(2,:))' + 1e-10);
        catch
            spectra.coherence.post_stim = nan(size(spectra.freq));
            spectra.psd_lfp.post_stim = nan(size(spectra.freq));
            spectra.psd_fiber.post_stim = nan(size(spectra.freq));
        end
    else
        spectra.coherence.post_stim = nan(size(spectra.freq));
        spectra.psd_lfp.post_stim = nan(size(spectra.freq));
        spectra.psd_fiber.post_stim = nan(size(spectra.freq));
    end
else
    spectra.coherence.post_stim = nan(size(spectra.freq));
    spectra.psd_lfp.post_stim = nan(size(spectra.freq));
    spectra.psd_fiber.post_stim = nan(size(spectra.freq));
end

end

%% ============================================================================
%  HELPER: Time-Resolved Coherence
%  ============================================================================
function time_resolved = compute_time_resolved_coherence(data, method, cfg)
%COMPUTE_TIME_RESOLVED_COHERENCE Sliding window coherence for heatmaps

fs = data.fs;
lfp = data.lfp;
gevi = data.gevi;
t = data.time;

time_resolved = struct();

% Use shorter window for stimulation data
window_sec = cfg.coherence.mscohere.time_window_sec;
step_sec = cfg.coherence.mscohere.time_step_sec;
window_samples = round(window_sec * fs);
step_samples = round(step_sec * fs);

n_samples = length(lfp);
n_windows = floor((n_samples - window_samples) / step_samples) + 1;

if n_windows < 1
    time_resolved.coherence = [];
    time_resolved.time = [];
    time_resolved.freq = [];
    return;
end

% Get frequency vector
segment_samples = round(cfg.coherence.mscohere.segment_sec * fs);
overlap_samples = round(cfg.coherence.mscohere.overlap_frac * segment_samples);
nfft = cfg.coherence.mscohere.nfft_factor * segment_samples;
freq_min = cfg.coherence.mscohere.freq_min;
freq_max = cfg.coherence.mscohere.freq_max;
window_hann = hanning(segment_samples);

% Compute the first window up front to get the frequency axis, then reuse
% its result instead of recomputing window 1 again inside the loop below.
[Cxy_first, F] = mscohere(lfp(1:window_samples), gevi(1:window_samples), ...
    window_hann, overlap_samples, nfft, fs);
freq_mask = (F >= freq_min) & (F <= freq_max);
F_out = F(freq_mask);

% Preallocate
coh_matrix = zeros(length(F_out), n_windows);
time_centers = zeros(1, n_windows);

coh_matrix(:, 1) = Cxy_first(freq_mask);
time_centers(1) = t(round(window_samples / 2));

for w = 2:n_windows
    start_idx = (w - 1) * step_samples + 1;
    end_idx = start_idx + window_samples - 1;
    
    if end_idx > n_samples
        break;
    end
    
    lfp_win = lfp(start_idx:end_idx);
    gevi_win = gevi(start_idx:end_idx);
    
    try
        [Cxy, ~] = mscohere(lfp_win, gevi_win, window_hann, overlap_samples, nfft, fs);
        coh_matrix(:, w) = Cxy(freq_mask);
    catch
        coh_matrix(:, w) = nan(length(F_out), 1);
    end
    
    % Time center of window
    center_idx = start_idx + window_samples / 2;
    time_centers(w) = t(round(center_idx));
end

time_resolved.coherence = coh_matrix;
time_resolved.time = time_centers;
time_resolved.freq = F_out;

end

%% ============================================================================
%  HELPER: Compute Spectrograms
%  ============================================================================
function [spec_lfp, spec_fiber] = compute_spectrograms(data, cfg)
%COMPUTE_SPECTROGRAMS Compute time-frequency spectrograms

fs = data.fs;
t = data.time;

spec_window_samples = round(cfg.spectrogram.window_sec * fs);
spec_overlap_samples = round(cfg.spectrogram.overlap_frac * spec_window_samples);

% Use nfft_mult from config if available, otherwise use default multiplier
if isfield(cfg.spectrogram, 'nfft_mult')
    spec_nfft = 2^nextpow2(spec_window_samples * cfg.spectrogram.nfft_mult);
else
    spec_nfft = 2^nextpow2(spec_window_samples * 2);  % Default: 2x window
end

freq_min = cfg.spectrogram.freq_min;
freq_max = cfg.spectrogram.freq_max;

% LFP spectrogram
[S_lfp, F_spec, T_spec] = spectrogram(data.lfp, spec_window_samples, ...
    spec_overlap_samples, spec_nfft, fs);
spec_power_lfp = abs(S_lfp).^2 / (fs * spec_window_samples);

% Fiber spectrogram
[S_fiber, ~, ~] = spectrogram(data.gevi, spec_window_samples, ...
    spec_overlap_samples, spec_nfft, fs);
spec_power_fiber = abs(S_fiber).^2 / (fs * spec_window_samples);

% Limit frequency range
freq_idx = (F_spec >= freq_min) & (F_spec <= freq_max);
F_out = F_spec(freq_idx);

% Apply smoothing
spec_power_lfp = smooth2a_radius_kernel(spec_power_lfp(freq_idx, :), ...
    cfg.spectrogram.smooth_freq, cfg.spectrogram.smooth_time);
spec_power_fiber = smooth2a_radius_kernel(spec_power_fiber(freq_idx, :), ...
    cfg.spectrogram.smooth_freq, cfg.spectrogram.smooth_time);

% Shift time
T_out = T_spec + t(1);

% Store
spec_lfp = struct();
spec_lfp.power = 10*log10(spec_power_lfp + 1e-10);
spec_lfp.freq = F_out;
spec_lfp.time = T_out;

spec_fiber = struct();
spec_fiber.power = 10*log10(spec_power_fiber + 1e-10);
spec_fiber.freq = F_out;
spec_fiber.time = T_out;

end

%% ============================================================================
%  HELPER: 2D Smoothing (from baseline pipeline)
%  ============================================================================
function out = smooth2a_radius_kernel(data, n_freq, n_time)
%SMOOTH2A_RADIUS_KERNEL 2D smoothing via a single (2*n+1)-per-side box kernel.
%   NOT the same algorithm as core/utils/smooth2a.m (separable per-axis
%   window-length convolution). Here n_freq/n_time are kernel RADII (kernel
%   size = 2*n+1), and smoothing is skipped only when BOTH radii are <= 0.
%   Deliberately named differently and kept local to this file so it is
%   never confused with, or silently merged into, the shared utility.

if n_freq <= 0 && n_time <= 0
    out = data;
    return;
end

% Create 2D averaging kernel
kernel_freq = max(1, 2*n_freq + 1);
kernel_time = max(1, 2*n_time + 1);
kernel = ones(kernel_freq, kernel_time) / (kernel_freq * kernel_time);

% Apply convolution
out = conv2(data, kernel, 'same');

end
