function cfg = stim_analysis_config()
%% ============================================================================
%  STIMULATION ANALYSIS CONFIGURATION
%  ============================================================================
%  Configuration file for stimulation trial spectral analysis.
%  This is SEPARATE from baseline analysis to avoid any interference.
%
%  USAGE:
%    cfg = stim_analysis_config();
%  ============================================================================

%% ============================================================================
%  STIMULATION DATABASE  -- EDIT THIS to your own cohort
%  ============================================================================
%  Define all animals and sessions with stimulation trials here.
%  Each session has specific stimulation parameters that will be stored
%  in the output for later comparison in Python.

cfg.stim_database = struct();

% -------------------------------------------------------------------------
% ANIMAL 1 -- example, 1s stimulation (bilateral)
% -------------------------------------------------------------------------
cfg.stim_database(1).mouse_id = 'Animal01';
cfg.stim_database(1).project = 'FiberVoltageImaging';
cfg.stim_database(1).data_root = fullfile(lab_paths().data_root, 'FiberVoltageImaging', 'Animal01', 'Fiber_Voltage_Processed');

% Session 1: 40Hz stimulation
cfg.stim_database(1).sessions(1).session_id = '01_09_25-R1';
cfg.stim_database(1).sessions(1).condition_name = '40Hz';
cfg.stim_database(1).sessions(1).num_trials = 10;
cfg.stim_database(1).sessions(1).stim_params.frequency_hz = 40;
cfg.stim_database(1).sessions(1).stim_params.pulse_width_us = 100;
cfg.stim_database(1).sessions(1).stim_params.voltage = 6;
cfg.stim_database(1).sessions(1).stim_params.location = 'bilateralANT';
cfg.stim_database(1).sessions(1).trial_folder_pattern = 'Trial%d_fov1_100us_40Hz_6V_1sec_bilateralANT_10_trials_1';
cfg.stim_database(1).sessions(1).mat_file_pattern = 'Animal01-%s_Trial%d_FiberPhotometry_Analysis.mat';

% Session 2: 135Hz stimulation
cfg.stim_database(1).sessions(2).session_id = '01_09_25-R2';
cfg.stim_database(1).sessions(2).condition_name = '135Hz';
cfg.stim_database(1).sessions(2).num_trials = 10;
cfg.stim_database(1).sessions(2).stim_params.frequency_hz = 135;
cfg.stim_database(1).sessions(2).stim_params.pulse_width_us = 100;
cfg.stim_database(1).sessions(2).stim_params.voltage = 6;
cfg.stim_database(1).sessions(2).stim_params.location = 'bilateralANT';
cfg.stim_database(1).sessions(2).trial_folder_pattern = 'Trial%d_fov1_100us_135Hz_6V_1sec_bilateralANT_10_trials_1';
cfg.stim_database(1).sessions(2).mat_file_pattern = 'Animal01-%s_Trial%d_FiberPhotometry_Analysis.mat';

% -------------------------------------------------------------------------
% ANIMAL 2 -- example, 1s DBS with amplitude- and energy-balanced conditions
% -------------------------------------------------------------------------
cfg.stim_database(2).mouse_id = 'Animal02';
cfg.stim_database(2).project = 'FiberVoltageImaging';
cfg.stim_database(2).data_root = fullfile(lab_paths().data_root, 'FiberVoltageImaging', 'Animal02', 'Fiber_Voltage_Processed');

% Session 1: 135Hz DBS (1 second)
cfg.stim_database(2).sessions(1).session_id = '01_02_26-R6';
cfg.stim_database(2).sessions(1).condition_name = '135Hz';
cfg.stim_database(2).sessions(1).num_trials = 10;
cfg.stim_database(2).sessions(1).stim_params.frequency_hz = 135;
cfg.stim_database(2).sessions(1).stim_params.pulse_width_us = 100;
cfg.stim_database(2).sessions(1).stim_params.voltage = 2.5;
cfg.stim_database(2).sessions(1).stim_params.location = 'ipsiHPDBS';
% Trial folder: Trial1_fov1_100us_135Hz_2.5V_1sec_ipsiHPDBS_HPipsicontra_1
cfg.stim_database(2).sessions(1).trial_folder_pattern = 'Trial%d_fov1_100us_135Hz_2.5V_1sec_ipsiHPDBS_HPipsicontra_%d';
cfg.stim_database(2).sessions(1).mat_file_pattern = 'Animal02-%s_Trial%d_FiberPhotometry_Analysis.mat';

% Session 2: 40Hz DBS - Amplitude Balanced (1 second) - 2.5V same as 135Hz
cfg.stim_database(2).sessions(2).session_id = '01_02_26-R9';
cfg.stim_database(2).sessions(2).condition_name = '40Hz_AmpBalanced';
cfg.stim_database(2).sessions(2).num_trials = 10;
cfg.stim_database(2).sessions(2).stim_params.frequency_hz = 40;
cfg.stim_database(2).sessions(2).stim_params.pulse_width_us = 100;
cfg.stim_database(2).sessions(2).stim_params.voltage = 2.5;
cfg.stim_database(2).sessions(2).stim_params.location = 'ipsiHPDBS';
cfg.stim_database(2).sessions(2).stim_params.balancing = 'amplitude';
% Trial folder: Trial1_fov1_100us_40Hz_2.5V_1sec_ipsiHPDBS_HPipsicontra_1
cfg.stim_database(2).sessions(2).trial_folder_pattern = 'Trial%d_fov1_100us_40Hz_2.5V_1sec_ipsiHPDBS_HPipsicontra_%d';
cfg.stim_database(2).sessions(2).mat_file_pattern = 'Animal02-%s_Trial%d_FiberPhotometry_Analysis.mat';

% Session 3: 40Hz DBS - Energy Balanced (1 second) - 4.5V for energy matching
cfg.stim_database(2).sessions(3).session_id = '01_02_26-R10';
cfg.stim_database(2).sessions(3).condition_name = '40Hz_EnergyBalanced';
cfg.stim_database(2).sessions(3).num_trials = 10;
cfg.stim_database(2).sessions(3).stim_params.frequency_hz = 40;
cfg.stim_database(2).sessions(3).stim_params.pulse_width_us = 100;
cfg.stim_database(2).sessions(3).stim_params.voltage = 4.5;
cfg.stim_database(2).sessions(3).stim_params.location = 'ipsiHPDBS';
cfg.stim_database(2).sessions(3).stim_params.balancing = 'energy';
% Trial folder: Trial1_fov1_100us_40Hz_4.5V_1sec_ipsiHPDBS_HPipsicontra_1
cfg.stim_database(2).sessions(3).trial_folder_pattern = 'Trial%d_fov1_100us_40Hz_4.5V_1sec_ipsiHPDBS_HPipsicontra_%d';
cfg.stim_database(2).sessions(3).mat_file_pattern = 'Animal02-%s_Trial%d_FiberPhotometry_Analysis.mat';

% -------------------------------------------------------------------------
% ANIMAL 3 -- example, 10s DBS (longer-duration protocol)
% -------------------------------------------------------------------------
cfg.stim_database(3).mouse_id = 'Animal03';
cfg.stim_database(3).project = 'FiberVoltageImaging';
cfg.stim_database(3).data_root = fullfile(lab_paths().data_root, 'FiberVoltageImaging', 'Animal03', 'Fiber_Voltage_Processed');

% Session 1: 135Hz DBS (10 seconds) - Reference for comparisons
cfg.stim_database(3).sessions(1).session_id = '01_03_26-R3';
cfg.stim_database(3).sessions(1).condition_name = '135Hz';
cfg.stim_database(3).sessions(1).num_trials = 10;
cfg.stim_database(3).sessions(1).stim_params.frequency_hz = 135;
cfg.stim_database(3).sessions(1).stim_params.pulse_width_us = 100;
cfg.stim_database(3).sessions(1).stim_params.voltage = 2.9;
cfg.stim_database(3).sessions(1).stim_params.duration_sec = 10;
cfg.stim_database(3).sessions(1).stim_params.location = 'ipsiANTDBS';
% Trial folder: Trial1_fov1_100us_135Hz_2.9V_10sec_ipsiANTDBS_ANT_1
cfg.stim_database(3).sessions(1).trial_folder_pattern = 'Trial%d_fov1_100us_135Hz_2.9V_10sec_ipsiANTDBS_ANT_%d';
cfg.stim_database(3).sessions(1).mat_file_pattern = 'Animal03-%s_Trial%d_FiberPhotometry_Analysis.mat';

% Session 2: 40Hz DBS - Amplitude Balanced (10 seconds) - 2.9V same as 135Hz
cfg.stim_database(3).sessions(2).session_id = '01_03_26-R5';
cfg.stim_database(3).sessions(2).condition_name = '40Hz_AmpBalanced';
cfg.stim_database(3).sessions(2).num_trials = 10;
cfg.stim_database(3).sessions(2).stim_params.frequency_hz = 40;
cfg.stim_database(3).sessions(2).stim_params.pulse_width_us = 100;
cfg.stim_database(3).sessions(2).stim_params.voltage = 2.9;
cfg.stim_database(3).sessions(2).stim_params.duration_sec = 10;
cfg.stim_database(3).sessions(2).stim_params.location = 'ipsiANTDBS';
cfg.stim_database(3).sessions(2).stim_params.balancing = 'amplitude';
% Trial folder: Trial1_fov1_100us_40Hz_2.9V_10sec_ipsiANTDBS_ANT_1
cfg.stim_database(3).sessions(2).trial_folder_pattern = 'Trial%d_fov1_100us_40Hz_2.9V_10sec_ipsiANTDBS_ANT_%d';
cfg.stim_database(3).sessions(2).mat_file_pattern = 'Animal03-%s_Trial%d_FiberPhotometry_Analysis.mat';

% Session 3: 40Hz DBS - Energy Balanced (10 seconds) - 4.5V for energy matching
cfg.stim_database(3).sessions(3).session_id = '01_03_26-R6';
cfg.stim_database(3).sessions(3).condition_name = '40Hz_EnergyBalanced';
cfg.stim_database(3).sessions(3).num_trials = 10;
cfg.stim_database(3).sessions(3).stim_params.frequency_hz = 40;
cfg.stim_database(3).sessions(3).stim_params.pulse_width_us = 100;
cfg.stim_database(3).sessions(3).stim_params.voltage = 4.5;
cfg.stim_database(3).sessions(3).stim_params.duration_sec = 10;
cfg.stim_database(3).sessions(3).stim_params.location = 'ipsiANTDBS';
cfg.stim_database(3).sessions(3).stim_params.balancing = 'energy';
% Trial folder: Trial1_fov1_100us_40Hz_4.5V_10sec_ipsiANTDBS_ANT_1
cfg.stim_database(3).sessions(3).trial_folder_pattern = 'Trial%d_fov1_100us_40Hz_4.5V_10sec_ipsiANTDBS_ANT_%d';
cfg.stim_database(3).sessions(3).mat_file_pattern = 'Animal03-%s_Trial%d_FiberPhotometry_Analysis.mat';

%% ============================================================================
%  STIMULATION TIMING PARAMETERS
%  ============================================================================
%  These define the periods for analysis.
%  Times are relative to recording start (stim onset is detected from data).
%
%  DEFAULT TIMING (1s-stimulation animals, e.g. Animal01/Animal02 above):
cfg.stim_timing.pre_stim_duration_sec = 4.0;      % Duration of pre-stim period (full 4s in recording)
cfg.stim_timing.stim_duration_sec = 1.0;           % Total stimulation duration
cfg.stim_timing.post_stim_duration_sec = 5.0;      % Duration of post-stim period (full 5s in recording)

%  10s-STIMULATION TIMING (e.g. Animal03 above):
%  To process a 10s-stim animal, uncomment these and comment out the defaults above:
% cfg.stim_timing.pre_stim_duration_sec = 10.0;     % Duration of pre-stim period (10s)
% cfg.stim_timing.stim_duration_sec = 10.0;         % Total stimulation duration (10s)
% cfg.stim_timing.post_stim_duration_sec = 10.0;    % Duration of post-stim period (10s)

% Stimulation sub-periods
% DEFAULT (1s stim): Short transient
cfg.stim_timing.transient_end_sec = 0.15;          % Transient period: 0 to 0.15s
cfg.stim_timing.sustained_start_sec = 0.15;        % Sustained period: 0.15s to stim_duration

% 10s stim: Longer transient to capture adaptation
% cfg.stim_timing.transient_end_sec = 1.0;          % Transient period: 0 to 1s
% cfg.stim_timing.sustained_start_sec = 1.0;        % Sustained period: 1s to stim_duration (9s)

% Period coherence analysis windows
% DEFAULT (1s stim): Use 1s windows for fair comparison
cfg.stim_timing.coherence_prestim_window_sec = 1.0;   % Use last 1s before stim
cfg.stim_timing.coherence_poststim_window_sec = 1.0;   % Use first 1s after stim

% 10s stim: Use larger windows for better frequency resolution
% cfg.stim_timing.coherence_prestim_window_sec = 8.0;   % Use last 8s before stim (skip first 2s)
% cfg.stim_timing.coherence_poststim_window_sec = 8.0;  % Use first 8s after stim

% Period-specific coherence parameters
% DEFAULT (1s stim): 0.5s segments for compatibility with short windows
cfg.coherence.period_segment_sec = 0.5;              % Shorter segment for period coherence
cfg.coherence.period_overlap_frac = 0.5;             % 50% overlap

% 10s stim: Longer segments for better frequency resolution (~1Hz bins)
% cfg.coherence.period_segment_sec = 1.0;              % 1s segment for better freq resolution
% cfg.coherence.period_overlap_frac = 0.5;             % 50% overlap for period coherence

% FieldTrip period-specific parameters
% DEFAULT (1s stim):
cfg.coherence.fieldtrip.period_pseudotrial_length_sec = 0.5;
cfg.coherence.fieldtrip.period_pseudotrial_overlap_sec = 0.25;

% 10s stim: Longer pseudotrials
% cfg.coherence.fieldtrip.period_pseudotrial_length_sec = 1.0;   % 1s pseudotrials
% cfg.coherence.fieldtrip.period_pseudotrial_overlap_sec = 0.5;  % 50% overlap

%% ============================================================================
%  SPECTRAL ANALYSIS PARAMETERS
%  ============================================================================
%  These match the baseline analysis config for consistency.

% Which methods to use
cfg.methods = {'mscohere', 'fieldtrip'};  % Options: 'mscohere', 'fieldtrip', or both

% Fiber/channel selection
cfg.fiber_index = 1;

%% ============================================================================
%  COHERENCE PARAMETERS - MSCOHERE
%  ============================================================================

cfg.coherence.mscohere.segment_sec = 1.0;     % Segment length (seconds)
cfg.coherence.mscohere.overlap_frac = 0.8;    % Overlap fraction
cfg.coherence.mscohere.freq_min = 2;          % Minimum frequency (Hz)
cfg.coherence.mscohere.freq_max = 150;        % Maximum frequency (Hz) - extended for 135Hz stim
cfg.coherence.mscohere.nfft_factor = 2;       % NFFT = factor × segment_samples

% Time-resolved coherence parameters - OPTIMIZED FOR SHORT TRIALS
cfg.coherence.mscohere.time_window_sec = 0.5;     % Window for time-resolved (s) - SHORTER for 10s trials
cfg.coherence.mscohere.time_step_sec = 0.05;      % Step for time-resolved (s) - FINER 50ms steps
cfg.coherence.mscohere.smooth_freq = 0;           % Frequency smoothing (bins) - reduce for sharpness
cfg.coherence.mscohere.smooth_time = 0;           % Time smoothing (bins) - reduce for sharpness

%% ============================================================================
%  COHERENCE PARAMETERS - FIELDTRIP
%  ============================================================================

cfg.coherence.fieldtrip.method = 'mtmfft';    % Method
cfg.coherence.fieldtrip.taper = 'hanning';   % Single Hanning taper (better for narrow-band stimulation frequencies)
cfg.coherence.fieldtrip.tapsmofrq = 2;        % Frequency smoothing (Hz) - only used with 'dpss'
cfg.coherence.fieldtrip.foi_min = 2;          % Minimum frequency (Hz)
cfg.coherence.fieldtrip.foi_max = 150;         % Maximum frequency (Hz) - extended for 135Hz stim
cfg.coherence.fieldtrip.foi_step = 0.5;       % Frequency resolution (Hz)

% Pseudo-trial parameters
cfg.coherence.fieldtrip.pseudotrial_length_sec = 1.0;    % seconds
cfg.coherence.fieldtrip.pseudotrial_overlap_sec = 0.8;   % seconds (80%)

% Time-resolved coherence parameters - OPTIMIZED FOR SHORT TRIALS
cfg.coherence.fieldtrip.time_window_epochs = 3;    % Epochs per time point (FEWER for short data)
cfg.coherence.fieldtrip.smooth_freq = 0;           % Frequency smoothing (bins) - reduce for sharpness
cfg.coherence.fieldtrip.smooth_time = 0;           % Time smoothing (bins) - reduce for sharpness

%% ============================================================================
%  SPECTROGRAM PARAMETERS - OPTIMIZED FOR SHORT (10s) TRIALS
%  ============================================================================
%  For 10s trials, we need balanced time-frequency resolution. Trade-offs:
%  - Shorter window = better time resolution, worse frequency resolution
%  - Longer window = better frequency resolution, worse time resolution
%  - Higher overlap = smoother appearance, but can blur sharp transitions

cfg.spectrogram.window_sec = 0.6;     % Window length (seconds) - OPTIMIZED: closer to preprocessing (0.96s) but adapted for short trials
cfg.spectrogram.overlap_frac = 0.88;   % Overlap fraction (88%) - OPTIMIZED: high overlap for smoothness (preprocessing uses 93.5%)
cfg.spectrogram.nfft_mult = 3;         % NFFT multiplier (NFFT = nfft_mult * window_samples) for freq resolution
cfg.spectrogram.freq_min = 2;         % Minimum frequency (Hz)
cfg.spectrogram.freq_max = 150;       % Maximum frequency (Hz) - extended for 135Hz stim
cfg.spectrogram.smooth_freq = 1;       % Frequency smoothing (bins) - light smoothing
cfg.spectrogram.smooth_time = 1;       % Time smoothing (bins) - light smoothing

%% ============================================================================
%  PSD PARAMETERS
%  ============================================================================

cfg.psd.window_sec = 1.0;             % Window length (seconds)
cfg.psd.overlap_frac = 0.8;           % Overlap fraction
cfg.psd.freq_min = 2;                 % Minimum frequency (Hz)
cfg.psd.freq_max = 150;               % Maximum frequency (Hz) - extended for 135Hz stim

%% ============================================================================
%  MOTION CONVERSION CONSTANTS
%  ============================================================================

cfg.motion.wheel_diameter_cm = 19.0;
cfg.motion.encoder_counts_per_rev = 1024;
cfg.motion.ephys_sampling_rate = 30000;
cfg.motion.smooth_samples = 10;

%% ============================================================================
%  PATH CONFIGURATION
%  ============================================================================

% Sourced from the centralised config (config/lab_paths.m). Override per machine
% via config/paths_local.m. See config/README.md.
cfg_dir   = fileparts(mfilename('fullpath'));            % .../spectral_analysis/config
repo_root = fileparts(fileparts(cfg_dir));               % project root
addpath(fullfile(repo_root, 'config'));
lp = lab_paths();

cfg.paths.output_root = fullfile(lp.figures_root, 'Stimulation_analysis', 'Spectral_data_outputs');
cfg.paths.code_root   = lp.spectral_code;

% Toolbox paths
cfg.paths.toolboxes   = lp.toolboxes;

%% ============================================================================
%  FIGURE PARAMETERS
%  ============================================================================

cfg.figure.dpi = 300;
cfg.figure.format = 'png';
cfg.figure.font_size_title = 14;
cfg.figure.font_size_label = 12;
cfg.figure.font_size_tick = 10;

% Colors for periods
cfg.figure.colors.pre_stim = [0.5, 0.5, 0.5];     % Grey
cfg.figure.colors.transient = [0.8, 0.2, 0.2];    % Red
cfg.figure.colors.sustained = [0.2, 0.6, 0.8];    % Blue
cfg.figure.colors.post_stim = [0.4, 0.7, 0.4];    % Green

%% ============================================================================
%  OUTPUT CONFIGURATION
%  ============================================================================

cfg.output.save_mat = true;
cfg.output.compress_mat = true;
cfg.output.save_raw_traces = true;    % Include raw LFP/fiber/motion in output

end
