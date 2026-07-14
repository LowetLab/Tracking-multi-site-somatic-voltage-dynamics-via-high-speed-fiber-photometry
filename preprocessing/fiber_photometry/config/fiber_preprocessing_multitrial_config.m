%% ============================================================================
%  FIBER PREPROCESSING -- MULTI-TRIAL RUN CONFIGURATION
%  ============================================================================
%  This is the ONE place to edit run parameters for the multi-trial pipeline.
%  It is a script-include: run_fiber_preprocessing_multitrial.m loads it
%  with run(), so every variable below is created directly in that workspace
%  with the same name it has always had (nothing else in the script changed).
%
%  Edit MOUSE_NAME / RECORDING_DATE / RECORDING_ID (and paths/options) for your
%  recording, then run the main script. base_folder and the STIM_* init values
%  are derived here from the parameters above them -- leave those as-is unless
%  you know why you are changing them.
%  ============================================================================

%% ANALYSIS MODE
% This script is ONLY for multi-trial processing
ANALYSIS_MODE = 'multi_trial';  % Fixed to multi_trial mode

%% ============================================================================
%  USER CONFIGURATION - MODIFY THESE PARAMETERS ONLY
%  ============================================================================

%% MOUSE AND RECORDING CONFIGURATION  -- EDIT THESE FOR YOUR RECORDING
MOUSE_NAME = 'Animal01';                    % Your animal/mouse identifier (matches its data folder name)
RECORDING_DATE = '01_01_25';                % Format: DD_MM_YY
RECORDING_ID = 'R1';                        % Format: R## (e.g., R14, R4)
EXPERIMENTER = 'YourInitials';              % Experimenter initials/name
EXPERIMENT_TYPE = 'DBS';                    % Session-type folder name used in your data structure
DATA_TYPE_IMAGING = 'ImagingData';
DATA_TYPE_OPEN_EPHYS = 'OpenEphys';

%% BASE PATH CONFIGURATION
% Base path up to the mouse folder -- EDIT THIS to your own data root, or
% (preferred) leave it derived from config/lab_paths.m below.
% BASE_PATH_ROOT = 'D:\Imaging_Data\';
BASE_PATH_ROOT = fullfile(lab_paths().data_root, 'FiberVoltageImaging');

% Construct base folder path automatically
base_folder = fullfile(BASE_PATH_ROOT, MOUSE_NAME, DATA_TYPE_IMAGING, EXPERIMENT_TYPE, RECORDING_DATE, RECORDING_ID);

%% PROCESSING OPTIONS
MOTION_CORRECTION = 0;           % 0 = off, 1 = on
CORRECTION_TYPE = 'rigid';       % 'rigid' or 'non-rigid'
PROCESS_FULL_FIELD = false;       % true = whole field average, false = manual ROI selection
INVERT_TRACE = false;            % Set true to invert fiber traces
ROI_SELECTION = 1;               % For multi-FOV: 0 = automatic, 1 = manual ROI selection

%% ENHANCED MULTI-FIBER VISUALIZATION (NEW)
GENERATE_MULTI_FIBER_PLOT = false;   % Generate combined multi-fiber visualization
MAX_FIBERS_DISPLAY = 6;            % Maximum fibers to display in combined plot

%% PHOTOBLEACHING CORRECTION
APPLY_PHOTOBLEACHING_CORRECTION = true;  % Enable photobleaching correction

%% SAMPLING RATES
EPHYS_FS = 30000;                % Open Ephys sampling rate (Hz) - fixed
IMAGING_FS = [];                 % Will be calculated from camera triggers

%% TIME PERIOD CONFIGURATION (FALLBACK ONLY - not used if automatic detection succeeds)
% These values are ONLY used as fallback if automatic stimulation detection fails
% For successful automatic detection, pre-stim period = [0, stim_onset] and
% stim period = [stim_onset, stim_onset + STIMULATION_DURATION_SEC]
PRE_STIM_PERIOD = [0, 10];        % [start, end] in seconds - FALLBACK baseline period
STIM_PERIOD = [10,20];            % [start, end] in seconds - FALLBACK stimulation period
% Post-stim period automatically extends from stim end to trial end

%% OPEN EPHYS CONFIGURATION
LOAD_EPHYS_DATA = true;          % Set false to skip Open Ephys data loading
EPHYS_FILE_PREFIX = 100;         % File prefix for Open Ephys files (usually 100)
LOAD_mPFC_LFP = true;        % Set true if mPFC LFP is also recorded as Ch1-Ch3 signal
LOAD_ipsiHP_LFP = true;      % Set true if ipsislateral HP LFP is also recorded as Ch2-Ch4 signal

%% FREQUENCY BAND DEFINITIONS
BAND_NAMES = {'Delta-Theta', 'Alpha', 'Beta', 'Low Gamma', 'High Gamma'};
BAND_RANGES = [1 8; 8 12; 13 30; 31 70; 71 185];  % [low high] in Hz
BAND_COLORS = [0 0 1; 0 1 1; 0 0.6 0; 1 0.5 0; 1 0 0];  % RGB colors

%% PLOT CONFIGURATION - Enable/disable specific plots
PLOTS = struct();
PLOTS.fiber_trace_comparison = true;      % Raw vs corrected fiber traces
PLOTS.fiber_spectrogram = true;           % Time-frequency spectrogram
PLOTS.fiber_band_power = true;            % Spectral band power over time
PLOTS.fiber_psd = true;                   % Power spectral density comparison
PLOTS.photobleaching_methods = true;      % Comparison of correction methods
PLOTS.exponential_fit_quality = true;     % Quality of exponential fit
PLOTS.lfp_fiber_spectrogram = true;       % Combined LFP and Fiber spectrograms
PLOTS.phase_locking_overall = true;       % Overall phase-locking vs frequency
PLOTS.phase_locking_behavior = true;      % Phase-locking during run vs rest
PLOTS.envelope_correlation = true;        % Theta envelope correlation

%% STIMULATION DETECTION
AUTO_DETECT_STIMULATION = true;  % Set false to use hardcoded values only
% BASELINE TRIAL SUPPORT:
% IS_BASELINE_TRIAL = []  -> Auto-detect (baseline if <5 pulses detected)
% IS_BASELINE_TRIAL = true -> Force baseline trial (skip stim detection)
% IS_BASELINE_TRIAL = false -> Force stimulation trial (use hardcoded periods if detection fails)
IS_BASELINE_TRIAL = [];           % [] = auto-detect, true = force baseline, false = force stim trial

% STIMULATION DURATION (for automatic offset calculation)
% This is used to calculate stimulation offset from detected onset
% Only used when automatic detection succeeds
STIMULATION_DURATION_SEC = 10.0;  % Duration of stimulation in seconds (used for offset calculation)

% Initialize stimulation detection variables
STIM_ONSET_EPHYS_SAMPLE = [];
STIM_OFFSET_EPHYS_SAMPLE = [];
STIM_PERIOD_HARDCODED = STIM_PERIOD;  % Store original
