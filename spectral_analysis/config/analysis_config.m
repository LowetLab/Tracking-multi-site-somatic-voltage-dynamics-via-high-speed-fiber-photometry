function cfg = analysis_config()
%% ============================================================================
%  MASTER CONFIGURATION FILE FOR SPECTRAL ANALYSIS PIPELINE
%  ============================================================================
%  This file centralizes ALL user-configurable parameters for the spectral
%  analysis pipeline (coherence, PSD, band power, behavioural classification).
%
%  USAGE:
%    cfg = analysis_config();
%    % Then pass cfg to pipeline functions
%  ============================================================================

%% ============================================================================
%  PIPELINE CONTROL - WHAT TO RUN
%  ============================================================================

% Which analysis levels to run
cfg.run_single_trial = false;         % Individual trials (generates heatmaps + spectra)
cfg.run_session_pooled = false;       % Pool trials within each session
cfg.run_animal_pooled = false;        % Pool all sessions per animal (averages session spectra)
cfg.run_animal_concatenated = false;  % Concatenate all raw data across sessions (computes spectra once)
cfg.run_group_statistics = true;      % Cross-animal statistics

% Which coherence methods to use
cfg.methods = {'mscohere', 'fieldtrip'};  % Options: 'mscohere', 'fieldtrip', or both

% Which animals to process ({} = all animals in config/animal_session_database.m)
cfg.animals_to_process = {};

% Figure generation
cfg.generate_figures = true;          % Set false to only compute data structs
cfg.save_figures = true;              % Save figures to disk
cfg.close_figures_after_save = true;  % Close figures after saving (memory)

%% ============================================================================
%  REST/RUN CLASSIFICATION MODE
%  ============================================================================
%
%  MODE 1: 'standard' (Original method)
%    - RUN:  speed > run_threshold (default 3 cm/s)
%    - REST: everything else (speed <= run_threshold)
%    - Short bouts are merged into surrounding state
%    - Ensures REST + RUN = 100% of data
%
%  MODE 2: 'clear' (Strict classification - RECOMMENDED)
%    - RUN:  speed > run_threshold AND bout duration >= min_run_bout_sec
%    - REST: speed < rest_threshold AND bout duration >= min_rest_bout_sec
%    - EXCLUDED: intermediate speeds OR short bouts
%    - More conservative but cleaner separation
%    - REST + RUN + EXCLUDED = 100% (intermediate periods excluded)
%
%  ============================================================================

cfg.behavior.classification_mode = 'clear';     % 'standard' or 'clear'

% Speed thresholds (cm/s)
cfg.behavior.run_threshold = 2.0;     % Speed above this = RUN (default 2.0 cm/s)
cfg.behavior.rest_threshold = 0.1;    % Speed below this = REST (only used in 'clear' mode)

% Bout duration filtering - SEPARATE for REST and RUN
% In 'standard' mode: min_bout_duration_sec applies to both
% In 'clear' mode: use min_run_bout_sec and min_rest_bout_sec separately
cfg.behavior.min_bout_duration_sec = 0.5;    % Legacy: shared minimum bout (seconds)
cfg.behavior.min_run_bout_sec = 0.5;         % Minimum RUN bout duration (seconds)
cfg.behavior.min_rest_bout_sec = 0.5;        % Minimum REST bout duration (seconds)
cfg.behavior.apply_bout_filter = true;       % Apply bout duration filter

% Motion smoothing
cfg.behavior.motion_smooth_samples = 10;     % Samples for speed smoothing [LEGACY: 10]

%% ============================================================================
%  MOTION CONVERSION CONSTANTS (from legacy scripts)
%  ============================================================================
%  These values MUST match the hardware setup used during recording.
%  CRITICAL: WHEEL_DIAMETER_CM was 19.0 in legacy, NOT 20.0!

cfg.motion.wheel_diameter_cm = 19.0;             % [LEGACY: 19.0] Wheel diameter in cm
cfg.motion.encoder_counts_per_rev = 1024;        % [LEGACY: 1024] Yumo E6B2 encoder
cfg.motion.ephys_sampling_rate = 30000;          % [LEGACY: 30000] Open Ephys sampling rate

%% ============================================================================
%  ARTIFACT HANDLING
%  ============================================================================
%  How to handle artifacts detected by artifact_removal_lfp_multisession.m.
%  Dependent fields (output_folder_suffix, artifact_exclusion) are DERIVED from
%  cfg.artifact.mode in run_spectral_pipeline.m.
%
%    'none'    - Use all data as-is.
%    'exclude' - TRIAL EXCLUSION: skip whole trials with > threshold contamination.
%    'clean'   - DATA CLEANING: excise artifact segments from within each trial.
cfg.artifact.mode = 'clean';            % 'none' | 'exclude' | 'clean'
cfg.artifact.threshold = 0.10;          % 'exclude' mode: skip trials with > this fraction contaminated

% Cleaning parameters (only used in 'clean' mode)
cfg.artifact.pre_pad_sec = 0.150;       % Remove this many seconds before each artifact
cfg.artifact.post_pad_sec = 0.150;      % Remove this many seconds after each artifact
cfg.artifact.smooth_window_sec = 0.100; % Merge artifacts closer than this (0 = disabled)

%% ============================================================================
%  SPECTROGRAM PARAMETERS [LEGACY: figure2_coherence_mscohere.m lines 130-137]
%  ============================================================================

cfg.spectrogram.window_sec = 0.75;    % [LEGACY: 0.75] Window length (seconds)
cfg.spectrogram.overlap_frac = 0.9;   % [LEGACY: 0.9] Overlap fraction (0-1)
cfg.spectrogram.freq_min = 2;         % [LEGACY: 2] Minimum frequency (Hz)
cfg.spectrogram.freq_max = 70;        % [LEGACY: 70] Maximum frequency (Hz)
cfg.spectrogram.smooth_freq = 1;      % [LEGACY: 1] Frequency smoothing (bins)
cfg.spectrogram.smooth_time = 5;      % [LEGACY: 5] Time smoothing (bins)

%% ============================================================================
%  COHERENCE PARAMETERS - MSCOHERE [LEGACY: figure2_coherence_mscohere.m lines 147-154]
%  ============================================================================

cfg.coherence.mscohere.segment_sec = 1.0;     % [LEGACY: 1.0] Segment length (seconds)
cfg.coherence.mscohere.overlap_frac = 0.8;    % [LEGACY: 0.8] Overlap fraction
cfg.coherence.mscohere.freq_min = 2;          % [LEGACY: 2] Minimum frequency (Hz)
cfg.coherence.mscohere.freq_max = 70;         % [LEGACY: 70] Maximum frequency (Hz)
cfg.coherence.mscohere.nfft_factor = 2;       % NFFT = factor × segment_samples

% Time-resolved coherence parameters
cfg.coherence.mscohere.time_window_sec = 5.0;    % [LEGACY: 5.0] Window for time-resolved (s)
cfg.coherence.mscohere.time_step_sec = 0.25;     % [LEGACY: 0.25] Step for time-resolved (s)
cfg.coherence.mscohere.smooth_freq = 1;          % [LEGACY: 1] Frequency smoothing (bins)
cfg.coherence.mscohere.smooth_time = 1;          % [LEGACY: 1] Time smoothing (bins)

%% ============================================================================
%  COHERENCE PARAMETERS - FIELDTRIP [LEGACY: figure2_coherence_fieldtrip.m lines 488-587]
%  ============================================================================
%
%  TAPER OPTIONS:
%    'dpss' (default) - Multi-taper with Slepian sequences
%                     - Smoother estimates, lower variance
%                     - Frequency smoothing controlled by tapsmofrq
%                     - Good for broadband characterization
%
%    'hanning'        - Single Hanning window taper
%                     - Sharper frequency resolution (no smoothing)
%                     - Higher variance but better for detecting narrow-band peaks
%                     - Recommended for stimulation analysis (e.g., 40Hz, 135Hz)
%                     - tapsmofrq is ignored when using hanning
%
%  ============================================================================

cfg.coherence.fieldtrip.method = 'mtmfft';    % [LEGACY: 'mtmfft'] Method
cfg.coherence.fieldtrip.taper = 'hanning';    % 'dpss' (multi-taper) or 'hanning' (single taper)
cfg.coherence.fieldtrip.tapsmofrq = 2;        % Frequency smoothing (Hz) - only used with 'dpss'
cfg.coherence.fieldtrip.foi_min = 2;          % [LEGACY: 2] Minimum frequency (Hz)
cfg.coherence.fieldtrip.foi_max = 70;         % [LEGACY: 70] Maximum frequency (Hz)
cfg.coherence.fieldtrip.foi_step = 0.5;       % Frequency resolution (Hz) - higher than legacy 1 Hz

% Pseudo-trial parameters - CRITICAL: MUST MATCH LEGACY!
% Legacy (lines 488-489): EPOCH_LENGTH_SEC=1.0, EPOCH_OVERLAP_SEC=0.8
% This produces ~296 pseudo-trials for 60s recording (vs ~59 with 2s/50%)
cfg.coherence.fieldtrip.pseudotrial_length_sec = 1.0;    % [LEGACY: 1.0] seconds
cfg.coherence.fieldtrip.pseudotrial_overlap_sec = 0.8;   % [LEGACY: 0.8] seconds (80%)

% Time-resolved coherence parameters
% Legacy (lines 586-587): time_window_epochs=10, loops by 1 (not 2!)
cfg.coherence.fieldtrip.time_window_epochs = 10;    % [LEGACY: 10] Epochs per time point
cfg.coherence.fieldtrip.smooth_freq = 1;            % [LEGACY: 1] Frequency smoothing (bins)
cfg.coherence.fieldtrip.smooth_time = 1;            % [LEGACY: 1] Time smoothing (bins)

%% ============================================================================
%  PSD PARAMETERS
%  ============================================================================

cfg.psd.window_sec = 1.0;             % Window length (seconds)
cfg.psd.overlap_frac = 0.8;           % Overlap fraction
cfg.psd.freq_min = 2;                 % Minimum frequency (Hz)
cfg.psd.freq_max = 70;                % Maximum frequency (Hz)

%% ============================================================================
%  FREQUENCY BANDS FOR STATISTICS
%  ============================================================================

cfg.freq_bands.theta = [4, 8];        % Theta band (Hz)
cfg.freq_bands.alpha = [8, 12];       % Alpha band (Hz)
cfg.freq_bands.beta = [12, 30];       % Beta band (Hz)
cfg.freq_bands.gamma = [30, 70];      % Gamma band (Hz)

%% ============================================================================
%  PATH CONFIGURATION
%  ============================================================================

% Base paths -- sourced from the centralised config (config/lab_paths.m).
% Override per machine via config/paths_local.m. See config/README.md.
cfg_dir   = fileparts(mfilename('fullpath'));            % .../spectral_analysis/config
repo_root = fileparts(fileparts(cfg_dir));               % project root
addpath(fullfile(repo_root, 'config'));
lp = lab_paths();

cfg.paths.data_root   = lp.data_root;
cfg.paths.output_root = lp.spectral_output_root;
cfg.paths.code_root   = lp.spectral_code;

% Toolbox paths (added at runtime by run_spectral_pipeline)
cfg.paths.toolboxes   = lp.toolboxes;

%% ============================================================================
%  FIGURE PARAMETERS
%  ============================================================================

cfg.figure.dpi = 300;                 % Figure resolution
cfg.figure.format = 'png';            % Output format ('png', 'pdf', 'svg', 'fig')
cfg.figure.font_size_title = 14;      % Title font size
cfg.figure.font_size_label = 12;      % Axis label font size
cfg.figure.font_size_tick = 10;       % Tick label font size

% Colors
cfg.figure.colors.rest = [0.2, 0.4, 0.6];        % REST color (blue)
cfg.figure.colors.run = [0.8, 0.4, 0.2];         % RUN color (orange)
cfg.figure.colors.overall = [0.3, 0.3, 0.3];     % Overall color (gray)

%% ============================================================================
%  GROUP STATISTICS PARAMETERS
%  ============================================================================

cfg.stats.alpha_level = 0.05;         % Significance threshold
cfg.stats.use_fdr_correction = true;  % Apply FDR correction

% Cluster-based permutation (FieldTrip)
%
% CLUSTER_ALPHA: Cluster-forming threshold. Higher values are MORE SENSITIVE
%   because more frequency bins can join clusters. The cluster-level test
%   still properly controls false positives. 0.15-0.20 recommended for small N.
%
% TAIL: For directional hypotheses (e.g., RUN > REST for theta coherence),
%   use tail = 1 (one-tailed). This doubles statistical power.
%   Use tail = 0 (two-tailed) for non-directional hypotheses.
%
cfg.stats.cluster.alpha = 0.15;           % Cluster-forming threshold (more lenient for small N)
cfg.stats.cluster.pval = 0.05;            % Cluster significance threshold
cfg.stats.cluster.num_randomizations = 'all';  % 'all' for exact test with small N
cfg.stats.cluster.tail = 1;               % 1 = one-tailed (RUN > REST), 0 = two-tailed

%% ============================================================================
%  FIBER/CHANNEL SELECTION
%  ============================================================================

cfg.fiber_index = 1;                  % Which fiber to analyze (1-based)

%% ============================================================================
%  OUTPUT DATA STRUCTURE CONFIGURATION
%  ============================================================================

cfg.output.save_mat = true;           % Save .mat data files
cfg.output.save_csv = false;          % Save summary CSV (not implemented)
cfg.output.compress_mat = true;       % Use compression for .mat files

end

