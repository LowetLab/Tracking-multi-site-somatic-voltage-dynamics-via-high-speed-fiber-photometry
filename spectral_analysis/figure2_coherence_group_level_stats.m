%% ============================================================================
%  GROUP-LEVEL STATISTICS: LFP-GEVI Coherence and PSD Across Animals
%  ============================================================================
%  This script computes group-level statistics by loading POOLED data from
%  multiple animals and computing mean ± SEM with paired t-tests.
%
%  POOLING STRATEGIES:
%    1. Session-pooled: Use one representative session's pooled data per animal
%    2. Animal-level (configurable via GROUP_POOLING_LEVEL):
%       - 'animal_pooled': Spectra computed per session, then averaged
%       - 'animal_concatenated': All raw data concatenated, spectra computed once
%
%  TWO COHERENCE METHODS:
%    1. mscohere: Fourier-based coherence
%       - Standard paired t-tests + FDR correction
%       - Band-averaged statistics
%
%    2. FieldTrip: Multi-taper coherence
%       - Standard paired t-tests + FDR correction  
%       - Band-averaged statistics
%       - FieldTrip cluster-based permutation testing
%
%  Output saved to:
%    Spectral_data_outputs/{BEHAVIOR_MODE}/group_level/data/
%      mscohere_session_pooled.mat
%      mscohere_{GROUP_POOLING_LEVEL}.mat
%      fieldtrip_session_pooled.mat          (standard stats)
%      fieldtrip_{GROUP_POOLING_LEVEL}.mat   (standard stats)
%      fieldtrip_cluster_session_pooled.mat  (cluster stats)
%      fieldtrip_cluster_{GROUP_POOLING_LEVEL}.mat (cluster stats)
%  ============================================================================

% NOTE: Only clear if NOT called from pipeline (to preserve cfg variable)
if ~exist('cfg', 'var')
    close all; clear; clc;
else
    % Called from pipeline - preserve cfg, just close figures
    close all; clc;
    fprintf('  Inheriting configuration from pipeline...\n');
end

%% ============================================================================
%  ADD REQUIRED PATHS (for FieldTrip cluster statistics and animal database)
%  ============================================================================

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'config'));   % animal_session_database.m
addpath(fullfile(script_dir, 'core'));     % if needed by pipeline
% External toolboxes via the centralised config (config/lab_paths.m).
addpath(fullfile(fileparts(script_dir), 'config'));  % project config/
setup_lab_paths();

% Check for FieldTrip (needed for cluster stats)
if exist('ft_defaults', 'file')
    ft_defaults;
    HAS_FIELDTRIP = true;
    fprintf('FieldTrip available - cluster statistics will be computed.\n');
else
    HAS_FIELDTRIP = false;
    fprintf('WARNING: FieldTrip not found - cluster statistics will be skipped.\n');
end

%% ============================================================================
%  USER CONFIGURATION
%  ============================================================================

BASE_OUTPUT_DIR = fullfile(lab_paths().figures_root, 'Spectral_data_outputs');

% Behavior classification mode: 'standard' or 'clear'
% If called from pipeline, inherit from cfg; otherwise use this default
if exist('cfg', 'var') && isfield(cfg, 'behavior') && isfield(cfg.behavior, 'classification_mode')
    BEHAVIOR_MODE = cfg.behavior.classification_mode;
else
    BEHAVIOR_MODE = 'clear';  % ← CHANGE THIS to match your analysis
end

% Output folder suffix (for artifact handling modes)
% Applied to ROOT directory: Spectral_data_outputs → Spectral_data_outputs{SUFFIX}
%
% ARTIFACT HANDLING MODES (set in run_spectral_pipeline.m):
%   ''                   - No artifact handling (original data)
%   '_artifact_excluded' - Trial exclusion: Trials with >30% contamination skipped
%   '_artifact_cleaned'  - Data cleaning: Artifact segments removed from all trials
%
% Inherit from cfg if available, otherwise use this default
if exist('cfg', 'var') && isfield(cfg, 'output_folder_suffix')
    OUTPUT_SUFFIX = cfg.output_folder_suffix;
else
    OUTPUT_SUFFIX = '_artifact_cleaned';  % ← CHANGE: '', '_artifact_excluded', or '_artifact_cleaned'
end

% Apply suffix to ROOT directory (not to behavior mode)
% Structure: {BASE_OUTPUT_DIR}{SUFFIX}/{BEHAVIOR_MODE}/...
% Examples: 
%   Spectral_data_outputs/clear/group_level/data/
%   Spectral_data_outputs_artifact_excluded/clear/group_level/data/
%   Spectral_data_outputs_artifact_cleaned/clear/group_level/data/
if ~isempty(OUTPUT_SUFFIX)
    EFFECTIVE_OUTPUT_DIR = [BASE_OUTPUT_DIR, OUTPUT_SUFFIX];
    fprintf('  Using output suffix: %s\n', OUTPUT_SUFFIX);
else
    EFFECTIVE_OUTPUT_DIR = BASE_OUTPUT_DIR;
end

GROUP_OUTPUT_DIR = fullfile(EFFECTIVE_OUTPUT_DIR, BEHAVIOR_MODE, 'group_level', 'data');

% Coherence methods to process
METHODS = {'mscohere', 'fieldtrip'};  % Can be {'mscohere'}, {'fieldtrip'}, or both

% ==============================================================================
%  GROUP-LEVEL POOLING CONFIGURATION
% ==============================================================================
% Which animal-level pooling to use for group statistics:
%   'animal_pooled'       - Spectra computed per session, then averaged (default)
%   'animal_concatenated' - All raw data concatenated across sessions, spectra computed once
%
% animal_concatenated provides a single spectrum from more data points but requires
% proper edge artifact handling for FieldTrip analysis. Both levels compute
% one data point per animal for group statistics.
GROUP_POOLING_LEVEL = 'animal_concatenated';  % 'animal_pooled' or 'animal_concatenated'

% Animal configuration: load from animal_session_database (single source of truth)
% ANIMALS is built as {mouse_id, representative_session} for each animal.
% Representative session for session-pooled stats is the FIRST session in the database.
% When called from run_spectral_pipeline.m, ANIMALS_TO_PROCESS restricts to that subset.
animals_db = animal_session_database();
valid = arrayfun(@(x) ~isempty(x.sessions), animals_db);
animals_db = animals_db(valid);
ANIMALS = cell(length(animals_db), 2);
for a = 1:length(animals_db)
    ANIMALS{a, 1} = animals_db(a).mouse_id;
    % First session as representative for session-pooled group-level stats
    ANIMALS{a, 2} = animals_db(a).sessions(1).session_id;
end

% If pipeline passed ANIMALS_TO_PROCESS, use only those animals (so group-level N matches pipeline)
if exist('ANIMALS_TO_PROCESS', 'var') && ~isempty(ANIMALS_TO_PROCESS)
    keep = cellfun(@(x) ismember(x, ANIMALS_TO_PROCESS), ANIMALS(:,1));
    ANIMALS = ANIMALS(keep, :);
    fprintf('  Group-level: using %d animals from pipeline: %s\n', size(ANIMALS,1), strjoin(ANIMALS_TO_PROCESS, ', '));
end

num_animals = size(ANIMALS, 1);
if num_animals == 0
    error('Group-level: no animals to process. Check ANIMALS_TO_PROCESS matches animal_session_database mouse_id entries.');
end

% Statistical parameters
FREQ_BANDS = struct();
FREQ_BANDS.theta = [5, 9];
FREQ_BANDS.alpha = [9, 12];
FREQ_BANDS.beta = [15, 30];
FREQ_BANDS.gamma = [30, 60];

ALPHA_LEVEL = 0.05;
USE_FDR_CORRECTION = true;
FREQ_MIN = 2;
FREQ_MAX = 70;

% Cluster-based permutation parameters (FieldTrip only)
% -------------------------------------------------------------------------
% OPTIMIZATION NOTES for small N (≤7 animals):
%   - clusteralpha: Higher values (0.15-0.20) are MORE SENSITIVE to detect
%     effects because more frequency bins can join clusters. The cluster-level
%     alpha still controls false positives correctly.
%   - Use 'all' permutations when N is small for exact p-values
%   - ONE-TAILED TEST: Use when you have a directional hypothesis (e.g., 
%     theta coherence is HIGHER during RUN than REST). This doubles power.
%     tail=1 tests for RUN > REST, tail=-1 tests for REST > RUN, tail=0 is two-tailed.
% -------------------------------------------------------------------------
CLUSTER_ALPHA = 0.15;           % More liberal cluster-forming threshold (recommended for small N)
CLUSTER_PVAL = 0.05;            % Significance threshold for clusters
NUM_RANDOMIZATIONS = 'all';     % Use ALL permutations for exact p-values (important for small N)
CLUSTER_TAIL = 1;               % 1 = one-tailed (RUN > REST), 0 = two-tailed, -1 = REST > RUN

%% ============================================================================
%  CREATE OUTPUT DIRECTORY
%  ============================================================================

if ~exist(GROUP_OUTPUT_DIR, 'dir')
    mkdir(GROUP_OUTPUT_DIR);
    fprintf('Created: %s\n', GROUP_OUTPUT_DIR);
end

%% ============================================================================
%  PROCESS EACH METHOD
%  ============================================================================

for method_idx = 1:length(METHODS)
    
    METHOD = METHODS{method_idx};
    
    fprintf('\n');
    fprintf('================================================================\n');
    fprintf('  PROCESSING METHOD: %s\n', upper(METHOD));
    fprintf('================================================================\n\n');
    
    %% ========================================================================
    %  LOAD SESSION-POOLED DATA
    %  ========================================================================
    
    fprintf('--- Loading Session-Pooled Data (%s) ---\n', METHOD);
    
    session_data = struct();
    valid_session = 0;
    
    for a = 1:num_animals
        mouse_id = ANIMALS{a, 1};
        session_id = ANIMALS{a, 2};
        
        % Path structure: {BASE_OUTPUT_DIR}{SUFFIX}/{BEHAVIOR_MODE}/session_pooled/{mouse_id}/{session_id}/data/{method}.mat
        data_path = fullfile(EFFECTIVE_OUTPUT_DIR, BEHAVIOR_MODE, 'session_pooled', mouse_id, session_id, 'data', ...
            sprintf('%s.mat', METHOD));
        
        if ~exist(data_path, 'file')
            fprintf('  WARNING: Not found: %s/%s\n', mouse_id, session_id);
            continue;
        end
        
        data = load(data_path);
        
        % Handle field name inconsistency: pooledtrials uses 'coh_spectrum', animalpooled uses 'coh_overall'
        if isfield(data, 'coh_spectrum') && ~isfield(data, 'coh_overall')
            data.coh_overall = data.coh_spectrum;
        end
        
        if ~isfield(data, 'coh_overall') || ~isfield(data, 'freq')
            fprintf('  WARNING: Missing fields in %s\n', data_path);
            continue;
        end
        
        valid_session = valid_session + 1;
        session_data(valid_session).mouse_id = mouse_id;
        session_data(valid_session).session_id = session_id;
        session_data(valid_session).data = data;
        
        % Handle missing fields gracefully
        pct_rest = 0; pct_run = 0;
        if isfield(data, 'pct_rest'), pct_rest = data.pct_rest; end
        if isfield(data, 'pct_run'), pct_run = data.pct_run; end
        
        fprintf('  Loaded %s/%s: %.1f%% rest, %.1f%% run\n', ...
            mouse_id, session_id, pct_rest, pct_run);
    end
    
    fprintf('  Valid: %d / %d animals\n\n', valid_session, num_animals);
    
    %% ========================================================================
    %  LOAD ANIMAL-LEVEL DATA (animal_pooled OR animal_concatenated)
    %  ========================================================================
    
    fprintf('--- Loading %s Data (%s) ---\n', upper(strrep(GROUP_POOLING_LEVEL, '_', '-')), METHOD);
    
    animal_data = struct();
    valid_animal = 0;
    
    for a = 1:num_animals
        mouse_id = ANIMALS{a, 1};
        
        % Path structure: {BASE_OUTPUT_DIR}{SUFFIX}/{BEHAVIOR_MODE}/{GROUP_POOLING_LEVEL}/{mouse_id}/data/{method}.mat
        data_path = fullfile(EFFECTIVE_OUTPUT_DIR, BEHAVIOR_MODE, GROUP_POOLING_LEVEL, mouse_id, 'data', ...
            sprintf('%s.mat', METHOD));
        
        if ~exist(data_path, 'file')
            fprintf('  WARNING: Not found: %s\n', mouse_id);
            continue;
        end
        
        data = load(data_path);
        
        % Handle field name inconsistency: pooledtrials uses 'coh_spectrum', animalpooled uses 'coh_overall'
        if isfield(data, 'coh_spectrum') && ~isfield(data, 'coh_overall')
            data.coh_overall = data.coh_spectrum;
        end
        
        if ~isfield(data, 'coh_overall') || ~isfield(data, 'freq')
            fprintf('  WARNING: Missing fields\n');
            continue;
        end
        
        valid_animal = valid_animal + 1;
        animal_data(valid_animal).mouse_id = mouse_id;
        animal_data(valid_animal).data = data;
        
        % Handle missing fields gracefully
        if isfield(data, 'num_sessions')
            n_sess = data.num_sessions;
        else
            n_sess = NaN;  % Unknown - older data files may not have this
        end
        pct_rest = 0; pct_run = 0;
        if isfield(data, 'pct_rest'), pct_rest = data.pct_rest; end
        if isfield(data, 'pct_run'), pct_run = data.pct_run; end
        
        fprintf('  Loaded %s: %d sessions, %.1f%% rest, %.1f%% run\n', ...
            mouse_id, n_sess, pct_rest, pct_run);
    end
    
    fprintf('  Valid: %d / %d animals\n\n', valid_animal, num_animals);
    
    %% ========================================================================
    %  COMPILE DATA MATRICES AND COMPUTE STATISTICS
    %  ========================================================================
    
    % Shared config for compute_and_save_pooled_stats (below): identical for
    % both the session-pooled and animal-pooled calls in this METHOD iteration.
    pooled_stats_cfg = struct();
    pooled_stats_cfg.METHOD = METHOD;
    pooled_stats_cfg.HAS_FIELDTRIP = HAS_FIELDTRIP;
    pooled_stats_cfg.ALPHA_LEVEL = ALPHA_LEVEL;
    pooled_stats_cfg.USE_FDR_CORRECTION = USE_FDR_CORRECTION;
    pooled_stats_cfg.FREQ_BANDS = FREQ_BANDS;
    pooled_stats_cfg.CLUSTER_ALPHA = CLUSTER_ALPHA;
    pooled_stats_cfg.CLUSTER_PVAL = CLUSTER_PVAL;
    pooled_stats_cfg.CLUSTER_TAIL = CLUSTER_TAIL;
    pooled_stats_cfg.NUM_RANDOMIZATIONS = NUM_RANDOMIZATIONS;
    pooled_stats_cfg.FREQ_MIN = FREQ_MIN;
    pooled_stats_cfg.FREQ_MAX = FREQ_MAX;
    pooled_stats_cfg.GROUP_OUTPUT_DIR = GROUP_OUTPUT_DIR;

    % --- SESSION-POOLED STATISTICS ---
    compute_and_save_pooled_stats(session_data, valid_session, 'session_pooled', pooled_stats_cfg, true);

    % --- ANIMAL-POOLED STATISTICS ---
    compute_and_save_pooled_stats(animal_data, valid_animal, GROUP_POOLING_LEVEL, pooled_stats_cfg, false);
end

fprintf('\n================================================================\n');
fprintf('  GROUP-LEVEL STATISTICS COMPLETE\n');
fprintf('================================================================\n');
fprintf('Output directory: %s\n', GROUP_OUTPUT_DIR);

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================

function compute_and_save_pooled_stats(pooled_data, valid_count, pooling_type, cfg, save_session_ids)
%COMPUTE_AND_SAVE_POOLED_STATS  Compute + save REST/RUN coherence & PSD stats
%(standard, plus FieldTrip cluster-permutation stats when applicable) for one
%pooling level, and write the corresponding .mat file(s).
%
%   pooled_data       : struct array with fields .data, .mouse_id
%                       (session-level entries also need .session_id when
%                       save_session_ids is true)
%   valid_count       : number of valid entries in pooled_data
%   pooling_type      : 'session_pooled', or GROUP_POOLING_LEVEL's value for
%                       the animal-level call -- used for BOTH output
%                       filenames and status messages
%   cfg               : struct with METHOD, HAS_FIELDTRIP, ALPHA_LEVEL,
%                       USE_FDR_CORRECTION, FREQ_BANDS, CLUSTER_ALPHA,
%                       CLUSTER_PVAL, CLUSTER_TAIL, NUM_RANDOMIZATIONS,
%                       FREQ_MIN, FREQ_MAX, GROUP_OUTPUT_DIR
%   save_session_ids  : true to also save out.session_ids (session-pooled
%                       call only)
%
%   Replaces what used to be two ~250-line byte-for-byte-parallel blocks
%   (session-pooled vs animal-pooled), one per pooling level.

label_display = upper(strrep(pooling_type, '_', '-'));

if valid_count < 2
    fprintf('  Insufficient data for %s statistics\n\n', strrep(pooling_type, '_', '-'));
    return;
end

fprintf('--- Computing %s Statistics ---\n', label_display);

% Get frequency axis
freq = pooled_data(1).data.freq(:);
num_freq = length(freq);

% Check for PSD
has_psd = isfield(pooled_data(1).data, 'psd_freq') && ~isempty(pooled_data(1).data.psd_freq);
if has_psd
    freq_psd = pooled_data(1).data.psd_freq(:);
    num_freq_psd = length(freq_psd);
end

% Compile matrices
coh_overall = zeros(num_freq, valid_count);
coh_rest = zeros(num_freq, valid_count);
coh_run = zeros(num_freq, valid_count);
has_restrun = false(1, valid_count);

if has_psd
    psd_lfp_rest = zeros(num_freq_psd, valid_count);
    psd_lfp_run = zeros(num_freq_psd, valid_count);
    psd_gevi_rest = zeros(num_freq_psd, valid_count);
    psd_gevi_run = zeros(num_freq_psd, valid_count);
end

for a = 1:valid_count
    d = pooled_data(a).data;

    coh_overall(:, a) = interp_to_freq(d.coh_overall, d.freq, freq);

    % Check dimensions match before interpolation
    has_rest = ~isempty(d.coh_rest) && length(d.coh_rest) == length(d.freq);
    has_run = ~isempty(d.coh_run) && length(d.coh_run) == length(d.freq);
    if has_rest && has_run
        has_restrun(a) = true;
        % Use d.freq for both rest and run (same frequency axis)
        coh_rest(:, a) = interp_to_freq(d.coh_rest, d.freq, freq);
        coh_run(:, a) = interp_to_freq(d.coh_run, d.freq, freq);
    end

    if has_psd && isfield(d, 'psd_freq') && ~isempty(d.psd_freq)
        % Use d.psd_freq for all PSD (same frequency axis)
        % Check each field separately - some may be empty if no REST or RUN
        if isfield(d, 'psd_lfp_rest') && ~isempty(d.psd_lfp_rest) && length(d.psd_lfp_rest) == length(d.psd_freq)
            psd_lfp_rest(:, a) = interp_to_freq(d.psd_lfp_rest, d.psd_freq, freq_psd);
        end
        if isfield(d, 'psd_lfp_run') && ~isempty(d.psd_lfp_run) && length(d.psd_lfp_run) == length(d.psd_freq)
            psd_lfp_run(:, a) = interp_to_freq(d.psd_lfp_run, d.psd_freq, freq_psd);
        end
        if isfield(d, 'psd_gevi_rest') && ~isempty(d.psd_gevi_rest) && length(d.psd_gevi_rest) == length(d.psd_freq)
            psd_gevi_rest(:, a) = interp_to_freq(d.psd_gevi_rest, d.psd_freq, freq_psd);
        end
        if isfield(d, 'psd_gevi_run') && ~isempty(d.psd_gevi_run) && length(d.psd_gevi_run) == length(d.psd_freq)
            psd_gevi_run(:, a) = interp_to_freq(d.psd_gevi_run, d.psd_freq, freq_psd);
        end
    end
end

% Compute statistics
idx_rr = find(has_restrun);

pooled_stats = struct();
pooled_stats.overall = compute_mean_sem(coh_overall);

if length(idx_rr) >= 2
    pooled_stats.coherence = compute_paired_stats(coh_rest(:, idx_rr), coh_run(:, idx_rr), cfg.ALPHA_LEVEL, cfg.USE_FDR_CORRECTION);
    pooled_stats.coherence.band_stats = compute_band_stats(coh_rest(:, idx_rr), coh_run(:, idx_rr), freq, cfg.FREQ_BANDS, cfg.ALPHA_LEVEL);
    fprintf('  Coherence: %d sig freqs\n', sum(pooled_stats.coherence.sig_mask));
end

if has_psd && length(idx_rr) >= 2
    pooled_stats.psd_lfp = compute_paired_stats(psd_lfp_rest(:, idx_rr), psd_lfp_run(:, idx_rr), cfg.ALPHA_LEVEL, cfg.USE_FDR_CORRECTION);
    pooled_stats.psd_lfp.band_stats = compute_band_stats(psd_lfp_rest(:, idx_rr), psd_lfp_run(:, idx_rr), freq_psd, cfg.FREQ_BANDS, cfg.ALPHA_LEVEL);
    pooled_stats.psd_gevi = compute_paired_stats(psd_gevi_rest(:, idx_rr), psd_gevi_run(:, idx_rr), cfg.ALPHA_LEVEL, cfg.USE_FDR_CORRECTION);
    pooled_stats.psd_gevi.band_stats = compute_band_stats(psd_gevi_rest(:, idx_rr), psd_gevi_run(:, idx_rr), freq_psd, cfg.FREQ_BANDS, cfg.ALPHA_LEVEL);
    fprintf('  PSD LFP: %d sig freqs, PSD GEVI: %d sig freqs\n', sum(pooled_stats.psd_lfp.sig_mask), sum(pooled_stats.psd_gevi.sig_mask));
end

% Save stats (standard)
out = create_output_struct(cfg.METHOD, pooling_type, valid_count, ...
    {pooled_data.mouse_id}, freq, pooled_stats, has_psd, freq_psd);
if save_session_ids
    out.session_ids = {pooled_data.session_id};
end

% Naming: {method}_{pooling_type}.mat
save_path = fullfile(cfg.GROUP_OUTPUT_DIR, sprintf('%s_%s.mat', cfg.METHOD, pooling_type));
save(save_path, '-struct', 'out', '-v7');
fprintf('  Saved: %s\n\n', save_path);

% =========================================================================
%  CLUSTER-BASED PERMUTATION STATS (FieldTrip only)
% =========================================================================
if strcmpi(cfg.METHOD, 'fieldtrip') && cfg.HAS_FIELDTRIP && length(idx_rr) >= 2
    fprintf('--- Computing Cluster-Based Permutation Statistics (%s) ---\n', label_display);
    fprintf('    Using %s test (tail=%d)\n', ternary(cfg.CLUSTER_TAIL==0, 'two-tailed', 'one-tailed'), cfg.CLUSTER_TAIL);

    cluster_out = compute_cluster_stats(coh_rest(:, idx_rr), coh_run(:, idx_rr), ...
        freq, 'coherence', cfg.CLUSTER_ALPHA, cfg.CLUSTER_PVAL, cfg.NUM_RANDOMIZATIONS, cfg.FREQ_MIN, cfg.FREQ_MAX, cfg.CLUSTER_TAIL);

    if has_psd
        cluster_out.psd_lfp_cluster = compute_cluster_stats(psd_lfp_rest(:, idx_rr), psd_lfp_run(:, idx_rr), ...
            freq_psd, 'psd_lfp', cfg.CLUSTER_ALPHA, cfg.CLUSTER_PVAL, cfg.NUM_RANDOMIZATIONS, cfg.FREQ_MIN, cfg.FREQ_MAX, cfg.CLUSTER_TAIL);
        cluster_out.psd_gevi_cluster = compute_cluster_stats(psd_gevi_rest(:, idx_rr), psd_gevi_run(:, idx_rr), ...
            freq_psd, 'psd_gevi', cfg.CLUSTER_ALPHA, cfg.CLUSTER_PVAL, cfg.NUM_RANDOMIZATIONS, cfg.FREQ_MIN, cfg.FREQ_MAX, cfg.CLUSTER_TAIL);
    end

    % Create cluster output structure
    cluster_save = struct();
    cluster_save.method = sprintf('fieldtrip_cluster_%s', pooling_type);
    cluster_save.stats_type = 'cluster_permutation';
    cluster_save.num_animals = length(idx_rr);
    cluster_save.animal_ids = {pooled_data(idx_rr).mouse_id};
    cluster_save.cluster_alpha = cfg.CLUSTER_ALPHA;
    cluster_save.cluster_pval = cfg.CLUSTER_PVAL;
    cluster_save.cluster_tail = cfg.CLUSTER_TAIL;  % 0=two-tailed, 1=RUN>REST, -1=REST>RUN
    cluster_save.num_randomizations = cfg.NUM_RANDOMIZATIONS;
    cluster_save.analysis_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
    cluster_save.freq = freq(:);

    % Copy overall coherence (for top row plots)
    cluster_save.coherence_overall = struct();
    cluster_save.coherence_overall.mean = pooled_stats.overall.mean(:);
    cluster_save.coherence_overall.sem = pooled_stats.overall.sem(:);

    % Copy means/SEMs from standard stats
    cluster_save.coherence = struct();
    cluster_save.coherence.REST = struct();
    cluster_save.coherence.REST.mean = pooled_stats.coherence.rest_mean(:);
    cluster_save.coherence.REST.sem = pooled_stats.coherence.rest_sem(:);
    cluster_save.coherence.RUN = struct();
    cluster_save.coherence.RUN.mean = pooled_stats.coherence.run_mean(:);
    cluster_save.coherence.RUN.sem = pooled_stats.coherence.run_sem(:);
    cluster_save.coherence.cluster_stats = cluster_out;
    if isfield(pooled_stats.coherence, 'band_stats')
        cluster_save.coherence.band_stats = pooled_stats.coherence.band_stats;
    end

    if has_psd
        cluster_save.freq_psd = freq_psd(:);
        cluster_save.psd_lfp = struct();
        cluster_save.psd_lfp.REST = struct();
        cluster_save.psd_lfp.REST.mean = pooled_stats.psd_lfp.rest_mean(:);
        cluster_save.psd_lfp.REST.sem = pooled_stats.psd_lfp.rest_sem(:);
        cluster_save.psd_lfp.RUN = struct();
        cluster_save.psd_lfp.RUN.mean = pooled_stats.psd_lfp.run_mean(:);
        cluster_save.psd_lfp.RUN.sem = pooled_stats.psd_lfp.run_sem(:);
        cluster_save.psd_lfp.cluster_stats = cluster_out.psd_lfp_cluster;
        cluster_save.psd_lfp.units = 'dB (re 1 uV^2/Hz)';
        if isfield(pooled_stats.psd_lfp, 'band_stats')
            cluster_save.psd_lfp.band_stats = pooled_stats.psd_lfp.band_stats;
        end

        cluster_save.psd_gevi = struct();
        cluster_save.psd_gevi.REST = struct();
        cluster_save.psd_gevi.REST.mean = pooled_stats.psd_gevi.rest_mean(:);
        cluster_save.psd_gevi.REST.sem = pooled_stats.psd_gevi.rest_sem(:);
        cluster_save.psd_gevi.RUN = struct();
        cluster_save.psd_gevi.RUN.mean = pooled_stats.psd_gevi.run_mean(:);
        cluster_save.psd_gevi.RUN.sem = pooled_stats.psd_gevi.run_sem(:);
        cluster_save.psd_gevi.cluster_stats = cluster_out.psd_gevi_cluster;
        cluster_save.psd_gevi.units = 'dB (re 1 (dF/F)^2/Hz)';
        if isfield(pooled_stats.psd_gevi, 'band_stats')
            cluster_save.psd_gevi.band_stats = pooled_stats.psd_gevi.band_stats;
        end
    end

    % Naming: fieldtrip_cluster_{pooling_type}.mat
    cluster_save_path = fullfile(cfg.GROUP_OUTPUT_DIR, sprintf('fieldtrip_cluster_%s.mat', pooling_type));
    save(cluster_save_path, '-struct', 'cluster_save', '-v7');
    fprintf('  Saved cluster stats: %s\n\n', cluster_save_path);
end
end

function y_interp = interp_to_freq(y, x_orig, x_target)
% Interpolate data to target frequency axis
y = y(:);
x_orig = x_orig(:);
x_target = x_target(:);

if length(y) == length(x_target)
    y_interp = y;
else
    y_interp = interp1(x_orig, y, x_target, 'linear', 'extrap');
end
end

function stats = compute_mean_sem(data_matrix)
% Compute mean and SEM across columns
stats.mean = mean(data_matrix, 2);
stats.std = std(data_matrix, 0, 2);
stats.sem = stats.std / sqrt(size(data_matrix, 2));
end

function stats = compute_paired_stats(rest_matrix, run_matrix, alpha, use_fdr)
% Compute paired t-test statistics per frequency
[num_freq, num_subj] = size(rest_matrix);

stats.rest_mean = mean(rest_matrix, 2);
stats.rest_std = std(rest_matrix, 0, 2);
stats.rest_sem = stats.rest_std / sqrt(num_subj);

stats.run_mean = mean(run_matrix, 2);
stats.run_std = std(run_matrix, 0, 2);
stats.run_sem = stats.run_std / sqrt(num_subj);

stats.pvals = zeros(num_freq, 1);
stats.tvals = zeros(num_freq, 1);

for f = 1:num_freq
    if num_subj >= 2
        [~, p, ~, stat] = ttest(rest_matrix(f, :), run_matrix(f, :));
        stats.pvals(f) = p;
        stats.tvals(f) = stat.tstat;
    else
        stats.pvals(f) = 1;
        stats.tvals(f) = 0;
    end
end

stats.pvals(isnan(stats.pvals)) = 1;

if use_fdr
    stats.pvals_fdr = fdr_bh(stats.pvals);
    stats.sig_mask = stats.pvals_fdr < alpha;
else
    stats.pvals_fdr = stats.pvals;
    stats.sig_mask = stats.pvals < alpha;
end
end

function pvals_fdr = fdr_bh(pvals)
% Benjamini-Hochberg FDR correction
n = length(pvals);
[sorted_pvals, sort_idx] = sort(pvals);

adjusted = zeros(n, 1);
adjusted(n) = sorted_pvals(n);

for i = (n-1):-1:1
    adjusted(i) = min(adjusted(i+1), sorted_pvals(i) * n / i);
end

[~, unsort_idx] = sort(sort_idx);
pvals_fdr = min(adjusted(unsort_idx), 1);
end

function band_stats = compute_band_stats(rest_matrix, run_matrix, freq, freq_bands, alpha)
% Compute band-averaged statistics with effect sizes
band_names = fieldnames(freq_bands);
band_stats = struct();

for b = 1:length(band_names)
    band = band_names{b};
    band_range = freq_bands.(band);
    band_idx = freq >= band_range(1) & freq <= band_range(2);
    
    if ~any(band_idx)
        continue;
    end
    
    rest_band = mean(rest_matrix(band_idx, :), 1);  % Average within band per subject
    run_band = mean(run_matrix(band_idx, :), 1);
    num_subj = size(rest_matrix, 2);
    
    band_stats.(band) = struct();
    band_stats.(band).freq_range = band_range;
    band_stats.(band).n_subjects = num_subj;
    
    % Descriptive statistics
    band_stats.(band).rest_mean = mean(rest_band);
    band_stats.(band).rest_std = std(rest_band);
    band_stats.(band).rest_sem = band_stats.(band).rest_std / sqrt(num_subj);
    band_stats.(band).run_mean = mean(run_band);
    band_stats.(band).run_std = std(run_band);
    band_stats.(band).run_sem = band_stats.(band).run_std / sqrt(num_subj);
    
    % Per-subject values (for Python plotting)
    band_stats.(band).rest_values = rest_band(:);
    band_stats.(band).run_values = run_band(:);
    
    % Paired t-test (two-tailed)
    [~, p_two, ~, stat] = ttest(rest_band, run_band);
    band_stats.(band).pval_twotailed = p_two;
    band_stats.(band).tval = stat.tstat;
    band_stats.(band).df = stat.df;
    band_stats.(band).sig_twotailed = p_two < alpha;
    
    % One-tailed p-value (for directional hypotheses)
    % If t > 0, REST > RUN; if t < 0, RUN > REST
    if stat.tstat > 0
        band_stats.(band).pval_onetailed_rest_greater = p_two / 2;
        band_stats.(band).pval_onetailed_run_greater = 1 - p_two / 2;
    else
        band_stats.(band).pval_onetailed_rest_greater = 1 - p_two / 2;
        band_stats.(band).pval_onetailed_run_greater = p_two / 2;
    end
    
    % Effect size: Cohen's d for paired samples
    % d = mean(difference) / std(difference)
    diff = rest_band - run_band;
    cohens_d = mean(diff) / std(diff);
    band_stats.(band).cohens_d = cohens_d;
    
    % Effect size interpretation
    % |d| < 0.2 = negligible, 0.2-0.5 = small, 0.5-0.8 = medium, > 0.8 = large
    if abs(cohens_d) < 0.2
        band_stats.(band).effect_size = 'negligible';
    elseif abs(cohens_d) < 0.5
        band_stats.(band).effect_size = 'small';
    elseif abs(cohens_d) < 0.8
        band_stats.(band).effect_size = 'medium';
    else
        band_stats.(band).effect_size = 'large';
    end
end
end

function out = create_output_struct(method, pooling_type, num_animals, animal_ids, freq, stats, has_psd, freq_psd)
% Create Python-compatible output structure

out = struct();

% Metadata
out.method = sprintf('%s_%s', method, pooling_type);
out.stats_type = 'paired_ttest_fdr';
out.num_animals = num_animals;
out.alpha_level = 0.05;
out.fdr_corrected = true;
out.analysis_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
out.animal_ids = animal_ids;

% Frequency axis
out.freq = freq(:);
if has_psd
    out.freq_psd = freq_psd(:);
end

% Overall coherence
out.coherence_overall = struct();
out.coherence_overall.mean = stats.overall.mean(:);
out.coherence_overall.sem = stats.overall.sem(:);

% Rest vs Run coherence
if isfield(stats, 'coherence')
    out.coherence = struct();
    out.coherence.REST = struct();
    out.coherence.REST.mean = stats.coherence.rest_mean(:);
    out.coherence.REST.sem = stats.coherence.rest_sem(:);
    out.coherence.RUN = struct();
    out.coherence.RUN.mean = stats.coherence.run_mean(:);
    out.coherence.RUN.sem = stats.coherence.run_sem(:);
    out.coherence.stats = struct();
    out.coherence.stats.pvals = stats.coherence.pvals(:);
    out.coherence.stats.pvals_fdr = stats.coherence.pvals_fdr(:);
    out.coherence.stats.sig_mask = stats.coherence.sig_mask(:);
    out.coherence.stats.tvals = stats.coherence.tvals(:);
    if isfield(stats.coherence, 'band_stats')
        out.coherence.band_stats = stats.coherence.band_stats;
    end
end

% PSD
if has_psd && isfield(stats, 'psd_lfp')
    out.psd_lfp = struct();
    out.psd_lfp.REST = struct();
    out.psd_lfp.REST.mean = stats.psd_lfp.rest_mean(:);
    out.psd_lfp.REST.sem = stats.psd_lfp.rest_sem(:);
    out.psd_lfp.RUN = struct();
    out.psd_lfp.RUN.mean = stats.psd_lfp.run_mean(:);
    out.psd_lfp.RUN.sem = stats.psd_lfp.run_sem(:);
    out.psd_lfp.stats = struct();
    out.psd_lfp.stats.pvals = stats.psd_lfp.pvals(:);
    out.psd_lfp.stats.pvals_fdr = stats.psd_lfp.pvals_fdr(:);
    out.psd_lfp.stats.sig_mask = stats.psd_lfp.sig_mask(:);
    out.psd_lfp.stats.tvals = stats.psd_lfp.tvals(:);
    out.psd_lfp.units = 'dB (re 1 uV^2/Hz)';
    if isfield(stats.psd_lfp, 'band_stats')
        out.psd_lfp.band_stats = stats.psd_lfp.band_stats;
    end
end

if has_psd && isfield(stats, 'psd_gevi')
    out.psd_gevi = struct();
    out.psd_gevi.REST = struct();
    out.psd_gevi.REST.mean = stats.psd_gevi.rest_mean(:);
    out.psd_gevi.REST.sem = stats.psd_gevi.rest_sem(:);
    out.psd_gevi.RUN = struct();
    out.psd_gevi.RUN.mean = stats.psd_gevi.run_mean(:);
    out.psd_gevi.RUN.sem = stats.psd_gevi.run_sem(:);
    out.psd_gevi.stats = struct();
    out.psd_gevi.stats.pvals = stats.psd_gevi.pvals(:);
    out.psd_gevi.stats.pvals_fdr = stats.psd_gevi.pvals_fdr(:);
    out.psd_gevi.stats.sig_mask = stats.psd_gevi.sig_mask(:);
    out.psd_gevi.stats.tvals = stats.psd_gevi.tvals(:);
    out.psd_gevi.units = 'dB (re 1 (dF/F)^2/Hz)';
    if isfield(stats.psd_gevi, 'band_stats')
        out.psd_gevi.band_stats = stats.psd_gevi.band_stats;
    end
end
end

function results = compute_cluster_stats(rest_matrix, run_matrix, freq, data_label, cluster_alpha, cluster_pval, num_rand, freq_min, freq_max, tail)
% COMPUTE_CLUSTER_STATS Perform FieldTrip cluster-based permutation testing
%
% This function performs cluster-based permutation testing using FieldTrip
% to compare REST vs RUN conditions in a within-subject (paired) design.
%
% Inputs:
%   rest_matrix    - [num_freq × num_subjects] REST data
%   run_matrix     - [num_freq × num_subjects] RUN data
%   freq           - frequency axis
%   data_label     - label for the data type (for FieldTrip structure)
%   cluster_alpha  - threshold for cluster formation
%   cluster_pval   - significance threshold for clusters
%   num_rand       - number of randomizations
%   freq_min/max   - frequency range for analysis
%   tail           - 0 = two-tailed, 1 = one-tailed (RUN > REST), -1 = one-tailed (REST > RUN)
%
% Output:
%   results - struct with cluster statistics

% Default to two-tailed if not specified
if nargin < 10 || isempty(tail)
    tail = 0;
end

results = struct();
[num_freq, num_subj] = size(rest_matrix);

% Construct FieldTrip structures
% -------------------------------------------------------------------------
% For frequency-domain cluster statistics, FieldTrip expects:
%   data.powspctrm - [Nsubj × Nchan × Nfreq]
%   data.freq      - frequency axis
%   data.label     - channel labels
%   data.dimord    - dimension order string
% -------------------------------------------------------------------------

ft_rest = struct();
ft_rest.label = {data_label};
ft_rest.freq = freq(:)';
ft_rest.dimord = 'subj_chan_freq';
ft_rest.powspctrm = zeros(num_subj, 1, num_freq);

ft_run = struct();
ft_run.label = {data_label};
ft_run.freq = freq(:)';
ft_run.dimord = 'subj_chan_freq';
ft_run.powspctrm = zeros(num_subj, 1, num_freq);

for m = 1:num_subj
    ft_rest.powspctrm(m, 1, :) = rest_matrix(:, m);
    ft_run.powspctrm(m, 1, :) = run_matrix(:, m);
end

% Design matrix for paired (within-subject) test
% -------------------------------------------------------------------------
% design = [ivar; uvar]
%   ivar (independent variable): condition label (1 or 2)
%   uvar (unit variable): subject label (1, 2, 3, ...)
%
% Example for 3 subjects:
%   design = [1 1 1 2 2 2;    % condition: 1=REST, 2=RUN
%             1 2 3 1 2 3]    % subject: 1, 2, 3
% -------------------------------------------------------------------------

design = zeros(2, 2 * num_subj);
design(1, :) = [ones(1, num_subj), 2 * ones(1, num_subj)];  % Condition
design(2, :) = [1:num_subj, 1:num_subj];                    % Subject

% Configure cluster-based permutation test
% -------------------------------------------------------------------------
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';      % Dependent (paired) samples T-test
cfg.correctm = 'cluster';           % Cluster-based correction
cfg.clusteralpha = cluster_alpha;   % Cluster-forming threshold (more liberal = more sensitive)
cfg.clusterstatistic = 'maxsum';    % Use maximum sum of t-values
cfg.alpha = cluster_pval;           % Significance threshold

% Tail configuration for one-tailed vs two-tailed tests
% -------------------------------------------------------------------------
%   tail = 0:  Two-tailed test (default) - detects any difference
%   tail = 1:  One-tailed test for RUN > REST (positive t-values)
%              RECOMMENDED when you expect theta coherence to INCREASE during running
%   tail = -1: One-tailed test for REST > RUN (negative t-values)
%
% One-tailed tests have ~2x the statistical power for directional hypotheses!
% -------------------------------------------------------------------------
cfg.tail = tail;
cfg.clustertail = tail;

% Handle 'all' permutations for exact p-values with small N
if ischar(num_rand) && strcmpi(num_rand, 'all')
    cfg.numrandomization = 'all';   % Exhaustive permutation (exact p-values)
else
    cfg.numrandomization = num_rand;
end

cfg.design = design;
cfg.ivar = 1;  % Row 1 is the independent variable (condition)
cfg.uvar = 2;  % Row 2 is the unit variable (subject)

% Frequency selection
cfg.frequency = [freq_min, freq_max];
cfg.avgoverfreq = 'no';

% Run cluster statistics
try
    stat = ft_freqstatistics(cfg, ft_rest, ft_run);
    
    % Extract results into simple structure
    results.success = true;
    results.tvals = squeeze(stat.stat(:));
    results.pvals = squeeze(stat.prob(:));
    results.freq = stat.freq(:);
    
    % Significance mask (ONLY significant clusters, p < cluster_pval)
    results.sig_mask = false(num_freq, 1);
    
    % Add significant positive clusters
    if isfield(stat, 'posclusters') && ~isempty(stat.posclusters) && isfield(stat, 'posclusterslabelmat')
        for c = 1:length(stat.posclusters)
            if stat.posclusters(c).prob < cluster_pval
                cluster_mask = squeeze(stat.posclusterslabelmat(:)) == c;
                results.sig_mask = results.sig_mask | cluster_mask;
            end
        end
    end
    
    % Add significant negative clusters
    if isfield(stat, 'negclusters') && ~isempty(stat.negclusters) && isfield(stat, 'negclusterslabelmat')
        for c = 1:length(stat.negclusters)
            if stat.negclusters(c).prob < cluster_pval
                cluster_mask = squeeze(stat.negclusterslabelmat(:)) == c;
                results.sig_mask = results.sig_mask | cluster_mask;
            end
        end
    end
    
    % Positive clusters
    results.num_pos_clusters = 0;
    results.num_sig_pos = 0;
    results.pos_cluster_pvals = [];
    results.pos_cluster_freqs = {};
    
    if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
        results.num_pos_clusters = length(stat.posclusters);
        results.pos_cluster_pvals = [stat.posclusters.prob];
        results.num_sig_pos = sum(results.pos_cluster_pvals < cluster_pval);
        
        if isfield(stat, 'posclusterslabelmat')
            for c = 1:results.num_pos_clusters
                cluster_mask = squeeze(stat.posclusterslabelmat(:)) == c;
                cluster_freqs = stat.freq(cluster_mask);
                results.pos_cluster_freqs{c} = cluster_freqs;
            end
        end
    end
    
    % Negative clusters
    results.num_neg_clusters = 0;
    results.num_sig_neg = 0;
    results.neg_cluster_pvals = [];
    results.neg_cluster_freqs = {};
    
    if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
        results.num_neg_clusters = length(stat.negclusters);
        results.neg_cluster_pvals = [stat.negclusters.prob];
        results.num_sig_neg = sum(results.neg_cluster_pvals < cluster_pval);
        
        if isfield(stat, 'negclusterslabelmat')
            for c = 1:results.num_neg_clusters
                cluster_mask = squeeze(stat.negclusterslabelmat(:)) == c;
                cluster_freqs = stat.freq(cluster_mask);
                results.neg_cluster_freqs{c} = cluster_freqs;
            end
        end
    end
    
    fprintf('    %s: %d pos clusters (%d sig), %d neg clusters (%d sig)\n', ...
        data_label, results.num_pos_clusters, results.num_sig_pos, ...
        results.num_neg_clusters, results.num_sig_neg);
    
catch ME
    fprintf('    WARNING: Cluster stats for %s failed: %s\n', data_label, ME.message);
    results.success = false;
    results.error = ME.message;
    results.sig_mask = false(num_freq, 1);
    results.tvals = zeros(num_freq, 1);
    results.pvals = ones(num_freq, 1);
    results.freq = freq(:);
end
end

%% ============================================================================
%  HELPER: Ternary operator
%  ============================================================================
function result = ternary(condition, true_val, false_val)
%TERNARY Inline conditional (like C's ?: operator)
if condition
    result = true_val;
else
    result = false_val;
end
end
