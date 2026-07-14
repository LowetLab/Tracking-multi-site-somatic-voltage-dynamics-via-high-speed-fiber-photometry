%% compute_phase_curve.m
%% ============================================================================
%  Phase-Curve Extraction
%  ============================================================================
%
%  DESCRIPTION:
%  Reads the phase-aligned spectrogram results from compute_state_spectrogram.m
%  (run_nonrun folder) and extracts two amplitude-vs-phase curves for each
%  signal type and behavioral state:
%    - Theta curve  : mean wavelet amplitude across THETA_BAND  (left y-axis)
%    - Gamma curve  : mean wavelet amplitude across GAMMA_BAND  (right y-axis)
%
%  Both curves are saved as vectors of length n_phase_bins (1-cycle) and
%  2*n_phase_bins (2-cycle) in a .mat file under Process/phase_curve/.
%
%  REQUIREMENTS:
%  - compute_phase_aligned_spectrogram_batch.m must have been run with
%    HighFreq_Range = [6, 100] so that theta frequencies are in freq_vector.
%  - compute_state_spectrogram.m must have been run to produce run_nonrun data.
%
%  OUTPUT (per signal type x state):
%    Process/phase_curve/{condition}/{animal}/
%      PhaseCurve_{SignalType}_{state}.mat
%        .animal_name       string
%        .condition         string
%        .signal_type       string  ('LFP' | 'Fiber1' | ...)
%        .state             string  ('running' | 'non_running' | 'all_epochs')
%        .theta_band        [1x2]   Hz
%        .gamma_band        [1x2]   Hz
%        .freq_vector       [1xF]   Hz
%        .phase_bins_deg    [1xN]   degrees, 1-cycle (-180 to 180)
%        .two_cycle_deg     [1x2N]  degrees, 2-cycle (0 to 720)
%        .theta_curve_1cyc  [1xN]   theta amplitude, 1 cycle
%        .gamma_curve_1cyc  [1xN]   gamma amplitude, 1 cycle
%        .theta_curve_2cyc  [1x2N]  theta amplitude, 2 cycles
%        .gamma_curve_2cyc  [1x2N]  gamma amplitude, 2 cycles
%        .n_epochs          scalar
%        .mean_velocity     scalar
%        .lfp_ephys_field   string  (LFP only; e.g. lfp_raw_aligned_HP — for plot axis labels)
%
%  ============================================================================

clear; clc; close all;

%% ============================================================================
%  CONFIGURATION
%  ============================================================================

THETA_BAND = [5, 9];     % Hz - low-frequency (carrier) band
GAMMA_BAND = [30, 60];   % Hz - low gamma band

STATES = {'running', 'non_running', 'all_epochs'};

% Velocity thresholds – must match compute_state_spectrogram.m
RUNNING_THRESHOLD     = 2.0;   % cm/s  – epoch running if ALL points > threshold
NON_RUNNING_THRESHOLD = 0.1;   % cm/s  – epoch rest   if ALL points < threshold
MIN_TIME_FRACTION     = 1.0;   % fraction of epoch that must satisfy criterion
EXCLUDE_ARTIFACT      = true;  % skip artifact-contaminated epochs

%% ============================================================================
%  SELECT ANIMAL FOLDER
%  ============================================================================

script_dir   = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);  % Project root = parent of PhaseCurve_Pipeline

default_path = fullfile(project_root, 'Process', 'phase_aligned_spectrogram');
if ~exist(default_path, 'dir')
    default_path = project_root;
end

animal_folder = uigetdir(default_path, 'Select Animal Folder (contains run_nonrun/)');
if isequal(animal_folder, 0)
    disp('Cancelled.'); return;
end

[~, animal_name] = fileparts(animal_folder);
condition_folder = fileparts(animal_folder);
[~, condition_name] = fileparts(condition_folder);

fprintf('\n=== Processing Animal: %s  Condition: %s ===\n', animal_name, condition_name);

run_nonrun_dir = fullfile(animal_folder, 'run_nonrun');
if ~exist(run_nonrun_dir, 'dir')
    error('run_nonrun folder not found: %s\nRun compute_state_spectrogram.m first.', run_nonrun_dir);
end

%% ============================================================================
%  OUTPUT DIRECTORY
%  ============================================================================

output_base = fullfile(project_root, 'Process', 'phase_curve', condition_name, animal_name);
if ~exist(output_base, 'dir')
    mkdir(output_base);
end
fprintf('Output: %s\n', output_base);

%% ============================================================================
%  PRE-COMPUTE PHASE-TRIGGERED WAVEFORMS & 95% CI FROM RAW EPOCHS
%  ============================================================================
%
%  Strategy: compute one mean value per phase bin PER EPOCH, then derive the
%  grand mean and 95% CI across epoch means (Welford-style online variance).
%
%  Theta waveform per sample:  theta_mag(t) × cos(phase(t))
%    = real part of analytic signal = bandpass-filtered 5-9 Hz signal
%
%  Gamma amplitude per sample: mean(spectrogram(gamma_mask, t))
%    = mean wavelet magnitude across GAMMA_BAND at time t
%
%  Both are binned by instantaneous phase and averaged per epoch, then
%  across epochs.  95% CI = 1.96 × SEM (across epoch means).

fprintf('\n--- Pre-computing phase-triggered curves with 95%% CI ---\n');

result_files = dir(fullfile(animal_folder, '**', 'PhaseAlignedSpectrogram_Results_*.mat'));
fprintf('Found %d PhaseAlignedSpectrogram_Results_*.mat files.\n', length(result_files));

% Determine phase bin edges from the first available Results file.
phase_bin_edges_rad = [];
for rfi = 1:length(result_files)
    fpath = fullfile(result_files(rfi).folder, result_files(rfi).name);
    tmp = load(fpath, 'results');
    if isfield(tmp, 'results') && isfield(tmp.results, 'metadata') && ...
            isfield(tmp.results.metadata, 'phase_bin_edges')
        phase_bin_edges_rad = double(tmp.results.metadata.phase_bin_edges(:)');
        break;
    end
end
if isempty(phase_bin_edges_rad)
    phase_bin_edges_rad = linspace(-pi, pi, 37);
    warning('compute_phase_curve:noBinEdges', ...
        'Could not find phase_bin_edges in Results files. Using 36-bin default.');
end
n_ptw_bins = length(phase_bin_edges_rad) - 1;

% Which ephys vector was used for LFP (batch writes results.metadata.lfp_ephys_field)
lfp_ephys_field_global = '';
for rfi = 1:length(result_files)
    fpath = fullfile(result_files(rfi).folder, result_files(rfi).name);
    tmp = load(fpath, 'results');
    if isfield(tmp, 'results') && isfield(tmp.results, 'metadata') && ...
            isfield(tmp.results.metadata, 'lfp_ephys_field')
        lfp_ephys_field_global = char(string(tmp.results.metadata.lfp_ephys_field));
        break;
    end
end

% Pre-allocate epoch-level accumulators (Welford online stats)
%   ep_sum  = sum of per-epoch bin means
%   ep_sq   = sum of squared per-epoch bin means
%   ep_cnt  = number of epochs that contributed data to each bin
SIGNAL_NAMES = {'LFP', 'Fiber1', 'Fiber2'};
for sni = 1:length(SIGNAL_NAMES)
    sn = SIGNAL_NAMES{sni};
    for sti = 1:length(STATES)
        st = STATES{sti};
        ptw_ep_sum.(sn).(st) = zeros(1, n_ptw_bins);
        ptw_ep_sq.(sn).(st)  = zeros(1, n_ptw_bins);
        ptw_ep_cnt.(sn).(st) = zeros(1, n_ptw_bins);
        gam_ep_sum.(sn).(st) = zeros(1, n_ptw_bins);
        gam_ep_sq.(sn).(st)  = zeros(1, n_ptw_bins);
        gam_ep_cnt.(sn).(st) = zeros(1, n_ptw_bins);
    end
end

% ---- 2-cycle accumulators (cycle1 | cycle2, truly independent) ----
n_2cyc = 2 * n_ptw_bins;
for sni = 1:length(SIGNAL_NAMES)
    sn = SIGNAL_NAMES{sni};
    for sti = 1:length(STATES)
        st = STATES{sti};
        ptw_ep_sum_2cyc.(sn).(st) = zeros(1, n_2cyc);
        ptw_ep_sq_2cyc.(sn).(st)  = zeros(1, n_2cyc);
        ptw_ep_cnt_2cyc.(sn).(st) = zeros(1, n_2cyc);
        gam_ep_sum_2cyc.(sn).(st) = zeros(1, n_2cyc);
        gam_ep_sq_2cyc.(sn).(st)  = zeros(1, n_2cyc);
        gam_ep_cnt_2cyc.(sn).(st) = zeros(1, n_2cyc);
    end
end

% Stream through all Results files
for rfi = 1:length(result_files)
    fpath = fullfile(result_files(rfi).folder, result_files(rfi).name);
    fprintf('  [%d/%d] %s\n', rfi, length(result_files), result_files(rfi).name);
    tmp_r = load(fpath, 'results');
    if ~isfield(tmp_r, 'results'); continue; end
    results = tmp_r.results;

    for sni = 1:length(SIGNAL_NAMES)
        sn = SIGNAL_NAMES{sni};
        if ~isfield(results, sn); continue; end
        sig_data = results.(sn);
        if ~isfield(sig_data, 'all_epochs_raw') || ~isfield(sig_data, 'freq_vector')
            continue;
        end

        freq_vec   = sig_data.freq_vector(:)';
        gamma_mask = freq_vec >= GAMMA_BAND(1) & freq_vec <= GAMMA_BAND(2);

        % Get sampling rate and design theta bandpass filter (once per signal per file)
        if isfield(results, 'metadata') && isfield(results.metadata, 'sampling_rate')
            fs_local = results.metadata.sampling_rate;
        else
            fs_local = 1000;  % fallback
            warning('compute_phase_curve:noFs', 'sampling_rate not found; assuming 1000 Hz.');
        end
        [b_theta, a_theta] = butter(4, THETA_BAND / (fs_local/2), 'bandpass');

        epochs = sig_data.all_epochs_raw;
        for ei = 1:length(epochs)
            epoch = epochs{ei};

            % Artifact rejection
            if EXCLUDE_ARTIFACT && isfield(epoch, 'has_artifact') && epoch.has_artifact
                continue;
            end

            % State classification (mirrors compute_state_spectrogram.m; velocity required)
            if ~isfield(epoch, 'velocity') || isempty(epoch.velocity)
                error('compute_phase_curve:MissingVelocity', ...
                    'Epoch %d (%s): missing velocity vector (required for classification).', ei, sn);
            end
            vel = epoch.velocity(:);
            vel = vel(~isnan(vel));
            nv = numel(vel);
            if nv == 0
                error('compute_phase_curve:EmptyVelocity', ...
                    'Epoch %d (%s): velocity is empty or all NaN (required for classification).', ei, sn);
            end
            run_frac    = sum(vel > RUNNING_THRESHOLD)    / nv;
            nonrun_frac = sum(vel < NON_RUNNING_THRESHOLD) / nv;
            is_running    = run_frac    >= MIN_TIME_FRACTION;
            is_non_running = nonrun_frac >= MIN_TIME_FRACTION;

            % ----------------------------------------------------------------
            % TRUE PHASE-TRIGGERED AVERAGE
            %   theta_wave  = bandpass-filtered (5-9 Hz) raw signal
            %               = the actual theta waveform, no cosine assumption
            %   gamma_mag   = mean wavelet magnitude across gamma band
            %               (no equivalent raw-signal version needed for amplitude)
            % ----------------------------------------------------------------
            ph = epoch.phase(:)';    % [1 x T]

            if isfield(epoch, 'theta_filtered') && ~isempty(epoch.theta_filtered)
                % Best path: use pre-filtered signal (filtered on the full trial
                % signal before epoch extraction → no boundary artefacts)
                theta_wave = double(epoch.theta_filtered(:))';   % [1 x T]
            elseif isfield(epoch, 'raw_signal') && ~isempty(epoch.raw_signal)
                % Legacy path: filter the epoch segment (slightly biased at edges)
                raw = double(epoch.raw_signal(:));
                if length(raw) < 3 * 8; continue; end
                theta_wave = filtfilt(b_theta, a_theta, raw)';   % [1 x T]
            else
                % Fallback: analytic reconstruction
                theta_mask_fb = freq_vec >= THETA_BAND(1) & freq_vec <= THETA_BAND(2);
                if ~any(theta_mask_fb); continue; end
                theta_mag  = mean(epoch.spectrogram(theta_mask_fb, :), 1);
                theta_wave = theta_mag .* cos(ph);
            end

            if any(gamma_mask)
                gamma_mag = mean(epoch.spectrogram(gamma_mask, :), 1); % [1 x T]
            else
                gamma_mag = zeros(size(theta_wave));
            end

            % Bin phase
            bin_idx = discretize(ph, phase_bin_edges_rad);
            valid   = ~isnan(bin_idx);
            if ~any(valid); continue; end

            % Per-epoch bin means (NaN if no samples in that bin)
            ep_theta = NaN(1, n_ptw_bins);
            ep_gamma = NaN(1, n_ptw_bins);
            bv = bin_idx(valid);
            for b = 1:n_ptw_bins
                mask_b = (bv == b);
                if any(mask_b)
                    tw = theta_wave(valid);
                    gw = gamma_mag(valid);
                    ep_theta(b) = mean(tw(mask_b));
                    ep_gamma(b) = mean(gw(mask_b));
                end
            end

            % ---- 2-cycle: split epoch at cycle_boundary, bin independently ----
            has_2cyc_ep = false;
            if isfield(epoch, 'cycle_boundary') && ~isempty(epoch.cycle_boundary)
                cb = epoch.cycle_boundary;
                T_ep = length(ph);
                if cb >= 2 && cb <= T_ep - 2
                    ep_theta_2c  = NaN(1, n_2cyc);
                    ep_gamma_2c  = NaN(1, n_2cyc);

                    for cyc_i = 1:2
                        if cyc_i == 1; idx_r = 1:cb;            boff = 0;
                        else;           idx_r = (cb+1):T_ep;    boff = n_ptw_bins;
                        end
                        ph_c  = ph(idx_r);
                        bi_c  = discretize(ph_c, phase_bin_edges_rad);
                        val_c = ~isnan(bi_c);
                        if ~any(val_c); continue; end
                        tw_c  = theta_wave(idx_r);
                        gm_c  = gamma_mag(idx_r);
                        for b = 1:n_ptw_bins
                            mask_b = val_c & (bi_c == b);
                            if any(mask_b)
                                ep_theta_2c(boff+b)  = mean(tw_c(mask_b));
                                ep_gamma_2c(boff+b)  = mean(gm_c(mask_b));
                            end
                        end
                    end
                    has_2cyc_ep = true;
                end
            end

            % Which states does this epoch belong to?
            states_to_update = {'all_epochs'};
            if is_running;     states_to_update{end+1} = 'running';     end
            if is_non_running; states_to_update{end+1} = 'non_running'; end

            % Accumulate epoch-level means into online stats
            for sti = 1:length(states_to_update)
                st = states_to_update{sti};
                for b = 1:n_ptw_bins
                    if ~isnan(ep_theta(b))
                        ptw_ep_sum.(sn).(st)(b) = ptw_ep_sum.(sn).(st)(b) + ep_theta(b);
                        ptw_ep_sq.(sn).(st)(b)  = ptw_ep_sq.(sn).(st)(b)  + ep_theta(b)^2;
                        ptw_ep_cnt.(sn).(st)(b) = ptw_ep_cnt.(sn).(st)(b) + 1;
                    end
                    if ~isnan(ep_gamma(b))
                        gam_ep_sum.(sn).(st)(b) = gam_ep_sum.(sn).(st)(b) + ep_gamma(b);
                        gam_ep_sq.(sn).(st)(b)  = gam_ep_sq.(sn).(st)(b)  + ep_gamma(b)^2;
                        gam_ep_cnt.(sn).(st)(b) = gam_ep_cnt.(sn).(st)(b) + 1;
                    end
                end

                % ---- 2-cycle accumulation ----
                if has_2cyc_ep
                    for b2 = 1:n_2cyc
                        if ~isnan(ep_theta_2c(b2))
                            ptw_ep_sum_2cyc.(sn).(st)(b2) = ptw_ep_sum_2cyc.(sn).(st)(b2) + ep_theta_2c(b2);
                            ptw_ep_sq_2cyc.(sn).(st)(b2)  = ptw_ep_sq_2cyc.(sn).(st)(b2)  + ep_theta_2c(b2)^2;
                            ptw_ep_cnt_2cyc.(sn).(st)(b2) = ptw_ep_cnt_2cyc.(sn).(st)(b2) + 1;
                        end
                        if ~isnan(ep_gamma_2c(b2))
                            gam_ep_sum_2cyc.(sn).(st)(b2) = gam_ep_sum_2cyc.(sn).(st)(b2) + ep_gamma_2c(b2);
                            gam_ep_sq_2cyc.(sn).(st)(b2)  = gam_ep_sq_2cyc.(sn).(st)(b2)  + ep_gamma_2c(b2)^2;
                            gam_ep_cnt_2cyc.(sn).(st)(b2) = gam_ep_cnt_2cyc.(sn).(st)(b2) + 1;
                        end
                    end
                end
            end
        end
    end
end

% Compute grand mean and 95% CI from epoch-level stats
ptw     = struct();   % mean theta waveform
ptw_ci  = struct();   % 95% CI half-width
gam_raw = struct();   % gamma mean from raw epochs
gam_ci  = struct();   % 95% CI for gamma

for sni = 1:length(SIGNAL_NAMES)
    sn = SIGNAL_NAMES{sni};
    for sti = 1:length(STATES)
        st  = STATES{sti};
        cnt_p = ptw_ep_cnt.(sn).(st);
        cnt_g = gam_ep_cnt.(sn).(st);

        % Theta waveform
        % – NaN for bins with no data (prevents false zero-pull)
        % – unbiased sample variance: (sum_sq - sum^2/n) / (n-1)
        m_p = ptw_ep_sum.(sn).(st) ./ cnt_p;          % NaN where cnt=0
        m_p(cnt_p == 0) = NaN;
        v_p_unbiased = (ptw_ep_sq.(sn).(st) - ptw_ep_sum.(sn).(st).^2 ./ max(cnt_p,1)) ...
                        ./ max(cnt_p - 1, 1);
        v_p_unbiased(cnt_p < 2) = NaN;                % need ≥2 epochs for variance
        ptw.(sn).(st)    = m_p;
        ptw_ci.(sn).(st) = 1.96 * sqrt(max(v_p_unbiased, 0) ./ cnt_p);
        ptw_ci.(sn).(st)(cnt_p < 2) = NaN;

        % Gamma amplitude
        m_g = gam_ep_sum.(sn).(st) ./ cnt_g;          % NaN where cnt=0
        m_g(cnt_g == 0) = NaN;
        v_g_unbiased = (gam_ep_sq.(sn).(st) - gam_ep_sum.(sn).(st).^2 ./ max(cnt_g,1)) ...
                        ./ max(cnt_g - 1, 1);
        v_g_unbiased(cnt_g < 2) = NaN;
        gam_raw.(sn).(st) = m_g;
        gam_ci.(sn).(st)  = 1.96 * sqrt(max(v_g_unbiased, 0) ./ cnt_g);
        gam_ci.(sn).(st)(cnt_g < 2) = NaN;
    end
end
% ---- 2-cycle mean + CI (same Welford formula, applied to 2N vectors) ----
ptw_2cyc_mean = struct();  ptw_2cyc_ci = struct();
gam_2cyc_mean = struct();  gam_2cyc_ci = struct();

for sni = 1:length(SIGNAL_NAMES)
    sn = SIGNAL_NAMES{sni};
    for sti = 1:length(STATES)
        st = STATES{sti};
        [ptw_2cyc_mean.(sn).(st), ptw_2cyc_ci.(sn).(st)] = welford_mean_ci( ...
            ptw_ep_sum_2cyc.(sn).(st), ptw_ep_sq_2cyc.(sn).(st), ptw_ep_cnt_2cyc.(sn).(st));
        [gam_2cyc_mean.(sn).(st), gam_2cyc_ci.(sn).(st)] = welford_mean_ci( ...
            gam_ep_sum_2cyc.(sn).(st), gam_ep_sq_2cyc.(sn).(st), gam_ep_cnt_2cyc.(sn).(st));
    end
end

fprintf('Phase-triggered curves with 95%% CI ready.\n');

%% ============================================================================
%  FIND SIGNAL FOLDERS
%  ============================================================================

signal_dirs = dir(run_nonrun_dir);
signal_dirs = signal_dirs([signal_dirs.isdir] & ~startsWith({signal_dirs.name}, '.'));

if isempty(signal_dirs)
    error('No signal folders found in: %s', run_nonrun_dir);
end

fprintf('\nFound signal folders: %s\n', strjoin({signal_dirs.name}, ', '));

%% ============================================================================
%  EXTRACT CURVES FOR EACH SIGNAL x STATE
%  ============================================================================

for si = 1:length(signal_dirs)
    signal_type = signal_dirs(si).name;   % e.g. 'LFP', 'Fiber1'
    signal_dir  = fullfile(run_nonrun_dir, signal_type);

    fprintf('\n--- Signal: %s ---\n', signal_type);

    for st = 1:length(STATES)
        state     = STATES{st};
        state_dir = fullfile(signal_dir, state);
        mat_file  = fullfile(state_dir, sprintf('PhaseAlignedSpectrogram_%s.mat', state));

        if ~exist(mat_file, 'file')
            fprintf('  [%s] Not found, skipping.\n', state);
            continue;
        end

        fprintf('  Loading [%s]... ', state);
        raw = load(mat_file);

        % ---- Extract required fields ----
        if ~isfield(raw, 'freq_vector') || ~isfield(raw, 'mean_spectrogram')
            fprintf('Missing fields, skipping.\n');
            continue;
        end

        freq_vector     = raw.freq_vector(:)';           % [1 x F]
        mean_spec       = raw.mean_spectrogram;          % [F x N_bins]
        n_epochs        = raw.n_epochs;
        mean_velocity   = raw.mean_velocity;

        % Phase bins (degrees, 1-cycle, -180 to 180)
        if isfield(raw, 'phase_bins_deg')
            phase_bins_deg = raw.phase_bins_deg(:)';
        else
            phase_bins_deg = linspace(-180, 180, size(mean_spec, 2) + 1);
            phase_bins_deg = phase_bins_deg(1:end-1) + diff(phase_bins_deg)/2;
        end

        % 2-cycle spectrogram
        if isfield(raw, 'two_cycle_spectrogram')
            spec_2cyc = raw.two_cycle_spectrogram;       % [F x 2N]
        elseif isfield(raw, 'cycle1_mean_spectrogram') && isfield(raw, 'cycle2_mean_spectrogram')
            spec_2cyc = [raw.cycle1_mean_spectrogram, raw.cycle2_mean_spectrogram];
        else
            spec_2cyc = [mean_spec, mean_spec];
        end

        % 2-cycle phase axis (0 to 720°)
        cycle1_deg  = phase_bins_deg + 180;              % shift -180:180 -> 0:360
        two_cycle_deg = [cycle1_deg, cycle1_deg + 360];  % 0 to 720

        % ---- Check frequency coverage ----
        if max(freq_vector) < THETA_BAND(1)
            fprintf('WARN: freq_vector max=%.1f Hz < theta band %.0f Hz.\n', ...
                    max(freq_vector), THETA_BAND(1));
            fprintf('      Re-run MATLAB with HighFreq_Range=[6,100].\n');
            continue;
        end

        % ---- Extract band curves ----
        theta_mask = freq_vector >= THETA_BAND(1) & freq_vector <= THETA_BAND(2);

        if ~any(theta_mask)
            fprintf('WARN: No freq bins in theta band [%.0f %.0f] Hz, skipping.\n', ...
                    THETA_BAND(1), THETA_BAND(2));
            continue;
        end

        theta_curve_1cyc = mean(mean_spec(theta_mask, :), 1);     % [1 x N]
        theta_curve_2cyc = mean(spec_2cyc(theta_mask, :), 1);     % [1 x 2N]

        % ---- Phase-triggered theta waveform + CI (from raw epochs) ----
        n_bins_state = size(mean_spec, 2);  % number of phase bins in state file

        if isfield(ptw, signal_type) && isfield(ptw.(signal_type), state)
            ptw_1cyc    = resample_vec(ptw.(signal_type).(state),    n_bins_state);
            ptw_ci_1cyc = resample_vec(ptw_ci.(signal_type).(state), n_bins_state);
        else
            ptw_1cyc    = zeros(1, n_bins_state);
            ptw_ci_1cyc = zeros(1, n_bins_state);
        end
        % True 2-cycle theta waveform from independent epoch splits
        has_ind_2cyc_ptw = false;
        if isfield(ptw_2cyc_mean, signal_type) && isfield(ptw_2cyc_mean.(signal_type), state)
            p2c = ptw_2cyc_mean.(signal_type).(state);
            if any(~isnan(p2c))
                ptw_2cyc    = resample_vec_2cyc(p2c, n_bins_state);
                ptw_ci_2cyc = resample_vec_2cyc(ptw_2cyc_ci.(signal_type).(state), n_bins_state);
                has_ind_2cyc_ptw = true;
            end
        end
        if ~has_ind_2cyc_ptw
            ptw_2cyc    = [ptw_1cyc,    ptw_1cyc];
            ptw_ci_2cyc = [ptw_ci_1cyc, ptw_ci_1cyc];
        end

        % ---- Gamma curve from spectrogram ----
        gamma_mask = freq_vector >= GAMMA_BAND(1) & freq_vector <= GAMMA_BAND(2);
        if ~any(gamma_mask)
            fprintf('WARN: No freq bins in gamma band [%.0f %.0f] Hz, skipping.\n', ...
                    GAMMA_BAND(1), GAMMA_BAND(2));
            continue;
        end

        gamma_curve_1cyc = mean(mean_spec(gamma_mask, :), 1);     % [1 x N]
        gamma_curve_2cyc = mean(spec_2cyc(gamma_mask, :), 1);     % [1 x 2N]

        fprintf('Done. (%d epochs, theta bins: %d, gamma bins: %d)\n', ...
                n_epochs, sum(theta_mask), sum(gamma_mask));

        % ---- Gamma curve + CI (from raw epochs) ----
        if isfield(gam_raw, signal_type) && isfield(gam_raw.(signal_type), state)
            gam_1cyc    = resample_vec(gam_raw.(signal_type).(state), n_bins_state);
            gam_ci_1cyc = resample_vec(gam_ci.(signal_type).(state),  n_bins_state);
        else
            gam_1cyc    = gamma_curve_1cyc;
            gam_ci_1cyc = zeros(1, n_bins_state);
        end
        % True 2-cycle gamma from independent epoch splits
        has_ind_2cyc_gam = false;
        if isfield(gam_2cyc_mean, signal_type) && isfield(gam_2cyc_mean.(signal_type), state)
            g2c = gam_2cyc_mean.(signal_type).(state);
            if any(~isnan(g2c))
                gam_2cyc    = resample_vec_2cyc(g2c, n_bins_state);
                gam_ci_2cyc = resample_vec_2cyc(gam_2cyc_ci.(signal_type).(state), n_bins_state);
                has_ind_2cyc_gam = true;
            end
        end
        if ~has_ind_2cyc_gam
            gam_2cyc    = [gam_1cyc,    gam_1cyc];
            gam_ci_2cyc = [gam_ci_1cyc, gam_ci_1cyc];
        end

        % ---- Save ----
        out = struct();
        out.animal_name               = animal_name;
        out.condition                 = condition_name;
        out.signal_type               = signal_type;
        out.state                     = state;
        out.theta_band                = THETA_BAND;
        out.gamma_band                = GAMMA_BAND;
        out.freq_vector               = freq_vector;
        out.phase_bins_deg            = phase_bins_deg;
        out.two_cycle_deg             = two_cycle_deg;
        out.theta_curve_1cyc          = theta_curve_1cyc;
        out.gamma_curve_1cyc          = gamma_curve_1cyc;
        out.theta_curve_2cyc          = theta_curve_2cyc;
        out.gamma_curve_2cyc          = gamma_curve_2cyc;
        out.theta_waveform_1cyc       = ptw_1cyc;
        out.theta_waveform_2cyc       = ptw_2cyc;
        out.theta_waveform_ci95_1cyc  = ptw_ci_1cyc;
        out.theta_waveform_ci95_2cyc  = ptw_ci_2cyc;
        out.gamma_raw_1cyc            = gam_1cyc;
        out.gamma_raw_2cyc            = gam_2cyc;
        out.gamma_ci95_1cyc           = gam_ci_1cyc;
        out.gamma_ci95_2cyc           = gam_ci_2cyc;
        out.n_epochs                  = n_epochs;
        out.mean_velocity             = mean_velocity;
        out.has_independent_2cyc      = has_ind_2cyc_ptw && has_ind_2cyc_gam;
        out.lfp_ephys_field           = '';
        if strcmp(signal_type, 'LFP')
            out.lfp_ephys_field = lfp_ephys_field_global;
        end

        out_filename = sprintf('PhaseCurve_%s_%s.mat', signal_type, state);
        save(fullfile(output_base, out_filename), '-struct', 'out', '-v7');
        fprintf('  Saved: %s\n', out_filename);
    end
end

fprintf('\n=== Done. Output: %s ===\n', output_base);

%% ============================================================================
%  LOCAL HELPER
%  ============================================================================

function v = resample_vec(x, n)
% RESAMPLE_VEC  Resample row vector x to length n using linear interpolation.
    x = x(:)';
    if length(x) == n
        v = x;
    else
        v = interp1(1:length(x), x, linspace(1, length(x), n), 'linear');
    end
end

function v = resample_vec_2cyc(x2n, n_target)
% RESAMPLE_VEC_2CYC  Resample a 2-cycle vector [c1|c2] (2N bins) so each
%   half maps to n_target bins, yielding a 2*n_target output.
    x2n = x2n(:)';
    half = floor(length(x2n) / 2);
    c1 = resample_vec(x2n(1:half),       n_target);
    c2 = resample_vec(x2n(half+1:end),   n_target);
    v  = [c1, c2];
end

function [m, ci] = welford_mean_ci(s, sq, cnt)
% WELFORD_MEAN_CI  Mean and 95% CI from Welford online accumulators.
    m = s ./ cnt;
    m(cnt == 0) = NaN;
    v = (sq - s.^2 ./ max(cnt,1)) ./ max(cnt - 1, 1);
    v(cnt < 2) = NaN;
    ci = 1.96 * sqrt(max(v, 0) ./ cnt);
    ci(cnt < 2) = NaN;
end
