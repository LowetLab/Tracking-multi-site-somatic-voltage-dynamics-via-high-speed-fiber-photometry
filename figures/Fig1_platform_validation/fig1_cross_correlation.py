"""
================================================================================
TIME-RESOLVED CROSS-CORRELATION ANALYSIS: GEVI vs LFP
================================================================================

Figure Layout (2x2):
  Top row:    Time-varying cross-correlation heatmaps (Locomotion | Rest)
  Bottom row: Mean cross-correlogram vs temporal offset (ms)

Processing:
  1. Concatenate 6 trials from one session
  2. Classify RUN / REST from speed
  3. Per-bout: bandpass -> trim edges -> sliding window xcorr
     (each window: z-score, Hann taper, normalised correlation)

Cross-correlation: correlate(GEVI, LFP)
  Positive offset  -> LFP lags GEVI
  Negative offset  -> LFP leads GEVI

USAGE:  python fig1_cross_correlation.py
================================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.gridspec import GridSpec
from scipy.signal import butter, filtfilt, correlate
from scipy.io import loadmat
import h5py
from pathlib import Path
import warnings

warnings.filterwarnings('ignore')

# Locate config/paths_config.py by walking up from this file (works regardless
# of how deep this script lives in the repo).
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from common import normalize_unc_path  # shared helper (was a local copy)

# ==============================================================================
#  USER CONFIGURATION -- EDIT THESE FOR YOUR RECORDING
# ==============================================================================

ANIMAL_ID     = 'Animal01'
SESSION_ID    = '01_09_25-R1'
NUM_TRIALS    = 6
FOLDER_SUFFIX = 'fov1_baselineRecording_60sec'

DATA_ROOT  = _LAB_DATA_ROOT / "FiberVoltageImaging"
OUTPUT_DIR = PROJECT_ROOT / "Figures" / "CrossCorrelation_figures"
FIGURE_NAME = f"cross_correlation_{ANIMAL_ID}_{SESSION_ID}"

# ==============================================================================
#  ANALYSIS PARAMETERS
# ==============================================================================

DEFAULT_FS = 500

RUN_THRESHOLD_CMS  = 2.0
REST_THRESHOLD_CMS = 0.1
MIN_BOUT_SEC       = 0.5

THETA_BAND = (5, 9)
DELTA_BAND = (0.5, 4)

WINDOW_SEC   = 2.0
STEP_SEC     = 0.25
MAX_LAG_MS   = 200
EDGE_TRIM_SEC = 0.5

# ==============================================================================
#  MULTI-SESSION SURROGATE TESTING
# ==============================================================================

ALL_SESSIONS = [
    {'session_id': '01_09_25-R1', 'num_trials': 6},
    {'session_id': '02_09_25-R1', 'num_trials': 6},
    {'session_id': '03_09_25-R1', 'num_trials': 6},
    {'session_id': '04_09_25-R1', 'num_trials': 2},
    {'session_id': '04_09_25-R2', 'num_trials': 2},
]

SURROGATE_BAND   = THETA_BAND
N_SURROGATES     = 500
MIN_SHIFT_CYCLES = 2
USE_FISHER_Z     = True

# ==============================================================================
#  FIGURE STYLING
# ==============================================================================

FIGURE_DPI     = 300
FIGURE_FORMATS = ['png', 'pdf']

FONT_SIZE_SUPTITLE   = 18
FONT_SIZE_TITLE      = 16
FONT_SIZE_LABEL      = 14
FONT_SIZE_TICK       = 12
FONT_SIZE_CBAR       = 11
FONT_SIZE_ANNOTATION = 10
FONT_SIZE_STATS      = 12

AXIS_LINEWIDTH = 2.0
TICK_WIDTH     = 1.8
TICK_LENGTH    = 6
LINE_WIDTH     = 2.5

COLOR_XCORR  = 'black'
CMAP_HEATMAP = 'inferno'

# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

def find_trial_path(animal_id, session_id, trial_num, folder_suffix=None):
    session_path = DATA_ROOT / animal_id / "Fiber_Voltage_Processed" / session_id
    if folder_suffix:
        trial_patterns = [
            f"Trial{trial_num}_{folder_suffix}*",
            f"Trial{trial_num}_{folder_suffix}",
        ]
    else:
        trial_patterns = [f"Trial{trial_num}_*", f"Trial{trial_num}", f"*Trial{trial_num}*"]

    for pattern in trial_patterns:
        trial_folders = list(session_path.glob(pattern))
        if trial_folders:
            trial_folder = trial_folders[0]
            mat_patterns = [
                f"{animal_id}-{session_id}_Trial{trial_num}_FiberPhotometry_Analysis.mat",
                f"{animal_id}*Trial{trial_num}*FiberPhotometry*.mat",
                "*FiberPhotometry_Analysis.mat",
            ]
            for mp in mat_patterns:
                mat_files = list(trial_folder.glob(mp))
                if mat_files:
                    return mat_files[0]
    return None


def load_trial_data(mat_path):
    mat_path_access = normalize_unc_path(str(mat_path), for_access=True)
    if not os.path.exists(mat_path_access):
        return None

    data = {}

    try:
        with h5py.File(mat_path_access, 'r') as f:
            root = f['FiberPhotometryAnalysis']
            if 'signals' in root:
                sig = root['signals']
                if 'deltaF_F_traces' in sig:
                    gevi_raw = np.array(sig['deltaF_F_traces'][()])
                    if gevi_raw.ndim == 2:
                        gevi = gevi_raw.T if gevi_raw.shape[0] < gevi_raw.shape[1] else gevi_raw
                    else:
                        gevi = gevi_raw.reshape(-1, 1)
                    data['gevi'] = gevi[:, 0]
            if 'ephys' in root:
                ephys = root['ephys']
                if 'lfp_raw_aligned_HP' in ephys:
                    data['lfp'] = np.array(ephys['lfp_raw_aligned_HP'][()]).flatten()
                elif 'lfp_z_HP' in ephys:
                    data['lfp'] = np.array(ephys['lfp_z_HP'][()]).flatten()
                if 'running_velocity_smooth' in ephys:
                    data['speed'] = np.array(ephys['running_velocity_smooth'][()]).flatten()
            if 'parameters' in root and 'sampling_rate' in root['parameters']:
                data['fs'] = float(np.array(root['parameters']['sampling_rate'][()]).flatten()[0])
            if 'fs' not in data:
                data['fs'] = DEFAULT_FS
            return data
    except Exception:
        pass

    try:
        mat = loadmat(mat_path_access, squeeze_me=True, struct_as_record=False)
        if 'FiberPhotometryAnalysis' in mat:
            fpa = mat['FiberPhotometryAnalysis']
            if hasattr(fpa, 'signals') and hasattr(fpa.signals, 'deltaF_F_traces'):
                gevi = np.atleast_2d(fpa.signals.deltaF_F_traces)
                if gevi.shape[0] < gevi.shape[1]:
                    gevi = gevi.T
                data['gevi'] = gevi[:, 0]
            if hasattr(fpa, 'ephys'):
                if hasattr(fpa.ephys, 'lfp_raw_aligned_HP'):
                    data['lfp'] = np.array(fpa.ephys.lfp_raw_aligned_HP).flatten()
                if hasattr(fpa.ephys, 'running_velocity_smooth'):
                    data['speed'] = np.array(fpa.ephys.running_velocity_smooth).flatten()
            if hasattr(fpa, 'parameters') and hasattr(fpa.parameters, 'sampling_rate'):
                data['fs'] = float(fpa.parameters.sampling_rate)
            else:
                data['fs'] = DEFAULT_FS
            return data
    except Exception:
        pass

    return None


def bandpass_filter(sig, fs, low_freq, high_freq, order=4):
    sig = np.asarray(sig, dtype=np.float64)
    sig = sig - np.nanmean(sig)
    sig = np.nan_to_num(sig, nan=0.0)
    nyq = fs / 2
    b, a = butter(order, [max(low_freq, 0.01) / nyq, min(high_freq, nyq * 0.95) / nyq], btype='band')
    return filtfilt(b, a, sig)


def classify_behavior(speed, fs, run_thresh, rest_thresh, min_bout_sec):
    """Simple behaviour classification — no boundary logic."""
    min_samples = int(min_bout_sec * fs)
    is_run_raw  = speed > run_thresh
    is_rest_raw = speed < rest_thresh

    def find_bouts(mask):
        bouts, in_bout, start = [], False, 0
        for i in range(len(mask)):
            if mask[i] and not in_bout:
                in_bout, start = True, i
            elif not mask[i] and in_bout:
                in_bout = False
                if i - start >= min_samples:
                    bouts.append((start, i))
        if in_bout and len(mask) - start >= min_samples:
            bouts.append((start, len(mask)))
        return bouts

    run_bouts  = find_bouts(is_run_raw)
    rest_bouts = find_bouts(is_rest_raw)

    is_run  = np.zeros(len(speed), dtype=bool)
    is_rest = np.zeros(len(speed), dtype=bool)
    for s, e in run_bouts:
        is_run[s:e] = True
    for s, e in rest_bouts:
        is_rest[s:e] = True

    return is_run, is_rest, run_bouts, rest_bouts


def compute_bout_xcorr(lfp_raw, gevi_raw, fs, bouts, freq_band,
                       window_sec, step_sec, max_lag_ms, edge_trim_sec):
    """Per-bout bandpass -> trim -> sliding-window normalised xcorr."""
    window_samples    = int(window_sec * fs)
    step_samples      = int(step_sec * fs)
    max_lag_samples   = int(max_lag_ms * fs / 1000)
    edge_trim_samples = int(edge_trim_sec * fs)
    min_bout_samples  = window_samples + 2 * edge_trim_samples
    hann = np.hanning(window_samples)

    all_xcorr, all_centers = [], []
    bouts_used = 0

    for bout_start, bout_end in bouts:
        if bout_end - bout_start < min_bout_samples:
            continue

        lfp_filt  = bandpass_filter(lfp_raw[bout_start:bout_end], fs, freq_band[0], freq_band[1])
        gevi_filt = bandpass_filter(gevi_raw[bout_start:bout_end], fs, freq_band[0], freq_band[1])

        lfp_trimmed  = lfp_filt[edge_trim_samples:-edge_trim_samples]
        gevi_trimmed = gevi_filt[edge_trim_samples:-edge_trim_samples]
        trimmed_start = bout_start + edge_trim_samples

        pos = 0
        while pos + window_samples <= len(lfp_trimmed):
            lw = lfp_trimmed[pos:pos + window_samples].copy()
            gw = gevi_trimmed[pos:pos + window_samples].copy()

            ls, gs = np.std(lw), np.std(gw)
            if ls == 0 or gs == 0:
                pos += step_samples
                continue
            lw = (lw - np.mean(lw)) / ls
            gw = (gw - np.mean(gw)) / gs

            lw *= hann
            gw *= hann

            xc = correlate(gw, lw, mode='full')
            norm = np.sqrt(np.sum(lw**2) * np.sum(gw**2))
            if norm > 0:
                xc /= norm

            c = len(xc) // 2
            all_xcorr.append(xc[c - max_lag_samples:c + max_lag_samples + 1])
            all_centers.append((trimmed_start + pos + window_samples / 2) / fs)
            pos += step_samples

        bouts_used += 1

    if not all_xcorr:
        return None, None, None

    xcorr_matrix   = np.column_stack(all_xcorr)
    window_centers = np.array(all_centers)
    lags_ms = np.arange(-max_lag_samples, max_lag_samples + 1) * 1000 / fs

    print(f"      {bouts_used}/{len(bouts)} bouts, {xcorr_matrix.shape[1]} windows")
    return xcorr_matrix, lags_ms, window_centers


# ==============================================================================
#  SURROGATE-TESTED CROSS-CORRELATION FUNCTIONS
# ==============================================================================

def mean_xcorr_windows(mat):
    """Mean cross-correlogram across sliding windows (axis=1), Fisher-z averaged."""
    if mat is None or mat.size == 0:
        return None
    if USE_FISHER_Z:
        clipped = np.clip(mat, -0.999999, 0.999999)
        return np.tanh(np.nanmean(np.arctanh(clipped), axis=1))
    return np.nanmean(mat, axis=1)


def get_peak_lag_limit(band):
    """Half-period at band centre frequency (ms)."""
    return 1000.0 / (2.0 * 0.5 * (band[0] + band[1]))


def peak_index_restricted(mean_xcorr_arr, lags_ms, band):
    """Index of peak |r| within the physiological lag window (± half theta period)."""
    limit = get_peak_lag_limit(band)
    mask = np.abs(lags_ms) <= limit
    cands = np.where(mask)[0]
    if len(cands) == 0:
        cands = np.arange(len(mean_xcorr_arr))
    return cands[np.argmax(np.abs(mean_xcorr_arr[cands]))], limit


def compute_trial_xcorr(lfp, gevi, fs, band,
                        window_sec=WINDOW_SEC, step_sec=STEP_SEC,
                        max_lag_ms=MAX_LAG_MS, edge_trim_sec=EDGE_TRIM_SEC):
    """Sliding-window normalised cross-correlation for one trial (no behaviour split).

    Returns dict with mat (n_lags x n_windows), lags_ms, centers, n_windows,
    or None if the trial is too short.
    """
    lfp_filt = bandpass_filter(lfp, fs, band[0], band[1])
    gevi_filt = bandpass_filter(gevi, fs, band[0], band[1])
    trim = int(edge_trim_sec * fs)
    if len(lfp_filt) <= 2 * trim + 10:
        return None
    lfp_filt = lfp_filt[trim:-trim]
    gevi_filt = gevi_filt[trim:-trim]
    n = len(lfp_filt)

    win = int(window_sec * fs)
    step = max(1, int(step_sec * fs))
    ml = int(max_lag_ms / 1000 * fs)
    hann = np.hanning(win)

    cols, ctrs = [], []
    for s0 in range(0, n - win + 1, step):
        lw = lfp_filt[s0:s0 + win].copy()
        gw = gevi_filt[s0:s0 + win].copy()
        lw -= np.mean(lw)
        gw -= np.mean(gw)
        ls, gs = np.std(lw), np.std(gw)
        if ls < 1e-10 or gs < 1e-10:
            continue
        lw /= ls
        gw /= gs
        lw *= hann
        gw *= hann
        xc = correlate(gw, lw, mode='full')
        norm = np.sqrt(np.sum(lw**2) * np.sum(gw**2))
        if norm > 0:
            xc /= norm
        mid = len(xc) // 2
        cols.append(xc[mid - ml:mid + ml + 1])
        ctrs.append(edge_trim_sec + (s0 + win / 2) / fs)

    if not cols:
        return None
    return {
        'mat': np.column_stack(cols),
        'lags_ms': np.arange(-ml, ml + 1) * 1000 / fs,
        'centers': np.array(ctrs),
        'n_windows': len(cols),
    }


def circular_shift_surrogate(lfp, gevi, fs, band, n_surrogates=N_SURROGATES):
    """Null distribution for peak |r| via circular time-shift of LFP.

    For each surrogate: circularly shift the filtered LFP by a random amount
    (>= MIN_SHIFT_CYCLES theta cycles), recompute sliding-window xcorr,
    extract peak |r| within the physiological lag window.
    """
    rng = np.random.default_rng(42)
    lfp_filt = bandpass_filter(lfp, fs, band[0], band[1])
    gevi_filt = bandpass_filter(gevi, fs, band[0], band[1])

    trim = int(EDGE_TRIM_SEC * fs)
    if len(lfp_filt) <= 2 * trim + 10:
        return np.array([])
    lfp_filt = lfp_filt[trim:-trim]
    gevi_filt = gevi_filt[trim:-trim]
    n = len(lfp_filt)

    f_mid = 0.5 * (band[0] + band[1])
    min_sh = int(MIN_SHIFT_CYCLES / f_mid * fs)
    max_sh = n - min_sh
    if max_sh <= min_sh:
        return np.array([])

    win = int(WINDOW_SEC * fs)
    step = max(1, int(STEP_SEC * fs))
    ml = int(MAX_LAG_MS / 1000 * fs)
    hann = np.hanning(win)
    starts = np.arange(0, n - win + 1, step)
    lags_arr = np.arange(-ml, ml + 1) * 1000 / fs
    limit = get_peak_lag_limit(band)

    null_peaks = np.full(n_surrogates, np.nan)
    for si in range(n_surrogates):
        shift = rng.integers(min_sh, max_sh)
        lfp_sh = np.roll(lfp_filt, shift)
        r_cols = []
        for s0 in starts:
            lw = lfp_sh[s0:s0 + win].copy()
            gw = gevi_filt[s0:s0 + win].copy()
            lw -= np.mean(lw)
            gw -= np.mean(gw)
            ls, gs = np.std(lw), np.std(gw)
            if ls < 1e-10 or gs < 1e-10:
                continue
            lw /= ls
            gw /= gs
            lw *= hann
            gw *= hann
            xc = correlate(gw, lw, mode='full')
            nv = np.sqrt(np.sum(lw**2) * np.sum(gw**2))
            if nv > 0:
                xc /= nv
            mid = len(xc) // 2
            r_cols.append(xc[mid - ml:mid + ml + 1])
        if not r_cols:
            continue
        mean_r = mean_xcorr_windows(np.column_stack(r_cols))
        if mean_r is None:
            continue
        mask = np.abs(lags_arr) <= limit
        sub = np.where(mask)[0]
        null_peaks[si] = float(np.nanmax(np.abs(mean_r[sub]))) if sub.size > 0 \
            else float(np.nanmax(np.abs(mean_r)))
    return null_peaks[np.isfinite(null_peaks)]


def collect_all_trials_xcorr(animal_id, sessions, folder_suffix, band):
    """Per-trial xcorr + circular-shift surrogates across all sessions.

    Returns dict with per-trial results, stacked matrix, grand mean +/- SEM,
    pooled null distribution, pooled p-value, and example trial data.
    """
    trials_out, all_null = [], []
    lags_ref = None
    failed = 0
    example_out, example_label = None, ""

    for sess in sessions:
        sid, nt = sess['session_id'], sess['num_trials']
        for t in range(1, nt + 1):
            mp = find_trial_path(animal_id, sid, t, folder_suffix)
            if mp is None:
                print(f"    {sid} T{t}: NOT FOUND")
                failed += 1
                continue
            d = load_trial_data(mp)
            if d is None or 'lfp' not in d or 'gevi' not in d:
                print(f"    {sid} T{t}: LOAD ERROR")
                failed += 1
                continue

            fs_val = d.get('fs', DEFAULT_FS)
            lfp_t = d['lfp'].astype(np.float64)
            gevi_t = d['gevi'].astype(np.float64)
            n = min(len(lfp_t), len(gevi_t))
            lfp_t, gevi_t = lfp_t[:n], gevi_t[:n]

            out = compute_trial_xcorr(lfp_t, gevi_t, fs_val, band)
            if out is None:
                print(f"    {sid} T{t}: TOO SHORT")
                failed += 1
                continue
            if example_out is None:
                example_out = out
                example_label = f"{sid} T{t}"

            lags_ms = out['lags_ms']
            if lags_ref is None:
                lags_ref = lags_ms

            mean_all = mean_xcorr_windows(out['mat'])
            if mean_all is None:
                failed += 1
                continue

            ir, _ = peak_index_restricted(mean_all, lags_ms, band)
            peak_r = float(np.abs(mean_all[ir]))
            peak_lag = float(lags_ms[ir])

            print(f"    {sid} T{t}: surrogates ({N_SURROGATES})...",
                  end="", flush=True)
            null_peaks = circular_shift_surrogate(lfp_t, gevi_t, fs_val, band)
            p_surr = np.nan
            if len(null_peaks) > 0:
                p_surr = float(np.mean(null_peaks >= peak_r))
                all_null.append(null_peaks)
            print(f" peak|r|={peak_r:.3f}, lag={peak_lag:.1f}ms, p={p_surr:.4f}")

            trials_out.append({
                'session': sid, 'trial_num': t,
                'mean_all': mean_all, 'peak_r': peak_r, 'peak_lag': peak_lag,
                'n_windows': out['n_windows'], 'p_surrogate': p_surr,
            })

    if not trials_out or lags_ref is None:
        return None

    mat_all = np.vstack([t['mean_all'] for t in trials_out])
    n_tr = len(trials_out)
    if USE_FISHER_Z:
        clipped = np.clip(mat_all, -0.999999, 0.999999)
        grand_mean = np.tanh(np.nanmean(np.arctanh(clipped), axis=0))
    else:
        grand_mean = np.nanmean(mat_all, axis=0)
    sem = np.nanstd(mat_all, axis=0, ddof=1) / np.sqrt(n_tr) if n_tr > 1 else None

    obs_peaks = np.array([t['peak_r'] for t in trials_out])
    pooled_null = np.concatenate(all_null) if all_null else np.array([])
    pooled_p = np.nan
    if len(pooled_null) > 0 and len(obs_peaks) > 0:
        median_obs = float(np.median(obs_peaks))
        pooled_p = float(np.mean(pooled_null >= median_obs))
    n_sig = sum(1 for t in trials_out if t['p_surrogate'] < 0.05)

    print(f"\n  SUMMARY: {n_tr} trials loaded, {failed} failed, "
          f"{n_sig}/{n_tr} significant (p<0.05)")

    return {
        'lags_ms': lags_ref, 'trials': trials_out,
        'mat_all': mat_all, 'grand_mean': grand_mean, 'sem': sem,
        'n_loaded': n_tr, 'n_failed': failed,
        'observed_peaks': obs_peaks, 'pooled_null': pooled_null,
        'pooled_p': pooled_p, 'n_sig_trials': n_sig,
        'example_out': example_out, 'example_label': example_label,
    }


# ==============================================================================
#  PLOTTING
# ==============================================================================

def style_axis(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis='both', which='major',
                   width=TICK_WIDTH, length=TICK_LENGTH, labelsize=FONT_SIZE_TICK)


def plot_heatmap(ax, xcorr_matrix, window_centers, lags_ms, cmap, vmin, vmax,
                 title, freq_label):
    """Straightforward imshow heatmap — no grid projection, no NaN tricks."""
    if xcorr_matrix is None or window_centers is None or len(window_centers) == 0:
        ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
        return

    extent = [window_centers[0], window_centers[-1], lags_ms[0], lags_ms[-1]]
    ax.imshow(xcorr_matrix, aspect='auto', origin='lower', extent=extent,
              cmap=cmap, vmin=vmin, vmax=vmax,
              interpolation='bilinear', rasterized=True)

    ax.axhline(y=0, color='white', linestyle='--', linewidth=1.5, alpha=0.8)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    sm.set_array([])
    cbar = ax.figure.colorbar(sm, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label('Cross-correlation\ncoefficient', fontsize=FONT_SIZE_CBAR)
    cbar.ax.tick_params(labelsize=FONT_SIZE_TICK - 1)

    ax.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
    ax.set_title(f'{title} ({freq_label})', fontsize=FONT_SIZE_TITLE, fontweight='bold', loc='left')
    ax.tick_params(axis='both', labelsize=FONT_SIZE_TICK, width=TICK_WIDTH, length=TICK_LENGTH)


def plot_mean_xcorr(ax, lags_ms, mean_xcorr, condition_label):
    """Plot the raw mean cross-correlogram with peak stats."""
    peak_offset, peak_val = None, None

    if mean_xcorr is not None:
        ax.plot(lags_ms, mean_xcorr, '-', color=COLOR_XCORR, linewidth=LINE_WIDTH)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)

        peak_idx    = np.argmax(np.abs(mean_xcorr))
        peak_offset = lags_ms[peak_idx]
        peak_val    = mean_xcorr[peak_idx]

        print(f"   {condition_label} peak: r={peak_val:.3f} at offset={peak_offset:.1f} ms")

        ax.text(0.98, 0.95, f'Peak r = {peak_val:.3f}\nPeak offset = {peak_offset:.0f} ms',
                transform=ax.transAxes, fontsize=FONT_SIZE_STATS, ha='right', va='top')

    ax.set_xlabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Cross-correlation coefficient', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim(-MAX_LAG_MS, MAX_LAG_MS)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    style_axis(ax)
    return peak_offset, peak_val


def plot_surrogate_test_panel(ax, all_trials_data):
    """Null distribution histogram with observed peak |r| markers and pooled p-value."""
    if all_trials_data is None:
        ax.text(0.5, 0.5, 'No data', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
        ax.set_title('Surrogate test', fontsize=FONT_SIZE_TITLE,
                     fontweight='bold', pad=8)
        style_axis(ax)
        return
    null_dist = all_trials_data['pooled_null']
    obs_peaks = all_trials_data['observed_peaks']
    pooled_p = all_trials_data['pooled_p']
    n_sig = all_trials_data['n_sig_trials']
    n_total = all_trials_data['n_loaded']
    if len(null_dist) == 0:
        ax.text(0.5, 0.5, 'No surrogates', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
        ax.set_title('Surrogate test', fontsize=FONT_SIZE_TITLE,
                     fontweight='bold', pad=8)
        style_axis(ax)
        return
    ax.hist(null_dist, bins=40, density=True, color='0.70',
            edgecolor='0.55', alpha=0.7, label='Null (circular shift)', zorder=1)
    if len(obs_peaks) > 0:
        for pk in obs_peaks:
            ax.axvline(pk, color='firebrick', ls='-', lw=1.0, alpha=0.4, zorder=2)
        median_obs = float(np.median(obs_peaks))
        ax.axvline(median_obs, color='firebrick', ls='-', lw=3.0, zorder=3,
                   label=f'Observed median |r| = {median_obs:.3f}')
        p95 = float(np.percentile(null_dist, 95))
        ax.axvline(p95, color='0.35', ls='--', lw=1.5, zorder=2,
                   label=f'Null 95th pctl = {p95:.3f}')
    ax.set_xlabel('Peak |r|', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Density', fontsize=FONT_SIZE_LABEL)
    ax.set_title(f'Circular-shift surrogate test (n={N_SURROGATES}/trial)',
                 fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=8)
    ax.legend(loc='upper right', fontsize=FONT_SIZE_ANNOTATION, frameon=True,
              framealpha=0.95, edgecolor='0.5')
    style_axis(ax)
    p_str = (f'p = {pooled_p:.4f}' if np.isfinite(pooled_p) and pooled_p >= 0.001
             else f'p = {pooled_p:.2e}' if np.isfinite(pooled_p) else 'p = n/a')
    stats_txt = (f'Trials significant: {n_sig}/{n_total}\n'
                 f'Pooled median vs null: {p_str}\n'
                 f'N surrogates/trial: {N_SURROGATES}')
    ax.text(0.02, 0.98, stats_txt, transform=ax.transAxes,
            fontsize=FONT_SIZE_ANNOTATION, va='top', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor='0.55', alpha=0.92))


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    print("=" * 70)
    print("TIME-RESOLVED CROSS-CORRELATION ANALYSIS")
    print(f"Animal: {ANIMAL_ID},  Session: {SESSION_ID}  ({NUM_TRIALS} trials)")
    print("=" * 70)

    # --- Load trials ---
    print(f"\n1. Loading {NUM_TRIALS} trials from {SESSION_ID}...")
    all_lfp, all_gevi, all_speed, all_time = [], [], [], []
    fs = None
    time_offset = 0
    total_loaded = 0

    for t in range(1, NUM_TRIALS + 1):
        mat_path = find_trial_path(ANIMAL_ID, SESSION_ID, t, FOLDER_SUFFIX)
        if mat_path is None:
            print(f"   Trial {t}: NOT FOUND"); continue
        data = load_trial_data(mat_path)
        if data is None or 'lfp' not in data or 'gevi' not in data:
            print(f"   Trial {t}: ERROR"); continue
        if fs is None:
            fs = data.get('fs', DEFAULT_FS)

        n = min(len(data['lfp']), len(data['gevi']))
        if 'speed' in data:
            n = min(n, len(data['speed']))
            all_speed.append(data['speed'][:n])
        all_lfp.append(data['lfp'][:n].astype(np.float64))
        all_gevi.append(data['gevi'][:n].astype(np.float64))
        all_time.append(np.arange(n) / fs + time_offset)
        time_offset = all_time[-1][-1] + 1 / fs
        total_loaded += 1
        print(f"   Trial {t}: {n/fs:.1f}s")

    if not all_lfp:
        print("ERROR: No valid trials"); return

    lfp   = np.concatenate(all_lfp)
    gevi  = np.concatenate(all_gevi)
    speed = np.concatenate(all_speed) if all_speed else np.zeros(len(lfp))
    print(f"\n   Concatenated: {total_loaded} trials, {len(lfp)/fs:.1f}s")

    # --- Classify behaviour ---
    print("\n2. Classifying behaviour...")
    _, _, run_bouts, rest_bouts = classify_behavior(
        speed, fs, RUN_THRESHOLD_CMS, REST_THRESHOLD_CMS, MIN_BOUT_SEC)
    print(f"   RUN:  {len(run_bouts)} bouts")
    print(f"   REST: {len(rest_bouts)} bouts")

    # --- Cross-correlation ---
    print(f"\n3. Computing GEVI-LFP cross-correlation...")
    print(f"   Window {WINDOW_SEC}s, step {STEP_SEC}s, max lag ±{MAX_LAG_MS}ms")

    print(f"\n   RUN (theta {THETA_BAND})...")
    xcorr_run, lags_ms, wc_run = compute_bout_xcorr(
        lfp, gevi, fs, run_bouts, THETA_BAND,
        WINDOW_SEC, STEP_SEC, MAX_LAG_MS, EDGE_TRIM_SEC)
    mean_run = np.mean(xcorr_run, axis=1) if xcorr_run is not None else None

    print(f"\n   REST (delta {DELTA_BAND})...")
    xcorr_rest, lags_rest, wc_rest = compute_bout_xcorr(
        lfp, gevi, fs, rest_bouts, DELTA_BAND,
        WINDOW_SEC, STEP_SEC, MAX_LAG_MS, EDGE_TRIM_SEC)
    mean_rest = np.mean(xcorr_rest, axis=1) if xcorr_rest is not None else None

    # --- Figure ---
    print("\n4. Creating figure...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.subplots_adjust(left=0.10, right=0.95, top=0.90, bottom=0.10,
                        wspace=0.30, hspace=0.40)

    vals = []
    if xcorr_run is not None:
        vals.extend(xcorr_run.flatten())
    if xcorr_rest is not None:
        vals.extend(xcorr_rest.flatten())
    if vals:
        vmax = min(max(np.ceil(np.percentile(np.abs(vals), 98) * 10) / 10, 0.2), 0.5)
    else:
        vmax = 0.3
    vmin = -vmax

    plot_heatmap(axes[0, 0], xcorr_run, wc_run, lags_ms, CMAP_HEATMAP, vmin, vmax,
                 'Locomotion', f'Theta {THETA_BAND[0]}-{THETA_BAND[1]} Hz')
    plot_heatmap(axes[0, 1], xcorr_rest, wc_rest, lags_rest, CMAP_HEATMAP, vmin, vmax,
                 'Rest', f'Delta {DELTA_BAND[0]}-{DELTA_BAND[1]} Hz')

    plot_mean_xcorr(axes[1, 0], lags_ms, mean_run, 'RUN')
    plot_mean_xcorr(axes[1, 1], lags_rest, mean_rest, 'REST')

    fig.suptitle(f'GEVI-LFP Cross-Correlation: {ANIMAL_ID} {SESSION_ID} '
                 f'({total_loaded} trials)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.97)

    # --- Save behaviour-separated figure ---
    print("\n5. Saving behaviour-separated figure...")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for fmt in FIGURE_FORMATS:
        out = OUTPUT_DIR / f"{FIGURE_NAME}.{fmt}"
        fig.savefig(str(out), format=fmt, dpi=FIGURE_DPI,
                    bbox_inches='tight', facecolor='white')
        print(f"   {out}")
    plt.close(fig)

    # ==========================================================================
    #  FIGURE 2: ALL DATA (no behaviour separation)
    # ==========================================================================
    print("\n6. Computing xcorr on ALL data (no RUN/REST split)...")
    all_bout = [(0, len(lfp))]

    print(f"\n   Theta {THETA_BAND}...")
    xcorr_theta, lags_theta, wc_theta = compute_bout_xcorr(
        lfp, gevi, fs, all_bout, THETA_BAND,
        WINDOW_SEC, STEP_SEC, MAX_LAG_MS, EDGE_TRIM_SEC)
    mean_theta = np.mean(xcorr_theta, axis=1) if xcorr_theta is not None else None

    print(f"\n   Delta {DELTA_BAND}...")
    xcorr_delta, lags_delta, wc_delta = compute_bout_xcorr(
        lfp, gevi, fs, all_bout, DELTA_BAND,
        WINDOW_SEC, STEP_SEC, MAX_LAG_MS, EDGE_TRIM_SEC)
    mean_delta = np.mean(xcorr_delta, axis=1) if xcorr_delta is not None else None

    print("\n7. Creating all-data figure...")
    fig2, ax2 = plt.subplots(2, 2, figsize=(14, 10))
    fig2.subplots_adjust(left=0.10, right=0.95, top=0.90, bottom=0.10,
                         wspace=0.30, hspace=0.40)

    vals2 = []
    if xcorr_theta is not None:
        vals2.extend(xcorr_theta.flatten())
    if xcorr_delta is not None:
        vals2.extend(xcorr_delta.flatten())
    if vals2:
        vmax2 = min(max(np.ceil(np.percentile(np.abs(vals2), 98) * 10) / 10, 0.2), 0.5)
    else:
        vmax2 = 0.3
    vmin2 = -vmax2

    plot_heatmap(ax2[0, 0], xcorr_theta, wc_theta, lags_theta, CMAP_HEATMAP,
                 vmin2, vmax2, 'Theta band', f'{THETA_BAND[0]}-{THETA_BAND[1]} Hz')
    plot_heatmap(ax2[0, 1], xcorr_delta, wc_delta, lags_delta, CMAP_HEATMAP,
                 vmin2, vmax2, 'Delta band', f'{DELTA_BAND[0]}-{DELTA_BAND[1]} Hz')

    plot_mean_xcorr(ax2[1, 0], lags_theta, mean_theta, 'THETA')
    plot_mean_xcorr(ax2[1, 1], lags_delta, mean_delta, 'DELTA')

    fig2.suptitle(f'GEVI-LFP Cross-Correlation: {ANIMAL_ID} {SESSION_ID} '
                  f'({total_loaded} trials, all data)',
                  fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.97)

    fig2_name = f"{FIGURE_NAME}_alldata"
    print("\n8. Saving all-data figure...")
    for fmt in FIGURE_FORMATS:
        out = OUTPUT_DIR / f"{fig2_name}.{fmt}"
        fig2.savefig(str(out), format=fmt, dpi=FIGURE_DPI,
                     bbox_inches='tight', facecolor='white')
        print(f"   {out}")
    plt.close(fig2)

    # ==========================================================================
    #  FIGURE 3: SURROGATE-TESTED GEVI-LFP CROSS-CORRELATION (ALL SESSIONS)
    # ==========================================================================
    print("\n" + "=" * 70)
    print("SURROGATE-TESTED CROSS-CORRELATION")
    print(f"Animal: {ANIMAL_ID}, Sessions: {len(ALL_SESSIONS)}, "
          f"Band: {SURROGATE_BAND[0]}-{SURROGATE_BAND[1]} Hz, "
          f"Surrogates: {N_SURROGATES}/trial")
    print("=" * 70)

    print(f"\n9.  Loading all trials across {len(ALL_SESSIONS)} sessions...")
    all_trials_data = collect_all_trials_xcorr(
        ANIMAL_ID, ALL_SESSIONS, FOLDER_SUFFIX, SURROGATE_BAND)

    if all_trials_data is not None:
        print("\n10. Creating surrogate-tested figure...")
        ex = all_trials_data['example_out']
        ex_label = all_trials_data['example_label']
        lags = all_trials_data['lags_ms']
        n_loaded = all_trials_data['n_loaded']

        fig3 = plt.figure(figsize=(20, 12))
        gs_outer = GridSpec(2, 1, figure=fig3, height_ratios=[1.0, 1.0],
                            hspace=0.35, left=0.06, right=0.96,
                            top=0.91, bottom=0.08)
        gs_top = gs_outer[0].subgridspec(
            1, 5, width_ratios=[1.2, 0.04, 0.08, 1.0, 0.04], wspace=0.12)
        gs_bot = gs_outer[1].subgridspec(
            1, 3, width_ratios=[1.0, 1.0, 0.85], wspace=0.28)

        # --- Top-left: example trial heatmap (time x lag) ---
        ax_hm_ex = fig3.add_subplot(gs_top[0])
        cax_ex = fig3.add_subplot(gs_top[1])
        if ex is not None and ex['mat'].size > 0:
            abs_vals = np.abs(ex['mat'].ravel())
            vmax_ex = min(0.5, max(0.2, float(
                np.ceil(np.percentile(abs_vals, 98) * 10) / 10)))
            extent_ex = [ex['centers'][0], ex['centers'][-1],
                         ex['lags_ms'][0], ex['lags_ms'][-1]]
            im_ex = ax_hm_ex.imshow(
                ex['mat'], aspect='auto', origin='lower', extent=extent_ex,
                cmap=CMAP_HEATMAP, vmin=-vmax_ex, vmax=vmax_ex,
                interpolation='bilinear', rasterized=True)
            ax_hm_ex.axhline(0, color='white', ls='--', lw=1.5, alpha=0.8)
            cb = fig3.colorbar(im_ex, cax=cax_ex)
            cb.set_label('r', fontsize=FONT_SIZE_CBAR, fontweight='bold')
            cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 1)
        else:
            cax_ex.set_visible(False)
        ax_hm_ex.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
        ax_hm_ex.set_ylabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
        ax_hm_ex.set_title(
            f'Example trial ({ex_label}, n={ex["n_windows"]} win)',
            fontsize=FONT_SIZE_TITLE, fontweight='bold', loc='left')
        ax_hm_ex.tick_params(labelsize=FONT_SIZE_TICK, width=TICK_WIDTH,
                             length=TICK_LENGTH)

        fig3.add_subplot(gs_top[2]).set_visible(False)

        # --- Top-right: all-trials heatmap (trial x lag) ---
        ax_hm_all = fig3.add_subplot(gs_top[3])
        cax_all = fig3.add_subplot(gs_top[4])
        mat_all = all_trials_data['mat_all']
        if mat_all is not None and mat_all.size > 0:
            abs_all = np.abs(mat_all.ravel())
            vmax_a = min(0.5, max(0.15, float(
                np.round(np.percentile(abs_all, 98), 2))))
            im_all = ax_hm_all.imshow(
                mat_all, aspect='auto', origin='lower',
                extent=[lags[0], lags[-1], 0.5, n_loaded + 0.5],
                cmap='RdBu_r', vmin=-vmax_a, vmax=vmax_a,
                interpolation='nearest', rasterized=True)
            ax_hm_all.axvline(0, color='black', ls='--', lw=1.5, alpha=0.6)
            cb2 = fig3.colorbar(im_all, cax=cax_all)
            cb2.set_label('r', fontsize=FONT_SIZE_CBAR, fontweight='bold')
            cb2.ax.tick_params(labelsize=FONT_SIZE_TICK - 1)
        else:
            cax_all.set_visible(False)
        ax_hm_all.set_xlabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
        ax_hm_all.set_ylabel('Trial', fontsize=FONT_SIZE_LABEL)
        ax_hm_all.set_title(f'All trials (n={n_loaded})',
                            fontsize=FONT_SIZE_TITLE, fontweight='bold',
                            loc='left')
        ax_hm_all.tick_params(labelsize=FONT_SIZE_TICK, width=TICK_WIDTH,
                              length=TICK_LENGTH)

        peak_lim = get_peak_lag_limit(SURROGATE_BAND)

        # --- Bottom-left: example trial mean correlogram ---
        ax_ex_line = fig3.add_subplot(gs_bot[0])
        if ex is not None and ex['mat'].size > 0:
            mean_ex = mean_xcorr_windows(ex['mat'])
            ax_ex_line.plot(ex['lags_ms'], mean_ex, '-', color=COLOR_XCORR,
                            lw=LINE_WIDTH)
            ir_ex, _ = peak_index_restricted(
                mean_ex, ex['lags_ms'], SURROGATE_BAND)
            ax_ex_line.text(
                0.02, 0.95,
                f'|r|={abs(mean_ex[ir_ex]):.3f} @ '
                f'{ex["lags_ms"][ir_ex]:.1f} ms'
                f' (n={ex["n_windows"]} win)',
                transform=ax_ex_line.transAxes,
                fontsize=FONT_SIZE_ANNOTATION,
                ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.6',
                          alpha=0.92))
        ax_ex_line.axvline(0, color='gray', ls='--', lw=1.5, alpha=0.7)
        ax_ex_line.axhline(0, color='gray', ls='-', lw=1, alpha=0.5)
        ax_ex_line.axvline(-peak_lim, color='0.65', ls=':', lw=1.0, alpha=0.7)
        ax_ex_line.axvline(peak_lim, color='0.65', ls=':', lw=1.0, alpha=0.7)
        ax_ex_line.set_xlim(-MAX_LAG_MS, MAX_LAG_MS)
        ax_ex_line.xaxis.set_major_locator(MultipleLocator(50))
        ax_ex_line.set_xlabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
        ax_ex_line.set_ylabel('Cross-correlation', fontsize=FONT_SIZE_LABEL)
        ax_ex_line.set_title(f'Example trial ({ex_label})',
                             fontsize=FONT_SIZE_TITLE, fontweight='bold',
                             pad=8)
        style_axis(ax_ex_line)

        # --- Bottom-center: grand mean +/- SEM ---
        ax_grand = fig3.add_subplot(gs_bot[1])
        gm = all_trials_data['grand_mean']
        sem_arr = all_trials_data['sem']
        if gm is not None:
            ax_grand.plot(lags, gm, '-', color=COLOR_XCORR,
                          lw=LINE_WIDTH + 0.5)
            if sem_arr is not None:
                ax_grand.fill_between(lags, gm - sem_arr, gm + sem_arr,
                                      color=COLOR_XCORR, alpha=0.2)
            ir_g, _ = peak_index_restricted(gm, lags, SURROGATE_BAND)
            ax_grand.text(
                0.02, 0.95,
                f'|r|={abs(gm[ir_g]):.3f} @ {lags[ir_g]:.1f} ms '
                f'(n={n_loaded})',
                transform=ax_grand.transAxes,
                fontsize=FONT_SIZE_ANNOTATION,
                ha='left', va='top',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='0.6',
                          alpha=0.92))
        ax_grand.axvline(0, color='gray', ls='--', lw=1.5, alpha=0.7)
        ax_grand.axhline(0, color='gray', ls='-', lw=1, alpha=0.5)
        ax_grand.axvline(-peak_lim, color='0.65', ls=':', lw=1.0, alpha=0.7)
        ax_grand.axvline(peak_lim, color='0.65', ls=':', lw=1.0, alpha=0.7)
        ax_grand.set_xlim(-MAX_LAG_MS, MAX_LAG_MS)
        ax_grand.xaxis.set_major_locator(MultipleLocator(50))
        ax_grand.set_xlabel('Temporal Offset (ms)', fontsize=FONT_SIZE_LABEL)
        ax_grand.set_ylabel('Cross-correlation', fontsize=FONT_SIZE_LABEL)
        ax_grand.set_title(
            f'All trials (n={n_loaded}): mean \u00b1 SEM',
            fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=8)
        style_axis(ax_grand)

        # --- Bottom-right: surrogate test panel ---
        ax_surr = fig3.add_subplot(gs_bot[2])
        plot_surrogate_test_panel(ax_surr, all_trials_data)

        fig3.suptitle(
            f'GEVI-LFP Theta Cross-Correlation \u2014 {ANIMAL_ID} '
            f'({n_loaded} trials, {len(ALL_SESSIONS)} sessions, '
            f'surrogate-tested)',
            fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.97)

        fig3_name = (f"cross_correlation_{ANIMAL_ID}_allSessions"
                     f"_surrogate_tested")
        print("\n11. Saving surrogate-tested figure...")
        for fmt in FIGURE_FORMATS:
            out_path = OUTPUT_DIR / f"{fig3_name}.{fmt}"
            fig3.savefig(str(out_path), format=fmt, dpi=FIGURE_DPI,
                         bbox_inches='tight', facecolor='white')
            print(f"    {out_path}")
        plt.close(fig3)

        # --- Print results summary for copy-paste ---
        pooled_p = all_trials_data['pooled_p']
        n_sig = all_trials_data['n_sig_trials']
        n_total = all_trials_data['n_loaded']
        median_r = float(np.median(all_trials_data['observed_peaks']))
        mean_lag = float(np.mean(
            [t['peak_lag'] for t in all_trials_data['trials']]))
        null_p95 = float(np.percentile(
            all_trials_data['pooled_null'], 95))
        p_str = (f"p = {pooled_p:.4f}"
                 if np.isfinite(pooled_p) and pooled_p >= 0.001
                 else f"p = {pooled_p:.2e}"
                 if np.isfinite(pooled_p) else "p = n/a")

        print("\n" + "=" * 70)
        print("SURROGATE TEST RESULTS")
        print("=" * 70)
        print(f"  Animal:             {ANIMAL_ID}")
        print(f"  Sessions:           {len(ALL_SESSIONS)}")
        print(f"  Trials loaded:      {n_total}")
        print(f"  Band:               {SURROGATE_BAND[0]}-"
              f"{SURROGATE_BAND[1]} Hz (theta)")
        print(f"  Surrogates/trial:   {N_SURROGATES}")
        print(f"  Median peak |r|:    {median_r:.3f}")
        print(f"  Mean peak lag:      {mean_lag:.1f} ms")
        print(f"  Trials significant: {n_sig}/{n_total} "
              f"(p<0.05, per-trial)")
        print(f"  Pooled {p_str}")
        print(f"  Null 95th pctl:     {null_p95:.3f}")
        print("-" * 70)
        print(f"RESULT: GEVI-LFP theta-band cross-correlation "
              f"({ANIMAL_ID}, n={n_total} trials across "
              f"{len(ALL_SESSIONS)} sessions): "
              f"median peak |r| = {median_r:.3f} at mean lag "
              f"{mean_lag:.1f} ms, "
              f"{n_sig}/{n_total} trials individually significant, "
              f"pooled {p_str} "
              f"(circular-shift surrogate test, "
              f"{N_SURROGATES} surrogates/trial)")
        print("=" * 70)
    else:
        print("\n  ERROR: No trials loaded — surrogate test skipped.")

    print("\nDONE")


if __name__ == "__main__":
    main()
