"""
Publication composite figures for DBS stimulation analysis.

Generates Figures 1-4 (a=Animal01, b=Animal02) for amp-balanced and energy-balanced
DBS conditions, following the same 6-panel format as the Striatum DBS composite.

Figures:
  1a/1b: 135 Hz condition (6 panels: A-F)
  2a/2b: 40 Hz Amp-balanced condition (6 panels: A-F)
  3a/3b: 40 Hz Energy-balanced condition (6 panels: A-F)
  4a/4b: All conditions composite (averaged fiber + trial heatmaps for all 3 conditions)
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import loadmat
from scipy import signal
from scipy.signal import butter, filtfilt
from scipy import stats
from scipy.ndimage import uniform_filter1d
from pathlib import Path
import warnings
import h5py
import os
import re
import sys

warnings.filterwarnings('ignore')

# common.py lives in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR = PROJECT_ROOT / "Figures" / "CompFig"

DPI = 300
FONT_SIZE = 14
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 22
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12

AXIS_LINEWIDTH = 2.5
TICK_WIDTH = 2.0
TICK_LENGTH = 8
LINE_WIDTH_TRACE = 1.2
LINE_WIDTH_THICK = 2.0

COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])
COLOR_LFP = np.array([0.35, 0.25, 0.45])
COLOR_MOTION = np.array([0.993248, 0.7, 0.4])
COLOR_STIM_PULSE = np.array([0.5, 0.1, 0.1])
COLOR_TRANSIENT = np.array([0.7, 0.15, 0.15])
COLOR_SUSTAINED = np.array([0.65, 0.25, 0.15])
COLOR_POST = np.array([0.15, 0.55, 0.55])

WHEEL_DIAMETER_CM = 19.0
WHEEL_CIRCUMFERENCE_CM = np.pi * WHEEL_DIAMETER_CM
ENCODER_COUNTS_PER_REV = 1024
EPHYS_SAMPLING_RATE = 30000
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM


from common import _infer_trial_from_name, create_parula_like_cmap, normalize_unc_path  # shared helpers (were local copies)


CMAP_PARULA_LIKE = create_parula_like_cmap()

# =============================================================================
# SESSION DEFINITIONS
# =============================================================================

SESSIONS = {
    'Animal03': {
        '135Hz': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal03" / "Fiber_Voltage_Processed",
            'session_id': "01_03_26-R3",
            'label': "135 Hz",
            'stim_freq': 135.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (130.0, 140.0),
        },
        '40Hz_Amp': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal03" / "Fiber_Voltage_Processed",
            'session_id': "01_03_26-R5",
            'label': "40 Hz (Amp-balanced)",
            'stim_freq': 40.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (35.0, 45.0),
        },
        '40Hz_Energy': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal03" / "Fiber_Voltage_Processed",
            'session_id': "01_03_26-R6",
            'label': "40 Hz (Energy-balanced)",
            'stim_freq': 40.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (35.0, 45.0),
        },
    },
    'Animal04': {
        '135Hz': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal04" / "Fiber_Voltage_Processed",
            'session_id': "01_04_26-R7",
            'label': "135 Hz",
            'stim_freq': 135.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (130.0, 140.0),
        },
        '40Hz_Amp': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal04" / "Fiber_Voltage_Processed",
            'session_id': "01_04_26-R10",
            'label': "40 Hz (Amp-balanced)",
            'stim_freq': 40.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (35.0, 45.0),
        },
        '40Hz_Energy': {
            'base_path': _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal04" / "Fiber_Voltage_Processed",
            'session_id': "01_04_26-R11",
            'label': "40 Hz (Energy-balanced)",
            'stim_freq': 40.0,
            'pre_sec': 10.0,
            'stim_sec': 10.0,
            'post_sec': 10.0,
            'band': (35.0, 45.0),
        },
    },
}

NUM_TRIALS = 10

# =============================================================================
# DATA LOADING
# =============================================================================



def load_session_trials(base_path, session_id, mouse_id, num_trials=NUM_TRIALS, fiber_index=0):
    """Load all trial .mat files for a session.
    
    Searches in two locations:
    1. Directly in the session folder (mat files at session level)
    2. Inside Trial{N}_* subdirectories (mat files in trial subfolders)
    """
    session_path = base_path / session_id
    session_str = normalize_unc_path(str(session_path))

    all_mats = []

    # Strategy 1: Try globbing directly in session folder
    for sp in [Path(normalize_unc_path(str(session_path), for_access=True)), Path(session_str)]:
        try:
            found = list(sp.glob("*_Trial*_FiberPhotometry_Analysis.mat"))
            if found:
                all_mats = sorted(found, key=lambda p: (_infer_trial_from_name(p.name) or 999))
                break
        except Exception:
            continue

    # Strategy 2: Search inside Trial{N}_* subdirectories
    if len(all_mats) == 0:
        print(f"    (Searching trial subdirectories for {session_id}...)")
        for sp in [Path(normalize_unc_path(str(session_path), for_access=True)), Path(session_str)]:
            try:
                if not sp.exists():
                    continue
                for trial_num in range(1, num_trials + 1):
                    trial_dirs = list(sp.glob(f"Trial{trial_num}_*"))
                    if not trial_dirs:
                        trial_dirs = list(sp.glob(f"Trial{trial_num}"))
                    for td in trial_dirs:
                        if td.is_dir():
                            mats_in_dir = list(td.glob("*_FiberPhotometry_Analysis.mat"))
                            all_mats.extend(mats_in_dir)
                if len(all_mats) > 0:
                    all_mats = sorted(all_mats, key=lambda p: (_infer_trial_from_name(p.name) or 999))
                    break
            except Exception as e:
                print(f"    (Subdir search error: {e})")
                continue

    # Strategy 3: Use os.listdir for robust directory listing (avoids glob issues on UNC)
    if len(all_mats) == 0:
        print(f"    (Trying os.listdir approach for {mouse_id}-{session_id}...)")
        for sp_str in [normalize_unc_path(str(session_path), for_access=True), session_str]:
            try:
                if not os.path.isdir(sp_str):
                    continue
                entries = os.listdir(sp_str)
                # Check for .mat files directly
                mat_files = [e for e in entries if e.endswith("_FiberPhotometry_Analysis.mat")]
                if mat_files:
                    all_mats = [Path(sp_str) / f for f in sorted(mat_files,
                                key=lambda x: (_infer_trial_from_name(x) or 999))]
                    break
                # Check inside Trial subdirectories
                trial_dirs = [e for e in entries if e.startswith("Trial") and os.path.isdir(os.path.join(sp_str, e))]
                for td_name in sorted(trial_dirs):
                    td_path = os.path.join(sp_str, td_name)
                    try:
                        sub_entries = os.listdir(td_path)
                        sub_mats = [f for f in sub_entries if f.endswith("_FiberPhotometry_Analysis.mat")]
                        for sm in sub_mats:
                            all_mats.append(Path(td_path) / sm)
                    except Exception:
                        continue
                if len(all_mats) > 0:
                    all_mats = sorted(all_mats, key=lambda p: (_infer_trial_from_name(p.name) or 999))
                    break
            except Exception as e:
                print(f"    (os.listdir error: {e})")
                continue

    all_mats = all_mats[:num_trials]

    trials = []
    for p in all_mats:
        try:
            d = _load_mat_file(str(p), fiber_index=fiber_index)
            d["trial_num"] = _infer_trial_from_name(p.name)
            trials.append(d)
            print(f"    Loaded: {p.name}")
        except Exception as e:
            print(f"    Warning: {p.name}: {e}")
    return trials


def _load_mat_file(mat_path, fiber_index=0):
    """Load one FiberPhotometry_Analysis.mat file."""
    mat_access = normalize_unc_path(str(mat_path), for_access=True)
    if not os.path.exists(mat_access):
        mat_access = normalize_unc_path(str(mat_path))
    if not os.path.exists(mat_access):
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    def _safe_scalar(arr, default=np.nan):
        try:
            arr = np.asarray(arr).ravel()
            return float(arr[0]) if arr.size else float(default)
        except Exception:
            return float(default)

    try:
        with h5py.File(mat_access, "r") as f:
            root = f["FiberPhotometryAnalysis"]
            t = np.array(root["time"]["time_vector_seconds"][()]).flatten()
            fs = _safe_scalar(root["time"]["sampling_rate"][()], default=np.nan)
            if not np.isfinite(fs) or fs <= 0:
                fs = 1.0 / np.nanmedian(np.diff(t))

            tr = np.array(root["signals"]["final_processed_traces"][()])
            if tr.ndim == 2:
                if tr.shape[0] < tr.shape[1]:
                    tr = tr.T
                fi = min(fiber_index, tr.shape[1] - 1)
                gevi = tr[:, fi]
            else:
                gevi = tr.flatten()

            ephys_grp = root["ephys"]
            if "lfp_raw_aligned_HP" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_HP"][()]).flatten()
            elif "lfp_raw_aligned_mPFC" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_mPFC"][()]).flatten()
            else:
                lfp = np.zeros_like(t)

            if "running_velocity_smooth" in ephys_grp:
                motion = np.array(ephys_grp["running_velocity_smooth"][()]).flatten()
            elif "running_velocity" in ephys_grp:
                motion = np.array(ephys_grp["running_velocity"][()]).flatten()
            else:
                motion = np.zeros_like(t)

            n = min(len(t), len(gevi), len(lfp), len(motion))
            return {
                "t": t[:n],
                "gevi": gevi[:n],
                "lfp": lfp[:n],
                "motion": motion[:n] * MOTION_TO_CM_PER_S,
                "fs": float(fs),
            }
    except Exception:
        pass

    mat = loadmat(mat_access, squeeze_me=False, struct_as_record=False)
    fp = mat["FiberPhotometryAnalysis"]
    t = np.asarray(fp.time.time_vector_seconds).flatten()
    fs = _safe_scalar(fp.time.sampling_rate, default=np.nan)
    if not np.isfinite(fs) or fs <= 0:
        fs = 1.0 / np.nanmedian(np.diff(t))

    tr = np.asarray(fp.signals.final_processed_traces)
    if tr.ndim == 2:
        if tr.shape[0] < tr.shape[1]:
            tr = tr.T
        fi = min(fiber_index, tr.shape[1] - 1)
        gevi = tr[:, fi]
    else:
        gevi = tr.flatten()

    ephys = fp.ephys
    if hasattr(ephys, 'lfp_raw_aligned_HP'):
        lfp = np.asarray(ephys.lfp_raw_aligned_HP).flatten()
    elif hasattr(ephys, 'lfp_raw_aligned_mPFC'):
        lfp = np.asarray(ephys.lfp_raw_aligned_mPFC).flatten()
    else:
        lfp = np.zeros_like(t)

    if hasattr(ephys, 'running_velocity_smooth'):
        motion = np.asarray(ephys.running_velocity_smooth).flatten()
    else:
        motion = np.asarray(ephys.running_velocity).flatten()

    n = min(len(t), len(gevi), len(lfp), len(motion))
    return {
        "t": t[:n],
        "gevi": gevi[:n],
        "lfp": lfp[:n],
        "motion": motion[:n] * MOTION_TO_CM_PER_S,
        "fs": float(fs),
    }


# =============================================================================
# ANALYSIS HELPERS
# =============================================================================

def generate_stim_pulses(t, stim_freq_hz, stim_duration_sec):
    """Generate square biphasic stimulation pulse train."""
    pulses = np.zeros_like(t)
    if stim_freq_hz <= 0 or stim_duration_sec <= 0:
        return pulses
    mask = (t >= 0.0) & (t <= stim_duration_sec)
    if not np.any(mask):
        return pulses
    period = 1.0 / stim_freq_hz
    t_stim = t[mask] - t[mask][0]
    phase = (t_stim % period) / period
    pulses[mask] = np.where(phase < 0.3, 1.0, np.where(phase < 0.6, -1.0, 0.0))
    return pulses


def _compute_trial_band_power(sig, fs, band_hz):
    """Compute band power in a frequency band using Welch's method."""
    x = np.asarray(sig, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 64:
        return np.nan
    nperseg = min(256, len(x))
    f, pxx = signal.welch(x, fs=fs, nperseg=nperseg, noverlap=nperseg // 2)
    mask = (f >= band_hz[0]) & (f <= band_hz[1])
    if not np.any(mask):
        return np.nan
    return float(np.trapz(pxx[mask], f[mask]))


def _plot_three_period_violin_with_stats(ax, data_dict, ylabel, title=None):
    """Plot violin/strip plot for Transient, Sustained, Post-stim periods with stats."""
    order = ["Transient", "Sustained", "Post-stim"]
    colors = {"Transient": COLOR_TRANSIENT, "Sustained": COLOR_SUSTAINED, "Post-stim": COLOR_POST}

    positions = []
    all_data = []
    labels = []
    for i, key in enumerate(order):
        vals = np.array(data_dict.get(key, []))
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            positions.append(i)
            all_data.append(vals)
            labels.append(key)

    if len(all_data) == 0:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center")
        return

    for i, (pos, vals, lab) in enumerate(zip(positions, all_data, labels)):
        color = colors[lab]
        if len(vals) >= 4:
            parts = ax.violinplot([vals], positions=[pos], showmeans=False, showextrema=False, widths=0.7)
            for pc in parts['bodies']:
                pc.set_facecolor(color)
                pc.set_alpha(0.3)
        jitter = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(vals))
        ax.scatter(np.full(len(vals), pos) + jitter, vals, color=color, s=30, alpha=0.7, edgecolors='none', zorder=3)
        ax.plot([pos - 0.15, pos + 0.15], [np.mean(vals), np.mean(vals)], color=color, linewidth=2.5, zorder=4)

    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8, alpha=0.6)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    # Stats: paired t-tests between adjacent periods
    fs_annot = FONT_SIZE - 2
    for i in range(len(all_data) - 1):
        if len(all_data[i]) >= 3 and len(all_data[i + 1]) >= 3:
            n_min = min(len(all_data[i]), len(all_data[i + 1]))
            _, p_val = stats.ttest_ind(all_data[i][:n_min], all_data[i + 1][:n_min])
            if p_val < 0.001:
                sig_text = "***"
            elif p_val < 0.01:
                sig_text = "**"
            elif p_val < 0.05:
                sig_text = "*"
            else:
                sig_text = "ns"
            y_max = max(np.max(all_data[i]), np.max(all_data[i + 1]))
            y_bar = y_max + 0.1 * (ax.get_ylim()[1] - ax.get_ylim()[0])
            ax.plot([positions[i], positions[i + 1]], [y_bar, y_bar], 'k-', linewidth=1.2)
            ax.text((positions[i] + positions[i + 1]) / 2, y_bar, sig_text,
                    ha='center', va='bottom', fontsize=fs_annot)


# =============================================================================
# FIGURE 1-3: 6-PANEL COMPOSITE (same format as Striatum DBS)
# =============================================================================

def create_six_panel_composite(trials, session_cfg, fig_label_prefix=""):
    """
    Create 6-panel composite figure for a single DBS condition.
    Panels: A=Single trial overview, B=Stim-onset zoom, C=Trial-averaged fiber,
            D=Trial-averaged spectrogram, E=Period violin (fiber), F=Period violin (band power)
    """
    pre_sec = session_cfg['pre_sec']
    stim_sec = session_cfg['stim_sec']
    post_sec = session_cfg['post_sec']
    band = session_cfg['band']
    stim_freq = session_cfg['stim_freq']
    transient_end = 1.0
    sustained_start = 1.0

    if len(trials) == 0:
        fig = plt.figure(figsize=(42, 30))
        fig.text(0.5, 0.5, "No data available", ha='center', va='center', fontsize=30)
        return fig

    rep = trials[0]
    for tr in trials:
        if tr.get("trial_num") == 1:
            rep = tr
            break

    fs_rep = rep["fs"]
    stim_on_idx = int(np.argmin(np.abs(rep["t"] - (rep["t"][0] + pre_sec))))
    t_rel = rep["t"] - rep["t"][stim_on_idx]
    fiber_pct = rep["gevi"] * 100.0
    motion = rep["motion"]
    stim_trace = generate_stim_pulses(t_rel, stim_freq, stim_sec)

    t_grid = np.linspace(-pre_sec, stim_sec + post_sec, 2500)
    fiber_all = []
    period_fiber = {"Transient": [], "Sustained": [], "Post-stim": []}
    period_band_db = {"Transient": [], "Sustained": [], "Post-stim": []}
    spec_stack = []

    for tr in trials:
        fs = tr["fs"]
        sidx = int(np.argmin(np.abs(tr["t"] - (tr["t"][0] + pre_sec))))
        tt = tr["t"] - tr["t"][sidx]
        yy = tr["gevi"] * 100.0
        keep = (tt >= -pre_sec) & (tt <= (stim_sec + post_sec))
        if np.sum(keep) < 100:
            continue
        tt = tt[keep]
        yy = yy[keep]
        yi = np.interp(t_grid, tt, yy, left=np.nan, right=np.nan)
        fiber_all.append(yi)

        bmask = (tt >= -pre_sec) & (tt < 0.0)
        baseline = np.nanmean(yy[bmask]) if np.any(bmask) else np.nan
        if np.isfinite(baseline):
            tmask = (tt >= 0.0) & (tt < transient_end)
            smask = (tt >= sustained_start) & (tt < stim_sec)
            pmask = (tt >= stim_sec) & (tt < (stim_sec + post_sec))
            if np.any(tmask):
                period_fiber["Transient"].append(float(np.nanmean(yy[tmask]) - baseline))
            if np.any(smask):
                period_fiber["Sustained"].append(float(np.nanmean(yy[smask]) - baseline))
            if np.any(pmask):
                period_fiber["Post-stim"].append(float(np.nanmean(yy[pmask]) - baseline))

        bpow = _compute_trial_band_power(yy[(tt >= -pre_sec) & (tt < 0.0)], fs, band)
        if np.isfinite(bpow) and bpow > 0:
            tpow = _compute_trial_band_power(yy[(tt >= 0.0) & (tt < transient_end)], fs, band)
            spow = _compute_trial_band_power(yy[(tt >= sustained_start) & (tt < stim_sec)], fs, band)
            ppow = _compute_trial_band_power(yy[(tt >= stim_sec) & (tt < (stim_sec + post_sec))], fs, band)
            if np.isfinite(tpow) and tpow > 0:
                period_band_db["Transient"].append(float(10 * np.log10(tpow / bpow)))
            if np.isfinite(spow) and spow > 0:
                period_band_db["Sustained"].append(float(10 * np.log10(spow / bpow)))
            if np.isfinite(ppow) and ppow > 0:
                period_band_db["Post-stim"].append(float(10 * np.log10(ppow / bpow)))

        nper = min(512, len(yy))
        if nper < 128:
            continue
        nover = min(int(0.875 * nper), nper - 1)
        nfft = max(2048, int(2 ** np.ceil(np.log2(nper))))
        f, ts, sxx = signal.spectrogram(
            np.nan_to_num(yy - np.nanmean(yy), nan=0.0),
            fs=fs, window="hann", nperseg=nper, noverlap=nover,
            nfft=nfft, detrend="linear", scaling="density", mode="psd",
        )
        ts_rel = ts + (tt[0] if len(tt) else 0.0)
        fmk = (f >= 1.0) & (f <= 100.0)
        if np.any(fmk):
            f2 = f[fmk]
            s2 = sxx[fmk, :]
            bmk = (ts_rel >= -pre_sec) & (ts_rel < -1.0)
            if np.sum(bmk) == 0:
                bmk = ts_rel < 0.0
            if np.sum(bmk) > 0:
                base = np.nanmean(s2[:, bmk], axis=1, keepdims=True)
                base = np.where(np.isfinite(base) & (base > 0), base, np.nan)
                s_rel = 10 * np.log10(np.maximum(s2, 1e-15) / np.maximum(base, 1e-15))
                if np.any(np.isfinite(s_rel)):
                    spec_stack.append((f2, ts_rel, s_rel))

    fiber_mat = np.asarray(fiber_all, dtype=float) if len(fiber_all) else np.empty((0, len(t_grid)))
    fiber_mu = np.nanmean(fiber_mat, axis=0) if fiber_mat.size else np.full(len(t_grid), np.nan)
    if fiber_mat.shape[0] > 1:
        fiber_sem = np.nanstd(fiber_mat, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(fiber_mat), axis=0))
    else:
        fiber_sem = np.full(len(t_grid), np.nan)

    spec_avg = None
    if len(spec_stack):
        f_ref = spec_stack[0][0]
        t_ref = spec_stack[0][1]
        mats = []
        for ff, tt_s, ss in spec_stack:
            if len(ff) != len(f_ref) or np.max(np.abs(ff - f_ref)) > 1e-9:
                continue
            arr = np.full((len(f_ref), len(t_ref)), np.nan)
            for i in range(len(f_ref)):
                vv = ss[i, :]
                mk = np.isfinite(tt_s) & np.isfinite(vv)
                if np.sum(mk) >= 2:
                    arr[i, :] = np.interp(t_ref, tt_s[mk], vv[mk], left=np.nan, right=np.nan)
            mats.append(arr)
        if len(mats):
            spec_avg = (f_ref, t_ref, np.nanmean(np.asarray(mats), axis=0))

    # Build figure
    fig = plt.figure(figsize=(42, 30))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.1, 1.1, 1.6],
                  width_ratios=[1.12, 1.0], hspace=0.56, wspace=0.30)
    fs_label = FONT_SIZE_LABEL + 28
    fs_tick = FONT_SIZE_TICK + 28

    x_lo = -min(pre_sec, 5.0)
    x_hi = stim_sec + min(post_sec, 5.0)

    # Panel A: Single trial overview
    axA = fig.add_subplot(gs[0, 0])
    mA = (t_rel >= x_lo) & (t_rel <= x_hi)
    tA = t_rel[mA]
    y_f = fiber_pct[mA] - np.nanmedian(fiber_pct[mA])
    y_m = motion[mA] - np.nanmedian(motion[mA]) if np.any(np.isfinite(motion[mA])) else np.zeros_like(tA)
    y_s = stim_trace[mA]
    nf = y_f / (np.nanpercentile(np.abs(y_f), 95) + 1e-9)
    nm = y_m / (np.nanpercentile(np.abs(y_m), 95) + 1e-9) if np.any(np.isfinite(y_m)) else np.zeros_like(tA)
    ns = y_s / (np.nanpercentile(np.abs(y_s), 95) + 1e-9) if np.any(np.isfinite(y_s)) else np.zeros_like(tA)
    off_f, off_m, off_s = 0.0, -2.8, 2.8
    axA.plot(tA, ns * 0.6 + off_s, color=COLOR_STIM_PULSE, linewidth=1.5)
    axA.plot(tA, nf * 0.9 + off_f, color=COLOR_GEVI, linewidth=2.1)
    motion_color_dark = np.clip(np.asarray(COLOR_MOTION) * 0.58, 0.0, 1.0)
    axA.plot(tA, nm * 0.9 + off_m, color=motion_color_dark, linewidth=1.6)
    axA.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axA.axvline(stim_sec, color="k", linestyle=":", linewidth=1.2, alpha=0.8)
    axA.set_xlim(x_lo, x_hi)
    axA.set_ylim(-4.3, 4.3)
    axA.set_yticks([off_m, off_f, off_s])
    axA.set_yticklabels(["Motion", "Fiber Vm", "Stim"], fontsize=fs_tick)
    axA.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")

    # Panel B: Stim-onset zoom
    axB = fig.add_subplot(gs[0, 1])
    zoom_pre, zoom_post = 0.5, 0.5
    mB = (t_rel >= -zoom_pre) & (t_rel <= zoom_post)
    tB = t_rel[mB]
    y_fb = fiber_pct[mB] - np.nanmedian(fiber_pct[mB])
    y_sb = stim_trace[mB]
    nfb = y_fb / (np.nanpercentile(np.abs(y_fb), 95) + 1e-9)
    nsb = y_sb / (np.nanpercentile(np.abs(y_sb), 95) + 1e-9) if np.any(np.isfinite(y_sb)) else np.zeros_like(tB)
    axB.plot(tB, nsb * 0.5 + 1.0, color=COLOR_STIM_PULSE, linewidth=1.4)
    axB.plot(tB, nfb * 0.9 - 0.8, color=COLOR_GEVI, linewidth=2.8)
    axB.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axB.set_xlim(-zoom_pre, zoom_post)
    axB.set_ylim(-2.0, 2.0)
    axB.set_yticks([1.0, -0.8])
    axB.set_yticklabels(["Stim", "Fiber Vm"], fontsize=fs_tick)
    axB.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")

    # Panel C: Trial-averaged fiber
    axC = fig.add_subplot(gs[1, 0])
    if np.any(np.isfinite(fiber_mu)):
        axC.plot(t_grid, fiber_mu, color=COLOR_GEVI, linewidth=2.2)
        if np.any(np.isfinite(fiber_sem)):
            axC.fill_between(t_grid, fiber_mu - fiber_sem, fiber_mu + fiber_sem,
                             color=COLOR_GEVI, alpha=0.28, linewidth=0)
    axC.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axC.axvline(stim_sec, color="k", linestyle=":", linewidth=1.2, alpha=0.8)
    axC.set_xlim(x_lo, x_hi)
    axC.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    axC.set_ylabel("Population Vm (\u0394F/F %)", fontsize=fs_label, fontweight="bold")

    # Panel D: Trial-averaged spectrogram
    axD = fig.add_subplot(gs[1, 1])
    if spec_avg is not None:
        ff, tt_s, ss = spec_avg
        ss_disp = uniform_filter1d(uniform_filter1d(ss, size=3, axis=0, mode="nearest"), size=5, axis=1, mode="nearest")
        vmax = np.nanpercentile(np.abs(ss_disp[np.isfinite(ss_disp)]), 95) if np.any(np.isfinite(ss_disp)) else 2.0
        im = axD.pcolormesh(tt_s, ff, ss_disp, cmap=CMAP_PARULA_LIKE, shading="auto", vmin=-vmax, vmax=vmax)
        axD.set_ylim(1, 100)
        cbar = fig.colorbar(im, ax=axD, pad=0.02)
        cbar.set_label("Power (dB)", fontsize=fs_tick)
        cbar.ax.tick_params(labelsize=max(fs_tick - 2, 1))
    else:
        axD.text(0.5, 0.5, "No spectrogram data", transform=axD.transAxes, ha="center", va="center", fontsize=fs_label)
    axD.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
    axD.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
    axD.set_xlim(x_lo, x_hi)
    axD.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    axD.set_ylabel("Frequency (Hz)", fontsize=fs_label, fontweight="bold")

    # Panel E: Period violin (fiber)
    axE = fig.add_subplot(gs[2, 0])
    _plot_three_period_violin_with_stats(axE, period_fiber, ylabel="Population Vm change (\u0394F/F %)")
    axE.set_ylabel("Population Vm change (\u0394F/F %)", fontsize=fs_label, fontweight="bold")

    # Panel F: Period violin (band power)
    axF = fig.add_subplot(gs[2, 1])
    band_label = f"{int(band[0])}-{int(band[1])} Hz"
    _plot_three_period_violin_with_stats(axF, period_band_db, ylabel=f"Relative {band_label} band power (dB)")
    axF.set_ylabel(f"Relative {band_label} band power (dB)", fontsize=fs_label, fontweight="bold")

    # Style all axes
    panel_labels = ["A", "B", "C", "D", "E", "F"]
    axes = [axA, axB, axC, axD, axE, axF]
    for ax, lbl in zip(axes, panel_labels):
        ax.tick_params(labelsize=fs_tick, width=TICK_WIDTH, length=TICK_LENGTH)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
        ax.text(-0.08, 1.05, f"{fig_label_prefix}{lbl}", transform=ax.transAxes,
                fontsize=fs_label + 8, fontweight='bold', va='top', ha='right')

    for ax in [axE, axF]:
        for txt in ax.texts:
            if txt.get_text() in ["***", "**", "*", "ns"]:
                txt.set_fontsize(txt.get_fontsize() + 8)

    fig.subplots_adjust(left=0.06, right=0.965, top=0.985, bottom=0.055)
    return fig


# =============================================================================
# FIGURE 4: ALL-CONDITIONS COMPOSITE (averaged fiber + trial heatmaps)
# =============================================================================

def create_all_conditions_composite(all_trials_by_condition, session_cfgs, mouse_label=""):
    """
    Create the all-conditions composite (Figure 4).
    6 rows: for each of 3 conditions, row1=averaged traces, row2=per-trial heatmaps.
    2 columns: col1=fiber signal, col2=stim band power.
    Rows 1,3,5 (averaged) have less height than rows 2,4,6 (heatmaps).
    """
    conditions_order = ['135Hz', '40Hz_Amp', '40Hz_Energy']
    n_cond = len(conditions_order)

    fig = plt.figure(figsize=(28, 38))
    gs = GridSpec(
        n_cond * 2, 2, figure=fig,
        height_ratios=[0.6, 1.0, 0.6, 1.0, 0.6, 1.0],
        width_ratios=[1.0, 1.0],
        hspace=0.35, wspace=0.30
    )

    panel_idx = 0
    panel_labels = "ABCDEFGHIJKL"
    fs_label = FONT_SIZE_LABEL + 6
    fs_tick = FONT_SIZE_TICK + 4
    fs_title = FONT_SIZE_TITLE + 2

    for cond_idx, cond_key in enumerate(conditions_order):
        cfg = session_cfgs[cond_key]
        trials = all_trials_by_condition[cond_key]
        pre_sec = cfg['pre_sec']
        stim_sec = cfg['stim_sec']
        post_sec = cfg['post_sec']
        band = cfg['band']
        stim_freq = cfg['stim_freq']
        label = cfg['label']

        x_lo = -min(pre_sec, 5.0)
        x_hi = stim_sec + min(post_sec, 5.0)
        t_grid = np.linspace(x_lo, x_hi, 2000)
        n_bins = 100

        fiber_all = []
        bandpower_all = []

        for tr in trials:
            fs = tr["fs"]
            sidx = int(np.argmin(np.abs(tr["t"] - (tr["t"][0] + pre_sec))))
            tt = tr["t"] - tr["t"][sidx]
            yy = tr["gevi"] * 100.0
            keep = (tt >= x_lo) & (tt <= x_hi)
            if np.sum(keep) < 100:
                continue
            tt_k = tt[keep]
            yy_k = yy[keep]

            # Baseline subtract
            bmask = tt_k < 0.0
            baseline = np.nanmean(yy_k[bmask]) if np.any(bmask) else 0.0
            yy_bs = yy_k - baseline

            # Interpolate to common grid
            yi = np.interp(t_grid, tt_k, yy_bs, left=np.nan, right=np.nan)
            fiber_all.append(yi)

            # Band power in sliding windows
            bin_edges = np.linspace(x_lo, x_hi, n_bins + 1)
            bp_trace = np.full(n_bins, np.nan)
            for bi in range(n_bins):
                win_mask = (tt_k >= bin_edges[bi]) & (tt_k < bin_edges[bi + 1])
                seg = yy_k[win_mask]
                if len(seg) >= 32:
                    bp_trace[bi] = _compute_trial_band_power(seg, fs, band)
            # Convert to dB relative to baseline bins
            bl_bins = bin_edges[:-1] < 0.0
            bl_vals = bp_trace[bl_bins]
            bl_vals = bl_vals[np.isfinite(bl_vals) & (bl_vals > 0)]
            if len(bl_vals) > 0:
                bp_ref = np.mean(bl_vals)
                bp_trace = 10 * np.log10(np.maximum(bp_trace, 1e-15) / bp_ref)
            bandpower_all.append(bp_trace)

        fiber_mat = np.asarray(fiber_all) if len(fiber_all) else np.full((1, len(t_grid)), np.nan)
        bp_mat = np.asarray(bandpower_all) if len(bandpower_all) else np.full((1, n_bins), np.nan)
        bin_centers = (np.linspace(x_lo, x_hi, n_bins + 1)[:-1] + np.linspace(x_lo, x_hi, n_bins + 1)[1:]) / 2

        # Row 1 (averaged traces): fiber avg (col1) + band power avg (col2)
        row_avg = cond_idx * 2
        row_heat = cond_idx * 2 + 1

        # Averaged fiber signal
        ax_avg_fiber = fig.add_subplot(gs[row_avg, 0])
        fiber_mu = np.nanmean(fiber_mat, axis=0)
        fiber_sem = np.nanstd(fiber_mat, axis=0, ddof=1) / np.sqrt(fiber_mat.shape[0]) if fiber_mat.shape[0] > 1 else np.zeros_like(fiber_mu)
        ax_avg_fiber.plot(t_grid, fiber_mu, color=COLOR_GEVI, linewidth=2.0)
        ax_avg_fiber.fill_between(t_grid, fiber_mu - fiber_sem, fiber_mu + fiber_sem,
                                  color=COLOR_GEVI, alpha=0.25, linewidth=0)
        ax_avg_fiber.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
        ax_avg_fiber.axvline(stim_sec, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
        ax_avg_fiber.set_xlim(x_lo, x_hi)
        ax_avg_fiber.set_ylabel("\u0394F/F (%)", fontsize=fs_label, fontweight="bold")
        ax_avg_fiber.set_title(f"{label} - Averaged Fiber Signal (n={fiber_mat.shape[0]})",
                               fontsize=fs_title, fontweight="bold", pad=8)
        ax_avg_fiber.text(-0.08, 1.12, panel_labels[panel_idx], transform=ax_avg_fiber.transAxes,
                          fontsize=fs_label + 8, fontweight='bold', va='top', ha='right')
        panel_idx += 1

        # Averaged band power
        ax_avg_bp = fig.add_subplot(gs[row_avg, 1])
        bp_mu = np.nanmean(bp_mat, axis=0)
        bp_sem = np.nanstd(bp_mat, axis=0, ddof=1) / np.sqrt(bp_mat.shape[0]) if bp_mat.shape[0] > 1 else np.zeros_like(bp_mu)
        ax_avg_bp.plot(bin_centers, bp_mu, color=COLOR_GEVI, linewidth=2.0)
        ax_avg_bp.fill_between(bin_centers, bp_mu - bp_sem, bp_mu + bp_sem,
                               color=COLOR_GEVI, alpha=0.25, linewidth=0)
        ax_avg_bp.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
        ax_avg_bp.axvline(stim_sec, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
        ax_avg_bp.axhline(0.0, color="grey", linestyle="--", linewidth=0.7, alpha=0.5)
        ax_avg_bp.set_xlim(x_lo, x_hi)
        band_label = f"{int(band[0])}-{int(band[1])}" if stim_freq <= 50 else f"{int(stim_freq)}\u00b15"
        ax_avg_bp.set_ylabel("Band Power (dB)", fontsize=fs_label, fontweight="bold")
        ax_avg_bp.set_title(f"{label} - {band_label} Hz Band Power (n={bp_mat.shape[0]})",
                            fontsize=fs_title, fontweight="bold", pad=8)
        ax_avg_bp.text(-0.08, 1.12, panel_labels[panel_idx], transform=ax_avg_bp.transAxes,
                       fontsize=fs_label + 8, fontweight='bold', va='top', ha='right')
        panel_idx += 1

        # Per-trial heatmap: fiber
        ax_heat_fiber = fig.add_subplot(gs[row_heat, 0])
        n_trials_disp = fiber_mat.shape[0]
        vmax_f = np.nanpercentile(np.abs(fiber_mat[np.isfinite(fiber_mat)]), 95) if np.any(np.isfinite(fiber_mat)) else 1.0
        im_f = ax_heat_fiber.imshow(
            fiber_mat, aspect='auto', origin='lower',
            extent=[t_grid[0], t_grid[-1], 0.5, n_trials_disp + 0.5],
            cmap=CMAP_PARULA_LIKE, vmin=-vmax_f, vmax=vmax_f, interpolation='nearest'
        )
        ax_heat_fiber.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
        ax_heat_fiber.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
        ax_heat_fiber.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
        ax_heat_fiber.set_ylabel("Trial", fontsize=fs_label, fontweight="bold")
        ax_heat_fiber.set_yticks(np.arange(1, n_trials_disp + 1))
        ax_heat_fiber.set_yticklabels([f"T{i}" for i in range(1, n_trials_disp + 1)])
        cbar_f = fig.colorbar(im_f, ax=ax_heat_fiber, pad=0.02, fraction=0.04)
        cbar_f.set_label("\u0394F/F (%)", fontsize=fs_tick)
        cbar_f.ax.tick_params(labelsize=fs_tick - 2)
        ax_heat_fiber.text(-0.08, 1.05, panel_labels[panel_idx], transform=ax_heat_fiber.transAxes,
                           fontsize=fs_label + 8, fontweight='bold', va='top', ha='right')
        panel_idx += 1

        # Per-trial heatmap: band power
        ax_heat_bp = fig.add_subplot(gs[row_heat, 1])
        n_trials_bp = bp_mat.shape[0]
        vmax_bp = np.nanpercentile(np.abs(bp_mat[np.isfinite(bp_mat)]), 95) if np.any(np.isfinite(bp_mat)) else 3.0
        im_bp = ax_heat_bp.imshow(
            bp_mat, aspect='auto', origin='lower',
            extent=[bin_centers[0], bin_centers[-1], 0.5, n_trials_bp + 0.5],
            cmap=CMAP_PARULA_LIKE, vmin=-vmax_bp, vmax=vmax_bp, interpolation='nearest'
        )
        ax_heat_bp.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
        ax_heat_bp.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
        ax_heat_bp.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
        ax_heat_bp.set_ylabel("Trial", fontsize=fs_label, fontweight="bold")
        ax_heat_bp.set_yticks(np.arange(1, n_trials_bp + 1))
        ax_heat_bp.set_yticklabels([f"T{i}" for i in range(1, n_trials_bp + 1)])
        cbar_bp = fig.colorbar(im_bp, ax=ax_heat_bp, pad=0.02, fraction=0.04)
        cbar_bp.set_label("Band Power (dB)", fontsize=fs_tick)
        cbar_bp.ax.tick_params(labelsize=fs_tick - 2)
        ax_heat_bp.text(-0.08, 1.05, panel_labels[panel_idx], transform=ax_heat_bp.transAxes,
                        fontsize=fs_label + 8, fontweight='bold', va='top', ha='right')
        panel_idx += 1

    # Style all axes
    for ax in fig.get_axes():
        if not hasattr(ax, 'images') or ax.__class__.__name__ != 'Axes':
            pass
        ax.tick_params(labelsize=fs_tick, width=TICK_WIDTH, length=TICK_LENGTH)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)

    fig.subplots_adjust(left=0.07, right=0.94, top=0.97, bottom=0.04)
    return fig


# =============================================================================
# FIGURE 5: 40Hz vs 135Hz COMPARISON (4 rows x 3 columns)
# =============================================================================

def _process_trials_for_comparison(trials, cfg):
    """Process trials and return aggregated data for comparison figure."""
    pre_sec = cfg['pre_sec']
    stim_sec = cfg['stim_sec']
    post_sec = cfg['post_sec']
    band = cfg['band']
    stim_freq = cfg['stim_freq']
    transient_end = 1.0
    sustained_start = 1.0

    x_lo = -min(pre_sec, 5.0)
    x_hi = stim_sec + min(post_sec, 5.0)
    t_grid = np.linspace(x_lo, x_hi, 2500)

    fiber_all = []
    period_fiber = {"Transient": [], "Sustained": [], "Post-stim": []}
    period_band_db = {"Transient": [], "Sustained": [], "Post-stim": []}
    spec_stack = []

    for tr in trials:
        fs = tr["fs"]
        sidx = int(np.argmin(np.abs(tr["t"] - (tr["t"][0] + pre_sec))))
        tt = tr["t"] - tr["t"][sidx]
        yy = tr["gevi"] * 100.0
        keep = (tt >= x_lo) & (tt <= x_hi)
        if np.sum(keep) < 100:
            continue
        tt = tt[keep]
        yy = yy[keep]
        yi = np.interp(t_grid, tt, yy, left=np.nan, right=np.nan)
        fiber_all.append(yi)

        bmask = (tt >= -pre_sec) & (tt < 0.0)
        baseline = np.nanmean(yy[bmask]) if np.any(bmask) else np.nan
        if np.isfinite(baseline):
            tmask = (tt >= 0.0) & (tt < transient_end)
            smask = (tt >= sustained_start) & (tt < stim_sec)
            pmask = (tt >= stim_sec) & (tt < (stim_sec + post_sec))
            if np.any(tmask):
                period_fiber["Transient"].append(float(np.nanmean(yy[tmask]) - baseline))
            if np.any(smask):
                period_fiber["Sustained"].append(float(np.nanmean(yy[smask]) - baseline))
            if np.any(pmask):
                period_fiber["Post-stim"].append(float(np.nanmean(yy[pmask]) - baseline))

        bpow = _compute_trial_band_power(yy[(tt >= -pre_sec) & (tt < 0.0)], fs, band)
        if np.isfinite(bpow) and bpow > 0:
            tpow = _compute_trial_band_power(yy[(tt >= 0.0) & (tt < transient_end)], fs, band)
            spow = _compute_trial_band_power(yy[(tt >= sustained_start) & (tt < stim_sec)], fs, band)
            ppow = _compute_trial_band_power(yy[(tt >= stim_sec) & (tt < (stim_sec + post_sec))], fs, band)
            if np.isfinite(tpow) and tpow > 0:
                period_band_db["Transient"].append(float(10 * np.log10(tpow / bpow)))
            if np.isfinite(spow) and spow > 0:
                period_band_db["Sustained"].append(float(10 * np.log10(spow / bpow)))
            if np.isfinite(ppow) and ppow > 0:
                period_band_db["Post-stim"].append(float(10 * np.log10(ppow / bpow)))

        nper = min(512, len(yy))
        if nper < 128:
            continue
        nover = min(int(0.875 * nper), nper - 1)
        nfft = max(2048, int(2 ** np.ceil(np.log2(nper))))
        f, ts, sxx = signal.spectrogram(
            np.nan_to_num(yy - np.nanmean(yy), nan=0.0),
            fs=fs, window="hann", nperseg=nper, noverlap=nover,
            nfft=nfft, detrend="linear", scaling="density", mode="psd",
        )
        ts_rel = ts + (tt[0] if len(tt) else 0.0)
        fmk = (f >= 1.0) & (f <= 150.0)
        if np.any(fmk):
            f2 = f[fmk]
            s2 = sxx[fmk, :]
            bmk = (ts_rel >= -pre_sec) & (ts_rel < -1.0)
            if np.sum(bmk) == 0:
                bmk = ts_rel < 0.0
            if np.sum(bmk) > 0:
                base_s = np.nanmean(s2[:, bmk], axis=1, keepdims=True)
                base_s = np.where(np.isfinite(base_s) & (base_s > 0), base_s, np.nan)
                s_rel = 10 * np.log10(np.maximum(s2, 1e-15) / np.maximum(base_s, 1e-15))
                if np.any(np.isfinite(s_rel)):
                    spec_stack.append((f2, ts_rel, s_rel))

    fiber_mat = np.asarray(fiber_all, dtype=float) if len(fiber_all) else np.empty((0, len(t_grid)))
    fiber_mu = np.nanmean(fiber_mat, axis=0) if fiber_mat.size else np.full(len(t_grid), np.nan)
    if fiber_mat.shape[0] > 1:
        fiber_sem = np.nanstd(fiber_mat, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(fiber_mat), axis=0))
    else:
        fiber_sem = np.full(len(t_grid), np.nan)

    spec_avg = None
    if len(spec_stack):
        f_ref = spec_stack[0][0]
        t_ref = spec_stack[0][1]
        mats = []
        for ff, tt_s, ss in spec_stack:
            if len(ff) != len(f_ref) or np.max(np.abs(ff - f_ref)) > 1e-9:
                continue
            arr = np.full((len(f_ref), len(t_ref)), np.nan)
            for i in range(len(f_ref)):
                vv = ss[i, :]
                mk = np.isfinite(tt_s) & np.isfinite(vv)
                if np.sum(mk) >= 2:
                    arr[i, :] = np.interp(t_ref, tt_s[mk], vv[mk], left=np.nan, right=np.nan)
            mats.append(arr)
        if len(mats):
            spec_avg = (f_ref, t_ref, np.nanmean(np.asarray(mats), axis=0))

    return {
        "t_grid": t_grid,
        "fiber_mu": fiber_mu,
        "fiber_sem": fiber_sem,
        "period_fiber": period_fiber,
        "period_band_db": period_band_db,
        "spec_avg": spec_avg,
        "n_trials": fiber_mat.shape[0],
    }


def _add_significance_bracket(ax, x1, x2, y, p_value, line_height=None, text_offset=None, fontsize=None):
    """Add significance bracket between two positions."""
    if line_height is None:
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        line_height = 0.015 * y_range
    if text_offset is None:
        y_range = ax.get_ylim()[1] - ax.get_ylim()[0]
        text_offset = 0.005 * y_range
    if fontsize is None:
        fontsize = FONT_SIZE_TICK

    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'

    ax.plot([x1, x2], [y, y], 'k-', linewidth=1.5, clip_on=False)
    ax.plot([x1, x1], [y - line_height, y], 'k-', linewidth=1.5, clip_on=False)
    ax.plot([x2, x2], [y - line_height, y], 'k-', linewidth=1.5, clip_on=False)
    ax.text((x1 + x2) / 2, y + text_offset, sig_text,
            ha='center', va='bottom', fontsize=fontsize, fontweight='bold')


def _perform_test(data1, data2, paired=False):
    """Perform t-test. Returns (p_value, test_name)."""
    d1 = np.array([v for v in data1 if v is not None and np.isfinite(v)])
    d2 = np.array([v for v in data2 if v is not None and np.isfinite(v)])
    if len(d1) < 3 or len(d2) < 3:
        return 1.0, 'insufficient'
    if paired:
        n = min(len(d1), len(d2))
        stat, p = stats.ttest_rel(d1[:n], d2[:n])
    else:
        stat, p = stats.ttest_ind(d1, d2)
    return float(p) if np.isfinite(p) else 1.0, 'paired t-test' if paired else 'unpaired t-test'


def _holm_correct(pvals):
    """Holm-Bonferroni correction."""
    n = len(pvals)
    if n <= 1:
        return list(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    corrected = [0.0] * n
    for rank, (orig_idx, p) in enumerate(indexed):
        corrected[orig_idx] = min(p * (n - rank), 1.0)
    for i in range(1, n):
        idx_cur = indexed[i][0]
        idx_prev = indexed[i - 1][0]
        corrected[idx_cur] = max(corrected[idx_cur], corrected[idx_prev])
    return corrected


def _plot_comparison_violin(ax, data_dict, colors_dict, ylabel, comparisons=None):
    """
    Plot publication-quality violin with individual dots, box, and significance brackets.
    No pre-stim category. Careful spacing to avoid overlaps.
    """
    positions = []
    data_list = []
    labels = []
    colors_list = []
    label_to_pos = {}

    for i, (label, values) in enumerate(data_dict.items()):
        valid = [v for v in values if v is not None and np.isfinite(v)]
        if len(valid) > 0:
            positions.append(i)
            data_list.append(np.array(valid))
            labels.append(label)
            colors_list.append(colors_dict.get(label, np.array([0.5, 0.5, 0.5])))
            label_to_pos[label] = i

    if len(data_list) == 0:
        ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes, fontsize=FONT_SIZE_LABEL)
        return

    parts = ax.violinplot(data_list, positions=positions, widths=0.7,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, color in zip(parts['bodies'], colors_list):
        pc.set_facecolor(color)
        pc.set_alpha(0.30)
        pc.set_edgecolor(np.clip(np.array(color) * 0.3, 0, 1))
        pc.set_linewidth(2.5)

    bp = ax.boxplot(data_list, positions=positions, widths=0.20, patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(np.clip(np.array(color) * 0.6 + 0.4, 0, 1))
        patch.set_edgecolor(np.clip(np.array(color) * 0.4, 0, 1))
        patch.set_linewidth(2.0)
        patch.set_alpha(0.9)
    for w in bp['whiskers']:
        w.set_color('black'); w.set_linewidth(1.5)
    for m in bp['medians']:
        m.set_color('black'); m.set_linewidth(2.5)
    for c in bp['caps']:
        c.set_color('black'); c.set_linewidth(1.5)

    rng = np.random.default_rng(42)
    for i, (pos, vals) in enumerate(zip(positions, data_list)):
        jitter = rng.uniform(-0.10, 0.10, len(vals))
        ax.scatter(pos + jitter, vals, color=colors_list[i], s=50, alpha=0.85,
                   zorder=10, edgecolors='black', linewidths=1.2)

    ax.axhline(0, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)

    if comparisons and len(comparisons) > 0:
        y_max = max(np.max(d) for d in data_list)
        y_min = min(np.min(d) for d in data_list)
        y_range = y_max - y_min if (y_max - y_min) != 0 else 1.0
        bracket_y_start = y_max + 0.10 * y_range
        bracket_spacing = 0.10 * y_range
        line_h = 0.02 * y_range
        txt_off = 0.008 * y_range

        raw_pvals = []
        valid_comps = []
        for (l1, l2, paired) in comparisons:
            if l1 in label_to_pos and l2 in label_to_pos:
                d1 = [v for v in data_dict[l1] if v is not None and np.isfinite(v)]
                d2 = [v for v in data_dict[l2] if v is not None and np.isfinite(v)]
                p, _ = _perform_test(d1, d2, paired=paired)
                raw_pvals.append(p)
                valid_comps.append((l1, l2, p))

        corrected = _holm_correct(raw_pvals) if len(raw_pvals) > 1 else raw_pvals

        for ci, (l1, l2, _) in enumerate(valid_comps):
            x1 = label_to_pos[l1]
            x2 = label_to_pos[l2]
            p_corr = corrected[ci]
            by = bracket_y_start + ci * bracket_spacing
            _add_significance_bracket(ax, x1, x2, by, p_corr,
                                      line_height=line_h, text_offset=txt_off,
                                      fontsize=FONT_SIZE_TICK)

        top = bracket_y_start + len(valid_comps) * bracket_spacing + 0.08 * y_range
        ax.set_ylim(ax.get_ylim()[0], top)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK, rotation=30, ha='right')
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(labelsize=FONT_SIZE_TICK, width=TICK_WIDTH, length=TICK_LENGTH)


# Frequency-based colors
COLOR_40HZ = np.array([0.15, 0.55, 0.55])
COLOR_40HZ_LIGHT = np.array([0.35, 0.70, 0.70])
COLOR_135HZ = np.array([0.95, 0.65, 0.45])
COLOR_135HZ_LIGHT = np.array([0.98, 0.80, 0.60])


def create_comparison_composite(trials_40, trials_135, cfg_40, cfg_135, comparison_label=""):
    """
    Create Figure 5: 40Hz vs 135Hz comparison.
    4 rows x 3 columns. Row 2 is shorter (averaged traces).

    Row 1: Example traces (col1=40Hz, col2=135Hz, col3=zoom split top/bottom)
    Row 2: Averaged fiber (col1=40Hz, col2=135Hz, col3=overlaid)
    Row 3: Violins fiber by period (col1=40Hz, col2=135Hz, col3=comparison transient/sustained)
    Row 4: Spectrograms (col1=40Hz, col2=135Hz, col3=band power violin comparison)
    """
    pre_sec = cfg_40['pre_sec']
    stim_sec = cfg_40['stim_sec']
    post_sec = cfg_40['post_sec']
    stim_freq_40 = cfg_40['stim_freq']
    stim_freq_135 = cfg_135['stim_freq']
    band_40 = cfg_40['band']
    band_135 = cfg_135['band']

    x_lo = -min(pre_sec, 5.0)
    x_hi = stim_sec + min(post_sec, 5.0)

    # Process both conditions
    res_40 = _process_trials_for_comparison(trials_40, cfg_40)
    res_135 = _process_trials_for_comparison(trials_135, cfg_135)

    # Get representative trials
    rep_40 = trials_40[0] if len(trials_40) > 0 else None
    rep_135 = trials_135[0] if len(trials_135) > 0 else None

    # Figure layout: 4 rows, 3 columns; row 2 shorter
    fig = plt.figure(figsize=(36, 40))
    gs = GridSpec(4, 3, figure=fig,
                  height_ratios=[1.2, 0.7, 1.2, 1.2],
                  width_ratios=[1.0, 1.0, 1.0],
                  hspace=0.40, wspace=0.32)

    fs_label = FONT_SIZE_LABEL + 8
    fs_tick = FONT_SIZE_TICK + 6
    fs_title = FONT_SIZE_TITLE + 4
    panel_idx = 0
    panel_labels = "ABCDEFGHIJKL"

    def _style_ax(ax):
        ax.tick_params(labelsize=fs_tick, width=TICK_WIDTH, length=TICK_LENGTH)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)

    def _add_panel_label(ax, label):
        ax.text(-0.10, 1.08, label, transform=ax.transAxes,
                fontsize=fs_label + 10, fontweight='bold', va='top', ha='right')

    # =========================================================================
    # ROW 1: Example traces
    # =========================================================================
    def _plot_example_traces(ax, trial, cfg_single, title_str):
        """Plot stim, LFP, fiber, motion for one example trial."""
        if trial is None:
            ax.text(0.5, 0.5, "No data", ha='center', va='center',
                    transform=ax.transAxes, fontsize=fs_label)
            return
        p = cfg_single['pre_sec']
        s = cfg_single['stim_sec']
        sidx = int(np.argmin(np.abs(trial["t"] - (trial["t"][0] + p))))
        t_rel = trial["t"] - trial["t"][sidx]
        fiber_pct = trial["gevi"] * 100.0
        lfp = trial["lfp"]
        motion = trial["motion"]
        stim = generate_stim_pulses(t_rel, cfg_single['stim_freq'], s)

        m = (t_rel >= x_lo) & (t_rel <= x_hi)
        t_p = t_rel[m]

        def _norm(y):
            p95 = np.nanpercentile(np.abs(y), 95)
            return y / (p95 + 1e-9) if p95 > 0 else y

        y_stim = _norm(stim[m]) * 0.5
        y_lfp = _norm(lfp[m] - np.nanmedian(lfp[m])) * 0.8
        y_fiber = _norm(fiber_pct[m] - np.nanmedian(fiber_pct[m])) * 0.8
        y_motion = _norm(motion[m] - np.nanmedian(motion[m])) * 0.7

        offsets = [4.5, 1.5, -1.5, -4.5]
        ax.plot(t_p, y_stim + offsets[0], color=COLOR_STIM_PULSE, linewidth=1.3)
        ax.plot(t_p, y_lfp + offsets[1], color=COLOR_LFP, linewidth=1.3)
        ax.plot(t_p, y_fiber + offsets[2], color=COLOR_GEVI, linewidth=1.8)
        motion_dark = np.clip(np.asarray(COLOR_MOTION) * 0.58, 0.0, 1.0)
        ax.plot(t_p, y_motion + offsets[3], color=motion_dark, linewidth=1.3)

        ax.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
        ax.axvline(s, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(-6.5, 6.5)
        ax.set_yticks(offsets)
        ax.set_yticklabels(["Stim", "LFP", "Fiber", "Motion"], fontsize=fs_tick)
        ax.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
        ax.set_title(title_str, fontsize=fs_title, fontweight='bold', pad=10)

    # Col 1: 40Hz example
    ax_r1c1 = fig.add_subplot(gs[0, 0])
    _plot_example_traces(ax_r1c1, rep_40, cfg_40, f"{cfg_40['label']} - Example Trial")
    _style_ax(ax_r1c1)
    _add_panel_label(ax_r1c1, panel_labels[panel_idx]); panel_idx += 1

    # Col 2: 135Hz example
    ax_r1c2 = fig.add_subplot(gs[0, 1])
    _plot_example_traces(ax_r1c2, rep_135, cfg_135, f"{cfg_135['label']} - Example Trial")
    _style_ax(ax_r1c2)
    _add_panel_label(ax_r1c2, panel_labels[panel_idx]); panel_idx += 1

    # Col 3: Zoom-in (split into two sub-rows via inner gridspec)
    ax_r1c3 = fig.add_subplot(gs[0, 2])
    ax_r1c3.set_visible(False)
    gs_zoom = gs[0, 2].subgridspec(2, 1, hspace=0.4)
    zoom_pre_s = 0.25
    zoom_post_s = 0.40

    def _plot_zoom(ax_z, trial, cfg_single, color_stim, label_str):
        if trial is None:
            ax_z.text(0.5, 0.5, "No data", ha='center', va='center',
                      transform=ax_z.transAxes, fontsize=fs_tick)
            return
        p = cfg_single['pre_sec']
        s = cfg_single['stim_sec']
        sidx = int(np.argmin(np.abs(trial["t"] - (trial["t"][0] + p))))
        t_rel = trial["t"] - trial["t"][sidx]
        fiber_pct = trial["gevi"] * 100.0
        stim = generate_stim_pulses(t_rel, cfg_single['stim_freq'], s)
        zm = (t_rel >= -zoom_pre_s) & (t_rel <= zoom_post_s)
        t_z = t_rel[zm]
        y_f = fiber_pct[zm] - np.nanmedian(fiber_pct[zm])
        y_s = stim[zm]
        nf = y_f / (np.nanpercentile(np.abs(y_f), 95) + 1e-9)
        ns = y_s / (np.nanpercentile(np.abs(y_s), 95) + 1e-9) if np.any(y_s != 0) else np.zeros_like(t_z)
        ax_z.plot(t_z, ns * 0.4 + 1.0, color=color_stim, linewidth=1.5)
        ax_z.plot(t_z, nf * 0.8 - 0.6, color=COLOR_GEVI, linewidth=2.2)
        ax_z.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
        ax_z.set_xlim(-zoom_pre_s, zoom_post_s)
        ax_z.set_ylim(-2.0, 2.0)
        ax_z.set_yticks([1.0, -0.6])
        ax_z.set_yticklabels(["Stim", "Fiber"], fontsize=fs_tick - 2)
        ax_z.set_title(label_str, fontsize=fs_tick, fontweight='bold', pad=4)
        _style_ax(ax_z)

    ax_zoom_40 = fig.add_subplot(gs_zoom[0])
    _plot_zoom(ax_zoom_40, rep_40, cfg_40, COLOR_STIM_PULSE, f"{cfg_40['label']} Zoom")
    ax_zoom_135 = fig.add_subplot(gs_zoom[1])
    _plot_zoom(ax_zoom_135, rep_135, cfg_135, COLOR_STIM_PULSE, f"{cfg_135['label']} Zoom")
    ax_zoom_135.set_xlabel("Time from stim onset (s)", fontsize=fs_label - 2, fontweight="bold")
    _add_panel_label(ax_zoom_40, panel_labels[panel_idx]); panel_idx += 1

    # =========================================================================
    # ROW 2: Averaged fiber signals
    # =========================================================================
    t_grid_40 = res_40["t_grid"]
    t_grid_135 = res_135["t_grid"]

    # Col 1: 40Hz averaged
    ax_r2c1 = fig.add_subplot(gs[1, 0])
    ax_r2c1.plot(t_grid_40, res_40["fiber_mu"], color=COLOR_40HZ, linewidth=2.0)
    ax_r2c1.fill_between(t_grid_40, res_40["fiber_mu"] - res_40["fiber_sem"],
                         res_40["fiber_mu"] + res_40["fiber_sem"],
                         color=COLOR_40HZ, alpha=0.25, linewidth=0)
    ax_r2c1.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
    ax_r2c1.axvline(stim_sec, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
    ax_r2c1.set_xlim(x_lo, x_hi)
    ax_r2c1.set_ylabel("\u0394F/F (%)", fontsize=fs_label, fontweight="bold")
    ax_r2c1.set_title(f"{cfg_40['label']} Averaged (n={res_40['n_trials']})",
                      fontsize=fs_title - 2, fontweight='bold', pad=6)
    _style_ax(ax_r2c1)
    _add_panel_label(ax_r2c1, panel_labels[panel_idx]); panel_idx += 1

    # Col 2: 135Hz averaged
    ax_r2c2 = fig.add_subplot(gs[1, 1])
    ax_r2c2.plot(t_grid_135, res_135["fiber_mu"], color=COLOR_135HZ, linewidth=2.0)
    ax_r2c2.fill_between(t_grid_135, res_135["fiber_mu"] - res_135["fiber_sem"],
                         res_135["fiber_mu"] + res_135["fiber_sem"],
                         color=COLOR_135HZ, alpha=0.25, linewidth=0)
    ax_r2c2.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
    ax_r2c2.axvline(stim_sec, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
    ax_r2c2.set_xlim(x_lo, x_hi)
    ax_r2c2.set_ylabel("\u0394F/F (%)", fontsize=fs_label, fontweight="bold")
    ax_r2c2.set_title(f"{cfg_135['label']} Averaged (n={res_135['n_trials']})",
                      fontsize=fs_title - 2, fontweight='bold', pad=6)
    _style_ax(ax_r2c2)
    _add_panel_label(ax_r2c2, panel_labels[panel_idx]); panel_idx += 1

    # Col 3: Overlaid
    ax_r2c3 = fig.add_subplot(gs[1, 2])
    ax_r2c3.plot(t_grid_40, res_40["fiber_mu"], color=COLOR_40HZ, linewidth=2.0,
                 label=f"{cfg_40['label']} (n={res_40['n_trials']})")
    ax_r2c3.fill_between(t_grid_40, res_40["fiber_mu"] - res_40["fiber_sem"],
                         res_40["fiber_mu"] + res_40["fiber_sem"],
                         color=COLOR_40HZ, alpha=0.20, linewidth=0)
    ax_r2c3.plot(t_grid_135, res_135["fiber_mu"], color=COLOR_135HZ, linewidth=2.0,
                 label=f"{cfg_135['label']} (n={res_135['n_trials']})")
    ax_r2c3.fill_between(t_grid_135, res_135["fiber_mu"] - res_135["fiber_sem"],
                         res_135["fiber_mu"] + res_135["fiber_sem"],
                         color=COLOR_135HZ, alpha=0.20, linewidth=0)
    ax_r2c3.axvline(0.0, color="k", linestyle="--", linewidth=1.0, alpha=0.7)
    ax_r2c3.axvline(stim_sec, color="k", linestyle=":", linewidth=1.0, alpha=0.7)
    ax_r2c3.set_xlim(x_lo, x_hi)
    ax_r2c3.set_ylabel("\u0394F/F (%)", fontsize=fs_label, fontweight="bold")
    ax_r2c3.set_title("Overlaid Comparison", fontsize=fs_title - 2, fontweight='bold', pad=6)
    ax_r2c3.legend(fontsize=fs_tick - 2, loc='upper right', framealpha=0.8)
    _style_ax(ax_r2c3)
    _add_panel_label(ax_r2c3, panel_labels[panel_idx]); panel_idx += 1

    # =========================================================================
    # ROW 3: Violin plots for fiber signal
    # =========================================================================
    # Col 1: 40Hz by period (Transient, Sustained, Post-stim only)
    ax_r3c1 = fig.add_subplot(gs[2, 0])
    violin_data_40 = {
        "Transient": res_40["period_fiber"]["Transient"],
        "Sustained": res_40["period_fiber"]["Sustained"],
        "Post-stim": res_40["period_fiber"]["Post-stim"],
    }
    violin_colors_40 = {
        "Transient": COLOR_TRANSIENT,
        "Sustained": COLOR_SUSTAINED,
        "Post-stim": COLOR_POST,
    }
    period_comps = [
        ("Transient", "Sustained", True),
        ("Sustained", "Post-stim", True),
    ]
    _plot_comparison_violin(ax_r3c1, violin_data_40, violin_colors_40,
                            "Fiber \u0394F/F change (%)", comparisons=period_comps)
    ax_r3c1.set_title(f"{cfg_40['label']} - Fiber by Period", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _add_panel_label(ax_r3c1, panel_labels[panel_idx]); panel_idx += 1

    # Col 2: 135Hz by period
    ax_r3c2 = fig.add_subplot(gs[2, 1])
    violin_data_135 = {
        "Transient": res_135["period_fiber"]["Transient"],
        "Sustained": res_135["period_fiber"]["Sustained"],
        "Post-stim": res_135["period_fiber"]["Post-stim"],
    }
    violin_colors_135 = {
        "Transient": COLOR_TRANSIENT,
        "Sustained": COLOR_SUSTAINED,
        "Post-stim": COLOR_POST,
    }
    _plot_comparison_violin(ax_r3c2, violin_data_135, violin_colors_135,
                            "Fiber \u0394F/F change (%)", comparisons=period_comps)
    ax_r3c2.set_title(f"{cfg_135['label']} - Fiber by Period", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _add_panel_label(ax_r3c2, panel_labels[panel_idx]); panel_idx += 1

    # Col 3: Cross-condition comparison (Transient and Sustained for 40 vs 135)
    ax_r3c3 = fig.add_subplot(gs[2, 2])
    comp_fiber_data = {
        "40Hz Trans": res_40["period_fiber"]["Transient"],
        "40Hz Sust": res_40["period_fiber"]["Sustained"],
        "135Hz Trans": res_135["period_fiber"]["Transient"],
        "135Hz Sust": res_135["period_fiber"]["Sustained"],
    }
    comp_fiber_colors = {
        "40Hz Trans": COLOR_40HZ,
        "40Hz Sust": COLOR_40HZ_LIGHT,
        "135Hz Trans": COLOR_135HZ,
        "135Hz Sust": COLOR_135HZ_LIGHT,
    }
    cross_comps = [
        ("40Hz Trans", "40Hz Sust", True),
        ("135Hz Trans", "135Hz Sust", True),
        ("40Hz Trans", "135Hz Trans", False),
        ("40Hz Sust", "135Hz Sust", False),
    ]
    _plot_comparison_violin(ax_r3c3, comp_fiber_data, comp_fiber_colors,
                            "Fiber \u0394F/F change (%)", comparisons=cross_comps)
    ax_r3c3.set_title("Fiber: 40Hz vs 135Hz", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _add_panel_label(ax_r3c3, panel_labels[panel_idx]); panel_idx += 1

    # =========================================================================
    # ROW 4: Spectrograms + band power violin
    # =========================================================================
    # Col 1: 40Hz spectrogram
    ax_r4c1 = fig.add_subplot(gs[3, 0])
    if res_40["spec_avg"] is not None:
        ff, tt_s, ss = res_40["spec_avg"]
        ss_d = uniform_filter1d(uniform_filter1d(ss, size=3, axis=0, mode="nearest"),
                                size=5, axis=1, mode="nearest")
        vmax = np.nanpercentile(np.abs(ss_d[np.isfinite(ss_d)]), 95) if np.any(np.isfinite(ss_d)) else 2.0
        im1 = ax_r4c1.pcolormesh(tt_s, ff, ss_d, cmap=CMAP_PARULA_LIKE, shading="auto",
                                  vmin=-vmax, vmax=vmax)
        ax_r4c1.set_ylim(1, 100)
        cbar1 = fig.colorbar(im1, ax=ax_r4c1, pad=0.02, fraction=0.04)
        cbar1.set_label("Power (dB)", fontsize=fs_tick)
        cbar1.ax.tick_params(labelsize=fs_tick - 2)
    else:
        ax_r4c1.text(0.5, 0.5, "No spectrogram", ha='center', va='center',
                     transform=ax_r4c1.transAxes, fontsize=fs_label)
    ax_r4c1.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
    ax_r4c1.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
    ax_r4c1.set_xlim(x_lo, x_hi)
    ax_r4c1.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    ax_r4c1.set_ylabel("Frequency (Hz)", fontsize=fs_label, fontweight="bold")
    ax_r4c1.set_title(f"{cfg_40['label']} Spectrogram", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _style_ax(ax_r4c1)
    _add_panel_label(ax_r4c1, panel_labels[panel_idx]); panel_idx += 1

    # Col 2: 135Hz spectrogram
    ax_r4c2 = fig.add_subplot(gs[3, 1])
    if res_135["spec_avg"] is not None:
        ff, tt_s, ss = res_135["spec_avg"]
        ss_d = uniform_filter1d(uniform_filter1d(ss, size=3, axis=0, mode="nearest"),
                                size=5, axis=1, mode="nearest")
        vmax = np.nanpercentile(np.abs(ss_d[np.isfinite(ss_d)]), 95) if np.any(np.isfinite(ss_d)) else 2.0
        im2 = ax_r4c2.pcolormesh(tt_s, ff, ss_d, cmap=CMAP_PARULA_LIKE, shading="auto",
                                  vmin=-vmax, vmax=vmax)
        ax_r4c2.set_ylim(1, 100)
        cbar2 = fig.colorbar(im2, ax=ax_r4c2, pad=0.02, fraction=0.04)
        cbar2.set_label("Power (dB)", fontsize=fs_tick)
        cbar2.ax.tick_params(labelsize=fs_tick - 2)
    else:
        ax_r4c2.text(0.5, 0.5, "No spectrogram", ha='center', va='center',
                     transform=ax_r4c2.transAxes, fontsize=fs_label)
    ax_r4c2.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
    ax_r4c2.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
    ax_r4c2.set_xlim(x_lo, x_hi)
    ax_r4c2.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    ax_r4c2.set_ylabel("Frequency (Hz)", fontsize=fs_label, fontweight="bold")
    ax_r4c2.set_title(f"{cfg_135['label']} Spectrogram", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _style_ax(ax_r4c2)
    _add_panel_label(ax_r4c2, panel_labels[panel_idx]); panel_idx += 1

    # Col 3: Band power violin (cross-condition comparison)
    ax_r4c3 = fig.add_subplot(gs[3, 2])
    comp_bp_data = {
        "40Hz Trans": res_40["period_band_db"]["Transient"],
        "40Hz Sust": res_40["period_band_db"]["Sustained"],
        "135Hz Trans": res_135["period_band_db"]["Transient"],
        "135Hz Sust": res_135["period_band_db"]["Sustained"],
    }
    comp_bp_colors = {
        "40Hz Trans": COLOR_40HZ,
        "40Hz Sust": COLOR_40HZ_LIGHT,
        "135Hz Trans": COLOR_135HZ,
        "135Hz Sust": COLOR_135HZ_LIGHT,
    }
    _plot_comparison_violin(ax_r4c3, comp_bp_data, comp_bp_colors,
                            "Relative band power (dB)", comparisons=cross_comps)
    ax_r4c3.set_title("Stim-Band Power: 40Hz vs 135Hz", fontsize=fs_title - 2, fontweight='bold', pad=8)
    _add_panel_label(ax_r4c3, panel_labels[panel_idx]); panel_idx += 1

    fig.subplots_adjust(left=0.06, right=0.96, top=0.96, bottom=0.04)
    return fig


# =============================================================================
# MAIN
# =============================================================================

def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 70)
    print("PUBLICATION COMPOSITE FIGURES")
    print(f"Output: {OUTPUT_DIR}")
    print("=" * 70)

    mice = ['Animal03', 'Animal04']
    mouse_suffix = {'Animal03': 'a', 'Animal04': 'b'}
    conditions_order = ['135Hz', '40Hz_Amp', '40Hz_Energy']
    fig_titles = {
        '135Hz': "135 Hz DBS",
        '40Hz_Amp': "40 Hz DBS (Amp-balanced)",
        '40Hz_Energy': "40 Hz DBS (Energy-balanced)",
    }

    for mouse in mice:
        suffix = mouse_suffix[mouse]
        print(f"\n{'='*70}")
        print(f"  MOUSE: {mouse} (Figure suffix: {suffix})")
        print(f"{'='*70}")

        all_trials = {}
        session_cfgs = SESSIONS[mouse]

        for cond_key in conditions_order:
            cfg = session_cfgs[cond_key]
            print(f"\n  Loading {cond_key} ({cfg['session_id']})...")
            trials = load_session_trials(
                cfg['base_path'], cfg['session_id'], mouse, num_trials=NUM_TRIALS
            )
            all_trials[cond_key] = trials
            print(f"    -> {len(trials)} trials loaded")

        # Figures 1-3: 6-panel composites per condition
        for fig_num, cond_key in enumerate(conditions_order, start=1):
            cfg = session_cfgs[cond_key]
            trials = all_trials[cond_key]
            print(f"\n  Creating Figure {fig_num}{suffix}: {fig_titles[cond_key]}...")

            fig = create_six_panel_composite(trials, cfg)
            fig.suptitle(f"Figure {fig_num}{suffix}: {fig_titles[cond_key]}",
                         fontsize=32, fontweight='bold', y=0.998)

            fname = f"F{fig_num}{suffix}"
            fig.savefig(str(OUTPUT_DIR / f"{fname}.pdf"), dpi=DPI, bbox_inches="tight")
            fig.savefig(str(OUTPUT_DIR / f"{fname}.png"), dpi=DPI, bbox_inches="tight")
            plt.close(fig)
            print(f"    Saved: {fname}.pdf/.png")

        # Figure 4: All-conditions composite
        print(f"\n  Creating Figure 4{suffix}: All conditions composite...")
        fig4 = create_all_conditions_composite(all_trials, session_cfgs, mouse_label=mouse)
        fig4.suptitle(f"Figure 4{suffix}: All DBS Conditions - Trial Heatmaps",
                      fontsize=28, fontweight='bold', y=0.998)

        fname4 = f"F4{suffix}"
        fig4.savefig(str(OUTPUT_DIR / f"{fname4}.pdf"), dpi=DPI, bbox_inches="tight")
        fig4.savefig(str(OUTPUT_DIR / f"{fname4}.png"), dpi=DPI, bbox_inches="tight")
        plt.close(fig4)
        print(f"    Saved: {fname4}.pdf/.png")

        # Figure 5: 40Hz vs 135Hz comparison (one per balancing condition)
        comparison_pairs = [
            ('40Hz_Amp', '135Hz', 'amp', '40 Hz Amp-balanced vs 135 Hz'),
            ('40Hz_Energy', '135Hz', 'nrg', '40 Hz Energy-balanced vs 135 Hz'),
        ]
        for cond_40_key, cond_135_key, comp_tag, comp_title in comparison_pairs:
            print(f"\n  Creating Figure 5{suffix}_{comp_tag}: {comp_title}...")
            fig5 = create_comparison_composite(
                all_trials[cond_40_key], all_trials[cond_135_key],
                session_cfgs[cond_40_key], session_cfgs[cond_135_key],
                comparison_label=comp_title
            )
            fig5.suptitle(f"Figure 5{suffix}: {comp_title}",
                          fontsize=28, fontweight='bold', y=0.995)

            fname5 = f"F5{suffix}_{comp_tag}"
            fig5.savefig(str(OUTPUT_DIR / f"{fname5}.pdf"), dpi=DPI, bbox_inches="tight")
            fig5.savefig(str(OUTPUT_DIR / f"{fname5}.png"), dpi=DPI, bbox_inches="tight")
            plt.close(fig5)
            print(f"    Saved: {fname5}.pdf/.png")

    print("\n" + "=" * 70)
    print("DONE! All publication composite figures generated.")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
