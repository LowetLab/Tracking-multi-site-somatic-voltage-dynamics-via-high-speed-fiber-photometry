"""
Plot Cellular DBS Comparison Figures (40 Hz vs 135 Hz).

Generates 8 figures for a pair of cellular-resolution DBS sessions:
1.  Trial-averaged traces — 40 Hz (neuron-avg → trial-avg: LFP + voltage ± SEM)
2.  Trial-averaged traces — 135 Hz (same layout)
3.  Overlaid trial-averaged LFP and voltage on shared axes (40 vs 135)
4.  Period violin: voltage signal by period — 40 Hz (pre, trans, sust, post)
5.  Period violin: voltage signal by period — 135 Hz
6.  40 vs 135 Hz comparison violin (trans/sust)
7.  Trial-averaged time–frequency spectrograms (LFP + neuron 1, baseline-normalised)
    — one figure per condition (40 Hz and 135 Hz), so 7a and 7b
8.  Stim-band relative power violin (trans vs sust, 40 vs 135)

Data comes from CellularAnalysis .mat files (HDF5 v7.3).
Styling follows stimulation_analysis.py (publication standard).
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy import signal, stats
from scipy.ndimage import uniform_filter1d
from pathlib import Path
import warnings
import sys
import h5py

warnings.filterwarnings("ignore")

# common.py lives in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# =============================================================================
# CONFIGURATION  -- EDIT THESE FOR YOUR OWN DATASET
# =============================================================================

BASE_PATH = _LAB_DATA_ROOT / "DBS" / "Animal01" / "BaselineDBS" / "CellularDataProcessed"
OUTPUT_DIR = PROJECT_ROOT / "Figures" / "Cellular_DBS_comparison"
MOUSE_NAME = "Animal01"

# Sessions to compare
SESSION_R1 = {"date": "04-06-25", "rec_id": "R13", "freq_hz": 40,
              "label": "40Hz (Energy-balanced)"}
SESSION_R2 = {"date": "05-06-25", "rec_id": "R10", "freq_hz": 135,
              "label": "130Hz (Energy-balanced)"}

# Time parameters (relative to stim onset = 0)
PRE_STIM_SEC = 3.0        # 3 s before stim onset
STIM_DURATION_SEC = 1.0   # 1 s stimulation
POST_STIM_SEC = 2.0       # 2 s after stim offset
TRANSIENT_END = 0.15      # first 150 ms of stim
SUSTAINED_START = 0.15    # rest of stim after transient

# Spectrogram baseline window [LO, HI) seconds relative to stim onset
# Matches fiber pipeline (stimulation_analysis.py): [-3, -1) avoids last second
# before onset and never includes stimulation.  Falls back to time < 0 if no bins match.
SPEC_BASELINE_T_LO = -3.0
SPEC_BASELINE_T_HI = -1.0

# Neuron index for spectrogram panels (0-based)
NEURON_IDX_SPEC = 0

# Stim-band frequency ranges for band-power violin
BAND_40HZ = (35, 45)
BAND_130HZ = (125, 135)

# Spectrogram computation — matched to MATLAB fiber pipeline (stim_analysis_config.m)
# Window 0.6 s, 95% overlap for dense time sampling, NFFT = nextpow2(3× window) for
# fine frequency grid.  shading="gouraud" in pcolormesh interpolates between bins.
SPEC_WINDOW_SEC = 0.6
SPEC_OVERLAP_FRAC = 0.95     # 95% → time step ≈ 30 ms → ~100 columns in 3 s
SPEC_NFFT_MULT = 3           # zero-pad to 3× window before rounding to next power of 2
SPEC_SMOOTH_FREQ = 3         # post-hoc smoothing in frequency bins
SPEC_SMOOTH_TIME = 3         # post-hoc smoothing in time bins
FREQ_RANGE = (1, 150)

# Spectrogram color scaling (shared across conditions, per modality)
# Use wider robust limits to avoid saturation of stimulation epochs.
SPEC_LFP_PCTL = (0.5, 99.7)
SPEC_NEUR_PCTL = (1.0, 99.0)

# Figure configuration
DPI = 300
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 22
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_SCALEBAR = 13

AXIS_LINEWIDTH = 2.5
TICK_WIDTH = 2.0
TICK_LENGTH = 8
LINE_WIDTH_TRACE = 1.2
LINE_WIDTH_THICK = 2.0

# Colors (consistent with fiber pipeline)
COLOR_LFP = np.array([0.35, 0.25, 0.45])
COLOR_NEURON = np.array([0.08, 0.45, 0.45])
COLOR_STIM_PULSE = np.array([0.4, 0.1, 0.1])

COLOR_PRE = np.array([0.5, 0.5, 0.5])
COLOR_TRANSIENT = np.array([0.7, 0.15, 0.15])
COLOR_SUSTAINED = np.array([0.65, 0.25, 0.15])
COLOR_POST = np.array([0.15, 0.55, 0.55])

COLOR_40HZ = np.array([0.15, 0.55, 0.55])
COLOR_40HZ_LIGHT = np.array([0.35, 0.70, 0.70])
COLOR_135HZ = np.array([0.95, 0.65, 0.45])
COLOR_135HZ_LIGHT = np.array([0.98, 0.80, 0.60])


# =============================================================================
# DATA LOADING  (reused from cellular_dbs_traces.py)
# =============================================================================

def load_cellular_data(session):
    """Load CellularAnalysis .mat (HDF5 v7.3) for one session."""
    date = session["date"]
    rec_id = session["rec_id"]
    session_folder = BASE_PATH / f"{date}-{rec_id}"
    mat_file = session_folder / f"{MOUSE_NAME}_{date}-{rec_id}_CellularAnalysis.mat"

    if not mat_file.exists():
        print(f"  WARNING: File not found: {mat_file}")
        return None

    print(f"  Loading: {mat_file}")

    try:
        with h5py.File(str(mat_file), "r") as f:
            cellular = f["CellularAnalysis"]
            num_trials = int(np.array(cellular["metadata"]["num_trials"]).flat[0])
            num_neurons = int(np.array(cellular["metadata"]["num_neurons"]).flat[0])
            print(f"    Found {num_trials} trials, {num_neurons} neurons")

            trials_refs = cellular["trials"][()]
            trials_data = []
            for trial_idx in range(num_trials):
                trial_grp = f[trials_refs.flat[trial_idx]]
                td = {}

                # Time vector
                try:
                    time_grp = trial_grp["time"]
                    tv_ref = time_grp["time_vector"][()]
                    if isinstance(tv_ref, np.ndarray) and tv_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        td["time_vector"] = f[tv_ref.flat[0]][()].flatten()
                    else:
                        td["time_vector"] = np.array(tv_ref).flatten()
                except Exception:
                    td["time_vector"] = np.array([])

                # Stimulus onset frame
                try:
                    stim_onset = time_grp["stimulus_onset_frame"][()]
                    if isinstance(stim_onset, np.ndarray) and stim_onset.dtype == h5py.special_dtype(ref=h5py.Reference):
                        td["stim_onset_frame"] = int(f[stim_onset.flat[0]][()].flat[0])
                    else:
                        td["stim_onset_frame"] = int(np.array(stim_onset).flat[0])
                except Exception:
                    td["stim_onset_frame"] = None

                # Sampling rate
                try:
                    params_grp = trial_grp["parameters"]
                    fs_data = params_grp["imaging_fs"][()]
                    if isinstance(fs_data, np.ndarray) and fs_data.dtype == h5py.special_dtype(ref=h5py.Reference):
                        td["fs"] = float(f[fs_data.flat[0]][()].flat[0])
                    else:
                        td["fs"] = float(np.array(fs_data).flat[0])
                except Exception:
                    td["fs"] = 1000.0

                # Fluorescence (frames × neurons)
                try:
                    signals_grp = trial_grp["signals"]
                    fluor_ref = signals_grp["fluorescence_corrected"][()]
                    if isinstance(fluor_ref, np.ndarray) and fluor_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        fluor_data = f[fluor_ref.flat[0]][()]
                    else:
                        fluor_data = np.array(fluor_ref)
                    if fluor_data.ndim == 2 and fluor_data.shape[0] < fluor_data.shape[1]:
                        fluor_data = fluor_data.T
                    td["fluorescence"] = fluor_data
                except Exception:
                    td["fluorescence"] = np.array([])

                # LFP
                try:
                    ephys_grp = trial_grp["ephys"]
                    lfp_ref = ephys_grp["lfp_raw_aligned"][()]
                    if isinstance(lfp_ref, np.ndarray) and lfp_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        td["lfp"] = f[lfp_ref.flat[0]][()].flatten()
                    else:
                        td["lfp"] = np.array(lfp_ref).flatten()
                except Exception:
                    td["lfp"] = None

                trials_data.append(td)
                print(f"      Trial {trial_idx+1}: {len(td['time_vector'])} samples, fs={td['fs']:.1f}Hz")

            return {
                "trials": trials_data,
                "num_trials": num_trials,
                "num_neurons": num_neurons,
                "session": session,
            }
    except Exception as e:
        print(f"    Error loading: {e}")
        import traceback; traceback.print_exc()
        return None


# =============================================================================
# HELPERS
# =============================================================================

def _window_mask(t, t_min, t_max):
    return (t >= t_min) & (t <= t_max)


def _collect_trial_traces(data):
    """
    Build aligned arrays across trials for the plotting window.

    Returns (t_common, lfp_all, voltage_all, fs)
        lfp_all     : (n_trials, n_samples)  — may have NaN rows if LFP missing
        voltage_all : (n_trials, n_samples)  — neuron-averaged per trial
    """
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    ref = data["trials"][0]
    fs = ref["fs"]
    mask0 = _window_mask(ref["time_vector"], t_min, t_max)
    t_common = ref["time_vector"][mask0]
    n_samp = len(t_common)

    lfp_list, volt_list = [], []
    for trial in data["trials"]:
        t_tr = trial["time_vector"]
        m = _window_mask(t_tr, t_min, t_max)
        if np.sum(m) != n_samp:
            continue
        if trial["lfp"] is not None:
            lfp_list.append(trial["lfp"][m])
        else:
            lfp_list.append(np.full(n_samp, np.nan))
        fl = trial["fluorescence"]
        if fl.ndim == 2 and fl.shape[0] >= n_samp:
            volt_list.append(np.mean(fl[m, :], axis=1))
        else:
            volt_list.append(np.full(n_samp, np.nan))

    return t_common, np.array(lfp_list), np.array(volt_list), fs


def _collect_neuron_trace(data, neuron_idx):
    """Single-neuron fluorescence across trials aligned to common grid."""
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    ref = data["trials"][0]
    mask0 = _window_mask(ref["time_vector"], t_min, t_max)
    t_common = ref["time_vector"][mask0]
    n_samp = len(t_common)
    traces = []
    for trial in data["trials"]:
        t_tr = trial["time_vector"]
        m = _window_mask(t_tr, t_min, t_max)
        if np.sum(m) != n_samp:
            continue
        fl = trial["fluorescence"]
        if fl.ndim == 2 and fl.shape[0] >= n_samp and neuron_idx < fl.shape[1]:
            traces.append(fl[m, neuron_idx])
        else:
            traces.append(np.full(n_samp, np.nan))
    return t_common, np.array(traces)


def _smooth(arr, w=5):
    return uniform_filter1d(arr, size=w)


def _mean_sem(arr2d, axis=0):
    m = np.nanmean(arr2d, axis=axis)
    s = np.nanstd(arr2d, axis=axis, ddof=1) / np.sqrt(np.sum(np.isfinite(arr2d), axis=axis).clip(1))
    return m, s


def _period_values(t, traces_2d):
    """
    Per-trial mean over 4 periods.
    traces_2d: (n_trials, n_time).
    Returns dict of lists per period.
    """
    masks = {
        "Pre-stim": t < 0,
        "Transient": (t >= 0) & (t < TRANSIENT_END),
        "Sustained": (t >= SUSTAINED_START) & (t < STIM_DURATION_SEC),
        "Post-stim": t >= STIM_DURATION_SEC,
    }
    out = {}
    for label, m in masks.items():
        vals = [float(np.nanmean(row[m])) if np.any(m) else np.nan for row in traces_2d]
        out[label] = vals
    return out


def _style_avg_axis(ax, show_x=True):
    ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
    ax.axvline(STIM_DURATION_SEC, color="black", linestyle=":", linewidth=1.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis="both", which="major", labelsize=FONT_SIZE_TICK,
                   width=TICK_WIDTH, length=TICK_LENGTH)
    if not show_x:
        ax.tick_params(labelbottom=False)


# --- Statistics & violin (identical logic to stimulation_analysis.py) ---

def _perform_stat(data1, data2, paired=False):
    """Returns (p_value, test_name, stat, n1, n2)."""
    d1 = np.array([v for v in data1 if v is not None and np.isfinite(v)])
    d2 = np.array([v for v in data2 if v is not None and np.isfinite(v)])
    if len(d1) < 3 or len(d2) < 3:
        return 1.0, "insufficient_data", float("nan"), len(d1), len(d2)
    try:
        _, p1 = stats.shapiro(d1) if len(d1) <= 5000 else (None, 0.05)
        _, p2 = stats.shapiro(d2) if len(d2) <= 5000 else (None, 0.05)
        normal = (p1 > 0.05) and (p2 > 0.05)
    except Exception:
        normal = False
    if paired:
        if normal and len(d1) == len(d2):
            stat, pv = stats.ttest_rel(d1, d2)
            return pv, "paired_ttest", stat, len(d1), len(d2)
        if len(d1) == len(d2):
            try:
                stat, pv = stats.wilcoxon(d1, d2, alternative="two-sided")
                return pv, "wilcoxon", stat, len(d1), len(d2)
            except Exception:
                pass
        if normal:
            stat, pv = stats.ttest_ind(d1, d2)
            return pv, "ttest_ind", stat, len(d1), len(d2)
        stat, pv = stats.mannwhitneyu(d1, d2, alternative="two-sided")
        return pv, "mannwhitney", stat, len(d1), len(d2)
    else:
        if normal:
            stat, pv = stats.ttest_ind(d1, d2)
            return pv, "ttest_ind", stat, len(d1), len(d2)
        stat, pv = stats.mannwhitneyu(d1, d2, alternative="two-sided")
        return pv, "mannwhitney", stat, len(d1), len(d2)


def _holm_bonferroni(p_values):
    """Holm-Bonferroni step-down correction."""
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    corrected = [None] * n
    cum_max = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        adj = min(p * (n - rank), 1.0)
        cum_max = max(cum_max, adj)
        corrected[orig_idx] = cum_max
    return corrected


def _perform_omnibus(data_dict):
    """Friedman test for repeated-measures (non-parametric). Returns (chi2, p, k, n) or None."""
    arrays = []
    for vals in data_dict.values():
        clean = [v for v in vals if v is not None and np.isfinite(v)]
        arrays.append(clean)
    lengths = [len(a) for a in arrays]
    if len(arrays) < 3 or min(lengths) < 3:
        return None
    min_len = min(lengths)
    trimmed = [a[:min_len] for a in arrays]
    try:
        chi2, pv = stats.friedmanchisquare(*trimmed)
        return chi2, pv, len(trimmed), min_len
    except Exception:
        return None


def _add_bracket(ax, x1, x2, y, pv, lh, to):
    txt = "***" if pv < 0.001 else ("**" if pv < 0.01 else ("*" if pv < 0.05 else "ns"))
    ax.plot([x1, x2], [y, y], "k-", lw=1.5, clip_on=False)
    ax.plot([x1, x1], [y - lh, y], "k-", lw=1.5, clip_on=False)
    ax.plot([x2, x2], [y - lh, y], "k-", lw=1.5, clip_on=False)
    ax.text((x1 + x2) / 2, y + to, txt, ha="center", va="bottom",
            fontsize=FONT_SIZE_TICK - 1, fontweight="bold")


def plot_violin_box(ax, data_dict, ylabel, colors_dict, title=None, comparisons=None,
                    correction="holm", omnibus=False, context_str=""):
    """Publication violin-box plot with verbose statistical output."""
    positions, data_list, labels, colors_list = [], [], [], []
    label_to_pos = {}
    for i, (label, values) in enumerate(data_dict.items()):
        valid = [v for v in values if v is not None and np.isfinite(v)]
        if len(valid) > 0:
            positions.append(i); data_list.append(valid); labels.append(label)
            colors_list.append(colors_dict.get(label, [0.5, 0.5, 0.5]))
            label_to_pos[label] = i
    if not data_list:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=FONT_SIZE_LABEL)
        if title:
            ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold")
        return
    parts = ax.violinplot(data_list, positions=positions, widths=0.75,
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, c in zip(parts["bodies"], colors_list):
        pc.set_facecolor(c); pc.set_alpha(0.35)
        pc.set_edgecolor(np.clip(np.array(c) * 0.3, 0, 1)); pc.set_linewidth(3.5)
    bp = ax.boxplot(data_list, positions=positions, widths=0.22,
                    patch_artist=True, showfliers=False)
    for patch, c in zip(bp["boxes"], colors_list):
        patch.set_facecolor(np.clip(np.array(c) * 0.6 + 0.4, 0, 1))
        patch.set_edgecolor(np.clip(np.array(c) * 0.4, 0, 1))
        patch.set_linewidth(2.5); patch.set_alpha(0.9)
    for w in bp["whiskers"]:
        w.set_color("black"); w.set_linewidth(2.0)
    for m in bp["medians"]:
        m.set_color("black"); m.set_linewidth(3.0)
    for c in bp["caps"]:
        c.set_color("black"); c.set_linewidth(2.0)
    for i, (pos, vals) in enumerate(zip(positions, data_list)):
        jitter = np.random.uniform(-0.10, 0.10, len(vals))
        ax.scatter(pos + jitter, vals, color=colors_list[i], s=70, alpha=0.9,
                   zorder=10, edgecolors="black", linewidths=2.0)
    panel_label = title if title else ylabel
    if omnibus and len(data_list) >= 3:
        omni = _perform_omnibus(data_dict)
        if omni is not None:
            chi2, p_omni, k, n = omni
            print(f"    OMNIBUS [{panel_label}]: Friedman chi2={chi2:.4f}, p={p_omni:.4e}, k={k} groups, n={n} per group")
        else:
            print(f"    OMNIBUS [{panel_label}]: Friedman test skipped (insufficient data)")

    if comparisons:
        y_max = max(max(v) for v in data_list)
        y_min = min(min(v) for v in data_list)
        y_range = max(y_max - y_min, 1e-6)
        by = y_max + 0.06 * y_range
        bs = 0.07 * y_range
        lh = 0.015 * y_range
        to = 0.005 * y_range

        raw_results = []
        for idx, (l1, l2, paired) in enumerate(comparisons):
            if l1 in label_to_pos and l2 in label_to_pos:
                pv, tname, stat, n1, n2 = _perform_stat(
                    data_dict.get(l1, []), data_dict.get(l2, []), paired)
                raw_results.append((idx, l1, l2, paired, pv, tname, stat, n1, n2))
            else:
                raw_results.append((idx, l1, l2, paired, None, "skipped", float("nan"), 0, 0))

        raw_pvals = [r[4] for r in raw_results if r[4] is not None]
        if correction == "holm" and len(raw_pvals) > 1:
            corrected_pvals = _holm_bonferroni(raw_pvals)
        else:
            corrected_pvals = list(raw_pvals)

        corr_idx = 0
        ctx = f" {context_str}" if context_str else ""
        print(f"    STATS [{panel_label}]{ctx} — {len(raw_results)} comparisons, correction={correction}:")
        for r in raw_results:
            idx_c, l1, l2, paired, p_raw, tname, stat, n1, n2 = r
            if p_raw is not None:
                p_corr = corrected_pvals[corr_idx]
                corr_idx += 1
                paired_str = "paired" if paired else "unpaired"
                print(f"      {l1} vs {l2}: {tname} ({paired_str}), stat={stat:.4f}, "
                      f"n1={n1}, n2={n2}, p_raw={p_raw:.4e}, p_corrected={p_corr:.4e}")
                p_use = p_corr
            else:
                print(f"      {l1} vs {l2}: skipped (label not in data)")
                p_use = 1.0

            if l1 in label_to_pos and l2 in label_to_pos:
                _add_bracket(ax, label_to_pos[l1], label_to_pos[l2],
                             by + idx_c * bs, p_use, lh, to)

        ax.set_ylim(ax.get_ylim()[0], by + len(raw_results) * bs + 0.05 * y_range)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK, rotation=45, ha="right")
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis="both", which="major", labelsize=FONT_SIZE_TICK,
                   width=TICK_WIDTH, length=TICK_LENGTH)
    if title:
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold")


# =============================================================================
# FIGURE 1 & 2: TRIAL-AVERAGED TRACES (per condition)
# =============================================================================

def fig_trial_avg_traces(data, condition_label):
    """LFP (trial-avg ± SEM) and neuron-avg voltage (trial-avg ± SEM), 2 panels."""
    t, lfp_all, volt_all, fs = _collect_trial_traces(data)
    n_trials = volt_all.shape[0]
    lfp_m, lfp_s = _mean_sem(lfp_all)
    volt_m, volt_s = _mean_sem(volt_all)
    lfp_m, lfp_s = _smooth(lfp_m), _smooth(lfp_s)
    volt_m, volt_s = _smooth(volt_m), _smooth(volt_s)

    fig, (ax_lfp, ax_v) = plt.subplots(2, 1, figsize=(16, 10), sharex=True,
                                        gridspec_kw={"hspace": 0.12})
    ax_lfp.fill_between(t, lfp_m - lfp_s, lfp_m + lfp_s, color=COLOR_LFP, alpha=0.3, lw=0)
    ax_lfp.plot(t, lfp_m, color=COLOR_LFP, lw=LINE_WIDTH_THICK)
    _style_avg_axis(ax_lfp, show_x=False)
    ax_lfp.set_ylabel("LFP (µV)", fontsize=FONT_SIZE_LABEL)
    ax_lfp.set_title(f"Trial-Averaged Traces (n={n_trials} trials) ± SEM",
                     fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)

    ax_v.fill_between(t, volt_m - volt_s, volt_m + volt_s, color=COLOR_NEURON, alpha=0.3, lw=0)
    ax_v.plot(t, volt_m, color=COLOR_NEURON, lw=LINE_WIDTH_THICK)
    _style_avg_axis(ax_v)
    ax_v.set_ylabel("Population Vm (ΔF/F)", fontsize=FONT_SIZE_LABEL)
    ax_v.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)

    fig.suptitle(f"{MOUSE_NAME} - {condition_label}",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.07)
    return fig


# =============================================================================
# FIGURE 3: OVERLAID TRIAL-AVERAGED TRACES (40 vs 135)
# =============================================================================

def fig_overlaid_avg(data_r1, data_r2, label_r1, label_r2):
    t1, lfp1, v1, _ = _collect_trial_traces(data_r1)
    t2, lfp2, v2, _ = _collect_trial_traces(data_r2)
    n1, n2 = v1.shape[0], v2.shape[0]

    # Use longest common time grid (both should be very similar)
    n_pts = min(len(t1), len(t2))
    t = t1[:n_pts]

    def _trim(a):
        return a[:, :n_pts]

    lfp1_m, lfp1_s = _mean_sem(_trim(lfp1)); lfp2_m, lfp2_s = _mean_sem(_trim(lfp2))
    v1_m, v1_s = _mean_sem(_trim(v1));       v2_m, v2_s = _mean_sem(_trim(v2))
    for arr in (lfp1_m, lfp1_s, lfp2_m, lfp2_s, v1_m, v1_s, v2_m, v2_s):
        arr[:] = _smooth(arr)

    fig, (ax_lfp, ax_v) = plt.subplots(2, 1, figsize=(16, 10), sharex=True,
                                        gridspec_kw={"hspace": 0.12})

    def _plot_pair(ax, m1, s1, m2, s2, ylabel, title):
        ax.fill_between(t, m1 - s1, m1 + s1, color=COLOR_40HZ, alpha=0.3, lw=0)
        ax.plot(t, m1, color=COLOR_40HZ, lw=LINE_WIDTH_THICK, label=f"{label_r1} (n={n1})")
        ax.fill_between(t, m2 - s2, m2 + s2, color=COLOR_135HZ, alpha=0.3, lw=0)
        ax.plot(t, m2, color=COLOR_135HZ, lw=LINE_WIDTH_THICK, label=f"{label_r2} (n={n2})")
        _style_avg_axis(ax, show_x=(ax is ax_v))
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)
        ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, framealpha=0.9)
        # Shared y from both envelopes
        lo = min(np.nanmin(m1 - s1), np.nanmin(m2 - s2))
        hi = max(np.nanmax(m1 + s1), np.nanmax(m2 + s2))
        span = max(hi - lo, 1e-6)
        ax.set_ylim(lo - 0.05 * span, hi + 0.05 * span)

    _plot_pair(ax_lfp, lfp1_m, lfp1_s, lfp2_m, lfp2_s,
               "LFP (µV)", f"Trial-Averaged LFP ± SEM (n={n1} vs {n2})")
    _plot_pair(ax_v, v1_m, v1_s, v2_m, v2_s,
               "Population Vm (ΔF/F)", f"Trial-Averaged Voltage ± SEM (n={n1} vs {n2})")
    ax_v.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)

    fig.suptitle(f"{MOUSE_NAME} - Trial-averaged overlaid ({label_r1} vs {label_r2})",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.07)
    return fig


# =============================================================================
# FIGURE 4 & 5: PERIOD VIOLIN (per condition)
# =============================================================================

def fig_period_violin(data, condition_label):
    t, _, volt_all, _ = _collect_trial_traces(data)
    n_trials = volt_all.shape[0]
    n_neurons = data["num_neurons"]
    pv = _period_values(t, volt_all)
    colors = {"Pre-stim": COLOR_PRE, "Transient": COLOR_TRANSIENT,
              "Sustained": COLOR_SUSTAINED, "Post-stim": COLOR_POST}
    comps = [
        ("Pre-stim", "Transient", True), ("Pre-stim", "Sustained", True),
        ("Pre-stim", "Post-stim", True), ("Transient", "Sustained", True),
        ("Transient", "Post-stim", True), ("Sustained", "Post-stim", True),
    ]
    print(f"    Population Vm: {n_neurons} neurons averaged per trial, {n_trials} trials")
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_violin_box(ax, pv, "Voltage Signal (ΔF/F)", colors,
                    "Voltage Signal by Period", comparisons=comps,
                    correction="holm", omnibus=True,
                    context_str=f"(n={n_trials} trials, {n_neurons} neurons)")
    fig.suptitle(f"{MOUSE_NAME} - {condition_label} Period Comparison",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
    return fig


# =============================================================================
# FIGURE 6: 40 vs 135 Hz TRANS / SUST VIOLIN
# =============================================================================

def fig_freq_comparison_violin(data_r1, data_r2, label_r1, label_r2):
    t1, _, v1, _ = _collect_trial_traces(data_r1)
    t2, _, v2, _ = _collect_trial_traces(data_r2)
    n1, n2 = v1.shape[0], v2.shape[0]
    nn1, nn2 = data_r1["num_neurons"], data_r2["num_neurons"]
    pv1 = _period_values(t1, v1)
    pv2 = _period_values(t2, v2)

    data_dict = {
        "40Hz Trans": pv1["Transient"], "40Hz Sust": pv1["Sustained"],
        "130Hz Trans": pv2["Transient"], "130Hz Sust": pv2["Sustained"],
    }
    colors = {"40Hz Trans": COLOR_40HZ, "40Hz Sust": COLOR_40HZ_LIGHT,
              "130Hz Trans": COLOR_135HZ, "130Hz Sust": COLOR_135HZ_LIGHT}
    comps = [
        ("40Hz Trans", "40Hz Sust", True),
        ("130Hz Trans", "130Hz Sust", True),
        ("40Hz Trans", "130Hz Trans", False),
        ("40Hz Sust", "130Hz Sust", False),
    ]
    print(f"    Comparison: 40Hz n={n1} trials/{nn1} neurons, 130Hz n={n2} trials/{nn2} neurons")
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_violin_box(ax, data_dict, "Voltage Signal (ΔF/F)", colors,
                    "Voltage Signal: 40Hz vs 130Hz", comparisons=comps,
                    correction="holm",
                    context_str=f"(40Hz: {n1} trials/{nn1} neur, 130Hz: {n2} trials/{nn2} neur)")
    fig.suptitle(f"{MOUSE_NAME} - 40Hz vs 130Hz Stimulation Comparison",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.15, right=0.95, top=0.90, bottom=0.15)
    return fig


# =============================================================================
# FIGURE 7a & 7b: TRIAL-AVERAGED SPECTROGRAMS (LFP + neuron 1)
# =============================================================================

from common import _next_pow2  # shared helpers (were local copies)


def _smooth2d(arr, n_freq=1, n_time=1):
    """Replicate MATLAB smooth2a: uniform-filter along each axis independently."""
    out = arr.copy()
    if n_freq > 1:
        out = uniform_filter1d(out, size=n_freq, axis=0, mode="nearest")
    if n_time > 1:
        out = uniform_filter1d(out, size=n_time, axis=1, mode="nearest")
    return out


def _compute_spectrogram(trace, fs):
    """
    Spectrogram matching MATLAB pipeline parameters:
    window = SPEC_WINDOW_SEC × fs, overlap = 88%, NFFT = nextpow2(3 × window),
    followed by 2-D smoothing.
    Returns (freq, time_bins, power) where power is smoothed |S|².
    """
    nperseg = round(SPEC_WINDOW_SEC * fs)
    nperseg = min(nperseg, len(trace))
    noverlap = round(SPEC_OVERLAP_FRAC * nperseg)
    nfft = _next_pow2(nperseg * SPEC_NFFT_MULT)

    freq, time_bins, Sxx = signal.spectrogram(
        trace, fs=fs, nperseg=nperseg, noverlap=noverlap, nfft=nfft,
        window="hann", scaling="density", mode="psd",
    )
    Sxx = _smooth2d(Sxx, SPEC_SMOOTH_FREQ, SPEC_SMOOTH_TIME)
    return freq, time_bins, Sxx


def _baseline_mask(t):
    m = (t >= SPEC_BASELINE_T_LO) & (t < SPEC_BASELINE_T_HI)
    if not np.any(m):
        m = t < 0.0
    return m


def _prepare_trialavg_spectrogram(data):
    """
    Prepare trial-averaged baseline-normalized spectrograms for LFP and neuron.
    Returns dict with freq/time axes and normalized matrices.
    """
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    ref = data["trials"][0]
    fs = ref["fs"]

    min_samples = round(SPEC_WINDOW_SEC * fs) + 1
    lfp_specs, neur_specs = [], []
    freq_out, time_out = None, None

    for trial in data["trials"]:
        t_tr = trial["time_vector"]
        m = _window_mask(t_tr, t_min, t_max)
        if np.sum(m) < min_samples:
            continue

        if trial["lfp"] is not None:
            lfp_seg = trial["lfp"][m]
            f, tb, Sxx = _compute_spectrogram(lfp_seg, fs)
            tb = tb + t_min
            lfp_specs.append(Sxx)
            if freq_out is None:
                freq_out, time_out = f, tb

        fl = trial["fluorescence"]
        if fl.ndim == 2 and fl.shape[0] >= np.sum(m) and NEURON_IDX_SPEC < fl.shape[1]:
            neur_seg = fl[m, NEURON_IDX_SPEC]
            f, tb, Sxx = _compute_spectrogram(neur_seg, fs)
            neur_specs.append(Sxx)

    if freq_out is None or len(lfp_specs) == 0:
        return None

    lfp_avg = np.nanmean(np.array(lfp_specs), axis=0)
    neur_avg = np.nanmean(np.array(neur_specs), axis=0) if neur_specs else None

    def _frac_change(power, t_axis):
        bm = _baseline_mask(t_axis)
        if not np.any(bm):
            return power - np.nanmean(power)
        bl = np.nanmean(power[:, bm], axis=1, keepdims=True)
        num = power - bl
        den = power + bl
        den[den == 0] = 1e-10
        return num / den

    lfp_norm = _frac_change(lfp_avg, time_out)
    has_neuron = neur_avg is not None
    neur_norm = _frac_change(neur_avg, time_out) if has_neuron else None
    fmask = (freq_out >= FREQ_RANGE[0]) & (freq_out <= FREQ_RANGE[1])

    return {
        "time_out": time_out,
        "freq_out": freq_out,
        "fmask": fmask,
        "lfp_norm": lfp_norm,
        "neur_norm": neur_norm,
        "has_neuron": has_neuron,
    }


def _compute_global_spec_limits(spec_prepared_list):
    """
    Compute global color limits across conditions, separately for LFP and neuron.
    """
    lfp_vals = []
    neur_vals = []
    for s in spec_prepared_list:
        if s is None:
            continue
        fmask = s["fmask"]
        lfp_vals.append(s["lfp_norm"][fmask, :].ravel())
        if s["neur_norm"] is not None:
            neur_vals.append(s["neur_norm"][fmask, :].ravel())

    if not lfp_vals:
        return None

    lfp_all = np.concatenate(lfp_vals)
    lfp_lo = float(np.nanpercentile(lfp_all, SPEC_LFP_PCTL[0]))
    lfp_hi = float(np.nanpercentile(lfp_all, SPEC_LFP_PCTL[1]))
    # Small expansion margin + enforce symmetric limits around zero.
    lfp_span = lfp_hi - lfp_lo
    if np.isfinite(lfp_span) and lfp_span > 1e-12:
        lfp_lo -= 0.05 * lfp_span
        lfp_hi += 0.08 * lfp_span
    lfp_abs = max(abs(lfp_lo), abs(lfp_hi))
    lfp_vmin, lfp_vmax = -lfp_abs, lfp_abs

    if neur_vals:
        neur_all = np.concatenate(neur_vals)
        neur_lo = float(np.nanpercentile(neur_all, SPEC_NEUR_PCTL[0]))
        neur_hi = float(np.nanpercentile(neur_all, SPEC_NEUR_PCTL[1]))
        neur_abs = max(abs(neur_lo), abs(neur_hi))
        neur_vmin, neur_vmax = -neur_abs, neur_abs
    else:
        neur_vmin, neur_vmax = lfp_vmin, lfp_vmax

    return {
        "lfp": (lfp_vmin, lfp_vmax),
        "neur": (neur_vmin, neur_vmax),
    }


def fig_spectrogram(data, condition_label, shared_limits=None):
    """
    Trial-averaged spectrogram for LFP and neuron NEURON_IDX_SPEC.
    Fractional change (P-B)/(P+B) in linear power.
    """
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    prepared = _prepare_trialavg_spectrogram(data)
    if prepared is None:
        fig, _ = plt.subplots(); return fig
    time_out = prepared["time_out"]
    freq_out = prepared["freq_out"]
    fmask = prepared["fmask"]
    lfp_norm = prepared["lfp_norm"]
    neur_norm = prepared["neur_norm"]
    has_neuron = prepared["has_neuron"]

    if shared_limits is None:
        lfp_lo = float(np.nanpercentile(lfp_norm[fmask, :], SPEC_LFP_PCTL[0]))
        lfp_hi = float(np.nanpercentile(lfp_norm[fmask, :], SPEC_LFP_PCTL[1]))
        lfp_span = lfp_hi - lfp_lo
        if np.isfinite(lfp_span) and lfp_span > 1e-12:
            lfp_lo -= 0.05 * lfp_span
            lfp_hi += 0.08 * lfp_span
        lfp_abs = max(abs(lfp_lo), abs(lfp_hi))
        lfp_vmin, lfp_vmax = -lfp_abs, lfp_abs
        if neur_norm is not None:
            neur_lo = float(np.nanpercentile(neur_norm[fmask, :], SPEC_NEUR_PCTL[0]))
            neur_hi = float(np.nanpercentile(neur_norm[fmask, :], SPEC_NEUR_PCTL[1]))
            neur_abs = max(abs(neur_lo), abs(neur_hi))
            neur_vmin, neur_vmax = -neur_abs, neur_abs
        else:
            neur_vmin, neur_vmax = lfp_vmin, lfp_vmax
    else:
        lfp_vmin, lfp_vmax = shared_limits["lfp"]
        neur_vmin, neur_vmax = shared_limits["neur"]

    n_rows = 2 if has_neuron else 1
    fig = plt.figure(figsize=(14, 5 * n_rows))
    gs = GridSpec(n_rows, 2, figure=fig,
                  height_ratios=([1] * n_rows),
                  width_ratios=[1, 0.03], hspace=0.25, wspace=0.05)

    # LFP spectrogram
    ax1 = fig.add_subplot(gs[0, 0])
    cax1 = fig.add_subplot(gs[0, 1])
    im1 = ax1.pcolormesh(time_out, freq_out[fmask], lfp_norm[fmask, :],
                         shading="gouraud", cmap="viridis",
                         vmin=lfp_vmin, vmax=lfp_vmax)
    ax1.axvline(0, color="white", ls="--", lw=1.5)
    ax1.axvline(STIM_DURATION_SEC, color="white", ls=":", lw=1.5)
    ax1.set_ylabel("LFP\nFrequency (Hz)", fontsize=FONT_SIZE_LABEL,
                   rotation=0, ha="right", va="center")
    if has_neuron:
        ax1.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK)
    else:
        ax1.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)
        ax1.tick_params(labelsize=FONT_SIZE_TICK)
    ax1.set_xlim(t_min, t_max)
    ax1.set_ylim(FREQ_RANGE)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax1.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    cb1 = fig.colorbar(im1, cax=cax1, orientation="vertical")
    cb1.set_label("Fractional Change\n(rel. baseline)", fontsize=FONT_SIZE_TICK)
    cb1.ax.tick_params(labelsize=FONT_SIZE_TICK - 2)

    # Neuron spectrogram
    if has_neuron:
        ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
        cax2 = fig.add_subplot(gs[1, 1])
        im2 = ax2.pcolormesh(time_out, freq_out[fmask], neur_norm[fmask, :],
                             shading="gouraud", cmap="viridis",
                             vmin=neur_vmin, vmax=neur_vmax)
        ax2.axvline(0, color="white", ls="--", lw=1.5)
        ax2.axvline(STIM_DURATION_SEC, color="white", ls=":", lw=1.5)
        ax2.set_ylabel(f"Neuron {NEURON_IDX_SPEC+1}\nFrequency (Hz)",
                       fontsize=FONT_SIZE_LABEL, rotation=0, ha="right", va="center")
        ax2.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)
        ax2.tick_params(labelsize=FONT_SIZE_TICK)
        ax2.set_ylim(FREQ_RANGE)
        for sp in ("top", "right"):
            ax2.spines[sp].set_visible(False)
        ax2.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax2.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
        cb2 = fig.colorbar(im2, cax=cax2, orientation="vertical")
        cb2.set_label("Fractional Change\n(rel. baseline)", fontsize=FONT_SIZE_TICK)
        cb2.ax.tick_params(labelsize=FONT_SIZE_TICK - 2)

    ax1.set_title(f"{MOUSE_NAME} - {condition_label} Spectral Analysis (trial-averaged)",
                  fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.14, right=0.88, top=0.92, bottom=0.08)
    return fig


def fig_spectrogram_fiber_only(data, condition_label, shared_limits=None):
    """
    Trial-averaged fiber-only spectrogram (same fiber panel used in fig_spectrogram).
    """
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    prepared = _prepare_trialavg_spectrogram(data)
    if prepared is None or prepared["neur_norm"] is None:
        fig, _ = plt.subplots()
        return fig

    time_out = prepared["time_out"]
    freq_out = prepared["freq_out"]
    fmask = prepared["fmask"]
    neur_norm = prepared["neur_norm"]

    if shared_limits is None:
        neur_lo = float(np.nanpercentile(neur_norm[fmask, :], SPEC_NEUR_PCTL[0]))
        neur_hi = float(np.nanpercentile(neur_norm[fmask, :], SPEC_NEUR_PCTL[1]))
        neur_abs = max(abs(neur_lo), abs(neur_hi))
        neur_vmin, neur_vmax = -neur_abs, neur_abs
    else:
        neur_vmin, neur_vmax = shared_limits["neur"]

    # Standalone fiber spectrogram: narrower and taller for publication layout.
    fig = plt.figure(figsize=(8.4, 7.2))
    gs_st = GridSpec(1, 2, figure=fig, width_ratios=[1, 0.03], wspace=0.05)
    ax = fig.add_subplot(gs_st[0, 0])
    cax = fig.add_subplot(gs_st[0, 1])

    im = ax.pcolormesh(
        time_out, freq_out[fmask], neur_norm[fmask, :],
        shading="gouraud", cmap="viridis", vmin=neur_vmin, vmax=neur_vmax
    )
    ax.axvline(0, color="white", ls="--", lw=1.5)
    ax.axvline(STIM_DURATION_SEC, color="white", ls=":", lw=1.5)
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(FREQ_RANGE)
    ax.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel(f"Neuron {NEURON_IDX_SPEC+1}\nFrequency (Hz)",
                  fontsize=FONT_SIZE_LABEL, rotation=0, ha="right", va="center")
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)

    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.set_label("Fractional Change\n(rel. baseline)", fontsize=FONT_SIZE_TICK)
    cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 2)

    ax.set_title(f"{MOUSE_NAME} - {condition_label} Fiber Spectrogram (trial-averaged)",
                 fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.14, right=0.90, top=0.90, bottom=0.13)
    return fig


def fig_spectrogram_single_trial(data, condition_label, trial_idx=0):
    """LFP spectrogram for a single trial (no averaging) — diagnostic figure."""
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    trial = data["trials"][trial_idx]
    fs = trial["fs"]
    t_tr = trial["time_vector"]
    m = _window_mask(t_tr, t_min, t_max)

    lfp_seg = trial["lfp"][m]
    freq, time_bins, Sxx = _compute_spectrogram(lfp_seg, fs)
    time_bins = time_bins + t_min

    bm = _baseline_mask(time_bins)
    if np.any(bm):
        bl = np.nanmean(Sxx[:, bm], axis=1, keepdims=True)
        den = Sxx + bl; den[den == 0] = 1e-10
        Sxx_norm = (Sxx - bl) / den
    else:
        Sxx_norm = Sxx - np.nanmean(Sxx)

    fmask = (freq >= FREQ_RANGE[0]) & (freq <= FREQ_RANGE[1])

    fig = plt.figure(figsize=(14, 5))
    gs_st = GridSpec(1, 2, figure=fig, width_ratios=[1, 0.03], wspace=0.05)
    ax = fig.add_subplot(gs_st[0, 0])
    cax = fig.add_subplot(gs_st[0, 1])

    im = ax.pcolormesh(time_bins, freq[fmask], Sxx_norm[fmask, :],
                       shading="gouraud", cmap="viridis")
    ax.axvline(0, color="white", ls="--", lw=1.5)
    ax.axvline(STIM_DURATION_SEC, color="white", ls=":", lw=1.5)
    ax.set_ylabel("LFP\nFrequency (Hz)", fontsize=FONT_SIZE_LABEL,
                  rotation=0, ha="right", va="center")
    ax.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(FREQ_RANGE)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
    ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    cb = fig.colorbar(im, cax=cax, orientation="vertical")
    cb.set_label("Fractional Change\n(rel. baseline)", fontsize=FONT_SIZE_TICK)
    cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 2)
    ax.set_title(f"{MOUSE_NAME} - {condition_label} LFP Spectrogram "
                 f"(single trial {trial_idx + 1})",
                 fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=10)
    fig.subplots_adjust(left=0.14, right=0.88, top=0.90, bottom=0.12)
    return fig


# =============================================================================
# FIGURE 8: STIM-BAND POWER VIOLIN (trans vs sust, 40 vs 135)
# =============================================================================

def _band_power_per_trial(data, freq_band):
    """
    Per-trial mean band power (dB) in four periods, computed from Welch spectrogram.
    Returns dict {period_label: [dB values per trial]}.
    """
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    ref = data["trials"][0]
    fs = ref["fs"]

    trial_powers = {"Pre-stim": [], "Transient": [], "Sustained": [], "Post-stim": []}

    min_samples = round(SPEC_WINDOW_SEC * fs) + 1

    for trial in data["trials"]:
        t_tr = trial["time_vector"]
        m = _window_mask(t_tr, t_min, t_max)
        if np.sum(m) < min_samples:
            for k in trial_powers:
                trial_powers[k].append(np.nan)
            continue

        fl = trial["fluorescence"]
        if fl.ndim == 2 and fl.shape[0] >= np.sum(m):
            seg = np.mean(fl[m, :], axis=1)
        else:
            for k in trial_powers:
                trial_powers[k].append(np.nan)
            continue

        freq, tb, Sxx = _compute_spectrogram(seg, fs)
        tb = tb + t_min
        fmask = (freq >= freq_band[0]) & (freq <= freq_band[1])
        band_power = np.mean(Sxx[fmask, :], axis=0)  # time series of band power
        band_db = 10 * np.log10(band_power + 1e-20)

        periods = {
            "Pre-stim": tb < 0,
            "Transient": (tb >= 0) & (tb < TRANSIENT_END),
            "Sustained": (tb >= SUSTAINED_START) & (tb < STIM_DURATION_SEC),
            "Post-stim": tb >= STIM_DURATION_SEC,
        }
        for label, pm in periods.items():
            if np.any(pm):
                trial_powers[label].append(float(np.nanmean(band_db[pm])))
            else:
                trial_powers[label].append(np.nan)

    return trial_powers


def fig_stim_band_violin(data_r1, data_r2, label_r1, label_r2):
    bp1 = _band_power_per_trial(data_r1, BAND_40HZ)
    bp2 = _band_power_per_trial(data_r2, BAND_130HZ)

    pre_mean_40 = np.nanmean(bp1["Pre-stim"])
    pre_mean_130 = np.nanmean(bp2["Pre-stim"])

    data_dict = {
        "40Hz Trans": [v - pre_mean_40 for v in bp1["Transient"]],
        "40Hz Sust": [v - pre_mean_40 for v in bp1["Sustained"]],
        "130Hz Trans": [v - pre_mean_130 for v in bp2["Transient"]],
        "130Hz Sust": [v - pre_mean_130 for v in bp2["Sustained"]],
    }
    colors = {"40Hz Trans": COLOR_40HZ, "40Hz Sust": COLOR_40HZ_LIGHT,
              "130Hz Trans": COLOR_135HZ, "130Hz Sust": COLOR_135HZ_LIGHT}
    comps = [
        ("40Hz Trans", "40Hz Sust", True),
        ("130Hz Trans", "130Hz Sust", True),
        ("40Hz Trans", "130Hz Trans", False),
        ("40Hz Sust", "130Hz Sust", False),
    ]
    n1, n2 = data_r1["num_trials"], data_r2["num_trials"]
    nn1, nn2 = data_r1["num_neurons"], data_r2["num_neurons"]
    print(f"    Stim-band power: 40Hz n={n1} trials/{nn1} neurons, 130Hz n={n2} trials/{nn2} neurons")
    fig, ax = plt.subplots(figsize=(10, 8))
    plot_violin_box(ax, data_dict, "Band Power (dB re: baseline)", colors,
                    f"Stim-Band Power\n(40Hz:{BAND_40HZ[0]}-{BAND_40HZ[1]}Hz, "
                    f"130Hz:{BAND_130HZ[0]}-{BAND_130HZ[1]}Hz)",
                    comparisons=comps, correction="holm",
                    context_str=f"(40Hz: {n1} trials/{nn1} neur, 130Hz: {n2} trials/{nn2} neur)")
    fig.suptitle(f"{MOUSE_NAME} - Stim-Band Relative Power Comparison",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.15, right=0.95, top=0.88, bottom=0.15)
    return fig


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("  CELLULAR DBS COMPARISON FIGURES")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"\n[1/2] Loading R1 ({SESSION_R1['label']})...")
    data_r1 = load_cellular_data(SESSION_R1)
    print(f"\n[2/2] Loading R2 ({SESSION_R2['label']})...")
    data_r2 = load_cellular_data(SESSION_R2)

    if data_r1 is None or data_r2 is None:
        print("ERROR: Could not load one or both sessions. Aborting.")
        return

    label_r1 = SESSION_R1["label"]
    label_r2 = SESSION_R2["label"]
    safe_r1 = label_r1.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")
    safe_r2 = label_r2.replace(" ", "_").replace("(", "").replace(")", "").replace("-", "")

    figures = {}

    print("\nCreating figures...")

    print("  Fig 1: Trial-avg traces R1")
    figures["avg_r1"] = fig_trial_avg_traces(data_r1, label_r1)

    print("  Fig 2: Trial-avg traces R2")
    figures["avg_r2"] = fig_trial_avg_traces(data_r2, label_r2)

    print("  Fig 3: Overlaid trial-avg traces")
    figures["overlay"] = fig_overlaid_avg(data_r1, data_r2, label_r1, label_r2)

    print("  Fig 4: Period violin R1")
    figures["period_r1"] = fig_period_violin(data_r1, label_r1)

    print("  Fig 5: Period violin R2")
    figures["period_r2"] = fig_period_violin(data_r2, label_r2)

    print("  Fig 6: 40 vs 135 Hz comparison violin")
    figures["freq_comp"] = fig_freq_comparison_violin(data_r1, data_r2, label_r1, label_r2)

    print("  Computing shared spectrogram limits across conditions (separate for LFP and neuron)...")
    spec_limits = _compute_global_spec_limits([
        _prepare_trialavg_spectrogram(data_r1),
        _prepare_trialavg_spectrogram(data_r2),
    ])

    print("  Fig 7a: Spectrogram R1")
    figures["spec_r1"] = fig_spectrogram(data_r1, label_r1, shared_limits=spec_limits)

    print("  Fig 7b: Spectrogram R2")
    figures["spec_r2"] = fig_spectrogram(data_r2, label_r2, shared_limits=spec_limits)

    print("  Fig 7e: Fiber-only spectrogram R1")
    figures["spec_fiber_r1"] = fig_spectrogram_fiber_only(data_r1, label_r1, shared_limits=spec_limits)

    print("  Fig 7f: Fiber-only spectrogram R2")
    figures["spec_fiber_r2"] = fig_spectrogram_fiber_only(data_r2, label_r2, shared_limits=spec_limits)

    print("  Fig 7c: Single-trial LFP spectrogram R1 (diagnostic)")
    figures["spec_single_r1"] = fig_spectrogram_single_trial(data_r1, label_r1, trial_idx=0)

    print("  Fig 7d: Single-trial LFP spectrogram R2 (diagnostic)")
    figures["spec_single_r2"] = fig_spectrogram_single_trial(data_r2, label_r2, trial_idx=0)

    print("  Fig 8: Stim-band power violin")
    figures["band_violin"] = fig_stim_band_violin(data_r1, data_r2, label_r1, label_r2)

    output_names = {
        "avg_r1":      f"{MOUSE_NAME}_01_trial_avg_{safe_r1}",
        "avg_r2":      f"{MOUSE_NAME}_02_trial_avg_{safe_r2}",
        "overlay":     f"{MOUSE_NAME}_03_trial_avg_overlay_{safe_r1}_vs_{safe_r2}",
        "period_r1":   f"{MOUSE_NAME}_04_period_violin_{safe_r1}",
        "period_r2":   f"{MOUSE_NAME}_05_period_violin_{safe_r2}",
        "freq_comp":   f"{MOUSE_NAME}_06_freq_comparison_violin",
        "spec_r1":        f"{MOUSE_NAME}_07a_spectrogram_{safe_r1}",
        "spec_r2":        f"{MOUSE_NAME}_07b_spectrogram_{safe_r2}",
        "spec_fiber_r1":  f"{MOUSE_NAME}_07e_spectrogram_fiberonly_{safe_r1}",
        "spec_fiber_r2":  f"{MOUSE_NAME}_07f_spectrogram_fiberonly_{safe_r2}",
        "spec_single_r1": f"{MOUSE_NAME}_07c_spectrogram_single_trial_{safe_r1}",
        "spec_single_r2": f"{MOUSE_NAME}_07d_spectrogram_single_trial_{safe_r2}",
        "band_violin":    f"{MOUSE_NAME}_08_stim_band_violin",
    }

    print("\nSaving figures...")
    for key, fig in figures.items():
        name = output_names[key]
        for ext in ("pdf", "png"):
            fig.savefig(str(OUTPUT_DIR / f"{name}.{ext}"), dpi=DPI, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {name}")

    print(f"\nDONE! Created {len(figures)} figures in {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
