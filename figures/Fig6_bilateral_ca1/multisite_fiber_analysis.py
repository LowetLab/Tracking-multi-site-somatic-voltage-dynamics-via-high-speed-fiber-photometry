"""
Multi-site Fiber Voltage Imaging Analysis — Right HP vs Left HP.

Generates publication-quality figures from preprocessed .mat files:

Figure 1: Representative traces (2-column layout)
    Col 1 (Right HP): LFP(ipsi) → GEVI → theta-filtered overlay → motion
    Col 2 (Left HP):  LFP(Ch11/contra) → GEVI → theta-filtered overlay → motion

Figure 2: Time-frequency spectrograms (2-column layout)
    Col 1 (Right HP): speed heatmap → LFP TFR → fiber TFR
    Col 2 (Left HP):  speed heatmap → LFP TFR → fiber TFR

Figure 3: Fiber-fiber theta-band cross-correlation (REST vs RUN)
    Top: single-trial heatmaps (inferno) + all-trials per-trial correlograms (RdBu_r)
    Bottom: single-trial mean correlogram + all-trials grand mean +/- SEM + violin (peak |r|)

Figure 4: Hippocampal theta phase-locking (Left − Right HP)
    Trial-equalized roses (dusty coral); paired half-violin/box/scatter for per-trial R.
    Optional pool all sessions (PHASE_AGGREGATE_ALL_SESSIONS). Example: 1 s θ traces + phase.

Figure 5: Bilateral fiber spectral analysis
    Row 1: Time-resolved coherence (coherogram) per trial for example session + speed bar.
    Row 2: Right HP PSD (single-trial | all-trials mean±SEM | theta prominence violin).
    Row 3: Left HP PSD (same layout).
    Row 4: Fiber-fiber coherence (single-trial | all-trials mean±SEM | theta coherence violin).

Data mapping (from FiberPhotometryAnalysis struct):
    Fiber 1 (Right HP) = signals.final_processed_traces(:, 1) or deltaF_F_traces(:, 1)
    Fiber 2 (Left HP)  = signals.final_processed_traces(:, 2) or deltaF_F_traces(:, 2)
    Right HP LFP (ipsi) = ephys.lfp_raw_aligned_ipsiHP
    Left HP LFP (contra/Ch11) = ephys.lfp_raw_aligned_HP
    Motion = ephys.running_velocity_smooth
    Time = time.time_vector_seconds
    fs = time.sampling_rate
"""

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MultipleLocator
from scipy import signal
from scipy.stats import gaussian_kde, shapiro, ttest_rel, ttest_1samp, wilcoxon
from scipy.ndimage import uniform_filter1d, label as ndimage_label
from scipy.io import loadmat
from pathlib import Path
from functools import lru_cache
import warnings
import h5py

warnings.filterwarnings("ignore")

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# common.py lives in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# =============================================================================
# CONFIGURATION -- EDIT THIS to your own cohort
# =============================================================================

BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging"

OUTPUT_DIR = PROJECT_ROOT / "Figures" / "Multisite_Fiber_Analysis"

RECORDINGS = {
    "Animal01": {
        "01_01_26-R2": {"n_trials": 5, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
        "02_01_26-R2": {"n_trials": 3, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
    },
    # Animal02 bilateral HP: update n_trials if folders are added/removed on disk.
    "Animal02": {
        "01_01_26-R1": {"n_trials": 5, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
        "02_01_26-R2": {"n_trials": 4, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
        "05_02_26-R2": {"n_trials": 3, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
        "16_02_26-R2": {"n_trials": 6, "suffix": "fov1_HPipsicontradual_baselineRecording_60sec"},
    },
}

EXAMPLE_ANIMAL = "Animal01"
EXAMPLE_SESSION = "02_01_26-R2"
EXAMPLE_TRIAL = 1

# Signal processing
THETA_BAND = (5, 9)
FILTER_ORDER = 3

# Cross-correlation
XCORR_WINDOW_SEC = 2.0
XCORR_STEP_SEC = 0.25
XCORR_MAX_LAG_MS = 200
XCORR_EDGE_TRIM_SEC = 0.5

# Peak lag on *mean* correlogram: argmax |r| only within |lag| <= this window (physiological).
# Use float, e.g. 50.0 for ±50 ms. Use "half_theta" for ±½ period at band centre: 1000/(2*f_mid).
XCORR_PEAK_LAG_LIMIT_MS = "half_theta"

# Average across sliding windows: Fisher z (atanh) then mean then tanh, so high-|r| windows
# do not dominate the arithmetic mean of correlation coefficients.
XCORR_MEAN_USE_FISHER_Z = True

# All-trials cross-correlation (Figure 3 expanded layout)
XCORR_AGGREGATE_ALL_SESSIONS = True   # pool all sessions for an animal (like phase figure)
XCORR_ALL_TRIALS_CMAP = "RdBu_r"     # diverging cmap for per-trial correlograms (complements inferno)
XCORR_ALL_TRIALS_NORMALIZE = "zscore" # per-trial normalization: "zscore", "peak", or None
XCORR_GRAND_MEAN_SEM_ALPHA = 0.25    # fill_between alpha for SEM shading
FIGSIZE_XCORR_INCH = (28, 14)        # wider figure for 4+3 panel layout

# Unsplit cross-correlation (no REST/RUN separation) + circular-shift surrogate test
XCORR_UNSPLIT_N_SURROGATES = 500     # number of circular-shift surrogates per trial
XCORR_UNSPLIT_MIN_SHIFT_CYCLES = 2   # min circular shift in theta cycles (avoids trivial shifts)
FIGSIZE_XCORR_UNSPLIT_INCH = (24, 14)

# Figure 5: Bilateral spectral analysis (PSD, coherence, coherogram)
SPEC5_PSD_WINDOW_SEC = 1.0           # Welch window length (≈ MATLAB pwelch default)
SPEC5_PSD_OVERLAP_FRAC = 0.5         # 50% overlap
SPEC5_PSD_NFFT_MULT = 2              # zero-padding factor
SPEC5_COHEROGRAM_WINDOW_SEC = 2.0    # sliding-window coherence window
SPEC5_COHEROGRAM_STEP_SEC = 0.25     # step size for coherogram
SPEC5_COHEROGRAM_CMAP = "inferno"    # time-resolved coherence colormap
SPEC5_PSD_FREQ_RANGE = (1, 50)       # display range for PSD/coherence plots
SPEC5_1F_LOW_FLANK = (2, 4)          # Hz — below theta, for 1/f linear fit
SPEC5_1F_HIGH_FLANK = (15, 25)       # Hz — above alpha, for 1/f linear fit
SPEC5_THETA_EXTRACTION_METHOD = "peak"  # 'peak', 'mean', or 'auc' for coherence violin
SPEC5_AGGREGATE_ALL_SESSIONS = True   # pool all sessions for rows 2-4
SPEC5_SEM_ALPHA = 0.25               # fill_between alpha for SEM shading
SPEC5_PSD_LINEWIDTH = 2.5            # PSD / coherence curve linewidth
SPEC5_EDGE_TRIM_SEC = 0.5            # trim trial edges before PSD/coherence
SPEC5_MIN_SEGMENT_SAMPLES = 64       # minimum REST/RUN segment for Welch
# Dusty coral family for synchrony measures (coherence / cross-corr / phase)
COLOR_SYNCHRONY_REST = np.array([0.82, 0.52, 0.48])   # same as COLOR_PHASE_ROSE_REST
COLOR_SYNCHRONY_RUN = np.array([0.55, 0.30, 0.27])    # same as COLOR_PHASE_ROSE_RUN
FIGSIZE_SPEC5_INCH = (32, 22)        # large 4-row figure

# REST/RUN behaviour — align with ../../spectral_analysis/config/analysis_config.m
# and core/classify_behavior.m (used by spectral_analysis.m).
BEHAVIOR_CLASSIFICATION_MODE = "clear"  # "standard" or "clear"
BEHAVIOR_RUN_THRESHOLD_CMS = 2.0
BEHAVIOR_REST_THRESHOLD_CMS = 0.1
BEHAVIOR_MIN_BOUT_DURATION_SEC = 2.0  # standard mode: merge bouts shorter than this
BEHAVIOR_MIN_RUN_BOUT_SEC = 0.3       # clear mode: min RUN bout
BEHAVIOR_MIN_REST_BOUT_SEC = 0.3      # clear mode: min REST bout
BEHAVIOR_APPLY_BOUT_FILTER = True
BEHAVIOR_MOTION_SMOOTH_SAMPLES = 10   # movmean window (samples at imaging fs)

# Spectrogram (matched to MATLAB fiber pipeline)
SPEC_WINDOW_SEC = 0.75
SPEC_OVERLAP_FRAC = 0.90
SPEC_NFFT_MULT = 2
SPEC_SMOOTH_FREQ = 1
SPEC_SMOOTH_TIME = 5
FREQ_RANGE = (2, 70)

# Wheel / motion conversion
WHEEL_DIAMETER_CM = 19.0
ENCODER_COUNTS_PER_REV = 1024
EPHYS_SAMPLING_RATE = 30000
DISTANCE_PER_EDGE_CM = (np.pi * WHEEL_DIAMETER_CM) / ENCODER_COUNTS_PER_REV
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM

# Figure styling
DPI = 300
FONT_SIZE_SUPTITLE = 22
FONT_SIZE_TITLE = 20
FONT_SIZE_LABEL = 18
FONT_SIZE_TICK = 16
FONT_SIZE_SCALEBAR = 14
FONT_SIZE_COLORBAR = 16
FONT_SIZE_LEGEND = 14
AXIS_LINEWIDTH = 2.0
TICK_WIDTH = 1.8
TICK_LENGTH = 6
LW_RAW = 0.8
LW_BAND = 0.8
SCALEBAR_LINEWIDTH = 3.5

# Colors
COLOR_RIGHT_HP = np.array([0.127568, 0.566949, 0.550556])   # teal (fiber 1)
# Left HP fiber: rusty terracotta (redder / less yellow than motion orange below)
COLOR_LEFT_HP = np.array([0.78, 0.34, 0.28])
# Distinct purple-grey shades for ipsi (right) vs Ch11/contra (left) LFP
COLOR_LFP_RIGHT = np.array([0.35, 0.25, 0.45])               # purple-grey (ipsi / right)
COLOR_LFP_LEFT = np.array([0.22, 0.32, 0.58])                # bluish-purple (Ch11 / left)
COLOR_MOTION = np.array([0.993248, 0.7, 0.4])                # orange
# Cross-correlogram traces: dark rusty coral between right (teal) and left (rust) fibers
COLOR_XCORR_RUST = np.array([0.55, 0.26, 0.22])
LINE_WIDTH_XCORR = 2.5  # match fig2_coherence coherence traces

# Raw-trace figure: show only this time window (seconds); bandpass uses full trial then slice
TRACES_TIME_T0_S = 8.0
TRACES_TIME_T1_S = 13.0
# Wide layout: large figure + margins so right-column scale bars are not clipped
FIGSIZE_TRACES_INCH = (30.0, 11.5)  # (width, height) — wide panels
TRACES_SUBPLOTS_WSPACE = 0.08      # horizontal gap between columns (smaller → wider panels)
# Leave extra space outside the right column for scale bars (x > 1 in axes coords)
TRACES_FIG_LEFT = 0.055
TRACES_FIG_RIGHT = 0.82
TRACES_FIG_TOP = 0.90
TRACES_FIG_BOTTOM = 0.065
TRACES_SAVE_PAD_INCH = 0.35        # bbox_inches='tight' padding for scale-bar labels
# Right-column scale bars: stay just outside axis edge (smaller x → less clipping)
TRACES_SCALEBAR_X_RIGHT = 1.01

# Phase-locking / Hilbert (Figure 4): aggregate across all trials in session
# θ band is THETA_BAND (5–9 Hz); prefilter removes HF before θ bandpass for cleaner Hilbert phase.
PHASE_LOWPASS_HZ = 50.0
PHASE_LOWPASS_ORDER = 4
PHASE_PREFILTER_MOVMEAN_SAMPLES = 5  # 1 = off; light temporal smooth after low-pass
PHASE_EDGE_TRIM_SEC = 0.5  # trim trial ends when pooling (reduces filter/Hilbert edge effects)
PHASE_ROSE_NBINS = 24  # wider bins → smoother rose (e.g. 20–24)
PHASE_MEAN_VECTOR_LW = 5.2
PHASE_MEAN_VECTOR_ARROW_MUTATION = 32  # arrowhead size (matplotlib mutation_scale)
# Polar rose publication styling (frame, grid, typography)
# Slightly conservative sizes/weights: very large bold + Agg PNG on Windows can hit FreeType raster overflow (0x62).
PHASE_ROSE_TITLE_FONTSIZE = 18
PHASE_ROSE_TITLE_PAD = 14
PHASE_ROSE_TITLE_FONTWEIGHT = "normal"
PHASE_ROSE_THETA_TICK_FONTSIZE = 18
PHASE_ROSE_POLAR_FRAME_LW = 2.5  # outer circular border
PHASE_ROSE_POLAR_FRAME_COLOR = "0.1"
PHASE_ROSE_GRID_LW = 1.0
PHASE_ROSE_GRID_ALPHA = 0.38
PHASE_ROSE_GRID_COLOR = "0.38"
PHASE_ROSE_BAR_EDGELW = 0.55
PHASE_ROSE_THETA_GRID_DEG = (0, 90, 180, 270)  # azimuth lines + labels
# RUN rose: keep samples with θ-envelope >= median within each contiguous RUN bout (top 50%)
PHASE_RUN_USE_ENVELOPE_GATE = True
# Figure 4: pool trials from every session in RECORDINGS[animal] (e.g. all 18 Animal02 trials)
PHASE_AGGREGATE_ALL_SESSIONS = True
# Per-trial inclusion for trial-equalized REST / RUN stats (no duration domination)
PHASE_TRIAL_MIN_REST_SAMPLES = 30
PHASE_TRIAL_MIN_RUN_BOUTS = 1       # contiguous gated-RUN segments in m_run
PHASE_TRIAL_MIN_RUN_SAMPLES = 40    # gated RUN samples; else trial omitted from RUN rose/stats
# When True, collect_trialwise_phase_difference prints each trial's R_rest, R_run and inclusion flags
PHASE_LOG_PER_TRIAL_STATS = True
# Example traces / sawtooth: 1 s window; time scale bar sized for short windows
PHASE_EXAMPLE_T0_S = 9.4
PHASE_EXAMPLE_T1_S = 10.2
PHASE_SCALEBAR_T_SECONDS = 0.25  # horizontal bar on θ-trace panel (not on phase panel)
FIGSIZE_PHASE_LOCKING_INCH = (19.5, 15.0)  # row0: REST rose | RUN rose | R violin
# Rose bars + R-panel violins (dusty coral family; RUN darker)
COLOR_PHASE_ROSE_REST = np.array([0.82, 0.52, 0.48])
COLOR_PHASE_ROSE_RUN = np.array([0.55, 0.30, 0.27])
# Half-violin layout (match fig5_theta / FOOOF group plots)
PHASE_VIOLIN_WIDTH = 0.28
PHASE_VIOLIN_BOX_WIDTH = 0.10
PHASE_VIOLIN_BOX_OFFSET = 0.08
PHASE_VIOLIN_DOT_OFFSET = 0.18
PHASE_VIOLIN_DOT_SIZE = 90
PHASE_VIOLIN_LINE_WIDTH = 1.5
PHASE_VIOLIN_LINE_ALPHA = 0.5
# Figure 4 PNG: lower dpi avoids rare Agg/FreeType glyph bitmap overflow; PDF keeps full DPI.
PHASE_FIG4_PNG_DPI = 200
PHASE_FIG4_PNG_DPI_FALLBACK = 120
PHASE_FIG4_SUPTITLE_FONTSIZE = 18
PHASE_FIG4_SUPTITLE_FONTWEIGHT = "normal"

# Figure 6: Fiber–LFP theta cross-correlation (4 ipsi/contra combinations)
FIG6_ANIMAL = "Animal01"
FIG6_EXAMPLE_SESSION = "02_01_26-R2"
FIG6_EXAMPLE_TRIAL = 1
FIG6_AGGREGATE_ALL_SESSIONS = True
FIGSIZE_FIG6_INCH = (34, 14)
FIG6_COMBINATIONS = [
    {"sig1_key": "fiber1", "sig2_key": "lfp_right",
     "label": "R Fiber \u2013 R LFP (ipsi)", "short": "R\u2013R"},
    {"sig1_key": "fiber1", "sig2_key": "lfp_left",
     "label": "R Fiber \u2013 L LFP (contra)", "short": "R\u2013L"},
    {"sig1_key": "fiber2", "sig2_key": "lfp_left",
     "label": "L Fiber \u2013 L LFP (ipsi)", "short": "L\u2013L"},
    {"sig1_key": "fiber2", "sig2_key": "lfp_right",
     "label": "L Fiber \u2013 R LFP (contra)", "short": "L\u2013R"},
]
FIG6_COLORS = {
    "R\u2013R": np.array([0.127568, 0.566949, 0.550556]),   # teal (right HP)
    "R\u2013L": np.array([0.42, 0.42, 0.65]),               # muted indigo (contra)
    "L\u2013L": np.array([0.78, 0.34, 0.28]),               # terracotta (left HP)
    "L\u2013R": np.array([0.58, 0.42, 0.58]),               # muted plum (contra)
}
FIG6_HEATMAP_CMAP = "inferno"
FIG6_SEM_ALPHA = 0.25
FIG6_CORR_LINE_WIDTH = 2.5

# Figure 7: Fiber–LFP phase-locking (4 combinations, all epochs)
FIG7_ANIMAL = "Animal02"
FIG7_AGGREGATE_ALL_SESSIONS = True
FIGSIZE_FIG7_INCH = (34, 10)
FIG7_ROSE_BAR_COLOR = "0.72"
FIG7_ROSE_BAR_ALPHA = 0.58

# Figure 6b: Fiber–LFP xcorr REST vs RUN (4 combos, behaviour-split)
FIG6B_ANIMAL = "Animal02"
FIG6B_AGGREGATE_ALL_SESSIONS = True
FIGSIZE_FIG6B_INCH = (30, 12)

# Figure 7b: Fiber–LFP phase-locking REST vs RUN (4 combos, behaviour-split)
FIG7B_ANIMAL = "Animal02"
FIG7B_AGGREGATE_ALL_SESSIONS = True
FIGSIZE_FIG7B_INCH = (34, 16)

# LFP–LFP coherence diagnostic
FIG_LFP_DIAG_ANIMAL = "Animal02"
FIG_LFP_DIAG_EXAMPLE_SESSION = "01_01_26-R1"
FIG_LFP_DIAG_EXAMPLE_TRIAL = 1
FIG_LFP_DIAG_AGGREGATE_ALL_SESSIONS = True
FIGSIZE_LFP_DIAG_INCH = (18, 7)
LFP_DIAG_FREQ_RANGE = (1, 80)
COLOR_LFP_COH = np.array([0.45, 0.28, 0.55])  # purple for LFP–LFP coherence

# Scale bars
SCALEBAR_T_SECONDS = 1.0
SCALEBAR_GEVI_PERCENT = 1.0
SCALEBAR_LFP_UV = 500.0
SCALEBAR_MOTION_CM_S = 15.0
SCALEBAR_GEVI_BAND_PERCENT = 0.5
SCALEBAR_LFP_BAND_UV = 200.0

# Speed colormap (monochromatic orange, matching run_all_plots.py)
def _create_speed_cmap():
    colors = [
        (0.25, 0.06, 0.00),
        (0.45, 0.12, 0.02),
        (0.65, 0.22, 0.04),
        (0.85, 0.38, 0.08),
        (0.98, 0.55, 0.18),
        (1.00, 0.75, 0.45),
        (1.00, 0.95, 0.82),
    ]
    return LinearSegmentedColormap.from_list("mono_orange", colors, N=256)

CMAP_SPEED = _create_speed_cmap()


def _apply_publication_rcparams():
    """Vector-friendly PDF text + white figures (aligned with other pipeline plot scripts)."""
    mpl.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


_apply_publication_rcparams()


# =============================================================================
# DATA LOADING
# =============================================================================

def _normalize_unc(path_str):
    """Prepend \\\\?\\UNC\\ for long UNC paths on Windows."""
    s = str(path_str)
    if s.startswith("\\\\") and not s.startswith("\\\\?\\"):
        return "\\\\?\\UNC\\" + s[2:]
    return s


def _trial_mat_path(animal, session, trial_num):
    """
    Resolve the .mat path for a trial. Two on-disk layouts exist:
      (a) session / Trial{N}_suffix_{N} / filename.mat   (older sessions)
      (b) session / filename.mat                          (newer sessions)
    Try (a) first; fall back to (b) if the subfolder doesn't exist.
    """
    info = RECORDINGS[animal][session]
    suffix = info["suffix"]
    filename = f"{animal}-{session}_Trial{trial_num}_FiberPhotometry_Analysis.mat"
    session_dir = BASE_PATH / animal / "Fiber_Voltage_Processed" / session

    subfolder_path = session_dir / f"Trial{trial_num}_{suffix}_{trial_num}" / filename
    if Path(_normalize_unc(subfolder_path)).exists():
        return subfolder_path

    flat_path = session_dir / filename
    if Path(_normalize_unc(flat_path)).exists():
        return flat_path

    return subfolder_path


def _hdf5_array(grp, key):
    """Extract a numpy float64 array from an HDF5 group, handling common shapes."""
    return np.array(grp[key][()], dtype=np.float64).squeeze()


@lru_cache(maxsize=None)
def load_trial(animal, session, trial_num):
    """
    Load one trial from preprocessed .mat file.
    Returns dict with keys: fiber1, fiber2, lfp_right, lfp_left, speed, time, fs.
    Tries h5py first (MATLAB v7.3), falls back to scipy.io.loadmat.

    Cached per (animal, session, trial_num): the same trial's file is loaded
    repeatedly across the different figure-generating sections in main()
    (up to ~4x per animal), and each load re-reads and re-parses the same
    .mat file from disk. The returned dict's array values are not mutated
    in place anywhere in this module, so sharing the cached dict across
    callers is safe.
    """
    mat_path = _trial_mat_path(animal, session, trial_num)
    access_path = _normalize_unc(mat_path)
    print(f"  Loading: {mat_path.name}")

    # ---- Try HDF5 (v7.3) first ----
    try:
        return _load_hdf5(access_path)
    except Exception:
        pass

    # ---- Fallback: scipy.io.loadmat ----
    try:
        return _load_scipy(access_path)
    except Exception as e:
        raise RuntimeError(f"Cannot load {mat_path}: {e}") from e


def _load_hdf5(path_str):
    with h5py.File(path_str, "r") as f:
        fpa = f["FiberPhotometryAnalysis"]

        fs = float(_hdf5_array(fpa["time"], "sampling_rate"))
        t = _hdf5_array(fpa["time"], "time_vector_seconds").reshape(-1)

        sig = fpa["signals"]
        for trace_key in ("final_processed_traces", "deltaF_F_traces"):
            if trace_key in sig:
                traces = _hdf5_array(sig, trace_key)
                break
        else:
            raise KeyError("No fiber traces found in signals group")

        if traces.ndim == 1:
            traces = traces.reshape(-1, 1)
        n_time = t.size
        if traces.shape[0] != n_time and traces.shape[1] == n_time:
            traces = traces.T
        fiber1 = traces[:, 0]
        fiber2 = traces[:, 1] if traces.shape[1] > 1 else np.full_like(fiber1, np.nan)

        eph = fpa["ephys"]

        lfp_right = _try_keys(eph, ("lfp_raw_aligned_ipsiHP", "lfp_z_ipsiHP"))
        lfp_left = _try_keys(eph, ("lfp_raw_aligned_HP", "lfp_z_HP"))
        speed_raw = _try_keys(eph, ("running_velocity_smooth", "running_velocity"))

    return _finalise(t, fiber1, fiber2, lfp_right, lfp_left, speed_raw, fs)


def _try_keys(grp, keys):
    for key in keys:
        if key in grp:
            try:
                return _hdf5_array(grp, key).reshape(-1)
            except Exception:
                continue
    return None


def _load_scipy(path_str):
    mat = loadmat(path_str, squeeze_me=True, struct_as_record=False)
    fpa = mat["FiberPhotometryAnalysis"]

    fs = float(fpa.time.sampling_rate)
    t = np.array(fpa.time.time_vector_seconds).reshape(-1)

    for attr in ("final_processed_traces", "deltaF_F_traces"):
        if hasattr(fpa.signals, attr):
            traces = np.atleast_2d(getattr(fpa.signals, attr))
            break
    else:
        raise KeyError("No fiber traces found")
    if traces.shape[0] < traces.shape[1]:
        traces = traces.T
    fiber1 = traces[:, 0]
    fiber2 = traces[:, 1] if traces.shape[1] > 1 else np.full_like(fiber1, np.nan)

    lfp_right = _try_attr(fpa.ephys, ("lfp_raw_aligned_ipsiHP", "lfp_z_ipsiHP"))
    lfp_left = _try_attr(fpa.ephys, ("lfp_raw_aligned_HP", "lfp_z_HP"))
    speed_raw = _try_attr(fpa.ephys, ("running_velocity_smooth", "running_velocity"))

    return _finalise(t, fiber1, fiber2, lfp_right, lfp_left, speed_raw, fs)


def _try_attr(obj, attrs):
    for a in attrs:
        if hasattr(obj, a):
            try:
                return np.array(getattr(obj, a), dtype=np.float64).reshape(-1)
            except Exception:
                continue
    return None


def _finalise(t, fiber1, fiber2, lfp_right, lfp_left, speed_raw, fs):
    n = len(t)
    fiber1 = fiber1[:n]
    fiber2 = fiber2[:n]
    if lfp_right is not None:
        lfp_right = lfp_right[:n]
    if lfp_left is not None:
        lfp_left = lfp_left[:n]
    if speed_raw is not None:
        speed = np.abs(speed_raw[:n]) * MOTION_TO_CM_PER_S
    else:
        speed = None

    return {
        "fiber1": fiber1, "fiber2": fiber2,
        "lfp_right": lfp_right, "lfp_left": lfp_left,
        "speed": speed, "time": t, "fs": fs,
    }


# =============================================================================
# SIGNAL PROCESSING HELPERS
# =============================================================================

def bandpass(sig, fs, lo, hi, order=FILTER_ORDER):
    nyq = fs / 2
    lo_n = max(0.01, lo) / nyq
    hi_n = min(0.95 * nyq, hi) / nyq
    b, a = signal.butter(order, [lo_n, hi_n], btype="band")
    out = sig - np.nanmean(sig)
    out = np.nan_to_num(out, nan=0.0)
    return signal.filtfilt(b, a, out)


def lowpass(sig, fs, cutoff_hz, order=4):
    """
    Butterworth low-pass (zero-phase filtfilt). Use before θ bandpass to attenuate
    high-frequency noise that degrades Hilbert phase (e.g. < 50 Hz at imaging fs).
    """
    nyq = fs / 2.0
    co = float(cutoff_hz)
    if co >= 0.95 * nyq:
        co = max(0.05 * nyq, min(co, 0.9 * nyq))
    wn = co / nyq
    wn = float(np.clip(wn, 0.01, 0.99))
    b, a = signal.butter(int(order), wn, btype="low")
    out = np.asarray(sig, dtype=float) - np.nanmean(sig)
    out = np.nan_to_num(out, nan=0.0)
    return signal.filtfilt(b, a, out)


def _smooth_speed_movmean(speed, win):
    """Moving average for behaviour classification (cf. classify_behavior.m)."""
    s = np.asarray(speed, dtype=float).ravel()
    if win <= 1:
        return s
    return uniform_filter1d(s, size=win, mode="nearest")


def _merge_short_bouts(is_running, min_samples):
    """
    Standard-mode bout merging (classify_behavior.m merge_short_bouts).
    Short RUN bouts become REST and vice versa by adopting the surrounding state.
    """
    x = np.asarray(is_running, dtype=bool).astype(np.int8).ravel()
    n = len(x)
    if n == 0 or min_samples <= 1:
        return x.astype(bool)
    out = x.copy()
    i = 0
    while i < n:
        j = i + 1
        while j < n and x[j] == x[i]:
            j += 1
        bout_len = j - i
        if bout_len < min_samples:
            if i > 0:
                surround = int(out[i - 1])
            elif j < n:
                surround = int(out[j])
            else:
                i = j
                continue
            out[i:j] = surround
        i = j
    return out.astype(bool)


def _filter_short_bouts_exclude(state, min_samples):
    """Clear mode: remove (set False) connected components shorter than min_samples."""
    state = np.asarray(state, dtype=bool).copy().ravel()
    if min_samples <= 1:
        return state
    labeled, nfeat = ndimage_label(state)
    out = state.copy()
    for k in range(1, nfeat + 1):
        mask = labeled == k
        if np.sum(mask) < min_samples:
            out[mask] = False
    return out


def classify_rest_run_masks(speed_cm_s, fs):
    """
    REST/RUN masks aligned with ../../spectral_analysis/core/classify_behavior.m
    and analysis_config.m defaults.
    Returns is_rest, is_run, is_excluded (excluded only in 'clear' mode).
    """
    speed_s = _smooth_speed_movmean(speed_cm_s, BEHAVIOR_MOTION_SMOOTH_SAMPLES)
    mode = str(BEHAVIOR_CLASSIFICATION_MODE).lower()
    min_bout = max(1, int(round(BEHAVIOR_MIN_BOUT_DURATION_SEC * fs)))
    min_run = max(1, int(round(BEHAVIOR_MIN_RUN_BOUT_SEC * fs)))
    min_rest = max(1, int(round(BEHAVIOR_MIN_REST_BOUT_SEC * fs)))

    if mode == "standard":
        is_run_raw = speed_s > BEHAVIOR_RUN_THRESHOLD_CMS
        if BEHAVIOR_APPLY_BOUT_FILTER and min_bout > 1:
            is_run = _merge_short_bouts(is_run_raw, min_bout)
        else:
            is_run = is_run_raw
        is_rest = ~is_run
        is_excluded = np.zeros_like(is_run, dtype=bool)
    elif mode == "clear":
        is_run_raw = speed_s > BEHAVIOR_RUN_THRESHOLD_CMS
        is_rest_raw = speed_s < BEHAVIOR_REST_THRESHOLD_CMS
        if BEHAVIOR_APPLY_BOUT_FILTER:
            is_run = _filter_short_bouts_exclude(is_run_raw, min_run)
            is_rest = _filter_short_bouts_exclude(is_rest_raw, min_rest)
        else:
            is_run, is_rest = is_run_raw, is_rest_raw
        is_excluded = ~is_run & ~is_rest
    else:
        raise ValueError(f"Unknown BEHAVIOR_CLASSIFICATION_MODE: {BEHAVIOR_CLASSIFICATION_MODE}")
    return is_rest, is_run, is_excluded


def _time_index_at_sec(time_vec, t_sec):
    """Map seconds to index into speed / behaviour arrays."""
    tv = np.asarray(time_vec, dtype=float).ravel()
    if tv.size == 0:
        return 0
    idx = int(np.searchsorted(tv, t_sec, side="right") - 1)
    return int(np.clip(idx, 0, tv.size - 1))


from common import _next_pow2  # shared helpers (were local copies)


def _smooth2d(arr, nf=1, nt=1):
    out = arr.copy()
    if nf > 1:
        out = uniform_filter1d(out, size=nf, axis=0, mode="nearest")
    if nt > 1:
        out = uniform_filter1d(out, size=nt, axis=1, mode="nearest")
    return out


def compute_spectrogram(trace, fs):
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


def wrap_to_pi(angles_rad):
    """Wrap angles to (−π, π]."""
    a = np.asarray(angles_rad, dtype=float)
    return (a + np.pi) % (2.0 * np.pi) - np.pi


def circular_mean_and_resultant_length(angles_rad):
    """
    Circular mean direction and mean resultant length R = |mean(exp(iφ))|.
    Returns (mean_phi_rad, R). Empty / all-NaN → (nan, nan).
    """
    x = np.asarray(angles_rad, dtype=float).ravel()
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    z = np.mean(np.exp(1j * x))
    R = float(np.abs(z))
    mean_phi = float(np.angle(z))
    return mean_phi, R


def preprocess_frac_for_theta(delta_f_frac, fs):
    """
    ΔF/F (fraction) → demean, % scale, low-pass (PHASE_LOWPASS_HZ), optional movmean.
    Output is wideband % signal fed into θ bandpass + Hilbert for Figure 4.
    """
    x = np.asarray(delta_f_frac, dtype=float).ravel()
    x = x - np.nanmean(x)
    x = np.nan_to_num(x, nan=0.0)
    x_pct = x * 100.0
    x_lp = lowpass(x_pct, fs, PHASE_LOWPASS_HZ, order=PHASE_LOWPASS_ORDER)
    win = int(PHASE_PREFILTER_MOVMEAN_SAMPLES)
    if win > 1:
        x_lp = uniform_filter1d(x_lp, size=win, mode="nearest")
    return x_lp


def theta_analytic_signal(preprocessed_pct, fs, band=THETA_BAND):
    """θ-band analytic signal (complex) from prefiltered % trace."""
    bp = bandpass(preprocessed_pct, fs, band[0], band[1])
    bp = np.nan_to_num(bp, nan=0.0)
    return signal.hilbert(bp)


def theta_analytic_from_raw(sig, fs, band=THETA_BAND):
    """
    Generic: raw signal → demean → low-pass → θ bandpass → Hilbert analytic signal.
    Works for both fiber ΔF/F and LFP signals (phase is scale-invariant).
    """
    x = np.asarray(sig, dtype=float).ravel()
    x = x - np.nanmean(x)
    x = np.nan_to_num(x, nan=0.0)
    x_lp = lowpass(x, fs, PHASE_LOWPASS_HZ, order=PHASE_LOWPASS_ORDER)
    win = int(PHASE_PREFILTER_MOVMEAN_SAMPLES)
    if win > 1:
        x_lp = uniform_filter1d(x_lp, size=win, mode="nearest")
    return theta_analytic_signal(x_lp, fs, band=band)


def instantaneous_theta_phase_hilbert(delta_f_frac, fs, band=THETA_BAND):
    """
    Low-pass (+ optional movmean) → θ-band (THETA_BAND) → Hilbert phase [−π, π].
    `delta_f_frac` unitless ΔF/F (same as struct), not ×100.
    """
    pre = preprocess_frac_for_theta(delta_f_frac, fs)
    z = theta_analytic_signal(pre, fs, band=band)
    return np.angle(z)


def phase_difference_left_minus_right(phi_left, phi_right):
    """Δφ = φ_LHP − φ_RHP, wrapped to [−π, π]."""
    return wrap_to_pi(np.asarray(phi_left, dtype=float) - np.asarray(phi_right, dtype=float))


def _envelope_top_half_per_run_bout(is_run, envelope):
    """
    Within each contiguous RUN bout, keep samples with envelope >= bout median
    (top 50% of envelope magnitude in that bout).
    """
    is_run = np.asarray(is_run, dtype=bool).ravel()
    env = np.asarray(envelope, dtype=float).ravel()
    out = np.zeros_like(is_run, dtype=bool)
    if is_run.size == 0 or env.size != is_run.size:
        return out
    labeled, nfeat = ndimage_label(is_run)
    for k in range(1, nfeat + 1):
        bm = labeled == k
        ev = env[bm]
        if ev.size == 0:
            continue
        thr = float(np.median(ev))
        out[bm] = env[bm] >= thr
    return out


def _n_contiguous_true_regions(mask_bool):
    """Number of separated True runs in a 1D mask (e.g. gated-RUN islands)."""
    m = np.asarray(mask_bool, dtype=bool).ravel()
    if m.size == 0 or not np.any(m):
        return 0
    _, nfeat = ndimage_label(m)
    return int(nfeat)


def collect_trialwise_phase_difference(
    animal,
    session=None,
    band=THETA_BAND,
    edge_trim_sec=PHASE_EDGE_TRIM_SEC,
):
    """
    One entry per trial: Δφ = φ_LHP − φ_RHP samples for REST and RUN (same preprocessing
    as Figure 4). `session=None` iterates every session in RECORDINGS[animal].

    Flags:
      use_for_rest_stats — trial has ≥ PHASE_TRIAL_MIN_REST_SAMPLES in REST.
      use_for_run_stats   — trial has ≥ PHASE_TRIAL_MIN_RUN_BOUTS contiguous gated-RUN
                            regions AND ≥ PHASE_TRIAL_MIN_RUN_SAMPLES gated RUN samples.

    Trial-equalized roses average each trial's *normalized* circular histogram and
    mean complex vector z_t = mean(exp(iΔφ)) across qualifying trials (equal weight).
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return {
            "trials": [],
            "n_loaded_ok": 0,
            "n_failed": 0,
            "n_rest_qualifying": 0,
            "n_run_qualifying": 0,
            "note": "unknown animal",
        }

    if session is not None:
        if session not in sessions_map:
            return {
                "trials": [],
                "n_loaded_ok": 0,
                "n_failed": 0,
                "n_rest_qualifying": 0,
                "n_run_qualifying": 0,
                "note": f"unknown session {session!r}",
            }
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    trials_out = []
    failed = 0
    note_parts = []

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception as e:
                failed += 1
                note_parts.append(f"{sess_key} T{trial_num}: load ({e})")
                continue

            f1 = np.asarray(data["fiber1"], dtype=float).ravel()
            f2 = np.asarray(data["fiber2"], dtype=float).ravel()
            if f2.size == 0 or not np.any(np.isfinite(f2)):
                failed += 1
                note_parts.append(f"{sess_key} T{trial_num}: no left HP")
                continue

            fs = float(data["fs"])
            speed = data["speed"]
            if speed is None:
                is_rest = np.ones(len(f1), dtype=bool)
                is_run = np.zeros_like(is_rest)
                is_excluded = np.zeros_like(is_rest)
            else:
                is_rest, is_run, is_excluded = classify_rest_run_masks(
                    np.asarray(speed, dtype=float).ravel(), fs
                )

            n = min(f1.size, f2.size, is_rest.size)
            f1, f2 = f1[:n], f2[:n]
            is_rest = is_rest[:n]
            is_run = is_run[:n]
            is_excluded = is_excluded[:n]

            trim = max(0, int(round(edge_trim_sec * fs)))
            if n <= 2 * trim + 8:
                failed += 1
                note_parts.append(f"{sess_key} T{trial_num}: too short")
                continue

            sl = slice(trim, n - trim)
            f1, f2 = f1[sl], f2[sl]
            is_rest = is_rest[sl]
            is_run = is_run[sl]
            is_excluded = is_excluded[sl]

            pre1 = preprocess_frac_for_theta(f1, fs)
            pre2 = preprocess_frac_for_theta(f2, fs)
            z_r = theta_analytic_signal(pre1, fs, band=band)
            z_l = theta_analytic_signal(pre2, fs, band=band)
            phi_r = np.angle(z_r)
            phi_l = np.angle(z_l)
            dphi = phase_difference_left_minus_right(phi_l, phi_r)
            env_mean = 0.5 * (np.abs(z_r) + np.abs(z_l))

            ok = np.isfinite(dphi) & (~is_excluded)
            m_rest = ok & is_rest
            if PHASE_RUN_USE_ENVELOPE_GATE:
                env_gate = _envelope_top_half_per_run_bout(is_run, env_mean)
                m_run = ok & is_run & env_gate
            else:
                m_run = ok & is_run

            phi_r_arr = dphi[m_rest]
            phi_n_arr = dphi[m_run]
            n_rest_s = int(phi_r_arr.size)
            n_run_s = int(phi_n_arr.size)
            n_run_bouts = _n_contiguous_true_regions(m_run)

            use_rest = n_rest_s >= PHASE_TRIAL_MIN_REST_SAMPLES
            use_run = (
                n_run_bouts >= PHASE_TRIAL_MIN_RUN_BOUTS
                and n_run_s >= PHASE_TRIAL_MIN_RUN_SAMPLES
            )

            # Per-trial circular stats (for tables / debugging)
            mu_r, r_r = circular_mean_and_resultant_length(phi_r_arr)
            mu_n, r_n = circular_mean_and_resultant_length(phi_n_arr)

            trials_out.append(
                {
                    "session": sess_key,
                    "trial_num": trial_num,
                    "phi_rest": phi_r_arr,
                    "phi_run": phi_n_arr,
                    "n_rest_samples": n_rest_s,
                    "n_run_samples": n_run_s,
                    "n_run_bouts": n_run_bouts,
                    "use_for_rest_stats": use_rest,
                    "use_for_run_stats": use_run,
                    "R_rest_trial": r_r,
                    "R_run_trial": r_n,
                    "mean_phi_rest_trial": mu_r,
                    "mean_phi_run_trial": mu_n,
                }
            )

    gate_txt = (
        "RUN gated: theta-env top 50%/bout"
        if PHASE_RUN_USE_ENVELOPE_GATE
        else "RUN: all run samples"
    )
    scope = "all sessions" if session is None else f"session {session}"
    note = (
        f"{scope}; loaded trials={len(trials_out)}, load_fail={failed}; "
        f"trial-equalized roses; REST min n>={PHASE_TRIAL_MIN_REST_SAMPLES}; "
        f"RUN min bouts>={PHASE_TRIAL_MIN_RUN_BOUTS}, min samples>={PHASE_TRIAL_MIN_RUN_SAMPLES}; "
        f"dphi=LHP-RHP; theta={band[0]}-{band[1]} Hz; LP<={PHASE_LOWPASS_HZ:.0f} Hz; "
        f"{gate_txt}; edge trim={edge_trim_sec}s"
    )
    if note_parts:
        note += " | " + "; ".join(note_parts[:6])
        if len(note_parts) > 6:
            note += "..."

    n_rest_qual = sum(1 for t in trials_out if t["use_for_rest_stats"])
    n_run_qual = sum(1 for t in trials_out if t["use_for_run_stats"])
    note += f" | qualifying: REST n_trials={n_rest_qual}, RUN n_trials={n_run_qual}"

    if PHASE_LOG_PER_TRIAL_STATS and trials_out:
        def _fmt_r(x):
            return f"{float(x):.4f}" if np.isfinite(x) else "nan"

        print(
            "\n  Per-trial dphi stats (circular R on within-trial samples; "
            "inclusion for trial-equalized roses):"
        )
        for t in trials_out:
            ur = "yes" if t["use_for_rest_stats"] else "no"
            un = "yes" if t["use_for_run_stats"] else "no"
            print(
                f"    {t['session']} trial {t['trial_num']}: "
                f"R_rest={_fmt_r(t['R_rest_trial'])} (REST incl. {ur}, n={t['n_rest_samples']}); "
                f"R_run={_fmt_r(t['R_run_trial'])} (RUN incl. {un}, n={t['n_run_samples']}, "
                f"bouts={t['n_run_bouts']})"
            )

    return {
        "trials": trials_out,
        "n_loaded_ok": len(trials_out),
        "n_failed": failed,
        "note": note,
        "n_rest_qualifying": n_rest_qual,
        "n_run_qualifying": n_run_qual,
    }


def _style_phase_rose_polar_axis(ax):
    """
    Publication-style polar axes: thick outer frame, azimuth grid only, legible θ labels.
    Call after y/r limits are set.
    """
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(PHASE_ROSE_POLAR_FRAME_LW)
        spine.set_edgecolor(PHASE_ROSE_POLAR_FRAME_COLOR)
        spine.set_zorder(15)
    ax.xaxis.grid(
        True,
        linestyle="-",
        linewidth=PHASE_ROSE_GRID_LW,
        color=PHASE_ROSE_GRID_COLOR,
        alpha=PHASE_ROSE_GRID_ALPHA,
        zorder=0.5,
    )
    ax.yaxis.grid(
        True,
        linestyle="-",
        linewidth=PHASE_ROSE_GRID_LW * 0.6,
        color=PHASE_ROSE_GRID_COLOR,
        alpha=PHASE_ROSE_GRID_ALPHA * 0.7,
        zorder=0.5,
    )
    ax.set_yticklabels([])
    ax.tick_params(
        axis="x",
        which="major",
        pad=8,
        colors=PHASE_ROSE_POLAR_FRAME_COLOR,
        width=TICK_WIDTH + 0.4,
        length=TICK_LENGTH + 3,
        labelsize=PHASE_ROSE_THETA_TICK_FONTSIZE,
    )
    deg = list(PHASE_ROSE_THETA_GRID_DEG)
    ax.set_thetagrids(
        deg,
        [f"${d}^\\circ$" for d in deg],
        fontsize=PHASE_ROSE_THETA_TICK_FONTSIZE,
        color=PHASE_ROSE_POLAR_FRAME_COLOR,
    )


def _trial_equalized_rose_histogram_and_vector(trial_phi_list, n_bins):
    """
    trial_phi_list: list of 1D arrays (one qualifying trial each).
    For each trial: normalized histogram (sums to 1); mean across trials → shape.
    Group R, μ: mean of per-trial complex mean vectors z_t = mean(exp(iΔφ)).

    Returns dict with theta, bar_heights (for polar radius), width, mu, R, n_trials, scale_max.
    """
    cleaned = []
    for a in trial_phi_list:
        if a is None:
            continue
        phi = np.asarray(a, dtype=float).ravel()
        phi = phi[np.isfinite(phi)]
        if phi.size > 0:
            cleaned.append(phi)
    n_tri = len(cleaned)
    if n_tri == 0:
        return None

    bins = np.linspace(-np.pi, np.pi, int(n_bins) + 1)
    theta = (bins[:-1] + bins[1:]) / 2.0
    width = 2.0 * np.pi / float(n_bins)
    h_norm = []
    z_list = []
    for phi in cleaned:
        c, _ = np.histogram(phi, bins=bins)
        s = float(np.sum(c))
        h_norm.append(c.astype(np.float64) / max(s, 1.0))
        z_list.append(np.mean(np.exp(1j * phi)))
    hist_mean = np.mean(h_norm, axis=0)
    Z = np.mean(z_list)
    mu = float(np.angle(Z))
    R = float(np.abs(Z))
    scale_max = float(np.max(hist_mean)) if hist_mean.size else 1.0
    if scale_max < 1e-15:
        scale_max = 1.0
    bar_h = hist_mean / scale_max * 10.0
    return {
        "theta": theta,
        "bar_heights": bar_h,
        "width": width,
        "mu": mu,
        "R": R,
        "n_trials": n_tri,
        "scale_max": scale_max,
    }


def _plot_phase_rose_polar_trial_equalized(
    ax,
    trial_phi_list,
    title,
    n_bins=PHASE_ROSE_NBINS,
    bar_color=None,
):
    """
    Rose from trial-mean normalized histograms; μ and R from mean of per-trial z-vectors.
    `bar_color`: RGB array or mpl color for histogram bars (default dusty coral REST tone).
    """
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    res = _trial_equalized_rose_histogram_and_vector(trial_phi_list, n_bins)
    if res is None:
        ax.set_title(
            title + "\n(no qualifying trials)",
            fontsize=PHASE_ROSE_TITLE_FONTSIZE,
            fontweight=PHASE_ROSE_TITLE_FONTWEIGHT,
            pad=PHASE_ROSE_TITLE_PAD,
        )
        _style_phase_rose_polar_axis(ax)
        return

    if bar_color is None:
        bar_color = tuple(np.clip(COLOR_PHASE_ROSE_REST, 0, 1))
    elif isinstance(bar_color, np.ndarray):
        bar_color = tuple(np.clip(bar_color, 0, 1))
    theta = res["theta"]
    bar_h = res["bar_heights"]
    width = res["width"]
    mu = res["mu"]
    R = res["R"]
    n_tri = res["n_trials"]
    r_max = float(np.max(bar_h)) if bar_h.size else 10.0

    ax.bar(
        theta,
        bar_h,
        width=width * 0.98,
        bottom=0.0,
        align="center",
        color=bar_color,
        edgecolor="white",
        linewidth=PHASE_ROSE_BAR_EDGELW,
        alpha=0.62,
    )

    if np.isfinite(mu) and np.isfinite(R):
        ax.annotate(
            "",
            xytext=(0.0, 0.0),
            xy=(mu, R * r_max),
            arrowprops=dict(
                arrowstyle="-|>",
                color="crimson",
                lw=PHASE_MEAN_VECTOR_LW,
                mutation_scale=PHASE_MEAN_VECTOR_ARROW_MUTATION,
                shrinkA=0,
                shrinkB=0,
            ),
            zorder=12,
        )

    ax.set_ylim(0.0, r_max * 1.12)
    _style_phase_rose_polar_axis(ax)
    deg = np.degrees(mu) if np.isfinite(mu) else float("nan")
    ttl = (
        f"{title}\n"
        f"n={n_tri} trials | phase={deg:.1f} deg | R={R:.3f}"
    )
    ax.set_title(
        ttl,
        fontsize=PHASE_ROSE_TITLE_FONTSIZE,
        fontweight=PHASE_ROSE_TITLE_FONTWEIGHT,
        pad=PHASE_ROSE_TITLE_PAD,
    )


def _paired_per_trial_R_rest_run(phase_data):
    """
    Trials that qualify for both REST and RUN: paired R_rest, R_run per trial.
    R = |⟨exp(iΔφ)⟩| within that trial (mean resultant length of phase difference).
    """
    r_rest, r_run = [], []
    for t in phase_data.get("trials", []):
        if not (t.get("use_for_rest_stats") and t.get("use_for_run_stats")):
            continue
        rr = t.get("R_rest_trial", np.nan)
        rn = t.get("R_run_trial", np.nan)
        if np.isfinite(rr) and np.isfinite(rn):
            r_rest.append(float(rr))
            r_run.append(float(rn))
    return np.asarray(r_rest, dtype=float), np.asarray(r_run, dtype=float)


def _phase_locking_plot_half_violin_R(ax, data_rest, data_run, positions, colors):
    """
    Half-violin + box + scatter + paired lines (same geometry as
    fig5_theta.plot_half_violin).
    """
    pos_left, pos_right = positions
    color_left = tuple(np.clip(colors[0], 0, 1))
    color_right = tuple(np.clip(colors[1], 0, 1))

    violin_base_left = pos_left - 0.05
    box_pos_left = pos_left + PHASE_VIOLIN_BOX_OFFSET
    dot_pos_left = pos_left + PHASE_VIOLIN_DOT_OFFSET

    violin_base_right = pos_right + 0.05
    box_pos_right = pos_right - PHASE_VIOLIN_BOX_OFFSET
    dot_pos_right = pos_right - PHASE_VIOLIN_DOT_OFFSET

    def _violin_path(data, position, side, width=PHASE_VIOLIN_WIDTH):
        data = np.asarray(data, dtype=float).ravel()
        if data.size < 2 or float(np.std(data)) < 1e-10:
            return None, None, None
        try:
            kde = gaussian_kde(data, bw_method=0.5)
            data_range = float(np.max(data) - np.min(data))
            padding = max(data_range * 0.2, float(np.std(data)) * 0.3)
            y_range = np.linspace(np.min(data) - padding, np.max(data) + padding, 100)
            density = kde(y_range)
            density = density / np.max(density) * width
            x_path = position - density if side == "left" else position + density
            return y_range, x_path, density
        except Exception:
            return None, None, None

    res_l = _violin_path(data_rest, violin_base_left, "left")
    if res_l[0] is not None:
        y_range, x_path, _ = res_l
        ax.fill_betweenx(y_range, violin_base_left, x_path, alpha=0.6, color=color_left)

    res_r = _violin_path(data_run, violin_base_right, "right")
    if res_r[0] is not None:
        y_range, x_path, _ = res_r
        ax.fill_betweenx(y_range, violin_base_right, x_path, alpha=0.6, color=color_right)

    if data_rest.size > 0:
        bp_l = ax.boxplot(
            [data_rest],
            positions=[box_pos_left],
            widths=PHASE_VIOLIN_BOX_WIDTH,
            patch_artist=True,
            showfliers=False,
            zorder=3,
        )
        for patch in bp_l["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor(color_left)
            patch.set_linewidth(2)
        for w in bp_l["whiskers"]:
            w.set_color(color_left)
            w.set_linewidth(2)
        for c in bp_l["caps"]:
            c.set_color(color_left)
            c.set_linewidth(2)
        for med in bp_l["medians"]:
            med.set_color("black")
            med.set_linewidth(2.5)

    if data_run.size > 0:
        bp_r = ax.boxplot(
            [data_run],
            positions=[box_pos_right],
            widths=PHASE_VIOLIN_BOX_WIDTH,
            patch_artist=True,
            showfliers=False,
            zorder=3,
        )
        for patch in bp_r["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor(color_right)
            patch.set_linewidth(2)
        for w in bp_r["whiskers"]:
            w.set_color(color_right)
            w.set_linewidth(2)
        for c in bp_r["caps"]:
            c.set_color(color_right)
            c.set_linewidth(2)
        for med in bp_r["medians"]:
            med.set_color("black")
            med.set_linewidth(2.5)

    if data_rest.size:
        ax.scatter(
            [dot_pos_left] * len(data_rest),
            data_rest,
            s=PHASE_VIOLIN_DOT_SIZE,
            c=[color_left],
            edgecolors="white",
            linewidths=1.5,
            zorder=5,
            alpha=0.9,
        )
    if data_run.size:
        ax.scatter(
            [dot_pos_right] * len(data_run),
            data_run,
            s=PHASE_VIOLIN_DOT_SIZE,
            c=[color_right],
            edgecolors="white",
            linewidths=1.5,
            zorder=5,
            alpha=0.9,
        )

    n_pair = min(len(data_rest), len(data_run))
    for i in range(n_pair):
        ax.plot(
            [dot_pos_left, dot_pos_right],
            [data_rest[i], data_run[i]],
            color="0.45",
            alpha=PHASE_VIOLIN_LINE_ALPHA,
            lw=PHASE_VIOLIN_LINE_WIDTH,
            zorder=4,
        )


def _perform_paired_test(rest_data, run_data):
    """
    Paired statistical test on REST vs RUN R-values (two-tailed).
    Shapiro-Wilk on differences → paired t-test if normal (or n >= 7), else Wilcoxon.
    Follows the pattern from stimulation_analysis.py / fig5_theta.py.
    """
    rest_data = np.asarray(rest_data, dtype=float)
    run_data = np.asarray(run_data, dtype=float)
    valid = ~(np.isnan(rest_data) | np.isnan(run_data))
    rest_data = rest_data[valid]
    run_data = run_data[valid]
    n = len(rest_data)
    if n < 3:
        return {"p_value": np.nan, "test_used": "insufficient_data", "n": n}

    differences = run_data - rest_data
    try:
        _, p_norm = shapiro(differences)
        is_normal = p_norm > 0.05
    except Exception:
        is_normal = False

    try:
        if is_normal or n >= 7:
            stat, p_value = ttest_rel(run_data, rest_data)
            test_used = "paired_ttest"
        else:
            stat, p_value = wilcoxon(run_data, rest_data, alternative="two-sided")
            test_used = "wilcoxon"
        mean_diff = float(np.mean(differences))
        std_diff = float(np.std(differences, ddof=1))
        effect_size = mean_diff / std_diff if std_diff > 0 else 0.0
        return {
            "p_value": float(p_value),
            "statistic": float(stat),
            "effect_size": effect_size,
            "n": n,
            "test_used": test_used,
            "is_normal": is_normal,
        }
    except Exception:
        return {"p_value": np.nan, "test_used": "error", "n": n}


def _add_significance_bracket(ax, x1, x2, y, p_value, line_height, text_offset):
    """
    Significance bracket with stars — consistent with stimulation_analysis.py.
    Draws horizontal bar + short vertical drops + star/ns annotation.
    """
    if p_value < 0.001:
        sig_text = "***"
    elif p_value < 0.01:
        sig_text = "**"
    elif p_value < 0.05:
        sig_text = "*"
    else:
        sig_text = "ns"

    ax.plot([x1, x2], [y, y], "k-", linewidth=1.5, clip_on=False)
    ax.plot([x1, x1], [y - line_height, y], "k-", linewidth=1.5, clip_on=False)
    ax.plot([x2, x2], [y - line_height, y], "k-", linewidth=1.5, clip_on=False)

    ax.text(
        (x1 + x2) / 2.0,
        y + text_offset,
        sig_text,
        ha="center",
        va="bottom",
        fontsize=FONT_SIZE_TICK - 1,
        fontweight="bold",
    )


def _plot_phase_locking_R_violin_panel(ax, phase_data):
    """
    Per-trial mean resultant length R of Δφ (REST vs RUN), paired across trials.
    R is the standard phase-locking / concentration metric for circular data (0 = uniform,
    1 = identical phases); here computed on left–right θ phase difference per trial.
    Includes paired statistical test with significance bracket.
    """
    r_rest, r_run = _paired_per_trial_R_rest_run(phase_data)
    positions = (0.7, 1.3)
    colors = (COLOR_PHASE_ROSE_REST, COLOR_PHASE_ROSE_RUN)
    labels = ("REST", "RUN")

    if r_rest.size == 0:
        ax.text(
            0.5,
            0.5,
            "Insufficient paired trials\n(need REST & RUN qualifying)",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=FONT_SIZE_LABEL,
            color="gray",
        )
        ax.set_title(
            "Per-trial R (mean resultant length)",
            fontsize=FONT_SIZE_TITLE,
            fontweight="normal",
            pad=8,
        )
        style_axis(ax)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
        ax.set_xlim(0.15, 1.85)
        ax.set_ylabel("R (mean resultant length)", fontsize=FONT_SIZE_LABEL)
        return

    _phase_locking_plot_half_violin_R(ax, r_rest, r_run, positions, colors)

    all_y = np.concatenate([r_rest, r_run])
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    y_range = y_max - y_min
    if y_range < 1e-10:
        y_range = max(abs(y_max), 0.05) * 0.3

    test_result = _perform_paired_test(r_rest, r_run)
    p_val = test_result.get("p_value", np.nan)
    test_name = test_result.get("test_used", "n/a")
    cohens_d = test_result.get("effect_size", np.nan)

    bracket_y = y_max + 0.06 * y_range
    line_height = 0.015 * y_range
    text_offset = 0.005 * y_range

    if np.isfinite(p_val):
        _add_significance_bracket(
            ax, positions[0], positions[1], bracket_y, p_val,
            line_height=line_height, text_offset=text_offset,
        )
        top_pad = bracket_y + 0.12 * y_range
    else:
        top_pad = y_max + 0.15 * y_range

    ax.set_ylim(y_min - y_range * 0.40, top_pad)

    mean_rest = float(np.mean(r_rest))
    mean_run = float(np.mean(r_run))
    ax.set_ylabel("R (mean resultant length)", fontsize=FONT_SIZE_LABEL)
    ax.set_title(
        "Per-trial R, paired REST vs RUN",
        fontsize=FONT_SIZE_TITLE - 1,
        fontweight="normal",
        pad=8,
    )
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
    ax.set_xlim(0.15, 1.85)
    style_axis(ax)

    p_str = (
        f"p = {p_val:.4f}" if np.isfinite(p_val) and p_val >= 0.001
        else f"p = {p_val:.2e}" if np.isfinite(p_val)
        else "p = n/a"
    )
    d_str = f"d = {cohens_d:.2f}" if np.isfinite(cohens_d) else "d = n/a"
    stats_txt = (
        f"Mean R  REST = {mean_rest:.3f}\n"
        f"Mean R  RUN  = {mean_run:.3f}\n"
        f"N paired = {len(r_rest)}\n"
        f"{test_name}, {p_str}\n"
        f"Cohen's {d_str}"
    )
    ax.text(
        0.02,
        0.98,
        stats_txt,
        transform=ax.transAxes,
        fontsize=FONT_SIZE_TICK - 1,
        va="top",
        ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.55", alpha=0.92),
    )


# =============================================================================
# AXIS STYLING
# =============================================================================

def style_axis(ax, show_spines=True):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if show_spines:
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis="both", which="major", labelsize=FONT_SIZE_TICK,
                   width=TICK_WIDTH, length=TICK_LENGTH)


def clean_axis(ax):
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])


def add_vscale_bar(ax, y_center, height, label, x_frac=-0.02, color="black",
                   fontsize=FONT_SIZE_SCALEBAR):
    """Vertical scale bar in axis-fraction coordinates (clip_off for publication margins)."""
    trans = ax.get_yaxis_transform()
    ax.plot([x_frac, x_frac], [y_center - height / 2, y_center + height / 2],
            transform=trans, color=color, lw=SCALEBAR_LINEWIDTH, clip_on=False, solid_capstyle="butt")
    # Right-side bars: place label slightly left of the line so it stays in the reserved margin
    if x_frac > 0.5:
        tx = x_frac - 0.018
        ha = "right"
    else:
        tx = x_frac - 0.01
        ha = "right"
    ax.text(
        tx, y_center, label, transform=trans, fontsize=fontsize,
        ha=ha, va="center", rotation=90, color=color, clip_on=False,
    )


# =============================================================================
# FIGURE 1: REPRESENTATIVE TRACES (2-column)
# =============================================================================

def fig_traces(data, animal, session, trial_num):
    """
    Two-column figure: Right HP | Left HP.
    Each column: LFP → GEVI → theta-filtered overlay → motion.
    Plots only TRACES_TIME_T0_S–TRACES_TIME_T1_S (bandpass on full trial, then slice).
    """
    t_full = np.asarray(data["time"], dtype=float)
    fs = data["fs"]

    win = (t_full >= TRACES_TIME_T0_S) & (t_full <= TRACES_TIME_T1_S)
    if not np.any(win):
        warnings.warn(
            f"No samples in [{TRACES_TIME_T0_S}, {TRACES_TIME_T1_S}] s; using full trace."
        )
        win = np.ones_like(t_full, dtype=bool)
    t = t_full[win]

    def _sl(a):
        if a is None:
            return None
        return np.asarray(a, dtype=float)[win]

    # Bandpass on full-length signals, then slice (reduces edge artifacts in window)
    f1_full = np.asarray(data["fiber1"], dtype=float)
    f2_full = np.asarray(data["fiber2"], dtype=float)
    gevi1_pct_full = f1_full * 100.0
    gevi2_pct_full = f2_full * 100.0
    gt1_full = bandpass(gevi1_pct_full, fs, *THETA_BAND)
    gt2_full = bandpass(gevi2_pct_full, fs, *THETA_BAND)

    lr_full = data["lfp_right"]
    ll_full = data["lfp_left"]
    ltr_full = bandpass(lr_full, fs, *THETA_BAND) if lr_full is not None else None
    ltl_full = bandpass(ll_full, fs, *THETA_BAND) if ll_full is not None else None

    speed_s = _sl(data["speed"])

    sides = [
        {
            "label": "Right Hippocampus",
            "fiber": _sl(f1_full),
            "lfp": _sl(lr_full),
            "gevi_theta": gt1_full[win],
            "lfp_theta": ltr_full[win] if ltr_full is not None else None,
            "color_fiber": COLOR_RIGHT_HP,
            "color_lfp": COLOR_LFP_RIGHT,
        },
        {
            "label": "Left Hippocampus",
            "fiber": _sl(f2_full),
            "lfp": _sl(ll_full),
            "gevi_theta": gt2_full[win],
            "lfp_theta": ltl_full[win] if ltl_full is not None else None,
            "color_fiber": COLOR_LEFT_HP,
            "color_lfp": COLOR_LFP_LEFT,
        },
    ]

    fig, axes = plt.subplots(
        4,
        2,
        figsize=FIGSIZE_TRACES_INCH,
        sharex=True,
        gridspec_kw={
            "hspace": 0.25,
            "wspace": TRACES_SUBPLOTS_WSPACE,
            "width_ratios": [1, 1],
        },
    )

    for col, side in enumerate(sides):
        fiber = side["fiber"]
        lfp = side["lfp"]
        c_fib = side["color_fiber"]
        c_lfp = side["color_lfp"]
        c_fib_t = tuple(np.clip(c_fib, 0, 1))
        c_lfp_t = tuple(np.clip(c_lfp, 0, 1))
        c_mot_t = tuple(np.clip(COLOR_MOTION, 0, 1))

        # --- Row 0: LFP ---
        ax = axes[0, col]
        if lfp is not None:
            ax.plot(t, lfp, color=c_lfp, lw=LW_RAW, alpha=0.9)
        clean_axis(ax)
        if lfp is not None:
            yc = float(np.nanmean(lfp))
            if col == 0:
                add_vscale_bar(ax, yc, SCALEBAR_LFP_UV, f"{SCALEBAR_LFP_UV:.0f} µV",
                               color=c_lfp_t)
            else:
                add_vscale_bar(ax, yc, SCALEBAR_LFP_UV, f"{SCALEBAR_LFP_UV:.0f} µV",
                               x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_lfp_t)
        ax.set_title(side["label"], fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)

        # --- Row 1: GEVI (ΔF/F) ---
        ax = axes[1, col]
        gevi_pct = fiber * 100.0  # ΔF/F to percent (fiber already ΔF/F in window)
        ax.plot(t, gevi_pct, color=c_fib, lw=LW_RAW, alpha=0.9)
        clean_axis(ax)
        yc_g = float(np.nanmean(gevi_pct))
        if col == 0:
            add_vscale_bar(ax, yc_g, SCALEBAR_GEVI_PERCENT,
                           f"{SCALEBAR_GEVI_PERCENT}%", color=c_fib_t)
        else:
            add_vscale_bar(ax, yc_g, SCALEBAR_GEVI_PERCENT,
                           f"{SCALEBAR_GEVI_PERCENT}%", x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_fib_t)

        # --- Row 2: Theta-filtered overlay (LFP + GEVI) ---
        ax = axes[2, col]
        gevi_theta = side["gevi_theta"]
        ax.plot(t, gevi_theta, color=c_fib, lw=LW_BAND, alpha=0.85, label="GEVI θ")
        if side["lfp_theta"] is not None:
            lfp_theta = side["lfp_theta"]
            ax2 = ax.twinx()
            ax2.plot(t, lfp_theta, color=c_lfp, lw=LW_BAND, alpha=0.85, label="LFP θ")
            clean_axis(ax2)
            if col == 0:
                add_vscale_bar(ax2, 0.0, SCALEBAR_LFP_BAND_UV,
                               f"{SCALEBAR_LFP_BAND_UV:.0f} µV",
                               x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_lfp_t)
            else:
                add_vscale_bar(ax2, 0.0, SCALEBAR_LFP_BAND_UV,
                               f"{SCALEBAR_LFP_BAND_UV:.0f} µV",
                               x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_lfp_t)
        clean_axis(ax)
        if col == 0:
            add_vscale_bar(ax, 0.0, SCALEBAR_GEVI_BAND_PERCENT,
                           f"{SCALEBAR_GEVI_BAND_PERCENT}%", color=c_fib_t)
        else:
            add_vscale_bar(ax, 0.0, SCALEBAR_GEVI_BAND_PERCENT,
                           f"{SCALEBAR_GEVI_BAND_PERCENT}%", x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_fib_t)

        # --- Row 3: Motion ---
        ax = axes[3, col]
        if speed_s is not None:
            ax.plot(t, speed_s, color=COLOR_MOTION, lw=LW_RAW, alpha=0.9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.set_yticks([])
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
        ax.tick_params(axis="x", labelsize=FONT_SIZE_TICK,
                       width=TICK_WIDTH, length=TICK_LENGTH)
        ax.set_xlabel("Time (s)", fontsize=FONT_SIZE_LABEL)
        if speed_s is not None:
            yc_m = float(np.nanmean(speed_s))
            if col == 0:
                add_vscale_bar(ax, yc_m, SCALEBAR_MOTION_CM_S,
                               f"{SCALEBAR_MOTION_CM_S:.0f} cm/s", color=c_mot_t)
            else:
                add_vscale_bar(ax, yc_m, SCALEBAR_MOTION_CM_S,
                               f"{SCALEBAR_MOTION_CM_S:.0f} cm/s",
                               x_frac=TRACES_SCALEBAR_X_RIGHT, color=c_mot_t)

    # Time scale bar (bottom, data coordinates)
    ax_bot = axes[3, 0]
    xlim = ax_bot.get_xlim()
    x0 = xlim[0] + 0.02 * (xlim[1] - xlim[0])
    y0 = ax_bot.get_ylim()[0]
    ax_bot.plot([x0, x0 + SCALEBAR_T_SECONDS], [y0, y0], color="black",
                lw=SCALEBAR_LINEWIDTH, clip_on=False, solid_capstyle="butt")
    ax_bot.text(x0 + SCALEBAR_T_SECONDS / 2, y0, f"{SCALEBAR_T_SECONDS:.0f} s",
                ha="center", va="top", fontsize=FONT_SIZE_SCALEBAR, color="black")

    fig.suptitle(
        f"{animal} - {session} Trial {trial_num} — Multi-site Fiber Traces "
        f"({TRACES_TIME_T0_S:.0f}–{TRACES_TIME_T1_S:.0f} s)",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98,
    )
    fig.subplots_adjust(
        left=TRACES_FIG_LEFT,
        right=TRACES_FIG_RIGHT,
        top=TRACES_FIG_TOP,
        bottom=TRACES_FIG_BOTTOM,
        hspace=0.25,
    )
    return fig


# =============================================================================
# FIGURE 2: TIME-FREQUENCY SPECTROGRAMS (2-column)
# =============================================================================

def fig_spectrograms(data, animal, session, trial_num):
    """
    Two-column: Right HP | Left HP.
    Each column: speed heatmap → LFP TFR → fiber TFR.

    Colour limits: vmin/vmax for LFP are computed from *both* LFP spectrograms together
    (2nd–98th percentile of log10 power); fiber limits likewise from both fibers.
    Left vs right in the same row therefore share one scale — apparent brightness
    differences reflect actual PSD differences on that scale, not per-panel rescaling.
    (LFP and fiber rows use separate colour bars because signals differ in physical units.)
    """
    t = data["time"]
    fs = data["fs"]

    sides = [
        {"label": "Right Hippocampus",
         "fiber": data["fiber1"], "lfp": data["lfp_right"]},
        {"label": "Left Hippocampus",
         "fiber": data["fiber2"], "lfp": data["lfp_left"]},
    ]

    fig = plt.figure(figsize=(18, 12))
    gs = GridSpec(3, 5, figure=fig,
                  height_ratios=[0.12, 1, 1],
                  width_ratios=[1, 0.025, 0.06, 1, 0.025],
                  hspace=0.25, wspace=0.06)

    # Compute common colour limits
    all_lfp_db, all_fib_db = [], []
    spec_cache = {}
    for col, side in enumerate(sides):
        for sig_type, sig_data in [("lfp", side["lfp"]), ("fiber", side["fiber"])]:
            if sig_data is None:
                continue
            freq, tb, Sxx = compute_spectrogram(sig_data, fs)
            tb_shifted = tb  # time already aligned to 0
            fmask = (freq >= FREQ_RANGE[0]) & (freq <= FREQ_RANGE[1])
            Sxx_db = 10 * np.log10(Sxx[fmask, :] + 1e-20)
            spec_cache[(col, sig_type)] = (freq[fmask], tb_shifted, Sxx_db)
            if sig_type == "lfp":
                all_lfp_db.append(Sxx_db)
            else:
                all_fib_db.append(Sxx_db)

    def _clim(arrs):
        if not arrs:
            return None, None
        combined = np.concatenate([a.ravel() for a in arrs])
        return np.nanpercentile(combined, 2), np.nanpercentile(combined, 98)

    vmin_lfp, vmax_lfp = _clim(all_lfp_db)
    vmin_fib, vmax_fib = _clim(all_fib_db)

    col_map = {0: (0, 1), 1: (3, 4)}  # (data_col, cbar_col)

    for col, side in enumerate(sides):
        dc, cc = col_map[col]

        # --- Row 0: Speed heatmap ---
        ax_sp = fig.add_subplot(gs[0, dc])
        if data["speed"] is not None:
            speed_2d = data["speed"].reshape(1, -1)
            vmax_s = np.nanpercentile(data["speed"], 99)
            ax_sp.imshow(speed_2d, aspect="auto", cmap=CMAP_SPEED,
                         extent=[t[0], t[-1], 0, 1], origin="lower",
                         interpolation="bilinear",
                         vmin=0, vmax=max(vmax_s, 1.0))
        ax_sp.set_ylim([0, 1])
        ax_sp.set_yticks([])
        ax_sp.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK)
        if col == 0:
            ax_sp.set_ylabel("Speed\n(cm/s)", fontsize=FONT_SIZE_LABEL,
                             rotation=0, ha="right", va="center")
        for sp in ("top", "right"):
            ax_sp.spines[sp].set_visible(False)
        ax_sp.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax_sp.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
        ax_sp.set_title(side["label"], fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)

        # Speed colorbar
        if data["speed"] is not None and len(ax_sp.images) > 0:
            cax_sp = fig.add_subplot(gs[0, cc])
            cb_sp = fig.colorbar(ax_sp.images[0], cax=cax_sp, orientation="vertical")
            cb_sp.set_label("cm/s", fontsize=FONT_SIZE_COLORBAR - 2, fontweight="bold")
            cb_sp.ax.tick_params(
                labelsize=FONT_SIZE_TICK - 2,
                width=TICK_WIDTH,
                length=max(TICK_LENGTH - 2, 4),
            )

        # --- Row 1: LFP TFR ---
        ax_lfp = fig.add_subplot(gs[1, dc])
        key_lfp = (col, "lfp")
        if key_lfp in spec_cache:
            fq, tb, Sxx_db = spec_cache[key_lfp]
            im_lfp = ax_lfp.pcolormesh(tb, fq, Sxx_db,
                                        shading="gouraud", cmap="viridis",
                                        vmin=vmin_lfp, vmax=vmax_lfp)
            cax_lfp = fig.add_subplot(gs[1, cc])
            cb_lfp = fig.colorbar(im_lfp, cax=cax_lfp, orientation="vertical")
            cb_lfp.set_label("Power (dB)", fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
            cb_lfp.ax.tick_params(
                labelsize=FONT_SIZE_TICK - 2,
                width=TICK_WIDTH,
                length=max(TICK_LENGTH - 2, 4),
            )
        ax_lfp.set_ylim(FREQ_RANGE)
        ax_lfp.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK)
        if col == 0:
            ax_lfp.set_ylabel("LFP\nFrequency (Hz)", fontsize=FONT_SIZE_LABEL)
        for sp in ("top", "right"):
            ax_lfp.spines[sp].set_visible(False)
        ax_lfp.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax_lfp.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)

        # --- Row 2: Fiber TFR ---
        ax_fib = fig.add_subplot(gs[2, dc], sharex=ax_lfp)
        key_fib = (col, "fiber")
        if key_fib in spec_cache:
            fq, tb, Sxx_db = spec_cache[key_fib]
            im_fib = ax_fib.pcolormesh(tb, fq, Sxx_db,
                                        shading="gouraud", cmap="viridis",
                                        vmin=vmin_fib, vmax=vmax_fib)
            cax_fib = fig.add_subplot(gs[2, cc])
            cb_fib = fig.colorbar(im_fib, cax=cax_fib, orientation="vertical")
            cb_fib.set_label("Power (dB)", fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
            cb_fib.ax.tick_params(
                labelsize=FONT_SIZE_TICK - 2,
                width=TICK_WIDTH,
                length=max(TICK_LENGTH - 2, 4),
            )
        ax_fib.set_ylim(FREQ_RANGE)
        ax_fib.set_xlabel("Time (s)", fontsize=FONT_SIZE_LABEL)
        ax_fib.tick_params(labelsize=FONT_SIZE_TICK)
        if col == 0:
            ax_fib.set_ylabel("Fiber\nFrequency (Hz)", fontsize=FONT_SIZE_LABEL)
        for sp in ("top", "right"):
            ax_fib.spines[sp].set_visible(False)
        ax_fib.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax_fib.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)

    # Hide spacer column
    for row in range(3):
        ax_gap = fig.add_subplot(gs[row, 2])
        ax_gap.set_visible(False)

    fig.suptitle(f"{animal} - {session} Trial {trial_num} — Time-Frequency Spectrograms",
                 fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98)
    fig.subplots_adjust(left=0.09, right=0.95, top=0.92, bottom=0.07)
    return fig


# =============================================================================
# FIGURE 3: FIBER-FIBER THETA CROSS-CORRELATION (REST VS RUN)
# =============================================================================

def compute_xcorr_fiber_by_behavior(
    fiber1,
    fiber2,
    fs,
    time_vec,
    speed_cm_s,
    band=THETA_BAND,
    window_sec=XCORR_WINDOW_SEC,
    step_sec=XCORR_STEP_SEC,
    max_lag_ms=XCORR_MAX_LAG_MS,
    trim_sec=XCORR_EDGE_TRIM_SEC,
):
    """
    Sliding-window cross-correlation; windows are split into REST vs RUN
    (classify_behavior.m / classify_rest_run_masks).

    Per window: segments are de-meaned, unit-variance scaled, Hann-tapered, then
    correlate(..., mode='full') with division by sqrt(Σs1²·Σs2²) so each curve is a
    normalised correlation in ~[-1, 1]. For the *mean* across windows, see
    mean_xcorr_over_windows() (Fisher-z average by default).
    """
    f1_filt = bandpass(fiber1, fs, band[0], band[1])
    f2_filt = bandpass(fiber2, fs, band[0], band[1])

    trim_samp = int(trim_sec * fs)
    if len(f1_filt) <= 2 * trim_samp + 10:
        raise ValueError("Trace too short for cross-correlation after edge trim.")
    f1_filt = f1_filt[trim_samp:-trim_samp]
    f2_filt = f2_filt[trim_samp:-trim_samp]
    n = len(f1_filt)

    win_samp = int(window_sec * fs)
    step_samp = max(1, int(step_sec * fs))
    max_lag_samp = int(max_lag_ms / 1000 * fs)

    starts = np.arange(0, n - win_samp + 1, step_samp)
    n_lags = 2 * max_lag_samp + 1
    lags_samp = np.arange(-max_lag_samp, max_lag_samp + 1)
    lags_ms = lags_samp / fs * 1000

    if speed_cm_s is None:
        warnings.warn(
            "Speed missing: cannot classify REST vs RUN; all windows assigned to REST.",
            UserWarning,
            stacklevel=2,
        )
        is_rest = np.ones(len(time_vec), dtype=bool)
        is_run = np.zeros_like(is_rest)
        is_excluded = np.zeros_like(is_rest)
        beh_txt = "no speed — all windows labelled REST"
    else:
        is_rest, is_run, is_excluded = classify_rest_run_masks(speed_cm_s, fs)
        n_map = min(
            len(np.asarray(time_vec, dtype=float).ravel()),
            len(is_rest),
            len(is_run),
            len(is_excluded),
        )
        is_rest = is_rest[:n_map]
        is_run = is_run[:n_map]
        is_excluded = is_excluded[:n_map]
        n_ex = int(np.sum(is_excluded))
        n_tot = len(is_rest)
        beh_txt = (
            f"mode={BEHAVIOR_CLASSIFICATION_MODE}, "
            f"REST {100 * np.mean(is_rest):.1f}%, RUN {100 * np.mean(is_run):.1f}%"
        )
        if BEHAVIOR_CLASSIFICATION_MODE.lower() == "clear" and n_tot:
            beh_txt += f", excluded {100 * n_ex / n_tot:.1f}%"

    cols_rest, cols_run, centers_rest, centers_run = [], [], [], []

    for s0 in starts:
        center_sec = trim_sec + (s0 + win_samp / 2) / fs
        ti = _time_index_at_sec(time_vec, center_sec)
        if ti < len(is_excluded) and bool(is_excluded[ti]):
            continue
        if ti < len(is_run) and bool(is_run[ti]):
            bucket = "run"
        elif ti < len(is_rest) and bool(is_rest[ti]):
            bucket = "rest"
        else:
            continue

        seg1 = f1_filt[s0 : s0 + win_samp].copy()
        seg2 = f2_filt[s0 : s0 + win_samp].copy()

        seg1 -= np.mean(seg1)
        seg2 -= np.mean(seg2)
        std1, std2 = np.std(seg1), np.std(seg2)
        if std1 < 1e-10 or std2 < 1e-10:
            continue

        seg1 /= std1
        seg2 /= std2

        hann = np.hanning(win_samp)
        seg1 *= hann
        seg2 *= hann

        full_xcorr = signal.correlate(seg1, seg2, mode="full")
        norm = np.sqrt(np.sum(seg2 ** 2) * np.sum(seg1 ** 2))
        if norm > 0:
            full_xcorr /= norm

        mid = len(full_xcorr) // 2
        col = full_xcorr[mid - max_lag_samp : mid + max_lag_samp + 1]

        if bucket == "run":
            cols_run.append(col)
            centers_run.append(center_sec)
        else:
            cols_rest.append(col)
            centers_rest.append(center_sec)

    mat_rest = np.column_stack(cols_rest) if cols_rest else None
    mat_run = np.column_stack(cols_run) if cols_run else None
    centers_rest = np.array(centers_rest) if centers_rest else np.array([])
    centers_run = np.array(centers_run) if centers_run else np.array([])

    return {
        "lags_ms": lags_ms,
        "mat_rest": mat_rest,
        "mat_run": mat_run,
        "centers_rest": centers_rest,
        "centers_run": centers_run,
        "behavior_note": beh_txt,
        "n_windows_rest": len(cols_rest),
        "n_windows_run": len(cols_run),
    }


def _peak_lag_search_half_width_ms(band):
    """
    Half-width (ms) for peak |r| search on mean correlogram.
    XCORR_PEAK_LAG_LIMIT_MS == 'half_theta' -> 1000/(2*f_mid) at centre of `band` (½ cycle).
    Otherwise treated as a float ±ms cap (e.g. 50.0).
    """
    lim = XCORR_PEAK_LAG_LIMIT_MS
    if isinstance(lim, str) and lim.lower() in ("half_theta", "half_theta_cycle"):
        f_mid = 0.5 * (float(band[0]) + float(band[1]))
        return 1000.0 / (2.0 * max(f_mid, 1e-6))
    return float(lim)


def mean_xcorr_over_windows(mat, use_fisher_z=None):
    """
    Mean cross-correlogram across sliding windows (axis=1).

    Each window is already a normalised correlation curve (zero-mean unit-variance
    segments + Hann taper + normalised correlate). To average coefficients without
    letting windows with larger |r| dominate, optionally use Fisher z: mean(atanh(r))
    then tanh (see XCORR_MEAN_USE_FISHER_Z).
    """
    if mat is None or mat.size == 0:
        return None
    if use_fisher_z is None:
        use_fisher_z = XCORR_MEAN_USE_FISHER_Z
    mat = np.asarray(mat, dtype=float)
    if use_fisher_z:
        m = np.clip(mat, -0.999999, 0.999999)
        with np.errstate(invalid="ignore"):
            z = np.arctanh(m)
        mean_z = np.nanmean(z, axis=1)
        out = np.tanh(mean_z)
        return out
    return np.nanmean(mat, axis=1)


def peak_lag_index_restricted(mean_r, lags_ms, band):
    """Index of max |mean_r| among lags with |lag| <= peak half-width; returns (idx, limit_ms)."""
    mean_r = np.asarray(mean_r, dtype=float)
    lags_ms = np.asarray(lags_ms, dtype=float)
    limit_ms = _peak_lag_search_half_width_ms(band)
    abs_r = np.abs(mean_r)
    mask = np.abs(lags_ms) <= limit_ms
    sub_idx = np.where(mask)[0]
    if sub_idx.size > 0 and np.any(np.isfinite(abs_r[sub_idx])):
        j = int(np.nanargmax(abs_r[sub_idx]))
        ir = int(sub_idx[j])
    else:
        ir = int(np.nanargmax(abs_r))
    return ir, limit_ms


def collect_trialwise_xcorr(animal, session=None, band=THETA_BAND):
    """
    Per-trial mean cross-correlograms (REST/RUN) across all sessions for an animal.

    Returns dict with lags_ms, per-trial dicts, stacked matrices for heatmaps,
    grand mean / SEM for line plots, and per-trial peak |r| for violin.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None

    if session is not None:
        if session not in sessions_map:
            return None
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    trials_out = []
    lags_ms_ref = None
    failed = 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                failed += 1
                continue

            f1 = np.asarray(data["fiber1"], dtype=float).ravel()
            f2 = np.asarray(data["fiber2"], dtype=float).ravel()
            if f2.size == 0 or not np.any(np.isfinite(f2)):
                failed += 1
                continue

            try:
                out = compute_xcorr_fiber_by_behavior(
                    f1, f2, float(data["fs"]), data["time"], data["speed"],
                    band=band,
                )
            except Exception:
                failed += 1
                continue

            lags_ms = out["lags_ms"]
            if lags_ms_ref is None:
                lags_ms_ref = lags_ms

            mean_rest = mean_xcorr_over_windows(out["mat_rest"])
            mean_run = mean_xcorr_over_windows(out["mat_run"])

            peak_r_rest, peak_lag_rest = np.nan, np.nan
            if mean_rest is not None:
                ir, _ = peak_lag_index_restricted(mean_rest, lags_ms, band)
                peak_r_rest = float(np.abs(mean_rest[ir]))
                peak_lag_rest = float(lags_ms[ir])

            peak_r_run, peak_lag_run = np.nan, np.nan
            if mean_run is not None:
                irn, _ = peak_lag_index_restricted(mean_run, lags_ms, band)
                peak_r_run = float(np.abs(mean_run[irn]))
                peak_lag_run = float(lags_ms[irn])

            trials_out.append({
                "session": sess_key,
                "trial_num": trial_num,
                "mean_rest": mean_rest,
                "mean_run": mean_run,
                "peak_r_rest": peak_r_rest,
                "peak_r_run": peak_r_run,
                "peak_lag_rest": peak_lag_rest,
                "peak_lag_run": peak_lag_run,
                "n_windows_rest": out["n_windows_rest"],
                "n_windows_run": out["n_windows_run"],
            })

    if not trials_out or lags_ms_ref is None:
        return None

    rest_rows = [t["mean_rest"] for t in trials_out if t["mean_rest"] is not None]
    run_rows = [t["mean_run"] for t in trials_out if t["mean_run"] is not None]
    mat_all_rest = np.vstack(rest_rows) if rest_rows else None
    mat_all_run = np.vstack(run_rows) if run_rows else None

    def _grand_mean_sem(mat, n):
        if mat is None or n < 1:
            return None, None
        if XCORR_MEAN_USE_FISHER_Z:
            m = np.clip(mat, -0.999999, 0.999999)
            z = np.arctanh(m)
            grand = np.tanh(np.nanmean(z, axis=0))
        else:
            grand = np.nanmean(mat, axis=0)
        sem = np.nanstd(mat, axis=0, ddof=1) / np.sqrt(n) if n > 1 else None
        return grand, sem

    grand_mean_rest, sem_rest = _grand_mean_sem(mat_all_rest, len(rest_rows))
    grand_mean_run, sem_run = _grand_mean_sem(mat_all_run, len(run_rows))

    print(
        f"  Cross-correlation all trials: {len(trials_out)} loaded, {failed} failed, "
        f"{len(rest_rows)} with REST, {len(run_rows)} with RUN"
    )

    return {
        "lags_ms": lags_ms_ref,
        "trials": trials_out,
        "mat_all_rest": mat_all_rest,
        "mat_all_run": mat_all_run,
        "grand_mean_rest": grand_mean_rest,
        "grand_mean_run": grand_mean_run,
        "sem_rest": sem_rest,
        "sem_run": sem_run,
        "n_loaded": len(trials_out),
        "n_failed": failed,
        "n_rest_trials": len(rest_rows),
        "n_run_trials": len(run_rows),
    }


# =============================================================================
# UNSPLIT CROSS-CORRELATION (all epochs) + CIRCULAR-SHIFT SURROGATE TEST
# =============================================================================

def compute_xcorr_fiber_unsplit(
    fiber1, fiber2, fs, band=THETA_BAND,
    window_sec=XCORR_WINDOW_SEC, step_sec=XCORR_STEP_SEC,
    max_lag_ms=XCORR_MAX_LAG_MS, trim_sec=XCORR_EDGE_TRIM_SEC,
):
    """
    Sliding-window cross-correlation over all epochs (no REST/RUN split).
    Returns dict with lags_ms, mat (n_lags x n_windows), centers_sec.
    """
    f1_filt = bandpass(fiber1, fs, band[0], band[1])
    f2_filt = bandpass(fiber2, fs, band[0], band[1])

    trim_samp = int(trim_sec * fs)
    if len(f1_filt) <= 2 * trim_samp + 10:
        raise ValueError("Trace too short for cross-correlation after edge trim.")
    f1_filt = f1_filt[trim_samp:-trim_samp]
    f2_filt = f2_filt[trim_samp:-trim_samp]
    n = len(f1_filt)

    win_samp = int(window_sec * fs)
    step_samp = max(1, int(step_sec * fs))
    max_lag_samp = int(max_lag_ms / 1000 * fs)

    starts = np.arange(0, n - win_samp + 1, step_samp)
    lags_samp = np.arange(-max_lag_samp, max_lag_samp + 1)
    lags_ms_arr = lags_samp / fs * 1000

    cols, centers = [], []
    for s0 in starts:
        seg1 = f1_filt[s0: s0 + win_samp].copy()
        seg2 = f2_filt[s0: s0 + win_samp].copy()
        seg1 -= np.mean(seg1)
        seg2 -= np.mean(seg2)
        std1, std2 = np.std(seg1), np.std(seg2)
        if std1 < 1e-10 or std2 < 1e-10:
            continue
        seg1 /= std1
        seg2 /= std2
        hann = np.hanning(win_samp)
        seg1 *= hann
        seg2 *= hann
        full_xcorr = signal.correlate(seg1, seg2, mode="full")
        norm = np.sqrt(np.sum(seg2 ** 2) * np.sum(seg1 ** 2))
        if norm > 0:
            full_xcorr /= norm
        mid = len(full_xcorr) // 2
        col = full_xcorr[mid - max_lag_samp: mid + max_lag_samp + 1]
        cols.append(col)
        centers.append(trim_sec + (s0 + win_samp / 2) / fs)

    mat = np.column_stack(cols) if cols else None
    return {
        "lags_ms": lags_ms_arr,
        "mat": mat,
        "centers": np.array(centers) if centers else np.array([]),
        "n_windows": len(cols),
    }


def _circular_shift_surrogate_peak_r(
    fiber1, fiber2, fs, band=THETA_BAND, n_surrogates=XCORR_UNSPLIT_N_SURROGATES,
):
    """
    Null distribution for peak |r| via circular time-shift surrogates.

    For each surrogate: circularly shift one filtered signal by a random amount
    (>= 2 theta cycles), recompute sliding-window xcorr, extract peak |r|.
    """
    rng = np.random.default_rng(42)
    f1_filt = bandpass(fiber1, fs, band[0], band[1])
    f2_filt = bandpass(fiber2, fs, band[0], band[1])

    trim_samp = int(XCORR_EDGE_TRIM_SEC * fs)
    if len(f1_filt) <= 2 * trim_samp + 10:
        return np.array([])
    f1_filt = f1_filt[trim_samp:-trim_samp]
    f2_filt = f2_filt[trim_samp:-trim_samp]
    n = len(f1_filt)

    f_mid = 0.5 * (band[0] + band[1])
    min_shift = int(XCORR_UNSPLIT_MIN_SHIFT_CYCLES / f_mid * fs)
    max_shift = n - min_shift
    if max_shift <= min_shift:
        return np.array([])

    win_samp = int(XCORR_WINDOW_SEC * fs)
    step_samp = max(1, int(XCORR_STEP_SEC * fs))
    max_lag_samp = int(XCORR_MAX_LAG_MS / 1000 * fs)
    starts = np.arange(0, n - win_samp + 1, step_samp)
    limit_ms = _peak_lag_search_half_width_ms(band)
    lags_samp = np.arange(-max_lag_samp, max_lag_samp + 1)
    lags_ms_arr = lags_samp / fs * 1000

    null_peaks = np.full(n_surrogates, np.nan)
    for si in range(n_surrogates):
        shift = rng.integers(min_shift, max_shift)
        f2_shifted = np.roll(f2_filt, shift)

        r_cols = []
        for s0 in starts:
            seg1 = f1_filt[s0: s0 + win_samp].copy()
            seg2 = f2_shifted[s0: s0 + win_samp].copy()
            seg1 -= np.mean(seg1)
            seg2 -= np.mean(seg2)
            std1, std2 = np.std(seg1), np.std(seg2)
            if std1 < 1e-10 or std2 < 1e-10:
                continue
            seg1 /= std1
            seg2 /= std2
            hann = np.hanning(win_samp)
            seg1 *= hann
            seg2 *= hann
            full_xcorr = signal.correlate(seg1, seg2, mode="full")
            norm_val = np.sqrt(np.sum(seg2 ** 2) * np.sum(seg1 ** 2))
            if norm_val > 0:
                full_xcorr /= norm_val
            mid = len(full_xcorr) // 2
            col = full_xcorr[mid - max_lag_samp: mid + max_lag_samp + 1]
            r_cols.append(col)

        if not r_cols:
            continue
        mean_r = mean_xcorr_over_windows(np.column_stack(r_cols))
        if mean_r is None:
            continue
        mask = np.abs(lags_ms_arr) <= limit_ms
        sub = np.where(mask)[0]
        if sub.size > 0:
            null_peaks[si] = float(np.nanmax(np.abs(mean_r[sub])))
        else:
            null_peaks[si] = float(np.nanmax(np.abs(mean_r)))

    return null_peaks[np.isfinite(null_peaks)]


def collect_trialwise_xcorr_unsplit(animal, session=None, band=THETA_BAND):
    """
    Per-trial unsplit (all-epoch) cross-correlograms + circular-shift surrogates.
    Returns dict with lags, per-trial mean correlograms, stacked matrix for heatmap,
    grand mean ± SEM, observed peak |r| per trial, and pooled null distribution.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    if session is not None:
        if session not in sessions_map:
            return None
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    trials_out = []
    lags_ms_ref = None
    failed = 0
    all_null_peaks = []

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                failed += 1
                continue

            f1 = np.asarray(data["fiber1"], dtype=float).ravel()
            f2 = np.asarray(data["fiber2"], dtype=float).ravel()
            if f2.size == 0 or not np.any(np.isfinite(f2)):
                failed += 1
                continue

            try:
                out = compute_xcorr_fiber_unsplit(f1, f2, float(data["fs"]), band=band)
            except Exception:
                failed += 1
                continue

            lags_ms = out["lags_ms"]
            if lags_ms_ref is None:
                lags_ms_ref = lags_ms

            mean_all = mean_xcorr_over_windows(out["mat"])
            if mean_all is None:
                failed += 1
                continue

            ir, _ = peak_lag_index_restricted(mean_all, lags_ms, band)
            peak_r = float(np.abs(mean_all[ir]))
            peak_lag = float(lags_ms[ir])

            print(f"    {sess_key} T{trial_num}: surrogates ({XCORR_UNSPLIT_N_SURROGATES})...",
                  end="", flush=True)
            null_peaks = _circular_shift_surrogate_peak_r(
                f1, f2, float(data["fs"]), band=band,
            )
            p_surr = np.nan
            if len(null_peaks) > 0:
                p_surr = float(np.mean(null_peaks >= peak_r))
                all_null_peaks.append(null_peaks)
            print(f" peak|r|={peak_r:.3f}, p_surr={p_surr:.4f}")

            trials_out.append({
                "session": sess_key,
                "trial_num": trial_num,
                "mean_all": mean_all,
                "peak_r": peak_r,
                "peak_lag": peak_lag,
                "n_windows": out["n_windows"],
                "p_surrogate": p_surr,
            })

    if not trials_out or lags_ms_ref is None:
        return None

    all_rows = [t["mean_all"] for t in trials_out]
    mat_all = np.vstack(all_rows) if all_rows else None

    def _grand_mean_sem(mat, n):
        if mat is None or n < 1:
            return None, None
        if XCORR_MEAN_USE_FISHER_Z:
            m = np.clip(mat, -0.999999, 0.999999)
            z = np.arctanh(m)
            grand = np.tanh(np.nanmean(z, axis=0))
        else:
            grand = np.nanmean(mat, axis=0)
        sem = np.nanstd(mat, axis=0, ddof=1) / np.sqrt(n) if n > 1 else None
        return grand, sem

    grand_mean, sem = _grand_mean_sem(mat_all, len(all_rows))

    observed_peaks = np.array([t["peak_r"] for t in trials_out])
    pooled_null = np.concatenate(all_null_peaks) if all_null_peaks else np.array([])

    pooled_p = np.nan
    if len(pooled_null) > 0 and len(observed_peaks) > 0:
        median_obs = float(np.median(observed_peaks))
        pooled_p = float(np.mean(pooled_null >= median_obs))

    n_sig = sum(1 for t in trials_out if t["p_surrogate"] < 0.05)

    print(
        f"  Unsplit xcorr all trials: {len(trials_out)} loaded, {failed} failed, "
        f"{n_sig}/{len(trials_out)} trials significant (p<0.05 surrogate)"
    )

    return {
        "lags_ms": lags_ms_ref,
        "trials": trials_out,
        "mat_all": mat_all,
        "grand_mean": grand_mean,
        "sem": sem,
        "n_loaded": len(trials_out),
        "n_failed": failed,
        "observed_peaks": observed_peaks,
        "pooled_null": pooled_null,
        "pooled_p": pooled_p,
        "n_sig_trials": n_sig,
    }


# ---------------------------------------------------------------------------
# Figure 6 helpers: fiber–LFP cross-correlation (4 combinations)
# ---------------------------------------------------------------------------

def collect_trialwise_fiber_lfp_xcorr(animal, session=None, band=THETA_BAND):
    """
    Per-trial unsplit theta cross-correlation for 4 fiber–LFP combinations
    (R–R ipsi, R–L contra, L–L ipsi, L–R contra).
    Pools across all sessions when *session* is None.

    Returns dict:
      combo_data  – per-combo aggregates (grand_mean, sem, peak_r_values, mat_all, ...)
      trial_peaks – list of dicts {session, trial_num, combo_short: peak_r, ...}
      lags_ms, n_loaded, n_failed
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    if session is not None:
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    combo_keys = [c["short"] for c in FIG6_COMBINATIONS]
    combo_data = {k: {"mean_rows": [], "peak_r_list": [], "peak_lag_list": []}
                  for k in combo_keys}
    trial_peaks = []
    lags_ms_ref = None
    n_loaded, n_failed = 0, 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue

            n_loaded += 1
            fs = float(data["fs"])
            tp = {"session": sess_key, "trial_num": trial_num}

            for combo in FIG6_COMBINATIONS:
                sig1 = np.asarray(data[combo["sig1_key"]], dtype=float).ravel()
                sig2 = data.get(combo["sig2_key"])
                if sig2 is None:
                    tp[combo["short"]] = np.nan
                    continue
                sig2 = np.asarray(sig2, dtype=float).ravel()
                if not np.any(np.isfinite(sig1)) or not np.any(np.isfinite(sig2)):
                    tp[combo["short"]] = np.nan
                    continue

                try:
                    out = compute_xcorr_fiber_unsplit(sig1, sig2, fs, band=band)
                except Exception:
                    tp[combo["short"]] = np.nan
                    continue

                if out["mat"] is None:
                    tp[combo["short"]] = np.nan
                    continue

                if lags_ms_ref is None:
                    lags_ms_ref = out["lags_ms"]

                mean_corr = mean_xcorr_over_windows(out["mat"])
                if mean_corr is None:
                    tp[combo["short"]] = np.nan
                    continue

                ir, _ = peak_lag_index_restricted(mean_corr, out["lags_ms"], band)
                peak_r = float(np.abs(mean_corr[ir]))
                peak_lag = float(out["lags_ms"][ir])

                cd = combo_data[combo["short"]]
                cd["mean_rows"].append(mean_corr)
                cd["peak_r_list"].append(peak_r)
                cd["peak_lag_list"].append(peak_lag)
                tp[combo["short"]] = peak_r

            trial_peaks.append(tp)
            combo_strs = "  ".join(
                f"{c['short']}={tp.get(c['short'], np.nan):.3f}"
                for c in FIG6_COMBINATIONS
            )
            print(f"    {sess_key} T{trial_num}: {combo_strs}")

    if lags_ms_ref is None:
        return None

    for key, cd in combo_data.items():
        rows = cd["mean_rows"]
        if rows:
            mat = np.vstack(rows)
            cd["mat_all"] = mat
            n = len(rows)
            if XCORR_MEAN_USE_FISHER_Z:
                m = np.clip(mat, -0.999999, 0.999999)
                z = np.arctanh(m)
                cd["grand_mean"] = np.tanh(np.nanmean(z, axis=0))
            else:
                cd["grand_mean"] = np.nanmean(mat, axis=0)
            cd["sem"] = (np.nanstd(mat, axis=0, ddof=1) / np.sqrt(n)
                         if n > 1 else None)
            cd["n_trials"] = n
            cd["peak_r_values"] = np.array(cd["peak_r_list"])
        else:
            cd["mat_all"] = None
            cd["grand_mean"] = None
            cd["sem"] = None
            cd["n_trials"] = 0
            cd["peak_r_values"] = np.array([])

    print(
        f"  Fiber-LFP xcorr: {n_loaded} trials loaded, {n_failed} failed"
    )
    return {
        "combo_data": combo_data,
        "trial_peaks": trial_peaks,
        "lags_ms": lags_ms_ref,
        "n_loaded": n_loaded,
        "n_failed": n_failed,
    }


def _plot_xcorr_heatmap_panel(ax, mat, centers, lags_ms, vmin, vmax, title):
    """Draw one REST or RUN heatmap; return mappable or None."""
    if mat is None or mat.size == 0 or len(centers) == 0:
        ax.text(
            0.5,
            0.5,
            "Insufficient data",
            transform=ax.transAxes,
            ha="center",
            va="center",
            fontsize=FONT_SIZE_LABEL,
            color="gray",
        )
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
        style_axis(ax)
        return None
    extent = [float(centers[0]), float(centers[-1]), float(lags_ms[0]), float(lags_ms[-1])]
    im = ax.imshow(
        mat,
        aspect="auto",
        origin="lower",
        interpolation="bilinear",
        cmap="inferno",
        rasterized=True,
        extent=extent,
        vmin=vmin,
        vmax=vmax,
    )
    ax.axhline(0, color="white", ls="--", lw=1.2, alpha=0.7)
    ax.set_ylabel("Lag (ms)", fontsize=FONT_SIZE_LABEL)
    ax.set_xlabel("Time (s)", fontsize=FONT_SIZE_LABEL)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
    style_axis(ax)
    return im


def _normalize_xcorr_rows(mat):
    """
    Per-trial normalization of correlogram rows (XCORR_ALL_TRIALS_NORMALIZE).
    'zscore': each row → (row - mean) / std  — equalises magnitude, preserves shape.
    'peak':   each row → row / max(|row|)    — scales to [-1, 1].
    None:     no change.
    """
    mode = XCORR_ALL_TRIALS_NORMALIZE
    if mode is None or mat is None or mat.size == 0:
        return mat
    out = mat.astype(float).copy()
    for i in range(out.shape[0]):
        row = out[i]
        if mode == "zscore":
            mu = np.nanmean(row)
            sd = np.nanstd(row)
            if sd > 1e-15:
                out[i] = (row - mu) / sd
        elif mode == "peak":
            mx = np.nanmax(np.abs(row))
            if mx > 1e-15:
                out[i] = row / mx
    return out


def _plot_xcorr_all_trials_heatmap(ax, mat, lags_ms, n_trials, vmin, vmax, cmap, title):
    """
    Heatmap where each row is one trial's mean correlogram.
    mat shape: (n_trials, n_lags). X = lag (ms), Y = trial index.
    Applies per-trial normalization if XCORR_ALL_TRIALS_NORMALIZE is set.
    """
    if mat is None or mat.size == 0:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
        style_axis(ax)
        return None

    mat_plot = _normalize_xcorr_rows(mat)

    if XCORR_ALL_TRIALS_NORMALIZE is not None:
        abs_max = float(np.nanpercentile(np.abs(mat_plot), 98))
        abs_max = max(abs_max, 0.5)
        vmin, vmax = -abs_max, abs_max

    extent = [float(lags_ms[0]), float(lags_ms[-1]), 0.5, n_trials + 0.5]
    im = ax.imshow(
        mat_plot, aspect="auto", origin="lower", interpolation="nearest",
        cmap=cmap, extent=extent, vmin=vmin, vmax=vmax,
    )
    ax.axvline(0, color="white", ls="--", lw=1.0, alpha=0.7)
    ax.set_ylabel("Trial", fontsize=FONT_SIZE_LABEL)
    ax.set_xlabel("Lag (ms)", fontsize=FONT_SIZE_LABEL)
    norm_lbl = f" [{XCORR_ALL_TRIALS_NORMALIZE}]" if XCORR_ALL_TRIALS_NORMALIZE else ""
    ax.set_title(title + norm_lbl, fontsize=FONT_SIZE_TITLE, fontweight="bold", pad=8)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    style_axis(ax)
    return im


def _paired_per_trial_peak_r_rest_run(all_trials_data):
    """Trials with both REST and RUN peak |r|: paired arrays."""
    r_rest, r_run = [], []
    if all_trials_data is None:
        return np.array([]), np.array([])
    for t in all_trials_data.get("trials", []):
        rr = t.get("peak_r_rest", np.nan)
        rn = t.get("peak_r_run", np.nan)
        if np.isfinite(rr) and np.isfinite(rn):
            r_rest.append(rr)
            r_run.append(rn)
    return np.asarray(r_rest, dtype=float), np.asarray(r_run, dtype=float)


def _plot_xcorr_peak_r_violin_panel(ax, all_trials_data):
    """
    Half-violin + box + scatter + paired lines + significance bracket for
    per-trial peak |r| (REST vs RUN). Reuses violin geometry from phase figure.
    """
    r_rest, r_run = _paired_per_trial_peak_r_rest_run(all_trials_data)
    positions = (0.7, 1.3)
    colors = (COLOR_PHASE_ROSE_REST, COLOR_PHASE_ROSE_RUN)
    labels = ("REST", "RUN")

    if r_rest.size == 0:
        ax.text(0.5, 0.5, "Insufficient paired trials\n(need REST & RUN)",
                transform=ax.transAxes, ha="center", va="center",
                fontsize=FONT_SIZE_LABEL, color="gray")
        ax.set_title("Peak |r| per trial", fontsize=FONT_SIZE_TITLE,
                      fontweight="normal", pad=8)
        style_axis(ax)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
        ax.set_xlim(0.15, 1.85)
        ax.set_ylabel("Peak |r|", fontsize=FONT_SIZE_LABEL)
        return

    _phase_locking_plot_half_violin_R(ax, r_rest, r_run, positions, colors)

    all_y = np.concatenate([r_rest, r_run])
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    y_range = y_max - y_min
    if y_range < 1e-10:
        y_range = max(abs(y_max), 0.05) * 0.3

    test_result = _perform_paired_test(r_rest, r_run)
    p_val = test_result.get("p_value", np.nan)
    test_name = test_result.get("test_used", "n/a")
    cohens_d = test_result.get("effect_size", np.nan)

    bracket_y = y_max + 0.06 * y_range
    line_height = 0.015 * y_range
    text_offset = 0.005 * y_range

    if np.isfinite(p_val):
        _add_significance_bracket(
            ax, positions[0], positions[1], bracket_y, p_val,
            line_height=line_height, text_offset=text_offset,
        )
        top_pad = bracket_y + 0.12 * y_range
    else:
        top_pad = y_max + 0.15 * y_range

    ax.set_ylim(y_min - y_range * 0.40, top_pad)

    mean_rest = float(np.mean(r_rest))
    mean_run = float(np.mean(r_run))
    ax.set_ylabel("Peak |r|", fontsize=FONT_SIZE_LABEL)
    ax.set_title("Peak |r|, paired REST vs RUN",
                  fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
    ax.set_xlim(0.15, 1.85)
    style_axis(ax)

    p_str = (
        f"p = {p_val:.4f}" if np.isfinite(p_val) and p_val >= 0.001
        else f"p = {p_val:.2e}" if np.isfinite(p_val)
        else "p = n/a"
    )
    d_str = f"d = {cohens_d:.2f}" if np.isfinite(cohens_d) else "d = n/a"
    stats_txt = (
        f"Mean  REST = {mean_rest:.3f}\n"
        f"Mean  RUN  = {mean_run:.3f}\n"
        f"N paired = {len(r_rest)}\n"
        f"{test_name}, {p_str}\n"
        f"Cohen's {d_str}"
    )
    ax.text(
        0.02, 0.98, stats_txt, transform=ax.transAxes,
        fontsize=FONT_SIZE_TICK - 1, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="0.55", alpha=0.92),
    )


def fig_xcorr(data, animal, session, trial_num, band=THETA_BAND):
    """
    Extended cross-correlation figure:
    Top row:    [Single REST hm | Single RUN hm | cb] [All-trials REST hm | All-trials RUN hm | cb]
    Bottom row: [Single-trial mean correlogram | All-trials mean+SEM | Peak |r| violin + stats]
    """
    # ---- single-trial computation (example trial) ----
    time_vec = data["time"]
    speed = data["speed"]
    out = compute_xcorr_fiber_by_behavior(
        data["fiber1"], data["fiber2"], data["fs"], time_vec, speed, band=band,
    )
    lags_ms = out["lags_ms"]
    mat_r = out["mat_rest"]
    mat_rn = out["mat_run"]
    c_r = out["centers_rest"]
    c_rn = out["centers_run"]

    # Colour limits for single-trial heatmaps (inferno, symmetric around 0)
    mats = [m for m in (mat_r, mat_rn) if m is not None and m.size > 0]
    if mats:
        combined_abs = np.abs(np.concatenate([m.ravel() for m in mats]))
        vmax_s = min(0.5, float(np.nanpercentile(combined_abs, 98)))
        vmax_s = max(0.2, np.round(vmax_s, 1))
        vmin_s = -vmax_s
    else:
        vmin_s, vmax_s = -0.3, 0.3

    # ---- all-trials computation ----
    xcorr_session = None if XCORR_AGGREGATE_ALL_SESSIONS else session
    print("  Collecting cross-correlations across all trials ...")
    all_trials = collect_trialwise_xcorr(animal, session=xcorr_session, band=band)

    # Colour limits for all-trials heatmaps (RdBu_r, symmetric around 0)
    vmin_a, vmax_a = -0.3, 0.3
    if all_trials is not None:
        all_mats = [
            m for m in (all_trials["mat_all_rest"], all_trials["mat_all_run"])
            if m is not None and m.size > 0
        ]
        if all_mats:
            all_abs = np.abs(np.concatenate([m.ravel() for m in all_mats]))
            vmax_a = min(0.5, float(np.nanpercentile(all_abs, 98)))
            vmax_a = max(0.15, np.round(vmax_a, 2))
            vmin_a = -vmax_a

    # ---- figure layout (nested GridSpec) ----
    fig = plt.figure(figsize=FIGSIZE_XCORR_INCH)
    gs_outer = GridSpec(
        2, 1, figure=fig, height_ratios=[1.25, 1.0], hspace=0.32,
        left=0.05, right=0.97, top=0.92, bottom=0.07,
    )
    gs_top = gs_outer[0, :].subgridspec(
        1, 7, width_ratios=[1, 1, 0.04, 0.10, 1, 1, 0.04], wspace=0.12,
    )
    gs_bot = gs_outer[1, :].subgridspec(
        1, 3, width_ratios=[1.0, 1.0, 0.82], wspace=0.28,
    )

    # ======== TOP ROW: 4 heatmaps ========

    # --- single-trial REST / RUN (inferno) ---
    ax_hm_sr = fig.add_subplot(gs_top[0])
    ax_hm_sn = fig.add_subplot(gs_top[1])
    cax_s = fig.add_subplot(gs_top[2])

    im_sr = _plot_xcorr_heatmap_panel(
        ax_hm_sr, mat_r, c_r, lags_ms, vmin_s, vmax_s,
        f"Example trial — REST (n={out['n_windows_rest']} win)",
    )
    im_sn = _plot_xcorr_heatmap_panel(
        ax_hm_sn, mat_rn, c_rn, lags_ms, vmin_s, vmax_s,
        f"Example trial — RUN (n={out['n_windows_run']} win)",
    )
    im_s_cb = im_sn if im_sn is not None else im_sr
    if im_s_cb is not None:
        cb_s = fig.colorbar(im_s_cb, cax=cax_s, orientation="vertical")
        cb_s.set_label("r", fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
        cb_s.ax.tick_params(labelsize=FONT_SIZE_TICK - 2, width=TICK_WIDTH,
                            length=max(TICK_LENGTH - 2, 4))
    else:
        cax_s.set_visible(False)

    # --- gap column (hidden) ---
    ax_gap = fig.add_subplot(gs_top[3])
    ax_gap.set_visible(False)

    # --- all-trials REST / RUN (RdBu_r) ---
    ax_hm_ar = fig.add_subplot(gs_top[4])
    ax_hm_an = fig.add_subplot(gs_top[5])
    cax_a = fig.add_subplot(gs_top[6])

    n_rest_tr = all_trials["n_rest_trials"] if all_trials else 0
    n_run_tr = all_trials["n_run_trials"] if all_trials else 0
    mat_ar = all_trials["mat_all_rest"] if all_trials else None
    mat_an = all_trials["mat_all_run"] if all_trials else None
    lags_a = all_trials["lags_ms"] if all_trials else lags_ms

    im_ar = _plot_xcorr_all_trials_heatmap(
        ax_hm_ar, mat_ar, lags_a, n_rest_tr, vmin_a, vmax_a,
        XCORR_ALL_TRIALS_CMAP, f"All trials — REST (n={n_rest_tr})",
    )
    im_an = _plot_xcorr_all_trials_heatmap(
        ax_hm_an, mat_an, lags_a, n_run_tr, vmin_a, vmax_a,
        XCORR_ALL_TRIALS_CMAP, f"All trials — RUN (n={n_run_tr})",
    )
    im_a_cb = im_an if im_an is not None else im_ar
    if im_a_cb is not None:
        cb_a = fig.colorbar(im_a_cb, cax=cax_a, orientation="vertical")
        cb_lbl = {"zscore": "z-score", "peak": "norm. r"}.get(XCORR_ALL_TRIALS_NORMALIZE, "r")
        cb_a.set_label(cb_lbl, fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
        cb_a.ax.tick_params(labelsize=FONT_SIZE_TICK - 2, width=TICK_WIDTH,
                            length=max(TICK_LENGTH - 2, 4))
    else:
        cax_a.set_visible(False)

    # ======== BOTTOM ROW ========

    rust = tuple(np.clip(COLOR_XCORR_RUST, 0, 1))
    peak_lim = _peak_lag_search_half_width_ms(band)

    # helper: style a correlogram axis
    def _style_corr_ax(ax):
        ax.axvline(0, color="grey", ls="--", lw=1)
        ax.axhline(0, color="grey", ls="-", lw=0.5, alpha=0.5)
        ax.set_xlim(-XCORR_MAX_LAG_MS, XCORR_MAX_LAG_MS)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.set_xlabel("Temporal offset (ms)", fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel("Cross-correlation", fontsize=FONT_SIZE_LABEL)
        style_axis(ax)
        ax.axvline(-peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7, zorder=0)
        ax.axvline(peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7, zorder=0)

    # --- Bottom-left: single-trial mean correlogram ---
    ax_line_s = fig.add_subplot(gs_bot[0])
    _style_corr_ax(ax_line_s)
    legend_s = []
    stats_s = []

    if mat_r is not None and mat_r.size > 0:
        mean_rest_s = mean_xcorr_over_windows(mat_r)
        ln, = ax_line_s.plot(lags_ms, mean_rest_s, "-", color=rust,
                             lw=LINE_WIDTH_XCORR, label="REST")
        legend_s.append(ln)
        ir_s, _ = peak_lag_index_restricted(mean_rest_s, lags_ms, band)
        stats_s.append(
            f"REST |r|={abs(mean_rest_s[ir_s]):.3f} @ {lags_ms[ir_s]:.1f} ms"
        )
    if mat_rn is not None and mat_rn.size > 0:
        mean_run_s = mean_xcorr_over_windows(mat_rn)
        ln, = ax_line_s.plot(lags_ms, mean_run_s, "--", color=rust,
                             lw=LINE_WIDTH_XCORR, label="RUN")
        legend_s.append(ln)
        irn_s, _ = peak_lag_index_restricted(mean_run_s, lags_ms, band)
        stats_s.append(
            f"RUN |r|={abs(mean_run_s[irn_s]):.3f} @ {lags_ms[irn_s]:.1f} ms"
        )
    if legend_s:
        ax_line_s.legend(handles=legend_s, loc="upper right", frameon=True,
                         fontsize=FONT_SIZE_LEGEND, framealpha=0.95, edgecolor="0.5")
    if stats_s:
        ax_line_s.text(
            0.02, 0.95, "\n".join(stats_s), transform=ax_line_s.transAxes,
            fontsize=FONT_SIZE_TICK - 2, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.92),
        )
    ax_line_s.set_title(
        f"Example trial ({session} T{trial_num}): REST solid / RUN dashed",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    # --- Bottom-center: all-trials grand mean + SEM ---
    ax_line_a = fig.add_subplot(gs_bot[1])
    _style_corr_ax(ax_line_a)
    legend_a = []
    stats_a = []

    if all_trials is not None:
        gm_rest = all_trials["grand_mean_rest"]
        gm_run = all_trials["grand_mean_run"]
        sem_rest = all_trials["sem_rest"]
        sem_run = all_trials["sem_run"]
        lag_a = all_trials["lags_ms"]

        if gm_rest is not None:
            ln, = ax_line_a.plot(lag_a, gm_rest, "-", color=rust,
                                 lw=LINE_WIDTH_XCORR + 0.5, label="REST")
            legend_a.append(ln)
            if sem_rest is not None:
                ax_line_a.fill_between(
                    lag_a, gm_rest - sem_rest, gm_rest + sem_rest,
                    color=rust, alpha=XCORR_GRAND_MEAN_SEM_ALPHA,
                )
            ir_a, _ = peak_lag_index_restricted(gm_rest, lag_a, band)
            stats_a.append(
                f"REST |r|={abs(gm_rest[ir_a]):.3f} @ {lag_a[ir_a]:.1f} ms "
                f"(n={all_trials['n_rest_trials']})"
            )

        if gm_run is not None:
            ln, = ax_line_a.plot(lag_a, gm_run, "--", color=rust,
                                 lw=LINE_WIDTH_XCORR + 0.5, label="RUN")
            legend_a.append(ln)
            if sem_run is not None:
                ax_line_a.fill_between(
                    lag_a, gm_run - sem_run, gm_run + sem_run,
                    color=rust, alpha=XCORR_GRAND_MEAN_SEM_ALPHA,
                )
            irn_a, _ = peak_lag_index_restricted(gm_run, lag_a, band)
            stats_a.append(
                f"RUN |r|={abs(gm_run[irn_a]):.3f} @ {lag_a[irn_a]:.1f} ms "
                f"(n={all_trials['n_run_trials']})"
            )

    if legend_a:
        ax_line_a.legend(handles=legend_a, loc="upper right", frameon=True,
                         fontsize=FONT_SIZE_LEGEND, framealpha=0.95, edgecolor="0.5")
    if stats_a:
        ax_line_a.text(
            0.02, 0.95, "\n".join(stats_a), transform=ax_line_a.transAxes,
            fontsize=FONT_SIZE_TICK - 2, ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.92),
        )
    n_all = all_trials["n_loaded"] if all_trials else 0
    ax_line_a.set_title(
        f"All trials (n={n_all}): mean +/- SEM",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    # --- Bottom-right: violin (peak |r|, REST vs RUN) ---
    ax_violin = fig.add_subplot(gs_bot[2])
    _plot_xcorr_peak_r_violin_panel(ax_violin, all_trials)

    # ---- suptitle ----
    n_loaded = all_trials["n_loaded"] if all_trials else 0
    scope_lbl = (
        f"example: {session} T{trial_num}; all trials: {n_loaded}"
        if XCORR_AGGREGATE_ALL_SESSIONS
        else f"{session} T{trial_num}"
    )
    fig.suptitle(
        f"{animal} — Fiber-Fiber theta cross-correlation ({scope_lbl})",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.97,
    )
    return fig


# =============================================================================
# FIGURE 3b: UNSPLIT CROSS-CORRELATION (all epochs) + SURROGATE TEST
# =============================================================================

def fig_xcorr_unsplit(data, animal, session, trial_num, band=THETA_BAND):
    """
    Cross-correlation figure without REST/RUN separation.
    Top row:  [Single-trial time×lag heatmap | cb] [All-trials trial×lag heatmap | cb]
    Bottom:   [Single-trial mean correlogram | All-trials mean±SEM | Surrogate test panel]
    """
    # ---- single-trial computation ----
    out = compute_xcorr_fiber_unsplit(
        data["fiber1"], data["fiber2"], data["fs"], band=band,
    )
    lags_ms = out["lags_ms"]
    mat = out["mat"]
    centers = out["centers"]

    mats = [mat] if mat is not None and mat.size > 0 else []
    if mats:
        combined_abs = np.abs(np.concatenate([m.ravel() for m in mats]))
        vmax_s = min(0.5, float(np.nanpercentile(combined_abs, 98)))
        vmax_s = max(0.2, np.round(vmax_s, 1))
        vmin_s = -vmax_s
    else:
        vmin_s, vmax_s = -0.3, 0.3

    # ---- all-trials + surrogates ----
    xcorr_sess = None if XCORR_AGGREGATE_ALL_SESSIONS else session
    print("  Collecting unsplit cross-correlations + surrogates across all trials ...")
    all_trials = collect_trialwise_xcorr_unsplit(animal, session=xcorr_sess, band=band)

    vmin_a, vmax_a = -0.3, 0.3
    if all_trials is not None and all_trials["mat_all"] is not None:
        all_abs = np.abs(all_trials["mat_all"].ravel())
        vmax_a = min(0.5, float(np.nanpercentile(all_abs, 98)))
        vmax_a = max(0.15, np.round(vmax_a, 2))
        vmin_a = -vmax_a

    # ---- figure layout ----
    fig = plt.figure(figsize=FIGSIZE_XCORR_UNSPLIT_INCH)
    gs_outer = GridSpec(
        2, 1, figure=fig, height_ratios=[1.25, 1.0], hspace=0.32,
        left=0.05, right=0.97, top=0.92, bottom=0.07,
    )

    gs_top = gs_outer[0, :].subgridspec(
        1, 5, width_ratios=[1.2, 0.04, 0.10, 1, 0.04], wspace=0.12,
    )
    gs_bot = gs_outer[1, :].subgridspec(
        1, 3, width_ratios=[1.0, 1.0, 0.82], wspace=0.28,
    )

    # ======== TOP ROW ========

    # --- single-trial heatmap (time × lag) ---
    ax_hm_s = fig.add_subplot(gs_top[0])
    cax_s = fig.add_subplot(gs_top[1])

    im_s = _plot_xcorr_heatmap_panel(
        ax_hm_s, mat, centers, lags_ms, vmin_s, vmax_s,
        f"Example trial ({session} T{trial_num}, n={out['n_windows']} win)",
    )
    if im_s is not None:
        cb_s = fig.colorbar(im_s, cax=cax_s, orientation="vertical")
        cb_s.set_label("r", fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
        cb_s.ax.tick_params(labelsize=FONT_SIZE_TICK - 2, width=TICK_WIDTH,
                            length=max(TICK_LENGTH - 2, 4))
    else:
        cax_s.set_visible(False)

    ax_gap = fig.add_subplot(gs_top[2])
    ax_gap.set_visible(False)

    # --- all-trials heatmap (trial × lag) ---
    ax_hm_a = fig.add_subplot(gs_top[3])
    cax_a = fig.add_subplot(gs_top[4])

    n_loaded = all_trials["n_loaded"] if all_trials else 0
    mat_all = all_trials["mat_all"] if all_trials else None
    lags_a = all_trials["lags_ms"] if all_trials else lags_ms

    im_a = _plot_xcorr_all_trials_heatmap(
        ax_hm_a, mat_all, lags_a, n_loaded, vmin_a, vmax_a,
        XCORR_ALL_TRIALS_CMAP, f"All trials (n={n_loaded})",
    )
    if im_a is not None:
        cb_a = fig.colorbar(im_a, cax=cax_a, orientation="vertical")
        cb_lbl = {"zscore": "z-score", "peak": "norm. r"}.get(XCORR_ALL_TRIALS_NORMALIZE, "r")
        cb_a.set_label(cb_lbl, fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
        cb_a.ax.tick_params(labelsize=FONT_SIZE_TICK - 2, width=TICK_WIDTH,
                            length=max(TICK_LENGTH - 2, 4))
    else:
        cax_a.set_visible(False)

    # ======== BOTTOM ROW ========

    rust = tuple(np.clip(COLOR_XCORR_RUST, 0, 1))
    peak_lim = _peak_lag_search_half_width_ms(band)

    def _style_corr_ax(ax):
        ax.axvline(0, color="grey", ls="--", lw=1)
        ax.axhline(0, color="grey", ls="-", lw=0.5, alpha=0.5)
        ax.set_xlim(-XCORR_MAX_LAG_MS, XCORR_MAX_LAG_MS)
        ax.xaxis.set_major_locator(MultipleLocator(50))
        ax.set_xlabel("Temporal offset (ms)", fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel("Cross-correlation", fontsize=FONT_SIZE_LABEL)
        style_axis(ax)
        ax.axvline(-peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7, zorder=0)
        ax.axvline(peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7, zorder=0)

    # --- single-trial mean correlogram ---
    ax_line_s = fig.add_subplot(gs_bot[0])
    _style_corr_ax(ax_line_s)

    if mat is not None and mat.size > 0:
        mean_all_s = mean_xcorr_over_windows(mat)
        ax_line_s.plot(lags_ms, mean_all_s, "-", color=rust, lw=LINE_WIDTH_XCORR)
        ir_s, _ = peak_lag_index_restricted(mean_all_s, lags_ms, band)
        ax_line_s.text(
            0.02, 0.95,
            f"|r|={abs(mean_all_s[ir_s]):.3f} @ {lags_ms[ir_s]:.1f} ms (n={out['n_windows']} win)",
            transform=ax_line_s.transAxes, fontsize=FONT_SIZE_TICK - 2,
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.92),
        )
    ax_line_s.set_title(
        f"Example trial ({session} T{trial_num})",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    # --- all-trials grand mean ± SEM ---
    ax_line_a = fig.add_subplot(gs_bot[1])
    _style_corr_ax(ax_line_a)

    if all_trials is not None and all_trials["grand_mean"] is not None:
        gm = all_trials["grand_mean"]
        sem_a = all_trials["sem"]
        lag_a = all_trials["lags_ms"]
        ax_line_a.plot(lag_a, gm, "-", color=rust, lw=LINE_WIDTH_XCORR + 0.5)
        if sem_a is not None:
            ax_line_a.fill_between(
                lag_a, gm - sem_a, gm + sem_a,
                color=rust, alpha=XCORR_GRAND_MEAN_SEM_ALPHA,
            )
        ir_a, _ = peak_lag_index_restricted(gm, lag_a, band)
        ax_line_a.text(
            0.02, 0.95,
            f"|r|={abs(gm[ir_a]):.3f} @ {lag_a[ir_a]:.1f} ms (n={n_loaded})",
            transform=ax_line_a.transAxes, fontsize=FONT_SIZE_TICK - 2,
            ha="left", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6", alpha=0.92),
        )
    ax_line_a.set_title(
        f"All trials (n={n_loaded}): mean +/- SEM",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    # --- surrogate test panel: observed vs null distribution ---
    ax_surr = fig.add_subplot(gs_bot[2])
    _plot_surrogate_test_panel(ax_surr, all_trials, band)

    # ---- suptitle ----
    scope_lbl = (
        f"example: {session} T{trial_num}; all trials: {n_loaded}"
        if XCORR_AGGREGATE_ALL_SESSIONS
        else f"{session} T{trial_num}"
    )
    fig.suptitle(
        f"{animal} — Fiber-Fiber theta cross-correlation, all epochs ({scope_lbl})",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.97,
    )
    return fig


def _plot_surrogate_test_panel(ax, all_trials, band):
    """
    Right-bottom panel: histogram of pooled null distribution (circular-shift surrogates)
    overlaid with observed per-trial peak |r| values as individual markers + median line.
    """
    if all_trials is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        ax.set_title("Surrogate test", fontsize=FONT_SIZE_TITLE - 1, pad=8)
        style_axis(ax)
        return

    null_dist = all_trials["pooled_null"]
    obs_peaks = all_trials["observed_peaks"]
    pooled_p = all_trials["pooled_p"]
    n_sig = all_trials["n_sig_trials"]
    n_total = all_trials["n_loaded"]

    if len(null_dist) == 0:
        ax.text(0.5, 0.5, "No surrogates", transform=ax.transAxes,
                ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        ax.set_title("Surrogate test", fontsize=FONT_SIZE_TITLE - 1, pad=8)
        style_axis(ax)
        return

    null_color = "0.70"
    obs_color = tuple(np.clip(COLOR_XCORR_RUST, 0, 1))

    ax.hist(
        null_dist, bins=40, density=True, color=null_color,
        edgecolor="0.55", alpha=0.7, label="Null (circular shift)", zorder=1,
    )

    if len(obs_peaks) > 0:
        for i, pk in enumerate(obs_peaks):
            ax.axvline(pk, color=obs_color, ls="-", lw=1.0, alpha=0.4, zorder=2)
        median_obs = float(np.median(obs_peaks))
        ax.axvline(
            median_obs, color=obs_color, ls="-", lw=3.0, zorder=3,
            label=f"Observed median |r| = {median_obs:.3f}",
        )

        p95 = float(np.percentile(null_dist, 95))
        ax.axvline(
            p95, color="0.35", ls="--", lw=1.5, zorder=2,
            label=f"Null 95th pctl = {p95:.3f}",
        )

    ax.set_xlabel("Peak |r|", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel("Density", fontsize=FONT_SIZE_LABEL)
    ax.set_title(
        f"Circular-shift surrogate test (n={XCORR_UNSPLIT_N_SURROGATES}/trial)",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )
    ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND - 1, frameon=True,
              framealpha=0.95, edgecolor="0.5")
    style_axis(ax)

    p_str = (
        f"p = {pooled_p:.4f}" if np.isfinite(pooled_p) and pooled_p >= 0.001
        else f"p = {pooled_p:.2e}" if np.isfinite(pooled_p)
        else "p = n/a"
    )
    stats_txt = (
        f"Trials significant: {n_sig}/{n_total}\n"
        f"Pooled median vs null: {p_str}\n"
        f"N surrogates/trial: {XCORR_UNSPLIT_N_SURROGATES}"
    )
    ax.text(
        0.02, 0.98, stats_txt, transform=ax.transAxes,
        fontsize=FONT_SIZE_TICK - 2, va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.55", alpha=0.92),
    )


# =============================================================================
# FIGURE 6: FIBER–LFP THETA CROSS-CORRELATION (4 COMBINATIONS)
# =============================================================================

def _plot_fig6_within_fiber_violin(ax, ipsi_arr, contra_arr, color_ipsi,
                                   color_contra, title, ylabel="Peak |r|"):
    """
    Paired half-violin/box/scatter for one fiber: ipsi vs contra.
    Reuses the same geometry as _phase_locking_plot_half_violin_R.
    """
    positions = (0.7, 1.3)
    colors = (color_ipsi, color_contra)
    labels = ("Ipsi", "Contra")

    if ipsi_arr.size == 0:
        ax.text(0.5, 0.5, "Insufficient data", transform=ax.transAxes,
                ha="center", va="center", fontsize=FONT_SIZE_LABEL - 2,
                color="gray")
        ax.set_title(title, fontsize=FONT_SIZE_TITLE - 2, fontweight="bold",
                     pad=6)
        style_axis(ax)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK - 1)
        ax.set_xlim(0.15, 1.85)
        return

    _phase_locking_plot_half_violin_R(ax, ipsi_arr, contra_arr, positions,
                                      colors)

    all_y = np.concatenate([ipsi_arr, contra_arr])
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    y_range = max(y_max - y_min, abs(y_max) * 0.05)

    test_result = _perform_paired_test(ipsi_arr, contra_arr)
    p_val = test_result.get("p_value", np.nan)
    cohens_d = test_result.get("effect_size", np.nan)
    test_name = test_result.get("test_used", "n/a")

    bracket_y = y_max + 0.08 * y_range
    if np.isfinite(p_val):
        _add_significance_bracket(
            ax, positions[0], positions[1], bracket_y, p_val,
            line_height=0.015 * y_range, text_offset=0.005 * y_range,
        )
        top_pad = bracket_y + 0.14 * y_range
    else:
        top_pad = y_max + 0.15 * y_range

    ax.set_ylim(y_min - y_range * 0.35, top_pad)

    p_str = (f"p={p_val:.4f}" if np.isfinite(p_val) and p_val >= 0.001
             else f"p={p_val:.2e}" if np.isfinite(p_val) else "p=n/a")
    d_str = f"d={cohens_d:.2f}" if np.isfinite(cohens_d) else "d=n/a"
    ax.text(
        0.03, 0.02,
        f"n={len(ipsi_arr)}, {test_name}\n{p_str}, Cohen's {d_str}",
        transform=ax.transAxes, fontsize=FONT_SIZE_TICK - 2,
        va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.55",
                  alpha=0.92),
    )

    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 2, fontweight="bold",
                 pad=6)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK - 1)
    ax.set_xlim(0.15, 1.85)
    style_axis(ax)


def _plot_fig6_contrast_index(ax, trial_peaks, ylabel="Contrast Index"):
    """
    Per-trial laterality contrast index for each fiber.

    CI = (ipsi_r - contra_r) / (ipsi_r + contra_r)

    Bounded [-1, +1].  CI > 0 ⇒ ipsilateral preference; CI < 0 ⇒ contralateral.
    CI = 0 ⇒ no laterality preference.  The denominator normalises out the overall
    magnitude, so CI is independent of each fiber's theta SNR.

    Two distributions are shown (R fiber, L fiber) and each is tested against zero
    (one-sample t-test or Wilcoxon signed-rank).
    """
    color_r = tuple(np.clip(FIG6_COLORS["R\u2013R"], 0, 1))
    color_l = tuple(np.clip(FIG6_COLORS["L\u2013L"], 0, 1))

    ci_right, ci_left = [], []
    for tp in trial_peaks:
        rr = tp.get("R\u2013R", np.nan)
        rl = tp.get("R\u2013L", np.nan)
        ll = tp.get("L\u2013L", np.nan)
        lr = tp.get("L\u2013R", np.nan)

        if np.isfinite(rr) and np.isfinite(rl) and (rr + rl) > 1e-12:
            ci_right.append((rr - rl) / (rr + rl))
        if np.isfinite(ll) and np.isfinite(lr) and (ll + lr) > 1e-12:
            ci_left.append((ll - lr) / (ll + lr))

    ci_r = np.array(ci_right)
    ci_l = np.array(ci_left)

    ax.axhline(0, color="0.55", ls="--", lw=1.2, zorder=0)

    positions = [0.7, 1.3]
    datasets = [ci_r, ci_l]
    colors = [color_r, color_l]
    labels = ["R Fiber", "L Fiber"]
    rng = np.random.default_rng(42)

    for pos, data_arr, color, label in zip(positions, datasets, colors, labels):
        if data_arr.size < 2:
            continue

        try:
            kde = gaussian_kde(data_arr, bw_method=0.5)
            padding = max(float(np.ptp(data_arr)) * 0.25,
                          float(np.std(data_arr)) * 0.4)
            y_range = np.linspace(np.min(data_arr) - padding,
                                  np.max(data_arr) + padding, 100)
            density = kde(y_range)
            density = density / np.max(density) * 0.22
            ax.fill_betweenx(y_range, pos - density, pos + density,
                             alpha=0.35, color=color, zorder=1)
        except Exception:
            pass

        bp = ax.boxplot(
            [data_arr], positions=[pos], widths=0.12,
            patch_artist=True, showfliers=False, zorder=3,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor(color)
            patch.set_linewidth(2)
        for w in bp["whiskers"]:
            w.set_color(color)
            w.set_linewidth(2)
        for c_el in bp["caps"]:
            c_el.set_color(color)
            c_el.set_linewidth(2)
        for med in bp["medians"]:
            med.set_color("black")
            med.set_linewidth(2.5)

        jitter = rng.uniform(-0.04, 0.04, len(data_arr))
        ax.scatter(pos + jitter, data_arr, s=50, c=[color],
                   edgecolors="white", linewidths=1.2, zorder=5, alpha=0.85)

    stats_lines = []
    for data_arr, label in zip(datasets, labels):
        if data_arr.size >= 3:
            try:
                _, p_norm = shapiro(data_arr)
                is_normal = p_norm > 0.05
            except Exception:
                is_normal = False

            if is_normal or data_arr.size >= 7:
                stat, p_val = ttest_1samp(data_arr, 0)
                test_name = "t-test"
            else:
                stat, p_val = wilcoxon(data_arr, alternative="two-sided")
                test_name = "wilcoxon"

            mean_ci = float(np.mean(data_arr))
            p_str = (f"p={p_val:.4f}" if p_val >= 0.001
                     else f"p={p_val:.2e}")
            stats_lines.append(
                f"{label}: CI={mean_ci:+.3f}, {test_name} {p_str}"
            )
        else:
            stats_lines.append(f"{label}: n<3")

    if stats_lines:
        ax.text(
            0.03, 0.02, "\n".join(stats_lines),
            transform=ax.transAxes, fontsize=FONT_SIZE_TICK - 3,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.55",
                      alpha=0.92),
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK - 1)
    ax.set_xlim(0.15, 1.85)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1)
    ax.set_title("Laterality CI\n(ipsi\u2212contra)/(ipsi+contra)",
                 fontsize=FONT_SIZE_TITLE - 3, fontweight="bold", pad=6)
    style_axis(ax)


def _plot_fiber_lfp_xcorr_comparison_panels(fig, gs_comp, all_trials, band):
    """
    Three-panel comparison for Figure 6:
      [0] R Fiber: ipsi vs contra (paired half-violin)
      [1] L Fiber: ipsi vs contra (paired half-violin)
      [2] Laterality Contrast Index (both fibers, tested against zero)
    """
    if all_trials is None:
        for i in range(3):
            ax = fig.add_subplot(gs_comp[i])
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                    color="gray")
            style_axis(ax)
        return

    trial_peaks = all_trials["trial_peaks"]

    def _paired_arrays(ipsi_key, contra_key):
        ipsi, contra = [], []
        for tp in trial_peaks:
            v_i = tp.get(ipsi_key, np.nan)
            v_c = tp.get(contra_key, np.nan)
            if np.isfinite(v_i) and np.isfinite(v_c):
                ipsi.append(v_i)
                contra.append(v_c)
        return np.array(ipsi), np.array(contra)

    ax0 = fig.add_subplot(gs_comp[0])
    r_ipsi, r_contra = _paired_arrays("R\u2013R", "R\u2013L")
    _plot_fig6_within_fiber_violin(
        ax0, r_ipsi, r_contra,
        color_ipsi=tuple(np.clip(FIG6_COLORS["R\u2013R"], 0, 1)),
        color_contra=tuple(np.clip(FIG6_COLORS["R\u2013L"], 0, 1)),
        title="R Fiber: Ipsi vs Contra",
    )

    ax1 = fig.add_subplot(gs_comp[1])
    l_ipsi, l_contra = _paired_arrays("L\u2013L", "L\u2013R")
    _plot_fig6_within_fiber_violin(
        ax1, l_ipsi, l_contra,
        color_ipsi=tuple(np.clip(FIG6_COLORS["L\u2013L"], 0, 1)),
        color_contra=tuple(np.clip(FIG6_COLORS["L\u2013R"], 0, 1)),
        title="L Fiber: Ipsi vs Contra",
    )

    ax2 = fig.add_subplot(gs_comp[2])
    _plot_fig6_contrast_index(ax2, trial_peaks)


def fig_fiber_lfp_xcorr(animal, session, trial_num, band=THETA_BAND):
    """
    Figure 6: Fiber–LFP theta cross-correlation for 4 combinations.

    Layout – 4 columns (R–R | R–L | L–L | L–R) + 1 comparison column:
      top row:    example-trial heatmap (Time on Y, Temporal offset on X)
      bottom row: all-trials grand mean ± SEM correlogram (shared X axis)
      right col:  grouped violin/box/scatter for peak |r| across combinations
    """
    data = load_trial(animal, session, trial_num)
    fs = float(data["fs"])

    single_trial = {}
    for combo in FIG6_COMBINATIONS:
        sig1 = np.asarray(data[combo["sig1_key"]], dtype=float).ravel()
        sig2 = data.get(combo["sig2_key"])
        if sig2 is None:
            single_trial[combo["short"]] = None
            continue
        sig2 = np.asarray(sig2, dtype=float).ravel()
        try:
            out = compute_xcorr_fiber_unsplit(sig1, sig2, fs, band=band)
            single_trial[combo["short"]] = out
        except Exception:
            single_trial[combo["short"]] = None

    valid_mats = [r["mat"] for r in single_trial.values()
                  if r is not None and r["mat"] is not None and r["mat"].size > 0]
    if valid_mats:
        combined_abs = np.abs(np.concatenate([m.ravel() for m in valid_mats]))
        vmax_s = min(0.5, float(np.nanpercentile(combined_abs, 98)))
        vmax_s = max(0.15, np.round(vmax_s, 2))
        vmin_s = -vmax_s
    else:
        vmin_s, vmax_s = -0.3, 0.3

    sess_arg = None if FIG6_AGGREGATE_ALL_SESSIONS else session
    print("  Collecting fiber\u2013LFP cross-correlations across all trials ...")
    all_trials = collect_trialwise_fiber_lfp_xcorr(animal, session=sess_arg, band=band)

    fig = plt.figure(figsize=FIGSIZE_FIG6_INCH)
    gs = GridSpec(
        2, 6, figure=fig,
        height_ratios=[1.3, 1.0],
        width_ratios=[1, 1, 1, 1, 0.04, 1.2],
        wspace=0.28, hspace=0.12,
        left=0.05, right=0.97, top=0.90, bottom=0.08,
    )

    peak_lim = _peak_lag_search_half_width_ms(band)
    im_ref = None

    for i, combo in enumerate(FIG6_COMBINATIONS):
        ax_hm = fig.add_subplot(gs[0, i])
        ax_corr = fig.add_subplot(gs[1, i], sharex=ax_hm)
        color = tuple(np.clip(FIG6_COLORS[combo["short"]], 0, 1))

        out = single_trial[combo["short"]]
        if out is not None and out["mat"] is not None and out["mat"].size > 0:
            mat_T = out["mat"].T
            extent = [float(out["lags_ms"][0]), float(out["lags_ms"][-1]),
                      float(out["centers"][0]), float(out["centers"][-1])]
            im = ax_hm.imshow(
                mat_T, aspect="auto", origin="lower", interpolation="bilinear",
                cmap=FIG6_HEATMAP_CMAP, rasterized=True, extent=extent,
                vmin=vmin_s, vmax=vmax_s,
            )
            ax_hm.axvline(0, color="white", ls="--", lw=1.2, alpha=0.7)
            if im_ref is None:
                im_ref = im
        else:
            ax_hm.text(0.5, 0.5, "No data", transform=ax_hm.transAxes,
                       ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                       color="gray")

        ax_hm.set_ylabel("Time (s)" if i == 0 else "", fontsize=FONT_SIZE_LABEL)
        ax_hm.set_title(combo["label"], fontsize=FONT_SIZE_TITLE - 2,
                        fontweight="bold", pad=8)
        style_axis(ax_hm)
        plt.setp(ax_hm.get_xticklabels(), visible=False)

        ax_corr.axvline(0, color="grey", ls="--", lw=1)
        ax_corr.axhline(0, color="grey", ls="-", lw=0.5, alpha=0.5)
        ax_corr.axvline(-peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7,
                        zorder=0)
        ax_corr.axvline(peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7,
                        zorder=0)

        if all_trials is not None:
            cd = all_trials["combo_data"][combo["short"]]
            gm = cd["grand_mean"]
            sem = cd["sem"]
            lags = all_trials["lags_ms"]

            if gm is not None:
                ax_corr.plot(lags, gm, "-", color=color, lw=FIG6_CORR_LINE_WIDTH)
                if sem is not None:
                    ax_corr.fill_between(lags, gm - sem, gm + sem,
                                         color=color, alpha=FIG6_SEM_ALPHA)
                ir, _ = peak_lag_index_restricted(gm, lags, band)
                ax_corr.text(
                    0.03, 0.95,
                    f"|r|={abs(gm[ir]):.3f}\n@ {lags[ir]:.1f} ms\n"
                    f"n={cd['n_trials']}",
                    transform=ax_corr.transAxes, fontsize=FONT_SIZE_TICK - 2,
                    ha="left", va="top",
                    bbox=dict(boxstyle="round,pad=0.3", fc="white",
                              ec="0.6", alpha=0.92),
                )

        ax_corr.set_xlim(-XCORR_MAX_LAG_MS, XCORR_MAX_LAG_MS)
        ax_corr.xaxis.set_major_locator(MultipleLocator(50))
        ax_corr.set_xlabel("Temporal offset (ms)", fontsize=FONT_SIZE_LABEL)
        ax_corr.set_ylabel("Cross-correlation" if i == 0 else "",
                           fontsize=FONT_SIZE_LABEL)
        style_axis(ax_corr)

    cax = fig.add_subplot(gs[0, 4])
    if im_ref is not None:
        cb = fig.colorbar(im_ref, cax=cax, orientation="vertical")
        cb.set_label("r", fontsize=FONT_SIZE_COLORBAR, fontweight="bold")
        cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 2, width=TICK_WIDTH,
                          length=max(TICK_LENGTH - 2, 4))
    else:
        cax.set_visible(False)

    ax_cax_bot = fig.add_subplot(gs[1, 4])
    ax_cax_bot.set_visible(False)

    gs_comp = gs[:, 5].subgridspec(3, 1, hspace=0.45)
    _plot_fiber_lfp_xcorr_comparison_panels(fig, gs_comp, all_trials, band)

    n_total = all_trials["n_loaded"] if all_trials else 0
    fig.suptitle(
        f"{animal} \u2014 Fiber\u2013LFP $\\theta$ cross-correlation "
        f"(example: {session} T{trial_num}, all trials: n={n_total})",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.97,
    )
    return fig


# =============================================================================
# FIGURE 7: FIBER–LFP PHASE-LOCKING (4 COMBINATIONS, ALL EPOCHS)
# =============================================================================

def collect_trialwise_fiber_lfp_phase(animal, session=None, band=THETA_BAND,
                                      edge_trim_sec=PHASE_EDGE_TRIM_SEC):
    """
    Per-trial phase-locking R for 4 fiber–LFP combinations (all epochs, no
    REST/RUN split).  For each combo the instantaneous phase difference
    Δφ = φ_fiber − φ_LFP is computed via Hilbert transform and the circular
    mean resultant length R is recorded per trial.

    Returns dict with:
      combo_data  – per-combo: phi_trials (list of arrays), R_trials, mu_trials
      trial_peaks – list of dicts {session, trial_num, combo_short: R, ...}
      n_loaded, n_failed
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    if session is not None:
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    combo_keys = [c["short"] for c in FIG6_COMBINATIONS]
    combo_data = {k: {"phi_trials": [], "R_trials": [], "mu_trials": []}
                  for k in combo_keys}
    trial_peaks = []
    n_loaded, n_failed = 0, 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue

            fs = float(data["fs"])
            n_loaded += 1
            tp = {"session": sess_key, "trial_num": trial_num}

            for combo in FIG6_COMBINATIONS:
                sig1 = data.get(combo["sig1_key"])
                sig2 = data.get(combo["sig2_key"])
                if sig1 is None or sig2 is None:
                    tp[combo["short"]] = np.nan
                    continue

                sig1 = np.asarray(sig1, dtype=float).ravel()
                sig2 = np.asarray(sig2, dtype=float).ravel()
                n = min(len(sig1), len(sig2))
                if n < 20:
                    tp[combo["short"]] = np.nan
                    continue

                sig1, sig2 = sig1[:n], sig2[:n]
                trim = max(0, int(round(edge_trim_sec * fs)))
                if n <= 2 * trim + 8:
                    tp[combo["short"]] = np.nan
                    continue
                sig1 = sig1[trim: n - trim]
                sig2 = sig2[trim: n - trim]

                try:
                    z1 = theta_analytic_from_raw(sig1, fs, band=band)
                    z2 = theta_analytic_from_raw(sig2, fs, band=band)
                except Exception:
                    tp[combo["short"]] = np.nan
                    continue

                dphi = wrap_to_pi(np.angle(z1) - np.angle(z2))
                ok = np.isfinite(dphi)
                dphi_clean = dphi[ok]

                mu, R = circular_mean_and_resultant_length(dphi_clean)

                cd = combo_data[combo["short"]]
                cd["phi_trials"].append(dphi_clean)
                cd["R_trials"].append(R)
                cd["mu_trials"].append(mu)
                tp[combo["short"]] = R

            trial_peaks.append(tp)
            combo_strs = "  ".join(
                f"{c['short']} R={tp.get(c['short'], np.nan):.4f}"
                for c in FIG6_COMBINATIONS
            )
            print(f"    {sess_key} T{trial_num}: {combo_strs}")

    if n_loaded == 0:
        return None

    for key, cd in combo_data.items():
        cd["R_values"] = np.array(cd["R_trials"])
        cd["n_trials"] = len(cd["R_trials"])
        cd["peak_r_values"] = cd["R_values"]

    print(f"  Fiber-LFP phase: {n_loaded} trials loaded, {n_failed} failed")
    return {
        "combo_data": combo_data,
        "trial_peaks": trial_peaks,
        "n_loaded": n_loaded,
        "n_failed": n_failed,
    }


def fig_fiber_lfp_phase_locking(animal, band=THETA_BAND):
    """
    Figure 7: Fiber–LFP phase-locking for 4 combinations (all epochs).

    Layout:
      4 rose plots (trial-equalized) + 3 comparison panels on the right
      (R-fiber ipsi vs contra | L-fiber ipsi vs contra | Contrast Index).
    """
    sess_arg = None if FIG7_AGGREGATE_ALL_SESSIONS else None
    print("  Collecting fiber\u2013LFP phase-locking across all trials ...")
    phase_data = collect_trialwise_fiber_lfp_phase(animal, session=sess_arg,
                                                   band=band)

    fig = plt.figure(figsize=FIGSIZE_FIG7_INCH)
    gs = GridSpec(
        1, 6, figure=fig,
        width_ratios=[1, 1, 1, 1, 0.05, 1.2],
        wspace=0.30,
        left=0.04, right=0.97, top=0.86, bottom=0.08,
    )

    for i, combo in enumerate(FIG6_COMBINATIONS):
        ax = fig.add_subplot(gs[0, i], projection="polar")
        color = tuple(np.clip(FIG6_COLORS[combo["short"]], 0, 1))

        if phase_data is not None:
            cd = phase_data["combo_data"][combo["short"]]
            phi_list = cd["phi_trials"]
        else:
            phi_list = []

        rose = _trial_equalized_rose_histogram_and_vector(
            phi_list, PHASE_ROSE_NBINS,
        )

        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        if rose is not None:
            ax.bar(
                rose["theta"], rose["bar_heights"],
                width=rose["width"] * 0.98, bottom=0.0, align="center",
                color=color, edgecolor="white",
                linewidth=PHASE_ROSE_BAR_EDGELW,
                alpha=FIG7_ROSE_BAR_ALPHA,
            )
            bar_max = float(np.max(rose["bar_heights"])) if rose["bar_heights"].size else 10.0
            if np.isfinite(rose["mu"]) and np.isfinite(rose["R"]):
                ax.annotate(
                    "", xytext=(0.0, 0.0),
                    xy=(rose["mu"], rose["R"] * bar_max),
                    arrowprops=dict(
                        arrowstyle="-|>", color="crimson",
                        lw=PHASE_MEAN_VECTOR_LW,
                        mutation_scale=PHASE_MEAN_VECTOR_ARROW_MUTATION,
                        shrinkA=0, shrinkB=0,
                    ),
                    zorder=6,
                )
            ax.set_ylim(0.0, bar_max * 1.12)

            deg = (np.degrees(rose["mu"]) if np.isfinite(rose["mu"])
                   else float("nan"))
            n_total = sum(len(p) for p in phi_list)
            ttl = (
                f"{combo['label']}\n"
                f"N={n_total}, {rose['n_trials']} trials\n"
                f"Mean={deg:.1f}$^\\circ$  R={rose['R']:.3f}"
            )
        else:
            ttl = f"{combo['label']}\n(no data)"

        _style_phase_rose_polar_axis(ax)
        ax.set_title(
            ttl, fontsize=PHASE_ROSE_TITLE_FONTSIZE - 1,
            fontweight=PHASE_ROSE_TITLE_FONTWEIGHT,
            pad=PHASE_ROSE_TITLE_PAD,
        )

    ax_spacer = fig.add_subplot(gs[0, 4])
    ax_spacer.set_visible(False)

    gs_comp = gs[0, 5].subgridspec(3, 1, hspace=0.50)

    if phase_data is not None:
        trial_peaks = phase_data["trial_peaks"]

        def _paired_arrays(ipsi_key, contra_key):
            ipsi, contra = [], []
            for tp in trial_peaks:
                v_i = tp.get(ipsi_key, np.nan)
                v_c = tp.get(contra_key, np.nan)
                if np.isfinite(v_i) and np.isfinite(v_c):
                    ipsi.append(v_i)
                    contra.append(v_c)
            return np.array(ipsi), np.array(contra)

        ax0 = fig.add_subplot(gs_comp[0])
        r_ipsi, r_contra = _paired_arrays("R\u2013R", "R\u2013L")
        _plot_fig6_within_fiber_violin(
            ax0, r_ipsi, r_contra,
            color_ipsi=tuple(np.clip(FIG6_COLORS["R\u2013R"], 0, 1)),
            color_contra=tuple(np.clip(FIG6_COLORS["R\u2013L"], 0, 1)),
            title="R Fiber: Ipsi vs Contra",
            ylabel="Phase-locking R",
        )

        ax1 = fig.add_subplot(gs_comp[1])
        l_ipsi, l_contra = _paired_arrays("L\u2013L", "L\u2013R")
        _plot_fig6_within_fiber_violin(
            ax1, l_ipsi, l_contra,
            color_ipsi=tuple(np.clip(FIG6_COLORS["L\u2013L"], 0, 1)),
            color_contra=tuple(np.clip(FIG6_COLORS["L\u2013R"], 0, 1)),
            title="L Fiber: Ipsi vs Contra",
            ylabel="Phase-locking R",
        )

        ax2 = fig.add_subplot(gs_comp[2])
        _plot_fig6_contrast_index(ax2, trial_peaks, ylabel="Laterality CI")
    else:
        for j in range(3):
            ax_e = fig.add_subplot(gs_comp[j])
            ax_e.text(0.5, 0.5, "No data", transform=ax_e.transAxes,
                      ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                      color="gray")
            style_axis(ax_e)

    n_total = phase_data["n_loaded"] if phase_data else 0
    scope = "all sessions" if FIG7_AGGREGATE_ALL_SESSIONS else ""
    fig.suptitle(
        f"{animal} \u2014 Fiber\u2013LFP $\\theta$ phase-locking, "
        f"all epochs ({scope}, n={n_total} trials)",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.97,
    )
    return fig


# =============================================================================
# LFP–LFP COHERENCE DIAGNOSTIC
# =============================================================================

def collect_trialwise_lfp_lfp_coherence(animal, session=None):
    """
    Per-trial coherence spectrum between R LFP and L LFP.
    Returns list of (freq, coh) tuples and the number of loaded/failed trials.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    if session is not None:
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    trials = []
    freq_ref = None
    n_failed = 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue

            lr = data.get("lfp_right")
            ll = data.get("lfp_left")
            if lr is None or ll is None:
                n_failed += 1
                continue

            lr = np.asarray(lr, dtype=float).ravel()
            ll = np.asarray(ll, dtype=float).ravel()
            fs = float(data["fs"])

            try:
                freq, coh = compute_coherence_welch(lr, ll, fs)
            except Exception:
                n_failed += 1
                continue

            if freq_ref is None:
                freq_ref = freq
            trials.append({
                "session": sess_key, "trial_num": trial_num,
                "freq": freq, "coh": coh,
            })

    if not trials or freq_ref is None:
        return None

    coh_mat = np.vstack([t["coh"] for t in trials])
    grand_mean = np.nanmean(coh_mat, axis=0)
    n = len(trials)
    sem = np.nanstd(coh_mat, axis=0, ddof=1) / np.sqrt(n) if n > 1 else None

    theta_mask = (freq_ref >= THETA_BAND[0]) & (freq_ref <= THETA_BAND[1])
    theta_coh_per_trial = [float(np.nanmean(t["coh"][theta_mask]))
                           for t in trials]

    return {
        "freq": freq_ref, "trials": trials, "coh_mat": coh_mat,
        "grand_mean": grand_mean, "sem": sem,
        "theta_coh_per_trial": np.array(theta_coh_per_trial),
        "n_loaded": n, "n_failed": n_failed,
    }


def fig_lfp_lfp_coherence_diagnostic(animal):
    """
    Diagnostic figure: R LFP – L LFP coherence spectrum.
    Panel 1: Single example trial.
    Panel 2: All-trials mean ± SEM with theta-band shading.
    """
    print("  Collecting LFP\u2013LFP coherence across all trials ...")
    sess_arg = (None if FIG_LFP_DIAG_AGGREGATE_ALL_SESSIONS
                else FIG_LFP_DIAG_EXAMPLE_SESSION)
    coh_data = collect_trialwise_lfp_lfp_coherence(animal, session=sess_arg)

    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE_LFP_DIAG_INCH,
                             gridspec_kw={"wspace": 0.30})

    coh_color = tuple(np.clip(COLOR_LFP_COH, 0, 1))

    def _style(ax):
        ax.set_xlim(*LFP_DIAG_FREQ_RANGE)
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL)
        ax.set_ylabel("Coherence", fontsize=FONT_SIZE_LABEL)
        ax.axvspan(THETA_BAND[0], THETA_BAND[1], alpha=0.10,
                   color="gold", zorder=0,
                   label=f"$\\theta$ ({THETA_BAND[0]}\u2013{THETA_BAND[1]} Hz)")
        ax.legend(fontsize=FONT_SIZE_LEGEND, loc="upper right")
        style_axis(ax)

    ax0, ax1 = axes

    if coh_data is None:
        for ax in axes:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                    color="gray")
            style_axis(ax)
        fig.suptitle(f"{animal} \u2014 R LFP vs L LFP coherence (no data)",
                     fontsize=FONT_SIZE_SUPTITLE, fontweight="bold")
        return fig

    ex = coh_data["trials"][0]
    ax0.plot(ex["freq"], ex["coh"], "-", color=coh_color, lw=2.5)
    _style(ax0)
    theta_mask = ((ex["freq"] >= THETA_BAND[0])
                  & (ex["freq"] <= THETA_BAND[1]))
    theta_val = float(np.nanmean(ex["coh"][theta_mask]))
    ax0.set_title(
        f"Example trial ({ex['session']} T{ex['trial_num']})\n"
        f"$\\theta$ coherence = {theta_val:.3f}",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    gm = coh_data["grand_mean"]
    sem = coh_data["sem"]
    freq = coh_data["freq"]
    n = coh_data["n_loaded"]

    ax1.plot(freq, gm, "-", color=coh_color, lw=2.5)
    if sem is not None:
        ax1.fill_between(freq, gm - sem, gm + sem, color=coh_color,
                         alpha=0.25)
    _style(ax1)
    theta_mask_g = (freq >= THETA_BAND[0]) & (freq <= THETA_BAND[1])
    theta_grand = float(np.nanmean(gm[theta_mask_g]))
    theta_vals = coh_data["theta_coh_per_trial"]
    ax1.set_title(
        f"All trials (n={n}): mean \u00b1 SEM\n"
        f"$\\theta$ coherence = {theta_grand:.3f} "
        f"(range {float(np.min(theta_vals)):.3f}\u2013"
        f"{float(np.max(theta_vals)):.3f})",
        fontsize=FONT_SIZE_TITLE - 1, fontweight="bold", pad=8,
    )

    fig.suptitle(
        f"{animal} \u2014 R LFP vs L LFP coherence diagnostic",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    return fig


# =============================================================================
# FIGURE 4: PHASE-LOCKING — ROSE (REST/RUN) + EXAMPLE ΔF/F & SAWTOOTH PHASE
# =============================================================================

def fig_phase_locking(data, animal, session, trial_num, band=THETA_BAND):
    """
    Top row: trial-equalized Δφ roses (mean normalized histogram + mean of per-trial z-vectors).

    Middle: θ-band ΔF/F (%, same pipeline as phase) — 1 s example; scale bars only (no y-axis ticks).

    Bottom: wrapped θ-phase (−π…π); y-axis tick labels only — no horizontal time scale bar (time from x-axis).
    """
    phase_session = None if PHASE_AGGREGATE_ALL_SESSIONS else session
    phase_data = collect_trialwise_phase_difference(
        animal,
        phase_session,
        band=band,
        edge_trim_sec=PHASE_EDGE_TRIM_SEC,
    )
    rest_trials = [t["phi_rest"] for t in phase_data["trials"] if t["use_for_rest_stats"]]
    run_trials = [t["phi_run"] for t in phase_data["trials"] if t["use_for_run_stats"]]

    t_full = np.asarray(data["time"], dtype=float)
    fs = float(data["fs"])
    win = (t_full >= PHASE_EXAMPLE_T0_S) & (t_full <= PHASE_EXAMPLE_T1_S)
    if not np.any(win):
        warnings.warn(
            f"No samples in phase example window [{PHASE_EXAMPLE_T0_S}, {PHASE_EXAMPLE_T1_S}] s; "
            "using full trace."
        )
        win = np.ones_like(t_full, dtype=bool)
    t = t_full[win]

    f1_full = np.asarray(data["fiber1"], dtype=float)
    f2_full = np.asarray(data["fiber2"], dtype=float)
    # Full-trial prefilter → θ analytic → phase (avoids Hilbert discontinuity at window edges)
    pre1_full = preprocess_frac_for_theta(f1_full, fs)
    pre2_full = preprocess_frac_for_theta(f2_full, fs)
    z1_full = theta_analytic_signal(pre1_full, fs, band=band)
    z2_full = theta_analytic_signal(pre2_full, fs, band=band)
    phi_r_full = np.angle(z1_full)
    phi_l_full = np.angle(z2_full)
    # Real part = θ-bandpassed ΔF/F (%) prior to Hilbert; matches phase extraction
    gt1_full = np.real(z1_full)
    gt2_full = np.real(z2_full)

    gt1_w = gt1_full[win]
    gt2_w = gt2_full[win]
    phi_r_w = phi_r_full[win]
    phi_l_w = phi_l_full[win]

    fig = plt.figure(figsize=FIGSIZE_PHASE_LOCKING_INCH)
    gs = GridSpec(
        3,
        3,
        figure=fig,
        height_ratios=[1.15, 0.95, 0.95],
        width_ratios=[1.0, 1.0, 1.02],
        hspace=0.38,
        wspace=0.35,
        left=0.06,
        right=0.97,
        top=0.90,
        bottom=0.07,
    )

    ax_pr = fig.add_subplot(gs[0, 0], projection="polar")
    ax_pn = fig.add_subplot(gs[0, 1], projection="polar")
    ax_rv = fig.add_subplot(gs[0, 2])
    _plot_phase_rose_polar_trial_equalized(
        ax_pr,
        rest_trials,
        "REST: dphi (trial-equalized)",
        n_bins=PHASE_ROSE_NBINS,
        bar_color=COLOR_PHASE_ROSE_REST,
    )
    run_rose_title = (
        "RUN: dphi (env-gated, trial-eq.)"
        if PHASE_RUN_USE_ENVELOPE_GATE
        else "RUN: dphi (trial-equalized)"
    )
    _plot_phase_rose_polar_trial_equalized(
        ax_pn,
        run_trials,
        run_rose_title,
        n_bins=PHASE_ROSE_NBINS,
        bar_color=COLOR_PHASE_ROSE_RUN,
    )
    _plot_phase_locking_R_violin_panel(ax_rv, phase_data)

    # --- Middle row: theta-band dF/F traces with proper axes ---
    trace_lw = 2.2
    ax_df = fig.add_subplot(gs[1, :])
    c_r = tuple(np.clip(COLOR_RIGHT_HP, 0, 1))
    c_l = tuple(np.clip(COLOR_LEFT_HP, 0, 1))
    ax_df.plot(
        t, gt1_w, color=c_r, lw=trace_lw, alpha=0.92,
        label=f"Right HP ({band[0]:.0f}-{band[1]:.0f} Hz)",
    )
    ax_df.plot(
        t, gt2_w, color=c_l, lw=trace_lw, alpha=0.92,
        label=f"Left HP ({band[0]:.0f}-{band[1]:.0f} Hz)",
    )
    ax_df.set_ylabel(r"$\theta$ $\Delta$F/F (%)", fontsize=FONT_SIZE_LABEL)
    ax_df.set_title(
        f"Example (trial {trial_num}): theta-filtered traces",
        fontsize=FONT_SIZE_TITLE, fontweight="normal", pad=6,
    )
    style_axis(ax_df)
    ax_df.tick_params(labelbottom=False)
    ax_df.legend(
        loc="upper right", fontsize=FONT_SIZE_LEGEND,
        frameon=True, framealpha=0.95, edgecolor="0.5",
    )

    # --- Bottom row: sawtooth phase with proper axes and pi symbols ---
    phase_lw = 2.4
    ax_ph = fig.add_subplot(gs[2, :], sharex=ax_df)
    ax_ph.plot(
        t, phi_r_w, color=c_r, ls="-", lw=phase_lw, alpha=0.95,
        label="Right HP (solid)",
    )
    ax_ph.plot(
        t, phi_l_w, color=c_l, ls="--", lw=phase_lw, alpha=0.95,
        label="Left HP (dashed)",
    )
    ax_ph.axhline(0.0, color="0.8", ls="-", lw=0.6, alpha=0.5)
    ax_ph.set_ylim(-np.pi * 1.08, np.pi * 1.08)
    ax_ph.set_yticks([-np.pi, 0.0, np.pi])
    ax_ph.set_yticklabels([r"$-\pi$", r"$0$", r"$\pi$"])
    ax_ph.set_ylabel("Phase (rad)", fontsize=FONT_SIZE_LABEL)
    ax_ph.set_xlabel("Time (s)", fontsize=FONT_SIZE_LABEL)
    style_axis(ax_ph)
    ax_ph.xaxis.set_major_locator(MultipleLocator(0.2))
    ax_ph.legend(
        loc="upper right", fontsize=FONT_SIZE_LEGEND,
        frameon=True, framealpha=0.95, edgecolor="0.5",
    )

    # --- Minimal suptitle (publication standard: no metadata dump) ---
    n_ld = phase_data["n_loaded_ok"]
    n_rq = phase_data["n_rest_qualifying"]
    n_nq = phase_data["n_run_qualifying"]
    fig.suptitle(
        f"{animal} — Bilateral hippocampal theta phase-locking"
        f" ({n_ld} trials, {n_rq} REST / {n_nq} RUN)",
        fontsize=PHASE_FIG4_SUPTITLE_FONTSIZE,
        fontweight=PHASE_FIG4_SUPTITLE_FONTWEIGHT,
        y=0.97,
    )
    return fig


# =============================================================================
# MAIN
# =============================================================================


def _save_phase_locking_figure(fig, path_no_ext):
    """
    Save Figure 4 as PDF + PNG, using a pure Agg canvas and cascading DPI
    fallbacks to work around FreeType raster-overflow (0x62) on Windows.

    Root cause: FreeType's smooth rasterizer can overflow its internal cell
    buffer when a glyph's x-origin is very large (long single-line text at
    high DPI).  Wrapping long text to < 90 chars/line is the primary fix;
    DPI cascade is the safety net.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg

    path_no_ext = Path(path_no_ext)
    base_kw = dict(bbox_inches="tight", pad_inches=0.25, facecolor="white", edgecolor="none")

    original_canvas = fig.canvas
    FigureCanvasAgg(fig)

    with mpl.rc_context({"text.hinting": "none"}):
        fig.savefig(str(path_no_ext.with_suffix(".pdf")), dpi=DPI, **base_kw)

        png_path = str(path_no_ext.with_suffix(".png"))
        png_saved = False
        for try_dpi in [PHASE_FIG4_PNG_DPI, 150, 100]:
            try:
                fig.savefig(png_path, dpi=try_dpi, **base_kw)
                png_saved = True
                break
            except RuntimeError as e:
                if "raster" in str(e).lower() or "glyph" in str(e).lower():
                    warnings.warn(
                        f"PNG at {try_dpi} dpi hit FreeType raster overflow; "
                        "trying lower dpi ..."
                    )
                else:
                    raise

        if not png_saved:
            warnings.warn(
                "PNG could not be saved (FreeType raster overflow at all DPIs). "
                "PDF was saved successfully."
            )

    fig.set_canvas(original_canvas)


# =============================================================================
# FIGURE 5: BILATERAL SPECTRAL ANALYSIS (PSD, COHERENCE, COHEROGRAM)
# =============================================================================

def compute_psd_welch_db(sig, fs,
                         window_sec=SPEC5_PSD_WINDOW_SEC,
                         overlap_frac=SPEC5_PSD_OVERLAP_FRAC,
                         nfft_mult=SPEC5_PSD_NFFT_MULT):
    """Welch PSD in dB (10·log10). Input assumed % or similar continuous signal."""
    sig = np.asarray(sig, dtype=float).ravel()
    sig = sig - np.nanmean(sig)
    sig = np.nan_to_num(sig, nan=0.0)
    nperseg = max(16, int(round(window_sec * fs)))
    nperseg = min(nperseg, len(sig))
    noverlap = int(round(overlap_frac * nperseg))
    nfft = _next_pow2(nperseg * nfft_mult)
    freq, pxx = signal.welch(sig, fs=fs, nperseg=nperseg, noverlap=noverlap,
                             nfft=nfft, window="hann", scaling="density")
    pxx_db = 10.0 * np.log10(np.maximum(pxx, 1e-30))
    return freq, pxx_db


def compute_coherence_welch(sig1, sig2, fs,
                            window_sec=SPEC5_PSD_WINDOW_SEC,
                            overlap_frac=SPEC5_PSD_OVERLAP_FRAC,
                            nfft_mult=SPEC5_PSD_NFFT_MULT):
    """Magnitude-squared coherence (0–1) via scipy.signal.coherence (Welch)."""
    s1 = np.asarray(sig1, dtype=float).ravel()
    s2 = np.asarray(sig2, dtype=float).ravel()
    n = min(len(s1), len(s2))
    s1, s2 = s1[:n] - np.nanmean(s1[:n]), s2[:n] - np.nanmean(s2[:n])
    s1, s2 = np.nan_to_num(s1), np.nan_to_num(s2)
    nperseg = max(16, int(round(window_sec * fs)))
    nperseg = min(nperseg, n)
    noverlap = int(round(overlap_frac * nperseg))
    nfft = _next_pow2(nperseg * nfft_mult)
    freq, coh = signal.coherence(s1, s2, fs=fs, nperseg=nperseg, noverlap=noverlap,
                                 nfft=nfft, window="hann")
    return freq, coh


def _concat_masked_segments(sig, mask):
    """Concatenate samples where mask is True, skipping gaps."""
    sig = np.asarray(sig, dtype=float).ravel()
    mask = np.asarray(mask, dtype=bool).ravel()
    n = min(len(sig), len(mask))
    return sig[:n][mask[:n]]


def compute_psd_rest_run(sig, speed, fs, time_vec):
    """Per-condition Welch PSD in dB. Returns (freq, psd_rest_db, psd_run_db)."""
    is_rest, is_run, _ = classify_rest_run_masks(speed, fs)
    n = min(len(sig), len(is_rest))
    trim = max(0, int(round(SPEC5_EDGE_TRIM_SEC * fs)))
    if n <= 2 * trim + SPEC5_MIN_SEGMENT_SAMPLES:
        return None, None, None
    sl = slice(trim, n - trim)
    sig_t = np.asarray(sig, dtype=float).ravel()[sl]
    is_rest_t, is_run_t = is_rest[sl], is_run[sl]

    psd_rest_db = psd_run_db = None
    seg_rest = _concat_masked_segments(sig_t, is_rest_t)
    if seg_rest.size >= SPEC5_MIN_SEGMENT_SAMPLES:
        freq, psd_rest_db = compute_psd_welch_db(seg_rest, fs)
    seg_run = _concat_masked_segments(sig_t, is_run_t)
    if seg_run.size >= SPEC5_MIN_SEGMENT_SAMPLES:
        freq_r, psd_run_db = compute_psd_welch_db(seg_run, fs)
        if psd_rest_db is None:
            freq = freq_r
    if psd_rest_db is None and psd_run_db is None:
        return None, None, None
    return freq, psd_rest_db, psd_run_db


def compute_coherence_rest_run(sig1, sig2, speed, fs, time_vec):
    """Per-condition Welch coherence. Returns (freq, coh_rest, coh_run)."""
    is_rest, is_run, _ = classify_rest_run_masks(speed, fs)
    n = min(len(sig1), len(sig2), len(is_rest))
    trim = max(0, int(round(SPEC5_EDGE_TRIM_SEC * fs)))
    if n <= 2 * trim + SPEC5_MIN_SEGMENT_SAMPLES:
        return None, None, None
    sl = slice(trim, n - trim)
    s1, s2 = np.asarray(sig1, dtype=float).ravel()[sl], np.asarray(sig2, dtype=float).ravel()[sl]
    is_rest_t, is_run_t = is_rest[sl], is_run[sl]

    coh_rest = coh_run = None
    sr1, sr2 = _concat_masked_segments(s1, is_rest_t), _concat_masked_segments(s2, is_rest_t)
    if sr1.size >= SPEC5_MIN_SEGMENT_SAMPLES:
        freq, coh_rest = compute_coherence_welch(sr1, sr2, fs)
    sn1, sn2 = _concat_masked_segments(s1, is_run_t), _concat_masked_segments(s2, is_run_t)
    if sn1.size >= SPEC5_MIN_SEGMENT_SAMPLES:
        freq_r, coh_run = compute_coherence_welch(sn1, sn2, fs)
        if coh_rest is None:
            freq = freq_r
    if coh_rest is None and coh_run is None:
        return None, None, None
    return freq, coh_rest, coh_run


def compute_coherogram(sig1, sig2, fs,
                       window_sec=SPEC5_COHEROGRAM_WINDOW_SEC,
                       step_sec=SPEC5_COHEROGRAM_STEP_SEC):
    """
    Sliding-window magnitude-squared coherence → (freq, time_centers, coh_matrix).
    coh_matrix shape: (n_freq, n_windows).
    """
    s1 = np.asarray(sig1, dtype=float).ravel()
    s2 = np.asarray(sig2, dtype=float).ravel()
    n = min(len(s1), len(s2))
    s1, s2 = s1[:n], s2[:n]
    win_samp = int(round(window_sec * fs))
    step_samp = max(1, int(round(step_sec * fs)))
    starts = np.arange(0, n - win_samp + 1, step_samp)
    if len(starts) == 0:
        return None, None, None

    cols = []
    freq_ref = None
    for s0 in starts:
        seg1 = s1[s0:s0 + win_samp] - np.mean(s1[s0:s0 + win_samp])
        seg2 = s2[s0:s0 + win_samp] - np.mean(s2[s0:s0 + win_samp])
        seg1 = np.nan_to_num(seg1)
        seg2 = np.nan_to_num(seg2)
        inner_nperseg = max(16, win_samp // 4)
        inner_noverlap = inner_nperseg // 2
        inner_nfft = _next_pow2(inner_nperseg * 2)
        f, c = signal.coherence(seg1, seg2, fs=fs, nperseg=inner_nperseg,
                                noverlap=inner_noverlap, nfft=inner_nfft, window="hann")
        if freq_ref is None:
            freq_ref = f
        cols.append(c)

    centers = (starts + win_samp / 2.0) / fs
    mat = np.column_stack(cols)
    return freq_ref, centers, mat


def compute_1f_deviation(freq, psd_db, theta_range=THETA_BAND):
    """
    Oscillatory theta above 1/f background (consistent with fig5_theta.py).
    Fit linear trend in log10(freq) vs PSD(dB) from flanking bands; return peak − predicted.
    """
    freq = np.asarray(freq, dtype=float).ravel()
    psd_db = np.asarray(psd_db, dtype=float).ravel()
    theta_mask = (freq >= theta_range[0]) & (freq <= theta_range[1])
    if not np.any(theta_mask):
        return np.nan
    theta_peak_db = float(np.nanmax(psd_db[theta_mask]))
    theta_peak_freq = float(freq[theta_mask][np.nanargmax(psd_db[theta_mask])])

    low_mask = (freq >= SPEC5_1F_LOW_FLANK[0]) & (freq <= SPEC5_1F_LOW_FLANK[1])
    high_mask = (freq >= SPEC5_1F_HIGH_FLANK[0]) & (freq <= SPEC5_1F_HIGH_FLANK[1])
    fit_mask = low_mask | high_mask
    if np.sum(fit_mask) < 4:
        return np.nan
    try:
        log_f = np.log10(freq[fit_mask])
        coeffs = np.polyfit(log_f, psd_db[fit_mask], 1)
        expected = coeffs[0] * np.log10(theta_peak_freq) + coeffs[1]
        return theta_peak_db - expected
    except Exception:
        return np.nan


def compute_theta_band_value(freq, spectrum, theta_range=THETA_BAND,
                             method=SPEC5_THETA_EXTRACTION_METHOD):
    """Extract a scalar theta-band value: peak, mean, or AUC."""
    freq = np.asarray(freq, dtype=float).ravel()
    spectrum = np.asarray(spectrum, dtype=float).ravel()
    mask = (freq >= theta_range[0]) & (freq <= theta_range[1])
    if not np.any(mask):
        return np.nan
    vals = spectrum[mask]
    if method == "peak":
        return float(np.nanmax(vals))
    elif method == "auc":
        df = float(np.mean(np.diff(freq[mask]))) if np.sum(mask) > 1 else 1.0
        return float(np.nansum(vals) * df)
    return float(np.nanmean(vals))


# ---- per-trial collection across sessions ----

def collect_trialwise_psd_and_coherence(animal, session=None, band=THETA_BAND):
    """
    Per-trial Welch PSD (both fibers) and fiber-fiber coherence, REST/RUN split.
    Returns dict with freq, per-trial arrays, grand means, SEMs, and per-trial
    theta values (1/f deviation for PSD, peak for coherence) for violins.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    if session is not None:
        if session not in sessions_map:
            return None
        session_items = [(session, sessions_map[session])]
    else:
        session_items = list(sessions_map.items())

    trials_out = []
    freq_psd_ref = freq_coh_ref = None
    failed = 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                failed += 1
                continue
            f1 = np.asarray(data["fiber1"], dtype=float).ravel()
            f2 = np.asarray(data["fiber2"], dtype=float).ravel()
            if f2.size == 0 or not np.any(np.isfinite(f2)):
                failed += 1
                continue
            fs = float(data["fs"])
            speed = data["speed"]
            if speed is None:
                failed += 1
                continue

            f1_pct = (f1 - np.nanmean(f1)) * 100.0
            f2_pct = (f2 - np.nanmean(f2)) * 100.0

            res_p1 = compute_psd_rest_run(f1_pct, speed, fs, data["time"])
            res_p2 = compute_psd_rest_run(f2_pct, speed, fs, data["time"])
            res_c = compute_coherence_rest_run(f1_pct, f2_pct, speed, fs, data["time"])

            if res_p1[0] is None and res_p2[0] is None and res_c[0] is None:
                failed += 1
                continue

            entry = {"session": sess_key, "trial_num": trial_num}

            if res_p1[0] is not None:
                freq_psd = res_p1[0]
                if freq_psd_ref is None:
                    freq_psd_ref = freq_psd
                entry["psd_r_rest"] = res_p1[1]
                entry["psd_r_run"] = res_p1[2]
                entry["theta_1f_r_rest"] = (
                    compute_1f_deviation(freq_psd, res_p1[1], band)
                    if res_p1[1] is not None else np.nan
                )
                entry["theta_1f_r_run"] = (
                    compute_1f_deviation(freq_psd, res_p1[2], band)
                    if res_p1[2] is not None else np.nan
                )

            if res_p2[0] is not None:
                freq_psd = res_p2[0]
                if freq_psd_ref is None:
                    freq_psd_ref = freq_psd
                entry["psd_l_rest"] = res_p2[1]
                entry["psd_l_run"] = res_p2[2]
                entry["theta_1f_l_rest"] = (
                    compute_1f_deviation(freq_psd, res_p2[1], band)
                    if res_p2[1] is not None else np.nan
                )
                entry["theta_1f_l_run"] = (
                    compute_1f_deviation(freq_psd, res_p2[2], band)
                    if res_p2[2] is not None else np.nan
                )

            if res_c[0] is not None:
                freq_coh = res_c[0]
                if freq_coh_ref is None:
                    freq_coh_ref = freq_coh
                entry["coh_rest"] = res_c[1]
                entry["coh_run"] = res_c[2]
                entry["theta_coh_rest"] = (
                    compute_theta_band_value(freq_coh, res_c[1], band)
                    if res_c[1] is not None else np.nan
                )
                entry["theta_coh_run"] = (
                    compute_theta_band_value(freq_coh, res_c[2], band)
                    if res_c[2] is not None else np.nan
                )

            trials_out.append(entry)

    if not trials_out:
        return None

    def _stack_and_stats(trials, key):
        rows = [t[key] for t in trials if key in t and t[key] is not None]
        if not rows:
            return None, None, None
        mat = np.vstack(rows)
        m = np.nanmean(mat, axis=0)
        sem = np.nanstd(mat, axis=0, ddof=1) / np.sqrt(len(rows)) if len(rows) > 1 else None
        return mat, m, sem

    result = {
        "trials": trials_out,
        "freq_psd": freq_psd_ref,
        "freq_coh": freq_coh_ref,
        "n_loaded": len(trials_out),
        "n_failed": failed,
    }

    for suffix in ("psd_r_rest", "psd_r_run", "psd_l_rest", "psd_l_run",
                    "coh_rest", "coh_run"):
        _, m, sem = _stack_and_stats(trials_out, suffix)
        result[f"mean_{suffix}"] = m
        result[f"sem_{suffix}"] = sem

    print(
        f"  Spectral all trials: {len(trials_out)} loaded, {failed} failed"
    )
    return result


# ---- plotting helpers for Figure 5 ----

def _plot_spec5_psd_panel(ax, freq, psd_rest, psd_run, color, title, show_legend=True):
    """PSD REST vs RUN on one axis (single-trial or mean+SEM)."""
    fmask = None
    if freq is not None:
        fmask = (freq >= SPEC5_PSD_FREQ_RANGE[0]) & (freq <= SPEC5_PSD_FREQ_RANGE[1])
    has_data = False
    if psd_rest is not None and freq is not None:
        ax.plot(freq[fmask], psd_rest[fmask], "-", color=color, lw=SPEC5_PSD_LINEWIDTH,
                alpha=0.9, label="REST")
        has_data = True
    if psd_run is not None and freq is not None:
        ax.plot(freq[fmask], psd_run[fmask], "--", color=color, lw=SPEC5_PSD_LINEWIDTH,
                alpha=0.9, label="RUN")
        has_data = True
    if not has_data:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                va="center", fontsize=FONT_SIZE_LABEL, color="gray")
    ax.set_xlabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel("PSD (dB)", fontsize=FONT_SIZE_LABEL)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.axvspan(THETA_BAND[0], THETA_BAND[1], alpha=0.08, color="grey")
    style_axis(ax)
    if show_legend and has_data:
        ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, frameon=True,
                  framealpha=0.95, edgecolor="0.5")


def _plot_spec5_psd_mean_sem(ax, freq, mean_rest, sem_rest, mean_run, sem_run,
                             color, title, n_rest=0, n_run=0):
    """Mean PSD ± SEM, REST vs RUN."""
    fmask = None
    if freq is not None:
        fmask = (freq >= SPEC5_PSD_FREQ_RANGE[0]) & (freq <= SPEC5_PSD_FREQ_RANGE[1])
    has_data = False
    if mean_rest is not None and freq is not None:
        f = freq[fmask]
        m = mean_rest[fmask]
        ax.plot(f, m, "-", color=color, lw=SPEC5_PSD_LINEWIDTH, label=f"REST (n={n_rest})")
        if sem_rest is not None:
            s = sem_rest[fmask]
            ax.fill_between(f, m - s, m + s, color=color, alpha=SPEC5_SEM_ALPHA)
        has_data = True
    if mean_run is not None and freq is not None:
        f = freq[fmask]
        m = mean_run[fmask]
        ax.plot(f, m, "--", color=color, lw=SPEC5_PSD_LINEWIDTH, label=f"RUN (n={n_run})")
        if sem_run is not None:
            s = sem_run[fmask]
            ax.fill_between(f, m - s, m + s, color=color, alpha=SPEC5_SEM_ALPHA)
        has_data = True
    if not has_data:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                va="center", fontsize=FONT_SIZE_LABEL, color="gray")
    ax.set_xlabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel("PSD (dB)", fontsize=FONT_SIZE_LABEL)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.axvspan(THETA_BAND[0], THETA_BAND[1], alpha=0.08, color="grey")
    style_axis(ax)
    if has_data:
        ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, frameon=True,
                  framealpha=0.95, edgecolor="0.5")


def _plot_spec5_coh_panel(ax, freq, coh_rest, coh_run, title, show_legend=True):
    """Coherence REST vs RUN, dusty coral colour scheme."""
    color = tuple(np.clip(COLOR_SYNCHRONY_REST, 0, 1))
    color_run = tuple(np.clip(COLOR_SYNCHRONY_RUN, 0, 1))
    fmask = None
    if freq is not None:
        fmask = (freq >= SPEC5_PSD_FREQ_RANGE[0]) & (freq <= SPEC5_PSD_FREQ_RANGE[1])
    has_data = False
    if coh_rest is not None and freq is not None:
        ax.plot(freq[fmask], coh_rest[fmask], "-", color=color, lw=SPEC5_PSD_LINEWIDTH,
                alpha=0.9, label="REST")
        has_data = True
    if coh_run is not None and freq is not None:
        ax.plot(freq[fmask], coh_run[fmask], "--", color=color_run, lw=SPEC5_PSD_LINEWIDTH,
                alpha=0.9, label="RUN")
        has_data = True
    if not has_data:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                va="center", fontsize=FONT_SIZE_LABEL, color="gray")
    ax.set_xlabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel("Coherence", fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.axvspan(THETA_BAND[0], THETA_BAND[1], alpha=0.08, color="grey")
    style_axis(ax)
    if show_legend and has_data:
        ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, frameon=True,
                  framealpha=0.95, edgecolor="0.5")


def _plot_spec5_coh_mean_sem(ax, freq, mean_rest, sem_rest, mean_run, sem_run,
                             title, n_rest=0, n_run=0):
    """Mean coherence ± SEM."""
    color = tuple(np.clip(COLOR_SYNCHRONY_REST, 0, 1))
    color_run = tuple(np.clip(COLOR_SYNCHRONY_RUN, 0, 1))
    fmask = None
    if freq is not None:
        fmask = (freq >= SPEC5_PSD_FREQ_RANGE[0]) & (freq <= SPEC5_PSD_FREQ_RANGE[1])
    has_data = False
    if mean_rest is not None and freq is not None:
        f = freq[fmask]
        m = mean_rest[fmask]
        ax.plot(f, m, "-", color=color, lw=SPEC5_PSD_LINEWIDTH, label=f"REST (n={n_rest})")
        if sem_rest is not None:
            s = sem_rest[fmask]
            ax.fill_between(f, m - s, m + s, color=color, alpha=SPEC5_SEM_ALPHA)
        has_data = True
    if mean_run is not None and freq is not None:
        f = freq[fmask]
        m = mean_run[fmask]
        ax.plot(f, m, "--", color=color_run, lw=SPEC5_PSD_LINEWIDTH, label=f"RUN (n={n_run})")
        if sem_run is not None:
            s = sem_run[fmask]
            ax.fill_between(f, m - s, m + s, color=color_run, alpha=SPEC5_SEM_ALPHA)
        has_data = True
    if not has_data:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center",
                va="center", fontsize=FONT_SIZE_LABEL, color="gray")
    ax.set_xlabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel("Coherence", fontsize=FONT_SIZE_LABEL)
    ax.set_ylim(0, 1)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.axvspan(THETA_BAND[0], THETA_BAND[1], alpha=0.08, color="grey")
    style_axis(ax)
    if has_data:
        ax.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, frameon=True,
                  framealpha=0.95, edgecolor="0.5")


def _plot_spec5_theta_violin(ax, trials, key_rest, key_run, ylabel, title, colors):
    """Half-violin + box + scatter + stats for per-trial theta metric (reuses existing geometry)."""
    rest_vals, run_vals = [], []
    for t in trials:
        rr = t.get(key_rest, np.nan)
        rn = t.get(key_run, np.nan)
        if np.isfinite(rr) and np.isfinite(rn):
            rest_vals.append(rr)
            run_vals.append(rn)
    r_rest = np.asarray(rest_vals, dtype=float)
    r_run = np.asarray(run_vals, dtype=float)
    positions = (0.7, 1.3)
    labels = ("REST", "RUN")

    if r_rest.size == 0:
        ax.text(0.5, 0.5, "Insufficient paired trials", transform=ax.transAxes,
                ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
        style_axis(ax)
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
        ax.set_xlim(0.15, 1.85)
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
        return

    _phase_locking_plot_half_violin_R(ax, r_rest, r_run, positions, colors)

    all_y = np.concatenate([r_rest, r_run])
    y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
    y_range = y_max - y_min
    if y_range < 1e-10:
        y_range = max(abs(y_max), 0.05) * 0.3

    test_result = _perform_paired_test(r_rest, r_run)
    p_val = test_result.get("p_value", np.nan)
    test_name = test_result.get("test_used", "n/a")
    cohens_d = test_result.get("effect_size", np.nan)

    bracket_y = y_max + 0.06 * y_range
    line_height = 0.015 * y_range
    text_offset = 0.005 * y_range
    if np.isfinite(p_val):
        _add_significance_bracket(ax, positions[0], positions[1], bracket_y, p_val,
                                  line_height=line_height, text_offset=text_offset)
        top_pad = bracket_y + 0.12 * y_range
    else:
        top_pad = y_max + 0.15 * y_range

    ax.set_ylim(y_min - y_range * 0.40, top_pad)
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
    ax.set_title(title, fontsize=FONT_SIZE_TITLE - 1, fontweight="normal", pad=8)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
    ax.set_xlim(0.15, 1.85)
    style_axis(ax)

    p_str = (
        f"p = {p_val:.4f}" if np.isfinite(p_val) and p_val >= 0.001
        else f"p = {p_val:.2e}" if np.isfinite(p_val)
        else "p = n/a"
    )
    d_str = f"d = {cohens_d:.2f}" if np.isfinite(cohens_d) else "d = n/a"
    stats_txt = (
        f"REST = {float(np.mean(r_rest)):.3f}\n"
        f"RUN  = {float(np.mean(r_run)):.3f}\n"
        f"N = {len(r_rest)}, {test_name}\n"
        f"{p_str}, Cohen's {d_str}"
    )
    ax.text(0.02, 0.98, stats_txt, transform=ax.transAxes, fontsize=FONT_SIZE_TICK - 2,
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="0.55", alpha=0.92))


# ---- master figure ----

def fig_spectral_bilateral(data, animal, session, trial_num, band=THETA_BAND):
    """
    Figure 5 — Bilateral fiber spectral analysis.
    Row 1: Coherogram (speed bar + per-trial heatmaps) for example session.
    Row 2: Right HP PSD (single-trial | all-trials mean±SEM | theta violin).
    Row 3: Left HP PSD (same layout).
    Row 4: Fiber-fiber coherence (same layout).
    """
    # ---- example-session trial data for Row 1 ----
    sess_info = RECORDINGS[animal][session]
    n_trials_sess = int(sess_info["n_trials"])
    session_trials = []
    for tn in range(1, n_trials_sess + 1):
        try:
            d = load_trial(animal, session, tn)
            session_trials.append((tn, d))
        except Exception:
            pass
    n_panels = len(session_trials)
    if n_panels == 0:
        n_panels = 1

    # ---- all-trials pooled data for Rows 2-4 ----
    pool_session = None if SPEC5_AGGREGATE_ALL_SESSIONS else session
    print("  Collecting bilateral spectral data across all trials ...")
    spec_data = collect_trialwise_psd_and_coherence(animal, session=pool_session, band=band)

    # ---- single-trial PSD/coherence for Rows 2-4 (from the provided `data`) ----
    fs = float(data["fs"])
    f1_pct = (np.asarray(data["fiber1"], dtype=float) - np.nanmean(data["fiber1"])) * 100.0
    f2_pct = (np.asarray(data["fiber2"], dtype=float) - np.nanmean(data["fiber2"])) * 100.0
    speed = data["speed"]
    st_psd_r = compute_psd_rest_run(f1_pct, speed, fs, data["time"])
    st_psd_l = compute_psd_rest_run(f2_pct, speed, fs, data["time"])
    st_coh = compute_coherence_rest_run(f1_pct, f2_pct, speed, fs, data["time"])

    # ---- figure layout ----
    fig = plt.figure(figsize=FIGSIZE_SPEC5_INCH)

    # Row 1: 2 sub-rows (speed bar + coherogram) per trial, laid out horizontally
    # Rows 2-4: 3 columns each (single trial | all trials | violin)
    n_coh_cols = max(n_panels, 1) + 1  # +1 for colorbar
    gs_outer = GridSpec(4, 1, figure=fig,
                        height_ratios=[1.2, 1.0, 1.0, 1.0],
                        hspace=0.35, left=0.05, right=0.97, top=0.94, bottom=0.05)

    # Row 1 sub-grid: speed bar row + coherogram row × n_panels + colorbar
    coh_w = [1.0] * n_panels + [0.04]
    gs_row1 = gs_outer[0].subgridspec(2, len(coh_w),
                                       height_ratios=[0.12, 1.0],
                                       width_ratios=coh_w, wspace=0.08, hspace=0.08)

    # Rows 2-4 sub-grids: 3 columns each
    gs_row2 = gs_outer[1].subgridspec(1, 3, width_ratios=[1, 1, 0.82], wspace=0.28)
    gs_row3 = gs_outer[2].subgridspec(1, 3, width_ratios=[1, 1, 0.82], wspace=0.28)
    gs_row4 = gs_outer[3].subgridspec(1, 3, width_ratios=[1, 1, 0.82], wspace=0.28)

    # ======== ROW 1: COHEROGRAMS ========
    coh_im_ref = None
    for i, (tn, d) in enumerate(session_trials):
        d_f1 = (np.asarray(d["fiber1"], dtype=float) - np.nanmean(d["fiber1"])) * 100.0
        d_f2 = (np.asarray(d["fiber2"], dtype=float) - np.nanmean(d["fiber2"])) * 100.0
        d_fs = float(d["fs"])
        d_time = d["time"]

        # Speed bar
        ax_sp = fig.add_subplot(gs_row1[0, i])
        if d["speed"] is not None:
            sp2d = np.asarray(d["speed"]).reshape(1, -1)
            vmax_sp = max(float(np.nanpercentile(d["speed"], 99)), 1.0)
            ax_sp.imshow(sp2d, aspect="auto", cmap=CMAP_SPEED,
                         extent=[d_time[0], d_time[-1], 0, 1], origin="lower",
                         interpolation="bilinear", vmin=0, vmax=vmax_sp)
        ax_sp.set_yticks([])
        ax_sp.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK - 2)
        for sp in ("top", "right", "left"):
            ax_sp.spines[sp].set_visible(False)
        ax_sp.spines["bottom"].set_linewidth(AXIS_LINEWIDTH * 0.5)
        ax_sp.set_title(f"Trial {tn}", fontsize=FONT_SIZE_TITLE - 2, fontweight="bold", pad=4)

        # Coherogram
        ax_cg = fig.add_subplot(gs_row1[1, i])
        freq_cg, t_cg, mat_cg = compute_coherogram(d_f1, d_f2, d_fs)
        if freq_cg is not None:
            fmask_cg = (freq_cg >= SPEC5_PSD_FREQ_RANGE[0]) & (freq_cg <= SPEC5_PSD_FREQ_RANGE[1])
            im = ax_cg.pcolormesh(t_cg, freq_cg[fmask_cg], mat_cg[fmask_cg, :],
                                  shading="gouraud", cmap=SPEC5_COHEROGRAM_CMAP,
                                  vmin=0, vmax=1.0, rasterized=True)
            coh_im_ref = im
        ax_cg.set_xlabel("Time (s)", fontsize=FONT_SIZE_LABEL - 2)
        if i == 0:
            ax_cg.set_ylabel("Frequency (Hz)", fontsize=FONT_SIZE_LABEL - 2)
        else:
            ax_cg.tick_params(labelleft=False)
        ax_cg.tick_params(labelsize=FONT_SIZE_TICK - 2)
        style_axis(ax_cg)

    # Colorbar for coherograms
    cax_cg = fig.add_subplot(gs_row1[1, n_panels])
    if coh_im_ref is not None:
        cb = fig.colorbar(coh_im_ref, cax=cax_cg, orientation="vertical")
        cb.set_label("Coherence", fontsize=FONT_SIZE_COLORBAR - 2, fontweight="bold")
        cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 2)
    else:
        cax_cg.set_visible(False)
    # Hide speed-bar colorbar slot
    ax_sp_cb = fig.add_subplot(gs_row1[0, n_panels])
    ax_sp_cb.set_visible(False)

    # ======== ROW 2: RIGHT HP PSD ========
    c_r = tuple(np.clip(COLOR_RIGHT_HP, 0, 1))
    ax_psd_r_st = fig.add_subplot(gs_row2[0])
    _plot_spec5_psd_panel(ax_psd_r_st, st_psd_r[0], st_psd_r[1], st_psd_r[2],
                          c_r, f"Right HP PSD — example trial")

    ax_psd_r_all = fig.add_subplot(gs_row2[1])
    if spec_data is not None:
        n_r_rest = sum(1 for t in spec_data["trials"] if t.get("psd_r_rest") is not None)
        n_r_run = sum(1 for t in spec_data["trials"] if t.get("psd_r_run") is not None)
        _plot_spec5_psd_mean_sem(
            ax_psd_r_all, spec_data["freq_psd"],
            spec_data["mean_psd_r_rest"], spec_data["sem_psd_r_rest"],
            spec_data["mean_psd_r_run"], spec_data["sem_psd_r_run"],
            c_r, "Right HP PSD — all trials", n_r_rest, n_r_run,
        )
    else:
        ax_psd_r_all.text(0.5, 0.5, "No data", transform=ax_psd_r_all.transAxes,
                          ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        style_axis(ax_psd_r_all)

    ax_psd_r_v = fig.add_subplot(gs_row2[2])
    if spec_data is not None:
        _plot_spec5_theta_violin(
            ax_psd_r_v, spec_data["trials"],
            "theta_1f_r_rest", "theta_1f_r_run",
            r"$\theta$ prominence (dB)", "Right HP theta",
            (COLOR_RIGHT_HP, COLOR_RIGHT_HP * 0.65),
        )
    else:
        style_axis(ax_psd_r_v)

    # ======== ROW 3: LEFT HP PSD ========
    c_l = tuple(np.clip(COLOR_LEFT_HP, 0, 1))
    ax_psd_l_st = fig.add_subplot(gs_row3[0])
    _plot_spec5_psd_panel(ax_psd_l_st, st_psd_l[0], st_psd_l[1], st_psd_l[2],
                          c_l, f"Left HP PSD — example trial")

    ax_psd_l_all = fig.add_subplot(gs_row3[1])
    if spec_data is not None:
        n_l_rest = sum(1 for t in spec_data["trials"] if t.get("psd_l_rest") is not None)
        n_l_run = sum(1 for t in spec_data["trials"] if t.get("psd_l_run") is not None)
        _plot_spec5_psd_mean_sem(
            ax_psd_l_all, spec_data["freq_psd"],
            spec_data["mean_psd_l_rest"], spec_data["sem_psd_l_rest"],
            spec_data["mean_psd_l_run"], spec_data["sem_psd_l_run"],
            c_l, "Left HP PSD — all trials", n_l_rest, n_l_run,
        )
    else:
        ax_psd_l_all.text(0.5, 0.5, "No data", transform=ax_psd_l_all.transAxes,
                          ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        style_axis(ax_psd_l_all)

    ax_psd_l_v = fig.add_subplot(gs_row3[2])
    if spec_data is not None:
        _plot_spec5_theta_violin(
            ax_psd_l_v, spec_data["trials"],
            "theta_1f_l_rest", "theta_1f_l_run",
            r"$\theta$ prominence (dB)", "Left HP theta",
            (COLOR_LEFT_HP, COLOR_LEFT_HP * 0.65),
        )
    else:
        style_axis(ax_psd_l_v)

    # ======== ROW 4: FIBER-FIBER COHERENCE ========
    ax_coh_st = fig.add_subplot(gs_row4[0])
    _plot_spec5_coh_panel(ax_coh_st, st_coh[0], st_coh[1], st_coh[2],
                          f"Fiber-fiber coherence — example trial")

    ax_coh_all = fig.add_subplot(gs_row4[1])
    if spec_data is not None:
        n_c_rest = sum(1 for t in spec_data["trials"] if t.get("coh_rest") is not None)
        n_c_run = sum(1 for t in spec_data["trials"] if t.get("coh_run") is not None)
        _plot_spec5_coh_mean_sem(
            ax_coh_all, spec_data["freq_coh"],
            spec_data["mean_coh_rest"], spec_data["sem_coh_rest"],
            spec_data["mean_coh_run"], spec_data["sem_coh_run"],
            "Fiber-fiber coherence — all trials", n_c_rest, n_c_run,
        )
    else:
        ax_coh_all.text(0.5, 0.5, "No data", transform=ax_coh_all.transAxes,
                        ha="center", va="center", fontsize=FONT_SIZE_LABEL, color="gray")
        style_axis(ax_coh_all)

    ax_coh_v = fig.add_subplot(gs_row4[2])
    if spec_data is not None:
        _plot_spec5_theta_violin(
            ax_coh_v, spec_data["trials"],
            "theta_coh_rest", "theta_coh_run",
            r"$\theta$ coherence", "Theta coherence",
            (COLOR_SYNCHRONY_REST, COLOR_SYNCHRONY_RUN),
        )
    else:
        style_axis(ax_coh_v)

    # ---- suptitle ----
    n_all = spec_data["n_loaded"] if spec_data else 0
    fig.suptitle(
        f"{animal} — Bilateral fiber spectral analysis "
        f"(coherogram: {session}; PSD/coherence: {n_all} trials pooled)",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.98,
    )
    return fig


# =============================================================================
# FIGURE 6b: FIBER–LFP THETA CROSS-CORRELATION – REST vs RUN (4 COMBOS)
# =============================================================================

def collect_trialwise_fiber_lfp_xcorr_by_behavior(animal, session=None,
                                                    band=THETA_BAND):
    """
    Per-trial REST/RUN split theta cross-correlation for the 4 fiber-LFP
    combinations.  Re-uses compute_xcorr_fiber_by_behavior() for the
    actual windowed correlogram computation.

    Returns dict with per-combo aggregated REST and RUN mean correlograms,
    peak |r| arrays, grand means / SEM, and per-trial peak dicts.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    session_items = ([(session, sessions_map[session])]
                     if session is not None else list(sessions_map.items()))

    combo_keys = [c["short"] for c in FIG6_COMBINATIONS]
    combo_data = {
        k: {
            "mean_rows_rest": [], "mean_rows_run": [],
            "peak_r_rest": [], "peak_r_run": [],
        }
        for k in combo_keys
    }
    trial_peaks = []
    lags_ms_ref = None
    n_loaded, n_failed = 0, 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue
            n_loaded += 1
            fs = float(data["fs"])
            tp = {"session": sess_key, "trial_num": trial_num}

            for combo in FIG6_COMBINATIONS:
                sig1 = np.asarray(data[combo["sig1_key"]], dtype=float).ravel()
                sig2 = data.get(combo["sig2_key"])
                if sig2 is None:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue
                sig2 = np.asarray(sig2, dtype=float).ravel()
                if not np.any(np.isfinite(sig1)) or not np.any(np.isfinite(sig2)):
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue

                try:
                    out = compute_xcorr_fiber_by_behavior(
                        sig1, sig2, fs, data["time"], data["speed"],
                        band=band,
                    )
                except Exception:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue

                if lags_ms_ref is None:
                    lags_ms_ref = out["lags_ms"]

                cd = combo_data[combo["short"]]

                for state in ("rest", "run"):
                    mat_key = f"mat_{state}"
                    mat = out[mat_key]
                    mean_corr = mean_xcorr_over_windows(mat)
                    if mean_corr is not None:
                        ir, _ = peak_lag_index_restricted(
                            mean_corr, out["lags_ms"], band)
                        peak_r = float(np.abs(mean_corr[ir]))
                        cd[f"mean_rows_{state}"].append(mean_corr)
                        cd[f"peak_r_{state}"].append(peak_r)
                        tp[combo["short"] + f"_{state}"] = peak_r
                    else:
                        tp[combo["short"] + f"_{state}"] = np.nan

            trial_peaks.append(tp)
            line_parts = []
            for c in FIG6_COMBINATIONS:
                vr = tp.get(c["short"] + "_rest", np.nan)
                vn = tp.get(c["short"] + "_run", np.nan)
                line_parts.append(
                    f"{c['short']} R={vr:.3f}/{vn:.3f}"
                )
            print(f"    {sess_key} T{trial_num}: {' | '.join(line_parts)}")

    if lags_ms_ref is None:
        return None

    def _agg(rows):
        if not rows:
            return None, None, 0, np.array([])
        mat = np.vstack(rows)
        n = len(rows)
        if XCORR_MEAN_USE_FISHER_Z:
            m = np.clip(mat, -0.999999, 0.999999)
            gm = np.tanh(np.nanmean(np.arctanh(m), axis=0))
        else:
            gm = np.nanmean(mat, axis=0)
        sem = (np.nanstd(mat, axis=0, ddof=1) / np.sqrt(n)
               if n > 1 else None)
        return gm, sem, n, mat

    for key, cd in combo_data.items():
        for state in ("rest", "run"):
            gm, sem, n, mat = _agg(cd[f"mean_rows_{state}"])
            cd[f"grand_mean_{state}"] = gm
            cd[f"sem_{state}"] = sem
            cd[f"n_trials_{state}"] = n
            cd[f"mat_all_{state}"] = mat
            cd[f"peak_r_{state}"] = np.array(cd[f"peak_r_{state}"])

    print(f"  Fiber-LFP xcorr REST/RUN: {n_loaded} loaded, {n_failed} failed")
    return {
        "combo_data": combo_data,
        "trial_peaks": trial_peaks,
        "lags_ms": lags_ms_ref,
        "n_loaded": n_loaded,
        "n_failed": n_failed,
    }


def fig_fiber_lfp_xcorr_rest_run(animal, band=THETA_BAND):
    """
    Figure 6b: Fiber–LFP theta cross-correlation REST vs RUN.

    Layout – 2 rows × (4 combo columns + 1 narrow spacer + 4 violin columns):
      top row:  all-trials grand-mean ± SEM correlogram per combo
                (REST solid, RUN dashed) — one panel per combo
      bottom row: half-violin/box/scatter REST vs RUN peak |r| per combo
    """
    sess_arg = None if FIG6B_AGGREGATE_ALL_SESSIONS else None
    print("  Collecting fiber–LFP xcorr by behaviour across all trials ...")
    all_data = collect_trialwise_fiber_lfp_xcorr_by_behavior(
        animal, session=sess_arg, band=band,
    )

    fig = plt.figure(figsize=FIGSIZE_FIG6B_INCH)
    gs = GridSpec(
        2, 4, figure=fig,
        height_ratios=[1.3, 1.0],
        wspace=0.30, hspace=0.35,
        left=0.05, right=0.97, top=0.88, bottom=0.08,
    )

    peak_lim = _peak_lag_search_half_width_ms(band)

    for i, combo in enumerate(FIG6_COMBINATIONS):
        color = tuple(np.clip(FIG6_COLORS[combo["short"]], 0, 1))
        color_dark = tuple(np.clip(FIG6_COLORS[combo["short"]] * 0.7, 0, 1))

        # --- top: correlogram ---
        ax_corr = fig.add_subplot(gs[0, i])
        ax_corr.axvline(0, color="grey", ls="--", lw=1)
        ax_corr.axhline(0, color="grey", ls="-", lw=0.5, alpha=0.5)
        ax_corr.axvline(-peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7,
                        zorder=0)
        ax_corr.axvline(peak_lim, color="0.65", ls=":", lw=1.0, alpha=0.7,
                        zorder=0)

        stats_lines = []
        if all_data is not None:
            cd = all_data["combo_data"][combo["short"]]
            lags = all_data["lags_ms"]
            for state, ls_style, lbl_st, col in [
                ("rest", "-", "REST", color),
                ("run", "--", "RUN", color_dark),
            ]:
                gm = cd[f"grand_mean_{state}"]
                sem = cd[f"sem_{state}"]
                n = cd[f"n_trials_{state}"]
                if gm is not None:
                    ax_corr.plot(lags, gm, ls_style, color=col,
                                lw=FIG6_CORR_LINE_WIDTH, label=lbl_st)
                    if sem is not None:
                        ax_corr.fill_between(
                            lags, gm - sem, gm + sem,
                            color=col, alpha=FIG6_SEM_ALPHA,
                        )
                    ir, _ = peak_lag_index_restricted(gm, lags, band)
                    stats_lines.append(
                        f"{lbl_st} |r|={abs(gm[ir]):.3f} (n={n})"
                    )

        if stats_lines:
            ax_corr.text(
                0.03, 0.97, "\n".join(stats_lines),
                transform=ax_corr.transAxes, fontsize=FONT_SIZE_TICK - 2,
                ha="left", va="top",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.6",
                          alpha=0.92),
            )

        ax_corr.set_xlim(-XCORR_MAX_LAG_MS, XCORR_MAX_LAG_MS)
        ax_corr.xaxis.set_major_locator(MultipleLocator(50))
        ax_corr.set_xlabel("Temporal offset (ms)", fontsize=FONT_SIZE_LABEL)
        ax_corr.set_ylabel("Cross-correlation" if i == 0 else "",
                           fontsize=FONT_SIZE_LABEL)
        ax_corr.set_title(combo["label"], fontsize=FONT_SIZE_TITLE - 2,
                          fontweight="bold", pad=8)
        if i == 0:
            ax_corr.legend(fontsize=FONT_SIZE_TICK - 2, loc="lower left",
                           framealpha=0.85)
        style_axis(ax_corr)

        # --- bottom: violin REST vs RUN ---
        ax_v = fig.add_subplot(gs[1, i])
        short = combo["short"]

        if all_data is not None:
            rest_vals, run_vals = [], []
            for tp in all_data["trial_peaks"]:
                vr = tp.get(short + "_rest", np.nan)
                vn = tp.get(short + "_run", np.nan)
                if np.isfinite(vr) and np.isfinite(vn):
                    rest_vals.append(vr)
                    run_vals.append(vn)
            r_rest_arr = np.array(rest_vals)
            r_run_arr = np.array(run_vals)
            n_pair = len(r_rest_arr)

            if n_pair >= 3:
                _phase_locking_plot_half_violin_R(
                    ax_v, r_rest_arr, r_run_arr,
                    positions=(0.7, 1.3),
                    colors=(color, color_dark),
                )
                all_y = np.concatenate([r_rest_arr, r_run_arr])
                y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
                y_range = max(y_max - y_min, abs(y_max) * 0.05)

                test_result = _perform_paired_test(r_rest_arr, r_run_arr)
                p_val = test_result.get("p_value", np.nan)
                cohens_d = test_result.get("effect_size", np.nan)
                test_name = test_result.get("test_used", "n/a")

                bracket_y = y_max + 0.08 * y_range
                if np.isfinite(p_val):
                    _add_significance_bracket(
                        ax_v, 0.7, 1.3, bracket_y, p_val,
                        line_height=0.015 * y_range,
                        text_offset=0.005 * y_range,
                    )
                    top_pad = bracket_y + 0.14 * y_range
                else:
                    top_pad = y_max + 0.15 * y_range

                ax_v.set_ylim(y_min - y_range * 0.35, top_pad)

                p_str = (f"p={p_val:.4f}" if np.isfinite(p_val) and p_val >= 0.001
                         else f"p={p_val:.2e}" if np.isfinite(p_val)
                         else "p=n/a")
                d_str = (f"d={cohens_d:.2f}" if np.isfinite(cohens_d)
                         else "d=n/a")
                ax_v.text(
                    0.03, 0.02,
                    f"n={n_pair}, {test_name}\n{p_str}, Cohen's {d_str}",
                    transform=ax_v.transAxes, fontsize=FONT_SIZE_TICK - 2,
                    va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec="0.55", alpha=0.92),
                )
            else:
                ax_v.text(0.5, 0.5, "Insufficient paired\nREST/RUN data",
                          transform=ax_v.transAxes, ha="center", va="center",
                          fontsize=FONT_SIZE_LABEL - 2, color="gray")
        else:
            ax_v.text(0.5, 0.5, "No data", transform=ax_v.transAxes,
                      ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                      color="gray")

        ax_v.set_ylabel("Peak |r|" if i == 0 else "",
                        fontsize=FONT_SIZE_LABEL - 1)
        ax_v.set_title(short + " REST vs RUN",
                       fontsize=FONT_SIZE_TITLE - 2, fontweight="bold", pad=6)
        ax_v.set_xticks([0.7, 1.3])
        ax_v.set_xticklabels(["REST", "RUN"], fontsize=FONT_SIZE_TICK - 1)
        ax_v.set_xlim(0.15, 1.85)
        style_axis(ax_v)

    n_total = all_data["n_loaded"] if all_data else 0
    scope = "all sessions" if FIG6B_AGGREGATE_ALL_SESSIONS else ""
    fig.suptitle(
        f"{animal} — Fiber–LFP $\\theta$ cross-correlation, "
        f"REST vs RUN ({scope}, n={n_total} trials)",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.96,
    )
    return fig


# =============================================================================
# FIGURE 7b: FIBER–LFP PHASE-LOCKING – REST vs RUN (4 COMBOS)
# =============================================================================

def collect_trialwise_fiber_lfp_phase_by_behavior(
        animal, session=None, band=THETA_BAND,
        edge_trim_sec=PHASE_EDGE_TRIM_SEC):
    """
    Per-trial REST/RUN split phase-locking R for the 4 fiber–LFP combinations.
    Phase difference Δφ = φ_fiber − φ_LFP is computed via Hilbert transform;
    samples are split by REST/RUN masks and R is computed per state per trial.
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    session_items = ([(session, sessions_map[session])]
                     if session is not None else list(sessions_map.items()))

    combo_keys = [c["short"] for c in FIG6_COMBINATIONS]
    combo_data = {
        k: {
            "phi_rest_trials": [], "phi_run_trials": [],
            "R_rest_trials": [], "R_run_trials": [],
            "mu_rest_trials": [], "mu_run_trials": [],
        }
        for k in combo_keys
    }
    trial_peaks = []
    n_loaded, n_failed = 0, 0

    for sess_key, info in session_items:
        n_tr = int(info["n_trials"])
        for trial_num in range(1, n_tr + 1):
            try:
                data = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue

            fs = float(data["fs"])
            n_loaded += 1
            tp = {"session": sess_key, "trial_num": trial_num}

            speed = data["speed"]
            if speed is None:
                is_rest = np.ones(1, dtype=bool)
                is_run = np.zeros(1, dtype=bool)
                is_excluded = np.zeros(1, dtype=bool)
            else:
                is_rest, is_run, is_excluded = classify_rest_run_masks(
                    np.asarray(speed, dtype=float).ravel(), fs)

            for combo in FIG6_COMBINATIONS:
                sig1 = data.get(combo["sig1_key"])
                sig2 = data.get(combo["sig2_key"])
                if sig1 is None or sig2 is None:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue

                sig1 = np.asarray(sig1, dtype=float).ravel()
                sig2 = np.asarray(sig2, dtype=float).ravel()
                n = min(len(sig1), len(sig2), len(is_rest))
                if n < 20:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue

                sig1, sig2 = sig1[:n], sig2[:n]
                m_rest = is_rest[:n]
                m_run = is_run[:n]
                m_excl = is_excluded[:n]

                trim = max(0, int(round(edge_trim_sec * fs)))
                if n <= 2 * trim + 8:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue
                sl = slice(trim, n - trim)
                sig1, sig2 = sig1[sl], sig2[sl]
                m_rest, m_run = m_rest[sl], m_run[sl]
                m_excl = m_excl[sl]

                try:
                    z1 = theta_analytic_from_raw(sig1, fs, band=band)
                    z2 = theta_analytic_from_raw(sig2, fs, band=band)
                except Exception:
                    tp[combo["short"] + "_rest"] = np.nan
                    tp[combo["short"] + "_run"] = np.nan
                    continue

                dphi = wrap_to_pi(np.angle(z1) - np.angle(z2))
                ok = np.isfinite(dphi) & (~m_excl)

                cd = combo_data[combo["short"]]

                for state, mask in [("rest", m_rest), ("run", m_run)]:
                    sel = ok & mask
                    phi_sel = dphi[sel]
                    if phi_sel.size >= 10:
                        mu, R = circular_mean_and_resultant_length(phi_sel)
                        cd[f"phi_{state}_trials"].append(phi_sel)
                        cd[f"R_{state}_trials"].append(R)
                        cd[f"mu_{state}_trials"].append(mu)
                        tp[combo["short"] + f"_{state}"] = R
                    else:
                        tp[combo["short"] + f"_{state}"] = np.nan

            trial_peaks.append(tp)
            line_parts = []
            for c in FIG6_COMBINATIONS:
                vr = tp.get(c["short"] + "_rest", np.nan)
                vn = tp.get(c["short"] + "_run", np.nan)
                line_parts.append(f"{c['short']} R={vr:.4f}/{vn:.4f}")
            print(f"    {sess_key} T{trial_num}: {' | '.join(line_parts)}")

    if n_loaded == 0:
        return None

    for key, cd in combo_data.items():
        for state in ("rest", "run"):
            cd[f"R_{state}_values"] = np.array(cd[f"R_{state}_trials"])
            cd[f"n_{state}_trials"] = len(cd[f"R_{state}_trials"])

    print(f"  Fiber-LFP phase REST/RUN: {n_loaded} loaded, {n_failed} failed")
    return {
        "combo_data": combo_data,
        "trial_peaks": trial_peaks,
        "n_loaded": n_loaded,
        "n_failed": n_failed,
    }


def fig_fiber_lfp_phase_rest_run(animal, band=THETA_BAND):
    """
    Figure 7b: Fiber–LFP phase-locking REST vs RUN for 4 combinations.

    Layout – 3 rows × 4 columns:
      row 0: REST rose plots (one per combo)
      row 1: RUN rose plots (one per combo)
      row 2: half-violin/box/scatter REST vs RUN R per combo
    """
    sess_arg = None if FIG7B_AGGREGATE_ALL_SESSIONS else None
    print("  Collecting fiber–LFP phase by behaviour across all trials ...")
    phase_data = collect_trialwise_fiber_lfp_phase_by_behavior(
        animal, session=sess_arg, band=band,
    )

    fig = plt.figure(figsize=FIGSIZE_FIG7B_INCH)
    gs = GridSpec(
        3, 4, figure=fig,
        height_ratios=[1.2, 1.2, 1.0],
        wspace=0.30, hspace=0.40,
        left=0.05, right=0.97, top=0.90, bottom=0.06,
    )

    for i, combo in enumerate(FIG6_COMBINATIONS):
        color = tuple(np.clip(FIG6_COLORS[combo["short"]], 0, 1))
        color_dark = tuple(np.clip(FIG6_COLORS[combo["short"]] * 0.7, 0, 1))

        for row_idx, (state, state_label, state_color) in enumerate([
            ("rest", "REST", color),
            ("run", "RUN", color_dark),
        ]):
            ax = fig.add_subplot(gs[row_idx, i], projection="polar")

            if phase_data is not None:
                cd = phase_data["combo_data"][combo["short"]]
                phi_list = cd[f"phi_{state}_trials"]
            else:
                phi_list = []

            rose = _trial_equalized_rose_histogram_and_vector(
                phi_list, PHASE_ROSE_NBINS,
            )

            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)

            if rose is not None:
                ax.bar(
                    rose["theta"], rose["bar_heights"],
                    width=rose["width"] * 0.98, bottom=0.0, align="center",
                    color=state_color, edgecolor="white",
                    linewidth=PHASE_ROSE_BAR_EDGELW,
                    alpha=FIG7_ROSE_BAR_ALPHA,
                )
                bar_max = (float(np.max(rose["bar_heights"]))
                           if rose["bar_heights"].size else 10.0)
                if np.isfinite(rose["mu"]) and np.isfinite(rose["R"]):
                    ax.annotate(
                        "", xytext=(0.0, 0.0),
                        xy=(rose["mu"], rose["R"] * bar_max),
                        arrowprops=dict(
                            arrowstyle="-|>", color="crimson",
                            lw=PHASE_MEAN_VECTOR_LW,
                            mutation_scale=PHASE_MEAN_VECTOR_ARROW_MUTATION,
                            shrinkA=0, shrinkB=0,
                        ),
                        zorder=6,
                    )
                ax.set_ylim(0.0, bar_max * 1.12)

                deg = (np.degrees(rose["mu"]) if np.isfinite(rose["mu"])
                       else float("nan"))
                n_total = sum(len(p) for p in phi_list)
                ttl = (
                    f"{combo['short']} {state_label}\n"
                    f"N={n_total}, {rose['n_trials']} trials\n"
                    f"Mean={deg:.1f}$^\\circ$  R={rose['R']:.3f}"
                )
            else:
                ttl = f"{combo['short']} {state_label}\n(no data)"

            _style_phase_rose_polar_axis(ax)
            ax.set_title(
                ttl, fontsize=PHASE_ROSE_TITLE_FONTSIZE - 2,
                fontweight=PHASE_ROSE_TITLE_FONTWEIGHT,
                pad=PHASE_ROSE_TITLE_PAD,
            )

        # --- bottom row: violin REST vs RUN ---
        ax_v = fig.add_subplot(gs[2, i])
        short = combo["short"]

        if phase_data is not None:
            rest_vals, run_vals = [], []
            for tp in phase_data["trial_peaks"]:
                vr = tp.get(short + "_rest", np.nan)
                vn = tp.get(short + "_run", np.nan)
                if np.isfinite(vr) and np.isfinite(vn):
                    rest_vals.append(vr)
                    run_vals.append(vn)
            r_rest_arr = np.array(rest_vals)
            r_run_arr = np.array(run_vals)
            n_pair = len(r_rest_arr)

            if n_pair >= 3:
                _phase_locking_plot_half_violin_R(
                    ax_v, r_rest_arr, r_run_arr,
                    positions=(0.7, 1.3),
                    colors=(color, color_dark),
                )
                all_y = np.concatenate([r_rest_arr, r_run_arr])
                y_min, y_max = float(np.min(all_y)), float(np.max(all_y))
                y_range = max(y_max - y_min, abs(y_max) * 0.05)

                test_result = _perform_paired_test(r_rest_arr, r_run_arr)
                p_val = test_result.get("p_value", np.nan)
                cohens_d = test_result.get("effect_size", np.nan)
                test_name = test_result.get("test_used", "n/a")

                bracket_y = y_max + 0.08 * y_range
                if np.isfinite(p_val):
                    _add_significance_bracket(
                        ax_v, 0.7, 1.3, bracket_y, p_val,
                        line_height=0.015 * y_range,
                        text_offset=0.005 * y_range,
                    )
                    top_pad = bracket_y + 0.14 * y_range
                else:
                    top_pad = y_max + 0.15 * y_range
                ax_v.set_ylim(y_min - y_range * 0.35, top_pad)

                p_str = (f"p={p_val:.4f}"
                         if np.isfinite(p_val) and p_val >= 0.001
                         else f"p={p_val:.2e}" if np.isfinite(p_val)
                         else "p=n/a")
                d_str = (f"d={cohens_d:.2f}" if np.isfinite(cohens_d)
                         else "d=n/a")
                ax_v.text(
                    0.03, 0.02,
                    f"n={n_pair}, {test_name}\n{p_str}, Cohen's {d_str}",
                    transform=ax_v.transAxes, fontsize=FONT_SIZE_TICK - 2,
                    va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec="0.55", alpha=0.92),
                )
            else:
                ax_v.text(0.5, 0.5, "Insufficient paired\nREST/RUN trials",
                          transform=ax_v.transAxes, ha="center", va="center",
                          fontsize=FONT_SIZE_LABEL - 2, color="gray")
        else:
            ax_v.text(0.5, 0.5, "No data", transform=ax_v.transAxes,
                      ha="center", va="center", fontsize=FONT_SIZE_LABEL,
                      color="gray")

        ax_v.set_ylabel("Phase-locking R" if i == 0 else "",
                        fontsize=FONT_SIZE_LABEL - 1)
        ax_v.set_title(combo["short"] + " REST vs RUN",
                       fontsize=FONT_SIZE_TITLE - 2, fontweight="bold", pad=6)
        ax_v.set_xticks([0.7, 1.3])
        ax_v.set_xticklabels(["REST", "RUN"], fontsize=FONT_SIZE_TICK - 1)
        ax_v.set_xlim(0.15, 1.85)
        style_axis(ax_v)

    n_total = phase_data["n_loaded"] if phase_data else 0
    scope = "all sessions" if FIG7B_AGGREGATE_ALL_SESSIONS else ""
    fig.suptitle(
        f"{animal} — Fiber–LFP $\\theta$ phase-locking, "
        f"REST vs RUN ({scope}, n={n_total} trials)",
        fontsize=FONT_SIZE_SUPTITLE, fontweight="bold", y=0.97,
    )
    return fig


RUN_ONLY_NEW_FIGURES = True  # set False to generate all figures (1-5 included)

def main():
    print("=" * 70)
    print("  MULTI-SITE FIBER ANALYSIS")
    print("=" * 70)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ----- Figure 6: Fiber–LFP theta cross-correlation (4 combos) -----
    animal6 = FIG6_ANIMAL
    session6 = FIG6_EXAMPLE_SESSION
    trial6 = FIG6_EXAMPLE_TRIAL
    print(f"\nFigure 6: Fiber\u2013LFP theta cross-correlation "
          f"({animal6}, example {session6} T{trial6})...")
    prefix6 = (
        f"{animal6}_allSessions_example_{session6}_Trial{trial6}"
        if FIG6_AGGREGATE_ALL_SESSIONS
        else f"{animal6}_{session6}_Trial{trial6}"
    )
    f6 = fig_fiber_lfp_xcorr(animal6, session6, trial6)
    _save_phase_locking_figure(
        f6, OUTPUT_DIR / f"{prefix6}_06_fiber_lfp_xcorr"
    )
    plt.close(f6)
    print("  Saved.")

    # ----- Figure 7: Fiber–LFP phase-locking (4 combos, all epochs) -----
    animal7 = FIG7_ANIMAL
    print(f"\nFigure 7: Fiber\u2013LFP theta phase-locking ({animal7})...")
    prefix7 = (
        f"{animal7}_allSessions"
        if FIG7_AGGREGATE_ALL_SESSIONS
        else animal7
    )
    f7 = fig_fiber_lfp_phase_locking(animal7)
    _save_phase_locking_figure(
        f7, OUTPUT_DIR / f"{prefix7}_07_fiber_lfp_phase_locking"
    )
    plt.close(f7)
    print("  Saved.")

    # ----- Figure 6b: Fiber–LFP xcorr REST vs RUN -----
    animal6b = FIG6B_ANIMAL
    print(f"\nFigure 6b: Fiber\u2013LFP theta xcorr REST vs RUN ({animal6b})...")
    prefix6b = (
        f"{animal6b}_allSessions"
        if FIG6B_AGGREGATE_ALL_SESSIONS
        else animal6b
    )
    f6b = fig_fiber_lfp_xcorr_rest_run(animal6b)
    _save_phase_locking_figure(
        f6b, OUTPUT_DIR / f"{prefix6b}_06b_fiber_lfp_xcorr_rest_vs_run"
    )
    plt.close(f6b)
    print("  Saved.")

    # ----- Figure 7b: Fiber–LFP phase REST vs RUN -----
    animal7b = FIG7B_ANIMAL
    print(f"\nFigure 7b: Fiber\u2013LFP theta phase-locking REST vs RUN ({animal7b})...")
    prefix7b = (
        f"{animal7b}_allSessions"
        if FIG7B_AGGREGATE_ALL_SESSIONS
        else animal7b
    )
    f7b = fig_fiber_lfp_phase_rest_run(animal7b)
    _save_phase_locking_figure(
        f7b, OUTPUT_DIR / f"{prefix7b}_07b_fiber_lfp_phase_rest_vs_run"
    )
    plt.close(f7b)
    print("  Saved.")

    # ----- LFP–LFP coherence diagnostic -----
    animal_diag = FIG_LFP_DIAG_ANIMAL
    print(f"\nLFP\u2013LFP coherence diagnostic ({animal_diag})...")
    prefix_diag = (
        f"{animal_diag}_allSessions"
        if FIG_LFP_DIAG_AGGREGATE_ALL_SESSIONS
        else animal_diag
    )
    f_diag = fig_lfp_lfp_coherence_diagnostic(animal_diag)
    _save_phase_locking_figure(
        f_diag, OUTPUT_DIR / f"{prefix_diag}_LFP_LFP_coherence_diagnostic"
    )
    plt.close(f_diag)
    print("  Saved.")

    if RUN_ONLY_NEW_FIGURES:
        print(f"\nDONE! Figures 6, 6b, 7, 7b, LFP diagnostic saved to "
              f"{OUTPUT_DIR}")
        return

    animal = EXAMPLE_ANIMAL
    session = EXAMPLE_SESSION
    trial = EXAMPLE_TRIAL

    print(f"\nLoading {animal} / {session} / Trial {trial}...")
    data = load_trial(animal, session, trial)
    print(f"  fs = {data['fs']:.1f} Hz, duration = {data['time'][-1]:.1f} s")
    print(f"  fiber1: {data['fiber1'].shape}, fiber2: {data['fiber2'].shape}")
    print(f"  lfp_right: {'OK' if data['lfp_right'] is not None else 'MISSING'}")
    print(f"  lfp_left: {'OK' if data['lfp_left'] is not None else 'MISSING'}")
    print(f"  speed: {'OK' if data['speed'] is not None else 'MISSING'}")

    prefix = f"{animal}_{session}_Trial{trial}"
    prefix_phase = (
        f"{animal}_allSessions_example_{session}_Trial{trial}"
        if PHASE_AGGREGATE_ALL_SESSIONS
        else prefix
    )

    # Figure 1: Traces
    print("\nFigure 1: Representative traces...")
    f1 = fig_traces(data, animal, session, trial)
    for ext in ("pdf", "png"):
        f1.savefig(
            str(OUTPUT_DIR / f"{prefix}_01_traces.{ext}"),
            dpi=DPI,
            bbox_inches="tight",
            pad_inches=TRACES_SAVE_PAD_INCH,
            facecolor="white",
            edgecolor="none",
        )
    plt.close(f1)
    print("  Saved.")

    # Figure 2: TFR spectrograms
    print("Figure 2: Time-frequency spectrograms...")
    f2 = fig_spectrograms(data, animal, session, trial)
    for ext in ("pdf", "png"):
        f2.savefig(str(OUTPUT_DIR / f"{prefix}_02_spectrograms.{ext}"),
                   dpi=DPI, bbox_inches="tight", facecolor="white", edgecolor="none")
    plt.close(f2)
    print("  Saved.")

    # Figure 3: Cross-correlation (REST vs RUN) — single-trial + all-trials layout
    print("Figure 3: Fiber-fiber theta cross-correlation (REST vs RUN)...")
    prefix_xcorr = (
        f"{animal}_allSessions_example_{session}_Trial{trial}"
        if XCORR_AGGREGATE_ALL_SESSIONS
        else prefix
    )
    if XCORR_AGGREGATE_ALL_SESSIONS:
        print(
            f"  Xcorr all-trials: all sessions in RECORDINGS[{animal!r}] "
            f"(example heatmap still {session} trial {trial})."
        )
    f3 = fig_xcorr(data, animal, session, trial)
    _save_phase_locking_figure(
        f3, OUTPUT_DIR / f"{prefix_xcorr}_03_xcorr_theta_rest_run"
    )
    plt.close(f3)
    print("  Saved.")

    # Figure 3b: Unsplit cross-correlation (all epochs) + surrogate test
    print("Figure 3b: Fiber-fiber theta cross-correlation (all epochs, unsplit)...")
    f3b = fig_xcorr_unsplit(data, animal, session, trial)
    _save_phase_locking_figure(
        f3b, OUTPUT_DIR / f"{prefix_xcorr}_03b_xcorr_theta_unsplit"
    )
    plt.close(f3b)
    print("  Saved.")

    # Figure 4: Phase-locking — trial-equalized roses + example ΔF/F & Hilbert phase traces
    print("Figure 4: Bilateral theta phase-locking (roses + example traces)...")
    if PHASE_AGGREGATE_ALL_SESSIONS:
        print(
            f"  Phase roses: all sessions in RECORDINGS[{animal!r}] "
            f"(example trace still {session} trial {trial})."
        )
    with mpl.rc_context({"text.hinting": "none"}):
        f4 = fig_phase_locking(data, animal, session, trial)
    _save_phase_locking_figure(
        f4, OUTPUT_DIR / f"{prefix_phase}_04_phase_locking_hilbert"
    )
    plt.close(f4)
    print("  Saved.")

    # Figure 5: Bilateral spectral analysis (PSD, coherence, coherogram)
    print("Figure 5: Bilateral spectral analysis (PSD / coherence / coherogram)...")
    prefix_spec = (
        f"{animal}_allSessions_example_{session}_Trial{trial}"
        if SPEC5_AGGREGATE_ALL_SESSIONS
        else prefix
    )
    if SPEC5_AGGREGATE_ALL_SESSIONS:
        print(
            f"  Spectral all-trials: all sessions in RECORDINGS[{animal!r}] "
            f"(coherogram + example trial: {session})."
        )
    f5 = fig_spectral_bilateral(data, animal, session, trial)
    _save_phase_locking_figure(
        f5, OUTPUT_DIR / f"{prefix_spec}_05_spectral_bilateral"
    )
    plt.close(f5)
    print("  Saved.")

    print(f"\nDONE! All figures saved to {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
