"""
Photobleaching control figure (Suppl. Fig. 1) -- signal-quality check across
repeated baseline trials for one animal.

Figure A (session-level, 3 x 6):
  Row 1: filtered fiber traces + fitted double-exponential decay (dotted)
  Row 2: instantaneous photobleaching rate (%/s)
  Row 3: theta-band (5-9 Hz) fiber-LFP coherence time series

Figure B (day-level drift summary):
  Across DAY1..DAY3, plot:
    - final cumulative bleaching (%)
    - mean theta coherence (5-9 Hz)

Inputs:
  - Preprocessed FiberPhotometryAnalysis trial files (for filtered traces)
  - FieldTrip single-trial outputs figure2_fieldtrip_trial*.mat (for coherence)
"""

from __future__ import annotations

import sys
from pathlib import Path
import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import h5py
from scipy.io import loadmat
from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# -----------------------------------------------------------------------------
# Configuration  -- EDIT THESE FOR YOUR OWN DATASET
# -----------------------------------------------------------------------------

BASE_DATA_ROOT = _LAB_DATA_ROOT / "FiberVoltageImaging"

COH_ROOT = (
    PROJECT_ROOT / "Figures" / "Spectral_data_outputs_artifact_cleaned"
    / "clear" / "single_trial"
)

OUTPUT_DIR = PROJECT_ROOT / "Figures" / "photobleaching_control"

ANIMAL = "Animal01"
SESSIONS = ["01_09_25-R1", "02_09_25-R1", "03_09_25-R1", "03_09_25-R2"]
DAY_TO_SESSIONS = {
    "DAY1": ["01_09_25-R1"],
    "DAY2": ["02_09_25-R1"],
    "DAY3": ["03_09_25-R1", "03_09_25-R2"],  # combine R1+R2 as requested
}
EXAMPLE_SESSION = "01_09_25-R1"
N_TRIALS_PER_SESSION = {
    "01_09_25-R1": 6,
    "02_09_25-R1": 6,
    "03_09_25-R1": 2,
    "03_09_25-R2": 2,
}
FIBER_INDEX = 0  # Python index (MATLAB fiber 1)

THETA_BAND = (5.0, 9.0)
THETA_COH_SMOOTH_SEC = 2.0

SAVE_DPI = 220


# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------

COLOR_FIBER = np.array([0.13, 0.57, 0.55])         # teal
COLOR_FIT = np.array([0.83, 0.45, 0.37])           # dusty coral
COLOR_RATE = np.array([0.45, 0.40, 0.74])          # purple-blue
COLOR_COH = np.array([0.86, 0.45, 0.17])           # orange

FONT_TITLE = 17
FONT_LABEL = 14
FONT_TICK = 12
DAY_BOX_WIDTH = 0.22
DAY_DOT_SIZE = 40


def style_axis(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.8)
    ax.spines["bottom"].set_linewidth(1.8)
    ax.tick_params(labelsize=FONT_TICK, width=1.5, length=5)


# -----------------------------------------------------------------------------
# Data loading helpers
# -----------------------------------------------------------------------------

def _double_exp(t: np.ndarray, a: float, tau1: float, b: float, tau2: float, c: float) -> np.ndarray:
    return a * np.exp(-t / tau1) + b * np.exp(-t / tau2) + c


def _trial_mat_path(animal: str, session: str, trial_num: int) -> Path:
    """
    Supports both layouts:
      1) session/TrialN_suffix_N/file.mat
      2) session/file.mat
    """
    sess_dir = BASE_DATA_ROOT / animal / "Fiber_Voltage_Processed" / session
    file_name = f"{animal}-{session}_Trial{trial_num}_FiberPhotometry_Analysis.mat"

    p_nested = sess_dir / f"Trial{trial_num}_fov1_baselineRecording_60sec_{trial_num}" / file_name
    if p_nested.exists():
        return p_nested

    p_direct = sess_dir / file_name
    if p_direct.exists():
        return p_direct

    return p_nested


def _coh_mat_path(animal: str, session: str, trial_num: int) -> Path:
    return COH_ROOT / animal / session / "data" / f"figure2_fieldtrip_trial{trial_num}.mat"


def _load_filtered_trace_h5(mat_path: Path, fiber_index: int = 0) -> tuple[np.ndarray, np.ndarray, float]:
    """FiberPhotometryAnalysis saved as MATLAB v7.3 (HDF5)."""
    with h5py.File(str(mat_path), "r") as f:
        if "FiberPhotometryAnalysis" not in f:
            raise KeyError(f"FiberPhotometryAnalysis not found in HDF5: {mat_path}")
        root = f["FiberPhotometryAnalysis"]
        t = np.asarray(root["time"]["time_vector_seconds"][()], dtype=float).reshape(-1)

        fs: float
        if "parameters" in root and "sampling_rate" in root["parameters"]:
            sr = np.asarray(root["parameters"]["sampling_rate"][()], dtype=float)
            fs = float(sr.item() if sr.size == 1 else sr.ravel()[0])
        elif "sampling_rate" in root["time"]:
            sr = np.asarray(root["time"]["sampling_rate"][()], dtype=float)
            fs = float(sr.item() if sr.size == 1 else sr.ravel()[0])
        else:
            fs = float(1.0 / np.median(np.diff(t))) if t.size > 1 else 500.0

        filt = np.asarray(root["signals"]["filtered_traces"][()], dtype=float)
        if filt.ndim == 1:
            y = filt.ravel()
        else:
            if filt.shape[0] == t.size:
                y = filt[:, fiber_index]
            elif filt.shape[1] == t.size:
                y = filt[fiber_index, :]
            else:
                raise ValueError(f"Unexpected filtered_traces shape {filt.shape} for t={t.size} (HDF5)")
        y = np.asarray(y, dtype=float).ravel()

    n = min(t.size, y.size)
    return t[:n], y[:n], fs


def load_filtered_trace(mat_path: Path, fiber_index: int = 0) -> tuple[np.ndarray, np.ndarray, float]:
    """
    Load time vector + filtered trace from FiberPhotometryAnalysis MAT.
    Classic .mat via scipy; v7.3 HDF5 via h5py.
    """
    try:
        d = loadmat(str(mat_path), simplify_cells=True)
    except NotImplementedError:
        return _load_filtered_trace_h5(mat_path, fiber_index)

    fpa = d.get("FiberPhotometryAnalysis", None)
    if fpa is None:
        raise KeyError(f"FiberPhotometryAnalysis not found: {mat_path}")

    t = np.asarray(fpa["time"]["time_vector_seconds"], dtype=float).ravel()
    fs = float(fpa["parameters"]["sampling_rate"]) if "parameters" in fpa else float(fpa["time"]["sampling_rate"])

    filt = np.asarray(fpa["signals"]["filtered_traces"], dtype=float)
    if filt.ndim == 1:
        y = filt.ravel()
    else:
        if filt.shape[0] == t.size:
            y = filt[:, fiber_index]
        elif filt.shape[1] == t.size:
            y = filt[fiber_index, :]
        else:
            raise ValueError(f"Unexpected filtered_traces shape {filt.shape} for t={t.size}")
    y = np.asarray(y, dtype=float).ravel()

    n = min(t.size, y.size)
    return t[:n], y[:n], fs


def _load_theta_coherence_h5(mat_path: Path, theta_band: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """figure2_fieldtrip_trial*.mat as HDF5 (v7.3)."""
    with h5py.File(str(mat_path), "r") as f:
        freq = np.asarray(f["freq"][()], dtype=float).ravel()
        time = np.asarray(f["time"][()], dtype=float).ravel()
        coh_key = "coh_fieldtrip" if "coh_fieldtrip" in f else ("coh_mscohere" if "coh_mscohere" in f else None)
        if coh_key is None:
            raise KeyError(f"No coherence matrix in HDF5: {mat_path.name}")
        coh = np.asarray(f[coh_key][()], dtype=float)
    return _collapse_theta_coh(freq, time, coh, theta_band)


def _collapse_theta_coh(
    freq: np.ndarray, time: np.ndarray, coh: np.ndarray, theta_band: tuple[float, float]
) -> tuple[np.ndarray, np.ndarray]:
    if coh.ndim != 2:
        raise ValueError(f"Unexpected coherence shape {coh.shape}")
    if coh.shape[0] == freq.size:
        coh_ft = coh
    elif coh.shape[1] == freq.size:
        coh_ft = coh.T
    else:
        raise ValueError(f"Cannot align coherence with freq: {coh.shape}, freq={freq.size}")
    mask = (freq >= theta_band[0]) & (freq <= theta_band[1])
    if not np.any(mask):
        raise ValueError("Theta mask empty for coherence frequencies.")
    theta_ts = np.nanmean(coh_ft[mask, :], axis=0)
    # Smooth coherence to reduce high-frequency fluctuations in the summary trace.
    if time.size > 1:
        dt = float(np.median(np.diff(time)))
        if dt > 0:
            w = max(3, int(round(THETA_COH_SMOOTH_SEC / dt)))
            theta_ts = uniform_filter1d(theta_ts, size=w, mode="nearest")
    n = min(theta_ts.size, time.size)
    return time[:n], theta_ts[:n]


def load_theta_coherence_timeseries(mat_path: Path, theta_band: tuple[float, float]) -> tuple[np.ndarray, np.ndarray]:
    """
    Load time-resolved coherence from figure2_fieldtrip_trial*.mat and collapse
    to theta-band mean coherence over time.
    """
    try:
        d = loadmat(str(mat_path), simplify_cells=True)
    except NotImplementedError:
        return _load_theta_coherence_h5(mat_path, theta_band)

    freq = np.asarray(d["freq"], dtype=float).ravel()
    time = np.asarray(d["time"], dtype=float).ravel()
    coh_key = "coh_fieldtrip" if "coh_fieldtrip" in d else ("coh_mscohere" if "coh_mscohere" in d else None)
    if coh_key is None:
        raise KeyError(f"No coherence matrix found in {mat_path.name}")
    coh = np.asarray(d[coh_key], dtype=float)
    return _collapse_theta_coh(freq, time, coh, theta_band)


# -----------------------------------------------------------------------------
# Metrics
# -----------------------------------------------------------------------------

def fit_double_exponential(t: np.ndarray, y: np.ndarray) -> dict:
    t0 = t - t[0]
    if t0.size < 10:
        raise ValueError("Trace too short for fit.")

    ymax = float(np.nanmax(y))
    ymin = float(np.nanmin(y))
    duration = max(float(t0[-1]), 1e-6)

    p0 = [0.5 * ymax, duration / 3.0, 0.5 * ymax, duration, ymin]
    bounds = ([0, 1e-6, 0, 1e-6, -np.inf], [np.inf, np.inf, np.inf, np.inf, np.inf])

    popt, _ = curve_fit(_double_exp, t0, y, p0=p0, bounds=bounds, maxfev=20000)
    yhat = _double_exp(t0, *popt)

    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = np.nan if ss_tot <= 1e-12 else 1.0 - ss_res / ss_tot
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))

    return {
        "params": popt,
        "fit_curve": yhat,
        "r2": r2,
        "rmse": rmse,
        "tau_fast": float(min(popt[1], popt[3])),
        "tau_slow": float(max(popt[1], popt[3])),
        "asymptote": float(popt[4]),
    }


def _double_exp_deriv(t: np.ndarray, a: float, tau1: float, b: float, tau2: float) -> np.ndarray:
    """Analytical derivative of a*exp(-t/tau1) + b*exp(-t/tau2) + c."""
    return -(a / tau1) * np.exp(-t / tau1) - (b / tau2) * np.exp(-t / tau2)


def compute_bleach_rate_metrics(t: np.ndarray, fit_params: np.ndarray) -> dict:
    """Compute bleaching metrics from the double-exponential fit.

    Parameters
    ----------
    t : array  – time vector (seconds, starting from 0)
    fit_params : array  – [a, tau1, b, tau2, c] from curve_fit
    """
    a, tau1, b, tau2, c = fit_params
    t0 = t - t[0]

    fit_curve = _double_exp(t0, *fit_params)
    fit_start = float(fit_curve[0])
    baseline = fit_start if abs(fit_start) > 1e-12 else 1.0

    rate_abs = _double_exp_deriv(t0, a, tau1, b, tau2)
    rate_pct = (rate_abs / baseline) * 100.0

    cumulative_pct = ((fit_curve - fit_start) / baseline) * 100.0
    return {
        "rate_pct": rate_pct,
        "cumulative_pct": cumulative_pct,
        "final_cumulative_pct": float(cumulative_pct[-1]),
        "mean_rate_pct_s": float(np.mean(rate_pct)),
        "baseline_ref": baseline,
    }


def collect_trial_bundle(animal: str, session: str, trial_num: int, fiber_index: int) -> dict | None:
    mat_path = _trial_mat_path(animal, session, trial_num)
    coh_path = _coh_mat_path(animal, session, trial_num)

    if not mat_path.exists():
        warnings.warn(f"Missing trial file: {mat_path}", stacklevel=2)
        return None
    if not coh_path.exists():
        warnings.warn(f"Missing coherence file: {coh_path}", stacklevel=2)
        return None

    t, y, fs = load_filtered_trace(mat_path, fiber_index=fiber_index)
    fit = fit_double_exponential(t, y)
    bleach = compute_bleach_rate_metrics(t, fit["params"])
    tc_t, tc = load_theta_coherence_timeseries(coh_path, THETA_BAND)

    theta_mean = float(np.nanmean(tc))

    return {
        "session": session,
        "trial": trial_num,
        "t": t,
        "filtered_trace": y,
        "fit_curve": fit["fit_curve"],
        "fit_r2": fit["r2"],
        "fit_rmse": fit["rmse"],
        "tau_fast": fit["tau_fast"],
        "tau_slow": fit["tau_slow"],
        "rate_pct": bleach["rate_pct"],
        "final_cumulative_pct": bleach["final_cumulative_pct"],
        "baseline_fluorescence": bleach["baseline_ref"],
        "tc_t": tc_t,
        "theta_coh_t": tc,
        "theta_coh_mean": theta_mean,
    }


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------

def _plot_summary_trend(ax: plt.Axes, values: list[float], color: np.ndarray,
                        ylabel: str, title: str) -> None:
    """Connected scatter showing metric trend across trials."""
    vals = np.asarray(values, dtype=float)
    x = np.arange(1, len(vals) + 1)
    edge_dark = np.clip(np.array(color) * 0.4, 0, 1)

    ax.plot(x, vals, color=edge_dark, lw=1.8, zorder=2)
    ax.scatter(x, vals, s=DAY_DOT_SIZE + 10, color=color, edgecolor="white",
               lw=1.0, zorder=3)
    ax.set_xlim(0.4, len(vals) + 0.6)
    ax.set_xticks(x)
    ax.set_xticklabels([str(i) for i in x], fontsize=FONT_TICK - 1)
    ax.set_xlabel("Trial", fontsize=FONT_TICK)
    ax.set_title(title, fontsize=FONT_TICK, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=FONT_TICK)
    style_axis(ax)


def fig_photobleach_session(trials: list[dict], animal: str, session: str) -> plt.Figure:
    from matplotlib.gridspec import GridSpec
    n_trials = len(trials)
    width_ratios = [4.2] * n_trials + [1.8]
    fig = plt.figure(figsize=(4.2 * n_trials + 2.4, 11.2))
    gs = GridSpec(3, n_trials + 1, figure=fig, hspace=0.38, wspace=0.35,
                  width_ratios=width_ratios)

    row_axes = {0: [], 1: [], 2: []}       # trial columns only

    for c, tr in enumerate(trials):
        # Row 1
        ax = fig.add_subplot(gs[0, c])
        ax.plot(tr["t"], tr["filtered_trace"], color=COLOR_FIBER, lw=1.6, label="Filtered trace")
        ax.plot(tr["t"], tr["fit_curve"], color=COLOR_FIT, lw=2.0, ls=":", label="Double-exp fit")
        style_axis(ax)
        ax.set_title(f"Trial {tr['trial']}", fontsize=FONT_TITLE, fontweight="bold")
        if c == 0:
            ax.set_ylabel("Fluorescence (a.u.)", fontsize=FONT_LABEL)
        txt = f"$R^2$={tr['fit_r2']:.3f}\n$\\tau_f$={tr['tau_fast']:.1f}s, $\\tau_s$={tr['tau_slow']:.1f}s"
        ax.text(
            0.02, 1.08, txt, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=FONT_TICK - 1, clip_on=False,
        )
        if c == n_trials - 1:
            ax.legend(loc="upper right", fontsize=FONT_TICK - 1, frameon=True, framealpha=0.95)
        row_axes[0].append(ax)

        # Row 2
        ax2 = fig.add_subplot(gs[1, c])
        ax2.plot(tr["t"], tr["rate_pct"], color=COLOR_RATE, lw=1.8)
        ax2.axhline(0.0, color="0.4", lw=1.0, ls="--", alpha=0.8)
        style_axis(ax2)
        if c == 0:
            ax2.set_ylabel("Bleaching rate (%/s)", fontsize=FONT_LABEL)
        ax2.text(
            0.02, 0.97, f"Final cumul. = {tr['final_cumulative_pct']:.2f}%",
            transform=ax2.transAxes, ha="left", va="top", fontsize=FONT_TICK - 1,
        )
        row_axes[1].append(ax2)

        # Row 3
        ax3 = fig.add_subplot(gs[2, c])
        ax3.plot(tr["tc_t"], tr["theta_coh_t"], color=COLOR_COH, lw=1.8)
        style_axis(ax3)
        ax3.set_ylim(0, 1)
        if c == 0:
            ax3.set_ylabel(r"Theta coherence (5-9 Hz)", fontsize=FONT_LABEL)
        ax3.set_xlabel("Time (s)", fontsize=FONT_LABEL)
        ax3.text(
            0.02, 0.97, f"Mean={tr['theta_coh_mean']:.3f}",
            transform=ax3.transAxes, ha="left", va="top", fontsize=FONT_TICK - 1,
        )
        row_axes[2].append(ax3)

        for ax_i in (ax, ax2, ax3):
            ax_i.xaxis.set_major_locator(MaxNLocator(5))

    # 7th column: summaries
    fluor_vals = [tr["baseline_fluorescence"] for tr in trials]
    cum_vals = [tr["final_cumulative_pct"] for tr in trials]
    coh_vals = [tr["theta_coh_mean"] for tr in trials]

    ax_s1 = fig.add_subplot(gs[0, n_trials])
    _plot_summary_trend(ax_s1, fluor_vals, COLOR_FIBER, "", "Summary")

    ax_s2 = fig.add_subplot(gs[1, n_trials])
    _plot_summary_trend(ax_s2, cum_vals, COLOR_RATE, "", "Summary")
    ax_s2.set_ylim(-5, 0)

    ax_s3 = fig.add_subplot(gs[2, n_trials])
    _plot_summary_trend(ax_s3, coh_vals, COLOR_COH, "", "Summary")

    # Unify y-axes: Rows 1 & 3 include summary; Row 2 autoscales per panel
    summary_axes = {0: ax_s1, 2: ax_s3}
    for row_idx in (0, 2):
        all_ax = row_axes[row_idx] + [summary_axes[row_idx]]
        ylims = [a.get_ylim() for a in all_ax]
        ymin = min(lo for lo, _ in ylims)
        ymax = max(hi for _, hi in ylims)
        for a in all_ax:
            a.set_ylim(ymin, ymax)

    fig.suptitle(
        f"{animal} {session} — Photobleaching control (filtered trace, bleaching rate, theta coherence)",
        fontsize=20, fontweight="bold", y=0.995
    )
    fig.subplots_adjust(top=0.93)
    return fig


def _plot_day_metric(ax: plt.Axes, by_day: dict[str, list[float]], ylabel: str, title: str, color: np.ndarray) -> None:
    days = list(by_day.keys())
    x = np.arange(len(days))

    day_medians = [float(np.nanmedian(by_day[d])) if len(by_day[d]) > 0 else np.nan for d in days]
    ax.plot(x, day_medians, color="0.65", lw=1.8, ls="-", alpha=0.6, zorder=1)

    for i, d in enumerate(days):
        vals = np.asarray(by_day[d], dtype=float)
        vals = vals[np.isfinite(vals)]
        if vals.size == 0:
            continue

        bp = ax.boxplot(
            [vals], positions=[i], widths=DAY_BOX_WIDTH,
            patch_artist=True, showfliers=False, zorder=2,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("white")
            patch.set_edgecolor(color)
            patch.set_linewidth(2.0)
        for whisker in bp["whiskers"]:
            whisker.set_color(color)
            whisker.set_linewidth(2.0)
        for cap in bp["caps"]:
            cap.set_color(color)
            cap.set_linewidth(2.0)
        for median in bp["medians"]:
            median.set_color("black")
            median.set_linewidth(2.5)

        jitter = (np.random.default_rng(100 + i).random(vals.size) - 0.5) * 0.14
        ax.scatter(
            np.full(vals.size, i) + jitter, vals,
            s=DAY_DOT_SIZE, color=color, edgecolor="white", lw=1.2, alpha=0.9, zorder=4
        )
        m = float(np.nanmean(vals))
        sem = float(np.nanstd(vals, ddof=1) / np.sqrt(vals.size)) if vals.size > 1 else 0.0
        ax.errorbar(i, m, yerr=sem, fmt="s", ms=6, color="black", capsize=3, lw=1.6, zorder=5)

    ax.set_xticks(x)
    ax.set_xticklabels(days, rotation=20, ha="right", fontsize=FONT_TICK)
    ax.set_ylabel(ylabel, fontsize=FONT_LABEL)
    ax.set_title(title, fontsize=FONT_TITLE, fontweight="bold")
    style_axis(ax)


def fig_day_drift(summary: list[dict], animal: str) -> plt.Figure:
    by_day_cum, by_day_theta, by_day_fluor = {}, {}, {}
    for day_label, sess_list in DAY_TO_SESSIONS.items():
        vals = [d for d in summary if d["session"] in sess_list]
        by_day_cum[day_label] = [v["final_cumulative_pct"] for v in vals]
        by_day_theta[day_label] = [v["theta_coh_mean"] for v in vals]
        by_day_fluor[day_label] = [v["baseline_fluorescence"] for v in vals]

    fig, axes = plt.subplots(1, 3, figsize=(18.5, 5.8))
    _plot_day_metric(
        axes[0], by_day_cum, "Final cumulative bleaching (%)",
        "End-of-trial cumulative bleaching across days", COLOR_RATE
    )
    _plot_day_metric(
        axes[1], by_day_theta, r"Mean theta coherence (5-9 Hz)",
        "Theta coherence across days", COLOR_COH
    )
    _plot_day_metric(
        axes[2], by_day_fluor, "Baseline fluorescence (a.u.)",
        "Absolute fluorescence across days", COLOR_FIBER
    )

    fig.suptitle(f"{animal} — Session-to-day photobleaching/coherence/fluorescence drift", fontsize=20, fontweight="bold")
    fig.tight_layout(rect=[0, 0.0, 1, 0.93])
    return fig


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Session figure (requested 6 trials of 16_09_25-R1)
    session_trials = []
    for tr in range(1, N_TRIALS_PER_SESSION[EXAMPLE_SESSION] + 1):
        out = collect_trial_bundle(ANIMAL, EXAMPLE_SESSION, tr, FIBER_INDEX)
        if out is not None:
            session_trials.append(out)
    if not session_trials:
        raise RuntimeError(f"No valid trials loaded for {ANIMAL} {EXAMPLE_SESSION}.")

    fig1 = fig_photobleach_session(session_trials, ANIMAL, EXAMPLE_SESSION)
    f1_base = OUTPUT_DIR / f"{ANIMAL}_{EXAMPLE_SESSION}_photobleaching_control"
    fig1.savefig(str(f1_base) + ".png", dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
    fig1.savefig(str(f1_base) + ".pdf", dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig1)

    # Day drift summary across requested sessions
    summary = []
    for sess in SESSIONS:
        n_trials = N_TRIALS_PER_SESSION.get(sess, 0)
        for tr in range(1, n_trials + 1):
            out = collect_trial_bundle(ANIMAL, sess, tr, FIBER_INDEX)
            if out is not None:
                summary.append(out)
    if not summary:
        raise RuntimeError("No valid trials loaded for day-drift summary.")

    fig2 = fig_day_drift(summary, ANIMAL)
    f2_base = OUTPUT_DIR / f"{ANIMAL}_session_to_day_drift_photobleaching"
    fig2.savefig(str(f2_base) + ".png", dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
    fig2.savefig(str(f2_base) + ".pdf", dpi=SAVE_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig2)

    print(f"Saved:\n  {f1_base}.png/.pdf\n  {f2_base}.png/.pdf")


if __name__ == "__main__":
    main()

