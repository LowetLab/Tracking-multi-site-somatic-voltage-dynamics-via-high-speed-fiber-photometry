"""
Phase-Amplitude Coupling (PAC) analysis for bilateral fiber–LFP combinations.

Computes LFP theta-phase-aligned wavelet spectrograms, phase curves
(theta waveform + low-gamma amplitude), and Tort modulation index (MI):
    R LFP θ → R Fiber γ  (ipsi)
    L LFP θ → R Fiber γ  (contra)
    L LFP θ → L Fiber γ  (ipsi)
    R LFP θ → L Fiber γ  (contra)

Core computation follows the ../../phase_amplitude_coupling/ MATLAB scripts; data loading and
behaviour classification reuse multisite_fiber_analysis.py.

Output per combination: 2-row figure (LFP row + Fiber row) each with
  spectrogram (RUN | REST) + phase curves (RUN | REST).
Plus a comparison/statistics figure (MI ipsi vs contra per fiber).
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.ticker import ScalarFormatter
from scipy import signal as sig
from scipy.stats import ttest_1samp
from pathlib import Path
import warnings

warnings.filterwarnings("ignore")

# multisite_fiber_analysis.py (and its RECORDINGS cohort config) lives in the
# Fig6_bilateral_ca1/ figure folder, not here -- point sys.path at it.
_MULTISITE_DIR = Path(__file__).resolve().parents[1] / "Fig6_bilateral_ca1"
if str(_MULTISITE_DIR) not in sys.path:
    sys.path.insert(0, str(_MULTISITE_DIR))

from multisite_fiber_analysis import (
    load_trial, RECORDINGS, THETA_BAND,
    style_axis, _save_phase_locking_figure, _perform_paired_test,
    _add_significance_bracket, _phase_locking_plot_half_violin_R,
    FIG6_COMBINATIONS, FIG6_COLORS,
    FONT_SIZE_TITLE, FONT_SIZE_SUPTITLE, FONT_SIZE_LABEL,
    FONT_SIZE_TICK, AXIS_LINEWIDTH,
    PROJECT_ROOT,
)

# =============================================================================
# CONFIGURATION
# =============================================================================

OUTPUT_DIR = PROJECT_ROOT / "Figures" / "Multisite_Fiber_Analysis" / "PAC"

PAC_ANIMAL = "Animal02"
PAC_AGGREGATE_ALL_SESSIONS = True

WAVELET_FREQ_RANGE = (5, 90)
WAVELET_N_FREQS = 86
WAVELET_CYCLES = 5

SPEC_DISPLAY_MIN = 25
SPEC_DISPLAY_MAX = 90

PHASE_CENTER_FREQ = 7.0
CYCLES_PER_EPOCH = 2
N_PHASE_BINS = 36

GAMMA_BAND = (30, 60)

EPOCH_RUN_THRESHOLD_CMS = 2.0
EPOCH_REST_THRESHOLD_CMS = 0.1
EPOCH_MIN_FRACTION = 1.0

MI_N_PHASE_BINS = 18
MI_N_SURROGATES = 200

FIGSIZE_COMBO_INCH = (24, 13)
FIGSIZE_COMPARISON_INCH = (16, 6)

COLOR_THETA_LFP = (0.50, 0.15, 0.55)
COLOR_THETA_FIBER = (0.15, 0.30, 0.70)
COLOR_GAMMA_LFP = (0.80, 0.52, 0.12)
COLOR_GAMMA_FIBER = (0.12, 0.62, 0.38)
LW_CURVE = 1.8
SEM_ALPHA = 0.18

CMAP_LFP_SPEC = "inferno"
CMAP_FIBER_SPEC = "viridis"


# =============================================================================
# WAVELET SPECTROGRAM  (port of MATLAB compute_wavelet_spectrogram)
# =============================================================================

def compute_morlet_spectrogram(signal_1d, fs, freq_vector, n_cycles=WAVELET_CYCLES):
    """Morlet wavelet magnitude spectrogram [n_freqs × n_samples]."""
    x = np.asarray(signal_1d, dtype=float).ravel()
    n = len(x)
    spec = np.empty((len(freq_vector), n), dtype=float)
    for fi, f in enumerate(freq_vector):
        sigma_t = n_cycles / (2 * np.pi * f)
        t_half = int(np.ceil(3 * sigma_t * fs))
        t_wav = np.arange(-t_half, t_half + 1) / fs
        wavelet = (np.exp(2j * np.pi * f * t_wav)
                   * np.exp(-t_wav ** 2 / (2 * sigma_t ** 2)))
        wavelet /= np.sum(np.abs(wavelet))
        spec[fi, :] = np.abs(np.convolve(x, wavelet, mode="same"))
    return spec


# =============================================================================
# PHASE EXTRACTION
# =============================================================================

def extract_theta_phase_wavelet(signal_1d, fs, center_freq=PHASE_CENTER_FREQ,
                                n_cycles=WAVELET_CYCLES):
    """Instantaneous theta phase via Morlet wavelet at centre frequency."""
    x = np.asarray(signal_1d, dtype=float).ravel()
    sigma_t = n_cycles / (2 * np.pi * center_freq)
    t_half = int(np.ceil(3 * sigma_t * fs))
    t_wav = np.arange(-t_half, t_half + 1) / fs
    wavelet = (np.exp(2j * np.pi * center_freq * t_wav)
               * np.exp(-t_wav ** 2 / (2 * sigma_t ** 2)))
    wavelet /= np.sum(np.abs(wavelet))
    return np.angle(np.convolve(x, wavelet, mode="same"))


# =============================================================================
# EPOCH DETECTION
# =============================================================================

def detect_theta_epochs(phase, cycles_per_epoch=CYCLES_PER_EPOCH, min_len=10):
    """Detect multi-cycle theta epochs from phase wrapping points."""
    dph = np.diff(phase)
    wrap_pts = np.where(dph < -np.pi)[0] + 1
    wrap_pts = wrap_pts[wrap_pts < len(phase)]
    n_cyc = len(wrap_pts)
    n_ep = n_cyc // cycles_per_epoch
    epochs = []
    for ei in range(n_ep):
        s = int(wrap_pts[ei * cycles_per_epoch])
        ei_end = (ei + 1) * cycles_per_epoch
        if ei_end < n_cyc:
            e = int(wrap_pts[ei_end]) - 1
        elif epochs:
            avg = int(np.mean([b - a for a, b in epochs]))
            e = min(s + avg, len(phase) - 1)
        else:
            e = len(phase) - 1
        if e - s >= min_len:
            epochs.append((s, e))
    return epochs


# =============================================================================
# EPOCH BEHAVIOUR CLASSIFICATION
# =============================================================================

def classify_epoch_state(speed_epoch):
    """Per-epoch: 'run', 'rest', or 'excluded' (strict all-sample criterion)."""
    if speed_epoch is None or len(speed_epoch) == 0:
        return "excluded"
    v = speed_epoch[~np.isnan(speed_epoch)]
    if len(v) == 0:
        return "excluded"
    if np.mean(v > EPOCH_RUN_THRESHOLD_CMS) >= EPOCH_MIN_FRACTION:
        return "run"
    if np.mean(v < EPOCH_REST_THRESHOLD_CMS) >= EPOCH_MIN_FRACTION:
        return "rest"
    return "excluded"


# =============================================================================
# PHASE BINNING + MI
# =============================================================================

def phase_bin_spectrogram(ep_spec, ep_phase, edges):
    """Bin [F×T] spectrogram by phase → [F × n_bins]."""
    n_f, n_b = ep_spec.shape[0], len(edges) - 1
    out = np.full((n_f, n_b), np.nan)
    for bi in range(n_b):
        lo, hi = edges[bi], edges[bi + 1]
        m = (ep_phase >= lo) & (ep_phase < hi) if bi < n_b - 1 else (ep_phase >= lo) & (ep_phase <= hi)
        if np.any(m):
            out[:, bi] = np.mean(ep_spec[:, m], axis=1)
    return out


def compute_tort_mi(phase, amplitude, n_bins=MI_N_PHASE_BINS):
    """Tort modulation index (KL divergence from uniform, normalised)."""
    edges = np.linspace(-np.pi, np.pi, n_bins + 1)
    mean_amp = np.zeros(n_bins)
    for bi in range(n_bins):
        m = (phase >= edges[bi]) & (phase < edges[bi + 1]) if bi < n_bins - 1 else (phase >= edges[bi]) & (phase <= edges[bi + 1])
        if np.any(m):
            mean_amp[bi] = np.mean(amplitude[m])
    total = np.sum(mean_amp)
    if total <= 0:
        return 0.0
    p = mean_amp / total
    nz = p > 0
    kl = np.sum(p[nz] * np.log(p[nz] * n_bins))
    return kl / np.log(n_bins)


# =============================================================================
# WELFORD ACCUMULATOR HELPERS
# =============================================================================

def _empty_accum():
    """Empty accumulator for one (combo, state, signal_role) cell."""
    return {
        "spec_sum": np.zeros((WAVELET_N_FREQS, N_PHASE_BINS)),
        "spec_cnt": np.zeros(N_PHASE_BINS),
        "theta_sum": np.zeros(N_PHASE_BINS),
        "theta_sq": np.zeros(N_PHASE_BINS),
        "theta_cnt": np.zeros(N_PHASE_BINS),
        "gamma_sum": np.zeros(N_PHASE_BINS),
        "gamma_sq": np.zeros(N_PHASE_BINS),
        "gamma_cnt": np.zeros(N_PHASE_BINS),
        "n_epochs": 0,
    }


def _accum_epoch(acc, binned_spec, ep_theta_binned, ep_gamma_binned,
                 ep_bin_ok):
    """Add one epoch's binned data to the accumulator."""
    for bi in range(N_PHASE_BINS):
        if not np.isnan(binned_spec[0, bi]):
            acc["spec_sum"][:, bi] += binned_spec[:, bi]
            acc["spec_cnt"][bi] += 1
        if ep_bin_ok[bi]:
            acc["theta_sum"][bi] += ep_theta_binned[bi]
            acc["theta_sq"][bi] += ep_theta_binned[bi] ** 2
            acc["theta_cnt"][bi] += 1
            acc["gamma_sum"][bi] += ep_gamma_binned[bi]
            acc["gamma_sq"][bi] += ep_gamma_binned[bi] ** 2
            acc["gamma_cnt"][bi] += 1
    acc["n_epochs"] += 1


def _finalise_accum(acc):
    """Compute mean spec, theta/gamma mean ± 95% CI from accumulators."""
    cnt = acc["spec_cnt"]
    acc["mean_spec"] = np.where(cnt[None, :] > 0,
                                acc["spec_sum"] / cnt[None, :], 0.0)
    for prefix in ("theta", "gamma"):
        s = acc[f"{prefix}_sum"]
        sq = acc[f"{prefix}_sq"]
        c = acc[f"{prefix}_cnt"]
        acc[f"{prefix}_mean"] = np.where(c > 0, s / c, 0.0)
        var = np.where(c > 1, (sq - s ** 2 / np.clip(c, 1, None))
                       / np.clip(c - 1, 1, None), 0.0)
        acc[f"{prefix}_ci"] = np.where(
            c > 1, 1.96 * np.sqrt(np.maximum(var, 0) / c), 0.0)


# =============================================================================
# MAIN DATA COLLECTION
# =============================================================================

def collect_pac_data(animal, session=None):
    """
    Collect phase-aligned spectrogram/phase-curve data and per-trial MI
    for all 4 combinations × 2 states (run/rest) × 2 signal roles (lfp/fiber).
    """
    sessions_map = RECORDINGS.get(animal, {})
    if not sessions_map:
        return None
    session_items = ([(session, sessions_map[session])]
                     if session else list(sessions_map.items()))

    freq_vector = np.linspace(*WAVELET_FREQ_RANGE, WAVELET_N_FREQS)
    phase_bin_edges = np.linspace(-np.pi, np.pi, N_PHASE_BINS + 1)
    phase_bin_centers = (phase_bin_edges[:-1] + phase_bin_edges[1:]) / 2
    gamma_mask = (freq_vector >= GAMMA_BAND[0]) & (freq_vector <= GAMMA_BAND[1])

    combos = FIG6_COMBINATIONS
    combo_keys = [c["short"] for c in combos]

    data_out = {}
    for k in combo_keys:
        data_out[k] = {}
        for state in ("run", "rest"):
            data_out[k][state] = {
                "lfp": _empty_accum(),
                "fiber": _empty_accum(),
            }

    trial_mi = []
    n_loaded, n_failed = 0, 0

    for sess_key, info in session_items:
        for trial_num in range(1, int(info["n_trials"]) + 1):
            try:
                trial = load_trial(animal, sess_key, trial_num)
            except Exception:
                n_failed += 1
                continue
            n_loaded += 1
            fs = float(trial["fs"])
            speed_arr = (np.asarray(trial["speed"], dtype=float).ravel()
                         if trial.get("speed") is not None else None)

            print(f"    {sess_key} T{trial_num}: wavelets", end="", flush=True)

            raw = {
                "fiber1": np.asarray(trial["fiber1"], dtype=float).ravel() * 100,
                "fiber2": (np.asarray(trial["fiber2"], dtype=float).ravel() * 100
                           if trial.get("fiber2") is not None else None),
                "lfp_right": (np.asarray(trial["lfp_right"], dtype=float).ravel()
                              if trial.get("lfp_right") is not None else None),
                "lfp_left": (np.asarray(trial["lfp_left"], dtype=float).ravel()
                             if trial.get("lfp_left") is not None else None),
            }

            specs, phases, theta_filts = {}, {}, {}
            for key, s in raw.items():
                if s is None or not np.any(np.isfinite(s)):
                    continue
                specs[key] = compute_morlet_spectrogram(s, fs, freq_vector)
                if key.startswith("lfp"):
                    phases[key] = extract_theta_phase_wavelet(s, fs)
                bp = sig.butter(4, [THETA_BAND[0] / (fs / 2),
                                    THETA_BAND[1] / (fs / 2)], btype="band")
                theta_filts[key] = sig.filtfilt(*bp, s - np.nanmean(s))

            print(" → epochs", end="", flush=True)

            tp = {"session": sess_key, "trial_num": trial_num}

            for combo in combos:
                fib_key = combo["sig1_key"]
                lfp_key = combo["sig2_key"]
                short = combo["short"]

                lfp_phase = phases.get(lfp_key)
                fib_spec = specs.get(fib_key)
                lfp_spec = specs.get(lfp_key)

                if lfp_phase is None or fib_spec is None:
                    tp[short + "_run"] = np.nan
                    tp[short + "_rest"] = np.nan
                    continue

                epochs = detect_theta_epochs(lfp_phase)
                if not epochs:
                    tp[short + "_run"] = np.nan
                    tp[short + "_rest"] = np.nan
                    continue

                mi_ph_run, mi_am_run = [], []
                mi_ph_rest, mi_am_rest = [], []

                for (es, ee) in epochs:
                    sl = slice(es, ee + 1)
                    ep_phase = lfp_phase[sl]
                    ep_fib_spec = fib_spec[:, sl]
                    ep_lfp_spec = lfp_spec[:, sl] if lfp_spec is not None else None
                    ep_gamma_fib = np.mean(ep_fib_spec[gamma_mask, :], axis=0)
                    ep_speed = (speed_arr[sl]
                                if speed_arr is not None and ee < len(speed_arr)
                                else None)
                    state = classify_epoch_state(ep_speed)
                    if state == "excluded":
                        continue

                    bin_idx = np.clip(
                        np.digitize(ep_phase, phase_bin_edges) - 1,
                        0, N_PHASE_BINS - 1)

                    lfp_tf_full = theta_filts.get(lfp_key)
                    fib_tf_full = theta_filts.get(fib_key)
                    lfp_tf_ep = (lfp_tf_full[sl] if lfp_tf_full is not None
                                 else np.zeros(ee - es + 1))
                    fib_tf_ep = (fib_tf_full[sl] if fib_tf_full is not None
                                 else np.zeros(ee - es + 1))

                    for role, ep_s, ep_tf in [
                        ("lfp", ep_lfp_spec, lfp_tf_ep),
                        ("fiber", ep_fib_spec, fib_tf_ep),
                    ]:
                        if ep_s is None:
                            continue
                        binned = phase_bin_spectrogram(ep_s, ep_phase,
                                                       phase_bin_edges)
                        ep_gam = np.mean(ep_s[gamma_mask, :], axis=0)
                        th_b = np.zeros(N_PHASE_BINS)
                        gm_b = np.zeros(N_PHASE_BINS)
                        ok = np.zeros(N_PHASE_BINS, dtype=bool)
                        for bi in range(N_PHASE_BINS):
                            m = bin_idx == bi
                            if np.any(m):
                                th_b[bi] = np.mean(ep_tf[m])
                                gm_b[bi] = np.mean(ep_gam[m])
                                ok[bi] = True
                        _accum_epoch(data_out[short][state][role],
                                     binned, th_b, gm_b, ok)

                    if state == "run":
                        mi_ph_run.append(ep_phase)
                        mi_am_run.append(ep_gamma_fib)
                    else:
                        mi_ph_rest.append(ep_phase)
                        mi_am_rest.append(ep_gamma_fib)

                tp[short + "_run"] = (
                    compute_tort_mi(np.concatenate(mi_ph_run),
                                    np.concatenate(mi_am_run))
                    if mi_ph_run and sum(len(p) for p in mi_ph_run) > 20
                    else np.nan)
                tp[short + "_rest"] = (
                    compute_tort_mi(np.concatenate(mi_ph_rest),
                                    np.concatenate(mi_am_rest))
                    if mi_ph_rest and sum(len(p) for p in mi_ph_rest) > 20
                    else np.nan)

            trial_mi.append(tp)
            mi_str = " | ".join(
                f"{c['short']}={tp.get(c['short']+'_run', np.nan):.4f}/"
                f"{tp.get(c['short']+'_rest', np.nan):.4f}"
                for c in combos)
            print(f" → {mi_str}")

    for k in combo_keys:
        for state in ("run", "rest"):
            for role in ("lfp", "fiber"):
                _finalise_accum(data_out[k][state][role])

    print(f"  Collection done: {n_loaded} loaded, {n_failed} failed")
    return {
        "results": data_out,
        "trial_mi": trial_mi,
        "freq_vector": freq_vector,
        "phase_bin_centers_rad": phase_bin_centers,
        "n_loaded": n_loaded,
        "n_failed": n_failed,
    }


# =============================================================================
# PLOT HELPERS
# =============================================================================

def _phase_deg_2cyc(pbc_rad):
    """1-cycle rad bins → 2-cycle degree axis (0–720)."""
    d = np.degrees(pbc_rad) + 180
    return np.concatenate([d, d + 360])


def _tile_2cyc(a):
    return np.concatenate([a, a])


def _disp_mask(fv):
    return (fv >= SPEC_DISPLAY_MIN) & (fv <= SPEC_DISPLAY_MAX)


def _format_dual_y(ax, curve, color, ylabel, symmetric=False, side="left"):
    """Style one y-axis of a dual-axis phase-curve panel."""
    ax.set_ylabel(ylabel, color=color, fontsize=FONT_SIZE_LABEL - 3)
    ax.tick_params(axis="y", colors=color, labelsize=FONT_SIZE_TICK - 3)
    ax.spines[side].set_color(color)
    lo, hi = float(np.nanmin(curve)), float(np.nanmax(curve))
    if symmetric:
        b = max(abs(lo), abs(hi))
        pad = b * 0.15
        ax.set_ylim(-b - pad, b + pad)
    else:
        pad = (hi - lo) * 0.12 if hi > lo else 0.01
        ax.set_ylim(lo - pad, hi + pad)
    mag = max(abs(lo), abs(hi))
    exp = int(np.floor(np.log10(mag))) if mag > 0 else 0
    sc = 10 ** exp
    ticks = np.linspace(*ax.get_ylim(), 5)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{t / sc:.1f}" for t in ticks],
                       fontsize=FONT_SIZE_TICK - 4)
    if exp != 0:
        ax.annotate(f"×10$^{{{exp}}}$",
                    xy=(0 if side == "left" else 1, 1),
                    xycoords="axes fraction",
                    xytext=(0, 4), textcoords="offset points",
                    fontsize=FONT_SIZE_TICK - 4, color=color,
                    ha="right" if side == "left" else "left", va="bottom")


def _plot_spectrogram_row(fig, gs, row, pac_res_state_dict, fv, pbc_rad,
                          cmap, row_label, xticks_2c):
    """Plot one spectrogram + phase-curve row (RUN | REST)."""
    dm = _disp_mask(fv)
    phase_2c = _phase_deg_2cyc(pbc_rad)
    is_lfp = "LFP" in row_label

    theta_col = COLOR_THETA_LFP if is_lfp else COLOR_THETA_FIBER
    gamma_col = COLOR_GAMMA_LFP if is_lfp else COLOR_GAMMA_FIBER
    unit = "mV" if is_lfp else "%"

    vmin, vmax = None, None
    for state in ("run", "rest"):
        ms = pac_res_state_dict[state]["mean_spec"][dm, :]
        pos = ms[ms > 0]
        if pos.size:
            lo, hi = float(pos.min()), float(ms.max())
            vmin = lo if vmin is None else min(vmin, lo)
            vmax = hi if vmax is None else max(vmax, hi)
    if vmin is None:
        vmin, vmax = 0, 1

    im_ref = None
    for ci, (state, st_lbl) in enumerate([("run", "RUN"), ("rest", "REST")]):
        ax = fig.add_subplot(gs[row, ci])
        ms = pac_res_state_dict[state]["mean_spec"]
        s2c = np.hstack([ms, ms])
        s_d = s2c[dm, :]
        fv_d = fv[dm]
        ph_edges = np.linspace(0, 720, s_d.shape[1] + 1)
        fr_edges = np.concatenate([
            [fv_d[0] - (fv_d[1] - fv_d[0]) / 2],
            (fv_d[:-1] + fv_d[1:]) / 2,
            [fv_d[-1] + (fv_d[-1] - fv_d[-2]) / 2],
        ])
        im = ax.pcolormesh(ph_edges, fr_edges, s_d, shading="flat",
                           cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
        if im_ref is None:
            im_ref = im
        ax.axvline(360, color="white", ls="--", lw=1.2, alpha=0.65)
        ax.set_xlim(0, 720)
        ax.set_xticks(xticks_2c)
        ax.set_xlabel("LFP Theta Phase (deg.)", fontsize=FONT_SIZE_LABEL - 3)
        if ci == 0:
            ax.set_ylabel(f"{row_label}\nFrequency (Hz)",
                          fontsize=FONT_SIZE_LABEL - 3)
        else:
            ax.set_yticklabels([])
        ax.set_title(st_lbl, fontsize=FONT_SIZE_TITLE - 3,
                     fontweight="bold", pad=5)
        ax.tick_params(labelsize=FONT_SIZE_TICK - 3)
        for sp in ax.spines.values():
            sp.set_linewidth(AXIS_LINEWIDTH * 0.6)
        n_ep = pac_res_state_dict[state]["n_epochs"]
        ax.text(0.97, 0.03, f"n = {n_ep}", transform=ax.transAxes,
                fontsize=FONT_SIZE_TICK - 4, ha="right", va="bottom",
                color="white")

    cax = fig.add_subplot(gs[row, 2])
    if im_ref is not None:
        cb = fig.colorbar(im_ref, cax=cax)
        lbl = f"Wavelet mag. ({unit})" if is_lfp else "Signal magnitude (%)"
        cb.set_label(lbl, fontsize=FONT_SIZE_TICK - 3, rotation=270,
                     labelpad=12)
        cb.ax.tick_params(labelsize=FONT_SIZE_TICK - 4)
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_scientific(True)
        fmt.set_powerlimits((-2, 2))
        cb.ax.yaxis.set_major_formatter(fmt)
    else:
        cax.set_visible(False)

    for ci, (state, st_lbl) in enumerate([("run", "RUN"), ("rest", "REST")]):
        ax_l = fig.add_subplot(gs[row, 3 + ci])
        ax_r = ax_l.twinx()
        d = pac_res_state_dict[state]
        tm = _tile_2cyc(d["theta_mean"])
        tc = _tile_2cyc(d["theta_ci"])
        gm = _tile_2cyc(d["gamma_mean"])
        gc = _tile_2cyc(d["gamma_ci"])

        ax_l.plot(phase_2c, tm, color=theta_col, lw=LW_CURVE, zorder=3)
        if np.any(tc > 0):
            ax_l.fill_between(phase_2c, tm - tc, tm + tc,
                              color=theta_col, alpha=SEM_ALPHA, lw=0, zorder=2)
        ax_l.axhline(0, color=theta_col, lw=0.5, ls=":", alpha=0.4)
        ax_r.plot(phase_2c, gm, color=gamma_col, lw=LW_CURVE, zorder=2)
        if np.any(gc > 0):
            ax_r.fill_between(phase_2c, gm - gc, gm + gc,
                              color=gamma_col, alpha=SEM_ALPHA, lw=0, zorder=1)
        ax_l.axvline(360, color="#999999", lw=0.6, ls="--", alpha=0.45)
        ax_l.set_xlim(0, 720)
        ax_l.set_xticks(xticks_2c)
        ax_l.set_xlabel("LFP Theta Phase (deg.)", fontsize=FONT_SIZE_LABEL - 3)
        ax_l.tick_params(axis="x", labelsize=FONT_SIZE_TICK - 3)

        t_full = np.concatenate([tm + tc, tm - tc])
        g_full = np.concatenate([gm + gc, gm - gc])
        _format_dual_y(ax_l, t_full, theta_col,
                       f"θ waveform ({unit})", symmetric=True, side="left")
        _format_dual_y(ax_r, g_full, gamma_col,
                       f"Low γ amplitude ({unit})", side="right")
        ax_l.spines["top"].set_visible(False)
        ax_r.spines["top"].set_visible(False)
        ax_r.spines["left"].set_visible(False)
        ax_l.spines["bottom"].set_linewidth(AXIS_LINEWIDTH * 0.6)
        ax_l.spines["left"].set_linewidth(AXIS_LINEWIDTH * 0.6)
        ax_r.spines["right"].set_linewidth(AXIS_LINEWIDTH * 0.6)
        ax_l.set_title(st_lbl, fontsize=FONT_SIZE_TITLE - 3,
                       fontweight="bold", pad=5)
        ax_l.text(0.97, 0.04, f"n = {d['n_epochs']}",
                  transform=ax_l.transAxes, fontsize=FONT_SIZE_TICK - 4,
                  ha="right", va="bottom", color="#666666")


# =============================================================================
# COMBINATION FIGURE (2 rows × 5 cols: spec_run, spec_rest, cbar, curve_run, curve_rest)
# =============================================================================

def fig_pac_combination(pac_data, combo, animal):
    """One combination figure matching the reference 2-row layout."""
    short = combo["short"]
    lfp_key = combo["sig2_key"]
    lfp_label = "R LFP" if lfp_key == "lfp_right" else "L LFP"
    fib_label = "R Fiber" if combo["sig1_key"] == "fiber1" else "L Fiber"
    res = pac_data["results"][short]
    fv = pac_data["freq_vector"]
    pbc = pac_data["phase_bin_centers_rad"]

    fig = plt.figure(figsize=FIGSIZE_COMBO_INCH)
    gs = GridSpec(2, 5, figure=fig,
                  width_ratios=[1, 1, 0.06, 1, 1],
                  hspace=0.38, wspace=0.35,
                  left=0.06, right=0.95, top=0.88, bottom=0.08)
    xticks = [0, 180, 360, 540, 720]

    lfp_data = {st: res[st]["lfp"] for st in ("run", "rest")}
    fib_data = {st: res[st]["fiber"] for st in ("run", "rest")}

    _plot_spectrogram_row(fig, gs, 0, lfp_data, fv, pbc,
                          CMAP_LFP_SPEC, lfp_label, xticks)
    _plot_spectrogram_row(fig, gs, 1, fib_data, fv, pbc,
                          CMAP_FIBER_SPEC, fib_label, xticks)

    fig.suptitle(
        f"{animal} — {combo['label']}\n"
        f"LFP $\\theta$ phase → wavelet spectrogram & phase curve",
        fontsize=FONT_SIZE_SUPTITLE - 2, fontweight="bold", y=0.97)
    return fig


# =============================================================================
# COMPARISON / STATISTICS FIGURE
# =============================================================================

def fig_pac_comparison(pac_data, animal):
    """
    3-panel comparison: R Fiber ipsi-vs-contra MI, L Fiber, Laterality CI.
    Uses RUN-epoch MI values.
    """
    trial_mi = pac_data["trial_mi"]

    fig = plt.figure(figsize=FIGSIZE_COMPARISON_INCH)
    gs = GridSpec(1, 3, figure=fig, wspace=0.42,
                  left=0.07, right=0.97, top=0.84, bottom=0.14)

    def _paired(ipsi_key, contra_key, sfx="_run"):
        ia, ca = [], []
        for tp in trial_mi:
            vi = tp.get(ipsi_key + sfx, np.nan)
            vc = tp.get(contra_key + sfx, np.nan)
            if np.isfinite(vi) and np.isfinite(vc):
                ia.append(vi)
                ca.append(vc)
        return np.array(ia), np.array(ca)

    def _violin(ax, a, b, col_a, col_b, title, ylabel="MI"):
        if len(a) >= 3:
            _phase_locking_plot_half_violin_R(
                ax, a, b, positions=(0.7, 1.3), colors=(col_a, col_b))
            yvals = np.concatenate([a, b])
            yr = max(float(np.ptp(yvals)), abs(float(np.max(yvals))) * 0.05)
            test = _perform_paired_test(a, b)
            p = test.get("p_value", np.nan)
            d = test.get("effect_size", np.nan)
            bk = float(np.max(yvals)) + 0.08 * yr
            if np.isfinite(p):
                _add_significance_bracket(ax, 0.7, 1.3, bk, p,
                                          line_height=0.015 * yr,
                                          text_offset=0.005 * yr)
                top = bk + 0.14 * yr
            else:
                top = float(np.max(yvals)) + 0.15 * yr
            ax.set_ylim(float(np.min(yvals)) - yr * 0.35, top)
            p_s = f"p={p:.4f}" if np.isfinite(p) and p >= .001 else (
                f"p={p:.2e}" if np.isfinite(p) else "p=n/a")
            d_s = f"d={d:.2f}" if np.isfinite(d) else "d=n/a"
            ax.text(0.03, 0.02,
                    f"n={len(a)}, {test.get('test_used','')}\n"
                    f"{p_s}, Cohen's {d_s}",
                    transform=ax.transAxes, fontsize=FONT_SIZE_TICK - 3,
                    va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.25",
                              fc="white", ec="0.55", alpha=0.92))
        else:
            ax.text(0.5, 0.5, "Insufficient\npaired data",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=FONT_SIZE_LABEL - 2, color="gray")
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL - 1)
        ax.set_title(title, fontsize=FONT_SIZE_TITLE - 2,
                     fontweight="bold", pad=6)
        ax.set_xticks([0.7, 1.3])
        ax.set_xticklabels(["Ipsi", "Contra"], fontsize=FONT_SIZE_TICK - 1)
        ax.set_xlim(0.15, 1.85)
        style_axis(ax)

    ax0 = fig.add_subplot(gs[0, 0])
    ri, rc = _paired("R\u2013R", "R\u2013L")
    _violin(ax0, ri, rc,
            tuple(np.clip(FIG6_COLORS["R\u2013R"], 0, 1)),
            tuple(np.clip(FIG6_COLORS["R\u2013L"], 0, 1)),
            "R Fiber: Ipsi vs Contra")

    ax1 = fig.add_subplot(gs[0, 1])
    li, lc = _paired("L\u2013L", "L\u2013R")
    _violin(ax1, li, lc,
            tuple(np.clip(FIG6_COLORS["L\u2013L"], 0, 1)),
            tuple(np.clip(FIG6_COLORS["L\u2013R"], 0, 1)),
            "L Fiber: Ipsi vs Contra")

    ax2 = fig.add_subplot(gs[0, 2])
    r_ci, l_ci = [], []
    for tp in trial_mi:
        rr = tp.get("R\u2013R_run", np.nan)
        rl = tp.get("R\u2013L_run", np.nan)
        ll = tp.get("L\u2013L_run", np.nan)
        lr = tp.get("L\u2013R_run", np.nan)
        if np.isfinite(rr) and np.isfinite(rl) and (rr + rl) > 0:
            r_ci.append((rr - rl) / (rr + rl))
        if np.isfinite(ll) and np.isfinite(lr) and (ll + lr) > 0:
            l_ci.append((ll - lr) / (ll + lr))
    r_ci, l_ci = np.array(r_ci), np.array(l_ci)

    cols = [tuple(np.clip(FIG6_COLORS["R\u2013R"], 0, 1)),
            tuple(np.clip(FIG6_COLORS["L\u2013L"], 0, 1))]
    pos = [0.7, 1.3]
    for arr, col, p_pos, lbl in zip([r_ci, l_ci], cols, pos,
                                     ["R Fiber", "L Fiber"]):
        if arr.size >= 3:
            parts = ax2.violinplot([arr], positions=[p_pos],
                                   showmeans=False, showextrema=False)
            for pc in parts["bodies"]:
                pc.set_facecolor(col)
                pc.set_alpha(0.5)
            bp = ax2.boxplot([arr], positions=[p_pos], widths=0.12,
                             patch_artist=True, showfliers=False, zorder=3)
            for patch in bp["boxes"]:
                patch.set_facecolor("white")
                patch.set_edgecolor(col)
                patch.set_linewidth(1.5)
            for w in bp["whiskers"] + bp["caps"]:
                w.set_color(col)
                w.set_linewidth(1.5)
            for m in bp["medians"]:
                m.set_color(col)
                m.set_linewidth(2)
            jit = np.random.uniform(-0.04, 0.04, len(arr))
            ax2.scatter(p_pos + jit, arr, s=25, color=col, alpha=0.7,
                        edgecolors="white", linewidths=0.5, zorder=4)
            _, p_val = ttest_1samp(arr, 0)
            star = ("***" if p_val < .001 else "**" if p_val < .01
                    else "*" if p_val < .05 else "n.s.")
            ymax = float(np.max(np.abs(arr)))
            ax2.text(p_pos, np.max(arr) + 0.06 * ymax, star,
                     ha="center", va="bottom",
                     fontsize=FONT_SIZE_TICK, fontweight="bold")
    ax2.axhline(0, color="grey", ls="--", lw=1)
    ax2.set_ylabel("Laterality CI (MI)", fontsize=FONT_SIZE_LABEL - 1)
    ax2.set_title("Laterality Contrast Index",
                  fontsize=FONT_SIZE_TITLE - 2, fontweight="bold", pad=6)
    ax2.set_xticks(pos)
    ax2.set_xticklabels(["R Fiber", "L Fiber"], fontsize=FONT_SIZE_TICK - 1)
    ax2.set_xlim(0.15, 1.85)
    style_axis(ax2)

    fig.suptitle(
        f"{animal} — LFP $\\theta$ → Fiber low-$\\gamma$ PAC (Tort MI), "
        f"RUN epochs (n={pac_data['n_loaded']} trials)",
        fontsize=FONT_SIZE_SUPTITLE - 2, fontweight="bold", y=0.96)
    return fig


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("  MULTI-SITE PAC ANALYSIS")
    print("=" * 70)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    animal = PAC_ANIMAL
    print(f"\nCollecting PAC data for {animal} ...")
    pac_data = collect_pac_data(animal)
    if pac_data is None:
        print("  No data.")
        return

    prefix = f"{animal}_allSessions" if PAC_AGGREGATE_ALL_SESSIONS else animal

    for combo in FIG6_COMBINATIONS:
        short = combo["short"]
        safe = short.replace("\u2013", "-")
        print(f"\n  Figure: {combo['label']} ...")
        f = fig_pac_combination(pac_data, combo, animal)
        _save_phase_locking_figure(f, OUTPUT_DIR / f"{prefix}_PAC_{safe}")
        plt.close(f)
        print("    Saved.")

    print("\n  Comparison figure ...")
    fc = fig_pac_comparison(pac_data, animal)
    _save_phase_locking_figure(fc, OUTPUT_DIR / f"{prefix}_PAC_comparison")
    plt.close(fc)
    print("    Saved.")

    print(f"\nDONE → {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
