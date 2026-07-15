"""
================================================================================
FIGURE 7 (panels E-H): dual-color, multi-animal fiber voltage imaging during
40 Hz DBS
================================================================================

Panels A-D (optical setup schematic, cage photo, FOV image, histology) were
assembled by hand and are not reproduced here.

  (E) GEVI (AceMNeon2) signals, raw (left) and 40 Hz-filtered (right), per
      fiber, during 40 Hz electrical stimulation in CA1.
  (F) Same as (E), for the simultaneously recorded mCherry reference signal.
  (G) Spectral power change relative to baseline for the GEVI signals (left:
      one line per fiber; right: 40 Hz peak quantification as a bar plot).
  (H) Same as (G), for the mCherry reference signal.

INPUT
-----
A single .mat file with one variable `all_traces`, shaped
(n_samples, n_channels, n_trials). n_channels = 2 * n_fibers: a GEVI channel
and a simultaneously-recorded reference (mCherry) channel per fiber.

Channel identity (GEVI vs reference) is NOT assumed from a fixed column
order. Instead it is auto-detected from the data itself: GEVI channels show
a rise in stimulation-frequency power during the stimulation window (that is
the entire point of the experiment -- somatic voltage entrainment to the
DBS pulse train), while the static reference fluorophore should not. This
keeps the script correct regardless of how the two color channels happen to
be interleaved/ordered in a given `all_traces` export -- see
`detect_gevi_channels()`. Fibers are then further ranked by entrainment
strength (`rank_fibers_by_entrainment()`) so the strongest (assumed CA1) are
grouped first and the weakest (assumed mPFC) last in panels G/H, matching
the manuscript's stated pattern -- again without assuming a fixed order.

USAGE:  python fig7_dualcolor_multianimal.py
================================================================================
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, welch

# common.py lives in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from common import load_matlab_struct  # noqa: E402

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT  # noqa: E402

# ==============================================================================
#  USER CONFIGURATION -- EDIT THESE FOR YOUR RECORDING
# ==============================================================================

# SET THIS: path to the .mat file containing the `all_traces` array
# (n_samples, n_channels, n_trials).
DATA_PATH = _LAB_DATA_ROOT / "Fig7_dual_color" / "all_traces_sens2.mat"

FS = 430.0                # imaging sampling rate (Hz) -- see Methods
STIM_FREQ = 40.0          # DBS stimulation frequency (Hz)
STIM_ONSET_SEC = 4.0      # SET THIS: stimulation onset within each trial (s)
STIM_DURATION_SEC = 1.0   # SET THIS: stimulation duration (s)
BAND_HALF_WIDTH_HZ = 2.0  # +/- Hz around STIM_FREQ for the filtered trace / band power

# Explicit 0-based channel indices for the GEVI and reference group. Leave as
# None to auto-split channels into two equal halves (first half = GEVI,
# second half = reference, matched fiber-for-fiber by within-half position)
# -- this matches how `all_traces_sens2.mat` is laid out. The auto-detection
# below (detect_gevi_channels) then VALIDATES that assumption against the
# data and swaps/warns if the "GEVI" half doesn't actually show stronger
# stimulus entrainment than the "reference" half.
GEVI_CHANNEL_IDX = None
REF_CHANNEL_IDX = None

# Per-fiber display labels, in the same order as GEVI_CHANNEL_IDX /
# REF_CHANNEL_IDX (or, if those are auto-detected, in entrainment-ranked
# order -- strongest first). SET THIS to match your recording; defaults to
# generic "Fiber N" labels. N_CA1_FIBERS controls the CA1/mPFC bar-plot
# grouping in panels G/H (first N_CA1_FIBERS fibers = "CA1", the rest =
# "mPFC").
FIBER_LABELS = None   # e.g. ['#1 CA1', '#2 CA1', '#3 CA1', '#1 mPFC', '#2 mPFC']
N_CA1_FIBERS = 3

FREQ_MIN, FREQ_MAX = 10, 70   # display range for panels G/H (Hz)

FIGURE_DPI = 300
FIGURE_FORMATS = ['png', 'pdf']
# NOTE: deliberately NOT named the same as this script's own source folder
# ("Fig7_dualcolor_multianimal") -- on case-insensitive filesystems (Windows,
# default macOS) "Figures/Fig7_dualcolor_multianimal" and this script's own
# "figures/Fig7_dualcolor_multianimal" source folder would resolve to the
# SAME directory, and generated output would land next to (and could
# overwrite) source files.
OUTPUT_DIR = PROJECT_ROOT / "Figures" / "Fig7_dualcolor_dbs"

COLOR_GEVI = (0.10, 0.45, 0.45)   # teal
COLOR_REF = (0.65, 0.10, 0.10)    # red
COLOR_STIM = (0.55, 0.30, 0.15)   # brown, matches the stim-pulse trace in the manuscript figure


# ==============================================================================
#  SIGNAL HELPERS
# ==============================================================================

def bandpass(x, fs, f_center, half_width, order=4):
    """Zero-phase Butterworth bandpass filter, f_center +/- half_width."""
    nyq = fs / 2.0
    lo = max(0.5, f_center - half_width) / nyq
    hi = min(nyq - 0.5, f_center + half_width) / nyq
    b, a = butter(order, [lo, hi], btype='band')
    return filtfilt(b, a, x, axis=0)


def welch_psd(x, fs, nperseg_sec=1.0):
    """Welch PSD of a 1D signal. Returns (freqs, psd)."""
    nperseg = min(len(x), int(round(fs * nperseg_sec)))
    nperseg = max(nperseg, 8)
    freqs, psd = welch(x, fs=fs, nperseg=nperseg)
    return freqs, psd


def band_power(freqs, psd, target_freq, half_width):
    band = (freqs >= target_freq - half_width) & (freqs <= target_freq + half_width)
    return float(np.mean(psd[band])) if np.any(band) else np.nan


def channel_entrainment_db(traces, channel_idx, fs, stim_freq, half_width,
                            baseline_slice, stim_slice):
    """
    Per-channel mean (across trials) stimulation-band power change, in dB:
    10*log10(P_stim / P_base). Used both to tell GEVI from reference
    channels and to rank fibers by entrainment strength (see
    detect_gevi_channels / rank_fibers_by_entrainment below).
    """
    gains = []
    for c in channel_idx:
        trial_db = []
        for t in range(traces.shape[2]):
            x = traces[:, c, t]
            f_base, p_base = welch_psd(x[baseline_slice], fs)
            f_stim, p_stim = welch_psd(x[stim_slice], fs)
            pb = band_power(f_base, p_base, stim_freq, half_width)
            ps = band_power(f_stim, p_stim, stim_freq, half_width)
            if pb > 0 and ps > 0:
                trial_db.append(10.0 * np.log10(ps / pb))
        gains.append(float(np.nanmean(trial_db)) if trial_db else np.nan)
    return np.array(gains)


def detect_gevi_channels(traces, fs, stim_freq, half_width, baseline_slice, stim_slice,
                          gevi_idx=None, ref_idx=None):
    """
    Split channels into a GEVI group and a reference group.

    If gevi_idx/ref_idx are given, use them directly. Otherwise split the
    channels into two equal halves (first half = GEVI-candidate, second
    half = reference-candidate) and validate the split using each channel's
    mean stimulation-band power change from baseline to stimulation, in dB
    (real GEVI entrainment; a static reference should show ~0 dB). If the
    "reference" half actually shows stronger average entrainment than the
    "GEVI" half, the two groups are swapped and a warning is printed, so the
    figure is labeled correctly regardless of how the channels happen to be
    ordered in `all_traces`.
    """
    n_channels = traces.shape[1]

    if gevi_idx is not None and ref_idx is not None:
        return list(gevi_idx), list(ref_idx)

    if n_channels % 2 != 0:
        raise ValueError(
            f"Expected an even number of channels (GEVI + reference pairs), "
            f"got {n_channels}. Set GEVI_CHANNEL_IDX / REF_CHANNEL_IDX explicitly."
        )
    n_fibers = n_channels // 2
    half_a = list(range(0, n_fibers))
    half_b = list(range(n_fibers, n_channels))

    gain_a = np.nanmean(channel_entrainment_db(traces, half_a, fs, stim_freq, half_width,
                                                baseline_slice, stim_slice))
    gain_b = np.nanmean(channel_entrainment_db(traces, half_b, fs, stim_freq, half_width,
                                                baseline_slice, stim_slice))
    print(f"  Channel-group stim-band entrainment: group A (idx {half_a}) = {gain_a:.2f} dB, "
          f"group B (idx {half_b}) = {gain_b:.2f} dB")

    if gain_b > gain_a:
        print("  NOTE: group B shows stronger stimulus entrainment than group A -- "
              "swapping so the more-entrained group is labeled GEVI.")
        return half_b, half_a
    return half_a, half_b


def rank_fibers_by_entrainment(traces, gevi_idx, ref_idx, fs, stim_freq, half_width,
                                baseline_slice, stim_slice):
    """
    Reorder (gevi_idx, ref_idx) together, strongest-to-weakest GEVI
    stimulation entrainment, so the manuscript's stated pattern (CA1 fibers
    show strong, consistent 40Hz entrainment; mPFC fibers show weaker
    entrainment) can be used to group bars as "CA1" (first N_CA1_FIBERS) vs
    "mPFC" (the rest) without assuming a specific channel order in the .mat
    file. This is a data-driven PROXY for anatomical identity, not ground
    truth -- if you know the true fiber order, set GEVI_CHANNEL_IDX /
    REF_CHANNEL_IDX / FIBER_LABELS explicitly instead and this ranking step
    is skipped.
    """
    gains = channel_entrainment_db(traces, gevi_idx, fs, stim_freq, half_width,
                                    baseline_slice, stim_slice)
    order = np.argsort(-gains)
    print(f"  Fiber entrainment ranking (dB, strongest first): "
          f"{[f'{g:.2f}' for g in gains[order]]}")
    gevi_sorted = [gevi_idx[i] for i in order]
    ref_sorted = [ref_idx[i] for i in order]
    return gevi_sorted, ref_sorted


# ==============================================================================
#  PLOTTING
# ==============================================================================

def plot_traces_panel(traces, channel_idx, fiber_labels, fs, onset_idx, offset_idx,
                       color, title, panel_letter, out_dir):
    """Panel E/F: raw (baseline-normalized) + stim-band-filtered traces, one row per fiber."""
    n_fibers = len(channel_idx)
    t = np.arange(traces.shape[0]) / fs

    fig, axes = plt.subplots(n_fibers, 2, figsize=(10, 1.6 * n_fibers), sharex=True)
    if n_fibers == 1:
        axes = axes.reshape(1, 2)

    for row, (ch, label) in enumerate(zip(channel_idx, fiber_labels)):
        trial_mean = np.mean(traces[:, ch, :], axis=1)
        baseline = trial_mean[:onset_idx]
        baseline_mean = np.mean(baseline)
        raw_norm = (trial_mean - baseline_mean) / baseline_mean
        filtered = bandpass(raw_norm, fs, STIM_FREQ, BAND_HALF_WIDTH_HZ)

        ax_raw, ax_filt = axes[row, 0], axes[row, 1]
        for ax, y in ((ax_raw, raw_norm), (ax_filt, filtered)):
            ax.plot(t, y, color=color, linewidth=1.0)
            ax.axvline(onset_idx / fs, color='k', linestyle='--', linewidth=0.8)
            ax.axvline(offset_idx / fs, color='k', linestyle='--', linewidth=0.8)
            ax.set_ylabel(label, rotation=0, ha='right', va='center', fontsize=9)
            ax.set_yticks([])
            for spine in ('top', 'right', 'left'):
                ax.spines[spine].set_visible(False)

        if row == 0:
            ax_raw.set_title('raw', fontsize=10, style='italic')
            ax_filt.set_title(f'filtered {int(STIM_FREQ)}Hz', fontsize=10, style='italic')

    axes[-1, 0].set_xlabel('Time (s)')
    axes[-1, 1].set_xlabel('Time (s)')
    fig.suptitle(f'({panel_letter}) {title}', fontsize=12, fontweight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, out_dir, f'Fig7{panel_letter}_traces')
    plt.close(fig)


def plot_spectral_panel(traces, channel_idx, fiber_labels, fs, onset_idx, offset_idx,
                         color, title, panel_letter, out_dir):
    """Panel G/H: spectral power change vs frequency (left) + 40Hz bar quantification (right)."""
    baseline_slice = slice(0, onset_idx)
    stim_slice = slice(onset_idx, offset_idx)
    n_fibers = len(channel_idx)

    freq_axis = None
    change_per_fiber = []
    peak_per_fiber = []
    for ch in channel_idx:
        trial_changes = []
        for t in range(traces.shape[2]):
            x = traces[:, ch, t]
            f_base, p_base = welch_psd(x[baseline_slice], fs)
            f_stim, p_stim = welch_psd(x[stim_slice], fs)
            # Interpolate the (shorter, stim-window) spectrum onto the
            # baseline frequency axis so trials/fibers can be averaged.
            p_stim_interp = np.interp(f_base, f_stim, p_stim)
            # dB change (10*log10(stim/baseline)), not a raw fractional
            # ratio: a fractional ratio blows up whenever baseline power at
            # a given frequency happens to be near zero, and dB is also
            # what Fig4/5 of this manuscript already use for stim-band
            # power quantification, so this keeps units consistent.
            with np.errstate(divide='ignore', invalid='ignore'):
                trial_db = 10.0 * np.log10(p_stim_interp / p_base)
            trial_changes.append(trial_db)
        change = np.nanmean(trial_changes, axis=0)
        if freq_axis is None:
            freq_axis = f_base
        change_per_fiber.append(change)
        peak_per_fiber.append(band_power(freq_axis, change, STIM_FREQ, BAND_HALF_WIDTH_HZ))

    fig, (ax_line, ax_bar) = plt.subplots(1, 2, figsize=(9, 3.5))

    freq_mask = (freq_axis >= FREQ_MIN) & (freq_axis <= FREQ_MAX)
    for change, label in zip(change_per_fiber, fiber_labels):
        ax_line.plot(freq_axis[freq_mask], change[freq_mask], color=color,
                     alpha=0.4 + 0.6 * np.random.default_rng(hash(label) % (2**31)).random(),
                     linewidth=1.5, label=label)
    ax_line.axvline(STIM_FREQ, color='0.4', linestyle=':', linewidth=1.0)
    ax_line.set_xlabel('Frequency (Hz)')
    ax_line.set_ylabel('Spectral Power Change (dB)')
    ax_line.set_title(f'({panel_letter}) {title}', fontsize=11, fontweight='bold')

    bar_colors = [color] * N_CA1_FIBERS + [tuple(0.6 * c + 0.4 for c in color)] * (n_fibers - N_CA1_FIBERS)
    ax_bar.bar(range(1, n_fibers + 1), peak_per_fiber, color=bar_colors[:n_fibers])
    ax_bar.axhline(0, color='k', linewidth=0.8)
    ax_bar.set_xticks(range(1, n_fibers + 1))
    ax_bar.set_ylabel(f'{int(STIM_FREQ)}Hz Power Change (dB)')
    n_mpfc = n_fibers - N_CA1_FIBERS
    if N_CA1_FIBERS > 0:
        ax_bar.text((1 + N_CA1_FIBERS) / 2, -0.05, 'CA1', ha='center', va='top',
                    transform=ax_bar.get_xaxis_transform())
    if n_mpfc > 0:
        ax_bar.text((N_CA1_FIBERS + 1 + n_fibers) / 2, -0.05, 'mPFC', ha='center', va='top',
                    transform=ax_bar.get_xaxis_transform())

    for ax in (ax_line, ax_bar):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    fig.tight_layout()
    _save(fig, out_dir, f'Fig7{panel_letter}_spectral')
    plt.close(fig)

    return freq_axis, change_per_fiber, peak_per_fiber


def _save(fig, out_dir, stem):
    out_dir.mkdir(parents=True, exist_ok=True)
    for fmt in FIGURE_FORMATS:
        path = out_dir / f'{stem}.{fmt}'
        fig.savefig(path, dpi=FIGURE_DPI, bbox_inches='tight')
    print(f'  Saved: {stem} ({", ".join(FIGURE_FORMATS)})')


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    print(f'Loading: {DATA_PATH}')
    data = load_matlab_struct(DATA_PATH)
    if data is None or 'all_traces' not in data:
        raise FileNotFoundError(
            f"Could not load variable 'all_traces' from {DATA_PATH}. "
            f"Set DATA_PATH at the top of this script."
        )

    traces = np.asarray(data['all_traces'], dtype=float)
    if traces.ndim != 3:
        raise ValueError(f"Expected all_traces to be 3D (samples x channels x trials), "
                          f"got shape {traces.shape}")
    n_samples, n_channels, n_trials = traces.shape
    print(f'  all_traces shape: {traces.shape} (samples x channels x trials)')

    onset_idx = int(round(STIM_ONSET_SEC * FS))
    offset_idx = int(round((STIM_ONSET_SEC + STIM_DURATION_SEC) * FS))
    if offset_idx >= n_samples or onset_idx <= 0:
        raise ValueError(
            f"STIM_ONSET_SEC ({STIM_ONSET_SEC}s) / STIM_DURATION_SEC ({STIM_DURATION_SEC}s) "
            f"don't fit inside the {n_samples / FS:.2f}s trial. Adjust these at the top of "
            f"this script to match your recording."
        )
    baseline_slice = slice(0, onset_idx)
    stim_slice = slice(onset_idx, offset_idx)

    gevi_idx, ref_idx = detect_gevi_channels(
        traces, FS, STIM_FREQ, BAND_HALF_WIDTH_HZ, baseline_slice, stim_slice,
        gevi_idx=GEVI_CHANNEL_IDX, ref_idx=REF_CHANNEL_IDX,
    )

    # Only re-rank fibers by entrainment strength when the channel groups
    # were auto-detected -- if GEVI_CHANNEL_IDX/REF_CHANNEL_IDX were set
    # explicitly, that order is assumed to already be correct (e.g. it may
    # reflect known anatomy) and is left untouched.
    if GEVI_CHANNEL_IDX is None or REF_CHANNEL_IDX is None:
        gevi_idx, ref_idx = rank_fibers_by_entrainment(
            traces, gevi_idx, ref_idx, FS, STIM_FREQ, BAND_HALF_WIDTH_HZ,
            baseline_slice, stim_slice,
        )

    print(f'  GEVI channels (ordered): {gevi_idx}')
    print(f'  Reference channels (ordered, paired with GEVI by position): {ref_idx}')

    n_fibers = len(gevi_idx)
    labels = FIBER_LABELS if FIBER_LABELS and len(FIBER_LABELS) == n_fibers \
        else [f'Fiber {i + 1}' for i in range(n_fibers)]

    plot_traces_panel(traces, gevi_idx, labels, FS, onset_idx, offset_idx,
                       COLOR_GEVI, 'GEVI (AceMNeon2)', 'E', OUTPUT_DIR)
    plot_traces_panel(traces, ref_idx, labels, FS, onset_idx, offset_idx,
                       COLOR_REF, 'REFERENCE (mCherry)', 'F', OUTPUT_DIR)
    plot_spectral_panel(traces, gevi_idx, labels, FS, onset_idx, offset_idx,
                         COLOR_GEVI, 'GEVI (AceMNeon2)', 'G', OUTPUT_DIR)
    plot_spectral_panel(traces, ref_idx, labels, FS, onset_idx, offset_idx,
                         COLOR_REF, 'REFERENCE (mCherry)', 'H', OUTPUT_DIR)

    print(f'\nDone. Figures saved to: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()
