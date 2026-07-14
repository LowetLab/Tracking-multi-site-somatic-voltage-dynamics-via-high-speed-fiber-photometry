"""
================================================================================
Phase-Aligned Spectrogram Plotting Script
================================================================================

This script loads the MATLAB-generated phase-aligned spectrogram results from 
the run_nonrun analysis and creates publication-quality figures.

Behavioral state (running / non_running) is **not** re-classified here: it must already
be encoded in the input .mat files from the upstream pipeline. The canonical rule
matches ``compute_state_spectrogram.m`` / Python ``test_phase_aligned_spectrogram_per_r.py``:
**running** only if every velocity sample in the epoch is **> 2** cm/s; **rest**
(non_running) only if every sample is **< 0.1** cm/s; intermediate epochs are excluded
from the running and rest averages; artifact epochs may be excluded upstream.

For each animal, it generates:
  - LFP running and non-running phase-aligned spectrograms
  - Fiber1 (and Fiber2 if available) running and non-running spectrograms

Input files:
  Located in: Process/phase_aligned_spectrogram/{condition}/{animal}/run_nonrun/
  - LFP/running/PhaseAlignedSpectrogram_running.mat
  - LFP/non_running/PhaseAlignedSpectrogram_non_running.mat
  - Fiber1/running/PhaseAlignedSpectrogram_running.mat
  - Fiber1/non_running/PhaseAlignedSpectrogram_non_running.mat
  (etc. for additional fibers)

Output:
  - Figures/phase_aligned_spectrogram/ — LFP purple-gold; Fiber viridis + pale-mint high end (per-signal colormaps)
  - Figures/phase_aligned_spectrogram_viridis/ — same layout, matplotlib viridis for all panels (classic)

Phase Definition:
  The phase is extracted from 5-9 Hz band using Morlet wavelet (center frequency 7 Hz).
  Cosine convention (angle of complex analytic signal):
  - 0° = Trough | 90° = Rising zero-crossing | 180° = Peak | 270° = Falling zero-crossing
  Epochs start at the trough (phase = –π wraps back to –π, displayed as 0°).
    (transition from negative to positive)
  - 90° corresponds to the PEAK of the oscillation
  - 180° corresponds to the FALLING ZERO-CROSSING (transition from positive to negative)
  - 270° (or -90°) corresponds to the TROUGH of the oscillation
  
  The 2-cycle view shows 0° to 720° (two complete oscillation cycles).

Usage:
    python plot_phase_aligned_spectrogram.py
================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter
from scipy.io import loadmat
from pathlib import Path
from datetime import datetime
import warnings

# Suppress scipy loadmat warnings
warnings.filterwarnings('ignore', category=UserWarning)


# ==============================================================================
#  CONFIGURATION
# ==============================================================================

# Base directories (relative to this script's location)
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent  # Go up from Code/ to project root

PROCESS_DIR = PROJECT_DIR / "Process" / "phase_aligned_spectrogram"
FIGURES_DIR = PROJECT_DIR / "Figures" / "phase_aligned_spectrogram"
# Second tree: classic matplotlib viridis for all panels (no LFP vs Fiber colormap split)
FIGURES_DIR_VIRIDIS = PROJECT_DIR / "Figures" / "phase_aligned_spectrogram_viridis"
CMAP_VIRIDIS = "viridis"

# Figure settings
FIGURE_DPI = 300
FIGURE_FORMAT = ['png', 'pdf', 'svg']  # Save in multiple formats

# Per-signal colormaps — LFP unchanged; Fiber uses full viridis then a short pale-mint tail at the top
CMAP_LFP = LinearSegmentedColormap.from_list('lfp_purple_gold', [
    (0.12, 0.02, 0.18),   # very dark purple
    (0.35, 0.08, 0.42),   # dark purple
    (0.55, 0.20, 0.58),   # medium purple
    (0.75, 0.42, 0.48),   # warm mauve
    (0.88, 0.65, 0.30),   # golden
    (0.97, 0.88, 0.50),   # bright gold
], N=256)


def _make_fiber_viridis_mint(n=256, tail_frac=0.12, mint_rgb=(0.88, 0.97, 0.94)):
    """
    Full matplotlib viridis ramp for most of the colormap; last ``tail_frac`` blends
    from viridis yellow to near-white mint (extends the top without losing viridis below).
    """
    try:
        vir = plt.colormaps['viridis']
    except (AttributeError, TypeError):
        vir = plt.cm.get_cmap('viridis')
    n_tail = max(1, int(round(n * tail_frac)))
    n_vir = n - n_tail
    colors = []
    for i in range(n_vir):
        t = i / (n_vir - 1) if n_vir > 1 else 0.0
        rgba = vir(t)
        colors.append((float(rgba[0]), float(rgba[1]), float(rgba[2])))
    v_end = np.asarray(vir(1.0)[:3], dtype=float)
    mint = np.asarray(mint_rgb, dtype=float)
    for j in range(n_tail):
        t = j / (n_tail - 1) if n_tail > 1 else 1.0
        rgb = v_end * (1.0 - t) + mint * t
        colors.append(tuple(rgb.tolist()))
    return LinearSegmentedColormap.from_list('fiber_viridis_mint', colors, N=n)


CMAP_FIBER = _make_fiber_viridis_mint()


def _cmap_for_signal(signal_type):
    """Return the colormap for a given signal type."""
    return CMAP_LFP if signal_type == 'LFP' else CMAP_FIBER


def _spectrogram_cmap(signal_type, cmap_mode="per_signal"):
    """
    cmap_mode:
      'per_signal' — LFP purple-gold, Fiber viridis + mint tail (FIGURES_DIR)
      'viridis'    — matplotlib viridis for all (FIGURES_DIR_VIRIDIS)
    """
    if cmap_mode == "viridis":
        return CMAP_VIRIDIS
    return _cmap_for_signal(signal_type)


def _spectrogram_cmap_lfp_fiber_pair(cmap_mode="per_signal"):
    """Return (cmap_lfp, cmap_fiber) for LFP vs Fiber comparison figure."""
    if cmap_mode == "viridis":
        return CMAP_VIRIDIS, CMAP_VIRIDIS
    return CMAP_LFP, CMAP_FIBER

# ------------------------------------------------------------------------------
# Y-axis display range and non-linear stretch
# FREQ_DISPLAY_MIN/MAX: only data within this range is shown
# FREQ_STRETCH_BREAK:   breakpoint (Hz) between the stretched lower and compressed upper segment
# FREQ_STRETCH_RATIO:   fraction of Y-axis height allocated to [FREQ_DISPLAY_MIN, FREQ_STRETCH_BREAK]
# FREQ_STRETCH_TICKS:   which Hz values to label on the Y-axis
# Set FREQ_STRETCH = False to revert to a plain linear Y-axis over the full data range
# ------------------------------------------------------------------------------
FREQ_STRETCH         = True   # Enable non-linear Y-axis
FREQ_DISPLAY_MIN     = 25     # Hz - lower clip of displayed frequency range
FREQ_DISPLAY_MAX     = 90     # Hz - upper clip of displayed frequency range
FREQ_STRETCH_BREAK   = 60     # Hz - breakpoint between stretched and compressed regions
FREQ_STRETCH_RATIO   = 0.70   # Fraction of Y-axis height for [FREQ_DISPLAY_MIN, FREQ_STRETCH_BREAK]
FREQ_STRETCH_TICKS   = [25, 60, 90]  # Y-axis tick labels (Hz)

# Signal type labels for simple figure (right corner annotation)
SIGNAL_LABELS = {
    'LFP': 'LFP',
    'Fiber1': 'Pyr',      # Pyramidal neurons
    'Fiber2': 'Pyr2',
}
PV_CONDITION = 'PV_Animals'
PV_SIGNAL_LABELS = {'Fiber1': 'PV+', 'Fiber2': 'PV+'}


def get_signal_label(signal_type, condition):
    """Return display label (Pyr vs PV by condition)."""
    if condition == PV_CONDITION and signal_type in PV_SIGNAL_LABELS:
        return PV_SIGNAL_LABELS[signal_type]
    return SIGNAL_LABELS.get(signal_type, signal_type)


def _decode_mat_char(val):
    """Normalize MATLAB char / bytes / 0-d array from loadmat for string fields."""
    if val is None:
        return ''
    if isinstance(val, bytes):
        return val.decode('utf-8', errors='replace').strip()
    if isinstance(val, np.ndarray):
        if val.size == 0:
            return ''
        x = val.flatten()[0]
        if isinstance(x, bytes):
            return x.decode('utf-8', errors='replace').strip()
        return str(x).strip()
    return str(val).strip()


def _lfp_wavelet_is_z_scored(lfp_ephys_field):
    """True only for z-scored ephys fields (lfp_z, lfp_z_HP, lfp_z_*). Raw / aligned first in batch → False."""
    s = _decode_mat_char(lfp_ephys_field).lower()
    if not s:
        return False
    if s in ('lfp_z_hp', 'lfp_z'):
        return True
    if s.startswith('lfp_z_'):
        return True
    return False


def lfp_spectrogram_colorbar_label(lfp_ephys_field):
    """Colorbar label for LFP wavelet magnitude (matches MATLAB lfp_wavelet_colorbar_ylabel)."""
    if _lfp_wavelet_is_z_scored(lfp_ephys_field):
        return 'Wavelet mag. (a.u., z-scored)'
    return 'Wavelet mag. (mV)'


def lfp_spectrogram_colorbar_label_title_case(lfp_ephys_field):
    if _lfp_wavelet_is_z_scored(lfp_ephys_field):
        return 'Wavelet Mag. (a.u., z-scored)'
    return 'Wavelet Mag. (mV)'


# ==============================================================================
#  PUBLICATION FONT SIZES - Following reference code style
# ==============================================================================
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 20
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_COLORBAR = 14
FONT_SIZE_ANNOTATION = 12

# Axis line widths
AXIS_LINEWIDTH = 2.0
TICK_WIDTH = 1.8
TICK_LENGTH = 7


# ==============================================================================
#  COLOR DEFINITIONS - Following reference code style
# ==============================================================================
# Colors consistent with plot_figure1_gevi_lfp.py
COLOR_LFP = np.array([0.35, 0.25, 0.45])               # Purple-grey
COLOR_FIBER = np.array([0.127568, 0.566949, 0.550556])  # Teal (from viridis)


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

def clip_freq_data(freq_vector, spectrogram):
    """
    Clip spectrogram rows to [FREQ_DISPLAY_MIN, FREQ_DISPLAY_MAX].
    Returns the clipped freq_vector and matching spectrogram rows.
    When FREQ_STRETCH is False the full data range is kept.
    """
    if not FREQ_STRETCH:
        return freq_vector, spectrogram
    mask = (freq_vector >= FREQ_DISPLAY_MIN) & (freq_vector <= FREQ_DISPLAY_MAX)
    return freq_vector[mask], spectrogram[mask, :]


def build_nonlinear_freq_edges(freq_edges):
    """
    Remap frequency edges for a non-linear Y-axis.
    The lower segment [FREQ_DISPLAY_MIN, FREQ_STRETCH_BREAK] is stretched to
    occupy FREQ_STRETCH_RATIO of the display span; the upper segment is
    compressed.  Uses the canonical FREQ_DISPLAY_MIN/MAX so that the mapping
    is identical to apply_nonlinear_yticks (tick positions match pixel grid).
    """
    if not FREQ_STRETCH:
        return freq_edges
    f_min   = FREQ_DISPLAY_MIN
    f_max   = FREQ_DISPLAY_MAX
    f_brk   = FREQ_STRETCH_BREAK
    ratio   = FREQ_STRETCH_RATIO
    span    = f_max - f_min
    lo_span = f_brk - f_min
    hi_span = f_max - f_brk
    display = np.where(
        freq_edges <= f_brk,
        f_min + (freq_edges - f_min) / lo_span * ratio * span,
        f_min + ratio * span + (freq_edges - f_brk) / hi_span * (1.0 - ratio) * span,
    )
    return display


def apply_nonlinear_yticks(ax, freq_vector):
    """
    After plotting with remapped freq_edges, place yticks at the correct
    display positions and label them with the original Hz values.
    Uses FREQ_DISPLAY_MIN/MAX as the canonical reference so that the
    endpoint ticks (e.g. 25 Hz, 90 Hz) are never filtered out by floating-
    point rounding of freq_vector edges.
    """
    if not FREQ_STRETCH:
        return
    f_min   = FREQ_DISPLAY_MIN
    f_max   = FREQ_DISPLAY_MAX
    f_brk   = FREQ_STRETCH_BREAK
    ratio   = FREQ_STRETCH_RATIO
    span    = f_max - f_min
    lo_span = f_brk - f_min
    hi_span = f_max - f_brk

    ticks_hz = [t for t in FREQ_STRETCH_TICKS if f_min <= t <= f_max]
    ticks_display = []
    for f in ticks_hz:
        if f <= f_brk:
            pos = f_min + (f - f_min) / lo_span * ratio * span
        else:
            pos = f_min + ratio * span + (f - f_brk) / hi_span * (1.0 - ratio) * span
        ticks_display.append(pos)

    ax.set_yticks(ticks_display)
    ax.set_yticklabels([str(int(f)) for f in ticks_hz])
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_clip_on(False)
        tick.tick1line.set_clip_on(False)


def _tighten_spectrogram_axes(ax, phase_edges, freq_edges):
    """Remove matplotlib's default ~5% margins so pcolormesh meets the spines flush.

    Edge yticks (25 Hz, 90 Hz) are snapped to the ylim boundaries so labels
    sit right on the axis frame with no gap.
    """
    ax.set_xlim(phase_edges[0], phase_edges[-1])
    y_lo, y_hi = freq_edges[0], freq_edges[-1]
    ax.set_ylim(y_lo, y_hi)
    ax.margins(x=0, y=0)

    ticks = list(ax.get_yticks())
    labels = [t.get_text() for t in ax.get_yticklabels()]
    if ticks and labels and len(ticks) == len(labels):
        if ticks[0] < y_lo:
            ticks[0] = y_lo
        if ticks[-1] > y_hi:
            ticks[-1] = y_hi
        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)


def style_axis_publication(ax, remove_top_right=True):
    """
    Apply publication-quality styling to an axis.
    Following reference code conventions.
    """
    if remove_top_right:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    
    ax.tick_params(axis='both', which='major', 
                   width=TICK_WIDTH, length=TICK_LENGTH,
                   labelsize=FONT_SIZE_TICK)
    ax.tick_params(axis='both', which='minor', 
                   width=TICK_WIDTH * 0.7, length=TICK_LENGTH * 0.6)
    
    ax.grid(False)


def load_spectrogram_data(mat_path):
    """
    Load phase-aligned spectrogram data from MATLAB .mat file.
    
    Returns
    -------
    data : dict
        Dictionary with keys:
            'phase_bins_deg': Phase bin centers in degrees
            'freq_vector': Frequency vector (Hz)
            'mean_spectrogram': 1-cycle averaged spectrogram
            'two_cycle_spectrogram': 2-cycle concatenated spectrogram
            'n_epochs': Number of epochs used
            'mean_velocity': Mean velocity of epochs
    """
    mat_path = Path(mat_path)
    if not mat_path.exists():
        raise FileNotFoundError(f"MAT file not found: {mat_path}")
    
    raw = loadmat(str(mat_path), squeeze_me=True, struct_as_record=False)
    
    data = {}
    
    # Phase bins (in degrees)
    if 'phase_bins_deg' in raw:
        data['phase_bins_deg'] = np.asarray(raw['phase_bins_deg']).flatten()
    elif 'phase_bins' in raw:
        # Convert from radians if necessary
        phase_bins = np.asarray(raw['phase_bins']).flatten()
        if np.max(np.abs(phase_bins)) <= np.pi + 0.1:  # Likely radians
            data['phase_bins_deg'] = np.rad2deg(phase_bins)
        else:
            data['phase_bins_deg'] = phase_bins
    else:
        raise KeyError("Could not find phase_bins or phase_bins_deg in MAT file")
    
    # Frequency vector
    data['freq_vector'] = np.asarray(raw['freq_vector']).flatten()
    
    # Spectrograms
    data['mean_spectrogram'] = np.asarray(raw['mean_spectrogram'])
    
    if 'two_cycle_spectrogram' in raw:
        data['two_cycle_spectrogram'] = np.asarray(raw['two_cycle_spectrogram'])
        data['has_two_cycle'] = True
    elif 'cycle1_mean_spectrogram' in raw and 'cycle2_mean_spectrogram' in raw:
        cycle1 = np.asarray(raw['cycle1_mean_spectrogram'])
        cycle2 = np.asarray(raw['cycle2_mean_spectrogram'])
        data['two_cycle_spectrogram'] = np.hstack([cycle1, cycle2])
        data['has_two_cycle'] = True
    else:
        data['two_cycle_spectrogram'] = None
        data['has_two_cycle'] = False
    
    # Metadata
    data['n_epochs'] = int(raw.get('n_epochs', 0))
    data['mean_velocity'] = float(raw.get('mean_velocity', 0))
    data['lfp_ephys_field'] = _decode_mat_char(raw.get('lfp_ephys_field', ''))
    
    return data


def create_phase_axis_2cycles(phase_bins_deg):
    """
    Create phase axis for 2-cycle view (0 to 720 degrees).
    
    The original phase bins are typically -180 to 180.
    For 2-cycle view, we want:
      Cycle 1: -180 to 180 -> mapped to 0 to 360
      Cycle 2: -180 to 180 -> mapped to 360 to 720
    
    Returns
    -------
    phase_axis : ndarray
        Phase values from 0 to 720 degrees
    """
    # Shift from -180~180 to 0~360 for cycle 1
    cycle1_phase = phase_bins_deg + 180  # Now 0 to 360
    # Cycle 2 is cycle 1 + 360
    cycle2_phase = cycle1_phase + 360    # Now 360 to 720
    
    return np.concatenate([cycle1_phase, cycle2_phase])


def _resolve_spectrogram_for_plot(data):
    """Pick 2-cycle or 1-cycle spectrogram depending on data availability.

    Returns
    -------
    spec : ndarray [F x M]
    phase_max : float  (360 or 720)
    xticks : list
    show_cycle_boundary : bool
    """
    if data.get('has_two_cycle', False) and data.get('two_cycle_spectrogram') is not None:
        return data['two_cycle_spectrogram'], 720.0, [0, 180, 360, 540, 720], True
    return data['mean_spectrogram'], 360.0, [0, 90, 180, 270, 360], False


def create_spectrogram_figure_full(
    data, signal_type, running_state, animal_name, condition_name, output_dir, cmap_mode="per_signal"
):
    """
    Create FULL version of phase-aligned spectrogram figure with all details.
    Includes phase definition, metadata, and comprehensive labeling.
    
    Parameters
    ----------
    data : dict
        Loaded spectrogram data from load_spectrogram_data()
    signal_type : str
        'LFP' or 'Fiber1', 'Fiber2', etc.
    running_state : str
        'running' or 'non_running'
    animal_name : str
        Name of the animal for title
    condition_name : str
        Experimental condition (e.g., 'contraHPfiber_contraHPLFP')
    output_dir : Path
        Directory to save figures
    cmap_mode : str
        'per_signal' (LFP/Fiber colormaps) or 'viridis' (canonical single map)
    """
    # LFP units follow results.metadata.lfp_ephys_field (raw HP → mV; z-scored → a.u.)
    if signal_type == 'LFP':
        cbar_label = lfp_spectrogram_colorbar_label_title_case(data.get('lfp_ephys_field', ''))
    else:
        cbar_label = 'Signal Magnitude (%)'
    
    freq_vector = data['freq_vector']
    spec, phase_max, xticks, show_boundary = _resolve_spectrogram_for_plot(data)
    freq_vector, spec = clip_freq_data(freq_vector, spec)

    # Figure setup - larger figure with space for annotations
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # Plot spectrogram using pcolormesh
    phase_edges = np.linspace(0, phase_max, spec.shape[1] + 1)
    freq_edges = np.concatenate([[freq_vector[0] - (freq_vector[1] - freq_vector[0])/2],
                                  (freq_vector[:-1] + freq_vector[1:]) / 2,
                                  [freq_vector[-1] + (freq_vector[-1] - freq_vector[-2])/2]])
    freq_edges = build_nonlinear_freq_edges(freq_edges)

    im = ax.pcolormesh(phase_edges, freq_edges, spec,
                       shading='flat', cmap=_spectrogram_cmap(signal_type, cmap_mode))
    apply_nonlinear_yticks(ax, freq_vector)
    _tighten_spectrogram_axes(ax, phase_edges, freq_edges)

    # Colorbar with scientific notation
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(cbar_label, fontsize=FONT_SIZE_LABEL)
    cbar.ax.tick_params(labelsize=FONT_SIZE_TICK)
    
    # Apply scientific notation to colorbar
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_scientific(True)
    formatter.set_powerlimits((-2, 2))
    cbar.ax.yaxis.set_major_formatter(formatter)
    cbar.ax.yaxis.get_offset_text().set_fontsize(FONT_SIZE_TICK)
    
    # Add vertical line at 360° (cycle boundary) only for real 2-cycle data
    if show_boundary:
        ax.axvline(x=360, color='white', linestyle='--', linewidth=2, alpha=0.8)
    
    # Axis labels
    ax.set_xlabel('LFP Theta Phase (deg.)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    
    ax.set_xticks(xticks)
    ax.set_xlim([0, phase_max])
    
    # Title with metadata
    state_map = {'running': 'Running', 'non_running': 'Non-Running', 'all_epochs': 'All Epochs'}
    state_label = state_map.get(running_state, running_state)
    title = f'{animal_name} - {signal_type} Phase-Aligned Spectrogram ({state_label})\n'
    title += f'Condition: {condition_name} | n = {data["n_epochs"]} epochs | Mean velocity = {data["mean_velocity"]:.2f} cm/s'
    ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Style the axis
    style_axis_publication(ax, remove_top_right=False)
    
    # Add phase definition annotation at bottom
    phase_def_text = (
        "Phase Definition (5-9 Hz Wavelet, cosine convention): "
        "0° = Trough | 90° = Rising zero-crossing | 180° = Peak | 270° = Falling zero-crossing"
    )
    fig.text(0.5, 0.02, phase_def_text, ha='center', va='bottom', 
             fontsize=FONT_SIZE_ANNOTATION, style='italic', color='gray')
    
    # Adjust layout to make room for annotation
    plt.subplots_adjust(bottom=0.12)
    
    # Save in multiple formats
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    base_name = f'{animal_name}_{signal_type}_{running_state}_full'
    
    for fmt in FIGURE_FORMAT:
        save_path = output_dir / f'{base_name}.{fmt}'
        # SVG keeps transparent background, others use white
        if fmt == 'svg':
            fig.savefig(save_path, dpi=FIGURE_DPI, bbox_inches='tight',
                       facecolor='none', edgecolor='none', transparent=True)
        else:
            fig.savefig(save_path, dpi=FIGURE_DPI, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
        print(f'    Saved (full): {save_path.name}')
    
    plt.close(fig)


def create_spectrogram_figure_simple(
    data, signal_type, running_state, animal_name, condition_name, output_dir, cmap_mode="per_signal"
):
    """
    Create SIMPLE version of phase-aligned spectrogram figure.
    Square plot with minimal labeling, suitable for publication panels.
    Matches the reference figure style exactly.
    
    Parameters
    ----------
    data : dict
        Loaded spectrogram data from load_spectrogram_data()
    signal_type : str
        'LFP' or 'Fiber1', 'Fiber2', etc.
    running_state : str
        'running' or 'non_running'
    animal_name : str
        Name of the animal
    output_dir : Path
        Directory to save figures
    cmap_mode : str
        'per_signal' or 'viridis'
    """
    if signal_type == 'LFP':
        cbar_label = lfp_spectrogram_colorbar_label(data.get('lfp_ephys_field', ''))
        corner_label = 'LFP'
    else:
        cbar_label = 'Signal magnitude (%)'
        corner_label = get_signal_label(signal_type, condition_name)
    
    freq_vector = data['freq_vector']
    spec, phase_max, xticks_ph, show_boundary = _resolve_spectrogram_for_plot(data)
    freq_vector, spec = clip_freq_data(freq_vector, spec)

    # Figure setup
    fig = plt.figure(figsize=(5.0, 4.2))
    
    # Main axes position [left, bottom, width, height] in figure coordinates
    ax = fig.add_axes([0.15, 0.15, 0.60, 0.72])  # Square-ish main plot
    
    # Colorbar axes - slightly shorter than main plot, positioned to the right
    cax = fig.add_axes([0.78, 0.18, 0.03, 0.66])  # Colorbar (shorter than 0.72)
    
    # Plot spectrogram using pcolormesh
    phase_edges = np.linspace(0, phase_max, spec.shape[1] + 1)
    freq_edges = np.concatenate([[freq_vector[0] - (freq_vector[1] - freq_vector[0])/2],
                                  (freq_vector[:-1] + freq_vector[1:]) / 2,
                                  [freq_vector[-1] + (freq_vector[-1] - freq_vector[-2])/2]])
    freq_edges = build_nonlinear_freq_edges(freq_edges)

    im = ax.pcolormesh(phase_edges, freq_edges, spec,
                       shading='flat', cmap=_spectrogram_cmap(signal_type, cmap_mode))
    apply_nonlinear_yticks(ax, freq_vector)
    _tighten_spectrogram_axes(ax, phase_edges, freq_edges)

    # Colorbar
    cbar = fig.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=9)
    
    # Colorbar label - single line, vertical
    cbar.set_label(cbar_label, fontsize=9, rotation=270, labelpad=12)
    
    # Get data range for colorbar ticks (bottom, 2-3 middle, top)
    vmin, vmax = im.get_clim()
    # Create 5 ticks: min, 25%, 50%, 75%, max
    cbar_ticks = np.linspace(vmin, vmax, 5)
    cbar.set_ticks(cbar_ticks)
    
    # Determine the scale factor for scientific notation
    max_abs = max(abs(vmin), abs(vmax))
    if max_abs > 0:
        exponent = int(np.floor(np.log10(max_abs)))
    else:
        exponent = 0
    scale_factor = 10 ** exponent
    
    # Scale the ticks and round to 1 decimal place
    scaled_ticks = cbar_ticks / scale_factor
    tick_labels = [f'{t:.1f}' for t in scaled_ticks]
    cbar.set_ticklabels(tick_labels)
    
    # Add exponent label at top of colorbar, left-aligned with colorbar left edge
    if exponent != 0:
        exp_label = f'×10$^{{{exponent}}}$'
        # Use figure coordinates: left edge of colorbar, higher y position
        cbar_pos = cax.get_position()
        fig.text(cbar_pos.x0, cbar_pos.y1 + 0.03, exp_label,
                 fontsize=9, ha='left', va='bottom')
    
    # Axis labels
    ax.set_xlabel('LFP Theta Phase (deg.)', fontsize=10)
    ax.set_ylabel('Frequency (Hz)', fontsize=10)
    
    ax.set_xticks(xticks_ph)
    ax.set_xlim([0, phase_max])
    ax.tick_params(labelsize=9)
    
    ax.text(0.02, 1.02, corner_label, transform=ax.transAxes,
            fontsize=13, fontweight='bold', ha='left', va='bottom')
    
    # Style the axis
    ax.spines['bottom'].set_linewidth(1.0)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['top'].set_linewidth(1.0)
    ax.spines['right'].set_linewidth(1.0)
    
    # Save in multiple formats
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    base_name = f'{animal_name}_{signal_type}_{running_state}_simple'
    
    for fmt in FIGURE_FORMAT:
        save_path = output_dir / f'{base_name}.{fmt}'
        # SVG keeps transparent background, others use white
        if fmt == 'svg':
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='none', edgecolor='none', transparent=True)
        else:
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='white', edgecolor='none')
        print(f'    Saved (simple): {save_path.name}')
    
    plt.close(fig)


def create_spectrogram_comparison_figure(
    data_running, data_rest, signal_type, animal_name, condition_name, output_dir, cmap_mode="per_signal"
):
    """
    Create comparison figure showing Running vs Rest side by side with shared colorbar.
    
    Parameters
    ----------
    data_running : dict
        Loaded spectrogram data for running state
    data_rest : dict
        Loaded spectrogram data for rest/non-running state
    signal_type : str
        'LFP' or 'Fiber1', 'Fiber2', etc.
    animal_name : str
        Name of the animal
    output_dir : Path
        Directory to save figures
    cmap_mode : str
        'per_signal' or 'viridis'
    """
    if signal_type == 'LFP':
        cbar_label = lfp_spectrogram_colorbar_label(data_running.get('lfp_ephys_field', ''))
        corner_label = 'LFP'
    else:
        cbar_label = 'Signal magnitude (%)'
        corner_label = get_signal_label(signal_type, condition_name)
    
    freq_vector = data_running['freq_vector']
    spec_run, phase_max, xticks_ph, _ = _resolve_spectrogram_for_plot(data_running)
    spec_rst, _, _, _ = _resolve_spectrogram_for_plot(data_rest)
    freq_vector, spec_run = clip_freq_data(freq_vector, spec_run)
    _, spec_rst = clip_freq_data(data_rest['freq_vector'], spec_rst)

    vmin = min(np.nanmin(spec_run), np.nanmin(spec_rst))
    vmax = max(np.nanmax(spec_run), np.nanmax(spec_rst))
    
    fig = plt.figure(figsize=(9.0, 4.2))
    ax1 = fig.add_axes([0.08, 0.15, 0.35, 0.72])
    ax2 = fig.add_axes([0.48, 0.15, 0.35, 0.72])
    cax = fig.add_axes([0.86, 0.18, 0.02, 0.66])
    
    phase_edges = np.linspace(0, phase_max, spec_run.shape[1] + 1)
    freq_edges = np.concatenate([[freq_vector[0] - (freq_vector[1] - freq_vector[0])/2],
                                  (freq_vector[:-1] + freq_vector[1:]) / 2,
                                  [freq_vector[-1] + (freq_vector[-1] - freq_vector[-2])/2]])
    freq_edges = build_nonlinear_freq_edges(freq_edges)

    sig_cmap = _spectrogram_cmap(signal_type, cmap_mode)
    im1 = ax1.pcolormesh(phase_edges, freq_edges, spec_run,
                         shading='flat', cmap=sig_cmap, vmin=vmin, vmax=vmax)
    ax1.set_xlabel('LFP Theta Phase (deg.)', fontsize=10)
    ax1.set_ylabel('Frequency (Hz)', fontsize=10)
    ax1.set_xticks(xticks_ph)
    ax1.set_xlim([0, phase_max])
    ax1.tick_params(labelsize=9)
    ax1.set_title('Running', fontsize=11, fontweight='bold')
    apply_nonlinear_yticks(ax1, freq_vector)
    _tighten_spectrogram_axes(ax1, phase_edges, freq_edges)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.0)

    im2 = ax2.pcolormesh(phase_edges, freq_edges, spec_rst,
                         shading='flat', cmap=sig_cmap, vmin=vmin, vmax=vmax)
    ax2.set_xlabel('LFP Theta Phase (deg.)', fontsize=10)
    ax2.set_ylabel('')
    ax2.set_yticks([])
    ax2.set_xticks(xticks_ph)
    ax2.set_xlim([0, phase_max])
    ax2.tick_params(labelsize=9)
    ax2.set_title('Rest', fontsize=11, fontweight='bold')
    _tighten_spectrogram_axes(ax2, phase_edges, freq_edges)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.0)
    
    # Shared colorbar
    cbar = fig.colorbar(im1, cax=cax)
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label(cbar_label, fontsize=9, rotation=270, labelpad=12)
    
    # Colorbar ticks with scientific notation
    cbar_ticks = np.linspace(vmin, vmax, 5)
    cbar.set_ticks(cbar_ticks)
    
    max_abs = max(abs(vmin), abs(vmax))
    if max_abs > 0:
        exponent = int(np.floor(np.log10(max_abs)))
    else:
        exponent = 0
    scale_factor = 10 ** exponent
    
    scaled_ticks = cbar_ticks / scale_factor
    tick_labels = [f'{t:.1f}' for t in scaled_ticks]
    cbar.set_ticklabels(tick_labels)
    
    if exponent != 0:
        exp_label = f'×10$^{{{exponent}}}$'
        cbar_pos = cax.get_position()
        fig.text(cbar_pos.x0, cbar_pos.y1 + 0.03, exp_label,
                 fontsize=9, ha='left', va='bottom')
    
    # Add corner label (signal type) at top left
    ax1.text(0.02, 1.02, corner_label, transform=ax1.transAxes,
             fontsize=13, fontweight='bold', ha='left', va='bottom')
    
    # Save in multiple formats
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    base_name = f'{animal_name}_{signal_type}_running_vs_rest'
    
    for fmt in FIGURE_FORMAT:
        save_path = output_dir / f'{base_name}.{fmt}'
        if fmt == 'svg':
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='none', edgecolor='none', transparent=True)
        else:
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='white', edgecolor='none')
        print(f'    Saved (comparison): {save_path.name}')
    
    plt.close(fig)


def create_lfp_fiber_comparison_figure(
    data_lfp, data_fiber, running_state, animal_name, condition_name, output_dir, cmap_mode="per_signal"
):
    """
    Create comparison figure showing LFP vs Fiber side by side with SEPARATE colorbars.
    (LFP uses mV or z-scored a.u. from lfp_ephys_field; Fiber uses %)
    
    Parameters
    ----------
    data_lfp : dict
        Loaded spectrogram data for LFP
    data_fiber : dict
        Loaded spectrogram data for Fiber
    running_state : str
        'running' or 'non_running'
    animal_name : str
        Name of the animal
    output_dir : Path
        Directory to save figures
    cmap_mode : str
        'per_signal' (LFP/Fiber maps) or 'viridis' (both panels use viridis)
    """
    cmap_lfp, cmap_fiber = _spectrogram_cmap_lfp_fiber_pair(cmap_mode)
    freq_vector = data_lfp['freq_vector']
    spec_lfp, phase_max, xticks_ph, _ = _resolve_spectrogram_for_plot(data_lfp)
    spec_fiber, _, _, _ = _resolve_spectrogram_for_plot(data_fiber)
    freq_vector, spec_lfp = clip_freq_data(freq_vector, spec_lfp)
    _, spec_fiber = clip_freq_data(data_fiber['freq_vector'], spec_fiber)

    fig = plt.figure(figsize=(11.0, 4.2))
    ax1 = fig.add_axes([0.06, 0.15, 0.32, 0.72])
    cax1 = fig.add_axes([0.39, 0.18, 0.015, 0.66])
    ax2 = fig.add_axes([0.52, 0.15, 0.32, 0.72])
    cax2 = fig.add_axes([0.85, 0.18, 0.015, 0.66])
    
    phase_edges = np.linspace(0, phase_max, spec_lfp.shape[1] + 1)
    freq_edges = np.concatenate([[freq_vector[0] - (freq_vector[1] - freq_vector[0])/2],
                                  (freq_vector[:-1] + freq_vector[1:]) / 2,
                                  [freq_vector[-1] + (freq_vector[-1] - freq_vector[-2])/2]])
    freq_edges = build_nonlinear_freq_edges(freq_edges)

    im1 = ax1.pcolormesh(phase_edges, freq_edges, spec_lfp,
                         shading='flat', cmap=cmap_lfp)
    ax1.set_xlabel('LFP Theta Phase (deg.)', fontsize=10)
    ax1.set_ylabel('Frequency (Hz)', fontsize=10)
    ax1.set_xticks(xticks_ph)
    ax1.set_xlim([0, phase_max])
    ax1.tick_params(labelsize=9)
    ax1.set_title('LFP', fontsize=11, fontweight='bold')
    apply_nonlinear_yticks(ax1, freq_vector)
    _tighten_spectrogram_axes(ax1, phase_edges, freq_edges)
    for spine in ax1.spines.values():
        spine.set_linewidth(1.0)

    cbar1 = fig.colorbar(im1, cax=cax1)
    cbar1.ax.tick_params(labelsize=8)
    cbar1.set_label(lfp_spectrogram_colorbar_label(data_lfp.get('lfp_ephys_field', '')),
                    fontsize=8, rotation=270, labelpad=10)
    _format_colorbar_scientific(cbar1, cax1, fig, im1.get_clim())

    im2 = ax2.pcolormesh(phase_edges, freq_edges, spec_fiber,
                         shading='flat', cmap=cmap_fiber)
    ax2.set_xlabel('LFP Theta Phase (deg.)', fontsize=10)
    ax2.set_ylabel('Frequency (Hz)', fontsize=10)
    ax2.set_xticks(xticks_ph)
    ax2.set_xlim([0, phase_max])
    ax2.tick_params(labelsize=9)
    fiber_label = get_signal_label('Fiber1', condition_name)
    ax2.set_title(fiber_label, fontsize=11, fontweight='bold')
    apply_nonlinear_yticks(ax2, freq_vector)
    _tighten_spectrogram_axes(ax2, phase_edges, freq_edges)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.0)
    
    # Fiber colorbar
    cbar2 = fig.colorbar(im2, cax=cax2)
    cbar2.ax.tick_params(labelsize=8)
    cbar2.set_label('Signal magnitude (%)', fontsize=8, rotation=270, labelpad=10)
    _format_colorbar_scientific(cbar2, cax2, fig, im2.get_clim())
    
    # Add state label as suptitle
    state_map = {'running': 'Running', 'non_running': 'Rest', 'all_epochs': 'All Epochs'}
    state_label = state_map.get(running_state, running_state)
    fig.suptitle(f'{animal_name} - {state_label}', fontsize=12, fontweight='bold', y=0.98)
    
    # Save in multiple formats
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    base_name = f'{animal_name}_LFP_vs_Fiber_{running_state}'
    
    for fmt in FIGURE_FORMAT:
        save_path = output_dir / f'{base_name}.{fmt}'
        if fmt == 'svg':
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='none', edgecolor='none', transparent=True)
        else:
            fig.savefig(save_path, dpi=FIGURE_DPI, 
                       facecolor='white', edgecolor='none')
        print(f'    Saved (LFP vs Fiber): {save_path.name}')
    
    plt.close(fig)


def _format_colorbar_scientific(cbar, cax, fig, clim):
    """Helper function to format colorbar with scientific notation."""
    vmin, vmax = clim
    cbar_ticks = np.linspace(vmin, vmax, 5)
    cbar.set_ticks(cbar_ticks)
    
    max_abs = max(abs(vmin), abs(vmax))
    if max_abs > 0:
        exponent = int(np.floor(np.log10(max_abs)))
    else:
        exponent = 0
    scale_factor = 10 ** exponent
    
    scaled_ticks = cbar_ticks / scale_factor
    tick_labels = [f'{t:.1f}' for t in scaled_ticks]
    cbar.set_ticklabels(tick_labels)
    
    if exponent != 0:
        exp_label = f'×10$^{{{exponent}}}$'
        cbar_pos = cax.get_position()
        fig.text(cbar_pos.x0, cbar_pos.y1 + 0.02, exp_label,
                 fontsize=8, ha='left', va='bottom')


def process_animal(animal_path):
    """
    Process all running/non-running spectrograms for a single animal.
    
    For each colormap mode (per-signal vs viridis), generates:
    1. Single spectrogram figures (simple + full) per signal/state
    2. Running vs Rest comparison per LFP and Fiber1
    3. LFP vs Fiber comparison per state (running, rest, all_epochs)
    Writes to FIGURES_DIR and FIGURES_DIR_VIRIDIS (mirrored folder layout).
    
    Parameters
    ----------
    animal_path : Path
        Path to the animal's run_nonrun folder
    """
    animal_name = animal_path.parent.name
    condition_name = animal_path.parent.parent.name
    
    print(f'\n{"="*60}')
    print(f'Processing animal: {animal_name}')
    print(f'Condition: {condition_name}')
    print(f'{"="*60}')
    
    # Output trees: per-signal colormaps vs canonical viridis (same folder layout)
    export_modes = [
        ("per_signal", FIGURES_DIR / condition_name / animal_name, "LFP/Fiber colormaps"),
        ("viridis", FIGURES_DIR_VIRIDIS / condition_name / animal_name, "viridis (canonical)"),
    ]
    
    # =========================================================================
    # STEP 1: Load all data first
    # =========================================================================
    # Structure: all_data[signal_type][state] = data_dict
    all_data = {}
    
    # Find all signal folders (LFP, Fiber1, Fiber2, etc.)
    signal_folders = [f for f in animal_path.iterdir() if f.is_dir()]
    
    for signal_folder in sorted(signal_folders):
        signal_type = signal_folder.name
        all_data[signal_type] = {}
        
        for state in ['running', 'non_running', 'all_epochs']:
            state_folder = signal_folder / state
            mat_file = state_folder / f'PhaseAlignedSpectrogram_{state}.mat'
            
            if not mat_file.exists():
                continue
            
            try:
                data = load_spectrogram_data(mat_file)
                all_data[signal_type][state] = data
                print(f'  Loaded: {signal_type}/{state} ({data["n_epochs"]} epochs)')
            except Exception as e:
                print(f'  Error loading {signal_type}/{state}: {e}')
    
    for cmap_mode, output_base, mode_label in export_modes:
        print(f'\n{"="*60}')
        print(f'Output: {mode_label}')
        print(f'Root: {output_base}')
        print(f'{"="*60}')
        
        # =========================================================================
        # STEP 2: Generate single spectrogram figures
        # =========================================================================
        print(f'\n--- Generating single spectrogram figures ({mode_label}) ---')
        
        for signal_type in all_data:
            for state in all_data[signal_type]:
                data = all_data[signal_type][state]
                output_dir = output_base / signal_type / state
                
                print(f'  [{signal_type}/{state}] Creating figures...')
                
                try:
                    create_spectrogram_figure_simple(
                        data, signal_type, state, animal_name, condition_name, output_dir,
                        cmap_mode=cmap_mode,
                    )
                except Exception as e:
                    print(f'    Error (simple): {e}')
                
                try:
                    create_spectrogram_figure_full(
                        data, signal_type, state, animal_name,
                        condition_name, output_dir, cmap_mode=cmap_mode,
                    )
                except Exception as e:
                    print(f'    Error (full): {e}')
        
        # =========================================================================
        # STEP 3: Running vs Rest comparison (LFP, Fiber1)
        # =========================================================================
        print(f'\n--- Generating Running vs Rest comparison ({mode_label}) ---')
        
        for signal_type in ['LFP', 'Fiber1']:
            if signal_type not in all_data:
                print(f'  [{signal_type}] Not available, skipping...')
                continue
            
            if 'running' not in all_data[signal_type] or 'non_running' not in all_data[signal_type]:
                print(f'  [{signal_type}] Missing running or rest data, skipping...')
                continue
            
            print(f'  [{signal_type}] Creating Running vs Rest comparison...')
            try:
                comparison_output_dir = output_base / signal_type
                create_spectrogram_comparison_figure(
                    all_data[signal_type]['running'],
                    all_data[signal_type]['non_running'],
                    signal_type,
                    animal_name,
                    condition_name,
                    comparison_output_dir,
                    cmap_mode=cmap_mode,
                )
            except Exception as e:
                print(f'    Error: {e}')
        
        # =========================================================================
        # STEP 4: LFP vs Fiber comparison
        # =========================================================================
        print(f'\n--- Generating LFP vs Fiber comparison ({mode_label}) ---')
        
        if 'LFP' not in all_data or 'Fiber1' not in all_data:
            print(f'  Missing LFP or Fiber1 data, skipping LFP vs Fiber comparisons...')
        else:
            for state in ['running', 'non_running', 'all_epochs']:
                state_map = {'running': 'Running', 'non_running': 'Rest', 'all_epochs': 'All Epochs'}
                state_label = state_map.get(state, state)
                
                if state not in all_data['LFP'] or state not in all_data['Fiber1']:
                    print(f'  [{state_label}] Missing LFP or Fiber1 data, skipping...')
                    continue
                
                print(f'  [{state_label}] Creating LFP vs Fiber comparison...')
                try:
                    comparison_output_dir = output_base / 'comparison'
                    create_lfp_fiber_comparison_figure(
                        all_data['LFP'][state],
                        all_data['Fiber1'][state],
                        state,
                        animal_name,
                        condition_name,
                        comparison_output_dir,
                        cmap_mode=cmap_mode,
                    )
                except Exception as e:
                    print(f'    Error: {e}')

    # =========================================================================
    # Write manifest file listing source trials and generation timestamp
    # =========================================================================
    _write_manifest(animal_path, all_data, animal_name, condition_name,
                    [FIGURES_DIR / condition_name / animal_name,
                     FIGURES_DIR_VIRIDIS / condition_name / animal_name])


def _write_manifest(animal_run_nonrun_path, all_data, animal_name, condition_name,
                    output_roots):
    """Write a manifest.txt in each output root with trials, parameters, and timestamps."""
    process_animal_dir = animal_run_nonrun_path.parent  # …/<AnimalID>
    trial_dirs = sorted([
        d.name for d in process_animal_dir.iterdir()
        if d.is_dir() and d.name != 'run_nonrun'
    ])

    # --- Collect per-state epoch info from all_data ---
    state_info = {}
    for sig_name, states in all_data.items():
        for state, d in states.items():
            if state not in state_info:
                state_info[state] = {}
            state_info[state][sig_name] = d.get('n_epochs', '?')
    first_data = None
    for sig in all_data.values():
        for d in sig.values():
            first_data = d
            break
        if first_data:
            break

    # --- Read run_nonrun .mat for velocity thresholds ---
    run_thresh, nonrun_thresh, lfp_field_str = '?', '?', '?'
    freq_range_str, n_freq_bins, n_phase_bins_str = '?', '?', '?'
    for mat_candidate in animal_run_nonrun_path.rglob('PhaseAlignedSpectrogram_*.mat'):
        try:
            md = loadmat(str(mat_candidate), squeeze_me=True)
            if 'running_threshold' in md:
                run_thresh = md['running_threshold']
            if 'non_running_threshold' in md:
                nonrun_thresh = md['non_running_threshold']
            if 'lfp_ephys_field' in md:
                lfp_field_str = _decode_mat_char(md['lfp_ephys_field'])
            fv = np.asarray(md.get('freq_vector', [])).flatten()
            if fv.size:
                freq_range_str = f"{fv[0]:.1f} – {fv[-1]:.1f} Hz"
                n_freq_bins = str(fv.size)
            pb = np.asarray(md.get('phase_bins', md.get('phase_bins_deg', []))).flatten()
            if pb.size:
                n_phase_bins_str = str(pb.size)
            break
        except Exception:
            pass

    # --- Read Results metadata for full analysis parameters ---
    analysis_params = {}
    for trial_dir in trial_dirs[:1]:
        results_glob = list((process_animal_dir / trial_dir).rglob(
            'PhaseAlignedSpectrogram_Results_*.mat'))
        if not results_glob:
            continue
        try:
            import mat73 as _m73
            res = _m73.loadmat(str(results_glob[0]))
            meta = (res.get('results') or res).get('metadata', {})
            if meta:
                analysis_params = {
                    'low_freq_band': meta.get('low_freq_band'),
                    'high_freq_range': meta.get('high_freq_range'),
                    'phase_method': meta.get('phase_method'),
                    'phase_source': meta.get('phase_source'),
                    'wavelet_cycles': meta.get('wavelet_cycles'),
                    'cycles_per_epoch': meta.get('cycles_per_epoch'),
                    'sampling_rate': meta.get('sampling_rate'),
                    'lfp_ephys_field': meta.get('lfp_ephys_field'),
                }
                if lfp_field_str == '?' and meta.get('lfp_ephys_field'):
                    lfp_field_str = str(meta['lfp_ephys_field'])
        except Exception:
            pass

    n_trials = len(trial_dirs) if trial_dirs else '?'

    lines = [
        "Phase-Aligned Spectrogram Figures — Generation Manifest",
        "========================================================",
        f"Animal:     {animal_name}",
        f"Condition:  {condition_name}",
        f"Generated:  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        f"Source data: Process/phase_aligned_spectrogram/{condition_name}/{animal_name}/run_nonrun/",
        "",
        f"Trials included ({n_trials}):",
    ]
    if trial_dirs:
        for t in trial_dirs:
            lines.append(f"  - {t}")
    else:
        lines.append("  (no per-trial folders found; data may be pre-merged)")

    lines += [
        "",
        "=== Analysis Parameters ===",
        f"LFP source field:          {lfp_field_str}",
        f"Sampling rate:             {analysis_params.get('sampling_rate', '?')} Hz",
        "",
        "Phase extraction:",
        f"  Low-frequency band:      {analysis_params.get('low_freq_band', '?')} Hz",
        f"  Method:                  {analysis_params.get('phase_method', '?')}",
        f"  Phase source:            {analysis_params.get('phase_source', '?')} "
        f"({'each signal uses LFP phase' if analysis_params.get('phase_source') == 'lfp' else 'each signal uses its own phase'})",
        "",
        "Wavelet spectrogram:",
        f"  Frequency range:         {freq_range_str}",
        f"  Number of freq bins:     {n_freq_bins}",
        f"  Morlet wavelet cycles:   {analysis_params.get('wavelet_cycles', '?')}",
        "",
        "Epoch & phase binning:",
        f"  Cycles per epoch:        {analysis_params.get('cycles_per_epoch', '?')}",
        f"  Number of phase bins:    {n_phase_bins_str}",
        "",
        "Velocity classification (strict per-sample):",
        f"  Running threshold:       > {run_thresh} cm/s (ALL samples in epoch)",
        f"  Rest threshold:          < {nonrun_thresh} cm/s (ALL samples in epoch)",
        f"  Intermediate epochs excluded from running/rest averages",
        "",
        "=== Epoch Counts Per State ===",
    ]
    for state in sorted(state_info):
        for sig, ne in sorted(state_info[state].items()):
            lines.append(f"  {state:15s} / {sig:8s}: {ne} epochs")

    lines += [
        "",
        f"Signals plotted: {sorted(all_data.keys())}",
        f"States plotted:  {sorted({s for sig in all_data.values() for s in sig})}",
        "",
        "Display settings:",
        f"  Freq display range:      {FREQ_DISPLAY_MIN} – {FREQ_DISPLAY_MAX} Hz",
        f"  Nonlinear Y-axis:        {'Yes' if FREQ_STRETCH else 'No'}",
    ]
    if FREQ_STRETCH:
        lines.append(f"  Stretch break:           {FREQ_STRETCH_BREAK} Hz")
        lines.append(f"  Stretch ratio:           {FREQ_STRETCH_RATIO}")
    lines += [
        "",
        f"Script: Code/plot_phase_aligned_spectrogram.py",
    ]

    content = "\n".join(lines) + "\n"
    for root in output_roots:
        root.mkdir(parents=True, exist_ok=True)
        manifest_path = root / "manifest.txt"
        manifest_path.write_text(content, encoding="utf-8")
        print(f"  Manifest written: {manifest_path}")


def find_all_animals(base_dir):
    """
    Find all animals with run_nonrun data across all conditions.
    
    Returns
    -------
    animal_paths : list
        List of Path objects pointing to run_nonrun folders
    """
    base_dir = Path(base_dir)
    animal_paths = []
    
    # Scan for conditions (e.g., contraHPfiber_contraHPLFP, ipsiHPfiber_contraHPLFP)
    for condition_folder in base_dir.iterdir():
        if not condition_folder.is_dir():
            continue
        
        # Scan for animals within each condition
        for animal_folder in condition_folder.iterdir():
            if not animal_folder.is_dir():
                continue
            
            # Check if this animal has run_nonrun data
            run_nonrun_path = animal_folder / 'run_nonrun'
            if run_nonrun_path.exists() and run_nonrun_path.is_dir():
                animal_paths.append(run_nonrun_path)
    
    return sorted(animal_paths)


def main():
    """Main entry point."""
    print("="*70)
    print("Phase-Aligned Spectrogram Plotting Script")
    print("="*70)
    
    # Print phase definition
    print("\n" + "-"*70)
    print("PHASE DEFINITION:")
    print("-"*70)
    print("Phase is extracted from 5-9 Hz band using Morlet wavelet (center frequency 7 Hz).")
    print("")
    print("  0° = RISING ZERO-CROSSING (transition from negative to positive)")
    print("  90° = PEAK of the low-frequency oscillation")
    print("  180° = FALLING ZERO-CROSSING (transition from positive to negative)")
    print("  270° = TROUGH of the low-frequency oscillation")
    print("")
    print("The 2-cycle view shows phases from 0° to 720° (two complete cycles).")
    print("The white dashed line at 360° marks the boundary between cycle 1 and cycle 2.")
    print("-"*70)
    
    # Check if Process directory exists
    if not PROCESS_DIR.exists():
        print(f"\nError: Process directory not found: {PROCESS_DIR}")
        print("Please ensure the MATLAB analysis has been run first.")
        return
    
    # Find all animals
    animal_paths = find_all_animals(PROCESS_DIR)
    
    if not animal_paths:
        print(f"\nNo animals with run_nonrun data found in: {PROCESS_DIR}")
        return
    
    print(f"\nFound {len(animal_paths)} animals with run_nonrun data:")
    for path in animal_paths:
        condition = path.parent.parent.name
        animal = path.parent.name
        print(f"  - {condition}/{animal}")
    
    # Process each animal
    for animal_path in animal_paths:
        process_animal(animal_path)
    
    print("\n" + "="*70)
    print("Processing complete!")
    print(f"Figures (per-signal colormaps): {FIGURES_DIR}")
    print(f"Figures (viridis only):         {FIGURES_DIR_VIRIDIS}")
    print("="*70)


if __name__ == '__main__':
    main()

