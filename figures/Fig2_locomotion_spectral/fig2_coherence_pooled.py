"""
================================================================================
FIGURE 2: Pooled Coherence & PSD Analysis - Python Plotting Script
================================================================================

This script creates publication-quality figures for:
  - SESSION-POOLED data (trials pooled within a recording session)
  - ANIMAL-POOLED data (all sessions pooled for each animal)

USAGE:
------
  Option 1: Run via master pipeline (recommended)
      python run_all_plots.py

  Option 2: Run standalone
      python fig2_coherence_pooled.py

================================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from pathlib import Path
import warnings

warnings.filterwarnings('ignore', category=UserWarning)

# common.py / plotting_config.py live in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# ==============================================================================
#  CONFIGURATION - Import from central config or use defaults
# ==============================================================================

try:
    from plotting_config import (
        BEHAVIOR_MODE, FIGURE_DPI, FIGURE_FORMATS, FREQ_MIN, FREQ_MAX,
        METHODS, ANIMALS, get_animals_to_process, get_session_ids,
        get_session_pooled_input_dir, get_session_pooled_output_dir,
        get_animal_pooled_input_dir, get_animal_pooled_output_dir,
        get_base_dir, FONT_SIZE_TITLE, FONT_SIZE_SUPTITLE, FONT_SIZE_LABEL,
        FONT_SIZE_TICK, FONT_SIZE_LEGEND, FONT_SIZE_STATS, FONT_SIZE_BAND,
        AXIS_LINEWIDTH, TICK_WIDTH, TICK_LENGTH, LINE_WIDTH_TRACE,
        LINE_WIDTH_DASHED, LINE_WIDTH_BAND, COLOR_REST, COLOR_RUN, COLOR_OVERALL,
        COLOR_GEVI, COLOR_LFP, BAND_LINE_COLORS, FREQ_BANDS, SEM_ALPHA,
        COLOR_COH_REST, COLOR_COH_RUN, COLOR_COH_OVERALL,
    )
    USING_CENTRAL_CONFIG = True
    BASE_DIR = get_base_dir()
except ImportError:
    USING_CENTRAL_CONFIG = False
    # Fallback defaults (only used if plotting_config.py can't be found)
    BASE_OUTPUT_DIR = Path(__file__).resolve().parents[2] / "Figures" / "Spectral_data_outputs"
    BEHAVIOR_MODE = 'clear'
    BASE_DIR = BASE_OUTPUT_DIR / BEHAVIOR_MODE
    FIGURE_DPI = 300
    FIGURE_FORMATS = ['png', 'pdf', 'svg']
    FREQ_MIN = 2
    FREQ_MAX = 70
    METHODS = ['mscohere', 'fieldtrip']
    
    # Fallback animal list -- EDIT THIS to your own cohort
    ANIMALS = [
        {'mouse_id': 'Animal01', 'sessions': [
            {'session_id': '01_09_25-R1'}, {'session_id': '02_09_25-R1'},
            {'session_id': '03_09_25-R1'}, {'session_id': '03_09_25-R2'},
        ]},
        {'mouse_id': 'Animal02', 'sessions': [
            {'session_id': '01_01_26-R1'}, {'session_id': '01_01_26-R3'},
        ]},
    ]
    
    def get_animals_to_process():
        return ANIMALS
    
    def get_session_ids(animal):
        return [s['session_id'] for s in animal['sessions']]

# ==============================================================================
#  PUBLICATION FONT SIZES - MATCHING SINGLE-TRIAL SCRIPT
# ==============================================================================
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 20
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_STATS = 13
FONT_SIZE_BAND = 13  # For frequency band labels

# Axis line widths
AXIS_LINEWIDTH = 2.0
TICK_WIDTH = 1.8
TICK_LENGTH = 7

# Plot line widths
LINE_WIDTH_TRACE = 2.5
LINE_WIDTH_DASHED = 2.5
LINE_WIDTH_BAND = 1.0  # For frequency band dotted lines

# ==============================================================================
#  COLOR DEFINITIONS - MATCHING SINGLE-TRIAL SCRIPT
# ==============================================================================
# Teal shades for rest vs run (PSD plots)
COLOR_REST = np.array([0.05, 0.35, 0.45])      # Darker teal for rest
COLOR_RUN = np.array([0.25, 0.65, 0.65])       # Lighter teal for run
COLOR_OVERALL = np.array([0.08, 0.45, 0.52])   # Mid teal for overall

# Coherence colors (Indigo/blue-violet blend - represents LFP-GEVI coupling)
COLOR_COH_REST = np.array([0.20, 0.25, 0.50])     # Dark indigo
COLOR_COH_RUN = np.array([0.45, 0.50, 0.70])      # Periwinkle/light indigo
COLOR_COH_OVERALL = np.array([0.30, 0.35, 0.58])  # Mid indigo for overall

# Colors for LFP and GEVI traces (consistent with Figure 1)
COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])  # Teal (from viridis)
COLOR_LFP = np.array([0.35, 0.25, 0.45])               # Purple-grey

# Frequency band line colors (specific colors for each band)
BAND_LINE_COLORS = {
    'theta': (0.4, 0.2, 0.6, 0.6),    # Purple
    'alpha': (0.2, 0.5, 0.5, 0.6),    # Teal
    'beta': (0.3, 0.6, 0.3, 0.6),     # Green
    'gamma': (0.7, 0.5, 0.2, 0.6),    # Orange/brown
}

# Frequency band boundaries (Hz)
FREQ_BANDS = {
    'theta': (4, 8, 'θ'),
    'alpha': (8, 12, 'α'),
    'beta': (12, 30, 'β'),
    'gamma': (30, 70, 'γ'),
}


# ==============================================================================
#  WINDOWS LONG PATH SUPPORT
# ==============================================================================

from common import load_matlab_struct as load_matlab_data, to_long_path  # shared helpers (were local copies)


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================


def get_field(data, field_name, default=None):
    """Safely get field from loaded MATLAB data."""
    if data is None:
        return default
    if field_name in data:
        val = data[field_name]
        if isinstance(val, np.ndarray):
            return val.flatten() if val.ndim == 1 or (val.ndim == 2 and 1 in val.shape) else val
        return val
    return default


def style_axis_publication(ax, remove_top_right=True):
    """
    Apply publication-quality styling to an axis.
    Matches the single-trial script styling.
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


def add_frequency_band_lines(ax, y_max=None, add_labels=True, labels_inside=True):
    """
    Add colored dotted vertical lines at frequency band boundaries.
    Each band has its own color scheme for visibility.
    Matches the single-trial script implementation.
    
    Parameters
    ----------
    ax : matplotlib axis
    y_max : float or None
        Y position for labels. If None, gets from axis ylim.
    add_labels : bool
        Whether to add frequency band labels (θ, α, β, γ)
    labels_inside : bool
        If True, place labels inside plot at top. If False, place above plot using transform.
    """
    # Band info: (start_freq, end_freq, label, color_key)
    band_info = [
        (4, 8, 'θ', 'theta'),
        (8, 12, 'α', 'alpha'),
        (12, 30, 'β', 'beta'),
        (30, 70, 'γ', 'gamma'),
    ]
    
    band_centers = [6, 10, 21, 45]  # Approximate centers for labels (gamma shifted left)
    
    # Get y_max from axis if not provided
    if y_max is None:
        _, y_max = ax.get_ylim()
    
    # Draw vertical lines at band boundaries with specific colors
    for (f_start, f_end, label, color_key), center in zip(band_info, band_centers):
        color = BAND_LINE_COLORS[color_key]
        
        # Draw start line
        if f_start >= FREQ_MIN and f_start <= FREQ_MAX:
            ax.axvline(x=f_start, color=color[:3], alpha=color[3],
                      linestyle=':', linewidth=LINE_WIDTH_BAND * 1.5)
        
        # Draw end line (if not already drawn by next band)
        if f_end >= FREQ_MIN and f_end <= FREQ_MAX and f_end == 70:
            ax.axvline(x=f_end, color=color[:3], alpha=color[3],
                      linestyle=':', linewidth=LINE_WIDTH_BAND * 1.5)
        
        # Add band label with matching color
        if add_labels and center <= FREQ_MAX and center >= FREQ_MIN:
            if labels_inside:
                # Place inside plot at top (for coherence plots with fixed 0-1 range)
                ax.text(center, y_max * 0.97, label, ha='center', va='top',
                       fontsize=FONT_SIZE_BAND, fontweight='bold', 
                       color=color[:3])
            else:
                # Place ABOVE the plot using axes transform (for PSD with variable range)
                x_data_to_axes = (center - FREQ_MIN) / (FREQ_MAX - FREQ_MIN)
                ax.text(x_data_to_axes, 1.02, label, ha='center', va='bottom',
                       fontsize=FONT_SIZE_BAND, fontweight='bold', 
                       color=color[:3], transform=ax.transAxes)


# ==============================================================================
#  PLOTTING FUNCTIONS
# ==============================================================================

def plot_coherence_spectrum(data, method, title_prefix, output_path):
    """
    Plot coherence spectrum with Overall and Rest vs Run panels.
    Publication-quality styling matching single-trial script.
    
    Creates a 1x2 figure:
      - Left: Overall coherence
      - Right: Rest vs Run coherence
    """
    freq = get_field(data, 'freq')
    # Check for both possible field names (pooledtrials uses coh_spectrum, animalpooled uses coh_overall)
    coh_overall = get_field(data, 'coh_overall')
    if coh_overall is None:
        coh_overall = get_field(data, 'coh_spectrum')  # Fallback for session-pooled data
    coh_rest = get_field(data, 'coh_rest')
    coh_run = get_field(data, 'coh_run')
    pct_rest = get_field(data, 'pct_rest', 0)
    pct_run = get_field(data, 'pct_run', 0)
    
    if freq is None or coh_overall is None:
        print(f"    Missing coherence data (no 'freq', 'coh_overall', or 'coh_spectrum'), skipping plot")
        return False
    
    # Check dimension match - critical for plotting
    if len(freq) != len(coh_overall):
        print(f"    Dimension mismatch: freq={len(freq)}, coh_overall={len(coh_overall)} - skipping plot")
        print(f"    (This may indicate an issue with the MATLAB data generation)")
        return False
    
    # Create figure with publication sizing
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.subplots_adjust(left=0.08, right=0.92, top=0.85, bottom=0.12, wspace=0.30)
    
    # === Left Panel: Overall Coherence ===
    ax1 = axes[0]
    
    # Add frequency band lines first (behind the data)
    add_frequency_band_lines(ax1, y_max=1.0, add_labels=True, labels_inside=True)
    
    # Plot overall coherence as LINE plot (no fill) - matching single-trial style
    ax1.plot(freq, coh_overall, '-', color=COLOR_COH_OVERALL, linewidth=LINE_WIDTH_TRACE)
    
    ax1.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax1.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
    ax1.set_xlim([FREQ_MIN, FREQ_MAX])
    ax1.set_ylim([0, 1])
    
    # Apply publication styling
    style_axis_publication(ax1)
    
    # Add statistics above plot (matching single-trial style)
    peak_idx = np.argmax(coh_overall)
    peak_coh = coh_overall[peak_idx]
    peak_freq = freq[peak_idx]
    theta_mask = (freq >= 4) & (freq <= 8)
    mean_theta = np.mean(coh_overall[theta_mask]) if np.any(theta_mask) else 0
    
    stats_text = f'Peak: {peak_coh:.2f} @ {peak_freq:.1f} Hz | θ-band: {mean_theta:.2f}'
    ax1.text(0.5, 1.06, stats_text, transform=ax1.transAxes,
            ha='center', va='bottom', fontsize=FONT_SIZE_STATS,
            fontweight='bold', color=COLOR_COH_OVERALL)
    
    # === Right Panel: Rest vs Run ===
    ax2 = axes[1]
    
    # Add frequency band lines first
    add_frequency_band_lines(ax2, y_max=1.0, add_labels=True, labels_inside=True)
    
    has_rest_run = (coh_rest is not None and coh_run is not None and 
                    len(coh_rest) > 0 and len(coh_run) > 0)
    
    legend_handles = []
    legend_labels = []
    stats_lines = []
    
    if has_rest_run:
        # Plot REST coherence (solid line - dark indigo)
        h_rest, = ax2.plot(freq, coh_rest, '-', color=COLOR_COH_REST, linewidth=LINE_WIDTH_TRACE)
        legend_handles.append(h_rest)
        legend_labels.append(f'Rest ({pct_rest:.0f}%)')
        
        # Get peak for rest
        peak_idx_rest = np.argmax(coh_rest)
        peak_coh_rest = coh_rest[peak_idx_rest]
        peak_freq_rest = freq[peak_idx_rest]
        stats_lines.append(f'Rest pk: {peak_coh_rest:.2f} @ {peak_freq_rest:.1f} Hz')
        
        # Plot RUN coherence (dashed line - light indigo)
        h_run, = ax2.plot(freq, coh_run, '--', color=COLOR_COH_RUN, linewidth=LINE_WIDTH_DASHED)
        legend_handles.append(h_run)
        legend_labels.append(f'Run ({pct_run:.0f}%)')
        
        # Get peak for run
        peak_idx_run = np.argmax(coh_run)
        peak_coh_run = coh_run[peak_idx_run]
        peak_freq_run = freq[peak_idx_run]
        stats_lines.append(f'Run pk: {peak_coh_run:.2f} @ {peak_freq_run:.1f} Hz')
        
        # Legend inside plot (upper right)
        ax2.legend(legend_handles, legend_labels, loc='upper right',
                  fontsize=FONT_SIZE_LEGEND, frameon=False)
    else:
        ax2.text(0.5, 0.5, 'Insufficient Rest/Run data', 
                transform=ax2.transAxes, ha='center', va='center',
                fontsize=FONT_SIZE_LABEL, color='gray')
    
    ax2.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax2.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
    ax2.set_xlim([FREQ_MIN, FREQ_MAX])
    ax2.set_ylim([0, 1])
    
    # Apply publication styling
    style_axis_publication(ax2)
    
    # Add peak stats above plot
    if stats_lines:
        stats_text = ' | '.join(stats_lines)
        ax2.text(0.5, 1.06, stats_text, transform=ax2.transAxes,
                ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                fontweight='bold', color=COLOR_COH_OVERALL)
    
    # Main title
    fig.suptitle(f'{title_prefix}\nCoherence Spectrum ({method.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)
    
    # Ensure output directory exists (with long path support)
    output_dir_long = to_long_path(str(output_path.parent))
    os.makedirs(output_dir_long, exist_ok=True)
    
    # Save figure
    for fmt in FIGURE_FORMATS:
        save_path = str(output_path.with_suffix(f'.{fmt}'))
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    
    plt.close(fig)
    print(f"    Saved: {output_path.stem}")
    return True


def plot_coherence_spectrum_logfreq(data, method, title_prefix, output_path):
    """
    Plot coherence spectrum with log-frequency x-axis.
    Same content as plot_coherence_spectrum, with vertical band lines retained.
    """
    freq = get_field(data, 'freq')
    coh_overall = get_field(data, 'coh_overall')
    if coh_overall is None:
        coh_overall = get_field(data, 'coh_spectrum')
    coh_rest = get_field(data, 'coh_rest')
    coh_run = get_field(data, 'coh_run')
    pct_rest = get_field(data, 'pct_rest', 0)
    pct_run = get_field(data, 'pct_run', 0)

    if freq is None or coh_overall is None:
        print(f"    Missing coherence data (no 'freq', 'coh_overall', or 'coh_spectrum'), skipping logfreq plot")
        return False
    if len(freq) != len(coh_overall):
        print(f"    Dimension mismatch: freq={len(freq)}, coh_overall={len(coh_overall)} - skipping logfreq plot")
        return False

    freq = np.asarray(freq, dtype=float).flatten()
    pos_mask = freq > 0
    if np.sum(pos_mask) < 2:
        print("    Not enough positive frequencies for log axis, skipping logfreq plot")
        return False

    f = freq[pos_mask]
    c_all = np.asarray(coh_overall).flatten()[pos_mask]
    c_rest = np.asarray(coh_rest).flatten()[pos_mask] if coh_rest is not None else None
    c_run = np.asarray(coh_run).flatten()[pos_mask] if coh_run is not None else None

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    fig.subplots_adjust(left=0.08, right=0.92, top=0.85, bottom=0.12, wspace=0.30)

    # Left: overall
    ax1 = axes[0]
    add_frequency_band_lines(ax1, y_max=1.0, add_labels=True, labels_inside=True)
    ax1.plot(f, c_all, '-', color=COLOR_COH_OVERALL, linewidth=LINE_WIDTH_TRACE)
    ax1.set_xscale('log')
    ax1.set_xlabel('Frequency (Hz, log scale)', fontsize=FONT_SIZE_LABEL)
    ax1.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
    ax1.set_xlim([max(FREQ_MIN, 1), FREQ_MAX])
    ax1.set_ylim([0, 1])
    style_axis_publication(ax1)

    peak_idx = np.argmax(c_all)
    peak_coh = c_all[peak_idx]
    peak_freq = f[peak_idx]
    theta_mask = (f >= 4) & (f <= 8)
    mean_theta = np.mean(c_all[theta_mask]) if np.any(theta_mask) else 0
    stats_text = f'Peak: {peak_coh:.2f} @ {peak_freq:.1f} Hz | θ-band: {mean_theta:.2f}'
    ax1.text(0.5, 1.06, stats_text, transform=ax1.transAxes,
             ha='center', va='bottom', fontsize=FONT_SIZE_STATS,
             fontweight='bold', color=COLOR_COH_OVERALL)

    # Right: rest/run
    ax2 = axes[1]
    add_frequency_band_lines(ax2, y_max=1.0, add_labels=True, labels_inside=True)
    has_rest_run = (c_rest is not None and c_run is not None and len(c_rest) > 0 and len(c_run) > 0)

    stats_lines = []
    if has_rest_run:
        h_rest, = ax2.plot(f, c_rest, '-', color=COLOR_COH_REST, linewidth=LINE_WIDTH_TRACE)
        h_run, = ax2.plot(f, c_run, '--', color=COLOR_COH_RUN, linewidth=LINE_WIDTH_DASHED)
        ax2.legend([h_rest, h_run], [f'Rest ({pct_rest:.0f}%)', f'Run ({pct_run:.0f}%)'],
                   loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
        i_rest = np.argmax(c_rest)
        i_run = np.argmax(c_run)
        stats_lines.append(f'Rest pk: {c_rest[i_rest]:.2f} @ {f[i_rest]:.1f} Hz')
        stats_lines.append(f'Run pk: {c_run[i_run]:.2f} @ {f[i_run]:.1f} Hz')
    else:
        ax2.text(0.5, 0.5, 'Insufficient Rest/Run data', transform=ax2.transAxes,
                 ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')

    ax2.set_xscale('log')
    ax2.set_xlabel('Frequency (Hz, log scale)', fontsize=FONT_SIZE_LABEL)
    ax2.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
    ax2.set_xlim([max(FREQ_MIN, 1), FREQ_MAX])
    ax2.set_ylim([0, 1])
    style_axis_publication(ax2)

    if stats_lines:
        ax2.text(0.5, 1.06, ' | '.join(stats_lines), transform=ax2.transAxes,
                 ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                 fontweight='bold', color=COLOR_COH_OVERALL)

    fig.suptitle(f'{title_prefix}\nCoherence Spectrum (log frequency, {method.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)

    output_dir_long = to_long_path(str(output_path.parent))
    os.makedirs(output_dir_long, exist_ok=True)
    for fmt in FIGURE_FORMATS:
        save_path = str(output_path.with_suffix(f'.{fmt}'))
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.stem}")
    return True


def plot_psd_spectrum(data, method, title_prefix, output_path):
    """
    Plot PSD spectrum for LFP and GEVI with Overall and Rest vs Run panels.
    Publication-quality styling matching single-trial script.
    
    Creates a 2x2 figure:
      - Top row: LFP PSD (Overall, Rest vs Run)
      - Bottom row: GEVI PSD (Overall, Rest vs Run)
    """
    freq = get_field(data, 'psd_freq')
    if freq is None:
        freq = get_field(data, 'freq')  # Fallback
    
    psd_lfp = get_field(data, 'psd_lfp')
    psd_gevi = get_field(data, 'psd_gevi')
    psd_lfp_rest = get_field(data, 'psd_lfp_rest')
    psd_lfp_run = get_field(data, 'psd_lfp_run')
    psd_gevi_rest = get_field(data, 'psd_gevi_rest')
    psd_gevi_run = get_field(data, 'psd_gevi_run')
    pct_rest = get_field(data, 'pct_rest', 0)
    pct_run = get_field(data, 'pct_run', 0)
    
    if freq is None or psd_lfp is None:
        print(f"    Missing PSD data, skipping plot")
        return False
    
    # Check dimension match
    if len(freq) != len(psd_lfp):
        print(f"    PSD dimension mismatch: freq={len(freq)}, psd_lfp={len(psd_lfp)} - skipping plot")
        return False
    
    # Create figure with publication sizing
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    fig.subplots_adjust(left=0.08, right=0.92, top=0.90, bottom=0.08, 
                        hspace=0.35, wspace=0.30)
    
    # === Top Left: LFP Overall ===
    ax_lfp_overall = axes[0, 0]
    
    if len(psd_lfp) > 0:
        ax_lfp_overall.plot(freq, psd_lfp, '-', color=COLOR_LFP, 
                           linewidth=LINE_WIDTH_TRACE, label='LFP')
        
        # Stats text above plot
        peak_idx = np.argmax(psd_lfp)
        peak_val = psd_lfp[peak_idx]
        peak_freq = freq[peak_idx]
        
        stats_text = f'Peak: {peak_val:.1f} dB @ {peak_freq:.1f} Hz'
        ax_lfp_overall.text(0.5, 1.08, stats_text, transform=ax_lfp_overall.transAxes,
                           ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                           fontweight='bold', color=COLOR_LFP)
    
    ax_lfp_overall.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax_lfp_overall.set_ylabel('LFP Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax_lfp_overall.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax_lfp_overall)
    ax_lfp_overall.yaxis.set_major_locator(MultipleLocator(5))
    ax_lfp_overall.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    add_frequency_band_lines(ax_lfp_overall, add_labels=True, labels_inside=False)
    
    # === Top Right: LFP Rest vs Run ===
    ax_lfp_rr = axes[0, 1]
    
    legend_handles_lfp = []
    legend_labels_lfp = []
    lfp_rest_peak = None
    lfp_run_peak = None
    
    has_lfp_rr = (psd_lfp_rest is not None and psd_lfp_run is not None and
                  len(psd_lfp_rest) > 0 and len(psd_lfp_run) > 0)
    
    if has_lfp_rr:
        # LFP Rest (solid)
        h_rest, = ax_lfp_rr.plot(freq, psd_lfp_rest, '-', color=COLOR_LFP, 
                                linewidth=LINE_WIDTH_TRACE, alpha=0.9)
        legend_handles_lfp.append(h_rest)
        legend_labels_lfp.append(f'Rest ({pct_rest:.0f}%)')
        peak_idx = np.argmax(psd_lfp_rest)
        lfp_rest_peak = (psd_lfp_rest[peak_idx], freq[peak_idx])
        
        # LFP Run (dashed)
        h_run, = ax_lfp_rr.plot(freq, psd_lfp_run, '--', color=COLOR_LFP, 
                               linewidth=LINE_WIDTH_DASHED, alpha=0.7)
        legend_handles_lfp.append(h_run)
        legend_labels_lfp.append(f'Run ({pct_run:.0f}%)')
        peak_idx = np.argmax(psd_lfp_run)
        lfp_run_peak = (psd_lfp_run[peak_idx], freq[peak_idx])
        
        ax_lfp_rr.legend(legend_handles_lfp, legend_labels_lfp, 
                        loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
    else:
        ax_lfp_rr.text(0.5, 0.5, 'Insufficient Rest/Run data', 
                      transform=ax_lfp_rr.transAxes, ha='center', va='center',
                      fontsize=FONT_SIZE_LABEL, color='gray')
    
    ax_lfp_rr.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax_lfp_rr.set_ylabel('LFP Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax_lfp_rr.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax_lfp_rr)
    ax_lfp_rr.yaxis.set_major_locator(MultipleLocator(5))
    ax_lfp_rr.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    add_frequency_band_lines(ax_lfp_rr, add_labels=True, labels_inside=False)
    
    # Peak stats above plot
    if lfp_rest_peak and lfp_run_peak:
        stats_text = f'Rest pk: {lfp_rest_peak[0]:.0f} dB | Run pk: {lfp_run_peak[0]:.0f} dB'
        ax_lfp_rr.text(0.5, 1.08, stats_text, transform=ax_lfp_rr.transAxes,
                      ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 2,
                      fontweight='bold', color=COLOR_LFP)
    
    # === Bottom Left: GEVI Overall ===
    ax_gevi_overall = axes[1, 0]
    
    if psd_gevi is not None and len(psd_gevi) > 0:
        ax_gevi_overall.plot(freq, psd_gevi, '-', color=COLOR_GEVI, 
                            linewidth=LINE_WIDTH_TRACE, label='GEVI')
        
        # Stats text above plot
        peak_idx = np.argmax(psd_gevi)
        peak_val = psd_gevi[peak_idx]
        peak_freq = freq[peak_idx]
        
        stats_text = f'Peak: {peak_val:.1f} dB @ {peak_freq:.1f} Hz'
        ax_gevi_overall.text(0.5, 1.08, stats_text, transform=ax_gevi_overall.transAxes,
                            ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                            fontweight='bold', color=COLOR_GEVI)
    else:
        ax_gevi_overall.text(0.5, 0.5, 'No GEVI PSD data', 
                            transform=ax_gevi_overall.transAxes, ha='center', va='center',
                            fontsize=FONT_SIZE_LABEL, color='gray')
    
    ax_gevi_overall.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax_gevi_overall.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax_gevi_overall.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax_gevi_overall)
    ax_gevi_overall.yaxis.set_major_locator(MultipleLocator(5))
    ax_gevi_overall.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    add_frequency_band_lines(ax_gevi_overall, add_labels=True, labels_inside=False)
    
    # === Bottom Right: GEVI Rest vs Run ===
    ax_gevi_rr = axes[1, 1]
    
    legend_handles_gevi = []
    legend_labels_gevi = []
    gevi_rest_peak = None
    gevi_run_peak = None
    
    has_gevi_rr = (psd_gevi_rest is not None and psd_gevi_run is not None and
                   len(psd_gevi_rest) > 0 and len(psd_gevi_run) > 0)
    
    if has_gevi_rr:
        # GEVI Rest (solid)
        h_rest, = ax_gevi_rr.plot(freq, psd_gevi_rest, '-', color=COLOR_GEVI, 
                                 linewidth=LINE_WIDTH_TRACE, alpha=0.9)
        legend_handles_gevi.append(h_rest)
        legend_labels_gevi.append(f'Rest ({pct_rest:.0f}%)')
        peak_idx = np.argmax(psd_gevi_rest)
        gevi_rest_peak = (psd_gevi_rest[peak_idx], freq[peak_idx])
        
        # GEVI Run (dashed)
        h_run, = ax_gevi_rr.plot(freq, psd_gevi_run, '--', color=COLOR_GEVI, 
                                linewidth=LINE_WIDTH_DASHED, alpha=0.7)
        legend_handles_gevi.append(h_run)
        legend_labels_gevi.append(f'Run ({pct_run:.0f}%)')
        peak_idx = np.argmax(psd_gevi_run)
        gevi_run_peak = (psd_gevi_run[peak_idx], freq[peak_idx])
        
        ax_gevi_rr.legend(legend_handles_gevi, legend_labels_gevi, 
                         loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
    else:
        ax_gevi_rr.text(0.5, 0.5, 'Insufficient Rest/Run data', 
                       transform=ax_gevi_rr.transAxes, ha='center', va='center',
                       fontsize=FONT_SIZE_LABEL, color='gray')
    
    ax_gevi_rr.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax_gevi_rr.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax_gevi_rr.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax_gevi_rr)
    ax_gevi_rr.yaxis.set_major_locator(MultipleLocator(5))
    ax_gevi_rr.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    add_frequency_band_lines(ax_gevi_rr, add_labels=True, labels_inside=False)
    
    # Peak stats above plot
    if gevi_rest_peak and gevi_run_peak:
        stats_text = f'Rest pk: {gevi_rest_peak[0]:.0f} dB | Run pk: {gevi_run_peak[0]:.0f} dB'
        ax_gevi_rr.text(0.5, 1.08, stats_text, transform=ax_gevi_rr.transAxes,
                       ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 2,
                       fontweight='bold', color=COLOR_GEVI)
    
    # Main title
    fig.suptitle(f'{title_prefix}\nPower Spectral Density – LFP and GEVI ({method.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)
    
    # Ensure output directory exists (with long path support)
    output_dir_long = to_long_path(str(output_path.parent))
    os.makedirs(output_dir_long, exist_ok=True)
    
    # Save figure
    for fmt in FIGURE_FORMATS:
        save_path = str(output_path.with_suffix(f'.{fmt}'))
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    
    plt.close(fig)
    print(f"    Saved: {output_path.stem}")
    return True


def plot_psd_spectrum_loglog(data, method, title_prefix, output_path):
    """
    Plot PSD spectrum on log-log scale for 1/f analysis.
    Log frequency axis shows 1/f slope as a straight line.
    """
    freq = get_field(data, 'psd_freq')
    if freq is None:
        freq = get_field(data, 'freq')
    
    psd_lfp = get_field(data, 'psd_lfp')
    psd_gevi = get_field(data, 'psd_gevi')
    psd_lfp_rest = get_field(data, 'psd_lfp_rest')
    psd_lfp_run = get_field(data, 'psd_lfp_run')
    psd_gevi_rest = get_field(data, 'psd_gevi_rest')
    psd_gevi_run = get_field(data, 'psd_gevi_run')
    pct_rest = get_field(data, 'pct_rest', 0)
    pct_run = get_field(data, 'pct_run', 0)
    
    if freq is None or psd_lfp is None:
        print(f"    Missing PSD data, skipping log-log plot")
        return False
    
    if len(freq) != len(psd_lfp):
        print(f"    PSD dimension mismatch - skipping log-log plot")
        return False
    
    # Filter to positive frequencies only (for log scale)
    freq_mask = freq > 0
    freq = freq[freq_mask]
    psd_lfp = psd_lfp[freq_mask]
    if psd_gevi is not None:
        psd_gevi = psd_gevi[freq_mask]
    if psd_lfp_rest is not None:
        psd_lfp_rest = psd_lfp_rest[freq_mask]
    if psd_lfp_run is not None:
        psd_lfp_run = psd_lfp_run[freq_mask]
    if psd_gevi_rest is not None:
        psd_gevi_rest = psd_gevi_rest[freq_mask]
    if psd_gevi_run is not None:
        psd_gevi_run = psd_gevi_run[freq_mask]
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 14))
    fig.subplots_adjust(left=0.08, right=0.92, top=0.90, bottom=0.08, 
                        hspace=0.35, wspace=0.30)
    
    # === Top Left: LFP Overall ===
    ax = axes[0, 0]
    if len(psd_lfp) > 0:
        ax.plot(freq, psd_lfp, '-', color=COLOR_LFP, linewidth=LINE_WIDTH_TRACE)
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('LFP Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    
    # === Top Right: LFP Rest vs Run ===
    ax = axes[0, 1]
    has_lfp_rr = (psd_lfp_rest is not None and psd_lfp_run is not None and
                  len(psd_lfp_rest) > 0 and len(psd_lfp_run) > 0)
    if has_lfp_rr:
        ax.plot(freq, psd_lfp_rest, '-', color=COLOR_LFP, linewidth=LINE_WIDTH_TRACE, 
                alpha=0.9, label=f'Rest ({pct_rest:.0f}%)')
        ax.plot(freq, psd_lfp_run, '--', color=COLOR_LFP, linewidth=LINE_WIDTH_DASHED, 
                alpha=0.7, label=f'Run ({pct_run:.0f}%)')
        ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('LFP Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    
    # === Bottom Left: GEVI Overall ===
    ax = axes[1, 0]
    if psd_gevi is not None and len(psd_gevi) > 0:
        ax.plot(freq, psd_gevi, '-', color=COLOR_GEVI, linewidth=LINE_WIDTH_TRACE)
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    
    # === Bottom Right: GEVI Rest vs Run ===
    ax = axes[1, 1]
    has_gevi_rr = (psd_gevi_rest is not None and psd_gevi_run is not None and
                   len(psd_gevi_rest) > 0 and len(psd_gevi_run) > 0)
    if has_gevi_rr:
        ax.plot(freq, psd_gevi_rest, '-', color=COLOR_GEVI, linewidth=LINE_WIDTH_TRACE, 
                alpha=0.9, label=f'Rest ({pct_rest:.0f}%)')
        ax.plot(freq, psd_gevi_run, '--', color=COLOR_GEVI, linewidth=LINE_WIDTH_DASHED, 
                alpha=0.7, label=f'Run ({pct_run:.0f}%)')
        ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
    ax.set_xscale('log')
    ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
    ax.set_xlim([FREQ_MIN, FREQ_MAX])
    style_axis_publication(ax)
    ax.yaxis.set_major_locator(MultipleLocator(5))
    
    fig.suptitle(f'{title_prefix}\nPower Spectral Density – Log-Log Scale ({method.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)
    
    output_dir_long = to_long_path(str(output_path.parent))
    os.makedirs(output_dir_long, exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = str(output_path.with_suffix(f'.{fmt}'))
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    
    plt.close(fig)
    print(f"    Saved: {output_path.stem}")
    return True


# ==============================================================================
#  MAIN PROCESSING
# ==============================================================================

def process_session_pooled(mouse_id, session_id, method):
    """Process session-pooled data for one session.
    
    New path structure:
        Input:  {BASE_DIR}/session_pooled/{mouse_id}/{session_id}/data/{method}.mat
        Output: {BASE_DIR}/session_pooled/{mouse_id}/{session_id}/figures/{method}_*.png
    """
    # Input file - new structure
    input_dir = BASE_DIR / 'session_pooled' / mouse_id / session_id / 'data'
    input_file = input_dir / f'{method}.mat'
    
    if not input_file.exists():
        print(f"    Input not found: {input_file}")
        return False
    
    # Output directory - sibling to data/
    output_dir = BASE_DIR / 'session_pooled' / mouse_id / session_id / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data = load_matlab_data(input_file)
    if data is None:
        return False
    
    title_prefix = f'{mouse_id} – {session_id} (Session Pooled)'
    
    # Plot coherence spectrum
    coh_output = output_dir / f'{method}_coherence'
    plot_coherence_spectrum(data, method, title_prefix, coh_output)
    
    # Plot PSD spectrum
    psd_output = output_dir / f'{method}_psd'
    plot_psd_spectrum(data, method, title_prefix, psd_output)
    
    return True


def process_animal_pooled(mouse_id, method):
    """Process animal-pooled data for one animal.
    
    New path structure:
        Input:  {BASE_DIR}/animal_pooled/{mouse_id}/data/{method}.mat
        Output: {BASE_DIR}/animal_pooled/{mouse_id}/figures/{method}_*.png
    """
    # Input file - new structure
    input_dir = BASE_DIR / 'animal_pooled' / mouse_id / 'data'
    input_file = input_dir / f'{method}.mat'
    
    if not input_file.exists():
        print(f"    Input not found: {input_file}")
        return False
    
    # Output directory - sibling to data/
    output_dir = BASE_DIR / 'animal_pooled' / mouse_id / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data = load_matlab_data(input_file)
    if data is None:
        return False
    
    title_prefix = f'{mouse_id} (Animal Pooled – All Sessions)'
    
    # Plot coherence spectrum
    coh_output = output_dir / f'{method}_coherence'
    plot_coherence_spectrum(data, method, title_prefix, coh_output)
    
    # Plot PSD spectrum
    psd_output = output_dir / f'{method}_psd'
    plot_psd_spectrum(data, method, title_prefix, psd_output)
    
    return True


def process_animal_concatenated(mouse_id, method):
    """Process animal-concatenated data for one animal.
    
    Animal-concatenated: All raw data from all sessions concatenated, spectra computed once.
    This differs from animal-pooled which averages spectra across sessions.
    
    Path structure:
        Input:  {BASE_DIR}/animal_concatenated/{mouse_id}/data/{method}.mat
        Output: {BASE_DIR}/animal_concatenated/{mouse_id}/figures/{method}_*.png
    """
    # Input file
    input_dir = BASE_DIR / 'animal_concatenated' / mouse_id / 'data'
    input_file = input_dir / f'{method}.mat'
    
    if not input_file.exists():
        print(f"    Input not found: {input_file}")
        return False
    
    # Output directory - sibling to data/
    output_dir = BASE_DIR / 'animal_concatenated' / mouse_id / 'figures'
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    data = load_matlab_data(input_file)
    if data is None:
        return False
    
    title_prefix = f'{mouse_id} (Animal Concatenated – All Sessions)'
    
    # Plot coherence spectrum
    coh_output = output_dir / f'{method}_coherence'
    plot_coherence_spectrum(data, method, title_prefix, coh_output)
    
    # Plot PSD spectrum
    psd_output = output_dir / f'{method}_psd'
    plot_psd_spectrum(data, method, title_prefix, psd_output)
    
    return True


def main():
    """Main function to process all animals and sessions."""
    print("=" * 80)
    print("POOLED COHERENCE & PSD PLOTTING (Publication Quality)")
    print("=" * 80)
    print(f"  Behavior Mode: {BEHAVIOR_MODE}")
    print(f"  Base Directory: {BASE_DIR}")
    print("=" * 80)
    
    total_session_plots = 0
    total_animal_plots = 0
    total_animal_concat_plots = 0
    
    # Get animals to process (from central config if available)
    animals_to_process = get_animals_to_process() if USING_CENTRAL_CONFIG else ANIMALS
    
    for animal in animals_to_process:
        mouse_id = animal['mouse_id']
        
        # Get session IDs (handle both old and new format)
        if USING_CENTRAL_CONFIG:
            sessions = get_session_ids(animal)
        else:
            # Old format compatibility
            if 'sessions' in animal and isinstance(animal['sessions'], list):
                if isinstance(animal['sessions'][0], dict):
                    sessions = [s['session_id'] for s in animal['sessions']]
                else:
                    sessions = animal['sessions']
            else:
                sessions = []
        
        print(f"\n{'='*60}")
        print(f"ANIMAL: {mouse_id}")
        print(f"{'='*60}")
        
        # Process session-pooled data
        print(f"\n  SESSION-POOLED DATA:")
        for session_id in sessions:
            print(f"\n  Session: {session_id}")
            for method in METHODS:
                print(f"    Method: {method}")
                if process_session_pooled(mouse_id, session_id, method):
                    total_session_plots += 1
        
        # Process animal-pooled data
        print(f"\n  ANIMAL-POOLED DATA:")
        for method in METHODS:
            print(f"    Method: {method}")
            if process_animal_pooled(mouse_id, method):
                total_animal_plots += 1
        
        # Process animal-concatenated data
        print(f"\n  ANIMAL-CONCATENATED DATA:")
        for method in METHODS:
            print(f"    Method: {method}")
            if process_animal_concatenated(mouse_id, method):
                total_animal_concat_plots += 1
    
    print(f"\n{'='*80}")
    print(f"COMPLETE!")
    print(f"  Session-pooled plot sets generated: {total_session_plots}")
    print(f"  Animal-pooled plot sets generated: {total_animal_plots}")
    print(f"  Animal-concatenated plot sets generated: {total_animal_concat_plots}")
    print(f"{'='*80}")


if __name__ == '__main__':
    main()
