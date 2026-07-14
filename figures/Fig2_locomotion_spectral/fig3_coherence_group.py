"""
================================================================================
FIGURE 3: Group-Level LFP-GEVI Coherence and PSD Analysis
================================================================================

FIGURE ORGANIZATION:
  - Coherence: 1 row x 2 columns (Session-Pooled | Animal-Pooled), REST vs RUN
  - PSD: 2 rows x 2 columns (LFP top, GEVI bottom; Session left, Animal right), REST vs RUN
  - Band-Averaged: Box plots for θ, α, β, γ bands

USAGE:
------
  Option 1: Run via master pipeline (recommended)
      python run_all_plots.py

  Option 2: Run standalone
      python fig3_coherence_group.py

================================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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
        get_group_level_input_dir, get_group_level_output_dir, get_animals_to_process,
        FONT_SIZE_SUPTITLE, FONT_SIZE_TITLE, FONT_SIZE_LABEL, FONT_SIZE_TICK,
        FONT_SIZE_LEGEND, FONT_SIZE_BAND, AXIS_LINEWIDTH, TICK_WIDTH, TICK_LENGTH,
        LINE_WIDTH_TRACE, COLOR_LFP_REST, COLOR_LFP_RUN, COLOR_GEVI_REST,
        COLOR_GEVI_RUN, COLOR_COH_REST, COLOR_COH_RUN, BAND_COLORS, SEM_ALPHA,
        GROUP_POOLING_LEVEL,
    )
    USING_CENTRAL_CONFIG = True
    INPUT_DIR = get_group_level_input_dir()
    OUTPUT_DIR = get_group_level_output_dir()
except ImportError:
    USING_CENTRAL_CONFIG = False
    # Fallback defaults (only used if plotting_config.py can't be found)
    BASE_OUTPUT_DIR = Path(__file__).resolve().parents[2] / "Figures" / "Spectral_data_outputs"
    BEHAVIOR_MODE = 'clear'
    # Group-level is NOT behavior-mode dependent
    INPUT_DIR = BASE_OUTPUT_DIR / "group_level_stats"
    OUTPUT_DIR = BASE_OUTPUT_DIR / "Spectral_python_figures" / "group_level"
    FIGURE_DPI = 300
    FIGURE_FORMATS = ['png', 'pdf', 'svg']
    FREQ_MIN = 2
    FREQ_MAX = 70
    GROUP_POOLING_LEVEL = 'animal_concatenated'  # Default fallback
    
    # Fallback for get_animals_to_process (empty list - FOOOF will be skipped)
    def get_animals_to_process():
        print("  WARNING: Central config not available, cannot get animal list")
        return []


# ==============================================================================
#  PUBLICATION FONT SIZES
# ==============================================================================
FONT_SIZE_SUPTITLE = 16
FONT_SIZE_TITLE = 14
FONT_SIZE_LABEL = 12
FONT_SIZE_TICK = 11
FONT_SIZE_LEGEND = 10
FONT_SIZE_BAND = 11

AXIS_LINEWIDTH = 2.0
TICK_WIDTH = 1.8
TICK_LENGTH = 7

LINE_WIDTH_TRACE = 2.5


# ==============================================================================
#  COLOR DEFINITIONS
# ==============================================================================

# LFP colors (PURPLE-GREY)
COLOR_LFP_REST = np.array([0.25, 0.18, 0.35])
COLOR_LFP_RUN = np.array([0.55, 0.45, 0.65])

# GEVI colors (TEAL)
COLOR_GEVI_REST = np.array([0.05, 0.35, 0.45])
COLOR_GEVI_RUN = np.array([0.25, 0.65, 0.65])

# Coherence colors - use from central config (Indigo/blue-violet blend)

SEM_ALPHA = 0.25
COLOR_SIG = 'red'
SIG_MARKER_SIZE = 10

# Band colors for box plots
BAND_COLORS = {
    'theta': np.array([0.4, 0.2, 0.6]),
    'alpha': np.array([0.2, 0.5, 0.5]),
    'beta': np.array([0.3, 0.6, 0.3]),
    'gamma': np.array([0.7, 0.5, 0.2]),
}


# ==============================================================================
#  WINDOWS LONG PATH SUPPORT
# ==============================================================================

from common import load_matlab_struct, to_long_path  # shared helpers (were local copies)


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================


def get_nested_field(obj, *keys, default=None):
    """Safely extract nested field from MATLAB data."""
    try:
        result = obj
        for key in keys:
            if result is None:
                return default
            if isinstance(result, dict) and key in result:
                result = result[key]
            elif hasattr(result, key):
                result = getattr(result, key)
            elif hasattr(result, '__dict__') and key in result.__dict__:
                result = result.__dict__[key]
            else:
                return default
        
        if result is None:
            return default
        if isinstance(result, np.ndarray):
            if result.ndim == 0:
                return result.item()
            return result.flatten() if result.ndim <= 2 and (result.ndim == 1 or 1 in result.shape) else result
        return result
    except:
        return default


def style_axis_publication(ax):
    """Apply publication-quality styling."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis='both', which='major', width=TICK_WIDTH, length=TICK_LENGTH, labelsize=FONT_SIZE_TICK)
    ax.grid(False)


def add_frequency_band_lines(ax):
    """Add colored dotted lines at frequency band boundaries."""
    bands = [(4, 'θ', 'theta'), (8, 'α', 'alpha'), (12, 'β', 'beta'), (30, 'γ', 'gamma')]
    centers = [6, 10, 21, 45]
    
    y_min, y_max = ax.get_ylim()
    
    for (f_start, label, key), center in zip(bands, centers):
        color = BAND_COLORS[key]
        if f_start >= FREQ_MIN and f_start <= FREQ_MAX:
            ax.axvline(x=f_start, color=color, alpha=0.5, linestyle=':', linewidth=1.5)
        if center <= FREQ_MAX:
            ax.text(center, y_max * 0.95, label, ha='center', va='top', fontsize=FONT_SIZE_BAND, fontweight='bold', color=color)


def plot_with_sem(ax, freq, mean, sem, color, linestyle='-', label=None):
    """Plot spectrum with SEM shading."""
    if mean is None or freq is None:
        return None
    if sem is not None:
        ax.fill_between(freq, mean - sem, mean + sem, color=color, alpha=SEM_ALPHA, linewidth=0)
    h, = ax.plot(freq, mean, linestyle, color=color, linewidth=LINE_WIDTH_TRACE, label=label)
    return h


def add_significance_markers(ax, freq, sig_mask):
    """Add red dots at significant frequencies."""
    if sig_mask is None:
        return
    sig_mask = np.asarray(sig_mask).flatten().astype(bool)
    if len(sig_mask) == len(freq) and np.any(sig_mask):
        y_min, y_max = ax.get_ylim()
        sig_y = y_min + (y_max - y_min) * 0.05
        sig_freqs = freq[sig_mask]
        ax.scatter(sig_freqs, np.full_like(sig_freqs, sig_y), c=COLOR_SIG, s=SIG_MARKER_SIZE, marker='o', zorder=5)


def get_sig_mask(data, measure_key, sig_source='standard'):
    """Extract significance mask from data."""
    if data is None or measure_key not in data:
        return None
    measure = data[measure_key]
    if sig_source == 'cluster':
        sig_mask = get_nested_field(measure, 'cluster_stats', 'sig_mask')
    else:
        sig_mask = get_nested_field(measure, 'stats', 'sig_mask')
    return sig_mask


def get_freq_axis(data, psd=False):
    """Get frequency axis from data."""
    if psd:
        freq = get_nested_field(data, 'freq_psd')
        if freq is None:
            freq = get_nested_field(data, 'freq')
    else:
        freq = get_nested_field(data, 'freq')
    return freq


# ==============================================================================
#  HELPER: Format pooling level name for display
# ==============================================================================

def get_pooling_level_label():
    """Get human-readable label for the current pooling level."""
    return GROUP_POOLING_LEVEL.replace('_', '-').title()


# ==============================================================================
#  COHERENCE FIGURE: Adapts layout based on available data
#  - If both session and animal data: 2 columns (Session | Animal-level)
#  - If only animal data: 1 column (Animal-level only)
# ==============================================================================

def create_coherence_figure(data_session, data_animal, output_path, method_title, sig_source='standard', show_significance=True):
    """
    Layout adapts based on data availability:
      - If data_session is provided: 1 row × 2 columns (Session-Pooled | Animal-level)
      - If only data_animal: 1 row × 1 column (Animal-level only)
    
    Parameters:
        show_significance: If True, display red significance markers. If False, omit them.
    """
    
    # Get the animal-level label based on config
    animal_level_label = get_pooling_level_label()
    
    # Determine layout: single panel (animal only) or dual panel (session + animal)
    single_panel = data_session is None
    
    if single_panel:
        fig, ax_single = plt.subplots(1, 1, figsize=(7, 5))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.85, bottom=0.12)
        data_list = [(data_animal, animal_level_label, ax_single)]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.subplots_adjust(left=0.10, right=0.95, top=0.85, bottom=0.12, wspace=0.25)
        data_list = [(data_session, 'Session-Pooled', axes[0]), (data_animal, animal_level_label, axes[1])]
    
    sig_label = 'Cluster' if sig_source == 'cluster' else 'FDR'
    
    for data, pool_label, ax in data_list:
        has_data = False
        if data is not None and 'coherence' in data:
            freq = get_freq_axis(data)
            rest_mean = get_nested_field(data, 'coherence', 'REST', 'mean')
            rest_sem = get_nested_field(data, 'coherence', 'REST', 'sem')
            run_mean = get_nested_field(data, 'coherence', 'RUN', 'mean')
            run_sem = get_nested_field(data, 'coherence', 'RUN', 'sem')
            sig_mask = get_sig_mask(data, 'coherence', sig_source) if show_significance else None
            
            if freq is not None and rest_mean is not None and run_mean is not None:
                has_data = True
                plot_with_sem(ax, freq, rest_mean, rest_sem, COLOR_COH_REST, '-', 'REST')
                plot_with_sem(ax, freq, run_mean, run_sem, COLOR_COH_RUN, '--', 'RUN')
                ax.set_xlim([FREQ_MIN, FREQ_MAX])
                ax.set_ylim([0, 1])
                add_frequency_band_lines(ax)
                if show_significance:
                    add_significance_markers(ax, freq, sig_mask)
                ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
                n = data.get('num_animals', get_nested_field(data, 'num_animals', default='?'))
                ax.set_title(f'{pool_label} (N={n})', fontsize=FONT_SIZE_TITLE, fontweight='bold')
        
        if not has_data:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
        
        ax.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
        ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        style_axis_publication(ax)
    
    # Title changes based on whether significance is shown
    if show_significance:
        fig.suptitle(f'LFP-GEVI Coherence: REST vs RUN\n{method_title} ({sig_label} sig = red)', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    else:
        fig.suptitle(f'LFP-GEVI Coherence: REST vs RUN\n{method_title}', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    
    # Ensure output directory exists (with long path support)
    os.makedirs(to_long_path(str(output_path.parent)), exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.name}")


def create_coherence_figure_logfreq(data_session, data_animal, output_path, method_title, sig_source='standard', show_significance=True):
    """
    Log-frequency variant of create_coherence_figure.
    Keeps the same REST/RUN means+SEM and frequency-band separators.
    """
    animal_level_label = get_pooling_level_label()
    single_panel = data_session is None

    if single_panel:
        fig, ax_single = plt.subplots(1, 1, figsize=(7, 5))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.85, bottom=0.12)
        data_list = [(data_animal, animal_level_label, ax_single)]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.subplots_adjust(left=0.10, right=0.95, top=0.85, bottom=0.12, wspace=0.25)
        data_list = [(data_session, 'Session-Pooled', axes[0]), (data_animal, animal_level_label, axes[1])]

    sig_label = 'Cluster' if sig_source == 'cluster' else 'FDR'

    for data, pool_label, ax in data_list:
        has_data = False
        if data is not None and 'coherence' in data:
            freq = get_freq_axis(data)
            rest_mean = get_nested_field(data, 'coherence', 'REST', 'mean')
            rest_sem = get_nested_field(data, 'coherence', 'REST', 'sem')
            run_mean = get_nested_field(data, 'coherence', 'RUN', 'mean')
            run_sem = get_nested_field(data, 'coherence', 'RUN', 'sem')
            sig_mask = get_sig_mask(data, 'coherence', sig_source) if show_significance else None

            if freq is not None and rest_mean is not None and run_mean is not None:
                freq = np.asarray(freq, dtype=float).flatten()
                pos = freq > 0
                if np.sum(pos) >= 2:
                    f = freq[pos]
                    rest_mean = np.asarray(rest_mean).flatten()[pos]
                    run_mean = np.asarray(run_mean).flatten()[pos]
                    rest_sem = np.asarray(rest_sem).flatten()[pos] if rest_sem is not None else None
                    run_sem = np.asarray(run_sem).flatten()[pos] if run_sem is not None else None
                    sig_mask = np.asarray(sig_mask).flatten()[pos] if sig_mask is not None else None

                    has_data = True
                    plot_with_sem(ax, f, rest_mean, rest_sem, COLOR_COH_REST, '-', 'REST')
                    plot_with_sem(ax, f, run_mean, run_sem, COLOR_COH_RUN, '--', 'RUN')
                    ax.set_xscale('log')
                    ax.set_xlim([max(FREQ_MIN, 1), FREQ_MAX])
                    ax.set_ylim([0, 1])
                    add_frequency_band_lines(ax)
                    if show_significance:
                        add_significance_markers(ax, f, sig_mask)
                    ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
                    n = data.get('num_animals', get_nested_field(data, 'num_animals', default='?'))
                    ax.set_title(f'{pool_label} (N={n})', fontsize=FONT_SIZE_TITLE, fontweight='bold')

        if not has_data:
            ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')

        ax.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
        ax.set_xlabel('Frequency (Hz, log scale)', fontsize=FONT_SIZE_LABEL)
        style_axis_publication(ax)

    if show_significance:
        fig.suptitle(f'LFP-GEVI Coherence: REST vs RUN (log frequency)\n{method_title} ({sig_label} sig = red)',
                     fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    else:
        fig.suptitle(f'LFP-GEVI Coherence: REST vs RUN (log frequency)\n{method_title}',
                     fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')

    os.makedirs(to_long_path(str(output_path.parent)), exist_ok=True)
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.name}")


# ==============================================================================
#  PSD FIGURE: Adapts layout based on available data
#  - If both session and animal data: 2 rows × 2 columns (LFP/GEVI × Session/Animal)
#  - If only animal data: 2 rows × 1 column (LFP/GEVI × Animal-level only)
# ==============================================================================

def create_psd_figure(data_session, data_animal, output_path, method_title, sig_source='standard', show_significance=True):
    """
    Layout adapts based on data availability:
      - If data_session is provided: 2 rows × 2 columns
      - If only data_animal: 2 rows × 1 column (Animal-level only)
    
    Parameters:
        show_significance: If True, display red significance markers. If False, omit them.
    """
    
    # Get the animal-level label based on config
    animal_level_label = get_pooling_level_label()
    
    # Determine layout: single column (animal only) or dual column (session + animal)
    single_column = data_session is None
    
    if single_column:
        fig, axes = plt.subplots(2, 1, figsize=(7, 10))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.08, hspace=0.30)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.subplots_adjust(left=0.10, right=0.95, top=0.90, bottom=0.08, wspace=0.25, hspace=0.30)
    
    sig_label = 'Cluster' if sig_source == 'cluster' else 'FDR'
    
    # Row 0: LFP PSD, Row 1: GEVI PSD
    for row_idx, (psd_key, signal_label, rest_color, run_color) in enumerate([
        ('psd_lfp', 'LFP', COLOR_LFP_REST, COLOR_LFP_RUN),
        ('psd_gevi', 'GEVI', COLOR_GEVI_REST, COLOR_GEVI_RUN)
    ]):
        if single_column:
            # Single column: only animal-level
            data_list = [(data_animal, animal_level_label, axes[row_idx])]
        else:
            # Dual column: session-pooled (col 0) and animal-level (col 1)
            data_list = [
                (data_session, 'Session-Pooled', axes[row_idx, 0]),
                (data_animal, animal_level_label, axes[row_idx, 1])
            ]
        
        for data, pool_label, ax in data_list:
            has_data = False
            if data is not None and psd_key in data:
                freq = get_freq_axis(data, psd=True)
                rest_mean = get_nested_field(data, psd_key, 'REST', 'mean')
                rest_sem = get_nested_field(data, psd_key, 'REST', 'sem')
                run_mean = get_nested_field(data, psd_key, 'RUN', 'mean')
                run_sem = get_nested_field(data, psd_key, 'RUN', 'sem')
                sig_mask = get_sig_mask(data, psd_key, sig_source) if show_significance else None
                
                if freq is not None and rest_mean is not None and run_mean is not None:
                    has_data = True
                    plot_with_sem(ax, freq, rest_mean, rest_sem, rest_color, '-', 'REST')
                    plot_with_sem(ax, freq, run_mean, run_sem, run_color, '--', 'RUN')
                    ax.set_xlim([FREQ_MIN, FREQ_MAX])
                    add_frequency_band_lines(ax)
                    if show_significance:
                        add_significance_markers(ax, freq, sig_mask)
                    ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
                    n = data.get('num_animals', get_nested_field(data, 'num_animals', default='?'))
                    ax.set_title(f'{signal_label} PSD – {pool_label} (N={n})', fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            if not has_data:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
                ax.set_title(f'{signal_label} PSD – {pool_label}', fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            ax.set_ylabel(f'{signal_label} Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
            ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
            style_axis_publication(ax)
    
    # Title changes based on whether significance is shown
    if show_significance:
        fig.suptitle(f'Power Spectral Density: REST vs RUN\n{method_title} ({sig_label} sig = red)', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    else:
        fig.suptitle(f'Power Spectral Density: REST vs RUN\n{method_title}', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    
    # Ensure output directory exists (with long path support)
    os.makedirs(to_long_path(str(output_path.parent)), exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.name}")


def create_psd_figure_loglog(data_session, data_animal, output_path, method_title, sig_source='standard', show_significance=True):
    """
    Log-log scale PSD figure for 1/f analysis.
    Same layout as create_psd_figure but with log frequency axis.
    """
    
    animal_level_label = get_pooling_level_label()
    single_column = data_session is None
    
    if single_column:
        fig, axes = plt.subplots(2, 1, figsize=(7, 10))
        fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.08, hspace=0.30)
    else:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.subplots_adjust(left=0.10, right=0.95, top=0.90, bottom=0.08, wspace=0.25, hspace=0.30)
    
    sig_label = 'Cluster' if sig_source == 'cluster' else 'FDR'
    
    for row_idx, (psd_key, signal_label, rest_color, run_color) in enumerate([
        ('psd_lfp', 'LFP', COLOR_LFP_REST, COLOR_LFP_RUN),
        ('psd_gevi', 'GEVI', COLOR_GEVI_REST, COLOR_GEVI_RUN)
    ]):
        if single_column:
            data_list = [(data_animal, animal_level_label, axes[row_idx])]
        else:
            data_list = [
                (data_session, 'Session-Pooled', axes[row_idx, 0]),
                (data_animal, animal_level_label, axes[row_idx, 1])
            ]
        
        for data, pool_label, ax in data_list:
            has_data = False
            if data is not None and psd_key in data:
                freq = get_freq_axis(data, psd=True)
                rest_mean = get_nested_field(data, psd_key, 'REST', 'mean')
                rest_sem = get_nested_field(data, psd_key, 'REST', 'sem')
                run_mean = get_nested_field(data, psd_key, 'RUN', 'mean')
                run_sem = get_nested_field(data, psd_key, 'RUN', 'sem')
                sig_mask = get_sig_mask(data, psd_key, sig_source) if show_significance else None
                
                if freq is not None and rest_mean is not None and run_mean is not None:
                    # Filter to positive frequencies for log scale
                    freq_mask = freq > 0
                    freq_plot = freq[freq_mask]
                    rest_mean_plot = rest_mean[freq_mask]
                    run_mean_plot = run_mean[freq_mask]
                    rest_sem_plot = rest_sem[freq_mask] if rest_sem is not None else None
                    run_sem_plot = run_sem[freq_mask] if run_sem is not None else None
                    sig_mask_plot = sig_mask[freq_mask] if sig_mask is not None else None
                    
                    has_data = True
                    plot_with_sem(ax, freq_plot, rest_mean_plot, rest_sem_plot, rest_color, '-', 'REST')
                    plot_with_sem(ax, freq_plot, run_mean_plot, run_sem_plot, run_color, '--', 'RUN')
                    ax.set_xscale('log')
                    ax.set_xlim([FREQ_MIN, FREQ_MAX])
                    if show_significance and sig_mask_plot is not None:
                        add_significance_markers(ax, freq_plot, sig_mask_plot)
                    ax.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
                    n = data.get('num_animals', get_nested_field(data, 'num_animals', default='?'))
                    ax.set_title(f'{signal_label} PSD – {pool_label} (N={n})', fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            if not has_data:
                ax.text(0.5, 0.5, 'No data', transform=ax.transAxes, ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
                ax.set_title(f'{signal_label} PSD – {pool_label}', fontsize=FONT_SIZE_TITLE, fontweight='bold')
            
            ax.set_ylabel(f'{signal_label} Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
            ax.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
            ax.yaxis.set_major_locator(MultipleLocator(5))
            ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
            style_axis_publication(ax)
    
    if show_significance:
        fig.suptitle(f'Power Spectral Density (Log-Log): REST vs RUN\n{method_title} ({sig_label} sig = red)', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    else:
        fig.suptitle(f'Power Spectral Density (Log-Log): REST vs RUN\n{method_title}', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    
    os.makedirs(to_long_path(str(output_path.parent)), exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.name}")


# ==============================================================================
#  BAND-AVERAGED BOX PLOT FIGURE
# ==============================================================================

def create_band_boxplot_figure(data_session, data_animal, output_path, method_title):
    """
    Create box plots comparing REST vs RUN for each frequency band.
    Uses animal-level data (animal_pooled or animal_concatenated) if available, else session-pooled.
    """
    
    data = data_animal if data_animal is not None else data_session
    if data is None:
        print("    WARNING: No data for band boxplot")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.subplots_adjust(left=0.08, right=0.95, top=0.85, bottom=0.12, wspace=0.30)
    
    bands = ['theta', 'alpha', 'beta', 'gamma']
    band_labels = ['θ\n4-8', 'α\n8-12', 'β\n12-30', 'γ\n30-70']
    
    n_animals = data.get('num_animals', get_nested_field(data, 'num_animals', default='?'))
    
    # Panel configs: (measure_key, title, y_label)
    panels = [
        ('coherence', 'LFP-GEVI Coherence', 'Coherence'),
        ('psd_lfp', 'LFP Power', 'Power (dB/Hz)'),
        ('psd_gevi', 'GEVI Power', 'Power (dB/Hz)'),
    ]
    
    for ax_idx, (measure_key, title, ylabel) in enumerate(panels):
        ax = axes[ax_idx]
        
        band_stats = get_nested_field(data, measure_key, 'band_stats')
        
        if band_stats is None:
            ax.text(0.5, 0.5, 'No band stats', transform=ax.transAxes, ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
            ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
            style_axis_publication(ax)
            continue
        
        # Collect data for box plots
        rest_data = []
        run_data = []
        pvals = []
        cohens_ds = []
        
        for band in bands:
            bs = get_nested_field(band_stats, band)
            if bs is not None:
                # Try to get individual subject values
                rest_vals = get_nested_field(bs, 'rest_values')
                run_vals = get_nested_field(bs, 'run_values')
                
                if rest_vals is not None and run_vals is not None:
                    rest_data.append(np.asarray(rest_vals).flatten())
                    run_data.append(np.asarray(run_vals).flatten())
                else:
                    # Fallback to means if individual values not available
                    rest_m = get_nested_field(bs, 'rest_mean', default=0)
                    run_m = get_nested_field(bs, 'run_mean', default=0)
                    rest_data.append([rest_m])
                    run_data.append([run_m])
                
                pvals.append(get_nested_field(bs, 'pval_twotailed', default=1.0))
                cohens_ds.append(get_nested_field(bs, 'cohens_d', default=0.0))
            else:
                rest_data.append([0])
                run_data.append([0])
                pvals.append(1.0)
                cohens_ds.append(0.0)
        
        # Create box plots with explicit z-order and thick black edges
        positions_rest = np.arange(len(bands)) * 2
        positions_run = positions_rest + 0.7
        
        bp_rest = ax.boxplot(rest_data, positions=positions_rest, widths=0.5, patch_artist=True,
                             boxprops=dict(facecolor=COLOR_COH_REST, edgecolor='black', linewidth=2),
                             medianprops=dict(color='white', linewidth=2.5),
                             whiskerprops=dict(color='black', linewidth=1.5),
                             capprops=dict(color='black', linewidth=2),
                             flierprops=dict(marker='o', markerfacecolor=COLOR_COH_REST, markeredgecolor='black', markersize=5))
        
        bp_run = ax.boxplot(run_data, positions=positions_run, widths=0.5, patch_artist=True,
                            boxprops=dict(facecolor=COLOR_COH_RUN, edgecolor='black', linewidth=2),
                            medianprops=dict(color='white', linewidth=2.5),
                            whiskerprops=dict(color='black', linewidth=1.5),
                            capprops=dict(color='black', linewidth=2),
                            flierprops=dict(marker='o', markerfacecolor=COLOR_COH_RUN, markeredgecolor='black', markersize=5))
        
        # CRITICAL: Set explicit x-axis limits FIRST to prevent edge clipping
        x_min = positions_rest[0] - 0.6  # Extra padding on left
        x_max = positions_run[-1] + 0.6   # Extra padding on right
        ax.set_xlim(x_min, x_max)
        
        # Disable clipping for ALL boxplot elements
        for bp in [bp_rest, bp_run]:
            for element_list in [bp['whiskers'], bp['caps'], bp['boxes'], bp['medians'], bp['fliers']]:
                for elem in element_list:
                    elem.set_clip_on(False)
            # Set z-order and linewidths
            for whisker in bp['whiskers']:
                whisker.set_zorder(3)
                whisker.set_linewidth(2)
                whisker.set_color('black')
            for cap in bp['caps']:
                cap.set_zorder(3)
                cap.set_linewidth(2.5)
                cap.set_color('black')
            for box in bp['boxes']:
                box.set_zorder(2)
                box.set_linewidth(2.5)
                box.set_edgecolor('black')
            for median in bp['medians']:
                median.set_zorder(4)
                median.set_linewidth(2.5)
            for flier in bp['fliers']:
                flier.set_zorder(5)
        
        # Draw explicit whisker lines on top of everything to guarantee visibility
        for i, (rest_d, run_d) in enumerate(zip(rest_data, run_data)):
            for data_arr, pos in [(rest_d, positions_rest[i]), (run_d, positions_run[i])]:
                if len(data_arr) > 1:
                    q1, q3 = np.percentile(data_arr, [25, 75])
                    iqr = q3 - q1
                    # Calculate whisker bounds (1.5 * IQR rule)
                    whisker_low = max(np.min(data_arr), q1 - 1.5 * iqr)
                    whisker_high = min(np.max(data_arr), q3 + 1.5 * iqr)
                    # Draw explicit whisker lines at highest z-order
                    ax.plot([pos, pos], [q1, whisker_low], color='black', linewidth=2, zorder=10, clip_on=False)
                    ax.plot([pos, pos], [q3, whisker_high], color='black', linewidth=2, zorder=10, clip_on=False)
                    # Draw caps
                    cap_width = 0.15
                    ax.plot([pos - cap_width, pos + cap_width], [whisker_low, whisker_low], color='black', linewidth=2.5, zorder=10, clip_on=False)
                    ax.plot([pos - cap_width, pos + cap_width], [whisker_high, whisker_high], color='black', linewidth=2.5, zorder=10, clip_on=False)
                    # Draw box outline
                    rect = Rectangle((pos - 0.25, q1), 0.5, q3 - q1, 
                                     fill=False, edgecolor='black', linewidth=2.5, zorder=10, clip_on=False)
                    ax.add_patch(rect)
        
        # Expand y-axis to include all data with padding
        all_data = [val for sublist in rest_data + run_data for val in sublist]
        data_min, data_max = min(all_data), max(all_data)
        data_range = data_max - data_min
        ax.set_ylim([data_min - 0.1 * data_range, data_max + 0.15 * data_range])
        
        # Add significance stars above the boxes
        y_min, y_max = ax.get_ylim()
        for i, (p, d) in enumerate(zip(pvals, cohens_ds)):
            x_pos = positions_rest[i] + 0.35
            
            # Significance stars only
            star = ''
            if p < 0.001:
                star = '***'
            elif p < 0.01:
                star = '**'
            elif p < 0.05:
                star = '*'
            
            if star:
                ax.text(x_pos, y_max * 0.95, star, ha='center', va='top', fontsize=14, fontweight='bold', color='black')
        
        # Formatting
        ax.set_xticks(positions_rest + 0.35)
        ax.set_xticklabels(band_labels, fontsize=FONT_SIZE_TICK)
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
        
        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=COLOR_COH_REST, edgecolor='black', linewidth=1.5, label='REST'),
                          Patch(facecolor=COLOR_COH_RUN, edgecolor='black', linewidth=1.5, label='RUN')]
        ax.legend(handles=legend_elements, loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=False)
        
        style_axis_publication(ax)
    
    fig.suptitle(f'Band-Averaged REST vs RUN (N={n_animals}) – {method_title}', 
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold')
    
    # Ensure output directory exists (with long path support)
    os.makedirs(to_long_path(str(output_path.parent)), exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {output_path.name}")


# ==============================================================================
#  MAIN EXECUTION
# ==============================================================================

def main():
    """Main execution function."""
    
    print("=" * 80)
    print("FIGURE 3: Group-Level LFP-GEVI Analysis")
    print("=" * 80)
    print(f"  Behavior Mode: {BEHAVIOR_MODE}")
    print(f"  Group Pooling Level: {GROUP_POOLING_LEVEL}")
    print(f"  Input Directory: {INPUT_DIR}")
    print(f"  Output Directory: {OUTPUT_DIR}")
    print("=" * 80)
    
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    
    # Define methods to process
    methods = [
        ('mscohere', 'mscohere', 'standard'),
        ('fieldtrip', 'FieldTrip', 'standard'),
        ('fieldtrip_cluster', 'FieldTrip Cluster', 'cluster'),
    ]
    
    for file_prefix, method_title, sig_source in methods:
        print(f"\n{'='*60}")
        print(f"  {method_title}")
        print(f"{'='*60}")
        
        # Load data - MATLAB naming convention: {method}_{level}.mat
        # Session-pooled: {method}_session_pooled.mat
        # Animal-level: {method}_{GROUP_POOLING_LEVEL}.mat (configurable)
        session_path = INPUT_DIR / f"{file_prefix}_session_pooled.mat"
        animal_path = INPUT_DIR / f"{file_prefix}_{GROUP_POOLING_LEVEL}.mat"
        
        data_session = load_matlab_struct(session_path)
        data_animal = load_matlab_struct(animal_path)
        
        if data_session:
            n = data_session.get('num_animals', get_nested_field(data_session, 'num_animals', default='?'))
            print(f"  Loaded session-pooled: N={n}")
        if data_animal:
            n = data_animal.get('num_animals', get_nested_field(data_animal, 'num_animals', default='?'))
            print(f"  Loaded {GROUP_POOLING_LEVEL.replace('_', '-')}: N={n}")
        
        if data_session is None and data_animal is None:
            print(f"  No data found, skipping...")
            continue
        
        # Create Coherence figure (1 row, REST vs RUN)
        print("\n  Creating Coherence Figure...")
        coh_path = OUTPUT_DIR / f"figure3_{file_prefix}_coherence"
        create_coherence_figure(data_session, data_animal, coh_path, method_title, sig_source)
        
        # Create PSD figure (2 rows: LFP/GEVI, 2 cols: Session/Animal, all REST vs RUN)
        print("  Creating PSD Figure...")
        psd_path = OUTPUT_DIR / f"figure3_{file_prefix}_psd"
        create_psd_figure(data_session, data_animal, psd_path, method_title, sig_source)
        
        # Create Band-Averaged Box Plot
        print("  Creating Band-Averaged Box Plot...")
        band_path = OUTPUT_DIR / f"figure3_{file_prefix}_band_boxplot"
        create_band_boxplot_figure(data_session, data_animal, band_path, method_title)
    
    # Run FOOOF analysis (Figure 4)
    # NOTE: FOOOF analyzes PSD (Power Spectral Density), NOT coherence.
    # Since PSD is computed identically for both mscohere and FieldTrip methods
    # (both use pwelch), we only need to run FOOOF once using mscohere data.
    print(f"\n{'='*60}")
    print("  FOOOF Analysis (Figure 4)")
    print(f"{'='*60}")
    try:
        from fig4_fooof import run_fooof_pipeline, FOOOF_AVAILABLE, get_pooling_level_label
        
        if FOOOF_AVAILABLE:
            pooling_label = get_pooling_level_label(GROUP_POOLING_LEVEL)
            print(f"  Running FOOOF Analysis ({pooling_label})...")
            print("  (FOOOF analyzes PSD which is identical for both methods)")
            animals = get_animals_to_process()
            
            # Run FOOOF on mscohere-derived PSD (same as fieldtrip PSD since both use pwelch)
            fooof_plots = run_fooof_pipeline(animals, OUTPUT_DIR, method='mscohere', 
                                             pooling_level=GROUP_POOLING_LEVEL)
            
            if fooof_plots > 0:
                print(f"  Generated {fooof_plots} FOOOF figure(s)")
        else:
            print("  Skipping FOOOF (not installed - pip install fooof)")
    except ImportError as e:
        print(f"  Skipping FOOOF (import error: {e})")
    except Exception as e:
        import traceback
        print(f"  FOOOF Error: {e}")
        traceback.print_exc()
    
    # Run Theta Band Analysis (Figure 5)
    print(f"\n{'='*60}")
    print("  Theta Band Analysis (Figure 5)")
    print(f"{'='*60}")
    try:
        from fig5_theta import run_theta_band_pipeline, get_pooling_level_label
        
        pooling_label = get_pooling_level_label(GROUP_POOLING_LEVEL)
        print(f"  Running Theta Band Analysis ({pooling_label})...")
        print(f"  Hypothesis: RUN > REST (one-tailed test)")
        animals = get_animals_to_process()
        
        theta_plots = run_theta_band_pipeline(animals, OUTPUT_DIR, method='mscohere',
                                              pooling_level=GROUP_POOLING_LEVEL)
        
        if theta_plots > 0:
            print(f"  Generated {theta_plots} Theta Band figure(s)")
    except ImportError as e:
        print(f"  Skipping Theta Band (import error: {e})")
    except Exception as e:
        import traceback
        print(f"  Theta Band Error: {e}")
        traceback.print_exc()
    
    print("\n" + "=" * 80)
    print("COMPLETE!")
    print("=" * 80)
    print(f"\nFigures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
