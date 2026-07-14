"""
================================================================================
FIGURE 2: LFP-GEVI Coherence Analysis - Python Plotting Script (Single-Trial)
================================================================================

This script loads the MATLAB-generated coherence analysis results and creates
publication-quality figures showing:
  - Running speed (as heatmap strip)
  - LFP spectrogram
  - LFP-GEVI coherence heatmap (time-frequency)
  - Coherence spectrum (coherence vs frequency) with rest/run breakdown
  - Power Spectral Density (PSD) for LFP and GEVI

USAGE:
------
  Option 1: Run via master pipeline (recommended)
      python run_all_plots.py

  Option 2: Run standalone
      1. Edit MOUSE_ID and SESSION_ID below
      2. python fig2_coherence.py

================================================================================
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, NullFormatter
from scipy.interpolate import RegularGridInterpolator
from pathlib import Path
import warnings

# Suppress scipy loadmat warnings about MATLAB structs
warnings.filterwarnings('ignore', category=UserWarning)

# common.py / plotting_config.py live in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# ==============================================================================
#  CONFIGURATION - Import from central config or use defaults
# ==============================================================================

try:
    from plotting_config import (
        BEHAVIOR_MODE, FIGURE_DPI, FIGURE_FORMATS, FREQ_MIN, FREQ_MAX,
        COH_VMIN, COH_VMAX, get_single_trial_input_dir, get_single_trial_output_dir,
        FONT_SIZE_TITLE, FONT_SIZE_SUPTITLE, FONT_SIZE_LABEL, FONT_SIZE_TICK,
        FONT_SIZE_LEGEND, FONT_SIZE_BAND, AXIS_LINEWIDTH, TICK_WIDTH, TICK_LENGTH,
        LINE_WIDTH_TRACE, COLOR_REST, COLOR_RUN, COLOR_OVERALL, COLOR_LFP, COLOR_GEVI,
        BAND_LINE_COLORS, SEM_ALPHA, BASE_OUTPUT_DIR,
    )
    USING_CENTRAL_CONFIG = True
    FIGURE_FORMAT = FIGURE_FORMATS
except ImportError:
    USING_CENTRAL_CONFIG = False
    # Fallback defaults (only used if plotting_config.py can't be found)
    BASE_OUTPUT_DIR = Path(__file__).resolve().parents[2] / "Figures" / "Spectral_data_outputs"
    BEHAVIOR_MODE = 'clear'
    FIGURE_DPI = 300
    FIGURE_FORMAT = ['png', 'pdf', 'svg']
    FREQ_MIN = 2
    FREQ_MAX = 70
    COH_VMIN = 0.0
    COH_VMAX = 1.0

# ------------------------------------------------------------------------------
#  STANDALONE MODE: Set these when running this script directly
# ------------------------------------------------------------------------------
MOUSE_ID = "Animal01"
SESSION_ID = "01_09_25-R1"

# Compute paths (will be overridden when called from master pipeline)
if USING_CENTRAL_CONFIG:
    INPUT_DIR = get_single_trial_input_dir(MOUSE_ID, SESSION_ID)
    OUTPUT_DIR = get_single_trial_output_dir(MOUSE_ID, SESSION_ID)
else:
    INPUT_DIR = BASE_OUTPUT_DIR / BEHAVIOR_MODE / "single_trial" / MOUSE_ID / SESSION_ID / "data"
    OUTPUT_DIR = BASE_OUTPUT_DIR / BEHAVIOR_MODE / "single_trial" / MOUSE_ID / SESSION_ID / "figures"

# Trial labels (auto-generated if None)
TRIAL_LABELS = None

# Colormaps
CMAP_SPECTROGRAM = 'viridis'

# ==============================================================================
#  PUBLICATION FONT SIZES - LARGER FOR CLARITY
# ==============================================================================
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 20
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_STATS = 13
FONT_SIZE_COLORBAR = 14
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
#  COLOR DEFINITIONS
# ==============================================================================
# Teal shades for rest vs run (different from Figure 1's teal which is [0.127568, 0.566949, 0.550556])
COLOR_REST = np.array([0.05, 0.35, 0.45])      # Darker teal for rest
COLOR_RUN = np.array([0.25, 0.65, 0.65])       # Lighter teal for run
COLOR_OVERALL = np.array([0.08, 0.45, 0.52])   # Mid teal for overall

# Colors for LFP and GEVI traces (consistent with Figure 1)
# Matched to fig1_gevi_lfp.py color scheme
COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])  # Teal (from viridis)
COLOR_LFP = np.array([0.35, 0.25, 0.45])               # Purple-grey

# Frequency band line colors (specific colors for each band)
BAND_LINE_COLORS = {
    'theta': (0.4, 0.2, 0.6, 0.6),    # Purple
    'alpha': (0.2, 0.5, 0.5, 0.6),    # Teal
    'beta': (0.3, 0.6, 0.3, 0.6),     # Green
    'gamma': (0.7, 0.5, 0.2, 0.6),    # Orange/brown
}

# Frequency band boundaries (Hz) - note: alpha ends at 12, beta starts at 12 (no gap)
FREQ_BANDS = {
    'theta': (4, 8, 'θ'),
    'alpha': (8, 12, 'α'),
    'beta': (12, 30, 'β'),   # Changed from 13 to 12 to remove gap
    'gamma': (30, 70, 'γ'),
}


# ==============================================================================
#  WINDOWS LONG PATH SUPPORT
# ==============================================================================

from common import load_matlab_struct as load_matlab_data, to_long_path  # shared helpers (were local copies)


# ==============================================================================
#  CUSTOM COLORMAPS
# ==============================================================================

from common import create_monochromatic_orange_cmap  # shared helpers (were local copies)

CMAP_SPEED = create_monochromatic_orange_cmap()


def create_dark_viridis_cmap():
    """
    Create a modified viridis colormap with darker low values.
    This enhances contrast for low coherence regions, making the dark hues darker.
    """
    # Get the standard viridis colormap
    viridis = plt.cm.viridis(np.linspace(0, 1, 256))
    
    # Darken the lower portion of the colormap (first ~30% of values)
    # Apply a power function to shift towards darker colors for low values
    n_dark = 80  # Number of color stops to modify (out of 256)
    darkening_factor = 0.5  # How much to darken (0 = black, 1 = original)
    
    for i in range(n_dark):
        # Gradual transition: more darkening at very low values
        blend = (i / n_dark) ** 0.7  # Power < 1 means more darkening retained longer
        viridis[i, :3] = viridis[i, :3] * (darkening_factor + (1 - darkening_factor) * blend)
    
    return LinearSegmentedColormap.from_list('dark_viridis', viridis, N=256)

# Use the dark viridis for coherence to make low values appear darker
CMAP_COHERENCE_DARK = create_dark_viridis_cmap()


# ==============================================================================
#  HELPER FUNCTIONS
# ==============================================================================

def ensure_2d(arr):
    """Ensure array is 2D (Nf × Nt)."""
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr


def style_axis_publication(ax, remove_top_right=True):
    """
    Apply publication-quality styling to an axis.
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
                # Convert center frequency to axes coordinates
                x_data_to_axes = (center - FREQ_MIN) / (FREQ_MAX - FREQ_MIN)
                ax.text(x_data_to_axes, 1.02, label, ha='center', va='bottom',
                       fontsize=FONT_SIZE_BAND, fontweight='bold', 
                       color=color[:3], transform=ax.transAxes)


# ==============================================================================
#  MAIN FIGURE: Speed + Spectrogram + Coherence Heatmap
# ==============================================================================

def _interpolate_coherence_to_fine_freq(coh_matrix, freq_orig, time_orig, target_nfreq=None):
    """
    Interpolate coherence matrix to a finer frequency grid for sharper display.
    
    The coherence is typically computed with fewer frequency bins than the LFP 
    spectrogram (especially for fiber data at ~500 Hz vs LFP at 1 kHz effective).
    This upsamples in the frequency dimension using linear interpolation.
    
    Parameters
    ----------
    coh_matrix : ndarray, shape (n_freq, n_time)
        Original coherence matrix
    freq_orig : ndarray
        Original frequency vector
    time_orig : ndarray
        Original time vector
    target_nfreq : int or None
        Target number of frequency bins. If None, uses 4x original resolution.
    
    Returns
    -------
    coh_interp : ndarray
        Interpolated coherence matrix
    freq_fine : ndarray
        New fine-grained frequency vector
    """
    n_freq_orig = len(freq_orig)
    if target_nfreq is None:
        target_nfreq = n_freq_orig * 4
    
    freq_fine = np.linspace(freq_orig[0], freq_orig[-1], target_nfreq)
    
    interp = RegularGridInterpolator(
        (freq_orig, time_orig), coh_matrix,
        method='linear', bounds_error=False, fill_value=None
    )
    
    freq_grid, time_grid = np.meshgrid(freq_fine, time_orig, indexing='ij')
    points = np.column_stack([freq_grid.ravel(), time_grid.ravel()])
    coh_interp = interp(points).reshape(len(freq_fine), len(time_orig))
    
    return coh_interp, freq_fine


def create_main_figure(data_list, trial_labels, method_name, output_path):
    """
    Create main figure with speed heatmap, spectrogram, and coherence heatmaps.
    All panels have consistent time axis and uniform colorbar limits across trials.
    
    Features:
    - Coherence heatmap is interpolated to finer frequency grid for sharper display
    - Logarithmic frequency axis emphasizes low-frequency activity
    - Colorbar limits have headroom beyond data extremes to avoid saturation
    
    Parameters
    ----------
    data_list : list of dict
        List of trial data dictionaries loaded from MATLAB .mat files
    trial_labels : list of str
        Labels for each trial (e.g., ['Trial 1', 'Trial 2', 'Trial 3'])
    method_name : str
        Analysis method name ('mscohere' or 'fieldtrip')
    output_path : Path
        Base path for saving output files
    """
    num_trials = len(data_list)
    
    # Check if motion data is available in any trial
    has_motion = any('motion' in data for data in data_list)
    
    # Get common time range across all trials for consistency
    time_ranges = []
    for data in data_list:
        if 'spec_time' in data:
            plot_time = np.asarray(data['spec_time']).flatten()
        else:
            plot_time = np.asarray(data['time']).flatten()
        time_ranges.append((plot_time[0], plot_time[-1]))
    t_min = max(tr[0] for tr in time_ranges)
    t_max = min(tr[1] for tr in time_ranges)
    
    # Get uniform spectrogram colorbar limits across all trials (with headroom)
    all_spec_p5 = []
    all_spec_p95 = []
    for data in data_list:
        spec = ensure_2d(np.asarray(data['spec_power']))
        all_spec_p5.append(np.nanpercentile(spec, 2))
        all_spec_p95.append(np.nanpercentile(spec, 98))
    spec_data_min = min(all_spec_p5)
    spec_data_max = max(all_spec_p95)
    spec_range = spec_data_max - spec_data_min
    spec_vmin = spec_data_min - 0.1 * spec_range
    spec_vmax = spec_data_max + 0.15 * spec_range
    
    # Coherence colorbar limits with headroom to avoid saturation
    # Compute data-driven limits from all trials
    all_coh_p2 = []
    all_coh_p98 = []
    for data in data_list:
        if method_name == 'mscohere' and 'coh_mscohere' in data:
            coh_tmp = ensure_2d(np.asarray(data['coh_mscohere']))
        elif method_name == 'fieldtrip' and 'coh_fieldtrip' in data:
            coh_tmp = ensure_2d(np.asarray(data['coh_fieldtrip']))
        else:
            for key in ['coh_mscohere', 'coh_fieldtrip']:
                if key in data:
                    coh_tmp = ensure_2d(np.asarray(data[key]))
                    break
        all_coh_p2.append(np.nanpercentile(coh_tmp, 1))
        all_coh_p98.append(np.nanpercentile(coh_tmp, 99))
    
    coh_data_min = min(all_coh_p2)
    coh_data_max = max(all_coh_p98)
    coh_range = coh_data_max - coh_data_min
    # Clamp to valid coherence range [0, 1] but add headroom
    coh_vmin = max(0.0, coh_data_min - 0.15 * coh_range)
    coh_vmax = min(1.0, coh_data_max + 0.20 * coh_range)
    
    # Get uniform speed colorbar limits across all trials
    speed_vmin = 0
    speed_vmax = 0
    if has_motion:
        for data in data_list:
            if 'motion' in data:
                motion = np.asarray(data['motion']).flatten()
                speed_vmax = max(speed_vmax, np.nanmax(motion))
    
    # Dynamic figure sizing based on number of trials
    fig_width = max(8, 5 + 5 * num_trials)
    fig_width = min(fig_width, 30)
    
    # Figure setup with dynamic layout
    if has_motion:
        fig = plt.figure(figsize=(fig_width, 12))
        gs = fig.add_gridspec(3, num_trials, height_ratios=[0.12, 1, 1], 
                              hspace=0.40, wspace=0.30,
                              left=0.06, right=0.92, top=0.92, bottom=0.08)
    else:
        fig = plt.figure(figsize=(fig_width, 10))
        gs = fig.add_gridspec(2, num_trials, height_ratios=[1, 1], 
                              hspace=0.40, wspace=0.30,
                              left=0.06, right=0.92, top=0.92, bottom=0.08)
    
    # Store axes for colorbar positioning
    speed_axes = []
    spec_axes = []
    coh_axes = []
    
    for col_idx, (data, trial_label) in enumerate(zip(data_list, trial_labels)):
        # Coherence time/freq axes
        time = np.asarray(data['time']).flatten()
        freq = np.asarray(data['freq']).flatten()
        
        # Spectrogram time/freq axes (may be different resolution)
        if 'spec_time' in data:
            spec_time = np.asarray(data['spec_time']).flatten()
        else:
            spec_time = time
        if 'spec_freq' in data:
            spec_freq = np.asarray(data['spec_freq']).flatten()
        else:
            spec_freq = freq
        
        row_offset = 0
        
        # ========== ROW 1: Speed Heatmap (thin strip) ==========
        if has_motion and 'motion' in data:
            ax_speed = fig.add_subplot(gs[0, col_idx])
            speed_axes.append(ax_speed)
            motion = np.asarray(data['motion']).flatten()
            
            motion_2d = motion.reshape(1, -1)
            
            if 'motion_time' in data:
                motion_time = np.asarray(data['motion_time']).flatten()
            else:
                motion_time = spec_time
            
            extent = [motion_time[0], motion_time[-1], 0, 1]
            im_speed = ax_speed.imshow(motion_2d, aspect='auto', cmap=CMAP_SPEED,
                                        extent=extent, origin='lower',
                                        interpolation='bilinear',
                                        vmin=speed_vmin, vmax=speed_vmax)
            
            ax_speed.set_xlim([t_min, t_max])
            ax_speed.set_ylim([0, 1])
            ax_speed.set_yticks([])
            ax_speed.set_xticklabels([])
            ax_speed.set_title(trial_label, fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=35)
            
            for spine in ax_speed.spines.values():
                spine.set_visible(False)
            
            cax = ax_speed.inset_axes([0.15, 1.20, 0.45, 0.30])
            cbar = fig.colorbar(im_speed, cax=cax, orientation='horizontal')
            cbar.ax.tick_params(labelsize=FONT_SIZE_COLORBAR)
            cbar.ax.xaxis.set_label_position('top')
            cbar.ax.xaxis.tick_top()
            
            ax_speed.text(0.72, 1.35, 'Speed (cm/s)', transform=ax_speed.transAxes,
                         fontsize=FONT_SIZE_COLORBAR + 2, fontweight='bold',
                         ha='left', va='center')
            
            row_offset = 1
        
        # ========== ROW 2: LFP Spectrogram ==========
        ax_spec = fig.add_subplot(gs[row_offset, col_idx])
        spec_axes.append(ax_spec)
        spec_power = ensure_2d(np.asarray(data['spec_power']))
        
        im_spec = ax_spec.pcolormesh(spec_time, spec_freq, spec_power,
                                      shading='gouraud', cmap=CMAP_SPECTROGRAM,
                                      vmin=spec_vmin, vmax=spec_vmax)
        
        ax_spec.set_ylabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_spec.set_yscale('log')
        ax_spec.set_ylim([FREQ_MIN, FREQ_MAX])
        ax_spec.set_xlim([t_min, t_max])
        ax_spec.set_xticklabels([])
        
        # Log-scale frequency ticks
        ax_spec.yaxis.set_major_locator(LogLocator(base=10, subs=[1.0]))
        ax_spec.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
        ax_spec.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax_spec.yaxis.set_minor_formatter(NullFormatter())
        ax_spec.set_yticks([2, 4, 8, 12, 30, 50, 70])
        ax_spec.set_yticklabels(['2', '4', '8', '12', '30', '50', '70'])
        
        if not has_motion:
            ax_spec.set_title(f'{trial_label}\nLFP Spectrogram', 
                             fontsize=FONT_SIZE_TITLE, fontweight='bold')
        else:
            ax_spec.set_title('LFP Spectrogram', fontsize=FONT_SIZE_TITLE - 2, pad=10)
        
        style_axis_publication(ax_spec, remove_top_right=False)
        ax_spec.tick_params(labelsize=FONT_SIZE_TICK)
        
        # ========== ROW 3: LFP-GEVI Coherence Heatmap ==========
        ax_coh = fig.add_subplot(gs[row_offset + 1, col_idx])
        coh_axes.append(ax_coh)
        
        # Get coherence matrix based on method
        if method_name == 'mscohere' and 'coh_mscohere' in data:
            coh_matrix = ensure_2d(np.asarray(data['coh_mscohere']))
        elif method_name == 'fieldtrip' and 'coh_fieldtrip' in data:
            coh_matrix = ensure_2d(np.asarray(data['coh_fieldtrip']))
        else:
            for key in ['coh_mscohere', 'coh_fieldtrip']:
                if key in data:
                    coh_matrix = ensure_2d(np.asarray(data[key]))
                    break
        
        # Interpolate coherence to finer frequency grid for sharper display
        coh_interp, freq_fine = _interpolate_coherence_to_fine_freq(
            coh_matrix, freq, time
        )
        
        im_coh = ax_coh.pcolormesh(time, freq_fine, coh_interp,
                                    shading='gouraud', cmap=CMAP_COHERENCE_DARK,
                                    vmin=coh_vmin, vmax=coh_vmax)
        
        ax_coh.set_xlabel('Time (s)', fontsize=FONT_SIZE_LABEL)
        ax_coh.set_ylabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_coh.set_yscale('log')
        ax_coh.set_ylim([FREQ_MIN, FREQ_MAX])
        ax_coh.set_xlim([t_min, t_max])
        ax_coh.set_title('LFP-GEVI Coherence', fontsize=FONT_SIZE_TITLE - 2, pad=10)
        
        # Log-scale frequency ticks (matching spectrogram)
        ax_coh.yaxis.set_major_locator(LogLocator(base=10, subs=[1.0]))
        ax_coh.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10) * 0.1))
        ax_coh.yaxis.set_major_formatter(FormatStrFormatter('%g'))
        ax_coh.yaxis.set_minor_formatter(NullFormatter())
        ax_coh.set_yticks([2, 4, 8, 12, 30, 50, 70])
        ax_coh.set_yticklabels(['2', '4', '8', '12', '30', '50', '70'])
        
        style_axis_publication(ax_coh, remove_top_right=False)
        ax_coh.tick_params(labelsize=FONT_SIZE_TICK)
    
    # Add colorbars for spectrogram and coherence (shared across trials)
    last_idx = num_trials - 1
    
    # Spectrogram colorbar
    cbar_spec_ax = fig.add_axes([0.94, spec_axes[last_idx].get_position().y0, 0.018, 
                                  spec_axes[last_idx].get_position().height])
    cbar_spec = fig.colorbar(im_spec, cax=cbar_spec_ax)
    cbar_spec.set_label('Power (dB/Hz)', fontsize=FONT_SIZE_COLORBAR, fontweight='bold')
    cbar_spec.ax.tick_params(labelsize=FONT_SIZE_COLORBAR)
    
    # Coherence colorbar
    cbar_coh_ax = fig.add_axes([0.94, coh_axes[last_idx].get_position().y0, 0.018, 
                                 coh_axes[last_idx].get_position().height])
    cbar_coh = fig.colorbar(im_coh, cax=cbar_coh_ax)
    cbar_coh.set_label('Coherence', fontsize=FONT_SIZE_COLORBAR, fontweight='bold')
    cbar_coh.ax.tick_params(labelsize=FONT_SIZE_COLORBAR)
    
    # Ensure output directory exists (with long path support)
    output_dir = Path(output_path).parent
    output_dir_long = to_long_path(str(output_dir))
    os.makedirs(output_dir_long, exist_ok=True)
    
    # Save in multiple formats - clean naming: {method}_heatmaps.{fmt}
    for fmt in FIGURE_FORMAT:
        save_path = str(Path(output_path).with_name(f"{Path(output_path).name}_heatmaps.{fmt}"))
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
        print(f"    Saved: {Path(save_path).name}")
    
    plt.close(fig)


# ==============================================================================
#  COHERENCE VS FREQUENCY FIGURE
# ==============================================================================

def create_coherence_spectrum_figure(data_list, trial_labels, method_name, output_path):
    """
    Create coherence vs frequency figure with overall and rest/run breakdown.
    
    Layout: 2 rows × N columns (dynamic based on number of trials)
        Row 1: Overall coherence vs frequency (line plot, no fill)
        Row 2: Rest (solid) vs Run (dashed) coherence
    
    Parameters
    ----------
    data_list : list of dict
        List of trial data dictionaries loaded from MATLAB .mat files
    trial_labels : list of str
        Labels for each trial
    method_name : str
        Analysis method name ('mscohere' or 'fieldtrip')
    output_path : Path
        Base path for saving output files
    
    Publication styling:
        - Colored dotted lines for frequency band boundaries (no shading)
        - No fill under curves - line plots only
        - Teal color shades for rest vs run
        - Stats text positioned to avoid overlap
    """
    num_trials = len(data_list)
    
    # Dynamic figure sizing - scale with number of trials
    fig_width = max(10, 6 * num_trials)
    fig_width = min(fig_width, 36)  # Cap at reasonable size
    
    fig, axes = plt.subplots(2, num_trials, figsize=(fig_width, 14))
    fig.subplots_adjust(left=0.08, right=0.92, top=0.90, bottom=0.08, 
                        hspace=0.50, wspace=0.25)  # Good spacing without subplot titles
    
    # Handle single trial case (axes won't be 2D)
    if num_trials == 1:
        axes = axes.reshape(2, 1)
    
    for col_idx, (data, trial_label) in enumerate(zip(data_list, trial_labels)):
        freq = np.asarray(data['freq']).flatten()
        coh_spectrum = np.asarray(data['coh_spectrum']).flatten()
        
        # ========== ROW 1: Overall coherence ==========
        ax_overall = axes[0, col_idx] if num_trials > 1 else axes[0]
        
        # Add colored dotted lines for frequency band boundaries
        add_frequency_band_lines(ax_overall, y_max=1.0, labels_inside=True)
        
        # Plot overall coherence as LINE plot (no fill)
        ax_overall.plot(freq, coh_spectrum, '-', color=COLOR_OVERALL, 
                       linewidth=LINE_WIDTH_TRACE)
        
        ax_overall.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_overall.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
        # No subplot title - info is in stats text
        ax_overall.set_xlim([FREQ_MIN, FREQ_MAX])
        ax_overall.set_ylim([0, 1])
        
        # Apply publication styling
        style_axis_publication(ax_overall)
        
        # Add statistics - positioned above plot (no box), consistent with PSD figure
        peak_idx = np.argmax(coh_spectrum)
        peak_coh = coh_spectrum[peak_idx]
        peak_freq = freq[peak_idx]
        theta_mask = (freq >= 4) & (freq <= 8)
        mean_theta = np.mean(coh_spectrum[theta_mask]) if np.any(theta_mask) else 0
        
        stats_text = f'Peak: {peak_coh:.2f} @ {peak_freq:.1f} Hz | θ-band: {mean_theta:.2f}'
        ax_overall.text(0.5, 1.06, stats_text, transform=ax_overall.transAxes,
                       ha='center', va='bottom', fontsize=FONT_SIZE_STATS,
                       fontweight='bold', color=COLOR_OVERALL)
        
        # ========== ROW 2: Rest vs Run coherence ==========
        ax_restrun = axes[1, col_idx] if num_trials > 1 else axes[1]
        
        # Add colored dotted lines for frequency band boundaries
        add_frequency_band_lines(ax_restrun, y_max=1.0, labels_inside=True)
        
        legend_handles = []
        legend_labels = []
        stats_lines = []  # For peak coherence stats
        
        # Plot REST coherence (solid line - darker teal)
        has_rest = 'coh_rest' in data and data['coh_rest'] is not None
        if has_rest:
            coh_rest = np.asarray(data['coh_rest']).flatten()
            freq_rest = np.asarray(data['coh_rest_freq']).flatten() if 'coh_rest_freq' in data else freq
            if len(coh_rest) > 0:
                h_rest, = ax_restrun.plot(freq_rest, coh_rest, '-', 
                                          color=COLOR_REST, linewidth=LINE_WIDTH_TRACE)
                pct_rest = data.get('pct_rest', 0)
                legend_handles.append(h_rest)
                legend_labels.append(f'Rest ({pct_rest:.0f}%)')
                
                # Get peak for rest
                peak_idx_rest = np.argmax(coh_rest)
                peak_coh_rest = coh_rest[peak_idx_rest]
                peak_freq_rest = freq_rest[peak_idx_rest]
                stats_lines.append((f'Rest peak: {peak_coh_rest:.2f} @ {peak_freq_rest:.1f} Hz', COLOR_REST))
        
        # Plot RUN coherence (dashed line - lighter teal)
        has_run = 'coh_run' in data and data['coh_run'] is not None
        if has_run:
            coh_run = np.asarray(data['coh_run']).flatten()
            freq_run = np.asarray(data['coh_run_freq']).flatten() if 'coh_run_freq' in data else freq
            if len(coh_run) > 0:
                h_run, = ax_restrun.plot(freq_run, coh_run, '--', 
                                         color=COLOR_RUN, linewidth=LINE_WIDTH_DASHED)
                pct_run = data.get('pct_run', 0)
                legend_handles.append(h_run)
                legend_labels.append(f'Run ({pct_run:.0f}%)')
                
                # Get peak for run
                peak_idx_run = np.argmax(coh_run)
                peak_coh_run = coh_run[peak_idx_run]
                peak_freq_run = freq_run[peak_idx_run]
                stats_lines.append((f'Run peak: {peak_coh_run:.2f} @ {peak_freq_run:.1f} Hz', COLOR_RUN))
        
        # Add legend - positioned INSIDE the plot (upper right) to avoid column overlap
        if legend_handles:
            leg = ax_restrun.legend(legend_handles, legend_labels, 
                                    loc='upper right',
                                    fontsize=FONT_SIZE_LEGEND, frameon=False)
        
        ax_restrun.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_restrun.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
        # No subplot title - info is in stats text
        ax_restrun.set_xlim([FREQ_MIN, FREQ_MAX])
        ax_restrun.set_ylim([0, 1])
        
        # Apply publication styling
        style_axis_publication(ax_restrun)
        
        # Add peak coherence stats above plot (no subplot label, consistent with PSD figure)
        if stats_lines:
            # Combine rest and run peaks into one line
            stats_text = ' | '.join([s[0].replace(' peak', ' pk') for s in stats_lines])
            ax_restrun.text(0.5, 1.06, stats_text, 
                           transform=ax_restrun.transAxes,
                           ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                           fontweight='bold', color=COLOR_OVERALL)
    
    # Main title - positioned at very top
    fig.suptitle(f'LFP-GEVI Coherence Spectrum ({method_name.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.99)
    
    # Ensure output directory exists (with long path support)
    output_dir = Path(output_path).parent
    output_dir_long = to_long_path(str(output_dir))
    os.makedirs(output_dir_long, exist_ok=True)
    
    # Save in multiple formats - clean naming: {method}_coherence.{fmt}
    for fmt in FIGURE_FORMAT:
        save_path = f"{output_path}_coherence.{fmt}"
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
        print(f"    Saved: {Path(save_path).name}")
    
    plt.close(fig)


# ==============================================================================
#  PSD FIGURE (Power Spectral Density)
# ==============================================================================

def create_psd_figure(data_list, trial_labels, method_name, output_path):
    """
    Create Power Spectral Density (PSD) figure for LFP and GEVI signals.
    
    Layout: 2 rows × (2 × num_trials) columns (dynamic)
        Columns: [LFP Trial1] [GEVI Trial1] [LFP Trial2] [GEVI Trial2] ...
        Row 1: Overall PSD
        Row 2: Rest (solid) vs Run (dashed) PSD
    
    Parameters
    ----------
    data_list : list of dict
        List of trial data dictionaries loaded from MATLAB .mat files
    trial_labels : list of str
        Labels for each trial
    method_name : str
        Analysis method name ('mscohere' or 'fieldtrip')
    output_path : Path
        Base path for saving output files
    
    Colors: 
        - LFP: purple-grey (matching Figure 1)
        - GEVI: teal (matching Figure 1)
        - Rest: solid line
        - Run: dashed line
    """
    num_trials = len(data_list)
    
    # Dynamic figure sizing - scale with number of trials
    fig_width = max(12, 6 * num_trials)
    fig_width = min(fig_width, 36)  # Cap at reasonable size
    
    # Create figure with 2 columns per trial (LFP, GEVI)
    fig, axes = plt.subplots(2, num_trials * 2, figsize=(fig_width, 14))
    fig.subplots_adjust(left=0.06, right=0.92, top=0.88, bottom=0.06, 
                        hspace=0.40, wspace=0.50)  # Dynamic spacing
    
    # Handle single trial case
    if num_trials == 1:
        axes = axes.reshape(2, 2)
    
    for trial_idx, (data, trial_label) in enumerate(zip(data_list, trial_labels)):
        # Check if PSD data is available
        if 'psd_freq' not in data or data['psd_freq'] is None:
            print(f"  Warning: No PSD data available for {trial_label}")
            continue
        
        freq = np.asarray(data['psd_freq']).flatten()
        psd_lfp = np.asarray(data.get('psd_lfp', [])).flatten()
        psd_gevi = np.asarray(data.get('psd_gevi', [])).flatten()
        
        # Column indices for this trial
        col_lfp = trial_idx * 2      # LFP column
        col_gevi = trial_idx * 2 + 1  # GEVI column
        
        # ========== ROW 1: Overall PSD ==========
        # --- LFP Overall ---
        ax_lfp_overall = axes[0, col_lfp]
        
        if len(psd_lfp) > 0:
            ax_lfp_overall.plot(freq, psd_lfp, '-', color=COLOR_LFP, 
                               linewidth=LINE_WIDTH_TRACE, label='LFP')
            
            # Stats text (no box)
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
        # Set y-axis ticks to increments of 5, no decimals
        ax_lfp_overall.yaxis.set_major_locator(MultipleLocator(5))
        ax_lfp_overall.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # Add band lines ABOVE plot (labels_inside=False to avoid overlap with traces)
        add_frequency_band_lines(ax_lfp_overall, add_labels=True, labels_inside=False)
        
        # --- GEVI Overall ---
        ax_gevi_overall = axes[0, col_gevi]
        
        if len(psd_gevi) > 0:
            ax_gevi_overall.plot(freq, psd_gevi, '-', color=COLOR_GEVI, 
                                linewidth=LINE_WIDTH_TRACE, label='GEVI')
            
            # Stats text (no box)
            peak_idx = np.argmax(psd_gevi)
            peak_val = psd_gevi[peak_idx]
            peak_freq = freq[peak_idx]
            
            stats_text = f'Peak: {peak_val:.1f} dB @ {peak_freq:.1f} Hz'
            ax_gevi_overall.text(0.5, 1.08, stats_text, transform=ax_gevi_overall.transAxes,
                                ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 1,
                                fontweight='bold', color=COLOR_GEVI)
        
        ax_gevi_overall.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_gevi_overall.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
        ax_gevi_overall.set_xlim([FREQ_MIN, FREQ_MAX])
        style_axis_publication(ax_gevi_overall)
        # Set y-axis ticks to increments of 5, no decimals
        ax_gevi_overall.yaxis.set_major_locator(MultipleLocator(5))
        ax_gevi_overall.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # Add band lines ABOVE plot
        add_frequency_band_lines(ax_gevi_overall, add_labels=True, labels_inside=False)
        
        # ========== ROW 2: Rest vs Run PSD ==========
        # --- LFP Rest/Run ---
        ax_lfp_rr = axes[1, col_lfp]
        
        legend_handles_lfp = []
        legend_labels_lfp = []
        
        # LFP Rest
        lfp_rest_peak = None
        lfp_run_peak = None
        has_lfp_rest = 'psd_lfp_rest' in data and data['psd_lfp_rest'] is not None
        if has_lfp_rest:
            psd_lfp_rest = np.asarray(data['psd_lfp_rest']).flatten()
            freq_rest = np.asarray(data.get('psd_rest_freq', freq)).flatten()
            if len(psd_lfp_rest) > 0:
                h_rest, = ax_lfp_rr.plot(freq_rest, psd_lfp_rest, '-', 
                                        color=COLOR_LFP, linewidth=LINE_WIDTH_TRACE,
                                        alpha=0.9)
                pct_rest = data.get('pct_rest', 0)
                legend_handles_lfp.append(h_rest)
                legend_labels_lfp.append(f'Rest ({pct_rest:.0f}%)')
                # Get peak for rest
                peak_idx = np.argmax(psd_lfp_rest)
                lfp_rest_peak = (psd_lfp_rest[peak_idx], freq_rest[peak_idx])
        
        # LFP Run
        has_lfp_run = 'psd_lfp_run' in data and data['psd_lfp_run'] is not None
        if has_lfp_run:
            psd_lfp_run = np.asarray(data['psd_lfp_run']).flatten()
            freq_run = np.asarray(data.get('psd_run_freq', freq)).flatten()
            if len(psd_lfp_run) > 0:
                h_run, = ax_lfp_rr.plot(freq_run, psd_lfp_run, '--', 
                                       color=COLOR_LFP, linewidth=LINE_WIDTH_DASHED,
                                       alpha=0.7)
                pct_run = data.get('pct_run', 0)
                legend_handles_lfp.append(h_run)
                legend_labels_lfp.append(f'Run ({pct_run:.0f}%)')
                # Get peak for run
                peak_idx = np.argmax(psd_lfp_run)
                lfp_run_peak = (psd_lfp_run[peak_idx], freq_run[peak_idx])
        
        # Legend positioned INSIDE the plot (upper right)
        if legend_handles_lfp:
            ax_lfp_rr.legend(legend_handles_lfp, legend_labels_lfp, 
                            loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
        
        ax_lfp_rr.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_lfp_rr.set_ylabel('LFP Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
        ax_lfp_rr.set_xlim([FREQ_MIN, FREQ_MAX])
        style_axis_publication(ax_lfp_rr)
        # Set y-axis ticks to increments of 5, no decimals
        ax_lfp_rr.yaxis.set_major_locator(MultipleLocator(5))
        ax_lfp_rr.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # Add band lines ABOVE plot (labels_inside=False)
        add_frequency_band_lines(ax_lfp_rr, add_labels=True, labels_inside=False)
        
        # Add peak stats text ABOVE band names (y=1.08)
        if lfp_rest_peak and lfp_run_peak:
            stats_text = f'Rest pk: {lfp_rest_peak[0]:.0f} dB | Run pk: {lfp_run_peak[0]:.0f} dB'
            ax_lfp_rr.text(0.5, 1.08, stats_text, transform=ax_lfp_rr.transAxes,
                          ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 2,
                          fontweight='bold', color=COLOR_LFP)
        
        # --- GEVI Rest/Run ---
        ax_gevi_rr = axes[1, col_gevi]
        
        legend_handles_gevi = []
        legend_labels_gevi = []
        gevi_rest_peak = None
        gevi_run_peak = None
        
        # GEVI Rest
        has_gevi_rest = 'psd_gevi_rest' in data and data['psd_gevi_rest'] is not None
        if has_gevi_rest:
            psd_gevi_rest = np.asarray(data['psd_gevi_rest']).flatten()
            freq_rest = np.asarray(data.get('psd_rest_freq', freq)).flatten()
            if len(psd_gevi_rest) > 0:
                h_rest, = ax_gevi_rr.plot(freq_rest, psd_gevi_rest, '-', 
                                         color=COLOR_GEVI, linewidth=LINE_WIDTH_TRACE,
                                         alpha=0.9)
                pct_rest = data.get('pct_rest', 0)
                legend_handles_gevi.append(h_rest)
                legend_labels_gevi.append(f'Rest ({pct_rest:.0f}%)')
                # Get peak for rest
                peak_idx = np.argmax(psd_gevi_rest)
                gevi_rest_peak = (psd_gevi_rest[peak_idx], freq_rest[peak_idx])
        
        # GEVI Run
        has_gevi_run = 'psd_gevi_run' in data and data['psd_gevi_run'] is not None
        if has_gevi_run:
            psd_gevi_run = np.asarray(data['psd_gevi_run']).flatten()
            freq_run = np.asarray(data.get('psd_run_freq', freq)).flatten()
            if len(psd_gevi_run) > 0:
                h_run, = ax_gevi_rr.plot(freq_run, psd_gevi_run, '--', 
                                        color=COLOR_GEVI, linewidth=LINE_WIDTH_DASHED,
                                        alpha=0.7)
                pct_run = data.get('pct_run', 0)
                legend_handles_gevi.append(h_run)
                legend_labels_gevi.append(f'Run ({pct_run:.0f}%)')
                # Get peak for run
                peak_idx = np.argmax(psd_gevi_run)
                gevi_run_peak = (psd_gevi_run[peak_idx], freq_run[peak_idx])
        
        # Legend positioned INSIDE the plot (upper right)
        if legend_handles_gevi:
            ax_gevi_rr.legend(legend_handles_gevi, legend_labels_gevi, 
                             loc='upper right', fontsize=FONT_SIZE_LEGEND - 1, frameon=False)
        
        ax_gevi_rr.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax_gevi_rr.set_ylabel('GEVI Power (dB/Hz)', fontsize=FONT_SIZE_LABEL)
        ax_gevi_rr.set_xlim([FREQ_MIN, FREQ_MAX])
        style_axis_publication(ax_gevi_rr)
        # Set y-axis ticks to increments of 5, no decimals
        ax_gevi_rr.yaxis.set_major_locator(MultipleLocator(5))
        ax_gevi_rr.yaxis.set_major_formatter(FormatStrFormatter('%d'))
        # Add band lines ABOVE plot (labels_inside=False)
        add_frequency_band_lines(ax_gevi_rr, add_labels=True, labels_inside=False)
        
        # Add peak stats text ABOVE band names (y=1.08)
        if gevi_rest_peak and gevi_run_peak:
            stats_text = f'Rest pk: {gevi_rest_peak[0]:.0f} dB | Run pk: {gevi_run_peak[0]:.0f} dB'
            ax_gevi_rr.text(0.5, 1.08, stats_text, transform=ax_gevi_rr.transAxes,
                           ha='center', va='bottom', fontsize=FONT_SIZE_STATS - 2,
                           fontweight='bold', color=COLOR_GEVI)
    
    # Main title - positioned at very top (above all peak stats at y=1.08)
    fig.suptitle(f'Power Spectral Density – LFP and GEVI ({method_name.upper()} Method)',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.95)
    
    # Ensure output directory exists (with long path support)
    output_dir = Path(output_path).parent
    output_dir_long = to_long_path(str(output_dir))
    os.makedirs(output_dir_long, exist_ok=True)
    
    # Save in multiple formats
    for fmt in FIGURE_FORMAT:
        save_path = f"{output_path}_psd.{fmt}"
        save_path_long = to_long_path(save_path)
        fig.savefig(save_path_long, dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
        print(f"    Saved: {Path(save_path).name}")
    
    plt.close(fig)


# ==============================================================================
#  HELPER: Discover Trial Files
# ==============================================================================

def discover_trial_files(input_dir, method):
    """
    Discover all available trial files for a given method.
    
    New naming convention: {method}_trial{N}.mat or {method}.mat (combined)
    
    Parameters
    ----------
    input_dir : Path
        Directory containing the trial .mat files (data/ folder)
    method : str
        Analysis method ('mscohere' or 'fieldtrip')
    
    Returns
    -------
    list of Path
        Sorted list of trial file paths
    """
    # Try new naming convention first
    pattern_new = f"{method}_trial*.mat"
    trial_files = sorted(input_dir.glob(pattern_new))
    
    # Fallback to old naming convention for backwards compatibility
    if not trial_files:
        pattern_old = f"figure2_{method}_trial*.mat"
        trial_files = sorted(input_dir.glob(pattern_old))
    
    # Also check for combined file (single mat with all trials)
    if not trial_files:
        combined_file = input_dir / f"{method}.mat"
        if combined_file.exists():
            trial_files = [combined_file]
    
    return trial_files


def generate_trial_labels(num_trials, custom_labels=None):
    """
    Generate trial labels for plotting.
    
    Parameters
    ----------
    num_trials : int
        Number of trials
    custom_labels : list or None
        Custom labels to use (if provided and matching length)
    
    Returns
    -------
    list of str
        Trial labels
    """
    if custom_labels is not None and len(custom_labels) >= num_trials:
        return custom_labels[:num_trials]
    
    # Auto-generate labels
    return [f'Trial {i+1}' for i in range(num_trials)]


# ==============================================================================
#  MAIN EXECUTION
# ==============================================================================

def main():
    """Main execution function."""
    
    print("=" * 70)
    print("FIGURE 2: LFP-GEVI Coherence Analysis - Python Plotting")
    print("=" * 70)
    
    # Print mouse and session info
    print(f"\nMouse ID:   {MOUSE_ID}")
    print(f"Session ID: {SESSION_ID}")
    print(f"Input Dir:  {INPUT_DIR}")
    
    # Check if input directory exists
    if not INPUT_DIR.exists():
        print(f"\nERROR: Input directory does not exist: {INPUT_DIR}")
        print("  Please verify MOUSE_ID and SESSION_ID are correct.")
        print(f"  Available subdirectories in {BASE_INPUT_DIR}:")
        if BASE_INPUT_DIR.exists():
            for subdir in sorted(BASE_INPUT_DIR.iterdir()):
                if subdir.is_dir():
                    print(f"    - {subdir.name}")
        return
    
    # Create output directory if it doesn't exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Output Dir: {OUTPUT_DIR}")
    
    # Process both methods
    for method in ['mscohere', 'fieldtrip']:
        print(f"\n{'='*50}")
        print(f"Processing {method.upper()} method...")
        print('='*50)
        
        # Discover all available trial files
        trial_files = discover_trial_files(INPUT_DIR, method)
        
        if len(trial_files) == 0:
            print(f"  WARNING: No trial files found for {method} method. Skipping.")
            continue
        
        print(f"  Found {len(trial_files)} trial file(s):")
        for tf in trial_files:
            print(f"    - {tf.name}")
        
        # Load all trial data
        data_list = []
        for trial_path in trial_files:
            print(f"  Loading: {trial_path.name}")
            data = load_matlab_data(trial_path)
            if data is None:
                print(f"  WARNING: Failed to load {trial_path.name}. Skipping.")
                continue
            data_list.append(data)

        if len(data_list) == 0:
            print(f"  WARNING: No trial data could be loaded for {method} method. Skipping.")
            continue

        # Generate trial labels
        trial_labels = generate_trial_labels(len(data_list), TRIAL_LABELS)
        
        # Print data summary for first trial
        data_trial1 = data_list[0]
        print(f"\n  Data summary (Trial 1 of {len(data_list)}):")
        print(f"    Time: {data_trial1['time'].shape} ({data_trial1['time'].min():.1f} - {data_trial1['time'].max():.1f} s)")
        print(f"    Freq: {data_trial1['freq'].shape} ({data_trial1['freq'].min():.1f} - {data_trial1['freq'].max():.1f} Hz)")
        print(f"    Spectrogram: {data_trial1['spec_power'].shape}")
        
        # Rest/Run info (if available)
        if 'pct_rest' in data_trial1:
            print(f"    Behavioral states: Rest {data_trial1['pct_rest']:.1f}%, Run {data_trial1['pct_run']:.1f}%")
        if 'coh_rest' in data_trial1 and data_trial1['coh_rest'] is not None:
            coh_rest = np.asarray(data_trial1['coh_rest']).flatten()
            print(f"    Rest coherence: {len(coh_rest)} freq points")
        if 'coh_run' in data_trial1 and data_trial1['coh_run'] is not None:
            coh_run = np.asarray(data_trial1['coh_run']).flatten()
            print(f"    Run coherence: {len(coh_run)} freq points")
        
        # PSD info (if available)
        if 'psd_freq' in data_trial1 and data_trial1['psd_freq'] is not None:
            psd_freq = np.asarray(data_trial1['psd_freq']).flatten()
            print(f"    PSD: {len(psd_freq)} freq points")
            if 'psd_lfp' in data_trial1:
                psd_lfp = np.asarray(data_trial1['psd_lfp']).flatten()
                print(f"    LFP PSD range: {psd_lfp.min():.1f} to {psd_lfp.max():.1f} dB")
            if 'psd_gevi' in data_trial1:
                psd_gevi = np.asarray(data_trial1['psd_gevi']).flatten()
                print(f"    GEVI PSD range: {psd_gevi.min():.1f} to {psd_gevi.max():.1f} dB")
        
        # Output path base - cleaner naming
        output_base = OUTPUT_DIR / method
        
        # Create figures
        print(f"\n  Generating figures for {len(data_list)} trials...")
        
        # Main figure (speed heatmap + spectrogram + coherence heatmaps)
        # Output: {method}_heatmaps.{fmt}
        create_main_figure(data_list, trial_labels, method, output_base)
        
        # Coherence spectrum figure (overall + rest vs run)
        # Output: {method}_coherence.{fmt}
        create_coherence_spectrum_figure(data_list, trial_labels, method, output_base)
        
        # PSD figure (LFP and GEVI power spectra)
        # Output: {method}_psd.{fmt}
        if 'psd_freq' in data_trial1 and data_trial1['psd_freq'] is not None:
            create_psd_figure(data_list, trial_labels, method, output_base)
        else:
            print("  Note: PSD data not available - run updated MATLAB scripts to generate")
    
    print("\n" + "=" * 70)
    print("FIGURE 2 PLOTTING COMPLETE")
    print("=" * 70)
    print(f"\nAll figures saved to: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
