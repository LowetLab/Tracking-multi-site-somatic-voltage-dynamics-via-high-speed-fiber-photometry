"""
================================================================================
FIGURE 5: Theta Band Analysis - Coherence & Power
================================================================================

Publication-ready figure showing theta band (5-9 Hz) comparisons for REST vs RUN:
  - LFP-GEVI Coherence in theta band
  - LFP Power in theta band
  - GEVI Power in theta band

Visualization: Half-violin plots with individual animal dots and paired lines.
Statistics: One-tailed paired t-test (directional hypothesis: RUN > REST)

================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.io import loadmat
from scipy.stats import wilcoxon, ttest_rel, shapiro, ttest_1samp
from pathlib import Path
import warnings
import sys

warnings.filterwarnings('ignore', category=UserWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# common.py / plotting_config.py live in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# ==============================================================================
#  CONFIGURATION
# ==============================================================================

try:
    from plotting_config import (
        BEHAVIOR_MODE, FIGURE_DPI, FIGURE_FORMATS, 
        get_animal_concatenated_input_dir, get_animal_pooled_input_dir,
        get_group_level_output_dir, GROUP_POOLING_LEVEL, get_animals_to_process,
        FONT_SIZE_SUPTITLE, FONT_SIZE_TITLE, FONT_SIZE_LABEL, FONT_SIZE_TICK,
        FONT_SIZE_LEGEND, AXIS_LINEWIDTH, TICK_WIDTH, TICK_LENGTH,
        COLOR_LFP_REST, COLOR_LFP_RUN, COLOR_GEVI_REST, COLOR_GEVI_RUN,
        COLOR_COH_REST, COLOR_COH_RUN,
    )
    USING_CENTRAL_CONFIG = True
except ImportError:
    USING_CENTRAL_CONFIG = False
    BEHAVIOR_MODE = 'clear'
    FIGURE_DPI = 300
    FIGURE_FORMATS = ['png', 'pdf', 'svg']
    FONT_SIZE_SUPTITLE = 20
    FONT_SIZE_TITLE = 16
    FONT_SIZE_LABEL = 14
    FONT_SIZE_TICK = 12
    FONT_SIZE_LEGEND = 12
    AXIS_LINEWIDTH = 2.0
    TICK_WIDTH = 1.8
    TICK_LENGTH = 7
    COLOR_LFP_REST = np.array([0.25, 0.18, 0.35])
    COLOR_LFP_RUN = np.array([0.55, 0.45, 0.65])
    COLOR_GEVI_REST = np.array([0.05, 0.35, 0.45])
    COLOR_GEVI_RUN = np.array([0.25, 0.65, 0.65])
    COLOR_COH_REST = np.array([0.05, 0.35, 0.45])
    COLOR_COH_RUN = np.array([0.25, 0.65, 0.65])
    GROUP_POOLING_LEVEL = 'animal_concatenated'
    
    def get_animals_to_process():
        return []

# Theta band definition
THETA_RANGE = [5, 9]  # Hz

# Spectral method to use: 'mscohere' or 'fieldtrip'
# IMPORTANT: Make sure this matches the method used for your visual coherence plots!
SPECTRAL_METHOD = 'fieldtrip'  # Changed from 'mscohere' to match FIELDTRIP plots

# Plot parameters (same as FOOOF)
VIOLIN_WIDTH = 0.28
BOX_WIDTH = 0.10
BOX_OFFSET = 0.08
DOT_OFFSET = 0.18
DOT_SIZE = 90
LINE_WIDTH = 1.5
LINE_ALPHA = 0.5

# Animal labels
SHOW_ANIMAL_LABELS = False

# ==============================================================================
#  WINDOWS LONG PATH SUPPORT
# ==============================================================================

from common import to_long_path  # shared helper (was a local copy)


# ==============================================================================
#  DATA LOADING
# ==============================================================================

def get_pooling_level_label(pooling_level):
    if pooling_level == 'animal_concatenated':
        return 'Animal-Concatenated'
    elif pooling_level == 'animal_pooled':
        return 'Animal-Pooled'
    return pooling_level.replace('_', '-').title()


# Method for extracting theta band value
# Options: 'peak' (maximum value), 'mean' (average), 'auc' (area under curve)
THETA_EXTRACTION_METHOD = 'peak'  # Changed to 'peak' for better oscillation detection

# Alpha band for theta/alpha ratio
ALPHA_RANGE = [9, 13]  # Hz


def compute_1f_deviation(freq, spectrum, theta_range=THETA_RANGE):
    """
    Compute how much theta deviates ABOVE the expected 1/f background.
    
    Method:
    1. Fit a linear trend (in log-freq space) to flanking bands (2-4 Hz and 15-25 Hz)
    2. Predict expected theta power from this 1/f trend
    3. Return actual_theta_peak - expected_theta
    
    This captures the "bump" in RUN that rises above the 1/f slope,
    vs the monotonic decrease in REST that follows the 1/f trend.
    
    Returns: deviation in dB (positive = theta bumps above 1/f trend)
    """
    freq = np.asarray(freq).flatten()
    spectrum = np.asarray(spectrum).flatten()
    
    # Get theta peak
    theta_mask = (freq >= theta_range[0]) & (freq <= theta_range[1])
    if not np.any(theta_mask):
        return np.nan
    
    theta_peak = np.nanmax(spectrum[theta_mask])
    theta_peak_freq = freq[theta_mask][np.nanargmax(spectrum[theta_mask])]
    
    # Get flanking bands for fitting (EXCLUDE theta and alpha to get clean 1/f)
    low_flank_mask = (freq >= 2) & (freq <= 4)
    high_flank_mask = (freq >= 15) & (freq <= 25)
    fit_mask = low_flank_mask | high_flank_mask
    
    if np.sum(fit_mask) < 4:
        return np.nan
    
    # Fit linear trend in log-frequency space (1/f appears linear in log-log)
    fit_freq = freq[fit_mask]
    fit_power = spectrum[fit_mask]
    
    # Use log(freq) for x-axis (1/f is linear in log-freq vs power in dB)
    log_freq = np.log10(fit_freq)
    
    # Linear fit: power = slope * log(freq) + intercept
    try:
        coeffs = np.polyfit(log_freq, fit_power, 1)
        slope, intercept = coeffs
        
        # Predict expected power at theta peak frequency
        expected_theta = slope * np.log10(theta_peak_freq) + intercept
        
        # Deviation = actual - expected (positive = bump above trend)
        deviation = theta_peak - expected_theta
        
        return deviation
    except:
        return np.nan


def compute_local_peak_index(freq, spectrum, theta_range=THETA_RANGE):
    """
    Compute theta spectral contrast: peak-to-trough difference.
    
    Method:
    Compare theta peak (5-9 Hz) to minimum power in post-theta trough (12-20 Hz).
    A distinct oscillation will show theta_peak >> trough_min.
    A monotonic 1/f decline will show smaller difference.
    
    Using minimum (not mean) in the trough band to capture the deepest point
    of the spectral "valley" after theta, which maximizes contrast detection.
    
    Returns: theta_peak - trough_min (in dB)
    """
    freq = np.asarray(freq).flatten()
    spectrum = np.asarray(spectrum).flatten()
    
    # Get theta peak
    theta_mask = (freq >= theta_range[0]) & (freq <= theta_range[1])
    if not np.any(theta_mask):
        return np.nan
    theta_peak = np.nanmax(spectrum[theta_mask])
    
    # Get post-theta trough (spectral valley after theta - widened range, use MIN)
    trough_mask = (freq >= 12) & (freq <= 20)
    if not np.any(trough_mask):
        return np.nan
    trough_min = np.nanmin(spectrum[trough_mask])  # Use MIN to capture deepest point
    
    # Spectral contrast = peak minus trough
    # Higher = more distinct peak, Lower = monotonic decline
    return theta_peak - trough_min


def compute_theta_band_value(freq, spectrum, theta_range=THETA_RANGE, method=None):
    """
    Compute theta band value using specified method.
    
    Methods:
        'peak': Maximum value in theta band (best for oscillation detection)
        'mean': Average across theta band (smooths over frequency variations)
        'auc': Area under curve (integrates total power/coherence)
    """
    if method is None:
        method = THETA_EXTRACTION_METHOD
    
    freq = np.asarray(freq).flatten()
    spectrum = np.asarray(spectrum).flatten()
    
    # Find theta frequency indices
    theta_mask = (freq >= theta_range[0]) & (freq <= theta_range[1])
    
    if not np.any(theta_mask):
        return np.nan
    
    theta_values = spectrum[theta_mask]
    
    if method == 'peak':
        return np.nanmax(theta_values)
    elif method == 'mean':
        return np.nanmean(theta_values)
    elif method == 'auc':
        # Area under curve (sum * frequency resolution)
        freq_res = np.mean(np.diff(freq[theta_mask])) if np.sum(theta_mask) > 1 else 1
        return np.nansum(theta_values) * freq_res
    else:
        return np.nanmean(theta_values)


def load_theta_band_data(animals_to_process, method='mscohere', pooling_level=None):
    """
    Load coherence and PSD data, extract theta band averages.
    
    Returns dict with:
        - animal_ids: list of animal IDs
        - coh_theta_rest, coh_theta_run: Coherence in theta band
        - lfp_theta_rest, lfp_theta_run: LFP power in theta band
        - gevi_theta_rest, gevi_theta_run: GEVI power in theta band
    """
    if pooling_level is None:
        pooling_level = GROUP_POOLING_LEVEL
    
    from plotting_config import get_base_dir
    base_path = get_base_dir() / pooling_level
    
    results = {
        'animal_ids': [],
        'coh_theta_rest': [],
        'coh_theta_run': [],
        'lfp_theta_rest': [],
        'lfp_theta_run': [],
        'gevi_theta_rest': [],
        'gevi_theta_run': [],
        # New metrics for peakiness
        'gevi_prominence_rest': [],
        'gevi_prominence_run': [],
        'gevi_theta_alpha_rest': [],
        'gevi_theta_alpha_run': [],
    }
    
    print(f"    Loading theta band data from {pooling_level}...")
    print(f"    Method: {method}")
    print(f"    Theta range: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz")
    print(f"    Extraction method: {THETA_EXTRACTION_METHOD} (peak=max value, mean=average)")
    
    for animal in animals_to_process:
        mouse_id = animal['mouse_id']
        data_path = base_path / mouse_id / 'data' / f'{method}.mat'
        data_path_long = to_long_path(str(data_path))
        
        if not os.path.exists(data_path_long):
            print(f"      [{mouse_id}] SKIP: File not found")
            continue
        
        try:
            mat = loadmat(data_path_long, squeeze_me=True)
            
            # Get frequency axis
            freq = mat.get('freq', mat.get('psd_freq', None))
            if freq is None:
                print(f"      [{mouse_id}] SKIP: No frequency axis")
                continue
            freq = np.asarray(freq).flatten()
            
            # Extract coherence (MATLAB saves as coh_rest, coh_run)
            coh_rest = mat.get('coh_rest', None)
            coh_run = mat.get('coh_run', None)
            
            # Extract PSD
            psd_lfp_rest = mat.get('psd_lfp_rest', None)
            psd_lfp_run = mat.get('psd_lfp_run', None)
            psd_gevi_rest = mat.get('psd_gevi_rest', None)
            psd_gevi_run = mat.get('psd_gevi_run', None)
            
            # Debug: show available keys
            available_keys = [k for k in mat.keys() if not k.startswith('_')]
            
            # Check required data - handle empty arrays from MATLAB
            def is_valid_array(arr):
                if arr is None:
                    return False
                arr = np.asarray(arr).flatten()
                return len(arr) > 0 and not np.all(np.isnan(arr))
            
            if not is_valid_array(coh_rest) or not is_valid_array(coh_run):
                print(f"      [{mouse_id}] SKIP: Missing/empty coherence data (coh_rest/coh_run)")
                print(f"                    Available keys: {available_keys}")
                continue
            if not is_valid_array(psd_lfp_rest) or not is_valid_array(psd_lfp_run):
                print(f"      [{mouse_id}] SKIP: Missing/empty LFP PSD data")
                continue
            
            # Flatten arrays
            coh_rest = np.asarray(coh_rest).flatten()
            coh_run = np.asarray(coh_run).flatten()
            psd_lfp_rest = np.asarray(psd_lfp_rest).flatten()
            psd_lfp_run = np.asarray(psd_lfp_run).flatten()
            
            # Compute theta band values (using THETA_EXTRACTION_METHOD: peak/mean/auc)
            coh_theta_r = compute_theta_band_value(freq, coh_rest)
            coh_theta_n = compute_theta_band_value(freq, coh_run)
            lfp_theta_r = compute_theta_band_value(freq, psd_lfp_rest)
            lfp_theta_n = compute_theta_band_value(freq, psd_lfp_run)
            
            # GEVI is optional
            gevi_theta_r = np.nan
            gevi_theta_n = np.nan
            gevi_prom_r = np.nan
            gevi_prom_n = np.nan
            gevi_ta_r = np.nan
            gevi_ta_n = np.nan
            
            if psd_gevi_rest is not None and psd_gevi_run is not None:
                psd_gevi_rest = np.asarray(psd_gevi_rest).flatten()
                psd_gevi_run = np.asarray(psd_gevi_run).flatten()
                if len(psd_gevi_rest) > 0 and len(psd_gevi_run) > 0:
                    gevi_theta_r = compute_theta_band_value(freq, psd_gevi_rest)
                    gevi_theta_n = compute_theta_band_value(freq, psd_gevi_run)
                    # Compute peakiness metrics (NEW: 1/f deviation and local peak index)
                    gevi_prom_r = compute_1f_deviation(freq, psd_gevi_rest)
                    gevi_prom_n = compute_1f_deviation(freq, psd_gevi_run)
                    gevi_ta_r = compute_local_peak_index(freq, psd_gevi_rest)
                    gevi_ta_n = compute_local_peak_index(freq, psd_gevi_run)
            
            # Store results
            results['animal_ids'].append(mouse_id)
            results['coh_theta_rest'].append(coh_theta_r)
            results['coh_theta_run'].append(coh_theta_n)
            results['lfp_theta_rest'].append(lfp_theta_r)
            results['lfp_theta_run'].append(lfp_theta_n)
            results['gevi_theta_rest'].append(gevi_theta_r)
            results['gevi_theta_run'].append(gevi_theta_n)
            results['gevi_prominence_rest'].append(gevi_prom_r)
            results['gevi_prominence_run'].append(gevi_prom_n)
            results['gevi_theta_alpha_rest'].append(gevi_ta_r)
            results['gevi_theta_alpha_run'].append(gevi_ta_n)
            
            print(f"      [{mouse_id}] OK: Coh={coh_theta_r:.3f}→{coh_theta_n:.3f}, "
                  f"LFP={lfp_theta_r:.1f}→{lfp_theta_n:.1f} dB, "
                  f"GEVI Prom={gevi_prom_r:.1f}→{gevi_prom_n:.1f}")
            
        except Exception as e:
            print(f"      [{mouse_id}] ERROR: {e}")
    
    # Convert to numpy arrays
    for key in results:
        if key != 'animal_ids':
            results[key] = np.array(results[key])
    
    print(f"    Loaded {len(results['animal_ids'])} animals")
    
    return results


# ==============================================================================
#  STATISTICAL TESTS
# ==============================================================================

def perform_one_tailed_test(rest_data, run_data, alternative='greater'):
    """
    Perform one-tailed paired test (hypothesis: RUN > REST).
    
    Uses paired t-test if data passes normality, otherwise Wilcoxon.
    'greater' means we test if RUN > REST (i.e., run_data - rest_data > 0)
    
    Returns dict with p_value, statistic, effect_size, n, test_used
    """
    rest_data = np.asarray(rest_data)
    run_data = np.asarray(run_data)
    
    # Remove NaN pairs
    valid = ~(np.isnan(rest_data) | np.isnan(run_data))
    rest_data = rest_data[valid]
    run_data = run_data[valid]
    
    n = len(rest_data)
    if n < 3:
        return {'p_value': np.nan, 'statistic': np.nan, 'effect_size': np.nan, 
                'n': n, 'test_used': 'insufficient_data'}
    
    # Calculate differences (RUN - REST)
    differences = run_data - rest_data
    
    # Check normality of differences
    try:
        _, p_normality = shapiro(differences)
        is_normal = p_normality > 0.05
    except:
        is_normal = False
    
    # Calculate effect size (Cohen's d for paired data)
    mean_diff = np.mean(differences)
    std_diff = np.std(differences, ddof=1)
    effect_size = mean_diff / std_diff if std_diff > 0 else 0
    
    try:
        if is_normal or n >= 7:  # Use t-test for normal data or larger samples
            # One-tailed paired t-test
            stat, p_two_tailed = ttest_rel(run_data, rest_data)
            # Convert to one-tailed
            if alternative == 'greater':
                p_value = p_two_tailed / 2 if stat > 0 else 1 - p_two_tailed / 2
            else:
                p_value = p_two_tailed / 2 if stat < 0 else 1 - p_two_tailed / 2
            test_used = 'paired_ttest_1tail'
        else:
            # Wilcoxon signed-rank test (one-tailed)
            stat, p_value = wilcoxon(run_data, rest_data, alternative=alternative)
            test_used = 'wilcoxon_1tail'
        
        return {
            'p_value': p_value,
            'statistic': stat,
            'effect_size': effect_size,
            'n': n,
            'test_used': test_used,
            'mean_diff': mean_diff,
        }
    except Exception as e:
        return {'p_value': np.nan, 'statistic': np.nan, 'effect_size': np.nan,
                'n': n, 'test_used': 'error'}


# ==============================================================================
#  PLOTTING FUNCTIONS
# ==============================================================================

def plot_half_violin(ax, data_left, data_right, positions, colors, labels, 
                     animal_ids=None, show_labels=False):
    """Plot half-violin with box and scatter (same style as FOOOF)."""
    from scipy.stats import gaussian_kde
    
    pos_left, pos_right = positions
    color_left, color_right = colors
    
    violin_base_left = pos_left - 0.05
    box_pos_left = pos_left + BOX_OFFSET
    dot_pos_left = pos_left + DOT_OFFSET
    
    violin_base_right = pos_right + 0.05
    box_pos_right = pos_right - BOX_OFFSET
    dot_pos_right = pos_right - DOT_OFFSET
    
    def get_violin_path(data, position, side='left', width=VIOLIN_WIDTH):
        if len(data) < 2 or np.std(data) < 1e-10:
            return None, None, None
        try:
            kde = gaussian_kde(data, bw_method=0.5)
            data_range = np.max(data) - np.min(data)
            padding = max(data_range * 0.2, np.std(data) * 0.3)
            y_range = np.linspace(np.min(data) - padding, np.max(data) + padding, 100)
            density = kde(y_range)
            density = density / np.max(density) * width
            x_path = position - density if side == 'left' else position + density
            return y_range, x_path, density
        except:
            return None, None, None
    
    # Half-violins
    result_left = get_violin_path(data_left, violin_base_left, 'left')
    if result_left[0] is not None:
        y_range, x_path, _ = result_left
        ax.fill_betweenx(y_range, violin_base_left, x_path, alpha=0.6, color=color_left)
    
    result_right = get_violin_path(data_right, violin_base_right, 'right')
    if result_right[0] is not None:
        y_range, x_path, _ = result_right
        ax.fill_betweenx(y_range, violin_base_right, x_path, alpha=0.6, color=color_right)
    
    # Box plots
    bp_left = ax.boxplot([data_left], positions=[box_pos_left], widths=BOX_WIDTH,
                         patch_artist=True, showfliers=False, zorder=3)
    bp_right = ax.boxplot([data_right], positions=[box_pos_right], widths=BOX_WIDTH,
                          patch_artist=True, showfliers=False, zorder=3)
    
    for bp, color in [(bp_left, color_left), (bp_right, color_right)]:
        for patch in bp['boxes']:
            patch.set_facecolor('white')
            patch.set_edgecolor(color)
            patch.set_linewidth(2)
        for whisker in bp['whiskers']:
            whisker.set_color(color)
            whisker.set_linewidth(2)
        for cap in bp['caps']:
            cap.set_color(color)
            cap.set_linewidth(2)
        for median in bp['medians']:
            median.set_color('black')
            median.set_linewidth(2.5)
    
    # Scatter points
    ax.scatter([dot_pos_left] * len(data_left), data_left, s=DOT_SIZE, c=[color_left],
               edgecolor='white', linewidth=1.5, zorder=5, alpha=0.9)
    ax.scatter([dot_pos_right] * len(data_right), data_right, s=DOT_SIZE, c=[color_right],
               edgecolor='white', linewidth=1.5, zorder=5, alpha=0.9)
    
    # Animal labels
    if show_labels and animal_ids is not None:
        for i, animal_id in enumerate(animal_ids):
            ax.annotate(animal_id, (dot_pos_left - 0.12, data_left[i]),
                       fontsize=5, ha='right', va='center', color='gray', alpha=0.9)
    
    # Paired lines
    for i in range(len(data_left)):
        ax.plot([dot_pos_left, dot_pos_right], [data_left[i], data_right[i]],
                color='gray', alpha=LINE_ALPHA, linewidth=LINE_WIDTH, zorder=4)


def add_significance_annotation(ax, p_value, y_min, y_max, positions):
    """Add significance stars above comparison."""
    pos_left, pos_right = positions
    x_center = (pos_left + pos_right) / 2
    
    y_range = y_max - y_min
    y_bracket = y_max + y_range * 0.08
    y_text = y_bracket + y_range * 0.05
    
    # Significance stars
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'n.s.'
    
    # Draw bracket
    bracket_width = 0.15
    ax.plot([pos_left + bracket_width, pos_right - bracket_width], 
            [y_bracket, y_bracket], color='black', linewidth=1.5)
    ax.plot([pos_left + bracket_width, pos_left + bracket_width],
            [y_bracket - y_range * 0.02, y_bracket], color='black', linewidth=1.5)
    ax.plot([pos_right - bracket_width, pos_right - bracket_width],
            [y_bracket - y_range * 0.02, y_bracket], color='black', linewidth=1.5)
    
    ax.text(x_center, y_text, sig_text, ha='center', va='bottom',
            fontsize=FONT_SIZE_TITLE, fontweight='bold',
            fontstyle='italic' if sig_text == 'n.s.' else 'normal')


# ==============================================================================
#  MAIN FIGURE CREATION
# ==============================================================================

def create_theta_band_figure(theta_data, output_path, method_title='mscohere'):
    """
    Create publication-ready theta band figure.
    
    Layout: 1 row × 5 columns
        Col 1: LFP-GEVI Coherence (theta band)
        Col 2: LFP Power (theta band)
        Col 3: GEVI Power (theta band)
        Col 4: GEVI Peak Prominence (peakiness metric 1)
        Col 5: GEVI Theta/Alpha Ratio (peakiness metric 2)
    """
    
    fig, axes = plt.subplots(1, 4, figsize=(20, 6))
    fig.subplots_adjust(left=0.06, right=0.96, top=0.82, bottom=0.12, wspace=0.35)
    
    positions = [0.7, 1.3]
    labels = ['REST', 'RUN']
    
    animal_ids = theta_data['animal_ids']
    
    # Colors for each panel
    # Coherence: Uses indigo/blue-violet from central config (blend of LFP purple and GEVI teal)
    colors_coh = [COLOR_COH_REST, COLOR_COH_RUN]
    
    colors_lfp = [COLOR_LFP_REST, COLOR_LFP_RUN]  # LFP: Purple
    colors_gevi = [COLOR_GEVI_REST, COLOR_GEVI_RUN]  # GEVI: Teal
    # Oscillatory Theta: Same colors as GEVI (panels 3 and 4 match)
    colors_osc = [COLOR_GEVI_REST, COLOR_GEVI_RUN]
    
    # Panel configurations (4 panels only)
    panels = [
        {
            'data_rest': theta_data['coh_theta_rest'],
            'data_run': theta_data['coh_theta_run'],
            'title': f'LFP-GEVI Coherence\n(Theta: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz)',
            'ylabel': 'Coherence',
            'colors': colors_coh,
        },
        {
            'data_rest': theta_data['lfp_theta_rest'],
            'data_run': theta_data['lfp_theta_run'],
            'title': f'LFP Power\n(Theta: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz)',
            'ylabel': 'Power (dB/Hz)',
            'colors': colors_lfp,
        },
        {
            'data_rest': theta_data['gevi_theta_rest'],
            'data_run': theta_data['gevi_theta_run'],
            'title': f'GEVI Power\n(Theta: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz)',
            'ylabel': 'Power (dB/Hz)',
            'colors': colors_gevi,
        },
        {
            'data_rest': theta_data['gevi_prominence_rest'],
            'data_run': theta_data['gevi_prominence_run'],
            'title': f'GEVI Oscillatory Theta\n(Aperiodic-Corrected)',
            'ylabel': 'Power above 1/f (dB)',
            'colors': colors_osc,
        },
    ]
    
    for col_idx, panel in enumerate(panels):
        ax = axes[col_idx]
        
        data_rest = np.asarray(panel['data_rest'])
        data_run = np.asarray(panel['data_run'])
        
        # Remove NaN values (keep pairs)
        valid = ~(np.isnan(data_rest) | np.isnan(data_run))
        data_rest_valid = data_rest[valid]
        data_run_valid = data_run[valid]
        animal_ids_valid = [animal_ids[i] for i in range(len(animal_ids)) if valid[i]]
        
        if len(data_rest_valid) < 2:
            ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                    ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
            ax.set_title(panel['title'], fontsize=FONT_SIZE_TITLE, fontweight='bold')
            continue
        
        # Plot
        plot_half_violin(ax, data_rest_valid, data_run_valid, positions, 
                        panel['colors'], labels, animal_ids_valid, SHOW_ANIMAL_LABELS)
        
        # Statistics (one-tailed: RUN > REST)
        stats = perform_one_tailed_test(data_rest_valid, data_run_valid, alternative='greater')
        
        # Y limits
        all_data = np.concatenate([data_rest_valid, data_run_valid])
        y_min, y_max = np.min(all_data), np.max(all_data)
        y_range = y_max - y_min
        if y_range < 1e-10:
            y_range = abs(y_max) * 0.2 if y_max != 0 else 1.0
        ax.set_ylim([y_min - y_range * 0.40, y_max + y_range * 0.35])
        
        # Significance annotation
        add_significance_annotation(ax, stats['p_value'], y_min, y_max, positions)
        
        # Styling
        ax.set_ylabel(panel['ylabel'], fontsize=FONT_SIZE_LABEL)
        ax.set_title(panel['title'], fontsize=FONT_SIZE_TITLE, fontweight='bold')
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
        ax.set_xlim([0.15, 1.85])
        
        # Stats text box
        test_abbrev = {'paired_ttest_1tail': 't-test (1-tail)', 
                       'wilcoxon_1tail': 'Wilcoxon (1-tail)'}.get(stats['test_used'], stats['test_used'])
        stats_text = f"p={stats['p_value']:.3f}\n{test_abbrev}\nN={stats['n']}\nd={stats['effect_size']:.2f}"
        ax.text(0.02, 0.02, stats_text, transform=ax.transAxes, fontsize=9,
                verticalalignment='bottom', fontfamily='monospace',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
        
        # Spines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
        ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
        ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICK,
                      width=TICK_WIDTH, length=TICK_LENGTH)
    
    # Legend (use GEVI colors as representative since 3/4 panels are GEVI-related)
    legend_elements = [
        Patch(facecolor=COLOR_GEVI_REST, alpha=0.6, label='REST', edgecolor='black', linewidth=1),
        Patch(facecolor=COLOR_GEVI_RUN, alpha=0.6, label='RUN', edgecolor='black', linewidth=1),
    ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=FONT_SIZE_LEGEND,
               frameon=True, bbox_to_anchor=(0.98, 0.98), framealpha=0.9, edgecolor='gray')
    
    # Suptitle
    pooling_label = get_pooling_level_label(GROUP_POOLING_LEVEL)
    method_label = {'peak': 'Peak', 'mean': 'Mean', 'auc': 'AUC'}.get(THETA_EXTRACTION_METHOD, THETA_EXTRACTION_METHOD)
    fig.suptitle(f'Theta Band Analysis ({method_label}): REST vs RUN\n{method_title} ({pooling_label}) - One-tailed test: RUN > REST',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)
    
    # Save
    output_dir = Path(output_path).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {Path(output_path).name}")


def run_theta_band_pipeline(animals_to_process, output_dir, method='mscohere', pooling_level=None):
    """Run theta band analysis pipeline."""
    
    if pooling_level is None:
        pooling_level = GROUP_POOLING_LEVEL
    
    pooling_label = get_pooling_level_label(pooling_level)
    
    print(f"\n    Loading theta band data ({pooling_label})...")
    theta_data = load_theta_band_data(animals_to_process, method, pooling_level)
    
    if len(theta_data['animal_ids']) == 0:
        print("    No data loaded, skipping theta band analysis")
        return 0
    
    # Print diagnostic table (4 metrics)
    print(f"\n    Theta Band Values (REST → RUN):")
    print(f"    " + "-" * 105)
    print(f"    {'Animal':<15s} | {'Coherence':<18s} | {'LFP Power':<18s} | {'GEVI Power':<18s} | {'Osc.Theta':<18s}")
    print(f"    " + "-" * 105)
    for i, animal_id in enumerate(theta_data['animal_ids']):
        coh_r, coh_n = theta_data['coh_theta_rest'][i], theta_data['coh_theta_run'][i]
        lfp_r, lfp_n = theta_data['lfp_theta_rest'][i], theta_data['lfp_theta_run'][i]
        gevi_r, gevi_n = theta_data['gevi_theta_rest'][i], theta_data['gevi_theta_run'][i]
        prom_r, prom_n = theta_data['gevi_prominence_rest'][i], theta_data['gevi_prominence_run'][i]
        
        coh_diff = coh_n - coh_r
        lfp_diff = lfp_n - lfp_r
        gevi_diff = gevi_n - gevi_r
        prom_diff = prom_n - prom_r
        
        coh_dir = "↑" if coh_diff > 0 else "↓"
        lfp_dir = "↑" if lfp_diff > 0 else "↓"
        gevi_dir = "↑" if gevi_diff > 0 else "↓"
        prom_dir = "↑" if prom_diff > 0 else "↓"
        
        print(f"    {animal_id:<15s} | {coh_r:.3f}→{coh_n:.3f} {coh_dir} | "
              f"{lfp_r:.1f}→{lfp_n:.1f} {lfp_dir} | "
              f"{gevi_r:.1f}→{gevi_n:.1f} {gevi_dir} | "
              f"{prom_r:.1f}→{prom_n:.1f} {prom_dir}")
    print(f"    " + "-" * 105)
    
    # Create figure
    level_abbrev = 'animal_concat' if pooling_level == 'animal_concatenated' else 'animal_pooled'
    output_path = output_dir / f"theta_band_{level_abbrev}"
    create_theta_band_figure(theta_data, output_path, method_title=method.upper())
    
    return 1


# ==============================================================================
#  MAIN
# ==============================================================================

def main():
    """Standalone execution."""
    print("=" * 70)
    print("  THETA BAND ANALYSIS")
    print("=" * 70)
    
    if not USING_CENTRAL_CONFIG:
        print("\nERROR: Central config not available")
        return
    
    animals = get_animals_to_process()
    output_dir = get_group_level_output_dir()
    pooling_level = GROUP_POOLING_LEVEL
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pooling_label = get_pooling_level_label(pooling_level)
    
    print(f"  Output: {output_dir}")
    print(f"  Animals: {[a['mouse_id'] for a in animals]}")
    print(f"  Pooling Level: {pooling_label}")
    print(f"  Spectral Method: {SPECTRAL_METHOD}")
    print(f"  Theta Range: {THETA_RANGE[0]}-{THETA_RANGE[1]} Hz")
    print(f"  1/f fit bands: 2-4 Hz, 15-25 Hz (flanks for 1/f deviation)")
    print(f"  Extraction: {THETA_EXTRACTION_METHOD} (peak=max, mean=average)")
    print(f"  Hypothesis: RUN > REST (one-tailed test)")
    print(f"  Oscillatory Theta: theta peak above fitted 1/f (aperiodic-corrected)")
    
    # Run analysis using configured spectral method
    run_theta_band_pipeline(animals, output_dir, method=SPECTRAL_METHOD, pooling_level=pooling_level)
    
    print("\n" + "=" * 70)
    print("  THETA BAND ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
