"""
================================================================================
FIGURE 4: FOOOF (Fitting Oscillations & One Over F) Analysis
================================================================================

Publication-ready figure showing oscillatory analysis with 1/f decomposition.
Compares REST vs RUN using FOOOF-derived measures:
  - Theta Peak Amplitude (oscillatory power above 1/f)
  - Aperiodic Exponent (1/f slope)
  - Aperiodic Offset (broadband power level)

Visualization: Half-violin plots with individual animal dots and paired lines.
Statistics: Wilcoxon signed-rank test (non-parametric paired test for small N).

NOTE: FOOOF analyzes PSD (Power Spectral Density), NOT coherence.
      PSD is computed identically for both mscohere and FieldTrip methods
      (both use MATLAB's pwelch function), so FOOOF only needs to run once.

================================================================================
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy.io import loadmat
from scipy.stats import wilcoxon, ttest_rel, shapiro
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
        get_animal_pooled_input_dir, get_animal_concatenated_input_dir,
        get_group_level_output_dir, GROUP_POOLING_LEVEL,
        FONT_SIZE_SUPTITLE, FONT_SIZE_TITLE, FONT_SIZE_LABEL, FONT_SIZE_TICK,
        FONT_SIZE_LEGEND, AXIS_LINEWIDTH, TICK_WIDTH, TICK_LENGTH,
        COLOR_LFP_REST, COLOR_LFP_RUN, COLOR_GEVI_REST, COLOR_GEVI_RUN,
        ANIMALS, get_animals_to_process,
    )
    USING_CENTRAL_CONFIG = True
except ImportError:
    USING_CENTRAL_CONFIG = False
    # Fallback defaults
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
    # LFP colors: Purple shades
    COLOR_LFP_REST = np.array([0.25, 0.18, 0.35])
    COLOR_LFP_RUN = np.array([0.55, 0.45, 0.65])
    # GEVI colors: Teal shades
    COLOR_GEVI_REST = np.array([0.05, 0.35, 0.45])
    COLOR_GEVI_RUN = np.array([0.25, 0.65, 0.65])
    GROUP_POOLING_LEVEL = 'animal_concatenated'  # Fallback default

# FOOOF parameters
FOOOF_FREQ_RANGE = [2, 70]  # Frequency range to fit
THETA_RANGE = [5, 9]        # Theta band for peak detection
APERIODIC_MODE = 'fixed'    # 'fixed' = offset + exponent, 'knee' = adds knee parameter

# Statistical test configuration
# For theta peak amplitude: directional hypothesis (RUN > REST) → use one-tailed test
# For exponent/offset: may be non-directional → use two-tailed test
USE_ONE_TAILED_THETA = True   # One-tailed test for theta (RUN > REST hypothesis)
USE_ONE_TAILED_EXPONENT = False  # Two-tailed for exponent (no directional hypothesis)
USE_ONE_TAILED_OFFSET = False   # Two-tailed for offset (no directional hypothesis)

# Test selection: 'auto' (choose based on normality), 'ttest' (paired t-test), 'wilcoxon' (Wilcoxon)
# 'auto' uses t-test if data passes normality (Shapiro-Wilk), otherwise Wilcoxon
STATISTICAL_TEST = 'auto'  # 'auto', 'ttest', or 'wilcoxon'

# Half-violin plot parameters (IMPROVED - separated layout)
VIOLIN_WIDTH = 0.28         # Narrower violin for cleaner look
BOX_WIDTH = 0.10            # Narrow box
BOX_OFFSET = 0.08           # Box shifted inward from violin position
DOT_OFFSET = 0.18           # Dots shifted further inward (between boxes)
DOT_SIZE = 90               # Dot size
LINE_WIDTH = 1.5
LINE_ALPHA = 0.5

# Animal label display
SHOW_ANIMAL_LABELS = True   # Set to True to show animal IDs next to data points

# Check if FOOOF is installed
try:
    from fooof import FOOOF
    FOOOF_AVAILABLE = True
except ImportError:
    FOOOF_AVAILABLE = False
    print("WARNING: FOOOF not installed. Install with: pip install fooof")


# ==============================================================================
#  WINDOWS LONG PATH SUPPORT
# ==============================================================================

from common import to_long_path  # shared helper (was a local copy)


def force_create_directory(dir_path):
    """
    Create directory, handling edge case where a file exists with the same name.
    This shouldn't normally happen, but can occur if something went wrong during a previous run.
    """
    import os
    import stat
    dir_path = Path(dir_path)
    dir_path_str = str(dir_path)
    dir_path_long = to_long_path(dir_path_str)
    
    # Check if it's already a directory (try both path formats)
    for path_to_check in [dir_path, Path(dir_path_long)]:
        try:
            if path_to_check.exists() and path_to_check.is_dir():
                return True
        except Exception:
            pass
    
    # Try to detect if a FILE exists (check both path formats)
    # Sometimes Path.exists() doesn't work properly, so we try to access it directly
    for path_to_check, path_str in [(dir_path, dir_path_str), (Path(dir_path_long), dir_path_long)]:
        try:
            # Try to stat it - this will tell us if it exists and what type it is
            stat_info = os.stat(path_str)
            if stat.S_ISDIR(stat_info.st_mode):
                # It's a directory - we're done
                return True
            elif stat.S_ISREG(stat_info.st_mode):
                # It's a regular file - remove it
                try:
                    os.remove(path_str)
                    print(f"      Removed file blocking directory: {dir_path_str}")
                except Exception as e:
                    print(f"      ERROR: Could not remove blocking file: {e}")
                    raise OSError(f"File exists at directory path and cannot be removed: {dir_path_str}")
        except FileNotFoundError:
            # Path doesn't exist - that's fine, we'll create it
            pass
        except OSError as e:
            # If it's error 183, it means a file exists but we can't stat it properly
            # Try to remove it anyway
            if hasattr(e, 'winerror') and e.winerror == 183:
                try:
                    os.remove(path_str)
                    print(f"      Removed file blocking directory (via error detection): {dir_path_str}")
                except Exception:
                    pass
        except Exception:
            # Other error - continue to next path format
            pass
    
    # Create the directory (try both path formats)
    for path_to_create, path_str in [(dir_path, dir_path_str), (Path(dir_path_long), dir_path_long)]:
        try:
            path_to_create.mkdir(parents=True, exist_ok=True)
            # Verify it was created
            try:
                stat_info = os.stat(path_str)
                if stat.S_ISDIR(stat_info.st_mode):
                    return True
            except Exception:
                pass
        except (FileExistsError, OSError) as e:
            # This is the error we're trying to fix - a file exists (error 183)
            # Try to remove it and recreate, even if we can't detect it
            winerror = getattr(e, 'winerror', None)
            if winerror == 183 or isinstance(e, FileExistsError):
                # Force remove the file (it might be hidden or not detectable)
                try:
                    # Try removing it - if it doesn't exist, that's fine
                    try:
                        os.remove(path_str)
                        print(f"      Removed blocking file (during creation): {dir_path_str}")
                    except FileNotFoundError:
                        # File doesn't exist - maybe it was already removed
                        pass
                    except Exception as remove_err:
                        print(f"      WARNING: Could not remove file: {remove_err}")
                    
                    # Now try creating the directory again
                    path_to_create.mkdir(parents=True, exist_ok=True)
                    # Verify it was created
                    try:
                        stat_info = os.stat(path_str)
                        if stat.S_ISDIR(stat_info.st_mode):
                            return True
                    except Exception:
                        pass
                except Exception as e2:
                    # Continue to next path format
                    continue
            else:
                # Other error - continue to next path format
                continue
        except Exception:
            # Try next path format
            continue
    
    # If all else fails, raise an error
    raise OSError(f"Could not create directory after all attempts: {dir_path_str}")


# ==============================================================================
#  DATA LOADING UTILITIES
# ==============================================================================

def load_animal_psd_data(animals_to_process, method='mscohere', pooling_level=None):
    """
    Load PSD data for each animal (REST and RUN separately).
    
    Parameters:
        animals_to_process: list of animal dictionaries
        method: 'mscohere' or 'fieldtrip'
        pooling_level: 'animal_pooled' or 'animal_concatenated' (defaults to GROUP_POOLING_LEVEL)
    
    Returns dict with animal_id keys, each containing:
        - freq: frequency axis
        - psd_lfp_rest, psd_lfp_run: LFP PSD for REST/RUN
        - psd_gevi_rest, psd_gevi_run: GEVI PSD for REST/RUN
    """
    animal_data = {}
    
    # Use the specified pooling level, defaulting to GROUP_POOLING_LEVEL
    if pooling_level is None:
        pooling_level = GROUP_POOLING_LEVEL
    
    # Use central config's path functions which handle suffix properly
    # Structure: {BASE_OUTPUT_DIR}{SUFFIX}/{BEHAVIOR_MODE}/{pooling_level}/{mouse_id}/data/
    from plotting_config import get_base_dir
    base_path = get_base_dir() / pooling_level
    
    print(f"    Attempting to load {len(animals_to_process)} animals:")
    print(f"    Expected: {[a['mouse_id'] for a in animals_to_process]}")
    print(f"    Pooling level: {pooling_level}")
    print(f"    Base path: {base_path}")
    
    for animal in animals_to_process:
        mouse_id = animal['mouse_id']
        data_path = base_path / mouse_id / 'data' / f'{method}.mat'
        data_path_long = to_long_path(str(data_path))
        
        if not os.path.exists(data_path_long):
            print(f"      [{mouse_id}] SKIP: File not found")
            print(f"                    Path: {data_path}")
            continue
        
        try:
            mat = loadmat(data_path_long, squeeze_me=True)
            
            # Debug: print available keys
            available_keys = [k for k in mat.keys() if not k.startswith('_')]
            
            # Extract frequency axis
            freq = mat.get('psd_freq', mat.get('freq', None))
            if freq is None:
                print(f"      [{mouse_id}] SKIP: No frequency axis (available: {available_keys})")
                continue
            freq = np.asarray(freq).flatten()
            
            # Extract PSD data (already in dB)
            psd_lfp_rest = mat.get('psd_lfp_rest', None)
            psd_lfp_run = mat.get('psd_lfp_run', None)
            psd_gevi_rest = mat.get('psd_gevi_rest', None)
            psd_gevi_run = mat.get('psd_gevi_run', None)
            
            # Check for empty arrays - LFP REST
            if psd_lfp_rest is not None:
                psd_lfp_rest = np.asarray(psd_lfp_rest).flatten()
                if len(psd_lfp_rest) == 0:
                    print(f"      [{mouse_id}] SKIP: psd_lfp_rest is empty array")
                    continue
                if np.all(np.isnan(psd_lfp_rest)):
                    print(f"      [{mouse_id}] SKIP: psd_lfp_rest is all NaN")
                    continue
            else:
                print(f"      [{mouse_id}] SKIP: psd_lfp_rest is None (available: {available_keys})")
                continue
            
            # Check for empty arrays - LFP RUN
            if psd_lfp_run is not None:
                psd_lfp_run = np.asarray(psd_lfp_run).flatten()
                if len(psd_lfp_run) == 0:
                    print(f"      [{mouse_id}] SKIP: psd_lfp_run is empty array")
                    continue
                if np.all(np.isnan(psd_lfp_run)):
                    print(f"      [{mouse_id}] SKIP: psd_lfp_run is all NaN")
                    continue
            else:
                print(f"      [{mouse_id}] SKIP: psd_lfp_run is None (available: {available_keys})")
                continue
            
            # Process GEVI (optional - don't skip if missing)
            if psd_gevi_rest is not None:
                psd_gevi_rest = np.asarray(psd_gevi_rest).flatten()
                if len(psd_gevi_rest) == 0 or np.all(np.isnan(psd_gevi_rest)):
                    psd_gevi_rest = None
            
            if psd_gevi_run is not None:
                psd_gevi_run = np.asarray(psd_gevi_run).flatten()
                if len(psd_gevi_run) == 0 or np.all(np.isnan(psd_gevi_run)):
                    psd_gevi_run = None
            
            animal_data[mouse_id] = {
                'freq': freq,
                'psd_lfp_rest': psd_lfp_rest,
                'psd_lfp_run': psd_lfp_run,
                'psd_gevi_rest': psd_gevi_rest,
                'psd_gevi_run': psd_gevi_run,
            }
            gevi_status = "GEVI OK" if (psd_gevi_rest is not None and psd_gevi_run is not None) else "GEVI missing"
            print(f"      [{mouse_id}] OK: freq={len(freq)}, LFP OK, {gevi_status}")
            
        except Exception as e:
            import traceback
            print(f"      [{mouse_id}] ERROR: {e}")
            traceback.print_exc()
    
    print(f"\n    Successfully loaded: {len(animal_data)} / {len(animals_to_process)} animals")
    if len(animal_data) < len(animals_to_process):
        missing = set([a['mouse_id'] for a in animals_to_process]) - set(animal_data.keys())
        print(f"    MISSING: {missing}")
    
    return animal_data


# ==============================================================================
#  FOOOF ANALYSIS
# ==============================================================================

def run_fooof_on_psd(freq, psd_db, freq_range=FOOOF_FREQ_RANGE):
    """
    Run FOOOF on a single PSD spectrum.
    
    Parameters:
        freq: Frequency axis (Hz)
        psd_db: Power spectrum in dB (will be converted to linear)
        freq_range: [min_freq, max_freq] for fitting
    
    Returns:
        dict with: aperiodic_offset, aperiodic_exponent, theta_peak_amplitude, theta_peak_cf
        Returns None if FOOOF fitting fails
    """
    if not FOOOF_AVAILABLE:
        return None
    
    # Convert dB to linear power (FOOOF expects linear)
    psd_linear = 10 ** (psd_db / 10)
    
    # Initialize FOOOF
    fm = FOOOF(
        peak_width_limits=[1, 8],      # Peak width range (Hz)
        max_n_peaks=6,                  # Max number of peaks to fit
        min_peak_height=0.1,            # Minimum peak height (in log power)
        peak_threshold=2.0,             # Threshold for peak detection
        aperiodic_mode=APERIODIC_MODE,  # 'fixed' or 'knee'
    )
    
    try:
        # Fit the model
        fm.fit(freq, psd_linear, freq_range)
        
        # Extract aperiodic parameters
        # For 'fixed' mode: [offset, exponent]
        # For 'knee' mode: [offset, knee, exponent]
        ap_params = fm.aperiodic_params_
        aperiodic_offset = ap_params[0]
        aperiodic_exponent = ap_params[-1]  # Last param is always exponent
        
        # Extract theta peak (if present)
        theta_peak_amplitude = 0.0
        theta_peak_cf = np.nan
        
        # fm.peak_params_ is array of [CF, PW, BW] for each peak
        # CF = center frequency, PW = power (height), BW = bandwidth
        if len(fm.peak_params_) > 0:
            for peak in fm.peak_params_:
                cf, pw, bw = peak
                if THETA_RANGE[0] <= cf <= THETA_RANGE[1]:
                    if pw > theta_peak_amplitude:
                        theta_peak_amplitude = pw
                        theta_peak_cf = cf
        
        return {
            'aperiodic_offset': aperiodic_offset,
            'aperiodic_exponent': aperiodic_exponent,
            'theta_peak_amplitude': theta_peak_amplitude,
            'theta_peak_cf': theta_peak_cf,
            'model_r2': fm.r_squared_,
            'model_error': fm.error_,
        }
        
    except Exception as e:
        print(f"      FOOOF fitting failed: {e}")
        return None


def run_fooof_analysis(animal_data, signal_type='lfp'):
    """
    Run FOOOF on all animals for REST and RUN conditions.
    
    Returns:
        dict with arrays for each measure:
            - theta_rest, theta_run
            - exponent_rest, exponent_run
            - offset_rest, offset_run
            - animal_ids (for labeling)
    """
    results = {
        'animal_ids': [],
        'theta_rest': [],
        'theta_run': [],
        'exponent_rest': [],
        'exponent_run': [],
        'offset_rest': [],
        'offset_run': [],
        'r2_rest': [],
        'r2_run': [],
    }
    
    psd_key_rest = f'psd_{signal_type}_rest'
    psd_key_run = f'psd_{signal_type}_run'
    
    print(f"      Running FOOOF on {len(animal_data)} animals for {signal_type.upper()}...")
    
    for animal_id, data in animal_data.items():
        freq = data['freq']
        psd_rest = data.get(psd_key_rest)
        psd_run = data.get(psd_key_run)
        
        if psd_rest is None or psd_run is None:
            print(f"        [{animal_id}] SKIP: Missing {signal_type.upper()} PSD (REST={psd_rest is not None}, RUN={psd_run is not None})")
            continue
        
        # Check for valid data
        if len(psd_rest) == 0 or len(psd_run) == 0:
            print(f"        [{animal_id}] SKIP: Empty PSD arrays")
            continue
        if np.all(np.isnan(psd_rest)) or np.all(np.isnan(psd_run)):
            print(f"        [{animal_id}] SKIP: All NaN values")
            continue
        
        print(f"        [{animal_id}] Fitting FOOOF...")
        
        # Run FOOOF on REST
        fooof_rest = run_fooof_on_psd(freq, psd_rest)
        if fooof_rest is None:
            print(f"        [{animal_id}] SKIP: FOOOF fitting failed for REST")
            continue
        
        # Run FOOOF on RUN
        fooof_run = run_fooof_on_psd(freq, psd_run)
        if fooof_run is None:
            print(f"        [{animal_id}] SKIP: FOOOF fitting failed for RUN")
            continue
        
        print(f"        [{animal_id}] OK: R²(REST)={fooof_rest['model_r2']:.3f}, R²(RUN)={fooof_run['model_r2']:.3f}")
        
        # Store results
        results['animal_ids'].append(animal_id)
        results['theta_rest'].append(fooof_rest['theta_peak_amplitude'])
        results['theta_run'].append(fooof_run['theta_peak_amplitude'])
        results['exponent_rest'].append(fooof_rest['aperiodic_exponent'])
        results['exponent_run'].append(fooof_run['aperiodic_exponent'])
        results['offset_rest'].append(fooof_rest['aperiodic_offset'])
        results['offset_run'].append(fooof_run['aperiodic_offset'])
        results['r2_rest'].append(fooof_rest['model_r2'])
        results['r2_run'].append(fooof_run['model_r2'])
    
    # Convert to numpy arrays
    for key in results:
        if key != 'animal_ids':
            results[key] = np.array(results[key])
    
    return results


# ==============================================================================
#  STATISTICAL TESTS
# ==============================================================================

def perform_statistical_test(rest_data, run_data, alternative='two-sided', test_type='auto'):
    """
    Perform paired statistical test with automatic test selection.
    
    Parameters:
        rest_data: array of REST condition values
        run_data: array of RUN condition values (paired with rest_data)
        alternative: 'two-sided', 'greater' (RUN > REST), or 'less' (REST > RUN)
        test_type: 'auto' (choose based on normality), 'ttest' (paired t-test), 'wilcoxon' (Wilcoxon)
    
    Returns:
        dict with: statistic, p_value, effect_size, n, test_used
    """
    n = len(rest_data)
    
    if n < 3:
        return {'statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan, 'n': n, 'test_used': 'none'}
    
    # Calculate differences for normality check
    diff = run_data - rest_data
    
    # Auto-select test based on normality
    if test_type == 'auto':
        try:
            # Shapiro-Wilk test for normality (only if n >= 3 and n <= 5000)
            if 3 <= n <= 5000:
                _, p_norm = shapiro(diff)
                use_ttest = p_norm > 0.05  # Use t-test if normal (p > 0.05)
            else:
                use_ttest = True  # For very small or very large N, default to t-test
        except:
            use_ttest = False  # If normality test fails, use Wilcoxon
    elif test_type == 'ttest':
        use_ttest = True
    else:  # 'wilcoxon'
        use_ttest = False
    
    try:
        if use_ttest:
            # Paired t-test (more powerful if data is normal)
            stat, p = ttest_rel(run_data, rest_data, alternative=alternative)
            
            # Effect size: Cohen's d for paired samples
            # d = mean_diff / std_diff
            mean_diff = np.mean(diff)
            std_diff = np.std(diff, ddof=1)
            if std_diff > 1e-10:
                d = mean_diff / std_diff
            else:
                d = 0.0
            
            test_used = 'paired_ttest'
            effect_size = d
            
        else:
            # Wilcoxon signed-rank test (non-parametric, robust)
            stat, p = wilcoxon(rest_data, run_data, alternative=alternative)
            
            # Effect size: matched-pairs rank-biserial correlation
            # r = 1 - (2*W) / (n*(n+1)/2)  where W is the smaller of W+ and W-
            ranks = np.abs(diff).argsort().argsort() + 1
            W_plus = np.sum(ranks[diff > 0])
            W_minus = np.sum(ranks[diff < 0])
            W = min(W_plus, W_minus)
            r = 1 - (2 * W) / (n * (n + 1) / 2)
            # Adjust sign based on direction
            if np.mean(run_data) < np.mean(rest_data):
                r = -r
            
            test_used = 'wilcoxon'
            effect_size = r
        
        return {
            'statistic': stat, 
            'p_value': p, 
            'effect_size': effect_size, 
            'n': n,
            'test_used': test_used
        }
        
    except Exception as e:
        print(f"      Statistical test failed: {e}")
        return {'statistic': np.nan, 'p_value': np.nan, 'effect_size': np.nan, 'n': n, 'test_used': 'error'}


def wilcoxon_test(rest_data, run_data, alternative='two-sided'):
    """
    Legacy wrapper for backward compatibility.
    Calls perform_statistical_test with test_type='wilcoxon'.
    """
    return perform_statistical_test(rest_data, run_data, alternative=alternative, test_type='wilcoxon')


# ==============================================================================
#  HALF-VIOLIN PLOTTING
# ==============================================================================

def plot_half_violin(ax, data_left, data_right, positions, colors, labels, animal_ids=None, show_labels=False):
    """
    Plot half-violin plots with SEPARATED layout: violin → box → dots
    
    Layout for REST (left side):
        [violin extends left] | [box shifted right] | [dots further right]
    
    Layout for RUN (right side):  
        [dots further left] | [box shifted left] | [violin extends right]
    
    Parameters:
        ax: matplotlib axis
        data_left: array for left condition (REST)
        data_right: array for right condition (RUN)
        positions: [left_pos, right_pos] x-positions (base positions)
        colors: [left_color, right_color]
        labels: [left_label, right_label]
        animal_ids: list of animal IDs for labeling points (optional)
        show_labels: if True, add animal ID labels next to REST points
    """
    from scipy.stats import gaussian_kde
    
    pos_left, pos_right = positions
    color_left, color_right = colors
    label_left, label_right = labels
    
    # Calculate positions for separated layout
    # REST: violin on left, box shifted right, dots further right
    violin_base_left = pos_left - 0.05  # Violin anchored here, extends left
    box_pos_left = pos_left + BOX_OFFSET  # Box shifted right of violin
    dot_pos_left = pos_left + DOT_OFFSET  # Dots further right
    
    # RUN: violin on right, box shifted left, dots further left  
    violin_base_right = pos_right + 0.05  # Violin anchored here, extends right
    box_pos_right = pos_right - BOX_OFFSET  # Box shifted left of violin
    dot_pos_right = pos_right - DOT_OFFSET  # Dots further left
    
    # Calculate KDE for violin shape
    def get_violin_path(data, position, side='left', width=VIOLIN_WIDTH):
        if len(data) < 2:
            return None, None, None
        if np.std(data) < 1e-10:  # Handle near-constant data
            return None, None, None
        
        try:
            kde = gaussian_kde(data, bw_method=0.5)
            data_range = np.max(data) - np.min(data)
            padding = max(data_range * 0.2, np.std(data) * 0.3)
            y_range = np.linspace(np.min(data) - padding, 
                                  np.max(data) + padding, 100)
            density = kde(y_range)
            density = density / np.max(density) * width
            
            if side == 'left':
                x_path = position - density
            else:
                x_path = position + density
            
            return y_range, x_path, density
        except:
            return None, None, None
    
    # Plot left half-violin (REST) - extends LEFT from violin_base_left
    result_left = get_violin_path(data_left, violin_base_left, 'left')
    if result_left[0] is not None:
        y_range, x_path, density = result_left
        ax.fill_betweenx(y_range, violin_base_left, x_path, alpha=0.6, color=color_left, label=label_left)
    
    # Plot right half-violin (RUN) - extends RIGHT from violin_base_right
    result_right = get_violin_path(data_right, violin_base_right, 'right')
    if result_right[0] is not None:
        y_range, x_path, density = result_right
        ax.fill_betweenx(y_range, violin_base_right, x_path, alpha=0.6, color=color_right, label=label_right)
    
    # Add box plots - SHIFTED positions (not overlapping violins)
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
    
    # Individual data points in NEAT VERTICAL STRIPS - between the boxes
    ax.scatter([dot_pos_left] * len(data_left), data_left, s=DOT_SIZE, c=[color_left], 
               edgecolor='white', linewidth=1.5, zorder=5, alpha=0.9)
    ax.scatter([dot_pos_right] * len(data_right), data_right, s=DOT_SIZE, c=[color_right],
               edgecolor='white', linewidth=1.5, zorder=5, alpha=0.9)
    
    # Add animal ID labels if provided and show_labels is True
    if show_labels and animal_ids is not None:
        for i, animal_id in enumerate(animal_ids):
            # Add full animal name label to the left of REST point
            ax.annotate(animal_id, (dot_pos_left - 0.12, data_left[i]),
                       fontsize=5, ha='right', va='center', color='gray', alpha=0.9)
    
    # Add paired connecting lines between REST and RUN dots
    for i in range(len(data_left)):
        ax.plot([dot_pos_left, dot_pos_right], 
                [data_left[i], data_right[i]], 
                color='gray', alpha=LINE_ALPHA, linewidth=LINE_WIDTH, zorder=4)


def add_significance_annotation(ax, p_value, y_min, y_max, positions):
    """Add significance stars above the comparison, ensuring no overlap with data."""
    pos_left, pos_right = positions
    x_center = (pos_left + pos_right) / 2
    
    # Position bracket well above the max data point
    y_range = y_max - y_min
    y_line = y_max + y_range * 0.08
    y_text = y_line + y_range * 0.02
    
    # Draw bracket (horizontal line with small vertical ticks)
    bracket_height = y_range * 0.02
    ax.plot([pos_left + 0.05, pos_left + 0.05, pos_right - 0.05, pos_right - 0.05], 
            [y_line - bracket_height, y_line, y_line, y_line - bracket_height],
            color='black', linewidth=1.5, clip_on=False)
    
    # Add stars based on p-value
    if np.isnan(p_value):
        stars = 'N/A'
    elif p_value < 0.001:
        stars = '***'
    elif p_value < 0.01:
        stars = '**'
    elif p_value < 0.05:
        stars = '*'
    else:
        stars = 'n.s.'
    
    ax.text(x_center, y_text, stars, ha='center', va='bottom', fontsize=12, fontweight='bold', 
            style='italic' if stars == 'n.s.' else 'normal')


# ==============================================================================
#  MAIN FIGURE CREATION
# ==============================================================================

def create_fooof_figure(fooof_results_lfp, fooof_results_gevi, output_path, method_title='mscohere'):
    """
    Create publication-ready FOOOF figure with half-violin + paired dot plots.
    
    Layout: 2 rows × 3 columns
        Row 1: LFP - Theta Peak | Aperiodic Exponent | Aperiodic Offset
        Row 2: GEVI - Theta Peak | Aperiodic Exponent | Aperiodic Offset
    """
    
    # LARGER figure size to prevent crowding
    fig, axes = plt.subplots(2, 3, figsize=(18, 14))
    fig.subplots_adjust(left=0.07, right=0.93, top=0.88, bottom=0.06, wspace=0.35, hspace=0.35)
    
    positions = [0.7, 1.3]  # Even wider spacing between REST and RUN
    labels = ['REST', 'RUN']
    
    # Signal-specific colors: LFP = Purple, GEVI = Teal
    colors_lfp = [COLOR_LFP_REST, COLOR_LFP_RUN]
    colors_gevi = [COLOR_GEVI_REST, COLOR_GEVI_RUN]
    
    # Column labels
    col_titles = ['Theta Peak Amplitude', 'Aperiodic Exponent', 'Aperiodic Offset']
    y_labels = ['Power (a.u.)', 'Exponent', 'Offset (log power)']
    
    # Process both signals
    for row_idx, (results, signal_name, colors) in enumerate([
        (fooof_results_lfp, 'LFP', colors_lfp),      # Purple for LFP
        (fooof_results_gevi, 'GEVI', colors_gevi)    # Teal for GEVI
    ]):
        if results is None or len(results['animal_ids']) == 0:
            for col_idx in range(3):
                axes[row_idx, col_idx].text(0.5, 0.5, 'No data', transform=axes[row_idx, col_idx].transAxes,
                                            ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
                axes[row_idx, col_idx].set_title(f'{signal_name}: {col_titles[col_idx]}', fontsize=FONT_SIZE_TITLE, fontweight='bold')
            continue
        
        # Get animal IDs for labeling
        animal_ids = results.get('animal_ids', [])
        
        # Helper to plot a single panel (using signal-specific colors)
        def plot_panel(ax, data_rest, data_run, title, ylabel, stats_result, panel_colors, animal_ids=None):
            if len(data_rest) == 0 or len(data_run) == 0:
                ax.text(0.5, 0.5, 'Insufficient data', transform=ax.transAxes,
                        ha='center', va='center', fontsize=FONT_SIZE_LABEL, color='gray')
                ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
                return
            
            plot_half_violin(ax, data_rest, data_run, positions, panel_colors, labels,
                           animal_ids=animal_ids, show_labels=SHOW_ANIMAL_LABELS)
            
            # Calculate y limits with GENEROUS padding for violin and annotation
            all_data = np.concatenate([data_rest, data_run])
            y_min = np.min(all_data)
            y_max = np.max(all_data)
            y_range = y_max - y_min
            if y_range < 1e-10:  # Handle near-constant data
                y_range = abs(y_max) * 0.2 if y_max != 0 else 1.0
            
            # Set y limits FIRST with generous padding
            # Bottom: 40% padding for violin tails
            # Top: 30% padding for significance annotation
            y_bottom = y_min - y_range * 0.40
            y_top = y_max + y_range * 0.30
            ax.set_ylim([y_bottom, y_top])
            
            # Add significance annotation AFTER setting y limits
            add_significance_annotation(ax, stats_result['p_value'], y_min, y_max, positions)
            
            ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
            ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
            ax.set_xticks(positions)
            ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK)
            ax.set_xlim([0.15, 1.85])  # Even wider x limits for violin tails
            
            # Stats annotation - BOTTOM LEFT corner to avoid violin overlap
            test_name = stats_result.get('test_used', 'unknown')
            test_abbrev = {'paired_ttest': 't-test', 'wilcoxon': 'Wilcoxon', 'none': 'N/A', 'error': 'Error'}.get(test_name, test_name)
            stats_text = f'p={stats_result["p_value"]:.3f}\n{test_abbrev}\nN={stats_result["n"]}'
            # Add effect size with appropriate label
            if test_name == 'paired_ttest':
                stats_text += f'\nd={stats_result["effect_size"]:.2f}'
            else:
                stats_text += f'\nr={stats_result["effect_size"]:.2f}'
            ax.text(0.03, 0.03, stats_text, transform=ax.transAxes, va='bottom', ha='left', fontsize=9, 
                    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, edgecolor='gray', linewidth=0.5))
        
        # Compute stats first so we can pass N to the panel
        # Use one-tailed test for theta (RUN > REST hypothesis) if configured
        theta_alternative = 'greater' if USE_ONE_TAILED_THETA else 'two-sided'
        exp_alternative = 'greater' if USE_ONE_TAILED_EXPONENT else 'two-sided'
        off_alternative = 'greater' if USE_ONE_TAILED_OFFSET else 'two-sided'
        
        stats_theta = perform_statistical_test(
            results['theta_rest'], results['theta_run'], 
            alternative=theta_alternative, test_type=STATISTICAL_TEST
        )
        stats_exp = perform_statistical_test(
            results['exponent_rest'], results['exponent_run'],
            alternative=exp_alternative, test_type=STATISTICAL_TEST
        )
        stats_off = perform_statistical_test(
            results['offset_rest'], results['offset_run'],
            alternative=off_alternative, test_type=STATISTICAL_TEST
        )
        
        # Column 0: Theta Peak Amplitude
        plot_panel(axes[row_idx, 0], results['theta_rest'], results['theta_run'],
                   f'{signal_name}: {col_titles[0]}', y_labels[0], stats_theta, colors, animal_ids)
        
        # Column 1: Aperiodic Exponent
        plot_panel(axes[row_idx, 1], results['exponent_rest'], results['exponent_run'],
                   f'{signal_name}: {col_titles[1]}', y_labels[1], stats_exp, colors, animal_ids)
        
        # Column 2: Aperiodic Offset
        plot_panel(axes[row_idx, 2], results['offset_rest'], results['offset_run'],
                   f'{signal_name}: {col_titles[2]}', y_labels[2], stats_off, colors, animal_ids)
    
    # Style all axes
    for ax in axes.flat:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
        ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
        ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICK,
                       width=TICK_WIDTH, length=TICK_LENGTH)
    
    # Add legend in upper right - show both signal types
    legend_elements = [
        Patch(facecolor=COLOR_LFP_REST, alpha=0.6, label='LFP REST', edgecolor='black', linewidth=1),
        Patch(facecolor=COLOR_LFP_RUN, alpha=0.6, label='LFP RUN', edgecolor='black', linewidth=1),
        Patch(facecolor=COLOR_GEVI_REST, alpha=0.6, label='GEVI REST', edgecolor='black', linewidth=1),
        Patch(facecolor=COLOR_GEVI_RUN, alpha=0.6, label='GEVI RUN', edgecolor='black', linewidth=1),
    ]
    fig.legend(handles=legend_elements, loc='upper right', fontsize=FONT_SIZE_LEGEND, 
               frameon=True, bbox_to_anchor=(0.99, 0.99), framealpha=0.9, edgecolor='gray', ncol=2)
    
    # Determine test name for title
    test_name = 'Auto-selected test' if STATISTICAL_TEST == 'auto' else \
                'Paired t-test' if STATISTICAL_TEST == 'ttest' else 'Wilcoxon signed-rank test'
    fig.suptitle(f'FOOOF Analysis: Oscillatory vs Aperiodic Components\n{method_title} ({test_name})',
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.95)
    
    # Ensure output directory exists (with long path support)
    output_dir = Path(output_path).parent
    force_create_directory(output_dir)
    
    # Save figure
    for fmt in FIGURE_FORMATS:
        save_path = f"{output_path}.{fmt}"
        fig.savefig(to_long_path(save_path), dpi=FIGURE_DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    Saved: {Path(output_path).name}")


def get_pooling_level_label(pooling_level):
    """Convert pooling level string to human-readable label."""
    if pooling_level == 'animal_concatenated':
        return 'Animal-Concatenated'
    elif pooling_level == 'animal_pooled':
        return 'Animal-Pooled'
    else:
        return pooling_level.replace('_', '-').title()


def run_fooof_pipeline(animals_to_process, output_dir, method='mscohere', pooling_level=None):
    """
    Run the complete FOOOF analysis pipeline.
    
    NOTE: The 'method' parameter specifies which data file to load (mscohere.mat or fieldtrip.mat).
    Since PSD is computed identically for both methods (using pwelch), the results will be the same.
    We default to 'mscohere' but either would work.
    
    Parameters:
        animals_to_process: list of animal dictionaries
        output_dir: Path object for output directory
        method: 'mscohere' or 'fieldtrip'
        pooling_level: 'animal_pooled' or 'animal_concatenated' (defaults to GROUP_POOLING_LEVEL)
    
    Returns number of figures generated.
    """
    if not FOOOF_AVAILABLE:
        print("    FOOOF not available. Install with: pip install fooof")
        return 0
    
    # Use the specified pooling level, defaulting to GROUP_POOLING_LEVEL
    if pooling_level is None:
        pooling_level = GROUP_POOLING_LEVEL
    
    pooling_label = get_pooling_level_label(pooling_level)
    
    print(f"\n    Loading animal PSD data ({pooling_label})...")
    animal_data = load_animal_psd_data(animals_to_process, method, pooling_level)
    
    if len(animal_data) == 0:
        print("    No animal data loaded, skipping FOOOF analysis")
        return 0
    
    print(f"\n    Running FOOOF analysis on LFP...")
    fooof_results_lfp = run_fooof_analysis(animal_data, signal_type='lfp')
    
    # Print LFP FOOOF values for diagnostics
    if fooof_results_lfp and len(fooof_results_lfp['animal_ids']) > 0:
        print(f"\n    LFP FOOOF Results (REST → RUN):")
        print(f"    " + "-" * 85)
        print(f"    {'Animal':<15s} | {'Theta Peak (a.u.)':<20s} | {'Aperiodic Offset':<20s} | {'Aperiodic Exponent':<18s}")
        print(f"    {'':<15s} | {'REST→RUN':<20s} | {'REST→RUN':<20s} | {'REST→RUN':<18s}")
        print(f"    " + "-" * 85)
        for i, animal_id in enumerate(fooof_results_lfp['animal_ids']):
            # Theta
            theta_rest = fooof_results_lfp['theta_rest'][i]
            theta_run = fooof_results_lfp['theta_run'][i]
            theta_diff = theta_run - theta_rest
            theta_dir = "↑" if theta_diff > 0 else "↓"
            
            # Offset (higher offset = higher overall power)
            off_rest = fooof_results_lfp['offset_rest'][i]
            off_run = fooof_results_lfp['offset_run'][i]
            off_diff = off_run - off_rest
            off_dir = "↑" if off_diff > 0 else "↓"
            
            # Exponent (higher exponent = steeper 1/f slope)
            exp_rest = fooof_results_lfp['exponent_rest'][i]
            exp_run = fooof_results_lfp['exponent_run'][i]
            exp_diff = exp_run - exp_rest
            exp_dir = "↑" if exp_diff > 0 else "↓"
            
            flag = " ***" if theta_diff < 0 else ""
            print(f"    {animal_id:<15s} | {theta_rest:.3f}→{theta_run:.3f} ({theta_diff:+.2f}{theta_dir}){flag:<3s} | "
                  f"{off_rest:.2f}→{off_run:.2f} ({off_diff:+.2f}{off_dir}) | "
                  f"{exp_rest:.2f}→{exp_run:.2f} ({exp_diff:+.2f}{exp_dir})")
        print(f"    " + "-" * 85)
        print(f"    Note: Theta Peak = oscillatory power ABOVE 1/f, not raw PSD power")
        print(f"          Higher offset during RUN but lower theta peak means 1/f shifted up")
        print(f"          but the actual theta oscillation (above 1/f) decreased.")
    
    print(f"\n    Running FOOOF analysis on GEVI...")
    fooof_results_gevi = run_fooof_analysis(animal_data, signal_type='gevi')
    
    # Create figure - use short pooling level abbreviation in filename
    # animal_concatenated -> animal_concat, animal_pooled -> animal_pooled
    level_abbrev = 'animal_concat' if pooling_level == 'animal_concatenated' else 'animal_pooled'
    output_path = output_dir / f"fooof_{level_abbrev}"
    create_fooof_figure(fooof_results_lfp, fooof_results_gevi, output_path, 
                        method_title=f'PSD Analysis ({pooling_label})')
    
    return 1


# ==============================================================================
#  STANDALONE EXECUTION
# ==============================================================================

def main():
    """Standalone execution for testing."""
    print("=" * 70)
    print("  FOOOF ANALYSIS - Oscillatory Decomposition")
    print("=" * 70)
    print("\n  NOTE: FOOOF analyzes PSD (computed identically for mscohere/FieldTrip)")
    print("        Running once is sufficient.\n")
    
    if not FOOOF_AVAILABLE:
        print("\nERROR: FOOOF not installed. Install with: pip install fooof")
        return
    
    # Get animals from config
    if USING_CENTRAL_CONFIG:
        animals = get_animals_to_process()
        output_dir = get_group_level_output_dir()
        pooling_level = GROUP_POOLING_LEVEL
    else:
        print("ERROR: Central config not available")
        return
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    pooling_label = get_pooling_level_label(pooling_level)
    
    print(f"  Output: {output_dir}")
    print(f"  Animals: {[a['mouse_id'] for a in animals]}")
    print(f"  Pooling Level: {pooling_label}")
    
    # Run FOOOF (only need to run once since PSD is identical for both methods)
    run_fooof_pipeline(animals, output_dir, method='mscohere', pooling_level=pooling_level)
    
    print("\n" + "=" * 70)
    print("  FOOOF ANALYSIS COMPLETE")
    print("=" * 70)


if __name__ == '__main__':
    main()
