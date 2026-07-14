"""
Plot comprehensive stimulation analysis figures.

Creates 12 separate figures:
1. Traces (R1): Single trial, zoomed (−0.5 s to stim+0.5 s), trial-averaged
2. Stim-onset zoom (R1): −0.5 s to +0.5 s from onset
3. Traces (R2): Same as (1) for second condition
4. Stim-onset zoom (R2)
5–6. Spectral heatmaps (R1 / R2): Speed, LFP spec, Fiber spec, Coherence
7–8. Period violin plots (R1 / R2): Fiber, LFP power, Fiber power, Coherence
9. 40Hz vs 135Hz comparison: Fiber, Fiber power, Coherence
10–11. Trial-by-trial heatmaps (R1 / R2): Baseline-subtracted and band-pass filtered power
12. Trial-averaged LFP and GEVI: R1 vs R2 overlaid on shared axes (same format as trace figure averages)

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from scipy.io import loadmat
from scipy import signal
from scipy.signal import butter, filtfilt
from scipy import stats
from scipy.ndimage import uniform_filter1d
from pathlib import Path
import warnings
import h5py
import os
import re
import sys
import argparse

warnings.filterwarnings('ignore')

# common.py lives in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# =============================================================================
# CONFIGURATION - CHANGE THESE FOR DIFFERENT COMPARISONS
# =============================================================================

# ============== COMPARISON SELECTOR ==============
# Set COMPARISON_MODE to select which analysis to run:
#   'Animal01_40vs135'          - Animal01: 40Hz vs 135Hz (original, 1s stim)
#   'Animal02_AmpBalanced'      - Animal02: 135Hz vs 40Hz (Amplitude balanced, 1s stim)
#   'Animal02_EnergyBalanced'   - Animal02: 135Hz vs 40Hz (Energy balanced, 1s stim)
#   'Animal03_AmpBalanced'      - Animal03: 135Hz vs 40Hz (Amplitude balanced, 10s stim)
#   'Animal03_EnergyBalanced'   - Animal03: 135Hz vs 40Hz (Energy balanced, 10s stim)
#   'Animal04_AmpBalanced'      - Animal04: 135Hz vs 40Hz (Amplitude balanced, 10s stim)
#   'Animal04_EnergyBalanced'   - Animal04: 135Hz vs 40Hz (Energy balanced, 10s stim)
COMPARISON_MODE = 'Animal04_EnergyBalanced'  # <-- CHANGE THIS TO SWITCH COMPARISONS

# ============== AUTO-CONFIGURED BASED ON COMPARISON_MODE ==============
if COMPARISON_MODE == 'Animal01_40vs135':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal01" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_09_25-R1"  # 40Hz stimulation
    SESSION_R2 = "01_09_25-R2"  # 135Hz stimulation
    MOUSE_ID = "Animal01"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0   # Stimulation frequency for R1
    FREQ_R2 = 135.0  # Stimulation frequency for R2
    STIM_DURATION_OVERRIDE = None  # Use default 1s
    
elif COMPARISON_MODE == 'Animal02_AmpBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal02" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_02_26-R9"   # 40Hz (Amplitude balanced)
    SESSION_R2 = "01_02_26-R6"   # 135Hz
    MOUSE_ID = "Animal02"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Amp-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = None  # Use default 1s
    
elif COMPARISON_MODE == 'Animal02_EnergyBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal02" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_02_26-R10"  # 40Hz (Energy balanced)
    SESSION_R2 = "01_02_26-R6"   # 135Hz
    MOUSE_ID = "Animal02"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Energy-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = None  # Use default 1s

elif COMPARISON_MODE == 'Animal03_AmpBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal03" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_03_26-R5"   # 40Hz @ 2.9V (Amplitude balanced)
    SESSION_R2 = "01_03_26-R3"   # 135Hz @ 2.9V
    MOUSE_ID = "Animal03"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Amp-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = 10.0  # 10 second stimulation

elif COMPARISON_MODE == 'Animal03_EnergyBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal03" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_03_26-R6"   # 40Hz @ 4.5V (Energy balanced)
    SESSION_R2 = "01_03_26-R3"   # 135Hz @ 2.9V
    MOUSE_ID = "Animal03"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Energy-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = 10.0  # 10 second stimulation

elif COMPARISON_MODE == 'Animal04_AmpBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal04" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_04_26-R10"   # 40Hz (Amplitude balanced)
    SESSION_R2 = "01_04_26-R7"    # 135Hz
    MOUSE_ID = "Animal04"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Amp-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = 10.0  # 10 second stimulation

elif COMPARISON_MODE == 'Animal04_EnergyBalanced':
    BASE_PATH = _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal04" / "Fiber_Voltage_Processed"
    SESSION_R1 = "01_04_26-R11"   # 40Hz (Energy balanced)
    SESSION_R2 = "01_04_26-R7"    # 135Hz
    MOUSE_ID = "Animal04"
    NUM_TRIALS = 10
    LABEL_R1 = "40Hz (Energy-balanced)"
    LABEL_R2 = "135Hz"
    FREQ_R1 = 40.0
    FREQ_R2 = 135.0
    STIM_DURATION_OVERRIDE = 10.0  # 10 second stimulation

else:
    raise ValueError(f"Unknown COMPARISON_MODE: {COMPARISON_MODE}")

# Spectral outputs from MATLAB pipeline
SPECTRAL_OUTPUT_ROOT = str(PROJECT_ROOT / "Figures" / "Stimulation_analysis" / "Spectral_data_outputs")
SPECTRAL_METHOD = "fieldtrip"  # or "mscohere"

# Stimulation parameters - adjusted based on comparison mode
if COMPARISON_MODE.startswith('Animal03') or COMPARISON_MODE.startswith('Animal04'):
    # These example animals used 10s pre-stim, 10s stim, 10s post-stim
    PRE_STIM_DURATION = 10.0   # seconds
    STIM_DURATION = 10.0       # seconds
    POST_STIM_DURATION = 10.0  # seconds
    STIM_TRANSIENT_END = 1.0   # seconds (transient period: 0-1s for longer stim)
    STIM_SUSTAINED_START = 1.0 # seconds (sustained period starts after transient)
else:
    # Default for Animal01, Animal02 (4s pre, 1s stim, 5s post)
    PRE_STIM_DURATION = 4.0    # seconds
    STIM_DURATION = STIM_DURATION_OVERRIDE if STIM_DURATION_OVERRIDE is not None else 1.0  # seconds
    POST_STIM_DURATION = 5.0   # seconds
    STIM_TRANSIENT_END = 0.15  # seconds (transient period: 0-0.15s)
    STIM_SUSTAINED_START = 0.15  # seconds (sustained period starts after transient)

# USER CONFIGURABLE: Skip first N seconds of pre-stim (e.g., if noisy)
# Set to 0 to use full pre-stim period, set to 1.0 to skip first second
if COMPARISON_MODE.startswith('Animal03') or COMPARISON_MODE.startswith('Animal04'):
    PRE_STIM_SKIP_SEC = 0.0  # Animal03/Animal04: use full pre-stim (longer baseline available)
else:
    PRE_STIM_SKIP_SEC = 1.0  # Animal01/Animal02: Skip first 1 second (noise issue)

# Effective pre-stim duration for plotting
PRE_STIM_EFFECTIVE = PRE_STIM_DURATION - PRE_STIM_SKIP_SEC

# Representative trial for single-trial plots
REPRESENTATIVE_TRIAL = 1

# Figure configuration - INCREASED for publication
DPI = 300
FONT_SIZE = 14
FONT_SIZE_TITLE = 18
FONT_SIZE_SUPTITLE = 22
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_SCALEBAR = 13

# Axis styling - INCREASED for publication
AXIS_LINEWIDTH = 2.5
TICK_WIDTH = 2.0
TICK_LENGTH = 8
LINE_WIDTH_TRACE = 1.2
LINE_WIDTH_THICK = 2.0

# Colors (consistent with fig1_gevi_lfp.py)
COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])  # teal
COLOR_LFP = np.array([0.35, 0.25, 0.45])  # purple-grey
COLOR_MOTION = np.array([0.993248, 0.7, 0.4])  # orange

# Period colors - UPDATED per user request
COLOR_PRE = np.array([0.5, 0.5, 0.5])  # grey
COLOR_TRANSIENT = np.array([0.7, 0.15, 0.15])  # deep red
COLOR_SUSTAINED = np.array([0.65, 0.25, 0.15])  # maroon-orange
COLOR_POST = np.array([0.15, 0.55, 0.55])  # teal-ish
COLOR_STIM_PULSE = np.array([0.5, 0.1, 0.1])  # dark red for stimulation pulses

# Frequency comparison colors - UPDATED per user request
COLOR_40HZ = np.array([0.15, 0.55, 0.55])  # teal-ish for 40Hz
COLOR_40HZ_LIGHT = np.array([0.35, 0.70, 0.70])  # lighter teal
COLOR_135HZ = np.array([0.95, 0.65, 0.45])  # peach for 135Hz
COLOR_135HZ_LIGHT = np.array([0.98, 0.80, 0.60])  # lighter peach

# Scale bars
SCALEBAR_GEVI_PERCENT = 1.0
SCALEBAR_LFP_UV = 500.0
SCALEBAR_MOTION_CM_S = 5.0  # Reduced from 15 to avoid overshadowing fiber scale bar
SCALEBAR_LINEWIDTH = 3.5

# Separate stim-onset zoom figures only: time from stim onset (s), symmetric window around t=0
STIM_ONSET_ZOOM_PRE_S = 0.25   # seconds before stim onset
STIM_ONSET_ZOOM_POST_S = 0.4  # seconds after stim onset

# LFP already in µV (not Volts)
LFP_IN_VOLTS = False

# Motion conversion
WHEEL_DIAMETER_CM = 19.0
WHEEL_CIRCUMFERENCE_CM = np.pi * WHEEL_DIAMETER_CM
ENCODER_COUNTS_PER_REV = 1024
EPHYS_SAMPLING_RATE = 30000
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM

# Frequency range for spectrograms
FREQ_RANGE = (1, 150)  # Extended to capture 135 Hz stimulation

# Spectral heatmaps: baseline = mean over spectrogram time bins in [LO, HI) seconds relative to stim onset.
# Default [-3, -1): avoids the last second before onset (possible edge/leakage/anticipation) and never includes
# stimulation (HI < 0). If no bins fall in this window, code falls back to all pre-onset bins (time < 0).
SPECTROGRAM_BASELINE_T_LO = -3.0
SPECTROGRAM_BASELINE_T_HI = -1.0

# Custom colormap for speed heatmap (consistent with run_all_plots.py)
from common import _infer_trial_from_name, create_monochromatic_orange_cmap, create_parula_like_cmap, normalize_unc_path  # shared helpers (were local copies)

CMAP_SPEED = create_monochromatic_orange_cmap()


def create_yellow_purple_diverging_cmap():
    """
    Create a smooth yellow-purple diverging colormap with no white at center.
    Yellow for positive values, purple for negative values.
    """
    colors = [
        (0.4, 0.0, 0.5),       # Dark purple (negative extreme)
        (0.5, 0.2, 0.6),       # Purple
        (0.6, 0.4, 0.7),       # Light purple
        (0.7, 0.6, 0.8),       # Lavender
        (0.75, 0.7, 0.85),     # Very light purple (near center, no white)
        (0.8, 0.75, 0.7),      # Very light yellow (near center, no white)
        (0.85, 0.8, 0.5),      # Light yellow
        (0.9, 0.85, 0.3),      # Yellow
        (0.95, 0.9, 0.1),      # Bright yellow
        (1.0, 0.95, 0.0),      # Pure yellow (positive extreme)
    ]
    return LinearSegmentedColormap.from_list('yellow_purple', colors, N=256)

CMAP_YELLOW_PURPLE = create_yellow_purple_diverging_cmap()



CMAP_PARULA_LIKE = create_parula_like_cmap()

# Output
OUTPUT_DIR = PROJECT_ROOT / "Figures" / "Stimulation_analysis"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def get_condition_name(session_id):
    """Get condition name for a session based on COMPARISON_MODE."""
    # Map session IDs to condition names based on comparison mode
    if COMPARISON_MODE == 'Animal01_40vs135':
        if session_id == '01_09_25-R1':
            return '40Hz'
        elif session_id == '01_09_25-R2':
            return '135Hz'
    elif COMPARISON_MODE == 'Animal02_AmpBalanced':
        if session_id == '01_02_26-R6':
            return '135Hz'
        elif session_id == '01_02_26-R9':
            return '40Hz_AmpBalanced'
    elif COMPARISON_MODE == 'Animal02_EnergyBalanced':
        if session_id == '01_02_26-R6':
            return '135Hz'
        elif session_id == '01_02_26-R10':
            return '40Hz_EnergyBalanced'
    elif COMPARISON_MODE == 'Animal03_AmpBalanced':
        if session_id == '01_03_26-R3':
            return '135Hz'
        elif session_id == '01_03_26-R5':
            return '40Hz_AmpBalanced'
    elif COMPARISON_MODE == 'Animal03_EnergyBalanced':
        if session_id == '01_03_26-R3':
            return '135Hz'
        elif session_id == '01_03_26-R6':
            return '40Hz_EnergyBalanced'
    elif COMPARISON_MODE == 'Animal04_AmpBalanced':
        if session_id == '01_04_26-R7':
            return '135Hz'
        elif session_id == '01_04_26-R10':
            return '40Hz_AmpBalanced'
    elif COMPARISON_MODE == 'Animal04_EnergyBalanced':
        if session_id == '01_04_26-R7':
            return '135Hz'
        elif session_id == '01_04_26-R11':
            return '40Hz_EnergyBalanced'
    
    # Fallback
    return 'unknown'


def build_spectral_path(session_id, trial_num, method):
    """Build path to spectral results file."""
    condition = get_condition_name(session_id)
    
    spectral_dir = os.path.join(
        normalize_unc_path(SPECTRAL_OUTPUT_ROOT), 
        MOUSE_ID, 
        f"{session_id}_{condition}"
    )
    filename = f"{MOUSE_ID}_{session_id}_Trial{trial_num}_{method}_SpectralResults.mat"
    return os.path.join(spectral_dir, filename)


# =============================================================================
# DATA LOADING FUNCTIONS
# =============================================================================

def matlab_struct_to_dict(matobj):
    """Recursively convert MATLAB structs to nested Python dicts."""
    from scipy.io.matlab import mat_struct
    
    # Handle mat_struct objects (from struct_as_record=False)
    if isinstance(matobj, mat_struct):
        out = {}
        for field_name in matobj._fieldnames:
            out[field_name] = matlab_struct_to_dict(getattr(matobj, field_name))
        return out
    
    if isinstance(matobj, np.ndarray) and matobj.dtype == np.object_ and matobj.size == 1:
        matobj = matobj.item()
    if isinstance(matobj, np.void):
        out = {}
        for field_name in matobj.dtype.names:
            out[field_name] = matlab_struct_to_dict(matobj[field_name])
        return out
    if isinstance(matobj, np.ndarray) and matobj.dtype == np.object_:
        return [matlab_struct_to_dict(el) for el in matobj]
    return matobj


def load_raw_trial_data(session_id, trial_num, fiber_index=0):
    """Load raw preprocessed trial data for traces."""
    base_path = BASE_PATH / session_id
    
    # Find trial directory
    trial_dirs = list(base_path.glob(f"Trial{trial_num}_*"))
    if not trial_dirs:
        raise FileNotFoundError(f"Trial {trial_num} directory not found for {session_id}")
    
    trial_dir = trial_dirs[0]
    mat_file = trial_dir / f"{MOUSE_ID}-{session_id}_Trial{trial_num}_FiberPhotometry_Analysis.mat"
    mat_file_str = normalize_unc_path(str(mat_file))
    
    if not os.path.exists(mat_file_str):
        raise FileNotFoundError(f"MAT file not found: {mat_file}")
    
    # Try HDF5 first (MATLAB -v7.3 format) - more reliable
    try:
        with h5py.File(mat_file_str, "r") as f:
            root = f["FiberPhotometryAnalysis"]
            
            # Time vector - handle both row and column vectors
            t_raw = np.array(root["time"]["time_vector_seconds"][()])
            t = t_raw.flatten()
            
            sr = np.array(root["time"]["sampling_rate"][()])
            fs = float(sr.item() if sr.size == 1 else sr.ravel()[0])
            
            # GEVI traces - handle transposed arrays
            gevi_all = np.array(root["signals"]["final_processed_traces"][()])
            if gevi_all.ndim == 2:
                # MATLAB stores column-major, so [N x fibers] in MATLAB becomes [fibers x N] in Python
                if gevi_all.shape[0] < gevi_all.shape[1]:
                    gevi_all = gevi_all.T  # Transpose to [N x fibers]
                gevi = gevi_all[:, fiber_index] if gevi_all.shape[1] > fiber_index else gevi_all.flatten()
            else:
                gevi = gevi_all.flatten()
            
            ephys_grp = root["ephys"]
            if "lfp_raw_aligned_HP" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_HP"][()]).flatten()
            elif "lfp_raw_aligned_mPFC" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_mPFC"][()]).flatten()
            else:
                raise KeyError("No LFP trace found")
            
            if "running_velocity_smooth" in ephys_grp:
                motion = np.array(ephys_grp["running_velocity_smooth"][()]).flatten()
            else:
                motion = np.array(ephys_grp["running_velocity"][()]).flatten()
                
    except (OSError, KeyError) as e:
        # Fall back to scipy loadmat for older MAT format
        mat = loadmat(mat_file_str, squeeze_me=False, struct_as_record=False)
        fp_struct = mat["FiberPhotometryAnalysis"]
        
        # Extract with careful handling of nested structs
        t = np.asarray(fp_struct.time.time_vector_seconds).flatten()
        try:
            fs = float(np.asarray(fp_struct.time.sampling_rate).flatten()[0])
        except:
            fs = 1.0 / np.median(np.diff(t))
        
        gevi_all = np.asarray(fp_struct.signals.final_processed_traces)
        if gevi_all.ndim == 1:
            gevi = gevi_all.flatten()
        else:
            if gevi_all.shape[0] < gevi_all.shape[1]:
                gevi_all = gevi_all.T
            gevi = gevi_all[:, fiber_index] if gevi_all.ndim > 1 and gevi_all.shape[1] > fiber_index else gevi_all.flatten()
        
        ephys = fp_struct.ephys
        if hasattr(ephys, 'lfp_raw_aligned_HP'):
            lfp = np.asarray(ephys.lfp_raw_aligned_HP).flatten()
        elif hasattr(ephys, 'lfp_raw_aligned_mPFC'):
            lfp = np.asarray(ephys.lfp_raw_aligned_mPFC).flatten()
        else:
            raise KeyError("No LFP trace found")
        
        if hasattr(ephys, 'running_velocity_smooth'):
            motion = np.asarray(ephys.running_velocity_smooth).flatten()
        else:
            motion = np.asarray(ephys.running_velocity).flatten()
    
    n = min(len(t), len(gevi), len(lfp), len(motion))
    
    # Validate data
    if n < 100:
        raise ValueError(f"Data too short: only {n} samples. t:{len(t)}, gevi:{len(gevi)}, lfp:{len(lfp)}, motion:{len(motion)}")
    
    # Debug: print loaded data info
    print(f"      [Loaded: {n} samples, t=[{t[0]:.2f}, {t[-1]:.2f}]s, fs={fs:.1f}Hz]")
    
    return {
        't': t[:n], 
        'gevi': gevi[:n], 
        'lfp': lfp[:n],
        'motion': motion[:n] * MOTION_TO_CM_PER_S, 
        'fs': fs
    }


def load_spectral_results(session_id, trial_num, method=SPECTRAL_METHOD):
    """Load spectral results from MATLAB pipeline output."""
    spectral_file_raw = build_spectral_path(session_id, trial_num, method)
    # Convert to extended path format to bypass 260 char limit
    spectral_file = normalize_unc_path(spectral_file_raw, for_access=True)
    spectral_dir = normalize_unc_path(os.path.dirname(spectral_file_raw), for_access=True)
    filename = os.path.basename(spectral_file_raw)
    
    def load_hdf5_group(grp):
        """Recursively load HDF5 group into dict."""
        out = {}
        for key in grp.keys():
            item = grp[key]
            if isinstance(item, h5py.Group):
                out[key] = load_hdf5_group(item)
            elif isinstance(item, h5py.Dataset):
                data = item[()]
                if isinstance(data, bytes):
                    out[key] = data.decode('utf-8')
                elif isinstance(data, np.ndarray) and data.dtype.kind == 'S':
                    out[key] = ''.join([chr(c) for c in data.flatten()])
                elif data.ndim == 0:
                    out[key] = data.item()
                else:
                    out[key] = np.array(data).squeeze()
            else:
                out[key] = item
        return out
    
    # Try HDF5 first (MATLAB -v7.3)
    try:
        with h5py.File(spectral_file, "r") as f:
            root = f["StimSpectralResults"]
            return load_hdf5_group(root)
    except OSError as hdf5_err:
        pass  # Not HDF5 format or file not found
    except Exception as hdf5_err:
        pass
    
    # Try scipy.io.loadmat (MATLAB -v7)
    try:
        mat = loadmat(spectral_file, squeeze_me=True, struct_as_record=False)
        return matlab_struct_to_dict(mat["StimSpectralResults"])
    except FileNotFoundError:
        pass
    except Exception as scipy_err:
        pass
    
    # Both failed - provide helpful error message
    try:
        existing_files = os.listdir(spectral_dir)
        if filename in existing_files:
            raise IOError(f"File exists but could not be read: {spectral_file_raw}\n  (Path length: {len(spectral_file_raw)} chars)")
        raise FileNotFoundError(f"Spectral results not found: {spectral_file_raw}\n  Dir has: {existing_files[:5]}...")
    except FileNotFoundError:
        raise
    except IOError:
        raise
    except Exception as e:
        raise FileNotFoundError(f"Could not access spectral directory: {spectral_dir} - {e}")


def load_all_spectral_results(session_id, num_trials=NUM_TRIALS, method=SPECTRAL_METHOD):
    """Load spectral results for all trials of a session."""
    results = []
    for trial_num in range(1, num_trials + 1):
        try:
            result = load_spectral_results(session_id, trial_num, method)
            results.append(result)
            print(f"    Loaded spectral: {session_id} Trial {trial_num}")
        except Exception as e:
            print(f"    Warning: {session_id} Trial {trial_num}: {e}")
            results.append(None)
    return results


def load_all_raw_trials(session_id, num_trials=NUM_TRIALS):
    """Load all raw trial data for a session."""
    trials = []
    for trial_num in range(1, num_trials + 1):
        try:
            data = load_raw_trial_data(session_id, trial_num)
            trials.append(data)
            print(f"    Loaded raw: {session_id} Trial {trial_num}")
        except Exception as e:
            print(f"    Warning: {session_id} Trial {trial_num}: {e}")
            trials.append(None)
    return trials


# =============================================================================
# PLOTTING HELPERS
# =============================================================================

def vbar_proportional(ax, x_frac, height_value, label, side="left", color="k"):
    """Draw a vertical scale bar with publication-ready formatting.
    
    Scale bar is positioned OUTSIDE the plot area on the left side with label to its left.
    """
    y_min, y_max = ax.get_ylim()
    y_range = y_max - y_min
    if y_range == 0:
        return
    
    # Calculate bar height as fraction of y-range
    bar_frac = height_value / y_range
    # Center the bar vertically
    y0 = y_min + (0.5 - bar_frac/2) * y_range
    y1 = y_min + (0.5 + bar_frac/2) * y_range
    
    x_min, x_max = ax.get_xlim()
    x_range = x_max - x_min
    
    # Position scale bar OUTSIDE the plot on the left (in data coordinates but clipped off)
    # Using negative offset from x_min to place it outside
    x_bar = x_min - 0.02 * x_range  # Bar position outside left edge
    x_text = x_min - 0.04 * x_range  # Text to the left of bar
    
    # Draw the vertical bar
    ax.plot([x_bar, x_bar], [y0, y1], color=color, linewidth=SCALEBAR_LINEWIDTH, 
            clip_on=False, solid_capstyle='butt')
    
    # Add label to the left of the bar, rotated 90 degrees
    ax.text(x_text, (y0 + y1)/2, label, ha="right", va="center", 
            fontsize=FONT_SIZE_SCALEBAR, rotation=90, color=color, clip_on=False,
            fontweight='bold')


def add_significance_bracket(ax, x1, x2, y, p_value, line_height=0.015, text_offset=0.005):
    """
    Add significance bracket with asterisks between two groups.
    
    Style: Compact brackets with stars closely positioned above the horizontal line.
    Consistent style for all violin plots.
    
    Args:
        ax: Matplotlib axis
        x1, x2: X positions of the two groups
        y: Y position for the bracket (in data coordinates)
        p_value: P-value for significance
        line_height: Height of bracket vertical lines (in data coordinates)
        text_offset: Offset for text above bracket (in data coordinates)
    """
    # Determine significance level
    if p_value < 0.001:
        sig_text = '***'
    elif p_value < 0.01:
        sig_text = '**'
    elif p_value < 0.05:
        sig_text = '*'
    else:
        sig_text = 'ns'
    
    # Draw horizontal line (thinner for cleaner look)
    ax.plot([x1, x2], [y, y], 'k-', linewidth=1.5, clip_on=False)
    
    # Draw short vertical lines at ends
    ax.plot([x1, x1], [y - line_height, y], 'k-', linewidth=1.5, clip_on=False)
    ax.plot([x2, x2], [y - line_height, y], 'k-', linewidth=1.5, clip_on=False)
    
    # Add text - positioned CLOSELY above the horizontal line (smaller font for compact look)
    ax.text((x1 + x2) / 2, y + text_offset, sig_text, 
            ha='center', va='bottom', fontsize=FONT_SIZE_TICK - 1, fontweight='bold')


def perform_statistical_test(data1, data2, paired=False):
    """
    Perform appropriate statistical test between two groups.
    
    Returns:
        (p_value, test_name, stat, n1, n2)
    """
    d1 = np.array([v for v in data1 if v is not None and np.isfinite(v)])
    d2 = np.array([v for v in data2 if v is not None and np.isfinite(v)])
    
    if len(d1) < 3 or len(d2) < 3:
        return 1.0, 'insufficient_data', float('nan'), len(d1), len(d2)
    
    try:
        _, p_norm1 = stats.shapiro(d1) if len(d1) <= 5000 else (None, 0.05)
        _, p_norm2 = stats.shapiro(d2) if len(d2) <= 5000 else (None, 0.05)
        normal = (p_norm1 > 0.05) and (p_norm2 > 0.05)
    except:
        normal = False
    
    if paired:
        if normal and len(d1) == len(d2):
            stat, p_value = stats.ttest_rel(d1, d2)
            test_name = 'paired_ttest'
        else:
            if len(d1) == len(d2):
                stat, p_value = stats.wilcoxon(d1, d2, alternative='two-sided')
                test_name = 'wilcoxon'
            else:
                if normal:
                    stat, p_value = stats.ttest_ind(d1, d2)
                    test_name = 'ttest_ind'
                else:
                    stat, p_value = stats.mannwhitneyu(d1, d2, alternative='two-sided')
                    test_name = 'mannwhitney'
    else:
        if normal:
            stat, p_value = stats.ttest_ind(d1, d2)
            test_name = 'ttest_ind'
        else:
            stat, p_value = stats.mannwhitneyu(d1, d2, alternative='two-sided')
            test_name = 'mannwhitney'
    
    return p_value, test_name, stat, len(d1), len(d2)


def holm_bonferroni_correction(p_values):
    """Apply Holm-Bonferroni step-down correction to a list of p-values."""
    n = len(p_values)
    if n == 0:
        return []
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    corrected = [None] * n
    cumulative_max = 0.0
    for rank, (orig_idx, p) in enumerate(indexed):
        adjusted = p * (n - rank)
        adjusted = min(adjusted, 1.0)
        cumulative_max = max(cumulative_max, adjusted)
        corrected[orig_idx] = cumulative_max
    return corrected


def perform_omnibus_repeated_measures(data_dict):
    """
    Friedman test for repeated-measures design (non-parametric).
    Expects data_dict with >=3 groups of equal length.
    Returns (chi2, p_value, k, n) or None on failure.
    """
    arrays = []
    for vals in data_dict.values():
        clean = [v for v in vals if v is not None and np.isfinite(v)]
        arrays.append(clean)
    lengths = [len(a) for a in arrays]
    if len(arrays) < 3 or min(lengths) < 3:
        return None
    min_len = min(lengths)
    trimmed = [a[:min_len] for a in arrays]
    try:
        chi2, p_value = stats.friedmanchisquare(*trimmed)
        return chi2, p_value, len(trimmed), min_len
    except Exception:
        return None


def plot_violin_box(ax, data_dict, ylabel, colors_dict, title=None, comparisons=None,
                    correction='holm', omnibus=False):
    """
    Create publication-ready full violin-box plot with individual trial dots.
    Prints verbose statistical output to terminal.
    """
    positions = []
    data_list = []
    labels = []
    colors_list = []
    label_to_pos = {}  # Map label to position for comparisons
    
    for i, (label, values) in enumerate(data_dict.items()):
        valid_values = [v for v in values if v is not None and np.isfinite(v)]
        if len(valid_values) > 0:
            positions.append(i)
            data_list.append(valid_values)
            labels.append(label)
            colors_list.append(colors_dict.get(label, [0.5, 0.5, 0.5]))
            label_to_pos[label] = i
    
    if len(data_list) == 0:
        ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes, fontsize=FONT_SIZE_LABEL)
        if title:
            ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')
        return
    
    # Full violin plots with transparency and PROMINENT edges
    parts = ax.violinplot(data_list, positions=positions, widths=0.75, 
                          showmeans=False, showmedians=False, showextrema=False)
    for pc, color in zip(parts['bodies'], colors_list):
        pc.set_facecolor(color)
        pc.set_alpha(0.35)
        # Darker, thicker edge for better visibility
        edge_color = np.clip(np.array(color) * 0.3, 0, 1)  # Even darker edge
        pc.set_edgecolor(edge_color)
        pc.set_linewidth(3.5)  # Thicker edge
    
    # Box plots - FILLED with lighter color, dark edge
    bp = ax.boxplot(data_list, positions=positions, widths=0.22, patch_artist=True, showfliers=False)
    for patch, color in zip(bp['boxes'], colors_list):
        # Lighter fill color
        light_color = np.clip(np.array(color) * 0.6 + 0.4, 0, 1)
        patch.set_facecolor(light_color)
        patch.set_edgecolor(np.clip(np.array(color) * 0.4, 0, 1))  # Dark edge
        patch.set_linewidth(2.5)
        patch.set_alpha(0.9)
    for whisker in bp['whiskers']:
        whisker.set_color('black')
        whisker.set_linewidth(2.0)
    for median in bp['medians']:
        median.set_color('black')
        median.set_linewidth(3.0)
    for cap in bp['caps']:
        cap.set_color('black')
        cap.set_linewidth(2.0)
    
    # Individual trial dots - REDUCED SIZE with PROMINENT black borders
    for i, (pos, values) in enumerate(zip(positions, data_list)):
        jitter = np.random.uniform(-0.10, 0.10, len(values))
        ax.scatter(pos + jitter, values, 
                   color=colors_list[i], s=70, alpha=0.9, zorder=10,  # Reduced from 100 to 70
                   edgecolors='black', linewidths=2.0)
    
    # Omnibus test (Friedman) if requested
    panel_label = title if title else ylabel
    if omnibus and len(data_list) >= 3:
        omni = perform_omnibus_repeated_measures(data_dict)
        if omni is not None:
            chi2, p_omni, k, n = omni
            print(f"    OMNIBUS [{panel_label}]: Friedman chi2={chi2:.4f}, p={p_omni:.4e}, k={k} groups, n={n} per group")
        else:
            print(f"    OMNIBUS [{panel_label}]: Friedman test skipped (insufficient data)")

    # Add statistical comparisons with COMPACT, CONSISTENT style
    if comparisons:
        y_max = max([max(v) for v in data_list if len(v) > 0])
        y_min = min([min(v) for v in data_list if len(v) > 0])
        y_range = y_max - y_min
        
        bracket_y_start = y_max + 0.06 * y_range
        bracket_spacing = 0.07 * y_range
        line_height = 0.015 * y_range
        text_offset = 0.005 * y_range
        
        raw_results = []
        for idx, (label1, label2, paired) in enumerate(comparisons):
            if label1 in label_to_pos and label2 in label_to_pos:
                data1 = data_dict.get(label1, [])
                data2 = data_dict.get(label2, [])
                p_value, test_name, stat, n1, n2 = perform_statistical_test(data1, data2, paired=paired)
                raw_results.append((idx, label1, label2, paired, p_value, test_name, stat, n1, n2))
            else:
                raw_results.append((idx, label1, label2, paired, None, 'skipped', float('nan'), 0, 0))
        
        raw_pvals = [r[4] for r in raw_results if r[4] is not None]
        if correction == 'holm' and len(raw_pvals) > 1:
            corrected_pvals = holm_bonferroni_correction(raw_pvals)
        else:
            corrected_pvals = list(raw_pvals)
        
        corr_idx = 0
        print(f"    STATS [{panel_label}] — {len(raw_results)} comparisons, correction={correction}:")
        for r in raw_results:
            idx_c, l1, l2, paired, p_raw, tname, stat, n1, n2 = r
            if p_raw is not None:
                p_corr = corrected_pvals[corr_idx]
                corr_idx += 1
                paired_str = "paired" if paired else "unpaired"
                print(f"      {l1} vs {l2}: {tname} ({paired_str}), stat={stat:.4f}, "
                      f"n1={n1}, n2={n2}, p_raw={p_raw:.4e}, p_corrected={p_corr:.4e}")
                p_use = p_corr
            else:
                print(f"      {l1} vs {l2}: skipped (label not in data)")
                p_use = 1.0

            if l1 in label_to_pos and l2 in label_to_pos:
                x1 = label_to_pos[l1]
                x2 = label_to_pos[l2]
                bracket_y = bracket_y_start + idx_c * bracket_spacing
                add_significance_bracket(ax, x1, x2, bracket_y, p_use,
                                        line_height=line_height, text_offset=text_offset)
        
        current_ylim = ax.get_ylim()
        top_padding = bracket_y_start + len(raw_results) * bracket_spacing + 0.05 * y_range
        ax.set_ylim(current_ylim[0], top_padding)
    
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FONT_SIZE_TICK, rotation=45, ha='right')
    ax.set_ylabel(ylabel, fontsize=FONT_SIZE_LABEL)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICK, width=TICK_WIDTH, length=TICK_LENGTH)
    
    if title:
        ax.set_title(title, fontsize=FONT_SIZE_TITLE, fontweight='bold')


def extract_band_power(spectral_results, period, signal_type, freq_band, debug=False):
    """
    Extract power in a specific frequency band from PSD.
    
    Args:
        spectral_results: List of spectral result dicts
        period: 'pre_stim', 'transient', 'sustained', 'post_stim'
        signal_type: 'lfp' or 'fiber'
        freq_band: tuple (low_freq, high_freq) in Hz
        debug: Print debug info
    
    Returns:
        List of mean band power values (one per trial), in dB
    """
    period_keys = {
        'pre_stim': ['pre_stim', 'prestim', 'pre'],
        'transient': ['transient', 'stim_transient', 'trans'],
        'sustained': ['sustained', 'stim_sustained', 'sust'],
        'post_stim': ['post_stim', 'poststim', 'post']
    }
    
    values = []
    for i, result in enumerate(spectral_results):
        if result is None:
            values.append(None)
            continue
        try:
            # Get frequency axis
            freq = None
            if 'freq' in result:
                freq = np.array(result['freq']).flatten()
            elif 'psd_lfp' in result and isinstance(result['psd_lfp'], dict):
                # Try to find freq in coherence or psd sub-dict
                for sub_key in ['overall', 'pre_stim']:
                    if sub_key in result.get('coherence', {}):
                        # freq should be at top level
                        break
            
            if freq is None:
                if debug and i == 0:
                    print(f"      DEBUG: No frequency axis found in result")
                values.append(None)
                continue
            
            # Get PSD data
            psd_key = f'psd_{signal_type}'
            psd_dict = result.get(psd_key, {})
            
            if not isinstance(psd_dict, dict):
                values.append(None)
                continue
            
            # Find period data
            period_data = None
            for key in period_keys.get(period, [period]):
                if key in psd_dict:
                    period_data = np.array(psd_dict[key]).flatten()
                    break
            
            if period_data is None or period_data.size == 0:
                values.append(None)
                continue
            
            # Extract band power
            low_freq, high_freq = freq_band
            band_mask = (freq >= low_freq) & (freq <= high_freq)
            
            if np.sum(band_mask) > 0 and len(period_data) == len(freq):
                band_power = period_data[band_mask]
                finite_mask = np.isfinite(band_power)
                if np.sum(finite_mask) > 0:
                    # Mean power in band (already in dB)
                    mean_band_power = np.nanmean(band_power[finite_mask])
                    values.append(mean_band_power)
                    if debug and i == 0:
                        print(f"      DEBUG: {signal_type} band power [{low_freq}-{high_freq}Hz] for {period}: {mean_band_power:.2f} dB")
                else:
                    values.append(None)
            else:
                if debug and i == 0:
                    print(f"      DEBUG: Band mask empty or size mismatch. freq len={len(freq)}, psd len={len(period_data)}, band_mask sum={np.sum(band_mask)}")
                values.append(None)
        except Exception as e:
            if debug:
                print(f"      DEBUG: Error extracting band power: {e}")
            values.append(None)
    
    return values


def extract_spectral_metric(spectral_results, period, metric_type, debug=False):
    """
    Extract a metric from spectral results.
    
    Period names: 'pre_stim', 'transient', 'sustained', 'post_stim'
    Metric types: 'fiber_mean', 'coherence_mean', 'psd_lfp_mean', 'psd_fiber_mean'
    """
    # Map period names to potential MATLAB key names
    period_keys = {
        'pre_stim': ['pre_stim', 'prestim', 'pre'],
        'transient': ['transient', 'stim_transient', 'trans'],
        'sustained': ['sustained', 'stim_sustained', 'sust'],
        'post_stim': ['post_stim', 'poststim', 'post']
    }
    
    values = []
    for i, result in enumerate(spectral_results):
        if result is None:
            values.append(None)
            continue
        try:
            if metric_type == 'fiber_mean':
                # Extract from raw traces by computing period mean
                traces = result.get('traces', {})
                if isinstance(traces, dict) and 'fiber' in traces and 'time_sec' in traces:
                    fiber = np.array(traces.get('fiber', [])).flatten()
                    t = np.array(traces.get('time_sec', [])).flatten()
                    if period == 'pre_stim':
                        mask = (t >= -PRE_STIM_DURATION) & (t < 0)
                    elif period == 'transient':
                        mask = (t >= 0) & (t < STIM_TRANSIENT_END)
                    elif period == 'sustained':
                        mask = (t >= STIM_SUSTAINED_START) & (t < STIM_DURATION)
                    elif period == 'post_stim':
                        mask = (t >= STIM_DURATION) & (t < STIM_DURATION + POST_STIM_DURATION)
                    else:
                        mask = np.ones(len(t), dtype=bool)
                    if np.sum(mask) > 0 and len(fiber) == len(t):
                        values.append(np.nanmean(fiber[mask]))
                    else:
                        values.append(None)
                else:
                    values.append(None)
            
            elif metric_type in ['coherence_mean', 'psd_lfp_mean', 'psd_fiber_mean']:
                # Try multiple key names for the metric
                if metric_type == 'coherence_mean':
                    data_dict = result.get('coherence', {})
                elif metric_type == 'psd_lfp_mean':
                    data_dict = result.get('psd_lfp', {})
                else:  # psd_fiber_mean
                    data_dict = result.get('psd_fiber', {})
                
                # If not a dict, try to access as nested structure
                if not isinstance(data_dict, dict):
                    if debug and i == 0:
                        print(f"      DEBUG: {metric_type} is not a dict, type={type(data_dict)}")
                    values.append(None)
                    continue
                
                # Try various period key names
                found = False
                for key in period_keys.get(period, [period]):
                    if key in data_dict:
                        period_data = np.array(data_dict[key]).flatten()
                        if period_data.size > 0:
                            finite_mask = np.isfinite(period_data)
                            n_finite = np.sum(finite_mask)
                            if n_finite > 0:
                                values.append(np.nanmean(period_data[finite_mask]))
                                found = True
                                if debug and i == 0:
                                    print(f"      DEBUG: {metric_type}/{period} found (key='{key}'), size={period_data.size}, finite={n_finite}, mean={np.nanmean(period_data[finite_mask]):.3f}")
                                break
                            else:
                                if debug and i == 0:
                                    print(f"      DEBUG: {metric_type}/{period} found (key='{key}') but all NaN, size={period_data.size}")
                        else:
                            if debug and i == 0:
                                print(f"      DEBUG: {metric_type}/{period} found (key='{key}') but empty array")
                
                if not found:
                    # Debug: show available keys for first trial
                    if debug and i == 0:
                        avail = list(data_dict.keys()) if isinstance(data_dict, dict) else []
                        print(f"      DEBUG: {metric_type}/{period} not found. Available keys: {avail}")
                    values.append(None)
            else:
                values.append(None)
        except Exception as e:
            if debug:
                print(f"      DEBUG: Error extracting {metric_type}/{period}: {e}")
            values.append(None)
    return values


# =============================================================================
# FIGURE 1 & 2: TRACES
# =============================================================================

def generate_stim_pulses(t, stim_freq_hz, stim_duration_sec):
    """
    Generate square biphasic stimulation pulses.
    
    Parameters:
    -----------
    t : array
        Time vector (relative to stim onset, so 0 = stim start)
    stim_freq_hz : float
        Stimulation frequency in Hz (e.g., 40 or 135)
    stim_duration_sec : float
        Duration of stimulation in seconds
    
    Returns:
    --------
    pulses : array
        Pulse signal (1 = positive phase, -1 = negative phase, 0 = no pulse)
    """
    pulses = np.zeros_like(t)
    
    # Only generate pulses during stimulation period (0 to stim_duration_sec)
    stim_mask = (t >= 0) & (t <= stim_duration_sec)
    if not np.any(stim_mask):
        return pulses
    
    # Period of one complete biphasic pulse (positive + negative phase)
    period = 1.0 / stim_freq_hz  # seconds
    half_period = period / 2.0  # Each phase is half the period
    
    # Get time values during stimulation
    t_stim = t[stim_mask]
    
    # Generate square biphasic pulses
    # For each time point, determine which phase of the pulse cycle we're in
    cycle_times = t_stim % period
    
    # Positive phase: 0 to half_period
    # Negative phase: half_period to period
    positive_phase = cycle_times < half_period
    negative_phase = ~positive_phase
    
    # Set pulse values
    stim_indices = np.where(stim_mask)[0]
    pulses[stim_indices[positive_phase]] = 1.0
    pulses[stim_indices[negative_phase]] = -1.0
    
    return pulses


def plot_trace_panel(ax, t, data, color, ylabel, scalebar_val, scalebar_label, show_stim=False, show_xlabel=False, is_motion=False):
    """Plot a single trace panel with ONLY scale bars (no axis lines)."""
    if is_motion:
        data_smooth = uniform_filter1d(data, size=5)
        linewidth = LINE_WIDTH_TRACE * 1.5
    else:
        data_smooth = data
        linewidth = LINE_WIDTH_TRACE

    ax.plot(t, data_smooth, color=color, linewidth=linewidth)

    data_range = np.nanmax(data_smooth) - np.nanmin(data_smooth)
    margin = max(data_range * 0.1, 0.1)
    ax.set_ylim(np.nanmin(data_smooth) - margin, np.nanmax(data_smooth) + margin)
    ax.set_xlim(t.min(), t.max())

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis='both', which='both', length=0, labelleft=False, labelbottom=False)
    ax.set_xticks([])
    ax.set_yticks([])

    vbar_proportional(ax, 0.0, scalebar_val, scalebar_label, "left", color)


def _avg_grid_num_points(raw_trials, stim_onset_time):
    """Number of time samples after pre-stim skip (for common grid resolution)."""
    for trial_data in raw_trials:
        if trial_data is not None and len(trial_data["t"]) > 10:
            t_trial_raw = trial_data["t"] - stim_onset_time
            trial_keep = t_trial_raw >= -PRE_STIM_EFFECTIVE
            return int(max(np.sum(trial_keep), 2))
    return 1000


def compute_trial_averaged_lfp_gevi(raw_trials, stim_onset_time, n_points):
    """
    Mean ± SEM of LFP (µV) and GEVI (ΔF/F %) on a uniform time grid, with the same
    smoothing as create_traces_figure (5-sample moving average on mean and SEM).

    Returns
    -------
    tuple or None
        (t_common, lfp_avg_s, lfp_sem_s, gevi_avg_s, gevi_sem_s, n_trials) or None if no data.
    """
    t_min = -PRE_STIM_EFFECTIVE
    t_max = STIM_DURATION + POST_STIM_DURATION
    t_common = np.linspace(t_min, t_max, max(int(n_points), 2))

    gevi_all = []
    lfp_all = []

    for trial_data in raw_trials:
        if trial_data is not None and len(trial_data["t"]) > 10:
            t_trial_raw = trial_data["t"] - stim_onset_time
            trial_keep = t_trial_raw >= -PRE_STIM_EFFECTIVE
            t_trial = t_trial_raw[trial_keep]
            gevi_trial = (trial_data["gevi"] * 100)[trial_keep]
            lfp_trial = (trial_data["lfp"] if not LFP_IN_VOLTS else trial_data["lfp"] * 1e6)[trial_keep]

            if len(t_trial) > 10:
                gevi_interp = np.interp(t_common, t_trial, gevi_trial)
                lfp_interp = np.interp(t_common, t_trial, lfp_trial)
                gevi_all.append(gevi_interp)
                lfp_all.append(lfp_interp)

    if len(gevi_all) == 0:
        return None

    gevi_all = np.array(gevi_all)
    lfp_all = np.array(lfp_all)

    gevi_avg = np.nanmean(gevi_all, axis=0)
    gevi_sem = np.nanstd(gevi_all, axis=0) / np.sqrt(gevi_all.shape[0])
    lfp_avg = np.nanmean(lfp_all, axis=0)
    lfp_sem = np.nanstd(lfp_all, axis=0) / np.sqrt(lfp_all.shape[0])

    smooth_window = 5
    gevi_avg_s = uniform_filter1d(gevi_avg, size=smooth_window)
    lfp_avg_s = uniform_filter1d(lfp_avg, size=smooth_window)
    gevi_sem_s = uniform_filter1d(gevi_sem, size=smooth_window)
    lfp_sem_s = uniform_filter1d(lfp_sem, size=smooth_window)

    return t_common, lfp_avg_s, lfp_sem_s, gevi_avg_s, gevi_sem_s, int(gevi_all.shape[0])


def create_overlaid_avg_traces_figure(raw_trials_r1, raw_trials_r2, stim_onset_time, label_r1, label_r2):
    """
    Trial-averaged LFP and Population Vm (GEVI) for both conditions on shared axes,
    matching styling of the trial-averaged panels in create_traces_figure (SEM shading,
    stim onset/offset lines, spine/tick styling). Uses COLOR_40HZ / COLOR_135HZ for R1/R2.
    """
    n_points = max(
        _avg_grid_num_points(raw_trials_r1, stim_onset_time),
        _avg_grid_num_points(raw_trials_r2, stim_onset_time),
    )

    res1 = compute_trial_averaged_lfp_gevi(raw_trials_r1, stim_onset_time, n_points)
    res2 = compute_trial_averaged_lfp_gevi(raw_trials_r2, stim_onset_time, n_points)

    fig, (ax_lfp, ax_gevi) = plt.subplots(
        2, 1, figsize=(16, 10), sharex=True, gridspec_kw={"hspace": 0.12}
    )

    if res1 is None or res2 is None:
        fig.suptitle(
            f"{MOUSE_ID} — Trial-averaged overlay: insufficient data "
            f"({label_r1 if res1 is None else 'ok'}, {label_r2 if res2 is None else 'ok'})",
            fontsize=FONT_SIZE_SUPTITLE,
            fontweight="bold",
            y=0.98,
        )
        fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.08)
        return fig

    t_c, lfp1, lfp1_sem, g1, g1_sem, n1 = res1
    _, lfp2, lfp2_sem, g2, g2_sem, n2 = res2

    t_min = -PRE_STIM_EFFECTIVE
    t_max = STIM_DURATION + POST_STIM_DURATION

    def style_avg_axis(ax, show_x_labels=True):
        ax.axvline(0, color="black", linestyle="--", linewidth=1.5)
        ax.axvline(STIM_DURATION, color="black", linestyle=":", linewidth=1.5)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
        ax.tick_params(
            axis="both",
            which="major",
            labelsize=FONT_SIZE_TICK,
            width=TICK_WIDTH,
            length=TICK_LENGTH,
        )
        if not show_x_labels:
            ax.tick_params(labelbottom=False)

    # LFP — R1 then R2 (SEM under lines for readability)
    ax_lfp.fill_between(
        t_c, lfp1 - lfp1_sem, lfp1 + lfp1_sem, color=COLOR_40HZ, alpha=0.3, linewidth=0
    )
    ax_lfp.plot(
        t_c,
        lfp1,
        color=COLOR_40HZ,
        linewidth=LINE_WIDTH_THICK,
        label=f"{label_r1} (n={n1})",
    )
    ax_lfp.fill_between(
        t_c, lfp2 - lfp2_sem, lfp2 + lfp2_sem, color=COLOR_135HZ, alpha=0.3, linewidth=0
    )
    ax_lfp.plot(
        t_c,
        lfp2,
        color=COLOR_135HZ,
        linewidth=LINE_WIDTH_THICK,
        label=f"{label_r2} (n={n2})",
    )
    style_avg_axis(ax_lfp, show_x_labels=False)
    ax_lfp.set_ylabel("LFP (µV)", fontsize=FONT_SIZE_LABEL)
    ax_lfp.set_title(
        f"Trial-Averaged LFP ± SEM (overlaid, n={n1} vs {n2} trials)",
        fontsize=FONT_SIZE_TITLE,
        fontweight="bold",
        pad=10,
    )
    ax_lfp.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, framealpha=0.9)
    ax_lfp.set_xlim(t_min, t_max)

    # GEVI / Population Vm
    ax_gevi.fill_between(
        t_c, g1 - g1_sem, g1 + g1_sem, color=COLOR_40HZ, alpha=0.3, linewidth=0
    )
    ax_gevi.plot(
        t_c,
        g1,
        color=COLOR_40HZ,
        linewidth=LINE_WIDTH_THICK,
        label=f"{label_r1} (n={n1})",
    )
    ax_gevi.fill_between(
        t_c, g2 - g2_sem, g2 + g2_sem, color=COLOR_135HZ, alpha=0.3, linewidth=0
    )
    ax_gevi.plot(
        t_c,
        g2,
        color=COLOR_135HZ,
        linewidth=LINE_WIDTH_THICK,
        label=f"{label_r2} (n={n2})",
    )
    style_avg_axis(ax_gevi, show_x_labels=True)
    ax_gevi.set_ylabel("Population Vm (ΔF/F %)", fontsize=FONT_SIZE_LABEL)
    ax_gevi.set_xlabel("Time from stim onset (s)", fontsize=FONT_SIZE_LABEL)
    ax_gevi.legend(loc="upper right", fontsize=FONT_SIZE_LEGEND, framealpha=0.9)
    ax_gevi.set_xlim(t_min, t_max)

    # Shared y-limits so both conditions use the same scale per modality
    def pair_ylim(lower_envelope, upper_envelope):
        lo = np.nanmin(lower_envelope)
        hi = np.nanmax(upper_envelope)
        span = hi - lo
        m = max(span * 0.05, 0.01) if np.isfinite(span) else 1.0
        return lo - m, hi + m

    ax_lfp.set_ylim(
        pair_ylim(
            np.minimum(lfp1 - lfp1_sem, lfp2 - lfp2_sem),
            np.maximum(lfp1 + lfp1_sem, lfp2 + lfp2_sem),
        )
    )
    ax_gevi.set_ylim(
        pair_ylim(
            np.minimum(g1 - g1_sem, g2 - g2_sem),
            np.maximum(g1 + g1_sem, g2 + g2_sem),
        )
    )

    fig.suptitle(
        f"{MOUSE_ID} — Trial-averaged traces overlaid ({label_r1} vs {label_r2})",
        fontsize=FONT_SIZE_SUPTITLE,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(left=0.12, right=0.95, top=0.92, bottom=0.07)

    return fig


def create_traces_figure(session_id, condition_name, raw_trials, stim_onset_time):
    """
    Create figure with stacked LFP, Fiber, Motion traces (NOT overlaid).

    Layout: Stim pulses, single-trial full view, zoom (−0.5 s to stim+0.5 s), trial-averaged.
    Each single/zoom section has Stim, LFP, Fiber, Motion. Trial-averaged: LFP and Fiber ± SEM.

    Applies PRE_STIM_SKIP_SEC to exclude noisy initial data if configured.
    """
    # Extract stimulation frequency from condition_name
    if '40' in condition_name or '40Hz' in condition_name:
        stim_freq = 40
    elif '135' in condition_name or '135Hz' in condition_name:
        stim_freq = 135
    else:
        stim_freq = 40  # Default
    
    # 11 rows: full (4) + zoom title + zoom traces (4) + averaged (2)
    fig = plt.figure(figsize=(16, 26))
    gs = GridSpec(11, 1, figure=fig, height_ratios=[0.4, 1, 1, 0.8, 0.25, 0.4, 1, 1, 0.8, 1.2, 1.2], hspace=0.12)
    
    # Get representative trial
    trial = None
    for t in raw_trials:
        if t is not None:
            trial = t
            break
    
    if trial is None:
        fig.suptitle(f"{MOUSE_ID} - {condition_name} - No data available", fontsize=FONT_SIZE_SUPTITLE)
        return fig
    
    # Get time relative to stim onset
    t_raw = trial['t'] - stim_onset_time
    
    # Apply pre-stim skip: exclude data before -PRE_STIM_EFFECTIVE (i.e., skip first PRE_STIM_SKIP_SEC)
    keep_mask = t_raw >= -PRE_STIM_EFFECTIVE
    t_full = t_raw[keep_mask]
    gevi_percent = (trial['gevi'] * 100)[keep_mask]
    lfp_uv = (trial['lfp'] if not LFP_IN_VOLTS else trial['lfp'] * 1e6)[keep_mask]
    motion_cm_s = trial['motion'][keep_mask]
    
    # === SECTION 0: Stimulation Pulses (row 0) ===
    ax_stim = fig.add_subplot(gs[0])
    stim_pulses = generate_stim_pulses(t_full, stim_freq, STIM_DURATION)
    ax_stim.plot(t_full, stim_pulses, color=COLOR_STIM_PULSE, linewidth=LINE_WIDTH_TRACE)
    ax_stim.set_ylim(-1.2, 1.2)
    full_xlim = (t_full.min(), t_full.max())  # Store for all full-view axes
    ax_stim.set_xlim(full_xlim)
    # Clean axis - no lines, no ticks
    for spine in ax_stim.spines.values():
        spine.set_visible(False)
    ax_stim.tick_params(axis='both', which='both', length=0, labelleft=False, labelbottom=False)
    ax_stim.set_xticks([])
    ax_stim.set_yticks([])
    
    # === SECTION 1: Single Trial - Full view (rows 1-3) ===
    # These share x-axis with each other for alignment, but NOT with zoomed section
    ax_lfp_full = fig.add_subplot(gs[1], sharex=ax_stim)
    ax_gevi_full = fig.add_subplot(gs[2], sharex=ax_stim)
    ax_motion_full = fig.add_subplot(gs[3], sharex=ax_stim)

    # Plot single trial traces
    plot_trace_panel(ax_lfp_full, t_full, lfp_uv, COLOR_LFP, 'LFP (µV)', 
                     SCALEBAR_LFP_UV, f"{int(SCALEBAR_LFP_UV)} µV")
    plot_trace_panel(ax_gevi_full, t_full, gevi_percent, COLOR_GEVI, 'Fiber (ΔF/F %)', 
                     SCALEBAR_GEVI_PERCENT, f"{SCALEBAR_GEVI_PERCENT:.1f}%")
    plot_trace_panel(ax_motion_full, t_full, motion_cm_s, COLOR_MOTION, 'Speed (cm/s)', 
                     SCALEBAR_MOTION_CM_S, f"{int(SCALEBAR_MOTION_CM_S)} cm/s", show_xlabel=False, is_motion=True)
    
    ax_stim.set_title(f"Stimulation Pulses ({stim_freq}Hz) - Single Trial Traces (Trial {REPRESENTATIVE_TRIAL})", 
                      fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    
    # Add horizontal time scale bar to bottom of single trial section
    add_time_scalebar(ax_motion_full, 1.0, "1 s")
    
    # === SECTION 2: Zoomed (-0.5s to stim+0.5s) (rows 4-7) ===
    # NOTE: Zoomed section has INDEPENDENT x-axis from full section!
    zoom_start = -0.5
    zoom_end = STIM_DURATION + 0.5
    zoom_mask = (t_full >= zoom_start) & (t_full <= zoom_end)
    t_zoom = t_full[zoom_mask]
    
    if len(t_zoom) > 0:
        # Title for zoomed section (row 4)
        ax_zoom_title = fig.add_subplot(gs[4])
        ax_zoom_title.axis('off')
        ax_zoom_title.text(0.5, 0.5, f"Zoomed View ({zoom_start:.1f}s to {zoom_end:.1f}s)", 
                          ha='center', va='center', fontsize=FONT_SIZE_TITLE, fontweight='bold',
                          transform=ax_zoom_title.transAxes)
        
        # Stim pulses in zoomed view (row 5) - NO sharex with full plots!
        ax_stim_zoom = fig.add_subplot(gs[5])
        stim_pulses_zoom = stim_pulses[zoom_mask]
        ax_stim_zoom.plot(t_zoom, stim_pulses_zoom, color=COLOR_STIM_PULSE, linewidth=LINE_WIDTH_TRACE)
        ax_stim_zoom.set_ylim(-1.2, 1.2)
        ax_stim_zoom.set_xlim(zoom_start, zoom_end)  # Explicitly set zoom limits
        for spine in ax_stim_zoom.spines.values():
            spine.set_visible(False)
        ax_stim_zoom.tick_params(axis='both', which='both', length=0, labelleft=False, labelbottom=False)
        ax_stim_zoom.set_xticks([])
        ax_stim_zoom.set_yticks([])
        
        # Zoomed traces share x with each other (within zoom section), NOT with full section
        ax_lfp_zoom = fig.add_subplot(gs[6], sharex=ax_stim_zoom)
        ax_gevi_zoom = fig.add_subplot(gs[7], sharex=ax_stim_zoom)
        ax_motion_zoom = fig.add_subplot(gs[8], sharex=ax_stim_zoom)
        
        # Extract zoomed data
        lfp_zoom_data = lfp_uv[zoom_mask]
        gevi_zoom_data = gevi_percent[zoom_mask]
        motion_zoom_data = motion_cm_s[zoom_mask]
        
        plot_trace_panel(ax_lfp_zoom, t_zoom, lfp_zoom_data, COLOR_LFP, 'LFP (µV)', 
                         SCALEBAR_LFP_UV/2, f"{int(SCALEBAR_LFP_UV/2)} µV")
        plot_trace_panel(ax_gevi_zoom, t_zoom, gevi_zoom_data, COLOR_GEVI, 'Fiber (ΔF/F %)', 
                         SCALEBAR_GEVI_PERCENT/2, f"{SCALEBAR_GEVI_PERCENT/2:.2f}%")
        plot_trace_panel(ax_motion_zoom, t_zoom, motion_zoom_data, COLOR_MOTION, 'Speed (cm/s)', 
                         SCALEBAR_MOTION_CM_S/2, f"{int(SCALEBAR_MOTION_CM_S/2)} cm/s", show_xlabel=False, is_motion=True)
        
        add_time_scalebar(ax_motion_zoom, 0.2, "0.2 s")

    # === SECTION 3: Trial-Averaged with SEM (rows 9–10) ===
    # Create a common time vector for all trials (using the effective time range)
    t_min = -PRE_STIM_EFFECTIVE
    t_max = STIM_DURATION + POST_STIM_DURATION
    n_points = len(t_full)
    t_common = np.linspace(t_min, t_max, n_points)
    
    # Collect all trials for averaging and SEM
    gevi_all = []
    lfp_all = []
    
    for trial_data in raw_trials:
        if trial_data is not None and len(trial_data['t']) > 10:
            # Get trial time relative to stim onset and apply skip
            t_trial_raw = trial_data['t'] - stim_onset_time
            trial_keep = t_trial_raw >= -PRE_STIM_EFFECTIVE
            t_trial = t_trial_raw[trial_keep]
            gevi_trial = (trial_data['gevi'] * 100)[trial_keep]
            lfp_trial = (trial_data['lfp'] if not LFP_IN_VOLTS else trial_data['lfp'] * 1e6)[trial_keep]
            
            # Interpolate to common time grid
            if len(t_trial) > 10:
                gevi_interp = np.interp(t_common, t_trial, gevi_trial)
                lfp_interp = np.interp(t_common, t_trial, lfp_trial)
                gevi_all.append(gevi_interp)
                lfp_all.append(lfp_interp)
    
    if len(gevi_all) > 0:
        gevi_all = np.array(gevi_all)
        lfp_all = np.array(lfp_all)
        
        # Compute mean and SEM across trials (axis=0 is trials, axis=1 is time)
        gevi_avg = np.nanmean(gevi_all, axis=0)
        gevi_sem = np.nanstd(gevi_all, axis=0) / np.sqrt(gevi_all.shape[0])
        lfp_avg = np.nanmean(lfp_all, axis=0)
        lfp_sem = np.nanstd(lfp_all, axis=0) / np.sqrt(lfp_all.shape[0])
        
        # Apply light smoothing to averaged traces for cleaner appearance
        from scipy.ndimage import uniform_filter1d
        smooth_window = 5  # 5-sample moving average
        gevi_avg_smooth = uniform_filter1d(gevi_avg, size=smooth_window)
        lfp_avg_smooth = uniform_filter1d(lfp_avg, size=smooth_window)
        gevi_sem_smooth = uniform_filter1d(gevi_sem, size=smooth_window)
        lfp_sem_smooth = uniform_filter1d(lfp_sem, size=smooth_window)
        
        ax_lfp_avg = fig.add_subplot(gs[9])
        ax_gevi_avg = fig.add_subplot(gs[10], sharex=ax_lfp_avg)
        
        # LFP averaged with SEM - NO red shading, just stim onset line
        ax_lfp_avg.fill_between(t_common, lfp_avg_smooth - lfp_sem_smooth, lfp_avg_smooth + lfp_sem_smooth, 
                                 color=COLOR_LFP, alpha=0.3, linewidth=0)
        ax_lfp_avg.plot(t_common, lfp_avg_smooth, color=COLOR_LFP, linewidth=LINE_WIDTH_THICK)
        ax_lfp_avg.axvline(0, color='black', linestyle='--', linewidth=1.5, label='Stim onset')
        ax_lfp_avg.axvline(STIM_DURATION, color='black', linestyle=':', linewidth=1.5, label='Stim offset')
        ax_lfp_avg.spines['top'].set_visible(False)
        ax_lfp_avg.spines['right'].set_visible(False)
        ax_lfp_avg.spines['left'].set_linewidth(AXIS_LINEWIDTH)
        ax_lfp_avg.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
        ax_lfp_avg.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICK, width=TICK_WIDTH, length=TICK_LENGTH)
        ax_lfp_avg.set_ylabel('LFP (µV)', fontsize=FONT_SIZE_LABEL)
        ax_lfp_avg.tick_params(labelbottom=False)
        ax_lfp_avg.set_title(f'Trial-Averaged Traces (n={gevi_all.shape[0]} trials) ± SEM', 
                             fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
        
        # Fiber averaged with SEM - NO red shading
        ax_gevi_avg.fill_between(t_common, gevi_avg_smooth - gevi_sem_smooth, gevi_avg_smooth + gevi_sem_smooth, 
                                  color=COLOR_GEVI, alpha=0.3, linewidth=0)
        ax_gevi_avg.plot(t_common, gevi_avg_smooth, color=COLOR_GEVI, linewidth=LINE_WIDTH_THICK)
        ax_gevi_avg.axvline(0, color='black', linestyle='--', linewidth=1.5)
        ax_gevi_avg.axvline(STIM_DURATION, color='black', linestyle=':', linewidth=1.5)
        ax_gevi_avg.spines['top'].set_visible(False)
        ax_gevi_avg.spines['right'].set_visible(False)
        ax_gevi_avg.spines['left'].set_linewidth(AXIS_LINEWIDTH)
        ax_gevi_avg.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
        ax_gevi_avg.tick_params(axis='both', which='major', labelsize=FONT_SIZE_TICK, width=TICK_WIDTH, length=TICK_LENGTH)
        ax_gevi_avg.set_ylabel('Population Vm (ΔF/F %)', fontsize=FONT_SIZE_LABEL)
        ax_gevi_avg.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL)
    
    fig.suptitle(f'{MOUSE_ID} - {condition_name} Traces', fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.995)
    
    # Adjust layout
    fig.subplots_adjust(left=0.12, right=0.95, top=0.95, bottom=0.04)
    
    return fig


def add_time_scalebar(ax, duration_sec, label):
    """Add horizontal time scale bar below the axis, positioned at bottom left."""
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    y_range = ylim[1] - ylim[0]
    
    # Position at bottom left of the plot area
    x_start = xlim[0]  # Start at left edge of data
    x_end = x_start + duration_sec
    y_pos = ylim[0] - y_range * 0.15  # Below the axis
    
    # Draw horizontal bar
    ax.plot([x_start, x_end], [y_pos, y_pos], color='black', linewidth=SCALEBAR_LINEWIDTH, clip_on=False)
    # Label centered below the bar
    ax.text((x_start + x_end) / 2, y_pos - y_range * 0.06, label, 
            ha='center', va='top', fontsize=FONT_SIZE_SCALEBAR, fontweight='bold', clip_on=False)


def create_traces_stim_onset_zoom_figure(session_id, condition_name, raw_trials, stim_onset_time):
    """
    Separate figure only: stimulation pulses, LFP, GEVI, and motion on a fixed window
    from STIM_ONSET_ZOOM_PRE_S before stim onset through STIM_ONSET_ZOOM_POST_S after
    (default −0.5 s to +0.5 s). Same representative trial and visual style as the main
    traces zoom panels (halved vertical scale bars).
    """
    if '40' in condition_name or '40Hz' in condition_name:
        stim_freq = 40
    elif '135' in condition_name or '135Hz' in condition_name:
        stim_freq = 135
    else:
        stim_freq = 40

    fig = plt.figure(figsize=(16, 11))
    gs = GridSpec(5, 1, figure=fig, height_ratios=[0.3, 0.4, 1, 1, 0.8], hspace=0.12)

    trial = None
    for t in raw_trials:
        if t is not None:
            trial = t
            break

    if trial is None:
        fig.suptitle(f"{MOUSE_ID} - {condition_name} - No data available", fontsize=FONT_SIZE_SUPTITLE)
        return fig

    t_raw = trial['t'] - stim_onset_time
    keep_mask = t_raw >= -PRE_STIM_EFFECTIVE
    t_full = t_raw[keep_mask]
    gevi_percent = (trial['gevi'] * 100)[keep_mask]
    lfp_uv = (trial['lfp'] if not LFP_IN_VOLTS else trial['lfp'] * 1e6)[keep_mask]
    motion_cm_s = trial['motion'][keep_mask]

    z0 = -float(STIM_ONSET_ZOOM_PRE_S)
    z1 = float(STIM_ONSET_ZOOM_POST_S)
    zoom_mask = (t_full >= z0) & (t_full <= z1)
    t_zoom = t_full[zoom_mask]
    stim_pulses = generate_stim_pulses(t_full, stim_freq, STIM_DURATION)

    if len(t_zoom) == 0:
        fig.suptitle(
            f"{MOUSE_ID} - {condition_name} - No samples in [{z0:.2f}, {z1:.2f}] s",
            fontsize=FONT_SIZE_SUPTITLE,
        )
        fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.08)
        return fig

    ax_title = fig.add_subplot(gs[0])
    ax_title.axis("off")
    ax_title.text(
        0.5,
        0.5,
        f"Stim onset zoom ({z0:.1f} s to {z1:.1f} s)",
        ha="center",
        va="center",
        fontsize=FONT_SIZE_TITLE,
        fontweight="bold",
        transform=ax_title.transAxes,
    )

    ax_stim = fig.add_subplot(gs[1])
    ax_stim.plot(t_zoom, stim_pulses[zoom_mask], color=COLOR_STIM_PULSE, linewidth=LINE_WIDTH_TRACE)
    ax_stim.set_ylim(-1.2, 1.2)
    ax_stim.set_xlim(z0, z1)
    for spine in ax_stim.spines.values():
        spine.set_visible(False)
    ax_stim.tick_params(axis="both", which="both", length=0, labelleft=False, labelbottom=False)
    ax_stim.set_xticks([])
    ax_stim.set_yticks([])

    ax_lfp = fig.add_subplot(gs[2], sharex=ax_stim)
    ax_gevi = fig.add_subplot(gs[3], sharex=ax_stim)
    ax_motion = fig.add_subplot(gs[4], sharex=ax_stim)

    plot_trace_panel(
        ax_lfp, t_zoom, lfp_uv[zoom_mask], COLOR_LFP, "LFP (µV)",
        SCALEBAR_LFP_UV / 2, f"{int(SCALEBAR_LFP_UV / 2)} µV",
    )
    plot_trace_panel(
        ax_gevi, t_zoom, gevi_percent[zoom_mask], COLOR_GEVI, "Fiber (ΔF/F %)",
        SCALEBAR_GEVI_PERCENT / 2, f"{SCALEBAR_GEVI_PERCENT / 2:.2f}%",
    )
    plot_trace_panel(
        ax_motion, t_zoom, motion_cm_s[zoom_mask], COLOR_MOTION, "Speed (cm/s)",
        SCALEBAR_MOTION_CM_S / 2, f"{int(SCALEBAR_MOTION_CM_S / 2)} cm/s",
        show_xlabel=False, is_motion=True,
    )

    for ax_z in (ax_stim, ax_lfp, ax_gevi, ax_motion):
        ax_z.set_xlim(z0, z1)

    add_time_scalebar(ax_motion, 0.2, "0.2 s")

    fig.suptitle(
        f"{MOUSE_ID} - {condition_name} — Stimulation onset ({stim_freq} Hz, Trial {REPRESENTATIVE_TRIAL})",
        fontsize=FONT_SIZE_SUPTITLE,
        fontweight="bold",
        y=0.98,
    )
    fig.subplots_adjust(left=0.12, right=0.95, top=0.90, bottom=0.06)

    return fig


# =============================================================================
# FIGURE 3 & 4: SPECTRAL HEATMAPS
# =============================================================================

def _spectrogram_baseline_time_mask(time_arr):
    """
    Columns of the spectrogram used for per-frequency baseline mean (fractional-change norm).
    Primary window: [SPECTROGRAM_BASELINE_T_LO, SPECTROGRAM_BASELINE_T_HI) relative to stim onset.
    If empty (e.g. sparse time grid), fall back to all pre-onset bins (time < 0).
    """
    t = np.asarray(time_arr).squeeze()
    m = (t >= SPECTROGRAM_BASELINE_T_LO) & (t < SPECTROGRAM_BASELINE_T_HI)
    if not np.any(m):
        m = t < 0.0
    return m


def create_heatmaps_figure(session_id, condition_name, spectral_results, raw_trial, stim_onset_time):
    """
    Create figure with spectrograms and coherence spectrum.

    Layout: Speed (small), LFP spectrogram, Fiber spectrogram, Coherence SPECTRUM (not heatmap)
    - Spectrograms: fractional change (P-B)/(P+B) in linear power; B = per-frequency mean over bins in
      [SPECTROGRAM_BASELINE_T_LO, SPECTROGRAM_BASELINE_T_HI), else fallback to time < 0.
    - Colorbars with proper units (placed outside plot area for equal subplot widths)
    - Coherence as line plot with 3 periods (pre, stim, post) + SEM shading
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    
    # Apply pre-stim skip
    t_min_plot = -PRE_STIM_EFFECTIVE
    
    # Create figure with unequal heights: speed much smaller
    # Use GridSpec with TWO columns: main plot (wider) and colorbar (narrow)
    fig = plt.figure(figsize=(14, 14))
    gs = GridSpec(4, 2, figure=fig, height_ratios=[0.3, 1, 1, 1.2], 
                  width_ratios=[1, 0.03], hspace=0.25, wspace=0.05)
    
    # Get first valid spectral result for spectrogram display
    spectral_example = None
    for s in spectral_results:
        if s is not None:
            spectral_example = s
            break
    
    # Get time axis from raw trial (apply skip)
    if raw_trial is not None:
        t_raw = raw_trial['t'] - stim_onset_time
        keep_mask = t_raw >= t_min_plot
        t_full = t_raw[keep_mask]
        motion_cm_s = raw_trial['motion'][keep_mask]
    else:
        t_full = np.linspace(t_min_plot, STIM_DURATION + POST_STIM_DURATION, 1000)
        motion_cm_s = np.zeros_like(t_full)
    
    # === Speed heatmap (compact) ===
    ax0 = fig.add_subplot(gs[0, 0])  # Row 0, Column 0 (main plot)
    cax0 = fig.add_subplot(gs[0, 1])  # Row 0, Column 1 (colorbar)
    
    motion_2d = motion_cm_s.reshape(1, -1)
    vmax_motion = np.percentile(motion_cm_s[np.isfinite(motion_cm_s)], 99) if np.any(np.isfinite(motion_cm_s)) else 1.0
    im0 = ax0.imshow(motion_2d, aspect='auto', cmap=CMAP_SPEED,
                     extent=[t_full.min(), t_full.max(), 0, 1], origin='lower',
                     vmin=0, vmax=max(vmax_motion, 1.0))
    ax0.set_ylim([0, 1])
    ax0.set_yticks([])
    ax0.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK)
    ax0.set_ylabel('Speed\n(cm/s)', fontsize=FONT_SIZE_LABEL, rotation=0, ha='right', va='center')
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)
    ax0.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax0.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    # Add vertical lines for stim onset/offset (no red shading)
    ax0.axvline(0, color='black', linestyle='--', linewidth=1.5)
    ax0.axvline(STIM_DURATION, color='black', linestyle=':', linewidth=1.5)
    
    cbar0 = fig.colorbar(im0, cax=cax0, orientation='vertical')
    cbar0.set_label('cm/s', fontsize=FONT_SIZE_TICK)
    cbar0.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    
    ax0.set_title(f'{MOUSE_ID} - {condition_name} Spectral Analysis', 
                  fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    
    # === LFP spectrogram ===
    ax1 = fig.add_subplot(gs[1, 0])  # Row 1, Column 0 (main plot)
    cax1 = fig.add_subplot(gs[1, 1])  # Row 1, Column 1 (colorbar)
    if spectral_example and isinstance(spectral_example, dict) and 'spectrogram_lfp' in spectral_example:
        spec = spectral_example['spectrogram_lfp']
        if isinstance(spec, dict):
            power = np.array(spec.get('power', [])).squeeze()
            freq = np.array(spec.get('freq', [])).squeeze()
            time = np.array(spec.get('time', [])).squeeze() - stim_onset_time
            
            if power.size > 0 and power.ndim == 2:
                # Filter to effective time range
                time_mask = time >= t_min_plot
                time = time[time_mask]
                power = power[:, time_mask]
                
                # Baseline normalize using fractional change method (like MATLAB)
                # Formula: (P - P_baseline) / (P + P_baseline) for symmetric, bounded [-1, +1] range
                pre_mask = _spectrogram_baseline_time_mask(time)
                if np.any(pre_mask):
                    # Convert dB back to linear power for the formula
                    power_linear = 10**(power / 10.0)  # Convert dB to linear
                    baseline_linear = np.nanmean(power_linear[:, pre_mask], axis=1, keepdims=True)
                    
                    # Apply fractional change formula: (signal - baseline) / (signal + baseline)
                    numerator = power_linear - baseline_linear
                    denominator = power_linear + baseline_linear
                    # Avoid division by zero
                    denominator[denominator == 0] = 1e-10
                    power_norm = numerator / denominator  # Bounded [-1, +1]
                else:
                    # No baseline available, use zero-centered dB
                    power_norm = power - np.nanmean(power)
                
                # Use viridis colormap with NATURAL data scaling (no manipulation)
                # Let matplotlib determine the color range from the actual data
                im1 = ax1.pcolormesh(time, freq, power_norm, shading='auto', cmap='viridis')
                ax1.set_ylim(FREQ_RANGE)
                
                cbar1 = fig.colorbar(im1, cax=cax1, orientation='vertical')
                cbar1.set_label('Fractional Change\n(rel. baseline)', fontsize=FONT_SIZE_TICK)
                cbar1.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    
    ax1.axvline(0, color='white', linestyle='--', linewidth=1.5)
    ax1.axvline(STIM_DURATION, color='white', linestyle=':', linewidth=1.5)
    ax1.set_ylabel('LFP\nFrequency (Hz)', fontsize=FONT_SIZE_LABEL, rotation=0, ha='right', va='center')
    ax1.tick_params(labelbottom=False, labelsize=FONT_SIZE_TICK)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax1.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax1.set_xlim(t_min_plot, STIM_DURATION + POST_STIM_DURATION)
    
    # === Fiber spectrogram ===
    ax2 = fig.add_subplot(gs[2, 0], sharex=ax1)  # Row 2, Column 0 (main plot)
    cax2 = fig.add_subplot(gs[2, 1])  # Row 2, Column 1 (colorbar)
    if spectral_example and isinstance(spectral_example, dict) and 'spectrogram_fiber' in spectral_example:
        spec = spectral_example['spectrogram_fiber']
        if isinstance(spec, dict):
            power = np.array(spec.get('power', [])).squeeze()
            freq = np.array(spec.get('freq', [])).squeeze()
            time = np.array(spec.get('time', [])).squeeze() - stim_onset_time
            
            if power.size > 0 and power.ndim == 2:
                time_mask = time >= t_min_plot
                time = time[time_mask]
                power = power[:, time_mask]
                
                # Baseline normalize using fractional change method (like MATLAB)
                # Formula: (P - P_baseline) / (P + P_baseline) for symmetric, bounded [-1, +1] range
                pre_mask = _spectrogram_baseline_time_mask(time)
                if np.any(pre_mask):
                    # Convert dB back to linear power for the formula
                    power_linear = 10**(power / 10.0)  # Convert dB to linear
                    baseline_linear = np.nanmean(power_linear[:, pre_mask], axis=1, keepdims=True)
                    
                    # Apply fractional change formula: (signal - baseline) / (signal + baseline)
                    numerator = power_linear - baseline_linear
                    denominator = power_linear + baseline_linear
                    # Avoid division by zero
                    denominator[denominator == 0] = 1e-10
                    power_norm = numerator / denominator  # Bounded [-1, +1]
                else:
                    # No baseline available, use zero-centered dB
                    power_norm = power - np.nanmean(power)
                
                # Use viridis colormap with NATURAL data scaling (no manipulation)
                # Let matplotlib determine the color range from the actual data
                im2 = ax2.pcolormesh(time, freq, power_norm, shading='auto', cmap='viridis')
                ax2.set_ylim(FREQ_RANGE)
                
                cbar2 = fig.colorbar(im2, cax=cax2, orientation='vertical')
                cbar2.set_label('Fractional Change\n(rel. baseline)', fontsize=FONT_SIZE_TICK)
                cbar2.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    
    ax2.axvline(0, color='white', linestyle='--', linewidth=1.5)
    ax2.axvline(STIM_DURATION, color='white', linestyle=':', linewidth=1.5)
    ax2.set_ylabel('Fiber\nFrequency (Hz)', fontsize=FONT_SIZE_LABEL, rotation=0, ha='right', va='center')
    ax2.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL)
    ax2.tick_params(labelsize=FONT_SIZE_TICK)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax2.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    
    # === Coherence SPECTRUM (3 period lines with SEM) ===
    ax3 = fig.add_subplot(gs[3, :])  # Row 3, span both columns (no colorbar needed)
    
    # Collect period-based coherence across all trials
    coh_pre_all, coh_stim_all, coh_post_all = [], [], []
    freq_common = None
    
    # Debug: print structure of first result
    if spectral_results and spectral_results[0] is not None:
        first_result = spectral_results[0]
        print(f"    DEBUG coherence: Top-level keys = {list(first_result.keys()) if isinstance(first_result, dict) else type(first_result)}")
        if isinstance(first_result, dict) and 'coherence' in first_result:
            coh_top = first_result['coherence']
            print(f"    DEBUG coherence: coherence keys = {list(coh_top.keys()) if isinstance(coh_top, dict) else type(coh_top)}")
            # Print actual coherence data shapes
            if isinstance(coh_top, dict):
                for key in ['pre_stim', 'transient', 'sustained', 'post_stim']:
                    if key in coh_top:
                        coh_val = np.array(coh_top[key])
                        print(f"    DEBUG coherence: {key} shape = {coh_val.shape}, size = {coh_val.size}, any finite = {np.any(np.isfinite(coh_val))}")
        if isinstance(first_result, dict) and 'freq' in first_result:
            freq_test = np.array(first_result['freq']).flatten()
            print(f"    DEBUG coherence: freq found at top level, size = {freq_test.size}, range = [{freq_test.min() if freq_test.size > 0 else 'N/A'}, {freq_test.max() if freq_test.size > 0 else 'N/A'}]")
    
    for result in spectral_results:
        if result is None:
            continue
        
        coh_dict = result.get('coherence', {})
        if not isinstance(coh_dict, dict):
            continue
        
        # Get frequency axis - check top level first (MATLAB stores it there)
        if freq_common is None:
            # Try top-level freq first (most common in MATLAB output)
            freq_common = np.array(result.get('freq', [])).flatten()
            if freq_common.size == 0:
                # Try inside coherence dict
                freq_common = np.array(coh_dict.get('freq', [])).flatten()
            if freq_common.size == 0:
                # Try psd_lfp freq (should match coherence freq)
                psd_lfp = result.get('psd_lfp', {})
                if isinstance(psd_lfp, dict):
                    freq_common = np.array(psd_lfp.get('freq', [])).flatten()
            if freq_common.size == 0:
                # Try time_resolved
                tr = result.get('time_resolved', {})
                if isinstance(tr, dict):
                    freq_common = np.array(tr.get('freq', [])).flatten()
            
            if freq_common.size > 0:
                print(f"    DEBUG coherence: Using freq_common with {freq_common.size} points, range [{freq_common.min():.1f}, {freq_common.max():.1f}] Hz")
        
        # Try all possible period key names (MATLAB uses: pre_stim, transient, sustained, post_stim)
        # For stim period: Use sustained (more reliable than transient which is very short)
        # If sustained is missing, fall back to transient
        period_mappings = [
            (['pre_stim', 'prestim', 'pre'], coh_pre_all),
            # For stim: prefer sustained (longer, more reliable), fallback to transient
            (['sustained', 'stim_sustained', 'sust'], coh_stim_all),
            (['transient', 'stim_transient', 'trans'], coh_stim_all),  # Fallback if sustained missing
            (['post_stim', 'poststim', 'post'], coh_post_all)
        ]
        
        stim_found = False  # Track if we found stim period data
        for key_variants, storage in period_mappings:
            for key in key_variants:
                if key in coh_dict:
                    coh_data = np.array(coh_dict[key]).flatten()
                    # Remove NaN values and check if we have valid data
                    coh_data_clean = coh_data[np.isfinite(coh_data)]
                    if coh_data_clean.size > 0:
                        # Check if length matches freq_common
                        if freq_common is not None and len(coh_data) == len(freq_common):
                            # For stim period, only add if we haven't found it yet (prefer sustained over transient)
                            if storage is coh_stim_all and stim_found:
                                print(f"    DEBUG coherence: Skipped {key} - already have stim period data")
                                break
                            storage.append(coh_data)
                            if storage is coh_stim_all:
                                stim_found = True
                            print(f"    DEBUG coherence: Added {key} coherence ({coh_data.size} points), mean={np.nanmean(coh_data):.4f}, std={np.nanstd(coh_data):.4f}, min={np.nanmin(coh_data):.4f}, max={np.nanmax(coh_data):.4f}")
                        elif freq_common is None:
                            # Store anyway, we'll figure out freq later
                            storage.append(coh_data)
                            if storage is coh_stim_all:
                                stim_found = True
                        else:
                            print(f"    DEBUG coherence: Skipped {key} - length mismatch: coh={len(coh_data)}, freq={len(freq_common)}")
                    else:
                        print(f"    DEBUG coherence: Skipped {key} - no finite values")
                    break  # Only use first matching key
    
    print(f"    DEBUG coherence: freq_common size = {freq_common.size if freq_common is not None else 0}")
    print(f"    DEBUG coherence: pre={len(coh_pre_all)}, stim={len(coh_stim_all)}, post={len(coh_post_all)}")
    
    # Plot coherence spectra if we have data
    has_data = False
    if freq_common is not None and freq_common.size > 0:
        freq_mask = (freq_common >= FREQ_RANGE[0]) & (freq_common <= FREQ_RANGE[1])
        freq_plot = freq_common[freq_mask]
        
        def plot_period_coherence(coh_list, color, label):
            nonlocal has_data
            if len(coh_list) == 0:
                return
            # Ensure all arrays have same length as freq_common
            valid_coh = []
            for c in coh_list:
                if len(c) == len(freq_common):
                    valid_coh.append(c[freq_mask])
                elif len(c) == len(freq_plot):
                    valid_coh.append(c)
            
            if len(valid_coh) == 0:
                print(f"    DEBUG: No valid coherence data for {label}")
                return
            
            coh_arr = np.array(valid_coh)
            # Debug: Check if all values are the same (would cause straight line)
            if coh_arr.size > 0:
                unique_vals = np.unique(coh_arr[~np.isnan(coh_arr)])
                if len(unique_vals) <= 2:
                    print(f"    WARNING: {label} coherence has very few unique values: {unique_vals[:10]}")
                    print(f"    WARNING: {label} coherence shape: {coh_arr.shape}, mean per freq: {np.nanmean(coh_arr, axis=0)[:5]}")
            
            coh_mean = np.nanmean(coh_arr, axis=0)
            coh_sem = np.nanstd(coh_arr, axis=0) / np.sqrt(coh_arr.shape[0])
            
            # Apply smoothing to coherence mean and SEM for cleaner appearance
            # Use a moderate window (7 points) to reduce spikiness while preserving major features
            COHERENCE_SMOOTH_WINDOW = 7  # Adjust this value: larger = smoother
            if len(coh_mean) > COHERENCE_SMOOTH_WINDOW:
                coh_mean_smooth = uniform_filter1d(coh_mean, size=COHERENCE_SMOOTH_WINDOW, mode='nearest')
                coh_sem_smooth = uniform_filter1d(coh_sem, size=COHERENCE_SMOOTH_WINDOW, mode='nearest')
            else:
                coh_mean_smooth = coh_mean
                coh_sem_smooth = coh_sem
            
            # Debug: Print summary stats
            if label == 'Stim':
                print(f"    DEBUG stim coherence: mean across all freqs = {np.nanmean(coh_mean_smooth):.4f}, std = {np.nanstd(coh_mean_smooth):.4f}")
                print(f"    DEBUG stim coherence: min = {np.nanmin(coh_mean_smooth):.4f}, max = {np.nanmax(coh_mean_smooth):.4f}")
                print(f"    DEBUG stim coherence: first 5 values = {coh_mean_smooth[:5]}")
                print(f"    DEBUG stim coherence: last 5 values = {coh_mean_smooth[-5:]}")
            
            ax3.fill_between(freq_plot, coh_mean_smooth - coh_sem_smooth, coh_mean_smooth + coh_sem_smooth, 
                             color=color, alpha=0.3, linewidth=0)
            ax3.plot(freq_plot, coh_mean_smooth, color=color, linewidth=LINE_WIDTH_THICK, label=label)
            has_data = True
        
        plot_period_coherence(coh_pre_all, COLOR_PRE, 'Pre-stim')
        plot_period_coherence(coh_stim_all, COLOR_TRANSIENT, 'Stim')
        plot_period_coherence(coh_post_all, COLOR_POST, 'Post-stim')
        
        if has_data:
            ax3.set_xlim(FREQ_RANGE)
            ax3.set_ylim(0, 1)
            ax3.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
            ax3.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
            ax3.tick_params(labelsize=FONT_SIZE_TICK)
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            ax3.spines['left'].set_linewidth(AXIS_LINEWIDTH)
            ax3.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
            ax3.legend(loc='upper right', fontsize=FONT_SIZE_LEGEND, frameon=True, framealpha=0.9)
            ax3.set_title('Coherence Spectrum by Period (Mean ± SEM)', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    
    if not has_data:
        ax3.text(0.5, 0.5, 'No coherence data found\n(check MATLAB output structure)', 
                 ha='center', va='center', transform=ax3.transAxes, fontsize=FONT_SIZE_LABEL)
        ax3.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
        ax3.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
        ax3.set_xlim(FREQ_RANGE)
        ax3.set_ylim(0, 1)
    
    ax3.set_xlabel('Frequency (Hz)', fontsize=FONT_SIZE_LABEL)
    ax3.set_ylabel('LFP-GEVI Coherence', fontsize=FONT_SIZE_LABEL)
    ax3.set_title('Coherence Spectrum by Period (Mean ± SEM)', fontsize=FONT_SIZE_TITLE, fontweight='bold')
    ax3.tick_params(labelsize=FONT_SIZE_TICK)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax3.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax3.tick_params(width=TICK_WIDTH, length=TICK_LENGTH)
    
    # Adjust layout
    fig.subplots_adjust(left=0.12, right=0.88, top=0.94, bottom=0.06)
    
    return fig


def create_fiber_spectrogram_only_figure(session_id, condition_name, spectral_results, stim_onset_time):
    """
    Standalone fiber spectrogram panel (same normalization/format as heatmap figure).
    """
    fig = plt.figure(figsize=(8.4, 7.2))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1, 0.03], wspace=0.05)
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[0, 1])

    t_min_plot = -PRE_STIM_EFFECTIVE
    t_max_plot = 3.0  # stim onset + 3 seconds

    spectral_example = None
    for s in spectral_results:
        if s is not None:
            spectral_example = s
            break

    im = None
    if spectral_example and isinstance(spectral_example, dict) and 'spectrogram_fiber' in spectral_example:
        spec = spectral_example['spectrogram_fiber']
        if isinstance(spec, dict):
            power = np.array(spec.get('power', [])).squeeze()
            freq = np.array(spec.get('freq', [])).squeeze()
            time = np.array(spec.get('time', [])).squeeze() - stim_onset_time

            if power.size > 0 and power.ndim == 2:
                time_mask = (time >= t_min_plot) & (time <= t_max_plot)
                time = time[time_mask]
                power = power[:, time_mask]

                pre_mask = _spectrogram_baseline_time_mask(time)
                if np.any(pre_mask):
                    power_linear = 10**(power / 10.0)
                    baseline_linear = np.nanmean(power_linear[:, pre_mask], axis=1, keepdims=True)
                    numerator = power_linear - baseline_linear
                    denominator = power_linear + baseline_linear
                    denominator[denominator == 0] = 1e-10
                    power_norm = numerator / denominator
                else:
                    power_norm = power - np.nanmean(power)

                im = ax.pcolormesh(time, freq, power_norm, shading='auto', cmap='viridis')
                ax.set_ylim(FREQ_RANGE)

    if im is not None:
        cbar = fig.colorbar(im, cax=cax, orientation='vertical')
        cbar.set_label('Fractional Change\n(rel. baseline)', fontsize=FONT_SIZE_TICK)
        cbar.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    else:
        ax.text(0.5, 0.5, 'No fiber spectrogram data', transform=ax.transAxes,
                ha='center', va='center', fontsize=FONT_SIZE_LABEL)
        cax.set_visible(False)

    ax.axvline(0, color='white', linestyle='--', linewidth=1.5)
    ax.axvline(STIM_DURATION, color='white', linestyle=':', linewidth=1.5)
    ax.set_ylabel('Fiber\nFrequency (Hz)', fontsize=FONT_SIZE_LABEL, rotation=0, ha='right', va='center')
    ax.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL)
    ax.tick_params(labelsize=FONT_SIZE_TICK)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax.set_xlim(t_min_plot, t_max_plot)
    ax.set_title(f'{MOUSE_ID} - {condition_name} Fiber Spectrogram',
                 fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)

    fig.subplots_adjust(left=0.14, right=0.90, top=0.90, bottom=0.13)
    return fig


# =============================================================================
# FIGURE 5 & 6: PERIOD VIOLIN PLOTS
# =============================================================================

def create_period_violins_figure(session_id, condition_name, spectral_results):
    """
    Create figure with 4 violin plots for period comparison.
    
    Colors: grey (pre-stim), deep red (transient), maroon-orange (sustained), teal (post-stim)
    """
    # Increased spacing to accommodate significance brackets and prevent title overlap
    fig, axes = plt.subplots(2, 2, figsize=(16, 14), gridspec_kw={'hspace': 0.50, 'wspace': 0.35})
    
    # Color scheme per user request
    colors_period = {
        'Pre-stim': COLOR_PRE,        # grey
        'Transient': COLOR_TRANSIENT, # deep red
        'Sustained': COLOR_SUSTAINED, # maroon-orange
        'Post-stim': COLOR_POST       # teal-ish
    }
    
    # Subplot 1: Fiber signal (from raw traces)
    fiber_data = {
        'Pre-stim': extract_spectral_metric(spectral_results, 'pre_stim', 'fiber_mean', debug=True),
        'Transient': extract_spectral_metric(spectral_results, 'transient', 'fiber_mean'),
        'Sustained': extract_spectral_metric(spectral_results, 'sustained', 'fiber_mean'),
        'Post-stim': extract_spectral_metric(spectral_results, 'post_stim', 'fiber_mean')
    }
    # Statistical comparisons for period violins (paired - same trials across periods)
    period_comparisons = [
        ('Pre-stim', 'Transient', True),   # Paired: same trials
        ('Pre-stim', 'Sustained', True),
        ('Pre-stim', 'Post-stim', True),
        ('Transient', 'Sustained', True),
        ('Transient', 'Post-stim', True),
        ('Sustained', 'Post-stim', True)
    ]
    plot_violin_box(axes[0, 0], fiber_data, 'Fiber Signal (ΔF/F %)', colors_period, 
                    'Fiber Signal by Period', comparisons=period_comparisons,
                    correction='holm', omnibus=True)
    
    # Subplot 2: LFP power (baseline normalized)
    # PSD is stored in dB, so normalization is subtraction, not division
    pre_psd_lfp = extract_spectral_metric(spectral_results, 'pre_stim', 'psd_lfp_mean', debug=True)
    pre_mean_lfp = np.nanmean([v for v in pre_psd_lfp if v is not None]) if any(v is not None for v in pre_psd_lfp) else 0.0
    
    # Extract all periods with debug for sustained
    transient_lfp = extract_spectral_metric(spectral_results, 'transient', 'psd_lfp_mean', debug=False)
    sustained_lfp = extract_spectral_metric(spectral_results, 'sustained', 'psd_lfp_mean', debug=True)  # Debug to see why missing
    post_lfp = extract_spectral_metric(spectral_results, 'post_stim', 'psd_lfp_mean', debug=False)
    
    lfp_power_data = {
        'Pre-stim': [v - pre_mean_lfp if v is not None else None for v in pre_psd_lfp],  # dB subtraction
        'Transient': [v - pre_mean_lfp if v is not None else None for v in transient_lfp],  # dB subtraction
        'Sustained': [v - pre_mean_lfp if v is not None else None for v in sustained_lfp],  # dB subtraction
        'Post-stim': [v - pre_mean_lfp if v is not None else None for v in post_lfp]  # dB subtraction
    }
    plot_violin_box(axes[0, 1], lfp_power_data, 'LFP Power (dB re: baseline)', colors_period, 
                    'LFP Power (Baseline Normalized)', comparisons=period_comparisons,
                    correction='holm', omnibus=True)
    axes[0, 1].set_title('LFP Power (Baseline Normalized)', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Subplot 3: Fiber power (baseline normalized)
    # PSD is stored in dB, so normalization is subtraction, not division
    pre_psd_fiber = extract_spectral_metric(spectral_results, 'pre_stim', 'psd_fiber_mean', debug=True)
    pre_mean_fiber = np.nanmean([v for v in pre_psd_fiber if v is not None]) if any(v is not None for v in pre_psd_fiber) else 0.0
    
    # Extract all periods with debug for sustained
    transient_fiber = extract_spectral_metric(spectral_results, 'transient', 'psd_fiber_mean', debug=False)
    sustained_fiber = extract_spectral_metric(spectral_results, 'sustained', 'psd_fiber_mean', debug=True)  # Debug to see why missing
    post_fiber = extract_spectral_metric(spectral_results, 'post_stim', 'psd_fiber_mean', debug=False)
    
    fiber_power_data = {
        'Pre-stim': [v - pre_mean_fiber if v is not None else None for v in pre_psd_fiber],  # dB subtraction
        'Transient': [v - pre_mean_fiber if v is not None else None for v in transient_fiber],  # dB subtraction
        'Sustained': [v - pre_mean_fiber if v is not None else None for v in sustained_fiber],  # dB subtraction
        'Post-stim': [v - pre_mean_fiber if v is not None else None for v in post_fiber]  # dB subtraction
    }
    plot_violin_box(axes[1, 0], fiber_power_data, 'Fiber Power (dB re: baseline)', colors_period, 
                    'Fiber Power (Baseline Normalized)', comparisons=period_comparisons,
                    correction='holm', omnibus=True)
    axes[1, 0].set_title('Fiber Power (Baseline Normalized)', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Subplot 4: Coherence
    coh_data = {
        'Pre-stim': extract_spectral_metric(spectral_results, 'pre_stim', 'coherence_mean', debug=True),
        'Transient': extract_spectral_metric(spectral_results, 'transient', 'coherence_mean'),
        'Sustained': extract_spectral_metric(spectral_results, 'sustained', 'coherence_mean'),
        'Post-stim': extract_spectral_metric(spectral_results, 'post_stim', 'coherence_mean')
    }
    plot_violin_box(axes[1, 1], coh_data, 'Coherence', colors_period, 
                    'LFP-GEVI Coherence by Period', comparisons=period_comparisons,
                    correction='holm', omnibus=True)
    axes[1, 1].set_title('LFP-GEVI Coherence by Period', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    fig.suptitle(f'{MOUSE_ID} - {condition_name} Period Comparison', 
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.995)
    
    # Adjust layout with more top padding to prevent overlap with significance brackets
    fig.subplots_adjust(left=0.08, right=0.95, top=0.94, bottom=0.08)
    
    return fig


# =============================================================================
# FIGURE 7: 40Hz vs 135Hz COMPARISON
# =============================================================================

def create_comparison_figure(spectral_r1, spectral_r2):
    """
    Create figure comparing 40Hz vs 135Hz stimulation.
    
    5 subplots: Fiber signal, LFP power (overall), Fiber power (overall), 
                Stim-band power (40Hz band vs 135Hz band), Coherence
    Colors: 40Hz = teal shades, 135Hz = peach shades
    """
    # 2 rows x 3 columns to accommodate 5 subplots (one empty spot)
    fig, axes = plt.subplots(2, 3, figsize=(22, 14), gridspec_kw={'hspace': 0.50, 'wspace': 0.35})
    
    # Color scheme: 40Hz = teal, 135Hz = peach
    colors_freq = {
        '40Hz Trans': COLOR_40HZ,
        '40Hz Sust': COLOR_40HZ_LIGHT,
        '135Hz Trans': COLOR_135HZ,
        '135Hz Sust': COLOR_135HZ_LIGHT
    }
    
    # Subplot 1: Fiber signal
    fiber_data = {
        '40Hz Trans': extract_spectral_metric(spectral_r1, 'transient', 'fiber_mean', debug=True),
        '40Hz Sust': extract_spectral_metric(spectral_r1, 'sustained', 'fiber_mean'),
        '135Hz Trans': extract_spectral_metric(spectral_r2, 'transient', 'fiber_mean'),
        '135Hz Sust': extract_spectral_metric(spectral_r2, 'sustained', 'fiber_mean')
    }
    # Statistical comparisons for 40Hz vs 135Hz (unpaired - different sessions)
    # Also compare Trans vs Sust within each frequency (paired - same trials)
    freq_comparisons = [
        ('40Hz Trans', '40Hz Sust', True),      # Paired: same trials, different periods
        ('135Hz Trans', '135Hz Sust', True),    # Paired: same trials, different periods
        ('40Hz Trans', '135Hz Trans', False),   # Unpaired: different sessions
        ('40Hz Sust', '135Hz Sust', False)      # Unpaired: different sessions
    ]
    plot_violin_box(axes[0, 0], fiber_data, 'Fiber Signal (ΔF/F %)', colors_freq, 
                    'Fiber Signal: 40Hz vs 135Hz', comparisons=freq_comparisons,
                    correction='holm')
    
    # Subplot 2: LFP power (baseline normalized) - NEW
    # PSD is stored in dB, so normalization is subtraction, not division
    pre_lfp_40 = extract_spectral_metric(spectral_r1, 'pre_stim', 'psd_lfp_mean')
    pre_lfp_135 = extract_spectral_metric(spectral_r2, 'pre_stim', 'psd_lfp_mean')
    pre_mean_lfp_40 = np.nanmean([v for v in pre_lfp_40 if v is not None]) if any(v is not None for v in pre_lfp_40) else 0.0
    pre_mean_lfp_135 = np.nanmean([v for v in pre_lfp_135 if v is not None]) if any(v is not None for v in pre_lfp_135) else 0.0
    
    trans_lfp_40 = extract_spectral_metric(spectral_r1, 'transient', 'psd_lfp_mean')
    sust_lfp_40 = extract_spectral_metric(spectral_r1, 'sustained', 'psd_lfp_mean', debug=True)
    trans_lfp_135 = extract_spectral_metric(spectral_r2, 'transient', 'psd_lfp_mean')
    sust_lfp_135 = extract_spectral_metric(spectral_r2, 'sustained', 'psd_lfp_mean', debug=True)
    
    lfp_power_data = {
        '40Hz Trans': [v - pre_mean_lfp_40 if v is not None else None for v in trans_lfp_40],  # dB subtraction
        '40Hz Sust': [v - pre_mean_lfp_40 if v is not None else None for v in sust_lfp_40],  # dB subtraction
        '135Hz Trans': [v - pre_mean_lfp_135 if v is not None else None for v in trans_lfp_135],  # dB subtraction
        '135Hz Sust': [v - pre_mean_lfp_135 if v is not None else None for v in sust_lfp_135]  # dB subtraction
    }
    plot_violin_box(axes[0, 1], lfp_power_data, 'LFP Power (norm.)', colors_freq, 
                    None, comparisons=freq_comparisons, correction='holm')
    axes[0, 1].set_title('LFP Power: 40Hz vs 135Hz', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Subplot 3: Fiber power (baseline normalized)
    pre_fiber_40 = extract_spectral_metric(spectral_r1, 'pre_stim', 'psd_fiber_mean')
    pre_fiber_135 = extract_spectral_metric(spectral_r2, 'pre_stim', 'psd_fiber_mean')
    pre_mean_fiber_40 = np.nanmean([v for v in pre_fiber_40 if v is not None]) if any(v is not None for v in pre_fiber_40) else 0.0
    pre_mean_fiber_135 = np.nanmean([v for v in pre_fiber_135 if v is not None]) if any(v is not None for v in pre_fiber_135) else 0.0
    
    trans_fiber_40 = extract_spectral_metric(spectral_r1, 'transient', 'psd_fiber_mean')
    sust_fiber_40 = extract_spectral_metric(spectral_r1, 'sustained', 'psd_fiber_mean', debug=True)
    trans_fiber_135 = extract_spectral_metric(spectral_r2, 'transient', 'psd_fiber_mean')
    sust_fiber_135 = extract_spectral_metric(spectral_r2, 'sustained', 'psd_fiber_mean', debug=True)
    
    fiber_power_data = {
        '40Hz Trans': [v - pre_mean_fiber_40 if v is not None else None for v in trans_fiber_40],  # dB subtraction
        '40Hz Sust': [v - pre_mean_fiber_40 if v is not None else None for v in sust_fiber_40],  # dB subtraction
        '135Hz Trans': [v - pre_mean_fiber_135 if v is not None else None for v in trans_fiber_135],  # dB subtraction
        '135Hz Sust': [v - pre_mean_fiber_135 if v is not None else None for v in sust_fiber_135]  # dB subtraction
    }
    plot_violin_box(axes[0, 2], fiber_power_data, 'Fiber Power (norm.)', colors_freq, 
                    None, comparisons=freq_comparisons, correction='holm')
    axes[0, 2].set_title('Fiber Power (Overall): 40Hz vs 135Hz', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Subplot 4: FIBER STIM-BAND SPECIFIC POWER
    # Extract FIBER power specifically in the stimulation frequency band
    # 40Hz stim: extract 35-45 Hz band power from fiber
    # 135Hz stim: extract 130-140 Hz band power from fiber
    BAND_40HZ = (35, 45)
    BAND_135HZ = (130, 140)
    
    # Get baseline FIBER power for normalization
    pre_band_fiber_40 = extract_band_power(spectral_r1, 'pre_stim', 'fiber', BAND_40HZ, debug=True)
    pre_band_fiber_135 = extract_band_power(spectral_r2, 'pre_stim', 'fiber', BAND_135HZ, debug=True)
    pre_mean_band_40 = np.nanmean([v for v in pre_band_fiber_40 if v is not None]) if any(v is not None for v in pre_band_fiber_40) else 0.0
    pre_mean_band_135 = np.nanmean([v for v in pre_band_fiber_135 if v is not None]) if any(v is not None for v in pre_band_fiber_135) else 0.0
    
    # Extract FIBER stim-band power for transient and sustained periods
    trans_band_40 = extract_band_power(spectral_r1, 'transient', 'fiber', BAND_40HZ)
    sust_band_40 = extract_band_power(spectral_r1, 'sustained', 'fiber', BAND_40HZ, debug=True)
    trans_band_135 = extract_band_power(spectral_r2, 'transient', 'fiber', BAND_135HZ)
    sust_band_135 = extract_band_power(spectral_r2, 'sustained', 'fiber', BAND_135HZ, debug=True)
    
    stim_band_data = {
        '40Hz Trans': [v - pre_mean_band_40 if v is not None else None for v in trans_band_40],
        '40Hz Sust': [v - pre_mean_band_40 if v is not None else None for v in sust_band_40],
        '135Hz Trans': [v - pre_mean_band_135 if v is not None else None for v in trans_band_135],
        '135Hz Sust': [v - pre_mean_band_135 if v is not None else None for v in sust_band_135]
    }
    plot_violin_box(axes[1, 0], stim_band_data, 'Band Power (dB re: baseline)', colors_freq, 
                    None, comparisons=freq_comparisons, correction='holm')
    axes[1, 0].set_title('Fiber Stim-Band Power\n(40Hz:35-45Hz, 135Hz:130-140Hz)', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Subplot 5: Coherence
    coh_data = {
        '40Hz Trans': extract_spectral_metric(spectral_r1, 'transient', 'coherence_mean'),
        '40Hz Sust': extract_spectral_metric(spectral_r1, 'sustained', 'coherence_mean'),
        '135Hz Trans': extract_spectral_metric(spectral_r2, 'transient', 'coherence_mean'),
        '135Hz Sust': extract_spectral_metric(spectral_r2, 'sustained', 'coherence_mean')
    }
    plot_violin_box(axes[1, 1], coh_data, 'Coherence', colors_freq, 
                    None, comparisons=freq_comparisons, correction='holm')
    axes[1, 1].set_title('LFP-GEVI Coherence: 40Hz vs 135Hz', fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=15)
    
    # Hide the 6th subplot (axes[1, 2]) - not needed
    axes[1, 2].axis('off')
    
    fig.suptitle(f'{MOUSE_ID} - 40Hz vs 135Hz Stimulation Comparison', 
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.995)
    
    # Adjust layout with more top padding to prevent overlap with significance brackets
    fig.subplots_adjust(left=0.06, right=0.96, top=0.94, bottom=0.08)
    
    return fig


# =============================================================================
# FIGURE 8: TRIAL-BY-TRIAL HEATMAPS
# =============================================================================

def create_trials_heatmap_figure(session_id, condition_name, raw_trials, stim_onset_time, stim_freq, spectral_results=None):
    """
    Create trial-by-trial heatmap figure with THREE subplots per session:
    1. Baseline-subtracted fiber signal (binned, like MATLAB RMS approach)
    2. Band power in stimulation frequency band (from spectrogram data)
    3. Time course of averaged band power with SEM shading (NEW!)
    
    Uses spectrogram time resolution for intelligent binning of band power.
    
    Args:
        session_id: Session identifier (e.g., "04_09_25-R1")
        condition_name: Condition name for title (e.g., "40Hz")
        raw_trials: List of trial data dictionaries
        stim_onset_time: Stimulation onset time (seconds)
        stim_freq: Stimulation frequency (Hz) for band power
        spectral_results: List of spectral result dictionaries (optional, for band power)
    """
    fig = plt.figure(figsize=(16, 14))  # Taller to accommodate 3 subplots
    gs = GridSpec(3, 1, figure=fig, height_ratios=[1, 1, 0.7], hspace=0.30)
    
    # Filter valid trials
    valid_trials = [t for t in raw_trials if t is not None]
    if not valid_trials:
        print(f"    Warning: No valid trials for {session_id}")
        return fig
    
    num_trials = len(valid_trials)
    print(f"    Processing {num_trials} trials for {session_id}")
    
    # Define time range for plotting (skip first second if configured)
    t_min_plot = -PRE_STIM_EFFECTIVE if PRE_STIM_SKIP_SEC > 0 else -PRE_STIM_DURATION
    t_max_plot = STIM_DURATION + POST_STIM_DURATION
    
    # Create time bins (0.1s bins like MATLAB)
    bin_size = 0.1
    time_bins = np.arange(t_min_plot, t_max_plot + bin_size, bin_size)
    n_bins = len(time_bins) - 1
    time_bin_centers = time_bins[:-1] + bin_size / 2  # Bin centers for x-axis
    
    # Initialize matrices for heatmaps
    # Shape: [num_trials, n_bins]
    signal_matrix = np.full((num_trials, n_bins), np.nan)
    bandpower_matrix = np.full((num_trials, n_bins), np.nan)
    
    # Define frequency band for band power (same for all trials)
    # Band: stim_freq ± 5 Hz
    lowcut = max(1.0, stim_freq - 5.0)
    highcut = stim_freq + 5.0  # Will be clipped per trial based on Nyquist
    
    # Process each trial
    for trial_idx, trial_data in enumerate(valid_trials):
        gevi = trial_data['gevi']
        t_trial = trial_data['t']
        fs_trial = trial_data['fs']
        
        # Clip highcut to Nyquist for this trial
        highcut_trial = min(fs_trial / 2.1, highcut)
        
        # Align time to stim onset
        t_aligned = t_trial - stim_onset_time
        
        # Get time mask for this trial
        time_mask_trial = (t_aligned >= t_min_plot) & (t_aligned <= t_max_plot)
        gevi_plot = gevi[time_mask_trial]
        t_trial_plot = t_aligned[time_mask_trial]
        
        # === HEATMAP 1: Baseline-Subtracted Signal (binned) ===
        # Calculate baseline as mean of pre-stim period (excluding first second if configured)
        pre_mask = t_trial_plot < 0
        if np.any(pre_mask):
            baseline = np.nanmean(gevi_plot[pre_mask])
        else:
            baseline = np.nanmean(gevi_plot)
        
        # Bin the signal: compute mean per bin (like RMS in MATLAB)
        for b in range(n_bins):
            bin_start = time_bins[b]
            bin_end = time_bins[b + 1]
            
            # Find samples in this bin
            bin_mask = (t_trial_plot >= bin_start) & (t_trial_plot < bin_end)
            if np.any(bin_mask):
                # Compute mean of signal in this bin, then baseline subtract
                bin_signal = np.nanmean(gevi_plot[bin_mask])
                signal_matrix[trial_idx, b] = bin_signal - baseline
        
        # === HEATMAP 2: Band Power (from spectrogram) ===
        # Use spectrogram data if available, otherwise compute from raw signal
        if spectral_results and trial_idx < len(spectral_results) and spectral_results[trial_idx] is not None:
            spec_data = spectral_results[trial_idx].get('spectrogram_fiber', {})
            if isinstance(spec_data, dict):
                spec_power = np.array(spec_data.get('power', [])).squeeze()
                spec_freq = np.array(spec_data.get('freq', [])).squeeze()
                spec_time = np.array(spec_data.get('time', [])).squeeze() - stim_onset_time
                
                if spec_power.size > 0 and spec_power.ndim == 2 and len(spec_freq) > 0 and len(spec_time) > 0:
                    # Filter to time range
                    time_mask_spec = (spec_time >= t_min_plot) & (spec_time <= t_max_plot)
                    if np.any(time_mask_spec):
                        spec_power_plot = spec_power[:, time_mask_spec]
                        spec_time_plot = spec_time[time_mask_spec]
                        
                        # Extract power in frequency band
                        freq_mask_band = (spec_freq >= lowcut) & (spec_freq <= highcut_trial)
                        if np.any(freq_mask_band):
                            # Average power across frequency band (mean across frequencies in band)
                            # Power is already in dB from MATLAB
                            band_power_time = np.nanmean(spec_power_plot[freq_mask_band, :], axis=0)
                            
                            # Baseline normalize (subtract in dB space)
                            pre_mask_spec = spec_time_plot < 0
                            if np.any(pre_mask_spec):
                                baseline_power = np.nanmean(band_power_time[pre_mask_spec])
                            else:
                                baseline_power = np.nanmean(band_power_time)
                            
                            # Baseline normalize (subtract in dB space)
                            band_power_norm = band_power_time - baseline_power
                            
                            # Determine intelligent bin size based on spectrogram time resolution
                            if len(spec_time_plot) > 1:
                                # Use spectrogram's native time resolution
                                spec_time_step = np.nanmean(np.diff(spec_time_plot))
                                # Use bin size that's 2-3x the spectrogram time step for good visualization
                                # But not too fine (min 0.05s) or too coarse (max 0.2s)
                                intelligent_bin_size = max(0.05, min(0.2, spec_time_step * 2.5))
                            else:
                                intelligent_bin_size = bin_size  # Fallback to 0.1s
                            
                            # Interpolate to bin centers using intelligent binning
                            for b in range(n_bins):
                                bin_center = time_bin_centers[b]
                                # Find time points within this bin (using intelligent bin size)
                                bin_start = bin_center - intelligent_bin_size / 2
                                bin_end = bin_center + intelligent_bin_size / 2
                                bin_mask_spec = (spec_time_plot >= bin_start) & (spec_time_plot < bin_end)
                                
                                if np.any(bin_mask_spec):
                                    # Average power within this bin
                                    bandpower_matrix[trial_idx, b] = np.nanmean(band_power_norm[bin_mask_spec])
                                elif len(spec_time_plot) > 0:
                                    # Fallback: use nearest neighbor
                                    nearest_idx = np.argmin(np.abs(spec_time_plot - bin_center))
                                    if abs(spec_time_plot[nearest_idx] - bin_center) < intelligent_bin_size:
                                        bandpower_matrix[trial_idx, b] = band_power_norm[nearest_idx]
        
        # Fallback: compute from raw signal if spectrogram not available
        if np.all(np.isnan(bandpower_matrix[trial_idx, :])):
            for b in range(n_bins):
                bin_start = time_bins[b]
                bin_end = time_bins[b + 1]
                
                # Find samples in this bin
                bin_mask = (t_trial_plot >= bin_start) & (t_trial_plot < bin_end)
                if np.any(bin_mask) and np.sum(bin_mask) >= 32:  # Need enough samples for PSD
                    bin_data = gevi_plot[bin_mask]
                    
                    # Compute PSD using Welch's method for this bin
                    try:
                        # Use shorter window for short bins
                        nperseg = min(64, len(bin_data) // 2)
                        if nperseg < 16:
                            nperseg = len(bin_data)
                        
                        freqs, psd = signal.welch(bin_data, fs=fs_trial, nperseg=nperseg, 
                                                 noverlap=nperseg//2, window='hanning')
                        
                        # Extract power in the frequency band
                        freq_mask = (freqs >= lowcut) & (freqs <= highcut_trial)
                        if np.any(freq_mask):
                            # Integrate power in the band (sum of PSD)
                            band_power = np.trapz(psd[freq_mask], freqs[freq_mask])
                            # Convert to dB and baseline normalize
                            band_power_db = 10 * np.log10(band_power + 1e-10)  # Add small value to avoid log(0)
                            # Baseline: use pre-stim bins
                            pre_bins = [i for i in range(n_bins) if time_bin_centers[i] < 0]
                            if pre_bins:
                                baseline_power = np.nanmean([bandpower_matrix[trial_idx, i] for i in pre_bins if not np.isnan(bandpower_matrix[trial_idx, i])])
                                if not np.isnan(baseline_power):
                                    band_power_db = band_power_db - baseline_power
                            bandpower_matrix[trial_idx, b] = band_power_db
                    except:
                        # If PSD fails, leave as NaN
                        pass
    
    # === SUBPLOT 1: Baseline-Subtracted Signal ===
    ax1 = fig.add_subplot(gs[0])
    
    # Use symmetric colormap for signed data
    # Use parula-like colormap for signal strength heatmap (matches example figure style)
    # Note: values are baseline-subtracted, so range can be negative
    vmax1 = np.nanpercentile(np.abs(signal_matrix), 95)
    vmin1 = -vmax1  # Symmetric around zero
    
    im1 = ax1.imshow(signal_matrix, aspect='auto', origin='lower',
                     extent=[time_bin_centers[0], time_bin_centers[-1], 0.5, num_trials + 0.5],
                     cmap=CMAP_PARULA_LIKE, vmin=vmin1, vmax=vmax1, interpolation='nearest')
    
    # Mark stimulation period
    ax1.axvline(0, color='white', linestyle='--', linewidth=2, alpha=0.8)
    ax1.axvline(STIM_DURATION, color='white', linestyle=':', linewidth=2, alpha=0.8)
    
    ax1.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax1.set_ylabel('Trial Number', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax1.set_title(f'{condition_name} - Baseline-Subtracted Fiber Signal (ΔF/F)', 
                  fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    ax1.set_yticks(np.arange(1, num_trials + 1))
    ax1.set_yticklabels([f'T{i}' for i in range(1, num_trials + 1)])
    ax1.tick_params(labelsize=FONT_SIZE_TICK)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax1.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    
    cbar1 = fig.colorbar(im1, ax=ax1, orientation='vertical', pad=0.02)
    cbar1.set_label('ΔF/F (rel. baseline)', fontsize=FONT_SIZE_TICK)
    cbar1.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    
    # === SUBPLOT 2: Band Power ===
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    
    # Use parula-like colormap for band power heatmap (matches example figure style)
    # Power is in dB relative to baseline, so can be positive or negative
    vmax2 = np.nanpercentile(np.abs(bandpower_matrix[np.isfinite(bandpower_matrix)]), 95) if np.any(np.isfinite(bandpower_matrix)) else 1.0
    vmin2 = -vmax2  # Symmetric around zero
    
    im2 = ax2.imshow(bandpower_matrix, aspect='auto', origin='lower',
                     extent=[time_bin_centers[0], time_bin_centers[-1], 0.5, num_trials + 0.5],
                     cmap=CMAP_PARULA_LIKE, vmin=vmin2, vmax=vmax2, interpolation='nearest')
    
    # Mark stimulation period
    ax2.axvline(0, color='white', linestyle='--', linewidth=2, alpha=0.8)
    ax2.axvline(STIM_DURATION, color='white', linestyle=':', linewidth=2, alpha=0.8)
    
    ax2.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_ylabel('Trial Number', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax2.set_title(f'{condition_name} - Band Power ({lowcut:.0f}-{stim_freq+5:.0f} Hz)', 
                  fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    ax2.set_yticks(np.arange(1, num_trials + 1))
    ax2.set_yticklabels([f'T{i}' for i in range(1, num_trials + 1)])
    ax2.tick_params(labelsize=FONT_SIZE_TICK)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax2.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    
    cbar2 = fig.colorbar(im2, ax=ax2, orientation='vertical', pad=0.02)
    cbar2.set_label('Band Power (dB re: baseline)', fontsize=FONT_SIZE_TICK)
    cbar2.ax.tick_params(labelsize=FONT_SIZE_TICK-2)
    
    # === SUBPLOT 3: Time Course of Averaged Band Power with SEM ===
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    
    # Compute mean and SEM across trials for band power
    # bandpower_matrix shape: [num_trials, n_bins]
    if np.any(np.isfinite(bandpower_matrix)):
        # Compute statistics across trials
        mean_bandpower = np.nanmean(bandpower_matrix, axis=0)
        # SEM = std / sqrt(n), using valid (non-NaN) trials per bin
        n_valid = np.sum(np.isfinite(bandpower_matrix), axis=0)
        std_bandpower = np.nanstd(bandpower_matrix, axis=0)
        sem_bandpower = std_bandpower / np.sqrt(np.maximum(n_valid, 1))
        
        # Plot mean with SEM shading
        ax3.fill_between(time_bin_centers, mean_bandpower - sem_bandpower, mean_bandpower + sem_bandpower,
                         color=COLOR_GEVI, alpha=0.3, linewidth=0)
        ax3.plot(time_bin_centers, mean_bandpower, color=COLOR_GEVI, linewidth=LINE_WIDTH_THICK, label='Mean ± SEM')
        
        # Add reference line at 0 (baseline)
        ax3.axhline(0, color='gray', linestyle='--', linewidth=1.5, alpha=0.7)
        
        # Mark stimulation period
        ax3.axvline(0, color='black', linestyle='--', linewidth=2, label='Stim onset')
        ax3.axvline(STIM_DURATION, color='black', linestyle=':', linewidth=2, label='Stim offset')
    
    ax3.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax3.set_ylabel('Band Power\n(dB re: baseline)', fontsize=FONT_SIZE_LABEL, fontweight='bold')
    ax3.set_title(f'Average Band Power ({lowcut:.0f}-{stim_freq+5:.0f} Hz) - Mean ± SEM (n={num_trials} trials)', 
                  fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    ax3.tick_params(labelsize=FONT_SIZE_TICK)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_linewidth(AXIS_LINEWIDTH)
    ax3.spines['bottom'].set_linewidth(AXIS_LINEWIDTH)
    ax3.legend(loc='upper right', fontsize=FONT_SIZE_TICK-1)
    
    fig.suptitle(f'{MOUSE_ID} - {condition_name} Trial-by-Trial Analysis', 
                 fontsize=FONT_SIZE_SUPTITLE, fontweight='bold', y=0.98)
    
    fig.subplots_adjust(left=0.08, right=0.92, top=0.94, bottom=0.06)
    
    return fig


# =============================================================================
# STRIATUM DBS COMPOSITE (single-session mode)
# =============================================================================



def _load_preprocessed_mat_file(mat_path, fiber_index=0):
    """Generic loader for one FiberPhotometry_Analysis.mat from absolute path."""
    mat_path = Path(mat_path)
    mat_access = normalize_unc_path(str(mat_path), for_access=True)
    if not os.path.exists(mat_access):
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    def _safe_scalar(x, default=np.nan):
        try:
            arr = np.asarray(x).flatten()
            return float(arr[0]) if arr.size else float(default)
        except Exception:
            return float(default)

    try:
        with h5py.File(mat_access, "r") as f:
            root = f["FiberPhotometryAnalysis"]
            t = np.array(root["time"]["time_vector_seconds"][()]).flatten()
            fs = _safe_scalar(root["time"]["sampling_rate"][()], default=np.nan)
            if not np.isfinite(fs) or fs <= 0:
                fs = 1.0 / np.nanmedian(np.diff(t))

            tr = np.array(root["signals"]["final_processed_traces"][()])
            if tr.ndim == 2 and tr.shape[0] < tr.shape[1]:
                tr = tr.T
            if tr.ndim == 1:
                gevi = tr.flatten()
            else:
                fi = min(fiber_index, tr.shape[1] - 1)
                gevi = tr[:, fi]

            ephys = root.get("ephys", None)
            lfp = None
            motion = None
            stim_pulses = None
            if ephys is not None:
                for key in ("lfp_raw_aligned_HP", "lfp_raw_aligned_mPFC", "lfp_raw_aligned"):
                    if key in ephys:
                        lfp = np.array(ephys[key][()]).flatten()
                        break
                if "running_velocity_smooth" in ephys:
                    motion = np.array(ephys["running_velocity_smooth"][()]).flatten()
                elif "running_velocity" in ephys:
                    motion = np.array(ephys["running_velocity"][()]).flatten()
                if "stim_pulses" in ephys:
                    stim_pulses = np.array(ephys["stim_pulses"][()]).flatten()

            n = len(t)
            if lfp is None:
                lfp = np.full(n, np.nan, dtype=float)
            if motion is None:
                motion = np.full(n, np.nan, dtype=float)
            if stim_pulses is None:
                stim_pulses = np.zeros(n, dtype=float)
            n = min(len(t), len(gevi), len(lfp), len(motion), len(stim_pulses))
            return {
                "t": t[:n],
                "gevi": np.asarray(gevi[:n], dtype=float),
                "lfp": np.asarray(lfp[:n], dtype=float),
                "motion": np.asarray(motion[:n], dtype=float) * MOTION_TO_CM_PER_S,
                "stim_pulses": np.asarray(stim_pulses[:n], dtype=float),
                "fs": float(fs),
            }
    except Exception:
        pass

    mat = loadmat(mat_access, squeeze_me=False, struct_as_record=False)
    fp = mat["FiberPhotometryAnalysis"]
    t = np.asarray(fp.time.time_vector_seconds).flatten()
    fs = _safe_scalar(fp.time.sampling_rate, default=np.nan)
    if not np.isfinite(fs) or fs <= 0:
        fs = 1.0 / np.nanmedian(np.diff(t))
    tr = np.asarray(fp.signals.final_processed_traces)
    if tr.ndim == 2 and tr.shape[0] < tr.shape[1]:
        tr = tr.T
    if tr.ndim == 1:
        gevi = tr.flatten()
    else:
        fi = min(fiber_index, tr.shape[1] - 1)
        gevi = tr[:, fi]

    ephys = fp.ephys if hasattr(fp, "ephys") else None
    lfp = None
    motion = None
    stim_pulses = None
    if ephys is not None:
        for key in ("lfp_raw_aligned_HP", "lfp_raw_aligned_mPFC", "lfp_raw_aligned"):
            if hasattr(ephys, key):
                lfp = np.asarray(getattr(ephys, key)).flatten()
                break
        if hasattr(ephys, "running_velocity_smooth"):
            motion = np.asarray(ephys.running_velocity_smooth).flatten()
        elif hasattr(ephys, "running_velocity"):
            motion = np.asarray(ephys.running_velocity).flatten()
        if hasattr(ephys, "stim_pulses"):
            stim_pulses = np.asarray(ephys.stim_pulses).flatten()

    n = len(t)
    if lfp is None:
        lfp = np.full(n, np.nan, dtype=float)
    if motion is None:
        motion = np.full(n, np.nan, dtype=float)
    if stim_pulses is None:
        stim_pulses = np.zeros(n, dtype=float)
    n = min(len(t), len(gevi), len(lfp), len(motion), len(stim_pulses))
    return {
        "t": t[:n],
        "gevi": np.asarray(gevi[:n], dtype=float),
        "lfp": np.asarray(lfp[:n], dtype=float),
        "motion": np.asarray(motion[:n], dtype=float) * MOTION_TO_CM_PER_S,
        "stim_pulses": np.asarray(stim_pulses[:n], dtype=float),
        "fs": float(fs),
    }


def _compute_trial_band_power(sig, fs, band_hz):
    x = np.asarray(sig, dtype=float)
    x = x[np.isfinite(x)]
    # Allow short transient windows (e.g., 150 ms) while still requiring enough samples.
    if len(x) < max(16, int(0.08 * fs)):
        return np.nan
    nper = min(256, len(x))
    nover = min(nper // 2, nper - 1)
    f, pxx = signal.welch(x, fs=fs, window="hann", nperseg=nper, noverlap=nover, detrend="linear")
    mk = np.isfinite(f) & np.isfinite(pxx) & (f >= band_hz[0]) & (f <= band_hz[1])
    if not np.any(mk):
        return np.nan
    return float(np.nanmean(pxx[mk]))


def _plot_three_period_violin_with_stats(ax, data_dict, ylabel, title=None):
    order = ["Transient", "Sustained", "Post-stim"]
    colors = {"Transient": COLOR_TRANSIENT, "Sustained": COLOR_SUSTAINED, "Post-stim": COLOR_POST}
    data_use = {k: data_dict.get(k, []) for k in order}
    plot_violin_box(
        ax,
        data_use,
        ylabel=ylabel,
        colors_dict=colors,
        title=title,
        comparisons=None,
        correction="holm",
        omnibus=True,
    )
    ax.tick_params(axis="x", labelrotation=0)
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(0)
        lbl.set_horizontalalignment("center")
    ax.axhline(0, color="black", linestyle=":", linewidth=1.5, alpha=0.85, zorder=0)

    # Baseline-vs-zero stars/ns above each violin (first layer).
    def _p_to_stars_local(p):
        if p < 0.001:
            return "***"
        if p < 0.01:
            return "**"
        if p < 0.05:
            return "*"
        return "ns"

    vals_all = []
    for lab in order:
        vals = [v for v in data_use.get(lab, []) if v is not None and np.isfinite(v)]
        vals_all.extend(vals)
    if len(vals_all) > 0:
        y_min = float(np.nanmin(vals_all))
        y_max = float(np.nanmax(vals_all))
        y_rng = max(y_max - y_min, 1e-9)
        base_y = y_max + 0.08 * y_rng
        for xi, lab in enumerate(order):
            vals = np.asarray([v for v in data_use.get(lab, []) if v is not None and np.isfinite(v)], dtype=float)
            if len(vals) < 3:
                continue
            try:
                _, p_norm = stats.shapiro(vals) if len(vals) <= 5000 else (None, 0.05)
                if p_norm > 0.05:
                    _, p_zero = stats.ttest_1samp(vals, popmean=0.0, nan_policy="omit")
                else:
                    _, p_zero = stats.wilcoxon(vals, alternative="two-sided")
            except Exception:
                p_zero = 1.0
            stars = _p_to_stars_local(float(p_zero))
            if stars != "ns":
                y_here = base_y + (0.045 * y_rng) * (xi % 2)
                ax.text(xi, y_here, stars, ha="center", va="bottom", fontsize=FONT_SIZE_TICK + 28, fontweight="bold")

        # Pairwise brackets and stars/ns (second layer above baseline symbols).
        comps = [("Transient", "Sustained"), ("Transient", "Post-stim"), ("Sustained", "Post-stim")]
        raw = []
        for a, b in comps:
            p, _, _, _, _ = perform_statistical_test(data_use.get(a, []), data_use.get(b, []), paired=True)
            raw.append((a, b, float(p)))
        corr = holm_bonferroni_correction([r[2] for r in raw]) if len(raw) > 1 else [r[2] for r in raw]
        pair_base = y_max + 0.34 * y_rng
        pair_step = 0.20 * y_rng
        line_h = 0.03 * y_rng
        for i, (a, b, _) in enumerate(raw):
            x1 = order.index(a)
            x2 = order.index(b)
            p_adj = corr[i]
            txt = _p_to_stars_local(p_adj)
            y = pair_base + i * pair_step
            ax.plot([x1, x2], [y, y], "k-", linewidth=3.2, clip_on=False)
            ax.plot([x1, x1], [y - line_h, y], "k-", linewidth=3.2, clip_on=False)
            ax.plot([x2, x2], [y - line_h, y], "k-", linewidth=3.2, clip_on=False)
            ax.text((x1 + x2) / 2.0, y + 0.028 * y_rng, txt, ha="center", va="bottom", fontsize=FONT_SIZE_TICK + 26, fontweight="bold")

        ax.set_ylim(y_min - 0.12 * y_rng, y_max + 1.18 * y_rng)


def create_striatum_dbs_composite_figure(trial1_mat_path, num_trials=10, invert_fiber=True):
    """
    Composite Striatum DBS figure (A-F) for a single 40Hz session.
    """
    trial1 = Path(trial1_mat_path)
    if not trial1.exists():
        raise FileNotFoundError(f"Trial1 mat not found: {trial1_mat_path}")
    session_dir = trial1.parent
    all_mats = sorted(
        session_dir.glob("*_Trial*_FiberPhotometry_Analysis.mat"),
        key=lambda p: (_infer_trial_from_name(p.name) if _infer_trial_from_name(p.name) is not None else 10**9),
    )
    if len(all_mats) == 0:
        raise FileNotFoundError(f"No trial mats found in session folder: {session_dir}")
    all_mats = all_mats[: int(num_trials)]

    trials = []
    for p in all_mats:
        try:
            d = _load_preprocessed_mat_file(str(p), fiber_index=0)
            d["trial_file"] = p.name
            d["trial_num"] = _infer_trial_from_name(p.name)
            trials.append(d)
            print(f"  Loaded: {p.name}")
        except Exception as e:
            print(f"  Warning loading {p.name}: {e}")
    if len(trials) == 0:
        raise RuntimeError("No valid trials loaded.")

    pre_sec = 5.0
    stim_sec = 1.0
    post_sec = 5.0
    transient_end = 0.15
    sustained_start = 0.15
    band_40 = (35.0, 45.0)

    rep = trials[0]
    for tr in trials:
        if tr.get("trial_num", None) == 1:
            rep = tr
            break

    fs_rep = rep["fs"]
    fiber_sign = -1.0 if invert_fiber else 1.0
    stim_on_idx = int(np.argmin(np.abs(rep["t"] - (rep["t"][0] + pre_sec))))
    t_rel = rep["t"] - rep["t"][stim_on_idx]
    fiber_pct = fiber_sign * rep["gevi"] * 100.0
    motion = rep["motion"]
    stim_trace = generate_stim_pulses(t_rel, 40.0, stim_sec)

    t_grid = np.linspace(-pre_sec, stim_sec + post_sec, 2500)
    fiber_all = []
    period_fiber = {"Transient": [], "Sustained": [], "Post-stim": []}
    period_band_db = {"Transient": [], "Sustained": [], "Post-stim": []}
    spec_stack = []

    for tr in trials:
        fs = tr["fs"]
        sidx = int(np.argmin(np.abs(tr["t"] - (tr["t"][0] + pre_sec))))
        tt = tr["t"] - tr["t"][sidx]
        yy = fiber_sign * tr["gevi"] * 100.0
        keep = (tt >= -pre_sec) & (tt <= (stim_sec + post_sec))
        if np.sum(keep) < 100:
            continue
        tt = tt[keep]
        yy = yy[keep]
        yi = np.interp(t_grid, tt, yy, left=np.nan, right=np.nan)
        fiber_all.append(yi)

        bmask = (tt >= -pre_sec) & (tt < 0.0)
        baseline = np.nanmean(yy[bmask]) if np.any(bmask) else np.nan
        if np.isfinite(baseline):
            tmask = (tt >= 0.0) & (tt < transient_end)
            smask = (tt >= sustained_start) & (tt < stim_sec)
            pmask = (tt >= stim_sec) & (tt < (stim_sec + post_sec))
            if np.any(tmask):
                period_fiber["Transient"].append(float(np.nanmean(yy[tmask]) - baseline))
            if np.any(smask):
                period_fiber["Sustained"].append(float(np.nanmean(yy[smask]) - baseline))
            if np.any(pmask):
                period_fiber["Post-stim"].append(float(np.nanmean(yy[pmask]) - baseline))

        bpow = _compute_trial_band_power(yy[(tt >= -pre_sec) & (tt < 0.0)], fs, band_40)
        if np.isfinite(bpow) and bpow > 0:
            tpow = _compute_trial_band_power(yy[(tt >= 0.0) & (tt < transient_end)], fs, band_40)
            spow = _compute_trial_band_power(yy[(tt >= sustained_start) & (tt < stim_sec)], fs, band_40)
            ppow = _compute_trial_band_power(yy[(tt >= stim_sec) & (tt < (stim_sec + post_sec))], fs, band_40)
            if np.isfinite(tpow) and tpow > 0:
                period_band_db["Transient"].append(float(10 * np.log10(tpow / bpow)))
            if np.isfinite(spow) and spow > 0:
                period_band_db["Sustained"].append(float(10 * np.log10(spow / bpow)))
            if np.isfinite(ppow) and ppow > 0:
                period_band_db["Post-stim"].append(float(10 * np.log10(ppow / bpow)))

        # Higher spectral resolution than default, then light display smoothing.
        nper = min(512, len(yy))
        if nper < 128:
            continue
        nover = min(int(0.875 * nper), nper - 1)
        nfft = max(2048, int(2 ** np.ceil(np.log2(nper))))
        f, ts, sxx = signal.spectrogram(
            np.nan_to_num(yy - np.nanmean(yy), nan=0.0),
            fs=fs,
            window="hann",
            nperseg=nper,
            noverlap=nover,
            nfft=nfft,
            detrend="linear",
            scaling="density",
            mode="psd",
        )
        ts_rel = ts + (tt[0] if len(tt) else 0.0)
        fmk = (f >= 1.0) & (f <= 100.0)
        if np.any(fmk):
            f2 = f[fmk]
            s2 = sxx[fmk, :]
            bmk = (ts_rel >= -pre_sec) & (ts_rel < -1.0)
            if np.sum(bmk) == 0:
                bmk = ts_rel < 0.0
            if np.sum(bmk) > 0:
                base = np.nanmean(s2[:, bmk], axis=1, keepdims=True)
                base = np.where(np.isfinite(base) & (base > 0), base, np.nan)
                s_rel = 10 * np.log10(np.maximum(s2, 1e-15) / np.maximum(base, 1e-15))
                if np.any(np.isfinite(s_rel)):
                    spec_stack.append((f2, ts_rel, s_rel))

    fiber_mat = np.asarray(fiber_all, dtype=float) if len(fiber_all) else np.empty((0, len(t_grid)))
    fiber_mu = np.nanmean(fiber_mat, axis=0) if fiber_mat.size else np.full(len(t_grid), np.nan)
    if fiber_mat.shape[0] > 1:
        fiber_sem = np.nanstd(fiber_mat, axis=0, ddof=1) / np.sqrt(np.sum(np.isfinite(fiber_mat), axis=0))
    else:
        fiber_sem = np.full(len(t_grid), np.nan)

    spec_avg = None
    if len(spec_stack):
        f_ref = spec_stack[0][0]
        t_ref = spec_stack[0][1]
        mats = []
        for ff, tt, ss in spec_stack:
            if len(ff) != len(f_ref) or np.max(np.abs(ff - f_ref)) > 1e-9:
                continue
            arr = np.full((len(f_ref), len(t_ref)), np.nan)
            for i in range(len(f_ref)):
                vv = ss[i, :]
                mk = np.isfinite(tt) & np.isfinite(vv)
                if np.sum(mk) >= 2:
                    arr[i, :] = np.interp(t_ref, tt[mk], vv[mk], left=np.nan, right=np.nan)
            mats.append(arr)
        if len(mats):
            spec_avg = (f_ref, t_ref, np.nanmean(np.asarray(mats), axis=0))

    fig = plt.figure(figsize=(42, 30))
    gs = GridSpec(3, 2, figure=fig, height_ratios=[1.1, 1.1, 1.6], width_ratios=[1.12, 1.0], hspace=0.56, wspace=0.30)
    fs_label = FONT_SIZE_LABEL + 28
    fs_tick = FONT_SIZE_TICK + 28

    # Panel A
    axA = fig.add_subplot(gs[0, 0])
    mA = (t_rel >= -pre_sec) & (t_rel <= (stim_sec + post_sec))
    tA = t_rel[mA]
    y_f = fiber_pct[mA] - np.nanmedian(fiber_pct[mA])
    if np.any(np.isfinite(motion[mA])):
        y_m = motion[mA] - np.nanmedian(motion[mA])
    else:
        y_m = np.zeros_like(tA)
    y_s = stim_trace[mA]
    nf = y_f / (np.nanpercentile(np.abs(y_f), 95) + 1e-9)
    nm = y_m / (np.nanpercentile(np.abs(y_m), 95) + 1e-9) if np.any(np.isfinite(y_m)) else np.zeros_like(tA)
    ns = y_s / (np.nanpercentile(np.abs(y_s), 95) + 1e-9) if np.any(np.isfinite(y_s)) else np.zeros_like(tA)
    # Larger inter-trace vertical spacing (Stim/Fiber/Motion).
    off_f, off_m, off_s = 0.0, -2.8, 2.8
    axA.plot(tA, ns * 0.6 + off_s, color=COLOR_STIM_PULSE, linewidth=1.5)
    axA.plot(tA, nf * 0.9 + off_f, color=COLOR_GEVI, linewidth=2.1)
    motion_color_dark = np.clip(np.asarray(COLOR_MOTION) * 0.58, 0.0, 1.0)
    axA.plot(tA, nm * 0.9 + off_m, color=motion_color_dark, linewidth=1.6)
    axA.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axA.axvline(stim_sec, color="k", linestyle=":", linewidth=1.2, alpha=0.8)
    axA.set_xlim(-3.0, 4.0)
    axA.set_ylim(-4.3, 4.3)
    axA.set_yticks([off_m, off_f, off_s])
    axA.set_yticklabels(["Motion", "Fiber Vm", "Stim"], fontsize=fs_tick)
    axA.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")

    # Panel B
    axB = fig.add_subplot(gs[0, 1])
    mB = (t_rel >= -0.3) & (t_rel <= 0.3)
    tB = t_rel[mB]
    y_fb = fiber_pct[mB] - np.nanmedian(fiber_pct[mB])
    y_sb = stim_trace[mB]
    nfb = y_fb / (np.nanpercentile(np.abs(y_fb), 95) + 1e-9)
    nsb = y_sb / (np.nanpercentile(np.abs(y_sb), 95) + 1e-9) if np.any(np.isfinite(y_sb)) else np.zeros_like(tB)
    axB.plot(tB, nsb * 0.5 + 1.0, color=COLOR_STIM_PULSE, linewidth=1.4)
    axB.plot(tB, nfb * 0.9 - 0.8, color=COLOR_GEVI, linewidth=2.8)
    axB.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axB.set_xlim(-0.3, 0.3)
    axB.set_ylim(-2.0, 2.0)
    axB.set_yticks([1.0, -0.8])
    axB.set_yticklabels(["Stim", "Fiber Vm"], fontsize=fs_tick)
    axB.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")

    # Panel C
    axC = fig.add_subplot(gs[1, 0])
    if np.any(np.isfinite(fiber_mu)):
        axC.plot(t_grid, fiber_mu, color=COLOR_GEVI, linewidth=2.2)
        if np.any(np.isfinite(fiber_sem)):
            axC.fill_between(t_grid, fiber_mu - fiber_sem, fiber_mu + fiber_sem, color=COLOR_GEVI, alpha=0.28, linewidth=0)
    axC.axvline(0.0, color="k", linestyle="--", linewidth=1.2, alpha=0.8)
    axC.axvline(stim_sec, color="k", linestyle=":", linewidth=1.2, alpha=0.8)
    axC.set_xlim(-3.0, 4.0)
    axC.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    axC.set_ylabel("Population Vm (ΔF/F %)", fontsize=fs_label, fontweight="bold")

    # Panel D
    axD = fig.add_subplot(gs[1, 1])
    if spec_avg is not None:
        ff, tt, ss = spec_avg
        ss_disp = uniform_filter1d(uniform_filter1d(ss, size=3, axis=0, mode="nearest"), size=5, axis=1, mode="nearest")
        vmax = np.nanpercentile(np.abs(ss_disp[np.isfinite(ss_disp)]), 95) if np.any(np.isfinite(ss_disp)) else 2.0
        im = axD.pcolormesh(tt, ff, ss_disp, cmap=CMAP_PARULA_LIKE, shading="auto", vmin=-vmax, vmax=vmax)
        axD.set_ylim(1, 100)
        cbar = fig.colorbar(im, ax=axD, pad=0.02)
        cbar.set_label("Power (dB)", fontsize=fs_tick)
        cbar.ax.tick_params(labelsize=max(fs_tick - 2, 1))
    else:
        axD.text(0.5, 0.5, "No spectrogram data", transform=axD.transAxes, ha="center", va="center", fontsize=fs_label)
    axD.axvline(0.0, color="w", linestyle="--", linewidth=1.2, alpha=0.9)
    axD.axvline(stim_sec, color="w", linestyle=":", linewidth=1.2, alpha=0.9)
    axD.set_xlim(-3.0, 4.0)
    axD.set_xlabel("Time from stim onset (s)", fontsize=fs_label, fontweight="bold")
    axD.set_ylabel("Frequency (Hz)", fontsize=fs_label, fontweight="bold")

    # Panel E/F
    axE = fig.add_subplot(gs[2, 0])
    _plot_three_period_violin_with_stats(
        axE,
        period_fiber,
        ylabel="Population Vm change (ΔF/F %)",
        title=None,
    )
    axE.set_ylabel("Population Vm change (ΔF/F %)", fontsize=fs_label, fontweight="bold")
    axF = fig.add_subplot(gs[2, 1])
    _plot_three_period_violin_with_stats(
        axF,
        period_band_db,
        ylabel="Relative 40Hz band power (dB)",
        title=None,
    )
    axF.set_ylabel("Relative 40Hz band power (dB)", fontsize=fs_label, fontweight="bold")

    for ax in [axA, axB, axC, axD, axE, axF]:
        ax.tick_params(labelsize=fs_tick, width=TICK_WIDTH, length=TICK_LENGTH)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(AXIS_LINEWIDTH)
        ax.spines["bottom"].set_linewidth(AXIS_LINEWIDTH)
    for ax in [axE, axF]:
        for txt in ax.texts:
            txt.set_fontsize(txt.get_fontsize() + 8)
    fig.subplots_adjust(left=0.06, right=0.965, top=0.985, bottom=0.055)
    return fig


def run_striatum_dbs_composite(trial1_mat_path, num_trials=10, invert_fiber=True):
    print("=" * 70)
    print("STRIATUM DBS COMPOSITE FIGURE")
    print("=" * 70)
    fig = create_striatum_dbs_composite_figure(
        trial1_mat_path=trial1_mat_path,
        num_trials=num_trials,
        invert_fiber=invert_fiber,
    )
    out_base = OUTPUT_DIR / "StriatumDBS"
    fig.savefig(str(out_base) + ".pdf", dpi=DPI, bbox_inches="tight")
    fig.savefig(str(out_base) + ".png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_base}.pdf/.png")

# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 70)
    print("STIMULATION ANALYSIS - MULTIPLE FIGURES")
    print("=" * 70)
    
    stim_onset_time = PRE_STIM_DURATION
    
    # Load data
    print("\n[1/4] Loading R1 (40Hz) raw trials...")
    raw_trials_r1 = load_all_raw_trials(SESSION_R1)
    
    print("\n[2/4] Loading R2 (135Hz) raw trials...")
    raw_trials_r2 = load_all_raw_trials(SESSION_R2)
    
    print("\n[3/4] Loading R1 (40Hz) spectral results...")
    spectral_r1 = load_all_spectral_results(SESSION_R1)
    
    print("\n[4/4] Loading R2 (135Hz) spectral results...")
    spectral_r2 = load_all_spectral_results(SESSION_R2)
    
    # Get representative trials for heatmaps
    raw_trial_r1 = next((t for t in raw_trials_r1 if t is not None), None)
    raw_trial_r2 = next((t for t in raw_trials_r2 if t is not None), None)
    
    # Create figures
    print("\n" + "=" * 70)
    print("CREATING FIGURES")
    print("=" * 70)
    
    figures = {}
    
    print(f"\n  Figure 1: R1 ({LABEL_R1}) Traces...")
    figures['traces_r1'] = create_traces_figure(SESSION_R1, LABEL_R1, raw_trials_r1, stim_onset_time)

    print(f"  Figure 2: R1 ({LABEL_R1}) Stim-onset zoom (−{STIM_ONSET_ZOOM_PRE_S} to +{STIM_ONSET_ZOOM_POST_S} s)...")
    figures['stim_zoom_r1'] = create_traces_stim_onset_zoom_figure(
        SESSION_R1, LABEL_R1, raw_trials_r1, stim_onset_time
    )

    print(f"  Figure 3: R2 ({LABEL_R2}) Traces...")
    figures['traces_r2'] = create_traces_figure(SESSION_R2, LABEL_R2, raw_trials_r2, stim_onset_time)

    print(f"  Figure 4: R2 ({LABEL_R2}) Stim-onset zoom (−{STIM_ONSET_ZOOM_PRE_S} to +{STIM_ONSET_ZOOM_POST_S} s)...")
    figures['stim_zoom_r2'] = create_traces_stim_onset_zoom_figure(
        SESSION_R2, LABEL_R2, raw_trials_r2, stim_onset_time
    )

    print(f"  Figure 5: R1 ({LABEL_R1}) Spectral Heatmaps...")
    figures['heatmaps_r1'] = create_heatmaps_figure(SESSION_R1, LABEL_R1, spectral_r1, raw_trial_r1, stim_onset_time)
    
    print(f"  Figure 6: R2 ({LABEL_R2}) Spectral Heatmaps...")
    figures['heatmaps_r2'] = create_heatmaps_figure(SESSION_R2, LABEL_R2, spectral_r2, raw_trial_r2, stim_onset_time)
    
    print(f"  Figure 7: R1 ({LABEL_R1}) Period Violins...")
    figures['violins_r1'] = create_period_violins_figure(SESSION_R1, LABEL_R1, spectral_r1)
    
    print(f"  Figure 8: R2 ({LABEL_R2}) Period Violins...")
    figures['violins_r2'] = create_period_violins_figure(SESSION_R2, LABEL_R2, spectral_r2)
    
    print(f"  Figure 9: {LABEL_R1} vs {LABEL_R2} Comparison...")
    figures['comparison'] = create_comparison_figure(spectral_r1, spectral_r2)
    
    print(f"  Figure 10: R1 ({LABEL_R1}) Trial-by-Trial Heatmaps...")
    figures['trials_r1'] = create_trials_heatmap_figure(SESSION_R1, LABEL_R1, raw_trials_r1, stim_onset_time, stim_freq=FREQ_R1, spectral_results=spectral_r1)
    
    print(f"  Figure 11: R2 ({LABEL_R2}) Trial-by-Trial Heatmaps...")
    figures['trials_r2'] = create_trials_heatmap_figure(SESSION_R2, LABEL_R2, raw_trials_r2, stim_onset_time, stim_freq=FREQ_R2, spectral_results=spectral_r2)

    print(f"  Figure 12: Trial-averaged LFP/GEVI overlaid ({LABEL_R1} vs {LABEL_R2})...")
    figures["avg_overlay"] = create_overlaid_avg_traces_figure(
        raw_trials_r1, raw_trials_r2, stim_onset_time, LABEL_R1, LABEL_R2
    )

    print(f"  Figure 13: R1 ({LABEL_R1}) Fiber Spectrogram Only...")
    figures['fiber_spec_r1'] = create_fiber_spectrogram_only_figure(
        SESSION_R1, LABEL_R1, spectral_r1, stim_onset_time
    )

    print(f"  Figure 14: R2 ({LABEL_R2}) Fiber Spectrogram Only...")
    figures['fiber_spec_r2'] = create_fiber_spectrogram_only_figure(
        SESSION_R2, LABEL_R2, spectral_r2, stim_onset_time
    )

    # Save figures
    print("\n" + "=" * 70)
    print("SAVING FIGURES")
    print("=" * 70)
    
    # Create sanitized labels for filenames (remove spaces and special chars)
    label_r1_safe = LABEL_R1.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')
    label_r2_safe = LABEL_R2.replace(' ', '_').replace('(', '').replace(')', '').replace('-', '')
    
    output_names = {
        'traces_r1': f'{MOUSE_ID}_01_traces_{label_r1_safe}',
        'stim_zoom_r1': f'{MOUSE_ID}_02_stim_onset_zoom_{label_r1_safe}',
        'traces_r2': f'{MOUSE_ID}_03_traces_{label_r2_safe}',
        'stim_zoom_r2': f'{MOUSE_ID}_04_stim_onset_zoom_{label_r2_safe}',
        'heatmaps_r1': f'{MOUSE_ID}_05_heatmaps_{label_r1_safe}',
        'heatmaps_r2': f'{MOUSE_ID}_06_heatmaps_{label_r2_safe}',
        'violins_r1': f'{MOUSE_ID}_07_violins_{label_r1_safe}',
        'violins_r2': f'{MOUSE_ID}_08_violins_{label_r2_safe}',
        'comparison': f'{MOUSE_ID}_09_comparison_{label_r1_safe}_vs_{label_r2_safe}',
        'trials_r1': f'{MOUSE_ID}_10_trials_{label_r1_safe}',
        'trials_r2': f'{MOUSE_ID}_11_trials_{label_r2_safe}',
        'avg_overlay': f'{MOUSE_ID}_12_trial_avg_overlay_{label_r1_safe}_vs_{label_r2_safe}',
        'fiber_spec_r1': f'{MOUSE_ID}_13_fiber_spectrogram_only_{label_r1_safe}',
        'fiber_spec_r2': f'{MOUSE_ID}_14_fiber_spectrogram_only_{label_r2_safe}',
    }
    
    for key, fig in figures.items():
        output_name = output_names[key]
        output_path = OUTPUT_DIR / output_name
        
        fig.savefig(str(output_path) + '.pdf', dpi=DPI, bbox_inches='tight')
        fig.savefig(str(output_path) + '.png', dpi=DPI, bbox_inches='tight')
        print(f"  Saved: {output_name}.pdf/.png")
        plt.close(fig)
    
    print("\n" + "=" * 70)
    print("DONE! Created 14 figures.")
    print("=" * 70)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Stimulation analysis and Striatum DBS composite plotting.")
    parser.add_argument(
        "--mode",
        choices=["default", "striatum-dbs"],
        default="default",
        help="Run full default stimulation-analysis suite, or the Striatum DBS composite figure.",
    )
    parser.add_argument(
        "--trial1-mat",
        type=str,
        default=None,
        help="Absolute path to Trial1 *_FiberPhotometry_Analysis.mat for Striatum DBS composite mode.",
    )
    parser.add_argument(
        "--num-trials",
        type=int,
        default=10,
        help="Number of trial mats to include from session folder in Striatum DBS mode.",
    )
    args = parser.parse_args()

    if args.mode == "striatum-dbs":
        if not args.trial1_mat:
            raise ValueError("--trial1-mat is required when --mode striatum-dbs")
        run_striatum_dbs_composite(args.trial1_mat, num_trials=args.num_trials)
    else:
        main()
