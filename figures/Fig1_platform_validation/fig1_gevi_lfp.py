"""
Plot Figure 1 for GEVI–LFP–motion recordings from MATLAB FiberPhotometryAnalysis .mat files.

This script:
- Loads a single-trial FiberPhotometryAnalysis .mat file
- Extracts:
    - Time vector (seconds)
    - ΔF/F GEVI trace (per-fiber; default: first fiber)
    - Aligned LFP trace
    - Motion trace
- Creates a 5-panel figure with shared time axis:
    1. Raw ΔF/F GEVI
    2. Raw LFP
    3. Motion
    4. Theta-band (5–9 Hz) GEVI & LFP (z-scored within panel)
    5. Beta-band (15–30 Hz) GEVI & LFP (z-scored within panel)
- Saves a vector PDF into a dedicated figure folder inside the project root.

Dependencies:
- numpy
- matplotlib
- scipy (scipy.io, scipy.signal, scipy.stats)

Typical usage:

1) Simple (recommended for daily use) – configure inside this script and just run:

    - Set DEFAULT_MAT_PATH, DEFAULT_FIBER_INDEX, DEFAULT_FIGURE_NAME
      in the USER CONFIGURATION section below.
    - Then run from the terminal:

        python fig1_gevi_lfp.py

2) Advanced – override from the command line:

        python fig1_gevi_lfp.py \\
            --mat-path "PATH/TO/FiberPhotometry_Analysis.mat" \\
            --fiber-index 0 \\
            --figure-name "MyCustomFigureName"

If --mat-path is omitted, the DEFAULT_MAT_PATH from the USER CONFIGURATION
section (or the built-in example) is used.
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy import signal
from scipy.stats import zscore
import h5py

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# =============================================================================
# USER CONFIGURATION (EDIT THESE FOR YOUR OWN DATASET)
# =============================================================================

# Default MAT file path -- EDIT THIS to point at your own recording.
DEFAULT_MAT_PATH = str(
    _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal01" / "Fiber_Voltage_Processed"
    / "01_01_25-R1" / "Trial1_fov1_HPmPFCdual_baselineRecording_60sec_1"
    / "Animal01-01_01_25-R1_Trial1_FiberPhotometry_Analysis.mat"
)


# By default, always use the first fiber (index 0) unless overridden.
DEFAULT_FIBER_INDEX = 0

# Default figure name (without extension); used when --figure-name is not given.
DEFAULT_FIGURE_NAME = "Fig_head_fixed_traces_Animal01_01-01-25-R1_rest"

# Optional: custom output directory. If None, a "python_figures" folder will be
# created inside the project root.
DEFAULT_OUTPUT_DIR = str(PROJECT_ROOT / "Figures" / "Traces_figures")
# Optional: time window [t_start, t_end] in seconds for plotting.
# If set to None, the full recording is used.
DEFAULT_TIME_WINDOW = (6.0, 6.5)   # example: show 6.0-6.5 s

# LFP unit conversion: Set to True if the MATLAB file stores LFP in Volts (needs *1e6).
# Set to False if LFP is already stored in µV (no conversion needed).
LFP_IN_VOLTS = False  # <-- CHANGE THIS based on your data!

# Scale bar configuration (in physical units)
SCALEBAR_T_SECONDS = 0.5  # horizontal scale bar duration

# Raw trace vertical scale bars
SCALEBAR_GEVI_PERCENT = 1.0      # 1% ΔF/F for raw GEVI
SCALEBAR_LFP_UV = 500.0          # 500 µV for raw LFP

# Band-limited trace scale bars (slightly smaller)
SCALEBAR_GEVI_BAND_PERCENT = 0.5  # 0.5% ΔF/F for filtered GEVI
SCALEBAR_LFP_BAND_UV = 200.0      # 200 µV for filtered LFP

# =============================================================================
# MOTION/RUNNING SPEED CONVERSION PARAMETERS
# =============================================================================
# The MATLAB code stores motion as: running_velocity_smooth
# Processing in MATLAB:
#   1. running_velocity = diff(running_wheel_raw) > 1  (rising edges at 30 kHz)
#   2. fastsmooth(running_velocity, 30, 1, 1)  → edge density (30-sample = 1 ms window)
#   3. Aligned to imaging rate via camera triggers
#   4. running_velocity_aligned × 1000, then fastsmooth(..., 40, 1, 1)
#
# To convert back to cm/s:
#   - The stored value represents (edges_per_ms × 33.33) approximately
#   - speed (cm/s) = stored_value × (EPHYS_FS / 1000) × (circumference / counts_per_rev)

# Wheel parameters (Yumo E6B2 encoder)
WHEEL_DIAMETER_CM = 19.0                           # Wheel diameter in cm
WHEEL_CIRCUMFERENCE_CM = np.pi * WHEEL_DIAMETER_CM # ≈ 59.69 cm
ENCODER_COUNTS_PER_REV = 1024                      # Yumo E6B2: 1024 pulses/revolution
EPHYS_SAMPLING_RATE = 30000                        # Open Ephys sampling rate (Hz)

# Distance traveled per detected rising edge
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV  # ≈ 0.0583 cm

# Conversion factor from MATLAB's running_velocity_smooth to cm/s
# Based on MATLAB processing: value = (edges_per_30samples / 30) × 1000 = edges_per_ms × 33.33
# To get edges/s: value × (1000 / 33.33) = value × 30
# To get cm/s: edges/s × distance_per_edge = value × 30 × distance_per_edge
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM  # ≈ 1.749

# Optional smoothing for the converted speed trace (in samples at imaging rate)
# Set to 0 or None to disable additional smoothing
# Note: MATLAB already applies ~80 ms smoothing
# Adjust this value to taste:
#   - 0 = no extra smoothing (spiky)
#   - 20 = ~40 ms additional smoothing (light)
#   - 50 = ~100 ms additional smoothing (heavy)
MOTION_SMOOTH_SAMPLES = 10  # Light smoothing (~40 ms at ~500 Hz)

# Motion scale bar (now in cm/s)
SCALEBAR_MOTION_CM_S = 15.0  # 15 cm/s scale bar

# Publication-ready font sizes
FONT_SIZE_SCALEBAR = 11      # Scale bar labels
FONT_SIZE_TITLE = 14         # Figure title
FONT_SIZE_AXIS_LABEL = 12    # Axis labels (if any)

# Scale bar line thickness
SCALEBAR_LINEWIDTH = 3.5     # Thicker for publication visibility


def matlab_struct_to_dict(matobj):
    """
    Recursively convert MATLAB structs (from scipy.io.loadmat) to nested Python dicts.

    This handles objects of type numpy.void (MATLAB struct) and nested fields.
    """
    if isinstance(matobj, np.ndarray) and matobj.dtype == np.object_ and matobj.size == 1:
        # Sometimes structs come as 1x1 object arrays
        matobj = matobj.item()

    if isinstance(matobj, np.void):
        out = {}
        for field_name in matobj.dtype.names:
            out[field_name] = matlab_struct_to_dict(matobj[field_name])
        return out

    if isinstance(matobj, np.ndarray) and matobj.dtype == np.object_:
        # Object array: convert each element
        return [matlab_struct_to_dict(el) for el in matobj]

    return matobj


def bandpass_filter_trace(x, fs, low_hz, high_hz, order=3):
    """
    Apply a zero-phase bandpass filter to a 1D trace using scipy.signal.

    Parameters
    ----------
    x : array-like
        Input signal (1D).
    fs : float
        Sampling rate in Hz.
    low_hz : float
        Low cutoff frequency in Hz.
    high_hz : float
        High cutoff frequency in Hz.
    order : int
        Filter order (Butterworth).

    Returns
    -------
    y : ndarray
        Bandpass-filtered signal.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 1:
        raise ValueError("bandpass_filter_trace expects a 1D array")

    nyq = 0.5 * fs
    low = low_hz / nyq
    high = high_hz / nyq
    if low <= 0 or high >= 1 or low >= high:
        raise ValueError(
            f"Invalid band ({low_hz}, {high_hz}) for sampling rate {fs} Hz"
        )

    b, a = signal.butter(order, [low, high], btype="bandpass")

    # filtfilt requires the input length to be greater than padlen (~3 * max(len(a), len(b)))
    padlen = 3 * max(len(a), len(b))
    if x.size <= padlen:
        # For very short signals, skip bandpass filtering and return a demeaned version
        # so that downstream z-scoring and plotting still work without error.
        # This is only for visualization; it avoids filtfilt edge errors.
        # (We keep the scale comparable by subtracting the mean only.)
        return x - np.nanmean(x)

    y = signal.filtfilt(b, a, x, method="pad")
    return y


def convert_motion_to_speed(motion_raw, smooth_samples=None):
    """
    Convert MATLAB's running_velocity_smooth to true running speed in cm/s.

    The MATLAB code processes the wheel encoder signal as follows:
    1. running_velocity = diff(running_wheel_raw) > 1  (detects rising edges at 30 kHz)
    2. fastsmooth(running_velocity, 30, 1, 1)  → 30-sample moving average (1 ms window)
    3. Aligned to imaging rate via camera triggers
    4. running_velocity_smooth = fastsmooth(running_velocity_aligned × 1000, 40, 1, 1)

    The stored value is approximately: (edges_per_ms) × 33.33
    To convert to cm/s:
        speed = stored_value × MOTION_TO_CM_PER_S

    Parameters
    ----------
    motion_raw : array-like
        The running_velocity_smooth trace from MATLAB (arbitrary units).
    smooth_samples : int or None
        Optional: apply additional moving average smoothing with this window size.
        Set to None or 0 to skip.

    Returns
    -------
    speed_cm_s : ndarray
        Running speed in cm/s.
    """
    motion = np.asarray(motion_raw, dtype=float)
    
    # Convert to cm/s using the derived conversion factor
    speed_cm_s = motion * MOTION_TO_CM_PER_S
    
    # Optional smoothing
    if smooth_samples is not None and smooth_samples > 1:
        # Simple moving average using convolution
        kernel = np.ones(smooth_samples) / smooth_samples
        # Use 'same' mode and handle edges
        speed_cm_s = np.convolve(speed_cm_s, kernel, mode='same')
    
    return speed_cm_s


def _load_fiberphotometry_hdf5(mat_path, fiber_index=0):
    """
    Load a MATLAB v7.3 (HDF5) FiberPhotometryAnalysis file using h5py.

    This is used as a fallback when scipy.io.loadmat cannot read the file.
    """
    mat_path = Path(mat_path)
    if not mat_path.is_file():
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    with h5py.File(mat_path, "r") as f:
        if "FiberPhotometryAnalysis" not in f:
            raise KeyError(
                "FiberPhotometryAnalysis group not found in HDF5 MAT file. "
                "Please ensure you are using the MATLAB output from the fiber_processing_* scripts."
            )

        root = f["FiberPhotometryAnalysis"]

        # --- Time vector & sampling rate ---
        try:
            t = np.array(root["time"]["time_vector_seconds"][()]).reshape(-1)
        except KeyError as exc:
            raise KeyError(
                "Could not find time/time_vector_seconds in FiberPhotometryAnalysis (HDF5)"
            ) from exc

        if "sampling_rate" in root["time"]:
            # Extract scalar safely to avoid deprecation warnings
            sr = np.array(root["time"]["sampling_rate"][()])
            fs_gevi = float(sr.item() if sr.size == 1 else sr.ravel()[0])
        else:
            if t.size < 2:
                raise ValueError("Time vector too short to infer sampling rate (HDF5)")
            dt = np.median(np.diff(t))
            fs_gevi = 1.0 / dt

        # --- GEVI ΔF/F trace ---
        try:
            gevi_all = np.array(root["signals"]["final_processed_traces"][()])
        except KeyError as exc:
            raise KeyError(
                "Could not find signals/final_processed_traces (ΔF/F) in FiberPhotometryAnalysis (HDF5)"
            ) from exc

        # Normalize shape/orientation.
        # We expect [time x fibers]. Handle common cases like [1 x N] or [N x 1].
        if gevi_all.ndim == 1:
            gevi_df = gevi_all.reshape(-1)
        elif gevi_all.ndim == 2 and 1 in gevi_all.shape:
            # Single fiber stored as row/column vector
            gevi_df = gevi_all.reshape(-1)
        else:
            n_time = t.size
            # If first dimension is not time but second is, transpose
            if gevi_all.shape[0] != n_time and gevi_all.shape[1] == n_time:
                gevi_all = gevi_all.T
            if gevi_all.shape[0] != n_time:
                raise ValueError(
                    f"Unexpected ΔF/F array shape {gevi_all.shape} (HDF5); "
                    f"could not align with time vector of length {n_time}."
                )
            # Now gevi_all is [time x fibers]
            if fiber_index < 0 or fiber_index >= gevi_all.shape[1]:
                raise IndexError(
                    f"fiber_index {fiber_index} out of range for ΔF/F traces with shape {gevi_all.shape} (HDF5)"
                )
            gevi_df = gevi_all[:, fiber_index]

        # --- LFP trace ---
        lfp = None
        if "ephys" in root:
            ephys_grp = root["ephys"]
            if "lfp_raw_aligned_HP" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_HP"][()]).reshape(-1)
            elif "lfp_raw_aligned_mPFC" in ephys_grp:
                lfp = np.array(ephys_grp["lfp_raw_aligned_mPFC"][()]).reshape(-1)

        if lfp is None:
            raise KeyError(
                "Could not find an aligned LFP trace in FiberPhotometryAnalysis/ephys (HDF5) "
                "(expected 'lfp_raw_aligned_HP' or 'lfp_raw_aligned_mPFC')."
            )

        # --- Motion trace ---
        motion = None
        if "ephys" in root:
            ephys_grp = root["ephys"]
            if "running_velocity_smooth" in ephys_grp:
                motion = np.array(ephys_grp["running_velocity_smooth"][()]).reshape(-1)
            elif "running_velocity" in ephys_grp:
                motion = np.array(ephys_grp["running_velocity"][()]).reshape(-1)

        if motion is None:
            raise KeyError(
                "Could not find motion trace in FiberPhotometryAnalysis/ephys (HDF5) "
                "(expected 'running_velocity_smooth' or 'running_velocity')."
            )

    # --- Sanity: align lengths ---
    n = min(len(t), len(gevi_df), len(lfp), len(motion))
    if n == 0:
        raise ValueError("One or more signals are empty after loading (HDF5).")

    t = t[:n]
    gevi_df = gevi_df[:n]
    lfp = lfp[:n]
    motion = motion[:n]

    fs_lfp = fs_gevi

    return {
        "t": t,
        "gevi_df": gevi_df,
        "lfp": lfp,
        "motion": motion,
        "fs_gevi": fs_gevi,
        "fs_lfp": fs_lfp,
    }


def load_fiberphotometry_mat(mat_path, fiber_index=0):
    """
    Load a FiberPhotometryAnalysis .mat file and extract relevant signals.

    Parameters
    ----------
    mat_path : str or Path
        Path to the .mat file containing FiberPhotometryAnalysis.
    fiber_index : int
        Zero-based index of the fiber to use (default: 0, i.e., first fiber).

    Returns
    -------
    data : dict
        Dictionary with keys:
            't'          : time vector in seconds (1D)
            'gevi_df'    : ΔF/F GEVI trace for the selected fiber (1D)
            'lfp'        : aligned LFP trace (1D)
            'motion'     : motion / running-velocity trace (1D)
            'fs_gevi'    : sampling rate of GEVI (float)
            'fs_lfp'     : sampling rate of LFP (float; equals fs_gevi after alignment)
    """
    mat_path = Path(mat_path)
    if not mat_path.is_file():
        raise FileNotFoundError(f"MAT file not found: {mat_path}")

    # First try the standard MATLAB reader (for v7.2 and below).
    # If the file is v7.3 (HDF5-based), fall back to h5py loader.
    try:
        mat = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    except NotImplementedError:
        # v7.3 file – use HDF5 reader instead
        return _load_fiberphotometry_hdf5(mat_path, fiber_index=fiber_index)
    if "FiberPhotometryAnalysis" not in mat:
        raise KeyError(
            "FiberPhotometryAnalysis struct not found in MAT file. "
            "Please ensure you are using the MATLAB output from the fiber_processing_* scripts."
        )

    fp_struct = mat["FiberPhotometryAnalysis"]
    fp = matlab_struct_to_dict(fp_struct)

    # --- Time vector & sampling rate ---
    try:
        t = np.asarray(fp["time"]["time_vector_seconds"]).reshape(-1)
    except KeyError as exc:
        raise KeyError("Could not find time.time_vector_seconds in FiberPhotometryAnalysis") from exc

    try:
        fs_gevi = float(fp["time"]["sampling_rate"])
    except KeyError:
        # Fallback: infer from time vector
        if t.size < 2:
            raise ValueError("Time vector too short to infer sampling rate")
        dt = np.median(np.diff(t))
        fs_gevi = 1.0 / dt

    # --- GEVI ΔF/F trace ---
    try:
        gevi_all = np.asarray(fp["signals"]["final_processed_traces"])
    except KeyError as exc:
        raise KeyError("Could not find signals.final_processed_traces (ΔF/F) in FiberPhotometryAnalysis") from exc

    # Normalize shape/orientation for non-HDF case as well.
    if gevi_all.ndim == 1:
        gevi_df = gevi_all.reshape(-1)
    elif gevi_all.ndim == 2 and 1 in gevi_all.shape:
        # Single fiber as row/column
        gevi_df = gevi_all.reshape(-1)
    else:
        n_time = t.size
        if gevi_all.shape[0] != n_time and gevi_all.shape[1] == n_time:
            gevi_all = gevi_all.T
        if gevi_all.shape[0] != n_time:
            raise ValueError(
                f"Unexpected ΔF/F array shape {gevi_all.shape}; "
                f"could not align with time vector of length {n_time}."
            )
        if fiber_index < 0 or fiber_index >= gevi_all.shape[1]:
            raise IndexError(
                f"fiber_index {fiber_index} out of range for ΔF/F traces with shape {gevi_all.shape}"
            )
        gevi_df = gevi_all[:, fiber_index]

    # --- LFP trace ---
    lfp = None
    ephys = fp.get("ephys", {})
    if isinstance(ephys, dict):
        # Prefer HP LFP (Ch11) if present
        if "lfp_raw_aligned_HP" in ephys:
            lfp = np.asarray(ephys["lfp_raw_aligned_HP"]).reshape(-1)
        elif "lfp_raw_aligned_mPFC" in ephys:
            lfp = np.asarray(ephys["lfp_raw_aligned_mPFC"]).reshape(-1)

    if lfp is None:
        raise KeyError(
            "Could not find an aligned LFP trace in FiberPhotometryAnalysis.ephys "
            "(expected 'lfp_raw_aligned_HP' or 'lfp_raw_aligned_mPFC')."
        )

    # --- Motion trace ---
    motion = None
    if isinstance(ephys, dict):
        if "running_velocity_smooth" in ephys:
            motion = np.asarray(ephys["running_velocity_smooth"]).reshape(-1)
        elif "running_velocity" in ephys:
            motion = np.asarray(ephys["running_velocity"]).reshape(-1)

    if motion is None:
        raise KeyError(
            "Could not find motion trace in FiberPhotometryAnalysis.ephys "
            "(expected 'running_velocity_smooth' or 'running_velocity')."
        )

    # --- Sanity: align lengths ---
    n = min(len(t), len(gevi_df), len(lfp), len(motion))
    if n == 0:
        raise ValueError("One or more signals are empty after loading.")

    t = t[:n]
    gevi_df = gevi_df[:n]
    lfp = lfp[:n]
    motion = motion[:n]

    # After alignment, LFP sampling rate equals GEVI sampling rate
    fs_lfp = fs_gevi

    return {
        "t": t,
        "gevi_df": gevi_df,
        "lfp": lfp,
        "motion": motion,
        "fs_gevi": fs_gevi,
        "fs_lfp": fs_lfp,
    }


def create_figure_1(data, figure_name, output_dir):
    """
    Create Figure 1 as 5 vertically stacked traces (no boxed panels),
    each with its own vertical scale bar, plus a single horizontal scale bar
    at the bottom. Save in multiple vector/raster formats.

    For band-filtered rows (theta, beta), GEVI and LFP are plotted on separate
    y-axes (twinx) so both signals are visually comparable despite very different
    physical-unit amplitudes.

    Parameters
    ----------
    data : dict
        Output of load_fiberphotometry_mat().
    figure_name : str
        Base name for the figure file (without extension).
    output_dir : str or Path
        Directory where the PDF will be saved.
    """
    t_full = data["t"]
    gevi_df_full = data["gevi_df"]          # ΔF/F (fraction)
    lfp_full = data["lfp"]                  # Volts (from Open Ephys)
    motion_full = data["motion"]            # Arbitrary units
    fs = data["fs_gevi"]                    # GEVI & LFP are aligned

    # Apply optional time window
    if DEFAULT_TIME_WINDOW is not None:
        t_start, t_end = DEFAULT_TIME_WINDOW
        mask = (t_full >= t_start) & (t_full <= t_end)
        if not np.any(mask):
            raise ValueError(
                f"Time window {DEFAULT_TIME_WINDOW} does not overlap recording "
                f"(t in [{t_full.min():.3f}, {t_full.max():.3f}])."
            )
        t = t_full[mask]
        gevi_df = gevi_df_full[mask]
        lfp = lfp_full[mask]
        motion = motion_full[mask]
    else:
        t = t_full
        gevi_df = gevi_df_full
        lfp = lfp_full
        motion = motion_full

    # Convert to physical units for plotting
    gevi_percent = gevi_df * 100.0           # ΔF/F in percent
    
    # LFP conversion depends on source units
    if LFP_IN_VOLTS:
        lfp_uV = lfp * 1e6                   # volts -> microvolts
    else:
        lfp_uV = lfp                         # already in µV, no conversion
    
    # Debug: print data ranges to help diagnose scale bar issues
    print(f"[DEBUG] GEVI range: {np.nanmin(gevi_percent):.2f} to {np.nanmax(gevi_percent):.2f} %")
    print(f"[DEBUG] LFP range: {np.nanmin(lfp_uV):.2f} to {np.nanmax(lfp_uV):.2f} µV")
    print(f"[DEBUG] GEVI scale bar: {SCALEBAR_GEVI_PERCENT}% -> fraction of range: {SCALEBAR_GEVI_PERCENT / (np.nanmax(gevi_percent) - np.nanmin(gevi_percent)):.3f}")
    print(f"[DEBUG] LFP scale bar: {SCALEBAR_LFP_UV} µV -> fraction of range: {SCALEBAR_LFP_UV / (np.nanmax(lfp_uV) - np.nanmin(lfp_uV)):.3f}")

    # Band-limited versions (kept in physical units)
    gevi_theta = bandpass_filter_trace(gevi_percent, fs, 5.0, 9.0)
    lfp_theta = bandpass_filter_trace(lfp_uV, fs, 5.0, 9.0)
    gevi_beta = bandpass_filter_trace(gevi_percent, fs, 30.0, 60.0)
    lfp_beta = bandpass_filter_trace(lfp_uV, fs, 30.0, 60.0)

    # Convert motion to true running speed (cm/s)
    motion_cm_s = convert_motion_to_speed(motion, smooth_samples=MOTION_SMOOTH_SAMPLES)
    
    # Debug: print motion speed range
    print(f"[DEBUG] Motion speed range: {np.nanmin(motion_cm_s):.2f} to {np.nanmax(motion_cm_s):.2f} cm/s")
    print(f"[DEBUG] Motion scale bar: {SCALEBAR_MOTION_CM_S} cm/s -> fraction of range: {SCALEBAR_MOTION_CM_S / (np.nanmax(motion_cm_s) - np.nanmin(motion_cm_s) + 1e-9):.3f}")

    # Prepare figure
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    base_path = output_dir / figure_name

    COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])  # teal
    COLOR_LFP = np.array([0.35, 0.25, 0.45])              # purple-grey (more distinct from teal)
    COLOR_MOTION = np.array([0.993248, 0.7, 0.4])         # orange

    # Thinner line width for cleaner look
    LW_RAW = 0.8
    LW_BAND = 0.8

    # Transparency for overlaid filtered traces
    ALPHA_BAND = 0.85

    fig, axes = plt.subplots(
        5, 1, figsize=(8.0, 8.0), sharex=False,
        gridspec_kw={"hspace": 0.25}
    )

    t_min, t_max = t.min(), t.max()
    t_range = t_max - t_min

    # Padding fractions for x-axis limits (to leave space for scale bars)
    PAD_LEFT = 0.06   # fraction of t_range
    PAD_RIGHT = 0.06  # fraction of t_range

    def vbar_proportional(ax, x_frac, height_value, label, side="left", color="k"):
        """
        Draw a vertical scale bar using axis-relative positioning.
        
        Parameters
        ----------
        ax : matplotlib axis
        x_frac : float
            Position as fraction of axis width (0=left edge, 1=right edge).
            Use negative for left of plot, >1 for right of plot.
        height_value : float
            The scale bar height in DATA units (e.g., µV or %).
        label : str
            Text label for the scale bar.
        side : str
            'left' or 'right' - determines text placement relative to bar.
        color : color
            Color for both bar and text.
        """
        # Get current data limits
        y_min, y_max = ax.get_ylim()
        y_range = y_max - y_min
        
        # Calculate bar height as TRUE fraction of visible range
        # The bar accurately represents the scale value (no min/max clamping)
        bar_frac = height_value / y_range if y_range != 0 else 0.3
        
        # Center the bar vertically in the middle of the axis
        y_center_frac = 0.5
        y0_frac = y_center_frac - bar_frac / 2
        y1_frac = y_center_frac + bar_frac / 2
        
        # Convert to data coordinates
        y0 = y_min + y0_frac * y_range
        y1 = y_min + y1_frac * y_range
        
        # x position in data coordinates
        x_min, x_max = ax.get_xlim()
        x_range = x_max - x_min
        x_pos = x_min + x_frac * x_range
        
        # Draw the bar (use configurable linewidth for publication quality)
        ax.plot([x_pos, x_pos], [y0, y1], color=color, linewidth=SCALEBAR_LINEWIDTH, 
                clip_on=False, solid_capstyle='butt')
        
        # Text position
        text_offset_frac = 0.025  # Slightly more offset for thicker bars
        if side == "left":
            tx_frac = x_frac - text_offset_frac
            ha = "right"
        else:
            tx_frac = x_frac + text_offset_frac
            ha = "left"
        tx = x_min + tx_frac * x_range
        
        ax.text(
            tx,
            (y0 + y1) / 2.0,
            label,
            ha=ha,
            va="center",
            fontsize=FONT_SIZE_SCALEBAR,
            rotation=90,
            color=color,
            clip_on=False,
        )

    # Helper to remove all spines, ticks, labels
    def clean_axis(ax, remove_xlabel=True):
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.tick_params(axis="both", which="both", length=0, labelleft=False, labelbottom=False)
        ax.set_yticks([])
        ax.set_xticks([])
        if remove_xlabel:
            ax.set_xlabel("")
        ax.set_ylabel("")

    # Set x limits with padding for scale bars
    xlim_left = t_min - PAD_LEFT * t_range * 2.0
    xlim_right = t_max + PAD_RIGHT * t_range * 2.0

    # Fractional positions for scale bars (negative = left of plot, >1 = right of plot)
    SCALEBAR_X_LEFT = -0.02   # just outside left edge
    SCALEBAR_X_RIGHT = 1.02   # just outside right edge

    # ========== Row 1: raw GEVI ΔF/F (%) ==========
    ax_gevi = axes[0]
    ax_gevi.plot(t, gevi_percent, color=COLOR_GEVI, linewidth=LW_RAW)
    clean_axis(ax_gevi)
    ax_gevi.set_xlim(xlim_left, xlim_right)
    vbar_proportional(ax_gevi, SCALEBAR_X_LEFT, SCALEBAR_GEVI_PERCENT,
                      f"{SCALEBAR_GEVI_PERCENT:.1f}%", side="left", color=COLOR_GEVI)

    # ========== Row 2: raw LFP (µV) ==========
    ax_lfp = axes[1]
    ax_lfp.plot(t, lfp_uV, color=COLOR_LFP, linewidth=LW_RAW)
    clean_axis(ax_lfp)
    ax_lfp.set_xlim(xlim_left, xlim_right)
    vbar_proportional(ax_lfp, SCALEBAR_X_LEFT, SCALEBAR_LFP_UV,
                      f"{SCALEBAR_LFP_UV/1000.0:.1f} mV", side="left", color=COLOR_LFP)

    # ========== Row 3: Motion / Running Speed (cm/s) ==========
    # Converted from MATLAB's running_velocity_smooth using wheel parameters:
    #   - Wheel diameter: 19 cm → circumference ≈ 59.69 cm
    #   - Encoder: Yumo E6B2, 1024 counts/revolution
    #   - Distance per edge: ~0.0583 cm
    #   - Conversion factor: MOTION_TO_CM_PER_S ≈ 1.749
    ax_motion = axes[2]
    ax_motion.plot(t, motion_cm_s, color=COLOR_MOTION, linewidth=LW_RAW)
    clean_axis(ax_motion)
    ax_motion.set_xlim(xlim_left, xlim_right)
    vbar_proportional(ax_motion, SCALEBAR_X_LEFT, SCALEBAR_MOTION_CM_S, 
                      f"{int(SCALEBAR_MOTION_CM_S)} cm/s", side="left", color=COLOR_MOTION)

    # ========== Row 4: Theta band (5–9 Hz) - GEVI on left axis, LFP on right axis ==========
    ax_theta_gevi = axes[3]
    ax_theta_gevi.plot(t, gevi_theta, color=COLOR_GEVI, linewidth=LW_BAND, alpha=ALPHA_BAND)
    clean_axis(ax_theta_gevi)
    ax_theta_gevi.set_xlim(xlim_left, xlim_right)

    ax_theta_lfp = ax_theta_gevi.twinx()
    ax_theta_lfp.plot(t, lfp_theta, color=COLOR_LFP, linewidth=LW_BAND, alpha=ALPHA_BAND)
    clean_axis(ax_theta_lfp)
    ax_theta_lfp.set_xlim(xlim_left, xlim_right)

    # Scale bars for theta (left: GEVI, right: LFP)
    vbar_proportional(ax_theta_gevi, SCALEBAR_X_LEFT, SCALEBAR_GEVI_BAND_PERCENT,
                      f"{SCALEBAR_GEVI_BAND_PERCENT:.1f}% ΔF/F", side="left", color=COLOR_GEVI)
    vbar_proportional(ax_theta_lfp, SCALEBAR_X_RIGHT, SCALEBAR_LFP_BAND_UV,
                      f"{int(SCALEBAR_LFP_BAND_UV)} µV", side="right", color=COLOR_LFP)

    # ========== Row 5: Beta band (15–30 Hz) - GEVI on left axis, LFP on right axis ==========
    ax_beta_gevi = axes[4]
    ax_beta_gevi.plot(t, gevi_beta, color=COLOR_GEVI, linewidth=LW_BAND, alpha=ALPHA_BAND)
    clean_axis(ax_beta_gevi, remove_xlabel=False)
    ax_beta_gevi.set_xlim(xlim_left, xlim_right)

    ax_beta_lfp = ax_beta_gevi.twinx()
    ax_beta_lfp.plot(t, lfp_beta, color=COLOR_LFP, linewidth=LW_BAND, alpha=ALPHA_BAND)
    clean_axis(ax_beta_lfp)
    ax_beta_lfp.set_xlim(xlim_left, xlim_right)

    # Scale bars for beta (left: GEVI, right: LFP)
    vbar_proportional(ax_beta_gevi, SCALEBAR_X_LEFT, SCALEBAR_GEVI_BAND_PERCENT,
                      f"{SCALEBAR_GEVI_BAND_PERCENT:.1f}% ΔF/F", side="left", color=COLOR_GEVI)
    vbar_proportional(ax_beta_lfp, SCALEBAR_X_RIGHT, SCALEBAR_LFP_BAND_UV,
                      f"{int(SCALEBAR_LFP_BAND_UV)} µV", side="right", color=COLOR_LFP)

    # ========== Horizontal time scale bar (as separate annotation below figure) ==========
    # Use figure coordinates to place the horizontal scale bar reliably
    # First, get the position of the bottom axis in figure coordinates
    fig.canvas.draw()  # needed to get accurate positions
    bbox = ax_beta_gevi.get_position()
    
    # Scale bar in figure coordinates (x as fraction of figure width)
    # Place it at the right side of the bottom panel, below
    scalebar_y = bbox.y0 - 0.06  # below the bottom axis
    scalebar_x_end = bbox.x1 - 0.02
    # Calculate width in figure coords based on time
    time_to_fig = (bbox.x1 - bbox.x0) / t_range
    scalebar_width = SCALEBAR_T_SECONDS * time_to_fig
    scalebar_x_start = scalebar_x_end - scalebar_width
    
    # Draw using figure coordinates (thicker line for publication)
    from matplotlib.lines import Line2D
    line = Line2D([scalebar_x_start, scalebar_x_end], [scalebar_y, scalebar_y],
                  transform=fig.transFigure, color='k', linewidth=SCALEBAR_LINEWIDTH)
    fig.add_artist(line)
    
    # Text below the scale bar (publication font size)
    fig.text((scalebar_x_start + scalebar_x_end) / 2.0, scalebar_y - 0.025,
             f"{SCALEBAR_T_SECONDS:.1f} s", ha="center", va="top", fontsize=FONT_SIZE_SCALEBAR)

    # Title (publication font size)
    fig.suptitle("Figure 1: GEVI, LFP, Motion and Band-Limited Activity", 
                 fontsize=FONT_SIZE_TITLE, y=0.97)

    # Use more bottom margin to accommodate scale bar
    fig.tight_layout(rect=[0.10, 0.10, 0.90, 0.95])

    # Save in multiple formats: PNG, PDF, SVG
    for ext in ("png", "pdf", "svg"):
        out_path = base_path.with_suffix(f".{ext}")
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(f"Saved Figure 1 to: {out_path}")

    plt.close(fig)


def get_default_mat_path():
    """
    Default MAT file path.

    Priority:
    1) If DEFAULT_MAT_PATH is set in the USER CONFIGURATION section, use that.
    2) Otherwise, fall back to a generic example path under DATA_ROOT -- edit
       DEFAULT_MAT_PATH above to point at a real recording of yours.
    """
    if DEFAULT_MAT_PATH:
        return DEFAULT_MAT_PATH

    return str(
        _LAB_DATA_ROOT / "FiberVoltageImaging" / "Animal01" / "Fiber_Voltage_Processed"
        / "01_01_25-R1" / "Animal01-01_01_25-R1_FiberPhotometry_Analysis.mat"
    )


def get_default_output_dir():
    """
    Determine a dedicated figure folder inside the project root, unless a
    custom DEFAULT_OUTPUT_DIR is specified in USER CONFIGURATION.
    """
    if DEFAULT_OUTPUT_DIR:
        return Path(DEFAULT_OUTPUT_DIR)

    return PROJECT_ROOT / "python_figures"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate Figure 1 (GEVI, LFP, motion) from a FiberPhotometryAnalysis .mat file."
    )
    parser.add_argument(
        "--mat-path",
        type=str,
        default=None,  # if None, we use get_default_mat_path()
        help="Path to FiberPhotometry_Analysis.mat (default: from USER CONFIG / built-in example).",
    )
    parser.add_argument(
        "--fiber-index",
        type=int,
        default=None,  # if None, we use DEFAULT_FIBER_INDEX
        help="Zero-based fiber index to plot (default: from USER CONFIG; usually 0 for first fiber).",
    )
    parser.add_argument(
        "--figure-name",
        type=str,
        default=None,  # if None, we use DEFAULT_FIGURE_NAME
        help="Base name for the output PDF figure (without extension; default from USER CONFIG).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=None,  # if None, we use get_default_output_dir()
        help="Directory to save the figure (default: 'python_figures' inside the project root, or DEFAULT_OUTPUT_DIR).",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Resolve configuration with priority:
    # 1) Command-line arguments (if provided)
    # 2) USER CONFIGURATION defaults (for fiber index, figure name, output dir)
    # 3) Built-in example path (for MAT file) if DEFAULT_MAT_PATH is not set

    mat_path = args.mat_path or get_default_mat_path()
    fiber_index = DEFAULT_FIBER_INDEX if args.fiber_index is None else args.fiber_index
    figure_name = DEFAULT_FIGURE_NAME if args.figure_name is None else args.figure_name
    output_dir = args.output_dir or get_default_output_dir()

    print(f"Loading MAT file: {mat_path}")
    data = load_fiberphotometry_mat(mat_path, fiber_index=fiber_index)

    print(f"Creating Figure 1 (fiber index {fiber_index})...")
    create_figure_1(data, figure_name=figure_name, output_dir=output_dir)


if __name__ == "__main__":
    main()


