# Cellular Voltage Imaging Preprocessing Pipeline

## Overview

This pipeline processes cellular-resolution voltage imaging data (single-neuron ROIs) and
aligns it with LFP recordings from Open Ephys. It is designed for DBS (Deep Brain Stimulation)
analysis with support for comparing different stimulation frequencies. It produces the
single-cell recordings behind **Figure 5** of the manuscript.

## Directory Structure

```
cellular_imaging/
├── README.md                              # This file
├── cellular_processing_multitrial_DBS.m   # Main processing script
├── run_cellular_batch_processing.m        # Batch processing wrapper
├── extract_spike_ml_data.m                # Export a concise spike-centric struct
├── validate_spike_ml_data.m               # Visual sanity check of that export
├── spike_detect_SNR_sim3_SC.m             # Spike detection (SNR-based, default method)
└── config/
    └── cellular_session_database.m        # Session metadata database (edit with your cohort)
```

## Quick Start

### Single Session Processing

1. Open `cellular_processing_multitrial_DBS.m`
2. Modify the configuration section at the top:
   ```matlab
   MOUSE_NAME = 'Animal01';
   RECORDING_DATE = '01_06_25';
   RECORDING_ID = 'R2';
   DBS_FREQUENCY_HZ = 40;
   DBS_COMPARISON_TYPE = 'EnergyBalanced';
   ```
3. Run the script
4. Select ROIs when prompted (motion correction ROI, then neuron ROIs)

### Batch Processing

1. Open `run_cellular_batch_processing.m`
2. Edit `SESSIONS_TO_PROCESS` to list sessions:
   ```matlab
   SESSIONS_TO_PROCESS = {
       'Animal01', '01_06_25', 'R2';    % 40Hz Energy Balanced
       'Animal01', '01_06_25', 'R9';    % 130Hz
   };
   ```
3. Run the script

## Features

- **Multi-trial processing**: Automatically detects trial folders with `_1`, `_2`, `_3` suffixes
- **Multi-TIFF support**: Concatenates multiple TIFF files per trial
- **Automatic path construction**: Finds Open Ephys data automatically
- **ROI sharing**: Neuron ROIs selected on first trial are applied to all trials
- **Motion correction**: Optional rigid motion correction using NoRMCorre
- **Photobleaching correction**: Linear detrending + double exponential fitting
- **Spike detection**: SNR-based (default) or fixed baseline-window threshold
- **DBS metadata**: Tracks stimulation frequency and comparison type

## Data Organization

### Expected Input Structure

```
CellularVoltageImaging/
└── <MOUSE_NAME>/
    └── DBS/
        └── <RECORDING_DATE>/
            └── <RECORDING_ID>/
                ├── trial_1/
                │   ├── recording_001.tif
                │   └── recording_002.tif
                ├── trial_2/
                │   └── recording.tif
                └── ...
```

### Output Structure

```
Preprocessed_Data/Cellular/
└── <MOUSE_NAME>/
    └── <RECORDING_DATE>-<RECORDING_ID>/
        ├── <MOUSE_NAME>_<RECORDING_DATE>-<RECORDING_ID>_CellularAnalysis.mat
        └── figures/
            ├── <MOUSE_NAME>_<RECORDING_DATE>-<RECORDING_ID>_ROIs.png
            └── <MOUSE_NAME>_<RECORDING_DATE>-<RECORDING_ID>_ROIs.fig
```

## Session Database

The session database (`config/cellular_session_database.m`) contains metadata for all
recordings. This enables:

- Automatic parameter lookup
- Comparison pair identification (40Hz vs 130Hz)
- Processing status tracking

Replace the example entries with your own cohort.

### Adding New Sessions

Edit `config/cellular_session_database.m`:

```matlab
idx = idx + 1;
sessions(idx).mouse_name = 'NewMouse';
sessions(idx).recording_date = 'DD_MM_YY';
sessions(idx).recording_id = 'R1';
sessions(idx).experiment_type = 'DBS';
sessions(idx).dbs_frequency_hz = 40;
sessions(idx).dbs_comparison_type = 'EnergyBalanced';
sessions(idx).num_trials_expected = 5;
sessions(idx).num_neurons_expected = NaN;
sessions(idx).notes = 'Description';
sessions(idx).status = 'pending';
```

## DBS Comparison Types

For DBS experiments, two comparison types are supported:

1. **Energy Balanced**: 40Hz vs 130Hz - same total energy delivered
2. **Amplitude Balanced**: 40Hz vs 130Hz - same pulse amplitude

## Output Data Structure

The main output file (`*_CellularAnalysis.mat`) contains:

```matlab
CellularAnalysis
├── metadata
│   ├── mouse_name
│   ├── recording_date
│   ├── recording_id
│   ├── experiment_type
│   ├── dbs_frequency_hz
│   ├── dbs_comparison_type
│   ├── num_trials
│   └── num_neurons
├── shared_rois
│   ├── masks           % Cell array of neuron ROI masks
│   └── centroids       % Nx2 array of ROI centers
└── trials              % Cell array with per-trial data
    └── trial{n}
        ├── metadata
        ├── parameters
        ├── time
        │   ├── time_vector
        │   └── stimulus_onset_frame
        ├── signals
        │   ├── lfp_aligned
        │   ├── stimulus_aligned
        │   ├── fluorescence_raw
        │   ├── fluorescence_detrended
        │   ├── fluorescence_corrected
        │   └── subthreshold
        ├── rois
        ├── spikes
        │   ├── spike_raster
        │   ├── firing_rates_hz
        │   └── detection_results
        └── dbs
            ├── frequency_hz
            ├── comparison_type
            └── duration_sec
```

## Spike Detection

`spike_detect_SNR_sim3_SC.m` (bundled in this folder) is the default method
(`SPIKE_DETECTION_TYPE = 'baselineSD'`): SNR-thresholded upward-deflection detection with
optional spike-waveform extraction. See its own header comment for the full input/output
struct reference. The alternative fixed-threshold method
(`SPIKE_DETECTION_TYPE = 'thresholdSD'`) calls `spike_detect_baseline_threshold_SC.m`, which is
**not bundled** -- see below.

## Required Dependencies

The following external functions must be on the MATLAB path (not bundled with this repo --
point `config/lab_paths.m`'s `p.toolboxes` at your local copies, see `config/README.md`):

- `load_open_ephys_data.m` - Open Ephys data loading
- `NoRMCorreSetParms.m` / `normcorre.m` - Motion correction (NoRMCorre)
- `spike_detect_baseline_threshold_SC.m` - Spike detection (fixed baseline-window threshold; alternative to the bundled SNR-based method)
- `fastsmooth.m` - Fast smoothing function (also bundled in `../fiber_photometry/core/utils/`)
- `nanfastsmooth.m` - NaN-aware smoothing, used by `spike_detect_SNR_sim3_SC.m`

## Comparison with Fiber Pipeline

This pipeline is modeled after `../fiber_photometry/run_fiber_preprocessing_multitrial.m` but
adapted for cellular voltage imaging:

| Feature | Fiber Pipeline | Cellular Pipeline |
|---------|---------------|-------------------|
| ROI Selection | Pre-defined fiber positions | Interactive polygon drawing |
| Signal Type | Photometry (2 fibers) | Fluorescence (N neurons) |
| Motion Correction | Not needed | Optional (NoRMCorre) |
| Spike Detection | Not applicable | SNR/Baseline threshold |
| Photobleaching | Linear detrend | Linear + Double exponential |

## Troubleshooting

### "No trial folders found"
- Check that trial folders have `_1`, `_2`, `_3` suffixes
- Verify the data path is correct

### "No TIFF files found"
- Check that TIFF files are directly in trial folders
- Verify file extensions (.tif or .tiff)

### "Could not find OpenEphys data"
- The script tries multiple path patterns
- If automatic detection fails, it will prompt for manual selection

### Low computed frame rate
- If frame triggers are unreliable, the script uses manual frame rate (650 Hz)
- Adjust `MANUAL_FRAME_RATE` if your camera runs at different speed
