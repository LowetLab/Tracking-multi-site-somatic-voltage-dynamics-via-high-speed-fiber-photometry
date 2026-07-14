# `*_FiberPhotometry_Analysis.mat` — output schema (contract)

The preprocessing scripts save one struct, `FiberPhotometryAnalysis`, per trial
(`-v7.3`). This is the **frozen interface** the downstream pipelines read
(`../../spectral_analysis/`, `../../phase_amplitude_coupling/`, and the
`../../figures/` plotting scripts). Treat field **names, shapes, and meanings as
a contract**: add fields freely, but do not rename/repurpose/remove existing
ones without updating every consumer.

This document is hand-derived from `SECTION 10` of
`run_fiber_preprocessing_multitrial.m` (the single-trial
script writes the same groups). `tests/check_output_struct.m` asserts the
load-bearing fields are present, non-empty and finite after a real run.

> Conventions below: **F** = number of frames (imaging samples), **N** = number
> of fibers/ROIs. Traces are `F × N` unless noted. "cond:" marks fields written
> only when a condition holds.

## Top-level groups
`metadata · parameters · time_periods · time · signals · ephys ·
photobleaching · rois · spectral · phase_locking · theta_analysis ·
lfp_comparison · plot_config · summary_stats`

### `metadata`
`mouse_name, experimenter, recording_date, recording_id, analysis_date,
analysis_mode, base_folder, save_directory, trial_name, trial_folder,
trial_index, total_trials`.

### `parameters`
`sampling_rate` (imaging Hz), `ephys_sampling_rate`, `num_frames` (F),
`num_fibers` (N), `recording_duration_sec`, `motion_correction`,
`correction_type`, `process_full_field`, `invert_trace`, `ephys_loaded`,
`polarity`, `is_baseline_trial`.

### `time_periods`
`is_baseline_trial`, `detection_method`, `detection_metadata`.
- **stim trial:** `pre_stim_period`, `stim_period` ([start end] s),
  `stim_onset_frame`, `stim_offset_frame`; cond: `hardcoded_period`.
- **baseline trial:** the four above are `[]`; instead `early_period`,
  `middle_period`, `late_period` (thirds of the recording).

### `time`
`time_vector_seconds` (F×1), `sampling_rate`.

### `signals`  (the core fiber traces — all `F × N`)
`raw_traces` → `filtered_traces` → `corrected_traces` (photobleach-corrected)
→ `deltaF_F_traces` (ΔF/F, static baseline) → `zscored_traces` →
`final_processed_traces` (= ΔF/F). Plus `deltaF_F_baseline_window`
([startFrame endFrame]), `deltaF_F_baseline_time`, `deltaF_F_method`,
`F0_values` (1×N). cond (photobleaching on): `detrended_traces`,
`exp_corrected_traces`.

### `ephys`  (cond: `ephys_loaded`)
HP (always, when ephys loaded): `lfp_HP_full_original` (cond), `lfp_raw_aligned_HP`,
`lfp_z_HP`, `lfp_sampling_rate`, `lfp_original_sampling_rate`.
Behaviour/stim: `running_velocity`, `running_velocity_smooth`, `stim_pulses`,
`stim_onset`, `camera_trigger_indices`, `theta_envelope_correlation`.
- cond mPFC (`LOAD_mPFC_LFP & loaded`): `mPFC_loaded`, `lfp_mPFC_full_original`
  (cond), `lfp_raw_aligned_mPFC`, `lfp_z_mPFC`, `mPFC_computation`.
- cond ipsiHP (`LOAD_ipsiHP_LFP & loaded`): `ipsiHP_loaded`,
  `lfp_ipsiHP_full_original` (cond), `lfp_raw_aligned_ipsiHP`, `lfp_z_ipsiHP`,
  `ipsiHP_computation`.

### `photobleaching`
`correction_applied`. cond (applied): `methods_used`, `fit_success`, and (on a
successful double-exp fit) `fit_params`, `pre_stim_trace`, `time_pre_stim`,
`time_full`, `fitted_curve`. Per-fiber `fiber(i).bleaching_rate`,
`fiber(i).bleach_trend`.

### `rois`
`type` (`'full_field'` | `'manual'`), `roi_data` (the ROI masks/coords).

### `spectral`
`band_names`, `band_ranges` (nBands×2), `band_colors`. Per-fiber struct array
`fiber(i)` with the fiber spectrogram, freqs, times and `band_power`
(nBands×T). cond (ephys): `lfp(...)`; cond mPFC/ipsiHP: `lfp_mPFC`,
`lfp_ipsiHP`. (Exact leaf names: see SECTION 10.)

### `phase_locking`  (cond: ephys)
`freq_vector`, `time_vector`, `phase_fiber`, `phase_lfp`, `phase_difference`,
`plv_all`; behaviour split: `running_threshold`, `running_periods`,
`rest_periods`, `plv_running`, `plv_rest`.

### `theta_analysis`  (cond: ephys)
`lfp_filtered_theta`, `lfp_envelope`, `lfp_envelope_smooth`, `filter_settings`,
and per-fiber `fiber(i)` (theta-band fiber envelope / correlation).

### `lfp_comparison`  (cond: mPFC loaded)
`mPFC.(...)` — HP-vs-mPFC comparison metrics.

### `plot_config`
`colors`, `band_colors`, `viridis_period_colors` — colours used by the figures
(saved so figures can be reproduced).

### `summary_stats`
Per-fiber RMS by period and ratios (pre/stim/post, or early/middle/late for
baseline trials).

---
**If you change this struct**, update: `tests/check_output_struct.m`, the MATLAB
spectral loader (`../../spectral_analysis/core/spectral_analysis.m`), and the
Python readers in `../../figures/`.
