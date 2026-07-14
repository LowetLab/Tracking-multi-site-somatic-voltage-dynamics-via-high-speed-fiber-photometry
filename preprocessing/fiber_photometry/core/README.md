# `core/` — extracted fiber-preprocessing functions

Reusable, unit-tested functions pulled out of the two monolithic preprocessing
scripts (`run_fiber_preprocessing_{multitrial,singletrial}.m`)
as part of the Tier B decomposition. Every function here is **behaviour-preserving**:
it reproduces, bit-for-bit, the inline code it replaced (proven by the tests in
`tests/` and `utils/tests/`).

Both main scripts put this folder on the path **on top** of everything else
(`addpath(genpath(.../core))`), so these versions win over any external or
built-in functions of the same name.

## `core/` — pipeline computations

| Function | What it does |
|----------|--------------|
| `detect_trial_folders` | Find + numerically sort the per-trial subfolders (`..._1, _2, _10`). |
| `correct_photobleaching` | Linear detrend + double-exponential `lsqcurvefit` correction. Parameterised (`opts.BaselineClampMin`, `opts.ShortBaselineGuard`) so it reproduces **both** scripts' variants; returns a `fit_diag` struct the scripts fold back into the output. |
| `detect_stim_onset` | First rising edge of a pulse/trigger trace (`diff > threshold`), optionally the first onset after a trial start. |
| `compute_deltaF_F` | Static-baseline `ΔF/F = (F_corr − F0)/F0` per fiber, plus the z-scored version and `F0_values`. |
| `compute_band_power` | Mean band power `mean(|S|.^2)` per frequency band from a complex spectrogram. |
| `zscore_smooth_bands` | Per-band z-score then `fastsmooth` of a band-power matrix. |

## `core/utils/` — generic signal/display helpers

| Function | What it does |
|----------|--------------|
| `smooth2a` | Separable 2-D moving-average smoothing. |
| `fastsmooth` | 1-D smoothing (flat / triangular / Gaussian; two edge modes). |
| `viridis`, `turbo` | Perceptually-uniform colormaps (own copies; beat MATLAB's built-in `turbo`). |
| `generate_biphasic_pulses` | Synthetic biphasic stim pulse train for a time vector. |
| `replace_outliers_with_median` | Clamp `|z| > threshold` samples of a vector to its median (LFP cleaning; default 10). |
| `clean_display_frame` | Time-average an image stack and knock down bright (`z > threshold`) pixels for the ROI-selection background (display only; default 15). |

## `config/` (sibling folder, not under `core/`)

Run parameters live in `config/fiber_preprocessing_{multitrial,singletrial}_config.m`.
Each main script loads its config with `run(...)` (a script-include that
preserves every variable name). **Edit the config file, not the script top.**

## Tests

- Per-function unit tests: `core/tests/test_*.m`, `core/utils/tests/test_core_utils.m`.
- One-shot runner (all offline suites): from `preprocessing/fiber_photometry/`,
  ```matlab
  addpath('tests'); run_all_tests
  ```
- Offline refactor check: `tests/validate_refactor.m`.
- Post-run output check (needs a real `*_FiberPhotometry_Analysis.mat`):
  `tests/check_output_struct.m`. See `../VALIDATION.md`.

These are all "Level 1" checks (no data needed). The full pipeline still
requires a real recording + interactive ROI selection; a real-data run is the
final sign-off for any stage that touches the saved output.
