# Figure 1 — Platform Validation

Reproduces the manuscript's Figure 1: representative GEVI/LFP/motion traces from the CA1
fiber-photometry platform, and the time-resolved GEVI-LFP theta cross-correlation that
validates the recording (panels C-D).

## Inputs

`*_FiberPhotometry_Analysis.mat` files produced by
[`../../preprocessing/fiber_photometry/`](../../preprocessing/fiber_photometry/).

## Scripts

| Script | Produces | Run |
|--------|----------|-----|
| `fig1_gevi_lfp.py` | Panel C: raw ΔF/F GEVI, LFP, motion, theta/beta-filtered overlays for one trial | `python fig1_gevi_lfp.py --mat-path <path> --fiber-index 0` (or edit `DEFAULT_MAT_PATH` at the top and run with no args) |
| `fig1_cross_correlation.py` | Panel D: time-resolved GEVI-LFP cross-correlation heatmaps (RUN vs REST) and grand-mean correlogram, pooled across a session's trials, with circular-shift surrogate significance testing | `python fig1_cross_correlation.py` (edit `ANIMAL_ID`/`SESSION_ID`/`NUM_TRIALS` at the top) |
| `print_xcorr_stats.py` | Prints the same surrogate-test statistics as text (no figure) for the **bilateral** cross-correlation used in Figure 6 | `python print_xcorr_stats.py [animal_name]` |

`print_xcorr_stats.py` imports from `multisite_fiber_analysis.py` in
[`../Fig6_bilateral_ca1/`](../Fig6_bilateral_ca1/) — it reports statistics for the same
right-HP/left-HP dataset used there, not for `fig1_cross_correlation.py`'s single-site data.

## Configuration

Edit the `USER CONFIGURATION` block at the top of each script: `ANIMAL_ID`, `SESSION_ID`,
`NUM_TRIALS`, `FOLDER_SUFFIX`. Paths are derived automatically from
[`../../config/paths_config.py`](../../config/paths_config.py) — set your data root there (or
in a local `paths_local.py` override) once, rather than per-script.

## Output

`Figures/CrossCorrelation_figures/` and `Figures/Traces_figures/` under the project root
(PNG + PDF, or PNG/PDF/SVG depending on script).
