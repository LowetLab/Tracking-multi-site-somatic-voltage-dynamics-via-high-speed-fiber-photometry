# Figure 2 — Locomotion-Dependent Spectral & Coherence Dynamics

Reproduces the manuscript's Figure 2: time-frequency structure, spectral power, and
LFP-GEVI coherence during rest vs. running, at the single-trial, session-pooled, and
group (cross-animal) levels, plus the FOOOF aperiodic-corrected oscillatory power panel and
the theta-band summary violin plots.

## Inputs

Coherence/PSD `.mat` outputs from
[`../../spectral_analysis/`](../../spectral_analysis/) (`run_spectral_pipeline.m`), read via
`../common/plotting_config.py`'s path helpers.

## Scripts

| Script | Manuscript panel | Level |
|--------|------|-------|
| `fig2_coherence.py` | 2A (spectrogram/coherence heatmaps) + 2B/2D left column (example session) | Single-trial |
| `fig2_coherence_pooled.py` | Session- and animal-pooled coherence/PSD | Session-pooled, animal-pooled |
| `fig3_coherence_group.py` | 2B/2C/2D middle+right columns (population average + band-averaged stats) | Group-level (cross-animal) |
| `fig4_fooof.py` | 2C (aperiodic-corrected oscillatory GEVI power) | FOOOF 1/f decomposition of PSD |
| `fig5_theta.py` | 2B/2C/2D right column (theta-band half-violin summary, paired stats) | Group-level |

Despite the `fig2`/`fig3`/`fig4`/`fig5` filenames (an internal numbering left over from
development), **all five scripts feed manuscript Figure 2** — none of them correspond to
paper Figures 3, 4, or 5.

## Configuration

All five read shared settings from [`../common/plotting_config.py`](../common/plotting_config.py)
(behaviour mode, artifact-handling suffix, animal/session cohort, styling). Each script also
has a `USING_CENTRAL_CONFIG` fallback block used only if `plotting_config.py` can't be found —
you shouldn't need it if the repo layout is intact.

Recommended: run via the orchestrator, `run_all_plots.py`, which drives all single-trial,
session/animal-pooled, and group-level figures (plus FOOOF) in one pass, honoring the
`PLOT_*` flags and `ANIMALS_TO_PROCESS` in `../common/plotting_config.py`:
```bash
python run_all_plots.py
```
Or run any script standalone:
```bash
python fig2_coherence.py
python fig2_coherence_pooled.py
python fig3_coherence_group.py
python fig4_fooof.py
python fig5_theta.py
```

## Output

`Figures/Spectral_data_outputs{_artifact_suffix}/{behaviour_mode}/...` under the project root
(see `../../spectral_analysis/README.md` for the full output-tree layout).
