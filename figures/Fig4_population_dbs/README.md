# Figure 4 — Population DBS Response (+ Suppl. Fig. 2: Striatal DBS)

Reproduces the manuscript's Figure 4: population (fiber-photometry) GEVI and LFP responses to
40 Hz vs 135 Hz electrical deep-brain stimulation — trial-averaged traces, pre/transient/
sustained/post-stim epoch quantification, time-frequency spectrograms, and stimulation-band
power. `stimulation_analysis.py --mode striatum-dbs` produces the striatal-DBS composite
(Suppl. Fig. 2).

## Inputs

`*_FiberPhotometry_Analysis.mat` files from
[`../../preprocessing/fiber_photometry/`](../../preprocessing/fiber_photometry/) (stimulation
trials), and the corresponding coherence spectral outputs from
[`../../spectral_analysis/run_stim_spectral_pipeline.m`](../../spectral_analysis/run_stim_spectral_pipeline.m).

## Scripts

### `stimulation_analysis.py`

The main comparison figures (12 panel types: traces, stim-onset zoom, spectral heatmaps,
period violins, 40 vs 135 Hz comparison, trial-by-trial heatmaps). Select which
animal/condition pair to plot via `COMPARISON_MODE` at the top:
```python
COMPARISON_MODE = 'Animal04_EnergyBalanced'  # <-- CHANGE THIS
```
See the comment block above it for all available modes (`Animal01_40vs135`,
`Animal02_AmpBalanced`/`Animal02_EnergyBalanced` for 1s stimulation,
`Animal03_*`/`Animal04_*` for 10s stimulation). Edit `SESSION_R1`/`SESSION_R2` and `BASE_PATH`
within each mode's block to point at your own sessions.

Run:
```bash
python stimulation_analysis.py                          # default suite, uses COMPARISON_MODE above
python stimulation_analysis.py --mode striatum-dbs \
    --trial1-mat "<path to Trial1_..._FiberPhotometry_Analysis.mat>" \
    --num-trials 10                                       # Suppl. Fig. 2 composite
```

### `publication_composites.py`

6-panel (A-F) publication composite figures, one per stimulation condition, for two example
animals (`SESSIONS` dict at the top — edit `base_path`/`session_id` for your own data), plus an
all-conditions composite.

Run:
```bash
python publication_composites.py
```

## Configuration

Both scripts derive `DATA_ROOT`/`PROJECT_ROOT` from
[`../../config/paths_config.py`](../../config/paths_config.py) — set your data root there once.

## Output

`Figures/Stimulation_analysis/` and `Figures/CompFig/` under the project root.
