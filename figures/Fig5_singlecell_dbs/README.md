# Figure 5 — Single-Cell DBS Response

Reproduces the manuscript's Figure 5: cellular-resolution (single-neuron) GEVI voltage
responses to 40 Hz vs 130 Hz DBS — per-neuron and population-averaged traces, spike rasters,
period-quantification violins, and time-frequency spectrograms.

## Inputs

`*_CellularAnalysis.mat` files from
[`../../preprocessing/cellular_imaging/`](../../preprocessing/cellular_imaging/).

## Scripts

| Script | Produces |
|--------|----------|
| `cellular_dbs_traces.py` | Trial-averaged and single-trial example traces (stimulation, LFP, population/individual-neuron voltage) with spike markers, across a session list (`SESSIONS`) |
| `cellular_dbs_comparison.py` | 8-figure 40Hz-vs-135Hz comparison (traces, period violins, spectrograms, stim-band power) for one pair of sessions |

## Configuration

Both derive `DATA_ROOT`/`PROJECT_ROOT` from
[`../../config/paths_config.py`](../../config/paths_config.py). Edit at the top of each script:
```python
MOUSE_NAME = "Animal01"
```
`cellular_dbs_traces.py`'s `SESSIONS` list and `cellular_dbs_comparison.py`'s
`SESSION_R1`/`SESSION_R2` need your own session dates/IDs.

## Run

```bash
python cellular_dbs_traces.py
python cellular_dbs_comparison.py
```

## Output

`Figures/Cellular_DBS_traces/` and `Figures/Cellular_DBS_comparison/` under the project root.
