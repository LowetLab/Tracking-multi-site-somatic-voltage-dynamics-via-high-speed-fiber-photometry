# Figure 6 — Bilateral CA1 Coupling

Reproduces the manuscript's Figure 6: simultaneous right/left CA1 GEVI+LFP recordings,
inter-hemispheric fiber-fiber theta cross-correlation (with circular-shift surrogate testing),
and theta phase-locking (rose plots + resultant vector length R, rest vs. running).

`multisite_fiber_analysis.py` is also the shared data-loading/styling module imported by
[`../Fig3_theta_gamma_pac/multisite_pac.py`](../Fig3_theta_gamma_pac/multisite_pac.py) — don't
rename or move it without updating that import.

## Inputs

`*_FiberPhotometry_Analysis.mat` files (bilateral fiber sessions) from
[`../../preprocessing/fiber_photometry/`](../../preprocessing/fiber_photometry/).

## What `python multisite_fiber_analysis.py` generates

Running the script end-to-end (`main()`) produces, in order:

| # | Figure | Manuscript mapping |
|---|--------|------|
| 1 | Representative traces (right/left HP LFP, GEVI, theta overlay, motion) | Fig 6C-D |
| 2 | Time-frequency spectrograms (speed, LFP TFR, fiber TFR, both sites) | — |
| 3 | Fiber-fiber theta cross-correlation (single-trial + all-trials, surrogate test) | **Fig 6E-G** |
| 3b | Same, unsplit (no REST/RUN separation) | — |
| 4 | Hippocampal theta phase-locking (rose plots, R violin) | **Fig 6H-J** |
| 5 | Bilateral spectral analysis (coherogram, PSD, fiber-fiber coherence) | — |
| 6/7/6b/7b | Fiber-LFP cross-correlation & phase-locking, 4 ipsi/contra combinations | Supports Fig 3 PAC analysis |
| LFP diag | LFP-LFP coherence diagnostic | — |

## Configuration

Edit the `RECORDINGS` dict and `EXAMPLE_ANIMAL`/`EXAMPLE_SESSION`/`EXAMPLE_TRIAL` (Figures 1-5)
and `FIG6_ANIMAL`/`FIG7_ANIMAL`/etc. (Figures 6/7/6b/7b/LFP-diag) near the top of the file —
these must be keys/sessions present in `RECORDINGS`. `BASE_PATH`/`OUTPUT_DIR` derive from
[`../../config/paths_config.py`](../../config/paths_config.py).

## Run

```bash
python multisite_fiber_analysis.py
```

Or import and call an individual figure function (e.g. `fig_xcorr(...)`, `fig_phase_locking(...)`)
from a Python session if you only need one panel — see the `def fig_*` functions for signatures.

## Output

`Figures/Multisite_Fiber_Analysis/` under the project root.
