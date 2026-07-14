# Figure 3 — Theta-Gamma Phase-Amplitude Coupling

Reproduces the manuscript's Figure 3: LFP theta-phase-aligned wavelet spectrograms, phase-
resolved gamma-amplitude curves, and Tort's Modulation Index (MI) for LFP and GEVI signals
during rest vs. running.

## Compute (MATLAB) → Plot (Python)

The heavy computation (phase-aligned wavelet spectrograms, MI + IAAFT surrogate testing) is
MATLAB, in [`../../phase_amplitude_coupling/`](../../phase_amplitude_coupling/) — **read that
folder's README first**; it documents the full compute pipeline (Steps 1-6) and output layout.

`multisite_pac.py` here is the plotting layer: it loads the MATLAB PAC outputs and produces
per-combination figures (LFP row + Fiber row, each with spectrogram + phase curves, RUN vs
REST) plus a summary statistics figure comparing MI across combinations.

`multisite_pac.py` imports data-loading and styling helpers from `multisite_fiber_analysis.py`
in [`../Fig6_bilateral_ca1/`](../Fig6_bilateral_ca1/) (the bilateral-recording dataset is
shared between Figures 3 and 6) — this is resolved automatically via `sys.path`, no action
needed on your part.

## Configuration

Edit the top of `multisite_pac.py`:
```python
PAC_ANIMAL = "Animal02"           # must be a key in RECORDINGS (Fig6_bilateral_ca1/multisite_fiber_analysis.py)
PAC_AGGREGATE_ALL_SESSIONS = True
WAVELET_FREQ_RANGE = (5, 90)
```

## Run

```bash
python multisite_pac.py
```

## Output

`Figures/Multisite_Fiber_Analysis/PAC/` under the project root.
