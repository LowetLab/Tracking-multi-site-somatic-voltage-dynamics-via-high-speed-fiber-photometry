# Figure 7 (panels E-H) — Dual-Color, Multi-Animal Fiber Voltage Imaging During DBS

Reproduces panels E-H of the manuscript's Figure 7: raw and 40Hz-filtered GEVI vs. mCherry
reference traces during 40Hz electrical stimulation, and the corresponding spectral-power-change
quantification (line plot + 40Hz bar graph), for 5 simultaneously-recorded fibers across 3
PV-Cre mice.

Panels A-D (optical setup schematic, cage photo, FOV image, histology) were assembled by hand
and are not reproduced here.

## Input

A single `.mat` file with one variable `all_traces`, shaped `(n_samples, n_channels, n_trials)`.
`n_channels = 2 * n_fibers`: a GEVI channel and a simultaneously-recorded mCherry reference
channel per fiber.

## Important: channel identity is auto-detected, not assumed

The script does **not** assume a fixed channel order for "which half is GEVI" or "which fiber is
CA1 vs. mPFC". Instead:

- **GEVI vs. reference** is detected from the data itself: GEVI channels show a rise in
  stimulation-frequency (40Hz) power during the stimulation window (the entire point of the
  experiment), while the static mCherry reference should not. See `detect_gevi_channels()`.
- **CA1 vs. mPFC grouping** (for the bar plots in G/H) is a data-driven proxy: fibers are ranked
  by entrainment strength, strongest first, since the manuscript states CA1 shows strong,
  consistent 40Hz entrainment in all three mice while mPFC shows weaker entrainment. See
  `rank_fibers_by_entrainment()`.

If you know the true channel order for your own data, set `GEVI_CHANNEL_IDX` / `REF_CHANNEL_IDX`
/ `FIBER_LABELS` explicitly at the top of the script — the auto-detection and ranking steps are
skipped whenever those are set.

## Configuration

Edit the block marked `USER CONFIGURATION` near the top of `fig7_dualcolor_multianimal.py`:

- `DATA_PATH` — where your `all_traces` `.mat` file lives.
- `FS`, `STIM_FREQ` — sampling rate and DBS frequency (defaults match the manuscript: 430Hz, 40Hz).
- `STIM_ONSET_SEC` / `STIM_DURATION_SEC` — stimulation timing within each trial. **Set these to
  match your recording** — the defaults (4s onset, 1s duration) are a reasonable placeholder
  based on the manuscript's general baseline-window convention, not a verified value for this
  specific dataset.
- `N_CA1_FIBERS` — how many of the ranked fibers to group as "CA1" vs. "mPFC" in panels G/H.

## Run

```bash
python fig7_dualcolor_multianimal.py
```

## Output

`Figures/Fig7_dualcolor_dbs/` under the project root (`Fig7E_traces`, `Fig7F_traces`,
`Fig7G_spectral`, `Fig7H_spectral`, each as `.png` and `.pdf`).
