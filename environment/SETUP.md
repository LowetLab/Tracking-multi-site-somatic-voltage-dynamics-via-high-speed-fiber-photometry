# Setup

How to get this pipeline running on a machine. Paths are centralised in
[`config/`](../config/README.md) — you should never have to edit a hardcoded absolute path
inside an analysis script; set your data root once (see step 4) and every script derives its
paths from it.

## 1. Get the code

```bash
git clone https://github.com/LowetLab/Tracking-multi-site-somatic-voltage-dynamics-via-high-speed-fiber-photometry.git
```

## 2. Python environment

The Python code (figure plotting) needs Python 3.10+.

```bash
# conda (recommended):
conda env create -f environment.yml
conda activate fiber-pipeline

# or pip:
pip install -r requirements.txt
```

Packages: numpy, scipy, matplotlib, h5py, fooof (optional: mat73, for large MATLAB v7.3 files
in `phase_amplitude_coupling/`).

## 3. MATLAB

- **MATLAB R2016b+** (developed/tested on R2024b).
- Toolboxes: **Signal Processing Toolbox**; **Parallel Computing Toolbox** (optional — only
  for `parfor` in the PAC MI surrogates; replace with `for` if absent).
- **FieldTrip** — required for the spectral pipeline's multitaper coherence and
  cluster-based permutation stats. Download from
  <https://www.fieldtriptoolbox.org/download/> and add to the MATLAB path.

### External dependencies (not bundled)

The entry scripts call **`setup_lab_paths()`** (from `config/`), which adds toolbox folders
listed in `config/lab_paths.m`'s `p.toolboxes` to the MATLAB path; missing ones are skipped
silently (or with a warning if you call `setup_lab_paths(true)`). Install these separately and
point `p.toolboxes` at your local copies:

| Provides | Where to get it |
|----------|------------------|
| `load_open_ephys_data.m` | [Open Ephys](https://open-ephys.org/) MATLAB tools |
| NoRMCorre (`normcorre`, `NoRMCorreSetParms`) | [NoRMCorre](https://github.com/flatironinstitute/NoRMCorre) (optional motion correction) |
| `spike_detect_baseline_threshold_SC.m` | your own spike-detection script, if you use the alternative fixed-threshold method (see `preprocessing/cellular_imaging/README.md`) -- the default SNR-based method is bundled |
| `fastsmooth.m` | already bundled in `preprocessing/fiber_photometry/core/utils/` |
| `nanfastsmooth.m` | NaN-aware smoothing, used by the bundled spike-detection script |

## 4. Set your data root (do this first)

Defaults are placeholders (`C:\PATH\TO\YOUR\DATA_SHARE`) that won't resolve on your machine.
Either:
- **Quick**: edit the default directly in `config/lab_paths.m` (`p.lab_root`) and
  `config/paths_config.py` (`LAB_ROOT`), or
- **Preferred**: copy `config/paths_local.example.m` → `config/paths_local.m` and
  `config/paths_local.example.py` → `config/paths_local.py`, then edit only the fields that
  differ on your machine.

Both `paths_local.*` files are gitignored and applied last, so local settings always win and
never get committed. See [`config/README.md`](../config/README.md).

## 5. Where things go

- Project (code + figures): self-located from the repo root.
- Raw/processed **data**: under whatever you set as your data root in step 4.
- **Figure / spectral outputs**: under the repo's `Figures/` (gitignored).
