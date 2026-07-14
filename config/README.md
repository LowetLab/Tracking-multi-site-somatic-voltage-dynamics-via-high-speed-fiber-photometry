# `config/` — centralised paths (single source of truth)

This folder defines **all project paths in one place** for both MATLAB and Python,
so the pipeline is portable and you never have to edit a hardcoded absolute path
inside an analysis script.

## Design

Two anchors:

- **`project_root`** — *self-located* from the config file. Works wherever you clone
  or copy this repo. Code and figure outputs live under here.
- **`lab_root`** (MATLAB) / **`LAB_ROOT`** (Python) — an *explicit* default that you
  must edit for your own machine (`config/lab_paths.m` line ~14, `config/paths_config.py`
  line ~13). Your raw imaging/ephys data and any external MATLAB toolboxes live
  **outside** the repo, so they cannot be self-located — everything external hangs
  off this one anchor.

## Files

| File | Purpose |
|------|---------|
| `lab_paths.m`            | MATLAB: returns a struct `p` of all paths. |
| `setup_lab_paths.m`      | MATLAB: `setup_lab_paths()` adds all external toolboxes to the path (use instead of `addpath(genpath(...))` blocks). |
| `paths_config.py`        | Python: module exposing `PROJECT_ROOT`, `DATA_ROOT`, `SPECTRAL_OUTPUT_ROOT`, … |
| `paths_local.example.m`  | Template for per-machine MATLAB overrides. |
| `paths_local.example.py` | Template for per-machine Python overrides. |

## Usage

**MATLAB** — at the top of an entry-point script:

```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'config'));  % reach config/
setup_lab_paths();          % add external toolboxes to the path
p = lab_paths();            % get the paths struct
data_dir = fullfile(p.data_root, 'YourAnimalID', 'YourSessionID');
```

**Python** — from any script, regardless of depth:

```python
from paths_config import add_config_to_path
add_config_to_path()                       # ensures config/ is importable
from paths_config import DATA_ROOT, SPECTRAL_OUTPUT_ROOT
```

(or locate the repo without importing: `paths_config.find_project_root(__file__)`).

## Setting up your own paths (do this first)

1. **Quick start**: edit the default in `lab_paths.m` (`p.lab_root`) and
   `paths_config.py` (`LAB_ROOT`) directly to point at your own data location.
2. **Preferred** (keeps tracked files untouched, useful if you contribute changes
   back): copy `paths_local.example.m` → `paths_local.m` and/or
   `paths_local.example.py` → `paths_local.py`, and override only the fields that
   differ, e.g. `p.data_root = 'D:\Imaging_Data';`. These two files are
   **gitignored** and applied **last**, so your local settings always win and never
   get committed.
3. **External toolboxes**: `p.toolboxes` / the equivalent Python setup lists the
   external dependencies this pipeline expects on the path (FieldTrip, optionally
   NoRMCorre, an Open Ephys MATLAB loader, your own spike-detection scripts — see
   `environment/SETUP.md` and `preprocessing/cellular_imaging/README.md`). None of
   these are bundled with this repo; install them separately and point at your
   local copies.
