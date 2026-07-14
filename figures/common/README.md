# `figures/common/` — shared plotting infrastructure

Imported by every script under `figures/Fig*/` and `figures/Supplementary/`. Each figure
script locates this folder at runtime via `sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))`
(and, separately, `../../config/paths_config.py` for data/output roots) — you don't need to
add anything to your `PYTHONPATH` manually.

## Files

| File | Purpose |
|------|---------|
| `common.py` | De-duplicated helpers previously copy-pasted across figure scripts: `to_long_path` (Windows long-path handling), colormap builders, small numeric utilities. Unit-tested in `tests/test_common.py`. |
| `plotting_config.py` | Central config for the **spectral figures** (Fig1, Fig2): behaviour mode, artifact-handling suffix, the animal/session cohort database (mirrors `../../spectral_analysis/config/animal_session_database.m`), figure styling constants (fonts, colors, line widths), and path-derivation helpers (`get_single_trial_input_dir`, etc.). |

## Usage

Any figure script does:
```python
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from common import to_long_path
from plotting_config import BEHAVIOR_MODE, get_single_trial_input_dir, ...
```

If you're editing the cohort, edit `plotting_config.py`'s `ANIMALS` list (and keep it in sync
with the MATLAB `animal_session_database.m` -- they're maintained separately but describe the
same underlying sessions).

## Tests

```bash
python tests/test_common.py
```
