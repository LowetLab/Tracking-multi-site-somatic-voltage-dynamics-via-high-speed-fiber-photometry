r"""
Canonical paths for this pipeline (single source of truth, Python side).
================================================================================

All paths are derived by self-location from this file, so moving / copying the
project tree (preserving its internal structure) just works -- nothing is
hardcoded to a specific drive.

Two anchors:
    * PROJECT_ROOT -- SELF-LOCATED from this file (works from a clone anywhere
                      on disk). Code and figure outputs live here.
    * LAB_ROOT     -- an EXPLICIT default you must edit (see below). Your data
                      lives OUTSIDE the repo, so it cannot be derived from
                      PROJECT_ROOT. Change it per machine via paths_local.py.

Usage:
    # from any script, regardless of depth:
    from paths_config import PROJECT_ROOT, DATA_ROOT, SPECTRAL_OUTPUT_ROOT
    # (ensure <project_root>/config is on sys.path -- see add_config_to_path below)

Per-machine overrides:
    Copy paths_local.example.py to paths_local.py (gitignored) and redefine any
    of the names below. If present it is imported last and wins.
"""

from __future__ import annotations

from pathlib import Path

# --- Self-located project root -----------------------------------------------
CONFIG_DIR = Path(__file__).resolve().parent          # <project_root>/config
PROJECT_ROOT = CONFIG_DIR.parent                      # <project_root>

# --- Explicit data-share anchor (override in paths_local.py per machine) ------
# Data lives outside the repo -- EDIT THIS to point at your own data location,
# or (preferred) leave it and set DATA_ROOT in a local paths_local.py override.
LAB_ROOT = Path(r"C:\PATH\TO\YOUR\DATA_SHARE")

# --- Data & outputs ----------------------------------------------------------
DATA_ROOT = LAB_ROOT / "Data"          # raw imaging + ephys data root
FIGURES_ROOT = PROJECT_ROOT / "Figures"
SPECTRAL_OUTPUT_ROOT = FIGURES_ROOT / "Spectral_data_outputs"


def add_config_to_path() -> None:
    """Put <project_root>/config on sys.path so `import paths_config` works
    from any script regardless of its directory depth."""
    import sys

    cfg = str(CONFIG_DIR)
    if cfg not in sys.path:
        sys.path.insert(0, cfg)


def find_project_root(start) -> Path:
    """Walk up from `start` until a repo marker (.git / README.md) is found.
    Fallback locator for scripts that prefer not to import this module."""
    p = Path(start).resolve()
    for d in (p, *p.parents):
        if (d / ".git").exists() or (d / "README.md").exists():
            return d
    return PROJECT_ROOT


# --- Per-machine overrides (optional, gitignored) ----------------------------
try:  # pragma: no cover - environment specific
    from paths_local import *  # noqa: F401,F403
except ImportError:
    pass
