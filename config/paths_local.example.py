"""
Per-machine path overrides (Python) -- EXAMPLE TEMPLATE.
================================================================================

Copy this file to `paths_local.py` (same folder) and edit. `paths_local.py` is
gitignored, so each machine can point at local data without touching tracked
code. paths_config.py imports it LAST, so anything you set here wins.

Only redefine what differs on this machine; delete the rest.
"""

from pathlib import Path

# Example: data lives on a fast local drive on the acquisition PC.
# DATA_ROOT = Path(r"D:\Imaging_Data")

# Example: write figures to a local scratch folder instead of the network share.
# FIGURES_ROOT = Path(r"C:\fiber_outputs\Figures")
# SPECTRAL_OUTPUT_ROOT = FIGURES_ROOT / "Spectral_data_outputs"
