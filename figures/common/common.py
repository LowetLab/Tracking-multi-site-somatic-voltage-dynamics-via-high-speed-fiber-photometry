"""
common.py -- shared helpers for the publication-figure scripts in figures/.

Single home for logic that was previously copy-pasted (and silently drifting)
across the individual ``fig*`` / analysis scripts. Import from here instead of
re-defining locally, e.g.::

    from common import to_long_path

Grow this module by section (paths, IO, signal, plotting) as duplication is
consolidated. Behaviour-preserving: each helper reproduces the canonical copy it
replaced. Pure helpers are unit-tested in ``tests/test_common.py``.
"""
from __future__ import annotations

import os
import re

import numpy as np
from scipy.io import loadmat


# --------------------------------------------------------------------------- #
#  Paths
# --------------------------------------------------------------------------- #
def to_long_path(path_str):
    r"""
    Convert a path to Windows extended-length format to bypass MAX_PATH (260 char) limit.
    Works for both UNC paths (\\server\share) and local paths (C:\folder).
    """
    path_str = str(path_str)

    # Already extended-length format
    if path_str.startswith('\\\\?\\'):
        return path_str

    # UNC path (network share)
    if path_str.startswith('\\\\'):
        # \\server\share -> \\?\UNC\server\share
        return '\\\\?\\UNC\\' + path_str[2:]

    # Local path with drive letter
    if len(path_str) >= 2 and path_str[1] == ':':
        return '\\\\?\\' + path_str

    # Return as-is for relative paths
    return path_str


def normalize_unc_path(path_str, for_access=False):
    r"""
    Canonicalize a path, optionally converting it to extended-length format.

    First strips any existing extended-length prefix (``\\?\UNC\`` or ``\\?\``)
    back to a plain UNC/local path. If ``for_access`` is True, the plain path
    is then converted to extended-length format via ``to_long_path`` so file
    I/O bypasses the Windows 260-character MAX_PATH limit.

    Use ``for_access=False`` (default) to get a normalized plain path for
    display/logging/dict-key purposes; use ``for_access=True`` right before
    a file-system call (``os.path.exists``, ``open``, ``loadmat``, ...).
    """
    path_str = str(path_str)

    if path_str.startswith('\\\\?\\UNC\\'):
        path_str = '\\\\' + path_str[8:]
    elif path_str.startswith('\\\\?\\'):
        path_str = path_str[4:]

    if for_access:
        path_str = to_long_path(path_str)

    return path_str


# --------------------------------------------------------------------------- #
#  IO
# --------------------------------------------------------------------------- #
def load_matlab_struct(filepath):
    """
    Load a MATLAB .mat file and return its variables as a dict.

    Uses the Windows extended-length path prefix to bypass the 260-character
    MAX_PATH limit, strips MATLAB's internal ``__header__``/``__version__``/
    ``__globals__`` keys, and returns None (after printing a message) if the
    file is missing or fails to load, rather than raising.
    """
    filepath_str = str(filepath)
    filepath_long = to_long_path(filepath_str)

    if not os.path.exists(filepath_long):
        return None
    try:
        data = loadmat(filepath_long, squeeze_me=True, struct_as_record=False)
        return {k: v for k, v in data.items() if not k.startswith('__')}
    except Exception as e:
        print(f"  ERROR loading {filepath}: {e}")
        return None


# --------------------------------------------------------------------------- #
#  Math / signal
# --------------------------------------------------------------------------- #
def _next_pow2(n):
    return int(2 ** np.ceil(np.log2(n)))


# --------------------------------------------------------------------------- #
#  Naming
# --------------------------------------------------------------------------- #
def _infer_trial_from_name(file_name):
    m = re.search(r"_Trial(\d+)_FiberPhotometry_Analysis\.mat$", file_name, re.IGNORECASE)
    return int(m.group(1)) if m else None


# --------------------------------------------------------------------------- #
#  Colormaps
# --------------------------------------------------------------------------- #
def create_monochromatic_orange_cmap():
    """
    Create a high-contrast monochromatic orange colormap for the speed heatmap.
    Goes from very dark orange (low speed / rest) to bright warm cream (high speed / running).

    Orange is complementary to the teal tones in viridis, creating visual contrast
    while maintaining aesthetic harmony.
    """
    from matplotlib.colors import LinearSegmentedColormap
    colors = [
        (0.25, 0.06, 0.00),    # Very dark orange/rust (low speed / rest)
        (0.45, 0.12, 0.02),    # Dark orange
        (0.65, 0.22, 0.04),    # Deep orange
        (0.85, 0.38, 0.08),    # Rich orange
        (0.98, 0.55, 0.18),    # Vivid orange
        (1.00, 0.75, 0.45),    # Warm golden orange
        (1.00, 0.95, 0.82),    # Warm cream (high speed / running)
    ]
    return LinearSegmentedColormap.from_list('mono_orange', colors, N=256)


def create_parula_like_cmap():
    from matplotlib.colors import LinearSegmentedColormap
    colors = [
        (0.208, 0.166, 0.529),
        (0.212, 0.325, 0.616),
        (0.192, 0.475, 0.635),
        (0.165, 0.600, 0.580),
        (0.220, 0.690, 0.490),
        (0.420, 0.760, 0.380),
        (0.660, 0.810, 0.320),
        (0.880, 0.850, 0.290),
        (0.976, 0.910, 0.145),
    ]
    return LinearSegmentedColormap.from_list('parula_like', colors, N=256)
