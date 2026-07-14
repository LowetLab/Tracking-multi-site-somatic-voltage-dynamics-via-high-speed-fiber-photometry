"""
================================================================================
CENTRAL CONFIGURATION FOR PYTHON PLOTTING PIPELINE
================================================================================

This module provides centralized configuration for all Python plotting scripts,
matching the MATLAB spectral analysis pipeline structure.

Usage:
    from plotting_config import CONFIG, ANIMALS
================================================================================
"""

from pathlib import Path
import numpy as np

# ==============================================================================
#  USER-CONFIGURABLE PARAMETERS
# ==============================================================================

# Behavior classification mode: 'standard' or 'clear'
# Must match the MATLAB pipeline's BEHAVIOR_MODE setting
BEHAVIOR_MODE = 'clear'

# Output folder suffix (to match MATLAB's OUTPUT_FOLDER_SUFFIX)
# ARTIFACT HANDLING MODES:
#   ''                   - No artifact handling (Spectral_data_outputs/)
#   '_artifact_excluded' - Trial exclusion: Trials with >30% contamination skipped
#   '_artifact_cleaned'  - Data cleaning: Artifact segments removed from all trials
#
# IMPORTANT DISTINCTION:
#   'excluded' = entire problematic TRIALS are skipped (coarse)
#   'cleaned'  = artifact TIME SEGMENTS removed from within trials (fine-grained)
OUTPUT_FOLDER_SUFFIX = '_artifact_cleaned'  # Set to match MATLAB pipeline

# Animals to process (set to None or [] to process all)
# Example: ['Animal01', 'Animal02'] or None for all
ANIMALS_TO_PROCESS = ['Animal01']

# Sessions to process (set to None or [] to process all sessions for selected animals)
# Example: ['01_09_25-R1'] to only process that session
SESSIONS_TO_PROCESS = ['01_09_25-R1']

# Analysis levels to plot
PLOT_SINGLE_TRIAL = True
PLOT_SESSION_POOLED = False
PLOT_ANIMAL_POOLED = False
PLOT_ANIMAL_CONCATENATED = False # Concatenated raw data across all sessions
PLOT_GROUP_LEVEL = False

# Optional additional figures
PLOT_PSD_LOGLOG = True  # Log-log scale PSD figure (for 1/f analysis)

# Methods to process (for single-trial, session-pooled, animal-pooled levels)
METHODS = ['mscohere', 'fieldtrip']

# ==============================================================================
#  GROUP-LEVEL FIGURE OPTIONS
# ==============================================================================
# Which group-level methods to generate figures for:
#   - 'fieldtrip_cluster': FieldTrip with cluster-based permutation statistics (RECOMMENDED)
#   - 'mscohere': mscohere with standard statistics 
#   - 'fieldtrip': FieldTrip with standard statistics (not cluster-corrected)
#
# NOTE: Cluster-based statistics are the gold standard for neural data.
# The non-cluster methods are optional and generate additional figures.

GROUP_LEVEL_METHODS = ['fieldtrip_cluster']  # Primary method (cluster-corrected)

# Set to True to also generate non-cluster FieldTrip and mscohere figures
# (These use standard statistics, not cluster-corrected - generates many more figures)
INCLUDE_NON_CLUSTER_METHODS = False

# Generate band-averaged boxplots?
# - Always generated for cluster method
# - Also generated for non-cluster methods IF INCLUDE_NON_CLUSTER_METHODS = True
GENERATE_BAND_BOXPLOTS = True

# Show session-pooled data alongside animal-pooled?
# - False (DEFAULT): Only show animal-pooled (the statistically correct unit of replication)
# - True: Show both session-pooled (left panel) and animal-pooled (right panel)
# NOTE: Session-pooled is PSEUDOREPLICATION - sessions from same animal are not independent.
#       Only use for supplementary/exploratory analysis.
SHOW_SESSION_POOLED = False

# ==============================================================================
#  GROUP-LEVEL POOLING CONFIGURATION
# ==============================================================================
# Which animal-level pooling to use for group statistics:
#   'animal_pooled'      - Spectra computed per session, then averaged (default)
#   'animal_concatenated' - All raw data concatenated, spectra computed once
#
# animal_concatenated provides a single spectrum from more data points but requires
# proper edge artifact handling for FieldTrip analysis.
GROUP_POOLING_LEVEL = 'animal_concatenated'  # 'animal_pooled' or 'animal_concatenated'

# Figure output formats
FIGURE_FORMATS = ['png', 'pdf', 'svg']
FIGURE_DPI = 300

# Frequency range for display (Hz)
FREQ_MIN = 2
FREQ_MAX = 70

# Coherence colorbar limits
COH_VMIN = 0.0
COH_VMAX = 1.0


# ==============================================================================
#  PATH CONFIGURATION
# ==============================================================================

# Base output directory (same as MATLAB) -- from the centralised config
# (config/paths_config.py). Override per machine via config/paths_local.py.
import sys
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import SPECTRAL_OUTPUT_ROOT
BASE_OUTPUT_DIR = SPECTRAL_OUTPUT_ROOT

# Derived paths based on behavior mode and output folder suffix
def get_base_dir():
    """Get base directory for current behavior mode.
    
    Structure: {BASE_OUTPUT_DIR}{SUFFIX}/{BEHAVIOR_MODE}/
    Example: Spectral_data_outputs_artifact_excluded/clear/
    """
    # Apply suffix to ROOT directory (not behavior mode)
    effective_root = Path(str(BASE_OUTPUT_DIR) + OUTPUT_FOLDER_SUFFIX)
    return effective_root / BEHAVIOR_MODE


# ==============================================================================
#  ANIMAL DATABASE  -- EDIT THIS to your own cohort (mirrors
#  ../../spectral_analysis/config/animal_session_database.m)
# ==============================================================================

ANIMALS = [
    {
        'mouse_id': 'Animal01',
        'project': 'FiberVoltageImaging',
        'sessions': [
            {'session_id': '01_09_25-R1', 'num_trials': 6},
            {'session_id': '02_09_25-R1', 'num_trials': 6},
            {'session_id': '03_09_25-R1', 'num_trials': 2},
            {'session_id': '03_09_25-R2', 'num_trials': 2},
        ],
        # Sessions to pool together for session-pooled analysis
        'session_pooled_groups': [
            ['03_09_25-R1', '03_09_25-R2'],  # Combined as '03_09_25-R1-combined'
        ],
    },
    {
        'mouse_id': 'Animal02',
        'project': 'FiberVoltageImaging',
        'sessions': [
            {'session_id': '01_01_26-R1', 'num_trials': 5},
            {'session_id': '01_01_26-R3', 'num_trials': 6},
            # DBS Stimulation sessions
            {'session_id': '01_02_26-R6', 'num_trials': 10, 'type': 'stim', 'stim_freq': '135Hz'},
            {'session_id': '01_02_26-R9', 'num_trials': 10, 'type': 'stim', 'stim_freq': '40Hz_AmpBalanced'},
            {'session_id': '01_02_26-R10', 'num_trials': 10, 'type': 'stim', 'stim_freq': '40Hz_EnergyBalanced'},
        ],
        'session_pooled_groups': [],
    },
]


def get_animals_to_process():
    """
    Get list of animals to process based on ANIMALS_TO_PROCESS and 
    SESSIONS_TO_PROCESS settings.
    
    Returns
    -------
    list
        List of animal dictionaries to process (with sessions filtered if configured)
    """
    if ANIMALS_TO_PROCESS is None or len(ANIMALS_TO_PROCESS) == 0:
        animals = ANIMALS
    else:
        animals = [a for a in ANIMALS if a['mouse_id'] in ANIMALS_TO_PROCESS]
    
    # Filter sessions if SESSIONS_TO_PROCESS is set
    if SESSIONS_TO_PROCESS is not None and len(SESSIONS_TO_PROCESS) > 0:
        filtered = []
        for a in animals:
            a_copy = dict(a)
            a_copy['sessions'] = [s for s in a['sessions'] 
                                  if s['session_id'] in SESSIONS_TO_PROCESS]
            if a_copy['sessions']:
                filtered.append(a_copy)
        return filtered
    
    return animals


def get_session_ids(animal):
    """
    Get all session IDs for an animal.
    
    Parameters
    ----------
    animal : dict
        Animal dictionary from ANIMALS
    
    Returns
    -------
    list
        List of session ID strings
    """
    return [s['session_id'] for s in animal['sessions']]


def get_trial_count(animal, session_id):
    """
    Get number of trials for a specific session.
    
    Parameters
    ----------
    animal : dict
        Animal dictionary from ANIMALS
    session_id : str
        Session identifier
    
    Returns
    -------
    int
        Number of trials, or 0 if session not found
    """
    for sess in animal['sessions']:
        if sess['session_id'] == session_id:
            return sess['num_trials']
    return 0


# ==============================================================================
#  PUBLICATION STYLING CONSTANTS
# ==============================================================================

# Font sizes
FONT_SIZE_SUPTITLE = 20
FONT_SIZE_TITLE = 18
FONT_SIZE_LABEL = 16
FONT_SIZE_TICK = 14
FONT_SIZE_LEGEND = 12
FONT_SIZE_STATS = 13
FONT_SIZE_BAND = 13

# Axis styling
AXIS_LINEWIDTH = 2.0
TICK_WIDTH = 1.8
TICK_LENGTH = 7

# Line widths
LINE_WIDTH_TRACE = 2.5
LINE_WIDTH_DASHED = 2.5
LINE_WIDTH_BAND = 1.0


# ==============================================================================
#  COLOR DEFINITIONS
# ==============================================================================

# REST vs RUN colors (Teal shades)
COLOR_REST = np.array([0.05, 0.35, 0.45])      # Darker teal for REST
COLOR_RUN = np.array([0.25, 0.65, 0.65])       # Lighter teal for RUN
COLOR_OVERALL = np.array([0.08, 0.45, 0.52])   # Mid teal for overall

# Signal type colors
COLOR_GEVI = np.array([0.127568, 0.566949, 0.550556])  # Teal (from viridis)
COLOR_LFP = np.array([0.35, 0.25, 0.45])               # Purple-grey

# LFP REST/RUN colors
COLOR_LFP_REST = np.array([0.25, 0.18, 0.35])
COLOR_LFP_RUN = np.array([0.55, 0.45, 0.65])

# GEVI REST/RUN colors
COLOR_GEVI_REST = np.array([0.05, 0.35, 0.45])
COLOR_GEVI_RUN = np.array([0.25, 0.65, 0.65])

# Coherence colors (Indigo/blue-violet blend - represents LFP-GEVI coupling)
# This is a blend between LFP purple and GEVI teal
COLOR_COH_REST = np.array([0.20, 0.25, 0.50])  # Dark indigo
COLOR_COH_RUN = np.array([0.45, 0.50, 0.70])   # Periwinkle/light indigo
COLOR_COH_OVERALL = np.array([0.30, 0.35, 0.58])  # Mid indigo for overall

# Frequency band line colors (with alpha)
BAND_LINE_COLORS = {
    'theta': (0.4, 0.2, 0.6, 0.6),    # Purple
    'alpha': (0.2, 0.5, 0.5, 0.6),    # Teal
    'beta': (0.3, 0.6, 0.3, 0.6),     # Green
    'gamma': (0.7, 0.5, 0.2, 0.6),    # Orange/brown
}

# Frequency band colors for box plots
BAND_COLORS = {
    'theta': np.array([0.4, 0.2, 0.6]),
    'alpha': np.array([0.2, 0.5, 0.5]),
    'beta': np.array([0.3, 0.6, 0.3]),
    'gamma': np.array([0.7, 0.5, 0.2]),
}

# Frequency band boundaries (Hz)
FREQ_BANDS = {
    'theta': (4, 8, 'θ'),
    'alpha': (8, 12, 'α'),
    'beta': (12, 30, 'β'),
    'gamma': (30, 70, 'γ'),
}

# SEM shading alpha
SEM_ALPHA = 0.25

# Significance marker color
COLOR_SIG = 'red'
SIG_MARKER_SIZE = 10


# ==============================================================================
#  PATH HELPER FUNCTIONS
# ==============================================================================

def get_single_trial_input_dir(mouse_id, session_id):
    """Get input directory for single-trial data."""
    return get_base_dir() / 'single_trial' / mouse_id / session_id / 'data'


def get_single_trial_output_dir(mouse_id, session_id):
    """Get output directory for single-trial figures."""
    return get_base_dir() / 'single_trial' / mouse_id / session_id / 'figures'


def get_session_pooled_input_dir(mouse_id, session_id):
    """Get input directory for session-pooled data."""
    return get_base_dir() / 'session_pooled' / mouse_id / session_id / 'data'


def get_session_pooled_output_dir(mouse_id, session_id):
    """Get output directory for session-pooled figures."""
    return get_base_dir() / 'session_pooled' / mouse_id / session_id / 'figures'


def get_animal_pooled_input_dir(mouse_id):
    """Get input directory for animal-pooled data."""
    return get_base_dir() / 'animal_pooled' / mouse_id / 'data'


def get_animal_pooled_output_dir(mouse_id):
    """Get output directory for animal-pooled figures."""
    return get_base_dir() / 'animal_pooled' / mouse_id / 'figures'


def get_animal_concatenated_input_dir(mouse_id):
    """Get input directory for animal-concatenated data."""
    return get_base_dir() / 'animal_concatenated' / mouse_id / 'data'


def get_animal_concatenated_output_dir(mouse_id):
    """Get output directory for animal-concatenated figures."""
    return get_base_dir() / 'animal_concatenated' / mouse_id / 'figures'


def get_group_level_input_dir():
    """Get input directory for group-level data.
    
    MATLAB saves to: {BASE_OUTPUT_DIR}/{BEHAVIOR_MODE}/group_level/data/
    """
    return get_base_dir() / 'group_level' / 'data'


def get_group_level_output_dir():
    """Get output directory for group-level figures."""
    return get_base_dir() / 'group_level' / 'figures'


# ==============================================================================
#  CONFIGURATION SUMMARY
# ==============================================================================

def print_config_summary():
    """Print a summary of current configuration."""
    animals = get_animals_to_process()
    animal_names = [a['mouse_id'] for a in animals]
    
    print("=" * 70)
    print("PYTHON PLOTTING PIPELINE CONFIGURATION")
    print("=" * 70)
    print(f"  Behavior Mode:     {BEHAVIOR_MODE}")
    print(f"  Base Directory:    {get_base_dir()}")
    print(f"  Methods:           {', '.join(METHODS)}")
    print(f"  Figure Formats:    {', '.join(FIGURE_FORMATS)}")
    print(f"  DPI:               {FIGURE_DPI}")
    print(f"  Frequency Range:   {FREQ_MIN}-{FREQ_MAX} Hz")
    print("-" * 70)
    print(f"  Animals to Plot:   {', '.join(animal_names)}")
    print("-" * 70)
    print("  Analysis Levels:")
    print(f"    Single-Trial:        {PLOT_SINGLE_TRIAL}")
    print(f"    Session-Pooled:      {PLOT_SESSION_POOLED}")
    print(f"    Animal-Pooled:       {PLOT_ANIMAL_POOLED}")
    print(f"    Animal-Concatenated: {PLOT_ANIMAL_CONCATENATED}")
    print(f"    Group-Level:         {PLOT_GROUP_LEVEL}")
    print(f"    Group Pooling Level: {GROUP_POOLING_LEVEL}")
    print("-" * 70)
    print("  Group-Level Options:")
    print(f"    Primary Methods:         {', '.join(GROUP_LEVEL_METHODS)}")
    print(f"    Include Non-Cluster:     {INCLUDE_NON_CLUSTER_METHODS}")
    print(f"    Generate Band Boxplots:  {GENERATE_BAND_BOXPLOTS}")
    print(f"    Show Session-Pooled:     {SHOW_SESSION_POOLED}")
    print("=" * 70)


# Create a CONFIG object for easy import
CONFIG = {
    'behavior_mode': BEHAVIOR_MODE,
    'animals_to_process': ANIMALS_TO_PROCESS,
    'plot_single_trial': PLOT_SINGLE_TRIAL,
    'plot_session_pooled': PLOT_SESSION_POOLED,
    'plot_animal_pooled': PLOT_ANIMAL_POOLED,
    'plot_animal_concatenated': PLOT_ANIMAL_CONCATENATED,
    'plot_group_level': PLOT_GROUP_LEVEL,
    'methods': METHODS,
    'figure_formats': FIGURE_FORMATS,
    'figure_dpi': FIGURE_DPI,
    'freq_min': FREQ_MIN,
    'freq_max': FREQ_MAX,
    'coh_vmin': COH_VMIN,
    'coh_vmax': COH_VMAX,
    # Group-level options
    'group_level_methods': GROUP_LEVEL_METHODS,
    'include_non_cluster_methods': INCLUDE_NON_CLUSTER_METHODS,
    'generate_band_boxplots': GENERATE_BAND_BOXPLOTS,
    'show_session_pooled': SHOW_SESSION_POOLED,
    'group_pooling_level': GROUP_POOLING_LEVEL,
}


if __name__ == '__main__':
    print_config_summary()
    print("\nAnimal Database:")
    for animal in ANIMALS:
        print(f"\n  {animal['mouse_id']} ({animal['project']}):")
        for sess in animal['sessions']:
            print(f"    - {sess['session_id']}: {sess['num_trials']} trials")
