"""
================================================================================
MASTER PLOTTING PIPELINE - Automated Figure Generation
================================================================================

This script automates the generation of all publication-quality figures at
all analysis levels:
  - Single-Trial:     Per-trial spectrograms, coherence heatmaps, spectra
  - Session-Pooled:   Coherence and PSD spectra per session
  - Animal-Pooled:    Coherence and PSD spectra per animal
  - Group-Level:      Group statistics with REST vs RUN comparisons

USAGE:
------
1. Edit plotting_config.py to set:
   - BEHAVIOR_MODE: 'standard' or 'clear'
   - ANIMALS_TO_PROCESS: None (all) or ['Animal1', 'Animal2']
   - PLOT_*: True/False for each analysis level

2. Run this script:
   python run_all_plots.py
================================================================================
"""

import sys
import os
from pathlib import Path
import warnings
import importlib.util

# Add current directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))
# common.py / plotting_config.py live in ../common/ (shared across all figures)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))

# Import configuration
from plotting_config import (
    CONFIG, ANIMALS, METHODS,
    get_animals_to_process, get_session_ids, get_trial_count,
    get_single_trial_input_dir, get_single_trial_output_dir,
    get_session_pooled_input_dir, get_session_pooled_output_dir,
    get_animal_pooled_input_dir, get_animal_pooled_output_dir,
    get_animal_concatenated_input_dir, get_animal_concatenated_output_dir,
    get_group_level_input_dir, get_group_level_output_dir,
    print_config_summary, FIGURE_FORMATS, FIGURE_DPI,
    PLOT_SINGLE_TRIAL, PLOT_SESSION_POOLED, PLOT_ANIMAL_POOLED, 
    PLOT_ANIMAL_CONCATENATED, PLOT_GROUP_LEVEL,
)

# Import optional figure flags (with defaults)
try:
    from plotting_config import PLOT_PSD_LOGLOG
except ImportError:
    PLOT_PSD_LOGLOG = False

# Import group-level specific options (with defaults if not present)
try:
    from plotting_config import (
        GROUP_LEVEL_METHODS, INCLUDE_NON_CLUSTER_METHODS, GENERATE_BAND_BOXPLOTS,
        SHOW_SESSION_POOLED, GROUP_POOLING_LEVEL
    )
except ImportError:
    GROUP_LEVEL_METHODS = ['fieldtrip_cluster']
    INCLUDE_NON_CLUSTER_METHODS = False
    GENERATE_BAND_BOXPLOTS = True
    SHOW_SESSION_POOLED = False
    GROUP_POOLING_LEVEL = 'animal_concatenated'

warnings.filterwarnings('ignore', category=UserWarning)


# ==============================================================================
#  IMPORT PLOTTING MODULES
# ==============================================================================

def import_plotting_modules():
    """Import the individual plotting modules."""
    modules = {}
    
    script_dir = Path(__file__).parent
    
    # Single-trial plotting
    try:
        spec = importlib.util.spec_from_file_location(
            "plot_single_trial", 
            script_dir / "fig2_coherence.py"
        )
        modules['single_trial'] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(modules['single_trial'])
    except Exception as e:
        print(f"Warning: Could not import single-trial module: {e}")
        modules['single_trial'] = None
    
    # Pooled plotting
    try:
        spec = importlib.util.spec_from_file_location(
            "plot_pooled",
            script_dir / "fig2_coherence_pooled.py"
        )
        modules['pooled'] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(modules['pooled'])
    except Exception as e:
        print(f"Warning: Could not import pooled module: {e}")
        modules['pooled'] = None
    
    # Group-level plotting
    try:
        spec = importlib.util.spec_from_file_location(
            "plot_group",
            script_dir / "fig3_coherence_group.py"
        )
        modules['group'] = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(modules['group'])
    except Exception as e:
        print(f"Warning: Could not import group-level module: {e}")
        modules['group'] = None
    
    return modules


from common import to_long_path  # shared helper (was a local copy)


def path_exists(path):
    """Check if path exists, using long path format for Windows."""
    return os.path.exists(to_long_path(str(path)))


def ensure_directory_exists(dir_path):
    """
    Create directory with robust handling for network paths and long paths.
    Uses Windows extended-length path prefix to bypass MAX_PATH limit.
    Handles edge case where a file (not directory) exists at the path.
    """
    from pathlib import Path
    dir_path = Path(dir_path)
    dir_path_str = str(dir_path)
    long_path = to_long_path(dir_path_str)
    
    # Check if directory already exists
    if dir_path.exists() and dir_path.is_dir():
        return True
    
    # Edge case: if a file exists with this name, remove it
    if dir_path.exists() and dir_path.is_file():
        try:
            dir_path.unlink()
            print(f"      Removed file blocking directory: {dir_path_str}")
        except Exception as e:
            print(f"      ERROR: Could not remove file at directory path: {e}")
            print(f"      Path: {dir_path_str}")
            return False
    
    # Try to create with long path first (for network paths > 260 chars)
    try:
        os.makedirs(long_path, exist_ok=True)
        if os.path.exists(long_path) and os.path.isdir(long_path):
            return True
    except Exception:
        pass
    
    # Fallback to regular path or Path.mkdir
    try:
        dir_path.mkdir(parents=True, exist_ok=True)
        return True
    except Exception as e:
        print(f"      ERROR: Could not create directory: {e}")
        print(f"      Path: {dir_path_str}")
        return False


# ==============================================================================
#  SINGLE-TRIAL PLOTTING
# ==============================================================================

def run_single_trial_plotting(modules):
    """Run single-trial plotting for all selected animals."""
    print("\n" + "=" * 70)
    print("  LEVEL 1: SINGLE-TRIAL ANALYSIS")
    print("=" * 70)
    
    if modules.get('single_trial') is None:
        print("  Single-trial module not available, skipping...")
        return 0
    
    mod = modules['single_trial']
    animals = get_animals_to_process()
    total_plots = 0
    
    for animal in animals:
        mouse_id = animal['mouse_id']
        print(f"\n  Animal: {mouse_id}")
        print("-" * 50)
        
        for session in animal['sessions']:
            session_id = session['session_id']
            num_trials = session['num_trials']
            
            print(f"\n    Session: {session_id} ({num_trials} trials)")
            
            input_dir = get_single_trial_input_dir(mouse_id, session_id)
            output_dir = get_single_trial_output_dir(mouse_id, session_id)
            
            if not path_exists(input_dir):
                print(f"      Input not found: {input_dir}")
                continue
            
            # Create output directory for figures
            print(f"      Creating output dir: {output_dir}")
            if not ensure_directory_exists(output_dir):
                print(f"      ERROR: Failed to create output directory, skipping session")
                continue
            print(f"      Output dir exists: {path_exists(output_dir)}")
            
            for method in METHODS:
                print(f"      Method: {method}")
                
                # Find trial files
                trial_files = sorted(input_dir.glob(f"figure2_{method}_trial*.mat"))
                
                if not trial_files:
                    print(f"        No trial files found")
                    continue
                
                print(f"        Found {len(trial_files)} trial files")
                
                # Load trial data
                data_list = []
                trial_labels = []
                for i, tf in enumerate(trial_files, 1):
                    try:
                        data = mod.load_matlab_data(tf)
                        if data is not None:
                            data_list.append(data)
                            trial_labels.append(f'Trial {i}')
                    except Exception as e:
                        print(f"        Error loading {tf.name}: {e}")
                
                if not data_list:
                    print(f"        No valid data loaded")
                    continue
                
                # Generate figures
                # Use short method names to avoid long path issues
                short_method = 'ft' if method == 'fieldtrip' else 'ms'
                
                try:
                    # Main figure (spectrograms + coherence heatmaps)
                    # Output: {short_method}_heatmaps.{fmt}
                    main_output = output_dir / short_method
                    mod.create_main_figure(data_list, trial_labels, method, main_output)
                    total_plots += 1
                    
                    # Coherence spectrum figure
                    # Output: {short_method}_coherence.{fmt}  (function appends _coherence)
                    coh_output = output_dir / short_method
                    mod.create_coherence_spectrum_figure(data_list, trial_labels, method, coh_output)
                    total_plots += 1
                    
                    # PSD figure
                    # Output: {short_method}_psd.{fmt}  (function appends _psd)
                    psd_output = output_dir / short_method
                    mod.create_psd_figure(data_list, trial_labels, method, psd_output)
                    total_plots += 1
                    
                    print(f"        Generated 3 figure sets ({short_method}_*)")
                except Exception as e:
                    import traceback
                    print(f"        Error generating figures: {e}")
                    traceback.print_exc()
    
    return total_plots


# ==============================================================================
#  SESSION-POOLED PLOTTING
# ==============================================================================

def run_session_pooled_plotting(modules):
    """Run session-pooled plotting for all selected animals."""
    print("\n" + "=" * 70)
    print("  LEVEL 2: SESSION-POOLED ANALYSIS")
    print("=" * 70)
    
    if modules.get('pooled') is None:
        print("  Pooled module not available, skipping...")
        return 0
    
    mod = modules['pooled']
    animals = get_animals_to_process()
    total_plots = 0
    
    for animal in animals:
        mouse_id = animal['mouse_id']
        print(f"\n  Animal: {mouse_id}")
        print("-" * 50)
        
        # Build list of session IDs to process:
        # 1. Individual sessions that are NOT part of a pooled group
        # 2. Combined sessions from session_pooled_groups
        sessions_to_plot = []
        
        # Get session IDs that are part of pooled groups
        sessions_in_groups = set()
        session_groups = animal.get('session_pooled_groups', [])
        
        for group in session_groups:
            if len(group) >= 2:
                # Add all sessions in this group to the exclusion set
                for sess_id in group:
                    sessions_in_groups.add(sess_id)
                # Add the combined session name (MATLAB convention: first_session-combined)
                combined_name = f"{group[0]}-combined"
                sessions_to_plot.append({'session_id': combined_name, 'is_combined': True})
        
        # Add individual sessions not in any group
        for session in animal['sessions']:
            if session['session_id'] not in sessions_in_groups:
                sessions_to_plot.append({'session_id': session['session_id'], 'is_combined': False})
        
        for session_info in sessions_to_plot:
            session_id = session_info['session_id']
            is_combined = session_info.get('is_combined', False)
            label = "(combined)" if is_combined else ""
            print(f"\n    Session: {session_id} {label}")
            
            input_dir = get_session_pooled_input_dir(mouse_id, session_id)
            output_dir = get_session_pooled_output_dir(mouse_id, session_id)
            
            if not path_exists(input_dir):
                print(f"      Input not found: {input_dir}")
                continue
            
            ensure_directory_exists(output_dir)
            
            for method in METHODS:
                print(f"      Method: {method}")
                
                input_file = input_dir / f"{method}.mat"
                
                if not path_exists(input_file):
                    print(f"        File not found: {input_file.name}")
                    continue
                
                # Load data
                data = mod.load_matlab_data(input_file)
                if data is None:
                    print(f"        Failed to load data")
                    continue
                
                title_prefix = f'{mouse_id} – {session_id} (Session Pooled)'
                
                try:
                    # Coherence spectrum
                    coh_output = output_dir / f'{method}_coherence'
                    mod.plot_coherence_spectrum(data, method, title_prefix, coh_output)
                    
                    # PSD spectrum
                    psd_output = output_dir / f'{method}_psd'
                    mod.plot_psd_spectrum(data, method, title_prefix, psd_output)
                    
                    total_plots += 2
                    print(f"        Generated 2 figures")
                except Exception as e:
                    print(f"        Error: {e}")
    
    return total_plots


# ==============================================================================
#  ANIMAL-POOLED PLOTTING
# ==============================================================================

def run_animal_pooled_plotting(modules):
    """Run animal-pooled plotting for all selected animals."""
    print("\n" + "=" * 70)
    print("  LEVEL 3: ANIMAL-POOLED ANALYSIS")
    print("=" * 70)
    
    if modules.get('pooled') is None:
        print("  Pooled module not available, skipping...")
        return 0
    
    mod = modules['pooled']
    animals = get_animals_to_process()
    total_plots = 0
    
    for animal in animals:
        mouse_id = animal['mouse_id']
        print(f"\n  Animal: {mouse_id}")
        print("-" * 50)
        
        input_dir = get_animal_pooled_input_dir(mouse_id)
        output_dir = get_animal_pooled_output_dir(mouse_id)
        
        if not path_exists(input_dir):
            print(f"    Input not found: {input_dir}")
            continue
        
        ensure_directory_exists(output_dir)
        
        for method in METHODS:
            print(f"    Method: {method}")
            
            input_file = input_dir / f"{method}.mat"
            
            if not path_exists(input_file):
                print(f"      File not found: {input_file.name}")
                continue
            
            # Load data
            data = mod.load_matlab_data(input_file)
            if data is None:
                print(f"      Failed to load data")
                continue
            
            title_prefix = f'{mouse_id} (Animal Pooled – All Sessions)'
            
            try:
                # Coherence spectrum
                coh_output = output_dir / f'{method}_coherence'
                mod.plot_coherence_spectrum(data, method, title_prefix, coh_output)
                
                # PSD spectrum
                psd_output = output_dir / f'{method}_psd'
                mod.plot_psd_spectrum(data, method, title_prefix, psd_output)
                
                total_plots += 2
                print(f"      Generated 2 figures")
            except Exception as e:
                print(f"      Error: {e}")
    
    return total_plots


# ==============================================================================
#  ANIMAL-CONCATENATED PLOTTING
# ==============================================================================

def run_animal_concatenated_plotting(modules):
    """Run animal-concatenated plotting for all selected animals.
    
    Animal-concatenated: All raw data from all sessions concatenated, spectra computed once.
    This differs from animal-pooled which averages spectra across sessions.
    """
    print("\n" + "=" * 70)
    print("  LEVEL 3b: ANIMAL-CONCATENATED ANALYSIS")
    print("=" * 70)
    
    if modules.get('pooled') is None:
        print("  Pooled module not available, skipping...")
        return 0
    
    mod = modules['pooled']
    animals = get_animals_to_process()
    total_plots = 0
    
    for animal in animals:
        mouse_id = animal['mouse_id']
        print(f"\n  Animal: {mouse_id}")
        print("-" * 50)
        
        input_dir = get_animal_concatenated_input_dir(mouse_id)
        output_dir = get_animal_concatenated_output_dir(mouse_id)
        
        if not path_exists(input_dir):
            print(f"    Input not found: {input_dir}")
            continue
        
        ensure_directory_exists(output_dir)
        
        for method in METHODS:
            print(f"    Method: {method}")
            
            input_file = input_dir / f"{method}.mat"
            
            if not path_exists(input_file):
                print(f"      File not found: {input_file.name}")
                continue
            
            # Load data
            data = mod.load_matlab_data(input_file)
            if data is None:
                print(f"      Failed to load data")
                continue
            
            title_prefix = f'{mouse_id} (Animal Concatenated – All Sessions)'
            
            try:
                # Coherence spectrum
                coh_output = output_dir / f'{method}_coherence'
                mod.plot_coherence_spectrum(data, method, title_prefix, coh_output)

                # Additional coherence spectrum with log-frequency x-axis
                coh_log_output = output_dir / f'{method}_coherence_logfreq'
                if hasattr(mod, 'plot_coherence_spectrum_logfreq'):
                    mod.plot_coherence_spectrum_logfreq(data, method, title_prefix, coh_log_output)
                
                # PSD spectrum
                psd_output = output_dir / f'{method}_psd'
                mod.plot_psd_spectrum(data, method, title_prefix, psd_output)
                
                n_plots = 3
                
                # Optional: Log-log PSD spectrum (for 1/f analysis)
                if PLOT_PSD_LOGLOG:
                    psd_log_output = output_dir / f'{method}_psd_log'
                    if hasattr(mod, 'plot_psd_spectrum_loglog'):
                        mod.plot_psd_spectrum_loglog(data, method, title_prefix, psd_log_output)
                        n_plots += 1
                
                total_plots += n_plots
                print(f"      Generated {n_plots} figures")
            except Exception as e:
                print(f"      Error: {e}")
    
    return total_plots


# ==============================================================================
#  GROUP-LEVEL PLOTTING
# ==============================================================================

def run_group_level_plotting(modules):
    """Run group-level plotting.
    
    Group-level statistics explanation:
    ----------------------------------
    - SESSION-POOLED: Each session is an observation (N = total sessions across animals)
      WARNING: This is PSEUDOREPLICATION - sessions from same animal are not independent!
      Use for exploratory analysis / showing session variability only.
      
    - ANIMAL-POOLED: Each animal is an observation (N = number of animals)
      This is the STATISTICALLY CORRECT approach for publication.
      The independent unit of replication is the animal.
    
    Both are plotted side-by-side for comparison, but ANIMAL-POOLED should be
    the primary result reported in publications.
    """
    print("\n" + "=" * 70)
    print("  LEVEL 4: GROUP-LEVEL STATISTICS")
    print("=" * 70)
    
    if modules.get('group') is None:
        print("  Group-level module not available, skipping...")
        return 0
    
    mod = modules['group']
    
    input_dir = get_group_level_input_dir()
    output_dir = get_group_level_output_dir()
    
    if not path_exists(input_dir):
        print(f"  Input directory not found: {input_dir}")
        print("  (Group-level data must be generated from MATLAB first)")
        return 0
    
    ensure_directory_exists(output_dir)
    
    print(f"  Input:  {input_dir}")
    print(f"  Output: {output_dir}")
    print(f"\n  Configuration:")
    print(f"    Primary methods: {GROUP_LEVEL_METHODS}")
    print(f"    Include non-cluster methods: {INCLUDE_NON_CLUSTER_METHODS}")
    print(f"    Generate band boxplots: {GENERATE_BAND_BOXPLOTS}")
    print(f"    Show session-pooled: {SHOW_SESSION_POOLED}")
    
    total_plots = 0
    
    # Build list of methods to process based on configuration
    # Primary method is always FieldTrip Cluster (cluster-based permutation statistics)
    methods_config = []
    
    # Always include configured primary methods
    if 'fieldtrip_cluster' in GROUP_LEVEL_METHODS:
        methods_config.append(('fieldtrip_cluster', 'FieldTrip Cluster', 'cluster'))
    
    # Optionally include non-cluster methods
    if INCLUDE_NON_CLUSTER_METHODS:
        if 'mscohere' in GROUP_LEVEL_METHODS or INCLUDE_NON_CLUSTER_METHODS:
            methods_config.append(('mscohere', 'mscohere', 'standard'))
        if 'fieldtrip' in GROUP_LEVEL_METHODS or INCLUDE_NON_CLUSTER_METHODS:
            methods_config.append(('fieldtrip', 'FieldTrip', 'standard'))
    
    # Get human-readable label for pooling level
    pooling_label = GROUP_POOLING_LEVEL.replace('_', '-')
    
    for file_prefix, method_title, sig_source in methods_config:
        print(f"\n    Method: {method_title}")
        
        # Load data - MATLAB naming convention: {method}_{level}.mat
        # Uses configurable GROUP_POOLING_LEVEL for animal-level data
        session_path = input_dir / f"{file_prefix}_session_pooled.mat"
        animal_path = input_dir / f"{file_prefix}_{GROUP_POOLING_LEVEL}.mat"
        
        # Only load session data if user wants to show it
        data_session = None
        if SHOW_SESSION_POOLED and path_exists(session_path):
            data_session = mod.load_matlab_struct(session_path)
        
        data_animal = mod.load_matlab_struct(animal_path) if path_exists(animal_path) else None
        
        if data_animal is None:
            print(f"      No {pooling_label} data found, skipping...")
            continue
        
        if data_session:
            n_sess = data_session.get('num_sessions', '?')
            print(f"      Loaded session-pooled data (N={n_sess} sessions)")
        
        n_anim = data_animal.get('num_animals', '?')
        print(f"      Loaded {pooling_label} data (N={n_anim} animals) [PRIMARY]")
        
        try:
            # Coherence figure - WITH significance markers
            coh_path = output_dir / f"figure3_{file_prefix}_coherence"
            mod.create_coherence_figure(data_session, data_animal, coh_path, method_title, sig_source, show_significance=True)
            total_plots += 1
            
            # Coherence figure - WITHOUT significance markers
            coh_path_nosig = output_dir / f"figure3_{file_prefix}_coherence_no_sig"
            mod.create_coherence_figure(data_session, data_animal, coh_path_nosig, method_title, sig_source, show_significance=False)
            total_plots += 1

            # Additional log-frequency coherence figure - WITH significance
            coh_log_path = output_dir / f"figure3_{file_prefix}_coherence_logfreq"
            if hasattr(mod, 'create_coherence_figure_logfreq'):
                mod.create_coherence_figure_logfreq(
                    data_session, data_animal, coh_log_path, method_title, sig_source, show_significance=True
                )
                total_plots += 1

            # Additional log-frequency coherence figure - WITHOUT significance
            coh_log_path_nosig = output_dir / f"figure3_{file_prefix}_coherence_logfreq_no_sig"
            if hasattr(mod, 'create_coherence_figure_logfreq'):
                mod.create_coherence_figure_logfreq(
                    data_session, data_animal, coh_log_path_nosig, method_title, sig_source, show_significance=False
                )
                total_plots += 1
            
            # PSD figure - WITH significance markers
            psd_path = output_dir / f"figure3_{file_prefix}_psd"
            mod.create_psd_figure(data_session, data_animal, psd_path, method_title, sig_source, show_significance=True)
            total_plots += 1
            
            # PSD figure - WITHOUT significance markers
            psd_path_nosig = output_dir / f"figure3_{file_prefix}_psd_no_sig"
            mod.create_psd_figure(data_session, data_animal, psd_path_nosig, method_title, sig_source, show_significance=False)
            total_plots += 1
            
            # Optional: Log-log PSD figure (for 1/f analysis)
            n_extra = 0
            if PLOT_PSD_LOGLOG and hasattr(mod, 'create_psd_figure_loglog'):
                psd_log_path = output_dir / f"figure3_{file_prefix}_psd_log"
                mod.create_psd_figure_loglog(data_session, data_animal, psd_log_path, method_title, sig_source, show_significance=False)
                total_plots += 1
                n_extra += 1
            
            # Band-averaged box plot
            # - Always for cluster method if GENERATE_BAND_BOXPLOTS is True
            # - Also for non-cluster methods if INCLUDE_NON_CLUSTER_METHODS is True
            if GENERATE_BAND_BOXPLOTS:
                band_path = output_dir / f"figure3_{file_prefix}_band_boxplot"
                mod.create_band_boxplot_figure(data_session, data_animal, band_path, method_title)
                total_plots += 1
                print(f"      Generated {7 + n_extra} figures (coh+coh_logfreq/psd with/without sig + band boxplot{' + psd_log' if n_extra else ''})")
            else:
                print(f"      Generated {6 + n_extra} figures (coh+coh_logfreq/psd with/without significance{' + psd_log' if n_extra else ''})")
                
        except Exception as e:
            print(f"      Error: {e}")
    
    # Run FOOOF analysis
    # NOTE: FOOOF analyzes PSD (Power Spectral Density), NOT coherence.
    # Since PSD is computed identically for both mscohere and FieldTrip methods
    # (both use pwelch), we only need to run FOOOF once using mscohere data.
    try:
        from fig4_fooof import run_fooof_pipeline, FOOOF_AVAILABLE, get_pooling_level_label
        
        if FOOOF_AVAILABLE:
            pooling_label = get_pooling_level_label(GROUP_POOLING_LEVEL)
            print(f"\n    Running FOOOF Analysis ({pooling_label})...")
            print("    (FOOOF analyzes PSD which is identical for both methods)")
            animals = get_animals_to_process()
            
            # Run FOOOF on mscohere-derived PSD (same as fieldtrip PSD since both use pwelch)
            # Uses GROUP_POOLING_LEVEL from config by default
            fooof_plots = run_fooof_pipeline(animals, output_dir, method='mscohere', 
                                             pooling_level=GROUP_POOLING_LEVEL)
            total_plots += fooof_plots
            
            if fooof_plots > 0:
                print(f"\n      Generated {fooof_plots} FOOOF figure(s)")
        else:
            print("\n    Skipping FOOOF (not installed - pip install fooof)")
    except ImportError as e:
        print(f"\n    Skipping FOOOF (import error: {e})")
    except Exception as e:
        import traceback
        print(f"\n    FOOOF Error: {e}")
        traceback.print_exc()
    
    return total_plots


# ==============================================================================
#  MAIN EXECUTION
# ==============================================================================

def main():
    """Main execution function."""
    
    # Print configuration summary
    print_config_summary()
    
    # Import plotting modules
    print("\nLoading plotting modules...")
    modules = import_plotting_modules()
    
    # Track statistics
    stats = {
        'single_trial': 0,
        'session_pooled': 0,
        'animal_pooled': 0,
        'animal_concatenated': 0,
        'group_level': 0,
    }
    
    # Run each level based on configuration
    if PLOT_SINGLE_TRIAL:
        stats['single_trial'] = run_single_trial_plotting(modules)
    else:
        print("\n  Skipping Single-Trial (disabled in config)")
    
    if PLOT_SESSION_POOLED:
        stats['session_pooled'] = run_session_pooled_plotting(modules)
    else:
        print("\n  Skipping Session-Pooled (disabled in config)")
    
    if PLOT_ANIMAL_POOLED:
        stats['animal_pooled'] = run_animal_pooled_plotting(modules)
    else:
        print("\n  Skipping Animal-Pooled (disabled in config)")
    
    if PLOT_ANIMAL_CONCATENATED:
        stats['animal_concatenated'] = run_animal_concatenated_plotting(modules)
    else:
        print("\n  Skipping Animal-Concatenated (disabled in config)")
    
    if PLOT_GROUP_LEVEL:
        stats['group_level'] = run_group_level_plotting(modules)
    else:
        print("\n  Skipping Group-Level (disabled in config)")
    
    # Print summary
    total = sum(stats.values())
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    print(f"  Single-Trial figures:       {stats['single_trial']}")
    print(f"  Session-Pooled figures:     {stats['session_pooled']}")
    print(f"  Animal-Pooled figures:      {stats['animal_pooled']}")
    print(f"  Animal-Concatenated figures: {stats['animal_concatenated']}")
    print(f"  Group-Level figures:        {stats['group_level']}")
    print("-" * 70)
    print(f"  TOTAL FIGURES GENERATED:    {total}")
    print("=" * 70)


if __name__ == '__main__':
    main()
