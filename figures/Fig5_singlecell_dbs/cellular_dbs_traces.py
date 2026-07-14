"""
Plot Cellular Voltage Imaging DBS Traces

Creates two types of figures for each session:
1. All-trial averaged: Stimulation, LFP (avg), all-neuron avg trace, individual neuron traces (avg)
   - SEM shading across trials
2. Single-trial example: Same layout but for one trial
   - SEM shading across neurons for averaged trace
   - Spike markers for individual neuron traces

Sessions:
- Animal01 01-06-25-R1: 40Hz AmpBalanced
- Animal01 01-06-25-R2: 40Hz EnergyBalanced
- Animal01 01-06-25-R10: 130Hz EnergyBalanced
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.io import loadmat
from scipy import stats
from pathlib import Path
import warnings
import sys
import h5py

warnings.filterwarnings('ignore')

# Locate config/paths_config.py by walking up from this file.
for _d in Path(__file__).resolve().parents:
    if (_d / "config" / "paths_config.py").exists():
        sys.path.insert(0, str(_d / "config"))
        break
from paths_config import DATA_ROOT as _LAB_DATA_ROOT, PROJECT_ROOT

# =============================================================================
# CONFIGURATION  -- EDIT THESE FOR YOUR OWN DATASET
# =============================================================================

# Base path to cellular processed data
BASE_PATH = _LAB_DATA_ROOT / "DBS" / "Animal01" / "BaselineDBS" / "CellularDataProcessed"

# Output folder for figures
OUTPUT_FOLDER = PROJECT_ROOT / "Figures" / "Cellular_DBS_traces"

# Session configuration
# Note: 05-06-25 R1 and R2 are low quality, using 04-06-25 R3 and R13 instead
SESSIONS = [
    # 05-06-25 sessions (original, lower quality)
    {'date': '05-06-25', 'rec_id': 'R1', 'freq_hz': 40, 'comparison': 'AmpBalanced', 'label': '40Hz Amp-Balanced (05-06)'},
    {'date': '05-06-25', 'rec_id': 'R2', 'freq_hz': 40, 'comparison': 'EnergyBalanced', 'label': '40Hz Energy-Balanced (05-06)'},
    {'date': '05-06-25', 'rec_id': 'R10', 'freq_hz': 130, 'comparison': 'EnergyBalanced', 'label': '130Hz Energy-Balanced (05-06)'},
    # 04-06-25 sessions (better quality replacements for 05-06-25)
    {'date': '04-06-25', 'rec_id': 'R3', 'freq_hz': 40, 'comparison': 'AmpBalanced', 'label': '40Hz Amp-Balanced (04-06)'},
    {'date': '04-06-25', 'rec_id': 'R13', 'freq_hz': 40, 'comparison': 'EnergyBalanced', 'label': '40Hz Energy-Balanced (04-06)'},
    # 18-06-25 sessions (alternative set with more trials)
    {'date': '18-06-25', 'rec_id': 'R2', 'freq_hz': 40, 'comparison': 'AmpBalanced', 'label': '40Hz Amp-Balanced (18-06)'},
    {'date': '18-06-25', 'rec_id': 'R4', 'freq_hz': 130, 'comparison': 'AmpBalanced', 'label': '130Hz Amp-Balanced (18-06)'},
    {'date': '18-06-25', 'rec_id': 'R5', 'freq_hz': 130, 'comparison': 'EnergyBalanced', 'label': '130Hz Energy-Balanced (18-06)'},
]

# Set which sessions to process (comment/uncomment as needed)
# Examples:
# SESSIONS_TO_PROCESS = ['04-06-25-R3', '04-06-25-R13', '05-06-25-R10']  # Original good quality set
# SESSIONS_TO_PROCESS = ['18-06-25-R2', '18-06-25-R4', '18-06-25-R5']    # Alternative set (18-06)
SESSIONS_TO_PROCESS = None  # None = process all sessions

MOUSE_NAME = 'Animal01'

# Time window configuration (relative to stim onset at 0)
PRE_STIM_SEC = 1.0    # Show 1 second before stim onset
STIM_DURATION_SEC = 1.0  # Stimulation lasts 1 second
POST_STIM_SEC = 1.0   # Show 1 second after stim end
# Total window: -1 to 2 seconds

# Representative trial for single-trial plots
REPRESENTATIVE_TRIAL = 1  # Trial 1 (0-indexed in Python)

# Figure configuration
DPI = 300
FONT_SIZE_TITLE = 16
FONT_SIZE_LABEL = 14
FONT_SIZE_TICK = 12
FONT_SIZE_LEGEND = 10
FONT_SIZE_SCALEBAR = 11

# Line widths
LINE_WIDTH_TRACE = 1.0
LINE_WIDTH_THICK = 1.5
LINE_WIDTH_AVG = 2.0
AXIS_LINEWIDTH = 1.5
SCALEBAR_LINEWIDTH = 3.0

# Colors (fiber pipeline style)
COLOR_LFP = np.array([0.35, 0.25, 0.45])  # purple-grey
COLOR_NEURON_AVG = np.array([0.08, 0.45, 0.45])  # dark teal for averaged
COLOR_NEURON_IND = np.array([0.127568, 0.566949, 0.550556])  # teal for individual
COLOR_SPIKE = np.array([0.8, 0.1, 0.1])  # red for spike markers
COLOR_STIM_PULSE = np.array([0.4, 0.1, 0.1])  # dark red for stim pulses
COLOR_SEM = 0.3  # alpha for SEM shading
COLOR_STIM_SHADE = np.array([1.0, 0.9, 0.9])  # light red for stim period background

# Trace offsets for stacking
TRACE_SPACING = 1.0  # Spacing between individual neuron traces (normalized units)
TRACE_SCALE = 5.0    # Scale factor to amplify small fluorescence signals for visibility


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def get_hdf5_value(f, dataset):
    """
    Safely extract value from HDF5 dataset, handling references.
    """
    data = dataset[()]
    # If it's a reference, dereference it
    if isinstance(data, np.ndarray) and data.dtype == h5py.special_dtype(ref=h5py.Reference):
        if data.size == 1:
            return f[data.flat[0]][()]
        else:
            return np.array([f[ref][()] for ref in data.flat])
    return data


def load_cellular_data(session):
    """
    Load cellular analysis data for a session.
    
    Returns dict with trials data or None if loading fails.
    """
    date = session['date']
    rec_id = session['rec_id']
    
    session_folder = BASE_PATH / f"{date}-{rec_id}"
    mat_file = session_folder / f"{MOUSE_NAME}_{date}-{rec_id}_CellularAnalysis.mat"
    
    if not mat_file.exists():
        print(f"  WARNING: File not found: {mat_file}")
        return None
    
    print(f"  Loading: {mat_file}")
    
    try:
        # Load with h5py (MATLAB v7.3)
        with h5py.File(str(mat_file), 'r') as f:
            cellular = f['CellularAnalysis']
            
            # Get metadata - direct access like other scripts do
            num_trials = int(np.array(cellular['metadata']['num_trials']).flat[0])
            num_neurons = int(np.array(cellular['metadata']['num_neurons']).flat[0])
            
            print(f"    Found {num_trials} trials, {num_neurons} neurons")
            
            # Get trials references - this is a cell array stored as references
            trials_refs = cellular['trials'][()]
            
            trials_data = []
            for trial_idx in range(num_trials):
                # Dereference to get the actual trial group
                trial_ref = trials_refs.flat[trial_idx]
                trial_grp = f[trial_ref]
                
                trial_dict = {}
                
                # Time vector
                try:
                    time_grp = trial_grp['time']
                    time_vec_ref = time_grp['time_vector'][()]
                    if isinstance(time_vec_ref, np.ndarray) and time_vec_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        trial_dict['time_vector'] = f[time_vec_ref.flat[0]][()].flatten()
                    else:
                        trial_dict['time_vector'] = np.array(time_vec_ref).flatten()
                except Exception as e:
                    print(f"      Warning: Could not load time_vector for trial {trial_idx+1}: {e}")
                    trial_dict['time_vector'] = np.array([])
                
                # Stimulus onset frame
                try:
                    stim_onset = time_grp['stimulus_onset_frame'][()]
                    if isinstance(stim_onset, np.ndarray) and stim_onset.dtype == h5py.special_dtype(ref=h5py.Reference):
                        trial_dict['stim_onset_frame'] = int(f[stim_onset.flat[0]][()].flat[0])
                    else:
                        trial_dict['stim_onset_frame'] = int(np.array(stim_onset).flat[0])
                except:
                    trial_dict['stim_onset_frame'] = None
                
                # Sampling rate
                try:
                    params_grp = trial_grp['parameters']
                    fs_data = params_grp['imaging_fs'][()]
                    if isinstance(fs_data, np.ndarray) and fs_data.dtype == h5py.special_dtype(ref=h5py.Reference):
                        trial_dict['fs'] = float(f[fs_data.flat[0]][()].flat[0])
                    else:
                        trial_dict['fs'] = float(np.array(fs_data).flat[0])
                except:
                    trial_dict['fs'] = 1000.0
                
                # Fluorescence traces
                try:
                    signals_grp = trial_grp['signals']
                    fluor_ref = signals_grp['fluorescence_corrected'][()]
                    if isinstance(fluor_ref, np.ndarray) and fluor_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        fluor_data = f[fluor_ref.flat[0]][()]
                    else:
                        fluor_data = np.array(fluor_ref)
                    
                    # Ensure frames x neurons orientation
                    if fluor_data.ndim == 2 and fluor_data.shape[0] < fluor_data.shape[1]:
                        fluor_data = fluor_data.T
                    trial_dict['fluorescence'] = fluor_data
                except Exception as e:
                    print(f"      Warning: Could not load fluorescence for trial {trial_idx+1}: {e}")
                    trial_dict['fluorescence'] = np.array([])
                
                # LFP trace
                try:
                    ephys_grp = trial_grp['ephys']
                    lfp_ref = ephys_grp['lfp_raw_aligned'][()]
                    if isinstance(lfp_ref, np.ndarray) and lfp_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        trial_dict['lfp'] = f[lfp_ref.flat[0]][()].flatten()
                    else:
                        trial_dict['lfp'] = np.array(lfp_ref).flatten()
                except:
                    trial_dict['lfp'] = None
                
                # Stimulus trace
                try:
                    stim_ref = ephys_grp['stimulus_aligned'][()]
                    if isinstance(stim_ref, np.ndarray) and stim_ref.dtype == h5py.special_dtype(ref=h5py.Reference):
                        trial_dict['stimulus'] = f[stim_ref.flat[0]][()].flatten()
                    else:
                        trial_dict['stimulus'] = np.array(stim_ref).flatten()
                except:
                    trial_dict['stimulus'] = None
                
                # Spike times per neuron (cell array)
                try:
                    spikes_grp = trial_grp['spikes']
                    spike_times_refs = spikes_grp['spike_times_seconds'][()]
                    spike_times_list = []
                    
                    for n_idx in range(num_neurons):
                        try:
                            st_ref = spike_times_refs.flat[n_idx]
                            st_data = f[st_ref][()]
                            spike_times_list.append(np.array(st_data).flatten())
                        except:
                            spike_times_list.append(np.array([]))
                    
                    trial_dict['spike_times'] = spike_times_list
                except:
                    trial_dict['spike_times'] = [np.array([]) for _ in range(num_neurons)]
                
                trials_data.append(trial_dict)
                print(f"      Trial {trial_idx+1}: {len(trial_dict['time_vector'])} samples, fs={trial_dict['fs']:.1f}Hz")
            
            return {
                'trials': trials_data,
                'num_trials': num_trials,
                'num_neurons': num_neurons,
                'session': session
            }
            
    except Exception as e:
        print(f"    Error loading with h5py: {e}")
        import traceback
        traceback.print_exc()
        return None


def generate_stim_pulses(t, stim_freq_hz, stim_duration_sec):
    """
    Generate square biphasic stimulation pulses.
    
    Parameters:
    -----------
    t : array
        Time vector (relative to stim onset, so 0 = stim start)
    stim_freq_hz : float
        Stimulation frequency in Hz
    stim_duration_sec : float
        Duration of stimulation in seconds
    
    Returns:
    --------
    pulses : array
        Pulse signal (1 = positive phase, -1 = negative phase, 0 = no pulse)
    """
    pulses = np.zeros_like(t)
    
    # Only generate pulses during stimulation period (0 to stim_duration_sec)
    stim_mask = (t >= 0) & (t < stim_duration_sec)
    
    if not np.any(stim_mask):
        return pulses
    
    # Period of stimulation
    period = 1.0 / stim_freq_hz
    pulse_width = period * 0.4  # 40% duty cycle per phase
    
    t_stim = t[stim_mask]
    
    # Phase within each cycle
    phase = np.mod(t_stim, period)
    
    # Biphasic pulse: positive phase then negative phase
    pos_phase = phase < pulse_width
    neg_phase = (phase >= period/2) & (phase < period/2 + pulse_width)
    
    pulses_stim = np.zeros_like(t_stim)
    pulses_stim[pos_phase] = 1
    pulses_stim[neg_phase] = -1
    
    pulses[stim_mask] = pulses_stim
    
    return pulses


def compute_sem(data, axis=0):
    """Compute standard error of the mean along axis."""
    n = data.shape[axis]
    return np.std(data, axis=axis, ddof=1) / np.sqrt(n)


def get_time_window_mask(time_vec, t_min, t_max):
    """Get boolean mask for time window."""
    return (time_vec >= t_min) & (time_vec <= t_max)


def add_scale_bar(ax, x_pos, y_pos, value, label, orientation='vertical', color='black'):
    """Add a scale bar to axis."""
    if orientation == 'vertical':
        ax.plot([x_pos, x_pos], [y_pos, y_pos + value], color=color, linewidth=SCALEBAR_LINEWIDTH, solid_capstyle='butt')
        ax.text(x_pos + 0.05, y_pos + value/2, label, fontsize=FONT_SIZE_SCALEBAR, va='center', ha='left')
    else:
        ax.plot([x_pos, x_pos + value], [y_pos, y_pos], color=color, linewidth=SCALEBAR_LINEWIDTH, solid_capstyle='butt')
        ax.text(x_pos + value/2, y_pos - 0.1, label, fontsize=FONT_SIZE_SCALEBAR, va='top', ha='center')


# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

def plot_all_trial_averaged(data, output_path):
    """
    Create all-trial averaged figure.
    
    Layout:
    - Row 0: Stimulation pulses
    - Row 1: LFP trace (trial-averaged with SEM)
    - Row 2: All-neuron averaged trace (trial-averaged with SEM)
    - Rows 3+: Individual neuron traces (trial-averaged with SEM)
    """
    session = data['session']
    num_trials = data['num_trials']
    num_neurons = data['num_neurons']
    trials = data['trials']
    
    print(f"  Creating all-trial averaged figure...")
    
    # Time window
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    
    # Get reference time vector from first trial
    ref_time = trials[0]['time_vector']
    fs = trials[0]['fs']
    
    # Create common time vector for window
    time_mask = get_time_window_mask(ref_time, t_min, t_max)
    t_plot = ref_time[time_mask]
    n_samples = np.sum(time_mask)
    
    # Collect data across trials
    lfp_all = []
    fluor_all = []  # trials x neurons x time
    
    for trial in trials:
        t_trial = trial['time_vector']
        mask = get_time_window_mask(t_trial, t_min, t_max)
        
        if np.sum(mask) != n_samples:
            # Interpolate to common time base if needed
            continue
        
        if trial['lfp'] is not None:
            lfp_all.append(trial['lfp'][mask])
        
        fluor_trial = trial['fluorescence'][mask, :]
        fluor_all.append(fluor_trial)
    
    if len(fluor_all) == 0:
        print("    No valid trials found!")
        return
    
    # Stack arrays
    fluor_all = np.array(fluor_all)  # trials x time x neurons
    if len(lfp_all) > 0:
        lfp_all = np.array(lfp_all)  # trials x time
        lfp_mean = np.mean(lfp_all, axis=0)
        lfp_sem = compute_sem(lfp_all, axis=0)
        has_lfp = True
    else:
        has_lfp = False
    
    # Compute averages and SEM
    # All-neuron average per trial, then average across trials
    neuron_avg_per_trial = np.mean(fluor_all, axis=2)  # trials x time
    neuron_avg_mean = np.mean(neuron_avg_per_trial, axis=0)  # time
    neuron_avg_sem = compute_sem(neuron_avg_per_trial, axis=0)  # time
    
    # Individual neuron averages across trials
    neuron_means = np.mean(fluor_all, axis=0)  # time x neurons
    neuron_sems = compute_sem(fluor_all, axis=0)  # time x neurons
    
    # Create figure
    n_rows = 3 + num_neurons  # stim + lfp + avg + individual neurons
    height_ratios = [0.3, 1, 1] + [1] * num_neurons
    
    fig = plt.figure(figsize=(12, 2 + num_neurons * 0.8))
    gs = GridSpec(n_rows, 1, figure=fig, height_ratios=height_ratios, hspace=0.1)
    
    # Row 0: Stimulation pulses
    ax_stim = fig.add_subplot(gs[0])
    stim_pulses = generate_stim_pulses(t_plot, session['freq_hz'], STIM_DURATION_SEC)
    ax_stim.plot(t_plot, stim_pulses, color=COLOR_STIM_PULSE, linewidth=LINE_WIDTH_TRACE)
    ax_stim.set_ylim(-1.5, 1.5)
    ax_stim.set_xlim(t_min, t_max)
    ax_stim.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_stim.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_stim.set_ylabel('Stim', fontsize=FONT_SIZE_LABEL)
    ax_stim.set_xticks([])
    ax_stim.set_yticks([])
    for spine in ax_stim.spines.values():
        spine.set_visible(False)
    ax_stim.set_title(f"{MOUSE_NAME} {session['date']}-{session['rec_id']}: {session['label']} - Trial Averaged (n={num_trials})",
                      fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    
    # Row 1: LFP trace
    ax_lfp = fig.add_subplot(gs[1], sharex=ax_stim)
    if has_lfp:
        # Z-score LFP for display
        lfp_z = (lfp_mean - np.mean(lfp_mean)) / np.std(lfp_mean)
        lfp_sem_z = lfp_sem / np.std(lfp_mean)
        
        ax_lfp.fill_between(t_plot, lfp_z - lfp_sem_z, lfp_z + lfp_sem_z, 
                           color=COLOR_LFP, alpha=COLOR_SEM)
        ax_lfp.plot(t_plot, lfp_z, color=COLOR_LFP, linewidth=LINE_WIDTH_AVG)
    ax_lfp.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_lfp.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_lfp.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_lfp.set_ylabel('LFP (z)', fontsize=FONT_SIZE_LABEL)
    ax_lfp.set_xticks([])
    ax_lfp.spines['top'].set_visible(False)
    ax_lfp.spines['right'].set_visible(False)
    ax_lfp.spines['bottom'].set_visible(False)
    
    # Row 2: All-neuron averaged trace
    ax_avg = fig.add_subplot(gs[2], sharex=ax_stim)
    ax_avg.fill_between(t_plot, neuron_avg_mean - neuron_avg_sem, neuron_avg_mean + neuron_avg_sem,
                        color=COLOR_NEURON_AVG, alpha=COLOR_SEM)
    ax_avg.plot(t_plot, neuron_avg_mean, color=COLOR_NEURON_AVG, linewidth=LINE_WIDTH_AVG)
    ax_avg.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_avg.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_avg.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_avg.set_ylabel(f'Avg (n={num_neurons})\nΔF/F', fontsize=FONT_SIZE_LABEL)
    ax_avg.set_xticks([])
    ax_avg.spines['top'].set_visible(False)
    ax_avg.spines['right'].set_visible(False)
    ax_avg.spines['bottom'].set_visible(False)
    
    # Rows 3+: Individual neuron traces (stacked)
    ax_neurons = fig.add_subplot(gs[3:])
    
    for n_idx in range(num_neurons):
        offset = (num_neurons - 1 - n_idx) * TRACE_SPACING
        
        # Z-score or scale the trace for visibility, then add offset
        trace_raw = neuron_means[:, n_idx]
        trace_scaled = (trace_raw - np.mean(trace_raw)) * TRACE_SCALE + offset
        sem_scaled = neuron_sems[:, n_idx] * TRACE_SCALE
        
        # SEM shading
        ax_neurons.fill_between(t_plot, trace_scaled - sem_scaled, trace_scaled + sem_scaled,
                               color=COLOR_NEURON_IND, alpha=0.2)
        ax_neurons.plot(t_plot, trace_scaled, color=COLOR_NEURON_IND, linewidth=LINE_WIDTH_TRACE)
        
        # Neuron label on the left
        ax_neurons.text(t_min - 0.05, offset, f'N{n_idx+1}', fontsize=FONT_SIZE_TICK, 
                       va='center', ha='right')
    
    ax_neurons.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_neurons.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_neurons.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_neurons.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL)
    ax_neurons.set_ylabel('Individual Neurons', fontsize=FONT_SIZE_LABEL)
    ax_neurons.set_xlim(t_min, t_max)
    ax_neurons.set_ylim(-0.5, num_neurons * TRACE_SPACING + 0.5)
    ax_neurons.spines['top'].set_visible(False)
    ax_neurons.spines['right'].set_visible(False)
    ax_neurons.spines['left'].set_visible(False)
    ax_neurons.set_yticks([])  # Hide y-axis ticks
    ax_neurons.tick_params(labelsize=FONT_SIZE_TICK)
    
    plt.tight_layout()
    
    # Save figure
    fig_name = f"{MOUSE_NAME}_{session['date']}-{session['rec_id']}_TrialAveraged"
    fig.savefig(output_path / f"{fig_name}.png", dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path / f"{fig_name}.pdf", bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    print(f"    Saved: {fig_name}")


def plot_single_trial(data, output_path, trial_idx=0):
    """
    Create single-trial figure with spike markers.
    
    Layout:
    - Row 0: Stimulation pulses
    - Row 1: LFP trace
    - Row 2: All-neuron averaged trace (with SEM across neurons)
    - Rows 3+: Individual neuron traces with spike markers
    """
    session = data['session']
    num_neurons = data['num_neurons']
    trials = data['trials']
    
    if trial_idx >= len(trials):
        print(f"    Trial {trial_idx+1} not available!")
        return
    
    trial = trials[trial_idx]
    print(f"  Creating single-trial figure (Trial {trial_idx+1})...")
    
    # Time window
    t_min = -PRE_STIM_SEC
    t_max = STIM_DURATION_SEC + POST_STIM_SEC
    
    t_trial = trial['time_vector']
    time_mask = get_time_window_mask(t_trial, t_min, t_max)
    t_plot = t_trial[time_mask]
    
    # Get fluorescence data
    fluor = trial['fluorescence'][time_mask, :]  # time x neurons
    
    # Compute all-neuron average with SEM across neurons
    neuron_avg = np.mean(fluor, axis=1)  # time
    neuron_sem = compute_sem(fluor.T, axis=0)  # SEM across neurons
    
    # Get LFP
    if trial['lfp'] is not None:
        lfp = trial['lfp'][time_mask]
        has_lfp = True
    else:
        has_lfp = False
    
    # Get spike times
    spike_times = trial['spike_times']
    
    # Create figure
    n_rows = 3 + num_neurons
    height_ratios = [0.3, 1, 1] + [1] * num_neurons
    
    fig = plt.figure(figsize=(12, 2 + num_neurons * 0.8))
    gs = GridSpec(n_rows, 1, figure=fig, height_ratios=height_ratios, hspace=0.1)
    
    # Row 0: Stimulation pulses
    ax_stim = fig.add_subplot(gs[0])
    stim_pulses = generate_stim_pulses(t_plot, session['freq_hz'], STIM_DURATION_SEC)
    ax_stim.plot(t_plot, stim_pulses, color=COLOR_STIM_PULSE, linewidth=LINE_WIDTH_TRACE)
    ax_stim.set_ylim(-1.5, 1.5)
    ax_stim.set_xlim(t_min, t_max)
    ax_stim.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_stim.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_stim.set_ylabel('Stim', fontsize=FONT_SIZE_LABEL)
    ax_stim.set_xticks([])
    ax_stim.set_yticks([])
    for spine in ax_stim.spines.values():
        spine.set_visible(False)
    ax_stim.set_title(f"{MOUSE_NAME} {session['date']}-{session['rec_id']}: {session['label']} - Trial {trial_idx+1}",
                      fontsize=FONT_SIZE_TITLE, fontweight='bold', pad=10)
    
    # Row 1: LFP trace
    ax_lfp = fig.add_subplot(gs[1], sharex=ax_stim)
    if has_lfp:
        lfp_z = (lfp - np.mean(lfp)) / np.std(lfp)
        ax_lfp.plot(t_plot, lfp_z, color=COLOR_LFP, linewidth=LINE_WIDTH_TRACE)
    ax_lfp.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_lfp.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_lfp.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_lfp.set_ylabel('LFP (z)', fontsize=FONT_SIZE_LABEL)
    ax_lfp.set_xticks([])
    ax_lfp.spines['top'].set_visible(False)
    ax_lfp.spines['right'].set_visible(False)
    ax_lfp.spines['bottom'].set_visible(False)
    
    # Row 2: All-neuron averaged trace with SEM across neurons
    ax_avg = fig.add_subplot(gs[2], sharex=ax_stim)
    ax_avg.fill_between(t_plot, neuron_avg - neuron_sem, neuron_avg + neuron_sem,
                        color=COLOR_NEURON_AVG, alpha=COLOR_SEM)
    ax_avg.plot(t_plot, neuron_avg, color=COLOR_NEURON_AVG, linewidth=LINE_WIDTH_AVG)
    ax_avg.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_avg.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_avg.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_avg.set_ylabel(f'Avg (n={num_neurons})\nΔF/F', fontsize=FONT_SIZE_LABEL)
    ax_avg.set_xticks([])
    ax_avg.spines['top'].set_visible(False)
    ax_avg.spines['right'].set_visible(False)
    ax_avg.spines['bottom'].set_visible(False)
    
    # Rows 3+: Individual neuron traces with spike markers
    ax_neurons = fig.add_subplot(gs[3:])
    
    for n_idx in range(num_neurons):
        offset = (num_neurons - 1 - n_idx) * TRACE_SPACING
        
        # Scale the trace for visibility, then add offset
        trace_raw = fluor[:, n_idx]
        trace_scaled = (trace_raw - np.mean(trace_raw)) * TRACE_SCALE + offset
        
        # Plot trace
        ax_neurons.plot(t_plot, trace_scaled, color=COLOR_NEURON_IND, linewidth=LINE_WIDTH_TRACE)
        
        # Plot spike markers (vertical lines)
        if n_idx < len(spike_times) and len(spike_times[n_idx]) > 0:
            spikes = spike_times[n_idx]
            # Filter spikes within time window
            spike_mask = (spikes >= t_min) & (spikes <= t_max)
            spikes_in_window = spikes[spike_mask]
            
            for spike_t in spikes_in_window:
                # Find trace value at spike time
                spike_idx = np.argmin(np.abs(t_plot - spike_t))
                spike_y = trace_scaled[spike_idx]
                
                # Draw small vertical line above spike
                spike_height = TRACE_SPACING * 0.3  # Proportional to spacing
                ax_neurons.plot([spike_t, spike_t], [spike_y, spike_y + spike_height], 
                               color=COLOR_SPIKE, linewidth=1.0, solid_capstyle='butt')
        
        # Neuron label on the left
        ax_neurons.text(t_min - 0.05, offset, f'N{n_idx+1}', fontsize=FONT_SIZE_TICK, 
                       va='center', ha='right')
    
    ax_neurons.axvline(0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_neurons.axvline(STIM_DURATION_SEC, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
    ax_neurons.axvspan(0, STIM_DURATION_SEC, color=COLOR_STIM_SHADE, alpha=0.3, zorder=0)
    ax_neurons.set_xlabel('Time from stim onset (s)', fontsize=FONT_SIZE_LABEL)
    ax_neurons.set_ylabel('Individual Neurons', fontsize=FONT_SIZE_LABEL)
    ax_neurons.set_xlim(t_min, t_max)
    ax_neurons.set_ylim(-0.5, num_neurons * TRACE_SPACING + 0.5)
    ax_neurons.spines['top'].set_visible(False)
    ax_neurons.spines['right'].set_visible(False)
    ax_neurons.spines['left'].set_visible(False)
    ax_neurons.set_yticks([])  # Hide y-axis ticks
    ax_neurons.tick_params(labelsize=FONT_SIZE_TICK)
    
    plt.tight_layout()
    
    # Save figure
    fig_name = f"{MOUSE_NAME}_{session['date']}-{session['rec_id']}_Trial{trial_idx+1}"
    fig.savefig(output_path / f"{fig_name}.png", dpi=DPI, bbox_inches='tight', facecolor='white')
    fig.savefig(output_path / f"{fig_name}.pdf", bbox_inches='tight', facecolor='white')
    plt.close(fig)
    
    print(f"    Saved: {fig_name}")


# =============================================================================
# MAIN
# =============================================================================

def main():
    """Main function to generate all figures."""
    
    print("=" * 70)
    print("  CELLULAR DBS TRACE PLOTTING")
    print("=" * 70)
    
    # Create output folder
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput folder: {OUTPUT_FOLDER}")
    
    # Process each session
    for session in SESSIONS:
        session_id = f"{session['date']}-{session['rec_id']}"
        
        # Skip if not in SESSIONS_TO_PROCESS (when specified)
        if SESSIONS_TO_PROCESS is not None and session_id not in SESSIONS_TO_PROCESS:
            print(f"\n  Skipping: {session_id} (not in SESSIONS_TO_PROCESS)")
            continue
        print(f"\n{'='*70}")
        print(f"  Processing: {session['date']}-{session['rec_id']} ({session['label']})")
        print(f"{'='*70}")
        
        # Load data
        data = load_cellular_data(session)
        
        if data is None:
            print(f"  SKIPPED: Could not load data")
            continue
        
        # Generate all-trial averaged figure
        plot_all_trial_averaged(data, OUTPUT_FOLDER)
        
        # Generate single-trial figure (Trial 1)
        plot_single_trial(data, OUTPUT_FOLDER, trial_idx=REPRESENTATIVE_TRIAL - 1)
    
    print("\n" + "=" * 70)
    print("  PLOTTING COMPLETE")
    print("=" * 70)


if __name__ == "__main__":
    main()
