function animals = animal_session_database()
%% ============================================================================
%  ANIMAL AND SESSION DATABASE  -- EDIT THIS to your own cohort
%  ============================================================================
%  Centralized database of all animals, sessions, and trial paths for the
%  spectral analysis pipeline. Replace the example animals below with your
%  own -- copy one ANIMAL block per animal, incrementing animal_idx.
%
%  STRUCTURE:
%    animals(i).mouse_id      - Animal identifier
%    animals(i).project       - Project/experiment name (your data-root subfolder)
%    animals(i).sessions      - Array of session structs with:
%                                .session_id    - Session identifier (e.g., '01_01_25-R1')
%                                .date          - Recording date
%                                .num_trials    - Number of trials
%                                .trial_paths   - Cell array of trial MAT file paths
%                                .trial_labels  - Cell array of trial labels
%                                .notes         - Any special notes
%    animals(i).session_pooled_groups - Groups of sessions to pool together
%                                       (indices into the sessions array; use
%                                       this when a recording was split across
%                                       two runs on the same day, see Animal01
%                                       below)
%  ============================================================================

% Base data path -- EDIT THIS (or leave derived from config/lab_paths.m)
BASE_PATH = fullfile(lab_paths().data_root, 'FiberVoltageImaging');

animals = struct('mouse_id', {}, 'project', {}, 'sessions', {}, 'session_pooled_groups', {});

%% ============================================================================
%  ANIMAL 1 -- example with a pooled (split) session and a stimulation session
%  ============================================================================
animal_idx = 1;
animals(animal_idx).mouse_id = 'Animal01';
animals(animal_idx).project = 'FiberVoltageImaging';
base_path = fullfile(BASE_PATH, 'Animal01', 'Fiber_Voltage_Processed');

% --- Session 1: baseline recording (6 trials) ---
s = 1;
animals(animal_idx).sessions(s).session_id = '01_01_25-R1';
animals(animal_idx).sessions(s).date = '01/01/2025';
animals(animal_idx).sessions(s).num_trials = 6;
animals(animal_idx).sessions(s).trial_paths = generate_trial_paths(base_path, '01_01_25-R1', ...
    'Animal01', 6, 'fov1_baselineRecording_60sec');
animals(animal_idx).sessions(s).trial_labels = generate_trial_labels(6);
animals(animal_idx).sessions(s).notes = '';

% --- Session 2: split recording, part A (2 trials) ---
s = 2;
animals(animal_idx).sessions(s).session_id = '02_01_25-R1';
animals(animal_idx).sessions(s).date = '02/01/2025';
animals(animal_idx).sessions(s).num_trials = 2;
animals(animal_idx).sessions(s).trial_paths = generate_trial_paths(base_path, '02_01_25-R1', ...
    'Animal01', 2, 'fov1_baselineRecording_60sec');
animals(animal_idx).sessions(s).trial_labels = generate_trial_labels(2);
animals(animal_idx).sessions(s).notes = 'Pool with 02_01_25-R2 for session-level analysis';

% --- Session 3: split recording, part B (2 trials) ---
s = 3;
animals(animal_idx).sessions(s).session_id = '02_01_25-R2';
animals(animal_idx).sessions(s).date = '02/01/2025';
animals(animal_idx).sessions(s).num_trials = 2;
animals(animal_idx).sessions(s).trial_paths = generate_trial_paths(base_path, '02_01_25-R2', ...
    'Animal01', 2, 'fov1_baselineRecording_60sec');
animals(animal_idx).sessions(s).trial_labels = generate_trial_labels(2);
animals(animal_idx).sessions(s).notes = 'Pool with 02_01_25-R1 for session-level analysis';

% --- Session 4: stimulation session (10 trials) ---
s = 4;
animals(animal_idx).sessions(s).session_id = '03_01_25-R1';
animals(animal_idx).sessions(s).date = '03/01/2025';
animals(animal_idx).sessions(s).num_trials = 10;
animals(animal_idx).sessions(s).trial_paths = generate_trial_paths(base_path, '03_01_25-R1', ...
    'Animal01', 10, 'fov1_100us_40Hz_6V_1sec_bilateralANT_10_trials');
animals(animal_idx).sessions(s).trial_labels = generate_trial_labels(10);
animals(animal_idx).sessions(s).notes = 'Stimulation session: 40Hz';

% Session pooling groups (indices into sessions array)
% Group 1: Sessions 2 and 3 should be pooled as one combined session
animals(animal_idx).session_pooled_groups = {[2, 3]};

%% ============================================================================
%  ANIMAL 2 -- example with no pooling
%  ============================================================================
animal_idx = 2;
animals(animal_idx).mouse_id = 'Animal02';
animals(animal_idx).project = 'FiberVoltageImaging';
base_path = fullfile(BASE_PATH, 'Animal02', 'Fiber_Voltage_Processed');

% --- Session 1: baseline recording (5 trials) ---
s = 1;
animals(animal_idx).sessions(s).session_id = '10_01_25-R1';
animals(animal_idx).sessions(s).date = '10/01/2025';
animals(animal_idx).sessions(s).num_trials = 5;
animals(animal_idx).sessions(s).trial_paths = generate_trial_paths(base_path, '10_01_25-R1', ...
    'Animal02', 5, 'fov1_HPmPFCdual_baselineRecording_60sec');
animals(animal_idx).sessions(s).trial_labels = generate_trial_labels(5);
animals(animal_idx).sessions(s).notes = 'Baseline recording';

animals(animal_idx).session_pooled_groups = {};  % No special pooling

end

%% ============================================================================
%  HELPER FUNCTION: Generate Trial Paths
%  ============================================================================
function paths = generate_trial_paths(base_path, session_id, mouse_id, num_trials, folder_suffix)
%GENERATE_TRIAL_PATHS Generate standard trial paths
%
%  Path format:
%    base_path/session_id/Trial{N}_{folder_suffix}_{N}/
%        {mouse_id}-{session_id}_Trial{N}_FiberPhotometry_Analysis.mat

paths = cell(1, num_trials);

for t = 1:num_trials
    trial_folder = sprintf('Trial%d_%s_%d', t, folder_suffix, t);
    trial_filename = sprintf('%s-%s_Trial%d_FiberPhotometry_Analysis.mat', ...
        mouse_id, session_id, t);
    paths{t} = fullfile(base_path, session_id, trial_folder, trial_filename);
end

end

%% ============================================================================
%  HELPER FUNCTION: Generate Trial Labels
%  ============================================================================
function labels = generate_trial_labels(num_trials)
%GENERATE_TRIAL_LABELS Generate standard trial labels

labels = cell(1, num_trials);
for t = 1:num_trials
    labels{t} = sprintf('Trial %d', t);
end

end
