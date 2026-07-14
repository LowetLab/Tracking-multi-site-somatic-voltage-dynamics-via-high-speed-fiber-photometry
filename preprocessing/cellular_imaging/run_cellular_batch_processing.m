%% ============================================================================
%  BATCH CELLULAR PREPROCESSING - WRAPPER SCRIPT
%  ============================================================================
%  Processes multiple cellular imaging sessions automatically.
%  Uses the session database for configuration.
%
%  USAGE:
%  1. Edit SESSIONS_TO_PROCESS below to specify which sessions to process
%  2. Run this script
%  3. ROI selection is done once per session (shared across trials)
%
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  CONFIGURATION
%% ============================================================================

% Add required paths
% External toolboxes via the centralised config (config/lab_paths.m).
% Override per machine via config/paths_local.m. See config/README.md.
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'config'));
setup_lab_paths();
addpath(fullfile(fileparts(mfilename('fullpath')), 'config'));  % local config (cellular_session_database.m)

% Specify sessions to process  -- EDIT THIS to your own sessions
% Format: {'mouse_name', 'recording_date', 'recording_id'}
SESSIONS_TO_PROCESS = {
    'Animal01', '01_06_25', 'R2';    % e.g. 40Hz Energy Balanced
    % 'Animal01', '01_06_25', 'R10'; % e.g. 130Hz (uncomment to process)
    % 'Animal01', '01_06_25', 'R1';  % e.g. 40Hz Amp Balanced (uncomment to process)
};

% Processing options
SKIP_IF_ALREADY_PROCESSED = true;  % Skip sessions that already have output files
INTERACTIVE_ROI_SELECTION = true;  % If false, will try to load ROIs from previous session

%% ============================================================================
%  LOAD SESSION DATABASE
%% ============================================================================

fprintf('Loading session database...\n');
all_sessions = cellular_session_database();
print_session_summary(all_sessions);

%% ============================================================================
%  PROCESS EACH SESSION
%% ============================================================================

num_sessions_to_process = size(SESSIONS_TO_PROCESS, 1);
processing_results = cell(num_sessions_to_process, 1);

fprintf('\n');
fprintf('================================================================\n');
fprintf('  STARTING BATCH PROCESSING: %d sessions\n', num_sessions_to_process);
fprintf('================================================================\n\n');

for session_idx = 1:num_sessions_to_process
    mouse_name = SESSIONS_TO_PROCESS{session_idx, 1};
    recording_date = SESSIONS_TO_PROCESS{session_idx, 2};
    recording_id = SESSIONS_TO_PROCESS{session_idx, 3};
    
    fprintf('\n');
    fprintf('****************************************************************\n');
    fprintf('  SESSION %d/%d: %s %s-%s\n', session_idx, num_sessions_to_process, ...
        mouse_name, recording_date, recording_id);
    fprintf('****************************************************************\n');
    
    % Find session in database
    session_match = find(strcmp({all_sessions.mouse_name}, mouse_name) & ...
                        strcmp({all_sessions.recording_date}, recording_date) & ...
                        strcmp({all_sessions.recording_id}, recording_id));
    
    if isempty(session_match)
        warning('Session not found in database. Using default parameters.');
        session_info = struct();
        session_info.mouse_name = mouse_name;
        session_info.recording_date = recording_date;
        session_info.recording_id = recording_id;
        session_info.experiment_type = 'Unknown';
        session_info.dbs_frequency_hz = NaN;
        session_info.dbs_comparison_type = '';
    else
        session_info = all_sessions(session_match);
    end
    
    % Check if already processed
    output_folder = fullfile(lab_paths().project_root, 'Preprocessed_Data', 'Cellular', ...
        mouse_name, sprintf('%s-%s', recording_date, recording_id));
    output_file = fullfile(output_folder, sprintf('%s_%s-%s_CellularAnalysis.mat', ...
        mouse_name, recording_date, recording_id));
    
    if SKIP_IF_ALREADY_PROCESSED && exist(output_file, 'file')
        fprintf('Session already processed. Skipping.\n');
        fprintf('Output file: %s\n', output_file);
        processing_results{session_idx} = struct('status', 'skipped', 'reason', 'already_processed');
        continue;
    end
    
    % Process the session
    try
        % Call the main processing function with session parameters
        result = process_cellular_session(session_info, INTERACTIVE_ROI_SELECTION);
        processing_results{session_idx} = struct('status', 'success', 'result', result);
        fprintf('\nSession %s-%s processed successfully.\n', recording_date, recording_id);
    catch ME
        warning('Failed to process session %s-%s: %s', recording_date, recording_id, ME.message);
        processing_results{session_idx} = struct('status', 'failed', 'error', ME);
    end
end

%% ============================================================================
%  SUMMARY
%% ============================================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('  BATCH PROCESSING COMPLETE\n');
fprintf('================================================================\n\n');

num_success = sum(cellfun(@(x) strcmp(x.status, 'success'), processing_results));
num_skipped = sum(cellfun(@(x) strcmp(x.status, 'skipped'), processing_results));
num_failed = sum(cellfun(@(x) strcmp(x.status, 'failed'), processing_results));

fprintf('  Successful: %d\n', num_success);
fprintf('  Skipped:    %d\n', num_skipped);
fprintf('  Failed:     %d\n', num_failed);
fprintf('\n');

% List any failures
if num_failed > 0
    fprintf('Failed sessions:\n');
    for i = 1:length(processing_results)
        if strcmp(processing_results{i}.status, 'failed')
            fprintf('  - %s %s-%s: %s\n', ...
                SESSIONS_TO_PROCESS{i, 1}, SESSIONS_TO_PROCESS{i, 2}, ...
                SESSIONS_TO_PROCESS{i, 3}, processing_results{i}.error.message);
        end
    end
end

fprintf('================================================================\n');


%% ============================================================================
%  HELPER FUNCTION: Process a single session
%% ============================================================================

function result = process_cellular_session(session_info, interactive_rois)
    % PROCESS_CELLULAR_SESSION - Process a single cellular imaging session
    %
    % This function wraps the main processing logic from 
    % cellular_processing_multitrial_DBS.m for use in batch processing.
    
    % Extract session parameters
    mouse_name = session_info.mouse_name;
    recording_date = session_info.recording_date;
    recording_id = session_info.recording_id;
    experiment_type = session_info.experiment_type;
    
    if isfield(session_info, 'dbs_frequency_hz') && ~isnan(session_info.dbs_frequency_hz)
        dbs_frequency = session_info.dbs_frequency_hz;
    else
        dbs_frequency = 0;
    end
    
    if isfield(session_info, 'dbs_comparison_type')
        dbs_comparison = session_info.dbs_comparison_type;
    else
        dbs_comparison = '';
    end
    
    % Set up paths
    base_path_data = fullfile(lab_paths().data_root, 'CellularVoltageImaging');
    base_path_output = fullfile(lab_paths().project_root, 'Preprocessed_Data', 'Cellular');
    
    data_folder = fullfile(base_path_data, mouse_name, experiment_type, recording_date, recording_id);
    output_folder = fullfile(base_path_output, mouse_name, sprintf('%s-%s', recording_date, recording_id));
    
    % Verify data exists
    if ~exist(data_folder, 'dir')
        % Try alternative path structure
        data_folder = fullfile(base_path_data, mouse_name, recording_date, recording_id);
        if ~exist(data_folder, 'dir')
            error('Data folder not found: %s', data_folder);
        end
    end
    
    % Create output folder
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end
    
    % Detect trial folders
    folder_contents = dir(data_folder);
    folder_contents = folder_contents([folder_contents.isdir]);
    folder_contents = folder_contents(~ismember({folder_contents.name}, {'.', '..'}));
    
    trial_folders = {};
    trial_numbers = [];
    
    for i = 1:length(folder_contents)
        folder_name = folder_contents(i).name;
        match = regexp(folder_name, '_(\d+)$', 'tokens');
        if ~isempty(match)
            trial_num = str2double(match{1}{1});
            trial_folders{end+1} = folder_name;
            trial_numbers(end+1) = trial_num;
        end
    end
    
    [trial_numbers, sort_idx] = sort(trial_numbers);
    trial_folders = trial_folders(sort_idx);
    num_trials = length(trial_folders);
    
    fprintf('Found %d trials\n', num_trials);
    
    % Return basic result (full processing would call the main script)
    result = struct();
    result.data_folder = data_folder;
    result.output_folder = output_folder;
    result.num_trials = num_trials;
    result.trial_folders = trial_folders;
    result.dbs_frequency = dbs_frequency;
    result.dbs_comparison = dbs_comparison;
    
    % Note: For full processing, you would call the main processing functions here
    % or run the main script with these parameters. The current implementation
    % is a framework that demonstrates the batch processing structure.
    
    fprintf('Session parameters configured. Ready for full processing.\n');
    fprintf('  Data: %s\n', data_folder);
    fprintf('  Output: %s\n', output_folder);
    fprintf('  Trials: %d\n', num_trials);
    if dbs_frequency > 0
        fprintf('  DBS: %d Hz (%s)\n', dbs_frequency, dbs_comparison);
    end
    
    % To run full processing, uncomment and adapt:
    % run('cellular_processing_multitrial_DBS.m');
    
end


%% Helper function from database file
function print_session_summary(sessions)
    fprintf('\n=== CELLULAR SESSION DATABASE ===\n\n');
    fprintf('%-12s %-12s %-6s %-12s %-8s %-15s %-10s\n', ...
        'Mouse', 'Date', 'RecID', 'Experiment', 'DBS Hz', 'Comparison', 'Status');
    fprintf('%s\n', repmat('-', 1, 85));
    
    for i = 1:length(sessions)
        s = sessions(i);
        if isnan(s.dbs_frequency_hz)
            freq_str = '-';
        else
            freq_str = sprintf('%d', s.dbs_frequency_hz);
        end
        fprintf('%-12s %-12s %-6s %-12s %-8s %-15s %-10s\n', ...
            s.mouse_name, s.recording_date, s.recording_id, ...
            s.experiment_type, freq_str, s.dbs_comparison_type, s.status);
    end
    fprintf('\n');
end
