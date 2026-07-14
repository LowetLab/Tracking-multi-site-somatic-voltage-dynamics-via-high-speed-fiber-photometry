%% ============================================================================
%  STIMULATION TRIAL SPECTRAL ANALYSIS PIPELINE
%  ============================================================================
%  This script analyzes stimulation trials and saves spectral data for
%  Python visualization. It processes each trial independently (no pooling).
%
%  OUTPUT STRUCTURE:
%    For each trial, saves a MAT file containing:
%    - Coherence (overall, pre-stim, transient, sustained, post-stim)
%    - PSD LFP (overall, pre-stim, transient, sustained, post-stim)
%    - PSD Fiber (overall, pre-stim, transient, sustained, post-stim)
%    - Time-resolved coherence (for heatmaps)
%    - Spectrograms (LFP and Fiber)
%    - Raw traces (LFP, Fiber, Motion)
%    - Metadata (animal, session, condition, timing)
%
%  USAGE:
%    1. Modify the QUICK CONFIG section below
%    2. Run: run_stim_spectral_pipeline
%
%  NOTE: This pipeline is SEPARATE from the baseline analysis pipeline.
%        It does NOT affect baseline analysis behavior.
%
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  QUICK CONFIG - MODIFY THESE FOR YOUR RUN
%  ============================================================================

% Which methods to use?
METHODS_TO_RUN = {'mscohere', 'fieldtrip'};  % {'mscohere'}, {'fieldtrip'}, or both

% Which animals to process? (empty = all animals in stim_database)
% Options: {'Animal01'}, {'Animal02'}, {'Animal03'}, or {} for all -- see config/stim_analysis_config.m
ANIMALS_TO_PROCESS = {'Animal02'};

% Which sessions to process? (empty = all sessions for each animal)
% Format: session_id strings
%
% For Animal02-style 1s-stim DBS comparisons:
%   Amplitude Balanced: {'01_02_26-R6', '01_02_26-R9'}   (135Hz vs 40Hz AmpBalanced)
%   Energy Balanced:    {'01_02_26-R6', '01_02_26-R10'}  (135Hz vs 40Hz EnergyBalanced)
%   All DBS sessions:   {'01_02_26-R6', '01_02_26-R9', '01_02_26-R10'}
%
% For Animal03-style 10s-stim DBS comparisons:
%   Amplitude Balanced: {'01_03_26-R3', '01_03_26-R5'}   (135Hz vs 40Hz AmpBalanced @ 2.9V)
%   Energy Balanced:    {'01_03_26-R3', '01_03_26-R6'}   (135Hz vs 40Hz EnergyBalanced @ 4.5V)
%   All sessions:       {'01_03_26-R3', '01_03_26-R5', '01_03_26-R6'}
%
SESSIONS_TO_PROCESS = {'01_02_26-R6', '01_02_26-R10'};  % Example: Amplitude Balanced comparison

% Override default fiber index?
FIBER_INDEX = 1;  % 1-based index

% Figure generation (for debugging)
GENERATE_DEBUG_FIGURES = false;

%% ============================================================================
%  ADD PATHS
%  ============================================================================

% Get the directory containing this script
script_dir = fileparts(mfilename('fullpath'));

% Add subdirectories to path
addpath(fullfile(script_dir, 'config'));
addpath(fullfile(script_dir, 'core'));

% Load configuration
cfg = stim_analysis_config();

% Override with quick config
cfg.methods = METHODS_TO_RUN;
cfg.fiber_index = FIBER_INDEX;

% Add toolbox paths
for i = 1:length(cfg.paths.toolboxes)
    if exist(cfg.paths.toolboxes{i}, 'dir')
        addpath(genpath(cfg.paths.toolboxes{i}));
    end
end

% Check FieldTrip
if any(strcmpi(METHODS_TO_RUN, 'fieldtrip'))
    if ~exist('ft_defaults', 'file')
        warning('FieldTrip not found. Removing fieldtrip from methods.');
        cfg.methods(strcmpi(cfg.methods, 'fieldtrip')) = [];
    else
        ft_defaults;
        fprintf('FieldTrip available.\n');
    end
end

%% ============================================================================
%  LOAD STIMULATION DATABASE
%  ============================================================================

stim_database = cfg.stim_database;
num_animals_total = length(stim_database);

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  STIMULATION TRIAL SPECTRAL ANALYSIS PIPELINE\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  Methods: %s\n', strjoin(cfg.methods, ', '));
fprintf('  Fiber index: %d\n', cfg.fiber_index);
fprintf('  Total animals in database: %d\n', num_animals_total);
fprintf('  Output root: %s\n', cfg.paths.output_root);
fprintf('════════════════════════════════════════════════════════════════════════\n\n');

% Filter animals if specified
if ~isempty(ANIMALS_TO_PROCESS)
    animal_names = {stim_database.mouse_id};
    animal_mask = ismember(animal_names, ANIMALS_TO_PROCESS);
    stim_database = stim_database(animal_mask);
    fprintf('Processing %d selected animals: %s\n\n', ...
        length(stim_database), strjoin({stim_database.mouse_id}, ', '));
else
    fprintf('Processing all %d animals.\n\n', num_animals_total);
end

% Create output directory
if ~exist(cfg.paths.output_root, 'dir')
    mkdir(cfg.paths.output_root);
    fprintf('Created output directory: %s\n', cfg.paths.output_root);
end

%% ============================================================================
%  MAIN PROCESSING LOOP
%  ============================================================================

total_trials_processed = 0;
total_trials_failed = 0;

for animal_idx = 1:length(stim_database)
    animal = stim_database(animal_idx);
    
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  ANIMAL %d/%d: %s\n', animal_idx, length(stim_database), animal.mouse_id);
    fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');
    fprintf('  Project: %s\n', animal.project);
    fprintf('  Sessions: %d\n', length(animal.sessions));
    
    % Create animal output directory
    animal_output_dir = fullfile(cfg.paths.output_root, animal.mouse_id);
    if ~exist(animal_output_dir, 'dir')
        mkdir(animal_output_dir);
    end
    
    % =========================================================================
    %  PROCESS EACH SESSION
    % =========================================================================
    for sess_idx = 1:length(animal.sessions)
        session = animal.sessions(sess_idx);
        
        % Skip if not in sessions to process
        if ~isempty(SESSIONS_TO_PROCESS)
            if ~ismember(session.session_id, SESSIONS_TO_PROCESS)
                fprintf('  Skipping session %s (not in SESSIONS_TO_PROCESS)\n', ...
                    session.session_id);
                continue;
            end
        end
        
        fprintf('\n  ─── Session: %s (%s) ───\n', ...
            session.session_id, session.condition_name);
        fprintf('      Trials: %d\n', session.num_trials);
        fprintf('      Stim params: %d Hz, %d µs, %d V\n', ...
            session.stim_params.frequency_hz, ...
            session.stim_params.pulse_width_us, ...
            session.stim_params.voltage);
        
        % Create session output directory
        session_output_dir = fullfile(animal_output_dir, ...
            sprintf('%s_%s', session.session_id, session.condition_name));
        if ~exist(session_output_dir, 'dir')
            mkdir(session_output_dir);
        end
        
        % =====================================================================
        %  PROCESS EACH TRIAL
        % =====================================================================
        for trial_num = 1:session.num_trials
            fprintf('\n    Trial %d/%d:\n', trial_num, session.num_trials);
            
            % -----------------------------------------------------------------
            %  PROCESS WITH EACH METHOD
            % -----------------------------------------------------------------
            for method_idx = 1:length(cfg.methods)
                method = cfg.methods{method_idx};
                fprintf('      [%s] ', upper(method));
                
                % Run analysis
                try
                    results = stim_spectral_analysis(animal, session, trial_num, method, cfg);
                    
                    % Validate results before saving
                    if ~isfield(results, 'success') || ~results.success
                        fprintf('FAILED (analysis returned unsuccessful)\n');
                        total_trials_failed = total_trials_failed + 1;
                        continue;
                    end
                    
                    % Additional validation: check if critical fields exist
                    if ~isfield(results, 'coherence') || ~isfield(results, 'psd_lfp') || ~isfield(results, 'freq')
                        fprintf('FAILED (results structure incomplete)\n');
                        total_trials_failed = total_trials_failed + 1;
                        continue;
                    end
                    
                    % Save results
                    output_filename = sprintf('%s_%s_Trial%d_%s_SpectralResults.mat', ...
                        animal.mouse_id, session.session_id, trial_num, method);
                    output_path = fullfile(session_output_dir, output_filename);
                    
                    % Ensure output directory exists
                    if ~exist(session_output_dir, 'dir')
                        mkdir(session_output_dir);
                    end
                    
                    % Robustly delete existing file if it exists
                    if exist(output_path, 'file')
                        max_retries = 3;
                        deleted = false;
                        for retry = 1:max_retries
                            try
                                % Try MATLAB delete first
                                delete(output_path);
                                pause(0.2); % Wait for OS to release file handle
                                
                                % Verify deletion
                                if ~exist(output_path, 'file')
                                    deleted = true;
                                    fprintf('      Deleted existing file: %s\n', output_filename);
                                    break;
                                end
                            catch del_err
                                if retry < max_retries
                                    pause(0.3 * retry); % Exponential backoff
                                else
                                    % Last resort: try to rename it
                                    try
                                        backup_name = [output_path, '.old_', datestr(now, 'HHMMSS')];
                                        if exist(backup_name, 'file')
                                            delete(backup_name);
                                        end
                                        movefile(output_path, backup_name);
                                        fprintf('      Renamed existing file to: %s\n', [output_filename, '.old']);
                                        deleted = true;
                                    catch
                                        warning('Could not delete or rename existing file: %s (delete error: %s)', ...
                                            output_path, del_err.message);
                                    end
                                end
                            end
                        end
                    end
                    
                    % Save with compression
                    StimSpectralResults = results;
                    try
                        % Clean the struct of any problematic fields before saving
                        StimSpectralResults = clean_struct_for_saving(StimSpectralResults);
                        
                        % Diagnostic: Check struct size
                        try
                            whos_info = whos('StimSpectralResults');
                            struct_size_mb = whos_info.bytes / (1024^2);
                            if struct_size_mb > 100
                                fprintf('      Warning: Struct is very large (%.1f MB)\n', struct_size_mb);
                            end
                        catch
                            % Ignore whos errors
                        end
                        
                        % ===============================================================
                        % CRITICAL: HDF5 (-v7.3) has issues writing to network paths
                        % Solution: Save to LOCAL temp file first, then COPY to network
                        % ===============================================================
                        
                        % Create truly local temp file (not on network!)
                        local_temp_dir = tempdir;  % System temp directory (e.g., C:\Users\...\Temp)
                        temp_filename = sprintf('stim_spectral_%s.mat', datestr(now, 'yyyymmdd_HHMMSS_FFF'));
                        local_temp_path = fullfile(local_temp_dir, temp_filename);
                        
                        max_save_retries = 3;
                        save_success = false;
                        
                        for save_retry = 1:max_save_retries
                            try
                                % Delete local temp file if it exists
                                if exist(local_temp_path, 'file')
                                    delete(local_temp_path);
                                    pause(0.05);
                                end
                                
                                % Save to LOCAL temp file first (avoids network HDF5 issues)
                                % Use -v7 by default (more reliable for network paths)
                                % Fall back to -v7.3 only if struct is very large
                                use_v73 = false;
                                try
                                    whos_check = whos('StimSpectralResults');
                                    if whos_check.bytes > 2e9  % Only use v7.3 for >2GB
                                        use_v73 = true;
                                    end
                                catch
                                end
                                
                                if use_v73
                                    fprintf('      Saving to local temp (v7.3 - large file): %s\n', local_temp_path);
                                    save(local_temp_path, 'StimSpectralResults', '-v7.3');
                                else
                                    fprintf('      Saving to local temp (v7): %s\n', local_temp_path);
                                    save(local_temp_path, 'StimSpectralResults', '-v7');
                                end
                                
                                % Verify local temp file was created and has content
                                pause(0.1); % Brief pause for file system
                                temp_info = dir(local_temp_path);
                                if isempty(temp_info)
                                    error('Local temp file was not created');
                                end
                                if temp_info.bytes == 0
                                    error('Local temp file is empty (0 bytes)');
                                end
                                
                                fprintf('      Local temp created: %.1f KB\n', temp_info.bytes / 1024);
                                
                                % Now COPY to network location (copyfile is more reliable than movefile for network)
                                fprintf('      Copying to network: %s\n', output_path);
                                
                                % Ensure target directory exists
                                if ~exist(session_output_dir, 'dir')
                                    mkdir(session_output_dir);
                                    pause(0.2);
                                end
                                
                                % Delete existing target file
                                if exist(output_path, 'file')
                                    delete(output_path);
                                    pause(0.1);
                                end
                                
                                % Copy to network (more reliable than movefile for network paths)
                                [copy_status, copy_msg] = copyfile(local_temp_path, output_path, 'f');
                                if ~copy_status
                                    error('Failed to copy to network: %s', copy_msg);
                                end
                                
                                % Verify final file on network
                                pause(0.2);
                                file_info = dir(output_path);
                                if isempty(file_info)
                                    error('Network file was not created');
                                end
                                if file_info.bytes == 0
                                    error('Network file is empty (0 bytes)');
                                end
                                
                                % Clean up local temp
                                try
                                    delete(local_temp_path);
                                catch
                                    % Ignore cleanup errors
                                end
                                
                                save_success = true;
                                break;
                                
                            catch save_retry_err
                                % Clean up local temp file on error
                                if exist(local_temp_path, 'file')
                                    try
                                        delete(local_temp_path);
                                    catch
                                    end
                                end
                                
                                fprintf('      Save attempt %d failed: %s\n', save_retry, save_retry_err.message);
                                
                                if save_retry < max_save_retries
                                    pause(0.5 * save_retry);
                                    continue;
                                else
                                    rethrow(save_retry_err);
                                end
                            end
                        end
                        
                        if save_success
                            fprintf('      Saved: %s (%d bytes)\n', output_filename, file_info.bytes);
                            total_trials_processed = total_trials_processed + 1;
                        else
                            error('Failed to save after %d attempts', max_save_retries);
                        end
                    catch save_err
                        fprintf('      ERROR saving file: %s\n', save_err.message);
                        if ~isempty(save_err.stack)
                            fprintf('      Error location: %s (line %d)\n', ...
                                save_err.stack(1).name, save_err.stack(1).line);
                        end
                        total_trials_failed = total_trials_failed + 1;
                    end
                    
                catch ME
                    fprintf('ERROR: %s\n', ME.message);
                    fprintf('      Stack trace: %s\n', getReport(ME, 'extended'));
                    total_trials_failed = total_trials_failed + 1;
                end
            end
        end
        
        fprintf('\n      Session %s completed.\n', session.session_id);
    end
    
    fprintf('\n  ✓ Animal %s completed.\n', animal.mouse_id);
end

%% ============================================================================
%  SUMMARY
%  ============================================================================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  STIMULATION PIPELINE COMPLETED\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  Animals processed: %d\n', length(stim_database));
fprintf('  Methods: %s\n', strjoin(cfg.methods, ', '));
fprintf('  Total trials processed: %d\n', total_trials_processed);
fprintf('  Total trials failed: %d\n', total_trials_failed);
fprintf('  Output directory: %s\n', cfg.paths.output_root);
fprintf('════════════════════════════════════════════════════════════════════════\n');

%% ============================================================================
%  GENERATE SUMMARY FILE
%  ============================================================================

% Create a summary file for Python to easily find all outputs
summary = struct();
summary.pipeline_version = '1.0';
summary.run_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
summary.methods = cfg.methods;
summary.animals = {};

for animal_idx = 1:length(stim_database)
    animal = stim_database(animal_idx);
    animal_summary = struct();
    animal_summary.mouse_id = animal.mouse_id;
    animal_summary.sessions = {};
    
    for sess_idx = 1:length(animal.sessions)
        session = animal.sessions(sess_idx);
        
        % Skip if not processed
        if ~isempty(SESSIONS_TO_PROCESS)
            if ~ismember(session.session_id, SESSIONS_TO_PROCESS)
                continue;
            end
        end
        
        sess_summary = struct();
        sess_summary.session_id = session.session_id;
        sess_summary.condition = session.condition_name;
        sess_summary.stim_params = session.stim_params;
        sess_summary.num_trials = session.num_trials;
        sess_summary.output_dir = fullfile(cfg.paths.output_root, animal.mouse_id, ...
            sprintf('%s_%s', session.session_id, session.condition_name));
        
        animal_summary.sessions{end+1} = sess_summary;
    end
    
    summary.animals{end+1} = animal_summary;
end

summary.timing_params = cfg.stim_timing;
summary.output_root = cfg.paths.output_root;

% Save summary
summary_path = fullfile(cfg.paths.output_root, 'pipeline_summary.mat');
save(summary_path, 'summary');
fprintf('\nSummary saved to: %s\n', summary_path);

fprintf('\nDone!\n');

%% ============================================================================
%  HELPER FUNCTION: Clean struct for saving
%  ============================================================================
function cleaned = clean_struct_for_saving(s)
    %CLEAN_STRUCT_FOR_SAVING Remove problematic fields that can't be saved
    
    cleaned = s;
    
    % Recursively clean struct fields
    if isstruct(cleaned)
        fn = fieldnames(cleaned);
        for i = 1:length(fn)
            field_val = cleaned.(fn{i});
            
            % Remove function handles
            if isa(field_val, 'function_handle')
                cleaned = rmfield(cleaned, fn{i});
                continue;
            end
            
            % Recursively clean nested structs
            if isstruct(field_val)
                cleaned.(fn{i}) = clean_struct_for_saving(field_val);
            elseif iscell(field_val)
                % Clean cell arrays
                for j = 1:numel(field_val)
                    if isstruct(field_val{j})
                        field_val{j} = clean_struct_for_saving(field_val{j});
                    elseif isa(field_val{j}, 'function_handle')
                        field_val{j} = []; % Replace function handles with empty
                    end
                end
                cleaned.(fn{i}) = field_val;
            end
        end
    end
end
