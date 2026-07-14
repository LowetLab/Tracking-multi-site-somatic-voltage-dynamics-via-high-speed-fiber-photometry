%% ============================================================================
%  MASTER SPECTRAL ANALYSIS PIPELINE
%  ============================================================================
%  This script orchestrates the complete spectral analysis workflow:
%    1. Single-trial analysis (heatmaps + spectra)
%    2. Session-pooled analysis (spectra only)
%    3. Animal-pooled analysis (spectra only)
%    4. Group-level statistics (optional)
%
%  USAGE:
%    1. Modify the QUICK CONFIG section below
%    2. Run: run_spectral_pipeline
%
%  All parameters are loaded from config/analysis_config.m
%  Animal/session data is loaded from config/animal_session_database.m
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  PER-RUN SETTINGS  ->  config/analysis_config.m
%  ============================================================================
%  All run settings now live in ONE place (single source of truth):
%      config/analysis_config.m
%  Edit there to change a run. The most-edited knobs are grouped at the top of
%  that file: analysis levels, coherence methods, animals, REST/RUN behaviour
%  mode + thresholds, FieldTrip taper, figure generation, and artifact mode.
%
%  For a quick ONE-OFF change without editing analysis_config.m, use the
%  OPTIONAL PER-RUN OVERRIDES block further below (after cfg is loaded).
%  ============================================================================

%% ============================================================================
%  ADD PATHS
%  ============================================================================

% Get the directory containing this script
script_dir = fileparts(mfilename('fullpath'));

% Add subdirectories to path
addpath(fullfile(script_dir, 'config'));
addpath(fullfile(script_dir, 'core'));

% Load configuration
cfg = analysis_config();

% Add toolbox paths
for i = 1:length(cfg.paths.toolboxes)
    if exist(cfg.paths.toolboxes{i}, 'dir')
        addpath(genpath(cfg.paths.toolboxes{i}));
    end
end

%% ============================================================================
%  OPTIONAL PER-RUN OVERRIDES
%  ============================================================================
%  cfg now holds every setting from analysis_config.m (the single source of
%  truth). To change one thing for THIS run only -- without editing
%  analysis_config.m -- uncomment and edit below. Examples:
%     cfg.animals_to_process = {'Animal01'};
%     cfg.behavior.classification_mode = 'standard';
%     cfg.methods = {'mscohere'};
%     cfg.artifact.mode = 'none';

% (no overrides by default)

%% ============================================================================
%  DERIVE ARTIFACT-DEPENDENT SETTINGS  (from cfg.artifact.mode)
%  ============================================================================
switch cfg.artifact.mode
    case 'none'
        cfg.output_folder_suffix = '';
    case 'exclude'
        cfg.output_folder_suffix = '_artifact_excluded';
    case 'clean'
        cfg.output_folder_suffix = '_artifact_cleaned';
    otherwise
        error('Invalid cfg.artifact.mode: %s. Use ''none'', ''exclude'', or ''clean''', cfg.artifact.mode);
end

% Legacy fields for backward compatibility
cfg.artifact_exclusion.enabled = strcmpi(cfg.artifact.mode, 'exclude');
cfg.artifact_exclusion.threshold = cfg.artifact.threshold;

% Apply output-folder suffix to the ROOT directory (not the behavior-mode level)
%   {output_root}{suffix}/{behavior_mode}/...
if ~isempty(cfg.output_folder_suffix)
    cfg.paths.output_root = [cfg.paths.output_root, cfg.output_folder_suffix];
    fprintf('  Artifact mode: %s\n', upper(cfg.artifact.mode));
    fprintf('  Output folder suffix: %s\n', cfg.output_folder_suffix);
    fprintf('  Data will be saved to: %s/%s/\n', cfg.paths.output_root, cfg.behavior.classification_mode);
else
    fprintf('  Artifact mode: NONE (using all data)\n');
end

%% ============================================================================
%  CHECK FIELDTRIP  (after cfg.methods is final)
%  ============================================================================
if any(strcmpi(cfg.methods, 'fieldtrip'))
    if ~exist('ft_defaults', 'file')
        warning('FieldTrip not found. Removing fieldtrip from methods.');
        cfg.methods(strcmpi(cfg.methods, 'fieldtrip')) = [];
    else
        ft_defaults;
        fprintf('FieldTrip available.\n');
    end
end

%% ============================================================================
%  LOAD ANIMAL DATABASE
%  ============================================================================

animals = animal_session_database();
num_animals_total = length(animals);

fprintf('\n');
fprintf('========================================================================\n');
fprintf('  SPECTRAL ANALYSIS PIPELINE\n');
fprintf('========================================================================\n');
fprintf('  Behavior Mode: %s\n', cfg.behavior.classification_mode);
fprintf('  RUN threshold: %.1f cm/s\n', cfg.behavior.run_threshold);
if strcmpi(cfg.behavior.classification_mode, 'clear')
    fprintf('  REST threshold: %.1f cm/s\n', cfg.behavior.rest_threshold);
end
fprintf('  Min bout duration: %.1f s\n', cfg.behavior.min_bout_duration_sec);
fprintf('  Methods: %s\n', strjoin(cfg.methods, ', '));
if cfg.artifact_exclusion.enabled
    fprintf('  Artifact exclusion: ENABLED (threshold: %.0f%%)\n', cfg.artifact_exclusion.threshold * 100);
else
    fprintf('  Artifact exclusion: DISABLED\n');
end
fprintf('  Total animals in database: %d\n', num_animals_total);
fprintf('  Output root: %s\n', cfg.paths.output_root);
fprintf('========================================================================\n\n');

% Filter animals if specified
if ~isempty(cfg.animals_to_process)
    animal_mask = ismember({animals.mouse_id}, cfg.animals_to_process);
    animals = animals(animal_mask);
    fprintf('Processing %d selected animals: %s\n\n', ...
        length(animals), strjoin({animals.mouse_id}, ', '));
else
    fprintf('Processing all %d animals.\n\n', num_animals_total);
end

%% ============================================================================
%  MAIN PROCESSING LOOP
%  ============================================================================

for animal_idx = 1:length(animals)
    animal = animals(animal_idx);
    
    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════════════════╗\n');
    fprintf('║  ANIMAL %d/%d: %s\n', animal_idx, length(animals), animal.mouse_id);
    fprintf('╚══════════════════════════════════════════════════════════════════════╝\n');
    fprintf('  Project: %s\n', animal.project);
    fprintf('  Sessions: %d\n', length(animal.sessions));
    
    % =========================================================================
    %  LEVEL 1: SINGLE-TRIAL ANALYSIS
    % =========================================================================
    if cfg.run_single_trial
        fprintf('\n  ─── LEVEL 1: Single-Trial Analysis ───\n');
        
        for method_idx = 1:length(cfg.methods)
            method = cfg.methods{method_idx};
            fprintf('    Method: %s\n', upper(method));
            
            for sess_idx = 1:length(animal.sessions)
                session = animal.sessions(sess_idx);
                fprintf('      Session %s (%d trials)...\n', ...
                    session.session_id, session.num_trials);
                
                % Call single-trial analysis function
                run_single_trial_analysis(animal, session, method, cfg);
            end
        end
    end
    
    % =========================================================================
    %  LEVEL 2: SESSION-POOLED ANALYSIS
    % =========================================================================
    if cfg.run_session_pooled
        fprintf('\n  ─── LEVEL 2: Session-Pooled Analysis ───\n');
        
        for method_idx = 1:length(cfg.methods)
            method = cfg.methods{method_idx};
            fprintf('    Method: %s\n', upper(method));
            
            % Get session groups (handle special pooling cases)
            session_groups = get_session_groups(animal);
            
            for group_idx = 1:length(session_groups)
                group = session_groups{group_idx};
                group_name = group.name;
                fprintf('      Session group: %s...\n', group_name);
                
                % Call session-pooled analysis function
                run_session_pooled_analysis(animal, group, method, cfg);
            end
        end
    end
    
    % =========================================================================
    %  LEVEL 3: ANIMAL-POOLED ANALYSIS (averages session spectra)
    % =========================================================================
    if cfg.run_animal_pooled
        fprintf('\n  ─── LEVEL 3: Animal-Pooled Analysis ───\n');
        
        for method_idx = 1:length(cfg.methods)
            method = cfg.methods{method_idx};
            fprintf('    Method: %s\n', upper(method));
            
            % Call animal-pooled analysis function
            run_animal_pooled_analysis(animal, method, cfg);
        end
    end
    
    % =========================================================================
    %  LEVEL 4: ANIMAL-CONCATENATED ANALYSIS (concatenates all raw data)
    % =========================================================================
    if cfg.run_animal_concatenated
        fprintf('\n  ─── LEVEL 4: Animal-Concatenated Analysis ───\n');
        fprintf('    (Concatenates all trials from all sessions, computes spectra once)\n');
        
        for method_idx = 1:length(cfg.methods)
            method = cfg.methods{method_idx};
            fprintf('    Method: %s\n', upper(method));
            
            % Call animal-concatenated analysis function
            run_animal_concatenated_analysis(animal, method, cfg);
        end
    end
    
    fprintf('\n  ✓ Animal %s completed.\n', animal.mouse_id);
end

%% ============================================================================
%  GROUP-LEVEL STATISTICS
%  ============================================================================
if cfg.run_group_statistics
    fprintf('\n');
    fprintf('════════════════════════════════════════════════════════════════════════\n');
    fprintf('  GROUP-LEVEL STATISTICS\n');
    fprintf('════════════════════════════════════════════════════════════════════════\n');
    
    % Run group statistics script
    run(fullfile(script_dir, 'figure2_coherence_group_level_stats.m'));
end

%% ============================================================================
%  SUMMARY
%  ============================================================================
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  PIPELINE COMPLETED\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
if exist('animals', 'var')
    fprintf('  Animals processed: %d\n', length(animals));
end
if exist('cfg', 'var')
    fprintf('  Methods: %s\n', strjoin(cfg.methods, ', '));
    fprintf('  Behavior mode: %s\n', cfg.behavior.classification_mode);
    if cfg.artifact_exclusion.enabled
        fprintf('  Artifact exclusion: ENABLED\n');
    end
    fprintf('  Output directory: %s\n', cfg.paths.output_root);
else
    fprintf('  (Group-level only - no config loaded)\n');
end
fprintf('════════════════════════════════════════════════════════════════════════\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================

function session_groups = get_session_groups(animal)
%GET_SESSION_GROUPS Organize sessions into groups for pooled analysis
%  Handles special cases where multiple sessions should be pooled together

session_groups = {};

% Track which sessions are in special groups
sessions_in_groups = [];
if ~isempty(animal.session_pooled_groups)
    for g = 1:length(animal.session_pooled_groups)
        group_indices = animal.session_pooled_groups{g};
        sessions_in_groups = [sessions_in_groups, group_indices]; %#ok<AGROW>
        
        % Create combined group manually to avoid struct() cell array expansion
        combined_sessions = animal.sessions(group_indices);
        group_name = sprintf('%s-combined', combined_sessions(1).session_id);
        
        group = struct();
        group.name = group_name;
        group.sessions = combined_sessions;  % Keep as struct array
        group.is_combined = true;
        session_groups{end+1} = group; %#ok<AGROW>
    end
end

% Add remaining sessions as individual groups
for s = 1:length(animal.sessions)
    if ~ismember(s, sessions_in_groups)
        % Create group struct manually to avoid struct() cell array expansion
        group = struct();
        group.name = animal.sessions(s).session_id;
        group.sessions = animal.sessions(s);  % Keep as 1x1 struct
        group.is_combined = false;
        session_groups{end+1} = group; %#ok<AGROW>
    end
end

end

function run_single_trial_analysis(animal, session, method, cfg)
%RUN_SINGLE_TRIAL_ANALYSIS Execute single-trial analysis for one session
%  Calls the unified spectral_analysis function

try
    spectral_analysis('single_trial', method, animal, session, cfg);
catch ME
    warning('Single-trial analysis failed for %s/%s: %s', ...
        animal.mouse_id, session.session_id, ME.message);
end

end

function run_session_pooled_analysis(animal, session_group, method, cfg)
%RUN_SESSION_POOLED_ANALYSIS Execute session-pooled analysis
%  Calls the unified spectral_analysis function

try
    spectral_analysis('session_pooled', method, animal, session_group, cfg);
catch ME
    warning('Session-pooled analysis failed for %s/%s: %s', ...
        animal.mouse_id, session_group.name, ME.message);
end

end

function run_animal_pooled_analysis(animal, method, cfg)
%RUN_ANIMAL_POOLED_ANALYSIS Execute animal-pooled analysis
%  Calls the unified spectral_analysis function

try
    spectral_analysis('animal_pooled', method, animal, [], cfg);
catch ME
    warning('Animal-pooled analysis failed for %s: %s', ...
        animal.mouse_id, ME.message);
end

end

function run_animal_concatenated_analysis(animal, method, cfg)
%RUN_ANIMAL_CONCATENATED_ANALYSIS Execute animal-concatenated analysis
%  Concatenates all raw data from all sessions and computes spectra once.
%  This differs from animal_pooled which averages spectra across sessions.

try
    spectral_analysis('animal_concatenated', method, animal, [], cfg);
catch ME
    warning('Animal-concatenated analysis failed for %s: %s', ...
        animal.mouse_id, ME.message);
end

end

