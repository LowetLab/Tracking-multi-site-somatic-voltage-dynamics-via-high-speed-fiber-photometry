function sessions = cellular_session_database()
% CELLULAR_SESSION_DATABASE - Central database for all cellular imaging sessions
%
% This file contains metadata for all voltage imaging recordings that need to be
% processed. Edit this file to add new sessions or modify existing ones.
%
% Fields for each session:
%   mouse_name          - Animal identifier (e.g., 'Animal01')
%   recording_date      - Date in DD_MM_YY format
%   recording_id        - Recording ID (e.g., 'R1', 'R2', 'R9')
%   experiment_type     - Type of experiment ('DBS', 'Baseline')
%   dbs_frequency_hz    - Stimulation frequency (40, 130, etc.), NaN for non-DBS
%   dbs_comparison_type - 'EnergyBalanced', 'AmpBalanced', or '' for non-DBS
%   num_trials_expected - Expected number of trials (for validation)
%   num_neurons_expected- Expected number of neurons (for validation, can be NaN)
%   notes               - Any additional notes
%   status              - Processing status ('pending', 'processed', 'needs_review')
%
% Usage:
%   sessions = cellular_session_database();
%   animal_sessions = sessions(strcmp({sessions.mouse_name}, 'Animal01'));

    sessions = struct([]);
    idx = 0;

    %% ========================================================================
    %  EDIT THIS: add your own animal's DBS sessions below
    %% ========================================================================

    % Example: 40Hz Amplitude Balanced
    idx = idx + 1;
    sessions(idx).mouse_name = 'Animal01';
    sessions(idx).recording_date = '01_06_25';
    sessions(idx).recording_id = 'R1';
    sessions(idx).experiment_type = 'DBS';
    sessions(idx).dbs_frequency_hz = 40;
    sessions(idx).dbs_comparison_type = 'AmpBalanced';
    sessions(idx).num_trials_expected = 4;
    sessions(idx).num_neurons_expected = 11;
    sessions(idx).notes = '40Hz Amplitude Balanced - compare with R9 for Amp Balanced analysis';
    sessions(idx).status = 'pending';

    % Example: 40Hz Energy Balanced
    idx = idx + 1;
    sessions(idx).mouse_name = 'Animal01';
    sessions(idx).recording_date = '01_06_25';
    sessions(idx).recording_id = 'R2';
    sessions(idx).experiment_type = 'DBS';
    sessions(idx).dbs_frequency_hz = 40;
    sessions(idx).dbs_comparison_type = 'EnergyBalanced';
    sessions(idx).num_trials_expected = 5;
    sessions(idx).num_neurons_expected = 11;
    sessions(idx).notes = '40Hz Energy Balanced - compare with R9 for Energy Balanced analysis';
    sessions(idx).status = 'pending';

    % Example: 130Hz (comparison session for both 40Hz conditions)
    idx = idx + 1;
    sessions(idx).mouse_name = 'Animal01';
    sessions(idx).recording_date = '01_06_25';
    sessions(idx).recording_id = 'R9';
    sessions(idx).experiment_type = 'DBS';
    sessions(idx).dbs_frequency_hz = 130;
    sessions(idx).dbs_comparison_type = 'Reference';  % Used as reference for both comparisons
    sessions(idx).num_trials_expected = 5;
    sessions(idx).num_neurons_expected = 11;
    sessions(idx).notes = '130Hz Reference - paired with R1 (Amp) and R2 (Energy) for comparisons';
    sessions(idx).status = 'pending';

    %% ========================================================================
    %  TEMPLATE FOR NEW SESSIONS
    %% ========================================================================
    % Copy and modify this template to add new sessions:
    %
    % idx = idx + 1;
    % sessions(idx).mouse_name = 'MouseName';
    % sessions(idx).recording_date = 'DD_MM_YY';
    % sessions(idx).recording_id = 'RX';
    % sessions(idx).experiment_type = 'DBS';  % or 'Baseline', etc.
    % sessions(idx).dbs_frequency_hz = XX;     % NaN for non-DBS
    % sessions(idx).dbs_comparison_type = '';  % 'EnergyBalanced', 'AmpBalanced', or ''
    % sessions(idx).num_trials_expected = X;
    % sessions(idx).num_neurons_expected = NaN;  % NaN if unknown
    % sessions(idx).notes = '';
    % sessions(idx).status = 'pending';
    
end

%% Helper functions for session database queries
function selected = get_sessions_by_mouse(sessions, mouse_name)
    % GET_SESSIONS_BY_MOUSE - Filter sessions by mouse name
    selected = sessions(strcmp({sessions.mouse_name}, mouse_name));
end

function selected = get_sessions_by_experiment(sessions, experiment_type)
    % GET_SESSIONS_BY_EXPERIMENT - Filter sessions by experiment type
    selected = sessions(strcmp({sessions.experiment_type}, experiment_type));
end

function selected = get_dbs_comparison_pair(sessions, mouse_name, date, comparison_type)
    % GET_DBS_COMPARISON_PAIR - Get the two sessions needed for a DBS comparison
    %
    % Returns a struct with:
    %   .low_freq   - The 40Hz session
    %   .high_freq  - The 130Hz session
    
    mouse_sessions = get_sessions_by_mouse(sessions, mouse_name);
    date_sessions = mouse_sessions(strcmp({mouse_sessions.recording_date}, date));
    
    if strcmp(comparison_type, 'EnergyBalanced')
        % For energy balanced: R2 (40Hz) vs R9 (130Hz)
        selected.low_freq = date_sessions(strcmp({date_sessions.recording_id}, 'R2'));
        selected.high_freq = date_sessions(strcmp({date_sessions.recording_id}, 'R9'));
    elseif strcmp(comparison_type, 'AmpBalanced')
        % For amplitude balanced: R1 (40Hz) vs R9 (130Hz)
        selected.low_freq = date_sessions(strcmp({date_sessions.recording_id}, 'R1'));
        selected.high_freq = date_sessions(strcmp({date_sessions.recording_id}, 'R9'));
    else
        error('Unknown comparison type: %s', comparison_type);
    end
end

function print_session_summary(sessions)
    % PRINT_SESSION_SUMMARY - Display summary table of all sessions
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
