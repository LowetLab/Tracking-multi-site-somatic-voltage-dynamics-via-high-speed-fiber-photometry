function results = spectral_analysis(analysis_level, method, animal, session_data, cfg)
%% ============================================================================
%  SPECTRAL_ANALYSIS - Unified spectral analysis function
%  ============================================================================
%  Core analysis function that handles all analysis levels and methods.
%
%  SYNTAX:
%    results = spectral_analysis(analysis_level, method, animal, session_data, cfg)
%
%  INPUTS:
%    analysis_level - 'single_trial', 'session_pooled', or 'animal_pooled'
%    method         - 'mscohere' or 'fieldtrip'
%    animal         - Animal struct from animal_session_database()
%    session_data   - Session struct(s) to analyze
%    cfg            - Configuration from analysis_config()
%
%  OUTPUT:
%    results - Struct containing all computed spectra and metadata
%  ============================================================================

%% Validate inputs
valid_levels = {'single_trial', 'session_pooled', 'animal_pooled', 'animal_concatenated'};
valid_methods = {'mscohere', 'fieldtrip'};

if ~ismember(lower(analysis_level), valid_levels)
    error('Invalid analysis_level. Use: %s', strjoin(valid_levels, ', '));
end
if ~ismember(lower(method), valid_methods)
    error('Invalid method. Use: %s', strjoin(valid_methods, ', '));
end

%% Initialize results
results = struct();
results.analysis_level = analysis_level;
results.method = method;
results.mouse_id = animal.mouse_id;
results.analysis_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');
results.config = cfg;

%% Route to appropriate analysis function
switch lower(analysis_level)
    case 'single_trial'
        results = analyze_single_trial(results, method, animal, session_data, cfg);
        
    case 'session_pooled'
        results = analyze_session_pooled(results, method, animal, session_data, cfg);
        
    case 'animal_pooled'
        results = analyze_animal_pooled(results, method, animal, cfg);
        
    case 'animal_concatenated'
        results = analyze_animal_concatenated(results, method, animal, cfg);
end

end

%% ============================================================================
%  SINGLE TRIAL ANALYSIS
%  ============================================================================
function results = analyze_single_trial(results, method, animal, session, cfg)
%ANALYZE_SINGLE_TRIAL Process individual trials within a session

results.session_id = session.session_id;
results.num_trials = session.num_trials;
results.trials = [];

fprintf('      Processing %d trials...\n', session.num_trials);

valid_trial_count = 0;

for trial_idx = 1:session.num_trials
    trial_path = session.trial_paths{trial_idx};
    
    % Check if file exists
    if ~exist(trial_path, 'file')
        warning('Trial file not found: %s', trial_path);
        continue;
    end
    
    % CRITICAL: Extract actual trial number from filename for artifact lookup
    % Artifact masks are stored by actual trial number (e.g., Trial1, Trial5), not loop index
    actual_trial_num = extract_trial_number_from_path(trial_path);
    if isempty(actual_trial_num)
        actual_trial_num = trial_idx;  % Fallback to loop index
    end
    
    % Load trial data
    fprintf('        Trial %d (mask idx %d): Loading...', trial_idx, actual_trial_num);
    trial_data = load_trial_data(trial_path, cfg.fiber_index);
    
    if isempty(trial_data)
        fprintf(' FAILED\n');
        continue;
    end
    
    % Store original length for verification
    original_trial_samples = length(trial_data.lfp);
    
    % Apply artifact cleaning if in 'clean' mode
    % ARTIFACT MODES:
    %   'none'    - Use all data as-is
    %   'exclude' - Skip entire trials (handled elsewhere)
    %   'clean'   - Remove artifact segments from within trials
    if isfield(cfg, 'artifact') && isfield(cfg.artifact, 'mode') && ...
            strcmpi(cfg.artifact.mode, 'clean')
        [trial_data, cleaning_info] = apply_artifact_cleaning(trial_data, trial_path, actual_trial_num, cfg);
        if cleaning_info.applied && cleaning_info.removed_pct > 0
            fprintf(' [cleaned: %.1f%% removed, %d->%d samples]', ...
                cleaning_info.removed_pct, original_trial_samples, length(trial_data.lfp));
            % Verify the data length actually changed
            if length(trial_data.lfp) == original_trial_samples
                warning('single_trial:cleaningNotApplied', ...
                    'Trial %d: Cleaning reported %.1f%% removed but data length unchanged!', ...
                    trial_idx, cleaning_info.removed_pct);
            end
        elseif ~cleaning_info.applied
            fprintf(' [NO MASK for idx %d]', actual_trial_num);
        end
    end
    
    % Classify behavior
    [is_rest, is_run, behavior_info] = classify_behavior(trial_data.speed, ...
        trial_data.fs, cfg);
    
    % Compute all spectral analyses (spectrograms, coherence, PSD)
    switch lower(method)
        case 'mscohere'
            trial_results = compute_single_trial_mscohere(trial_data, is_rest, is_run, cfg, ...
                session.trial_labels{trial_idx}, trial_idx, behavior_info);
        case 'fieldtrip'
            trial_results = compute_single_trial_fieldtrip(trial_data, is_rest, is_run, cfg, ...
                session.trial_labels{trial_idx}, trial_idx, behavior_info);
    end
    
    if isempty(trial_results)
        fprintf(' FAILED\n');
        continue;
    end
    
    % Store results - handle struct array assignment properly
    valid_trial_count = valid_trial_count + 1;
    if valid_trial_count == 1
        results.trials = trial_results;  % First assignment initializes struct array
    else
        results.trials(valid_trial_count) = trial_results;  % Subsequent assignments
    end
    
    fprintf(' Done (REST: %.1f%%, RUN: %.1f%%)\n', ...
        behavior_info.pct_rest, behavior_info.pct_run);
    
    % Save per-trial file (legacy format for Python compatibility)
    save_single_trial_file(trial_results, animal, session.session_id, method, trial_idx, cfg);
end

% Check if we have any valid trials
if valid_trial_count == 0
    warning('No valid trials processed for %s/%s/%s', ...
        animal.mouse_id, session.session_id, method);
    return;
end

results.num_valid_trials = valid_trial_count;

% Note: Per-trial files already saved via save_single_trial_file()
% No combined session file needed - avoids redundancy

end

%% ============================================================================
%  SESSION POOLED ANALYSIS
%  ============================================================================
function results = analyze_session_pooled(results, method, animal, session_group, cfg)
%ANALYZE_SESSION_POOLED Pool trials within a session and compute spectra

results.session_id = session_group.name;
results.is_combined = session_group.is_combined;

% Collect all trial paths from session group
all_trial_paths = {};
all_trial_labels = {};

% Handle both struct array and single struct cases
sessions = session_group.sessions;
if isstruct(sessions) && ~isempty(sessions)
    num_sessions = numel(sessions);
    for s = 1:num_sessions
        if num_sessions == 1
            sess = sessions;  % Single struct (not array)
        else
            sess = sessions(s);  % Struct array element
        end
        
        % Safely get trial paths
        if isfield(sess, 'trial_paths')
            tp = sess.trial_paths;
            if iscell(tp)
                % Convert to row cell array and add each path
                for p = 1:length(tp)
                    all_trial_paths{end+1} = tp{p}; %#ok<AGROW>
                end
            else
                warning('trial_paths is not a cell array in session %d', s);
                continue;
            end
        else
            warning('trial_paths field missing in session %d', s);
            continue;
        end
        
        % Safely get trial labels  
        if isfield(sess, 'trial_labels') && iscell(sess.trial_labels)
            tl = sess.trial_labels;
            for p = 1:length(tl)
                all_trial_labels{end+1} = tl{p}; %#ok<AGROW>
            end
        else
            % Create empty labels for paths we just added
            n_added = length(sess.trial_paths);
            for p = 1:n_added
                all_trial_labels{end+1} = ''; %#ok<AGROW>
            end
        end
    end
else
    warning('session_group.sessions is not a valid struct');
    return;
end

results.num_trials = length(all_trial_paths);
fprintf('      Pooling %d trials from session group "%s"...\n', ...
    results.num_trials, session_group.name);

% Debug: print first trial path to verify correctness
if ~isempty(all_trial_paths)
    fprintf('        First trial: %s\n', all_trial_paths{1});
end

% Load and concatenate trial data
all_lfp = [];
all_gevi = [];
all_speed = [];
all_time = [];
fs = [];
excluded_trials = [];
included_trials = [];
trial_boundaries = [0];  % Track trial boundaries for FieldTrip edge artifact avoidance

for trial_idx = 1:length(all_trial_paths)
    trial_path = all_trial_paths{trial_idx};
    fprintf('        Loading trial %d/%d... ', trial_idx, length(all_trial_paths));
    
    if ~exist(trial_path, 'file')
        fprintf('NOT FOUND\n');
        continue;
    end
    
    % CRITICAL: Extract actual trial number from filename for artifact lookup
    % The loop trial_idx counts continuously across combined sessions (1,2,3,4...)
    % but artifact info uses per-session indices (1,2 for each session).
    % Extract actual trial number from filename pattern like "_Trial1_" or "_Trial2_"
    actual_trial_num = extract_trial_number_from_path(trial_path);
    if isempty(actual_trial_num)
        actual_trial_num = trial_idx;  % Fallback to loop index if parsing fails
        warning('Could not extract trial number from path: %s', trial_path);
    end
    
    % Handle artifact based on mode:
    % 'exclude' mode: Skip entire trials with high contamination
    % 'clean' mode:   Remove artifact segments, keep clean portions
    artifact_mode = 'none';
    if isfield(cfg, 'artifact') && isfield(cfg.artifact, 'mode')
        artifact_mode = cfg.artifact.mode;
    elseif isfield(cfg, 'artifact_exclusion') && cfg.artifact_exclusion.enabled
        artifact_mode = 'exclude';  % Legacy compatibility
    end
    
    if strcmpi(artifact_mode, 'exclude')
        [should_exclude, artifact_pct] = check_artifact_exclusion(trial_path, actual_trial_num, cfg);
        if should_exclude
            fprintf('EXCLUDED (%.1f%% artifacts)\n', artifact_pct);
            excluded_trials = [excluded_trials, trial_idx]; %#ok<AGROW>
            continue;
        end
    end
    
    trial_data = load_trial_data(trial_path, cfg.fiber_index);
    
    if isempty(trial_data)
        fprintf('LOAD FAILED\n');
        continue;
    end
    
    % Apply artifact cleaning if in 'clean' mode
    if strcmpi(artifact_mode, 'clean')
        [trial_data, cleaning_info] = apply_artifact_cleaning(trial_data, trial_path, actual_trial_num, cfg);
        if cleaning_info.applied && cleaning_info.removed_pct > 0
            fprintf('[cleaned: %.1f%% removed (orig: %.1f%%)] ', ...
                cleaning_info.removed_pct, cleaning_info.removed_pct_original);
        elseif cleaning_info.applied && cleaning_info.removed_pct == 0
            fprintf('[cleaned: 0%% removed - mask had no artifacts] ');
        elseif ~cleaning_info.applied
            fprintf('[no artifact info found] ');
        end
    end
    
    included_trials = [included_trials, trial_idx]; %#ok<AGROW>
    fprintf('OK (%.1f s)\n', trial_data.duration);
    
    if isempty(fs)
        fs = trial_data.fs;
    end
    
    % Concatenate
    all_lfp = [all_lfp; trial_data.lfp(:)]; %#ok<AGROW>
    all_gevi = [all_gevi; trial_data.gevi(:)]; %#ok<AGROW>
    all_speed = [all_speed; trial_data.speed(:)]; %#ok<AGROW>
    
    % Track trial boundary (end of this trial in concatenated indexing)
    % This is used by FieldTrip to avoid pseudo-trials crossing trial boundaries
    trial_boundaries = [trial_boundaries; length(all_lfp)]; %#ok<AGROW>
    
    % Calculate time offset safely
    if isempty(all_time)
        time_offset = 0;
    else
        time_offset = all_time(end) + 1/fs;
    end
    all_time = [all_time; trial_data.time(:) + time_offset]; %#ok<AGROW>
end

if isempty(all_lfp)
    warning('No valid trial data loaded for session %s', session_group.name);
    return;
end

% Log exclusion summary
if ~isempty(excluded_trials)
    fprintf('        ARTIFACT EXCLUSION: %d trials excluded (trials: %s)\n', ...
        length(excluded_trials), strjoin(arrayfun(@num2str, excluded_trials, 'UniformOutput', false), ', '));
end
fprintf('        Total pooled: %d trials, %d samples (%.1f s)\n', ...
    length(included_trials), length(all_lfp), length(all_lfp)/fs);

% Store exclusion info in results
results.excluded_trials = excluded_trials;
results.included_trials = included_trials;

% Create pooled data struct
pooled_data = struct();
pooled_data.lfp = all_lfp;
pooled_data.gevi = all_gevi;
pooled_data.speed = all_speed;
pooled_data.time = all_time;
pooled_data.fs = fs;
pooled_data.trial_boundaries = trial_boundaries;  % For FieldTrip edge artifact avoidance

% Classify behavior on pooled data
fprintf('        Classifying behavior...\n');
[is_rest, is_run, behavior_info] = classify_behavior(pooled_data.speed, fs, cfg);

fprintf('        Behavior: REST %.1f%%, RUN %.1f%%', ...
    behavior_info.pct_rest, behavior_info.pct_run);
if strcmpi(cfg.behavior.classification_mode, 'clear')
    fprintf(', EXCLUDED %.1f%%', behavior_info.pct_excluded);
end
fprintf('\n');

% Compute spectra
switch lower(method)
    case 'mscohere'
        spectra = compute_spectra_mscohere(pooled_data, is_rest, is_run, cfg);
    case 'fieldtrip'
        % Use boundary-aware FieldTrip computation to avoid edge artifacts
        spectra = compute_spectra_fieldtrip_with_boundaries(pooled_data, is_rest, is_run, cfg);
end

results.behavior = behavior_info;
results.spectra = spectra;
results.total_duration_sec = length(all_lfp) / fs;

% Save results
save_results(results, animal, session_group.name, 'session_pooled', method, cfg);

end

%% ============================================================================
%  ANIMAL POOLED ANALYSIS
%  ============================================================================
function results = analyze_animal_pooled(results, method, animal, cfg)
%ANALYZE_ANIMAL_POOLED Pool all sessions for an animal by AVERAGING SPECTRA
%  IMPORTANT: We compute spectra per session, then AVERAGE (not concatenate raw signals)

results.num_sessions = length(animal.sessions);
fprintf('      Averaging spectra across %d sessions...\n', results.num_sessions);

% Storage for per-session spectra
session_spectra = struct([]);
valid_sessions = 0;

for sess_idx = 1:length(animal.sessions)
    session = animal.sessions(sess_idx);
    fprintf('        Session %s: ', session.session_id);
    
    % Load and concatenate trials within this session
    sess_lfp = [];
    sess_gevi = [];
    sess_speed = [];
    fs = [];
    
    % Determine artifact handling mode
    artifact_mode = 'none';
    if isfield(cfg, 'artifact') && isfield(cfg.artifact, 'mode')
        artifact_mode = cfg.artifact.mode;
    elseif isfield(cfg, 'artifact_exclusion') && cfg.artifact_exclusion.enabled
        artifact_mode = 'exclude';  % Legacy compatibility
    end
    
    excluded_in_session = 0;
    cleaned_samples_session = 0;
    original_samples_session = 0;
    
    for trial_idx = 1:session.num_trials
        trial_path = session.trial_paths{trial_idx};
        
        if ~exist(trial_path, 'file')
            continue;
        end
        
        % CRITICAL: Extract actual trial number from filename for artifact lookup
        % Artifact masks are stored by the actual trial number in the filename (e.g., Trial1, Trial5)
        % NOT by the loop index. This matches the single-trial analysis behavior.
        actual_trial_num = extract_trial_number_from_path(trial_path);
        if isempty(actual_trial_num)
            actual_trial_num = trial_idx;  % Fallback to loop index if parsing fails
            warning('session_pooled:trialNumParseFailed', ...
                'Could not extract trial number from path: %s. Using loop index %d.', ...
                trial_path, trial_idx);
        end
        
        % Diagnostic: warn if actual trial number doesn't match loop index
        if actual_trial_num ~= trial_idx
            fprintf('\n        [NOTE: Filename has Trial%d but loop index is %d] ', ...
                actual_trial_num, trial_idx);
        end
        
        % 'exclude' mode: Skip entire trials with high contamination
        if strcmpi(artifact_mode, 'exclude')
            [should_exclude, artifact_pct] = check_artifact_exclusion(trial_path, actual_trial_num, cfg);
            if should_exclude
                excluded_in_session = excluded_in_session + 1;
                continue;  % Skip this trial
            end
        end
        
        trial_data = load_trial_data(trial_path, cfg.fiber_index);
        
        if isempty(trial_data)
            continue;
        end
        
        % Store original length for verification
        original_trial_samples = length(trial_data.lfp);
        
        % 'clean' mode: Remove artifact segments from within trials
        if strcmpi(artifact_mode, 'clean')
            original_samples_session = original_samples_session + length(trial_data.time);
            [trial_data, cleaning_info] = apply_artifact_cleaning(trial_data, trial_path, actual_trial_num, cfg);
            cleaned_samples_session = cleaned_samples_session + cleaning_info.removed_samples;
            
            % DIAGNOSTIC: Verify cleaning was applied
            if cleaning_info.applied && cleaning_info.removed_pct > 0
                fprintf('        [Trial %d (mask idx %d): %.1f%% removed, %d->%d samples] ', ...
                    trial_idx, actual_trial_num, cleaning_info.removed_pct, ...
                    cleaning_info.original_samples, cleaning_info.cleaned_samples);
                
                % CRITICAL: Verify the data length actually changed
                cleaned_trial_samples = length(trial_data.lfp);
                if cleaned_trial_samples == original_trial_samples
                    warning('session_pooled:cleaningNotApplied', ...
                        'Trial %d: Cleaning reported %.1f%% removed but data length unchanged (%d samples)!', ...
                        trial_idx, cleaning_info.removed_pct, original_trial_samples);
                end
            elseif ~cleaning_info.applied
                fprintf('        [Trial %d: NO MASK FOUND for actual_trial_num=%d] ', ...
                    trial_idx, actual_trial_num);
            end
        end
        
        if isempty(fs)
            fs = trial_data.fs;
        end
        
        % CRITICAL: Verify cleaned data lengths match before concatenation
        if length(trial_data.lfp) ~= length(trial_data.gevi) || ...
           length(trial_data.lfp) ~= length(trial_data.speed)
            warning('Length mismatch in cleaned trial %d: lfp=%d, gevi=%d, speed=%d', ...
                trial_idx, length(trial_data.lfp), length(trial_data.gevi), length(trial_data.speed));
        end
        
        sess_lfp = [sess_lfp; trial_data.lfp(:)]; %#ok<AGROW>
        sess_gevi = [sess_gevi; trial_data.gevi(:)]; %#ok<AGROW>
        sess_speed = [sess_speed; trial_data.speed(:)]; %#ok<AGROW>
    end
    
    % Report artifact handling results
    if excluded_in_session > 0
        fprintf('[%d trials excluded] ', excluded_in_session);
    end
    if strcmpi(artifact_mode, 'clean') && original_samples_session > 0
        pct_cleaned = 100 * cleaned_samples_session / original_samples_session;
        if pct_cleaned > 0
            fprintf('[%.1f%% cleaned] ', pct_cleaned);
        end
    end
    
    if isempty(sess_lfp)
        fprintf('No data\n');
        continue;
    end
    
    % Create session data struct
    sess_data = struct('lfp', sess_lfp, 'gevi', sess_gevi, ...
        'speed', sess_speed, 'fs', fs);
    
    % DIAGNOSTIC: Verify cleaned data lengths match
    if strcmpi(artifact_mode, 'clean') && original_samples_session > 0
        total_cleaned_samples = length(sess_lfp);
        expected_cleaned_samples = original_samples_session - cleaned_samples_session;
        if abs(total_cleaned_samples - expected_cleaned_samples) > 10
            warning('CLEANED DATA MISMATCH: Expected %d cleaned samples, got %d (diff: %d)', ...
                expected_cleaned_samples, total_cleaned_samples, ...
                abs(total_cleaned_samples - expected_cleaned_samples));
        end
        fprintf('        [Session total: %d original -> %d cleaned (%.1f%% removed)]\n', ...
            original_samples_session, total_cleaned_samples, ...
            100 * cleaned_samples_session / original_samples_session);
        
        % CRITICAL: Store cleaning metadata for verification in output
        cleaning_metadata.original_samples = original_samples_session;
        cleaning_metadata.cleaned_samples = total_cleaned_samples;
        cleaning_metadata.removed_samples = cleaned_samples_session;
        cleaning_metadata.removed_pct = 100 * cleaned_samples_session / original_samples_session;
        cleaning_metadata.data_checksum = sum(sess_lfp(1:min(1000, length(sess_lfp))));  % Quick checksum
    else
        cleaning_metadata.original_samples = length(sess_lfp);
        cleaning_metadata.cleaned_samples = length(sess_lfp);
        cleaning_metadata.removed_samples = 0;
        cleaning_metadata.removed_pct = 0;
        cleaning_metadata.data_checksum = sum(sess_lfp(1:min(1000, length(sess_lfp))));
    end
    
    % Classify behavior
    [is_rest, is_run, behavior_info] = classify_behavior(sess_speed, fs, cfg);
    
    % Compute spectra for this session
    switch lower(method)
        case 'mscohere'
            spectra = compute_spectra_mscohere(sess_data, is_rest, is_run, cfg);
        case 'fieldtrip'
            spectra = compute_spectra_fieldtrip(sess_data, is_rest, is_run, cfg);
    end
    
    if ~isempty(spectra) && isfield(spectra, 'coherence') && ~isempty(spectra.coherence.overall)
        valid_sessions = valid_sessions + 1;
        session_spectra(valid_sessions).session_id = session.session_id;
        session_spectra(valid_sessions).spectra = spectra;
        session_spectra(valid_sessions).behavior = behavior_info;
        session_spectra(valid_sessions).cleaning = cleaning_metadata;  % Store cleaning info
        fprintf('REST %.1f%%, RUN %.1f%% [data: %d samples, checksum: %.2f]\n', ...
            behavior_info.pct_rest, behavior_info.pct_run, ...
            cleaning_metadata.cleaned_samples, cleaning_metadata.data_checksum);
    else
        fprintf('Spectra computation failed\n');
    end
end

if valid_sessions == 0
    warning('No valid sessions for animal %s', animal.mouse_id);
    return;
end

% Average spectra across sessions
results.num_valid_sessions = valid_sessions;
results.session_spectra = session_spectra;
results.spectra = average_spectra(session_spectra);

% Average behavior stats (use arrayfun for proper struct array handling)
% CRITICAL: Store as results.behavior.pct_rest/run to match save_results expectation
avg_pct_rest = mean(arrayfun(@(x) x.behavior.pct_rest, session_spectra));
avg_pct_run = mean(arrayfun(@(x) x.behavior.pct_run, session_spectra));
results.behavior.pct_rest = avg_pct_rest;
results.behavior.pct_run = avg_pct_run;

fprintf('        Average: REST %.1f%%, RUN %.1f%% (across %d sessions)\n', ...
    avg_pct_rest, avg_pct_run, valid_sessions);

% Save results
save_results(results, animal, 'all_sessions', 'animal_pooled', method, cfg);

end

%% ============================================================================
%  ANIMAL CONCATENATED ANALYSIS
%  ============================================================================
function results = analyze_animal_concatenated(results, method, animal, cfg)
%ANALYZE_ANIMAL_CONCATENATED Concatenate all raw data across sessions, compute spectra once
%  This differs from animal_pooled which averages spectra computed per-session.
%
%  IMPORTANT: This approach concatenates raw LFP, GEVI, and speed across ALL
%  sessions for an animal. To avoid edge artifacts with FieldTrip:
%    - Trial boundaries AND session boundaries are tracked
%    - Pseudo-trials crossing any boundary are skipped
%    - REST/RUN classification uses pseudo-trial selection (not raw indexing)
%
%  CAUTION: Different sessions may have different DC offsets and amplitudes.
%  This approach assumes such differences are acceptable for the analysis.

results.num_sessions = length(animal.sessions);
fprintf('      Concatenating raw data from %d sessions...\n', results.num_sessions);

% Storage for concatenated data
all_lfp = [];
all_gevi = [];
all_speed = [];
all_time = [];
fs = [];

% Track ALL boundaries: both trial boundaries within sessions and session boundaries
% trial_boundaries: [0, trial1_end, trial2_end, ..., trialN_end]
% session_boundaries: subset of trial_boundaries marking session ends
trial_boundaries = [0];
session_boundaries = [0];

% Track metadata for each session
session_info = struct([]);
total_trials = 0;
total_excluded_trials = 0;

for sess_idx = 1:length(animal.sessions)
    session = animal.sessions(sess_idx);
    fprintf('        Session %d/%d (%s): ', sess_idx, length(animal.sessions), session.session_id);
    
    session_start_sample = length(all_lfp) + 1;  % Track where this session starts
    trials_in_session = 0;
    excluded_in_session = 0;
    cleaned_samples_session = 0;
    original_samples_session = 0;
    
    % Determine artifact handling mode
    artifact_mode = 'none';
    if isfield(cfg, 'artifact') && isfield(cfg.artifact, 'mode')
        artifact_mode = cfg.artifact.mode;
    elseif isfield(cfg, 'artifact_exclusion') && cfg.artifact_exclusion.enabled
        artifact_mode = 'exclude';
    end
    
    for trial_idx = 1:session.num_trials
        trial_path = session.trial_paths{trial_idx};
        
        if ~exist(trial_path, 'file')
            continue;
        end
        
        % Extract actual trial number from filename for artifact lookup
        actual_trial_num = extract_trial_number_from_path(trial_path);
        if isempty(actual_trial_num)
            actual_trial_num = trial_idx;
        end
        
        % 'exclude' mode: Skip entire trials with high contamination
        if strcmpi(artifact_mode, 'exclude')
            [should_exclude, ~] = check_artifact_exclusion(trial_path, actual_trial_num, cfg);
            if should_exclude
                excluded_in_session = excluded_in_session + 1;
                continue;
            end
        end
        
        trial_data = load_trial_data(trial_path, cfg.fiber_index);
        
        if isempty(trial_data)
            continue;
        end
        
        % 'clean' mode: Remove artifact segments
        if strcmpi(artifact_mode, 'clean')
            original_samples_session = original_samples_session + length(trial_data.time);
            [trial_data, cleaning_info] = apply_artifact_cleaning(trial_data, trial_path, actual_trial_num, cfg);
            cleaned_samples_session = cleaned_samples_session + cleaning_info.removed_samples;
        end
        
        if isempty(fs)
            fs = trial_data.fs;
        end
        
        % Concatenate this trial's data
        all_lfp = [all_lfp; trial_data.lfp(:)]; %#ok<AGROW>
        all_gevi = [all_gevi; trial_data.gevi(:)]; %#ok<AGROW>
        all_speed = [all_speed; trial_data.speed(:)]; %#ok<AGROW>
        
        % Track trial boundary
        trial_boundaries = [trial_boundaries; length(all_lfp)]; %#ok<AGROW>
        
        % Calculate time offset
        if isempty(all_time)
            time_offset = 0;
        else
            time_offset = all_time(end) + 1/fs;
        end
        all_time = [all_time; trial_data.time(:) + time_offset]; %#ok<AGROW>
        
        trials_in_session = trials_in_session + 1;
        total_trials = total_trials + 1;
    end
    
    % Mark session boundary (the last trial boundary is also a session boundary)
    if ~isempty(all_lfp)
        session_boundaries = [session_boundaries; length(all_lfp)]; %#ok<AGROW>
    end
    
    total_excluded_trials = total_excluded_trials + excluded_in_session;
    
    % Store session metadata
    session_info(sess_idx).session_id = session.session_id;
    session_info(sess_idx).trials_included = trials_in_session;
    session_info(sess_idx).trials_excluded = excluded_in_session;
    session_info(sess_idx).start_sample = session_start_sample;
    session_info(sess_idx).end_sample = length(all_lfp);
    session_info(sess_idx).duration_sec = (session_info(sess_idx).end_sample - session_start_sample + 1) / fs;
    
    if strcmpi(artifact_mode, 'clean') && original_samples_session > 0
        session_info(sess_idx).cleaned_pct = 100 * cleaned_samples_session / original_samples_session;
        fprintf('%d trials, %.1f%% cleaned\n', trials_in_session, session_info(sess_idx).cleaned_pct);
    else
        session_info(sess_idx).cleaned_pct = 0;
        if excluded_in_session > 0
            fprintf('%d trials (%d excluded)\n', trials_in_session, excluded_in_session);
        else
            fprintf('%d trials\n', trials_in_session);
        end
    end
end

if isempty(all_lfp)
    warning('No valid data for animal %s', animal.mouse_id);
    return;
end

% Summary statistics
total_duration_sec = length(all_lfp) / fs;
fprintf('      Total concatenated: %d trials, %d samples (%.1f min)\n', ...
    total_trials, length(all_lfp), total_duration_sec/60);
fprintf('      Trial boundaries: %d, Session boundaries: %d\n', ...
    length(trial_boundaries)-1, length(session_boundaries)-1);
if total_excluded_trials > 0
    fprintf('      Total trials excluded: %d\n', total_excluded_trials);
end

% Create concatenated data struct with boundary information
concat_data = struct();
concat_data.lfp = all_lfp;
concat_data.gevi = all_gevi;
concat_data.speed = all_speed;
concat_data.time = all_time;
concat_data.fs = fs;
concat_data.trial_boundaries = trial_boundaries;
concat_data.session_boundaries = session_boundaries;

% Classify behavior on concatenated data
fprintf('      Classifying behavior...\n');
[is_rest, is_run, behavior_info] = classify_behavior(concat_data.speed, fs, cfg);

fprintf('      Behavior: REST %.1f%%, RUN %.1f%%', ...
    behavior_info.pct_rest, behavior_info.pct_run);
if strcmpi(cfg.behavior.classification_mode, 'clear')
    fprintf(', EXCLUDED %.1f%%', behavior_info.pct_excluded);
end
fprintf('\n');

% Compute spectra
fprintf('      Computing spectra (%s)...\n', upper(method));
switch lower(method)
    case 'mscohere'
        spectra = compute_spectra_mscohere(concat_data, is_rest, is_run, cfg);
    case 'fieldtrip'
        % Use boundary-aware FieldTrip to avoid edge artifacts
        spectra = compute_spectra_fieldtrip_with_boundaries(concat_data, is_rest, is_run, cfg);
end

% Store results
results.num_trials = total_trials;
results.num_excluded_trials = total_excluded_trials;
results.session_info = session_info;
results.trial_boundaries = trial_boundaries;
results.session_boundaries = session_boundaries;
results.behavior = behavior_info;
results.spectra = spectra;
results.total_duration_sec = total_duration_sec;

% Save results
save_results(results, animal, 'all_sessions', 'animal_concatenated', method, cfg);

fprintf('      ✓ Animal concatenated analysis complete.\n');

end

%% ============================================================================
%  HELPER: Load Trial Data
%  ============================================================================
function data = load_trial_data(trial_path, fiber_index)
%LOAD_TRIAL_DATA Cached wrapper around load_trial_data_uncached.
%   spectral_analysis() is called once per (analysis_level, method) pair
%   (up to 4 levels x 2 methods = 8 calls), and each call reloads the same
%   trial .mat files from disk. Caching here (keyed by trial_path and
%   fiber_index) avoids that redundant I/O within a MATLAB session.
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
cache_key = sprintf('%s|%d', trial_path, fiber_index);
if isKey(cache, cache_key)
    data = cache(cache_key);
    return;
end
data = load_trial_data_uncached(trial_path, fiber_index);
cache(cache_key) = data;
end

function data = load_trial_data_uncached(trial_path, fiber_index)
%LOAD_TRIAL_DATA_UNCACHED Load and extract data from trial MAT file
%  EXACT match to legacy figure2_coherence_mscohere.m data loading

data = [];

% Motion conversion constants - MUST MATCH LEGACY (19.0 cm, NOT 20.0!)
% [LEGACY: figure2_coherence_mscohere.m lines 162-174]
WHEEL_DIAMETER_CM = 19.0;  % CRITICAL: Legacy uses 19.0, not 20.0!
WHEEL_CIRCUMFERENCE_CM = pi * WHEEL_DIAMETER_CM;
ENCODER_COUNTS_PER_REV = 1024;    % [LEGACY: 1024] Yumo E6B2 encoder
EPHYS_SAMPLING_RATE = 30000;       % [LEGACY: 30000] Open Ephys
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV;
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM;
MOTION_SMOOTH_SAMPLES = 10;        % [LEGACY: 10]

try
    loaded = load(trial_path);
    
    % LEGACY: Expect FiberPhotometryAnalysis struct
    if ~isfield(loaded, 'FiberPhotometryAnalysis')
        warning('FiberPhotometryAnalysis structure not found in %s', trial_path);
        return;
    end
    
    FPA = loaded.FiberPhotometryAnalysis;
    
    % Extract time vector and sampling rate - LEGACY EXACT
    if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
        t = FPA.time.time_vector_seconds;
        if isstruct(t)
            warning('time_vector_seconds is a struct, not a vector in %s', trial_path);
            return;
        end
        t = t(:);
    else
        warning('Time vector not found in %s', trial_path);
        return;
    end
    
    % Sampling rate - LEGACY EXACT
    if isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
        fs = FPA.parameters.sampling_rate;
    elseif isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
        fs = FPA.time.sampling_rate;
    else
        % Infer from time vector
        fs = 1 / median(diff(t));
    end
    
    % Make sure fs is scalar
    if ~isscalar(fs) || isstruct(fs)
        warning('Sampling rate is not a scalar in %s', trial_path);
        return;
    end
    
    % Extract ΔF/F fiber trace - LEGACY EXACT
    if isfield(FPA, 'signals') && isfield(FPA.signals, 'final_processed_traces')
        fiber_all = FPA.signals.final_processed_traces;
        if isstruct(fiber_all)
            warning('final_processed_traces is a struct in %s', trial_path);
            return;
        end
        if size(fiber_all, 2) >= fiber_index
            fiber_trace = fiber_all(:, fiber_index);
        elseif size(fiber_all, 1) >= fiber_index
            fiber_trace = fiber_all(fiber_index, :)';
        else
            warning('Fiber index %d out of range in %s', fiber_index, trial_path);
            return;
        end
        fiber_trace = fiber_trace(:);
    elseif isfield(FPA, 'signals') && isfield(FPA.signals, 'deltaF_F_traces')
        fiber_all = FPA.signals.deltaF_F_traces;
        if isstruct(fiber_all)
            warning('deltaF_F_traces is a struct in %s', trial_path);
            return;
        end
        fiber_trace = fiber_all(:, min(fiber_index, size(fiber_all, 2)));
        fiber_trace = fiber_trace(:);
    else
        warning('Fiber trace not found in %s', trial_path);
        return;
    end
    
    % Extract LFP trace - LEGACY EXACT
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_HP')
        lfp_trace = FPA.ephys.lfp_raw_aligned_HP;
        if isstruct(lfp_trace)
            warning('lfp_raw_aligned_HP is a struct in %s', trial_path);
            return;
        end
        lfp_trace = lfp_trace(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'lfp_raw_aligned_mPFC')
        lfp_trace = FPA.ephys.lfp_raw_aligned_mPFC;
        if isstruct(lfp_trace)
            warning('lfp_raw_aligned_mPFC is a struct in %s', trial_path);
            return;
        end
        lfp_trace = lfp_trace(:);
    else
        warning('LFP trace not found in %s', trial_path);
        return;
    end
    
    % Extract motion trace - LEGACY EXACT
    if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
        motion_raw = FPA.ephys.running_velocity_smooth;
        if isstruct(motion_raw)
            warning('running_velocity_smooth is a struct in %s', trial_path);
            return;
        end
        motion_raw = motion_raw(:);
    elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
        motion_raw = FPA.ephys.running_velocity;
        if isstruct(motion_raw)
            warning('running_velocity is a struct in %s', trial_path);
            return;
        end
        motion_raw = motion_raw(:);
    else
        warning('Motion trace not found in %s', trial_path);
        return;
    end
    
    % Align lengths (in case of minor mismatches) - LEGACY EXACT
    n = min([length(t), length(fiber_trace), length(lfp_trace), length(motion_raw)]);
    t = t(1:n);
    fiber_trace = fiber_trace(1:n);
    lfp_trace = lfp_trace(1:n);
    motion_raw = motion_raw(1:n);
    
    % Convert motion to speed (cm/s) - LEGACY EXACT
    speed_cm_s = motion_raw * MOTION_TO_CM_PER_S;
    
    % Optional smoothing - LEGACY EXACT
    if MOTION_SMOOTH_SAMPLES > 1
        kernel = ones(MOTION_SMOOTH_SAMPLES, 1) / MOTION_SMOOTH_SAMPLES;
        speed_cm_s = conv(speed_cm_s, kernel, 'same');
    end
    
    % Store data (matching legacy field names where applicable)
    data.time = t;
    data.fs = fs;
    data.gevi = fiber_trace;  % Called 'fiber' in legacy, 'gevi' here for clarity
    data.lfp = lfp_trace;
    data.speed = speed_cm_s;
    data.duration = t(end) - t(1);
    
catch ME
    warning('Failed to load %s: %s', trial_path, ME.message);
    data = [];
end

end

%% ============================================================================
%  SINGLE-TRIAL COMPUTATION (with spectrograms and time-resolved coherence)
%  ============================================================================

function trial_out = compute_single_trial_mscohere(data, is_rest, is_run, cfg, trial_label, trial_idx, behavior_info)
%COMPUTE_SINGLE_TRIAL_MSCOHERE Complete single-trial analysis with spectrograms
%  Matches legacy script output structure exactly

trial_out = struct();
fs = data.fs;

% Spectrogram parameters
spec_window_samples = round(cfg.spectrogram.window_sec * fs);
spec_overlap_samples = round(cfg.spectrogram.overlap_frac * spec_window_samples);
spec_nfft = 2^nextpow2(spec_window_samples * 2);

% Coherence parameters
coh_segment_samples = round(cfg.coherence.mscohere.segment_sec * fs);
coh_overlap_samples = round(cfg.coherence.mscohere.overlap_frac * coh_segment_samples);
coh_nfft = cfg.coherence.mscohere.nfft_factor * coh_segment_samples;
freq_min = cfg.coherence.mscohere.freq_min;
freq_max = cfg.coherence.mscohere.freq_max;

% ============================================================================
%  LFP SPECTROGRAM
% ============================================================================
[S, F_spec, T_spec] = spectrogram(data.lfp, spec_window_samples, ...
    spec_overlap_samples, spec_nfft, fs);

% Convert to power spectral density (µV²/Hz)
spec_power = abs(S).^2 / (fs * spec_window_samples);

% Limit frequency range
freq_idx = (F_spec >= freq_min) & (F_spec <= freq_max);
F_spec = F_spec(freq_idx);
spec_power = spec_power(freq_idx, :);

% Apply smoothing (2D smoothing function)
spec_power_smoothed = smooth2a(spec_power, cfg.spectrogram.smooth_freq, cfg.spectrogram.smooth_time);

% Shift time to align with original time vector
T_spec = T_spec + data.time(1);

trial_out.spec_power = spec_power_smoothed;  % [freq × time]
trial_out.spec_freq = F_spec(:);              % [freq × 1]
trial_out.spec_time = T_spec(:)';             % [1 × time]

% ============================================================================
%  OVERALL COHERENCE SPECTRUM
% ============================================================================
window = hanning(coh_segment_samples);
[Cxy, F_coh] = mscohere(data.lfp, data.gevi, window, coh_overlap_samples, coh_nfft, fs);
[Pxy, ~] = cpsd(data.lfp, data.gevi, window, coh_overlap_samples, coh_nfft, fs);
coh_phase = angle(Pxy);

% Limit frequency range
freq_idx_coh = (F_coh >= freq_min) & (F_coh <= freq_max);
F_coh = F_coh(freq_idx_coh);
Cxy = Cxy(freq_idx_coh);
coh_phase = coh_phase(freq_idx_coh);

trial_out.coh_mag = Cxy(:);        % [freq × 1]
trial_out.coh_phase = coh_phase(:); % [freq × 1]
trial_out.coh_freq = F_coh(:);      % [freq × 1]

% ============================================================================
%  TIME-RESOLVED COHERENCE (sliding window)
% ============================================================================
time_coh_window_sec = 5.0;   % 5-second windows
time_coh_step_sec = 0.25;    % 0.25-second step

time_coh_window_samples = round(time_coh_window_sec * fs);
time_coh_step_samples = round(time_coh_step_sec * fs);

n_samples = length(data.lfp);
n_time_windows = floor((n_samples - time_coh_window_samples) / time_coh_step_samples) + 1;

coh_time_resolved = zeros(length(F_coh), n_time_windows);
coh_time_centers = zeros(1, n_time_windows);

for win_idx = 1:n_time_windows
    start_sample = (win_idx - 1) * time_coh_step_samples + 1;
    end_sample = start_sample + time_coh_window_samples - 1;
    
    lfp_win = data.lfp(start_sample:end_sample);
    gevi_win = data.gevi(start_sample:end_sample);
    
    [Cxy_win, ~] = mscohere(lfp_win, gevi_win, window, coh_overlap_samples, coh_nfft, fs);
    coh_time_resolved(:, win_idx) = Cxy_win(freq_idx_coh);
    coh_time_centers(win_idx) = data.time(1) + (start_sample + time_coh_window_samples/2 - 1) / fs;
end

% Apply smoothing (matches legacy: minimal smoothing in freq and time)
coh_time_resolved_smoothed = smooth2a(coh_time_resolved, 1, 1);

trial_out.coh_time_resolved = coh_time_resolved_smoothed;  % [freq × time]
trial_out.coh_time = coh_time_centers(:)';                  % [1 × time]

% ============================================================================
%  REST vs RUN COHERENCE SPECTRA
% ============================================================================
min_coh_bout_samples = round(0.1 * fs);  % Minimum 0.1s for coherence

% Find continuous bouts
rest_idx = find(is_rest);
rest_bouts = find_continuous_bouts(rest_idx, min_coh_bout_samples);

coh_rest = [];
coh_rest_freq = [];
if ~isempty(rest_bouts)
    lfp_rest = [];
    gevi_rest = [];
    for bout = 1:size(rest_bouts, 1)
        lfp_rest = [lfp_rest; data.lfp(rest_bouts(bout,1):rest_bouts(bout,2))]; %#ok<AGROW>
        gevi_rest = [gevi_rest; data.gevi(rest_bouts(bout,1):rest_bouts(bout,2))]; %#ok<AGROW>
    end
    if length(lfp_rest) >= coh_segment_samples * 2
        [Cxy_rest, F_rest] = mscohere(lfp_rest, gevi_rest, window, coh_overlap_samples, coh_nfft, fs);
        coh_rest = Cxy_rest(freq_idx_coh);
        coh_rest_freq = F_rest(freq_idx_coh);
    end
end

run_idx = find(is_run);
run_bouts = find_continuous_bouts(run_idx, min_coh_bout_samples);

coh_run = [];
coh_run_freq = [];
if ~isempty(run_bouts)
    lfp_run = [];
    gevi_run = [];
    for bout = 1:size(run_bouts, 1)
        lfp_run = [lfp_run; data.lfp(run_bouts(bout,1):run_bouts(bout,2))]; %#ok<AGROW>
        gevi_run = [gevi_run; data.gevi(run_bouts(bout,1):run_bouts(bout,2))]; %#ok<AGROW>
    end
    if length(lfp_run) >= coh_segment_samples * 2
        [Cxy_run, F_run] = mscohere(lfp_run, gevi_run, window, coh_overlap_samples, coh_nfft, fs);
        coh_run = Cxy_run(freq_idx_coh);
        coh_run_freq = F_run(freq_idx_coh);
    end
end

trial_out.coh_rest = coh_rest(:);
trial_out.coh_rest_freq = coh_rest_freq(:);
trial_out.coh_run = coh_run(:);
trial_out.coh_run_freq = coh_run_freq(:);

% ============================================================================
%  PSD SPECTRA
% ============================================================================
psd_segment_samples = round(cfg.psd.window_sec * fs);
psd_overlap_samples = round(cfg.psd.overlap_frac * psd_segment_samples);
psd_nfft = 2^nextpow2(psd_segment_samples * 2);

[Pxx_lfp, F_psd] = pwelch(data.lfp, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
[Pxx_gevi, ~] = pwelch(data.gevi, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);

freq_idx_psd = (F_psd >= freq_min) & (F_psd <= freq_max);
F_psd = F_psd(freq_idx_psd);
psd_lfp = 10*log10(Pxx_lfp(freq_idx_psd) + eps);
psd_gevi = 10*log10(Pxx_gevi(freq_idx_psd) + eps);

trial_out.psd_freq = F_psd(:);
trial_out.psd_lfp = psd_lfp(:);
trial_out.psd_gevi = psd_gevi(:);

% REST PSD
psd_lfp_rest = [];
psd_gevi_rest = [];
psd_rest_freq = [];
if ~isempty(rest_bouts) && length(lfp_rest) >= psd_segment_samples * 2
    [Pxx_lfp_rest, ~] = pwelch(lfp_rest, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    [Pxx_gevi_rest, ~] = pwelch(gevi_rest, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    psd_lfp_rest = 10*log10(Pxx_lfp_rest(freq_idx_psd) + eps);
    psd_gevi_rest = 10*log10(Pxx_gevi_rest(freq_idx_psd) + eps);
    psd_rest_freq = F_psd;
end

% RUN PSD
psd_lfp_run = [];
psd_gevi_run = [];
psd_run_freq = [];
if ~isempty(run_bouts) && length(lfp_run) >= psd_segment_samples * 2
    [Pxx_lfp_run, ~] = pwelch(lfp_run, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    [Pxx_gevi_run, ~] = pwelch(gevi_run, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    psd_lfp_run = 10*log10(Pxx_lfp_run(freq_idx_psd) + eps);
    psd_gevi_run = 10*log10(Pxx_gevi_run(freq_idx_psd) + eps);
    psd_run_freq = F_psd;
end

trial_out.psd_lfp_rest = psd_lfp_rest(:);
trial_out.psd_gevi_rest = psd_gevi_rest(:);
trial_out.psd_rest_freq = psd_rest_freq(:);
trial_out.psd_lfp_run = psd_lfp_run(:);
trial_out.psd_gevi_run = psd_gevi_run(:);
trial_out.psd_run_freq = psd_run_freq(:);

% ============================================================================
%  METADATA
% ============================================================================
trial_out.t = data.time(:)';  % [1 × time] - row vector
trial_out.fs = fs;
trial_out.speed = data.speed(:)';  % [1 × time] - row vector
trial_out.duration = data.time(end) - data.time(1);
trial_out.method = 'mscohere';
trial_out.trial_label = trial_label;
trial_out.trial_id = trial_idx;
trial_out.behavior = behavior_info;  % Include behavior classification info

end

function trial_out = compute_single_trial_fieldtrip(data, is_rest, is_run, cfg, trial_label, trial_idx, behavior_info)
%COMPUTE_SINGLE_TRIAL_FIELDTRIP Complete single-trial analysis with FieldTrip
%  Uses FieldTrip's ft_freqanalysis and ft_connectivityanalysis for coherence
%  Matches legacy script approach exactly

trial_out = struct();
fs = data.fs;

% Check if FieldTrip is available
if ~exist('ft_defaults', 'file')
    warning('FieldTrip not found, falling back to mscohere for trial %d', trial_idx);
    trial_out = compute_single_trial_mscohere(data, is_rest, is_run, cfg, trial_label, trial_idx, behavior_info);
    trial_out.method = 'fieldtrip';
    return;
end

ft_defaults;

% Spectrogram parameters (same as mscohere for consistency)
spec_window_samples = round(cfg.spectrogram.window_sec * fs);
spec_overlap_samples = round(cfg.spectrogram.overlap_frac * spec_window_samples);
spec_nfft = 2^nextpow2(spec_window_samples * 2);

% Coherence parameters
freq_min = cfg.coherence.fieldtrip.foi_min;
freq_max = cfg.coherence.fieldtrip.foi_max;
tapsmofrq = cfg.coherence.fieldtrip.tapsmofrq;

% Pseudo-trial parameters (for FieldTrip)
epoch_length_sec = cfg.coherence.fieldtrip.pseudotrial_length_sec;
epoch_overlap_sec = cfg.coherence.fieldtrip.pseudotrial_overlap_sec;

% ============================================================================
%  LFP SPECTROGRAM (MATLAB native, same as mscohere for consistency)
% ============================================================================
[S, F_spec, T_spec] = spectrogram(data.lfp, spec_window_samples, ...
    spec_overlap_samples, spec_nfft, fs);

spec_power = abs(S).^2 / (fs * spec_window_samples);
freq_idx = (F_spec >= freq_min) & (F_spec <= freq_max);
F_spec = F_spec(freq_idx);
spec_power = spec_power(freq_idx, :);
spec_power_smoothed = smooth2a(spec_power, cfg.spectrogram.smooth_freq, cfg.spectrogram.smooth_time);
T_spec = T_spec + data.time(1);

trial_out.spec_power = spec_power_smoothed;
trial_out.spec_freq = F_spec(:);
trial_out.spec_time = T_spec(:)';

% ============================================================================
%  FIELDTRIP COHERENCE (Pseudo-Trials Approach)
% ============================================================================
n_samples = length(data.lfp);
epoch_samples = round(epoch_length_sec * fs);
epoch_step = round((epoch_length_sec - epoch_overlap_sec) * fs);
n_epochs = floor((n_samples - epoch_samples) / epoch_step) + 1;

% Create FieldTrip data structure with pseudo-trials
ft_data = [];
ft_data.label = {'LFP'; 'Fiber'};
ft_data.fsample = fs;
ft_data.trial = cell(1, n_epochs);
ft_data.time = cell(1, n_epochs);
ft_data.sampleinfo = zeros(n_epochs, 2);

for ep = 1:n_epochs
    start_sample = (ep - 1) * epoch_step + 1;
    end_sample = start_sample + epoch_samples - 1;
    
    ft_data.trial{ep} = [data.lfp(start_sample:end_sample)'; ...
                         data.gevi(start_sample:end_sample)'];
    ft_data.time{ep} = (0:epoch_samples-1) / fs;
    ft_data.sampleinfo(ep, :) = [start_sample, end_sample];
end

% Time-frequency decomposition - taper is user-configurable
cfg_freq = [];
cfg_freq.method = 'mtmfft';
cfg_freq.output = 'fourier';
cfg_freq.taper = cfg.coherence.fieldtrip.taper;  % 'dpss' or 'hanning'
cfg_freq.foi = freq_min:cfg.coherence.fieldtrip.foi_step:freq_max;
cfg_freq.keeptrials = 'yes';
cfg_freq.pad = 'nextpow2';

% tapsmofrq only applies to 'dpss' (multi-taper)
if strcmpi(cfg.coherence.fieldtrip.taper, 'dpss')
    cfg_freq.tapsmofrq = tapsmofrq;
end

freq_result = ft_freqanalysis(cfg_freq, ft_data);

% Compute coherence
cfg_conn = [];
cfg_conn.method = 'coh';
cfg_conn.complex = 'abs';

conn_result = ft_connectivityanalysis(cfg_conn, freq_result);

% Extract LFP-Fiber coherence
coh_matrix = conn_result.cohspctrm;
coh_lfp_fiber = squeeze(coh_matrix(1, 2, :));
F_coh = conn_result.freq(:);

trial_out.coh_mag = coh_lfp_fiber(:);
trial_out.coh_freq = F_coh(:);

% ============================================================================
%  TIME-RESOLVED COHERENCE (Sliding Window FieldTrip)
% ============================================================================
% Use same parameters as legacy code for full time coverage
time_window_epochs = 10;  % Number of consecutive epochs per time point (10 × 1s = 10s effective window)

% IMPORTANT: Legacy code steps by 1 epoch (not 2), giving full time coverage
% This produces ~290 time points for a 60s recording vs ~145 with step=2
n_time_points = n_epochs - time_window_epochs + 1;
coh_time_resolved = zeros(length(F_coh), n_time_points);
coh_time_centers = zeros(1, n_time_points);

for t_idx = 1:n_time_points
    epoch_start = t_idx;  % Step by 1 epoch (matches legacy)
    epoch_end = t_idx + time_window_epochs - 1;
    
    % Select subset of trials
    cfg_sel = [];
    cfg_sel.trials = epoch_start:epoch_end;
    ft_data_subset = ft_selectdata(cfg_sel, ft_data);
    
    % Frequency analysis on subset
    freq_subset = ft_freqanalysis(cfg_freq, ft_data_subset);
    
    % Coherence on subset
    conn_subset = ft_connectivityanalysis(cfg_conn, freq_subset);
    
    % Store
    coh_time_resolved(:, t_idx) = squeeze(conn_subset.cohspctrm(1, 2, :));
    
    % Time center (middle of the window)
    middle_epoch = epoch_start + floor(time_window_epochs / 2);
    center_sample = ft_data.sampleinfo(middle_epoch, 1) + epoch_samples / 2;
    coh_time_centers(t_idx) = data.time(1) + (center_sample - 1) / fs;
end

% Apply smoothing (matches legacy: minimal smoothing in freq and time)
coh_time_resolved_smoothed = smooth2a(coh_time_resolved, 1, 1);

trial_out.coh_time_resolved = coh_time_resolved_smoothed;
trial_out.coh_time = coh_time_centers(:)';

% ============================================================================
%  REST vs RUN COHERENCE (using mscohere like legacy script)
% ============================================================================
% Legacy script uses mscohere for REST/RUN separation for consistency
% We'll do the same here

min_coh_bout_samples = round(0.1 * fs);
coh_segment_samples = round(cfg.coherence.mscohere.segment_sec * fs);
coh_overlap_samples = round(cfg.coherence.mscohere.overlap_frac * coh_segment_samples);
coh_nfft = cfg.coherence.mscohere.nfft_factor * coh_segment_samples;
window = hanning(coh_segment_samples);

% REST coherence
rest_idx = find(is_rest);
rest_bouts = find_continuous_bouts(rest_idx, min_coh_bout_samples);

coh_rest = [];
coh_rest_freq = [];
if ~isempty(rest_bouts)
    lfp_rest = [];
    gevi_rest = [];
    for bout = 1:size(rest_bouts, 1)
        lfp_rest = [lfp_rest; data.lfp(rest_bouts(bout,1):rest_bouts(bout,2))]; %#ok<AGROW>
        gevi_rest = [gevi_rest; data.gevi(rest_bouts(bout,1):rest_bouts(bout,2))]; %#ok<AGROW>
    end
    if length(lfp_rest) >= coh_segment_samples * 2
        [Cxy_rest, F_rest] = mscohere(lfp_rest, gevi_rest, window, coh_overlap_samples, coh_nfft, fs);
        freq_mask = (F_rest >= freq_min) & (F_rest <= freq_max);
        coh_rest = Cxy_rest(freq_mask);
        coh_rest_freq = F_rest(freq_mask);
    end
end

% RUN coherence
run_idx = find(is_run);
run_bouts = find_continuous_bouts(run_idx, min_coh_bout_samples);

coh_run = [];
coh_run_freq = [];
if ~isempty(run_bouts)
    lfp_run = [];
    gevi_run = [];
    for bout = 1:size(run_bouts, 1)
        lfp_run = [lfp_run; data.lfp(run_bouts(bout,1):run_bouts(bout,2))]; %#ok<AGROW>
        gevi_run = [gevi_run; data.gevi(run_bouts(bout,1):run_bouts(bout,2))]; %#ok<AGROW>
    end
    if length(lfp_run) >= coh_segment_samples * 2
        [Cxy_run, F_run] = mscohere(lfp_run, gevi_run, window, coh_overlap_samples, coh_nfft, fs);
        freq_mask = (F_run >= freq_min) & (F_run <= freq_max);
        coh_run = Cxy_run(freq_mask);
        coh_run_freq = F_run(freq_mask);
    end
end

trial_out.coh_rest = coh_rest(:);
trial_out.coh_rest_freq = coh_rest_freq(:);
trial_out.coh_run = coh_run(:);
trial_out.coh_run_freq = coh_run_freq(:);

% ============================================================================
%  PSD SPECTRA (same as mscohere)
% ============================================================================
psd_segment_samples = round(cfg.psd.window_sec * fs);
psd_overlap_samples = round(cfg.psd.overlap_frac * psd_segment_samples);
psd_nfft = 2^nextpow2(psd_segment_samples * 2);

[Pxx_lfp, F_psd] = pwelch(data.lfp, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
[Pxx_gevi, ~] = pwelch(data.gevi, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);

freq_idx_psd = (F_psd >= freq_min) & (F_psd <= freq_max);
F_psd = F_psd(freq_idx_psd);
psd_lfp = 10*log10(Pxx_lfp(freq_idx_psd) + eps);
psd_gevi = 10*log10(Pxx_gevi(freq_idx_psd) + eps);

trial_out.psd_freq = F_psd(:);
trial_out.psd_lfp = psd_lfp(:);
trial_out.psd_gevi = psd_gevi(:);

% REST PSD
psd_lfp_rest = [];
psd_gevi_rest = [];
psd_rest_freq = [];
if ~isempty(rest_bouts) && length(lfp_rest) >= psd_segment_samples * 2
    [Pxx_lfp_rest, ~] = pwelch(lfp_rest, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    [Pxx_gevi_rest, ~] = pwelch(gevi_rest, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    psd_lfp_rest = 10*log10(Pxx_lfp_rest(freq_idx_psd) + eps);
    psd_gevi_rest = 10*log10(Pxx_gevi_rest(freq_idx_psd) + eps);
    psd_rest_freq = F_psd;
end

% RUN PSD
psd_lfp_run = [];
psd_gevi_run = [];
psd_run_freq = [];
if ~isempty(run_bouts) && length(lfp_run) >= psd_segment_samples * 2
    [Pxx_lfp_run, ~] = pwelch(lfp_run, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    [Pxx_gevi_run, ~] = pwelch(gevi_run, hanning(psd_segment_samples), psd_overlap_samples, psd_nfft, fs);
    psd_lfp_run = 10*log10(Pxx_lfp_run(freq_idx_psd) + eps);
    psd_gevi_run = 10*log10(Pxx_gevi_run(freq_idx_psd) + eps);
    psd_run_freq = F_psd;
end

trial_out.psd_lfp_rest = psd_lfp_rest(:);
trial_out.psd_gevi_rest = psd_gevi_rest(:);
trial_out.psd_rest_freq = psd_rest_freq(:);
trial_out.psd_lfp_run = psd_lfp_run(:);
trial_out.psd_gevi_run = psd_gevi_run(:);
trial_out.psd_run_freq = psd_run_freq(:);

% ============================================================================
%  METADATA
% ============================================================================
trial_out.t = data.time(:)';
trial_out.fs = fs;
trial_out.speed = data.speed(:)';
trial_out.duration = data.time(end) - data.time(1);
trial_out.method = 'fieldtrip';
trial_out.trial_label = trial_label;
trial_out.trial_id = trial_idx;
trial_out.behavior = behavior_info;

end

%% ============================================================================
%  HELPER: Compute Spectra (mscohere) - for pooled analysis
%  ============================================================================
function spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg)
%COMPUTE_SPECTRA_MSCOHERE Compute coherence and PSD using mscohere method

spectra = struct();
fs = data.fs;

% Parameters
segment_samples = round(cfg.coherence.mscohere.segment_sec * fs);
overlap_samples = round(cfg.coherence.mscohere.overlap_frac * segment_samples);
nfft = cfg.coherence.mscohere.nfft_factor * segment_samples;
freq_min = cfg.coherence.mscohere.freq_min;
freq_max = cfg.coherence.mscohere.freq_max;

% Window
window = hanning(segment_samples);

% Overall coherence
[coh_overall, freq] = mscohere(data.lfp, data.gevi, window, overlap_samples, nfft, fs);

% Frequency mask
freq_mask = (freq >= freq_min) & (freq <= freq_max);
freq = freq(freq_mask);
coh_overall = coh_overall(freq_mask);

% REST coherence
if sum(is_rest) > segment_samples
    lfp_rest = data.lfp(is_rest);
    gevi_rest = data.gevi(is_rest);
    coh_rest = mscohere(lfp_rest, gevi_rest, window, overlap_samples, nfft, fs);
    coh_rest = coh_rest(freq_mask);
else
    coh_rest = [];
end

% RUN coherence
if sum(is_run) > segment_samples
    lfp_run = data.lfp(is_run);
    gevi_run = data.gevi(is_run);
    coh_run = mscohere(lfp_run, gevi_run, window, overlap_samples, nfft, fs);
    coh_run = coh_run(freq_mask);
else
    coh_run = [];
end

% Store coherence
spectra.coherence.overall = coh_overall(:);
spectra.coherence.rest = coh_rest(:);
spectra.coherence.run = coh_run(:);
spectra.coherence.freq = freq(:);

% PSD - LFP
[psd_lfp_overall, freq_psd] = pwelch(data.lfp, window, overlap_samples, nfft, fs);
freq_mask_psd = (freq_psd >= freq_min) & (freq_psd <= freq_max);
freq_psd = freq_psd(freq_mask_psd);
psd_lfp_overall = 10*log10(psd_lfp_overall(freq_mask_psd));

% PSD - GEVI
psd_gevi_overall = pwelch(data.gevi, window, overlap_samples, nfft, fs);
psd_gevi_overall = 10*log10(psd_gevi_overall(freq_mask_psd));

% REST PSD
if sum(is_rest) > segment_samples
    psd_lfp_rest = pwelch(data.lfp(is_rest), window, overlap_samples, nfft, fs);
    psd_lfp_rest = 10*log10(psd_lfp_rest(freq_mask_psd));
    psd_gevi_rest = pwelch(data.gevi(is_rest), window, overlap_samples, nfft, fs);
    psd_gevi_rest = 10*log10(psd_gevi_rest(freq_mask_psd));
else
    psd_lfp_rest = [];
    psd_gevi_rest = [];
end

% RUN PSD
if sum(is_run) > segment_samples
    psd_lfp_run = pwelch(data.lfp(is_run), window, overlap_samples, nfft, fs);
    psd_lfp_run = 10*log10(psd_lfp_run(freq_mask_psd));
    psd_gevi_run = pwelch(data.gevi(is_run), window, overlap_samples, nfft, fs);
    psd_gevi_run = 10*log10(psd_gevi_run(freq_mask_psd));
else
    psd_lfp_run = [];
    psd_gevi_run = [];
end

% Store PSD
spectra.psd_lfp.overall = psd_lfp_overall(:);
spectra.psd_lfp.rest = psd_lfp_rest(:);
spectra.psd_lfp.run = psd_lfp_run(:);
spectra.psd_gevi.overall = psd_gevi_overall(:);
spectra.psd_gevi.rest = psd_gevi_rest(:);
spectra.psd_gevi.run = psd_gevi_run(:);
spectra.psd_freq = freq_psd(:);

end

%% ============================================================================
%  HELPER: Compute Spectra (FieldTrip)
%  ============================================================================
function spectra = compute_spectra_fieldtrip(data, is_rest, is_run, cfg)
%COMPUTE_SPECTRA_FIELDTRIP Compute coherence and PSD using FieldTrip
%  CRITICAL: FieldTrip coherence requires multiple pseudo-trials.
%  This function segments continuous data into overlapping epochs,
%  matching the legacy approach in figure2_coherence_fieldtrip.m

spectra = struct();

% Check FieldTrip availability
if ~exist('ft_freqanalysis', 'file')
    warning('FieldTrip not available, falling back to mscohere');
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

fs = data.fs;

% Pseudo-trial parameters - matching legacy script
EPOCH_LENGTH_SEC = cfg.coherence.fieldtrip.pseudotrial_length_sec;  % 2.0s default
EPOCH_OVERLAP_SEC = cfg.coherence.fieldtrip.pseudotrial_overlap_sec; % 1.0s default

% -------------------------------------------------------------------------
% HELPER: Create FieldTrip data structure with pseudo-trials
% -------------------------------------------------------------------------
    function ft_data = create_ft_data_pseudotrials(lfp_signal, gevi_signal, fs_local, epoch_len_sec, epoch_ovlp_sec)
        n_samples_local = length(lfp_signal);
        epoch_samples = round(epoch_len_sec * fs_local);
        epoch_step = round((epoch_len_sec - epoch_ovlp_sec) * fs_local);
        
        n_epochs = floor((n_samples_local - epoch_samples) / epoch_step) + 1;
        
        if n_epochs < 1
            ft_data = [];
            return;
        end
        
        ft_data = struct();
        ft_data.label = {'LFP'; 'GEVI'};
        ft_data.fsample = fs_local;
        ft_data.trial = cell(1, n_epochs);
        ft_data.time = cell(1, n_epochs);
        ft_data.sampleinfo = zeros(n_epochs, 2);
        
        for ep = 1:n_epochs
            start_sample = (ep - 1) * epoch_step + 1;
            end_sample = start_sample + epoch_samples - 1;
            
            if end_sample > n_samples_local
                break;
            end
            
            ft_data.trial{ep} = [lfp_signal(start_sample:end_sample)'; ...
                                 gevi_signal(start_sample:end_sample)'];
            ft_data.time{ep} = (0:epoch_samples-1) / fs_local;
            ft_data.sampleinfo(ep, :) = [start_sample, end_sample];
        end
    end

% -------------------------------------------------------------------------
% OVERALL coherence (all data)
% -------------------------------------------------------------------------
ft_data_overall = create_ft_data_pseudotrials(data.lfp, data.gevi, fs, ...
    EPOCH_LENGTH_SEC, EPOCH_OVERLAP_SEC);

if isempty(ft_data_overall)
    warning('Insufficient data for FieldTrip pseudo-trials, falling back to mscohere');
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

% Frequency analysis config - CRITICAL: keeptrials='yes', output='fourier'
ft_cfg = [];
ft_cfg.method = cfg.coherence.fieldtrip.method;  % 'mtmfft'
ft_cfg.taper = cfg.coherence.fieldtrip.taper;    % 'dpss'
ft_cfg.foi = cfg.coherence.fieldtrip.foi_min:cfg.coherence.fieldtrip.foi_step:cfg.coherence.fieldtrip.foi_max;
ft_cfg.keeptrials = 'yes';   % MUST keep trials for coherence
ft_cfg.output = 'fourier';   % Complex Fourier spectra for coherence
ft_cfg.pad = 'nextpow2';

if strcmpi(cfg.coherence.fieldtrip.taper, 'dpss')
    ft_cfg.tapsmofrq = cfg.coherence.fieldtrip.tapsmofrq;
end

% Coherence config
coh_cfg = [];
coh_cfg.method = 'coh';
coh_cfg.complex = 'abs';

try
    freq_overall = ft_freqanalysis(ft_cfg, ft_data_overall);
    coh_overall = ft_connectivityanalysis(coh_cfg, freq_overall);
    
    % cohspctrm is [chan × chan × freq] after trial averaging
    spectra.coherence.overall = squeeze(coh_overall.cohspctrm(1,2,:));
    spectra.coherence.freq = coh_overall.freq(:);
    
    % PSD from power spectrum (average over trials)
    % powspctrm dims: [ntrials × nchan × nfreq] when keeptrials='yes'
    % Need to first get power, then average
    ft_cfg_pow = ft_cfg;
    ft_cfg_pow.output = 'pow';
    ft_cfg_pow.keeptrials = 'no';  % Average for PSD
    freq_pow = ft_freqanalysis(ft_cfg_pow, ft_data_overall);
    
    spectra.psd_lfp.overall = 10*log10(squeeze(freq_pow.powspctrm(1,:)))';
    spectra.psd_gevi.overall = 10*log10(squeeze(freq_pow.powspctrm(2,:)))';
    spectra.psd_freq = freq_pow.freq(:);
    
catch ME
    warning('FieldTrip analysis failed: %s. Using mscohere.', ME.message);
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

% -------------------------------------------------------------------------
% REST coherence
% -------------------------------------------------------------------------
rest_samples = sum(is_rest);
if rest_samples > fs * EPOCH_LENGTH_SEC * 2  % At least 2 epochs worth
    ft_data_rest = create_ft_data_pseudotrials(data.lfp(is_rest), ...
        data.gevi(is_rest), fs, EPOCH_LENGTH_SEC, EPOCH_OVERLAP_SEC);
    
    if ~isempty(ft_data_rest) && length(ft_data_rest.trial) >= 2
        try
            freq_rest = ft_freqanalysis(ft_cfg, ft_data_rest);
            coh_rest = ft_connectivityanalysis(coh_cfg, freq_rest);
            spectra.coherence.rest = squeeze(coh_rest.cohspctrm(1,2,:));
            
            % PSD for rest
            freq_pow_rest = ft_freqanalysis(ft_cfg_pow, ft_data_rest);
            spectra.psd_lfp.rest = 10*log10(squeeze(freq_pow_rest.powspctrm(1,:)))';
            spectra.psd_gevi.rest = 10*log10(squeeze(freq_pow_rest.powspctrm(2,:)))';
        catch
            spectra.coherence.rest = [];
            spectra.psd_lfp.rest = [];
            spectra.psd_gevi.rest = [];
        end
    else
        spectra.coherence.rest = [];
        spectra.psd_lfp.rest = [];
        spectra.psd_gevi.rest = [];
    end
else
    spectra.coherence.rest = [];
    spectra.psd_lfp.rest = [];
    spectra.psd_gevi.rest = [];
end

% -------------------------------------------------------------------------
% RUN coherence
% -------------------------------------------------------------------------
run_samples = sum(is_run);
if run_samples > fs * EPOCH_LENGTH_SEC * 2  % At least 2 epochs worth
    ft_data_run = create_ft_data_pseudotrials(data.lfp(is_run), ...
        data.gevi(is_run), fs, EPOCH_LENGTH_SEC, EPOCH_OVERLAP_SEC);
    
    if ~isempty(ft_data_run) && length(ft_data_run.trial) >= 2
        try
            freq_run = ft_freqanalysis(ft_cfg, ft_data_run);
            coh_run = ft_connectivityanalysis(coh_cfg, freq_run);
            spectra.coherence.run = squeeze(coh_run.cohspctrm(1,2,:));
            
            % PSD for run
            freq_pow_run = ft_freqanalysis(ft_cfg_pow, ft_data_run);
            spectra.psd_lfp.run = 10*log10(squeeze(freq_pow_run.powspctrm(1,:)))';
            spectra.psd_gevi.run = 10*log10(squeeze(freq_pow_run.powspctrm(2,:)))';
        catch
            spectra.coherence.run = [];
            spectra.psd_lfp.run = [];
            spectra.psd_gevi.run = [];
        end
    else
        spectra.coherence.run = [];
        spectra.psd_lfp.run = [];
        spectra.psd_gevi.run = [];
    end
else
    spectra.coherence.run = [];
    spectra.psd_lfp.run = [];
    spectra.psd_gevi.run = [];
end

end

%% ============================================================================
%  HELPER: Compute Spectra (FieldTrip) WITH BOUNDARY AWARENESS
%  ============================================================================
function spectra = compute_spectra_fieldtrip_with_boundaries(data, is_rest, is_run, cfg)
%COMPUTE_SPECTRA_FIELDTRIP_WITH_BOUNDARIES Compute coherence/PSD using FieldTrip
%  with proper handling of trial/session boundaries to avoid edge artifacts.
%
%  This function implements the boundary-aware approach from the legacy script
%  figure2_coherence_fieldtrip_pooledtrials.m:
%    1. Creates pseudo-trials from concatenated data
%    2. SKIPS pseudo-trials that cross any trial/session boundary
%    3. For REST/RUN: selects pseudo-trials by behavioral state (not raw indexing)
%
%  INPUTS:
%    data     - Struct with fields: lfp, gevi, speed, time, fs, trial_boundaries
%    is_rest  - Logical vector for REST classification
%    is_run   - Logical vector for RUN classification
%    cfg      - Configuration struct
%
%  OUTPUT:
%    spectra  - Struct with coherence and PSD results

spectra = struct();

% Check FieldTrip availability
if ~exist('ft_freqanalysis', 'file')
    warning('FieldTrip not available, falling back to mscohere');
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

fs = data.fs;
n_samples = length(data.lfp);

% Get trial boundaries (default to no boundaries if not provided)
if isfield(data, 'trial_boundaries') && ~isempty(data.trial_boundaries)
    trial_boundaries = data.trial_boundaries(:);
else
    trial_boundaries = [0; n_samples];  % No internal boundaries
end

% Pseudo-trial parameters
EPOCH_LENGTH_SEC = cfg.coherence.fieldtrip.pseudotrial_length_sec;
EPOCH_OVERLAP_SEC = cfg.coherence.fieldtrip.pseudotrial_overlap_sec;

epoch_samples = round(EPOCH_LENGTH_SEC * fs);
epoch_step = round((EPOCH_LENGTH_SEC - EPOCH_OVERLAP_SEC) * fs);

% Helper function to check if a range crosses any trial boundary
crosses_boundary = @(start_idx, end_idx) any(trial_boundaries > start_idx & trial_boundaries < end_idx);

% -------------------------------------------------------------------------
% First pass: identify valid pseudo-trials (those not crossing boundaries)
% -------------------------------------------------------------------------
n_potential_epochs = floor((n_samples - epoch_samples) / epoch_step) + 1;
valid_pt_indices = [];

for pt = 1:n_potential_epochs
    start_idx = (pt - 1) * epoch_step + 1;
    end_idx = start_idx + epoch_samples - 1;
    
    if end_idx <= n_samples && ~crosses_boundary(start_idx, end_idx)
        valid_pt_indices = [valid_pt_indices; pt]; %#ok<AGROW>
    end
end

n_valid_epochs = length(valid_pt_indices);
n_skipped = n_potential_epochs - n_valid_epochs;

fprintf('          FieldTrip: %d pseudo-trials (%.1fs each, %.0f%% overlap)\n', ...
    n_valid_epochs, EPOCH_LENGTH_SEC, (EPOCH_OVERLAP_SEC/EPOCH_LENGTH_SEC)*100);
if n_skipped > 0
    fprintf('          Skipped %d pseudo-trials crossing trial boundaries\n', n_skipped);
end

if n_valid_epochs < 2
    warning('Insufficient valid pseudo-trials for FieldTrip coherence, falling back to mscohere');
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

% -------------------------------------------------------------------------
% Build FieldTrip data structure with valid pseudo-trials only
% -------------------------------------------------------------------------
ft_data = struct();
ft_data.label = {'LFP'; 'GEVI'};
ft_data.fsample = fs;
ft_data.trial = cell(1, n_valid_epochs);
ft_data.time = cell(1, n_valid_epochs);
ft_data.sampleinfo = zeros(n_valid_epochs, 2);

for i = 1:n_valid_epochs
    pt = valid_pt_indices(i);
    start_idx = (pt - 1) * epoch_step + 1;
    end_idx = start_idx + epoch_samples - 1;
    
    ft_data.trial{i} = [data.lfp(start_idx:end_idx)'; ...
                        data.gevi(start_idx:end_idx)'];
    ft_data.time{i} = (0:epoch_samples-1) / fs;
    ft_data.sampleinfo(i, :) = [start_idx, end_idx];
end

% -------------------------------------------------------------------------
% FieldTrip frequency analysis configuration
% -------------------------------------------------------------------------
ft_cfg = [];
ft_cfg.method = cfg.coherence.fieldtrip.method;
ft_cfg.taper = cfg.coherence.fieldtrip.taper;
ft_cfg.foi = cfg.coherence.fieldtrip.foi_min:cfg.coherence.fieldtrip.foi_step:cfg.coherence.fieldtrip.foi_max;
ft_cfg.keeptrials = 'yes';
ft_cfg.output = 'fourier';
ft_cfg.pad = 'nextpow2';

if strcmpi(cfg.coherence.fieldtrip.taper, 'dpss')
    ft_cfg.tapsmofrq = cfg.coherence.fieldtrip.tapsmofrq;
end

% Coherence config
coh_cfg = [];
coh_cfg.method = 'coh';
coh_cfg.complex = 'abs';

% PSD config
ft_cfg_pow = ft_cfg;
ft_cfg_pow.output = 'pow';
ft_cfg_pow.keeptrials = 'no';

% -------------------------------------------------------------------------
% OVERALL coherence (all valid pseudo-trials)
% -------------------------------------------------------------------------
try
    freq_overall = ft_freqanalysis(ft_cfg, ft_data);
    coh_overall = ft_connectivityanalysis(coh_cfg, freq_overall);
    
    spectra.coherence.overall = squeeze(coh_overall.cohspctrm(1,2,:));
    spectra.coherence.freq = coh_overall.freq(:);
    
    freq_pow = ft_freqanalysis(ft_cfg_pow, ft_data);
    spectra.psd_lfp.overall = 10*log10(squeeze(freq_pow.powspctrm(1,:)))';
    spectra.psd_gevi.overall = 10*log10(squeeze(freq_pow.powspctrm(2,:)))';
    spectra.psd_freq = freq_pow.freq(:);
    
catch ME
    warning('FieldTrip overall analysis failed: %s. Using mscohere.', ME.message);
    spectra = compute_spectra_mscohere(data, is_rest, is_run, cfg);
    return;
end

% -------------------------------------------------------------------------
% Classify pseudo-trials by behavioral state (REST/RUN)
% Instead of extracting and concatenating raw samples, we select existing
% pseudo-trials that fall predominantly within one behavioral state.
% This avoids creating new discontinuities.
% -------------------------------------------------------------------------
REST_THRESHOLD = 0.90;  % Require 90% of pseudo-trial to be in REST
RUN_THRESHOLD = 0.90;   % Require 90% of pseudo-trial to be in RUN

rest_pseudo_idx = [];
run_pseudo_idx = [];

for i = 1:n_valid_epochs
    start_idx = ft_data.sampleinfo(i, 1);
    end_idx = ft_data.sampleinfo(i, 2);
    
    % Check behavioral state of this pseudo-trial
    pt_is_rest = is_rest(start_idx:end_idx);
    pt_is_run = is_run(start_idx:end_idx);
    
    pct_rest = sum(pt_is_rest) / length(pt_is_rest);
    pct_run = sum(pt_is_run) / length(pt_is_run);
    
    if pct_rest >= REST_THRESHOLD
        rest_pseudo_idx = [rest_pseudo_idx; i]; %#ok<AGROW>
    elseif pct_run >= RUN_THRESHOLD
        run_pseudo_idx = [run_pseudo_idx; i]; %#ok<AGROW>
    end
end

fprintf('          REST pseudo-trials: %d, RUN pseudo-trials: %d\n', ...
    length(rest_pseudo_idx), length(run_pseudo_idx));

% -------------------------------------------------------------------------
% REST coherence (select REST pseudo-trials)
% -------------------------------------------------------------------------
if length(rest_pseudo_idx) >= 2
    ft_rest = struct();
    ft_rest.label = ft_data.label;
    ft_rest.fsample = fs;
    ft_rest.trial = ft_data.trial(rest_pseudo_idx);
    ft_rest.time = ft_data.time(rest_pseudo_idx);
    ft_rest.sampleinfo = ft_data.sampleinfo(rest_pseudo_idx, :);
    
    try
        freq_rest = ft_freqanalysis(ft_cfg, ft_rest);
        coh_rest = ft_connectivityanalysis(coh_cfg, freq_rest);
        spectra.coherence.rest = squeeze(coh_rest.cohspctrm(1,2,:));
        
        freq_pow_rest = ft_freqanalysis(ft_cfg_pow, ft_rest);
        spectra.psd_lfp.rest = 10*log10(squeeze(freq_pow_rest.powspctrm(1,:)))';
        spectra.psd_gevi.rest = 10*log10(squeeze(freq_pow_rest.powspctrm(2,:)))';
    catch
        spectra.coherence.rest = [];
        spectra.psd_lfp.rest = [];
        spectra.psd_gevi.rest = [];
    end
else
    spectra.coherence.rest = [];
    spectra.psd_lfp.rest = [];
    spectra.psd_gevi.rest = [];
end

% -------------------------------------------------------------------------
% RUN coherence (select RUN pseudo-trials)
% -------------------------------------------------------------------------
if length(run_pseudo_idx) >= 2
    ft_run = struct();
    ft_run.label = ft_data.label;
    ft_run.fsample = fs;
    ft_run.trial = ft_data.trial(run_pseudo_idx);
    ft_run.time = ft_data.time(run_pseudo_idx);
    ft_run.sampleinfo = ft_data.sampleinfo(run_pseudo_idx, :);
    
    try
        freq_run = ft_freqanalysis(ft_cfg, ft_run);
        coh_run = ft_connectivityanalysis(coh_cfg, freq_run);
        spectra.coherence.run = squeeze(coh_run.cohspctrm(1,2,:));
        
        freq_pow_run = ft_freqanalysis(ft_cfg_pow, ft_run);
        spectra.psd_lfp.run = 10*log10(squeeze(freq_pow_run.powspctrm(1,:)))';
        spectra.psd_gevi.run = 10*log10(squeeze(freq_pow_run.powspctrm(2,:)))';
    catch
        spectra.coherence.run = [];
        spectra.psd_lfp.run = [];
        spectra.psd_gevi.run = [];
    end
else
    spectra.coherence.run = [];
    spectra.psd_lfp.run = [];
    spectra.psd_gevi.run = [];
end

end

%% ============================================================================
%  HELPER: Average Spectra Across Sessions
%  ============================================================================
function avg_spectra = average_spectra(session_spectra)
%AVERAGE_SPECTRA Average spectral results across sessions

avg_spectra = struct();
n_sessions = length(session_spectra);

if n_sessions == 0
    return;
end

% Get frequency axes from first session
avg_spectra.coherence.freq = session_spectra(1).spectra.coherence.freq;
avg_spectra.psd_freq = session_spectra(1).spectra.psd_freq;

% Collect spectra
coh_overall = [];
coh_rest = [];
coh_run = [];
psd_lfp_overall = [];
psd_lfp_rest = [];
psd_lfp_run = [];
psd_gevi_overall = [];
psd_gevi_rest = [];
psd_gevi_run = [];

for s = 1:n_sessions
    sp = session_spectra(s).spectra;
    
    if ~isempty(sp.coherence.overall)
        coh_overall = [coh_overall, sp.coherence.overall(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.coherence.rest)
        coh_rest = [coh_rest, sp.coherence.rest(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.coherence.run)
        coh_run = [coh_run, sp.coherence.run(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_lfp.overall)
        psd_lfp_overall = [psd_lfp_overall, sp.psd_lfp.overall(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_lfp.rest)
        psd_lfp_rest = [psd_lfp_rest, sp.psd_lfp.rest(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_lfp.run)
        psd_lfp_run = [psd_lfp_run, sp.psd_lfp.run(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_gevi.overall)
        psd_gevi_overall = [psd_gevi_overall, sp.psd_gevi.overall(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_gevi.rest)
        psd_gevi_rest = [psd_gevi_rest, sp.psd_gevi.rest(:)]; %#ok<AGROW>
    end
    if ~isempty(sp.psd_gevi.run)
        psd_gevi_run = [psd_gevi_run, sp.psd_gevi.run(:)]; %#ok<AGROW>
    end
end

% Compute means - handle empty arrays properly
avg_spectra.coherence.overall = safe_mean(coh_overall);
avg_spectra.coherence.rest = safe_mean(coh_rest);
avg_spectra.coherence.run = safe_mean(coh_run);
avg_spectra.psd_lfp.overall = safe_mean(psd_lfp_overall);
avg_spectra.psd_lfp.rest = safe_mean(psd_lfp_rest);
avg_spectra.psd_lfp.run = safe_mean(psd_lfp_run);
avg_spectra.psd_gevi.overall = safe_mean(psd_gevi_overall);
avg_spectra.psd_gevi.rest = safe_mean(psd_gevi_rest);
avg_spectra.psd_gevi.run = safe_mean(psd_gevi_run);

% Store SEM for error bars (only if we have multiple samples)
if ~isempty(coh_overall) && size(coh_overall, 2) > 1
    avg_spectra.coherence.overall_sem = std(coh_overall, 0, 2) / sqrt(size(coh_overall, 2));
end
if ~isempty(coh_rest) && size(coh_rest, 2) > 1
    avg_spectra.coherence.rest_sem = std(coh_rest, 0, 2) / sqrt(size(coh_rest, 2));
end
if ~isempty(coh_run) && size(coh_run, 2) > 1
    avg_spectra.coherence.run_sem = std(coh_run, 0, 2) / sqrt(size(coh_run, 2));
end

end

%% ============================================================================
%  HELPER: Save Single-Trial File (Legacy Format)
%  ============================================================================
function save_single_trial_file(trial_data, animal, session_id, method, trial_idx, cfg)
%SAVE_SINGLE_TRIAL_FILE Save per-trial file in legacy format for Python compatibility
%  Matches exact structure from legacy figure2_coherence_mscohere.m script

if ~cfg.output.save_mat
    return;
end

% Build output path
% Note: output_folder_suffix is already applied to cfg.paths.output_root in pipeline
behavior_mode = cfg.behavior.classification_mode;
output_dir = fullfile(cfg.paths.output_root, behavior_mode, 'single_trial', ...
    animal.mouse_id, session_id, 'data');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Create legacy-format output struct (exact match to legacy script)
fig2_out = struct();

% --- Axes (ensure correct orientation) ---
% Time: 1 × Nt row vector (coherence time axis)
fig2_out.time = trial_data.coh_time(:)';

% Frequency: Nf × 1 column vector (coherence frequency axis)
fig2_out.freq = trial_data.coh_freq(:);

% --- Spectrogram ---
% Power in dB/Hz (Nf × Nt matrix)
% Note: spec_power is already in linear units, convert to dB
spec_db = 10 * log10(trial_data.spec_power + eps);

% Save spectrogram on NATIVE time axis (not interpolated to coherence axis)
% This preserves full resolution regardless of coherence method
fig2_out.spec_power = spec_db;                    % Nf_spec × Nt_spec, dB/Hz
fig2_out.spec_time = trial_data.spec_time(:)';    % 1 × Nt_spec
fig2_out.spec_freq = trial_data.spec_freq(:);     % Nf_spec × 1

% --- Coherence (method-specific field name) ---
% Magnitude-squared coherence (Nf × Nt matrix, 0-1)
if strcmpi(method, 'mscohere')
    fig2_out.coh_mscohere = trial_data.coh_time_resolved;  % Nf × Nt
else
    fig2_out.coh_fieldtrip = trial_data.coh_time_resolved;  % Nf × Nt
end

% Overall coherence spectrum (Nf × 1)
fig2_out.coh_spectrum = trial_data.coh_mag(:);

% --- Rest vs Run Coherence ---
% Behavioral state-dependent coherence spectra
fig2_out.pct_rest = trial_data.behavior.pct_rest;
fig2_out.pct_run = trial_data.behavior.pct_run;
fig2_out.run_threshold_cms = cfg.behavior.run_threshold;
fig2_out.min_bout_sec = cfg.behavior.min_bout_duration_sec;

if ~isempty(trial_data.coh_rest)
    fig2_out.coh_rest = trial_data.coh_rest(:);
    fig2_out.coh_rest_freq = trial_data.coh_rest_freq(:);
else
    fig2_out.coh_rest = [];
    fig2_out.coh_rest_freq = [];
end

if ~isempty(trial_data.coh_run)
    fig2_out.coh_run = trial_data.coh_run(:);
    fig2_out.coh_run_freq = trial_data.coh_run_freq(:);
else
    fig2_out.coh_run = [];
    fig2_out.coh_run_freq = [];
end

% --- Motion (if available) ---
% Running speed in cm/s, interpolated to SPECTROGRAM time axis (not coherence)
% This ensures consistent motion display regardless of coherence method
if isfield(trial_data, 'speed') && ~isempty(trial_data.speed)
    % Interpolate to spectrogram time axis for method-consistent display
    motion_interp = interp1(trial_data.t, trial_data.speed, fig2_out.spec_time, 'linear', NaN);
    fig2_out.motion = motion_interp(:)';
    fig2_out.motion_time = fig2_out.spec_time;  % Motion uses spectrogram time axis
else
    fig2_out.motion = [];
    fig2_out.motion_time = [];
end

% --- Power Spectral Density (PSD) ---
% Overall PSD for LFP and GEVI (dB units)
fig2_out.psd_freq = trial_data.psd_freq(:);           % Frequency axis (Nf × 1)
fig2_out.psd_lfp = trial_data.psd_lfp(:);             % LFP PSD (dB re 1 µV²/Hz)
fig2_out.psd_gevi = trial_data.psd_gevi(:);           % GEVI PSD (dB re 1 (ΔF/F)²/Hz)

% Rest vs Run PSD
if ~isempty(trial_data.psd_lfp_rest)
    fig2_out.psd_lfp_rest = trial_data.psd_lfp_rest(:);
    fig2_out.psd_gevi_rest = trial_data.psd_gevi_rest(:);
    fig2_out.psd_rest_freq = trial_data.psd_rest_freq(:);
else
    fig2_out.psd_lfp_rest = [];
    fig2_out.psd_gevi_rest = [];
    fig2_out.psd_rest_freq = [];
end

if ~isempty(trial_data.psd_lfp_run)
    fig2_out.psd_lfp_run = trial_data.psd_lfp_run(:);
    fig2_out.psd_gevi_run = trial_data.psd_gevi_run(:);
    fig2_out.psd_run_freq = trial_data.psd_run_freq(:);
else
    fig2_out.psd_lfp_run = [];
    fig2_out.psd_gevi_run = [];
    fig2_out.psd_run_freq = [];
end

% --- Metadata (simple scalars/strings only) ---
fig2_out.fs = trial_data.fs;
fig2_out.duration = trial_data.duration;
fig2_out.method = method;
fig2_out.trial_label = trial_data.trial_label;
fig2_out.analysis_date = datestr(now, 'yyyy-mm-dd');

% --- Units documentation ---
fig2_out.units_spec = 'dB/Hz (re 1 uV^2/Hz)';
fig2_out.units_coh = 'magnitude-squared coherence (0-1)';
fig2_out.units_psd_lfp = 'dB (re 1 uV^2/Hz)';
fig2_out.units_psd_gevi = 'dB (re 1 (dF/F)^2/Hz)';
fig2_out.units_time = 'seconds';
fig2_out.units_freq = 'Hz';
fig2_out.units_motion = 'cm/s';

% Save as struct (Python-compatible format)
% Legacy naming: figure2_{method}_trial{idx}.mat
filename = sprintf('figure2_%s_trial%d.mat', method, trial_idx);
filepath = fullfile(output_dir, filename);

save(filepath, '-struct', 'fig2_out', '-v7');  % -v7 for scipy compatibility

end

%% ============================================================================
%  HELPER: Safe Mean Function (handles empty arrays)
%  ============================================================================
function result = safe_mean(data)
%SAFE_MEAN Compute mean along columns, returning [] if input is empty
if isempty(data)
    result = [];
else
    result = mean(data, 2);
end
end

%% ============================================================================
%  HELPER: 2D Smoothing Function
%  ============================================================================
function smoothed = smooth2a(data, row_window, col_window)
%SMOOTH2A 2D smoothing using moving average
%  Applies smoothing in both dimensions (rows=freq, cols=time)

if nargin < 3
    col_window = row_window;
end

[rows, cols] = size(data);
smoothed = data;

% Smooth along rows (frequency)
if row_window > 1
    kernel_row = ones(row_window, 1) / row_window;
    for c = 1:cols
        smoothed(:, c) = conv(data(:, c), kernel_row, 'same');
    end
end

% Smooth along columns (time)
if col_window > 1
    kernel_col = ones(1, col_window) / col_window;
    for r = 1:rows
        smoothed(r, :) = conv(smoothed(r, :), kernel_col, 'same');
    end
end

end

%% ============================================================================
%  HELPER: Find Continuous Bouts
%  ============================================================================
function bouts = find_continuous_bouts(indices, min_length)
%FIND_CONTINUOUS_BOUTS Identify continuous bouts from sparse indices

if isempty(indices)
    bouts = [];
    return;
end

indices = indices(:);
diff_idx = diff(indices);
break_points = find(diff_idx > 1);

bout_starts = [indices(1); indices(break_points + 1)];
bout_ends = [indices(break_points); indices(end)];

bout_lengths = bout_ends - bout_starts + 1;
valid = bout_lengths >= min_length;

bouts = [bout_starts(valid), bout_ends(valid)];

end

%% ============================================================================
%  HELPER: Save Results
%  ============================================================================
function save_results(results, animal, session_id, level, method, cfg)
%SAVE_RESULTS Save results to MAT file in Python-compatible format
%
%  Output structure:
%    {output_root}/{behavior_mode}/{level}/{MouseID}/{SessionID}/data/{method}.mat
%
%  For animal_pooled or animal_concatenated (no session):
%    {output_root}/{behavior_mode}/{level}/{MouseID}/data/{method}.mat
%
%  Fields are flattened for Python compatibility (scipy.io.loadmat)

if ~cfg.output.save_mat
    return;
end

% Get behavior mode for folder structure
% Note: output_folder_suffix is already applied to cfg.paths.output_root in pipeline
behavior_mode = cfg.behavior.classification_mode;

% Build output path based on level
% animal_pooled and animal_concatenated don't have session subdirectories
if strcmp(level, 'animal_pooled') || strcmp(level, 'animal_concatenated')
    output_dir = fullfile(cfg.paths.output_root, behavior_mode, level, animal.mouse_id, 'data');
else
    output_dir = fullfile(cfg.paths.output_root, behavior_mode, level, animal.mouse_id, session_id, 'data');
end

filename = sprintf('%s.mat', method);

% Create directory if needed
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

filepath = fullfile(output_dir, filename);

% Create Python-compatible flat structure
out = struct();

% For pooled data, flatten the spectra structure
if isfield(results, 'spectra') && ~isempty(results.spectra)
    sp = results.spectra;
    
    % Coherence
    if isfield(sp, 'coherence')
        out.freq = sp.coherence.freq(:);
        out.coh_overall = sp.coherence.overall(:);
        if isfield(sp.coherence, 'rest') && ~isempty(sp.coherence.rest)
            out.coh_rest = sp.coherence.rest(:);
        else
            out.coh_rest = [];
        end
        if isfield(sp.coherence, 'run') && ~isempty(sp.coherence.run)
            out.coh_run = sp.coherence.run(:);
        else
            out.coh_run = [];
        end
    end
    
    % PSD frequency axis
    if isfield(sp, 'psd_freq')
        out.psd_freq = sp.psd_freq(:);
    end
    
    % LFP PSD
    if isfield(sp, 'psd_lfp')
        out.psd_lfp = sp.psd_lfp.overall(:);
        if isfield(sp.psd_lfp, 'rest') && ~isempty(sp.psd_lfp.rest)
            out.psd_lfp_rest = sp.psd_lfp.rest(:);
        else
            out.psd_lfp_rest = [];
        end
        if isfield(sp.psd_lfp, 'run') && ~isempty(sp.psd_lfp.run)
            out.psd_lfp_run = sp.psd_lfp.run(:);
        else
            out.psd_lfp_run = [];
        end
    end
    
    % GEVI PSD
    if isfield(sp, 'psd_gevi')
        out.psd_gevi = sp.psd_gevi.overall(:);
        if isfield(sp.psd_gevi, 'rest') && ~isempty(sp.psd_gevi.rest)
            out.psd_gevi_rest = sp.psd_gevi.rest(:);
        else
            out.psd_gevi_rest = [];
        end
        if isfield(sp.psd_gevi, 'run') && ~isempty(sp.psd_gevi.run)
            out.psd_gevi_run = sp.psd_gevi.run(:);
        else
            out.psd_gevi_run = [];
        end
    end
end

% Behavior info
if isfield(results, 'behavior')
    out.pct_rest = results.behavior.pct_rest;
    out.pct_run = results.behavior.pct_run;
    out.run_threshold_cms = cfg.behavior.run_threshold;
    out.min_bout_sec = cfg.behavior.min_bout_duration_sec;
end

% Number of sessions (for animal-pooled level)
if isfield(results, 'num_sessions')
    out.num_sessions = results.num_sessions;
end

% Metadata
out.mouse_id = animal.mouse_id;
out.method = method;
out.level = level;
out.behavior_mode = behavior_mode;
out.analysis_date = datestr(now, 'yyyy-mm-dd');

if ~strcmp(level, 'animal_pooled')
    out.session_id = session_id;
end

% Save as flat struct for Python compatibility
save(filepath, '-struct', 'out', '-v7');  % -v7 for scipy compatibility

fprintf('        Saved: %s\n', filepath);

end

%% ============================================================================
%  ARTIFACT EXCLUSION HELPER
%  ============================================================================
function [should_exclude, artifact_pct] = check_artifact_exclusion(trial_path, trial_idx, cfg)
%CHECK_ARTIFACT_EXCLUSION Check if a trial should be excluded based on artifact info
%
%  Looks for artifact info file in the same directory as the trial file
%  and checks if the trial is marked for exclusion.
%
%  INPUTS:
%    trial_path - Full path to the trial data file
%    trial_idx  - Trial index within the session
%    cfg        - Configuration struct with artifact_exclusion settings
%
%  OUTPUTS:
%    should_exclude - true if trial should be excluded
%    artifact_pct   - percentage of trial that is artifacts (for logging)

should_exclude = false;
artifact_pct = 0;

% Get the directory containing the trial file
trial_dir = fileparts(trial_path);

% Look for artifact info file - could be in trial folder or session folder
% Pattern: {MouseName}-{SessionID}_artifact_removal.mat or similar

% First try: look in the session folder (parent of trial folder)
session_dir = fileparts(trial_dir);
artifact_files = dir(fullfile(session_dir, '*_artifact_removal.mat'));

% Second try: look directly in trial folder
if isempty(artifact_files)
    artifact_files = dir(fullfile(trial_dir, '*_artifact_removal.mat'));
end

if isempty(artifact_files)
    % No artifact info file found - include trial by default
    return;
end

% Load the artifact info
try
    artifact_path = fullfile(artifact_files(1).folder, artifact_files(1).name);
    loaded = load(artifact_path, 'ArtifactInfo');
    
    if ~isfield(loaded, 'ArtifactInfo')
        return;
    end
    
    info = loaded.ArtifactInfo;
    
    % Find this trial in the artifact info
    if isfield(info, 'trials') && length(info.trials) >= trial_idx
        trial_info = info.trials(trial_idx);
        
        % Get artifact percentage
        if isfield(trial_info, 'artifact_pct')
            artifact_pct = trial_info.artifact_pct;
        end
        
        % Get threshold from config (prefer new field, fall back to legacy)
        % Threshold is stored as fraction (0.08 = 8%), but artifact_pct is percentage (0-100)
        threshold = 0.30;  % Default fallback (30%)
        if isfield(cfg, 'artifact') && isfield(cfg.artifact, 'threshold')
            threshold = cfg.artifact.threshold;
        elseif isfield(cfg, 'artifact_exclusion') && isfield(cfg.artifact_exclusion, 'threshold')
            threshold = cfg.artifact_exclusion.threshold;
        end
        
        % CRITICAL: Compare artifact_pct against the threshold from config
        % The artifact removal tool's final_decision uses its own threshold (30% default),
        % but we need to use the threshold specified in the spectral analysis pipeline.
        % artifact_pct is stored as percentage (0-100), threshold is fraction (0-1)
        % Convert threshold to percentage for comparison
        threshold_pct = threshold * 100;
        if artifact_pct > threshold_pct
            should_exclude = true;
        end
        
        % Note: We ignore the pre-made final_decision from artifact removal tool
        % because it was made with a different threshold. We make our own decision
        % based on the pipeline's threshold.
    end
    
catch ME
    warning('artifact_exclusion:loadFailed', ...
        'Failed to load artifact info: %s', ME.message);
end

end

%% ============================================================================
%  ARTIFACT CLEANING FUNCTIONS
%  ============================================================================
%  These functions implement the 'clean' artifact mode, which removes artifact
%  time segments from within trials (as opposed to 'exclude' mode which skips
%  entire trials).
%
%  TERMINOLOGY:
%    - 'exclude' mode: Skip entire trials with >threshold% contamination
%    - 'clean' mode:   Remove artifact segments, keep clean portions of ALL trials
%  ============================================================================

function [cleaned_data, cleaning_info] = apply_artifact_cleaning(data, trial_path, trial_idx, cfg)
%APPLY_ARTIFACT_CLEANING Remove artifact segments from trial data
%
%  Loads artifact info from the corresponding *_artifact_removal.mat file
%  and removes flagged time segments from LFP, GEVI, speed, and time data.
%
%  INPUTS:
%    data       - Struct with fields: time, fs, lfp, gevi, speed
%    trial_path - Full path to the trial data file
%    trial_idx  - Trial index within the session
%    cfg        - Configuration struct
%
%  OUTPUTS:
%    cleaned_data  - Data struct with artifact segments removed
%    cleaning_info - Struct with cleaning statistics:
%                    .original_samples   - Original number of samples
%                    .cleaned_samples    - Samples after cleaning
%                    .removed_samples    - Number of samples removed
%                    .removed_pct        - Percentage of data removed
%                    .num_segments       - Number of artifact segments
%                    .artifact_mask      - Original artifact mask (for reference)

cleaning_info = struct();
cleaning_info.original_samples = length(data.time);
cleaning_info.cleaned_samples = cleaning_info.original_samples;
cleaning_info.removed_samples = 0;
cleaning_info.removed_pct = 0;
cleaning_info.removed_pct_original = 0;  % Before intelligent processing
cleaning_info.num_segments = 0;
cleaning_info.num_segments_original = 0;
cleaning_info.num_segments_processed = 0;
cleaning_info.artifact_mask_original = [];
cleaning_info.artifact_mask_processed = [];
cleaning_info.applied = false;

% Default: return data unchanged
cleaned_data = data;

% Load artifact mask for this trial
[artifact_mask, mask_info] = load_artifact_mask(trial_path, trial_idx);

if isempty(artifact_mask)
    % No artifact info found - return data unchanged
    warning('artifact_cleaning:noMask', ...
        'No artifact mask found for trial %d (path: %s). Data returned unchanged.', ...
        trial_idx, trial_path);
    return;
end

% Diagnostic: Check if mask has any artifacts
n_artifacts_original = sum(artifact_mask);
if n_artifacts_original == 0
    warning('artifact_cleaning:noArtifacts', ...
        'Artifact mask for trial %d contains no artifacts (all false). Data returned unchanged.', ...
        trial_idx);
    % Still process to update cleaning_info, but no actual cleaning needed
    cleaning_info.artifact_mask_original = artifact_mask;
    cleaning_info.num_segments_original = 0;
    cleaning_info.removed_pct_original = 0;
    cleaning_info.applied = false;
    return;
end

% Store original mask info
cleaning_info.artifact_mask_original = artifact_mask;
cleaning_info.num_segments_original = mask_info.num_segments;

% Ensure mask length matches data length
n_data = length(data.time);
n_mask = length(artifact_mask);

if n_mask ~= n_data
    % Try to handle length mismatch
    if n_mask > n_data
        artifact_mask = artifact_mask(1:n_data);
        warning('artifact_cleaning:lengthMismatch', ...
            'Artifact mask longer than data (%d vs %d). Truncating mask.', n_mask, n_data);
    else
        % Pad mask with false (assume extra samples are clean)
        artifact_mask = [artifact_mask; false(n_data - n_mask, 1)];
        warning('artifact_cleaning:lengthMismatch', ...
            'Artifact mask shorter than data (%d vs %d). Padding with clean.', n_mask, n_data);
    end
end

% ============================================================================
% INTELLIGENT ARTIFACT MASK PROCESSING
% ============================================================================
% Apply intelligent processing to the artifact mask:
%   1. Optional smoothing to merge nearby artifacts
%   2. Ensure contiguous blocks are merged
%   3. Add padding around artifacts to remove potentially contaminated data
% ============================================================================

% Get cleaning parameters from config (with defaults)
pre_pad_sec = 0.100;   % Default: 100ms
post_pad_sec = 0.100;  % Default: 100ms
smooth_window_sec = 0.050;  % Default: 50ms

if isfield(cfg, 'artifact')
    if isfield(cfg.artifact, 'pre_pad_sec')
        pre_pad_sec = cfg.artifact.pre_pad_sec;
    end
    if isfield(cfg.artifact, 'post_pad_sec')
        post_pad_sec = cfg.artifact.post_pad_sec;
    end
    if isfield(cfg.artifact, 'smooth_window_sec')
        smooth_window_sec = cfg.artifact.smooth_window_sec;
    end
end

fs = data.fs;

% Step 1: Optional smoothing to merge nearby artifacts
% This helps merge artifacts that are very close together (within smooth_window)
% Uses morphological dilation: if any sample in a window around a point is artifact,
% mark that point as artifact. This effectively extends each artifact by half the
% window size on each side, merging nearby artifacts.
if smooth_window_sec > 0
    smooth_samples = max(1, round(smooth_window_sec * fs));
    if smooth_samples > 1
        % Efficient dilation using convolution with a ones kernel
        % This is equivalent to: for each point, check if any sample in window is artifact
        half_win = floor(smooth_samples / 2);
        artifact_mask_dilated = false(size(artifact_mask));
        
        % Use a more efficient approach: find all artifact positions and dilate them
        artifact_indices = find(artifact_mask);
        for idx = artifact_indices(:)'
            start_idx = max(1, idx - half_win);
            end_idx = min(length(artifact_mask), idx + half_win);
            artifact_mask_dilated(start_idx:end_idx) = true;
        end
        
        artifact_mask = artifact_mask_dilated;
    end
end

% Step 2: Ensure contiguous blocks are merged (should already be done, but ensure)
% Find all artifact segments and merge any that are adjacent
artifact_segments = find_artifact_segments(artifact_mask);
if ~isempty(artifact_segments)
    % Reconstruct mask from segments (this ensures contiguous blocks)
    artifact_mask_merged = false(size(artifact_mask));
    for s = 1:size(artifact_segments, 1)
        start_idx = artifact_segments(s, 1);
        end_idx = artifact_segments(s, 2);
        artifact_mask_merged(start_idx:end_idx) = true;
    end
    artifact_mask = artifact_mask_merged;
end

% Step 3: Add padding around artifacts
% Expand each artifact segment by pre_pad before and post_pad after
pre_pad_samples = round(pre_pad_sec * fs);
post_pad_samples = round(post_pad_sec * fs);

if pre_pad_samples > 0 || post_pad_samples > 0
    artifact_mask_expanded = artifact_mask;
    
    % Find artifact segments
    segments = find_artifact_segments(artifact_mask);
    
    for s = 1:size(segments, 1)
        seg_start = segments(s, 1);
        seg_end = segments(s, 2);
        
        % Expand segment with padding
        expanded_start = max(1, seg_start - pre_pad_samples);
        expanded_end = min(length(artifact_mask), seg_end + post_pad_samples);
        
        % Mark expanded region as artifact
        artifact_mask_expanded(expanded_start:expanded_end) = true;
    end
    
    artifact_mask = artifact_mask_expanded;
end

% Store processed mask
cleaning_info.artifact_mask_processed = artifact_mask;
cleaning_info.num_segments_processed = size(find_artifact_segments(artifact_mask), 1);

% Create clean mask (inverse of processed artifact mask)
clean_mask = ~artifact_mask;

% Count statistics (using processed mask)
n_clean = sum(clean_mask);
n_artifact = sum(artifact_mask);
n_original_artifact = sum(cleaning_info.artifact_mask_original);

cleaning_info.cleaned_samples = n_clean;
cleaning_info.removed_samples = n_artifact;
cleaning_info.removed_pct = 100 * n_artifact / n_data;
cleaning_info.removed_pct_original = 100 * n_original_artifact / n_data;  % Before processing
cleaning_info.num_segments = cleaning_info.num_segments_processed;
cleaning_info.applied = true;

% Diagnostic warning if no artifacts remain after processing (shouldn't happen, but check)
if n_artifact == 0 && n_original_artifact > 0
    warning('artifact_cleaning:processingRemovedAll', ...
        'WARNING: All artifacts were removed during processing for trial %d! Original: %d artifacts, Processed: 0. This may indicate a bug.', ...
        trial_idx, n_original_artifact);
end

% Diagnostic warning if very little data remains
if cleaning_info.removed_pct > 50
    warning('artifact_cleaning:excessiveRemoval', ...
        'WARNING: Trial %d has %.1f%% of data removed by artifact cleaning. This may affect analysis quality.', ...
        trial_idx, cleaning_info.removed_pct);
end

% Apply cleaning - HARD REMOVAL (samples dropped entirely, not NaN-masked)
% IMPORTANT: This creates discontinuities in time, but spectral analysis
% methods (PSD via pwelch, coherence via mscohere) work on the concatenated
% clean segments. These functions handle gaps by processing segments independently
% and then averaging, so hard removal is appropriate and more efficient than NaN masking.
%
% The cleaned data will have:
%   - Shorter time vector (only clean samples)
%   - Discontinuous time (gaps where artifacts were removed)
%   - All signals (LFP, GEVI, speed) aligned to the same indices
cleaned_data.time = data.time(clean_mask);
cleaned_data.lfp = data.lfp(clean_mask);
cleaned_data.gevi = data.gevi(clean_mask);
cleaned_data.speed = data.speed(clean_mask);
cleaned_data.fs = data.fs;  % Sampling rate unchanged
cleaned_data.duration = length(cleaned_data.time) / data.fs;  % Updated duration (actual clean data duration)

end

function [artifact_mask, mask_info] = load_artifact_mask(trial_path, trial_idx)
%LOAD_ARTIFACT_MASK Load artifact mask for a specific trial
%
%  Searches for *_artifact_removal.mat file in session folder and extracts
%  the artifact mask for the specified trial.
%
%  INPUTS:
%    trial_path - Full path to the trial data file
%    trial_idx  - Trial index within the session
%
%  OUTPUTS:
%    artifact_mask - Boolean vector (true = artifact, false = clean)
%    mask_info     - Struct with: .num_segments, .artifact_pct

artifact_mask = [];
mask_info = struct('num_segments', 0, 'artifact_pct', 0);

% Get the directory containing the trial file
trial_dir = fileparts(trial_path);

% Look for artifact info file
% Pattern: {MouseName}-{SessionID}_artifact_removal.mat
% Can be in session folder (parent of trial folder) or trial folder

% First try: session folder (parent of trial folder)
session_dir = fileparts(trial_dir);
artifact_files = dir(fullfile(session_dir, '*_artifact_removal.mat'));

% Second try: directly in trial folder (for alternate structures)
if isempty(artifact_files)
    artifact_files = dir(fullfile(trial_dir, '*_artifact_removal.mat'));
end

% Third try: in the same folder as trial (for combined sessions)
if isempty(artifact_files)
    % Try parent's parent for deeply nested structures
    parent_dir = fileparts(session_dir);
    artifact_files = dir(fullfile(parent_dir, '*_artifact_removal.mat'));
end

if isempty(artifact_files)
    % No artifact info file found - this is expected for some sessions
    % Don't warn here as it's normal for sessions without artifact removal
    return;
end

% Load the artifact info
try
    artifact_path = fullfile(artifact_files(1).folder, artifact_files(1).name);
    loaded = load(artifact_path, 'ArtifactInfo');
    
    if ~isfield(loaded, 'ArtifactInfo')
        warning('artifact_cleaning:invalidFile', ...
            'ArtifactInfo struct not found in %s', artifact_path);
        return;
    end
    
    info = loaded.ArtifactInfo;
    
    % Find this trial in the artifact info
    if ~isfield(info, 'trials')
        warning('artifact_cleaning:noTrialsField', ...
            'ArtifactInfo has no "trials" field in %s', artifact_path);
        return;
    end
    
    if length(info.trials) < trial_idx
        warning('artifact_cleaning:trialNotFound', ...
            'Trial %d not found in artifact info (has %d trials, file: %s). Check trial numbering.', ...
            trial_idx, length(info.trials), artifact_path);
        return;
    end
    
    trial_info = info.trials(trial_idx);
    
    % Extract artifact mask
    if isfield(trial_info, 'artifact_mask')
        artifact_mask = trial_info.artifact_mask(:);  % Ensure column vector
        
        % Check if mask is valid
        if isempty(artifact_mask)
            warning('artifact_cleaning:emptyMask', ...
                'Artifact mask for trial %d is empty in %s', trial_idx, artifact_path);
            return;
        end
        
        if ~islogical(artifact_mask)
            % Convert to logical if needed
            artifact_mask = logical(artifact_mask);
        end
        
        mask_info.artifact_pct = trial_info.artifact_pct;
        
        if isfield(trial_info, 'num_artifacts')
            mask_info.num_segments = trial_info.num_artifacts;
        elseif isfield(trial_info, 'artifact_segments')
            mask_info.num_segments = size(trial_info.artifact_segments, 1);
        else
            % Count segments manually if not stored
            segments = find_artifact_segments(artifact_mask);
            mask_info.num_segments = size(segments, 1);
        end
        
        % Diagnostic: Log mask info
        n_artifacts = sum(artifact_mask);
        if n_artifacts > 0
            fprintf('        [Artifact mask loaded: %d artifacts (%.1f%%), %d segments]\n', ...
                n_artifacts, mask_info.artifact_pct, mask_info.num_segments);
        end
    else
        warning('artifact_cleaning:noMask', ...
            'No artifact_mask field in trial %d of %s', trial_idx, artifact_path);
    end
    
catch ME
    warning('artifact_cleaning:loadFailed', ...
        'Failed to load artifact mask: %s', ME.message);
end

end

%% ============================================================================
%  HELPER: Find Artifact Segments
%  ============================================================================
function segments = find_artifact_segments(artifact_mask)
%FIND_ARTIFACT_SEGMENTS Find contiguous artifact segments in a binary mask
%
%  Identifies all contiguous blocks of artifact samples (true values) in the mask.
%  Returns a matrix where each row is [start_idx, end_idx] for one segment.
%
%  INPUT:
%    artifact_mask - Logical vector (true = artifact, false = clean)
%
%  OUTPUT:
%    segments - [N × 2] matrix where each row is [start_idx, end_idx]
%               Empty matrix if no artifacts found

segments = [];

if isempty(artifact_mask) || ~any(artifact_mask)
    return;
end

% Find transitions: where mask changes from false to true (start) or true to false (end)
mask_diff = diff([false; artifact_mask(:); false]);
starts = find(mask_diff == 1);   % Transitions from false to true
ends = find(mask_diff == -1) - 1; % Transitions from true to false (adjust for diff offset)

% Combine into segments matrix
if ~isempty(starts) && length(starts) == length(ends)
    segments = [starts, ends];
end

end

%% ============================================================================
%  HELPER: Extract Trial Number from File Path
%  ============================================================================
function trial_num = extract_trial_number_from_path(trial_path)
%EXTRACT_TRIAL_NUMBER_FROM_PATH Extract trial number from filename
%
%  Parses trial number from filename patterns like:
%    - "..._Trial1_..." → 1
%    - "..._Trial2_..." → 2
%    - "MouseName-SessionID_Trial3_FiberPhotometry_Analysis.mat" → 3
%
%  This is CRITICAL for artifact cleaning in combined sessions where
%  the loop index doesn't match the actual trial numbering in artifact info.
%
%  INPUT:
%    trial_path - Full path to trial file
%
%  OUTPUT:
%    trial_num - Extracted trial number (integer), or [] if not found

trial_num = [];

% Get just the filename
[~, filename, ~] = fileparts(trial_path);

% Try to match pattern "_Trial{N}_" in filename
% Common patterns:
%   MouseName-SessionID_Trial1_FiberPhotometry_Analysis
%   Trial1_fov1_baselineRecording_60sec_1
tokens = regexp(filename, '_Trial(\d+)_', 'tokens');

if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Try pattern "Trial{N}_" at start of filename (for trial folder names)
tokens = regexp(filename, '^Trial(\d+)_', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Try to extract from folder path (e.g., ".../Trial1_xxx/...")
tokens = regexp(trial_path, '[/\\]Trial(\d+)_', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
    return;
end

% Fallback: try any pattern with Trial followed by a number
tokens = regexp(filename, 'Trial(\d+)', 'tokens');
if ~isempty(tokens)
    trial_num = str2double(tokens{1}{1});
end

end
