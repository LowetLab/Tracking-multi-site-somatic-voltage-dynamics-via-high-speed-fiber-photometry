function [is_rest, is_run, behavior_info] = classify_behavior(speed, fs, cfg)
%% ============================================================================
%  CLASSIFY_BEHAVIOR - Classify behavioral states (REST vs RUN) from speed
%  ============================================================================
%
%  SYNTAX:
%    [is_rest, is_run, behavior_info] = classify_behavior(speed, fs, cfg)
%
%  INPUTS:
%    speed   - Speed time series (cm/s), [1 × N] or [N × 1]
%    fs      - Sampling frequency (Hz)
%    cfg     - Configuration struct from analysis_config() containing:
%              cfg.behavior.classification_mode: 'standard' or 'clear'
%              cfg.behavior.run_threshold: Speed threshold for RUN (cm/s)
%              cfg.behavior.rest_threshold: Speed threshold for REST (cm/s)
%              cfg.behavior.min_bout_duration_sec: Minimum bout duration (s)
%              cfg.behavior.apply_bout_filter: Apply bout duration filtering
%              cfg.behavior.motion_smooth_samples: Samples for smoothing
%
%  OUTPUTS:
%    is_rest       - Logical array, true where state = REST [N × 1]
%    is_run        - Logical array, true where state = RUN [N × 1]
%    behavior_info - Struct with classification details:
%                    .pct_rest, .pct_run, .pct_excluded
%                    .n_rest_samples, .n_run_samples, .n_excluded_samples
%                    .mode, .run_threshold, .rest_threshold
%                    .min_bout_sec, .n_rest_bouts, .n_run_bouts
%
%  CLASSIFICATION MODES:
%
%    'standard' (Original method):
%      - RUN:  speed > run_threshold
%      - REST: speed <= run_threshold
%      - Short bouts merged into surrounding state
%      - REST + RUN = 100%
%
%    'clear' (Strict classification):
%      - RUN:  speed > run_threshold
%      - REST: speed < rest_threshold
%      - EXCLUDED: rest_threshold <= speed <= run_threshold
%      - REST + RUN + EXCLUDED = 100%
%      - Short bouts filtered (not merged, just excluded)
%  ============================================================================

%% Input validation
if nargin < 3
    error('classify_behavior requires 3 inputs: speed, fs, cfg');
end

speed = speed(:);  % Ensure column vector
n_samples = length(speed);

% Extract configuration
mode = cfg.behavior.classification_mode;
run_threshold = cfg.behavior.run_threshold;
rest_threshold = cfg.behavior.rest_threshold;
apply_bout_filter = cfg.behavior.apply_bout_filter;
smooth_samples = cfg.behavior.motion_smooth_samples;

% Get bout duration parameters - support both legacy and new separate parameters
min_bout_sec = cfg.behavior.min_bout_duration_sec;  % Legacy: shared for standard mode

% For 'clear' mode, use separate min bout durations if available
if isfield(cfg.behavior, 'min_run_bout_sec')
    min_run_bout_sec = cfg.behavior.min_run_bout_sec;
else
    min_run_bout_sec = min_bout_sec;  % Fallback to legacy
end

if isfield(cfg.behavior, 'min_rest_bout_sec')
    min_rest_bout_sec = cfg.behavior.min_rest_bout_sec;
else
    min_rest_bout_sec = min_bout_sec;  % Fallback to legacy
end

min_bout_samples = round(min_bout_sec * fs);
min_run_bout_samples = round(min_run_bout_sec * fs);
min_rest_bout_samples = round(min_rest_bout_sec * fs);

%% Smooth speed signal (optional)
if smooth_samples > 1
    speed_smooth = movmean(speed, smooth_samples);
else
    speed_smooth = speed;
end

%% Classification based on mode
switch lower(mode)
    case 'standard'
        % Standard mode: binary classification with bout merging
        % Uses shared min_bout_samples for both states
        [is_rest, is_run, behavior_info] = classify_standard(speed_smooth, ...
            run_threshold, min_bout_samples, apply_bout_filter, fs);
        
    case 'clear'
        % Clear mode: three-state classification (REST, RUN, EXCLUDED)
        % Uses SEPARATE min bout durations for REST and RUN
        [is_rest, is_run, behavior_info] = classify_clear(speed_smooth, ...
            run_threshold, rest_threshold, min_run_bout_samples, min_rest_bout_samples, ...
            apply_bout_filter, fs);
        
    otherwise
        error('Unknown classification mode: %s. Use ''standard'' or ''clear''.', mode);
end

%% Add configuration info to output
behavior_info.mode = mode;
behavior_info.run_threshold = run_threshold;
behavior_info.rest_threshold = rest_threshold;
behavior_info.min_bout_sec = min_bout_sec;  % Legacy parameter
behavior_info.min_run_bout_sec = min_run_bout_sec;
behavior_info.min_rest_bout_sec = min_rest_bout_sec;
behavior_info.apply_bout_filter = apply_bout_filter;
behavior_info.fs = fs;
behavior_info.total_samples = n_samples;
behavior_info.total_duration_sec = n_samples / fs;

end

%% ============================================================================
%  HELPER FUNCTION: Standard Classification
%  ============================================================================
function [is_rest, is_run, info] = classify_standard(speed, run_threshold, ...
    min_bout_samples, apply_bout_filter, fs)
%CLASSIFY_STANDARD Binary classification with bout merging
%  Everything is either REST or RUN. Short bouts are merged into surrounding.

n_samples = length(speed);

% Initial classification
is_running_raw = speed > run_threshold;

if apply_bout_filter && min_bout_samples > 1
    % Merge short bouts into surrounding state
    is_running_filtered = merge_short_bouts(is_running_raw, min_bout_samples);
else
    is_running_filtered = is_running_raw;
end

% Final assignment
is_run = is_running_filtered;
is_rest = ~is_running_filtered;

% Compute statistics
n_rest = sum(is_rest);
n_run = sum(is_run);

info.pct_rest = 100 * n_rest / n_samples;
info.pct_run = 100 * n_run / n_samples;
info.pct_excluded = 0;  % No exclusion in standard mode

info.n_rest_samples = n_rest;
info.n_run_samples = n_run;
info.n_excluded_samples = 0;

info.duration_rest_sec = n_rest / fs;
info.duration_run_sec = n_run / fs;
info.duration_excluded_sec = 0;

% Count bouts
[info.n_rest_bouts, info.n_run_bouts] = count_bouts(is_rest, is_run);

end

%% ============================================================================
%  HELPER FUNCTION: Clear Classification
%  ============================================================================
function [is_rest, is_run, info] = classify_clear(speed, run_threshold, ...
    rest_threshold, min_run_bout_samples, min_rest_bout_samples, apply_bout_filter, fs)
%CLASSIFY_CLEAR Three-state classification with exclusion zone
%  REST: speed < rest_threshold AND bout duration >= min_rest_bout_samples
%  RUN:  speed > run_threshold AND bout duration >= min_run_bout_samples
%  EXCLUDED: intermediate speeds OR short bouts
%
%  This mode uses SEPARATE minimum bout durations for REST and RUN,
%  allowing more flexible behavioral classification.

n_samples = length(speed);

% Initial three-state classification based on speed thresholds only
is_running_raw = speed > run_threshold;
is_resting_raw = speed < rest_threshold;
is_ambiguous = ~is_running_raw & ~is_resting_raw;  % Intermediate zone

if apply_bout_filter
    % For clear mode: filter short bouts by excluding them (not merging)
    % Use SEPARATE min bout durations for REST and RUN
    if min_run_bout_samples > 1
        is_running_filtered = filter_short_bouts_exclude(is_running_raw, min_run_bout_samples);
    else
        is_running_filtered = is_running_raw;
    end
    
    if min_rest_bout_samples > 1
        is_resting_filtered = filter_short_bouts_exclude(is_resting_raw, min_rest_bout_samples);
    else
        is_resting_filtered = is_resting_raw;
    end
else
    is_running_filtered = is_running_raw;
    is_resting_filtered = is_resting_raw;
end

% Final assignment
is_run = is_running_filtered;
is_rest = is_resting_filtered;
is_excluded = ~is_run & ~is_rest;

% Compute statistics
n_rest = sum(is_rest);
n_run = sum(is_run);
n_excluded = sum(is_excluded);

info.pct_rest = 100 * n_rest / n_samples;
info.pct_run = 100 * n_run / n_samples;
info.pct_excluded = 100 * n_excluded / n_samples;

info.n_rest_samples = n_rest;
info.n_run_samples = n_run;
info.n_excluded_samples = n_excluded;

info.duration_rest_sec = n_rest / fs;
info.duration_run_sec = n_run / fs;
info.duration_excluded_sec = n_excluded / fs;

% Count bouts
[info.n_rest_bouts, info.n_run_bouts] = count_bouts(is_rest, is_run);

end

%% ============================================================================
%  HELPER FUNCTION: Merge Short Bouts (for Standard Mode)
%  ============================================================================
function is_running_filtered = merge_short_bouts(is_running, min_samples)
%MERGE_SHORT_BOUTS Merge short bouts into surrounding state
%  Short REST bouts surrounded by RUN → converted to RUN
%  Short RUN bouts surrounded by REST → converted to REST

is_running_filtered = is_running(:);
n_samples = length(is_running);

% Find state transitions
diff_state = [0; diff(double(is_running_filtered))];
transitions = find(diff_state ~= 0);

if isempty(transitions)
    return;  % Single state throughout
end

% Add boundaries
transitions = [1; transitions; n_samples + 1];

% Process each bout
for i = 1:(length(transitions) - 1)
    bout_start = transitions(i);
    bout_end = transitions(i + 1) - 1;
    bout_length = bout_end - bout_start + 1;
    
    if bout_length < min_samples
        % Merge into surrounding state
        if bout_start > 1
            surrounding_state = is_running_filtered(bout_start - 1);
        elseif bout_end < n_samples
            surrounding_state = is_running_filtered(bout_end + 1);
        else
            continue;  % Cannot determine surrounding state
        end
        is_running_filtered(bout_start:bout_end) = surrounding_state;
    end
end

end

%% ============================================================================
%  HELPER FUNCTION: Filter Short Bouts by Exclusion (for Clear Mode)
%  ============================================================================
function state_filtered = filter_short_bouts_exclude(state, min_samples)
%FILTER_SHORT_BOUTS_EXCLUDE Remove short bouts by setting them to false
%  Unlike merge_short_bouts, this just excludes short bouts entirely

state_filtered = state(:);
n_samples = length(state);

% Find bout boundaries using connected components approach
labeled = bwlabel(state_filtered);
n_bouts = max(labeled);

for bout_id = 1:n_bouts
    bout_mask = (labeled == bout_id);
    bout_length = sum(bout_mask);
    
    if bout_length < min_samples
        state_filtered(bout_mask) = false;
    end
end

end

%% ============================================================================
%  HELPER FUNCTION: Count Bouts
%  ============================================================================
function [n_rest_bouts, n_run_bouts] = count_bouts(is_rest, is_run)
%COUNT_BOUTS Count the number of contiguous bouts for each state

% Count REST bouts
rest_diff = diff([0; is_rest(:); 0]);
n_rest_bouts = sum(rest_diff == 1);

% Count RUN bouts
run_diff = diff([0; is_run(:); 0]);
n_run_bouts = sum(run_diff == 1);

end

