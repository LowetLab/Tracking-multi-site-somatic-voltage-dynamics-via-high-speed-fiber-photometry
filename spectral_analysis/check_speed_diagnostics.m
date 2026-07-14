%% ============================================================================
%  SPEED DIAGNOSTICS SCRIPT
%  ============================================================================
%  This script checks the actual speed values in a trial to diagnose why
%  RUN classification might be showing 0% even when motion is visible.
%
%  USAGE:
%    1. Modify the QUICK CONFIG section below
%    2. Run: check_speed_diagnostics
%
%  ============================================================================

close all; clear; clc;

%% ============================================================================
%  QUICK CONFIG - MODIFY THESE
%  ============================================================================

% Animal and session to check
ANIMAL_ID = 'Animal01';
SESSION_ID = '01_01_25-R1';
TRIAL_NUM = 1;  % Which trial to check (1-based)

% RUN threshold (should match run_spectral_pipeline.m)
RUN_THRESHOLD_CMS = 2.0;  % cm/s

% Motion conversion constants (MUST MATCH pipeline)
WHEEL_DIAMETER_CM = 19.0;
ENCODER_COUNTS_PER_REV = 1024;
EPHYS_SAMPLING_RATE = 30000;
MOTION_SMOOTH_SAMPLES = 10;

%% ============================================================================
%  LOAD ANIMAL DATABASE
%  ============================================================================

addpath(fullfile(fileparts(mfilename('fullpath')), 'config'));
animals = animal_session_database();

% Find the animal
animal_idx = find(strcmp({animals.mouse_id}, ANIMAL_ID));
if isempty(animal_idx)
    error('Animal "%s" not found in database', ANIMAL_ID);
end
animal = animals(animal_idx);

% Find the session
session_idx = find(strcmp({animal.sessions.session_id}, SESSION_ID));
if isempty(session_idx)
    error('Session "%s" not found for animal "%s"', SESSION_ID, ANIMAL_ID);
end
session = animal.sessions(session_idx);

% Get trial path
if TRIAL_NUM > length(session.trial_paths)
    error('Trial %d not found. Session has %d trials.', TRIAL_NUM, length(session.trial_paths));
end
trial_path = session.trial_paths{TRIAL_NUM};

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  SPEED DIAGNOSTICS\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  Animal: %s\n', ANIMAL_ID);
fprintf('  Session: %s\n', SESSION_ID);
fprintf('  Trial: %d\n', TRIAL_NUM);
fprintf('  File: %s\n', trial_path);
fprintf('  RUN Threshold: %.1f cm/s\n', RUN_THRESHOLD_CMS);
fprintf('════════════════════════════════════════════════════════════════════════\n\n');

%% ============================================================================
%  LOAD TRIAL DATA
%  ============================================================================

if ~exist(trial_path, 'file')
    error('Trial file not found: %s', trial_path);
end

fprintf('Loading trial data...\n');
loaded = load(trial_path);

if ~isfield(loaded, 'FiberPhotometryAnalysis')
    error('FiberPhotometryAnalysis structure not found in %s', trial_path);
end

FPA = loaded.FiberPhotometryAnalysis;

% Extract time vector
if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
    t = FPA.time.time_vector_seconds(:);
else
    error('Time vector not found');
end

% Extract sampling rate
if isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
    fs = FPA.time.sampling_rate;
elseif isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
    fs = FPA.parameters.sampling_rate;
else
    fs = 1 / median(diff(t));
    fprintf('  Warning: Sampling rate not found, inferred as %.1f Hz\n', fs);
end

% Extract motion trace
if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
    motion_raw = FPA.ephys.running_velocity_smooth(:);
    motion_source = 'running_velocity_smooth';
elseif isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity')
    motion_raw = FPA.ephys.running_velocity(:);
    motion_source = 'running_velocity';
else
    error('Motion trace not found');
end

fprintf('  Motion source: %s\n', motion_source);
fprintf('  Sampling rate: %.1f Hz\n', fs);
fprintf('  Duration: %.1f seconds\n', t(end) - t(1));
fprintf('  Total samples: %d\n', length(t));

%% ============================================================================
%  CONVERT MOTION TO SPEED (cm/s)
%  ============================================================================

fprintf('\nConverting motion to speed (cm/s)...\n');

% Motion conversion formula (EXACT match to pipeline)
WHEEL_CIRCUMFERENCE_CM = pi * WHEEL_DIAMETER_CM;
DISTANCE_PER_EDGE_CM = WHEEL_CIRCUMFERENCE_CM / ENCODER_COUNTS_PER_REV;
MOTION_TO_CM_PER_S = (EPHYS_SAMPLING_RATE / 1000) * DISTANCE_PER_EDGE_CM;

fprintf('  Wheel diameter: %.1f cm\n', WHEEL_DIAMETER_CM);
fprintf('  Wheel circumference: %.2f cm\n', WHEEL_CIRCUMFERENCE_CM);
fprintf('  Encoder counts per rev: %d\n', ENCODER_COUNTS_PER_REV);
fprintf('  Distance per edge: %.6f cm\n', DISTANCE_PER_EDGE_CM);
fprintf('  Conversion factor: %.6f cm/s per encoder count\n', MOTION_TO_CM_PER_S);

% Convert to speed
speed_cm_s = motion_raw * MOTION_TO_CM_PER_S;

% Apply smoothing (if enabled)
if MOTION_SMOOTH_SAMPLES > 1
    kernel = ones(MOTION_SMOOTH_SAMPLES, 1) / MOTION_SMOOTH_SAMPLES;
    speed_cm_s = conv(speed_cm_s, kernel, 'same');
    fprintf('  Applied smoothing: %d samples\n', MOTION_SMOOTH_SAMPLES);
end

%% ============================================================================
%  COMPUTE STATISTICS
%  ============================================================================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  SPEED STATISTICS\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');

% Basic statistics
fprintf('\nBasic Statistics:\n');
fprintf('  Minimum speed:     %.3f cm/s\n', min(speed_cm_s));
fprintf('  Maximum speed:     %.3f cm/s\n', max(speed_cm_s));
fprintf('  Mean speed:        %.3f cm/s\n', mean(speed_cm_s));
fprintf('  Median speed:      %.3f cm/s\n', median(speed_cm_s));
fprintf('  Std deviation:     %.3f cm/s\n', std(speed_cm_s));

% Percentiles
fprintf('\nPercentiles:\n');
percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99];
for p = percentiles
    val = prctile(speed_cm_s, p);
    fprintf('  %2dth percentile:   %.3f cm/s\n', p, val);
end

% Classification statistics
fprintf('\nClassification (RUN threshold = %.1f cm/s):\n', RUN_THRESHOLD_CMS);
n_total = length(speed_cm_s);
n_above_threshold = sum(speed_cm_s > RUN_THRESHOLD_CMS);
n_below_threshold = sum(speed_cm_s <= RUN_THRESHOLD_CMS);
pct_above = 100 * n_above_threshold / n_total;
pct_below = 100 * n_below_threshold / n_total;

fprintf('  Samples above threshold: %d (%.1f%%)\n', n_above_threshold, pct_above);
fprintf('  Samples below threshold: %d (%.1f%%)\n', n_below_threshold, pct_below);

% Time above threshold
duration_above_sec = n_above_threshold / fs;
duration_below_sec = n_below_threshold / fs;
fprintf('  Time above threshold:    %.2f seconds\n', duration_above_sec);
fprintf('  Time below threshold:    %.2f seconds\n', duration_below_sec);

% Speed distribution
fprintf('\nSpeed Distribution:\n');
speed_ranges = [0, 0.1, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, inf];
for i = 1:(length(speed_ranges) - 1)
    mask = (speed_cm_s >= speed_ranges(i)) & (speed_cm_s < speed_ranges(i+1));
    n_in_range = sum(mask);
    pct_in_range = 100 * n_in_range / n_total;
    if speed_ranges(i+1) == inf
        fprintf('  %.1f - inf cm/s:      %d samples (%.1f%%)\n', speed_ranges(i), n_in_range, pct_in_range);
    else
        fprintf('  %.1f - %.1f cm/s:      %d samples (%.1f%%)\n', speed_ranges(i), speed_ranges(i+1), n_in_range, pct_in_range);
    end
end

%% ============================================================================
%  PLOT SPEED TRACE
%  ============================================================================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  GENERATING PLOT\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');

figure('Position', [100, 100, 1200, 600]);

% Main speed trace
subplot(2, 1, 1);
plot(t, speed_cm_s, 'b-', 'LineWidth', 1);
hold on;
yline(RUN_THRESHOLD_CMS, 'r--', 'LineWidth', 2, 'Label', sprintf('RUN Threshold (%.1f cm/s)', RUN_THRESHOLD_CMS));
xlabel('Time (seconds)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Speed (cm/s)', 'FontSize', 12, 'FontWeight', 'bold');
title(sprintf('Speed Trace: %s %s Trial %d', ANIMAL_ID, SESSION_ID, TRIAL_NUM), ...
    'FontSize', 14, 'FontWeight', 'bold');
grid on;
legend('Speed', 'RUN Threshold', 'Location', 'best');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Histogram
subplot(2, 1, 2);
histogram(speed_cm_s, 100, 'FaceColor', 'b', 'EdgeColor', 'none');
hold on;
xline(RUN_THRESHOLD_CMS, 'r--', 'LineWidth', 2, 'Label', sprintf('RUN Threshold (%.1f cm/s)', RUN_THRESHOLD_CMS));
xlabel('Speed (cm/s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Frequency', 'FontSize', 12, 'FontWeight', 'bold');
title('Speed Distribution', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

sgtitle(sprintf('Speed Diagnostics: %s %s Trial %d', ANIMAL_ID, SESSION_ID, TRIAL_NUM), ...
    'FontSize', 16, 'FontWeight', 'bold');

fprintf('  Plot generated.\n');

%% ============================================================================
%  SUMMARY AND RECOMMENDATIONS
%  ============================================================================

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  SUMMARY AND RECOMMENDATIONS\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');

if pct_above < 0.1
    fprintf('\n⚠️  WARNING: Very little RUN data detected (%.1f%%)\n', pct_above);
    fprintf('   This explains why RUN classification shows 0.0%%.\n');
    fprintf('\n   Possible reasons:\n');
    fprintf('   1. Speed values are genuinely below %.1f cm/s\n', RUN_THRESHOLD_CMS);
    fprintf('   2. Motion conversion parameters may be incorrect\n');
    fprintf('   3. Motion trace may need different processing\n');
    fprintf('\n   Recommendations:\n');
    fprintf('   - Check if max speed (%.3f cm/s) is reasonable for this animal\n', max(speed_cm_s));
    fprintf('   - Consider lowering RUN threshold if motion is clearly visible\n');
    fprintf('   - Verify wheel diameter and encoder settings match your setup\n');
else
    fprintf('\n✓ RUN data detected: %.1f%% of samples above threshold\n', pct_above);
    fprintf('  If classification still shows 0%%, check:\n');
    fprintf('  - Minimum bout duration filtering (may be too strict)\n');
    fprintf('  - Behavior classification mode settings\n');
end

fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
fprintf('  DIAGNOSTICS COMPLETE\n');
fprintf('════════════════════════════════════════════════════════════════════════\n');
