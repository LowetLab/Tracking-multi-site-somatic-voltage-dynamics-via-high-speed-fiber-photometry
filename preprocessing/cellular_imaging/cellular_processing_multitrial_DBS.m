%% ============================================================================
%  COMPREHENSIVE CELLULAR VOLTAGE IMAGING AND LFP ANALYSIS (MULTI-TRIAL VERSION)
%  ============================================================================
%  Processes cellular voltage imaging data with neuron ROI selection and aligns
%  with LFP recordings from Open Ephys. Designed for DBS (Deep Brain Stimulation)
%  analysis with support for comparing different stimulation frequencies.
%
%  FEATURES:
%  - Multi-trial processing (automatically detects _1, _2, _3 trial folders)
%  - Multi-TIFF support (concatenates multiple TIFF files per trial)
%  - Automatic path construction for imaging and Open Ephys data
%  - Neuron ROI selection with polygon drawing
%  - Motion correction (optional)
%  - Photobleaching correction (linear detrending + double exponential)
%  - Spike detection
%  - DBS metadata tracking (frequency, comparison type)
%  - Organized output folder structure matching fiber pipeline
%
%  BASED ON: run_fiber_preprocessing_multitrial.m
%  ADAPTED FOR: Cellular voltage imaging with neuron ROI selection
%  ============================================================================

close all; clear; clc;

%% Add required paths
% External toolboxes via the centralised config (config/lab_paths.m).
% Override per machine via config/paths_local.m. See config/README.md.
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'config'));
setup_lab_paths();

%% ============================================================================
%  SECTION 1: USER CONFIGURATION - MODIFY THESE PARAMETERS ONLY
%  ============================================================================

%% MOUSE AND RECORDING CONFIGURATION  -- EDIT THESE FOR YOUR RECORDING
MOUSE_NAME = 'Animal01';                % Your animal/mouse identifier
RECORDING_DATE = '01-01-25';            % Format: DD-MM-YY (with dashes to match folder)
RECORDING_ID = 'R1';                    % Recording ID (R1, R2, R9, etc.)
EXPERIMENTER = 'YourInitials';          % Experimenter initials/name
EXPERIMENT_TYPE = '';                    % Empty: no experiment type subfolder in this path

%% INDICATOR CONFIGURATION
INDICATOR_PROMOTER = 'CamKII';          % Promoter targeting pyramidal neurons
INDICATOR_POLARITY = 'positive';        % 'positive' = depolarization increases fluorescence

%% DBS STIMULATION PARAMETERS (Baseline - no stimulation)
DBS_FREQUENCY_HZ = 0;                   % No stimulation (baseline recording)
DBS_COMPARISON_TYPE = '';               % N/A for baseline
DBS_DURATION_SEC = 0;                   % No stimulation
DBS_VOLTAGE = 0;                        % No stimulation
DBS_PULSE_WIDTH_US = 0;                 % No stimulation
DBS_CURRENT_UA = [];                    % No stimulation

%% SESSION DATABASE  -- EDIT THIS to your own cohort
% Reference table for all of this animal's DBS sessions (for your own
% bookkeeping -- not read automatically by the script below).
% Format: {Date, RecID, StimFreq, NumTrials, NumNeurons, ComparisonType}
ANIMAL_SESSIONS = {
    '01-06-25', 'R1',  40,   4, 11, 'AmpBalanced';     % 40Hz Amp balanced
    '01-06-25', 'R2',  40,   5, 11, 'EnergyBalanced';  % 40Hz Energy balanced
    '01-06-25', 'R10', 130,  5, 11, 'EnergyBalanced';  % 130Hz (comparison for both)
};
% Note: comparisons are typically run within the SAME amp-/energy-balanced
% condition across stimulation frequencies (e.g. 40Hz vs 130Hz).

%% BASE PATH CONFIGURATION
% Data location - actual folder structure:
% <BASE_PATH_DATA>\{Mouse}\{ExperimentType}\Voltage_Imaging\{Date}\{RecID}\{TrialFolder}\*.tif
BASE_PATH_DATA = fullfile(lab_paths().data_root, 'DBS');

% Construct DATA_FOLDER path automatically
% Path: BASE_PATH_DATA / MOUSE_NAME / EXPERIMENT_TYPE / Voltage_Imaging / RECORDING_DATE / RECORDING_ID
DATA_FOLDER = fullfile(BASE_PATH_DATA, MOUSE_NAME, EXPERIMENT_TYPE, 'Voltage_Imaging', RECORDING_DATE, RECORDING_ID);

% Output location (constructed relative to data source, like fiber code)
% Creates CellularDataProcessed folder at same level as Voltage_Imaging
% Path: BASE_PATH_DATA / MOUSE_NAME / EXPERIMENT_TYPE / CellularDataProcessed / DATE-RECID
experiment_base_path = fullfile(BASE_PATH_DATA, MOUSE_NAME, EXPERIMENT_TYPE);
processed_base = fullfile(experiment_base_path, 'CellularDataProcessed');
session_folder_name = sprintf('%s-%s', RECORDING_DATE, RECORDING_ID);
OUTPUT_FOLDER = fullfile(processed_base, session_folder_name);

%% PROCESSING OPTIONS
MOTION_CORRECTION_ENABLED = true;    % Enable/disable motion correction
MOTION_CORRECTION_TYPE = 'rigid';    % 'rigid' or 'non-rigid'
APPLY_PHOTOBLEACHING_CORRECTION = true;  % Enable photobleaching correction

%% SAMPLING RATES
EPHYS_FS = 30000;                    % Open Ephys sampling rate (Hz) - fixed hardware
IMAGING_FS = [];                     % Will be calculated from camera triggers

%% OPEN EPHYS CONFIGURATION
LOAD_EPHYS_DATA = true;              % Set false to skip ephys loading
EPHYS_FILE_PREFIX = 100;             % File prefix for Open Ephys files (usually 100)

%% CHANNEL DEFINITIONS FOR OPEN EPHYS
CHANNELS = struct();
CHANNELS.LFP = 11;                   % Channel 11: LFP data
CHANNELS.CAMERA_TRIGGER = 1;         % ADC1: Camera frame triggers
CHANNELS.MOTION_X = 2;               % ADC2: Motion data (X-axis)
CHANNELS.MOTION_Y = 3;               % ADC3: Motion data (Y-axis)
CHANNELS.MOTION_Z = 4;               % ADC4: Motion data (Z-axis)
CHANNELS.STIMULUS_OUTPUT = 6;        % ADC6: Stimulus pulses (biphasic)
CHANNELS.TRIAL_ILLUMINATION = 7;     % ADC7: Trial start/end illumination marker

%% DETECTION THRESHOLDS
CAMERA_TRIGGER_THRESHOLD = 0.5;      % Threshold for camera triggers
STIMULUS_TRIGGER_THRESHOLD = 0.1;    % Threshold for stimulus pulses
ILLUMINATION_THRESHOLD = 0.5;        % Threshold for illumination detection

%% SPIKE DETECTION PARAMETERS
SPIKE_DETECTION_TYPE = 'baselineSD'; % 'baselineSD' (SNR-based) or 'thresholdSD' (fixed baseline-window threshold)
SPIKE_PARAMS = struct();
SPIKE_PARAMS.up_threshold = 3.5;       % Threshold for upward deflections
SPIKE_PARAMS.down_threshold = 3.5;     % Threshold for downward deflections
SPIKE_PARAMS.smoothing_factor = 5;   % Smoothing for spike detection
SPIKE_PARAMS.waveform_window = [-2, 6]; % Sample window around spike peak

%% FIGURE GENERATION OPTIONS
GENERATE_FIGURES = true;
SAVE_FIGURES = true;

%% METADATA STRUCTURE
METADATA = struct();
METADATA.mouse_name = MOUSE_NAME;
METADATA.experimenter = EXPERIMENTER;
METADATA.recording_date = RECORDING_DATE;
METADATA.recording_id = RECORDING_ID;
METADATA.experiment_type = EXPERIMENT_TYPE;
METADATA.indicator_promoter = INDICATOR_PROMOTER;
METADATA.indicator_polarity = INDICATOR_POLARITY;
METADATA.dbs_frequency_hz = DBS_FREQUENCY_HZ;
METADATA.dbs_comparison_type = DBS_COMPARISON_TYPE;
METADATA.dbs_duration_sec = DBS_DURATION_SEC;
METADATA.is_baseline = (DBS_FREQUENCY_HZ == 0);
METADATA.processing_date = datestr(now, 'yyyy-mm-dd HH:MM:SS');

%% ============================================================================
%  SECTION 2: PATH VERIFICATION AND TRIAL DETECTION
%  ============================================================================

fprintf('\n');
fprintf('============================================================================\n');
fprintf('  CELLULAR VOLTAGE IMAGING PREPROCESSING PIPELINE\n');
fprintf('============================================================================\n');
fprintf('  Mouse: %s\n', MOUSE_NAME);
fprintf('  Date: %s\n', RECORDING_DATE);
fprintf('  Recording: %s\n', RECORDING_ID);
fprintf('  Experiment Type: %s\n', EXPERIMENT_TYPE);
fprintf('  DBS Frequency: %d Hz\n', DBS_FREQUENCY_HZ);
fprintf('  Comparison Type: %s\n', DBS_COMPARISON_TYPE);
fprintf('============================================================================\n\n');

% Verify data folder exists
if ~exist(DATA_FOLDER, 'dir')
    error('Data folder does not exist: %s\nPlease check configuration parameters.', DATA_FOLDER);
end

fprintf('Data folder: %s\n', DATA_FOLDER);

% Create output folder structure (like fiber code)
% First create CellularDataProcessed folder if needed
if ~exist(processed_base, 'dir')
    mkdir(processed_base);
    fprintf('Created processed data folder: %s\n', processed_base);
end

% Then create session-specific folder
if ~exist(OUTPUT_FOLDER, 'dir')
    mkdir(OUTPUT_FOLDER);
    fprintf('Created session folder: %s\n', OUTPUT_FOLDER);
else
    fprintf('Output folder: %s\n', OUTPUT_FOLDER);
end

%% Detect trial folders
fprintf('\nSearching for trial folders...\n');

folder_contents = dir(DATA_FOLDER);
folder_contents = folder_contents([folder_contents.isdir]);
folder_contents = folder_contents(~ismember({folder_contents.name}, {'.', '..'}));

if isempty(folder_contents)
    error('No trial folders found in: %s', DATA_FOLDER);
end

% Extract trial numbers from folder names (look for _1, _2, _3, etc.)
trial_folders = {};
trial_numbers = [];

for i = 1:length(folder_contents)
    folder_name = folder_contents(i).name;
    % Look for pattern ending with _ followed by digits
    match = regexp(folder_name, '_(\d+)$', 'tokens');
    if ~isempty(match)
        trial_num = str2double(match{1}{1});
        trial_folders{end+1} = folder_name;
        trial_numbers(end+1) = trial_num;
    end
end

if isempty(trial_folders)
    error('No trial folders with _1, _2, _3 suffixes found in: %s', DATA_FOLDER);
end

% Sort by trial number
[trial_numbers, sort_idx] = sort(trial_numbers);
trial_folders = trial_folders(sort_idx);
num_trials = length(trial_folders);

fprintf('Found %d trial folders:\n', num_trials);
for i = 1:num_trials
    fprintf('  Trial %d: %s\n', trial_numbers(i), trial_folders{i});
end

%% ============================================================================
%  SECTION 3: OPEN EPHYS PATH CONSTRUCTION
%  ============================================================================

fprintf('\n=== CONSTRUCTING OPEN EPHYS PATH ===\n');

% Open Ephys data location - try multiple common patterns
% Pattern 1: Parallel to Voltage_Imaging: .../OpenEphys/DATE/RECID/
% Pattern 2: Inside same folder structure
% Pattern 3: In Data/OpenEphys/ base folder
EPHYS_ALT_PATHS = {
    fullfile(BASE_PATH_DATA, MOUSE_NAME, EXPERIMENT_TYPE, 'Open_Ephys', RECORDING_DATE, RECORDING_ID),
    fullfile(BASE_PATH_DATA, MOUSE_NAME, EXPERIMENT_TYPE, RECORDING_DATE, RECORDING_ID, 'OpenEphys'),
    fullfile(DATA_FOLDER, 'OpenEphys'),
    fullfile(lab_paths().data_root, 'OpenEphys', MOUSE_NAME, EXPERIMENT_TYPE, RECORDING_DATE, RECORDING_ID),
    fullfile(lab_paths().data_root, 'OpenEphys', MOUSE_NAME, RECORDING_DATE, RECORDING_ID),
};

% Find the correct ephys path
ephys_found = false;
for path_idx = 1:length(EPHYS_ALT_PATHS)
    test_path = EPHYS_ALT_PATHS{path_idx};
    if exist(test_path, 'dir')
        % Check if it contains .continuous files
        continuous_files = dir(fullfile(test_path, '*.continuous'));
        if ~isempty(continuous_files)
            EPHYS_DATA_FOLDER = test_path;
            ephys_found = true;
            fprintf('Found OpenEphys data at: %s\n', EPHYS_DATA_FOLDER);
            break;
        end
        % Check subdirectories (Record Node pattern)
        subdirs = dir(test_path);
        for sd = 1:length(subdirs)
            if subdirs(sd).isdir && ~ismember(subdirs(sd).name, {'.', '..'})
                subpath = fullfile(test_path, subdirs(sd).name);
                continuous_files = dir(fullfile(subpath, '**', '*.continuous'));
                if ~isempty(continuous_files)
                    % Use the folder containing the .continuous files directly
                    EPHYS_DATA_FOLDER = continuous_files(1).folder;
                    ephys_found = true;
                    fprintf('Found OpenEphys data at: %s\n', EPHYS_DATA_FOLDER);
                    break;
                end
            end
        end
        if ephys_found, break; end
    end
end

if ~ephys_found && LOAD_EPHYS_DATA
    warning('Could not find OpenEphys data automatically. Will prompt for manual selection.');
    [~, EPHYS_DATA_FOLDER] = uigetfile('*.continuous', 'Select any Open Ephys .continuous file');
    if isequal(EPHYS_DATA_FOLDER, 0)
        error('No ephys folder selected. Cannot proceed.');
    end
    EPHYS_DATA_FOLDER = fileparts(EPHYS_DATA_FOLDER);
end

%% ============================================================================
%  SECTION 4: PROCESS EACH TRIAL
%  ============================================================================

% Store results for all trials
all_trial_results = cell(num_trials, 1);

% ROI masks will be shared across trials (selected from first trial)
shared_neuron_roi_masks = [];
shared_neuron_centroids = [];
shared_motion_roi_position = [];
num_neurons = 0;

for trial_idx = 1:num_trials
    fprintf('\n');
    fprintf('########################################################################\n');
    fprintf('  PROCESSING TRIAL %d/%d: %s\n', trial_idx, num_trials, trial_folders{trial_idx});
    fprintf('########################################################################\n');
    
    trial_folder_path = fullfile(DATA_FOLDER, trial_folders{trial_idx});
    
    %% ========================================================================
    %  LOAD VOLTAGE IMAGING DATA
    %% ========================================================================
    
    fprintf('\n=== LOADING VOLTAGE IMAGING DATA ===\n');
    
    % Find all TIFF files in the trial folder
    tiff_files = dir(fullfile(trial_folder_path, '*.tif'));
    if isempty(tiff_files)
        tiff_files = dir(fullfile(trial_folder_path, '*.tiff'));
    end
    
    if isempty(tiff_files)
        warning('No TIFF files found in trial folder: %s. Skipping trial.', trial_folder_path);
        continue;
    end
    
    % Sort TIFF files by name (to ensure correct concatenation order)
    [~, sort_idx_tiff] = sort({tiff_files.name});
    tiff_files = tiff_files(sort_idx_tiff);
    
    fprintf('Found %d TIFF file(s) in trial folder\n', length(tiff_files));
    
    % Load and concatenate all TIFF files
    voltage_imaging_stack = [];
    
    for tiff_idx = 1:length(tiff_files)
        tiff_path = fullfile(tiff_files(tiff_idx).folder, tiff_files(tiff_idx).name);
        fprintf('  Loading: %s\n', tiff_files(tiff_idx).name);
        
        % Get TIFF info
        tiff_info = imfinfo(tiff_path);
        num_frames_this_file = numel(tiff_info);
        frame_height = tiff_info(1).Height;
        frame_width = tiff_info(1).Width;
        
        fprintf('    Frames: %d, Size: %dx%d\n', num_frames_this_file, frame_height, frame_width);
        
        % Load frames
        tiff_file_obj = Tiff(tiff_path, 'r');
        temp_stack = zeros(frame_height, frame_width, num_frames_this_file, 'uint16');
        
        warning('off', 'all');
        for frame_idx = 1:num_frames_this_file
            if mod(frame_idx, 500) == 0
                fprintf('    Frame %d/%d\n', frame_idx, num_frames_this_file);
            end
            tiff_file_obj.setDirectory(frame_idx);
            temp_stack(:,:,frame_idx) = tiff_file_obj.read();
        end
        tiff_file_obj.close();
        warning('on', 'all');
        
        % Concatenate to main stack
        if isempty(voltage_imaging_stack)
            voltage_imaging_stack = temp_stack;
        else
            % Verify dimensions match
            if size(temp_stack,1) ~= size(voltage_imaging_stack,1) || ...
               size(temp_stack,2) ~= size(voltage_imaging_stack,2)
                error('TIFF dimension mismatch between files');
            end
            voltage_imaging_stack = cat(3, voltage_imaging_stack, temp_stack);
        end
    end
    
    num_imaging_frames = size(voltage_imaging_stack, 3);
    frame_height = size(voltage_imaging_stack, 1);
    frame_width = size(voltage_imaging_stack, 2);
    
    fprintf('Total imaging stack: %d frames, %dx%d pixels\n', num_imaging_frames, frame_height, frame_width);
    
    %% ========================================================================
    %  MOTION CORRECTION ROI SELECTION (First trial only)
    %% ========================================================================
    
    if trial_idx == 1
        fprintf('\n=== MOTION CORRECTION ROI SELECTION ===\n');
        
        % Create reference frame
        reference_frame_start = min(100, num_imaging_frames);
        reference_frame_end = min(190, num_imaging_frames);
        reference_frame = nanmean(voltage_imaging_stack(:,:,reference_frame_start:reference_frame_end), 3);
        
        figure('Name', 'Motion Correction ROI Selection', 'Position', [100, 100, 800, 600]);
        imagesc(reference_frame);
        axis image; axis off; colormap(gray);
        title('Select ROI for Motion Correction (drag rectangle)');
        drawnow;
        
        motion_correction_roi_handle = imrect;
        if isvalid(motion_correction_roi_handle)
            shared_motion_roi_position = round(getPosition(motion_correction_roi_handle));
            close(gcf);
        else
            error('Motion correction ROI not selected properly.');
        end
        
        fprintf('Motion ROI: [%d,%d] size [%d,%d]\n', ...
            shared_motion_roi_position(1), shared_motion_roi_position(2), ...
            shared_motion_roi_position(3), shared_motion_roi_position(4));
    end
    
    %% ========================================================================
    %  APPLY MOTION CORRECTION
    %% ========================================================================
    
    % Extract motion ROI from stack
    motion_roi_stack = voltage_imaging_stack(...
        shared_motion_roi_position(2):min(shared_motion_roi_position(2)+shared_motion_roi_position(4), frame_height), ...
        shared_motion_roi_position(1):min(shared_motion_roi_position(1)+shared_motion_roi_position(3), frame_width), :);
    
    if MOTION_CORRECTION_ENABLED
        fprintf('\n=== APPLYING MOTION CORRECTION (Trial %d/%d) ===\n', trial_idx, num_trials);
        
        motion_corrected_stack = single(motion_roi_stack);
        
        % 3D smoothing
        smoothed_stack = imboxfilt3(motion_corrected_stack, [1 1 11]);
        
        % High-pass filter for motion detection
        gaussian_small = fspecial('gaussian', 50, 1);
        gaussian_large = fspecial('gaussian', 50, 25);
        motion_detection_filter = gaussian_small - gaussian_large;
        filtered_stack = gather(imfilter(smoothed_stack, motion_detection_filter, 'replicate', 'same'));
        
        % NoRMCorre parameters
        motion_correction_options = NoRMCorreSetParms(...
            'd1', size(motion_corrected_stack,1), ...
            'd2', size(motion_corrected_stack,2), ...
            'bin_width', 10, ...
            'max_shift', 50, ...
            'us_fac', 1);
        
        % Run motion estimation
        [~, motion_shifts, ~] = normcorre(filtered_stack, motion_correction_options);
        
        % Apply shifts to each frame
        corrected_imaging_stack = zeros(size(motion_corrected_stack), 'single');
        for frame_idx = 1:size(motion_corrected_stack, 3)
            corrected_imaging_stack(:,:,frame_idx) = circshift(...
                motion_corrected_stack(:,:,frame_idx), motion_shifts(frame_idx).shifts);
        end
        corrected_imaging_stack = uint16(corrected_imaging_stack);
        
        % Store motion shifts for this trial
        trial_motion_shifts = motion_shifts;
        
        fprintf('Motion correction applied.\n');
    else
        fprintf('Skipping motion correction.\n');
        corrected_imaging_stack = uint16(motion_roi_stack);
        trial_motion_shifts = [];
    end
    
    % Transpose for proper orientation
    corrected_imaging_stack = permute(corrected_imaging_stack, [2 1 3]);
    fprintf('Processed stack: %dx%dx%d\n', size(corrected_imaging_stack));
    
    %% ========================================================================
    %  NEURON ROI SELECTION (First trial only)
    %% ========================================================================
    
    if trial_idx == 1
        fprintf('\n=== NEURON ROI SELECTION ===\n');
        
        average_frame_for_rois = mean(corrected_imaging_stack, 3);
        
        % Remove outliers for visualization
        outlier_threshold = 15;
        outlier_mask = abs(zscore(average_frame_for_rois(:))) > outlier_threshold;
        average_frame_clean = average_frame_for_rois;
        average_frame_clean(outlier_mask) = NaN;
        
        display_min = min(average_frame_clean(:));
        display_max = prctile(average_frame_clean(:), 99.9);
        
        average_frame_clean(isnan(average_frame_clean)) = median(average_frame_clean(:), 'omitnan');
        
        % Initialize ROI collection
        shared_neuron_roi_masks = {};
        combined_roi_mask = zeros(size(average_frame_for_rois));
        
        figure('Name', 'Neuron ROI Selection', 'Position', [100, 100, 900, 700]);
        imagesc(average_frame_clean, [display_min, display_max]);
        axis image; axis off; colormap(gray);
        colorbar;
        title('Draw polygons around neurons (press Enter when done)');
        drawnow;
        
        roi_polygon = drawpolygon;
        roi_counter = 0;
        
        while isvalid(roi_polygon)
            roi_counter = roi_counter + 1;
            
            current_roi_mask = createMask(roi_polygon);
            shared_neuron_roi_masks{end+1} = current_roi_mask;
            combined_roi_mask = double(combined_roi_mask | current_roi_mask);
            
            imagesc(average_frame_clean .* (1 - combined_roi_mask * 0.3), [display_min, display_max]);
            axis image; axis off; colormap(gray);
            colorbar;
            title(sprintf('ROIs selected: %d (draw next or press Enter when done)', roi_counter));
            drawnow;
            
            roi_polygon = drawpolygon;
        end
        
        close(gcf);
        
        num_neurons = length(shared_neuron_roi_masks);
        fprintf('Selected %d neuron ROIs\n', num_neurons);
        
        if num_neurons == 0
            error('No neuron ROIs selected. Cannot proceed.');
        end
        
        % Calculate centroids
        shared_neuron_centroids = zeros(num_neurons, 2);
        for neuron_idx = 1:num_neurons
            roi_properties = regionprops(shared_neuron_roi_masks{neuron_idx}, 'Centroid');
            shared_neuron_centroids(neuron_idx,:) = roi_properties.Centroid;
        end
        
        % Save ROI reference frame
        roi_reference_frame = average_frame_clean;
    end
    
    %% ========================================================================
    %  EXTRACT FLUORESCENCE TRACES
    %% ========================================================================
    
    fprintf('\n=== EXTRACTING FLUORESCENCE TRACES ===\n');
    
    neuron_fluorescence_traces = zeros(num_imaging_frames, num_neurons);
    
    for neuron_idx = 1:num_neurons
        masked_stack = corrected_imaging_stack .* uint16(shared_neuron_roi_masks{neuron_idx});
        roi_pixel_count = sum(shared_neuron_roi_masks{neuron_idx}(:));
        neuron_fluorescence_traces(:, neuron_idx) = squeeze(sum(masked_stack, [1, 2])) / roi_pixel_count;
    end
    
    % Background trace
    background_mask = true(size(corrected_imaging_stack, 1), size(corrected_imaging_stack, 2));
    for neuron_idx = 1:num_neurons
        background_mask = background_mask & ~shared_neuron_roi_masks{neuron_idx};
    end
    masked_background = corrected_imaging_stack .* uint16(background_mask);
    background_pixel_count = sum(background_mask(:));
    background_trace = squeeze(sum(masked_background, [1, 2])) / background_pixel_count;
    
    fprintf('Extracted traces for %d neurons\n', num_neurons);
    
    %% ========================================================================
    %  LOAD OPEN EPHYS DATA
    %% ========================================================================
    
    if LOAD_EPHYS_DATA
        fprintf('\n=== LOADING OPEN EPHYS DATA ===\n');
        
        % Change to ephys directory (like fiber code - load_open_ephys_data works better this way)
        original_dir = pwd;
        
        % Verify ephys folder exists
        if ~exist(EPHYS_DATA_FOLDER, 'dir')
            error('EPHYS_DATA_FOLDER does not exist: %s', EPHYS_DATA_FOLDER);
        end
        
        cd(EPHYS_DATA_FOLDER);
        fprintf('Changed to ephys folder: %s\n', EPHYS_DATA_FOLDER);
        
        % Construct filenames (not full paths, since we're in the directory)
        lfp_filename = sprintf('%d_RhythmData_Ch%d.continuous', EPHYS_FILE_PREFIX, CHANNELS.LFP);
        camera_trigger_filename = sprintf('%d_RhythmData_ADC%d.continuous', EPHYS_FILE_PREFIX, CHANNELS.CAMERA_TRIGGER);
        trial_illumination_filename = sprintf('%d_RhythmData_ADC%d.continuous', EPHYS_FILE_PREFIX, CHANNELS.TRIAL_ILLUMINATION);
        stimulus_output_filename = sprintf('%d_RhythmData_ADC%d.continuous', EPHYS_FILE_PREFIX, CHANNELS.STIMULUS_OUTPUT);
        
        % Verify files exist before loading
        if ~exist(lfp_filename, 'file')
            cd(original_dir);
            error('LFP file not found: %s in %s', lfp_filename, EPHYS_DATA_FOLDER);
        end
        
        % Load channels
        fprintf('Loading LFP (Ch%d)...\n', CHANNELS.LFP);
        [lfp_data_raw, lfp_timestamps, ~] = load_open_ephys_data(lfp_filename);
        
        fprintf('Loading Camera Triggers (ADC%d)...\n', CHANNELS.CAMERA_TRIGGER);
        [camera_trigger_raw, ~, ~] = load_open_ephys_data(camera_trigger_filename);
        
        fprintf('Loading Trial Illumination (ADC%d)...\n', CHANNELS.TRIAL_ILLUMINATION);
        [trial_illumination_raw, ~, ~] = load_open_ephys_data(trial_illumination_filename);
        
        fprintf('Loading Stimulus Output (ADC%d)...\n', CHANNELS.STIMULUS_OUTPUT);
        [stimulus_output_raw, ~, ~] = load_open_ephys_data(stimulus_output_filename);
        
        % Return to original directory
        cd(original_dir);
        
        % Align to minimum length
        min_ephys_samples = min([length(lfp_data_raw), length(camera_trigger_raw), ...
                                length(trial_illumination_raw), length(stimulus_output_raw)]);
        
        lfp_data_raw = lfp_data_raw(1:min_ephys_samples);
        camera_trigger_raw = camera_trigger_raw(1:min_ephys_samples);
        trial_illumination_raw = trial_illumination_raw(1:min_ephys_samples);
        stimulus_output_raw = stimulus_output_raw(1:min_ephys_samples);
        
        fprintf('Ephys data aligned: %d samples\n', min_ephys_samples);
        
        %% Detect camera triggers and compute frame rate (fiber code approach)
        fprintf('\n=== DETECTING CAMERA TRIGGERS ===\n');
        
        % Use positive-going edge detection (like fiber code)
        all_camera_triggers = find(diff(camera_trigger_raw) > CAMERA_TRIGGER_THRESHOLD);
        fprintf('Total trigger edges detected: %d\n', length(all_camera_triggers));
        
        % Filter out spurious triggers (minimum spacing filter, like fiber code)
        % At 650 Hz, minimum interval is ~46 samples at 30kHz
        min_trigger_spacing = 8;  % Samples at EPHYS_FS
        trigger_intervals_raw = diff(all_camera_triggers);
        valid_trigger_mask = [true; trigger_intervals_raw > min_trigger_spacing];
        camera_frame_indices = all_camera_triggers(valid_trigger_mask);
        
        fprintf('Camera frames detected (after filtering): %d\n', length(camera_frame_indices));
        
        % Compute frame rate from trigger intervals
        if length(camera_frame_indices) > 1
            trigger_intervals = diff(camera_frame_indices);
            trigger_intervals_sec = trigger_intervals / EPHYS_FS;
            instantaneous_frame_rates = 1 ./ trigger_intervals_sec;
            
            median_frame_rate = median(instantaneous_frame_rates);
            mean_frame_rate = mean(instantaneous_frame_rates);
            std_frame_rate = std(instantaneous_frame_rates);
            
            fprintf('  Median frame rate: %.2f Hz\n', median_frame_rate);
            fprintf('  Mean frame rate: %.2f Hz\n', mean_frame_rate);
            fprintf('  Std deviation: %.2f Hz\n', std_frame_rate);
            
            IMAGING_FS = median_frame_rate;
            
            if std_frame_rate > 50.0
                warning('High frame rate variability detected (std = %.2f Hz)', std_frame_rate);
            end
        else
            % Fallback to manual frame rate
            MANUAL_FRAME_RATE = 650;
            fprintf('Insufficient triggers. Using manual frame rate: %.1f Hz\n', MANUAL_FRAME_RATE);
            IMAGING_FS = MANUAL_FRAME_RATE;
            
            illumination_on = trial_illumination_raw > ILLUMINATION_THRESHOLD;
            illumination_start = find(illumination_on, 1, 'first');
            illumination_end = find(illumination_on, 1, 'last');
            
            if ~isempty(illumination_start)
                illumination_duration_sec = (illumination_end - illumination_start) / EPHYS_FS;
                expected_frames = round(illumination_duration_sec * MANUAL_FRAME_RATE);
                camera_frame_indices = round(linspace(illumination_start, illumination_end, expected_frames));
            end
        end
        
        %% Detect trial boundaries and stimulus (fiber code approach)
        fprintf('\n=== DETECTING STIMULATION PERIOD ===\n');
        
        trial_illumination_on = find(diff(trial_illumination_raw) > ILLUMINATION_THRESHOLD);
        trial_illumination_off = find(diff(trial_illumination_raw) < -ILLUMINATION_THRESHOLD);
        
        if ~isempty(trial_illumination_on) && ~isempty(trial_illumination_off)
            trial_start_ephys = trial_illumination_on(1);
            trial_end_ephys = trial_illumination_off(end);
            fprintf('Trial period: samples %d to %d (%.2f to %.2f sec)\n', ...
                trial_start_ephys, trial_end_ephys, ...
                trial_start_ephys/EPHYS_FS, trial_end_ephys/EPHYS_FS);
        else
            trial_start_ephys = 1;
            trial_end_ephys = length(trial_illumination_raw);
        end
        
        % Detect stimulus onset from pulse channel (like fiber code)
        STIM_ONSET_THRESHOLD = 0.1;  % Match fiber code threshold
        stimulus_diff = diff(stimulus_output_raw);
        onset_candidates = find(stimulus_diff > STIM_ONSET_THRESHOLD);
        
        if ~isempty(onset_candidates)
            % Filter to only include onsets within trial period
            valid_onsets = onset_candidates(onset_candidates > trial_start_ephys & ...
                                            onset_candidates < trial_end_ephys);
            if ~isempty(valid_onsets)
                stimulus_onset_ephys = valid_onsets(1);
                stim_onset_time = stimulus_onset_ephys / EPHYS_FS;
                
                % Calculate stimulus offset using duration
                stimulus_offset_ephys = stimulus_onset_ephys + round(DBS_DURATION_SEC * EPHYS_FS);
                stim_offset_time = stimulus_offset_ephys / EPHYS_FS;
                
                fprintf('  Stimulus onset detected at sample: %d (%.3f sec)\n', stimulus_onset_ephys, stim_onset_time);
                fprintf('  Stimulus offset set to sample: %d (%.3f sec)\n', stimulus_offset_ephys, stim_offset_time);
                fprintf('  Stimulus duration: %.2f sec\n', DBS_DURATION_SEC);
            else
                stimulus_onset_ephys = [];
                stimulus_offset_ephys = [];
                fprintf('  No stimulus detected within trial period\n');
            end
        else
            stimulus_onset_ephys = [];
            stimulus_offset_ephys = [];
            fprintf('  No stimulus pulses detected\n');
        end
        
        %% Align imaging with ephys
        fprintf('\n=== ALIGNING DATA ===\n');
        
        frames_to_align = min(length(camera_frame_indices), num_imaging_frames);
        aligned_camera_indices = camera_frame_indices(1:frames_to_align);
        
        lfp_aligned = lfp_data_raw(aligned_camera_indices);
        stimulus_output_smoothed = fastsmooth(stimulus_output_raw, 31, 1, 1);
        stimulus_aligned = stimulus_output_smoothed(aligned_camera_indices);
        
        neuron_traces_aligned = neuron_fluorescence_traces(1:frames_to_align, :);
        
        % Convert stimulus onset/offset to frame indices (like fiber code)
        if ~isempty(stimulus_onset_ephys)
            [~, stimulus_onset_frame] = min(abs(aligned_camera_indices - stimulus_onset_ephys));
            time_reference_frame = stimulus_onset_frame;
            time_reference_label = 'stimulus onset';
            
            % Also convert stimulus offset
            if exist('stimulus_offset_ephys', 'var') && ~isempty(stimulus_offset_ephys)
                [~, stimulus_offset_frame] = min(abs(aligned_camera_indices - stimulus_offset_ephys));
            else
                stimulus_offset_frame = stimulus_onset_frame + round(DBS_DURATION_SEC * IMAGING_FS);
            end
            
            fprintf('  Stimulus onset frame: %d (%.3f sec)\n', stimulus_onset_frame, stimulus_onset_frame/IMAGING_FS);
            fprintf('  Stimulus offset frame: %d (%.3f sec)\n', stimulus_offset_frame, stimulus_offset_frame/IMAGING_FS);
        else
            stimulus_onset_frame = [];
            stimulus_offset_frame = [];
            time_reference_frame = 1;
            time_reference_label = 'recording start';
        end
        
        time_vector = ((1:frames_to_align) - time_reference_frame) / IMAGING_FS;
        
        fprintf('Aligned %d frames\n', frames_to_align);
        fprintf('Time reference: %s\n', time_reference_label);
        
    else
        % No ephys - create simple time vector
        IMAGING_FS = 650;  % Default
        frames_to_align = num_imaging_frames;
        time_vector = (1:frames_to_align) / IMAGING_FS;
        lfp_aligned = [];
        stimulus_aligned = [];
        stimulus_onset_frame = [];
        neuron_traces_aligned = neuron_fluorescence_traces;
    end
    
    %% ========================================================================
    %  PHOTOBLEACHING CORRECTION
    %% ========================================================================
    
    if APPLY_PHOTOBLEACHING_CORRECTION
        fprintf('\n=== APPLYING PHOTOBLEACHING CORRECTION (Trial %d/%d) ===\n', trial_idx, num_trials);
        
        % Use pre-calculated stimulus onset/offset frames
        if ~isempty(stimulus_onset_frame) && exist('stimulus_offset_frame', 'var') && ~isempty(stimulus_offset_frame)
            stim_onset = stimulus_onset_frame;
            stim_offset = stimulus_offset_frame;
            fprintf('  Stimulus period: frames %d to %d\n', stim_onset, stim_offset);
        elseif ~isempty(stimulus_onset_frame)
            stim_onset = stimulus_onset_frame;
            stim_offset = stim_onset + round(DBS_DURATION_SEC * IMAGING_FS);
            fprintf('  Stimulus period: frames %d to %d (offset from duration)\n', stim_onset, stim_offset);
        else
            % No stimulus detected - use near-end of recording for baseline fitting
            stim_onset = round(frames_to_align * 0.95);
            stim_offset = round(frames_to_align * 0.99);
            fprintf('  No stimulus detected. Using frames %d to %d as "stimulus period"\n', ...
                    stim_onset, stim_offset);
        end
        
        fluorescence_detrended = neuron_traces_aligned;
        fluorescence_exp_corrected = neuron_traces_aligned;
        
        % METHOD 1: Linear detrending using pre-stimulus period only
        fprintf('  Applying linear detrending correction...\n');
        for neuron_idx = 1:num_neurons
            trace = neuron_traces_aligned(:, neuron_idx);
            
            % Use only pre-stimulus period for trend estimation
            pre_stim_idx = 1:stim_onset;
            pre_stim_values = trace(pre_stim_idx);
            
            % Fit linear trend to pre-stimulus data
            trend_coef = polyfit(pre_stim_idx, pre_stim_values, 1);
            
            % Apply detrending to entire trace
            trend_line = polyval(trend_coef, 1:frames_to_align);
            detrended = trace - trend_line';
            
            % Normalize to baseline
            baseline = mean(detrended(1:min(stim_onset-1, 100)));
            if baseline ~= 0
                fluorescence_detrended(:, neuron_idx) = detrended / baseline;
            else
                fluorescence_detrended(:, neuron_idx) = detrended;
            end
        end
        
        % METHOD 2: Double exponential fitting using pre-stimulus period
        fprintf('  Applying double exponential correction...\n');
        double_exp_func = @(p, t) p(1) * exp(-t/p(2)) + p(3) * exp(-t/p(4)) + p(5);
        exp_fit_failures = 0;
        
        % Store F0 values and fitted parameters for each neuron
        F0_values = zeros(1, num_neurons);
        exp_fit_params = cell(num_neurons, 1);
        
        for neuron_idx = 1:num_neurons
            trace = neuron_traces_aligned(:, neuron_idx);
            pre_stim_trace = trace(1:stim_onset);
            time_pre = (0:length(pre_stim_trace)-1)' / IMAGING_FS;
            time_full = (0:frames_to_align-1)' / IMAGING_FS;
            
            % Initial parameter guess [A1, tau1, A2, tau2, offset]
            init_params = [0.5*max(pre_stim_trace), time_pre(end)/3, ...
                          0.5*max(pre_stim_trace), time_pre(end), min(pre_stim_trace)];
            lb = [0, 0, 0, 0, 0];
            ub = [Inf, Inf, Inf, Inf, Inf];
            opts = optimoptions('lsqcurvefit', 'Display', 'off');
            
            try
                fitted_params = lsqcurvefit(double_exp_func, init_params, time_pre, pre_stim_trace, lb, ub, opts);
                exp_fit_params{neuron_idx} = fitted_params;
                fitted_curve = double_exp_func(fitted_params, time_full);
                
                % Correct trace and compute ΔF/F
                corrected = trace ./ fitted_curve;
                baseline_f0 = mean(corrected(1:min(stim_onset-1, 100)));
                F0_values(neuron_idx) = baseline_f0;
                
                if baseline_f0 ~= 0
                    delta_f = (corrected - baseline_f0) / baseline_f0;
                    fluorescence_exp_corrected(:, neuron_idx) = delta_f;
                else
                    fluorescence_exp_corrected(:, neuron_idx) = fluorescence_detrended(:, neuron_idx);
                end
            catch ME
                exp_fit_failures = exp_fit_failures + 1;
                % Use detrended trace as fallback
                fluorescence_exp_corrected(:, neuron_idx) = fluorescence_detrended(:, neuron_idx);
                F0_values(neuron_idx) = NaN;
            end
        end
        
        % Store baseline window used for ΔF/F
        deltaF_F_baseline_window = [1, min(stim_onset-1, 100)];
        deltaF_F_baseline_time = deltaF_F_baseline_window / IMAGING_FS;
        
        if exp_fit_failures > 0
            warning('Double exponential fitting failed for %d/%d neurons, using linear detrending as fallback', ...
                    exp_fit_failures, num_neurons);
        end
        
        % Compute z-scored traces (for compatibility with fiber pipeline)
        fluorescence_zscored = zeros(size(fluorescence_exp_corrected));
        for neuron_idx = 1:num_neurons
            fluorescence_zscored(:, neuron_idx) = zscore(fluorescence_exp_corrected(:, neuron_idx));
        end
        
        fprintf('  Photobleaching correction completed for %d neurons\n', num_neurons);
        fprintf('  ΔF/F baseline window: frames %d to %d (%.3f to %.3f s)\n', ...
                deltaF_F_baseline_window(1), deltaF_F_baseline_window(2), ...
                deltaF_F_baseline_time(1), deltaF_F_baseline_time(2));
    else
        fprintf('\n=== SKIPPING PHOTOBLEACHING CORRECTION (Trial %d/%d) ===\n', trial_idx, num_trials);
        fluorescence_exp_corrected = neuron_traces_aligned;
        fluorescence_detrended = neuron_traces_aligned;
        fluorescence_zscored = zeros(size(neuron_traces_aligned));
        for neuron_idx = 1:num_neurons
            fluorescence_zscored(:, neuron_idx) = zscore(neuron_traces_aligned(:, neuron_idx));
        end
        F0_values = NaN(1, num_neurons);
        exp_fit_params = cell(num_neurons, 1);
        deltaF_F_baseline_window = [];
        deltaF_F_baseline_time = [];
    end
    
    %% ========================================================================
    %  SPIKE DETECTION
    %% ========================================================================
    
    fprintf('\n=== DETECTING SPIKES (Trial %d/%d) ===\n', trial_idx, num_trials);
    
    spike_detection_results = cell(num_neurons, 1);
    
    % Define baseline window for spike detection
    % For DBS trials: use pre-stimulus period as baseline (more accurate)
    % Fallback to first 30% only if no stimulus detected
    if ~isempty(stimulus_onset_frame) && stimulus_onset_frame > 30
        % Use pre-stimulus period as baseline (leave small buffer before stim)
        baseline_window = [1, stimulus_onset_frame - 10];
        fprintf('PRE-STIM BASELINE: Using frames 1 to %d (%.2f sec before stim onset)\n', ...
            baseline_window(2), baseline_window(2)/IMAGING_FS);
    else
        baseline_window = [1, round(frames_to_align * 0.3)];
        fprintf('BASELINE MODE: Using first 30%% as baseline (frames %d to %d)\n', ...
            baseline_window(1), baseline_window(2));
    end
    
    % Detect spikes in each neuron
    for neuron_idx = 1:num_neurons
        if mod(neuron_idx, 5) == 0 || neuron_idx == 1
            fprintf('  Detecting spikes in neuron %d/%d...\n', neuron_idx, num_neurons);
        end
        
        trace = fluorescence_exp_corrected(:, neuron_idx);
        
        try
            if strcmp(SPIKE_DETECTION_TYPE, 'thresholdSD')
                spike_detection_results{neuron_idx} = spike_detect_baseline_threshold_SC(...
                    trace, baseline_window, ...
                    SPIKE_PARAMS.up_threshold, SPIKE_PARAMS.down_threshold, ...
                    SPIKE_PARAMS.smoothing_factor, SPIKE_PARAMS.waveform_window);
            else
                spike_detection_results{neuron_idx} = spike_detect_SNR_sim3_SC(...
                    trace, SPIKE_PARAMS.up_threshold, SPIKE_PARAMS.down_threshold, ...
                    SPIKE_PARAMS.smoothing_factor, SPIKE_PARAMS.waveform_window);
            end
        catch ME
            warning('Spike detection failed for neuron %d: %s', neuron_idx, ME.message);
            spike_detection_results{neuron_idx} = struct('spike_idx', {{}}, 'spike_snr', {{}}, ...
                'roaster', [], 'roaster2', [], 'trace_ws', [], 'orig_trace', [], ...
                'denoise_trace', [], 'hp_trace', [], 'trace_noise', NaN);
        end
    end
    
    % Pre-allocate signal extraction arrays
    original_traces = zeros(frames_to_align, num_neurons);
    denoised_traces = zeros(frames_to_align, num_neurons);
    highpass_filtered_traces = zeros(frames_to_align, num_neurons);
    subthreshold_traces = zeros(frames_to_align, num_neurons);
    spike_raster = false(frames_to_align, num_neurons);
    spike_raster_alt = false(frames_to_align, num_neurons);
    
    % Extract all signal types from spike detection results
    for neuron_idx = 1:num_neurons
        if ~isempty(spike_detection_results{neuron_idx})
            result = spike_detection_results{neuron_idx};
            
            % Extract original trace
            if isfield(result, 'orig_trace') && ~isempty(result.orig_trace)
                original_traces(:, neuron_idx) = result.orig_trace(1, :)';
            else
                original_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            end
            
            % Extract denoised trace
            if isfield(result, 'denoise_trace') && ~isempty(result.denoise_trace)
                denoised_traces(:, neuron_idx) = result.denoise_trace(1, :)';
            else
                denoised_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            end
            
            % Extract high-pass filtered trace
            if isfield(result, 'hp_trace') && ~isempty(result.hp_trace)
                highpass_filtered_traces(:, neuron_idx) = result.hp_trace(1, :)';
            end
            
            % Extract subthreshold trace (with spikes removed)
            if isfield(result, 'trace_ws') && ~isempty(result.trace_ws)
                subthreshold_traces(:, neuron_idx) = result.trace_ws(1, :)';
            else
                subthreshold_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            end
            
            % Extract spike rasters
            if isfield(result, 'roaster') && ~isempty(result.roaster)
                spike_raster(:, neuron_idx) = result.roaster(1, :)';
            end
            if isfield(result, 'roaster2') && ~isempty(result.roaster2)
                spike_raster_alt(:, neuron_idx) = result.roaster2(1, :)';
            end
        else
            % Fallback if detection returned empty
            original_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            denoised_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            subthreshold_traces(:, neuron_idx) = fluorescence_exp_corrected(:, neuron_idx);
            warning('Spike detection returned empty for neuron %d', neuron_idx);
        end
    end
    
    % Calculate detailed spike statistics
    recording_duration = frames_to_align / IMAGING_FS;
    total_spikes = sum(spike_raster(:));
    firing_rates = sum(spike_raster, 1)' / recording_duration;
    neuron_noise_levels = zeros(num_neurons, 1);
    mean_spike_snr = zeros(num_neurons, 1);
    mean_spike_amplitude = zeros(num_neurons, 1);
    
    for neuron_idx = 1:num_neurons
        result = spike_detection_results{neuron_idx};
        if ~isempty(result)
            % Extract noise level
            if isfield(result, 'trace_noise') && ~isempty(result.trace_noise)
                neuron_noise_levels(neuron_idx) = result.trace_noise(1);
            end
            
            % Extract spike SNR
            if isfield(result, 'spike_snr') && ~isempty(result.spike_snr) && ...
               iscell(result.spike_snr) && ~isempty(result.spike_snr{1})
                mean_spike_snr(neuron_idx) = mean(result.spike_snr{1});
            end
            
            % Extract spike amplitude
            if isfield(result, 'spike_amplitude') && ~isempty(result.spike_amplitude) && ...
               iscell(result.spike_amplitude) && ~isempty(result.spike_amplitude{1})
                mean_spike_amplitude(neuron_idx) = mean(result.spike_amplitude{1});
            end
        end
    end
    
    % Print spike detection summary
    fprintf('\nSpike detection summary:\n');
    fprintf('  Total spikes detected: %d\n', total_spikes);
    fprintf('  Average spikes per neuron: %.1f\n', total_spikes/num_neurons);
    fprintf('  Average firing rate: %.2f Hz\n', mean(firing_rates));
    if any(mean_spike_snr > 0)
        fprintf('  Mean SNR across neurons: %.2f\n', mean(mean_spike_snr(mean_spike_snr > 0)));
    end
    fprintf('  Recording duration: %.1f seconds\n', recording_duration);
    
    %% ========================================================================
    %  STORE TRIAL RESULTS (Comprehensive - matching fiber pipeline)
    %% ========================================================================
    
    trial_result = struct();
    
    %% METADATA
    trial_result.metadata = METADATA;
    trial_result.metadata.trial_number = trial_numbers(trial_idx);
    trial_result.metadata.trial_folder = trial_folders{trial_idx};
    trial_result.metadata.trial_index = trial_idx;
    trial_result.metadata.total_trials = num_trials;
    
    %% PARAMETERS
    trial_result.parameters.imaging_fs = IMAGING_FS;
    trial_result.parameters.ephys_fs = EPHYS_FS;
    trial_result.parameters.num_frames = frames_to_align;
    trial_result.parameters.num_neurons = num_neurons;
    trial_result.parameters.recording_duration_sec = recording_duration;
    trial_result.parameters.motion_correction_enabled = MOTION_CORRECTION_ENABLED;
    trial_result.parameters.motion_correction_type = MOTION_CORRECTION_TYPE;
    trial_result.parameters.photobleaching_correction_enabled = APPLY_PHOTOBLEACHING_CORRECTION;
    trial_result.parameters.spike_detection_params = SPIKE_PARAMS;
    
    %% TIME VECTORS AND STIMULUS TIMING
    trial_result.time.time_vector = time_vector;
    trial_result.time.sampling_rate = IMAGING_FS;
    trial_result.time.stimulus_onset_frame = stimulus_onset_frame;
    if exist('stim_offset', 'var')
        trial_result.time.stimulus_offset_frame = stim_offset;
    end
    
    %% TIME PERIODS (matching fiber pipeline structure)
    trial_result.time_periods.stimulus_onset_frame = stimulus_onset_frame;
    if exist('stim_offset', 'var')
        trial_result.time_periods.stimulus_offset_frame = stim_offset;
        trial_result.time_periods.stim_duration_sec = (stim_offset - stimulus_onset_frame) / IMAGING_FS;
    end
    if ~isempty(stimulus_onset_frame)
        trial_result.time_periods.pre_stim_period = [0, (stimulus_onset_frame-1)/IMAGING_FS];
        trial_result.time_periods.stim_period = [stimulus_onset_frame/IMAGING_FS, ...
            (stimulus_onset_frame + DBS_DURATION_SEC*IMAGING_FS)/IMAGING_FS];
    else
        % For trials without stimulus, divide into early/middle/late
        total_time = frames_to_align / IMAGING_FS;
        trial_result.time_periods.early_period = [0, total_time/3];
        trial_result.time_periods.middle_period = [total_time/3, 2*total_time/3];
        trial_result.time_periods.late_period = [2*total_time/3, total_time];
    end
    
    %% SIGNALS - COMPREHENSIVE FLUORESCENCE TRACE STORAGE
    % Store all processing stages for complete traceability (matching fiber pipeline)
    
    % Raw and basic processing
    trial_result.signals.fluorescence_raw = neuron_traces_aligned;           % Raw from ROI extraction
    trial_result.signals.fluorescence_detrended = fluorescence_detrended;    % After linear detrending
    trial_result.signals.fluorescence_corrected = fluorescence_exp_corrected; % After exp correction (ΔF/F)
    trial_result.signals.fluorescence_zscored = fluorescence_zscored;        % Z-scored version
    trial_result.signals.final_processed_traces = fluorescence_exp_corrected; % Final output (ΔF/F)
    trial_result.signals.background = background_trace(1:frames_to_align);
    
    % From spike detection
    trial_result.signals.original_traces = original_traces;
    trial_result.signals.denoised_traces = denoised_traces;
    trial_result.signals.highpass_filtered = highpass_filtered_traces;
    trial_result.signals.subthreshold = subthreshold_traces;
    
    % ΔF/F metadata (matching fiber pipeline)
    trial_result.signals.deltaF_F_method = 'double_exponential_correction';
    if exist('deltaF_F_baseline_window', 'var') && ~isempty(deltaF_F_baseline_window)
        trial_result.signals.deltaF_F_baseline_window = deltaF_F_baseline_window;
        trial_result.signals.deltaF_F_baseline_time = deltaF_F_baseline_time;
    end
    if exist('F0_values', 'var')
        trial_result.signals.F0_values = F0_values;  % F0 value for each neuron
    end
    if exist('exp_fit_params', 'var')
        trial_result.signals.exp_fit_params = exp_fit_params;  % Exponential fit parameters
    end
    
    %% EPHYS SIGNALS
    if LOAD_EPHYS_DATA && ~isempty(lfp_aligned)
        trial_result.ephys.lfp_raw_aligned = lfp_aligned;
        trial_result.ephys.lfp_zscored = zscore(lfp_aligned);
        trial_result.ephys.lfp_sampling_rate = IMAGING_FS;  % After downsampling
        trial_result.ephys.lfp_original_sampling_rate = EPHYS_FS;
        trial_result.ephys.stimulus_aligned = stimulus_aligned;
        trial_result.ephys.channel_definitions = CHANNELS;
        
        % Store camera trigger indices for alignment verification
        if exist('aligned_camera_indices', 'var')
            trial_result.ephys.camera_trigger_indices = aligned_camera_indices;
        end
        
        % Frame rate statistics
        if exist('computed_frame_rate', 'var') && exist('frame_rate_variability', 'var')
            trial_result.frame_rates.computed_hz = computed_frame_rate;
            trial_result.frame_rates.variability_percent = frame_rate_variability;
            trial_result.frame_rates.used_manual = (computed_frame_rate < MIN_ACCEPTABLE_FRAME_RATE);
        end
    end
    
    %% MOTION CORRECTION
    trial_result.motion_correction.enabled = MOTION_CORRECTION_ENABLED;
    trial_result.motion_correction.roi_position = shared_motion_roi_position;
    if exist('trial_motion_shifts', 'var') && ~isempty(trial_motion_shifts)
        % Extract shifts as array for easier analysis
        shifts_array = zeros(length(trial_motion_shifts), 2);
        for s_idx = 1:length(trial_motion_shifts)
            if isstruct(trial_motion_shifts(s_idx)) && isfield(trial_motion_shifts(s_idx), 'shifts')
                shifts_array(s_idx, :) = trial_motion_shifts(s_idx).shifts;
            end
        end
        trial_result.motion_correction.shifts = shifts_array;
        trial_result.motion_correction.max_shift_x = max(abs(shifts_array(:,1)));
        trial_result.motion_correction.max_shift_y = max(abs(shifts_array(:,2)));
        trial_result.motion_correction.mean_shift = mean(sqrt(sum(shifts_array.^2, 2)));
    end
    
    %% ROI INFO
    trial_result.rois.masks = shared_neuron_roi_masks;
    trial_result.rois.centroids = shared_neuron_centroids;
    trial_result.rois.num_neurons = num_neurons;
    if trial_idx == 1
        trial_result.rois.reference_frame = roi_reference_frame;
    end
    
    %% SPIKES - Detection results and statistics
    trial_result.spikes.detection_results = spike_detection_results;
    trial_result.spikes.spike_raster = spike_raster;
    trial_result.spikes.spike_raster_alt = spike_raster_alt;
    trial_result.spikes.firing_rates_hz = firing_rates;
    trial_result.spikes.total_spikes = total_spikes;
    trial_result.spikes.mean_firing_rate = mean(firing_rates);
    trial_result.spikes.neuron_noise_levels = neuron_noise_levels;
    trial_result.spikes.mean_spike_snr = mean_spike_snr;
    trial_result.spikes.mean_spike_amplitude = mean_spike_amplitude;
    trial_result.spikes.detection_type = SPIKE_DETECTION_TYPE;
    trial_result.spikes.baseline_window = baseline_window;
    
    % Per-neuron spike data (for detailed analysis)
    spike_times_per_neuron = cell(num_neurons, 1);
    spike_indices_per_neuron = cell(num_neurons, 1);
    spike_snr_per_neuron = cell(num_neurons, 1);
    spike_amplitudes_per_neuron = cell(num_neurons, 1);
    spike_waveforms_per_neuron = cell(num_neurons, 1);
    mean_waveforms_per_neuron = cell(num_neurons, 1);
    
    for n_idx = 1:num_neurons
        result = spike_detection_results{n_idx};
        if ~isempty(result) && isfield(result, 'spike_idx') && iscell(result.spike_idx) && ~isempty(result.spike_idx{1})
            valid_indices = result.spike_idx{1};
            valid_indices = valid_indices(valid_indices <= frames_to_align);
            
            spike_indices_per_neuron{n_idx} = valid_indices;
            spike_times_per_neuron{n_idx} = time_vector(valid_indices);
            
            if isfield(result, 'spike_snr') && iscell(result.spike_snr)
                spike_snr_per_neuron{n_idx} = result.spike_snr{1};
            end
            if isfield(result, 'spike_amplitude') && iscell(result.spike_amplitude)
                spike_amplitudes_per_neuron{n_idx} = result.spike_amplitude{1};
            end
            if isfield(result, 'spike_waveforms') && iscell(result.spike_waveforms)
                spike_waveforms_per_neuron{n_idx} = result.spike_waveforms{1};
            end
            if isfield(result, 'mean_waveform') && iscell(result.mean_waveform)
                mean_waveforms_per_neuron{n_idx} = result.mean_waveform{1};
            end
        end
    end
    
    trial_result.spikes.spike_times_seconds = spike_times_per_neuron;
    trial_result.spikes.spike_indices = spike_indices_per_neuron;
    trial_result.spikes.spike_snr_values = spike_snr_per_neuron;
    trial_result.spikes.spike_amplitudes = spike_amplitudes_per_neuron;
    trial_result.spikes.individual_waveforms = spike_waveforms_per_neuron;
    trial_result.spikes.mean_waveforms = mean_waveforms_per_neuron;
    
    %% CONNECTIVITY DATA - Ready-to-use data for connectivity analysis
    trial_result.connectivity_data.spike_trains = spike_raster;
    trial_result.connectivity_data.subthreshold_signals = subthreshold_traces;
    trial_result.connectivity_data.full_signals = fluorescence_exp_corrected;
    trial_result.connectivity_data.denoised_signals = denoised_traces;
    trial_result.connectivity_data.highpass_signals = highpass_filtered_traces;
    trial_result.connectivity_data.time_vector = time_vector;
    trial_result.connectivity_data.sampling_rate = IMAGING_FS;
    
    %% DBS INFO
    trial_result.dbs.frequency_hz = DBS_FREQUENCY_HZ;
    trial_result.dbs.comparison_type = DBS_COMPARISON_TYPE;
    trial_result.dbs.duration_sec = DBS_DURATION_SEC;
    trial_result.dbs.voltage = DBS_VOLTAGE;
    trial_result.dbs.pulse_width_us = DBS_PULSE_WIDTH_US;
    trial_result.dbs.current_ua = DBS_CURRENT_UA;
    
    all_trial_results{trial_idx} = trial_result;
    
    fprintf('\nTrial %d processing complete.\n', trial_idx);
end

%% ============================================================================
%  SECTION 5: SAVE RESULTS
%  ============================================================================

fprintf('\n');
fprintf('============================================================================\n');
fprintf('  SAVING RESULTS\n');
fprintf('============================================================================\n');

% Create comprehensive output structure
CellularAnalysis = struct();
CellularAnalysis.metadata = METADATA;
CellularAnalysis.metadata.num_trials = num_trials;
CellularAnalysis.metadata.num_neurons = num_neurons;
CellularAnalysis.trials = all_trial_results;

% Shared ROI info (same across trials)
CellularAnalysis.shared_rois.masks = shared_neuron_roi_masks;
CellularAnalysis.shared_rois.centroids = shared_neuron_centroids;

% Save main data file
output_filename = sprintf('%s_%s-%s_CellularAnalysis.mat', ...
    MOUSE_NAME, RECORDING_DATE, RECORDING_ID);
output_filepath = fullfile(OUTPUT_FOLDER, output_filename);

fprintf('Saving to: %s\n', output_filepath);
save(output_filepath, 'CellularAnalysis', '-v7.3');

%% ============================================================================
%  SECTION 6: VISUALIZATION
%  ============================================================================
% Generate comprehensive figures for quality control and analysis overview.
% Each trial gets its own set of figures saved to the output folder.

if GENERATE_FIGURES
    fprintf('\n');
    fprintf('============================================================================\n');
    fprintf('  GENERATING VISUALIZATION FIGURES\n');
    fprintf('============================================================================\n');
    
    % Create figures folder
    figures_folder = fullfile(OUTPUT_FOLDER, 'figures');
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end
    
    % Define consistent styling
    STYLE = struct();
    STYLE.title_fontsize = 14;
    STYLE.label_fontsize = 12;
    STYLE.axis_fontsize = 10;
    STYLE.linewidth = 1.2;
    STYLE.colors.lfp = [0.25, 0.25, 0.25];           % Dark gray for LFP
    STYLE.colors.neuron = [1, 0.60, 0.20];           % Orange for neurons
    STYLE.colors.spike = [0.6, 0, 0];                % Dark red for spikes
    STYLE.colors.stim_shade = [1, 0.85, 0.85];       % Soft red for stim period
    STYLE.colors.roi_outline = [1, 0, 0];            % Red for ROI outlines
    
    % Define frequency bands for spectral analysis
    FREQ_BANDS = {
        'Delta',  [1, 4],   [0.2, 0.4, 0.8];
        'Theta',  [4, 8],   [0.3, 0.6, 0.3];
        'Alpha',  [8, 13],  [0.8, 0.6, 0.2];
        'Beta',   [13, 30], [0.6, 0.3, 0.6];
        'Gamma',  [30, 80], [0.8, 0.2, 0.2];
    };
    
    %% ------------------------------------------------------------------------
    %  FIGURE 1: ROI VISUALIZATION
    %  ------------------------------------------------------------------------
    fprintf('  Generating ROI visualization...\n');
    
    fig_roi = figure('Name', 'Neuron ROIs', 'Color', 'w', 'Position', [100, 100, 900, 700]);
    
    imagesc(roi_reference_frame);
    axis image; axis off; colormap(gray);
    colorbar('FontSize', STYLE.axis_fontsize);
    hold on;
    
    % Draw ROI outlines and labels
    for n = 1:num_neurons
        [boundaries, ~] = bwboundaries(shared_neuron_roi_masks{n}, 'noholes');
        if ~isempty(boundaries)
            b = boundaries{1};
            plot(b(:,2), b(:,1), 'Color', STYLE.colors.roi_outline, 'LineWidth', 2);
            text(shared_neuron_centroids(n,1), shared_neuron_centroids(n,2), ...
                sprintf('%d', n), 'Color', STYLE.colors.roi_outline, ...
                'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    
    title(sprintf('%s %s-%s: Neuron ROIs (n=%d)', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, num_neurons), ...
          'FontSize', STYLE.title_fontsize);
    
    if SAVE_FIGURES
        save_figure(fig_roi, figures_folder, sprintf('%s_%s-%s_01_ROIs', MOUSE_NAME, RECORDING_DATE, RECORDING_ID));
    end
    
    %% ------------------------------------------------------------------------
    %  FIGURES 2-N: PER-TRIAL VISUALIZATIONS
    %  ------------------------------------------------------------------------
    
    for trial_idx = 1:num_trials
        fprintf('  Generating figures for Trial %d/%d...\n', trial_idx, num_trials);
        
        % Get trial data - skip if empty (trial was skipped during processing)
        trial_data = all_trial_results{trial_idx};
        if isempty(trial_data) || ~isstruct(trial_data)
            fprintf('    Trial %d was skipped (no data). Skipping visualization.\n', trial_idx);
            continue;
        end
        
        trial_name = sprintf('Trial%d', trial_data.metadata.trial_number);
        
        % Extract key signals
        time_vec = trial_data.time.time_vector;
        fluor_corrected = trial_data.signals.fluorescence_corrected;
        fluor_zscored = trial_data.signals.fluorescence_zscored;
        spike_raster_trial = trial_data.spikes.spike_raster;
        stim_onset = trial_data.time.stimulus_onset_frame;
        fs = trial_data.parameters.imaging_fs;
        
        % Get LFP if available
        if isfield(trial_data, 'ephys') && isfield(trial_data.ephys, 'lfp_zscored')
            lfp_z = trial_data.ephys.lfp_zscored;
            has_lfp = true;
        else
            has_lfp = false;
        end
        
        %% FIGURE 2: Combined LFP and Voltage Traces with Spike Rasters
        fig_traces = figure('Name', sprintf('%s - Neural Activity', trial_name), ...
                           'Color', 'w', 'Position', [100, 100, 1200, 800]);
        
        % Select neurons to display (up to 8, or all if fewer)
        max_display = min(8, num_neurons);
        selected_neurons = 1:max_display;
        trace_spacing = 6;  % Vertical spacing between traces
        
        axes('Position', [0.08, 0.1, 0.82, 0.85]);
        hold on;
        
        % Calculate offsets
        num_traces = length(selected_neurons);
        if has_lfp
            lfp_offset = num_traces * trace_spacing + 4;
        end
        
        % Plot LFP at top
        if has_lfp
            plot(time_vec, lfp_z + lfp_offset, 'Color', STYLE.colors.lfp, ...
                 'LineWidth', STYLE.linewidth);
            text(max(time_vec)*1.01, lfp_offset, 'LFP', 'Color', STYLE.colors.lfp, ...
                 'FontSize', STYLE.axis_fontsize, 'FontWeight', 'bold');
        end
        
        % Plot neuron traces
        for i = 1:length(selected_neurons)
            neuron_id = selected_neurons(i);
            offset = (i-1) * trace_spacing;
            
            % Use z-scored traces for display
            trace_display = fluor_zscored(:, neuron_id);
            
            plot(time_vec, trace_display + offset, 'Color', STYLE.colors.neuron, ...
                 'LineWidth', STYLE.linewidth);
            
            % Add neuron label
            text(max(time_vec)*1.01, offset, sprintf('N%d', neuron_id), ...
                 'Color', 'k', 'FontSize', STYLE.axis_fontsize-1);
        end
        
        % Plot spike markers
        for i = 1:length(selected_neurons)
            neuron_id = selected_neurons(i);
            offset = (i-1) * trace_spacing + trace_spacing * 0.7;
            
            spike_indices = find(spike_raster_trial(:, neuron_id));
            if ~isempty(spike_indices)
                spike_times = time_vec(spike_indices);
                for s = 1:length(spike_times)
                    line([spike_times(s), spike_times(s)], offset + [0, 1], ...
                         'Color', STYLE.colors.spike, 'LineWidth', 0.5);
                end
            end
        end
        
        % Add stimulus onset line
        if ~isempty(stim_onset)
            xline(0, '--', 'Color', [0.4, 0.4, 0.4], 'LineWidth', 1.5);
            
            % Add stimulus period shading
            stim_end_time = DBS_DURATION_SEC;
            yl = ylim;
            patch([0, stim_end_time, stim_end_time, 0], [yl(1), yl(1), yl(2), yl(2)], ...
                  STYLE.colors.stim_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.4);
            uistack(gca, 'bottom');
        end
        
        % Formatting
        title(sprintf('%s %s-%s %s: Neural Activity', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name), ...
              'FontSize', STYLE.title_fontsize);
        xlabel('Time from stimulus onset (s)', 'FontSize', STYLE.label_fontsize);
        ylabel('Z-scored signal', 'FontSize', STYLE.label_fontsize);
        set(gca, 'FontSize', STYLE.axis_fontsize, 'YTick', []);
        xlim([min(time_vec), max(time_vec)*1.1]);
        
        if SAVE_FIGURES
            save_figure(fig_traces, figures_folder, ...
                sprintf('%s_%s-%s_%s_02_traces', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name));
        end
        
        %% FIGURE 3: Spike Waveforms
        spike_results = trial_data.spikes.detection_results;
        
        % Count neurons with detected waveforms
        neurons_with_spikes = sum(cellfun(@(x) ~isempty(x) && isfield(x, 'mean_waveform') && ...
                                  ~isempty(x.mean_waveform) && iscell(x.mean_waveform) && ...
                                  ~isempty(x.mean_waveform{1}), spike_results));
        
        if neurons_with_spikes > 0
            fig_waveforms = figure('Name', sprintf('%s - Spike Waveforms', trial_name), ...
                                  'Color', 'w', 'Position', [100, 100, 1000, 800]);
            
            % Calculate subplot layout
            n_rows = ceil(sqrt(neurons_with_spikes));
            n_cols = ceil(neurons_with_spikes / n_rows);
            
            % Calculate waveform time vector (ms)
            waveform_window = SPIKE_PARAMS.waveform_window;
            waveform_time_ms = 1000 * (waveform_window(1):waveform_window(2)) / fs;
            
            subplot_idx = 0;
            for neuron_idx = 1:num_neurons
                result = spike_results{neuron_idx};
                if isempty(result) || ~isfield(result, 'mean_waveform') || ...
                   ~iscell(result.mean_waveform) || isempty(result.mean_waveform{1})
                    continue;
                end
                
                subplot_idx = subplot_idx + 1;
                subplot(n_rows, n_cols, subplot_idx);
                hold on;
                
                % Plot individual waveforms (semi-transparent)
                if isfield(result, 'spike_waveforms') && iscell(result.spike_waveforms) && ...
                   ~isempty(result.spike_waveforms{1})
                    waveforms = result.spike_waveforms{1};
                    n_to_plot = min(size(waveforms, 1), 50);
                    for w = 1:n_to_plot
                        plot(waveform_time_ms, waveforms(w, :), 'Color', [0.7, 0.7, 0.7, 0.3]);
                    end
                end
                
                % Plot mean waveform
                mean_wf = result.mean_waveform{1};
                plot(waveform_time_ms, mean_wf, 'k-', 'LineWidth', 2);
                
                % Add spike peak marker
                xline(0, 'r--', 'LineWidth', 1);
                
                % Format subplot
                n_spikes = sum(spike_raster_trial(:, neuron_idx));
                title(sprintf('N%d (n=%d)', neuron_idx, n_spikes), 'FontSize', STYLE.axis_fontsize);
                
                if subplot_idx > (n_rows-1) * n_cols
                    xlabel('Time (ms)', 'FontSize', STYLE.axis_fontsize-1);
                end
                if mod(subplot_idx-1, n_cols) == 0
                    ylabel('ΔF/F', 'FontSize', STYLE.axis_fontsize-1);
                end
                
                xlim([waveform_time_ms(1), waveform_time_ms(end)]);
                set(gca, 'FontSize', STYLE.axis_fontsize-1);
                grid on;
            end
            
            sgtitle(sprintf('%s %s-%s %s: Spike Waveforms', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name), ...
                    'FontSize', STYLE.title_fontsize);
            
            if SAVE_FIGURES
                save_figure(fig_waveforms, figures_folder, ...
                    sprintf('%s_%s-%s_%s_03_waveforms', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name));
            end
        end
        
        %% FIGURE 4: Spectral Analysis
        if has_lfp
            fig_spectral = figure('Name', sprintf('%s - Spectral Analysis', trial_name), ...
                                 'Color', 'w', 'Position', [100, 100, 1400, 900]);
            
            % Spectrogram parameters
            spec_window = min(380, round(length(time_vec) / 5));
            spec_overlap = spec_window - 1;
            freq_vec = 1:min(100, fs/2);
            
            % Compute spectrograms
            [S_lfp, F_lfp, T_lfp] = spectrogram(trial_data.ephys.lfp_raw_aligned, ...
                spec_window, spec_overlap, freq_vec, fs);
            
            % Use first neuron for voltage imaging spectrogram
            [S_neuron, F_neuron, T_neuron] = spectrogram(fluor_corrected(:,1), ...
                spec_window, spec_overlap, freq_vec, fs);
            
            % Align time to stimulus onset
            if ~isempty(stim_onset)
                T_lfp_aligned = T_lfp - stim_onset/fs;
                T_neuron_aligned = T_neuron - stim_onset/fs;
            else
                T_lfp_aligned = T_lfp;
                T_neuron_aligned = T_neuron;
            end
            
            % Subplot 1: LFP Spectrogram
            subplot(2, 3, 1);
            imagesc(T_lfp_aligned, F_lfp, abs(S_lfp) .* repmat(sqrt(F_lfp), 1, size(S_lfp, 2)));
            axis xy; colormap(gca, hot); colorbar;
            title('LFP Spectrogram', 'FontSize', STYLE.label_fontsize);
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            if ~isempty(stim_onset), hold on; xline(0, 'w--', 'LineWidth', 1.5); end
            set(gca, 'FontSize', STYLE.axis_fontsize);
            
            % Subplot 2: Neuron Spectrogram
            subplot(2, 3, 2);
            imagesc(T_neuron_aligned, F_neuron, abs(S_neuron) .* repmat(sqrt(F_neuron), 1, size(S_neuron, 2)));
            axis xy; colormap(gca, parula); colorbar;
            title('Voltage Imaging Spectrogram (N1)', 'FontSize', STYLE.label_fontsize);
            xlabel('Time (s)'); ylabel('Frequency (Hz)');
            if ~isempty(stim_onset), hold on; xline(0, 'w--', 'LineWidth', 1.5); end
            set(gca, 'FontSize', STYLE.axis_fontsize);
            
            % Subplot 3: LFP Band Power
            subplot(2, 3, 3);
            hold on;
            for b = 1:size(FREQ_BANDS, 1)
                band_range = FREQ_BANDS{b, 2};
                band_color = FREQ_BANDS{b, 3};
                freq_idx = F_lfp >= band_range(1) & F_lfp <= band_range(2);
                band_power = mean(abs(S_lfp(freq_idx, :)), 1);
                plot(T_lfp_aligned, band_power, 'Color', band_color, 'LineWidth', 1.5);
            end
            if ~isempty(stim_onset), xline(0, 'r--', 'LineWidth', 1.5); end
            title('LFP Band Power', 'FontSize', STYLE.label_fontsize);
            xlabel('Time (s)'); ylabel('Power');
            legend(FREQ_BANDS(:,1), 'Location', 'best', 'FontSize', 8);
            set(gca, 'FontSize', STYLE.axis_fontsize);
            grid on;
            
            % Subplot 4: Neuron Band Power
            subplot(2, 3, 4);
            hold on;
            for b = 1:size(FREQ_BANDS, 1)
                band_range = FREQ_BANDS{b, 2};
                band_color = FREQ_BANDS{b, 3};
                freq_idx = F_neuron >= band_range(1) & F_neuron <= band_range(2);
                band_power = mean(abs(S_neuron(freq_idx, :)), 1);
                plot(T_neuron_aligned, band_power, 'Color', band_color, 'LineWidth', 1.5);
            end
            if ~isempty(stim_onset), xline(0, 'r--', 'LineWidth', 1.5); end
            title('Voltage Band Power (N1)', 'FontSize', STYLE.label_fontsize);
            xlabel('Time (s)'); ylabel('Power');
            legend(FREQ_BANDS(:,1), 'Location', 'best', 'FontSize', 8);
            set(gca, 'FontSize', STYLE.axis_fontsize);
            grid on;
            
            % Subplot 5-6: Power Spectral Density
            subplot(2, 3, [5, 6]);
            lfp_psd = mean(abs(S_lfp), 2);
            neuron_psd = mean(abs(S_neuron), 2);
            
            yyaxis left;
            plot(F_lfp, lfp_psd, 'k-', 'LineWidth', 1.5);
            ylabel('LFP Power', 'Color', 'k');
            set(gca, 'YColor', 'k');
            
            yyaxis right;
            plot(F_neuron, neuron_psd, 'Color', STYLE.colors.neuron, 'LineWidth', 1.5);
            ylabel('Voltage Imaging Power', 'Color', STYLE.colors.neuron);
            set(gca, 'YColor', STYLE.colors.neuron);
            
            title('Power Spectral Density', 'FontSize', STYLE.label_fontsize);
            xlabel('Frequency (Hz)');
            legend('LFP', 'Voltage Imaging', 'Location', 'best');
            set(gca, 'FontSize', STYLE.axis_fontsize);
            grid on;
            
            sgtitle(sprintf('%s %s-%s %s: Spectral Analysis', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name), ...
                    'FontSize', STYLE.title_fontsize);
            
            if SAVE_FIGURES
                save_figure(fig_spectral, figures_folder, ...
                    sprintf('%s_%s-%s_%s_04_spectral', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name));
            end
        end
        
        %% FIGURE 5: Spike Statistics Summary
        fig_stats = figure('Name', sprintf('%s - Spike Statistics', trial_name), ...
                          'Color', 'w', 'Position', [100, 100, 1000, 600]);
        
        % Subplot 1: Firing rates bar plot
        subplot(2, 2, 1);
        bar(trial_data.spikes.firing_rates_hz, 'FaceColor', STYLE.colors.neuron);
        xlabel('Neuron ID', 'FontSize', STYLE.label_fontsize);
        ylabel('Firing Rate (Hz)', 'FontSize', STYLE.label_fontsize);
        title('Firing Rates', 'FontSize', STYLE.label_fontsize);
        set(gca, 'FontSize', STYLE.axis_fontsize);
        
        % Subplot 2: Spike SNR
        subplot(2, 2, 2);
        snr_data = trial_data.spikes.mean_spike_snr;
        snr_data(snr_data == 0) = NaN;
        bar(snr_data, 'FaceColor', [0.3, 0.6, 0.3]);
        xlabel('Neuron ID', 'FontSize', STYLE.label_fontsize);
        ylabel('Mean Spike SNR', 'FontSize', STYLE.label_fontsize);
        title('Spike Quality (SNR)', 'FontSize', STYLE.label_fontsize);
        set(gca, 'FontSize', STYLE.axis_fontsize);
        
        % Subplot 3: Raster plot
        subplot(2, 2, [3, 4]);
        hold on;
        for n = 1:num_neurons
            spike_times = time_vec(spike_raster_trial(:, n));
            if ~isempty(spike_times)
                plot(spike_times, n * ones(size(spike_times)), 'k.', 'MarkerSize', 2);
            end
        end
        if ~isempty(stim_onset)
            xline(0, 'r--', 'LineWidth', 1.5);
            patch([0, DBS_DURATION_SEC, DBS_DURATION_SEC, 0], [0, 0, num_neurons+1, num_neurons+1], ...
                  STYLE.colors.stim_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end
        xlabel('Time (s)', 'FontSize', STYLE.label_fontsize);
        ylabel('Neuron ID', 'FontSize', STYLE.label_fontsize);
        title('Population Spike Raster', 'FontSize', STYLE.label_fontsize);
        ylim([0, num_neurons+1]);
        xlim([min(time_vec), max(time_vec)]);
        set(gca, 'FontSize', STYLE.axis_fontsize);
        
        sgtitle(sprintf('%s %s-%s %s: Spike Statistics (Total: %d spikes, Mean rate: %.2f Hz)', ...
                MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name, ...
                trial_data.spikes.total_spikes, trial_data.spikes.mean_firing_rate), ...
                'FontSize', STYLE.title_fontsize);
        
        if SAVE_FIGURES
            save_figure(fig_stats, figures_folder, ...
                sprintf('%s_%s-%s_%s_05_spikestats', MOUSE_NAME, RECORDING_DATE, RECORDING_ID, trial_name));
        end
        
    end  % End trial loop for figures
    
    fprintf('  Figure generation complete.\n');
    
    if ~SAVE_FIGURES
        fprintf('  NOTE: Figures displayed but not saved (SAVE_FIGURES = false)\n');
    end
end

fprintf('\n');
fprintf('============================================================================\n');
fprintf('  PROCESSING COMPLETE\n');
fprintf('============================================================================\n');
fprintf('  Mouse: %s\n', MOUSE_NAME);
fprintf('  Recording: %s-%s\n', RECORDING_DATE, RECORDING_ID);
fprintf('  Trials processed: %d\n', num_trials);
fprintf('  Neurons: %d\n', num_neurons);
fprintf('  DBS: %d Hz (%s)\n', DBS_FREQUENCY_HZ, DBS_COMPARISON_TYPE);
fprintf('  Output: %s\n', output_filepath);
fprintf('============================================================================\n');

%% ============================================================================
%  LOCAL HELPER FUNCTIONS
%  ============================================================================

function save_figure(fig_handle, folder, filename_base)
    % SAVE_FIGURE - Save figure in multiple formats (PNG, FIG)
    %
    % Saves the figure to the specified folder with the given base filename.
    % Creates both a PNG (for quick viewing) and a MATLAB FIG file (for editing).
    %
    % Args:
    %   fig_handle: Handle to the figure to save
    %   folder: Directory path to save the figure
    %   filename_base: Base filename without extension
    
    try
        % Save as PNG (high resolution for publications)
        png_path = fullfile(folder, [filename_base, '.png']);
        saveas(fig_handle, png_path);
        
        % Save as MATLAB figure (for later editing)
        fig_path = fullfile(folder, [filename_base, '.fig']);
        saveas(fig_handle, fig_path);
        
        % Close figure after saving
        close(fig_handle);
        
    catch ME
        warning('Failed to save figure %s: %s', filename_base, ME.message);
    end
end
