%% ============================================================================
%  COMPREHENSIVE FIBER PHOTOMETRY AND LFP ANALYSIS (ENHANCED WITH MULTI-FIBER)
%  ============================================================================
%  Combines voltage imaging (fiber photometry) with local field potential (LFP)
%  recordings from Open Ephys. Includes photobleaching correction, motion
%  correction, time-frequency analysis, and phase-locking analysis.
%
%  ENHANCED FEATURES:
%  - Original single fiber analysis preserved
%  - Added multi-fiber combined visualization
%  - Individual fiber plots remain unchanged
%  - Enhanced multi-fiber trace comparison
%  - All original functionality maintained
%  ============================================================================

close all; clear; clc;

%% Add required paths
% External toolboxes via the centralised config (config/lab_paths.m).
% Override per machine via config/paths_local.m. See config/README.md.
addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'config'));
setup_lab_paths();
% This pipeline's own helpers (core/, core/utils/) ON TOP of the path, so our
% smooth2a/fastsmooth/viridis/turbo are used rather than any external or
% built-in versions.
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), 'core')));

%% ============================================================================
%  SECTION 1: CONFIGURATION PARAMETERS
%  ============================================================================

%% RUN CONFIGURATION
% All run parameters (mouse/date/paths/options/bands/plots/...) now live in
% config/fiber_preprocessing_singletrial_config.m. Edit THAT file for your
% recording. It is loaded here as a script-include, so every variable it sets
% lands in this workspace with exactly the same name as before.
run(fullfile(fileparts(mfilename('fullpath')), 'config', 'fiber_preprocessing_singletrial_config.m'));

%% METADATA
METADATA = struct();
METADATA.mouse_name = MOUSE_NAME;
METADATA.experimenter = EXPERIMENTER;
METADATA.recording_date = RECORDING_DATE;
METADATA.recording_id = RECORDING_ID;

fprintf('Analysis metadata:\n');
fprintf('  Mouse: %s\n', METADATA.mouse_name);
fprintf('  Date: %s\n', METADATA.recording_date);
fprintf('  Recording ID: %s\n', METADATA.recording_id);
fprintf('  Experimenter: %s\n', METADATA.experimenter);

%% AUTOMATIC TRIAL DETECTION
if strcmp(ANALYSIS_MODE, 'single_trial')
    PROCESS_SINGLE_TRIAL = true;

    % Verify base folder exists
    if ~exist(base_folder, 'dir')
        error('Base folder does not exist: %s\nPlease check MOUSE_NAME, RECORDING_DATE, and RECORDING_ID parameters.', base_folder);
    end

    % Automatically find trial folder inside base_folder
    fprintf('Searching for trial folders in: %s\n', base_folder);

    % Look for folders (excluding . and ..)
    folder_contents = dir(base_folder);
    folder_contents = folder_contents([folder_contents.isdir]);  % Keep only directories
    folder_contents = folder_contents(~ismember({folder_contents.name}, {'.', '..'}));  % Remove . and ..

    if isempty(folder_contents)
        error('No trial folders found in: %s', base_folder);
    end

    % If multiple folders found, list them and use the first one
    if length(folder_contents) > 1
        fprintf('Multiple folders found:\n');
        for i = 1:length(folder_contents)
            fprintf('  %d. %s\n', i, folder_contents(i).name);
        end

        % Automatically use the first folder
        trial_name = folder_contents(1).name;
        fprintf('Using first folder: %s\n', trial_name);
    else
        % Only one folder found - use it automatically
        trial_name = folder_contents(1).name;
        fprintf('Found trial folder: %s\n', trial_name);
    end

    % Construct full trial folder path
    trial_folder = fullfile(base_folder, trial_name);
    fprintf('Processing trial: %s\n', trial_folder);

else
    % Multi-FOV mode: process all fov folders in base directory
    PROCESS_SINGLE_TRIAL = false;

    if ~exist(base_folder, 'dir')
        error('Base folder does not exist: %s', base_folder);
    end

    fprintf('Multi-FOV mode: Will process all fov* folders in %s\n', base_folder);
end

fprintf('=== COMPREHENSIVE FIBER PHOTOMETRY & LFP ANALYSIS (ENHANCED) ===\n');
fprintf('Analysis Mode: %s\n', ANALYSIS_MODE);
if PROCESS_SINGLE_TRIAL
    fprintf('Processing Trial: %s\n', trial_name);
end
fprintf('Motion Correction: %s\n', char(string(logical(MOTION_CORRECTION))));
fprintf('Photobleaching Correction: %s\n', char(string(APPLY_PHOTOBLEACHING_CORRECTION)));
fprintf('Multi-Fiber Plot: %s\n', char(string(GENERATE_MULTI_FIBER_PLOT)));
fprintf('Load Open Ephys Data: %s\n', char(string(LOAD_EPHYS_DATA)));
fprintf('\n');

%% ============================================================================
%  SECTION 2: LOAD AND PROCESS IMAGING DATA
%  ============================================================================

fprintf('=== LOADING FIBER PHOTOMETRY DATA ===\n');

if PROCESS_SINGLE_TRIAL
    %% Single trial processing
    trial_folder = fullfile(base_folder, trial_name);
    if ~exist(trial_folder, 'dir')
        error('Trial folder not found: %s', trial_folder);
    end

    % % Find TIFF file
    % ome_file = dir(fullfile(trial_folder, '*.ome.tif'));
    % if isempty(ome_file)
    %     tif_file = dir(fullfile(trial_folder, '*.tif'));
    %     if isempty(tif_file)
    %         error('No TIFF file found in %s', trial_folder);
    %     end
    %     ome_file = tif_file;
    % end
    % 
    % fullpath = fullfile(trial_folder, ome_file(1).name);
    % fprintf('Loading: %s\n', ome_file(1).name);
    % 
    % % Load TIFF stack
    % info = imfinfo(fullpath);
    % numFrames = numel(info);
    % fprintf('  Total frames: %d\n', numFrames);
    % 
    % t = Tiff(fullpath, 'r');
    % stack = zeros(info(1).Height, info(1).Width, numFrames, 'uint16');
    % for k = 1:numFrames
    %     if mod(k, 500) == 0
    %         fprintf('  Loading frame %d/%d\n', k, numFrames);
    %     end
    %     t.setDirectory(k);
    %     stack(:,:,k) = t.read();
    % end
    % t.close();

    % Find TIFF files (handle multiple OME-TIFF parts)
    ome_files = dir(fullfile(trial_folder, '*.ome.tif'));
    if isempty(ome_files)
        tif_files = dir(fullfile(trial_folder, '*.tif'));
        if isempty(tif_files)
            error('No TIFF file found in %s', trial_folder);
        end
        ome_files = tif_files;
    end

    % Sort files by name to ensure correct order (_1, _2, etc.)
    [~, idx] = sort({ome_files.name});
    ome_files = ome_files(idx);

    % Display found files
    fprintf('Found %d TIFF file(s):\n', length(ome_files));
    for i = 1:length(ome_files)
        fprintf('  %s\n', ome_files(i).name);
    end

    % Count total number of frames across all files
    numTotalFrames = 0;
    frameCounts = zeros(length(ome_files), 1);
    firstInfo = [];

    for i = 1:length(ome_files)
        fullpath = fullfile(trial_folder, ome_files(i).name);
        info = imfinfo(fullpath);
        frameCounts(i) = numel(info);
        numTotalFrames = numTotalFrames + frameCounts(i);

        if i == 1
            firstInfo = info(1);
            fprintf('  Frame dimensions: %d x %d\n', firstInfo.Height, firstInfo.Width);
        end
    end

    fprintf('  Total frames across all files: %d\n', numTotalFrames);

    % Preallocate full stack
    stack = zeros(firstInfo.Height, firstInfo.Width, numTotalFrames, 'uint16');

    % Load all TIFF parts into a single concatenated stack
    frameIdx = 1;
    for i = 1:length(ome_files)
        fullpath = fullfile(trial_folder, ome_files(i).name);
        fprintf('Loading file %d/%d: %s (%d frames)\n', i, length(ome_files), ome_files(i).name, frameCounts(i));

        t = Tiff(fullpath, 'r');
        for k = 1:frameCounts(i)
            if mod(frameIdx, 1000) == 0
                fprintf('  Loading frame %d/%d\n', frameIdx, numTotalFrames);
            end
            t.setDirectory(k);
            stack(:,:,frameIdx) = t.read();
            frameIdx = frameIdx + 1;
        end
        t.close();
    end

    % Motion correction ROI selection
    if MOTION_CORRECTION
        disp('Select ROI for motion correction.');
        % avgFrame = mean(stack(:,:,100:min(190, numFrames)), 3);
        avgFrame = mean(stack(:,:,100:min(190, numTotalFrames)), 3);
        figure('Name', 'Motion_Correction_ROI_Selection');
        imagesc(avgFrame); axis image off; colormap(gray);
        title('Draw motion correction ROI');
        ROI_rect = imrect;
        pos = round(getPosition(ROI_rect));
        close;
        roiWindow = stack(pos(2):(pos(2)+pos(4)), pos(1):(pos(1)+pos(3)), :);
    else
        roiWindow = stack;
    end

    % Apply motion correction
    if MOTION_CORRECTION
        fprintf('  Applying motion correction...\n');
        Y1 = single(roiWindow);
        if strcmp(CORRECTION_TYPE, 'rigid')
            options_rigid = NoRMCorreSetParms('d1', size(Y1,1), 'd2', size(Y1,2), ...
                'bin_width', 10, 'max_shift', 50, 'us_fac', 1, 'init_batch', 100);
            [Y_corrected, ~, ~] = normcorre(Y1, options_rigid);
        else
            error('Only rigid correction currently supported.');
        end
        vall = uint16(Y_corrected);
    else
        vall = roiWindow;
    end

    vall = permute(vall, [2 1 3]);

    % Store in cell array for consistency with multi-FOV processing
    all_image_stacks = {vall};
    num_fovs = 1;

else
    %% Multi-FOV processing
    cd(base_folder)
    fov_folders = dir('fov*');
    num_fovs = length(fov_folders);

    fprintf('Found %d FOV folders to process\n', num_fovs);

    all_image_stacks = cell(num_fovs, 1);

    for fov_idx = 1:num_fovs
        fprintf('\n--- Processing FOV %d of %d ---\n', fov_idx, num_fovs);

        current_fov_path = fullfile(fov_folders(fov_idx).folder, fov_folders(fov_idx).name);
        cd(current_fov_path)

        % Find TIFF file
        tif_file = dir('*.tif');
        if isempty(tif_file)
            warning('No TIFF file found in %s, skipping...', current_fov_path);
            continue;
        end

        tif_filename = fullfile(current_fov_path, tif_file.name);
        fprintf('Loading: %s\n', tif_file.name);

        % Load TIFF stack
        info = imfinfo(tif_filename);
        numFrames = numel(info);

        imageStack = zeros(info(1).Height, info(1).Width, numFrames, 'uint16');
        t = Tiff(tif_filename, 'r');
        warning off
        for frame_idx = 1:numFrames
            if mod(frame_idx, 500) == 0
                fprintf('  Loading frame %d/%d\n', frame_idx, numFrames);
            end
            t.setDirectory(frame_idx);
            imageStack(:,:,frame_idx) = t.read();
        end
        t.close();
        warning on

        % Motion correction
        if MOTION_CORRECTION
            fprintf('  Applying motion correction...\n');
            imageStack_single = single(imageStack);
            smoothed_images = imboxfilt3(imageStack_single, [1 1 11]);
            h_highpass = fspecial('gaussian', 50, 1) - fspecial('gaussian', 50, 25);
            filtered_images = imfilter(smoothed_images, h_highpass, 'replicate', 'same');
            clear smoothed_images

            options_rigid = NoRMCorreSetParms('d1', size(imageStack_single,1), ...
                'd2', size(imageStack_single,2), 'bin_width', 10, ...
                'max_shift', 50, 'us_fac', 1);
            [~, shifts, ~] = normcorre(filtered_images, options_rigid);
            clear filtered_images

            corrected_stack = zeros(size(imageStack_single), 'uint16');
            for frame_idx = 1:size(imageStack_single, 3)
                corrected_stack(:,:,frame_idx) = circshift(imageStack_single(:,:,frame_idx), ...
                    shifts(frame_idx).shifts);
            end
            imageStack = corrected_stack;
        end

        % Reorient data
        imageStack = permute(imageStack, [2 1 3]);
        all_image_stacks{fov_idx} = imageStack;
    end
end

%% ============================================================================
%  SECTION 3: EXTRACT FIBER TRACES FROM IMAGING DATA
%  ============================================================================

fprintf('\n=== EXTRACTING FIBER TRACES ===\n');

all_traces = [];
all_ROIs = {};

for fov_idx = 1:num_fovs
    fprintf('Processing FOV %d...\n', fov_idx);

    imageStack = all_image_stacks{fov_idx};

    if PROCESS_FULL_FIELD
        % Full field processing
        N = 1;
        ROIs = {ones(size(imageStack,1), size(imageStack,2))};
        traces = squeeze(nanmean(nanmean(imageStack, 1), 2));
        fprintf('  Full field average extracted\n');
    else
        % ROI-based processing
        averageFrame = clean_display_frame(imageStack, 15);  % ROI-selection background (see core/utils)

        mask = zeros(size(averageFrame));
        ROIs = {};

        cmin = min(averageFrame(:));
        cmax = prctile(averageFrame(:), 99.9);

        figure('Name', sprintf('ROI Selection - FOV %d', fov_idx));
        imagesc(averageFrame, [cmin, cmax]);
        axis image off; colormap(gray); colorbar;
        title(sprintf('Draw ROIs for FOV %d - Close polygon when done', fov_idx));

        ROI = drawpolygon;
        while isvalid(ROI)
            bw = createMask(ROI);
            ROIs{end+1} = bw;
            mask = mask | bw;
            imagesc(averageFrame .* (1 - mask * 0.1), [cmin, cmax]);
            axis image off; colormap(gray); colorbar;
            ROI = drawpolygon;
        end
        close;

        N = length(ROIs);
        traces = zeros(size(imageStack, 3), N);

        fprintf('  Extracting traces from %d ROIs\n', N);
        for roi_idx = 1:N
            neuron_signal = imageStack .* uint16(ROIs{roi_idx});
            traces(:, roi_idx) = squeeze(sum(neuron_signal, [1,2])) / sum(ROIs{roi_idx}(:));
        end
    end

    all_traces = [all_traces, traces];
    all_ROIs{fov_idx} = ROIs;
end

num_fibers = size(all_traces, 2);
fprintf('Total traces extracted: %d fibers\n', num_fibers);

%% ============================================================================
%  SECTION 4: ARTIFACT REMOVAL AND FILTERING
%  ============================================================================

fprintf('\n=== APPLYING ARTIFACT REMOVAL ===\n');

% Placeholder sampling rate (will be updated after ephys alignment)
if isempty(IMAGING_FS)
    IMAGING_FS = 500;  % Temporary default
end

% Remove end artifact + DBS harmonics (120-124, 130-132 Hz). See core/remove_stim_artifacts.m.
filtered_traces = remove_stim_artifacts(all_traces, IMAGING_FS);

fprintf('Artifact removal complete\n');

processed_traces = filtered_traces;  % Temporary - will be replaced after photobleaching correction
traces_detrended = filtered_traces;  % Temporary placeholder
traces_exp_corrected = filtered_traces;  % Temporary placeholder


%% ============================================================================
%  SECTION 5: LOAD OPEN EPHYS DATA AND ALIGN WITH IMAGING
%  ============================================================================

if LOAD_EPHYS_DATA
    fprintf('\n=== LOADING OPEN EPHYS DATA ===\n');

    %% Automatically construct Open Ephys path
    % Expected structure: <BASE_PATH_ROOT>\<MOUSE_NAME>\Open_Ephys\<DD-MM-YY>\[folder_with_R##]\Record Node 10X\

    % Build base Open Ephys path
    open_ephys_base = fullfile(BASE_PATH_ROOT, MOUSE_NAME, DATA_TYPE_OPEN_EPHYS, EXPERIMENT_TYPE, RECORDING_DATE);

    if ~exist(open_ephys_base, 'dir')
        warning('Open Ephys base folder does not exist: %s', open_ephys_base);
        [ephys_path, LOAD_EPHYS_DATA] = prompt_manual_ephys_file_selection();
        if ~LOAD_EPHYS_DATA
            EPHYS_LOADED = false;
        end
    else
        % Search for folder containing the Recording ID
        fprintf('Searching for recording folder containing "%s" in: %s\n', RECORDING_ID, open_ephys_base);

        all_folders = dir(open_ephys_base);
        all_folders = all_folders([all_folders.isdir]);
        all_folders = all_folders(~ismember({all_folders.name}, {'.', '..'}));

        % Find folder with Recording ID in name (using pattern to avoid substring matches)
        % e.g., "R1" should match "R1" but not "R10", "R11", etc.
        % Pattern ensures RECORDING_ID is not preceded or followed by a digit
        recording_folder = [];
        escaped_id = regexptranslate('escape', RECORDING_ID);
        % Pattern: RECORDING_ID must be at start/end or preceded/followed by non-digit
        % This matches "R1" in "_R1", "R1_", "R1-", etc., but not "R10", "R11"
        pattern = ['(^|[^0-9])' escaped_id '([^0-9]|$)'];
        for i = 1:length(all_folders)
            if ~isempty(regexp(all_folders(i).name, pattern, 'once'))
                recording_folder = all_folders(i).name;
                fprintf('Found recording folder: %s\n', recording_folder);
                break;
            end
        end

        if isempty(recording_folder)
            warning('Could not find folder containing "%s" in %s', RECORDING_ID, open_ephys_base);
            [ephys_path, LOAD_EPHYS_DATA] = prompt_manual_ephys_file_selection();
            if ~LOAD_EPHYS_DATA
                EPHYS_LOADED = false;
            end
        else
            % Now look for Record Node folder (try 103 first, then 104)
            recording_path = fullfile(open_ephys_base, recording_folder);

            record_node_path = [];
            if exist(fullfile(recording_path, 'Record Node 103'), 'dir')
                record_node_path = fullfile(recording_path, 'Record Node 103');
                fprintf('Found: Record Node 103\n');
            elseif exist(fullfile(recording_path, 'Record Node 104'), 'dir')
                record_node_path = fullfile(recording_path, 'Record Node 104');
                fprintf('Found: Record Node 104\n');
            else
                warning('Could not find "Record Node 103" or "Record Node 104" in %s', recording_path);
                [ephys_path, LOAD_EPHYS_DATA] = prompt_manual_ephys_file_selection();
                if ~LOAD_EPHYS_DATA
                    EPHYS_LOADED = false;
                end
            end

            if ~isempty(record_node_path)
                % Verify that .continuous files exist
                test_file = fullfile(record_node_path, sprintf('%d_RhythmData_Ch11.continuous', EPHYS_FILE_PREFIX));
                if exist(test_file, 'file')
                    ephys_path = record_node_path;
                    fprintf('Open Ephys path confirmed: %s\n', ephys_path);
                else
                    warning('Could not find expected .continuous files in %s', record_node_path);
                    [ephys_path, LOAD_EPHYS_DATA] = prompt_manual_ephys_file_selection();
                    if ~LOAD_EPHYS_DATA
                        EPHYS_LOADED = false;
                    end
                end
            end
        end
    end

    %% Load Open Ephys data if path was found
    if LOAD_EPHYS_DATA && exist('ephys_path', 'var')
        cd(ephys_path)

        % Load all channels
        fprintf('Loading Open Ephys channels...\n');

        try
            [lfp_data, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_Ch11.continuous', EPHYS_FILE_PREFIX));
            fprintf('  Ch11 (LFP): %d samples\n', length(lfp_data));
            if LOAD_mPFC_LFP

                fprintf('Loading mPFC differential LFP channels---\n');
                [lfp_data2, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_Ch1.continuous', EPHYS_FILE_PREFIX));
                fprintf('  Ch1 (mPFC+): %d samples\n', length(lfp_data2));

                [lfp_data3, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_Ch3.continuous', EPHYS_FILE_PREFIX));
                fprintf('  Ch3 (mPFC-): %d samples\n', length(lfp_data3));

                % Compute differential mPFC LFP signal
                lfp_data_mPFC = lfp_data2 - lfp_data3;
                fprintf('mPFC LFP computed (Ch1 - Ch3): %d samples\n', length(lfp_data_mPFC));
                % Verify lengths match
                if length(lfp_data2) ~= length(lfp_data3)
                    warning('Channel 1 and Channel 3 have different lengths!');
                    fprintf('  Ch1 length: %d\n', length(lfp_data2));
                    fprintf('  Ch3 length: %d\n', length(lfp_data3));
                    % Trim to shorter length
                    min_len = min(length(lfp_data2), length(lfp_data3));
                    lfp_data_mPFC = lfp_data2(1:min_len) - lfp_data3(1:min_len);
                    fprintf('  Trimmed to: %d samples\n', min_len);
                end

                mPFC_LFP_LOADED = true;
            else
                fprintf('mPFC LFP loading disabled (LOAD_mPFC_LFP = false)\n');
                mPFC_LFP_LOADED = false;
            end

            [camera_triggers, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_ADC1.continuous', EPHYS_FILE_PREFIX));
            fprintf('  ADC1 (Camera triggers): %d samples\n', length(camera_triggers));

            [trial_markers, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_ADC7.continuous', EPHYS_FILE_PREFIX));
            fprintf('  ADC7 (Trial markers): %d samples\n', length(trial_markers));

            [stim_channel_5, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_ADC5.continuous', EPHYS_FILE_PREFIX));
            [stim_channel_6, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_ADC6.continuous', EPHYS_FILE_PREFIX));
            fprintf('  ADC5/6 (Stim channels): loaded\n');

            [running_wheel_raw, ~, ~] = load_open_ephys_data(sprintf('%d_RhythmData_ADC4.continuous', EPHYS_FILE_PREFIX));
            fprintf('  ADC4 (Running wheel): %d samples\n', length(running_wheel_raw));

            %% Determine which channel contains stimulation pulses
            fprintf('\nDetermining stimulation pulse channel...\n');
            stim5_pulses = find(diff(abs(stim_channel_5)) > 0.02);
            stim6_pulses = find(diff(abs(stim_channel_6)) > 0.02);

            num_pulses_ch5 = length(stim5_pulses(diff(stim5_pulses) > 1));
            num_pulses_ch6 = length(stim6_pulses(diff(stim6_pulses) > 1));

            fprintf('  ADC5 detected pulses: %d\n', num_pulses_ch5);
            fprintf('  ADC6 detected pulses: %d\n', num_pulses_ch6);

            %% BASELINE TRIAL DETECTION
            % Check if user has manually set IS_BASELINE_TRIAL
            if isempty(IS_BASELINE_TRIAL)
                % Auto-detect: baseline trial if very few pulses detected in both channels
                % Threshold: < 5 pulses total suggests baseline trial
                BASELINE_PULSE_THRESHOLD = 5;
                total_pulses = num_pulses_ch5 + num_pulses_ch6;
                
                if total_pulses < BASELINE_PULSE_THRESHOLD
                    IS_BASELINE_TRIAL = true;
                    fprintf('\n=== BASELINE TRIAL DETECTED ===\n');
                    fprintf('  Total pulses detected: %d (threshold: %d)\n', total_pulses, BASELINE_PULSE_THRESHOLD);
                    fprintf('  No stimulation period will be defined\n');
                else
                    IS_BASELINE_TRIAL = false;
                    fprintf('\n=== STIMULATION TRIAL DETECTED ===\n');
                    fprintf('  Total pulses detected: %d\n', total_pulses);
                end
            else
                % User override
                if IS_BASELINE_TRIAL
                    fprintf('\n=== BASELINE TRIAL (USER OVERRIDE) ===\n');
                else
                    fprintf('\n=== STIMULATION TRIAL (USER OVERRIDE) ===\n');
                end
            end

            if ~IS_BASELINE_TRIAL
                % Only determine pulse channels if this is a stimulation trial
                if num_pulses_ch5 > num_pulses_ch6
                    stim_pulses = stim_channel_5;
                    stim_onset_trigger = stim_channel_6;
                    fprintf('  -> ADC5 contains stimulation pulses\n');
                else
                    stim_pulses = stim_channel_6;
                    stim_onset_trigger = stim_channel_5;
                    fprintf('  -> ADC6 contains stimulation pulses\n');
                end
            else
                % For baseline trials, still assign channels for compatibility
                % but they won't be used for detection
                if num_pulses_ch5 > num_pulses_ch6
                    stim_pulses = stim_channel_5;
                    stim_onset_trigger = stim_channel_6;
                else
                    stim_pulses = stim_channel_6;
                    stim_onset_trigger = stim_channel_5;
                end
            end


            %% Automatic stimulation period detection
            if LOAD_EPHYS_DATA && AUTO_DETECT_STIMULATION && ~IS_BASELINE_TRIAL

                fprintf('\n=== DETECTING STIMULATION PERIOD AUTOMATICALLY ===\n');

                % Define thresholds
                STIM_ONSET_THRESHOLD = 0.1;  % Adjust based on your trigger voltage
                STIM_PULSE_THRESHOLD = 0.1; % For detecting individual pulses

                % First, align stimulation channels to imaging timebase
                stim_onset_ephys_full = stim_onset_trigger;  % The channel with onset marker
                stim_pulses_ephys_full = stim_pulses;        % The channel with pulse train

                % Detect stimulus ONSET from the pulse train channel (not trigger channel).
                % Detection kernel extracted to core/detect_stim_onset.m (unit-tested);
                % prefer the first onset within the trial period if markers exist.
                fprintf('Detecting stimulation onset from pulse train...\n');
                if exist('trial_starts', 'var') && ~isempty(trial_starts)
                    stim_onset_sample = detect_stim_onset(stim_pulses_ephys_full, STIM_ONSET_THRESHOLD, trial_starts(1));
                else
                    stim_onset_sample = detect_stim_onset(stim_pulses_ephys_full, STIM_ONSET_THRESHOLD);
                end

                if ~isempty(stim_onset_sample)
                    fprintf('  Stimulus onset detected at Open Ephys sample: %d\n', stim_onset_sample);

                    % Convert to time
                    stim_onset_time = stim_onset_sample / EPHYS_FS;
                    fprintf('  Stimulus onset time: %.3f seconds\n', stim_onset_time);

                else
                    warning('Could not detect stimulus onset from trigger channel');
                    stim_onset_sample = [];
                    stim_onset_time = [];
                end

                %% Set stimulus OFFSET using fixed duration (like voltage imaging code)
                fprintf('Setting stimulation offset using fixed duration...\n');

                if ~isempty(stim_onset_sample)
                    % Use fixed 10-second duration (adjust as needed for your experiments)
                    STIMULATION_DURATION = 10.0;  % seconds
                    stim_offset_sample = stim_onset_sample + round(STIMULATION_DURATION * EPHYS_FS);
                    stim_offset_time = stim_onset_time + STIMULATION_DURATION;

                    fprintf('  Stimulus offset set to Open Ephys sample: %d\n', stim_offset_sample);
                    fprintf('  Stimulus offset time: %.3f seconds\n', stim_offset_time);
                else
                    warning('No stimulus onset detected. Cannot set offset.');
                    stim_offset_sample = [];
                    stim_offset_time = [];
                end

                % Update STIM_PERIOD if detection was successful
                if ~isempty(stim_onset_time) && ~isempty(stim_offset_time)
                    % Store original hardcoded values for reference
                    STIM_PERIOD_HARDCODED = STIM_PERIOD;

                    % Convert ephys time to imaging time for visualization
                    % We need to wait until camera triggers are processed to do this conversion
                    % For now, store the ephys sample indices for later conversion
                    STIM_ONSET_EPHYS_SAMPLE = stim_onset_sample;
                    STIM_OFFSET_EPHYS_SAMPLE = stim_offset_sample;

                    % Don't update STIM_PERIOD here - it will be updated after camera alignment

                    fprintf('\n=== STIMULATION PERIOD UPDATED ===\n');
                    fprintf('  Original (hardcoded): [%.2f, %.2f] seconds\n', ...
                        STIM_PERIOD_HARDCODED(1), STIM_PERIOD_HARDCODED(2));
                    fprintf('  Detected (automatic): [%.2f, %.2f] seconds\n', ...
                        STIM_PERIOD(1), STIM_PERIOD(2));
                    fprintf('  Duration: %.2f seconds\n', STIM_PERIOD(2) - STIM_PERIOD(1));

                    % Calculate offset from hardcoded values
                    onset_diff = STIM_PERIOD(1) - STIM_PERIOD_HARDCODED(1);
                    offset_diff = STIM_PERIOD(2) - STIM_PERIOD_HARDCODED(2);

                    if abs(onset_diff) > 1.0 || abs(offset_diff) > 1.0
                        warning('Large difference from hardcoded values detected!');
                        fprintf('  Onset difference: %.2f seconds\n', onset_diff);
                        fprintf('  Offset difference: %.2f seconds\n', offset_diff);
                    end

                    % Store detection metadata
                    STIM_DETECTION = struct();
                    STIM_DETECTION.method = 'automatic';
                    STIM_DETECTION.onset_sample = stim_onset_sample;
                    STIM_DETECTION.offset_sample = stim_offset_sample;
                    STIM_DETECTION.onset_threshold = STIM_ONSET_THRESHOLD;
                    STIM_DETECTION.pulse_threshold = STIM_PULSE_THRESHOLD;
                    STIM_DETECTION.onset_diff_from_hardcoded = onset_diff;
                    STIM_DETECTION.offset_diff_from_hardcoded = offset_diff;

                else
                    warning('Automatic detection failed. Using hardcoded STIM_PERIOD values.');
                    fprintf('  Using: [%.2f, %.2f] seconds\n', STIM_PERIOD(1), STIM_PERIOD(2));

                    STIM_DETECTION = struct();
                    STIM_DETECTION.method = 'hardcoded';
                    STIM_DETECTION.reason = 'automatic_detection_failed';
                end

                fprintf('\n');
            end

            %% Ensure STIM_DETECTION exists even if detection was skipped
            if ~exist('STIM_DETECTION', 'var')
                if IS_BASELINE_TRIAL
                    STIM_DETECTION = struct();
                    STIM_DETECTION.method = 'baseline_trial';
                    STIM_DETECTION.reason = 'no_stimulation_pulses_detected';
                    fprintf('\n=== BASELINE TRIAL: No stimulation detection performed ===\n');
                else
                    STIM_DETECTION = struct();
                    STIM_DETECTION.method = 'skipped';
                    STIM_DETECTION.reason = 'ephys_not_loaded_or_detection_disabled';
                end
            end


            %% Convert stimulation period from ephys time to imaging frame indices
            % This will be used later after camera triggers are aligned
            % Store the ephys sample indices for later conversion
            if exist('stim_onset_sample', 'var') && ~isempty(stim_onset_sample)
                STIM_ONSET_EPHYS_SAMPLE = stim_onset_sample;
                STIM_OFFSET_EPHYS_SAMPLE = stim_offset_sample;
            else
                STIM_ONSET_EPHYS_SAMPLE = [];
                STIM_OFFSET_EPHYS_SAMPLE = [];
            end
            %% Identify valid trial periods
            fprintf('\nIdentifying valid trial periods...\n');
            TRIAL_THRESHOLD = 0.5;
            trial_active = trial_markers > TRIAL_THRESHOLD;
            trial_starts = find(diff([0; trial_active]) == 1);
            trial_stops = find(diff([trial_active; 0]) == -1);
            fprintf('  Found %d trial periods\n', length(trial_starts));

            %% Extract and filter camera triggers
            fprintf('\nExtracting camera trigger timestamps...\n');
            all_camera_triggers = find(diff(camera_triggers) > 0.5);
            trigger_intervals = diff(all_camera_triggers);
            valid_trigger_spacing = [true; trigger_intervals > 8];
            camera_trigger_indices = all_camera_triggers(valid_trigger_spacing);

            fprintf('  Total camera triggers detected: %d\n', length(camera_trigger_indices));

            % Filter to valid trials
            if ~isempty(trial_starts)
                fprintf('  Filtering triggers to valid trial periods...\n');
                triggers_in_trials = false(size(camera_trigger_indices));
                for trial_idx = 1:length(trial_starts)
                    in_this_trial = (camera_trigger_indices >= trial_starts(trial_idx)) & ...
                        (camera_trigger_indices <= trial_stops(trial_idx));
                    triggers_in_trials = triggers_in_trials | in_this_trial;
                end
                camera_trigger_indices = camera_trigger_indices(triggers_in_trials);
            end

            fprintf('  Valid camera triggers: %d\n', length(camera_trigger_indices));

            %% Calculate actual imaging sampling rate
            fprintf('\nCalculating actual imaging sampling rate...\n');
            trigger_intervals = diff(camera_trigger_indices);
            trigger_intervals_sec = trigger_intervals / EPHYS_FS;
            instantaneous_frame_rates = 1 ./ trigger_intervals_sec;

            median_frame_rate = median(instantaneous_frame_rates);
            mean_frame_rate = mean(instantaneous_frame_rates);
            std_frame_rate = std(instantaneous_frame_rates);

            fprintf('  Median frame rate: %.2f Hz\n', median_frame_rate);
            fprintf('  Mean frame rate: %.2f Hz\n', mean_frame_rate);
            fprintf('  Std deviation: %.2f Hz\n', std_frame_rate);

            IMAGING_FS = median_frame_rate;

            if std_frame_rate > 1.0
                warning('Frame rate variability detected (std = %.2f Hz)', std_frame_rate);
            end

            %% Verify alignment
            num_fiber_frames = size(processed_traces, 1);
            num_trigger_frames = length(camera_trigger_indices);

            fprintf('  Imaging frames: %d\n', num_fiber_frames);
            fprintf('  Camera triggers: %d\n', num_trigger_frames);

            if abs(num_trigger_frames - num_fiber_frames) > 10
                warning('Significant mismatch between triggers and imaging frames!');
                fprintf('  Difference: %d frames\n', abs(num_trigger_frames - num_fiber_frames));
            end

            % Align to shorter length
            min_frames = min(num_fiber_frames, num_trigger_frames);
            camera_trigger_indices = camera_trigger_indices(1:min_frames);
            processed_traces = processed_traces(1:min_frames, :);
            filtered_traces = filtered_traces(1:min_frames, :);
            traces_detrended = traces_detrended(1:min_frames, :);
            traces_exp_corrected = traces_exp_corrected(1:min_frames, :);

            %% Align ephys signals to imaging timebase
            fprintf('\nAligning ephys data to imaging frames...\n');
            lfp_aligned = lfp_data(camera_trigger_indices);

            % Conditionally align mPFC LFP
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                lfp_aligned_mPFC = lfp_data_mPFC(camera_trigger_indices);
                fprintf('  mPFC LFP aligned: %d samples\n', length(lfp_aligned_mPFC));
            end

            % Process running signal
            running_velocity = diff(running_wheel_raw) > 1;
            running_velocity_smooth = fastsmooth(running_velocity, 30, 1, 1);
            running_velocity_aligned = running_velocity_smooth(camera_trigger_indices);

            % Align stimulation signals
            stim_pulses_aligned = stim_pulses(camera_trigger_indices);
            stim_onset_aligned = stim_onset_trigger(camera_trigger_indices);

            fprintf('  LFP downsampled from %d Hz to %.2f Hz\n', EPHYS_FS, IMAGING_FS);
            fprintf('  Alignment complete\n');

            % Remove outliers from LFP
            lfp_aligned = replace_outliers_with_median(lfp_aligned, 10);

            % Remove outliers from mPFC LFP
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                lfp_aligned_mPFC = replace_outliers_with_median(lfp_aligned_mPFC, 10);
                fprintf('  mPFC LFP outliers removed\n');
            end

            % Smooth running velocity
            running_velocity_smooth = fastsmooth(running_velocity_aligned * 1000, 40, 1, 1);

            EPHYS_LOADED = true;

            % Update STIM_PERIOD with properly aligned timing (only for stimulation trials)
            if IS_BASELINE_TRIAL
                % Baseline trial: set STIM_PERIOD to empty to indicate no stimulation
                STIM_PERIOD = [];
                fprintf('Baseline trial: STIM_PERIOD set to empty (no stimulation period)\n');
            elseif ~isempty(STIM_ONSET_EPHYS_SAMPLE) && exist('camera_trigger_indices', 'var')
                % Find closest camera trigger to detected stim onset/offset
                [~, stim_onset_frame] = min(abs(camera_trigger_indices - STIM_ONSET_EPHYS_SAMPLE));
                [~, stim_offset_frame] = min(abs(camera_trigger_indices - STIM_OFFSET_EPHYS_SAMPLE));

                % Convert to imaging time
                stim_onset_time_aligned = stim_onset_frame / IMAGING_FS;
                stim_offset_time_aligned = stim_offset_frame / IMAGING_FS;

                % NOW update STIM_PERIOD for visualization
                STIM_PERIOD = [stim_onset_time_aligned, stim_offset_time_aligned];

                fprintf('STIM_PERIOD updated for visualization: [%.2f, %.2f] seconds\n', ...
                    STIM_PERIOD(1), STIM_PERIOD(2));
            else
                fprintf('Using original STIM_PERIOD values for visualization\n');
            end

        catch ME
            warning('Error loading Open Ephys data: %s', ME.message);
            fprintf('Continuing without ephys analysis.\n');
            LOAD_EPHYS_DATA = false;
            EPHYS_LOADED = false;
        end
    end
else
    EPHYS_LOADED = false;
    fprintf('Skipping Open Ephys data loading (disabled in configuration)\n');
end

%% Handle baseline trial detection when ephys is not loaded
if ~EPHYS_LOADED && isempty(IS_BASELINE_TRIAL)
    % If ephys not loaded and user hasn't specified, default to stimulation trial
    % (maintains backward compatibility)
    IS_BASELINE_TRIAL = false;
    fprintf('Ephys not loaded: Assuming stimulation trial (use IS_BASELINE_TRIAL = true to override)\n');
elseif ~EPHYS_LOADED && IS_BASELINE_TRIAL
    fprintf('Baseline trial mode: User override (ephys not loaded)\n');
    % Initialize STIM_DETECTION for baseline trial
    if ~exist('STIM_DETECTION', 'var')
        STIM_DETECTION = struct();
        STIM_DETECTION.method = 'baseline_trial';
        STIM_DETECTION.reason = 'user_override_no_ephys';
    end
end

%% ============================================================================
%  SECTION 6: PHOTOBLEACHING CORRECTION
%  ============================================================================

if APPLY_PHOTOBLEACHING_CORRECTION
    fprintf('\n=== APPLYING PHOTOBLEACHING CORRECTION ===\n');
    fprintf('Note: Working with ephys-aligned data (length: %d frames)\n', size(filtered_traces, 1));

    % Determine stimulus period for correction using detected values
    if PROCESS_SINGLE_TRIAL
        if IS_BASELINE_TRIAL
            % Baseline trial: use entire 60s trace for photobleaching correction
            stim_onset_frame = size(filtered_traces, 1);  % Use entire trace
            stim_offset_frame = size(filtered_traces, 1);
            fprintf('Baseline trial: Using entire 60s trace for photobleaching correction\n');
            fprintf('  Total frames: %d (%.2f seconds)\n', stim_onset_frame, stim_onset_frame/IMAGING_FS);
        elseif ~isempty(STIM_ONSET_EPHYS_SAMPLE) && exist('camera_trigger_indices', 'var')
            % Use automatically detected stimulation period aligned to imaging frames
            % Find closest camera trigger to detected stim onset
            [~, stim_onset_frame] = min(abs(camera_trigger_indices - STIM_ONSET_EPHYS_SAMPLE));
            [~, stim_offset_frame] = min(abs(camera_trigger_indices - STIM_OFFSET_EPHYS_SAMPLE));

            fprintf('Using automatically detected stimulation period:\n');
            fprintf('  Onset: frame %d (%.2f seconds)\n', stim_onset_frame, stim_onset_frame/IMAGING_FS);
            fprintf('  Offset: frame %d (%.2f seconds)\n', stim_offset_frame, stim_offset_frame/IMAGING_FS);
        else
            % Fallback to time-based calculation using updated STIM_PERIOD
            if ~isempty(STIM_PERIOD)
                stim_onset_frame = round(STIM_PERIOD(1) * IMAGING_FS) + 1;
                stim_offset_frame = round(STIM_PERIOD(2) * IMAGING_FS) + 1;

                fprintf('Using time-based stimulation period (fallback):\n');
                fprintf('  Onset: frame %d (%.2f seconds)\n', stim_onset_frame, stim_onset_frame/IMAGING_FS);
                fprintf('  Offset: frame %d (%.2f seconds)\n', stim_offset_frame, stim_offset_frame/IMAGING_FS);
            else
                % If STIM_PERIOD is empty and no detection, use first 80% as baseline
                stim_onset_frame = round(size(filtered_traces, 1) * 0.8);
                stim_offset_frame = size(filtered_traces, 1);
                fprintf('No stimulation period available: Using first 80%% as baseline\n');
            end
        end
    else
        % For multi-FOV, use first 20% of recording as baseline
        stim_onset_frame = round(size(filtered_traces, 1) * 0.2);
        stim_offset_frame = size(filtered_traces, 1);
        fprintf('Multi-FOV mode: Using first 20%% as baseline\n');
    end

    % Ensure frames are within bounds
    stim_onset_frame = max(1, min(stim_onset_frame, size(filtered_traces, 1)));
    stim_offset_frame = max(stim_onset_frame, min(stim_offset_frame, size(filtered_traces, 1)));

    fprintf('Using frames 1-%d as baseline for photobleaching correction\n', stim_onset_frame);

    % Per-trace photobleaching correction (linear detrend + double-exponential).
    % Extracted to core/correct_photobleaching.m (unit-tested in core/tests/).
    % Single-trial keeps its historical behaviour: NO baseline clamp and NO
    % short-baseline guard (the multi-trial script uses the guarded defaults).
    pb_opts = struct('BaselineClampMin', false, 'ShortBaselineGuard', false);
    [traces_detrended, traces_exp_corrected, fit_diag] = correct_photobleaching( ...
        filtered_traces, stim_onset_frame, IS_BASELINE_TRIAL, IMAGING_FS, pb_opts);

    % Restore the fit "leftover" variables that later diagnostic figures and the
    % saved data struct read (isfield reproduces the original exist('var') guards).
    double_exp_function = fit_diag.double_exp_function;
    pre_stim_trace = fit_diag.pre_stim_trace;
    if isfield(fit_diag, 'time_pre_stim'), time_pre_stim = fit_diag.time_pre_stim; end
    if isfield(fit_diag, 'time_full'),     time_full     = fit_diag.time_full;     end
    if isfield(fit_diag, 'fitted_params'), fitted_params = fit_diag.fitted_params; end

    % Photobleaching correction outputs F_corr(t) (corrected fluorescence, not ΔF/F)
    % ΔF/F will be computed in the next section
    % Do NOT set processed_traces here - it will be set after ΔF/F computation

    fprintf('Photobleaching correction completed\n');
else
    % If no photobleaching correction, use filtered traces as F_corr(t)
    % (they will still go through ΔF/F computation)
    traces_detrended = filtered_traces;
    traces_exp_corrected = filtered_traces;
    fprintf('Photobleaching correction disabled: Using filtered traces as F_corr(t)\n');
end

%% ============================================================================
%  SECTION 6B: ΔF/F COMPUTATION (STATIC BASELINE METHOD)
%  ============================================================================

fprintf('\n=== COMPUTING ΔF/F (STATIC BASELINE METHOD) ===\n');

% Determine baseline window for F0 calculation
if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
    % For stimulation trials: use [-60, 0] seconds before stim onset
    baseline_start_time = STIM_PERIOD(1) - 60;  % 60 seconds before stim
    baseline_end_time = STIM_PERIOD(1);          % Up to stim onset
    
    % Convert to frame indices (ensure within bounds)
    baseline_start_frame = max(1, round(baseline_start_time * IMAGING_FS) + 1);
    baseline_end_frame = min(round(baseline_end_time * IMAGING_FS) + 1, size(traces_exp_corrected, 1));
    
    % Ensure valid window
    if baseline_start_frame >= baseline_end_frame
        warning('Baseline window invalid, using first 60 seconds instead');
        baseline_start_frame = 1;
        baseline_end_frame = min(round(60 * IMAGING_FS) + 1, size(traces_exp_corrected, 1));
    end
    
    fprintf('Stimulation trial: Using baseline window [%.1f, %.1f] seconds before stim\n', ...
        max(0, baseline_start_time), baseline_end_time);
    fprintf('  Baseline frames: %d to %d (actual time: [%.2f, %.2f] s)\n', ...
        baseline_start_frame, baseline_end_frame, ...
        (baseline_start_frame-1)/IMAGING_FS, (baseline_end_frame-1)/IMAGING_FS);
else
    % For baseline trials: use entire 60s trace for F0 calculation
    % (baseline trials are 60s total, so use full trace)
    baseline_start_frame = 1;
    baseline_end_frame = size(traces_exp_corrected, 1);  % Use entire trace
    
    fprintf('Baseline trial: Using entire 60s trace for F0 calculation\n');
    fprintf('  Baseline frames: %d to %d (%.2f seconds)\n', ...
        baseline_start_frame, baseline_end_frame, ...
        (baseline_end_frame-1)/IMAGING_FS);
end

% Compute ΔF/F (static baseline) plus z-scored ΔF/F for the chosen window.
% Window SELECTION is above; this just applies it. See core/compute_deltaF_F.m.
[traces_deltaF_F, F0_values, processed_traces_zscored] = ...
    compute_deltaF_F(traces_exp_corrected, baseline_start_frame, baseline_end_frame);

% For visualization, use ΔF/F instead of z-scored
processed_traces = traces_deltaF_F;

fprintf('Final processed traces prepared (length: %d frames)\n', size(processed_traces, 1));
fprintf('  Using ΔF/F for visualization (z-scored version also available)\n');

%% Verify all traces have same length after correction
if EPHYS_LOADED
    expected_length = length(camera_trigger_indices);
    actual_length = size(processed_traces, 1);

    if expected_length ~= actual_length
        warning('Length mismatch after photobleaching correction!');
        fprintf('  Expected (from camera triggers): %d\n', expected_length);
        fprintf('  Actual (processed traces): %d\n', actual_length);
        fprintf('  Trimming to match...\n');

        min_len = min(expected_length, actual_length);
        processed_traces = processed_traces(1:min_len, :);
        traces_detrended = traces_detrended(1:min_len, :);
        traces_exp_corrected = traces_exp_corrected(1:min_len, :);
        filtered_traces = filtered_traces(1:min_len, :);
    else
        fprintf('✓ All traces properly aligned (%d frames)\n', actual_length);
    end
end

%% ============================================================================
%  SECTION 7: TIME-FREQUENCY ANALYSIS
%  ============================================================================

fprintf('\n=== PERFORMING TIME-FREQUENCY ANALYSIS ===\n');

% Prepare voltage signal
voltage_signal = processed_traces(:);  % Use first trace or concatenate all
tvec = (1:length(voltage_signal)) / IMAGING_FS;

% if EPHYS_LOADED
%     %% Combined Fiber + LFP analysis using FieldTrip
%     fprintf('Creating FieldTrip data structure...\n');
%
%     % Fix: Use only the first fiber trace for FieldTrip analysis instead of concatenating all
%     voltage_signal_single = processed_traces(:, 1);  % Use first trace only
%     tvec = (1:length(voltage_signal_single)) / IMAGING_FS;
%
%     lfp_ft = [];
%     lfp_ft.trial{1}(1,:) = voltage_signal_single;
%     lfp_ft.trial{1}(2,:) = lfp_aligned;
%     lfp_ft.time{1} = tvec;
%     lfp_ft.label{1} = 'Fiber';
%     lfp_ft.label{2} = 'LFP';
%
%     % Configure time-frequency analysis
%     cfg = [];
%     cfg.method = 'mtmconvol';
%     cfg.output = 'fourier';
%     cfg.taper = 'hanning';
%     cfg.keeptapers = 'yes';
%     cfg.keeptrials = 'yes';
%     cfg.tapsmofrq = 5;
%     cfg.channel = 'all';
%     cfg.foi = 2:1:70;
%     cfg.toi = lfp_ft.time{1}(1:1:end);
%     cfg.width = 8;
%     cfg.t_ftimwin = ones(1, length(cfg.foi)) * 0.5;
%
%     fprintf('Running FieldTrip time-frequency analysis...\n');
%     freq_result = ft_freqanalysis(cfg, lfp_ft);
%     fprintf('Analysis complete\n');
%
%     % Extract spectrograms
%     spectrogram_fiber = squeeze(abs(freq_result.fourierspctrm(1, 1, :, :)));
%     spectrogram_lfp = squeeze(abs(freq_result.fourierspctrm(1, 2, :, :)));
%     time_vector = freq_result.time;
%     freq_vector = freq_result.freq';
%
%     % Extract phase information
%     phase_fiber = squeeze(angle(freq_result.fourierspctrm(1, 1, :, :)));
%     phase_lfp = squeeze(angle(freq_result.fourierspctrm(1, 2, :, :)));
%
% else
%     %% Fiber-only analysis
%     fprintf('Performing fiber-only time-frequency analysis...\n');
%
%     % For each trace, compute spectrogram
%     num_traces = size(processed_traces, 2);
%     spectrograms_all = cell(num_traces, 1);
%
%     for trace_idx = 1:num_traces
%         Vx = processed_traces(:, trace_idx);
%         [s, w, t] = spectrogram(Vx, 480, 449, 1:1:150, IMAGING_FS);
%         spectrograms_all{trace_idx} = struct('s', s, 'w', w, 't', t);
%     end
%
%     % Use first trace for main analysis
%     spectrogram_fiber = abs(spectrograms_all{1}.s);
%     freq_vector = spectrograms_all{1}.w;
%     time_vector = spectrograms_all{1}.t;
% end

%% ============================================================================
%  CREATE SAVE DIRECTORY
%  ============================================================================

fprintf('\n=== CREATING SAVE DIRECTORY ===\n');

% Mouse base path: <BASE_PATH_ROOT>\<MOUSE_NAME>\
mouse_base_path = fullfile(BASE_PATH_ROOT, MOUSE_NAME);

if ~exist(mouse_base_path, 'dir')
    error('Mouse base path does not exist: %s', mouse_base_path);
end

% Create Fiber_Voltage_Processed folder if it doesn't exist
processed_base = fullfile(mouse_base_path, 'Fiber_Voltage_Processed');
if ~exist(processed_base, 'dir')
    mkdir(processed_base);
    fprintf('Created directory: %s\n', processed_base);
end

% Create session-specific folder: Date-RecordingID (e.g., 23-09-25-R14)
session_folder_name = sprintf('%s-%s', METADATA.recording_date, METADATA.recording_id);
save_directory = fullfile(processed_base, session_folder_name);

if ~exist(save_directory, 'dir')
    mkdir(save_directory);
    fprintf('Created session directory: %s\n', save_directory);
end

fprintf('All results will be saved to: %s\n', save_directory);

%% Visualize stimulation detection (diagnostic plot)
if exist('STIM_DETECTION', 'var') && EPHYS_LOADED && ~IS_BASELINE_TRIAL
    if strcmp(STIM_DETECTION.method, 'automatic')
    figure('Name', 'Stimulation Detection Verification', 'Color', 'w', 'Position', [100, 100, 1400, 500]);

    % Get the pulse channel data (the one we actually used for detection)
    if num_pulses_ch5 > num_pulses_ch6
        pulse_channel_data = stim_channel_5;
        channel_name = 'ADC5';
    else
        pulse_channel_data = stim_channel_6;
        channel_name = 'ADC6';
    end

    time_ephys = (1:length(pulse_channel_data)) / EPHYS_FS;

    % Single clean plot showing only the pulse train with detection markers
    hold on;

    % Plot the pulse train
    plot(time_ephys, pulse_channel_data, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Stimulation Pulses');

    % Mark detected onset
    if ~isempty(STIM_ONSET_EPHYS_SAMPLE)
        xline(STIM_ONSET_EPHYS_SAMPLE/EPHYS_FS, 'g-', 'LineWidth', 3, 'DisplayName', 'Detected Onset');
    end

    % Mark calculated offset (onset + 10s)
    if ~isempty(STIM_OFFSET_EPHYS_SAMPLE)
        xline(STIM_OFFSET_EPHYS_SAMPLE/EPHYS_FS, 'r-', 'LineWidth', 3, 'DisplayName', 'Calculated Offset (+10s)');
    end

    % Add hardcoded references for comparison (if they differ significantly)
    if exist('STIM_PERIOD_HARDCODED', 'var')
        onset_diff = abs(STIM_ONSET_EPHYS_SAMPLE/EPHYS_FS - STIM_PERIOD_HARDCODED(1));
        if onset_diff > 1.0  % Only show if difference > 1 second
            xline(STIM_PERIOD_HARDCODED(1), 'g:', 'LineWidth', 2, 'DisplayName', 'Hardcoded Onset');
            xline(STIM_PERIOD_HARDCODED(2), 'r:', 'LineWidth', 2, 'DisplayName', 'Hardcoded Offset');
        end
    end

    % Add shaded stimulation period
    if ~isempty(STIM_ONSET_EPHYS_SAMPLE) && ~isempty(STIM_OFFSET_EPHYS_SAMPLE)
        y_limits = ylim;
        patch([STIM_ONSET_EPHYS_SAMPLE/EPHYS_FS, STIM_OFFSET_EPHYS_SAMPLE/EPHYS_FS, ...
            STIM_OFFSET_EPHYS_SAMPLE/EPHYS_FS, STIM_ONSET_EPHYS_SAMPLE/EPHYS_FS], ...
            [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [1, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'Stimulation Period');
    end

    % Formatting
    title(sprintf('Stimulation Detection from %s (Biphasic Pulse Train)', channel_name), ...
        'FontSize', 16, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Voltage (V)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 12);
    grid on;
    set(gca, 'FontSize', 12, 'LineWidth', 1.5);

    % Add detection summary text
    if ~isempty(STIM_ONSET_EPHYS_SAMPLE)
        detection_text = sprintf('Onset: %.2f s | Duration: %.1f s | Method: First pulse + fixed duration', ...
            STIM_ONSET_EPHYS_SAMPLE/EPHYS_FS, ...
            (STIM_OFFSET_EPHYS_SAMPLE - STIM_ONSET_EPHYS_SAMPLE)/EPHYS_FS);
        text(0.02, 0.95, detection_text, 'Units', 'normalized', 'FontSize', 11, ...
            'BackgroundColor', 'white', 'EdgeColor', 'black', 'Margin', 5);
    end

    % Save diagnostic figure
    diag_fig_filename = sprintf('StimDetection_Diagnostic_%s', session_folder_name);
    saveas(gcf, fullfile(save_directory, [diag_fig_filename '.fig']));
    saveas(gcf, fullfile(save_directory, [diag_fig_filename '.png']));
    fprintf('Stimulation detection diagnostic saved: %s\n', diag_fig_filename);
    end
elseif IS_BASELINE_TRIAL && EPHYS_LOADED
    % Create a simple diagnostic plot for baseline trials
    figure('Name', 'Baseline Trial Detection', 'Color', 'w', 'Position', [100, 100, 800, 400]);
    text(0.5, 0.5, sprintf('BASELINE TRIAL DETECTED\n\nNo stimulation pulses found.\nTotal pulses: %d + %d = %d', ...
        num_pulses_ch5, num_pulses_ch6, num_pulses_ch5 + num_pulses_ch6), ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
        'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.2 0.6 0.2]);
    axis off;
    diag_fig_filename = sprintf('BaselineTrial_Diagnostic_%s', session_folder_name);
    saveas(gcf, fullfile(save_directory, [diag_fig_filename '.fig']));
    saveas(gcf, fullfile(save_directory, [diag_fig_filename '.png']));
    fprintf('Baseline trial diagnostic saved: %s\n', diag_fig_filename);
end
%% ============================================================================
%  SECTION 8: ENHANCED VISUALIZATION WITH MULTI-FIBER SUPPORT
%  ============================================================================

fprintf('\n=== GENERATING VISUALIZATIONS ===\n');

% Pre-compute z-scored signals once to avoid repeated computation
if EPHYS_LOADED
    lfp_z_scored = zscore(lfp_aligned);
    motion_z_scored = zscore(running_velocity_smooth);
    if LOAD_mPFC_LFP && mPFC_LFP_LOADED
        mPFC_z_scored = zscore(lfp_aligned_mPFC);
    end
end

% Define color schemes
viridis_colors = [
    0.267004, 0.004874, 0.329415;  % Dark purple (Pre-stim)
    0.127568, 0.566949, 0.550556;  % Teal (Stim)
    0.993248, 0.906157, 0.143936   % Yellow (Post-stim)
    ];

% Define consistent colors for different signals
COLOR_LFP = [0.2, 0.2, 0.25];           % Dark grey-black for LFP
COLOR_FIBER = [0.127568, 0.566949, 0.550556];  % Teal for Fiber (from viridis)
COLOR_MOTION = [0.993248, 0.7, 0.4];    % Orange-pink for motion
COLOR_STIM_SHADE = [1, 0.9, 0.9];       % Light red for stim period shading

% Enhanced multi-fiber colors
if num_fibers > 1
    % Generate distinct colors for multiple fibers
    fiber_colors = lines(min(num_fibers, MAX_FIBERS_DISPLAY));
    % Ensure good contrast
    if num_fibers > 6
        fiber_colors = [fiber_colors; jet(num_fibers - 6)];
    end
else
    fiber_colors = COLOR_FIBER;
end

% Viridis-compatible colors for correction methods
COLOR_LINEAR = [0.27, 0.49, 0.77];    % soft steel blue
COLOR_EXPONENTIAL = [0.83, 0.45, 0.37];    % muted coral red
COLOR_mPFC = [0.8, 0.3, 0.3];  % Reddish for mPFC

% Colors for behavioral states
COLOR_RUNNING = [0.163625, 0.471133, 0.558148];  % Greenish-blue from viridis
COLOR_REST = [0.741388, 0.173449, 0.149561];     % Darkish red

polarity = -1 * INVERT_TRACE + 1 * ~INVERT_TRACE;

% Set default figure properties
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultAxesFontWeight', 'bold');
set(0, 'DefaultAxesLineWidth', 1.5);
set(0, 'DefaultLineLineWidth', 2);

%% ========================================================================
%  NEW: MULTI-FIBER COMBINED VISUALIZATION - SIMPLIFIED
%  ========================================================================

if GENERATE_MULTI_FIBER_PLOT && num_fibers > 1
    fprintf('\n=== GENERATING MULTI-FIBER COMBINED VISUALIZATION ===\n');

    % Prepare time vector for all fibers
    tvec_fiber = (1:size(processed_traces, 1)) / IMAGING_FS;

    % Calculate proper spacing for trace overlay
    fiber_display_count = min(num_fibers, MAX_FIBERS_DISPLAY);
    signal_ranges = zeros(fiber_display_count, 1);
    for i = 1:fiber_display_count
        signal_ranges(i) = range(processed_traces(:, i));
    end
    max_signal_range = max(signal_ranges);
    trace_spacing = max_signal_range * 1.8;
    fiber_spacing = max_signal_range * 1.3;  % Tighter spacing between fibers
    lfp_fiber_gap = max_signal_range * 0.5;  % Extra gap between LFP and first fiber

    %% FIGURE 1: Multi-Fiber Overview (Full Recording) with mPFC logic
    fig_overview = figure('Name', 'Multi_Fiber_Overview', ...
        'Color', 'w', 'Position', [50, 50, 1400, 800]);

    % Add stimulation period shading if applicable (skip for baseline trials)
    % We'll set y-limits after plotting all signals based on actual data ranges
    if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
        % Use temporary y-limits for shading, will be updated after plotting
        if EPHYS_LOADED
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                temp_y_top = (fiber_display_count + 2) * trace_spacing;
            else
                temp_y_top = (fiber_display_count + 1) * trace_spacing;
            end
            temp_y_bottom = -trace_spacing;
        else
            temp_y_top = fiber_display_count * trace_spacing;
            temp_y_bottom = 0;
        end
        patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
            [temp_y_bottom temp_y_bottom temp_y_top temp_y_top], ...
            COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        hold on;
    else
        hold on;
    end

    % Plot signals based on configuration
    current_offset = fiber_display_count + 1;

    % Plot LFP signals at top (if available)
    if EPHYS_LOADED
        lfp_display = lfp_z_scored * max_signal_range * 0.2;

        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % Plot both LFP1 and mPFC
            mPFC_display = mPFC_z_scored * max_signal_range * 0.2;
            y_offset_mPFC = (current_offset + 1) * trace_spacing;
            y_offset_lfp = current_offset * trace_spacing;

            plot(tvec_fiber, mPFC_display + y_offset_mPFC, 'Color', COLOR_mPFC, ...
                'LineWidth', 1, 'DisplayName', 'mPFC LFP');
            plot(tvec_fiber, lfp_display + y_offset_lfp, 'Color', COLOR_LFP, ...
                'LineWidth', 1, 'DisplayName', 'LFP1');

            current_offset = current_offset - 1; % Adjust for next signals
        else
            % Plot only LFP1
            y_offset_lfp = current_offset * trace_spacing;
            plot(tvec_fiber, lfp_display + y_offset_lfp, 'Color', COLOR_LFP, ...
                'LineWidth', 1, 'DisplayName', 'LFP1');
        end

        current_offset = current_offset - 1; % Adjust for fibers
        % Add extra gap between LFP and first fiber to prevent overlap
        current_offset = current_offset - lfp_fiber_gap / trace_spacing;
    end

    % Plot each fiber with different colors and offsets (use tighter spacing)
    for i = 1:fiber_display_count
        fiber_signal = processed_traces(:, i) * polarity * 2.0;
        y_offset = (current_offset + 1 - i) * fiber_spacing;

        if size(fiber_colors, 1) >= i
            color = fiber_colors(i, :);
        else
            color = COLOR_FIBER;
        end

        plot(tvec_fiber, fiber_signal + y_offset, 'Color', color, ...
            'LineWidth', 1, 'DisplayName', sprintf('Fiber %d', i));
    end

    % Plot motion signal at bottom (if available)
    % Position motion using trace_spacing system, scale similar to LFP
    if EPHYS_LOADED
        motion_display = motion_z_scored * max_signal_range * 0.2; % Same scaling as LFP
        % Position motion closer to fibers (reduce whitespace)
        motion_y_offset = -fiber_spacing * 0.8; % Tighter spacing for motion
        plot(tvec_fiber, motion_display + motion_y_offset, 'Color', COLOR_MOTION, ...
            'LineWidth', 1, 'DisplayName', 'Motion');
    end
    
    % Calculate y-axis limits based on actual plotted data (tighter margins)
    all_children = get(gca, 'Children');
    all_y_data = [];
    for i = 1:length(all_children)
        if strcmp(get(all_children(i), 'Type'), 'line')
            y_data = get(all_children(i), 'YData');
            all_y_data = [all_y_data, y_data];
        end
    end
    if ~isempty(all_y_data)
        data_min = min(all_y_data);
        data_max = max(all_y_data);
        data_range = data_max - data_min;
        y_bottom = data_min - data_range * 0.02; % Reduced to 2% margin below
        y_top = data_max + data_range * 0.02; % Reduced to 2% margin above
    else
        % Fallback if no data found
        if EPHYS_LOADED
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                y_top = (fiber_display_count + 2) * trace_spacing + trace_spacing * 0.1;
            else
                y_top = (fiber_display_count + 1) * trace_spacing + trace_spacing * 0.1;
            end
            y_bottom = -fiber_spacing * 0.8 - fiber_spacing * 0.1;
        else
            y_top = fiber_display_count * fiber_spacing + fiber_spacing * 0.1;
            y_bottom = -fiber_spacing * 0.1;
        end
    end

    % Add labels on the right side
    label_positions = [];
    if EPHYS_LOADED
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % Labels for: mPFC, LFP1, Fibers, Motion
            total_signals = 3 + fiber_display_count; % mPFC + LFP1 + fibers + motion
            text(1.02, 0.92, 'mPFC', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_mPFC, 'VerticalAlignment', 'middle');
            text(1.02, 0.88, 'LFP1', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');

            % Fiber labels
            for i = 1:fiber_display_count
                y_pos = 0.84 - (i * 0.65 / (total_signals));
                if size(fiber_colors, 1) >= i
                    color = fiber_colors(i, :);
                else
                    color = COLOR_FIBER;
                end
                text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
            end

            text(1.02, 0.08, 'Motion', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');
        else
            % Labels for: LFP1, Fibers, Motion
            total_signals = 2 + fiber_display_count; % LFP1 + fibers + motion
            text(1.02, 0.9, 'LFP1', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');

            % Fiber labels
            for i = 1:fiber_display_count
                y_pos = 0.85 - (i * 0.7 / (total_signals));
                if size(fiber_colors, 1) >= i
                    color = fiber_colors(i, :);
                else
                    color = COLOR_FIBER;
                end
                text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
            end

            text(1.02, 0.1, 'Motion', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');
        end
    else
        % Only fiber labels
        for i = 1:fiber_display_count
            y_pos = 0.9 - (i * 0.8 / fiber_display_count);
            if size(fiber_colors, 1) >= i
                color = fiber_colors(i, :);
            else
                color = COLOR_FIBER;
            end
            text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
        end
    end

    % Add scale bar (5% ΔF/F)
    scale_value_deltaF_F = 0.05;  % 5% ΔF/F
    x_scale = max(tvec_fiber) * 0.85;
    % Position scale bar on top fiber trace
    y_scale = (fiber_display_count - 1) * fiber_spacing + fiber_spacing * 0.1;
    plot([x_scale x_scale], [y_scale y_scale + scale_value_deltaF_F], 'k-', 'LineWidth', 4);
    text(x_scale + max(tvec_fiber)*0.02, y_scale + scale_value_deltaF_F/2, ...
        '5% ΔF/F', 'FontSize', 12, 'FontWeight', 'bold');

    xlim([0, max(tvec_fiber)]);
    ylim([y_bottom, y_top]); % Set y-limits to match calculated range
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Multi-Fiber Recording Overview (%d fibers, ΔF/F)', num_fibers), ...
        'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'YTick', []);
    set(gca, 'FontSize', 14);
    grid off; box on;

    % Save overview figure
    fig_overview_filename = sprintf('Multi_Fiber_Overview_%s', session_folder_name);
    saveas(fig_overview, fullfile(save_directory, [fig_overview_filename '.fig']));
    saveas(fig_overview, fullfile(save_directory, [fig_overview_filename '.png']));
    fprintf('Multi-fiber overview saved: %s\n', fig_overview_filename);

    %% FIGURE 2: Multi-Fiber Zoomed View with mPFC logic
    fig_zoomed = figure('Name', 'Multi_Fiber_Zoomed', ...
        'Color', 'w', 'Position', [100, 100, 1400, 800]);

    % Define zoom window
    if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
        zoom_start = STIM_PERIOD(1) - 5;
        zoom_end = min(STIM_PERIOD(2) + 15, max(tvec_fiber));
    else
        % For baseline trials or multi-FOV, use middle portion
        zoom_start = max(tvec_fiber) * 0.4;
        zoom_end = max(tvec_fiber) * 0.6;
    end
    zoom_idx = tvec_fiber >= zoom_start & tvec_fiber <= zoom_end;

    % Add stimulation shading for zoom (skip for baseline trials)
    % We'll set y-limits after plotting all signals based on actual data ranges
    if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
        % Use temporary y-limits for shading, will be updated after plotting
        if EPHYS_LOADED
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                temp_y_top_zoom = (fiber_display_count + 2) * trace_spacing;
            else
                temp_y_top_zoom = (fiber_display_count + 1) * trace_spacing;
            end
            temp_y_bottom_zoom = -trace_spacing;
        else
            temp_y_top_zoom = fiber_display_count * trace_spacing;
            temp_y_bottom_zoom = 0;
        end
        patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
            [temp_y_bottom_zoom temp_y_bottom_zoom temp_y_top_zoom temp_y_top_zoom], ...
            COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
        hold on;
    else
        hold on;
    end

    % Plot zoomed signals with same logic as overview
    current_offset = fiber_display_count + 1;

    if EPHYS_LOADED
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % Plot both LFP signals
            mPFC_display = mPFC_z_scored * max_signal_range * 0.2;
            y_offset_mPFC = (current_offset + 1) * trace_spacing;
            y_offset_lfp = current_offset * trace_spacing;

            plot(tvec_fiber(zoom_idx), mPFC_display(zoom_idx) + y_offset_mPFC, ...
                'Color', COLOR_mPFC, 'LineWidth', 1);
            plot(tvec_fiber(zoom_idx), lfp_display(zoom_idx) + y_offset_lfp, ...
                'Color', COLOR_LFP, 'LineWidth', 1);

            current_offset = current_offset - 1;
        else
            % Plot only LFP1
            y_offset_lfp = current_offset * trace_spacing;
            plot(tvec_fiber(zoom_idx), lfp_display(zoom_idx) + y_offset_lfp, ...
                'Color', COLOR_LFP, 'LineWidth', 1);
        end

        current_offset = current_offset - 1;
        % Add extra gap between LFP and first fiber to prevent overlap
        current_offset = current_offset - lfp_fiber_gap / trace_spacing;
    end

    % Plot fibers (use tighter spacing)
    for i = 1:fiber_display_count
        fiber_signal = processed_traces(:, i) * polarity * 2.0;
        y_offset = (current_offset + 1 - i) * fiber_spacing;

        if size(fiber_colors, 1) >= i
            color = fiber_colors(i, :);
        else
            color = COLOR_FIBER;
        end

        plot(tvec_fiber(zoom_idx), fiber_signal(zoom_idx) + y_offset, ...
            'Color', color, 'LineWidth', 1);
    end

    % Plot motion (using same y_offset and scaling as overview)
    if EPHYS_LOADED
        motion_display_zoom = motion_z_scored * max_signal_range * 0.2; % Same scaling as LFP
        motion_y_offset = -fiber_spacing * 0.8; % Tighter spacing for motion (same as overview)
        plot(tvec_fiber(zoom_idx), motion_display_zoom(zoom_idx) + motion_y_offset, ...
            'Color', COLOR_MOTION, 'LineWidth', 1);
    end
    
    % Calculate y-axis limits for zoomed view based on actual plotted data (tighter margins)
    all_children_zoom = get(gca, 'Children');
    all_y_data_zoom = [];
    for i = 1:length(all_children_zoom)
        if strcmp(get(all_children_zoom(i), 'Type'), 'line')
            y_data_zoom = get(all_children_zoom(i), 'YData');
            all_y_data_zoom = [all_y_data_zoom, y_data_zoom];
        end
    end
    if ~isempty(all_y_data_zoom)
        data_min_zoom = min(all_y_data_zoom);
        data_max_zoom = max(all_y_data_zoom);
        data_range_zoom = data_max_zoom - data_min_zoom;
        y_bottom_zoom = data_min_zoom - data_range_zoom * 0.02; % Reduced to 2% margin below
        y_top_zoom = data_max_zoom + data_range_zoom * 0.02; % Reduced to 2% margin above
    else
        % Fallback if no data found
        if EPHYS_LOADED
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                y_top_zoom = (fiber_display_count + 2) * trace_spacing + trace_spacing * 0.1;
            else
                y_top_zoom = (fiber_display_count + 1) * trace_spacing + trace_spacing * 0.1;
            end
            y_bottom_zoom = -fiber_spacing * 0.8 - fiber_spacing * 0.1;
        else
            y_top_zoom = fiber_display_count * fiber_spacing + fiber_spacing * 0.1;
            y_bottom_zoom = -fiber_spacing * 0.1;
        end
    end

    % Add labels (same as overview)
    if EPHYS_LOADED
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            text(1.02, 0.92, 'mPFC', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_mPFC, 'VerticalAlignment', 'middle');
            text(1.02, 0.88, 'LFP1', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');

            for i = 1:fiber_display_count
                y_pos = 0.84 - (i * 0.65 / (3 + fiber_display_count));
                if size(fiber_colors, 1) >= i
                    color = fiber_colors(i, :);
                else
                    color = COLOR_FIBER;
                end
                text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
            end

            text(1.02, 0.08, 'Motion', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');
        else
            text(1.02, 0.9, 'LFP1', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');

            for i = 1:fiber_display_count
                y_pos = 0.85 - (i * 0.7 / (2 + fiber_display_count));
                if size(fiber_colors, 1) >= i
                    color = fiber_colors(i, :);
                else
                    color = COLOR_FIBER;
                end
                text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                    'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
            end

            text(1.02, 0.1, 'Motion', 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');
        end
    else
        for i = 1:fiber_display_count
            y_pos = 0.9 - (i * 0.8 / fiber_display_count);
            if size(fiber_colors, 1) >= i
                color = fiber_colors(i, :);
            else
                color = COLOR_FIBER;
            end
            text(1.02, y_pos, sprintf('Fiber %d', i), 'Units', 'normalized', 'FontSize', 14, ...
                'FontWeight', 'bold', 'Color', color, 'VerticalAlignment', 'middle');
        end
    end

    xlim([zoom_start, zoom_end]);
    ylim([y_bottom_zoom, y_top_zoom]); % Set y-limits to match calculated range
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Multi-Fiber Zoomed View (%.1f-%.1f s, ΔF/F)', zoom_start, zoom_end), ...
        'FontSize', 18, 'FontWeight', 'bold');
    set(gca, 'YTick', []);
    set(gca, 'FontSize', 14);
    grid off; box on;

    % Save zoomed figure
    fig_zoomed_filename = sprintf('Multi_Fiber_Zoomed_%s', session_folder_name);
    saveas(fig_zoomed, fullfile(save_directory, [fig_zoomed_filename '.fig']));
    saveas(fig_zoomed, fullfile(save_directory, [fig_zoomed_filename '.png']));
    fprintf('Multi-fiber zoomed view saved: %s\n', fig_zoomed_filename);

    fprintf('Multi-fiber visualization completed - 2 clean figures generated\n');
end

%% ========================================================================
%  ORIGINAL: INDIVIDUAL FIBER ANALYSIS (PRESERVED)
%  ========================================================================

% Process each fiber individually (original functionality preserved)
for fiber_idx = 1:num_fibers

    fprintf('Creating visualizations for fiber %d...\n', fiber_idx);

    Vx = processed_traces(:, fiber_idx);
    tvec_fiber = (1:length(Vx)) / IMAGING_FS;

    %% Create FieldTrip structure for THIS specific fiber
    if EPHYS_LOADED
        fprintf('Running FieldTrip analysis for fiber %d...\n', fiber_idx);
        lfp_ft = [];

        % Determine number of channels based on mPFC loading
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % 3-channel analysis: Fiber, LFP1, mPFC LFP
            lfp_ft.trial{1}(1,:) = Vx;                    % Fiber
            lfp_ft.trial{1}(2,:) = lfp_aligned;           % LFP1 (Ch11)
            lfp_ft.trial{1}(3,:) = lfp_aligned_mPFC;      % mPFC LFP (Ch1 - Ch3)
            lfp_ft.time{1} = tvec_fiber;
            lfp_ft.label{1} = sprintf('Fiber_%d', fiber_idx);
            lfp_ft.label{2} = 'LFP1';
            lfp_ft.label{3} = 'mPFC_LFP';
            fprintf('  Using 3-channel analysis (Fiber + LFP1 + mPFC)\n');
        else
            % 2-channel analysis: Fiber, LFP1 only
            lfp_ft.trial{1}(1,:) = Vx;                    % Fiber
            lfp_ft.trial{1}(2,:) = lfp_aligned;           % LFP1 (Ch11)
            lfp_ft.time{1} = tvec_fiber;
            lfp_ft.label{1} = sprintf('Fiber_%d', fiber_idx);
            lfp_ft.label{2} = 'LFP1';
            fprintf('  Using 2-channel analysis (Fiber + LFP1)\n');
        end

        % Configure time-frequency analysis
        cfg = [];
        cfg.method = 'mtmconvol';
        cfg.output = 'fourier';
        cfg.taper = 'hanning';
        cfg.keeptapers = 'yes';
        cfg.keeptrials = 'yes';
        cfg.tapsmofrq = 5;
        cfg.channel = 'all';
        cfg.foi = 2:1:70;
        cfg.toi = lfp_ft.time{1}(1:1:end);
        cfg.width = 8;
        cfg.t_ftimwin = ones(1, length(cfg.foi)) * 0.5;

        fprintf('Running FieldTrip time-frequency analysis...\n');
        freq_result = ft_freqanalysis(cfg, lfp_ft);
        fprintf('Analysis complete\n');

        % Extract phase information
        phase_fiber = squeeze(angle(freq_result.fourierspctrm(1, 1, :, :)));
        phase_lfp = squeeze(angle(freq_result.fourierspctrm(1, 2, :, :)));

        % Conditionally extract mPFC phase
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            phase_lfp_mPFC = squeeze(angle(freq_result.fourierspctrm(1, 3, :, :)));
            fprintf('  mPFC LFP phase extracted\n');
        end

        freq_vector = freq_result.freq';
        time_vector = freq_result.time;
        fprintf('FieldTrip analysis complete for fiber %d\n', fiber_idx);
    end

    %% Compute local spectrogram for fiber
    [s_fiber, w_fiber, t_fiber] = spectrogram(Vx, 480, 449, 1:1:150, IMAGING_FS);
    fprintf('Fiber %d: Spectrogram computed, size of s_fiber: [%d, %d]\n', fiber_idx, size(s_fiber, 1), size(s_fiber, 2));

    %% Compute band powers for fiber
    band_power_fiber = compute_band_power(s_fiber, w_fiber, BAND_RANGES);

    % Z-score and smooth fiber band powers
    band_power_fiber = zscore_smooth_bands(band_power_fiber, 20);

    if EPHYS_LOADED
        %% Compute local spectrogram for LFP
        [s_lfp, w_lfp, t_lfp] = spectrogram(lfp_aligned, 480, 449, 1:1:150, IMAGING_FS);

        %% Conditionally compute spectrogram for mPFC LFP
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            [s_lfp_mPFC, w_lfp_mPFC, t_lfp_mPFC] = spectrogram(lfp_aligned_mPFC, 480, 449, 1:1:150, IMAGING_FS);
            fprintf('mPFC LFP spectrogram computed\n');
        end

        %% Compute band powers for LFP
        band_power_lfp = compute_band_power(s_lfp, w_lfp, BAND_RANGES);

        %% Conditionally compute band powers for mPFC LFP
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            band_power_lfp_mPFC = compute_band_power(s_lfp_mPFC, w_lfp_mPFC, BAND_RANGES);
            fprintf('mPFC LFP band powers computed\n');
        end

        % Z-score and smooth LFP band powers (per band, each trace independent)
        band_power_lfp = zscore_smooth_bands(band_power_lfp, 20);
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            band_power_lfp_mPFC = zscore_smooth_bands(band_power_lfp_mPFC, 20);
        end
    end

    % Define viridis-compatible colors for band power
    viridis_band_colors = [
        0.267004, 0.004874, 0.329415;  % Dark purple - Delta-Theta
        0.253935, 0.265254, 0.529983;  % Blue - Alpha
        0.127568, 0.566949, 0.550556;  % Teal - Beta
        0.477504, 0.821444, 0.318195;  % Yellow-green - Low Gamma
        0.993248, 0.906157, 0.143936   % Yellow - High Gamma
        ];

    %% FIGURE 1: MAIN ANALYSIS - Conditional mPFC LFP
    if EPHYS_LOADED
        % Determine figure height based on number of LFPs
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            fig_height = 1600;  % 4 rows for both LFPs
            num_rows = 4;
        else
            fig_height = 1200;  % 3 rows for single LFP
            num_rows = 3;
        end

        fig1 = figure('Name', sprintf('Figure1_MainAnalysis_Fiber%d', fiber_idx), ...
            'Color', 'w', 'Position', [50, 50, 1700, fig_height]);

        %% Define color for mPFC LFP
        COLOR_mPFC = [0.8, 0.3, 0.3];  % Reddish for mPFC

        %% Subplot 1: Combined traces
        subplot(num_rows, 6, [1 3]);

        % Determine number of signals to plot
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            num_signals = 4;  % LFP1, mPFC, Fiber, Motion
        else
            num_signals = 3;  % LFP1, Fiber, Motion
        end

        % Normalize traces (use pre-computed z-scored signals)
        lfp_z = lfp_z_scored;
        fiber_z = polarity * Vx;
        motion_z = motion_z_scored;

        % Adjust spacing
        lfp_range = max(lfp_z) - min(lfp_z);
        fiber_range = max(fiber_z) - min(fiber_z);
        motion_range = max(motion_z) - min(motion_z);
        max_range = max([lfp_range, fiber_range, motion_range]);
        trace_spacing = max_range * 1.8;

        lfp_normalized = (lfp_z - mean(lfp_z)) / (max(abs(lfp_z - mean(lfp_z))) + eps) * (trace_spacing * 0.4);
        fiber_normalized = (fiber_z - mean(fiber_z)) / (max(abs(fiber_z - mean(fiber_z))) + eps) * (trace_spacing * 0.4);
        motion_normalized = (motion_z - mean(motion_z)) / (max(abs(motion_z - mean(motion_z))) + eps) * (trace_spacing * 0.4);

        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            lfp_z_mPFC = mPFC_z_scored;
            lfp_normalized_mPFC = (lfp_z_mPFC - mean(lfp_z_mPFC)) / (max(abs(lfp_z_mPFC - mean(lfp_z_mPFC))) + eps) * (trace_spacing * 0.4);
        end

        % Plot with stimulation shading (skip for baseline trials)
        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                y_bottom = -trace_spacing*0.5;
                y_top = 3*trace_spacing + trace_spacing*0.5;
            else
                y_bottom = -trace_spacing*0.5;
                y_top = 2*trace_spacing + trace_spacing*0.5;
            end

            patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                [y_bottom, y_bottom, y_top, y_top], ...
                COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            hold on;
        else
            hold on;
        end

        % Plot signals based on configuration
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % 4 signals: LFP1, mPFC, Fiber, Motion
            plot(tvec_fiber, lfp_normalized + 3*trace_spacing, 'Color', COLOR_LFP, 'LineWidth', 1.5);
            plot(tvec_fiber, lfp_normalized_mPFC + 2*trace_spacing, 'Color', COLOR_mPFC, 'LineWidth', 1.5);
            plot(tvec_fiber, fiber_normalized + trace_spacing, 'Color', COLOR_FIBER, 'LineWidth', 1.5);
            plot(tvec_fiber, motion_normalized, 'Color', COLOR_MOTION, 'LineWidth', 1.5);

            % Labels
            text(1.02, 0.88, 'LFP1', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');
            text(1.02, 0.66, 'mPFC', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_mPFC, 'VerticalAlignment', 'middle');
            text(1.02, 0.44, 'Fiber', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_FIBER, 'VerticalAlignment', 'middle');
            text(1.02, 0.22, 'Motion', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');

            ylim([-trace_spacing*0.7, 3*trace_spacing + trace_spacing*0.7]);
            title('LFP1, mPFC, Fiber, and Motion', 'FontSize', 14, 'FontWeight', 'bold');
        else
            % 3 signals: LFP1, Fiber, Motion
            plot(tvec_fiber, lfp_normalized + 2*trace_spacing, 'Color', COLOR_LFP, 'LineWidth', 1.5);
            plot(tvec_fiber, fiber_normalized + trace_spacing, 'Color', COLOR_FIBER, 'LineWidth', 1.5);
            plot(tvec_fiber, motion_normalized, 'Color', COLOR_MOTION, 'LineWidth', 1.5);

            % Labels
            text(1.02, 0.83, 'LFP1', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');
            text(1.02, 0.50, 'Fiber', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_FIBER, 'VerticalAlignment', 'middle');
            text(1.02, 0.17, 'Motion', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');

            ylim([-trace_spacing*0.7, 2*trace_spacing + trace_spacing*0.7]);
            title('LFP, Fiber, and Motion', 'FontSize', 14, 'FontWeight', 'bold');
        end

        xlim([0, max(tvec_fiber)]);
        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Normalized Signal (LFP: z-score, Fiber: ΔF/F, Motion: z-score)', 'FontSize', 13, 'FontWeight', 'bold');
        set(gca, 'YTick', []);
        grid off; box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 2: Zoomed traces - SPANS 2 COLUMNS with mPFC logic
        ax = subplot(num_rows, 6, [4 6]); % Use num_rows instead of hardcoded 3
        pos = get(ax, 'Position');
        pos(1) = pos(1) + 0.02;
        pos(3) = pos(3) - 0.02;
        set(ax, 'Position', pos);

        % Define zoom window (skip stim period for baseline trials)
        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            zoom_start = STIM_PERIOD(1)-5;
            zoom_end = min(STIM_PERIOD(2) + 40, max(tvec_fiber));
        else
            % For baseline trials or multi-FOV, use middle portion
            zoom_start = max(tvec_fiber) * 0.4;
            zoom_end = max(tvec_fiber) * 0.6;
        end
        zoom_idx = tvec_fiber >= zoom_start & tvec_fiber <= zoom_end;

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            % Determine y-limits based on mPFC configuration
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                y_bottom = -trace_spacing*0.5;
                y_top = 3*trace_spacing + trace_spacing*0.5;
            else
                y_bottom = -trace_spacing*0.5;
                y_top = 2*trace_spacing + trace_spacing*0.5;
            end

            % Add stimulation shading patches for each signal
            if LOAD_mPFC_LFP && mPFC_LFP_LOADED
                % 4 signals: LFP1, mPFC, Fiber, Motion
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [3*trace_spacing - trace_spacing*0.5, 3*trace_spacing - trace_spacing*0.5, ...
                    3*trace_spacing + trace_spacing*0.5, 3*trace_spacing + trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                hold on;
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [2*trace_spacing - trace_spacing*0.5, 2*trace_spacing - trace_spacing*0.5, ...
                    2*trace_spacing + trace_spacing*0.5, 2*trace_spacing + trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [trace_spacing - trace_spacing*0.5, trace_spacing - trace_spacing*0.5, ...
                    trace_spacing + trace_spacing*0.5, trace_spacing + trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [-trace_spacing*0.5, -trace_spacing*0.5, trace_spacing*0.5, trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            else
                % 3 signals: LFP1, Fiber, Motion
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [2*trace_spacing - trace_spacing*0.5, 2*trace_spacing - trace_spacing*0.5, ...
                    2*trace_spacing + trace_spacing*0.5, 2*trace_spacing + trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                hold on;
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [trace_spacing - trace_spacing*0.5, trace_spacing - trace_spacing*0.5, ...
                    trace_spacing + trace_spacing*0.5, trace_spacing + trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [-trace_spacing*0.5, -trace_spacing*0.5, trace_spacing*0.5, trace_spacing*0.5], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            end
        else
            hold on;
        end

        % Plot signals based on configuration
        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            % 4 signals: LFP1, mPFC, Fiber, Motion
            plot(tvec_fiber(zoom_idx), lfp_normalized(zoom_idx) + 3*trace_spacing, 'Color', COLOR_LFP, 'LineWidth', 1.5);
            plot(tvec_fiber(zoom_idx), lfp_normalized_mPFC(zoom_idx) + 2*trace_spacing, 'Color', COLOR_mPFC, 'LineWidth', 1.5);
            plot(tvec_fiber(zoom_idx), fiber_normalized(zoom_idx) + trace_spacing, 'Color', COLOR_FIBER, 'LineWidth', 1.5);
            plot(tvec_fiber(zoom_idx), motion_normalized(zoom_idx), 'Color', COLOR_MOTION, 'LineWidth', 1.5);

            % Add labels outside
            text(1.02, 0.88, 'LFP1', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');
            text(1.02, 0.66, 'mPFC', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_mPFC, 'VerticalAlignment', 'middle');
            text(1.02, 0.44, 'Fiber', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_FIBER, 'VerticalAlignment', 'middle');
            text(1.02, 0.22, 'Motion', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');

            ylim([-trace_spacing*0.7, 3*trace_spacing + trace_spacing*0.7]);
        else
            % 3 signals: LFP1, Fiber, Motion
            plot(tvec_fiber(zoom_idx), lfp_normalized(zoom_idx) + 2*trace_spacing, 'Color', COLOR_LFP, 'LineWidth', 1.5);
            plot(tvec_fiber(zoom_idx), fiber_normalized(zoom_idx) + trace_spacing, 'Color', COLOR_FIBER, 'LineWidth', 1.5);
            plot(tvec_fiber(zoom_idx), motion_normalized(zoom_idx), 'Color', COLOR_MOTION, 'LineWidth', 1.5);

            % Add labels outside
            text(1.02, 0.83, 'LFP1', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_LFP, 'VerticalAlignment', 'middle');
            text(1.02, 0.50, 'Fiber', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_FIBER, 'VerticalAlignment', 'middle');
            text(1.02, 0.17, 'Motion', 'Units', 'normalized', 'FontSize', 11, ...
                'FontWeight', 'bold', 'Color', COLOR_MOTION, 'VerticalAlignment', 'middle');

            ylim([-trace_spacing*0.7, 2*trace_spacing + trace_spacing*0.7]);
        end

        xlim([zoom_start, zoom_end]);
        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Zoomed View (%.0f-%.0f s)', zoom_start, zoom_end), 'FontSize', 14, 'FontWeight', 'bold');
        set(gca, 'YTick', []);
        grid off; box on;
        set(gca, 'LineWidth', 1.5);

        %% ROW 2: FIBER ANALYSIS

        %% Subplot 3: Fiber Spectrogram
        subplot(num_rows, 6, [7 10]);

        spec_scaled = log10(abs(s_fiber) + eps);
        imagesc(t_fiber, w_fiber, spec_scaled);

        axis xy;
        colormap(gca, viridis());
        ylabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        title('Fiber Spectrogram', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [mean(get(gca,'XLim')), 75, 0]);
        c = colorbar;
        c.Label.String = 'log_{10}(Power)';
        c.Label.FontSize = 12;
        c.Label.FontWeight = 'bold';

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            hold on;
            xline([STIM_PERIOD(1) STIM_PERIOD(2)], 'Color', [1 1 1], 'LineWidth', 2);
        end

        ylim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 4: Fiber Band Power
        subplot(num_rows, 6, 11);
        hold on;
        for b = 1:length(BAND_NAMES)
            plot(t_fiber, band_power_fiber(b, :) + 5 * (b-1), ...
                'Color', viridis_band_colors(b,:), 'LineWidth', 2);
        end

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            y_limits = ylim;  % get current y-axis limits
            patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                [0.8 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end


        ylim([-2, 5 * length(BAND_NAMES)]);
        yticks(5 * (0:length(BAND_NAMES)-1));
        yticklabels(BAND_NAMES);
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
        title('Fiber Band Power', 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 5: Fiber PSD
        subplot(num_rows, 6, 12);

        if PROCESS_SINGLE_TRIAL
            if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
                % Baseline trial: compute PSD over entire recording
                [psd_full, f_full] = pwelch(Vx, [], [], [], IMAGING_FS);
                hold on;
                loglog(f_full, psd_full, 'Color', COLOR_FIBER, 'LineWidth', 2, 'DisplayName', 'Baseline');
                legend({'Baseline'}, 'Location', 'northeast', 'FontSize', 9);
            else
                pre_stim_idx = tvec_fiber >= PRE_STIM_PERIOD(1) & tvec_fiber < PRE_STIM_PERIOD(2);
                stim_idx = tvec_fiber >= STIM_PERIOD(1) & tvec_fiber < STIM_PERIOD(2);
                post_stim_idx = tvec_fiber >= STIM_PERIOD(2);

                % Check if we have enough data for pwelch (need at least ~50 samples for reliable PSD)
                MIN_SAMPLES_FOR_PSD = 50;

                % Pre-stimulation PSD
                if sum(pre_stim_idx) >= MIN_SAMPLES_FOR_PSD
                    [psd_pre, f_pre] = pwelch(Vx(pre_stim_idx), [], [], [], IMAGING_FS);
                    plot_pre = true;
                else
                    fprintf('  Warning: Pre-stim period too short (%d samples) for PSD computation\n', sum(pre_stim_idx));
                    plot_pre = false;
                end

                % Stimulation PSD
                if sum(stim_idx) >= MIN_SAMPLES_FOR_PSD
                    [psd_stim, f_stim] = pwelch(Vx(stim_idx), [], [], [], IMAGING_FS);
                    plot_stim = true;
                else
                    fprintf('  Warning: Stim period too short (%d samples) for PSD computation\n', sum(stim_idx));
                    plot_stim = false;
                end

                % Post-stimulation PSD
                if sum(post_stim_idx) >= MIN_SAMPLES_FOR_PSD
                    [psd_post, f_post] = pwelch(Vx(post_stim_idx), [], [], [], IMAGING_FS);
                    plot_post = true;
                else
                    fprintf('  Warning: Post-stim period too short (%d samples) for PSD computation\n', sum(post_stim_idx));
                    plot_post = false;
                end

                % Plot only the PSDs that we successfully computed
                hold on;
                if plot_pre
                    loglog(f_pre, psd_pre, 'Color', viridis_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Pre');
                end
                if plot_stim
                    loglog(f_stim, psd_stim, 'Color', viridis_colors(2,:), 'LineWidth', 2, 'DisplayName', 'Stim');
                end
                if plot_post
                    loglog(f_post, psd_post, 'Color', viridis_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Post');
                end

                % Create legend only for available data
                legend_entries = {};
                if plot_pre, legend_entries{end+1} = 'Pre'; end
                if plot_stim, legend_entries{end+1} = 'Stim'; end
                if plot_post, legend_entries{end+1} = 'Post'; end

                if ~isempty(legend_entries)
                    legend(legend_entries, 'Location', 'northeast', 'FontSize', 9);
                else
                    % If no PSDs could be computed, show a message
                    text(0.5, 0.5, 'Insufficient data for PSD computation', ...
                        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                        'FontSize', 12, 'Color', 'red');
                end
            end

        else
            [psd, f] = pwelch(Vx, [], [], [], IMAGING_FS);
            loglog(f, psd, 'Color', COLOR_FIBER, 'LineWidth', 2);
        end

        xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Fiber PSD', 'FontSize', 11, 'FontWeight', 'bold');
        xlim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% ROW 3: LFP ANALYSIS

        %% Subplot 6: LFP Spectrogram
        subplot(num_rows, 6, [13 16]);

        if exist('spectrogram_lfp', 'var')
            spec_scaled = log10(spectrogram_lfp + eps);
            imagesc(time_vector, freq_vector, spec_scaled);
        else
            spec_scaled = log10(abs(s_lfp) + eps);
            imagesc(t_lfp, w_lfp, spec_scaled);
        end

        axis xy;
        colormap(gca, viridis());
        ylabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        title('LFP Spectrogram', 'FontSize', 14, 'FontWeight', 'bold', 'Position', [mean(get(gca,'XLim')), 75, 0]);
        c = colorbar;
        c.Label.String = 'log_{10}(Power)';
        c.Label.FontSize = 12;
        c.Label.FontWeight = 'bold';

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            hold on;
            xline([STIM_PERIOD(1) STIM_PERIOD(2)], 'Color', [1 1 1], 'LineWidth', 2);
        end

        ylim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 7: LFP Band Power
        subplot(num_rows, 6, 17);
        hold on;
        for b = 1:length(BAND_NAMES)
            plot(t_lfp, band_power_lfp(b, :) + 5 * (b-1), ...
                'Color', viridis_band_colors(b,:), 'LineWidth', 2);
        end

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            y_limits = ylim;  % get current y-axis limits
            patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                [0.8 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        end

        ylim([-2, 5 * length(BAND_NAMES)]);
        yticks(5 * (0:length(BAND_NAMES)-1));
        yticklabels(BAND_NAMES);
        xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
        title('LFP Band Power', 'FontSize', 12, 'FontWeight', 'bold');
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 8: LFP PSD
        subplot(num_rows, 6, 18);

        if PROCESS_SINGLE_TRIAL
            if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
                % Baseline trial: compute PSD over entire recording
                [psd_full, f_full] = pwelch(lfp_aligned, [], [], [], IMAGING_FS);
                hold on;
                loglog(f_full, psd_full, 'Color', COLOR_LFP, 'LineWidth', 2, 'DisplayName', 'Baseline');
                legend({'Baseline'}, 'Location', 'northeast', 'FontSize', 9);
            else
                pre_stim_idx = tvec_fiber >= PRE_STIM_PERIOD(1) & tvec_fiber < PRE_STIM_PERIOD(2);
                stim_idx = tvec_fiber >= STIM_PERIOD(1) & tvec_fiber < STIM_PERIOD(2);
                post_stim_idx = tvec_fiber >= STIM_PERIOD(2);

                [psd_pre, f_pre] = pwelch(lfp_aligned(pre_stim_idx), [], [], [], IMAGING_FS);
                [psd_stim, f_stim] = pwelch(lfp_aligned(stim_idx), [], [], [], IMAGING_FS);
                [psd_post, f_post] = pwelch(lfp_aligned(post_stim_idx), [], [], [], IMAGING_FS);

                hold on;
                loglog(f_pre, psd_pre, 'Color', viridis_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Pre');
                loglog(f_stim, psd_stim, 'Color', viridis_colors(2,:), 'LineWidth', 2, 'DisplayName', 'Stim');
                loglog(f_post, psd_post, 'Color', viridis_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Post');

                legend('Location', 'northeast', 'FontSize', 9);
            end
        else
            [psd, f] = pwelch(lfp_aligned, [], [], [], IMAGING_FS);
            loglog(f, psd, 'Color', COLOR_LFP, 'LineWidth', 2);
        end

        xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('LFP PSD', 'FontSize', 11, 'FontWeight', 'bold');
        xlim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% ROW 4: mPFC LFP ANALYSIS (CONDITIONAL)

        if LOAD_mPFC_LFP && mPFC_LFP_LOADED
            %% Subplot: mPFC Spectrogram
            subplot(num_rows, 6, [19 22]);
            spec_scaled_mPFC = log10(abs(s_lfp_mPFC) + eps);
            imagesc(t_lfp_mPFC, w_lfp_mPFC, spec_scaled_mPFC);
            axis xy;
            colormap(gca, viridis());
            ylabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
            xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
            title('mPFC LFP Spectrogram (Ch1 - Ch3)', 'FontSize', 14, 'FontWeight', 'bold');
            c = colorbar;
            c.Label.String = 'log_{10}(Power)';
            c.Label.FontSize = 12;
            c.Label.FontWeight = 'bold';

            if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
                hold on;
                xline([STIM_PERIOD(1) STIM_PERIOD(2)], 'Color', [1 1 1], 'LineWidth', 2);
            end

            ylim([2 70]);
            grid off; box on;
            set(gca, 'LineWidth', 1.5);

            %% Subplot: mPFC Band Power
            subplot(num_rows, 6, 23);
            hold on;
            for b = 1:length(BAND_NAMES)
                plot(t_lfp_mPFC, band_power_lfp_mPFC(b, :) + 5 * (b-1), ...
                    'Color', viridis_band_colors(b,:), 'LineWidth', 2);
            end

            if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
                y_limits = ylim;
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [y_limits(1) y_limits(1) y_limits(2) y_limits(2)], ...
                    [0.8 0.3 0.3], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            end

            ylim([-2, 5 * length(BAND_NAMES)]);
            yticks(5 * (0:length(BAND_NAMES)-1));
            yticklabels(BAND_NAMES);
            xlabel('Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
            title('mPFC Band Power', 'FontSize', 12, 'FontWeight', 'bold');
            grid off; box on;
            set(gca, 'LineWidth', 1.5);

            %% Subplot: mPFC PSD
            subplot(num_rows, 6, 24);
            if PROCESS_SINGLE_TRIAL
                if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
                    % Baseline trial: compute PSD over entire recording
                    [psd_full, f_full] = pwelch(lfp_aligned_mPFC, [], [], [], IMAGING_FS);
                    hold on;
                    loglog(f_full, psd_full, 'Color', COLOR_mPFC, 'LineWidth', 2, 'DisplayName', 'Baseline');
                    legend({'Baseline'}, 'Location', 'northeast', 'FontSize', 9);
                else
                    pre_stim_idx = tvec_fiber >= PRE_STIM_PERIOD(1) & tvec_fiber < PRE_STIM_PERIOD(2);
                    stim_idx = tvec_fiber >= STIM_PERIOD(1) & tvec_fiber < STIM_PERIOD(2);
                    post_stim_idx = tvec_fiber >= STIM_PERIOD(2);

                    [psd_pre, f_pre] = pwelch(lfp_aligned_mPFC(pre_stim_idx), [], [], [], IMAGING_FS);
                    [psd_stim, f_stim] = pwelch(lfp_aligned_mPFC(stim_idx), [], [], [], IMAGING_FS);
                    [psd_post, f_post] = pwelch(lfp_aligned_mPFC(post_stim_idx), [], [], [], IMAGING_FS);

                    hold on;
                    loglog(f_pre, psd_pre, 'Color', viridis_colors(1,:), 'LineWidth', 2, 'DisplayName', 'Pre');
                    loglog(f_stim, psd_stim, 'Color', viridis_colors(2,:), 'LineWidth', 2, 'DisplayName', 'Stim');
                    loglog(f_post, psd_post, 'Color', viridis_colors(3,:), 'LineWidth', 2, 'DisplayName', 'Post');
                    legend('Location', 'northeast', 'FontSize', 9);
                end
            else
                [psd, f] = pwelch(lfp_aligned_mPFC, [], [], [], IMAGING_FS);
                loglog(f, psd, 'Color', COLOR_mPFC, 'LineWidth', 2);
            end

            xlabel('Frequency (Hz)', 'FontSize', 11, 'FontWeight', 'bold');
            ylabel('mPFC PSD', 'FontSize', 11, 'FontWeight', 'bold');
            xlim([2 70]);
            grid off; box on;
            set(gca, 'LineWidth', 1.5);
        end

        % Save Figure 1
        fig1_filename = sprintf('Figure1_MainAnalysis_Fiber%d_%s', fiber_idx, session_folder_name);
        saveas(fig1, fullfile(save_directory, [fig1_filename '.fig']));
        saveas(fig1, fullfile(save_directory, [fig1_filename '.png']));
        fprintf('  Figure 1 saved: %s\n', fig1_filename);
    end

    %% ========================================================================
    %  FIGURE 2: PHOTOBLEACHING CORRECTION ANALYSIS
    %% ========================================================================

    if APPLY_PHOTOBLEACHING_CORRECTION
        fig2 = figure('Name', sprintf('Figure2_Photobleaching_Fiber%d', fiber_idx), ...
            'Color', 'w', 'Position', [100, 100, 1500, 500]);

        %% Subplot 1: Correction Methods Comparison
        subplot(1, 3, 1);

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            % Add shaded box for stimulation
            patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                [min(get(gca,'YLim')), min(get(gca,'YLim')), max(get(gca,'YLim')), max(get(gca,'YLim'))], ...
                COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            hold on;
        else
            hold on;
        end

        % Compute ΔF/F for linear method (exponential method already computed in traces_deltaF_F)
        % Use existing baseline window from main analysis to ensure consistency
        F0_lin = mean(traces_detrended(baseline_start_frame:baseline_end_frame, fiber_idx));
        deltaF_F_lin = (traces_detrended(:, fiber_idx) - F0_lin) / F0_lin;
        
        % Use existing exponential ΔF/F (already computed in traces_deltaF_F)
        deltaF_F_exp = traces_deltaF_F(:, fiber_idx);
        
        h_lin = plot(tvec_fiber, deltaF_F_lin, ...
            'Color', COLOR_LINEAR, 'LineWidth', 2, 'DisplayName', 'Linear');
        hold on;
        h_exp = plot(tvec_fiber, deltaF_F_exp, ...
            'Color', COLOR_EXPONENTIAL, 'LineWidth', 2, 'DisplayName', 'Exponential');

        title('Correction Methods (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('ΔF/F', 'FontSize', 13, 'FontWeight', 'bold');
        legend([h_lin, h_exp], {'Linear Fit', 'Exponential Fit'}, ...
            'Location', 'best', 'FontSize', 11, 'Box', 'off');
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 2: Exponential Fit Quality
        subplot(1, 3, 2);
        if exist('fitted_params', 'var') && exist('time_pre_stim', 'var')
            plot(time_pre_stim, pre_stim_trace, 'k-', 'LineWidth', 2.5, 'DisplayName', 'Raw');
            hold on;
            fitted_baseline = double_exp_function(fitted_params, time_pre_stim);
            plot(time_pre_stim, fitted_baseline, 'Color', COLOR_EXPONENTIAL, ...
                'LineWidth', 2.5, 'LineStyle', '--', 'DisplayName', 'Fit');

            title('Exponential Fit Quality', 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
            ylabel('Fluorescence', 'FontSize', 13, 'FontWeight', 'bold');
            legend('Location', 'best', 'FontSize', 11);
        else
            text(0.5, 0.5, 'Fit data not available', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 12);
        end
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 3: Photobleaching Rate + Cumulative Change
        subplot(1, 3, 3);

        window_size = round(IMAGING_FS * 10);  % 5-second windows
        if window_size < length(filtered_traces(:, fiber_idx))
            raw_trace = filtered_traces(:, fiber_idx);

            % Compute running mean (bleaching trend)
            bleach_trend = movmean(raw_trace, window_size);

            % Use the first few seconds as baseline reference
            baseline_duration = 10;  % seconds
            baseline_samples = round(baseline_duration * IMAGING_FS);
            baseline_value = mean(raw_trace(1:min(baseline_samples, length(raw_trace))));

            % Calculate instantaneous rate in %/s
            bleach_rate_abs = [0; diff(bleach_trend)] * IMAGING_FS;
            bleach_rate_percent = (bleach_rate_abs / baseline_value) * 100;
            bleach_rate_smooth = smooth(bleach_rate_percent, round(IMAGING_FS));

            % Calculate cumulative change (total photobleaching over time)
            cumulative_change = ((bleach_trend - baseline_value) / baseline_value) * 100;

            % Calculate final total photobleaching
            total_photobleaching_percent = cumulative_change(end);

            % Add stimulation shading BEFORE plotting (skip for baseline trials)
            if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
                % Get y-limits for both axes to create proper shading
                rate_range = [min(bleach_rate_smooth), max(bleach_rate_smooth)];
                cumul_range = [min(cumulative_change), max(cumulative_change)];

                % Create shading patch (will be behind both plots)
                patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                    [rate_range(1)*1.2, rate_range(1)*1.2, rate_range(2)*1.2, rate_range(2)*1.2], ...
                    COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.6);
                hold on;
            end

            % LEFT Y-AXIS: Instantaneous Rate (%/s)
            yyaxis left
            plot(tvec_fiber, bleach_rate_smooth, 'Color',[0.2, 0.4, 0.8], 'LineWidth', 2.5);
            yline(0, '--', 'Color',[0.2, 0.4, 0.8], 'LineWidth', 1, 'Alpha', 0.7);
            ylabel('Instantaneous Rate (%/s)', 'FontSize', 13, 'FontWeight', 'bold', 'Color',[0.2, 0.4, 0.8]);
            set(gca, 'YColor', [0.2, 0.4, 0.8]);

            % RIGHT Y-AXIS: Cumulative Change (%)
            yyaxis right
            plot(tvec_fiber, cumulative_change, 'Color',[0.8, 0.3, 0.2],  'LineWidth', 2.5);
            ylabel('Cumulative Change (%)', 'FontSize', 13, 'FontWeight', 'bold', 'Color',[0.8, 0.3, 0.2]);
            set(gca, 'YColor', [0.8, 0.3, 0.2]);

            % Add zero line for cumulative (if it crosses zero)
            if min(cumulative_change) < 0 && max(cumulative_change) > 0
                yline(0, '--', 'Color',[0.8, 0.3, 0.2], 'LineWidth', 1, 'Alpha', 0.7);
            end

            % Formatting
            title('Photobleaching: Rate & Cumulative', 'FontSize', 14, 'FontWeight', 'bold');
            xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');

            % Add legend
            % legend({'Rate (%/s)', 'Cumulative (%)'}, 'Location', 'best', 'FontSize', 10);

            % Add information text box
            info_text = sprintf(['Rate: d(trend)/dt\n' ...
                'Total: %.2f%%\n' ...
                'Mean rate: %.3f %%/s'], ...
                total_photobleaching_percent, mean(bleach_rate_smooth));

            text(0.02, 0.98, info_text, 'Units', 'normalized', ...
                'FontSize', 10, 'VerticalAlignment', 'top', ...
                'BackgroundColor', [1 1 1 0.85], 'EdgeColor', [0.3 0.3 0.3], ...
                'Margin', 4);

        else
            text(0.5, 0.5, 'Insufficient data for rate calculation', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 12);
        end

        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        % Save Figure 2
        fig2_filename = sprintf('Figure2_Photobleaching_Fiber%d_%s', fiber_idx, session_folder_name);
        saveas(fig2, fullfile(save_directory, [fig2_filename '.fig']));
        saveas(fig2, fullfile(save_directory, [fig2_filename '.png']));
        fprintf('  Figure 2 saved: %s\n', fig2_filename);
    end

    %% ========================================================================
    %  FIGURE 3: PHASE-LOCKING AND CORRELATION ANALYSIS
    %% ========================================================================

    if EPHYS_LOADED
        fig3 = figure('Name', sprintf('Figure3_PhaseLocking_Fiber%d', fiber_idx), ...
            'Color', 'w', 'Position', [150, 150, 1400, 900]);

        % Compute phase-locking using FieldTrip if available
        if exist('freq_result', 'var')
            phase_fiber = squeeze(angle(freq_result.fourierspctrm(1, 1, :, :)));
            phase_lfp = squeeze(angle(freq_result.fourierspctrm(1, 2, :, :)));
            freq_vector = freq_result.freq';
            time_vector = freq_result.time;
        end

        %% Subplot 1: Overall Phase-Locking
        subplot(2, 2, 1);
        phase_difference = circ_dist(phase_fiber, phase_lfp);
        plv_all = compute_plv(phase_difference);

        plot(freq_vector, plv_all, 'Color', COLOR_FIBER, 'LineWidth', 2.5);
        axis tight;
        xlabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Phase Locking Value', 'FontSize', 13, 'FontWeight', 'bold');
        title('LFP-Fiber Phase Locking', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 2: Phase-Locking by Behavior
        subplot(2, 2, 2);

        running_threshold = 0.1;
        running_periods = running_velocity_smooth > running_threshold;
        rest_periods = running_velocity_smooth < running_threshold;

        plv_running = compute_plv(phase_difference(:, running_periods));
        plv_rest = compute_plv(phase_difference(:, rest_periods));

        plot(freq_vector, plv_running, 'Color', COLOR_RUNNING, 'LineWidth', 2.5, 'DisplayName', 'Running');
        hold on;
        plot(freq_vector, plv_rest, 'Color', COLOR_REST, 'LineWidth', 2.5, 'DisplayName', 'Rest');

        axis tight;
        legend('Location', 'best', 'FontSize', 11);
        xlabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Phase Locking Value', 'FontSize', 13, 'FontWeight', 'bold');
        title('Phase Locking: Run vs Rest', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([2 70]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 3: Theta Envelope Correlation - Time Series
        subplot(2, 2, 3);

        % Compute theta envelope correlation
        lowCutoff = 5;
        highCutoff = 10;
        filterOrder = 3;

        voltage_signal = processed_traces(:, fiber_idx);
        tvec = (1:length(voltage_signal)) / IMAGING_FS;

        % Theta-band LFP envelope (bandpass + Hilbert + smooth). See core/compute_band_envelope.m.
        [lfp_filtered_theta, lfp_envelope, lfp_envelope_smooth] = ...
            compute_band_envelope(lfp_aligned, IMAGING_FS, [lowCutoff, highCutoff], filterOrder, 90);
        lfp_envelope_z_scored = zscore(lfp_envelope_smooth);  % Compute once and reuse

        fiber_smooth = fastsmooth(voltage_signal, 300, 1, 1);

        % Compute correlation: LFP envelope (z-scored) vs Fiber (ΔF/F)
        % Correlation is scale-invariant, so we can use z-scored LFP and ΔF/F fiber
        [correlation_matrix, pval] = corrcoef(lfp_envelope_z_scored, fiber_smooth);
        theta_envelope_correlation = correlation_matrix(1, 2);

        if PROCESS_SINGLE_TRIAL && ~IS_BASELINE_TRIAL && ~isempty(STIM_PERIOD)
            % Add shaded box for stimulation
            y_limits_temp = [min([lfp_envelope_z_scored; fiber_smooth]), ...
                max([lfp_envelope_z_scored; fiber_smooth])];
            patch([STIM_PERIOD(1) STIM_PERIOD(2) STIM_PERIOD(2) STIM_PERIOD(1)], ...
                [y_limits_temp(1)*1.1, y_limits_temp(1)*1.1, y_limits_temp(2)*1.1, y_limits_temp(2)*1.1], ...
                COLOR_STIM_SHADE, 'EdgeColor', 'none', 'FaceAlpha', 0.8);
            hold on;
        else
            hold on;
        end

        % Plot: LFP envelope (z-scored) and Fiber (ΔF/F)
        h1 = plot(tvec, lfp_envelope_z_scored, 'Color', COLOR_LFP, 'LineWidth', 2);
        hold on;
        h2 = plot(tvec, fiber_smooth, 'Color', COLOR_FIBER, 'LineWidth', 2);

        xlabel('Time (s)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Signal (LFP: z-score, Fiber: ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Theta Envelope (r = %.3f)', theta_envelope_correlation), 'FontSize', 14, 'FontWeight', 'bold');

        % Explicit legend using plot handles (prevents "data1")
        legend([h1, h2], {'LFP Theta Envelope', 'Fiber Signal'}, 'Location', 'best', 'FontSize', 11);

        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        %% Subplot 4: Theta Envelope Correlation - Scatter
        subplot(2, 2, 4);
        scatter(lfp_envelope_z_scored, fiber_smooth, 15, COLOR_FIBER, ...
            'filled', 'MarkerFaceAlpha', 0.4);
        hold on;

        % Add regression line
        p = polyfit(lfp_envelope_z_scored, fiber_smooth, 1);
        x_fit = linspace(min(lfp_envelope_z_scored), max(lfp_envelope_z_scored), 100);
        y_fit = polyval(p, x_fit);
        plot(x_fit, y_fit, 'Color', [0 0 0], 'LineWidth', 2.5);

        xlabel('LFP Theta Envelope (z-score)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Fiber Signal (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
        title(sprintf('Correlation: r = %.3f, p = %.4f', theta_envelope_correlation, pval(1,2)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        text(0.05, 0.95, 'Black line: Linear fit', 'Units', 'normalized', ...
            'FontSize', 10, 'VerticalAlignment', 'top', 'Color', [0.5 0.5 0.5]);
        grid off;
        box on;
        set(gca, 'LineWidth', 1.5);

        % Save Figure 3
        fig3_filename = sprintf('Figure3_PhaseLocking_Fiber%d_%s', fiber_idx, session_folder_name);
        saveas(fig3, fullfile(save_directory, [fig3_filename '.fig']));
        saveas(fig3, fullfile(save_directory, [fig3_filename '.png']));
        fprintf('  Figure 3 saved: %s\n', fig3_filename);
    end
end

fprintf('All visualizations complete\n');

%% Store theta envelope correlation for summary
if EPHYS_LOADED && exist('theta_envelope_correlation', 'var')
    fprintf('\nTheta envelope correlation: r = %.3f\n', theta_envelope_correlation);
end

%% ========================================================================
%  FIGURE 4: CONDITIONAL HP-mPFC LFP Coherence and Cross-Correlation
%% ========================================================================

if EPHYS_LOADED && LOAD_mPFC_LFP && mPFC_LFP_LOADED
    fig_lfp_comparison = figure('Name', sprintf('LFP1_mPFC_Comparison_Fiber%d', fiber_idx), ...
        'Color', 'w', 'Position', [200, 200, 1400, 600]);

    %% Subplot 1: Coherence
    subplot(1, 3, 1);
    [coh, f_coh] = mscohere(lfp_aligned, lfp_aligned_mPFC, [], [], [], IMAGING_FS);
    plot(f_coh, coh, 'Color', [0.4, 0.2, 0.6], 'LineWidth', 2.5);
    xlabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Coherence', 'FontSize', 13, 'FontWeight', 'bold');
    title('LFP1-mPFC Coherence', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([2 70]);
    ylim([0 1]);
    grid on; box on;

    %% Subplot 2: Cross-Correlation
    subplot(1, 3, 2);
    [xcorr_vals, lags] = xcorr(zscore(lfp_aligned), zscore(lfp_aligned_mPFC), 'coeff');
    lags_time = lags / IMAGING_FS;
    plot(lags_time, xcorr_vals, 'Color', [0.4, 0.2, 0.6], 'LineWidth', 2);
    xlabel('Time Lag (s)', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Cross-Correlation', 'FontSize', 13, 'FontWeight', 'bold');
    title('LFP1-mPFC Cross-Correlation', 'FontSize', 14, 'FontWeight', 'bold');
    xlim([-1 1]);
    grid on; box on;

    %% Subplot 3: Phase Locking
    subplot(1, 3, 3);
    if exist('phase_lfp', 'var') && exist('phase_lfp_mPFC', 'var')
        phase_diff_lfps = circ_dist(phase_lfp, phase_lfp_mPFC);
        plv_lfps = compute_plv(phase_diff_lfps);
        plot(freq_vector, plv_lfps, 'Color', [0.4, 0.2, 0.6], 'LineWidth', 2.5);
        xlabel('Frequency (Hz)', 'FontSize', 13, 'FontWeight', 'bold');
        ylabel('Phase Locking Value', 'FontSize', 13, 'FontWeight', 'bold');
        title('LFP1-mPFC Phase Locking', 'FontSize', 14, 'FontWeight', 'bold');
        xlim([2 70]);
        grid on; box on;
    end

    % Save figure
    fig_lfp_comp_filename = sprintf('LFP1_mPFC_Comparison_Fiber%d_%s', fiber_idx, session_folder_name);
    saveas(fig_lfp_comparison, fullfile(save_directory, [fig_lfp_comp_filename '.fig']));
    saveas(fig_lfp_comparison, fullfile(save_directory, [fig_lfp_comp_filename '.png']));
    fprintf('  LFP1-mPFC comparison figure saved: %s\n', fig_lfp_comp_filename);
end

%% ============================================================================
%  SECTION 9: SUMMARY STATISTICS
%  ============================================================================

fprintf('\n=== SUMMARY STATISTICS ===\n');
fprintf('Number of fibers analyzed: %d\n', size(processed_traces, 2));
fprintf('Recording duration: %.2f seconds\n', length(processed_traces)/IMAGING_FS);
fprintf('Sampling rate: %.2f Hz\n', IMAGING_FS);

if PROCESS_SINGLE_TRIAL
    for fiber_idx = 1:size(processed_traces, 2)
        Vx = processed_traces(:, fiber_idx);
        tvec_fiber = (1:length(Vx)) / IMAGING_FS;

        if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
            % Baseline trial: split into early, middle, late periods
            total_time = max(tvec_fiber);
            pre_stim_idx = tvec_fiber < total_time/3;
            stim_idx = tvec_fiber >= total_time/3 & tvec_fiber < 2*total_time/3;
            post_stim_idx = tvec_fiber >= 2*total_time/3;
            
            pre_rms = rms(Vx(pre_stim_idx));
            stim_rms = rms(Vx(stim_idx));
            post_rms = rms(Vx(post_stim_idx));

            fprintf('\nFiber %d (Baseline Trial):\n', fiber_idx);
            fprintf('  Early period RMS: %.4f\n', pre_rms);
            fprintf('  Middle period RMS: %.4f\n', stim_rms);
            fprintf('  Late period RMS: %.4f\n', post_rms);
            fprintf('  Middle/Early ratio: %.4f\n', stim_rms/pre_rms);
            fprintf('  Late/Early ratio: %.4f\n', post_rms/pre_rms);
        else
            pre_stim_idx = tvec_fiber >= PRE_STIM_PERIOD(1) & tvec_fiber < PRE_STIM_PERIOD(2);
            stim_idx = tvec_fiber >= STIM_PERIOD(1) & tvec_fiber < STIM_PERIOD(2);
            post_stim_idx = tvec_fiber >= STIM_PERIOD(2);

            pre_rms = rms(Vx(pre_stim_idx));
            stim_rms = rms(Vx(stim_idx));
            post_rms = rms(Vx(post_stim_idx));

            fprintf('\nFiber %d:\n', fiber_idx);
            fprintf('  Pre-stim RMS: %.4f\n', pre_rms);
            fprintf('  Stimulation RMS: %.4f\n', stim_rms);
            fprintf('  Post-stim RMS: %.4f\n', post_rms);
            fprintf('  Stim/Pre ratio: %.4f\n', stim_rms/pre_rms);
            fprintf('  Post/Pre ratio: %.4f\n', post_rms/pre_rms);
        end
    end
end

%% ============================================================================
%  SECTION 10: SAVE COMPREHENSIVE ANALYSIS DATA
%  ============================================================================

fprintf('\n=== SAVING ANALYSIS RESULTS ===\n');
% Create comprehensive data structure
FiberPhotometryAnalysis = struct();

%% METADATA
FiberPhotometryAnalysis.metadata = METADATA;
FiberPhotometryAnalysis.metadata.analysis_date = datestr(now);
FiberPhotometryAnalysis.metadata.analysis_mode = ANALYSIS_MODE;
FiberPhotometryAnalysis.metadata.base_folder = base_folder;
FiberPhotometryAnalysis.metadata.save_directory = save_directory;
if PROCESS_SINGLE_TRIAL
    FiberPhotometryAnalysis.metadata.trial_name = trial_name;
    FiberPhotometryAnalysis.metadata.trial_folder = trial_folder;
end

%% EXPERIMENTAL PARAMETERS
FiberPhotometryAnalysis.parameters.sampling_rate = IMAGING_FS;
FiberPhotometryAnalysis.parameters.ephys_sampling_rate = EPHYS_FS;
FiberPhotometryAnalysis.parameters.num_frames = size(processed_traces, 1);
FiberPhotometryAnalysis.parameters.num_fibers = size(processed_traces, 2);
FiberPhotometryAnalysis.parameters.recording_duration_sec = size(processed_traces, 1) / IMAGING_FS;
FiberPhotometryAnalysis.parameters.motion_correction = MOTION_CORRECTION;
FiberPhotometryAnalysis.parameters.correction_type = CORRECTION_TYPE;
FiberPhotometryAnalysis.parameters.process_full_field = PROCESS_FULL_FIELD;
FiberPhotometryAnalysis.parameters.invert_trace = INVERT_TRACE;
FiberPhotometryAnalysis.parameters.ephys_loaded = EPHYS_LOADED;
FiberPhotometryAnalysis.parameters.polarity = polarity;
FiberPhotometryAnalysis.parameters.is_baseline_trial = IS_BASELINE_TRIAL;

%% TIME PERIOD CONFIGURATION
if PROCESS_SINGLE_TRIAL
    FiberPhotometryAnalysis.time_periods.is_baseline_trial = IS_BASELINE_TRIAL;
    if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
        FiberPhotometryAnalysis.time_periods.pre_stim_period = [];  % Not applicable
        FiberPhotometryAnalysis.time_periods.stim_period = [];  % Not applicable
        FiberPhotometryAnalysis.time_periods.stim_onset_frame = [];
        FiberPhotometryAnalysis.time_periods.stim_offset_frame = [];
        % For baseline trials, record the period splits used
        total_time = size(processed_traces, 1) / IMAGING_FS;
        FiberPhotometryAnalysis.time_periods.early_period = [0, total_time/3];
        FiberPhotometryAnalysis.time_periods.middle_period = [total_time/3, 2*total_time/3];
        FiberPhotometryAnalysis.time_periods.late_period = [2*total_time/3, total_time];
    else
        FiberPhotometryAnalysis.time_periods.pre_stim_period = PRE_STIM_PERIOD;
        FiberPhotometryAnalysis.time_periods.stim_period = STIM_PERIOD;
        FiberPhotometryAnalysis.time_periods.stim_onset_frame = stim_onset_frame;
        FiberPhotometryAnalysis.time_periods.stim_offset_frame = stim_offset_frame;
    end

    % Add stimulation detection metadata
    if exist('STIM_DETECTION', 'var')
        FiberPhotometryAnalysis.time_periods.detection_method = STIM_DETECTION.method;
        FiberPhotometryAnalysis.time_periods.detection_metadata = STIM_DETECTION;
        if exist('STIM_PERIOD_HARDCODED', 'var')
            FiberPhotometryAnalysis.time_periods.hardcoded_period = STIM_PERIOD_HARDCODED;
        end
    end
end

%% TIME VECTORS
FiberPhotometryAnalysis.time.time_vector_seconds = (1:size(processed_traces, 1))' / IMAGING_FS;
FiberPhotometryAnalysis.time.sampling_rate = IMAGING_FS;

%% SIGNALS - COMPREHENSIVE FIBER TRACE STORAGE
% Store all processing stages for complete traceability
FiberPhotometryAnalysis.signals.raw_traces = all_traces;                    % Raw from ROI extraction
FiberPhotometryAnalysis.signals.filtered_traces = filtered_traces;           % After artifact removal (120-124 Hz, 130-132 Hz)
FiberPhotometryAnalysis.signals.corrected_traces = traces_exp_corrected;    % After photobleaching correction (exponential method)
FiberPhotometryAnalysis.signals.deltaF_F_traces = traces_deltaF_F;         % ΔF/F computed using static baseline method
FiberPhotometryAnalysis.signals.zscored_traces = processed_traces_zscored;  % Z-scored version (for backward compatibility)
FiberPhotometryAnalysis.signals.final_processed_traces = processed_traces;  % Currently set to ΔF/F for visualization

% Store F0 values and baseline window used for ΔF/F computation
if exist('baseline_start_frame', 'var') && exist('baseline_end_frame', 'var')
    FiberPhotometryAnalysis.signals.deltaF_F_baseline_window = [baseline_start_frame, baseline_end_frame];
    FiberPhotometryAnalysis.signals.deltaF_F_baseline_time = ...
        [(baseline_start_frame-1)/IMAGING_FS, (baseline_end_frame-1)/IMAGING_FS];
    FiberPhotometryAnalysis.signals.deltaF_F_method = 'static_baseline';
    if exist('F0_values', 'var')
        FiberPhotometryAnalysis.signals.F0_values = F0_values;  % F0 value for each fiber
    end
end

if APPLY_PHOTOBLEACHING_CORRECTION
    FiberPhotometryAnalysis.signals.detrended_traces = traces_detrended;    % Linear detrending method
    FiberPhotometryAnalysis.signals.exp_corrected_traces = traces_exp_corrected;  % Exponential correction method
end

%% EPHYS SIGNALS - COMPREHENSIVE LFP STORAGE
if EPHYS_LOADED
    % Primary LFP (Ch11/HP) - Store both full original and aligned versions
    if exist('lfp_data', 'var')
        FiberPhotometryAnalysis.ephys.lfp_HP_full_original = lfp_data;      % Full original HP (Ch11) data (30 kHz, not downsampled)
    end
    FiberPhotometryAnalysis.ephys.lfp_raw_aligned_HP = lfp_aligned;         % HP LFP aligned to camera triggers (30x downsampled)
    FiberPhotometryAnalysis.ephys.lfp_z_HP = zscore(lfp_aligned);            % Z-scored HP LFP
    FiberPhotometryAnalysis.ephys.lfp_sampling_rate = IMAGING_FS;            % Effective sampling rate after downsampling
    FiberPhotometryAnalysis.ephys.lfp_original_sampling_rate = EPHYS_FS;    % Original ephys sampling rate (30 kHz)
    
    % Conditionally save mPFC LFP
    if LOAD_mPFC_LFP && mPFC_LFP_LOADED
        FiberPhotometryAnalysis.ephys.mPFC_loaded = true;
        % Store only full original mPFC differential (30 kHz, not downsampled, not aligned)
        if exist('lfp_data_mPFC', 'var')
            FiberPhotometryAnalysis.ephys.lfp_mPFC_full_original = lfp_data_mPFC;  % Full original mPFC (Ch1 - Ch3, 30 kHz, not downsampled)
        end
        % Store aligned/downsampled mPFC version only
        FiberPhotometryAnalysis.ephys.lfp_raw_aligned_mPFC = lfp_aligned_mPFC;      % mPFC LFP aligned to camera triggers (30x downsampled)
        FiberPhotometryAnalysis.ephys.lfp_z_mPFC = zscore(lfp_aligned_mPFC);       % Z-scored mPFC LFP
        FiberPhotometryAnalysis.ephys.mPFC_computation = 'Ch1 - Ch3';
    else
        FiberPhotometryAnalysis.ephys.mPFC_loaded = false;
    end
    
    % Behavioral and stimulation signals
    FiberPhotometryAnalysis.ephys.running_velocity = running_velocity_aligned;
    FiberPhotometryAnalysis.ephys.running_velocity_smooth = running_velocity_smooth;
    FiberPhotometryAnalysis.ephys.stim_pulses = stim_pulses_aligned;
    FiberPhotometryAnalysis.ephys.stim_onset = stim_onset_aligned;
    FiberPhotometryAnalysis.ephys.camera_trigger_indices = camera_trigger_indices;

    if exist('theta_envelope_correlation', 'var')
        FiberPhotometryAnalysis.ephys.theta_envelope_correlation = theta_envelope_correlation;
    end
end

%% PHOTOBLEACHING CORRECTION
FiberPhotometryAnalysis.photobleaching.correction_applied = APPLY_PHOTOBLEACHING_CORRECTION;
if APPLY_PHOTOBLEACHING_CORRECTION
    FiberPhotometryAnalysis.photobleaching.methods_used = {'linear_detrend', 'double_exponential'};
    if exist('fitted_params', 'var')
        FiberPhotometryAnalysis.photobleaching.fit_params = fitted_params;
        FiberPhotometryAnalysis.photobleaching.fit_success = true;

        if exist('pre_stim_trace', 'var')
            FiberPhotometryAnalysis.photobleaching.pre_stim_trace = pre_stim_trace;
        end
        if exist('time_pre_stim', 'var')
            FiberPhotometryAnalysis.photobleaching.time_pre_stim = time_pre_stim;
        end
        if exist('time_full', 'var')
            FiberPhotometryAnalysis.photobleaching.time_full = time_full;
            double_exp_function = @(params, t) params(1) * exp(-t/params(2)) + ...
                params(3) * exp(-t/params(4)) + params(5);
            FiberPhotometryAnalysis.photobleaching.fitted_curve = double_exp_function(fitted_params, time_full);
        end
    else
        FiberPhotometryAnalysis.photobleaching.fit_success = false;
    end

    for fiber_idx = 1:size(filtered_traces, 2)
        window_size = round(IMAGING_FS * 5);
        if window_size < length(filtered_traces(:, fiber_idx))
            raw_trace = filtered_traces(:, fiber_idx);
            bleach_trend = movmean(raw_trace, window_size);
            bleach_rate = [0; diff(bleach_trend)] * IMAGING_FS;
            FiberPhotometryAnalysis.photobleaching.fiber(fiber_idx).bleaching_rate = bleach_rate;
            FiberPhotometryAnalysis.photobleaching.fiber(fiber_idx).bleach_trend = bleach_trend;
        end
    end
end

%% ROI INFORMATION
if PROCESS_FULL_FIELD
    FiberPhotometryAnalysis.rois.type = 'full_field';
else
    FiberPhotometryAnalysis.rois.type = 'manual_selection';
end
FiberPhotometryAnalysis.rois.roi_data = all_ROIs;

%% SPECTRAL ANALYSIS
FiberPhotometryAnalysis.spectral.band_names = BAND_NAMES;
FiberPhotometryAnalysis.spectral.band_ranges = BAND_RANGES;
FiberPhotometryAnalysis.spectral.band_colors = BAND_COLORS;

% Save spectral data that was already computed during visualization
for fiber_idx = 1:size(processed_traces, 2)
    if exist('s_fiber', 'var')
        FiberPhotometryAnalysis.spectral.fiber(fiber_idx).spectrogram = abs(s_fiber);
        FiberPhotometryAnalysis.spectral.fiber(fiber_idx).frequencies = w_fiber;
        FiberPhotometryAnalysis.spectral.fiber(fiber_idx).time = t_fiber;
    end
    if exist('band_power_fiber', 'var')
        FiberPhotometryAnalysis.spectral.fiber(fiber_idx).band_power = band_power_fiber;
    end
    FiberPhotometryAnalysis.spectral.fiber(fiber_idx).trace_for_plot = polarity * processed_traces(:, fiber_idx);
end

if EPHYS_LOADED
    if exist('s_lfp', 'var')
        FiberPhotometryAnalysis.spectral.lfp.spectrogram = abs(s_lfp);
        FiberPhotometryAnalysis.spectral.lfp.frequencies = w_lfp;
        FiberPhotometryAnalysis.spectral.lfp.time = t_lfp;
    end
    if exist('band_power_lfp', 'var')
        FiberPhotometryAnalysis.spectral.lfp.band_power = band_power_lfp;
    end

    if exist('freq_result', 'var')
        FiberPhotometryAnalysis.phase_locking.freq_vector = freq_result.freq';
        FiberPhotometryAnalysis.phase_locking.time_vector = freq_result.time;
        FiberPhotometryAnalysis.phase_locking.phase_fiber = phase_fiber;
        FiberPhotometryAnalysis.phase_locking.phase_lfp = phase_lfp;

        if exist('phase_difference', 'var')
            FiberPhotometryAnalysis.phase_locking.phase_difference = phase_difference;
            FiberPhotometryAnalysis.phase_locking.plv_all = plv_all;
        end

        if exist('plv_running', 'var')
            FiberPhotometryAnalysis.phase_locking.running_threshold = running_threshold;
            FiberPhotometryAnalysis.phase_locking.running_periods = running_periods;
            FiberPhotometryAnalysis.phase_locking.rest_periods = rest_periods;
            FiberPhotometryAnalysis.phase_locking.plv_running = plv_running;
            FiberPhotometryAnalysis.phase_locking.plv_rest = plv_rest;
        end
    end

    if exist('lfp_envelope_smooth', 'var')
        FiberPhotometryAnalysis.theta_analysis.lfp_filtered_theta = lfp_filtered_theta;
        FiberPhotometryAnalysis.theta_analysis.lfp_envelope = lfp_envelope;
        FiberPhotometryAnalysis.theta_analysis.lfp_envelope_smooth = lfp_envelope_smooth;
        FiberPhotometryAnalysis.theta_analysis.filter_settings.low_cutoff = 5;
        FiberPhotometryAnalysis.theta_analysis.filter_settings.high_cutoff = 10;
        FiberPhotometryAnalysis.theta_analysis.filter_settings.filter_order = 3;

        if exist('fiber_smooth', 'var')
            FiberPhotometryAnalysis.theta_analysis.fiber(1).fiber_smooth = fiber_smooth;
            FiberPhotometryAnalysis.theta_analysis.fiber(1).correlation = theta_envelope_correlation;
            if exist('pval', 'var')
                FiberPhotometryAnalysis.theta_analysis.fiber(1).pvalue = pval(1, 2);
            end
            if exist('p', 'var')
                FiberPhotometryAnalysis.theta_analysis.fiber(1).regression_coeffs = p;
            end
        end
    end

    % Save mPFC spectral data
    if LOAD_mPFC_LFP && mPFC_LFP_LOADED
        if exist('s_lfp_mPFC', 'var')
            FiberPhotometryAnalysis.spectral.lfp_mPFC.spectrogram = abs(s_lfp_mPFC);
            FiberPhotometryAnalysis.spectral.lfp_mPFC.frequencies = w_lfp_mPFC;
            FiberPhotometryAnalysis.spectral.lfp_mPFC.time = t_lfp_mPFC;
        end

        if exist('band_power_lfp_mPFC', 'var')
            FiberPhotometryAnalysis.spectral.lfp_mPFC.band_power = band_power_lfp_mPFC;
        end

        % Save LFP1-mPFC coherence
        if exist('coh', 'var')
            FiberPhotometryAnalysis.lfp_comparison.coherence = coh;
            FiberPhotometryAnalysis.lfp_comparison.coherence_freqs = f_coh;
        end

        if exist('xcorr_vals', 'var')
            FiberPhotometryAnalysis.lfp_comparison.cross_correlation = xcorr_vals;
            FiberPhotometryAnalysis.lfp_comparison.cross_corr_lags = lags_time;
        end

        if exist('plv_lfps', 'var')
            FiberPhotometryAnalysis.lfp_comparison.phase_locking = plv_lfps;
            FiberPhotometryAnalysis.lfp_comparison.phase_locking_freqs = freq_vector;
        end
    end

end

%% PLOTTING CONFIGURATION
FiberPhotometryAnalysis.plot_config.colors.lfp = [0.2, 0.2, 0.25];
FiberPhotometryAnalysis.plot_config.colors.fiber = [0.127568, 0.566949, 0.550556];
FiberPhotometryAnalysis.plot_config.colors.motion = [0.993248, 0.7, 0.4];
FiberPhotometryAnalysis.plot_config.colors.stim_shade = [1, 0.9, 0.9];
FiberPhotometryAnalysis.plot_config.colors.linear = [0.27, 0.49, 0.77];
FiberPhotometryAnalysis.plot_config.colors.exponential = [0.83, 0.45, 0.37];
FiberPhotometryAnalysis.plot_config.colors.running = [0.163625, 0.471133, 0.558148];
FiberPhotometryAnalysis.plot_config.colors.rest = [0.741388, 0.173449, 0.149561];

viridis_band_colors = [
    0.267004, 0.004874, 0.329415;
    0.253935, 0.265254, 0.529983;
    0.127568, 0.566949, 0.550556;
    0.477504, 0.821444, 0.318195;
    0.993248, 0.906157, 0.143936
    ];
FiberPhotometryAnalysis.plot_config.band_colors = viridis_band_colors;

viridis_colors = [
    0.267004, 0.004874, 0.329415;
    0.127568, 0.566949, 0.550556;
    0.993248, 0.906157, 0.143936
    ];
FiberPhotometryAnalysis.plot_config.viridis_period_colors = viridis_colors;

%% SUMMARY STATISTICS
if PROCESS_SINGLE_TRIAL
    for fiber_idx = 1:size(processed_traces, 2)
        Vx = processed_traces(:, fiber_idx);
        tvec_fiber = (1:length(Vx))' / IMAGING_FS;

        if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
            % Baseline trial: split into early, middle, late periods
            total_time = max(tvec_fiber);
            pre_stim_idx = tvec_fiber < total_time/3;
            stim_idx = tvec_fiber >= total_time/3 & tvec_fiber < 2*total_time/3;
            post_stim_idx = tvec_fiber >= 2*total_time/3;
            
            FiberPhotometryAnalysis.summary_stats(fiber_idx).fiber_id = fiber_idx;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).is_baseline_trial = true;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).early_rms = rms(Vx(pre_stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).middle_rms = rms(Vx(stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).late_rms = rms(Vx(post_stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).middle_early_ratio = ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).middle_rms / ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).early_rms;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).late_early_ratio = ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).late_rms / ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).early_rms;
        else
            pre_stim_idx = tvec_fiber >= PRE_STIM_PERIOD(1) & tvec_fiber < PRE_STIM_PERIOD(2);
            stim_idx = tvec_fiber >= STIM_PERIOD(1) & tvec_fiber < STIM_PERIOD(2);
            post_stim_idx = tvec_fiber >= STIM_PERIOD(2);

            FiberPhotometryAnalysis.summary_stats(fiber_idx).fiber_id = fiber_idx;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).is_baseline_trial = false;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).pre_stim_rms = rms(Vx(pre_stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).stim_rms = rms(Vx(stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).post_stim_rms = rms(Vx(post_stim_idx));
            FiberPhotometryAnalysis.summary_stats(fiber_idx).stim_pre_ratio = ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).stim_rms / ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).pre_stim_rms;
            FiberPhotometryAnalysis.summary_stats(fiber_idx).post_pre_ratio = ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).post_stim_rms / ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).pre_stim_rms;
        end
        FiberPhotometryAnalysis.summary_stats(fiber_idx).mean_signal = mean(Vx);
        FiberPhotometryAnalysis.summary_stats(fiber_idx).std_signal = std(Vx);
    end
end

%% Save data
fprintf('Saving data structure...\n');
save_filename = sprintf('%s-%s-%s_FiberPhotometry_Analysis.mat', ...
    METADATA.mouse_name, METADATA.recording_date, METADATA.recording_id);
save_filepath = fullfile(save_directory, save_filename);
save(save_filepath, 'FiberPhotometryAnalysis', '-v7.3');

fprintf('Analysis data saved to: %s\n', save_filepath);

%% Save summary report
fprintf('Generating summary report...\n');
summary_filename = sprintf('%s-%s-%s_Summary.txt', ...
    METADATA.mouse_name, METADATA.recording_date, METADATA.recording_id);
summary_filepath = fullfile(save_directory, summary_filename);

fid = fopen(summary_filepath, 'w');
fprintf(fid, '=== COMPREHENSIVE FIBER PHOTOMETRY & LFP ANALYSIS SUMMARY (ENHANCED) ===\n');
fprintf(fid, 'Analysis Date: %s\n', datestr(now));
fprintf(fid, 'Mouse: %s\n', METADATA.mouse_name);
fprintf(fid, 'Recording Date: %s\n', METADATA.recording_date);
fprintf(fid, 'Recording ID: %s\n', METADATA.recording_id);
fprintf(fid, 'Analysis Mode: %s\n', ANALYSIS_MODE);
if PROCESS_SINGLE_TRIAL
    fprintf(fid, 'Trial: %s\n', trial_name);
end

fprintf(fid, '\n=== RECORDING PARAMETERS ===\n');
fprintf(fid, 'Sampling Rate: %.2f Hz\n', IMAGING_FS);
fprintf(fid, 'Total Frames: %d\n', size(processed_traces, 1));
fprintf(fid, 'Recording Duration: %.2f seconds\n', size(processed_traces, 1)/IMAGING_FS);
fprintf(fid, 'Number of Fibers: %d\n', size(processed_traces, 2));
fprintf(fid, 'Processing Mode: %s\n', char(string(PROCESS_FULL_FIELD)));
fprintf(fid, 'Motion Correction: %s\n', char(string(logical(MOTION_CORRECTION))));
fprintf(fid, 'Multi-Fiber Visualization: %s\n', char(string(GENERATE_MULTI_FIBER_PLOT)));
fprintf(fid, 'Baseline Trial: %s\n', char(string(IS_BASELINE_TRIAL)));

if PROCESS_SINGLE_TRIAL
    fprintf(fid, '\n=== TIME PERIODS ===\n');
    if IS_BASELINE_TRIAL || isempty(STIM_PERIOD)
        fprintf(fid, 'Trial Type: BASELINE (No stimulation)\n');
        total_time = size(processed_traces, 1) / IMAGING_FS;
        fprintf(fid, 'Early period: 0.0 - %.1f seconds\n', total_time/3);
        fprintf(fid, 'Middle period: %.1f - %.1f seconds\n', total_time/3, 2*total_time/3);
        fprintf(fid, 'Late period: %.1f - %.1f seconds\n', 2*total_time/3, total_time);
    else
        fprintf(fid, 'Trial Type: STIMULATION\n');
        fprintf(fid, 'Pre-stimulus: %.1f - %.1f seconds\n', PRE_STIM_PERIOD(1), PRE_STIM_PERIOD(2));
        fprintf(fid, 'Stimulation: %.1f - %.1f seconds\n', STIM_PERIOD(1), STIM_PERIOD(2));
        fprintf(fid, 'Post-stimulus: %.1f seconds to end\n', STIM_PERIOD(2));
    end

    % Add detection method information
    if exist('STIM_DETECTION', 'var')
        fprintf(fid, '\n=== STIMULATION DETECTION ===\n');
        fprintf(fid, 'Detection Method: %s\n', STIM_DETECTION.method);
        if strcmp(STIM_DETECTION.method, 'baseline_trial')
            fprintf(fid, 'Baseline trial detected: No stimulation pulses found\n');
            if exist('num_pulses_ch5', 'var') && exist('num_pulses_ch6', 'var')
                fprintf(fid, 'ADC5 pulses: %d\n', num_pulses_ch5);
                fprintf(fid, 'ADC6 pulses: %d\n', num_pulses_ch6);
            end
        elseif strcmp(STIM_DETECTION.method, 'automatic')
            fprintf(fid, 'Onset Sample: %d\n', STIM_DETECTION.onset_sample);
            fprintf(fid, 'Offset Sample: %d\n', STIM_DETECTION.offset_sample);
            if exist('STIM_PERIOD_HARDCODED', 'var')
                fprintf(fid, 'Difference from hardcoded onset: %.2f seconds\n', STIM_DETECTION.onset_diff_from_hardcoded);
                fprintf(fid, 'Difference from hardcoded offset: %.2f seconds\n', STIM_DETECTION.offset_diff_from_hardcoded);
            end
        else
            fprintf(fid, 'Reason: %s\n', STIM_DETECTION.reason);
        end
    end
end

fprintf(fid, '\n=== PHOTOBLEACHING CORRECTION ===\n');
fprintf(fid, 'Correction Applied: %s\n', char(string(APPLY_PHOTOBLEACHING_CORRECTION)));
if APPLY_PHOTOBLEACHING_CORRECTION
    fprintf(fid, 'Methods: Linear detrending, Double exponential\n');
    if exist('fitted_params', 'var')
        fprintf(fid, 'Exponential Fit Success: Yes\n');
    else
        fprintf(fid, 'Exponential Fit Success: No\n');
    end
end

fprintf(fid, '\n=== OPEN EPHYS DATA ===\n');
fprintf(fid, 'Ephys Data Loaded: %s\n', char(string(EPHYS_LOADED)));
if EPHYS_LOADED
    fprintf(fid, 'Ephys Sampling Rate: %d Hz\n', EPHYS_FS);
    fprintf(fid, 'Aligned Samples: %d\n', length(lfp_aligned));
    if exist('theta_envelope_correlation', 'var')
        fprintf(fid, 'Theta Envelope Correlation: %.4f\n', theta_envelope_correlation);
    end
end

if PROCESS_SINGLE_TRIAL && isfield(FiberPhotometryAnalysis, 'summary_stats')
    for fiber_idx = 1:length(FiberPhotometryAnalysis.summary_stats)
        fprintf(fid, '\n=== FIBER %d STATISTICS ===\n', fiber_idx);
        if isfield(FiberPhotometryAnalysis.summary_stats(fiber_idx), 'is_baseline_trial') && ...
                FiberPhotometryAnalysis.summary_stats(fiber_idx).is_baseline_trial
            fprintf(fid, 'Trial Type: BASELINE\n');
            fprintf(fid, 'Early period RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).early_rms);
            fprintf(fid, 'Middle period RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).middle_rms);
            fprintf(fid, 'Late period RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).late_rms);
            fprintf(fid, 'Middle/Early Ratio: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).middle_early_ratio);
            fprintf(fid, 'Late/Early Ratio: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).late_early_ratio);
        else
            fprintf(fid, 'Trial Type: STIMULATION\n');
            fprintf(fid, 'Pre-stim RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).pre_stim_rms);
            fprintf(fid, 'Stimulation RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).stim_rms);
            fprintf(fid, 'Post-stim RMS: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).post_stim_rms);
            fprintf(fid, 'Stim/Pre Ratio: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).stim_pre_ratio);
            fprintf(fid, 'Post/Pre Ratio: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).post_pre_ratio);
        end
        fprintf(fid, 'Mean Signal: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).mean_signal);
        fprintf(fid, 'Std Signal: %.4f\n', FiberPhotometryAnalysis.summary_stats(fiber_idx).std_signal);
    end
end

fprintf(fid, '\n=== ENHANCED FEATURES ===\n');
if GENERATE_MULTI_FIBER_PLOT && num_fibers > 1
    fprintf(fid, 'Multi-fiber combined visualization: Generated\n');
    fprintf(fid, 'Inter-fiber correlation analysis: Included\n');
    fprintf(fid, 'Multi-fiber power spectral comparison: Included\n');
else
    fprintf(fid, 'Multi-fiber visualization: Skipped (single fiber or disabled)\n');
end

fprintf(fid, '\n=== SAVED DATA STRUCTURE ===\n');
fprintf(fid, 'All plots can be regenerated from saved .mat file\n');
fprintf(fid, 'Enhanced multi-fiber capabilities integrated\n');

fclose(fid);
fprintf('Summary report saved to: %s\n', summary_filepath);

%% ============================================================================
%  COMPLETION MESSAGE
%  ============================================================================

fprintf('\n=== ENHANCED FIBER PHOTOMETRY ANALYSIS COMPLETE ===\n');
fprintf('Results saved to: %s\n', save_directory);
fprintf('Number of fibers processed: %d\n', num_fibers);
if EPHYS_LOADED
    fprintf('LFP integration: Successful\n');
end
if GENERATE_MULTI_FIBER_PLOT && num_fibers > 1
    fprintf('Multi-fiber combined visualization: Generated\n');
end
fprintf('All original functionality preserved\n');
fprintf('Enhanced features successfully integrated\n');

%% ============================================================================
%  HELPER FUNCTIONS
%  ============================================================================
%  Moved to core/utils/ so the single- and multi-trial scripts share ONE copy:
%      smooth2a, fastsmooth, viridis, turbo (+ generate_biphasic_pulses)
%  They are added to the MATLAB path (on top, so these versions are used rather
%  than any external or built-in ones) by the path-setup block near the top of
%  this script. See core/utils/ and core/utils/tests/test_core_utils.m.