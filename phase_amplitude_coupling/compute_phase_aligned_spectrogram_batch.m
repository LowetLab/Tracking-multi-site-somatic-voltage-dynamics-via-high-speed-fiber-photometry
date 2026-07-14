%% compute_phase_aligned_spectrogram_batch.m
%% ============================================================================
%  BATCH Phase-Aligned Wavelet Spectrogram Analysis
%  ============================================================================
%
%  DESCRIPTION:
%  This script batch processes all *_FiberPhotometry_Analysis.mat files in
%  BaselineData folder and computes phase-aligned wavelet spectrograms for
%  LFP and Fiber signals.
%
%  Output is saved to Process/phase_aligned_spectrogram/ with matching folder
%  structure. No GUI display - figures are saved directly to disk.
%  ============================================================================

clear; clc; close all;

%% ==================== Configuration Parameters ====================
% Low frequency band for phase extraction (carrier frequency)
LowFreq_Band = [5, 9];           % Hz - band for phase extraction
Phase_Method = 'wavelet';        % 'wavelet' or 'hilbert' - method for phase extraction

% High frequency range for wavelet analysis
% Extended down to 2 Hz so theta-band (5-9 Hz) amplitude is available for
% phase-curve plots. 95 bins gives exactly 1 Hz resolution over [6, 100].
HighFreq_Range = [2, 100];       % Hz - range for wavelet spectrogram
HighFreq_NumFreqs = 95;          % Number of frequency bins (1 Hz resolution: (100-6)/(95-1)=1 Hz)

% Epoch and phase bin settings
Cycles_Per_Epoch = 2;            % Number of low-frequency cycles per epoch
Num_Phase_Bins = 36;             % Number of phase bins (-π to π), each bin = 10°

% Phase source configuration
Phase_Source = 'lfp';            % Options: 'self' or 'lfp' (using LFP phase for Fiber)

% Wavelet parameters
Wavelet_Cycles = 5;              % Number of cycles for Morlet wavelet
% NOTE: Wavelet output is MAGNITUDE (not power), colorbar shows signal amplitude
% NOTE: Fiber signals are multiplied by 100 (percentage representation)

% Batch processing settings
SKIP_EXISTING = false;            % Set to false to reprocess all files
SAVE_FIGURES = true;             % Save figures to disk
FIGURE_FORMAT = {'png', 'fig'}; % Output formats

% ---------- Animal / path filter (optional) ----------
% Set to '' to process ALL files; otherwise only files whose full path
% contains this substring will be processed.
% Examples:
%   ANIMAL_FILTER = 'Animal01';                   % one animal
%   ANIMAL_FILTER = 'PV_Animals';                  % one condition group
%   ANIMAL_FILTER = 'Animal01\27_01_26-R1';        % one specific recording
%   ANIMAL_FILTER = '';                             % process everything
ANIMAL_FILTER = '';

% Artifact mask settings
% Artifact masks are loaded from Artifact_masks/ in the project root.
% Each epoch will receive has_artifact (logical) and artifact_fraction (0-1).
% If no mask file is found for a trial, has_artifact = false for all epochs.
ARTIFACT_MASKS_DIR = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Artifact_masks');

%% ==================== Add Required Paths ====================
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);

% Add tortlab path for eegfilt
tortlab_path = fullfile(script_dir, 'sample', 'tortlab_phase_amplitude_coupling_master', 'phase-amplitude-coupling-master');
if exist(tortlab_path, 'dir')
    addpath(tortlab_path);
    fprintf('Added tortlab path.\n');
else
    error('Tortlab path not found: %s', tortlab_path);
end

%% ==================== Find All Analysis Files ====================
baseline_data_path = fullfile(project_root, 'BaselineData');
if ~exist(baseline_data_path, 'dir')
    error('BaselineData folder not found: %s', baseline_data_path);
end

fprintf('\n=== Scanning for analysis files in BaselineData ===\n');

% Recursively find all *_FiberPhotometry_Analysis.mat files
all_files = dir(fullfile(baseline_data_path, '**', '*_FiberPhotometry_Analysis.mat'));

if isempty(all_files)
    error('No *_FiberPhotometry_Analysis.mat files found in BaselineData folder.');
end

fprintf('Found %d files total.\n', length(all_files));

% Apply ANIMAL_FILTER if set
if ~isempty(ANIMAL_FILTER)
    keep = false(length(all_files), 1);
    for fi = 1:length(all_files)
        keep(fi) = contains(fullfile(all_files(fi).folder, all_files(fi).name), ANIMAL_FILTER);
    end
    all_files = all_files(keep);
    if isempty(all_files)
        error('ANIMAL_FILTER = ''%s'' matched 0 files. Check spelling / path.', ANIMAL_FILTER);
    end
    fprintf('ANIMAL_FILTER = ''%s'' → %d file(s) selected.\n', ANIMAL_FILTER, length(all_files));
end

%% ==================== Prepare Parameters ====================
% Define frequency vector for wavelet analysis
freq_vector = linspace(HighFreq_Range(1), HighFreq_Range(2), HighFreq_NumFreqs);

% Define phase bins
phase_bin_edges = linspace(-pi, pi, Num_Phase_Bins + 1);
phase_bin_centers = (phase_bin_edges(1:end-1) + phase_bin_edges(2:end)) / 2;

%% ==================== Batch Processing Loop ====================
fprintf('\n=== Starting Batch Processing ===\n');
total_files = length(all_files);
processed_count = 0;
skipped_count = 0;
error_count = 0;
error_files = {};

tic;  % Start timer

for file_idx = 1:total_files
    current_file = all_files(file_idx);
    filepath = current_file.folder;
    filename = current_file.name;
    full_mat_path = fullfile(filepath, filename);
    [~, base_filename_full, ~] = fileparts(filename);
    base_filename = strrep(base_filename_full, '_FiberPhotometry_Analysis', '');
    
    fprintf('\n[%d/%d] Processing: %s\n', file_idx, total_files, base_filename);
    
    %% Construct output path
    baseline_keyword = 'BaselineData';
    idx = strfind(filepath, baseline_keyword);
    if ~isempty(idx)
        relative_path = filepath(idx(1) + length(baseline_keyword) + 1 : end);
    else
        path_parts = strsplit(filepath, filesep);
        relative_path = fullfile(path_parts{max(1,end-3):end});
    end
    
    output_dir = fullfile(project_root, 'Process', 'phase_aligned_spectrogram', relative_path);
    
    % Check if already processed
    results_filename = sprintf('PhaseAlignedSpectrogram_Results_%s.mat', base_filename);
    results_path = fullfile(output_dir, results_filename);
    
    if SKIP_EXISTING && exist(results_path, 'file')
        fprintf('  -> Already processed, skipping.\n');
        skipped_count = skipped_count + 1;
        continue;
    end
    
    try
        %% Load data
        data = load(full_mat_path);
        
        if ~isfield(data, 'FiberPhotometryAnalysis')
            warning('  FiberPhotometryAnalysis not found, skipping.');
            error_count = error_count + 1;
            error_files{end+1} = filename;
            continue;
        end
        
        FPA = data.FiberPhotometryAnalysis;
        
        % Get sampling rate
        if isfield(FPA, 'time') && isfield(FPA.time, 'sampling_rate')
            fs = FPA.time.sampling_rate;
        elseif isfield(FPA, 'parameters') && isfield(FPA.parameters, 'sampling_rate')
            fs = FPA.parameters.sampling_rate;
        else
            warning('  Sampling rate not found, skipping.');
            error_count = error_count + 1;
            error_files{end+1} = filename;
            continue;
        end
        
        % Get time vector
        if isfield(FPA, 'time') && isfield(FPA.time, 'time_vector_seconds')
            time_vector = FPA.time.time_vector_seconds;
        else
            warning('  Time vector not found, skipping.');
            error_count = error_count + 1;
            error_files{end+1} = filename;
            continue;
        end
        
        % LFP: 顺位读取，列表靠前优先；前若干项为未 z-score 的电压迹，其后为 z-score 回退
        LFP_EPHYS_CANDIDATES = { ...
            'lfp_raw_aligned_HP', 'lfp_raw_aligned_mPFC', 'lfp_aligned', ...
            'lfp_z_HP', 'lfp_z' };
        lfp_data = [];
        lfp_field_name = '';
        if isfield(FPA, 'ephys')
            for ci = 1:numel(LFP_EPHYS_CANDIDATES)
                fn = LFP_EPHYS_CANDIDATES{ci};
                if isfield(FPA.ephys, fn)
                    lfp_data = FPA.ephys.(fn);
                    lfp_field_name = fn;
                    break;
                end
            end
        end
        
        if isempty(lfp_data)
            warning('  LFP data not found, skipping.');
            error_count = error_count + 1;
            error_files{end+1} = filename;
            continue;
        end
        lfp_data = lfp_data(:);
        fprintf('  LFP source: FPA.ephys.%s\n', lfp_field_name);
        lfp_mag_cbar_yl = lfp_wavelet_colorbar_ylabel(lfp_field_name);
        fiber_mag_cbar_yl = 'Mag. (%)';
        
        % Get Fiber data
        if ~isfield(FPA, 'signals') || ~isfield(FPA.signals, 'final_processed_traces')
            warning('  Fiber traces not found, skipping.');
            error_count = error_count + 1;
            error_files{end+1} = filename;
            continue;
        end
        fiber_data = FPA.signals.final_processed_traces;
        num_fibers = size(fiber_data, 2);
        
        % Get running velocity
        running_velocity = [];
        if isfield(FPA, 'ephys') && isfield(FPA.ephys, 'running_velocity_smooth')
            running_velocity = FPA.ephys.running_velocity_smooth;
            running_velocity = running_velocity(:);
        end
        
        % Build fiber-to-ROI mapping
        fiber_roi_map = build_fiber_roi_map(FPA);

        % Load artifact mask for this trial (if available)
        artifact_mask = load_artifact_mask(ARTIFACT_MASKS_DIR, base_filename, filename);
        if ~isempty(artifact_mask)
            fprintf('  Artifact mask loaded: %.1f%% of samples flagged.\n', ...
                    100 * mean(double(artifact_mask)));
        else
            fprintf('  No artifact mask found for this trial.\n');
        end

        % Create output directory
        if ~exist(output_dir, 'dir')
            mkdir(output_dir);
        end
        
        fprintf('  Sampling rate: %.1f Hz | Fibers: %d\n', fs, num_fibers);
        
        %% Initialize results structure
        results = struct();
        results.metadata.source_file = full_mat_path;
        results.metadata.base_filename = base_filename;
        results.metadata.sampling_rate = fs;
        results.metadata.time_vector = time_vector;
        results.metadata.low_freq_band = LowFreq_Band;
        results.metadata.phase_method = Phase_Method;
        results.metadata.high_freq_range = HighFreq_Range;
        results.metadata.freq_vector = freq_vector;
        results.metadata.phase_bin_centers = phase_bin_centers;
        results.metadata.phase_bin_edges = phase_bin_edges;
        results.metadata.num_phase_bins = Num_Phase_Bins;
        results.metadata.cycles_per_epoch = Cycles_Per_Epoch;
        results.metadata.phase_source = Phase_Source;
        results.metadata.wavelet_cycles = Wavelet_Cycles;
        results.metadata.artifact_mask_available = ~isempty(artifact_mask);
        results.metadata.lfp_ephys_field = lfp_field_name;

        %% Process LFP
        fprintf('  Processing LFP...\n');
        
        % Extract LFP phase
        lfp_phase = extract_instantaneous_phase(lfp_data, fs, LowFreq_Band, Phase_Method, Wavelet_Cycles);
        
        % Compute wavelet spectrogram
        fprintf('    Computing wavelet spectrogram...');
        [lfp_spectrogram, ~] = compute_wavelet_spectrogram(lfp_data, fs, freq_vector, Wavelet_Cycles);
        fprintf(' Done.\n');
        
        % Detect epochs
        [epoch_starts, epoch_ends] = detect_epochs(lfp_phase, Cycles_Per_Epoch);
        fprintf('    Found %d epochs\n', length(epoch_starts));
        
        % Pre-filter full LFP signal for phase-triggered average.
        % Filtering on the FULL signal (then slicing by epoch) avoids the
        % systematic filtfilt boundary artefact that occurs when each epoch
        % is filtered independently (epochs always start at the trough).
        [b_pta, a_pta] = butter(4, LowFreq_Band / (fs/2), 'bandpass');
        lfp_theta_filt  = filtfilt(b_pta, a_pta, double(lfp_data));

        % Align epochs by phase
        [mean_spectrogram, all_epochs_raw, epoch_spectrograms, cycle1_mean_spectrogram, cycle2_mean_spectrogram] = ...
            align_epochs_by_phase(lfp_spectrogram, lfp_phase, lfp_data, lfp_theta_filt, epoch_starts, epoch_ends, ...
                                  phase_bin_edges, freq_vector, running_velocity);
        
        % Annotate epochs with artifact information
        all_epochs_raw = annotate_epochs_artifact(all_epochs_raw, artifact_mask);

        % Compute epoch locomotion summary
        fprintf('    Summarizing locomotion...');
        epoch_mean_velocity = zeros(length(epoch_starts), 1);
        for ei = 1:length(epoch_starts)
            epoch_mean_velocity(ei) = all_epochs_raw{ei}.mean_velocity;
        end
        fprintf(' Done.\n');
        
        % Store LFP results
        fprintf('    Storing LFP results...');
        results.LFP.phase_bins = phase_bin_centers;
        results.LFP.freq_vector = freq_vector;
        results.LFP.mean_spectrogram = mean_spectrogram;
        results.LFP.cycle1_mean_spectrogram = cycle1_mean_spectrogram;
        results.LFP.cycle2_mean_spectrogram = cycle2_mean_spectrogram;
        results.LFP.all_epochs_raw = all_epochs_raw;
        results.LFP.epoch_spectrograms = epoch_spectrograms;
        results.LFP.n_epochs = length(epoch_starts);
        results.LFP.instantaneous_phase = lfp_phase;
        results.LFP.full_spectrogram = lfp_spectrogram;
        results.LFP.epoch_mean_velocity = epoch_mean_velocity;
        fprintf(' Done.\n');
        
        %% Generate and save LFP figure
        if SAVE_FIGURES
            fprintf('    Generating LFP figure...');
            fig_lfp = figure('Position', [100, 100, 1400, 1000], 'Visible', 'off');
            
            % Row 1: Wavelet spectrogram
            ax1 = subplot(4, 2, [1, 2]);
            imagesc(time_vector, freq_vector, results.LFP.full_spectrogram);
            axis xy;
            cb_lfp = colorbar;
            ylabel(cb_lfp, lfp_mag_cbar_yl, 'FontSize', 8);
            ylabel('Frequency (Hz)');
            title(sprintf('LFP Wavelet Spectrogram - %s', strrep(base_filename, '_', '\_')));
            colormap(gca, 'jet');
            set(gca, 'XTickLabel', []);
            xlim([time_vector(1), time_vector(end)]);
            
            % Row 2: Running velocity
            ax2 = subplot(4, 2, [3, 4]);
            if ~isempty(running_velocity)
                plot(time_vector, running_velocity(1:length(time_vector)), 'k-', 'LineWidth', 0.8);
            end
            ylabel('Velocity');
            title('Running Velocity (smoothed)');
            xlim([time_vector(1), time_vector(end)]);
            set(gca, 'XTickLabel', []);
            grid on; box off;
            
            % Row 3: Instantaneous phase
            ax3 = subplot(4, 2, [5, 6]);
            plot(time_vector, results.LFP.instantaneous_phase, 'b-', 'LineWidth', 0.5);
            xlabel('Time (s)');
            ylabel('Phase (rad)');
            title(sprintf('LFP Instantaneous Phase (%d-%d Hz)', LowFreq_Band(1), LowFreq_Band(2)));
            ylim([-pi, pi]);
            yticks([-pi, -pi/2, 0, pi/2, pi]);
            yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
            xlim([time_vector(1), time_vector(end)]);
            grid on;
            
            % linkaxes([ax1, ax2, ax3], 'x');  % Disabled for batch mode
            
            % Row 4 Left: 1 cycle
            subplot(4, 2, 7);
            imagesc(rad2deg(phase_bin_centers), freq_vector, results.LFP.mean_spectrogram);
            axis xy;
            cb_lfp2 = colorbar;
            ylabel(cb_lfp2, lfp_mag_cbar_yl, 'FontSize', 8);
            xlabel('Phase (degrees)');
            ylabel('Frequency (Hz)');
            title('Phase-Aligned (1 cycle)');
            colormap(gca, 'jet');
            xticks([-180, -90, 0, 90, 180]);
            
            % Row 4 Right: 2 cycles
            subplot(4, 2, 8);
            two_cycle_spectrogram = [results.LFP.cycle1_mean_spectrogram, results.LFP.cycle2_mean_spectrogram];
            two_cycle_phase = [rad2deg(phase_bin_centers), rad2deg(phase_bin_centers) + 360];
            imagesc(two_cycle_phase, freq_vector, two_cycle_spectrogram);
            axis xy;
            cb_lfp3 = colorbar;
            ylabel(cb_lfp3, lfp_mag_cbar_yl, 'FontSize', 8);
            xlabel('Phase (degrees)');
            ylabel('Frequency (Hz)');
            title('Phase-Aligned (2 cycles)');
            colormap(gca, 'jet');
            xticks([-180, 0, 180, 360, 540]);
            hold on;
            plot([180, 180], [freq_vector(1), freq_vector(end)], 'w--', 'LineWidth', 1.5);
            hold off;
            
            fprintf(' Saving...');
            % Save figure
            lfp_fig_name = sprintf('PhaseAlignedSpectrogram_LFP_%s', base_filename);
            for fmt_idx = 1:length(FIGURE_FORMAT)
                fmt = FIGURE_FORMAT{fmt_idx};
                switch fmt
                    case 'fig'
                        saveas(fig_lfp, fullfile(output_dir, [lfp_fig_name, '.fig']));
                    case 'png'
                        print(fig_lfp, fullfile(output_dir, [lfp_fig_name, '.png']), '-dpng', '-r150');
                end
            end
            close(fig_lfp);
            fprintf(' Done.\n');
        end
        
        %% Process each Fiber
        for fiber_idx = 1:num_fibers
            fprintf('  Processing Fiber %d/%d...\n', fiber_idx, num_fibers);
            
            fiber_signal = fiber_data(:, fiber_idx);
            
            % Convert to percentage (multiply by 100)
            % This makes the wavelet magnitude represent % signal change
            fiber_signal = fiber_signal * 100;
            
            % PV voltage-imaging indicator: signal polarity is inverted.
            % Flip along Y-axis (y = -y) so fiber matches conventional interpretation.
            if contains(relative_path, 'PV_Animals')
                fiber_signal = -fiber_signal;
                fprintf('    [PV animal] Fiber signal flipped (Y-axis).\n');
            end
            
            % Determine phase source
            if strcmp(Phase_Source, 'lfp')
                current_phase = lfp_phase;
            else
                current_phase = extract_instantaneous_phase(fiber_signal, fs, LowFreq_Band, Phase_Method, Wavelet_Cycles);
            end
            
            % Compute wavelet spectrogram
            fprintf('    Computing wavelet spectrogram...');
            [fiber_spectrogram, ~] = compute_wavelet_spectrogram(fiber_signal, fs, freq_vector, Wavelet_Cycles);
            fprintf(' Done.\n');
            
            % Detect epochs
            [epoch_starts, epoch_ends] = detect_epochs(current_phase, Cycles_Per_Epoch);
            
            % Pre-filter full Fiber signal for phase-triggered average
            [b_pta_f, a_pta_f] = butter(4, LowFreq_Band / (fs/2), 'bandpass');
            fiber_theta_filt = filtfilt(b_pta_f, a_pta_f, double(fiber_signal));

            % Align epochs by phase
            [mean_spectrogram, all_epochs_raw, epoch_spectrograms, cycle1_mean_spectrogram, cycle2_mean_spectrogram] = ...
                align_epochs_by_phase(fiber_spectrogram, current_phase, fiber_signal, fiber_theta_filt, epoch_starts, epoch_ends, ...
                                      phase_bin_edges, freq_vector, running_velocity);
            
            % Annotate epochs with artifact information
            all_epochs_raw = annotate_epochs_artifact(all_epochs_raw, artifact_mask);

            % Compute epoch locomotion summary
            epoch_mean_velocity = zeros(length(epoch_starts), 1);
            for ei = 1:length(epoch_starts)
                epoch_mean_velocity(ei) = all_epochs_raw{ei}.mean_velocity;
            end
            
            % Get ROI info
            roi_info = fiber_roi_map{fiber_idx};
            
            % Store Fiber results
            field_name = sprintf('Fiber%d', fiber_idx);
            results.(field_name).phase_bins = phase_bin_centers;
            results.(field_name).freq_vector = freq_vector;
            results.(field_name).mean_spectrogram = mean_spectrogram;
            results.(field_name).cycle1_mean_spectrogram = cycle1_mean_spectrogram;
            results.(field_name).cycle2_mean_spectrogram = cycle2_mean_spectrogram;
            results.(field_name).all_epochs_raw = all_epochs_raw;
            results.(field_name).epoch_spectrograms = epoch_spectrograms;
            results.(field_name).n_epochs = length(epoch_starts);
            results.(field_name).instantaneous_phase = current_phase;
            results.(field_name).full_spectrogram = fiber_spectrogram;
            results.(field_name).roi_info = roi_info;
            results.(field_name).epoch_mean_velocity = epoch_mean_velocity;
            
            %% Generate and save Fiber figure
            if SAVE_FIGURES
                fig_fiber = figure('Position', [100, 100, 1400, 1000], 'Visible', 'off');
                
                % Row 1: Wavelet spectrogram
                ax1 = subplot(4, 2, [1, 2]);
                imagesc(time_vector, freq_vector, results.(field_name).full_spectrogram);
                axis xy;
                cb_f = colorbar;
                ylabel(cb_f, fiber_mag_cbar_yl, 'FontSize', 8);
                ylabel('Frequency (Hz)');
                title(sprintf('Fiber%d (%s) Wavelet Spectrogram - %s', fiber_idx, roi_info, strrep(base_filename, '_', '\_')));
                colormap(gca, 'jet');
                set(gca, 'XTickLabel', []);
                xlim([time_vector(1), time_vector(end)]);
                
                % Row 2: Running velocity
                ax2 = subplot(4, 2, [3, 4]);
                if ~isempty(running_velocity)
                    plot(time_vector, running_velocity(1:length(time_vector)), 'k-', 'LineWidth', 0.8);
                end
                ylabel('Velocity');
                title('Running Velocity (smoothed)');
                xlim([time_vector(1), time_vector(end)]);
                set(gca, 'XTickLabel', []);
                grid on; box off;
                
                % Row 3: Instantaneous phase
                ax3 = subplot(4, 2, [5, 6]);
                plot(time_vector, results.(field_name).instantaneous_phase, 'g-', 'LineWidth', 0.5);
                xlabel('Time (s)');
                ylabel('Phase (rad)');
                title(sprintf('Fiber%d Instantaneous Phase (%d-%d Hz)', fiber_idx, LowFreq_Band(1), LowFreq_Band(2)));
                ylim([-pi, pi]);
                yticks([-pi, -pi/2, 0, pi/2, pi]);
                yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
                xlim([time_vector(1), time_vector(end)]);
                grid on;
                
                % linkaxes([ax1, ax2, ax3], 'x');  % Disabled for batch mode
                
                % Row 4 Left: 1 cycle
                subplot(4, 2, 7);
                imagesc(rad2deg(phase_bin_centers), freq_vector, results.(field_name).mean_spectrogram);
                axis xy;
                cb_f2 = colorbar;
                ylabel(cb_f2, fiber_mag_cbar_yl, 'FontSize', 8);
                xlabel('Phase (degrees)');
                ylabel('Frequency (Hz)');
                title('Phase-Aligned (1 cycle)');
                colormap(gca, 'jet');
                xticks([-180, -90, 0, 90, 180]);
                
                % Row 4 Right: 2 cycles
                subplot(4, 2, 8);
                two_cycle_spectrogram = [results.(field_name).cycle1_mean_spectrogram, results.(field_name).cycle2_mean_spectrogram];
                two_cycle_phase = [rad2deg(phase_bin_centers), rad2deg(phase_bin_centers) + 360];
                imagesc(two_cycle_phase, freq_vector, two_cycle_spectrogram);
                axis xy;
                cb_f3 = colorbar;
                ylabel(cb_f3, fiber_mag_cbar_yl, 'FontSize', 8);
                xlabel('Phase (degrees)');
                ylabel('Frequency (Hz)');
                title('Phase-Aligned (2 cycles)');
                colormap(gca, 'jet');
                xticks([-180, 0, 180, 360, 540]);
                hold on;
                plot([180, 180], [freq_vector(1), freq_vector(end)], 'w--', 'LineWidth', 1.5);
                hold off;
                
                % Save figure
                fiber_fig_name = sprintf('PhaseAlignedSpectrogram_Fiber%d_%s_%s', fiber_idx, roi_info, base_filename);
                for fmt_idx = 1:length(FIGURE_FORMAT)
                    fmt = FIGURE_FORMAT{fmt_idx};
                    switch fmt
                        case 'fig'
                            saveas(fig_fiber, fullfile(output_dir, [fiber_fig_name, '.fig']));
                        case 'png'
                            print(fig_fiber, fullfile(output_dir, [fiber_fig_name, '.png']), '-dpng', '-r150');
                    end
                end
                close(fig_fiber);
            end
        end
        
        %% Save results
        save(fullfile(output_dir, results_filename), 'results', '-v7.3');
        
        % Save Python-friendly version
        python_results = struct();
        python_results.metadata = results.metadata;
        python_results.phase_bins_rad = phase_bin_centers;
        python_results.phase_bins_deg = rad2deg(phase_bin_centers);
        python_results.freq_vector = freq_vector;
        python_results.two_cycle_phase_deg = [rad2deg(phase_bin_centers), rad2deg(phase_bin_centers) + 360];
        
        if isfield(results, 'LFP')
            python_results.LFP_mean_spectrogram = results.LFP.mean_spectrogram;
            python_results.LFP_cycle1_mean_spectrogram = results.LFP.cycle1_mean_spectrogram;
            python_results.LFP_cycle2_mean_spectrogram = results.LFP.cycle2_mean_spectrogram;
            python_results.LFP_two_cycle_spectrogram = [results.LFP.cycle1_mean_spectrogram, results.LFP.cycle2_mean_spectrogram];
            python_results.LFP_n_epochs = results.LFP.n_epochs;
            python_results.LFP_epoch_spectrograms = results.LFP.epoch_spectrograms;
            python_results.LFP_epoch_mean_velocity = results.LFP.epoch_mean_velocity;
        end
        
        for fiber_idx = 1:num_fibers
            field_name = sprintf('Fiber%d', fiber_idx);
            python_results.(sprintf('Fiber%d_mean_spectrogram', fiber_idx)) = results.(field_name).mean_spectrogram;
            python_results.(sprintf('Fiber%d_cycle1_mean_spectrogram', fiber_idx)) = results.(field_name).cycle1_mean_spectrogram;
            python_results.(sprintf('Fiber%d_cycle2_mean_spectrogram', fiber_idx)) = results.(field_name).cycle2_mean_spectrogram;
            python_results.(sprintf('Fiber%d_two_cycle_spectrogram', fiber_idx)) = [results.(field_name).cycle1_mean_spectrogram, results.(field_name).cycle2_mean_spectrogram];
            python_results.(sprintf('Fiber%d_n_epochs', fiber_idx)) = results.(field_name).n_epochs;
            python_results.(sprintf('Fiber%d_epoch_spectrograms', fiber_idx)) = results.(field_name).epoch_spectrograms;
            python_results.(sprintf('Fiber%d_roi_info', fiber_idx)) = results.(field_name).roi_info;
            python_results.(sprintf('Fiber%d_epoch_mean_velocity', fiber_idx)) = results.(field_name).epoch_mean_velocity;
        end
        
        python_filename = sprintf('PhaseAlignedSpectrogram_ForPython_%s.mat', base_filename);
        save(fullfile(output_dir, python_filename), '-struct', 'python_results', '-v7');
        
        fprintf('  -> Saved to: %s\n', output_dir);
        processed_count = processed_count + 1;
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        error_count = error_count + 1;
        error_files{end+1} = filename;
    end
end

%% ==================== Summary ====================
elapsed_time = toc;

fprintf('\n========================================\n');
fprintf('       BATCH PROCESSING COMPLETE       \n');
fprintf('========================================\n');
fprintf('Total files:     %d\n', total_files);
fprintf('Processed:       %d\n', processed_count);
fprintf('Skipped:         %d\n', skipped_count);
fprintf('Errors:          %d\n', error_count);
fprintf('Elapsed time:    %.1f seconds (%.1f min)\n', elapsed_time, elapsed_time/60);

if ~isempty(error_files)
    fprintf('\nFiles with errors:\n');
    for i = 1:length(error_files)
        fprintf('  - %s\n', error_files{i});
    end
end

fprintf('========================================\n');

%% ==================== Helper Functions ====================

function phase = extract_instantaneous_phase(signal, fs, freq_band, method, n_cycles)
    if nargin < 4
        method = 'hilbert';
    end
    if nargin < 5
        n_cycles = 5;
    end
    
    signal = signal(:)';
    
    if strcmpi(method, 'wavelet')
        % Wavelet-based phase extraction
        % Use center frequency of the band
        center_freq = mean(freq_band);
        
        % Morlet wavelet at center frequency
        sigma_t = n_cycles / (2 * pi * center_freq);
        t_wavelet = -3*sigma_t : 1/fs : 3*sigma_t;
        wavelet = exp(2*1i*pi*center_freq*t_wavelet) .* exp(-t_wavelet.^2 / (2*sigma_t^2));
        wavelet = wavelet / sum(abs(wavelet));
        
        % Convolve and extract phase
        convolved = conv(signal, wavelet, 'same');
        phase = angle(convolved);
    else
        % Traditional Hilbert-based phase extraction
        filtered_signal = eegfilt(signal, fs, freq_band(1), freq_band(2));
        analytic_signal = hilbert(filtered_signal);
        phase = angle(analytic_signal);
    end
    
    phase = phase(:);
end

function [spectrogram_matrix, time_out] = compute_wavelet_spectrogram(signal, fs, freq_vector, n_cycles)
    % Compute wavelet spectrogram using Morlet wavelet
    % Output: magnitude (not power) - i.e., abs(convolved) instead of abs(convolved)^2
    signal = signal(:);
    n_samples = length(signal);
    n_freqs = length(freq_vector);
    
    time_out = (0:n_samples-1) / fs;
    spectrogram_matrix = zeros(n_freqs, n_samples);
    
    for fi = 1:n_freqs
        f = freq_vector(fi);
        sigma_t = n_cycles / (2 * pi * f);
        t_wavelet = -3*sigma_t : 1/fs : 3*sigma_t;
        wavelet = exp(2*1i*pi*f*t_wavelet) .* exp(-t_wavelet.^2 / (2*sigma_t^2));
        wavelet = wavelet / sum(abs(wavelet));
        convolved = conv(signal, wavelet, 'same');
        spectrogram_matrix(fi, :) = abs(convolved);  % Magnitude (not power)
    end
end

function [epoch_starts, epoch_ends] = detect_epochs(phase, cycles_per_epoch)
    phase_diff = diff(phase);
    wrap_points = find(phase_diff < -pi);
    cycle_starts = wrap_points + 1;
    cycle_starts = cycle_starts(cycle_starts <= length(phase));
    
    n_cycles = length(cycle_starts);
    n_epochs = floor(n_cycles / cycles_per_epoch);
    
    if n_epochs == 0
        epoch_starts = [];
        epoch_ends = [];
        return;
    end
    
    epoch_starts = zeros(n_epochs, 1);
    epoch_ends = zeros(n_epochs, 1);
    
    for ei = 1:n_epochs
        start_cycle_idx = (ei - 1) * cycles_per_epoch + 1;
        end_cycle_idx = ei * cycles_per_epoch;
        
        epoch_starts(ei) = cycle_starts(start_cycle_idx);
        
        if end_cycle_idx < n_cycles
            epoch_ends(ei) = cycle_starts(end_cycle_idx + 1) - 1;
        else
            if ei > 1
                avg_epoch_length = mean(epoch_ends(1:ei-1) - epoch_starts(1:ei-1));
                epoch_ends(ei) = min(epoch_starts(ei) + round(avg_epoch_length), length(phase));
            else
                epoch_ends(ei) = length(phase);
            end
        end
    end
    
    min_epoch_length = 10;
    valid = (epoch_ends - epoch_starts) >= min_epoch_length;
    epoch_starts = epoch_starts(valid);
    epoch_ends = epoch_ends(valid);
end

function [mean_spectrogram, all_epochs_raw, epoch_spectrograms, cycle1_mean_spectrogram, cycle2_mean_spectrogram] = ...
    align_epochs_by_phase(full_spectrogram, phase, raw_signal, theta_filtered, epoch_starts, epoch_ends, phase_bin_edges, freq_vector, running_velocity)
    
    n_freqs = length(freq_vector);
    n_bins = length(phase_bin_edges) - 1;
    n_epochs = length(epoch_starts);
    
    has_velocity = ~isempty(running_velocity);
    
    all_epochs_raw = cell(n_epochs, 1);
    epoch_spectrograms = zeros(n_freqs, n_bins, n_epochs);
    
    bin_power_sum = zeros(n_freqs, n_bins);
    bin_count = zeros(1, n_bins);
    
    cycle1_power_sum = zeros(n_freqs, n_bins);
    cycle1_count = zeros(1, n_bins);
    cycle2_power_sum = zeros(n_freqs, n_bins);
    cycle2_count = zeros(1, n_bins);
    
    % Progress display
    fprintf('    Aligning %d epochs by phase...', n_epochs);
    
    for ei = 1:n_epochs
        start_idx = epoch_starts(ei);
        end_idx = epoch_ends(ei);
        
        epoch_spec = full_spectrogram(:, start_idx:end_idx);
        epoch_phase = phase(start_idx:end_idx);
        epoch_len = length(epoch_phase);
        
        if has_velocity && end_idx <= length(running_velocity)
            epoch_velocity = running_velocity(start_idx:end_idx);
            epoch_mean_velocity = mean(epoch_velocity, 'omitnan');
        else
            epoch_velocity = [];
            epoch_mean_velocity = NaN;
        end
        
        phase_diff = diff(epoch_phase);
        wrap_points = find(phase_diff < -pi);
        
        if ~isempty(wrap_points)
            cycle_boundary = wrap_points(1);
        else
            cycle_boundary = round(epoch_len / 2);
        end
        
        cycle1_indices = 1:cycle_boundary;
        cycle2_indices = (cycle_boundary+1):epoch_len;
        
        all_epochs_raw{ei}.spectrogram  = epoch_spec;
        all_epochs_raw{ei}.phase        = epoch_phase;
        all_epochs_raw{ei}.raw_signal       = raw_signal(start_idx:end_idx);
        all_epochs_raw{ei}.theta_filtered   = theta_filtered(start_idx:end_idx);  % pre-filtered on full signal
        all_epochs_raw{ei}.time_indices     = [start_idx, end_idx];
        all_epochs_raw{ei}.cycle_boundary = cycle_boundary;
        all_epochs_raw{ei}.velocity = epoch_velocity;
        all_epochs_raw{ei}.mean_velocity = epoch_mean_velocity;
        
        epoch_binned = zeros(n_freqs, n_bins);
        epoch_bin_count = zeros(1, n_bins);
        
        for bi = 1:n_bins
            in_bin = epoch_phase >= phase_bin_edges(bi) & epoch_phase < phase_bin_edges(bi+1);
            if any(in_bin)
                epoch_binned(:, bi) = mean(epoch_spec(:, in_bin), 2);
                epoch_bin_count(bi) = sum(in_bin);
            end
        end
        
        epoch_spectrograms(:, :, ei) = epoch_binned;
        
        for bi = 1:n_bins
            if epoch_bin_count(bi) > 0
                bin_power_sum(:, bi) = bin_power_sum(:, bi) + epoch_binned(:, bi);
                bin_count(bi) = bin_count(bi) + 1;
            end
        end
        
        if ~isempty(cycle1_indices)
            cycle1_phase = epoch_phase(cycle1_indices);
            cycle1_spec = epoch_spec(:, cycle1_indices);
            for bi = 1:n_bins
                in_bin = cycle1_phase >= phase_bin_edges(bi) & cycle1_phase < phase_bin_edges(bi+1);
                if any(in_bin)
                    cycle1_power_sum(:, bi) = cycle1_power_sum(:, bi) + mean(cycle1_spec(:, in_bin), 2);
                    cycle1_count(bi) = cycle1_count(bi) + 1;
                end
            end
        end
        
        if ~isempty(cycle2_indices)
            cycle2_phase = epoch_phase(cycle2_indices);
            cycle2_spec = epoch_spec(:, cycle2_indices);
            for bi = 1:n_bins
                in_bin = cycle2_phase >= phase_bin_edges(bi) & cycle2_phase < phase_bin_edges(bi+1);
                if any(in_bin)
                    cycle2_power_sum(:, bi) = cycle2_power_sum(:, bi) + mean(cycle2_spec(:, in_bin), 2);
                    cycle2_count(bi) = cycle2_count(bi) + 1;
                end
            end
        end
    end
    
    mean_spectrogram = zeros(n_freqs, n_bins);
    for bi = 1:n_bins
        if bin_count(bi) > 0
            mean_spectrogram(:, bi) = bin_power_sum(:, bi) / bin_count(bi);
        end
    end
    
    cycle1_mean_spectrogram = zeros(n_freqs, n_bins);
    for bi = 1:n_bins
        if cycle1_count(bi) > 0
            cycle1_mean_spectrogram(:, bi) = cycle1_power_sum(:, bi) / cycle1_count(bi);
        end
    end
    
    cycle2_mean_spectrogram = zeros(n_freqs, n_bins);
    for bi = 1:n_bins
        if cycle2_count(bi) > 0
            cycle2_mean_spectrogram(:, bi) = cycle2_power_sum(:, bi) / cycle2_count(bi);
        end
    end
    
    fprintf(' Done.\n');
end

function fiber_roi_map = build_fiber_roi_map(FPA)
    fiber_roi_map = {};
    
    if ~isfield(FPA, 'rois') || ~isfield(FPA.rois, 'roi_data')
        if isfield(FPA, 'parameters') && isfield(FPA.parameters, 'num_fibers')
            num_fibers = FPA.parameters.num_fibers;
        elseif isfield(FPA, 'signals') && isfield(FPA.signals, 'final_processed_traces')
            num_fibers = size(FPA.signals.final_processed_traces, 2);
        else
            num_fibers = 1;
        end
        
        for k = 1:num_fibers
            fiber_roi_map{k} = sprintf('FOV1_ROI%d', k);
        end
        return;
    end
    
    roi_data = FPA.rois.roi_data;
    fiber_idx = 1;
    
    for fov_idx = 1:length(roi_data)
        fov_rois = roi_data{fov_idx};
        if isempty(fov_rois)
            continue;
        end
        
        for roi_idx = 1:length(fov_rois)
            fiber_roi_map{fiber_idx} = sprintf('FOV%d_ROI%d', fov_idx, roi_idx);
            fiber_idx = fiber_idx + 1;
        end
    end
end


function artifact_mask = load_artifact_mask(artifact_masks_dir, base_filename, filename)
% LOAD_ARTIFACT_MASK  Find and load the artifact mask for a single trial.
%
%   base_filename : e.g. 'Animal01-01_09_25-R1_Trial1'
%   filename      : e.g. 'Animal01-01_09_25-R1_Trial1_FiberPhotometry_Analysis.mat'
%
%   Matching logic:
%     1. Strip '_Trial\d+' suffix from base_filename to get session prefix.
%     2. Look for <session_prefix>_artifact_removal.mat in artifact_masks_dir.
%     3. Inside the file, find the trial whose .filename matches <filename>.
%     4. Return that trial's artifact_mask as a logical column vector.
%        If no match is found at any step, return [].

    artifact_mask = [];

    if ~exist(artifact_masks_dir, 'dir')
        return;
    end

    % Derive session prefix: 'Animal01-01_09_25-R1_Trial1' -> 'Animal01-01_09_25-R1'
    session_prefix = regexprep(base_filename, '_Trial\d+$', '');
    mask_file = fullfile(artifact_masks_dir, [session_prefix '_artifact_removal.mat']);

    if ~exist(mask_file, 'file')
        return;
    end

    try
        art_data = load(mask_file, 'ArtifactInfo');
        if ~isfield(art_data, 'ArtifactInfo')
            return;
        end
        AI = art_data.ArtifactInfo;
        for ti = 1:length(AI.trials)
            if strcmp(AI.trials(ti).filename, filename)
                artifact_mask = logical(AI.trials(ti).artifact_mask(:));
                return;
            end
        end
    catch ME
        warning('compute_phase_aligned_spectrogram_batch:artifactLoad', ...
                'Could not load artifact mask from %s: %s', mask_file, ME.message);
    end
end


function lbl = lfp_wavelet_colorbar_ylabel(lfp_field_name)
%LFP_WAVELET_COLORBAR_YLABEL  Y-label for colorbars on LFP wavelet magnitude panels.
    if isempty(lfp_field_name)
        lbl = 'Wavelet mag. (mV)';
        return;
    end
    if strcmp(lfp_field_name, 'lfp_z_HP') || strcmp(lfp_field_name, 'lfp_z') ...
            || strncmp(lfp_field_name, 'lfp_z_', 6)
        lbl = 'Wavelet mag. (a.u., z-scored)';
    else
        lbl = 'Wavelet mag. (mV)';
    end
end


function all_epochs_raw = annotate_epochs_artifact(all_epochs_raw, artifact_mask)
% ANNOTATE_EPOCHS_ARTIFACT  Add has_artifact and artifact_fraction fields to
% each epoch in all_epochs_raw using the epoch's time_indices.
%
%   has_artifact      - true if ANY sample within the epoch is flagged
%   artifact_fraction - fraction (0-1) of samples flagged as artifact
%
% If artifact_mask is empty, both fields are set to false/0.

    for ei = 1:length(all_epochs_raw)
        if isempty(artifact_mask)
            all_epochs_raw{ei}.has_artifact      = false;
            all_epochs_raw{ei}.artifact_fraction = 0;
        else
            si  = all_epochs_raw{ei}.time_indices(1);
            ei_ = all_epochs_raw{ei}.time_indices(2);
            % Guard against artifact_mask being shorter than signal
            if ei_ <= length(artifact_mask)
                epoch_art = artifact_mask(si:ei_);
            elseif si <= length(artifact_mask)
                epoch_art = artifact_mask(si:end);
            else
                epoch_art = false(1,1);
            end
            all_epochs_raw{ei}.has_artifact      = any(epoch_art);
            all_epochs_raw{ei}.artifact_fraction = mean(double(epoch_art));
        end
    end
end

