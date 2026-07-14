%% VALIDATE_SPIKE_ML_DATA
%  Quick visual sanity check of the SpikeMLData structure.
%  Generates three figures:
%    1. Grid of mean waveforms per neuron
%    2. Individual spike overlay for a selected neuron
%    3. Full trace with spike locations marked for a selected neuron
%
%  USAGE:
%    1. Set ML_DATA_FILE below
%    2. Run the script
%
clear; clc; close all;

%% Configuration  -- EDIT THIS FOR YOUR SESSION
ML_DATA_FILE = 'C:\PATH\TO\YOUR\CellularDataProcessed\01-01-25-R1\Animal01_01-01-25-R1_SpikeMLData.mat';

fprintf('Loading SpikeMLData from:\n  %s\n', ML_DATA_FILE);
S = load(ML_DATA_FILE, 'SpikeMLData');
D = S.SpikeMLData;

info    = D.info;
trial   = D.trials{1};
neurons = trial.neurons;
pooled  = D.snippets_pooled;

num_neurons    = info.num_neurons;
snippet_time   = info.snippet_time_axis_ms;
fs             = info.sampling_rate_hz;

fprintf('  %d neurons, %d total spikes, snippet = %d samples (%.0f ms)\n', ...
    num_neurons, info.total_spikes, info.snippet_length, ...
    1000 * info.snippet_length / fs);

%% ========================================================================
%  FIGURE 1: Mean waveform per neuron (grid)
%  ========================================================================

n_cols = ceil(sqrt(num_neurons));
n_rows = ceil(num_neurons / n_cols);

fig1 = figure('Name', 'Mean Waveforms', 'Color', 'w', 'Position', [80, 80, 1100, 750]);

for n = 1:num_neurons
    subplot(n_rows, n_cols, n);
    hold on;

    if neurons(n).num_spikes > 0
        % Light overlay of individual snippets (up to 30)
        snips = neurons(n).spike_snippets;
        n_show = min(size(snips, 1), 30);
        for k = 1:n_show
            plot(snippet_time, snips(k, :), 'Color', [0.75, 0.75, 0.75, 0.4], 'LineWidth', 0.5);
        end
        % Mean waveform
        plot(snippet_time, neurons(n).mean_waveform, 'k-', 'LineWidth', 1.8);
    end

    xline(0, 'r--', 'LineWidth', 1);
    title(sprintf('N%d (%d)', n, neurons(n).num_spikes), 'FontSize', 9);
    if n > (n_rows-1)*n_cols, xlabel('ms'); end
    if mod(n-1, n_cols) == 0, ylabel('\DeltaF/F'); end
    set(gca, 'FontSize', 8);
    xlim([snippet_time(1), snippet_time(end)]);
    grid on; box on;
end

sgtitle(sprintf('%s %s-%s: Mean Spike Waveforms (total %d spikes)', ...
    info.mouse_name, info.recording_date, info.recording_id, info.total_spikes), ...
    'FontSize', 13, 'FontWeight', 'bold');

%% ========================================================================
%  FIGURE 2: Detailed single-neuron spike overlay
%  Pick the neuron with the most spikes for a rich display
%  ========================================================================

[~, best_neuron] = max([neurons.num_spikes]);
bn = neurons(best_neuron);

fig2 = figure('Name', sprintf('Neuron %d Spike Overlay', best_neuron), ...
    'Color', 'w', 'Position', [100, 100, 800, 500]);

snips = bn.spike_snippets;
n_spikes = size(snips, 1);

% Color-code by SNR
snr_vals = bn.spike_snr;
if all(isnan(snr_vals))
    colors = repmat([0.6, 0.6, 0.6], n_spikes, 1);
    has_snr = false;
else
    cmap = parula(256);
    snr_norm = (snr_vals - min(snr_vals)) / max(1e-9, max(snr_vals) - min(snr_vals));
    color_idx = max(1, min(256, round(snr_norm * 255) + 1));
    colors = cmap(color_idx, :);
    has_snr = true;
end

hold on;
for k = 1:n_spikes
    plot(snippet_time, snips(k, :), 'Color', [colors(k,:), 0.5], 'LineWidth', 0.8);
end
plot(snippet_time, bn.mean_waveform, 'k-', 'LineWidth', 2.5);
xline(0, 'r--', 'LineWidth', 1.5);

xlabel('Time relative to spike peak (ms)', 'FontSize', 12);
ylabel('\DeltaF/F', 'FontSize', 12);
title(sprintf('Neuron %d: %d spikes overlaid (%.1f Hz)', best_neuron, n_spikes, bn.firing_rate_hz), ...
    'FontSize', 13);
set(gca, 'FontSize', 11);
grid on; box on;

if has_snr
    colormap(parula);
    cb = colorbar;
    caxis([min(snr_vals), max(snr_vals)]);
    cb.Label.String = 'Spike SNR';
    cb.Label.FontSize = 11;
end

%% ========================================================================
%  FIGURE 3: Full trace with spike locations for the same neuron
%  ========================================================================

fig3 = figure('Name', sprintf('Neuron %d Full Trace', best_neuron), ...
    'Color', 'w', 'Position', [120, 120, 1200, 400]);

trace = trial.traces.corrected(:, best_neuron);
time_sec = trial.traces.time_vector;

plot(time_sec, trace, 'Color', [0.3, 0.3, 0.3], 'LineWidth', 0.8);
hold on;

spike_idx = bn.spike_indices;
if ~isempty(spike_idx)
    plot(time_sec(spike_idx), trace(spike_idx), 'rv', ...
        'MarkerSize', 6, 'MarkerFaceColor', [1, 0.3, 0.3]);
end

xlabel('Time (s)', 'FontSize', 12);
ylabel('\DeltaF/F', 'FontSize', 12);
title(sprintf('Neuron %d: Full corrected trace with %d detected spikes', best_neuron, n_spikes), ...
    'FontSize', 13);
set(gca, 'FontSize', 11);
grid on; box on;
legend('Corrected trace', 'Detected spike peaks', 'Location', 'best');

%% ========================================================================
%  FIGURE 4: Pooled snippet matrix overview
%  ========================================================================

fig4 = figure('Name', 'Pooled Snippets Heatmap', 'Color', 'w', 'Position', [140, 140, 900, 500]);

% Sort pooled snippets by neuron id then by time
[~, sort_order] = sortrows([pooled.neuron_id, pooled.spike_time_sec]);
sorted_data = pooled.data(sort_order, :);
sorted_nids = pooled.neuron_id(sort_order);

imagesc(snippet_time, 1:size(sorted_data,1), sorted_data);
colormap(parula); colorbar;
xlabel('Time relative to peak (ms)', 'FontSize', 12);
ylabel('Spike index (sorted by neuron)', 'FontSize', 12);
title(sprintf('All %d spike snippets (sorted by neuron)', info.total_spikes), 'FontSize', 13);
set(gca, 'FontSize', 11);

% Draw horizontal lines at neuron boundaries
hold on;
for n = 1:num_neurons-1
    boundary = find(sorted_nids == n, 1, 'last');
    if ~isempty(boundary)
        yline(boundary + 0.5, 'w-', 'LineWidth', 1);
    end
end

fprintf('\nValidation complete. Four figures generated.\n');
fprintf('  Fig 1: Mean waveforms grid\n');
fprintf('  Fig 2: Single-neuron spike overlay (N%d, %d spikes)\n', best_neuron, n_spikes);
fprintf('  Fig 3: Full trace with spike markers (N%d)\n', best_neuron);
fprintf('  Fig 4: Pooled snippets heatmap\n');
