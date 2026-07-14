function [traces_deltaF_F, F0_values, processed_traces_zscored] = compute_deltaF_F(traces_exp_corrected, baseline_start_frame, baseline_end_frame)
%COMPUTE_DELTAF_F  Static-baseline ΔF/F (and z-scored ΔF/F) per fiber.
%   [traces_deltaF_F, F0_values, processed_traces_zscored] =
%   COMPUTE_DELTAF_F(traces_exp_corrected, baseline_start_frame, baseline_end_frame)
%
%   For each column (fiber) of the photobleaching-corrected fluorescence
%   F_corr(t), F0 is the mean over the baseline window
%   [baseline_start_frame, baseline_end_frame], and
%
%       ΔF/F = (F_corr - F0) / F0
%
%   The z-scored ΔF/F is also returned (kept for backward compatibility).
%
%   Inputs
%   ------
%   traces_exp_corrected : frames x fibers matrix of corrected fluorescence.
%   baseline_start_frame : first frame of the F0 baseline window (inclusive).
%   baseline_end_frame   : last frame of the F0 baseline window (inclusive).
%
%   Outputs
%   -------
%   traces_deltaF_F         : frames x fibers ΔF/F traces.
%   F0_values               : 1 x fibers baseline F0 per fiber.
%   processed_traces_zscored: frames x fibers z-scored ΔF/F.
%
%   The CHOICE of baseline window (baseline trial / automatic stim onset /
%   fallback / 80%) stays in the calling script; this function only applies a
%   window it is handed. Extracted (behaviour-preserving) from the fiber
%   preprocessing scripts, where the per-fiber loop was inlined identically in
%   both the multi- and single-trial versions. Unit-tested in
%   core/tests/test_compute_deltaF_F.m.

% Compute ΔF/F for each fiber
% Input: traces_exp_corrected contains F_corr(t) (photobleaching-corrected fluorescence)
% Output: traces_deltaF_F contains ΔF/F = (F_corr - F0) / F0
traces_deltaF_F = zeros(size(traces_exp_corrected));
F0_values = zeros(1, size(traces_exp_corrected, 2));  % Store F0 values for each fiber
for fiber_idx = 1:size(traces_exp_corrected, 2)
    fiber_corrected = traces_exp_corrected(:, fiber_idx);  % F_corr(t)

    % Calculate F0 as mean of baseline window from corrected trace
    baseline_window = fiber_corrected(baseline_start_frame:baseline_end_frame);
    F0 = mean(baseline_window);
    F0_values(fiber_idx) = F0;  % Store F0 value

    % Compute ΔF/F = (F_corr - F0) / F0
    traces_deltaF_F(:, fiber_idx) = (fiber_corrected - F0) / F0;

    fprintf('  Fiber %d: F0 = %.4f, ΔF/F range = [%.4f, %.4f]\n', ...
        fiber_idx, F0, min(traces_deltaF_F(:, fiber_idx)), max(traces_deltaF_F(:, fiber_idx)));
end

fprintf('ΔF/F computation complete\n');

% Z-score normalization (preserved for backward compatibility)
% Z-score the ΔF/F traces, not F_corr
processed_traces_zscored = zeros(size(traces_deltaF_F));
for i = 1:size(traces_deltaF_F, 2)
    processed_traces_zscored(:, i) = zscore(traces_deltaF_F(:, i));
end
end
