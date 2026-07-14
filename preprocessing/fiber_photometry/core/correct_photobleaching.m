function [traces_detrended, traces_exp_corrected, fit_diag] = correct_photobleaching( ...
        filtered_traces, stim_onset_frame, is_baseline_trial, imaging_fs, opts)
%CORRECT_PHOTOBLEACHING  Two photobleaching corrections for fiber traces.
%   [traces_detrended, traces_exp_corrected, fit_diag] = CORRECT_PHOTOBLEACHING( ...
%       filtered_traces, stim_onset_frame, is_baseline_trial, imaging_fs)
%   corrects each column (trace) of FILTERED_TRACES for photobleaching using
%   the pre-stimulation baseline, by two methods:
%
%     Method 1 (linear)    : subtract a straight-line fit of the baseline
%                            -> traces_detrended  (F_corr, NOT yet dF/F).
%     Method 2 (double-exp): divide by a double-exponential fit of the baseline
%                            -> traces_exp_corrected (F_corr, NOT yet dF/F).
%
%   The baseline period is frames 1..stim_onset_frame (whole trace if
%   is_baseline_trial is true) or 1..stim_onset_frame-1 otherwise.
%
%   Inputs
%   ------
%   filtered_traces   : frames x nTraces matrix of (motion/spike) filtered traces.
%   stim_onset_frame  : first stimulation frame (baseline = frames before it).
%   is_baseline_trial : true to treat the whole trace as baseline.
%   imaging_fs        : imaging sampling rate (Hz) for the exponential time axis.
%
%   Outputs
%   -------
%   traces_detrended      : linear-detrended traces (same size as input).
%   traces_exp_corrected  : double-exponential-corrected traces (same size).
%   fit_diag              : struct of "last successful iteration" fit values that
%                           the original inline code left in the workspace and
%                           that downstream figures + the saved data struct read.
%                           Fields (mirroring the original leftover semantics):
%                             .double_exp_function : the model handle (always)
%                             .pre_stim_trace      : last trace's baseline (always)
%                             .time_pre_stim       : present iff a baseline >= 3
%                             .time_full           : present iff a baseline >= 3
%                             .fitted_params       : present iff a fit succeeded
%                           Restore these into the workspace after calling, e.g.
%                             double_exp_function = fit_diag.double_exp_function;
%                             pre_stim_trace      = fit_diag.pre_stim_trace;
%                             if isfield(fit_diag,'fitted_params'), fitted_params = fit_diag.fitted_params; end
%                           Using isfield reproduces the original exist('var')
%                           guards exactly.
%
%   Notes
%   -----
%   * dF/F is NOT computed here -- these are corrected fluorescence (F_corr).
%   * If a baseline is shorter than 3 frames, or the exponential fit fails, that
%     trace falls back to the linear-detrended result (with a warning), exactly
%     as the original inline code did.
%
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts.
%   Unit-tested in core/tests/test_correct_photobleaching.m.

% Optional behaviour flags. DEFAULTS reproduce the MULTI-TRIAL script exactly;
% the single-trial script passes both as false (its historical behaviour).
%   opts.BaselineClampMin   (default true)  -> Method-2 baseline index uses
%                                              max(1, stim_onset_frame-1).
%   opts.ShortBaselineGuard (default true)  -> a baseline < 3 frames falls back
%                                              to linear detrend (skips exp fit).
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'BaselineClampMin'),   opts.BaselineClampMin = true;   end
if ~isfield(opts, 'ShortBaselineGuard'), opts.ShortBaselineGuard = true; end

% Initialize corrected traces (start from the filtered traces)
traces_detrended = filtered_traces;
traces_exp_corrected = filtered_traces;

% Double exponential model: A1*exp(-t/tau1) + A2*exp(-t/tau2) + offset
double_exp_function = @(params, t) params(1) * exp(-t/params(2)) + ...
    params(3) * exp(-t/params(4)) + params(5);

% fit_diag carries the loop's "leftover" variables (see help). Only the model
% handle is guaranteed; the rest are added under the same conditions the
% original inline code would have set them, so isfield(...) mirrors exist(...).
fit_diag = struct('double_exp_function', double_exp_function);

for trace_idx = 1:size(filtered_traces, 2)
    current_trace = filtered_traces(:, trace_idx);

    %% METHOD 1: Linear detrending
    % For photobleaching correction, use pre-stimulation period only
    if is_baseline_trial
        pre_stim_indices = 1:stim_onset_frame;  % Use entire trace for baseline
    else
        pre_stim_indices = 1:stim_onset_frame-1;  % Exclude stim onset frame
    end
    pre_stim_time = pre_stim_indices;
    pre_stim_values = current_trace(pre_stim_indices);

    trend_coefficients = polyfit(pre_stim_time, pre_stim_values, 1);
    full_time_indices = 1:size(filtered_traces, 1);
    trend_line = polyval(trend_coefficients, full_time_indices);
    detrended_trace = current_trace - trend_line';

    % Store detrended trace WITHOUT normalization (F_corr for linear method)
    % Normalization will be done in dF/F computation step
    traces_detrended(:, trace_idx) = detrended_trace;

    %% METHOD 2: Double exponential fitting
    % For photobleaching correction, use pre-stimulation period only.
    % For baseline trials, stim_onset_frame = end of trace, so use entire trace.
    if is_baseline_trial
        pre_stim_trace = current_trace(1:stim_onset_frame);  % Use entire trace for baseline
    elseif opts.BaselineClampMin
        pre_stim_trace = current_trace(1:max(1, stim_onset_frame-1));  % Exclude stim onset frame, ensure at least 1 frame
    else
        pre_stim_trace = current_trace(1:stim_onset_frame-1);  % Exclude stim onset frame (single-trial: no clamp)
    end
    fit_diag.pre_stim_trace = pre_stim_trace;  % set every iteration (as inline code did)

    % Safety check: ensure we have enough data points
    if opts.ShortBaselineGuard && length(pre_stim_trace) < 3
        warning('Pre-stim trace too short (%d frames) for double exponential fit. Using linear detrending only.', ...
            length(pre_stim_trace));
        % Skip double exponential fitting, use only linear detrending
        traces_exp_corrected(:, trace_idx) = detrended_trace;
        continue;
    end

    time_pre_stim = (0:length(pre_stim_trace)-1)' / imaging_fs;
    time_full = (0:size(filtered_traces, 1)-1)' / imaging_fs;
    fit_diag.time_pre_stim = time_pre_stim;
    fit_diag.time_full = time_full;

    initial_params = [0.5*max(pre_stim_trace), time_pre_stim(end)/3, ...
        0.5*max(pre_stim_trace), time_pre_stim(end), min(pre_stim_trace)];
    lower_bounds = [0, 0, 0, 0, 0];
    upper_bounds = [Inf, Inf, Inf, Inf, Inf];
    fit_options = optimoptions('lsqcurvefit', 'Display', 'off');

    try
        fitted_params = lsqcurvefit(double_exp_function, initial_params, ...
            time_pre_stim, pre_stim_trace, ...
            lower_bounds, upper_bounds, fit_options);

        fitted_curve_full = double_exp_function(fitted_params, time_full);
        % Photobleaching correction: divide by fitted curve to get F_corr(t)
        % Do NOT compute dF/F here - that will be done in the next step
        corrected_trace = current_trace ./ fitted_curve_full;

        traces_exp_corrected(:, trace_idx) = corrected_trace;
        fit_diag.fitted_params = fitted_params;  % set only on a successful fit
    catch ME
        warning('Exponential fitting failed for trace %d: %s', trace_idx, ME.message);
        traces_exp_corrected(:, trace_idx) = traces_detrended(:, trace_idx);
    end
end
end
