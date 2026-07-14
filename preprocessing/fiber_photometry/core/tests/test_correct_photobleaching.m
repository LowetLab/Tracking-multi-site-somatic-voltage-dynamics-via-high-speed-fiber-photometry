function test_correct_photobleaching()
%TEST_CORRECT_PHOTOBLEACHING  Equivalence test for core/correct_photobleaching.
%   Compares the extracted function against a VERBATIM inline copy of the
%   original loop on synthetic photobleached traces. Because both run identical
%   code (including the lsqcurvefit path, or its catch-fallback if the
%   Optimization Toolbox is absent), the corrected traces must be bit-for-bit
%   identical. The test ALSO checks that fit_diag reproduces exactly the
%   "leftover" variables the original loop left in scope (which downstream
%   figures and the saved data struct read), by capturing them in the reference
%   with the same exist(...) guards the main script uses.
%
%   Run:  test_correct_photobleaching

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
rng(7);
nfail = 0;

nframes = 600;
t = (0:nframes-1)';
make = @(A,tau,off,nz) A*exp(-t/tau) + off + nz*randn(nframes,1);
F = [make(100,150,20,0.5), make(60,300,15,0.5), make(140,90,25,0.5)];
imaging_fs = 650;

% --- case A: stimulation trial -------------------------------------------
[d,e,g]    = correct_photobleaching(F, 400, false, imaging_fs);
[dr,er,gr] = ref(F, 400, false, imaging_fs);
nfail = nfail + chk('stim: detrended', d, dr) + chk('stim: exp', e, er) + chkdiag('stim: diag', g, gr);

% --- case B: baseline trial ----------------------------------------------
[d,e,g]    = correct_photobleaching(F, nframes, true, imaging_fs);
[dr,er,gr] = ref(F, nframes, true, imaging_fs);
nfail = nfail + chk('baseline: detrended', d, dr) + chk('baseline: exp', e, er) + chkdiag('baseline: diag', g, gr);

% --- case C: very short baseline -> linear-only fallback -----------------
ws = warning('off', 'all');
[d,e,g]    = correct_photobleaching(F, 3, false, imaging_fs);
[dr,er,gr] = ref(F, 3, false, imaging_fs);
warning(ws);
nfail = nfail + chk('short: detrended', d, dr) + chk('short: exp', e, er) + chkdiag('short: diag', g, gr);

% --- single-trial variant: no baseline clamp, no short-baseline guard ----
opts_s = struct('BaselineClampMin', false, 'ShortBaselineGuard', false);
[d,e,g]    = correct_photobleaching(F, 400, false, imaging_fs, opts_s);
[dr,er,gr] = ref_single(F, 400, false, imaging_fs);
nfail = nfail + chk('single normal: detrended', d, dr) + chk('single normal: exp', e, er) + chkdiag('single normal: diag', g, gr);

ws2 = warning('off', 'all');
[d,e,g]    = correct_photobleaching(F, 2, false, imaging_fs, opts_s);   % 1-frame baseline, NO guard
[dr,er,gr] = ref_single(F, 2, false, imaging_fs);
warning(ws2);
nfail = nfail + chk('single short(no guard): detrended', d, dr) + chk('single short(no guard): exp', e, er) + chkdiag('single short(no guard): diag', g, gr);

if nfail == 0
    fprintf('\nALL correct_photobleaching TESTS PASSED\n');
else
    error('test_correct_photobleaching:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = chk(name, got, want)
if isequaln(got, want)
    nf = 0;
else
    fprintf('  FAIL %s\n', name); nf = 1;
end
end

function nf = chkdiag(name, got, want)
nf = 0;
fg = sort(fieldnames(got)); fw = sort(fieldnames(want));
if ~isequal(fg, fw)
    fprintf('  FAIL %s: fields {%s} vs {%s}\n', name, strjoin(fg,','), strjoin(fw,','));
    nf = 1; return;
end
for i = 1:numel(fg)
    f = fg{i};
    if isa(got.(f), 'function_handle')
        oktype = strcmp(func2str(got.(f)), func2str(want.(f)));
    else
        oktype = isequaln(got.(f), want.(f));
    end
    if ~oktype
        fprintf('  FAIL %s: field %s differs\n', name, f);
        nf = nf + 1;
    end
end
end

% =========================================================================
% Reference: VERBATIM original loop, then capture the genuine leftover state
% exactly as the main script reads it (exist guards) -> ground-truth fit_diag.
function [traces_detrended, traces_exp_corrected, refdiag] = ref(filtered_traces, stim_onset_frame, is_baseline_trial, imaging_fs)
traces_detrended = filtered_traces;
traces_exp_corrected = filtered_traces;
double_exp_function = @(params, t) params(1) * exp(-t/params(2)) + ...
    params(3) * exp(-t/params(4)) + params(5);
for trace_idx = 1:size(filtered_traces, 2)
    current_trace = filtered_traces(:, trace_idx);
    if is_baseline_trial
        pre_stim_indices = 1:stim_onset_frame;
    else
        pre_stim_indices = 1:stim_onset_frame-1;
    end
    pre_stim_time = pre_stim_indices;
    pre_stim_values = current_trace(pre_stim_indices);
    trend_coefficients = polyfit(pre_stim_time, pre_stim_values, 1);
    full_time_indices = 1:size(filtered_traces, 1);
    trend_line = polyval(trend_coefficients, full_time_indices);
    detrended_trace = current_trace - trend_line';
    traces_detrended(:, trace_idx) = detrended_trace;
    if is_baseline_trial
        pre_stim_trace = current_trace(1:stim_onset_frame);
    else
        pre_stim_trace = current_trace(1:max(1, stim_onset_frame-1));
    end
    if length(pre_stim_trace) < 3
        traces_exp_corrected(:, trace_idx) = detrended_trace;
        continue;
    end
    time_pre_stim = (0:length(pre_stim_trace)-1)' / imaging_fs;
    time_full = (0:size(filtered_traces, 1)-1)' / imaging_fs;
    initial_params = [0.5*max(pre_stim_trace), time_pre_stim(end)/3, ...
        0.5*max(pre_stim_trace), time_pre_stim(end), min(pre_stim_trace)];
    lower_bounds = [0, 0, 0, 0, 0];
    upper_bounds = [Inf, Inf, Inf, Inf, Inf];
    fit_options = optimoptions('lsqcurvefit', 'Display', 'off');
    try
        fitted_params = lsqcurvefit(double_exp_function, initial_params, ...
            time_pre_stim, pre_stim_trace, lower_bounds, upper_bounds, fit_options);
        fitted_curve_full = double_exp_function(fitted_params, time_full);
        traces_exp_corrected(:, trace_idx) = current_trace ./ fitted_curve_full;
    catch
        traces_exp_corrected(:, trace_idx) = traces_detrended(:, trace_idx);
    end
end
% Capture leftover state the same way the main script does (exist guards):
refdiag = struct('double_exp_function', double_exp_function);
if exist('pre_stim_trace', 'var'), refdiag.pre_stim_trace = pre_stim_trace; end
if exist('time_pre_stim', 'var'),  refdiag.time_pre_stim  = time_pre_stim;  end
if exist('time_full', 'var'),      refdiag.time_full      = time_full;      end
if exist('fitted_params', 'var'),  refdiag.fitted_params  = fitted_params;  end
end

% =========================================================================
% Reference: VERBATIM original SINGLE-TRIAL loop (no baseline clamp, no
% short-baseline guard), then capture genuine leftover state via exist guards.
function [traces_detrended, traces_exp_corrected, refdiag] = ref_single(filtered_traces, stim_onset_frame, is_baseline_trial, imaging_fs)
traces_detrended = filtered_traces;
traces_exp_corrected = filtered_traces;
double_exp_function = @(params, t) params(1) * exp(-t/params(2)) + ...
    params(3) * exp(-t/params(4)) + params(5);
for trace_idx = 1:size(filtered_traces, 2)
    current_trace = filtered_traces(:, trace_idx);
    if is_baseline_trial
        pre_stim_indices = 1:stim_onset_frame;
    else
        pre_stim_indices = 1:stim_onset_frame-1;
    end
    pre_stim_time = pre_stim_indices;
    pre_stim_values = current_trace(pre_stim_indices);
    trend_coefficients = polyfit(pre_stim_time, pre_stim_values, 1);
    full_time_indices = 1:size(filtered_traces, 1);
    trend_line = polyval(trend_coefficients, full_time_indices);
    detrended_trace = current_trace - trend_line';
    traces_detrended(:, trace_idx) = detrended_trace;
    if is_baseline_trial
        pre_stim_trace = current_trace(1:stim_onset_frame);
    else
        pre_stim_trace = current_trace(1:stim_onset_frame-1);   % no max(1,...) clamp
    end
    time_pre_stim = (0:length(pre_stim_trace)-1)' / imaging_fs;  % no <3 guard
    time_full = (0:size(filtered_traces, 1)-1)' / imaging_fs;
    initial_params = [0.5*max(pre_stim_trace), time_pre_stim(end)/3, ...
        0.5*max(pre_stim_trace), time_pre_stim(end), min(pre_stim_trace)];
    lower_bounds = [0, 0, 0, 0, 0];
    upper_bounds = [Inf, Inf, Inf, Inf, Inf];
    fit_options = optimoptions('lsqcurvefit', 'Display', 'off');
    try
        fitted_params = lsqcurvefit(double_exp_function, initial_params, ...
            time_pre_stim, pre_stim_trace, lower_bounds, upper_bounds, fit_options);
        fitted_curve_full = double_exp_function(fitted_params, time_full);
        traces_exp_corrected(:, trace_idx) = current_trace ./ fitted_curve_full;
    catch
        traces_exp_corrected(:, trace_idx) = traces_detrended(:, trace_idx);
    end
end
refdiag = struct('double_exp_function', double_exp_function);
if exist('pre_stim_trace', 'var'), refdiag.pre_stim_trace = pre_stim_trace; end
if exist('time_pre_stim', 'var'),  refdiag.time_pre_stim  = time_pre_stim;  end
if exist('time_full', 'var'),      refdiag.time_full      = time_full;      end
if exist('fitted_params', 'var'),  refdiag.fitted_params  = fitted_params;  end
end
