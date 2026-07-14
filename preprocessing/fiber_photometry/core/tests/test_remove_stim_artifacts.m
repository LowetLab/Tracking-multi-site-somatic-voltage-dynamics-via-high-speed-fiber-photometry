function test_remove_stim_artifacts()
%TEST_REMOVE_STIM_ARTIFACTS  Unit test for core/remove_stim_artifacts.
%   Covers: bit-identical equivalence to the original inline loop, output
%   shape, multi-fiber independence, the >90-sample end-replacement guard
%   (short traces are left un-replaced), and that in-band (120-132 Hz) power
%   is reduced.
%
%   Run:  test_remove_stim_artifacts

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
nfail = 0;
fs = 500;

% --- equivalence with the verbatim inline loop (multi-fiber) ------------
rng(13);
t = (0:1/fs:6)';
traces = 0.5*randn(numel(t), 3) ...
    + sin(2*pi*122*t) + sin(2*pi*131*t);     % 122 & 131 Hz artifacts in-band
got = remove_stim_artifacts(traces, fs);
ref = ref_inline(traces, fs);
nfail = nfail + expect_eq('matches inline (bit-identical)', got, ref);
nfail = nfail + expect_eq('output shape', size(got), size(traces));

% --- per-fiber independence ---------------------------------------------
nfail = nfail + expect_eq('fiber 2 independent', ...
    remove_stim_artifacts(traces(:,2), fs), ref_inline(traces(:,2), fs));

% --- short trace (<=90 samples): end-replacement guard skips ------------
short = randn(50, 1);
gs = remove_stim_artifacts(short, fs);
nfail = nfail + expect_eq('short trace matches inline', gs, ref_inline(short, fs));

% --- artifact power is actually reduced ---------------------------------
pure = sin(2*pi*122*t);
cleaned = remove_stim_artifacts(pure, fs);
nfail = nfail + expect_true('122 Hz power reduced', rms(cleaned(100:end-100)) < 0.3*rms(pure));

if nfail == 0
    fprintf('\nALL remove_stim_artifacts TESTS PASSED\n');
else
    error('test_remove_stim_artifacts:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect_eq(name, got, want)
if isequal(got, want), nf = 0; else
    fprintf('  FAIL %s\n', name); nf = 1;
end
end

function nf = expect_true(name, cond)
if cond, nf = 0; else
    fprintf('  FAIL %s\n', name); nf = 1;
end
end

% Verbatim original inline SECTION 4 loop.
function filtered_traces = ref_inline(all_traces, IMAGING_FS)
filtered_traces = zeros(size(all_traces));
for trace_idx = 1:size(all_traces, 2)
    current_trace = all_traces(:, trace_idx);
    if size(current_trace, 1) > 90
        current_trace(end-80:end) = median(current_trace(1:end-90));
    end
    lowCutoff = 120; highCutoff = 124; filterOrder = 3;
    [b, a] = butter(filterOrder, [lowCutoff, highCutoff] / (IMAGING_FS / 2), 'bandpass');
    artifact_120Hz = filtfilt(b, a, current_trace);
    current_trace = current_trace - artifact_120Hz;
    lowCutoff = 130; highCutoff = 132;
    [b, a] = butter(filterOrder, [lowCutoff, highCutoff] / (IMAGING_FS / 2), 'bandpass');
    artifact_130Hz = filtfilt(b, a, current_trace);
    current_trace = current_trace - artifact_130Hz;
    filtered_traces(:, trace_idx) = current_trace;
end
end
