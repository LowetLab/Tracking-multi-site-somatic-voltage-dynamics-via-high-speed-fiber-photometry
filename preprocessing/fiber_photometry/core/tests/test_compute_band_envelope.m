function test_compute_band_envelope()
%TEST_COMPUTE_BAND_ENVELOPE  Unit test for core/compute_band_envelope.
%   Covers: bit-identical equivalence to the original inline block
%   (butter -> filtfilt -> abs(hilbert) -> fastsmooth), output shapes, and a
%   sanity check that an in-band tone yields a roughly constant envelope.
%
%   Run:  test_compute_band_envelope

addpath(fileparts(fileparts(mfilename('fullpath'))));                       % core/
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));    % fastsmooth
nfail = 0;

fs = 500; band = [5 10]; order = 3; win = 90;
t = (0:1/fs:20)';
rng(9);
sig = sin(2*pi*7*t) + 0.5*randn(size(t));   % 7 Hz (in band) + noise

% --- equivalence with the verbatim inline block -------------------------
[f1, e1, s1] = compute_band_envelope(sig, fs, band, order, win);
[f0, e0, s0] = ref_inline(sig, fs, band, order, win);
nfail = nfail + expect_eq('filtered matches inline', f1, f0);
nfail = nfail + expect_eq('envelope matches inline', e1, e0);
nfail = nfail + expect_eq('smoothed envelope matches inline', s1, s0);

% --- shapes preserved ----------------------------------------------------
nfail = nfail + expect_eq('filtered shape', size(f1), size(sig));
nfail = nfail + expect_eq('smoothed shape', size(s1), size(sig));

% --- sanity: in-band tone -> roughly constant envelope (~amplitude 1) ---
pure = sin(2*pi*7*t);
[~, ~, env_s] = compute_band_envelope(pure, fs, band, order, win);
mid = env_s(round(0.2*end):round(0.8*end));   % drop edge transients
nfail = nfail + expect_true('in-band envelope near 1', abs(median(mid) - 1) < 0.15);

if nfail == 0
    fprintf('\nALL compute_band_envelope TESTS PASSED\n');
else
    error('test_compute_band_envelope:fail', '%d test(s) FAILED', nfail);
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

% Verbatim original inline theta-envelope block.
function [filtered, env, env_smooth] = ref_inline(signal, fs, band, filter_order, smooth_window)
[b, a] = butter(filter_order, band / (fs / 2), 'bandpass');
filtered = filtfilt(b, a, signal);
env = abs(hilbert(filtered));
env_smooth = fastsmooth(env, smooth_window, 1, 1);
end
