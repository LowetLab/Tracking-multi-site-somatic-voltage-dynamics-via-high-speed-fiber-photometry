function test_compute_deltaF_F()
%TEST_COMPUTE_DELTAF_F  Unit test for core/compute_deltaF_F.
%   Covers: the ΔF/F definition (F_corr-F0)/F0 against a hand value, F0 taken
%   over the requested window only, multi-fiber independence, the z-scored
%   output, and bit-identical equivalence to the original inlined loop (the
%   block that was shared verbatim by the multi- and single-trial scripts).
%
%   Run:  test_compute_deltaF_F

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
nfail = 0;

% --- hand-checked single fiber, full-trace baseline ----------------------
% F_corr = [1 2 3 4], F0 = mean = 2.5, ΔF/F = (x-2.5)/2.5
x = [1; 2; 3; 4];
[dff, F0, z] = compute_deltaF_F(x, 1, 4);
nfail = nfail + expect_eq('F0 full window', F0, 2.5, 0);
nfail = nfail + expect_close('dF/F full window', dff, (x - 2.5) / 2.5, 1e-12);
nfail = nfail + expect_close('zscore full window', z, zscore(x), 1e-12);

% --- F0 from a sub-window only ------------------------------------------
% baseline = frames 1:2 -> F0 = 1.5; rest of trace must NOT affect F0.
[dff2, F02] = compute_deltaF_F(x, 1, 2);
nfail = nfail + expect_eq('F0 sub-window', F02, 1.5, 0);
nfail = nfail + expect_close('dF/F sub-window', dff2, (x - 1.5) / 1.5, 1e-12);

% --- multi-fiber independence + shapes ----------------------------------
M = [10 100; 20 200; 30 300; 40 400];   % 4 frames x 2 fibers
[dffM, F0M, zM] = compute_deltaF_F(M, 1, 4);
nfail = nfail + expect_eq('F0 per fiber', F0M, [25 250], 1e-12);
nfail = nfail + expect_eq('dF/F shape', size(dffM), size(M), 0);
nfail = nfail + expect_close('fiber1 indep of fiber2', dffM(:,1), (M(:,1)-25)/25, 1e-12);
nfail = nfail + expect_close('fiber2 indep of fiber1', dffM(:,2), (M(:,2)-250)/250, 1e-12);
nfail = nfail + expect_close('zscore per fiber', zM(:,2), zscore(M(:,2)), 1e-12);

% --- equivalence with the original inlined loop --------------------------
rng(7);
R = 100 + 5 * randn(500, 3);   % realistic-ish corrected fluorescence
bs = 50; be = 380;
[dffR, F0R, zR] = compute_deltaF_F(R, bs, be);
[dffRef, F0Ref, zRef] = ref_inline(R, bs, be);
nfail = nfail + expect_eq('matches inline dF/F (bit-identical)', dffR, dffRef, 0);
nfail = nfail + expect_eq('matches inline F0 (bit-identical)', F0R, F0Ref, 0);
nfail = nfail + expect_eq('matches inline zscore (bit-identical)', zR, zRef, 0);

if nfail == 0
    fprintf('\nALL compute_deltaF_F TESTS PASSED\n');
else
    error('test_compute_deltaF_F:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect_eq(name, got, want, ~)
if isequal(got, want)
    nf = 0;
else
    fprintf('  FAIL %s: got %s want %s\n', name, mat2str(got), mat2str(want));
    nf = 1;
end
end

function nf = expect_close(name, got, want, tol)
if isequal(size(got), size(want)) && all(abs(got(:) - want(:)) <= tol)
    nf = 0;
else
    fprintf('  FAIL %s: max abs diff %g (tol %g)\n', name, max(abs(got(:) - want(:))), tol);
    nf = 1;
end
end

% Verbatim original inlined ΔF/F loop (as it appeared in BOTH scripts).
function [traces_deltaF_F, F0_values, processed_traces_zscored] = ref_inline(traces_exp_corrected, baseline_start_frame, baseline_end_frame)
traces_deltaF_F = zeros(size(traces_exp_corrected));
F0_values = zeros(1, size(traces_exp_corrected, 2));
for fiber_idx = 1:size(traces_exp_corrected, 2)
    fiber_corrected = traces_exp_corrected(:, fiber_idx);
    baseline_window = fiber_corrected(baseline_start_frame:baseline_end_frame);
    F0 = mean(baseline_window);
    F0_values(fiber_idx) = F0;
    traces_deltaF_F(:, fiber_idx) = (fiber_corrected - F0) / F0;
end
processed_traces_zscored = zeros(size(traces_deltaF_F));
for i = 1:size(traces_deltaF_F, 2)
    processed_traces_zscored(:, i) = zscore(traces_deltaF_F(:, i));
end
end
