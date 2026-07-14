function test_compute_plv()
%TEST_COMPUTE_PLV  Unit test for core/compute_plv.
%   Covers: perfect locking -> 1, uniform spread -> ~0, NaN handling
%   (nanmean), per-row independence/shape, column-subset use, and
%   bit-identical equivalence to the original inline one-liner.
%
%   Run:  test_compute_plv

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
nfail = 0;

% --- perfect locking: constant phase diff per row -> plv = 1 -------------
pd = repmat([0; 0.7; -1.3], 1, 200);
nfail = nfail + expect_close('constant phase diff -> 1', compute_plv(pd), ones(3,1), 1e-12);

% --- uniform spread over a full circle -> ~0 ----------------------------
row = linspace(0, 2*pi, 100000);
row = row(1:end-1);                       % avoid duplicate 0/2pi endpoint
nfail = nfail + expect_close('uniform spread -> ~0', compute_plv(row), 0, 1e-3);

% --- NaN handling: nanmean ignores NaNs ---------------------------------
pdn = [0 0 0 NaN 0];                      % all defined samples identical
nfail = nfail + expect_close('NaN ignored', compute_plv(pdn), 1, 1e-12);

% --- shape + column subset ----------------------------------------------
rng(5); P = randn(8, 300);
nfail = nfail + expect_eq('output shape', size(compute_plv(P)), [8 1]);
idx = 1:50;
nfail = nfail + expect_close('column subset == manual', ...
    compute_plv(P(:, idx)), abs(nanmean(exp(1i .* P(:, idx)), 2)), 1e-12);

% --- equivalence with the original inline one-liner ---------------------
nfail = nfail + expect_eq('matches inline (bit-identical)', ...
    compute_plv(P), abs(nanmean(exp(1i .* P), 2)));

if nfail == 0
    fprintf('\nALL compute_plv TESTS PASSED\n');
else
    error('test_compute_plv:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect_eq(name, got, want)
if isequal(got, want), nf = 0; else
    fprintf('  FAIL %s\n', name); nf = 1;
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
