function ok = validate_refactor()
%VALIDATE_REFACTOR  One-command sanity check for the preprocessing refactor.
%   Run this from MATLAB after pulling a new refactor stage. It needs NO data
%   and NO ROI selection -- it checks the things that can be checked offline:
%
%     1. The shared helpers resolve to core/utils (path priority is correct).
%     2. Those helpers are bit-for-bit equivalent to the originals
%        (runs core/utils/tests/test_core_utils.m).
%     3. The output comparator (compare_fiber_datastructs) itself works
%        (self-test on synthetic data).
%
%   Prints a clear PASS / FAIL and what to do next.
%
%   Usage:
%       cd preprocessing/fiber_photometry
%       addpath(genpath('core')); addpath('tests');
%       validate_refactor
%
%   For the FULL data-level check (after stages that change pipeline logic),
%   see VALIDATION.md -- you process one recording before & after and compare
%   the two output .mat files with compare_fiber_datastructs.

here = fileparts(mfilename('fullpath'));            % .../fiber_photometry/tests
root = fileparts(here);                             % .../fiber_photometry
addpath(genpath(fullfile(root, 'core')));
addpath(here);

fprintf('\n========================================================\n');
fprintf('  PREPROCESSING REFACTOR -- OFFLINE VALIDATION\n');
fprintf('========================================================\n');

pass = true;

% --- 1 & 2: helper equivalence + path priority ---------------------------
fprintf('\n[1/3] Shared helper equivalence + path priority...\n');
try
    test_core_utils();   % errors on any mismatch; prints its own summary
catch err
    pass = false;
    fprintf('  FAILED: %s\n', err.message);
end

% --- 3: comparator self-test --------------------------------------------
fprintf('\n[3/3] Output comparator self-test...\n');
try
    pass = comparator_selftest() && pass;
catch err
    pass = false;
    fprintf('  FAILED: %s\n', err.message);
end

% --- verdict -------------------------------------------------------------
fprintf('\n========================================================\n');
if pass
    fprintf('  OFFLINE VALIDATION: PASS\n');
    fprintf('  Next: do ONE real-recording check (see VALIDATION.md) so the\n');
    fprintf('  full pipeline is confirmed on your machine + data.\n');
else
    fprintf('  OFFLINE VALIDATION: FAIL  -- paste the output back, do not run\n');
    fprintf('  real data until this is green.\n');
end
fprintf('========================================================\n\n');

if nargout > 0, ok = pass; end
end

% =========================================================================
function ok = comparator_selftest()
% Confirm compare_fiber_datastructs flags differences and passes identicals.
s1.meta.mouse = 'Animal01';
s1.signals.lfp = (1:1000)' + 0.0;
s1.trials = {struct('x', randn(50,1)), struct('x', randn(10,3))};

% identical copy -> must PASS
r_same = compare_fiber_datastructs(s1, s1, 'MaxReport', 5);

% perturb one element -> must report exactly one mismatch
s2 = s1;
s2.signals.lfp(7) = s2.signals.lfp(7) + 1e-3;
r_diff = compare_fiber_datastructs(s1, s2, 'MaxReport', 5);

ok = r_same.pass && ~r_diff.pass && r_diff.n_mismatch == 1;
if ok
    fprintf('  comparator self-test PASSED (detects a 1e-3 change, ignores identicals)\n');
else
    fprintf('  comparator self-test FAILED (same.pass=%d diff.mismatch=%d)\n', ...
            r_same.pass, r_diff.n_mismatch);
end
end
