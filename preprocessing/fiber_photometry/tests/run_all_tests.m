function run_all_tests()
%RUN_ALL_TESTS  Run every offline test for the fiber preprocessing refactor.
%   One entry point for the whole test suite. Runs:
%     * validate_refactor          - core/utils equivalence + path priority +
%                                    output-comparator self-test (no data)
%     * core/utils/tests/...       - shared helper equivalence tests
%     * core/tests/...             - extracted core-function unit tests
%   and prints a PASS/FAIL summary. Returns a non-zero MATLAB error if any
%   suite fails, so it can gate a commit.
%
%   None of these need raw data, hardware, or interactive input -- they are the
%   "Level 1" checks. The real-recording ("Level 2") check is separate; see
%   VALIDATION.md and tests/check_output_struct.m.
%
%   Run:  cd .../preprocessing/fiber_photometry; addpath('tests'); run_all_tests

here = fileparts(mfilename('fullpath'));               % .../fiber_photometry/tests
root = fileparts(here);                                % .../fiber_photometry
addpath(here);
addpath(genpath(fullfile(root, 'core')));

suites = {
    'validate_refactor',          @() run_quiet(@validate_refactor)
    'test_core_utils',            @() runfn(fullfile(root,'core','utils','tests'), @test_core_utils)
    'test_detect_trial_folders',  @() runfn(fullfile(root,'core','tests'), @test_detect_trial_folders)
    'test_correct_photobleaching',@() runfn(fullfile(root,'core','tests'), @test_correct_photobleaching)
    'test_detect_stim_onset',     @() runfn(fullfile(root,'core','tests'), @test_detect_stim_onset)
    'test_compute_deltaF_F',      @() runfn(fullfile(root,'core','tests'), @test_compute_deltaF_F)
    'test_compute_band_power',    @() runfn(fullfile(root,'core','tests'), @test_compute_band_power)
    'test_zscore_smooth_bands',   @() runfn(fullfile(root,'core','tests'), @test_zscore_smooth_bands)
    'test_compute_plv',           @() runfn(fullfile(root,'core','tests'), @test_compute_plv)
    'test_compute_band_envelope', @() runfn(fullfile(root,'core','tests'), @test_compute_band_envelope)
    'test_remove_stim_artifacts', @() runfn(fullfile(root,'core','tests'), @test_remove_stim_artifacts)
    };

fprintf('\n========================================================\n');
fprintf('  FIBER PREPROCESSING -- FULL OFFLINE TEST SUITE\n');
fprintf('========================================================\n');

n = size(suites, 1);
results = cell(n, 1);
nfail = 0;
for k = 1:n
    name = suites{k, 1};
    try
        suites{k, 2}();
        results{k} = 'PASS';
    catch ME
        results{k} = sprintf('FAIL (%s)', ME.message);
        nfail = nfail + 1;
    end
end

fprintf('\n--------------------- SUMMARY --------------------------\n');
for k = 1:n
    fprintf('  %-30s %s\n', suites{k, 1}, results{k});
end
fprintf('--------------------------------------------------------\n');

if nfail == 0
    fprintf('  ALL %d SUITES PASSED\n', n);
    fprintf('========================================================\n');
else
    error('run_all_tests:fail', '%d of %d suites FAILED (see summary above)', nfail, n);
end
end

% =========================================================================
function runfn(test_dir, fn)
% Run a test function from its own directory (tests addpath their own deps).
old = cd(test_dir);
cleanup = onCleanup(@() cd(old));
fn();
end

function run_quiet(fn)
% validate_refactor lives in tests/ (already on path); just call it.
fn();
end
