function test_detect_trial_folders()
%TEST_DETECT_TRIAL_FOLDERS  Unit test for core/detect_trial_folders.
%   Builds a temporary folder tree and checks discovery, numeric sorting, and
%   the three error conditions. Prints a summary; errors on any failure.

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/ (detect_trial_folders)
nfail = 0;

% --- build a temp recording with mixed subfolders ------------------------
base = fullfile(tempdir, sprintf('dtf_test_%d', randi(1e9)));
cleanup = onCleanup(@() rmdir(base, 's'));
mkdir(base);
mk(base, 'recording_2');
mk(base, 'recording_10');   % must sort AFTER 2 (numeric, not lexicographic)
mk(base, 'recording_1');
mk(base, 'notes');          % no trailing _<number> -> ignored
mk(base, 'OpenEphys');      % ignored

% --- discovery + numeric sort --------------------------------------------
[folders, numbers, n] = detect_trial_folders(base);
nfail = nfail + expect('num_trials', n, 3);
nfail = nfail + expect('trial_numbers sorted', numbers, [1 2 10]);
nfail = nfail + expect('folders sorted', folders, {'recording_1','recording_2','recording_10'});

% --- error: no suffixed folders ------------------------------------------
base2 = fullfile(tempdir, sprintf('dtf_test2_%d', randi(1e9)));
c2 = onCleanup(@() rmdir(base2, 's'));
mkdir(base2); mk(base2, 'notes'); mk(base2, 'OpenEphys');
nfail = nfail + expect_error('no suffix -> error', @() detect_trial_folders(base2));

% --- error: missing base folder ------------------------------------------
nfail = nfail + expect_error('missing base -> error', ...
    @() detect_trial_folders(fullfile(tempdir, 'definitely_not_here_8417')));

if nfail == 0
    fprintf('\nALL detect_trial_folders TESTS PASSED\n');
else
    error('test_detect_trial_folders:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function mk(parent, name)
mkdir(fullfile(parent, name));
end

function nf = expect(name, got, want)
if isequal(got, want)
    nf = 0;
else
    fprintf('  FAIL %s\n', name);
    disp(got);
    nf = 1;
end
end

function nf = expect_error(name, fn)
threw = false;
try
    fn();
catch
    threw = true;
end
if threw
    nf = 0;
else
    fprintf('  FAIL %s (expected an error, none thrown)\n', name);
    nf = 1;
end
end
