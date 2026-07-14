function test_detect_stim_onset()
%TEST_DETECT_STIM_ONSET  Unit test for core/detect_stim_onset.
%   Covers: basic rising-edge detection, no-onset -> [], the trial_start
%   preference, the fallback when no onset is after trial_start, and that the
%   two call forms reproduce the original two inline detection sites.
%
%   Run:  test_detect_stim_onset

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
nfail = 0;

% A step up at index 100 (so diff>thr at sample 100) and another at 500.
p = zeros(1000, 1);
p(101:end) = p(101:end) + 1;     % rising edge -> diff(100) = 1
p(501:end) = p(501:end) + 1;     % rising edge -> diff(500) = 1
thr = 0.1;

% --- basic: first onset, no trial filter --------------------------------
nfail = nfail + expect('basic first onset', detect_stim_onset(p, thr), 100);
nfail = nfail + expect('omitted trial_start == []', detect_stim_onset(p, thr, []), 100);

% --- no onset -> [] ------------------------------------------------------
nfail = nfail + expect('flat -> empty', detect_stim_onset(zeros(50,1), thr), []);
nfail = nfail + expect('below threshold -> empty', detect_stim_onset(p, 5), []);

% --- trial_start preference ---------------------------------------------
nfail = nfail + expect('prefer onset after trial_start', detect_stim_onset(p, thr, 200), 500);
% no onset after trial_start -> fall back to first overall
nfail = nfail + expect('fallback when none after trial_start', detect_stim_onset(p, thr, 900), 100);
% trial_start before all -> first overall
nfail = nfail + expect('trial_start before all', detect_stim_onset(p, thr, 10), 100);

% --- equivalence with the original inline logic --------------------------
nfail = nfail + expect('matches site-1 inline (with trial_start)', ...
    detect_stim_onset(p, thr, 200), ref_site1(p, thr, 200));
nfail = nfail + expect('matches site-2 inline (no trial_start)', ...
    detect_stim_onset(p, thr), ref_site2(p, thr));

if nfail == 0
    fprintf('\nALL detect_stim_onset TESTS PASSED\n');
else
    error('test_detect_stim_onset:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect(name, got, want)
if isequal(got, want)
    nf = 0;
else
    fprintf('  FAIL %s: got %s want %s\n', name, mat2str(got), mat2str(want));
    nf = 1;
end
end

% Verbatim original site-1 onset logic (with trial-period filtering).
function s = ref_site1(stim_pulses_for_period, STIM_ONSET_THRESHOLD, trial_start)
trial_starts = trial_start;
stim_onset_diff = diff(stim_pulses_for_period);
onset_candidates = find(stim_onset_diff > STIM_ONSET_THRESHOLD);
if ~isempty(onset_candidates)
    if exist('trial_starts', 'var') && ~isempty(trial_starts)
        valid_onsets = onset_candidates(onset_candidates > trial_starts(1));
        if ~isempty(valid_onsets), s = valid_onsets(1); else, s = onset_candidates(1); end
    else
        s = onset_candidates(1);
    end
else
    s = [];
end
end

% Verbatim original site-2 onset logic (no trial filtering).
function s = ref_site2(trial_stim_for_period, STIM_ONSET_THRESHOLD)
stim_onset_diff = diff(trial_stim_for_period);
onset_candidates = find(stim_onset_diff > STIM_ONSET_THRESHOLD);
if ~isempty(onset_candidates)
    s = onset_candidates(1);
else
    s = [];
end
end
