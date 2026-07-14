function onset_sample = detect_stim_onset(stim_pulses, threshold, trial_start)
%DETECT_STIM_ONSET  First stimulation onset sample from a pulse/trigger trace.
%   onset_sample = DETECT_STIM_ONSET(stim_pulses, threshold) returns the index
%   of the first rising edge where diff(stim_pulses) > threshold, or [] if none.
%
%   onset_sample = DETECT_STIM_ONSET(stim_pulses, threshold, trial_start)
%   prefers the first onset that occurs AFTER trial_start (sample index); if no
%   onset occurs after trial_start it falls back to the first onset overall.
%   Pass trial_start = [] (or omit it) to disable this filtering.
%
%   Inputs
%   ------
%   stim_pulses : vector trigger/pulse trace (e.g. an Open Ephys ADC channel).
%   threshold   : rising-edge threshold on diff(stim_pulses) (e.g. 0.1).
%   trial_start : optional sample index; prefer the first onset after it.
%
%   Output
%   ------
%   onset_sample : index into diff(stim_pulses) of the detected onset, or [].
%
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts, where
%   it was inlined at two stimulation-detection sites. Unit-tested in
%   core/tests/test_detect_stim_onset.m.

if nargin < 3
    trial_start = [];
end

stim_onset_diff = diff(stim_pulses);
onset_candidates = find(stim_onset_diff > threshold);

if isempty(onset_candidates)
    onset_sample = [];
    return;
end

% Prefer the first onset within the trial period if a trial start is given.
if ~isempty(trial_start)
    valid_onsets = onset_candidates(onset_candidates > trial_start);
    if ~isempty(valid_onsets)
        onset_sample = valid_onsets(1);
    else
        onset_sample = onset_candidates(1);
    end
else
    onset_sample = onset_candidates(1);
end
end
