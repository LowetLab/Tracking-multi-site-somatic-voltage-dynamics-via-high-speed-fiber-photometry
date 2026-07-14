function filtered_traces = remove_stim_artifacts(traces, fs)
%REMOVE_STIM_ARTIFACTS  Clean DBS-harmonic artifacts from fiber traces.
%   filtered_traces = REMOVE_STIM_ARTIFACTS(traces, fs) processes each column
%   (fiber) of `traces` independently:
%     1. replaces the final 81 samples with the median of the trace (excluding
%        its last 90 samples) to knock down an end-of-recording artifact, and
%     2. subtracts the 120-124 Hz and 130-132 Hz band-pass components
%        (zero-phase Butterworth, order 3) to remove stimulation harmonics.
%
%   Inputs
%   ------
%   traces : frames x fibers matrix of raw (ROI-extracted) fluorescence.
%   fs     : imaging sampling rate (Hz).
%
%   Output
%   ------
%   filtered_traces : same size as `traces`, artifact-removed.
%
%   The end-artifact replacement only runs when a trace has more than 90
%   samples. Requires the Signal Processing Toolbox (butter/filtfilt).
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts'
%   SECTION 4, where this loop was repeated verbatim. Unit-tested in
%   core/tests/test_remove_stim_artifacts.m.

filtered_traces = zeros(size(traces));

for trace_idx = 1:size(traces, 2)
    current_trace = traces(:, trace_idx);

    % Replace end artifacts
    if size(current_trace, 1) > 90
        current_trace(end-80:end) = median(current_trace(1:end-90));
    end

    % Remove 120-124 Hz artifact
    lowCutoff = 120;
    highCutoff = 124;
    filterOrder = 3;
    [b, a] = butter(filterOrder, [lowCutoff, highCutoff] / (fs / 2), 'bandpass');
    artifact_120Hz = filtfilt(b, a, current_trace);
    current_trace = current_trace - artifact_120Hz;

    % Remove 130-132 Hz artifact
    lowCutoff = 130;
    highCutoff = 132;
    [b, a] = butter(filterOrder, [lowCutoff, highCutoff] / (fs / 2), 'bandpass');
    artifact_130Hz = filtfilt(b, a, current_trace);
    current_trace = current_trace - artifact_130Hz;

    filtered_traces(:, trace_idx) = current_trace;
end
end
