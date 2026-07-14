function plv = compute_plv(phase_difference)
%COMPUTE_PLV  Phase-locking value over time from a phase-difference matrix.
%   plv = COMPUTE_PLV(phase_difference) returns, for each row (frequency) of
%   phase_difference, the phase-locking value across the columns (time):
%
%       plv = |mean_t( exp(1i * phase_difference) )|
%
%   computed with nanmean so NaN samples are ignored. plv is in [0, 1]: 1 means
%   the phase difference is constant over time (perfectly locked), ~0 means it
%   is uniformly distributed (no locking).
%
%   Inputs
%   ------
%   phase_difference : nFreq x nTime matrix of phase differences (radians),
%                      e.g. circ_dist(phase_fiber, phase_lfp). Pass a column
%                      subset (phase_difference(:, idx)) to restrict to a set of
%                      time points (e.g. running vs rest).
%
%   Output
%   ------
%   plv : nFreq x 1 vector of phase-locking values.
%
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts, where
%   this one-liner was repeated for plv_all / plv_running / plv_rest / plv_lfps.
%   The phase-difference itself (e.g. via circ_dist) is computed by the caller.
%   Unit-tested in core/tests/test_compute_plv.m.

plv = abs(nanmean(exp(1i .* phase_difference), 2));
end
