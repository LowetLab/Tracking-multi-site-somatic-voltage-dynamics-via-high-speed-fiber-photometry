function pulses = generate_biphasic_pulses(t, stim_onset_sec, stim_freq_hz, stim_duration_sec, ~)
%GENERATE_BIPHASIC_PULSES  Synthetic square biphasic stimulation pulse train.
%   Uses a 50% duty cycle for visualization (matches the Python implementation).
%
%   pulses = generate_biphasic_pulses(t, stim_onset_sec, stim_freq_hz, ...
%                                     stim_duration_sec, pulse_width_us)
%
%   Inputs
%   ------
%   t                 : time vector in seconds (absolute, from recording start)
%   stim_onset_sec    : stimulation onset time (s)
%   stim_freq_hz      : stimulation frequency (Hz), e.g. 20, 40, 135
%   stim_duration_sec : duration of stimulation (s)
%   pulse_width_us    : ignored -- kept for API compatibility
%
%   Output
%   ------
%   pulses : same size as t; +1 = positive phase, -1 = negative phase, 0 = none

pulses = zeros(size(t));

% Calculate stim end time
stim_end_sec = stim_onset_sec + stim_duration_sec;

% Only generate pulses during stimulation period
stim_mask = (t >= stim_onset_sec) & (t <= stim_end_sec);
if ~any(stim_mask)
    return;
end

% Period of one complete biphasic pulse (positive + negative phase)
period = 1.0 / stim_freq_hz;  % seconds
half_period = period / 2.0;   % Each phase is half the period (50% duty cycle)

% Time relative to stim onset
t_relative = t(stim_mask) - stim_onset_sec;

% For each time point, determine which phase of the pulse cycle we're in
cycle_times = mod(t_relative, period);

% Positive phase: 0 to half_period (50% of cycle)
% Negative phase: half_period to period (50% of cycle)
positive_phase = cycle_times < half_period;
negative_phase = ~positive_phase;  % Everything else is negative

% Set pulse values
stim_indices = find(stim_mask);
pulses(stim_indices(positive_phase)) = 1.0;
pulses(stim_indices(negative_phase)) = -1.0;

end
