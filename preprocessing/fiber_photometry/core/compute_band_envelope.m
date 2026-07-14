function [filtered, env, env_smooth] = compute_band_envelope(signal, fs, band, filter_order, smooth_window)
%COMPUTE_BAND_ENVELOPE  Band-pass a signal and take its smoothed Hilbert envelope.
%   [filtered, env, env_smooth] = COMPUTE_BAND_ENVELOPE(signal, fs, band,
%   filter_order, smooth_window) zero-phase band-pass filters `signal`
%   (Butterworth `filtfilt`), takes the analytic-signal amplitude envelope
%   (|hilbert|), and smooths the envelope with a flat moving average.
%
%   Inputs
%   ------
%   signal        : vector time series (e.g. aligned LFP).
%   fs            : sampling rate (Hz).
%   band          : [low high] pass band in Hz (e.g. [5 10] for theta).
%   filter_order  : Butterworth order (e.g. 3).
%   smooth_window : fastsmooth window for the envelope (e.g. 90).
%
%   Outputs
%   -------
%   filtered   : band-pass filtered signal.
%   env        : |hilbert(filtered)| amplitude envelope.
%   env_smooth : fastsmooth(env, smooth_window, 1, 1).
%
%   Requires the Signal Processing Toolbox (butter/filtfilt/hilbert).
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts'
%   theta-envelope-correlation block (band [5 10], order 3, smooth 90).
%   Unit-tested in core/tests/test_compute_band_envelope.m.

[b, a] = butter(filter_order, band / (fs / 2), 'bandpass');
filtered = filtfilt(b, a, signal);
env = abs(hilbert(filtered));
env_smooth = fastsmooth(env, smooth_window, 1, 1);
end
