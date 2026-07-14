function band_power = compute_band_power(spectrogram_complex, freqs, band_ranges)
%COMPUTE_BAND_POWER  Mean power in frequency bands from a complex spectrogram.
%   band_power = COMPUTE_BAND_POWER(spectrogram_complex, freqs, band_ranges)
%   returns, for each frequency band and each time bin, the mean power
%   (|S|.^2 averaged over the band's frequency bins) of a short-time Fourier
%   spectrogram S.
%
%   Inputs
%   ------
%   spectrogram_complex : nFreq x nTime complex spectrogram (e.g. the first
%                         output of MATLAB's spectrogram()).
%   freqs               : nFreq vector of frequencies (the spectrogram's w).
%   band_ranges         : nBands x 2 matrix of [low high] band edges (Hz),
%                         inclusive on both ends (e.g. BAND_RANGES).
%
%   Output
%   ------
%   band_power : nBands x nTime matrix; row b is the mean power in
%                band_ranges(b,:) over time. No z-scoring or smoothing is
%                applied here (the scripts do that separately).
%
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts,
%   where this loop was repeated for the fiber, HP-LFP, mPFC and ipsiHP
%   spectrograms. Unit-tested in core/tests/test_compute_band_power.m.

n_bands = size(band_ranges, 1);
band_power = zeros(n_bands, size(spectrogram_complex, 2));
for b = 1:n_bands
    f_idx = freqs >= band_ranges(b,1) & freqs <= band_ranges(b,2);
    band_power(b, :) = mean(abs(spectrogram_complex(f_idx, :)).^2, 1);
end
end
