function band_power = zscore_smooth_bands(band_power, smooth_window)
%ZSCORE_SMOOTH_BANDS  Per-band z-score then moving-average smooth.
%   band_power = ZSCORE_SMOOTH_BANDS(band_power, smooth_window) takes an
%   nBands x nTime band-power matrix (e.g. from COMPUTE_BAND_POWER) and, for
%   each band (row) independently, z-scores it over time and then smooths it
%   with fastsmooth(..., smooth_window, 1, 1) (flat moving average).
%
%   Inputs
%   ------
%   band_power   : nBands x nTime matrix.
%   smooth_window: fastsmooth window length (the scripts use 20).
%
%   Output
%   ------
%   band_power : same size, each row z-scored then smoothed.
%
%   Bands are processed independently, so this reproduces the original inline
%   loops exactly whether they ran over one matrix (fiber) or several at once
%   (HP/mPFC/ipsiHP in a shared loop). Extracted (behaviour-preserving) from
%   the fiber preprocessing scripts. Unit-tested in
%   core/tests/test_zscore_smooth_bands.m.

for b = 1:size(band_power, 1)
    band_power(b, :) = zscore(band_power(b, :));
    band_power(b, :) = fastsmooth(band_power(b, :), smooth_window, 1, 1);
end
end
