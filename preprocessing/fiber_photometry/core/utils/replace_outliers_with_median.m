function y = replace_outliers_with_median(x, z_threshold)
%REPLACE_OUTLIERS_WITH_MEDIAN  Clamp extreme samples to the signal median.
%   y = REPLACE_OUTLIERS_WITH_MEDIAN(x, z_threshold) returns a copy of x in
%   which every sample whose absolute z-score exceeds z_threshold is replaced
%   by median(x). Both the z-scores and the median are taken from the ORIGINAL
%   x (so the replacement value does not depend on which samples are replaced).
%
%   y = REPLACE_OUTLIERS_WITH_MEDIAN(x) uses z_threshold = 10 (the value the
%   fiber preprocessing scripts use for cleaning aligned LFP traces).
%
%   Inputs
%   ------
%   x           : numeric vector (e.g. an aligned LFP trace).
%   z_threshold : absolute z-score above which a sample is an outlier (def 10).
%
%   Output
%   ------
%   y : x with outliers replaced by median(x); same size/orientation as x.
%
%   Extracted (behaviour-preserving) from the fiber preprocessing scripts,
%   where the one-liner `x(abs(zscore(x))>10)=median(x)` was repeated at every
%   LFP-alignment site (HP / mPFC / ipsiHP). Unit-tested in
%   core/utils/tests/test_core_utils.m.

if nargin < 2
    z_threshold = 10;
end

y = x;
y(abs(zscore(x)) > z_threshold) = median(x);
end
