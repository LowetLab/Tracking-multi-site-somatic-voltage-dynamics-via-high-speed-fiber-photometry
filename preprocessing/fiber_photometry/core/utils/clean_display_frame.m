function average_frame = clean_display_frame(image_stack, z_threshold)
%CLEAN_DISPLAY_FRAME  Time-averaged frame with bright outliers knocked down.
%   average_frame = CLEAN_DISPLAY_FRAME(image_stack, z_threshold) averages the
%   image stack over time (dim 3) and then replaces unusually bright pixels --
%   those whose z-score across the frame exceeds z_threshold -- so they do not
%   blow out the contrast when the frame is shown for ROI selection.
%
%   average_frame = CLEAN_DISPLAY_FRAME(image_stack) uses z_threshold = 15
%   (the value the fiber preprocessing scripts use for the ROI background).
%
%   This is a DISPLAY/diagnostic frame only -- it is not used for trace
%   extraction. Behaviour is preserved exactly from the original inline code,
%   including the two-step "set outliers to NaN, then fill NaNs with the frame
%   median" (one-sided, > z_threshold). Unit-tested in
%   core/utils/tests/test_core_utils.m.

if nargin < 2
    z_threshold = 15;
end

average_frame = mean(image_stack, 3);
average_frame(zscore(average_frame(:)) > z_threshold) = NaN;
average_frame(isnan(average_frame)) = median(average_frame(:));
end
