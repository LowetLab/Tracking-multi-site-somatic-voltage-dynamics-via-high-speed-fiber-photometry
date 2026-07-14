function smoothed = smooth2a(data, row_window, col_window)
%SMOOTH2A  Separable 2-D moving-average smoothing.
%   smoothed = SMOOTH2A(data, row_window, col_window) smooths a 2-D matrix by
%   a moving average of length row_window down the columns (dim 1) and
%   col_window across the rows (dim 2). If col_window is omitted it defaults to
%   row_window. A window of <= 1 leaves that dimension unchanged.
%
%   Convolution uses 'same' so the output is the same size as the input (edges
%   are implicitly zero-padded). NaNs propagate.
%
%   NOTE: this is the pipeline's own lightweight version. It deliberately keeps
%   identical behaviour to the copy that used to live inside the fiber
%   preprocessing scripts -- do not "upgrade" it (e.g. to a NaN-aware File
%   Exchange smooth2a) without re-validating downstream outputs.
%
%   NOTE: spectral_analysis/core/stim_spectral_analysis.m has its OWN local
%   function, smooth2a_radius_kernel, using a different (radius-based, single
%   2D box kernel) algorithm. The two are not interchangeable -- do not merge
%   them without re-validating downstream figure output.
%
%   See also CONV, FASTSMOOTH.

if nargin < 3
    col_window = row_window;
end

[rows, cols] = size(data);
smoothed = data;

% Smooth along rows
if row_window > 1
    kernel_row = ones(row_window, 1) / row_window;
    for c = 1:cols
        smoothed(:, c) = conv(data(:, c), kernel_row, 'same');
    end
end

% Smooth along columns
if col_window > 1
    kernel_col = ones(1, col_window) / col_window;
    for r = 1:rows
        smoothed(r, :) = conv(smoothed(r, :), kernel_col, 'same');
    end
end
end
