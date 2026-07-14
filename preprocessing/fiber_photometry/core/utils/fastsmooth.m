function smoothed = fastsmooth(data, window_size, smooth_type, edge_mode)
%FASTSMOOTH  1-D smoothing by moving-average, triangular, or Gaussian kernel.
%   smoothed = FASTSMOOTH(data, window_size) applies a moving-average filter of
%   length window_size and returns a column vector.
%
%   smoothed = FASTSMOOTH(data, window_size, smooth_type, edge_mode)
%     smooth_type : 1 = moving average (default), 2 = triangular, 3 = Gaussian
%     edge_mode   : 1 = zero padding (default), 2 = replicate edges
%
%   A window_size of <= 1 returns the input unchanged.
%
%   IMPORTANT: this is the pipeline's OWN fastsmooth and is intentionally NOT
%   the same as O'Haver's widely-distributed fastsmooth.m (different kernels and
%   edge handling). Keep this version on the TOP of the MATLAB path for the
%   fiber preprocessing scripts so results stay reproducible.
%
%   See also CONV, SMOOTH2A.

if nargin < 3
    smooth_type = 1;
end
if nargin < 4
    edge_mode = 1;
end

if window_size <= 1
    smoothed = data;
    return;
end

data = data(:);  % Ensure column vector
n = length(data);

% Create kernel based on smooth_type
switch smooth_type
    case 1  % Moving average
        kernel = ones(window_size, 1) / window_size;
    case 2  % Triangular
        half_width = floor(window_size / 2);
        kernel = [1:half_width, half_width+1, half_width:-1:1]';
        kernel = kernel / sum(kernel);
    case 3  % Gaussian
        sigma = window_size / 6;
        x = -floor(window_size/2):floor(window_size/2);
        kernel = exp(-(x.^2) / (2 * sigma^2));
        kernel = kernel' / sum(kernel);
    otherwise
        kernel = ones(window_size, 1) / window_size;
end

% Apply convolution with edge handling
if edge_mode == 2
    % Replicate edges
    pad_size = floor(length(kernel) / 2);
    padded_data = [repmat(data(1), pad_size, 1); data; repmat(data(end), pad_size, 1)];
    smoothed_padded = conv(padded_data, kernel, 'same');
    smoothed = smoothed_padded(pad_size+1:end-pad_size);
else
    % Zero padding (default)
    smoothed = conv(data, kernel, 'same');
end
end
