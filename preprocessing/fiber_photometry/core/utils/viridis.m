function cmap = viridis()
%VIRIDIS  256x3 viridis colormap (interpolated from 11 anchor colours).
%   cmap = VIRIDIS() returns a 256-by-3 RGB colormap approximating the viridis
%   perceptually-uniform colormap, built by pchip interpolation of 11 anchors.
%
%   This in-repo approximation is kept so figures reproduce exactly regardless
%   of MATLAB version or whether color_maps/viridis.m is on the path. Keep it on
%   the TOP of the path for the fiber preprocessing scripts.
%
%   See also TURBO, COLORMAP.

n = 256;
values = [
    0.267004, 0.004874, 0.329415;
    0.282623, 0.140926, 0.457517;
    0.253935, 0.265254, 0.529983;
    0.206756, 0.371758, 0.553117;
    0.163625, 0.471133, 0.558148;
    0.127568, 0.566949, 0.550556;
    0.134692, 0.658636, 0.517649;
    0.266941, 0.748751, 0.440573;
    0.477504, 0.821444, 0.318195;
    0.741388, 0.873449, 0.149561;
    0.993248, 0.906157, 0.143936
    ];

% Interpolate to create full colormap
xi = linspace(1, size(values, 1), n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(1:size(values, 1), values(:, i), xi, 'pchip');
end
end
