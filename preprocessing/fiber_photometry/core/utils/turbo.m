function cmap = turbo()
%TURBO  256x3 turbo colormap (interpolated from 11 anchor colours).
%   cmap = TURBO() returns a 256-by-3 RGB colormap approximating the turbo
%   colormap, built by pchip interpolation of 11 anchors.
%
%   This in-repo approximation is kept so figures reproduce exactly and so the
%   scripts do not depend on MATLAB's built-in turbo (R2020a+), which differs
%   slightly. Keep it on the TOP of the path for the fiber preprocessing scripts.
%
%   See also VIRIDIS, COLORMAP.

n = 256;
values = [
    0.18995, 0.07176, 0.23217;
    0.25107, 0.25237, 0.63374;
    0.27628, 0.48555, 0.85658;
    0.25862, 0.67862, 0.89715;
    0.32778, 0.84556, 0.79041;
    0.54658, 0.95717, 0.60574;
    0.76279, 0.97649, 0.42830;
    0.93717, 0.89854, 0.28334;
    0.98447, 0.71862, 0.14951;
    0.90006, 0.49541, 0.13068;
    0.70004, 0.24514, 0.10015
    ];

% Interpolate to create full colormap
xi = linspace(1, size(values, 1), n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(1:size(values, 1), values(:, i), xi, 'pchip');
end
end
