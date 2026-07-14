function test_core_utils()
%TEST_CORE_UTILS  Equivalence tests for the extracted core/utils helpers.
%   Confirms each extracted function reproduces, bit-for-bit, the original
%   implementation that used to live inside the fiber preprocessing scripts.
%   The reference implementations below are verbatim copies of the originals.
%
%   Run:  test_core_utils
%   Passes silently (prints a summary); errors on any mismatch.

rng(42);
nfail = 0;

% --- smooth2a -------------------------------------------------------------
M = randn(20, 15);
cases = {{3,5}, {1,4}, {4,1}, {1,1}, {6}};
for k = 1:numel(cases)
    a = cases{k};
    if numel(a) == 1
        got = smooth2a(M, a{1});         ref = ref_smooth2a(M, a{1});
    else
        got = smooth2a(M, a{1}, a{2});   ref = ref_smooth2a(M, a{1}, a{2});
    end
    nfail = nfail + check(sprintf('smooth2a case %d', k), got, ref);
end

% --- fastsmooth -----------------------------------------------------------
v = randn(100, 1);
for st = 1:3
    for em = 1:2
        got = fastsmooth(v, 7, st, em);  ref = ref_fastsmooth(v, 7, st, em);
        nfail = nfail + check(sprintf('fastsmooth type%d edge%d', st, em), got, ref);
    end
end
nfail = nfail + check('fastsmooth default args', fastsmooth(v, 5), ref_fastsmooth(v, 5));
nfail = nfail + check('fastsmooth window<=1', fastsmooth(v, 1), ref_fastsmooth(v, 1));

% --- viridis / turbo ------------------------------------------------------
nfail = nfail + check('viridis', viridis(), ref_viridis());
nfail = nfail + check('turbo',   turbo(),   ref_turbo());

% --- generate_biphasic_pulses --------------------------------------------
t = 0:0.001:5;
got = generate_biphasic_pulses(t, 1, 40, 2, 100);
ref = ref_genpulses(t, 1, 40, 2, 100);
nfail = nfail + check('generate_biphasic_pulses', got, ref);

% --- replace_outliers_with_median ----------------------------------------
% Equivalence to the inline LFP cleaning one-liner, on a trace with spikes.
lfp = randn(2000, 1);
lfp([37 500 1234]) = [80; -120; 95];   % extreme samples (|z| >> 10)
nfail = nfail + check('replace_outliers default thr=10', ...
    replace_outliers_with_median(lfp), ref_replace_outliers(lfp, 10));
nfail = nfail + check('replace_outliers explicit thr=10', ...
    replace_outliers_with_median(lfp, 10), ref_replace_outliers(lfp, 10));
% A clean trace (no |z|>thr) must be returned unchanged.
clean = randn(500, 1);
nfail = nfail + check('replace_outliers no-op on clean', ...
    replace_outliers_with_median(clean, 10), clean);
% Row vector keeps its orientation.
nfail = nfail + check('replace_outliers row vector', ...
    replace_outliers_with_median(lfp', 10), ref_replace_outliers(lfp', 10));

% --- clean_display_frame --------------------------------------------------
% Equivalence to the inline ROI-background cleaning. Inject a few very bright
% pixels into the stack so the time-averaged frame has clear outliers.
stack = 100 + 5 * randn(40, 50, 12);
stack(7, 9, :) = 9000;  stack(30, 41, :) = 12000;   % bright pixels -> z > 15
% NaN-aware compare: the median-fill leaves NaNs in place when outliers
% exist (preserved behaviour), and isequal treats NaN ~= NaN.
nfail = nfail + checkn('clean_display_frame default thr=15', ...
    clean_display_frame(stack), ref_clean_frame(stack, 15));
nfail = nfail + checkn('clean_display_frame explicit thr=8', ...
    clean_display_frame(stack, 8), ref_clean_frame(stack, 8));

% --- path-priority report -------------------------------------------------
fprintf('\nResolved on path:\n');
for f = {'smooth2a','fastsmooth','viridis','turbo','generate_biphasic_pulses','replace_outliers_with_median','clean_display_frame'}
    fprintf('   %-26s -> %s\n', f{1}, which(f{1}));
end

if nfail == 0
    fprintf('\nALL CORE/UTILS EQUIVALENCE TESTS PASSED\n');
else
    error('test_core_utils:fail', '%d core/utils test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = check(name, got, ref)
if isequal(got, ref)
    nf = 0;
else
    d = max(abs(double(got(:)) - double(ref(:))));
    fprintf('  FAIL %-28s max|diff|=%.3g\n', name, d);
    nf = 1;
end
end

% NaN-aware equality (NaN positions must match; NaN treated as equal).
function nf = checkn(name, got, ref)
if isequaln(got, ref)
    nf = 0;
else
    fprintf('  FAIL %-28s (isequaln mismatch)\n', name);
    nf = 1;
end
end

% =========================================================================
% Reference implementations -- VERBATIM copies of the original local functions.
% =========================================================================
function smoothed = ref_smooth2a(data, row_window, col_window)
if nargin < 3, col_window = row_window; end
[rows, cols] = size(data);
smoothed = data;
if row_window > 1
    kernel_row = ones(row_window, 1) / row_window;
    for c = 1:cols
        smoothed(:, c) = conv(data(:, c), kernel_row, 'same');
    end
end
if col_window > 1
    kernel_col = ones(1, col_window) / col_window;
    for r = 1:rows
        smoothed(r, :) = conv(smoothed(r, :), kernel_col, 'same');
    end
end
end

function smoothed = ref_fastsmooth(data, window_size, smooth_type, edge_mode)
if nargin < 3, smooth_type = 1; end
if nargin < 4, edge_mode = 1; end
if window_size <= 1, smoothed = data; return; end
data = data(:);
switch smooth_type
    case 1
        kernel = ones(window_size, 1) / window_size;
    case 2
        half_width = floor(window_size / 2);
        kernel = [1:half_width, half_width+1, half_width:-1:1]';
        kernel = kernel / sum(kernel);
    case 3
        sigma = window_size / 6;
        x = -floor(window_size/2):floor(window_size/2);
        kernel = exp(-(x.^2) / (2 * sigma^2));
        kernel = kernel' / sum(kernel);
    otherwise
        kernel = ones(window_size, 1) / window_size;
end
if edge_mode == 2
    pad_size = floor(length(kernel) / 2);
    padded_data = [repmat(data(1), pad_size, 1); data; repmat(data(end), pad_size, 1)];
    smoothed_padded = conv(padded_data, kernel, 'same');
    smoothed = smoothed_padded(pad_size+1:end-pad_size);
else
    smoothed = conv(data, kernel, 'same');
end
end

function cmap = ref_viridis()
n = 256;
values = [
    0.267004, 0.004874, 0.329415; 0.282623, 0.140926, 0.457517;
    0.253935, 0.265254, 0.529983; 0.206756, 0.371758, 0.553117;
    0.163625, 0.471133, 0.558148; 0.127568, 0.566949, 0.550556;
    0.134692, 0.658636, 0.517649; 0.266941, 0.748751, 0.440573;
    0.477504, 0.821444, 0.318195; 0.741388, 0.873449, 0.149561;
    0.993248, 0.906157, 0.143936 ];
xi = linspace(1, size(values, 1), n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(1:size(values, 1), values(:, i), xi, 'pchip');
end
end

function cmap = ref_turbo()
n = 256;
values = [
    0.18995, 0.07176, 0.23217; 0.25107, 0.25237, 0.63374;
    0.27628, 0.48555, 0.85658; 0.25862, 0.67862, 0.89715;
    0.32778, 0.84556, 0.79041; 0.54658, 0.95717, 0.60574;
    0.76279, 0.97649, 0.42830; 0.93717, 0.89854, 0.28334;
    0.98447, 0.71862, 0.14951; 0.90006, 0.49541, 0.13068;
    0.70004, 0.24514, 0.10015 ];
xi = linspace(1, size(values, 1), n);
cmap = zeros(n, 3);
for i = 1:3
    cmap(:, i) = interp1(1:size(values, 1), values(:, i), xi, 'pchip');
end
end

function y = ref_replace_outliers(x, z_threshold)
% Verbatim original inline behaviour: x(abs(zscore(x))>thr) = median(x).
y = x;
y(abs(zscore(x)) > z_threshold) = median(x);
end

function af = ref_clean_frame(image_stack, z_threshold)
% Verbatim original inline ROI-background cleaning.
af = mean(image_stack, 3);
af(zscore(af(:)) > z_threshold) = NaN;
af(isnan(af)) = median(af(:));
end

function pulses = ref_genpulses(t, stim_onset_sec, stim_freq_hz, stim_duration_sec, ~)
pulses = zeros(size(t));
stim_end_sec = stim_onset_sec + stim_duration_sec;
stim_mask = (t >= stim_onset_sec) & (t <= stim_end_sec);
if ~any(stim_mask), return; end
period = 1.0 / stim_freq_hz;
half_period = period / 2.0;
t_relative = t(stim_mask) - stim_onset_sec;
cycle_times = mod(t_relative, period);
positive_phase = cycle_times < half_period;
negative_phase = ~positive_phase;
stim_indices = find(stim_mask);
pulses(stim_indices(positive_phase)) = 1.0;
pulses(stim_indices(negative_phase)) = -1.0;
end
