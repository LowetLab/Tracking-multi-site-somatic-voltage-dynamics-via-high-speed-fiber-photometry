function test_compute_band_power()
%TEST_COMPUTE_BAND_POWER  Unit test for core/compute_band_power.
%   Covers: output shape, a hand-checked single-band value, inclusive band
%   edges, and bit-identical equivalence to the original inlined loop (the
%   block repeated for the fiber / LFP / mPFC / ipsiHP spectrograms).
%
%   Run:  test_compute_band_power

addpath(fileparts(fileparts(mfilename('fullpath'))));   % core/
nfail = 0;

% --- hand-checked: constant-magnitude spectrogram -----------------------
% |S|=2 everywhere -> power = 4 in every band/time bin.
w = (1:10)';
S = 2 * ones(10, 7);                 % 10 freqs x 7 time bins
ranges = [1 3; 4 6; 7 10];
bp = compute_band_power(S, w, ranges);
nfail = nfail + expect_eq('shape', size(bp), [3 7]);
nfail = nfail + expect_close('constant power = 4', bp, 4 * ones(3, 7), 1e-12);

% --- band selection picks the right rows --------------------------------
% Put distinct magnitudes in distinct freq rows, one freq per band.
w2 = [1; 2; 3];
S2 = [ (1+0i)*ones(1,4); (3)*ones(1,4); (5)*ones(1,4) ];  % |S| = 1,3,5
r2 = [1 1; 2 2; 3 3];               % each band = exactly one freq row
bp2 = compute_band_power(S2, w2, r2);
nfail = nfail + expect_close('band1 power=1', bp2(1,:), ones(1,4), 1e-12);
nfail = nfail + expect_close('band2 power=9', bp2(2,:), 9*ones(1,4), 1e-12);
nfail = nfail + expect_close('band3 power=25', bp2(3,:), 25*ones(1,4), 1e-12);

% --- inclusive edges: a band covering 2-3 averages those two rows -------
bp3 = compute_band_power(S2, w2, [2 3]);
nfail = nfail + expect_close('inclusive edges mean(9,25)=17', bp3, 17*ones(1,4), 1e-12);

% --- equivalence with the original inlined loop -------------------------
rng(11);
Sr = randn(150, 40) + 1i*randn(150, 40);   % complex spectrogram
wr = (1:150)';
BAND_RANGES = [1 8; 8 12; 13 30; 31 70; 71 185];
got = compute_band_power(Sr, wr, BAND_RANGES);
ref = ref_inline(Sr, wr, BAND_RANGES, 5, size(Sr,2));
nfail = nfail + expect_eq('matches inline loop (bit-identical)', got, ref);

if nfail == 0
    fprintf('\nALL compute_band_power TESTS PASSED\n');
else
    error('test_compute_band_power:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect_eq(name, got, want)
if isequal(got, want)
    nf = 0;
else
    fprintf('  FAIL %s: got %s want %s\n', name, mat2str(got), mat2str(want));
    nf = 1;
end
end

function nf = expect_close(name, got, want, tol)
if isequal(size(got), size(want)) && all(abs(got(:) - want(:)) <= tol)
    nf = 0;
else
    fprintf('  FAIL %s: max abs diff %g (tol %g)\n', name, max(abs(got(:) - want(:))), tol);
    nf = 1;
end
end

% Verbatim original inlined band-power loop (as it appeared at every site).
% n_band_names mirrors length(BAND_NAMES); n_t mirrors length(t_*).
function band_power = ref_inline(s, w, BAND_RANGES, n_band_names, n_t)
band_power = zeros(n_band_names, n_t);
for b = 1:size(BAND_RANGES, 1)
    f_idx = w >= BAND_RANGES(b,1) & w <= BAND_RANGES(b,2);
    band_power(b, :) = mean(abs(s(f_idx, :)).^2, 1);
end
end
