function test_zscore_smooth_bands()
%TEST_ZSCORE_SMOOTH_BANDS  Unit test for core/zscore_smooth_bands.
%   Covers: per-band independence, output shape, and bit-identical
%   equivalence to the original inline loop (both the single-matrix fiber
%   form and the shared-loop HP/mPFC/ipsiHP form, which are equivalent
%   because each band row is processed independently).
%
%   Run:  test_zscore_smooth_bands

addpath(fileparts(fileparts(mfilename('fullpath'))));            % core/
addpath(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'utils'));  % fastsmooth
nfail = 0;
win = 20;

% --- equivalence with the inline single-matrix loop ---------------------
rng(3);
bp = randn(5, 300);
nfail = nfail + expect_eq('matches inline (fiber form)', ...
    zscore_smooth_bands(bp, win), ref_inline(bp, win));
nfail = nfail + expect_eq('shape preserved', size(zscore_smooth_bands(bp, win)), [5 300]);

% --- per-band independence: shared-loop form == 3 separate calls --------
% Original LFP block ran one loop updating 3 matrices row-by-row; that equals
% calling the helper on each matrix separately.
A = randn(5, 200); B = randn(5, 200); C = randn(5, 200);
[A2, B2, C2] = ref_inline_shared(A, B, C, win);
nfail = nfail + expect_eq('shared-loop A == helper A', zscore_smooth_bands(A, win), A2);
nfail = nfail + expect_eq('shared-loop B == helper B', zscore_smooth_bands(B, win), B2);
nfail = nfail + expect_eq('shared-loop C == helper C', zscore_smooth_bands(C, win), C2);

% --- a single-row matrix still works ------------------------------------
one = randn(1, 120);
nfail = nfail + expect_eq('single band row', zscore_smooth_bands(one, win), ref_inline(one, win));

if nfail == 0
    fprintf('\nALL zscore_smooth_bands TESTS PASSED\n');
else
    error('test_zscore_smooth_bands:fail', '%d test(s) FAILED', nfail);
end
end

% =========================================================================
function nf = expect_eq(name, got, want)
if isequal(got, want)
    nf = 0;
else
    fprintf('  FAIL %s\n', name);
    nf = 1;
end
end

% Verbatim original inline loop (single matrix).
function bp = ref_inline(bp, smooth_window)
for b = 1:size(bp, 1)
    bp(b, :) = zscore(bp(b, :));
    bp(b, :) = fastsmooth(bp(b, :), smooth_window, 1, 1);
end
end

% Verbatim original inline shared loop (three matrices in one b-loop).
function [bp1, bp2, bp3] = ref_inline_shared(bp1, bp2, bp3, smooth_window)
for b = 1:size(bp1, 1)
    bp1(b, :) = zscore(bp1(b, :));
    bp1(b, :) = fastsmooth(bp1(b, :), smooth_window, 1, 1);
    bp2(b, :) = zscore(bp2(b, :));
    bp2(b, :) = fastsmooth(bp2(b, :), smooth_window, 1, 1);
    bp3(b, :) = zscore(bp3(b, :));
    bp3(b, :) = fastsmooth(bp3(b, :), smooth_window, 1, 1);
end
end
