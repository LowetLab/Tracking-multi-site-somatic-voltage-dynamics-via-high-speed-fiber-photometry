function stars = p_to_stars(p_value)
%P_TO_STARS  Convert a p-value to a conventional significance-star string.
%   stars = P_TO_STARS(p_value) returns '***' (p<0.001), '**' (p<0.01),
%   '*' (p<0.05), or 'n.s.' (not significant, or non-finite p_value).
%
%   Shared by the three compute_mi_from_wavelet_*.m surrogate-method scripts,
%   which used to each carry a byte-identical copy of this function.

    if ~isfinite(p_value)
        stars = 'n.s.';
    elseif p_value < 0.001
        stars = '***';
    elseif p_value < 0.01
        stars = '**';
    elseif p_value < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end
