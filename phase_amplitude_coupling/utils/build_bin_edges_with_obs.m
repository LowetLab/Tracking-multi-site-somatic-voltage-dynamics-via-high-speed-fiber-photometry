function edges = build_bin_edges_with_obs(dist, obs, bin_width)
%BUILD_BIN_EDGES_WITH_OBS  Histogram edges spanning a surrogate distribution and the observed value.
%   edges = BUILD_BIN_EDGES_WITH_OBS(dist, obs, bin_width) returns edges from
%   0 up to (at least) the larger of max(dist) and obs, in steps of bin_width,
%   with one extra trailing bin so the rightmost value is never on the edge.
%
%   Shared by the three compute_mi_from_wavelet_*.m surrogate-method scripts,
%   which used to each carry a byte-identical copy of this function.

    all_vals = [dist(:); obs];
    all_vals = all_vals(isfinite(all_vals));
    dmax = max(all_vals);
    if dmax == 0
        dmax = bin_width;
    end
    start_edge = 0;
    end_edge = ceil(dmax / bin_width) * bin_width;
    end_edge = end_edge + bin_width;
    edges = start_edge:bin_width:end_edge;
end
