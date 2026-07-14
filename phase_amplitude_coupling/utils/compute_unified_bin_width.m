function bin_width = compute_unified_bin_width(all_values, bin_width_sci, target_bins, min_bins)
%COMPUTE_UNIFIED_BIN_WIDTH  Pick a "nice" histogram bin width for surrogate/observed MI plots.
%   bin_width = COMPUTE_UNIFIED_BIN_WIDTH(all_values, bin_width_sci, target_bins, min_bins)
%   scales bin_width_sci (a bin width already expressed on a 0-10 mantissa)
%   to the order of magnitude of all_values, then widens it (snapping to a
%   {0.1, 0.2, 0.5, 1.0} mantissa) if that would otherwise produce fewer than
%   min_bins bins over the data range.
%
%   Shared by the three compute_mi_from_wavelet_*.m surrogate-method scripts,
%   which used to each carry a byte-identical copy of this function.

    all_values = all_values(isfinite(all_values));
    if isempty(all_values)
        bin_width = 1;
        return;
    end
    dmin = max(0, min(all_values));
    dmax = max(all_values);
    if dmax <= dmin
        dmax = dmin + 1;
    end
    data_range = dmax - dmin;
    scale = max(abs([dmin, dmax]));
    if scale <= 0
        bin_width = data_range / target_bins;
    else
        exponent = floor(log10(scale));
        bin_width = bin_width_sci * 10^exponent;
    end
    expected_nbins = ceil(data_range / bin_width);
    if expected_nbins < min_bins
        bin_width = data_range / min_bins;
        if bin_width > 0
            exp_new = floor(log10(bin_width));
            mantissa = bin_width / 10^exp_new;
            nice_vals = [0.1, 0.2, 0.5, 1.0];
            [~, idx] = min(abs(nice_vals - mantissa));
            bin_width = nice_vals(idx) * 10^exp_new;
        end
    end
end
