function apply_simple_ticks(ax, edges, ~, unified_exp)
%APPLY_SIMPLE_TICKS  Set round-number x-tick marks (scaled by 10^unified_exp) on a histogram axis.
%   APPLY_SIMPLE_TICKS(ax, edges, ~, unified_exp) sets ax.XLim to span edges
%   and picks an integer tick step (in units of 10^unified_exp) so there are
%   roughly 5 ticks across the axis, labeling each with the unscaled integer.
%
%   Shared by the three compute_mi_from_wavelet_*.m surrogate-method scripts,
%   which used to each carry a byte-identical copy of this function.

    if isempty(edges)
        return;
    end
    xmin = edges(1);
    xmax = edges(end);
    ax.XLim = [xmin, xmax];

    scale = 10^unified_exp;
    max_scaled = xmax / scale;

    if max_scaled <= 5
        tick_step = 1;
    elseif max_scaled <= 10
        tick_step = 2;
    else
        tick_step = ceil(max_scaled / 5);
    end

    tick_vals_int = 0:tick_step:ceil(max_scaled);
    tick_vals = tick_vals_int * scale;
    tick_vals = tick_vals(tick_vals <= xmax * 1.01);

    ax.XTick = tick_vals;
    ax.XTickLabel = arrayfun(@num2str, tick_vals_int(1:length(tick_vals)), 'UniformOutput', false);
end
