function mi = compute_tort_mi(phase, amp, nbin)
%COMPUTE_TORT_MI  Phase-amplitude coupling Modulation Index (Tort et al., 2010).
%   mi = COMPUTE_TORT_MI(phase, amp, nbin) computes the KL-divergence-based MI
%   between an instantaneous phase (radians, [-pi, pi]) and instantaneous
%   amplitude, binning phase into nbin equal-width bins.
%
%   Shared by the three compute_mi_from_wavelet_*.m surrogate-method scripts
%   (results / epoch_shuffle / iaaft_per_epoch), which used to each carry a
%   byte-identical copy of this function.

    phase = phase(:)';
    amp = amp(:)';

    % Create phase bins
    bin_edges = linspace(-pi, pi, nbin + 1);

    % Compute mean amplitude per bin
    mean_amp = zeros(1, nbin);
    for bi = 1:nbin
        if bi < nbin
            idx = (phase >= bin_edges(bi)) & (phase < bin_edges(bi+1));
        else
            idx = (phase >= bin_edges(bi)) & (phase <= bin_edges(bi+1));
        end
        if sum(idx) > 0
            mean_amp(bi) = mean(amp(idx));
        end
    end

    % Normalize to probability distribution
    if sum(mean_amp) > 0
        p = mean_amp / sum(mean_amp);
    else
        mi = 0;
        return;
    end

    % Uniform distribution
    q = ones(1, nbin) / nbin;

    % KL divergence (only sum over non-zero bins, since 0*log(0) -> 0)
    nonzero_idx = (p > 0);
    kl_div = sum(p(nonzero_idx) .* log(p(nonzero_idx) ./ q(nonzero_idx)));

    % Normalize by log(nbin) to get MI in [0, 1]
    mi = kl_div / log(nbin);
end
