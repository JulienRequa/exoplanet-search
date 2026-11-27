"""v3 candidate: red-noise-aware detection SNR (time-averaging method; Pont, Zucker & Queloz 2006,
as implemented by Winn et al. 2008). NOT part of v2 (frozen); pipeline.py is unchanged.

Idea (same spirit as the v2 odd/even fix, PAPER_NOTES M2): the white-noise SNR assumes the scatter
of single points averages down as 1/sqrt(N). Measure how the OUT-OF-TRANSIT residuals actually
average down over transit-length bins, and inflate the noise by the ratio.

    beta_L = sigma_L(measured) / [ sigma_1 / sqrt(n_L) * sqrt(M / (M - 1)) ]

  sigma_1  = scatter of single out-of-transit points
  sigma_L  = scatter of the means of out-of-transit bins of length L (time-contiguous)
  n_L      = mean number of points per bin, M = number of bins
  L        = 0.75, 1.0, 1.25 × T14 (Winn et al. 2008 average over bin sizes near the duration)
  beta     = max(1, median of beta_L)      (never makes the SNR larger)

  SNR_red  = SNR_white / beta

Out of transit = more than 1 × T14 from the mid-transit (so the signal cannot leak into beta).
Bins with fewer than 80% of the expected points (gaps) are dropped. Robust scatter: 1.4826 × MAD.
"""
import numpy as np


def _rstd(x):
    x = np.asarray(x, float)
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def beta_time_averaging(t, f, t0, P, T14, factors=(0.75, 1.0, 1.25), min_fill=0.8):
    t = np.asarray(t, float); f = np.asarray(f, float)
    ok = np.isfinite(t) & np.isfinite(f)
    t, f = t[ok], f[ok]
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    oot = np.abs(ph) * P > T14
    t, r = t[oot], f[oot] - np.median(f[oot])
    if r.size < 100:
        return 1.0, {}
    dt = np.median(np.diff(t))
    s1 = _rstd(r)
    betas = {}
    for fac in factors:
        L = fac * T14
        idx = np.floor((t - t[0]) / L).astype(int)
        u, inv, cnt = np.unique(idx, return_inverse=True, return_counts=True)
        need = min_fill * L / dt
        good = cnt >= need
        if good.sum() < 10:
            continue
        means = np.bincount(inv, weights=r)[good] / cnt[good]
        M, n = good.sum(), cnt[good].mean()
        expected = s1 / np.sqrt(n) * np.sqrt(M / (M - 1))
        betas[fac] = _rstd(means) / expected
    if not betas:
        return 1.0, {}
    return max(1.0, float(np.median(list(betas.values())))), betas


def snr_white(t, f, model):
    """v2's detection SNR (pipeline.vet_star), reproduced for comparison."""
    intr = model < 1
    if not intr.any():
        return 0.0
    return (1 - model.min()) / (np.std(f[~intr]) / np.sqrt(intr.sum()))


def snr_rednoise(t, f, model, t0, P):
    """Returns (snr_white, beta, snr_red). T14 taken from the model itself."""
    t = np.asarray(t, float); f = np.asarray(f, float); model = np.asarray(model, float)
    sw = snr_white(t, f, model)
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    intr = model < 1
    T14 = 2 * np.max(np.abs(ph[intr])) * P if intr.any() else 0.0
    if T14 <= 0:
        return sw, 1.0, sw
    beta, _ = beta_time_averaging(t, f, t0, P, T14)
    return sw, beta, sw / beta
