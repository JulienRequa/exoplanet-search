"""Difference-image centroid on SIMULATED pixel data (PROTOTYPE, not in pipeline).
Can a PSF fit to (out-of-transit − in-transit) locate an eclipsing neighbour < 1 px away?
Run from the repo root: python sims/diffimage_sim.py"""
import sys
import numpy as np
from scipy.optimize import curve_fit

NPIX, SIG = 11, 0.75                       # stamp size; Gaussian PSF σ (FWHM ≈ 1.8 px)
yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)

def psf(x0, y0, sig=SIG):
    g = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sig ** 2))
    return g / g.sum()

def fit_source(img, sig_fit, pix_err=None):
    f = lambda xy, A, x0, y0, c: (A * psf(x0, y0, sig_fit) + c).ravel()
    p0 = [img.sum(), *np.unravel_index(np.argmax(img), img.shape)[::-1], 0.0]
    p, cov = curve_fit(f, None, img.ravel(), p0=p0, maxfev=5000,
                       sigma=None if pix_err is None else pix_err.ravel(),
                       absolute_sigma=pix_err is not None)
    return (p[1], p[2]) if pix_err is None else (p[1], p[2], np.sqrt(cov[1, 1]), np.sqrt(cov[2, 2]))

def simulate(rng, sep, fratio, on_target=False, sig_fit=SIG, tmag=11.0,
             depth_obs=0.005, n_tr=6, cad_per_tr=90, bkg=200.0, method="per-transit"):
    Ft = 15000 * 10 ** (-0.4 * (tmag - 10)) * 120            # e- per 2-min cadence
    Fn = fratio * Ft
    tgt = np.array([5.0, 5.0]) + rng.uniform(-0.5, 0.5, 2)
    ang = rng.uniform(0, 2 * np.pi)
    nb = tgt + sep * np.array([np.cos(ang), np.sin(ang)])
    # eclipse depth of the eclipsing star needed for the observed (diluted) depth
    if on_target:
        d_t, d_n = depth_obs * (Ft + Fn) / Ft, 0.0
    else:
        d_t, d_n = 0.0, depth_obs * (Ft + Fn) / Fn
        if d_n > 1:
            return None                                     # neighbour too faint to make the dip
    base_t, base_n = Ft * psf(*tgt), Fn * psf(*nb)
    def frames(n, dt, dn):
        mean = base_t * (1 - dt) + base_n * (1 - dn) + bkg
        return mean + rng.normal(0, 1, (n, NPIX, NPIX)) * np.sqrt(mean)
    if method == "stacked":
        # B: one stacked difference image; per-pixel noise measured from out-of-transit frames
        out_fr = frames(2 * cad_per_tr * n_tr, 0, 0)
        in_fr = frames(cad_per_tr * n_tr, d_t, d_n)
        d_img = out_fr.mean(0) - in_fr.mean(0)
        pix_sd = out_fr.std(0, ddof=1)                      # empirical (would include jitter)
        pix_err = pix_sd * np.sqrt(1 / len(out_fr) + 1 / len(in_fr))
        x, y, ex, ey = fit_source(d_img, sig_fit, pix_err)
        off = np.array([x, y]) - tgt
        from scipy.stats import chi2, norm
        z = norm.isf(chi2.sf((off[0] / ex) ** 2 + (off[1] / ey) ** 2, 2) / 2)
        true_off = (nb - tgt) if not on_target else np.zeros(2)
        return np.hypot(*off), np.hypot(*true_off), z
    pos = []
    for _ in range(n_tr):                                   # A: per-transit difference images
        d_img = frames(2 * cad_per_tr, 0, 0).mean(0) - frames(cad_per_tr, d_t, d_n).mean(0)
        pos.append(fit_source(d_img, sig_fit))
    pos = np.array(pos)
    src = pos.mean(0)
    off = src - tgt
    # Hotelling's T² (small-sample test of a 2-D mean offset): error estimated from only
    # n_tr per-transit positions, so the F distribution (not χ²) gives the right false-alarm rate
    from scipy.stats import f as fdist, norm
    n, p = len(pos), 2
    S = np.cov(pos.T, ddof=1)
    T2 = n * off @ np.linalg.solve(S, off)
    pval = fdist.sf((n - p) / (p * (n - 1)) * T2, p, n - p)
    z = norm.isf(pval / 2)
    true_off = (nb - tgt) if not on_target else np.zeros(2)
    return np.hypot(*off), np.hypot(*true_off), z

if __name__ == "__main__":
    rng = np.random.default_rng(1)
    N = 200
    print(f"{'case':38s} {'detect >3σ':>10s} {'measured sep (median)':>22s}")
    setups = [(SIG, "Tmag 11, depth 0.5%, 6 transits, PSF exact", {}),
              (1.15 * SIG, "same, PSF model 15% too wide", {}),
              (1.15 * SIG, "HARD: Tmag 13, depth 0.2%, 3 transits, PSF 15% wide",
               dict(tmag=13.0, depth_obs=0.002, n_tr=3)),
              (1.15 * SIG, "HARD, method B (stacked image, empirical pixel noise)",
               dict(tmag=13.0, depth_obs=0.002, n_tr=3, method="stacked"))]
    if len(sys.argv) > 1 and sys.argv[1] == "hard":
        setups = setups[2:]
    for sig_fit, label, kw in setups:
        print(f"--- {label}")
        res = [simulate(rng, 0.0, 0.2, on_target=True, sig_fit=sig_fit, **kw) for _ in range(N)]
        z = np.array([r[2] for r in res])
        print(f"{'NULL: eclipse on target (20% nbr)':38s} {np.mean(z > 3):10.0%} {np.median([r[0] for r in res]):17.3f} px")
        for fr in (0.2, 0.05):
            for sep in (0.3, 0.5, 0.8, 1.2):
                res = [r for r in (simulate(rng, sep, fr, sig_fit=sig_fit, **kw) for _ in range(N // 4)) if r]
                if not res:
                    print(f"neighbour {fr:.0%} flux at {sep} px".ljust(38), "  (too faint to make the dip)"); continue
                z = np.array([r[2] for r in res])
                print(f"{('neighbour %d%% flux at %.1f px' % (fr * 100, sep)):38s} {np.mean(z > 3):10.0%} "
                      f"{np.median([r[0] for r in res]):17.3f} px")
