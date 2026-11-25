"""Gaia neighbour-depth check on SIMULATED scenes (PROTOTYPE, not in pipeline).
Each hypothesis "star j is the one eclipsing" must (1) be able to make the observed depth
(d_j = δ_obs·S/S_j ≤ 1, S = aperture flux, S_j = star j's flux inside the aperture) and
(2) predict the measured in/out-of-transit moment-centroid shift Δ. Δ_pred is computed from a
PSF × aperture model of the catalogue stars, so aperture truncation is included.
Rule: FAIL if the target hypothesis is rejected (χ², 2 dof, p < 0.27%) AND some neighbour is
consistent; WEAK if the target is rejected and no neighbour fits; else PASS.
Criteria: planets on the target (incl. crowded fields like HATS-55 b, CROWDSAP ~0.4) are not
FAILed; eclipsing neighbours (incl. < 1 px) are caught.
Run from the repo root: python sims/neighbor_depth_sim.py"""
import numpy as np
from scipy.stats import chi2

NPIX, SIG_TRUE = 13, 0.75
yy, xx = np.mgrid[0:NPIX, 0:NPIX].astype(float)
CRIT = chi2.isf(0.0027, 2)                                       # 11.83

def psf(x0, y0, sig):
    g = np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * sig ** 2))
    return g / g.sum()

def counts(tmag):
    return 15000 * 10 ** (-0.4 * (tmag - 10)) * 120                # e- per 2-min cadence

def scene(rng, crowded, extra=None):
    """Target near the centre + neighbours within 2.5 px. Returns list of (x, y, tmag)."""
    tgt = (6.0 + rng.uniform(-0.5, 0.5), 6.0 + rng.uniform(-0.5, 0.5), 12.0)
    n = rng.poisson(5 if crowded else 1)
    stars = [tgt]
    for _ in range(n):
        r, a = 2.5 * np.sqrt(rng.uniform(0.04, 1)), rng.uniform(0, 2 * np.pi)
        stars.append((tgt[0] + r * np.cos(a), tgt[1] + r * np.sin(a),
                      12.0 + (rng.uniform(0, 5) if crowded else rng.uniform(1, 7))))
    if extra:
        stars.append(extra(tgt))
    return stars

def aperture(stars, sig):
    tp = psf(stars[0][0], stars[0][1], sig)
    return tp >= 0.02 * tp.max() * 2                                # SPOC-like, ~9 px

def model_terms(stars, ap, sig):
    """Per star: S_j (flux in aperture) and Sx_j, Sy_j (flux-weighted coordinates)."""
    out = []
    for x, y, m in stars:
        img = counts(m) * psf(x, y, sig) * ap
        out.append((img.sum(), (img * xx).sum(), (img * yy).sum()))
    return np.array(out)

def predict(terms, j, d_obs):
    S, Sx, Sy = terms.sum(0)
    Sj, Sxj, Syj = terms[j]
    dj = d_obs * S / Sj
    if dj > 1:
        return None                                                 # can't make the depth
    return np.array([(Sx - dj * Sxj) / (S - dj * Sj) - Sx / S, (Sy - dj * Syj) / (S - dj * Sj) - Sy / S])

def observe(rng, stars, ap, ecl, d_obs, n_in=540, n_out=1080, bkg=200.0, jitter=0.003):
    """Mean in/out images with photon noise and pointing jitter → δ_obs, Δ (measured)."""
    S_true = model_terms(stars, ap, SIG_TRUE)
    dj = d_obs * S_true.sum(0)[0] / S_true[ecl][0]
    def mean_img(n, dim, shift):
        img = sum(counts(m) * psf(x + shift[0], y + shift[1], SIG_TRUE) * (1 - (dim if i == ecl else 0))
                  for i, (x, y, m) in enumerate(stars)) + bkg
        return img + rng.normal(0, 1, img.shape) * np.sqrt(img / n)
    def cen(img):
        a = (img - bkg) * ap
        return a.sum(), np.array([(a * xx).sum(), (a * yy).sum()]) / a.sum()
    Fo, co = cen(mean_img(n_out, 0, rng.normal(0, jitter / np.sqrt(n_out / 30), 2)))
    Fi, ci = cen(mean_img(n_in, dj, rng.normal(0, jitter / np.sqrt(n_in / 30), 2)))
    return 1 - Fi / Fo, ci - co

def check(rng, stars, ecl, d_obs=0.005, sig_model=SIG_TRUE, wcs_err=0.3, mag_err=0.05):
    ap = aperture(stars, SIG_TRUE)
    d_meas, delta = observe(rng, stars, ap, ecl, d_obs)
    # null (no eclipse anywhere) → empirical σ of Δ, as the pipeline gets from fake transits
    null = np.array([observe(rng, stars, ap, 0, 0.0)[1] for _ in range(30)])
    sx, sy = null.std(0, ddof=1)
    # catalogue: common WCS offset + small per-star errors; model PSF possibly wrong
    off = wcs_err * np.array([np.cos(a := rng.uniform(0, 2 * np.pi)), np.sin(a)])
    cat = [(x + off[0] + rng.normal(0, 0.02), y + off[1] + rng.normal(0, 0.02), m + rng.normal(0, mag_err))
           for x, y, m in stars]
    terms = model_terms(cat, ap, sig_model)
    chis = []
    for j in range(len(cat)):
        pr = predict(terms, j, d_meas)
        if pr is None:
            chis.append(np.inf); continue
        sys2 = (0.1 * np.abs(pr)) ** 2 + 0.0005 ** 2                # model systematic floor
        chis.append(np.sum((delta - pr) ** 2 / (np.array([sx, sy]) ** 2 + sys2)))
    tgt_rej = chis[0] > CRIT
    nbr_ok = any(c <= CRIT for c in chis[1:])
    return "FAIL" if tgt_rej and nbr_ok else ("WEAK" if tgt_rej else "PASS")

if __name__ == "__main__":
    rng = np.random.default_rng(7)
    print("PLANETS on the target (should not FAIL):")
    for label, crowded, sigm in [("normal field", False, SIG_TRUE), ("crowded field (HATS-55-like)", True, SIG_TRUE),
                                 ("crowded field, PSF model 15% wide", True, 1.15 * SIG_TRUE)]:
        res = [check(rng, scene(rng, crowded), 0, sig_model=sigm) for _ in range(200)]
        print(f"  {label:36s} N=200: FAIL {res.count('FAIL')/2:.1f}%  WEAK {res.count('WEAK')/2:.1f}%")
    print("ECLIPSING NEIGHBOUR (should FAIL), normal field + one neighbour:")
    for dmag in (2.0, 4.0):
        for sep in (0.3, 0.5, 1.0, 2.0):
            res = []
            for _ in range(50):
                a = rng.uniform(0, 2 * np.pi)
                st = scene(rng, False, extra=lambda t: (t[0] + sep * np.cos(a), t[1] + sep * np.sin(a), 12.0 + dmag))
                res.append(check(rng, st, len(st) - 1, sig_model=1.15 * SIG_TRUE))
            print(f"  Δmag {dmag:.0f}, {sep:.1f} px: FAIL {res.count('FAIL')/50:4.0%}  WEAK {res.count('WEAK')/50:4.0%}"
                  f"  PASS {res.count('PASS')/50:4.0%}")
