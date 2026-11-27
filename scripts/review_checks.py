"""Expert-review checks for the 4 strong candidates (diagnostic; v2 unchanged). 2026-11-26.

    python scripts/review_checks.py OUTDIR

1. Density / duration. The transit gives ρ* = 3π (a/Rs)³ / (G P²), but a/Rs and b are degenerate
   (v2's covariance errors are huge here). So: refit the v2-cleaned light curve (v2 sectors) with
   b fixed on a grid (0–1.2; t0, Rp/Rs, a/Rs free; limb darkening fixed at the v2 values) and
   record χ² and ρ*(b). Compare with the catalogue density (below) and report Δχ² at the b where
   the transit density equals the catalogue density. Also the expected T14 at the best-fit b if the
   star had the catalogue density.
   Catalogue densities (ρ = 1.41 M/R³ g/cm³):
     3200.01: Gaia DR3 GSP-Phot R 1.43 (1.34–1.67), FLAME R 1.81 / M 1.06 → 0.25–0.51;
              asteroseismic Δν = 99 μHz (Karim et al. 2026) → ρ = 1.41 (Δν/135.1)² = 0.76.
              Range used: 0.25–0.76.
     3265.01: no Gaia/TIC parameters (no Gaia parallax). TIC "DWARF", Teff ≈ 6,125 K (ExoFOP) →
              F-dwarf range R 1.1–1.4, M 1.05–1.25 → 0.45–0.95.
     4842.01: TIC R 3.01 (giant); Gaia parallax 1.454 mas, G 12.49 → L ≈ 4.5 L☉ without
              extinction → R ≈ 2.4 R☉ at Teff 5,466 K; mass assumed 1.0–1.3 M☉ → 0.05–0.13.
     4694.01: Gaia DR3 R 1.40 (GSP-Phot) / 1.55 (FLAME), M 1.08 → 0.41–0.55.
2. Ellipsoidal variation. Per sector (SPOC 2-min PDCSAP, normalised): mask transit and secondary
   windows (±1 T14), fit 1 + quadratic trend + cos/sin(2πφ) + cos/sin(4πφ). Ellipsoidal light has
   minima at phase 0 and 0.5 → negative cos(4πφ) coefficient, sin(4πφ) ≈ 0. Errors: white-noise
   formal errors × √(points per 30 min) (a simple correlated-noise inflation). Expected amplitude
   A = α q (R1/a)³ sin²i with α ≈ 1.4 (u ≈ 0.6, τ ≈ 0.4), a/R1 from the catalogue density:
   planet q ≈ 0.001 (1 M_Jup) vs star q ≈ 0.1–1.
3. Is the secondary physically consistent with an eclipsing companion of R2 = k·R1 on the target?
   Secondary depth ≈ k² × S2/S1 (S = TESS-band surface brightness). Solve for the companion Teff
   (Planck, 600–1000 nm) and compare with the main-sequence radius at that Teff
   (≈ 0.15 R☉ at 3,000 K, 0.35 at 3,400, 0.6 at 3,900, 0.7 at 4,400, 0.8 at 5,000, 0.9 at 5,500,
   1.1 at 6,000 K). A required Teff whose MS radius is much larger than k·R1 means the pair cannot
   be the target plus a small companion: the eclipse is diluted, i.e. a blended EB.
"""
import os, sys, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "."); sys.path.insert(0, "scripts")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, brentq
import pipeline as pl
from report_cards import t14

G = 6.674e-8
RHO = {3200.01: (0.25, 0.76), 3265.01: (0.45, 0.95), 4842.01: (0.05, 0.13), 4694.01: (0.41, 0.55)}
TEFF = {3200.01: 5818, 3265.01: 6125, 4842.01: 5466, 4694.01: 6165}
R1 = {3200.01: (1.43, 1.81), 3265.01: (1.1, 1.4), 4842.01: (2.4, 3.0), 4694.01: (1.40, 1.55)}
SEC = {3265.01: (1666, 209), 3200.01: (824, 155), 4842.01: (690, 135), 4694.01: (2447, 576)}  # all-sector combined
MS = [(3000, .15), (3400, .35), (3900, .6), (4400, .7), (5000, .8), (5500, .9), (6000, 1.1), (6500, 1.3)]


def rho_of(a, P):
    return 3 * np.pi * a ** 3 / (G * (P * 86400) ** 2)


def a_of(rho, P):
    return (rho * G * (P * 86400) ** 2 / (3 * np.pi)) ** (1 / 3)


def sb_ratio(t2, t1):
    hP, kB, c = 6.626e-34, 1.381e-23, 3e8
    lam = np.linspace(600e-9, 1000e-9, 200)
    B = lambda T: 1 / (np.exp(hP * c / (lam * kB * T)) - 1) / lam ** 5
    return np.trapz(B(t2), lam) / np.trapz(B(t1), lam)


def load(toi):
    h = pd.read_csv("results/hunt_v2_2026-11-25.csv")
    r = h[h.TOI.round(2) == round(toi, 2)].iloc[0]
    tic = int(r.star.split()[-1])
    D = t14(r.period, r["Rp/Rs"], r["a/Rs"], r.b, r.inc)
    lcs = pl.download(f"TIC {tic}", 2)
    return r, D, lcs


def density_grid(toi, r, D, lcs):
    flat = pl.clean_all(lcs, period=r.period, t0=r.t0, duration=D)
    t = np.asarray(flat.time.value, float); f = np.asarray(flat.flux.value, float)
    ok = np.isfinite(f); t, f = t[ok], f[ok]
    ph = np.abs(((t - r.t0) / r.period + 0.5) % 1 - 0.5) * r.period
    w = ph < 3 * D; t, f = t[w], f[w]
    sig = np.std(f[ph[w] > D])
    rows = []
    for b in [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2]:
        model = lambda tt, t0, k, a: pl.transit_model(tt, t0, r.period, k, a, b, r.u1, r.u2)
        best = None
        for k0 in (0.05, 0.1, 0.2, 0.35):
            for a0 in (3, 8, 20):
                try:
                    p, _ = curve_fit(model, t, f, p0=[r.t0, k0, a0],
                                     bounds=([r.t0 - 0.05, 0.005, 1.2], [r.t0 + 0.05, 1.0, 150]), maxfev=20000)
                    chi = np.sum(((f - model(t, *p)) / sig) ** 2)
                    if best is None or chi < best[0]:
                        best = (chi, p)
                except Exception:
                    pass
        if best:
            chi, p = best
            rows.append(dict(TOI=toi, b=b, chi2=chi, k=p[1], a=p[2], rho=rho_of(p[2], r.period)))
    g = pd.DataFrame(rows); g["dchi2"] = g.chi2 - g.chi2.min()
    lo, hi = RHO[toi]
    best = g.loc[g.chi2.idxmin()]
    inside = g[(g.rho >= lo) & (g.rho <= hi)]
    # interpolate Δχ² where ρ(b) crosses the catalogue range
    cross = inside.dchi2.min() if len(inside) else np.nan
    if not len(inside):
        gg = g.sort_values("b")
        for i in range(len(gg) - 1):
            r1, r2 = gg.rho.iloc[i], gg.rho.iloc[i + 1]
            for target in (lo, hi):
                if (r1 - target) * (r2 - target) < 0:
                    x = (target - r1) / (r2 - r1)
                    cross = np.nanmin([cross, gg.dchi2.iloc[i] + x * (gg.dchi2.iloc[i + 1] - gg.dchi2.iloc[i])])
    a_cat = a_of(np.sqrt(lo * hi), r.period)
    s = np.sqrt(max((1 + best.k) ** 2 - best.b ** 2, 1e-6)) / a_cat
    t14_cat = r.period / np.pi * np.arcsin(min(s, 1)) * 24
    summ = dict(TOI=toi, best_b=best.b, best_rho=best.rho, rho_cat_lo=lo, rho_cat_hi=hi,
                ratio_best=best.rho / np.sqrt(lo * hi), dchi2_at_catalog_rho=cross,
                b_where_catalog=",".join(f"{x:.2f}" for x in inside.b) or "none on grid",
                t14_obs_h=D * 24, t14_expected_at_best_b_h=t14_cat)
    return g, summ


def ellipsoidal(toi, r, D, lcs):
    out = []
    for lc in lcs:
        lc = lc.remove_nans().normalize()
        t = np.asarray(lc.time.value, float); f = np.asarray(lc.flux.value, float)
        phi = ((t - r.t0) / r.period) % 1
        dist0 = np.minimum(phi, 1 - phi) * r.period; dist5 = np.abs(phi - 0.5) * r.period
        m = (dist0 > D) & (dist5 > D) & np.isfinite(f)
        t, f, phi = t[m], f[m], phi[m]
        x = (t - t.mean()) / (t.max() - t.min())
        A = np.vstack([np.ones_like(t), x, x ** 2, np.cos(2 * np.pi * phi), np.sin(2 * np.pi * phi),
                       np.cos(4 * np.pi * phi), np.sin(4 * np.pi * phi)]).T
        c, *_ = np.linalg.lstsq(A, f, rcond=None)
        res = f - A @ c
        cov = np.linalg.inv(A.T @ A) * np.var(res)
        infl = np.sqrt(30 / (np.median(np.diff(t)) * 1440))
        e = np.sqrt(np.diag(cov)) * infl
        out.append(dict(TOI=toi, sector=lc.meta.get("SECTOR"), cos4=c[5] * 1e6, cos4_err=e[5] * 1e6,
                        sin4=c[6] * 1e6, sin4_err=e[6] * 1e6, cos2=c[3] * 1e6, sin2=c[4] * 1e6))
    df = pd.DataFrame(out)
    lo, hi = RHO[toi]
    aR = a_of(np.sqrt(lo * hi), r.period)
    exp_planet = 1.4 * 1e-3 * aR ** -3 * 1e6
    exp_star = (1.4 * 0.1 * aR ** -3 * 1e6, 1.4 * 1.0 * aR ** -3 * 1e6)
    w = 1 / df.cos4_err ** 2
    mean = (df.cos4 * w).sum() / w.sum(); err = w.sum() ** -0.5
    chi = np.sum(((df.cos4 - mean) / df.cos4_err) ** 2)
    return df, dict(TOI=toi, cos4_mean=mean, cos4_err=err, cos4_chi2=chi, n_sectors=len(df),
                    a_over_R_catalog=aR, expected_planet_ppm=exp_planet,
                    expected_star_ppm=f"{exp_star[0]:.0f}-{exp_star[1]:.0f}")


def secondary_physics(toi, k):
    sec, e = SEC[toi]
    need = sec * 1e-6 / k ** 2
    T1 = TEFF[toi]
    try:
        T2 = brentq(lambda T: sb_ratio(T, T1) - need, 1500, 15000)
    except ValueError:
        T2 = np.nan
    ms_r = np.interp(T2, [m[0] for m in MS], [m[1] for m in MS]) if np.isfinite(T2) else np.nan
    r2 = [k * x for x in R1[toi]]
    return dict(TOI=toi, k=k, sec_ppm=sec, required_SB_ratio=need, required_T2=T2,
                MS_radius_at_T2=ms_r, companion_radius_kR1=f"{r2[0]:.2f}-{r2[1]:.2f}",
                expected_sec_if_MS_dwarf_of_that_size_ppm=k ** 2 * sb_ratio(
                    np.interp(np.mean(r2), [m[1] for m in MS], [m[0] for m in MS]), T1) * 1e6)


if __name__ == "__main__":
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    S1, S2, S3, G1, E1 = [], [], [], [], []
    for toi in [3265.01, 3200.01, 4842.01, 4694.01]:
        r, D, lcs = load(toi)
        g, s = density_grid(toi, r, D, lcs); G1.append(g); S1.append(s)
        e, s2 = ellipsoidal(toi, r, D, lcs); E1.append(e); S2.append(s2)
        kbest = g.loc[g.chi2.idxmin()].k
        S3.append(secondary_physics(toi, kbest))
        print(f"done {toi}", flush=True)
    pd.concat(G1).to_csv(f"{out}/density_bgrid.csv", index=False)
    pd.concat(E1).to_csv(f"{out}/ellipsoidal_per_sector.csv", index=False)
    for name, S in (("density", S1), ("ellipsoidal", S2), ("secondary_physics", S3)):
        d = pd.DataFrame(S); d.to_csv(f"{out}/{name}_summary.csv", index=False)
        print(f"== {name}"); print(d.round(3).to_string())
