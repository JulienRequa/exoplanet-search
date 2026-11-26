"""Targeted caveat checks for three flagged TOIs (diagnostic; v2 unchanged). 2026-11-26.

    python scripts/caveat_checks.py OUTDIR

1. TOI 4694.01: the v2 fit sits at the grazing bound (b = 1.20, Rp/Rs = 0.36), which inflates the
   largest secondary a planet could make. Refit the same v2-cleaned light curve with b fixed on a
   grid (0.0 to 1.2; Rp/Rs, a/Rs, t0 free; limb darkening fixed at the v2 values) and report χ²
   per b. Recompute the planet secondary limit (albedo 1 + dayside thermal) for each b.
2. TOI 4855.01: ExoFOP comment "slight depth aperture correlation". Measure the transit depth in
   circular apertures of radius 1, 1.5, 2, 3 px centred on the target (TESS pixel files; local
   out-of-transit baseline per transit), and compare with how the depth would change if the
   signal came from the target and the extra light were just the other stars (TIC Tmag, Gaussian
   PSF σ = 0.75 px).
3. TOI 4273.01: the size (2.68 R_Jup) depends on R* (TIC 2.23 R☉) and on SPOC's crowding
   correction (CROWDSAP 0.30). Fetch Gaia DR3 radii (GSP-Phot and FLAME) and measure the depth in
   small apertures (less contaminated) as in 2. Recompute Rp for the radius choices.
"""
import os, sys, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "."); sys.path.insert(0, "scripts")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
import lightkurve as lk
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
import pipeline as pl
from report_cards import t14
from diffimage_check import download_tpfs

HUNT = "results/hunt_v2_2026-11-25.csv"
RJ_PER_RSUN = 9.731


def planet_limit(k, a, teff):
    hP, kB, c = 6.626e-34, 1.381e-23, 3e8
    lam = np.linspace(600e-9, 1000e-9, 200)
    B = lambda T: 1 / (np.exp(hP * c / (lam * kB * T)) - 1) / lam ** 5
    tday = teff * np.sqrt(1 / a) * (2 / 3) ** 0.25
    return ((k / a) ** 2 + k ** 2 * np.trapz(B(tday), lam) / np.trapz(B(teff), lam)) * 1e6


def row(toi):
    h = pd.read_csv(HUNT)
    return h[h.TOI.round(2) == round(toi, 2)].iloc[0]


def check_4694(out):
    h = row(4694.01); tic = int(h.star.split()[-1])
    D = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
    lcs = pl.download(f"TIC {tic}", 2)
    flat = pl.clean_all(lcs, period=h.period, t0=h.t0, duration=D)
    t = np.asarray(flat.time.value, float); f = np.asarray(flat.flux.value, float)
    ok = np.isfinite(f); t, f = t[ok], f[ok]
    w = np.abs(((t - h.t0) / h.period + 0.5) % 1 - 0.5) * h.period < 3 * D
    t, f = t[w], f[w]
    sig = np.std(f[np.abs(((t - h.t0) / h.period + 0.5) % 1 - 0.5) * h.period > D])
    res = []
    for b in [0.0, 0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]:
        model = lambda tt, t0, k, a: pl.transit_model(tt, t0, h.period, k, a, b, h.u1, h.u2)
        best = None
        for k0 in (0.08, 0.15, 0.3):
            try:
                p, _ = curve_fit(model, t, f, p0=[h.t0, k0, max(h["a/Rs"], 3)],
                                 bounds=([h.t0 - 0.05, 0.01, 1.5], [h.t0 + 0.05, 1.0, 100]), maxfev=20000)
                chi = np.sum(((f - model(t, *p)) / sig) ** 2)
                if best is None or chi < best[0]:
                    best = (chi, p)
            except Exception:
                pass
        if best:
            chi, p = best
            res.append(dict(TOI=4694.01, b=b, chi2=chi, k=p[1], a=p[2],
                            rp_jup_TIC=p[1] * 1.53 * RJ_PER_RSUN, planet_limit_ppm=planet_limit(p[1], p[2], 5885)))
    df = pd.DataFrame(res); df["dchi2"] = df.chi2 - df.chi2.min()
    df.to_csv(f"{out}/toi4694_b_grid.csv", index=False)
    print(df.round(3).to_string()); return df


def aperture_depths(toi, radii=(1.0, 1.5, 2.0, 3.0)):
    h = row(toi); tic = int(h.star.split()[-1])
    D = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
    from astroquery.mast import Catalogs
    out = []
    for tpf in download_tpfs(tic):
        hd = tpf.hdu[0].header
        c0 = SkyCoord(hd["RA_OBJ"], hd["DEC_OBJ"], unit="deg")
        tx, ty = [float(v) for v in tpf.wcs.world_to_pixel(c0)]
        tt = Catalogs.query_region(c0, radius=150 * u.arcsec, catalog="TIC")
        stars = []
        for r in tt:
            cx, cy = [float(v) for v in tpf.wcs.world_to_pixel(SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg"))]
            stars.append((int(r["ID"]), cx, cy, 10 ** (-0.4 * float(r["Tmag"]))))
        t = np.asarray(tpf.time.value, float); cube = np.asarray(tpf.flux.value, float)
        good = (tpf.quality == 0) & np.isfinite(t); t, cube = t[good], cube[good]
        ny, nx = cube.shape[1:]; yy, xx = np.mgrid[0:ny, 0:nx]
        n = np.round((t - h.t0) / h.period); dt = t - (h.t0 + n * h.period)
        for rad in radii:
            ap = np.hypot(xx - tx, yy - ty) <= rad
            lc = np.nansum(cube[:, ap], axis=1)
            ds = []
            for k in np.unique(n[np.abs(dt) < 0.35 * D]):
                m = n == k
                inn, oo = m & (np.abs(dt) < 0.35 * D), m & (np.abs(dt) > 0.75 * D) & (np.abs(dt) < 2 * D)
                if inn.sum() > 5 and oo.sum() > 10:
                    ds.append(1 - np.nanmean(lc[inn]) / np.nanmean(lc[oo]))
            # expected fraction of aperture light from the target (Gaussian PSF, σ 0.75 px)
            fr = []
            for sid, cx, cy, fl in stars:
                psf = np.exp(-((xx - cx) ** 2 + (yy - cy) ** 2) / (2 * 0.75 ** 2)); psf /= psf.sum()
                fr.append((sid, fl * psf[ap].sum()))
            tot = sum(v for _, v in fr); targ = sum(v for s, v in fr if s == tic)
            out.append(dict(TOI=toi, sector=tpf.sector, radius_px=rad, n_transits=len(ds),
                            depth_ppm=np.mean(ds) * 1e6 if ds else np.nan,
                            depth_err_ppm=np.std(ds) / np.sqrt(len(ds)) * 1e6 if len(ds) > 1 else np.nan,
                            target_light_fraction=targ / tot if tot else np.nan))
    df = pd.DataFrame(out)
    # if the signal is on the target, depth × (1/target fraction) should be ~constant
    df["undiluted_if_on_target_ppm"] = df.depth_ppm / df.target_light_fraction
    return df


def gaia_radii(toi):
    from astroquery.vizier import Vizier
    h = row(toi); tic = int(h.star.split()[-1])
    from astroquery.mast import Catalogs
    r = Catalogs.query_criteria(catalog="TIC", ID=tic)[0]
    c = SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")
    v = Vizier(columns=["**", "+_r"], row_limit=5, vizier_server="vizier.cfa.harvard.edu")
    out = dict(TIC_rad=float(r["rad"]) if r["rad"] else np.nan, GAIA=str(r["GAIA"]))
    try:
        p = v.query_region(c, radius=3 * u.arcsec, catalog="I/355/paramp")[0][0]
        for col in ("Rad", "Rad-Flame", "Teff", "logg", "Lum-Flame"):
            if col in p.colnames:
                out[col] = float(p[col]) if p[col] is not None else np.nan
    except Exception as e:
        out["gaia_error"] = type(e).__name__
    return out


if __name__ == "__main__":
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    print("== 4694.01 b grid"); check_4694(out)
    print("== 4855.01 aperture depths")
    a = aperture_depths(4855.01); a.to_csv(f"{out}/toi4855_apertures.csv", index=False); print(a.round(3).to_string())
    print("== 4273.01 aperture depths + Gaia radii")
    b = aperture_depths(4273.01); b.to_csv(f"{out}/toi4273_apertures.csv", index=False); print(b.round(3).to_string())
    g = gaia_radii(4273.01); print(g)
    pd.DataFrame([g]).to_csv(f"{out}/toi4273_radii.csv", index=False)
