"""Targeted difference-image check for report-card survivors (diagnostic; not a v2 change).

    python scripts/diffimage_check.py LIST.csv OUTDIR      (LIST has a TOI column)

Where does the dimming happen? For each sector of TESS 2-min SPOC pixel data (same first two
sectors as the hunt):
 - in-transit cadences: |t − tc| < 0.35·T14; local out-of-transit: 0.75·T14 < |t − tc| < 1.75·T14,
   per transit (each transit uses its own neighbouring baseline, so slow trends cancel).
 - difference image = mean(out) − mean(in), averaged over transits; pixel errors from the
   transit-to-transit scatter.
 - source position: Gaussian-PSF (σ 0.75 px) + constant fit to the difference image; the same fit
   to the mean out-of-transit image gives the light centroid of the field (mostly the target).
 - offsets: source − out-of-transit centroid (relative, cancels WCS errors) and source − target
   position (TIC coordinates moved for proper motion to the sector's mid-time).
 - error: scatter of per-transit source positions / √n (bootstrap-free; needs ≥ 3 transits).
A sector whose difference image has S/N < 5 at its peak pixel is "no detection".
Result per sector: "on target" if the source is < 3σ and < 0.5 px from the target position,
"off target" if > 3σ and > 0.5 px (nearest TIC star ≤ 2.5 px named), else "inconclusive".
Caveat: the real-data positional floor of this method is not yet measured (that is the v3
go/no-go test on known planets); treat results as indicative.
"""
import os, sys, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.time import Time
from scipy.optimize import curve_fit

HUNT = "results/hunt_v2_2026-11-25.csv"
SIG = 0.75


def t14(P, k, a, b, inc):
    s = np.sqrt(max((1 + k) ** 2 - b ** 2, 1e-6)) / (a * np.sin(np.radians(inc)))
    return P / np.pi * np.arcsin(min(s, 1.0))


def fit_psf(img, guess=None):
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    def f(_, A, x0, y0, c):
        return (A * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * SIG ** 2)) + c).ravel()
    if guess is None:
        y0, x0 = np.unravel_index(np.nanargmax(img), img.shape)
    else:
        x0, y0 = guess
    good = np.isfinite(img)
    im = np.where(good, img, np.nanmedian(img))
    p, _ = curve_fit(f, None, im.ravel(), p0=[np.nanmax(im) - np.nanmedian(im), x0, y0, np.nanmedian(im)],
                     bounds=([0, -0.5, -0.5, -np.inf], [np.inf, nx - 0.5, ny - 0.5, np.inf]), maxfev=5000)
    return p[1], p[2]


def download_tpfs(tic, n=2, attempts=3):
    for a in range(attempts):
        try:
            sr = lk.search_targetpixelfile(f"TIC {tic}", mission="TESS", author="SPOC", exptime=120)
            seen, keep = set(), []
            for i, m in enumerate(sr.table["mission"]):
                if str(m) not in seen:
                    seen.add(str(m)); keep.append(i)
            return list(sr[keep[:n]].download_all())
        except Exception:
            if a == attempts - 1:
                raise


def sector_result(tpf, tic, t0, P, D, nb):
    t = np.asarray(tpf.time.value, float)
    cube = np.asarray(tpf.flux.value, float)
    good = (tpf.quality == 0) & np.isfinite(t) & np.all(np.isfinite(cube.reshape(len(t), -1)) | True, axis=1)
    t, cube = t[good], cube[good]
    n = np.round((t - t0) / P)
    dt = t - (t0 + n * P)
    diffs, oots = [], []
    for k in np.unique(n[np.abs(dt) < 0.35 * D]):
        m = n == k
        inn = m & (np.abs(dt) < 0.35 * D)
        out = m & (np.abs(dt) > 0.75 * D) & (np.abs(dt) < 1.75 * D)
        if inn.sum() >= 5 and out.sum() >= 10:
            o = np.nanmean(cube[out], axis=0)
            diffs.append(o - np.nanmean(cube[inn], axis=0)); oots.append(o)
    if len(diffs) < 3:
        return None
    diffs, oots = np.array(diffs), np.array(oots)
    dimg, oimg = np.nanmean(diffs, axis=0), np.nanmean(oots, axis=0)
    dsnr = np.nanmax(dimg) / (np.nanstd(diffs, axis=0)[np.unravel_index(np.nanargmax(dimg), dimg.shape)] / np.sqrt(len(diffs)))
    sx, sy = fit_psf(dimg)
    ox, oy = fit_psf(oimg)
    per = []
    for d in diffs:
        try:
            per.append(fit_psf(d, guess=(sx, sy)))
        except Exception:
            pass
    per = np.array(per)
    per = per[(np.abs(per[:, 0] - sx) < 3) & (np.abs(per[:, 1] - sy) < 3)] if len(per) else per
    err = np.std(per, axis=0) / np.sqrt(len(per)) if len(per) >= 3 else np.array([np.nan, np.nan])
    # target position with proper motion
    h = tpf.hdu[0].header
    ra, dec = h["RA_OBJ"], h["DEC_OBJ"]
    pmra, pmdec = h.get("PMRA") or 0.0, h.get("PMDEC") or 0.0
    yrs = (Time(np.median(t) + 2457000, format="jd").jyear - 2000.0)
    dec2 = dec + pmdec / 3.6e6 * yrs
    ra2 = ra + pmra / 3.6e6 * yrs / np.cos(np.radians(dec))
    col, row = tpf.wcs.world_to_pixel(SkyCoord(ra2, dec2, unit="deg"))
    tx, ty = float(col), float(row)
    off = np.hypot(sx - tx, sy - ty)
    eoff = np.hypot(err[0], err[1])
    eoff = eoff if eoff > 1e-3 else np.nan          # identical per-transit fits = failed fits
    sig = off / eoff if np.isfinite(eoff) else np.nan
    rel = np.hypot(sx - ox, sy - oy)
    nearest = ""
    best = 9e9
    for s in nb:
        c, r = tpf.wcs.world_to_pixel(SkyCoord(s["ra"], s["dec"], unit="deg"))
        dd = np.hypot(sx - float(c), sy - float(r))
        s["px"] = (float(c), float(r))
        if dd < best:
            best, nearest = dd, f"TIC {s['tic']} (Tmag {s['tmag']:.1f}, {dd:.2f} px from source)"
    if not (np.isfinite(dsnr) and dsnr >= 5):
        res = "no detection"
    elif np.isfinite(sig) and sig < 3 and off < 0.5:
        res = "on target"
    elif np.isfinite(sig) and sig >= 3 and off >= 0.5:
        res = "off target"
    else:
        res = "inconclusive"
    return dict(sector=tpf.sector, n_transits=len(diffs), diff_snr=dsnr, src=(sx, sy), oot=(ox, oy),
                tgt=(tx, ty), off_px=off, off_err=eoff, off_sigma=sig, rel_oot_px=rel,
                nearest=nearest if best < 2.5 else "", result=res, dimg=dimg, oimg=oimg)


if __name__ == "__main__":
    lst, outdir = sys.argv[1:3]
    os.makedirs(outdir, exist_ok=True)
    hunt = pd.read_csv(HUNT)
    rows = []
    out_csv = f"{outdir}/diffimage.csv"
    done = set(pd.read_csv(out_csv).TOI.round(2)) if os.path.exists(out_csv) else set()
    from astroquery.mast import Catalogs
    for toi in pd.read_csv(lst).TOI:
        if round(toi, 2) in done:
            continue
        h = hunt[hunt.TOI == toi].iloc[0]
        tic = int(str(h.star).split()[-1])
        D = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
        recs = []
        try:
            tpfs = download_tpfs(tic)
            c0 = SkyCoord(tpfs[0].hdu[0].header["RA_OBJ"], tpfs[0].hdu[0].header["DEC_OBJ"], unit="deg")
            tt = Catalogs.query_region(c0, radius=60 * u.arcsec, catalog="TIC")
            nb = [dict(tic=int(r["ID"]), ra=float(r["ra"]), dec=float(r["dec"]), tmag=float(r["Tmag"])) for r in tt]
            for tpf in tpfs:
                r = sector_result(tpf, tic, h.t0, h.period, D, [dict(s) for s in nb])
                if r:
                    recs.append(r)
        except Exception as e:
            pd.DataFrame([dict(TOI=toi, tic=tic, result="error", note=f"{type(e).__name__}: {e}"[:200])]).to_csv(
                out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
            print(f"TOI {toi:.2f}: error {type(e).__name__}", flush=True)
            continue
        if not recs:
            pd.DataFrame([dict(TOI=toi, tic=tic, result="no usable sector", note="< 3 usable transits")]).to_csv(
                out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
            print(f"TOI {toi:.2f}: no usable sector", flush=True)
            continue
        fig, ax = plt.subplots(len(recs), 2, figsize=(10, 4.6 * len(recs)), squeeze=False)
        for i, r in enumerate(recs):
            for j, (img, title) in enumerate(((r["oimg"], "mean out-of-transit"), (r["dimg"], "difference (out − in)"))):
                ax[i, j].imshow(img, origin="lower", cmap="viridis")
                ax[i, j].plot(*r["tgt"], "w*", ms=14, mec="k", label="target (PM-corrected)")
                ax[i, j].plot(*r["src"], "rx", ms=12, mew=2, label="dimming source")
                ax[i, j].set_title(f"S{r['sector']} {title}")
            ax[i, 1].legend(loc="upper right", fontsize=7)
            ax[i, 1].set_xlabel(f"offset {r['off_px']:.2f}±{r['off_err']:.2f} px ({r['off_sigma']:.1f}σ) → {r['result']}")
        fig.suptitle(f"TOI {toi:.2f} difference image (diagnostic)")
        fig.tight_layout(); fig.savefig(f"{outdir}/TOI_{toi:.2f}_diffimage.png", dpi=100); plt.close(fig)
        res = [r["result"] for r in recs if r["result"] != "no detection"]
        overall = ("off target" if "off target" in res and "on target" not in res else
                   "on target" if "on target" in res and "off target" not in res else
                   "mixed" if "on target" in res else "inconclusive" if res else "no detection")
        row = dict(TOI=toi, tic=tic, overall=overall,
                   sectors=" | ".join(f"S{r['sector']}: {r['result']}, {r['off_px']:.2f}±{r['off_err']:.2f} px "
                                      f"({r['off_sigma']:.1f}σ), rel. OOT {r['rel_oot_px']:.2f} px, n={r['n_transits']}, "
                                      f"diff S/N {r['diff_snr']:.0f}" + (f", nearest {r['nearest']}" if r['nearest'] else "")
                                      for r in recs), result=overall)
        pd.DataFrame([row]).to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
        print(f"TOI {toi:.2f}: {overall}", flush=True)
