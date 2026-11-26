"""Does the secondary eclipse repeat in other TESS sectors? (diagnostic; v2 unchanged)

    python scripts/extra_sector_secondary.py OUTDIR TOI [TOI ...]

For every TESS sector of the star (any cadence), one light curve per sector is chosen in this
order: SPOC 2-min > TESS-SPOC FFI > QLP. Eleanor, T16 and TGLC are only used if nothing else
exists for that sector. Each sector is measured on its own:
 - flux: PDCSAP (SPOC / TESS-SPOC) or SAP (QLP, the raw aperture flux, so all products get the
   same detrending); quality-flagged cadences removed (for QLP: its own QUALITY flag must be 0 —
   added after a first run showed scattered-light glitches of 10⁴–10⁵ ppm in S34/S35); upward
   outliers > 4σ clipped; downward points > 10σ removed (gross glitches only; far above any
   eclipse depth here).
 - detrending: the v2 method (lightkurve flatten, transits masked over 1.5 × T14) with the window
   set to ~1 day in that sector's cadence (v2 uses 721 two-minute points ≈ 1 day).
 - transit timing: per sector, the time shift (±0.45 P around the ExoFOP ephemeris, 2-min
   steps) that gives the deepest local transit, so the secondary window follows the observed
   transits. A sector whose recovered transit is < 50% of the median over the star's sectors is
   marked "transit not recovered" and excluded from the combined secondary.
 - secondary: the report-card method — depth at phase 0.5 (±0.4·T14) against its local baseline
   (0.75–2·T14), significance against a red-noise null of 200 random phases (seed 0).
 - also: the same local depth scanned over phase 0.30–0.70 (in case the orbit is eccentric).
Light curves are binned to 10 min when the cadence is shorter, so every sector is treated alike
(this changes no depth, only the noise per point).
"""
import os, sys, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "."); sys.path.insert(0, "scripts")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import lightkurve as lk
import pipeline as pl
from report_cards import t14, secondary_local, window_depth, centres

HUNT = "results/hunt_v2_2026-11-25.csv"
TOIS = "data/toi_2026-11-26.csv"
RANK = {"SPOC": 0, "TESS-SPOC": 1, "QLP": 2}


def pick_products(tic):
    sr = lk.search_lightcurve(f"TIC {tic}", mission="TESS")
    best = {}
    for i, r in enumerate(sr.table):
        sec = int(str(r["mission"]).split()[-1])
        rank = RANK.get(str(r["author"]), 9)
        if sec not in best or rank < best[sec][0]:
            best[sec] = (rank, i, str(r["author"]), float(r["exptime"]))
    return sr, dict(sorted(best.items()))


def load(sr, i, author):
    lc = sr[i].download(quality_bitmask="default")
    if author == "QLP" and "quality" in lc.colnames:
        lc = lc[np.asarray(lc["quality"].value if hasattr(lc["quality"], "value") else lc["quality"]) == 0]
    if author == "QLP":
        col = "sap_flux" if "sap_flux" in lc.colnames else "flux"
        lc = lk.LightCurve(time=lc.time, flux=lc[col], flux_err=lc["flux_err"] if "flux_err" in lc.colnames else None)
    lc = lc.remove_nans().normalize()
    lc = lc.remove_outliers(sigma_upper=4, sigma_lower=np.inf)
    return lc.remove_outliers(sigma_upper=np.inf, sigma_lower=10)   # gross glitches only


def transit_shift(t, f, popt, P, D):
    """Time shift (±0.45 P, 2-min steps) that maximises the mean local transit depth.
    (First version used a χ² grid over ±0.3 d; ExoFOP ephemerides drift by more than that for
    some stars, e.g. TOI 4694.01 by ~0.8 d at sector 98.)"""
    best, bs = -np.inf, 0.0
    for s in np.arange(-0.45 * P, 0.45 * P, 2 / 1440):
        d = window_depth(t, f, centres(t, popt[0] + s, P, 0.0), 0.35 * D, 0.75 * D, 2.0 * D)
        if len(d) and d.mean() > best:
            best, bs = d.mean(), s
    return bs


def measure(toi, outdir):
    hunt, tois = pd.read_csv(HUNT), pd.read_csv(TOIS)
    h = hunt[hunt.TOI.round(2) == round(toi, 2)].iloc[0]
    x = tois[tois.TOI.round(2) == round(toi, 2)].iloc[0]
    tic = int(str(h.star).split()[-1])
    P = float(x["Period (days)"])
    t0 = float(x["Epoch (BJD)"]) - 2457000
    D = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
    popt = [t0, P, h["Rp/Rs"], h["a/Rs"], h.b, h.u1, h.u2]
    sr, prods = pick_products(tic)
    rows, curves = [], []
    for sec, (rank, i, author, exp) in prods.items():
        try:
            lc = load(sr, i, author)
        except Exception as e:
            rows.append(dict(TOI=toi, sector=sec, product=author, error=type(e).__name__)); continue
        if exp < 600:
            lc = lc.bin(time_bin_size=10 / 1440).remove_nans()
        cad = np.median(np.diff(lc.time.value))
        win = int(round(1.0 / cad)) | 1
        t = np.asarray(lc.time.value, float)
        ph = ((t - t0) / P + 0.5) % 1 - 0.5
        mask = np.abs(ph * P) < 0.75 * D
        flat = lc.flatten(window_length=max(win, 7), mask=mask)
        f = np.asarray(flat.flux.value, float); ok = np.isfinite(f); t, f = t[ok], f[ok]
        s = transit_shift(t, f, popt, P, D)
        tc0 = t0 + s
        rng = np.random.default_rng(0)
        tr = window_depth(t, f, centres(t, tc0, P, 0.0), 0.4 * D, 0.75 * D, 2.0 * D)
        sec_d, sec_sd, sec_sig = secondary_local(t, f, tc0, P, D, rng)
        scan = []
        for phs in np.arange(0.30, 0.705, 0.01):
            d = window_depth(t, f, centres(t, tc0, P, phs * P), 0.4 * D, 0.75 * D, 2.0 * D)
            scan.append(d.mean() if len(d) else np.nan)
        rows.append(dict(TOI=toi, tic=tic, sector=sec, product=author, cadence_s=exp,
                         used_by_v2=(author == "SPOC" and exp == 120), shift_min=s * 1440,
                         n_transits=len(tr), transit_ppm=tr.mean() * 1e6 if len(tr) else np.nan,
                         sec_ppm=sec_d * 1e6, sec_err_ppm=sec_sd * 1e6, sec_sigma=sec_sig,
                         ratio=sec_d / tr.mean() if len(tr) and tr.mean() > 0 else np.nan))
        curves.append((sec, author, np.arange(0.30, 0.705, 0.01), np.array(scan), sec_sd))
        print(f"TOI {toi:.2f} S{sec} {author}: transit {rows[-1]['transit_ppm']:.0f} ppm ({len(tr)} tr), "
              f"secondary {sec_d*1e6:.0f} ± {sec_sd*1e6:.0f} ppm ({sec_sig:.1f}σ), shift {s*1440:.0f} min", flush=True)
    if curves:
        fig, ax = plt.subplots(figsize=(9, 5))
        for sec, au, phs, sc, sd in curves:
            ax.plot(phs, sc * 1e6, "o-", ms=3, label=f"S{sec} {au} (null σ {sd*1e6:.0f} ppm)")
        ax.axvline(0.5, color="k", ls=":"); ax.axhline(0, color="gray", lw=0.5)
        ax.set_xlabel("orbital phase"); ax.set_ylabel("local depth (ppm)")
        ax.set_title(f"TOI {toi:.2f}: local depth vs phase, each sector separately (diagnostic)")
        ax.legend(fontsize=7); fig.tight_layout()
        fig.savefig(f"{outdir}/TOI_{toi:.2f}_sectors.png", dpi=110); plt.close(fig)
    return rows


if __name__ == "__main__":
    outdir = sys.argv[1]; os.makedirs(outdir, exist_ok=True)
    out = f"{outdir}/extra_sectors.csv"
    for toi in [float(x) for x in sys.argv[2:]]:
        try:
            rows = measure(toi, outdir)
        except Exception as e:
            rows = [dict(TOI=toi, error=f"{type(e).__name__}: {e}"[:200])]
            print(f"TOI {toi:.2f}: error {rows[0]['error']}", flush=True)
        df = pd.DataFrame(rows)
        if "transit_ppm" in df:
            med = df.transit_ppm.median()
            df["transit_recovered"] = df.transit_ppm >= 0.5 * med
            ok = df[df.transit_recovered & df.sec_err_ppm.notna() & (df.sec_err_ppm > 0)]
            for label, sub in (("all", ok), ("extra", ok[~ok.used_by_v2.astype(bool)])):
                w = 1 / sub.sec_err_ppm ** 2
                if len(sub):
                    m, e = (sub.sec_ppm * w).sum() / w.sum(), w.sum() ** -0.5
                    print(f"TOI {toi:.2f} combined {label} sectors ({len(sub)}): {m:.0f} ± {e:.0f} ppm ({m/e:.1f}σ)", flush=True)
        df.to_csv(out, mode="a", header=not os.path.exists(out), index=False)
