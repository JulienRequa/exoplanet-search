"""Measure the real red-noise level on DEV-set light curves (never blind_v2) and translate it into
the simulation's red-noise amplitude. Diagnostic for v3; v2 frozen. Julien's instruction 2026-11-27.

    python scripts/measure_rednoise_dev.py OUTDIR            (resumable)
    python scripts/measure_rednoise_dev.py OUTDIR --calibrate

Stars: dev rows of results/devfull200_multistart_2026-11-25.csv (dev3) and
results/devfull93_multistart_2026-11-25.csv (dev2 / v1-dev / known planets) where the v2 fit found
the catalog period ("P match" = yes), so the fitted ephemeris masks the real transits. Any TIC in
data/blind_v2.csv is excluded (checked). Same sectors as those dev runs (first N SPOC 2-min sectors),
cleaned exactly as the search sees them (pipeline.clean_all with the transits masked).
Masked before measuring: primary and secondary eclipse phases, ±1 × T14.

Measured per star (scripts/rednoise_snr.py statistic, Winn et al. 2008 / Pont et al. 2006):
  beta at fixed bin lengths 1 h, 2 h, 3 h, and at the star's own T14; sigma_1 (point scatter, ppm).

Translation to the simulation (--calibrate): pure-noise simulated stars (2 sectors, white noise
3–5 ppt as in the threshold simulation, AR(1) with tau 1 h and amplitude alpha × white; new seeds
30000+i), cleaned the same way, 60 stars per alpha on a grid 0–1.0 → median beta(2 h) versus alpha.

PRE-REGISTERED CHOICE OF LEVELS (written 2026-11-27, before measuring):
  "typical" level = alpha whose simulated median beta(2 h) equals the dev median beta(2 h);
  "pessimistic" level = alpha matching the dev 90th-percentile beta(2 h).
  (2 h ≈ the median transit duration in the threshold simulation's P 2–9 d range.)
  Both are reported; the threshold simulation is rerun at both. Also reported (not used for the
  choice): planet hosts vs FP stars, and stars with sigma_1 in 2,000–7,000 ppm (closest to the
  simulated faint stars).
"""
import os, sys, io, contextlib, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "."); sys.path.insert(0, "scripts"); sys.path.insert(0, "sims")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
import lightkurve as lk
import pipeline as pl
from rednoise_snr import _rstd
from report_cards import t14

BINS_H = (1.0, 2.0, 3.0)


def beta_at(t, r, L, min_fill=0.8):
    dt = np.median(np.diff(t)); s1 = _rstd(r)
    idx = np.floor((t - t[0]) / L).astype(int)
    u, inv, cnt = np.unique(idx, return_inverse=True, return_counts=True)
    good = cnt >= min_fill * L / dt
    if good.sum() < 10:
        return np.nan
    means = np.bincount(inv, weights=r)[good] / cnt[good]
    M, n = good.sum(), cnt[good].mean()
    return _rstd(means) / (s1 / np.sqrt(n) * np.sqrt(M / (M - 1)))


def measure(t, f, t0, P, T14):
    t = np.asarray(t, float); f = np.asarray(f, float)
    ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]
    ph = ((t - t0) / P) % 1
    dist = np.minimum.reduce([ph, 1 - ph, np.abs(ph - 0.5)]) * P
    keep = dist > T14
    t, r = t[keep], f[keep] - np.median(f[keep])
    out = {f"beta_{h:g}h": beta_at(t, r, h / 24) for h in BINS_H}
    out["beta_T14"] = beta_at(t, r, T14)
    out["sigma1_ppm"] = _rstd(r) * 1e6
    out["n_points"] = r.size
    return out


def stars():
    blind = set(pd.read_csv("data/blind_v2.csv")["TIC ID"].astype(int))
    rows = []
    for f in ("results/devfull200_multistart_2026-11-25.csv", "results/devfull93_multistart_2026-11-25.csv"):
        d = pd.read_csv(f)
        if "P match" in d:
            d = d[d["P match"] == "yes"]
        d = d[d.t0.notna()]
        rows.append(d)
    d = pd.concat(rows)
    d["tic"] = pd.to_numeric(d.star.astype(str).str.replace("TIC ", "", regex=False), errors="coerce")
    btoi = set(pd.read_csv("data/blind_v2.csv").TOI.round(2))
    assert not (set(d.tic.dropna().astype(int)) & blind), "blind star in dev list"
    assert not (set(d.TOI.dropna().round(2)) & btoi), "blind TOI in dev list"
    return d.drop_duplicates("star")


def run(out):
    os.makedirs(out, exist_ok=True); csv = f"{out}/dev_rednoise.csv"
    done = set(pd.read_csv(csv).star) if os.path.exists(csv) else set()
    for _, h in stars().iterrows():
        if h.star in done:
            continue
        rec = dict(TOI=h.TOI, star=h.star, tic=h.tic, kind=h.kind, period=h.period)
        try:
            T14 = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
            n = int(h.sectors_requested) if np.isfinite(h.sectors_requested) else 2
            with contextlib.redirect_stdout(io.StringIO()):
                lcs = pl.download(str(h.star), n)
                cl = pl.clean_all(lcs, period=h.period, t0=h.t0, duration=T14)
            rec.update(T14_h=T14 * 24, n_sectors=len(lcs),
                       **measure(cl.time.value, cl.flux.value, h.t0, h.period, T14))
        except Exception as e:
            rec.update(error=f"{type(e).__name__}: {e}"[:120])
        pd.DataFrame([rec]).to_csv(csv, mode="a", header=not os.path.exists(csv), index=False)
        print(f"{h.kind:13s} {h.star}: beta2h {rec.get('beta_2h', np.nan):.2f}", flush=True)


def calibrate(out, alphas=(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0), n=60):
    from rednoise_snr_sim import make
    os.makedirs(out, exist_ok=True)
    res = []
    for al in alphas:
        for i in range(n):
            seed = 30000 + i
            lcs, P, _ = make(seed, "noise", 2, al > 0, al, 1.0)
            with contextlib.redirect_stdout(io.StringIO()):
                cl = pl.clean_all(lcs)
            # mask a fake "transit" ephemeris of 2 h like the dev measurement
            res.append(dict(alpha=al, seed=seed, **measure(cl.time.value, cl.flux.value, 0.7, P, 2 / 24)))
        print(f"alpha {al}: median beta2h {np.median([r['beta_2h'] for r in res if r['alpha'] == al]):.2f}", flush=True)
    pd.DataFrame(res).to_csv(f"{out}/sim_alpha_calibration.csv", index=False)


if __name__ == "__main__":
    out = sys.argv[1]
    if "--calibrate" in sys.argv:
        calibrate(out)
    else:
        run(out)
