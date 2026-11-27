"""Validation of the red-noise-aware detection SNR (scripts/rednoise_snr.py) on the multi-sector
threshold simulation, and the threshold re-run at a measured red-noise level. SIMULATION ONLY;
v2 frozen. Julien's instruction 2026-11-27.

    python sims/rednoise_snr_sim.py OUT.csv --alpha A [--tau-h 1] [--summary]

Same simulated stars as sims/multisector_threshold_sim.py (same seeds, same draw order; with
--alpha 0.5 --tau-h 1 the stars are identical, so snr_white must reproduce that run's snr; this
is checked in the summary). Red noise: AR(1), correlation time tau, amplitude alpha × white noise.
Per configuration (2 / 4 sectors × white / white+red): 1,000 noise-only stars, 250 faint planets.
For each star: snr_white (v2-style), beta (time-averaging), snr_red = snr_white / beta.

PRE-REGISTERED VALIDATION CRITERIA (written 2026-11-27, before the validation run at alpha 0.5):
  V1 (no penalty on white noise): in the white-noise configurations, median beta ≤ 1.10, and
     faint-planet recovery at the FAR-1% threshold with snr_red is within 3 percentage points of
     the recovery with snr_white.
  V2 (threshold no longer depends on the noise colour): for each sector count, the FAR-1%
     threshold of snr_red in white+red noise is within 25% of the one in white noise
     (snr_white: 8.4 vs 28.6, i.e. +240%).
  V3 (transfer): applying the white-noise FAR-1% snr_red threshold to the white+red noise stars
     gives a false-alarm rate ≤ 3%.
  PASS = V1 and V2 and V3. If FAIL, the method is not adopted as is; any change is a new,
  separately committed version, validated on new seeds.
Which false-alarm rate to adopt is Julien's decision; this script only measures.
"""
import os, sys, io, contextlib, time, argparse
sys.path.insert(0, "."); sys.path.insert(0, "sims"); sys.path.insert(0, "scripts")
import numpy as np
import pandas as pd
import lightkurve as lk
import pipeline as pl
_argv, sys.argv = sys.argv, sys.argv[:1]   # the prototype reads sys.argv[1] at import time
from multisector_search_sim import multi, sector, STARTS
sys.argv = _argv
from rednoise_snr import snr_rednoise
from multisector_threshold_sim import _binom_ci

N_NOISE, N_PLANET = 1000, 250


def make(seed, kind, n_sec, red, alpha, tau_h):
    """Identical draw order to multisector_threshold_sim.make."""
    rng = np.random.default_rng(seed)
    P = rng.uniform(2, 9); t0 = rng.uniform(0.3, 1.5); a = 4.2 * P ** (2 / 3)
    noise = rng.uniform(3e-3, 5e-3)
    k, b = rng.uniform(0.03, 0.065), rng.uniform(0, 0.6)
    lcs = []
    for s0 in STARTS[n_sec]:
        T = sector(s0)
        m = pl.transit_model(T, t0, P, k, a, b, 0.4, 0.2) if kind == "planet" else np.ones(T.size)
        f = m + rng.normal(0, noise, T.size)
        if red:
            phi = np.exp(-(2 / 1440) / (tau_h / 24))
            e = rng.normal(0, alpha * noise * np.sqrt(1 - phi ** 2), T.size)
            r = np.empty(T.size); r[0] = e[0] / np.sqrt(1 - phi ** 2)
            for i in range(1, T.size):
                r[i] = phi * r[i - 1] + e[i]
            f = f + r
        lc = lk.LightCurve(time=T, flux=f, flux_err=np.full(T.size, noise))
        lc.meta["SECTOR"] = len(lcs) + 1
        lcs.append(lc)
    return lcs, P, k


def one(seed, kind, n_sec, red, alpha, tau_h):
    lcs, P, k = make(seed, kind, n_sec, red, alpha, tau_h)
    t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            popt, cl = multi(lcs, P)
        tt = np.asarray(cl.time.value, float); ff = np.asarray(cl.flux.value, float)
        sw, beta, sr = snr_rednoise(tt, ff, pl.transit_model(tt, *popt), popt[0], popt[1])
        good_P = abs(popt[1] - P) / P < 0.01
    except Exception:
        sw, beta, sr, good_P = 0.0, np.nan, 0.0, False
    return dict(kind=kind, n_sec=n_sec, red=red, alpha=alpha if red else 0.0, tau_h=tau_h, seed=seed,
                P=P, k=k if kind == "planet" else np.nan, snr_white=sw, beta=beta, snr_red=sr,
                correct_period=good_P, runtime_s=time.time() - t)


def thr_rec(nz, pl_, col, far):
    thr = np.quantile(nz, 1 - far)
    return thr, np.mean((pl_[col] > thr) & pl_.correct_period)


def summary(csv, ref="results/multisector_threshold_sim_2026-11-27.csv"):
    d = pd.read_csv(csv)
    L = [f"# Results: sims/rednoise_snr_sim.py, {os.path.basename(csv)} ({time.strftime('%Y-%m-%d')})", ""]
    if os.path.exists(ref) and (d[d.red].alpha.round(3) == 0.5).all() and (d.tau_h == 1).all():
        r = pd.read_csv(ref).merge(d, on=["kind", "n_sec", "red", "seed"])
        L += [f"Reproducibility vs {ref}: {len(r)} matched stars, max |snr_white − snr| = "
              f"{np.nanmax(np.abs(r.snr_white - r.snr)):.2e}", ""]
    thr1 = {}
    for (n_sec, red), g in d.groupby(["n_sec", "red"]):
        nz = g[g.kind == "noise"]; pl_ = g[g.kind == "planet"]
        L.append(f"## {n_sec} sectors, {'white + red (alpha %.2f, tau %.1f h)' % (g.alpha.iloc[0], g.tau_h.iloc[0]) if red else 'white'} noise "
                 f"({len(nz)} noise stars, {len(pl_)} faint planets)")
        L.append(f"- beta on noise stars: median {nz.beta.median():.2f} (16–84%: "
                 f"{nz.beta.quantile(.16):.2f}–{nz.beta.quantile(.84):.2f}); on planets median {pl_.beta.median():.2f}")
        for col in ("snr_white", "snr_red"):
            L.append(f"- {col}: FAR at 7.1 = {np.mean(nz[col] > 7.1):.1%}; "
                     f"planets recovered at 7.1 = {np.mean((pl_[col] > 7.1) & pl_.correct_period):.0%}")
            for far in (0.05, 0.01, 0.0025):
                thr, rec = thr_rec(nz[col].values, pl_, col, far)
                n_fa = int(np.sum(nz[col] > thr)); lo, hi = [x / len(nz) for x in _binom_ci(n_fa, len(nz))]
                L.append(f"  - FAR {far:.2%}: threshold {thr:.2f} (FAR 95% CI {lo:.2%}–{hi:.2%}); "
                         f"faint planets recovered {rec:.0%}")
                if far == 0.01:
                    thr1[(n_sec, red, col)] = (thr, rec)
        L.append("")
    # validation criteria
    L.append("## Pre-registered validation")
    ok = True
    for n_sec in sorted(d.n_sec.unique()):
        gw = d[(d.n_sec == n_sec) & ~d.red]; gr = d[(d.n_sec == n_sec) & d.red]
        if not len(gw) or not len(gr):
            continue
        mb = gw[gw.kind == "noise"].beta.median()
        dr = thr1[(n_sec, False, "snr_red")][1] - thr1[(n_sec, False, "snr_white")][1]
        v1 = mb <= 1.10 and abs(dr) <= 0.03
        tw, tr = thr1[(n_sec, False, "snr_red")][0], thr1[(n_sec, True, "snr_red")][0]
        v2 = abs(tr / tw - 1) <= 0.25
        far3 = np.mean(gr[gr.kind == "noise"].snr_red > tw)
        v3 = far3 <= 0.03
        ok &= v1 and v2 and v3
        L.append(f"- {n_sec} sectors: V1 median beta (white) {mb:.2f}, recovery change {dr:+.0%} → "
                 f"{'pass' if v1 else 'FAIL'}; V2 thresholds white {tw:.2f} vs red {tr:.2f} "
                 f"({tr / tw - 1:+.0%}) → {'pass' if v2 else 'FAIL'}; V3 FAR on red stars at the white "
                 f"threshold {far3:.1%} → {'pass' if v3 else 'FAIL'}")
    L += ["", f"**Validation: {'PASS' if ok else 'FAIL'}**"]
    return "\n".join(L)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out"); ap.add_argument("--alpha", type=float, required=True)
    ap.add_argument("--tau-h", type=float, default=1.0); ap.add_argument("--summary", action="store_true")
    ap.add_argument("--configs", default="2w,2r,4w,4r")
    a = ap.parse_args()
    if a.summary:
        print(summary(a.out)); sys.exit()
    done = set()
    if os.path.exists(a.out):
        d = pd.read_csv(a.out); done = set(zip(d.kind, d.n_sec, d.red, d.seed))
    for c in a.configs.split(","):
        n_sec, red = int(c[0]), c[1] == "r"
        off = (n_sec * 10 + int(red)) * 100000
        jobs = [("noise", off + 10000 + i) for i in range(N_NOISE)] + \
               [("planet", off + 20000 + i) for i in range(N_PLANET)]
        for kind, seed in jobs:
            if (kind, n_sec, red, seed) in done:
                continue
            r = one(seed, kind, n_sec, red, a.alpha, a.tau_h)
            pd.DataFrame([r]).to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
        print(f"done {n_sec} sectors red={red}", flush=True)
