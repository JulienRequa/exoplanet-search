"""Detection threshold for the v3 multi-sector search, from the simulated noise false-alarm rate
(Julien's rule, 2026-11-25). SIMULATION ONLY; not in the pipeline; v2 frozen.

    python sims/multisector_threshold_sim.py OUT.csv          (resumable)

Uses the prototype search from sims/multisector_search_sim.py (`multi`, stitched BLS within ±2%
of a catalog period, then the v2 fit) and its detection SNR (`snr`, as in v2's vet_star).

Configurations: 2 and 4 sectors (two pairs 700 d apart), × white noise only / white + red noise
(AR(1), 1-hour correlation time, amplitude 50% of the white noise — similar to the red noise in
the odd/even simulations). Per configuration:
  - 1,000 noise-only stars with a random "catalog" period 2–9 d → distribution of the SNR the
    search returns. That distribution is the false-alarm curve.
  - 250 faint planets (Rp/Rs 0.03–0.065, b 0–0.6, P 2–9 d) → SNR at the correct period, to show
    what each threshold costs in recovered planets.
Seeds fixed: noise stars 10000+i, planets 20000+i (per configuration offset).
Summary (sims/multisector_threshold_sim_results_<date>.md, written by --summary): for each
configuration, the SNR threshold giving a false-alarm rate of 5%, 1% and 0.25% (with a binomial
95% interval on the rate at that threshold), and the fraction of faint planets recovered at each.
Which false-alarm rate to adopt is Julien's decision; this script only measures.
"""
import os, sys, io, contextlib, time
sys.path.insert(0, "."); sys.path.insert(0, "sims")
import numpy as np
import pandas as pd
import lightkurve as lk
import pipeline as pl
_argv, sys.argv = sys.argv, sys.argv[:1]   # the prototype reads sys.argv[1] at import time
from multisector_search_sim import multi, snr, sector, STARTS
sys.argv = _argv

N_NOISE, N_PLANET = 1000, 250


def make(seed, kind, n_sec, red):
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
            phi = np.exp(-(2 / 1440) / (1 / 24))           # AR(1), 1 h correlation
            e = rng.normal(0, 0.5 * noise * np.sqrt(1 - phi ** 2), T.size)
            r = np.empty(T.size); r[0] = e[0] / np.sqrt(1 - phi ** 2)
            for i in range(1, T.size):
                r[i] = phi * r[i - 1] + e[i]
            f = f + r
        lc = lk.LightCurve(time=T, flux=f, flux_err=np.full(T.size, noise))
        lc.meta["SECTOR"] = len(lcs) + 1
        lcs.append(lc)
    return lcs, P, k


def one(seed, kind, n_sec, red):
    lcs, P, k = make(seed, kind, n_sec, red)
    t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            popt, cl = multi(lcs, P)
        s = snr(cl, popt)
        good_P = abs(popt[1] - P) / P < 0.01
    except Exception:
        s, good_P = 0.0, False
    return dict(kind=kind, n_sec=n_sec, red=red, seed=seed, P=P, k=k if kind == "planet" else np.nan,
                snr=s, correct_period=good_P, runtime_s=time.time() - t)


def summary(csv):
    d = pd.read_csv(csv)
    lines = [f"# Results: sims/multisector_threshold_sim.py ({time.strftime('%Y-%m-%d')})", ""]
    for (n_sec, red), g in d.groupby(["n_sec", "red"]):
        nz = g[g.kind == "noise"].snr.values; pl_ = g[g.kind == "planet"]
        lines.append(f"## {n_sec} sectors, {'white + red' if red else 'white'} noise "
                     f"({len(nz)} noise stars, {len(pl_)} faint planets)")
        lines.append(f"- false-alarm rate at SNR 7.1 (v2 threshold): {np.mean(nz > 7.1):.1%}")
        for far in (0.05, 0.01, 0.0025):
            thr = np.quantile(nz, 1 - far)
            n_fa = int(np.sum(nz > thr)); n = len(nz)
            lo, hi = [x / n for x in _binom_ci(n_fa, n)]
            rec = np.mean((pl_.snr > thr) & pl_.correct_period)
            lines.append(f"- FAR {far:.2%}: threshold SNR {thr:.2f} (FAR 95% CI {lo:.2%}–{hi:.2%}); "
                         f"faint planets recovered {rec:.0%}")
        rec71 = np.mean((pl_.snr > 7.1) & pl_.correct_period)
        lines.append(f"- faint planets recovered at SNR 7.1: {rec71:.0%}; median runtime "
                     f"{g.runtime_s.median():.1f} s/star")
        lines.append("")
    return "\n".join(lines)


def _binom_ci(k, n, z=1.96):
    p = k / n; den = 1 + z ** 2 / n
    c = (p + z ** 2 / (2 * n)) / den; h = z * np.sqrt(p * (1 - p) / n + z ** 2 / (4 * n ** 2)) / den
    return (max(0, c - h) * n, min(1, c + h) * n)


if __name__ == "__main__":
    out = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == "--summary":
        print(summary(out)); sys.exit()
    done = set()
    if os.path.exists(out):
        d = pd.read_csv(out); done = set(zip(d.kind, d.n_sec, d.red, d.seed))
    for n_sec in (2, 4):
        for red in (False, True):
            off = (n_sec * 10 + int(red)) * 100000
            jobs = [("noise", off + 10000 + i) for i in range(N_NOISE)] + \
                   [("planet", off + 20000 + i) for i in range(N_PLANET)]
            for kind, seed in jobs:
                if (kind, n_sec, red, seed) in done:
                    continue
                r = one(seed, kind, n_sec, red)
                pd.DataFrame([r]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
            print(f"done {n_sec} sectors red={red}", flush=True)
