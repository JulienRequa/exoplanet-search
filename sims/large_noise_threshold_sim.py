"""Large noise-only simulation for the v3 multi-sector detection threshold at the 0.25% false-alarm
rate. SIMULATION ONLY; v2 frozen; nothing here changes the pipeline. Julien's decision 2026-11-27:
"Provisional v3 multi-sector threshold: 1% false-alarm rate (SNR ≈ 9.1–9.2 at the typical measured
noise level). Target: 0.25%. Before adopting it, run a larger noise-only simulation (≥10,000 stars
per setup, new seeds, rules committed first)."

    python sims/large_noise_threshold_sim.py OUT.csv [--workers 6]     (resumable)
    python sims/large_noise_threshold_sim.py OUT.csv --summary

Same star model, search and SNR as sims/rednoise_snr_sim.py (`make`, `one`; prototype search
`multi` from sims/multisector_search_sim.py, v2-style snr_white). No code in those files changes.

PRE-REGISTERED RULES (written 2026-11-27, committed before the run):
  Setups: 2 and 4 sectors × (a) white noise only, (b) white + red noise at the measured TYPICAL level
    (AR(1), tau 1 h, alpha 0.091; results/rednoise_dev_2026-11-27/SUMMARY.md). The typical-level
    setups decide; the white setups are a reference only.
  Per setup: 10,000 noise-only stars + 1,000 faint planets (same distributions as before:
    P 2–9 d, Rp/Rs 0.03–0.065, b 0–0.6, white noise 3–5 ppt).
  New seeds: 50,000,000 + (10 × n_sectors + red) × 1,000,000 + i (noise) and + 500,000 + i
    (planets); disjoint from every earlier simulation (earlier seeds are all < 5,000,000; asserted).
  Threshold at FAR f (f = 5%, 1%, 0.25%): the empirical (1 − f) quantile of the noise-star SNR
    (numpy default, linear interpolation). Uncertainty: 95% bootstrap interval of that quantile
    (2,000 resamples, rng seed 12345) and the binomial 95% interval of the false-alarm rate at it.
  Faint-planet recovery: fraction with SNR above the threshold AND the correct period (±1%).
  Enough false alarms for 0.25% (the target), per typical-level setup: at least 20 noise stars
    above the threshold AND bootstrap 95% interval no wider than 1.0 SNR unit. If either fails,
    the result is reported as "not enough" and the run is extended by another 10,000 noise stars
    per setup (seeds continue at i = 10,000…); no other change.
  Proposed v3 number at each FAR: the HIGHER of the 2- and 4-sector typical-level thresholds,
    rounded UP to 0.1 (one threshold for all stars). Reported, not adopted: adoption is Julien's.
  Also reported (not a gate): agreement with the 1,000-star run (sims/rednoise_threshold_measured_
    levels_2026-11-27.md) and the beta distribution on the noise stars (input to the two-tier design).
"""
import os, sys, time, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
sys.path.insert(0, "."); sys.path.insert(0, "sims"); sys.path.insert(0, "scripts")
import numpy as np
import pandas as pd
from multiprocessing import Pool

ALPHA_TYPICAL, TAU_H = 0.091, 1.0
N_NOISE, N_PLANET = 10_000, 1_000
SETUPS = [(2, False), (2, True), (4, False), (4, True)]
FARS = (0.05, 0.01, 0.0025)


def seed_of(kind, n_sec, red, i):
    base = 50_000_000 + (10 * n_sec + int(red)) * 1_000_000
    return base + (500_000 + i if kind == "planet" else i)


def jobs(n_noise=N_NOISE, n_planet=N_PLANET):
    out = []
    for n_sec, red in SETUPS:
        out += [("noise", n_sec, red, seed_of("noise", n_sec, red, i)) for i in range(n_noise)]
        out += [("planet", n_sec, red, seed_of("planet", n_sec, red, i)) for i in range(n_planet)]
    assert min(j[3] for j in out) >= 50_000_000 and len({j[3] for j in out}) == len(out)
    return out


def _run(job):
    import io, contextlib
    from rednoise_snr_sim import one
    kind, n_sec, red, seed = job
    with contextlib.redirect_stdout(io.StringIO()):
        return one(seed, kind, n_sec, red, ALPHA_TYPICAL, TAU_H)


def boot_q(x, q, n=2000, seed=12345):
    rng = np.random.default_rng(seed)
    qs = [np.quantile(rng.choice(x, x.size), q) for _ in range(n)]
    return np.quantile(qs, [0.025, 0.975])


def summary(csv):
    from multisector_threshold_sim import _binom_ci
    d = pd.read_csv(csv)
    L = [f"# Results: sims/large_noise_threshold_sim.py ({time.strftime('%Y-%m-%d')})", "",
         f"Rows: {len(d)}; failed searches (SNR 0): {int((d.snr_white == 0).sum())}", ""]
    prop = {f: [] for f in FARS}; enough = True
    for (n_sec, red), g in d.groupby(["n_sec", "red"]):
        nz = g[g.kind == "noise"].snr_white.values; pl_ = g[g.kind == "planet"]
        lab = f"white + red (typical, alpha {ALPHA_TYPICAL})" if red else "white only (reference)"
        L.append(f"## {n_sec} sectors, {lab}: {nz.size} noise stars, {len(pl_)} faint planets")
        L.append(f"- false alarms at SNR 7.1 (v2): {np.mean(nz > 7.1):.1%}; beta on noise stars: median "
                 f"{g[g.kind == 'noise'].beta.median():.2f}, 90th pct {g[g.kind == 'noise'].beta.quantile(.9):.2f}")
        for f in FARS:
            thr = np.quantile(nz, 1 - f); lo, hi = boot_q(nz, 1 - f)
            k = int(np.sum(nz > thr)); flo, fhi = [x / nz.size for x in _binom_ci(k, nz.size)]
            rec = np.mean((pl_.snr_white > thr) & pl_.correct_period)
            line = (f"- FAR {f:.2%}: threshold {thr:.2f} (bootstrap 95% {lo:.2f}–{hi:.2f}); {k} noise stars above "
                    f"(FAR 95% CI {flo:.2%}–{fhi:.2%}); faint planets kept {rec:.1%}")
            if red and f == 0.0025:
                ok = k >= 20 and (hi - lo) <= 1.0; enough &= ok
                line += f" → enough false alarms: {'yes' if ok else 'NO'}"
            L.append(line)
            if red:
                prop[f].append(thr)
        L.append("")
    L.append("## Proposed v3 threshold (higher of 2 / 4 sectors at the typical level, rounded up to 0.1)")
    for f in FARS:
        if prop[f]:
            L.append(f"- FAR {f:.2%}: SNR {np.ceil(max(prop[f]) * 10) / 10:.1f}")
    L += ["", f"**0.25% rests on enough false alarms: {'YES' if enough else 'NO — extend by 10,000 noise stars per setup'}**",
          "", "Adoption is Julien's decision."]
    return "\n".join(L)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out"); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--summary", action="store_true"); ap.add_argument("--n-noise", type=int, default=N_NOISE)
    a = ap.parse_args()
    if a.summary:
        print(summary(a.out)); sys.exit()
    done = set()
    if os.path.exists(a.out):
        d = pd.read_csv(a.out); done = set(zip(d.kind, d.n_sec, d.red, d.seed))
    todo = [j for j in jobs(a.n_noise) if j not in done]
    print(f"{len(todo)} stars to run ({len(done)} done), {a.workers} workers", flush=True)
    t0 = time.time()
    with Pool(a.workers) as pool:
        for n, r in enumerate(pool.imap_unordered(_run, todo, chunksize=8), 1):
            pd.DataFrame([r]).to_csv(a.out, mode="a", header=not os.path.exists(a.out), index=False)
            if n % 500 == 0:
                el = time.time() - t0
                print(f"{n}/{len(todo)}  {el / 3600:.2f} h elapsed, ~{el / n * (len(todo) - n) / 3600:.1f} h left", flush=True)
    print("done", flush=True)
