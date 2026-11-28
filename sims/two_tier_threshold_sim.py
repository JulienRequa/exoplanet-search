"""Two-tier (three-bucket) detection threshold for the v3 multi-sector search: calibration and
validation. SIMULATION ONLY; v2 frozen. Design: docs/design_two_tier_threshold.md (approved by
Julien 2026-11-27 with two additions: conservative tier assignment using the uncertainty of the
red-noise factor; tier-3 signals listed separately, never claimed).

    python sims/two_tier_threshold_sim.py calibrate CAL.csv [--workers 6]       (resumable)
    python sims/two_tier_threshold_sim.py validate  VAL.csv --cal CAL.csv [--workers 6]
    python sims/two_tier_threshold_sim.py summary   VAL.csv --cal CAL.csv

Star model, search and SNR as sims/rednoise_snr_sim.py (`make`; prototype `multi`; v2-style SNR).

PRE-REGISTERED RULES (written 2026-11-27, committed before any run):
 Red-noise factor of a star (beta_2h): as scripts/measure_rednoise_dev.py `measure`/`beta_at`, on the
   cleaned light curve returned by the search, masking the FOUND signal's primary and secondary
   ±1 × T14 (T14 from the fitted model's in-transit points). M = number of 2-h bins used.
 Uncertainty: sigma_beta = beta × sqrt(1.36 / M) (asymptotic variance of the normalised MAD scale
   estimator, 1.36 sigma^2 / M; checked before writing: predicts 5.0% for 2-sector simulated stars,
   empirical 16–84% half-width 4.2–5.7% across alpha 0–1.0).
 Tier (conservative): beta_hi = beta + sigma_beta.  beta_hi ≤ 1.15 → normal; 1.15 < beta_hi ≤ 2.0 →
   strict; beta_hi > 2.0 (or beta not measurable) → tier 3. A star whose interval crosses a cutoff
   therefore goes to the stricter tier.
 Tier 3: searched like any star; a signal is LISTED separately (low confidence, follow-up only,
   never claimed, never counted in any rate or detection total) if its SNR ≥ the strict threshold.
   Its listing rate on pure noise is reported, not a gate.

 CALIBRATION (seeds 60,000,000 + (10·n_sec + tier)·1,000,000 + i; planets + 500,000 + i):
   Edge alphas by linear interpolation of the existing calibration table
   (results/rednoise_dev_2026-11-27/sim_alpha_calibration.csv, median beta_2h vs alpha) at
   beta 1.15 (normal) and 2.0 (strict). Per tier × 2 and 4 sectors: 10,000 noise stars + 1,000
   faint planets at the edge alpha (all stars used, whatever their measured beta).
   Tier threshold at FAR f (1% and 0.25%): (1 − f) quantile of the noise SNR, the higher of
   2 / 4 sectors, rounded UP to 0.1.

 VALIDATION (new seeds 70,000,000 + n_sec·1,000,000 + i; planets + 500,000 + i):
   Mixed population like the dev stars: each star's alpha from a random dev beta_2h
   (results/rednoise_dev_2026-11-27/dev_rednoise.csv, rng seeded by the star seed) mapped to alpha
   through the same table (beta < 0.90 → 0, > 3.61 → 1.0). 10,000 noise + 2,000 faint planets per
   sector count (2, 4). Every planet is also run as its noise-only twin (same seed: the noise draws
   are identical), to see whether the planet changes its own tier.
   Each star is tiered by its OWN measured beta_hi (not its true alpha), and a noise star is a
   false alarm if SNR > its tier's threshold (tier 3: never claimed).
   Criteria, for each f in (1%, 0.25%) and each sector count:
     T1 overall: claimed false alarms / all noise stars — binomial 95% interval reaches ≤ f.
     T2 per tier (normal, strict): false-alarm rate in the tier ≤ 2f.
     T3 worth it: share of faint planets claimed (correct period, normal or strict tier, above
        threshold) is HIGHER than with one single threshold on all stars set at the same overall
        false-alarm rate on this mixed population.
   Reported, not gates: tier shares (dev: ~74 / 18 / 8% before the conservative rule); planets
   whose tier differs from their noise-only twin; tier-3 listing rate on noise; sigma_beta formula
   vs the empirical spread of beta at fixed alpha in the calibration run.
   PASS = T1 and T2 and T3 everywhere. FAIL → not adopted as is; any change is a new committed
   version validated on new seeds. Adoption is Julien's decision.
"""
import os, sys, time, argparse
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
sys.path.insert(0, "."); sys.path.insert(0, "sims"); sys.path.insert(0, "scripts")
import numpy as np
import pandas as pd
from multiprocessing import Pool

CUT_NORMAL, CUT_STRICT = 1.15, 2.0
TAU_H = 1.0
FARS = (0.01, 0.0025)
CALTAB = "results/rednoise_dev_2026-11-27/sim_alpha_calibration.csv"
DEVTAB = "results/rednoise_dev_2026-11-27/dev_rednoise.csv"


def _table():
    c = pd.read_csv(CALTAB).groupby("alpha").beta_2h.median()
    return c.values, c.index.values            # beta (increasing), alpha


def beta_to_alpha(b):
    bb, aa = _table()
    return float(np.interp(b, bb, aa, left=0.0, right=1.0))




def edge_alpha(tier):
    return beta_to_alpha(CUT_NORMAL if tier == "normal" else CUT_STRICT)


def beta2h(t, f, t0, P, T14, L=2 / 24, min_fill=0.8):
    """Same computation as measure_rednoise_dev.measure + beta_at, also returning M."""
    from rednoise_snr import _rstd
    t = np.asarray(t, float); f = np.asarray(f, float)
    ok = np.isfinite(t) & np.isfinite(f); t, f = t[ok], f[ok]
    ph = ((t - t0) / P) % 1
    dist = np.minimum.reduce([ph, 1 - ph, np.abs(ph - 0.5)]) * P
    keep = dist > T14
    t, r = t[keep], f[keep] - np.median(f[keep])
    dt = np.median(np.diff(t)); s1 = _rstd(r)
    idx = np.floor((t - t[0]) / L).astype(int)
    u, inv, cnt = np.unique(idx, return_inverse=True, return_counts=True)
    good = cnt >= min_fill * L / dt
    if good.sum() < 10:
        return np.nan, 0
    means = np.bincount(inv, weights=r)[good] / cnt[good]
    M, n = int(good.sum()), cnt[good].mean()
    return _rstd(means) / (s1 / np.sqrt(n) * np.sqrt(M / (M - 1))), M


def tier_of(beta, M):
    if not np.isfinite(beta) or M < 10:
        return "tier3"
    hi = beta * (1 + np.sqrt(1.36 / M))
    return "normal" if hi <= CUT_NORMAL else ("strict" if hi <= CUT_STRICT else "tier3")


def run_star(job):
    import io, contextlib
    import pipeline as pl
    from rednoise_snr_sim import make
    from rednoise_snr import snr_white
    _a, sys.argv = sys.argv, sys.argv[:1]
    from multisector_search_sim import multi
    sys.argv = _a
    kind, n_sec, alpha, seed, extra = job
    lcs, P, k = make(seed, kind, n_sec, alpha > 0, alpha, TAU_H)
    rec = dict(kind=kind, n_sec=n_sec, alpha=alpha, seed=seed, P=P, k=k if kind == "planet" else np.nan, **extra)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            popt, cl = multi(lcs, P)
        tt = np.asarray(cl.time.value, float); ff = np.asarray(cl.flux.value, float)
        model = pl.transit_model(tt, *popt)
        sw = snr_white(tt, ff, model)
        intr = model < 1
        ph = ((tt - popt[0]) / popt[1] + 0.5) % 1 - 0.5
        T14 = 2 * np.max(np.abs(ph[intr])) * popt[1] if intr.any() else 2 / 24
        b, M = beta2h(tt, ff, popt[0], popt[1], T14)
        rec.update(snr=sw, correct_period=abs(popt[1] - P) / P < 0.01, beta2h=b, M=M)
    except Exception as e:
        rec.update(snr=0.0, correct_period=False, beta2h=np.nan, M=0, error=type(e).__name__)
    rec["tier"] = tier_of(rec["beta2h"], rec["M"])
    return rec


def cal_jobs():
    out = []
    for ti, tier in enumerate(("normal", "strict")):
        a = edge_alpha(tier)
        for n_sec in (2, 4):
            base = 60_000_000 + (10 * n_sec + ti) * 1_000_000
            out += [("noise", n_sec, a, base + i, {"cal_tier": tier}) for i in range(10_000)]
            out += [("planet", n_sec, a, base + 500_000 + i, {"cal_tier": tier}) for i in range(1_000)]
    return out


def val_jobs():
    devb = pd.read_csv(DEVTAB).beta_2h.dropna().values
    out = []
    for n_sec in (2, 4):
        base = 70_000_000 + n_sec * 1_000_000
        for kind, off, n in (("noise", 0, 10_000), ("planet", 500_000, 2_000)):
            for i in range(n):
                s = base + off + i
                b = devb[np.random.default_rng([s, 7]).integers(devb.size)]
                a = beta_to_alpha(b)
                out.append((kind, n_sec, a, s, {"dev_beta": b, "twin": False}))
                if kind == "planet":
                    out.append(("noise", n_sec, a, s, {"dev_beta": b, "twin": True}))
    return out


def thresholds(cal):
    d = pd.read_csv(cal); thr = {}
    for tier in ("normal", "strict"):
        for f in FARS:
            q = [np.quantile(g.snr.values, 1 - f) for _, g in d[(d.cal_tier == tier) & (d.kind == "noise")].groupby("n_sec")]
            thr[(tier, f)] = np.ceil(max(q) * 10) / 10
    return thr


def summary(val, cal):
    from multisector_threshold_sim import _binom_ci
    c = pd.read_csv(cal); v = pd.read_csv(val); thr = thresholds(cal)
    L = [f"# Results: sims/two_tier_threshold_sim.py ({time.strftime('%Y-%m-%d')})", "",
         f"Edge alphas: normal {edge_alpha('normal'):.3f}, strict {edge_alpha('strict'):.3f}", "", "## Calibration"]
    for (tier, n_sec), g in c.groupby(["cal_tier", "n_sec"]):
        nz = g[g.kind == "noise"]; pl_ = g[g.kind == "planet"]
        emp = (nz.beta2h.quantile(.84) - nz.beta2h.quantile(.16)) / 2 / nz.beta2h.median()
        L.append(f"- {tier}, {n_sec} sectors: median beta2h {nz.beta2h.median():.2f}; sigma_beta formula "
                 f"{np.sqrt(1.36 / nz.M.median()):.1%} vs empirical {emp:.1%}; " +
                 "; ".join(f"FAR {f:.2%} quantile {np.quantile(nz.snr, 1 - f):.2f} (planets kept "
                           f"{np.mean((pl_.snr > np.quantile(nz.snr, 1 - f)) & pl_.correct_period):.0%})" for f in FARS))
    L += ["", "Tier thresholds (higher of 2 / 4 sectors, rounded up): " +
          ", ".join(f"{t} FAR {f:.2%} = {thr[(t, f)]:.1f}" for (t, f) in thr), "", "## Validation"]
    ok = True
    for n_sec, g in v.groupby("n_sec"):
        nz = g[(g.kind == "noise") & ~g.twin]; pl_ = g[g.kind == "planet"]; tw = g[(g.kind == "noise") & g.twin]
        sh = nz.tier.value_counts(normalize=True)
        chg = pl_.merge(tw[["seed", "tier"]], on="seed", suffixes=("", "_twin"))
        L.append(f"### {n_sec} sectors ({len(nz)} noise, {len(pl_)} planets)")
        L.append(f"- tier shares (noise stars): " + ", ".join(f"{t} {sh.get(t, 0):.1%}" for t in ("normal", "strict", "tier3")) +
                 f"; planets whose tier differs from their noise twin: {np.mean(chg.tier != chg.tier_twin):.1%}")
        for f in FARS:
            def claimed(df):
                t = df.tier.map({"normal": thr[("normal", f)], "strict": thr[("strict", f)]})
                return (df.tier != "tier3") & (df.snr > t)
            fa = claimed(nz); k = int(fa.sum()); lo, hi = [x / len(nz) for x in _binom_ci(k, len(nz))]
            t1 = lo <= f
            rates = {t: fa[nz.tier == t].mean() if (nz.tier == t).any() else 0.0 for t in ("normal", "strict")}
            t2 = all(r <= 2 * f for r in rates.values())
            rec2 = np.mean(claimed(pl_) & pl_.correct_period)
            single = np.quantile(nz.snr, 1 - k / len(nz)) if k else np.inf
            rec1 = np.mean((pl_.snr > single) & pl_.correct_period)
            t3 = rec2 > rec1
            list3 = np.mean((nz.tier == "tier3") & (nz.snr >= thr[("strict", f)]))
            ok &= t1 and t2 and t3
            L.append(f"- FAR {f:.2%}: T1 overall {k / len(nz):.2%} (95% {lo:.2%}–{hi:.2%}) → {'pass' if t1 else 'FAIL'}; "
                     f"T2 normal {rates['normal']:.2%}, strict {rates['strict']:.2%} → {'pass' if t2 else 'FAIL'}; "
                     f"T3 planets claimed {rec2:.1%} vs single threshold {single:.2f} {rec1:.1%} → {'pass' if t3 else 'FAIL'}; "
                     f"tier-3 noise listed {list3:.2%}")
    L += ["", f"**Validation: {'PASS' if ok else 'FAIL'}** (adoption is Julien's decision)"]
    return "\n".join(L)


def run(jobs, out, workers):
    done = set()
    if os.path.exists(out):
        d = pd.read_csv(out); done = set(zip(d.kind, d.n_sec, d.seed, d.get("twin", pd.Series([False] * len(d)))))
    todo = [j for j in jobs if (j[0], j[1], j[3], j[4].get("twin", False)) not in done]
    print(f"{len(todo)} stars to run ({len(done)} done)", flush=True)
    t0 = time.time()
    with Pool(workers) as pool:
        for n, r in enumerate(pool.imap_unordered(run_star, todo, chunksize=8), 1):
            pd.DataFrame([r]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
            if n % 1000 == 0:
                el = time.time() - t0
                print(f"{n}/{len(todo)} ~{el / n * (len(todo) - n) / 3600:.1f} h left", flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["calibrate", "validate", "summary"]); ap.add_argument("out")
    ap.add_argument("--cal"); ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    if a.mode == "calibrate":
        run(cal_jobs(), a.out, a.workers)
    elif a.mode == "validate":
        assert a.cal and os.path.exists(a.cal), "calibration first"
        print("thresholds:", thresholds(a.cal), flush=True)
        run(val_jobs(), a.out, a.workers)
    else:
        print(summary(a.out, a.cal))
