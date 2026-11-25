"""Period-window rule for period_hint (PROTOTYPE, not in pipeline).
Current rule: highest BLS power among 1×, ½×, 2× windows. Proposed: keep 1× unless
power(alt)/power(1×) > R*, with R* = 99.7th percentile of that ratio for simulated TRUE-1×
planets on spotted stars (calibration seeds 0–99, never real targets).
Tests on fresh seeds: planets on spotted stars; EBs whose true period is 2× the catalog one.
Run from the repo root: python sims/period_window_sim.py"""
import sys, io, contextlib
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl

cad = 2 / 1440
T = np.r_[np.arange(0, 12.5, cad), np.arange(13.5, 26, cad)]          # one sector

def spots(rng):
    """Quasi-periodic spot modulation: rotation 0.5–5 d, 0.5–3% amplitude, evolving."""
    Prot, A = rng.uniform(0.5, 5), rng.uniform(0.005, 0.03)
    env = 1 + 0.3 * np.sin(2 * np.pi * T / rng.uniform(5, 20) + rng.uniform(0, 6.3))
    s = sum(A / h * env * np.sin(2 * np.pi * h * T / Prot + rng.uniform(0, 6.3)) for h in (1, 2, 3))
    return s, Prot, A

def make(seed, kind):
    rng = np.random.default_rng(seed)
    Pcat = rng.uniform(2, 9); t0 = rng.uniform(0.3, 1.5)
    a = 4.2 * Pcat ** (2 / 3)
    if kind == "planet":
        m = pl.transit_model(T, t0, Pcat, 0.07, a, rng.uniform(0, 0.6), 0.4, 0.2)
    else:                                   # EB: true period 2×Pcat, secondary = ratio × primary
        ratio = {"eb60": 0.6, "eb90": 0.9}[kind]
        prim = pl.transit_model(T, t0, 2 * Pcat, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
        sec = pl.transit_model(T, t0 + Pcat, 2 * Pcat, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
        m = prim - ratio * (1 - sec)
    s, Prot, A = spots(rng)
    f = (m + s) + rng.normal(0, 1e-3, T.size)
    lc = lk.LightCurve(time=T, flux=f, flux_err=np.full(T.size, 1e-3))
    return lc, Pcat, t0

def windows(lc, hint, tol=0.02, n_periods=2000):
    t = np.asarray(lc.time.value, float); half = (t.max() - t.min()) / 2
    out = {}
    for fct, lab in [(1, "1×"), (0.5, "½×"), (2, "2×")]:
        P = fct * hint
        if P * (1 + tol) > half:
            continue
        pg = lc.to_periodogram(method="bls", period=np.linspace(P * (1 - tol), P * (1 + tol), n_periods),
                               duration=[0.05, 0.08, 0.1, 0.15, 0.2])
        out[lab] = (pl.plain(pg.max_power), pg)
    return out

def choose(w, Rstar):
    """current: max power. proposed: 1× unless an alternative's power ratio > R*."""
    cur = max(w, key=lambda k: w[k][0])
    alts = {k: v[0] / w["1×"][0] for k, v in w.items() if k != "1×"}
    best_alt = max(alts, key=alts.get) if alts else None
    prop = best_alt if (best_alt and alts[best_alt] > Rstar) else "1×"
    return cur, prop, (max(alts.values()) if alts else 0)

def caught(lc, w, lab):
    pg = w[lab][1]
    with contextlib.redirect_stdout(io.StringIO()):
        popt, *_ = pl.fit_transit(lc, pg.period_at_max_power, pg.transit_time_at_max_power,
                                  pg.duration_at_max_power, pg.depth_at_max_power)
        _, _, oe = pl.odd_even(lc, popt)
        _, _, ss, sr = pl.secondary_eclipse(lc, popt)
    return (oe is not None and oe >= 3) or (ss > 3 and sr > 0.1)

if __name__ == "__main__":
    def prepared(seed, kind):
        lc, Pcat, t0 = make(seed, kind)
        with contextlib.redirect_stdout(io.StringIO()):
            c = pl.clean(lc)                                   # the pipeline's own cleaning
        return c, Pcat, windows(c, Pcat)

    # 1. calibration: ratio distribution for TRUE-1× planets on spotted stars
    ratios = []
    for s in range(100):
        c, Pcat, w = prepared(s, "planet")
        ratios.append(choose(w, np.inf)[2])
    ratios = np.array(ratios)
    Rstar = float(np.percentile(ratios, 99.7))
    print(f"CALIBRATION (100 planets, spotted): ratio alt/1× median {np.median(ratios):.2f}, "
          f"max {ratios.max():.2f}, 99.7th pct R* = {Rstar:.2f}  "
          f"(genuine ½×/2× alternative expected ≈ 2)", flush=True)

    # 2. fresh planets
    res = []
    for s in range(100, 160):
        c, Pcat, w = prepared(s, "planet")
        cur, prop, _ = choose(w, Rstar)
        res.append((cur == "1×", prop == "1×"))
    r = np.array(res)
    print(f"FRESH PLANETS (60, spotted): keep 1× — current {r[:,0].mean():.0%}, proposed {r[:,1].mean():.0%}", flush=True)

    # 3. EBs with true period 2× catalog
    for kind in ("eb60", "eb90"):
        res = []
        for s in range(200, 230):
            c, Pcat, w = prepared(s, kind)
            cur, prop, _ = choose(w, Rstar)
            res.append((caught(c, w, cur), caught(c, w, prop), cur, prop))
        cc = np.mean([x[0] for x in res]); cp = np.mean([x[1] for x in res])
        picks = {k: sum(1 for x in res if x[3] == k) for k in ("1×", "½×", "2×")}
        print(f"EB {kind} (30, true P = 2× catalog): caught — current {cc:.0%}, proposed {cp:.0%}; "
              f"proposed picks {picks}", flush=True)
