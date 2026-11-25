"""Period window on ORDINARY (quiet) stars: current max-power rule vs always-1×.
Criteria (Julien, 2026-11-24): injected planets recovered at the correct period (final fit
within 1%); EBs whose true period is 2× the catalog value still caught (odd/even or secondary).
One sector, white noise 1000–4000 ppm per 2-min point (faint stars), Rp/Rs 0.05–0.10.
Run from the repo root: python sims/period_window_ordinary_sim.py"""
import sys, io, contextlib
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl
from sims.period_window_sim import T, windows, caught

def make(seed, kind):
    rng = np.random.default_rng(seed)
    Pcat = rng.uniform(2, 9); t0 = rng.uniform(0.3, 1.5); a = 4.2 * Pcat ** (2 / 3)
    noise = rng.uniform(1e-3, 4e-3)
    if kind == "planet":
        m = pl.transit_model(T, t0, Pcat, rng.uniform(0.05, 0.10), a, rng.uniform(0, 0.6), 0.4, 0.2)
    else:
        ratio = {"eb60": 0.6, "eb90": 0.9}[kind]
        prim = pl.transit_model(T, t0, 2 * Pcat, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
        sec = pl.transit_model(T, t0 + Pcat, 2 * Pcat, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
        m = prim - ratio * (1 - sec)
    lc = lk.LightCurve(time=T, flux=m + rng.normal(0, noise, T.size), flux_err=np.full(T.size, noise))
    return lc, Pcat

def fitted_ok(c, w, lab, Pcat):
    pg = w[lab][1]
    try:
        popt, *_ = pl.fit_transit(c, pg.period_at_max_power, pg.transit_time_at_max_power,
                                  pg.duration_at_max_power, pg.depth_at_max_power)
        return abs(popt[1] - Pcat) / Pcat < 0.01
    except Exception:
        return False

if __name__ == "__main__":
    res = []
    for s in range(300, 400):
        lc, Pcat = make(s, "planet")
        with contextlib.redirect_stdout(io.StringIO()):
            c = pl.clean(lc); w = windows(c, Pcat)
            cur = max(w, key=lambda k: w[k][0])
            res.append((cur == "1×", fitted_ok(c, w, cur, Pcat), fitted_ok(c, w, "1×", Pcat)))
    r = np.array(res)
    print(f"PLANETS on quiet stars (100): window = 1× — current {r[:,0].mean():.0%}, always-1× 100%; "
          f"final period within 1% — current {r[:,1].mean():.0%}, always-1× {r[:,2].mean():.0%}", flush=True)
    for kind in ("eb60", "eb90"):
        cc = []
        for s in range(500, 530):
            lc, Pcat = make(s, kind)
            with contextlib.redirect_stdout(io.StringIO()):
                c = pl.clean(lc); w = windows(c, Pcat)
                cur = max(w, key=lambda k: w[k][0])
                cc.append((caught(c, w, cur), caught(c, w, "1×")))
        cc = np.array(cc)
        print(f"EB {kind} (30, true P = 2× catalog, quiet star): caught — current {cc[:,0].mean():.0%}, "
              f"always-1× {cc[:,1].mean():.0%}", flush=True)
