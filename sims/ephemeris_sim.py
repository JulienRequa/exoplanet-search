"""Ephemeris drift across long gaps: current multi-sector refinement vs per-sector transit
times + straight-line ephemeris (PROTOTYPE, not in pipeline).
Run from the repo root: python sims/ephemeris_sim.py"""
import sys, io, contextlib
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl

cad = 2 / 1440
def sector(t_start):
    return np.r_[np.arange(t_start, t_start + 12.5, cad), np.arange(t_start + 13.5, t_start + 26, cad)]

from sims.oddeven_sim import red

def make(seed, gap, Pmin=2, Pmax=9, k=0.08, red_amp=0.0):
    rng = np.random.default_rng(seed)
    P = rng.uniform(Pmin, Pmax); t0 = rng.uniform(0.5, 2)
    p = [t0, P, k, 4.2 * P ** (2 / 3) * rng.uniform(0.9, 1.1), rng.uniform(0, 0.7), 0.4, 0.2]
    T1, T2 = sector(0.0), sector(gap)
    mk = lambda T: lk.LightCurve(time=T, flux=pl.transit_model(T, *p) + rng.normal(0, 1e-3, T.size)
                                 + (red(rng, T.size, red_amp) if red_amp else 0))
    return p, mk(T1), mk(T2)

def per_sector_ephemeris(lc1, lc2, popt1, D, dep):
    """Transit time in each sector, then one straight line through them."""
    P1 = popt1[1]
    P_, T0_, D_, dep_, _ = pl.search_near(lc2, P1, tol=0.02)            # find sector-2 transits
    popt2, *_ = pl.fit_transit(lc2, P_, T0_, D, dep_)
    t1, t2 = popt1[0], popt2[0]
    N = np.round((t2 - t1) / P1)
    return t1, (t2 - t1) / N

def miss_hours(p, t0, P, t_check):
    """|predicted − true| transit time near t_check, in hours"""
    k = np.round((t_check - t0) / P)
    k_true = np.round((t_check - p[0]) / p[1])
    return abs((t0 + k * P) - (p[0] + k_true * p[1])) * 24

def T14_hours(p):
    g = np.linspace(-0.5, 0.5, 20001)
    return np.mean(pl.transit_model(p[0] + g * p[1], *p) < 1) * p[1] * 24

if __name__ == "__main__":
    # "easy": P 2–9 d, Rp/Rs 0.08, white noise, 700 d gap
    # "hard": TOI-2207-like — P 6–9 d (2–4 transits/sector), Rp/Rs 0.067, red noise = white, 1100 d gap
    case = sys.argv[1] if len(sys.argv) > 1 else "easy"
    cfg = {"easy": dict(gap=700.0, Pmin=2, Pmax=9, k=0.08, red_amp=0.0),
           "hard": dict(gap=1100.0, Pmin=6, Pmax=9, k=0.067, red_amp=1e-3)}[case]
    gap, N = cfg.pop("gap"), 30
    print("case:", case, cfg, "gap", gap)
    rows = []
    for s in range(N):
        p, lc1, lc2 = make(s, gap, **cfg)
        both = lc1.append(lc2)
        with contextlib.redirect_stdout(io.StringIO()):
            P_, T0_, D_, dep_ = pl.search(lc1)
            popt1, errs1, *_ = pl.fit_transit(lc1, P_, T0_, D_, dep_)
            cur, *_ = pl.fit_transit(both, P_, T0_, D_, dep_, start=popt1, start_errs=errs1)
            t0n, Pn = per_sector_ephemeris(lc1, lc2, popt1, pl.plain(D_), pl.plain(dep_))
            new, *_ = pl.fit_transit(both, Pn, t0n, D_, dep_, start=np.r_[t0n, Pn, popt1[2:]],
                                     start_errs=[0.002, 1e-5])
        half = T14_hours(p) / 2
        mc, mn = miss_hours(p, cur[0], cur[1], gap + 13), miss_hours(p, new[0], new[1], gap + 13)
        drift1 = miss_hours(p, popt1[0], popt1[1], gap + 13)
        rows.append((mc > half, mn > half, drift1 > half, drift1))
        print(f"seed {s:2d} P={p[1]:.3f}  1-sector drift {drift1:6.2f} h | current miss {mc:6.2f} h | "
              f"per-sector {mn:6.2f} h | half-duration {half:.2f} h", flush=True)
    r = np.array(rows)
    print(f"\nGap {gap:.0f} d, N={N}: sector-2 transits LOST — current {r[:,0].mean():.0%}, "
          f"per-sector ephemeris {r[:,1].mean():.0%}  (1-sector ephemeris alone drifts > ½ duration in "
          f"{r[:,2].mean():.0%}; median drift {np.median(r[:,3]):.1f} h)")
