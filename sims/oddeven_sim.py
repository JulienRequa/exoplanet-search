"""Simulated planets / EBs: odd/even FAIL rate, pipeline.odd_even vs prototype odd_even_v3.
Run from the repo root: python sims/oddeven_sim.py   (results: sims/oddeven_sim_results_2026-11-24.md)"""
import sys, io, contextlib
sys.path.insert(0, ".")  # run from the repo root: python sims/oddeven_sim.py
import numpy as np, lightkurve as lk, pipeline as pl
from sims.oddeven_prototypes import odd_even_v3  # prototype; pipeline.odd_even is the same algorithm

# two TESS-like sectors (each 2 × 12.5 d orbits with a 1 d gap), 2-min cadence
cad = 2 / 1440
T = np.r_[np.arange(0, 12.5, cad), np.arange(13.5, 26, cad), np.arange(55, 67.5, cad), np.arange(68.5, 81, cad)]

def red(rng, n, amp, tau_min=60):
    phi = np.exp(-2 / tau_min); x = np.zeros(n); e = rng.normal(0, amp * np.sqrt(1 - phi**2), n)
    for i in range(1, n): x[i] = phi * x[i - 1] + e[i]
    return x

def one(seed, red_amp, odd_scale):
    rng = np.random.default_rng(seed)
    P, t0 = 3.3 + rng.uniform(0, 1), 1.0 + rng.uniform(0, 1)
    p = [t0, P, 0.1, 9.0, rng.uniform(0, 0.7), 0.4, 0.2]
    m = pl.transit_model(T, *p)
    if odd_scale != 1:                                    # EB-like: odd eclipses shallower
        odd = np.round((T - t0) / P) % 2 == 1
        m = np.where(odd, 1 - odd_scale * (1 - m), m)
    f = m + rng.normal(0, 1e-3, T.size) + (red(rng, T.size, red_amp) if red_amp else 0)
    lc = lk.LightCurve(time=T, flux=f)
    with contextlib.redirect_stdout(io.StringIO()):
        popt, *_ = pl.fit_transit(lc, P, t0, 0.12, 0.01)
        _, _, s_old = pl.odd_even(lc, popt)
    s_new = odd_even_v3(lc, popt)["sigma"]
    return s_old >= 3, (s_new is not None and s_new >= 3)

if __name__ == "__main__":
    for label, red_amp, odd_scale, N in [("planet, white noise", 0, 1, 60), ("planet, red noise (1 h, = white)", 1e-3, 1, 60),
                                         ("EB odd 15% shallower, white", 0, 0.85, 30), ("EB odd 15% shallower, red", 1e-3, 0.85, 30)]:
        res = np.array([one(s, red_amp, odd_scale) for s in range(N)])
        print(f"{label:34s} N={N}: FAIL rate current {res[:,0].mean():5.1%}   v3 {res[:,1].mean():5.1%}", flush=True)
