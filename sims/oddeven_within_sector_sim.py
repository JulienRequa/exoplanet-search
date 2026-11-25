"""Within-sector odd/even (PROTOTYPE, not in pipeline) vs pipeline odd/even v3.
Physics: in crowded fields the crowding (dilution) correction differs per sector, so the transit
depth can differ between sectors. If odd and even transits come unequally from the two sectors,
that sector difference looks like an odd/even difference. Comparing odd vs even INSIDE each
sector and combining the per-sector differences removes this.
Criteria: planets with sector depth differences are not falsely failed; EBs still caught.
Run from the repo root: python sims/oddeven_within_sector_sim.py"""
import sys, io, contextlib
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl
from sims.oddeven_sim import T, red

SECTOR = np.where(T < 40, 1, 2)                     # two sectors (see oddeven_sim.T)

def odd_even_ws(lc, popt, sector, n_null=300, seed=0):
    """Odd − even depth difference within each sector (v3 β and full-transit rule), combined
    across sectors with inverse-variance weights. Returns sigma (None if no sector has both)."""
    t = np.asarray(lc.time.value, float); f = np.asarray(lc.flux.value, float)
    t0, P = popt[0], popt[1]
    g = np.linspace(-0.5, 0.5, 20001); mg = pl.transit_model(t0 + g * P, *popt)
    T14 = np.mean(mg < 1) * P; core = 0.5 * np.mean(mg < 1 - 0.8 * (1 - mg.min())) * P
    rng = np.random.default_rng(seed); frac = T14 / P
    allowed = np.r_[np.linspace(2 * frac, 0.5 - 2 * frac, 400), np.linspace(0.5 + 2 * frac, 1 - 2 * frac, 400)]
    z = []
    for sh in rng.choice(allowed, n_null) * P:                      # same β as pipeline v3
        d = pl._single_depths(t, f, t0, P, T14, core, shift=sh)
        if len(d):
            z.extend((d[:, 1] - np.median(d[:, 1])) / d[:, 2])
    beta = max(1.0, float(np.std(z))) if len(z) > 30 else 1.0
    num = den = 0.0
    for s in np.unique(sector):
        m = sector == s
        d = pl._single_depths(t[m], f[m], t0, P, T14, core)
        if len(d) == 0 or len(np.unique(d[:, 0])) < 2:
            continue
        dep, err = {}, {}
        for p in (0, 1):
            gp = d[d[:, 0] == p]
            dep[p] = np.sum(gp[:, 1] * gp[:, 3]) / gp[:, 3].sum()
            err[p] = beta * np.sqrt(1 / np.sum(1 / gp[:, 2] ** 2))
        diff, e2 = dep[0] - dep[1], err[0] ** 2 + err[1] ** 2
        num += diff / e2; den += 1 / e2
    return None if den == 0 else abs(num / den) * np.sqrt(den)

def one(seed, case):
    rng = np.random.default_rng(seed)
    P, t0 = 3.3 + rng.uniform(0, 1), 1.0 + rng.uniform(0, 1)
    p = [t0, P, 0.1, 9.0, rng.uniform(0, 0.7), 0.4, 0.2]
    m = pl.transit_model(T, *p)
    n = np.round((T - t0) / P)
    if case.startswith("eb"):                                       # odd eclipses 15% shallower
        m = np.where(n % 2 == 1, 1 - 0.85 * (1 - m), m)
    if "crowd" in case:                                             # sector-2 depth 10–30% different
        scale = np.where(SECTOR == 2, 1 + rng.choice([-1, 1]) * rng.uniform(0.1, 0.3), 1.0)
        m = 1 - scale * (1 - m)
    keep = np.ones(T.size, bool)                                    # random missing transits (gaps)
    for k in np.unique(n):
        if rng.uniform() < 0.25:
            keep &= ~((n == k) & (np.abs(T - (t0 + k * P)) < 0.3))
    f = m + rng.normal(0, 1e-3, T.size) + red(rng, T.size, 5e-4)
    lc = lk.LightCurve(time=T[keep], flux=f[keep])
    with contextlib.redirect_stdout(io.StringIO()):
        popt, *_ = pl.fit_transit(lc, P, t0, 0.12, 0.01)
        s_v3 = pl.odd_even(lc, popt)[2]
    s_ws = odd_even_ws(lc, popt, SECTOR[keep])
    fail = lambda s: s is not None and s >= 3
    return fail(s_v3), fail(s_ws)

if __name__ == "__main__":
    for case, N in [("planet", 60), ("planet crowd", 60), ("eb", 30), ("eb crowd", 30)]:
        r = np.array([one(s, case) for s in range(N)])
        print(f"{case:13s} N={N}: FAIL rate — pipeline v3 {r[:,0].mean():5.1%}   within-sector {r[:,1].mean():5.1%}", flush=True)
