"""Ephemeris: pipeline.per_sector_ephemeris (timing fix) vs the old refinement.
Criteria: final period correct (|ΔP| within 3σ and sector-2 transits hit); EBs with true
period = 2× catalog still caught. Run from the repo root: python sims/ephemeris_proto_sim.py"""
import sys, io, contextlib
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl
from sims.oddeven_sim import red

cad = 2 / 1440
def sector(t0):
    return np.r_[np.arange(t0, t0 + 12.5, cad), np.arange(t0 + 13.5, t0 + 26, cad)]
GAP = 1100.0

def lcs_for(flux_fn, rng, red_amp):
    out = []
    for s0, sec in [(0.0, 1), (GAP, 2)]:
        T = sector(s0)
        f = flux_fn(T) + rng.normal(0, 1e-3, T.size) + (red(rng, T.size, red_amp) if red_amp else 0)
        lc = lk.LightCurve(time=T, flux=f, flux_err=np.full(T.size, 1e-3)); lc.meta["SECTOR"] = sec
        out.append(lc)
    return out

def analyse(lcs, hint, use_proto):
    """The pipeline's own steps (minus download): 1-sector search + fit, then refinement."""
    with contextlib.redirect_stdout(io.StringIO()):
        c0 = pl.clean(lcs[0])
        P, T0, D, dep, _ = pl.search_near(c0, hint)
        c1 = pl.clean(lcs[0], period=P, t0=T0, duration=D)
        popt, errs, *_ = pl.fit_transit(c1, P, T0, D, dep)
        if use_proto:
            T0e, Pe, sT0, sP, _ = pl.per_sector_ephemeris(lcs, popt, errs, D, dep)
            start, serr = np.r_[T0e, Pe, popt[2:]], [sT0, sP]
        else:
            start, serr, T0e, Pe = popt, errs, popt[0], popt[1]
        cl = pl.clean_all(lcs, period=Pe, t0=T0e, duration=D)
        pf, ef, *_ = pl.fit_transit(cl, P, T0, D, dep, start=start, start_errs=serr)
        _, _, oe = pl.odd_even(cl, pf)
        _, _, ss, sr = pl.secondary_eclipse(cl, pf)
    return pf, ef, (oe is not None and oe >= 3) or (ss > 3 and sr > 0.1)

if __name__ == "__main__":
    print("PLANETS (hard: P 6–9 d, Rp/Rs 0.067, red = white noise, 1100 d gap), N=30")
    ok = {False: [], True: []}
    for s in range(30):
        rng = np.random.default_rng(s)
        P = rng.uniform(6, 9); t0 = rng.uniform(0.5, 2)
        p = [t0, P, 0.067, 4.2 * P ** (2 / 3) * rng.uniform(0.9, 1.1), rng.uniform(0, 0.7), 0.4, 0.2]
        lcs = lcs_for(lambda T: pl.transit_model(T, *p), rng, 1e-3)
        g = np.linspace(-0.5, 0.5, 20001); half = np.mean(pl.transit_model(t0 + g * P, *p) < 1) * P / 2
        for proto in (False, True):
            pf, ef, _ = analyse(lcs, P, proto)
            k = np.round((GAP + 13 - pf[0]) / pf[1]); kt = np.round((GAP + 13 - t0) / P)
            hit = abs((pf[0] + k * pf[1]) - (t0 + kt * P)) < half
            within = abs(pf[1] - P) < 3 * ef[1]
            ok[proto].append((hit and within, abs(pf[1] - P)))
    for proto, lab in [(False, "current"), (True, "prototype")]:
        a = np.array(ok[proto])
        print(f"  {lab:9s}: correct period (sector-2 hit AND |ΔP| < 3σ) {a[:,0].mean():.0%};"
              f" median |ΔP| {np.median(a[:,1])*86400:.1f} s")
    print("EBs (true P = 2× catalog, secondary 60% / 90% of primary, white noise, 1100 d gap), N=20 each")
    for ratio in (0.6, 0.9):
        c = {False: 0, True: 0}
        for s in range(100, 120):
            rng = np.random.default_rng(s)
            Pc = rng.uniform(2, 9); t0 = rng.uniform(0.3, 1.5); a = 1.6 * 4.2 * Pc ** (2 / 3)
            flux = lambda T: (pl.transit_model(T, t0, 2 * Pc, 0.15, a, 0.3, 0.4, 0.2)
                              - ratio * (1 - pl.transit_model(T, t0 + Pc, 2 * Pc, 0.15, a, 0.3, 0.4, 0.2)))
            lcs = lcs_for(flux, rng, 0)
            for proto in (False, True):
                try:
                    c[proto] += analyse(lcs, Pc, proto)[2]
                except Exception:
                    pass
        print(f"  secondary {ratio:.0%}: caught — current {c[False]/20:.0%}, prototype {c[True]/20:.0%}")
