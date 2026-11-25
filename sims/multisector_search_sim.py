"""Multi-sector search with a period hint (PROTOTYPE, not in pipeline) vs the current
single-sector search. Physics: a planet too faint to find in one sector (SNR ∝ √N_transits)
can be found when all sectors are stacked. Criteria: faint planets recovered at the correct
period more often; EBs at 2× catalog still caught; noise-only stars do NOT pass the detection
check more often (no new false passes); runtime per star.
Run from the repo root: python sims/multisector_search_sim.py"""
import sys, io, contextlib, time
sys.path.insert(0, ".")
import numpy as np, lightkurve as lk, pipeline as pl

cad = 2 / 1440
def sector(t0):
    return np.r_[np.arange(t0, t0 + 12.5, cad), np.arange(t0 + 13.5, t0 + 26, cad)]
STARTS = {2: [0.0, 700.0], 4: [0.0, 27.0, 700.0, 727.0]}

def make(seed, kind, n_sec):
    rng = np.random.default_rng(seed)
    P = rng.uniform(2, 9); t0 = rng.uniform(0.3, 1.5); a = 4.2 * P ** (2 / 3)
    noise = rng.uniform(3e-3, 5e-3)                       # faint star
    lcs = []
    for s0 in STARTS[n_sec]:
        T = sector(s0)
        if kind == "planet":
            m = pl.transit_model(T, t0, P, rng.uniform(0.045, 0.065), a, rng.uniform(0, 0.6), 0.4, 0.2)
        elif kind == "eb":                                 # true period 2× catalog, secondary 60%
            prim = pl.transit_model(T, t0, 2 * P, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
            sec = pl.transit_model(T, t0 + P, 2 * P, 0.15, 1.6 * a, 0.3, 0.4, 0.2)
            m = prim - 0.6 * (1 - sec)
        else:
            m = np.ones(T.size)                            # noise only
        lc = lk.LightCurve(time=T, flux=m + rng.normal(0, noise, T.size), flux_err=np.full(T.size, noise))
        lc.meta["SECTOR"] = len(lcs) + 1
        lcs.append(lc)
    return lcs, P

def snr(lc, popt):
    """Detection SNR as in vet_star."""
    t = np.asarray(lc.time.value, float); f = np.asarray(lc.flux.value, float)
    m = pl.transit_model(t, *popt); intr = m < 1
    return (1 - m.min()) / (np.std(f[~intr]) / np.sqrt(max(intr.sum(), 1))) if intr.any() else 0.0

def single(lcs, P):
    c0 = pl.clean(lcs[0]); Ph, T0, D, dep, _ = pl.search_near(c0, P)
    c1 = pl.clean(lcs[0], period=Ph, t0=T0, duration=D)
    popt, errs, *_ = pl.fit_transit(c1, Ph, T0, D, dep)
    if len(lcs) > 1:
        T0e, Pe, sT0, sP, _ = pl.per_sector_ephemeris(lcs, popt, errs, D, dep)
        cl = pl.clean_all(lcs, period=Pe, t0=T0e, duration=D)
        popt, *_ = pl.fit_transit(cl, Ph, T0, D, dep, start=np.r_[T0e, Pe, popt[2:]], start_errs=[sT0, sP])
        return popt, cl
    return popt, c1

def multi(lcs, P, tol=0.02):
    st = pl.clean_all(lcs)
    t = np.asarray(st.time.value, float); span = t.max() - t.min()
    dP = P * 0.05 / (3 * span)                             # grid fine enough for the full span
    n = int(min(200000, 2 * tol * P / dP))
    pg = st.to_periodogram(method="bls", period=np.linspace(P * (1 - tol), P * (1 + tol), n),
                           duration=[0.05, 0.08, 0.1, 0.15, 0.2])
    Pb, T0, D, dep = pg.period_at_max_power, pg.transit_time_at_max_power, pg.duration_at_max_power, pg.depth_at_max_power
    cl = pl.clean_all(lcs, period=Pb, t0=T0, duration=D)
    popt, *_ = pl.fit_transit(cl, Pb, T0, D, dep, start=[pl.plain(T0), pl.plain(Pb), np.sqrt(max(pl.plain(dep), 1e-5)),
                                                         4.2 * pl.plain(Pb) ** (2 / 3), 0.3, 0.3, 0.2],
                              start_errs=[0.004, 3 * dP])
    return popt, cl

def run(fn, lcs, P):
    t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            popt, cl = fn(lcs, P)
        return popt, cl, time.time() - t
    except Exception:
        return None, None, time.time() - t

if __name__ == "__main__":
    for n_sec in (2, 4):
        out = {"single": [], "multi": []}
        for s in range(40):
            lcs, P = make(s, "planet", n_sec)
            for name, fn in (("single", single), ("multi", multi)):
                popt, cl, dt = run(fn, lcs, P)
                ok = popt is not None and abs(popt[1] - P) / P < 0.01 and snr(cl, popt) > 7.1
                out[name].append((ok, dt))
        print(f"FAINT PLANETS, {n_sec} sectors (40): recovered at correct period AND SNR > 7.1 — "
              f"single {np.mean([o[0] for o in out['single']]):.0%}, multi {np.mean([o[0] for o in out['multi']]):.0%}; "
              f"runtime/star single {np.median([o[1] for o in out['single']]):.0f} s, multi {np.median([o[1] for o in out['multi']]):.0f} s", flush=True)
        eb = {"single": [], "multi": []}
        for s in range(100, 120):
            lcs, P = make(s, "eb", n_sec)
            for name, fn in (("single", single), ("multi", multi)):
                popt, cl, _ = run(fn, lcs, P)
                c = False
                if popt is not None:
                    with contextlib.redirect_stdout(io.StringIO()):
                        oe = pl.odd_even(cl, popt)[2]; _, _, ss, sr = pl.secondary_eclipse(cl, popt)
                    c = (oe is not None and oe >= 3) or (ss > 3 and sr > 0.1)
                eb[name].append(c)
        print(f"EBs at 2× catalog, {n_sec} sectors (20): caught — single {np.mean(eb['single']):.0%}, multi {np.mean(eb['multi']):.0%}", flush=True)
        nz = {"single": [], "multi": []}
        for s in range(200, 240):
            lcs, P = make(s, "noise", n_sec)
            for name, fn in (("single", single), ("multi", multi)):
                popt, cl, _ = run(fn, lcs, P)
                nz[name].append(popt is not None and snr(cl, popt) > 7.1)
        print(f"NOISE ONLY, {n_sec} sectors (40): falsely pass the detection check — single {np.mean(nz['single']):.0%}, "
              f"multi {np.mean(nz['multi']):.0%}", flush=True)
