"""PROTOTYPE odd/even v2 (not in pipeline): (b) red-noise inflation β from the scatter of
single-transit depths *within* each parity group; (d) drop transits not fully covered."""
import numpy as np, pipeline as pl

def odd_even_v2(lc, popt):
    t = np.asarray(lc.time.value, float); f = np.asarray(lc.flux.value, float)
    t0, P = popt[0], popt[1]
    n = np.round((t - t0) / P)
    model = pl.transit_model(t, *popt)
    core = model < 1 - 0.8 * (1 - model.min())
    grid = np.linspace(-0.5, 0.5, 20001)
    T14 = np.mean(pl.transit_model(t0 + grid * P, *popt) < 1) * P       # days
    cad = np.median(np.diff(t))
    per = []                                       # (parity, depth, white err, n_pts)
    dropped = 0
    for k in np.unique(n[core]):
        tc = t0 + k * P
        w = np.abs(t - tc) <= T14 / 2
        tw = t[w]
        full = (w.sum() > 0 and tw.min() <= tc - T14 / 2 + 2 * cad and tw.max() >= tc + T14 / 2 - 2 * cad
                and (np.diff(tw).max() if w.sum() > 1 else np.inf) < 10 * cad)
        m = core & (n == k)
        if not full or m.sum() < 3:
            dropped += 1
            continue
        per.append((int(k % 2), 1 - f[m].mean(), f[m].std() / np.sqrt(m.sum()), m.sum()))
    per = np.array(per, dtype=float)
    out = {"dropped": dropped}
    if len(per) == 0 or len(np.unique(per[:, 0])) < 2:
        out.update(sigma=None, beta=None); return out
    chi2, dof, dep, err = 0.0, 0, {}, {}
    for p in (0, 1):
        g = per[per[:, 0] == p]
        wts = g[:, 3]
        dep[p] = np.sum(g[:, 1] * wts) / wts.sum()                 # = mean of all core points
        err[p] = np.sqrt(1 / np.sum(1 / g[:, 2] ** 2))             # white-noise error of the group
        if len(g) > 1:
            chi2 += np.sum(((g[:, 1] - dep[p]) / g[:, 2]) ** 2); dof += len(g) - 1
    beta = max(1.0, np.sqrt(chi2 / dof)) if dof > 0 else 1.0
    sigma = abs(dep[0] - dep[1]) / (beta * np.sqrt(err[0] ** 2 + err[1] ** 2))
    out.update(sigma=sigma, beta=beta, dof=dof, n_even=int((per[:, 0] == 0).sum()), n_odd=int((per[:, 0] == 1).sum()))
    return out


def odd_even_v3(lc, popt, n_null=300, seed=0):
    """(b') red-noise factor β from FAKE transits out of transit (signal can't leak in);
    (d) drop transits not fully covered. Group depths/white errors as before."""
    t = np.asarray(lc.time.value, float); f = np.asarray(lc.flux.value, float)
    t0, P = popt[0], popt[1]
    grid = np.linspace(-0.5, 0.5, 20001)
    mg = pl.transit_model(t0 + grid * P, *popt)
    T14 = np.mean(mg < 1) * P
    core_half = 0.5 * np.mean(mg < 1 - 0.8 * (1 - mg.min())) * P      # days
    cad = np.median(np.diff(t))

    def transit_depths(shift):
        """single-'transit' depths of core windows centred at t0 + shift + k·P"""
        n = np.round((t - t0 - shift) / P)
        dt = t - (t0 + shift + n * P)
        res = []
        for k in np.unique(n[np.abs(dt) <= core_half]):
            tc = t0 + shift + k * P
            w = np.abs(t - tc) <= T14 / 2
            tw = t[w]
            full = (w.sum() > 1 and tw.min() <= tc - T14 / 2 + 2 * cad
                    and tw.max() >= tc + T14 / 2 - 2 * cad and np.diff(tw).max() < 10 * cad)
            m = (n == k) & (np.abs(dt) <= core_half)
            if full and m.sum() >= 3:
                res.append((int(k % 2), 1 - f[m].mean(), f[m].std() / np.sqrt(m.sum()), m.sum()))
        return np.array(res, float)

    per = transit_depths(0.0)
    n_all = len(np.unique(np.round((t - t0) / P)[np.abs(t - t0 - np.round((t - t0) / P) * P) <= core_half]))
    out = {"dropped": n_all - len(per)}
    if len(per) == 0 or len(np.unique(per[:, 0])) < 2:
        out.update(sigma=None, beta=None); return out
    # null: fake transits away from the real transit and phase 0.5
    rng = np.random.default_rng(seed)
    frac = T14 / P
    ok = np.r_[np.linspace(2 * frac, 0.5 - 2 * frac, 400), np.linspace(0.5 + 2 * frac, 1 - 2 * frac, 400)]
    z = []
    for sh in rng.choice(ok, n_null) * P:
        d = transit_depths(sh)
        if len(d):
            z.extend((d[:, 1] - np.median(d[:, 1])) / d[:, 2])   # normalised single-depth residuals
    beta = max(1.0, float(np.std(z))) if len(z) > 30 else 1.0
    dep, err = {}, {}
    for p in (0, 1):
        g = per[per[:, 0] == p]
        dep[p] = np.sum(g[:, 1] * g[:, 3]) / g[:, 3].sum()
        err[p] = np.sqrt(1 / np.sum(1 / g[:, 2] ** 2))
    sigma = abs(dep[0] - dep[1]) / (beta * np.sqrt(err[0] ** 2 + err[1] ** 2))
    out.update(sigma=sigma, beta=beta, n_even=int((per[:, 0] == 0).sum()), n_odd=int((per[:, 0] == 1).sum()))
    return out
