"""Fast tests (no internet except the TIC test, which is skipped offline).
Run: pytest -q
"""
import numpy as np
import lightkurve as lk
import pytest
import pipeline as pl

T = np.arange(0, 27, 2 / 1440)            # one TESS sector at 2-min cadence


def fake_lc(params, noise, seed=0):
    rng = np.random.default_rng(seed)
    f = pl.transit_model(T, *params) + rng.normal(0, noise, T.size)
    return lk.LightCurve(time=T, flux=f, flux_err=np.full(T.size, noise))


def fit(lc, ld_prior=None):
    P, T0, D, dep = pl.search(lc)
    popt, errs, inc, ie, flags, pcov = pl.fit_transit(lc, P, T0, D, dep, return_pcov=True,
                                                      ld_prior=ld_prior)
    r = {n: (v, e) for n, v, e in zip(pl.PARAM_NAMES, popt, errs)}
    r["cov_rp_b"], r["flags"], r["lc"], r["ld_prior"] = pcov[2, 4], flags, lc, ld_prior
    return r


def test_injected_planet_recovered():
    r = fit(fake_lc([1.3, 3.0, 0.1, 9.0, 0.3, 0.4, 0.2], 1e-3))
    assert abs(r["period"][0] - 3.0) < 1e-3
    assert abs(r["Rp/Rs"][0] - 0.1) < 3 * r["Rp/Rs"][1] + 1e-3


def test_deep_eclipse_does_not_crash():      # bug 1: guess sqrt(depth) > old 0.30 cap
    r = fit(fake_lc([1.3, 3.0, 0.5, 9.0, 0.2, 0.4, 0.2], 1e-3))
    assert 0.4 < r["Rp/Rs"][0] < 0.6


@pytest.mark.parametrize("hint, choice", [(3.0, "1×"), (6.0, "½×"), (1.5, "2×")])
def test_period_hint_picks_right_multiple(hint, choice):   # bug 3
    lc = fake_lc([1.3, 3.0, 0.1, 9.0, 0.3, 0.4, 0.2], 1e-3)
    P, T0, D, dep, ch = pl.search_near(lc, hint)
    assert abs(float(P.value) - 3.0) < 0.01 and ch == choice


def test_star_info_header_used_and_nan_is_missing(monkeypatch):   # bug 2
    lc = lk.LightCurve(time=np.arange(5.0), flux=np.ones(5))
    lc.meta.update(TICID=None, TEFF=5800.0, LOGG=4.4, RADIUS=1.0)
    info, src = pl.get_star_info(lc)
    assert info["RADIUS"] == 1.0 and set(src.values()) == {"header"}
    lc.meta.update(RADIUS=np.nan)            # TICID=None → lookup fails → stays missing
    info, src = pl.get_star_info(lc)
    assert info["RADIUS"] is None and src["RADIUS"] == "missing"


@pytest.mark.parametrize("prior", [None, (0.4, 0.2, 0.1)])
def test_normal_transit_not_flagged_grazing(prior):
    r = fit(fake_lc([1.3, 3.0, 0.15, 9.0, 0.3, 0.4, 0.2], 1e-3, seed=1), prior)
    g, sigma = pl.grazing_check(r)
    assert sigma < 3


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_grazing_transit_flagged(seed):
    """b + Rp/Rs = 1.17, 200 ppm, LD prior off by 0.1 (a realistic table error).
    Profile method + prior: flagged in all 3 seeds in the 2026-11-24 simulations."""
    r = fit(fake_lc([1.3, 3.0, 0.15, 9.0, 1.02, 0.4, 0.2], 2e-4, seed=seed), (0.5, 0.1, 0.1))
    g, sigma = pl.grazing_check(r)
    assert sigma > 3


def test_limb_darkening_always_physical():      # Kipping q1, q2: u1 + u2 ≤ 1, u1 ≥ 0, u1 + 2u2 ≥ 0
    r = fit(fake_lc([1.3, 3.0, 0.5, 9.0, 0.2, 0.4, 0.2], 1e-3))
    u1, u2 = r["u1"][0], r["u2"][0]
    assert u1 >= -1e-9 and u1 + u2 <= 1 + 1e-9 and u1 + 2 * u2 >= -1e-9


def test_q_u_roundtrip():
    for u1, u2 in [(0.4, 0.2), (0.1, 0.5), (0.9, 0.05)]:
        assert np.allclose(pl.q_to_u(*pl.u_to_q(u1, u2)), (u1, u2))


def test_ld_prior_pulls_fit():
    lc = fake_lc([1.3, 3.0, 0.1, 9.0, 0.3, 0.4, 0.2], 3e-3)
    P, T0, D, dep = pl.search(lc)
    popt, *_ = pl.fit_transit(lc, P, T0, D, dep, ld_prior=(0.6, 0.1, 0.02))
    assert abs(popt[5] - 0.6) < 0.06 and abs(popt[6] - 0.1) < 0.06


def test_ld_prior_lookup():
    prior, src = pl.ld_prior_for({"TEFF": 5750.0, "LOGG": 4.5})     # exact grid point
    assert prior is not None and abs(prior[0] - 0.3075) < 1e-3 and abs(prior[1] - 0.2743) < 1e-3
    assert pl.ld_prior_for({"TEFF": 6000.0, "LOGG": 8.0})[0] is None   # white dwarf: off table
    assert pl.ld_prior_for({"TEFF": None, "LOGG": 4.4})[0] is None


def _two_sector_lc(seed, odd_scale=1.0, red_amp=0.0):
    from sims.oddeven_sim import T, red
    rng = np.random.default_rng(seed)
    p = [1.5, 3.7, 0.1, 9.0, 0.3, 0.4, 0.2]
    m = pl.transit_model(T, *p)
    if odd_scale != 1:
        m = np.where(np.round((T - 1.5) / 3.7) % 2 == 1, 1 - odd_scale * (1 - m), m)
    f = m + rng.normal(0, 1e-3, T.size) + (red(rng, T.size, red_amp) if red_amp else 0)
    return lk.LightCurve(time=T, flux=f), p


def test_odd_even_planet_in_red_noise_passes():      # v3: old rule failed ~57% of these
    lc, p = _two_sector_lc(1, red_amp=1e-3)
    assert pl.odd_even(lc, p)[2] < 3


def test_odd_even_catches_eb():
    lc, p = _two_sector_lc(2, odd_scale=0.85)
    assert pl.odd_even(lc, p)[2] > 3


def test_download_retries_network_errors_only(monkeypatch):
    calls = {"n": 0}

    class FakeSearch(list):
        def __getitem__(self, s):
            return self
        def download_all(self):
            return ["lc1", "lc2"]

    def flaky(*a, **k):
        calls["n"] += 1
        if calls["n"] == 1:
            raise TimeoutError("Timeout limit of 600 exceeded.")
        return FakeSearch([1, 2])

    monkeypatch.setattr(pl.lk, "search_lightcurve", flaky)
    monkeypatch.setattr(pl.time, "sleep", lambda s: None)
    assert pl.download("X", 2) == ["lc1", "lc2"] and calls["n"] == 2     # retried once

    calls["n"] = 0
    monkeypatch.setattr(pl.lk, "search_lightcurve", lambda *a, **k: (calls.__setitem__("n", calls["n"] + 1), FakeSearch())[1])
    with pytest.raises(ValueError):                                       # no data: not retried
        pl.download("Y", 2)
    assert calls["n"] == 1
