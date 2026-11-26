"""Diagnostic report cards for flagged hunt candidates (strong fakes). Diagnostic only:
reads the frozen v2 hunt results and the cached light curves; changes nothing in v2.

    python scripts/report_cards.py LIST.csv OUTDIR      (LIST has a TOI column)

Per target: a PNG card + one row in OUTDIR/report_cards.csv. Questions answered:
 (a) Could stellar variability or a phase curve explain the flag?
     - secondary eclipse: depth at phase 0.5 measured against its LOCAL baseline, and its
       significance against a red-noise null (the same window placed at 200 random phases away
       from transit and phase 0.5, seed 0), in the v2-cleaned light curve.
     - odd/even: each full transit's depth against its own local baseline; odd − even
       significance against the same kind of null.
     - Lomb–Scargle periodogram of the un-flattened light curve (transits masked), with P, P/2
       and 2P marked; phase-curve amplitudes (cos 2πφ and cos 4πφ terms).
 (b) Could a nearby eclipsing binary explain it?
     - TIC stars within 2′ (Tmag, separation, eclipse depth a star would need to make the
       observed dip if fully blended: depth × 10^(0.4 ΔTmag)).
     - Known EBs / variables within 2′: VSX (B/vsx), Gaia DR3 eclipsing binaries (I/358/veb),
       TESS EB catalog (J/ApJS/258/16), via the VizieR mirror vizier.cfa.harvard.edu.
Assessment rules (written before running on any target; diagnostic thresholds, not v2):
 - A secondary flag survives (a) if local-baseline depth > 0 at >= 3σ against the red-noise null.
 - An odd/even flag survives (a) if the local-baseline odd − even difference is >= 3σ against
   the red-noise null.
 - Centroid / neighbour-depth / size flags are not affected by (a).
 - (b) fails if a known EB/variable (other than the target) lies within 2′, is bright enough
   (needed depth <= 100%), and its period is within 1% of P, 2P or P/2 ("known NEB"), or its
   period is unknown/different but it could still contribute ("possible contamination").
 - The target itself listed as an EB/variable is reported (it means the flag is already public).
Sort: survives both / variability / contamination / unclear.
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.timeseries import LombScargle
import pipeline as pl

HUNT = "results/hunt_v2_2026-11-25.csv"
STRONG = "results/hunt_v2_compare_2026-11-26/strong_fakes_followup.csv"
RNG_SEED, N_NULL, RADIUS = 0, 200, 2 * u.arcmin


def t14(P, k, a, b, inc_deg):
    s = np.sqrt(max((1 + k) ** 2 - b ** 2, 1e-6)) / (a * np.sin(np.radians(inc_deg)))
    return P / np.pi * np.arcsin(min(s, 1.0))


def window_depth(t, f, tc_list, half, lo, hi):
    """Mean depth inside |t-tc|<half vs local baseline lo<|t-tc|<hi, per centre."""
    out = []
    for tc in tc_list:
        d = np.abs(t - tc)
        w, b = d < half, (d > lo) & (d < hi)
        if w.sum() >= 3 and b.sum() >= 6:
            out.append(np.mean(f[b]) - np.mean(f[w]))
    return np.array(out)


def centres(t, t0, P, shift):
    n = np.arange(np.floor((t.min() - t0) / P) - 1, np.ceil((t.max() - t0) / P) + 2)
    return t0 + shift + n * P


def secondary_local(t, f, t0, P, D, rng):
    half, lo, hi = 0.4 * D, 0.75 * D, 2.0 * D
    obs = window_depth(t, f, centres(t, t0, P, 0.5 * P), half, lo, hi)
    if len(obs) < 2:
        return np.nan, np.nan, np.nan
    m = obs.mean()
    null = []
    phases = rng.uniform(0.1, 0.9, N_NULL * 3)
    phases = phases[np.abs(phases - 0.5) > 0.1][:N_NULL]
    for ph in phases:
        x = window_depth(t, f, centres(t, t0, P, ph * P), half, lo, hi)
        if len(x) >= 2:
            null.append(x.mean())
    sd = np.std(null)
    return m, sd, m / sd if sd > 0 else np.nan


def oddeven_local(t, f, t0, P, D, rng):
    half, lo, hi = 0.25 * D, 0.75 * D, 2.0 * D
    def diff(shift):
        cs = centres(t, t0, P, shift)
        k = np.round((cs - t0 - shift) / P).astype(int)
        d = window_depth(t, f, cs, half, lo, hi)
        ok = np.array([len(window_depth(t, f, [c], half, lo, hi)) == 1 for c in cs])
        kk = k[ok]
        if len(d) < 4 or len(set(kk % 2)) < 2:
            return np.nan, np.nan, np.nan
        return d[kk % 2 == 1].mean() - d[kk % 2 == 0].mean(), d[kk % 2 == 1].mean(), d[kk % 2 == 0].mean()
    obs, dodd, deven = diff(0.0)
    if np.isnan(obs):
        return np.nan, np.nan, np.nan, np.nan, np.nan
    null = []
    for ph in rng.uniform(0.15, 0.85, N_NULL):
        if abs(ph - 0.5) < 0.1:
            continue
        x = diff(ph * P)[0]
        if np.isfinite(x):
            null.append(x)
    sd = np.std(null)
    return obs, sd, abs(obs) / sd if sd > 0 else np.nan, dodd, deven


def known_variables(coord, P, tic):
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**", "+_r"], row_limit=200, vizier_server="vizier.cfa.harvard.edu")
    out, err = [], []
    for cat in ["B/vsx/vsx", "I/358/veb", "J/ApJS/258/16"]:
        try:
            res = v.query_region(coord, radius=RADIUS, catalog=cat)
        except Exception as e:
            err.append(f"{cat}: {type(e).__name__}")
            continue
        for tab in res:
            for r in tab:
                sep = float(r["_r"]) * 60 if "_r" in tab.colnames else np.nan  # arcmin → arcsec
                per, typ, name = np.nan, "", ""
                if cat == "B/vsx/vsx":
                    per = float(r["Period"]) if "Period" in tab.colnames and r["Period"] else np.nan
                    typ, name = str(r["Type"]), str(r["Name"])
                elif cat == "I/358/veb":
                    fq = float(r["Freq"]) if r["Freq"] else np.nan
                    per, typ, name = (1 / fq if fq else np.nan), "EB (Gaia DR3)", f"Gaia {r['Source']}"
                else:
                    for c in ("Per", "period", "P"):
                        if c in tab.colnames:
                            per = float(r[c]); break
                    typ = "EB (TESS EB cat.)"
                    name = f"TIC {r['TIC']}" if "TIC" in tab.colnames else tab.meta.get("name", "")
                out.append(dict(cat=cat, name=name, type=typ, period=per, sep_arcsec=sep))
    return out, err


def period_match(p, P):
    return np.isfinite(p) and any(abs(p / (m * P) - 1) < 0.01 for m in (0.5, 1, 2))


def card(toi, hunt, strong, outdir, rng):
    h = hunt[hunt.TOI == toi].iloc[0]
    s = strong[strong.TOI == toi].iloc[0]
    tic = int(str(h.star).split()[-1])
    t0, P, k, a, b, inc = h.t0, h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc
    D = t14(P, k, a, b, inc)
    popt = [t0, P, k, a, b, h.u1, h.u2]
    depth = 1 - pl.transit_model(np.array([t0]), *popt)[0]

    lcs = pl.download(f"TIC {tic}", 2)
    flat = pl.clean_all(lcs, period=P, t0=t0, duration=D)
    t = np.asarray(flat.time.value, float); f = np.asarray(flat.flux.value, float)
    ok = np.isfinite(f); t, f = t[ok], f[ok]

    raw_t, raw_f = [], []
    for lc in lcs:
        lc = lc.remove_nans()
        x = np.asarray(lc.flux.value, float); x = x / np.nanmedian(x)
        raw_t.append(np.asarray(lc.time.value, float)); raw_f.append(x)
    rt, rf = np.concatenate(raw_t), np.concatenate(raw_f)
    rph = ((rt - t0) / P + 0.5) % 1 - 0.5
    oot = np.abs(rph) > 0.75 * D / P
    good = oot & np.isfinite(rf)
    clip = np.abs(rf[good] - np.median(rf[good])) < 5 * np.std(rf[good])
    lt, lf = rt[good][clip], rf[good][clip]

    sec, sec_sd, sec_sig = secondary_local(t, f, t0, P, D, rng)
    oe, oe_sd, oe_sig, d_odd, d_even = oddeven_local(t, f, t0, P, D, rng)

    freq = np.linspace(1 / 30, 1 / 0.1, 30000)
    ls = LombScargle(lt, lf)
    pw = ls.power(freq)
    ipk = np.argmax(pw); p_peak = 1 / freq[ipk]
    amp_peak = np.sqrt(pw[ipk] * 2) * np.std(lf)  # rough semi-amplitude
    near = [m for m in (0.5, 1, 2) if abs(p_peak / (m * P) - 1) < 0.02]
    phi = ((lt - t0) / P) % 1
    A = np.vstack([np.ones_like(phi), np.cos(2 * np.pi * phi), np.sin(2 * np.pi * phi),
                   np.cos(4 * np.pi * phi), np.sin(4 * np.pi * phi)]).T
    coef = np.linalg.lstsq(A, lf, rcond=None)[0]
    amp1, amp2 = np.hypot(coef[1], coef[2]), np.hypot(coef[3], coef[4])

    coord = SkyCoord(float(lcs[0].meta["RA_OBJ"]), float(lcs[0].meta["DEC_OBJ"]), unit="deg")
    tmag_t = float(h.TESSMAG) if pd.notna(h.TESSMAG) else np.nan
    nb = []
    try:
        from astroquery.mast import Catalogs
        tt = Catalogs.query_region(coord, radius=RADIUS, catalog="TIC")
        for r in tt:
            if int(r["ID"]) == tic:
                if not np.isfinite(tmag_t): tmag_t = float(r["Tmag"])
                continue
            sep = coord.separation(SkyCoord(float(r["ra"]), float(r["dec"]), unit="deg")).arcsec
            nb.append(dict(tic=int(r["ID"]), sep=sep, tmag=float(r["Tmag"]),
                           dra=(float(r["ra"]) - coord.ra.deg) * np.cos(coord.dec.radian) * 3600,
                           ddec=(float(r["dec"]) - coord.dec.deg) * 3600))
        tic_err = ""
    except Exception as e:
        tic_err = type(e).__name__
    for n in nb:
        n["needed_depth"] = depth * 10 ** (0.4 * (n["tmag"] - tmag_t)) if np.isfinite(tmag_t) else np.nan
    could = [n for n in nb if n["needed_depth"] <= 1]
    kv, kv_err = known_variables(coord, P, tic)
    target_listed = [v for v in kv if v["sep_arcsec"] < 10]
    others = [v for v in kv if v["sep_arcsec"] >= 10]
    eb_words = ("EA", "EB", "EW", "ELL", "EB (")

    def is_eb(v): return any(w in v["type"] for w in eb_words)
    neb_match = [v for v in others if is_eb(v) and period_match(v["period"], P)]
    neb_other = [v for v in others if is_eb(v) and not period_match(v["period"], P)]

    reasons = s.clear_reason
    a_parts, a_ok = [], True
    if "secondary" in reasons:
        surv = np.isfinite(sec_sig) and sec_sig >= 3
        a_parts.append(f"secondary local {sec*1e6:.0f}±{sec_sd*1e6:.0f} ppm ({sec_sig:.1f}σ red-noise) → "
                       + ("survives" if surv else "NOT robust"))
        a_ok &= surv
    if "odd/even" in reasons:
        surv = np.isfinite(oe_sig) and oe_sig >= 3
        a_parts.append(f"odd/even local {oe*1e6:.0f}±{oe_sd*1e6:.0f} ppm ({oe_sig:.1f}σ red-noise) → "
                       + ("survives" if surv else "NOT robust"))
        a_ok &= surv
    other_reason = any(w in reasons for w in ("centroid", "neighbour", "size"))
    if not a_parts:
        a_parts.append("n/a (centroid/neighbour/size flag)")
    if near:
        a_parts.append(f"variability peak {p_peak:.3f} d ≈ {near[0]}×P, ~{amp_peak*1e6:.0f} ppm")
    # A flag survives (a) overall if every variability-sensitive reason survives, or if a
    # variability-insensitive reason (centroid/neighbour/size) is also present.
    a_overall = a_ok or other_reason

    b_parts, b_ok = [], True
    if neb_match:
        b_ok = False
        b_parts.append("known EB at matching period: " + ", ".join(f"{v['name']} ({v['sep_arcsec']:.0f}″, P={v['period']:.4f})" for v in neb_match))
    if neb_other:
        b_parts.append("known EB nearby, other period: " + ", ".join(f"{v['name']} ({v['sep_arcsec']:.0f}″, P={v['period']:.3f})" for v in neb_other))
        if not ("centroid" in reasons or "neighbour" in reasons):
            b_ok = False
    if target_listed:
        b_parts.append("TARGET itself listed: " + ", ".join(f"{v['cat']} {v['type']} P={v['period']:.4f}" for v in target_listed))
    b_parts.append(f"{len(could)} TIC stars within 2′ bright enough to mimic (needed depth ≤ 100%)")
    unclear = (not np.isfinite(sec_sig) and "secondary" in reasons and not other_reason) or \
              (not np.isfinite(oe_sig) and "odd/even" in reasons and not other_reason) or bool(kv_err) or bool(tic_err)
    if unclear:
        cat = "unclear"
    elif a_overall and b_ok:
        cat = "survives both"
    elif not a_overall:
        cat = "variability"
    else:
        cat = "contamination"

    # ---- figure
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    ph = ((t - t0) / P + 0.5) % 1 - 0.5
    n_orb = np.round((t - t0) / P).astype(int)
    w = np.abs(ph * P) < 2 * D
    for par, col, lab in ((0, "tab:blue", "even"), (1, "tab:red", "odd")):
        m = w & (n_orb % 2 == par)
        ax[0, 0].plot(ph[m] * P * 24, f[m], ".", ms=1, alpha=0.3, color=col)
        bins = np.linspace(-2 * D, 2 * D, 50) * 24
        x = ph[m] * P * 24
        bm = [np.mean(f[m][(x >= lo) & (x < hi)]) if ((x >= lo) & (x < hi)).sum() else np.nan for lo, hi in zip(bins[:-1], bins[1:])]
        ax[0, 0].plot(0.5 * (bins[1:] + bins[:-1]), bm, "o-", ms=3, color=col, label=lab)
    ax[0, 0].set_title("Folded transit (odd vs even)"); ax[0, 0].set_xlabel("hours from mid-transit"); ax[0, 0].legend()
    ph2 = ((t - t0) / P) % 1
    m2 = np.abs(ph2 - 0.5) < 0.2
    ax[0, 1].plot(ph2[m2], f[m2], ".", ms=1, alpha=0.2, color="gray")
    bins = np.linspace(0.3, 0.7, 41)
    bm = [np.mean(f[m2][(ph2[m2] >= lo) & (ph2[m2] < hi)]) for lo, hi in zip(bins[:-1], bins[1:])]
    ax[0, 1].plot(0.5 * (bins[1:] + bins[:-1]), bm, "o-", color="red", ms=3)
    ax[0, 1].axvspan(0.5 - 0.5 * D / P, 0.5 + 0.5 * D / P, alpha=0.15)
    ax[0, 1].set_ylim(np.nanmin(bm) - 3 * np.nanstd(bm) - 1e-4, np.nanmax(bm) + 3 * np.nanstd(bm) + 1e-4)
    ax[0, 1].set_title(f"Phase 0.5 window — v2 {h.sec_sigma:.1f}σ; local {sec_sig:.1f}σ (red-noise)")
    phr = ((rt - t0) / P) % 1
    bins = np.linspace(0, 1, 101)
    okr = np.isfinite(rf)
    bm = [np.median(rf[okr][(phr[okr] >= lo) & (phr[okr] < hi)]) for lo, hi in zip(bins[:-1], bins[1:])]
    ax[0, 2].plot(phr[okr], rf[okr], ".", ms=1, alpha=0.1, color="gray")
    ax[0, 2].plot(0.5 * (bins[1:] + bins[:-1]), bm, "-", color="k")
    lo_, hi_ = np.nanpercentile(rf[okr], [1, 99]); ax[0, 2].set_ylim(lo_, hi_)
    ax[0, 2].set_title(f"Un-flattened phase curve: cos2πφ {amp1*1e6:.0f} ppm, cos4πφ {amp2*1e6:.0f} ppm")
    ax[1, 0].plot(1 / freq, pw, lw=0.6)
    for m_, c_ in ((0.5, "g"), (1, "r"), (2, "m")):
        ax[1, 0].axvline(m_ * P, color=c_, ls="--", lw=1, label=f"{m_}×P")
    ax[1, 0].set_xscale("log"); ax[1, 0].set_xlabel("period (d)"); ax[1, 0].legend()
    ax[1, 0].set_title(f"Lomb–Scargle: peak {p_peak:.3f} d")
    ax[1, 1].scatter([0], [0], marker="*", s=200, color="k")
    for n in nb:
        ax[1, 1].scatter(n["dra"], n["ddec"], s=max(5, 120 - 12 * (n["tmag"] - tmag_t)),
                         color="tab:orange" if n["needed_depth"] <= 1 else "lightgray")
    for v in others:
        ax[1, 1].annotate(v["type"][:6], (0, 0), xytext=(5, 5), textcoords="offset points", color="red") if False else None
    ax[1, 1].add_patch(plt.Circle((0, 0), 21, fill=False, ls=":"))
    ax[1, 1].set_xlim(130, -130); ax[1, 1].set_ylim(-130, 130); ax[1, 1].set_aspect("equal")
    ax[1, 1].set_title("TIC stars ≤2′ (orange = bright enough; dotted = 1 px)"); ax[1, 1].set_xlabel("ΔRA″"); ax[1, 1].set_ylabel("ΔDec″")
    ax[1, 2].axis("off")
    txt = [f"TOI {toi:.2f}  TIC {tic}  P = {P:.5f} d  depth ≈ {depth*1e6:.0f} ppm  T14 = {D*24:.1f} h",
           f"v2 flag: {reasons}", f"group: {s.group}", "", "(a) variability / phase curve:"] + \
          ["   " + x for x in a_parts] + ["", "(b) nearby EB:"] + ["   " + x for x in b_parts] + \
          ["", f"ASSESSMENT: {cat.upper()}"] + ([f"catalog errors: {kv_err} {tic_err}"] if (kv_err or tic_err) else [])
    import textwrap
    wrapped = []
    for line in txt:
        wrapped += textwrap.wrap(line, 70, subsequent_indent="      ") or [""]
    ax[1, 2].text(0, 1, "\n".join(wrapped), va="top", family="monospace", fontsize=8)
    fig.suptitle(f"Report card — TOI {toi:.2f} (diagnostic only; frozen v2 results unchanged)")
    fig.tight_layout()
    fig.savefig(f"{outdir}/TOI_{toi:.2f}.png", dpi=110); plt.close(fig)

    return dict(TOI=toi, tic=tic, period=P, depth_ppm=depth * 1e6, t14_h=D * 24, flag=reasons,
                group=s.group, sec_local_ppm=sec * 1e6, sec_local_sigma=sec_sig,
                oe_local_ppm=oe * 1e6, oe_local_sigma=oe_sig, ls_peak_d=p_peak,
                ls_peak_near=",".join(f"{m}xP" for m in near), amp_cos1_ppm=amp1 * 1e6, amp_cos2_ppm=amp2 * 1e6,
                n_tic_2arcmin=len(nb), n_could_mimic=len(could),
                known_neb_match="; ".join(v["name"] for v in neb_match),
                known_eb_other="; ".join(v["name"] for v in neb_other),
                target_listed="; ".join(f"{v['cat']}:{v['type']}" for v in target_listed),
                a_note=" | ".join(a_parts), b_note=" | ".join(b_parts), assessment=cat,
                errors=";".join(kv_err + ([tic_err] if tic_err else [])))


if __name__ == "__main__":
    lst, outdir = sys.argv[1:3]
    os.makedirs(outdir, exist_ok=True)
    hunt, strong = pd.read_csv(HUNT), pd.read_csv(STRONG)
    out_csv = f"{outdir}/report_cards.csv"
    done = set(pd.read_csv(out_csv).TOI.round(2)) if os.path.exists(out_csv) else set()
    for toi in pd.read_csv(lst).TOI:
        if round(toi, 2) in done:
            continue
        rng = np.random.default_rng(RNG_SEED)
        try:
            row = card(toi, hunt, strong, outdir, rng)
        except Exception as e:
            row = dict(TOI=toi, assessment="unclear", errors=f"{type(e).__name__}: {e}"[:200])
        pd.DataFrame([row]).to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
        print(f"TOI {toi:.2f}: {row['assessment']}  {row.get('errors', '')}", flush=True)
