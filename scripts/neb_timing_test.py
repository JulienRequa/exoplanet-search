"""Nearby-binary timing test for report-card "contamination" cases (Julien, 2026-11-26).

    python scripts/neb_timing_test.py CARDS_DIR

Question: could a catalogued nearby eclipsing binary (EB) whose period is NOT P, 2P or P/2 put
its eclipses at our observed secondary-eclipse times (phase 0.5) in BOTH sectors?

For every TOI whose card lists "known EB nearby, other period":
 - Our secondary windows: t0 + (m + ½)P with ≥ 3 data points within ±0.4·T14, per sector.
 - EB ephemeris from the catalogue:
     Gaia DR3 (I/358/veb): P = 1/Freq, reference BJD = TimeRef + 2455197.5, primary and
       secondary eclipse phases PhaseGauss1/2, half-width 2·sigPhaseGauss·P (≈95% of the dip).
     VSX (B/vsx): Period, Epoch (HJD ≈ BJD for this purpose); primary only plus a secondary at
       phase 0.5 (assumed); half-width 0.08·P (assumed; wide, i.e. conservative).
   Timing error at TESS epochs = cycles × σ_P (+ epoch error); σ_P from e_Freq (Gaia) or half
   the last quoted digit (VSX).
 - A window is "hit" if any EB eclipse centre lies within 0.4·T14 + half-width + 2σ_t.
 - "with catalogue epoch": fraction of windows hit per sector (only if 2σ_t < P_EB/4).
 - "worst case": the largest fraction over all possible EB epochs (200 shifts), per sector —
   valid even when the epoch is unknown.
Decision (written before running): the EB CAN explain the secondary if the hit fraction is
≥ 50% in every sector (catalogue-epoch fraction if available, else worst case). Otherwise the
star moves to "survives" for (b). An EB with no catalogued period cannot be tested → stays.
Writes CARDS_DIR/neb_timing.csv.
"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import numpy as np
import pandas as pd
import pipeline as pl

BTJD0, GAIA0 = 2457000.0, 2455197.5


def eb_ephemeris(name):
    from astroquery.vizier import Vizier
    v = Vizier(columns=["**"], row_limit=5, vizier_server="vizier.cfa.harvard.edu")
    if name.startswith("Gaia"):
        src = name.split()[-1]
        r = v.query_constraints(catalog="I/358/veb", Source=src)[0][0]
        P = 1 / float(r["Freq"]); sP = float(r["e_Freq"]) * P ** 2
        tref = float(r["TimeRef"]) + GAIA0 - BTJD0
        ecl = []
        for k in ("1", "2"):
            ph, sg = r[f"PhaseGauss{k}"], r[f"sigPhaseGauss{k}"]
            if ph is not None and np.isfinite(float(ph)):
                ecl.append((float(ph), 2 * float(sg) * P))
        if not ecl:
            ph = float(r["PhaseE1"]); ecl = [(ph, 0.08 * P), (ph + 0.5, 0.08 * P)]
        return dict(P=P, sP=sP, tref=tref, s_ep=0.0, ecl=ecl, src="Gaia DR3 veb")
    r = v.query_constraints(catalog="B/vsx/vsx", Name=name)[0][0]
    if not r["Period"] or not np.isfinite(float(r["Period"])):
        return None
    P = float(r["Period"]); digits = len(str(r["Period"]).split(".")[-1])
    sP = 0.5 * 10 ** (-digits)
    ep = float(r["Epoch"]) - BTJD0 if r["Epoch"] else np.nan
    return dict(P=P, sP=sP, tref=ep, s_ep=0.01 * P, ecl=[(0.0, 0.08 * P), (0.5, 0.08 * P)], src="VSX")


def windows(lcs, t0, P, half):
    out = []
    for lc in lcs:
        t = np.asarray(lc.remove_nans().time.value, float)
        m = np.arange(np.floor((t.min() - t0) / P) - 1, np.ceil((t.max() - t0) / P) + 1)
        c = t0 + (m + 0.5) * P
        c = [x for x in c if (np.abs(t - x) < half).sum() >= 3]
        out.append((int(lc.meta.get("SECTOR") or 0), np.array(c)))
    return out


def hit_fraction(centres, eb, tol_extra, shift=0.0):
    hits = 0
    for tc in centres:
        ok = False
        for ph, hw in eb["ecl"]:
            n = np.round((tc - eb["tref"] - shift - ph * eb["P"]) / eb["P"])
            te = eb["tref"] + shift + (n + ph) * eb["P"]
            if abs(te - tc) < hw + tol_extra:
                ok = True
        hits += ok
    return hits / len(centres) if len(centres) else np.nan


if __name__ == "__main__":
    d = sys.argv[1]
    cards = pd.read_csv(f"{d}/report_cards.csv")
    hunt = pd.read_csv("results/hunt_v2_2026-11-25.csv")
    rows = []
    for _, c in cards[cards.known_eb_other.notna()].iterrows():
        h = hunt[hunt.TOI == c.TOI].iloc[0]
        half = 0.4 * c.t14_h / 24
        lcs = pl.download(f"TIC {int(c.tic)}", 2)
        wins = windows(lcs, h.t0, h.period, half)
        for name in [x.strip() for x in c.known_eb_other.split(";")]:
            try:
                eb = eb_ephemeris(name)
            except Exception as e:
                rows.append(dict(TOI=c.TOI, eb=name, result="not testable", note=f"{type(e).__name__}")); continue
            if eb is None:
                rows.append(dict(TOI=c.TOI, eb=name, result="not testable", note="no catalogued period")); continue
            per_sector, can_all = [], True
            for sec, cs in wins:
                cyc = np.abs(np.mean(cs) - eb["tref"]) / eb["P"] if len(cs) else 0
                st = cyc * eb["sP"] + eb["s_ep"]
                use_epoch = np.isfinite(eb["tref"]) and 2 * st < eb["P"] / 4
                f_ep = hit_fraction(cs, eb, half + 2 * st) if use_epoch else np.nan
                f_wc = max(hit_fraction(cs, eb, half, s) for s in np.linspace(0, eb["P"], 200, endpoint=False))
                f = f_ep if use_epoch else f_wc
                can_all &= bool(np.isfinite(f) and f >= 0.5)
                per_sector.append(f"s{sec}: n={len(cs)}, epoch {f_ep:.2f}, worst {f_wc:.2f}, σt={st*24:.1f} h")
            rows.append(dict(TOI=c.TOI, eb=name, eb_source=eb["src"], eb_period=eb["P"], target_period=h.period,
                             result="can explain" if can_all else "cannot explain", note=" | ".join(per_sector)))
            print(c.TOI, name, rows[-1]["result"], rows[-1]["note"], flush=True)
    pd.DataFrame(rows).to_csv(f"{d}/neb_timing.csv", index=False)
