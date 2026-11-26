"""Follow-up status and originality check for the hunt's strong fakes (2026-11-26).

    python scripts/strong_fakes_check.py CACHE_DIR OUTDIR

CACHE_DIR holds what was downloaded on 2026-11-26 (raw files not committed; sources below):
  exofop_json/<TIC>.json   ExoFOP target page JSON (target.php?id=TIC&json)
  notes_text.json          ExoFOP Observing Notes text per TIC (edit_obsnotes.php?id=TIC)
  tce/*.sel                SPOC DV TCE statistics rows for these TICs, all sectors and multi-sector
                           runs (archive.stsci.edu/missions/tess/catalogs/tce/*_dvr-tcestats.csv)
  raven/*.csv              RAVEN catalog, Lafarga et al. 2026 (MNRAS; arXiv:2603.22597),
                           Zenodo record 19661443
Inputs from the repo: results/hunt_v2_compare_2026-11-26/strong_fakes.csv.

Follow-up (ground) = any TFOP time-series or spectroscopy entry, or an Observing Note.
FOLLOWUP_FINDING: read by hand from the ExoFOP notes and SG1 result files (2026-11-26).
  "resolved" = follow-up already identified a false-positive cause (retired APC/BEB, SB2, EB,
  NEB on a neighbour, or the event ruled out on target); "partial" = points that way but the
  observers call it not conclusive.

SPOC DV flags (standard DV report diagnostics; DV is SPOC's own vetting, not ours):
  secondary: weak-secondary max MES (tce_ws_maxmes) >= 7.1
  centroid:  difference-image offset from the TIC position (tce_ditco_msky) >= 3 sigma
  odd/even:  odd-even depth statistic (tce_bin_oedp_stat) >= 3 sigma
  A TCE matches when its period is within 1% of the TOI period (or ½×, 2×).
  QLP has no public DV reports, so only SPOC TCEs can be checked.
RAVEN: a match (same TIC, period within 1% or ½×/2×) that is in the NSFP >= 0.9 table but not
  in the vetted table counts as "RAVEN did not pass it"; a match in the vetted/validated tables
  means RAVEN judged it planet-like (a disagreement with v2, reported separately).

Groups:
  already resolved        FOLLOWUP_FINDING == resolved
  independently flagged   not resolved, but flagged by SPOC DV, by RAVEN, by partial follow-up
                          evidence, or by a paper cited in the ExoFOP notes
  new                     none of the above
"""
import glob, json, os, re, sys
import numpy as np
import pandas as pd

FOLLOWUP_FINDING = {
    # resolved: follow-up identified the cause
    6847.01: ("resolved", "retired APC/BEB on ExoFOP (chromatic transit)"),
    3384.01: ("resolved", "retired as BEB on ExoFOP (chromatic event, MuSCAT2)"),
    3668.01: ("resolved", "retired APC/BEB on ExoFOP (blended EB, chromatic)"),
    5546.01: ("resolved", "double-lined spectrum (SB2), TRES / ExoFOP-K2 note"),
    6022.02: ("resolved", "EB at 1.93 d per Glusman et al. (cited in Observing Notes)"),
    3518.01: ("resolved", "SG1: probable binary-star eclipse"),
    6379.01: ("resolved", "SG1: NEB detected in TIC 450308270, the same neighbour v2 found"),
    4200.01: ("resolved", "SG1: event ruled out on target; T1/T2 NEB-like or variable"),
    859.01:  ("resolved", "SG1: event ruled out on target twice (ASTEP)"),
    3076.01: ("resolved", "SG1: 17.6 ppt event ruled out on target"),
    2974.01: ("resolved", "SG1: 3.6 ppt event ruled out on target"),
    4998.01: ("resolved", "SG1: ingress ruled out on target"),
    941.01:  ("resolved", "SG1: no transit on target or T2; all neighbours cleared"),
    4728.01: ("resolved", "SG1: V-shaped, chromatic event on target"),
    # partial: points to a false positive, not conclusive
    3192.01: ("partial", "SG1: depth too deep for a planet; NEB-like event on T2; 'if confirmed, retire as EB'"),
    3743.01: ("partial", "SG1: chromatic egress; full transit needed"),
    2057.01: ("partial", "SG1: event likely ruled out on and off target; TIC match to a nearby star"),
    3035.01: ("partial", "SG1: event possibly ruled out on target"),
    5217.01: ("partial", "SG1: event possibly on target or on T2 (ambiguous)"),
}


def ratio_match(p1, p2):
    return any(abs(p1 / (k * p2) - 1) < 0.01 for k in (0.5, 1, 2))


def dv_flags(tce, tic, period):
    rows = tce[(tce.ticid == tic)]
    rows = rows[[ratio_match(p, period) for p in rows.tce_period]] if len(rows) else rows
    if not len(rows):
        return None, ""
    f = set()
    for _, r in rows.iterrows():
        if pd.notna(r.tce_ws_maxmes) and r.tce_ws_maxmes >= 7.1:
            f.add(f"secondary (MES {r.tce_ws_maxmes:.1f})")
        if pd.notna(r.tce_ditco_msky_err) and r.tce_ditco_msky_err > 0 and r.tce_ditco_msky / r.tce_ditco_msky_err >= 3:
            f.add(f"centroid offset ({r.tce_ditco_msky / r.tce_ditco_msky_err:.1f}σ)")
        if pd.notna(r.tce_bin_oedp_stat) and r.tce_bin_oedp_stat >= 3:
            f.add(f"odd/even ({r.tce_bin_oedp_stat:.1f}σ)")
    return len(rows), "; ".join(sorted(f))


if __name__ == "__main__":
    cache, outdir = sys.argv[1:3]
    s = pd.read_csv(f"{outdir}/strong_fakes.csv")
    s["tic"] = s.star.str.split().str[-1].astype(int)
    notes = json.load(open(f"{cache}/notes_text.json"))
    tce = pd.concat([pd.read_csv(f) for f in glob.glob(f"{cache}/tce/*.sel") if os.path.getsize(f)],
                    ignore_index=True)
    tce.columns = [c.strip() for c in tce.columns]
    tce = tce[pd.to_numeric(tce.ticid, errors="coerce").notna()].copy()
    for c in ["ticid", "tce_period", "tce_ws_maxmes", "tce_ditco_msky", "tce_ditco_msky_err", "tce_bin_oedp_stat"]:
        tce[c] = pd.to_numeric(tce[c], errors="coerce")
    nsfp = pd.read_csv(f"{cache}/raven/nsfp09_table.csv")
    vet = pd.read_csv(f"{cache}/raven/vet_table.csv")
    val = pd.read_csv(f"{cache}/raven/val_table.csv")

    rows = []
    for _, r in s.iterrows():
        j = json.load(open(f"{cache}/exofop_json/{r.tic}.json"))
        n_ts, n_spec, n_img = (len(j.get(k) or []) for k in ("time_series", "spectroscopy", "imaging"))
        obs = re.sub(r"\d+\.\d+ PC Active planet candidate \S+ \S+ \S+", "", notes.get(str(r.tic), "")).strip()
        has_fu = (n_ts + n_spec) > 0 or bool(obs)
        fu_kind, fu_text = FOLLOWUP_FINDING.get(round(r.TOI, 2), ("", ""))
        n_tce, dv = dv_flags(tce, r.tic, r.period)
        rv = nsfp[(nsfp.ticid == r.tic).to_numpy() & np.array([ratio_match(p, r.period) for p in nsfp.bls_per], dtype=bool)]
        raven = ""
        if len(rv):
            if r.tic in set(val.ticid):
                raven = "RAVEN validated as planet (disagrees)"
            elif r.tic in set(vet.ticid):
                raven = "RAVEN vetted as planet candidate (disagrees)"
            else:
                raven = "RAVEN did not pass it (NSFP>0.9 but not vetted)"
        raven_flag = raven.startswith("RAVEN did not")
        if fu_kind == "resolved":
            group = "already resolved"
        elif dv or raven_flag or fu_kind == "partial":
            group = "independently flagged"
        else:
            group = "new"
        rows.append(dict(TOI=r.TOI, tic=r.tic, clear_reason=r.clear_reason, strength=r.strength,
                         period=r.period, rp_jup=r.rp_jup, group=group,
                         has_ground_followup=has_fu, n_timeseries=n_ts, n_spectra=n_spec, n_imaging=n_img,
                         observing_note=bool(obs), followup_finding=fu_text,
                         spoc_tce_matches=n_tce or 0, spoc_dv_flags=dv, raven=raven,
                         exofop_comment=r.Comments))
    out = pd.DataFrame(rows)
    out.to_csv(f"{outdir}/strong_fakes_followup.csv", index=False)
    print(out.group.value_counts().to_dict())
    print("ground follow-up:", out.has_ground_followup.sum(), " none:", (~out.has_ground_followup).sum())
    print("SPOC TCE matched:", (out.spoc_tce_matches > 0).sum(), " DV-flagged:", (out.spoc_dv_flags != "").sum())
    print("RAVEN:", out.raven.value_counts().to_dict())
