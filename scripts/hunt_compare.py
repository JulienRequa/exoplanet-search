"""Compare hunt verdicts with the ExoFOP TOI table (status + public TOI comment).

    python scripts/hunt_compare.py HUNT.csv TOI_SNAPSHOT.csv OUTDIR

Writes OUTDIR/strong_fakes.csv, clean_candidates.csv, oddeven_only.csv and prints a summary.
Selection rules for the two short lists (ranking/listing only; vetting thresholds untouched):

Strong fake  = LIKELY FALSE POSITIVE with at least one "clear" FAIL, and an ExoFOP comment that
               raises no false-positive flag. "Clear" = a FAIL that no known planet in the dev
               (124) or blind (99) samples ever got:
                 - secondary eclipse
                 - centroid
                 - neighbour depth
                 - odd/even at >= 6 sigma (known planets lost up to 5.2 sigma, TOI-837 b)
                 - planet size with b below the grazing bound (both size losses had b = 1.2)
Clean        = PLANET CANDIDATE (all checks PASS: no FAIL, WARN, WEAK or N/A), and an ExoFOP
               comment with no false-positive flag.
"""
import re, sys, os
import pandas as pd

FLAG = re.compile(r"v-?\s?shape|\bEB\b|\bNEB\b|\bBEB\b|eclips|binary|secondary|odd|even|"
                  r"centroid|offset|nearby|contamin|blend|\bSB[12]?\b|companion|large ?Rp|"
                  r"too large|synchron|half phase|phase curve|false|\bFP\b|ellipsoidal|"
                  r"off-?target|true source|same pixel|stars in pixel|crowd", re.I)
CHECKS = ["Detection SNR", "Odd/even depths", "Secondary eclipse", "Density", "Planet size",
          "Grazing", "Centroid", "Neighbor depth"]


def clear_reasons(r):
    out = []
    if r["check: Secondary eclipse"] == "FAIL":
        out.append(f"secondary eclipse {r['sec_sigma']:.1f}σ, {100*r['sec_ratio']:.0f}% of depth")
    if r["check: Odd/even depths"] == "FAIL" and r["oe_sigma"] >= 6:
        out.append(f"odd/even {r['oe_sigma']:.1f}σ")
    if r["check: Centroid"] == "FAIL":
        out.append(f"centroid {r['centroid_sigma']:.1f}σ, source {r['centroid_dist']:.2f} px off")
    if r["check: Neighbor depth"] == "FAIL":
        out.append(f"signal fits neighbour ({r['nd_best']})")
    if r["check: Planet size"] == "FAIL" and r["b"] < 1.15:
        out.append(f"size {r['rp_jup']:.1f} R_Jup (b {r['b']:.2f})")
    return "; ".join(out)


def strength(r):
    s = [r["sec_sigma"] if r["check: Secondary eclipse"] == "FAIL" else 0,
         r["oe_sigma"] if r["check: Odd/even depths"] == "FAIL" else 0,
         r["centroid_sigma"] if r["check: Centroid"] == "FAIL" else 0]
    return max([x for x in s if pd.notna(x)] + [0])


if __name__ == "__main__":
    hunt, snap, outdir = sys.argv[1:4]
    os.makedirs(outdir, exist_ok=True)
    h = pd.read_csv(hunt)
    t = pd.read_csv(snap)[["TOI", "TFOPWG Disposition", "Comments"]]
    m = h.merge(t, on="TOI", how="left")
    m["Comments"] = m["Comments"].fillna("")
    m["exofop_flag"] = m["Comments"].str.contains(FLAG)
    tested = m[m["error"].isna()]
    fails = tested[[f"check: {c}" for c in CHECKS]].eq("FAIL")
    nfail = fails.sum(axis=1)

    fp = tested[tested["verdict"] == "LIKELY FALSE POSITIVE"].copy()
    fp["clear_reason"] = fp.apply(clear_reasons, axis=1)
    fp["strength"] = fp.apply(strength, axis=1)
    strong = fp[(fp["clear_reason"] != "") & ~fp["exofop_flag"]].sort_values("strength", ascending=False)

    clean = tested[(tested["verdict"] == "PLANET CANDIDATE") & ~tested["exofop_flag"]]
    clean = clean.sort_values("snr", ascending=False)

    oe_only = tested[(nfail == 1) & fails["check: Odd/even depths"]]

    cols = ["TOI", "star", "period", "rp_jup", "b", "snr", "TFOPWG Disposition", "Comments"]
    strong[["TOI", "star", "clear_reason", "strength"] + cols[2:]].to_csv(f"{outdir}/strong_fakes.csv", index=False)
    clean[cols + ["density_ratio", "oe_sigma"]].to_csv(f"{outdir}/clean_candidates.csv", index=False)
    oe_only[cols + ["oe_sigma"]].to_csv(f"{outdir}/oddeven_only.csv", index=False)
    # Second ranking (Julien 2026-11-26): smallest radius first; flag, don't remove, Rp > 1.8 R_Jup.
    # Targets without a stellar radius (rp_jup missing) go last.
    by_r = clean.sort_values(["rp_jup", "snr"], ascending=[True, False], na_position="last").copy()
    by_r["flag"] = ""
    by_r.loc[by_r["rp_jup"] > 1.8, "flag"] = "large for a planet (possible brown dwarf / low-mass star)"
    by_r[cols + ["density_ratio", "flag"]].to_csv(f"{outdir}/clean_candidates_by_radius.csv", index=False)

    print(f"targets {len(m)}, tested {len(tested)}, no data {m['error'].notna().sum()}")
    print(f"current ExoFOP status of all targets: {m['TFOPWG Disposition'].value_counts().to_dict()}")
    print(f"verdicts: {tested['verdict'].value_counts().to_dict()}")
    print(f"FP verdicts {len(fp)}: with a clear reason {(fp['clear_reason'] != '').sum()}, "
          f"of which ExoFOP comment already flags {((fp['clear_reason'] != '') & fp['exofop_flag']).sum()}; "
          f"strong fakes {len(strong)}")
    print(f"PLANET CANDIDATE (all pass) {(tested['verdict'] == 'PLANET CANDIDATE').sum()}, "
          f"ExoFOP comment flags {((tested['verdict'] == 'PLANET CANDIDATE') & tested['exofop_flag']).sum()}; "
          f"clean {len(clean)}")
    print(f"failed only odd/even: {len(oe_only)}")
