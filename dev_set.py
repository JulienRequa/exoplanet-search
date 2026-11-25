"""Development set for vetting rules v2 — DEVELOPMENT ONLY (never reported as an improvement
over v1, never used to tune thresholds by example).

Targets: 13 known planets (pipeline.TARGETS) + the 20 v1 FPs (data/fp_sample_v1.csv)
+ the dev_v2 extension (data/dev_v2.csv: 40 FP + 20 CP/KP, random_state=3), if present.
Every run uses the published / TOI period as period_hint (as the blind test will).

    python dev_set.py run OUT.csv [STAR ...]     run all (or only the named stars)
    python dev_set.py summary OUT.csv            planets kept, FPs caught, which check caught each
"""
import sys
import pandas as pd
import pipeline as pl
from regression import row_for


def targets():
    t = [(s, P, "known planet", None, "") for s, P in pl.TARGETS.items()]
    samp = pd.read_csv("data/fp_sample_v1.csv")
    t += [(f"TIC {int(r['TIC ID'])}", r["Period (days)"], "v1 FP", r["TOI"], "")
          for _, r in samp.iterrows()]
    try:
        dev = pd.read_csv("data/dev_v2.csv")
        t += [(f"TIC {int(r['TIC ID'])}", r["Period (days)"],
               "dev FP" if r["group"] == "FP" else "dev planet", r["TOI"], r["Comments"])
              for _, r in dev.iterrows()]
    except FileNotFoundError:
        pass
    return t


def run(out, only=None):
    rows = []
    for star, P, kind, toi, comment in targets():
        if only and star not in only:
            continue
        row = {"kind": kind, "TOI": toi, **row_for(star, 2, P, hint=True)}
        row["Comments"] = comment
        rows.append(row)
        print(f"{kind:12s} {star:14s} {row.get('verdict', row.get('error'))}", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("saved", out)


def summary(csv):
    d = pd.read_csv(csv, keep_default_na=False, na_values=[""])
    if "error" not in d:
        d["error"] = pd.NA
    chk = [c for c in d.columns if c.startswith("check:")]
    planets = d[d.kind.isin(["known planet", "dev planet"])]
    fps = d[d.kind.isin(["v1 FP", "dev FP"])]

    def split(g):
        err = g["error"].notna()
        nodata = err & g["error"].fillna("").str.contains("No TESS 2-minute")
        mism = ~err & ~g["P match"].isin(["yes", "½×", "2×"])
        return err, nodata, mism, g[~err & ~mism]

    err, nodata, mism, tp = split(planets)
    kept = ~tp["verdict"].str.startswith("LIKELY FALSE")
    print(f"PLANETS ({len(planets)}): no data {nodata.sum()}, other errors {(err & ~nodata).sum()}, "
          f"period mismatch {mism.sum()}, testable {len(tp)} → kept {kept.sum()} ({kept.mean():.0%})")
    print("   verdicts:", tp["verdict"].value_counts().to_dict())
    for _, r in tp[~kept].iterrows():
        print(f"   LOST {r['star']} ({r['kind']}): FAIL in", [c[7:] for c in chk if r[c] == "FAIL"])

    err, nodata, mism, tf = split(fps)
    caught = tf["verdict"].str.startswith("LIKELY FALSE")
    print(f"FPs ({len(fps)}): no data {nodata.sum()}, other errors {(err & ~nodata).sum()}, "
          f"period mismatch {mism.sum()}, testable {len(tf)} → caught {caught.sum()} ({caught.mean():.0%})")
    print("   FAILs per check (testable FPs):", {c[7:]: int((tf[c] == "FAIL").sum()) for c in chk})
    only1 = (tf[chk] == "FAIL").sum(axis=1).eq(1)
    print("   caught by ONLY that check:", {c[7:]: int(((tf[c] == "FAIL") & only1).sum()) for c in chk})
    print("   WARN (any row):", int((d[chk] == "WARN").any(axis=1).sum()),
          "| errors other than no-data:",
          d.loc[d["error"].notna() & ~d["error"].fillna("").str.contains("No TESS"), "star"].tolist())
    return d


if __name__ == "__main__":
    if sys.argv[1] == "run":
        run(sys.argv[2], only=sys.argv[3:] or None)
    else:
        summary(sys.argv[2])
