"""Development set for vetting rules v2: 13 known planets + the 20 v1 false positives.
All runs use the published / TOI period as period_hint (as the v2 blind test will).
These results are for DEVELOPMENT ONLY: v1 is never re-scored as an improvement.

    python dev_set.py results/devset_<label>_<date>.csv
"""
import sys
import pandas as pd
import pipeline as pl
from regression import row_for

def main(out):
    v1 = pd.read_csv("results/fp_test_v1.csv")
    targets = [(s, P, "known planet") for s, P in pl.TARGETS.items()]
    targets += [(r["star"], r["published P"], "v1 FP") for _, r in v1.iterrows()]
    rows = []
    for star, P, kind in targets:
        row = {"kind": kind, **row_for(star, 2, P, hint=True)}
        if kind == "v1 FP":
            row["TOI"] = v1.set_index("star").loc[star, "TOI"]
        rows.append(row)
        print(f"{kind:12s} {star:14s} {row.get('verdict', row.get('error'))}", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("saved", out)

if __name__ == "__main__":
    main(sys.argv[1])
