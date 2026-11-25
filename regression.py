"""Regression run: the 13 known planets (2 sectors) + WASP-18 (4 sectors).

    python regression.py run  results/known_<label>_<date>.csv
    python regression.py compare OLD.csv NEW.csv

`run` saves every fitted number and vetting result so later code changes can be checked.
`compare` lists every value that changed, in units of its error bar.
"""
import sys, io, contextlib, time
import numpy as np
import pandas as pd
import pipeline as pl

RUNS = [(star, 2, P) for star, P in pl.TARGETS.items()] + [("WASP-18", 4, 0.9415)]


def row_for(star, n_sectors, P_pub, hint=False):
    row = {"star": star, "sectors_requested": n_sectors}
    t = time.time()
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            r = pl.analyze_star(star, n_sectors=n_sectors, verbose=False, vet=True,
                                period_hint=P_pub if hint else None)
    except Exception as e:
        row["error"] = f"{type(e).__name__}: {e}"[:200]
        return row
    row["n_sectors"] = r["n_sectors"]
    for name in pl.PARAM_NAMES:
        row[name], row[name + "_err"] = r[name]
    row["inc"], row["inc_err"] = r["inc"]
    row["density"], row["density_err"] = r["density"]
    row["flags"] = ", ".join(r["flags"])
    row["b_upper"] = r.get("b_upper", "")
    row["inc_lower"] = r.get("inc_lower", "")
    row["P match"] = pl.period_match(r["period"][0], P_pub)
    for k in ("star_source", "period_choice", "ld_source", "ephem_notes"):          # added by later versions
        if k in r:
            row[k] = r[k]
    v = r["vetting"]
    row["snr"], row["rp_jup"] = v["snr"], v["rp_jup"]
    for name, status, detail in v["checks"]:
        row["check: " + name] = status
        row["detail: " + name] = detail
    row["verdict"], row["n_pass"], row["n_run"] = v["verdict"], v["n_pass"], v["n_run"]
    c = v.get("centroid")
    if c:
        row["centroid_sigma"], row["centroid_dist"], row["centroid_dist_err"] = c["sigma"], c["dist"], c["dist_err"]
    for k in ("CROWDSAP", "TESSMAG"):
        row[k] = r["raw_lcs"][0].meta.get(k) if "raw_lcs" in r else None
    row["time_s"] = round(time.time() - t)
    return row


def run(out):
    rows = []
    for star, n, P in RUNS:
        rows.append(row_for(star, n, P))
        print(f"{star} ({n} sectors): {rows[-1].get('verdict', rows[-1].get('error'))}", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("saved", out)


def compare(old_csv, new_csv):
    a, b = pd.read_csv(old_csv), pd.read_csv(new_csv)
    key = ["star", "sectors_requested"]
    m = a.merge(b, on=key, suffixes=("_old", "_new"))
    nums = pl.PARAM_NAMES + ["inc", "density"]
    worst = 0
    for _, r in m.iterrows():
        notes = []
        for p in nums:
            d = r[p + "_new"] - r[p + "_old"]
            err = r[p + "_err_old"]
            flagged = p in str(r["flags_old"]).split(", ") or (p == "inc" and "b" in str(r["flags_old"]))
            if d != 0 and (flagged or not err > 0):  # at a bound: its ± is meaningless
                if abs(d) > 1e-4:
                    notes.append(f"{p} {r[p+'_old']:.5g}→{r[p+'_new']:.5g} (at a bound, ± unreliable)")
                continue
            if d != 0:
                s = abs(d) / err
                worst = max(worst, s)
                if s > 0.1:
                    notes.append(f"{p} {r[p+'_old']:.5g}→{r[p+'_new']:.5g} ({s:.2f}σ)")
        for c in [c for c in a.columns if c.startswith("check:") or c in ("verdict", "flags", "P match")]:
            o, n = r.get(c + "_old"), r.get(c + "_new")
            if not (pd.isna(o) and pd.isna(n)) and o != n:
                notes.append(f"{c}: {o} → {n}")
        print(f"{r['star']:10s} ({r['sectors_requested']}): " + ("; ".join(notes) if notes else "unchanged (all < 0.1σ)"))
    print(f"\nLargest parameter change: {worst:.3f}σ")


if __name__ == "__main__":
    {"run": lambda: run(sys.argv[2]), "compare": lambda: compare(sys.argv[2], sys.argv[3])}[sys.argv[1]]()
