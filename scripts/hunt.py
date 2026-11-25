"""v2 hunt: vet every TOI in data/hunt_v2_targets_2026-11-25.csv with the frozen v2 rules.

    python scripts/hunt.py OUT.csv          (resumes: TOIs already in OUT are skipped)

pipeline.py must be identical to tag v2-freeze (checked at start). Saves, for every target, the
fit, every check's status and text, AND every check's underlying numbers (odd/even depths,
secondary-eclipse depth, density ratio, grazing σ, centroid and neighbour-depth details), so v3
can re-vet later without rerunning. The numbers are recomputed with the same frozen functions
(fixed seeds → identical values) and cross-checked against the vetting text. Full nested
details (per-sector centroids, neighbour fit) go to OUT.jsonl, one line per TOI.
"""
import sys, os, json, subprocess, contextlib, io, warnings
sys.path.insert(0, ".")
warnings.filterwarnings("ignore")
import numpy as np
import pandas as pd
import pipeline as pl
from regression import row_for

TARGETS = "data/hunt_v2_targets_2026-11-25.csv"
TARGETS_MD5 = "684506ef9e006f93f19e50b695e3986f"


def _num(x):
    try:
        x = float(x)
        return x if np.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def details(r):
    """Never stops the hunt: a failure here is recorded in details_error."""
    try:
        return _details(r)
    except Exception as e:
        return {"details_error": f"{type(e).__name__}: {e}"[:200]}


def _details(r):
    """Underlying numbers of every check (same frozen functions as vet_star)."""
    lc, popt = r["lc"], pl.get_popt(r)
    with contextlib.redirect_stdout(io.StringIO()):
        d_odd, d_even, oe_sig = pl.odd_even(lc, popt)
        sec_d, sec_e, sec_sig, sec_ratio = pl.secondary_eclipse(lc, popt)
        dens = pl.density_check(r)
    g, g_sig = pl.grazing_check(r)
    v = r["vetting"]
    checks = {n: (s, d) for n, s, d in v["checks"]}
    if oe_sig is not None:                      # cross-check: same σ as in the verdict text
        assert checks["Odd/even depths"][1] == f"{oe_sig:.1f}σ difference", r["star"]
    out = {"oe_depth_odd": _num(d_odd), "oe_depth_even": _num(d_even), "oe_sigma": _num(oe_sig),
           "sec_depth": _num(sec_d), "sec_err": _num(sec_e), "sec_sigma": _num(sec_sig),
           "sec_ratio": _num(sec_ratio), "grazing_g": _num(g), "grazing_sigma": _num(g_sig)}
    if dens is not None:
        out.update(density_cat=_num(dens[0]), density_ratio=_num(dens[1]),
                   density_sigma=_num(dens[2]), density_weak=bool(dens[3]))
    nd = v.get("neighbor_depth")
    if nd:
        out.update(nd_chi2_target=nd["chi2_target"], nd_crit=nd["crit"],
                   nd_n_neighbours=nd["n_neighbours"], nd_n_possible=nd["n_possible"],
                   nd_best=json.dumps(nd["best_neighbour"]))
    c = v.get("centroid")
    out["_json"] = json.dumps({"centroid": None if not c else
                               {k: c[k] for k in ("sigma", "D2", "dof", "dist", "dist_err")}
                               | {"sectors": c["sectors"]},
                               "neighbor_depth": nd, "checks": v["checks"]}, default=float)
    return out


def main(out):
    frozen = subprocess.run(["git", "diff", "--quiet", "v2-freeze", "--", "pipeline.py"])
    assert frozen.returncode == 0, "pipeline.py differs from v2-freeze — refusing to hunt"
    md5 = subprocess.run(["md5", "-q", TARGETS], capture_output=True, text=True).stdout.strip()
    assert md5 == TARGETS_MD5, f"target list changed ({md5})"
    tg = pd.read_csv(TARGETS)
    rows = pd.read_csv(out).to_dict("records") if os.path.exists(out) else []
    done = {r["TOI"] for r in rows}
    jl = open(out.replace(".csv", ".jsonl"), "a")
    for _, t in tg.iterrows():
        if t["TOI"] in done:
            continue
        star = f"TIC {int(t['TIC ID'])}"
        row = {"TOI": t["TOI"], **row_for(star, 2, t["Period (days)"], hint=True, extra=details)}
        js = row.pop("_json", None)
        jl.write(json.dumps({"TOI": t["TOI"], "star": star, "details": js}) + "\n"); jl.flush()
        rows.append(row)
        pd.DataFrame(rows).to_csv(out + ".tmp", index=False)
        os.replace(out + ".tmp", out)
        print(f"TOI {t['TOI']:<9} {star:14s} {row.get('verdict', row.get('error'))}", flush=True)
    print("saved", out)


if __name__ == "__main__":
    main(sys.argv[1])
