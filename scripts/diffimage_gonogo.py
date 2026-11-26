"""v3 go/no-go test for the difference-image method: positional accuracy on real data.
DEV DATA ONLY (dev sample 2 = data/dev_v3.csv); blind_v2 is not touched. Diagnostic; v2 unchanged.

    python scripts/diffimage_gonogo.py OUTDIR

Samples (from results/devfull200_multistart_2026-11-25.csv, rows where v2 found the catalog
period):
  - planets: dev3 confirmed/known planets → the dimming IS on the target.
  - NEB controls: dev3 TFOP false positives whose ExoFOP comment says NEB / off-target / centroid
    offset → the dimming is (mostly) NOT on the target.
Method: exactly scripts/diffimage_check.py (same functions), first SPOC 2-min sector only
(to limit downloads). Sectors with difference-image S/N < 5 count as "no detection".

PRE-REGISTERED CRITERIA (written 2026-11-26, before running):
  Floor  = 95th percentile of |source − target| (px) among detected planets.
  Accuracy = fraction of detected planets within 0.5 px of the target.
  Error calibration = fraction of detected planets with offset/σ > 3 (≈ 0.3% if errors are right).
  GO (an "on target" result counts as evidence that the dimming is within the floor of the
     target) if Accuracy ≥ 90% AND Error calibration ≤ 10%.
  PARTIAL GO if Accuracy ≥ 90% but Error calibration > 10%: use the empirical floor (offset
     in px) instead of σ.
  NO-GO otherwise: "on target" results are not used as evidence.
  Sensitivity (reported, not part of GO) = fraction of detected NEB controls classified "off
  target" by offset > floor.
Outputs: OUTDIR/gonogo.csv (per star), OUTDIR/GONOGO.md (numbers + verdict).
"""
import os, sys, socket, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "."); sys.path.insert(0, "scripts")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.coordinates import SkyCoord
from diffimage_check import download_tpfs, sector_result
from report_cards import t14

DEV = "results/devfull200_multistart_2026-11-25.csv"


def samples():
    d = pd.read_csv(DEV)
    d = d[(d["P match"] == "yes") & d.t0.notna()]
    meta = pd.read_csv("data/dev_v3.csv")
    neb = meta[(meta.group == "FP") & meta.Comments.fillna("").str.contains(
        r"NEB|off-target|off target|offset", case=False)].TOI
    pl_ = d[d.kind == "dev3 planet"].assign(sample="planet")
    nb = d[(d.kind == "dev3 FP") & d.TOI.isin(neb)].assign(sample="NEB control")
    return pd.concat([pl_, nb])


if __name__ == "__main__":
    out = sys.argv[1]; os.makedirs(out, exist_ok=True)
    csv = f"{out}/gonogo.csv"
    done = set(pd.read_csv(csv).TOI.round(2)) if os.path.exists(csv) else set()
    from astroquery.mast import Catalogs
    for _, h in samples().iterrows():
        if round(h.TOI, 2) in done:
            continue
        tic = int(str(h.star).split()[-1])
        D = t14(h.period, h["Rp/Rs"], h["a/Rs"], h.b, h.inc)
        rec = dict(TOI=h.TOI, tic=tic, sample=h["sample"])
        try:
            tpf = download_tpfs(tic, n=1)[0]
            hd = tpf.hdu[0].header
            c0 = SkyCoord(hd["RA_OBJ"], hd["DEC_OBJ"], unit="deg")
            tt = Catalogs.query_region(c0, radius=60 * u.arcsec, catalog="TIC")
            nb = [dict(tic=int(r["ID"]), ra=float(r["ra"]), dec=float(r["dec"]), tmag=float(r["Tmag"])) for r in tt]
            r = sector_result(tpf, tic, h.t0, h.period, D, nb)
            if r is None:
                rec.update(result="too few transits")
            else:
                rec.update(sector=r["sector"], n_transits=r["n_transits"], diff_snr=r["diff_snr"],
                           off_px=r["off_px"], off_err=r["off_err"], off_sigma=r["off_sigma"],
                           result=r["result"], nearest=r["nearest"])
        except Exception as e:
            rec.update(result="error", note=f"{type(e).__name__}: {e}"[:150])
        pd.DataFrame([rec]).to_csv(csv, mode="a", header=not os.path.exists(csv), index=False)
        print(f"{h['sample']:12s} TOI {h.TOI:.2f}: {rec.get('result')} {rec.get('off_px', np.nan):.2f} px", flush=True)

    g = pd.read_csv(csv)
    det = g[g.result.isin(["on target", "off target", "inconclusive"])]
    p = det[det["sample"] == "planet"]
    floor = np.nanpercentile(p.off_px, 95)
    acc = (p.off_px <= 0.5).mean()
    cal = (p.off_sigma > 3).mean()
    verdict = ("GO" if acc >= 0.9 and cal <= 0.1 else "PARTIAL GO (use empirical floor)" if acc >= 0.9 else "NO-GO")
    n = det[det["sample"] == "NEB control"]
    sens = (n.off_px > floor).mean() if len(n) else np.nan
    lines = [f"# Difference-image go/no-go (dev data only) — {pd.Timestamp.now():%Y-%m-%d}", "",
             f"- Planets: {len(g[g['sample']=='planet'])} tried, {len(p)} detected (diff S/N ≥ 5).",
             f"- Accuracy (within 0.5 px): {acc:.1%}  (criterion ≥ 90%)",
             f"- Error calibration (offset > 3σ): {cal:.1%}  (criterion ≤ 10%)",
             f"- Floor (95th percentile offset): {floor:.2f} px = {floor*21:.1f}″",
             f"- Median planet offset: {np.nanmedian(p.off_px):.2f} px",
             f"- NEB controls: {len(g[g['sample']=='NEB control'])} tried, {len(n)} detected; "
             f"offset > floor: {sens:.1%}",
             "", f"**Verdict: {verdict}**"]
    open(f"{out}/GONOGO.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
