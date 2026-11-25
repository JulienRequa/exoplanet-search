"""Pre-cache light curves (DOWNLOAD ONLY) for blind_v2 and the v2 hunting targets.
No pipeline runs and no per-star output: only counts are logged, so nothing about the blind
sample is revealed. First 2 sectors per star (what the pipeline uses), with download retries.
Stops if free disk space drops below MIN_FREE_GB. Restartable (lightkurve skips cached files).

    python scripts/precache.py LOGFILE
"""
import sys, os, shutil, time, warnings
sys.path.insert(0, ".")
warnings.filterwarnings("ignore")
import pandas as pd
import pipeline as pl

MIN_FREE_GB = 8.0

def targets():
    toi = pd.read_csv("data/toi_2026-11-24.csv")
    pc = toi[(toi["TFOPWG Disposition"] == "PC") & toi["Period (days)"].between(0.5, 10)
             & (toi["Depth (ppm)"] > 2000)]
    blind = list(pd.read_csv("data/blind_v2.csv")["TIC ID"])
    hunt = [t for t in pc.drop_duplicates("TIC ID")["TIC ID"] if t not in set(blind)]
    return [("blind", t) for t in blind] + [("hunt", t) for t in hunt]

def free_gb():
    return shutil.disk_usage(os.path.expanduser("~")).free / 1e9

if __name__ == "__main__":
    log = open(sys.argv[1], "a")
    tg = targets()
    n = {"blind": [0, 0, 0], "hunt": [0, 0, 0]}          # downloaded, no 2-min data, errors
    t0 = time.time()
    for i, (grp, tic) in enumerate(tg):
        if i % 25 == 0 and free_gb() < MIN_FREE_GB:
            log.write(f"STOPPED at {i}/{len(tg)}: free disk {free_gb():.1f} GB < {MIN_FREE_GB} GB\n"); break
        try:
            pl.download(f"TIC {int(tic)}", 2)
            n[grp][0] += 1
        except ValueError:
            n[grp][1] += 1
        except Exception:
            n[grp][2] += 1
        if (i + 1) % 100 == 0 or i + 1 == len(tg):
            log.write(f"{time.strftime('%H:%M')} {i+1}/{len(tg)} | blind done/nodata/err {n['blind']} | "
                      f"hunt {n['hunt']} | free {free_gb():.1f} GB | {(time.time()-t0)/60:.0f} min\n"); log.flush()
    log.write("FINISHED\n"); log.close()
