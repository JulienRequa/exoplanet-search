"""Download-only pass over the v2 hunting stars with the current (duplicate-safe) download.
Already-cached files are skipped by lightkurve. Logs counts and the TICs that still fail.

    python scripts/prefetch_hunt.py LOGFILE
"""
import sys, warnings
sys.path.insert(0, ".")
warnings.filterwarnings("ignore")
import pandas as pd
import pipeline as pl

if __name__ == "__main__":
    log = open(sys.argv[1], "a")
    tics = pd.read_csv("data/hunt_v2_targets_2026-11-25.csv")["TIC ID"].drop_duplicates()
    ok, nodata, failed = 0, 0, []
    for i, tic in enumerate(tics):
        try:
            pl.download(f"TIC {int(tic)}", 2); ok += 1
        except ValueError:
            nodata += 1
        except Exception as e:
            failed.append((int(tic), type(e).__name__))
        if (i + 1) % 250 == 0 or i + 1 == len(tics):
            log.write(f"{i+1}/{len(tics)} ok {ok} nodata {nodata} failed {len(failed)}\n"); log.flush()
    log.write(f"FAILED: {failed}\nFINISHED\n"); log.close()
