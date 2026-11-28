"""DOWNLOAD ONLY: 2-minute SPOC light curves (no pixel files) for the sibling-search candidate pool.
No analysis, no search; only counts are logged, nothing about individual stars is printed.
Julien's instruction 2026-11-27. Pool definition: docs/sibling_pool_definition.md (saved first).

    python scripts/prefetch_sibling_pool.py define      → data/sibling_pool_2026-11-27.csv
    python scripts/prefetch_sibling_pool.py fetch LOG   (resumable; run under the watchdog)

Pool (written before any download):
  Stars in data/toi_2026-11-26.csv with at least one TOI whose TFOPWG disposition is CP, KP or PC,
  and at least 6 TESS sectors listed in the snapshot's "Sectors" column (pre-filter, any cadence).
  Excluded: every blind_v2 and blind_v3 star (blind samples stay untouched by v3 work).
  Order (Julien 2026-11-27): confirmed-planet hosts (any CP/KP TOI) before PC-only hosts; within
  each, most listed sectors first; then TIC. A partial download thus covers the best targets.
  The exact rule (≥ 6 sectors of 2-min SPOC data) is applied at download time from the MAST
  listing: a star with fewer is recorded as "below 6" and nothing is downloaded for it.
Download: every 2-min SPOC light-curve file of the star (one product per sector), one file at a
  time, 120 s socket time limit per file, 3 attempts per file (waits 30 s, 90 s). Already-cached
  files are skipped by lightkurve. Stops when free disk < 10 GB (checked before each star; was 20 GB, lowered by Julien 2026-11-28).
Records per star go to results/sibling_pool_prefetch_2026-11-27/status.csv (TIC, status, number of
  2-min sectors, files ok / failed); the log has counts only.
"""
import os, sys, time, shutil, socket, warnings
warnings.filterwarnings("ignore")
socket.setdefaulttimeout(120)
import numpy as np
import pandas as pd

SNAP = "data/toi_2026-11-26.csv"
POOL = "data/sibling_pool_2026-11-27.csv"
OUT = "results/sibling_pool_prefetch_2026-11-27"
MIN_FREE_GB, MIN_SECTORS = 10.0, 6


def free_gb():
    return shutil.disk_usage(os.path.expanduser("~")).free / 1e9


def define():
    toi = pd.read_csv(SNAP)
    toi = toi[toi["TFOPWG Disposition"].isin(["CP", "KP", "PC"])].copy()
    toi["n_listed"] = toi["Sectors"].fillna("").astype(str).apply(
        lambda s: len({x.strip() for x in s.split(",") if x.strip()}))
    toi["confirmed"] = toi["TFOPWG Disposition"].isin(["CP", "KP"])
    st = toi.groupby("TIC ID").agg(n_listed=("n_listed", "max"), confirmed=("confirmed", "any")).reset_index()
    st["host"] = np.where(st.confirmed, "CP/KP", "PC")
    excl = set(pd.read_csv("data/blind_v2.csv")["TIC ID"]) | set(pd.read_csv("data/blind_v3.csv")["TIC ID"])
    st = st[(st.n_listed >= MIN_SECTORS) & ~st["TIC ID"].isin(excl)]
    st = st.sort_values(["confirmed", "n_listed", "TIC ID"], ascending=[False, False, True])
    st.to_csv(POOL, index=False)
    print(f"pool: {len(st)} stars ({st.confirmed.sum()} CP/KP hosts, {(~st.confirmed).sum()} PC-only) (listed sectors ≥ {MIN_SECTORS}); median listed sectors "
          f"{st.n_listed.median():.0f}; blind stars excluded")


def fetch(log):
    import lightkurve as lk
    os.makedirs(OUT, exist_ok=True); stat = f"{OUT}/status.csv"
    done = set(pd.read_csv(stat)["TIC ID"]) if os.path.exists(stat) else set()
    pool = [t for t in pd.read_csv(POOL)["TIC ID"] if t not in done]
    L = open(log, "a"); t0 = time.time(); c = {}
    for i, tic in enumerate(pool):
        if free_gb() < MIN_FREE_GB:
            L.write(f"STOPPED: free disk {free_gb():.1f} GB < {MIN_FREE_GB} GB\n"); L.flush(); break
        rec = {"TIC ID": int(tic), "n2min": 0, "ok": 0, "failed": 0}
        try:
            for a in range(3):
                try:
                    sr = lk.search_lightcurve(f"TIC {int(tic)}", mission="TESS", author="SPOC", exptime=120); break
                except Exception:
                    if a == 2: raise
                    time.sleep((30, 90)[a])
            seen, keep = set(), []
            for j, m in enumerate(sr.table["mission"]):
                if str(m) not in seen:
                    seen.add(str(m)); keep.append(j)
            rec["n2min"] = len(keep)
            if len(keep) < MIN_SECTORS:
                rec["status"] = "below6"
            else:
                for j in keep:
                    for a in range(3):
                        try:
                            if sr[j].download() is None:
                                raise ConnectionError("empty")
                            rec["ok"] += 1; break
                        except Exception:
                            if a == 2: rec["failed"] += 1
                            else: time.sleep((30, 90)[a])
                rec["status"] = "done" if rec["failed"] == 0 else "partial"
        except Exception as e:
            rec["status"] = "search_failed"
        pd.DataFrame([rec]).to_csv(stat, mode="a", header=not os.path.exists(stat), index=False)
        c[rec["status"]] = c.get(rec["status"], 0) + 1
        if (i + 1) % 25 == 0 or i + 1 == len(pool):
            L.write(f"{time.strftime('%H:%M')} {i + 1 + len(done)}/{len(pool) + len(done)} {c} | free "
                    f"{free_gb():.1f} GB | {(time.time() - t0) / 60:.0f} min\n"); L.flush()
    else:
        L.write("FINISHED\n")
    L.close()


if __name__ == "__main__":
    define() if sys.argv[1] == "define" else fetch(sys.argv[2])
