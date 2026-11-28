"""Reproducible TOI samples from a dated ExoFOP snapshot (data/toi_<date>.csv).

    python sampling.py blind    → data/blind_v2.csv  (100 FP + 100 CP/KP, random_state=2)
    python sampling.py dev      → data/dev_v2.csv    (40 FP + 20 CP/KP, random_state=3)
    python sampling.py blind3   → data/blind_v3.csv  (100 FP + 100 CP/KP, random_state=5, snapshot
                                  2026-11-26; excludes every star used so far; never printed)
    python sampling.py hunt     → data/hunt_v2_targets_2026-11-25.csv (all open PC TOIs in the
                                  v2 tested range; every TOI kept, even several per star)

Same cuts as v1: period 0.5–10 d, depth > 2000 ppm. Exclusions are by TIC ID (one star can
host several TOIs), and each pool keeps one TOI per star (lowest TOI number) so no star is
analysed twice. Order of operations matters: blind is drawn first, then dev excludes it.
The blind sample is only counted here, never printed.
"""
import sys
import pandas as pd

SNAPSHOT = "data/toi_2026-11-24.csv"
V1_SAMPLE = "data/fp_sample_v1.csv"
HUNT_LIST = "data/hunt_v2_targets_2026-11-25.csv"
# Hosts of the 13 known-planet regression targets (pipeline.TARGETS)
KNOWN_TICS = {"WASP-18": 100100827, "WASP-121": 22529346, "WASP-19": 35516889,
              "WASP-43": 36734222, "HD 189733": 256364928, "HD 209458": 420814525,
              "KELT-9": 16740101, "WASP-33": 129979528, "HAT-P-7": 424865156,
              "WASP-62": 149603524, "WASP-100": 38846515, "WASP-126": 25155310,
              "pi Men": 261136679}   # MAST name → nearest TIC (all 0.0″), 2026-11-24


def pool(toi, dispositions, exclude_tics):
    p = toi[toi["TFOPWG Disposition"].isin(dispositions)
            & toi["Period (days)"].between(0.5, 10)
            & (toi["Depth (ppm)"] > 2000)
            & ~toi["TIC ID"].isin(exclude_tics)]
    return p.sort_values("TOI").drop_duplicates("TIC ID", keep="first")


def draw(toi, n_fp, n_pl, seed, exclude_tics):
    fp = pool(toi, ["FP"], exclude_tics).sample(n_fp, random_state=seed)
    pl = pool(toi, ["CP", "KP"], exclude_tics).sample(n_pl, random_state=seed)
    out = pd.concat([fp.assign(group="FP"), pl.assign(group="CP/KP")])
    return out[["group", "TOI", "TIC ID", "TFOPWG Disposition", "Period (days)",
                "Depth (ppm)", "TESS Mag", "Comments"]]


def main(which):
    toi = pd.read_csv(SNAPSHOT)
    exclude = set(pd.read_csv(V1_SAMPLE)["TIC ID"]) | set(KNOWN_TICS.values())
    if which == "blind":
        s = draw(toi, 100, 100, 2, exclude)
        s.to_csv("data/blind_v2.csv", index=False)
        print("blind_v2:", s["group"].value_counts().to_dict(), "| unique TICs:", s["TIC ID"].nunique())
    elif which == "dev":
        exclude |= set(pd.read_csv("data/blind_v2.csv")["TIC ID"])
        s = draw(toi, 40, 20, 3, exclude)
        s.to_csv("data/dev_v2.csv", index=False)
        print("dev_v2:", s["group"].value_counts().to_dict(), "| unique TICs:", s["TIC ID"].nunique())
    elif which == "dev3":           # "dev sample 2" (2026-11-24): finding problems, not reporting
        exclude |= set(pd.read_csv("data/blind_v2.csv")["TIC ID"])
        exclude |= set(pd.read_csv("data/dev_v2.csv")["TIC ID"])
        s = draw(toi, 100, 100, 4, exclude)
        s.to_csv("data/dev_v3.csv", index=False)
        print("dev_v3:", s["group"].value_counts().to_dict(), "| unique TICs:", s["TIC ID"].nunique())
    elif which == "blind3":         # v3 blind sample (Julien 2026-11-27), locked before further v3 work
        blind3()
        return
    elif which == "hunt":           # v2 hunting list (2026-11-25): open candidates, no sampling
        for f in ("data/blind_v2.csv", "data/dev_v2.csv", "data/dev_v3.csv"):
            exclude |= set(pd.read_csv(f)["TIC ID"])
        h = toi[(toi["TFOPWG Disposition"] == "PC")
                & toi["Period (days)"].between(0.5, 10)
                & (toi["Depth (ppm)"] > 2000)
                & ~toi["TIC ID"].isin(exclude)].sort_values("TOI")
        h = h[["TOI", "TIC ID", "TFOPWG Disposition", "Period (days)", "Depth (ppm)",
               "TESS Mag", "Comments"]]
        h.to_csv(HUNT_LIST, index=False)
        print("hunt_v2:", len(h), "TOIs on", h["TIC ID"].nunique(), "stars")
        return
    # how many matching FPs/planets are still unused by any sample
    used = exclude | (set(pd.read_csv("data/dev_v3.csv")["TIC ID"]) if which == "dev3" else set())
    print("unused matching stars — FP:", len(pool(toi, ["FP"], used)),
          "| CP/KP:", len(pool(toi, ["CP", "KP"], used)))


BLIND3_SNAPSHOT = "data/toi_2026-11-26.csv"   # latest dated snapshot
BLIND3_SEED = 5                               # new: 2 (blind_v2), 3 (dev_v2), 4 (dev_v3) used before


def used_tics():
    """Every star used so far, by TIC: v1, the 13 known hosts, blind_v2, dev_v2, dev_v3, the whole
    v2 hunt list (all 2,896 TOIs were run; includes every hunt flag reviewed in detail), and any star
    named in any results/ or private/ table (TIC columns, "TIC n" star names, or TOI numbers mapped to
    TIC through both snapshots). Deliberately broader than needed."""
    import glob, re
    toi = pd.concat([pd.read_csv("data/toi_2026-11-24.csv"), pd.read_csv(BLIND3_SNAPSHOT)])
    toi2tic = dict(zip(toi["TOI"].round(2), toi["TIC ID"]))
    used = set(KNOWN_TICS.values())
    files = [f for f in glob.glob("data/*.csv") if "toi_" not in f and "claret" not in f]
    files += glob.glob("results/**/*.csv", recursive=True) + glob.glob("private/**/*.csv", recursive=True)
    for f in files:
        try:
            d = pd.read_csv(f, low_memory=False)
        except Exception:
            continue
        for c in d.columns:
            cl = c.strip().lower()
            if cl in ("tic id", "tic"):
                used |= set(pd.to_numeric(d[c], errors="coerce").dropna().astype(int))
            elif cl == "star":
                m = d[c].astype(str).str.extract(r"TIC\s*(\d+)")[0].dropna()
                used |= set(m.astype(int))
            elif cl == "toi":
                t = pd.to_numeric(d[c], errors="coerce").dropna().round(2)
                used |= {int(toi2tic[x]) for x in t if x in toi2tic}
    return used


def blind3():
    import hashlib, os
    out = "data/blind_v3.csv"
    assert not os.path.exists(out), "blind_v3 already locked; never redraw"
    toi = pd.read_csv(BLIND3_SNAPSHOT)
    used = used_tics()
    s = draw(toi, 100, 100, BLIND3_SEED, used)
    assert not (set(s["TIC ID"]) & used) and s["TIC ID"].nunique() == 200
    s.to_csv(out, index=False)
    raw = open(out, "rb").read()
    fp = (f"file: {out}\nsnapshot: {BLIND3_SNAPSHOT}\nseed: {BLIND3_SEED}\nrows: {len(s)} "
          f"({(s.group == 'FP').sum()} FP, {(s.group == 'CP/KP').sum()} CP/KP)\n"
          f"md5: {hashlib.md5(raw).hexdigest()}\nsha256: {hashlib.sha256(raw).hexdigest()}\n"
          f"excluded stars (used so far): {len(used)}\n")
    open("data/blind_v3.fingerprint.txt", "w").write(fp)
    print(fp, end="")   # counts and hashes only; the stars are never printed
    left = used | set(s["TIC ID"])
    print("unused matching stars after blind_v3 — FP:", len(pool(toi, ["FP"], left)),
          "| CP/KP:", len(pool(toi, ["CP", "KP"], left)))


if __name__ == "__main__":
    main(sys.argv[1])
