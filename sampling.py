"""Reproducible TOI samples from a dated ExoFOP snapshot (data/toi_<date>.csv).

    python sampling.py blind    → data/blind_v2.csv  (100 FP + 100 CP/KP, random_state=2)
    python sampling.py dev      → data/dev_v2.csv    (40 FP + 20 CP/KP, random_state=3)

Same cuts as v1: period 0.5–10 d, depth > 2000 ppm. Exclusions are by TIC ID (one star can
host several TOIs), and each pool keeps one TOI per star (lowest TOI number) so no star is
analysed twice. Order of operations matters: blind is drawn first, then dev excludes it.
The blind sample is only counted here, never printed.
"""
import sys
import pandas as pd

SNAPSHOT = "data/toi_2026-11-24.csv"
V1_SAMPLE = "data/fp_sample_v1.csv"
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


if __name__ == "__main__":
    main(sys.argv[1])
