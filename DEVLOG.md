# Development log

## 2026-11-24 — False-positive test v1 (official blind result, rules frozen 2026-11-24)

Sample: 20 TFOPWG "FP" TOIs (period 0.5–10 d, depth > 2000 ppm, `random_state=1`).

- 3 had no 2-min data; 2 crashed; 3 search found a different signal (period mismatch)
- 12 testable → 5 caught as LIKELY FALSE POSITIVE (≈42%), 7 passed as planet candidates

**Do NOT re-score this sample with new code and report it as an improvement.**
Improvements are evaluated only on a fresh sample.

### Bugs found (fixes are bug fixes, not threshold changes)

1. Crash "Initial guess is outside of provided bounds" in `fit_transit`: for very deep
   eclipses, the guess Rp/Rs = sqrt(depth) exceeds the upper bound 0.30. Fix: raise the
   Rp/Rs upper bound to 1.0 (eclipsing binaries can be star-sized) and clip every initial
   guess to lie inside its bounds. One TIC also hit exactly Rp/Rs = 0.3000, the cap.
2. Missing stellar radius: some targets have no RADIUS/LOGG in the light-curve header, so
   the size and density checks become N/A. Fix: fall back to the TIC catalog (radius,
   log g, Teff) only when header values are missing; record which source was used.
3. Period hint: optional `period_hint` for `analyze_star` (and batch targets). When given,
   restrict BLS to ±2% of the hint (also try ½× and 2×, report which matched), so we vet
   the TOI's actual signal instead of a different one.

### Next

- Proposed (not yet approved): grazing / V-shape check (b + Rp/Rs > 1).
- Rerun 13 known-planet regression (n_sectors=2): expect 13/13, same verdicts.
- Freeze vetting rules v2 (VETTING_RULES.md), snapshot ExoFOP TOI table, then blind test
  on a fresh sample (`random_state=2`, excluding v1): 100 FP + 100 CP/KP TOIs.

## 2026-11-24 — Review of v1 misses (from TOI Comments)

Of the 7 FPs that passed as planet candidates in v1:
- 4 nearby-star contamination (centroid offset / signal on a neighbor)
- 2 V-shaped
- 1 needed a stellar radius and more sectors
- 1 unknown

### Decisions
- Bug 1 accepted, including its side effect: the "at a boundary" margin for Rp/Rs grows
  from 0.0029 to 0.0099 (1% of the new 0.01–1.0 range). `profile_b`'s Rp/Rs cap is also
  raised to 1.0 to match `fit_transit`.
- Bug 2: NaN header values count as missing (a NaN radius used to make the size check
  FAIL and a NaN log g made the density check PASS).
- Bug 3: with a period hint, BLS is run within ±2% of 1×, ½× and 2× the hint; the one with
  the highest BLS power wins and is recorded. Candidates longer than half the span are skipped.
- Grazing check APPROVED as a WARNING (never a FAIL): grazing if g = b + Rp/Rs − 1 > 0 by
  more than 3σ, σ_g from the fit covariance (including the b–Rp/Rs covariance term).
  T14/T23 shape comparison skipped for now.
- Verdict becomes three levels: PLANET CANDIDATE (all pass) / CANDIDATE WITH WARNINGS
  (no FAIL, ≥1 warning) / LIKELY FALSE POSITIVE (any FAIL).
- Edge case: white-dwarf hosts can give Rp/Rs >> 0.3; the size check must report this sensibly.
- Centroid check (target pixel files, in- vs out-of-transit centroid shift): design to be
  shown before building.
- Order: setup → baseline → bug fixes (one commit each) → regression → grazing + 3-level
  verdict → regression → centroid design (approval) → freeze v2 → blind test
  (100 FP + 100 CP/KP, random_state=2, excluding the 20 v1 TICs).

### v1 sample (exclude from all future samples)
TIC (TOI): 29336129 (2933.01), 269450900 (1127.01), 263203229 (6741.01), 55390413 (6994.01),
290348382 (1113.01), 28873762 (3684.01), 68322189 (5016.01), 66356824 (315.01),
308396509 (2625.01), 144043410 (2554.01), 686072378 (2438.01), 146846569 (734.01),
441738827 (2084.02), 231076038 (2618.01), 285233023 (6298.01), 294315741 (4985.01),
146261607 (591.01), 140706664 (1020.01), 306547699 (3216.01), 279769094 (4047.01)

## 2026-11-24 — Moving out of Colab
- Cloned the GitHub repo to ~/code/exoplanet-search; working on branch `v2-vetting`.
- Added `pipeline.py` (Colab export of pipeline.ipynb + the v1 FP-test cells) and CLAUDE.md.
