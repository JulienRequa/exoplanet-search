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

## 2026-11-24 — Environment, baseline, bug fixes, regression
- conda env `exoplanet` (environment.yml; exact versions in environment.lock.yml).
- pipeline.py made importable (no function changed). `regression.py run/compare` added.
- Local baseline (pre-fix code): 14/14 PLANET CANDIDATE, all periods match; WASP-18 and
  WASP-121 match Colab within 0.25σ. pi Men density PASS (σρ/ρ = 41%) vs WEAK in Colab.
- Bug fixes committed separately: (1) Rp/Rs cap 1.0 + guess clipping, also in profile_b;
  (2) TIC fallback for missing/NaN Teff, log g, radius (`star_source` recorded);
  (3) `period_hint` → `search_near` (±2% of 1×/½×/2×, highest BLS power, `period_choice`).
- Regression after fixes: 13/13 recovered, every period matches; 12 planets 5/5, pi Men 4/4
  with density WEAK — identical to the Colab regression. Largest change in any fitted value
  that is not at a bound: 0.015σ. WASP-18 (4 sectors) unchanged.
- Note for the paper: pi Men's fit has b, u1, u2 at bounds, so its covariance is nearly
  singular. Best-fit values don't move, but the a/Rs error swings (±2.2 → ±6.3) with
  microscopic changes in where the fit lands; its density WEAK/PASS status is fragile.
  Bug fix 1's side effect also flags pi Men's Rp/Rs (0.0166, within 0.0099 of the 0.01 floor).

## 2026-11-24 — Grazing check, v1 record, GitHub sync
- Grazing check + three-level verdict added. Regression: 13/13, all fitted values identical,
  no grazing warnings on known planets (b + Rp/Rs between 0.02 and 0.85).
- **Grazing check has almost no power** (simulations, 3 seeds × 2 noise levels): 0/6
  simulated grazing transits (b = 1.02, Rp/Rs = 0.15) were flagged, even at 200 ppm. With free
  limb darkening, a high-b non-grazing solution fits as well or better (χ² 19143.5 vs 19146.7
  for the truth), sometimes with confident-looking but wrong errors (b+k = 0.993 ± 0.010).
  Kept as tests/test_fixes.py::test_grazing_transit_flagged (xfail). Decide before v2 freeze.
- results/fp_test_v1.csv committed as the official v1 record (from Colab).
- Pulled final Colab pipeline.ipynb (commit 3ed189d): only adds the FP-test v1 cells, no
  function changes. Its saved outputs record the v1 TOI table: 8175 TOIs, 627 FPs passing
  the cuts. fp_test_v1() in pipeline.py now also merges in the Comments, as the notebook does.
- WD 1856+534 b (white-dwarf host, P = 1.4079 d, hint → 1×): Rp/Rs hit the 1.0 cap; TIC has
  no Teff/log g/radius (all NaN) → size and density N/A; verdict PLANET CANDIDATE (some checks
  weak). Note: fit gave u1 = u2 = 1.0 (u1 + u2 = 2, physically impossible limb darkening →
  model depth 103.9%); current LD bounds allow it — for Julien to decide.
- Centroid-check design written: docs/centroid_check_design.md (awaiting approval).

## 2026-11-24 — Limb darkening (q1/q2 + priors) and centroid check
- Decisions (Julien): Kipping q1, q2 + Claret Gaussian priors (σ = 0.1); size lower-limit
  rule approved; centroid check approved (FAIL if shift > 3σ AND source > 1 px).
- Implemented q1/q2 in fit_transit and profile_b. Priors enter curve_fit as two extra
  residuals (u − u_table)/σ, flux points weighted by their scatter; profile_b's χ² includes
  them. Claret (2017) TESS table NOT yet downloaded: CDS (VizieR/FTP/TAP) was down all
  session → every star currently runs with free q1, q2 (`ld_source` = "free (no LD table)").
- Regression (q1/q2, no prior): 13/13, largest change 0.9σ (WASP-33: old fit had
  unphysical u2 = −0.33). pi Men: a/Rs error 6.3 → 2.5, density WEAK → PASS; still at an LD
  corner (u1 = 0, u2 = 0.99) — needs the prior to judge whether it stabilizes.
- Grazing simulations (b + k = 1.17; 2 noise levels × 3 seeds):
  covariance σ (approved rule): free 0/6, prior 0/6 (best +2.5σ).
  Profile-likelihood σ (refit with b + k ≤ 1, σ = √Δχ²; prototype, NOT in pipeline):
  free 2/6, prior 4/6 (all three at 200 ppm, 5–7σ), table off by 0.1: same 4/6.
  No false alarms in any setup (b = 0.3). → proposal for Julien.
- Centroid check built as approved (light-curve centroid columns, null from 200 fake
  transits, source offset vs the target's WCS position). WCS-vs-centroid systematic ≈ 0.3–0.4 px.
  Dev FPs with neighbour Comments: shifts 2.8σ, 1.5σ, 3.2σ, 7.8σ, but the named neighbours
  are at 0.60, 0.47, 1.22, ~0.5–0.6 px → 3 of 4 are inside the 1 px limit by construction;
  the 1.22 px one is measured at 0.53 px (aperture truncation biases moment centroids low).

## 2026-11-24 — Dev set (13 known planets + 20 v1 FPs), draft v2 rules, no LD prior yet
DEVELOPMENT ONLY — v1 is not re-scored as an improvement. All runs with period_hint.
- Known planets: 13/13 kept (12 PLANET CANDIDATE, pi Men PC with Density + Centroid WEAK).
  No grazing WARN, no centroid FAIL. Centroid shifts on planets are often highly significant
  (up to 13σ, HD 189733) but the implied source stays within 0.39 px of the target.
- v1 FPs: 3 no 2-min data, 1 error (TIC 146846569, 23% deep: curve_fit max evaluations
  exceeded), 1 period mismatch (TIC 441738827), 15 testable → 6 caught (40%).
  FAILs: odd/even 3, size 3, SNR 1, secondary 1, centroid 1 (TIC 269450900, 1.24 px).
- Density is WEAK for 9 of 15 testable FPs (free LD → large a/Rs errors): the LD prior
  should matter a lot here.
- 4 neighbour-contaminated dev FPs: none FAIL the centroid check (3 neighbours < 1 px;
  measured 0.28–0.66 px vs planets ≤ 0.39 px). Two V-shaped FPs: no grazing WARN.
- TIC has no log g / radius for 4 FPs (so size/density N/A) — the fallback worked; the TIC
  simply lacks them. TOI 2438.01 (Comment: Jupiter around a WD?) hits Rp/Rs = 1.0, no radius.
- data/fp_sample_v1.csv: v1 sample TOI periods/depths from pipeline.ipynb's saved output.

## 2026-11-24 — Plan change: improve v2 before the freeze
- Blind sample LOCKED first (data/blind_v2.csv, 100 FP + 100 CP/KP, random_state=2, from
  data/toi_2026-11-24.csv; excludes v1 + 13 known hosts by TIC; not run or inspected).
- Dev extension data/dev_v2.csv (40 FP + 20 CP/KP, random_state=3; excludes blind + v1).
  Full dev set = 93 targets (13 known + 20 v1 FP + 60 new).
- Baseline (code 35797d7, no LD prior): planets 28/32 kept, FPs 19/45 caught (42%).
  4 dev planets lost: 3 to odd/even, 1 to density. 32 testable planets (1 period mismatch).
- Fixes committed separately: (2a) grazing → profile likelihood; (2b) curve_fit
  maxfev 20000; (2c) Claret priors (VizieR mirror vizier.cfa.harvard.edu; CDS down).
- Design notes (not built): docs/design_difference_image.md, design_neighbor_depth.md,
  design_per_transit_depth.md.
- Mac note: the "restarts" were automatic log-outs after 5 min idle (setting now off).
  dev_set.py now saves after every star and resumes.

### Diagnosis plan for the 4 lost planets (before any fix)
1. Per planet: detected vs catalog period; odd and even depths ± errors; number of transits
   in each group; sectors used; for the density loss, which numbers drove it. Explain why
   32 planets were testable, not 33.
2. Test causes: (a) period alias (½× / 2×); (b) underestimated depth errors (odd_even's
   white-noise error vs transit-to-transit scatter / red noise); (c) depth differences
   between sectors (crowding/dilution); (d) detrending distortion near data gaps.
   Report which cause explains each case.
3. Any fix: physically justified, checked on simulated planets, and its effect on the 19
   currently caught FPs reported. No threshold tuning.
4. Report how many of the 93 dev targets get Claret priors vs fall back to free LD.

### "Ready for blind test" bar (set by Julien, 2026-11-24)
- ≥ 95% of real dev planets kept, with every loss explained
- 0 crashes
- each new check behaves as its design note predicts on simulations
- FP catch rate clearly above the 42% baseline (no fixed target)

### Push plan (later, not yet)
After the freeze: set up GitHub login, confirm commits use the no-reply email, and push
the frozen rules BEFORE running the blind sample (public timestamp). Nothing pushed yet.

### Diagnosis results (baseline code; per-transit depths from the cleaned light curve)
- 32 not 33: HIP 67522 b (TIC 166527623, TOI 6551.01, KP, P = 6.959 d). With the hint, BLS chose
  the ½× window and found 3.545 d (1.9% off ½×, at the window edge) — not the planet. Very
  young, spotted star: variability beats the transit in BLS → period mismatch, untestable.
- HATS-55 b (TIC 114299824, odd/even 3.0σ) → cause (b) + (c). Single-transit depths scatter
  8.8× (reduced χ²) more than odd_even's white-noise errors; scatter-based error gives 1.3σ.
  Sector 61 is 29% deeper than sector 34 (19 600 vs 15 200 ppm) in a very crowded aperture
  (CROWDSAP 0.39); Rp/Rs 0.135 vs published 0.114. Parity is balanced per sector, so (c)
  adds scatter rather than creating the difference. Not (a).
- TOI-2207 b (TIC 90850770, 11.4σ) → none of (a)–(d): NEW cause (e) ephemeris drift. The
  two sector-67 "transits" have negative depth (no transit at the predicted times). Fitted
  P = 7.99905 d vs TOI 8.00049 vs published 8.00197 d: over 136 cycles (S27 → S67) that is
  4.7–9.5 h, so the transits are missed; both S67 epochs are "even". The multi-sector
  refit's period window (5× the 1-sector error) did not contain the true period. Side
  effect: S67's real transits were not masked during flattening (d).
- TOI-2364 b (TIC 39414571, 3.5σ, 1 sector) → cause (d). Epoch 0 sits at a gap edge
  (0.00 d) and has negative depth (−2769 ppm): a partial transit distorted by trim/ramp.
  The other four transits agree (8548/9783 odd, 9322/9747 even).
- TOI-2416 b (TIC 237922465, density 2.78×, 5.1σ) → eccentricity (photoeccentric effect).
  Published e = 0.32; the density check assumes a circular orbit. ρ_obs/ρ_true =
  [(1 + e sinω)/√(1 − e²)]³ reaches 2.7 for e = 0.32; our a/Rs 17.4 vs published 14.4.
  (Published values from the NASA Exoplanet Archive pscomppars, 2026-11-24.)

### Odd/even prototypes (NOT in pipeline)
- v2: β from scatter of single-transit depths within each parity group + drop partial
  transits → rescues all 3 planets but un-catches 7/9 odd/even FPs (β up to 55: within-group
  scatter absorbs ephemeris failures and real depth variation). REJECTED.
- v3: β from ~300 fake transits outside transit (signal cannot leak in) + drop transits not
  fully covered → FPs: 7/9 still FAIL (loses TOI 1127.01, still caught by 2 other checks, and
  TOI 2554.01, β = 4.5, Comment "stellar variability…" — net 19 → 18 caught);
  planets: HATS-55 b 2.6σ PASS, TOI-2364 b 0.7σ PASS, TOI-2207 b 5.8σ still FAIL (needs (e)).
