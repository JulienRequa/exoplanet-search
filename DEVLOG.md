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
- Simulations (sims/oddeven_sim.py; 2 TESS-like sectors, Rp/Rs 0.1, 1000 ppm white noise;
  red noise = AR(1), 1 h correlation, same amplitude):
  | case | current FAIL | v3 FAIL |
  | planet, white (60) | 0% | 0% |
  | planet, red (60) | **57%** | 0% |
  | EB odd 15% shallower, white (30) | 100% | 100% |
  | EB odd 15% shallower, red (30) | 100% | 63% |
  → the current odd/even rule assumes white noise and falsely fails most planets on stars
  with correlated noise. v3 fixes this; the EB power lost in red noise is real (the data
  can't prove a 15% difference there), not a tuning effect.
- Density: the circular-orbit assumption can bias ρ by [(1 ± e)/√(1 − e²)]³ (0.37–2.7 for
  e = 0.32; 0.13–8 for e = 0.6), which covers the whole FAIL region (outside 0.5–2×).
  Baseline: density caught 0 of the 45 testable FPs on its own or otherwise.

### Dev-set runs after each fix (93 targets; each run from a clean checkout of its commit)
| run | planets kept | FPs caught | notes |
| baseline (35797d7) | 28/32 (88%) | 19/45 (42%) | 1 crash (TOI 734.01) |
| 2a grazing profile (ab9d789) | 28/32 | 19/45 | WARN on 3 FPs (all already caught), 0 planets |
| 2b maxfev 20000 (8d7b7b1) | 28/32 | 19/46 (41%) | 0 crashes; 80 converged fits bit-identical (max diff 0) |
| 2c Claret priors (51edd55) | 29/32 (91%) | 22/46 (48%) | density catches 3 FPs (2 only by density); grazing WARN on 6 FPs, 0 testable planets |
- Claret coverage (2c): 70 of 93 targets got a prior; 10 free (no Teff/log g); 1 free (outside
  table); 12 not fitted (no 2-min data).
- CORRECTION to the TOI-2416 b diagnosis: with priors b moves 0 → 0.64, a/Rs 17.4 → 13.6
  (published 14.4) and density PASSES (1.32, 0.5σ). The main cause was the b–a/Rs
  degeneracy with b pinned at 0 under free limb darkening; eccentricity (e = 0.32) is a
  smaller contributor. Proposal "density → WARN" WITHDRAWN (it would now cost 2 catches).
- pi Men did NOT stabilise with the prior: a/Rs ± 2.7 → ± 5.2, density still WEAK
  (b ≈ 0.02 at its bound; shallow transit).
- Remaining planet losses with 2c: HATS-55 b, TOI-2364 b (odd/even v3 proposal would rescue
  both), TOI-2207 b (needs the ephemeris fix (e)).

## 2026-11-24 — Odd/even v3 implemented (approved)
- pipeline.odd_even = prototype v3 exactly (β from 300 out-of-transit fake transits, seed 0;
  only fully covered transits). Pipeline σ = prototype σ to 4 decimals. Simulation script and
  results kept in sims/ (oddeven_sim.py, oddeven_sim_results_2026-11-24.md).
- **Honest loss:** TOI 2554.01 (TIC 144043410, v1 FP) was caught only by odd/even (9.5σ);
  v3 gives 2.2σ (β = 4.5: strong correlated noise; its Comment mentions stellar variability
  on short timescales). It is now expected to pass as a candidate. We accept this: the old
  9.5σ came from treating correlated noise as white.
- **Known limitation:** pi Men's density stays WEAK even with the LD prior (a/Rs ± 5.2):
  very shallow transit (Rp/Rs ≈ 0.017) and b stuck at 0, so a/Rs (and ρ ∝ (a/Rs)³) are
  poorly constrained. Density stays a FAIL-capable check (density → WARN withdrawn).

## 2026-11-24 — Designs + simulations for the remaining fixes (none in the pipeline yet)
- Ephemeris drift (docs/design_ephemeris.md, sims/ephemeris_sim*.py/md): per-sector transit
  times + straight-line ephemeris. Hard TOI-2207-like case (P 6–9 d, Rp/Rs 0.067, red noise,
  1100 d gap): sector-2 transits lost — current 40%, per-sector 0%. Easy case: 0% / 0%.
- Period window (docs/design_period_window.md, sims/period_window_sim*.py/md): pre-registered
  rule "1× unless alt/1× BLS power > R*", R* from simulated spotted-star planets, valid only if
  R* < 2 (a genuine alternative gives ≈ 2). Result R* = 3.73 → always use 1× with a hint.
  Spotted planets kept at 1×: current 75%, proposed 100%. EBs with true P = 2× catalog:
  caught 70% / 50% (secondary 60% / 90% depth) under both rules. Not tuned on HIP 67522 b.
- Difference-image centroid prototype (sims/diffimage_sim*.py/md, simulated pixels):
  first error estimate (per-transit scatter as exact χ²) gave 4–8% null false alarms →
  fixed (Hotelling T²; or stacked image with empirical pixel noise = method B). Method B in a
  hard case (Tmag 13, 0.2% dip, 3 transits, PSF 15% wrong): neighbour at 0.5 px detected
  96–100%, 0.3 px 34%, 0% false alarms; ~10% distance bias from the PSF error. Idealised:
  no jitter, PRF variation or WCS/Gaia error — the real-data systematic floor is unknown.

## 2026-11-24 (evening) — "dev sample 2" (data/dev_v3.csv)
- 100 FP + 100 CP/KP, random_state=4, same cuts; excludes blind_v2, v1, dev_v2 and the 13 known
  hosts by TIC (overlaps: 0). Purpose: FINDING PROBLEMS before the freeze — not for reporting.
- blind_v2 untouched: md5 unchanged, only commit = lock commit 0217ed9, none of its stars in
  any results file.
- Unused matching stars after dev_v3: 367 FP, 506 CP/KP (enough for a ~100 FP v3 blind sample).
- Code for the dev_v3 run: odd/even v3 + Claret priors (commit after 927087a); period-window and
  ephemeris fixes NOT included (known failure modes).

## 2026-11-24 (evening) — Decisions (Julien)
**v2 scope, final:** odd/even v3 (in); per-sector ephemeris fix (if its simulations meet the
criteria: injected planets at the correct period, EBs at 2× catalog still caught, 13-planet
regression unchanged); Gaia neighbour-depth check (validated on simulations; must not flag
real planets, including crowded fields like HATS-55 b). Then ONE full rerun of all dev sets
→ freeze → push → blind_v2. No other additions to v2.

**Period window: not in v2** — it's mostly a detrending problem. Known v2 limitation:
HIP 67522 b and other very active stars. With a ~1-day detrending window, spot residuals on
fast rotators (P_rot 0.5–5 d, 0.5–3% amplitude) can be larger than the transit, and BLS locks
onto them even inside the ±2% 1× window. Simulation (60 spotted planets): fitted period
within 1% of the truth — current max-power rule 43%, always-1× 53%.

**v3 plan (after the blind_v2 run, with its own fresh blind sample), in order:**
1. Difference-image centroid — go/no-go first: measure the real-data positional floor on the
   known planets (~½ day). Build the full check only if the floor is below ~0.3 px;
   otherwise stop and log the result.
2. Detrending for fast-rotating spotted stars (own design note + simulations; needs a full
   dev rerun because it changes cleaning for every star).
3. Per-transit depth, only if time allows (needs red-noise errors and crowding handling).
v2 hunting on open PC TOIs can start right after the blind run.

**GitHub:** gh CLI from the official release into ~/.local/bin; browser login. Check the
no-reply format; push v2-vetting as-is tonight (with OK) to publicly timestamp that blind_v2
was locked (commit 0217ed9) before the v2 improvements. Later pushes only at: the freeze
(tag v2-freeze), right after the blind run, and before hunting — each with Julien's OK.

## 2026-11-24 — GitHub
- gh 2.105.0 (official macOS arm64 release, SHA-256 verified against the release checksums,
  code-signed) in ~/.local/bin; config in ~/.local/share/gh because ~/.config is owned by
  root (fish install, Jul 2024). Logged in as JulienRequa via browser (device code).
- Account created 2026-11-16 → no-reply format 357599494+JulienRequa@users.noreply.github.com,
  set for FUTURE commits in this repo. Existing 38 commits keep JulienRequa@users.noreply.github.com
  (private, not rewritten, so 0217ed9 and all DEVLOG hashes stay valid).
- Pushed v2-vetting at dce7e0b (2026-11-24 17:56 PST): public record that blind_v2 was locked
  (0217ed9) before the v2 improvements. Next pushes only at: freeze (tag v2-freeze), right after
  the blind run, before hunting — each with Julien's OK.

## 2026-11-24 (evening) — Dev set (93) with odd/even v3 (commit 927087a + priors)
- Planets: 30/32 kept (94%) — up from 29/32. Rescued: TOI-2364 b (partial transit dropped).
  Still lost: TOI-2207 b (ephemeris drift → per-sector ephemeris fix) and HATS-55 b
  (odd/even 3.2σ; β ≈ 1.15 so red noise is small; likely cause (c): sector 61 is 29% deeper
  than sector 34 in a very crowded aperture, and with a partial transit dropped the odd and
  even groups are no longer balanced across sectors, so the sector depth difference leaks into
  odd/even). Possible physical fix for discussion: compare odd vs even within each sector.
  The prototype's 2.6σ for HATS-55 b was measured with the free-LD fit; with priors it is 3.2σ.
- FPs: 21/46 caught (46%), 12 of 60 no 2-min data, 2 period mismatches. Only change vs the
  priors run: TOI 2554.01 caught → missed (odd/even 9.6σ → 2.2σ), as predicted.
- Two planets (TOI 255.01, TOI 115.01) first hit MAST download timeouts (600 s) while several
  jobs ran at once; rerun → both PLANET CANDIDATE. Not code crashes; network.
