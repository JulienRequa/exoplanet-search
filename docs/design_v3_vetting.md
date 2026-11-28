# Design notes: v3 vetting (a)–(f) — PROPOSALS, NOT BUILT, NOT RUN (2026-11-28, for Julien's review)

**v3's goal: two outputs from one run.**
- **(A)** Protect real planets and find planet candidates.
- **(B)** Classify every false positive (FP) by type and collect new eclipsing binaries (EBs).

v2 stays frozen. v3 never changes paper 1's numbers.

## Shared rules for all six notes
- **Data.**
  - Development uses only dev data: v1 FPs, the 13 known hosts, dev_v2, dev_v3, and the
    sibling-pool light curves for simulation backgrounds.
  - blind_v3 (data/blind_v3.csv, md5 1cfefce0…, locked 3fea7cb) is run once, at the end, with the
    frozen v3.
  - blind_v2 is not used, unless Julien decides to make it v3 development data (see decision D1).
- **Order for each check.** Physics argument → validation simulation with new seeds and
  pass/fail rules committed before running → dev-set check → freeze.
- **Statuses.** Every check reports exactly one of PASS / WEAK / FAIL / CANNOT_RUN, with a reason.
  Nothing is skipped silently.
- **Noise model in every simulation.**
  - White noise 3–5 ppt.
  - AR(1) red noise at the measured dev distribution of beta_2h: each simulated star's alpha is
    drawn as in sims/two_tier_threshold_sim.py.
  - 2 and 4 sectors, from the same sector templates as before.
- **Significance.** Every sigma below is measured against an **empirical null**: the same
  statistic computed at ~200 fake event times away from transit and phase 0.5 on the same star.
  This is the M2 / report-card method, so red noise and variability are inside the error bar.
- **Multiple checks.** About 8 checks are applied to every planet. If each falsely fails a
  planet with probability p, about 8p of planets are lost. So a "strong" failure is set at a
  false-failure rate of about 0.1% per check (≈ 5σ for Gaussian nulls). Then the planets lost
  to strong failures stay ≤ ~1%, which fits the 95%-kept bar.

---

## (a) Odd/even depth test for variable stars

**Physics.**
- An EB found at twice its true period shows alternating eclipse depths, unless its two stars
  are identical. A planet's transits are all the same depth.
- On a spotted star, the brightness at each transit differs. If the flattening (detrending)
  window is long compared with the spot timescale, each transit's depth inherits a baseline
  error.
- If the rotation period is near 2P (or another commensurate ratio), these errors alternate
  between odd and even transits. That fakes an odd/even difference.
- v2's red-noise fix (M2) inflates the errors but still measures depths against a global
  baseline. The report cards measured each transit against its own local baseline instead.

**Method.**
- **Per-transit depth.** Fit a local polynomial baseline to 1–3 × T14 on each side of each fully
  covered transit (order 1, or order 2 if it improves χ² by more than 10). Measure the depth
  against it, with the transit shape fixed from the global fit.
- **Statistic.** D = mean(odd) − mean(even), with significance from the empirical null. The null
  uses the same local fits at fake event pairs spaced P apart.
- **Timing drift.** Re-centre on the fitted multi-year ephemeris (M4), since timing drift was
  the main cause of false odd/even rejections (R2).
- **Variability flag (reported, not a failure).** A Lomb–Scargle peak at P, 2P or P/2 with
  amplitude ≥ 0.3 × depth sets "variability near orbital period". This feeds (f).

**Validation simulation.**
- **Stars:**
  - spotted stars: rotation period 0.5–20 d, amplitude 0.1–3%, spots that evolve (sum of
    sinusoids whose phase drifts on 1–3 rotation periods);
  - a resonance subset with rotation period = 2P ± 5%;
  - quiet stars;
  - noise as in the shared rules.
- **Signals:** 2,000 planets (equal depths); 2,000 EBs at 2P with odd/even depth differences of
  5–30%.
- **Reference:** v2's odd/even test run on the same light curves.

**Pass/fail rules (pre-registered).**
1. Planets falsely FAIL at ≤ 1% overall, ≤ 3% in every rotation-period bin, and ≤ 3% in the
   resonance subset.
2. Quiet stars: catches EBs with ≥ 10% depth difference at a rate within 5 percentage points (pp)
   of v2's.
3. Spotted stars: catches ≥ 50% of EBs with ≥ 20% depth difference. This reports an honest loss of
   power instead of false catches.
4. Dev set: no dev known planet newly lost. The change in dev FP catches is reported.

PASS requires 1–4.

---

## (b) Secondary-eclipse check robust to stellar variability

**Physics.**
- A planet's secondary-eclipse depth is at most about (Rp/Rs)² × (planet/star surface
  brightness), from thermal emission plus reflected light. That is usually ≤ a few hundred ppm.
  A deeper eclipse at phase 0.5 means a self-luminous companion.
- The upper limit, `planet_limit` in scripts/caveat_checks.py, uses albedo 1 plus full
  redistribution-free thermal emission, so it is generous.
- Phase curves (reflection at P, ellipsoidal variation at P/2, beaming) and spots change the
  flux near phase 0.5 relative to a wide baseline. That fakes a secondary. In the report cards,
  9 of 13 v2 secondary flags disappeared against a local baseline (PAPER_NOTES L5).

**Method (report-card method, made a vetting step).**
- **Depth.** A window of width T14 at phase 0.5, measured against flanks 1–3 × T14 on both sides,
  with a local linear baseline, event by event, then averaged.
- **Significance.** From the empirical null.
- **Eccentric orbits.** An optional scan of phases 0.35–0.65. When used, the same scan is applied
  to every null window, so the extra trials are paid for in the null.
- **Two outcomes.**
  - **Detected:** ≥ 3σ.
  - **Too deep for a planet:** depth − 2σ > planet_limit. This is what makes FAIL.
- **Ellipsoidal amplitude.** Measured and reported (feeds (f) and the companion-mass estimate).

**Validation simulation.**
- **Planets:** 2,000, each with a phase curve (reflection 0–300 ppm, ellipsoidal 0–200 ppm,
  beaming 0–20 ppm) and a true secondary at or below the planet limit, on quiet and spotted
  stars.
- **EBs:** 2,000 with secondary depths of 300–20,000 ppm, some eccentric (e up to 0.3).
- **Reference:** v2's secondary test on the same stars.

**Pass/fail rules.**
1. Planets FAIL at ≤ 1% overall, and ≤ 2% for planets with phase curves of 200–500 ppm.
2. EBs whose secondary is ≥ 2 × planet_limit and ≥ 5σ in white noise are caught at ≥ 90%.
3. Eccentric EBs: catch rate reported, not a gate.
4. Dev set: no dev known planet newly lost. The change in dev FP catches is reported. The hot
   Jupiters with real measured secondaries (e.g. WASP-18, KELT-9, WASP-33) must PASS, not FAIL:
   their measured secondary is below planet_limit.

PASS requires 1, 2 and 4.

---

## (c) The pixel (difference-image) check as a full vetting step, 4″ floor

**Physics.**
- In-transit minus out-of-transit images cancel every constant star. What is left is the
  point-spread function (PSF) of the source that dimmed.
- On dev data, 95% of 79 detected known planets were within **0.19 px (≈ 4″)** of the target
  (M9). The formal errors were too small, so the floor is the honest accuracy.

**Method.**
- **Per sector:** a difference image from in-transit cadences (the deep core) versus local
  out-of-transit cadences at 1–3 × T14. A Gaussian or PRF fit locates the dimming source. The
  target position is Gaia, moved to the sector epoch with proper motion (M5). A sector is
  usable if the difference-image S/N is ≥ 5.
- **Offset per sector:**
  - "off target" if the offset is > 0.19 px and > 3σ_stat;
  - "on target" if the offset is ≤ 0.19 px.
- **Statuses:**
  - **FAIL:** off target in ≥ 2 usable sectors, or in 1 of 1 usable sector with offset ≥ 0.5 px.
  - **WEAK:** off target in 1 sector out of ≥ 2.
  - **PASS:** on target in every usable sector.
  - **CANNOT_RUN:** no usable sector. The reason is given: S/N < 5, no TPF, or a crowded stamp.
- **Source match (feeds (f)).** Match the measured source position to Gaia neighbours: the
  nearest neighbour whose brightness allows the observed depth, using the neighbour-depth
  relation δ_n = δ_obs / f_n ≤ 1.
- **Data cost.** About 10–50 MB of pixel files per star per sector. These downloads are
  separate and need Julien's OK for each batch.

**Validation simulation.**
- **Synthetic pixel stamps:** 11×11 px, PRF-like PSFs, Gaia-like neighbour fields, pointing
  jitter of 0.01–0.05 px, and background noise.
- **Signals:** 1,000 planets on the target; 1,000 neighbour EBs (NEBs) at separations of
  0.1–2 px and brightness ratios ΔT of 0–6 mag.

**Pass/fail rules.**
1. Planets on the target FAIL at ≤ 1%.
2. NEBs at ≥ 0.4 px with difference-image S/N ≥ 10 are caught at ≥ 85%.
3. Source match: the correct neighbour is picked in ≥ 80% of caught NEBs.
4. Dev re-check (not independent, since the floor came from these stars; reported only): known
   planets FAIL at ≤ 3%, and the NEB controls are caught at ≥ 80%.

PASS requires 1–3. The independent test is blind_v3.

---

## (d) Gaia upgrades: RUWE, binary catalogues, Gaia radius and mass, and explicit CANNOT_RUN

**Physics.**
- **RUWE.** RUWE > 1.4 means the single-star astrometric model fits poorly. That usually
  points to an unresolved companion or an orbital wobble (Lindegren et al. 2021; Belokurov et
  al. 2020).
  - It does not prove the signal is an EB: many planet hosts have wide companions.
  - So RUWE is a **flag**, used by (e) and (f). It is never a FAIL on its own.
- **Gaia DR3 binary catalogues.**
  - Two-body orbit solutions (`nss_two_body_orbit`: SB1, SB2, astrometric).
  - Eclipsing binaries (`vari_eclipsing_binary`, Mowlavi et al. 2023).
  - An EB solution with a period within 1% of P, 2P or P/2 is direct evidence of an EB at the
    target position.
  - An SB1/SB2 orbit with a period near P, and a semi-amplitude K too large for a planet (K >
    1 km/s for P < 10 d), is the same evidence.
- **Radius and mass.** The density and size checks need stellar values. 16% of tested hunt
  targets lacked them (v3_list).
  - Fallback order: TIC → Gaia FLAME (radius_flame, mass_flame) → Gaia GSP-Phot radius with a
    mass from the main-sequence relation.
  - The source used is recorded for every star.

**Method.**
- Query once per star and save a dated snapshot of each catalogue result (CLAUDE.md rule 6).
- **Statuses:**
  - CANNOT_RUN, with the reason, for: no Gaia match within 2″, Gaia values flagged as
    unreliable, or no source for R or M.
  - Density and size each say which source was used.
- **Consistency flag:** when both TIC and Gaia exist and |ln(R_Gaia/R_TIC)| > 0.3.

**Validation (catalogue checks, not a light-curve simulation).**
1. On dev stars that have both TIC and Gaia radii: median |ln ratio| ≤ 0.10, and 90th
   percentile ≤ 0.30.
2. With the fallback on, dev planets lost to density or size are ≤ v2's dev count + 1.
3. Coverage: the share of dev stars where density can run rises. This is reported (expected
   84% → > 95%).
4. RUWE and catalogue flags: the rates among dev CP/KP and dev FP are reported. A flag is useful
   if it is ≥ 2× more common among FPs. That is reported, not a gate.
5. Unit test: every check returns one of the four statuses on every dev star (no blanks).

PASS requires 1, 2 and 5.

---

## (e) Three-way verdict: planet candidate / needs follow-up / likely false positive

**Why.** v2 calls any single FAIL a false positive. That merges strong evidence (a 10σ
secondary) with marginal evidence (a 3σ odd/even difference on a noisy star). The marginal cases
cost real planets: the 8 blind_v2 losses were mostly 3–5σ.

**Rule.** Strong-failure thresholds are set by the multiple-checks argument in the shared rules
(about 0.1% false rate per check), using each check's simulated null distribution.
- **Likely false positive:** at least one strong failure. Strong failures are:
  - secondary too deep for a planet at ≥ 5σ;
  - odd/even ≥ 5σ;
  - pixel FAIL in ≥ 2 sectors;
  - a Gaia EB or SB orbit at P, 2P or P/2;
  - size > 2.5 R_Jup with b well measured (2σ upper limit on b < 0.9);
  - density mismatch ≥ 5σ with stellar values from TIC or FLAME.
- **Planet candidate:** no FAIL or WEAK. The three key checks (secondary, odd/even, pixel) all
  ran. The star is in the normal or strict noise tier.
- **Needs follow-up:** everything else. Examples:
  - 3–5σ failures;
  - a single off-target sector;
  - a key check that could not run;
  - RUWE > 1.4 with otherwise clean checks;
  - tier-3 (too-noisy) stars, whose signals are listed as low confidence, never claimed.

**Validation.**
- **Simulation:** a combined population. Planets, EBs, blended EBs and NEBs from (a)–(c), plus
  variable stars and systematics, all through every check. Pass rules:
  1. Planets → likely FP at ≤ 2%.
  2. Planets → planet candidate at ≥ 80%.
  3. EBs with strong signatures → likely FP at ≥ 85%.
- **Dev:**
  4. Dev known planets → likely FP at ≤ 3% (v2 dev: 2.4% lost).
  5. Dev FPs → likely FP at a rate ≥ v2's dev catch rate − 5 pp.
  6. "Needs follow-up" ≤ 35% of all dev stars, so the middle class doesn't swallow everything.

PASS requires 1–6.

**blind_v3 metrics (pre-registered now, reported as run).** The share of planets and of FPs in
each of the three classes, with binomial 95% intervals.

---

## (f) False-positive type labels and the new-EB list

**Labels** are assigned only to likely-FP and needs-follow-up stars. The rules apply in order,
and the first one that matches wins:
1. **Neighbour EB.** Pixel check off target, with a matched Gaia neighbour; or a catalogued EB
   or variable within 2′ whose period is within 1% of P, 2P or P/2 and which is bright enough
   (the report-card NEB timing test).
2. **EB on target.** Dimming on target, plus a strong EB signature: a too-deep secondary,
   odd/even, a Gaia EB or SB orbit at the period, size > 2.5 R_Jup, or ellipsoidal amplitude
   implying a stellar-mass companion.
3. **Blended EB.** Dimming on target within 4″, but the transit shape or density says the eclipsed
   star is not the target. Examples: transit density inconsistent with the host by ≥ 3σ; a
   secondary too bright for a companion of the size the transit allows (the 3265.01 argument); a
   known close companion (RUWE, or high-resolution imaging on ExoFOP).
4. **Variability.** No robust eclipse signature, but a Lomb–Scargle peak at P, 2P or P/2 with
   amplitude ≥ the depth, or local-baseline tests that collapse versus wide-baseline tests.
5. **Systematic.** The signal is in a single sector or a single orbit, clusters at
   momentum-dump times or sector edges, or has a period within 2% of a TESS orbit harmonic.
6. **Unclassified.** Nothing above matches.

**Measured and reported for every labelled star:**
- period, from the fit;
- secondary depth ± error, from (b), with its phase;
- odd/even depths;
- estimated companion size and temperature. The companion temperature comes from the
  secondary/primary surface-brightness ratio, using the Planck ratio in the TESS band (as in
  scripts/review_checks.py sb_ratio). The radius comes from that temperature through a
  main-sequence table (Pecaut & Mamajek). The transit-depth radius is given alongside, with a
  dilution correction if the source is a neighbour.

**New-EB list.**
- Every EB-labelled star (labels 1–3) is cross-matched within 2′ against VSX (B/vsx), Gaia DR3
  eclipsing binaries (I/358/veb) and the TESS EB catalogue (Prša et al. 2022, J/ApJS/258/16),
  through the VizieR mirror.
- A catalogue period within 1% of P, 2P or P/2 marks the star "known EB".
- No match makes it "possibly new EB". The list goes to output (B), with the measured
  quantities above.
- Catalogue queries are saved as dated snapshots.

**Validation.**
- **Simulation.** A labelled population: on-target EBs, blended EBs (an EB diluted by a brighter
  star in the same pixel), NEBs (pixel level, from (c)), spotted stars with no eclipse, and
  systematics (single-sector ramps, momentum-dump dips). 1,000 of each. Pass rules:
  1. Correct label for ≥ 80% of NEBs and ≥ 80% of on-target EBs.
  2. ≥ 60% of blended EBs, variables and systematics.
  3. Companion temperature within 20% of truth for ≥ 70% of EBs with a ≥ 5σ secondary.
- **Unit test.**
  4. No star that matches a catalogued EB within 2′ at the period is ever listed "possibly
     new". Checked with injected catalogue rows.
- **Dev.**
  5. Dev FPs whose ExoFOP comments name a type ("NEB", "EB", "SB2", "BEB", "systematic",
     "variable") are used as labels. Label accuracy ≥ 60%. These comments are noisy labels, so
     the bar is modest. The comments are read for dev stars only.

PASS requires 1–5.

---

## Decisions for Julien before any of this is built
- **D1. blind_v2's status.**
  - Keep it untouched forever, or make it v3 dev data so its 8 planet losses can be studied?
  - The losses were mostly marginal odd/even (5) and density (2), exactly what (a), (d) and (e)
    target.
  - blind_v3 is now locked, so converting blind_v2 is allowed by the 2026-11-25 rule. Either way
    it is never used for accuracy claims.
- **D2. Build order.** Proposed: (d) first (catalogue work, cheap), then (b), (a), (c), (e), (f).
  (e) and (f) depend on the others.
- **D3. The strong-failure level.** About 0.1% false rate per check (≈ 5σ), from the
  multiple-checks argument. Or another level.
- **D4. Pixel data.** (c) needs pixel files for the dev stars, and later for blind_v3 and the
  hunt: about 10–50 MB per star per sector. Downloads in batches, with OK.
