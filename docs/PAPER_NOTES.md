# Paper notes — running list of paper-worthy findings

Separate from DEVLOG (which records *what we did*). Each entry: the finding, key numbers,
where the evidence lives, and the paper section it fits. Newest at the bottom of each
section. Numbers from the **development sets are for method development only**; the
headline performance numbers will come from the frozen v2 blind test (blind_v2).

---

## Methods

### M1. Development / blind / freeze design (pre-registration)
**Finding.** Vetting rules are developed only on development samples; a blind sample is
drawn and locked *before* any v2 change, never run or inspected until the rules are frozen,
and the lock is publicly timestamped on GitHub.
**Numbers.** v1 official blind test: 20 TFOPWG FP TOIs (random_state=1). blind_v2: 100 FP +
100 CP/KP (random_state=2) from a dated ExoFOP snapshot (8175 TOIs, 2026-11-24); cuts
P = 0.5–10 d, depth > 2000 ppm; one TOI per star; excludes v1 and the 13 regression hosts.
Dev samples: dev_v2 (40 FP + 20 CP/KP, rs=3), dev_v3 "dev sample 2" (100 + 100, rs=4).
After dev_v3, 367 unused matching FPs remain for a v3 blind sample.
**Evidence.** data/toi_2026-11-24.csv, data/blind_v2.csv (lock commit 0217ed9, pushed in
v2-vetting at dce7e0b on 2026-11-24), sampling.py, data/dev_v2.csv, data/dev_v3.csv, DEVLOG.
**Section.** methods.

### M2. Red-noise-aware odd/even test
**Finding.** The standard odd/even depth test assumes white noise; on stars with correlated
("red") noise it falsely rejects most real planets. Inflating the errors by a red-noise
factor measured from ~300 fake transits *outside* transit (so the signal cannot leak in),
and using only fully covered transits, removes the false rejections while still catching
eclipsing binaries.
**Numbers.** Simulations (2 sectors, Rp/Rs 0.1, 1000 ppm white noise; red noise AR(1), 1 h,
1000 ppm): planets in red noise falsely failed **56.7% → 0%**; planets in white noise 0% → 0%;
EBs (odd eclipses 15% shallower) caught 100% → 100% (white), 100% → 63.3% (red; honest
loss of power). An earlier variant (β from scatter *between* transits) was rejected: it
hid 7 of 9 dev FPs. Dev set: rescues HATS-55 b and TOI-2364 b; costs one FP (TOI 2554.01,
correlated noise).
**Evidence.** commit 927087a; sims/oddeven_sim.py, sims/oddeven_sim_results_2026-11-24.md,
sims/oddeven_prototypes.py; DEVLOG "Diagnosis results".
**Section.** methods (+ results for the dev-set effect).

### M3. Physical limb darkening + theoretical priors un-stick the impact parameter
**Finding.** Fitting limb darkening freely lets b pin at 0, which inflates a/Rs and hence
the transit-derived stellar density ∝ (a/Rs)³, so real planets can fail the density check.
Kipping (2013) q1, q2 (physical profiles only) plus Gaussian priors from Claret (2017)
TESS tables (σ = 0.1) fix this.
**Numbers.** TOI-2416 b: b 0 → 0.64, a/Rs 17.4 → 13.6 (published 14.4), density
2.78× catalog (5.1σ, FAIL) → 1.32× (0.5σ, PASS). Dev set: planets kept 28/32 → 29/32;
FPs caught 19/46 → 22/46 (density becomes a working check). 70 of 93 dev targets got a
prior (10 no Teff/log g, 1 outside the table, 12 no data).
**Evidence.** commits 7b463b2 (q1/q2), 51edd55 (priors);
data/claret2017_tess_quadratic_2026-11-24.csv (ATLAS table25, Z = 0, ξ = 2 km/s, from the
VizieR mirror vizier.cfa.harvard.edu); results/devset93_baseline_2026-11-24.csv vs
results/devset93_c_priors_2026-11-24.csv.
**Section.** methods (+ results).

### M4. Ephemeris drift across multi-year sector gaps
**Finding.** Refining a period from one sector inside a narrow window fails when sectors
are years apart: small period errors accumulate to hours and later transits are missed
(then odd/even compares real transits with empty windows). Measuring a transit time per
sector and fitting one straight-line ephemeris fixes it.
**Numbers.** TOI-2207 b (P = 8.0 d, sectors 27 → 67, 136 orbits): fitted P off by 0.0014 d →
4.7–9.5 h drift, sector-67 transits missed, false odd/even FAIL. Simulation, hard case
(P 6–9 d, Rp/Rs 0.067, red noise, 1100 d gap): correct period 33% → 70% (prototype),
median |ΔP| **126 s → 0.7 s**; EBs at 2× catalog caught 95%/100% → 100%/100%.
Earlier stand-in: sector-2 transits lost 40% → 0%. (The remaining 30% still to be split
into "missed transits" vs "error bar too small" before the fix goes in.)
**Evidence.** docs/design_ephemeris.md; sims/ephemeris_sim.py + results (106b3f3);
prototype simulation sims/ephemeris_proto_sim.py (to be committed with the fix).
**Status.** prototype; enters v2 only if the criteria are met.
**Section.** methods.

### M5. Proper motion must be applied before comparing centroids with the target position
**Finding.** Using the catalog (J2000) position for fast-moving nearby stars makes the
target look offset from its own light, so the centroid check wrongly blames a neighbour.
**Numbers.** GJ 436 (1.21″/yr): moved 1.16 px by TESS sector 22; centroid − J2000 position
1.18 px vs centroid − propagated position 0.15 px. LTT 1445 A: 0.65 → 0.36 px. WASP-18
(control): 0.03 px of motion.
**Evidence.** DEVLOG 2026-11-24 (night); results/devset_v3sample_2026-11-24.csv.
**Status.** bug found in dev sample 2; fix pending approval.
**Section.** methods (+ a caution for other TESS vetting pipelines).

### M6. Period window with a catalog hint: always use 1×
**Finding.** When vetting a known TOI, searching only near the catalog period (not ½× or 2×)
keeps the analysis on the signal that was actually reported.
**Numbers.** Simulations: safe but neutral — quiet stars: final period correct 100% under both
rules, EBs at 2× catalog caught 100% under both; spotted stars: 43% → 53%. Real dev data:
rescues LTT 1445 A b and TOI-2266 (both had been analysed at 2× / ½×). Faint planets
(Kepler-44 b, TOI-4306, K2-415 b) are NOT rescued: too weak to find in one sector.
**Evidence.** sims/period_window_sim*.py/md, sims/period_window_ordinary_sim*.py/md;
results/devset_v3sample_window_rerun_2026-11-24.csv; DEVLOG.
**Section.** methods (honest: benefit shown on real data, not in simulation).

### M7. Multi-start fit against the b = 0 local minimum
**Finding.** Transit fits can stick at impact parameter b = 0 (b–a/Rs degeneracy), inflating
a/Rs and the transit density; LD priors (M3) reduce but do not remove it. Fitting from a
second start (b = 0.6) and keeping the lower-χ² solution fixes the remaining cases.
**Numbers.** HATS-28 b: b 0.00 → 0.57, a/Rs 12.65 → 10.37, density 2.38× (4.3σ, FAIL) →
1.31× (0.7σ, PASS). Same pattern earlier: TOI-2416 b, pi Men c. 13-planet regression: no
verdict change. Full dev rerun: planets kept unchanged (121/124); cost: Kepler-74 b (faint) jumps
to a grazing solution at the b bound (Rp 5 R_Jup → size FAIL), and 3 borderline FPs pass
(69 → 67/127 caught). Honest trade-off: multi-start finds the lower χ², which on noisy data can
be the grazing branch of the b–Rp/Rs degeneracy.
**Evidence.** commit e216703; DEVLOG 2026-11-24/25.
**Section.** methods.

---

## Results

### R1. Light-curve-only false-positive catch rate
**Finding.** With checks that use only the light curve and its centroids, a bit under half of
deep, short-period TESS false positives are caught.
**Numbers.** v1 official blind: 5/12 testable (≈42%). Development: 19/45 (42%) baseline →
22/46 (48%) with priors (dev, not a reported result). Planets kept on dev: 88% → 91%.
The v2 blind number will replace these.
**Evidence.** results/fp_test_v1.csv; results/devset93_*_2026-11-24.csv; dev_set.py summary.
**Section.** results.

### R2. Timing drift is the main cause of false odd/even rejections in multi-year data
**Finding.** In dev sample 2, most real planets rejected by odd/even were bright, clean hot
Jupiters whose sectors are years apart; a slightly wrong period misaligns later transits.
**Numbers.** 16 planets wrongly rejected; per-sector ephemeris rescues 8 (odd/even 3.3–9.7σ →
PASS), e.g. HATS-59 b (period 0.06% off). Remaining: 5 period-search, 2 proper-motion bug, 1 open.
**Evidence.** results/devset_v3sample_2026-11-24.csv;
results/devset_v3sample_lost_rerun_timingfix_2026-11-24.csv; DEVLOG.
**Section.** results (development) / methods motivation.

---

## Limitations

### L1. Neighbours closer than ~1 TESS pixel (21″)
**Finding.** Most contaminating neighbours in our dev sample are inside 1 pixel, where the
light-curve (moment) centroid cannot separate them; moment centroids also underestimate
offsets of sources at the aperture edge.
**Numbers.** 4 dev FPs with neighbour Comments: neighbours at 0.60, 0.47, 1.22, ~0.5–0.6 px
(3 of 4 < 1 px); measured source offsets 0.28, 0.09, 0.53, 0.66 px (the 1.22 px one at
0.53 px); centroid shifts 2.8σ, 1.5σ, 3.2σ, 7.8σ. Known planets: shifts up to 13σ but
source ≤ 0.39 px from target. Difference-image prototype on simulated pixels: 0.5 px
neighbours found 96–100%, 0.3 px 34%, 0% false alarms (idealised; real-data floor unknown).
**Evidence.** commit 3058d01 (centroid check); results/devset_v2draft_2026-11-24.csv;
sims/diffimage_sim.py + results (dfcfffb); DEVLOG.
**Section.** limitations (→ v3 difference-image work).

### L2. Very active, fast-rotating stars
**Finding.** With a ~1-day detrending window, spot residuals on fast rotators can be larger
than the transit, and BLS locks onto them even when restricted to ±2% of the catalog
period. Mostly a detrending problem, not a window-choice problem.
**Numbers.** 60 simulated spotted-star planets (P_rot 0.5–5 d, 0.5–3% amplitude): fitted
period within 1% of truth — max-power window 43%, always-1× 53%. Window choice: 75% → 100%
keep 1×. Real case: HIP 67522 b (untestable in dev, wrong period).
**Evidence.** docs/design_period_window.md; sims/period_window_sim.py + results (7f54085);
fitted-period simulation sims/period_window_fit_sim.py (prototype checkout, to be committed).
**Section.** limitations (→ v3 detrending).

### L3. Shallow transits with b at 0: weak density test
**Finding.** pi Men c (Rp/Rs ≈ 0.017): even with the LD prior, a/Rs is poorly constrained
(± 5.2), so the density check is WEAK.
**Evidence.** results/devset93_c_priors_2026-11-24.csv; DEVLOG.
**Section.** limitations.

### L4. Planets with strong transit-timing variations (TTVs)
**Finding.** The pipeline fits one straight-line ephemeris (fixed period). Planets whose
transit times wander (gravitational tugs from a neighbour planet) are not modelled, which can
misalign transits and push the odd/even test over its threshold.
**Numbers.** K2-19 b (near 3:2 resonance with K2-19 c; TTVs of hours over years): odd/even
3.0σ → FAIL in the full v2 dev rerun (honest threshold loss, threshold unchanged). Sectors
45/46 times agree with the linear ephemeris within errors, so the effect here is small/
borderline. Related: HATS-64 b also lost at 3.1σ (no TTVs known; statistical/systematic).
**Evidence.** results/devfull200_v2final_2026-11-24.csv; DEVLOG 2026-11-24/25.
**Section.** limitations (→ future: per-transit times, which also feed the sibling search).
