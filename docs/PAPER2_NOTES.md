# Paper 2 notes (sibling search with v3) — living document

Results and decisions that belong to paper 2 (v3, multi-sector / sibling search). Paper 1 (frozen
v2) notes are in docs/PAPER_NOTES.md; nothing here changes paper 1's numbers.

## R1. A red-noise-aware detection SNR (time-averaging beta) failed its pre-registered validation (2026-11-27)

**What was tried.** v3 candidate scripts/rednoise_snr.py: the time-averaging method (Pont, Zucker &
Queloz 2006; as implemented by Winn et al. 2008). Out-of-transit residuals (> 1 × T14 from
mid-transit) are binned at 0.75, 1.0 and 1.25 × T14; beta = measured bin scatter ÷ white-noise
expectation, median over the three bin lengths, never below 1. SNR_red = SNR_white / beta.

**Rules (committed before the run; sims/rednoise_snr_sim.py header).** On the 5,000 stars of the
multi-sector threshold simulation (2 and 4 sectors, white vs white + red noise, AR(1) τ 1 h,
amplitude 50% of white):
- V1, no penalty on white noise: median beta ≤ 1.10 and faint-planet recovery at the FAR-1%
  threshold within 3 percentage points of the white-noise SNR's.
- V2, threshold no longer depends on the noise colour: FAR-1% threshold with red noise within 25%
  of the white-noise one.
- V3, transfer: the white-noise FAR-1% threshold gives ≤ 3% false alarms on red-noise stars.
- PASS = all three; if FAIL, not adopted, and any change is a new committed version on new seeds.

**Verdict: FAIL** (sims/rednoise_snr_validation_results_2026-11-27.md). The white-noise SNR in this
run reproduced the earlier simulation exactly (5,000 stars, max difference 0).

| | 2 sectors | 4 sectors |
|---|---|---|
| V1: median beta on white noise; recovery change | 1.00; +0 pp → pass | 1.00; +0 pp → pass |
| V2: FAR-1% threshold white vs red (white-noise SNR was 8.4 vs 28.6, +240%) | 8.38 vs 15.30 (+83%) → FAIL | 8.15 vs 14.58 (+79%) → FAIL |
| V3: false alarms on red stars at the white threshold | 30.4% → FAIL | 36.9% → FAIL |

**What it fixed.** It removed about two-thirds of the red/white threshold gap (+240% → about +80%)
at no cost on white noise. At the measured pessimistic level (alpha 0.292) the FAR-1% threshold
fell from 16.3 / 17.2 to 11.6 / 11.3, but faint-planet recovery was not better (51% / 70% vs
60% / 74%), because beta also deflates the planets' SNR.

**Likely causes (not tested).**
1. The search keeps the deepest dip over many trial phases and periods (look-elsewhere effect).
   Red noise fattens the tail of that maximum more than it raises the typical bin scatter that
   beta measures, so a correction from average bins under-corrects the tail that sets a
   false-alarm threshold.
2. The estimator reads low on white noise: before the floor at 1, simulated white noise gives
   beta ≈ 0.90 at 2 h bins (sim calibration, results/rednoise_dev_2026-11-27/), so weak red
   noise is partly hidden under the floor. At the measured typical level (alpha 0.091) beta is
   1.00 on noise stars, yet the FAR-1% threshold still rises from 8.4 to 9.2.

**Consequence.** Not adopted. v3 keeps the white-noise-style SNR with a threshold calibrated on
simulated noise at the measured red-noise level (Julien's decision 2026-11-27: provisional FAR
1%, target 0.25%, pending a 10,000-star-per-setup run, sims/large_noise_threshold_sim.py). The
noisiest stars are to be handled by a two-tier threshold (docs/design_two_tier_threshold.md).

**For the paper.** Report as a negative methods result: a standard red-noise correction from the
transit-fitting literature is not sufficient to calibrate a search threshold, because a search
statistic is a maximum over many trials. Evidence: scripts/rednoise_snr.py,
sims/rednoise_snr_sim.py, results/rednoise_snr_validation_2026-11-27.csv.
