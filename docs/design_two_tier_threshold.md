# Design: two-tier detection threshold for noisy stars (v3 multi-sector search) — FROZEN 2026-11-27

Status: approved by Julien 2026-11-27 (three buckets as proposed, plus additions A1 and A2 below).
Implemented in sims/two_tier_threshold_sim.py; its header holds the binding pre-registered rules.
Runs after sims/large_noise_threshold_sim.py finishes (one CPU job at a time).

Julien's decision 2026-11-27: "design a two-tier threshold (each star's measured red-noise factor,
a cutoff fixed in advance, the stricter threshold above it) with its own pre-registered simulation
and validation. Design note first; show me before running." v2 is frozen and not affected.

## Problem
One SNR threshold cannot serve all stars. At the typical measured red-noise level the 1%
false-alarm threshold is SNR ≈ 9.2, but at the dev 90th-percentile level it is ≈ 17
(sims/rednoise_threshold_measured_levels_2026-11-27.md). A single threshold of 17 would throw
away many planets on quiet stars; a single threshold of 9.2 lets noise through on noisy stars.
The red-noise-aware SNR did not solve this (PAPER2_NOTES R1).

## The red-noise factor of a star
beta_2h, exactly as measured on the dev stars (scripts/measure_rednoise_dev.py `measure`): on the
cleaned multi-sector light curve, mask the found signal's primary and secondary ±1 × T14, bin the
remaining points in 2-hour bins (bins < 80% full dropped), beta = robust scatter of the bin means
÷ (single-point scatter / √n · √(M/(M−1))). No floor at 1.
Why this one: it is the statistic whose real distribution we have measured (250 dev stars:
median 1.00, 75th pct 1.16, 90th 1.61, 95th 2.8, 99th 5.5) and whose simulation calibration
exists (alpha → median beta_2h: 0 → 0.90, 0.1 → 1.02, 0.15 → 1.14, 0.2 → 1.30, 0.3 → 1.64,
0.4 → 1.99, 0.5 → 2.32). 2 h is about the median transit duration in the simulated range.

## Principle (physics justification)
In the simulation the false-alarm rate at a fixed threshold rises monotonically with the red-noise
amplitude. So if each tier's threshold is calibrated on simulated noise at the **noisiest edge of
that tier**, every star in the tier has a false-alarm rate at or below the target. This is
conservative by construction; nothing is tuned on real stars.

## Tiers (approved)
| Tier | Star's beta_2h | Calibrated at | Share of dev stars |
|---|---|---|---|
| Normal | ≤ 1.15 | alpha whose median beta_2h = 1.15 (≈ 0.15, interpolated from the existing calibration) | 74% |
| Strict | 1.15 – 2.0 | alpha whose median beta_2h = 2.0 (≈ 0.40) | 18% |
| Too noisy (no claim) | > 2.0 | — no detection claimed; listed as "too noisy for this search" | 8% |

Why 1.15: just above the spread of the typical level (alpha 0.091 → beta 16–84%: 0.97–1.06 in the
calibration) so ordinary stars are not pushed into the strict tier by measurement scatter, and
near the dev 75th percentile. Why 2.0: above it the strict threshold would exceed ~20 and keep
few faint planets; the dev tail above 2.0 (20 of 250 stars, up to beta ≈ 5) cannot be covered by
one more number without simulating far beyond the measured calibration.
Julien approved the third bucket (2026-11-27). The alternative (no upper limit, strict tier
calibrated at the dev 99th percentile, beta ≈ 5.5) was not chosen.

## A1. Conservative tier assignment (Julien, 2026-11-27)
If a star's red-noise factor uncertainty crosses a cutoff, the stricter tier is used.
- Uncertainty: sigma_beta = beta × sqrt(1.36 / M), M = number of 2-h bins (the asymptotic variance
  of the normalised MAD scale estimator is 1.36 sigma²/M). Checked against the existing simulation
  calibration before freezing: it predicts 5.0% for 2-sector stars; the measured 16–84% half-width
  of beta at fixed alpha is 4.2–5.7% across alpha 0–1.0.
- Tier is decided by beta_hi = beta + sigma_beta: ≤ 1.15 normal; ≤ 2.0 strict; above → tier 3.
  A star whose beta cannot be measured (< 10 full bins) goes to tier 3.
- Effect: shifts some stars just under 1.15 (about beta 1.10–1.15) into the strict tier; the
  validation reports the resulting shares.

## A2. Tier-3 stars are searched, never claimed (Julien, 2026-11-27)
Tier-3 stars go through the same search. A signal with SNR ≥ the strict-tier threshold is listed
in a separate low-confidence list, for follow-up only. It is never claimed as a detection and
never counted in any rate or total. The validation reports how often pure noise gets listed, so
readers know how unreliable that list is.

## Simulation 1: calibration (new seeds, rules committed first)
Same star model and search as sims/large_noise_threshold_sim.py (white noise 3–5 ppt, AR(1) τ 1 h).
- Edge alphas: obtained by linear interpolation of the existing calibration table
  (results/rednoise_dev_2026-11-27/sim_alpha_calibration.csv) at beta 1.15 and 2.0; no re-tuning.
- Per tier × 2 and 4 sectors: 10,000 noise stars + 1,000 faint planets at the edge alpha.
- Tier threshold at FAR f (1% and 0.25%): (1 − f) quantile of the noise SNR; with 2 / 4 sectors,
  take the higher, rounded up to 0.1 (same rule as the large run).
- Recorded per star: SNR, found ephemeris, beta_2h measured as above.
- Edge alphas, computed before freezing: 0.153 (beta 1.15) and 0.402 (beta 2.0).
- Seeds: 80,000,000–101,500,999 (calibration), 72,000,000–74,501,999 (validation); disjoint from all earlier runs.

## Simulation 2: validation (new seeds, separate from calibration)
A mixed population that looks like the real dev stars: each simulated star's alpha is drawn by
picking a random dev star's beta_2h and mapping it to alpha through the calibration (beta below
0.90 → alpha 0; above 3.6 → alpha 1.0). 10,000 noise stars + 2,000 faint planets, 2 and 4
sectors. Each star is assigned to a tier by **its own measured beta_hi** (A1; not its true alpha),
as a real star would be, and judged against that tier's threshold. Every planet is also run as its
noise-only twin (same seed, identical noise) to see whether the planet moves its own star's tier.

Pre-registered criteria (for each FAR f in 1%, 0.25%; each sector count):
- T1 (overall): claimed false alarms ÷ all noise stars (tier 3 never claims) ≤ f, judged by the
  binomial 95% interval reaching a value ≤ f.
- T2 (per tier, catches stars misclassified near the cutoff): false-alarm rate in each tier ≤ 2f.
- T3 (worth it): two-tier keeps more faint planets than a single threshold calibrated on the same
  mixed population at the same overall false-alarm rate (computed in the same run).
- Reported, not a gate: share of stars in each tier vs the dev shares (74 / 18 / 8%); share of
  planet stars whose tier changes because of the planet (should be ~0 with the masking).
- PASS = T1, T2, T3. FAIL → not adopted as is; any change is a new committed version on new seeds.

## Runtime and data
Calibration 44,000 stars, validation 28,000 (with the twins): about 4–6 hours with 6 CPU workers, overnight,
no downloads. Real data: none used. After a PASS, the tier shares on real stars are reported on
dev stars only (never blind_v2).

## What it does not fix
- Red noise with a different time scale than 1 h (the AR(1) model) is not tested; the dev
  measurement at 1 h and 3 h bins gave similar medians (1.02, 0.97).
- Stars in the "too noisy" bucket get no multi-sector detection claim; that is a stated loss.
