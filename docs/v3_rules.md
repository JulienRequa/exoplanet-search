# v3 rules — adopted decisions (v3 is NOT frozen; v2 rules stay in VETTING_RULES.md, unchanged)

Each rule: what it is, who decided, when, and the evidence. Candidate ideas not yet decided are in
docs/v3_list.md. Nothing here changes v1/v2 results or paper 1. blind_v2 is v3 dev data since 2026-11-28 (D1);
blind_v3 (data/blind_v3.csv) is the v3 blind test, run once with the frozen v3.

## R1. Multi-sector detection threshold: SNR 10.0 (false-alarm rate 0.25%) — SUPERSEDED 2026-11-28 by R2 for the normal tier
- **Rule.** In the v3 multi-sector search, a signal is detected if its detection SNR (v2-style:
  depth ÷ (out-of-transit scatter / √N in-transit points)) is at least **10.0**. This is the
  normal-tier threshold.
- **Decided by** Julien, 2026-11-27 (after the restart): "adopt the 0.25% false-alarm rate. v3
  multi-sector detection threshold = SNR 10.0 (normal tier)."
- **Evidence.** sims/large_noise_threshold_sim.py (rules committed before the run, b7786d5);
  sims/large_noise_threshold_results_2026-11-27.md (5574053). Measured typical red-noise level
  (alpha 0.091, results/rednoise_dev_2026-11-27/), 10,000 noise stars per setup: FAR-0.25%
  threshold 9.96 (2 sectors) / 9.95 (4 sectors); pre-registered rule → higher, rounded up = 10.0.
  Faint planets kept 87% / 97%. "Enough false alarms" test: yes (25 noise stars above in each).
- **Open point (for Julien, after the two-tier calibration).** The frozen two-tier design
  (docs/design_two_tier_threshold.md, 413f7bb) calibrates the normal tier at its noisiest edge
  (beta 1.15, alpha 0.153), not at the typical level (alpha 0.091). That calibration will very likely
  give a normal-tier number somewhat above 10.0. Which number applies to the normal tier, 10.0 or
  the edge-calibrated one, is Julien's decision once both are measured. Until then R1 stands as
  decided. The strict tier and tier 3 follow the two-tier validation.

## R2. Two-tier detection thresholds (false-alarm rate 0.25%), supersedes R1's normal-tier number
- **Rule.** Each star is tiered by its measured red-noise factor (beta_2h + its uncertainty, rules in
  sims/two_tier_threshold_sim.py / docs/design_two_tier_threshold.md):
  - normal tier (beta_hi ≤ 1.15): detected (planet candidate eligible) if SNR ≥ **11.6**;
  - strict tier (1.15 < beta_hi ≤ 2.0): detected if SNR ≥ **25.5**;
  - tier 3 (beta_hi > 2.0): searched; signals with SNR ≥ 25.5 listed separately as low-confidence,
    follow-up only, never claimed.
- **Normal-tier signals with 10.0 ≤ SNR < 11.6** go to "needs follow-up" (design (e)): not claimed as
  candidates, not discarded.
- **Decided by** Julien, 2026-11-28: "Normal tier: 11.6 (edge-calibrated, as the frozen design
  specified) … Strict tier: 25.5 (0.25%)."
- **Evidence.** Rules 413f7bb (pushed before the run); results 32b186a
  (sims/two_tier_threshold_results_2026-11-27.md): validation PASS (T1–T3) at 0.25% and 1%.
