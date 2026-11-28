# v3 rules — adopted decisions (v3 is NOT frozen; v2 rules stay in VETTING_RULES.md, unchanged)

Each rule: what it is, who decided, when, and the evidence. Candidate ideas not yet decided are in
docs/v3_list.md. Nothing here is ever applied to v1/v2 results or to blind_v2.

## R1. Multi-sector detection threshold: SNR 10.0 (false-alarm rate 0.25%)
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
