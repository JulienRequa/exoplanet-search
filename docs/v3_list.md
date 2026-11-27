# v3 list (living document; v2 is frozen — see VETTING_RULES.md)

Each item needs a physics or simulation justification and is developed on dev data only
(blind_v2 becomes v3 dev data if its losses are investigated; v3 gets a fresh blind sample).
v3 never changes paper 1's numbers.

From VETTING_RULES.md v2 (frozen list):
- Difference-image centroid. Go/no-go done 2026-11-26 (dev data): PARTIAL GO, positions accurate
  (100% of 79 known planets within 0.5 px), formal errors too small → use the empirical floor
  0.19 px (≈ 4″). PAPER_NOTES M9.
- Detrending for spotted / fast-rotating stars.
- Multi-sector search (stricter threshold from the simulated noise false-alarm rate).
- Degeneracy-aware size check (b–Rp/Rs; Kepler-74 b, TOI-1408 b).
- Per-transit depth.

Added 2026-11-26 (Julien):
- Secondary-eclipse check robust to phase curves and stellar variability: local baseline +
  red-noise null (report cards: 9 of 13 v2 secondary flags not robust; PAPER_NOTES L5).
- Variability near the orbital period (P, P/2) flagged explicitly.
- Nearby eclipsing binaries: catalogue cross-match (VSX, Gaia DR3 EBs, TESS EB catalogue) plus
  the period/epoch timing test.
- Density check: (1) report explicitly when it can't run (missing catalogue values) with a status
  of its own (v2 leaves the status blank and writes "no catalog data" in the detail; 388 of 2,397
  tested hunt targets, 16%), and (2) fall back to Gaia DR3 (GSP-Phot / FLAME radius and mass) or
  TIC radius/mass where available. Also: use a b-profile rather than covariance errors when b is
  unconstrained (the 4 strong candidates had density errors of 100–1,000%).

Added 2026-11-27 (Julien):
- Detection SNR red-noise aware: Pont et al. (2006) time-averaging (noise measured in
  transit-length bins), validated on the threshold simulation; then measure the real red-noise
  level on dev-set light curves (never blind_v2) and rerun the threshold simulation at that level.
  Julien chooses the false-alarm rate from those results. (PAPER_NOTES L6.)
  Outcome 2026-11-27: the red-noise-aware SNR FAILED validation, not adopted (PAPER2_NOTES R1).
  Julien: provisional multi-sector threshold at FAR 1% (SNR ≈ 9.1–9.2, typical level); target
  0.25% after a 10,000-star-per-setup run (sims/large_noise_threshold_sim.py); noisiest stars by a
  two-tier threshold on each star's measured beta (docs/design_two_tier_threshold.md, design only).
