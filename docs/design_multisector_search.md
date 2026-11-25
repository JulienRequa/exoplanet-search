# Design: multi-sector search with a period hint (proposal, 2026-11-24 — simulated, NOT in pipeline)

## Problem
The v2 search finds the period in ONE sector (the first), then refines with the others.
A planet too faint to stand out in one sector is then found at a wrong period, even inside
the ±2% 1× window, and fails the detection-SNR check. Dev sample 2 lost 4 real planets
this way: Kepler-44 b, TOI-4306, K2-415 b (SNR 2–3) and HAT-P-37 b (crowded).

## Physics
Transit SNR grows like √(number of transits). Stacking 2 sectors gives √2 ≈ 1.4× the SNR of
one; 4 sectors give 2×. A signal at SNR 4–5 in one sector can pass 7.1 with 4 sectors.

## Method
1. Clean every sector (as now) and stitch them.
2. BLS on the stitched light curve within ±2% of the hint (1× only), with a period grid fine
   enough for the full time span: spacing ≈ P · (duration/P) / (3 · span) — e.g. ~1,700 trial
   periods for 2 sectors two years apart (no slower than today's 2,000-point grid per window).
3. Fit the stitched data starting from the BLS answer (tight window), then vet as usual.
4. Optionally use more than 2 sectors (n_sectors > 2) — a separate decision (more downloads).

## What could go wrong (and what the simulation checks)
- More data also means more chances for a spurious "signal" in noise → measure how often
  noise-only stars pass the detection check (must not increase).
- EBs whose true period is 2× catalog must still be caught.
- Runtime per star (the hunting run has ~2,900 stars).

## Simulation (sims/multisector_search_sim.py)
Faint planets (Rp/Rs 0.045–0.065, 3000–5000 ppm noise per 2-min point), 2 and 4 sectors with a
700-day gap; EBs at 2× catalog; noise-only stars with a random hint. Compared: current
single-sector search (+ timing fix) vs stitched search. Results:
sims/multisector_search_sim_results_2026-11-24.md (overnight run).
