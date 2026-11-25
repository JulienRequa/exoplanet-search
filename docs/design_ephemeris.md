# Design: per-sector ephemeris (proposal, 2026-11-24 — simulated, NOT in pipeline)

## Problem (TOI-2207 b)
The multi-sector refit starts from a 1-sector period and may only move it by 5× the
1-sector error. With few transits per sector, red noise (which makes curve_fit's error too
small) and a long gap, the predicted sector-2 transit times drift by hours and the refit
cannot recover: the later sector's transits are missed (TOI-2207 b: 4.7–9.5 h over 136 cycles).

## Physics / method
A transit time is measured far more precisely than a period from one sector. So:
1. Fit sector 1 → t₁, P₁ (as now).
2. Find sector 2's transits on their own: BLS within ±2% of P₁ on sector 2 only, then fit → t₂.
3. Count whole orbits: N = round((t₂ − t₁)/P₁); period P = (t₂ − t₁)/N (straight-line
   ephemeris; with ≥ 3 sectors, a least-squares line through all transit times).
4. Final joint fit starts from (t₁, P) with a tight window.
Cycle-count safety: N is wrong only if P₁'s drift exceeds P/2 (days); simulated drifts are
hours (median 3.5 h in the hard case), a margin of ~20–60×. A check: the residual
|t₂ − (t₁ + N·P₁)| must be < P/4, otherwise status "ephemeris ambiguous" (no vetting verdict).

## Simulation (sims/ephemeris_sim.py; results in sims/ephemeris_sim_results_2026-11-24.md)
- easy (P 2–9 d, white noise, 700 d gap): lost — current 0%, per-sector 0%.
- hard (TOI-2207-like: P 6–9 d, Rp/Rs 0.067, red = white noise, 1100 d gap):
  lost — **current 40%, per-sector 0%**.

## Effect on vetting
When the ephemeris is lost, odd/even compares real transits with empty windows (false FAIL)
and the later sector's transits are not masked when flattening. Fixing the ephemeris
removes that failure mode; it doesn't touch any threshold.
