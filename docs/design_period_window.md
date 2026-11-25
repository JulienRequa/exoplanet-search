# Design: period window with a catalog hint (proposal, 2026-11-24 — simulated, NOT in pipeline)

## Problem (HIP 67522 b)
With a period hint, search_near takes the window (1×, ½×, 2×) with the highest BLS power.
On a spotted star, rotational variability that survives detrending makes box-like dips that
can beat the transit in a narrow window → wrong period → untestable.

## Physics
BLS power here is a log-likelihood ∝ depth² × number of in-transit points. A GENUINE
alternative beats the catalog window by ≈ 2× (e.g. true period = ½ catalog: the ½× fold
stacks all transits, the 1× fold only half of them at full depth). Spot residuals have no
reason to prefer any window, and can be large.

## Pre-registered rule and result
Keep 1× unless power(alt)/power(1×) > R*, R* = 99.7th percentile for simulated true-1×
planets on spotted stars; usable only if R* < 2. Simulation: **R* = 3.73 > 2** → spot
artifacts can exceed any genuine alternative, so **always use the 1× window when a hint is
given** (record the ½×/2× powers as diagnostics only).
- Fresh spotted planets kept at 1×: current 75%, proposed 100%.
- EBs with true P = 2× catalog: caught 70% (secondary 60% depth) / 50% (90%) under both
  rules — limited by odd/even power on spotted stars; the 1× fold is what lets odd/even
  compare primaries with secondaries.
Not tuned on HIP 67522 b (its outcome was not used to set anything).
