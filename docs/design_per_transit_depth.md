# Design: per-transit depth consistency (proposal, 2026-11-24 — NOT built)

## Why it works (physics)
A planet blocks the same fraction of the same star every orbit, so every transit has the
same depth (to within noise; starspot crossings add small bumps). Depths that change are a
sign that something other than a planet is going on:
1. **Signal on a neighbour, different apertures per sector.** SPOC corrects each sector's
   flux for crowding assuming the variability is on the target (CROWDSAP). If the dip is
   really on a neighbour, that correction is wrong by a different amount in each sector
   (different aperture, different neighbour fraction) → **depth changes sector to sector**.
2. **False alarms / systematics.** If BLS stacked a few real-looking events with
   momentum-dump or scattered-light glitches, individual "transits" disagree wildly and
   some are missing.
3. **Variable dilution** by a nearby variable star.
(Unequal alternating depths are already covered by odd/even.)

## Method
1. Shape fixed from the global fit; for each transit fit only a depth scale factor
   s_i plus a local straight-line baseline → depth d_i = s_i·δ ± σ_i.
2. Red noise: σ_i inflated by a factor measured the same way as the centroid null — fit
   the same model at ~200 fake transit times; the scatter of those fake depths is the real
   per-transit noise.
3. Statistics: χ² = Σ ((d_i − d̄)/σ_i)² with N − 1 dof (all transits), and a separate
   per-sector comparison (sector means). Proposed: WARN at p < 0.27% for the all-transit
   test; FAIL for a significant **sector-to-sector** difference (physics reason 1).
   Status N/A if < 3 transits.

## Confirming it doesn't flag real planets
1. **Injection–recovery:** inject constant-depth transits into real light curves of quiet
   stars (several depths, periods, magnitudes). Flag rate must be ≈ 0.27% — if it's
   higher, the red-noise inflation is wrong and must be fixed before use.
2. **Spotted stars:** the known-planet set includes active hosts (WASP-19, WASP-43,
   HD 189733). They must not be flagged; if they are, spot noise needs to be modelled
   (e.g. use in-transit residual scatter) — not a threshold change.
3. **Power curve:** inject depth changes of 5–50% between sectors → detection rate vs size.

## Cost
Light: uses the cleaned light curve we already have; ~N small fits per star.
