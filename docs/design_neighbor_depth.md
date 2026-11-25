# Design: neighbour-depth check (proposal, 2026-11-24 — NOT built)

## Why it works (physics)
Suppose a neighbour contributes a fraction f_n of the light in our aperture. If *it* were
the star eclipsing, its own eclipse depth would have to be

    δ_n = δ_obs / f_n          (δ_obs = depth we measure)

A star cannot lose more than all its light, so **any neighbour with f_n < δ_obs cannot be
the source**. Example: a 2% dip can't come from a neighbour supplying 1% of the light.
That alone rules out most faint neighbours with no fitting at all.

For neighbours that *could* make the depth, the moment centroid we already measure makes a
**specific prediction**: if source n eclipses, the centroid moves by

    Δ_pred(n) = −δ_obs · (x_n − c_out) / (1 − δ_obs)     (vector, pixels)

(same algebra as the current check, run forwards). Each hypothesis — "target", "neighbour
1", "neighbour 2", … — predicts a direction *and* size of shift. We compare each with the
measured Δ ± σ (σ from the fake-transit null we already compute):

    χ²(n) = Σ_sectors [(Δx − Δx_pred)/σx]² + [(Δy − Δy_pred)/σy]²

This does **not** need the source to be > 1 px away: a 0.5 px neighbour in a known
direction predicts a known shift, and we test that exact prediction. That is why it
targets the < 1 px problem, using only data we already have.

## Method
1. Neighbours: TIC v8.2 (Gaia-based positions + TESS magnitudes) within 2.5 px (52″),
   proper motion propagated to each sector's epoch; position in pixels via the file WCS.
2. f_n: flux of each star falling in the SPOC aperture, from its TESS magnitude and a PRF
   model (Gaussian first), divided by the total. Cross-check: our predicted target fraction
   should match SPOC's CROWDSAP; if not, flag WEAK.
3. Feasible sources: target + neighbours with δ_obs / f_n ≤ 1.
4. χ² for each feasible hypothesis.
   - Target hypothesis rejected (p < 0.27%) AND some neighbour hypothesis acceptable
     → **FAIL** (signal on neighbour n — report which TIC).
   - Target rejected and no neighbour fits → WEAK (uncatalogued source / systematics).
   - Target acceptable → PASS (record which neighbours are also still possible).
   This **replaces** the current "> 1 px" rule in the centroid check (they use the same Δ).

## Confirming it doesn't flag real planets
1. **Simulated scenes:** a target + catalogue-like neighbours, Gaussian PRF, SPOC-like
   aperture, real noise (taken from quiet stars' centroid series). Eclipse the target →
   FAIL rate must be ≈ 0.27%. Eclipse a neighbour at 0.3–2 px, various Δmag → power curve.
2. **Systematic floor:** the WCS-vs-centroid offset (≈ 0.3–0.4 px seen on WASP-18) enters
   Δ_pred through c_out vs x_n. We measure this floor on the 13 known planets and the
   simulations and add it as σ_sys — never tuned on FPs.
3. **13 known planets:** target hypothesis must be accepted for all.

## Cost
Light: one TIC cone query per star; no new downloads; no new packages.
