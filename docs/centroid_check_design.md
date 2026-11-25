# Centroid check — design proposal (2026-11-24, NOT built yet)

**Question:** when the star dims, does the light's "center of mass" on the detector move?
If the dip happens on a neighboring star blended into our pixels, the center of light
shifts *away from* the neighbor during transit (the neighbor's share of the light drops).

## 1. Data: use the centroids already in the light-curve files (recommended for v2)
SPOC 2-min light-curve files already contain a flux-weighted centroid for every cadence
(`MOM_CENTR1/2` → lightkurve `lc.centroid_col`, `lc.centroid_row`, in pixels), computed
over the same aperture as the flux. So **no extra downloads**.
Target pixel files (TPFs) would allow difference images (which pixel the dip comes from),
but they are ~50 MB per sector: ~20 GB for the 200-star test. Upgrade for later.

## 2. Measurement (per sector, because pixels/orientation change between sectors)
1. Take centroid_col and centroid_row; remove the same bad points as `clean()`.
2. Detrend each with a running median (same ~24 h window as the flux), **masking transits**
   (spacecraft pointing drifts move centroids slowly).
3. In-transit points: the deep core of each transit (model < 1 − 0.5·depth).
   Out-of-transit points: a local baseline from 1× to 3× the transit duration on both sides.
4. Shift: Δx = mean(in) − mean(out), Δy likewise.

## 3. Significance (the important part)
Centroid noise is *not* white: pointing jitter makes neighbouring points correlated, so the
textbook σ = √(var_in/N_in + var_out/N_out) underestimates the noise.
**Empirical null:** put the same in/out windows at ~200 random fake transit times
(avoiding real transits and phase 0.5), measure Δx, Δy there → their spread σx, σy is the
real noise, including red noise and jitter.
- Per sector: D² = (Δx/σx)² + (Δy/σy)², which has a χ² distribution with 2 degrees of freedom.
- Combine sectors by adding D² (2 dof per sector). "Significant" = probability of chance
  < 0.27% (the 3σ level): D² > 11.8 for one sector, > 16.8 for two.

## 4. Where is the source? (avoids false alarms on deep, real transits)
Even an on-target transit moves the centroid a little if other stars add light to the
aperture. So significance alone isn't enough. Following the Kepler method
(Bryson et al. 2013), the dimming source sits at about
  **source offset ≈ −Δ · (1 − δ) / δ**  (δ = transit depth; offset measured from the out-of-transit centroid).
Compare that with where the target star itself sits (its TIC position converted to pixels).
The offset error comes from σx, σy scaled the same way.

## 5. Proposed rule (YOUR thresholds to choose; numbers below are placeholders)
- **FAIL** (signal on a neighbour): shift significant (3σ, as above) AND implied source is
  > 1 TESS pixel (21″) from the target at > 3σ.
- **WEAK**: offset uncertainty > 1 pixel (too noisy to say).
- **PASS**: otherwise. **N/A**: no centroid columns / too few in-transit points.
- Optional: significant shift but < 1 pixel → WARN (possible blend very close to target).

## 6. Known limits
- Neighbours closer than ~1 pixel (21″) can't be separated: low power for tight blends.
- Very bright/saturated stars (bleed columns) have unreliable centroids → should be WEAK.
- If the contaminating star is outside the aperture, little of its dimming reaches the flux at all.

## 7. Validation BEFORE the v2 freeze (no peeking at test data)
1. 13 known planets: all should PASS (or WEAK), none FAIL.
2. Injection: shift the real centroid series by a fake transit-shaped offset → check it's
   recovered at the expected significance.
3. A small **development set** of ~10 FP TOIs whose ExoFOP Comments mention a centroid
   offset, NOT in v1, and **excluded from the v2 random sample**, so the blind test stays blind.
