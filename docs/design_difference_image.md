# Design: difference-image centroid (proposal, 2026-11-24 — NOT built)

## Why it works (physics)
Every pixel's brightness is a sum of all stars' light spread by the telescope's point
spread function (PSF; for TESS the pixel response function, PRF, ~1–2 px wide).
Average the target pixel file (TPF) images **in transit** and **just outside transit** and
subtract: every star that stays constant cancels exactly, and what is left is the PRF of
**the one source that dimmed**, scaled by how much light it lost. Fitting a PRF to this
difference image gives the dimming source's position directly.

Why it beats the moment centroid we have now:
- No aperture truncation: it uses the whole 11×11+ pixel stamp, so a source at the aperture
  edge is not pulled toward the centre (the 1.22 px neighbour measured at 0.53 px).
- A PRF fit localises a bright signal to a small fraction of a pixel
  (σ_pos ≈ PRF width / SNR of the difference image), so neighbours at 0.3–0.6 px become
  resolvable when the transit is well detected. This is the Kepler/TESS-SPOC DV method
  (Bryson et al. 2013; Twicken et al. 2018).

## Method
1. Download SPOC 2-min TPFs for the same sectors (lightkurve `search_targetpixelfile`).
2. Per sector: in-transit cadences = deep core of the fitted model; out-of-transit = 1–3
   durations either side (as in the current check). Difference image D = mean(out) − mean(in).
3. Fit a PRF model to D → (x_src, y_src). Fit the same model to mean(out) for reference.
4. Target position: Gaia DR3 RA/Dec propagated to the sector epoch with proper motion,
   converted to pixels with the TPF's WCS.
5. Uncertainty: fit each transit's own difference image; the scatter of per-transit
   positions / √N is the statistical error (includes red noise). Add a systematic floor
   σ_sys (WCS + PRF model error), measured on injections and known planets (below).
6. Offset significance: |src − target| / √(σ_stat² + σ_sys²). Proposed: FAIL if > 3σ.
   No fixed 1-px distance cut — significance alone, because σ is now honest.

PRF model choice (your decision): (a) a 2D Gaussian fitted per sector (no new packages,
small known bias), or (b) the official TESS PRF models via a package such as `tess-prf`
(would need installing — I'd ask first).

## Confirming it doesn't flag real planets
1. **Injection on-target:** take TPFs of quiet stars, multiply the target's PRF-weighted
   pixels by a transit model (signal truly on target), run the check. Over many injections
   the FAIL rate must be ≈ 0.27% (the 3σ false-alarm rate) — this measures σ_sys honestly.
2. **Injection on a synthetic neighbour** at 0.3, 0.5, 1, 2 px and several brightness
   ratios: measure recovered position vs truth → power curve and any bias.
3. **13 known planets:** offsets must be consistent with 0. Any that aren't reveal a
   systematic we haven't modelled (e.g. saturated pi Men → mark WEAK, don't fit).

## Cost
TPFs are ~30–50 MB per sector: ~7 GB for the dev set (93 stars × 2 sectors), ~16 GB for
the blind set. Fitting is fast. This is the biggest piece of work of the three.
