# Design: hidden-sibling search (proposal, 2026-11-24 — NOT built)

**Goal.** Find *new* transiting planet candidates in systems that already have a known
planet or TOI, by removing the known signal and searching the remaining multi-sector data
in the regime where single-run pipeline searches are weakest: near-threshold depths that
only become significant after stacking many sectors, and longer periods with few transits
per sector. Built on the frozen, blind-tested v2 vetting (its accuracy is the paper's
credibility section).

## Why this regime can still hold new candidates
SPOC and QLP already re-search each light curve after removing detected signals, so
obvious siblings are found. What remains: (1) signals just below the detection threshold in
any single sector or search run that become significant when many sectors (often across
several years) are stacked; (2) periods ~10–50 d with few transits per sector; (3) systems
where the known planet's residuals or stellar activity hid a weaker signal.

## Targets
Stars with a known planet or TOI and **many sectors of 2-min SPOC data** (priority: near the
ecliptic poles / TESS continuous viewing zones, ≥ 10 sectors), excluding v1, dev and blind
samples' known FPs. Start with ~100 systems; record the selection query and a dated snapshot.

## Search
1. Clean each sector (as v2), mask every known transit (all TOIs/planets on the star, using
   published or our refit ephemerides, mask width 1.5–2 × duration).
2. Stitch all sectors; BLS on a period grid fine enough for the full baseline
   (Δf ≈ duration / (3 · baseline²) spacing in frequency), periods 0.5–50 d, require ≥ 3
   transits with data.
3. Signal threshold: SNR ≥ 7.1 on the stitched data AND the signal present (consistent depth,
   within errors) in at least two independent halves of the data (first vs second half of the
   sectors) — systematics rarely repeat across years.
4. Iterate: mask the strongest new signal, search again (up to 3 new signals per star).

## Injection–recovery (required before trusting any detection)
Standard practice (e.g. Christiansen et al. 2013, ApJS 207, 35 for Kepler). For each target,
inject batman transits into the *real* light curves (before cleaning) on a grid of
depth (200–3000 ppm) × period (0.5–50 d) × random phase, ~200 injections per star, run the
identical search + vetting, and record recovered (period within 1%, epoch within ½ duration).
Products: recovery fraction vs depth and period (a sensitivity map), the false-alarm rate
from the uninjected data, and the fraction of injected planets that v2 vetting wrongly
rejects in the shallow regime (v2 is only blind-tested for depth > 2000 ppm). The paper's
non-detections become quantitative: "in these systems we exclude additional planets with
depth > X ppm and P < Y d at Z% completeness".

## Vetting of a new signal
All v2 checks, plus:
- **Multi-sector consistency:** depth consistent between halves of the data; transit times
  follow a straight-line ephemeris (per-sector times, as in the ephemeris fix).
- **(a) Background eclipsing binaries.** Shallow signals are MORE vulnerable to diluted
  background EBs: a faint neighbour only needs a modest eclipse to mimic a few hundred ppm
  (δ_n = δ_obs / f_n stays ≤ 1 for much fainter stars than for a deep signal). So every
  candidate must pass the **Gaia neighbour-depth check**, and ideally the
  **difference-image centroid** (v3), before it is submitted as a CTOI. Without them a
  shallow candidate is reported as "unvetted for blends".
- **Statistical validation (optional, if time):** TRICERATOPS (Giacalone et al. 2021, AJ
  161, 24), which uses TESS pixel data and nearby stars to estimate false-positive probability.

## (b) The multiplicity argument
In Kepler data, candidates in systems with more than one transiting candidate are far more
often real planets, because false positives rarely line up on the same star by chance
(Lissauer et al. 2012, ApJ 750, 112, "Almost all of Kepler's multiple-planet candidates are
planets"; used for statistical validation by Lissauer et al. 2014, ApJ 784, 44 and Rowe et al.
2014, ApJ 784, 45). This gives a sibling candidate a higher prior probability of being a
planet. Caveats we must state: the boost was derived for Kepler's 4″ pixels; with TESS's 21″
pixels more background stars share the aperture, so the boost is smaller and has to be
argued, not assumed; and it applies only if the known signal is itself a planet on the same
star.

## (c) What counts as "new"
A signal is new only if ALL of these hold on the date of the check (recorded, dated):
1. **Not a TOI or CTOI** on the star at a matching period (within 1%, or at ½×/2×/3×)
   — ExoFOP TOI and CTOI tables (dated snapshots).
2. **Not a SPOC or QLP threshold-crossing event (TCE)** at a matching period — SPOC
   multi-sector TCE lists / DV reports on MAST (Twicken et al. 2018, PASP 130, 064502) and
   QLP candidate products where available (Huang et al. 2020, RNAAS 4, 204).
3. **Not in the literature** — NASA Exoplanet Archive (confirmed planets and candidates),
   plus an ADS/arXiv search on the star's names (TIC, TOI, HD/HIP/Gaia ID).
4. Passes all vetting above.
If 1–3 hold but vetting is incomplete, it is reported as "new signal, unvetted", not as a
candidate.

## Output and next step for a surviving candidate
Per candidate: ephemeris, depth, Rp, SNR, per-sector depths and times, vetting report,
neighbour-depth table, sensitivity context, and novelty evidence (1–3 above with dates).
A surviving candidate is submitted as a CTOI on ExoFOP (with Julien's approval) so the TFOP
community can follow it up.

## Validation summary (what must be true before running on real targets)
- Injection–recovery on ≥ 20 stars shows the expected recovery and a measured false-alarm rate.
- Known multi-planet systems (TOIs .01/.02 on the same star) are recovered when only .01 is masked.
- v2 checks don't reject injected shallow planets beyond the measured rate.
