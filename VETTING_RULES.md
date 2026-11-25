# Vetting rules v2 — FROZEN 2026-11-25

These are the rules the blind test (data/blind_v2.csv) and the first hunting run use.
They are frozen at git tag `v2-freeze`. Nothing here changes after the blind results are seen;
anything new goes to v3. Developed only on development samples (DEVLOG; docs/PAPER_NOTES.md M1).

## Data
- TESS 2-minute SPOC light curves (lightkurve), first **2 sectors** in archive order; each sector
  used once (duplicate archive listings skipped).
- Stellar Teff, log g, radius from the light-curve header; missing values from the TIC.
- Cleaning per sector: flatten with a 721-point (~24 h) window, transits masked when known;
  upward outliers > 4σ removed.

## Search and fit
- **Period:** BLS within ±2% of the catalog (TOI) period, 2,000 trial periods, durations
  0.05–0.2 d; always the 1× window when it fits in half the data span (½× or 2× only otherwise).
  Search on the first sector.
- **Timing (multi-sector):** a transit time per sector, whole-orbit count, weighted straight-line
  ephemeris; a sector > P/4 from a whole orbit is dropped and noted.
- **Model:** batman, quadratic limb darkening; fit t0, P, Rp/Rs, a/Rs, b, q1, q2 (Kipping 2013,
  q ∈ [0, 1]); bounds Rp/Rs 0.01–1.0, a/Rs 1.5–100, b 0–1.2; curve_fit (trf), maxfev 20,000.
- **Limb-darkening prior:** Claret (2017) TESS table (ATLAS, Z = 0, ξ = 2 km/s), Gaussian σ = 0.1
  on u1, u2; free if Teff/log g missing or off-table.
- **Multi-start:** every fit also starts from b = 0.6; the lower-χ² solution is kept.
- **b upper limit:** profile scan when b is at a bound (Δχ² = 4).

## Checks (thresholds frozen)
| Check | FAIL if | Other outcomes |
|---|---|---|
| Detection SNR | SNR ≤ 7.1 | — |
| Odd/even depths | difference ≥ 3σ (errors × red-noise factor from 300 out-of-transit fake transits, seed 0; fully covered transits only; per-sector differences combined) | N/A if too few transits |
| Secondary eclipse | > 3σ **and** > 10% of transit depth | — |
| Density | transit/catalog ratio < 0.5 or > 2 **and** > 3σ (catalog ±30%) | WEAK if transit density error > 50% |
| Planet size | Rp ≥ 2.5 R_Jup | if Rp/Rs at its 1.0 cap: lower limit → FAIL if ≥ 2.5 R_Jup, else WEAK |
| Grazing | — (never FAIL) | WARN if grazing preferred at > 3σ (profile likelihood, b = s(1 − k)) |
| Centroid | shift > 3σ (null from 200 fake transits, seed 0) **and** source > 1 px from target (proper motion applied) | WEAK if source position error > 1 px |
| Neighbor depth | target rejected (p < 0.0027) **and** a TIC neighbour (≤ 2.5 px, ≤ 7.5 mag fainter, depth ≤ 100%) fits | WEAK if target rejected and no neighbour fits. PSF σ 0.75 px, real SPOC aperture, systematic 10% + 0.0005 px |

## Verdict
- Any FAIL → **LIKELY FALSE POSITIVE**
- else any WARN → **CANDIDATE WITH WARNINGS**
- else any WEAK / N/A → **PLANET CANDIDATE (some checks weak or not run)**
- else → **PLANET CANDIDATE**

## Tested range
P 0.5–10 d, catalog depth > 2000 ppm (the dev and blind samples' cuts). Outside it, results are
untested.

## Development results (method development only — the headline numbers come from blind_v2)
293 dev stars (93-star set + dev sample 2), code e216703: real planets kept **121/124 (97.6%)**,
false positives caught **67/127 (53%)** (v1 baseline 42%), **0 crashes**. 13-planet regression:
no verdict change. Lost: HATS-64 b, K2-19 b (odd/even just over 3σ), Kepler-74 b (grazing branch
of the b–Rp/Rs degeneracy on a faint star).

## Known limitations
- Neighbours closer than ~1 TESS pixel are not separated (centroid / neighbour-depth).
- Fast-rotating, spotted stars: detrending can leave residuals larger than the transit.
- Shallow transits: b and a/Rs poorly constrained (density often WEAK).
- Faint, noisy stars: grazing and central solutions indistinguishable; size check can FAIL a
  real planet (Kepler-74 b).
- Transit-timing variations are not modelled (straight-line ephemeris; K2-19 b).
- Faint planets undetectable in one sector are missed (search uses the first sector).
- About 11% of dev targets (32/293) have no 2-minute data (not testable).

## v3 list (not in v2)
Difference-image centroid (go/no-go first), detrending for spotted stars, multi-sector search
(stricter threshold from the simulated noise false-alarm rate), degeneracy-aware size check,
per-transit depth.
