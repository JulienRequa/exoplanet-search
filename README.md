# Exoplanet Search with TESS

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24571362.svg)](https://doi.org/10.5281/zenodo.24571362)

An open-source Python pipeline that checks NASA TESS planet candidates and flags the likely
false positives. For each star it downloads the TESS light curve, finds the transit near the
catalog period, fits a physical transit model, and runs eight vetting checks based on physics
(orbital mechanics, transit shape, stellar density, and where on the sky the signal comes from).
It returns a verdict: **planet candidate** or **likely false positive**.

## What it does
1. **Download** TESS 2-minute light curves (2 sectors) with `lightkurve`.
2. **Clean** each sector: remove slow trends with transits masked, drop outliers.
3. **Search** for the transit within ±2% of the catalog period (Box Least Squares).
4. **Time** the transit in each sector and fit one straight-line ephemeris, so transits years
   apart stay aligned.
5. **Fit** a limb-darkened transit model (batman), with theoretical limb-darkening priors and
   two starting points to avoid getting stuck at an edge-on solution.
6. **Vet** with eight checks: detection strength, odd/even depths, secondary eclipse, stellar
   density, planet size, grazing geometry, centroid shift, and neighbour-star depth.

All checks, thresholds and fit settings are listed in **[VETTING_RULES.md](VETTING_RULES.md)**.

## How it was validated
The rules were built on development samples and tested on a blind sample, like a
pre-registered experiment:

1. **Blind sample locked first.** 100 confirmed/known planets + 100 TFOPWG false positives
   were drawn from a dated ExoFOP snapshot (`data/blind_v2.csv`, `random_state=2`) and
   committed in [0217ed9](https://github.com/JulienRequa/exoplanet-search/commit/0217ed9)
   before any v2 change.
2. **Development only on other stars.** Every change was designed, simulated, and checked on
   293 development stars; thresholds were never tuned to rescue a particular star.
3. **Freeze.** Rules frozen at tag
   [`v2-freeze`](https://github.com/JulienRequa/exoplanet-search/tree/v2-freeze). The blind
   file was unchanged at that point: md5 `66805c4c73b34c6cef991a2e3fc7d2f1`.
4. **Blind run once**, on the frozen code, reported as-is. Changes after that go to v3.

The full record of every change, test and decision is in [DEVLOG.md](DEVLOG.md).

**Archived release.** The hunt results (frozen v2, 2,896 TOIs) are archived on Zenodo as
release `hunt-v2-results-2026-11-26`:
[doi:10.5281/zenodo.24571362](https://doi.org/10.5281/zenodo.24571362)
(zip md5 `f7f5351c6d39d1b42d8e2abbfa58fcb4`). Please cite this DOI when using these results.

## Results
| | Real planets kept | False positives caught | Crashes |
|---|---|---|---|
| v1 (blind, 20 FPs) | — | 5 / 12 testable (42%) | — |
| v2 development (293 stars) | 121 / 124 (97.6%) | 67 / 127 (53%) | 0 |
| **v2 blind (200 stars)** | **91 / 99 (92%)** | **43 / 85 (51%)** | **0** |

"Testable" = the star has TESS 2-minute data and the pipeline found the catalog period.
The blind result is reported exactly as it came out of the frozen code. Planets kept (92%) is
below the development figure: of the 8 real planets lost, 5 failed odd/even (3.1–5.2σ), 2 failed
density (3.0–3.2σ) and 1 failed planet size (fit at the grazing limit). These are being studied
for v3; v2 is not changed. Per-star results: `results/blind_v2_results_2026-11-25.csv`.

## Tested range and known limitations
Tested on orbital periods **0.5–10 days** and transit depths **> 2000 ppm**. Outside that
range the results are untested. Known limitations (details in
[VETTING_RULES.md](VETTING_RULES.md#known-limitations)):
- Neighbour stars closer than ~1 TESS pixel (21″) cannot be separated.
- Very active, fast-rotating stars can leave residuals larger than the transit.
- Shallow transits give a weak density test; on faint, noisy stars a grazing and a central
  transit can look the same.
- Transit-timing variations are not modelled.
- Planets too faint to see in one sector are missed.

## Quick start
```bash
conda env create -f environment.yml
conda activate exoplanet
```
```python
import pipeline as pl
results = pl.analyze_star("WASP-18", n_sectors=2, vet=True, period_hint=0.94145)
```
One call: download → clean → search → fit → vetting report with a verdict.

Other tools:
- `python regression.py run OUT.csv`: re-measure 13 known planets (checks that code changes
  don't change results)
- `python dev_set.py summary OUT.csv`: planets kept and false positives caught for a run
- `pytest`: unit tests

## What's next
- **Hunt: done.** Frozen v2 ran once on all 2,896 open TESS planet candidates in the tested
  range (target list pushed before the run). Results are reported as-is in release
  [`hunt-v2-results-2026-11-26`](https://github.com/JulienRequa/exoplanet-search/releases/tag/hunt-v2-results-2026-11-26),
  archived at [doi:10.5281/zenodo.24571362](https://doi.org/10.5281/zenodo.24571362).
- **Checking the flagged candidates:** diagnostic checks (variability, nearby eclipsing
  binaries, difference images) on the targets v2 flags as likely false positives. These are
  candidate flags only, not confirmed, and nothing is reported to ExoFOP before they are checked.
- **Paper 1:** the v2 tool, its blind-test accuracy, the hunt, and the flagged candidates.
- **v3:** difference-image centroids, detrending for spotted stars, multi-sector search, a
  degeneracy-aware size check, and a secondary-eclipse check robust to stellar variability;
  tested on a fresh blind sample. v3 never changes the v2 results above.
- **Hidden siblings (paper 2):** search the same stars for additional, smaller planets with v3,
  after injection-recovery tests.

## History
v1 was built and tested in Google Colab. The original notebooks (as of the v1 test):
- **pipeline.ipynb**: the v1 pipeline
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/3ed189d/pipeline.ipynb)
- **lab_notebook.ipynb**: lab notebook: building and testing the pipeline step by step,
  including every mistake and fix
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/3ed189d/lab_notebook.ipynb)

## Author
Julien Requa: Phillips Academy Andover '29
