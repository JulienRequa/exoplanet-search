# Exoplanet Search with TESS

A Python pipeline that searches NASA TESS light curves for transiting exoplanets,
measures them with a limb-darkened transit model, and vets them for false positives
using physics (orbital mechanics, transit geometry, stellar density).

## Quick start
```python
results = analyze_star("WASP-18", n_sectors=4, vet=True)
```
One call: download → clean → period search → transit fit → vetting report with a verdict.

## Status (in progress)
- [x] Download, clean, and fold multi-sector TESS light curves (`lightkurve`)
- [x] Trim spacecraft ramps around data gaps; mask transits during detrending
- [x] Box Least Squares period search
- [x] Limb-darkened transit fit (batman + SciPy) with fitted period, timing, and impact parameter
- [x] Automatic boundary flags + profile-scan upper limits for near edge-on orbits
- [x] Vetting checks: detection SNR, odd/even depths, secondary eclipse, stellar density, planet size
- [x] One-screen vetting report with verdict
- [ ] Batch run on many known planets → summary table
- [ ] Search stars without known planets for new candidates
- [ ] ML classifier (from scratch in NumPy) for planet vs. false positive

## Results
| | WASP-18 b (4 sectors) | WASP-121 b (2 sectors, blind test) |
|---|---|---|
| Period | 0.94145 d | 1.27492 d |
| Rp/R★ | 0.0992 ± 0.0002 | 0.1208 ± 0.0004 |
| Planet radius | 1.33 R_Jup | 1.83 R_Jup |
| a/R★ | 3.418 ± 0.016 | 3.828 ± 0.031 |
| Inclination | 83.17° ± 0.24° | > 87.1° (95%, profile scan) |
| Stellar density (transit / catalog) | 0.85 / 0.69 g/cm³ | 0.65 / 0.58 g/cm³ |
| Secondary eclipse (planet's own glow) | 487 ± 9 ppm | 379 ± 25 ppm |
| Verdict | Planet candidate (5/5) | Planet candidate (5/5) |

Both agree well with published values. Both secondary eclipses are detections of the
planets' thermal emission; depths are likely overestimated (ellipsoidal variation and
detrending; see notes).

## Notebooks
- **[pipeline.ipynb](pipeline.ipynb)**: the clean, reusable pipeline
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/main/pipeline.ipynb)
- **[lab_notebook.ipynb](lab_notebook.ipynb)**: lab notebook: building and testing the
  pipeline step by step, including every mistake and fix
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/main/lab_notebook.ipynb)

## Author
Julien Requa: Phillips Academy Andover '29
