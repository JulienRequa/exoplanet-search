# Exoplanet Search with TESS

Learning Python by building a pipeline that searches NASA TESS data for
transiting exoplanets, using physics (orbital mechanics, transit geometry,
limb darkening) to measure planets and tell real ones from false positives.

## Status (in progress)
- [x] Download, clean, and fold TESS light curves with `lightkurve`
- [x] Period search with Box Least Squares (recovers WASP-18 b with no input)
- [x] Mask transits during detrending; trim spacecraft "ramp" after data gaps
- [x] Limb-darkened transit model fit (batman + SciPy) with fitted period/timing
- [x] Stellar density from Kepler's 3rd law, with error propagation
- [x] Multi-sector fit (4 TESS sectors)
- [ ] Full `analyze_star()` pipeline tested on a second planet
- [ ] Vetting checks: odd/even depths, secondary eclipse, density consistency
- [ ] Search unflagged stars; ML classifier for planet vs. false positive

## Results: WASP-18 b (TESS sectors 2, 3, 29, 30)
| Parameter | Value |
|---|---|
| Period | 0.94145 d |
| Rp/R★ | 0.0991 ± 0.0002 |
| a/R★ | 3.421 ± 0.017 |
| Inclination | 83.2° ± 0.2° |
| Stellar density | 0.85 ± 0.01 g/cm³ |

## Notebook
- **[View the notebook on GitHub](exoplanet_search.ipynb)**: code, plots, and notes
- **Run it yourself:** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/main/exoplanet_search.ipynb)

## Author
Julien Requa: Phillips Academy Andover '29
