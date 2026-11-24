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
- [x] Full `analyze_star()` pipeline, blind-tested on WASP-121 b
- [x] Impact-parameter fit + profile-scan upper limits for near edge-on orbits
- [ ] Clean pipeline notebook with multi-sector support
- [ ] Vetting checks: odd/even depths, secondary eclipse, density consistency
- [ ] Search stars without known planets; ML classifier for planet vs. false positive

## Results: WASP-18 b (TESS sectors 2, 3, 29, 30)
| Parameter | Value |
|---|---|
| Period | 0.94145 d |
| Rp/R★ | 0.0991 ± 0.0002 |
| a/R★ | 3.421 ± 0.017 |
| Inclination | 83.2° ± 0.2° |
| Stellar density | 0.85 ± 0.01 g/cm³ |

## Blind test: WASP-121 b (TESS, 1 sector, no hints given)
| Parameter | Value |
|---|---|
| Period | 1.27493 d |
| Rp/R★ | 0.1217 ± 0.0006 |
| a/R★ | 3.82 ± 0.04 |
| Impact parameter | b < 0.22 (95%) → inclination > 86.7° |
| Stellar density | 0.65 g/cm³ |

Both planets agree well with published values.

## Next steps
1. Clean pipeline notebook with multi-sector support and automatic boundary flags
2. Vetting checks for false positives (odd/even depths, secondary eclipse, density)
3. Run on many known planets and compare to published values
4. Search stars without known planets for new candidates
5. Machine-learning classifier (built from scratch in NumPy) for planet vs. false positive

## Notebooks
- **[lab_notebook.ipynb](lab_notebook.ipynb)**: lab notebook: building and testing
  the pipeline on WASP-18 b and WASP-121 b, including every mistake and fix
  [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/main/lab_notebook.ipynb)
- **pipeline.ipynb**: the clean, reusable pipeline (coming next)

## Author
Julien Requa: Phillips Academy Andover '29
