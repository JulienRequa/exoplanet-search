# Exoplanet Search with TESS

Learning Python by searching NASA TESS data for transiting exoplanets,
and using physics (orbital mechanics, transit geometry) to tell real
planets from false positives.

## Status (in progress)
- [x] Download and plot TESS light curves with `lightkurve`
- [x] Recover the known hot Jupiter WASP-18 b and fold its transit
- [ ] Estimate planet size from transit depth: δ ≈ (Rp/Rs)²
- [ ] Orbital distance from Kepler's third law
- [ ] Own period search (box least squares) and transit model
- [ ] Vet unconfirmed TESS candidates

## First result: WASP-18 b (TESS Sector 2)
- Period ≈ 0.94 days
- Transit depth ≈ 1.05% → planet radius ≈ 0.10 × star radius
- Duration ≈ 2 hours; rounded bottom from limb darkening

## Notebook
- **[View the notebook on GitHub](exoplanet_search.ipynb)**: code, plots, and notes
- **Run it yourself:** [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/JulienRequa/exoplanet-search/blob/main/exoplanet_search.ipynb)

## Author
Julien Requa: Phillips Academy Andover '29
