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
`exoplanet_search.ipynb`: open with the "Open in Colab" badge to run it.

## Author
Julien Requa: Phillips Academy Andover '29
