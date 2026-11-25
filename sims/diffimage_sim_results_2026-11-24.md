# Difference-image centroid simulation (sims/diffimage_sim.py), 2026-11-24

11×11 px stamp, Gaussian PSF σ = 0.75 px; target + one neighbour; the neighbour eclipses (or the
target, for the NULL rows). Detection = offset from the target position > 3σ.
Method A: per-transit difference images, error from their scatter (Hotelling T²).
Method B: one stacked difference image, error from the PSF fit with per-pixel noise measured
from out-of-transit frames. Idealised: no jitter, no PRF variation, no WCS/Gaia error.

```
case                                   detect >3σ  measured sep (median)
--- Tmag 11, depth 0.5%, 6 transits, PSF exact
NULL: eclipse on target (20% nbr)              0%             0.014 px
neighbour 20% flux at 0.3 px                 100%             0.299 px
neighbour 20% flux at 0.5 px                 100%             0.500 px
neighbour 20% flux at 0.8 px                 100%             0.801 px
neighbour 20% flux at 1.2 px                 100%             1.199 px
neighbour 5% flux at 0.3 px                  100%             0.301 px
neighbour 5% flux at 0.5 px                  100%             0.499 px
neighbour 5% flux at 0.8 px                  100%             0.802 px
neighbour 5% flux at 1.2 px                  100%             1.199 px
--- same, PSF model 15% too wide
NULL: eclipse on target (20% nbr)              0%             0.014 px
neighbour 20% flux at 0.3 px                 100%             0.301 px
neighbour 20% flux at 0.5 px                 100%             0.498 px
neighbour 20% flux at 0.8 px                 100%             0.806 px
neighbour 20% flux at 1.2 px                 100%             1.199 px
neighbour 5% flux at 0.3 px                  100%             0.300 px
neighbour 5% flux at 0.5 px                  100%             0.500 px
neighbour 5% flux at 0.8 px                  100%             0.799 px
neighbour 5% flux at 1.2 px                  100%             1.198 px
--- HARD: Tmag 13, depth 0.2%, 3 transits, PSF 15% wide
NULL: eclipse on target (20% nbr)              0%             0.113 px
neighbour 20% flux at 0.3 px                   0%             0.290 px
neighbour 20% flux at 0.5 px                   0%             0.509 px
neighbour 20% flux at 0.8 px                   0%             0.820 px
neighbour 20% flux at 1.2 px                   2%             1.174 px
neighbour 5% flux at 0.3 px                    0%             0.312 px
neighbour 5% flux at 0.5 px                    0%             0.501 px
neighbour 5% flux at 0.8 px                    0%             0.789 px
neighbour 5% flux at 1.2 px                    2%             1.159 px

case                                   detect >3σ  measured sep (median)
--- HARD: Tmag 13, depth 0.2%, 3 transits, PSF 15% wide
NULL: eclipse on target (20% nbr)              0%             0.105 px
neighbour 20% flux at 0.3 px                   0%             0.326 px
neighbour 20% flux at 0.5 px                   0%             0.488 px
neighbour 20% flux at 0.8 px                   2%             0.787 px
neighbour 20% flux at 1.2 px                   4%             1.210 px
neighbour 5% flux at 0.3 px                    0%             0.335 px
neighbour 5% flux at 0.5 px                    0%             0.501 px
neighbour 5% flux at 0.8 px                    0%             0.778 px
neighbour 5% flux at 1.2 px                    2%             1.176 px
--- HARD, method B (stacked image, empirical pixel noise)
NULL: eclipse on target (20% nbr)              0%             0.081 px
neighbour 20% flux at 0.3 px                  34%             0.267 px
neighbour 20% flux at 0.5 px                 100%             0.441 px
neighbour 20% flux at 0.8 px                 100%             0.706 px
neighbour 20% flux at 1.2 px                 100%             1.115 px
neighbour 5% flux at 0.3 px                   34%             0.269 px
neighbour 5% flux at 0.5 px                   96%             0.430 px
neighbour 5% flux at 0.8 px                  100%             0.703 px
neighbour 5% flux at 1.2 px                  100%             1.078 px
```

First attempt (χ² with per-transit scatter, 6 transits) had a 4–8% NULL false-alarm rate:
small-sample error estimate treated as exact. Fixed with Hotelling T² (A) / method B.
