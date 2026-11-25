# Period-window simulation (sims/period_window_sim.py), 2026-11-24

One TESS sector; spotted star (rotation 0.5–5 d, 0.5–3% amplitude, evolving) + 1000 ppm white
noise; cleaned with pipeline.clean(); BLS power in ±2% windows around 1×, ½×, 2× the catalog period.
R* = 99.7th percentile of power(best alternative)/power(1×) over 100 calibration planets whose
true period IS the catalog period. Pre-registered criterion: a genuine alternative gives ≈ 2,
so the rule only makes sense if R* < 2.

```
CALIBRATION (100 planets, spotted): ratio alt/1× median 0.72, max 4.28, 99.7th pct R* = 3.73  (genuine ½×/2× alternative expected ≈ 2)
FRESH PLANETS (60, spotted): keep 1× — current 75%, proposed 100%
EB eb60 (30, true P = 2× catalog): caught — current 70%, proposed 70%; proposed picks {'1×': 30, '½×': 0, '2×': 0}
EB eb90 (30, true P = 2× catalog): caught — current 47%, proposed 50%; proposed picks {'1×': 30, '½×': 0, '2×': 0}
```
