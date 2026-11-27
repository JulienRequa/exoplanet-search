# Real red noise on dev light curves (2026-11-27; dev only, no blind_v2 star)

- Stars: 251 attempted, 250 measured, 0 errors. Median point scatter 3107 ppm.
- beta_1h: median 1.02, 16–84% 0.96–1.25, 90th pct 1.49
- beta_2h: median 1.00, 16–84% 0.91–1.29, 90th pct 1.61
- beta_3h: median 0.97, 16–84% 0.86–1.39, 90th pct 1.71
- beta_T14: median 0.98, 16–84% 0.85–1.34, 90th pct 1.56
- beta_2h by group: planet hosts median 0.99 (n=124); FP stars 1.00 (n=126); point scatter 2,000–7,000 ppm: median 0.98, 90th 1.30 (n=107) [reported only]

## Simulation calibration (median beta_2h of pure-noise simulated stars vs alpha; tau 1 h)

α 0: 0.90 · α 0.05: 0.93 · α 0.1: 1.02 · α 0.15: 1.14 · α 0.2: 1.30 · α 0.3: 1.64 · α 0.4: 1.99 · α 0.5: 2.32 · α 0.75: 3.05 · α 1: 3.61

Note: white noise alone gives beta_2h ≈ 0.90 with this estimator after detrending (the SNR uses max(1, beta)).

## Pre-registered levels
- typical (median): dev beta_2h 1.004 → alpha 0.091
- pessimistic (90th pct): dev beta_2h 1.614 → alpha 0.292

The earlier threshold simulation assumed alpha 0.5 (beta_2h ≈ 2.3): ~1.4× stronger than the dev 90th percentile.
