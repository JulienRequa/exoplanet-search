# Results: sims/two_tier_threshold_sim.py (2026-11-27)

Edge alphas: normal 0.153, strict 0.402

## Calibration
- normal, 2 sectors: median beta2h 1.13; sigma_beta formula 5.1% vs empirical 5.3%; FAR 1.00% quantile 10.76 (planets kept 86%); FAR 0.25% quantile 11.53 (planets kept 82%)
- normal, 4 sectors: median beta2h 1.14; sigma_beta formula 3.6% vs empirical 3.8%; FAR 1.00% quantile 10.76 (planets kept 95%); FAR 0.25% quantile 11.51 (planets kept 94%)
- strict, 2 sectors: median beta2h 1.95; sigma_beta formula 5.2% vs empirical 5.4%; FAR 1.00% quantile 23.00 (planets kept 37%); FAR 0.25% quantile 25.43 (planets kept 30%)
- strict, 4 sectors: median beta2h 1.96; sigma_beta formula 3.7% vs empirical 3.8%; FAR 1.00% quantile 22.46 (planets kept 62%); FAR 0.25% quantile 25.07 (planets kept 54%)

Tier thresholds (higher of 2 / 4 sectors, rounded up): normal FAR 1.00% = 10.8, normal FAR 0.25% = 11.6, strict FAR 1.00% = 23.0, strict FAR 0.25% = 25.5

## Validation
### 2 sectors (10000 noise, 2000 planets)
- tier shares (noise stars): normal 68.8%, strict 22.9%, tier3 8.3%; planets whose tier differs from their noise twin: 4.7%
- FAR 1.00%: T1 overall 0.10% (95% 0.05%–0.18%) → pass; T2 normal 0.13%, strict 0.04% → pass; T3 planets claimed 65.5% vs single threshold 43.67 2.1% → pass; tier-3 noise listed 4.71%
- FAR 0.25%: T1 overall 0.02% (95% 0.01%–0.07%) → pass; T2 normal 0.01%, strict 0.04% → pass; T3 planets claimed 61.1% vs single threshold 47.75 0.9% → pass; tier-3 noise listed 3.69%
### 4 sectors (10000 noise, 2000 planets)
- tier shares (noise stars): normal 70.5%, strict 20.9%, tier3 8.5%; planets whose tier differs from their noise twin: 2.5%
- FAR 1.00%: T1 overall 0.03% (95% 0.01%–0.09%) → pass; T2 normal 0.04%, strict 0.00% → pass; T3 planets claimed 80.1% vs single threshold 49.75 8.8% → pass; tier-3 noise listed 5.00%
- FAR 0.25%: T1 overall 0.00% (95% 0.00%–0.04%) → pass; T2 normal 0.00%, strict 0.00% → pass; T3 planets claimed 77.2% vs single threshold inf 0.0% → pass; tier-3 noise listed 3.97%

**Validation: PASS** (adoption is Julien's decision)
