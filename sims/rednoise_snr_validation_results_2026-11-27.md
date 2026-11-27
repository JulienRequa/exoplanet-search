# Results: sims/rednoise_snr_sim.py, rednoise_snr_validation_2026-11-27.csv (2026-11-27)

Reproducibility vs results/multisector_threshold_sim_2026-11-27.csv: 5000 matched stars, max |snr_white − snr| = 0.00e+00

## 2 sectors, white noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.00 (16–84%: 1.00–1.00); on planets median 1.00
- snr_white: FAR at 7.1 = 17.2%; planets recovered at 7.1 = 97%
  - FAR 5.00%: threshold 7.79 (FAR 95% CI 3.81%–6.53%); faint planets recovered 95%
  - FAR 1.00%: threshold 8.38 (FAR 95% CI 0.54%–1.83%); faint planets recovered 93%
  - FAR 0.25%: threshold 8.99 (FAR 95% CI 0.10%–0.88%); faint planets recovered 91%
- snr_red: FAR at 7.1 = 17.2%; planets recovered at 7.1 = 97%
  - FAR 5.00%: threshold 7.79 (FAR 95% CI 3.81%–6.53%); faint planets recovered 95%
  - FAR 1.00%: threshold 8.38 (FAR 95% CI 0.54%–1.83%); faint planets recovered 93%
  - FAR 0.25%: threshold 8.99 (FAR 95% CI 0.10%–0.88%); faint planets recovered 91%

## 2 sectors, white + red (alpha 0.50, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 2.24 (16–84%: 2.10–2.38); on planets median 2.27
- snr_white: FAR at 7.1 = 100.0%; planets recovered at 7.1 = 81%
  - FAR 5.00%: threshold 24.38 (FAR 95% CI 3.81%–6.53%); faint planets recovered 33%
  - FAR 1.00%: threshold 28.59 (FAR 95% CI 0.54%–1.83%); faint planets recovered 18%
  - FAR 0.25%: threshold 31.09 (FAR 95% CI 0.10%–0.88%); faint planets recovered 14%
- snr_red: FAR at 7.1 = 62.8%; planets recovered at 7.1 = 73%
  - FAR 5.00%: threshold 11.85 (FAR 95% CI 3.81%–6.53%); faint planets recovered 24%
  - FAR 1.00%: threshold 15.30 (FAR 95% CI 0.54%–1.83%); faint planets recovered 11%
  - FAR 0.25%: threshold 18.77 (FAR 95% CI 0.10%–0.88%); faint planets recovered 4%

## 4 sectors, white noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.00 (16–84%: 1.00–1.00); on planets median 1.00
- snr_white: FAR at 7.1 = 17.8%; planets recovered at 7.1 = 100%
  - FAR 5.00%: threshold 7.64 (FAR 95% CI 3.81%–6.53%); faint planets recovered 100%
  - FAR 1.00%: threshold 8.15 (FAR 95% CI 0.54%–1.83%); faint planets recovered 99%
  - FAR 0.25%: threshold 9.11 (FAR 95% CI 0.10%–0.88%); faint planets recovered 98%
- snr_red: FAR at 7.1 = 17.8%; planets recovered at 7.1 = 100%
  - FAR 5.00%: threshold 7.64 (FAR 95% CI 3.81%–6.53%); faint planets recovered 100%
  - FAR 1.00%: threshold 8.15 (FAR 95% CI 0.54%–1.83%); faint planets recovered 99%
  - FAR 0.25%: threshold 9.11 (FAR 95% CI 0.10%–0.88%); faint planets recovered 98%

## 4 sectors, white + red (alpha 0.50, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 2.27 (16–84%: 2.17–2.37); on planets median 2.29
- snr_white: FAR at 7.1 = 100.0%; planets recovered at 7.1 = 86%
  - FAR 5.00%: threshold 25.01 (FAR 95% CI 3.81%–6.53%); faint planets recovered 53%
  - FAR 1.00%: threshold 28.75 (FAR 95% CI 0.54%–1.83%); faint planets recovered 44%
  - FAR 0.25%: threshold 30.54 (FAR 95% CI 0.10%–0.88%); faint planets recovered 40%
- snr_red: FAR at 7.1 = 66.2%; planets recovered at 7.1 = 80%
  - FAR 5.00%: threshold 11.85 (FAR 95% CI 3.81%–6.53%); faint planets recovered 48%
  - FAR 1.00%: threshold 14.58 (FAR 95% CI 0.54%–1.83%); faint planets recovered 32%
  - FAR 0.25%: threshold 15.96 (FAR 95% CI 0.10%–0.88%); faint planets recovered 26%

## Pre-registered validation
- 2 sectors: V1 median beta (white) 1.00, recovery change +0% → pass; V2 thresholds white 8.38 vs red 15.30 (+83%) → FAIL; V3 FAR on red stars at the white threshold 30.4% → FAIL
- 4 sectors: V1 median beta (white) 1.00, recovery change +0% → pass; V2 thresholds white 8.15 vs red 14.58 (+79%) → FAIL; V3 FAR on red stars at the white threshold 36.9% → FAIL

**Validation: FAIL**
