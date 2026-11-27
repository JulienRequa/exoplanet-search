# Threshold simulation at the measured red-noise levels (2026-11-27)

White-noise rows are taken from results/rednoise_snr_validation_2026-11-27.csv (same seeds and code; they do not depend on
the red-noise level): 2 sectors, FAR 1% at SNR 8.38 (93% of faint planets kept); 4 sectors, 8.15 (99%).
Note: the red-noise-aware SNR (snr_red) FAILED its pre-registered validation (sims/rednoise_snr_validation_results_2026-11-27.md),
so its rows are shown for information only; it is not adopted. Levels from results/rednoise_dev_2026-11-27/SUMMARY.md.

# Level: typical (dev median, alpha 0.091)

## 2 sectors, white + red (alpha 0.09, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.00 (16–84%: 1.00–1.03); on planets median 1.00
- snr_white: FAR at 7.1 = 38.8%; planets recovered at 7.1 = 97%
  - FAR 5.00%: threshold 8.43 (FAR 95% CI 3.81%–6.53%); faint planets recovered 94%
  - FAR 1.00%: threshold 9.12 (FAR 95% CI 0.54%–1.83%); faint planets recovered 91%
  - FAR 0.25%: threshold 9.77 (FAR 95% CI 0.10%–0.88%); faint planets recovered 89%
- snr_red: FAR at 7.1 = 35.7%; planets recovered at 7.1 = 96%
  - FAR 5.00%: threshold 8.33 (FAR 95% CI 3.81%–6.53%); faint planets recovered 94%
  - FAR 1.00%: threshold 9.12 (FAR 95% CI 0.54%–1.83%); faint planets recovered 91%
  - FAR 0.25%: threshold 9.58 (FAR 95% CI 0.10%–0.88%); faint planets recovered 90%


## 4 sectors, white + red (alpha 0.09, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.00 (16–84%: 1.00–1.02); on planets median 1.00
- snr_white: FAR at 7.1 = 43.4%; planets recovered at 7.1 = 100%
  - FAR 5.00%: threshold 8.54 (FAR 95% CI 3.81%–6.53%); faint planets recovered 99%
  - FAR 1.00%: threshold 9.22 (FAR 95% CI 0.54%–1.83%); faint planets recovered 98%
  - FAR 0.25%: threshold 9.57 (FAR 95% CI 0.10%–0.88%); faint planets recovered 97%
- snr_red: FAR at 7.1 = 41.1%; planets recovered at 7.1 = 100%
  - FAR 5.00%: threshold 8.52 (FAR 95% CI 3.81%–6.53%); faint planets recovered 99%
  - FAR 1.00%: threshold 9.16 (FAR 95% CI 0.54%–1.83%); faint planets recovered 98%
  - FAR 0.25%: threshold 9.55 (FAR 95% CI 0.10%–0.88%); faint planets recovered 97%


# Level: pessimistic (dev 90th percentile, alpha 0.292)

## 2 sectors, white + red (alpha 0.29, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.55 (16–84%: 1.46–1.64); on planets median 1.58
- snr_white: FAR at 7.1 = 99.5%; planets recovered at 7.1 = 92%
  - FAR 5.00%: threshold 14.49 (FAR 95% CI 3.81%–6.53%); faint planets recovered 68%
  - FAR 1.00%: threshold 16.34 (FAR 95% CI 0.54%–1.83%); faint planets recovered 60%
  - FAR 0.25%: threshold 18.52 (FAR 95% CI 0.10%–0.88%); faint planets recovered 50%
- snr_red: FAR at 7.1 = 50.6%; planets recovered at 7.1 = 86%
  - FAR 5.00%: threshold 9.55 (FAR 95% CI 3.81%–6.53%); faint planets recovered 67%
  - FAR 1.00%: threshold 11.56 (FAR 95% CI 0.54%–1.83%); faint planets recovered 51%
  - FAR 0.25%: threshold 14.22 (FAR 95% CI 0.10%–0.88%); faint planets recovered 37%


## 4 sectors, white + red (alpha 0.29, tau 1.0 h) noise (1000 noise stars, 250 faint planets)
- beta on noise stars: median 1.57 (16–84%: 1.50–1.62); on planets median 1.58
- snr_white: FAR at 7.1 = 99.9%; planets recovered at 7.1 = 96%
  - FAR 5.00%: threshold 14.37 (FAR 95% CI 3.81%–6.53%); faint planets recovered 80%
  - FAR 1.00%: threshold 17.17 (FAR 95% CI 0.54%–1.83%); faint planets recovered 74%
  - FAR 0.25%: threshold 18.36 (FAR 95% CI 0.10%–0.88%); faint planets recovered 68%
- snr_red: FAR at 7.1 = 52.8%; planets recovered at 7.1 = 94%
  - FAR 5.00%: threshold 9.38 (FAR 95% CI 3.81%–6.53%); faint planets recovered 80%
  - FAR 1.00%: threshold 11.31 (FAR 95% CI 0.54%–1.83%); faint planets recovered 70%
  - FAR 0.25%: threshold 12.79 (FAR 95% CI 0.10%–0.88%); faint planets recovered 62%

