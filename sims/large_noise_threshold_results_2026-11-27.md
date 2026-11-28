# Results: sims/large_noise_threshold_sim.py (2026-11-27)

Rows: 44000; failed searches (SNR 0): 0

## 2 sectors, white only (reference): 10000 noise stars, 1000 faint planets
- false alarms at SNR 7.1 (v2): 16.7%; beta on noise stars: median 1.00, 90th pct 1.00
- FAR 5.00%: threshold 7.71 (bootstrap 95% 7.68–7.75); 500 noise stars above (FAR 95% CI 4.59%–5.44%); faint planets kept 96.0%
- FAR 1.00%: threshold 8.47 (bootstrap 95% 8.39–8.56); 100 noise stars above (FAR 95% CI 0.82%–1.21%); faint planets kept 93.7%
- FAR 0.25%: threshold 9.16 (bootstrap 95% 8.96–9.33); 25 noise stars above (FAR 95% CI 0.17%–0.37%); faint planets kept 90.9%

## 2 sectors, white + red (typical, alpha 0.091): 10000 noise stars, 1000 faint planets
- false alarms at SNR 7.1 (v2): 40.1%; beta on noise stars: median 1.00, 90th pct 1.04
- FAR 5.00%: threshold 8.43 (bootstrap 95% 8.38–8.48); 500 noise stars above (FAR 95% CI 4.59%–5.44%); faint planets kept 94.6%
- FAR 1.00%: threshold 9.22 (bootstrap 95% 9.12–9.34); 100 noise stars above (FAR 95% CI 0.82%–1.21%); faint planets kept 90.8%
- FAR 0.25%: threshold 9.96 (bootstrap 95% 9.70–10.17); 25 noise stars above (FAR 95% CI 0.17%–0.37%); faint planets kept 87.3% → enough false alarms: yes

## 4 sectors, white only (reference): 10000 noise stars, 1000 faint planets
- false alarms at SNR 7.1 (v2): 19.5%; beta on noise stars: median 1.00, 90th pct 1.00
- FAR 5.00%: threshold 7.79 (bootstrap 95% 7.75–7.82); 500 noise stars above (FAR 95% CI 4.59%–5.44%); faint planets kept 99.3%
- FAR 1.00%: threshold 8.45 (bootstrap 95% 8.36–8.54); 100 noise stars above (FAR 95% CI 0.82%–1.21%); faint planets kept 98.9%
- FAR 0.25%: threshold 9.08 (bootstrap 95% 8.95–9.18); 25 noise stars above (FAR 95% CI 0.17%–0.37%); faint planets kept 98.8%

## 4 sectors, white + red (typical, alpha 0.091): 10000 noise stars, 1000 faint planets
- false alarms at SNR 7.1 (v2): 43.2%; beta on noise stars: median 1.00, 90th pct 1.03
- FAR 5.00%: threshold 8.52 (bootstrap 95% 8.47–8.56); 500 noise stars above (FAR 95% CI 4.59%–5.44%); faint planets kept 99.1%
- FAR 1.00%: threshold 9.27 (bootstrap 95% 9.17–9.36); 100 noise stars above (FAR 95% CI 0.82%–1.21%); faint planets kept 98.3%
- FAR 0.25%: threshold 9.95 (bootstrap 95% 9.74–10.19); 25 noise stars above (FAR 95% CI 0.17%–0.37%); faint planets kept 97.3% → enough false alarms: yes

## Proposed v3 threshold (higher of 2 / 4 sectors at the typical level, rounded up to 0.1)
- FAR 5.00%: SNR 8.6
- FAR 1.00%: SNR 9.3
- FAR 0.25%: SNR 10.0

**0.25% rests on enough false alarms: YES**

Adoption is Julien's decision.
