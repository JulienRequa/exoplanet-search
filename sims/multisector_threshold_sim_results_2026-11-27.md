# Results: sims/multisector_threshold_sim.py (2026-11-27)

## 2 sectors, white noise (1000 noise stars, 250 faint planets)
- false-alarm rate at SNR 7.1 (v2 threshold): 17.2%
- FAR 5.00%: threshold SNR 7.79 (FAR 95% CI 3.81%–6.53%); faint planets recovered 95%
- FAR 1.00%: threshold SNR 8.38 (FAR 95% CI 0.54%–1.83%); faint planets recovered 93%
- FAR 0.25%: threshold SNR 8.99 (FAR 95% CI 0.10%–0.88%); faint planets recovered 91%
- faint planets recovered at SNR 7.1: 97%; median runtime 1.5 s/star

## 2 sectors, white + red noise (1000 noise stars, 250 faint planets)
- false-alarm rate at SNR 7.1 (v2 threshold): 100.0%
- FAR 5.00%: threshold SNR 24.38 (FAR 95% CI 3.81%–6.53%); faint planets recovered 33%
- FAR 1.00%: threshold SNR 28.59 (FAR 95% CI 0.54%–1.83%); faint planets recovered 18%
- FAR 0.25%: threshold SNR 31.09 (FAR 95% CI 0.10%–0.88%); faint planets recovered 14%
- faint planets recovered at SNR 7.1: 81%; median runtime 1.3 s/star

## 4 sectors, white noise (1000 noise stars, 250 faint planets)
- false-alarm rate at SNR 7.1 (v2 threshold): 17.8%
- FAR 5.00%: threshold SNR 7.64 (FAR 95% CI 3.81%–6.53%); faint planets recovered 100%
- FAR 1.00%: threshold SNR 8.15 (FAR 95% CI 0.54%–1.83%); faint planets recovered 99%
- FAR 0.25%: threshold SNR 9.11 (FAR 95% CI 0.10%–0.88%); faint planets recovered 98%
- faint planets recovered at SNR 7.1: 100%; median runtime 2.7 s/star

## 4 sectors, white + red noise (1000 noise stars, 250 faint planets)
- false-alarm rate at SNR 7.1 (v2 threshold): 100.0%
- FAR 5.00%: threshold SNR 25.01 (FAR 95% CI 3.81%–6.53%); faint planets recovered 53%
- FAR 1.00%: threshold SNR 28.75 (FAR 95% CI 0.54%–1.83%); faint planets recovered 44%
- FAR 0.25%: threshold SNR 30.54 (FAR 95% CI 0.10%–0.88%); faint planets recovered 40%
- faint planets recovered at SNR 7.1: 86%; median runtime 2.6 s/star

## Notes (2026-11-27)
- The search's SNR is the white-noise SNR used by v2 (depth ÷ out-of-transit point scatter × √N
  in-transit points). It does not account for time-correlated (red) noise.
- White noise: a fixed threshold works. SNR 7.1 gives ~17–18% false alarms; ~8.4 (2 sectors) or
  ~8.2 (4 sectors) gives 1%, losing only a few percent of faint planets.
- White + red noise: the threshold moves from ~8 to ~25–30, because the SNR above treats
  correlated noise as if it averaged down. The red noise simulated here (AR(1), 1 h, 50% of the
  white noise per point) is a deliberately strong, assumed level; it was not measured on real TESS
  light curves. So the red-noise thresholds are not yet usable as a v3 number.
- Implication: one fixed SNR number cannot serve both cases. Options (not chosen; Julien's call):
  (a) measure the real red-noise level on dev-data light curves and re-run with it;
  (b) make the v3 detection SNR red-noise aware (noise measured in transit-duration bins,
      Pont et al. 2006) and re-run this simulation with that SNR.
- Seeds and per-star results: results/multisector_threshold_sim_2026-11-27.csv (5,000 rows).
