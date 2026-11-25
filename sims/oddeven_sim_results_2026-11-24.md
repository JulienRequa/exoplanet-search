# Odd/even simulation results (sims/oddeven_sim.py), 2026-11-24

Two TESS-like sectors (2 × 12.5 d orbits each, 30 d apart), Rp/Rs = 0.1, a/Rs = 9,
P = 3.3–4.3 d, 1000 ppm white noise; red noise = AR(1), 1 h correlation, 1000 ppm.
FAIL = σ ≥ 3.

## Before (pipeline odd/even = original white-noise rule) vs prototype v3
| case | original | v3 |
|---|---|---|
| planet, white noise (N=60) | 0.0% | 0.0% |
| planet, red noise (N=60) | **56.7%** | 0.0% |
| EB, odd 15% shallower, white (N=30) | 100% | 100% |
| EB, odd 15% shallower, red (N=30) | 100% | 63.3% |

## After v3 was implemented (pipeline odd/even = v3): identical columns
| case | pipeline | prototype v3 |
|---|---|---|
| planet, white | 0.0% | 0.0% |
| planet, red | 0.0% | 0.0% |
| EB white | 100% | 100% |
| EB red | 63.3% | 63.3% |
Spot checks: pipeline σ = prototype σ to 4 decimals (0.4526, 19.7694, 0.1985).
