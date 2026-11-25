# Ephemeris simulation results (sims/ephemeris_sim.py), 2026-11-24

Two TESS-like sectors; 'lost' = final ephemeris misses sector-2 transits by more than half a
transit duration. Current = pipeline multi-sector refit (window 5× 1-sector error).
Per-sector = transit time per sector (sector 2 found with BLS ±2% around P₁) + straight line.

## easy
- Gap 700 d, N=30: sector-2 transits LOST — current 0%, per-sector ephemeris 0%
- 1-sector ephemeris alone drifts by more than half a duration in 0% (median 0.2 h)

## hard
- Gap 1100 d, N=30: sector-2 transits LOST — current 40%, per-sector ephemeris 0%
- 1-sector ephemeris alone drifts by more than half a duration in 73% (median 3.5 h)

easy: P 2–9 d, Rp/Rs 0.08, 1000 ppm white noise, 700 d gap
hard (TOI-2207-like): P 6–9 d, Rp/Rs 0.067, 1000 ppm white + 1000 ppm red (1 h), 1100 d gap
