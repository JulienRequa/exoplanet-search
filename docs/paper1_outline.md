# Paper 1 — outline (started 2026-11-26)

Working title: *A pre-registered, light-curve-only vetting pipeline for TESS planet candidates:
blind-test accuracy, a 2,896-TOI search, and previously unreported false-positive evidence*

Rule: paper 1 reports **frozen v2** results exactly. v3 never changes these numbers; a v3 re-run
of the hunt would be a separate result. Evidence locations are relative to the repo root.
Detailed findings live in docs/PAPER_NOTES.md (M = methods, R = results, L = limitations).

## 1. Introduction
- TESS has released >8,000 TOIs; most are unconfirmed; follow-up time is scarce, so cheap,
  reproducible vetting of light curves matters.
- What is new here: (1) a pre-registered dev → freeze → blind design for a vetting tool;
  (2) a blind accuracy measurement; (3) a complete run on all open PCs in the tested range;
  (4) candidates with strong false-positive evidence that has not been reported publicly.
- Key numbers: blind 92% planets kept, 51% FPs caught; 2,896 TOIs; FP fraction ~22%.

## 2. Related work
- TESS vetting tools: SPOC DV reports, QLP, DAVE, TRICERATOPS (Giacalone et al.),
  RAVEN (Hadjigeorghiou et al. 2025; Lafarga et al. 2026), LEO-vetter, CNN classifiers.
- TFOP follow-up (SG1–SG5) and ExoFOP dispositions.
- Kepler-era robovetter and false-positive rates (Fressin et al. 2013; Morton et al. 2016).
- TODO: literature pass (not yet done systematically).

## 3. Methods (the v2 pipeline)
- Data, cleaning, period search near the catalog period, per-sector ephemeris (M4, M6).
- Transit fit with Claret priors + Kipping q1/q2 (M3); multi-start (M7).
- Eight checks and verdict logic: VETTING_RULES.md (frozen); red-noise odd/even (M2);
  proper-motion centroid (M5); neighbour-depth check.
- Evidence: VETTING_RULES.md, pipeline.py at tag v2-freeze (46b253d), DEVLOG.

## 4. Validation
- Design: dev samples (293 stars) → freeze → blind (200 stars, locked in 0217ed9) (M1).
- Dev results (method development only): 121/124 kept, 67/127 caught (VETTING_RULES).
- **Blind (headline)**: 91/99 planets kept (92%), 43/85 FPs caught (51%), 0 crashes.
  Losses: 5 odd/even, 2 density, 1 size (DEVLOG blind_v2 entry).
- Evidence: results/blind_v2_results_2026-11-25.csv, data/blind_v2.csv (md5 66805c4c…).

## 5. Hunt results
- 2,896 PC TOIs (P 0.5–10 d, depth > 2000 ppm), list pushed before the run (md5 684506ef…).
- Tested 2,397; no 2-min data 499 (all rechecked); verdicts 998 / 984 / 414 / 1 (R3).
- FP fraction ≈ 22% (68% 14–26%), formula + assumptions (R4).
- 95 failed only odd/even (v3 recheck list).
- Clean candidates 854 (96 flagged > 1.8 R_Jup).
- Evidence: results/hunt_v2_2026-11-25.csv (md5 e4921968…), results/hunt_v2_compare_2026-11-26/.

## 6. Flagged candidates (headline)
- Listing rule written after the hunt, cutoffs only from dev + blind (M8): 105 strong fakes.
- Comparison with ExoFOP follow-up, SPOC DV and RAVEN: 14 already resolved, 47 independently
  flagged, 44 new; 11 with no ground follow-up.
- Diagnostic report cards (variability / phase curve; nearby EBs), difference images for
  survivors; final list = survivors with evidence summaries.
- Evidence: results/hunt_v2_compare_2026-11-26/strong_fakes_followup.csv,
  results/report_cards_2026-11-26/ (in progress).
- Note: a report-card result is a diagnostic, not a v2 change; frozen v2 verdicts are
  reported unchanged alongside.

## 7. Limitations
- L1–L4 (neighbours < 1 px, active stars, shallow transits, TTVs); faint planets missed;
  tested range only; the v2 secondary-eclipse test uses white-noise errors (report cards show
  variability can mimic it, e.g. TOI 5684.01); FP-fraction assumptions (R4).

## 8. AI-use disclosure
- Code written with an AI coding assistant (Claude) under the author's direction; the author
  set goals, physics rules, thresholds, and every decision; all decisions and changes are
  logged in DEVLOG with dates. Draft wording to be finalised.

## Data and code availability
- GitHub repo + Release `hunt-v2-results-2026-11-26` + Zenodo DOI (to add).
