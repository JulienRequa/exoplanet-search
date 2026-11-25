# Exoplanet search pipeline — notes for Claude

Owner: Julien (10th grade, learning Python). Explain every change in plain language.

## Ground rules
1. Do NOT change physics or math without asking first: `transit_model`, density
   calculations, the fit setup (`fit_transit` bounds/guesses), `profile_b`, and the vetting
   checks and their thresholds. Julien writes and owns these. If one looks wrong, explain
   why and let him decide.
2. Vetting thresholds are frozen (see VETTING_RULES.md once it exists). Never tune them to
   make a result look better. Never re-score an old blind sample with new code and call it
   an improvement — improvements are measured only on a fresh sample.
3. Refactoring must not change results. Verify against the regression numbers below.
4. Small steps. Explain each change, then commit ("Add X", "Fix Y").
5. Ask before installing anything, deleting files, or pushing to GitHub.
6. Reproducibility: fixed random seeds; save dated snapshots of downloaded catalogs
   (e.g. the TOI table in data/) that results depend on.
7. Keep DEVLOG.md updated each session (what we did and why).

## Environment
`conda env create -f environment.yml` then `conda activate exoplanet`.

## Regression results (published / Colab)
- Known-planet batch (n_sectors=2): 13/13 recovered, every period matches published,
  12 planets 5/5, pi Men 4/4 with density WEAK.
- WASP-18, 4 sectors: P 0.94145 d, Rp/Rs ≈ 0.0992, a/Rs ≈ 3.418, b ≈ 0.406,
  inc ≈ 83.17°, density ≈ 0.85 g/cm³, 5/5 PLANET CANDIDATE.
- WASP-121, 2 sectors: Rp/Rs ≈ 0.1208, b flagged → b < 0.19, inc > 87.1°.
- Local (this Mac) baseline of the 13-planet batch: results/baseline_known_2026-11-24.csv
