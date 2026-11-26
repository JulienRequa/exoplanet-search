# Report cards — first 21 flagged candidates (2026-11-26)

Diagnostic only: frozen v2 verdicts are unchanged. Script: scripts/report_cards.py (rules in its
header, written before running). Cards: TOI_*.png; numbers: report_cards.csv.
The 21 = 11 strong fakes with no ground follow-up + top 10 "new" by v2 strength.

| Assessment (rule) | TOIs |
|---|---|
| survives both | 6490.01, 3333.01, 5898.01, 2861.01, 2904.01, 5460.01, 6435.01 |
| contamination | 4802.01, 4842.01, 6569.01, 3115.01 |
| variability | 5060.01, 3079.01, 4976.01, 4868.01, 5684.01, 5871.01, 2738.01, 3188.01, 4871.01 |
| unclear | 1538.01 (too few full transits for the local odd/even test) |

Notes from looking at the cards (judgement, not rule changes):
- 4802.01: a Gaia DR3 eclipsing binary 11″ away has the same period (1.7171 d), so the signal very
  likely comes from that catalogued EB. The false-positive evidence is real, but it is already
  public (Gaia DR3).
- 4842.01, 6569.01, 3115.01: the secondary eclipse holds up against red noise (4.8σ, 5.9σ, 6.9σ
  locally). The rule flagged them only because a known EB with an unrelated period (0.38 d,
  0.73 d, 0.73 d) sits 53–118″ away. An EB at an unrelated period cannot add a dip at phase 0.5
  coherently across two sectors, so these look like survivors. Decision for Julien.
- 6435.01: the odd transits show no dip at all, and only the even transits are visible. That
  points to a true period twice the catalog value (12.23 d), not to an eclipsing binary. Its
  odd/even flag is therefore not false-positive evidence (it could still be a planet at 12.2 d).
  Move to "unclear (period alias)". Decision for Julien.
- Targets already listed in a variable-star catalog (so not "new"): 5684.01 and 4976.01 (VSX:
  RS CVn variables), 5871.01 (VSX: eclipsing cataclysmic variable, P = 0.1469 d), 5460.01 (VSX:
  pulsator; its flag is a centroid shift, which pulsation doesn't affect).
- Pattern: 9 of the 13 secondary-eclipse flags here do not survive the red-noise, local-baseline
  test. v2's secondary check uses white-noise errors and a wide baseline, and stellar variability
  near P or P/2 can mimic it. This goes into limitations (paper 1) and the v3 list.

## Update — all 51 cards + nearby-EB timing test (2026-11-26)
Cards for the remaining 30 "new" strong fakes were added (same script, same rules). Original
rule outcomes stay in report_cards.csv (`assessment`). Julien's decisions are applied in
report_cards_final.csv (`assessment_final`, `final_note`); test: scripts/neb_timing_test.py →
neb_timing.csv.

| Final | n | TOIs |
|---|---|---|
| survives both | 22 | 2344.01, 2861.01, 2904.01, 2990.01, 3169.01, 3200.01, 3265.01, 3333.01, 3395.01, 3709.01, 4161.01, 4273.01, 4694.01, 4714.01, 4842.01, 4855.01, 5460.01, 5898.01, 6045.01, 6490.01, 6527.01, 6569.01 |
| variability | 20 | see report_cards_final.csv |
| contamination | 6 | 4802.01 (Gaia EB, same period, 11″), 3115.01, 6583.01, 3179.01, 3183.01, 3210.01 |
| unclear | 3 | 1538.01, 6435.01 (likely period alias, P ≈ 12.23 d), 7594.01 |

Timing test results: 4842.01 and 6569.01 → survives (the nearby EB's catalogued ephemeris puts
its eclipses on only 20–50% of our phase-0.5 windows). 3115.01 stays contamination.
Caveats (judgement, for Julien):
- The test counts a window as "hit" whenever an EB eclipse overlaps it. For contact-type binaries
  Gaia's eclipse widths are wide (±0.17 in phase each), so "can explain" is easy to reach and is
  weak evidence. A sharper test would fold the EB's light-curve model at our period and check for
  a local dip at phase 0.5.
- 6583.01, 3179.01, 3183.01, 3210.01 are flagged on planet size, not on a secondary eclipse. The
  card rule labels any nearby EB with an unrelated period as contamination. The phase-0.5 timing
  test does not address size flags.
- Targets already in a variable-star catalog: 5684.01, 4976.01 (RS CVn), 5871.01 (eclipsing CV),
  5460.01 (pulsator; centroid flag), 3200.01 (VSX "ROT", a rotational variable; survives).
