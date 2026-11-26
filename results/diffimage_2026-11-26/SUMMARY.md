# Difference images — 22 report-card survivors (2026-11-26; diagnostic, v2 unchanged)

Script: scripts/diffimage_check.py (method and rules in its header). Per-star images:
TOI_*_diffimage.png; numbers: diffimage.csv. Caveat: the method's real-data positional floor is
not yet measured (v3 go/no-go). The "rel. OOT" offsets in diffimage.csv are not reliable: in
crowded stamps the out-of-transit fit often locks onto a brighter star. Use the offset from the
target position.

| v2 flag type | Result | TOIs | Reading |
|---|---|---|---|
| secondary eclipse | on target | 4842.01, 4694.01, 3200.01, 3265.01, 3709.01, 3395.01, 4855.01 | The dimming is on the target, so the secondary eclipse belongs to the target itself: evidence that the target is an eclipsing binary. Strongest group. |
| secondary eclipse | inconclusive | 6569.01 | Source 0.37 ± 0.12 px from the target (3.1σ); nearest star is the target. |
| planet size | on target | 4273.01, 6527.01, 2990.01, 3169.01 | On target; the size flag stands. But 6527.01, 2990.01 and 3169.01 have b = 1.05–1.09, close to the grazing degeneracy that cost v2 two real planets. Weaker group. |
| centroid shift | on target | 6490.01, 5460.01, 6045.01 | Disagrees with v2's centroid check, which put the source about 1 px off. Unresolved until the method's positional floor is known. |
| centroid shift | no detection | 5898.01 | Difference image too weak (S/N < 5). |
| neighbour depth | source on target and on v2's neighbour | 3333.01, 4161.01 | v2's suspected neighbour sits within ~0.2 px of the target, too close to separate. Not decisive. |
| neighbour depth | inconclusive / no detection | 2904.01; 2861.01, 2344.01, 4714.01 | Too weak or ambiguous. |
