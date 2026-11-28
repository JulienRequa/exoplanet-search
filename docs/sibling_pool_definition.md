# Sibling-search candidate pool (download only) — defined 2026-11-28, before any download

Julien's instruction 2026-11-27: light curves (2-minute SPOC, no pixel files) for a broad pool of
candidate systems for the sibling search: every TESS host star with a confirmed planet or PC that has
at least 6 sectors of 2-minute data. No analysis, no searching, nothing printed about individual stars.

- Source: data/toi_2026-11-26.csv (dated snapshot).
- Host: a star with at least one TOI whose TFOPWG disposition is CP, KP or PC (TOI hosts only; non-TOI
  confirmed planets, e.g. Kepler/K2 hosts outside the TOI list, are not included).
- Sectors: pre-filter ≥ 6 sectors in the snapshot's "Sectors" column (any cadence); the exact rule
  (≥ 6 sectors of 2-min SPOC data, one product per sector) is applied at download time from the MAST
  listing; stars with fewer are recorded "below6" and nothing is downloaded.
- Excluded: all blind_v2 and blind_v3 stars (blind samples untouched by v3 work).
- Order (Julien): confirmed-planet hosts (any CP/KP TOI) first, then PC-only hosts; within each, most
  listed sectors first; then TIC.
- List: data/sibling_pool_2026-11-27.csv — 1,356 stars (288 CP/KP hosts, 1,068 PC-only); median 9
  listed sectors; 19,826 listed sector-files in total (upper bound on downloads; 4,557 for CP/KP hosts).
- Download: scripts/prefetch_sibling_pool.py fetch — every 2-min SPOC light-curve file per qualifying
  star, 120 s time limit per file, 3 attempts per file, restarts by a watchdog if the log stalls 30 min;
  stops when free disk < 20 GB. Per-star status in results/sibling_pool_prefetch_2026-11-27/status.csv;
  the log holds counts only.
- This pool is a download list, not a target selection for paper 2: the sibling-search design
  (docs/design_sibling_search.md) still chooses its targets and rules before any search.
