# Brief A: mixed covers, cover side (heuristic), 2026-09-27

Working tree `~/math/square-packing/public/s12`.  Read `tasks/line-cover/README.md`, `FORMAT.md`, then
`search/S21_COVER.md` §0, `search/M4_MARGIN.md` §0, `search/S32_COVER.md` §0 and §3 (where the weight sits),
`search/S32_SHIFT.md` §0–1.  Tools to reuse/extend: `closed4.py` (LP, rows, colgen), `rung2_close.py`, `close_shifted.py`,
`s21_cover_eval.py` (`hardscan`, `family` = the strict protocol), `family_rows.py`.
**Resources:** `taskset -c 0-9,16-25`, RSS `<= 60 GB` total.  Budget: one working day wall-clock.

## Question

Does moving cover weight from discrete points to uniform densities on segments (mainly the grid lines) and convex
polygons shrink the validity overhead, i.e. the gap between the converged LP and the honest cost?  At `m = 5`, the
point-cover LP is `≈ 20.73` and the best honest cost `21.016` (`×1.0138`); we need `< 20.947`.

## Do

0. **First hour:** write `search/mixed_cover.py` per the contract in `FORMAT.md` (B imports it), commit it early.
1. Float machinery: row coefficient of a segment / polygon column = its mass in the pose's closed square (clipping;
   closed semantics: a segment on an edge of Q counts fully — at exact tile poses this matters).  Strict-protocol
   evaluator for mixed covers (a mixed-cover analogue of `hardscan` + `family`, same pitches, plus polish).
2. **m = 4 calibration.**  Column families to try (your judgement; start simple): grid-line segments of length
   `1/q` (piecewise-constant density, e.g. `q = 10–50`), optionally off-grid segments where point covers put
   curves, small polygons, plus the usual point columns.  LP + closing loop; report LP value, strict honest cost,
   where the worst poses are, vs the point-cover numbers (LP `12.34`, honest `12.58–12.70`).
   Key check: do the near-tile holes at `θ ≲ 1/D` disappear?
3. **m = 5**: same, warm-started from the best `m = 4` structure / the `m = 5` point LP (`S21_COVER.md`).  Close the loop
   with permanent dip rows (the `s(32)` recipe) until the strict minimum at pitch 0.002 is stable (`< 0.1 %` over
   3 rounds), then a pitch-0.001 confirmation scan.
4. Write the best `m = 5` cover as a mixed file, scaled to `1.0025 ×` over its confirmed minimum:
   `runs/line-cover_m5_candidate.txt` (and the `m = 4` one).  Also report the piece counts (vs 3–13k points).

Early stop: if at `m = 4` the mixed honest cost is not clearly better than the point cover's (`≈ 1.5 %+` over LP), say so
after a solid attempt and stop; that is a valid answer.

## Report

New note `search/LINE_COVER.md` (numbers first; tables of LP / honest cost / piece counts at `m = 4, 5`; where the
residual dips sit; gate verdict **GO / NO-GO**).  Commit your files (not pushed).  Final reply `<= 200` words, numbers first.
