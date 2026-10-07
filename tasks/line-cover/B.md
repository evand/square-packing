# Brief B: exact checker for mixed covers, 2026-09-27

Working tree `~/math/square-packing/public/s12`.  Read `tasks/line-cover/README.md`, `FORMAT.md`, then
`search/ZEROMARGIN.md`, `search/zeromargin.py` (CORE, P1, CHAIN/Lemma H, weighted `cert`, `--d4`), `search/S32_EXACT.md`
§1–2 and `search/S32_SHIFT.md` §1 (germ failures), `notes/s13-casefree.md` §1 (reduction), `search/M4_MARGIN.md` §1.
**Resources:** `taskset -c 10-14,26-30`, RSS `<= 30 GB`.  Budget: one working day wall-clock.

## Goal

`search/zm_mixed.py`: an exact checker for `FORMAT.md` files, **importing** `zeromargin.py` (do not edit it; its sha is
pinned).  Statement checked: every closed unit square in `[0,s]^2` has `μ(Q) >= 1`.  `search/mixed_cover.py` is written
by agent A (reader contract in `FORMAT.md`); until it lands, write against the contract.

## Design starting points (verify, don't trust)

* **Off the germs, segments and polygons are easy.**  Every point of `piece ∩ core(box)` lies in every square of the
  box, so `mass(piece ∩ core(box))` is a certified lower bound (core as in `zeromargin.py`; it is convex, an
  intersection of translates of the angle-bin core; replace its round part by an inscribed polygon to stay rational).
  Segment ∩ convex polygon is exact rational clipping; polygon ∩ polygon likewise.  The loss is `O(box size)`, and the
  weight is continuous in pose except where an edge of Q lies along a segment, so box refinement should terminate
  wherever there is margin.
* **The hard place is the axis-parallel family / tile germs.**  At `θ = 0` a grid-line segment on Q's edge counts in
  full; at `θ = 0+` its chord is a fraction depending on (edge offset)/θ.  So near germs the mass is not continuous in
  `(cx, cy, θ)`; expect to need a blow-up coordinate (offset/θ) or an exact treatment of `θ = 0` plus one-sided
  bounds, or the disjunctive/threshold reasoning of Lemma H applied to line densities (a segment is a continuum of
  pivots with no gaps between them, which is exactly what S32_SHIFT §1 defect 2 lacked).  Work out what is sound,
  write the argument down, then code it.
* Symmetry (`--d4` as in `zeromargin.py`) and multiprocessing: reuse.

## Tests (required)

* Point-only files: `zm_mixed.py` agrees with `zeromargin.py cert` on the shipped `certificates/rung2` `s(13)` cover
  (a short column sweep suffices; a full run if cheap).
* A mixed toy cover you construct at `m = 4` (e.g. the grid lines at a uniform density plus points; generous total is
  fine) that certifies, and **rejection tests**: covers with a known hole (remove / lighten a piece; shift a segment
  off a grid line) must fail at the right pose.  Float cross-check with `mixed_cover.mass_in_square_float`.
* Report box counts / CPU time vs the point checker on comparable covers, and specifically on a germ cell.
* When A has written an `m = 4` mixed cover (`runs/line-cover_m4_*`), try it (partial sweep ok).

## Report

New note `search/ZM_MIXED.md`: the soundness argument (lemmas, with proofs), what the checker does at germs, test
results, performance, known limitations.  Commit your files (not pushed).  Final reply `<= 200` words.
