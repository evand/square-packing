# rectfront: R_n = {(w, h) : n unit squares fit in w x h}, small n  (2026-10-20)

Code: `frontier.py` (fq rectangle mode, container [0, h + dx] x [0, h], min h; random / tilt-alphabet / band starts
+ kicks).  Runs: `runs/rectfront/f1` (n = 2..8, dx step 0.05, 40 + 40 quenches), `runs/rectfront/f3fine` (n = 3, 4,
dx 0.6..2.2 step 0.02, 150 + 150).  Lean: `lean/Sqpack/Rect.lean`.  Labels: [proved] Lean; [measured] float.

## Results

* **R_2** [proved, `rect2`]: two unit squares need max(w, h) >= 2.  Proof: with w, h < 2 every unit square
  contains the rectangle's centre in its interior (`centre_mem_sqInt`; `|dx|, |dy| < 1 - k`, `k = (|C|+|S|)/2`,
  and `2k(1 - k) <= 1/2`).
* **R_3 has a tilted piece** [measured]: two uprights in the bottom corners, a 45 deg square in the gap touching
  their inner top corners and the top wall: `h = 2 + sqrt2 - w/2`, below the grid for `2 sqrt2 < w < 3`
  (e.g. 2.942809 x 1.942809, = (2 + 2sqrt2/3) x (1 + 2sqrt2/3)).  Every sub-grid optimum the search found matches
  this formula to 1e-12; nothing else beat the grid.  Conjecture (w >= h): `h_3(w) = 2` on [2, 2sqrt2],
  `2 + sqrt2 - w/2` on [2sqrt2, 3), `1` from 3.
* **R_4** [measured]: two uprights at opposite corners and a tilted pair (39.3-43 deg, angle varying along the
  frontier, not 45) between them: below the grid for w in about (3.72, 4), e.g. 3.912287 x 1.912287 at 39.95 deg.
* **R_5** [measured]: the s(5) packing stretches: four corner uprights + a central 45 deg square touching the four
  inner corners: `w + h = 4 + sqrt2` for w in [2 + 1/sqrt2, 3) (h from 2.707 down to 2.414); at w = 3 the grid
  (h = 2) takes over.
* n = 6, 7, 8 [measured]: tilted pieces at h < 2 (strip-of-height-<2 family, x.91-ish) and n = 7 = R_5 piece + a
  column of two (w = 3.7071 at h = 2.7071, w + h = 6 + sqrt2 - ...); not yet read closely.

## Caveat
The search misses known optima often (n = 3: found h = 2 at dx = 0.84, 0.86, 0.92, 0.94 where the 45 deg family
gives less; n = 8 long thin cells worse than the trivial row).  Below-grid hits are real (checked feasible by fq);
"no hit" means little.  A continuation along the frontier (seed dx + 0.02 from dx's best) would fix most of it.

## Literature (scan `literature/rectangles-small-n-2026-10-20.md`)
No paper characterises R_n for n >= 2.  Yang (arXiv 2608.13595, G_3) settles R_3 by rescaling except 1.5 <= w < 2,
max(w, 3 - w/2) <= h < 3 (w short); the scan's "natural conjecture" (that region infeasible) is refuted by our 45 deg
family, which lies inside it.  R_4: Stenlund 2010 (Friedman's Squares in Dominoes) has 4 in 1.930 x 3.861, the
same shape as our tilted-pair family.  No rectangle result in any proof assistant.
