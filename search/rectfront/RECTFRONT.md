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

## R_3 lower bound: method and status (10-20/21)

Claim: three unit squares in [0,w] x [0,h], w < 2, h < 3 force 2w + h >= 4 + 2sqrt2.  With rect2 and the
constructions (rect3_wedge, grid, column) this gives R_3 exactly (w <= h): fits iff w >= 1 and (h >= 3 or w >= 2 or
2w + h >= 4 + 2sqrt2).

* Certificate (`r3bb.py`): SAT for squares makes each configuration one of 512 separation types (per pair, an edge
  normal of one square).  For weights lambda >= 0 on the pairs, with net pair forces G_K on the squares
  (sum G_K = 0), LP duality for the centres gives  sum lambda S <= sum_K |G_Kx|(w/2 - k_K) + |G_Ky|(h/2 - k_K).
  The right side is linear increasing in (w, h), so checking it fails at the two upper corners V1 = (2, 2sqrt2),
  V2 = ((1 + 2sqrt2)/2, 3) of the excluded region excludes the open region (margin >= 0 suffices).
* Tilt map (`tiltlp.py`, `tiltmap.py`, 5 deg grid, w0 = 1.95, 1.99): no tilt triple beats the boundary; the only
  exact tie is the wedge (two upright, one at 45).  But margin 0 holds on whole lines: two tilts = 0, third free
  (two upright squares side by side, tight at V1, w = 2).  So the tilt-space search has a tube around three lines.
* On the line, wedge types: weights (cos t, sin t) cancel the middle square's vertical force and give exactly
  w + ((p^2 - 1)/2) h >= p + p^2, p = cos t + sin t; at V1, V2 the margins p + p^2(1 - sqrt2) - 2 + sqrt2 and
  p - p^2/2 + 1 - sqrt2 are concave in p with zeros only at p = 1, sqrt2 (exact).
* Global search `runs/rectfront/r3bb2` (symmetry-reduced: labels sorted, mirror): 86,416 boxes certified, 423 tube
  boxes set aside, **no failure outside the tube**.  Float64 with an explicit 1e-9 margin over padded trig
  enclosures (exact-arithmetic replay still to do).
* Tube (`r3local.py`, `iad.py`): tilt-dependent weights (support and zero force components of the LP optimum at a
  probe; secondary weights proportional to the distance from the line), directed-rounding interval AD along rays
  from the line, sector/rho/beta subdivision.  Works on the side-by-side and wedge types tested; too slow on the
  raw tube boxes, so first `r3line1`: refine the tube, keeping only boxes that touch a zero line.

## 2026-10-21: chain reduction for R_3 (replaces the tube/origin interval search)

Tube run r3tube6 stopped at 42/423 (0 failures) in favour of this.  Scripts: `chain/`.

Setup: strip of height H < 2, length L; three unit squares sorted by centre x (a, b, c).
1. Wall clearance: x_a >= k_a, x_c <= L - k_c, k = (|cos| + |sin|)/2 (Lean: `half_extent`).
2. Pair separation: disjoint interiors give a separating edge normal u; u_x <= 0 would force |dy| >= 1,
   impossible since |dy| <= H - 1 < 1.  So dx >= (h_A(u) + h_B(u) - u_y dy)/u_x.
3. The resulting bound on L is linear in the centre heights, so b sits on the floor or ceiling
   (vertical mirror: floor).  The sides then decouple: L >= P(beta) + P(-beta), where P(beta) is the
   least x_b with b on the floor and some square to its left.
4. Numerically (lemP3.py) the neighbour can be taken upright, giving the closed form
   P(t) = min(1+k, U, V'),  P(-t) = min(1+k, V, U')  (t in [0,45], c=cos t, s=sin t, k=(c+s)/2):
     U  = 1/2 + (k+1/2 - s(k-1/2))/c        V  = 1/2 + (k+1/2 - s(H-1/2-k))/c
     U' = 1/2 + (k+1/2 - c(k-1/2))/s        V' = 1/2 + (k+1/2 - c(H-1/2-k))/s
   Closed form matches brute force over all neighbour tilts to 1.6e-15.
5. Target L >= 2 m(H), m = min(3/2, 2+sqrt2-H).  Splits into 9 term pairs, each linear in H on
   either side of the kink H = 1/2+sqrt2, so only H in {max(1,c+s), 1/2+sqrt2, 2} need checking
   (combos.py): all >= 0, zeros only at t=0 (all upright) and t=45 with pair (V', V) (the wedge).

Grid evidence only so far.  To prove: (i) the upright neighbour is worst (step 4); (ii) the 1-variable trig
inequalities in step 5; (iii) step 2 needs the separating-axis direction (disjoint => separating normal), or
in Lean its contrapositive with an explicit common interior point.
