# Mixed cover format v1 (points + segments + polygons)

A *mixed cover* is a finite measure on the container `[0,s]²` made of point masses, masses spread uniformly
along segments, and masses spread uniformly over convex polygons.  It extends the plain point format of
[`../FORMAT.md`](../FORMAT.md); a plain point file is a mixed file with no segments and no polygons.
`s21_mixed_cover_5.txt` is in this format (points and axis-parallel segments only).

Plain text, whitespace-separated integers, `#` comments to the end of a line.  Everything is exact.

```
mixed 1
s_num s_den           # container [0, s]^2 with s = s_num/s_den
D                     # coordinate denominator: every coordinate is X/D
W                     # mass denominator: every mass is w/W
np                    # number of point masses
X Y w                 # np lines: the point (X/D, Y/D) carries mass w/W
ns                    # number of segments
X0 Y0 X1 Y1 w         # ns lines: mass w/W spread uniformly (by length) on the closed segment (X0,Y0)-(X1,Y1)
npg                   # number of polygons
k w X1 Y1 ... Xk Yk   # npg lines: convex polygon, vertices counter-clockwise, mass w/W spread uniformly (by area)
```

Masses are **totals**, so the mass a square captures is rational: a segment contributes
`w/W · |seg ∩ Q| / |seg|` (a parametric fraction of the segment), a polygon `w/W · area(P ∩ Q) / area(P)`.
Degenerate pieces (a zero-length segment, a zero-area polygon) are invalid; use a point instead.  All pieces lie
in `[0,s]²` and all `w ≥ 0`.  The total mass of the file is the sum of all the `w/W`.

## What a file asserts

For every **closed** unit square `Q ⊆ [0,s]²`, at every centre and every angle, `μ(Q) ≥ 1`, where `μ` is the
sum of the point, segment and polygon measures.  Closed matters: a point on the boundary of `Q` counts, and a
segment lying along an edge of `Q` counts in full.

## Why that bounds `s(n)`

If the file's total mass is `< n`, then `s(n) ≥ s`.  Suppose `n` unit squares with disjoint interiors fit in a
square of side `s' < s`.  Scale the configuration about the container's centre by `s/s' > 1`: the squares become
squares of side `s/s' > 1` with disjoint interiors inside `[0,s]²`, and the concentric closed unit square inside
each of them is contained in its interior.  So we get `n` pairwise **disjoint** closed unit squares
`Q_1, …, Q_n` in `[0,s]²`, and

```
n ≤ Σ μ(Q_i) = μ(Q_1 ∪ … ∪ Q_n) ≤ μ([0,s]²) = total < n,
```

a contradiction.  (Closed squares are measurable; the argument needs nothing else about `μ`.)  This is proved in
Lean for arbitrary measures: `packing_le_measure` and `not_packs_of_measure` in `lean/Sqpack/MixedMeasure.lean`.

## Checkers

* `search/mixed_cover.py`: reader and well-formedness check (`load`, `validate`, `total`), plus a float evaluator.
* `search/zm_mixed.py`: exact checker (points, segments, polygons; `--cert-mode` refuses polygons).
* `verify2/src/bin/zmx2.rs` (`zmx2`): independent exact checker (points and axis-parallel segments; anything
  else is refused with `ERROR: unsupported`).
