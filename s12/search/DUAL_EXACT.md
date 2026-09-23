# Exact (rational) certificate:  `L(399/100) >= 12`

Code: `search/dual_exact.py`.  Certificate: `runs/dual_exact_3.99_support.txt`
(logs `runs/dual_exact_3.99.log`, `runs/dual_exact_3.99_check_F.log`, `runs/dual_exact_3.99_check_full.log`;
numbers in `runs/dual_exact_3.99.json`).  This makes the float-with-margins bound of `search/DUAL.md`
exact: every load-bearing quantity is a Python integer or `Fraction`; floats only *choose* the masses.

## Statement proved

Let `t = 399/100` and `C = [0, t]^2`.  There is an explicit measure `mu` on placements of **closed** unit
squares inside `C` — 156 poses `(p, q, cx, cy, m)` with rational rotation `theta = 2 arctan(p/q)` and
rational centre, each pose carrying `m/8` on each of its 8 dihedral images (1248 squares) — with

    mass(mu)  =  6004115377 / 500000000        =  12.008230754
    M         =  max_{x in C} mu({Q : x in Q})  =  7999999981 / 8000000000  =  0.999999997625   (exactly)
    L         =  mass / M                       =  96065846032 / 7999999981  =  12.00823078252...

Hence `mu / M` is a feasible fractional packing and the fractional packing number satisfies
`nu_f(t) >= L > 12` for every `t >= 399/100` (a packing measure for `C` is one for any larger container).
By weak LP duality (`search/CEILING.md`) **no weighted unavoidable set (LP cover) of total weight `< 12`
exists at any container side `>= 3.99`**, whatever point set, cell decomposition or verifier is used.
`L - 12 = 65846260 / 7999999981 = +8.23e-3`.

## What is exact, and how

**Poses.**  Each pose of `runs/dual_PA2_support.txt` (172 float poses, `DUAL.md`) was re-snapped:
`p/q = round(tan(theta/2) * 10^5) / 10^5` (reduced), so `cos = (q^2-p^2)/(q^2+p^2)`,
`sin = 2pq/(q^2+p^2)` are rational (`theta = 0` gives `p = 0`, exactly axis-aligned); the centre was
rounded to `10^-6` and then **clamped exactly** into `[w/2, t - w/2]` with `w/2 = (|cos|+|sin|)/2`
rational, so the closed square lies in the closed container.  Maximum displacement from the float pose:
`1.0e-5` rad in angle, `< 1e-6` in the centre.  The 8 dihedral images `(cx,cy,+theta), (t-cx,cy,-theta),
(cx,t-cy,-theta), (t-cx,t-cy,+theta), (cy,cx,-theta), (t-cy,cx,+theta), (cy,t-cx,+theta), (t-cy,t-cx,-theta)`
(the convention of `packing_dual.py` / `nu_f.max_coverage`) are formed exactly (`-theta` is `p -> -p`).
Admissibility is checked exactly for every image: the four corners `C + R(±1/2, ±1/2)`, written over the
common denominator `2 (q^2+p^2) Dk`, satisfy `0 <= X` and `100 X <= 399 D` (same for `Y`).

**Arrangement vertices.**  Each square contributes 4 closed edges with integer endpoints over its
denominator.  For every pair of non-parallel edges from different squares the intersection point of the
two lines is computed in integers (`t = cross(P2-P1, d2)/det`, `u = cross(P2-P1, d1)/det`), kept iff
`0 <= t <= 1` and `0 <= u <= 1` (closed segments; `det = 0` means parallel — no vertex, and the endpoints
of any collinear overlap are corners), and stored as a reduced triple `(X, Y, D)`.  Corners are added.
The only pruning is an exact grid prefilter (cells `floor(coord * 5)` computed by integer division of
the rational coordinate; two edges whose cell ranges are disjoint have disjoint bounding boxes), and each
pair is processed exactly once (in the cell `max` of the two min-cells, which lies in both ranges whenever
the ranges meet).  In the fundamental domain `F = {0 <= x <= y <= t/2}` the certificate has
**104,947 distinct vertices** (104,352 edge intersections + corners); over the whole container
**836,576** (831,816 intersections + 4,992 corners).

**Coverage.**  For a vertex `(X, Y, D)` and a square with centre `(CX, CY)/Dk`, `cos = a/r`, `sin = b/r`,
the point is in the closed square iff (with `dx = X Dk - CX D`, `dy = Y Dk - CY D`)

    2 |a dx + b dy| <= r D Dk     and     2 |a dy - b dx| <= r D Dk ,

i.e. `|u| <= 1/2`, `|v| <= 1/2` in the square's own frame, all in integers.  Candidate squares per vertex
come from an exact bounding-box cell prefilter (`floor(coord * 10)`).  `cov(v) = sum over the images
containing v of m_k / 8`; masses have the common denominator `10^9`, so `cov(v) * 8 * 10^9` is an integer
and `M` is the exact maximum.

**Why the maximum of `cov` over `C` is attained at an enumerated vertex.**  `cov` is a finite sum of
indicator functions of *closed* squares, so it is upper semicontinuous and takes finitely many values;
let `x*` attain the maximum and let `S` be the set of support squares containing `x*`.  Every point of
`P = ∩ S` lies in all squares of `S`, so `cov >= cov(x*)` on `P`, and `P` is a nonempty compact convex
polygon (an intersection of closed squares) — in particular it has a vertex `v`.  At a vertex of an
intersection of convex polygons two distinct supporting lines meet, each carrying an edge of some square
of `S` that contains `v` (closed edges).  Either both edges belong to the same square, and `v` is a corner
of that square, or they belong to two squares and `v` is the intersection point of two non-parallel
closed edge segments.  Both kinds are in the enumerated set, so `max_C cov = max_{enumerated} cov`.
Container corners are not needed (every point of `P` is inside a square; if `S` were empty the maximum
would be 0).  **D4 reduction.**  The support is exactly D4-invariant (each pose's mass is spread over its
full orbit), so `cov(g x) = cov(x)` for all `g` in the dihedral group of `C` and the arrangement is
invariant; every vertex is `g` of a vertex in `F`, hence it suffices to evaluate the vertices in `F`.
Both were run: `F` (104,947 vertices) and the whole container without the reduction (836,576 vertices)
give the same `M = 7999999981/8000000000`, attained at
`(889057260620813905630, 896900529698002207551)/898383317495909000000 ≈ (0.989619, 0.998349)`.

**Masses (the only heuristic step, not relied upon).**  On the snapped support the *float* masses of
`dual_PA2_support.txt` (rounded down) have exact `M = 4020298391/4000000000 = 1.00507`, i.e. only
`L = 11.9476` — the LP solution is degenerate on tight vertices, and a `1e-5` perturbation of the
geometry breaks the ties.  So the masses were re-optimised ("polished") once: the exact incidence
of all 126,078 vertices of the 172-pose arrangement in `F` gives 87,956 distinct incidence rows; the LP
`max sum m_k  s.t.  sum_k (count_{v,k}/8) m_k <= 1` (HiGHS, 14 s) gives mass `12.0082308`.  Its solution
was rounded **down** to multiples of `10^-9` (16 poses got mass 0) and re-checked exactly: `M < 1`
already, so no rescaling was needed (the code would otherwise multiply by `1/M` and round down again).
The exact check on the written support file (`check` mode: re-read, re-snap nothing, re-enumerate, no
LP) is what the statement rests on.

## Numbers and timings (4 processes, `nproc = 32` box shared with other jobs)

| step | count | time |
|---|---|---|
| `build`: snap 172 poses, 1376 images, admissibility | — | < 0.1 s |
| vertices in `F` of the 172-pose arrangement | 126,078 | 0.8 s |
| exact incidences | 15,122,067 (vertex, square) pairs | 4.7 s |
| LP (floats, HiGHS) | 87,956 rows x 172 cols, 6.6 M nz | 14.1 s |
| exact `M`, write support | 156 poses, mass `12.008230754` | < 1 s (total `build` 23.8 s) |
| `check` in `F` (from the file) | 104,947 vertices, 11.8 M pairs | 5.0 s |
| `check --full` (whole container, no D4 reduction) | 836,576 vertices, 94.2 M pairs | 38.0 s |

Cross-checks (not part of the proof): the integer intersection routine agrees with a float 2x2 solve on
300 random square pairs; the integer containment test agrees with a float frame test at 3000 random
rational points (no mismatches outside the `1e-7` ambiguity band); the exact coverage is D4-invariant
at 200 random points and their 8 images.

Comparison with `DUAL.md`: float certificate `12.00820` (`M <= 1 + 1.5e-12`, 1,007,508 float
candidates) versus exact `12.00823`; the slightly larger mass is because the LP was re-solved on the
perturbed arrangement with *all* vertex rows at once.

## Reproduce

```sh
python3 search/dual_exact.py build                     # -> runs/dual_exact_3.99_support.txt, .json, .log   (24 s)
python3 search/dual_exact.py check runs/dual_exact_3.99_support.txt          # fundamental domain (5 s)
python3 search/dual_exact.py check runs/dual_exact_3.99_support.txt --full   # whole container   (38 s)
```

`build` is deterministic given HiGHS; `check` depends only on the support file and Python integers.
Support file format: `pose p q cx cy mass` with `cx, cy, mass` as fractions, `t = 399/100`.
