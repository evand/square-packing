# The Helly gap: clique constraints the point-cover LP does not have (2026-08-28/29)

Code: `search/clique_check.py` (max-mass clique of a certified measure), `search/clique_lp.py`
(packing LP over a support file with lazily added clique cuts), `search/plus_pack.py` (heuristic
packer for the plus region).  Nothing in the certificate pipeline changed; `branch.py` gained
`--max-iters`, `--lam-hi` and a per-round dual dump (`runs/branch_TAG_dual.txt`).

## The observation

The cover LP is the fractional **point**-clique cover of the pose overlap graph: the poses through a
point pairwise overlap, so `sum y_S <= 1` over them, and the certificate's total weight bounds the
packing.  Axis-parallel boxes are a Helly family (pairwise intersecting => a common point), rotated
squares are not: three unit squares sitting on the outer sides of a small triangle pairwise overlap
with no common point.  So there are cliques `K` (pairwise interior-overlapping sets of poses) that
are contained in no point-clique, and for each of them `sum_{S in K} y_S <= 1` is a valid constraint
on every fractional packing that the point-cover LP does not impose.  In cover terms: a certificate
may consist of weighted *cliques* with every pose covered by weight `>= 1`, points being the special
case; the reduction to `n <= W` is the same one line (a packing has at most one square in a clique).

## Measured on the certified extremal measures

`clique_check.py` expands the D4 symmetry (mass `mu/8` per image), builds the interior-overlap graph
(separating-axis test, strict) and finds the maximum-mass clique exactly (branch and bound with a
greedy-colouring weight bound).  Input: the certified measures of `DUAL.md` (coverage `<= 1` at every
arrangement vertex).

| t | measure | certified mass | max clique mass | size | common point |
|---|---|---|---|---|---|
| 3.97 | `PD2` | 11.807 | **1.320** | 278 | none |
| 3.98 | `PD1` | 11.918 | **1.326** | 267 | none |
| 3.99 | `PA2` | 12.008 | **1.327** | 221 | none |
| 4.00 (closed) | `PC1` | 12.163 | **1.350** | 123 | none |

Independent check of the `PA2` clique: all 24,310 pairs overlap by polygon clipping (min intersection
area `1.7e-7`); no point of a 0.001 grid lies in more than 167 of the 221 squares.  Its anatomy: 163
members through `(0.99, 1.275)` with mass exactly 1.000 (the certified tight point: the wall square
`[0,1] x [1,2]` and its perturbations) plus 58 interior tilted poses of mass 0.327 (centres near
`(1.5-1.7, 1.5-1.7)` at 45 deg and `(0.65-0.7, 1.8)` at 22-37 deg) that overlap every one of them.  It
is *not* a "majority of a disc" clique (best disc: min fraction 0.30 < 1/2); it is the container-
dependent kind: near a wall every square through `p` is essentially the wall square, so interior
squares that overlap the wall square overlap all of them.  The heaviest non-Helly *triple* is only
0.27-0.89, so the violation is a property of the whole cluster, not of three poses.

So the extremal measures violate clique constraints by `~0.33` at every `t`, against a pure-method
excess of 0.008 at 3.99.

## How much the LP moves (fixed pose set)

`clique_lp.py` re-solves the packing LP over the support's poses (point constraints on a 0.01 grid of
the fundamental domain, so its round-0 value is a little above the certified one) and adds clique cuts
lazily: the global max-mass clique plus the best clique through each of the 40 heaviest poses per
round.  Mass migrates to neighbouring cliques — there is a continuum of them — so the descent is slow,
but it is monotone and it crosses 12:

| t | measure | LP over grid, no cuts | after cuts | cuts | max clique at the end |
|---|---|---|---|---|---|
| 3.98 | `PD1` | 11.952 | **11.816** | 2881 | 1.22 |
| 3.99 | `PA2` | 12.053 | **11.962** | 2409 | 1.28 |
| 4.00 (closed) | `PC1` | 12.062 | **11.754** | 2776 | 1.18 |

Caveats, both directions: the pose set is fixed (column generation over poses would raise the value
back somewhat — the extremal poses of the *point*-cover LP are not those of the clique LP), and the
final max clique is still `> 1` (more cuts would lower it further).  What is established is that on
the pure method's own extremal poses the clique-strengthened relaxation is below 12 at 3.99 and at
the limit 4.0, by 0.05 and 0.25 respectively (0.14 below the point-cover value at 3.98), with the
loops not yet converged.

## Reading

* The pure ceiling `s* in [3.9686, 3.99)` is the ceiling of the *Helly sub-relaxation*.  The clique
  relaxation's ceiling is higher; where it sits needs column generation on both sides (poses and
  cliques) — the pricing problem on the clique side is exactly `max_weight_clique`, which is seconds
  on ~1,500 poses.
* For a certificate, cliques must be verifiable pose *regions*.  Point cliques are `{S : p in S}`.
  The cliques found here have the form "all poses through `p` (near a wall) together with poses that
  overlap all of them"; a verifiable version is `{S : p in S} ∪ {S : S ⊇ A}` for a set `A`, valid
  when every admissible square through `p` meets `int A` — a non-avoidance lemma of the Bentz kind,
  checkable by the existing sweep.  Whether such structured cliques capture the 0.05-0.25 above is the
  next question; the finite-graph experiment says the mass is there.
* Framing: corner branching is a hand-picked lift; clique inequalities are the rank-1 strengthening of
  the same relaxation.  Both attack the same fractional mass (walls and interior), which is why the
  corner leaves an integral packing resembles gained nothing.

## Side results

**Plus region** `P_t = [0,t]^2` minus its four corner unit squares (`plus_pack.py`, penalty
minimisation with random restarts, bisection on `t`; heuristic, one-sided):

| squares | fits at | not found at |
|---|---|---|
| 8 | 4.00020 | 3.99990 (and 3.98 with 300 restarts: residual `1.3e-4`) |
| 7 | 3.92500 | 3.92422 |

So the only way to put 8 in the plus region is the grid at 4, and 7 already fit at 3.925: two
axis-parallel wall squares plus a tilted row of five at ~54 deg through the centre (Stromquist's
tilted-row trick).

**Correction to `BRANCH.md` ("the all-corners leaf is exactly 12").**  The `k = 4` loops never
produced a valid cover: their probe minimum was 0.27-0.94 in every one of 31 rounds (final rounds
0.0, 0.4, -0.0, 0.5 — the verifier's value is `covered - lambda*[box]`, pass is `>= 1`), and the
value pinned at `12(1+margin)` is the LP cycling on a pruned row set with a fixed column set.  The
80-point dyadic "cover" fails coverage by 0.75 at `N = 2000`, at 3.98, 3.995 and 4.0 alike.  Its
dual was never certified either (pricing at 1-3 %), so nothing is known about the leaf's value in
either direction; the statement "corner branching cannot beat the pure ceiling" rested on this.
A clean run (`t398hk4`: seeded from the `k = 2` certificate so that the warm start has 10,347 box
rows, no pruning, 60 rounds of column generation) climbed 11.889, 11.930, 11.940, 11.964, 11.971,
11.979, 11.979, 11.980, 11.984, 11.991, 11.993, 11.996, 11.996, 11.998, 11.999, 11.998 over rounds
0-15 with probe minima 0.73-0.98 (never a cover), and at round 16 — **without pruning** — jumped to the
same degenerate vertex, `W = 18.000024, lambda = 1.5`, value `12.000024`, 16 orbits.  So that vertex
is a genuine optimum of the row-sampled LP: the rows present admit a fractional packing of mass
exactly 12 over the current 6,095 columns (pricing gain 2 %, so uncertified).  Rounds 17-18 then
alternated 11.9995 (column step, probe 0.97) and 12.000024 (row step, the degenerate vertex again), at
which point the run was stopped: the leaf value at 3.98 is `12.000 +- 0.001`, straddled by the row and
column steps of the cutting-plane loop, and the corner level does not close it.  Its dual
(`runs/branch_t398hk4_dual_it16.txt`): corner mass 1.000 per box (0.845 on the axis-parallel corner
square, 0.15 on 8.5-deg tilted corner squares), 3.75 on the eight wall slots (axis-parallel wall
squares at 0.24 + 0.08 per slot and neighbours), 4.0 in the interior 1.0-1.6 from the walls, mostly
tilted 5-30 deg.  Its max clique has mass **1.49** (385 poses, no common point) — larger than the
1.33 of the pure measures.  Reading: the `k = 4` leaf at 3.98 sits at `12.00 +- 0.005` — the substance
of the withdrawn claim survives, its justification did not — and what keeps it there is a
wall-plus-interior cluster that clique constraints cut by 0.5.  (The `clique_lp.py --branch --kmass 4`
loop on this dual's 520 poses did *not* move from 12.000000 in 40+ rounds / 1,640 cuts: with or without
the corner constraint the 0.01-grid LP over these poses is exactly 12 — a finite instance tight at 12,
and too coarse here, its solution having coverage 1.06 on a 0.002 grid.  The pose set of a branch dual
is not a substitute for column generation.)  The multiplier cap is irrelevant to the
value: the same run with the cap raised from 1.5 to 4 (`t398ik4`) gave the identical round-0
objective 11.888517 with `W` higher by exactly `4 x 2.5` — the zero-cost direction (a four-point
diagonal corner orbit up, `lambda` up) — and was stopped.  If the leaf closes at 3.98, together with `k = 0, 1, 2` (certified) and
`1110` (11.97, rising) that is `s(12) >= 3.98`.

## Reproduce

```sh
python3 search/clique_check.py runs/dual_PA2_support.txt            # max clique, ~30 s
python3 search/clique_lp.py runs/dual_PA2_support.txt --rounds 60 --per-round 40
python3 search/plus_pack.py --m 8 --bisect 3.9 4.05 --restarts 120 --steps 9
python3 search/branch.py certificates/branch/s12_t3.98_corner_k2.txt t398hk4 --k 4 --r 1 --N 2000 \
    --colgen 60 --cg-want 300 --cols runs/branch_t398fk4_cols.txt --warm-thr 3 --prune-at 400000 --row-grid 0.002 0.005
```
