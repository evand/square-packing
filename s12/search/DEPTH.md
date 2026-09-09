# Depth test: branching the interior of the hardest `t = 4` leaf on its four quadrants

Task: `tasks/depth-test/README.md`.  Code: `search/t4leaf.py --quad` (a flag; off, the script is
bit for bit the one of `search/T4LEAF.md`, §6), `search/leaf_ceiling.py --quad` (still
self-contained).  Runs live in this worktree's `runs/tl_D*`, `runs/lc_D*` (gitignored, so every
number is quoted here).  Read against `search/T4LEAF.md` (the parent leaf and its instrument),
`search/T4SCREEN.md` §1 and §4.2 (semantics; `leaf value = 8 + interior`) and
`notes/branch-semantics.md` §3 (the max rule and its packing-side dual).

---

## 0. Verdict, up front

**Drop per level = 0, exactly — and it is a theorem, not a measurement.**  Branching the
interior `[1,3]^2` of leaf `01010101` on the four quadrant counts `q = (q_0, q_1, q_2, q_3)`,
`sum q = 4`, has a child — the balanced one, `1111` — whose cover-side bound is *identical* to
the parent's: the slot pattern is `C_4`-invariant, `C_4` permutes the quadrants transitively, and
the `C_4`-average of any feasible parent measure with interior mass `m` is a feasible measure of
the same mass with exactly `m/4` in each quadrant (§1.4).  So `max over children = parent`, the
level buys nothing at its hardest child, and by the brief's own rule (`<= 0.1` ⇒ new
mathematics) **this family needs a different branch or new mathematics, not one more level of
this one** (§5 says which branch).

The computation agrees, and adds what the argument does not give.  All at `t = 4`, closed
semantics, `r = 1`, corners `1111`, slots `01010101`, anchor cliques + chord rows + boundary
duplication on, cap relaxation `mu(Q_i) <= q_i` (§1.3), warm-started from `T4LEAF.md`'s
`A01010101` checkpoints plus the `C_4` images of their support:

| child `q` (Q0 Q1 Q2 Q3) | type | LP (float) | `M` | `kmax` | quadrant masses | exact certificate (`leaf_ceiling.py`) | drop vs parent `11.745` |
|---|---|---|---|---|---|---|---|
| **`1111`** (`D1111` st0 / st1) | balanced = parent | **`11.804463`** / `11.823552` | `1.000000000` / `1.0086` | `1.000000` / `1.0253` | `[0.959 0.942 0.960 0.943]`, every cap slack, cap duals `0` | coverage + regions exact at **`11.804462741`**; **not clique-feasible**: the exact scan finds an anchor clique of mass `1.165738` the separator missed (§4) | **`0`** (`-0.059`: the child is *above* the parent's recorded value, §2) |
| `0211`, `0121` (mirror pair) | 2 next to the empty quadrant | `11.000000` | `1.000000000` | `1.000000` | `[0 1 1 1]` — the cap `2` slack at `1` | **exact `11`, all conditions proved** (11 pairwise disjoint unit squares) | `0.745` |
| `0112` | 2 opposite the empty quadrant | `11.000000` | `1.000000000` | `1.000000` | `[0 1 1 1]` | exact `11`, all proved | `0.745` |
| `0022` | two 2's in one half | `10.000000` | `1.000000000` | `1.000000` | `[0 0 1 1]` | exact `10`, all proved | `1.745` |
| `0220` | two 2's diagonal | `10.000000` | `1.000000000` | `1.000000` | `[0 1 1 0]` | exact `10` (integral witness; the LP's own fractional optimum has an anchor clique of mass `3/2` the separator missed, §4) | `1.745` |
| `0013`, `0031` (mirror pair) | 3 next to the 1 | `10.000000` | `1.000000000` | `1.000000` | `[0 0 1 1]` — the cap `3` slack at `1` | exact `10`, all proved | `1.745` |
| `0130` | 3 opposite the 1 | `10.000000` | `1.000000000` | `1.000000` | `[0 1 1 0]` | exact `10` (integral witness) | `1.745` |
| equality child `1111` (`D1111EQ`, the brief's semantics) | — | **infeasible** at the first solve (9287 columns, 183 cliques) | — | — | — | — | (§1.2: `12` or nothing) |

Three things the table says:

* **The unbalanced children are trivial, and integral.**  Every one converges at its first stage
  to exactly `8 + (number of non-empty quadrants)`, with *one* unit square per non-empty quadrant
  and the caps of `2` and `3` slack at `1`; their optima are packings of 10 or 11 unit squares
  (exactly certified, pairwise disjoint, clique bound proved).  Emptying a quadrant does not let
  a neighbour hold more than one square: the frame leaves each quadrant room for one, and the
  parent's fractional excess `0.745` lives entirely in the **crowding at the centre `(2,2)`**,
  which only the balanced child keeps.  These children would close on the cover side with room
  to spare — but they are not where the difficulty is.
* **The balanced child is the parent, numerically as well as in principle.**  Its optimum has
  all four caps slack and all four cap duals exactly `0`, so by the dual argument of
  `T4LEAF.md` §1.3 the same measure is optimal for the parent LP on the same pose set: the child
  value *is* the parent value.  It came out `0.059` above `T4LEAF.md`'s `11.745193` because the
  warm start added the `C_4` images of the support (which the symmetric average needs) and the
  separator reported its cliques converged (`kmax = 1.000000` where `A` stage 2 still had
  `1.000859`).  The exact scan then shows the separator wrong by `0.17` (§4), so the honest
  reading is: balanced child `= parent` on every pose set (cap duals `0`), drop `0`; and the
  number `11.80` itself, like `T4LEAF.md`'s `11.745`, is a coverage-exact value that the full
  anchor-clique family does **not** admit — the clique-constrained leaf value is lower.
* **The brief's semantics is degenerate here (§1.2).**  With the frame at 8, pinning quadrants
  summing to 4 pins the total at 12: each equality child is infeasible or exactly `12`, and on
  the warm start every one is infeasible (`D1111EQ`) — which says only what `11.745 < 12` said.

---

## 1. What is measured

### 1.1 The LP

`T4LEAF.md` §1's LP for leaf `01010101` under corner leaf `1111` — closed unit squares in the
closed `[0,4]^2`, coverage `<= 1` at every point, four corner boxes pinned at 1, the eight slots
at `01010101`, anchor cliques `<= 1`, chord rows `<= 3`, the verifier-compatible boundary
duplication of `T4LEAF.md` §1.1 with `delta = 1e-6` — plus four **quadrant rows**

```
mu({S : centre(S) in Q_i})  (=  or  <=)  q_i ,       i = 0..3,
Q_0 = [1,2]x[1,2]   Q_1 = [2,3]x[1,2]   Q_2 = [1,2]x[2,3]   Q_3 = [2,3]x[2,3]
```

(bottom-left, bottom-right, top-left, top-right; column labels `12..15`).  Boundary semantics is
the verifier-compatible one: a pose within `delta` of `c_x = 2` or `c_y = 2` gets one column per
adjacent quadrant (identical coverage and clique coefficients, different label), a pose at
`(2, 2)` gets four, and a pose on the interior/slot boundary `c = 1` or `3` gets, as before, one
column per adjacent region on either side — a centre at `(1, 2)` has the four labels
`W_7, W_6, Q_0, Q_2`.  This is the packing-side dual of `notes/branch-semantics.md` §3.1's max
rule with the quadrants added to the region list, and `packing_le_weight_regions_choice` is the
theorem a child certificate would discharge.  Adding columns can only raise the value, so every
number here stays a lower bound on what a child certificate must beat.  (`x = 2`, `y = 2` were
already in `t4leaf.py`'s duplication list as the mid-wall slot boundaries, so the only new
duplicates are the interior poses on them: `319` boundary duplicates on the warm start against
`~210` without the quadrants.)

The pricer prices the sixteen labels separately (`t4screen.price` stratifies by region), the
reduced-cost check and the pose sifting see the quadrant multipliers through `lam_map`, the
matched variants (`pure`, `nochord`) carry the quadrant rows, and the per-stage record gains
`quads` and `quad_dual`.  With `--quad` absent nothing changes: `--threads 1`, the pre-change
script and this one produce identical logs and identical `hist` records on a 1 594-pose,
160-clique instance with cliques, chord rows, variants and boundary duplication on (§6).

### 1.2 The equality semantics is degenerate in an `m = 4` leaf

`T4SCREEN.md` §4.2: with the frame pinned, `leaf value = 8 + interior mass` exactly.  Pinning the
quadrants to `q` with `sum q = 4` pins the interior at 4 and hence the **total at 12** for every
feasible measure: the LP is a feasibility problem, its value is `12` when feasible and it has no
value when not, and by LP duality the cover side is then `12` or `-inf`.  There is no `11.7`.  In
particular

> the equality child `q` is feasible  **iff**  the leaf's true value is `>= 12` with a measure
> whose quadrant masses are exactly `q`;

and the parent's `11.745 < 12` on its pose set says the pinned-interior child (`q` free,
`sum = 4`) is infeasible there, hence every quadrant child is.  That is what `D1111EQ` reports —
HiGHS returns `kInfeasible` on the first solve, before any row or column generation — and it
measures no pruning power.  (The parent LP escapes this only because it leaves the interior
count free, which is what gives it a continuous value.)

### 1.3 The cap relaxation, and what its value means

Replace the equalities by `mu(Q_i) <= q_i` and maximise the mass as before (`--quad-mode le`).
Then

* `value(q) <= parent value` (the parent LP on the same pose set, minus constraints);
* `value(q) >= 12  iff  the equality child q is feasible` (if `>= 12` then the interior carries
  `>= 4` under caps summing to 4, so every cap is tight; conversely a feasible equality measure
  has mass 12 and satisfies the caps) — **the cap child closes iff the branch child closes**;
* cover side, it is the branch certificate restricted to `lambda_{Q_i} >= 0` — the quadrant
  multipliers may only *charge* an interior pose extra capture, never credit it — and the
  parent's certificate is the special case `lambda_Q = 0`.

`drop = parent - max_q value(q)` is then the honest form of the brief's number: it is `0` iff the
parent's optimum can be re-labelled as a child, and it is what one more level of the branch buys
at its hardest child.

### 1.4 The balanced child equals the parent — the symmetry argument

Let `rho` be the quarter turn about `(2,2)`.  The container, the corner boxes, the coverage and
clique constraints and the chord rows are `rho`-invariant; so is the slot pattern `01010101`
(`W_1 -> W_3 -> W_5 -> W_7 -> W_1`); and `rho` permutes the quadrants in the cycle
`Q_0 -> Q_1 -> Q_3 -> Q_2 -> Q_0`.  If `mu` is feasible for the parent with interior mass `m` and
quadrant masses `(a, b, c, d)`, then `mu' = (mu + rho mu + rho^2 mu + rho^3 mu) / 4` is feasible
for the parent (the constraints are linear and `rho`-invariant), has the same total mass, and has
quadrant masses `(m/4, m/4, m/4, m/4)`.  Hence

* **cap child `1111`:** `mu'` satisfies `mu(Q_i) <= 1` whenever `m <= 4`, so
  `value(1111) >= parent`, hence `value(1111) = parent` — on the continuum, and on every
  `C_4`-symmetric pose set;
* **equality child `1111`:** feasible iff the parent has a feasible measure of interior mass
  exactly 4, i.e. iff the parent's true value is `>= 12`.

Either way the balanced child is the parent under another name, and it is the hardest child
(every other `q` is a proper sub-cap of some rotation of the balanced situation, and §2 shows
the others are far below).  **The drop is zero before any computation.**  The runs only had to
confirm that the instrument agrees — it does, to the last digit: the child's optimum has all
four cap duals `0` — and to measure how cheap the unbalanced children are.

The same argument kills any `C_4`-invariant branch whose balanced child is the symmetric average
— a `k x k` interior grid with `k` even, "how many squares in each half", the quadrant counts of
any `C_4`-invariant region — under any `C_4`-invariant slot pattern (`01010101` and `10101010`
at `m = 4`).  It does not apply to the other two `m = 4` patterns (`D_2`-invariant only), nor to a
partition the symmetric average cannot land on (§5).

### 1.5 The bound `q_i <= 3` (and why not 2)

Two squares of a packing have side `L > 1` after rescaling to the certificate's container, and
their centres are at distance `>= L > 1`: the support function of a square of side `L` is
`>= L/2` in every direction, and interior-disjoint convex bodies are separated by a line, so the
centre difference projects to `>= L/2 + L/2` on its normal.  So the `q_i` centres in the closed
unit box `Q_i` are pairwise `> 1` apart, and

> **a closed unit square contains no four points with pairwise distances `> 1`.**

*Proof.*  Halve the square both ways: each quarter (diameter `sqrt 2 / 2 < 1`) holds at most one
point, so there is exactly one per quarter — `A` bottom-left, `B` bottom-right, `C` top-left,
`D` top-right.  `A, B` lie in the bottom strip of height `1/2`, so `|AB| > 1` forces
`B_x - A_x > sqrt(1 - 1/4)`, i.e. `A_x < delta_1 := 1 - sqrt(3)/2` and `B_x > 1 - delta_1`; the
same for the other three adjacent pairs puts each point within `delta_1` of its corner in both
coordinates.  Iterate: with both points of an adjacent pair within `delta_n` of the square's
edge in the transverse coordinate, `|AB| > 1` forces the longitudinal gap
`> sqrt(1 - delta_n^2)`, i.e. `delta_{n+1} = 1 - sqrt(1 - delta_n^2) < delta_n^2`.  So
`delta_n -> 0`, all four points sit at the four corners, and `|AB| = 1`, not `> 1`.  ∎

So `q_i <= 3`.  The brief's `q_i <= 2` ("diameter `sqrt 2`, so at most two centres more than 1
apart") is false for *points* — `(0,0), (1, 0.268), (0.268, 1)` are pairwise `1.035` apart — but
for *squares* the separation must be along an edge normal, which is stronger, and
`level2_capacity.py --t 4 --box 1 2 1 2 --k 3` finds nothing better than margin `+0.000000`
(three grid squares touching pairwise at corners, `(1,1), (2,1), (2,2)`), as for `k = 4`.  That is
one-sided evidence for `q_i <= 2`, not a proof; the enumeration uses the proved `q_i <= 3`, and
the table shows the `3`-children are worth exactly what the corresponding `2`-children are (the
third unit of cap is never used).

### 1.6 The children

`sum q = 4`, `0 <= q_i <= 3`: 31 count vectors.  The symmetry of the pattern is `C_4` (the
reflections map `01010101` to `10101010`, so they are not symmetries of the leaf); orbits under
the quadrant cycle `Q_0 -> Q_1 -> Q_3 -> Q_2`, written `(q_0 q_1 q_2 q_3)`:

```
q_i <= 3 : 9 children      (q_i <= 2 : the 6 without a 3)
  1111                       balanced                                        orbit 1
  0211   0121                2 next to the empty quadrant, a mirror pair      orbit 4 each
  0112                       2 opposite the empty quadrant                    orbit 4
  0022                       two 2's in one half                              orbit 4
  0220                       two 2's diagonal                                 orbit 2
  0013   0031                3 next to the 1, a mirror pair                   orbit 4 each
  0130                       3 opposite the 1                                 orbit 4
```

(Burnside: `(31 + 1 + 3 + 1) / 4 = 9`; with `q_i <= 2`, `(19 + 1 + 3 + 1) / 4 = 6`.)  The mirror
pairs are distinct children of *this* leaf; each is the mirror image of the other under the
mirror leaf `10101010`, so a tree over both leaves would run one of each pair.  All nine were
run.

---

## 2. The trajectory

Nine processes, all `t = 4.0`, corners `1111`, slots `01010101`, cliques (wall + interior, cap
700, age 5, want 200), chord rows, `bnd-delta 1e-6`, cap relaxation, `--variants
--variant-every 3`, `--pose-max 8000`, `--row-cap 3000`, `--rowloops 25`, HiGHS backend with
`--lp-tlim 20`.  Warm start: `T4LEAF.md`'s `tl_A01010101_{poses,dual,cliques}.txt` (8 050
columns, 315 binding rows, 162 cliques) plus the `C_4` images of that checkpoint's 394 support
poses (`--load-poses-c4`; 9 287 columns for 8 968 poses after duplication).  `D1111` ran with 2
threads and a long budget; the eight unbalanced children ran one after another with 1 thread
and 3 stages each (`runs/launchD.sh`, `runs/launchDseq.sh`).  A "stage" is one pricing round
after the row/clique loop and the settle solve, exactly as in `T4LEAF.md` §2; `conv` = `M <= 1 +
1e-9` and `kmax <= 1 + 1e-6` at the settle solve.

```
tag       q     st         LP           M     kmax      int                          quads  cols   rows   cq       pure    nochord    bnd conv   secs
D1111     1111   0  11.804463 1.000000000 1.000000 3.804463  [0.9594 0.9421 0.9597 0.9432]  9287  20335  197  11.932869  11.804463 0.3831 YES   3620
D1111     1111   1  11.823552 1.008595168 1.025307 3.823552  [0.9551 0.9568 0.9625 0.9491]  8077  27811  205          -  11.823552 0.4472 no   17160
D0211     0211   0  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  9287  10395  501  11.000000  11.000000 2.0000 YES    551
D0211     0211   1  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  8178  11283  405          -  11.000000 1.0000 YES    957
D0211     0211   2  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  8077  11283   40          -  11.000000 2.0000 YES    983
D0121     0121   0  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  9287  10775  283  11.000000  11.000000 1.0000 YES    352
D0121     0121   2  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  8066  11238  687          -  11.000000 2.0000 YES    618
D0112     0112   0  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  9287  10906  487  11.000000  11.000000 1.0000 YES    673
D0112     0112   2  11.000000 1.000000000 1.000000 3.000000  [0.0000 1.0000 1.0000 1.0000]  8068  11480  208          -  11.000000 1.0000 YES   1022
D0022     0022   0  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  9287   7245  202  10.000000  10.000000 1.0000 YES     29
D0022     0022   2  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  8100   8411  546          -  10.000000 2.0000 YES    373
D0220     0220   0  10.000000 1.000000000 1.000000 2.000000  [0.0000 1.0000 1.0000 0.0000]  9287   7385  200  10.000000  10.000000 2.5000 YES     22
D0220     0220   2  10.000000 1.000000000 1.000000 2.000000  [0.0000 1.0000 1.0000 0.0000]  8088   7614  392          -  10.000000 1.3333 YES     78
D0013     0013   0  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  9287   9113  598  10.000000  10.000000 1.0000 YES     47
D0013     0013   2  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  8085   9164  146          -  10.000000 1.0000 YES    154
D0031     0031   0  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  9287   7693  392  10.000000  10.000000 1.0000 YES    119
D0031     0031   2  10.000000 1.000000000 1.000000 2.000000  [0.0000 0.0000 1.0000 1.0000]  8102   8352   91          -  10.000000 1.0000 YES    508
D0130     0130   0  10.000000 1.000000000 1.000000 2.000000  [0.0000 1.0000 1.0000 0.0000]  9287   7481  628  10.000000  10.000000 1.0000 YES     55
D0130     0130   2  10.000000 1.000000000 1.000000 2.000000  [0.0000 1.0000 1.0000 0.0000]  8068   8031  370  10.000000  10.000000 1.0000 YES    406
```

(`search/depth_table.py`; the stage-1 rows of the unbalanced children are identical to their
stage-2 rows and are omitted.)  `int` is the interior mass, `quads` the mass per quadrant by
column label, `bnd` the mass within `1e-7` of a region boundary, `pure` / `nochord` the matched
variants of `T4LEAF.md` §1.3; `LP = 8 + int` holds to twelve digits in every row, with corners
`[1,1,1,1]` and slots `[0,1,0,1,0,1,0,1]` exactly.

**The balanced child.**  `D1111`'s stage 0 took the warm start from a row-starved `12.000000`
(`M = 1.28`) through 25 row/clique loops to `11.804463` with `M = 1.000000000` and `kmax =
1.000000` at the settle solve (one solve, iteration 12, fell into HiGHS's ipm-crossover
pathology of `T4LEAF.md` §1.3 and took 1 201 s; the rest were `85–115` s).  Its quadrant masses
`[0.959, 0.942, 0.960, 0.943]` are within `0.02` of the symmetric `3.804/4 = 0.951`, all four caps
are slack and all four cap duals are exactly `0.0`, so this measure is optimal for the parent LP
on the same rows and columns (`T4LEAF.md` §1.3's dual argument): **the balanced child's value is
the parent's value on this pose set.**  It sits `0.059` above `A01010101` stage 2's `11.745193`
because the pose set now contains the `C_4` images of the support and because the clique
separator converged (`A` stopped with one clique violated by `8.6e-4`).  Its matched no-clique
value is `11.932869` (`M_pure = 1.054`, not a bound, as always for the `pure` column).
Pricing after stage 0 found improving interior columns (best interior reduced cost `+0.337`),
and stage 1 reached `11.823552` (`+0.019`) with `M = 1.0086`, `kmax = 1.025` when its 25 loops
ran out — not converged, caps still slack (`[0.955 0.957 0.963 0.949]`).  Stage 1 took 3.8 h:
on the 8 077 x 27 811 model six solves fell into the ipm-crossover pathology (`1 134–3 439` s
each) against `85–115` s for the rest; the run was stopped in stage 2 at its 5-hour budget.  So
the parent keeps rising under pricing, `11.745 -> 11.804 -> 11.824`, exactly as `T4LEAF.md` §2
found — and the balanced child rises with it, because it is it.

**The unbalanced children** all converge at stage 0, in `22–673` s, to exactly
`8 + (number of non-empty quadrants)`, and two further pricing stages change nothing (interior
reduced cost `+0.000000`; the reported gaps of `+2` / `+3` are the frame-region artefact
`t4screen.price` documents).  The caps of `2` and `3` are never used: a quadrant next to an
empty one still carries exactly `1`.  The `pure` value equals the clique value in every case —
the anchor cliques do nothing once the centre is uncrowded.

---

## 3. Anatomy

### 3.1 The hardest child: the parent's "grid smeared", with the centre crowded

`runs/tl_D1111s0_measure.txt` (stage 0, 495 support columns, `11.804463`;
`t4leaf_table.py --anatomy --cliques`):

```
corners   4.000000   per box  1.00000 1.00000 1.00000 1.00000
slots     4.000000   per slot 0.00000 1.00000 0.00515 0.99485 0.11202 0.88798 0.08468 0.91532
interior  3.804463   per quadrant 0.9594 0.9421 0.9597 0.9432
wall strips 3 3 3 3;  mass on a region boundary (1e-7): 0.383106 on 35 poses
axis-parallel mass: 6.767503   tilted: 5.036960
      mu= 0.842569  ( 3.50000,  3.50000)     0.00 deg   C3          mu= 0.433904  ( 2.48771,  2.50786)    -0.96 deg   Q3 (19 poses)
      mu= 0.842118  ( 0.50000,  3.50000)     0.00 deg   C2          mu= 0.421605  ( 1.49038,  2.48653)    -0.96 deg   Q2 (18 poses)
      mu= 0.839792  ( 0.50007,  0.50007)     0.00 deg   C0          mu= 0.407214  ( 1.51437,  1.48892)    -0.96 deg   Q0 (20 poses)
      mu= 0.833777  ( 3.50000,  0.50000)     0.00 deg   C1          mu= 0.379437  ( 2.50776,  1.51177)    -0.96 deg   Q1 (20 poses)
      mu= 0.612239  ( 3.49904,  2.49484)     0.00 deg   W3          mu= 0.146189  ( 2.40455,  2.56316)     9.80 deg   Q3 (8 poses)
      mu= 0.579193  ( 0.50183,  1.50106)     0.00 deg   W7          mu= 0.141892  ( 2.56036,  1.58960)    10.00 deg   Q1 (10 poses)
      mu= 0.577781  ( 1.50485,  3.49932)     0.00 deg   W5          mu= 0.141378  ( 0.50719,  1.99873)     0.00 deg   W7 (8 poses)
      mu= 0.562583  ( 2.49874,  0.50072)     0.00 deg   W1          mu= 0.140134  ( 2.00107,  0.51077)     0.00 deg   W1 (9 poses)
```

(clusters of radius `0.06`, 6 degrees; 138 clusters in all.)  It is `T4LEAF.md` §4's object with
one slot per wall switched off: `0.83–0.84` on the four corner grid poses, `0.56–0.61` on the
four *occupied* wall grid poses, `0.38–0.43` on each of the four central grid poses — now all
tilted by `-0.96` degrees, a common small rotation that takes the four central squares off the
exact crossing `(2,2)` — and `0.14` on `10`-degree satellites of the two right-hand central
poses; `5.04` of the `11.80` is tilted, `0.38` sits on region boundaries (`0.14` at each of the
two mid-wall points `(0.5, 2)` and `(2, 0.5)`, booked to the occupied slot).  **The mass is grid
mass, and the crowding at the centre is where the fractional excess lives**: the four central
poses carry `1.64` between the four quadrants where a packing could place at most one square at
`(2,2)`-adjacent grid positions without them touching.  That is exactly what the unbalanced
children remove and what the balanced child keeps.

### 3.2 The unbalanced children: packings

`D0211`'s optimum (`runs/lc_D0211.txt`, exact):

```
  mu=1  (0.8087, 0.5000)   0.00 deg  C0       mu=1  (2.0000, 0.5000)   0.00 deg  W1 (on the mid-wall line, booked W1)
  mu=1  (3.4997, 0.5003)  -0.04 deg  C1       mu=1  (3.4975, 2.0000)   0.20 deg  W3 (on the mid-wall line, booked W3)
  mu=1  (0.5000, 3.5000)   0.00 deg  C2       mu=1  (1.5000, 3.5000)   0.00 deg  W5
  mu=1  (3.5000, 3.5000)   0.00 deg  C3       mu=1  (0.5086, 1.5120)  -1.00 deg  W7
  mu=1  (2.4592, 1.5000)   0.00 deg  Q1       mu=1  (2.5250, 2.7181) -39.10 deg  Q3       mu=1  (1.4628, 2.2786)  46.00 deg  Q2
```

Eleven unit squares, pairwise disjoint (the exact arrangement has **no** edge crossings), mass
1 each: a packing, not a fractional measure.  `D0112` and `D0121` are the same picture with the
tilted pair placed differently (`45` / `41` degrees); `D0022`, `D0013`, `D0031` are ten-square
packings with the two interior squares in the top half at `60` / `-37` degrees; `D0220` and
`D0130` return fractional vertices of value exactly `10` (masses `1/2` and `1/6`), for which the
ten-square packing `lc_D0211.txt` minus its `Q3` square is an exactly certified integral witness
of the same value.

---

## 4. Which numbers here are what

* **Exact.**  The eight unbalanced children: `leaf_ceiling.py snap --polish` (angle lattice
  `1e-8`, centres to `1e-8`, masses to `1e-9`) followed by `check --anchor clique` certifies, in
  integer/rational arithmetic, mass **`11`** (`0211`, `0121`, `0112`) and **`10`** (`0022`,
  `0013`, `0031`; `0220` and `0130` through the integral witness) with coverage `<= 1` at every
  arrangement vertex, the twelve frame equalities and the four caps exact under the choice
  semantics (each pose counted in the region the LP booked it to — `leaf_ceiling.py` now reads
  and honours those bookings, §6), the chord inequalities, and the anchor-clique family proved
  (the maximum-weight closed-intersection clique is `1`, i.e. the squares are pairwise
  disjoint).  These are theorem-grade lower bounds on the children's values, and they are
  attained by packings, so they are also upper bounds on nothing.
* **Exact coverage, cliques as the instrument says.**  The balanced child's stage-0 measure:
  rounded down and topped up without re-optimisation, mass **`11.804462741`** with exact
  `M = 1` and the regions and caps exact (`runs/lc_D1111s0_np.txt`); re-optimised on the complete
  exact arrangement rows of its support *without* the clique rows, mass **`11.872301651`**,
  exact `M = 1` (`runs/lc_D1111s0.txt`).  Coverage and regions are therefore certified for both.
  Anchor cliques are **violated**: the run's separator (wall pitch `0.02`, interior pitch
  `0.04`, 12 directions) reports `kmax = 1.000000` on the first, `t4leaf_table.py`'s finer scan
  (`0.01` / `0.02`, 24 directions) finds `1.024253`, and the exact scan of `leaf_ceiling.py`
  (`check --anchor scan`: float candidates at 200 arrangement vertices, the 60 best evaluated
  exactly, then the exact local ascent over the candidate set of `LEAF_CEILING.md` §2.3) finds
  **`mu(K(p, A)) = 1.165738071`** on the run's own measure — `p ~ (2.000001, 2.026658)`, `A` a
  short segment beside it, `|Cont| = 43`, `|Meet| = 90`, 86 members re-checked pairwise
  intersecting — and `1.211105988` (`p ~ (2.256, 2.033)`, 59 members) on the polished one.
  Both are certified lower bounds on the true maximum over the family (`scan` never proves
  `<= 1`; here it does not need to).  So `11.804463` is a value the poses attain for
  coverage and regions, and **not** a measure the anchor-clique family admits: the LP was short
  of the cliques that bite, and the clique-constrained value of the leaf on this pose set is
  strictly lower.  `T4LEAF.md`'s `11.745193` (`kmax = 1.000859` by the same separator) has the
  same status and had never been exactly scanned; its clique side is now known to be unreliable
  at the `0.1` level.  The clique-free `11.872302` is a genuine exact lower bound on the leaf's
  *no-clique* value (below `T4LEAF.md`'s float `11.920067`, because it lives on the clique run's
  support).
* **The separator is not a certificate — measured.**  On `D0220`'s 24-square fractional
  optimum the complete anchor enumeration of `LEAF_CEILING.md` §2.3 finishes and finds
  `mu(K(p, A)) = 3/2` — three half-mass squares `(3.054, 0.5)`, `(3.491, 0.509)`,
  `(2.762, 1.464, 41 deg)`, pairwise meeting with no common point, the non-Helly configuration of
  `CLIQUE_CONTINUUM.md` — where the run's separator had reported `kmax = 1.000000`.  Every
  `kmax` in this note and in `T4LEAF.md` is a statement about the separator's candidate grid,
  not about the family; on the leaf's 400–500-square measures the enumeration does not finish
  (`LEAF_CEILING.md` §2.5), so the clique side of `11.80` remains float.
* **Not bounds.**  Every `pure` column (`11.932869` for `D1111`; the row set is the clique
  run's and `M_pure = 1.054`).  `D1111EQ`'s "infeasible" is a statement about the warm-start
  pose set only.
* **Theorems.**  §1.2 (degeneracy), §1.3 (the cap child closes iff the equality child does),
  §1.4 (balanced child = parent, drop `0`), §1.5 (`q_i <= 3`), and `leaf value = 8 + interior`.
  The drop being zero does not depend on any number in the table.

---

## 5. What this says, and what to do instead

1. **The quadrant branch has no pruning power at `01010101`, exactly.**  Not "small", not "still
   rising": the hardest child is the parent, by symmetry, and its value rises with the parent's
   (`11.745 -> 11.804` here, `T4LEAF.md` §7 item 1's "converge the clique leaf" is what the
   balanced child's run actually did).  The unbalanced children are integral packings `1–2` units
   below and would close trivially.  Under the brief's own rule (`<= 0.1` ⇒ new mathematics) the
   corner + slot + anchor-clique + quadrant family is not a finite engineering project.
2. **Why: the fractional excess is symmetric and central.**  The parent's optimum puts `0.4`
   on each of the four central grid poses, `1.64` where a packing can put at most one square
   that touches the crossing; the `C_4`-symmetric fractional optimum is itself a feasible point
   of the balanced child.  Any `C_4`-invariant count branch has this property.
3. **What would have pruning power.**  A branch the symmetric average *cannot* land on:
   * an **odd** interior grid — `3 x 3` cells of side `2/3` on `[1,3]^2`.  The four squares
     among nine cells, up to `C_4`: the symmetric fractional optimum spreads over all nine, but a
     `C_4`-invariant integer child is `(1,1,1,1 | 0,0,0,0 | 0)` (one per corner cell) or
     `(0,0,0,0 | 1,1,1,1 | 0)` (one per edge cell); no child inherits the parent's optimum, and the
     centre cell (which holds the crowding) can carry at most one.  The number of children is
     `(number of 4-subsets/multisets of 9 cells with the cell capacities) / 4`, a few dozen; the
     same `--quad` machinery with a `3 x 3` label map runs it unchanged;
   * a branch on **which square sits at the centre** — "exactly one centre in `[1.5, 2.5]^2`,
     and its angle bin" — which breaks the symmetry by hand and pins the one square whose
     fractional smearing is the whole excess;
   * or the other two `m = 4` patterns first (`D_2`-invariant, where the quadrant branch is not
     killed by the argument), to learn whether the family closes anything at `t = 4`.
4. **The instrument.**  (a) The separator misses anchor cliques of mass `3/2` on a 24-square
   measure and `1.166` on the leaf's own 495-square measure while reporting `kmax = 1.000000`;
   every `kmax` in `T4LEAF.md` and here should be read as "at the candidate grid's resolution",
   and `leaf_ceiling.py check --anchor scan` (3.5 h single-threaded on 495 squares) should be
   run on every measure whose clique value is quoted — better, its exact witnesses should be
   fed back as rows (`--cq-load` accepts anchors), which is the cheapest way to make the clique
   leaf honest.  (b) `snap --polish` re-optimises without the
   clique rows and can therefore *raise* a clique-run measure into clique-infeasible territory
   (`11.804 -> 11.872`); to certify a clique run's own measure, snap it without `--polish`.  (c)
   `leaf_ceiling.py` now honours a `t4leaf` measure file's region bookings (choice semantics) —
   with the sum semantics the mid-wall poses every optimum here uses make the exact region
   check infeasible.  (d) With `--threads > 1` the search is not reproducible: the LP is
   degenerate and a parallel solve lands on a different optimal vertex; single-threaded it is
   bit for bit.
5. **The decisive computation is still the cover side** (`T4LEAF.md` §7 item 2), and nothing
   here changes that; what changes is that the next level of *this* branch is known not to help
   before that computation is attempted.

---

## 6. Reproduce

```sh
# regression: the flag off reproduces T4LEAF.md's script bit for bit (single-threaded; with
# --threads 2 the SAME script diverges from itself at the second inner iteration -- the LP is
# degenerate and the parallel solve picks a different optimal vertex, hence different rows)
python3 search/t4leaf.py 4.0 REG --corners 1111 --patterns 01010101 --price-pattern 01010101 \
    --cliques --cq-wall --cq-interior --cq-want 40 --cq-max 200 --cq-age 5 --chord --variants \
    --bnd-delta 1e-6 --stages 2 --rowloops 4 --pose-max 2000 --cg-want 150 --seed-pitch 0.5 \
    --seed-dth 30 --row-pitch 0.16 --corner-rows 300 --load-poses runs/reg_poses.txt \
    --load-rows runs/tl_A01010101_dual.txt --cq-load runs/tl_A01010101_cliques.txt \
    --time 100000 --threads 1 --lp-tlim 60 --ckpt-secs 100000
# runs/reg_poses.txt = the 394 support poses of tl_A01010101_poses.txt + 1200 random others
# (random.seed(1)); against the pre-change script the logs differ only in the args line and
# timings, and runs/tl_REG.json `hist` agrees on LP, M, kmax, cols, rows, cq, interior,
# corners, slots, strips, bnd, LP_pure, M_pure.  Smoke: --quad 1111 --quad-mode eq on the same
# instance is infeasible; --quad-mode le gives 11.754 with quads ~0.94 and cap duals 0.

# the children (runs/launchD.sh TAG QUAD MODE THREADS SECONDS [STAGES]; runs/launchDseq.sh)
sh runs/launchD.sh D1111 1111 le 2 18000
for Q in 0211 0022 0112 0220 0121 0013 0031 0130; do sh runs/launchD.sh D$Q $Q le 1 2400 3; done
# = python3 search/t4leaf.py 4.0 TAG --corners 1111 --patterns 01010101 --price-pattern 01010101
#     --quad QUAD --quad-mode le --cliques --cq-wall --cq-interior --cq-max 700 --cq-age 5
#     --cq-want 200 --chord --variants --variant-every 3 --bnd-delta 1e-6 --lp-tlim 20
#     --load-poses runs/tl_A01010101_poses.txt --load-poses-c4 runs/tl_A01010101_poses.txt
#     --load-rows runs/tl_A01010101_dual.txt --cq-load runs/tl_A01010101_cliques.txt
#     --pose-max 8000 --row-cap 3000 --rowloops 25 --stages 400 --time SECONDS --threads N
# the literal equality child (infeasible at the first solve): --quad 1111 --quad-mode eq --stages 1
python3 search/depth_table.py                       # the table of section 2

# exact certification of a child's measure (runs/snapall.sh TAG QUAD MODE [PROCS] [CLIQUE_SECS])
python3 search/leaf_ceiling.py snap runs/tl_D0211_measure.txt --t 4 --out runs/lc_D0211.txt \
    --corners 1111 --patterns 01010101 --quad 0211 --quad-mode le --polish --procs 1 \
    --Q 100000000 --Dc 100000000
python3 search/leaf_ceiling.py check runs/lc_D0211.txt --t 4 --corners 1111 --patterns 01010101 \
    --quad 0211 --quad-mode le --chord --anchor clique
# the balanced child, without re-optimisation (its own measure) and with (no clique rows):
python3 search/leaf_ceiling.py snap runs/tl_D1111s0_measure.txt --t 4 --out runs/lc_D1111s0_np.txt \
    --corners 1111 --patterns 01010101 --quad 1111 --quad-mode le --Q 100000000 --Dc 100000000
python3 search/leaf_ceiling.py snap ... --polish --rows0 50000 --polish-rounds 60 --out runs/lc_D1111s0.txt
python3 search/leaf_ceiling.py check runs/lc_D1111s0_np.txt ... --anchor scan
# the complete anchor enumeration on a small measure (finds the 3/2 clique on D0220's optimum)
python3 search/leaf_ceiling.py check runs/lc_D0220.txt --t 4 --corners 1111 --patterns 01010101 \
    --quad 0220 --quad-mode le --chord --anchor all --budget 900

# the child count (section 1.6) and the q_i bound witness (section 1.5)
python3 search/level2_capacity.py --t 4.0 --box 1 2 1 2 --k 3 --starts 200
```

`leaf_ceiling.py` changes (all behind `--quad` or a new file feature; `selftest` unchanged and
passing): quadrant boxes `Q0..Q3`; `--quad` / `--quad-mode` on `check` and `snap` (caps as
`A_ub` rows in the polish, as upper bounds in the exact assignment search); `# region N`
bookings of `t4leaf` measure files read and honoured in both commands (`--no-labels` restores
the sum semantics), a booked centre within `--bnd-delta` of a boundary moved onto it exactly;
`pose` lines may carry a trailing `# LABEL`.  The angle lattice must be fine (`--Q 1e8`) for
measures that are packings with `1e-6` gaps: the default `1e-5` closes the gaps and reads
`M = 2`.

Warm starts are `T4LEAF.md`'s `A01010101` checkpoints in
`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a8e11c0e17344cba4/runs/`
(read-only, copied into `runs/`).
