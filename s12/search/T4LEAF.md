# The hardest level-2 leaf at `t = 4`, taken to convergence

Code: `search/t4leaf.py` (subclasses `search/t4screen.py`'s `CLeaf`), `search/t4leaf_table.py`.
Runs live in this worktree's `runs/tl_*` (gitignored, so every number is quoted here).
Nothing in `verify/`, `xcheck.py`, `lean/`, `certificates/` or `branch.py` is touched: this is a
packing-side instrument only.

Read against `search/T4SCREEN.md` (the screen this follows up — its §4.2, §4.3, §5, §7, §8 set the
task), `notes/level2-design.md` §2–§5 (the branch), `search/CLIQUE_CEILING.md` /
`search/CLIQUE_CONTINUUM.md` / `search/RECONCILE.md` (the cliques) and `search/DUAL.md` /
`search/CLOSED4.md` (the closed semantics).

---

## 0. Verdict, up front

**(c) still rising when the budget ran out.  Nothing reached 12: the leaf keeps `0.25` of margin
with the anchor cliques on and `0.08` without.  But the level-2 slot branch turns out to be worth
`0.08 – 0.25` at `t = 4`, not the `0.4 – 1.7` the screen measured, and the corner leaf above it is
already at `11.99` and climbing.**

The numbers, all at `t = 4.0`, closed semantics, `r = 1`, four corner boxes pinned at 1, with the
verifier-compatible boundary semantics of §1.1 and the chord rows of §1.2 present:

| what | value | `M` | `kmax` | interior | bound? |
|---|---|---|---|---|---|
| leaf `01010101`, anchor cliques on (`A01010101` st2) | **`11.745193`** | `1.000000000` | `1.000859` | `3.745193` | up to one clique violated by `8.6e-4` |
| leaf `01010101`, **no cliques** (`C01010101P` st2 / `E01010101NB` st0) | **`11.914756` / `11.920067`** | `1.000000000` | — | `3.914756` / `3.920067` | **yes** |
| corner leaf `k = 4`, slots free, cliques on (`B40K` st2) | **`11.988925`** | `1.000000000` | `1.000000` | `3.988925` | **yes** |

Each is a lower bound on its own quantity over its own pose set (§1), so the leaf's true value is
**at least** `11.745` with cliques and **at least** `11.920` without.  None of the three had
stopped moving: the last increments were `+0.014`, `+0.005` and `+0.007` per pricing stage.

What is new, and what matters for the family:

1. **The slot branch is worth far less than the screen said.**  Against the `k = 4` control on the
   same instrument the branch buys `11.9889 - 11.7452 = 0.244` with cliques and `~0.08` without,
   where `T4SCREEN.md` §4.1 reported `0.4 – 1.7` and `notes/level2-design.md` §4.3 measured `0.489`
   at `t = 3.98`.  The cause is the **pose set**, not the semantics: the screen's leaves ran on pose
   sets its own interior calibration showed to be short by `3.5`.
2. **`T4SCREEN.md` §5's boundary obstacle is real but free.**  `0.47 – 0.68` of the leaf optimum's
   mass sits within `1e-7` of a region boundary and the packing does use the freedom the
   verifier-compatible rule gives it — but the converged value is the same to `0.005` either way
   (§3.2).  So the branch can be made verifier-compatible at no cost in leaf value.
3. **The chord inequality `mu(wall strip) <= 3` buys nothing in an `m = 4` leaf, exactly** — its
   dual is `0` at every iteration of every leaf run, because the leaf's own region equalities
   already pin each wall strip at `3` (§3.1).  In the `k = 4` control it is worth `0.61` on a
   starved row set and `0` once the coverage rows converge.  The hypothesis being proved elsewhere
   is needed for exhaustiveness, but it is not a cut here.
4. **`T4SCREEN.md`'s "`k = 4` pure reads exactly `12.000000` six times" is not a value.**  Those
   matched-pure numbers share the clique run's row set, and the measure that attains them has
   certified maximum coverage `M = 1.5186` and a violated anchor clique `mu(K) = 1.2310` (§4).  It
   is an upper bound on the pure restricted-pose value and a bound on nothing else.  Its anatomy is
   worth having anyway, and it is not "the grid minus four": it is **the `4 x 4` grid smeared**
   (§4).

So the family is not refuted here, and it is not vindicated either.  What the run does establish is
that at `t = 4` the honest level-2 branch is a **small** cut on top of the corner branch, that
almost all of the leaf's remaining margin is bought by the anchor cliques rather than by the slot
branch, and that the corner branch alone is already at `11.99` and climbing towards 12.

---

## 1. What is measured, and the three things that changed

The LP is `t4screen.py`'s (`T4SCREEN.md` §1): over closed unit squares `S` contained in the closed
container `[0,4]^2`,

```
max  sum_S mu_S
s.t. coverage(p) = sum_S mu_S [p in S]        <= 1     for every point p of [0,4]^2
     mu({S : centre(S) in C_i}) = 1                    i = 0..3   (the corner boxes, r = 1)
     mu({S : centre(S) in W_j}) = K_j                  j = 0..7   (the eight wall slots)
     mu(K)        = sum_{S in K} mu_S         <= 1     for every anchor clique K = K(p, A)
     mu(strip_w)  <= 3                                 w = 0..3   (the chord rows, §1.2)
     mu >= 0
```

with `K = 01010101`, the hardest of the four `m = 4` leaves under the corner leaf `k = 4`.  By LP
duality this is the cover-side leaf value of a branch certificate with those regions and per-region
multipliers, strengthened by anchor-clique columns.  **A leaf closes iff this value is `< 12`.**
The semantics of the container, the squares, the incidence test, the angle lattice, the exact
arrangement-vertex certification and the anchor-clique membership test are all `T4SCREEN.md` §1.1's,
unchanged; so are the error directions of `T4SCREEN.md` §1.2 (poses restrict and LOWER, rows relax
and RAISE, cliques relax and RAISE), and so is the reading rule:

> **`LP` is a value the poses actually attain only when `M <= 1` and `kmax <= 1`.**  When both
> converge, `LP` is a genuine *lower* bound on the leaf's true value over all poses, and a
> converged value `>= 12` is a decisive no-go for this family.

In an `m = 4` leaf the arithmetic is exact and unchanged (`T4SCREEN.md` §4.2): the four corner boxes
carry 1 each and the four occupied slots carry 1 each, so

> **leaf value = 8 + (mass in the interior `[1,3]^2`)**, and the leaf fails to close iff the
> interior can carry `>= 4`.

Every leaf record below confirms it to twelve digits (corners `[1,1,1,1]`, slots
`[0,1,0,1,0,1,0,1]`).

Three things are new.

### 1.1 Boundary semantics: a pose on a region boundary may choose either side

`T4SCREEN.md` §5's obstacle: the packing optimum parks mass — in two leaves a whole unit at a single
pose — exactly on a slot boundary, `c_x` or `c_y` equal to `t/2 = 2` to within `1e-14`, and it lands
in the *occupied* slot rather than the empty one only because `level2_regions.classify` breaks the
tie that way.  No cell-based verifier can implement a `1e-14` tie-break: `verify/`'s region test is
"is this cell's bounding box inside the box?", exact for convex cells, and a cell straddling the
boundary is inside **neither** box.  So the branch is exhaustive only if a straddling cell is
required to satisfy **both** regions' thresholds.

The packing-side dual of that requirement is: **a pose within `delta` of a region boundary may be
assigned to either adjacent region, at the packing's choice.**  Cover side, the constraint on a
straddling pose becomes `capture(S) >= 1 + max(lam_j, lam_j')` over the adjacent regions `j, j'`;
that is two constraints, and the dual of two constraints is two columns.  `t4leaf.py` therefore
**duplicates such a pose, one column per adjacent region** — identical coverage coefficients,
identical clique coefficients, different region label — with `delta = 1e-6`.  The region boundaries
in the centre plane are `c in {r, t/2, t-r} = {1, 2, 3}` in each coordinate, and the label of each
copy is read off by evaluating `classify` at a point strictly inside the relevant side, so the
answer does not depend on `classify`'s own tie-break at all.

Adding columns can only **raise** the LP value, so the number stays a lower bound on what a
certificate must beat.  §3.2 measures what it costs the branch.

### 1.2 The chord inequality `mu(wall strip) <= 3`

`notes/level2-design.md` §2.3 (Nagamochi 2005 Lemma 7(i) / Stromquist's chord lemma): in a packing
of squares of side `L > 1` inside `[0,t]^2` with `t <= 4`, at most **3** squares have their centres
within 1 of a given wall.  It is **still a hypothesis in this repo** — it is a strictness statement,
false for closed unit squares, and is being proved separately — so it is behind a flag and every
value below is reported with and without it.

> **Correction (2026-09-22).**  No longer a hypothesis: `lean/Sqpack/Chord.lean` `wall_strip_le_three` proves it for closed unit squares pairwise disjoint as closed sets in `[0,t]^2`, `t <= 4` (the semantics used here; `notes/chord-lemma.md`). Only the `1e-9` coordinate tolerance in the rows below is argued outside Lean (the parenthesis after them).  Values reported "with the chord flag" are therefore values of a proved relaxation.

The four rows are built from the pose **coordinates**, not from the region labels:
`mu({c_y <= r + 1e-9}) <= 3` and its three images.  (Including a pose at `c_y = 1 + 1e-9` is safe:
the chord argument only needs `c_y <= 0.9 + w/2`, and `w >= 1` always.)  This matters, because from
the labels alone the inequality is vacuous in an `m = 4` leaf — see §3.1.

### 1.3 One long process, and a HiGHS backend so that it can converge

`T4SCREEN.md` §8: at 100–310 s per LP a 450 s chunk never completed a stage, so the follow-up run
never checkpointed and each chunk started over.  `t4leaf.py` runs as one detached process with
checkpoints (poses, cliques, point duals, the measure, and each variant's measure) at every stage
and at least every 240 s, and replaces scipy's `linprog` — which rebuilds and re-solves the whole
model every time — with a `highspy` model **kept alive across the row and clique generation loop**:
interior point with crossover for the first solve after a rebuild, then warm-started simplex for
every incremental solve (rows appended, cliques added or deleted).  Row ageing is **off**.
`--load-rows` seeds a restart's row set from the previous run's point-dual checkpoint, which is what
makes a restart cheap: only the rows carrying dual mass matter (`270` of `22 796` on the leaf).

Same LP, same answer.  Cross-checked against the scipy path on identical models — leaf and control,
cliques on/off, chord on/off, and after a clique deletion plus a row addition — agreement to
`1e-11`:

```
01010101 cliques+chord    highs 12.000000000  scipy 12.000000000   diff +6.4e-14
01010101 chord, no cq     highs 12.000000000  scipy 12.000000000   diff +6.4e-14
........ cliques+chord    highs 12.105336845  scipy 12.105336845   diff -1.1e-14
........ cliques, nochord highs 12.719713277  scipy 12.719713277   diff +3.9e-14
+300 cliques, cq+chord    highs 11.814161688  scipy 11.814161688   diff +3.3e-13
after -300 cq, +500 rows  highs 12.000000000  scipy 12.000000000   diff +6.8e-14
```

Two further changes affect what the reported number *means*.

* **A settle solve.**  `t4screen.py`'s inner loop calls the clique separator **after** its last
  solve, so the record it reports is a value of a model that no longer exists (its clique rows have
  grown).  `t4leaf.py` ends each stage with a solve of the model as it then stands, recomputes `M`
  and `kmax` against that measure (separating with `want = 0`, i.e. finding the worst clique but
  adding nothing), and values the matched variants on exactly those rows and columns.  One stage
  row is one self-consistent LP.
* **The matched variants are solved on a one-shot model.**  Widening the clique rows' bounds in
  place leaves the warm basis holding nonbasic-at-bound logicals whose bounds no longer exist, and
  HiGHS spends longer repairing that than solving the model outright (26 minutes on a `12.8k x 22.8k`
  instance whose base model solves in 4 s).  Variants build a fresh model and use ipm+crossover.
  The no-chord variant needs no solve at all when the chord rows carry zero dual: that dual solution
  is feasible for the model without them and has the same objective, so `LP_nochord = LP` exactly.

---

## 2. The trajectory

Five processes, all `t = 4.0`, closed semantics, `r = 1`, all four corner boxes pinned at 1, row
ageing off, one long detached process each (`runs/launch{A,B,C,D,E}.sh`, logs `runs/tl_*.log`,
checkpoints `runs/tl_*_{poses,dual,cliques,measure}.txt`):

| tag | pattern | cliques | chord rows | boundary duplication |
|---|---|---|---|---|
| `A01010101` | `01010101` | yes (wall + interior) | yes | yes |
| `C01010101P` | `01010101` | **no** | yes | yes |
| `E01010101NB` | `01010101` | **no** | yes | **no** (`t4screen.py`'s tie-break) |
| `B40K` | `........` (`k = 4` control) | yes | yes | yes |
| `D40KNC` | `........` | yes | **no** | yes |

`A` and `B` were warm-started from `T4SCREEN.md`'s refined interior poses (`t4_INTP_poses.txt` —
the `34–60 deg` interior configurations a seed grid misses) plus the screen's own leaf / control
pose and clique checkpoints; `C` from `A`'s, `E` from `C`'s, `D` from `B`'s.  A "stage" is one
pricing round, taken after the row and clique generation loop has run to its limit and the model
has been re-solved (the settle solve of §1.3).  `conv` = `M <= 1 + 1e-9` **and** `kmax <= 1 + 1e-6`
at that settle solve, i.e. the value is then a genuine lower bound on the leaf over its own pose
set.

```
tag          st pattern          LP         M     kmax       pure    nochord      int          strips        cols  poses   rows   cq     bnd conv
C01010101P   -1 01010101  11.874968  1.000000        -          -  11.874968  3.87497 [3.000 3.000 3.000 3.000] 12828  12684  11411    0       -  YES
A01010101    -1 01010101  11.708426  1.018857  1.00605  11.883716  11.708426  3.70843 [3.000 3.000 3.000 3.000] 12828  12684  30090  206  0.4676  no
A01010101     0 01010101  11.713896  1.001972  1.01299  11.903293  11.713896  3.71390 [3.000 3.000 3.000 3.000] 10348  10156  15919  228  0.5096  no
A01010101     1 01010101  11.730895  1.020709  1.00385          -  11.730895  3.73090 [3.000 3.000 3.000 3.000]  8037   7819  23184  315  0.6771  no
A01010101     2 01010101  11.745193  1.000000  1.00086          -  11.745193  3.74519 [3.000 3.000 3.000 3.000]  8050   7808  27251  162  0.6222  no*
C01010101P    0 01010101  11.894316  1.000000        -          -  11.894316  3.89432 [3.000 3.000 3.000 3.000]  8497   8316  12877    0  0.5440  YES
C01010101P    1 01010101  11.909254  1.000000        -          -  11.909254  3.90925 [3.000 3.000 3.000 3.000]  8117   7910  17373    0  0.5524  YES
C01010101P    2 01010101  11.914756  1.000000        -          -  11.914756  3.91476 [3.000 3.000 3.000 3.000]  8132   7915  19543    0  0.5734  YES
C01010101P    3 01010101  11.920049  1.024053        -          -          -  3.92005 [3.000 3.000 3.000 3.000]  8071   7857  24256    0  0.4996  no
E01010101NB   0 01010101  11.920067  1.000000        -          -  11.920067  3.92007 [3.000 3.000 3.000 3.000]  8341   8341  13289    0  0.2810  YES
E01010101NB   1 01010101  11.922396  1.019849        -          -  11.922396  3.92240 [3.000 3.000 3.000 3.000]  8000   8000  19375    0  0.2775  no
B40K          0 ........  11.962637  1.001354  1.00040  12.000000  11.962637  3.96264 [3.000 3.000 3.000 3.000] 13006  12872  25836  293  0.1042  no
B40K          1 ........  11.977428  1.000000  1.00000          -          -  3.97743 [3.000 3.000 3.000 3.000]  8006   7875  42294  232  0.0840  YES
B40K          2 ........  11.988925  1.000000  1.00000          -          -  3.98893 [3.000 3.000 3.000 3.000]  8018   7882  55580  341  0.0708  YES
B40K          3 ........  11.995681  1.000508  1.00000  12.000000  11.995681  3.99568 [3.000 3.000 3.000 3.000]  8017   7905  63926  201  0.0641  no
D40KNC        0 ........  11.992299  1.000201  1.00067  12.000000     (none)  3.99230 [3.000 3.000 3.000 3.000]  8670   8517  18327  285  0.0668  no
```

`st = -1` is the stage of the run before the last restart, on the larger (12 828-column) pose set.
`bnd` is the mass whose centre lies within `1e-7` of a region boundary; `strips` is the mass of
centres within 1 of each wall (the chord quantity); `int` is the interior mass.  `LP = 8 + int`
holds to twelve digits in every leaf row, with region split `corners [1,1,1,1]` and
`slots [0,1,0,1,0,1,0,1]` exactly.  `*` : `A` stage 2 has `M = 1.000000000` exactly and
`kmax = 1.000859`, i.e. the coverage rows are converged and one separated anchor clique is still
violated by `8.6e-4`.

**Where it stopped.**  All three main runs were still rising when their wall-clock budget ran out:

```
A (leaf, cliques)              11.7139 -> 11.7309 -> 11.7452              +0.0170 +0.0143
C (leaf, no cliques)           11.8943 -> 11.9093 -> 11.9148 -> 11.9200  +0.0150 +0.0055 +0.0053  (11.8750 before)
E (leaf, no cliques, no bnd)   11.9201 -> 11.9224                         +0.0023
B (k = 4 control)              11.9626 -> 11.9774 -> 11.9889 -> 11.9957   +0.0148 +0.0115 +0.0068
```

Compute actually spent: about `5.5` hours of wall clock on `<= 8` threads (three to five
single-threaded HiGHS processes at a time, plus their multithreaded geometry kernels).  Restarts:
`A` three times, `C` twice, each from its own pose + point-dual + clique checkpoint, for the two
solver pathologies of §1.3.  Nothing was lost — a restart reloads the pose set and the binding
rows, so the trajectory above is continuous across them.

---

## 3. Reading

### 3.1 The chord inequality buys nothing in an `m = 4` leaf, and nothing in the control either

In every leaf record the chord duals are `[0, 0, 0, 0]` and the wall-strip masses are
`[3, 3, 3, 3]` exactly, so `LP_nochord = LP` exactly (§1.3's dual argument).  The reason is
arithmetic, not numerical: with `r = 1`, `t = 4` the bottom strip `{c_y < 1}` **is**
`C_0 u C_1 u W_0 u W_1` as a set of centres, and the leaf pins
`mu(C_0) + mu(C_1) + mu(W_0) + mu(W_1) = 1 + 1 + K_0 + K_1 = 3`, because `K_0 + K_1 = 1` in every
`m = 4` pattern.  The only way the row can bite is mass sitting exactly on `c_y = 1` and labelled
interior or side wall, which the boundary duplication of §1.1 makes possible — and the LP never
uses it.

In the `k = 4` control the slots are free, so the inequality is not implied, and on a **starved**
row set it is worth a lot: on the calibration instance of §1.3, `12.105337` with the chord rows
against `12.719713` without — a cut of `0.614`.  Once the coverage rows converge that cut is gone:
`B40K` stages 0 and 3 both give `LP_nochord = LP` to all printed digits (`11.962637`, `11.995681`),
and the control's own optimum puts exactly `3` in every wall strip unprompted.  `D40KNC` is the
same control run with the chord rows removed from the model altogether, warm-started from `B`'s
checkpoint: it reaches `11.992299` (`M = 1.000201`, `kmax = 1.000666`) on a third of `B`'s rows,
i.e. the same place `B` is, from below.  **At `t = 4` the
chord bound is implied by the closed-square coverage constraint at the optimum.**  It remains what
makes the leaf *set* exhaustive (`notes/level2-design.md` §2.3) — that is not an LP statement — but
as a cut it is worth zero here.

### 3.2 The verifier-compatible boundary semantics costs the leaf almost nothing

`T4SCREEN.md` §5 found the leaf optima parking `1 – 11 %` of their mass — in two converged leaves a
whole unit at a single pose — exactly on a slot boundary, with the branch decided by
`level2_regions.classify`'s `1e-14` tie-break.  That is a real obstacle for a cell-based verifier
(§1.1), and the mass is still there: `0.47 – 0.68` of the leaf optimum sits within `1e-7` of a
region boundary, and the duplication is used — in `C`'s stage-2 optimum `0.18` of mass sits at
`c_x = 2` labelled `W_5`, where `classify` would have said `W_4` (empty), and `0.095` likewise at
the left wall.  (`classify` happens to favour the *occupied* slot on the bottom and right walls and
the *empty* one on the top and left walls, so only two of the four walls are affected.)

But it costs the value nothing measurable.  `E01010101NB` is the same no-clique leaf run with
duplication switched off, warm-started from `C`'s pose set:

```
C (duplication on)   stage 2   11.914756   converged, boundary mass 0.5734
E (duplication off)  stage 0   11.920067   converged, boundary mass 0.2810
```

The two agree to `0.005`, and `E` — the *restricted* semantics — is the higher of the two.  (Their
pose and row sets differ, so this is not a matched pair; on an identical pose set duplication can
only raise the value.)  With a rich pose set the LP simply uses a near-neighbour pose instead.
**So `T4SCREEN.md` §5's obstacle is an implementation problem for the verifier, not a source of leaf
value: making the branch verifier-compatible is free.**

### 3.3 What the slot branch is worth at `t = 4`, honestly

| | with anchor cliques | without cliques |
|---|---|---|
| corner leaf `k = 4`, slots free, best **certified** (`B40K` st2) | `11.988925` (`M = 1`, `kmax = 1`) | `12.000000` (**not** a bound — §4) |
| corner leaf `k = 4`, last stage (`B40K` st3) | `11.995681` (`M = 1.0005`) | `12.000000` (not a bound) |
| leaf `01010101`, best (`A` st2 / `E` st0) | `11.745193` (`M = 1.000000`, `kmax = 1.00086`) | `11.920067` (`M = 1.000000`) |
| **slot branch is worth** | **`0.244` – `0.251`** | **`~0.08`** |

`T4SCREEN.md` §4.1 put the same branch at `0.4 – 1.7`, and `notes/level2-design.md` §4.3 measured
`0.489` for this leaf at `t = 3.98`.  The collapse is a pose-set effect, not a semantics effect
(§3.2): the screen's leaf runs sat `0.4 – 1.0` below 12 on `3.4k–6.6k` poses that under-measured
the interior by `3.5` on the screen's own calibration.  Feeding the leaf the refined interior poses
moves the leaf up by `~0.15` and the *control* up by `~0.10`, and what survives as a difference is
`0.25`.

The matched clique gain — same rows, same columns — is `+0.19` (`A` stage 0: `11.903293` pure
against `11.713896`), against `+0.02 – 0.08` for the corner leaf and `+0.25` in the screen.  **At
`t = 4` the anchor cliques now do more of the work on this leaf than the slot branch does.**

### 3.4 The threshold, and how far from it

`leaf value = 8 + interior mass`, so the leaf fails to close iff the interior `[1,3]^2` carries
`>= 4`:

```
interior mass, leaf 01010101, t = 4, closed semantics, corner leaf k = 4
  with anchor cliques   3.745193   (coverage certified; one clique violated by 8.6e-4)   deficit 0.2548
  without cliques       3.920067   (certified)                                           deficit 0.0799
for reference, k = 4 with the slots free   3.995681   (M = 1.0005)                       deficit 0.0043
```

`T4SCREEN.md` §4.2's calibration disaster — the interior sub-problem reading a "converged"
`4.000000` on a seed pose set and `7.53` after pricing — does not recur here in the same form:
these runs *start* from that pricing's output, the interior mass is `3.75 – 3.92` rather than
`2.3 – 3.5`, and each stage's rise is `0.005 – 0.017` rather than `0.05`.  But the direction is the
same and it has not stopped.

### 3.5 Which numbers here are what

* **Bounds.**  `C01010101P` stages 0–2, `E01010101NB` stage 0 and `B40K` stages 1–2 have
  `M = 1.000000000` and (where cliques are on) `kmax = 1.000000`: genuine lower bounds on their
  leaf's value over their own pose sets, hence on the true leaf value.  `A01010101` stage 2 has
  `M = 1.000000000` and `kmax = 1.000859`, so it is a bound up to one anchor clique violated by
  `8.6e-4`.
* **Not bounds.**  Every `pure` column — including the `12.000000` of `B40K` — shares the clique
  run's row set and its measure is not certified (§4).  `B40K` stage 3 and `E` stage 1 have `M`
  just above 1.
* **Float, not exact.**  Everything above is a float LP.  What is exact is the arithmetic
  `leaf value = 8 + interior` and `T4SCREEN.md` §4.2's argument that no *integral* witness of mass
  12 can exist (it would refute `s(12) = 4` outright, and `s(13) = 4` is a theorem).

---

## 4. The `k = 4` measure at exactly `12.000000`

`T4SCREEN.md` §8 records the corner leaf's *matched pure* value landing on exactly `12.000000` in
six independent runs.  It does so here too, twice more — `B40K` stages 0 and 3, both to `1e-14`,
with interior mass `4.000000000000002` and `4.000000000000132`.  Two things about it.

**It is not a value.**  The matched pure solve uses the clique run's row set, and the measure that
attains `12.000000` has certified maximum coverage `M = 1.5186` and a violated anchor clique
`mu(K) = 1.231016` — anchor `p = (2.0, 1.1)`, segment endpoints `(2.06125, 1.325)`–`(1.93875,
1.325)`; the point clique at the same anchor carries only `1.019233`, so this is exactly the
non-Helly mass `CLIQUE_CONTINUUM.md` describes.  So it is an upper bound on the pure restricted-pose
value and a bound on nothing else, and the same caveat applies to every `LP_pure` column in
`T4SCREEN.md` §3.

**Its anatomy: the `4 x 4` grid, smeared.**  `runs/tl_B40K_measure_pure_nochord.txt` — 288 support
poses, `mass = 12.000000000`, `4.000000 / 4.000000 / 4.000000` over corners / slots / interior,
`0.475 – 0.583` per slot, exactly `3.000000` in every wall strip:

```
mu= 0.988174  (3.50000, 3.50000)    0.00 deg  C3      mu= 0.240856  (2.49737, 1.50129)   0.00 deg  I
mu= 0.984833  (0.50000, 0.50000)    0.00 deg  C0      mu= 0.211444  (2.49516, 2.49763)   0.00 deg  I
mu= 0.967053  (3.49456, 0.50179)    0.00 deg  C1      mu= 0.187861  (1.50079, 2.49873)   0.00 deg  I
mu= 0.927020  (0.50329, 3.50000)    0.00 deg  C2      mu= 0.156118  (1.50100, 1.50834)   0.00 deg  I
mu= 0.367064  (3.49916, 1.50073)    0.00 deg  W2      mu= 0.140872  (2.75213, 2.55748) -32.85 deg  I
mu= 0.339476  (2.50000, 3.50000)    0.00 deg  W4      mu= 0.139160  (0.50000, 2.13294)   0.00 deg  W6
mu= 0.338047  (0.50153, 2.49947)    0.00 deg  W6      mu= 0.123274  (1.87027, 0.50145)  -0.20 deg  W0
mu= 0.310452  (2.50000, 0.50000)    0.00 deg  W1      mu= 0.116875  (1.42127, 1.28710)  27.60 deg  I
mu= 0.306086  (1.51851, 0.50051)    0.00 deg  W0      mu= 0.109734  (1.50490, 1.42690)  -8.35 deg  I
mu= 0.305526  (1.50080, 3.50000)    0.00 deg  W5      mu= 0.101605  (2.71995, 1.42521)  26.00 deg  I
mu= 0.303743  (3.49981, 2.49921)    0.00 deg  W3      mu= 0.099680  (1.42505, 2.72765) -31.00 deg  I
mu= 0.302911  (0.50000, 1.50000)    0.00 deg  W7
```

`7.744652` of the mass is axis-parallel and `4.255348` tilted; the sixteen `4 x 4` grid poses carry
about `7.2` of the `12` between them and the rest is spread over near-grid and `26 – 33 deg` tilted
poses.  **It is not the grid minus four.**  It is the whole grid with the mass graded by how
crowded each grid square is: a corner grid square touches two neighbours and gets `~0.98`, a wall
grid square touches three and gets `~0.31`, a central grid square touches four and gets `~0.20`,
and the `4.76` the crowding removes is redistributed onto tilted interior and wall poses.  (Mass 1
on each of the sixteen — or on each of a "grid minus four" twelve — is *infeasible* for this LP:
two closed unit squares sharing an edge give coverage 2 on that edge.  That is exactly why `t = 4`
is non-trivial and `nu_f(4) = 12.163`.)

The leaf's own converged optimum (`runs/tl_C01010101P_measure.txt`, `11.914756`, 431 support poses)
has the same shape with one slot per wall switched off: `0.87 – 0.99` at the four corner grid
poses, `0.49 – 0.54` at the four *occupied* wall grid poses, `0.21` at each of the four central grid
poses, and `4.8` of tilted mass.

**Clique feasibility.**  The `12.000000` measure is **not** clique-feasible (`mu(K) = 1.231016`) and
not coverage-feasible (`M = 1.5186`).  It is the same object `CLIQUE_CEILING.md` and `RECONCILE.md`
chase — the frame driven to `4 + 4` with the residue in the interior — and the clique family does
cut it: by `0.037` on `B40K` stage 0 (`12.000000 -> 11.962637`) and by `0.19` on the leaf.

---

## 5. Verdict

**(c) still rising when the budget ran out.**  The hardest level-2 leaf at `t = 4`, pattern
`01010101` under the corner leaf `k = 4`, with anchor cliques, the chord rows and the
verifier-compatible boundary semantics, reached **`11.745193` with the coverage rows certified
exactly (`M = 1.000000000`) and one anchor clique still violated by `8.6e-4`** — interior mass
`3.745193`, deficit `0.2548` against the threshold of 4 — and was rising by `+0.014` per pricing
stage after three stages, the increments decaying from `+0.017`.  Without the cliques the same leaf
reaches a fully certified `11.920067` (deficit `0.0799`), rising `+0.002 – 0.006`.  The corner leaf
`k = 4` alone, on the same instrument, reaches a certified `11.988925` and is at `11.995681` and
still climbing: on the packing side the corner leaf really is sitting on 12 at `t = 4`, exactly as
the level-1 tree says.  **So the value does not reach 12 here, and the leaf keeps about a quarter of
a unit of margin with the cliques on; but it is not "clearly below and converged" either — the
margin has shrunk monotonically under every refinement, and the no-clique leaf value is already
within `0.08` of 12.**  On this evidence the corner + slot + anchor-clique family is **not refuted**
at its hardest leaf and remains plausible, with essentially all of the remaining margin coming from
the anchor cliques (`+0.19`) rather than from the slot branch (`+0.08` on the no-clique side); the
level-2 branch, once its region test is made verifier-compatible, is a much smaller cut at `t = 4`
than the screen suggested, and the chord lemma — the thing that makes the fifteen leaves
exhaustive — is worth exactly zero as an inequality here.

---

## 6. Reproduce

```sh
# the hardest leaf: anchor cliques + chord rows + boundary semantics
python3 search/t4leaf.py 4.0 A01010101 --corners 1111 --patterns 01010101 \
    --price-pattern 01010101 --cliques --cq-wall --cq-interior --cq-max 700 --cq-age 5 \
    --cq-want 200 --chord --variants --variant-every 3 --bnd-delta 1e-6 --lp-tlim 20 \
    --load-poses POSES.txt --load-rows DUAL.txt --cq-load CLIQUES.txt \
    --pose-max 8000 --row-cap 3000 --rowloops 25 --stages 400 --time 6300 --threads 4

# the same leaf without cliques (C: drop --cliques ...), and with t4screen.py's original
# region tie-break instead of the boundary semantics (E: --bnd-delta -1)
# the k = 4 control (B): --patterns "........" ;  without the chord rows (D): drop --chord

python3 search/t4leaf_table.py --table runs/tl_*.json
python3 search/t4leaf_table.py --anatomy runs/tl_B40K_measure_pure_nochord.txt --cliques
```

Warm starts come from `T4SCREEN.md`'s checkpoints in
`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a29754bcd81341ece/runs/`
(`t4_INTP_poses.txt`, `t4_L01010101_poses.txt`, `t4_K40_poses.txt` and the matching
`_cliques.txt`), read-only, copied into this worktree's `runs/`.

A single `t4leaf.py` invocation is a pure function of its arguments except for the `--time` and
`--lp-tlim` cuts, both of which land mid-loop; a run with both large and `--stages` / `--rowloops`
as the only stopping rule reproduces bit for bit.  `--backend scipy` reproduces `t4screen.py`'s
solver path and agrees with the HiGHS path to `1e-11` on identical models.

---

## 7. What a follow-up would need

1. **Converge the clique leaf.**  `A01010101` stopped at `11.745193` after three pricing stages
   with increments still `+0.014`.  Another `10 – 20` stages would settle whether it flattens near
   `11.85` or keeps going.  The cost is bounded now: with the HiGHS backend a stage is `10 – 25`
   minutes at `8k` columns and `25k` rows, against a stage that never completed at all in
   `T4SCREEN.md` §8.
2. **The cover side.**  Nothing on this side can *prove* a leaf closes; a value `< 12` here is the
   optimistic direction.  The decisive computation is `branch.py --cq-interior` at `t = 4` with the
   twelve region multipliers and the verifier as separation oracle
   (`notes/level2-design.md` §8: `13 – 27` h per leaf).
3. **The other three `m = 4` leaves, and the tree under `k < 4`.**  `01010101` was the hardest of
   the four on the screen's numbers; whether it still is once every leaf gets the refined interior
   poses is unmeasured.  The tree under `k = 0,1,2,3` is untouched.
4. **The region test in `verify/`.**  §3.2 says the boundary semantics is free in value, so the
   remaining work is purely to implement it: a cell straddling a slot boundary must be required to
   meet both slots' thresholds.  That is a small change to the cell-vs-box test and it makes the
   branch exhaustive without a tie-break.
5. **The chord lemma** is still needed for exhaustiveness (`notes/level2-design.md` §2.3, §6) even
   though §3.1 shows it is worth nothing as a cut.
