# unavoid13-no: is `h(F) >= 14` reachable for closed unit squares in `[0,4]^2`?  (2026-09-22)

Task brief: `tasks/unavoid13-no/README.md`.  Previous note: `notes/unavoid13.md`.  Code (new, all
importing `search/unavoid13_lib.py` unmodified): `search/unavoid13no_lib.py` (columns, exact/lenient
incidence, D4 action, HiGHS wrappers), `search/unavoid13no_cells.c` (+ `.so`, maximal arrangement cells
in C/OpenMP), `search/unavoid13no_loop.py` (family loop), `search/unavoid13no_*` (see each section).
Runs: `runs/unavoid13no_*`.  Certificates: none — no `certificates/unavoid13/no13_*` was earned.

## 0. Verdict

**Not reached.**  No finite family `F` with `h(F) >= 14` was found; `h(F) = 13` at every round of every
loop, up to `|F| = 1,761` squares (`≈ 34k` exact maximal cells, LP `12.25`) — except on the two largest
families, where no solver returned any verdict (below).  What this task establishes:

* [proved] the finite cover LP over `tiles ∪ grid2 ∪ COVER4-support` with **all** arrangement cells is
  `12.268804`, equal to the certified packing mass — the LP side is exhausted at the support; every
  family here has LP `12.20–12.27` and the integer question is a pure integrality-gap question `12.2 → 14`;
* [measured] the loop's 13-sets have two recurring shapes (§3) and their worst violation converges to
  the previous task's maximin level `−0.026 … −0.03` within ~10 rounds, after which each round removes
  one more "hardest" 13-set of a continuum of them — the trend does not show how many rounds remain;
* [measured] the practical wall is the finder: exact cells grow by `≈ 60` per added killer square and
  HiGHS's time to find a 13-set goes from seconds at `6k` cells to `> 20–40 min` at `27–30k` (§4);
  a purpose-built B&B (§5) proves a genuine `T = 3` infeasibility in 85 s (HiGHS: 49 s) but is a worse
  finder than HiGHS above `12k` cells;
* [proved, tooling] `search/unavoid13no_cells.c` computes the dominance-reduced cell set of a
  3,000-square family in 2 min (the previous task's pass could not go beyond ~1,000 squares) and
  `search/unavoid13no_recheck.py` re-solves an instance exactly (Fractions) with candidates verified
  complete by an exact witness check, validated on the `T = 3` theorem (`7,361` exact candidates, three
  solvers, `k = 6` infeasible);
* [proved] the brief's LP test cannot detect an incomplete candidate set (§1), and the previous task's
  lenient tolerance `1e-7` silently costs `0.46` of LP once the support is in `F` (§1) — every loop of
  `notes/unavoid13.md` at `T = 4` ran with a weaker LP than it needed to.

No 13-set was certified either (none had a nonnegative float margin; the checker was never reached).
The symmetric-exclusion runs of the previous task (`Rx`, `Rd`) did not finish during this task (§6).

**The two largest families are undecided, not feasible [measured].**  On the union of all killers over the
`support:80` bulk (`2,809` squares, `132,301` cells, LP `12.2495`) HiGHS ran **2 h without finding a
13-set and without excluding one** (dual bound `13.0`, 2 nodes); on L2e's `1,617`-square family the
4-process portfolio found **no 13-set in 1 h** (§4).  Every smaller family had a 13-set within minutes.
This is consistent with `h = 14` on those two families and equally with HiGHS running out of steam at
`30k–130k` columns; nothing is proved.  If the exact re-solve ever returns `INFEASIBLE` on one of them,
`runs/unavoid13no_L2e/F_round01.txt` / `runs/unavoid13no_union/u1_family.txt` are the candidates.

**Update (2026-09-22, 21:40).**  The 7,200 s portfolio on L2e's `1,617`-square family **found a 13-set** (optimisation
form, seed 4, `3,088 s`, LP `12.2057`; `runs/unavoid13no_L2e/round_log.txt` round 2), so that family has `h = 13` and is
no longer a candidate.  L3d's `1,761`-square family also gave a 13-set (`2,661 s`).  Both loops continued (1,697 and
1,841 squares).  Only the `2,809`-square union remains undecided.

**State at hand-off** (detached, still running; do not restart on top of them): `L2e` (`runs/unavoid13no_L2e`,
round 2 = a 7,200 s portfolio on the 1,617-square family; a `13-set by IP` line means feasible,
`INFEASIBLE` means the dump `final_*` exists and must go through `search/unavoid13no_recheck.py`) and
`L3d` (`runs/unavoid13no_L3d`, round 5, 1,761 squares).  Each has 4 spawned HiGHS workers
(`multiprocessing.spawn`), all 1-thread.  The previous task's `symRx` and `t4L` are untouched.

**What would be needed** (§7): a finder/prover that does not pay for every cell, i.e. either (i) a
row-pruned lazy loop keeping the IP family at `≤ 800` rows (`≤ 12k` cells, minutes per round) with the
full `F` used only for the checks, run for days, plus one final long HiGHS/B&B proof on the accumulated
family, or (ii) a genuinely different prover (geometric region branching with column generation).  The
outcome of either is uncertain: the loops at `T = 4` have now consumed `~2,300` killer squares in this
task (`≈ 940` in the L2 lineage, `≈ 760` in L3, `≈ 600` in L1) on top of the previous task's 700, against
200 at `T = 3`, without `h` provably moving.

## 1. Method and soundness

**Objects, lower bound, candidates** — as in `notes/unavoid13.md` §1: closed unit square `Q(cx, cy, u)`,
`u = tan(θ/2)` rational; `h(F) ≤ p(4)` for every finite family `F` of admissible squares; the hitting-set
IP over the vertices of the arrangement of `F` computes `h(F)` exactly, and candidates with identical
incidence sets can be merged, dominated ones dropped.

**Maximal cells in C [proved criterion, measured agreement].**  The previous task's bottleneck was the
candidate set: all arrangement vertices, then an `O(K · posting)` dominance pass (11k candidates from 900
squares, 10+ min IPs).  `search/unavoid13no_cells.c` streams every arrangement vertex `v` and keeps it iff
it is *maximal*, decided locally: with `S(v)` the squares containing `v` and `P(v) = ∩ S(v)`,

> `v` is dominated by some arrangement vertex  ⇔  some square `Q ∉ S(v)` meets `P(v)`.

(⇐: `Q ∩ P(v)` is a convex polygon inside every square of `S(v)` and inside `Q`; its vertices are
arrangement vertices with incidence `⊇ S(v) ∪ {Q}`.  ⇒: a dominating `v'` lies in `P(v)` and in some
`Q ∉ S(v)`.)  `P(v)` is computed by half-plane clipping, the test by a bounding-box pre-filter and a clip.
In floats the incidence is lenient (`tol = 1e-7`, squares enlarged by `tol` in the clips), which can only
keep more vertices — the safe direction.  Survivors are deduped on `S(v)` and passed once more through the
posting-list filter (cheap at that size).  [measured] On the `T = 3` final family (575 squares) it returns
**7,361** candidates in 0.7 s — exactly the exact recheck's count (`notes/unavoid13.md` §2.3) — and
`h = 7`; on `runs/unavoid13_t4b/F_round05.txt` (889 squares) 10,533 candidates (the previous float code:
10,673; same LP `11.3939`), `h = 13`.

**Box relaxation — tried and rejected [measured].**  A first design used a D4-symmetric grid of closed
boxes as columns with `inc(box, Q) = [box ∩ Q ≠ ∅]` (a sound relaxation: every point lies in a box whose
incidence contains its own).  At pitch `0.05` (6,400 boxes, 2,216 after dominance) over the 3,201-square
family of §2 the LP is **9.82** and **11 boxes** hit everything.  The reason is that the problem is sharp at
side exactly 4: `ν_f^closed(3.99) = 12.008` against `ν_f^closed(4) = 12.269` (`DUAL_EXACT.md`, `COVER4.md`),
so any leniency of order `δ` costs of order `100 δ` in the bound.  Exact cells are essential; boxes are
only usable as a coarse first stage far from the tile skeleton, and were not used further.

**The brief's "LP-completeness test" is a check of incidence, not of completeness [proved].**  For any
candidate set `C` the hitting-set LP over `F ⊇ COVER4-support` has dual `max Σ y_Q s.t. Σ_{Q ∋ c} y_Q ≤ 1
(c ∈ C)`, and the certified measure `μ/M` (`COVER4.md`) is dual-feasible for *every* `C` (its coverage is
`≤ 1` at every point of the plane), so `LP ≥ 12.2688038611` holds *whatever* `C` is — fewer candidates
make the LP larger, not smaller.  An LP below `12.2688` can therefore only mean wrong (too lenient)
incidence or a mis-read support; it cannot detect a missing candidate.  The value that *is* informative is
the LP with *all* cells — the finite cover LP over the support — and how far above `12.2688` it sits (§2).

**The lenient tolerance `1e-7` is a trap once the support is in `F` [measured, then explained].**
With `tol = 1e-7` (the previous task's setting) the LP over all cells of `tiles ∪ grid2 ∪ support`
(2,393 squares, 126,649 cells) is **11.806**, below the theorem's `12.2688`.  Coverage of the certified
measure on those cells: maximum **1.541** at the lattice point `(1,1)` (exact coverage there: exactly
`1.000`; 12 support squares pass within `1e-7` of the point without containing it).  The reason is in
`COVER4.md`: the support's centres were snapped to `1e-8` and *clamped*, so many squares deliberately
stop `1e-8` short of a tile corner — a lenient `1e-7` incidence puts the corner inside them, the
columns become too strong and the LP collapses by `0.46`.  All cell computations here therefore use
`tol = 1e-11` (float errors are `~1e-15`, so the incidence is still a superset of the truth).  Anything
in `notes/unavoid13.md` computed with `tol = 1e-7` on a family containing near-coincident support
squares is *sound* (a lower bound stays a lower bound) but *weaker than it needed to be*; its LP values
`11.3–11.6` for families with `support:40`/`support:120` are partly this effect.

## 2. Step 1: the right `F` and its LP  [measured; the inequality is proved]

`python3 search/unavoid13no_lp.py --init tiles,grid2,support --out runs/unavoid13no_lp/support11`
(`runs/unavoid13no_lp/support11.log`): `tiles ∪ grid2 ∪ support` = **2,393** squares (the support's 294
poses give 2,348 distinct squares under D4 — orbits of the corner tile and of poses on symmetry axes are
smaller than 8 — plus 45 axis-parallel grid squares not in the support); 3,682,040 arrangement vertices
streamed, **127,573** maximal cells (`tol = 1e-11`, D4-closed as computed, 0 images missing).

> **LP over all cells = 12.268804** (HiGHS, 1,415 s; dual support 1,991 rows).

This equals the certified mass `24537607710/1999999999 = 12.2688038611` to the printed precision: the finite
cover LP over the support family is *tight* at the theorem's floor — the certified measure is an optimal
dual for this `F` and every arrangement vertex is (weakly) covered.  The test passes; as §1 explains, it
tests incidence, not completeness.  With `tol = 1e-7` the same instance gives `11.806` (§1).

**Cheaper bulks [measured]** (`runs/unavoid13no_test/bulks.log`; `support:N` = the `N` heaviest poses):

| family | squares | maximal cells | LP |
|---|---|---|---|
| tiles ∪ grid2 ∪ support:40 | 361 | 1,521 | 11.981 |
| tiles ∪ grid2 ∪ support:80 | 681 | 5,305 | **12.197** |
| tiles ∪ grid2 ∪ support:120 | 1,001 | 15,957 | 12.249 |
| tiles ∪ grid2 ∪ support:200 | 1,641 | 50,709 | 12.267 |
| A1 (grid 1/5) ∪ support:80 | 1,371 | 10,313 | 12.206 |
| A1 (grid 1/10) ∪ support:80 | 3,051 | 16,969 | 12.219 |
| tiles ∪ grid2 ∪ support (all) | 2,393 | 127,573 | 12.269 |

The 80 heaviest poses carry `12.20` of the `12.27`; the structured `45°` bulk adds only `0.01–0.02` of LP
for `5–12k` extra cells.  The loop's bulks are therefore `support:80` (L2) and `support:120` (L3).

## 3. Step 3: structure — `h` on structured subfamilies  [measured]

`search/unavoid13no_struct.py --grid 10` (`runs/unavoid13no_struct/`).  Each stage is D4-closed; centres
on the `1/10` lattice (clamped exactly into the admissible range); `h` by `highspy` (gap 0).

| stage | squares | maximal cells | LP | `h` | IP time | optimal set |
|---|---|---|---|---|---|---|
| A0: axis-parallel (contains the tiles and every shifted `4×3`, `3×4`, `3×3` tiling at shift `1/10`) | 961 | **441** | 9.000 | **9** | 0 s | `{1, 2.1, 3}²`-like: 4 lattice points, 5 on tile edges |
| A1: + `45°` (`u = 5/12, 7/17`: `45.24°`, `44.76°`) | 2,419 | **1,285** | 10.500 | **11** | 3 s | corners `(1,1),(1,3),(3,1),(3,3)`; `(0.904,2),(2,0.904),(2,3.096),(3.1,2)`; `(2,1.5),(2.5,2.1),(1.452,2.148)` |
| A2: + near-axis `2.0°, 4.1°, 6.0°` (grid 1/5 only: 2,473 sq.) | 2,473 | **63,616** | 10.000 | **10** | 252 s | (grid-1/5 A1 alone: 739 sq., 580 cells, `h = 10`) |
| A3: + Pythagorean `u ∈ {1/8,1/7,1/6,1/5,1/4,1/3,2/5,1/2,3/5,2/3}` | *(pending, grid 1/5)* | | | | | |

Near-axis classes are the worst of both worlds: their edge lines are all distinct (`cos, sin` with
denominator `3250`), so cells explode (`580 → 63,616` at grid 1/5), and they add nothing to `h`
(`10 → 10`); the previous note's "wall squares" matter for *certifying* an upper bound, not for the lower
bound.  The `T = 4` stages A0 → A1 are `h = 9 → 11`; the rest of the way to 13 comes from the general
angles of the support (§2: `support:80` alone has LP `12.20`).

**The shape of the 13-sets [measured, `search/unavoid13no_skeleton.py`].**  Over the 12 polished 13-sets
of the loops L1–L3 (§4) two skeletons recur, D4-images of each other:

* **(A) rows `4 + 3 + 4` on `y ≈ 1, 2, 3` plus 2 free points** at `y ≈ 1.4–1.6` and `2.4–2.5` (x-pattern
  `4+1+3+1+4` at `x ≈ 1, 1.5, 2, 2.5, 3`): 11 of the 13 points within `0.12` of the lines `y ∈ {1,2,3}`.
  This is Friedman's `4+3+3+4` with the two middle rows (`y = 1.8, 2.2`) merged into one row of 3 at
  `y = 2` and the two saved points spent as free points in the gaps.
* **(B) columns `3 + 4 + 3 + 3`** at `x ≈ 0.9, 1.75, 2.55, 3.1`: three columns with points at `y = 1, 2, 3`
  and one column of 4 at `y ≈ 1.0, 1.55, 2.45, 3.0` (the previous task's `t4L` sets have this shape).

Both have `≥ 10` points on the tile lines `x ∈ {1,2,3}` or `y ∈ {1,2,3}`; the free points sit near the
tile centres' vertical/horizontal lines.  No rigid skeleton shared by *all* sets was found (L3 round 0
has neither), so a skeleton/no-skeleton split was not attempted; what the two shapes say is that a
13-set needs 9 points doing the axis-parallel work and only 4 for everything rotated, and the killers
found each round are rotated squares (`20°–70°`) through the gaps between the free points.

**Why structured families are cheap [measured, explained].**  With centres on a lattice and few angles,
edges of different squares are collinear (axis-parallel: only `2 × 31` distinct edge lines), so the
arrangement has few distinct vertices and almost every vertex is dominated: 961 axis-parallel squares
give 441 maximal cells, while 889 squares of the previous task's family (arbitrary angles) give 10,533
and 2,393 squares of the support give 126,649.  Adding arbitrary-angle squares to A1 costs about
**150 cells per square** (A1 + 80 arbitrary squares from the `t4L` family, D4-closed to 2,899 squares:
14,757 cells; + 240: 36,669 cells).  So the plan is a structured bulk plus a few hundred arbitrary
"killer" squares — not the support wholesale.

## 4. Step 5: the family loop on exact cells  [measured]

`search/unavoid13no_loop.py`.  Per round: maximal cells of `F` (C helper, `tol = 1e-11`, D4-closed) →
13-column hitting set by the B&B dive (§5) with HiGHS as fallback → polish (max-margin LP over `F`) →
float violation search (`0.01 × 0.5°` grid, Nelder–Mead, distinct local minima `≥ 0.05` apart) → the
worst `add` poses with depth `≤ −min_viol` rationalised (`den 1000`) and added with their D4 images.

| run | bulk | `add` / `min_viol` | rounds done | `|F|` | cells | LP | IP time (finder) | worst violation per round |
|---|---|---|---|---|---|---|---|---|
| L1 | A1 (grid 1/10) ∪ support:40 (snap 1/20) | 24 / 0.003 | 5 | 2,731 → 3,347 | 6,193 → 40,749 | 12.057 → 12.114 | 43 → 842 s | −0.049, −0.045, −0.037, −0.033, −0.057 |
| L2 → L2b → L2e | support:80 | 12 → 10 / 0.003 → 0.005 | 12 (+1 unresolved) | 681 → 1,617 | 5,305 → 33,353 | 12.197 → 12.206 | 2 s (dive) → 85–272 s → 716 s (portfolio) → **> 3,600 s, no set** | −0.189, −0.096, −0.093, −0.076, −0.055, −0.040, −0.050, −0.033, −0.026, −0.040 |
| L3 → L3b → L3d | support:120 | 8 → 10 / 0.01 → 0.005 | 11+ | 1,001 → 1,761 | 15,957 → 34k | 12.249 → 12.249 | 22 s (dive) → 87 s → 127–2,495 s (portfolio) | −0.175, −0.088, −0.126, −0.061, −0.054, −0.021, −0.028, −0.033, −0.027, −0.016 |

(L1 was stopped at round 5: the A1 bulk buys no LP and its cells make every IP slow.  L2/L3 were
restarted from their own families (`--resume`) three times: with the short B&B budget and the feasibility
form (L2c/L3b), then with the 4-process portfolio finder (L2d/L3c — these had the fork bug of §5 and did
one round), then with the fixed portfolio (L2e/L3d), which is what is running at hand-off.  Several
in-between single HiGHS calls ran `> 45 min` without a verdict and were discarded.)

**What the numbers say [measured].**
* `h(F) = 13` at every round of every run — **no family reached 14**.
* Cells grow by `≈ 60` per added killer square (`|F| ≈ 1,000 → 12k cells, 1,500 → 30k`), and HiGHS's time
  to *find* a 13-set grows from seconds (`≤ 6k` cells) to minutes (`15–20k`) to `> 20` minutes (`27–30k`):
  the loop's per-round cost is dominated by the finder, exactly as in the previous task, only with the
  wall moved from `10k` cells to `30k` by the C cells and the feasibility form.
* The worst violation of the polished 13-set falls from `−0.19` (round 0) through `−0.05` (rounds 4–7)
  to `−0.026` (L2b round 3) — i.e. to the level of the previous task's best maximin margin (`−0.028`).
  From there on every new 13-set is one of the "hardest" ones and each round removes only its `0.03`-
  neighbourhood in a continuum of such sets, so the number of rounds still needed is not visible from
  the trend.
* The LP barely moves (`12.20 → 12.21`, `12.25 → 12.25`): killers are chosen to cut integer sets, not
  fractional ones, and the LP is pinned by the support's measure.  The integrality gap the prover must
  close is therefore `12.2 → 14`, i.e. `1.8`; at `T = 3` HiGHS closed `5.53 → 7` (`1.47`) at the root.

## 5. Step 4: the purpose-built branch-and-bound  [measured; rules proved]

`search/unavoid13no_bb.py` (class `BB`, also used inside the loop).  Node = chosen cells; branch on the
*unhit row with the fewest cells* (every hitting set contains one of its cells), children ordered by the
node LP's values; bound = LP relaxation of the residual cover over **all** cells (HiGHS simplex,
warm-started by toggling row bounds; no presolve), pruned when `depth + ⌈LP − 1e-6⌉ > k`; plus a free
greedy packing bound, reduced-cost pruning of children (`⌈LP + rc(c)⌉ > k − depth`), and residual
dominance among children (drop `c` if its unhit-row incidence is contained in a kept sibling's).

**D4 root rule [proved].**  `F` is D4-closed and the cell set is made exactly D4-closed by adding any
missing images (`d4_close_cells`, matched through incidence sets, not coordinates — the float cells of
these families were already closed: 0 images added in every instance).  Let `Q0` be the axis-parallel
square centred at `(2,2)` (D4-invariant, in every family).  A hitting set `X` contains a cell `c ∈ Q0`;
the cells of `Q0` form a D4-invariant set; pick `g` with `gc` the orbit representative; `gX` is a hitting
set of the same size (`inc(gc', gQ) = inc(c', Q)`).  So at the root it suffices to branch on one cell per
orbit of `Q0`'s cells (`537` of `537 × ~8` for L1 round 0), and infeasibility with the rule implies
infeasibility without it.  In the HiGHS runs the same fact is the single row `Σ_{reps} x ≥ 1`.  Nothing
about symmetry is assumed deeper in the tree (a hitting set has no symmetry in general — the previous
note proved that a 13-set has none).

| instance | verdict | B&B nodes / LPs / time | HiGHS (same instance, same threads) |
|---|---|---|---|
| `T = 3` final family, `k = 6` (LP 5.53, `h = 7`) | INFEASIBLE | 606 / 606 / **85 s** (max depth 2) | 27 nodes, 49 s |
| L1 round 0 (2,731 sq., 6,193 cells), `k = 13` | FEASIBLE (dive) | 15 / 14 / **5 s** | root, 43 s |
| L1 round 0, `k = 12` (LP 12.06) | INFEASIBLE | 1 / 1 / 3 s (root LP) | 7 s |
| L2 round 4 (993 sq., 12,153 cells), `k = 13` | FEASIBLE | 333 / 332 / 130 s (712 / 270 s before the rc/dominance pruning) | root, 145 s (round 5, 14.8k cells) |
| L2 round 5 (14,797 cells), `k = 13` | undecided in 600 s (1,094 nodes) | | root, 145 s |
| L2b rounds 0–2 (20–23k cells), `k = 13` | undecided in 60 s | | root, 59–272 s |

So: as a *prover* the B&B is within 2× of HiGHS on a genuine infeasibility (`T = 3`, gap `5.53 → 7`),
which makes it an acceptable independent second solver; as a *finder* it is better than HiGHS on
`≤ 6k`-cell instances and worse above `~12k` cells, where HiGHS's root heuristics find the 13-set after
its cut loop while the LP-guided dive wanders (LP solutions of these instances are highly fractional —
`x ≈ 1/2` on many cells around the tile lines — so the "LP-value-first" child order carries little
information).  The loop therefore runs the dive with a 60 s budget and falls back to HiGHS.
HiGHS's pure feasibility form (no objective) and the min-`Σx`-with-cutoff form: `106` vs `481–510` s on
L2b round 0 with 1 thread (the D4 row and `mip_detect_symmetry` change nothing; `mip_heuristic_effort
0.3` does not help).

**The HiGHS stalls are variance, not hardness [measured] — and the fix is a portfolio.**  On L3b's
round-0 instance (1,337 squares, 26,993 cells) the loop's single HiGHS call (feasibility form, 3 threads)
ran **> 48 min without a verdict**, while the same instance with the optimisation form and `random_seed 7`
(1 thread) found a 13-set in **263 s**, and with seed 2 in 413 s.  HiGHS's MIP is effectively
single-threaded here (`~100 %` of one core at `threads = 3`), so the right finder is `solve_ip_race`
(`unavoid13no_lib.py`): four 1-thread HiGHS processes (feasibility/optimisation forms × seeds), first
verdict wins, the rest are killed — a verdict from any one run is a verdict for the instance.  The loops
L3c/L2d run with `--race 4`: `≈ 8` min per round at `27–30k` cells instead of `> 45` min.

## 6. Adversarial check of `notes/unavoid13.md` §4.2 and §2.2  [checked, no flaw found]

**§4.2 (symmetric exclusions).**  Replacement argument re-derived: for a `G`-symmetric hitting set `P`
and `p ∈ P` off the fixed set of `G`, a dominating arrangement vertex `v` has `|orbit(v)| ≤ |orbit(p)| =
|G|` and, `F` being `G`-closed, `S(gv) = g·S(v) ⊇ g·S(p) = S(gp)`, so `orbit(v)` replaces `orbit(p)`
at no extra cost; for `p` on an axis `ℓ` (the only way to have a smaller orbit is a nontrivial stabiliser,
i.e. a reflection axis or the centre), `P(p) ∩ ℓ` is a segment whose endpoints lie on square edges, hence
are edge–axis intersections, cost 1 — exactly the candidates (a) adds; (b) "no dominance pruning" is
necessary and present.  The orbit sizes in the float loop are rounded at `1e-6` (an off-axis vertex
`1e-8` from the axis counts as on-axis, i.e. *cheaper* — the safe direction for an infeasibility claim)
and exact in `unavoid13_symcheck.py`.  Lenient incidence is likewise the safe direction.  The theorem
for `C2, C4, V, D4` rests on the exact symcheck with two solvers.  **No gap found.**  `Rx` is still
running (`runs/unavoid13_symRx`, round 18 at 1,393 squares, 272k orbit candidates, ~5,000 s/round
at hand-off time); `Rd` never started (it was queued behind `Rx`), so the "no symmetry at all" theorem
did **not** land during this task.

**§2.2 (SEG lemma).**  Statement and proof re-derived (if `p ∉ Q` with `|x'_p| ≤ ½`, `y'_p ≤ ½` then
`y'_p < −½`, so `y'_q < −½ + |q − p| ≤ ½`, and with `y'_q ≥ −½`, `|x'_q| ≤ ½`, `q ∈ Q`; closed squares
throughout).  Box exactness: `max_{c, θ} ⟨p − c, e(θ)⟩ = max over the 4 rectangle corners of
max over θ`, and for fixed `v` the max of `|v| cos(θ − φ)` over a bin of width `< π` is at an endpoint
unless `φ` is inside (then `|v|`) — the cone test in `_le_half` is two cross products with the bin's
endpoint directions; the `≥ −½` conditions are the same test on `−v`.  Implementation checked against
the statement (`_seg_ok`: both points in the slab, `free(p) ≤ ½`, `free(q) ≥ −½`, four variants,
float pre-test then exact).  **Correct.**  One remark for the record: the lemma needs unit weights on
both points (`W[i], W[j] ≥ 1`), which the code enforces, so it is a *pure*-set primitive and would not
transfer to weighted covers without a weighted restatement.

## 7. "Cannot reach" — the sequence, and what would be needed  [measured / heuristic]

*(final numbers below are as of the last read of the detached loops; they keep running in
`runs/unavoid13no_L2c`, `runs/unavoid13no_L3b`; their logs are `stdout.log` and `round_log.txt`.)*

**The sequence of 13-sets** (polished; worst violation against the continuum in parentheses) is in
`runs/unavoid13no_L*/round_log.txt`, one line per round, and `P_round*.npy`; `search/unavoid13no_skeleton.py L2 L2b L2c L3 L3b`
prints their column/row skeletons.  Representative members:

| where | shape | points | worst violation |
|---|---|---|---|
| L2 round 0 | rows 4+3+4 + 2 free | `(0.89,{1,2,3})…` | −0.189 |
| L2 round 5 | columns 3+3+2+3 + 2 free at x≈2.2 | `(0.82,1) (0.75,2) (0.77,3) (1.25,3) (1.25,1) (1.29,2) (2.07,1) (2.02,3) (3.10,1) (3.13,3) (3.14,2) (2.20,1.59) (2.16,2.49)` | −0.055 |
| L2b round 3 | columns 3+4+3+3 (x≈0.93, 1.78, 2.6, 3.1) | `(0.95,1) (0.90,2) (0.95,3) (2.54,2) (2.65,1) (2.62,3) (3.10,1) (1.77,0.94) (1.77,3.06) (1.79,1.51) (3.09,2) (3.05,3) (1.78,2.46)` | −0.026 |
| skeleton IP, ≤ 9 near-line | columns 3+3+4+3 (x≈0.87, 1.41, 2.19, 3.18), free column at y = 0.74, 1.52, 2.50, 3.32 | | (hits F) |
| L3d round 0 (1,385 sq.) | columns 3+3+4+3 (x≈0.94, 1.40, 2.24, 3.06) | `(0.91,1) (0.90,2) (0.96,3) (1.49,3) (1.33,1) (1.38,2) (2.32,3) (3.05,1) (3.10,2) (3.05,3) (2.17,1.56) (2.24,2.42) (2.20,0.98)` | −0.021 |
| L3d round 1 (1,449 sq.) | | | −0.028 |
| L3d round 2 (1,529 sq.) | | | −0.033 |
| L2e round 0 (1,537 sq.) | columns 3+3+4+3 (x≈0.90, 1.39, 2.26, 3.10) | | −0.040 |

**First finder failure [measured].**  On L2e's round-1 family (`1,617` squares, `33,353` cells, LP
`12.206`) the 4-process portfolio (feasibility/optimisation × seeds 1–4, 1 thread each) found **no
13-set in 3,600 s** — the first instance of the day on which every finder failed; the loop restarted it
with a 7,200 s limit, and `search/unavoid13no_bb.py` (D4 root rule, LP bounds, 2 threads, no time limit)
was started on the same family as a prover (`runs/unavoid13no_L2e/bb_round01.log`; `Q0` has 8,405 cells,
`< 200` nodes in its first 15 min).  The B&B was stopped after 30 min below 200 nodes (each node LP has 33k columns; the
prover of §5 does not scale to this size either).  At hand-off the loop's 7,200 s race on this family
was still running; if an optimisation-form HiGHS run ends `INFEASIBLE`, the family `runs/unavoid13no_L2e/F_round01.txt` is the candidate certificate
and `python3 search/unavoid13no_recheck.py runs/unavoid13no_L2e/F_round01.txt --k 13` is the exact
re-solve (expected cost: exact vertices of 1,617 squares in Fractions, ~1 h; then three solvers).

**One mode, up to D4 [measured].**  From L2b round 3 on, *every* 13-set found by either lineage
(L2b–L2e, L3b–L3d: 8 consecutive sets) is the same configuration up to a symmetry of the container:
three columns of three points on the tile lines `y = 1, 2, 3` at `x ≈ 0.90, 1.40, 3.10` and one column
of four at `x ≈ 2.25` with `y ≈ 0.95, 1.52, 2.45, 3.02` — or its mirror image `x → 4 − x`
(columns at `0.90, 2.60, 3.10`, four-column at `1.75`), or the transpose.  The loops have converged onto
this single mode and each round removes a `0.02–0.04`-neighbourhood of one of its positions; the worst
violations in the last five rounds are `−0.021, −0.028, −0.033, −0.027` (L3d) and `−0.040` (L2e).  A
targeted attack on this mode was tried (`search/unavoid13no_mode.py`: random mode parameters → max-margin
LP polish against all of `F` → harvest violated squares if the polished set hits `F`): in 17 min on each
of the L2e and L3d families **not one** of the random mode sets could be polished into an `F`-hitting
set [measured] — the `F`-hitting positions of the mode are already too sparse to be sampled, which is
why only the exact IP finds them, and why each IP round now costs 5–40 min.

**Structural lemma attempted and refuted [measured].**  "Every 13-set hitting `F` has `≥ 10` points within
`0.12` of the tile lines `x, y ∈ {1,2,3}`" is *false* for the L2 round-6 family: the feasibility IP with the
row `Σ_{near-line cells} x ≤ 9` is feasible (44 s), the set above.  So the rigid part of a 13-set is at
most the three `{y = 1, 2, 3}` columns (9 points) and even those slide (`x ≈ 0.75–0.95`, `1.25–1.45` or
`2.5–2.7`, `3.1–3.25`).

**Killer angles [measured].**  The violated poses found by the loops lie at `20°–70°` (peaks near `28°`,
`52°`, `60°`; complementary pairs from the D4 closure), never near-axis: the near-axis "wall" squares of
the previous note matter for certifying an upper bound, not for cutting 13-sets.

**Budget arithmetic [heuristic].**  Killers per round `64–96` (8–12 poses × 8 images), `≈ 60` cells per
killer, and HiGHS time roughly `∝ cells^2` beyond `15k`: from `30k` cells (`1,500` squares) each further
round costs `≥ 30` min and the family can grow by at most `~15` rounds a day.  If the number of "hardest"
13-set modes to be cut is in the hundreds (the T = 3 loop needed 200 killers *in total*; here the last
10 rounds spent 800 killers on sets with violation `0.03–0.05`), the loop as built needs weeks, and the
final proof (HiGHS on `50–100k` cells with gap `1.8`) is itself an open question — at `T = 3` HiGHS
closed a gap of `1.47` at the root, which is the one piece of evidence that it might.

**What would change the picture.**  (i) *Row-pruned lazy IP*: keep the IP family to the bulk plus the
killers violated by the last few 13-sets (`≤ 800` rows, `≤ 12k` cells, minutes per round), check every
13-set against all of `F`, re-activate what it misses; the union `F` keeps growing but the IP does not.
(ii) *Geometric branching with column generation* (§5): branch "the point hitting `Q` is in the left /
right half of `Q`", bound by the LP over all cells via pricing on the dual's `≈ 300` support squares —
no `30k`-column LP per node.  (iii) *A better killer choice*: the deepest violation per 13-set kills a
`0.03`-neighbourhood; diverse shallow ones kill little; `min_viol` should rise with the rounds.
None of these was reached within the day.

**Stopped 2026-09-22, 22:35 (the question is parked; `TODO`).**  Last states, none with a verdict:
`L2e` round 4 solved (13-set, LP `12.2095`, `1,762 s`; worst violation `-0.024`), stopped in round 5 at `|F| = 1,857`;
`L3d` round 6 solved (13-set, LP `12.2491`, `3,156 s`; worst `-0.032`), stopped in round 7 at `|F| = 1,921`;
`symRx` (reflection `x -> 4 - x`, `notes/unavoid13.md`) round 19 at `|F| = 1,465` **feasible** (a symmetric 13-set, 8 orbits,
IP `18,566 s`; worst violation `-0.076`), stopped in round 20 at `|F| = 1,545`.  `symRd` never ran past round 1.
Families and logs: `runs/unavoid13no_L2e`, `runs/unavoid13no_L3d`, `runs/unavoid13_symRx`.

## 8. What in the brief or in `notes/unavoid13.md` turned out wrong

1. **The brief's step-1 test is not a completeness test** (§1): `LP ≥ 12.2688` holds for *any* candidate
   subset, by duality with the certified measure; it can only catch over-lenient incidence.  The brief has
   the direction backwards ("if your LP is below that, the candidate set is incomplete").
2. **The brief's "with `F ⊇` the full 2,352-square support"**: the support gives 2,348 distinct squares
   (corner-tile and on-axis poses have small orbits), and including it wholesale costs `127k` cells — it
   is the wrong `F` for the IP; the 80 heaviest poses carry `12.20` of the `12.27` at `5k` cells (§2).
3. **`notes/unavoid13.md`'s lenient tolerance `1e-7`** is unsound-in-the-useful-direction once the support
   is in `F`: sound for lower bounds, but it hands the columns squares they do not contain (coverage
   `1.54` at `(1,1)` against a true `1.000`) and costs `0.46` of LP (§1).  The previous note's LP values
   `11.3–11.6` and its "the LP over every family only 11.3–11.6 although `ν_f ≥ 12.27`" puzzle are
   partly this artefact, partly the missing support.
4. **The brief's "structure first" stages**: near-axis (`< 5°`) squares are the *worst* thing to add
   (`63k` cells for nothing, §3); the `45°` class is cheap and raises `h` from 9 to 11; the rest of the
   way to 13 needs the support's general angles.  The "9 near-lattice points forced, 4 free" picture is
   right in spirit but the 9 are not at lattice points — three sliding columns/rows of three on the tile
   lines — and no rigid skeleton is forced (§7).
5. **The brief's step 4 expectation** ("branch on the square with the fewest candidate cells … LP bound at
   each node"): built and calibrated (§5).  It is fine as a prover and poor as a finder because these
   instances' LP solutions are half-integral on the tile lines and carry no branching information;
   HiGHS's root heuristics do better.  "If HiGHS … proves infeasibility … acceptable proof provided the
   instance is dumped and re-solved by a second method": prepared (`unavoid13no_recheck.py`, three
   solvers, validated on `T = 3`) but never needed.
6. **`notes/unavoid13.md` §4.2 / §2.2**: checked line by line (§6) — nothing wrong found.  The hand-off's
   `Rx`/`Rd` promise did not materialise: `Rx` is at round 19 (`~5,000 s`/round, 305k orbit candidates)
   and `Rd` was queued behind it and never started.
7. **`notes/unavoid13.md` §4.3's lazy-row loop (`t4L`)** is stuck: its round-1 inner step 10 hit a
   4,852 s time limit without an incumbent on 11.8k candidates and doubled to 9,600 s — the same wall,
   at a third of the cell count, because of the `1e-7` tolerance and the optimisation form.
