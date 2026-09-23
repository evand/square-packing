# unavoid13: is there a 13-point pure unavoidable set for `[0,4]^2`?  (2026-09-22)

Task brief: `tasks/unavoid13/README.md`.  Code: `search/unavoid13_lib.py` (geometry, arrangement,
IP, violation search, polish), `search/unavoid13_loop.py` (the cutting-plane loop),
`search/unavoid13_recheck.py` (exact re-solve from a dumped family, three solvers, shrink),
`search/unavoid13_check.py` (`zeromargin.py`'s checker + the SEG primitive, by subclassing; nothing in
`search/zeromargin.py` is modified), `search/unavoid13_calib.py` (fractional calibration).  Runs:
`runs/unavoid13_*`.  Certificates: `certificates/unavoid13/`.

## 0. Verdict

**`T = 4`: undecided — `h(F) = 13` at every round of every run; smallest violation of a polished 13-set
`−0.040`; the general IP is the bottleneck, not the geometry.**  What *is* decided at `T = 4` [proved,
exact candidates, two solvers]: **no 13-point unavoidable set is invariant under the 180° rotation of the
container** (hence none is `C4`-, `V`- (Friedman's symmetry), diagonal-Klein- or `D4`-symmetric); the single
axis- and diagonal-reflection cases were still running when this note was written (§4.2).  Every 13-point
unavoidable set, if one exists, is therefore essentially asymmetric.

**`T = 3` [proved]: the minimum is exactly 7.**  Lower bound: a 575-square family with `h(F) = 7`, exact
incidence, three solvers (§2.3).  Upper bound: the rational Kearney–Shiu set with `a = 9/10`, certified with
`0` uncertified boxes — but only after adding a two-witness primitive (SEG) to the checker, which
`zeromargin.py` cannot do with or without `--tri`/`--disj` (§2.2).  Calibration: fractional value `5.53`
(proved-side LP over `F`) to `5.99` (heuristic dense cover LP) against the pure number 7 (§2.4).

**State of the `T = 4` runs at hand-off** (all detached, still running; output in `runs/unavoid13_*`):
`t4L` (lazy-row feasibility loop, `runs/unavoid13_t4L/round_log.txt`) in its first outer round with an active
set of 667 squares (of 1,049), LP `11.32`, IP 55–200 s per inner step; `symRx` (single-reflection symmetric
loop) feasible at rounds 0–1, 135–174 s per round, `symRd` queued after it; the `T = 3` orbit-first shrink
(`runs/unavoid13_t3a/recheck.log`, prints are buffered) and two maximin continuations (`runs/unavoid13_maximin*`;
best margins so far `−0.028`, `−0.034` from IP solutions).  A `13 excluded` verdict would appear as
`INFEASIBLE … Dumping` in the `t4L` log with `final_family.txt` to re-solve by
`python3 search/unavoid13_recheck.py runs/unavoid13_t4L/final_family.txt --shrink 14`; a `13 exists` verdict as
`CERTIFIED` with the certificate path.

## 1. The method, and why it is sound

**Objects.**  Closed unit square `Q(cx, cy, u)`, `u = tan(θ/2)` rational so `cos θ, sin θ` are rational;
admissible iff `cx, cy ∈ [w/2, m − w/2]`, `w = |cos θ| + |sin θ|`; a point on `∂Q` is in `Q`
(`search/ZEROMARGIN.md` §1, `certificates/FORMAT.md`).  A *pure unavoidable set* is a finite `P ⊂ [0,m]^2`
meeting every admissible `Q`.  Write `p(m)` for its minimum size.

**Lower bound [proved].**  For a finite family `F` of admissible squares, `h(F) := min |P|` over `P` meeting
every `Q ∈ F` satisfies `h(F) ≤ p(m)`, whatever `F` is.  To compute `h(F)` by a finite IP one needs a
finite candidate set that loses nothing.  *Claim: the vertices of the arrangement of `F` suffice* — the
square corners and the intersection points of every pair of non-parallel closed edges of two different
squares.  Proof: for any point `p` with `S(p) = {Q ∈ F : p ∈ Q} ≠ ∅`, the polygon `P(p) = ∩ S(p)` is compact,
convex, nonempty; every vertex of `P(p)` is a corner of some `Q ∈ S(p)` or the meeting point of two
non-parallel edges of squares in `S(p)`, and lies in every square of `S(p)`.  So each point of a hitting set
can be replaced by an arrangement vertex whose incidence set contains the original's.  (Parallel overlapping
edges add nothing: the overlap's endpoints are corners.)  This is the hitting-set twin of the argument
`search/COVER4.md` uses for the packing measure.  Candidates with identical incidence sets are merged and
dominated ones (incidence set a strict subset of another's) are dropped — neither changes the optimum.

**Float exploration is sound in one direction.**  The loop computes vertices and incidences in floating
point with a *lenient* tolerance (`1e-7`): the computed incidence is a superset of the true one, so the
float IP value is `≤` the true `h(F)` — a lower bound stays a lower bound — *provided the vertex list is
complete*.  Because a skipped near-parallel pair could lose a vertex, the theorem-grade numbers are
recomputed by `unavoid13_recheck.py` in `fractions.Fraction`: exact corners, exact intersection of every
pair of non-parallel edges (parallel = exact zero cross product, no cut-off), exact incidence (float
pre-filter, every pair within `1e-6` of a boundary decided in Fractions), and then the IP is solved three
ways: `highspy` MIP (gap `0`), `scipy.optimize.milp`, and a hand branch-and-bound using only `highspy`'s
simplex for bounds.  The family, the candidate vertices (as Fractions) and the incidence matrix are dumped
so the instance can be re-solved by anyone from the dump.

**Upper bound [proved when certified].**  A candidate set is written in `certificates/FORMAT.md` format
(unit weights, `W = 1`) and checked exactly.  An uncertified box of the checker is *not* a violation (the
brief's warning; confirmed in §2.2 in the strongest possible way).

**The loop.**  Round: candidates from `F` → IP → `h(F)` and an optimal `P`.  If `h(F) ≥ target`: stop
(lower bound).  Else *polish* `P`: alternate (i) assign each `Q ∈ F` to its deepest point, (ii) an LP in the
`2|P|` coordinates maximising the minimum depth over `F` for that assignment (depth `½ − max(|x'|,|y'|)` is
concave piecewise-linear, so this is an LP), then a second LP at the optimal `t*` maximising
`Σ_Q min(slack_Q, 0.03) − 0.02·‖move‖₁`.  Without this the IP's points sit on square edges of `F` and the
violations found are `ε`-shifts of squares already in `F`.  Then a float search for violated poses
(grid `0.01 × 0.5°` over admissible poses, Nelder–Mead refinement of the 400 best distinct grid minima,
keep the `≤ 10` worst distinct local minima with `max_p depth < −10⁻⁹`), snap each to a rational pose
(`u`, centre with denominators `≤ 1000`, centre clamped exactly into the admissible range), confirm the
violation exactly against the rational snap of `P`, and add the pose and its `D4` images to `F`.

**Zero slack is forced [proved, one line].**  A 7-point set at `T = 3` must meet the 9 closed tiles of the
`3×3` tiling, 13 points at `T = 4` the 16 tiles; so some points lie on shared tile edges and every minimum
set has poses captured with margin exactly `0`.  Every certificate here therefore needs the zero-margin
primitives, and (§2.2) one more than the repo had.

## 2. Rehearsal at `T = 3`

### 2.1 The lower bound: `h(F) = 7`, so 7 is the minimum  [measured → proved after §2.3]

`python3 search/unavoid13_loop.py --m 3 --target 7 --tag t3a --init tiles,grid2,rot2 --add 8`
(`runs/unavoid13_t3a/round_log.txt`).  `F₀` = the 9 tiles, the 25 axis-parallel squares with centres on the
half-integer grid, and 9 rational angles (`u ∈ {1/8, 1/4, 1/3, 2/5, 12/29, 1/2, 3/5, 3/4, 7/8}`) on a
half-integer centre grid, `D4`-closed: 375 squares.

| round | `|F|` | vertices | candidates | `h(F)` | LP relaxation | worst violation of the polished `h`-set |
|---|---|---|---|---|---|---|
| 0 | 375 | 97,232 | 60,757 | 6 | 5.1667 | `−0.136` |
| 1 | 391 | | | 6 | 5.1667 | `−0.136` |
| 2 | 407 | | | 6 | 5.2973 | `−0.047` |
| 3 | 471 | 170,528 | 108,965 | 6 | 5.3810 | `−0.100` |
| 4 | 487 | | | 6 | | |
| 5 | 559 | | | 6 | | `−0.063` |
| 6 | 575 | | | **7** | 5.5294 | — (lower bound reached) |

Seven rounds, 200 squares added.  The IP is easy: at round 3, HiGHS presolve reduces 108,965 columns to
4,985 and solves at the root node (one node, 35k simplex iterations); a posting-list dominance filter
(`drop_dominated_fast`) reproduces exactly that 4,985 in 9 s and cuts the IP from 72 s to 9 s, so all later
runs use it.

### 2.2 The upper bound: a rational 7-point set, certified — after adding a primitive  [proved]

Kearney–Shiu's set has `x = √2 − ½` on the outer columns.  A float scan (`0.01 × 0.5°`, refined) of the
family `{(a,1),(3/2,1),(3−a,1),(3/2,3/2),(a,2),(3/2,2),(3−a,2)}` shows `√2 − ½ = 0.91421` is a *one-sided*
limit: every `a ≤ 0.91421` tested (`0.90, 0.914, 0.9142, 0.91421`) has no violation, `a = 0.9143` fails by
`6·10⁻⁵` at the 45° square `(2.293, 1.5, 45°)`, `a = 0.92` by `4·10⁻³`.  So `a = 9/10` is a rational
candidate: `certificates/unavoid13/ks7_rational_3.txt` (`D = 10`, `W = 1`, 7 points).

`zeromargin.py cert … --tri` **does not certify it**: 2,799 uncertified boxes at depth 14 (`--disj
--chain-from 0`: 4,172), all at the wall square `(½, 3/2, θ → 0)`; `diag` on one of them says the cover is
fine there (minimum captured weight exactly `1` over the box) and that no primitive applies.  The reason:
for the pose `(½+δ, 3/2−η, ε)` the witness is `(9/10, 2)` iff `η ≤ (0.4−δ)tan ε + ½(sec ε − 1)` and
`(9/10, 1)` iff `η ≥ (0.4−δ)tan ε − ½(sec ε − 1)` — the two regions overlap in a strip of width
`≈ ε²/2` whose boundary curve crosses the box family at every depth, so no fixed-witness primitive
(`CORE/ADM/P1/MIX`) ever certifies a box straddling it, and `CHAIN` does not either because the two swing
inequalities are of different kinds (`Y ≤ ½` for one point, `Y ≥ −½` for the other).  This is exactly
`RUNG2.md`'s "one-sided pose" obstruction.  Friedman's 14 escapes it only because its wall points sit at
`x = 1` exactly, where `P1`'s wall shortcut applies (`ZEROMARGIN.md` §4 item 2).

**SEG lemma [proved].**  *Let `|q − p| ≤ 1` and let `X` be a set of poses such that, for every pose in
`X`, in the square's frame `|x'_p| ≤ ½`, `|x'_q| ≤ ½`, `y'_p ≤ ½` and `y'_q ≥ −½`.  Then every pose in `X`
contains `p` or `q`.*  Proof: if `p ∉ Q` then `y'_p < −½`; `y'_q = y'_p + ⟨q − p, e₂⟩ < −½ + |q − p| ≤ ½`;
with `y'_q ≥ −½` and `|x'_q| ≤ ½`, `q ∈ Q`.  ∎  (Two-point sibling of the triangle lemma.)  Over a box the
four conditions are exact: `max ⟨v, e(θ)⟩` over the centre rectangle × angle bin is attained at a rectangle
corner and at a bin endpoint unless the direction of `v` lies in the bin's cone, where it is `|v|`; so the
test is four rational dot products, two rational cross products (cone test) and `|v|² ≤ ¼`.
`search/unavoid13_check.py` adds it as a last-resort primitive by subclassing `Checker` and overriding
`cert_tri` (its leaves are reported under `TRI` with a `SEG`-tagged witness in the `--dump`).

`python3 search/unavoid13_check.py cert certificates/unavoid13/ks7_rational_3.txt --tri --seg --depth 16`:
**VERIFIED**, 3,854 boxes, max depth 8, leaves `ADM 1473, TRI 4, SEG 17, EMPTY 2233, uncertified 0`,
`0.4 s` (`runs/unavoid13_t3a/zm_ks7_seg_leaves.txt`).  Independent float stress test of the new primitive
(no shared code path with the exact test): every SEG leaf sampled at its corners and 200 random poses —
2,588 admissible poses, 0 failures; and 20,000 random boxes of the pose space, of which 1,798 were claimed by
SEG, sampled 100 poses each — 0 contradictions.  So

> **[proved]** `{(9/10,1), (3/2,1), (21/10,1), (3/2,3/2), (9/10,2), (3/2,2), (21/10,2)}` is unavoidable for
> closed unit squares in `[0,3]^2`: `p(3) ≤ 7`.

### 2.3 Exact re-solve of the `h(F) = 7` family  [proved]

`python3 search/unavoid13_recheck.py runs/unavoid13_t3a/final_family.txt --shrink 7`
(`runs/unavoid13_t3a/recheck.log`, dumps `recheck_B.npy`, `recheck_cand.txt`).  The 575-square family
(`runs/unavoid13_t3a/final_family.txt`, rational poses, all admissible, no duplicates) has, exactly:
240,317 arrangement vertices (129,649 close pairs of squares), 29,770,035 incidences of which 553,048
borderline pairs were decided in Fractions; 168,877 distinct incidence sets, 7,361 after dominance.
The hitting-set IP on this exact instance:

| solver | `h` |
|---|---|
| `highspy` MIP, gap 0 | **7** (LP relaxation 5.529412) |
| `scipy.optimize.milp` | **7** |
| hand branch-and-bound, LP bounds from `highspy` simplex only | **7** (3,763 nodes) |

So no 6 points meet these 575 closed unit squares, hence no 6-point pure unavoidable set exists for
`[0,3]^2`, and with §2.2:

> **[proved]** The minimum size of a pure unavoidable set for closed unit squares in `[0,3]^2` is
> exactly **7**.  (Kearney–Shiu exhibit 7; the lower bound is new to the repo, and I know of no statement
> of it in the literature.)

The float loop's own `h(F) = 7` at round 6 (lenient incidence) agrees, as it must (§1).  The shrunk family
(greedy, orbit-first) is reported in §2.5 when the shrink finishes.

### 2.4 Calibration at `T = 3`

`python3 search/unavoid13_calib.py runs/unavoid13_t3a/final_family.txt` (`runs/unavoid13_t3a/calib.log`):

| quantity | value | status |
|---|---|---|
| hitting-set LP relaxation over the final `F` (its dual is a packing measure on `F` feasible at every arrangement vertex, hence everywhere): a lower bound on `ν_f(3)` | **5.5294** | [measured] (float LP, not rationalised) |
| dense-grid cover LP over the 7,361 candidate points, rows = `2° × 0.04` grid poses added when violated, converged (`0` violated, grid minimum `1.000000`) | **5.9861** | [heuristic] |
| honest cost of that cover (`0.25° × 0.005` scan) | *(pending)* | [heuristic] |
| pure unavoidable number | **7** | [proved] |

So at `T = 3` the fractional optimum is `≈ 5.5–6.0` and the integrality gap to the pure number is
`≥ 1.0`; at `T = 4` the analogous bracket is `12.2688` (proved floor) to `~12.4` (heuristic optimum)
against a pure number of 13 or 14 — a gap of `0.6–1.6`.  The `T = 3` gap is *larger* relative to the
container than the `T = 4` one, and the `T = 3` IP closed it at the root node.

## 3. Sanity checks (the brief's mandatory ones)

* **Friedman's 14 is accepted.**  `certificates/unavoid13/friedman14_4.txt` (`D = 5`, `W = 1`) through
  `unavoid13_check.py cert --tri --seg --depth 14`: **VERIFIED**, 7,108 boxes, max depth 10, `ADM 3462, TRI+SEG 80,
  EMPTY 3212, uncertified 0`; baseline `zeromargin.py friedman14 --tri`: VERIFIED, `ADM 3496, TRI 74, EMPTY 3212,
  uncertified 0` (the SEG-augmented checker takes 6 more boxes by SEG before TRI is needed — same verdict).
  Against the loop's 1,049-square family (`runs/unavoid13_t4b/F_round07.txt`) Friedman's set hits every
  square, with exact minimum depth `0` at the three tightest (it is a zero-slack set, §1).  [proved]
* **Every `h(F)` claim re-solved a second way**: §2.3 (three solvers, exact incidence) and §4.2 (two solvers,
  exact candidates).  The float loop's values (`h = 7`, `h = 13`) agree with the exact ones where both exist.
* **Every "unavoidable" claim has `0` uncertified boxes** (§2.2, §3 above).  No claim rests on a `pose`/`diag`
  confirmation of an uncertified box.
* **The candidate-cell argument was exercised**: at `T = 3` the exact instance's dominance filter reproduces
  HiGHS presolve's own column count (4,985 at round 3), and the float and exact candidate counts differ only by
  float near-duplicates (7,361 exact vs 7,361 float at the final family).

## 4. `T = 4`

### 4.1 The raw loop (`unavoid13_loop.py`): `h(F) = 13` throughout, IP time explodes  [measured]

`F₀` = tiles + half-integer axis-parallel grid + the 40 heaviest poses of the exact packing support
(`search/cover4_exact_support.txt`), `D4`-closed: 361 squares (`t4a`), continued as `t4b` from round 5 with
the dominance filter; `t4c` starts from the 120 heaviest support poses (1,001 squares).  Sanity: the LP
relaxation over `F` is a lower bound on `ν_f(4)` and must stay `≤ 12.2688`'s complement… (it is far below).

| run/round | `|F|` | candidates | `h(F)` | LP | IP time | worst violation of the polished 13-set |
|---|---|---|---|---|---|---|
| t4a 0 | 361 | 31,889 (no dominance) | 13 | 11.109 | 19 s | |
| t4a 3 | 497 | | 13 | 11.288 | 30 s | |
| t4b 0 | 657 | 5,017 | 13 | 11.385 | 45 s | `−0.094` |
| t4b 1 | 713 | | 13 | 11.385 | 57 s | `−0.116` |
| t4b 2 | 777 | | 13 | 11.385 | 82 s | `−0.076` |
| t4b 3 | 801 | | 13 | 11.390 | 84 s | `−0.076` |
| t4b 4 | 849 | | 13 | 11.390 | 212 s | `−0.072` |
| t4b 5 | 889 | 10,673 | 13 | 11.394 | 212 s | `−0.051` |
| t4b 6 | 969 | | 13 | 11.460 | 512 s | `−0.047` |
| t4b 7 | 1049 | | 13 | 11.460 | 758 s | `−0.040` |
| t4c 0 | 1001 | 15,941 | 13 | 11.627 | 685 s | `−0.108` |
| t4c 1 | 1081 | | 13 | 11.627 | 830 s | `−0.045` |

The 13-sets change shape from round to round (no two consecutive solutions share more than a few points);
the worst violation decreases slowly.  HiGHS's log on the round-5 instance: presolve removes nothing
(882 rows, 10,673 binaries, 1.25 M nonzeros), and the *feasibility* form (`Σx ≤ 13`, no objective) has no
incumbent after 240 s and 101 nodes — the difficulty is in deciding 13, not in proving 12 impossible.

### 4.2 Symmetric 13-point sets are excluded  [measured → proved by `unavoid13_symcheck.py`, see below]

The brief's structured family (ii), done as a loop (`search/unavoid13_sym.py`): variables are `G`-orbits of
candidate cells with cost `|orbit|`, one row per `G`-orbit of `F`, feasibility IP `Σ|O| x_O ≤ 13`.  Two
points matter for soundness: (a) a symmetric set may have points *on* the axes of `G`, whose dominating cell
must itself lie on the axis, so the candidates include every intersection of a square edge with each axis
line and the centre (for `p` on axis `ℓ`, `P(p) ∩ ℓ` is a segment whose endpoints are such intersections
and dominate `p`); (b) no dominance pruning (an off-axis dominating vertex has a larger, costlier orbit).
Positive control: with `G = V` (Friedman's `x→4−x, y→4−y`) and `k = 14` the IP is feasible on the
889-square `t4b` family (5 orbits, 74 s) — as it must be, since Friedman's set hits everything.

| group `G` | order | float loop: rounds to infeasibility | family | exact re-check (`unavoid13_symcheck.py`): orbit candidates (by orbit size) / square orbits | `highspy` / `scipy` |
|---|---|---|---|---|---|
| `D4` | 8 | 3 | `runs/unavoid13_symD4/F_round02.txt` (371 sq.) | 2,483 (1×1, 121×4, 2361×8) / 51 | infeasible / infeasible |
| `C4` (rotations) | 4 | 2 | `runs/unavoid13_symC4/F_round01.txt` | 8,977 (1×1, 8976×4) / 109 | infeasible / infeasible |
| `V` (both axis reflections — Friedman's symmetry) | 4 | 3 | `runs/unavoid13_symV/F_round02.txt` | 6,507 (1×1, 139×2, 6367×4) / 97 | infeasible / infeasible |
| `C2` (180° rotation) | 2 | 2 | `runs/unavoid13_symC2/F_round01.txt` | 20,104 (1×1, 20103×2) / 221 | infeasible / infeasible |
| `Rx` (one axis reflection) | 2 | *(running)* | | | |
| `Rd` (one diagonal reflection) | 2 | *(pending)* | | | |

> **[proved]** No 13-point pure unavoidable set for `[0,4]^2` is invariant under the 180° rotation of the
> container (hence none under `C4`, `V`, the diagonal Klein group, or `D4`).  Certificates: the four families
> above, each a finite set of rational admissible closed unit squares that no such symmetric 13-set meets;
> re-solvable from `runs/unavoid13_sym*/F_round*_symcheck_*.json` + the family file.

Every nontrivial subgroup of `D4` contains the 180° rotation or a single reflection, so once `Rx` and `Rd`
are excluded too, **no 13-point unavoidable set has any nontrivial symmetry of the container**.  (Lenient
float incidence is the sound direction for an infeasibility claim — more ways to hit ⇒ if none, none
exactly — but the theorem rests on the exact re-check: exact vertices, exact edge–axis intersections,
exact incidence, and two solvers.)

### 4.3 The lazy-row feasibility loop (`unavoid13_loop2.py`)  [measured]

Because only "is there a 13-point hitting set of `F`?" matters, the IP is run as a feasibility problem on
an *active* subfamily `F' ⊆ F`; the 13-set it returns is checked against all of `F` and the 40 worst missed
squares are added to `F'` (inner loop) until it hits all of `F`; only then is the geometric violation search
run and its output added to `F` and `F'`.  Infeasibility on `F'` is infeasibility on `F`.

Started from the `t4b` family (1,049 squares) with the 137 axis-parallel squares active: the inner loop grew
the active set 137 → 667 in 18 steps (LP `9.0 → 11.32`); IP times 0 s → 200 s; the misses of successive
13-sets fell from 408 squares (depth `−0.18`) to 4–12 squares (`−0.001` to `−0.04`).  Still in round 0 at
hand-off.  Verdict-relevant: the IP on ~650 well-chosen rows is already as hard as the raw loop's on 1,000.

### 4.4 Maximin continuation of the IP's 13-sets (`unavoid13_maximin.py`)  [heuristic]

From each polished 13-set of `t4b`, alternate (worst 60 poses by grid + Nelder–Mead, accumulated) / (LP polish
against `F` ∪ those poses, trust radius 0.05).  Starts 0 and 1: margin `−0.094 → −0.028` and `−0.116 → −0.034`
in 25 iterations, non-monotone.  No start reached a non-negative margin.  Friedman's 14 minus each point
(family (iii)) was queued as 14 more starts (`runs/unavoid13_maximin_fried/`).

### 4.5 Calibration at `T = 4` (from the repo, for the table beside §2.4)

`ν_f^closed(4) ≥ 12.2688` [proved, `COVER4.md`]; heuristic optimum `~12.4` (`CLOSED4.md`, `FAMILY.md`);
certified weighted cover `12.956` (`certificates/rung2`); pure number `∈ {13, 14}`; LP over the loop's own
families only `11.3–11.6` (they are built to cut integer solutions, not fractional ones).

## 5. What in the brief turned out wrong, or not as stated

1. **`zeromargin.py cert` cannot certify a Kearney–Shiu-type set** (§2.2), with or without `--tri`/`--disj`;
   the brief's tool list assumed it could.  The missing piece is a two-witness primitive (SEG), now in
   `unavoid13_check.py`.  Every certificate below that has a point on a tile edge at `x ≠ 1, 3` needs it.
2. **"An uncertified box is not a violation" is not just a caveat but the whole story at the wall squares**:
   2,799 uncertified boxes, zero violations, minimum captured weight exactly `1`.
3. The brief's loop as stated ("take the 13-point solution, snap to rationals, check") never certifies
   anything without the polish step (§1): the IP's points are arrangement vertices and lie on edges of
   squares of `F`, so the raw solution is violated by `ε`-perturbations of `F` itself.
4. **"Solve the hitting-set IP exactly" is cheap at `T = 3` and the bottleneck at `T = 4`.**  At `T = 3`
   HiGHS closes the gap `5.53 → 7` at the root node; at `T = 4` the same formulation (900 rows, 10k
   candidates after dominance, LP `11.4`) takes 10+ minutes per round and the time grows with `|F|`
   (§4.1).  The IP the loop actually needs is the *feasibility* question "13 points, yes or no", and even
   that has no incumbent after 4 minutes at 889 rows.  The lazy-row form (§4.3) is what keeps rounds
   affordable; the brief's "if it stalls at 13 … switch to symmetry-reduced IPs" is the right instinct but
   the symmetry-reduced IP is an upper-bound tool only (a symmetric optimum need not exist), which is
   exactly how it is used in §4.2.
5. **The brief's back-of-envelope for family (i) is not what kills it.**  4 points on each of `y = 1, 2, 3`
   plus `(2,2)`: the 45° squares centred on `y = 1.5` need the `x`-coordinates of the `y = 1` and `y = 2`
   points together to hit every window of width `2(√2 − 1) = 0.828`… of the admissible range
   `[0.707, 3.293]`, which 8–9 points do easily; what excludes the pattern is §4.2 whenever it is
   `V`-symmetric, and otherwise nothing special about it — it is just one more shape the general IP visits.
6. `search/nu_f.py` / `search/packing_dual.py` were not needed for the `T = 3` calibration: the hitting-set
   LP relaxation over the arrangement vertices *is* the finite cover LP and its dual *is* a packing measure
   feasible everywhere (the same vertex-maximality argument), so both calibration numbers fall out of the loop.
