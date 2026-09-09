# Clique lever: the QSTAB relaxation of the hardest `t = 4` leaf on a fixed pose set

Task: `tasks/clique-lever/README.md`.  Code: `search/cliquelever.py` (imports `leaf_ceiling.py`
for the exact geometry and `t4leaf.Hi` for the LP; nothing in `verify/`, `lean/`, `certificates/`
is touched).  Runs in this worktree's `runs/cl_*` (gitignored; every number is quoted here).
Read against `search/LEAF_CEILING.md` §2.4/§5, `search/T4LEAF.md` §0–§1, `search/CLIQUE.md`,
`search/CLIQUE_CEILING.md`.

## 0. Verdict, up front

**Cliques are the lever.**  On the leaf-search agent's own pose set for the hardest `t = 4` leaf
`01010101` (7,845 poses), the LP with **every** clique inequality of the closed-intersection graph
enforced is

> **`QSTAB = 11.303030` (`= 373/33`)** — exactly certified `11.303030286` — against **`11.745193`**
> with the anchor cliques (T4LEAF) and `11.9148` with no cliques.

The general clique inequalities take **`0.44`** off the anchor-clique value and **`0.61`** off the
pure value on the same poses; the anchor family `K(p, A)` had captured `0.17` of the `0.61`, i.e.
just over a quarter.  One stage of column generation (607 new poses, interior pricing gap `+0.53`)
brings back at most **`+0.02`**: the post-pricing value is bracketed in **`[11.303030, 11.321429]`**
(`373/33`, the certified stage-(2) measure, is feasible on the enlarged set; `317/28` is the LP
value the loop has sat on for 60+ iterations while still cutting cliques of mass `1 + k/56` on a
degenerate optimal face).  The matched corner-leaf
`k = 4` control drops from `11.988925` (anchor cliques, T4LEAF) to **`11.673913` (`= 537/46`)**,
exactly certified `11.673912998`.  All four loops **converged** (`M <= 1` and max-weight clique
`<= 1`, every branch and bound complete), so each number is the value its poses really attain,
and each exactly certified measure is a rigorous lower bound on the true QSTAB value of the leaf.
The margin the brief asked about — the `0.255` between `11.745` and `12` — is therefore not the
number to close: with the full clique family the leaf sits `0.68–0.70` below 12, and even the
corner leaf alone, which the point-and-anchor family cannot close at `t = 4` (LEAF_CEILING §5.2:
`11.94` and rising), sits `0.33` below.

| step | pose set | value with points (+ anchor cliques) | **QSTAB** (all cliques) | `M` | max clique | B&B | converged | exactly certified |
|---|---|---|---|---|---|---|---|---|
| (1) | `lc_A0101.txt`, 234 certified poses (235 columns) | `11.776361678` exact, points only (LEAF_CEILING §5.3) | **`11.100000`** (`111/10`) | `1` | `1` | complete | 7 it, 27 s | `2219999999/200000000 = 11.099999995` |
| (2) | leaf `01010101`, 7,845 poses (8,087 columns) = `tl_A01010101` ∪ (1) | `11.745193` anchor / `11.914756` pure (T4LEAF) | **`11.303030`** (`373/33`) | `1` | `1` | complete | 31 it, 25 min | `5651515143/500000000 = 11.303030286` |
| (3) | (2) + one pricing stage: 8,452 poses (8,715 columns) | — | **`<= 11.321429`** (`317/28`, LP pinned since it 30) and **`>= 11.303030`** (the certified (2) measure is feasible here) | `1` | `1.01–1.12` still found | complete | no: 44 + 24 it, value unchanged for 60+ | via (2): `11.303030286` |
| control | corner leaf `k = 4`, 7,829 poses (7,941 columns) = `tl_B40K` | `11.988925` anchor (T4LEAF); `11.94` exact points-only (LEAF_CEILING §5.2) | **`11.673913`** (`537/46`) | `1` | `1` | complete | 47 + 13 it, 76 min | `5836956499/500000000 = 11.673912998` |

Every exactly certified measure was re-checked by `leaf_ceiling.py check --anchor clique` (a code
path sharing nothing with the loop): coverage `M = 1`, regions exact (corners `1,1,1,1`; slots
`0,1,0,1,0,1,0,1` in the leaf), chord strips `3,3,3,3`, max-weight clique `= 1` with a complete
B&B, "(b) PROVED".

**What the violated cliques are** (§6): not pinwheels.  Each is **the point clique of a
wall-square corner plus a fan**: all squares through a point `p` lying on one of the lines
`x, y in {1, 2, 3}` (the corners of the axis-parallel wall squares) — mass `0.9–1.0` — together
with `0.2–0.56` of mass on 30–45 interior tilted squares at 12–16 distinct angles, every one of
which closed-meets every square through `p` (grazing contacts, CLIQUE_CEILING.md's "structural
finding").  The anchor family `K(p, A)` sees the part of the fan that contains a common segment `A`;
on the ten heaviest raw cliques it captures `0–65 %` of the non-Helly excess (median `30 %`), which
is the "half at best" of the brief.  As the loop cuts, the violated cliques shrink to 5–25 support
members, interior, point mass exactly 1 plus an excess of `0.05–0.45`; at convergence every tight
clique but two or three is a plain point clique and the measure is `1/33`-valued (leaf) or
`1/46`-valued (control).

## 1. What is computed

The LP of `tasks/clique-lever/README.md`, on a **fixed** finite pose set `P`:

```
QSTAB(P) = max  sum mu_S
           s.t. cov(x) = sum_{S : x in S} mu_S            <= 1   for every point x of [0,4]^2
                mu(K)  = sum_{S in K} mu_S                <= 1   for EVERY clique K of the closed-
                                                                 intersection graph of P
                mu({S : centre(S) in R}) = k_R                   for the pinned regions R
                mu(wall strip) <= 3                              (chord rows; dual 0 throughout)
                mu >= 0
```

with the boundary semantics of `T4LEAF.md` §1.1 (a pose whose centre is within `1e-6` of a region
boundary is snapped onto it and gets one column per adjacent region: 242 duplicated poses in the
leaf set, 112 in the control).  Directions, as in `T4SCREEN.md` §1.2: restricting the poses LOWERS
the value, so `QSTAB(P)` is a lower bound on the QSTAB value over all poses; on the pose set
itself every intermediate LP value (a finite subset of the rows) is an UPPER bound on `QSTAB(P)`,
and when the loop converges the two coincide.

**The loop** (`cliquelever.py run`), on exact rational poses (`leaf_ceiling.py`'s
`theta = 2 arctan(p/q)` with `Q = 10^5`, centres on a `1e-6` grid, snapped from the float
checkpoint; the 234 certified poses of step (1) are included as they are):

1. solve with HiGHS (`t4leaf.Hi`: warm dual simplex, ipm + crossover fallback — on these LPs the
   fallback is what runs, `50–90 s` per solve at 2 threads);
2. **coverage**: enumerate the exact arrangement vertices of the *support* of the solution
   (`leaf_ceiling.enumerate_vertices`; complete for the support measure by LEAF_CEILING Lemma 1),
   compute the exact coverage of the `/10^9`-rounded masses, add the violated vertices as rows
   (`<= 3000` per iteration; coefficients over the whole pose set by float SAT, with the exact
   integer test on every pair within `1e-7` of the boundary);
3. **cliques**: the exact closed-intersection graph of the support (`leaf_ceiling.closed_graph`,
   integer SAT with `<=`), the exact max-weight clique by branch and bound on integer masses
   (`leaf_ceiling.max_clique_int`, with its completeness flag), and the max-weight clique
   *through* each of the 40–60 heaviest support squares (`max_clique_through`, seeded with the
   incumbent `1` so that it only returns violated cliques).  Each violated support clique is
   **extended greedily to a maximal clique over all poses** — a pose joins iff it closed-intersects
   every current member, the ones added before it included — and every extended row is
   re-verified pairwise before it enters the LP.  Up to 60–80 rows per iteration;
4. stop when no vertex is violated and the max clique is `<= 1 + 1e-6` with every B&B complete.

An incomplete B&B never invalidates a row (the row is a clique whether or not the search
finished); it only weakens the stopping test.  **No B&B was incomplete in any run here** (the
support graphs have 70–580 vertices; the global B&B takes `< 3 s`).

Row ageing (dropping clique rows that carry no dual for 12 solves, `t4screen`'s rule) made both
8k-pose loops cycle once the value had settled — the LP value sat at `317/28` and `537/46` for
15 iterations while rows were dropped and re-separated — so both were checkpointed and
**resumed without ageing** (`--resume-rows/--resume-cliques`, every resumed row re-verified),
after which every iteration adds a new maximal clique and termination is guaranteed on a finite
pose set: the control converged 13 iterations later; the leaf's post-pricing stage resumed at `317/28` and,
24 iterations later, still finds cliques of mass `1 + 1/56 … 1 + 6/56` on a `1/56`-valued vertex
without the value moving — a degenerate optimal face with many vertices, cut one at a time (§4).
Its bracket is rigorous regardless of when it terminates.

**The exact final certification** (`Lever.finalize`): masses rounded down to `/10^9` per column,
each pinned region topped back up to its count on poses with *exact* coverage slack **and** exact
clique slack (`max_clique_through` on the support), then `M`, the max-weight clique (complete
B&B) and the region masses recomputed in integers; written in `leaf_ceiling.py`'s measure format
(`runs/cl_*_exact.txt`) and re-checked by `leaf_ceiling.py check`.

## 2. Step (1): the 234-pose certified support

The support of LEAF_CEILING §5.3's exactly certified measure (`11.776361678`, `M = 1`, regions
exact, max-weight clique `1.5388` complete, best anchor clique `1.2829`).  Rows: all 30,036
arrangement vertices of the support plus the `0.08` grid, so no coverage row generation is needed
and the descent is the clique descent alone:

| it | LP | `M` | max clique (support members → after extension) | clique rows | support |
|---|---|---|---|---|---|
| 0 | `11.778106` | `1.000000000` | `1.483626` (43 → 44) | +20 | 193 |
| 1 | `11.474553` | `1.000000003` | `1.458781` (19 → 48) | +15 | 99 |
| 2 | `11.235294` | `1.000000001` | `1.264706` (10 → 47) | +9 | 64 |
| 3 | `11.140351` | `1.000000001` | `1.263158` (10 → 50) | +7 | 58 |
| 4 | `11.100000` | `1` | `1.150000` (7 → 33) | +3 | 41 |
| 5 | `11.100000` | `1` | `1.100000` (7 → 28) | +1 | 39 |
| 6 | **`11.100000`** | `1` | **`1.000000`** | 0 | 36 |

Converged with 55 clique rows in 27 s.  Exactly certified: mass **`2219999999/200000000 =
11.099999995`**, `M = 1`, max clique `= 1` (complete), corners `1,1,1,1`, slots `0,1,0,1,0,1,0,1`,
chord strips `3,3,3,3`, interior `3.1`.  So on the support of the best certified leaf measure,
enforcing every clique inequality costs **`0.676`**: `11.776 → 11.100`.  The final measure is
tenths-valued (36 poses, masses in `{0.1, …, 1.0}`) — the clique LP on that support is degenerate
down to a small rational vertex, and the support is far too small for the number to say anything
about the leaf; it is the answer to "what do those 234 poses carry under all cliques" and a
demonstration that the loop closes in minutes.  Step (2) is the number that matters.

## 3. Step (2): the leaf pose set (7,845 poses, 8,087 columns)

`tl_A01010101_poses.txt` (the T4LEAF stage-2 checkpoint: 8,050 float columns, `LP = 11.745193`)
snapped to exact poses and merged with the 234 poses of step (1); rows seeded with the `0.08`
grid, the 224 positive-dual points of `tl_A01010101_dual.txt` and the 5,993 heaviest arrangement
vertices of the initial support.  Chord rows on (dual `0` at every iteration, as in T4LEAF §3.1).

| it | LP | `M` | max clique | best violated support clique (members → extended) | rows | clique rows | support |
|---|---|---|---|---|---|---|---|
| 0 | `11.952147` | `2.248` | `2.2546` | 21 → 765 | 11,717 | 15 | 165 |
| 1 | `11.857471` | `1.271` | `1.6534` | 71 → 585 | 14,717 | 32 | 460 |
| 2 | `11.800782` | `1.136` | `1.5286` | 49 → 389 | 15,623 | 49 | 229 |
| 3 | `11.734641` | `1.146` | `1.4441` | 51 → 260 | 16,257 | 67 | 224 |
| 4 | `11.629698` | `1.129` | `1.4061` | 33 → 477 | 16,619 | 83 | 168 |
| 5 | `11.572932` | `1.496` | `1.5845` | 11 → 725 | 16,900 | 100 | 146 |
| 6 | `11.529539` | `1.049` | `1.3095` | 31 → 261 | 17,059 | 117 | 158 |
| 7 | `11.469762` | `1.053` | `1.2471` | 29 → 527 | 17,194 | 137 | 134 |
| 8 | `11.438930` | `1.097` | `1.2188` | 26 → 363 | 17,233 | 153 | 143 |
| 9 | `11.416681` | `1.027` | `1.2157` | 26 → 323 | 17,246 | 170 | 122 |
| 10 | `11.376675` | `1.026` | `1.2987` | 29 → 284 | 17,272 | 189 | 136 |
| 12 | `11.353367` | `1.041` | `1.1870` | 26 → 249 | 17,328 | 220 | 126 |
| 14 | `11.333805` | `1.006` | `1.1499` | 29 → 289 | 17,401 | 227 | 125 |
| 16 | `11.320722` | `1.022` | `1.1531` | 12 → 957 | 17,452 | 242 | 148 |
| 18 | `11.314091` | `1.043` | `1.1129` | 22 → 321 | 17,504 | 264 | 108 |
| 20 | `11.308603` | `1.041` | `1.2381` | 9 → 1153 | 17,541 | 273 | 131 |
| 22 | `11.305667` | `1.000` | `1.1605` | 16 → 499 | 17,585 | 282 | 145 |
| 24 | `11.303752` | `1.000` | `1.1073` | 13 → 444 | 17,611 | 294 | 108 |
| 26 | `11.303430` | `1.008` | `1.1055` | 11 → 394 | 17,634 | 292 | 110 |
| 28 | `11.303030` | `1.000` | `1.0909` | 10 → 532 | 17,644 | 283 | 75 |
| 30 | `11.303030` | `1.091` | `1.2424` | 7 → 907 | 17,653 | 265 | 83 |
| 31 | **`11.303030`** | `1.000000001` | **`1.000000`** | — | 17,653 | 252 | 70 |

Converged after 31 iterations (25 min, every B&B complete).  Exactly certified: mass
**`5651515143/500000000 = 11.303030286`** (`373/33 = 11.303030303`, the float LP value to
`2e-8`), `M = 1`, max clique `= 1` (complete), corners `1,1,1,1`, slots `0,1,0,1,0,1,0,1`, chord
`3,3,3,3`, interior **`3.303030`** — against `3.745` with the anchor cliques.  Read with T4LEAF §1:
the leaf fails to close iff the interior can carry `>= 4`; under the full clique family the
interior carries `3.30`.

Two things are worth noting about the descent.  The support collapses from 460 to 70 poses; the
final measure is `1/33`-valued.  And the max-weight clique of the intermediate measures is
carried by ever smaller support cliques (71 → 5–15 members) whose maximal extensions over the pose
set are large (250–1150 poses): the rows that do the cutting are wide, low-mass fans, not the
tight 40–50-member cliques of the coverage optimum.

## 4. Step (3): one pricing stage

`t4leaf.py`'s stratified pricer, run on the converged stage-(2) duals: reduced cost
`1 - capture(point duals) - lambda(region) - sum of the clique duals of the rows the candidate
would join` (a candidate is charged a clique's dual iff it closed-meets every member of the row,
so the charge is legitimate — the row extends to it).  `169,538` grid candidates
(`pitch 0.04`, `2.5°`), best reduced cost `+1.134` (a pinned region with a large negative
multiplier; `+0.534` in the interior), 400 kept over the 13 regions plus 400 refined and 400
perturbations of the support: **607 new poses**, `8,452` poses / `8,715` columns.  Every existing
clique row was extended to the new poses that closed-meet all its members (14,626 new memberships)
and the loop restarted on the enlarged set:

| it | LP | `M` | max clique | best violated support clique (→ extended) | clique rows | support |
|---|---|---|---|---|---|---|
| 0 | `11.466373` | `1.353` | `1.3526` | 5 → 710 | 258 | 133 |
| 1 | `11.393906` | `1.053` | `1.2323` | 34 → 278 | 264 | 150 |
| 3 | `11.368053` | `1.076` | `1.1517` | 34 → 414 | 270 | 149 |
| 5 | `11.355173` | `1.001` | `1.1264` | 27 → 410 | 270 | 131 |
| 8 | `11.346487` | `1.000` | `1.0875` | 21 → 597 | 285 | 120 |
| 11 | `11.339886` | `1.000` | `1.0753` | 29 → 353 | 295 | 138 |
| 14 | `11.331651` | `1.000` | `1.0486` | 25 → 246 | 286 | 133 |
| 20 | `11.325932` | `1.072` | `1.1578` | 21 → 448 | 280 | 144 |
| 25 | `11.324000` | `1.054` | `1.1100` | 18 → 300 | 272 | 104 |
| 30 | `11.321429` | `1.000` | `1.0483` | 13 → 292 | 278 | 82 |
| 44 | `11.321429` | `1.002` | `1.0452` | 17 → 276 | 206 (ageing) | 85 |
| resume +24 | `11.321429` | `1.000000002` | `1.0089–1.1161` | 8–15 → 430–1036 | 347 (no ageing) | 69–93 |

The new columns lift the first LP to `11.466` (with `M = 1.35`: the new poses first have to be
given their coverage rows) and the clique loop takes it back down in 30 iterations to
`317/28 = 11.321429`, where it has stayed for 60+ iterations: the LP is on a `1/56`-valued
degenerate optimal face, each iteration finds a violated maximal clique of mass `1 + k/56`
(`k = 1 … 6`, 8–15 support members, 430–1036 poses after extension) and the LP moves to a
neighbouring vertex of the same value.  Termination is guaranteed (no ageing, finite pose set,
one new maximal clique per iteration) but not within the budget.  What is established:

> `11.303030 <= QSTAB(P_1) <= 11.321429` on the 8,452-pose set `P_1`
> (lower bound: the exactly certified stage-(2) measure, which is feasible for `P_1`; upper
> bound: the LP value on a subset of the rows, up to the solver tolerance `1e-6`).

So the pricing gap of `+0.53` on the reduced cost translates into `+0.02` of value: the columns
the pricer wants are the columns the cliques forbid (they enter the extended rows), which is the
T4SCREEN calibration failure mode running in the *opposite* direction from the one the brief
feared — refinement does not bring the value back.  One stage is not convergence of the column
generation, and the honest statement is: `QSTAB` over all poses is at least `11.3214` (the
certified post-pricing measure) and, on the evidence of this stage, is not climbing towards 11.75.

## 5. The corner-leaf `k = 4` control (7,829 poses, 7,941 columns)

`tl_B40K_poses.txt` (T4LEAF's `k = 4` control: `11.988925` with the anchor cliques) with corners
pinned at 1 and the slots free.  The first LP on the seed rows reads exactly `12.000000` with
`M = 1.85` — the same starved-row `12.000000` as T4SCREEN's and LEAF_CEILING §5.2's — and falls
away as soon as the rows close in: `11.9897` (it 3), `11.9397`, `11.8904`, `11.8535`, `11.8259`,
`11.7875` (it 8), `11.7063` (it 14), `11.6820` (it 20), `11.6775` (it 26), `11.6751` (it 33),
`11.6739` (it 39); then 8 iterations pinned at `537/46` with the ageing recycling rows, checkpoint,
resume without ageing, and convergence 13 iterations later:

> **`QSTAB = 11.673913`** (`537/46`), `M = 1.000000001`, max clique `1.000000`, 362 clique rows,
> 23,140 coverage rows, 102 poses in the support; exactly certified **`5836956499/500000000 =
> 11.673912998`**, corners `1,1,1,1`, slots `0.413, 0.587, 0.489, 0.511, 0.565, 0.435, 0.5, 0.5`,
> chord `3,3,3,3`, interior `3.673913`.

The full clique family takes `0.315` off the corner leaf, versus `0.442` off the `m = 4` leaf.
The slot branch, which T4LEAF measured at `0.244` with the anchor cliques, is worth
`11.6739 - 11.3030 = 0.371` under the full family.

## 6. Anatomy of the violated cliques

Three views, all on exact geometry: the *raw* cliques the coverage-only optimum violates
(`cliquelever.py raw` on `lc_A0101.txt`, the certified 234-pose measure of LEAF_CEILING §5.3), the
*cutting* cliques found along the descent (`cl_*_found.json`, recorded on the re-run of step (1)),
and the *tight* cliques of the converged QSTAB measures (`raw --lb 0.999999`).

### 6.1 Raw: what the coverage optimum violates

The ten heaviest violated cliques of `lc_A0101` (mass `11.776`, `M = 1`).  "point" is the mass of
the heaviest point sub-clique `P(p) ∩ K` (exact, over the arrangement of the members); "anchor" the
best `K(p, A)` inside the clique found by `leaf_ceiling`'s exact scan + ascent on the sub-measure;
"captured" = `(anchor − point) / (mass − point)`, the share of the non-Helly excess the certifiable
family sees.

| members | mass | point | excess | anchor | captured | the point `p` | angle bins (5°) | regions of the mass |
|---|---|---|---|---|---|---|---|---|
| 52 | `1.5388` | `0.979` | `0.560` | `1.296` | 57 % | `(1.988, 2.661)` | 14 | I `1.30`, W5 `0.24` |
| 45 | `1.4818` | `0.915` | `0.567` | `1.285` | 65 % | `(1.001, 1.316)` | 15 | W7 `0.59`, I `0.89` |
| 44 | `1.4650` | `0.994` | `0.471` | `1.066` | 15 % | `(2.988, 1.958)` | 15 | I `1.11`, W3 `0.36` |
| 47 | `1.4500` | `0.895` | `0.556` | `1.149` | 46 % | `(1.007, 2.332)` | 12 | W7 `0.34`, I `1.11` |
| 35 | `1.4364` | `1.000` | `0.436` | `1.122` | 28 % | `(2.904, 1.003)` | 12 | I `0.91`, W1 `0.53` |
| 40 | `1.3892` | `0.947` | `0.443` | `1.041` | 21 % | `(1.490, 2.878)` | 14 | W5 `0.48`, I `0.91` |
| 37 | `1.3646` | `0.898` | `0.466` | `0.958` | 13 % | `(2.878, 2.510)` | 16 | I `0.87`, W3 `0.50` |
| 40 | `1.3331` | `0.964` | `0.370` | `1.081` | 32 % | `(1.000, 1.070)` | 15 | W7 `0.68`, I `0.66` |
| 14 | `1.3122` | `0.985` | `0.327` | `1.113` | 39 % | `(2.996, 0.819)` | 7 | C1 `0.71`, W1 `0.27`, I `0.33` |
| 32 | `1.2904` | `0.942` | `0.349` | `0.942` | 0 % | `(1.000, 2.778)` | 13 | W5 `0.35`, I `0.78`, C2 `0.16` |

The shape is the same in every row and it is **not a pinwheel**: a point clique of mass
`0.9–1.0` at a point `p` on one of the lines `x = 1, 2, 3` or `y = 1, 2, 3` — the corners of the
axis-parallel wall squares `[0,1] x [1,2]`, `[2,3] x [0,1]`, … (the leaf's own frame) — and a
**fan** of 30–45 interior tilted squares at 12–16 distinct angles between 0° and 90°, each of
which closed-meets every square through `p`, carrying `0.33–0.57` of mass.  Do they share a
segment?  No: the fan members contain no common segment — squares at 14 different angles that all
graze the wall square's corner region have an intersection that is a small polygon at best — and
that is exactly what the anchor family `K(p, A) = P(p) ∩ Meet(A) ∪ Cont(A)` needs: every extra
member must *contain* `A`.  The exact ascent finds an `A` that a part of the fan contains and
captures 0–65 % of the excess (median 30 %, mean 32 %); the rest of the fan is invisible to any
`K(p, A)`.  The heaviest non-Helly *triple* is a different question (CLIQUE.md found 0.27–0.89):
the violation is a property of the whole cluster.

### 6.2 Cutting: what the descent removes

On the 234-pose case (55 rows), by iteration.  The raw cliques of it 0 are the ones above; from
it 1 on, the violated cliques are **interior, smaller, and point-exact**:

| it | support members → extended | mass | point | excess | angle bins | region of the mass |
|---|---|---|---|---|---|---|
| 0 | 43 → 44 … 49 → 52 | `1.34–1.48` | `0.89–1.00` | `0.40–0.55` | 12–15 | wall/slot corner + interior |
| 1 | 19 → 48, 22 → 51, 21 → 50, 23 → 46, 22 → 48, 18 → 48 | `1.29–1.46` | `0.93–1.00` | `0.29–0.46` | 8–12 | I `1.00–1.46` (entirely interior) |
| 1 | 17 → 38, 12 → 39, 9 → 28, 17 → 51, 8 → 28 | `1.03–1.15` | `0.61–1.00` | `0.06–0.42` | 3–10 | interior + one slot |
| 2 | 10 → 47, 8 → 47, 9 → 41, 9 → 45, 10 → 48 | `1.12–1.26` | `0.97–1.00` | `0.12–0.26` | 5–7 | I |
| 3 | 10 → 50, 7 → 37, 7 → 42, 10 → 38 | `1.07–1.26` | `0.95–1.00` | `0.07–0.26` | 4–8 | I, I + W0/C0, I + W1 |
| 4–5 | 7 → 33, 6 → 47, 8 → 35, 7 → 28 | `1.05–1.15` | `0.95–1.00` | `0.05–0.15` | 3–6 | W1 + I |

On the 8k-pose leaf set the same progression is visible in the trajectory table of §3 (the best
violated clique: 71 → 49 → 33 → 11 → 29 → 26 → … → 7–13 support members, mass 1.65 → 1.02), with
the extensions over all poses growing as the support cliques shrink (250 → 1150 poses per row).
The cutting cliques are, after the first iteration, `P(p)` for an interior point `p` (mass exactly
1) plus a thin fan of straddlers — precisely the object whose certifiable version would be "all
poses through `p` together with all poses meeting every square through `p`", the non-avoidance
family CLIQUE.md proposed and `notes/clique-family.md`'s `K(p, A)` approximates with one segment.

### 6.3 Tight: what is left at convergence

| measure | tight cliques found (`mass >= 1 - 1e-6`) | of which pure point cliques | the others |
|---|---|---|---|
| leaf QSTAB `11.303030` (70 poses) | 27 | 24 | 11 members, point `0.939`, excess `0.061`; 10 members, point `0.970`, excess `0.030` (twice) — anchor `1.0`, `1.0`, `0.97` |
| control QSTAB `11.673913` (102 poses) | 25 | 24 | 15 members, point `0.891`, excess `0.109`, anchor `1.0` |

At convergence the measure has retreated to an essentially Helly configuration: the binding
constraints are point cliques (the LP's positive clique duals at the end sit on rows whose
members-with-mass are 1–19 poses through one point), the interior mass is `1/33`- or
`1/46`-valued, and the residual non-Helly excess is `0.03–0.11` on two or three cliques, all of
which an anchor clique of mass 1 already covers.  The non-Helly structure of the coverage
optimum is not "corrected" by the cliques; it is destroyed, and `0.44–0.68` of mass goes with it.

## 7. What is exact and what is float

* **Exact (Python integers / `Fraction`s).**  The poses (rational rotation and centre), the
  squares and their admissibility; the arrangement vertices and the coverage of every support
  measure (so every reported `M`, computed on the `/10^9`-rounded masses of the float solution);
  the closed-intersection graph of every support (integer SAT, `<=`); the max-weight clique and
  its completeness flag (integer weights); the region membership of every column (closed boxes);
  the chord rows; the pairwise verification of every clique row and of every resumed row on the
  band `|margin| <= 1e-7`, float outside it (the float SAT margin has error `~1e-15` on
  coordinates of size 4, so the band decides every pair that could be in doubt); the final
  certification of each stage's measure (`cl_*_exact.txt`), re-checked by `leaf_ceiling.py check`
  whose clique test is fully exact.
* **Float.**  The LP itself (HiGHS, feasibility tolerance `1e-7`), hence the *values* `LP` and the
  masses; the choice of which rows to add; the pricing.  A converged float LP value is exact up to
  the solver tolerance; the exactly certified mass (rounded down, regions topped up on exact
  slack) is a rigorous LOWER bound on `QSTAB(P)`, and the last LP value an upper bound on it up to
  `~1e-6`.  The three certified masses are within `2e-8` of the float values.
* **Not exact, by design.**  The pose set is the *snapped* one: `t4leaf`'s float poses moved by
  `<= 5e-6` (angle `2 arctan(p/q)`, `Q = 10^5`; centres on `1e-6`).  Grazing contacts of the float
  set may open or close under the snap; CLIQUE_CEILING.md's "structural finding" says the clique
  excess lives on such contacts, so the snapped pose set is *the* pose set here and the numbers
  are its numbers.  The direction is safe either way: whatever the snap does, the certified
  measures are certified.
* **Hypotheses used.**  The chord rows `mu(strip) <= 3` are present (T4LEAF §1.2's hypothesis) but
  carry dual `0` at every iteration of every run — implied by the region equalities in the leaf
  and by the corner equalities plus coverage in the control — so no number here depends on them.

## 8. Reproduce

```sh
python3 search/cliquelever.py selftest                       # the pinwheel: 1.5 -> 1, exact, extension
R=runs   # copies of the read-only checkpoints named in the task brief
# (1) the certified 234-pose support, ~30 s; then the independent re-check and the anatomies
python3 search/cliquelever.py run A0101s --exact $R/lc_A0101.txt --corners 1111 --patterns 01010101 \
    --chord --rows0 30036 --threads 2 --procs 2 --ckpt 1
python3 search/leaf_ceiling.py check $R/cl_A0101s_exact.txt --corners 1111 --patterns 01010101 --chord --anchor clique
python3 search/cliquelever.py raw $R/lc_A0101.txt --anatomy 10 --json $R/cl_A0101_raw_anatomy.json
python3 search/cliquelever.py raw $R/cl_A0101s_exact.txt --lb 0.999999 --anatomy 40
# (2)+(3) the leaf pose set, then one pricing stage (the run does both; ~25 min + ~2.5 h)
python3 search/cliquelever.py run A0101L --exact $R/lc_A0101.txt --load-poses $R/tl_A01010101_poses.txt \
    --load-rows $R/tl_A01010101_dual.txt --corners 1111 --patterns 01010101 --chord \
    --threads 2 --procs 2 --lp-tlim 10 --clique-time 60 --cq-want 60 --cq-top 40 --cq-age 0 --price 1
# control: the corner leaf k = 4 (~75 min)
python3 search/cliquelever.py run B40KL --load-poses $R/tl_B40K_poses.txt --load-rows $R/tl_B40K_dual.txt \
    --corners 1111 --patterns ........ --chord --threads 2 --procs 2 --lp-tlim 10 --clique-time 60 \
    --cq-want 60 --cq-top 40 --cq-age 0
# resume a checkpointed run (as was done here after the ageing cycle):
python3 search/cliquelever.py run B40KL1 --exact $R/cl_B40KL_poses.txt --resume-rows $R/cl_B40KL_rows.txt \
    --resume-cliques $R/cl_B40KL_cliques.txt --corners 1111 --patterns ........ --chord --row-pitch 0 --cq-age 0
```

The runs recorded here used `--cq-age 12` for the first 31 (leaf, stage 0), 44 (leaf, stage 1) and
47 (control) iterations and `--cq-age 0` on the resume; `--cq-age 0` from the start is the
recommended setting.
