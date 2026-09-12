# rank-diag: what non-clique structure carries the residual excess of the QSTAB measures?

Task: `tasks/rank-diag/README.md`.  Code: `search/rankdiag.py` (imports `leaf_ceiling.py` for the
exact geometry; nothing in `verify/`, `lean/`, `certificates/` is touched).  Runs in this
worktree's `runs/` (gitignored; every number is quoted here).  Read against
`search/CLIQUELEVER.md` §6, `notes/clique-family.md` §0–§1, `search/BOXCLIQUE.md`.

## 0. Verdict, up front

**Build it — but the object that bites is not an odd hole of poses, it is a pentagon of anchors.**

Three findings, in order of how much they change the plan.

1. **`alpha = 11` on every converged QSTAB support, exactly, with a complete branch and bound.**
   The non-clique gap `mu(V) - alpha(G)` is **`0.32` / `0.44` / `0.67` / `0.70` / `0.79`** on the
   five measures.  Each support really does carry 11 pairwise-disjoint closed unit squares and
   provably not 12, so the gap is exactly the amount of mass that the *rank* inequality
   `mu(X) <= alpha(G[X])` takes off and no clique inequality can.

2. **The odd-cycle family is worth exactly zero.**  Over all five supports there is **not one**
   chordless `C5` or `C7` with mass above its rank bound, **not one** odd antihole on 7 vertices
   above 2, and **not one** 5-wheel above 2 — and it is not close: the heaviest chordless `C5`
   anywhere reaches `1.46` against the bound 2, the heaviest `C7` `1.85` against 3, the heaviest
   7-antihole `1.21` against 2, the heaviest 5-wheel `1.47` against 2.  On the leaf that is `0.00`
   of the `0.32` gap, against the brief's `>= 0.1` threshold.  The mass sits on 71–152 poses with
   a median of `0.03–0.07` each; no five or seven of them ever add up to 2 or 3.  A separator built
   on holes, antiholes and wheels of *poses* would never fire.

3. **A five-anchor rank-2 family does bite, on all five measures, and it is certifiable with no
   geometry beyond convexity.**  Take five points `A_0..A_4`, the five closed segments
   `E_i = [A_i, A_{i+1 mod 5}]` of the pentagon they span, and the pieces
   `X_i = {S admissible : E_i ⊆ S}`.  Any two members of the same piece both contain `A_i`; any
   member of `X_i` and any member of `X_{i+1}` both contain `A_{i+1}`; so an independent set of
   `∪ X_i` injects into an independent set of `C5` and

   > **`mu(X_0 ∪ … ∪ X_4) <= alpha(C5) = 2`** — valid for every packing of disjoint closed unit
   > squares, hence a valid cut for the LP, which the measures below violate.

   Hill-climbing the five anchors over the exact arrangement vertices of each support finds

   | measure | pentagon `mu` | excess over 2 | share of the gap | poses hit |
   |---|---|---|---|---|
   | `cl_A0101L1` | **`2.285714`** | `+0.285714` | **89 %** | 19 |
   | `cl_B40KL1` | **`2.271739`** | `+0.271739` | 40 % | 28 |
   | `cl_A0101L2` | **`2.435484`** | `+0.435484` | **100 %** | 32 |
   | `cl_B40KL2` | **`2.213253`** | `+0.213253` | 27 % | 45 |
   | `cl_PUREQ` | **`2.076164`** | `+0.076164` | 11 % | 34 |

   Each of these was re-derived from the exact rational anchors by a second, independent code
   path (`sq_contains` on both endpoints of each segment) — membership **matches** — and each
   witness set had its `alpha(G[X]) = 2` re-confirmed by a fresh unseeded complete branch and
   bound, as Lemma 0 predicts.  `rankdiag.py --recheck runs/rankdiag.json` then rebuilds all five
   from the measure files and the stored rational anchors and checks Lemma 0 by brute force —
   **every one of the `C(|X|,3)` triples (up to 14,190 of them) has a meeting pair**:
   `ALL PENTAGONS VERIFIED`.  On the leaf the single cut is `0.286` against a `0.321` gap; on
   `cl_A0101L2` it is the entire gap.

**Why this is the family to build.**  In the notation of `notes/clique-family.md` §1 this is the
`m`-anchor family `K(A_1,…,A_m)` with `m = 5` and guaranteed-meet graph `H = C5`, in the
*special case where consecutive anchors meet* (they share an endpoint), so **no "meets" filter is
needed at all** — the pieces are plain containment sets `{S : E_i ⊆ S}`.  That is precisely a box
clique of `search/BOXCLIQUE.md` with one rule relaxed: a box clique is a union of boxes whose
*cores* (the rectangles every member square contains) **pairwise** meet, credited weight 1; the
pentagon is a union of five boxes whose cores meet **cyclically**, credited weight 2.  The
verifier's `cores_meet` (exact separating-axis on rotated rectangles in `i128`) is already the
test; what changes is which pairs it is asked about and what right-hand side the row carries.
Nothing new has to be trusted, no new geometric lemma is needed, and the validity proof is the
four-line argument above plus `alpha(C5) = 2`.  It would be the first inequality in this project
that is neither a coverage row nor a clique row.

**The next experiment, concretely.**  `rankdiag.py --step 4` already *is* a separator: given a
measure it returns a violated pentagon in under a minute.  Dropping it into `cliquelever.py`'s
loop as a third separation stage (after coverage and cliques), with the pentagon rows extended
over the whole pose set the way clique rows already are, measures the thing this document does not:
how far the QSTAB value actually falls.  The upper limit is `mass - alpha = 0.32–0.79` per measure
— that is all the rank family of *any* kind can ever be worth on these supports, and it is
`11.32 - 11 = 0.32` on the leaf, so even a perfect rank separator leaves the leaf at 11 and the
`t = 4` ceiling question unchanged.  What it would change is the shape of the argument: a
certified bound of 11 on a leaf's own pose set, from rows a verifier can check.

**Where the excess lives.**  Every winning pentagon is anchored on the wall-square corner lines
`x ∈ {1,2,3}` or `y ∈ {1,2,3}` and wraps a *grid vertex* — `(3,2)` for both `A0101L1` and
`B40KL1`, `(2,3)` for `A0101L2`, `(1,2)` for `B40KL2` and `PUREQ` — with three or four anchors on
the line and one or two just inside.  The mass it collects is split between the interior region
`I` and the adjacent wall region (`W3`, `W5`, `W6/W7`), it spans 14–32 distinct angles, and no
point is common to all members.  This is the same object `CLIQUELEVER.md` §6.1 found as the
violated *clique*: a wall-square corner plus a fan of tilted interior squares.  The clique family
sees only the part of the fan through one point; the pentagon sees a fan that wraps the corner and
pays for it with a right-hand side of 2 instead of 1.

**What is NOT there.**  There is no small combinatorial witness of the classical kind.  Minimal
violated rank subsets (no single pose removable) have **56–103 members** and an excess of one mass
quantum (`0.035714 = 1/28`, `0.010870`, `0.016129`, `0.001446`, `0.025388`); the exhaustive search
over sets of `<= 8` poses with `alpha <= 2` finds nothing on any measure.  The diffuseness is real
— it is just that diffuseness of the *support* does not imply diffuseness of the *family*: the
pentagon witness is 19–45 poses, but it is generated by five points, and it was the anchor
geometry, not the subset search, that found it (on `cl_B40KL2` and `cl_PUREQ` the subset search
found no violated `alpha = 2` set at all, while the pentagon search found ones of mass `2.213` and
`2.076`).

## 1. What is computed

`search/rankdiag.py` reads a `leaf_ceiling.py` measure file, rebuilds the exact squares
(`read_measure` + `build_squares`; masses become integers over the common denominator
`DM = 10^9`), and builds the exact closed-intersection graph `G` of the support
(`leaf_ceiling.closed_graph`: integer separating-axis SAT with `<=`, so touching counts as
meeting — the semantics of `notes/clique-family.md` §0).

* **step 1 `alpha`.**  `alpha(G)` = maximum clique of the complement, by a bitset Tomita branch
  and bound with a greedy-colouring bound, exact and with a completeness flag.  Reports
  `mu(V) - alpha(G)`.
* **step 2 `cycles`.**  Every chordless cycle of length 5 and 7 in `G` with mass `> (k-1)/2`;
  every chordless 7-cycle of the *complement* (= odd antihole on 7 vertices) with mass `> 2`;
  every 5-wheel (a vertex adjacent to all of a chordless `C5`) with hub + rim mass `> 2`.  The
  enumeration is canonical (the heaviest member is the root, one traversal direction) and pruned
  by "the partial path can gain at most `(k - len) * w[root]`", which is exact and, because the
  root is the heaviest member, complete.  When nothing is violated the threshold is walked down
  until something is found, so the *heaviest* structure of each kind is reported.
* **step 3 `search`.**  A heuristic hunt for `X` with `mu(X) - alpha(G[X])` large: (a) an ascent
  that deletes the lightest pose whose deletion lowers `alpha`; (b) minimalisation to sets from
  which no single pose can be removed; (c) an `alpha`-level sweep — for each `k`, greedy growth
  of a heaviest `X` with `alpha(G[X]) <= k`, seeded from the levels below; (d) an exhaustive
  search over sets of `<= 8` poses with `alpha <= 2`.  Everything any of these touches is offered
  to a size-vs-excess Pareto front, which is what gets reported.  One economy makes this cheap: if
  `W` is a maximum independent set of `G[X]` and `v ∉ W` then `alpha(X\v) = alpha(X)`, so only the
  `alpha` members of `W` can lower it — each round costs one `alpha` call plus `<= 11` tests, not
  `|X|`.
* **step 4 `anchors`.**  The separation problem for the pentagon family of §0.3: enumerate the
  exact arrangement vertices of the support (`leaf_ceiling.enumerate_vertices` +
  `incidences`), reduce them to distinct point cliques, keep the heaviest 500, and hill-climb five
  anchor indices by coordinate descent from many random restarts, maximising
  `mu(∪_i (Cont(A_i) ∩ Cont(A_{i+1})))` — a square is convex, so it contains the segment
  `[A_i, A_{i+1}]` iff it contains both endpoints.

`python3 search/rankdiag.py --selftest` checks the machinery on two synthetic measures: a ring of
five unit squares that meet only consecutively (`alpha = 2`, mass `2.5`, one violated chordless
`C5` of excess `0.5`, no `C7`, no antihole, and the exhaustive small search finds it) and a ring of
seven (`alpha = 3`, one violated chordless `C7` of mass `3.5`, no `C5`).

## 2. Step 1: `alpha` and the gap

| measure | poses | mass | `alpha(G)` | complete B&B | **gap `mu - alpha`** | edges (density) |
|---|---|---|---|---|---|---|
| `cl_A0101L1` (leaf `01010101`) | 71 | `11.321429` | **11** | yes | **`0.321429`** | 834 (`0.336`) |
| `cl_A0101L2` (leaf `01010101`, later run) | 94 | `11.435484` | **11** | yes | **`0.435484`** | 1506 (`0.345`) |
| `cl_B40KL1` (corner leaf `k = 4`) | 102 | `11.673913` | **11** | yes | **`0.673913`** | 1836 (`0.356`) |
| `cl_B40KL2` (corner leaf, later run) | 149 | `11.785783` | **11** | yes | **`0.785783`** | 4383 (`0.398`) |
| `cl_PUREQ` | 152 | `11.702398` | **11** | yes | **`0.702398`** | 3850 (`0.335`) |

(Pose counts are the nonzero-mass columns the reader rebuilds; `sym = 1` on all five, so one
square per pose.  The brief's "70 / 102 / 152" refer to the same files.)

Every `alpha` took under `0.05 s`.  `alpha = 11` on all five is the expected answer and the
reassuring one: the supports are finite pose sets on which the LP sits at `11.32–11.79`, and an
independent set of 12 would have been 12 pairwise-disjoint closed unit squares in `[0,4]²` — a
disproof of the conjecture.  Each maximum independent set found is the four corner squares (or
three of them plus a near-corner) together with seven interior/wall squares; the corner squares
carry mass `0.87–1.00` and the rest `0.005–0.39`, which is why the gap is small: the measure puts
its `11.3–11.8` on a support whose integral optimum is 11, and the surplus is spread over the
light poses.

## 3. Step 2: the odd-cycle family finds nothing

| measure | violated `C5` (`>2`) | violated `C7` (`>3`) | violated 7-antihole (`>2`) | violated 5-wheel (`>2`) |
|---|---|---|---|---|
| `cl_A0101L1` | **0** | **0** | **0** | **0** |
| `cl_B40KL1` | **0** | **0** | **0** | **0** |
| `cl_A0101L2` | **0** | **0** | **0** | **0** |
| `cl_B40KL2` | **0** | **0** | **0** | **0** |
| `cl_PUREQ` | **0** | **0** | **0** | **0** |

Not marginally zero.  The heaviest structure of each kind, anywhere:

| measure | heaviest chordless `C5` (bound 2) | heaviest `C7` (bound 3) | heaviest 7-antihole (bound 2) | heaviest 5-wheel (bound 2) |
|---|---|---|---|---|
| `cl_A0101L1` | `1.392857` | `1.750000` | `1.500000` | `1.464286` |
| `cl_B40KL1` | `1.119565` | `1.456522` | `1.097826` | `1.108696` |
| `cl_A0101L2` | `0.967742` | `1.322581` | `1.209677` | — (no `C5` with rim `> 1`) |
| `cl_B40KL2` | `1.027711` | `1.220241` | `1.087711` | `1.053494` |
| `cl_PUREQ` | `1.458392` | `1.853315` | `1.187588` | `1.466855` |

The arithmetic reason is blunt.  A violated `C5` needs its heaviest member to carry more than
`2/5 = 0.4`; only 4–7 poses per measure do, and on four of the five they are exactly the four
corner squares (mass `0.87–1.00`), which meet almost nothing and cannot sit in a chordless 5-cycle
with four more heavy poses.  A
violated `C7` needs `> 3/7 = 0.43` on its heaviest member; a violated 7-antihole `> 2/7 = 0.29`.
At convergence the mass is spread over 71–152 poses with individual masses from `1.00` down to
`0.0005` and a median of `0.03–0.07` (`CLIQUELEVER.md` §6.3: "essentially Helly"),
and five or seven such poses simply do not add up to 2 or 3.

**The answer to the brief's question, stated plainly: the odd-cycle family accounts for `0.000` of
the `0.32 / 0.44 / 0.67 / 0.70 / 0.79` gap on every measure.**  The brief's `>= 0.1` threshold is
missed by the whole gap.

## 4. Step 3: where the violated rank structure actually is

Violated rank inequalities are everywhere — they are just large.  The size-vs-excess front (the
largest `mu(X) - alpha(G[X])` found at each `|X|`, over everything the search touched):

| measure | smallest `\|X\|` with excess `>= 0.1` | smallest with excess `>=` half the gap | smallest with excess `>=` the gap | best excess found (at `\|X\|`) |
|---|---|---|---|---|
| `cl_A0101L1` | **17** (`alpha = 2`, `mu = 2.107143`) | 19 | 25 | `+0.785714` (66) |
| `cl_B40KL1` | **32** (`alpha = 2`, `mu = 2.108696`) | 48 | 69 | `+0.858696` (89) |
| `cl_A0101L2` | **46** (`alpha = 5`, `mu = 5.112903`) | 49 | 53 | `+1.000000` (87) |
| `cl_B40KL2` | **51** (`alpha = 4`, `mu = 4.105783`) | 88 | 132 | `+0.785783` (137) |
| `cl_PUREQ` | **37** (`alpha = 4`, `mu = 4.100141`) | 108 | 144 | `+0.702398` (144) |

Every set on those fronts was re-verified after the fact — mass resummed, `alpha` recomputed by a
fresh unseeded complete branch and bound — with **0 mismatches** (25 / 58 / 36 / 78 / 69 sets).

Three things are worth separating.

* **Minimal violated subsets are useless as witnesses.**  Minimalising (delete poses while
  `mu(X) > alpha(G[X])` survives) from the full support, greedily and with 20 randomised restarts,
  lands on sets of **56 / 65 / 69 / 86 / 103** poses whose excess is one mass quantum
  (`0.035714 = 1/28`, `0.010870`, `0.016129`, `0.001446`, `0.025388`).  That is the granularity
  floor, not structure: the minimalisation stops when the next removal would cross zero.  These
  are the "large `X`, no small certifiable witness" sets the brief warned about, and they are what
  a naive minimal-witness search would return.
* **The excess is *not* tied to large `alpha`.**  On `cl_A0101L1` the front reaches the full gap
  `0.321429` at `alpha = 3` on `|X| = 25`, and `0.178571` (56 % of the gap) at `alpha = 2` on 19
  poses; on `cl_B40KL1` and `cl_A0101L2` the subset search also found violated `alpha = 2` sets
  (`mu = 2.119565` on 33, `mu = 2.096774` on 32).  That is what made step 4 worth doing: an
  `alpha = 2` set of 20–40 poses is exactly the shape a five-piece anchor family can cover.  On
  `cl_B40KL2` and `cl_PUREQ` the whole front bottoms out at `alpha = 4` — the subset search found
  no violated set of rank 3 or less at all — and step 4 found a violated rank-2 one on both, which
  is the measure of how weak the subset heuristics are relative to searching the anchor geometry
  directly.
* **Nothing small.**  The exhaustive search over every `X` with `|X| <= 8` and `alpha(G[X]) <= 2`
  (grown so the complement stays triangle-free, mass-pruned) returns **0** on all five measures.
  There is no 5-, 6-, 7- or 8-pose certificate of any rank-2 kind.

## 5. Step 4: the pentagon anchor family, and its anatomy

**The inequality.**  For points `A_0..A_4`, segments `E_i = [A_i, A_{i+1 mod 5}]`, pieces
`X_i = {S admissible : E_i ⊆ S}`:

* two members of the same `X_i` both contain `A_i`, so they closed-intersect;
* `S ∈ X_i` and `S' ∈ X_{i+1}` both contain `A_{i+1}`, so they closed-intersect;
* hence any pairwise-disjoint subfamily of `∪_i X_i` meets each `X_i` at most once and never two
  consecutive ones, i.e. it injects into an independent set of `C5`;
* `alpha(C5) = 2`, so **any packing of pairwise-disjoint closed unit squares in `[0,4]²` places
  at most 2 of its squares in `X_0 ∪ … ∪ X_4`** — i.e. `mu(∪ X_i) <= 2` is valid for the integer
  hull and may be added to the LP as a cut, exactly as a clique inequality is.  (Fractional `mu`
  need not satisfy it; that is the point.)

This is Lemma 0 of `notes/clique-family.md` with `m = 5` pieces and `H = C5`, in the sub-case where
consecutive anchors intersect so the "meets" filter degenerates.  No non-avoidance lemma, no
wall lemma, no geometry beyond "a square is convex".  It generalises verbatim: for any `m` points
`A_0..A_{m-1}` in a cycle with `m` odd, the same argument gives `mu(∪ X_i) <= alpha(C_m) =
(m-1)/2`; `m = 5` is the first case that is not a clique.  Note that the `C_m` here is a cycle of
*anchors*, of which there are five — not a cycle of *poses*, of which §3 shows there are none.
That distinction is the whole content of this document.

**The separation results** (hill climb over the top 500 point cliques of each support's
arrangement, 2000 random restarts plus heavy-point seeds; well under a minute per measure —
`192 s` for all four steps on all five measures, single-threaded):

| measure | gap | pentagon `mu` | excess | share | `\|X\|` | `alpha(G[X])` re-checked | re-derived from exact anchors |
|---|---|---|---|---|---|---|---|
| `cl_A0101L1` | `0.321429` | `2.285714` | `+0.285714` | 89 % | 19 | `2` (complete) | matches |
| `cl_B40KL1` | `0.673913` | `2.271739` | `+0.271739` | 40 % | 28 | `2` (complete) | matches |
| `cl_A0101L2` | `0.435484` | `2.435484` | `+0.435484` | 100 % | 32 | `2` (complete) | matches |
| `cl_B40KL2` | `0.785783` | `2.213253` | `+0.213253` | 27 % | 45 | `2` (complete) | matches |
| `cl_PUREQ` | `0.702398` | `2.076164` | `+0.076164` | 11 % | 34 | `2` (complete) | matches |

**Anatomy** — where they are, and whether they wrap a grid vertex:

| measure | anchors `A_0..A_4` (decimal) | grid vertex | regions of the mass | distinct angles |
|---|---|---|---|---|
| `cl_A0101L1` | `(3.004,2.986) (3.178,2.494) (3.015,1.929) (2.786,1.988) (2.470,2.152)` | **`(3,2)`**, 5 of 19 members contain it | `I 1.286`, `W3 0.679`, `W2/W3 0.321` | 14 |
| `cl_B40KL1` | `(2.284,1.723) (3,1.082) (3,1.999) (3,2.384) (2.340,1.982)` | **`(3,2)`** 7, `(2,2)` 6 of 28 | `I 1.272`, `W3 0.511`, `W2 0.489` | 17 |
| `cl_A0101L2` | `(1.879,2.519) (2.000,2.135) (2.077,3) (1.5,3.993) (1.014,3)` | **`(2,3)`**, 11 of 32 | `I 1.435`, `W5 0.694`, `W4/W5 0.306` | 14 |
| `cl_B40KL2` | `(1,2.001) (1,1.029) (1.521,1.855) (1.484,1.998) (1,2.297)` | **`(1,2)`**, 11 of 45 | `I 1.213`, `W6 0.555`, `W7 0.445` | 32 |
| `cl_PUREQ` | `(1.656,1.767) (1,1.160) (0.851,2.001) (1.000,2.787) (1.329,2.284)` | `(1,2)` 3, `(2,2)` 2 of 34 | `I 1.451`, `W7 0.315`, `W6 0.310` | 28 |

Three or four of the five anchors lie exactly on one of the wall-square corner lines
`x ∈ {1,2,3}` / `y ∈ {1,2,3}`, and the remaining one or two sit just inside it.  The pentagon is
therefore a thin sliver wrapped around a grid vertex on a wall boundary — not `(2,2)` itself in any
of the five, but `(3,2)`, `(2,3)` or `(1,2)`, i.e. the corner of an axis-parallel wall square.  The
mass it collects is split roughly half interior / half wall, and no point lies in every member
(`common grid vertex: none` in all five), which is exactly why no *clique* inequality catches it.
This is the same geometry `CLIQUELEVER.md` §6.1 reported for the violated cliques of the coverage
optimum; the difference is that the fan wraps the corner rather than passing through one point, so
it takes a rank-2 right-hand side rather than a rank-1 one to bound.

**What this does and does not claim.**  It claims: on each converged QSTAB measure there is a
single valid inequality, outside the clique family, certifiable from convexity alone, that the
measure violates by `0.08–0.44`.  It does not claim that adding it drops the LP value by that
much — the LP will redistribute, and one cut is one cut.  The pentagon rows are valid over the
whole continuum of admissible poses (the pieces are defined for all `S`, not just the support), so
they are genuine rows for the pricing LP; measuring what a separation loop over them is worth is
the next experiment, not this one.

## 6. What is exact and what is float

* **Exact (Python integers / `Fraction`s).**  The poses and squares (`read_measure`,
  `build_squares`); the closed-intersection graph (`closed_graph`, integer separating-axis SAT with
  `<=`); every mass comparison — all masses are integers over `DM = 10^9`, so `mu(X) > alpha`,
  `mu(C) > (k-1)/2` and `mu(pentagon) > 2` are integer comparisons with no tolerance anywhere;
  `alpha` and its completeness flag (integer branch and bound on the complement); the arrangement
  vertices (exact reduced rational triples) and the point-in-square incidence
  (`sq_contains`, integers); the region membership (closed boxes, exact `Fraction` comparisons);
  the grid-vertex containment tests.
* **Float, and only in ways that cannot change a verdict.**  Which subsets the heuristics try
  (the ascent, minimalisation, `alpha`-level sweep and pentagon hill climb are all heuristics —
  every set they return is then verified exactly); the decimals printed; the angles in degrees;
  and the centre-distance prefilter inside `leaf_ceiling.closed_graph`, which skips the exact
  pairwise test when the float squared centre distance exceeds `2 + 1e-6` — safe, because two unit
  squares have circumradius `sqrt(2)/2` each and cannot meet beyond centre distance `sqrt 2`, and
  the float error on coordinates of size 4 is `~1e-15`.
* **Heuristic, and labelled as such.**  Everything reported as a *maximum* over subsets — the
  size-vs-excess front, the `alpha`-level sweep, the "heaviest `alpha = 2` subset found", the best
  pentagon — is a search result, hence a **lower** bound on the true maximum.  The exhaustive
  statements are the ones that are complete: `alpha` itself (complete B&B), the odd-cycle /
  antihole / wheel enumerations at their rank thresholds (complete, by the mass-descending root
  prune), and the `|X| <= 8`, `alpha <= 2` search.  The pentagon search repeatedly beat the
  `alpha`-level sweep's own best `alpha = 2` set, which is a reminder of how loose the sweep is.
* **Not exact, by design.**  The measures are the snapped pose sets of `CLIQUELEVER.md` §7 (angle
  `2 arctan(p/q)`, `Q = 10^5`; centres on `1e-6`).  Grazing contacts may open or close under the
  snap; as there, the snapped set is *the* pose set here and the numbers are its numbers.

## 7. Reproduce

`~3.5 min` total for all five measures, single-threaded, no detached runs.

```sh
python3 search/rankdiag.py --selftest          # the 5-ring and the 7-ring: alpha, holes, exhaustive
R=runs                                         # copies of the read-only measures in the task brief
python3 search/rankdiag.py --step 4 --sweep 8 --restarts 20 --pent 2000 --cands 500 \
    --json $R/rankdiag.json \
    $R/cl_A0101L1_exact.txt $R/cl_B40KL1_exact.txt $R/cl_A0101L2_exact.txt \
    $R/cl_B40KL2_exact.txt $R/cl_PUREQ_exact.txt
# re-verify every pentagon from the JSON alone: rebuild the measures, take the exact rational
# anchors, recompute membership, and check Lemma 0 by brute force over all triples (~20 s)
python3 search/rankdiag.py --recheck $R/rankdiag.json
# step 1 only (alpha and the gap, < 1 s for all five):
python3 search/rankdiag.py --step 1 $R/cl_*_exact.txt
```

`--step 1` is `alpha`, `2` adds the odd-cycle enumerations, `3` adds the subset search, `4` adds
the pentagon separation.  The JSON carries, per measure, the maximum independent set, every
reported structure with its members, the size-vs-excess front, and the pentagon's anchors as exact
rationals (`anchors_exact`) so the inequality can be rebuilt without rerunning the search.

---

# Round 2: the family in the loop

Code: `search/rankfamily.py` (all of it) and a hook in `search/cliquelever.py` behind `--pent`
(imports, `Lever.add_pgon` plus one LP block, one separator call in `run_stage`, two guard lines
in `finalize`, the CLI flags).  With `--pent` empty — the default — the loop is bit-for-bit what
it was.

## 8. The generalised family: odd polygons of `k` anchors

§5's pentagon is the case `k = 5` of

> **Definition.**  For `k` odd and points `A_0 … A_{k-1}`, let `E_i = [A_i, A_{i+1 mod k}]` be the
> `k` sides of the polygon they span and `X_i = { S admissible : E_i ⊆ S }`.
>
> **Lemma.**  `mu(X_0 ∪ … ∪ X_{k-1}) <= alpha(C_k) = (k-1)/2` for every packing of pairwise
> disjoint closed unit squares.
>
> *Proof.*  Two members of `X_i` both contain `A_i`.  A member of `X_i` and one of `X_{i+1}` both
> contain `A_{i+1} ∈ E_i ∩ E_{i+1}`.  So a pairwise-disjoint subfamily of the union hits each
> piece at most once and never two cyclically consecutive pieces: choosing one index per member
> injects it into an independent set of `C_k`, whose size is at most `⌊k/2⌋ = (k-1)/2`. ∎

`k = 3` is a clique row (`(k-1)/2 = 1`); `k = 5, 7, 9` are the new ones.  **The only
well-formedness condition is `E_i ∩ E_{i+1} ≠ ∅`** — the sides need not be the sides of a convex
polygon, need not be short, and the `A_i` need not be distinct; taking them to be the sides of a
closed walk on `k` points makes the condition automatic, which is why the separator is
parameterised that way.  Degenerate choices are harmless: `A_i = A_{i+1}` makes `X_i` a point
clique and the row weaker, never invalid.

**Separation** (`rankfamily.separate`, called where the clique separator is called).  The
candidate anchors are the exact arrangement vertices of the current support, reduced to distinct
point cliques and ranked by point mass (top `--pent-cands`, default 240).  For each `k` in
`--pent`, coordinate descent on the `k` anchor indices maximises the support mass of the union;
it is fully vectorised (for one anchor position, the two sides that touch it are swept over all
candidates at once as a boolean matrix product), so a sweep is five to nine numpy operations and
the whole separation costs `0.1–0.5 s` per iteration.  Restarts are random plus the previous
iteration's best anchors (`--pent-restarts`, default 60; `--pent-time` caps the budget).

**Rows are maximal over the whole pose set**, as clique rows are: a pose is in the row iff its
square contains some `E_i`, i.e. contains both endpoints of that side (a square is convex).  The
membership test is `leaf_ceiling.sq_contains` in integers, with a float prefilter trusted only
outside a `1e-7` band.  Before a row enters the LP it is (a) **re-derived exactly** — every listed
member must pass the integer containment test on some side — and (b) its support members' exact
independence number is computed by a complete branch and bound and checked against `(k-1)/2`.  A
row failing either test is dropped with a log line.  The hill climb is a heuristic and can only
fail to find a row, never produce an invalid one.  `finalize`'s region top-up now also respects
polygon slack, and the final exact certification re-derives every polygon row from its stored
rational anchors and reports its exact mass and slack.

## 10. The number: what the family is worth in the loop

Two experiments.  The **support experiment** puts the family on the pose set the converged QSTAB
measure actually uses (94 and 149 poses) — small LPs, minutes, and the question is clean because
`alpha` of that very pose set is known exactly (§2: 11 in both cases).  The **loop experiment**
puts it on the full recorded pose sets of `CLIQUELEVER.md` §3 (10,175 and 9,551 columns), resumed
exactly as `runs/launch_2026-09-11.sh` resumes them, on the restricted master of
`search/CLMASTER.md` §4.

### 10.1 Support experiment: the family closes the whole gap on the leaf

| support | poses | `alpha` | QSTAB (cliques only) | **QSTAB + odd polygons** | drop | share of `mass - alpha` |
|---|---|---|---|---|---|---|
| leaf `01010101` (`cl_A0101L2`) | 94 | **11** | `11.435484` (`SUPA0`) | **`11.000000`** (`SUPA`) | `0.435484` | **100 %** |
| corner `k = 4` (`cl_B40KL2`) | 149 | **11** | `11.785783` (`SUPB0`) | **`11.470839`** (`SUPB2`) | `0.314944` | **40 %** |

Both `SUPA` and `SUPB2` **converged** (no violated coverage vertex, no clique of mass `> 1` with a
complete B&B, and no violated polygon the separator can find) and both final measures were exactly
certified and re-checked by `leaf_ceiling.py check --anchor clique`:

* **`SUPA`**: `mass = 10999999997/1000000000 = 10.999999997`, `M = 1`, max clique `= 1` (complete),
  regions `C=1,1,1,1`, slots `0,1,0,1,0,1,0,1`, chord strips `3,3,3,3` — **on 11 poses whose
  closed-intersection graph has 0 edges**.  The LP has walked all the way down to an *integral*
  packing: eleven pairwise-disjoint unit squares in the leaf's own region pattern.  43 polygon rows,
  all satisfied, 27 of them tight, every one re-derived from its exact rational anchors.
* **`SUPB2`**: `11.470839`, `M = 1`, max clique `= 1` (complete), regions OK; 362 polygon rows, all
  satisfied, every one re-derived.  `alpha` of its own support is again 11, so `0.47` of gap is
  left — the separator stops finding violated polygons before the LP reaches the integer optimum.

The separator's restart budget matters a great deal: on the corner support, 300 restarts per `k`
converge at `11.519638`, 2,500 restarts at `11.470839`.  Everything reported here is therefore an
UPPER bound on what the family can do — a better separator can only push it lower.

### 10.2 Anatomy of the binding polygons at convergence

`rankdiag.py --pgons` rebuilds every row from its stored exact rational anchors and re-checks it.
On `SUPB` (124 rows, 38 with positive dual): **0 violated, 0 with `alpha(G[X]) > (k-1)/2`**.
Every binding row is a `k = 5` polygon, and every one of them sits exactly where §5 said the
round-1 pentagons sat:

| binding row | dual | `\|X\|` | `mu` | `alpha` | anchors | grid vertex | regions of the mass | angles |
|---|---|---|---|---|---|---|---|---|
| 1 | `0.322` | 30 | `2.000000` | 2 | four on the line `x = 1`, one interior | **`(1,2)`** (9 of 30) | `I 1.000`, `W6 0.505`, `W7 0.495` | 18 |
| 2 | `0.289` | 28 | `2.000000` | 2 | three on `x = 3`, two interior | **`(3,2)`** (9 of 28) | `I 1.000`, `W2 0.569`, `W3 0.431` | 18 |
| 3 | `0.210` | 28 | `2.000000` | 2 | two on `y = 1`, three interior | **`(2,1)`** (9 of 28) | `I 1.000`, `W1 0.552`, `W0 0.448` | 18 |

Three things are worth naming.  First, **the binding rows are tight to the last digit**: `mu`
is exactly `2.000000` on each, so the right-hand side `2` is doing real work — these are facets of
the relaxation as the LP sees it, not slack decoration.  Second, **the mass splits `1.000` interior
and `1.000` wall**, and the wall half splits again across the two wall regions the grid vertex
separates (`W6 + W7`, `W2 + W3`, `W0 + W1`).  A polygon wrapped around the corner of an
axis-parallel wall square is exactly an object that says "the interior square and the wall square
here cannot both be replaced by two" — which no clique row can say, because no point lies in all
28–30 members.  Third, **`k = 5` is what binds**.  The separator finds and adds plenty of longer
polygons — `SUPB2` ends with 291 `k = 5`, 38 `k = 7` and 33 `k = 9` rows, `PGA` with 100 / 65 / 76
— and many of them are tight at intermediate solves, but at convergence the rows carrying positive
dual are 50 of `k = 5` and exactly **one** of `k = 7` (`SUPB2`), none of `k = 9`; on `SUPA`,
`SUPAP`, `SUPB` and `PGA` every dual-carrying row is `k = 5`.  The extra freedom of a longer cycle
buys almost nothing here: the geometry being cut is a wall-square corner, and a corner needs five
anchors.  A verifier that only ever learned the `k = 5` block would lose very little.

### 10.3 The pricing stage brings back nothing

`CLIQUELEVER.md` §0's clique result had a sting in it: one pricing stage on the leaf brought the
value back up by `+0.02` and the loop then sat on a degenerate optimal face for 60 iterations
without moving.  The same experiment with the rank family (`SUPAP`: `SUPA` plus `--price 1`):

| stage | poses | columns | LP | `M` | max clique | converged | exactly certified |
|---|---|---|---|---|---|---|---|
| 0 | 94 | 94 | **`11.000000`** | `1` | `1` (complete) | yes, 9 it | `10.999999998` |
| pricing | +559 new poses (lattice gap `+1.000000`); clique rows extended by 3,084 memberships; **the 131 polygon rows re-derived over the new pose set, +14,865 memberships** |
| 1 | 653 | 670 | **`11.000000`** | `1` | `1` (complete) | yes, 10 it | `10.999999997` |

**The value does not climb.**  Stage 1 converges at exactly `11.000000` again, on a support of 11
pairwise-disjoint squares, with 273 polygon rows (140 tight, none violated, every one re-derived
from its exact rational anchors) and `leaf_ceiling.py check --anchor clique` giving
`coverage OK, regions OK, anchor cliques PROVED <= 1`.  The pricing delta is **`+0.000000`**,
against `+0.02` for the clique family.

That is the qualitative difference between the two families on this leaf.  Clique rows cut the
non-Helly mass down to `11.30–11.79` and then the pricer finds new poses that rebuild a fractional
optimum just above it.  The polygon rows cut to the *integer* optimum, and there is nothing for the
pricer to rebuild: any new pose is either already in some polygon row (the rows are re-derived
maximal over the enlarged set, which is what `rankfamily.regrow` is for) or does not help.  On the
leaf's own pose set, `alpha = 11` is not merely the floor — it is where the relaxation lands.

## 9. Spec: what the verifier and Lean would need

No verifier or Lean change has been made.  This is the spec for when the family earns it.

### 9.1 `verify/`: a cyclic block next to `cliques`

`search/BOXCLIQUE.md` already has every primitive.  A **box clique** is a union of boxes; the
verifier computes each box's **core** (the concentric shrunk rectangle every unit square with a
pose in that box contains — the shrink lemma), refuses the block unless **every** pair of cores
meets (`cores_meet`, an exact separating-axis test between two rotated rectangles in `i128`),
credits a swept cell the weight `w_K` iff the cell lies inside one of the boxes, and adds `w_K`
once to the bound total.

A **cycle block** is the same object with two rules changed:

| | `cliques` block (today) | `cycle` block (proposed) |
|---|---|---|
| boxes | any number, unordered | `k` boxes in a **stated cyclic order**, `k` odd, `k >= 5` |
| core test | `cores_meet(i, j)` for **every** pair `i <= j` | `cores_meet(i, i+1 mod k)` for the `k` consecutive pairs, plus each core non-empty |
| sweep credit | cell inside some box gets `w` | unchanged |
| bound total | `+ w` | `+ w * (k-1)/2` |

Everything else — the parser, the empty-core refusal, the "a cell that merely *meets* a box gets
nothing and its LP witness is placed outside the box" rule, the no-symmetry-images rule that forces
a `[0°, 90°]` sweep — carries over unchanged, because none of it depends on how many pairs of cores
were required to meet.  The rejection tests gain the mirror images of the existing ones: a cycle
whose `cores_meet(i, i+1)` fails for one `i`, an even `k`, a `k < 5`, a repeated box index, and a
cycle credited `(k+1)/2` instead of `(k-1)/2`.  `certificates/FORMAT.md` gains one block type whose
body is the existing box list plus the cyclic order.

The separator here emits point anchors, i.e. degenerate cores; the verifier's cores are inward
rounded rectangles, so a certificate would use short segment/rectangle cores around each `A_i` —
`cores_meet` on consecutive pairs is then exactly `E_i ∩ E_{i+1} ≠ ∅` with slack.

### 9.2 `lean/`: the sibling of `clique_of_cores`

Two lemmas, both short.  First the combinatorics, which is where all the new content is:

```lean
/-- An independent set of the cycle `C_k` has at most `k / 2` elements: the translates
    `{s, s+1}` for `s ∈ S` are pairwise disjoint (a shared element would force two members
    of `S` to be equal or cyclically adjacent), each has 2 elements, and all lie in `ZMod k`. -/
lemma card_le_of_cycle_independent {k : ℕ} (hk : 2 ≤ k) (S : Finset (ZMod k))
    (hS : ∀ s ∈ S, ∀ s' ∈ S, s ≠ s' → s' ≠ s + 1) : 2 * S.card ≤ k
```

Then the geometric sibling of `clique_of_cores`, with `card_filter_clique_le_one` generalised from
`1` to `(k-1)/2`:

```lean
/-- **Cyclic core families.**  Pieces `B : Fin k → pose → Prop` with cores `core i` contained in
    every closed unit square of piece `i` (`hcore`, the verifier's shrink lemma), and consecutive
    cores meeting (`hmeet i : (core i ∩ core (i+1)).Nonempty`).  In a packing at most `k / 2`
    squares have a pose in the union. -/
lemma card_filter_cycle_le {n k : ℕ} (L : ℝ) (hL : 1 < L) (ctr ang) (hdisj …)
    (B : Fin k → ℝ × ℝ → ℝ → Prop) (core : Fin k → Set (ℝ × ℝ))
    (hcore : ∀ i c θ, B i c θ → core i ⊆ sq c θ 1)
    (hmeet : ∀ i : ZMod k, (core i ∩ core (i + 1)).Nonempty) :
    ((univ : Finset (Fin n)).filter (fun i => ∃ j, B j (ctr i) (ang i))).card ≤ k / 2
```

*Proof.*  Choose for each packed square in the union a piece index (`Finset.exists_of_filter` /
choice).  The index map is injective on the filtered set and its image is cycle-independent: two
squares with equal or cyclically adjacent indices both contain a common point of the corresponding
cores (exactly the two cases of `clique_of_cores`' one-line proof), so their `L`-interiors would
meet, contradicting `hdisj` via `unit_subset_interior`.  Apply `card_le_of_cycle_independent`. ∎

Finally `packing_le_weight_cliques` generalises by replacing the scalar `1` in its `hcard` step
with a per-family cap `r j`:

```lean
theorem packing_le_weight_ranks
    (A w hw C) (m) (K : Fin m → ℝ × ℝ → ℝ → Prop) (v : Fin m → ℝ) (hv : ∀ j, 0 ≤ v j)
    (r : Fin m → ℕ)
    (hK : ∀ j, ((univ : Finset (Fin n)).filter (fun i => K j (ctr i) (ang i))).card ≤ r j)
    (hcover : …) : (n : ℝ) ≤ ∑ a ∈ A, w a + ∑ j, v j * r j
```

The existing proof goes through verbatim: `hsum` is unchanged, and `hcard` becomes
`v j * card ≤ v j * r j` by the same `nlinarith`.  `packing_le_weight_cliques` is the case
`r = 1` with `hK` discharged by `card_filter_clique_le_one`; the cycle block is the case
`r j = k_j / 2` with `hK` discharged by `card_filter_cycle_le`.  No existing lemma is touched.

## 11. Reproduce (round 2)

```sh
python3 search/rankfamily.py                  # the k-ring selftest for k = 5, 7, 9
python3 search/cliquelever.py selftest        # the pinwheel, old path AND --master, with --pent on
R=runs
# the support experiment (minutes): QSTAB alone, then QSTAB + the rank family, on the pose set
# the converged measure itself uses
sh search/rankfamily_launch.sh support        # SUPA0/SUPA (leaf), SUPB0/SUPB (corner)
# one pricing stage on top of the leaf's converged 11.000000
sh search/rankfamily_launch.sh price          # SUPAP
# the loop experiment on the full recorded pose sets, restricted master (CLMASTER.md 4)
T=9000 sh search/rankfamily_launch.sh loop    # PGA (leaf, --price 1), PGB (corner)

# independent re-checks of any converged measure
python3 search/leaf_ceiling.py check $R/cl_SUPA_exact.txt --corners 1111 \
    --patterns 01010101 --chord --anchor clique
python3 search/rankdiag.py --pgons $R/cl_SUPA_pgons.json $R/cl_SUPA_exact.txt --anat 6
python3 search/rankdiag.py --step 1 $R/cl_SUPA_exact.txt      # alpha of the new support
```

`--pent 5,7,9` turns the family on; `--pent-restarts` is the one setting that matters (300 per `k`
converges the corner support at `11.519638`, 2,500 at `11.470839`), `--pent-cands` the number of
candidate anchor points, `--pent-time` the per-iteration budget, `--pent-want` the rows added per
`k` per iteration.  With `--pent` absent the loop is exactly the one `CLMASTER.md` describes.

Every converged run writes `runs/cl_<TAG>_pgons.json`: for each polygon row its `k`, its members,
its dual, its mass at the last solve, and **its anchors as exact rationals**, which is all
`rankdiag.py --pgons` needs to rebuild and re-check the row from scratch.
