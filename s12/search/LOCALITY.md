# locality-ceiling: how much can ANY local family ever prove?  (`LOC_D`, 2026-09-11)

Task: `tasks/locality-ceiling/README.md`.  Code: `search/locality.py` (imports `leaf_ceiling.py`
for the exact geometry and `cliquelever.Poses` for the exact pose set / boundary semantics; nothing
in `verify/`, `lean/`, `certificates/` is touched).  Runs in this worktree's `runs/loc_*`
(gitignored; every number is quoted here).  Read against `search/CLIQUELEVER.md` §0–§1,
`search/RANKDIAG.md` §0, `search/CLMASTER.md` §1.

## 0. Verdict, up front

**The "global counting is unavoidable" branch of the brief is ruled out, and the obstruction that
remains is the pentagon — measured, exactly, as a factor `5/4`.**

Three findings, in order of how much they change the plan.

1. **`LOC_D <= QSTAB(P)` for every `D >= 2`, so `LOC_3 >= 11.5` is impossible.**  Every clique of
   the closed-intersection graph has centre-extent at most `sqrt 2` (two closed unit squares meet
   only if their centres are within `sqrt 2`), and a set of centre-extent `s` fits inside some
   window of the grid as soon as `s <= D - pitch = 3D/4`, i.e. `D >= 4 sqrt(2)/3 = 1.8856`.  So at
   `D = 2, 2.5, 3, 3.5` **every** clique inequality — and every coverage row, a point clique being
   a clique — is *implied* by window consistency, with no separation at all.  With
   `CLIQUELEVER.md`'s converged values that pins the whole upper half of the table:

   > `LOC_D(leaf) <= 11.435484` and `LOC_D(corner k = 4) <= 11.785783` for every `D >= 2`.

   The brief's decision rule was `LOC_3 >= 11.5 ⇒ no local family can close the leaf and the
   endgame needs global counting`.  **That cannot happen**: `LOC_3 <= LOC_2 <= 11.4355 < 11.5` on
   the leaf's own pose set.  Whatever the exact value is, the locality ceiling on this pose set
   lies in `[alpha, 11.4355] = [11, 11.4355]` — a window of `0.44`, not of `1.0`.

2. **`D = 3.5` is already the whole container, and `D = 3` is nearly so.**  Every admissible centre
   lies in `[1/2, 7/2]^2`, and the box `[0, 3.5]^2` is a member of the `D = 3.5`, pitch-`0.875`
   grid, so at `D = 3.5` the inclusion-maximal window family is a **single** window containing all
   `10,175` poses: `LOC_3.5 = LOC_4 = max{mass : mu in conv(packings of P), regions pinned}`, which
   is at most `alpha(P)`.  At `D = 3` only four windows survive, of `5,240`–`6,132` poses.  The table's interesting
   rows are therefore `D = 1.5, 2, 2.5` — and, as (3) shows, they all carry the *same* obstruction.

3. **The obstruction is the pentagon, and it is worth exactly `5/4`.**  `locality.py check` takes a
   measure and computes, exactly, the Dantzig–Wolfe convexity multiplier each window needs,

   > `z_j = min { sum_I lam_I : sum_{I ni i} lam_I >= mu_i, lam >= 0, I independent in G[P_j] }`,
   > and `mu|_{P_j} in Pi_j` iff `z_j <= 1`,

   by column generation with exact max-weight-independent-set pricing, and then recomputes the
   optimum by an exact rational simplex over the generated columns.  On `CLIQUELEVER.md`'s
   **converged, exactly certified** QSTAB measures:

   | measure | mass | `D = 1.5` | `D = 2` | `D = 2.5` | worst window | its mass |
   |---|---|---|---|---|---|---|
   | leaf `01010101` (`cl_A0101L2_exact`) | `11.435484` | **`1.249999995`** | **`1.249999995`** | **`1.249999995`** | `[1, 3] x [1.5, 3.5]` (D=2) | `3.548387` |
   | corner `k = 4` (`cl_B40KL2_exact`) | `11.785783` | **`1.161807219`** | `>= 1.233674684` | (running) | `[2.25, 3.75] x [1.125, 2.625]` (D=1.5) | `2.550602` |

   (the corner `D = 2` entry is a partial sweep: five of the nine windows had been certified when
   the time ran out, at `z = 1.1597, 1.1599, 1.2307, 1.2316, 1.2337` — every one of them `> 1`, so
   the conclusion "not window-consistent" is already settled there too, only the maximum is not.)

   `249999999/200000000 = 1.249999995` is `5/4` to the last digit the `/10^9` measure can carry.
   The QSTAB optimum is therefore **not** window-consistent at any of these `D`, and by a wide
   margin: the whole measure would have to be scaled by `1/1.25` to fit.  The `D = 1.5` witness
   window is `[1.125, 2.625] x [2.25, 3.75]`, `35` support poses of mass `2.564516` — **the box
   that wraps the grid vertex `(2, 3)`**, which is exactly where `RANKDIAG.md` §0 found this
   measure's violated pentagon (`mu = 2.435484` on 32 poses, excess `+0.435484` = 100 % of the
   non-clique gap).  `5/4` is the fractional cover number of a `C5` at mass `1/2` per class —
   `chi_f(C5)/2 = 5/2 / 2`.  The window subproblem, which knows nothing about pentagons, rediscovers
   the pentagon as **the** local obstruction and prices it: enforcing local consistency costs the
   QSTAB measure a factor `5/4` in one `1.5 x 1.5` box.

**What this means for the plan.**  Local families are *not* the binding constraint at `t = 4`: the
ceiling of everything local is at most `11.44` on the leaf and at most `11.79` on the corner
control, against the `12` that has to be excluded.  The gap that has to be closed is the
`0.44`/`0.79` between the all-clique value and `alpha = 11`, and item (3) says the pentagon rank
inequality of `RANKDIAG.md` is the object that carries it — which is a *local* inequality, of
centre-extent `<= 1.5`.  **Global counting is not forced by these measurements; a rank family at
window `1.5` is what is missing**, and `rankdiag.py --step 4` dropped into `cliquelever.py`'s loop
(RANKDIAG §"The next experiment") is the concrete next step.  The one thing not settled here is
whether the ceiling is exactly `11` at some `D <= 3`: see §4.

## 1. What is computed, and what is a bound

`LOC_D(P) = max sum mu` subject to `cov(x) <= 1` at every arrangement vertex, the leaf's pinned
region counts (boundary poses free to declare either side, `cliquelever`'s semantics), the chord
rows `mu(wall strip) <= 3`, and `mu|_{P_j} in Pi_j` for every window, with
`P_j = {S : centre(S) in W_j}` and `Pi_j = conv{1_I : I independent in G[P_j]}` the packing
polytope of the exact closed-intersection graph inside the window.  Windows are closed boxes of
side `D` (or closed discs of diameter `D`) with centres on a pitch-`D/4` grid, clipped to the
container, reduced to the inclusion-**maximal** pose sets (a window whose poses are contained in
another's imposes a weaker constraint; two windows with the same poses impose the same one).

`mu|_{P_j} in Pi_j` implies *every* valid inequality of *every* local family supported in `W_j` —
point rows, anchor cliques, general cliques, the pentagon, anything not yet invented — so `LOC_D`
is the ceiling of local reasoning at window diameter `D`.

Directions, as `T4SCREEN.md` §1.2:

* restricting the poses LOWERS the value, so every number here is a lower bound on the continuum
  `LOC_D`;
* the grid family brackets the continuum one: every continuum window of side `D` sits inside some
  grid window of side `D + pitch`, so `LOC^{grid}_D <= LOC^{cont}_{3D/4}`, and `LOC^{grid}_D >=
  LOC^{cont}_D` because the grid family is a subfamily.  A grid number at `D` is therefore a
  statement about the continuum between `3D/4` and `D`;
* in the **cutting-plane** form the loop runs (§2), every row is valid, so every intermediate LP
  value is an UPPER bound on `LOC_D(P)`; at convergence it is the value;
* the exactly certified final measure is a rigorous LOWER bound on `LOC_D(P)`.

## 2. The method

The brief's Dantzig–Wolfe master is

    conv_j:      sum_I lam_{j,I} <= 1                              (one row per window)
    link_{j,i}:  mu_i <= sum_{I ni i} lam_{j,I}                    (one row per (j, S))

— the linking rows may be inequalities because the stable-set polytope is down-monotone.  Pricing
is a maximum-weight independent set in `G[P_j]` under the linking duals, solved exactly as a
max-weight clique of the complement (`leaf_ceiling.max_clique_int`; the complement of a unit-square
intersection graph inside a window is sparse).  `locality.py` implements that master
(`Model`, the linking/convexity blocks) but runs it in its equivalent **cutting-plane** form,
because on this instance the primal master climbs far too slowly (§3):

1. the outer LP is in the `mu` variables alone (coverage rows, region pins, chord rows, window
   cuts), so it descends from `12.000000` and every value is an upper bound;
2. the window subproblem for the current `mu` is the Dantzig–Wolfe master of window `j` restricted
   to the **support** of `mu` — a pose with `mu_i = 0` constrains nothing, and an independent set
   of `G[P_j]` meets the support in an independent set of the support subgraph, so nothing is lost
   and the subproblem has 5–60 rows instead of 4,251.  Its value is `z_j`; its optimal dual `pi`
   satisfies `sum_{i in I} pi_i <= 1` for every independent set — which is exactly what the pricing
   max-weight independent set *proves* on the iteration that ends the column generation — so
   `sum_i pi_i mu_i <= 1` is a valid packing inequality violated by `z_j - 1`;
3. that cut is extended over the whole window: a pose outside the support may carry
   `1 - (heaviest independent set of the support it can join)`, and the extended vector is
   renormalised by its own exact max-weight independent set, which leaves it valid whatever
   happened.  When `pi` is the indicator of a clique this is exactly `cliquelever`'s greedy maximal
   clique extension;
4. convergence: no violated arrangement vertex and no window with `z_j > 1`.

**Exact certification.**  `window_certificates` re-does every window subproblem in exact rational
arithmetic: the columns come from the float column generation (whose last pricing B&B proves none
is missing), the optimum is recomputed by a Fraction tableau simplex with Bland's rule
(`exact_cover_lp`), and every independent set in the certificate's support is re-verified pairwise
disjoint by the exact integer separating-axis test (`leaf_ceiling.sq_meets_sq`).  The final measure
is rounded down to `/10^9`, the pinned regions are topped back up on exact *window* slack, and `M`,
the region counts, the max-weight clique and `alpha` of the support are recomputed in integers.

**`selftest`** covers the exact simplex (a triangle at `1/3` → `1`, a disjoint pair at `1/2` →
`1/2`, `C5` at `1/2` → `5/4 = chi_f(C5)/2`), the maximal-independent-set enumeration on `C5`, and
the exact max-weight independent set.  `5/5`.

## 3. The window families, and why the coverage rows are free

Both pose sets, inclusion-maximal windows, `pitch = D/4`, boxes:

| `D` | pitch | leaf: windows | largest `P_j` | linking rows | corner: windows | largest `P_j` |
|---|---|---|---|---|---|---|
| `1.5` | `0.375` | `36` | `2,308` | `53,423` | `36` | `2,269` |
| `2` | `0.5` | `9` | `4,251` | `27,899` | `9` | `4,458` |
| `2.5` | `0.625` | `9` | `4,516` | `36,102` | `9` | `4,809` |
| `3` | `0.75` | `4` | `6,132` | `22,512` | `4` | `5,240` |
| `3.5` | `0.875` | **`1`** | `10,175` | `10,175` | **`1`** | `9,546` |

(leaf `10,175` poses / `10,463` columns; corner `9,546` / `9,669` — the post-pricing `cliquelever`
checkpoints, values `11.435484` and `11.785783`.)  The largest independent set inside a `D = 2`
window of the leaf pose set is **4–5** (20 randomised greedy maximal sets per window), and the
certified QSTAB measure puts `3.44`–`4.98` of mass in each of them: the windows are tight objects,
not slack ones.

`D(1 - pitch) = 3D/4` versus `sqrt 2 = 1.4142`: `1.125` at `D = 1.5` (not implied), `1.5`, `1.875`,
`2.25`, `2.625` at `D = 2, 2.5, 3, 3.5` (implied).  For `D >= 2` the loop therefore carries **no**
coverage rows beyond a `0.08` bounding grid — they are by far the densest rows in the model
(a generic container point lies in 500–2,000 of the 10,175 poses, so the `25,012` recorded rows of
`cl_A0101L2_rows.txt` are ~6M nonzeros) and removing them took the LP solve from `76 s` to `2 s`.

## 4. What did not finish

**The outer cutting-plane loop did not descend from `12.000000` in the time available.**  On the
leaf at `D = 2`, 200 iterations (about 5 s each: LP `2 s`, all nine window subproblems `1 s`) added
`18` cuts apiece — the raw cut and its extension per window — and the LP value stayed at `12`
throughout while `zmax` wandered between `1.84` and `2.57` and the support moved between `12` and
`123` poses.  This is the same phenomenon `CLIQUELEVER.md` §3 records for the clique loop, only
worse: with `10,463` columns and a 4-dimensional pose lattice, a cut supported on `~150` poses is
easy for the LP to walk around, and `cliquelever` needed `17,653` rows and 31 iterations from a
seed of `11,717` rows to move `11.95 → 11.30`.  The honest reading is that the descent needs the
recorded coverage/clique row bases as a seed (they are valid for `LOC_D` too, being implied by it)
and a run of hours, not minutes.  The primal Dantzig–Wolfe master — which is what the brief
specifies — was tried first and is worse on this instance: seeded with `22,581`/`12,704`
independent-set columns it started at `9.00`/`10.06` and climbed `+0.02` per iteration.

So **the `LOC_D` table is not filled in with converged values**, and no `loc_*_exact.txt` measure
was produced.  What replaces it is the pair of rigorous statements above: the *upper* half of the
table is pinned by §0(1) with no LP at all,

| pose set | `D = 1.5` | `D = 2` | `D = 2.5` | `D = 3` | `D = 3.5` | `D = 4` |
|---|---|---|---|---|---|---|
| leaf `01010101` | `<= 11.914756` | **`<= 11.435484`** | **`<= 11.435484`** | **`<= 11.435484`** | `= LOC_4` | `= LOC_4` |
| corner `k = 4` | `<= 11.988925` | **`<= 11.785783`** | **`<= 11.785783`** | **`<= 11.785783`** | `= LOC_4` | `= LOC_4` |

The `D >= 2` entries are §0(1): the clique family is implied, so `LOC_D <= QSTAB(P)`, and `QSTAB(P)`
is `CLIQUELEVER.md`'s converged value on each pose set.  At `D = 1.5` only cliques of centre-extent
`<= 1.125` are implied, so the bound falls back to the pure coverage+region LP (`T4LEAF.md`
`11.914756` on the leaf; `11.988925`, the anchor-clique value, on the corner control).  `D = 3.5`
and `D = 4` are literally the same relaxation, both being the single window that is the whole pose
set.  The lower end, `LOC_D >= alpha(P)`, is `11` if `alpha(P) = 11` as the brief states (this task
did not recompute `alpha` over the full 10,175-pose set; `RANKDIAG.md` §0 proves it on the QSTAB
supports, and note that it is the REGION PINS, not the container, that make `12` hard here — a
`4 x 4` box holds 16 axis-parallel unit squares, but not with corner counts `1,1,1,1` and slot
counts `0,1,0,1,0,1,0,1`).  And the *discriminating* measurement — is the QSTAB optimum itself
window-consistent? — is answered exactly, and negatively, by `5/4`.

**The smallest `D` at which the leaf reaches `alpha = 11` is therefore known only as `<= 3.5`**, and
at `D = 3.5` it is a degeneracy of the container rather than a fact about locality (the window is
everything).  Whether it already happens at `D = 2` — the brief's `LOC_2 ≈ 11` branch — is the one
number this task owes.  Discs were not run.

**To finish it** (in order): seed the outer LP with `cl_A0101L2_rows.txt` (25,012 exact coverage
rows) and `cl_A0101L2_cliques.txt` (1,128 verified clique rows), both valid for `LOC_D`, so that
the descent starts at `11.4355` instead of `12`; then only the window cuts have to do work, and the
question becomes how far below `11.4355` they push.  `locality.py` already has `--row-pitch` and
would need a `--resume-rows/--resume-cliques` pair copied from `cliquelever.Lever`.

## 6. Round 2: the seeded `LOC_2` descent, and the cut the window subproblem actually writes down

### 6.1 Seeding the outer LP (`--resume-rows`, `--resume-cliques`, `--tight`)

§4's route, implemented.  A coverage row is part of the definition of `LOC_D`, and a clique row
`mu(K) <= 1` is the `pi = 1_K` case of a window cut — implied by `LOC_D` whenever the clique fits
in a window, which at `D = 2` is always (§0(1)).  Both are therefore valid for `LOC_D` and seeding
with them only starts the descent lower; they cannot change the value.  `--resume-cliques`
re-verifies every row pairwise (exact integer SAT on the borderline pairs) before it enters the LP.

The full `25,012`-row seed is not affordable: a generic container point lies in 500–2,000 of the
`10,175` poses, so the row block is ~6M nonzeros and one ipm+crossover on `26,140 x 10,463` ran for
**22 minutes** without finishing.  Three economies, in the order they matter:

* `--tight 1e-6` keeps only the rows the seed measure already saturates — **`4,952` of `25,012`**
  on the leaf, **`4,811` of `32,191`** on the corner control.  Same warm start at a fifth of the
  cost: the dropped rows carry no dual at the seed measure, they are still implied by the windows,
  and the coverage separation would re-add any that bite.
* **no crossover.**  The cutting-plane loop never reads the outer LP's duals — the window
  subproblem separates on `mu` alone — so the crossover is pure cost on a `6k x 10k` coverage LP
  with a large degenerate optimal face.  Off by default now (`--crossover` restores it).
* `--lp-tlim` adds `t4leaf.Hi`'s rule (warm dual simplex after the first solve, ipm fallback).  On
  these LPs it is `CLMASTER.md` §2.1 again — the simplex burned the whole 45 s budget and the ipm
  ran anyway, every single iteration — so `--lp-tlim 0` (always ipm) is the right setting.

**The seed does exactly what it was meant to do: iteration 0 of the seeded run is
`LP = 11.435484`**, the converged QSTAB value of `CLIQUELEVER.md` reproduced to the digit, with
`zmax = 1.532258` — i.e. before any window cut the outer LP *is* the QSTAB LP, as §0(1) predicts,
and the vertex it sits on is even further from window-consistent than the certified measure.

### 6.2 The window subproblem's optimal dual, named

`locality.py check` now reports the *shape* of the separating inequality, not just its violation:
the dual `pi >= 0` of the window subproblem, its distinct values, the independence number of its
support (complete B&B), and the centroid / angle range / region split of each `pi`-class.  The
reading is mechanical — `pi = 1` on a clique is a clique row; `pi = 1/r` on a support with
`alpha = r` is a rank inequality `mu(X) <= r`, and with `2r+1` classes it is the odd polygon of
`RANKDIAG.md`; anything else is new.

**On the leaf's certified QSTAB measure (`11.435484`), at `D = 1.5` and at `D = 2`, the cut is
literally the same object:**

> `pi = 1/2` on **33 poses**, `alpha` of that support `= 2` (complete branch and bound),
> `mu(X) = 2.500000` exactly, so the inequality is the **rank inequality `mu(X) <= 2`**, violated by
> `0.5`, and `z = 2.5/2 = 5/4`.  Centroid `(1.831, 2.859)`, angles spanning `0.00–89.96` degrees,
> regions `I: 22, W5: 10, W4: 1`.

That is `RANKDIAG.md` §0's pentagon in its exact-rank closure: same place (the box wrapping the
grid vertex `(2,3)`), same regions (interior plus the adjacent wall region `W5`), 33 poses against
the pentagon's 32, and `2.500000` of mass against the pentagon's `2.435484`.  The pentagon
*undercounts the obstruction by `0.065`*; the rank inequality on the window subproblem's support is
the sharp form.

### 6.3 What the pentagon family leaves behind

`cl_PGA_measure.txt`, the rank agent's leaf run with QSTAB + odd polygons at `LP = 11.313432`
(iteration 33, not yet converged), `230` support poses:

| measure | mass | `z_max` at `D = 1.5` | witness window | its mass / support |
|---|---|---|---|---|
| `cl_A0101L2_exact` (QSTAB) | `11.435484` | `1.249999995` | `[1.125, 2.625] x [2.25, 3.75]` | `2.564516` / 35 |
| `cl_PGA_measure` (QSTAB + polygons) | `11.313432` | **`1.240874607`** | `[0.375, 1.875] x [1.125, 2.625]` | `2.651778` / 82 |

**The odd-polygon family removes `0.122` of mass and essentially none of the window
inconsistency**: `z_max` moves `1.250 -> 1.241`.  Local reasoning is *not* exhausted by the polygon
family at `D = 1.5`.

And the cut that now separates is **not** an anchor polygon.  It has two levels:

> `pi = 1` on a **24-pose clique** (mass `0.530137`, centroid `(1.460, 1.741)`, regions `I: 22,
> W7: 2`) **plus** `pi = 1/2` on a **47-pose rank-2 set** (mass `1.434513`, centroid
> `(1.077, 1.883)`, regions `I: 28, W7: 17, W6: 2`), with `alpha` of the whole 71-pose support
> `= 2` (complete B&B).  `pi . mu = 1.247394`.

Read as an inequality, `mu(K) + (1/2) mu(X) <= 1`: an independent set that uses a member of the
clique `K` can use nothing of `X`, and one that avoids `K` can use at most two of `X`.  It is
*strictly stronger than rank on the same support* — `mu(K union X) = 1.9647 < 2 = alpha`, so the
plain rank inequality `mu <= 2` does **not** fire there — and it is not an odd hole, wheel or
antihole of poses either (`RANKDIAG.md` §0 item 2 already found those worth zero).  It is the
join of a clique and a rank-2 set, which is the next family to build if the polygon family is to
be pushed further.  Geometrically it is again the interior against a wall region, this time `W7`
and the left wall, one window over from the pentagon.

### 6.4 `LOC_2`, seeded: what the runs show

Both seeded runs reproduce their pose set's QSTAB value at iteration 0 and then cut:

| run | pose set | it 0 | it 1 | it 2 | `zmax` | rows / cuts at it 0 | s / it |
|---|---|---|---|---|---|---|---|
| `A2S` | leaf `01010101` | `11.435484` | `11.435484` | — | `1.532258 -> 1.645161` | `4,952` / `1,146` | `61` |
| `B2S` | corner `k = 4` | `11.785783` | `11.785783` | `11.785783` | `1.239277 -> 1.272771` | `4,811` / `2,663` | `85` |

`LOC_2` is therefore **not yet pinned numerically**: the value sits on a wide degenerate optimal
face — `11.435484` and `11.785783` unchanged over the first iterations while `zmax` moves, and in
fact RISES, because each cut pushes the LP onto a neighbouring vertex of the same value that is
even less window-consistent (`1.53 -> 1.65` on the leaf) — which is the behaviour `CLIQUELEVER.md` §4 records for the
clique loop on the same instances (`317/28` and `537/46` held for 15–60 iterations while rows were
cut one vertex at a time).  Each iteration adds two cuts per window (the raw dual and its maximal
extension), 18 in all.

**What is still missing is a restricted master over the COLUMNS.**  `CLMASTER.md` §2.2 is the fix:
these LPs carry `10,463` columns of which `1,500`–`3,800` are ever active, and the ipm is
`t ~ n^1.5`.  `locality.py` does not have it (the cutting-plane loop was written not to need the
duals, which is exactly what pricing loaded columns requires), and adding it is the next piece of
work on this file.

## 7. Round 3: fitting a certifiable anchor family to the exact window cut (`search/anchorfit.py`)

### 7.1 The family, and why it is certifiable

The window subproblem's dual `pi` is the exact local facet on the support, but a verifier cannot
check "this list of 33 poses has independence number 2".  The certifiable form is the ANCHOR GRAPH.
Let `H` be a graph whose vertices carry anchors `A_v` (a point or a closed segment), let
`X_v = {S admissible : A_v subset S}`, and **assign every square to the LEAST `v` with
`A_v subset S`**, so the pieces are pairwise disjoint by a rule a verifier can apply.  Then for any
weights `pi_v >= 0` and any packing `P`:

* two squares of `P` cannot both lie in `X_v` — both contain `A_v`, so they meet;
* if `uv` is an edge of `H`, certified by `A_u cap A_v != {}`, a square of `P` in `X_u` and one in
  `X_v` both contain a common point, so they meet; being in disjoint pieces they are distinct,
  which is impossible.

So `S -> (its piece)` injects `P cap (union X_v)` into an independent set of `H`, each square paying
its own piece's weight, and

> **`sum_v pi_v mu(X_v) <= alpha_pi(H) := max{ sum_{v in I} pi_v : I independent in H }`**

is valid for every packing of disjoint closed unit squares, with no geometry beyond "a convex set
contains a segment iff it contains both endpoints" and the exact intersection test on `H`'s edges —
`cores_meet` on `H`'s edges plus the credit `alpha_pi(H)`, exactly as the brief's Lemma 0 predicts.
`H = C_k`, `pi = 1` is `rankfamily.py`'s odd polygon (`alpha_pi = (k-1)/2`); `H = K1 join C_5` with
`pi = 1` on the hub and `1/2` on the rim is the odd wheel, `alpha_pi = 1`.  The disjoint-assignment
rule is what makes the WEIGHTED form valid: with overlapping pieces a square would pay several
weights and the injection argument bounds only the union.

`anchorfit.py fit` takes a cut from `locality.py check`, hill-climbs anchors on the arrangement
vertices of the cut's support (`rankfamily._cand_masks` / `_climb`, imported, not forked), builds
`H` from the EXACT anchor intersections — a non-edge is never assumed, and an unplanned edge only
lowers `alpha_pi` — computes `alpha_pi(H)` by brute force over the independent sets of a graph with
at most 8 vertices, extends the pieces over the WHOLE loaded pose set (maximal rows), and reports

    captured = (sum_v pi_v mu(X_v) - alpha_pi(H)) / (pi . mu - 1).

`selftest` covers the exact anchor-intersection predicate (crossing, touching, collinear overlap,
point-on-segment) and `alpha_pi` on `C5`, `C7`, `W5`, `W7`: **10/10**.

### 7.2 The table: cut -> fitted `H` -> captured fraction

The leaf's certified QSTAB measure (`cl_A0101L2_exact`, `11.435484`).  **Every one of the nine
`D = 2` windows is violated**, and so are eight of the 36 at `D = 1.5`; every cut but one is a pure
rank-2 facet, `pi = 1/2` on 26–33 poses with `alpha = 2` under a complete branch and bound.

| cut | `z` | violation | `pi` | fitted | excess | **captured** |
|---|---|---|---|---|---|---|
| `D=2` w0 | `1.217742` | `0.217742` | `1/2` | `C5` | `+0.435484` | **`200.0 %`** |
| `D=2` w1 | `1.096774` | `0.096774` | `1/2` | `C5` | `+0.129032` | `133.3 %` |
| `D=2` w2 | `1.217742` | `0.217742` | `1/2` | `C5` | `+0.387097` | `177.8 %` |
| `D=2` w3 | `1.250000` | `0.250000` | `1/2` | `C5` | `+0.435484` | `174.2 %` |
| `D=2` w4 | `1.217742` | `0.217742` | `1/2` | `C5` | `+0.306452` | `140.7 %` |
| `D=2` w5 | `1.161290` | `0.161290` | `1/2` | `C5` | `+0.322581` | `200.0 %` |
| `D=2` w6 | `1.193548` | `0.193548` | `1/2` | `C5` | `+0.177419` | `91.7 %` |
| `D=2` w7 | `1.209677` | `0.209677` | `1/2` | `C5` | `+0.274194` | `130.8 %` |
| `D=2` w8 | `1.217742` | `0.217742` | `1/2` | `C5` | `+0.290323` | `133.3 %` |
| `D=1.5` w2 | `1.217742` | `0.217742` | `1/2` | `C5` | `+0.387097` | `177.8 %` |
| `D=1.5` w12 | `1.161290` | `0.161290` | `1/2` | `C5` | `+0.306452` | `190.0 %` |
| `D=1.5` w13 | `1.250000` | `0.250000` | `1/2` | `C5` | `+0.435484` | `174.2 %` |
| `D=1.5` w15 | `1.137097` | `0.137097` | `1/2` | `C5` | `+0.145161` | `105.9 %` |
| `D=1.5` w18 | `1.193548` | `0.193548` | `1/2` | `C5` | `+0.177419` | `91.7 %` |
| `D=1.5` w32 | `1.209677` | `0.209677` | `1, 1/2` | `C5` | `+0.274194` | `130.8 %` |
| `D=1.5` w34 | `1.185484` | `0.185484` | `1, 1/2` | `C5` | `+0.209677` | `113.0 %` |
| `D=1.5` w28 | `1.024194` | `0.024194` | `1, 1/2` | **none** | `-0.500000` | **`< 0`** |

**16 of 17 cuts are captured by a plain `C5`, at 92–200 %.**  Over 100 % is not a rounding artefact:
the cut is confined to the window's support, while the fitted polygon row is *maximal over all
10,175 poses*, so the certified inequality is **stronger** than the exact local facet that produced
it.  The best fit on `D=1.5` window 13 — the `5/4` window at the grid vertex `(2,3)` — is
`mu(X) = 2.435484` against `alpha(C5) = 2`, which is `RANKDIAG.md` §0's pentagon for this measure
**to the last digit**: the fitter, starting only from the window subproblem's dual, rediscovers the
hand-found pentagon.

### 7.3 What does not fit, and why the wheel is not realisable

`C7`, `K1 join C5` and `K1 join C7` **never fire**: on every one of the 17 cuts their excess is
negative.  For `C7` the reason is arithmetic — the fitted seven pieces reach `2.44–2.56` of mass
against `alpha(C7) = 3`.  For the wheels the reason is geometric and worth recording:

> the fitted hub anchor never meets the rim segments, so `build_H` certifies no hub–rim edge and
> `alpha_pi` comes out `2` (hub `1` plus two opposite rim vertices at `1/2`) instead of the wheel's
> `1`.  A genuine `W5` needs a point that lies on **all five** sides of the pentagon; but the rim
> edges of `W5` are the cyclic ones, so the five anchors must also meet consecutively — and five
> segments through one common point form a CLIQUE, not a cycle.  **`W5` has no anchor realisation
> with point/segment anchors in this configuration**, which is why the two-level cut
> `mu(K) + (1/2) mu(X) <= 1` found in round 2 has no small anchor description.

Where the measure's cut is two-level (`pi in {1, 1/2}`, windows 32 and 34 at `D = 1.5`) a plain
`C5` fitted to the half-level still captures `113–131 %`, so nothing is lost there.  The one cut
that defeats the menu is `D = 1.5` window 28, the *smallest* violation in the table
(`z = 1.024194`, `0.0242`): its best `C5` is `1.500000` against `alpha = 2`, i.e. `-0.50`.

### 7.4 Reading

**The anchor-graph family is the certifiable closure of local reasoning on this measure, and the
odd pentagon alone is very nearly all of it.**  Every substantial local facet the window
subproblem produces — violations from `0.097` to `0.250` — is matched, and usually beaten, by a
single `C5` anchor row that the verifier can already check with `cores_meet` and a right-hand side
of `2`.  The residue is one cut of violation `0.024` with no small anchor description; whether that
is a fitting failure or a genuinely new family is the open question, and it is a `0.024` question,
not a `0.4` one.

### 7.5 The same instrument on the pentagon-loop measure: nothing fits

`cl_PGA_measure.txt` — the rank agent's leaf loop with QSTAB **plus** odd polygons, `LP = 11.313432`
(iteration 33; PGA/PGB were at iteration 44, `11.306280` / `11.716220`, with no `_exact.txt` yet, so
this is a float checkpoint and `check --fast` was used: float cover value, exact anatomy).

The cuts change character completely.  On the QSTAB measure every one of the 17 cuts was a rank-2
facet with a single weight `1/2`.  Here, of the 31 violated windows (22 with violation `>= 0.02`):

| | leaf QSTAB measure | pentagon-loop measure |
|---|---|---|
| `z_max`, `D = 2` | `1.250000` | `1.343676` |
| `z_max`, `D = 1.5` | `1.250000` | `1.247394` |
| violated windows, `D = 2` | 9 of 9 | 9 of 9 |
| violated windows, `D = 1.5` | 8 of 36 | 22 of 36 |
| `alpha` of the cut support | `2` on every cut | `1, 2, 3, 4, 5` |
| distinct `pi` values per cut | `1` (or `2`) | `1` to **`35`** |
| **captured by the anchor menu** | **16 of 17, at 92–200 %** | **0 of 22** |

Not one cut fits.  The best of all 22 is `D = 1.5` window 32, a `C5` at excess `-0.056982`
(`-44.6 %`); the rest run from `-100 %` to `-3700 %`.  The `pi` spectra say why: at `D = 2` they are
ladders — `0.56, 0.52, 0.48, … , 0.04` (step `1/25`), `7/9, 6/9, … , 1/9`, `1, 5/6, 4/6, … , 1/6` —
on supports of 67–105 poses with `alpha = 3, 4, 5` under a complete branch and bound.  A `C_k` or a
wheel has one or two weight levels and `alpha <= 3`; these facets have a dozen levels and
`alpha = 5`.  The fitted pieces come back lopsided (one large piece, two or three empty), which is
the geometric statement that no five or seven anchors carry the mass.

Three of the 22 (`D = 1.5` windows 17, 23, 30; `alpha = 1`, single weight `1`, violations
`0.023–0.050`) are plain CLIQUE rows — residue of a loop that has not finished its own clique
separation, not new families.  The other nineteen are not.

**Reading.**  The anchor-graph family is the certifiable closure of local reasoning *on the QSTAB
measure* — the pentagon is essentially all of it, and `anchorfit` recovers `RANKDIAG.md`'s pentagon
from the window dual alone.  But it is **not** the closure of local reasoning full stop: once the
polygon family has been applied, the local facets that remain are higher-rank (`alpha` up to `5`)
and many-levelled, and no `C5`, `C7`, `W5`, `W7` or `K join C5` describes any of them.  Either a
genuinely wider certifiable family is needed for the last `0.3` — the ladder spectra suggest
antiwebs or comb inequalities rather than polygons — or that last `0.3` is where global counting
finally has to enter.  This is the first measurement in the project that separates "local" from
"certifiable-local".

### 7.6 The converged pentagon measures: the polygon family is SATURATED, and `0.15` of local
violation survives it

`cl_PGA_exact.txt` and `cl_PGB_exact.txt` landed (the rank agent's QSTAB + odd-polygon loops,
`RESULT tag=PGA s0=LP11.303071/M1.000000/K1.000000`): the leaf at **`11.303071`** on 201 poses,
exactly certified `M = 1`, max clique `= 1` with a complete branch and bound; the corner control at
`11.713225` on 380 poses (its max clique certifies at `1.006327`, so that one is not clique-clean
and its numbers below are read with that caveat).

| measure | mass | `z_max`, `D = 1.5` | `z_max`, `D = 2` | violated windows `D=1.5 / D=2` |
|---|---|---|---|---|
| leaf QSTAB `cl_A0101L2_exact` | `11.435484` | `1.250000` | `1.250000` | 8/36, 9/9 |
| leaf **+ polygons** `cl_PGA_exact` | `11.303071` | **`1.170998`** | **`1.290264`** | 19/36, 9/9 |
| corner QSTAB `cl_B40KL2_exact` | `11.785783` | `1.161807` | `>= 1.233675` | —, 9/9 |
| corner **+ polygons** `cl_PGB_exact` | `11.713225` | `1.256771` | `1.664832` | —, 9/9 |

Two things happen at once.  At `D = 1.5` the polygon family **does** reduce the local inconsistency,
`1.250000 -> 1.170998` — a third of the way from `5/4` to `1`, in the same window 13 that wraps the
grid vertex `(2,3)`.  At `D = 2` it makes it **worse**, `1.250000 -> 1.290264`: the LP, cut by rows
the polygon family can write, moves mass into a configuration those rows do not see but the exact
`D = 2` window facet does.

And the fit says the family is spent.  Of the 19 cuts with violation `>= 0.02`, **none is captured**
— but the near misses are now very tight:

| cut | `z` | violation | best fit | its value vs `alpha_pi` | excess | captured |
|---|---|---|---|---|---|---|
| `D=1.5` w13 | `1.170998` | `0.170998` | `C5` | `1.985575` vs `2` | `-0.014425` | `-8.4 %` |
| `D=1.5` w2 | `1.168754` | `0.168754` | `C5` | `1.957039` vs `2` | `-0.042961` | `-25.5 %` |
| `D=1.5` w18 | `1.154157` | `0.154157` | `C5` | — | `-0.083260` | `-54.0 %` |
| `D=1.5` w32 | `1.146976` | `0.146976` | `C5` | — | `-0.086941` | `-59.2 %` |
| `D=2` w1 | `1.205023` | `0.205023` | `C5` | — | `-0.154676` | `-75.4 %` |
| … 14 more | `1.022`–`1.290` | — | `C5` / `K1+C5` | — | `-0.17` … `-0.77` | `-75 %` … `-2440 %` |

**The best pentagon anywhere on this measure reaches `1.985575` against its right-hand side `2`.**
That is the definition of a saturated family: the loop has pushed every fittable `C5` to within
`0.014` of its bound, exactly as a converged polygon separator should, and the exact local facet is
*still* violated by `0.171` in the same window.  The residual is not a pentagon the search missed;
there is no pentagon left to find.

What the residual is instead: `alpha` of the cut support runs `2, 3, 4, 5` (complete B&B) and the
dual takes up to **30 distinct values**, in ladders — `1, 6/7, 5/7, … , 1/7` on `D=1.5` w15,
`1, 12/13, 10/13, … , 1/13` on `D=2` w4, a `1/48` ladder on `D=2` w7.  Those are antiweb- or
comb-shaped weights, not polygon-shaped ones.

**Reading, final.**  On the QSTAB measure the anchor-graph family *is* the certifiable closure of
local reasoning and the pentagon is essentially all of it (§7.2).  Run the pentagon family to
convergence and it takes the leaf from `11.435484` to `11.303071` — and then it is exhausted, with
`0.15–0.17` of local violation per window left over in facets with no small anchor description.
The next certifiable family, if there is one, is an antiweb/comb family at window `1.5`; the
alternative is that the last stretch to `alpha = 11` is where global counting has to enter.

## 5. Reproduce

```
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-11/* runs/

python3 search/locality.py selftest

# the discriminating measurement: is a given measure window-consistent, exactly?
python3 search/locality.py check QA --poses runs/cl_A0101L2_poses.txt \
    --measure runs/cl_A0101L2_exact.txt --D 1.5 2 2.5 3 3.5
python3 search/locality.py check QB --poses runs/cl_B40KL2_poses.txt \
    --measure runs/cl_B40KL2_exact.txt --D 1.5 2 2.5 3 3.5
python3 search/locality.py check PGAm --poses runs/cl_PGA_poses.txt \
    --measure runs/cl_PGA_measure.txt --D 1.5 2        # a pentagon-cut measure

# the LOC_D loop, seeded with the recorded rows and clique rows (6.1)
python3 search/locality.py run A2S --poses runs/cl_A0101L2_poses.txt \
    --resume-rows runs/cl_A0101L2_rows.txt --resume-cliques runs/cl_A0101L2_cliques.txt \
    --tight 1e-6 --corners 1111 --patterns 01010101 --chord --D 2 --shape box \
    --threads 2 --procs 2 --lp-tlim 45
bash search/locality_launch.sh A 2400            # the whole D row, detached, 2 threads each
```

Logs `runs/loc_<TAG>.log`, `runs/loc_<TAG>_check.log`; JSON `runs/loc_<TAG>{,_check}.json`.
