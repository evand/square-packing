# bentz-incidence: Bentz's case analysis as incidence-pattern branching, measured at `t = 4` (2026-09-13)

Task: `tasks/bentz-incidence/README.md` (+ its "Addendum 2026-09-13: priority order").
Code: `search/bentz.py` (new) and a `--pose-filter` hook in `search/cliquelever.py` (~25 lines;
inert without the flag).  Nothing in `verify/`, `verify2/`, `lean/`, `certificates/` is touched:
this is a measurement.  Runs live in this worktree's `runs/` (gitignored), so every number is
quoted here.

Read against `notes/proof-anatomy.md` §2.3 and §7.2, `notes/review-2026-09-12.md`,
`search/HONEST.md` §0, `search/ALLMEET.md` §0/§5, `search/T4LEAF.md` §0, `search/T4SCREEN.md` §0–§1,
`runs/launch_2026-09-12_pgonly.sh` (the `PGCORN` configuration every run below extends).

---

## 0. Verdict, up front

**A — the fully pinned pattern leaf does not die, fractionally or otherwise: it is exactly where
the `PGCORN` optimum already lives.**  The certified `PGCORN` measure
(`search/pgonly_corner_exact.txt`, mass `5999999963/500000000 = 11.999999926`) satisfies **every**
constraint of the fully pinned leaf: all 162 of its poses have one of the twelve admitted patterns,
and the twelve count rows `mu(R_pi) = 1` hold with deficits between `1e-9` and `1.3e-8` — the
round-down of `finalize`, and nothing else.  So

> **leaf `(AB)^4`, fully pinned, `t = 4`, no cliques:  `11.999999926 <= value <= 12`, exactly.**

The lower bound is an exact certified measure (`M = 1` exactly over all 14,805 arrangement
vertices, all four wall strips exactly `3`, all four corner regions `1.000000`, all 1,919 polygon
rows satisfied).  The upper bound is not a measurement at all: the twelve admitted patterns are
pairwise distinct, every admitted pose has exactly one of them, and `mu(R_pi) <= mu(R_p) <= 1` for
any `p in pi` by the clique lemma, so the leaf's objective is `<= 12` on **any** pose set.  There is
no pose-set-richness caveat on either side.  The leaf is closed as a question: **pattern branching
on Bentz's 16 points has zero power at `t = 4`.**

Per-pattern masses of that measure, exactly (`python3 search/bentz.py patterns runs/pgonly_corner_exact.txt`):

| pattern | mass | `1 - mass` | poses carrying it |
|---|---|---|---|
| `{a0,b0}` | `1` | `0` | 1 |
| `{a1,b1}` | `999999999/1000000000` | `1e-9` | 4 |
| `{a2,b2}` | `1` | `0` | 2 |
| `{a3,b3}` | `1` | `0` | 1 |
| `{c0}` | `199999999/200000000` | `5e-9` | 14 |
| `{c1}` | `124999999/125000000` | `8e-9` | 14 |
| `{c2}` | `199999999/200000000` | `5e-9` | 10 |
| `{c3}` | `249999999/250000000` | `4e-9` | 10 |
| `{d0}` | `999999987/1000000000` | `1.3e-8` | 27 |
| `{d1}` | `249999997/250000000` | `1.2e-8` | 24 |
| `{d2}` | `249999997/250000000` | `1.2e-8` | 26 |
| `{d3}` | `499999993/500000000` | `1.4e-8` | 29 |
| **total** | `5999999963/500000000` | `7.4e-8` | 162 |

The addendum's expectation was that leaf A "deletes poses".  It deletes `1,746` of the `14,835`
poses of the corner instance `cl_E2Bg_poses.txt` — **12 %** — and **none** of the 162 poses that
carry the optimum.  That is the whole finding: the leaf whose integral version dies instantly by
`s(4) = 2` is a leaf whose fractional version the LP was already sitting inside.

**C1 — axis-parallel only: `9.000000` exactly, converged, integral.**  The corner `k = 4` leaf at
`t = 4` with the pose set restricted to `theta = 0 mod 90 deg` **exactly** converges in 71 s to
`LP = 9.000000`, support 9 poses each of mass exactly `1`, `M = 1` exactly, corner regions
`[1,1,1,1]`, all 1,872 polygon rows satisfied.  That is the `3 x 3` staggered packing, and it is
also what the continuum gives (uniform density `1` over centres in `[1/2, 7/2]^2` has coverage
exactly `1` and mass `9`).  Relaxing to `|theta mod 90 deg| <= 1.0000034 deg` lifts it to at least
`10.75` (still running, §4).  **So the tilt is essential, and much more essential than the
`~4` units of visibly-tilted mass in the smeared optimum suggested: `3` of the `12` cannot be
reached without tilt at all, and `1.25` of the remainder needs more than a degree of it.**  This
answers the addendum's alternative in the negative: the obstruction is *not* the pure axis-parallel
smear, and an inequality tested on the axis-parallel case first would be testing a case that is
already `3` below `12`.

**C2 — corner `k = 3`: still rising, `>= 11.883193` after two iterations (running).**  See §4 for
the live value and §7 for the detached-run record.

**What this says the next branching variable must be** (§6): not incidence with any finite point
set.  The leaf's mass is carried by `10` to `29` poses per pattern that differ by `1e-3` in centre
and by tens of degrees in angle; a count row over a pattern region cannot see the difference
between one square and a cloud of near-squares, and `s(4) = 2` — the statement that kills the leaf
integrally — is a statement about *four pairwise-disjoint* squares in a region of side `< 2`, i.e.
about a clique-like object of four poses, not about point membership.

---

## 1. What the branching object is, and the semantics

Fix a finite `P0 subset [0,4]^2`.  The **pattern** of a pose `S` is `pi(S) = P0 ∩ S` with **closed**
containment (a point on `dS` counts, `search/ZEROMARGIN.md` §1).  `R_pi = {S admissible : pi(S) = pi}`.

Semantics, unchanged from the brief: `t = 4`, closed unit squares, closed containment, closed
container; a packing is a family of closed unit squares **pairwise disjoint as closed sets** (the
dilation argument of `notes/s13-casefree.md` §1).  Consequences used below:

* **(P1)** two packing squares never contain a common point;
* **(P2)** hence for `pi != empty` and `p in pi`, `R_pi ⊆ R_p = {S : p in S}` and `R_p` is a clique,
  so **a packing has at most one square in `R_pi`** — in Lean this is
  `card_filter_clique_le_one` (`lean/Sqpack/Basic.lean:249`) applied to the predicate
  `K c θ := p ∈ sq c θ 1`, whose clique hypothesis `hK` is discharged by `p` itself;
* **(P3)** the patterns used by the squares of one packing are pairwise disjoint subsets of `P0`;
* **(P4)** at most one packing square has its centre in a corner box `[0,1]^2` (`BRANCH.md`), which
  is what makes the corner `k = 4` branch a partition in the first place.

`P0` is Bentz's 16 points at `m = 4` (`notes/proof-anatomy.md` §2.3): `A(1, 457/500)`,
`C(457/500, 2)`, `D(33/20, 33/20)` and their `D4` images — exactly, as rationals, with Bentz's own
decimal `0.914 = 457/500` (not `sqrt2 - 1/2 = 0.9142…`; he rounds down so the lemma hypotheses
hold, and we keep his number).  The orbit structure is `8 + 4 + 4`:

* `a_i, b_i` (`i = 0..3`) — the eight corner points, one pair per corner box `C_i`; `a_i` sits on a
  vertical grid line `x in {1,3}`, `b_i` on a horizontal one.  They form **one** `D4` orbit of
  size 8 (the diagonal reflections swap `a <-> b`).
* `c_0..c_3` — the mid-wall points, the orbit of `C(0.914, 2)`.
* `d_0..d_3` — the interior points, the orbit of `D(1.65, 1.65)`.

`python3 search/bentz.py points` prints all sixteen with their exact coordinates and the induced
`D4` permutations.

---

## 2. The exhaustiveness argument (written before the runs)

**Claim.**  Under the corner `k = 4` branch, the 256 tuples `(sigma_0, …, sigma_3)` with
`sigma_i in {{a_i,b_i}, {a_i}, {b_i}, empty}` partition the packings of 12, and up to `D4` there
are **43 classes**.

*Proof.*  Four steps.

1. **The pattern is a function of the pose.**  `pi(S) = P0 ∩ S` is determined by `S`, so
   `{R_pi}_{pi ⊆ P0}` is a partition of pose space.  (Exhaustiveness at this level is free; it is
   the *refinement* below that needs an argument.)
2. **Corner lemma: a square centred in `C_i` sees only its own pair.**  A closed unit square is
   contained in the closed disc of radius `sqrt(2)/2` about its centre.  Every point of `P0` other
   than `a_i, b_i` is at squared distance `>= 169/200 = 0.845 > 1/2` from the box `C_i` (the
   closest is `d_0 = (1.65,1.65)` from `C_0 = [0,1]^2`; checked exactly by
   `bentz.corner_reach_check`, printed by `python3 search/bentz.py leaves`).  Hence
   `pi(T_i) ⊆ {a_i, b_i}` for the corner square `T_i`, and `sigma_i := pi(T_i)` really does range
   over the four listed sets.
3. **The four choices are exhaustive and mutually exclusive** for each `i` by definition of
   `sigma_i`, and by **(P4)** the corner square `T_i` is unique in the `k = 4` branch.  So the map
   (packing) `-> (sigma_0, …, sigma_3)` is total and single-valued, and the 256 fibres partition
   the `k = 4` packings.
4. **`D4` reduction.**  `P0` is `D4`-invariant (checked: `bentz.d4_perm` is a permutation of the
   sixteen indices for each of the eight group elements), the container is `D4`-invariant, and the
   corner-count vector `(1,1,1,1)` is `D4`-invariant.  So `D4` acts on the 256 tuples — permuting
   corners and, on the diagonal reflections, swapping `A <-> B` within each corner — and the value
   of a leaf is constant on an orbit.  The orbit count is **43** (orbit sizes summing to 256;
   `python3 search/bentz.py leaves` prints the representatives, orbit sizes and Bentz's counting
   bound per class).  ∎

**Bentz's counting at 12 boxes, and why `(AB)^4` is already fully pinned.**  With 12 squares and
16 points, patterns pairwise disjoint (**P3**) and non-empty (**§3**), write
`K = Σ_S (|pi(S)| - 1)` and `u = #{points in no square}`.  Then `Σ_S |pi(S)| = 16 - u` and
`Σ_S |pi(S)| = 12 + K`, so

> **`K + u = 4`**

— Bentz's `k + u <= 3` for 13 boxes (`notes/proof-anatomy.md` §2.3, §7.2), one unit looser, exactly
as the addendum says.  On the class `(AB)^4` the four corner squares already contribute `K >= 4`,
so `K = 4` and `u = 0`: **every other square has a singleton pattern, and the eight remaining
points `c_0..c_3, d_0..d_3` are used once each by the eight non-corner squares.**  So on this class
the "corner-pattern level" and the "fully pinned leaf" of addendum A are the *same leaf* — the
induced pinning is forced, not an extra choice.  The other 42 classes have `K >= #{i : sigma_i = AB}`
and `u <= 4 - K`, and their set-packings of the remaining points are a genuine further branch; §5
says why measuring them is beside the point once `(AB)^4` reads 12.

**What a Lean statement of the leaf reduction would be.**
`packing_le_weight_regions_choice` (`lean/Sqpack/Basic.lean:207`) is the model and it already has
the right shape: its regions are arbitrary predicates `R : Fin m → (ℝ×ℝ) → ℝ → Prop` on poses, with
neither disjointness nor exhaustiveness required, and its conclusion `n + Σ_i lam (σ i) ≤ W` is
relative to an arbitrary assignment `σ` of each square to a region containing its pose.  A pattern
region is the predicate

```
R_pi c θ  :=  ∀ p ∈ P0, (p ∈ sq c θ 1 ↔ p ∈ pi)
```

and the leaf's `<= 1` is **not** a diameter bound but `card_filter_clique_le_one` at the point
`p ∈ pi` (the clique hypothesis `hK` is witnessed by `p`), where `packing_le_weight_regions_cliques`
uses a box clique.  So the two Lean ingredients exist; what is Python and not Lean is exactly what
`notes/branch-semantics.md` §5 already flags — the *enumeration* (that the 43 classes are all of
them, i.e. steps 2 and 4 above), which here is `bentz.corner_reach_check` plus an orbit count.

---

## 3. Is `P0` a closed cover of `[0,4]^2`?

The brief's item 1 and the addendum both ask this before `k_empty = 0` may be used.  **The answer
is: yes on every test I can run, but `zeromargin.py` cannot certify it, and the reason is
instructive.**

* `python3 search/zeromargin.py cert runs/bentz16_cover.txt --tri --disj --depth 18` (unit weights,
  the exact adaptive box sweep; the point set is symmetric so the reduced domain applies):
  **`NOT VERIFIED`** — `31,110` boxes, `ADM 10963 / CHAIN 489 / TRI 12 / EMPTY 3775 /
  UNCERTIFIED 3516`, 15 s (`runs/bentz16_cover.log`).
* But none of the uncertified boxes contains a counterexample.  All `94,933` oracle poses of the
  uncertified boxes have a **non-empty** pattern (exact test).  A fine float scan of the
  uncertified band (`cx in [w/2, 0.80] x cy in [1.40, 1.60] x theta in [40°, 44°]`, pitch `2e-4` /
  `0.005°`) bottoms out at containment margin **`+0.000222503`** at `(0.706996, 1.494800, 43.985°)`.
  A global scan (`theta in [0,45°]` step `0.1°`, centre pitch `0.004`, the `D4` fundamental domain)
  finds **no** admissible pose with a negative margin; the global minimum is exactly `0`, attained
  at the sixteen squares of the `4 x 4` tiling, where the points lie *on* an edge.
* So `P0` is a closed cover but a **tight** one — margin `0` at the tiling and `2.2e-4` in the
  `44°` band — which is precisely why the box checker's `CORE`/`ADM`/`P1` primitives fail on it:
  they need a positive core over a whole angle bin.  Bentz's own proof of unavoidability at these
  poses is Lemma 5 (Stromquist's `f(a)`, a one-parameter trigonometric inequality) and Corollary 3,
  which `notes/proof-anatomy.md` §7.4 already lists as *not* box-certifiable.

**Status used below:** `k_empty = 0` is a hypothesis with strong numerical evidence, not a
certified fact.  It is not load-bearing for §0: leaf A's lower bound is an explicit measure, and
its upper bound `<= 12` needs only that the twelve admitted patterns are distinct.

---

## 4. The table of leaf values

All at `t = 4`, closed semantics, `r = 1`, four corner boxes pinned at `1` unless said otherwise,
`--cq-want 0` (no general cliques, no number below is carried by one), points + odd polygons +
region equalities + chord — the `PGCORN` configuration of `runs/launch_2026-09-12_pgonly.sh`.

| run | leaf | LP | exact measure | `M` | count rows | interior | note |
|---|---|---|---|---|---|---|---|
| `PGCORN` (2026-09-12, quoted) | corner `k = 4`, no patterns | `12.000000` | `5999999963/500000000 = 11.999999926` | `1` exact | — | `4.000000` | the baseline |
| **A** | **`(AB)^4` fully pinned** | **`[11.999999926, 12]` exactly** | the `PGCORN` measure, unchanged | `1` exact | **all 12 hold to `1.3e-8`** | `4.000000` | §0; the leaf deletes `1746/14835` poses and `0` of the optimum's |
| `BZC1B` | corner `k = 4`, `theta = 0 mod 90°` exactly | **`9.000000`** converged | `9` exactly (9 poses, mass `1` each) | `1` exact | — | `1.000000` | integral; 71 s; slots `[0,1,0,1,1,0,0,1]` |
| `BZC1A` | corner `k = 4`, `|theta mod 90°| <= 1.0000034°` | `>= 10.750000` | — | `1` exact at s0.0 | — | `2.750000` | running, §7 |
| `BZC2` | corner `k = 3` (`--corners 1110`) | `>= 11.883193` | — | `1.372` at s0.1 | — | `3.171161` | running, §7 |
| `BZA` | leaf A re-measured with the 16 points as coverage rows | — | — | — | — | — | running, §7 (confirmatory) |
| `BZAC` | leaf A with the 12 count rows as equalities | — | — | — | — | — | running, §7 (confirmatory) |

Notes on the table.

* **A is not a run.**  It is the observation that the certified `PGCORN` measure already satisfies
  the leaf, plus the two-line upper bound of §0.  `BZA`/`BZAC` re-derive it from the LP; they are
  confirmation, not the evidence.
* The `PGCORN` measure's exact re-check in this worktree (`runs/bentz_A_check.log`,
  `python3 search/leaf_ceiling.py check … --corners 1111 --chord --anchor none`) reproduces
  `MASS = 5999999963/500000000`, `M = max cov = 1` over `14,805` arrangement vertices, wall strips
  `S0..S3 = 3.000000` each `<= 3 OK`, and corner masses `C0..C3 = 1.000000` with
  `I = 4.000000`.  The region **target** check prints `FAIL` for one reason only, the one the
  2026-09-12 review already recorded: a single pose of mass `0.024234` has its centre on a region
  boundary (`c_x = 1`) and no assignment of that one pose makes all four corner equalities hold
  *exactly* after the round-down.  Coverage, chord and the polygon rows are unconditional.
* **Chord duals** are `0` in every run, as in `T4LEAF.md` §0 item 3: the corner equalities already
  pin each wall strip at `3`.
* **Count-row duals.**  With `--pose-filter` carrying `"counts"`, the twelve equality rows enter
  after the region rows and their duals are read back as `lever.count_dual` and printed on the
  `STAGE` line.  They are the multipliers `lam_pi` a branch certificate would carry.  In that mode
  the lattice pricer's reduced cost does **not** charge them (it charges region, clique and polygon
  duals only), so "no improving column" is not a convergence claim there; that is why `BZA` (no
  equalities, the caps implied by the coverage rows at `P0`) is the primary run and `BZAC` the
  cross-check.

---

## 5. The anatomy of the non-closing leaf, and why level B does not need measuring

Three things about the optimum, all exact.

**(i) The leaf forbids the `4 x 4` tiling and forces the nudge — which the optimum already used.**
Patterns of the sixteen exact tiling squares (`c = (i+1/2, j+1/2)`, `theta = 0`):

```
 {a2,b2}  {a2,c2}  {a3,c2}  {a3,b3}
 {b2,c0}   {d1}     {d3}    {b3,c3}
 {b0,c0}   {d0}     {d2}    {b1,c3}
 {a0,b0}  {a0,c1}  {a1,c1}  {a1,b1}
```

The four corner squares and the four interior squares are admitted; **all eight wall squares are
forbidden**, because each contains a corner point *and* a mid-wall point — a doubled pattern that
`K + u = 4` has no budget for once the corners are doubled.  Nudge a wall square by `1e-3`
*towards the wall's centre* and its pattern collapses to the singleton `{c_j}`, which is admitted;
nudge it *away* and it collapses to `{a_i}`, a corner point without its partner, which is not.  The
`PGCORN` optimum's wall mass is exactly the admitted nudge (`0.25–0.30` per wall square at `1e-3`
off the grid lines — `notes/review-2026-09-12.md`), i.e. the item-3 family of
`search/family_rows.py` that `HONEST.md` §0 item 3 found the pricing lattice hiding.  **The leaf
does not delete the optimum; it deletes the optimum's un-nudged idealisation.**

**(ii) Each unit of pattern mass is a cloud, not a square.**  `{d0}` is carried by 27 poses,
`{d3}` by 29, `{c0}` by 14 — centres within about `1e-3` of each other in the axis-parallel part
and angles spread over `33°, 39°, 42.5°, 75°, 80°` in the tilted part.  A count row
`mu(R_pi) = 1` is satisfied by any such cloud.  The integral statement that kills the leaf —
between four pinned corner squares with positive gaps at most one wall square per wall fits, and
then four interior squares would have to fit in a region of side `< 2`, i.e. `s(4) = 2` — is a
statement about *four pairwise-disjoint squares*, and the only object in the repo that expresses
"these four cannot coexist" is a clique row, which `ALLMEET.md` §0 has already shown is not soundly
creditable where the mass is.

**(iii) So level B is decided by one leaf.**  A branch tree closes only if **every** leaf closes.
`(AB)^4` is one of the 43 classes and it is at `12.000` with an exactly certified measure, so the
corner-pattern level does not close, and neither does any refinement of it that keeps `(AB)^4`
whole.  The remaining 42 classes are therefore measured (if at all) for tree-size bookkeeping, not
for the verdict; §2 records their Bentz budget `K >= #{sigma_i = AB}`, `u <= 4 - K`, which is the
information a tree design would need.  The coarsening I would use if it were worth doing is the one
§2 proves is free: `(AB)^4` needs no sub-branch at all, because `K + u = 4` forces the rest.

---

## 6. What the measurement says the next branching variable must be

Three constraints on it, from the numbers above:

1. **Not incidence with a finite point set.**  The one pattern leaf whose integral version is
   hopeless is at `12.000` fractionally, and it is not close: the LP loses nothing at all on it.
   Enlarging `P0` cannot help in the direction that matters — more points make more patterns, and
   the optimum's clouds just acquire finer patterns; what a leaf would have to forbid is a *cloud*,
   and a pattern region is a union of clouds by construction.
2. **The axis-parallel case is not where to test.**  `C1` puts the axis-parallel corner `k = 4`
   leaf at exactly `9`, `3` below `12`, and the `1°` band at `>= 10.75`.  Every one of the missing
   `3` units is bought by tilt, and `1.25` of them by tilt of more than a degree.  Any new
   inequality has to be a statement about tilted squares; a purely axis-parallel argument is
   answering a question that is already answered.
3. **It has to be a statement about a few squares at once, with positive volume in pose space.**
   That is `ALLMEET.md` §5's conjecture from the other side: the `s(4) = 2` step is a
   4-square statement, the clique family is the repo's 4-square object, and the clique family is
   uncreditable.  The gap between "`s(4) = 2` kills leaf A integrally" and "the LP sits at 12 on
   it" is exactly one unit of `notes/proof-anatomy.md` §7.2's deficit, and it is now localised to a
   single, completely explicit object: the four singleton patterns `{d_0}, {d_1}, {d_2}, {d_3}`
   carrying mass `1` each on `24–29` mutually overlapping poses inside `[1,3]^2`.  A certifiable
   rule that says "these four cannot each carry `1`" is the whole remaining content of the Bentz
   template for `n = 12`.

---

## 7. Runs still going, and how to read them

Launch scripts: `runs/launch_bentz_C1.sh` (C1), `runs/launch_bentz_C2.sh` (C2),
`runs/launch_bentz_A.sh` (A), `runs/launch_bentz_cover.sh` (§3), `runs/launch_bentz_check.sh` (the
exact re-check).  Each `cliquelever` run writes `runs/<TAG>.out`, `runs/cl_<TAG>.json` and
checkpoints `runs/cl_<TAG>_{measure,poses,rows,pgons}.*`.

_(live status recorded in §7.1 below at the end of the session)_

---

## 8. What is exact and what is float

**Exact** (integers or `Fraction`, no float decides anything):

* **Pattern membership.**  `bentz.pattern_of` is `leaf_ceiling.sq_contains` — four integer
  half-plane tests in the square's own frame, on a point `(X, Y, D)` with `D = 500` and a pose
  `(p, q, cx, cy)` with `theta = 2 arctan(p/q)`.  No tolerance, no band.  The `--pose-filter`
  admission test, the count-row index of a pose, and the census of §0 all use it.
* **The angle filter of C1.**  `min(|a|, |b|) * den <= num * r` with
  `(a, b, r) = (q^2 - p^2, 2pq, q^2 + p^2)`; `num/den = 17453/1000000` for the `1°` band
  (`|theta mod 90°| <= arcsin(0.017453) = 1.0000034°`) and `num = 0` for exactly axis-parallel.
* **The corner lemma** of §2 step 2: squared distances in `Fraction`, `169/200 > 1/2`.
* **The `D4` orbit computation** and the 43-class count: exact rational coordinates, so the orbit
  of a point is decided by equality of `Fraction`s.
* **The measure check.**  `leaf_ceiling.py check` — mass, the maximum coverage `M` over *every*
  vertex of the arrangement of the support squares, region masses, wall-strip masses: all exact.
  `cliquelever`'s `FINAL EXACT` line likewise (it rounds the LP solution down to `1/10^9` and
  re-certifies).
* **The `P0` cover file** `runs/bentz16_cover.txt` and `zeromargin.py`'s sweep over it.

**Float**:

* **Every LP value** (`LP=…` in the logs) is HiGHS in double precision.  A value is turned into a
  claim only through the rounded-down exact measure and its exact re-check.
* **The pricing reduced costs and the polygon separator** are float and heuristic, as in
  `RANKDIAG.md` / `HONEST.md`; "converged" means no violated row and no improving column *on the
  lattice the pricer looks at*, which is evidence, not a theorem.  On top of that, with
  `--pose-filter` the pricer ranks candidates *before* the filter rejects them, so the printed
  pricing `gap` is the unfiltered maximum and overstates what the filtered run could still gain;
  the honest statement for `BZC1B` is "two successive full lattice passes contributed 2 and then 0
  admitted poses, and the LP over the admitted columns is at its exact optimum".
* **The §3 cover scans** (global and local) are float; they are the reason for saying "no
  counterexample found", not "certified".  The only exact statement in §3 is the negative one:
  `zeromargin.py` leaves 3,516 boxes uncertified.

**Not quoted as a bound anywhere**: any number carried by a general clique row.  Every run above is
`--cq-want 0` and `--ktol 1000`, so the `kmax` field in the logs is decorative and the `certified`
flag's clique component is meaningless here, exactly as `runs/launch_2026-09-12_pgonly.sh` says.

---

## 9. Reproduce

```sh
# the point set, the D4 action, the 43 corner-pattern classes and the corner lemma
python3 search/bentz.py points
python3 search/bentz.py leaves

# the verdict of section 0: the pattern census of the certified PGCORN measure
cp /path/to/s12/search/pgonly_corner_exact.txt runs/
python3 search/bentz.py patterns runs/pgonly_corner_exact.txt

# its exact re-check (M, regions, chord)
sh runs/launch_bentz_check.sh              # -> runs/bentz_A_check.log

# is P0 a closed cover?
sh runs/launch_bentz_cover.sh              # -> runs/bentz16_cover.log

# the leaf specs, the pose seeds and the runs
python3 search/bentz.py spec runs/specs/leafA.json        --leaf AB,AB,AB,AB --full
python3 search/bentz.py spec runs/specs/leafA_counts.json --leaf AB,AB,AB,AB --full --counts
python3 search/bentz.py seed runs/seed_leafA.txt  --pitch 0.02 --spec runs/specs/leafA.json
python3 search/bentz.py seed runs/seed_axis1.txt  --pitch 0.02 --spec runs/specs/axis1deg.json
python3 search/bentz.py seed runs/seed_axis0.txt  --pitch 0.01 --spec runs/specs/axis0.json
sh runs/launch_bentz_C1.sh                 # BZC1A, BZC1B
sh runs/launch_bentz_C2.sh                 # BZC2
sh runs/launch_bentz_A.sh                  # BZA, BZAC
```

Inputs used, read-only from the main checkout and copied into `runs/`:
`runs/inputs-2026-09-12/cl_E2Bg_{poses,rows,pgons}.*` and `search/pgonly_corner_exact.txt`.
