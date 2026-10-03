# bentz-incidence: Bentz's case analysis as an LP-certified branch tree at `t = 4` (2026-09-12, launched 2026-09-13)

**Why.**  The cover route for `n = 12` is closed on every certifiable family measured so far
(`search/HONEST.md`, `search/ALLMEET.md`, `runs/PGCORN.out`): points + odd polygons + centre-box
regions + chord sit at **exactly `12.000` on the corner `k = 4` leaf** (`11.999999926` certified,
`search/pgonly_corner_exact.txt`), at `>= 12.17` on the pure instance, and *no* sound clique rule can
recover the fan mass (`ALLMEET.md` §0: any sound credited set is `⊆ A(K)`, hence dominated by `meet`,
hence `honest >= 18`).  The slot branch (centre boxes on the walls) is worth only `0.08` without
cliques (`search/T4LEAF.md` §0).  The extremal object on the corner leaf is always the **`4 × 4`
grid smeared**: corner mass `1,1,0.95,0.93` at `θ = 0`, wall squares `0.25–0.30` nudged `1e-3` off
the grid lines, interior grid centres `0.13–0.16`, the rest tilted (`notes/review-2026-09-12.md`).

`TODO.md` critical-path item 3 is the fallback that needs no new family: **Bentz's template**.  His
`s(13) = 4` (`notes/proof-anatomy.md` §2.3) is a 16-point set plus a 6-leaf case analysis whose
branching variable is *which points each box contains* — "a box that covers `A` but not `B`, `C`,
`D`" — and whose leaves are refuted by counting.  Every counting step is an integrality step the
fractional LP cannot make (§7.2: for 12 boxes all eight corner points can sit in four double boxes,
which is the trivial packing minus four, and that is precisely the smeared-tiling measure the LP
puts at 12).  The one object in the repo that is both **certifiable and a clique** is the *point
region* `R_p = {S : p ∈ S}` (`notes/level2-design.md` §2.2: `card_filter_clique_le_one` already in
Lean, Helly core `{p}`, positive volume) — Bentz's branching variable is exactly membership in
point regions.  So the question this task measures is:

> Branching on **incidence patterns** — for a small point set `P₀`, which points each packing square
> contains — with the counts pinned per leaf and everything else the certifiable family already has
> (points, polygons, centre-box regions, chord; **no general cliques**), do the leaves under the
> corner `k = 4` branch go below 12 at `t = 4`?  If yes, how big is the tree; if no, what does the
> extremal measure of a non-closing leaf look like?

**Semantics** (do not change): `t = 4`, closed unit squares, closed containment, a point on `∂Q`
counts (`search/ZEROMARGIN.md` §1); a packing is a family of closed unit squares that are pairwise
**disjoint as closed sets** (the dilation argument, `notes/s13-casefree.md` §1, `notes/chord-lemma.md`
B1/B2).  Consequences you may use: two packing squares cannot both contain a point; at most 3 packing
squares have centres in any wall strip `[0,4] × [0,1]` (chord lemma, Lean); at most one packing
square has its centre in a corner box `[0,1]²` (`BRANCH.md`).

**The branching object.**  Fix a finite `P₀ ⊂ [0,4]²` (candidates below).  The *pattern* of a pose `S`
is `π(S) = P₀ ∩ S` (closed membership; exact in rationals — the four half-plane tests of
`leaf_ceiling.py` / `honestcost.py`).  `R_π = {S admissible : π(S) = π}` partitions pose space.  For
`π ≠ ∅`, `R_π ⊆ R_p` for any `p ∈ π`, so **a packing has at most one square in `R_π`**, and two
patterns used by the same packing are disjoint subsets of `P₀`.  A *leaf* is a set of linear count
constraints on the pattern regions — `Σ_{S ∈ R_π} y_S = k_π` with `k_π ∈ {0, 1}` for `π ≠ ∅`,
`Σ_π k_π = 12` — and a *tree* is a finite family of leaves whose constraint sets are exhaustive: every
packing of 12 satisfies at least one leaf's constraints.  Coarse leaves are allowed and preferred
(e.g. "four squares each contain an `{A_i, B_i}` corner pair" is `Σ_{S ∈ R_{A_i} ∩ R_{B_i}} y_S = 1`
for `i = 0..3`, a union of many fine leaves); exhaustiveness has to be argued for the tree you
build, in the style of `notes/branch-semantics.md` §4, and checked by a script that enumerates
patterns up to `D4`.

**Candidate `P₀`** (try in this order; report each):
1. **Bentz's 16 points** at `m = 4`: `A(1, 0.914), B(0.914, 1), C(0.914, 2), D(1.65, 1.65)` and their
   `D4` images (`proof-anatomy.md` §2.3).  First check, with `search/zeromargin.py` or
   `verify2/zmcheck` on a unit-weight file, whether this set is a closed cover of `[0,4]²` (every
   closed unit square contains `>= 1`); if it is, `k_∅ = 0` on every leaf and Bentz's counting
   `k + u <= 4` applies verbatim.  If it is not (closed vs open box semantics can differ at a few
   poses), find the failing poses and either nudge the set or keep `k_∅` as a leaf variable.
2. The **heavy atoms of the rung-2 cover** `certificates/rung2/s13_closed_cover_4.txt` — the corner
   and grid-line atoms carrying most of its `12.956` — as `P₀`, with the cover's own weights kept as
   the cover-side objective's warm start.
3. The **8 corner points alone** (`{A_i, B_i}`), for the cheapest possible first measurement.

**What to measure**, packing side, `t = 4`, `--cq-want 0` (no general cliques), points + polygons +
centre-box regions + chord as in `runs/launch_2026-09-12_pgonly.sh` (the `PGCORN` configuration),
extended with pattern-region count rows.  Inputs: `runs/inputs-2026-09-12/cl_E2Bg_*` (the corner
instance's pose set and polygon rows) and `search/pgonly_corner_exact.txt` (its certified measure);
copy into your worktree's `runs/`.  Every leaf value is to be reported the way `PGCORN` was: the LP
value, the exact rounded-down measure with `M`, and whether the count rows are satisfied exactly.

1. **The AB-corner leaf** (the smeared tiling's home, Bentz's "adjacent" case four times over):
   `Σ_{R_{A_i} ∩ R_{B_i}} y = 1` for each corner, total 12.  This is the single most informative
   number in the task.  If it is `< 12` by a margin, the corner `k = 4` branch closes along this
   leaf and the remaining corner-pattern leaves (a corner square containing `A_i` but not `B_i`,
   neither, …) are the tree.  If it is still `12.000`, look at the measure: what patterns does the
   wall/interior mass use, and does pinning *those* (the next level) move it?
2. **The full corner-pattern level** under `k = 4`: for each corner, the pattern of the corner-box
   square with respect to `{A_i, B_i}` is one of `{A,B}, {A}, {B}, ∅` — `4⁴ / D4` leaves; measure all
   of them (most should be cheap — Bentz's Lemma 10/11 say a box with `A` and not `B` is pushed
   inward and contains other points; the LP should see the same).
3. **Whether the tree needs the wall/interior level at all**: for every leaf of (2) that does not
   close, the pattern census of its extremal measure with respect to the *whole* `P₀`, and the value
   after pinning the top pattern it uses.
4. A table like `T4LEAF.md` §0: leaf, value, `M`, interior mass, chord duals, count-row duals
   (the count rows' duals are the multipliers `λ_π` a branch certificate would carry).

**The exhaustiveness argument** must be written before the runs: the pattern of a square is a
function of its pose; patterns of a packing are pairwise disjoint; the corner box holds at most one
centre; so the corner-pattern level is a partition of the packings with `k = 4`.  Say exactly which
lemmas it uses (chord for the counts elsewhere, `card_filter_clique_le_one` for `k_π <= 1`), and what
a Lean statement of the leaf reduction would be (`packing_le_weight_regions_choice` is the model;
point regions replace centre boxes and the `<= 1` is the clique lemma, not a diameter bound).

**Do not**: build a verifier or edit `verify/`, `verify2/`, `lean/`, `certificates/` — this is a
measurement; V1-type work waits for a leaf that closes.  Do not quote any number carried by general
cliques as a bound.  Do not edit `TODO.md`, `README.md`, or other tasks' files.  No runs over 10 min
in the foreground (`setsid nohup … > runs/… 2>&1 &`); at most 16 threads; leave `.claude/worktrees/`
alone.  Exact where the family is exact (pattern membership, the final measure check via
`leaf_ceiling.py check`); say what is float.

**Deliverables.**  `search/bentz.py` (pattern enumeration up to `D4`, the leaf constraint builder,
the leaf runner; reproducible commands), `search/BENTZ.md` — verdict up front (does the AB-corner
leaf go below 12; if yes, which leaves remain and the tree size; if no, the extremal measure's
anatomy and what it says the next branching variable must be), the exhaustiveness argument, the
`P₀` closed-cover check, the table of leaf values, what is exact vs float, reproduce.  Commit on your
worktree branch.  Budget: the AB-corner leaf within ~2 h; the corner-pattern level by ~5 h; write
the note with what you have at ~6 h even if partial.

---

## Addendum 2026-09-13: priority order (read before starting)

The review of 2026-09-13 changed the emphasis of "What to measure".  Do these in this order.

**A. The fully pinned pattern leaf is the decisive number, not the corner-count pins.**  Item 1 above
(`Σ_{R_{A_i} ∩ R_{B_i}} y = 1` per corner) will almost certainly read `12.000` again — the smeared
tiling already has corner mass `≈ 1` at `θ = 0` containing both `A_i` and `B_i`, so a count pin on it
changes nothing.  The power of pattern branching is that a leaf **deletes poses**, not that it counts
them.  So the first leaf to measure is the one Bentz's counting forces once the four corners hold
`{A_i, B_i}`: with `P₀` a closed cover, `16 − 8 = 8` points remain for `8` squares, so **every
non-corner square has a singleton pattern from `{C₁..C₄, D₁..D₄}`**, one square per point.
Constraints: `Σ_{R_{{A_i,B_i}}} y = 1` (4 rows), `Σ_{R_{{C_j}}} y = 1` and `Σ_{R_{{D_j}}} y = 1`
(8 rows), and **`y_S = 0` for every pose whose pattern is anything else** (a pose containing two
non-corner points, a pose containing a corner point without its partner, a pose with pattern `∅`).
Integrally this leaf dies by `s(4) = 2` on the interior (at most one wall square per wall fits between
two pinned corners with positive gaps, then four interior squares must fit in a region of side `< 2`).
The question is whether it dies *fractionally*: report the LP value, the exact rounded-down measure
with `M`, and whether the twelve count rows hold exactly, as for `PGCORN`.  If this leaf is `< 12` by
a margin, the template has integrality power and the rest of the brief is the tree.  If it is still
`12.000`, look at the measure — which poses carry it and how they evade the point set — because that
says pattern branching on *any* `P₀` of this kind is dead and the note must say so.

**B. Then the full corner-pattern level** (item 2 above), each leaf **with the induced pattern
pinning on the remaining squares** (points freed by a corner that takes only one of `{A_i, B_i}`
make that leaf's set-packing looser: with `12` squares and `16` points, `k` doubled patterns and `u`
uncovered points satisfy `k + u ≤ 4`, Bentz's counting one unit looser; enumerate the set-packings
up to `D4`, coarsen where the LP cannot tell leaves apart, and say which coarsening you used).

**C. Two cheap diagnostics, in parallel with A, on the corner `k = 4` instance without patterns:**
1. **Axis-parallel only.**  The corner-`k = 4` value at `t = 4` with the pose set restricted to
   `|θ| ≤ 1°` (and again at exactly `θ = 0` if the tooling allows).  Nobody has measured this.  If it
   is `< 12`, the `≈ 4` units of tilted mass (`33°, 42.5°, 75°, 80°`) in the smeared measure are
   *essential*; if it is `12.000`, the obstruction is the pure axis-parallel smear and the
   axis-parallel case is where any new inequality should be tested first.
2. **Corner `k = 3` at `t = 4`, no cliques.**  `T4SCREEN.md` measured only `k = 4` and the slot
   leaves at `t = 4`; the `1110` leaf was `11.75–11.97` and rising at `t = 3.98`.  One run in the
   `PGCORN` configuration with the corner rows `1,1,1,0` (up to `D4`).  If it also sits at `12.000`,
   the corner branch has two non-closing leaves and the tree design must know that.

**D. Optional, only if A closes:** the same fully pinned leaf at `t = 3.99` in the `verify/` (sweep
verifier) setting, where margins are positive and the existing branch trailer handles regions — a
cover certificate for it there would be the first end-to-end test of pattern-leaf certification with
tools that already exist.  Do not build a verifier for `t = 4`; that is V1 and waits for the
measurement.

**Semantics reminder for `P₀`.**  Bentz's points sit *on* the grid lines (`A = (1, 0.914)`), and in
closed semantics a square with an edge on `x = 1` contains `A`; a square nudged `1e-3` off does not.
Pattern membership must be exact (rational point coordinates, the four half-plane tests as integers)
and the pose set must include the item-3 family poses of `search/family_rows.py` (nudged grid
squares), or the lattice will hide the worst poses by the factor `HONEST.md` §0 item 3 measured.
State whether Bentz's 16 points are a closed cover of `[0,4]²` (brief item 1 of "Candidate `P₀`")
before relying on `k_∅ = 0`; if a handful of poses escape, add the atoms of the rung-2 cover needed
to close them and report the enlarged `P₀`.

Report the answer to A first in `search/BENTZ.md` §0, with the raw LP/measure numbers, before
anything else in the note.
