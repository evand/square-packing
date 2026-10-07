# helly-allmeet: can any positive-volume, sound rule credit the fan cliques? (2026-09-12)

**Why.**  `search/HONEST.md` (task `honest-cost`) showed that every clique row carrying dual on the
best pure `t = 4` instances is **non-Helly with an empty core** (`638/638` on `E2Pg`, `317/317` on
`SUPEPGF`; max-margin deficits `-0.0007 … -0.042`), so the verifier's sound rule `clique_of_cores`
credits them nothing off the pose set, and the pricer's rule ("credit `z_K` to any pose that
closed-meets every member") is unsound (304/317 rows admit two disjoint credited squares).  The
`62–64 %` of `Θ` that sits on these rows is what makes the packing LP read `11.80–11.87` instead of
`≥ 12`; without cliques the certifiable family does not go below 12 (`runs/PGPURE.out`,
`runs/PGCORN.out`, `notes/review-2026-09-12.md`).  `TODO.md` critical-path item 2 asks the sharp
question, **first as a measurement**:

> For each dual-carrying clique row `K`, is there a **positive-volume box `B` of poses** such that
> `K ∪ B` is pairwise closed-meeting (so that "credit `z_K` to every pose in `K ∪ B`" is a *sound*
> clique rule with positive-volume membership), and how much does such a rule recover?

If the answer is "no such `B` for the rows that carry the mass" — or "yes, but the honest cost with
the sound rule stays far above 12" — the cover route for `n = 12` is closed at `≥ 12` and the
Helly-type statement becomes the theorem to state and prove.  If the answer is "yes, and the honest
cost drops to `≈ 12`", V1 is back on the path with a certifiable object.

**Semantics** (do not change): `t = 4`; closed unit squares; closed containment and closed meeting
(a shared boundary point counts); a pose `(c_x, c_y, θ)` is admissible iff its closed square lies in
`[0,4]²`; angles as `u = tan(θ/2)` rational where exactness is needed (`search/ZEROMARGIN.md` §1,
`search/RUNG2.md` §1).  Two closed squares meet iff no one of their 8 edge normals separates them
(`leaf_ceiling.sq_meets_sq`, exact integers).  A **box** of poses is `[x0,x1] × [y0,y1] × [u0,u1]`
with `x1 > x0, y1 > y0, u1 > u0` — positive volume means all three strictly.

**Inputs** (read-only; copy what you need into your worktree's `runs/`):
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-12/` — `cl_E2Pg_{poses,rows,cliques,measure}.txt`,
`cl_E2Pg.json`, `cl_E2Pg_pgons.json` (the full 14,611-pose instance, LP `11.864926892`) and the
`cl_SUPEPGF_*` files (the certified 956-pose support, LP `11.800418089`; **start here, it is 15×
smaller**).  The dual `z_K` per clique row is obtained exactly as `HONEST.md` §6 describes
(`search/honestcost.py`: `load_state` / `solve_state` / `dump_dual`, and `cliquelever.py --dump-dual`);
`honestcost.clique_core` is the existing empty-core test and `honestcost.cmd_sound` / `pair_meets`
already find poses that meet every member of a row.  `search/leaf_ceiling.py` has the exact
square/segment/square meet tests and the `CHAIN` machinery of `zeromargin.py` is the place to
see separation conditions written as affine-in-centre, quadratic-in-`u` polynomials — the
box-vs-square "all poses in `B` meet square `S`" test you need is of that class (a sign condition
on a bounded-degree polynomial over a box; certify it exactly with an interval/Bernstein bound or
by the corner argument where the dependence is affine, and say which).

**What to compute, per dual-carrying row `K`** (all `317` on `SUPEPGF`, then all `638` on `E2Pg`):

1. `A(K) = { admissible poses S : S closed-meets every member of K }` — a float picture first: does
   it have non-empty interior?  Seed from `honestcost`'s credited poses and from `K`'s own members,
   compute for each candidate pose the *signed* meet margin against every member (the least
   separating-axis slack; `> 0` robust overlap, `= 0` grazing, `< 0` disjoint), and report the
   distribution of `min-over-members margin` over `A(K)`: if it is `≤ 0` everywhere `A(K)` has empty
   interior.  Then the largest box `B ⊂ A(K)` you can inscribe by local search (grow from the
   best-margin pose), **certified exactly** by the box-vs-square all-meet test against every member.
2. Whether `B` is itself pairwise meeting (any two poses of `B` meet) — true for a small enough
   box, but check it exactly: the worst pair is at opposite corners; state the criterion.
3. The **sound rule** `credit z_K iff S ∈ K ∪ B_K` and, with it, the honest cost of the dual
   (`honestcost.py`'s `capture`/`honest` with this rule in place of `meet`/`core`) — this is the
   headline number.  Also the fraction of clique dual carried by rows with `vol(B_K) > 0`, and by
   rows with `vol(B_K) = 0`.
4. For the rows with empty `A(K)`-interior: what is the geometry — how many members are pairwise
   grazing (`|margin| < 1e-6`), and do the grazing pairs form the "fan" (one square touched by many)
   the review describes?  A table of the top 10 rows by `z_K` with `|K|`, core deficit, number of
   grazing pairs, `vol(B_K)`, `min margin` on `B_K`.

**A caution for step 4 → conjecture.**  The Helly-type statement in `TODO.md` ("every positive-measure
pairwise-intersecting family of closed unit squares in `[0,4]²` has a common point") may be false as
literally stated: three unit squares with robustly overlapping pairwise intersections and empty
triple intersection, each thickened by a small ball in pose space, are a positive-measure pairwise
meeting family with no common point — check whether such a robust non-Helly triple exists at `t = 4`
(construct one, or show it cannot fit).  If it exists, the statement worth proving is about the **LP
mass**, not the family: e.g. that a clique row with a positive-volume pairwise-meeting extension
captures no more packing mass than its point clique, or that the packing LP's optimum over
positive-volume-creditable rows is `≥ 12`.  Say which statement the data supports, precisely, and
what a proof would need.  Do not claim a theorem; state a conjecture with the measurement that
motivates it.

**Exactness.**  Float search is fine for finding boxes; every `vol(B_K) > 0` you report must be
backed by an exact (integer / `Fraction`) certificate of "every pose in `B_K` meets every member of
`K` and every pose in `B_K` meets every other" — write the checker as a separate function with its
own docstring, and run it on every reported box.  Sound-rule honest cost is a float measurement like
`HONEST.md`'s; say so.

**Deliverables.**  `search/allmeet.py` (CLI: `boxes`, `honest`, `report`; reproducible with the
commands in the note), `search/ALLMEET.md` — verdict up front (one paragraph: does a sound
positive-volume rule recover the fan mass, yes/no, with the honest-cost number), then the tables,
the top-10 anatomy, the conjecture paragraph, what is exact vs float, reproduce.  Commit on your
worktree branch.  Do not edit `TODO.md`, `README.md`, `HONEST.md`, or other tasks' files.

**Do not**: touch pids 1109693–1109695 or anything under `.claude/worktrees/`; run anything over
10 min in the foreground (detach with `setsid nohup … &` into your worktree's `runs/`); use more
than 16 threads.  Budget: `SUPEPGF` numbers within ~2 h; `E2Pg` after; write the note with what you
have at ~5 h even if `E2Pg` is partial (say what is partial).
