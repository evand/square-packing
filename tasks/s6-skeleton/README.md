# s6-skeleton: re-prove `s(6) = 3` with a separation / transversal skeleton, as the rehearsal for `s(12) = 4`  (2026-09-20)

**Why.**  `search/CENSUS.md`: at `n = 12`, `t = 4` the point-pattern tree is blind to the obstruction
(`>= 2,779` zero-margin leaves); what holds the plateau shut is a wall-to-wall row of four (99 %) or a
corner-contact "pinwheel" core (axis-parallel only).  The proposed next step is a proof skeleton that branches
on separation structure and kills zero-margin regions by *exact* arguments.  Before building that in 36
dimensions, build it where the answer is known: **six closed unit squares, pairwise disjoint as closed sets, do
not fit in `[0,3]^2`** (`s(6) = 3`, Kearney–Shiu 2002, `notes/proof-anatomy.md` §5.2).  It is the `n^2 - n` case
at `n = 3`, 18 dimensions, margin zero (tiling minus 3), has the same two core types (row of three; tiling minus
a permutation = pinwheels, every row and column with slack 1), and the same tilted competitor one square down
(`s(5) = 2 + 1/sqrt 2 < 3`).  The `n = 11` rank-10 rehearsal is margin-positive and does not exercise any of this.
If the skeleton cannot close `s(6)` mechanically it will not close `s(12)`; if it can, its leaf counts are the
scaling estimate.

**Sanity filter (binding).**  `s(n^2 - n) = n` is **false** for large `n`: Cleemann packs `272 = 17^2 - 17` in
side `< 17` (Friedman DS7).  So any kill that works uniformly in `n` is wrong.  For every exact lemma you use,
say where the chain length (`3` here, `4` at `n = 12`) enters quantitatively, and what it would say for a chain
of 17.  The axis-parallel Dilworth argument (`CENSUS.md` §3) is uniform in `n` and is fine only because it
assumes no tilted-only pair; the content is whatever handles those.

**Build, in order.**
1. **Instrument.**  Closed gap for six squares in `[0,3]^2` (reuse `search/rank8.py` `pair_gap`,
   `search/chains.py analyse`).  Reproduce the picture first: margin exactly `0`; the zero set; the core type at
   each stopping point (expect rows of three and pinwheels); is any core ever tilted?
2. **Skeleton.**  Branch-and-bound over configuration boxes (centres, angles mod 90°, `D4` and label symmetry
   quotiented), with three kinds of kill, counted separately:
   * **(N) negative margin**: some pair overlaps, or a square leaves the container, on the whole box (interval
     bound).  Ordinary; cannot terminate near the zero set.
   * **(K1) transversal chord lemma** (candidate answer to `CENSUS.md` §3's "what replaces width `>= 1`"): a line
     crossing a unit square through two *opposite* sides has chord `1 / cos(angle) >= 1` at any tilt; three
     pairwise-disjoint closed chords of length `>= 1` do not fit on a closed segment of length `3`.  For a square
     tilted `theta in [0, 45°]` the horizontal lines that cross opposite sides form the band
     `|y - y_c| <= (cos theta - sin theta) / 2`.  "The bands of three squares share a height" is an **open**
     condition (checkable on a whole box with slack) giving an **exact** kill with no axis-parallel hypothesis.
     Same for vertical lines.  This should eat neighbourhoods of the row-of-three plateau.
   * **(K2) whatever kills the pinwheels.**  Unknown; this is the point of the task.  At the axis-parallel
     pinwheel, pushing a square to open the corner contact puts it on a line already carrying two full chords
     (K1 fires); the escape is a *tilted* square poking a corner across that line, a corner-cut chord of length
     `d (tan theta + cot theta)` at depth `d`.  Find the exact statement (chord sums on every axis-parallel line
     `< 3`?  Kearney–Shiu's Lemma 2/3 inequalities?  something else) and whether it is an open condition.
3. **Report what is left.**  Boxes killed by N / K1 / K2, the residual set described geometrically (not just
   counted), and whether the residue shrinks under refinement or is a genuine third mechanism.  "Did not close,
   and here is the residue" is a fine result; a vague "mostly closes" is not.
4. **Transfer.**  One section: which of this carries to `n = 12` verbatim (K1 with `4` for `3`), what changes
   (eleven-square sub-families are realisable with margin `3.8e-6`, `RANK8.md` §2), and an estimate of the box
   count at 36 dimensions with the `P0` pattern leaves as the outer layer.

**Semantics.**  Closed squares, closed container, packing = pairwise disjoint as closed sets
(`notes/s13-casefree.md` §1).  Floating point is fine for the instrument and for exploring; every kill must be
*stated* as an exact lemma with hypotheses, and if you claim a box count is rigorous, say which arithmetic made
it so.  Do not claim a proof of `s(6) = 3` unless every box is closed by a stated lemma.

**Housekeeping.**  Deliverable `search/S6_SKELETON.md` (verdict first, every number quoted in the file since
`runs/` is gitignored) + code `search/s6skel.py`.  Do not edit `TODO.md`; do not touch `verify*/`,
`certificates/`, `lean/`, or existing `search/*.py` (import them).  Anything over 10 minutes: detached
(`setsid nohup`), 16 threads is the machine's width, never `pkill`, and leave no `until … sleep` / `tail -f`
watchers behind.  Do not commit.
