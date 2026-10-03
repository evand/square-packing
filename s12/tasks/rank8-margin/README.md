# rank8-margin: the margin and the minimal core of the missing unit  (2026-09-13)

**Why.**  `search/BENTZ.md` §0, §5(iii), `notes/status.md`: at `t = 4` every degree-1 family sits at 12 on
the fully pinned leaf A (four corner squares holding `{a_i, b_i}`, eight squares holding one of `c_j, d_j`
each), while integrally the leaf admits 11 — the eight non-corner singletons admit only seven together.
That rank-8 statement is the thing a proof needs.  Since tiling-minus-4 is a 12-packing with *touching*
squares, its margin under closed semantics is probably exactly zero, approached from the tiling; whether
that is so, and what the smallest sub-configuration carrying the obstruction is, decides what kind of
argument (linearised rigidity vs. interval checks) can prove it.

**Measure, in order.**
1. **Sup of the min pairwise closed gap over leaf A**: 12 labelled squares, each constrained to its pattern
   region (containment of its points, non-containment of the other `P0` points, exact half-plane tests as
   in `search/bentz.py`), objective `max min_{i<j} gap(S_i, S_j)` with `gap` the separating-axis signed
   distance (negative = overlap).  Multistart local optimisation (hundreds of starts, including the nudged
   tiling and the 162 support poses of `search/pgonly_corner_exact.txt`).  Report the best value, whether
   the maximisers converge to the tiling, and whether any tilted family gets within `1e-3` of `0` far from it.
2. **Minimal infeasible sub-configurations**: for subsets of the eight singletons (with the four corners
   always present), which are jointly realisable and which are not; list the minimal infeasible ones up to
   `D4`.  If a rank-5 or rank-6 statement on half the container already carries the obstruction, say so —
   that is the simplification we want.
3. **First-order system around the tiling**: linearise the pairwise-disjointness constraints of the twelve
   perturbed tiles (centre offsets, tilts; bounding-box cost `O(|θ|)`), test infeasibility as an LP, and
   state the neighbourhood (in centre offset and angle) on which the linearisation is valid.
4. (30 min, side deliverable) **SA level-2 on the 162-pose support**: pair marginals for the twelve
   labelled squares supported on disjoint pose pairs; the value, and which pairs carry the dual.

**Semantics.**  `t = 4`, closed squares, closed containment, packing = pairwise disjoint as closed sets
(`notes/s13-casefree.md` §1).  Floating-point is fine for 1–3 (this is a measurement of a margin, not a
certificate); say so in the write-up.

**Deliverable.**  `search/RANK8.md` with the four numbers/lists and the extremal configurations in `runs/`;
code under `search/`.  Do not edit `TODO.md`; do not touch `verify*/`, `certificates/`, `lean/`.
