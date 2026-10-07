# F. Level-2 branching design (wall-slot occupancy)

**When.** Only if task A says the clique relaxation alone stays `≥ 12` at `t = 4`.

**Goal.** A design, not code: from the `k = 4` leaf's dual (`runs/branch_t398hk4_dual_it16.txt`:
corners 1.0 each, 3.75 on eight wall slots, 4.0 in a tilted ring 1.0–1.6 from the walls), choose
the level-2 regions (eight wall-slot boxes by centre? by "contains a wall point"?), the multiplier
scheme (per-region `λ`, `packing_le_weight_regions` already covers families), the leaf count up to
symmetry, and the expected leaf values from the dual measure.  Include which capacity/rank columns
(`α(R) > 1`: e.g. "≤ 3 squares centred in a wall strip", from the chord lemma and disjointness) are
worth adding to the same certificate format as box cliques.

**Done when.** `notes/level2-design.md` with the tree, expected values, and the LP cost per leaf
under the backend from task C.

## Update 2026-08-30 (Phase 1 merged)

**When** is now: the clique go/no-go (`tasks/clique-continuum/`) runs in parallel; this design is the
fallback if it fails, and the two must be comparable.  New inputs since this brief was written:

* `search/CLIQUE_CEILING.md`: under clique constraints the certified 3.99 measure has corners
  **4.000** and wall slots **4.000** (integral: `(0.5, 1.5, 0°)` at 3.273 + a 7.5° tilt at 0.727 per
  slot pair) and 3.364 interior on seven tilted poses; the residual excess over 11 is entirely
  interior.  The pure measure at the same `t` has corners 3.40, walls ~4.3, interior ~4.3.
* `search/CLIQUE.md` "Plus region": 8 squares fit in `[0,t]²` minus its corner unit squares only at
  `t >= 4.0002` (heuristic), 7 fit at 3.925.  So with the corners integral the whole question is
  whether the eight non-corner squares can be shown to need `t >= 4`.
* `search/LPSPEED.md`: a cover-side leaf costs ≈ 20 min/round with `BRANCH_SOLVER=restricted`
  (~9× extrapolated at full size), so a tree of a few dozen leaves is now affordable; a tree of
  hundreds is not.
* `search/ZEROMARGIN.md`: any leaf that has to be closed *at* `t = 4` is verified by the exact
  closed-container checker; wall/corner tight poses are handled by the 2×2-box lemma, tilted tight
  families need the triangle lemma.  Design leaves so that their tight poses are of those kinds.

**Rule of thumb: explore on the packing side.**  For every leaf you propose, compute its value from
the packing side (`search/packing_dual.py` / `search/clique_ceiling.py --branch` style LP with the
leaf's occupancy constraints; ~1–3k poses, seconds to minutes) before costing it on the cover side.
Deliver: the tree (regions, multipliers, leaf count up to symmetry), each leaf's packing-side value at
`t = 3.98` and `t = 4.0` (heuristic, say in which direction it errs), which leaves already close,
which need cliques, which need `t = 4` zero-margin verification, and the cover-side cost estimate.
`<= 2` cores.  `runs/` is gitignored: read `/home/evand/math/square-packing/s12/runs/` by absolute
path.  Commit in your worktree; do not touch `TODO.md` or other tasks' files.
