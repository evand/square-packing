# s17-core-isolation: is there an s(11)-style endgame for s(17)?  (brief written 2026-09-29; ON HOLD)

**Status: ON HOLD. Do not start until Evan says so.**  There are two reasons.  (1) The Bentz / k²−3 → k²−4 work
has the CPU and the priority.  (2) Courtesy: the s(11) method is Queuingtheorydotcom's
(github.com/Queuingtheorydotcom/11SquaresOptimal, posted 2026-09-28).  We are waiting to see how quickly and how
far they extend it, and they keep priority on that.  Before starting, check whether they (or Mira-acc) have
posted s(17) work, and ask Evan.

## Background

Their s(11) proof (PROOF.md; they describe it as unverified, one implementation, no fresh public replay) has
four parts:

1. **Case split.**  A closed cover of the centre square by 16 cells, each of diameter < 1/(U−1), so no cell can
   hold two centres.  This gives 2,184 cell-pattern cases up to the half-turn.
2. **Case exclusions.**  A pose-space branch-and-bound excludes 2,180 of the cases at a rational side U slightly
   above T.  For each square it tracks the set of reachable poses (outer domain) and a region surely inside the
   square (owned hull); the owned hulls rule out poses of the other squares.
3. **Symmetry.**  An exact D4 argument reduces the 4 surviving cases to one.
4. **Capture and local isolation.**  That case is captured into a box around Trump's packing.  Inside the box a
   local theorem holds: in the fixed container [0,T]², the only feasible perturbation is zero.  The proof takes
   dual certificates on the tied contact rows and adds Taylor curvature bounds (PROOF.md §7;
   `src/evidence/research/global-math/FOCUSED-LOCAL-RECTANGLE.md`).

For s(17) (Bidwell 1997, T = 4.675530093604551…), the quantity that decides whether the endgame exists is
whether step 4 goes through.  The packing is not rigid.  Our site analysis (float, `site/www/data/p/square-17.json`,
0-based labels) finds:

* squares 5 and 12 can move alone;
* squares 4 and 10 can move only together with those (4 moving squares in all);
* the other 13 have no first-order motion.

The local theorem does not need every square pinned.  Dropping constraints only weakens the system, so drop the
4 moving squares and every constraint that mentions them.  If the 13-square core is isolated in [0,T]², any
packing near Bidwell's has its core at Bidwell's core, and therefore spans T.  This brief tests exactly that.
The global part (steps 1–3) is **out of scope**; only sizing it is in scope (§5 below).

Context: `notes/n12-gap.md` explains why this route does not apply to s(12) (the extremal set there is a
positive-dimensional zero-margin set, not an isolated point).  The s(17) lower-bound state is in `TODO.md`
(Secondary, s(17)).

## Tasks

0. **Exact packing.**
   * Derive an exact description of Bidwell's packing, starting from the site's float data and its contact list.
     Solve the tied system symbolically in tan-half-angle variables and get the minimal polynomial of T plus an
     isolating interval.  If a published exact form exists (Friedman's page, Gensane, jlevy/squares), cross-check
     against it; don't just adopt it.
   * Verify exactly in the number field (not by tolerance) that it is a packing: all 136 pairs are weakly
     separated, all squares are contained, and it spans T in both axes.
   * Record the tied rows.
1. **Core system.**
   * Take the 13 core squares plus the walls, which gives 39 perturbation coordinates.
   * For each contacting core pair, enumerate the separation features: owner of the axis × which of its two axes
     × which side.  Sort them exactly into available and unavailable at the packing, as in their 112 → 24/88 split.
   * Enumerate the branch matrices this produces, merging aliases that have identical gradients.
2. **First-order rigidity.**
   * For every branch matrix A and every signed coordinate ±e_j, find a rational λ ≥ 0 with λᵀA = ±e_jᵀ.  Exact
     equality is preferred; otherwise bound the residual ε_j rigorously.
   * If some j fails, extract the infinitesimal motion and **stop and report**.  That result is the finding: the
     endgame would need a second-order argument.
   * Watch for this failure mode: the core's rigidity may rely on rows through a moving square.  Such rows can
     still carry λ if they stay tight during the motion.  If so, retry with 4 and/or 10 put back into the core, and
     report which set is the minimal isolated core.
3. **Curvature, radii, margin.**
   * Derive your own second-derivative bounds K for the corner/axis gap functions over a coordinate box.  Use
     their formula (PROOF.md §7) only as a cross-check afterwards.
   * Choose per-coordinate radii r_j to maximise the worst ratio.  Verify, exactly:
     * every unavailable feature stays strictly negative on the box: g(0) + Σ|∂g| r + K/2 < 0;
     * for all branches, coordinates and signs: M_j < 2(r_j − ε_j R).
   * Write out the Taylor / τ-saturation soundness argument yourself, in the note.  Deliverables: the worst ratio
     (theirs was 0.68 at n = 11) and the radii.  **The radii matter most**: they set how fine the global capture
     would have to be.
4. **Negative controls.**  Perturb the packing so it is infeasible, perturb a radius beyond the verified limit,
   and drop a tied row.  The checker must reject each.
5. **Sizing the global part (estimate only).**
   * Find the minimum number of closed cells of diameter < 1/(U−1) (U just above T) that cover [0,1]², or give a
     good constructive upper bound plus a lower bound.  My guess is 26–30, against their 16.
   * Report C(cells, 17) up to symmetry.
   * List the known competitor packings with their gaps above T (site data: 4.68013 `_sym`, 4.69036 `_r2`,
     4.70711 `_r1*`).  The exclusions must kill every competitor at about 0.1% margin.
   * One paragraph: is the global case analysis plausible on 16 cores, and roughly how long would it take?

## Rules

* **Their code is not ours.**  Read their PROOF.md and FOCUSED-LOCAL-RECTANGLE.md for the mathematics, and cite
  them for the method.  Do not execute or copy their code (`~/math/_untrusted-third-party/CLAUDE.md`).  If you need
  a local copy, it goes under `_untrusted-third-party/` with a PROVENANCE entry; ask Evan first, since that tree is
  read-only.
* Checker: exact rationals / number-field arithmetic in the accept path; floats only to propose λ and r.
* Files: `search/s17_core_iso.py` (or a small package) plus the note `search/S17_CORE_ISOLATION.md`.  Commit on a
  worktree branch; no push, no merge.
* **Compute:** small; 1–2 physical cores at most.  Check what the Bentz / k²−4 runs are pinned to and stay off
  those cores.
* Report: the exact T and its polynomial; the core set; rigid yes/no (and the motion, if no); the worst ratio and
  radii; the cell-count estimate; anything in their §7 argument that looked wrong or needed repair (that is a
  finding for s(11) too).
