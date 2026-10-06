# n = 12: lemmas from the group-chat 3.99 push (2026-10-05)

Source: `~/math/square-packing/n12_progress_2026-10-04.md` (group-chat checkpoint, not in this repo; its evidence
tarball `n12_progress_evidence_2026-10-04.tar.gz` is not on this machine).  Credit them before anything here goes public.

Their target is `L = 3.9901`, just past our pure-cover ceiling (`search/DUAL_EXACT.md`: no cover < 12 at `t >= 3.99`),
by a pose-state case tree.  5 residual roots open (the permutation-hole families).  Not a route to `t = 4` (it needs
slack); the lemmas are.

## Lemmas worth having, and whether they hold at t = 4

* **Wall strip, h = 603553/500000 (just under (1+√2)/2):** at most 3 centres within h of a wall, every side `L < 4`.
  Ours (`lean/Sqpack/Chord.lean`) has h = 1.  Proof: certificate, 8,155,478 leaves.  Transfers to closed `t = 4` by the
  dilation of `notes/s13-casefree.md` §1 (centres move toward the walls).
* **Corner–edge:** `U + V >= (391/100) c - 191/100` (their `prove12_local/corner_edge_unconditional.md`).  No container
  side in it, so it holds at 4.  Analytic + one rational Bernstein quartic: Lean-friendly.  4,000 adversarial local
  optimisations here (scratch, 10-05): no violation, slack -> 0 only as the tilt -> 0.
* Chord band `B(θ)` = our K1; window capacity = our K2/total band.  Max-tilt branching and per-family point
  certificates are tree machinery: need slack, so `t < 4` only.
* Sanity filter: six equally tilted squares realise their "L-shaped" local pattern at side 3.949751.

## The counting sketch (10-05, unchecked beyond this)

Regions: corner `K_i` (centre within h of two walls), edge `E_j` (one wall), centre `Z` (none).  Summing the four
strip rows: `ΣE + 2ΣK <= 12`, so

    N = ΣK + ΣE + Z  <=  12 + Z − ΣK.

`K_i <= 1` (centre domain ⊂ [1/2, h]², diagonal ~1).  The tiling minus a permutation hole set with the central 2×2
intact has `ΣK = Z = 4`, `N = 12` — touching only.  So `s(12) = 4` along this line is **`Z <= ΣK − 1` for closed
packings**: the equality case (4 corners, each strip saturated, central 2×2 full) is a rank-8-type rigidity statement,
and their cyclic corner→edge→…→corner theorem is aimed at exactly it (gap: central contact ownership).  **`Z <= 4` is false as a local statement** at any h <= h*: five
squares (all tilted ~26°) fit with centres in the box of side `3 − √2`, min gap `+0.0156` (float witness, scratch
`zwit.py`, 10-05; an earlier "−0.28" was a penalty-optimiser failure).  Locally `Z <= 4` needs box side ~1.5, h ~1.25.
So the count needs a *joint* lemma: a central pinwheel reaches toward the walls and should cost E/K squares.
Other cases `ΣK <= 3` need `Z <= ΣK − 1` too: corner-occupancy lemmas.

## h vs leaves

* h* = (1+√2)/2 exactly: three floor squares + a 45° square dipping its vertex into a gap g -> 1 has top centre at
  1 − g/2 + √2/2 -> h* as L -> 4.  Their h is 7e-7 below h*: that's why 8.2M leaves.
* **Lines alone: one line, at y = √2 − 1/2, proves the strip lemma for h <= (3√2 − 2)/2 = 1.12132** (every pose with
  centre in [hw, h] has chord >= 1 there, by `B(θ) >= (√2−1)/2`); the 1-D density LP over all line weightings gives the
  same threshold (scratch `rho.py`, `thr.py`).  Short, Lean-able (K1 band + our B1).
* Between 1.1213 and h*: needs branching (or 2-D).  Which h a consumer needs is now open again (`Z <= 4` gone).

## Lean: single-line strip lemma

`lean/Sqpack/ChordLine.lean` (2026-10-05), imported from `Sqpack.lean`, axiom lines added to `lean/Axioms.lean`.

* **`SquarePacking.wall_strip_le_three_line`**: for `t <= 4`, if `N` closed unit squares `sq (ctr i) (ang i) 1`
  lie in `box t = [0,t]^2`, are **pairwise disjoint as closed sets**, and each has centre height
  `(ctr i).2 <= (3√2 − 2)/2` (≈ 1.12132), then `N <= 3`.
* **`wall_strip_le_three_line_of_packing`**: the same for side-`L > 1` squares with disjoint interiors (B2 form).
* Key lemma `exists_chord_line`: such a square contains a closed segment of length 1 on `y = √2 − 1/2`.
* Boundary semantics: same as `wall_strip_le_three`, i.e. `t <= 4` with closed disjointness (no `t < 4` version needed).
* Proof reuses `exists_unit_interval_subset_slab`, `half_height_le_centre`, the slab sign symmetries and
  `no_four_spread` from `Chord.lean`.  The only new piece is the inequality.  With `u = |c|+|s| ∈ [1,√2]` and
  `2|c||s| = u² − 1`, it reduces to `(√2 − u)(u − (2 − √2)) >= 0` (below the line) and `(√2 − u)(√2 − 1 + u) >= 0`
  (above), both closed by `nlinarith`.
* `lake build Sqpack.ChordLine` is OK (~7 s warm), 0 sorries, `#print axioms` gives `[propext, Classical.choice, Quot.sound]`
  for both theorems.  The full `lake build` was not run because the tree is shared.

## Outcome, 2026-10-06 (three agents; details in `census/RESULTS.md`, `strip-demo/PROOF.md`)

* **Lean:** `wall_strip_le_three_line` (`lean/Sqpack/ChordLine.lean`), h = (3√2−2)/2, closed t <= 4.  Full `lake build`
  not yet run (shared tree).
* **Branching demo:** h = 1.13 by per-pair lines ("chain of skewed widths"): 26 boxes, 234 exact leaves, min chain
  1000357/250000.  Rechecked + controls rerun 10-06.  Hand lemmas (ordering, chain sum) checked by reading, not Lean.
  Size curve: 1.15 → 462 leaves, 1.18 → 4,505, 1.19 → ~1e5 (unmerged); exact counterexample at 1.2081, so h* is sharp.
* **Census [measured]: the counting route partitions the zero set, it doesn't shrink it.**  With 12 squares the
  count forces z >= k, so it prunes nothing by itself.  18 occupancy classes reach margin 0 at every h (rule: z <= 4 and
  <= 2 edge squares per wall).  Two candidate joint lemmas with room: (A) z >= 5 ⇒ margin <= −0.034; (B) 3 edge squares
  on one wall ⇒ <= −0.050.  Neither touches the 18 zero-margin classes, which each need the leaf-A rigidity argument.
  Witness margins rechecked independently 10-06.
* 10-06 follow-up [measured]: all 80 margin-0 census witnesses (18 classes) contain a wall-to-wall touching path of
  exactly axis-parallel squares (tilted squares up to 44.7° ride along).  Same picture as `BANDCUT_K.md` §0 and
  `RANK8.md`: occupancy class does not change the per-configuration problem (chain existence at δ >= 0 = obstruction X,
  `notes/n12-gap.md` §4.11).  Only 4 witnesses per class except the equality class (12): optimiser outputs, not a
  plateau sample.
