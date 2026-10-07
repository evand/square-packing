# Lean: uniform-area (Lebesgue) mass primitives — Lemma U and Lemma K

Why now: every insertable certificate must be Lebesgue on its slab-crossing square (S2 lemma, `search/S2_INSERTABLE.md`,
in progress), and `ValidTilt7/9` (the last unformalised steps of k²−3 / k²−4) need the same primitives.  So Lean
leaf tests for tilted squares meeting a uniform-area region are on every route.

Read: `search/QUADRANT_EXACT.md` §4.1–4.4 (LEB = Lemma U, CAP = Lemma K incl. the width lemma and the tangent-cap
refinement, chord-end formulas §4.3; Lemma E is a later batch), `search/qx2_zm.py` (~l.939–1060, the tests as run),
`lean/LADDER.md` (ValidSplit section, segments batch 1), `notes/lean-bentz-reduction.md`, `notes/lean-segments.md`,
`lean/Sqpack/` (`MixedCover`, `CovM`, `segMeasure`, the box/pose-bin machinery, Lemmas P/S/T/L).

## Steps
1. **Design note first** (`notes/lean-leb-mass.md`): how a pose box (centre box × angle bin, with the `ŵ` bound on
   `cos θ + sin θ`) is represented in the existing ZMTree/BoxTree leaf machinery; exact Lean statements of U and K as
   leaf tests; what is reused (segment lemmas, chord ends) vs new (width lemma, area of `Q ∩ {y < a}` bound,
   `λ(Q) = 1` for a rotated unit square); estimate of size.  Stop and report if the design shows it is much bigger
   than a few sessions.
2. **Lemma U** (square inside the Lebesgue square ⇒ mass ≥ 1) as a sound leaf test over a pose box.
3. **Lemma K** (caps): the width lemma, then the crude test, then the tangent-cap refinement if time allows.
4. A smoke test: a small toy certificate whose leaves close by U/K (e.g. a few boxes of the k²−3 box-7 run
   `search/qx2_data/L4_k02_box7.txt`), so the tests are exercised by `decide`/kernel reduction, not only stated.

## Rules
Kernel reduction only (no `native_decide`, no new axioms, no `sorry` in delivered theorems); `#print axioms` = the
standard three; add to `lean/Axioms.lean`; default build must still pass.  CPU: no core pinning; keep the machine's
total busy processes ≲ 16 (check `uptime` / load before big builds; another agent is running LPs).  Commit on your
worktree branch (message ends with the session's Co-Authored-By line); don't push or merge.  Update `lean/LADDER.md`
(a section).  Report: exact statements, axioms, timings, branch + commits, design-note conclusions, size estimate for
the rest (Lemma E, ValidTilt7).
