# Lean: s(6) = 3 (Kearney–Shiu 2002)

Goal: `minSide 6 = 3` in Lean, hypothesis-free, in the repo's spec (`lean/Sqpack`, `minSide`, `Packs`, `Spec.lean`
bridges).  It is the last k of `s(k²−3) = k` not in Lean (k = 4: `s13_eq_4`; k = 5 conditional on `S21CheckerCover`;
k ≥ 6 conditional on `ValidTilt7`).  Upper bound `s(6) ≤ 3` is the grid (`packs_grid`).

## Read first
- Kearney, Shiu, *Efficient packing of unit squares in a square*, EJC 9 (2002) #R14 (open access).  Our summary:
  `notes/proof-anatomy.md` §5.2 (7-point "green/red lattices", Lemma 1 arc argument, s(6) §3 with Lemmas 2–3 and
  inequalities (7)).  Conventions: their boxes are open, `notes/proof-anatomy.md` header and §7; ours: `notes/branch-semantics.md`,
  `lean/Sqpack/Spec.lean`.
- `certificates/unavoid13/ks7_rational_3.txt`: a rational 7-point unavoidable set for `[0,3]²` (p(3) ≤ 7), certified
  by `search/unavoid13_check.py`; the Lean ZMTree machinery (`notes/s13-casefree.md`, LADDER rung 2) checks
  "every closed unit square in the box contains a point" statements.  Use it if the K–S accounting tolerates
  rational points; K–S's own points involve √2.
- Also: `search/S6_SKELETON.md`, `notes/t3-*.md` (our failed independent routes; context only).

## Steps
1. Check external Lean work first (formal-conjectures, wand125, chelokot, Queuingtheorydotcom repos; mathlib):
   if s(6) = 3 is already formalised, report and stop.
2. Plan note `notes/lean-s6.md`: the K–S proof broken into Lean lemmas (geometric lemmas with u = tan(θ/2)
   parametrisation as elsewhere in the repo; which steps are certificate checks vs hand proofs), conventions
   reconciled, size estimate.  Stop and report if it is much bigger than a few sessions.
3. Formalise.  Reuse existing lemmas (chord lemma, ZMTree, `not_packs_of_measure` style reductions).

## Rules
Kernel reduction only (no `native_decide`, no new axioms, no `sorry` in delivered theorems); `#print axioms` = the
standard three; add to `lean/Axioms.lean`; default build must pass.  CPU: no pinning; keep total busy processes ≲ 16.
Commit on your worktree branch (messages end with the session's Co-Authored-By line); don't push or merge.  Update
`lean/LADDER.md`.  Report: statements, axioms, timings, branch + commits, anything in K–S that needed repair.
