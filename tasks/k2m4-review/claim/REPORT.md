# k2m4-review / claim — REPORT (transcribed by the session from the agent's final message; the agent could not write files)

BREAKS: none.

GAP (Lean coverage, not truth)
1. `Valid9` quantifies over every pose in [0,9]²; the run certifies the D4 region (u ∈ (0,1/2]) + Lemma Z (θ = 0).  The step
   is argued only in a docstring (lean/Sqpack/Bentz4.lean:46–52); same for `Valid7` (Bentz.lean:64).  Lean pieces exist:
   MixedMeasure.d4_reduction_measure_u (l.164), d4InvM (l.470), pattern S21.lean:233–292.  Until closed, say "the run
   certifies Valid9 via the cover's D4 symmetry".  (Symmetry checked exactly by the reviewer.)  → task lean-valid-split Part A.

MINOR
2. "Friedman's F₄" is our notation; DS7 Conjecture 1: "if s(n²−k)=n then s((n+1)²−k)=n+1" (his k = our c).  Say "the
   deficit-4 case of Friedman's Conjecture 1".  Holds for all n (n = 3 hypothesis false: s(5) < 3; n = 4 conclusion
   s(21) = 5 true; s(12) irrelevant).  FRIEDMAN.md:20's "⇔" is only given s(21) = 5.
3. State trust per k: k = 5 Lean conditional on a hypothesis certified by two checkers; k = 6 Lean hypothesis-free; k = 7 two
   checkers, no Lean; k ≥ 8 Lean conditional on Valid9 (single implementation).  Second routes: k = 8 s(60), k = 9 wand125.
4. Cite wand125's s(77) = 9 (2026-10-01, certificates/k2m4_n77_L9; extends our s(60) cover; our checkers; not independently
   reviewed per its README).

Checked OK: Lean statement semantics vs FORMAT.md (closed squares/container, parametric segment fraction, density-1
polygon); no sorry/native_decide/axiom; Lean data = file's 2076 segments in order, all positive-mass unit grid segments off
the walls, Lebesgue [14/5,31/5]², total 81 − 4D < 77, exactly 16 doubled segments, D4-invariant (exact, check_claim.py);
checkers read doubled segments additively (zm_mixed.py:79–102, qx2_zm.py:268–276, zm_mixed.py:350–371); k = 5, 6, 7 match
the s21/s32/s45 bundles.
Literature: no s(k²−4) with k ≥ 4 known before 2026-09-26.  s(5) = 2 + 1/√2 < 3 (Göbel) so the family starts at k = 5;
Nagamochi only a general bound (proof incomplete, Karakuş arXiv:2609.37410); Kearney–Shiu, Bentz: k²−3 only; Karakuş:
k²−1 only; Ellsworth marks none of 21/32/45/60/77 proved; jlevy credits 21/32/45 to us.
Suggested LIT text: "No exact value of s(k²−4) for k ≥ 4 was known before 2026; ours: s(32) = 6 (09-26), s(21) = 5 and
s(45) = 7 (09-27), s(60) = 8 (09-28); wand125 proved s(77) = 9 (10-01, extending our s(60) cover with our checkers); the
deficit-4 case of Friedman's Conjecture 1 follows."
Not checked: run leaves, Lean rebuild, certificates/k2m4/, wand125 s(77) replay, Bentz preprint status.
