# Brief: audit zm_mixed.py and the s(21) run (adversarial), 2026-09-27

Read `tasks/s21-finish/README.md`, `tasks/line-cover/FORMAT.md`, `search/ZM_MIXED.md` (all), `search/zm_mixed.py`,
`search/zm_mixed_test.py`, and as needed `search/ZEROMARGIN.md`, `search/zeromargin.py`, `search/RUNG2.md` §3–7, `S32_EXACT.md` §7
(D4 argument).  Cores `taskset -c 10-11` only.  Budget: one working day.  **Your job is to find what is wrong.**

1. **Proofs.**  Check every lemma (A/B/C use, S, T and Corollaries T/T′, L, L′, Corollary L, R/SPLIT, V, P, the D4 and
   `--full` symmetry arguments, the reduction in FORMAT.md) line by line.  Closed semantics, `θ = 0` vs `θ > 0`, admissibility
   clipping, open/closed interval endpoints, double counting between groups / lines / points / phantom, withheld points, cache
   inheritance of `L` and of the witness cache across children (can a child inherit a stale phantom weight?), rounding of `L`.
2. **Code vs proof.**  For each lemma, find the code, confirm it implements the stated hypotheses (e.g. Lemma L's untyped-interval
   containment, Bernstein positivity of `Den`, the sub-bin split of `w/2`, `σ` slacks, breakpoint enumeration incl. tails).  Any float
   that can certify (not just reject) is a bug.
3. **Adversarial tests.**  Write new ones the author did not: targeted covers designed to break a specific lemma (a hole at a
   germ pivot pose; a hole exactly at a sub-bin boundary; segments meeting at grid crossings; points on group lines at the
   crossing; tiny tilts `u ~ 1e-9`); float re-sampling of certified leaves of the real m5 run if a leaf dump exists (or re-run a
   few roots with dumps).
4. **The run.**  Confirm the 40,000 roots tile the whole D4 region `[0, 5/2]² × u ∈ [0, 1/2]`, the checker sha / settings used
   match the committed file, the cover file checked is the one whose total is claimed, and D4 invariance of the measure is exact.

Report: `search/ZM_MIXED_AUDIT.md`: findings graded **must-fix** (soundness) / **should-fix** / **nit**, each with evidence; a verdict
on the s(21) claim.  Do not edit `zm_mixed.py` yourself (propose patches in the note).  Commit your note/tests only.
