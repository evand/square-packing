# k2m4-review: adversarial review of the s(k²−4) = k certificate (2026-10-03)

**Claim under review** (paths relative to `~/math/square-packing/public/s12`): the box cover
`search/qx2_data/K4_k008_box9.txt` (2076 unit segments on the 1/5-grid + Lebesgue `[14/5, 31/5]²`, total `81 − 4D`,
`D = 214770225571/200000000000`) satisfies **Valid9**: every closed unit square in `[0,9]²`, any centre, any angle, has mass
≥ 1.  Then Lean `bentz4_of_valid9` (`lean/Sqpack/Bentz4.lean`, kernel-checked) gives s(k² − 4) = k for all k ≥ 8, and our
bundles s21, s32, s45 give k = 5, 6, 7 (Friedman's F₄).  Certificate of record: `qx2_zm.py` (sha 6294052a…, the same code
as the reviewed k2m3 certificate) run `runs/qx2_k4x_k008/qxzm_full.{out,jsonl}` (VERIFIED-D4, u > 0, 16,200 roots, 0
uncertified); θ = 0 by Lemma Z (`runs/qx2_k4x_k008_xloop.out`, axis face min 1); exact germ scan (`germscan.out`).
Write-up: `search/K2M4_MARGIN.md` §4, `search/QUADRANT_EXACT.md` (the k2m3 method), `notes/lean-k2m4-reduction.md`.
The k2m3 review (`~/math/square-packing/private/s12/tasks/k2m3-review/*/REPORT.md`) covered the checker code; **focus on
what is new here**, but anything in the code that only bites at R = 3 / box 9 / this cover is in scope.

**What is new vs k2m3:** R = w = 3, box 9; margin settings κ = 0.08, θ₀ = 3°; the corner module puts mass on the band's
end lines x = 3, y = 3 where the profile's phase-0 line also lies (16 unit segments listed twice in the box file); far more
Lemma E leaves (EXACT 20,577, EXACT0 1,424, EXACT45 1,392; PIECE 84,152); 4.3× the expected CPU, concentrated in roots
5,400–12,000.

**Your job: try to break it.** Hostile referee; default posture: wrong somewhere.  Concrete findings are worth most (an
exact pose with mass < 1; a leaf whose certification is unjustified; a case not covered).  "No problem found" is fine if
you say exactly what you checked and how, and what you did not.  Derive/check statements yourself before reading proofs;
design your own tests.  Known historical bug class: tile-germ limits (θ → 0⁺ with lines on the square's edges) slightly
below 1, invisible to float sampling.  Tiny tilts, one-sided limits, degenerate boxes, U's boundary lines, D4 bookkeeping,
duplicated segments, "dropped because float says empty".

**Rules.** Read-only on `public/` (no edits, no commits, no state-changing git).  Scratch + report in
`~/math/square-packing/private/s12/tasks/k2m4-review/<area>/`.  Exact arithmetic (`fractions.Fraction`) for any violation
claim.  Only your assigned cores (`taskset -c`); ≤ 2 CPU-h total unless a concrete lead needs more; ≤ 8 GB RAM.  Never
execute anything under `~/math/_untrusted-third-party/`.  Report `<area>/REPORT.md` and final message ≤ 400 words: findings
ranked **BREAKS / GAP / MINOR** with evidence and file:line; then "checked, OK" with how; then "not checked".

**Areas.**
* `break-it` (cores 4–5): hunt for violating poses: exact mass at germs/tiny tilts near the doubled end lines, the corner
  module ↔ band junction, the Lebesgue boundary at 14/5, the wall; your own exact searches (not the shipped germ scan).
* `leaves` (core 6): the record.  Root set vs the D4 region of [0,9]² (is every pose of Valid9 reached after the D4 and
  θ ↦ π/2 − θ reductions?); sample leaves of each kind (bias to EXACT/PIECE and the slow roots 5,400–12,000) and
  re-justify them independently (recompute the exact bound your own way where feasible); Lemma E's hypotheses at these
  slopes (u ≤ 1/2, exact-from depth 3); the axis (θ = 0) face and its one-sided limits for this cover.
* `claim` (core 7): the statement chain.  Lean `bentz4_of_valid9` statement and `Valid9` definition vs the box file and
  FORMAT.md semantics (closed squares, mass by parametric fraction / area); whether k = 5, 6, 7 really follow from the
  s21/s32/s45 bundles as stated; F₄ wording vs Friedman's conjecture; **literature**: what was known for s(k² − 4) = k
  before (Bentz, Kearney–Shiu, Nagamochi (note its proof gap), Friedman's page, Karakuş, jlevy/squares, wand125's s(77)),
  with sources; anything we'd overclaim.  Light compute.
