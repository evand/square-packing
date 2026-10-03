# Brief: Lean for s(21) = 5 (mixed-measure reduction), 2026-09-27

Read `tasks/s21-finish/README.md`, `tasks/line-cover/FORMAT.md`, `notes/lean-s32.md`, `lean/Sqpack/Basic.lean`
(`packing_le_weight*`), `lean/Sqpack/D4.lean`, `lean/Sqpack/S32.lean`.  Cores `taskset -c 12-13`.  Budget: one working day.

1. Generalise the reduction from finite weighted point sets to the mixed measures of FORMAT.md: finitely many point masses,
   segments carrying mass uniformly by length, convex polygons carrying mass uniformly by area (a clean route: any finite
   Borel measure `μ` with `μ(sq c θ 1) ≥ 1` for every admissible closed unit square gives `n ≤ μ(C)` for a packing at side `L > 1`
   scaled into `C`, by disjointness of the shrunken closed squares; then instantiate `μ` as the sum of Dirac, 1-D Hausdorff /
   pushforward-of-Lebesgue on segments, and restricted Lebesgue on polygons, and show its total is the file's total).
2. D4 symmetry for such measures (the analogue of what `D4.lean` does for point sets): if `μ` is invariant under the square's
   symmetries, the cover hypothesis on the D4 fundamental region implies it everywhere.
3. `s21_eq_five_of_checker`: `s(21) = 5` from one clearly stated computational hypothesis (the checker's statement on the D4
   region for this cover), mirroring `s32_eq_six_of_checker`, plus `s(21) <= 5` from the 5×5 tiling.  No `sorry`; list axioms.
Note: `notes/lean-s21.md`.  Commit own files, never push.
