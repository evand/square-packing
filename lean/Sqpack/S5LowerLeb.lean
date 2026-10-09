import Sqpack.CovMP

/-!
# `s(5) ≥ 2` from a Lebesgue cover: the toy end-to-end run of `CovMP`

The cover: Lebesgue measure on the container `[0,2]²` (one rectangle entry of density `1`), total
`4 < 5`.  Every closed unit square inside the container has mass exactly `1` (`LEB`, one leaf), so
five unit squares do not fit in a square of side `< 2`.  The bound is trivial (`s(5) = 2.707…`);
the point is the pipeline: `CovMP.of_leb`, `rectD4Check`, `le_minSide_mixedP`.
-/

namespace SquarePacking

open ZMTreeM

/-- **`s(5) ≥ 2`**, from the Lebesgue measure on the container, by the `LEB` leaf. -/
theorem s5_ge_2_leb : (2 : ℝ) ≤ minSide 5 := by
  have hcov : CovMP 1 1 2 2 1 PTree.leaf.toList STree.leaf.toList [(0, 0, 2, 2, 1)] 0
      (2 * 1 - 2 * 1 / 2) 0 (2 * 1 - 2 * 1 / 2) 0 1 :=
    CovMP.of_leb (by norm_num) (by norm_num) (by norm_num) (List.mem_singleton_self _)
      (by decide +kernel)
  have h := le_minSide_mixedP 1 1 2 2 1 1 PTree.leaf STree.leaf [(0, 0, 2, 2, 1)]
    (by norm_num) (by norm_num) (by norm_num) (by norm_num) (by decide +kernel)
    (by decide +kernel) (by decide +kernel) (by decide +kernel) (by decide +kernel) hcov 5 3
    (by decide +kernel) (by norm_num)
  have e : ((2 : ℕ) : ℝ) / ((1 : ℕ) : ℝ) = 2 := by norm_num
  rw [e] at h
  exact h

end SquarePacking
