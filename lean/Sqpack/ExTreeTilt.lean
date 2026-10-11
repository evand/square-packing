import Sqpack.LebMassBridge

/-!
# `ValidTilt` from a pose tree

The tilted region `ValidTilt` (centres in `[0, m/2]²`, `0 < u`, `u² + 2u ≤ 1`) is what a `CovT` root
over `[0, m/2]² × (0, Um/R]` (`Um/R ≥ 1/2`) states: through `CovT.toLebMass` (`LebMassBridge.lean`) it
is a `LebMass.CovT` box over the root, and `LebMass.validTilt_of_covT` closes it.  With main's Lemma Z
(`Bentz.validAxis7`, `Bentz.valid7_of_tilt`) this replaces `cover_of_tree` and the generated face
`θ = 0` (the former `V7/Face.lean`).
-/

open MeasureTheory

namespace SquarePacking

namespace ZMTreeM

open BoxTree ZMTree

/-- **`ValidTilt` from a pose tree**: a `CovT` root over the D4 root box `[0, Mq/(2D)]² × (0, Um/R]`
with `Um/R ≥ 1/2` gives `ValidTilt (Mq/D)` for the mixed cover. -/
theorem validTilt_of_tree (D S Mq R Um W : ℕ) (pts : PTree) (segs : STree) (rects : List RectE)
    (hD : 0 < D) (hS : 0 < S) (hR : 0 < R) (hUm : R ≤ 2 * Um)
    (hcov : CovT D S Mq R W pts.toList segs.toList rects 0 (Mq * S - Mq * S / 2) 0
      (Mq * S - Mq * S / 2) 0 Um)
    (hW : 0 < W) :
    ValidSplit.ValidTilt ((Mq : ℝ) / D) (mcoverP D W pts segs rects).measure := by
  have h := hcov.toLebMass hW
  have hDr : (0 : ℝ) < D := by exact_mod_cast hD
  have hSr : (0 : ℝ) < S := by exact_mod_cast hS
  have hRr : (0 : ℝ) < R := by exact_mod_cast hR
  have hhalf : (Mq : ℝ) / D / 2 ≤ (((Mq * S - Mq * S / 2 : ℕ) : ℝ)) / (D * S) := by
    have h2 : Mq * S ≤ 2 * (Mq * S - Mq * S / 2) := by omega
    have h2r : ((Mq * S : ℕ) : ℝ) ≤ 2 * ((Mq * S - Mq * S / 2 : ℕ) : ℝ) := by exact_mod_cast h2
    rw [div_div, div_le_div_iff₀ (by positivity) (by positivity)]
    push_cast at h2r
    nlinarith
  have hU : (1 : ℝ) / 2 ≤ (Um : ℝ) / R := by
    rw [div_le_div_iff₀ (by norm_num) hRr]
    have : (R : ℝ) ≤ 2 * Um := by exact_mod_cast hUm
    linarith
  apply LebMass.validTilt_of_covT
  intro c u hx0 hx1 hy0 hy1 hu0 hu1 hu hu' hsub
  exact h c u (by simpa using hx0) (le_trans hx1 hhalf) (by simpa using hy0) (le_trans hy1 hhalf)
    hu0 (le_trans hu1 hU) hu hu' hsub

end ZMTreeM

end SquarePacking
