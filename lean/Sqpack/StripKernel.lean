import Mathlib

/-!
# The strip slope-walk kernel: every step costs at least the zigzag rate

Strip of height `H < 2`, separating lines `x + m y = s` between consecutive squares (walls: `m = 0`).
A unit square with frame `(C, S)` (`C, S ≥ 0`, `C² + S² = 1`, half-extent `k = (C + S)/2`) between lines of
slopes `m` (left) and `m'` (right) contributes at least
`κ = (g m + g m')/2 - |m - m'| (H/2 - k)`,  `g m = |C + m S| + |m C - S|`
to the length (`search/rectfront/RECTFRONT.md`, 2026-10-22).  `kappa_ge_zigzag`: for `2√2 - 1 ≤ H ≤ 2` and
`C + S ≤ H` (the square fits), `κ ≥ 2√2 - H`, the cost per square of the 45° zigzag.  Equality only at 45°.

Proof: for `m ≥ m'`, bound the four absolute values by linear forms whose `m`, `m'` coefficients cancel
(weights `(2c - S)/C` and `(C - 2c)/S`, `c = H/2 - k`, both in `[-1, 1]`), leaving `2SC κ ≥ 2p - H`,
`p = C + S`; then `2p - H - (2√2 - H)(p² - 1) = (2√2 - H)(√2 - p)(p - p₂)` with `p₂ ≤ 1`.  Upright squares
(`SC = 0`) cost `≥ 1 ≥ 2√2 - H`.
-/

namespace SquarePacking

/-- The step cost of one square between separating lines of slopes `m` and `m'`. -/
noncomputable def kap (C S H m m' : ℝ) : ℝ :=
  (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) / 2 - |m - m'| * (H / 2 - (C + S) / 2)

lemma kap_symm (C S H m m' : ℝ) : kap C S H m m' = kap C S H m' m := by
  simp only [kap, abs_sub_comm m m']; ring

/-- Linear minorant of `|x|` with weight `w`, `|w| ≤ W`, scaled by `W`. -/
lemma lin_min {w W x : ℝ} (h : |w| ≤ W) : w * x ≤ W * |x| :=
  calc w * x ≤ |w * x| := le_abs_self _
    _ = |w| * |x| := abs_mul _ _
    _ ≤ W * |x| := mul_le_mul_of_nonneg_right h (abs_nonneg _)

/-- Near-upright case `H ≤ 2S`: weights `(C ± 2c)/S` give `κ ≥ 1/S`. -/
lemma kap_ge_one_of_le {C S H m m' : ℝ} (hC : 0 ≤ C) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH2 : H ≤ 2) (hfit : C + S ≤ H) (hHS : H ≤ 2 * S) (hmm : m' ≤ m) :
    1 ≤ kap C S H m m' := by
  have hS : 0 < S := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  have hC1 : 2 * C ≤ H := by nlinarith
  unfold kap
  rw [abs_of_nonneg (by linarith : 0 ≤ m - m')]
  set c := H / 2 - (C + S) / 2 with hc
  have b1 := lin_min (x := C + m * S) (abs_le.2 ⟨by linarith, by linarith⟩ : |C + 2 * c| ≤ S)
  have b3 := lin_min (x := C + m' * S) (abs_le.2 ⟨by linarith, by linarith⟩ : |C - 2 * c| ≤ S)
  have b2 : -(m * C - S) ≤ |m * C - S| := neg_le_abs _
  have b4 : -(m' * C - S) ≤ |m' * C - S| := neg_le_abs _
  -- S (|..|+..) - 2 S (m - m') c ≥ 2 (C² + S²) = 2
  have key : 2 ≤ S * (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) - 2 * S * ((m - m') * c) := by
    have e2 := mul_le_mul_of_nonneg_left b2 hS.le
    have e4 := mul_le_mul_of_nonneg_left b4 hS.le
    have : (C + 2 * c) * (C + m * S) + (C - 2 * c) * (C + m' * S) + S * -(m * C - S) + S * -(m' * C - S)
        - 2 * S * ((m - m') * c) = 2 := by linear_combination 2 * hCS
    linarith
  set K := (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) / 2 - (m - m') * c
  by_contra hK; push Not at hK
  have : S * K < S * 1 := mul_lt_mul_of_pos_left hK hS
  have : S * K * 2 = S * (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) - 2 * S * ((m - m') * c) := by
    simp only [K]; ring
  linarith

/-- The ordered case `m' ≤ m`. -/
lemma kap_ge_zigzag_of_le {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2) (hfit : C + S ≤ H) (hmm : m' ≤ m) :
    2 * Real.sqrt 2 - H ≤ kap C S H m m' := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2p : 0 ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
  have r2g : 1.414 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hp2 : C + S ≤ Real.sqrt 2 := (Real.le_sqrt (by linarith) (by norm_num)).2 (by nlinarith [sq_nonneg (C - S)])
  rcases le_or_gt H (2 * S) with hHS | hHS
  · linarith [kap_ge_one_of_le hC hCS hH2 hfit hHS hmm]
  rcases le_or_gt H (2 * C) with hHC | hHC
  · -- mirror: swap C, S and reflect slopes m ↦ -m (order reverses)
    have e : kap C S H m m' = kap S C H (-m') (-m) := by
      unfold kap
      rw [show -m' - -m = m - m' by ring, show |S + -m' * C| = |m' * C - S| by
            rw [← abs_neg]; ring_nf,
          show |-m' * S - C| = |C + m' * S| by rw [← abs_neg]; ring_nf,
          show |S + -m * C| = |m * C - S| by rw [← abs_neg]; ring_nf,
          show |-m * S - C| = |C + m * S| by rw [← abs_neg]; ring_nf]
      ring
    rw [e]
    linarith [kap_ge_one_of_le hS (by linarith) hH2 (by linarith) hHC (by linarith : -m ≤ -m')]
  -- tilted: 2S < H, 2C < H; weights (2c - S)/C and (C - 2c)/S
  have hSp : 0 < S := by
    by_contra h; have : S = 0 := by linarith
    subst this
    have : C < 1 := by linarith
    nlinarith
  have hCp : 0 < C := by
    by_contra h; have : C = 0 := by linarith
    subst this
    have : S < 1 := by linarith
    nlinarith
  unfold kap
  rw [abs_of_nonneg (by linarith : 0 ≤ m - m')]
  set c := H / 2 - (C + S) / 2 with hc
  have a1 : C + m * S ≤ |C + m * S| := le_abs_self _
  have a4 : -(m' * C - S) ≤ |m' * C - S| := neg_le_abs _
  have b2 := lin_min (x := m * C - S) (abs_le.2 ⟨by linarith, by linarith⟩ : |2 * c - S| ≤ C)
  have b3 := lin_min (x := C + m' * S) (abs_le.2 ⟨by linarith, by linarith⟩ : |C - 2 * c| ≤ S)
  have key : 2 * (C + S) - H ≤ S * C * (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) -
      2 * S * C * ((m - m') * c) := by
    have e1 := mul_le_mul_of_nonneg_left a1 (mul_nonneg hSp.le hCp.le)
    have e4 := mul_le_mul_of_nonneg_left a4 (mul_nonneg hSp.le hCp.le)
    have e2 := mul_le_mul_of_nonneg_left b2 hSp.le
    have e3 := mul_le_mul_of_nonneg_left b3 hCp.le
    nlinarith
  have quad : (2 * Real.sqrt 2 - H) * (2 * S * C) ≤ 2 * (C + S) - H := by
    have hSC : 2 * S * C = (C + S) ^ 2 - 1 := by nlinarith
    rw [hSC]
    have f1 : 0 ≤ (Real.sqrt 2 - (C + S)) * ((2 * Real.sqrt 2 - H) * (C + S) - (Real.sqrt 2 * H - 2)) := by
      apply mul_nonneg (by linarith)
      have hA : 0 ≤ 2 * Real.sqrt 2 - H := by linarith
      have h1 := mul_le_mul_of_nonneg_left hp1 hA
      have h2 := mul_le_mul_of_nonneg_right hH2 (by linarith : (0:ℝ) ≤ 1 + Real.sqrt 2)
      linarith
    linear_combination f1 + (2 * (C + S) - H) * r2
  have hSC0 : 0 < S * C := mul_pos hSp hCp
  set K := (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) / 2 - (m - m') * c
  by_contra hK; push Not at hK
  have h3 : S * C * K < S * C * (2 * Real.sqrt 2 - H) := mul_lt_mul_of_pos_left hK hSC0
  have h4 : S * C * K * 2 = S * C * (|C + m * S| + |m * C - S| + |C + m' * S| + |m' * C - S|) -
      2 * S * C * ((m - m') * c) := by simp only [K]; ring
  linarith

/-- **Every step costs at least the zigzag rate**: `κ ≥ 2√2 - H`. -/
theorem kap_ge_zigzag {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2) (hfit : C + S ≤ H) :
    2 * Real.sqrt 2 - H ≤ kap C S H m m' := by
  rcases le_total m' m with h | h
  · exact kap_ge_zigzag_of_le hC hS hCS hH1 hH2 hfit h
  · rw [kap_symm]; exact kap_ge_zigzag_of_le hC hS hCS hH1 hH2 hfit h

end SquarePacking
