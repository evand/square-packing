import Sqpack.StripKernel

/-!
# Strip slope walk, odd `n`: per-square inequalities behind the prefix potentials

Notation (strip height `H`, `2√2 - 1 ≤ H ≤ 2`): `λ = 2√2 - H` (zigzag cost per square), `β = 1 - H/2`,
`a = 1 - λ` (excess of a grid step).  Potentials (`search/rectfront/RECTFRONT.md`, 2026-10-22):
`u₁ m = a + β|m|`, `u₀ m = min (2a + β|m|) (a + β(2 - |m|))`.

* `line_T`: the per-line part `T x = g x / 2 - |x| (1 - k)` of a square's step cost is smallest when the line runs
  along the square's edge; for the two edge lines `T (tan t) + T (-cot t) = 2/(1 + cos t + sin t)`.
* `kap_B`: every step has `kap ≥ 2(√2 - 1) + β(|m| + |m'|)`, i.e. excess `≥ β(|m| + |m'| - 2)`.
-/

namespace SquarePacking

/-- **Per-line bound** (`x ≥ 0`): `2C · T x ≥ 1 - 2S + SC + S²`, equality for the line along the edge
(`x = tan t`). -/
lemma line_T {C S x : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hx : 0 ≤ x) :
    1 - 2 * S + S * C + S ^ 2 ≤
      2 * C * ((|C + x * S| + |x * C - S|) / 2 - x * (1 - (C + S) / 2)) := by
  have hS1 : S ≤ 1 := by nlinarith
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  rw [abs_of_nonneg (by positivity : 0 ≤ C + x * S)]
  have b := lin_min (x := x * C - S) (abs_le.2 ⟨by linarith, by linarith⟩ : |2 - C - 2 * S| ≤ C)
  nlinarith

/-- The per-line bound for a line of negative slope (mirror). -/
lemma line_T_neg {C S x : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hx : x ≤ 0) :
    1 - 2 * C + S * C + C ^ 2 ≤
      2 * S * ((|C + x * S| + |x * C - S|) / 2 - (-x) * (1 - (C + S) / 2)) := by
  have := line_T (x := -x) hS hC (by linarith) (by linarith)
  rw [show |S + -x * C| = |x * C - S| by rw [← abs_neg]; ring_nf,
    show |-x * S - C| = |C + x * S| by rw [← abs_neg]; ring_nf] at this
  linarith

/-- **A square between its two edge lines**: the two per-line parts add up to at least `2(√2 - 1)`
(exactly `2/(1 + cos t + sin t)` when both lines run along edges). -/
lemma two_lines {C S m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hm : 0 ≤ m) (hm' : m' ≤ 0) :
    2 * (Real.sqrt 2 - 1) ≤ ((|C + m * S| + |m * C - S|) / 2 - m * (1 - (C + S) / 2)) +
      ((|C + m' * S| + |m' * C - S|) / 2 - (-m') * (1 - (C + S) / 2)) := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2p : 0 ≤ Real.sqrt 2 := Real.sqrt_nonneg 2
  have hp2 : C + S ≤ Real.sqrt 2 := (Real.le_sqrt (by linarith) (by norm_num)).2 (by nlinarith [sq_nonneg (C - S)])
  set T1 := (|C + m * S| + |m * C - S|) / 2 - m * (1 - (C + S) / 2)
  set T2 := (|C + m' * S| + |m' * C - S|) / 2 - (-m') * (1 - (C + S) / 2)
  have h1 : 1 - 2 * S + S * C + S ^ 2 ≤ 2 * C * T1 := line_T hC hS hCS hm
  have h2 : 1 - 2 * C + S * C + C ^ 2 ≤ 2 * S * T2 := line_T_neg hC hS hCS hm'
  -- (√2 - 1)(1 + p) ≤ 1
  have hq : (Real.sqrt 2 - 1) * (1 + (C + S)) ≤ 1 := by nlinarith
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · -- upright: T1 + T2 = 1
    subst hS0
    have hC1 : C = 1 := by nlinarith
    subst hC1
    have e1 : T1 = 1 / 2 := by
      simp only [T1]; rw [abs_of_nonneg (by linarith), abs_of_nonneg (by linarith)]; ring
    have e2 : T2 = 1 / 2 := by
      simp only [T2]; rw [abs_of_nonneg (by linarith), abs_of_nonpos (by linarith)]; ring
    nlinarith
  rcases eq_or_lt_of_le hC with hC0 | hCp
  · subst hC0
    have hS1 : S = 1 := by nlinarith
    subst hS1
    have e1 : T1 = 1 / 2 := by
      simp only [T1]; rw [abs_of_nonneg (by linarith), abs_of_nonpos (by linarith)]; ring
    have e2 : T2 = 1 / 2 := by
      simp only [T2]; rw [abs_of_nonpos (by linarith), abs_of_nonpos (by linarith)]; ring
    nlinarith
  -- (1 + p)(S A + C B) = 4 S C
  have id : (1 + (C + S)) * (S * (1 - 2 * S + S * C + S ^ 2) + C * (1 - 2 * C + S * C + C ^ 2)) = 4 * S * C := by
    linear_combination (C + S) * (C + S - 1) * hCS
  have hSC : 0 < S * C := mul_pos hSp hCp
  -- 2SC (T1 + T2) ≥ S A + C B
  have h3 : S * (1 - 2 * S + S * C + S ^ 2) + C * (1 - 2 * C + S * C + C ^ 2) ≤ 2 * S * C * (T1 + T2) := by
    nlinarith [mul_le_mul_of_nonneg_left h1 hSp.le, mul_le_mul_of_nonneg_left h2 hCp.le]
  -- so (1 + p)(T1 + T2) ≥ 2
  have h4 : 2 ≤ (1 + (C + S)) * (T1 + T2) := by
    have := mul_le_mul_of_nonneg_left h3 (by linarith : (0:ℝ) ≤ 1 + (C + S))
    rw [id] at this
    by_contra hc; push Not at hc
    have : 2 * S * C * ((1 + (C + S)) * (T1 + T2)) < 2 * S * C * 2 := mul_lt_mul_of_pos_left hc (by linarith)
    nlinarith
  by_contra hc; push Not at hc
  have : (1 + (C + S)) * (T1 + T2) < (1 + (C + S)) * (2 * (Real.sqrt 2 - 1)) :=
    mul_lt_mul_of_pos_left hc (by linarith)
  nlinarith

/-- `F1`: `1 - 2(√2-1) C - (2 - H) S > 0` (`(2(√2-1))² + (3 - 2√2)² < 1`). -/
lemma F1 {C S H : ℝ} (_hC : 0 ≤ C) (_hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H)
    (hH2 : H ≤ 2) : 0 < 1 - 2 * (Real.sqrt 2 - 1) * C - (2 - H) * S := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  set r := 2 * (Real.sqrt 2 - 1); set b := 2 - H
  have hb0 : 0 ≤ b := by simp only [b]; linarith
  have hb : b ≤ 3 - 2 * Real.sqrt 2 := by simp only [b]; linarith
  have hr : 0.8284 < r ∧ r < 0.8286 := ⟨by simp only [r]; linarith, by simp only [r]; linarith⟩
  -- (rC + bS)² ≤ (r² + b²)(C² + S²) < 1
  have cs : (r * C + b * S) ^ 2 ≤ r ^ 2 + b ^ 2 := by nlinarith [sq_nonneg (r * S - b * C)]
  have hsum : r ^ 2 + b ^ 2 < 1 := by nlinarith
  nlinarith [sq_nonneg (r * C + b * S - 1)]

/-- `F2`: `1 + SC - S - 2(√2-1) C ≥ 0`. -/
lemma F2 {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) :
    0 ≤ 1 + S * C - S - 2 * (Real.sqrt 2 - 1) * C := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  rcases le_total (2 * (Real.sqrt 2 - 1)) S with h | h
  · nlinarith [mul_nonneg hC (sub_nonneg.2 h)]
  · nlinarith [mul_le_mul_of_nonneg_right hC1 (sub_nonneg.2 h)]

/-- Same-sign slopes `0 ≤ m' ≤ m`: `kap ≥ 2(√2 - 1) + β(m + m')`.  By the kink of each line (`m C ≷ S`) the
excess is linear in `m, m'`; pushed to the vertices `(0,0)`, `(tan t, 0)`, `(tan t, tan t)` it is `F1` or `F2`. -/
lemma kap_B_same {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2) (hfit : C + S ≤ H) (hm' : 0 ≤ m') (hmm : m' ≤ m) :
    2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (m + m') ≤ kap C S H m m' := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hm : 0 ≤ m := hm'.trans hmm
  have f1 := F1 hC hS hCS hH1 hH2
  have f2 := F2 hC hS hCS
  set r := 2 * (Real.sqrt 2 - 1)
  unfold kap
  rw [abs_of_nonneg (by linarith : 0 ≤ m - m'), abs_of_nonneg (by nlinarith : 0 ≤ C + m * S),
    abs_of_nonneg (by nlinarith : 0 ≤ C + m' * S)]
  -- E = kap - r - β(m + m'), as a function of the kink cases
  have Cpos : ∀ x, S ≤ x * C → 0 < C := fun x h => by
    rcases eq_or_lt_of_le hC with h0 | h0
    · subst h0; have : S = 1 := by nlinarith
      linarith
    · exact h0
  rcases le_total S (m' * C) with h2 | h2
  · -- both lines above the kink: C E ≥ 1 - rC - (2-H) S
    have h1 : S ≤ m * C := h2.trans (mul_le_mul_of_nonneg_right hmm hC)
    rw [abs_of_nonneg (by linarith : 0 ≤ m * C - S), abs_of_nonneg (by linarith : 0 ≤ m' * C - S)]
    have hCp := Cpos m' h2
    have a1 : 0 ≤ (m - m') * C * (C + S - 1) := mul_nonneg (mul_nonneg (by linarith) hC) (by linarith)
    have a2 : 0 ≤ (m' * C - S) * (C + S + H - 2) := mul_nonneg (by linarith) (by linarith)
    have key : 0 ≤ C * ((C + m * S + (m * C - S) + (C + m' * S) + (m' * C - S)) / 2 -
        (m - m') * (H / 2 - (C + S) / 2) - r - (1 - H / 2) * (m + m')) := by nlinarith
    have := (mul_nonneg_iff_of_pos_left hCp).1 key
    linarith
  rcases le_total S (m * C) with h1 | h1
  · -- m above, m' below the kink
    rw [abs_of_nonneg (by linarith : 0 ≤ m * C - S), abs_of_nonpos (by linarith : m' * C - S ≤ 0)]
    have hCp := Cpos m h1
    have a1 : 0 ≤ (m * C - S) * (C + S - 1) := mul_nonneg (by linarith) (by linarith)
    rcases le_total 0 (H - 1 - C) with hh | hh
    · have a2 : 0 ≤ m' * C * (H - 1 - C) := mul_nonneg (mul_nonneg hm' hC) hh
      have key : 0 ≤ C * ((C + m * S + (m * C - S) + (C + m' * S) + -(m' * C - S)) / 2 -
          (m - m') * (H / 2 - (C + S) / 2) - r - (1 - H / 2) * (m + m')) := by nlinarith
      have := (mul_nonneg_iff_of_pos_left hCp).1 key
      linarith
    · have a2 : 0 ≤ (S - m' * C) * (C + 1 - H) := mul_nonneg (by linarith) (by linarith)
      have a3 : 0 ≤ (m * C - S) * (C + S + H - 2) := mul_nonneg (by linarith) (by linarith)
      have key : 0 ≤ C * ((C + m * S + (m * C - S) + (C + m' * S) + -(m' * C - S)) / 2 -
          (m - m') * (H / 2 - (C + S) / 2) - r - (1 - H / 2) * (m + m')) := by nlinarith
      have := (mul_nonneg_iff_of_pos_left hCp).1 key
      linarith
  · -- both below the kink
    rw [abs_of_nonpos (by linarith : m * C - S ≤ 0), abs_of_nonpos (by linarith : m' * C - S ≤ 0)]
    rcases eq_or_lt_of_le hC with hC0 | hCp
    · subst hC0
      have : S = 1 := by nlinarith
      subst this
      nlinarith
    have a1 : 0 ≤ (S - m * C) * (1 - S) := mul_nonneg (by linarith) (by linarith)
    rcases le_total 0 (H - 1 - C) with hh | hh
    · have a2 : 0 ≤ m' * C * (H - 1 - C) := mul_nonneg (mul_nonneg hm' hC) hh
      have key : 0 ≤ C * ((C + m * S + -(m * C - S) + (C + m' * S) + -(m' * C - S)) / 2 -
          (m - m') * (H / 2 - (C + S) / 2) - r - (1 - H / 2) * (m + m')) := by nlinarith
      have := (mul_nonneg_iff_of_pos_left hCp).1 key
      linarith
    · have a2 : 0 ≤ (S - m' * C) * (C + 1 - H) := mul_nonneg (by linarith) (by linarith)
      have key : 0 ≤ C * ((C + m * S + -(m * C - S) + (C + m' * S) + -(m' * C - S)) / 2 -
          (m - m') * (H / 2 - (C + S) / 2) - r - (1 - H / 2) * (m + m')) := by nlinarith
      have := (mul_nonneg_iff_of_pos_left hCp).1 key
      linarith

/-- **(B)**: every step has `kap ≥ 2(√2 - 1) + β(|m| + |m'|)`, i.e. excess `kap - λ ≥ β(|m| + |m'| - 2)`. -/
theorem kap_B {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2) (hfit : C + S ≤ H) :
    2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (|m| + |m'|) ≤ kap C S H m m' := by
  -- opposite signs, m ≥ 0 ≥ m'
  have opp : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H → 0 ≤ m → m' ≤ 0 →
      2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (|m| + |m'|) ≤ kap C S H m m' := by
    intro C S m m' hC hS hCS hfit hm hm'
    have := two_lines hC hS hCS hm hm'
    unfold kap
    rw [abs_of_nonneg hm, abs_of_nonpos hm', abs_of_nonneg (by linarith : 0 ≤ m - m')]
    linarith
  -- same sign, both ≥ 0
  have same : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H → 0 ≤ m → 0 ≤ m' →
      2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (|m| + |m'|) ≤ kap C S H m m' := by
    intro C S m m' hC hS hCS hfit hm hm'
    rw [abs_of_nonneg hm, abs_of_nonneg hm']
    rcases le_total m' m with h | h
    · exact kap_B_same hC hS hCS hH1 hH2 hfit hm' h
    · rw [kap_symm, add_comm m m']; exact kap_B_same hC hS hCS hH1 hH2 hfit hm h
  -- the mirror turns (m, m') into (-m', -m) and swaps the frame
  have mir : ∀ {m m' : ℝ}, 2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (|-m'| + |-m|) ≤ kap S C H (-m') (-m) →
      2 * (Real.sqrt 2 - 1) + (1 - H / 2) * (|m| + |m'|) ≤ kap C S H m m' := by
    intro m m' h; rw [kap_mirror]; rw [abs_neg, abs_neg] at h; linarith
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · exact same hC hS hCS hfit hm hm'
  · exact opp hC hS hCS hfit hm hm'
  · -- m ≤ 0 ≤ m': swap the two lines
    rw [kap_symm, add_comm |m| |m'|]; exact opp hC hS hCS hfit hm' hm
  · exact mir (same hS hC hSC (by linarith) (by linarith) (by linarith))

end SquarePacking
