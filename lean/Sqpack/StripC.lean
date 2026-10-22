import Sqpack.StripS
import Sqpack.StripOdd
import Sqpack.StripZigzag

/-!
# (C') and the odd-`n` strip theorem for `1.87 ≤ H < 2`

`kap_C`: `kap ≥ 2λ - 1 + β||m| - |m'||` (`λ = 2√2 - H`, `β = 1 - H/2`) for `1.87 ≤ H < 2`.  Same-sign slopes: `kap_same`
(`2λ - 1 ≤ 1`).  Opposite slopes: per lines, `T m + T'(-y) ≥ 2λ - 1`; with both lines along edges this is
`SC (W + V_y) = p + C² - S² - C² H ≥ (4√2 - 1 - 2H) SC`, linear in `H`; at `H = 2` it is `p - 1 ≥ (4√2 - 5) SC`, at
`H = 1.87` a one-variable inequality with slack ≈ 6e-4 (at `tan t = 2`), certified by six product terms
(`chainn/certlp.py`).  It fails below `H ≈ 1.8685` (3-step prefixes beat the one-step potential there).

`strip_odd_iff`: for odd `n ≥ 3` and `1.87 ≤ H < 2`, `n` unit squares fit in `L × H` iff `L ≥ min(n, F_n)`,
`F_n = (n - 2)(2√2 - H) + 4 - H` (the 45° zigzag).
-/

namespace SquarePacking

/-- The one-variable inequality at `H = 1.87` (`4√2 - 4.74 ≤ 0.9172`). -/
lemma C_187 {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) :
    0 ≤ C + S - 87 / 100 * C ^ 2 - S ^ 2 - 2293 / 2500 * (S * C) := by
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have a1 := mul_nonneg hC (sq_nonneg (S - 2 * C))
  have a2 := mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hp1)
  have a3 := mul_nonneg (sub_nonneg.2 hS1) (sq_nonneg (S - 2 * C))
  have a4 := mul_nonneg (sub_nonneg.2 hp1) (sq_nonneg (C - 4472 / 10000))
  have a5 := mul_nonneg (sub_nonneg.2 hp1) (sq_nonneg (S - 8944 / 10000))
  have e1 : C * (C ^ 2 + S ^ 2) = C := by rw [hCS]; ring
  have e2 : S * (C ^ 2 + S ^ 2) = S := by rw [hCS]; ring
  nlinarith

/-- Both lines along edges: `(4√2 - 1 - 2H) SC ≤ p + C² - S² - C² H` for `1.87 ≤ H ≤ 2`. -/
lemma C_aligned {C S H : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hH1 : 1.87 ≤ H) (hH2 : H ≤ 2) :
    (4 * Real.sqrt 2 - 1 - 2 * H) * (S * C) ≤ C + S + C ^ 2 - S ^ 2 - C ^ 2 * H := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hSC : 0 ≤ S * C := mul_nonneg hS hC
  have hC1 : C ≤ 1 := by nlinarith
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  -- at H = 1.87
  have g0 : (4 * Real.sqrt 2 - 1 - 2 * 1.87) * (S * C) ≤ C + S + C ^ 2 - S ^ 2 - C ^ 2 * 1.87 := by
    have := C_187 hC hS hCS
    nlinarith [mul_le_mul_of_nonneg_right (by linarith : 4 * Real.sqrt 2 - 4.74 ≤ 2293 / 2500) hSC]
  -- at H = 2: p - 1 = 2SC/(p + 1) ≥ (4√2 - 5) SC
  have g2 : (4 * Real.sqrt 2 - 1 - 2 * 2) * (S * C) ≤ C + S + C ^ 2 - S ^ 2 - C ^ 2 * 2 := by
    have hq : (C + S - 1) * (C + S + 1) = 2 * (S * C) := by linear_combination hCS
    have hp2 : C + S ≤ 1.4143 := by nlinarith [sq_nonneg (C - S)]
    nlinarith [mul_le_mul_of_nonneg_left (by linarith : C + S + 1 ≤ 2.4143) (by linarith : (0:ℝ) ≤ C + S - 1)]
  -- linear in H
  have e : ∀ h : ℝ, C + S + C ^ 2 - S ^ 2 - C ^ 2 * h - (4 * Real.sqrt 2 - 1 - 2 * h) * (S * C) =
      (C + S + C ^ 2 - S ^ 2 - (4 * Real.sqrt 2 - 1) * (S * C)) + h * (2 * (S * C) - C ^ 2) := fun h => by ring
  have l0 := e 1.87; have l2 := e 2; have lH := e H
  have : 0 ≤ (2 - H) * (C + S + C ^ 2 - S ^ 2 - C ^ 2 * 1.87 - (4 * Real.sqrt 2 - 1 - 2 * 1.87) * (S * C)) +
      (H - 1.87) * (C + S + C ^ 2 - S ^ 2 - C ^ 2 * 2 - (4 * Real.sqrt 2 - 1 - 2 * 2) * (S * C)) :=
    add_nonneg (mul_nonneg (by linarith) (by linarith)) (mul_nonneg (by linarith) (by linarith))
  nlinarith

/-- **(C')** for `1.87 ≤ H < 2`. -/
theorem kap_C {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 1.87 ≤ H) (hH2 : H < 2) (hfit : C + S ≤ H) :
    2 * (2 * Real.sqrt 2 - H) - 1 + (1 - H / 2) * |(|m| - |m'|)| ≤ kap C S H m m' := by
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hlam : 2 * (2 * Real.sqrt 2 - H) - 1 ≤ 1 := by linarith
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  -- same sign (both ≥ 0)
  have same : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → 0 ≤ m → 0 ≤ m' →
      2 * (2 * Real.sqrt 2 - H) - 1 + (1 - H / 2) * |(|m| - |m'|)| ≤ kap C S H m m' := by
    intro C S m m' hC hS hCS hm hm'
    rw [abs_of_nonneg hm, abs_of_nonneg hm']
    rcases le_total m' m with h | h
    · rw [abs_of_nonneg (by linarith)]; have := kap_same (H := H) hC hS hCS hm' h; linarith
    · rw [abs_of_nonpos (by linarith), kap_symm]; have := kap_same (H := H) hC hS hCS hm h; linarith
  -- opposite, m ≥ 0 ≥ m' = -y
  have opp : ∀ {C S m y : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H → 0 ≤ y → y ≤ m →
      2 * (2 * Real.sqrt 2 - H) - 1 + (1 - H / 2) * (m - y) ≤ kap C S H m (-y) := by
    intro C S m y hC hS hCS hfit hy hym
    have hm : 0 ≤ m := hy.trans hym
    have hC1 : C ≤ 1 := by nlinarith
    have hS1 : S ≤ 1 := by nlinarith
    have hW0 := line_W hC hS hCS hm
    have hLy0 := line_Ly hC hS hCS hH2.le hy
    obtain ⟨k, hk⟩ : ∃ k : ℝ, k = (C + S) / 2 := ⟨_, rfl⟩
    obtain ⟨Tm, hTm⟩ : ∃ Tm : ℝ, Tm = (C + m * S + |m * C - S|) / 2 - m * (1 - (C + S) / 2) := ⟨_, rfl⟩
    obtain ⟨Ty, hTy⟩ : ∃ Ty : ℝ, Ty = (|C - y * S| + (y * C + S)) / 2 - y * (H - 1 - (C + S) / 2) := ⟨_, rfl⟩
    rw [← hTm, ← hk] at hW0
    rw [← hTy, ← hk] at hLy0
    have hkap : kap C S H m (-y) = Tm + Ty + (1 - H / 2) * (m - y) := by
      unfold kap
      rw [show m - -y = m + y by ring, abs_of_nonneg (by linarith : 0 ≤ m + y),
        abs_of_nonneg (by positivity : 0 ≤ C + m * S), show -y * C - S = -(y * C + S) by ring, abs_neg,
        abs_of_nonneg (by positivity : 0 ≤ y * C + S), show C + -y * S = C - y * S by ring, hTm, hTy]
      ring
    rw [hkap]
    -- it suffices: T m + T'(-y) ≥ 2λ - 1 and β|m - y| = β(m - y) up to the symmetric case
    have key : 2 * (2 * Real.sqrt 2 - H) - 1 ≤ Tm + Ty := by
      rcases eq_or_lt_of_le hS with hS0 | hSp
      · subst hS0
        have : C = 1 := by nlinarith
        subst this
        rw [hTm, hTy, abs_of_nonneg (by linarith : (0:ℝ) ≤ m * 1 - 0), abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - y * 0)]
        nlinarith
      rcases eq_or_lt_of_le hC with hC0 | hCp
      · subst hC0
        have : S = 1 := by nlinarith
        subst this
        rw [hTm, hTy, abs_of_nonpos (by linarith : m * 0 - 1 ≤ (0:ℝ)), abs_of_nonpos (by linarith : (0:ℝ) - y * 1 ≤ 0)]
        nlinarith
      rcases hLy0 with hy1 | hy1
      · have hg : C * k + S * (S - 1) + C * k - C = (1 - C) * (1 - S) := by rw [hk]; linear_combination hCS
        have h1 := mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)
        have h2 := mul_le_mul_of_nonneg_left hy1 hCp.le
        have : C * 1 ≤ C * (Tm + Ty) := by linarith
        have := (mul_le_mul_iff_of_pos_left hCp).1 this
        linarith
      · have hid : S * (C * k + S * (S - 1)) + C * (S * k + C * (C + 1 - H)) = C + S + C ^ 2 - S ^ 2 - C ^ 2 * H := by
          rw [hk]; linear_combination (C + S) * hCS
        have hWS := mul_le_mul_of_nonneg_left hW0 hSp.le
        have hyC := mul_le_mul_of_nonneg_left hy1 hCp.le
        have ca := C_aligned hC hS hCS hH1 hH2.le
        have : S * C * (2 * (2 * Real.sqrt 2 - H) - 1) ≤ S * C * (Tm + Ty) := by
          have e : S * C * (Tm + Ty) = S * (C * Tm) + C * (S * Ty) := by ring
          rw [e]; linarith
        exact (mul_le_mul_iff_of_pos_left (mul_pos hSp hCp)).1 this
    linarith
  -- opposite slopes in either order (mirror), then all sign cases
  have opp' : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H → 0 ≤ m → m' ≤ 0 →
      2 * (2 * Real.sqrt 2 - H) - 1 + (1 - H / 2) * |(|m| - |m'|)| ≤ kap C S H m m' := by
    intro C S m m' hC hS hCS hfit hm hm'
    rw [abs_of_nonneg hm, abs_of_nonpos hm']
    rcases le_total (-m') m with h | h
    · rw [abs_of_nonneg (by linarith)]
      have := opp (y := -m') hC hS hCS hfit (by linarith) h
      rwa [neg_neg] at this
    · rw [abs_of_nonpos (by linarith), kap_mirror]
      have := opp (m := -m') (y := m) hS hC (by linarith) (by linarith) hm h
      linarith
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · exact same hC hS hCS hm hm'
  · exact opp' hC hS hCS hfit hm hm'
  · rw [kap_symm, abs_sub_comm]; exact opp' hC hS hCS hfit hm' hm
  · rw [kap_mirror]
    have := same (m := -m') (m' := -m) hS hC hSC (by linarith) (by linarith)
    rwa [abs_neg, abs_neg, abs_sub_comm] at this

/-- `(C')` in the one-sided form used by the walk argument. -/
theorem kapC_of {H : ℝ} (hH1 : 1.87 ≤ H) (hH2 : H < 2) : KapC H := by
  intro C S m m' hC hS hCS hfit
  have := kap_C (m := m) (m' := m') hC hS hCS hH1 hH2 hfit
  have hb : 0 ≤ 1 - H / 2 := by linarith
  have : (1 - H / 2) * (|m'| - |m|) ≤ (1 - H / 2) * |(|m| - |m'|)| :=
    mul_le_mul_of_nonneg_left (by rw [abs_sub_comm]; exact le_abs_self _) hb
  linarith

/-- **Odd `n` in a strip, `1.87 ≤ H < 2`, lower bound**: `L ≥ min(n, F_n)`, `F_n = (n - 2)(2√2 - H) + 4 - H`. -/
theorem strip_odd_lower {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hH1 : 1.87 ≤ H) (hH2 : H < 2)
    (hp : RectPacks n L H) : min (n : ℝ) (((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H) ≤ L := by
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have := strip_odd_of_C (kapC_of hH1 hH2) (by linarith) hH2 hn hp
  have e : (n : ℝ) * (2 * Real.sqrt 2 - H) + 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) =
      ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H := by ring
  rwa [e] at this

/-- **Odd `n` in a strip, exactly** (`n ≥ 3` odd, `1.87 ≤ H < 2`): `n` unit squares fit in `[0,L] × [0,H]` iff
`L ≥ n` or `L ≥ (n - 2)(2√2 - H) + 4 - H` (the row, or the 45° zigzag). -/
theorem strip_odd_iff {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hn3 : 3 ≤ n) (hH1 : 1.87 ≤ H) (hH2 : H < 2) :
    RectPacks n L H ↔ (n : ℝ) ≤ L ∨ ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H ≤ L := by
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  constructor
  · intro hp; exact min_le_iff.1 (strip_odd_lower hn hH1 hH2 hp)
  · rintro (h | h)
    · exact strip_row h (by linarith)
    · exact strip_zigzag hn hn3 (by linarith) hH2.le h

end SquarePacking
