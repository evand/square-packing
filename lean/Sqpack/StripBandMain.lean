import Sqpack.StripBandCells
import Sqpack.StripC

/-!
# Odd `n` in a strip, `H*(9) ≤ H < 2`

`strip_odd_band`: on the band `H*(9) = (14√2 - 5)/8 ≤ H ≤ 1.87` the tangent potentials `bandU0`, `bandU1`
(`StripBand`) are admissible (`band_potOK`), so `strip_odd_pot` gives `L ≥ min(n, F_n)` for every odd `n`.
With `strip_odd_lower` (`H ≥ 1.87`): `strip_odd_lower'` on `H*(9) ≤ H < 2`.  Since `H*(7)`, `H*(9) ≥ H*(9)`,
monotonicity in `H` completes `n = 7, 9` (`strip_odd_iff_le9`: every odd `n ≤ 9`, every `1 ≤ H < 2`).

Per-square facts: opposite slopes from the generated cores `band_Sa`, `band_Sb`, `band_C0`, `band_C1`, `band_N`
(`StripBandCells`); same-sign slopes from `kap_same` and the Lipschitz bounds of the potentials.
-/

namespace SquarePacking

section lip
variable {H : ℝ}

lemma bandU0_neg (m : ℝ) : bandU0 H (-m) = bandU0 H m := by unfold bandU0; rw [abs_neg]
lemma bandU1_neg (m : ℝ) : bandU1 H (-m) = bandU1 H m := by unfold bandU1; rw [abs_neg]

/-- The lines of the potentials have slopes `β`, `-β`, `T₁` with `|T₁| ≤ β`, so the potentials are
`β`-Lipschitz in `|m|`. -/
lemma bandL_lip (hb : InBand H) (u v : ℝ) :
    bandL H u ≤ bandL H v + (1 - H / 2) * |u - v| := by
  obtain ⟨h9, h87⟩ := hb
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hb0 : 0 ≤ 1 - H / 2 := by linarith
  have d1 : (1 - H / 2) * (u - v) ≤ (1 - H / 2) * |u - v| := mul_le_mul_of_nonneg_left (le_abs_self _) hb0
  have d2 : (1 - H / 2) * (v - u) ≤ (1 - H / 2) * |u - v| := by
    rw [abs_sub_comm]; exact mul_le_mul_of_nonneg_left (le_abs_self _) hb0
  have d3 : (H / 16 - 1 / 10) * (u - v) ≤ (1 - H / 2) * |u - v| := by
    have : |H / 16 - 1 / 10| ≤ 1 - H / 2 := abs_le.2 ⟨by linarith, by linarith⟩
    exact lin_min this
  unfold bandL bandLz bandT
  apply max_le
  · have : (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - u) ≤
        ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - v)) + (1 - H / 2) * |u - v| := by linarith
    exact this.trans (by gcongr; exact le_max_left _ _)
  · have : (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * (u - 8 / 7) ≤
        ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * (v - 8 / 7)) +
          (1 - H / 2) * |u - v| := by linarith
    exact this.trans (by gcongr; exact le_max_right _ _)

lemma bandU0_lip (hb : InBand H) (s t : ℝ) : bandU0 H s ≤ bandU0 H t + (1 - H / 2) * |s - t| := by
  have hb0 : 0 ≤ 1 - H / 2 := by linarith [hb.2]
  have hst : (1 - H / 2) * |(|s| - |t|)| ≤ (1 - H / 2) * |s - t| :=
    mul_le_mul_of_nonneg_left (abs_abs_sub_abs_le_abs_sub s t) hb0
  have d1 : (1 - H / 2) * (|s| - |t|) ≤ (1 - H / 2) * |(|s| - |t|)| :=
    mul_le_mul_of_nonneg_left (le_abs_self _) hb0
  have hL := bandL_lip hb |s| |t|
  unfold bandU0
  rcases min_choice (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |t|) (bandL H |t|) with h | h <;> rw [h]
  · exact le_trans (min_le_left _ _) (by linarith)
  · exact le_trans (min_le_right _ _) (by linarith)

lemma bandU1_lip (hb : InBand H) (s t : ℝ) : bandU1 H s ≤ bandU1 H t + (1 - H / 2) * |s - t| := by
  have hb0 : 0 ≤ 1 - H / 2 := by linarith [hb.2]
  have hst : (1 - H / 2) * |(|s| - |t|)| ≤ (1 - H / 2) * |s - t| :=
    mul_le_mul_of_nonneg_left (abs_abs_sub_abs_le_abs_sub s t) hb0
  have d1 : (1 - H / 2) * (|s| - |t|) ≤ (1 - H / 2) * |(|s| - |t|)| :=
    mul_le_mul_of_nonneg_left (le_abs_self _) hb0
  have hL := bandL_lip hb |s| |t|
  unfold bandU1
  rcases min_choice ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |t|) ((1 - (2 * Real.sqrt 2 - H)) + bandL H |t|)
    with h | h <;> rw [h]
  · exact le_trans (min_le_left _ _) (by linarith)
  · exact le_trans (min_le_right _ _) (by linarith)

lemma band_a_nonneg (hb : InBand H) : 0 ≤ 1 - (2 * Real.sqrt 2 - H) := by
  have := sqrt2_hi'; linarith [hb.1]

lemma bandU1_le (hb : InBand H) (m : ℝ) : bandU1 H m ≤ bandU0 H m + (1 - (2 * Real.sqrt 2 - H)) := by
  have ha := band_a_nonneg hb
  unfold bandU1 bandU0
  rw [← min_add_add_right]
  exact le_min (le_trans (min_le_left _ _) (by linarith)) (le_trans (min_le_right _ _) (by linarith))

lemma bandU0_le (hb : InBand H) (m : ℝ) : bandU0 H m ≤ bandU1 H m + (1 - (2 * Real.sqrt 2 - H)) := by
  have ha := band_a_nonneg hb
  unfold bandU1 bandU0
  rw [← min_add_add_right]
  exact le_min (le_trans (min_le_left _ _) (by linarith)) (le_trans (min_le_right _ _) (by linarith))

end lip

/-- Lines leaning the same way: a square costs at least a grid unit plus `β|m - m'|`. -/
lemma kap_same' {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hs : 0 ≤ m * m') :
    1 + (1 - H / 2) * |m - m'| ≤ kap C S H m m' := by
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  have pos : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → 0 ≤ m → 0 ≤ m' →
      1 + (1 - H / 2) * |m - m'| ≤ kap C S H m m' := by
    intro C S m m' hC hS hCS hm hm'
    rcases le_total m' m with h | h
    · rw [abs_of_nonneg (by linarith)]; exact kap_same hC hS hCS hm' h
    · rw [abs_of_nonpos (by linarith), kap_symm]
      have := kap_same (H := H) hC hS hCS hm h; linarith
  have neg : m ≤ 0 → m' ≤ 0 → 1 + (1 - H / 2) * |m - m'| ≤ kap C S H m m' := by
    intro hm hm'
    rw [kap_mirror, show m - m' = -m' - -m by ring, abs_sub_comm]
    have := pos (m := -m') (m' := -m) hS hC hSC (by linarith) (by linarith)
    rwa [abs_sub_comm] at this
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · exact pos hC hS hCS hm hm'
  · rcases mul_eq_zero.1 (le_antisymm (mul_nonpos_of_nonneg_of_nonpos hm hm') hs) with h | h
    · exact neg (by linarith) hm'
    · exact pos hC hS hCS hm (by linarith)
  · rcases mul_eq_zero.1 (le_antisymm (mul_nonpos_of_nonpos_of_nonneg hm hm') hs) with h | h
    · exact pos hC hS hCS (by linarith) hm'
    · exact neg hm (by linarith)
  · exact neg hm hm'

/-! ## The band potentials are admissible -/

section adm
variable {H : ℝ}

lemma band_step01 (hb : InBand H) {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (m m' : ℝ) :
    bandU1 H m' ≤ bandU0 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) := by
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · have k := kap_same' (H := H) hC hS hCS (mul_nonneg hm hm')
    have l := bandU1_lip hb m' m; have u := bandU1_le hb m
    rw [abs_sub_comm] at l; linarith
  · have := band_Sa (x := m) (y := -m') hC hS hCS hb hm (by linarith)
    rwa [neg_neg, bandU1_neg] at this
  · have := band_Sa (x := -m) (y := m') hS hC hSC hb (by linarith) hm'
    rwa [bandU0_neg, ← kap_symm, ← kap_mirror] at this
  · have k := kap_same' (H := H) hC hS hCS (mul_nonneg_of_nonpos_of_nonpos hm hm')
    have l := bandU1_lip hb m' m; have u := bandU1_le hb m
    rw [abs_sub_comm] at l; linarith

lemma band_step10 (hb : InBand H) {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (m m' : ℝ) :
    bandU0 H m' ≤ bandU1 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) := by
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · have k := kap_same' (H := H) hC hS hCS (mul_nonneg hm hm')
    have l := bandU0_lip hb m' m; have u := bandU0_le hb m
    rw [abs_sub_comm] at l; linarith
  · have := band_Sb (x := m) (y := -m') hC hS hCS hb hm (by linarith)
    rwa [neg_neg, bandU0_neg] at this
  · have := band_Sb (x := -m) (y := m') hS hC hSC hb (by linarith) hm'
    rwa [bandU1_neg, ← kap_symm, ← kap_mirror] at this
  · have k := kap_same' (H := H) hC hS hCS (mul_nonneg_of_nonpos_of_nonpos hm hm')
    have l := bandU0_lip hb m' m; have u := bandU0_le hb m
    rw [abs_sub_comm] at l; linarith

/-- A cheap step, `m ≥ 0 ≥ m'`. -/
lemma band_cheap_pos (hb : InBand H) {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) {m m' : ℝ}
    (hm : 0 ≤ m) (hm' : m' ≤ 0)
    (hch : kap C S H m m' - (2 * Real.sqrt 2 - H) < 1 - (2 * Real.sqrt 2 - H)) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU0 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + bandU0 H m' ∧
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU1 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + bandU1 H m' := by
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  have e : m' = -(-m') := by ring
  rcases le_total m (1 / 2) with hx | hx
  · have := band_N (x := m) (y := -m') hC hS hCS hb hm (by linarith) hx
    rw [← e] at this; linarith
  rcases le_total (-m') (1 / 2) with hy | hy
  · have := band_N (x := -m') (y := m) hS hC hSC hb (by linarith) hm hy
    rw [← kap_mirror] at this; linarith
  have c0 := band_C0 (x := m) (y := -m') hC hS hCS hb hm (by linarith) hx hy
  have c1 := band_C1 (x := m) (y := -m') hC hS hCS hb hm (by linarith) hx hy
  rw [← e, bandU0_neg] at c0; rw [← e, bandU1_neg] at c1
  exact ⟨c0, c1⟩

lemma band_cheap (hb : InBand H) {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (m m' : ℝ)
    (hch : kap C S H m m' - (2 * Real.sqrt 2 - H) < 1 - (2 * Real.sqrt 2 - H)) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU0 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + bandU0 H m' ∧
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU1 H m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + bandU1 H m' := by
  have hb0 : 0 ≤ 1 - H / 2 := by linarith [hb.2]
  have same : 0 ≤ m * m' → False := fun hs => by
    have k := kap_same' (H := H) hC hS hCS hs
    have : 0 ≤ (1 - H / 2) * |m - m'| := mul_nonneg hb0 (abs_nonneg _)
    linarith
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · exact (same (mul_nonneg hm hm')).elim
  · exact band_cheap_pos hb hC hS hCS hm hm' hch
  · rw [kap_symm] at hch ⊢
    have := band_cheap_pos hb hC hS hCS hm' hm hch
    constructor <;> linarith [this.1, this.2]
  · exact (same (mul_nonneg_of_nonpos_of_nonpos hm hm')).elim

theorem band_potOK (hb : InBand H) : PotOK H (bandU0 H) (bandU1 H) where
  wall := fun C S m hC hS hCS _ => by
    have := kap_wall (H := H) (m := m) hC hS hCS
    have : bandU1 H m ≤ (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m| := by unfold bandU1; exact min_le_left _ _
    linarith
  step10 := fun C S m m' hC hS hCS _ => band_step10 hb hC hS hCS m m'
  step01 := fun C S m m' hC hS hCS _ => band_step01 hb hC hS hCS m m'
  cheap := fun C S m m' hC hS hCS _ hch => band_cheap hb hC hS hCS m m' hch

end adm

/-! ## Odd `n` in a strip down to `H*(9)` -/

/-- **Odd `n`, `H*(9) ≤ H ≤ 1.87`**: `L ≥ min(n, F_n)`. -/
theorem strip_odd_band {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hb : InBand H) (hp : RectPacks n L H) :
    min (n : ℝ) (((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H) ≤ L := by
  have := strip_odd_pot (band_potOK hb) (by linarith [hb.2]) hn hp
  have e : (n : ℝ) * (2 * Real.sqrt 2 - H) + 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) =
      ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H := by ring
  rwa [e] at this

/-- **Odd `n`, `H*(9) ≤ H < 2`**: `L ≥ min(n, F_n)`, `F_n = (n - 2)(2√2 - H) + 4 - H`. -/
theorem strip_odd_lower' {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hH1 : (14 * Real.sqrt 2 - 5) / 8 ≤ H) (hH2 : H < 2)
    (hp : RectPacks n L H) : min (n : ℝ) (((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H) ≤ L := by
  rcases le_total H (187 / 100) with h | h
  · exact strip_odd_band hn ⟨hH1, h⟩ hp
  · exact strip_odd_lower hn (by norm_num at h ⊢; linarith) hH2 hp

/-- **Odd `n`, `H*(9) ≤ H < 2`, exactly**: `n` unit squares fit in `L × H` iff `L ≥ n` or `L ≥ F_n`. -/
theorem strip_odd_iff' {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hn3 : 3 ≤ n) (hH1 : (14 * Real.sqrt 2 - 5) / 8 ≤ H)
    (hH2 : H < 2) : RectPacks n L H ↔ (n : ℝ) ≤ L ∨ ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H ≤ L := by
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  constructor
  · intro hp; exact min_le_iff.1 (strip_odd_lower' hn hH1 hH2 hp)
  · rintro (h | h)
    · exact strip_row h (by linarith)
    · exact strip_zigzag hn hn3 (by linarith) hH2.le h

/-- **Odd `n ≤ 9`, every strip `1 ≤ H < 2`, exactly** (`n = 3, 5, 7, 9`): `n` unit squares fit in `L × H` iff
`L ≥ n` or `L ≥ F_n`.  Below `H*(9)` the strip lies inside the strip of height `H*(9)`, where `F_n - n = (9 - n) a ≥ 0`. -/
theorem strip_odd_iff_le9 {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hn3 : 3 ≤ n) (hn9 : n ≤ 9) (hH1 : 1 ≤ H)
    (hH2 : H < 2) : RectPacks n L H ↔ (n : ℝ) ≤ L ∨ ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H ≤ L := by
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hn9' : (n : ℝ) ≤ 9 := by exact_mod_cast hn9
  have hn3' : (3 : ℝ) ≤ n := by exact_mod_cast hn3
  rcases le_total ((14 * Real.sqrt 2 - 5) / 8) H with h | h
  · exact strip_odd_iff' hn hn3 h hH2
  constructor
  · intro hp
    have hp9 : RectPacks n L ((14 * Real.sqrt 2 - 5) / 8) := (RectPacks.mono hp.swap h).swap
    have := strip_odd_lower' hn le_rfl (by linarith) hp9
    left
    rcases min_le_iff.1 this with h1 | h1
    · exact h1
    · -- at H*(9): F_n - n = (9 - n)(3 - 2√2)/8 ≥ 0
      nlinarith
  · rintro (h1 | h1)
    · exact strip_row h1 hH1
    · -- below H*(9), F_n ≥ n, so the row fits
      exact strip_row (by nlinarith) hH1

end SquarePacking
