import Sqpack.StripWalk
import Sqpack.StripPot
import Sqpack.StripS

/-!
# Odd `n` in a strip: the walk argument, given the per-square inequalities `(S)` and `(C')`

Strip height `H`, `λ = 2√2 - H`, `a = 1 - λ`, `β = 1 - H/2`, `g = a + β`.  For odd `n` the target is
`L ≥ min n (n λ + 2g)` (`= min(n, F_n)`, `F_n = (n - 2) λ + 4 - H`, the 45° zigzag).

Per-square inequalities (step cost `κ = kap C S H m m'`, excess `c = κ - λ`):
* `kap_wall` (proved): `κ(0, m) ≥ 1 + β|m|`;
* `kap_B` (proved): `κ ≥ 2(√2 - 1) + β(|m| + |m'|)`;
* `KapS H` (hypothesis here): `κ ≥ min (1 + β||m| - |m'||) (Z - β(|m| + |m'|))`, `Z = 2 + 2√2 - 2H`;
* `KapC H` (hypothesis here): `κ ≥ 2λ - 1 + β(|m'| - |m|)`.

Walk argument (`strip_odd_of`): prefix potentials by line-index parity, `u₁ x = a + β|x|` (odd),
`u₀ x = min (2a + β|x|) (a + β(2 - |x|))` (even); `prefix_bound` by induction, `suffix_bound` backwards (`kap_symm`).
If no step is cheap (`c ≥ a`) the sum is `≥ n a`; a cheap step is never at a wall (`kap_wall`), so it has
interior lines on both sides, of the same parity (n odd), and prefix + step + suffix `≥ 2g`.
-/

open Finset

namespace SquarePacking

/-- `(S)`: a square costs a grid unit or the zigzag amount. -/
def KapS (H : ℝ) : Prop := ∀ C S m m' : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
  min (1 + (1 - H / 2) * |(|m| - |m'|)|) (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) ≤
    kap C S H m m'

/-- `(C')`. -/
def KapC (H : ℝ) : Prop := ∀ C S m m' : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
  2 * (2 * Real.sqrt 2 - H) - 1 + (1 - H / 2) * (|m'| - |m|) ≤ kap C S H m m'

section
variable (H : ℝ)

/-- The parity potentials. -/
noncomputable def uOdd (x : ℝ) : ℝ := (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |x|
noncomputable def uEven (x : ℝ) : ℝ :=
  min (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |x|) ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - |x|))
noncomputable def uPot (j : ℕ) (x : ℝ) : ℝ := if j % 2 = 1 then uOdd H x else uEven H x
end

/-- One step of the potentials: `u_j(m) + c(m, m') ≥ u_{j+1}(m')` (either parity). -/
lemma pot_step {H C S m m' : ℝ} (j : ℕ) (hS : KapS H) (hC : KapC H) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2)
    (hC0 : 0 ≤ C) (hS0 : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hfit : C + S ≤ H) :
    uPot H (j + 1) m' ≤ uPot H j m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) := by
  have hb : 0 ≤ 1 - H / 2 := by linarith
  have B := kap_B (m := m) (m' := m') hC0 hS0 hCS hH1 hH2 hfit
  have Cc := hC C S m m' hC0 hS0 hCS hfit
  have Sx := hS C S m m' hC0 hS0 hCS hfit
  unfold uPot
  rcases Nat.mod_two_eq_zero_or_one j with hj | hj
  · -- even → odd: u1(m') ≤ u0(m) + c  (both routes of u0)
    have : (j + 1) % 2 = 1 := by omega
    simp only [this, hj, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false, uOdd, uEven]
    rcases min_choice (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|)
      ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - |m|)) with h | h <;> rw [h] <;> nlinarith
  · -- odd → even: u0(m') ≤ u1(m) + c  via (S)
    have : (j + 1) % 2 = 0 := by omega
    simp only [this, hj, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false, uOdd, uEven]
    have hab : (1 - H / 2) * |(|m| - |m'|)| ≥ (1 - H / 2) * (|m'| - |m|) :=
      mul_le_mul_of_nonneg_left (by rw [abs_sub_comm]; exact le_abs_self _) hb
    rcases min_choice (1 + (1 - H / 2) * |(|m| - |m'|)|)
      (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) with h | h <;> rw [h] at Sx
    · exact (min_le_left _ _).trans (by linarith)
    · exact (min_le_right _ _).trans (by linarith)

/-- A cheap step (`c < a`) between lines of the same parity class pays the two zigzag ends:
`u_j(m) + c + u_j(m') ≥ 2g`. -/
lemma pot_cheap {H C S m m' : ℝ} (j : ℕ) (hS : KapS H) (hC : KapC H) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H)
    (hH2 : H ≤ 2) (hC0 : 0 ≤ C) (hS0 : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hfit : C + S ≤ H)
    (hcheap : kap C S H m m' - (2 * Real.sqrt 2 - H) < 1 - (2 * Real.sqrt 2 - H)) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
      uPot H j m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + uPot H j m' := by
  have hb : 0 ≤ 1 - H / 2 := by linarith
  have ha : 0 ≤ 1 - (2 * Real.sqrt 2 - H) := by linarith
  have B := kap_B (m := m) (m' := m') hC0 hS0 hCS hH1 hH2 hfit
  have C1 := hC C S m m' hC0 hS0 hCS hfit
  have C2 := hC C S m' m hC0 hS0 hCS hfit
  rw [kap_symm] at C2
  -- cheap ⇒ the zigzag branch of (S)
  have Z : 2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|) ≤ kap C S H m m' := by
    have Sx := hS C S m m' hC0 hS0 hCS hfit
    rcases min_choice (1 + (1 - H / 2) * |(|m| - |m'|)|)
      (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) with h | h <;> rw [h] at Sx
    · nlinarith [abs_nonneg (|m| - |m'|)]
    · exact Sx
  unfold uPot
  rcases Nat.mod_two_eq_zero_or_one j with hj | hj
  · simp only [hj, show (0 : ℕ) ≠ 1 by norm_num, if_false, uEven]
    rcases min_choice (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|)
      ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - |m|)) with h | h <;> rw [h] <;>
    rcases min_choice (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m'|)
      ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - |m'|)) with h' | h' <;> rw [h'] <;> linarith
  · simp only [hj, if_true, uOdd]; linarith

/-- **Prefixes**: the excess up to line `j` (`1 ≤ j < n`) is at least the potential of line `j`. -/
lemma prefix_bound {H : ℝ} {n : ℕ} {P Q m : ℕ → ℝ} (hS : KapS H) (hC : KapC H)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2)
    (hfr : ∀ i < n, 0 ≤ P i ∧ 0 ≤ Q i ∧ P i ^ 2 + Q i ^ 2 = 1 ∧ P i + Q i ≤ H) (hm0 : m 0 = 0) :
    ∀ j, 1 ≤ j → j < n → uPot H j (m j) ≤
      ∑ l ∈ range j, (kap (P l) (Q l) H (m l) (m (l + 1)) - (2 * Real.sqrt 2 - H)) := by
  intro j hj1 hjn
  induction j, hj1 using Nat.le_induction with
  | base =>
    obtain ⟨h1, h2, h3, -⟩ := hfr 0 (by omega)
    have := kap_wall (H := H) (m := m 1) h1 h2 h3
    simp only [range_one, sum_singleton, zero_add, hm0, uPot, show (1 : ℕ) % 2 = 1 from rfl, if_true, uOdd]
    linarith
  | succ j hj ih =>
    obtain ⟨h1, h2, h3, h4⟩ := hfr j (by omega)
    rw [sum_range_succ]
    have := pot_step (m := m j) (m' := m (j + 1)) j hS hC hH1 hH2 h1 h2 h3 h4
    linarith [ih (by omega)]

/-- **Suffixes**: the excess from line `j` (`1 ≤ j < n`) to the right wall is at least the potential of class
`n - j`. -/
lemma suffix_bound {H : ℝ} {n : ℕ} {P Q m : ℕ → ℝ} (hS : KapS H) (hC : KapC H)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H ≤ 2)
    (hfr : ∀ i < n, 0 ≤ P i ∧ 0 ≤ Q i ∧ P i ^ 2 + Q i ^ 2 = 1 ∧ P i + Q i ≤ H) (hmn : m n = 0) :
    ∀ d, 1 ≤ d → ∀ j, j + d = n → 1 ≤ j → uPot H d (m j) ≤
      ∑ l ∈ Ico j n, (kap (P l) (Q l) H (m l) (m (l + 1)) - (2 * Real.sqrt 2 - H)) := by
  intro d hd
  induction d, hd using Nat.le_induction with
  | base =>
    intro j hj _
    obtain ⟨h1, h2, h3, -⟩ := hfr j (by omega)
    have := kap_wall (H := H) (m := m j) h1 h2 h3
    have hjn : j + 1 = n := hj
    rw [show Ico j n = {j} by rw [← hjn]; exact Nat.Ico_succ_singleton j, sum_singleton, hjn, hmn, kap_symm]
    simp only [uPot, show (1 : ℕ) % 2 = 1 from rfl, if_true, uOdd]
    linarith
  | succ d hd ih =>
    intro j hj hj1
    obtain ⟨h1, h2, h3, h4⟩ := hfr j (by omega)
    rw [← Finset.insert_Ico_add_one_left_eq_Ico (by omega : j < n), sum_insert (by simp)]
    have := pot_step (m := m (j + 1)) (m' := m j) d hS hC hH1 hH2 h1 h2 h3 h4
    rw [kap_symm] at this
    linarith [ih (j + 1) (by omega) (by omega)]

/-- **Odd `n` in a strip, given `(S)` and `(C')`**: `L ≥ min(n, n λ + 2g)`. -/
theorem strip_odd_of {n : ℕ} {L H : ℝ} (hS : KapS H) (hC : KapC H) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H < 2)
    (hn : n % 2 = 1) (hp : RectPacks n L H) :
    min (n : ℝ) (n * (2 * Real.sqrt 2 - H) + 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2))) ≤ L := by
  obtain ⟨P, Q, m, hfr, hm0, hmn, hsum⟩ := strip_walk (by omega) hp hH2
  set lam := 2 * Real.sqrt 2 - H
  set a := 1 - lam
  set c : ℕ → ℝ := fun l => kap (P l) (Q l) H (m l) (m (l + 1)) - lam with hc
  have hsplit : ∑ l ∈ range n, kap (P l) (Q l) H (m l) (m (l + 1)) = n * lam + ∑ l ∈ range n, c l := by
    simp only [c, sum_sub_distrib, sum_const, card_range, nsmul_eq_mul]; ring
  have ha : 0 ≤ a := by simp only [a, lam]; linarith
  -- walls are never cheap
  have wall0 : a ≤ c 0 := by
    obtain ⟨h1, h2, h3, -⟩ := hfr 0 (by omega)
    have := kap_wall (H := H) (m := m 1) h1 h2 h3
    simp only [c, a, hm0]; nlinarith [abs_nonneg (m 1)]
  have wallN : a ≤ c (n - 1) := by
    obtain ⟨h1, h2, h3, -⟩ := hfr (n - 1) (by omega)
    have := kap_wall (H := H) (m := m (n - 1)) h1 h2 h3
    simp only [c, a, show n - 1 + 1 = n by omega, hmn, kap_symm _ _ _ _ (0 : ℝ)]; nlinarith [abs_nonneg (m (n - 1))]
  by_cases hall : ∀ l < n, a ≤ c l
  · have : n * a ≤ ∑ l ∈ range n, c l := by
      have := sum_le_sum (fun l (hl : l ∈ range n) => hall l (mem_range.1 hl))
      simpa [mul_comm] using this
    rw [min_le_iff]; left
    have : (n : ℝ) = n * lam + n * a := by simp only [a]; ring
    linarith
  push Not at hall
  obtain ⟨i, hin, hci⟩ := hall
  have hi0 : i ≠ 0 := fun h => by subst h; linarith
  have hiN : i ≠ n - 1 := fun h => by subst h; linarith
  have hpre := prefix_bound hS hC hH1 hH2.le hfr hm0 i (by omega) (by omega)
  have hsuf := suffix_bound hS hC hH1 hH2.le hfr hmn (n - (i + 1)) (by omega) (i + 1) (by omega) (by omega)
  rw [show n - (i + 1) = (n - (i + 1)) by rfl] at hsuf
  have hpar : uPot H (n - (i + 1)) = uPot H i := by
    unfold uPot; ext x
    have : (n - (i + 1)) % 2 = i % 2 := by omega
    rw [this]
  rw [hpar] at hsuf
  obtain ⟨h1, h2, h3, h4⟩ := hfr i hin
  have hch := pot_cheap (m := m i) (m' := m (i + 1)) i hS hC hH1 hH2.le h1 h2 h3 h4 hci
  have hdecomp : ∑ l ∈ range n, c l = ∑ l ∈ range i, c l + c i + ∑ l ∈ Ico (i + 1) n, c l := by
    rw [← sum_range_succ, sum_range_add_sum_Ico _ (by omega : i + 1 ≤ n)]
  rw [min_le_iff]; right
  simp only [c] at hdecomp hpre hsuf
  linarith

/-- `(S)` holds on the whole range `2√2 - 1 ≤ H < 2` (`kap_S`). -/
theorem kapS_of {H : ℝ} (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H < 2) : KapS H :=
  fun _ _ _ _ hC hS hCS hfit => kap_S hC hS hCS hH1 hH2 hfit

/-- **Odd `n` in a strip, given `(C')`** (`(S)` is proved): `L ≥ min(n, F_n)`. -/
theorem strip_odd_of_C {n : ℕ} {L H : ℝ} (hC : KapC H) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H < 2)
    (hn : n % 2 = 1) (hp : RectPacks n L H) :
    min (n : ℝ) (n * (2 * Real.sqrt 2 - H) + 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2))) ≤ L :=
  strip_odd_of (kapS_of hH1 hH2) hC hH1 hH2 hn hp

end SquarePacking
