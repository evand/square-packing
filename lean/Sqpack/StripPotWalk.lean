import Sqpack.StripWalk
import Sqpack.StripPot

/-!
# Odd `n` in a strip: the walk argument for any parity potentials

`strip_odd_pot`: let `u₀, u₁ : ℝ → ℝ` be potentials on even / odd lines (`uP u₀ u₁ j` picks by the parity of `j`).
Per square, with frame `(C, S)` and excess `c = kap - λ` (`λ = 2√2 - H`, `a = 1 - λ`, `g = a + (1 - H/2)`):
* (wall) `u₁ m ≤ c(0, m)`;
* (step) `u₁ m' ≤ u₀ m + c(m, m')` and `u₀ m' ≤ u₁ m + c(m, m')`;
* (cheap) `c(m, m') < a ⇒ 2g ≤ u_P m + c + u_P m'` for both parities `P`.
Then for odd `n`, `L ≥ min(n, nλ + 2g) = min(n, F_n)`.  Same proof as `strip_odd_of` (`StripOdd`): prefix and
suffix bounds by induction; a cheap step is never at a wall (`kap_wall`), and for odd `n` the prefix and suffix around
it end on lines of the same parity.
-/

open Finset

namespace SquarePacking

/-- The potential of line `j`: `u₁` on odd lines, `u₀` on even lines. -/
def uP (u₀ u₁ : ℝ → ℝ) (j : ℕ) : ℝ → ℝ := if j % 2 = 1 then u₁ else u₀

/-- The per-square hypotheses on a pair of parity potentials. -/
structure PotOK (H : ℝ) (u₀ u₁ : ℝ → ℝ) : Prop where
  wall : ∀ C S m : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
    u₁ m ≤ kap C S H 0 m - (2 * Real.sqrt 2 - H)
  step10 : ∀ C S m m' : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
    u₀ m' ≤ u₁ m + (kap C S H m m' - (2 * Real.sqrt 2 - H))
  step01 : ∀ C S m m' : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
    u₁ m' ≤ u₀ m + (kap C S H m m' - (2 * Real.sqrt 2 - H))
  cheap : ∀ C S m m' : ℝ, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H →
    kap C S H m m' - (2 * Real.sqrt 2 - H) < 1 - (2 * Real.sqrt 2 - H) →
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ u₀ m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + u₀ m' ∧
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ u₁ m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + u₁ m'

lemma PotOK.step {H : ℝ} {u₀ u₁ : ℝ → ℝ} (hu : PotOK H u₀ u₁) (j : ℕ) {C S m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S)
    (hCS : C ^ 2 + S ^ 2 = 1) (hfit : C + S ≤ H) :
    uP u₀ u₁ (j + 1) m' ≤ uP u₀ u₁ j m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) := by
  unfold uP
  rcases Nat.mod_two_eq_zero_or_one j with hj | hj
  · have : (j + 1) % 2 = 1 := by omega
    simp only [this, hj, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false]
    exact hu.step01 C S m m' hC hS hCS hfit
  · have : (j + 1) % 2 = 0 := by omega
    simp only [this, hj, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false]
    exact hu.step10 C S m m' hC hS hCS hfit

lemma PotOK.cheap' {H : ℝ} {u₀ u₁ : ℝ → ℝ} (hu : PotOK H u₀ u₁) (j : ℕ) {C S m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S)
    (hCS : C ^ 2 + S ^ 2 = 1) (hfit : C + S ≤ H)
    (hch : kap C S H m m' - (2 * Real.sqrt 2 - H) < 1 - (2 * Real.sqrt 2 - H)) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
      uP u₀ u₁ j m + (kap C S H m m' - (2 * Real.sqrt 2 - H)) + uP u₀ u₁ j m' := by
  have := hu.cheap C S m m' hC hS hCS hfit hch
  unfold uP
  split_ifs
  · exact this.2
  · exact this.1

/-- **Prefixes**: the excess up to line `j` (`1 ≤ j < n`) is at least the potential of line `j`. -/
lemma prefix_bound_pot {H : ℝ} {u₀ u₁ : ℝ → ℝ} (hu : PotOK H u₀ u₁) {n : ℕ} {P Q m : ℕ → ℝ}
    (hfr : ∀ i < n, 0 ≤ P i ∧ 0 ≤ Q i ∧ P i ^ 2 + Q i ^ 2 = 1 ∧ P i + Q i ≤ H) (hm0 : m 0 = 0) :
    ∀ j, 1 ≤ j → j < n → uP u₀ u₁ j (m j) ≤
      ∑ l ∈ range j, (kap (P l) (Q l) H (m l) (m (l + 1)) - (2 * Real.sqrt 2 - H)) := by
  intro j hj1 hjn
  induction j, hj1 using Nat.le_induction with
  | base =>
    obtain ⟨h1, h2, h3, h4⟩ := hfr 0 (by omega)
    have := hu.wall (P 0) (Q 0) (m 1) h1 h2 h3 h4
    simp only [range_one, sum_singleton, zero_add, hm0, uP, show (1 : ℕ) % 2 = 1 from rfl, if_true]
    exact this
  | succ j hj ih =>
    obtain ⟨h1, h2, h3, h4⟩ := hfr j (by omega)
    rw [sum_range_succ]
    have := hu.step (m := m j) (m' := m (j + 1)) j h1 h2 h3 h4
    linarith [ih (by omega)]

/-- **Suffixes**: the excess from line `j` to the right wall (`d = n - j ≥ 1` squares) is at least `uP d (m j)`. -/
lemma suffix_bound_pot {H : ℝ} {u₀ u₁ : ℝ → ℝ} (hu : PotOK H u₀ u₁) {n : ℕ} {P Q m : ℕ → ℝ}
    (hfr : ∀ i < n, 0 ≤ P i ∧ 0 ≤ Q i ∧ P i ^ 2 + Q i ^ 2 = 1 ∧ P i + Q i ≤ H) (hmn : m n = 0) :
    ∀ d, 1 ≤ d → ∀ j, j + d = n → 1 ≤ j → uP u₀ u₁ d (m j) ≤
      ∑ l ∈ Ico j n, (kap (P l) (Q l) H (m l) (m (l + 1)) - (2 * Real.sqrt 2 - H)) := by
  intro d hd
  induction d, hd using Nat.le_induction with
  | base =>
    intro j hj _
    obtain ⟨h1, h2, h3, h4⟩ := hfr j (by omega)
    have := hu.wall (P j) (Q j) (m j) h1 h2 h3 h4
    have hjn : j + 1 = n := hj
    rw [show Ico j n = {j} by rw [← hjn]; exact Nat.Ico_succ_singleton j, sum_singleton, hjn, hmn, kap_symm]
    simp only [uP, show (1 : ℕ) % 2 = 1 from rfl, if_true]
    exact this
  | succ d hd ih =>
    intro j hj hj1
    obtain ⟨h1, h2, h3, h4⟩ := hfr j (by omega)
    rw [← Finset.insert_Ico_add_one_left_eq_Ico (by omega : j < n), sum_insert (by simp)]
    have := hu.step (m := m (j + 1)) (m' := m j) d h1 h2 h3 h4
    rw [kap_symm] at this
    linarith [ih (j + 1) (by omega) (by omega)]

/-- **Odd `n` in a strip, for any admissible parity potentials**: `L ≥ min(n, n λ + 2g)`. -/
theorem strip_odd_pot {n : ℕ} {L H : ℝ} {u₀ u₁ : ℝ → ℝ} (hu : PotOK H u₀ u₁) (hH2 : H < 2)
    (hn : n % 2 = 1) (hp : RectPacks n L H) :
    min (n : ℝ) (n * (2 * Real.sqrt 2 - H) + 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2))) ≤ L := by
  obtain ⟨P, Q, m, hfr, hm0, hmn, hsum⟩ := strip_walk (by omega) hp hH2
  set lam := 2 * Real.sqrt 2 - H
  set a := 1 - lam
  set c : ℕ → ℝ := fun l => kap (P l) (Q l) H (m l) (m (l + 1)) - lam with hc
  have hsplit : ∑ l ∈ range n, kap (P l) (Q l) H (m l) (m (l + 1)) = n * lam + ∑ l ∈ range n, c l := by
    simp only [c, sum_sub_distrib, sum_const, card_range, nsmul_eq_mul]; ring
  have hb : 0 ≤ 1 - H / 2 := by linarith
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
  have hpre := prefix_bound_pot hu hfr hm0 i (by omega) (by omega)
  have hsuf := suffix_bound_pot hu hfr hmn (n - (i + 1)) (by omega) (i + 1) (by omega) (by omega)
  have hpar : uP u₀ u₁ (n - (i + 1)) = uP u₀ u₁ i := by
    unfold uP
    have : (n - (i + 1)) % 2 = i % 2 := by omega
    rw [this]
  rw [hpar] at hsuf
  obtain ⟨h1, h2, h3, h4⟩ := hfr i hin
  have hch := hu.cheap' (m := m i) (m' := m (i + 1)) i h1 h2 h3 h4 hci
  have hdecomp : ∑ l ∈ range n, c l = ∑ l ∈ range i, c l + c i + ∑ l ∈ Ico (i + 1) n, c l := by
    rw [← sum_range_succ, sum_range_add_sum_Ico _ (by omega : i + 1 ≤ n)]
  rw [min_le_iff]; right
  simp only [c] at hdecomp hpre hsuf
  linarith

end SquarePacking
