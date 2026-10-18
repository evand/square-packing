import Sqpack.RyuStair

/-!
# The staircase end: the parameters and the count

The construction of Ryu's Section 6 (S1)–(S6), with two simplifications that keep the proof short: any tilt `α` with
`m cos α + sin α ≤ g₀` is allowed (no rational grid), and the layers above the last column run up to
`⌊g(ξ_I − ḡ tan α)⌋`.  For an end region with ceiling `g(X) = y + (X − σ) μ`, `μ = σ/γ`:

* `m = ⌊ḡ⌋ + 1` (`ḡ = g(w)`), columns `i < I = ⌊(w − mτ)/ς⌋` (`τ = tan α`, `ς = sec α`) in the slabs
  `ξ_i = mτ + iς`, with `J_i = ⌊(g(iς) − H₀)/(1 − cos α)⌋` (`H₀ = m cos α + sin α`) and `m − J_i` squares;
* layers `j < J₁` of `⌊w − a_j⌋` squares, `a_j = ξ_{p(j)} − jτ`, `p(j) = #{i < I : J_i ≤ j}`;
* wedge rows `r < ⌊g₀⌋` of `⌊(m − r − 1)τ⌋` squares.

`EndP.packs`: these pack `endReg`.  `EndP.count`: the number of squares is at least `EndP.lb`, a closed form
(the floors lose at most `1` each).
-/

namespace UnitSquarePacking.Ryu

open Set UnitSquarePacking.Conjectures

/-- The data of a staircase end: the end region (`γ, σ, w, y`) and the column tilt (`ca, sa`). -/
structure EndP where
  (γ σ w y ca sa : ℝ)

namespace EndP

variable (P : EndP)

noncomputable def μ : ℝ := P.σ / P.γ
noncomputable def g (X : ℝ) : ℝ := P.y + (X - P.σ) * P.μ
noncomputable def g0 : ℝ := P.g 0
noncomputable def gb : ℝ := P.g P.w
noncomputable def m : ℕ := ⌊P.gb⌋₊ + 1
noncomputable def τ : ℝ := P.sa / P.ca
noncomputable def ς : ℝ := 1 / P.ca
noncomputable def H0 : ℝ := P.m * P.ca + P.sa
noncomputable def I : ℕ := ⌊(P.w - P.m * P.τ) / P.ς⌋₊
noncomputable def J (i : ℕ) : ℕ := ⌊(P.g (i * P.ς) - P.H0) / (1 - P.ca)⌋₊
noncomputable def xi (i : ℕ) : ℝ := P.m * P.τ + i * P.ς
noncomputable def J1 : ℕ := max ((Finset.range P.I).sup P.J) ⌊P.g (P.xi P.I - P.gb * P.τ)⌋₊
noncomputable def p (j : ℕ) : ℕ := ((Finset.range P.I).filter (fun i => P.J i ≤ j)).card
noncomputable def a (j : ℕ) : ℝ := P.xi (P.p j) - j * P.τ
noncomputable def L (j : ℕ) : ℕ := ⌊P.w - P.a j⌋₊
noncomputable def R0 : ℕ := ⌊P.g0⌋₊
noncomputable def l (r : ℕ) : ℕ := ⌊((P.m : ℝ) - r - 1) * P.τ⌋₊
/-- The number of squares. -/
noncomputable def N : ℕ :=
  (∑ i ∈ Finset.range P.I, (P.m - P.J i)) + (∑ j ∈ Finset.range P.J1, P.L j) + ∑ r ∈ Finset.range P.R0, P.l r

/-- The hypotheses: the band side (`γ > 0`, `σ ≥ 0`), the tilt (`0 < ca < 1`, `sa > 0`, `ca² + sa² = 1`,
`m ca + sa ≤ g₀`) and `mτ ≤ w`. -/
structure Hyp : Prop where
  hγ : 0 < P.γ
  hσ : 0 ≤ P.σ
  hca : 0 < P.ca
  hca1 : P.ca < 1
  hsa : 0 < P.sa
  hu : P.ca ^ 2 + P.sa ^ 2 = 1
  htilt : P.H0 ≤ P.g0
  hw : P.m * P.τ ≤ P.w

variable {P}

section facts
variable (h : P.Hyp)
include h

lemma μ_nonneg : 0 ≤ P.μ := div_nonneg h.hσ h.hγ.le
lemma τ_pos : 0 < P.τ := div_pos h.hsa h.hca
lemma ς_pos : 0 < P.ς := div_pos one_pos h.hca
lemma ca_τ : P.ca * P.τ = P.sa := by unfold τ; rw [mul_div_assoc', mul_div_cancel_left₀ _ h.hca.ne']
lemma ca_ς : P.ca * P.ς = 1 := by unfold ς; rw [mul_one_div, div_self h.hca.ne']
lemma one_le_ς : 1 ≤ P.ς := by unfold ς; rw [le_div_iff₀ h.hca]; linarith [h.hca1]
lemma sa_le_τ : P.sa ≤ P.τ := by
  unfold τ; rw [le_div_iff₀ h.hca]; nlinarith [h.hca1, h.hsa]

lemma g_mono {X X' : ℝ} (hX : X ≤ X') : P.g X ≤ P.g X' := by
  unfold g; nlinarith [μ_nonneg h]

lemma g0_le_gb (hw0 : 0 ≤ P.w) : P.g0 ≤ P.gb := g_mono h hw0

lemma gb_lt_m : P.gb < P.m := by
  unfold m; push_cast; exact Nat.lt_floor_add_one _

lemma xi_I_le : P.xi P.I ≤ P.w := by
  unfold xi
  have hP : 0 ≤ P.w - P.m * P.τ := by linarith [h.hw]
  have := Nat.floor_le (div_nonneg hP (ς_pos h).le)
  have e : ((P.w - P.m * P.τ) / P.ς) * P.ς = P.w - P.m * P.τ := by field_simp [(ς_pos h).ne']
  have := mul_le_mul_of_nonneg_right this (ς_pos h).le
  unfold I; nlinarith

lemma w_nonneg : 0 ≤ P.w := le_trans (mul_nonneg (by positivity) (τ_pos h).le) h.hw

lemma xi_le {i : ℕ} (hi : i ≤ P.I) : P.xi i ≤ P.xi P.I := by
  unfold xi
  have : (i : ℝ) ≤ P.I := by exact_mod_cast hi
  nlinarith [ς_pos h]

/-- `H₀ + J_i (1 − ca) ≤ g(iς)`: the column top is under the ceiling. -/
lemma J_spec (i : ℕ) (hi : i ≤ P.I) : P.H0 + P.J i * (1 - P.ca) ≤ P.g (i * P.ς) := by
  have hc : 0 < 1 - P.ca := by linarith [h.hca1]
  have hg : P.g0 ≤ P.g (i * P.ς) := g_mono h (by positivity [ς_pos h])
  have h0 : 0 ≤ (P.g (i * P.ς) - P.H0) / (1 - P.ca) := div_nonneg (by linarith [h.htilt]) hc.le
  have := Nat.floor_le h0
  unfold J
  have e := mul_le_mul_of_nonneg_right this hc.le
  rw [div_mul_cancel₀ _ hc.ne'] at e
  linarith

lemma J_mono {i i' : ℕ} (hii : i ≤ i') : P.J i ≤ P.J i' := by
  unfold J
  apply Nat.floor_mono
  apply div_le_div_of_nonneg_right _ (by linarith [h.hca1])
  have : (i : ℝ) * P.ς ≤ i' * P.ς := mul_le_mul_of_nonneg_right (by exact_mod_cast hii) (ς_pos h).le
  linarith [g_mono h this]

lemma iς_le_w {i : ℕ} (hi : i ≤ P.I) : (i : ℝ) * P.ς ≤ P.w := by
  have := xi_le h hi
  have := xi_I_le h
  unfold xi at *
  nlinarith [τ_pos h]

/-- `J_i < m`. -/
lemma J_lt_m (i : ℕ) (hi : i ≤ P.I) : P.J i < P.m := by
  have hc : 0 < 1 - P.ca := by linarith [h.hca1]
  have hgi : P.g (i * P.ς) ≤ P.gb := g_mono h (iς_le_w h hi)
  have hm := gb_lt_m h
  unfold J
  rw [Nat.floor_lt' (by have := Nat.succ_pos ⌊P.gb⌋₊; unfold m; omega), div_lt_iff₀ hc]
  unfold H0
  nlinarith [h.hsa]

lemma p_le (j : ℕ) : P.p j ≤ P.I := by
  unfold p
  exact (Finset.card_filter_le _ _).trans (by simp)

/-- If `J_i ≤ j` then `i < p(j)`. -/
lemma lt_p {i j : ℕ} (hi : i < P.I) (hJ : P.J i ≤ j) : i < P.p j := by
  unfold p
  have : Finset.range (i + 1) ⊆ (Finset.range P.I).filter (fun i' => P.J i' ≤ j) := by
    intro i' hi'
    simp only [Finset.mem_range] at hi'
    simp only [Finset.mem_filter, Finset.mem_range]
    exact ⟨by omega, (J_mono h (by omega)).trans hJ⟩
  have := Finset.card_le_card this
  simp only [Finset.card_range] at this
  omega

/-- If `p(j) ≤ i < I` then `j < J_i`. -/
lemma lt_J {i j : ℕ} (hi : i < P.I) (hp : P.p j ≤ i) : j < P.J i := by
  by_contra hc
  have := lt_p h hi (by omega : P.J i ≤ j)
  omega

/-- If `j < J_i` then `p(j) ≤ i`. -/
lemma p_le_of {i j : ℕ} (hi : i < P.I) (hj : j < P.J i) : P.p j ≤ i := by
  unfold p
  have : (Finset.range P.I).filter (fun i' => P.J i' ≤ j) ⊆ Finset.range i := by
    intro i' hi'
    simp only [Finset.mem_filter, Finset.mem_range] at hi' ⊢
    by_contra hc
    have := J_mono h (by omega : i ≤ i')
    omega
  have := Finset.card_le_card this
  simpa using this

lemma J_le_J1 {i : ℕ} (hi : i < P.I) : P.J i ≤ P.J1 :=
  le_max_of_le_left (Finset.le_sup (Finset.mem_range.2 hi))

lemma H0_nonneg : 0 ≤ P.H0 := add_nonneg (mul_nonneg (Nat.cast_nonneg _) h.hca.le) h.hsa.le
lemma g0_nonneg : 0 ≤ P.g0 := le_trans (H0_nonneg h) h.htilt
lemma gb_nonneg : 0 ≤ P.gb := le_trans (g0_nonneg h) (g0_le_gb h (w_nonneg h))

lemma J1_lt_m : P.J1 < P.m := by
  have hm : P.m = ⌊P.gb⌋₊ + 1 := rfl
  unfold J1
  apply max_lt
  · have : (Finset.range P.I).sup P.J ≤ P.m - 1 := Finset.sup_le fun i hi => by
      have := J_lt_m h i (by simp at hi; omega); omega
    omega
  · have hgb : P.g (P.xi P.I - P.gb * P.τ) ≤ P.gb := by
      apply g_mono h
      have : 0 ≤ P.gb * P.τ := mul_nonneg (gb_nonneg h) (τ_pos h).le
      linarith [xi_I_le h]
    have := Nat.floor_le_floor hgb
    omega

lemma xi_mono {i i' : ℕ} (hii : i ≤ i') : P.xi i ≤ P.xi i' := by
  unfold xi
  have : (i : ℝ) ≤ i' := by exact_mod_cast hii
  nlinarith [ς_pos h]

lemma a_le_w (j : ℕ) : P.a j ≤ P.w := by
  unfold a
  have := (xi_mono h (p_le h j)).trans (xi_I_le h)
  have : 0 ≤ (j : ℝ) * P.τ := mul_nonneg (Nat.cast_nonneg _) (τ_pos h).le
  linarith

lemma a_ge (j : ℕ) : P.m * P.τ - j * P.τ ≤ P.a j := by
  unfold a xi
  have : 0 ≤ ((P.p j : ℕ) : ℝ) * P.ς := mul_nonneg (Nat.cast_nonneg _) (ς_pos h).le
  linarith

lemma L_le (j : ℕ) : (P.L j : ℝ) ≤ P.w - P.a j := Nat.floor_le (by linarith [a_le_w h j])

/-- `H(J) = H₀ + J(1 − ca) ≥ J` for `J ≤ m`. -/
lemma H_ge {Jc : ℕ} (hJ : Jc ≤ P.m) : (Jc : ℝ) ≤ P.H0 + Jc * (1 - P.ca) := by
  unfold H0
  have : (Jc : ℝ) ≤ P.m := by exact_mod_cast hJ
  nlinarith [h.hca.le, h.hsa.le]

/-- **The layers fit**: `j + 1 ≤ g(a_j)` for `j < J₁`. -/
lemma layer_fit {j : ℕ} (hj : j < P.J1) : (j : ℝ) + 1 ≤ P.g (P.a j) := by
  have hjm : j < P.m := hj.trans (J1_lt_m h)
  have hjm' : (j : ℝ) ≤ P.m := by exact_mod_cast hjm.le
  rcases Nat.lt_or_ge (P.p j) P.I with hp | hp
  · -- a column `q = p(j)` with `J_q > j`
    set q := P.p j with hq
    have hJq : j < P.J q := lt_J h hp le_rfl
    have h1 : P.xi q - P.m * P.τ ≤ P.a j := by
      unfold a; rw [← hq]; nlinarith [τ_pos h]
    have h2 : P.xi q - P.m * P.τ = q * P.ς := by unfold xi; ring
    have h3 := J_spec h q hp.le
    have h4 := H_ge h (J_lt_m h q hp.le).le
    have h5 : (j : ℝ) + 1 ≤ P.J q := by exact_mod_cast hJq
    have := g_mono h h1
    rw [h2] at this
    linarith
  · -- right of all columns: `j < ⌊g(ξ_I − ḡτ)⌋`
    have hpI : P.p j = P.I := le_antisymm (p_le h j) hp
    have hjB : j < ⌊P.g (P.xi P.I - P.gb * P.τ)⌋₊ := by
      rcases lt_max_iff.1 hj with hA | hB
      · obtain ⟨i, hi, hji⟩ := Finset.lt_sup_iff.1 hA
        have := p_le_of h (Finset.mem_range.1 hi) hji
        have := Finset.mem_range.1 hi
        omega
      · exact hB
    have hjg : (j : ℝ) + 1 ≤ ⌊P.g (P.xi P.I - P.gb * P.τ)⌋₊ := by exact_mod_cast hjB
    have hx1 : 1 ≤ P.g (P.xi P.I - P.gb * P.τ) := Nat.floor_pos.1 (by omega)
    have hfl := Nat.floor_le (show 0 ≤ P.g (P.xi P.I - P.gb * P.τ) by linarith)
    have hjgb : (j : ℝ) ≤ P.gb := by
      have := Nat.floor_le_floor (g_mono h (show P.xi P.I - P.gb * P.τ ≤ P.w by
        have : 0 ≤ P.gb * P.τ := mul_nonneg (gb_nonneg h) (τ_pos h).le
        linarith [xi_I_le h]))
      have : (j : ℝ) + 1 ≤ ⌊P.gb⌋₊ := le_trans hjg (by exact_mod_cast this)
      linarith [Nat.floor_le (gb_nonneg h)]
    have h1 : P.xi P.I - P.gb * P.τ ≤ P.a j := by
      unfold a; rw [hpI]; nlinarith [τ_pos h]
    linarith [g_mono h h1]

/-- The ceiling, `γ`-scaled: `γ(Y) − σX ≤ γy − σ²` iff `Y ≤ g(X)`. -/
lemma ceil_iff (X Y : ℝ) (hY : Y ≤ P.g X) : P.γ * Y - P.σ * X ≤ P.γ * P.y - P.σ ^ 2 := by
  have hγ := h.hγ
  have e : P.γ * P.g X = P.γ * P.y + P.σ * X - P.σ ^ 2 := by
    unfold g μ; field_simp; ring
  nlinarith [mul_le_mul_of_nonneg_left hY hγ.le]

/-- **The staircase packs the end region.** -/
theorem packs : PacksIn P.N (endReg P.γ P.σ P.w P.y) := by
  obtain ⟨α, hc, hs⟩ := exists_angle h.hu
  have hca := h.hca
  have hsa := h.hsa
  have hτ := τ_pos h
  have hς := ς_pos h
  have cτ := ca_τ h
  have cς := ca_ς h
  have hς1 := one_le_ς h
  have hsaτ := sa_le_τ h
  have castn : ∀ i ≤ P.I, (((P.m - P.J i : ℕ)) : ℝ) = P.m - P.J i := fun i hi => by
    rw [Nat.cast_sub (J_lt_m h i hi).le]
  have hU : ∀ i : ℕ, P.ca * (P.xi i - P.J i * P.τ) + P.sa * P.J i = P.m * P.sa + i := by
    intro i; unfold xi; linear_combination (P.m - P.J i : ℝ) * cτ + (i : ℝ) * cς
  have hUa : ∀ j : ℕ, P.ca * P.a j + P.sa * j = P.m * P.sa + P.p j := by
    intro j; unfold a xi; linear_combination (P.m - j : ℝ) * cτ + (P.p j : ℝ) * cς
  have hgb := gb_lt_m h
  have hR0 : ∀ r < P.R0, (r : ℝ) + 1 ≤ P.g0 := fun r hr => by
    have : (r : ℝ) + 1 ≤ P.R0 := by exact_mod_cast hr
    exact this.trans (Nat.floor_le (g0_nonneg h))
  have hg0m : P.g0 < P.m := lt_of_le_of_lt (g0_le_gb h (w_nonneg h)) hgb
  have hl : ∀ r < P.R0, (P.l r : ℝ) ≤ ((P.m : ℝ) - r - 1) * P.τ := fun r hr => by
    have := hR0 r hr
    exact Nat.floor_le (mul_nonneg (by linarith) hτ.le)
  have hp := staircase (γ := P.γ) (σ := P.σ) (w := P.w) (y := P.y) h.hγ.le h.hσ (α := α) (by rw [hc]; exact hca.le)
    (by rw [hs]; exact hsa.le) P.I P.J1 P.R0 (fun i => P.xi i - P.J i * P.τ) P.J (fun i => P.m - P.J i) P.L P.l
    (fun i hi => by
      rw [hs, castn i hi.le]; unfold xi
      have := J_lt_m h i hi.le
      have : (P.J i : ℝ) ≤ P.m := by exact_mod_cast this.le
      nlinarith [mul_le_mul_of_nonneg_left hsaτ (by linarith : (0 : ℝ) ≤ P.m - P.J i),
        mul_nonneg (Nat.cast_nonneg i) hς.le, mul_nonneg (Nat.cast_nonneg (P.J i)) hτ.le])
    (fun i hi => by
      rw [hc]
      have := (xi_mono h (show i + 1 ≤ P.I by omega)).trans (xi_I_le h)
      unfold xi at this ⊢; push_cast at this
      nlinarith [mul_nonneg (Nat.cast_nonneg (P.J i)) hτ.le, h.hca1])
    (fun i hi => by
      rw [hc, hs, castn i hi.le]
      have e1 := J_spec h i hi.le
      have e2 : (i : ℝ) * P.ς ≤ P.xi i - P.J i * P.τ - (P.m - P.J i) * P.sa := by
        unfold xi
        have := J_lt_m h i hi.le
        have : (P.J i : ℝ) ≤ P.m := by exact_mod_cast this.le
        nlinarith [mul_le_mul_of_nonneg_left hsaτ (by linarith : (0 : ℝ) ≤ P.m - P.J i)]
      have e3 : P.J i + P.sa + (P.m - P.J i) * P.ca = P.H0 + P.J i * (1 - P.ca) := by unfold H0; ring
      rw [e3]
      have := ceil_iff h (i * P.ς) _ e1
      nlinarith [mul_le_mul_of_nonneg_left e2 h.hσ])
    (fun i hi i' hi' hii => by
      rw [hc, hs, hU, hU]
      have : (i : ℝ) + 1 ≤ i' := by exact_mod_cast hii
      linarith)
    (fun j hj => by linarith [L_le h j, a_ge h j, show (j : ℝ) ≤ P.m by exact_mod_cast (hj.trans (J1_lt_m h)).le,
      mul_le_mul_of_nonneg_right (show (j : ℝ) ≤ P.m by exact_mod_cast (hj.trans (J1_lt_m h)).le) hτ.le])
    (fun j hj => by
      have h1 := layer_fit h hj
      have h2 : P.g (P.a j) ≤ P.g (P.w - P.L j) := g_mono h (by linarith [L_le h j])
      exact ceil_iff h _ _ (by push_cast; linarith))
    (fun r hr => by
      have := hl r hr
      have : ((P.m : ℝ) - r - 1) * P.τ ≤ P.m * P.τ := by nlinarith [(Nat.cast_nonneg r : (0 : ℝ) ≤ r)]
      linarith [h.hw])
    (fun r hr => by
      have := hR0 r hr
      have e : P.γ * P.g0 = P.γ * P.y - P.σ ^ 2 := by unfold g0 g μ; field_simp [h.hγ.ne']; ring
      nlinarith [mul_le_mul_of_nonneg_left this h.hγ.le])
    (fun i hi j hj => by
      rw [hc, hs]
      by_cases hJ : P.J i ≤ j
      · left
        have h1 := lt_p h hi hJ
        have h1' : (i : ℝ) + 1 ≤ P.p j := by exact_mod_cast h1
        rw [hU]
        have := mul_le_mul_of_nonneg_left (show P.a j ≤ P.w - P.L j by linarith [L_le h j]) hca.le
        linarith [hUa j]
      · right; push_cast; exact_mod_cast (by omega : j + 1 ≤ P.J i))
    (fun i hi r hr => by
      rw [hc, hs, hU]
      have := mul_le_mul_of_nonneg_left (hl r hr) hca.le
      have e : P.ca * (((P.m : ℝ) - r - 1) * P.τ) = ((P.m : ℝ) - r - 1) * P.sa := by
        rw [← cτ]; ring
      nlinarith [(Nat.cast_nonneg i : (0 : ℝ) ≤ i)])
    (fun j hj r hr => by
      rw [hc, hs]
      have := mul_le_mul_of_nonneg_left (hl r hr) hca.le
      have e : P.ca * (((P.m : ℝ) - r - 1) * P.τ) = ((P.m : ℝ) - r - 1) * P.sa := by
        rw [← cτ]; ring
      have := mul_le_mul_of_nonneg_left (show P.a j ≤ P.w - P.L j by linarith [L_le h j]) hca.le
      nlinarith [hUa j, (Nat.cast_nonneg (P.p j) : (0 : ℝ) ≤ P.p j)])
  exact hp

lemma sum_p : ∑ j ∈ Finset.range P.J1, (P.p j : ℝ) = ∑ i ∈ Finset.range P.I, ((P.J1 : ℝ) - P.J i) := by
  have hnat : ∑ j ∈ Finset.range P.J1, P.p j = ∑ i ∈ Finset.range P.I, (P.J1 - P.J i) := by
    unfold p
    simp only [Finset.card_filter]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i hi => ?_
    rw [← Finset.card_filter]
    have : (Finset.range P.J1).filter (fun j => P.J i ≤ j) = Finset.Ico (P.J i) P.J1 := by
      ext j; simp only [Finset.mem_filter, Finset.mem_range, Finset.mem_Ico]; omega
    rw [this, Nat.card_Ico]
  have := congrArg (fun n : ℕ => (n : ℝ)) hnat
  simp only [Nat.cast_sum] at this
  rw [this]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Nat.cast_sub (J_le_J1 h (Finset.mem_range.1 hi))]

end facts

lemma sum_range_cast (n : ℕ) : ∑ i ∈ Finset.range n, (i : ℝ) = n * (n - 1) / 2 := by
  induction n with
  | zero => simp
  | succ n ih => rw [Finset.sum_range_succ, ih]; push_cast; ring

/-- The closed form of the count. -/
noncomputable def lb (P : EndP) : ℝ :=
  P.I * P.m + P.J1 * (P.w - P.m * P.τ - P.I * P.ς - 1) + P.τ * P.J1 * (P.J1 - 1) / 2 +
    (P.ς - 1) * ((P.I * (P.g0 - P.H0) + P.ς * P.μ * (P.I * (P.I - 1) / 2)) / (1 - P.ca) - P.I) +
    (P.R0 * (P.m - 1) * P.τ - P.τ * P.R0 * (P.R0 - 1) / 2 - P.R0)

section count
variable (h : P.Hyp)
include h

/-- **The count**: `N ≥ lb`. -/
theorem count : P.lb ≤ P.N := by
  have hτ := τ_pos h
  have hς := ς_pos h
  have hς1 := one_le_ς h
  have hc : 0 < 1 - P.ca := by linarith [h.hca1]
  -- the three sums, as reals
  have e1 : ((∑ i ∈ Finset.range P.I, (P.m - P.J i) : ℕ) : ℝ) = ∑ i ∈ Finset.range P.I, ((P.m : ℝ) - P.J i) := by
    push_cast
    refine Finset.sum_congr rfl fun i hi => ?_
    rw [Nat.cast_sub (J_lt_m h i (Finset.mem_range.1 hi).le).le]
  have e2 : ∑ j ∈ Finset.range P.J1, (P.w - P.a j - 1) ≤ ((∑ j ∈ Finset.range P.J1, P.L j : ℕ) : ℝ) := by
    push_cast
    exact Finset.sum_le_sum fun j _ => by have := Nat.lt_floor_add_one (P.w - P.a j); unfold L; linarith
  have e3 : ∑ r ∈ Finset.range P.R0, (((P.m : ℝ) - r - 1) * P.τ - 1) ≤
      ((∑ r ∈ Finset.range P.R0, P.l r : ℕ) : ℝ) := by
    push_cast
    exact Finset.sum_le_sum fun r _ => by
      have := Nat.lt_floor_add_one (((P.m : ℝ) - r - 1) * P.τ); unfold l; linarith
  -- the layer sum
  have e4 : ∑ j ∈ Finset.range P.J1, (P.w - P.a j - 1) =
      P.J1 * (P.w - P.m * P.τ - 1) - P.ς * ∑ j ∈ Finset.range P.J1, (P.p j : ℝ) + P.τ * (P.J1 * (P.J1 - 1) / 2) := by
    unfold a xi
    simp only [Finset.sum_sub_distrib, Finset.sum_add_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul,
      ← Finset.sum_mul, ← Finset.mul_sum, sum_range_cast]
    ring
  rw [sum_p h] at e4
  -- the column sum and `Σ J_i`
  have e5 : ∑ i ∈ Finset.range P.I, ((P.m : ℝ) - P.J i) = P.I * P.m - ∑ i ∈ Finset.range P.I, (P.J i : ℝ) := by
    rw [Finset.sum_sub_distrib]; simp
  have e6 : ∑ i ∈ Finset.range P.I, ((P.J1 : ℝ) - P.J i) = P.I * P.J1 - ∑ i ∈ Finset.range P.I, (P.J i : ℝ) := by
    rw [Finset.sum_sub_distrib]; simp
  have hJ : ∑ i ∈ Finset.range P.I, ((P.g0 - P.H0 + i * P.ς * P.μ) / (1 - P.ca) - 1) ≤
      ∑ i ∈ Finset.range P.I, (P.J i : ℝ) := by
    refine Finset.sum_le_sum fun i hi => ?_
    have := Nat.lt_floor_add_one ((P.g (i * P.ς) - P.H0) / (1 - P.ca))
    have eg : P.g (i * P.ς) = P.g0 + i * P.ς * P.μ := by unfold g0 g; ring
    unfold J; rw [eg] at this ⊢
    rw [show P.g0 + ↑i * P.ς * P.μ - P.H0 = P.g0 - P.H0 + ↑i * P.ς * P.μ by ring] at this ⊢
    linarith
  have eJ : ∑ i ∈ Finset.range P.I, ((P.g0 - P.H0 + i * P.ς * P.μ) / (1 - P.ca) - 1) =
      (P.I * (P.g0 - P.H0) + P.ς * P.μ * (P.I * (P.I - 1) / 2)) / (1 - P.ca) - P.I := by
    simp only [Finset.sum_sub_distrib, ← Finset.sum_div, Finset.sum_add_distrib, Finset.sum_const, Finset.card_range,
      nsmul_eq_mul, ← Finset.sum_mul, sum_range_cast]
    ring
  rw [eJ] at hJ
  have e7 : ∑ r ∈ Finset.range P.R0, (((P.m : ℝ) - r - 1) * P.τ - 1) =
      P.R0 * (P.m - 1) * P.τ - P.τ * P.R0 * (P.R0 - 1) / 2 - P.R0 := by
    simp only [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul, ← Finset.sum_mul,
      sum_range_cast]
    ring
  rw [e7] at e3
  unfold N lb
  push_cast at e1 e2 e3 ⊢
  rw [e1] at *
  have := mul_le_mul_of_nonneg_left hJ (by linarith : (0 : ℝ) ≤ P.ς - 1)
  nlinarith

end count


end EndP

end UnitSquarePacking.Ryu
