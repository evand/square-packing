import Sqpack.Rect5Mid

/-!
# Five unit squares in a strip of height `< 2`

(work in progress)  `rect5_strip_lower_of`: the chain argument for five squares, reduced to the real
inequality `Final5` (still to be proved).  The chain: sort `a, b, c, d, e` by centre `x`; left shelf on `b`
(from `a`), right shelf on `d` (from `e`); lines `b | c` and `c | d`, normalised to unit directions
`(Cᵢ, Sᵢ)` with a sign; the one-square lemma on `b` and on `d` (`osl_set*`), the middle square `c` (`mid_set`).
Target (Friedman's pieces 4.1, 4.2 of `S₄`): `L ≥ 5` or `L + 4H ≥ 4 + 6√2`.
-/

namespace SquarePacking

/-- The frame aligned with a line of direction `(C, S)` and slope sign `ε` (`true`: normal `(S, C)`). -/
def alFr (ε : Bool) (C S : ℝ) : ℝ × ℝ := if ε then (S, C) else (C, S)

/-- The signed `y`-coefficient of a line's normal. -/
def sgnC (ε : Bool) (C : ℝ) : ℝ := if ε then C else -C

/-- One-square bound for `X`, in direction `(C, S)`. -/
def OSB (C S H X : ℝ) : Prop := 2 * S - C * (H - 1) ≤ X ∨ C * (S - C) + (1 + C + S - H) ≤ C * X

/-- **The remaining real inequality** for five squares in a strip.  Line 1 (`b | c`) and line 2 (`c | d`), with
`b`'s bound on `s₁` (shifted by `C₁ H` when `ε₁`), `d`'s bound on `S₂ L - s₂` (shifted by `C₂ H` when `¬ε₂`),
and the middle square's width at one of its four candidate frames. -/
def Final5 : Prop :=
  ∀ (ε₁ ε₂ : Bool) (C₁ S₁ C₂ S₂ H L s₁ s₂ : ℝ) (x : ℝ × ℝ), 0 ≤ C₁ → 0 < S₁ → C₁ ^ 2 + S₁ ^ 2 = 1 →
    0 ≤ C₂ → 0 < S₂ → C₂ ^ 2 + S₂ ^ 2 = 1 → 1 ≤ H → H < 2 →
    OSB C₁ S₁ H (if ε₁ then s₁ - C₁ * H else s₁) →
    OSB C₂ S₂ H (if ε₂ then S₂ * L - s₂ else S₂ * L - C₂ * H - s₂) →
    (x = (1, 0) ∨ x = alFr ε₁ C₁ S₁ ∨ x = alFr ε₂ C₂ S₂ ∨ x = (0, 1)) →
    midW (sgnC ε₂ C₂ / S₂ - sgnC ε₁ C₁ / S₁) (1 / (2 * S₁)) (1 / (2 * S₂)) (alFr ε₁ C₁ S₁) (alFr ε₂ C₂ S₂) x
      (if 0 ≤ sgnC ε₂ C₂ / S₂ - sgnC ε₁ C₁ / S₁ then (x.1 + x.2) / 2 else H - (x.1 + x.2) / 2) ≤
      s₂ / S₂ - s₁ / S₁ →
    5 ≤ L ∨ 4 + 6 * Real.sqrt 2 ≤ L + 4 * H

/-- The aligned frame of `alFr`. -/
lemma alFr_aligned {ε : Bool} {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) :
    Aligned S (sgnC ε C) (alFr ε C S) := by
  cases ε
  · simpa [alFr, sgnC] using aligned_neg hC hS hCS
  · simpa [alFr, sgnC] using aligned_pos hC hS hCS

/-- A separating line `u₁ x + u₂ y = s` (`u₁ > 0`) as a unit direction with a sign. -/
lemma unit_line {u₁ u₂ : ℝ} (hu : 0 < u₁) : ∃ (ε : Bool) (C S r : ℝ), 0 ≤ C ∧ 0 < S ∧ C ^ 2 + S ^ 2 = 1 ∧
    0 < r ∧ ∀ p : ℝ × ℝ, S * p.1 + sgnC ε C * p.2 = (u₁ * p.1 + u₂ * p.2) / r := by
  set r := Real.sqrt (u₁ ^ 2 + u₂ ^ 2)
  have hr0 : 0 < r := Real.sqrt_pos.2 (by positivity)
  have hr2 : r ^ 2 = u₁ ^ 2 + u₂ ^ 2 := Real.sq_sqrt (by positivity)
  have hCS : (|u₂| / r) ^ 2 + (u₁ / r) ^ 2 = 1 := by rw [div_pow, div_pow, sq_abs]; field_simp; linarith
  rcases le_total u₂ 0 with h | h
  · refine ⟨false, |u₂| / r, u₁ / r, r, by positivity, div_pos hu hr0, hCS, hr0, fun p => ?_⟩
    rw [show sgnC false (|u₂| / r) = -(|u₂| / r) from rfl, abs_of_nonpos h]; field_simp
  · refine ⟨true, |u₂| / r, u₁ / r, r, by positivity, div_pos hu hr0, hCS, hr0, fun p => ?_⟩
    rw [show sgnC true (|u₂| / r) = |u₂| / r from rfl, abs_of_nonneg h]; field_simp

/-- **Five squares in a strip**, given `Final5`. -/
theorem rect5_strip_lower_of (hF : Final5) {L H : ℝ} (hp : RectPacks 5 L H) (hH : H < 2) :
    5 ≤ L ∨ 4 + 6 * Real.sqrt 2 ≤ L + 4 * H := by
  have hH1 : 1 ≤ H := (one_le_sides (hp.choose_spec.choose_spec.1 0)).2
  obtain ⟨ctr, ang, hin, hd⟩ := hp
  set σ := Tuple.sort (fun i => (ctr i).1)
  have hm := Tuple.monotone_sort (fun i => (ctr i).1)
  have hne : ∀ i j : Fin 5, i ≠ j → σ i ≠ σ j := fun i j h => σ.injective.ne h
  have h01 : (ctr (σ 0)).1 ≤ (ctr (σ 1)).1 := hm (by decide)
  have h12 : (ctr (σ 1)).1 ≤ (ctr (σ 2)).1 := hm (by decide)
  have h23 : (ctr (σ 2)).1 ≤ (ctr (σ 3)).1 := hm (by decide)
  have h34 : (ctr (σ 3)).1 ≤ (ctr (σ 4)).1 := hm (by decide)
  obtain ⟨lo, l0, l1, hl⟩ := left_shelf (hin (σ 0)) (hin (σ 1)) (hd _ _ (hne 0 1 (by decide))) h01 hH
  obtain ⟨lo2, r0, r1, hr⟩ :=
    right_shelf (hin (σ 3)) (hin (σ 4)) (hd _ _ (hne 3 4 (by decide))) h34 hH
  obtain ⟨u₁, u₂, s, hu, hb, hc⟩ :=
    separate_lr (hin (σ 1)) (hin (σ 2)) (hd _ _ (hne 1 2 (by decide))) h12 hH
  obtain ⟨w₁, w₂, s', hw, hc', hdd⟩ :=
    separate_lr (hin (σ 2)) (hin (σ 3)) (hd _ _ (hne 2 3 (by decide))) h23 hH
  obtain ⟨ε₁, C₁, S₁, r, hC₁, hS₁, hCS₁, hrp, e₁⟩ := unit_line (u₂ := u₂) hu
  obtain ⟨ε₂, C₂, S₂, r', hC₂, hS₂, hCS₂, hrp', e₂⟩ := unit_line (u₂ := w₂) hw
  have scale : ∀ {a b ρ : ℝ}, 0 < ρ → a ≤ b → a / ρ ≤ b / ρ := fun hρ h => div_le_div_of_nonneg_right h hρ.le
  -- line 1: `b` left of it, `c` right of it; line 2: `c` left, `d` right
  have B : ∀ p ∈ sq (ctr (σ 1)) (ang (σ 1)) 1, S₁ * p.1 + sgnC ε₁ C₁ * p.2 ≤ s / r :=
    fun p hp => by rw [e₁]; exact scale hrp (hb p hp)
  have C1 : ∀ p ∈ sq (ctr (σ 2)) (ang (σ 2)) 1, s / r ≤ S₁ * p.1 + sgnC ε₁ C₁ * p.2 :=
    fun p hp => by rw [e₁]; exact scale hrp (hc p hp)
  have C2 : ∀ p ∈ sq (ctr (σ 2)) (ang (σ 2)) 1, S₂ * p.1 + sgnC ε₂ C₂ * p.2 ≤ s' / r' :=
    fun p hp => by rw [e₂]; exact scale hrp' (hc' p hp)
  have D : ∀ p ∈ sq (ctr (σ 3)) (ang (σ 3)) 1, s' / r' ≤ S₂ * p.1 + sgnC ε₂ C₂ * p.2 :=
    fun p hp => by rw [e₂]; exact scale hrp' (hdd p hp)
  -- the one-square bounds
  have ob : OSB C₁ S₁ H (if ε₁ then s / r - C₁ * H else s / r) := by
    cases ε₁
    · exact osl_set (hin (σ 1)) hH l0 l1 hl hC₁ hS₁ hCS₁ (fun p hp => by simpa [sgnC] using B p hp)
    · exact osl_set_up (hin (σ 1)) hH l0 l1 hl hC₁ hS₁ hCS₁ (fun p hp => by simpa [sgnC] using B p hp)
  have od : OSB C₂ S₂ H (if ε₂ then S₂ * L - s' / r' else S₂ * L - C₂ * H - s' / r') := by
    cases ε₂
    · exact osl_set_rot (hin (σ 3)) hH r0 r1 hr hC₂ hS₂ hCS₂ (fun p hp => by simpa [sgnC] using D p hp)
    · exact osl_set_mir (hin (σ 3)) hH r0 r1 hr hC₂ hS₂ hCS₂ (fun p hp => by simpa [sgnC] using D p hp)
  obtain ⟨x, hx, hmid⟩ := mid_set (hin (σ 2)) hS₁ hS₂ (alFr_aligned hC₁ hS₁.le hCS₁)
    (alFr_aligned hC₂ hS₂.le hCS₂) C1 C2
  exact hF ε₁ ε₂ C₁ S₁ C₂ S₂ H L (s / r) (s' / r') x hC₁ hS₁ hCS₁ hC₂ hS₂ hCS₂ hH1 hH ob od hx hmid

end SquarePacking
