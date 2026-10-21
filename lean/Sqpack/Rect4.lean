import Sqpack.RectChain
import Sqpack.Rect4Core

/-!
# Four unit squares in a rectangle

(work in progress)  The one-square lemma `one_square`: a square `b` with a left shelf in a strip of height
`H < 2`, and a unit direction `(C, S)` (`S > 0`, `C ≥ 0`, normal `(S, -C)`).  Its support
`T = max_{p ∈ b} (S p₁ - C p₂)` satisfies `T ≥ 2S - C(H - 1)` (an upright square at the top beside the
shelf) or `C T ≥ C(S - C) + 1 + C + S - H` (the square at tilt `t` on the ceiling against the shelf's corner,
`(C, S) = (cos t, sin t)`: Friedman's `R₄` piece).  By witness case: `V` (side vertex) and `U` (edge above
it) give the upright bound; `D` (edge below it) gives one of the two, through `core_D1` / `core_D2`.
-/

namespace SquarePacking

section one
variable {P Q x y lo H C S : ℝ}

/-- Case `U`, `σ ≤ 0`. -/
lemma os_U_neg (hQ : 0 ≤ Q) (hC : 0 ≤ C) (hS : 0 < S) (hlo' : lo ≤ H - 1)
    (hQ1 : Q ≤ 1) (l1 : y + (P - Q) / 2 < lo)
    (l3 : 0 ≤ Q * (x - 1 - (P + Q) / 2) + P * (lo - (y + (P - Q) / 2)))
    (hs : S * P + -C * Q ≤ 0)
    (rep : 2 * (S + C) ≤ (P + Q) * (S + C) + (-(S * P + -C * Q)) + (S * Q + C * P)) :
    0 ≤ Q * (S * x - C * y + (-(S * P + -C * Q) + (S * Q + C * P)) / 2 - (2 * S - C * (H - 1))) := by
  linarith [mul_nonneg hS.le l3, mul_nonneg (mul_nonneg hQ hC) (by linarith : 0 ≤ H - 1 - lo),
    mul_nonneg (by linarith : 0 ≤ lo - (y + (P - Q) / 2)) (by linarith : 0 ≤ -(S * P + -C * Q)),
    mul_nonneg hQ (by linarith : 0 ≤ (P + Q) * (S + C) + (-(S * P + -C * Q)) + (S * Q + C * P) - 2 * (S + C)),
    mul_nonneg (mul_nonneg hQ hC) (by linarith : 0 ≤ 1 - Q)]

/-- Case `U`, `σ ≥ 0`. -/
lemma os_U_pos (hP : 0 ≤ P) (hQ : 0 ≤ Q) (hC : 0 ≤ C) (hS : 0 < S) (hPQ : P ^ 2 + Q ^ 2 = 1)
    (hy1 : (P + Q) / 2 ≤ y) (hlo' : lo ≤ H - 1) (hH : H < 2) (hP1 : P ≤ 1) (hQ1 : Q ≤ 1)
    (l3 : 0 ≤ Q * (x - 1 - (P + Q) / 2) + P * (lo - (y + (P - Q) / 2)))
    (hs : 0 ≤ S * P + -C * Q) :
    0 ≤ Q * (S * x - C * y + ((S * P + -C * Q) + (S * Q + C * P)) / 2 - (2 * S - C * (H - 1))) := by
  have e1 : S * P ^ 2 + S * Q ^ 2 = S := by linear_combination S * hPQ
  linarith [mul_nonneg hS.le l3, mul_nonneg hs (by linarith : 0 ≤ y - (P + Q) / 2),
    mul_nonneg (mul_nonneg hS.le hP) (by linarith : 0 ≤ H - 1 - lo),
    mul_nonneg hs (by linarith : 0 ≤ 2 - H),
    mul_nonneg (by linarith : 0 ≤ 1 - Q)
      (add_nonneg (mul_nonneg hC hQ) (mul_nonneg hS.le (by linarith : 0 ≤ 1 - P)))]

/-- Case `D`: `P T ≥ P (S - C) + P σ⁺ + γ (1 + P + Q - H)`, for `s = |σ|`. -/
lemma os_D_bound (hQ : 0 ≤ Q) (hC : 0 ≤ C) (hS : 0 < S) (hP : 0 ≤ P) (hlo : 0 ≤ lo)
    (hy2 : y ≤ H - (P + Q) / 2)
    (l3 : 0 ≤ P * (x - 1 - (P + Q) / 2) + Q * (y + (P - Q) / 2 - lo - 1)) (s : ℝ)
    (hs : s = S * P + -C * Q ∨ s = -(S * P + -C * Q)) :
    P * (S - C) + P * ((s + (S * P + -C * Q)) / 2) + (S * Q + C * P) * (1 + P + Q - H) ≤
      P * (S * x - C * y + (s + (S * Q + C * P)) / 2) := by
  have hγ0 : 0 ≤ S * Q + C * P := add_nonneg (mul_nonneg hS.le hQ) (mul_nonneg hC hP)
  rcases hs with rfl | rfl <;>
  linarith [mul_nonneg hS.le l3, mul_nonneg (mul_nonneg hS.le hQ) hlo,
    mul_nonneg hγ0 (by linarith : 0 ≤ H - (P + Q) / 2 - y)]

/-- Case `D`, `σ ≥ 0` (`core_D1`). -/
lemma os_D1 (hP0 : 0 < P) (hQ : 0 ≤ Q) (hC : 0 ≤ C) (hS : 0 < S) (hPQ : P ^ 2 + Q ^ 2 = 1)
    (hCS : C ^ 2 + S ^ 2 = 1) (hs : 0 ≤ S * P + -C * Q) {T : ℝ}
    (b : P * (S - C) + P * (((S * P + -C * Q) + (S * P + -C * Q)) / 2) +
      (S * Q + C * P) * (1 + P + Q - H) ≤ P * T) :
    2 * S - C * (H - 1) ≤ T ∨ C * (S - C) + (1 + C + S - H) ≤ C * T := by
  have core := core_D1 hP0.le hQ hC hS.le hPQ hCS (by linarith)
  have idD1 : S * (P * (S - C) + P * (S * P + -C * Q) + (S * Q + C * P) * (1 + P + Q - H)) =
      Q * (C * (S - C) + (1 + C + S - H)) + (S * P + -C * Q) * (2 * S - C * (H - 1)) +
        polyD1 P Q C S := by
    simp only [polyD1]; linear_combination S * (C + S) * hPQ + (-Q * (H - 2)) * hCS
  have hSP : 0 < S * P := mul_pos hS hP0
  have b' := mul_le_mul_of_nonneg_left b hS.le
  rcases le_total (C * (S - C) + (1 + C + S - H)) (C * (2 * S - C * (H - 1))) with hψ | hψ
  · right
    have h : 0 ≤ S * P * (C * T - (C * (S - C) + (1 + C + S - H))) := by
      nlinarith [mul_nonneg hs (sub_nonneg.2 hψ), mul_nonneg hC core]
    have := (mul_nonneg_iff_of_pos_left hSP).1 h; linarith
  · left
    have h : 0 ≤ S * P * (T - (2 * S - C * (H - 1))) := by
      nlinarith [mul_le_mul_of_nonneg_left hψ hQ]
    have := (mul_nonneg_iff_of_pos_left hSP).1 h; linarith

/-- Case `D`, `σ ≤ 0` (`core_D2`). -/
lemma os_D2 (hP0 : 0 < P) (hQ : 0 ≤ Q) (hC : 0 ≤ C) (hS : 0 < S) (hPQ : P ^ 2 + Q ^ 2 = 1)
    (hCS : C ^ 2 + S ^ 2 = 1) (hH : H < 2) (hs : S * P + -C * Q ≤ 0) {T : ℝ}
    (b : P * (S - C) + P * ((-(S * P + -C * Q) + (S * P + -C * Q)) / 2) +
      (S * Q + C * P) * (1 + P + Q - H) ≤ P * T) :
    C * (S - C) + (1 + C + S - H) ≤ C * T := by
  have core := core_D2 hP0.le hQ hC hS.le hPQ hCS (by linarith)
  have idD2 : C * (P * (S - C) + (S * Q + C * P) * (1 + P + Q - H)) - P * (C * (S - C) + (1 + C + S - H)) =
      polyD2 P Q C S + (2 - H) * (-S * (S * P + -C * Q)) := by
    simp only [polyD2]; linear_combination (-P * (H - 2)) * hCS
  have h : 0 ≤ P * (C * T - (C * (S - C) + (1 + C + S - H))) := by
    nlinarith [mul_le_mul_of_nonneg_left b hC, mul_nonneg (by linarith : 0 ≤ 2 - H)
      (mul_nonneg hS.le (by linarith : 0 ≤ -(S * P + -C * Q)))]
  have := (mul_nonneg_iff_of_pos_left hP0).1 h; linarith

end one

/-- The support of the framed square `(P, Q)`, centre `(x, y)`, in the direction `(S, -C)`. -/
noncomputable def suppT (C S P Q x y : ℝ) : ℝ :=
  S * x - C * y + (|S * P + -C * Q| + |-S * Q + -C * P|) / 2

/-- **One square.** -/
theorem one_square {P Q x y lo H C S : ℝ} (hP : 0 ≤ P) (hQ : 0 ≤ Q) (hPQ : P ^ 2 + Q ^ 2 = 1)
    (hy1 : (P + Q) / 2 ≤ y) (hy2 : y ≤ H - (P + Q) / 2) (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hH : H < 2)
    (hw : Wit P Q y lo (x - 1 - (P + Q) / 2)) (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1) :
    2 * S - C * (H - 1) ≤ suppT C S P Q x y ∨
      C * (S - C) + (1 + C + S - H) ≤ C * suppT C S P Q x y := by
  unfold suppT
  have hγ : |-S * Q + -C * P| = S * Q + C * P := by
    rw [abs_of_nonpos (by nlinarith [mul_nonneg hS.le hQ, mul_nonneg hC hP])]; ring
  have rep := upright_repl P Q S (-C) hP hQ hPQ hS.le
  rw [hγ, abs_neg, abs_of_nonneg hC] at rep
  rw [hγ]
  have hk : 1 ≤ P + Q := by nlinarith
  have hP1 : P ≤ 1 := by nlinarith
  have hQ1 : Q ≤ 1 := by nlinarith
  rcases hw with hA | ⟨l1, l2, l3⟩ | ⟨l1, l2, l3⟩
  · -- V: the side vertex clears the wall square
    left
    linarith [mul_le_mul_of_nonneg_left hA hS.le, mul_le_mul_of_nonneg_left hy2 hC]
  · -- U: the edge above the side vertex meets the band's bottom
    left
    have hQ0 : 0 < Q := by linarith
    rcases le_total (S * P + -C * Q) 0 with hs | hs
    · rw [abs_of_nonpos hs] at rep ⊢
      have := (mul_nonneg_iff_of_pos_left hQ0).1 (os_U_neg hQ hC hS hlo' hQ1 l1 l3 hs rep)
      linarith
    · rw [abs_of_nonneg hs]
      have := (mul_nonneg_iff_of_pos_left hQ0).1 (os_U_pos hP hQ hC hS hPQ hy1 hlo' hH hP1 hQ1 l3 hs)
      linarith
  · -- D: the edge below the side vertex meets the band's top
    have hP0 : 0 < P := by linarith
    rcases le_total 0 (S * P + -C * Q) with hs | hs
    · rw [abs_of_nonneg hs]
      exact os_D1 hP0 hQ hC hS hPQ hCS hs (os_D_bound hQ hC hS hP hlo hy2 l3 _ (Or.inl rfl))
    · rw [abs_of_nonpos hs]
      exact Or.inr (os_D2 hP0 hQ hC hS hPQ hCS hH hs (os_D_bound hQ hC hS hP hlo hy2 l3 _ (Or.inr rfl)))

/-! ## Squares against a shelf, in the four orientations -/

/-- **One square, set form.**  A square with a left shelf in the half-plane `S x - C y ≤ s`. -/
theorem osl_set {c : ℝ × ℝ} {θ L H lo C S s : ℝ} (hin : sq c θ 1 ⊆ rbox L H) (hH : H < 2)
    (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hl : ∀ p ∈ sq c θ 1, lo ≤ p.2 → p.2 ≤ lo + 1 → 1 ≤ p.1)
    (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hs : ∀ p ∈ sq c θ 1, S * p.1 + -C * p.2 ≤ s) :
    2 * S - C * (H - 1) ≤ s ∨ C * (S - C) + (1 + C + S - H) ≤ C * s := by
  obtain ⟨P, Q, hP, hQ, hPQ, -, -, hy1, hy2, hpt⟩ := framed hin
  have hk : 1 ≤ P + Q := by nlinarith
  have hw := witness (x := c.1) (by linarith) (by linarith)
    (fun u v hu hv h1 h2 => hl _ (hpt u v hu hv) h1 h2)
  have hT : suppT C S P Q c.1 c.2 ≤ s := by
    have := support_le hpt hs; unfold suppT; linarith
  rcases one_square hP hQ hPQ hy1 hy2 hlo hlo' hH hw hC hS hCS with h | h
  · left; linarith
  · right; nlinarith [mul_le_mul_of_nonneg_left hT hC]

/-- A square with a left shelf in the half-plane `S x + C y ≤ s` (vertical mirror). -/
theorem osl_set_up {c : ℝ × ℝ} {θ L H lo C S s : ℝ} (hin : sq c θ 1 ⊆ rbox L H) (hH : H < 2)
    (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hl : ∀ p ∈ sq c θ 1, lo ≤ p.2 → p.2 ≤ lo + 1 → 1 ≤ p.1)
    (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hs : ∀ p ∈ sq c θ 1, S * p.1 + C * p.2 ≤ s) :
    2 * S - C * (H - 1) ≤ s - C * H ∨ C * (S - C) + (1 + C + S - H) ≤ C * (s - C * H) := by
  refine osl_set (sq_vmir_subset hin) hH (lo := H - 1 - lo) (by linarith) (by linarith)
    (fun q hq h1 h2 => ?_) hC hS hCS (fun q hq => ?_)
  · rw [← vmir_vmir H q, mem_sq_vmir] at hq
    have := hl (vmir H q) hq (by simp only [vmir]; linarith) (by simp only [vmir]; linarith)
    simpa only [vmir] using this
  · rw [← vmir_vmir H q, mem_sq_vmir] at hq
    have := hs _ hq; simp only [vmir] at this; linarith

/-- A square with a right shelf in the half-plane `S x - C y ≥ s` (half turn). -/
theorem osl_set_rot {c : ℝ × ℝ} {θ L H lo C S s : ℝ} (hin : sq c θ 1 ⊆ rbox L H) (hH : H < 2)
    (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hr : ∀ p ∈ sq c θ 1, lo ≤ p.2 → p.2 ≤ lo + 1 → p.1 ≤ L - 1)
    (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hs : ∀ p ∈ sq c θ 1, s ≤ S * p.1 + -C * p.2) :
    2 * S - C * (H - 1) ≤ S * L - C * H - s ∨
      C * (S - C) + (1 + C + S - H) ≤ C * (S * L - C * H - s) := by
  refine osl_set (sq_rot_subset hin) hH (lo := H - 1 - lo) (by linarith) (by linarith)
    (fun q hq h1 h2 => ?_) hC hS hCS (fun q hq => ?_)
  · rw [← rot_rot L H q, mem_sq_rot] at hq
    have := hr _ hq (by simp only [rot]; linarith) (by simp only [rot]; linarith)
    simp only [rot] at this; linarith
  · rw [← rot_rot L H q, mem_sq_rot] at hq
    have := hs _ hq; simp only [rot] at this; linarith

/-- A square with a right shelf in the half-plane `S x + C y ≥ s` (horizontal mirror). -/
theorem osl_set_mir {c : ℝ × ℝ} {θ L H lo C S s : ℝ} (hin : sq c θ 1 ⊆ rbox L H) (hH : H < 2)
    (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hr : ∀ p ∈ sq c θ 1, lo ≤ p.2 → p.2 ≤ lo + 1 → p.1 ≤ L - 1)
    (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hs : ∀ p ∈ sq c θ 1, s ≤ S * p.1 + C * p.2) :
    2 * S - C * (H - 1) ≤ S * L - s ∨ C * (S - C) + (1 + C + S - H) ≤ C * (S * L - s) := by
  refine osl_set (sq_mir_subset hin) hH (lo := lo) hlo hlo' (fun q hq h1 h2 => ?_) hC hS hCS
    (fun q hq => ?_)
  · rw [← mir_mir L q, mem_sq_mir] at hq
    have := hr _ hq h1 h2; simp only [mir] at this; linarith
  · rw [← mir_mir L q, mem_sq_mir] at hq
    have := hs _ hq; simp only [mir] at this; linarith

/-- Adding the two halves. -/
theorem combine4 {C S H L X Y : ℝ} (hC : 0 ≤ C) (hS : 0 < S) (hCS : C ^ 2 + S ^ 2 = 1) (hH : H < 2)
    (hXY : X + Y = S * L - C * H)
    (hX : 2 * S - C * (H - 1) ≤ X ∨ C * (S - C) + (1 + C + S - H) ≤ C * X)
    (hY : 2 * S - C * (H - 1) ≤ Y ∨ C * (S - C) + (1 + C + S - H) ≤ C * Y) :
    4 ≤ L ∨ (0 < C ∧ 2 * C * (S - C) + 2 * (1 + C + S - H) + C ^ 2 * H ≤ S * C * L) := by
  have hS1 : S ≤ 1 := by nlinarith
  -- `C = 0` leaves only the upright branch
  rcases eq_or_lt_of_le hC with h0 | hC0
  · subst h0
    have hS1' : S = 1 := by nlinarith
    subst hS1'
    left
    rcases hX with hX | hX <;> rcases hY with hY | hY <;> nlinarith
  rcases le_total (C * (S - C) + (1 + C + S - H)) (C * (2 * S - C * (H - 1))) with hψ | hψ
  · right
    refine ⟨hC0, ?_⟩
    have hx : C * (S - C) + (1 + C + S - H) ≤ C * X := by
      rcases hX with hX | hX
      · nlinarith [mul_le_mul_of_nonneg_left hX hC]
      · exact hX
    have hy : C * (S - C) + (1 + C + S - H) ≤ C * Y := by
      rcases hY with hY | hY
      · nlinarith [mul_le_mul_of_nonneg_left hY hC]
      · exact hY
    nlinarith
  · left
    have hx : 2 * S - C * (H - 1) ≤ X := by
      rcases hX with hX | hX
      · exact hX
      · exact le_of_mul_le_mul_left (by linarith) hC0
    have hy : 2 * S - C * (H - 1) ≤ Y := by
      rcases hY with hY | hY
      · exact hY
      · exact le_of_mul_le_mul_left (by linarith) hC0
    have : 4 * S ≤ S * L := by nlinarith
    nlinarith

/-! ## Four squares -/

/-- **Four unit squares, lower bound.**  If four unit squares pack `[0,L] × [0,H]` with `H < 2`, then
`L ≥ 4`, or `S C L ≥ 2C(S - C) + 2(1 + C + S - H) + C² H` for some unit `(C, S)`, `C, S > 0`: that is,
`L ≥ L(t, H) = 2 + (1 - 2(H - k - 1) sin t)/cos t + (1 - (H - 2k) cos t)/sin t`, `k = (cos t + sin t)/2`,
Friedman's piece of `S₃` (mostly Morandi).  `(C, S)` is the direction of a line separating the middle two
squares. -/
theorem rect4_lower {L H : ℝ} (hp : RectPacks 4 L H) (hH : H < 2) :
    4 ≤ L ∨ ∃ C S : ℝ, 0 < C ∧ 0 < S ∧ C ^ 2 + S ^ 2 = 1 ∧
      2 * C * (S - C) + 2 * (1 + C + S - H) + C ^ 2 * H ≤ S * C * L := by
  obtain ⟨ctr, ang, hin, hd⟩ := hp
  set σ := Tuple.sort (fun i => (ctr i).1)
  have hm := Tuple.monotone_sort (fun i => (ctr i).1)
  have hne : ∀ i j : Fin 4, i ≠ j → σ i ≠ σ j := fun i j h => σ.injective.ne h
  have h01 : (ctr (σ 0)).1 ≤ (ctr (σ 1)).1 := hm (by decide)
  have h12 : (ctr (σ 1)).1 ≤ (ctr (σ 2)).1 := hm (by decide)
  have h23 : (ctr (σ 2)).1 ≤ (ctr (σ 3)).1 := hm (by decide)
  obtain ⟨lo, l0, l1, hl⟩ := left_shelf (hin (σ 0)) (hin (σ 1)) (hd _ _ (hne 0 1 (by decide))) h01 hH
  obtain ⟨lo2, r0, r1, hr⟩ :=
    right_shelf (hin (σ 2)) (hin (σ 3)) (hd _ _ (hne 2 3 (by decide))) h23 hH
  obtain ⟨u1, u2, s, hu1, hb, hc⟩ :=
    separate_lr (hin (σ 1)) (hin (σ 2)) (hd _ _ (hne 1 2 (by decide))) h12 hH
  set r := Real.sqrt (u1 ^ 2 + u2 ^ 2)
  have hr0 : 0 < r := Real.sqrt_pos.2 (by positivity)
  have hr2 : r ^ 2 = u1 ^ 2 + u2 ^ 2 := Real.sq_sqrt (by positivity)
  have hS : 0 < u1 / r := div_pos hu1 hr0
  have hC : 0 ≤ |u2| / r := div_nonneg (abs_nonneg _) hr0.le
  have hCS : (|u2| / r) ^ 2 + (u1 / r) ^ 2 = 1 := by
    rw [div_pow, div_pow, sq_abs]; field_simp; linarith
  have scale : ∀ a b : ℝ, a ≤ b → a / r ≤ b / r := fun a b h => div_le_div_of_nonneg_right h hr0.le
  have fin := fun X Y hXY hX hY => combine4 (L := L) (X := X) (Y := Y) hC hS hCS hH hXY hX hY
  have out : ∀ C S : ℝ, (4 ≤ L ∨ (0 < C ∧ 2 * C * (S - C) + 2 * (1 + C + S - H) + C ^ 2 * H ≤ S * C * L)) →
      0 < S → C ^ 2 + S ^ 2 = 1 → (4 ≤ L ∨ ∃ C S : ℝ, 0 < C ∧ 0 < S ∧ C ^ 2 + S ^ 2 = 1 ∧
        2 * C * (S - C) + 2 * (1 + C + S - H) + C ^ 2 * H ≤ S * C * L) := by
    rintro C S (h | ⟨h1, h2⟩) hS hCS
    · exact Or.inl h
    · exact Or.inr ⟨C, S, h1, hS, hCS, h2⟩
  refine out _ _ ?_ hS hCS
  rcases le_total u2 0 with hu2 | hu2
  · have e : ∀ p : ℝ × ℝ, u1 / r * p.1 + -(|u2| / r) * p.2 = (u1 * p.1 + u2 * p.2) / r := by
      intro p; rw [abs_of_nonpos hu2]; field_simp
    have hX := osl_set (hin (σ 1)) hH l0 l1 hl hC hS hCS (s := s / r)
      (fun p hp => by rw [e]; exact scale _ _ (hb p hp))
    have hY := osl_set_rot (hin (σ 2)) hH r0 r1 hr hC hS hCS (s := s / r)
      (fun p hp => by rw [e]; exact scale _ _ (hc p hp))
    exact fin _ _ (by ring) hX hY
  · have e : ∀ p : ℝ × ℝ, u1 / r * p.1 + |u2| / r * p.2 = (u1 * p.1 + u2 * p.2) / r := by
      intro p; rw [abs_of_nonneg hu2]; field_simp
    have hX := osl_set_up (hin (σ 1)) hH l0 l1 hl hC hS hCS (s := s / r)
      (fun p hp => by rw [e]; exact scale _ _ (hb p hp))
    have hY := osl_set_mir (hin (σ 2)) hH r0 r1 hr hC hS hCS (s := s / r)
      (fun p hp => by rw [e]; exact scale _ _ (hc p hp))
    exact fin _ _ (by ring) hX hY

end SquarePacking
