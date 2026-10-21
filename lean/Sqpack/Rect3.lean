import Sqpack.Rect

/-!
# Three unit squares in a rectangle

`rect3_lower`: if three unit squares (any angles) pack `[0,L] × [0,H]` with `H < 2`, then
`L ≥ 3` or `L + 2H ≥ 4 + 2√2`.  With the constructions of `Rect.lean` this gives `R₃` exactly
(`rect3_iff`).  The boundary (with the wedge) is the boundary of `S₂` in Friedman's Math Magic, June 2010,
given there by constructions; this file proves the lower bound.

The proof is a *chain* argument, no case enumeration over separating axes:
1. Sort the squares by the x-coordinate of their centres: `a`, `b`, `c`.
2. `left_shelf`: Hahn–Banach separates `a` from `b` by a line with normal `u`, `u₁ > 0` (the vertical
   centre offset is `< 1`).  An *upright* unit square against the left wall at a suitable height lies
   on `a`'s side as well (`upright_repl`, a polynomial inequality), so every point of `b` at the heights
   of that upright square has `x ≥ 1`.  `right_shelf` is the mirror image.
3. `middle`: one square `b` with two such shelves.  Its witness point on each side is its side vertex
   or a point on an adjacent edge at a shelf's top or bottom height; the resulting `3 × 3` real
   inequalities are `case_*`.  Equality: all squares upright (`L = 3`), or the wedge (`L + 2H = 4 + 2√2`,
   `b` at `45°`).
-/

namespace SquarePacking

/-! ## Real inequalities -/

/-- Limit along a segment: `a + t b < s` for all `t ∈ [0,1)` gives `a + b ≤ s`. -/
theorem lim_le (a b s : ℝ) (h : ∀ t : ℝ, 0 ≤ t → t < 1 → a + t * b < s) : a + b ≤ s := by
  by_contra hc; push Not at hc
  have h0 := h 0 le_rfl one_pos; simp at h0
  have hb : 0 < b := by linarith
  set d := a + b - s
  by_cases ht : d ≤ 2 * b
  · have := h (1 - d / (2 * b)) (by rw [sub_nonneg, div_le_one (by linarith)]; exact ht)
      (by have : 0 < d / (2 * b) := div_pos (by linarith) (by linarith); linarith)
    have e : (1 - d / (2 * b)) * b = b - d / 2 := by field_simp
    rw [e] at this; linarith
  · linarith

/-- Upright replacement, core case `C, S, u₂ ≥ 0`. -/
theorem upright_core (C S u1 u2 : ℝ) (hCS : C^2 + S^2 = 1) (hu1 : 0 ≤ u1) (hC : 0 ≤ C) (hS : 0 ≤ S)
    (hu2 : 0 ≤ u2) :
    2 * (u1 + u2) ≤ (C + S) * (u1 + u2) + |u1 * C + u2 * S| + |-u1 * S + u2 * C| := by
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  rw [abs_of_nonneg (by positivity : 0 ≤ u1 * C + u2 * S)]
  rcases le_total 0 (-u1 * S + u2 * C) with hY | hY
  · rw [abs_of_nonneg hY]
    have h1 : 0 ≤ (1 + C) * (u2 * (C + S - 1) - u1 * (1 - C)) := by
      have e : (1 + C) * (u2 * (C + S - 1) - u1 * (1 - C)) = S * (u2 * (1 + C - S) - u1 * S) := by
        linear_combination (u1 + u2) * hCS
      rw [e]; apply mul_nonneg hS; nlinarith [mul_nonneg hu2 (sub_nonneg.2 hS1)]
    have h2 : 0 ≤ u2 * (C + S - 1) - u1 * (1 - C) := by
      by_contra h; push Not at h; nlinarith
    nlinarith
  · rw [abs_of_nonpos hY]
    have h1 : 0 ≤ (1 + S) * (u1 * (C + S - 1) - u2 * (1 - S)) := by
      have e : (1 + S) * (u1 * (C + S - 1) - u2 * (1 - S)) = C * (u1 * (1 + S - C) - u2 * C) := by
        linear_combination (u1 + u2) * hCS
      rw [e]; apply mul_nonneg hC; nlinarith [mul_nonneg hu1 (sub_nonneg.2 hC1)]
    have h2 : 0 ≤ u1 * (C + S - 1) - u2 * (1 - S) := by
      by_contra h; push Not at h; nlinarith
    nlinarith

/-- **Upright replacement.**  For a unit square with axis `(C, S)` and a direction `u` with `u₁ ≥ 0`:
`2 (u₁ + |u₂|) ≤ (|C| + |S|)(u₁ + |u₂|) + |u·e| + |u·e⊥|`.  (An upright square pushed to the wall is
on the same side of any line `u·p = s`, `u₁ > 0`, as the tilted one.) -/
theorem upright_repl (C S u1 u2 : ℝ) (hCS : C^2 + S^2 = 1) (hu1 : 0 ≤ u1) :
    2 * (u1 + |u2|) ≤ (|C| + |S|) * (u1 + |u2|) + |u1 * C + u2 * S| + |-u1 * S + u2 * C| := by
  -- reduce to `u₂ ≥ 0` (reflect `S`), then to `C, S ≥ 0` (quarter turns)
  have red : ∀ C S u2 : ℝ, C^2 + S^2 = 1 → 0 ≤ u2 →
      2 * (u1 + u2) ≤ (|C| + |S|) * (u1 + u2) + |u1 * C + u2 * S| + |-u1 * S + u2 * C| := by
    intro C S u2 hCS hu2
    rcases le_total 0 C with hC | hC <;> rcases le_total 0 S with hS | hS
    · rw [abs_of_nonneg hC, abs_of_nonneg hS]; exact upright_core C S u1 u2 hCS hu1 hC hS hu2
    · have := upright_core (-S) C u1 u2 (by linear_combination hCS) hu1 (by linarith) hC hu2
      rw [abs_of_nonneg hC, abs_of_nonpos hS]
      have e1 : |u1 * (-S) + u2 * C| = |-u1 * S + u2 * C| := by congr 1; ring
      have e2 : |-u1 * C + u2 * (-S)| = |u1 * C + u2 * S| := by
        rw [show -u1 * C + u2 * (-S) = -(u1 * C + u2 * S) by ring, abs_neg]
      linarith
    · have := upright_core S (-C) u1 u2 (by linear_combination hCS) hu1 hS (by linarith) hu2
      rw [abs_of_nonpos hC, abs_of_nonneg hS]
      have e1 : |u1 * S + u2 * (-C)| = |-u1 * S + u2 * C| := by
        rw [show u1 * S + u2 * (-C) = -(-u1 * S + u2 * C) by ring, abs_neg]
      have e2 : |-u1 * (-C) + u2 * S| = |u1 * C + u2 * S| := by congr 1; ring
      linarith
    · have := upright_core (-C) (-S) u1 u2 (by linear_combination hCS) hu1 (by linarith) (by linarith) hu2
      rw [abs_of_nonpos hC, abs_of_nonpos hS]
      have e1 : |u1 * (-C) + u2 * (-S)| = |u1 * C + u2 * S| := by
        rw [show u1 * (-C) + u2 * (-S) = -(u1 * C + u2 * S) by ring, abs_neg]
      have e2 : |-u1 * (-S) + u2 * (-C)| = |-u1 * S + u2 * C| := by
        rw [show -u1 * (-S) + u2 * (-C) = -(-u1 * S + u2 * C) by ring, abs_neg]
      linarith
  rcases le_total 0 u2 with hu2 | hu2
  · rw [abs_of_nonneg hu2]; exact red C S u2 hCS hu2
  · have := red C (-S) (-u2) (by linear_combination hCS) (by linarith)
    rw [abs_of_nonpos hu2]
    have e0 : |-S| = |S| := abs_neg S
    have e1 : |u1 * C + -u2 * -S| = |u1 * C + u2 * S| := by congr 1; ring
    have e2 : |-u1 * -S + -u2 * C| = |-u1 * S + u2 * C| := by
      rw [show -u1 * -S + -u2 * C = -(-u1 * S + u2 * C) by ring, abs_neg]
    rw [e0, e1, e2] at this; exact this

/-- `|u₁|, |u₂| ≤ |u·e| + |u·e⊥|` for a unit vector `e = (C, S)`. -/
theorem abs_le_supp (C S u1 u2 : ℝ) (hCS : C^2 + S^2 = 1) :
    |u1| ≤ |u1 * C + u2 * S| + |-u1 * S + u2 * C| ∧ |u2| ≤ |u1 * C + u2 * S| + |-u1 * S + u2 * C| := by
  have hC1 : |C| ≤ 1 := by rw [abs_le]; constructor <;> nlinarith
  have hS1 : |S| ≤ 1 := by rw [abs_le]; constructor <;> nlinarith
  set X := u1 * C + u2 * S; set Y := -u1 * S + u2 * C
  have e1 : u1 = C * X - S * Y := by simp only [X, Y]; linear_combination (-u1) * hCS
  have e2 : u2 = S * X + C * Y := by simp only [X, Y]; linear_combination (-u2) * hCS
  have hX := abs_nonneg X; have hY := abs_nonneg Y
  constructor
  · rw [e1]
    calc |C * X - S * Y| ≤ |C| * |X| + |S| * |Y| := by rw [← abs_mul, ← abs_mul]; exact abs_sub _ _
      _ ≤ |X| + |Y| := by nlinarith [mul_le_mul_of_nonneg_right hC1 hX, mul_le_mul_of_nonneg_right hS1 hY]
  · rw [e2]
    calc |S * X + C * Y| ≤ |S| * |X| + |C| * |Y| := by rw [← abs_mul, ← abs_mul]; exact abs_add_le _ _
      _ ≤ |X| + |Y| := by nlinarith [mul_le_mul_of_nonneg_right hC1 hY, mul_le_mul_of_nonneg_right hS1 hX]

/-! ### The final `3 × 3` cases

Variables: the middle square's frame `(C, S)` (`C, S ≥ 0`), centre height `y`, shelf bottoms `lo`
(left) and `lo2` (right), and `A = x - 1 - k`, `B = (L - x) - 1 - k` with `k = (C + S)/2`, so that
`L = A + B + 2 + C + S`.  Left: vertex (`V`), edge above it (`U`), edge below it (`D`); same on the right. -/

section cases
set_option linter.unusedSectionVars false
variable (C S H y lo lo2 A B s : ℝ) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C^2 + S^2 = 1)
    (hH : H < 2) (hy : (C+S)/2 ≤ y) (hy2 : y ≤ H - (C+S)/2) (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1)
    (hlo2 : 0 ≤ lo2) (hlo2' : lo2 ≤ H - 1) (hs : 0 < s) (hs2 : s^2 = 2)
include hC hS hCS hH hy hy2 hlo hlo' hlo2 hlo2' hs hs2

/-- The goal of every case. -/
abbrev Goal3 := 3 ≤ A + B + 2 + C + S ∨ 4 + 2 * s ≤ A + B + 2 + C + S + 2 * H

/-- Core of the two wedge cases: `C S (A + B) ≥ C³ + S³ - G` gives `A + B ≥ 2√2 - 2G - C - S`. -/
theorem wedge_core (G : ℝ) (hC0 : 0 < C) (hS0 : 0 < S) (hG : G < 1)
    (h : C ^ 3 + S ^ 3 - G ≤ C * S * (A + B)) : 2 * s - 2 * G - C - S ≤ A + B := by
  have hp : C + S ≤ s := by nlinarith [sq_nonneg (C - S), sq_nonneg (C + S - s)]
  have hq : 0 ≤ 1 - 2 * (C * S) := by nlinarith [sq_nonneg (C - S)]
  have key : 1 - 2 * (C * S) ≤ C + S - 2 * s * (C * S) := by
    have e : 2 * (C * S) = (C + S) ^ 2 - 1 := by nlinarith
    have : 0 ≤ (C + S - 1) * (1 - (s - 1) * (C + S + 1)) := by
      apply mul_nonneg (by nlinarith [mul_pos hC0 hS0])
      nlinarith
    nlinarith
  have hcs : 0 < C * S := mul_pos hC0 hS0
  have h3 : C ^ 3 + S ^ 3 = (C + S) * (1 - C * S) := by nlinarith
  have : C * S * (2 * s - 2 * G - C - S) ≤ C * S * (A + B) := by
    nlinarith [mul_le_mul_of_nonneg_right hG.le hq]
  exact le_of_mul_le_mul_left this hcs

theorem case_VV (hA : 0 ≤ A) (hB : 0 ≤ B) : Goal3 C S H A B s := by
  left; nlinarith

theorem case_VU (hA : 0 ≤ A) (hR1 : y + (S - C)/2 < lo2) (hR2 : lo2 ≤ y + (C+S)/2)
    (hB : 0 ≤ C * B + S * (lo2 - (y + (S - C)/2))) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

theorem case_VD (hA : 0 ≤ A) (hR1 : lo2 + 1 < y + (S - C)/2) (hR2 : y - (C+S)/2 ≤ lo2 + 1)
    (hB : 0 ≤ S * B + C * (y + (S - C)/2 - lo2 - 1)) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

theorem case_UV (hL1 : y + (C - S)/2 < lo) (hL2 : lo ≤ y + (C+S)/2)
    (hA : 0 ≤ S * A + C * (lo - (y + (C - S)/2))) (hB : 0 ≤ B) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

theorem case_DV (hL1 : lo + 1 < y + (C - S)/2) (hL2 : y - (C+S)/2 ≤ lo + 1)
    (hA : 0 ≤ C * A + S * (y + (C - S)/2 - lo - 1)) (hB : 0 ≤ B) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

theorem case_UD (hL1 : y + (C - S)/2 < lo) (hL2 : lo ≤ y + (C+S)/2)
    (hA : 0 ≤ S * A + C * (lo - (y + (C - S)/2)))
    (_hR1 : lo2 + 1 < y + (S - C)/2) (_hR2 : y - (C+S)/2 ≤ lo2 + 1)
    (hB : 0 ≤ S * B + C * (y + (S - C)/2 - lo2 - 1)) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

theorem case_DU (hL1 : lo + 1 < y + (C - S)/2) (hL2 : y - (C+S)/2 ≤ lo + 1)
    (hA : 0 ≤ C * A + S * (y + (C - S)/2 - lo - 1))
    (hR1 : y + (S - C)/2 < lo2) (_hR2 : lo2 ≤ y + (C+S)/2)
    (hB : 0 ≤ C * B + S * (lo2 - (y + (S - C)/2))) : Goal3 C S H A B s := by
  by_contra hc; push Not at hc
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  nlinarith [mul_nonneg hC hS, mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)]

/-- The wedge case: `b` hangs from the top with both shelves above its side vertices. -/
theorem case_UU (hL1 : y + (C - S)/2 < lo) (hL2 : lo ≤ y + (C+S)/2)
    (hA : 0 ≤ S * A + C * (lo - (y + (C - S)/2)))
    (hR1 : y + (S - C)/2 < lo2) (hR2 : lo2 ≤ y + (C+S)/2)
    (hB : 0 ≤ C * B + S * (lo2 - (y + (S - C)/2))) : Goal3 C S H A B s := by
  right
  have hX : lo - (y + (C - S)/2) ≤ H - 1 - C := by linarith
  have hY : lo2 - (y + (S - C)/2) ≤ H - 1 - S := by linarith
  have := wedge_core C S H y lo lo2 A B s hC hS hCS hH hy hy2 hlo hlo' hlo2 hlo2' hs hs2 (H - 1)
    (by linarith) (by linarith) (by linarith) (by
      nlinarith [mul_le_mul_of_nonneg_left hX (sq_nonneg C), mul_le_mul_of_nonneg_left hY (sq_nonneg S)])
  linarith

/-- The mirrored wedge case. -/
theorem case_DD (hL1 : lo + 1 < y + (C - S)/2) (hL2 : y - (C+S)/2 ≤ lo + 1)
    (hA : 0 ≤ C * A + S * (y + (C - S)/2 - lo - 1))
    (hR1 : lo2 + 1 < y + (S - C)/2) (hR2 : y - (C+S)/2 ≤ lo2 + 1)
    (hB : 0 ≤ S * B + C * (y + (S - C)/2 - lo2 - 1)) : Goal3 C S H A B s := by
  right
  have hX : y + (C - S)/2 - lo - 1 ≤ H - 1 - S := by linarith
  have hY : y + (S - C)/2 - lo2 - 1 ≤ H - 1 - C := by linarith
  have := wedge_core C S H y lo lo2 A B s hC hS hCS hH hy hy2 hlo hlo' hlo2 hlo2' hs hs2 (H - 1)
    (by linarith) (by linarith) (by linarith) (by
      nlinarith [mul_le_mul_of_nonneg_left hX (sq_nonneg S), mul_le_mul_of_nonneg_left hY (sq_nonneg C)])
  linarith

end cases

/-! ## Squares: interiors, separation, support -/

lemma coord_shrink (c p : ℝ × ℝ) (θ t : ℝ) :
    coord c θ (c.1 + t * (p.1 - c.1), c.2 + t * (p.2 - c.2)) =
      (t * (coord c θ p).1, t * (coord c θ p).2) := by
  simp only [coord]; ext <;> ring

lemma sqInt_shrink {c p : ℝ × ℝ} {θ t : ℝ} (hp : p ∈ sq c θ 1) (ht0 : 0 ≤ t) (ht1 : t < 1) :
    (c.1 + t * (p.1 - c.1), c.2 + t * (p.2 - c.2)) ∈ sqInt c θ 1 := by
  obtain ⟨h1, h2⟩ := hp
  simp only [sqInt, Set.mem_ofPred_eq, coord_shrink, abs_mul, abs_of_nonneg ht0]
  constructor <;> nlinarith

/-- A linear bound on the open square extends to the closed square. -/
lemma lin_le_of_sqInt {c : ℝ × ℝ} {θ u1 u2 s : ℝ}
    (h : ∀ q ∈ sqInt c θ 1, u1 * q.1 + u2 * q.2 < s) : ∀ p ∈ sq c θ 1, u1 * p.1 + u2 * p.2 ≤ s := by
  intro p hp
  have := lim_le (u1 * c.1 + u2 * c.2) (u1 * (p.1 - c.1) + u2 * (p.2 - c.2)) s (fun t ht0 ht1 => by
    have := h _ (sqInt_shrink hp ht0 ht1); simp only at this; linarith)
  linarith

lemma centre_mem_sqInt' (c : ℝ × ℝ) (θ : ℝ) : c ∈ sqInt c θ 1 := by
  simp [sqInt, coord]

lemma sqInt_isOpen (c : ℝ × ℝ) (θ : ℝ) : IsOpen (sqInt c θ 1) := by
  have h1 : Continuous fun p : ℝ × ℝ => |(coord c θ p).1| := by unfold coord; fun_prop
  have h2 : Continuous fun p : ℝ × ℝ => |(coord c θ p).2| := by unfold coord; fun_prop
  exact (isOpen_lt h1 continuous_const).inter (isOpen_lt h2 continuous_const)

lemma sqInt_convex (c : ℝ × ℝ) (θ : ℝ) : Convex ℝ (sqInt c θ 1) := by
  intro p hp q hq a b ha hb hab
  obtain ⟨hp1, hp2⟩ := hp; obtain ⟨hq1, hq2⟩ := hq
  have e : coord c θ (a • p + b • q) = (a * (coord c θ p).1 + b * (coord c θ q).1,
      a * (coord c θ p).2 + b * (coord c θ q).2) := by
    simp only [coord, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd, smul_eq_mul]
    ext
    · simp only; linear_combination (c.1 * Real.cos θ + c.2 * Real.sin θ) * hab
    · simp only; linear_combination (-c.1 * Real.sin θ + c.2 * Real.cos θ) * hab
  have key : ∀ x y : ℝ, |x| < 1 / 2 → |y| < 1 / 2 → |a * x + b * y| < 1 / 2 := by
    intro x y hx hy
    calc |a * x + b * y| ≤ a * |x| + b * |y| := by
          calc |a * x + b * y| ≤ |a * x| + |b * y| := abs_add_le _ _
            _ = a * |x| + b * |y| := by rw [abs_mul, abs_mul, abs_of_nonneg ha, abs_of_nonneg hb]
      _ < 1 / 2 := by
          rcases eq_or_lt_of_le ha with h0 | h0
          · subst h0; simp at hab; subst hab; simpa using hy
          · nlinarith [abs_nonneg x, abs_nonneg y]
  simp only [sqInt, Set.mem_ofPred_eq, e]
  exact ⟨key _ _ hp1 hq1, key _ _ hp2 hq2⟩

/-- **Separation.**  Two unit squares with disjoint interiors lie on opposite sides of a line. -/
theorem separate {ca cb : ℝ × ℝ} {θa θb : ℝ} (hd : Disjoint (sqInt ca θa 1) (sqInt cb θb 1)) :
    ∃ u1 u2 s : ℝ, (u1 ≠ 0 ∨ u2 ≠ 0) ∧ (∀ p ∈ sq ca θa 1, u1 * p.1 + u2 * p.2 ≤ s) ∧
      (∀ p ∈ sq cb θb 1, s ≤ u1 * p.1 + u2 * p.2) := by
  obtain ⟨f, s, hf1, hf2⟩ := geometric_hahn_banach_open_open (sqInt_convex ca θa) (sqInt_isOpen ca θa)
    (sqInt_convex cb θb) (sqInt_isOpen cb θb) hd
  have hf : ∀ p : ℝ × ℝ, f p = f (1, 0) * p.1 + f (0, 1) * p.2 := by
    intro p
    have : p = p.1 • ((1 : ℝ), (0 : ℝ)) + p.2 • ((0 : ℝ), (1 : ℝ)) := by ext <;> simp
    conv_lhs => rw [this]
    rw [map_add, map_smul, map_smul, smul_eq_mul, smul_eq_mul]; ring
  refine ⟨f (1, 0), f (0, 1), s, ?_, ?_, ?_⟩
  · by_contra h0; push Not at h0
    have a := hf1 ca (centre_mem_sqInt' ca θa); have b := hf2 cb (centre_mem_sqInt' cb θb)
    rw [hf, h0.1, h0.2] at a b; simp at a b; linarith
  · exact lin_le_of_sqInt (fun q hq => by rw [← hf]; exact hf1 q hq)
  · intro p hp
    have := lin_le_of_sqInt (u1 := -f (1, 0)) (u2 := -f (0, 1)) (s := -s)
      (fun q hq => by have := hf2 q hq; rw [hf] at this; linarith) p hp
    linarith

/-- `|a|/2` as `±a/2` with the sign chosen by `a`. -/
lemma sgn_half' (a : ℝ) : ∃ t : ℝ, |t| = 1 / 2 ∧ t * a = |a| / 2 := by
  rcases le_or_gt 0 a with h | h
  · exact ⟨1 / 2, by norm_num, by rw [abs_of_nonneg h]; ring⟩
  · exact ⟨-(1 / 2), by norm_num, by rw [abs_of_neg h]; ring⟩

/-- **Support.**  A closed unit square below the line `u·p = s` has
`u·c + (|u·e| + |u·e⊥|)/2 ≤ s`. -/
theorem support_le {c : ℝ × ℝ} {θ u1 u2 s : ℝ} (h : ∀ p ∈ sq c θ 1, u1 * p.1 + u2 * p.2 ≤ s) :
    u1 * c.1 + u2 * c.2 + (|u1 * Real.cos θ + u2 * Real.sin θ| +
      |-u1 * Real.sin θ + u2 * Real.cos θ|) / 2 ≤ s := by
  obtain ⟨a, ha, hax⟩ := sgn_half' (u1 * Real.cos θ + u2 * Real.sin θ)
  obtain ⟨b, hb, hby⟩ := sgn_half' (-u1 * Real.sin θ + u2 * Real.cos θ)
  have := h _ (pt_mem_sq c θ a b (le_of_eq ha) (le_of_eq hb))
  simp only at this
  nlinarith [hax, hby]

theorem support_ge {c : ℝ × ℝ} {θ u1 u2 s : ℝ} (h : ∀ p ∈ sq c θ 1, s ≤ u1 * p.1 + u2 * p.2) :
    s ≤ u1 * c.1 + u2 * c.2 - (|u1 * Real.cos θ + u2 * Real.sin θ| +
      |-u1 * Real.sin θ + u2 * Real.cos θ|) / 2 := by
  have := support_le (c := c) (θ := θ) (u1 := -u1) (u2 := -u2) (s := -s)
    (fun p hp => by have := h p hp; linarith)
  rw [show -u1 * Real.cos θ + -u2 * Real.sin θ = -(u1 * Real.cos θ + u2 * Real.sin θ) by ring,
    show - -u1 * Real.sin θ + -u2 * Real.cos θ = -(-u1 * Real.sin θ + u2 * Real.cos θ) by ring,
    abs_neg, abs_neg] at this
  linarith

/-- The separating normal points right when the centres are ordered in `x` and the vertical
offset is below `1`. -/
theorem u1_pos {u1 u2 xa xb ya yb Sa Sb s : ℝ} (hu : u1 ≠ 0 ∨ u2 ≠ 0) (hx : xa ≤ xb)
    (hdy : |yb - ya| < 1)
    (ha : u1 * xa + u2 * ya + Sa / 2 ≤ s) (hb : s ≤ u1 * xb + u2 * yb - Sb / 2)
    (ha1 : |u1| ≤ Sa) (ha2 : |u2| ≤ Sa) (hb1 : |u1| ≤ Sb) (hb2 : |u2| ≤ Sb) : 0 < u1 := by
  by_contra hc; push Not at hc
  set dx := xb - xa; set dy := yb - ya
  have h : (Sa + Sb) / 2 ≤ u1 * dx + u2 * dy := by simp only [dx, dy]; linarith
  have h1 : u1 * dx ≤ 0 := mul_nonpos_of_nonpos_of_nonneg hc (by simp only [dx]; linarith)
  have h2 : u2 * dy ≤ |u2| * |dy| := by rw [← abs_mul]; exact le_abs_self _
  rcases hu with hu | hu
  · have : 0 < |u1| := abs_pos.2 hu
    nlinarith [abs_nonneg dy, abs_nonneg u2]
  · have : 0 < |u2| := abs_pos.2 hu
    nlinarith [abs_nonneg dy]

/-! ## Shelves -/

/-- The half extent `k ≥ 1/2`. -/
lemma half_le_k (θ : ℝ) : 1 / 2 ≤ (|Real.cos θ| + |Real.sin θ|) / 2 := by
  have := one_le_abs_cos_add_abs_sin θ; linarith

/-- **Left shelf.**  If `a` is left of `b` (by centre) in a strip of height `H < 2`, some height band
`[lo, lo + 1] ⊆ [0, H]` has every point of `b` at those heights at `x ≥ 1`. -/
theorem left_shelf {ca cb : ℝ × ℝ} {θa θb L H : ℝ} (hA : sq ca θa 1 ⊆ rbox L H)
    (hB : sq cb θb 1 ⊆ rbox L H) (hd : Disjoint (sqInt ca θa 1) (sqInt cb θb 1))
    (hx : ca.1 ≤ cb.1) (hH : H < 2) :
    ∃ lo, 0 ≤ lo ∧ lo ≤ H - 1 ∧ ∀ p ∈ sq cb θb 1, lo ≤ p.2 → p.2 ≤ lo + 1 → 1 ≤ p.1 := by
  obtain ⟨u1, u2, s, hu, hsa, hsb⟩ := separate hd
  have sa := support_le hsa
  have sb := support_ge hsb
  obtain ⟨a1, a2, a3, a4⟩ := half_extent hA
  obtain ⟨b1, b2, b3, b4⟩ := half_extent hB
  have ka := half_le_k θa; have kb := half_le_k θb
  have hCSa : Real.cos θa ^ 2 + Real.sin θa ^ 2 = 1 := Real.cos_sq_add_sin_sq θa
  have hCSb : Real.cos θb ^ 2 + Real.sin θb ^ 2 = 1 := Real.cos_sq_add_sin_sq θb
  obtain ⟨la1, la2⟩ := abs_le_supp _ _ u1 u2 hCSa
  obtain ⟨lb1, lb2⟩ := abs_le_supp _ _ u1 u2 hCSb
  have hu1 : 0 < u1 := u1_pos hu hx (by rw [abs_lt]; constructor <;> linarith) sa sb la1 la2 lb1 lb2
  have rep := upright_repl _ _ u1 u2 hCSa hu1.le
  set k := (|Real.cos θa| + |Real.sin θa|) / 2 with hk
  -- the upright square `[0,1] × [lo, lo+1]` is below the line
  have key : ∀ lo : ℝ, (0 ≤ u2 → lo = ca.2 - k) → (u2 < 0 → lo = ca.2 + k - 1) →
      ∀ y, lo ≤ y → y ≤ lo + 1 → u1 + u2 * y ≤ s := by
    intro lo h1 h2 y hy1 hy2
    have hxa : u1 * k ≤ u1 * ca.1 := mul_le_mul_of_nonneg_left a1 hu1.le
    rcases le_or_gt 0 u2 with hu2 | hu2
    · have := h1 hu2; subst this
      rw [abs_of_nonneg hu2] at rep
      nlinarith [mul_le_mul_of_nonneg_left hy2 hu2]
    · have := h2 hu2; subst this
      rw [abs_of_neg hu2] at rep
      nlinarith [mul_le_mul_of_nonpos_left hy1 hu2.le]
  have fin : ∀ lo : ℝ, (∀ y, lo ≤ y → y ≤ lo + 1 → u1 + u2 * y ≤ s) →
      ∀ p ∈ sq cb θb 1, lo ≤ p.2 → p.2 ≤ lo + 1 → 1 ≤ p.1 := by
    intro lo hl p hp h1 h2
    have := hl p.2 h1 h2; have := hsb p hp
    by_contra hc; push Not at hc; nlinarith
  rcases le_or_gt 0 u2 with hu2 | hu2
  · exact ⟨ca.2 - k, by linarith, by linarith, fin _ (key _ (fun _ => rfl) (fun h => absurd h (not_lt.2 hu2)))⟩
  · exact ⟨ca.2 + k - 1, by linarith, by linarith,
      fin _ (key _ (fun h => absurd h (not_le.2 hu2)) (fun _ => rfl))⟩

/-- **Right shelf**, the mirror image: points of `b` at the band's heights have `x ≤ L - 1`. -/
theorem right_shelf {cb cc : ℝ × ℝ} {θb θc L H : ℝ} (hB : sq cb θb 1 ⊆ rbox L H)
    (hC : sq cc θc 1 ⊆ rbox L H) (hd : Disjoint (sqInt cb θb 1) (sqInt cc θc 1))
    (hx : cb.1 ≤ cc.1) (hH : H < 2) :
    ∃ lo, 0 ≤ lo ∧ lo ≤ H - 1 ∧ ∀ p ∈ sq cb θb 1, lo ≤ p.2 → p.2 ≤ lo + 1 → p.1 ≤ L - 1 := by
  obtain ⟨u1, u2, s, hu, hsb, hsc⟩ := separate hd
  have sb := support_le hsb
  have sc := support_ge hsc
  obtain ⟨b1, b2, b3, b4⟩ := half_extent hB
  obtain ⟨c1, c2, c3, c4⟩ := half_extent hC
  have kb := half_le_k θb; have kc := half_le_k θc
  have hCSb : Real.cos θb ^ 2 + Real.sin θb ^ 2 = 1 := Real.cos_sq_add_sin_sq θb
  have hCSc : Real.cos θc ^ 2 + Real.sin θc ^ 2 = 1 := Real.cos_sq_add_sin_sq θc
  obtain ⟨lb1, lb2⟩ := abs_le_supp _ _ u1 u2 hCSb
  obtain ⟨lc1, lc2⟩ := abs_le_supp _ _ u1 u2 hCSc
  have hu1 : 0 < u1 := u1_pos hu hx (by rw [abs_lt]; constructor <;> linarith) sb sc lb1 lb2 lc1 lc2
  have rep := upright_repl _ _ u1 u2 hCSc hu1.le
  set k := (|Real.cos θc| + |Real.sin θc|) / 2 with hk
  have key : ∀ lo : ℝ, (0 ≤ u2 → lo = cc.2 + k - 1) → (u2 < 0 → lo = cc.2 - k) →
      ∀ y, lo ≤ y → y ≤ lo + 1 → s ≤ u1 * (L - 1) + u2 * y := by
    intro lo h1 h2 y hy1 hy2
    have hxc : u1 * cc.1 ≤ u1 * (L - k) := mul_le_mul_of_nonneg_left c2 hu1.le
    rcases le_or_gt 0 u2 with hu2 | hu2
    · have := h1 hu2; subst this
      rw [abs_of_nonneg hu2] at rep
      nlinarith [mul_le_mul_of_nonneg_left hy1 hu2]
    · have := h2 hu2; subst this
      rw [abs_of_neg hu2] at rep
      nlinarith [mul_le_mul_of_nonpos_left hy2 hu2.le]
  have fin : ∀ lo : ℝ, (∀ y, lo ≤ y → y ≤ lo + 1 → s ≤ u1 * (L - 1) + u2 * y) →
      ∀ p ∈ sq cb θb 1, lo ≤ p.2 → p.2 ≤ lo + 1 → p.1 ≤ L - 1 := by
    intro lo hl p hp h1 h2
    have := hl p.2 h1 h2; have := hsb p hp
    by_contra hc; push Not at hc; nlinarith
  rcases le_or_gt 0 u2 with hu2 | hu2
  · exact ⟨cc.2 + k - 1, by linarith, by linarith,
      fin _ (key _ (fun _ => rfl) (fun h => absurd h (not_lt.2 hu2)))⟩
  · exact ⟨cc.2 - k, by linarith, by linarith,
      fin _ (key _ (fun h => absurd h (not_le.2 hu2)) (fun _ => rfl))⟩

/-! ## The middle square -/

/-- A frame `(P, Q)`, `P, Q ≥ 0`, of the square: its points are `c + u (P, Q) + v (-Q, P)`. -/
theorem frame (c : ℝ × ℝ) (θ : ℝ) : ∃ P Q : ℝ, 0 ≤ P ∧ 0 ≤ Q ∧ P ^ 2 + Q ^ 2 = 1 ∧
    P + Q = |Real.cos θ| + |Real.sin θ| ∧ ∀ u v : ℝ, |u| ≤ 1 / 2 → |v| ≤ 1 / 2 →
      ((c.1 + (u * P - v * Q), c.2 + (u * Q + v * P)) : ℝ × ℝ) ∈ sq c θ 1 := by
  set C := Real.cos θ; set S := Real.sin θ
  have hCS : C ^ 2 + S ^ 2 = 1 := Real.cos_sq_add_sin_sq θ
  have hn : ∀ x : ℝ, |x| ≤ 1 / 2 → |-x| ≤ 1 / 2 := fun x h => by rwa [abs_neg]
  rcases le_total 0 C with hC | hC <;> rcases le_total 0 S with hS | hS
  · refine ⟨C, S, hC, hS, hCS, by rw [abs_of_nonneg hC, abs_of_nonneg hS], fun u v hu hv => ?_⟩
    exact pt_mem_sq c θ u v hu hv
  · -- (P, Q) = (-S, C):  u' = -v, v' = u
    refine ⟨-S, C, by linarith, hC, by linear_combination hCS,
      by rw [abs_of_nonneg hC, abs_of_nonpos hS]; ring, fun u v hu hv => ?_⟩
    have := pt_mem_sq c θ (-v) u (hn v hv) hu
    convert this using 2 <;> simp only [C, S] <;> ring
  · -- (P, Q) = (S, -C):  u' = v, v' = -u
    refine ⟨S, -C, hS, by linarith, by linear_combination hCS,
      by rw [abs_of_nonpos hC, abs_of_nonneg hS]; ring, fun u v hu hv => ?_⟩
    have := pt_mem_sq c θ v (-u) hv (hn u hu)
    convert this using 2 <;> simp only [C, S] <;> ring
  · refine ⟨-C, -S, by linarith, by linarith, by linear_combination hCS,
      by rw [abs_of_nonpos hC, abs_of_nonpos hS], fun u v hu hv => ?_⟩
    have := pt_mem_sq c θ (-u) (-v) (hn u hu) (hn v hv)
    convert this using 2 <;> simp only [C, S] <;> ring

/-- **The middle square** with a left shelf and a right shelf. -/
theorem middle {cb : ℝ × ℝ} {θb L H lo lo2 : ℝ} (hB : sq cb θb 1 ⊆ rbox L H) (hH : H < 2)
    (hlo : 0 ≤ lo) (hlo' : lo ≤ H - 1) (hlo2 : 0 ≤ lo2) (hlo2' : lo2 ≤ H - 1)
    (hl : ∀ p ∈ sq cb θb 1, lo ≤ p.2 → p.2 ≤ lo + 1 → 1 ≤ p.1)
    (hr : ∀ p ∈ sq cb θb 1, lo2 ≤ p.2 → p.2 ≤ lo2 + 1 → p.1 ≤ L - 1) :
    3 ≤ L ∨ 4 + 2 * Real.sqrt 2 ≤ L + 2 * H := by
  obtain ⟨P, Q, hP, hQ, hPQ, hk, hpt⟩ := frame cb θb
  obtain ⟨-, -, b3, b4⟩ := half_extent hB
  rw [← hk] at b3 b4
  set x := cb.1; set y := cb.2
  have hs : 0 < Real.sqrt 2 := by positivity
  have hs2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hk1 : 1 ≤ P + Q := by nlinarith
  set A := x - 1 - (P + Q) / 2
  set B := (L - x) - 1 - (P + Q) / 2
  have hL : L = A + B + 2 + P + Q := by simp only [A, B]; ring
  rw [hL]
  have h2 : (1 / 2 : ℝ) ≤ 1 / 2 := le_rfl
  have hh : |(1 / 2 : ℝ)| ≤ 1 / 2 := by norm_num
  have hmh : |(-(1 / 2) : ℝ)| ≤ 1 / 2 := by norm_num
  -- left witness
  have left : 0 ≤ A ∨ (y + (P - Q)/2 < lo ∧ lo ≤ y + (P+Q)/2 ∧ 0 ≤ Q * A + P * (lo - (y + (P - Q)/2)))
      ∨ (lo + 1 < y + (P - Q)/2 ∧ y - (P+Q)/2 ≤ lo + 1 ∧ 0 ≤ P * A + Q * (y + (P - Q)/2 - lo - 1)) := by
    have top : lo ≤ y + (P+Q)/2 := by linarith
    have bot : y - (P+Q)/2 ≤ lo + 1 := by linarith
    rcases lt_or_ge (y + (P - Q)/2) lo with h | h
    · -- edge from the left vertex up to the top vertex, at height `lo`
      have hQ0 : 0 < Q := by linarith
      set u := (lo - y - P / 2) / Q
      have hu : |u| ≤ 1 / 2 := by
        rw [abs_le]; constructor
        · rw [le_div_iff₀ hQ0]; linarith
        · rw [div_le_iff₀ hQ0]; linarith
      have hm := hpt u (1 / 2) hu hh
      have huQ : u * Q = lo - y - P / 2 := by simp only [u]; field_simp
      have := hl _ hm (by simp only; linarith) (by simp only; linarith)
      simp only at this
      right; left
      refine ⟨h, top, ?_⟩
      have : Q * (x - 1 - (P + Q) / 2) + P * (lo - (y + (P - Q) / 2)) = Q * (x + u * P - 1 / 2 * Q - 1) := by
        linear_combination (-P) * huQ
      simp only [A]; rw [this]; exact mul_nonneg hQ0.le (by linarith)
    · rcases le_or_gt (y + (P - Q)/2) (lo + 1) with h' | h'
      · -- the left vertex itself
        have hm := hpt (-(1 / 2)) (1 / 2) hmh hh
        have := hl _ hm (by simp only; linarith) (by simp only; linarith)
        simp only at this
        left; simp only [A]; linarith
      · -- edge from the left vertex down to the bottom vertex, at height `lo + 1`
        have hP0 : 0 < P := by linarith
        set v := (lo + 1 - y + Q / 2) / P
        have hv : |v| ≤ 1 / 2 := by
          rw [abs_le]; constructor
          · rw [le_div_iff₀ hP0]; linarith
          · rw [div_le_iff₀ hP0]; linarith
        have hm := hpt (-(1 / 2)) v hmh hv
        have hvP : v * P = lo + 1 - y + Q / 2 := by simp only [v]; field_simp
        have := hl _ hm (by simp only; linarith) (by simp only; linarith)
        simp only at this
        right; right
        refine ⟨h', bot, ?_⟩
        have : P * (x - 1 - (P + Q) / 2) + Q * (y + (P - Q) / 2 - lo - 1) = P * (x + (-(1 / 2) * P - v * Q) - 1) := by
          linear_combination Q * hvP
        simp only [A]; rw [this]; exact mul_nonneg hP0.le (by linarith)
  -- right witness
  have right : 0 ≤ B ∨ (y + (Q - P)/2 < lo2 ∧ lo2 ≤ y + (P+Q)/2 ∧ 0 ≤ P * B + Q * (lo2 - (y + (Q - P)/2)))
      ∨ (lo2 + 1 < y + (Q - P)/2 ∧ y - (P+Q)/2 ≤ lo2 + 1 ∧ 0 ≤ Q * B + P * (y + (Q - P)/2 - lo2 - 1)) := by
    have top : lo2 ≤ y + (P+Q)/2 := by linarith
    have bot : y - (P+Q)/2 ≤ lo2 + 1 := by linarith
    rcases lt_or_ge (y + (Q - P)/2) lo2 with h | h
    · -- edge from the right vertex up to the top vertex (`u = 1/2`), at height `lo2`
      have hP0 : 0 < P := by linarith
      set v := (lo2 - y - Q / 2) / P
      have hv : |v| ≤ 1 / 2 := by
        rw [abs_le]; constructor
        · rw [le_div_iff₀ hP0]; linarith
        · rw [div_le_iff₀ hP0]; linarith
      have hm := hpt (1 / 2) v hh hv
      have hvP : v * P = lo2 - y - Q / 2 := by simp only [v]; field_simp
      have := hr _ hm (by simp only; linarith) (by simp only; linarith)
      simp only at this
      right; left
      refine ⟨h, top, ?_⟩
      have : P * (L - x - 1 - (P + Q) / 2) + Q * (lo2 - (y + (Q - P) / 2)) = P * (L - 1 - (x + (1 / 2 * P - v * Q))) := by
        linear_combination (-Q) * hvP
      simp only [B]; rw [this]; exact mul_nonneg hP0.le (by linarith)
    · rcases le_or_gt (y + (Q - P)/2) (lo2 + 1) with h' | h'
      · have hm := hpt (1 / 2) (-(1 / 2)) hh hmh
        have := hr _ hm (by simp only; linarith) (by simp only; linarith)
        simp only at this
        left; simp only [B]; linarith
      · -- edge from the right vertex down to the bottom vertex (`v = -1/2`), at height `lo2 + 1`
        have hQ0 : 0 < Q := by linarith
        set u := (lo2 + 1 - y + P / 2) / Q
        have hu : |u| ≤ 1 / 2 := by
          rw [abs_le]; constructor
          · rw [le_div_iff₀ hQ0]; linarith
          · rw [div_le_iff₀ hQ0]; linarith
        have hm := hpt u (-(1 / 2)) hu hmh
        have huQ : u * Q = lo2 + 1 - y + P / 2 := by simp only [u]; field_simp
        have := hr _ hm (by simp only; linarith) (by simp only; linarith)
        simp only at this
        right; right
        refine ⟨h', bot, ?_⟩
        have : Q * (L - x - 1 - (P + Q) / 2) + P * (y + (Q - P) / 2 - lo2 - 1) = Q * (L - 1 - (x + (u * P - -(1 / 2) * Q))) := by
          linear_combination P * huQ
        simp only [B]; rw [this]; exact mul_nonneg hQ0.le (by linarith)
  have hy1 : (P + Q) / 2 ≤ y := b3
  have hy2 : y ≤ H - (P + Q) / 2 := b4
  have c := fun hA hB => case_VV P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 hA hB
  rcases left with hA | ⟨l1, l2, l3⟩ | ⟨l1, l2, l3⟩ <;> rcases right with hB | ⟨r1, r2, r3⟩ | ⟨r1, r2, r3⟩
  · exact c hA hB
  · exact case_VU P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 hA r1 r2 r3
  · exact case_VD P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 hA r1 r2 r3
  · exact case_UV P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 hB
  · exact case_UU P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 r1 r2 r3
  · exact case_UD P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 r1 r2 r3
  · exact case_DV P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 hB
  · exact case_DU P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 r1 r2 r3
  · exact case_DD P Q H y lo lo2 A B _ hP hQ hPQ hH hy1 hy2 hlo hlo' hlo2 hlo2' hs hs2 l1 l2 l3 r1 r2 r3

/-! ## Three squares -/

/-- The chain for three squares ordered by centre `x`. -/
theorem chain3 {ca cb cc : ℝ × ℝ} {θa θb θc L H : ℝ} (hA : sq ca θa 1 ⊆ rbox L H)
    (hB : sq cb θb 1 ⊆ rbox L H) (hC : sq cc θc 1 ⊆ rbox L H)
    (hab : Disjoint (sqInt ca θa 1) (sqInt cb θb 1)) (hbc : Disjoint (sqInt cb θb 1) (sqInt cc θc 1))
    (h1 : ca.1 ≤ cb.1) (h2 : cb.1 ≤ cc.1) (hH : H < 2) :
    3 ≤ L ∨ 4 + 2 * Real.sqrt 2 ≤ L + 2 * H := by
  obtain ⟨lo, l0, l1, hl⟩ := left_shelf hA hB hab h1 hH
  obtain ⟨lo2, r0, r1, hr⟩ := right_shelf hB hC hbc h2 hH
  exact middle hB hH l0 l1 r0 r1 hl hr

/-- **Three unit squares, lower bound.**  If three unit squares pack `[0,L] × [0,H]` and `H < 2`, then
`L ≥ 3` or `L + 2H ≥ 4 + 2√2`. -/
theorem rect3_lower {L H : ℝ} (hp : RectPacks 3 L H) (hH : H < 2) :
    3 ≤ L ∨ 4 + 2 * Real.sqrt 2 ≤ L + 2 * H := by
  obtain ⟨ctr, ang, hin, hd⟩ := hp
  have go : ∀ i j k : Fin 3, i ≠ j → j ≠ k → (ctr i).1 ≤ (ctr j).1 → (ctr j).1 ≤ (ctr k).1 →
      3 ≤ L ∨ 4 + 2 * Real.sqrt 2 ≤ L + 2 * H := fun i j k hij hjk h1 h2 =>
    chain3 (hin i) (hin j) (hin k) (hd i j hij) (hd j k hjk) h1 h2 hH
  rcases le_total (ctr 0).1 (ctr 1).1 with a | a <;>
  rcases le_total (ctr 1).1 (ctr 2).1 with b | b <;>
  rcases le_total (ctr 0).1 (ctr 2).1 with c | c
  all_goals first
    | exact go 0 1 2 (by decide) (by decide) a b
    | exact go 0 2 1 (by decide) (by decide) c b
    | exact go 2 0 1 (by decide) (by decide) c a
    | exact go 1 0 2 (by decide) (by decide) a c
    | exact go 1 2 0 (by decide) (by decide) b c
    | exact go 2 1 0 (by decide) (by decide) b a

/-! ## Orientation and the exact region -/

lemma coord_swap (c p : ℝ × ℝ) (θ : ℝ) :
    coord (c.2, c.1) (Real.pi / 2 - θ) p =
      ((coord c θ (p.2, p.1)).1, -(coord c θ (p.2, p.1)).2) := by
  simp only [coord, Real.cos_pi_div_two_sub, Real.sin_pi_div_two_sub]; ext <;> simp only <;> ring

lemma mem_sq_swap (c p : ℝ × ℝ) (θ : ℝ) :
    p ∈ sq (c.2, c.1) (Real.pi / 2 - θ) 1 ↔ (p.2, p.1) ∈ sq c θ 1 := by
  simp only [sq, Set.mem_ofPred_eq, coord_swap, abs_neg]

lemma mem_sqInt_swap (c p : ℝ × ℝ) (θ : ℝ) :
    p ∈ sqInt (c.2, c.1) (Real.pi / 2 - θ) 1 ↔ (p.2, p.1) ∈ sqInt c θ 1 := by
  simp only [sqInt, Set.mem_ofPred_eq, coord_swap, abs_neg]

/-- Packings are symmetric in the two sides. -/
theorem RectPacks.swap {n : ℕ} {w h : ℝ} (hp : RectPacks n w h) : RectPacks n h w := by
  obtain ⟨ctr, ang, hin, hd⟩ := hp
  refine ⟨fun i => ((ctr i).2, (ctr i).1), fun i => Real.pi / 2 - ang i, fun i p hp => ?_,
    fun i j hij => ?_⟩
  · have := hin i ((mem_sq_swap _ _ _).1 hp)
    simp only [rbox, Set.mem_ofPred_eq] at this ⊢
    exact ⟨this.2.2.1, this.2.2.2, this.1, this.2.1⟩
  · rw [Set.disjoint_left]
    intro p h1 h2
    exact Set.disjoint_left.1 (hd i j hij) ((mem_sqInt_swap _ _ _).1 h1) ((mem_sqInt_swap _ _ _).1 h2)

/-- Three upright squares from their centres, when the centres are pairwise `1` apart in some axis. -/
lemma upright3 {w h : ℝ} (c : Fin 3 → ℝ × ℝ) (hin : ∀ i, sq (c i) 0 1 ⊆ rbox w h)
    (hsep : ∀ i j, i ≠ j → 1 ≤ |(c i).1 - (c j).1| ∨ 1 ≤ |(c i).2 - (c j).2|) :
    RectPacks 3 w h := by
  refine ⟨c, fun _ => 0, hin, fun i j hij => ?_⟩
  rw [Set.disjoint_left]
  intro p hp hq
  rw [sqInt_zero] at hp hq
  rcases hsep i j hij with h | h
  · have := abs_sub_le ((c i).1) p.1 ((c j).1)
    rw [abs_sub_comm ((c i).1) p.1] at this; linarith [hp.1, hq.1]
  · have := abs_sub_le ((c i).2) p.2 ((c j).2)
    rw [abs_sub_comm ((c i).2) p.2] at this; linarith [hp.2, hq.2]

/-- Three in a row. -/
theorem rect3_row {w h : ℝ} (hw : 3 ≤ w) (hh : 1 ≤ h) : RectPacks 3 w h := by
  apply upright3 ![(1 / 2, 1 / 2), (3 / 2, 1 / 2), (5 / 2, 1 / 2)]
  · intro i; fin_cases i <;> exact sq_zero_subset (by norm_num) (by simp; linarith)
      (by norm_num) (by simp; linarith)
  · intro i j hij; left
    fin_cases i <;> fin_cases j <;> simp at hij ⊢ <;> norm_num [abs_of_pos, abs_of_neg]

/-- An L of three. -/
theorem rect3_ell {w h : ℝ} (hw : 2 ≤ w) (hh : 2 ≤ h) : RectPacks 3 w h := by
  apply upright3 ![(1 / 2, 1 / 2), (3 / 2, 1 / 2), (1 / 2, 3 / 2)]
  · intro i; fin_cases i <;> exact sq_zero_subset (by norm_num) (by simp; linarith)
      (by norm_num) (by simp; linarith)
  · intro i j hij
    fin_cases i <;> fin_cases j <;> simp at hij ⊢ <;> norm_num [abs_of_pos, abs_of_neg]

/-- **`R₃` exactly.**  Three unit squares (any angles) fit in `[0,w] × [0,h]` iff both sides are at
least `1` and: both are at least `2`, or one is at least `3`, or the short side `s < 2` and the long
side `ℓ` satisfy `ℓ + 2s ≥ 4 + 2√2` (the wedge). -/
theorem rect3_iff {w h : ℝ} : RectPacks 3 w h ↔ 1 ≤ w ∧ 1 ≤ h ∧
    ((2 ≤ w ∧ 2 ≤ h) ∨ 3 ≤ w ∨ 3 ≤ h ∨ (h < 2 ∧ 4 + 2 * Real.sqrt 2 ≤ w + 2 * h) ∨
      (w < 2 ∧ 4 + 2 * Real.sqrt 2 ≤ h + 2 * w)) := by
  have hs : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs1 : 1 < Real.sqrt 2 := by nlinarith [Real.sqrt_nonneg 2]
  constructor
  · intro hp
    obtain ⟨hw, hh⟩ := one_le_sides (hp.choose_spec.choose_spec.1 0)
    refine ⟨hw, hh, ?_⟩
    by_cases hw2 : 2 ≤ w <;> by_cases hh2 : 2 ≤ h
    · exact Or.inl ⟨hw2, hh2⟩
    · push Not at hh2
      rcases rect3_lower hp hh2 with h3 | h3
      · exact Or.inr (Or.inl h3)
      · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨hh2, h3⟩)))
    · push Not at hw2
      rcases rect3_lower hp.swap hw2 with h3 | h3
      · exact Or.inr (Or.inr (Or.inl h3))
      · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨hw2, h3⟩)))
    · push Not at hh2
      rcases rect3_lower hp hh2 with h3 | h3
      · exact Or.inr (Or.inl h3)
      · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨hh2, h3⟩)))
  · rintro ⟨hw, hh, h2 | h3 | h3 | ⟨hh2, hwe⟩ | ⟨hw2, hwe⟩⟩
    · exact rect3_ell h2.1 h2.2
    · exact rect3_row h3 hh
    · exact (rect3_row h3 hw).swap
    · rcases le_or_gt 3 w with h3 | h3
      · exact rect3_row h3 hh
      · exact rect3_wedge (by nlinarith) h3.le (by linarith)
    · rcases le_or_gt 3 h with h3 | h3
      · exact (rect3_row h3 hw).swap
      · exact (rect3_wedge (by nlinarith) h3.le (by linarith)).swap

end SquarePacking
