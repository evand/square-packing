import Sqpack.Rect5Pts
import Sqpack.Rect4
import Sqpack.Rect5Pack

/-!
# Five unit squares, short side `≥ 2`: the rest of the boundary (cover C)

`rect5_wide_lower`: if `2 ≤ a, b < 3` and `a + b < 4 + √2`, five unit squares do not fit in `[0,a] × [0,b]`
(Friedman's pieces 4.3–4.5 of `S₄`: `b = 3` for `2 ≤ a ≤ 1 + √2`, `a + b = 4 + √2` up to `2 + 1/√2`).

Reduction: widen the box to `a, b ≥ 1 + √2` (still `a + b < 4 + √2`), then scale it up onto the line
`a + b = 4 + √2`.  On the line, for `a ≥ 4 - √2` the four points of `rect5_cover` work; for `1 + √2 ≤ a ≤ 4 - √2`
we use **cover C**: ten points, weights `2/3` and `1/3`, total `14/3 < 5`
* `(c₀, 1)` and its three mirror images, weight `2/3` (`c₀ = 91/100`);
* `(1, 1)` and its mirror images (the four points of `rect5_cover`), weight `1/3`;
* `(1, b/2)`, `(a - 1, b/2)` on the midline, weight `1/3`.
Every closed unit square contains a **block** of weight exactly `1`: a row pair `{(c₀, 1), (1, 1)}`, a wall block
`{(c₀, 1 or b - 1), (1, b/2)}`, or a centre block `{(1, 1 or b - 1), (1, b/2), (a - 1, b/2)}` (up to mirrors).
Numerics and paper proofs: `search/rectfront/r5pts/R5PTS.md` §9–10.
-/

open Finset
open scoped Classical

namespace SquarePacking

/-! ## Weighted point counting with indexed points -/

/-- **Pigeonhole with weights**, for an indexed family of points (repetitions allowed): if every closed unit
square in `C` catches weight `≥ 1`, then `n` squares of side `L > 1` in `C` with disjoint interiors have
`n ≤` total weight. -/
theorem packing_le_idx {N : ℕ} (pt : Fin N → ℝ × ℝ) (w : Fin N → ℝ) (hw : ∀ j, 0 ≤ w j)
    (C : Set (ℝ × ℝ))
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C → 1 ≤ ∑ j ∈ univ.filter (fun j => pt j ∈ sq c θ 1), w j)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) ≤ ∑ j, w j := by
  have hIntSub : ∀ (c : ℝ × ℝ) (θ : ℝ), sqInt c θ L ⊆ sq c θ L := by
    rintro c θ p ⟨h1, h2⟩; exact ⟨le_of_lt h1, le_of_lt h2⟩
  set S : Fin n → Set (ℝ × ℝ) := fun i => sq (ctr i) (ang i) 1 with hSdef
  have hU : ∀ i, S i ⊆ C := fun i =>
    subset_trans (subset_trans (unit_subset_interior hL) (hIntSub _ _)) (hin i)
  have hSdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j) := fun i j hij =>
    Set.disjoint_of_subset (unit_subset_interior hL) (unit_subset_interior hL) (hdisj i j hij)
  have h1 : (n : ℝ) ≤ ∑ i : Fin n, ∑ j ∈ univ.filter (fun j => pt j ∈ S i), w j := by
    calc (n : ℝ) = ∑ _i : Fin n, (1 : ℝ) := by simp
      _ ≤ _ := Finset.sum_le_sum fun i _ => hcover (ctr i) (ang i) (hU i)
  have h2 : ∑ i : Fin n, ∑ j ∈ univ.filter (fun j => pt j ∈ S i), w j
      = ∑ j, ∑ i : Fin n, (if pt j ∈ S i then w j else 0) := by
    simp only [Finset.sum_filter]; exact Finset.sum_comm
  have h3 : ∀ j, ∑ i : Fin n, (if pt j ∈ S i then w j else 0) ≤ w j := by
    intro j
    have hcard : (((univ : Finset (Fin n)).filter (fun i => pt j ∈ S i)).card : ℝ) ≤ 1 := by
      exact_mod_cast card_filter_le_one S hSdisj (pt j)
    have hEq : ∑ i : Fin n, (if pt j ∈ S i then w j else 0)
        = (((univ : Finset (Fin n)).filter (fun i => pt j ∈ S i)).card : ℝ) * w j := by
      rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
    rw [hEq]; nlinarith [hw j]
  linarith [Finset.sum_le_sum fun j (_ : j ∈ univ) => h3 j]

/-- A block `J` of points in `X` with weight `≥ 1` gives capture `≥ 1`. -/
lemma cap_of_block {N : ℕ} {pt : Fin N → ℝ × ℝ} {w : Fin N → ℝ} (hw : ∀ j, 0 ≤ w j) {X : Set (ℝ × ℝ)}
    (J : Finset (Fin N)) (hJ : ∀ j ∈ J, pt j ∈ X) (h1 : 1 ≤ ∑ j ∈ J, w j) :
    1 ≤ ∑ j ∈ univ.filter (fun j => pt j ∈ X), w j :=
  h1.trans (Finset.sum_le_sum_of_subset_of_nonneg (fun j hj => Finset.mem_filter.2 ⟨mem_univ _, hJ j hj⟩)
    (fun j _ _ => hw j))

/-- **Scaling.**  A weighted cover of total `< 5` for the box `[0, l a] × [0, l b]`, `l > 1`, rules out five
unit squares in `[0,a] × [0,b]`. -/
theorem no_pack_of_cover {N : ℕ} (pt : Fin N → ℝ × ℝ) (w : Fin N → ℝ) (hw : ∀ j, 0 ≤ w j)
    {a b l : ℝ} (hl : 1 < l)
    (hcov : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ rbox (l * a) (l * b) →
      1 ≤ ∑ j ∈ univ.filter (fun j => pt j ∈ sq c θ 1), w j)
    (hW : ∑ j, w j < 5) : ¬ RectPacks 5 a b := by
  rintro ⟨ctr, ang, hin, hd⟩
  have l0 : 0 < l := by linarith
  have key := packing_le_idx pt w hw (rbox (l * a) (l * b)) hcov 5 l hl
    (fun i => (l * (ctr i).1, l * (ctr i).2)) ang
    (fun i p hp => by
      have := hin i ((mem_sq_scale l0).1 hp)
      simp only [rbox, Set.mem_ofPred_eq] at this ⊢
      obtain ⟨e1, e2, e3, e4⟩ := this
      refine ⟨?_, ?_, ?_, ?_⟩
      · have := mul_nonneg l0.le e1; field_simp at this; linarith
      · have := mul_le_mul_of_nonneg_left e2 l0.le; field_simp at this; linarith
      · have := mul_nonneg l0.le e3; field_simp at this; linarith
      · have := mul_le_mul_of_nonneg_left e4 l0.le; field_simp at this; linarith)
    (fun i j hij => by
      rw [Set.disjoint_left]; intro p h1 h2
      exact Set.disjoint_left.1 (hd i j hij) ((mem_sqInt_scale l0).1 h1) ((mem_sqInt_scale l0).1 h2))
  norm_num at key; linarith

/-! ## Small facts -/

lemma sqrt2_bounds : 1.4142 < Real.sqrt 2 ∧ Real.sqrt 2 < 1.4143 := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  constructor <;> nlinarith [Real.sqrt_nonneg 2]

lemma frame_more {c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) :
    2 * (c * s) = (c + s) ^ 2 - 1 ∧ (c + s) ^ 2 ≤ 2 ∧ 2 * (c * s) ≤ c + s := by
  refine ⟨by linear_combination -hcs, by nlinarith [sq_nonneg (c - s)], ?_⟩
  nlinarith [sq_nonneg (c - s), mul_nonneg hc hs]

lemma mp_E1 {c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) :
    0 ≤ c * (c - s) / 2 - (19 / 10 - Real.sqrt 2) * c + (2 - Real.sqrt 2) * s := by
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  nlinarith [mul_nonneg hc hs, sq_nonneg (c - 1), sq_nonneg (s - 1 / 4),
    mul_nonneg hs (by nlinarith : (0:ℝ) ≤ 1 - c)]

lemma mp_E2 {c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) :
    0 ≤ c * (c - s) / 2 - (Real.sqrt 2 - 11 / 10) * c + s / 2 := by
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  nlinarith [mul_nonneg hc hs, sq_nonneg (c - 1 / 3), mul_nonneg hs (by nlinarith : (0:ℝ) ≤ 1 - c)]

/-! ## The row lemma -/

/-- `0.39 c + (1 - k) s ≤ 1/2` on the quarter circle (tight at `θ = 90°`). -/
lemma row_key {c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) :
    39 / 100 * s + (1 - (c + s) / 2) * c ≤ 1 / 2 := by
  nlinarith [sq_nonneg (c - 48 / 100), sq_nonneg (s - 88 / 100), mul_nonneg hc hs, sq_nonneg (c + s - 1),
    mul_nonneg (mul_nonneg hc hs) hs, mul_nonneg (mul_nonneg hc hs) hc]

/-- **Row lemma.**  A point within `0.39` of the centre horizontally and `1 - k` vertically is in the square. -/
lemma In_row {dx dy c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hx : |dx| ≤ 39 / 100) (hy : |dy| ≤ 1 - (c + s) / 2) : In dx dy c s := by
  have k1 := row_key hc hs hcs
  have k2 := row_key hs hc (by linarith)
  have b1 := mul_le_mul_of_nonneg_right hx hc; have b2 := mul_le_mul_of_nonneg_right hy hs
  have b3 := mul_le_mul_of_nonneg_right hx hs; have b4 := mul_le_mul_of_nonneg_right hy hc
  constructor
  · calc |dx * c + dy * s| ≤ |dx * c| + |dy * s| := abs_add_le _ _
      _ = |dx| * c + |dy| * s := by rw [abs_mul, abs_mul, abs_of_nonneg hc, abs_of_nonneg hs]
      _ ≤ 1 / 2 := by linarith
  · calc |-dx * s + dy * c| ≤ |-dx * s| + |dy * c| := abs_add_le _ _
      _ = |dx| * s + |dy| * c := by rw [abs_mul, abs_mul, abs_neg, abs_of_nonneg hc, abs_of_nonneg hs]
      _ ≤ 1 / 2 := by linarith

/-! ## Side lemmas with the centre past the point row -/

/-- Core for a centre at depth `β ≥ 0` beyond the point row (offsets `(A, -β)`, `(-A', -β)`), `s ≤ c`: it
suffices that `c s X + β ≤ (c + s)/2`. -/
lemma neg_core {A A' β c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) (hsc : s ≤ c)
    (hA : 0 ≤ A) (hA' : 0 ≤ A') (hβ : 0 ≤ β) (hβ2 : β ≤ 1 / 2) (hX : A + A' ≤ 1)
    (hK : c * s * (A + A') + β ≤ (c + s) / 2) : In A (-β) c s ∨ In (-A') (-β) c s := by
  obtain ⟨p1, p2, c1, s1⟩ := frame_facts hc hs hcs
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  by_contra hno
  rw [not_or] at hno
  obtain ⟨hQ, hP⟩ := hno
  have hc0 : 0 < c := by nlinarith
  have f1 := mul_nonneg hA hc; have f2 := mul_nonneg hA hs; have f3 := mul_nonneg hA' hc
  have f4 := mul_nonneg hA' hs; have f5 := mul_nonneg hβ hc; have f6 := mul_nonneg hβ hs
  have g1 : β * c ≤ 1 / 2 := by nlinarith [mul_le_mul_of_nonneg_left c1 hβ]
  have g2 : β * s ≤ 1 / 2 := by nlinarith [mul_le_mul_of_nonneg_left s1 hβ]
  have gX : (A + A') * c ≤ 1 := by nlinarith [mul_le_mul_of_nonneg_left c1 (by linarith : 0 ≤ A + A')]
  have gXs : (A + A') * s ≤ 1 := by nlinarith [mul_le_mul_of_nonneg_left s1 (by linarith : 0 ≤ A + A')]
  have e1 : β * c ^ 2 + β * s ^ 2 = β := by linear_combination β * hcs
  have gcs : c * s * (A + A') ≤ c * s := by nlinarith [mul_nonneg hc hs]
  have gp : 2 * (c * s) ≤ c + s := by nlinarith
  have wQ : 1 / 2 < A * c + -β * s ∨ -A * s + -β * c < -(1 / 2) := by
    rcases not_In hQ with h | h | h | h
    · exact Or.inl h
    · linarith
    · linarith
    · exact Or.inr h
  have wP : -A' * c + -β * s < -(1 / 2) ∨ 1 / 2 < -(-A') * s + -β * c := by
    rcases not_In hP with h | h | h | h
    · linarith
    · exact Or.inl h
    · exact Or.inr h
    · linarith
  rcases wQ with q | q <;> rcases wP with r | r
  · linarith
  · -- weights `(s, c)`: `c s X - β > (c + s)/2`
    have k1 := mul_le_mul_of_nonneg_left q.le hs
    have k2 := mul_lt_mul_of_pos_left r hc0
    nlinarith
  · -- weights `(c, s)`: `c s X + β > (c + s)/2`
    have k1 := mul_lt_mul_of_pos_left q hc0
    have k2 := mul_le_mul_of_nonneg_left r.le hs
    nlinarith
  · linarith

/-- The case `B ≥ 0` of `In_side2_le`. -/
lemma side2_pos {A A' B y₀ c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) (hsc : s ≤ c)
    (hA : 0 ≤ A) (hA' : 0 ≤ A') (hX : A + A' ≤ 1) (hy₀ : y₀ ≤ 1) (hXy : (A + A') / 2 + y₀ ≤ Real.sqrt 2)
    (hB : B ≤ y₀ - (c + s) / 2) (hB0 : 0 ≤ B) : In A B c s ∨ In (-A') B c s := by
  obtain ⟨p1, p2, c1, s1⟩ := frame_facts hc hs hcs
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2g : 1 < Real.sqrt 2 := by nlinarith [Real.sqrt_nonneg 2]
  by_contra hno
  rw [not_or] at hno
  obtain ⟨hQ, hP⟩ := hno
  have hc0 : 0 < c := by
    by_contra h; push Not at h
    have : c = 0 := le_antisymm h hc
    subst this; nlinarith
  have f1 := mul_nonneg hA hc; have f2 := mul_nonneg hA hs; have f3 := mul_nonneg hA' hc
  have f4 := mul_nonneg hA' hs; have f5 := mul_nonneg hB0 hc; have f6 := mul_nonneg hB0 hs
  have g1 : A * s ≤ A * c := mul_le_mul_of_nonneg_left hsc hA
  have g2 : B * c ≤ B := by nlinarith [mul_le_mul_of_nonneg_left c1 hB0]
  have g3 : B * s ≤ B := by nlinarith [mul_le_mul_of_nonneg_left s1 hB0]
  have g4 : (A + A') * c ≤ A + A' := by nlinarith [mul_le_mul_of_nonneg_left c1 (by linarith : 0 ≤ A + A')]
  have hu : 1 / 2 < A * c + B * s := by
    rcases not_In hQ with h | h | h | h
    · exact h
    · linarith
    · linarith
    · linarith
  rcases not_In hP with h | h | h | h
  · linarith
  · linarith
  · -- `c s X + B > (c + s)/2`, against the chord of `X (p² - 1)/2 + y₀ - p`
    have k1 := mul_le_mul_of_nonneg_left hu.le hs
    have k2 := mul_lt_mul_of_pos_left h hc0
    have e1 : B * c ^ 2 + B * s ^ 2 = B := by linear_combination B * hcs
    have e2 : (A + A') * c ^ 2 + (A + A') * s ^ 2 = A + A' := by linear_combination (A + A') * hcs
    have m1 := mul_nonneg (mul_nonneg (by linarith : 0 ≤ A + A') (by linarith : 0 ≤ c + s - 1))
      (by linarith : 0 ≤ Real.sqrt 2 - (c + s))
    have m2 := mul_nonneg (by linarith : 0 ≤ Real.sqrt 2 - (c + s)) (by linarith : 0 ≤ 1 - y₀)
    have m3 := mul_nonneg (by linarith : 0 ≤ c + s - 1) (by linarith : 0 ≤ Real.sqrt 2 - (A + A') / 2 - y₀)
    have hL : (Real.sqrt 2 - 1) * ((A + A') * (Real.sqrt 2 + 1) * (c + s - 1) / 2 + y₀ - (c + s)) =
        (Real.sqrt 2 - (c + s)) * (y₀ - 1) + (c + s - 1) * ((A + A') / 2 + y₀ - Real.sqrt 2) := by
      linear_combination ((A + A') * (c + s - 1) / 2) * r2
    have hLpos : 0 < (A + A') * (Real.sqrt 2 + 1) * (c + s - 1) / 2 + y₀ - (c + s) := by linarith
    have := mul_pos (by linarith : 0 < Real.sqrt 2 - 1) hLpos
    linarith
  · linarith

/-- **Side lemma, general distance** (DS7 Lemma 2 with the point row at distance `y₀` from the wall; the centre
may be up to `(√2 - 1)/2` past the row). -/
lemma In_side2_le {A A' B y₀ c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1) (hsc : s ≤ c)
    (hA : 0 ≤ A) (hA' : 0 ≤ A') (hX : A + A' ≤ 1) (hy₀ : y₀ ≤ 1) (hXy : (A + A') / 2 + y₀ ≤ Real.sqrt 2)
    (hB : B ≤ y₀ - (c + s) / 2) (hB' : -((Real.sqrt 2 - 1) / 2) ≤ B) : In A B c s ∨ In (-A') B c s := by
  rcases le_total 0 B with hB0 | hB0
  · exact side2_pos hc hs hcs hsc hA hA' hX hy₀ hXy hB hB0
  · obtain ⟨p1, p2, c1, s1⟩ := frame_facts hc hs hcs
    obtain ⟨e1, -, -⟩ := frame_more hc hs hcs
    obtain ⟨r2g, r2l⟩ := sqrt2_bounds
    have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    have h1 := mul_le_mul_of_nonneg_left hX (mul_nonneg hc hs)
    -- `(c + s - √2)(c + s + √2 - 1) ≤ 0`
    have hm := mul_nonneg (by linarith : 0 ≤ Real.sqrt 2 - (c + s)) (by linarith : 0 ≤ c + s + Real.sqrt 2 - 1)
    have := neg_core hc hs hcs hsc hA hA' (β := -B) (by linarith) (by linarith) hX (by nlinarith)
    rwa [neg_neg] at this

lemma In_side2 {A A' B y₀ c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hA : 0 ≤ A) (hA' : 0 ≤ A') (hX : A + A' ≤ 1) (hy₀ : y₀ ≤ 1) (hXy : (A + A') / 2 + y₀ ≤ Real.sqrt 2)
    (hB : B ≤ y₀ - (c + s) / 2) (hB' : -((Real.sqrt 2 - 1) / 2) ≤ B) : In A B c s ∨ In (-A') B c s := by
  rcases le_total s c with h | h
  · exact In_side2_le hc hs hcs h hA hA' hX hy₀ hXy hB hB'
  · have n : ∀ x y : ℝ, In (-x) y s c → In x y c s := fun x y h => In_negx.1 h
    have n' : ∀ x y : ℝ, In x y s c → In (-x) y c s := fun x y h => In_negx.1 (by rwa [neg_neg])
    rcases In_side2_le hs hc (by linarith) h hA' hA (by linarith) hy₀ (by linarith) (by linarith) hB'
      with h1 | h1
    · exact Or.inr (n' _ _ h1)
    · exact Or.inl (n _ _ h1)

/-- **Far side** (no wall): offsets `(A, -β)`, `(-A', -β)` with `√2 - 1 ≤ X ≤ 1`, `X + 2β ≤ √2`; tight for the
`45°` square at the centre of the box on the line. -/
lemma In_far {A A' β c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hA : 0 ≤ A) (hA' : 0 ≤ A') (hβ : 0 ≤ β) (hX1 : Real.sqrt 2 - 1 ≤ A + A') (hX : A + A' ≤ 1)
    (hXβ : A + A' + 2 * β ≤ Real.sqrt 2) : In A (-β) c s ∨ In (-A') (-β) c s := by
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1415 / 1000 := by nlinarith [Real.sqrt_nonneg 2]
  have key : ∀ {c s : ℝ}, 0 ≤ c → 0 ≤ s → c ^ 2 + s ^ 2 = 1 →
      c * s * (A + A') + β ≤ (c + s) / 2 := by
    intro c s hc hs hcs
    obtain ⟨p1, p2, -, -⟩ := frame_facts hc hs hcs
    have f := mul_nonneg (by linarith : 0 ≤ Real.sqrt 2 - (c + s))
      (by nlinarith : 0 ≤ (A + A') * (Real.sqrt 2 + (c + s)) - 1)
    have e : 2 * (c * s) = (c + s) ^ 2 - 1 := by linear_combination -hcs
    nlinarith
  rcases le_total s c with h | h
  · exact neg_core hc hs hcs h hA hA' hβ (by linarith) hX (key hc hs hcs)
  · have n : ∀ x y : ℝ, In (-x) y s c → In x y c s := fun x y h => In_negx.1 h
    have n' : ∀ x y : ℝ, In x y s c → In (-x) y c s := fun x y h => In_negx.1 (by rwa [neg_neg])
    rcases neg_core hs hc (by linarith) h hA' hA hβ (by linarith) (by linarith)
      (by have := key hs hc (by linarith); linarith) with h1 | h1
    · exact Or.inr (n' _ _ h1)
    · exact Or.inl (n _ _ h1)

/-- **The far midline point** `m′ = (a - 1, b/2)` from a centre right of `x = 1.1` (cell 3). -/
lemma In_mp {μ ν a b c s : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hab : a + b = 4 + Real.sqrt 2) (hb1 : 2 * Real.sqrt 2 ≤ b) (hb2 : b ≤ 3)
    (hμ0 : 0 ≤ μ) (hμ : μ ≤ a - 21 / 10) (hν0 : 0 ≤ ν) (hν : ν ≤ b / 2 - 2 + (c + s) / 2) : In μ ν c s := by
  obtain ⟨p1, p2, c1, s1⟩ := frame_facts hc hs hcs
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  have E1 := mp_E1 hc hs hcs
  have E2 := mp_E2 hc hs hcs
  -- `H` is linear in `b/2 ∈ [√2, 3/2]`
  have H : 0 ≤ c * (c - s) / 2 - (a - 21 / 10) * c + (2 - b / 2) * s := by
    have q1 := mul_nonneg (by linarith : 0 ≤ 3 / 2 - b / 2) E1
    have q2 := mul_nonneg (by linarith : 0 ≤ b / 2 - Real.sqrt 2) E2
    have hpos : 0 < 3 / 2 - Real.sqrt 2 := by linarith
    have e : (3 / 2 - Real.sqrt 2) * (c * (c - s) / 2 - (a - 21 / 10) * c + (2 - b / 2) * s) =
        (3 / 2 - b / 2) * (c * (c - s) / 2 - (19 / 10 - Real.sqrt 2) * c + (2 - Real.sqrt 2) * s) +
        (b / 2 - Real.sqrt 2) * (c * (c - s) / 2 - (Real.sqrt 2 - 11 / 10) * c + s / 2) := by
      have ha : a = 4 + Real.sqrt 2 - b := by linarith
      rw [ha]; ring
    by_contra hneg; push Not at hneg
    have := mul_neg_of_pos_of_neg hpos hneg
    linarith
  have hμ2 : μ ≤ 1 / 2 := by linarith
  have hν2 : ν ≤ 1 / 2 := by linarith
  have u1 := mul_le_mul_of_nonneg_right hμ hc; have u2 := mul_le_mul_of_nonneg_right hν hs
  have u3 := mul_nonneg hμ0 hc; have u4 := mul_nonneg hν0 hs
  have v1 := mul_le_mul_of_nonneg_right hμ2 hs; have v2 := mul_le_mul_of_nonneg_right hν2 hc
  have v3 := mul_nonneg hμ0 hs; have v4 := mul_nonneg hν0 hc
  constructor
  · rw [abs_le]; constructor
    · linarith
    · linarith
  · rw [abs_le]; constructor
    · linarith
    · linarith

/-! ## Cover C -/

/-- The ten points of cover C. -/
noncomputable def ptC (a b : ℝ) : Fin 10 → ℝ × ℝ :=
  ![(91 / 100, 1), (a - 91 / 100, 1), (91 / 100, b - 1), (a - 91 / 100, b - 1),
    (1, 1), (a - 1, 1), (1, b - 1), (a - 1, b - 1), (1, b / 2), (a - 1, b / 2)]

/-- Their weights. -/
noncomputable def wC : Fin 10 → ℝ := ![2 / 3, 2 / 3, 2 / 3, 2 / 3, 1 / 3, 1 / 3, 1 / 3, 1 / 3, 1 / 3, 1 / 3]

lemma wC_nonneg (j : Fin 10) : 0 ≤ wC j := by fin_cases j <;> simp [wC] <;> norm_num

lemma wC_sum : ∑ j, wC j = 14 / 3 := by simp [Fin.sum_univ_succ, wC]; norm_num

/-- The line `a + b = 4 + √2` with `1 + √2 ≤ a ≤ 4 - √2`. -/
def OnLine (a b : ℝ) : Prop := 1 + Real.sqrt 2 ≤ a ∧ a ≤ 4 - Real.sqrt 2 ∧ a + b = 4 + Real.sqrt 2

/-- A square with frame `(c, s)` and centre `(x, y)` in the lower-left quarter of the box. -/
def QCell (a b c s x y : ℝ) : Prop :=
  0 ≤ c ∧ 0 ≤ s ∧ c ^ 2 + s ^ 2 = 1 ∧ (c + s) / 2 ≤ x ∧ x ≤ a / 2 ∧ (c + s) / 2 ≤ y ∧ y ≤ b / 2

/-- The blocks a lower-left square can contain (offsets from its centre). -/
def QB (a b c s x y : ℝ) : Prop :=
  (In (91 / 100 - x) (1 - y) c s ∧ In (1 - x) (1 - y) c s) ∨
  (In (91 / 100 - x) (1 - y) c s ∧ In (1 - x) (b / 2 - y) c s) ∨
  (In (91 / 100 - x) (b - 1 - y) c s ∧ In (1 - x) (b / 2 - y) c s) ∨
  (In (1 - x) (1 - y) c s ∧ In (1 - x) (b / 2 - y) c s ∧ In (a - 1 - x) (b / 2 - y) c s) ∨
  (In (1 - x) (b - 1 - y) c s ∧ In (1 - x) (b / 2 - y) c s ∧ In (a - 1 - x) (b / 2 - y) c s)

/-- **Every lower-left square contains a block.**  Three cells: the row (`y ≤ 2 - k`: `(c₀, 1)`, `(1, 1)`),
the wall band (`x ≤ 1.1`: `(1, b/2)` and `(c₀, 1)` or `(c₀, b - 1)`), the centre (`x ≥ 1.1`: both midline points
and `(1, 1)` or `(1, b - 1)`). -/
theorem quadC {a b c s x y : ℝ} (hL : OnLine a b) (hq : QCell a b c s x y) : QB a b c s x y := by
  obtain ⟨ha1, ha2, hab⟩ := hL
  obtain ⟨hc, hs, hcs, hx1, hx2, hy1, hy2⟩ := hq
  have hsc : s ^ 2 + c ^ 2 = 1 := by linarith
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1415 / 1000 := by nlinarith [Real.sqrt_nonneg 2]
  have r2g : 1414 / 1000 < Real.sqrt 2 := by nlinarith [Real.sqrt_nonneg 2]
  obtain ⟨p1, p2, -, -⟩ := frame_facts hc hs hcs
  have ab : ∀ {z : ℝ}, -(1 - (c + s) / 2) ≤ z → z ≤ 1 - (c + s) / 2 → |z| ≤ 1 - (c + s) / 2 :=
    fun h h' => abs_le.2 ⟨h, h'⟩
  rcases le_total y (2 - (c + s) / 2) with y1 | y1
  · -- the row
    refine Or.inl ⟨?_, In_corner hc hs hcs (ab (by linarith) (by linarith)) (ab (by linarith) (by linarith))⟩
    rcases le_total x 1 with x1 | x1
    · exact In_corner hc hs hcs (ab (by linarith) (by linarith)) (ab (by linarith) (by linarith))
    · exact In_row hc hs hcs (abs_le.2 ⟨by linarith, by linarith⟩) (ab (by linarith) (by linarith))
  · have hm : In (1 - x) (b / 2 - y) c s := by
      rcases le_total x 1 with x1 | x1 <;>
        exact In_corner hc hs hcs (ab (by linarith) (by linarith)) (ab (by linarith) (by linarith))
    rcases le_total x (11 / 10) with x1 | x1
    · -- the wall band
      rcases In_side2 hs hc hsc (A := b - 1 - y) (A' := y - 1) (B := 91 / 100 - x) (y₀ := 91 / 100)
          (by linarith) (by linarith) (by linarith) (by norm_num) (by linarith) (by linarith) (by linarith)
        with h | h
      · exact Or.inr (Or.inr (Or.inl ⟨In_swap.1 h, hm⟩))
      · exact Or.inr (Or.inl ⟨In_congr (In_swap.1 h) rfl (by ring), hm⟩)
    · -- the centre
      have hm' : In (a - 1 - x) (b / 2 - y) c s :=
        In_mp hc hs hcs hab (by linarith) (by linarith) (by linarith) (by linarith) (by linarith) (by linarith)
      rcases In_far hs hc hsc (A := b - 1 - y) (A' := y - 1) (β := x - 1) (by linarith) (by linarith)
          (by linarith) (by linarith) (by linarith) (by linarith) with h | h
      · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨In_congr (In_swap.1 h) (by ring) rfl, hm, hm'⟩)))
      · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨In_congr (In_swap.1 h) (by ring) (by ring), hm, hm'⟩)))

/-- Mirrors of offsets. -/
lemma In_mx {q dy x c s a : ℝ} (h : In (q - (a - x)) dy s c) : In (a - q - x) dy c s :=
  In_negx.1 (In_congr h (by ring) rfl)

lemma In_my {q dx y c s b : ℝ} (h : In dx (q - (b - y)) s c) : In dx (b - q - y) c s :=
  In_negx.1 (In_neg.1 (In_congr h (by ring) (by ring)))

lemma In_mxy {q q' x y c s a b : ℝ} (h : In (q - (a - x)) (q' - (b - y)) c s) : In (a - q - x) (b - q' - y) c s :=
  In_neg.1 (In_congr h (by ring) (by ring))

/-- **Cover C** catches weight `≥ 1` from every closed unit square in the box. -/
theorem coverC {a b : ℝ} (hL : OnLine a b) {m : ℝ × ℝ} {θ : ℝ}
    (hin : sq m θ 1 ⊆ rbox a b) : 1 ≤ ∑ j ∈ univ.filter (fun j => ptC a b j ∈ sq m θ 1), wC j := by
  obtain ⟨P, Q, hP, hQ, hPQ, h1, h2, h3, h4, hf⟩ := framed hin
  have hQP : Q ^ 2 + P ^ 2 = 1 := by linarith
  -- a block of indices, all of whose points satisfy `In`, gives capture `≥ 1`
  have B : ∀ J : Finset (Fin 10), (∀ j ∈ J, In ((ptC a b j).1 - m.1) ((ptC a b j).2 - m.2) P Q) →
      1 ≤ ∑ j ∈ J, wC j → 1 ≤ ∑ j ∈ univ.filter (fun j => ptC a b j ∈ sq m θ 1), wC j :=
    fun J hJ hw => cap_of_block wC_nonneg J (fun j hj => mem_of_In hf hPQ (hJ j hj)) hw
  set x := m.1; set y := m.2
  rcases le_total x (a / 2) with hx | hx <;> rcases le_total y (b / 2) with hy | hy
  · -- lower left
    rcases quadC hL ⟨hP, hQ, hPQ, h1, hx, h3, hy⟩ with ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2, i3⟩ |
      ⟨i1, i2, i3⟩
    · exact B {0, 4} (by simp [ptC]; exact ⟨i1, i2⟩) (by simp [wC]; norm_num)
    · exact B {0, 8} (by simp [ptC]; exact ⟨i1, i2⟩) (by simp [wC]; norm_num)
    · exact B {2, 8} (by simp [ptC]; exact ⟨i1, i2⟩) (by simp [wC]; norm_num)
    · exact B {4, 8, 9} (by simp [ptC]; exact ⟨i1, i2, i3⟩) (by simp [wC]; norm_num)
    · exact B {6, 8, 9} (by simp [ptC]; exact ⟨i1, i2, i3⟩) (by simp [wC]; norm_num)
  · -- upper left: mirror in `y`
    rcases quadC hL (a := a) (b := b) (c := Q) (s := P) (x := x) (y := b - y)
        ⟨hQ, hP, hQP, by linarith, hx, by linarith, by linarith⟩ with
      ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2, i3⟩ | ⟨i1, i2, i3⟩
    · exact B {2, 6} (by simp [ptC]; exact ⟨In_congr (In_my i1) (by ring) (by ring),
        In_congr (In_my i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {2, 8} (by simp [ptC]; exact ⟨In_congr (In_my i1) (by ring) (by ring),
        In_congr (In_my i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {0, 8} (by simp [ptC]; exact ⟨In_congr (In_my i1) (by ring) (by ring),
        In_congr (In_my i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {6, 8, 9} (by simp [ptC]; exact ⟨In_congr (In_my i1) (by ring) (by ring),
        In_congr (In_my i2) (by ring) (by ring), In_congr (In_my i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)
    · exact B {4, 8, 9} (by simp [ptC]; exact ⟨In_congr (In_my i1) (by ring) (by ring),
        In_congr (In_my i2) (by ring) (by ring), In_congr (In_my i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)
  · -- lower right: mirror in `x`
    rcases quadC hL (a := a) (b := b) (c := Q) (s := P) (x := a - x) (y := y)
        ⟨hQ, hP, hQP, by linarith, by linarith, by linarith, hy⟩ with
      ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2, i3⟩ | ⟨i1, i2, i3⟩
    · exact B {1, 5} (by simp [ptC]; exact ⟨In_congr (In_mx i1) (by ring) (by ring),
        In_congr (In_mx i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {1, 9} (by simp [ptC]; exact ⟨In_congr (In_mx i1) (by ring) (by ring),
        In_congr (In_mx i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {3, 9} (by simp [ptC]; exact ⟨In_congr (In_mx i1) (by ring) (by ring),
        In_congr (In_mx i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {5, 9, 8} (by simp [ptC]; exact ⟨In_congr (In_mx i1) (by ring) (by ring),
        In_congr (In_mx i2) (by ring) (by ring), In_congr (In_mx i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)
    · exact B {7, 9, 8} (by simp [ptC]; exact ⟨In_congr (In_mx i1) (by ring) (by ring),
        In_congr (In_mx i2) (by ring) (by ring), In_congr (In_mx i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)
  · -- upper right: the half turn
    rcases quadC hL (a := a) (b := b) (c := P) (s := Q) (x := a - x) (y := b - y)
        ⟨hP, hQ, hPQ, by linarith, by linarith, by linarith, by linarith⟩ with
      ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2⟩ | ⟨i1, i2, i3⟩ | ⟨i1, i2, i3⟩
    · exact B {3, 7} (by simp [ptC]; exact ⟨In_congr (In_mxy i1) (by ring) (by ring),
        In_congr (In_mxy i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {3, 9} (by simp [ptC]; exact ⟨In_congr (In_mxy i1) (by ring) (by ring),
        In_congr (In_mxy i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {1, 9} (by simp [ptC]; exact ⟨In_congr (In_mxy i1) (by ring) (by ring),
        In_congr (In_mxy i2) (by ring) (by ring)⟩) (by simp [wC]; norm_num)
    · exact B {7, 9, 8} (by simp [ptC]; exact ⟨In_congr (In_mxy i1) (by ring) (by ring),
        In_congr (In_mxy i2) (by ring) (by ring), In_congr (In_mxy i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)
    · exact B {5, 9, 8} (by simp [ptC]; exact ⟨In_congr (In_mxy i1) (by ring) (by ring),
        In_congr (In_mxy i2) (by ring) (by ring), In_congr (In_mxy i3) (by ring) (by ring)⟩)
        (by simp [wC]; norm_num)

/-! ## Assembly -/

/-- The four points of `rect5_cover`, indexed, weight `1`. -/
noncomputable def ptD (a b : ℝ) : Fin 4 → ℝ × ℝ := ![(1, 1), (a - 1, 1), (1, b - 1), (a - 1, b - 1)]

lemma coverD {a b : ℝ} (ha : 2 ≤ a) (ha' : a ≤ 2 * Real.sqrt 2) (hb' : b ≤ 2 * Real.sqrt 2)
    (hab : a + b ≤ 4 + Real.sqrt 2) {m : ℝ × ℝ} {θ : ℝ} (hin : sq m θ 1 ⊆ rbox a b) :
    1 ≤ ∑ j ∈ univ.filter (fun j => ptD a b j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j := by
  have B : ∀ j : Fin 4, ptD a b j ∈ sq m θ 1 →
      1 ≤ ∑ j ∈ univ.filter (fun j => ptD a b j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j :=
    fun j hj => cap_of_block (fun _ => zero_le_one) {j} (by simpa using hj) (by simp)
  rcases rect5_cover ha ha' hb' hab hin with h | h | h | h
  · exact B 0 (by simpa [ptD] using h)
  · exact B 1 (by simpa [ptD] using h)
  · exact B 2 (by simpa [ptD] using h)
  · exact B 3 (by simpa [ptD] using h)

lemma rect5_wide_aux {a b : ℝ} (ha : 2 ≤ a) (hb3 : b < 3) (hab' : a ≤ b)
    (hab : a + b < 4 + Real.sqrt 2) : ¬ RectPacks 5 a b := by
  intro hp
  have r2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2l : Real.sqrt 2 < 1415 / 1000 := by nlinarith [Real.sqrt_nonneg 2]
  have r2g : 1414 / 1000 < Real.sqrt 2 := by nlinarith [Real.sqrt_nonneg 2]
  -- widen to `a₁, b₁ ≥ 1 + √2`
  set a₁ := max a (1 + Real.sqrt 2) with ha₁
  set b₁ := max b (1 + Real.sqrt 2) with hb₁
  have hp₁ : RectPacks 5 a₁ b₁ :=
    (RectPacks.mono (RectPacks.mono hp (le_max_left a _)).swap (le_max_left b _)).swap
  have e1 : 1 + Real.sqrt 2 ≤ a₁ := le_max_right _ _
  have e2 : a₁ ≤ b₁ := max_le_max hab' le_rfl
  have e3 : a₁ + b₁ < 4 + Real.sqrt 2 := by
    rcases max_cases a (1 + Real.sqrt 2) with ⟨h1, -⟩ | ⟨h1, -⟩ <;>
      rcases max_cases b (1 + Real.sqrt 2) with ⟨h2, -⟩ | ⟨h2, -⟩ <;>
      rw [← ha₁] at h1 <;> rw [← hb₁] at h2 <;> rw [h1, h2] <;> linarith
  -- scale onto the line
  set l := (4 + Real.sqrt 2) / (a₁ + b₁) with hl
  have l1 : 1 < l := (one_lt_div (by linarith)).2 e3
  have hsum : l * a₁ + l * b₁ = 4 + Real.sqrt 2 := by
    rw [← mul_add, hl, div_mul_cancel₀ _ (by linarith : a₁ + b₁ ≠ 0)]
  have f1 : 1 + Real.sqrt 2 ≤ l * a₁ := by nlinarith
  have f2 : l * a₁ ≤ l * b₁ := mul_le_mul_of_nonneg_left e2 (by linarith)
  rcases le_total (l * a₁) (4 - Real.sqrt 2) with h | h
  · exact no_pack_of_cover (ptC (l * a₁) (l * b₁)) wC wC_nonneg l1
      (fun c θ hc => coverC ⟨f1, h, hsum⟩ hc) (by rw [wC_sum]; norm_num) hp₁
  · exact no_pack_of_cover (ptD (l * a₁) (l * b₁)) (fun _ => 1) (fun _ => zero_le_one) l1
      (fun c θ hc => coverD (by linarith) (by linarith) (by linarith) (by linarith) hc) (by simp; norm_num) hp₁

/-- **`R₅` with both sides `≥ 2`, lower bound** : five unit squares
need a side `≥ 3` or `a + b ≥ 4 + √2` (Friedman's pieces 4.3–4.5 of `S₄`). -/
theorem rect5_wide_lower {a b : ℝ} (ha : 2 ≤ a) (hb : 2 ≤ b) (ha3 : a < 3)
    (hb3 : b < 3) (hab : a + b < 4 + Real.sqrt 2) : ¬ RectPacks 5 a b := by
  rcases le_total a b with h | h
  · exact rect5_wide_aux ha hb3 h hab
  · exact fun hp => rect5_wide_aux hb ha3 h (by linarith) hp.swap

/-- **`R₅` with both sides `≥ 2`, exactly**: five unit squares fit in `[0,a] × [0,b]`, `a, b ≥ 2`, iff one side
is at least `3` or `a + b ≥ 4 + √2` (Friedman's pieces 4.3–4.5 of `S₄`, Math Magic June 2010; at `a = b`,
`s(5) = 2 + 1/√2`, Göbel 1979). -/
theorem rect5_wide_iff {a b : ℝ} (ha : 2 ≤ a) (hb : 2 ≤ b) :
    RectPacks 5 a b ↔ 3 ≤ a ∨ 3 ≤ b ∨ 4 + Real.sqrt 2 ≤ a + b := by
  constructor
  · intro hp
    by_contra h
    push Not at h
    exact rect5_wide_lower ha hb h.1 h.2.1 h.2.2 hp
  · rintro (h | h | h)
    · exact (rect5_grid hb h).swap
    · exact rect5_grid ha h
    · exact rect5_stretch ha hb h

end SquarePacking
