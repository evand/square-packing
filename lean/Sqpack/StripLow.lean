import Sqpack.Rect5Low

/-!
# Low strips, every `n`: `n - 1` points on the midline

`strip_low`: in a strip of height `1 ≤ H < 2√2 - 1 ≈ 1.828`, `n` unit squares need length `L ≥ n`.

The proof of `rect5_low` (`Rect5Low.lean`) never used `n = 5`: in `[0,n] × [0,H]`, `H ≤ 2√2 - 1`, every closed
unit square contains one of the midline points `(i, H/2)`, `i = 1, …, n - 1` (`stripN_cover`).  A centre within
`1` of an end wall gets the end point (`In_corner`); a centre between two points gets one of them (`In_side2`).
Widen the strip, scale it onto `L = n` (`no_pack_of_cover_n`), and `n - 1` points cannot hold `n` disjoint
squares.  The bound is sharp as `n → ∞`: the tilted pieces of the strip frontier start just above
`H = 2√2 - 1` and approach it (Friedman, Math Magic June 2010, for `n ≤ 10`).
-/

open Finset
open scoped Classical

namespace SquarePacking

/-- **Scaling, any `n`.**  A weighted cover of total `< n` for the box `[0, l a] × [0, l b]`, `l > 1`, rules out
`n` unit squares in `[0,a] × [0,b]`.  (`no_pack_of_cover` is the case `n = 5`.) -/
theorem no_pack_of_cover_n {N : ℕ} (n : ℕ) (pt : Fin N → ℝ × ℝ) (w : Fin N → ℝ) (hw : ∀ j, 0 ≤ w j)
    {a b l : ℝ} (hl : 1 < l)
    (hcov : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ rbox (l * a) (l * b) →
      1 ≤ ∑ j ∈ univ.filter (fun j => pt j ∈ sq c θ 1), w j)
    (hW : ∑ j, w j < n) : ¬ RectPacks n a b := by
  rintro ⟨ctr, ang, hin, hd⟩
  have l0 : 0 < l := by linarith
  have key := packing_le_idx pt w hw (rbox (l * a) (l * b)) hcov n l hl
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
  linarith

/-- The lower half of `[0,n] × [0,H]`: some midline point `(i, H/2)`, `1 ≤ i ≤ n - 1`, is in the square. -/
lemma stripN_cell {n : ℕ} (hn : 2 ≤ n) {H c s x y : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hH : H ≤ 2 * Real.sqrt 2 - 1) (hx1 : (c + s) / 2 ≤ x) (hx2 : x ≤ n - (c + s) / 2)
    (hy1 : (c + s) / 2 ≤ y) (hy2 : y ≤ H / 2) :
    ∃ i : ℕ, 1 ≤ i ∧ i + 1 ≤ n ∧ In ((i : ℝ) - x) (H / 2 - y) c s := by
  obtain ⟨p1, p2, -, -⟩ := frame_facts hc hs hcs
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  have hn' : (2 : ℝ) ≤ n := by exact_mod_cast hn
  have ab : ∀ {z : ℝ}, -(1 - (c + s) / 2) ≤ z → z ≤ 1 - (c + s) / 2 → |z| ≤ 1 - (c + s) / 2 :=
    fun h h' => abs_le.2 ⟨h, h'⟩
  have corner : ∀ i : ℕ, 1 ≤ i → i + 1 ≤ n → |(i : ℝ) - x| ≤ 1 - (c + s) / 2 →
      ∃ i : ℕ, 1 ≤ i ∧ i + 1 ≤ n ∧ In ((i : ℝ) - x) (H / 2 - y) c s :=
    fun i h1 h2 hx => ⟨i, h1, h2, In_corner hc hs hcs hx (ab (by linarith) (by linarith))⟩
  -- near an end wall: the end point
  rcases le_total x 1 with x1 | x1
  · exact corner 1 le_rfl hn (by push_cast; exact ab (by linarith) (by linarith))
  rcases le_or_gt ((n : ℝ) - 1) x with xn | xn
  · refine corner (n - 1) (by omega) (by omega) ?_
    rw [Nat.cast_sub (by omega)]; push_cast; exact ab (by linarith) (by linarith)
  -- between the points `i = ⌊x⌋₊` and `i + 1`
  set i := ⌊x⌋₊ with hi
  have hix : (i : ℝ) ≤ x := Nat.floor_le (by linarith)
  have hxi : x < i + 1 := Nat.lt_floor_add_one x
  have i1 : 1 ≤ i := Nat.le_floor (by exact_mod_cast x1)
  have i2 : i + 2 ≤ n := by
    have : (i : ℝ) + 1 < n := by linarith
    exact_mod_cast (show ((i + 1 : ℕ) : ℝ) < n by push_cast; linarith)
  rcases In_side2 hc hs hcs (A := i + 1 - x) (A' := x - i) (B := H / 2 - y) (y₀ := H / 2)
    (by linarith) (by linarith) (by linarith) (by linarith) (by linarith) (by linarith) (by linarith)
    with h | h
  · exact ⟨i + 1, by omega, by omega, In_congr h (by push_cast; ring) rfl⟩
  · exact ⟨i, i1, by omega, In_congr h (by ring) rfl⟩

/-- The `n - 1` midline points `(j + 1, H/2)`. -/
noncomputable def ptN (n : ℕ) (H : ℝ) : Fin (n - 1) → ℝ × ℝ := fun j => ((j : ℝ) + 1, H / 2)

/-- **The midline cover** of `[0,n] × [0,H]`, `H ≤ 2√2 - 1`. -/
theorem stripN_cover {n : ℕ} (hn : 2 ≤ n) {H : ℝ} (hH : H ≤ 2 * Real.sqrt 2 - 1) {m : ℝ × ℝ} {θ : ℝ}
    (hin : sq m θ 1 ⊆ rbox n H) :
    1 ≤ ∑ j ∈ univ.filter (fun j => ptN n H j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j := by
  obtain ⟨P, Q, hP, hQ, hPQ, h1, h2, h3, h4, hf⟩ := framed hin
  have hQP : Q ^ 2 + P ^ 2 = 1 := by linarith
  have B : ∀ i : ℕ, 1 ≤ i → i + 1 ≤ n → ((i : ℝ), H / 2) ∈ sq m θ 1 →
      1 ≤ ∑ j ∈ univ.filter (fun j => ptN n H j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j := by
    intro i hi1 hi2 hmem
    have hj : ptN n H ⟨i - 1, by omega⟩ = ((i : ℝ), H / 2) := by
      simp only [ptN, Prod.mk.injEq, and_true]; rw [Nat.cast_sub hi1]; ring
    exact cap_of_block (fun _ => zero_le_one) {⟨i - 1, by omega⟩} (by simpa [hj] using hmem) (by simp)
  rcases le_total m.2 (H / 2) with hy | hy
  · obtain ⟨i, hi1, hi2, h⟩ := stripN_cell hn hP hQ hPQ hH h1 h2 h3 hy
    exact B i hi1 hi2 (mem_of_In hf hPQ h)
  · obtain ⟨i, hi1, hi2, h⟩ := stripN_cell (x := m.1) (y := H - m.2) hn hQ hP hQP hH (by linarith)
      (by linarith) (by linarith) (by linarith)
    exact B i hi1 hi2 (mem_of_In hf hPQ (In_congr (In_my (b := H) h) rfl (by ring)))

/-- **Low strips, every `n`**: `n` unit squares in `[0,L] × [0,H]` with `H < 2√2 - 1` need `L ≥ n`. -/
theorem strip_low {n : ℕ} (hn1 : 1 ≤ n) {L H : ℝ} (hp : RectPacks n L H) (hH : H < 2 * Real.sqrt 2 - 1) :
    (n : ℝ) ≤ L := by
  rcases Nat.lt_or_ge n 2 with hn | hn
  · obtain rfl : n = 1 := by omega
    push_cast; exact (one_le_sides (hp.choose_spec.choose_spec.1 0)).1
  by_contra hL; push Not at hL
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  have hn' : (2 : ℝ) ≤ n := by exact_mod_cast hn
  -- widen to `L'` with `n H / L' < 2√2 - 1`, then scale by `n / L'`
  obtain ⟨M, hM1, hM2, hM5⟩ : ∃ M, L ≤ M ∧ n * H / (2 * Real.sqrt 2 - 1) ≤ M ∧ M < n :=
    ⟨_, le_max_left _ _, le_max_right _ _, max_lt hL (by rw [div_lt_iff₀ (by linarith)]; nlinarith)⟩
  obtain ⟨L', hL', hL'5, hL'H⟩ : ∃ L', L ≤ L' ∧ L' < n ∧ n * H / (2 * Real.sqrt 2 - 1) < L' :=
    ⟨(M + n) / 2, by linarith, by linarith, by linarith⟩
  have hH0 : 0 ≤ H := by
    have := (one_le_sides (hp.choose_spec.choose_spec.1 ⟨0, by omega⟩)).2; linarith
  have hL'0 : 0 < L' := by
    have : 0 ≤ n * H / (2 * Real.sqrt 2 - 1) := div_nonneg (by positivity) (by linarith)
    linarith
  have hp' := RectPacks.mono hp hL'
  set l := (n : ℝ) / L'
  have l1 : 1 < l := (one_lt_div hL'0).2 hL'5
  have e5 : l * L' = n := by simp only [l]; field_simp
  have eH : l * H ≤ 2 * Real.sqrt 2 - 1 := by
    have := (div_lt_iff₀ (by linarith : (0:ℝ) < 2 * Real.sqrt 2 - 1)).1 hL'H
    simp only [l]; rw [div_mul_eq_mul_div, div_le_iff₀ hL'0]; nlinarith
  refine no_pack_of_cover_n n (ptN n (l * H)) (fun _ => 1) (fun _ => zero_le_one) l1
    (fun c θ hc => stripN_cover hn eH (by rwa [e5] at hc)) ?_ hp'
  simp only [sum_const, card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
  rw [Nat.cast_sub (by omega)]; push_cast; linarith

end SquarePacking
