import Sqpack.Rect5Line

/-!
# Five unit squares in a low strip: four points on the midline

`rect5_low`: in a strip of height `H < 2√2 - 1 ≈ 1.828`, five unit squares need length `L ≥ 5`.
Cover: in `[0,5] × [0,H]` with `H ≤ 2√2 - 1`, every closed unit square contains one of `(i, H/2)`, `i = 1, …, 4`
(`strip_cover`): a centre within `1` of an end wall gets the end point (`In_corner`), a centre between two
points gets one of them (`In_side2`: gap `1`, distance `H/2` from the floor, `1/2 + H/2 ≤ √2`).  Widen the
strip, scale it onto `L = 5`, and four points cannot hold five disjoint squares.
-/

open Finset
open scoped Classical

namespace SquarePacking

/-- The lower half of the strip: some midline point `(i, H/2)` is in the square. -/
lemma strip_cell {H c s x y : ℝ} (hc : 0 ≤ c) (hs : 0 ≤ s) (hcs : c ^ 2 + s ^ 2 = 1)
    (hH : H ≤ 2 * Real.sqrt 2 - 1) (hx1 : (c + s) / 2 ≤ x) (hx2 : x ≤ 5 - (c + s) / 2)
    (hy1 : (c + s) / 2 ≤ y) (hy2 : y ≤ H / 2) :
    In (1 - x) (H / 2 - y) c s ∨ In (2 - x) (H / 2 - y) c s ∨ In (3 - x) (H / 2 - y) c s ∨
      In (4 - x) (H / 2 - y) c s := by
  obtain ⟨p1, p2, -, -⟩ := frame_facts hc hs hcs
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  have ab : ∀ {z : ℝ}, -(1 - (c + s) / 2) ≤ z → z ≤ 1 - (c + s) / 2 → |z| ≤ 1 - (c + s) / 2 :=
    fun h h' => abs_le.2 ⟨h, h'⟩
  have side : ∀ i : ℝ, i ≤ x → x ≤ i + 1 → In (i + 1 - x) (H / 2 - y) c s ∨ In (i - x) (H / 2 - y) c s :=
    fun i h1 h2 => by
      rcases In_side2 hc hs hcs (A := i + 1 - x) (A' := x - i) (B := H / 2 - y) (y₀ := H / 2)
        (by linarith) (by linarith) (by linarith) (by linarith) (by linarith) (by linarith) (by linarith)
        with h | h
      · exact Or.inl h
      · exact Or.inr (In_congr h (by ring) rfl)
  rcases le_total x 1 with x1 | x1
  · exact Or.inl (In_corner hc hs hcs (ab (by linarith) (by linarith)) (ab (by linarith) (by linarith)))
  rcases le_total x 2 with x2 | x2
  · rcases side 1 x1 (by linarith) with h | h
    · exact Or.inr (Or.inl (In_congr h (by ring) rfl))
    · exact Or.inl h
  rcases le_total x 3 with x3 | x3
  · rcases side 2 x2 (by linarith) with h | h
    · exact Or.inr (Or.inr (Or.inl (In_congr h (by ring) rfl)))
    · exact Or.inr (Or.inl h)
  rcases le_total x 4 with x4 | x4
  · rcases side 3 x3 (by linarith) with h | h
    · exact Or.inr (Or.inr (Or.inr (In_congr h (by ring) rfl)))
    · exact Or.inr (Or.inr (Or.inl h))
  · exact Or.inr (Or.inr (Or.inr
      (In_corner hc hs hcs (ab (by linarith) (by linarith)) (ab (by linarith) (by linarith)))))

/-- The four midline points. -/
noncomputable def ptS (H : ℝ) : Fin 4 → ℝ × ℝ := ![(1, H / 2), (2, H / 2), (3, H / 2), (4, H / 2)]

/-- **The midline cover** of `[0,5] × [0,H]`, `H ≤ 2√2 - 1`. -/
theorem strip_cover {H : ℝ} (hH : H ≤ 2 * Real.sqrt 2 - 1) {m : ℝ × ℝ} {θ : ℝ} (hin : sq m θ 1 ⊆ rbox 5 H) :
    1 ≤ ∑ j ∈ univ.filter (fun j => ptS H j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j := by
  obtain ⟨P, Q, hP, hQ, hPQ, h1, h2, h3, h4, hf⟩ := framed hin
  have hQP : Q ^ 2 + P ^ 2 = 1 := by linarith
  have B : ∀ j : Fin 4, In ((ptS H j).1 - m.1) ((ptS H j).2 - m.2) P Q →
      1 ≤ ∑ j ∈ univ.filter (fun j => ptS H j ∈ sq m θ 1), (fun _ => (1 : ℝ)) j :=
    fun j hj => cap_of_block (fun _ => zero_le_one) {j} (by simpa using mem_of_In hf hPQ hj) (by simp)
  rcases le_total m.2 (H / 2) with hy | hy
  · rcases strip_cell hP hQ hPQ hH h1 h2 h3 hy with h | h | h | h
    · exact B 0 (by simpa [ptS] using h)
    · exact B 1 (by simpa [ptS] using h)
    · exact B 2 (by simpa [ptS] using h)
    · exact B 3 (by simpa [ptS] using h)
  · rcases strip_cell (x := m.1) (y := H - m.2) hQ hP hQP hH (by linarith) (by linarith) (by linarith)
      (by linarith) with h | h | h | h
    · exact B 0 (by simp [ptS]; exact In_congr (In_my (b := H) h) rfl (by ring))
    · exact B 1 (by simp [ptS]; exact In_congr (In_my (b := H) h) rfl (by ring))
    · exact B 2 (by simp [ptS]; exact In_congr (In_my (b := H) h) rfl (by ring))
    · exact B 3 (by simp [ptS]; exact In_congr (In_my (b := H) h) rfl (by ring))

/-- **Low strips**: five unit squares in `[0,L] × [0,H]` with `H < 2√2 - 1` need `L ≥ 5`. -/
theorem rect5_low {L H : ℝ} (hp : RectPacks 5 L H) (hH : H < 2 * Real.sqrt 2 - 1) : 5 ≤ L := by
  by_contra hL; push Not at hL
  obtain ⟨r2g, r2l⟩ := sqrt2_bounds
  have hH1 : 1 ≤ H := (one_le_sides (hp.choose_spec.choose_spec.1 0)).2
  have hL1 : 1 ≤ L := (one_le_sides (hp.choose_spec.choose_spec.1 0)).1
  -- widen to `L'` with `5 H / L' < 2√2 - 1`, then scale by `5 / L'`
  obtain ⟨M, hM1, hM2, hM5⟩ : ∃ M, L ≤ M ∧ 5 * H / (2 * Real.sqrt 2 - 1) ≤ M ∧ M < 5 :=
    ⟨_, le_max_left _ _, le_max_right _ _, max_lt hL (by rw [div_lt_iff₀ (by linarith)]; nlinarith)⟩
  obtain ⟨L', hL', hL'5, hL'H⟩ : ∃ L', L ≤ L' ∧ L' < 5 ∧ 5 * H / (2 * Real.sqrt 2 - 1) < L' :=
    ⟨(M + 5) / 2, by linarith, by linarith, by linarith⟩
  have hL'0 : 0 < L' := by linarith
  have hp' := RectPacks.mono hp hL'
  set l := 5 / L'
  have l1 : 1 < l := (one_lt_div hL'0).2 hL'5
  have e5 : l * L' = 5 := by simp only [l]; field_simp
  have eH : l * H ≤ 2 * Real.sqrt 2 - 1 := by
    have := (div_lt_iff₀ (by linarith : (0:ℝ) < 2 * Real.sqrt 2 - 1)).1 hL'H
    simp only [l]; rw [div_mul_eq_mul_div, div_le_iff₀ hL'0]; nlinarith
  exact no_pack_of_cover (ptS (l * H)) (fun _ => 1) (fun _ => zero_le_one) l1
    (fun c θ hc => strip_cover eH (by rwa [e5] at hc)) (by simp; norm_num) hp'

end SquarePacking
