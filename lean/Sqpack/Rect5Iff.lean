import Sqpack.Rect5
import Sqpack.Rect5Low

/-!
# `R₅`: five unit squares in a rectangle, exactly

`rect5_iff_of`: given the real inequality `Final5Hi` (strips of height in `[2√2 - 1, 2)`), five unit squares fit in
`[0,w] × [0,h]` iff both sides are at least `1` and one of
* a side is at least `5` (the row);
* both sides are at least `2` and a side is at least `3` (the `2 × 3` grid) or `w + h ≥ 4 + √2` (the stretched
  `s(5)` packing; at `w = h` this is `s(5) = 2 + 1/√2`, Göbel 1979);
* a side `s < 2` and the other `ℓ` with `ℓ + 4 s ≥ 4 + 6√2` (the zigzag).
Boundary: Friedman, Math Magic June 2010 (pieces 4.1–4.5 of `S₄`).  Lower bounds: low strips by four midline
points (`rect5_low`), high strips by the chain (`rect5_strip_lower_hi_of`), both sides `≥ 2` by weighted point
covers (`rect5_wide_iff`).
-/

namespace SquarePacking

/-- The strip part: in `[0,L] × [0,H]` with `H < 2`, five squares need `L ≥ 5` or `L + 4H ≥ 4 + 6√2`. -/
theorem rect5_strip_of (hF : Final5Hi) {L H : ℝ} (hp : RectPacks 5 L H) (hH : H < 2) :
    5 ≤ L ∨ 4 + 6 * Real.sqrt 2 ≤ L + 4 * H := by
  rcases lt_or_ge H (2 * Real.sqrt 2 - 1) with h | h
  · exact Or.inl (rect5_low hp h)
  · exact rect5_strip_lower_hi_of hF hp hH h

/-- The condition for five unit squares to fit in `[0,w] × [0,h]`. -/
def Fits5 (w h : ℝ) : Prop :=
  1 ≤ w ∧ 1 ≤ h ∧ (5 ≤ w ∨ 5 ≤ h ∨ (2 ≤ w ∧ 2 ≤ h ∧ (3 ≤ w ∨ 3 ≤ h ∨ 4 + Real.sqrt 2 ≤ w + h)) ∨
    (h < 2 ∧ 4 + 6 * Real.sqrt 2 ≤ w + 4 * h) ∨ (w < 2 ∧ 4 + 6 * Real.sqrt 2 ≤ h + 4 * w))

/-- **`R₅` exactly**, given `Final5Hi`. -/
theorem rect5_iff_of (hF : Final5Hi) {w h : ℝ} : RectPacks 5 w h ↔ Fits5 w h := by
  constructor
  · intro hp
    obtain ⟨hw, hh⟩ := one_le_sides (hp.choose_spec.choose_spec.1 0)
    refine ⟨hw, hh, ?_⟩
    rcases lt_or_ge h 2 with h2 | h2
    · rcases rect5_strip_of hF hp h2 with h5 | hz
      · exact Or.inl h5
      · exact Or.inr (Or.inr (Or.inr (Or.inl ⟨h2, hz⟩)))
    rcases lt_or_ge w 2 with w2 | w2
    · rcases rect5_strip_of hF hp.swap w2 with h5 | hz
      · exact Or.inr (Or.inl h5)
      · exact Or.inr (Or.inr (Or.inr (Or.inr ⟨w2, hz⟩)))
    exact Or.inr (Or.inr (Or.inl ⟨w2, h2, (rect5_wide_iff w2 h2).1 hp⟩))
  · rintro ⟨hw, hh, h5 | h5 | ⟨w2, h2, hc⟩ | ⟨h2, hz⟩ | ⟨w2, hz⟩⟩
    · exact rect5_row h5 hh
    · exact (rect5_row h5 hw).swap
    · exact (rect5_wide_iff w2 h2).2 hc
    · obtain ⟨r2g, r2l⟩ := sqrt2_bounds
      rcases lt_or_ge h (Real.sqrt 2) with hs | hs
      · exact rect5_row (by linarith) hh
      · exact rect5_zigzag hs h2.le hz
    · obtain ⟨r2g, r2l⟩ := sqrt2_bounds
      rcases lt_or_ge w (Real.sqrt 2) with hs | hs
      · exact (rect5_row (by linarith) hw).swap
      · exact (rect5_zigzag hs w2.le hz).swap

end SquarePacking
