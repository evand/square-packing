import Sqpack.StripPotWalk

/-!
# Odd `n` in a strip, `1.849 ≤ H < 1.87`: tangent potentials (work in progress)

Notation: `λ = 2√2 - H`, `a = 1 - λ`, `β = 1 - H/2`, `g = a + β`.  For a line of slope `m`, `x = |m|`.

* Opposite slopes have a closed form (`kap_opp`): for `x, y ≥ 0`,
  `kap C S H x (-y) = C + S + (S - H/2) x + (C - H/2) y + (C x - S)⁺ + (S y - C)⁺`.
* Potentials (`bandU0`, `bandU1`): with `A0 = 2a + βx`, `Lz = a + β(2 - x)` (tangent at `x = 1` of the
  aligned-route cost `E`), `T` = the tangent of `E` at `x = 4/3` (exact: `√(1 + 16/9) = 5/3`), `L = max Lz T`:
  `u₀ = min A0 L`, `u₁ = min (a + βx) (a + L)`.
  They are the parity distances of the walk (prefixes: two uprights; one upright then a square aligned with two
  perpendicular lines) with `E` replaced by two tangents; numerics in `search/rectfront/chainn/` (`pfix.py`,
  `lincheck.py`, `vert7.py`, `band_dec.py`).
-/

namespace SquarePacking

/-- **Opposite slopes in closed form.** -/
lemma kap_opp {C S H x y : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hx : 0 ≤ x) (hy : 0 ≤ y) :
    kap C S H x (-y) = C + S + (S - H / 2) * x + (C - H / 2) * y + max (C * x - S) 0 + max (S * y - C) 0 := by
  unfold kap
  have e1 : |C + x * S| = C + x * S := abs_of_nonneg (by positivity)
  have e2 : |C + -y * S| = |S * y - C| := by rw [← abs_neg]; ring_nf
  have e3 : |-y * C - S| = y * C + S := by
    rw [show -y * C - S = -(y * C + S) by ring, abs_neg]; exact abs_of_nonneg (by positivity)
  have e4 : |x - -y| = x + y := by rw [sub_neg_eq_add]; exact abs_of_nonneg (by linarith)
  have m1 : |x * C - S| = 2 * max (C * x - S) 0 - (C * x - S) := by
    rcases le_total 0 (C * x - S) with h | h
    · rw [max_eq_left h, abs_of_nonneg (by linarith)]; ring
    · rw [max_eq_right h, abs_of_nonpos (by linarith)]; ring
  have m2 : |S * y - C| = 2 * max (S * y - C) 0 - (S * y - C) := by
    rcases le_total 0 (S * y - C) with h | h
    · rw [max_eq_left h, abs_of_nonneg h]; ring
    · rw [max_eq_right h, abs_of_nonpos h]; ring
  rw [e1, e2, e3, e4, m1, m2]
  ring

section pots
variable (H : ℝ)

/-- `T`: the tangent at `x = 4/3` of the aligned-route cost `E` (`E(4/3) = 7H/12 + 14/3 - 4√2`,
`E'(4/3) = H/16 - 1/10`). -/
noncomputable def bandT (x : ℝ) : ℝ :=
  7 * H / 12 + 14 / 3 - 4 * Real.sqrt 2 + (H / 16 - 1 / 10) * (x - 4 / 3)

/-- `Lz`: the tangent of `E` at `x = 1`. -/
noncomputable def bandLz (x : ℝ) : ℝ := (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - x)

noncomputable def bandL (x : ℝ) : ℝ := max (bandLz H x) (bandT H x)

/-- Even lines: two uprights, or one upright and an aligned square. -/
noncomputable def bandU0 (m : ℝ) : ℝ :=
  min (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|) (bandL H |m|)

/-- Odd lines: one upright, or two uprights and an aligned square. -/
noncomputable def bandU1 (m : ℝ) : ℝ :=
  min ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|) ((1 - (2 * Real.sqrt 2 - H)) + bandL H |m|)

end pots

end SquarePacking
