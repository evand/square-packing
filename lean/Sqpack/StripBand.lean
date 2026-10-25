import Sqpack.StripPotWalk

/-!
# Odd `n` in a strip, `H*(9) ≤ H ≤ 1.87`: the tangent potentials

Notation: `λ = 2√2 - H`, `a = 1 - λ`, `β = 1 - H/2`, `g = a + β`.  For a line of slope `m`, `x = |m|`.

* Opposite slopes have a closed form (`kap_opp`): for `x, y ≥ 0`,
  `kap C S H x (-y) = C + S + (S - H/2) x + (C - H/2) y + (C x - S)⁺ + (S y - C)⁺`.
* Potentials (`bandU0`, `bandU1`): with `A0 = 2a + βx`, `Lz = a + β(2 - x)` (the tangent at `x = 1` of the
  aligned-route cost `E`) and `T` (through `(8/7, Lz(8/7))` with `E`'s slope at `x = 4/3`), `L = max Lz T`:
  `u₀ = min A0 L`, `u₁ = min (a + βx) (a + L)`.
  They are the parity distances of the walk (prefixes: two uprights; one upright then a square aligned with two
  perpendicular lines) with `E` replaced by two tangents.

The proof that they are admissible: `StripBandLines` (one square against one line), `StripBandPairs` (the pair
table), `StripBandCore` (assembly), `StripBandMain` (the theorems).  This file: definitions, the pieces of the
potentials, and the three hand-proved polynomial inequalities where the table is tight or nearly so.
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

/-- `Lz`: the tangent at `x = 1` of the aligned-route cost `E`. -/
noncomputable def bandLz (x : ℝ) : ℝ := (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - x)

/-- `T`: through `(8/7, Lz(8/7))` with the slope `H/16 - 1/10` of `E` at `4/3` (close to `E`'s tangent there;
the fixed crossing `8/7` keeps every breakpoint free of square roots). -/
noncomputable def bandT (x : ℝ) : ℝ :=
  (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * (x - 8 / 7)

noncomputable def bandL (x : ℝ) : ℝ := max (bandLz H x) (bandT H x)

/-- Even lines: two uprights, or one upright and an aligned square. -/
noncomputable def bandU0 (m : ℝ) : ℝ :=
  min (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|) (bandL H |m|)

/-- Odd lines: one upright, or two uprights and an aligned square. -/
noncomputable def bandU1 (m : ℝ) : ℝ :=
  min ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * |m|) ((1 - (2 * Real.sqrt 2 - H)) + bandL H |m|)

end pots

/-! ## The band and its constants -/

/-- The band: `H*(9) = (14√2 - 5)/8 ≤ H ≤ 1.87` (`H*(9)`: where the zigzag for nine squares has length 9). -/
def InBand (H : ℝ) : Prop := (14 * Real.sqrt 2 - 5) / 8 ≤ H ∧ H ≤ 187 / 100

lemma sqrt2_lo' : (141421356 : ℝ) / 100000000 ≤ Real.sqrt 2 := by
  rw [Real.le_sqrt (by norm_num) (by norm_num)]; norm_num

lemma sqrt2_hi' : Real.sqrt 2 ≤ (141421357 : ℝ) / 100000000 := by
  rw [Real.sqrt_le_left (by norm_num)]; norm_num

lemma sqrt2_sq : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)

/-! ## The potentials piece by piece (`x ≥ 0`)
Breakpoints: `x₂ = (2b - a)/(2b)` (`A0 = Lz`), `8/7` (`Lz = T`), `y₁ = T₀/(b - T₁)` (`a + bx = a + T`), where
`a = 1 - λ`, `b = 1 - H/2`, `T = T₀ + T₁ x`; `y₁ ≥ 8/7` exactly when `H ≥ H*(9)`. -/

section pieces
variable {H x : ℝ}

lemma u0_A0 (hb : InBand H) (hx : 0 ≤ x)
    (h : x * (2 * (1 - H / 2)) ≤ 2 * (1 - H / 2) - (1 - (2 * Real.sqrt 2 - H))) :
    bandU0 H x = 2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * x := by
  unfold bandU0 bandL bandLz
  rw [abs_of_nonneg hx]
  apply min_eq_left
  exact le_max_of_le_left (by nlinarith)

lemma u0_Lz (hb : InBand H) (hx : 0 ≤ x)
    (h1 : 2 * (1 - H / 2) - (1 - (2 * Real.sqrt 2 - H)) ≤ x * (2 * (1 - H / 2))) (h2 : 7 * x ≤ 8) :
    bandU0 H x = bandLz H x := by
  obtain ⟨h9, h87⟩ := hb
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  unfold bandU0 bandL
  rw [abs_of_nonneg hx]
  have hTL : bandT H x ≤ bandLz H x := by
    unfold bandT bandLz; nlinarith
  rw [max_eq_left hTL]
  apply min_eq_right
  unfold bandLz; nlinarith

lemma u0_T (hb : InBand H) (hx : 8 ≤ 7 * x) : bandU0 H x = bandT H x := by
  obtain ⟨h9, h87⟩ := hb
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  unfold bandU0 bandL
  rw [abs_of_nonneg (by linarith)]
  have hLT : bandLz H x ≤ bandT H x := by
    unfold bandT bandLz; nlinarith
  rw [max_eq_right hLT]
  apply min_eq_right
  unfold bandT; nlinarith

lemma u1_O (hb : InBand H) (hx : 0 ≤ x)
    (h : x * ((1 - H / 2) - (H / 16 - 1 / 10)) ≤
      (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) - (H / 16 - 1 / 10) * (8 / 7)) :
    bandU1 H x = (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * x := by
  unfold bandU1 bandL
  rw [abs_of_nonneg hx]
  apply min_eq_left
  have : (1 - H / 2) * x ≤ bandT H x := by unfold bandT; nlinarith
  linarith [le_max_right (bandLz H x) (bandT H x)]

lemma u1_aT (hb : InBand H)
    (h : (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) - (H / 16 - 1 / 10) * (8 / 7) ≤
      x * ((1 - H / 2) - (H / 16 - 1 / 10))) :
    bandU1 H x = (1 - (2 * Real.sqrt 2 - H)) + bandT H x := by
  obtain ⟨h9, h87⟩ := hb
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hd : 0 < (1 - H / 2) - (H / 16 - 1 / 10) := by linarith
  -- y₁ ≥ 8/7 on the band, so x ≥ 8/7
  have hx : 8 ≤ 7 * x := by
    by_contra hc; push Not at hc
    have : x * ((1 - H / 2) - (H / 16 - 1 / 10)) < 8 / 7 * ((1 - H / 2) - (H / 16 - 1 / 10)) :=
      mul_lt_mul_of_pos_right (by linarith) hd
    nlinarith
  unfold bandU1 bandL
  rw [abs_of_nonneg (by linarith)]
  have hLT : bandLz H x ≤ bandT H x := by unfold bandT bandLz; nlinarith
  rw [max_eq_right hLT]
  apply min_eq_right
  unfold bandT; nlinarith

end pieces

/-! ## The hand-proved leaves: the 45° zigzag (`tight_P1`, `tight_P2`; `p = C + S`, `r = √2`) and 3-4-5 (`tight_P3`)
`P₁ = (p - 1)(r - 1)(r - p)`, `P₂ = (r - p)((p + 1)(r + 1) - H(r + p))`. -/

lemma tight_P1 {C S r : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hr2 : r ^ 2 = 2) (hr : 0 ≤ r) :
    0 ≤ C ^ 3 + C ^ 2 * S - C ^ 2 + C * S ^ 2 - 2 * C * S * r + 2 * C * S + S ^ 3 - S ^ 2 := by
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hpr : C + S ≤ r := by nlinarith [sq_nonneg (C - S), sq_nonneg (C + S - r)]
  have hr1 : 1 ≤ r := by nlinarith
  have e : C ^ 3 + C ^ 2 * S - C ^ 2 + C * S ^ 2 - 2 * C * S * r + 2 * C * S + S ^ 3 - S ^ 2 =
      (C + S - 1) * (r - 1) * (r - (C + S)) := by
    linear_combination (C + S + r - 2) * hCS + (1 - C - S) * hr2
  rw [e]
  exact mul_nonneg (mul_nonneg (by linarith) (by linarith)) (by linarith)

lemma tight_P2 {C S H r : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hr2 : r ^ 2 = 2) (hr : 0 ≤ r)
    (hH : H ≤ 2) :
    0 ≤ C ^ 3 - C ^ 2 * H + C ^ 2 * S + C ^ 2 + 2 * C * H * S + C * S ^ 2 - 2 * C * S * r - 2 * C * S - H * S ^ 2 +
      S ^ 3 + S ^ 2 := by
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hpr : C + S ≤ r := by nlinarith [sq_nonneg (C - S), sq_nonneg (C + S - r)]
  have hr1 : 1 ≤ r := by nlinarith
  have e : C ^ 3 - C ^ 2 * H + C ^ 2 * S + C ^ 2 + 2 * C * H * S + C * S ^ 2 - 2 * C * S * r - 2 * C * S - H * S ^ 2 +
      S ^ 3 + S ^ 2 = (r - (C + S)) * ((C + S + 1) * (r + 1) - H * (r + (C + S))) := by
    linear_combination (C + S + r + 2 - 2 * H) * hCS + (H - C - S - 1) * hr2
  rw [e]
  apply mul_nonneg (by linarith)
  nlinarith [mul_nonneg (sub_nonneg.2 hp1) (sub_nonneg.2 hr1)]

/-- **The 3-4-5 case.**  When the even potential's straight branch meets the odd potential's tangent `T`, the gap
vanishes to second order at the 3-4-5 angle (`tan θ = 3/4`, where `T` was tangent to the aligned-route cost) up to
the Pell gap `99/70 - √2`:
`(5 + 4C + 3S) · 560 P₃ = 7 (3C - 4S)² (40 + 48C + 40S - 5H(5 + 4C + 3S)) + 16 (5 + 4C + 3S) C S (99 - 70 r)`
on the circle. -/
lemma tight_P3 {C S H r : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hr : r ≤ 99 / 70)
    (hH : H ≤ 187 / 100) :
    0 ≤ (560 * C ^ 3 - 315 * C ^ 2 * H + 560 * C ^ 2 * S + 56 * C ^ 2 + 840 * C * H * S + 560 * C * S ^ 2 -
      1120 * C * S * r - 544 * C * S - 560 * H * S ^ 2 + 560 * S ^ 3 + 560 * S ^ 2) / 560 := by
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hw : 0 < 5 + 4 * C + 3 * S := by linarith
  have hf : 0 ≤ 40 + 48 * C + 40 * S - 5 * H * (5 + 4 * C + 3 * S) := by nlinarith
  have e : (5 + 4 * C + 3 * S) * ((560 * C ^ 3 - 315 * C ^ 2 * H + 560 * C ^ 2 * S + 56 * C ^ 2 + 840 * C * H * S +
      560 * C * S ^ 2 - 1120 * C * S * r - 544 * C * S - 560 * H * S ^ 2 + 560 * S ^ 3 + 560 * S ^ 2) / 560) =
      (7 * (3 * C - 4 * S) ^ 2 * (40 + 48 * C + 40 * S - 5 * H * (5 + 4 * C + 3 * S)) +
        16 * (5 + 4 * C + 3 * S) * (C * S) * (99 - 70 * r)) / 560 := by
    linear_combination (C + S) * (4 * C + 3 * S) * hCS
  have : 0 ≤ (5 + 4 * C + 3 * S) * ((560 * C ^ 3 - 315 * C ^ 2 * H + 560 * C ^ 2 * S + 56 * C ^ 2 + 840 * C * H * S +
      560 * C * S ^ 2 - 1120 * C * S * r - 544 * C * S - 560 * H * S ^ 2 + 560 * S ^ 3 + 560 * S ^ 2) / 560) := by
    rw [e]; apply div_nonneg _ (by norm_num)
    have : 0 ≤ 99 - 70 * r := by linarith
    positivity
  exact nonneg_of_mul_nonneg_right this hw

end SquarePacking
