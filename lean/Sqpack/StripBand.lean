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

/-! ## A convex piece is nonnegative if it is at its ends and its kink

`f(x) = A + B x + max (C x - S) 0` is convex with one kink at `S / C`; on `[ln/ld, hn/hd]` its least value is at an
end or at the kink.  With a positive multiplier `m` and an added constant `E` (for the nested use: first in `x`,
then in `y`).  Ends are fractions so that no division appears. -/

/-- An affine function nonnegative at two fractions `n₁/d₁ ≤ x ≤ n₂/d₂` is nonnegative at `x`. -/
lemma aff_frac {α β x n₁ d₁ n₂ d₂ : ℝ} (hd₁ : 0 < d₁) (hd₂ : 0 < d₂) (h₁ : n₁ ≤ x * d₁) (h₂ : x * d₂ ≤ n₂)
    (e₁ : 0 ≤ α * d₁ + β * n₁) (e₂ : 0 ≤ α * d₂ + β * n₂) : 0 ≤ α + β * x := by
  rcases le_total 0 β with hβ | hβ
  · have : 0 ≤ (α + β * x) * d₁ := by nlinarith [mul_le_mul_of_nonneg_left h₁ hβ]
    exact nonneg_of_mul_nonneg_left this hd₁
  · have : 0 ≤ (α + β * x) * d₂ := by nlinarith [mul_le_mul_of_nonpos_left h₂ hβ]
    exact nonneg_of_mul_nonneg_left this hd₂

/-- Nondecreasing affine function, nonnegative at `n₁/d₁ ≤ x`. -/
lemma aff_lo {α β x n₁ d₁ : ℝ} (hd₁ : 0 < d₁) (h₁ : n₁ ≤ x * d₁) (hβ : 0 ≤ β) (e₁ : 0 ≤ α * d₁ + β * n₁) :
    0 ≤ α + β * x := by
  have : 0 ≤ (α + β * x) * d₁ := by nlinarith [mul_le_mul_of_nonneg_left h₁ hβ]
  exact nonneg_of_mul_nonneg_left this hd₁

/-- **A convex piece** `m (A + B x + max (C x - S) 0) + E` on `[ln/ld, hn/hd]` is nonnegative if it is at both ends
and, when the kink `S/C` lies in the interval, at the kink. -/
lemma plx {m E A B C S x ln ld hn hd F : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hm : 0 < m) (hld : 0 < ld)
    (hhd : 0 < hd) (hl : ln ≤ x * ld) (hh : x * hd ≤ hn)
    (h1 : 0 ≤ m * (A * ld + B * ln + max (C * ln - S * ld) 0) + E * ld)
    (h2 : 0 ≤ m * (A * hd + B * hn + max (C * hn - S * hd) 0) + E * hd)
    (hk : 0 < C → ln * C ≤ S * ld → S * hd ≤ hn * C → 0 ≤ m * (A * C + B * S) + E * C)
    (hF : F = m * (A + B * x + max (C * x - S) 0) + E) : 0 ≤ F := by
  rw [hF]
  rcases le_total (C * x) S with hx | hx
  · rw [max_eq_right (by linarith)]
    have hl0 : C * ln ≤ S * ld := by nlinarith [mul_le_mul_of_nonneg_left hl hC]
    rw [max_eq_right (by linarith)] at h1
    have e₁ : 0 ≤ (m * A + E) * ld + m * B * ln := by linarith
    rcases le_or_gt (C * hn) (S * hd) with hh0 | hh0
    · rw [max_eq_right (by linarith)] at h2
      have := aff_frac (α := m * A + E) (β := m * B) hld hhd hl hh e₁ (by linarith)
      linarith
    · have hCp : 0 < C := by
        rcases eq_or_lt_of_le hC with h0 | h0
        · rw [← h0] at hh0; nlinarith
        · exact h0
      have ek := hk hCp (by linarith) (by linarith)
      have := aff_frac (α := m * A + E) (β := m * B) (n₂ := S) (d₂ := C) hld hCp hl (by linarith) e₁ (by linarith)
      linarith
  · rw [max_eq_left (by linarith)]
    have hh0 : S * hd ≤ C * hn := by nlinarith [mul_le_mul_of_nonneg_left hh hC]
    rw [max_eq_left (by linarith)] at h2
    have e₂ : 0 ≤ (m * A - m * S + E) * hd + m * (B + C) * hn := by linarith
    rcases le_or_gt (S * ld) (C * ln) with hl0 | hl0
    · rw [max_eq_left (by linarith)] at h1
      have := aff_frac (α := m * A - m * S + E) (β := m * (B + C)) hld hhd hl hh (by linarith) e₂
      linarith
    · have hCp : 0 < C := by
        rcases eq_or_lt_of_le hC with h0 | h0
        · rw [← h0] at hx hl0; nlinarith
        · exact h0
      have ek := hk hCp (by linarith) (by nlinarith)
      have := aff_frac (α := m * A - m * S + E) (β := m * (B + C)) (n₁ := S) (d₁ := C) hCp hhd (by linarith) hh (by linarith) e₂
      linarith

/-- **A convex piece on a half-line** `[ln/ld, ∞)`: nonnegative at the start, at the kink (if after the start), and
nondecreasing at infinity (`B + C ≥ 0`). -/
lemma plx_inf {m E A B C S x ln ld F : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hm : 0 < m) (hld : 0 < ld)
    (hl : ln ≤ x * ld)
    (h1 : 0 ≤ m * (A * ld + B * ln + max (C * ln - S * ld) 0) + E * ld)
    (hs : 0 ≤ B + C)
    (hk : 0 < C → ln * C ≤ S * ld → 0 ≤ m * (A * C + B * S) + E * C)
    (hF : F = m * (A + B * x + max (C * x - S) 0) + E) : 0 ≤ F := by
  rw [hF]
  rcases le_total (C * x) S with hx | hx
  · rw [max_eq_right (by linarith)]
    have hl0 : C * ln ≤ S * ld := by nlinarith [mul_le_mul_of_nonneg_left hl hC]
    rw [max_eq_right (by linarith)] at h1
    have e₁ : 0 ≤ (m * A + E) * ld + m * B * ln := by linarith
    rcases le_or_gt 0 B with hB | hB
    · have := aff_lo (α := m * A + E) (β := m * B) hld hl (by positivity) e₁
      linarith
    · have hCp : 0 < C := by linarith
      have ek := hk hCp (by nlinarith)
      have := aff_frac (α := m * A + E) (β := m * B) (n₂ := S) (d₂ := C) hld hCp hl (by linarith) e₁ (by linarith)
      linarith
  · rw [max_eq_left (by linarith)]
    have hmBC : 0 ≤ m * (B + C) := by positivity
    rcases le_or_gt (S * ld) (C * ln) with hl0 | hl0
    · rw [max_eq_left (by linarith)] at h1
      have := aff_lo (α := m * A - m * S + E) (β := m * (B + C)) hld hl hmBC (by linarith)
      linarith
    · have hCp : 0 < C := by
        rcases eq_or_lt_of_le hC with h0 | h0
        · rw [← h0] at hx hl0; nlinarith
        · exact h0
      have ek := hk hCp (by linarith)
      have := aff_lo (x := x) (α := m * A - m * S + E) (β := m * (B + C)) (n₁ := S) (d₁ := C) hCp (by linarith) hmBC (by linarith)
      linarith

/-! ## The two tight leaves (equality at the 45° zigzag), `p = C + S`, `r = √2`
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

end SquarePacking
