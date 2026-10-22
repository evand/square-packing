import Sqpack.StripPot

/-!
# (S): every square costs a grid unit or the zigzag amount

`kap_S`: `kap ≥ min (1 + β ||m| - |m'||) (Z - β(|m| + |m'|))`, `β = 1 - H/2`, `Z = 2 + 2√2 - 2H`, for
`2√2 - 1 ≤ H < 2`.  Tight for upright squares (first branch) and the 45° zigzag (second).

* Lines leaning the same way: `kap_same` (always the first branch).
* Opposite slopes `m ≥ y ≥ 0`, `m' = -y`: `kap = T m + T' (-y) + β(m - y)` with per-line parts
  `T x = g x/2 - |x|(1 - k)`, `T' x = g x/2 - |x|(H - 1 - k)`, `T' = T + 2β|x|`; the branches differ by
  `2βm - (Z - 1)`, so the threshold is `m ≷ μ* = (Z - 1)/(2β)`.  Per line, the minimum is at the vertical line or at the
  line along an edge; everything then reduces to
  - `G2`: `W + k ≥ 1` (`(1 - C)(1 - S) ≥ 0`),
  - `Z1` (`aligned`): both lines along edges, `V_m + V_y = (p + 1 - H)/(SC) ≥ Z`,
  - `X3`: the `m`-line at the threshold `μ* ≤ tan t`, the `y`-line along an edge.  In terms of `μ = μ*`
    (`H = (1 + 2√2 - 2μ)/(2 - μ)`), `(2 - μ) S X3 = F μ`, a convex quadratic in `μ` on `[0, min 1 (tan t)]`, decreasing at
    the right end; the end values reduce to `P1`, `P2`, `P3` (one-variable, `P1`, `P2` tight at 45°).
-/

namespace SquarePacking

/-- `P1` (`t ≥ 45°`): `(1 - S)² + C(S + 2 - 2√2) ≥ 0`. -/
lemma P1 {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hCS' : C ≤ S) :
    0 ≤ (1 - S) ^ 2 + C * (S + 2 - 2 * Real.sqrt 2) := by
  set r := Real.sqrt 2
  have r2 : r ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2g : 1.4142 < r := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have r2l : r < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hS1 : S ≤ 1 := by nlinarith
  have hSr : r / 2 ≤ S := by nlinarith [sq_nonneg (S - r / 2), sq_nonneg (S + r / 2)]
  rcases le_total (2 * r - 2) S with h | h
  · nlinarith [mul_nonneg hC (by linarith : (0:ℝ) ≤ S + 2 - 2 * r)]
  -- (1-S)³ ≥ (1+S)(2r-2-S)²  ⇐  -2(S - r/2) Q₁(S) ≥ 0, Q₁ ≤ 0 on [r/2, 2r-2]
  have Q1 : S ^ 2 + S * (1 - 3 * r / 2) - 11 * r / 2 + 8 ≤ 0 := by
    nlinarith [mul_nonneg (sub_nonneg.2 hSr) (sub_nonneg.2 h)]
  have cube : (1 + S) * (2 * r - 2 - S) ^ 2 ≤ (1 - S) ^ 3 := by
    have e : (1 - S) ^ 3 - (1 + S) * (2 * r - 2 - S) ^ 2 =
        -2 * (S - r / 2) * (S ^ 2 + S * (1 - 3 * r / 2) - 11 * r / 2 + 8) := by
      linear_combination (-(5 * S - 3) / 2) * r2
    nlinarith [mul_nonneg (sub_nonneg.2 hSr) (neg_nonneg.2 Q1)]
  -- so (1-S)⁴ ≥ C²(2r-2-S)², and both bases are ≥ 0
  have h4 : (C * (2 * r - 2 - S)) ^ 2 ≤ ((1 - S) ^ 2) ^ 2 := by
    have : C ^ 2 = (1 - S) * (1 + S) := by nlinarith
    nlinarith [mul_le_mul_of_nonneg_left cube (by linarith : (0:ℝ) ≤ 1 - S)]
  have := (sq_le_sq₀ (mul_nonneg hC (by linarith)) (by positivity)).1 h4
  nlinarith

/-- `P2` (`t ≤ 45°`): `2 + S - 2S² + C(1 - 2√2 + 2S - 2√2 S) ≥ 0`; squared, it is
`(16 - 8√2)(S - √2/2)²(S² + (1 + √2)S/2 - 1/4 + 3√2/8) ≥ 0`. -/
lemma P2 {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hSC : S ≤ C) :
    0 ≤ 2 + S - 2 * S ^ 2 + C * (1 - 2 * Real.sqrt 2 + 2 * S - 2 * Real.sqrt 2 * S) := by
  set r := Real.sqrt 2
  have r2 : r ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2g : 1.4142 < r := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have r2l : r < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hS2 : S ^ 2 ≤ 1 / 2 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  have hA : 0 ≤ 2 + S - 2 * S ^ 2 := by nlinarith
  have hB : 0 ≤ 2 * r - 1 + 2 * (r - 1) * S := by nlinarith
  have e : (2 + S - 2 * S ^ 2) ^ 2 - (1 - S ^ 2) * (2 * r - 1 + 2 * (r - 1) * S) ^ 2 =
      (16 - 8 * r) * (S - r / 2) ^ 2 * (S ^ 2 + S * (1 / 2 + r / 2) - 1 / 4 + 3 * r / 8) := by
    linear_combination ((16 * S ^ 4 + 16 * S ^ 3 - 8 * S ^ 2 * r + 12 * S ^ 2 + 4 * S * r ^ 2 - 16 * S * r +
      3 * r ^ 2 - 8 * r - 6) / 4) * r2
  have hq : 0 ≤ S ^ 2 + S * (1 / 2 + r / 2) - 1 / 4 + 3 * r / 8 := by nlinarith
  have h4 : (C * (2 * r - 1 + 2 * (r - 1) * S)) ^ 2 ≤ (2 + S - 2 * S ^ 2) ^ 2 := by
    have : C ^ 2 = 1 - S ^ 2 := by linarith
    have : 0 ≤ (16 - 8 * r) * (S - r / 2) ^ 2 * (S ^ 2 + S * (1 / 2 + r / 2) - 1 / 4 + 3 * r / 8) :=
      mul_nonneg (mul_nonneg (by linarith) (sq_nonneg _)) hq
    nlinarith
  have := (sq_le_sq₀ (mul_nonneg hC hB) hA).1 h4
  nlinarith

/-- `P3` (`t ≤ 45°`): `C(2S + 1) ≥ 1 + S²`. -/
lemma P3 {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hSC : S ≤ C) :
    1 + S ^ 2 ≤ C * (2 * S + 1) := by
  have hS2 : S ^ 2 ≤ 1 / 2 := by nlinarith
  have hS3 : S ≤ 0.7072 := by nlinarith
  have h4 : (1 + S ^ 2) ^ 2 ≤ (C * (2 * S + 1)) ^ 2 := by
    have e : (C * (2 * S + 1)) ^ 2 - (1 + S ^ 2) ^ 2 = S * (S + 4 - 4 * S ^ 2 - 5 * S ^ 3) := by
      linear_combination (2 * S + 1) ^ 2 * hCS
    have : 0 ≤ S + 4 - 4 * S ^ 2 - 5 * S ^ 3 := by nlinarith [mul_nonneg hS (sq_nonneg S)]
    nlinarith [mul_nonneg hS this]
  exact (sq_le_sq₀ (by positivity) (by positivity)).1 h4

/-- **`X3` in the variable `μ`**: `F μ ≥ 0` on `0 ≤ μ ≤ 1`, `μ C ≤ S`.  `F` is convex in `μ` and decreasing at the right
end `min 1 (tan t)`, where it is `P1` (`t ≥ 45°`) or `(1 - S) P2 / C²` (`t ≤ 45°`). -/
lemma X3F {C S u : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (_hu0 : 0 ≤ u) (hu1 : u ≤ 1)
    (huC : u * C ≤ S) :
    0 ≤ -C ^ 2 * u + 2 * C ^ 2 - C * S * u + 2 * C * S - 2 * C * Real.sqrt 2 + C * u + C - S ^ 2 * u ^ 2 +
      S ^ 2 * u + 2 * S ^ 2 + S * u ^ 2 - S * u - 2 * S := by
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  rcases le_total C S with h | h
  · have p1 := P1 hC hS hCS h
    have e : -C ^ 2 * u + 2 * C ^ 2 - C * S * u + 2 * C * S - 2 * C * Real.sqrt 2 + C * u + C - S ^ 2 * u ^ 2 +
        S ^ 2 * u + 2 * S ^ 2 + S * u ^ 2 - S * u - 2 * S =
        ((1 - S) ^ 2 + C * (S + 2 - 2 * Real.sqrt 2)) + (1 - C) * (1 - S) * (1 - u) + S * (1 - S) * (u - 1) ^ 2 := by
      linear_combination (2 - u) * hCS
    rw [e]
    have := mul_nonneg (mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)) (sub_nonneg.2 hu1)
    have := mul_nonneg (mul_nonneg hS (sub_nonneg.2 hS1)) (sq_nonneg (u - 1))
    linarith
  · have p2 := P2 hC hS hCS h
    have p3 := P3 hC hS hCS h
    have hCp : 0 < C := by nlinarith
    set F := -C ^ 2 * u + 2 * C ^ 2 - C * S * u + 2 * C * S - 2 * C * Real.sqrt 2 + C * u + C - S ^ 2 * u ^ 2 +
      S ^ 2 * u + 2 * S ^ 2 + S * u ^ 2 - S * u - 2 * S
    have e : C ^ 2 * F = (1 - S) * (2 + S - 2 * S ^ 2 + C * (1 - 2 * Real.sqrt 2 + 2 * S - 2 * Real.sqrt 2 * S)) -
        (1 - S) * (2 * C * S + C - S ^ 2 - 1) * (u * C - S) + S * (1 - S) * (u * C - S) ^ 2 := by
      simp only [F]
      linear_combination (-C ^ 2 * u + 2 * C ^ 2 - C * S * u + 2 * C * S - 2 * C * Real.sqrt 2 + C * u + C -
        2 * S + 2) * hCS
    have t1 := mul_nonneg (sub_nonneg.2 hS1) p2
    have t2 : 0 ≤ -(1 - S) * (2 * C * S + C - S ^ 2 - 1) * (u * C - S) := by
      have := mul_nonneg (mul_nonneg (sub_nonneg.2 hS1) (by linarith : (0:ℝ) ≤ 2 * C * S + C - S ^ 2 - 1))
        (sub_nonneg.2 huC)
      nlinarith
    have t3 := mul_nonneg (mul_nonneg hS (sub_nonneg.2 hS1)) (sq_nonneg (u * C - S))
    have : 0 ≤ C ^ 2 * F := by rw [e]; linarith
    exact (mul_nonneg_iff_of_pos_left (by positivity)).1 this

/-- **`Z1`, both lines along edges**: `Z · SC ≤ C + S + 1 - H` (that is `V_m + V_y ≥ Z`);
`(C + S + 1 - H) - Z(p² - 1)/2 = (√2 - p)((1 + √2 - H)p + 1 + √2 - √2 H) ≥ 0`. -/
lemma aligned {C S H : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (_hH1 : 2 * Real.sqrt 2 - 1 ≤ H)
    (hH2 : H ≤ 2) : (2 + 2 * Real.sqrt 2 - 2 * H) * (S * C) ≤ C + S + 1 - H := by
  set r := Real.sqrt 2
  have r2 : r ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have r2g : 1.4142 < r := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have hp2 : C + S ≤ r := (Real.le_sqrt (by linarith) (by norm_num)).2 (by nlinarith [sq_nonneg (C - S)])
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have h2 : 0 ≤ (1 + r - H) * (C + S) + 1 + r - r * H := by nlinarith
  have := mul_nonneg (sub_nonneg.2 hp2) h2
  have e : C + S + 1 - H - (2 + 2 * r - 2 * H) * (S * C) =
      (r - (C + S)) * ((1 + r - H) * (C + S) + 1 + r - r * H) := by
    linear_combination (1 + r - H) * hCS - (C + S + 1 - H) * r2
  linarith

end SquarePacking
