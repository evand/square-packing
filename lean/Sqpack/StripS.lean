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

/-- Per line (`m ≥ 0`): `T m ≥ W`, the value on the edge line (`C W = C k + S(S - 1)`). -/
lemma line_W {C S m : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (_hm : 0 ≤ m) :
    C * ((C + S) / 2) + S * (S - 1) ≤ C * ((C + m * S + |m * C - S|) / 2 - m * (1 - (C + S) / 2)) := by
  have hS1 : S ≤ 1 := by nlinarith
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  rcases le_total S (m * C) with h | h
  · rw [abs_of_nonneg (by linarith)]; nlinarith [mul_nonneg (sub_nonneg.2 h) (by linarith : (0:ℝ) ≤ C + S - 1)]
  · rw [abs_of_nonpos (by linarith)]; nlinarith [mul_nonneg (sub_nonneg.2 h) (sub_nonneg.2 hS1)]

/-- Per line (`y ≥ 0`, slope `-y`): `T'(-y) ≥ min k V_y` (`S V_y = S k + C(C + 1 - H)`). -/
lemma line_Ly {C S H y : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hH : H ≤ 2) (hy : 0 ≤ y) :
    (C + S) / 2 ≤ (|C - y * S| + (y * C + S)) / 2 - y * (H - 1 - (C + S) / 2) ∨
      S * ((C + S) / 2) + C * (C + 1 - H) ≤ S * ((|C - y * S| + (y * C + S)) / 2 - y * (H - 1 - (C + S) / 2)) := by
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  rcases le_total (y * S) C with h | h
  · rw [abs_of_nonneg (by linarith)]
    rcases le_total 0 (C + 1 - H) with h2 | h2
    · left; nlinarith [mul_nonneg hy h2]
    · right; nlinarith [mul_nonneg (sub_nonneg.2 h) (neg_nonneg.2 h2)]
  · rw [abs_of_nonpos (by linarith)]
    right; nlinarith [mul_nonneg (sub_nonneg.2 h) (by linarith : (0:ℝ) ≤ C + S + 1 - H)]

/-- **`X3`** (the `m`-line at the threshold `μ*`, the `y`-line along an edge), in the strip variables:
`(2 - H) μ = 1 + 2√2 - 2H`, `0 ≤ μ`, `μ C ≤ S` ⇒ `S (p - 1 - μ(1 - S)) + C (C + 1 - H) ≥ 0`.
Multiplied by `2 - μ` it is `X3F μ`. -/
lemma X3 {C S H μ : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hH1 : 2 * Real.sqrt 2 - 1 ≤ H)
    (hH2 : H < 2) (hμ : (2 - H) * μ = 1 + 2 * Real.sqrt 2 - 2 * H) (hμ0 : 0 ≤ μ) (hμC : μ * C ≤ S) :
    0 ≤ S * (C + S - 1 - μ * (1 - S)) + C * (C + 1 - H) := by
  have hμ1 : μ ≤ 1 := by
    by_contra h; push Not at h
    have : (2 - H) * 1 < (2 - H) * μ := mul_lt_mul_of_pos_left h (by linarith)
    linarith
  have F := X3F hC hS hCS hμ0 hμ1 hμC
  -- (2 - μ) · goal = F μ, using H (2 - μ) = 1 + 2√2 - 2μ
  have hH : H * (2 - μ) = 1 + 2 * Real.sqrt 2 - 2 * μ := by linear_combination hμ
  have e : (2 - μ) * (S * (C + S - 1 - μ * (1 - S)) + C * (C + 1 - H)) =
      -C ^ 2 * μ + 2 * C ^ 2 - C * S * μ + 2 * C * S - 2 * C * Real.sqrt 2 + C * μ + C - S ^ 2 * μ ^ 2 +
        S ^ 2 * μ + 2 * S ^ 2 + S * μ ^ 2 - S * μ - 2 * S := by
    linear_combination (-C) * hH
  have : 0 ≤ (2 - μ) * (S * (C + S - 1 - μ * (1 - S)) + C * (C + 1 - H)) := by rw [e]; exact F
  exact (mul_nonneg_iff_of_pos_left (by linarith)).1 this

/-- **(S), opposite slopes** `m ≥ y ≥ 0`, `m' = -y`. -/
lemma kap_S_opp {C S H m y : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H < 2) (hfit : C + S ≤ H) (hy : 0 ≤ y) (hym : y ≤ m) :
    min (1 + (1 - H / 2) * (m - y)) (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (m + y)) ≤ kap C S H m (-y) := by
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hC1 : C ≤ 1 := by nlinarith
  have hS1 : S ≤ 1 := by nlinarith
  have hm : 0 ≤ m := hy.trans hym
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  -- per-line parts as plain variables
  have hW0 := line_W hC hS hCS hm
  have hLy0 := line_Ly hC hS hCS hH2.le hy
  obtain ⟨k, hk⟩ : ∃ k : ℝ, k = (C + S) / 2 := ⟨_, rfl⟩
  obtain ⟨Tm, hTm⟩ : ∃ Tm : ℝ, Tm = (C + m * S + |m * C - S|) / 2 - m * (1 - (C + S) / 2) := ⟨_, rfl⟩
  obtain ⟨Ty, hTy⟩ : ∃ Ty : ℝ, Ty = (|C - y * S| + (y * C + S)) / 2 - y * (H - 1 - (C + S) / 2) := ⟨_, rfl⟩
  rw [← hTm, ← hk] at hW0
  rw [← hTy, ← hk] at hLy0
  have hkap : kap C S H m (-y) = Tm + Ty + (1 - H / 2) * (m - y) := by
    unfold kap
    rw [show m - -y = m + y by ring, abs_of_nonneg (by linarith : 0 ≤ m + y),
      abs_of_nonneg (by positivity : 0 ≤ C + m * S), show -y * C - S = -(y * C + S) by ring, abs_neg,
      abs_of_nonneg (by positivity : 0 ≤ y * C + S), show C + -y * S = C - y * S by ring, hTm, hTy]
    ring
  rw [hkap]
  -- upright squares: the grid branch (T m = T'(-y) = 1/2)
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · subst hS0
    have hC1' : C = 1 := by nlinarith
    subst hC1'
    apply min_le_of_left_le
    rw [hTm, hTy, abs_of_nonneg (by linarith : (0:ℝ) ≤ m * 1 - 0), abs_of_nonneg (by linarith : (0:ℝ) ≤ 1 - y * 0)]
    nlinarith
  rcases eq_or_lt_of_le hC with hC0 | hCp
  · subst hC0
    have hS1' : S = 1 := by nlinarith
    subst hS1'
    apply min_le_of_left_le
    rw [hTm, hTy, abs_of_nonpos (by linarith : m * 0 - 1 ≤ (0:ℝ)), abs_of_nonpos (by linarith : (0:ℝ) - y * 1 ≤ 0)]
    nlinarith
  have Z1 := aligned hC hS hCS hH1 hH2.le
  -- S (C V_m) + C (S V_y) = p + 1 - H
  have hid : S * (C * k + S * (S + 1 - H)) + C * (S * k + C * (C + 1 - H)) = C + S + 1 - H := by
    rw [hk]; linear_combination (C + S + 1 - H) * hCS
  have hTpm : Tm + Ty + (1 - H / 2) * (m - y) = (Tm + (2 - H) * m) + Ty - (1 - H / 2) * (m + y) := by ring
  rcases hLy0 with hy1 | hy1
  · -- the y-line is at least k: grid branch from W + k ≥ 1, i.e. (1 - C)(1 - S) ≥ 0
    apply min_le_of_left_le
    have hg : C * k + S * (S - 1) + C * k - C = (1 - C) * (1 - S) := by rw [hk]; linear_combination hCS
    have h1 := mul_nonneg (sub_nonneg.2 hC1) (sub_nonneg.2 hS1)
    have h2 := mul_le_mul_of_nonneg_left hy1 hCp.le
    have : C * 1 ≤ C * (Tm + Ty) := by linarith
    have := (mul_le_mul_iff_of_pos_left hCp).1 this
    linarith
  -- the y-line is at least V_y
  have hyC := mul_le_mul_of_nonneg_left hy1 hCp.le
  rcases le_total S (m * C) with hmk | hmk
  · -- m past its kink: T'(m) ≥ V_m, then Z1: the zigzag branch
    apply min_le_of_right_le
    have hVm : C * k + S * (S + 1 - H) ≤ C * (Tm + (2 - H) * m) := by
      rw [hTm, hk, abs_of_nonneg (by linarith : 0 ≤ m * C - S)]
      nlinarith [mul_nonneg (sub_nonneg.2 hmk) (by linarith : (0:ℝ) ≤ C + S + 1 - H)]
    have hVmS := mul_le_mul_of_nonneg_left hVm hSp.le
    have : S * C * (2 + 2 * Real.sqrt 2 - 2 * H) ≤ S * C * ((Tm + (2 - H) * m) + Ty) := by
      have e : S * C * ((Tm + (2 - H) * m) + Ty) = S * (C * (Tm + (2 - H) * m)) + C * (S * Ty) := by ring
      rw [e]; linarith
    have := (mul_le_mul_iff_of_pos_left (mul_pos hSp hCp)).1 this
    rw [hTpm]; linarith
  -- m before its kink: T m = k + m(S - 1)
  have hTm' : Tm = k + m * (S - 1) := by
    rw [hTm, hk, abs_of_nonpos (by linarith : m * C - S ≤ 0)]; ring
  have hSk : S * k = S * (C + S) / 2 := by rw [hk]; ring
  obtain ⟨μ, hμdef⟩ : ∃ μ : ℝ, μ = (1 + 2 * Real.sqrt 2 - 2 * H) / (2 - H) := ⟨_, rfl⟩
  have hμ : (2 - H) * μ = 1 + 2 * Real.sqrt 2 - 2 * H := by
    rw [hμdef]; field_simp [(by linarith : (0:ℝ) < 2 - H).ne']
  have hμS : S * ((2 - H) * μ) = S * (1 + 2 * Real.sqrt 2 - 2 * H) := by rw [hμ]
  rcases le_total μ m with hmμ | hmμ
  · -- zigzag branch
    apply min_le_of_right_le
    rw [hTpm]
    rcases le_total (S + 1 - H) 0 with hc | hc
    · -- T'(m) decreasing in m: ≥ V_m
      have hVm : C * k + S * (S + 1 - H) ≤ C * (Tm + (2 - H) * m) := by
        rw [hTm']; nlinarith [mul_nonneg (sub_nonneg.2 hmk) (neg_nonneg.2 hc)]
      have hVmS := mul_le_mul_of_nonneg_left hVm hSp.le
      have : S * C * (2 + 2 * Real.sqrt 2 - 2 * H) ≤ S * C * ((Tm + (2 - H) * m) + Ty) := by
        have e : S * C * ((Tm + (2 - H) * m) + Ty) = S * (C * (Tm + (2 - H) * m)) + C * (S * Ty) := by ring
        rw [e]; linarith
      have := (mul_le_mul_iff_of_pos_left (mul_pos hSp hCp)).1 this
      linarith
    rcases le_total 0 μ with hμ0 | hμ0
    · -- T'(m) ≥ T'(μ), then X3
      have x3 := X3 hC hS hCS hH1 hH2 hμ hμ0 ((mul_le_mul_of_nonneg_right hmμ hC).trans hmk)
      have hprod := mul_nonneg (sub_nonneg.2 hmμ) (mul_nonneg hSp.le hc)
      have : S * (2 + 2 * Real.sqrt 2 - 2 * H) ≤ S * ((Tm + (2 - H) * m) + Ty) := by
        rw [hTm']; linarith
      have := (mul_le_mul_iff_of_pos_left hSp).1 this
      linarith
    · -- μ ≤ 0, so Z ≤ 1: T'(m) ≥ k and k + V_y ≥ 1 ≥ Z
      have hZ : 2 + 2 * Real.sqrt 2 - 2 * H ≤ 1 := by
        have : (2 - H) * μ ≤ 0 := mul_nonpos_of_nonneg_of_nonpos (by linarith) hμ0
        linarith
      have hprod := mul_nonneg hm (mul_nonneg hSp.le hc)
      have hprod2 := mul_nonneg hCp.le hc
      have hid2 : S * (C + S) + C * (C + 1 - H) = S + (1 - S) + C * (S + 1 - H) := by linear_combination hCS
      have : S * 1 ≤ S * ((Tm + (2 - H) * m) + Ty) := by
        rw [hTm']; linarith
      have := (mul_le_mul_iff_of_pos_left hSp).1 this
      linarith
  · -- grid branch
    apply min_le_of_left_le
    rcases le_total (μ * C) S with hμk | hμk
    · -- T m ≥ T μ, then X3
      have x3 := X3 hC hS hCS hH1 hH2 hμ (hm.trans hmμ) hμk
      have hprod := mul_nonneg (sub_nonneg.2 hmμ) (mul_nonneg hSp.le (sub_nonneg.2 hS1))
      have : S * 1 ≤ S * (Tm + Ty) := by
        rw [hTm']; linarith
      have := (mul_le_mul_iff_of_pos_left hSp).1 this
      linarith
    · -- tan t ≤ μ: T m ≥ W = V_m - 2β tan ≥ V_m - (Z - 1), then Z1
      have h2b : (2 - H) * S ≤ (1 + 2 * Real.sqrt 2 - 2 * H) * C :=
        calc (2 - H) * S ≤ (2 - H) * (μ * C) := mul_le_mul_of_nonneg_left hμk (by linarith)
          _ = (1 + 2 * Real.sqrt 2 - 2 * H) * C := by rw [← hμ]; ring
      have hWS := mul_le_mul_of_nonneg_left hW0 hSp.le
      have h2bS := mul_le_mul_of_nonneg_left h2b hSp.le
      have : S * C * 1 ≤ S * C * (Tm + Ty) := by
        have e : S * C * (Tm + Ty) = S * (C * Tm) + C * (S * Ty) := by ring
        have hid3 : S * (C * k + S * (S - 1)) + C * (S * k + C * (C + 1 - H)) = C + S + 1 - H - (2 - H) * S ^ 2 := by
          rw [hk]; linear_combination (C + S + 1 - H) * hCS
        rw [e]; linarith
      have := (mul_le_mul_iff_of_pos_left (mul_pos hSp hCp)).1 this
      linarith

/-- **(S)**: every square costs a grid unit or the zigzag amount,
`kap ≥ min (1 + β||m| - |m'||) (Z - β(|m| + |m'|))`. -/
theorem kap_S {C S H m m' : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)
    (hH1 : 2 * Real.sqrt 2 - 1 ≤ H) (hH2 : H < 2) (hfit : C + S ≤ H) :
    min (1 + (1 - H / 2) * |(|m| - |m'|)|) (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) ≤
      kap C S H m m' := by
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  have hfit' : S + C ≤ H := by linarith
  -- both slopes ≥ 0: the first branch, by kap_same
  have same : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → 0 ≤ m → 0 ≤ m' →
      min (1 + (1 - H / 2) * |(|m| - |m'|)|) (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) ≤
        kap C S H m m' := by
    intro C S m m' hC hS hCS hm hm'
    apply min_le_of_left_le
    rw [abs_of_nonneg hm, abs_of_nonneg hm']
    rcases le_total m' m with h | h
    · rw [abs_of_nonneg (by linarith)]; exact kap_same hC hS hCS hm' h
    · rw [abs_of_nonpos (by linarith), kap_symm]
      have := kap_same (H := H) hC hS hCS hm h; linarith
  -- opposite slopes m ≥ 0 ≥ m'
  have opp : ∀ {C S m m' : ℝ}, 0 ≤ C → 0 ≤ S → C ^ 2 + S ^ 2 = 1 → C + S ≤ H → 0 ≤ m → m' ≤ 0 →
      min (1 + (1 - H / 2) * |(|m| - |m'|)|) (2 + 2 * Real.sqrt 2 - 2 * H - (1 - H / 2) * (|m| + |m'|)) ≤
        kap C S H m m' := by
    intro C S m m' hC hS hCS hfit hm hm'
    rw [abs_of_nonneg hm, abs_of_nonpos hm']
    rcases le_total (-m') m with h | h
    · rw [abs_of_nonneg (by linarith)]
      have := kap_S_opp (y := -m') hC hS hCS hH1 hH2 hfit (by linarith) h
      rwa [neg_neg] at this
    · -- mirror: (m, m') ↦ (-m', -m) with the frame swapped
      rw [abs_of_nonpos (by linarith), kap_mirror]
      have := kap_S_opp (m := -m') (y := m) hS hC (by linarith) hH1 hH2 (by linarith) hm h
      rw [show m + -m' = -m' + m by ring, show -(m - -m') = -m' - m by ring]
      exact this
  rcases le_total 0 m with hm | hm <;> rcases le_total 0 m' with hm' | hm'
  · exact same hC hS hCS hm hm'
  · exact opp hC hS hCS hfit hm hm'
  · rw [kap_symm, abs_sub_comm, add_comm |m| |m'|]; exact opp hC hS hCS hfit hm' hm
  · rw [kap_mirror]
    have := same (m := -m') (m' := -m) hS hC hSC (by linarith) (by linarith)
    rwa [abs_neg, abs_neg, abs_sub_comm, add_comm |m'| |m|] at this

end SquarePacking
