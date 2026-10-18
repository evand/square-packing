import Sqpack.RyuStairCount

/-!
# The staircase end: the end cost `≤ (2√2 + 3/2) v² + O(v)`

For a band with `b ≤ v³` squares per row (`v ≥ 10`, `b ≥ v³/2`) and a lift `y ∈ [3v²/4, 3v²/4 + 2]`, the staircase
of `RyuStairCount` with the column tilt `tan(α/2) = (√2 + 1/v)/v` packs the end region with end cost
`E(y) = A₀(y) − N ≤ 4.33 v² + 11 v + 30`, where `A₀(y) = yw + γσ(1 + b²)/2` (Ryu's Proposition 5.2).

* `lb_identity`: `w g₀ + μ w²/2 − lb` as a sum of six terms (an identity of rational functions).
* `term_a`, …, `term_f`: the six terms, bounded from the parameter bounds; `tilt_ok`: the tilt condition.
* `end_bound`: the end packing and its cost.

The leading terms are Ryu's: `w sin α ≈ 2√2 v²` (the triangles beside the columns) and `J₁ + R₀ ≈ 3v²/2` (the ends
of the layers and of the wedge rows).
-/

namespace UnitSquarePacking.Ryu

open UnitSquarePacking.Conjectures

namespace EndP

/-- `area − lb`, with `area = w g₀ + μ w²/2`, `Pw = w − mτ`, `ρ = Pw − Iς`. -/
lemma lb_identity (P : EndP) (hca : P.ca ≠ 0) (hca1 : 1 - P.ca ≠ 0) :
    P.w * P.g0 + P.μ * P.w ^ 2 / 2 - P.lb =
      (P.w - P.m * P.τ) * P.sa + (P.w - P.m * P.τ) * (1 - P.ca) + P.J1 + P.R0 +
      P.μ * ((P.w - P.m * P.τ) * P.m * P.τ + (P.w - P.m * P.τ) * (P.w - P.m * P.τ - P.I * P.ς) +
        (P.w - P.m * P.τ) * P.ς / 2 + P.m ^ 2 * P.τ ^ 2 / 2 - (P.w - P.m * P.τ - P.I * P.ς) ^ 2 / 2 -
        (P.w - P.m * P.τ - P.I * P.ς) * P.ς / 2) +
      (P.w - P.m * P.τ - P.I * P.ς) * (P.g0 - P.sa - P.J1 - (1 - P.ca)) +
      P.τ * (-(P.J1 : ℝ) ^ 2 / 2 + P.J1 / 2 + (P.R0 : ℝ) ^ 2 / 2 - P.R0 * P.m + P.R0 / 2 + P.g0 * P.m) := by
  unfold lb τ ς H0
  field_simp
  ring

end EndP

section terms
variable {v s t : ℝ}

/-- Common facts on `v ≥ 10`, `s = √2`, `t = (s + 1/v)/v`. -/
structure TV (v s t : ℝ) : Prop where
  hv : 10 ≤ v
  hs : s ^ 2 = 2
  hs0 : 0 ≤ s
  htv : t * v = s + 1 / v

lemma TV.v0 (h : TV v s t) : 0 < v := by linarith [h.hv]
lemma TV.s1 (h : TV v s t) : s ≤ 1415 / 1000 := by nlinarith [h.hs, h.hs0]
lemma TV.s2 (h : TV v s t) : 1414 / 1000 ≤ s := by nlinarith [h.hs, h.hs0]
lemma TV.iv (h : TV v s t) : 1 / v ≤ 1 / 10 := by
  rw [div_le_div_iff₀ h.v0 (by norm_num)]; linarith [h.hv]
lemma TV.t0 (h : TV v s t) : 0 < t := by
  have := h.htv; have := h.v0; have : 0 < s + 1 / v := by have := h.s2; positivity
  nlinarith
lemma TV.t1 (h : TV v s t) : t ≤ 1515 / 10000 := by
  have := h.htv; have := h.s1; have := h.iv; have := h.v0; have := h.hv
  nlinarith
lemma TV.v3 (h : TV v s t) : 1000 ≤ v ^ 3 := by
  have := h.hv; have h2 : 100 ≤ v ^ 2 := by nlinarith
  nlinarith

/-- (a) `P sa ≤ 2 s v² + 2 v`. -/
lemma term_a (h : TV v s t) {P sa : ℝ} (hP : P ≤ v ^ 3) (hsa0 : 0 ≤ sa) (hsa : sa ≤ 2 * t) :
    P * sa ≤ 2 * s * v ^ 2 + 2 * v := by
  have h1 : P * sa ≤ v ^ 3 * sa := mul_le_mul_of_nonneg_right hP hsa0
  have h2 : v ^ 3 * sa ≤ v ^ 3 * (2 * t) := mul_le_mul_of_nonneg_left hsa (by have := h.v0; positivity)
  have h3 : v ^ 3 * (2 * t) = 2 * s * v ^ 2 + 2 * v := by
    have e : v ^ 3 * (2 * t) = 2 * v ^ 2 * (t * v) := by ring
    rw [e, h.htv]; have := h.v0; field_simp
  linarith

/-- (b) `P (1 − ca) ≤ 4 v + 4 s + 2/v`. -/
lemma term_b (h : TV v s t) {P ca : ℝ} (hP : P ≤ v ^ 3) (hca0 : 0 ≤ 1 - ca) (hca : 1 - ca ≤ 2 * t ^ 2) :
    P * (1 - ca) ≤ 4 * v + 4 * s + 2 / v := by
  have h1 : P * (1 - ca) ≤ v ^ 3 * (1 - ca) := mul_le_mul_of_nonneg_right hP hca0
  have h2 : v ^ 3 * (1 - ca) ≤ v ^ 3 * (2 * t ^ 2) := mul_le_mul_of_nonneg_left hca (by have := h.v0; positivity)
  have h3 : v ^ 3 * (2 * t ^ 2) = 4 * v + 4 * s + 2 / v := by
    have e : v ^ 3 * (2 * t ^ 2) = 2 * v * (t * v) ^ 2 := by ring
    rw [e, h.htv]; have := h.v0; field_simp; linear_combination (2 * v ^ 2) * h.hs
  linarith

/-- `X τ ≤ 2.177 v + 3.1` for `0 ≤ X ≤ 3v²/4 + 5.001`. -/
lemma term_Xτ (h : TV v s t) {X τ : ℝ} (hX0 : 0 ≤ X) (hX : X ≤ 3 / 4 * v ^ 2 + 5 + 1 / 1000) (hτ0 : 0 ≤ τ)
    (hτ : τ ≤ 41 / 20 * t) : X * τ ≤ 2177 / 1000 * v + 31 / 10 := by
  have h1 : X * τ ≤ (3 / 4 * v ^ 2 + 5 + 1 / 1000) * (41 / 20 * t) :=
    mul_le_mul hX hτ hτ0 (by have := h.v0; positivity)
  have h2 : (3 / 4 * v ^ 2 + 5 + 1 / 1000) * (41 / 20 * t) =
      41 / 20 * (3 / 4 * v * (t * v) + (5 + 1 / 1000) * t) := by ring
  have h4 : 3 / 4 * v * (s + 1 / v) = 3 / 4 * s * v + 3 / 4 := by have := h.v0; field_simp
  rw [h2, h.htv, h4] at h1
  have := h.t1; have := h.s1; have := h.v0
  nlinarith

/-- (d) the `μ` terms. -/
lemma term_d (h : TV v s t) {μ w m τ ς ρ P : ℝ} (hμ0 : 0 ≤ μ) (hμ : μ ≤ 5 / v ^ 3) (hΔ2 : w * μ ≤ 2 + 1 / 1000)
    (hw1 : 0 ≤ w) (hP0 : 0 ≤ P) (hPw : P ≤ w) (hρ0 : 0 ≤ ρ) (hρ : ρ ≤ ς) (hς : ς ≤ 21 / 20)
    (hmτ0 : 0 ≤ m * τ) (hmτ : m * τ ≤ 2177 / 1000 * v + 31 / 10) :
    μ * (P * m * τ + P * ρ + P * ς / 2 + m ^ 2 * τ ^ 2 / 2 - ρ ^ 2 / 2 - ρ * ς / 2) ≤ 436 / 100 * v + 12 := by
  have hv0 := h.v0
  have hv := h.hv
  have h1 : P * m * τ + P * ρ + P * ς / 2 + m ^ 2 * τ ^ 2 / 2 - ρ ^ 2 / 2 - ρ * ς / 2 ≤
      w * (m * τ) + w * ς + w * ς / 2 + (m * τ) ^ 2 / 2 := by
    have e1 : P * (m * τ) ≤ w * (m * τ) := mul_le_mul_of_nonneg_right hPw hmτ0
    have e2 : P * ρ ≤ w * ς := mul_le_mul hPw hρ hρ0 hw1
    have e3 : P * ς ≤ w * ς := mul_le_mul_of_nonneg_right hPw (by linarith)
    have e4 : 0 ≤ ρ ^ 2 / 2 + ρ * ς / 2 := by have := mul_nonneg hρ0 (by linarith : (0 : ℝ) ≤ ς); positivity
    have e5 : m ^ 2 * τ ^ 2 = (m * τ) ^ 2 := by ring
    have e6 : P * m * τ = P * (m * τ) := by ring
    rw [e5, e6]; linarith
  have h2 := mul_le_mul_of_nonneg_left h1 hμ0
  have h3 : μ * (w * (m * τ) + w * ς + w * ς / 2 + (m * τ) ^ 2 / 2) =
      (w * μ) * (m * τ) + (w * μ) * (3 / 2 * ς) + μ * (m * τ) ^ 2 / 2 := by ring
  have h4 : (w * μ) * (m * τ) ≤ (2 + 1 / 1000) * (2177 / 1000 * v + 31 / 10) :=
    mul_le_mul hΔ2 hmτ hmτ0 (by norm_num)
  have h5 : (w * μ) * (3 / 2 * ς) ≤ (2 + 1 / 1000) * (3 / 2 * (21 / 20)) :=
    mul_le_mul hΔ2 (by linarith) (by nlinarith) (by norm_num)
  have h6 : μ * (m * τ) ^ 2 / 2 ≤ 2 := by
    have e1 : (m * τ) ^ 2 ≤ (2177 / 1000 * v + 31 / 10) ^ 2 := pow_le_pow_left₀ hmτ0 hmτ 2
    have e2 : μ * (m * τ) ^ 2 ≤ 5 / v ^ 3 * (2177 / 1000 * v + 31 / 10) ^ 2 :=
      mul_le_mul hμ e1 (sq_nonneg _) (by positivity)
    have e3 : 5 / v ^ 3 * (2177 / 1000 * v + 31 / 10) ^ 2 ≤ 4 := by
      rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
      nlinarith [mul_nonneg (by linarith : (0 : ℝ) ≤ v - 10) (sq_nonneg v)]
    linarith
  nlinarith

/-- (f) the `τ` term, from `J₁ ≥ g₀ + c`, `R₀ ∈ [g₀ − 1, g₀]`. -/
lemma term_f (h : TV v s t) {τ J1 R0 m g0 Δ ε : ℝ} (hτ0 : 0 ≤ τ) (hτ : τ ≤ 41 / 20 * t)
    (hΔ1 : 2 ≤ Δ) (hΔ2 : Δ ≤ 2 + 1 / 1000) (hε0 : 0 ≤ ε) (hε : ε ≤ 1 / 2) (hgε : g0 * ε ≤ 11)
    (hg0 : 0 ≤ g0) (hJl : g0 + (Δ - 1 - ε) ≤ J1) (hm : m ≤ g0 + Δ + 1) (hmg : g0 < m)
    (hR1 : g0 - 1 ≤ R0) (hR2 : R0 ≤ g0) :
    τ * (-J1 ^ 2 / 2 + J1 / 2 + R0 ^ 2 / 2 - R0 * m + R0 / 2 + g0 * m) ≤ 5 := by
  set c := Δ - 1 - ε with hcd
  have hc1 : 1 / 2 ≤ c := by rw [hcd]; linarith
  have q1 : -J1 ^ 2 / 2 + J1 / 2 ≤ -(g0 + c) ^ 2 / 2 + (g0 + c) / 2 := by
    have := mul_nonneg (by linarith : (0 : ℝ) ≤ J1 - (g0 + c)) (by linarith : (0 : ℝ) ≤ J1 + (g0 + c) - 1)
    nlinarith
  have q2 : R0 ^ 2 / 2 - R0 * m + R0 / 2 ≤ (g0 - 1) ^ 2 / 2 - (g0 - 1) * m + (g0 - 1) / 2 := by
    have := mul_nonneg (by linarith : (0 : ℝ) ≤ R0 - (g0 - 1)) (by linarith : (0 : ℝ) ≤ 2 * m - R0 - g0)
    nlinarith
  have hQ : -J1 ^ 2 / 2 + J1 / 2 + R0 ^ 2 / 2 - R0 * m + R0 / 2 + g0 * m ≤ 15 := by
    have e : -(g0 + c) ^ 2 / 2 + (g0 + c) / 2 + ((g0 - 1) ^ 2 / 2 - (g0 - 1) * m + (g0 - 1) / 2) + g0 * m =
        m - c * g0 - c * (c - 1) / 2 := by ring
    have h1 : 1 - c ≤ ε := by rw [hcd]; linarith
    have h2 : g0 * (1 - c) ≤ g0 * ε := mul_le_mul_of_nonneg_left h1 hg0
    nlinarith [sq_nonneg (c - 1 / 2)]
  have hτ' : τ * 15 ≤ 5 := by
    have := h.t1; nlinarith
  rcases le_total 0 (-J1 ^ 2 / 2 + J1 / 2 + R0 ^ 2 / 2 - R0 * m + R0 / 2 + g0 * m) with hq | hq
  · nlinarith [mul_le_mul_of_nonneg_left hQ hτ0]
  · nlinarith [mul_nonpos_of_nonneg_of_nonpos hτ0 hq]

/-- **The tilt condition** `(G + Δ + 1) ca + sa ≤ G` at `tan(α/2) = t`, for `G ≥ 3v²/4 − 1/v³`, `Δ ≤ 2 + 1/v³`. -/
lemma tilt_ok (h : TV v s t) {G Δ : ℝ} (hG : 3 / 4 * v ^ 2 - 1 / v ^ 3 ≤ G) (hΔ0 : 0 ≤ Δ)
    (hΔ : Δ ≤ 2 + 1 / v ^ 3) :
    (G + Δ + 1) * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2) ≤ G := by
  have hv0 := h.v0
  have ht0 := h.t0
  have ht1 := h.t1
  have hs2 := h.s2
  have hp : 0 < 1 + t ^ 2 := by positivity
  have hq : 0 ≤ 1 - t ^ 2 := by nlinarith
  set a := 1 / v with had
  have ha0 : 0 < a := by positivity
  have ha1 : a ≤ 1 / 10 := h.iv
  have hva : v * a = 1 := by rw [had]; field_simp
  have hta : t = (s + a) * a := by
    have := h.htv; rw [← had] at this
    have e : t = t * v * a := by rw [mul_assoc, hva, mul_one]
    rw [e, this]
  have ha3 : 1 / v ^ 3 = a ^ 3 := by rw [had]; field_simp
  rw [ha3] at hG hΔ
  have hGt : 3 / 4 * (s + a) ^ 2 - a ^ 3 * t ^ 2 ≤ G * t ^ 2 := by
    have e1 := mul_le_mul_of_nonneg_right hG (sq_nonneg t)
    have e2 : v ^ 2 * t ^ 2 = (s + a) ^ 2 := by
      have := h.htv; rw [← had] at this; rw [← this]; ring
    nlinarith
  have e : (G + Δ + 1) * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2) =
      ((G + Δ + 1) * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) := by field_simp
  rw [e, div_le_iff₀ hp]
  have k1 : (Δ + 1) * (1 - t ^ 2) ≤ (3 + a ^ 3) * (1 - t ^ 2) := mul_le_mul_of_nonneg_right (by linarith) hq
  have k2 : a ^ 3 * (1 - t ^ 2) ≤ a ^ 3 := by nlinarith [pow_pos ha0 3]
  have k3 : (s + a) ^ 2 = 2 + 2 * s * a + a ^ 2 := by nlinarith [h.hs]
  have k4 : a ^ 3 + 1 / 2 * a ^ 2 ≤ s * a := by nlinarith
  have k5 : 2 * a ^ 3 * t ^ 2 ≤ 3 * t ^ 2 := by nlinarith [sq_nonneg t, pow_pos ha0 3]
  rw [hta] at k1 k2 k5 hGt ⊢
  nlinarith

end terms


set_option maxHeartbeats 4000000 in
/-- **The staircase end.**  For a band side with `b ≤ v³` squares per row (`v ≥ 10`) and a lift
`y ∈ [3v²/4, 3v²/4 + 2]`: a packing of `endReg γ σ w y` with end cost `A₀(y) − n ≤ 4.33 v² + 11 v + 30`. -/
theorem end_bound {b : ℕ} {v γ σ w y : ℝ} (hv : 10 ≤ v) (hγ : 0 < γ) (hσ : 0 < σ) (hu : γ ^ 2 + σ ^ 2 = 1)
    (hwd : w = b * γ + σ) (hΔ1 : 2 ≤ w * (σ / γ)) (hΔ2 : w * (σ / γ) ≤ 2 + 1 / v ^ 3)
    (hμ : σ / γ ≤ 5 / v ^ 3) (hσμ : σ * (σ / γ) ≤ 1 / v ^ 3) (hw1 : v ^ 3 / 3 ≤ w) (hw2 : w ≤ v ^ 3)
    (hy1 : 3 / 4 * v ^ 2 ≤ y) (hy2 : y ≤ 3 / 4 * v ^ 2 + 2) :
    ∃ n : ℕ, PacksIn n (endReg γ σ w y) ∧ y * w + γ * σ * (1 + b ^ 2) / 2 - n ≤ 433 / 100 * v ^ 2 + 11 * v + 30 := by
  have hv0 : 0 < v := by linarith
  set s := Real.sqrt 2 with hsd
  have hs : s ^ 2 = 2 := Real.sq_sqrt (by norm_num)
  have hs0 : 0 ≤ s := Real.sqrt_nonneg 2
  have hs1 : s ≤ 1415 / 1000 := by nlinarith only [hs, hs0]
  have hs2 : 1414 / 1000 ≤ s := by nlinarith only [hs, hs0]
  have hiv : 1 / v ≤ 1 / 10 := by rw [div_le_div_iff₀ hv0 (by norm_num)]; linarith
  have hiv0 : 0 < 1 / v := by positivity
  have hv21 : 10 * v ≤ v ^ 2 := by nlinarith only [hv]
  have hv2 : 100 ≤ v ^ 2 := by nlinarith only [hv]
  have hv32 : 10 * v ^ 2 ≤ v ^ 3 := by nlinarith only [hv, hv2]
  have hv3 : 1000 ≤ v ^ 3 := by linarith only [hv2, hv32]
  have hv3v : 100 * v ≤ v ^ 3 := by linarith only [hv21, hv32]
  have hiv3 : 1 / v ^ 3 ≤ 1 / 1000 := by rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith
  have hσμ0 : 0 ≤ σ * (σ / γ) := by positivity
  set t := (s + 1 / v) / v with ht
  have htv : t * v = s + 1 / v := by rw [ht]; field_simp
  have ht0 : 0 < t := by positivity
  have ht1 : t ≤ 1515 / 10000 := by rw [ht, div_le_iff₀ hv0]; linarith
  have ht2 : t ^ 2 ≤ 1 / 41 := by nlinarith only [ht0, ht1]
  have hp : 0 < 1 + t ^ 2 := by positivity
  have hq : 0 < 1 - t ^ 2 := by linarith only [ht2]
  set ca := (1 - t ^ 2) / (1 + t ^ 2) with hcad
  set sa := 2 * t / (1 + t ^ 2) with hsad
  have hca : 0 < ca := div_pos hq hp
  have hca1 : ca < 1 := by rw [hcad, div_lt_one hp]; nlinarith only [ht0]
  have hsa : 0 < sa := by positivity
  have hu' : ca ^ 2 + sa ^ 2 = 1 := by rw [hcad, hsad]; field_simp; ring
  set P : EndP := ⟨γ, σ, w, y, ca, sa⟩ with hPd
  have hτ : P.τ = 2 * t / (1 - t ^ 2) := by
    show sa / ca = _; rw [hcad, hsad]; field_simp
  have hς : P.ς = (1 + t ^ 2) / (1 - t ^ 2) := by
    show 1 / ca = _; rw [hcad]; field_simp
  have hτ0 : 0 ≤ P.τ := by rw [hτ]; positivity
  have hτb : P.τ ≤ 41 / 20 * t := by
    rw [hτ, div_le_iff₀ hq]; nlinarith only [ht2, ht0]
  have hς1 : 1 ≤ P.ς := by rw [hς, le_div_iff₀ hq]; nlinarith only [ht0]
  have hςb : P.ς ≤ 21 / 20 := by rw [hς, div_le_iff₀ hq]; linarith only [ht2]
  have hsab : sa ≤ 2 * t := by rw [hsad, div_le_iff₀ hp]; nlinarith only [ht0]
  have hcab : 1 - ca ≤ 2 * t ^ 2 := by
    have : 1 - ca = 2 * t ^ 2 / (1 + t ^ 2) := by rw [hcad]; field_simp; ring
    rw [this, div_le_iff₀ hp]; nlinarith only [ht0]
  have hμd : P.μ = σ / γ := rfl
  have hμ0 : 0 ≤ P.μ := div_nonneg hσ.le hγ.le
  have hg0 : P.g0 = y - σ * (σ / γ) := by show y + (0 - σ) * (σ / γ) = _; ring
  have hgb : P.gb = P.g0 + w * (σ / γ) := by
    show y + (w - σ) * (σ / γ) = y + (0 - σ) * (σ / γ) + w * (σ / γ); ring
  have hgb0 : 0 ≤ P.gb := by rw [hgb, hg0]; linarith only [hy1, hσμ, hiv3, hΔ1, hv2]
  have hmR : (P.m : ℝ) ≤ P.gb + 1 := by
    show ((⌊P.gb⌋₊ + 1 : ℕ) : ℝ) ≤ _; push_cast; linarith [Nat.floor_le hgb0]
  have hm0 : (0 : ℝ) ≤ P.m := Nat.cast_nonneg _
  have hmX : (P.m : ℝ) ≤ 3 / 4 * v ^ 2 + 5 + 1 / 1000 := by
    rw [hgb, hg0] at hmR; linarith only [hmR, hy2, hΔ2, hσμ0, hiv3]
  -- the tilt condition `m ca + sa ≤ g₀`
  have htilt : P.H0 ≤ P.g0 := by
    show (P.m : ℝ) * ca + sa ≤ P.g0
    have e1 : (P.m : ℝ) * ca ≤ (P.g0 + w * (σ / γ) + 1) * ca := by
      apply mul_le_mul_of_nonneg_right _ hca.le; rw [hgb] at hmR; exact hmR
    have key := tilt_ok ⟨hv, hs, hs0, htv⟩ (G := P.g0) (Δ := w * (σ / γ)) (by rw [hg0]; linarith)
      (by linarith) hΔ2
    rw [← hcad, ← hsad] at key
    linarith only [e1, key]
  have hT : TV v s t := ⟨hv, hs, hs0, htv⟩
  -- `mτ ≤ w`
  have hmτ' := term_Xτ hT hm0 hmX hτ0 hτb
  have hwmτ : (P.m : ℝ) * P.τ ≤ w := by linarith only [hmτ', hw1, hv3v, hv]
  have h : P.Hyp := ⟨hγ, hσ.le, hca, hca1, hsa, hu', htilt, hwmτ⟩
  refine ⟨P.N, EndP.packs h, ?_⟩
  have hc := EndP.count h
  have hid := EndP.lb_identity P hca.ne' (by linarith)
  -- `A₀ − area = μ/2`
  have hA : y * w + γ * σ * (1 + b ^ 2) / 2 - (w * P.g0 + P.μ * w ^ 2 / 2) = P.μ / 2 := by
    have key : y * w + γ * σ * (1 + b ^ 2) / 2 - (w * (y - σ * (σ / γ)) + σ / γ * w ^ 2 / 2) - σ / γ / 2 =
        σ * (γ ^ 2 + σ ^ 2 - 1) / (2 * γ) := by rw [hwd]; field_simp; ring
    rw [hu, sub_self, mul_zero, zero_div, sub_eq_zero] at key
    rw [hg0, hμd]; exact key
  -- the parameter bounds
  have hP0 : 0 ≤ P.w - P.m * P.τ := by show 0 ≤ w - _; linarith
  have hI2 : P.w - P.m * P.τ - P.I * P.ς ≤ P.ς := by
    have hς0 : 0 < P.ς := by linarith
    have := Nat.lt_floor_add_one ((P.w - P.m * P.τ) / P.ς)
    rw [div_lt_iff₀ hς0] at this
    have e : (P.I : ℝ) = ⌊(P.w - P.m * P.τ) / P.ς⌋₊ := rfl
    rw [e]; linarith only [this]
  have hI1 : 0 ≤ P.w - P.m * P.τ - P.I * P.ς := by
    have hς0 : 0 < P.ς := by linarith
    have := Nat.floor_le (div_nonneg hP0 hς0.le)
    have e := mul_le_mul_of_nonneg_right this hς0.le
    rw [div_mul_cancel₀ _ hς0.ne'] at e
    have e2 : (P.I : ℝ) = ⌊(P.w - P.m * P.τ) / P.ς⌋₊ := rfl
    rw [e2]; linarith only [e]
  have hJ1a : (P.J1 : ℝ) ≤ P.g0 + w * (σ / γ) := by
    have := EndP.J1_lt_m h
    have : P.J1 ≤ ⌊P.gb⌋₊ := by have e : P.m = ⌊P.gb⌋₊ + 1 := rfl; omega
    have : (P.J1 : ℝ) ≤ ⌊P.gb⌋₊ := by exact_mod_cast this
    rw [← hgb]; linarith [Nat.floor_le hgb0]
  have hJ1b : P.g0 + w * (σ / γ) - 1 - σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ) ≤ P.J1 := by
    have h1 : (⌊P.g (P.xi P.I - P.gb * P.τ)⌋₊ : ℝ) ≤ P.J1 := by
      have : ⌊P.g (P.xi P.I - P.gb * P.τ)⌋₊ ≤ P.J1 := le_max_right _ _
      exact_mod_cast this
    have h2 := Nat.lt_floor_add_one (P.g (P.xi P.I - P.gb * P.τ))
    have h3 : P.g (P.xi P.I - P.gb * P.τ) = P.gb - (P.w - P.xi P.I + P.gb * P.τ) * (σ / γ) := by
      show y + (P.xi P.I - P.gb * P.τ - σ) * (σ / γ) = y + (w - σ) * (σ / γ) - _; ring
    have h4 : P.w - P.xi P.I ≤ P.ς := by unfold EndP.xi; linarith
    have h5 : (P.w - P.xi P.I + P.gb * P.τ) * (σ / γ) ≤ (P.ς + P.gb * P.τ) * (σ / γ) :=
      mul_le_mul_of_nonneg_right (by linarith) hμ0
    rw [← hgb]
    linarith only [h1, h2, h3, h5]
  have hR1 : P.g0 - 1 ≤ P.R0 := by have := Nat.lt_floor_add_one P.g0; show _ ≤ ((⌊P.g0⌋₊ : ℕ) : ℝ); linarith
  have hg00 : 0 ≤ P.g0 := by rw [hg0]; linarith only [hy1, hσμ, hiv3, hv2]
  have hR2 : (P.R0 : ℝ) ≤ P.g0 := Nat.floor_le hg00
  -- the six terms
  have hPw : P.w - P.m * P.τ ≤ v ^ 3 := by show w - _ ≤ _; linarith only [hw2, mul_nonneg hm0 hτ0]
  have ta := term_a hT hPw hsa.le hsab
  have tb := term_b hT hPw (by linarith) hcab
  have hΔ2' : w * (σ / γ) ≤ 2 + 1 / 1000 := by linarith only [hΔ2, hiv3]
  have td := term_d hT hμ0 hμ hΔ2' (by linarith) hP0 (by show w - _ ≤ w; linarith only [mul_nonneg hm0 hτ0])
    hI1 hI2 hςb (mul_nonneg hm0 hτ0) hmτ'
  have hgb0' : 0 ≤ P.g0 + w * (σ / γ) := by linarith
  have hgbX : P.g0 + w * (σ / γ) ≤ 3 / 4 * v ^ 2 + 5 + 1 / 1000 := by
    rw [hg0]; linarith only [hy2, hσμ0, hΔ2, hiv3]
  have hgbτ := term_Xτ hT hgb0' hgbX hτ0 hτb
  have hgτ0 := mul_nonneg hgb0' hτ0
  have hε0 : 0 ≤ σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ) :=
    mul_nonneg hμ0 (by linarith only [hς1, hgτ0])
  have hεb : σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ) ≤ 5 / v ^ 3 * (21 / 20 + (2177 / 1000 * v + 31 / 10)) :=
    mul_le_mul hμ (by linarith only [hςb, hgbτ]) (by linarith only [hς1, hgτ0]) (by positivity)
  have hε : σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ) ≤ 1 / 2 := by
    have : 5 / v ^ 3 * (21 / 20 + (2177 / 1000 * v + 31 / 10)) ≤ 1 / 2 := by
      rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]; linarith only [hv3v, hv]
    linarith
  have hgε : P.g0 * (σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ)) ≤ 11 := by
    have h1 : P.g0 * (σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ)) ≤
        (3 / 4 * v ^ 2 + 2) * (5 / v ^ 3 * (21 / 20 + (2177 / 1000 * v + 31 / 10))) :=
      mul_le_mul (by rw [hg0]; linarith only [hy2, hσμ0]) hεb hε0 (by positivity)
    have h2 : (3 / 4 * v ^ 2 + 2) * (5 / v ^ 3 * (21 / 20 + (2177 / 1000 * v + 31 / 10))) ≤ 11 := by
      rw [show (3 / 4 * v ^ 2 + 2) * (5 / v ^ 3 * (21 / 20 + (2177 / 1000 * v + 31 / 10))) =
        (3 / 4 * v ^ 2 + 2) * (5 * (21 / 20 + (2177 / 1000 * v + 31 / 10))) / v ^ 3 by ring,
        div_le_iff₀ (by positivity)]
      nlinarith only [hv32, hv21, hv, hv2]
    linarith
  have hmg : P.g0 < P.m := by
    have := Nat.lt_floor_add_one P.gb
    have e : (P.m : ℝ) = ⌊P.gb⌋₊ + 1 := by show ((⌊P.gb⌋₊ + 1 : ℕ) : ℝ) = _; push_cast; ring
    rw [e]; linarith only [this, hgb, hΔ1]
  have hmR' : (P.m : ℝ) ≤ P.g0 + w * (σ / γ) + 1 := by rw [hgb] at hmR; exact hmR
  have tf := term_f (τ := P.τ) (J1 := P.J1) (R0 := P.R0) (m := P.m) (g0 := P.g0) (Δ := w * (σ / γ))
    (ε := σ / γ * (P.ς + (P.g0 + w * (σ / γ)) * P.τ)) hT hτ0 hτb hΔ1 hΔ2' hε0 hε hgε hg00
    (by linarith only [hJ1b]) hmR' hmg hR1 hR2
  have te : (P.w - P.m * P.τ - P.I * P.ς) * (P.g0 - sa - P.J1 - (1 - ca)) ≤ 0 :=
    mul_nonpos_of_nonneg_of_nonpos hI1 (by linarith only [hJ1b, hΔ1, hε, hsa, hca1])
  have hμ2 : σ / γ / 2 ≤ 1 / 100 := by
    have : 5 / v ^ 3 ≤ 5 / 1000 := by rw [div_le_div_iff₀ (by positivity) (by norm_num)]; linarith only [hv3]
    linarith
  have h2v : 2 / v ≤ 1 / 5 := by rw [div_le_div_iff₀ hv0 (by norm_num)]; linarith
  have hTc : (P.J1 : ℝ) + P.R0 ≤ 3 / 2 * v ^ 2 + 6 + 1 / 1000 := by
    rw [hg0] at hJ1a hR2; linarith only [hJ1a, hR2, hy2, hσμ0, hΔ2, hiv3]
  -- assemble
  have hPca : P.ca = ca := rfl
  have hPsa : P.sa = sa := rfl
  rw [hPca, hPsa] at hid
  rw [hμd] at hA hid
  have hs15 : 2 * s + 3 / 2 ≤ 433 / 100 := by linarith
  have hlead : (2 * s + 3 / 2) * v ^ 2 ≤ 433 / 100 * v ^ 2 := mul_le_mul_of_nonneg_right hs15 (sq_nonneg v)
  have hN : (P.lb : ℝ) ≤ P.N := hc
  have hA' : y * P.w + γ * σ * (1 + b ^ 2) / 2 - (P.w * P.g0 + σ / γ * P.w ^ 2 / 2) = σ / γ / 2 := hA
  -- `A₀ − N ≤ μ/2 + (area − lb)`
  have key : y * P.w + γ * σ * (1 + b ^ 2) / 2 - P.N ≤ σ / γ / 2 + (P.w * P.g0 + σ / γ * P.w ^ 2 / 2 - P.lb) := by
    linarith only [hA', hN]
  rw [hid] at key
  have d1 : (P.w - P.m * P.τ) * sa + (P.w - P.m * P.τ) * (1 - ca) ≤ 2 * s * v ^ 2 + 2 * v + (4 * v + 4 * s + 2 / v) := by
    linarith only [ta, tb]
  have d2 : (P.J1 : ℝ) + P.R0 ≤ 3 / 2 * v ^ 2 + 6 + 1 / 1000 := hTc
  rw [hμd] at td
  show y * P.w + γ * σ * (1 + b ^ 2) / 2 - (P.N : ℝ) ≤ _
  linarith only [key, d1, d2, td, te, tf, hμ2, h2v, hlead, hs1, hv0]

end UnitSquarePacking.Ryu
