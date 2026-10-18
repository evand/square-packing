import Sqpack.RyuStairEnd
import Sqpack.RyuLayers

/-!
# Staircase ends: `c*(k) = O(k^{2/5})` (after Ryu's Theorem 1.2)

S. Ryu, *k² − c upper bounds* (squarepacker/k2-minus-c-upper v1.1; jlevy/squares#471, #486), Theorem 1.2:
`c*(k) < 20.668 k^{2/5} + 0.623` for `k ≥ 1.6 · 10⁷`.  Here, with the staircase ends of `RyuStairEnd`:

* `packs_stair`: for `v ≥ 10`, `v³/2 ≤ b ≤ v³` and `k ≥ b + 2v² + 10`, `k² − c` unit squares pack in a square of
  side `< k` whenever `c ≥ 4k/b + 4(4.33 v² + 11 v + 30)`.  The L-shape of `RyuBand.lshape`, both bands lifted by
  `3v²/4` (`RyuEnds.band_ends`), all four ends packed by staircases.
* `packs_stair_k25`: `c*(k) < 18.9 k^{2/5} + 35.7 k^{1/5} + 122` for every `k ≥ 3 · 10⁵`.
* `packs_ryu_thm12`: Ryu's `c*(k) < 20.668 k^{2/5} + 0.623`, for every `k ≥ 8 · 10⁶`.

The count: `k² − N = (k² − S²) + Σ_bands (R − 1) tan θ + Σ_ends E(y)` (Ryu's Proposition 5.2, `ident_band`), with
`E(y) ≤ 4.33 v² + 11 v + 30` at every end, `(R − 1) tan θ ≤ Lσ`, and a band angle `tan(θ/2) = 1/b + 1/(8kb³)`.
-/

namespace UnitSquarePacking.Ryu

open UnitSquarePacking.Conjectures

/-- The end-cost bound. -/
noncomputable def Eub (v : ℝ) : ℝ := 433 / 100 * v ^ 2 + 11 * v + 30

/-! ## The band at `tan(θ/2) = t ∈ (1/b, 1/b + 1/(8b⁴)]` -/

section band
variable {b t v : ℝ}

/-- The facts on the band that `end_bound` and the count use. -/
structure BandFacts (b t v : ℝ) : Prop where
  hγ : 0 < (1 - t ^ 2) / (1 + t ^ 2)
  hσ : 0 < 2 * t / (1 + t ^ 2)
  hσt : 2 * t / (1 + t ^ 2) ≤ 2 * t
  hiγ : 1 / ((1 - t ^ 2) / (1 + t ^ 2)) ≤ 2
  hh : b * (2 * t / (1 + t ^ 2)) + (1 - t ^ 2) / (1 + t ^ 2) ≤ 4
  hΔ1 : 2 ≤ (b * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2)) *
    (2 * t / (1 + t ^ 2) / ((1 - t ^ 2) / (1 + t ^ 2)))
  hΔ2 : (b * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2)) *
    (2 * t / (1 + t ^ 2) / ((1 - t ^ 2) / (1 + t ^ 2))) ≤ 2 + 1 / v ^ 3
  hμ : 2 * t / (1 + t ^ 2) / ((1 - t ^ 2) / (1 + t ^ 2)) ≤ 5 / v ^ 3
  hσμ : 2 * t / (1 + t ^ 2) * (2 * t / (1 + t ^ 2) / ((1 - t ^ 2) / (1 + t ^ 2))) ≤ 1 / v ^ 3
  hw1 : v ^ 3 / 3 ≤ b * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2)
  hw2 : b * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2) ≤ v ^ 3

set_option maxHeartbeats 1000000 in
theorem bandFacts (hb : 500 ≤ b) (ht0 : 1 / b < t) (hbt : b * t ≤ 1 + 1 / (8 * b ^ 2)) (hv : 10 ≤ v)
    (hb1 : b ≤ v ^ 3) (hb2 : v ^ 3 ≤ 2 * b) : BandFacts b t v := by
  have hb0 : 0 < b := by linarith
  have hv2 : 100 ≤ v ^ 2 := by nlinarith only [hv]
  have hv3 : 1000 ≤ v ^ 3 := by nlinarith only [hv, hv2]
  have hv3p : 0 < v ^ 3 := by linarith
  have htpos : 0 < t := lt_trans (by positivity) ht0
  have hbt1 : 1 ≤ b * t := by rw [div_lt_iff₀ hb0] at ht0; linarith
  have h8 : 1 / (8 * b ^ 2) ≤ 1 / 1000 := by rw [div_le_div_iff₀ (by positivity) (by norm_num)]; nlinarith only [hb]
  have hbt2 : b * t ≤ 1001 / 1000 := by linarith
  have ht1 : t ≤ 1 / 250 := by nlinarith only [hbt2, hb, htpos]
  have ht2 : t ^ 2 ≤ 1 / 62500 := by nlinarith only [ht1, htpos]
  have hp : 0 < 1 + t ^ 2 := by positivity
  have hq : 0 < 1 - t ^ 2 := by linarith only [ht2]
  -- closed forms
  have eμ : 2 * t / (1 + t ^ 2) / ((1 - t ^ 2) / (1 + t ^ 2)) = 2 * t / (1 - t ^ 2) := by field_simp
  have ew : b * ((1 - t ^ 2) / (1 + t ^ 2)) + 2 * t / (1 + t ^ 2) = (b * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) := by
    field_simp
  have eΔ : (b * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) * (2 * t / (1 - t ^ 2)) =
      (2 * b * t * (1 - t ^ 2) + 4 * t ^ 2) / ((1 + t ^ 2) * (1 - t ^ 2)) := by field_simp; ring
  have hbt3 : b * t ^ 2 ≤ 1001 / 1000 * t := by nlinarith only [hbt2, htpos]
  have hiv3 : 1 / v ^ 3 ≥ 1 / (2 * b) := by rw [ge_iff_le, div_le_div_iff₀ (by positivity) hv3p]; linarith
  have hμb : 2 * t / (1 - t ^ 2) ≤ 201 / 100 / b := by
    rw [div_le_div_iff₀ hq hb0]; nlinarith only [hbt2, ht2]
  have hib : 1 / b ≤ 2 / v ^ 3 := by rw [div_le_div_iff₀ hb0 hv3p]; linarith
  have hwb : (b * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) ≤ b := by
    rw [div_le_iff₀ hp]; nlinarith only [hbt1, htpos, hb]
  have hwb' : b - 1 ≤ (b * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) := by
    rw [le_div_iff₀ hp]; nlinarith only [hbt1, hbt3, htpos, hb, ht2]
  refine ⟨div_pos hq hp, by positivity, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · rw [div_le_iff₀ hp]; nlinarith only [htpos]
  · rw [one_div_div, div_le_iff₀ hq]; nlinarith only [ht2]
  · have : b * (2 * t / (1 + t ^ 2)) ≤ 2 * (1001 / 1000) := by
      rw [mul_div_assoc', div_le_iff₀ hp]; nlinarith only [hbt2, htpos]
    have : (1 - t ^ 2) / (1 + t ^ 2) ≤ 1 := by rw [div_le_one hp]; nlinarith only [sq_nonneg t]
    linarith
  · rw [eμ, ew, eΔ, le_div_iff₀ (by positivity)]
    nlinarith only [hbt1, hbt2, htpos, ht2, mul_nonneg (by linarith : (0 : ℝ) ≤ b * t - 1) (sq_nonneg t)]
  · rw [eμ, ew]
    have h1 : (b * (1 - t ^ 2) + 2 * t) / (1 + t ^ 2) * (2 * t / (1 - t ^ 2)) ≤ b * (2 * t / (1 - t ^ 2)) :=
      mul_le_mul_of_nonneg_right hwb (by positivity)
    have h2 : b * (2 * t / (1 - t ^ 2)) ≤ 2 + 1 / (2 * b) := by
      rw [mul_div_assoc', div_le_iff₀ hq]
      have e : (2 + 1 / (2 * b)) * (1 - t ^ 2) = 2 - 2 * t ^ 2 + 1 / (2 * b) - t ^ 2 / (2 * b) := by ring
      rw [e]
      have h3 : t ^ 2 / (2 * b) ≤ t ^ 2 := by rw [div_le_iff₀ (by positivity)]; nlinarith only [hb, sq_nonneg t]
      have h4 : 1 / (8 * b ^ 2) * 2 + 3 * t ^ 2 ≤ 1 / (2 * b) := by
        have : t ^ 2 ≤ (1001 / 1000) ^ 2 / b ^ 2 := by rw [le_div_iff₀ (by positivity)]; nlinarith only [hbt2, htpos, hbt1]
        have e2 : 1 / (8 * b ^ 2) * 2 + 3 * ((1001 / 1000) ^ 2 / b ^ 2) ≤ 1 / (2 * b) := by
          rw [show 1 / (8 * b ^ 2) * 2 + 3 * ((1001 / 1000) ^ 2 / b ^ 2) = (1 / 4 + 3 * (1001 / 1000) ^ 2) / b ^ 2 by
            field_simp; ring, div_le_div_iff₀ (by positivity) (by positivity)]
          nlinarith only [hb]
        linarith
      nlinarith only [hbt, h3, h4]
    linarith
  · rw [eμ]
    calc 2 * t / (1 - t ^ 2) ≤ 201 / 100 / b := hμb
      _ = 201 / 100 * (1 / b) := by ring
      _ ≤ 201 / 100 * (2 / v ^ 3) := by gcongr
      _ ≤ 5 / v ^ 3 := by rw [mul_div_assoc', div_le_div_iff₀ hv3p hv3p]; nlinarith only [hv3p]
  · rw [eμ]
    have h1 : 2 * t / (1 + t ^ 2) * (2 * t / (1 - t ^ 2)) ≤ 2 * t * (201 / 100 / b) :=
      mul_le_mul (by rw [div_le_iff₀ hp]; nlinarith only [htpos]) hμb (by positivity) (by positivity)
    have h2 : 2 * t * (201 / 100 / b) ≤ 1 / (2 * b) := by
      rw [mul_div_assoc', div_le_div_iff₀ hb0 (by positivity)]; nlinarith only [ht1, hb, htpos]
    linarith
  · rw [ew]; linarith
  · rw [ew]; linarith

end band

/-! ## Rows and the waste identity -/

/-- Rows `R = ⌊(L − h − 2Y) γ⌋ + 1`: the top lift `L − H(R) − Y` is in `[Y, Y + 2]`, and `(R − 1) μ ≤ Lσ`. -/
lemma rows_exist {γ σ L Y h : ℝ} (hγ : 0 < γ) (hiγ : 1 / γ ≤ 2) (hσ : 0 ≤ σ) (hh0 : 0 ≤ h) (hY : 0 ≤ Y)
    (hL : h + 2 * Y ≤ L) :
    ∃ R : ℕ, 1 ≤ R ∧ Y ≤ L - (((R : ℝ) - 1) / γ + h) - Y ∧ L - (((R : ℝ) - 1) / γ + h) - Y ≤ Y + 2 ∧
      ((R : ℝ) - 1) * (σ / γ) ≤ L * σ := by
  have hpos : 0 ≤ (L - h - 2 * Y) * γ := mul_nonneg (by linarith) hγ.le
  refine ⟨⌊(L - h - 2 * Y) * γ⌋₊ + 1, by omega, ?_, ?_, ?_⟩ <;> push_cast <;> simp only [add_sub_cancel_right]
  · have h1 := Nat.floor_le hpos
    have : (⌊(L - h - 2 * Y) * γ⌋₊ : ℝ) / γ ≤ L - h - 2 * Y := by rw [div_le_iff₀ hγ]; linarith
    linarith
  · have h1 := Nat.lt_floor_add_one ((L - h - 2 * Y) * γ)
    have : L - h - 2 * Y - 1 / γ ≤ (⌊(L - h - 2 * Y) * γ⌋₊ : ℝ) / γ := by
      rw [le_div_iff₀ hγ, sub_mul, one_div_mul_cancel hγ.ne']; linarith
    linarith
  · have h1 := Nat.floor_le hpos
    have e : (⌊(L - h - 2 * Y) * γ⌋₊ : ℝ) * (σ / γ) ≤ (L - h - 2 * Y) * γ * (σ / γ) :=
      mul_le_mul_of_nonneg_right h1 (div_nonneg hσ hγ.le)
    have e2 : (L - h - 2 * Y) * γ * (σ / γ) = (L - h - 2 * Y) * σ := by field_simp
    have e3 : (L - h - 2 * Y) * σ ≤ L * σ := mul_le_mul_of_nonneg_right (by linarith) hσ
    linarith

/-- **The waste identity** (Ryu's Proposition 5.2, one band): `bR + A₀(y₀) + A₀(y₁) = wL − (R − 1) tan θ` for
`y₁ = L − H(R) − y₀`. -/
lemma ident_band {b γ σ L Y : ℝ} (R : ℝ) (hγ : γ ≠ 0) (hu : γ ^ 2 + σ ^ 2 = 1) :
    R * b + (Y * (b * γ + σ) + γ * σ * (1 + b ^ 2) / 2) +
      ((L - ((R - 1) / γ + b * σ + γ) - Y) * (b * γ + σ) + γ * σ * (1 + b ^ 2) / 2) =
      (b * γ + σ) * L - (R - 1) * (σ / γ) := by
  have key : R * b + (Y * (b * γ + σ) + γ * σ * (1 + b ^ 2) / 2) +
      ((L - ((R - 1) / γ + b * σ + γ) - Y) * (b * γ + σ) + γ * σ * (1 + b ^ 2) / 2) -
      ((b * γ + σ) * L - (R - 1) * (σ / γ)) = -b * (γ ^ 2 + σ ^ 2 - 1) := by
    field_simp; ring
  rw [hu, sub_self, mul_zero, sub_eq_zero] at key
  exact key

/-! ## The packing -/

set_option maxHeartbeats 2000000 in
/-- **The L-shape with staircase ends.** -/
theorem packs_stair (k b : ℕ) (v : ℝ) (hv : 10 ≤ v) (hb1 : (b : ℝ) ≤ v ^ 3) (hb2 : v ^ 3 ≤ 2 * b)
    (hk : (b : ℝ) + 2 * v ^ 2 + 10 ≤ k) (c : ℕ) (hc : 4 * k / b + 4 * Eub v ≤ c) :
    ∃ s < (k : ℝ), Packs (k * k - c) s := by
  have hv0 : 0 < v := by linarith
  have hv2 : 100 ≤ v ^ 2 := by nlinarith only [hv]
  have hv3 : 1000 ≤ v ^ 3 := by nlinarith only [hv, hv2]
  have hb : (500 : ℝ) ≤ b := by linarith
  have hb0 : (0 : ℝ) < b := by linarith
  have hK : (0 : ℝ) < k := by linarith
  have hbk : (b : ℝ) ≤ k := by linarith
  -- the band angle `t = 1/b + δ`, `δ = 1/(8kb³)`
  set δ : ℝ := 1 / (8 * k * b ^ 3) with hδd
  have hδ0 : 0 < δ := by positivity
  set t : ℝ := 1 / b + δ with ht
  have ht0 : 1 / (b : ℝ) < t := by rw [ht]; linarith
  have hbδ : (b : ℝ) * δ ≤ 1 / (8 * b ^ 2) := by
    rw [hδd, mul_one_div, div_le_div_iff₀ (by positivity) (by positivity)]
    have hk1 : (1 : ℝ) ≤ k := by linarith
    have hb3 : (0 : ℝ) < b ^ 3 := by positivity
    nlinarith only [hk1, hb3]
  have hbt : (b : ℝ) * t ≤ 1 + 1 / (8 * b ^ 2) := by
    rw [ht, mul_add, mul_one_div_cancel hb0.ne']; linarith
  have BF := bandFacts hb ht0 hbt hv hb1 hb2
  have hp : 0 < 1 + t ^ 2 := by positivity
  have hu : ((1 - t ^ 2) / (1 + t ^ 2)) ^ 2 + (2 * t / (1 + t ^ 2)) ^ 2 = 1 := by field_simp; ring
  obtain ⟨θ, hc', hs'⟩ := exists_angle hu
  -- plain names for `cos θ`, `sin θ` and the band width
  obtain ⟨γ, hγc⟩ : ∃ γ, Real.cos θ = γ := ⟨_, rfl⟩
  obtain ⟨σ, hσs⟩ : ∃ σ, Real.sin θ = σ := ⟨_, rfl⟩
  rw [hγc] at hc'
  rw [hσs] at hs'
  obtain ⟨w, hwd⟩ : ∃ w, (b : ℝ) * γ + σ = w := ⟨_, rfl⟩
  have hγ : 0 < γ := by rw [hc']; exact BF.hγ
  have hσ : 0 < σ := by rw [hs']; exact BF.hσ
  have huθ : γ ^ 2 + σ ^ 2 = 1 := by rw [hc', hs']; exact hu
  have ew : w = ww b t := by rw [← hwd, hc', hs', ww]; field_simp
  have hwb : w < b := by rw [ew]; exact w_lt hb0 ht0
  have hwge : (b : ℝ) - 2 * b * δ ≤ w := by
    have := w_ge hb0 ht0.le; rw [← ew, ht, add_sub_cancel_left] at this; linarith
  have hσt : σ ≤ 2 * t := by rw [hs']; exact BF.hσt
  have hiγ : 1 / γ ≤ 2 := by rw [hc']; exact BF.hiγ
  have hh : (b : ℝ) * σ + γ ≤ 4 := by rw [hc', hs']; exact BF.hh
  have hw1 : v ^ 3 / 3 ≤ w := by rw [← hwd, hc', hs']; exact BF.hw1
  have hw0 : 0 ≤ w := by linarith
  -- the end packings
  have E := fun (y : ℝ) (h1 : 3 / 4 * v ^ 2 ≤ y) (h2 : y ≤ 3 / 4 * v ^ 2 + 2) =>
    end_bound (b := b) (γ := γ) (σ := σ) (w := w) hv hγ hσ huθ hwd.symm
      (by rw [← hwd, hc', hs']; exact BF.hΔ1) (by rw [← hwd, hc', hs']; exact BF.hΔ2) (by rw [hc', hs']; exact BF.hμ)
      (by rw [hc', hs']; exact BF.hσμ) hw1 (by rw [← hwd, hc', hs']; exact BF.hw2) h1 h2
  -- the grid and the two bands
  obtain ⟨mg, hmg⟩ : ∃ mg : ℕ, mg = k - b := ⟨_, rfl⟩
  have hmgR : (mg : ℝ) = k - b := by rw [hmg, Nat.cast_sub (by exact_mod_cast hbk)]
  obtain ⟨Y, hY⟩ : ∃ Y : ℝ, Y = 3 / 4 * v ^ 2 := ⟨_, rfl⟩
  have hY0 : 0 ≤ Y := by rw [hY]; positivity
  have hLh : (b : ℝ) * σ + γ + 2 * Y ≤ mg := by rw [hmgR, hY]; linarith
  have hh0 : 0 ≤ (b : ℝ) * σ + γ := by positivity
  obtain ⟨R₁, hR₁, hy₁a, hy₁b, hR₁μ⟩ := rows_exist (L := (mg : ℝ) + w) (h := (b : ℝ) * σ + γ) hγ hiγ hσ.le hh0 hY0
    (by linarith)
  obtain ⟨R₂, hR₂, hy₂a, hy₂b, hR₂μ⟩ := rows_exist (L := (mg : ℝ)) (h := (b : ℝ) * σ + γ) hγ hiγ hσ.le hh0 hY0 hLh
  obtain ⟨n0, P0, hn0⟩ := E Y hY.symm.le (by rw [hY]; linarith)
  obtain ⟨n1, P1, hn1⟩ := E ((mg : ℝ) + w - (((R₁ : ℝ) - 1) / γ + ((b : ℝ) * σ + γ)) - Y)
    (by rw [← hY]; exact hy₁a) (by rw [← hY]; exact hy₁b)
  obtain ⟨n2, P2, hn2⟩ := E ((mg : ℝ) - (((R₂ : ℝ) - 1) / γ + ((b : ℝ) * σ + γ)) - Y)
    (by rw [← hY]; exact hy₂a) (by rw [← hY]; exact hy₂b)
  have hγθ : 0 < Real.cos θ := by rw [hγc]; exact hγ
  have hσθ : 0 ≤ Real.sin θ := by rw [hσs]; exact hσ.le
  have B1 := band_ends (θ := θ) (L := (mg : ℝ) + w) hγθ hσθ R₁ b hR₁ hY0 (le_trans hY0 hy₁a)
    (by rw [hγc, hσs]; exact le_of_eq (by ring)) (by rw [hγc, hσs, hwd]; exact P0) (by rw [hγc, hσs, hwd]; exact P1)
  have B2 := band_ends (θ := θ) (L := (mg : ℝ)) hγθ hσθ R₂ b hR₂ hY0 (le_trans hY0 hy₂a)
    (by rw [hγc, hσs]; exact le_of_eq (by ring)) (by rw [hγc, hσs, hwd]; exact P0) (by rw [hγc, hσs, hwd]; exact P2)
  rw [hγc, hσs, hwd] at B1 B2
  have P := lshape hw0 B1 B2
  refine ⟨(mg : ℝ) + w, by rw [hmgR]; linarith, packs_of_le P ?_⟩
  -- the count
  have i1 := ident_band (b := (b : ℝ)) (L := (mg : ℝ) + w) (Y := Y) (R₁ : ℝ) hγ.ne' huθ
  have i2 := ident_band (b := (b : ℝ)) (L := (mg : ℝ)) (Y := Y) (R₂ : ℝ) hγ.ne' huθ
  rw [hwd] at i1 i2
  have hS : (k : ℝ) * k - ((mg : ℝ) * mg + w * ((mg : ℝ) + w) + w * mg) ≤ 2 * k * (2 * b * δ) := by
    rw [hmgR]
    have e : (k : ℝ) * k - ((k - b) * (k - b) + w * ((k - b) + w) + w * (k - b)) = (b - w) * (2 * k - b + w) := by ring
    rw [e]
    have h1 : (b : ℝ) - w ≤ 2 * b * δ := by linarith
    have h2 : 0 ≤ (b : ℝ) - w := by linarith
    have h3 : 2 * (k : ℝ) - b + w ≤ 2 * k := by linarith
    have h4 : 0 ≤ 2 * (k : ℝ) - b + w := by linarith
    calc ((b : ℝ) - w) * (2 * k - b + w) ≤ (2 * b * δ) * (2 * k) := mul_le_mul h1 h3 h4 (by positivity)
      _ = 2 * k * (2 * b * δ) := by ring
  have hLσ : ((mg : ℝ) + w) * σ + mg * σ ≤ 2 * k * (2 * t) := by
    rw [hmgR]
    have h1 : ((k : ℝ) - b + w) + (k - b) ≤ 2 * k := by linarith
    have h2 : (((k : ℝ) - b + w) + (k - b)) * σ ≤ 2 * k * σ := mul_le_mul_of_nonneg_right h1 hσ.le
    have h3 : 2 * (k : ℝ) * σ ≤ 2 * k * (2 * t) := mul_le_mul_of_nonneg_left hσt (by positivity)
    linarith
  have h2kt : 2 * (k : ℝ) * (2 * t) = 4 * k / b + 4 * k * δ := by rw [ht]; field_simp; norm_num
  have hsmall : 2 * (k : ℝ) * (2 * b * δ) + 4 * k * δ < 1 := by
    have e : 2 * (k : ℝ) * (2 * b * δ) + 4 * k * δ = ((b : ℝ) + 1) / (2 * b ^ 3) := by
      rw [hδd]; field_simp; ring
    rw [e, div_lt_one (by positivity)]
    have hb2' : (250000 : ℝ) ≤ b ^ 2 := by nlinarith only [hb]
    nlinarith only [hb, hb2']
  have hE4 : (k : ℝ) * k - c < (mg : ℝ) * mg + (n0 + R₁ * b + n1) + (n0 + R₂ * b + n2) + 1 := by
    unfold Eub at hc
    linarith only [i1, i2, hn0, hn1, hn2, hR₁μ, hR₂μ, hS, hLσ, h2kt, hsmall, hc]
  have hNn : (k * k - c : ℕ) ≤ mg * mg + (n0 + R₁ * b + n1) + (n0 + R₂ * b + n2) := by
    have : ((k * k : ℕ) : ℝ) < (mg * mg + (n0 + R₁ * b + n1) + (n0 + R₂ * b + n2) + c + 1 : ℕ) := by
      push_cast; linarith
    have := Nat.cast_lt.1 this
    omega
  exact hNn

/-! ## Closed forms -/

lemma rpow_fifth (k : ℕ) : ((k : ℝ) ^ ((1 : ℝ) / 5)) ^ 5 = k := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul (Nat.cast_nonneg k)]; norm_num

lemma rpow_two_fifths (k : ℕ) : (k : ℝ) ^ ((2 : ℝ) / 5) = ((k : ℝ) ^ ((1 : ℝ) / 5)) ^ 2 := by
  rw [← Real.rpow_natCast, ← Real.rpow_mul (Nat.cast_nonneg k)]; norm_num

/-- From `v = 0.81 z` (`z⁵ = k`, `b = ⌊v³⌋`) to a packing, given the polynomial inequality. -/
lemma packs_of_z (k : ℕ) {z X : ℝ} (hz5 : z ^ 5 = k) (hz : 124 / 10 ≤ z)
    (hX : 4 * z ^ 5 / ((81 / 100 * z) ^ 3 - 1) + 4 * Eub (81 / 100 * z) ≤ X) :
    ∃ s < (k : ℝ), Packs (k * k - ⌈X⌉₊) s := by
  have hz0 : 0 ≤ z := by linarith
  set v := 81 / 100 * z with hvd
  have hv : 10 ≤ v := by rw [hvd]; linarith
  have hvz : v ≤ z := by rw [hvd]; linarith
  have hv0 : 0 ≤ v := by linarith
  have hv2 : 100 ≤ v ^ 2 := by nlinarith only [hv]
  have hv3 : 1000 ≤ v ^ 3 := by nlinarith only [hv, hv2]
  set b := ⌊v ^ 3⌋₊ with hbd
  have hb1 : (b : ℝ) ≤ v ^ 3 := Nat.floor_le (by positivity)
  have hb0 : v ^ 3 - 1 < b := by have := Nat.lt_floor_add_one (v ^ 3); rw [← hbd] at this; linarith
  have hb2 : v ^ 3 ≤ 2 * b := by linarith
  have hz2 : 153 ≤ z ^ 2 := by nlinarith only [hz]
  have hz53 : 153 * z ^ 3 ≤ z ^ 5 := by nlinarith only [hz2, pow_nonneg hz0 3]
  have hv3z : v ^ 3 ≤ z ^ 3 := pow_le_pow_left₀ hv0 hvz 3
  have hv2z : v ^ 2 ≤ z ^ 2 := pow_le_pow_left₀ hv0 hvz 2
  have hz32 : 12 * z ^ 2 ≤ z ^ 3 := by nlinarith only [hz, sq_nonneg z]
  have hk : (b : ℝ) + 2 * v ^ 2 + 10 ≤ k := by
    rw [← hz5]; nlinarith only [hb1, hv3z, hv2z, hz53, hz32, hz2]
  refine packs_stair k b v hv hb1 hb2 hk ⌈X⌉₊ ?_
  have h1 : 4 * (k : ℝ) / b ≤ 4 * z ^ 5 / (v ^ 3 - 1) := by
    rw [← hz5]; apply div_le_div_of_nonneg_left (by positivity) (by linarith) hb0.le
  exact le_trans (by linarith) (Nat.le_ceil X)

/-- **`c*(k) < 18.9 k^{2/5} + 35.7 k^{1/5} + 122`** for every `k ≥ 3 · 10⁵`. -/
theorem packs_stair_k25 (k : ℕ) (hk : 300000 ≤ k) :
    ∃ s < (k : ℝ), Packs (k * k - ⌈189 / 10 * (k : ℝ) ^ ((2 : ℝ) / 5) + 357 / 10 * (k : ℝ) ^ ((1 : ℝ) / 5) + 122⌉₊)
      s := by
  set z := (k : ℝ) ^ ((1 : ℝ) / 5) with hzd
  have hz5 : z ^ 5 = k := rpow_fifth k
  have hz0 : 0 ≤ z := by positivity
  have hz : 124 / 10 ≤ z := by
    by_contra hc
    push_neg at hc
    have : z ^ 5 < (124 / 10) ^ 5 := pow_lt_pow_left₀ hc hz0 (by norm_num)
    have : (300000 : ℝ) ≤ k := by exact_mod_cast hk
    linarith
  rw [rpow_two_fifths, ← hzd]
  apply packs_of_z k hz5 hz
  unfold Eub
  have hden : 0 < (81 / 100 * z) ^ 3 - 1 := by nlinarith only [hz, pow_le_pow_left₀ (by norm_num) hz 3]
  rw [div_add' _ _ _ hden.ne', div_le_iff₀ hden]
  have hz2 : 153 ≤ z ^ 2 := by nlinarith only [hz]
  nlinarith only [hz, hz2, pow_nonneg hz0 3, pow_nonneg hz0 4, mul_nonneg (pow_nonneg hz0 3) (by linarith : (0 : ℝ) ≤ z ^ 2 - 153)]

/-- **Ryu's Theorem 1.2**: `c*(k) < 20.668 k^{2/5} + 0.623`, here for every `k ≥ 8 · 10⁶` (Ryu: `k ≥ 1.6 · 10⁷`). -/
theorem packs_ryu_thm12 (k : ℕ) (hk : 8000000 ≤ k) :
    ∃ s < (k : ℝ), Packs (k * k - ⌈20668 / 1000 * (k : ℝ) ^ ((2 : ℝ) / 5) + 623 / 1000⌉₊) s := by
  set z := (k : ℝ) ^ ((1 : ℝ) / 5) with hzd
  have hz5 : z ^ 5 = k := rpow_fifth k
  have hz0 : 0 ≤ z := by positivity
  have hz : 24 ≤ z := by
    by_contra hc
    push_neg at hc
    have : z ^ 5 < 24 ^ 5 := pow_lt_pow_left₀ hc hz0 (by norm_num)
    have : (8000000 : ℝ) ≤ k := by exact_mod_cast hk
    linarith
  rw [rpow_two_fifths, ← hzd]
  apply packs_of_z k hz5 (by linarith)
  unfold Eub
  have hden : 0 < (81 / 100 * z) ^ 3 - 1 := by nlinarith only [hz, pow_le_pow_left₀ (by norm_num) hz 3]
  rw [div_add' _ _ _ hden.ne', div_le_iff₀ hden]
  have hz2 : 576 ≤ z ^ 2 := by nlinarith only [hz]
  nlinarith only [hz, hz2, pow_nonneg hz0 3, pow_nonneg hz0 4,
    mul_nonneg (pow_nonneg hz0 3) (by linarith : (0 : ℝ) ≤ z ^ 2 - 576),
    mul_nonneg (pow_nonneg hz0 3) (by linarith : (0 : ℝ) ≤ z - 24)]

end UnitSquarePacking.Ryu
