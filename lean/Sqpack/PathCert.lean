import Sqpack.LocalMin

/-!
# Path certificates: two packings joined below a level

`search/packer/analysis/pathcert.py` decides that two packings are the same local minimum by a chain of
configurations `A₀, A₁, …` (centres, angles, box side) such that, for every step `A → B` and every `t ∈ [0, 1]`,
the linear interpolation `(1 - t) A + t B`, scaled about the origin to side `L`, is a packing.  It checks each step
with one inequality per constraint (a wall, or a near pair with a chosen separating side line):

  `(1 - t) g_A + t g_B - t (1 - t) q / 2 + k (1 - s(t) / L) ≥ 0`   for all `t ∈ [0, 1]`,

where `g` is the constraint's margin, `s(t)` the interpolated side, `k = 1` for pairs and `1/2` for walls, and `q`
an explicit bound on the margin's second derivative.  This file proves that check sound.

* §1 `chord_lower`: `f'' ≤ q` on `[0, 1]` gives `f t ≥ (1 - t) f 0 + t f 1 - t (1 - t) q / 2`.
* §2 the margins along a step: the half-width `hw φ = (|cos φ| + |sin φ|) / 2` of a unit square seen at relative
  angle `φ` (`hw_chord`, as a maximum of four cosines, so no kink argument), the pair margin with
  `q = W Δθ_o² + 2 D |Δθ_o| + (√2/2) Δφ²` (`pairGap_chord`), and the wall margin with `q = (√2/2) Δθ²`
  (`wallGap_chord`).
* §3 one configuration: a nonnegative pair margin separates (`disjoint_of_pairGap`, through `ExactPack`'s
  `SepSide`), nonnegative wall margins put the square in the box, centres `√2` apart never overlap.
* §4 scaling: with `ρ = L / s ≥ 1`, margin `g` becomes at least `ρ g + (ρ - 1)` (pairs) or `ρ m + (ρ - 1)/2` (walls),
  so the check's inequality at `t` makes the scaled configuration a packing (`step_isPacking`).
* §5 the path: `step_joined` (each step is a path inside the packings of side `≤ L`) and `chain_joined`.
* §6 `qmin_nonneg`, `cond_eq`: the closed-form minimum the script evaluates (`_cond`, `_qmin`) covers `[0, 1]`.
* §7 `joined_of_chain`: a certified chain whose ends are packings joins them by a path of packings, all in boxes of
  side `≤ L`: the two packings are in one component of the sublevel set, i.e. the same basin at level `L`.

What stays outside: the script evaluates the margins at the nodes and the minima in 50-digit floating point with a
`1e-30` safety margin (not interval arithmetic), and its float pre-filter decides which pairs are "far".
-/

namespace UnitSquarePacking

open Set Real

/-! ## 1. The chord bound -/

/-- `f'' ≤ q` on `[0, 1]` bounds `f` below by its chord minus `t (1 - t) q / 2`. -/
lemma chord_lower {f f' f'' : ℝ → ℝ} {q : ℝ} (hf : ∀ x, HasDerivAt f (f' x) x)
    (hf' : ∀ x, HasDerivAt f' (f'' x) x) (hq : ∀ x ∈ Icc (0 : ℝ) 1, f'' x ≤ q) {t : ℝ}
    (ht : t ∈ Icc (0 : ℝ) 1) :
    (1 - t) * f 0 + t * f 1 - t * (1 - t) * q / 2 ≤ f t := by
  set u : ℝ → ℝ := fun x => f x - q / 2 * (x * x) with hu
  have hu' : ∀ x, HasDerivAt u (f' x - q * x) x := fun x =>
    ((hf x).sub (((hasDerivAt_id x).mul (hasDerivAt_id x)).const_mul (q / 2))).congr_deriv
      (by simp only [id]; ring)
  have hu'' : ∀ x, HasDerivAt (fun x => f' x - q * x) (f'' x - q) x := fun x =>
    ((hf' x).sub ((hasDerivAt_id x).const_mul q)).congr_deriv (by simp)
  have hc : ConcaveOn ℝ (Icc (0 : ℝ) 1) u :=
    concaveOn_of_hasDerivWithinAt2_nonpos (convex_Icc 0 1)
      (fun x _ => (hu' x).continuousAt.continuousWithinAt)
      (fun x _ => (hu' x).hasDerivWithinAt) (fun x _ => (hu'' x).hasDerivWithinAt)
      (fun x hx => by have := hq x (interior_subset hx); linarith)
  have h := hc.2 (left_mem_Icc.2 zero_le_one) (right_mem_Icc.2 zero_le_one) (sub_nonneg.2 ht.2) ht.1
    (sub_add_cancel 1 t)
  simp only [hu, smul_eq_mul] at h
  have e : (1 - t) * 0 + t * 1 = t := by ring
  rw [e] at h
  nlinarith

/-! ## 2. Margins along a step -/

/-- Half the width of a unit square, measured across a direction at angle `φ` to its sides. -/
noncomputable def hw (φ : ℝ) : ℝ := (|cos φ| + |sin φ|) / 2

lemma hw_nonneg (φ : ℝ) : 0 ≤ hw φ := by unfold hw; positivity

/-- `a x + b y ≤ W` when `a² + b² ≤ 1` and `x² + y² ≤ W²`, `W ≥ 0` (Cauchy–Schwarz). -/
lemma dot_le {a b x y W : ℝ} (hab : a ^ 2 + b ^ 2 ≤ 1) (hxy : x ^ 2 + y ^ 2 ≤ W ^ 2) (hW : 0 ≤ W) :
    a * x + b * y ≤ W := by
  have h1 : (a * x + b * y) ^ 2 ≤ W ^ 2 := by
    nlinarith [sq_nonneg (a * y - b * x), sq_nonneg a, sq_nonneg b, sq_nonneg x, sq_nonneg y]
  exact (abs_le_of_sq_le_sq' h1 hW).2

/-- A sign `σ = ±1` with `σ c = |c|`. -/
lemma exists_sign (c : ℝ) : ∃ σ : ℝ, σ * c = |c| ∧ σ ^ 2 = 1 ∧ ∀ d, σ * d ≤ |d| := by
  rcases le_total 0 c with h | h
  · exact ⟨1, by rw [abs_of_nonneg h]; ring, by norm_num, fun d => by linarith [le_abs_self d]⟩
  · exact ⟨-1, by rw [abs_of_nonpos h]; ring, by norm_num, fun d => by linarith [neg_abs_le d]⟩

lemma hasDerivAt_affine (a d x : ℝ) : HasDerivAt (fun x => a + x * d) d x :=
  ((hasDerivAt_id x).mul_const d).const_add a |>.congr_deriv (by simp)

/-- The half-width along a linear change of angle stays below its chord plus `t (1 - t) (√2/2) Δφ² / 2`. -/
lemma hw_chord (φ₀ dφ : ℝ) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    hw (φ₀ + t * dφ) ≤ (1 - t) * hw φ₀ + t * hw (φ₀ + 1 * dφ) + t * (1 - t) * (√2 / 2 * dφ ^ 2) / 2 := by
  obtain ⟨σ, hσ, hσ2, hσd⟩ := exists_sign (cos (φ₀ + t * dφ))
  obtain ⟨τ, hτ, hτ2, hτd⟩ := exists_sign (sin (φ₀ + t * dφ))
  -- the piece p = (σ cos + τ sin) / 2 equals hw at t and lies below hw everywhere
  set p : ℝ → ℝ := fun x => -((σ * cos (φ₀ + x * dφ) + τ * sin (φ₀ + x * dφ)) / 2) with hp
  set p' : ℝ → ℝ := fun x => -((σ * (-sin (φ₀ + x * dφ) * dφ) + τ * (cos (φ₀ + x * dφ) * dφ)) / 2)
  set p'' : ℝ → ℝ := fun x => -((σ * (-(cos (φ₀ + x * dφ) * dφ) * dφ) + τ * (-sin (φ₀ + x * dφ) * dφ * dφ)) / 2)
  have d1 : ∀ x, HasDerivAt p (p' x) x := fun x =>
    (((hasDerivAt_affine φ₀ dφ x).cos.const_mul σ).add ((hasDerivAt_affine φ₀ dφ x).sin.const_mul τ)).div_const 2
      |>.neg
  have d2 : ∀ x, HasDerivAt p' (p'' x) x := fun x =>
    ((((hasDerivAt_affine φ₀ dφ x).sin.neg.mul_const dφ).const_mul σ).add
      (((hasDerivAt_affine φ₀ dφ x).cos.mul_const dφ).const_mul τ)).div_const 2 |>.neg |>.congr_deriv
      (by simp only [p''])
  have hq : ∀ x ∈ Icc (0 : ℝ) 1, p'' x ≤ √2 / 2 * dφ ^ 2 := by
    intro x _
    have hs : (σ * cos (φ₀ + x * dφ) + τ * sin (φ₀ + x * dφ)) ≤ √2 := by
      have hc := cos_sq_add_sin_sq (φ₀ + x * dφ)
      have h2 : (σ * cos (φ₀ + x * dφ) + τ * sin (φ₀ + x * dφ)) ^ 2 ≤ √2 ^ 2 := by
        rw [Real.sq_sqrt (by norm_num : (0 : ℝ) ≤ 2)]
        nlinarith [sq_nonneg (σ * sin (φ₀ + x * dφ) - τ * cos (φ₀ + x * dφ))]
      exact abs_le_of_sq_le_sq' h2 (Real.sqrt_nonneg 2) |>.2
    have : p'' x = (σ * cos (φ₀ + x * dφ) + τ * sin (φ₀ + x * dφ)) / 2 * dφ ^ 2 := by simp only [p'']; ring
    rw [this]
    nlinarith [sq_nonneg dφ]
  have h := chord_lower d1 d2 hq ht
  simp only [hp, zero_mul, add_zero] at h
  have h0 : σ * cos φ₀ + τ * sin φ₀ ≤ 2 * hw φ₀ := by unfold hw; linarith [hσd (cos φ₀), hτd (sin φ₀)]
  have h1 : σ * cos (φ₀ + 1 * dφ) + τ * sin (φ₀ + 1 * dφ) ≤ 2 * hw (φ₀ + 1 * dφ) := by
    unfold hw; linarith [hσd (cos (φ₀ + 1 * dφ)), hτd (sin (φ₀ + 1 * dφ))]
  have ht' : hw (φ₀ + t * dφ) = (σ * cos (φ₀ + t * dφ) + τ * sin (φ₀ + t * dφ)) / 2 := by
    unfold hw; rw [hσ, hτ]
  rw [ht']
  nlinarith [ht.1, ht.2]

/-- The margin of a pair along a step: the owner's side normal at angle `α₀ + t dα`, the centre offset (other minus
owner) `(wx + t dwx, wy + t dwy)`, the other square's angle relative to the owner `φ₀ + t dφ`.  `≥ 0` means the other
square lies beyond that side line (§3). -/
noncomputable def pairGap (α₀ dα wx dwx wy dwy φ₀ dφ t : ℝ) : ℝ :=
  cos (α₀ + t * dα) * (wx + t * dwx) + sin (α₀ + t * dα) * (wy + t * dwy) - 1 / 2 - hw (φ₀ + t * dφ)

/-- The smooth part of the pair margin, `n(t) · w(t)`, has second derivative `≤ W dα² + 2 D |dα|`. -/
lemma dotPart_chord {α₀ dα wx dwx wy dwy W D : ℝ}
    (hW : ∀ x ∈ Icc (0 : ℝ) 1, (wx + x * dwx) ^ 2 + (wy + x * dwy) ^ 2 ≤ W ^ 2) (hW0 : 0 ≤ W)
    (hD : dwx ^ 2 + dwy ^ 2 ≤ D ^ 2) (hD0 : 0 ≤ D) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    let F := fun x => cos (α₀ + x * dα) * (wx + x * dwx) + sin (α₀ + x * dα) * (wy + x * dwy)
    (1 - t) * F 0 + t * F 1 - t * (1 - t) * (W * dα ^ 2 + 2 * D * |dα|) / 2 ≤ F t := by
  intro F
  set F' : ℝ → ℝ := fun x => -sin (α₀ + x * dα) * dα * (wx + x * dwx) + cos (α₀ + x * dα) * dwx +
    (cos (α₀ + x * dα) * dα * (wy + x * dwy) + sin (α₀ + x * dα) * dwy)
  set F'' : ℝ → ℝ := fun x => -(dα ^ 2) * (cos (α₀ + x * dα) * (wx + x * dwx) + sin (α₀ + x * dα) * (wy + x * dwy)) +
    2 * dα * (-sin (α₀ + x * dα) * dwx + cos (α₀ + x * dα) * dwy)
  have d1 : ∀ x, HasDerivAt F (F' x) x := fun x =>
    (((hasDerivAt_affine α₀ dα x).cos.mul (hasDerivAt_affine wx dwx x)).add
      ((hasDerivAt_affine α₀ dα x).sin.mul (hasDerivAt_affine wy dwy x)))
  have d2 : ∀ x, HasDerivAt F' (F'' x) x := fun x =>
    (((((hasDerivAt_affine α₀ dα x).sin.neg.mul_const dα).mul (hasDerivAt_affine wx dwx x)).add
      ((hasDerivAt_affine α₀ dα x).cos.mul_const dwx)).add
      ((((hasDerivAt_affine α₀ dα x).cos.mul_const dα).mul (hasDerivAt_affine wy dwy x)).add
      ((hasDerivAt_affine α₀ dα x).sin.mul_const dwy))).congr_deriv (by simp only [F'', Pi.neg_apply]; ring_nf)
  have hq : ∀ x ∈ Icc (0 : ℝ) 1, F'' x ≤ W * dα ^ 2 + 2 * D * |dα| := by
    intro x hx
    have hc := cos_sq_add_sin_sq (α₀ + x * dα)
    have e1 : -cos (α₀ + x * dα) * (wx + x * dwx) + -sin (α₀ + x * dα) * (wy + x * dwy) ≤ W :=
      dot_le (by nlinarith) (hW x hx) hW0
    have e2 : -sin (α₀ + x * dα) * dwx + cos (α₀ + x * dα) * dwy ≤ D := dot_le (by nlinarith) hD hD0
    have e3 : -(-sin (α₀ + x * dα) * dwx + cos (α₀ + x * dα) * dwy) ≤ D := by
      have := dot_le (a := sin (α₀ + x * dα)) (b := -cos (α₀ + x * dα)) (by nlinarith) hD hD0; linarith
    have e4 : dα * (-sin (α₀ + x * dα) * dwx + cos (α₀ + x * dα) * dwy) ≤ |dα| * D := by
      rcases le_total 0 dα with h | h
      · rw [abs_of_nonneg h]; exact mul_le_mul_of_nonneg_left e2 h
      · rw [abs_of_nonpos h]; nlinarith
    simp only [F'']
    nlinarith [sq_nonneg dα]
  exact chord_lower d1 d2 hq ht

/-- **Pair margin along a step.**  With `W ≥ |w(t)|` on `[0, 1]` and `D ≥ |Δw|`, the margin stays above its chord
minus `t (1 - t) q / 2`, `q = W dα² + 2 D |dα| + (√2/2) dφ²` (the script's `q` for a pair). -/
theorem pairGap_chord {α₀ dα wx dwx wy dwy φ₀ dφ W D : ℝ}
    (hW : ∀ x ∈ Icc (0 : ℝ) 1, (wx + x * dwx) ^ 2 + (wy + x * dwy) ^ 2 ≤ W ^ 2) (hW0 : 0 ≤ W)
    (hD : dwx ^ 2 + dwy ^ 2 ≤ D ^ 2) (hD0 : 0 ≤ D) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    (1 - t) * pairGap α₀ dα wx dwx wy dwy φ₀ dφ 0 + t * pairGap α₀ dα wx dwx wy dwy φ₀ dφ 1 -
        t * (1 - t) * (W * dα ^ 2 + 2 * D * |dα| + √2 / 2 * dφ ^ 2) / 2 ≤
      pairGap α₀ dα wx dwx wy dwy φ₀ dφ t := by
  have h1 := dotPart_chord (α₀ := α₀) (dα := dα) hW hW0 hD hD0 ht
  have h2 := hw_chord φ₀ dφ ht
  simp only [zero_mul, add_zero] at h1 h2 ⊢
  unfold pairGap
  simp only [zero_mul, add_zero, one_mul] at h1 h2 ⊢
  nlinarith [ht.1, ht.2]

/-- A wall margin along a step: `m₀ + t dm` (the centre's distance to the wall, linear) minus the half-width at
angle `θ₀ + t dθ`. -/
noncomputable def wallGap (m₀ dm θ₀ dθ t : ℝ) : ℝ := m₀ + t * dm - hw (θ₀ + t * dθ)

/-- **Wall margin along a step**: `q = (√2/2) dθ²`. -/
theorem wallGap_chord (m₀ dm θ₀ dθ : ℝ) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    (1 - t) * wallGap m₀ dm θ₀ dθ 0 + t * wallGap m₀ dm θ₀ dθ 1 - t * (1 - t) * (√2 / 2 * dθ ^ 2) / 2 ≤
      wallGap m₀ dm θ₀ dθ t := by
  have h2 := hw_chord θ₀ dθ ht
  unfold wallGap
  simp only [zero_mul, add_zero, one_mul] at h2 ⊢
  nlinarith [ht.1, ht.2]

/-- `|w(t)|² ≤ W²` on `[0, 1]` from the two ends (the squared norm is convex). -/
lemma norm_sq_le_of_ends {wx dwx wy dwy W : ℝ} (h0 : wx ^ 2 + wy ^ 2 ≤ W ^ 2)
    (h1 : (wx + dwx) ^ 2 + (wy + dwy) ^ 2 ≤ W ^ 2) :
    ∀ x ∈ Icc (0 : ℝ) 1, (wx + x * dwx) ^ 2 + (wy + x * dwy) ^ 2 ≤ W ^ 2 := by
  intro x ⟨hx0, hx1⟩
  have e : (wx + x * dwx) ^ 2 + (wy + x * dwy) ^ 2 =
      (1 - x) * (wx ^ 2 + wy ^ 2) + x * ((wx + dwx) ^ 2 + (wy + dwy) ^ 2) - x * (1 - x) * (dwx ^ 2 + dwy ^ 2) := by
    ring
  rw [e]
  nlinarith [mul_nonneg hx0 (sub_nonneg.2 hx1), sq_nonneg dwx, sq_nonneg dwy,
    mul_nonneg (mul_nonneg hx0 (sub_nonneg.2 hx1)) (add_nonneg (sq_nonneg dwx) (sq_nonneg dwy))]

/-! ## 3. One configuration -/

lemma half_mul_ge {a : ℝ} (ha : Half a) (u : ℝ) : -(|u| / 2) ≤ a * u := by
  rcases ha with rfl | rfl
  · linarith [neg_abs_le u]
  · linarith [le_abs_self u]

lemma one_half_le_hw (φ : ℝ) : 1 / 2 ≤ hw φ := by
  unfold hw
  have h := cos_sq_add_sin_sq φ
  have : 1 ≤ |cos φ| + |sin φ| := by
    nlinarith [abs_nonneg (cos φ), abs_nonneg (sin φ), sq_abs (cos φ), sq_abs (sin φ),
      mul_nonneg (abs_nonneg (cos φ)) (abs_nonneg (sin φ))]
  linarith

lemma hw_le_one (φ : ℝ) : hw φ ≤ 1 := by
  unfold hw; linarith [abs_cos_le_one φ, abs_sin_le_one φ]

/-- The static margin of square `u` against side `k` of its owner `o`. -/
noncomputable def gapAt (xo yo θo xu yu θu : ℝ) (k : Fin 4) : ℝ :=
  (nrm (cos θo) (sin θo) k).1 * (xu - xo) + (nrm (cos θo) (sin θo) k).2 * (yu - yo) - 1 / 2 - hw (θu - θo)

/-- A corner of square `u` seen along side `k` of `o` is at least `-hw` from `u`'s centre. -/
lemma corner_proj_ge (θo θu : ℝ) (k : Fin 4) {a b : ℝ} (ha : Half a) (hb : Half b) :
    -hw (θu - θo) ≤ (nrm (cos θo) (sin θo) k).1 * (cos θu * a - sin θu * b) +
      (nrm (cos θo) (sin θo) k).2 * (sin θu * a + cos θu * b) := by
  have hC : cos (θu - θo) = cos θu * cos θo + sin θu * sin θo := cos_sub θu θo
  have hS : sin (θu - θo) = sin θu * cos θo - cos θu * sin θo := sin_sub θu θo
  unfold hw
  have a1 := half_mul_ge ha (cos (θu - θo)); have a2 := half_mul_ge ha (sin (θu - θo))
  have a3 := half_mul_ge ha (-cos (θu - θo)); have a4 := half_mul_ge ha (-sin (θu - θo))
  have b1 := half_mul_ge hb (cos (θu - θo)); have b2 := half_mul_ge hb (sin (θu - θo))
  have b3 := half_mul_ge hb (-cos (θu - θo)); have b4 := half_mul_ge hb (-sin (θu - θo))
  simp only [abs_neg] at a3 a4 b3 b4
  fin_cases k <;> simp only [nrm] <;> rw [hC, hS] at * <;> nlinarith

/-- **A nonnegative pair margin separates.** -/
lemma disjoint_of_gapAt {xo yo θo xu yu θu : ℝ} {k : Fin 4} (h : 0 ≤ gapAt xo yo θo xu yu θu k) :
    Disjoint (interior (unitSq (xo, yo) θo)) (interior (unitSq (xu, yu) θu)) := by
  refine disjoint_of_sepSide rfl rfl rfl rfl (cos_sq_add_sin_sq θo) (cos_sq_add_sin_sq θu) (k := k) ?_
  intro a b ha hb
  have := corner_proj_ge θo θu k ha hb
  unfold gapAt at h
  nlinarith

/-- **Nonnegative wall margins put the square in the box.** -/
lemma subset_container_of_walls {S x y θ : ℝ} (h1 : 0 ≤ x - hw θ) (h2 : 0 ≤ S - x - hw θ)
    (h3 : 0 ≤ y - hw θ) (h4 : 0 ≤ S - y - hw θ) : unitSq (x, y) θ ⊆ container S := by
  refine unitSq_subset_container rfl rfl ?_
  intro a b ha hb
  unfold hw at h1 h2 h3 h4
  have := half_mul_ge ha (cos θ); have := half_mul_ge ha (-cos θ)
  have := half_mul_ge ha (sin θ); have := half_mul_ge ha (-sin θ)
  have := half_mul_ge hb (cos θ); have := half_mul_ge hb (-cos θ)
  have := half_mul_ge hb (sin θ); have := half_mul_ge hb (-sin θ)
  simp only [abs_neg] at *
  refine ⟨?_, ?_, ?_, ?_⟩ <;> nlinarith

/-- **Centres `√2` apart never overlap** (each square lies in the disc of radius `√2/2` about its centre). -/
lemma disjoint_of_far {xi yi θi xj yj θj : ℝ} (h : 2 ≤ (xj - xi) ^ 2 + (yj - yi) ^ 2) :
    Disjoint (interior (unitSq (xi, yi) θi)) (interior (unitSq (xj, yj) θj)) := by
  rw [Set.disjoint_left]
  intro p hpi hpj
  rw [interior_unitSq] at hpi hpj
  obtain ⟨h1, h2⟩ := hpi; obtain ⟨h3, h4⟩ := hpj
  simp only at h1 h2 h3 h4
  have r : ∀ (dx dy θ : ℝ), (dx * cos θ + dy * sin θ) ^ 2 + (-dx * sin θ + dy * cos θ) ^ 2 = dx ^ 2 + dy ^ 2 :=
    fun dx dy θ => by linear_combination (dx ^ 2 + dy ^ 2) * cos_sq_add_sin_sq θ
  have e1 := r (p.1 - xi) (p.2 - yi) θi
  have e2 := r (p.1 - xj) (p.2 - yj) θj
  have s1 : (p.1 - xi) ^ 2 + (p.2 - yi) ^ 2 < 1 / 2 := by
    rw [← e1]; nlinarith [abs_lt.1 h1, abs_lt.1 h2, sq_abs ((p.1 - xi) * cos θi + (p.2 - yi) * sin θi)]
  have s2 : (p.1 - xj) ^ 2 + (p.2 - yj) ^ 2 < 1 / 2 := by
    rw [← e2]; nlinarith [abs_lt.1 h3, abs_lt.1 h4]
  nlinarith [sq_nonneg (p.1 - xi + (p.1 - xj)), sq_nonneg (p.2 - yi + (p.2 - yj))]

/-! ## 4. Scaling about the origin -/

/-- Scaling the centres by `ρ ≥ 1` adds at least `(ρ - 1) (g + 1)` to a pair margin `g`. -/
lemma gapAt_scale {xo yo θo xu yu θu ρ : ℝ} {k : Fin 4} (hρ : 1 ≤ ρ) :
    ρ * (gapAt xo yo θo xu yu θu k + 1 - 1 / ρ) ≤ gapAt (ρ * xo) (ρ * yo) θo (ρ * xu) (ρ * yu) θu k := by
  have hρ0 : 0 < ρ := by linarith
  have hh := one_half_le_hw (θu - θo)
  unfold gapAt
  have e : ρ * (1 / ρ) = 1 := by field_simp
  set d := (nrm (cos θo) (sin θo) k).1 * (xu - xo) + (nrm (cos θo) (sin θo) k).2 * (yu - yo)
  have e2 : (nrm (cos θo) (sin θo) k).1 * (ρ * xu - ρ * xo) + (nrm (cos θo) (sin θo) k).2 * (ρ * yu - ρ * yo) =
      ρ * d := by simp only [d]; ring
  rw [e2]
  nlinarith

/-- The static certificate of a configuration, at the scale `ρ`: every wall margin `m` has
`m + (1 - 1/ρ) / 2 ≥ 0`; every pair is far apart or has a side line with `g + 1 - 1/ρ ≥ 0`. -/
def MarginsOK (n : ℕ) (S : ℝ) (x y θ : Fin n → ℝ) (ρ : ℝ) : Prop :=
  (∀ i, 0 ≤ x i - hw (θ i) + (1 - 1 / ρ) / 2 ∧ 0 ≤ S - x i - hw (θ i) + (1 - 1 / ρ) / 2 ∧
        0 ≤ y i - hw (θ i) + (1 - 1 / ρ) / 2 ∧ 0 ≤ S - y i - hw (θ i) + (1 - 1 / ρ) / 2) ∧
  ∀ i j, i < j → 2 ≤ (x j - x i) ^ 2 + (y j - y i) ^ 2 ∨
    (∃ k, 0 ≤ gapAt (x i) (y i) (θ i) (x j) (y j) (θ j) k + 1 - 1 / ρ) ∨
    (∃ k, 0 ≤ gapAt (x j) (y j) (θ j) (x i) (y i) (θ i) k + 1 - 1 / ρ)

/-- **The scaled configuration is a packing**: centres times `ρ ≥ 1`, box side `ρ S`. -/
theorem isPacking_scale {n : ℕ} {S ρ : ℝ} {x y θ : Fin n → ℝ} (hρ : 1 ≤ ρ) (h : MarginsOK n S x y θ ρ) :
    IsPacking n (ρ * S) (fun i => (ρ * x i, ρ * y i)) θ := by
  have hρ0 : 0 < ρ := by linarith
  have e : ρ * (1 / ρ) = 1 := by field_simp
  refine ⟨fun i => ?_, ?_⟩
  · obtain ⟨h1, h2, h3, h4⟩ := h.1 i
    have hh := one_half_le_hw (θ i)
    refine subset_container_of_walls ?_ ?_ ?_ ?_
    · nlinarith [mul_le_mul_of_nonneg_left h1 hρ0.le]
    · nlinarith [mul_le_mul_of_nonneg_left h2 hρ0.le]
    · nlinarith [mul_le_mul_of_nonneg_left h3 hρ0.le]
    · nlinarith [mul_le_mul_of_nonneg_left h4 hρ0.le]
  · have pos : ∀ i j, i < j → Disjoint (interior (unitSq (ρ * x i, ρ * y i) (θ i)))
        (interior (unitSq (ρ * x j, ρ * y j) (θ j))) := by
      intro i j hij
      rcases h.2 i j hij with hf | ⟨k, hk⟩ | ⟨k, hk⟩
      · refine disjoint_of_far ?_
        have : (ρ * x j - ρ * x i) ^ 2 + (ρ * y j - ρ * y i) ^ 2 = ρ ^ 2 * ((x j - x i) ^ 2 + (y j - y i) ^ 2) := by
          ring
        rw [this]; nlinarith
      · exact disjoint_of_gapAt ((mul_nonneg hρ0.le hk).trans (gapAt_scale hρ))
      · exact (disjoint_of_gapAt ((mul_nonneg hρ0.le hk).trans (gapAt_scale hρ))).symm
    intro i j hij
    rcases lt_or_gt_of_ne hij with hl | hl
    · exact pos i j hl
    · exact (pos j i hl).symm

/-! ## 5. Steps and chains -/

/-- A configuration: box side, centres `(x i, y i)`, angles `θ i`. -/
abbrev Node (n : ℕ) := ℝ × (Fin n → ℝ) × (Fin n → ℝ) × (Fin n → ℝ)

/-- The configurations that are packings in a box of side at most `L`. -/
def Packed (n : ℕ) (L : ℝ) : Set (Node n) :=
  {C | C.1 ≤ L ∧ IsPacking n C.1 (fun i => (C.2.1 i, C.2.2.1 i)) C.2.2.2}

/-- Linear interpolation of two configurations. -/
def lerpN {n : ℕ} (A B : Node n) (t : ℝ) : Node n :=
  (A.1 + t * (B.1 - A.1), fun i => A.2.1 i + t * (B.2.1 i - A.2.1 i),
    fun i => A.2.2.1 i + t * (B.2.2.1 i - A.2.2.1 i), fun i => A.2.2.2 i + t * (B.2.2.2 i - A.2.2.2 i))

/-- A configuration scaled about the origin to box side `L`. -/
noncomputable def scaleTo {n : ℕ} (C : Node n) (L : ℝ) : Node n :=
  (L, fun i => L / C.1 * C.2.1 i, fun i => L / C.1 * C.2.2.1 i, C.2.2.2)

/-- The script's wall test: `(1 - t) m_A + t m_B - t (1 - t) (√2/2) Δθ² / 2 + (1 - s(t)/L) / 2 ≥ 0` on `[0, 1]`,
`m = c - hw θ` with `c` the centre's distance to the wall. -/
def WallCert (c₀ c₁ θ₀ θ₁ s₀ s₁ L : ℝ) : Prop :=
  ∀ t ∈ Icc (0 : ℝ) 1, 0 ≤ (1 - t) * (c₀ - hw θ₀) + t * (c₁ - hw θ₁) - t * (1 - t) * (√2 / 2 * (θ₁ - θ₀) ^ 2) / 2 +
    (1 - (s₀ + t * (s₁ - s₀)) / L) / 2

/-- The script's pair test for owner `o`, other square `u`, side `k`: bounds `W ≥ |w|` at both ends and
`D ≥ |Δw|`, and `(1 - t) g_A + t g_B - t (1 - t) q / 2 + 1 - s(t)/L ≥ 0` on `[0, 1]`. -/
def AxisCert {n : ℕ} (A B : Node n) (L : ℝ) (o u : Fin n) (k : Fin 4) : Prop :=
  ∃ W D : ℝ, 0 ≤ W ∧ 0 ≤ D ∧
    (A.2.1 u - A.2.1 o) ^ 2 + (A.2.2.1 u - A.2.2.1 o) ^ 2 ≤ W ^ 2 ∧
    (B.2.1 u - B.2.1 o) ^ 2 + (B.2.2.1 u - B.2.2.1 o) ^ 2 ≤ W ^ 2 ∧
    ((B.2.1 u - B.2.1 o) - (A.2.1 u - A.2.1 o)) ^ 2 + ((B.2.2.1 u - B.2.2.1 o) - (A.2.2.1 u - A.2.2.1 o)) ^ 2 ≤
      D ^ 2 ∧
    ∀ t ∈ Icc (0 : ℝ) 1, 0 ≤ (1 - t) * gapAt (A.2.1 o) (A.2.2.1 o) (A.2.2.2 o) (A.2.1 u) (A.2.2.1 u) (A.2.2.2 u) k +
      t * gapAt (B.2.1 o) (B.2.2.1 o) (B.2.2.2 o) (B.2.1 u) (B.2.2.1 u) (B.2.2.2 u) k -
      t * (1 - t) * (W * (B.2.2.2 o - A.2.2.2 o) ^ 2 + 2 * D * |B.2.2.2 o - A.2.2.2 o| +
        √2 / 2 * ((B.2.2.2 u - B.2.2.2 o) - (A.2.2.2 u - A.2.2.2 o)) ^ 2) / 2 +
      (1 - (A.1 + t * (B.1 - A.1)) / L)

/-- **The certificate of one step `A → B` at level `L`**, as `pathcert.py`'s `step_ok` checks it: positive sides
`≤ L`, the four wall tests for every square, and for every pair either centres `≥ √2` apart throughout or a pair
test with one of the two squares as owner. -/
def StepCert {n : ℕ} (A B : Node n) (L : ℝ) : Prop :=
  0 < A.1 ∧ 0 < B.1 ∧ A.1 ≤ L ∧ B.1 ≤ L ∧
  (∀ i, WallCert (A.2.1 i) (B.2.1 i) (A.2.2.2 i) (B.2.2.2 i) A.1 B.1 L ∧
        WallCert (A.1 - A.2.1 i) (B.1 - B.2.1 i) (A.2.2.2 i) (B.2.2.2 i) A.1 B.1 L ∧
        WallCert (A.2.2.1 i) (B.2.2.1 i) (A.2.2.2 i) (B.2.2.2 i) A.1 B.1 L ∧
        WallCert (A.1 - A.2.2.1 i) (B.1 - B.2.2.1 i) (A.2.2.2 i) (B.2.2.2 i) A.1 B.1 L) ∧
  ∀ i j, i < j →
    (∀ t ∈ Icc (0 : ℝ) 1, 2 ≤ ((lerpN A B t).2.1 j - (lerpN A B t).2.1 i) ^ 2 +
      ((lerpN A B t).2.2.1 j - (lerpN A B t).2.2.1 i) ^ 2) ∨
    (∃ k, AxisCert A B L i j k) ∨ (∃ k, AxisCert A B L j i k)

/-- The side normals of a square at angle `θ` are the directions `θ + k π/2`. -/
lemma nrm_angle (θ : ℝ) (k : Fin 4) :
    nrm (cos θ) (sin θ) k = (cos (θ + (k : ℕ) * (π / 2)), sin (θ + (k : ℕ) * (π / 2))) := by
  fin_cases k
  · simp [nrm]
  · simp [nrm, cos_add_pi_div_two, sin_add_pi_div_two]
  · simp only [nrm, Nat.cast_ofNat]
    rw [show θ + 2 * (π / 2) = θ + π by ring, cos_add_pi, sin_add_pi]
  · simp only [nrm]
    rw [show θ + ((3 : ℕ) : ℝ) * (π / 2) = (θ + π) + π / 2 by push_cast; ring, cos_add_pi_div_two,
      sin_add_pi_div_two, cos_add_pi, sin_add_pi]
    simp

/-- The pair margin of the interpolated configuration is `pairGap`. -/
lemma gapAt_lerp (xo yo θo xu yu θu xo' yo' θo' xu' yu' θu' : ℝ) (k : Fin 4) (t : ℝ) :
    gapAt (xo + t * (xo' - xo)) (yo + t * (yo' - yo)) (θo + t * (θo' - θo)) (xu + t * (xu' - xu))
      (yu + t * (yu' - yu)) (θu + t * (θu' - θu)) k =
    pairGap (θo + (k : ℕ) * (π / 2)) (θo' - θo) (xu - xo) ((xu' - xo') - (xu - xo)) (yu - yo) ((yu' - yo') - (yu - yo))
      (θu - θo) ((θu' - θo') - (θu - θo)) t := by
  unfold gapAt pairGap
  rw [nrm_angle]
  simp only
  rw [show θo + t * (θo' - θo) + (k : ℕ) * (π / 2) = θo + (k : ℕ) * (π / 2) + t * (θo' - θo) by ring,
    show θu + t * (θu' - θu) - (θo + t * (θo' - θo)) = θu - θo + t * (θu' - θo' - (θu - θo)) by ring]
  ring

/-- **One step**: at every `t ∈ [0, 1]`, the interpolated configuration scaled to side `L` is a packing. -/
theorem step_isPacking {n : ℕ} {A B : Node n} {L : ℝ} (h : StepCert A B L) {t : ℝ} (ht : t ∈ Icc (0 : ℝ) 1) :
    scaleTo (lerpN A B t) L ∈ Packed n L := by
  obtain ⟨hA, hB, hAL, hBL, hwall, hpair⟩ := h
  set s := A.1 + t * (B.1 - A.1) with hs
  have hs0 : 0 < s := by
    rw [hs]; nlinarith [mul_nonneg ht.1 hB.le, mul_nonneg (sub_nonneg.2 ht.2) hA.le, ht.1, ht.2]
  have hsL : s ≤ L := by
    rw [hs]; nlinarith [mul_le_mul_of_nonneg_left hBL ht.1, mul_le_mul_of_nonneg_left hAL (sub_nonneg.2 ht.2)]
  have hL : 0 < L := by linarith
  set ρ := L / s with hρ
  have hρ1 : 1 ≤ ρ := by rw [hρ, le_div_iff₀ hs0]; linarith
  have hinv : 1 - 1 / ρ = 1 - s / L := by rw [hρ, one_div_div]
  have hρs : ρ * s = L := by rw [hρ]; field_simp
  -- a wall test gives the wall margin at t
  have wall : ∀ {c₀ c₁ θ₀ θ₁ : ℝ}, WallCert c₀ c₁ θ₀ θ₁ A.1 B.1 L →
      0 ≤ c₀ + t * (c₁ - c₀) - hw (θ₀ + t * (θ₁ - θ₀)) + (1 - 1 / ρ) / 2 := by
    intro c₀ c₁ θ₀ θ₁ hc
    have h1 := hc t ht
    have h2 := wallGap_chord c₀ (c₁ - c₀) θ₀ (θ₁ - θ₀) ht
    unfold wallGap at h2
    simp only [zero_mul, add_zero, one_mul, add_sub_cancel] at h2
    rw [hinv]
    linarith
  -- a pair test gives the pair margin at t
  have axis : ∀ {o u : Fin n} {k : Fin 4}, AxisCert A B L o u k →
      0 ≤ gapAt ((lerpN A B t).2.1 o) ((lerpN A B t).2.2.1 o) ((lerpN A B t).2.2.2 o) ((lerpN A B t).2.1 u)
        ((lerpN A B t).2.2.1 u) ((lerpN A B t).2.2.2 u) k + 1 - 1 / ρ := by
    intro o u k ⟨W, D, hW0, hD0, hWA, hWB, hD, hc⟩
    have h1 := hc t ht
    have hWB' : ((A.2.1 u - A.2.1 o) + ((B.2.1 u - B.2.1 o) - (A.2.1 u - A.2.1 o))) ^ 2 +
        ((A.2.2.1 u - A.2.2.1 o) + ((B.2.2.1 u - B.2.2.1 o) - (A.2.2.1 u - A.2.2.1 o))) ^ 2 ≤ W ^ 2 := by
      simpa only [add_sub_cancel] using hWB
    have hWt := norm_sq_le_of_ends hWA hWB'
    have h2 := pairGap_chord (α₀ := A.2.2.2 o + (k : ℕ) * (π / 2)) (dα := B.2.2.2 o - A.2.2.2 o)
      (φ₀ := A.2.2.2 u - A.2.2.2 o) (dφ := (B.2.2.2 u - B.2.2.2 o) - (A.2.2.2 u - A.2.2.2 o)) hWt hW0 hD hD0 ht
    have e0 := gapAt_lerp (A.2.1 o) (A.2.2.1 o) (A.2.2.2 o) (A.2.1 u) (A.2.2.1 u) (A.2.2.2 u) (B.2.1 o)
      (B.2.2.1 o) (B.2.2.2 o) (B.2.1 u) (B.2.2.1 u) (B.2.2.2 u) k 0
    have e1 := gapAt_lerp (A.2.1 o) (A.2.2.1 o) (A.2.2.2 o) (A.2.1 u) (A.2.2.1 u) (A.2.2.2 u) (B.2.1 o)
      (B.2.2.1 o) (B.2.2.2 o) (B.2.1 u) (B.2.2.1 u) (B.2.2.2 u) k 1
    have et := gapAt_lerp (A.2.1 o) (A.2.2.1 o) (A.2.2.2 o) (A.2.1 u) (A.2.2.1 u) (A.2.2.2 u) (B.2.1 o)
      (B.2.2.1 o) (B.2.2.2 o) (B.2.1 u) (B.2.2.1 u) (B.2.2.2 u) k t
    simp only [zero_mul, add_zero, one_mul, add_sub_cancel] at e0 e1
    simp only [lerpN]
    have hinv' : 1 / ρ = s / L := by rw [hρ, one_div_div]
    rw [et, hinv', hs]
    rw [← e0, ← e1] at h2
    linarith
  have hm : MarginsOK n s (lerpN A B t).2.1 (lerpN A B t).2.2.1 (lerpN A B t).2.2.2 ρ := by
    refine ⟨fun i => ?_, fun i j hij => ?_⟩
    · obtain ⟨w1, w2, w3, w4⟩ := hwall i
      refine ⟨wall w1, ?_, wall w3, ?_⟩
      · have := wall w2; simp only [lerpN, hs] at this ⊢; convert this using 2; ring
      · have := wall w4; simp only [lerpN, hs] at this ⊢; convert this using 2; ring
    · rcases hpair i j hij with hf | ⟨k, hk⟩ | ⟨k, hk⟩
      · exact Or.inl (hf t ht)
      · exact Or.inr (Or.inl ⟨k, axis hk⟩)
      · exact Or.inr (Or.inr ⟨k, axis hk⟩)
  have hp := isPacking_scale hρ1 hm
  rw [hρs] at hp
  exact ⟨le_refl L, hp⟩

/-- **One step is a path inside the packings of side `≤ L`**, from `A` scaled to side `L` to `B` scaled to side
`L`. -/
theorem step_joined {n : ℕ} {A B : Node n} {L : ℝ} (h : StepCert A B L) :
    JoinedIn (Packed n L) (scaleTo A L) (scaleTo B L) := by
  have hA := h.1; have hB := h.2.1
  have hpos : ∀ t ∈ unitInterval, A.1 + t * (B.1 - A.1) ≠ 0 := fun t ht => by
    have h1 := ht.1; have h2 := ht.2
    nlinarith [mul_nonneg h1 hB.le, mul_nonneg (sub_nonneg.2 h2) hA.le]
  refine JoinedIn.ofLine (f := fun t => scaleTo (lerpN A B t) L) ?_ ?_ ?_ ?_
  · refine continuousOn_const.prodMk (ContinuousOn.prodMk ?_ (ContinuousOn.prodMk ?_ ?_))
    · refine continuousOn_pi.2 fun i => ?_
      exact (continuousOn_const.div (by fun_prop) hpos).mul (by fun_prop)
    · refine continuousOn_pi.2 fun i => ?_
      exact (continuousOn_const.div (by fun_prop) hpos).mul (by fun_prop)
    · exact continuousOn_pi.2 fun i => by simp only [lerpN]; fun_prop
  · simp [scaleTo, lerpN]
  · simp only [scaleTo, lerpN, one_mul, add_sub_cancel]
  · rintro _ ⟨t, ht, rfl⟩
    exact step_isPacking h ht

/-- **A certified chain** `N 0 → N 1 → … → N m`: its ends, scaled to side `L`, are joined inside the packings of
side `≤ L`. -/
theorem chain_joined {n : ℕ} {L : ℝ} (m : ℕ) (N : ℕ → Node n) (h : ∀ i < m, StepCert (N i) (N (i + 1)) L)
    (h0 : scaleTo (N 0) L ∈ Packed n L) : JoinedIn (Packed n L) (scaleTo (N 0) L) (scaleTo (N m) L) := by
  induction m with
  | zero => exact JoinedIn.refl h0
  | succ m ih =>
    exact (ih fun i hi => h i (by omega)).trans (step_joined (h m (by omega)))

/-! ## 6. The closed-form check -/

/-- The script's `_cond`: the test's left side is the quadratic `c₀ + c₁ t + c₂ t²` with
`c₀ = g_A + k (1 - s_A/L)`, `c₁ = (g_B - g_A) - q/2 - k (s_B - s_A)/L`, `c₂ = q/2`. -/
lemma cond_eq (gA gB q k sA sB L t : ℝ) :
    (1 - t) * gA + t * gB - t * (1 - t) * q / 2 + k * (1 - (sA + t * (sB - sA)) / L) =
      (gA + k * (1 - sA / L)) + ((gB - gA) - q / 2 - k * (sB - sA) / L) * t + q / 2 * t ^ 2 := by
  ring

/-- The script's `_qmin`: a quadratic with `c₂ ≥ 0` that is `≥ 0` at `0`, at `1`, and at its vertex when the vertex
lies in `(0, 1)`, is `≥ 0` on `[0, 1]`. -/
lemma qmin_nonneg {c₀ c₁ c₂ : ℝ} (h2 : 0 ≤ c₂) (h0 : 0 ≤ c₀) (h1 : 0 ≤ c₀ + c₁ + c₂)
    (hv : 0 < c₂ → 0 < -c₁ / (2 * c₂) → -c₁ / (2 * c₂) < 1 →
      0 ≤ c₀ + c₁ * (-c₁ / (2 * c₂)) + c₂ * (-c₁ / (2 * c₂)) ^ 2) :
    ∀ t ∈ Icc (0 : ℝ) 1, 0 ≤ c₀ + c₁ * t + c₂ * t ^ 2 := by
  intro t ⟨ht0, ht1⟩
  rcases h2.lt_or_eq with hc | hc
  · set v := -c₁ / (2 * c₂) with hvdef
    have hc1 : c₁ = -2 * c₂ * v := by rw [hvdef]; field_simp
    -- `c₀ + c₁ t + c₂ t² = c₂ (t - v)² + (c₀ - c₂ v²)`
    have e : ∀ x, c₀ + c₁ * x + c₂ * x ^ 2 = c₂ * (x - v) ^ 2 + (c₀ - c₂ * v ^ 2) := fun x => by rw [hc1]; ring
    rw [e]
    by_cases hv0 : v ≤ 0
    · -- increasing on [0, 1]: at least the value at 0
      nlinarith [mul_nonneg hc.le (mul_nonneg ht0 (by linarith : 0 ≤ t - 2 * v))]
    by_cases hv1 : 1 ≤ v
    · -- decreasing on [0, 1]: at least the value at 1
      rw [hc1] at h1
      nlinarith [mul_nonneg hc.le (mul_nonneg (by linarith : 0 ≤ 1 - t) (by linarith : 0 ≤ 2 * v - 1 - t))]
    · have := hv hc (lt_of_not_ge hv0) (lt_of_not_ge hv1)
      rw [e] at this
      nlinarith [mul_nonneg hc.le (sq_nonneg (t - v))]
  · subst hc
    nlinarith [mul_nonneg ht0 (sub_nonneg.2 ht1)]

/-! ## 7. Ends: the packings themselves -/

lemma MarginsOK.mono {n : ℕ} {S : ℝ} {x y θ : Fin n → ℝ} {ρ ρ' : ℝ} (hρ : 0 < ρ) (hρρ' : ρ ≤ ρ')
    (h : MarginsOK n S x y θ ρ) : MarginsOK n S x y θ ρ' := by
  have hle : 1 / ρ' ≤ 1 / ρ := one_div_le_one_div_of_le hρ hρρ'
  refine ⟨fun i => ?_, fun i j hij => ?_⟩
  · obtain ⟨h1, h2, h3, h4⟩ := h.1 i
    exact ⟨by linarith, by linarith, by linarith, by linarith⟩
  · rcases h.2 i j hij with hf | ⟨k, hk⟩ | ⟨k, hk⟩
    · exact Or.inl hf
    · exact Or.inr (Or.inl ⟨k, by linarith⟩)
    · exact Or.inr (Or.inr ⟨k, by linarith⟩)

/-- A packing certified at its own side (margins `≥ 0`: `MarginsOK` at `ρ = 1`) is joined to itself scaled up to
side `L ≥` its side, inside the packings of side `≤ L`. -/
theorem scale_joined {n : ℕ} {P : Node n} {L : ℝ} (hP : 0 < P.1) (hPL : P.1 ≤ L)
    (hm : MarginsOK n P.1 P.2.1 P.2.2.1 P.2.2.2 1) : JoinedIn (Packed n L) P (scaleTo P L) := by
  have hr : 1 ≤ L / P.1 := by rw [le_div_iff₀ hP]; linarith
  refine JoinedIn.ofLine (f := fun t => ((1 + t * (L / P.1 - 1)) * P.1,
    fun i => (1 + t * (L / P.1 - 1)) * P.2.1 i, fun i => (1 + t * (L / P.1 - 1)) * P.2.2.1 i, P.2.2.2)) ?_ ?_ ?_ ?_
  · fun_prop
  · simp
  · simp only [scaleTo, one_mul, add_sub_cancel]
    rw [div_mul_cancel₀ L hP.ne']
  · rintro _ ⟨t, ⟨ht0, ht1⟩, rfl⟩
    have hμ : 1 ≤ 1 + t * (L / P.1 - 1) := by nlinarith
    refine ⟨?_, isPacking_scale hμ (hm.mono one_pos hμ)⟩
    have : (1 + t * (L / P.1 - 1)) * P.1 = P.1 + t * (L - P.1) := by field_simp
    show (1 + t * (L / P.1 - 1)) * P.1 ≤ L
    rw [this]; nlinarith

/-- **Same basin.**  A certified chain `N 0 → … → N m` at level `L` whose two ends are packings (certified at their
own sides) joins them by a path of packings, every one in a box of side `≤ L`. -/
theorem joined_of_chain {n : ℕ} {L : ℝ} (m : ℕ) (N : ℕ → Node n) (h : ∀ i < m, StepCert (N i) (N (i + 1)) L)
    (h0 : 0 < (N 0).1) (hm0 : 0 < (N m).1) (hL0 : (N 0).1 ≤ L) (hLm : (N m).1 ≤ L)
    (hP : MarginsOK n (N 0).1 (N 0).2.1 (N 0).2.2.1 (N 0).2.2.2 1)
    (hQ : MarginsOK n (N m).1 (N m).2.1 (N m).2.2.1 (N m).2.2.2 1) :
    JoinedIn (Packed n L) (N 0) (N m) := by
  have a := scale_joined h0 hL0 hP
  have b := scale_joined hm0 hLm hQ
  exact (a.trans (chain_joined m N h a.target_mem)).trans b.symm

end UnitSquarePacking
