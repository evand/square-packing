import Sqpack.StripBand

/-!
# One square, one line: where the extremes are

For opposite slopes the step cost splits into a left part and a right part (`kap_split`):
`kap C S H x (-y) = C + S + kh C S H x + kh S C H y`, where `kh c s H z = (s - H/2) z + (c z - s)⁺` is what one
square contributes against one line of slope `z ≥ 0`.  So each band condition is a statement about a left line and a
right line separately.

For one line, the potential plus `kh` is a minimum of convex functions with one kink each (and the potential minus
`kh` a maximum of concave ones).  So its extreme over the slope is at the start of the range, at a breakpoint of the
potential (`8/7`, `x₂`, `y₁`), or at the slope where `kh` has its kink: `z = s / c`, the line that runs **along a
side of the square** (the *aligned* line).  The line lemmas below state exactly this:

* `line_lo`: `min (p + βx) (q + L x) + kh` (both potentials are of this form) on `[X₀, ∞)`: least at `X₀`, `8/7`,
  or the aligned slope;
* `line_u1_hi`: `u₁ - kh` on `[0, ∞)`: greatest at `0`, `y₁`, or the aligned slope;
* `line_u0_hi`: `u₀ - kh` on `[0, ∞)`: greatest at `0`, `x₂`, `8/7`, or the aligned slope.

The aligned slope is passed as "any `t` with `c t = s`", so that no division appears and an upright square
(`c = 0`, no kink) needs no special case.
-/

namespace SquarePacking

/-- One square against one line of slope `z ≥ 0` (frame `(c, s)`, measured from the line's side). -/
noncomputable def kh (c s H z : ℝ) : ℝ := (s - H / 2) * z + max (c * z - s) 0

/-- **The step cost splits** into a left part and a right part. -/
lemma kap_split {C S H x y : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hx : 0 ≤ x) (hy : 0 ≤ y) :
    kap C S H x (-y) = C + S + kh C S H x + kh S C H y := by
  rw [kap_opp hC hS hx hy]; unfold kh; ring

/-! ## A function with one kink -/

/-- `A + B z + (c z - s)⁺` on `[lo, hi]` is least at an end or at its kink. -/
lemma kink_lo {A B c s lo hi x m : ℝ} (hc : 0 ≤ c) (hl : lo ≤ x) (hh : x ≤ hi)
    (h1 : m ≤ A + B * lo + max (c * lo - s) 0) (h2 : m ≤ A + B * hi + max (c * hi - s) 0)
    (hk : ∀ t, lo ≤ t → t ≤ hi → c * t = s → m ≤ A + B * t) :
    m ≤ A + B * x + max (c * x - s) 0 := by
  rcases le_total (c * x) s with hx | hx
  · rw [max_eq_right (by linarith)]
    rcases le_total 0 B with hB | hB
    · have : c * lo ≤ s := by nlinarith [mul_le_mul_of_nonneg_left hl hc]
      rw [max_eq_right (by linarith)] at h1; nlinarith [mul_le_mul_of_nonneg_left hl hB]
    · rcases le_or_gt (c * hi) s with hh' | hh'
      · rw [max_eq_right (by linarith)] at h2; nlinarith [mul_le_mul_of_nonpos_left hh hB]
      · have hc0 : 0 < c := by
          rcases eq_or_lt_of_le hc with h0 | h0
          · rw [← h0] at hx hh'; linarith
          · exact h0
        have ht := hk (s / c) (by rw [le_div_iff₀ hc0]; nlinarith [mul_le_mul_of_nonneg_left hl hc])
          (by rw [div_le_iff₀ hc0]; linarith) (by field_simp)
        have : x ≤ s / c := by rw [le_div_iff₀ hc0]; linarith
        nlinarith [mul_le_mul_of_nonpos_left this hB]
  · rw [max_eq_left (by linarith)]
    rcases le_total (B + c) 0 with hB | hB
    · have : s ≤ c * hi := by nlinarith [mul_le_mul_of_nonneg_left hh hc]
      rw [max_eq_left (by linarith)] at h2; nlinarith [mul_le_mul_of_nonpos_left hh hB]
    · rcases le_or_gt s (c * lo) with hl' | hl'
      · rw [max_eq_left (by linarith)] at h1; nlinarith [mul_le_mul_of_nonneg_left hl hB]
      · have hc0 : 0 < c := by
          rcases eq_or_lt_of_le hc with h0 | h0
          · rw [← h0] at hx hl'; linarith
          · exact h0
        have hsx : s / c ≤ x := by rw [div_le_iff₀ hc0]; linarith
        have ht := hk (s / c) (by rw [le_div_iff₀ hc0]; linarith) (by linarith) (by field_simp)
        have e : c * (s / c) = s := by field_simp
        nlinarith [mul_le_mul_of_nonneg_left hsx hB]

/-- `A + B z + (c z - s)⁺` on `[lo, ∞)` with `B + c ≥ 0`: least at `lo` or at its kink. -/
lemma kink_lo_inf {A B c s lo x m : ℝ} (hc : 0 ≤ c) (hl : lo ≤ x) (hs : 0 ≤ B + c)
    (h1 : m ≤ A + B * lo + max (c * lo - s) 0)
    (hk : ∀ t, lo ≤ t → c * t = s → m ≤ A + B * t) :
    m ≤ A + B * x + max (c * x - s) 0 := by
  rcases le_total (c * x) s with hx | hx
  · rcases le_total 0 B with hB | hB
    · rw [max_eq_right (by linarith)]
      have : c * lo ≤ s := by nlinarith [mul_le_mul_of_nonneg_left hl hc]
      rw [max_eq_right (by linarith)] at h1; nlinarith [mul_le_mul_of_nonneg_left hl hB]
    · -- `B ≤ 0 ≤ B + c` and `c x ≤ s`: use the kink `s / c ≥ x` (if `c = 0` then `B = 0`)
      rcases eq_or_lt_of_le hc with h0 | hc0
      · have hB0 : B = 0 := by linarith
        rw [← h0] at h1 ⊢; rw [hB0] at h1 ⊢; simpa using h1
      · exact kink_lo hc hl (le_refl _) h1 (by
          rw [max_eq_right (by linarith)]
          have ht := hk (s / c) (by rw [le_div_iff₀ hc0]; nlinarith [mul_le_mul_of_nonneg_left hl hc])
            (by field_simp)
          have : x ≤ s / c := by rw [le_div_iff₀ hc0]; linarith
          nlinarith [mul_le_mul_of_nonpos_left this hB])
          (fun t h1t h2t e => hk t h1t e)
  · exact kink_lo hc hl (le_refl _) h1 (by
      rcases le_or_gt s (c * lo) with hl' | hl'
      · rw [max_eq_left (by linarith)] at h1
        rw [max_eq_left (by linarith)]; nlinarith [mul_le_mul_of_nonneg_left hl hs]
      · have hc0 : 0 < c := by
          rcases eq_or_lt_of_le hc with h0 | h0
          · rw [← h0] at hx hl'; linarith
          · exact h0
        have hsx : s / c ≤ x := by rw [div_le_iff₀ hc0]; linarith
        have ht := hk (s / c) (by rw [le_div_iff₀ hc0]; linarith) (by field_simp)
        have e : c * (s / c) = s := by field_simp
        rw [max_eq_left (by linarith)]; nlinarith [mul_le_mul_of_nonneg_left hsx hs])
      (fun t h1t h2t e => hk t h1t e)

/-! ## The potentials, piece by piece -/

section pieces
variable {H : ℝ}

/-- `x₂`: where `A0 = 2a + βx` meets `Lz`. -/
noncomputable def bandX2 (H : ℝ) : ℝ :=
  (2 * (1 - H / 2) - (1 - (2 * Real.sqrt 2 - H))) / (2 * (1 - H / 2))

/-- `y₁`: where `a + βx` meets `a + T`; `y₁ ≥ 8/7` exactly when `H ≥ H*(9)`. -/
noncomputable def bandY1 (H : ℝ) : ℝ :=
  ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) - (H / 16 - 1 / 10) * (8 / 7)) /
    ((1 - H / 2) - (H / 16 - 1 / 10))

lemma bandX2_eq (hb : InBand H) :
    bandX2 H * (2 * (1 - H / 2)) = 2 * (1 - H / 2) - (1 - (2 * Real.sqrt 2 - H)) := by
  have : 2 * (1 - H / 2) ≠ 0 := by linarith [hb.2]
  unfold bandX2; exact div_mul_cancel₀ _ this

lemma bandY1_eq (hb : InBand H) :
    bandY1 H * ((1 - H / 2) - (H / 16 - 1 / 10)) =
      (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) - (H / 16 - 1 / 10) * (8 / 7) := by
  have : (1 - H / 2) - (H / 16 - 1 / 10) ≠ 0 := by linarith [hb.2]
  unfold bandY1; exact div_mul_cancel₀ _ this

/-- `T - Lz = (T₁ + β)(x - 8/7)` with `T₁ + β > 0`: below `8/7` the max is `Lz`, above it is `T`. -/
lemma bandT_sub (H x : ℝ) : bandT H x - bandLz H x = (9 / 10 - 7 * H / 16) * (x - 8 / 7) := by
  unfold bandT bandLz; ring

lemma bandL_lo (hb : InBand H) {x : ℝ} (h : x ≤ 8 / 7) : bandL H x = bandLz H x := by
  unfold bandL; apply max_eq_left
  have := bandT_sub H x
  nlinarith [mul_le_mul_of_nonneg_left (by linarith : x - 8 / 7 ≤ 0) (by linarith [hb.2] : (0 : ℝ) ≤ 9 / 10 - 7 * H / 16)]

lemma bandL_hi (hb : InBand H) {x : ℝ} (h : 8 / 7 ≤ x) : bandL H x = bandT H x := by
  unfold bandL; apply max_eq_right
  have := bandT_sub H x
  nlinarith [mul_le_mul_of_nonneg_left (by linarith : 0 ≤ x - 8 / 7) (by linarith [hb.2] : (0 : ℝ) ≤ 9 / 10 - 7 * H / 16)]


lemma bandX2_nonneg (hb : InBand H) : 0 ≤ bandX2 H := by
  have r1 := sqrt2_lo'
  have hd : 0 < 2 * (1 - H / 2) := by linarith [hb.2]
  have : 0 ≤ bandX2 H * (2 * (1 - H / 2)) := by rw [bandX2_eq hb]; linarith [hb.1, hb.2]
  exact nonneg_of_mul_nonneg_left this hd

lemma bandX2_le (hb : InBand H) : bandX2 H ≤ 8 / 7 := by
  have r1 := sqrt2_lo'
  have hd : 0 < 2 * (1 - H / 2) := by linarith [hb.2]
  have : bandX2 H * (2 * (1 - H / 2)) ≤ 8 / 7 * (2 * (1 - H / 2)) := by rw [bandX2_eq hb]; linarith [hb.1, hb.2]
  exact le_of_mul_le_mul_right this hd

lemma bandY1_nonneg (hb : InBand H) : 0 ≤ bandY1 H := by
  have r2 := sqrt2_hi'
  have hd : 0 < (1 - H / 2) - (H / 16 - 1 / 10) := by linarith [hb.2]
  have : 0 ≤ bandY1 H * ((1 - H / 2) - (H / 16 - 1 / 10)) := by rw [bandY1_eq hb]; linarith [hb.1, hb.2]
  exact nonneg_of_mul_nonneg_left this hd

lemma half_le_bandX2 (hb : InBand H) : 1 / 2 ≤ bandX2 H := by
  have r1 := sqrt2_lo'
  have hd : 0 < 2 * (1 - H / 2) := by linarith [hb.2]
  have : 1 / 2 * (2 * (1 - H / 2)) ≤ bandX2 H * (2 * (1 - H / 2)) := by rw [bandX2_eq hb]; linarith [hb.1, hb.2]
  exact le_of_mul_le_mul_right this hd

/-- `y₁ ≥ 8/7` exactly when `H ≥ H*(9)`: this is where the band starts. -/
lemma le_bandY1 (hb : InBand H) : 8 / 7 ≤ bandY1 H := by
  have hd : 0 < (1 - H / 2) - (H / 16 - 1 / 10) := by linarith [hb.2]
  have : 8 / 7 * ((1 - H / 2) - (H / 16 - 1 / 10)) ≤ bandY1 H * ((1 - H / 2) - (H / 16 - 1 / 10)) := by
    rw [bandY1_eq hb]; linarith [hb.1]
  exact le_of_mul_le_mul_right this hd

/-! The potentials on their pieces, in the form used by the pair table. -/

lemma bandU0_A0' (hb : InBand H) {z : ℝ} (hz : 0 ≤ z) (h : z ≤ bandX2 H) :
    bandU0 H z = 2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z :=
  u0_A0 hb hz (by rw [← bandX2_eq hb]; exact mul_le_mul_of_nonneg_right h (by linarith [hb.2]))

lemma bandU0_Lz' (hb : InBand H) {z : ℝ} (h : bandX2 H ≤ z) (h' : z ≤ 8 / 7) :
    bandU0 H z = (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - z) :=
  u0_Lz hb ((bandX2_nonneg hb).trans h) (by rw [← bandX2_eq hb]; exact mul_le_mul_of_nonneg_right h (by linarith [hb.2]))
    (by linarith)

lemma bandU0_T' (hb : InBand H) {z : ℝ} (h : 8 / 7 ≤ z) :
    bandU0 H z = (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * (z - 8 / 7) :=
  u0_T hb (by linarith)

lemma bandU1_O' (hb : InBand H) {z : ℝ} (hz : 0 ≤ z) (h : z ≤ bandY1 H) :
    bandU1 H z = (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z :=
  u1_O hb hz (by rw [← bandY1_eq hb]; exact mul_le_mul_of_nonneg_right h (by linarith [hb.2]))

lemma bandU1_aT' (hb : InBand H) {z : ℝ} (h : bandY1 H ≤ z) :
    bandU1 H z = (1 - (2 * Real.sqrt 2 - H)) +
      ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * (z - 8 / 7)) :=
  u1_aT hb (by rw [← bandY1_eq hb]; exact mul_le_mul_of_nonneg_right h (by linarith [hb.2]))

/-- At the aligned slope the kink term vanishes. -/
lemma kh_aligned {c s H z : ℝ} (h : c * z = s) : kh c s H z = (s - H / 2) * z := by
  unfold kh; rw [h, sub_self, max_self, add_zero]

lemma kh_pos {c s H z : ℝ} (h : s ≤ c * z) : kh c s H z = (s - H / 2) * z + (c * z - s) := by
  unfold kh; rw [max_eq_left (by linarith)]

lemma kh_neg {c s H z : ℝ} (h : c * z ≤ s) : kh c s H z = (s - H / 2) * z := by
  unfold kh; rw [max_eq_right (by linarith), add_zero]

end pieces

/-! ## The line lemmas -/

section lines
variable {H C S : ℝ}

/-- **Lower bound along one line.**  For `min (p + βx) (q + L x) + kh` on `[X₀, ∞)` (`u₀`: `p = 2a`, `q = 0`;
`u₁`: `p = q = a`), it suffices to check `X₀`, `8/7` and the aligned slope. -/
lemma line_lo (hb : InBand H) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) {p q X₀ m : ℝ}
    (hX : X₀ ≤ 8 / 7)
    (h0 : m ≤ min (p + (1 - H / 2) * X₀) (q + bandL H X₀) + kh C S H X₀)
    (h1 : m ≤ min (p + (1 - H / 2) * (8 / 7)) (q + bandL H (8 / 7)) + kh C S H (8 / 7))
    (ht : ∀ t, X₀ ≤ t → C * t = S → m ≤ min (p + (1 - H / 2) * t) (q + bandL H t) + kh C S H t)
    {x : ℝ} (hx : X₀ ≤ x) : m ≤ min (p + (1 - H / 2) * x) (q + bandL H x) + kh C S H x := by
  obtain ⟨h9, h87⟩ := hb
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  rw [← min_add_add_right]
  apply le_min
  · -- the straight branch `p + βx`: one kink
    have := kink_lo_inf (A := p) (B := (1 - H / 2) + (S - H / 2)) (c := C) (s := S) (lo := X₀) (m := m) hC hx
      (by linarith)
      (by have := h0.trans (add_le_add_left (min_le_left _ _) _); unfold kh at this; linarith)
      (fun t h1t e => by
        have := (ht t h1t e).trans (add_le_add_left (min_le_left _ _) _)
        unfold kh at this; rw [e, sub_self, max_self] at this; linarith)
    unfold kh; linarith
  · -- the branch `q + L x`: `L = Lz` below `8/7`, `T` above
    rcases le_total x (8 / 7) with hx8 | hx8
    · have := kink_lo (A := q + bandLz H 0) (B := -(1 - H / 2) + (S - H / 2)) (c := C) (s := S) (lo := X₀)
        (hi := 8 / 7) (m := m) hC hx hx8
        (by have := h0.trans (add_le_add_left (min_le_right _ _) _)
            rw [bandL_lo ⟨h9, h87⟩ hX] at this; unfold kh bandLz at this; unfold bandLz; linarith)
        (by have := h1.trans (add_le_add_left (min_le_right _ _) _)
            rw [bandL_lo ⟨h9, h87⟩ le_rfl] at this; unfold kh bandLz at this; unfold bandLz; linarith)
        (fun t h1t h2t e => by
          have := (ht t h1t e).trans (add_le_add_left (min_le_right _ _) _)
          rw [bandL_lo ⟨h9, h87⟩ h2t] at this; unfold kh bandLz at this; rw [e, sub_self, max_self] at this
          unfold bandLz; linarith)
      rw [bandL_lo ⟨h9, h87⟩ hx8]; unfold kh bandLz at *; linarith
    · have := kink_lo_inf (A := q + bandT H 0) (B := (H / 16 - 1 / 10) + (S - H / 2)) (c := C) (s := S)
        (lo := 8 / 7) (m := m) hC hx8 (by linarith)
        (by have := h1.trans (add_le_add_left (min_le_right _ _) _)
            rw [bandL_hi ⟨h9, h87⟩ le_rfl] at this; unfold kh bandT at this; unfold bandT; linarith)
        (fun t h1t e => by
          have := (ht t (hX.trans h1t) e).trans (add_le_add_left (min_le_right _ _) _)
          rw [bandL_hi ⟨h9, h87⟩ h1t] at this; unfold kh bandT at this; rw [e, sub_self, max_self] at this
          unfold bandT; linarith)
      rw [bandL_hi ⟨h9, h87⟩ hx8]; unfold kh bandT at *; linarith

/-- **Upper bound along one line, odd potential.**  `u₁ - kh` on `[0, ∞)` (`u₁ = a + βy` up to `y₁`, then `a + T`):
it suffices to check `0`, `y₁` and the aligned slope.  (Right-hand square: frame `(S, C)`.) -/
lemma line_u1_hi (hb : InBand H) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) {M : ℝ}
    (h0 : bandU1 H 0 - kh S C H 0 ≤ M) (h1 : bandU1 H (bandY1 H) - kh S C H (bandY1 H) ≤ M)
    (ht : ∀ t, 0 ≤ t → S * t = C → bandU1 H t - kh S C H t ≤ M)
    {y : ℝ} (hy : 0 ≤ y) : bandU1 H y - kh S C H y ≤ M := by
  obtain ⟨h9, h87⟩ := hb
  have hb := (⟨h9, h87⟩ : InBand H)
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hd : 0 < (1 - H / 2) - (H / 16 - 1 / 10) := by linarith
  have hY := bandY1_eq hb
  have hY0 : 0 ≤ bandY1 H := by
    have : 0 ≤ bandY1 H * ((1 - H / 2) - (H / 16 - 1 / 10)) := by rw [hY]; linarith
    exact nonneg_of_mul_nonneg_left this hd
  have pO : ∀ z, 0 ≤ z → z ≤ bandY1 H → bandU1 H z = (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z :=
    fun z hz h => u1_O hb hz (by rw [← hY]; exact mul_le_mul_of_nonneg_right h hd.le)
  have pT : ∀ z, bandY1 H ≤ z → bandU1 H z = (1 - (2 * Real.sqrt 2 - H)) + bandT H z :=
    fun z h => u1_aT hb (by rw [← hY]; exact mul_le_mul_of_nonneg_right h hd.le)
  rcases le_total y (bandY1 H) with hyY | hyY
  · have := kink_lo (A := M - (1 - (2 * Real.sqrt 2 - H))) (B := (C - H / 2) - (1 - H / 2)) (c := S) (s := C)
      (lo := 0) (hi := bandY1 H) (m := 0) hS hy hyY
      (by rw [pO 0 le_rfl hY0] at h0; unfold kh at h0; linarith)
      (by rw [pO _ hY0 le_rfl] at h1; unfold kh at h1; linarith)
      (fun t h1t h2t e => by
        have := ht t h1t e; rw [pO t h1t h2t] at this; unfold kh at this
        rw [e, sub_self, max_self] at this; linarith)
    rw [pO y hy hyY]; unfold kh; linarith
  · have := kink_lo_inf (A := M - (1 - (2 * Real.sqrt 2 - H)) - bandT H 0) (B := (C - H / 2) - (H / 16 - 1 / 10))
      (c := S) (s := C) (lo := bandY1 H) (m := 0) hS hyY (by linarith)
      (by rw [pT _ le_rfl] at h1; unfold kh bandT at h1; unfold bandT; linarith)
      (fun t h1t e => by
        have := ht t (hY0.trans h1t) e; rw [pT t h1t] at this; unfold kh bandT at this
        rw [e, sub_self, max_self] at this; unfold bandT; linarith)
    rw [pT y hyY]; unfold kh bandT at *; linarith

/-- **Upper bound along one line, even potential.**  `u₀ - kh` on `[0, ∞)` (`A0` up to `x₂`, `Lz` up to `8/7`,
then `T`): it suffices to check `0`, `x₂`, `8/7` and the aligned slope. -/
lemma line_u0_hi (hb : InBand H) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) {M : ℝ}
    (h0 : bandU0 H 0 - kh S C H 0 ≤ M) (h1 : bandU0 H (bandX2 H) - kh S C H (bandX2 H) ≤ M)
    (h2 : bandU0 H (8 / 7) - kh S C H (8 / 7) ≤ M)
    (ht : ∀ t, 0 ≤ t → S * t = C → bandU0 H t - kh S C H t ≤ M)
    {y : ℝ} (hy : 0 ≤ y) : bandU0 H y - kh S C H y ≤ M := by
  obtain ⟨h9, h87⟩ := hb
  have hb := (⟨h9, h87⟩ : InBand H)
  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  have hd : 0 < 2 * (1 - H / 2) := by linarith
  have hX := bandX2_eq hb
  have hX0 : 0 ≤ bandX2 H := by
    have : 0 ≤ bandX2 H * (2 * (1 - H / 2)) := by rw [hX]; linarith
    exact nonneg_of_mul_nonneg_left this hd
  have hX8 : bandX2 H ≤ 8 / 7 := by
    have : bandX2 H * (2 * (1 - H / 2)) ≤ 8 / 7 * (2 * (1 - H / 2)) := by rw [hX]; linarith
    exact le_of_mul_le_mul_right this hd
  have pA : ∀ z, 0 ≤ z → z ≤ bandX2 H → bandU0 H z = 2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z :=
    fun z hz h => u0_A0 hb hz (by rw [← hX]; exact mul_le_mul_of_nonneg_right h hd.le)
  have pL : ∀ z, bandX2 H ≤ z → z ≤ 8 / 7 → bandU0 H z = bandLz H z :=
    fun z h h' => u0_Lz hb (hX0.trans h) (by rw [← hX]; exact mul_le_mul_of_nonneg_right h hd.le) (by linarith)
  have pT : ∀ z, 8 / 7 ≤ z → bandU0 H z = bandT H z := fun z h => u0_T hb (by linarith)
  rcases le_total y (bandX2 H) with hy1 | hy1
  · have := kink_lo (A := M - 2 * (1 - (2 * Real.sqrt 2 - H))) (B := (C - H / 2) - (1 - H / 2)) (c := S) (s := C)
      (lo := 0) (hi := bandX2 H) (m := 0) hS hy hy1
      (by rw [pA 0 le_rfl hX0] at h0; unfold kh at h0; linarith)
      (by rw [pA _ hX0 le_rfl] at h1; unfold kh at h1; linarith)
      (fun t h1t h2t e => by
        have := ht t h1t e; rw [pA t h1t h2t] at this; unfold kh at this
        rw [e, sub_self, max_self] at this; linarith)
    rw [pA y hy hy1]; unfold kh; linarith
  rcases le_total y (8 / 7) with hy2 | hy2
  · have := kink_lo (A := M - bandLz H 0) (B := (C - H / 2) + (1 - H / 2)) (c := S) (s := C)
      (lo := bandX2 H) (hi := 8 / 7) (m := 0) hS hy1 hy2
      (by rw [pL _ le_rfl hX8] at h1; unfold kh bandLz at h1; unfold bandLz; linarith)
      (by rw [pL _ hX8 le_rfl] at h2; unfold kh bandLz at h2; unfold bandLz; linarith)
      (fun t h1t h2t e => by
        have := ht t (hX0.trans h1t) e; rw [pL t h1t h2t] at this; unfold kh bandLz at this
        rw [e, sub_self, max_self] at this; unfold bandLz; linarith)
    rw [pL y hy1 hy2]; unfold kh bandLz at *; linarith
  · have := kink_lo_inf (A := M - bandT H 0) (B := (C - H / 2) - (H / 16 - 1 / 10)) (c := S) (s := C)
      (lo := 8 / 7) (m := 0) hS hy2 (by linarith)
      (by rw [pT _ le_rfl] at h2; unfold kh bandT at h2; unfold bandT; linarith)
      (fun t h1t e => by
        have := ht t (by linarith) e; rw [pT t h1t] at this; unfold kh bandT at this
        rw [e, sub_self, max_self] at this; unfold bandT; linarith)
    rw [pT y hy2]; unfold kh bandT at *; linarith

end lines

end SquarePacking
