import Sqpack.Rect4

/-!
# Five unit squares in a strip: the middle square between two lines

For `R₅` with short side `< 2`, the chain `a, b, c, d, e` has a middle square `c` between two separating lines
(`b | c` and `c | d`).  Its width between them, `s₂/S₂ - s₁/S₁ ≥ μ y + (|cr a₁ e| + a₁·e)/(2S₁) +
(|cr a₂ e| + a₂·e)/(2S₂)` (`aᵢ` the frames aligned with the lines, `e` the frame of `c`, `y` its height),
is, after `y` goes to the floor or the ceiling, linear in `e` on each arc of the quarter circle between the
kinks `a₁`, `a₂`.  So its minimum is at `c` upright or aligned with one of the lines (`mid_set`).

* `arc_lin` — a linear form at least `m ≥ 0` at both ends of an arc of the quarter circle is at least `m` on it
  (the arc lies outside its chord; degenerate arcs by the frame decomposition `decomp`).
* `cr_trans` — the angle order (`cr x y ≥ 0`) is transitive on the quarter circle.
* `kink_min` — linear plus two kinks: at least the value at `(1, 0)`, `a₁`, `a₂` or `(0, 1)`.
* `mid_floor`, `mid_ceil`, `mid_set` — the middle square.
-/

namespace SquarePacking

/-- Cross product `x₁ y₂ - x₂ y₁` (`≥ 0` iff `y` is anticlockwise of `x`, on the quarter circle). -/
def cr (x y : ℝ × ℝ) : ℝ := x.1 * y.2 - x.2 * y.1

/-- Dot product. -/
def dt (x y : ℝ × ℝ) : ℝ := x.1 * y.1 + x.2 * y.2

/-- A unit vector in the closed first quadrant. -/
def QV (x : ℝ × ℝ) : Prop := 0 ≤ x.1 ∧ 0 ≤ x.2 ∧ x.1 ^ 2 + x.2 ^ 2 = 1

lemma QV.dt_nonneg {x y : ℝ × ℝ} (hx : QV x) (hy : QV y) : 0 ≤ dt x y := by
  obtain ⟨a, b, -⟩ := hx; obtain ⟨c, d, -⟩ := hy; unfold dt; positivity

lemma QV.dt_le_one {x y : ℝ × ℝ} (hx : QV x) (hy : QV y) : dt x y ≤ 1 := by
  obtain ⟨-, -, h1⟩ := hx; obtain ⟨-, -, h2⟩ := hy; unfold dt
  nlinarith [sq_nonneg (x.1 - y.1), sq_nonneg (x.2 - y.2)]

lemma QV.pyth {x y : ℝ × ℝ} (hx : QV x) (hy : QV y) : dt x y ^ 2 + cr x y ^ 2 = 1 := by
  obtain ⟨-, -, h1⟩ := hx; obtain ⟨-, -, h2⟩ := hy; unfold dt cr
  linear_combination (x.1 ^ 2 + x.2 ^ 2) * h2 + h1

lemma eq_one_of_sq {x : ℝ} (h0 : 0 ≤ x) (h : x ^ 2 = 1) : x = 1 := by
  have : (x - 1) * (x + 1) = 0 := by linear_combination h
  rcases mul_eq_zero.1 this with h1 | h1 <;> linarith

/-- Coordinates of `x` in the frame `(e, e⊥)` of a unit vector `e`. -/
lemma decomp {x e : ℝ × ℝ} (he : e.1 ^ 2 + e.2 ^ 2 = 1) :
    x.1 = dt x e * e.1 + cr x e * e.2 ∧ x.2 = dt x e * e.2 - cr x e * e.1 := by
  unfold dt cr; constructor
  · linear_combination (-x.1) * he
  · linear_combination (-x.2) * he

/-- The angle order is transitive on the quarter circle. -/
lemma cr_trans {x y z : ℝ × ℝ} (hx : QV x) (hy : QV y) (hz : QV z) (h1 : 0 ≤ cr x y) (h2 : 0 ≤ cr y z) :
    0 ≤ cr x z := by
  have e : cr x z = dt x y * cr y z + cr x y * dt y z := by
    obtain ⟨-, -, hy1⟩ := hy; unfold cr dt; linear_combination (-(x.1 * z.2 - x.2 * z.1)) * hy1
  rw [e]; have := hx.dt_nonneg hy; have := hy.dt_nonneg hz; positivity

/-- **Arc lemma.**  `e` between `e₁` and `e₂` on the quarter circle; a linear form at least `m ≥ 0` at both
ends is at least `m` at `e`. -/
theorem arc_lin {e₁ e₂ e : ℝ × ℝ} (h₁ : QV e₁) (h₂ : QV e₂) (he : QV e) (o₁ : 0 ≤ cr e₁ e) (o₂ : 0 ≤ cr e e₂)
    (α β m : ℝ) (hm : 0 ≤ m) (l₁ : m ≤ α * e₁.1 + β * e₁.2) (l₂ : m ≤ α * e₂.1 + β * e₂.2) :
    m ≤ α * e.1 + β * e.2 := by
  set s₁ := cr e₁ e; set s₂ := cr e e₂; set c₁ := dt e₁ e; set c₂ := dt e e₂
  have hc₁ : 0 ≤ c₁ := h₁.dt_nonneg he
  have hc₂ : 0 ≤ c₂ := he.dt_nonneg h₂
  -- `D e = s₂ e₁ + s₁ e₂` with `D = s₂ c₁ + s₁ c₂`
  have he1 : (s₂ * c₁ + s₁ * c₂) * e.1 = s₂ * e₁.1 + s₁ * e₂.1 := by
    obtain ⟨-, -, h⟩ := he; simp only [s₁, s₂, c₁, c₂, cr, dt]; linear_combination e.1 * (e₁.1 * e₂.2 - e₁.2 * e₂.1) * h
  have he2 : (s₂ * c₁ + s₁ * c₂) * e.2 = s₂ * e₁.2 + s₁ * e₂.2 := by
    obtain ⟨-, -, h⟩ := he; simp only [s₁, s₂, c₁, c₂, cr, dt]; linear_combination e.2 * (e₁.1 * e₂.2 - e₁.2 * e₂.1) * h
  have he' := he
  have hD : 0 ≤ s₂ * c₁ + s₁ * c₂ := by positivity
  rcases eq_or_lt_of_le hD with hD0 | hDp
  · -- degenerate: `e` is an end of the arc
    have p₁ : s₂ * c₁ = 0 := by nlinarith [mul_nonneg o₂ hc₁, mul_nonneg o₁ hc₂]
    have p₂ : s₁ * c₂ = 0 := by nlinarith [mul_nonneg o₂ hc₁, mul_nonneg o₁ hc₂]
    have y₁ := h₁.pyth he; have y₂ := he.pyth h₂
    obtain ⟨-, -, hee⟩ := he
    obtain ⟨x1, x2⟩ := decomp (x := e₁) hee
    change e₁.1 = c₁ * e.1 + s₁ * e.2 at x1
    change e₁.2 = c₁ * e.2 - s₁ * e.1 at x2
    obtain ⟨z1, z2⟩ := decomp (x := e₂) hee
    have cz : cr e₂ e = -s₂ := by simp only [s₂, cr]; ring
    have dz : dt e₂ e = c₂ := by simp only [c₂, dt]; ring
    rw [cz, dz] at z1 z2
    change c₁ ^ 2 + s₁ ^ 2 = 1 at y₁; change c₂ ^ 2 + s₂ ^ 2 = 1 at y₂
    rcases eq_or_lt_of_le hc₁ with hc | hc
    · -- `c₁ = 0`: `s₁ = 1`, `c₂ = 0`, `s₂ = 1`, so `e₁ = -e₂`: impossible on the quarter circle
      exfalso
      have hs₁ : s₁ = 1 := eq_one_of_sq o₁ (by rw [← hc] at y₁; linarith)
      have hc₂0 : c₂ = 0 := by rw [hs₁, one_mul] at p₂; exact p₂
      have hs₂ : s₂ = 1 := eq_one_of_sq o₂ (by rw [hc₂0] at y₂; linarith)
      rw [← hc, hs₁] at x1 x2; rw [hc₂0, hs₂] at z1 z2
      obtain ⟨a1, a2, -⟩ := h₁; obtain ⟨b1, b2, -⟩ := h₂
      have q1 : e.1 = 0 := by linarith
      have q2 : e.2 = 0 := by linarith
      rw [q1, q2] at hee; norm_num at hee
    · -- `c₁ > 0`: `s₂ = 0`, so `c₂ = 1`, `s₁ = 0`, `c₁ = 1`: `e = e₁`
      have hs₂ : s₂ = 0 := by
        rcases mul_eq_zero.1 p₁ with h | h
        · exact h
        · linarith
      have hc₂1 : c₂ = 1 := eq_one_of_sq hc₂ (by rw [hs₂] at y₂; linarith)
      have hs₁ : s₁ = 0 := by rw [hc₂1, mul_one] at p₂; exact p₂
      have hc₁1 : c₁ = 1 := eq_one_of_sq hc₁ (by rw [hs₁] at y₁; linarith)
      rw [hc₁1, hs₁] at x1 x2
      have d1 : e.1 = e₁.1 := by linarith
      have d2 : e.2 = e₁.2 := by linarith
      rw [d1, d2]; exact l₁
  · -- `s₁ + s₂ ≥ D`: the arc lies outside its chord
    have chord : s₂ * c₁ + s₁ * c₂ ≤ s₁ + s₂ := by
      nlinarith [mul_le_mul_of_nonneg_left (h₁.dt_le_one he) o₂, mul_le_mul_of_nonneg_left (he.dt_le_one h₂) o₁]
    have key : (s₂ * c₁ + s₁ * c₂) * (α * e.1 + β * e.2) = s₂ * (α * e₁.1 + β * e₁.2) + s₁ * (α * e₂.1 + β * e₂.2) := by
      linear_combination α * he1 + β * he2
    have : (s₂ * c₁ + s₁ * c₂) * m ≤ (s₂ * c₁ + s₁ * c₂) * (α * e.1 + β * e.2) := by
      rw [key]; nlinarith [mul_le_mul_of_nonneg_left l₁ o₂, mul_le_mul_of_nonneg_left l₂ o₁]
    exact le_of_mul_le_mul_left this hDp

/-- A function on the quarter circle that is linear plus two kinks `|cr aᵢ x|`. -/
def kinkG (α β w₁ w₂ : ℝ) (a₁ a₂ x : ℝ × ℝ) : ℝ :=
  α * x.1 + β * x.2 + w₁ * |cr a₁ x| + w₂ * |cr a₂ x|

/-- On an arc where both kinks keep their signs (`|cr aᵢ x| = tᵢ cr aᵢ x` at the ends and at `e`),
`kinkG` is linear, so it is at least its smaller end value (when that is `≥ 0`). -/
lemma arc_abs {a₁ a₂ e₁ e₂ e : ℝ × ℝ} (h₁ : QV e₁) (h₂ : QV e₂) (he : QV e) (o₁ : 0 ≤ cr e₁ e)
    (o₂ : 0 ≤ cr e e₂) (α β w₁ w₂ t₁ t₂ : ℝ)
    (s : ∀ x, x = e₁ ∨ x = e₂ ∨ x = e → |cr a₁ x| = t₁ * cr a₁ x ∧ |cr a₂ x| = t₂ * cr a₂ x)
    (hm : 0 ≤ min (kinkG α β w₁ w₂ a₁ a₂ e₁) (kinkG α β w₁ w₂ a₁ a₂ e₂)) :
    min (kinkG α β w₁ w₂ a₁ a₂ e₁) (kinkG α β w₁ w₂ a₁ a₂ e₂) ≤ kinkG α β w₁ w₂ a₁ a₂ e := by
  have lin : ∀ x, x = e₁ ∨ x = e₂ ∨ x = e → kinkG α β w₁ w₂ a₁ a₂ x =
      (α - w₁ * t₁ * a₁.2 - w₂ * t₂ * a₂.2) * x.1 + (β + w₁ * t₁ * a₁.1 + w₂ * t₂ * a₂.1) * x.2 := by
    intro x hx; obtain ⟨u₁, u₂⟩ := s x hx
    unfold kinkG; rw [u₁, u₂]; unfold cr; ring
  rw [lin e (Or.inr (Or.inr rfl))]
  apply arc_lin h₁ h₂ he o₁ o₂ _ _ _ hm
  · rw [← lin e₁ (Or.inl rfl)]; exact min_le_left _ _
  · rw [← lin e₂ (Or.inr (Or.inl rfl))]; exact min_le_right _ _

lemma QV.e0 : QV ((1, 0) : ℝ × ℝ) := ⟨by norm_num, by norm_num, by norm_num⟩
lemma QV.e90 : QV ((0, 1) : ℝ × ℝ) := ⟨by norm_num, by norm_num, by norm_num⟩

/-- **Kinks.**  `kinkG` (nonnegative on the quarter circle) is at least its least value at `(1, 0)`, `a₁`,
`a₂`, `(0, 1)`; here with `a₁` before `a₂`. -/
lemma kink_min_ord {a₁ a₂ e : ℝ × ℝ} (ha₁ : QV a₁) (ha₂ : QV a₂) (he : QV e) (ord : 0 ≤ cr a₁ a₂)
    (α β w₁ w₂ : ℝ) (hpos : ∀ x, QV x → 0 ≤ kinkG α β w₁ w₂ a₁ a₂ x) :
    min (min (kinkG α β w₁ w₂ a₁ a₂ (1, 0)) (kinkG α β w₁ w₂ a₁ a₂ a₁))
      (min (kinkG α β w₁ w₂ a₁ a₂ a₂) (kinkG α β w₁ w₂ a₁ a₂ (0, 1))) ≤ kinkG α β w₁ w₂ a₁ a₂ e := by
  have anti : ∀ x y : ℝ × ℝ, cr y x = -cr x y := fun x y => by unfold cr; ring
  have c0 : ∀ a : ℝ × ℝ, QV a → cr a (1, 0) ≤ 0 := fun a ha => by unfold cr; simp; exact ha.2.1
  have c90 : ∀ a : ℝ × ℝ, QV a → 0 ≤ cr a (0, 1) := fun a ha => by unfold cr; simp; exact ha.1
  have cself : ∀ a : ℝ × ℝ, cr a a = 0 := fun a => by unfold cr; ring
  have hp := fun x hx => hpos x hx
  by_cases i : 0 ≤ cr e a₁
  · -- `e` before `a₁`: arc `(1, 0) → a₁`, both kinks negative
    have e2 : 0 ≤ cr e a₂ := cr_trans he ha₁ ha₂ i ord
    refine le_trans (min_le_left _ _) (arc_abs QV.e0 ha₁ he (by unfold cr; simp; exact he.2.1) i α β w₁ w₂
      (-1) (-1) ?_ (le_min (hp _ QV.e0) (hp _ ha₁)))
    rintro x (rfl | rfl | rfl)
    · exact ⟨by rw [abs_of_nonpos (c0 _ ha₁)]; ring, by rw [abs_of_nonpos (c0 _ ha₂)]; ring⟩
    · exact ⟨by rw [cself, abs_zero]; ring, by rw [abs_of_nonpos (by rw [anti]; linarith)]; ring⟩
    · exact ⟨by rw [abs_of_nonpos (by rw [anti]; linarith)]; ring,
        by rw [abs_of_nonpos (by rw [anti]; linarith)]; ring⟩
  push Not at i
  have i' : 0 ≤ cr a₁ e := by rw [anti]; linarith
  by_cases j : 0 ≤ cr a₂ e
  · -- `e` after `a₂`: arc `a₂ → (0, 1)`, both kinks positive
    have e1 : 0 ≤ cr a₁ e := i'
    refine le_trans (min_le_right _ _) (arc_abs ha₂ QV.e90 he j (by unfold cr; simp; exact he.1) α β w₁ w₂
      1 1 ?_ (le_min (hp _ ha₂) (hp _ QV.e90)))
    rintro x (rfl | rfl | rfl)
    · exact ⟨by rw [abs_of_nonneg ord]; ring, by rw [cself, abs_zero]; ring⟩
    · exact ⟨by rw [abs_of_nonneg (c90 _ ha₁)]; ring, by rw [abs_of_nonneg (c90 _ ha₂)]; ring⟩
    · exact ⟨by rw [abs_of_nonneg e1]; ring, by rw [abs_of_nonneg j]; ring⟩
  push Not at j
  -- between: arc `a₁ → a₂`, first kink positive, second negative
  have j' : 0 ≤ cr e a₂ := by rw [anti]; linarith
  refine le_trans (le_min (le_trans (min_le_left _ _) (min_le_right _ _))
    (le_trans (min_le_right _ _) (min_le_left _ _)))
    (arc_abs ha₁ ha₂ he i' j' α β w₁ w₂ 1 (-1) ?_ (le_min (hp _ ha₁) (hp _ ha₂)))
  rintro x (rfl | rfl | rfl)
  · exact ⟨by rw [cself, abs_zero]; ring, by rw [abs_of_nonpos (by rw [anti]; linarith)]; ring⟩
  · exact ⟨by rw [abs_of_nonneg ord]; ring, by rw [cself, abs_zero]; ring⟩
  · exact ⟨by rw [abs_of_nonneg i']; ring, by rw [abs_of_nonpos j.le]; ring⟩

/-- **Kinks**, any order: `kinkG` is at least its value at one of `(1, 0)`, `a₁`, `a₂`, `(0, 1)`. -/
theorem kink_min {a₁ a₂ e : ℝ × ℝ} (ha₁ : QV a₁) (ha₂ : QV a₂) (he : QV e) (α β w₁ w₂ : ℝ)
    (hpos : ∀ x, QV x → 0 ≤ kinkG α β w₁ w₂ a₁ a₂ x) :
    kinkG α β w₁ w₂ a₁ a₂ (1, 0) ≤ kinkG α β w₁ w₂ a₁ a₂ e ∨ kinkG α β w₁ w₂ a₁ a₂ a₁ ≤ kinkG α β w₁ w₂ a₁ a₂ e ∨
      kinkG α β w₁ w₂ a₁ a₂ a₂ ≤ kinkG α β w₁ w₂ a₁ a₂ e ∨
      kinkG α β w₁ w₂ a₁ a₂ (0, 1) ≤ kinkG α β w₁ w₂ a₁ a₂ e := by
  rcases le_total 0 (cr a₁ a₂) with ord | ord
  · have := kink_min_ord ha₁ ha₂ he ord α β w₁ w₂ hpos
    simp only [min_le_iff] at this; tauto
  · have sw : ∀ x, kinkG α β w₁ w₂ a₁ a₂ x = kinkG α β w₂ w₁ a₂ a₁ x := fun x => by unfold kinkG; ring
    have ord' : 0 ≤ cr a₂ a₁ := by unfold cr at ord ⊢; linarith
    have := kink_min_ord ha₂ ha₁ he ord' α β w₂ w₁ (fun x hx => by rw [← sw]; exact hpos x hx)
    simp only [min_le_iff, ← sw] at this; tauto

/-- The width of the middle square between two lines, as a function of its frame `e` and centre height `y`:
`μ y + c₁ (|cr a₁ e| + a₁ · e) + c₂ (|cr a₂ e| + a₂ · e)` (`aᵢ` the frames aligned with the lines). -/
def midW (μ c₁ c₂ : ℝ) (a₁ a₂ e : ℝ × ℝ) (y : ℝ) : ℝ :=
  μ * y + c₁ * (|cr a₁ e| + dt a₁ e) + c₂ * (|cr a₂ e| + dt a₂ e)

/-- **Middle square, `μ ≥ 0`** (it sits lowest on the floor): the width is at least the width of a square on the
floor that is upright or aligned with one of the lines. -/
theorem mid_floor {a₁ a₂ e : ℝ × ℝ} (ha₁ : QV a₁) (ha₂ : QV a₂) (he : QV e) {μ c₁ c₂ y : ℝ} (hμ : 0 ≤ μ)
    (hc₁ : 0 ≤ c₁) (hc₂ : 0 ≤ c₂) (hy : (e.1 + e.2) / 2 ≤ y) :
    ∃ x, (x = (1, 0) ∨ x = a₁ ∨ x = a₂ ∨ x = (0, 1)) ∧
      midW μ c₁ c₂ a₁ a₂ x ((x.1 + x.2) / 2) ≤ midW μ c₁ c₂ a₁ a₂ e y := by
  have form : ∀ x : ℝ × ℝ, midW μ c₁ c₂ a₁ a₂ x ((x.1 + x.2) / 2) =
      kinkG (μ / 2 + c₁ * a₁.1 + c₂ * a₂.1) (μ / 2 + c₁ * a₁.2 + c₂ * a₂.2) c₁ c₂ a₁ a₂ x := by
    intro x; unfold midW kinkG dt; ring
  have hpos : ∀ x, QV x → 0 ≤ kinkG (μ / 2 + c₁ * a₁.1 + c₂ * a₂.1) (μ / 2 + c₁ * a₁.2 + c₂ * a₂.2)
      c₁ c₂ a₁ a₂ x := by
    intro x hx; rw [← form]; unfold midW
    have := ha₁.dt_nonneg hx; have := ha₂.dt_nonneg hx; obtain ⟨h1, h2, -⟩ := hx; positivity
  have hW : midW μ c₁ c₂ a₁ a₂ e ((e.1 + e.2) / 2) ≤ midW μ c₁ c₂ a₁ a₂ e y := by
    unfold midW; nlinarith [mul_le_mul_of_nonneg_left hy hμ]
  rcases kink_min ha₁ ha₂ he _ _ c₁ c₂ hpos with h | h | h | h <;> rw [← form, ← form] at h
  · exact ⟨_, Or.inl rfl, h.trans hW⟩
  · exact ⟨_, Or.inr (Or.inl rfl), h.trans hW⟩
  · exact ⟨_, Or.inr (Or.inr (Or.inl rfl)), h.trans hW⟩
  · exact ⟨_, Or.inr (Or.inr (Or.inr rfl)), h.trans hW⟩

/-- **Middle square, `μ ≤ 0`** (it sits highest against the ceiling). -/
theorem mid_ceil {a₁ a₂ e : ℝ × ℝ} (ha₁ : QV a₁) (ha₂ : QV a₂) (he : QV e) {μ c₁ c₂ y H : ℝ} (hμ : μ ≤ 0)
    (hc₁ : 0 ≤ c₁) (hc₂ : 0 ≤ c₂) (hy : y ≤ H - (e.1 + e.2) / 2) :
    ∃ x, (x = (1, 0) ∨ x = a₁ ∨ x = a₂ ∨ x = (0, 1)) ∧
      midW μ c₁ c₂ a₁ a₂ x (H - (x.1 + x.2) / 2) ≤ midW μ c₁ c₂ a₁ a₂ e y := by
  have form : ∀ x : ℝ × ℝ, midW μ c₁ c₂ a₁ a₂ x (H - (x.1 + x.2) / 2) = μ * H +
      kinkG (-μ / 2 + c₁ * a₁.1 + c₂ * a₂.1) (-μ / 2 + c₁ * a₁.2 + c₂ * a₂.2) c₁ c₂ a₁ a₂ x := by
    intro x; unfold midW kinkG dt; ring
  have hpos : ∀ x, QV x → 0 ≤ kinkG (-μ / 2 + c₁ * a₁.1 + c₂ * a₂.1) (-μ / 2 + c₁ * a₁.2 + c₂ * a₂.2)
      c₁ c₂ a₁ a₂ x := by
    intro x hx
    have e := form x; unfold midW at e
    have := ha₁.dt_nonneg hx; have := ha₂.dt_nonneg hx; obtain ⟨h1, h2, -⟩ := hx
    have : 0 ≤ -μ * ((x.1 + x.2) / 2) + c₁ * (|cr a₁ x| + dt a₁ x) + c₂ * (|cr a₂ x| + dt a₂ x) := by
      have : 0 ≤ -μ := by linarith
      positivity
    linarith
  have hW : midW μ c₁ c₂ a₁ a₂ e (H - (e.1 + e.2) / 2) ≤ midW μ c₁ c₂ a₁ a₂ e y := by
    unfold midW; nlinarith [mul_le_mul_of_nonpos_left hy hμ]
  rcases kink_min ha₁ ha₂ he _ _ c₁ c₂ hpos with h | h | h | h
  · exact ⟨_, Or.inl rfl, by rw [form]; linarith [form e]⟩
  · exact ⟨_, Or.inr (Or.inl rfl), by rw [form]; linarith [form e]⟩
  · exact ⟨_, Or.inr (Or.inr (Or.inl rfl)), by rw [form]; linarith [form e]⟩
  · exact ⟨_, Or.inr (Or.inr (Or.inr rfl)), by rw [form]; linarith [form e]⟩

/-! ## Lines and their aligned frames -/

/-- A line with normal `(S, v)`, `S > 0`, has aligned frame `a`: the support of a square with frame `e` in the
direction `(S, v)` is `|cr a e| + a · e`. -/
def Aligned (S v : ℝ) (a : ℝ × ℝ) : Prop :=
  QV a ∧ ∀ P Q : ℝ, 0 ≤ P → 0 ≤ Q → |S * P + v * Q| + |-S * Q + v * P| = |cr a (P, Q)| + dt a (P, Q)

/-- Normal `(S, -C)`: the aligned frame is `(C, S)`. -/
lemma aligned_neg {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) : Aligned S (-C) (C, S) := by
  refine ⟨⟨hC, hS, hCS⟩, fun P Q hP hQ => ?_⟩
  unfold cr dt; simp only
  rw [show S * P + -C * Q = -(C * Q - S * P) by ring, abs_neg,
    show -S * Q + -C * P = -(C * P + S * Q) by ring, abs_neg, abs_of_nonneg (by positivity : 0 ≤ C * P + S * Q)]

/-- Normal `(S, C)`: the aligned frame is `(S, C)`. -/
lemma aligned_pos {C S : ℝ} (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) : Aligned S C (S, C) := by
  refine ⟨⟨hS, hC, by linarith⟩, fun P Q hP hQ => ?_⟩
  unfold cr dt; simp only
  rw [abs_of_nonneg (by positivity : 0 ≤ S * P + C * Q),
    show -S * Q + C * P = -(S * Q - C * P) by ring, abs_neg]; ring

/-- **The middle square between two lines** (set form).  `c` lies right of line 1 (`S₁ x + v₁ y ≥ s₁`) and
left of line 2 (`S₂ x + v₂ y ≤ s₂`).  Then `s₂/S₂ - s₁/S₁ ≥ midW μ (1/(2S₁)) (1/(2S₂)) a₁ a₂ x y*` for some `x`
among `(1, 0)`, `a₁`, `a₂`, `(0, 1)`, with `μ = v₂/S₂ - v₁/S₁` and `y*` the floor (`μ ≥ 0`) or ceiling
(`μ ≤ 0`) height of a square with frame `x`. -/
theorem mid_set {c : ℝ × ℝ} {θ L H S₁ v₁ s₁ S₂ v₂ s₂ : ℝ} {a₁ a₂ : ℝ × ℝ} (hin : sq c θ 1 ⊆ rbox L H)
    (hS₁ : 0 < S₁) (hS₂ : 0 < S₂) (A₁ : Aligned S₁ v₁ a₁) (A₂ : Aligned S₂ v₂ a₂)
    (h₁ : ∀ p ∈ sq c θ 1, s₁ ≤ S₁ * p.1 + v₁ * p.2) (h₂ : ∀ p ∈ sq c θ 1, S₂ * p.1 + v₂ * p.2 ≤ s₂) :
    ∃ x, (x = (1, 0) ∨ x = a₁ ∨ x = a₂ ∨ x = (0, 1)) ∧
      midW (v₂ / S₂ - v₁ / S₁) (1 / (2 * S₁)) (1 / (2 * S₂)) a₁ a₂ x
        (if 0 ≤ v₂ / S₂ - v₁ / S₁ then (x.1 + x.2) / 2 else H - (x.1 + x.2) / 2) ≤ s₂ / S₂ - s₁ / S₁ := by
  obtain ⟨P, Q, hP, hQ, hPQ, -, -, hy1, hy2, hpt⟩ := framed hin
  have he : QV (P, Q) := ⟨hP, hQ, hPQ⟩
  have g₁ := support_ge hpt h₁
  have g₂ := support_le hpt h₂
  rw [A₁.2 P Q hP hQ] at g₁; rw [A₂.2 P Q hP hQ] at g₂
  have hW : midW (v₂ / S₂ - v₁ / S₁) (1 / (2 * S₁)) (1 / (2 * S₂)) a₁ a₂ (P, Q) c.2 ≤ s₂ / S₂ - s₁ / S₁ := by
    have b₂ : (S₂ * c.1 + v₂ * c.2 + (|cr a₂ (P, Q)| + dt a₂ (P, Q)) / 2) / S₂ ≤ s₂ / S₂ :=
      div_le_div_of_nonneg_right g₂ hS₂.le
    have b₁ : s₁ / S₁ ≤ (S₁ * c.1 + v₁ * c.2 - (|cr a₁ (P, Q)| + dt a₁ (P, Q)) / 2) / S₁ :=
      div_le_div_of_nonneg_right g₁ hS₁.le
    have e : midW (v₂ / S₂ - v₁ / S₁) (1 / (2 * S₁)) (1 / (2 * S₂)) a₁ a₂ (P, Q) c.2 =
        (S₂ * c.1 + v₂ * c.2 + (|cr a₂ (P, Q)| + dt a₂ (P, Q)) / 2) / S₂ -
          (S₁ * c.1 + v₁ * c.2 - (|cr a₁ (P, Q)| + dt a₁ (P, Q)) / 2) / S₁ := by
      unfold midW; field_simp; ring
    rw [e]; linarith
  have hc₁ : 0 ≤ 1 / (2 * S₁) := by positivity
  have hc₂ : 0 ≤ 1 / (2 * S₂) := by positivity
  split_ifs with hμ
  · obtain ⟨x, hx, hle⟩ := mid_floor A₁.1 A₂.1 he hμ hc₁ hc₂ (y := c.2) (by simpa using hy1)
    exact ⟨x, hx, hle.trans hW⟩
  · obtain ⟨x, hx, hle⟩ := mid_ceil A₁.1 A₂.1 he (le_of_lt (not_le.1 hμ)) hc₁ hc₂ (y := c.2) (H := H)
      (by simpa using hy2)
    exact ⟨x, hx, hle.trans hW⟩

end SquarePacking
