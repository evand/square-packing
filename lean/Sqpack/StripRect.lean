import Sqpack.ChordLine
import Sqpack.D4

/-!
# Rectangle and slab corollaries of the wall-strip lemma, parameterised by the strip height

`notes/strip-split.md` §1–3.  Everything here is stated for an abstract strip height `h` through
the predicate `WallStripBound h` (the bottom-wall statement of `wall_strip_le_three_line` with
`(3√2 - 2)/2` replaced by `h`), so that a stronger strip lemma upgrades every corollary.
Available in Lean: `wallStripBound_one` (`h = 1`, `Chord.lean`) and `wallStripBound_line`
(`h = (3√2 - 2)/2 ≈ 1.1213`, `ChordLine.lean`).

Semantics throughout: closed unit squares `sq c θ 1`, container `box t` with `t ≤ 4`, pairwise
disjoint as closed sets.  `vslab a b t = [a,b] × [0,t]`.

* `WallStripBound.card_le_three` — the Finset-indexed bottom-wall form.
* `slab_left_card_le_three`, `slab_right_card_le_three` — inside a slab `[a,b] × [0,t] ⊆ box t`,
  at most three squares have centre `x ≤ a + h`, and at most three have centre `x ≥ b - h`.
  Special cases: `left_wall_le_three` (box, `x ≤ h`) and `right_wall_le_three`
  (`[0,w] × [0,t]`, `x ≥ w - h`).
* `sq_subset_vslab_shift`, `sq_shift_eq` — translating a slab to `[0,w] × [0,t]`.
* Slab capacities (R1, R4, R7): `slab_R1`, `slab_R4`, `slab_R7` (and Finset forms
  `slab_one_le_width`, `slab_width_of_four`, `slab_width_of_seven`).
* `guillotine_false` — no `(k, 4, 8-k)` guillotine split of 12 squares when `h ≥ 7/6`.
* `line_sweep` — for `u ≤ h + 1/2` and `t - u ≤ 2h`, at least 3 squares meet the line `x = u`
  (are contained in neither side).
-/

open Finset
open scoped Classical

namespace SquarePacking

/-! ## 0. The strip predicate -/

/-- **The wall-strip bound at height `h`.**  In `[0,t]^2`, `t ≤ 4`, at most three closed unit
squares, pairwise disjoint as closed sets, have their centre at height at most `h`. -/
def WallStripBound (h : ℝ) : Prop :=
  ∀ (N : ℕ) (t : ℝ), t ≤ 4 → ∀ (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ),
    (∀ i, sq (ctr i) (ang i) 1 ⊆ box t) → (∀ i, (ctr i).2 ≤ h) →
    (∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) → N ≤ 3

theorem wallStripBound_line : WallStripBound ((3 * Real.sqrt 2 - 2) / 2) :=
  fun _ _ ht ctr ang hin hs hd => wall_strip_le_three_line ht ctr ang hin hs hd

theorem wallStripBound_one : WallStripBound 1 :=
  fun _ _ ht ctr ang hin hs hd => wall_strip_le_three ht ctr ang hin hs hd

theorem WallStripBound.mono {h h' : ℝ} (hh : h' ≤ h) (hW : WallStripBound h) :
    WallStripBound h' :=
  fun N t ht ctr ang hin hs hd => hW N t ht ctr ang hin (fun i => (hs i).trans hh) hd

/-! ## 1. The Finset-indexed form, and transport along a point map -/

/-- Transport: if each square of a finite subfamily is the preimage under one fixed map `f` of a
square `sq (ctr' i) (ang' i) 1` lying in the box with centre height `≤ h`, the subfamily has at
most three members. -/
lemma WallStripBound.card_le_three_of_map {h : ℝ} (hW : WallStripBound h) {t : ℝ} (ht : t ≤ 4)
    {ι : Type*} (s : Finset ι) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (ctr' : ι → ℝ × ℝ) (ang' : ι → ℝ) (f : ℝ × ℝ → ℝ × ℝ)
    (hf : ∀ i ∈ s, sq (ctr' i) (ang' i) 1 = f ⁻¹' sq (ctr i) (ang i) 1)
    (hin : ∀ i ∈ s, sq (ctr' i) (ang' i) 1 ⊆ box t)
    (hh : ∀ i ∈ s, (ctr' i).2 ≤ h)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    s.card ≤ 3 := by
  let e := s.equivFin
  refine hW s.card t ht (fun k => ctr' (e.symm k)) (fun k => ang' (e.symm k))
    (fun k => hin _ (e.symm k).2) (fun k => hh _ (e.symm k).2) ?_
  intro k l hkl
  have hne : (e.symm k : ι) ≠ e.symm l := fun h => hkl (e.symm.injective (Subtype.ext h))
  rw [hf _ (e.symm k).2, hf _ (e.symm l).2]
  exact (hdisj _ (e.symm k).2 _ (e.symm l).2 hne).preimage f

/-- **Bottom wall, Finset form.**  Any finite subfamily (of a family indexed by any type) of
closed, pairwise disjoint unit squares in `box t`, `t ≤ 4`, with centre heights `≤ h`, has at most
three members. -/
theorem WallStripBound.card_le_three {h : ℝ} (hW : WallStripBound h) {t : ℝ} (ht : t ≤ 4)
    {ι : Type*} (s : Finset ι) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ box t)
    (hh : ∀ i ∈ s, (ctr i).2 ≤ h)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    s.card ≤ 3 :=
  hW.card_le_three_of_map ht s ctr ang ctr ang id (fun _ _ => rfl) hin hh hdisj

/-! ## 2. Slabs -/

/-- The vertical slab `[a,b] × [0,t]`.  (`vslab 0 t t` is `box t`.) -/
def vslab (a b t : ℝ) : Set (ℝ × ℝ) := {p | a ≤ p.1 ∧ p.1 ≤ b ∧ 0 ≤ p.2 ∧ p.2 ≤ t}

lemma box_eq_vslab (t : ℝ) : box t = vslab 0 t t := rfl

/-- The support-function fact for a slab: a closed unit square lies in `[a,b] × [0,t]` exactly
when its centre lies in `[a + w/2, b - w/2] × [w/2, t - w/2]`, `w = |cos θ| + |sin θ|`. -/
theorem sq_subset_vslab_iff (a b t : ℝ) (c : ℝ × ℝ) (θ : ℝ) :
    sq c θ 1 ⊆ vslab a b t ↔
      a + wid θ / 2 ≤ c.1 ∧ c.1 ≤ b - wid θ / 2 ∧ wid θ / 2 ≤ c.2 ∧ c.2 ≤ t - wid θ / 2 := by
  constructor
  · intro hsub
    set sc : ℝ := if 0 ≤ Real.cos θ then 1 else -1 with hsc
    set ss : ℝ := if 0 ≤ Real.sin θ then 1 else -1 with hss
    have hscC : sc * Real.cos θ = |Real.cos θ| := by
      rw [hsc]; split_ifs with h
      · rw [abs_of_nonneg h]; ring
      · rw [abs_of_neg (not_le.mp h)]; ring
    have hssS : ss * Real.sin θ = |Real.sin θ| := by
      rw [hss]; split_ifs with h
      · rw [abs_of_nonneg h]; ring
      · rw [abs_of_neg (not_le.mp h)]; ring
    have habs : ∀ t : ℝ, (t = 1 ∨ t = -1) → |(-t) / 2| ≤ 1 / 2 ∧ |t / 2| ≤ 1 / 2 := by
      rintro t (rfl | rfl) <;> norm_num
    have hsc' : sc = 1 ∨ sc = -1 := by rw [hsc]; split_ifs <;> simp
    have hss' : ss = 1 ∨ ss = -1 := by rw [hss]; split_ifs <;> simp
    have hA := habs sc hsc'
    have hB := habs ss hss'
    refine ⟨?_, ?_, ?_, ?_⟩
    · have h := (hsub (mem_sq_offset c θ (-sc / 2) (ss / 2) hA.1 hB.2)).1
      simp only [wid] at *
      nlinarith [hscC, hssS, h]
    · have h := (hsub (mem_sq_offset c θ (sc / 2) (-ss / 2) hA.2 hB.1)).2.1
      simp only [wid] at *
      nlinarith [hscC, hssS, h]
    · have h := (hsub (mem_sq_offset c θ (-ss / 2) (-sc / 2) hB.1 hA.1)).2.2.1
      simp only [wid] at *
      nlinarith [hscC, hssS, h]
    · have h := (hsub (mem_sq_offset c θ (ss / 2) (sc / 2) hB.2 hA.2)).2.2.2
      simp only [wid] at *
      nlinarith [hscC, hssS, h]
  · rintro ⟨h1, h2, h3, h4⟩ p hp
    obtain ⟨hx, hy⟩ := abs_sub_le_wid hp
    rw [abs_le] at hx hy
    exact ⟨by linarith [hx.1], by linarith [hx.2], by linarith [hy.1], by linarith [hy.2]⟩

/-- The axis width of a unit square is at least `1`. -/
lemma one_le_wid (θ : ℝ) : 1 ≤ wid θ := by
  have hC := abs_nonneg (Real.cos θ)
  have hS := abs_nonneg (Real.sin θ)
  have hCS : |Real.cos θ| ^ 2 + |Real.sin θ| ^ 2 = 1 := by
    rw [sq_abs, sq_abs]; exact Real.cos_sq_add_sin_sq θ
  simp only [wid]
  nlinarith [mul_nonneg hC hS]

/-- A unit square in `[a,b] × [0,t]` has centre `x ∈ [a + 1/2, b - 1/2]`. -/
lemma centre_mem_of_sq_subset_vslab {a b t : ℝ} {c : ℝ × ℝ} {θ : ℝ}
    (hin : sq c θ 1 ⊆ vslab a b t) : a + 1 / 2 ≤ c.1 ∧ c.1 ≤ b - 1 / 2 := by
  obtain ⟨h1, h2, -, -⟩ := (sq_subset_vslab_iff a b t c θ).mp hin
  have := one_le_wid θ
  exact ⟨by linarith, by linarith⟩

/-! ### Translation -/

/-- Translating a square horizontally by `-u` is the preimage under `p ↦ p + (u,0)`. -/
lemma sq_shift_eq (c : ℝ × ℝ) (θ u : ℝ) :
    sq (c.1 - u, c.2) θ 1 = (fun p : ℝ × ℝ => (p.1 + u, p.2)) ⁻¹' sq c θ 1 := by
  ext p
  simp only [sq, coord, Set.mem_preimage, Set.mem_ofPred_eq]
  have e : p.1 + u - c.1 = p.1 - (c.1 - u) := by ring
  rw [e]

/-- **Translation of a slab.**  A unit square lies in `[u, u+w] × [0,t]` iff its translate by
`-u` lies in `[0,w] × [0,t]`. -/
lemma sq_subset_vslab_shift (u w t : ℝ) (c : ℝ × ℝ) (θ : ℝ) :
    sq c θ 1 ⊆ vslab u (u + w) t ↔ sq (c.1 - u, c.2) θ 1 ⊆ vslab 0 w t := by
  rw [sq_subset_vslab_iff, sq_subset_vslab_iff]
  simp only
  constructor
  · rintro ⟨h1, h2, h3, h4⟩; exact ⟨by linarith, by linarith, h3, h4⟩
  · rintro ⟨h1, h2, h3, h4⟩; exact ⟨by linarith, by linarith, h3, h4⟩

/-- Translation preserves (closed) disjointness. -/
lemma disjoint_sq_shift {c c' : ℝ × ℝ} {θ θ' : ℝ} (u : ℝ)
    (hd : Disjoint (sq c θ 1) (sq c' θ' 1)) :
    Disjoint (sq (c.1 - u, c.2) θ 1) (sq (c'.1 - u, c'.2) θ' 1) := by
  rw [sq_shift_eq, sq_shift_eq]; exact hd.preimage _

/-! ### The two walls of a slab -/

/-- The left wall `x = a`, turned into the bottom wall: `q ↦ (q.2 + a, q.1)`. -/
lemma sq_leftT_eq (a : ℝ) (c : ℝ × ℝ) (θ : ℝ) :
    sq (c.2, c.1 - a) (-θ) 1 = (fun q : ℝ × ℝ => (q.2 + a, q.1)) ⁻¹' sq c θ 1 := by
  ext q
  simp only [sq, Set.mem_preimage, Set.mem_ofPred_eq]
  have h1 : (coord (c.2, c.1 - a) (-θ) q).1 = (coord c θ (q.2 + a, q.1)).2 := by
    simp only [coord, Real.cos_neg, Real.sin_neg]; ring
  have h2 : (coord (c.2, c.1 - a) (-θ) q).2 = (coord c θ (q.2 + a, q.1)).1 := by
    simp only [coord, Real.cos_neg, Real.sin_neg]; ring
  rw [h1, h2, and_comm]

/-- The right wall `x = b`, turned into the bottom wall: `q ↦ (b - q.2, q.1)`. -/
lemma sq_rightT_eq (b : ℝ) (c : ℝ × ℝ) (θ : ℝ) :
    sq (c.2, b - c.1) θ 1 = (fun q : ℝ × ℝ => (b - q.2, q.1)) ⁻¹' sq c θ 1 := by
  ext q
  simp only [sq, Set.mem_preimage, Set.mem_ofPred_eq]
  have h1 : (coord (c.2, b - c.1) θ q).1 = (coord c θ (b - q.2, q.1)).2 := by
    simp only [coord]; ring
  have h2 : (coord (c.2, b - c.1) θ q).2 = -(coord c θ (b - q.2, q.1)).1 := by
    simp only [coord]; ring
  rw [h1, h2, abs_neg, and_comm]

/-- **Left wall of a slab.**  In `[a,b] × [0,t]` with `0 ≤ a`, `b ≤ t ≤ 4`, at most three closed
pairwise disjoint unit squares have centre `x ≤ a + h`. -/
theorem slab_left_card_le_three {h : ℝ} (hW : WallStripBound h) {a b t : ℝ}
    (ht : t ≤ 4) (ha : 0 ≤ a) (hb : b ≤ t)
    {ι : Type*} (s : Finset ι) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ vslab a b t)
    (hx : ∀ i ∈ s, (ctr i).1 ≤ a + h)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    s.card ≤ 3 := by
  refine hW.card_le_three_of_map ht s ctr ang (fun i => ((ctr i).2, (ctr i).1 - a))
    (fun i => -ang i) (fun q => (q.2 + a, q.1)) (fun i _ => sq_leftT_eq a _ _) ?_ ?_ hdisj
  · intro i hi
    obtain ⟨h1, h2, h3, h4⟩ := (sq_subset_vslab_iff a b t _ _).mp (hin i hi)
    rw [box_eq_vslab, sq_subset_vslab_iff, wid_neg]
    simp only
    exact ⟨by linarith, by linarith, by linarith, by linarith⟩
  · intro i hi; simp only; linarith [hx i hi]

/-- **Right wall of a slab.**  In `[a,b] × [0,t]` with `0 ≤ a`, `b ≤ t ≤ 4`, at most three closed
pairwise disjoint unit squares have centre `x ≥ b - h`. -/
theorem slab_right_card_le_three {h : ℝ} (hW : WallStripBound h) {a b t : ℝ}
    (ht : t ≤ 4) (ha : 0 ≤ a) (hb : b ≤ t)
    {ι : Type*} (s : Finset ι) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ vslab a b t)
    (hx : ∀ i ∈ s, b - h ≤ (ctr i).1)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    s.card ≤ 3 := by
  refine hW.card_le_three_of_map ht s ctr ang (fun i => ((ctr i).2, b - (ctr i).1))
    ang (fun q => (b - q.2, q.1)) (fun i _ => sq_rightT_eq b _ _) ?_ ?_ hdisj
  · intro i hi
    obtain ⟨h1, h2, h3, h4⟩ := (sq_subset_vslab_iff a b t _ _).mp (hin i hi)
    rw [box_eq_vslab, sq_subset_vslab_iff]
    simp only
    exact ⟨by linarith, by linarith, by linarith, by linarith⟩
  · intro i hi; simp only; linarith [hx i hi]

/-- **Left wall of the box** (`Fin N` form): at most three centres with `x ≤ h`. -/
theorem left_wall_le_three {h : ℝ} (hW : WallStripBound h) {N : ℕ} {t : ℝ} (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ box t)
    (hx : ∀ i, (ctr i).1 ≤ h)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    N ≤ 3 := by
  have := slab_left_card_le_three hW ht le_rfl le_rfl (univ : Finset (Fin N)) ctr ang
    (fun i _ => hin i) (fun i _ => by linarith [hx i]) (fun i _ j _ hij => hdisj i j hij)
  simpa using this

/-- **Right wall of a rectangle** `[0,w] × [0,t]`, `w ≤ t ≤ 4` (`Fin N` form): at most three
centres with `x ≥ w - h`. -/
theorem right_wall_le_three {h : ℝ} (hW : WallStripBound h) {N : ℕ} {w t : ℝ}
    (hwt : w ≤ t) (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab 0 w t)
    (hx : ∀ i, w - h ≤ (ctr i).1)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    N ≤ 3 := by
  have := slab_right_card_le_three hW ht le_rfl hwt (univ : Finset (Fin N)) ctr ang
    (fun i _ => hin i) (fun i _ => hx i) (fun i _ j _ hij => hdisj i j hij)
  simpa using this

/-! ## 3. Slab capacities -/

/-- **R1, Finset form.**  A slab holding a unit square has width at least `1`. -/
theorem slab_one_le_width {a b t : ℝ} {ι : Type*} (s : Finset ι) (hs : s.Nonempty)
    (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ vslab a b t) :
    1 ≤ b - a := by
  obtain ⟨i, hi⟩ := hs
  have := centre_mem_of_sq_subset_vslab (hin i hi)
  linarith [this.1, this.2]

/-- **R4, Finset form.**  Four closed pairwise disjoint unit squares in a slab
`[a,b] × [0,t] ⊆ box t`, `t ≤ 4`, need width `b - a > h + 1/2`. -/
theorem slab_width_of_four {h : ℝ} (hW : WallStripBound h) {a b t : ℝ}
    (ht : t ≤ 4) (ha : 0 ≤ a) (hb : b ≤ t)
    {ι : Type*} (s : Finset ι) (hs : 4 ≤ s.card) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ vslab a b t)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    h + 1 / 2 < b - a := by
  by_contra hw
  push Not at hw
  have := slab_left_card_le_three hW ht ha hb s ctr ang hin
    (fun i hi => by linarith [(centre_mem_of_sq_subset_vslab (hin i hi)).2]) hdisj
  omega

/-- **R7, Finset form.**  Seven closed pairwise disjoint unit squares in a slab
`[a,b] × [0,t] ⊆ box t`, `t ≤ 4`, need width `b - a > 2h`. -/
theorem slab_width_of_seven {h : ℝ} (hW : WallStripBound h) {a b t : ℝ}
    (ht : t ≤ 4) (ha : 0 ≤ a) (hb : b ≤ t)
    {ι : Type*} (s : Finset ι) (hs : 7 ≤ s.card) (ctr : ι → ℝ × ℝ) (ang : ι → ℝ)
    (hin : ∀ i ∈ s, sq (ctr i) (ang i) 1 ⊆ vslab a b t)
    (hdisj : ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    2 * h < b - a := by
  by_contra hw
  push Not at hw
  set L := s.filter (fun i => (ctr i).1 ≤ a + h)
  set R := s.filter (fun i => b - h ≤ (ctr i).1)
  have hL : L.card ≤ 3 := slab_left_card_le_three hW ht ha hb L ctr ang
    (fun i hi => hin i (mem_filter.mp hi).1) (fun i hi => (mem_filter.mp hi).2)
    (fun i hi j hj => hdisj i (mem_filter.mp hi).1 j (mem_filter.mp hj).1)
  have hR : R.card ≤ 3 := slab_right_card_le_three hW ht ha hb R ctr ang
    (fun i hi => hin i (mem_filter.mp hi).1) (fun i hi => (mem_filter.mp hi).2)
    (fun i hi j hj => hdisj i (mem_filter.mp hi).1 j (mem_filter.mp hj).1)
  have hsub : s ⊆ L ∪ R := by
    intro i hi
    rw [mem_union, mem_filter, mem_filter]
    by_cases hc : (ctr i).1 ≤ a + h
    · exact Or.inl ⟨hi, hc⟩
    · exact Or.inr ⟨hi, by linarith⟩
  have := (card_le_card hsub).trans (card_union_le L R)
  omega

/-- **R1.**  `N ≥ 1` unit squares in the slab `[u, u+w] × [0,t]` force `w ≥ 1`. -/
theorem slab_R1 {N : ℕ} (hN : 1 ≤ N) {u w t : ℝ}
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab u (u + w) t) :
    1 ≤ w := by
  have := slab_one_le_width (univ : Finset (Fin N)) ⟨⟨0, by omega⟩, mem_univ _⟩ ctr ang
    (fun i _ => hin i)
  linarith

/-- **R4.**  `N ≥ 4` closed pairwise disjoint unit squares in the slab `[u, u+w] × [0,t]`,
`0 ≤ u`, `u + w ≤ t ≤ 4`, force `w > h + 1/2`. -/
theorem slab_R4 {h : ℝ} (hW : WallStripBound h) {N : ℕ} (hN : 4 ≤ N) {u w t : ℝ}
    (hu : 0 ≤ u) (huw : u + w ≤ t) (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab u (u + w) t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    h + 1 / 2 < w := by
  have := slab_width_of_four hW ht hu huw (univ : Finset (Fin N)) (by simpa using hN) ctr ang
    (fun i _ => hin i) (fun i _ j _ hij => hdisj i j hij)
  linarith

/-- **R7.**  `N ≥ 7` closed pairwise disjoint unit squares in the slab `[u, u+w] × [0,t]`,
`0 ≤ u`, `u + w ≤ t ≤ 4`, force `w > 2h`. -/
theorem slab_R7 {h : ℝ} (hW : WallStripBound h) {N : ℕ} (hN : 7 ≤ N) {u w t : ℝ}
    (hu : 0 ≤ u) (huw : u + w ≤ t) (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab u (u + w) t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    2 * h < w := by
  have := slab_width_of_seven hW ht hu huw (univ : Finset (Fin N)) (by simpa using hN) ctr ang
    (fun i _ => hin i) (fun i _ j _ hij => hdisj i j hij)
  linarith

/-! ## 4. Guillotine splits -/

/-- **No guillotine `(k, 4, 8-k)` split.**  Twelve closed pairwise disjoint unit squares in
`[0,t]^2`, `t ≤ 4`; two cuts `0 ≤ u ≤ v ≤ t` such that every square lies in `[0,u] × [0,t]`,
`[u,v] × [0,t]` or `[v,t] × [0,t]`, with exactly four in the middle slab.  If the wall-strip bound
holds at some `h ≥ 7/6`, this is impossible. -/
theorem guillotine_false {h : ℝ} (hW : WallStripBound h) (hh : 7 / 6 ≤ h) {t u v : ℝ}
    (ht : t ≤ 4) (hu : 0 ≤ u) (huv : u ≤ v) (hv : v ≤ t)
    (ctr : Fin 12 → ℝ × ℝ) (ang : Fin 12 → ℝ)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1))
    (hcut : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab 0 u t ∨ sq (ctr i) (ang i) 1 ⊆ vslab u v t ∨
      sq (ctr i) (ang i) 1 ⊆ vslab v t t)
    (hmid : (univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab u v t)).card = 4) :
    False := by
  set L := univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab 0 u t)
  set M := univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab u v t)
  set R := univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab v t t)
  have hd : ∀ (s : Finset (Fin 12)), ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1) :=
    fun _ i _ j _ hij => hdisj i j hij
  have hinL : ∀ i ∈ L, sq (ctr i) (ang i) 1 ⊆ vslab 0 u t := fun i hi => (mem_filter.mp hi).2
  have hinM : ∀ i ∈ M, sq (ctr i) (ang i) 1 ⊆ vslab u v t := fun i hi => (mem_filter.mp hi).2
  have hinR : ∀ i ∈ R, sq (ctr i) (ang i) 1 ⊆ vslab v t t := fun i hi => (mem_filter.mp hi).2
  have hcount : 8 ≤ L.card + R.card := by
    have hsub : (univ : Finset (Fin 12)) ⊆ L ∪ M ∪ R := by
      intro i _
      simp only [L, M, R, mem_union, mem_filter, mem_univ, true_and]
      rcases hcut i with h | h | h
      · exact Or.inl (Or.inl h)
      · exact Or.inl (Or.inr h)
      · exact Or.inr h
    have h1 := card_le_card hsub
    have h2 := card_union_le (L ∪ M) R
    have h3 := card_union_le L M
    simp only [card_univ, Fintype.card_fin] at h1
    omega
  -- the middle slab
  have wM : h + 1 / 2 < v - u :=
    slab_width_of_four hW ht hu hv M (by omega) ctr ang hinM (hd M)
  -- left slab bounds
  have hL1 : 1 ≤ L.card → 1 ≤ u - 0 := fun hk =>
    slab_one_le_width L (card_pos.mp (by omega)) ctr ang hinL
  have hL4 : 4 ≤ L.card → h + 1 / 2 < u - 0 := fun hk =>
    slab_width_of_four hW ht le_rfl (huv.trans hv) L hk ctr ang hinL (hd L)
  have hL7 : 7 ≤ L.card → 2 * h < u - 0 := fun hk =>
    slab_width_of_seven hW ht le_rfl (huv.trans hv) L hk ctr ang hinL (hd L)
  have hR1 : 1 ≤ R.card → 1 ≤ t - v := fun hk =>
    slab_one_le_width R (card_pos.mp (by omega)) ctr ang hinR
  have hR4 : 4 ≤ R.card → h + 1 / 2 < t - v := fun hk =>
    slab_width_of_four hW ht (hu.trans huv) le_rfl R hk ctr ang hinR (hd R)
  have hR7 : 7 ≤ R.card → 2 * h < t - v := fun hk =>
    slab_width_of_seven hW ht (hu.trans huv) le_rfl R hk ctr ang hinR (hd R)
  rcases Nat.lt_or_ge L.card 1 with k0 | k1
  · have := hR7 (by omega); linarith
  rcases Nat.lt_or_ge R.card 1 with m0 | m1
  · have := hL7 (by omega); linarith
  rcases Nat.lt_or_ge L.card 4 with k3 | k4
  · have := hR4 (by omega); have := hL1 k1; linarith
  rcases Nat.lt_or_ge R.card 4 with m3 | m4
  · have := hL4 (by omega); have := hR1 m1; linarith
  · have := hL4 k4; have := hR4 m4; linarith

/-! ## 5. The line sweep -/

/-- **Line sweep.**  Twelve closed pairwise disjoint unit squares in `[0,t]^2`, `t ≤ 4`.  For any
`u` with `u ≤ h + 1/2` and `t - u ≤ 2h`, at least three of the squares lie neither in
`[0,u] × [0,t]` nor in `[u,t] × [0,t]` (i.e. they cross the line `x = u`). -/
theorem line_sweep {h : ℝ} (hW : WallStripBound h) {t u : ℝ} (ht : t ≤ 4)
    (hu1 : u ≤ h + 1 / 2) (hu2 : t - u ≤ 2 * h)
    (ctr : Fin 12 → ℝ × ℝ) (ang : Fin 12 → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ box t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    3 ≤ (univ.filter (fun i => ¬ sq (ctr i) (ang i) 1 ⊆ vslab 0 u t ∧
      ¬ sq (ctr i) (ang i) 1 ⊆ vslab u t t)).card := by
  set L := univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab 0 u t)
  set R := univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab u t t)
  set S := univ.filter (fun i => ¬ sq (ctr i) (ang i) 1 ⊆ vslab 0 u t ∧
      ¬ sq (ctr i) (ang i) 1 ⊆ vslab u t t)
  have hd : ∀ (s : Finset (Fin 12)), ∀ i ∈ s, ∀ j ∈ s, i ≠ j →
      Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1) :=
    fun _ i _ j _ hij => hdisj i j hij
  -- clamp the cut into `[0,t]` (the squares lie in the box anyway)
  have hinL : ∀ i ∈ L, sq (ctr i) (ang i) 1 ⊆ vslab 0 (min u t) t := by
    intro i hi p hp
    have h1 := (mem_filter.mp hi).2 hp
    have h2 := hin i hp
    exact ⟨h1.1, le_min h1.2.1 h2.2.1, h1.2.2⟩
  have hinR : ∀ i ∈ R, sq (ctr i) (ang i) 1 ⊆ vslab (max u 0) t t := by
    intro i hi p hp
    have h1 := (mem_filter.mp hi).2 hp
    have h2 := hin i hp
    exact ⟨max_le h1.1 h2.1, h1.2⟩
  have hL : L.card ≤ 3 := by
    by_contra hk
    have := slab_width_of_four hW ht le_rfl (min_le_right u t) L (by omega) ctr ang hinL (hd L)
    have := min_le_left u t
    linarith
  have hR : R.card ≤ 6 := by
    by_contra hk
    have := slab_width_of_seven hW ht (le_max_right u 0) le_rfl R (by omega) ctr ang hinR (hd R)
    have := le_max_left u 0
    linarith
  have hsub : (univ : Finset (Fin 12)) ⊆ L ∪ R ∪ S := by
    intro i _
    simp only [L, R, S, mem_union, mem_filter, mem_univ, true_and]
    tauto
  have h1 := card_le_card hsub
  have h2 := card_union_le (L ∪ R) S
  have h3 := card_union_le L R
  simp only [card_univ, Fintype.card_fin] at h1
  omega

/-! ## 6. At the height available in Lean, `h = (3√2 - 2)/2` -/

/-- R4 at `h = (3√2 - 2)/2`: four squares in a slab need width `> 3√2/2 - 1/2 ≈ 1.621`. -/
theorem slab_R4_line {N : ℕ} (hN : 4 ≤ N) {u w t : ℝ}
    (hu : 0 ≤ u) (huw : u + w ≤ t) (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab u (u + w) t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    (3 * Real.sqrt 2 - 2) / 2 + 1 / 2 < w :=
  slab_R4 wallStripBound_line hN hu huw ht ctr ang hin hdisj

/-- R7 at `h = (3√2 - 2)/2`: seven squares in a slab need width `> 3√2 - 2 ≈ 2.243`. -/
theorem slab_R7_line {N : ℕ} (hN : 7 ≤ N) {u w t : ℝ}
    (hu : 0 ≤ u) (huw : u + w ≤ t) (ht : t ≤ 4)
    (ctr : Fin N → ℝ × ℝ) (ang : Fin N → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab u (u + w) t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    3 * Real.sqrt 2 - 2 < w := by
  have := slab_R7 wallStripBound_line hN hu huw ht ctr ang hin hdisj
  linarith

/-- The line sweep at `h = (3√2 - 2)/2`. -/
theorem line_sweep_line {t u : ℝ} (ht : t ≤ 4)
    (hu1 : u ≤ (3 * Real.sqrt 2 - 2) / 2 + 1 / 2) (hu2 : t - u ≤ 3 * Real.sqrt 2 - 2)
    (ctr : Fin 12 → ℝ × ℝ) (ang : Fin 12 → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ box t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1)) :
    3 ≤ (univ.filter (fun i => ¬ sq (ctr i) (ang i) 1 ⊆ vslab 0 u t ∧
      ¬ sq (ctr i) (ang i) 1 ⊆ vslab u t t)).card :=
  line_sweep wallStripBound_line ht hu1 (by linarith) ctr ang hin hdisj

end SquarePacking
