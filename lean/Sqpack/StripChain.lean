import Sqpack.StripRect

/-!
# The wall-strip lemma past the one-line threshold: the hand-proved part

`tasks/n12-chat-lemmas/strip-demo/PROOF.md` §2–3.  Everything here is generic: no certificate data.

* §1  A closed unit square `sq c θ 1` is `sqCS x c C S` for a normalised rotation `C > 0`, `S ≥ 0`,
  `C² + S² = 1` (`exists_CS`: the four quarter-turns of `(cos θ, sin θ)`); in the box `[0,T]²` its
  centre satisfies `p ≤ c`, `p ≤ x ≤ T - p` with `p = (C+S)/2` (`sqCS_box`).
* §2  Half-chords.  On the line at height `c + e` (`|e| ≤ p`) the square's chord is
  `[x - hL C S e, x + hR C S e]` with (`PROOF.md` §2, `notes/chord-lemma.md` Lemma 1)
  `hR = min((1/2 - eS)/C, (eC + 1/2)/S)`, `hL = min((1/2 + eS)/C, (1/2 - eC)/S)` (both `1/2` when
  `S = 0`); the two endpoints lie in the square (`mem_right`, `mem_left`).
* §3  Lemma 1 (one order on every line): `order_transfer`, for any two disjoint convex sets.
  Lemma 2 (chain inequality): `chain_lt`, for four squares in order and any three separation heights
  at which every square is crossed.
* §4  The tan-half-angle parametrisation `t = S/(1+C) ∈ [0,1)` (`Ct`, `St`, `pT`, `DT`), node values
  `nodeT`, and Case A (`nodeT_line_ge_one`: on the line `y = a`, `a ≤ √2 - 1/2`, every non-H square
  has node value `≥ 1`).
* §5  `wallStrip_of_cert`: the wall-strip bound at height `h` from four abstract certificate facts
  (coverage, H-flags, and the chain condition with leaf bounds).  `Sqpack/Strip117.lean` supplies
  them for `h = 117/100`.
-/

open Finset

namespace SquarePacking
namespace StripChain

/-! ## 1. Normalised squares -/

/-- The closed unit square with centre `(x, c)` and rotation `(C, S)` (`C = cos θ`, `S = sin θ`). -/
def sqCS (x c C S : ℝ) : Set (ℝ × ℝ) :=
  {q | |(q.1 - x) * C + (q.2 - c) * S| ≤ 1 / 2 ∧ |-(q.1 - x) * S + (q.2 - c) * C| ≤ 1 / 2}

lemma sq_eq_sqCS_of (c : ℝ × ℝ) (θ C S : ℝ)
    (h : (C = Real.cos θ ∧ S = Real.sin θ) ∨ (C = Real.sin θ ∧ S = -Real.cos θ) ∨
      (C = -Real.cos θ ∧ S = -Real.sin θ) ∨ (C = -Real.sin θ ∧ S = Real.cos θ)) :
    sq c θ 1 = sqCS c.1 c.2 C S := by
  ext q
  simp only [sq, sqCS, coord, Set.mem_ofPred_eq]
  rcases h with ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩ | ⟨rfl, rfl⟩
  · rw [neg_mul]
  · rw [show (q.1 - c.1) * Real.sin θ + (q.2 - c.2) * -Real.cos θ =
        -(-(q.1 - c.1) * Real.sin θ + (q.2 - c.2) * Real.cos θ) by ring,
      show -(q.1 - c.1) * -Real.cos θ + (q.2 - c.2) * Real.sin θ =
        (q.1 - c.1) * Real.cos θ + (q.2 - c.2) * Real.sin θ by ring, abs_neg, and_comm]
  · rw [show (q.1 - c.1) * -Real.cos θ + (q.2 - c.2) * -Real.sin θ =
        -((q.1 - c.1) * Real.cos θ + (q.2 - c.2) * Real.sin θ) by ring,
      show -(q.1 - c.1) * -Real.sin θ + (q.2 - c.2) * -Real.cos θ =
        -(-(q.1 - c.1) * Real.sin θ + (q.2 - c.2) * Real.cos θ) by ring, abs_neg, abs_neg]
  · rw [show (q.1 - c.1) * -Real.sin θ + (q.2 - c.2) * Real.cos θ =
        -(q.1 - c.1) * Real.sin θ + (q.2 - c.2) * Real.cos θ by ring,
      show -(q.1 - c.1) * Real.cos θ + (q.2 - c.2) * -Real.sin θ =
        -((q.1 - c.1) * Real.cos θ + (q.2 - c.2) * Real.sin θ) by ring, abs_neg, and_comm]

/-- Every closed unit square is `sqCS` for a rotation `(C, S)` with `C > 0`, `S ≥ 0`. -/
lemma exists_CS (θ : ℝ) : ∃ C S : ℝ, 0 < C ∧ 0 ≤ S ∧ C ^ 2 + S ^ 2 = 1 ∧
    ∀ c : ℝ × ℝ, sq c θ 1 = sqCS c.1 c.2 C S := by
  have h1 : Real.cos θ ^ 2 + Real.sin θ ^ 2 = 1 := Real.cos_sq_add_sin_sq θ
  rcases lt_or_ge 0 (Real.cos θ) with hc | hc
  · rcases le_or_gt 0 (Real.sin θ) with hs | hs
    · exact ⟨_, _, hc, hs, h1, fun c => sq_eq_sqCS_of c θ _ _ (Or.inl ⟨rfl, rfl⟩)⟩
    · exact ⟨-Real.sin θ, Real.cos θ, by linarith, hc.le, by nlinarith,
        fun c => sq_eq_sqCS_of c θ _ _ (Or.inr (Or.inr (Or.inr ⟨rfl, rfl⟩)))⟩
  · rcases lt_or_ge 0 (Real.sin θ) with hs | hs
    · exact ⟨Real.sin θ, -Real.cos θ, hs, by linarith, by nlinarith,
        fun c => sq_eq_sqCS_of c θ _ _ (Or.inr (Or.inl ⟨rfl, rfl⟩))⟩
    · rcases lt_or_eq_of_le hc with hc' | hc'
      · exact ⟨-Real.cos θ, -Real.sin θ, by linarith, by linarith, by nlinarith,
          fun c => sq_eq_sqCS_of c θ _ _ (Or.inr (Or.inr (Or.inl ⟨rfl, rfl⟩)))⟩
      · -- `cos θ = 0`, `sin θ ≤ 0`: so `sin θ = -1`
        have hs1 : Real.sin θ < 0 := by
          rcases lt_or_eq_of_le hs with h | h
          · exact h
          · exfalso; rw [hc', h] at h1; norm_num at h1
        exact ⟨-Real.sin θ, Real.cos θ, by linarith, by linarith, by nlinarith,
          fun c => sq_eq_sqCS_of c θ _ _ (Or.inr (Or.inr (Or.inr ⟨rfl, rfl⟩)))⟩

/-- The point with local coordinates `(u, v)`. -/
lemma mem_sqCS_local {x c C S u v : ℝ} (h1 : C ^ 2 + S ^ 2 = 1) (hu : |u| ≤ 1 / 2)
    (hv : |v| ≤ 1 / 2) : (x + (u * C - v * S), c + (u * S + v * C)) ∈ sqCS x c C S := by
  simp only [sqCS, Set.mem_ofPred_eq]
  have e1 : (x + (u * C - v * S) - x) * C + (c + (u * S + v * C) - c) * S = u := by
    linear_combination u * h1
  have e2 : -(x + (u * C - v * S) - x) * S + (c + (u * S + v * C) - c) * C = v := by
    linear_combination v * h1
  rw [e1, e2]; exact ⟨hu, hv⟩

/-- A normalised square in `[0,T]²`: `p ≤ c`, `p ≤ x`, `x + p ≤ T`, `p = (C+S)/2`. -/
lemma sqCS_box {x c C S T : ℝ} (h1 : C ^ 2 + S ^ 2 = 1) (hin : sqCS x c C S ⊆ box T) :
    (C + S) / 2 ≤ c ∧ (C + S) / 2 ≤ x ∧ x + (C + S) / 2 ≤ T := by
  have hh : |(1 : ℝ) / 2| ≤ 1 / 2 := by rw [abs_of_pos (by norm_num)]
  have hh' : |-(1 : ℝ) / 2| ≤ 1 / 2 := by rw [neg_div, abs_neg, abs_of_pos (by norm_num)]
  have hlow := (hin (mem_sqCS_local (x := x) (c := c) h1 hh' hh')).2.2.1
  have hleft := (hin (mem_sqCS_local (x := x) (c := c) h1 hh' hh)).1
  have hright := (hin (mem_sqCS_local (x := x) (c := c) h1 hh hh')).2.1
  simp only at hlow hleft hright
  refine ⟨by linarith, by linarith, by linarith⟩

/-- Normalised squares are convex. -/
lemma sqCS_convex (x c C S : ℝ) : Convex ℝ (sqCS x c C S) := by
  intro p hp q hq a b ha hb hab
  simp only [sqCS, Set.mem_ofPred_eq, Prod.fst_add, Prod.snd_add, Prod.smul_fst, Prod.smul_snd,
    smul_eq_mul] at hp hq ⊢
  have e1 : (a * p.1 + b * q.1 - x) * C + (a * p.2 + b * q.2 - c) * S =
      a * ((p.1 - x) * C + (p.2 - c) * S) + b * ((q.1 - x) * C + (q.2 - c) * S) := by
    linear_combination (x * C + c * S) * hab
  have e2 : -(a * p.1 + b * q.1 - x) * S + (a * p.2 + b * q.2 - c) * C =
      a * (-(p.1 - x) * S + (p.2 - c) * C) + b * (-(q.1 - x) * S + (q.2 - c) * C) := by
    linear_combination (-(x * S) + c * C) * hab
  rw [e1, e2]
  obtain ⟨hp1, hp2⟩ := hp
  obtain ⟨hq1, hq2⟩ := hq
  rw [abs_le] at hp1 hp2 hq1 hq2 ⊢
  rw [abs_le]
  refine ⟨⟨?_, ?_⟩, ?_, ?_⟩ <;>
    nlinarith [mul_le_mul_of_nonneg_left hp1.1 ha, mul_le_mul_of_nonneg_left hp1.2 ha,
      mul_le_mul_of_nonneg_left hp2.1 ha, mul_le_mul_of_nonneg_left hp2.2 ha,
      mul_le_mul_of_nonneg_left hq1.1 hb, mul_le_mul_of_nonneg_left hq1.2 hb,
      mul_le_mul_of_nonneg_left hq2.1 hb, mul_le_mul_of_nonneg_left hq2.2 hb]

/-! ## 2. Half-chords -/

/-- Right half-chord on the line at offset `e` from the centre. -/
noncomputable def hR (C S e : ℝ) : ℝ :=
  if S = 0 then 1 / 2 else min ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)

/-- Left half-chord on the line at offset `e` from the centre. -/
noncomputable def hL (C S e : ℝ) : ℝ :=
  if S = 0 then 1 / 2 else min ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)

lemma C_eq_one {C S : ℝ} (hC : 0 < C) (h1 : C ^ 2 + S ^ 2 = 1) (hS : S = 0) : C = 1 := by
  subst hS; nlinarith

/-- The right end of the chord lies in the square. -/
lemma mem_right {x c C S e : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1)
    (he : |e| ≤ (C + S) / 2) : (x + hR C S e, c + e) ∈ sqCS x c C S := by
  simp only [sqCS, Set.mem_ofPred_eq, add_sub_cancel_left]
  rw [abs_le] at he
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · have hC1 := C_eq_one hC h1 hS0.symm
    rw [hR, if_pos hS0.symm, ← hS0, hC1]
    rw [← hS0, hC1] at he
    refine ⟨by norm_num, ?_⟩
    rw [abs_le]; constructor <;> linarith [he.1, he.2]
  rw [hR, if_neg hSp.ne']
  set ξ := min ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S) with hξ
  have hξ1 : ξ * C ≤ 1 / 2 - e * S := by
    have := min_le_left ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)
    rw [← hξ, le_div_iff₀ hC] at this; exact this
  have hξ2 : ξ * S ≤ e * C + 1 / 2 := by
    have := min_le_right ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)
    rw [← hξ, le_div_iff₀ hSp] at this; exact this
  rw [abs_le, abs_le]
  rcases min_choice ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S) with h | h <;> rw [← hξ] at h
  · have hx : ξ * C = 1 / 2 - e * S := by rw [h]; field_simp
    have key : (-ξ * S + e * C - 1 / 2) * C = e - (C + S) / 2 := by
      linear_combination (-S) * hx + e * h1
    have : -ξ * S + e * C - 1 / 2 ≤ 0 := by
      by_contra hc; push Not at hc; nlinarith
    exact ⟨⟨by linarith, by linarith⟩, by linarith, by linarith⟩
  · have hx : ξ * S = e * C + 1 / 2 := by rw [h]; field_simp
    have key : (ξ * C + e * S + 1 / 2) * S = e + (C + S) / 2 := by
      linear_combination C * hx + e * h1
    have : 0 ≤ ξ * C + e * S + 1 / 2 := by
      by_contra hc; push Not at hc; nlinarith
    exact ⟨⟨by linarith, by linarith⟩, by linarith, by linarith⟩

/-- The left end of the chord lies in the square. -/
lemma mem_left {x c C S e : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1)
    (he : |e| ≤ (C + S) / 2) : (x - hL C S e, c + e) ∈ sqCS x c C S := by
  simp only [sqCS, Set.mem_ofPred_eq, add_sub_cancel_left, sub_sub_cancel_left]
  rw [abs_le] at he
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · have hC1 := C_eq_one hC h1 hS0.symm
    rw [hL, if_pos hS0.symm, ← hS0, hC1]
    rw [← hS0, hC1] at he
    refine ⟨by norm_num, ?_⟩
    rw [abs_le]; constructor <;> linarith [he.1, he.2]
  rw [hL, if_neg hSp.ne']
  set ξ := min ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S) with hξ
  have hξ1 : ξ * C ≤ 1 / 2 + e * S := by
    have := min_le_left ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)
    rw [← hξ, le_div_iff₀ hC] at this; exact this
  have hξ2 : ξ * S ≤ 1 / 2 - e * C := by
    have := min_le_right ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)
    rw [← hξ, le_div_iff₀ hSp] at this; exact this
  rw [abs_le, abs_le]
  rcases min_choice ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S) with h | h <;> rw [← hξ] at h
  · have hx : ξ * C = 1 / 2 + e * S := by rw [h]; field_simp
    have key : (ξ * S + e * C + 1 / 2) * C = e + (C + S) / 2 := by
      linear_combination S * hx + e * h1
    have : 0 ≤ ξ * S + e * C + 1 / 2 := by
      by_contra hc; push Not at hc; nlinarith
    exact ⟨⟨by linarith, by linarith⟩, by linarith, by linarith⟩
  · have hx : ξ * S = 1 / 2 - e * C := by rw [h]; field_simp
    have key : (-ξ * C + e * S - 1 / 2) * S = e - (C + S) / 2 := by
      linear_combination (-C) * hx + e * h1
    have : -ξ * C + e * S - 1 / 2 ≤ 0 := by
      by_contra hc; push Not at hc; nlinarith
    exact ⟨⟨by linarith, by linarith⟩, by linarith, by linarith⟩

/-- `hR ≤ p`. -/
lemma hR_le {C S e : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1) :
    hR C S e ≤ (C + S) / 2 := by
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · rw [hR, if_pos hS0.symm, C_eq_one hC h1 hS0.symm, ← hS0]; norm_num
  rw [hR, if_neg hSp.ne']
  have a1 := min_le_left ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)
  have a2 := min_le_right ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)
  set ξ := min ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S)
  rw [le_div_iff₀ hC] at a1
  rw [le_div_iff₀ hSp] at a2
  nlinarith [mul_le_mul_of_nonneg_left a1 hC.le, mul_le_mul_of_nonneg_left a2 hS]

/-- `hL ≤ p`. -/
lemma hL_le {C S e : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1) :
    hL C S e ≤ (C + S) / 2 := by
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · rw [hL, if_pos hS0.symm, C_eq_one hC h1 hS0.symm, ← hS0]; norm_num
  rw [hL, if_neg hSp.ne']
  have a1 := min_le_left ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)
  have a2 := min_le_right ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)
  set ξ := min ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S)
  rw [le_div_iff₀ hC] at a1
  rw [le_div_iff₀ hSp] at a2
  nlinarith [mul_le_mul_of_nonneg_left a1 hC.le, mul_le_mul_of_nonneg_left a2 hS]

/-! ## 3. Lemma 1 (one order on every line) and Lemma 2 (the chain inequality) -/

/-- **Lemma 1.**  Two disjoint convex sets: if on one horizontal line a point of `P` lies left of a
point of `Q`, then on any other horizontal line every point of `P` lies left of every point of
`Q`. -/
lemma order_transfer {P Q : Set (ℝ × ℝ)} (hP : Convex ℝ P) (hQ : Convex ℝ Q) (hd : Disjoint P Q)
    {pa pz qa qz : ℝ × ℝ} (hpa : pa ∈ P) (hpz : pz ∈ P) (hqa : qa ∈ Q) (hqz : qz ∈ Q)
    (hya : pa.2 = qa.2) (hyz : pz.2 = qz.2) (hlt : pa.1 < qa.1) : pz.1 < qz.1 := by
  by_contra hge
  push Not at hge
  set α := qa.1 - pa.1 with hα
  set β := pz.1 - qz.1 with hβ
  have hα0 : 0 < α := by rw [hα]; linarith
  have hβ0 : 0 ≤ β := by rw [hβ]; linarith
  have hs : 0 < α + β := by linarith
  set s := α / (α + β) with hs_def
  have hs0 : 0 ≤ s := div_nonneg hα0.le hs.le
  have hs1 : 1 - s = β / (α + β) := by rw [hs_def]; field_simp; ring
  have hs1' : 0 ≤ 1 - s := by rw [hs1]; exact div_nonneg hβ0 hs.le
  have hX := hP hpa hpz hs1' hs0 (by ring)
  have hY := hQ hqa hqz hs1' hs0 (by ring)
  have heq : (1 - s) • pa + s • pz = (1 - s) • qa + s • qz := by
    ext
    · simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
      rw [hs1, hs_def]
      field_simp
      rw [hα, hβ]; ring
    · simp only [Prod.snd_add, Prod.smul_snd, smul_eq_mul]
      rw [hya, hyz]
  rw [heq] at hX
  exact Set.disjoint_left.mp hd hX hY

/-- The chord gap.  Two disjoint normalised squares `i`, `j`, both crossed by the lines `y = a` and
`y = z`; if on `y = a` the right end of `i`'s chord lies left of the right end of `j`'s, then on
`y = z` the right end of `i`'s chord lies strictly left of the left end of `j`'s. -/
lemma chord_gap {xi ci Ci Si xj cj Cj Sj a z : ℝ}
    (hCi : 0 < Ci) (hSi : 0 ≤ Si) (h1i : Ci ^ 2 + Si ^ 2 = 1)
    (hCj : 0 < Cj) (hSj : 0 ≤ Sj) (h1j : Cj ^ 2 + Sj ^ 2 = 1)
    (hd : Disjoint (sqCS xi ci Ci Si) (sqCS xj cj Cj Sj))
    (hai : |a - ci| ≤ (Ci + Si) / 2) (haj : |a - cj| ≤ (Cj + Sj) / 2)
    (hzi : |z - ci| ≤ (Ci + Si) / 2) (hzj : |z - cj| ≤ (Cj + Sj) / 2)
    (hord : xi + hR Ci Si (a - ci) < xj + hR Cj Sj (a - cj)) :
    xi + hR Ci Si (z - ci) < xj - hL Cj Sj (z - cj) := by
  have hpa := mem_right (x := xi) (c := ci) hCi hSi h1i hai
  have hpz := mem_right (x := xi) (c := ci) hCi hSi h1i hzi
  have hqa := mem_right (x := xj) (c := cj) hCj hSj h1j haj
  have hqz := mem_left (x := xj) (c := cj) hCj hSj h1j hzj
  exact order_transfer (sqCS_convex _ _ _ _) (sqCS_convex _ _ _ _) hd hpa hpz hqa hqz
    (by simp) (by simp) hord

/-- Every square is crossed by the line `y = z`. -/
def Crosses (Y C S : Fin 4 → ℝ) (z : ℝ) : Prop := ∀ k, |z - Y k| ≤ (C k + S k) / 2

/-- **Lemma 2 (chain inequality).**  Four pairwise disjoint normalised squares in `[0,T]²`,
`T ≤ 4`, ordered along the line `y = a`, and three separation heights `z₀, z₁, z₂`; all four
lines cross all four squares.  Then the node values sum to less than `4`. -/
theorem chain_lt (X Y C S : Fin 4 → ℝ) {T a z0 z1 z2 : ℝ}
    (hC : ∀ k, 0 < C k) (hS : ∀ k, 0 ≤ S k) (h1 : ∀ k, C k ^ 2 + S k ^ 2 = 1)
    (hx0 : (C 0 + S 0) / 2 ≤ X 0) (hx3 : X 3 + (C 3 + S 3) / 2 ≤ T) (hT : T ≤ 4)
    (hd : ∀ i j, i ≠ j → Disjoint (sqCS (X i) (Y i) (C i) (S i)) (sqCS (X j) (Y j) (C j) (S j)))
    (ha : Crosses Y C S a) (hz0 : Crosses Y C S z0) (hz1 : Crosses Y C S z1)
    (hz2 : Crosses Y C S z2)
    (o01 : X 0 + hR (C 0) (S 0) (a - Y 0) < X 1 + hR (C 1) (S 1) (a - Y 1))
    (o12 : X 1 + hR (C 1) (S 1) (a - Y 1) < X 2 + hR (C 2) (S 2) (a - Y 2))
    (o23 : X 2 + hR (C 2) (S 2) (a - Y 2) < X 3 + hR (C 3) (S 3) (a - Y 3)) :
    (C 0 + S 0) / 2 + hR (C 0) (S 0) (z0 - Y 0) +
      (hL (C 1) (S 1) (z0 - Y 1) + hR (C 1) (S 1) (z1 - Y 1)) +
      (hL (C 2) (S 2) (z1 - Y 2) + hR (C 2) (S 2) (z2 - Y 2)) +
      (hL (C 3) (S 3) (z2 - Y 3) + (C 3 + S 3) / 2) < 4 := by
  have g01 := chord_gap (hC 0) (hS 0) (h1 0) (hC 1) (hS 1) (h1 1) (hd 0 1 (by decide))
    (ha 0) (ha 1) (hz0 0) (hz0 1) o01
  have g12 := chord_gap (hC 1) (hS 1) (h1 1) (hC 2) (hS 2) (h1 2) (hd 1 2 (by decide))
    (ha 1) (ha 2) (hz1 1) (hz1 2) o12
  have g23 := chord_gap (hC 2) (hS 2) (h1 2) (hC 3) (hS 3) (h1 3) (hd 2 3 (by decide))
    (ha 2) (ha 3) (hz2 2) (hz2 3) o23
  linarith

/-! ## 4. The tan-half-angle parametrisation, node values, and Case A -/

/-- `cos θ` at `t = tan(θ/2)`. -/
noncomputable def Ct (t : ℝ) : ℝ := (1 - t ^ 2) / (1 + t ^ 2)
/-- `sin θ` at `t = tan(θ/2)`. -/
noncomputable def St (t : ℝ) : ℝ := 2 * t / (1 + t ^ 2)
/-- Half-height `p = (C+S)/2`. -/
noncomputable def pT (t : ℝ) : ℝ := (Ct t + St t) / 2
/-- Chord-band half-height `D = p - CS` (`notes/chord-lemma.md` Cor. 2). -/
noncomputable def DT (t : ℝ) : ℝ := pT t - Ct t * St t

/-- A pose at strip height `h`: `t ∈ [0,1)`, `p(t) ≤ c ≤ h`. -/
def Pose (h t c : ℝ) : Prop := 0 ≤ t ∧ t < 1 ∧ pT t ≤ c ∧ c ≤ h

/-- Node value of a square whose left separation is side `l` and right separation side `r`; side
`0` is the container wall (END), side `j ≥ 1` the line `y = zOf j`. -/
noncomputable def nodeT (zOf : ℕ → ℝ) (t c : ℝ) (l r : ℕ) : ℝ :=
  (if l = 0 then pT t else hL (Ct t) (St t) (zOf l - c)) +
    (if r = 0 then pT t else hR (Ct t) (St t) (zOf r - c))

/-- The normalised rotation as a tan-half-angle. -/
lemma tan_half {C S : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1) :
    0 ≤ S / (1 + C) ∧ S / (1 + C) < 1 ∧ Ct (S / (1 + C)) = C ∧ St (S / (1 + C)) = S := by
  have hpos : 0 < 1 + C := by linarith
  have hS1 : S ≤ 1 := by nlinarith
  refine ⟨div_nonneg hS hpos.le, (div_lt_one hpos).mpr (by linarith), ?_, ?_⟩
  · rw [Ct, div_eq_iff (by positivity)]
    field_simp
    linear_combination (-(1 + C)) * h1
  · rw [St, div_eq_iff (by positivity)]
    field_simp
    linear_combination (-S) * h1

lemma Ct_sq_add (t : ℝ) : Ct t ^ 2 + St t ^ 2 = 1 := by
  rw [Ct, St]; field_simp; ring

lemma Ct_pos {t : ℝ} (h0 : 0 ≤ t) (h1 : t < 1) : 0 < Ct t := by
  rw [Ct]; apply div_pos <;> nlinarith

lemma St_nonneg {t : ℝ} (h0 : 0 ≤ t) : 0 ≤ St t := by
  rw [St]; positivity

/-- **Case A.**  If `a ≤ √2 - 1/2` and the pose is not high (`c ≤ a + D`), the line `y = a` cuts a
chord of length `≥ 1`, and both end node values are `≥ 1`. -/
lemma line_ge_one {C S c a : ℝ} (hC : 0 < C) (hS : 0 ≤ S) (h1 : C ^ 2 + S ^ 2 = 1)
    (ha1 : 0 ≤ a + 1 / 2) (ha2 : (a + 1 / 2) ^ 2 ≤ 2)
    (hp : (C + S) / 2 ≤ c) (hH : c ≤ a + ((C + S) / 2 - C * S)) :
    1 ≤ hL C S (a - c) + hR C S (a - c) := by
  rcases eq_or_lt_of_le hS with hS0 | hSp
  · rw [hL, hR, if_pos hS0.symm, if_pos hS0.symm]; norm_num
  rw [hL, hR, if_neg hSp.ne', if_neg hSp.ne']
  set e := a - c
  -- the lower side of the band: `a - c ≤ a - p ≤ D`, from `a ≤ √2 - 1/2`
  have hu1 : 1 ≤ C + S := by nlinarith
  have hu2 : (C + S) ^ 2 ≤ 2 := by nlinarith [sq_nonneg (C - S)]
  have hlow : e ≤ (C + S) / 2 - C * S := by
    -- `a ≤ u - (u² - 1)/2` for `u ∈ [1, √2]`, given `(a + 1/2)² ≤ 2`
    have hCS : 2 * (C * S) = (C + S) ^ 2 - 1 := by nlinarith
    have hr := Real.sqrt_nonneg 2
    have hr2 : Real.sqrt 2 ^ 2 = 2 := Real.sq_sqrt (by norm_num)
    have har : a + 1 / 2 ≤ Real.sqrt 2 := by
      rw [← Real.sqrt_sq ha1]; exact Real.sqrt_le_sqrt ha2
    have hur : C + S ≤ Real.sqrt 2 := by nlinarith
    have hprod : 0 ≤ (Real.sqrt 2 - (C + S)) * ((C + S) - (2 - Real.sqrt 2)) :=
      mul_nonneg (by linarith) (by nlinarith)
    nlinarith
  have hup : -e ≤ (C + S) / 2 - C * S := by simp only [e]; linarith
  have hCS0 : 0 < C * S := mul_pos hC hSp
  -- four combinations
  have k11 : 1 ≤ (1 / 2 + e * S) / C + (1 / 2 - e * S) / C := by
    rw [← add_div, le_div_iff₀ hC]; nlinarith
  have k22 : 1 ≤ (1 / 2 - e * C) / S + (e * C + 1 / 2) / S := by
    rw [← add_div, le_div_iff₀ hSp]; nlinarith
  have k12 : 1 ≤ (1 / 2 + e * S) / C + (e * C + 1 / 2) / S := by
    rw [div_add_div _ _ hC.ne' hSp.ne', le_div_iff₀ (by positivity)]; nlinarith
  have k21 : 1 ≤ (1 / 2 - e * C) / S + (1 / 2 - e * S) / C := by
    rw [div_add_div _ _ hSp.ne' hC.ne', le_div_iff₀ (by positivity)]; nlinarith
  rcases min_choice ((1 / 2 + e * S) / C) ((1 / 2 - e * C) / S) with h | h <;> rw [h] <;>
    rcases min_choice ((1 / 2 - e * S) / C) ((e * C + 1 / 2) / S) with h' | h' <;> rw [h'] <;>
    linarith

/-! ## 5. The wall-strip bound from certificate facts -/

/-- **The wall-strip bound from a certificate.**  Abstract certificate data: side heights
`zOf` (side `0` = END), a rule `R b b'` (the separation side when box `b` is immediately left of
box `b'`), H-flags, and box membership `inB b t c`.  Facts:
* the ordering line `y = a` satisfies `h - 1/2 < a < 1` and `a ≤ √2 - 1/2`;
* every rule height lies in `(h - 1/2, 1)`;
* `hcover`: every pose lies in some box;
* `hflag`: in an unflagged box every pose has `c ≤ a + D(t)` (is not H);
* `hchain`: for every four boxes, one of them flagged, there are numbers `v k` summing to `≥ 4`, each
  a lower bound of the node value over its box with the separations the rule prescribes.
Then at most three squares fit in the strip of height `h`. -/
theorem wallStrip_of_cert {h a : ℝ} (m : ℕ) (zOf : ℕ → ℝ) (R : ℕ → ℕ → ℕ) (flag : ℕ → Bool)
    (inB : ℕ → ℝ → ℝ → Prop)
    (ha1 : 0 ≤ a + 1 / 2) (ha2 : (a + 1 / 2) ^ 2 ≤ 2) (ha3 : h - 1 / 2 < a) (ha4 : a < 1)
    (hR0 : ∀ b b', b < m → b' < m → R b b' ≠ 0)
    (hz : ∀ b b', b < m → b' < m → h - 1 / 2 < zOf (R b b') ∧ zOf (R b b') < 1)
    (hcover : ∀ t c, Pose h t c → ∃ b, b < m ∧ inB b t c)
    (hflag : ∀ b t c, b < m → Pose h t c → inB b t c → flag b = false → c ≤ a + DT t)
    (hchain : ∀ b : Fin 4 → ℕ, (∀ k, b k < m) → (∃ k, flag (b k) = true) → ∃ v : Fin 4 → ℝ,
      4 ≤ v 0 + v 1 + v 2 + v 3 ∧ ∀ t c, Pose h t c →
        (inB (b 0) t c → v 0 ≤ nodeT zOf t c 0 (R (b 0) (b 1))) ∧
        (inB (b 1) t c → v 1 ≤ nodeT zOf t c (R (b 0) (b 1)) (R (b 1) (b 2))) ∧
        (inB (b 2) t c → v 2 ≤ nodeT zOf t c (R (b 1) (b 2)) (R (b 2) (b 3))) ∧
        (inB (b 3) t c → v 3 ≤ nodeT zOf t c (R (b 2) (b 3)) 0)) :
    WallStripBound h := by
  intro N T hT ctr ang hin hs hd
  by_contra hN
  push Not at hN
  choose C S hC hS h1 hsq using fun i => exists_CS (ang i)
  have hin' : ∀ i, sqCS (ctr i).1 (ctr i).2 (C i) (S i) ⊆ box T := fun i => by
    rw [← hsq i (ctr i)]; exact hin i
  have hd' : ∀ i j, i ≠ j → Disjoint (sqCS (ctr i).1 (ctr i).2 (C i) (S i))
      (sqCS (ctr j).1 (ctr j).2 (C j) (S j)) := fun i j hij => by
    rw [← hsq i (ctr i), ← hsq j (ctr j)]; exact hd i j hij
  have hbox := fun i => sqCS_box (h1 i) (hin' i)
  -- the tan-half-angle of each square
  set tt : Fin N → ℝ := fun i => S i / (1 + C i) with htt
  have htan := fun i => tan_half (hC i) (hS i) (h1 i)
  have hpT : ∀ i, pT (tt i) = (C i + S i) / 2 := fun i => by
    simp only [pT, htt, (htan i).2.2.1, (htan i).2.2.2]
  have hpose : ∀ i, Pose h (tt i) (ctr i).2 := fun i =>
    ⟨(htan i).1, (htan i).2.1, by rw [hpT]; exact (hbox i).1, hs i⟩
  -- every line in `(h - 1/2, 1)` crosses every square
  have hcross : ∀ i (z : ℝ), h - 1 / 2 < z → z < 1 → |z - (ctr i).2| ≤ (C i + S i) / 2 := by
    intro i z hz1 hz2
    have hu : 1 ≤ C i + S i := by nlinarith [hC i, hS i, h1 i]
    have := (hbox i).1
    have := hs i
    rw [abs_le]; constructor <;> linarith
  -- four of the squares, sorted along the line `y = a`
  let e : Fin 4 → Fin N := fun k => ⟨k, by omega⟩
  have he : Function.Injective e := fun i j hij => by
    simp only [e, Fin.mk.injEq] at hij; exact Fin.ext hij
  let key : Fin 4 → ℝ := fun k =>
    (ctr (e k)).1 + hR (C (e k)) (S (e k)) (a - (ctr (e k)).2)
  let σ := Tuple.sort key
  have hmono : Monotone (key ∘ σ) := Tuple.monotone_sort key
  have hinj : Function.Injective (key ∘ σ) := by
    intro i k hik
    by_contra hne
    have hne' : e (σ i) ≠ e (σ k) := fun h => hne (σ.injective (he h))
    have m1 := mem_right (x := (ctr (e (σ i))).1) (c := (ctr (e (σ i))).2) (hC _) (hS _) (h1 _)
      (hcross (e (σ i)) a ha3 ha4)
    have m2 := mem_right (x := (ctr (e (σ k))).1) (c := (ctr (e (σ k))).2) (hC _) (hS _) (h1 _)
      (hcross (e (σ k)) a ha3 ha4)
    simp only [Function.comp, key] at hik
    rw [hik, add_sub_cancel] at m1
    rw [add_sub_cancel] at m2
    exact Set.disjoint_left.mp (hd' _ _ hne') m1 m2
  have hsm : StrictMono (key ∘ σ) := hmono.strictMono_of_injective hinj
  let j : Fin 4 → Fin N := fun k => e (σ k)
  have hj : Function.Injective j := he.comp σ.injective
  set X : Fin 4 → ℝ := fun k => (ctr (j k)).1
  set Y : Fin 4 → ℝ := fun k => (ctr (j k)).2
  set CC : Fin 4 → ℝ := fun k => C (j k)
  set SS : Fin 4 → ℝ := fun k => S (j k)
  have hdj : ∀ i k, i ≠ k → Disjoint (sqCS (X i) (Y i) (CC i) (SS i))
      (sqCS (X k) (Y k) (CC k) (SS k)) := fun i k hik => hd' _ _ (fun h => hik (hj h))
  have hcr : ∀ z, h - 1 / 2 < z → z < 1 → Crosses Y CC SS z := fun z h1 h2 k =>
    hcross (j k) z h1 h2
  have o01 : key (σ 0) < key (σ 1) := hsm (by decide)
  have o12 : key (σ 1) < key (σ 2) := hsm (by decide)
  have o23 : key (σ 2) < key (σ 3) := hsm (by decide)
  have hx0 : (CC 0 + SS 0) / 2 ≤ X 0 := (hbox (j 0)).2.1
  have hx3 : X 3 + (CC 3 + SS 3) / 2 ≤ T := (hbox (j 3)).2.2
  have hchainlt := fun z0 z1 z2 (hz0 : Crosses Y CC SS z0) (hz1 : Crosses Y CC SS z1)
      (hz2 : Crosses Y CC SS z2) =>
    chain_lt X Y CC SS (a := a) (fun k => hC _) (fun k => hS _) (fun k => h1 _) hx0 hx3 hT hdj
      (hcr a ha3 ha4) hz0 hz1 hz2 o01 o12 o23
  -- node values in the `t`-parametrisation
  have hnode : ∀ k l r, nodeT zOf (tt (j k)) (Y k) l r =
      (if l = 0 then (CC k + SS k) / 2 else hL (CC k) (SS k) (zOf l - Y k)) +
        (if r = 0 then (CC k + SS k) / 2 else hR (CC k) (SS k) (zOf r - Y k)) := by
    intro k l r
    have e1 : Ct (tt (j k)) = CC k := (htan _).2.2.1
    have e2 : St (tt (j k)) = SS k := (htan _).2.2.2
    simp only [nodeT, hpT, e1, e2]
    rfl
  by_cases hH : ∃ k : Fin 4, a + DT (tt (j k)) < Y k
  · -- Case B: the certificate
    choose bx hbm hbx using fun k : Fin 4 => hcover _ _ (hpose (j k))
    have hfl : ∃ k, flag (bx k) = true := by
      obtain ⟨k, hk⟩ := hH
      refine ⟨k, ?_⟩
      by_contra hf
      have := hflag (bx k) _ _ (hbm k) (hpose (j k)) (hbx k) (by simpa using hf)
      exact absurd hk (not_lt.mpr this)
    obtain ⟨v, hv4, hvle⟩ := hchain bx hbm hfl
    have v0 := ((hvle _ _ (hpose (j 0))).1 (hbx 0))
    have v1 := ((hvle _ _ (hpose (j 1))).2.1 (hbx 1))
    have v2 := ((hvle _ _ (hpose (j 2))).2.2.1 (hbx 2))
    have v3 := ((hvle _ _ (hpose (j 3))).2.2.2 (hbx 3))
    rw [hnode] at v0 v1 v2 v3
    simp only [hR0 _ _ (hbm _) (hbm _), ↓reduceIte] at v0 v1 v2 v3
    have hz01 := hz _ _ (hbm 0) (hbm 1)
    have hz12 := hz _ _ (hbm 1) (hbm 2)
    have hz23 := hz _ _ (hbm 2) (hbm 3)
    have := hchainlt (zOf (R (bx 0) (bx 1))) (zOf (R (bx 1) (bx 2))) (zOf (R (bx 2) (bx 3)))
      (hcr _ hz01.1 hz01.2) (hcr _ hz12.1 hz12.2) (hcr _ hz23.1 hz23.2)
    linarith
  · -- Case A: every separation on `y = a`
    push Not at hH
    have hA : ∀ k, 1 ≤ hL (CC k) (SS k) (a - Y k) + hR (CC k) (SS k) (a - Y k) := by
      intro k
      have hk := hH k
      have e1 : Ct (tt (j k)) = CC k := (htan _).2.2.1
      have e2 : St (tt (j k)) = SS k := (htan _).2.2.2
      have e3 : pT (tt (j k)) = (CC k + SS k) / 2 := hpT _
      rw [DT, e1, e2, e3] at hk
      exact line_ge_one (hC _) (hS _) (h1 _) ha1 ha2 (hbox (j k)).1 (by linarith)
    have hl0 : hL (CC 0) (SS 0) (a - Y 0) ≤ (CC 0 + SS 0) / 2 := hL_le (hC _) (hS _) (h1 _)
    have hr3 : hR (CC 3) (SS 3) (a - Y 3) ≤ (CC 3 + SS 3) / 2 := hR_le (hC _) (hS _) (h1 _)
    have := hchainlt a a a (hcr a ha3 ha4) (hcr a ha3 ha4) (hcr a ha3 ha4)
    have := hA 0
    have := hA 1
    have := hA 2
    have := hA 3
    linarith

end StripChain
end SquarePacking
