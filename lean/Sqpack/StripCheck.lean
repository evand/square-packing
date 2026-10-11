import Sqpack.StripChain
import Sqpack.BernsteinZ

/-!
# An exact checker for wall-strip certificates, with its soundness proof

`tasks/n12-chat-lemmas/strip-demo/PROOF.md` §4–5.  A certificate (`Cert`) lists boxes
`[T0/TR, T1/TR] × [C0/K, C1/K]` in the pose plane (`t = tan(θ/2)`, centre height `c`), a menu of
separation heights `Z/K`, a rule (the separation side used when box `b` is immediately left of box
`b'`), leaf claims `(l, r, Q)` ("the node value with sides `l`, `r` is `≥ Q/K` on the box"; side `0`
is the wall), H-flags, a coverage witness and min-plus tables for the chain condition.  Everything is
decided in integer arithmetic (`decide +kernel` in `Sqpack/Strip117.lean`); nothing in the data is
trusted.

* Leaves (K3).  With `N = 1 + t²`, `M = 1 - t²`, `W = 2t` we have `cos θ = M/N`, `sin θ = W/N`, and
  each half-chord is a minimum of two of `L1 = (N/2 + eW)/M`, `L2 = (N/2 - eM)/W` (left) and
  `U1 = (N/2 - eW)/M`, `U2 = (eM + N/2)/W` (right).  So a leaf claim is a conjunction of at most four
  "combination" inequalities; each, times a positive denominator, is `A(t) + K c · B(t) ≥ 0` with
  integer polynomials `A`, `B` (`A11` … `B2E`).  Linear in `c`, it holds on the box if it holds at
  `c = C1/K` and at a lower bound of `c`: either `C0/K` or the half-height `p(t)` (every pose has
  `c ≥ p(t)`; then the test is a degree-6 polynomial).  Each polynomial is checked on `[T0, T1]` by
  `BernZ.check` with bisection (`chkSub`).  `leafOK_sound`.
* Coverage (K2): `covAll_sound`.  Flags (K4): `flagOK_sound` (an unflagged box has
  `c ≤ a + D(t)` on it, a Bernstein test of a quartic).  Chain (K5): `chainOK_sound`, from three
  min-plus tables checked transition by transition.
* `Cert.wallStrip`: all checks pass ⟹ `WallStripBound (hN/K)`.
-/

namespace SquarePacking
namespace StripCheck

open BernZ StripChain

/-- Scale of heights and claims. -/
def K : ℕ := 1000000
/-- Scale of `t`. -/
def TR : ℕ := 80000

lemma K_pos : (0 : ℝ) < K := by norm_num [K]
lemma TR_pos : 0 < TR := by norm_num [TR]

/-! ## 1. Bernstein test with bisection -/

/-- `p ≥ 0` on `[U0/R, U1/R]`: Bernstein test, bisecting up to depth `d`. -/
def chkSub (p : Poly) (n : ℕ) : ℕ → ℕ → ℕ → ℕ → Bool
  | U0, U1, R, 0 => check p n U0 U1 R
  | U0, U1, R, d + 1 => check p n U0 U1 R ||
      (chkSub p n (2 * U0) (U0 + U1) (2 * R) d && chkSub p n (U0 + U1) (2 * U1) (2 * R) d)

theorem chkSub_sound (p : Poly) (n : ℕ) : ∀ (d U0 U1 R : ℕ), chkSub p n U0 U1 R d = true →
    0 < R → ∀ u : ℝ, (U0 : ℝ) / R ≤ u → u ≤ (U1 : ℝ) / R → 0 ≤ eval p u
  | 0, U0, U1, R, h, hR, u, h0, h1 => nonneg_of_check h hR h0 h1
  | d + 1, U0, U1, R, h, hR, u, h0, h1 => by
    simp only [chkSub, Bool.or_eq_true, Bool.and_eq_true] at h
    rcases h with h | ⟨ha, hb⟩
    · exact nonneg_of_check h hR h0 h1
    have hR' : (0 : ℝ) < R := by exact_mod_cast hR
    have hR2 : 0 < 2 * R := by omega
    rcases le_total u (((U0 : ℝ) + U1) / (2 * R)) with hu | hu
    · refine chkSub_sound p n d _ _ _ ha hR2 u ?_ ?_
      · push_cast; rw [mul_div_mul_left _ _ two_ne_zero]; exact h0
      · push_cast; exact hu
    · refine chkSub_sound p n d _ _ _ hb hR2 u ?_ ?_
      · push_cast; exact hu
      · push_cast; rw [mul_div_mul_left _ _ two_ne_zero]; exact h1

/-- Maximal bisection depth. -/
def depth : ℕ := 6

/-- `p ≥ 0` on `[T0/TR, T1/TR]`. -/
def chk (p : Poly) (T0 T1 : ℕ) : Bool := chkSub p (p.length - 1) T0 T1 TR depth

lemma chk_sound {p : Poly} {T0 T1 : ℕ} (h : chk p T0 T1 = true) {u : ℝ}
    (h0 : (T0 : ℝ) / TR ≤ u) (h1 : u ≤ (T1 : ℝ) / TR) : 0 ≤ eval p u :=
  chkSub_sound p _ depth T0 T1 TR h TR_pos u h0 h1

/-! ## 2. The basic polynomials -/

/-- `N = 1 + t²`. -/
def pN : Poly := [1, 0, 1]
/-- `M = 1 - t²`. -/
def pM : Poly := [1, 0, -1]
/-- `W = 2t`. -/
def pW : Poly := [0, 2]
/-- `M + W`. -/
def pMW : Poly := [1, 2, -1]

@[simp] lemma eval_pN (t : ℝ) : eval pN t = 1 + t ^ 2 := by simp [pN]; ring
@[simp] lemma eval_pM (t : ℝ) : eval pM t = 1 - t ^ 2 := by simp [pM]; ring
@[simp] lemma eval_pW (t : ℝ) : eval pW t = 2 * t := by simp [pW]; ring
@[simp] lemma eval_pMW (t : ℝ) : eval pMW t = 1 + 2 * t - t ^ 2 := by simp [pMW]; ring

lemma pT_eq (t : ℝ) : pT t = (1 + 2 * t - t ^ 2) / (2 * (1 + t ^ 2)) := by
  have : (1 + t ^ 2) ≠ 0 := by positivity
  rw [pT, Ct, St]; field_simp; ring

/-! ## 3. One combination inequality: linear in `c` -/

/-- `A(t) + K c B(t) ≥ 0` on the box (with `c ≥ p(t)`), from Bernstein tests. -/
def comboOK (A B : Poly) (T0 T1 C0 C1 : ℕ) : Bool :=
  match B with
  | [] => chk A T0 T1
  | b :: bs => chk (add A (smul C1 (b :: bs))) T0 T1 &&
      (chk (add A (smul C0 (b :: bs))) T0 T1 ||
        chk (add (mul [2, 0, 2] A) (mul (smul K pMW) (b :: bs))) T0 T1)

lemma comboOK_sound {A B : Poly} {T0 T1 C0 C1 : ℕ} (h : comboOK A B T0 T1 C0 C1 = true)
    {t c : ℝ} (ht0 : (T0 : ℝ) / TR ≤ t) (ht1 : t ≤ (T1 : ℝ) / TR) (hc0 : (C0 : ℝ) / K ≤ c)
    (hp : pT t ≤ c) (hc1 : c ≤ (C1 : ℝ) / K) : 0 ≤ eval A t + (K * c) * eval B t := by
  have hK := K_pos
  match B, h with
  | [], h =>
    simp only [comboOK] at h
    simpa using chk_sound h ht0 ht1
  | b :: bs, h =>
    simp only [comboOK, Bool.and_eq_true, Bool.or_eq_true] at h
    obtain ⟨hhi, hlo⟩ := h
    have fhi := chk_sound hhi ht0 ht1
    simp only [eval_add, eval_smul] at fhi
    push_cast at fhi
    have hx1 : K * c ≤ C1 := by rw [le_div_iff₀ hK] at hc1; linarith
    rcases le_total 0 (eval (b :: bs) t) with hB | hB
    · rcases hlo with hl | hl
      · have flo := chk_sound hl ht0 ht1
        simp only [eval_add, eval_smul] at flo
        push_cast at flo
        have hx0 : (C0 : ℝ) ≤ K * c := by rw [div_le_iff₀ hK] at hc0; linarith
        nlinarith [mul_le_mul_of_nonneg_right hx0 hB]
      · have fp := chk_sound hl ht0 ht1
        simp only [eval_add, eval_mul, eval_smul, eval_pMW] at fp
        have e2 : eval [2, 0, 2] t = 2 * (1 + t ^ 2) := by simp; ring
        rw [e2] at fp
        rw [pT_eq, div_le_iff₀ (by positivity)] at hp
        have hN : 0 < 2 * (1 + t ^ 2) := by positivity
        -- `2N (A + K c B) ≥ 2N A + K (M + W) B ≥ 0`
        have : 0 ≤ 2 * (1 + t ^ 2) * (eval A t + K * c * eval (b :: bs) t) := by
          have hmono : (K : ℝ) * (1 + 2 * t - t ^ 2) * eval (b :: bs) t ≤
              K * (c * (2 * (1 + t ^ 2))) * eval (b :: bs) t :=
            mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hp hK.le) hB
          push_cast at fp
          nlinarith
        exact (mul_nonneg_iff_of_pos_left hN).mp this
    · nlinarith [mul_le_mul_of_nonneg_right hx1 (neg_nonneg.mpr hB)]

/-! ## 4. Half-chords in the `t`-parametrisation -/

/-- `L1 = (N/2 + eW)/M`. -/
noncomputable def L1t (t e : ℝ) : ℝ := ((1 + t ^ 2) / 2 + e * (2 * t)) / (1 - t ^ 2)
/-- `L2 = (N/2 - eM)/W`. -/
noncomputable def L2t (t e : ℝ) : ℝ := ((1 + t ^ 2) / 2 - e * (1 - t ^ 2)) / (2 * t)
/-- `U1 = (N/2 - eW)/M`. -/
noncomputable def U1t (t e : ℝ) : ℝ := ((1 + t ^ 2) / 2 - e * (2 * t)) / (1 - t ^ 2)
/-- `U2 = (eM + N/2)/W`. -/
noncomputable def U2t (t e : ℝ) : ℝ := (e * (1 - t ^ 2) + (1 + t ^ 2) / 2) / (2 * t)

lemma St_ne_zero {t : ℝ} (ht : 0 < t) : St t ≠ 0 := by
  rw [St]; positivity

lemma hR_t {t e : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    hR (Ct t) (St t) e = min (U1t t e) (U2t t e) := by
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  have hM : (1 - t ^ 2) ≠ 0 := by nlinarith
  have hW : t ≠ 0 := ht0.ne'
  rw [hR, if_neg (St_ne_zero ht0)]
  congr 1
  · rw [U1t, Ct, St]; field_simp
  · rw [U2t, Ct, St]; field_simp

lemma hL_t {t e : ℝ} (ht0 : 0 < t) (ht1 : t < 1) :
    hL (Ct t) (St t) e = min (L1t t e) (L2t t e) := by
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  have hM : (1 - t ^ 2) ≠ 0 := by nlinarith
  have hW : t ≠ 0 := ht0.ne'
  rw [hL, if_neg (St_ne_zero ht0)]
  congr 1
  · rw [L1t, Ct, St]; field_simp
  · rw [L2t, Ct, St]; field_simp

lemma hR_t0 (e : ℝ) : hR (Ct 0) (St 0) e = U1t 0 e := by
  rw [hR, if_pos (by simp [St]), U1t]; norm_num

lemma hL_t0 (e : ℝ) : hL (Ct 0) (St 0) e = L1t 0 e := by
  rw [hL, if_pos (by simp [St]), L1t]; norm_num

/-! ## 5. The eight combination polynomials -/

/-- `K M (L1 + U1 - q)` (independent of `c`). -/
def A11 (Zl Zr Q : ℤ) : Poly := add (smul K pN) (add (smul (Zl - Zr) pW) (smul (-Q) pM))
/-- `K W (L2 + U2 - q)` (independent of `c`). -/
def A22 (Zl Zr Q : ℤ) : Poly := add (smul K pN) (add (smul (-(Zl - Zr)) pM) (smul (-Q) pW))
/-- `2 K M W (L1 + U2 - q)`. -/
def A12 (Zl Zr Q : ℤ) : Poly :=
  add (mul pW (add (smul K pN) (smul (2 * Zl) pW)))
    (add (mul pM (add (smul (2 * Zr) pM) (smul K pN))) (smul (-2 * Q) (mul pM pW)))
def B12 : Poly := smul (-2) (add (mul pW pW) (mul pM pM))
/-- `2 K M W (L2 + U1 - q)`. -/
def A21 (Zl Zr Q : ℤ) : Poly :=
  add (mul pM (add (smul K pN) (smul (-2 * Zl) pM)))
    (add (mul pW (add (smul K pN) (smul (-2 * Zr) pW))) (smul (-2 * Q) (mul pM pW)))
def B21 : Poly := smul 2 (add (mul pM pM) (mul pW pW))
/-- `2 K N M (p + U1 - q)`. -/
def AE1 (Zr Q : ℤ) : Poly :=
  add (smul K (mul pMW pM)) (add (mul pN (add (smul K pN) (smul (-2 * Zr) pW)))
    (smul (-2 * Q) (mul pN pM)))
def BE1 : Poly := smul 2 (mul pN pW)
/-- `2 K N W (p + U2 - q)`. -/
def AE2 (Zr Q : ℤ) : Poly :=
  add (smul K (mul pMW pW)) (add (mul pN (add (smul (2 * Zr) pM) (smul K pN)))
    (smul (-2 * Q) (mul pN pW)))
def BE2 : Poly := smul (-2) (mul pN pM)
/-- `2 K N M (L1 + p - q)`. -/
def A1E (Zl Q : ℤ) : Poly :=
  add (mul pN (add (smul K pN) (smul (2 * Zl) pW))) (add (smul K (mul pMW pM))
    (smul (-2 * Q) (mul pN pM)))
def B1E : Poly := smul (-2) (mul pN pW)
/-- `2 K N W (L2 + p - q)`. -/
def A2E (Zl Q : ℤ) : Poly :=
  add (mul pN (add (smul K pN) (smul (-2 * Zl) pM))) (add (smul K (mul pMW pW))
    (smul (-2 * Q) (mul pN pW)))
def B2E : Poly := smul 2 (mul pN pM)

section ids
variable (t c : ℝ) (Zl Zr Q : ℤ)

local notation "el" => (Zl : ℝ) / K - c
local notation "er" => (Zr : ℝ) / K - c
local notation "q" => (Q : ℝ) / K

lemma K_ne : (K : ℝ) ≠ 0 := K_pos.ne'

lemma id11 (hM : 1 - t ^ 2 ≠ 0) :
    eval (A11 Zl Zr Q) t + (K * c) * eval [] t = K * (1 - t ^ 2) * (L1t t el + U1t t er - q) := by
  have := K_ne
  simp only [A11, eval_add, eval_smul, eval_pN, eval_pM, eval_pW, eval_nil, L1t, U1t]
  push_cast; field_simp; ring

lemma id22 (hW : t ≠ 0) :
    eval (A22 Zl Zr Q) t + (K * c) * eval [] t = K * (2 * t) * (L2t t el + U2t t er - q) := by
  have := K_ne
  simp only [A22, eval_add, eval_smul, eval_pN, eval_pM, eval_pW, eval_nil, L2t, U2t]
  push_cast; field_simp; ring

lemma id12 (hM : 1 - t ^ 2 ≠ 0) (hW : t ≠ 0) :
    eval (A12 Zl Zr Q) t + (K * c) * eval B12 t =
      2 * K * (1 - t ^ 2) * (2 * t) * (L1t t el + U2t t er - q) := by
  have := K_ne
  simp only [A12, B12, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, L1t, U2t]
  push_cast; field_simp; ring

lemma id21 (hM : 1 - t ^ 2 ≠ 0) (hW : t ≠ 0) :
    eval (A21 Zl Zr Q) t + (K * c) * eval B21 t =
      2 * K * (1 - t ^ 2) * (2 * t) * (L2t t el + U1t t er - q) := by
  have := K_ne
  simp only [A21, B21, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, L2t, U1t]
  push_cast; field_simp; ring

lemma idE1 (hM : 1 - t ^ 2 ≠ 0) :
    eval (AE1 Zr Q) t + (K * c) * eval BE1 t =
      2 * K * (1 + t ^ 2) * (1 - t ^ 2) * (pT t + U1t t er - q) := by
  have := K_ne
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  simp only [AE1, BE1, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, eval_pMW, U1t,
    pT_eq]
  push_cast; field_simp; ring

lemma idE2 (hW : t ≠ 0) :
    eval (AE2 Zr Q) t + (K * c) * eval BE2 t =
      2 * K * (1 + t ^ 2) * (2 * t) * (pT t + U2t t er - q) := by
  have := K_ne
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  simp only [AE2, BE2, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, eval_pMW, U2t,
    pT_eq]
  push_cast; field_simp; ring

lemma id1E (hM : 1 - t ^ 2 ≠ 0) :
    eval (A1E Zl Q) t + (K * c) * eval B1E t =
      2 * K * (1 + t ^ 2) * (1 - t ^ 2) * (L1t t el + pT t - q) := by
  have := K_ne
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  simp only [A1E, B1E, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, eval_pMW, L1t,
    pT_eq]
  push_cast; field_simp; ring

lemma id2E (hW : t ≠ 0) :
    eval (A2E Zl Q) t + (K * c) * eval B2E t =
      2 * K * (1 + t ^ 2) * (2 * t) * (L2t t el + pT t - q) := by
  have := K_ne
  have hN : (1 + t ^ 2) ≠ 0 := by positivity
  simp only [A2E, B2E, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, eval_pMW, L2t,
    pT_eq]
  push_cast; field_simp; ring

end ids

/-! ## 6. Leaves -/

/-- Height numerator (over `K`) of side `j ≥ 1`. -/
def zval (menu : List ℕ) (j : ℕ) : ℕ := menu.getD (j - 1) 0

/-- Height of side `j ≥ 1`. -/
noncomputable def zOf (menu : List ℕ) (j : ℕ) : ℝ := (zval menu j : ℝ) / K

/-- The leaf claim "node value with sides `l`, `r` is `≥ Q/K` on the box `bx`". -/
def leafOK (menu : List ℕ) (bx : ℕ × ℕ × ℕ × ℕ) (l r Q : ℕ) : Bool :=
  let T0 := bx.1
  let T1 := bx.2.1
  let C0 := bx.2.2.1
  let C1 := bx.2.2.2
  let Zl : ℤ := zval menu l
  let Zr : ℤ := zval menu r
  if l = 0 then
    if r = 0 then false
    else comboOK (AE1 Zr Q) BE1 T0 T1 C0 C1 && comboOK (AE2 Zr Q) BE2 T0 T1 C0 C1
  else if r = 0 then
    comboOK (A1E Zl Q) B1E T0 T1 C0 C1 && comboOK (A2E Zl Q) B2E T0 T1 C0 C1
  else
    comboOK (A11 Zl Zr Q) [] T0 T1 C0 C1 && comboOK (A22 Zl Zr Q) [] T0 T1 C0 C1 &&
      comboOK (A12 Zl Zr Q) B12 T0 T1 C0 C1 && comboOK (A21 Zl Zr Q) B21 T0 T1 C0 C1

lemma le_of_poly {X q f D : ℝ} (hD : 0 < D) (hf : 0 ≤ f) (h : f = D * (X - q)) : q ≤ X := by
  rw [h] at hf
  linarith [(mul_nonneg_iff_of_pos_left hD).mp hf]

lemma le_min_add_min {q x1 x2 y1 y2 : ℝ} (h11 : q ≤ x1 + y1) (h12 : q ≤ x1 + y2)
    (h21 : q ≤ x2 + y1) (h22 : q ≤ x2 + y2) : q ≤ min x1 x2 + min y1 y2 := by
  rcases min_choice x1 x2 with h | h <;> rw [h] <;>
    rcases min_choice y1 y2 with h' | h' <;> rw [h'] <;> assumption

/-- **Leaf soundness.** -/
theorem leafOK_sound {menu : List ℕ} {bx : ℕ × ℕ × ℕ × ℕ} {l r Q : ℕ}
    (h : leafOK menu bx l r Q = true) {t c : ℝ} (ht0 : 0 ≤ t) (ht1 : t < 1)
    (hb0 : (bx.1 : ℝ) / TR ≤ t) (hb1 : t ≤ (bx.2.1 : ℝ) / TR)
    (hc0 : (bx.2.2.1 : ℝ) / K ≤ c) (hc1 : c ≤ (bx.2.2.2 : ℝ) / K) (hp : pT t ≤ c) :
    (Q : ℝ) / K ≤ nodeT (zOf menu) t c l r := by
  have hK := K_pos
  have hM : 0 < 1 - t ^ 2 := by nlinarith
  have hN : 0 < 1 + t ^ 2 := by positivity
  have hzl : zOf menu l - c = (((zval menu l : ℤ) : ℝ)) / K - c := by simp [zOf]
  have hzr : zOf menu r - c = (((zval menu r : ℤ) : ℝ)) / K - c := by simp [zOf]
  have hQ : (Q : ℝ) / K = ((Q : ℤ) : ℝ) / K := by simp
  have cb := fun {A B : Poly} (hh : comboOK A B bx.1 bx.2.1 bx.2.2.1 bx.2.2.2 = true) =>
    comboOK_sound hh hb0 hb1 hc0 hp hc1
  rw [nodeT, hQ]
  by_cases hl : l = 0 <;> by_cases hr : r = 0
  · simp [leafOK, hl, hr] at h
  · -- left end
    simp only [leafOK, hl, hr, ↓reduceIte, Bool.and_eq_true] at h
    rw [if_pos hl, if_neg hr, hzr]
    have k1 := le_of_poly (by positivity) (cb h.1) (idE1 t c _ _ hM.ne')
    rcases eq_or_lt_of_le ht0 with h0 | h0
    · subst h0; rw [hR_t0]; exact k1
    have k2 := le_of_poly (by positivity) (cb h.2) (idE2 t c _ _ h0.ne')
    rw [hR_t h0 ht1]
    rcases min_choice (U1t t ((zval menu r : ℤ) / K - c)) (U2t t ((zval menu r : ℤ) / K - c))
      with e | e <;> rw [e] <;> assumption
  · -- right end
    simp only [leafOK, hl, hr, ↓reduceIte, Bool.and_eq_true] at h
    rw [if_neg hl, if_pos hr, hzl]
    have k1 := le_of_poly (by positivity) (cb h.1) (id1E t c _ _ hM.ne')
    rcases eq_or_lt_of_le ht0 with h0 | h0
    · subst h0; rw [hL_t0]; exact k1
    have k2 := le_of_poly (by positivity) (cb h.2) (id2E t c _ _ h0.ne')
    rw [hL_t h0 ht1]
    rcases min_choice (L1t t ((zval menu l : ℤ) / K - c)) (L2t t ((zval menu l : ℤ) / K - c))
      with e | e <;> rw [e] <;> assumption
  · -- two lines
    simp only [leafOK, hl, hr, ↓reduceIte, Bool.and_eq_true] at h
    obtain ⟨⟨⟨h11, h22⟩, h12⟩, h21⟩ := h
    rw [if_neg hl, if_neg hr, hzl, hzr]
    have k11 := le_of_poly (by positivity) (cb h11) (id11 t c _ _ _ hM.ne')
    rcases eq_or_lt_of_le ht0 with h0 | h0
    · subst h0; rw [hL_t0, hR_t0]; exact k11
    have k22 := le_of_poly (by positivity) (cb h22) (id22 t c _ _ _ h0.ne')
    have k12 := le_of_poly (by positivity) (cb h12) (id12 t c _ _ _ hM.ne' h0.ne')
    have k21 := le_of_poly (by positivity) (cb h21) (id21 t c _ _ _ hM.ne' h0.ne')
    rw [hL_t h0 ht1, hR_t h0 ht1]
    exact le_min_add_min k11 k12 k21 k22

/-! ## 7. Certificates -/

/-- A wall-strip certificate.  Heights and claims are over `K`, `t`-values over `TR`. -/
structure Cert where
  /-- strip height `h = hN/K` -/
  hN : ℕ
  /-- ordering line `a = aN/K` -/
  aN : ℕ
  /-- separation heights; side `j ≥ 1` is `menu[j-1]/K`, side `0` the wall -/
  menu : List ℕ
  /-- boxes `(T0, T1, C0, C1)`: `[T0/TR, T1/TR] × [C0/K, C1/K]` -/
  boxes : List (ℕ × ℕ × ℕ × ℕ)
  /-- `rule[b][b']`: the separation side when box `b` is immediately left of box `b'` -/
  rule : List (List ℕ)
  /-- H-flags -/
  flags : List Bool
  /-- leaf claims `(l, r, Q)` per box -/
  claims : List (List (ℕ × ℕ × ℕ))
  /-- coverage: `t`-strips `[S0/TR, S1/TR]` and a chain of boxes covering `[p(t), h]` on each -/
  cover : List (ℕ × ℕ × List ℕ)
  /-- min-plus tables: `phi_k[b]` lists `(z, f, s)`, a lower bound `s` of the partial chain sum
  over the first `k` boxes when box `b` comes next with left side `z` and flag-so-far `f` -/
  phi1 : List (List (ℕ × Bool × ℕ))
  phi2 : List (List (ℕ × Bool × ℕ))
  phi3 : List (List (ℕ × Bool × ℕ))

namespace Cert

variable (C : Cert)

def m : ℕ := C.boxes.length
def box (b : ℕ) : ℕ × ℕ × ℕ × ℕ := C.boxes.getD b (0, 0, 0, 0)
def R (b b' : ℕ) : ℕ := (C.rule.getD b []).getD b' 0
def flag (b : ℕ) : Bool := C.flags.getD b false

/-- Pose `(t, c)` lies in box `b`. -/
def inB (b : ℕ) (t c : ℝ) : Prop :=
  ((C.box b).1 : ℝ) / TR ≤ t ∧ t ≤ ((C.box b).2.1 : ℝ) / TR ∧
    ((C.box b).2.2.1 : ℝ) / K ≤ c ∧ c ≤ ((C.box b).2.2.2 : ℝ) / K

def claim (b l r : ℕ) : Option ℕ :=
  ((C.claims.getD b []).find? (fun e => e.1 == l && e.2.1 == r)).map (fun e => e.2.2)

/-! ### Parameters and rule (K1) -/

def paramsOK : Bool :=
  decide ((2 * C.aN + K) ^ 2 ≤ 8 * K ^ 2) && decide (2 * C.hN < 2 * C.aN + K) && decide (C.aN < K)

def ruleOK : Bool :=
  (List.range C.m).all fun b => (List.range C.m).all fun b' =>
    decide (C.R b b' ≠ 0) && decide (2 * C.hN < 2 * zval C.menu (C.R b b') + K) &&
      decide (zval C.menu (C.R b b') < K)

/-! ### Coverage (K2) -/

/-- `C0/K ≤ p(t)` on `[S0/TR, S1/TR]`. -/
def pLow (C0 S0 S1 : ℕ) : Bool :=
  decide (2 * C0 ≤ K) || chk (add (smul K pMW) (smul (-2 * (C0 : ℤ)) pN)) S0 S1

def covChain (S0 S1 : ℕ) : ℕ → List ℕ → Bool
  | _, [] => false
  | cur, b :: bs => decide (b < C.m) && decide ((C.box b).1 ≤ S0) && decide (S1 ≤ (C.box b).2.1) &&
      decide ((C.box b).2.2.1 ≤ cur) &&
      (decide (C.hN ≤ (C.box b).2.2.2) || covChain S0 S1 (C.box b).2.2.2 bs)

def covEntry (S0 S1 : ℕ) : List ℕ → Bool
  | [] => false
  | b :: bs => decide (b < C.m) && decide ((C.box b).1 ≤ S0) && decide (S1 ≤ (C.box b).2.1) &&
      pLow (C.box b).2.2.1 S0 S1 &&
      (decide (C.hN ≤ (C.box b).2.2.2) || covChain C S0 S1 (C.box b).2.2.2 bs)

def covAll : ℕ → List (ℕ × ℕ × List ℕ) → Bool
  | s, [] => decide (TR ≤ s)
  | s, (S0, S1, bs) :: rest => decide (S0 ≤ s) && covEntry C S0 S1 bs && covAll S1 rest

def coverOK : Bool := covAll C 0 C.cover

/-! ### Flags (K4) -/

/-- `2 K N² (a - C1/K + D(t))`. -/
def flagPoly (aN C1 : ℕ) : Poly :=
  add (smul (2 * ((aN : ℤ) - C1)) (mul pN pN))
    (add (smul K (mul pN pMW)) (smul (-2 * (K : ℤ)) (mul pM pW)))

def flagsOK : Bool :=
  (List.range C.m).all fun b =>
    C.flag b || chk (flagPoly C.aN (C.box b).2.2.2) (C.box b).1 (C.box b).2.1

/-! ### Leaves (K3) -/

def boxLeavesOK (b : ℕ) : Bool :=
  (C.claims.getD b []).all fun e => leafOK C.menu (C.box b) e.1 e.2.1 e.2.2

/-! ### Chain (K5) -/

def lookup (φ : List (List (ℕ × Bool × ℕ))) (b z : ℕ) (f : Bool) : Option ℕ :=
  ((φ.getD b []).find? (fun e => e.1 == z && e.2.1 == f)).map (fun e => e.2.2)

def chain1 : Bool :=
  (List.range C.m).all fun b1 => (List.range C.m).all fun b2 =>
    match C.claim b1 0 (C.R b1 b2), lookup C.phi1 b2 (C.R b1 b2) (C.flag b1) with
    | some v, some s => decide (s ≤ v)
    | _, _ => false

def chainStep (φ ψ : List (List (ℕ × Bool × ℕ))) : Bool :=
  (List.range C.m).all fun b => (φ.getD b []).all fun e => (List.range C.m).all fun b' =>
    match C.claim b e.1 (C.R b b'), lookup ψ b' (C.R b b') (e.2.1 || C.flag b) with
    | some v, some s' => decide (s' ≤ e.2.2 + v)
    | _, _ => false

def chainEnd (φ : List (List (ℕ × Bool × ℕ))) : Bool :=
  (List.range C.m).all fun b => (φ.getD b []).all fun e =>
    !(e.2.1 || C.flag b) ||
      match C.claim b e.1 0 with
      | some v => decide (4 * K ≤ e.2.2 + v)
      | none => false

def chainOK : Bool :=
  C.chain1 && C.chainStep C.phi1 C.phi2 && C.chainStep C.phi2 C.phi3 && C.chainEnd C.phi3

/-! ### Soundness -/

lemma pT_ge_half {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1) : 1 / 2 ≤ pT t := by
  rw [pT_eq, le_div_iff₀ (by positivity)]; nlinarith

lemma pLow_sound {C0 S0 S1 : ℕ} (h : pLow C0 S0 S1 = true) {t : ℝ} (h0 : 0 ≤ t) (h1 : t ≤ 1)
    (hs0 : (S0 : ℝ) / TR ≤ t) (hs1 : t ≤ (S1 : ℝ) / TR) : (C0 : ℝ) / K ≤ pT t := by
  have hK := K_pos
  simp only [pLow, Bool.or_eq_true, decide_eq_true_eq] at h
  rcases h with h | h
  · have : (C0 : ℝ) / K ≤ 1 / 2 := by
      rw [div_le_iff₀ hK]; have : (2 * C0 : ℝ) ≤ K := by exact_mod_cast h
      linarith
    linarith [pT_ge_half h0 h1]
  · have f := chk_sound h hs0 hs1
    simp only [eval_add, eval_smul, eval_pMW, eval_pN] at f
    push_cast at f
    rw [pT_eq, div_le_div_iff₀ hK (by positivity)]
    nlinarith

lemma covChain_sound {S0 S1 : ℕ} : ∀ (cur : ℕ) (bs : List ℕ), C.covChain S0 S1 cur bs = true →
    ∀ t c : ℝ, (S0 : ℝ) / TR ≤ t → t ≤ (S1 : ℝ) / TR → (cur : ℝ) / K ≤ c →
      c ≤ (C.hN : ℝ) / K → ∃ b, b < C.m ∧ C.inB b t c
  | _, [], h => by simp [covChain] at h
  | cur, b :: bs, h => by
    intro t c hs0 hs1 hc hch
    have hK := K_pos
    have hTR : (0 : ℝ) < TR := by exact_mod_cast TR_pos
    simp only [covChain, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨⟨⟨hb, ht0⟩, ht1⟩, hc0⟩, hnext⟩ := h
    have it0 : ((C.box b).1 : ℝ) / TR ≤ t :=
      le_trans (div_le_div_of_nonneg_right (by exact_mod_cast ht0) hTR.le) hs0
    have it1 : t ≤ ((C.box b).2.1 : ℝ) / TR :=
      le_trans hs1 (div_le_div_of_nonneg_right (by exact_mod_cast ht1) hTR.le)
    have ic0 : ((C.box b).2.2.1 : ℝ) / K ≤ c :=
      le_trans (div_le_div_of_nonneg_right (by exact_mod_cast hc0) hK.le) hc
    rcases le_or_gt c (((C.box b).2.2.2 : ℝ) / K) with hcc | hcc
    · exact ⟨b, hb, it0, it1, ic0, hcc⟩
    · rcases hnext with hn | hn
      · exfalso
        have : (C.hN : ℝ) / K ≤ ((C.box b).2.2.2 : ℝ) / K :=
          div_le_div_of_nonneg_right (by exact_mod_cast hn) hK.le
        linarith
      · exact covChain_sound _ bs hn t c hs0 hs1 hcc.le hch

lemma covEntry_sound {S0 S1 : ℕ} {bs : List ℕ} (h : C.covEntry S0 S1 bs = true) {t c : ℝ}
    (h0 : 0 ≤ t) (h1 : t ≤ 1) (hs0 : (S0 : ℝ) / TR ≤ t) (hs1 : t ≤ (S1 : ℝ) / TR)
    (hp : pT t ≤ c) (hch : c ≤ (C.hN : ℝ) / K) : ∃ b, b < C.m ∧ C.inB b t c := by
  have hK := K_pos
  have hTR : (0 : ℝ) < TR := by exact_mod_cast TR_pos
  match bs, h with
  | [], h => simp [covEntry] at h
  | b :: bs, h =>
    simp only [covEntry, Bool.and_eq_true, Bool.or_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨⟨⟨hb, ht0⟩, ht1⟩, hpl⟩, hnext⟩ := h
    have it0 : ((C.box b).1 : ℝ) / TR ≤ t :=
      le_trans (div_le_div_of_nonneg_right (by exact_mod_cast ht0) hTR.le) hs0
    have it1 : t ≤ ((C.box b).2.1 : ℝ) / TR :=
      le_trans hs1 (div_le_div_of_nonneg_right (by exact_mod_cast ht1) hTR.le)
    have ic0 : ((C.box b).2.2.1 : ℝ) / K ≤ c := (pLow_sound hpl h0 h1 hs0 hs1).trans hp
    rcases le_or_gt c (((C.box b).2.2.2 : ℝ) / K) with hcc | hcc
    · exact ⟨b, hb, it0, it1, ic0, hcc⟩
    · rcases hnext with hn | hn
      · exfalso
        have : (C.hN : ℝ) / K ≤ ((C.box b).2.2.2 : ℝ) / K :=
          div_le_div_of_nonneg_right (by exact_mod_cast hn) hK.le
        linarith
      · exact C.covChain_sound _ bs hn t c hs0 hs1 hcc.le hch

lemma covAll_sound : ∀ (s : ℕ) (L : List (ℕ × ℕ × List ℕ)), C.covAll s L = true →
    ∀ t c : ℝ, (s : ℝ) / TR ≤ t → t < 1 → pT t ≤ c → c ≤ (C.hN : ℝ) / K →
      ∃ b, b < C.m ∧ C.inB b t c
  | s, [], h => by
    intro t c hs ht1
    have hTR : (0 : ℝ) < TR := by exact_mod_cast TR_pos
    simp only [covAll, decide_eq_true_eq] at h
    have : (1 : ℝ) ≤ s / TR := by rw [le_div_iff₀ hTR]; simpa using (by exact_mod_cast h : (TR : ℝ) ≤ s)
    linarith
  | s, (S0, S1, bs) :: rest, h => by
    intro t c hs ht1 hp hch
    have hTR : (0 : ℝ) < TR := by exact_mod_cast TR_pos
    simp only [covAll, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨hS0, he⟩, hrest⟩ := h
    have h0 : 0 ≤ t := le_trans (by positivity) hs
    rcases le_or_gt t ((S1 : ℝ) / TR) with h1 | h1
    · exact C.covEntry_sound he h0 ht1.le
        (le_trans (div_le_div_of_nonneg_right (by exact_mod_cast hS0) hTR.le) hs) h1 hp hch
    · exact covAll_sound S1 rest hrest t c h1.le ht1 hp hch

lemma flag_sound (hf : C.flagsOK = true) {b : ℕ} (hb : b < C.m) (hfl : C.flag b = false)
    {t c : ℝ} (h1 : C.inB b t c) : c ≤ (C.aN : ℝ) / K + DT t := by
  have hK := K_pos
  simp only [flagsOK, List.all_eq_true, List.mem_range, Bool.or_eq_true] at hf
  rcases hf b hb with h | h
  · rw [hfl] at h; exact absurd h (by simp)
  obtain ⟨ht0, ht1, -, hc1⟩ := h1
  have f := chk_sound h ht0 ht1
  have hN : (0 : ℝ) < 1 + t ^ 2 := by positivity
  have e : eval (flagPoly C.aN (C.box b).2.2.2) t =
      2 * K * (1 + t ^ 2) ^ 2 * ((C.aN : ℝ) / K - ((C.box b).2.2.2 : ℝ) / K + DT t) := by
    simp only [flagPoly, eval_add, eval_smul, eval_mul, eval_pN, eval_pM, eval_pW, eval_pMW, DT,
      pT_eq, Ct, St]
    push_cast; field_simp; ring
  rw [e] at f
  have := (mul_nonneg_iff_of_pos_left (by positivity)).mp f
  linarith

lemma mem_of_lookup {φ : List (List (ℕ × Bool × ℕ))} {b z : ℕ} {f : Bool} {s : ℕ}
    (h : lookup φ b z f = some s) : (z, f, s) ∈ φ.getD b [] := by
  simp only [lookup, Option.map_eq_some_iff] at h
  obtain ⟨e, he, rfl⟩ := h
  have hp := List.find?_some he
  simp only [Bool.and_eq_true, beq_iff_eq] at hp
  have hm := List.mem_of_find?_eq_some he
  obtain ⟨e1, e2, e3⟩ := e
  simp only at hp
  rw [hp.1, hp.2] at hm
  exact hm

lemma mem_of_claim {b l r v : ℕ} (h : C.claim b l r = some v) : (l, r, v) ∈ C.claims.getD b [] := by
  simp only [claim, Option.map_eq_some_iff] at h
  obtain ⟨e, he, rfl⟩ := h
  have hp := List.find?_some he
  simp only [Bool.and_eq_true, beq_iff_eq] at hp
  have hm := List.mem_of_find?_eq_some he
  obtain ⟨e1, e2, e3⟩ := e
  simp only at hp
  rw [hp.1, hp.2] at hm
  exact hm

/-- A claim that is present is a valid lower bound on its box. -/
lemma claim_sound (hl : ∀ b < C.m, C.boxLeavesOK b = true) {b l r v : ℕ} (hb : b < C.m)
    (h : C.claim b l r = some v) {t c : ℝ} (ht0 : 0 ≤ t) (ht1 : t < 1) (hp : pT t ≤ c)
    (hin : C.inB b t c) : (v : ℝ) / K ≤ nodeT (zOf C.menu) t c l r := by
  have hbl := hl b hb
  simp only [boxLeavesOK, List.all_eq_true] at hbl
  have := hbl _ (C.mem_of_claim h)
  exact leafOK_sound this ht0 ht1 hin.1 hin.2.1 hin.2.2.1 hin.2.2.2 hp

lemma chain_sound (hc : C.chainOK = true) (b : Fin 4 → ℕ) (hb : ∀ k, b k < C.m)
    (hf : ∃ k, C.flag (b k) = true) :
    ∃ v0 v1 v2 v3 : ℕ, C.claim (b 0) 0 (C.R (b 0) (b 1)) = some v0 ∧
      C.claim (b 1) (C.R (b 0) (b 1)) (C.R (b 1) (b 2)) = some v1 ∧
      C.claim (b 2) (C.R (b 1) (b 2)) (C.R (b 2) (b 3)) = some v2 ∧
      C.claim (b 3) (C.R (b 2) (b 3)) 0 = some v3 ∧ 4 * K ≤ v0 + v1 + v2 + v3 := by
  simp only [chainOK, Bool.and_eq_true] at hc
  obtain ⟨⟨⟨h1, h2⟩, h3⟩, he⟩ := hc
  -- position 1
  simp only [chain1, List.all_eq_true, List.mem_range] at h1
  have g1 := h1 (b 0) (hb 0) (b 1) (hb 1)
  split at g1 <;> try simp at g1
  rename_i v0 s1 hv0 hs1
  -- position 2
  have m1 := mem_of_lookup hs1
  simp only [chainStep, List.all_eq_true, List.mem_range] at h2
  have g2 := h2 (b 1) (hb 1) _ m1 (b 2) (hb 2)
  split at g2 <;> try simp at g2
  rename_i v1 s2 hv1 hs2
  have m2 := mem_of_lookup hs2
  simp only [chainStep, List.all_eq_true, List.mem_range] at h3
  have g3 := h3 (b 2) (hb 2) _ m2 (b 3) (hb 3)
  split at g3 <;> try simp at g3
  rename_i v2 s3 hv2 hs3
  have m3 := mem_of_lookup hs3
  simp only [chainEnd, List.all_eq_true, List.mem_range] at he
  have g4 := he (b 3) (hb 3) _ m3
  have hfl : ((C.flag (b 0) || C.flag (b 1) || C.flag (b 2)) || C.flag (b 3)) = true := by
    obtain ⟨k, hk⟩ := hf
    fin_cases k <;> simp_all
  simp only [hfl, Bool.not_true, Bool.false_or] at g4
  split at g4 <;> try simp at g4
  rename_i v3 hv3
  refine ⟨v0, v1, v2, v3, hv0, hv1, hv2, hv3, ?_⟩
  try simp only [decide_eq_true_eq] at g1 g2 g3 g4
  omega

lemma a_sq_le {aN k : ℕ} (hK : (0 : ℝ) < k) (h : (2 * aN + k) ^ 2 ≤ 8 * k ^ 2) :
    ((aN : ℝ) / k + 1 / 2) ^ 2 ≤ 2 := by
  rw [div_add' _ _ _ hK.ne', div_pow, div_le_iff₀ (by positivity)]
  have := (Nat.cast_le (α := ℝ)).mpr h
  push_cast at this
  nlinarith

/-- **The certificate theorem.** -/
theorem wallStrip (hp : C.paramsOK = true) (hr : C.ruleOK = true) (hcov : C.coverOK = true)
    (hf : C.flagsOK = true) (hl : ∀ b < C.m, C.boxLeavesOK b = true) (hc : C.chainOK = true) :
    WallStripBound ((C.hN : ℝ) / K) := by
  have hK := K_pos
  simp only [paramsOK, Bool.and_eq_true, decide_eq_true_eq] at hp
  obtain ⟨⟨hp1, hp2⟩, hp3⟩ := hp
  simp only [ruleOK, List.all_eq_true, List.mem_range, Bool.and_eq_true, decide_eq_true_eq] at hr
  refine wallStrip_of_cert (a := (C.aN : ℝ) / K) C.m (zOf C.menu) C.R C.flag C.inB
    (by positivity) ?_ ?_ ?_ (fun b b' hb hb' => (hr b hb b' hb').1.1) ?_ ?_ ?_ ?_
  · exact a_sq_le hK hp1
  · rw [div_sub' hK.ne', div_lt_div_iff_of_pos_right hK]
    have := (Nat.cast_lt (α := ℝ)).mpr hp2
    push_cast at this
    linarith
  · rw [div_lt_one hK]; exact_mod_cast hp3
  · intro b b' hb hb'
    obtain ⟨⟨-, h2⟩, h3⟩ := hr b hb b' hb'
    simp only [zOf]
    constructor
    · rw [div_sub' hK.ne', div_lt_div_iff_of_pos_right hK]
      have := (Nat.cast_lt (α := ℝ)).mpr h2
      push_cast at this
      linarith
    · rw [div_lt_one hK]; exact_mod_cast h3
  · intro t c ⟨h0, h1, hpc, hch⟩
    exact C.covAll_sound 0 C.cover hcov t c (by simp [h0]) h1 hpc hch
  · intro b t c hb _ hin hfl
    exact C.flag_sound hf hb hfl hin
  · intro b hb hfl
    obtain ⟨v0, v1, v2, v3, e0, e1, e2, e3, hsum⟩ := C.chain_sound hc b hb hfl
    refine ⟨![(v0 : ℝ) / K, (v1 : ℝ) / K, (v2 : ℝ) / K, (v3 : ℝ) / K], ?_, ?_⟩
    · simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.cons_val_two,
        Matrix.cons_val_three]
      have := (Nat.cast_le (α := ℝ)).mpr hsum
      push_cast at this
      have e : (v0 : ℝ) / K + (v1 : ℝ) / K + (v2 : ℝ) / K + (v3 : ℝ) / K =
          ((v0 : ℝ) + v1 + v2 + v3) / K := by ring
      exact le_of_le_of_eq (by rw [le_div_iff₀ hK]; linarith) e.symm
    · intro t c ⟨h0, h1, hpc, _⟩
      refine ⟨fun hin => ?_, fun hin => ?_, fun hin => ?_, fun hin => ?_⟩
      · exact C.claim_sound hl (hb 0) e0 h0 h1 hpc hin
      · exact C.claim_sound hl (hb 1) e1 h0 h1 hpc hin
      · exact C.claim_sound hl (hb 2) e2 h0 h1 hpc hin
      · exact C.claim_sound hl (hb 3) e3 h0 h1 hpc hin

end Cert

end StripCheck
end SquarePacking
