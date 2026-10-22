import Sqpack.Rect5Pack

/-!
# Strips: the constructions for every `n` (row) and every odd `n` (45° zigzag)

* `strip_row`: `n` upright squares in a row, `n ≤ L`, `1 ≤ H`.
* `strip_zigzag` (odd `n ≥ 3`, `√2 ≤ H ≤ 2`): uprights in both bottom corners, `n - 2` squares at `45°` alternating
  ceiling / floor (first and last on the ceiling, against the uprights' inner corners), consecutive ones edge to edge,
  `λ = 2√2 - H` apart.  Length `F_n = (n - 2)(2√2 - H) + 4 - H`.  For `n = 3, 5` these are Friedman's pieces
  (Math Magic, June 2010); for all odd `n` the same family (the boundary piece `b = 4 + (2n − 4)√2 − (n − 1)a`
  listed by wand125 `problems/rectangles` for `n ≤ 15`).

Every pair of squares is split by one line, found by `linarith` from the four linear bounds of each interior
(`up_int`, `dia_int`): an upright and its neighbouring diamond by `x + y = 2` (mirror: `x - y = L - 2`), an upright
and anything farther by `x = 1` (`x = L - 1`), consecutive diamonds by their common edge line, diamonds two or more
apart by a vertical line (`2λ ≥ √2`).
-/

namespace SquarePacking

/-- `n` upright squares in a row. -/
theorem strip_row {n : ℕ} {L H : ℝ} (hL : (n : ℝ) ≤ L) (hH : 1 ≤ H) : RectPacks n L H := by
  refine ⟨fun i => ((i : ℝ) + 1 / 2, 1 / 2), fun _ => 0, fun i => ?_, fun i j hij => ?_⟩
  · have : ((i : ℕ) : ℝ) + 1 ≤ n := by exact_mod_cast i.isLt
    have h0 : (0 : ℝ) ≤ ((i : ℕ) : ℝ) := Nat.cast_nonneg _
    exact sq_zero_subset (by simp only; linarith) (by simp only; linarith) (by norm_num) (by simp only; linarith)
  · rw [Set.disjoint_left]; intro p h1 h2
    obtain ⟨a1, a2, -, -⟩ := up_int h1 rfl; obtain ⟨b1, b2, -, -⟩ := up_int h2 rfl
    simp only at a1 a2 b1 b2
    have hne : (i : ℕ) ≠ j := fun h => hij (Fin.ext h)
    rcases Nat.lt_or_gt_of_ne hne with h | h
    · have : ((i : ℕ) : ℝ) + 1 ≤ j := by exact_mod_cast h
      linarith
    · have : ((j : ℕ) : ℝ) + 1 ≤ i := by exact_mod_cast h
      linarith

section zigzag
variable (n : ℕ) (L H : ℝ)

/-- Centre of square `i` of the zigzag. -/
noncomputable def zzCtr (i : ℕ) : ℝ × ℝ :=
  if i = 0 then (1 / 2, 1 / 2)
  else if i = n - 1 then (L - 1 / 2, 1 / 2)
  else (2 + Real.sqrt 2 - H + ((i : ℝ) - 1) * (2 * Real.sqrt 2 - H),
    if i % 2 = 1 then H - Real.sqrt 2 / 2 else Real.sqrt 2 / 2)

/-- Angle of square `i`: the two ends upright, the rest at `45°`. -/
noncomputable def zzAng (i : ℕ) : ℝ := if i = 0 ∨ i = n - 1 then 0 else Real.pi / 4

lemma zz_first : zzCtr n L H 0 = (1 / 2, 1 / 2) ∧ zzAng n 0 = 0 := by simp [zzCtr, zzAng]

lemma zz_last (hn : n - 1 ≠ 0) : zzCtr n L H (n - 1) = (L - 1 / 2, 1 / 2) ∧ zzAng n (n - 1) = 0 := by
  simp [zzCtr, zzAng, hn]

lemma zz_mid {i : ℕ} (h0 : i ≠ 0) (hN : i ≠ n - 1) :
    zzCtr n L H i = (2 + Real.sqrt 2 - H + ((i : ℝ) - 1) * (2 * Real.sqrt 2 - H),
      if i % 2 = 1 then H - Real.sqrt 2 / 2 else Real.sqrt 2 / 2) ∧ zzAng n i = Real.pi / 4 := by
  simp [zzCtr, zzAng, h0, hN]

end zigzag

/-- **The 45° zigzag** for odd `n ≥ 3`. -/
theorem strip_zigzag {n : ℕ} {L H : ℝ} (hn : n % 2 = 1) (hn3 : 3 ≤ n) (hH1 : Real.sqrt 2 ≤ H) (hH2 : H ≤ 2)
    (hL : ((n : ℝ) - 2) * (2 * Real.sqrt 2 - H) + 4 - H ≤ L) : RectPacks n L H := by
  have r2 : Real.sqrt 2 * Real.sqrt 2 = 2 := Real.mul_self_sqrt (by norm_num)
  have r2g : 1.4142 < Real.sqrt 2 := by rw [Real.lt_sqrt (by norm_num)]; norm_num
  have r2l : Real.sqrt 2 < 1.4143 := by rw [Real.sqrt_lt' (by norm_num)]; norm_num
  have hn' : (3 : ℝ) ≤ n := by exact_mod_cast hn3
  set lam := 2 * Real.sqrt 2 - H with hlam
  set x1 := 2 + Real.sqrt 2 - H
  have hlam0 : 0 < lam := by simp only [lam]; linarith
  -- the diamond centres' x as a function of the index
  have xmono : ∀ i j : ℕ, i ≤ j → x1 + ((i : ℝ) - 1) * lam ≤ x1 + ((j : ℝ) - 1) * lam := by
    intro i j h; have : (i : ℝ) ≤ j := by exact_mod_cast h
    nlinarith
  -- the last diamond is n - 2
  have xlast : ∀ i : ℕ, i ≤ n - 2 → x1 + ((i : ℝ) - 1) * lam + x1 ≤ L := by
    intro i hi
    have : (i : ℝ) ≤ n - 2 := by
      have : ((n - 2 : ℕ) : ℝ) = n - 2 := by rw [Nat.cast_sub (by omega)]; norm_num
      rw [← this]; exact_mod_cast hi
    simp only [x1, lam] at *; nlinarith
  have hL2 : 2 ≤ L := by nlinarith
  refine ⟨fun i => zzCtr n L H i, fun i => zzAng n i, fun i => ?_, fun i j hij => ?_⟩
  · -- containment
    have hi := i.isLt
    show sq (zzCtr n L H i) (zzAng n i) 1 ⊆ rbox L H
    by_cases h0 : (i : ℕ) = 0
    · rw [h0, (zz_first n L H).1, (zz_first n L H).2]
      exact sq_zero_subset (by norm_num) (by simp only; linarith) (by norm_num) (by simp only; linarith)
    by_cases hN : (i : ℕ) = n - 1
    · rw [hN, (zz_last n L H (by omega)).1, (zz_last n L H (by omega)).2]
      exact sq_zero_subset (by simp only; linarith) (by simp only; linarith) (by norm_num) (by simp only; linarith)
    rw [(zz_mid n L H h0 hN).1, (zz_mid n L H h0 hN).2]
    have hi1 : (1 : ℝ) ≤ (i : ℕ) := by exact_mod_cast (Nat.one_le_iff_ne_zero.2 h0)
    have := xlast i (by omega)
    apply dia_subset
    · simp only; nlinarith
    · simp only; nlinarith
    · simp only; split_ifs <;> linarith
    · simp only; split_ifs <;> linarith
  · -- disjointness
    show Disjoint (sqInt (zzCtr n L H i) (zzAng n i) 1) (sqInt (zzCtr n L H j) (zzAng n j) 1)
    -- interior bounds of each kind of square, by index
    have first : ∀ p ∈ sqInt (zzCtr n L H 0) (zzAng n 0) 1, p.1 < 1 ∧ p.2 < 1 := fun p hp => by
      rw [(zz_first n L H).1, (zz_first n L H).2] at hp
      obtain ⟨-, a2, -, a4⟩ := up_int hp rfl; simp only at a2 a4; constructor <;> linarith
    have last : ∀ p ∈ sqInt (zzCtr n L H (n - 1)) (zzAng n (n - 1)) 1, L - 1 < p.1 ∧ p.2 < 1 := fun p hp => by
      rw [(zz_last n L H (by omega)).1, (zz_last n L H (by omega)).2] at hp
      obtain ⟨a1, -, -, a4⟩ := up_int hp rfl; simp only at a1 a4; constructor <;> linarith
    have mid : ∀ k : ℕ, k ≠ 0 → k ≠ n - 1 → ∀ p ∈ sqInt (zzCtr n L H k) (zzAng n k) 1,
        let X := x1 + ((k : ℝ) - 1) * lam
        let Y := if k % 2 = 1 then H - Real.sqrt 2 / 2 else Real.sqrt 2 / 2
        X + Y - Real.sqrt 2 / 2 < p.1 + p.2 ∧ p.1 + p.2 < X + Y + Real.sqrt 2 / 2 ∧
          Y - X - Real.sqrt 2 / 2 < p.2 - p.1 ∧ p.2 - p.1 < Y - X + Real.sqrt 2 / 2 := by
      intro k h0 hN p hp
      rw [(zz_mid n L H h0 hN).1, (zz_mid n L H h0 hN).2] at hp
      exact dia_int hp rfl
    -- the generic pair argument for i < j
    have key : ∀ i j : ℕ, i < j → j < n → Disjoint (sqInt (zzCtr n L H i) (zzAng n i) 1)
        (sqInt (zzCtr n L H j) (zzAng n j) 1) := by
      intro i j hij hjn
      rw [Set.disjoint_left]; intro p h1 h2
      have hj1 : (1 : ℝ) ≤ (j : ℕ) := by exact_mod_cast (show 1 ≤ j by omega)
      have hij' : (i : ℝ) + 1 ≤ j := by exact_mod_cast hij
      by_cases hi0 : i = 0
      · subst hi0
        obtain ⟨a1, a2⟩ := first p h1
        by_cases hjN : j = n - 1
        · subst hjN; obtain ⟨b1, -⟩ := last p h2; linarith
        obtain ⟨b1, b2, b3, b4⟩ := mid j (by omega) hjN p h2
        by_cases hj1' : j = 1
        · -- the first diamond touches the corner (1,1): split by x + y = 2
          subst hj1'; simp only [show (1 : ℕ) % 2 = 1 from rfl, if_true, Nat.cast_one] at b1; linarith
        · -- farther diamonds: split by x = 1
          have : (2 : ℝ) ≤ (j : ℕ) := by exact_mod_cast (show 2 ≤ j by omega)
          split_ifs at b1 b2 b3 b4 <;> nlinarith
      have hi1 : (1 : ℝ) ≤ (i : ℕ) := by exact_mod_cast (show 1 ≤ i by omega)
      obtain ⟨a1, a2, a3, a4⟩ := mid i hi0 (by omega) p h1
      by_cases hjN : j = n - 1
      · subst hjN
        obtain ⟨b1, b2⟩ := last p h2
        have hxl := xlast i (by omega)
        by_cases hiN : i = n - 2
        · -- the last diamond touches the corner (L - 1, 1): split by x - y = L - 2
          have hodd : i % 2 = 1 := by omega
          simp only [hodd, if_true] at a1 a2 a3 a4; linarith
        · have : (i : ℝ) + 1 ≤ ((n - 2 : ℕ) : ℝ) := by exact_mod_cast (show i + 1 ≤ n - 2 by omega)
          have hx2 := xlast (n - 2) le_rfl
          split_ifs at a1 a2 a3 a4 <;> nlinarith
      obtain ⟨b1, b2, b3, b4⟩ := mid j (by omega) hjN p h2
      by_cases hsucc : j = i + 1
      · -- consecutive diamonds: common edge line
        subst hsucc
        push_cast at b1 b2 b3 b4
        rcases Nat.mod_two_eq_zero_or_one i with hp | hp
        · have : (i + 1) % 2 = 1 := by omega
          simp only [hp, this, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false] at a1 a2 a3 a4 b1 b2 b3 b4
          nlinarith
        · have : (i + 1) % 2 = 0 := by omega
          simp only [hp, this, if_true, show (0 : ℕ) ≠ 1 by norm_num, if_false] at a1 a2 a3 a4 b1 b2 b3 b4
          nlinarith
      · -- two or more apart: a vertical line (2λ ≥ √2)
        have h2 : (i : ℝ) + 2 ≤ j := by exact_mod_cast (show i + 2 ≤ j by omega)
        have hgap : 2 * lam ≤ (x1 + ((j : ℕ) - 1) * lam) - (x1 + ((i : ℕ) - 1) * lam) := by
          nlinarith [mul_le_mul_of_nonneg_right h2 hlam0.le]
        simp only [lam] at hgap
        linarith
    have hne : (i : ℕ) ≠ j := fun h => hij (Fin.ext h)
    rcases Nat.lt_or_gt_of_ne hne with h | h
    · exact key i j h j.isLt
    · exact (key j i h i.isLt).symm

end SquarePacking
