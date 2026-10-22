import Sqpack.RectChain
import Sqpack.StripKernel

/-!
# The slope walk: n unit squares in a strip of height `H < 2`

`strip_walk`: if `n ≥ 1` unit squares pack `[0,L] × [0,H]`, `H < 2`, then (after sorting the squares by centre `x`)
there are frames `(P i, Q i)` and line slopes `m 0 = 0, m 1, …, m n = 0` with
`∑_{i < n} kap (P i) (Q i) H (m i) (m (i+1)) ≤ L`.

* Lines: the walls `x = 0`, `x = L` (slope `0`) and, between consecutive squares, the separating line of
  `separate_lr` normalised to `x + m y = s` (its normal points right).
* `square_step`: a square between lines `(m, s)` (left) and `(m', s')` (right) has
  `kap + (m' - m) H/2 ≤ s' - s` (support function in direction `(1, m)`, then the centre height
  `y ∈ [k, H - k]`).
* Summing, `s` telescopes to `L` and `(m' - m) H/2` to `0`.

`search/rectfront/RECTFRONT.md` (2026-10-22): this is the exact relaxation (numerically) behind the strip frontier.
-/

open Finset

namespace SquarePacking

/-- **One square between two lines.** -/
theorem square_step {c : ℝ × ℝ} {θ L H m s m' s' : ℝ} (hin : sq c θ 1 ⊆ rbox L H)
    (hl : ∀ p ∈ sq c θ 1, s ≤ p.1 + m * p.2) (hr : ∀ p ∈ sq c θ 1, p.1 + m' * p.2 ≤ s') :
    ∃ P Q : ℝ, 0 ≤ P ∧ 0 ≤ Q ∧ P ^ 2 + Q ^ 2 = 1 ∧ P + Q ≤ H ∧
      kap P Q H m m' + (m' - m) * (H / 2) ≤ s' - s := by
  obtain ⟨P, Q, hP, hQ, hPQ, -, -, h3, h4, hf⟩ := framed hin
  refine ⟨P, Q, hP, hQ, hPQ, by linarith, ?_⟩
  have a := support_le (u1 := 1) (u2 := m') hf (fun p hp => by have := hr p hp; linarith)
  have b := support_ge (u1 := 1) (u2 := m) (s := s) hf (fun p hp => by have := hl p hp; linarith)
  -- (m' - m) y ≥ (m' - m) H/2 - |m' - m| (H/2 - k)
  have hy : (m' - m) * (H / 2) - |m - m'| * (H / 2 - (P + Q) / 2) ≤ (m' - m) * c.2 := by
    rw [abs_sub_comm]
    rcases le_total 0 (m' - m) with h | h
    · rw [abs_of_nonneg h]; nlinarith
    · rw [abs_of_nonpos h]; nlinarith
  unfold kap
  have e1 : |-1 * Q + m' * P| = |m' * P - Q| := by ring_nf
  have e2 : |-1 * Q + m * P| = |m * P - Q| := by ring_nf
  have e3 : |1 * P + m' * Q| = |P + m' * Q| := by ring_nf
  have e4 : |1 * P + m * Q| = |P + m * Q| := by ring_nf
  rw [e1, e3] at a; rw [e2, e4] at b
  linarith

/-- The squares sorted by centre `x`, as `ℕ`-indexed families (junk beyond `n`). -/
theorem RectPacks.sorted {n : ℕ} {L H : ℝ} (hp : RectPacks n L H) :
    ∃ (c : ℕ → ℝ × ℝ) (θ : ℕ → ℝ), (∀ i < n, sq (c i) (θ i) 1 ⊆ rbox L H) ∧
      (∀ i j, i < n → j < n → i ≠ j → Disjoint (sqInt (c i) (θ i) 1) (sqInt (c j) (θ j) 1)) ∧
      (∀ i j, i ≤ j → j < n → (c i).1 ≤ (c j).1) := by
  obtain ⟨ctr, ang, hin, hd⟩ := hp
  set σ := Tuple.sort (fun i => (ctr i).1)
  have hmono := Tuple.monotone_sort (fun i => (ctr i).1)
  refine ⟨fun i => if h : i < n then ctr (σ ⟨i, h⟩) else 0, fun i => if h : i < n then ang (σ ⟨i, h⟩) else 0,
    fun i hi => by simp only [hi, dite_true]; exact hin _, fun i j hi hj hij => ?_, fun i j hij hj => ?_⟩
  · simp only [hi, hj, dite_true]
    exact hd _ _ (fun h => hij (by simpa using σ.injective h))
  · have hi : i < n := lt_of_le_of_lt hij hj
    simp only [hi, hj, dite_true]
    exact hmono (show (⟨i, hi⟩ : Fin n) ≤ ⟨j, hj⟩ from hij)

/-- **The slope walk.**  `n ≥ 1` unit squares in `[0,L] × [0,H]`, `H < 2`: frames and slopes with
`m 0 = m n = 0` and `∑ kap ≤ L`. -/
theorem strip_walk {n : ℕ} {L H : ℝ} (hn : 1 ≤ n) (hp : RectPacks n L H) (hH : H < 2) :
    ∃ (P Q m : ℕ → ℝ), (∀ i < n, 0 ≤ P i ∧ 0 ≤ Q i ∧ P i ^ 2 + Q i ^ 2 = 1 ∧ P i + Q i ≤ H) ∧
      m 0 = 0 ∧ m n = 0 ∧ ∑ i ∈ range n, kap (P i) (Q i) H (m i) (m (i + 1)) ≤ L := by
  obtain ⟨c, θ, hin, hd, hx⟩ := hp.sorted
  -- lines 0..n: walls at 0 and n, separating lines between
  have hline : ∀ i : ℕ, ∃ ms : ℝ × ℝ, (i = 0 → ms = (0, 0)) ∧ (i = n → ms = (0, L)) ∧
      (0 < i → i ≤ n → ∀ p ∈ sq (c (i - 1)) (θ (i - 1)) 1, p.1 + ms.1 * p.2 ≤ ms.2) ∧
      (i < n → ∀ p ∈ sq (c i) (θ i) 1, ms.2 ≤ p.1 + ms.1 * p.2) := by
    intro i
    rcases Nat.eq_zero_or_pos i with rfl | hi0
    · refine ⟨(0, 0), fun _ => rfl, fun h => absurd h (by omega), fun h => absurd h (by omega), fun h p hp => ?_⟩
      have := hin 0 h hp; simp only [rbox, Set.mem_ofPred_eq] at this; simp; linarith
    rcases lt_or_ge i n with hin' | hin'
    · obtain ⟨u1, u2, s, hu, ha, hb⟩ := separate_lr (hin (i - 1) (by omega)) (hin i hin')
        (hd (i - 1) i (by omega) hin' (by omega)) (hx (i - 1) i (by omega) hin') hH
      refine ⟨(u2 / u1, s / u1), fun h => absurd h (by omega), fun h => absurd h (by omega), fun _ _ p hp => ?_, fun _ p hp => ?_⟩
      · have := ha p hp
        simp only; rw [le_div_iff₀ hu]; field_simp; linarith
      · have := hb p hp
        simp only; rw [div_le_iff₀ hu]; field_simp; linarith
    · rcases lt_or_eq_of_le hin' with hgt | hi
      · exact ⟨(0, 0), fun h => absurd h (by omega), fun h => absurd h (by omega),
          fun _ h => absurd h (by omega), fun h => absurd h (by omega)⟩
      subst hi
      refine ⟨(0, L), fun h => absurd h (by omega), fun _ => rfl, fun _ _ p hp => ?_, fun h => absurd h (by omega)⟩
      have := hin (n - 1) (by omega) hp; simp only [rbox, Set.mem_ofPred_eq] at this; simp; linarith
  choose ms hms0 hmsn hmsl hmsr using hline
  -- each square
  have hsq : ∀ i < n, ∃ PQ : ℝ × ℝ, 0 ≤ PQ.1 ∧ 0 ≤ PQ.2 ∧ PQ.1 ^ 2 + PQ.2 ^ 2 = 1 ∧ PQ.1 + PQ.2 ≤ H ∧
      kap PQ.1 PQ.2 H (ms i).1 (ms (i + 1)).1 + ((ms (i + 1)).1 - (ms i).1) * (H / 2) ≤
        (ms (i + 1)).2 - (ms i).2 := by
    intro i hi
    obtain ⟨P, Q, h1, h2, h3, h4, h5⟩ := square_step (hin i hi) (hmsr i hi)
      (by simpa using hmsl (i + 1) (by omega) (by omega))
    exact ⟨(P, Q), h1, h2, h3, h4, h5⟩
  choose! PQ hPQ using hsq
  refine ⟨fun i => (PQ i).1, fun i => (PQ i).2, fun i => (ms i).1,
    fun i hi => ⟨(hPQ i hi).1, (hPQ i hi).2.1, (hPQ i hi).2.2.1, (hPQ i hi).2.2.2.1⟩,
    by simp [hms0 0 rfl], by simp [hmsn n rfl], ?_⟩
  have hsum := Finset.sum_le_sum (fun i (hi : i ∈ range n) => (hPQ i (mem_range.1 hi)).2.2.2.2)
  rw [Finset.sum_range_sub (fun i => (ms i).2), Finset.sum_add_distrib, ← Finset.sum_mul,
    Finset.sum_range_sub (fun i => (ms i).1), hms0 0 rfl, hmsn n rfl] at hsum
  simpa using hsum

end SquarePacking
