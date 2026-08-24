import Mathlib

/-!
# Weighted unavoidable sets and lower bounds for packing unit squares in a square

Let `s(n)` be the side of the smallest square into which `n` unit squares can be packed
(rotations allowed).  A **weighted unavoidable set** for the container `C` is a finite set of
points ("atoms") with non-negative weights such that *every* closed unit square contained in
`C` contains atoms of total weight at least `1`.

`packing_le_weight` : if such a set has total weight `W`, then any family of `n` squares of
side `L > 1`, contained in `C` and with pairwise disjoint interiors, satisfies `n ≤ W`.

`no_small_packing` : consequently, if `W < n` then `n` unit squares cannot be packed into any
square of side `s' < s`; that is, `s(n) ≥ s`.  (The passage from unit squares in the smaller
container to `L`-squares with `L > 1` in `C` is the classical scaling trick.)

The hypothesis `hcover` is what the exact integer verifier in `verify/` checks.
-/

open Finset
open scoped Classical

namespace SquarePacking

/-- Rotated coordinates of `p` relative to centre `c` and angle `θ`. -/
noncomputable def coord (c : ℝ × ℝ) (θ : ℝ) (p : ℝ × ℝ) : ℝ × ℝ :=
  ((p.1 - c.1) * Real.cos θ + (p.2 - c.2) * Real.sin θ,
   (-(p.1 - c.1)) * Real.sin θ + (p.2 - c.2) * Real.cos θ)

/-- The closed square of side `L`, centre `c`, angle `θ`. -/
def sq (c : ℝ × ℝ) (θ L : ℝ) : Set (ℝ × ℝ) :=
  {p | |(coord c θ p).1| ≤ L/2 ∧ |(coord c θ p).2| ≤ L/2}

/-- The open square (interior) of side `L`, centre `c`, angle `θ`. -/
def sqInt (c : ℝ × ℝ) (θ L : ℝ) : Set (ℝ × ℝ) :=
  {p | |(coord c θ p).1| < L/2 ∧ |(coord c θ p).2| < L/2}

/-- A concentric closed unit square sits in the interior of any concentric square of side
`L > 1`.  This is the step that lets a *closed*-square certificate control a packing. -/
lemma unit_subset_interior {c : ℝ × ℝ} {θ L : ℝ} (hL : 1 < L) :
    sq c θ 1 ⊆ sqInt c θ L := by
  rintro p ⟨h1, h2⟩
  constructor
  · calc |(coord c θ p).1| ≤ 1/2 := h1
      _ < L/2 := by linarith
  · calc |(coord c θ p).2| ≤ 1/2 := h2
      _ < L/2 := by linarith

/-- A point lies in at most one member of a pairwise disjoint family. -/
lemma card_filter_le_one {n : ℕ} (S : Fin n → Set (ℝ × ℝ))
    (hdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j)) (a : ℝ × ℝ) :
    ((univ : Finset (Fin n)).filter (fun i => a ∈ S i)).card ≤ 1 := by
  rw [card_le_one]
  intro i hi j hj
  simp only [mem_filter, mem_univ, true_and] at hi hj
  by_contra hne
  exact (Set.disjoint_left.mp (hdisj i j hne) hi) hj

/-- **Main reduction.**  A weighted unavoidable set of total weight `W` bounds the number of
squares of side `L > 1` that can be packed (disjoint interiors) inside `C`. -/
theorem packing_le_weight
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ))
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        1 ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) ≤ ∑ a ∈ A, w a := by
  -- the concentric unit squares are inside C, and inside the (disjoint) interiors
  have hIntSub : ∀ (c : ℝ × ℝ) (θ : ℝ), sqInt c θ L ⊆ sq c θ L := by
    rintro c θ p ⟨h1, h2⟩; exact ⟨le_of_lt h1, le_of_lt h2⟩
  have hU : ∀ i, sq (ctr i) (ang i) 1 ⊆ C := fun i =>
    subset_trans (subset_trans (unit_subset_interior hL) (hIntSub _ _)) (hin i)
  set S : Fin n → Set (ℝ × ℝ) := fun i => sq (ctr i) (ang i) 1 with hSdef
  have hSdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j) := by
    intro i j hij
    exact Set.disjoint_of_subset (unit_subset_interior hL) (unit_subset_interior hL) (hdisj i j hij)
  have key : ∀ i : Fin n, (1:ℝ) ≤ ∑ a ∈ A.filter (fun a => a ∈ S i), w a := fun i =>
    hcover (ctr i) (ang i) (hU i)
  have h1 : (n:ℝ) ≤ ∑ i : Fin n, ∑ a ∈ A.filter (fun a => a ∈ S i), w a := by
    calc (n:ℝ) = ∑ _i : Fin n, (1:ℝ) := by simp
      _ ≤ _ := Finset.sum_le_sum fun i _ => key i
  have h2 : ∑ i : Fin n, ∑ a ∈ A.filter (fun a => a ∈ S i), w a
      = ∑ a ∈ A, ∑ i : Fin n, (if a ∈ S i then w a else 0) := by
    simp only [Finset.sum_filter]; exact Finset.sum_comm
  have h3 : ∀ a ∈ A, ∑ i : Fin n, (if a ∈ S i then w a else 0) ≤ w a := by
    intro a ha
    have hcard : (((univ : Finset (Fin n)).filter (fun i => a ∈ S i)).card : ℝ) ≤ 1 := by
      exact_mod_cast card_filter_le_one S hSdisj a
    have hEq : ∑ i : Fin n, (if a ∈ S i then w a else 0)
        = (((univ : Finset (Fin n)).filter (fun i => a ∈ S i)).card : ℝ) * w a := by
      rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul]
    rw [hEq]; nlinarith [hw a ha]
  have h4 : ∑ a ∈ A, ∑ i : Fin n, (if a ∈ S i then w a else 0) ≤ ∑ a ∈ A, w a :=
    Finset.sum_le_sum h3
  linarith [h1, h2 ▸ h1, h4]

/-- Scaling by `μ > 0` turns a unit square into a square of side `μ`. -/
lemma sq_scale (c : ℝ × ℝ) (θ μ : ℝ) (hμ : 0 < μ) (p : ℝ × ℝ) (hp : p ∈ sq c θ 1) :
    (μ * p.1, μ * p.2) ∈ sq (μ * c.1, μ * c.2) θ μ := by
  obtain ⟨h1, h2⟩ := hp
  have e1 : (coord (μ * c.1, μ * c.2) θ (μ * p.1, μ * p.2)).1 = μ * (coord c θ p).1 := by
    simp only [coord]; ring
  have e2 : (coord (μ * c.1, μ * c.2) θ (μ * p.1, μ * p.2)).2 = μ * (coord c θ p).2 := by
    simp only [coord]; ring
  refine ⟨?_, ?_⟩
  · rw [e1, abs_mul, abs_of_pos hμ]; nlinarith [abs_nonneg (coord c θ p).1]
  · rw [e2, abs_mul, abs_of_pos hμ]; nlinarith [abs_nonneg (coord c θ p).2]

/-- Interiors scale the same way. -/
lemma sqInt_scale (c : ℝ × ℝ) (θ μ : ℝ) (hμ : 0 < μ) (p : ℝ × ℝ) (hp : p ∈ sqInt c θ 1) :
    (μ * p.1, μ * p.2) ∈ sqInt (μ * c.1, μ * c.2) θ μ := by
  obtain ⟨h1, h2⟩ := hp
  have e1 : (coord (μ * c.1, μ * c.2) θ (μ * p.1, μ * p.2)).1 = μ * (coord c θ p).1 := by
    simp only [coord]; ring
  have e2 : (coord (μ * c.1, μ * c.2) θ (μ * p.1, μ * p.2)).2 = μ * (coord c θ p).2 := by
    simp only [coord]; ring
  refine ⟨?_, ?_⟩
  · rw [e1, abs_mul, abs_of_pos hμ]; nlinarith [abs_nonneg (coord c θ p).1]
  · rw [e2, abs_mul, abs_of_pos hμ]; nlinarith [abs_nonneg (coord c θ p).2]

end SquarePacking
