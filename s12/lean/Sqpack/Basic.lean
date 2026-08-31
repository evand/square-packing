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

/-- **Threshold form of the reduction.**  If every closed unit square inside `C` with centre `c`
and angle `θ` captures weight `≥ t c θ`, then a packing of `n` squares of side `L > 1` satisfies
`∑ i, t (ctr i) (ang i) ≤ W`.  `packing_le_weight` is the case `t = 1`. -/
theorem packing_le_weight_thresh
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (t : ℝ × ℝ → ℝ → ℝ)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        t c θ ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    ∑ i : Fin n, t (ctr i) (ang i) ≤ ∑ a ∈ A, w a := by
  have hIntSub : ∀ (c : ℝ × ℝ) (θ : ℝ), sqInt c θ L ⊆ sq c θ L := by
    rintro c θ p ⟨h1, h2⟩; exact ⟨le_of_lt h1, le_of_lt h2⟩
  have hU : ∀ i, sq (ctr i) (ang i) 1 ⊆ C := fun i =>
    subset_trans (subset_trans (unit_subset_interior hL) (hIntSub _ _)) (hin i)
  set S : Fin n → Set (ℝ × ℝ) := fun i => sq (ctr i) (ang i) 1 with hSdef
  have hSdisj : ∀ i j, i ≠ j → Disjoint (S i) (S j) := by
    intro i j hij
    exact Set.disjoint_of_subset (unit_subset_interior hL) (unit_subset_interior hL) (hdisj i j hij)
  have key : ∀ i : Fin n, t (ctr i) (ang i) ≤ ∑ a ∈ A.filter (fun a => a ∈ S i), w a := fun i =>
    hcover (ctr i) (ang i) (hU i)
  have h1 : ∑ i : Fin n, t (ctr i) (ang i) ≤ ∑ i : Fin n, ∑ a ∈ A.filter (fun a => a ∈ S i), w a :=
    Finset.sum_le_sum fun i _ => key i
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

/-- **Branch reduction.**  Let `R` be a region of poses (centre, angle) and `λ` a real number.
If every closed unit square inside `C` whose pose lies in `R` captures weight `≥ 1 + λ` and every
other one captures `≥ 1`, then a packing of `n` squares of side `L > 1` with exactly `k` squares
whose pose lies in `R` satisfies `n + λ k ≤ W`.  This is what a branch certificate
(`certificates/FORMAT.md`) asserts; with `λ = 0` it is `packing_le_weight`. -/
theorem packing_le_weight_region
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (R : ℝ × ℝ → ℝ → Prop) (lam : ℝ)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        (if R c θ then 1 + lam else 1) ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) + lam * (((univ : Finset (Fin n)).filter
        (fun i => R (ctr i) (ang i))).card : ℝ) ≤ ∑ a ∈ A, w a := by
  have h := packing_le_weight_thresh A w hw C (fun c θ => if R c θ then 1 + lam else 1) hcover
    n L hL ctr ang hin hdisj
  have e : ∀ i : Fin n, (if R (ctr i) (ang i) then (1:ℝ) + lam else 1)
      = 1 + lam * (if R (ctr i) (ang i) then 1 else 0) := by
    intro i; split_ifs <;> simp
  have hsum : ∑ i : Fin n, (if R (ctr i) (ang i) then (1:ℝ) + lam else 1)
      = (n : ℝ) + lam * (((univ : Finset (Fin n)).filter
          (fun i => R (ctr i) (ang i))).card : ℝ) := by
    simp_rw [e]
    rw [Finset.sum_add_distrib, ← Finset.mul_sum, Finset.sum_boole]
    simp
  simpa [hsum] using h

/-- **Several regions.**  Regions `R j` with multipliers `lam j`, `j : Fin m`: a square whose pose
lies in `R j` must capture `1 + lam j` more than a square outside every region (the
thresholds add if regions overlap); then `n + ∑ j, lam j · #{i | pose i ∈ R j} ≤ W`.  This is
the per-corner branch certificate (`lambda L1 L2 L3 L4 / k K1 K2 K3 K4`). -/
theorem packing_le_weight_regions
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (m : ℕ) (R : Fin m → ℝ × ℝ → ℝ → Prop) (lam : Fin m → ℝ)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        1 + ∑ j : Fin m, (if R j c θ then lam j else 0)
          ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) + ∑ j : Fin m, lam j * (((univ : Finset (Fin n)).filter
        (fun i => R j (ctr i) (ang i))).card : ℝ) ≤ ∑ a ∈ A, w a := by
  have h := packing_le_weight_thresh A w hw C
    (fun c θ => 1 + ∑ j : Fin m, (if R j c θ then lam j else 0)) hcover n L hL ctr ang hin hdisj
  have hsum : ∑ i : Fin n, (1 + ∑ j : Fin m, (if R j (ctr i) (ang i) then lam j else 0))
      = (n : ℝ) + ∑ j : Fin m, lam j * (((univ : Finset (Fin n)).filter
          (fun i => R j (ctr i) (ang i))).card : ℝ) := by
    rw [Finset.sum_add_distrib, Finset.sum_comm]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
    congr 1
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, mul_comm]
  simpa [hsum] using h

/-- **Cliques.**  A set of poses `K` is a *clique* if any two closed unit squares with poses in it
share a point.  In a packing (pairwise disjoint interiors of `L`-squares, `L > 1`) at most one
square has its pose in a clique: the concentric closed unit squares of two such squares would
share a point lying in both interiors. -/
lemma card_filter_clique_le_one {n : ℕ} (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L))
    (K : ℝ × ℝ → ℝ → Prop)
    (hK : ∀ (c c' : ℝ × ℝ) (θ θ' : ℝ), K c θ → K c' θ' → (sq c θ 1 ∩ sq c' θ' 1).Nonempty) :
    ((univ : Finset (Fin n)).filter (fun i => K (ctr i) (ang i))).card ≤ 1 := by
  rw [card_le_one]
  intro i hi j hj
  simp only [mem_filter, mem_univ, true_and] at hi hj
  by_contra hne
  obtain ⟨p, hpi, hpj⟩ := hK _ _ _ _ hi hj
  exact Set.disjoint_left.mp (hdisj i j hne) (unit_subset_interior hL hpi) (unit_subset_interior hL hpj)

/-- **Box cliques are cliques.**  If the poses of a family are split into pieces `B i` (the boxes),
each piece has a *core* `core i` contained in every closed unit square with a pose in the piece,
and the cores pairwise meet, then any two squares of the family share a point.  This is exactly
what the verifier checks for a `cliques` block (`certificates/FORMAT.md`): the core of a box is the
concentric shrunk square common to all its poses, and the pairwise test is an exact separating-axis
test.  The inclusion `core i ⊆ sq c θ 1` is the verifier's shrink lemma (a unit square at any
angle of a bin contains the concentric `σ_k`-square at the bin's first angle), stated here as the
hypothesis `hcore`. -/
lemma clique_of_cores {ι : Type*} (B : ι → ℝ × ℝ → ℝ → Prop) (core : ι → Set (ℝ × ℝ))
    (hcore : ∀ i (c : ℝ × ℝ) (θ : ℝ), B i c θ → core i ⊆ sq c θ 1)
    (hmeet : ∀ i j, (core i ∩ core j).Nonempty) :
    ∀ (c c' : ℝ × ℝ) (θ θ' : ℝ), (∃ i, B i c θ) → (∃ j, B j c' θ') →
      (sq c θ 1 ∩ sq c' θ' 1).Nonempty := by
  rintro c c' θ θ' ⟨i, hi⟩ ⟨j, hj⟩
  exact (hmeet i j).mono (Set.inter_subset_inter (hcore i c θ hi) (hcore j c' θ' hj))

/-- **Clique reduction.**  Points `A` with weights `w ≥ 0` and cliques `K j` with weights `v j ≥ 0`
(`j : Fin m`): if every closed unit square inside `C` is covered by weight `≥ 1` by the points it
contains together with the cliques it belongs to, then a packing of `n` squares of side `L > 1`
satisfies `n ≤ ∑ w + ∑ v`.  Points are the special case of cliques (`K = {poses containing p}`),
and `packing_le_weight` is the case `m = 0`. -/
theorem packing_le_weight_cliques
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (m : ℕ) (K : Fin m → ℝ × ℝ → ℝ → Prop) (v : Fin m → ℝ) (hv : ∀ j, 0 ≤ v j)
    (hK : ∀ j (c c' : ℝ × ℝ) (θ θ' : ℝ), K j c θ → K j c' θ' → (sq c θ 1 ∩ sq c' θ' 1).Nonempty)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        1 ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a
              + ∑ j : Fin m, (if K j c θ then v j else 0))
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) ≤ ∑ a ∈ A, w a + ∑ j : Fin m, v j := by
  have h := packing_le_weight_thresh A w hw C
    (fun c θ => 1 - ∑ j : Fin m, (if K j c θ then v j else 0))
    (fun c θ hc => by linarith [hcover c θ hc]) n L hL ctr ang hin hdisj
  have hsum : ∑ i : Fin n, (1 - ∑ j : Fin m, (if K j (ctr i) (ang i) then v j else 0))
      = (n : ℝ) - ∑ j : Fin m, v j * (((univ : Finset (Fin n)).filter
          (fun i => K j (ctr i) (ang i))).card : ℝ) := by
    rw [Finset.sum_sub_distrib, Finset.sum_comm]
    simp only [Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, mul_one]
    congr 1
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← Finset.sum_filter, Finset.sum_const, nsmul_eq_mul, mul_comm]
  have hcard : ∀ j : Fin m, v j * (((univ : Finset (Fin n)).filter
      (fun i => K j (ctr i) (ang i))).card : ℝ) ≤ v j := by
    intro j
    have h1 : (((univ : Finset (Fin n)).filter (fun i => K j (ctr i) (ang i))).card : ℝ) ≤ 1 := by
      exact_mod_cast card_filter_clique_le_one L hL ctr ang hdisj (K j) (hK j)
    nlinarith [hv j]
  have h2 : ∑ j : Fin m, v j * (((univ : Finset (Fin n)).filter
      (fun i => K j (ctr i) (ang i))).card : ℝ) ≤ ∑ j : Fin m, v j :=
    Finset.sum_le_sum fun j _ => hcard j
  rw [hsum] at h
  linarith

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
