import Sqpack.RyuEnds

/-!
# The staircase end: geometry (after Ryu, Section 6)

S. Ryu, *k² − c upper bounds* (squarepacker/k2-minus-c-upper v1.1; jlevy/squares#471, #486), Section 6, Lemma 6.3
(the idea of replacing the lower squares of tilted columns by horizontal layers is from Bui's Section 3, as Ryu
credits).  An end region `endReg γ σ w y` (`RyuEnds`) is packed by

* **columns**: for `i < I`, `n_i` squares at one tilt `α` (`ca = cos α`, `sa = sin α`), stacked along
  `v = (−sa, ca)` from the base point `(x_i, J_i)`;
* **layers**: for `j < J₁`, `L_j` axis squares `[w − L_j, w] × [j, j + 1]`;
* **wedge rows**: for `r < R₀`, `l_r` axis squares `[0, l_r] × [r, r + 1]`.

`staircase`: if each piece fits under the ceiling and the pieces are separated pairwise by the rotated coordinate
`U = ca X + sa Y` (columns occupy `U`-slabs of width `1`, ordered left to right) or by height, they pack the end
region.  The choice of the parameters, and the count, are in `RyuStairCount`.

`PacksIn.biUnion`: packings of finitely many closed regions with pairwise disjoint interiors combine.
-/

namespace UnitSquarePacking.Ryu

open Set UnitSquarePacking.Conjectures UnitSquarePacking.FangGlue

/-! ## Finite unions of closed regions -/

lemma disjoint_interior_biUnion {ι : Type*} (s : Finset ι) (S : ι → Set (ℝ × ℝ)) {U : Set (ℝ × ℝ)}
    (hU : IsOpen U) (hc : ∀ i ∈ s, IsClosed (S i)) (hd : ∀ i ∈ s, Disjoint U (interior (S i))) :
    Disjoint U (interior (⋃ i ∈ s, S i)) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih =>
    rw [Finset.set_biUnion_insert, Set.union_comm]
    exact disjoint_interior_union hU (hc a (Finset.mem_insert_self a s))
      (ih (fun i hi => hc i (Finset.mem_insert_of_mem hi)) (fun i hi => hd i (Finset.mem_insert_of_mem hi)))
      (hd a (Finset.mem_insert_self a s))

/-- Packings of finitely many closed regions with pairwise disjoint interiors combine. -/
lemma _root_.UnitSquarePacking.Conjectures.PacksIn.biUnion {ι : Type*} (s : Finset ι) (S : ι → Set (ℝ × ℝ))
    (n : ι → ℕ) (hc : ∀ i ∈ s, IsClosed (S i)) (hp : ∀ i ∈ s, PacksIn (n i) (S i))
    (hd : ∀ i ∈ s, ∀ j ∈ s, i ≠ j → Disjoint (interior (S i)) (interior (S j))) :
    PacksIn (∑ i ∈ s, n i) (⋃ i ∈ s, S i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simpa using PacksIn.zero _
  | insert a s ha ih =>
    rw [Finset.sum_insert ha, Finset.set_biUnion_insert]
    have ih' := ih (fun i hi => hc i (Finset.mem_insert_of_mem hi)) (fun i hi => hp i (Finset.mem_insert_of_mem hi))
      (fun i hi j hj hij => hd i (Finset.mem_insert_of_mem hi) j (Finset.mem_insert_of_mem hj) hij)
    refine (hp a (Finset.mem_insert_self a s)).union ih' ?_
    refine disjoint_interior_biUnion s S isOpen_interior (fun i hi => hc i (Finset.mem_insert_of_mem hi))
      (fun i hi => ?_)
    exact hd a (Finset.mem_insert_self a s) i (Finset.mem_insert_of_mem hi) (fun h => ha (h ▸ hi))

/-! ## One tilted column -/

/-- The region of a column: `U ∈ [U(P), U(P) + 1]`, `Y ∈ [J, J + sa + n ca]`, `X ∈ [x₀ − n sa, x₀ + ca]`, for the
base point `P = (x₀, J)`. -/
def colSet (ca sa x₀ J n : ℝ) : Set (ℝ × ℝ) :=
  {p | ca * x₀ + sa * J ≤ ca * p.1 + sa * p.2 ∧ ca * p.1 + sa * p.2 ≤ ca * x₀ + sa * J + 1 ∧
    J ≤ p.2 ∧ p.2 ≤ J + sa + n * ca ∧ x₀ - n * sa ≤ p.1 ∧ p.1 ≤ x₀ + ca}

lemma isClosed_colSet (ca sa x₀ J n : ℝ) : IsClosed (colSet ca sa x₀ J n) := by
  have cU : Continuous fun p : ℝ × ℝ => ca * p.1 + sa * p.2 := by fun_prop
  exact (isClosed_le continuous_const cU).inter ((isClosed_le cU continuous_const).inter
    ((isClosed_le continuous_const continuous_snd).inter ((isClosed_le continuous_snd continuous_const).inter
    ((isClosed_le continuous_const continuous_fst).inter (isClosed_le continuous_fst continuous_const)))))

/-- **A column**: `n` squares at angle `α`, stacked along `(−sin α, cos α)` from the base point `(x₀, J)`. -/
theorem column {α : ℝ} (hca : 0 ≤ Real.cos α) (hsa : 0 ≤ Real.sin α) (x₀ : ℝ) (J n : ℕ) :
    PacksIn n (colSet (Real.cos α) (Real.sin α) x₀ J n) := by
  set ca := Real.cos α with hcad
  set sa := Real.sin α with hsad
  have hu : ca ^ 2 + sa ^ 2 = 1 := by rw [hcad, hsad, add_comm]; exact Real.sin_sq_add_cos_sq α
  have hca1 : ca ≤ 1 := by nlinarith
  have hsa1 : sa ≤ 1 := by nlinarith
  set c : Fin n → ℝ × ℝ := fun e => (x₀ + ca / 2 - (e + 1 / 2) * sa, J + sa / 2 + (e + 1 / 2) * ca) with hc
  refine ⟨c, fun _ => α, fun e p hp => ?_, fun e e' hee => ?_⟩
  · have he : (e : ℝ) + 1 ≤ n := by exact_mod_cast e.isLt
    have he0 : (0 : ℝ) ≤ e := by positivity
    -- each of the six constraints, from the corners
    have K : ∀ a b : ℝ, Half a → Half b →
        (0 ≤ 1 / 2 + a ∧ 1 / 2 + a ≤ 1) ∧ (0 ≤ e + 1 / 2 + b ∧ e + 1 / 2 + b ≤ n) := by
      intro a b ha hb
      rcases ha with rfl | rfl <;> rcases hb with rfl | rfl <;> norm_num <;> constructor <;> linarith
    have cx : ∀ a b : ℝ, (c e).1 + ca * a - sa * b = x₀ + ca * (1 / 2 + a) - sa * (e + 1 / 2 + b) := by
      intro a b; simp only [hc]; ring
    have cy : ∀ a b : ℝ, (c e).2 + sa * a + ca * b = J + sa * (1 / 2 + a) + ca * (e + 1 / 2 + b) := by
      intro a b; simp only [hc]; ring
    have h1 := unitSq_subset_half (c₀ := c e) (θ := α) (α := -ca) (β := -sa) (c := -(ca * x₀ + sa * J)) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    have h2 := unitSq_subset_half (c₀ := c e) (θ := α) (α := ca) (β := sa) (c := ca * x₀ + sa * J + 1) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    have h3 := unitSq_subset_half (c₀ := c e) (θ := α) (α := 0) (β := -1) (c := -J) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    have h4 := unitSq_subset_half (c₀ := c e) (θ := α) (α := 0) (β := 1) (c := J + sa + n * ca) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    have h5 := unitSq_subset_half (c₀ := c e) (θ := α) (α := -1) (β := 0) (c := -(x₀ - n * sa)) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    have h6 := unitSq_subset_half (c₀ := c e) (θ := α) (α := 1) (β := 0) (c := x₀ + ca) (fun a b ha hb => by
      rw [← hcad, ← hsad, cx, cy]; obtain ⟨⟨k1, k2⟩, k3, k4⟩ := K a b ha hb; nlinarith) hp
    simp only [mem_ofPred_eq] at h1 h2 h3 h4 h5 h6
    exact ⟨by linarith, by linarith, by linarith, by linarith, by linarith, by linarith⟩
  · -- same angle, `V` differs by `e − e'`
    apply sep_V
    rw [← hcad, ← hsad]
    have hv : ∀ e : Fin n, V ca sa (c e) = -x₀ * sa + J * ca + (e + 1 / 2) := by
      intro e; simp only [V, hc]; linear_combination ((e : ℝ) + 1 / 2) * hu
    rw [hv, hv]
    have hne : (e : ℕ) ≠ e' := fun h => hee (Fin.ext h)
    rcases Nat.lt_or_gt_of_ne hne with h | h
    · have : (e : ℝ) + 1 ≤ e' := by exact_mod_cast h
      rw [abs_sub_comm, abs_of_nonneg (by linarith)]; linarith
    · have : (e' : ℝ) + 1 ≤ e := by exact_mod_cast h
      rw [abs_of_nonneg (by linarith)]; linarith

/-! ## The staircase -/

lemma isClosed_endReg (γ σ w y : ℝ) : IsClosed (endReg γ σ w y) := by
  have cV : Continuous fun p : ℝ × ℝ => γ * p.2 - σ * p.1 := by fun_prop
  exact (isClosed_le continuous_const continuous_fst).inter ((isClosed_le continuous_fst continuous_const).inter
    ((isClosed_le continuous_const continuous_snd).inter (isClosed_le cV continuous_const)))

/-- The pieces of the staircase, indexed by `Fin I ⊕ Fin J₁ ⊕ Fin R₀`. -/
def stairSet (ca sa w : ℝ) (I J₁ R₀ : ℕ) (x J : ℕ → ℝ) (n : ℕ → ℕ) (L l : ℕ → ℕ) :
    Fin I ⊕ Fin J₁ ⊕ Fin R₀ → Set (ℝ × ℝ)
  | .inl i => colSet ca sa (x i) (J i) (n i)
  | .inr (.inl j) => box (w - L j) (w - L j + L j) j (j + 1)
  | .inr (.inr r) => box 0 (0 + l r) r (r + 1)

/-- **The staircase** (Ryu's Lemma 6.3, packing part).  Columns `i < I` (base `(x_i, J_i)`, `n_i` squares at angle
`α`), layers `j < J₁` (`L_j` squares `[w − L_j, w] × [j, j + 1]`) and wedge rows `r < R₀` (`l_r` squares
`[0, l_r] × [r, r + 1]`) pack `endReg γ σ w y` under these conditions:

* columns: `x_i − n_i sa ≥ 0`, `x_i + ca ≤ w`, `J_i ≥ 0`, under the ceiling at their top-left extreme, and in
  increasing `U`-slabs `[U_i, U_i + 1]`, `U_i = ca x_i + sa J_i`;
* layers and wedge rows: inside the strip and under the ceiling;
* a layer is right of a column's slab or below the column; wedge rows are left of every column and layer. -/
theorem staircase {γ σ w y : ℝ} (hγ : 0 ≤ γ) (hσ : 0 ≤ σ) {α : ℝ} (hca : 0 ≤ Real.cos α) (hsa : 0 ≤ Real.sin α)
    (I J₁ R₀ : ℕ) (x : ℕ → ℝ) (J n L l : ℕ → ℕ)
    (hx0 : ∀ i < I, 0 ≤ x i - n i * Real.sin α) (hxw : ∀ i < I, x i + Real.cos α ≤ w)
    (hcc : ∀ i < I, γ * (J i + Real.sin α + n i * Real.cos α) - σ * (x i - n i * Real.sin α) ≤ γ * y - σ ^ 2)
    (hslab : ∀ i < I, ∀ i' < I, i < i' →
      Real.cos α * x i + Real.sin α * J i + 1 ≤ Real.cos α * x i' + Real.sin α * J i')
    (hLw : ∀ j < J₁, (L j : ℝ) ≤ w) (hLc : ∀ j < J₁, γ * (j + 1) - σ * (w - L j) ≤ γ * y - σ ^ 2)
    (hlw : ∀ r < R₀, (l r : ℝ) ≤ w) (hlc : ∀ r < R₀, γ * (r + 1) ≤ γ * y - σ ^ 2)
    (hLC : ∀ i < I, ∀ j < J₁, Real.cos α * x i + Real.sin α * J i + 1 ≤ Real.cos α * (w - L j) + Real.sin α * j ∨
      j + 1 ≤ J i)
    (hlC : ∀ i < I, ∀ r < R₀, Real.cos α * l r + Real.sin α * (r + 1) ≤ Real.cos α * x i + Real.sin α * J i)
    (hlL : ∀ j < J₁, ∀ r < R₀, Real.cos α * l r + Real.sin α * (r + 1) ≤ Real.cos α * (w - L j) + Real.sin α * j) :
    PacksIn ((∑ i ∈ Finset.range I, n i) + (∑ j ∈ Finset.range J₁, L j) + ∑ r ∈ Finset.range R₀, l r)
      (endReg γ σ w y) := by
  set ca := Real.cos α with hcad
  set sa := Real.sin α with hsad
  have hu : ca ^ 2 + sa ^ 2 = 1 := by rw [hcad, hsad, add_comm]; exact Real.sin_sq_add_cos_sq α
  have hαβ : ca ^ 2 + sa ^ 2 ≠ 0 := by rw [hu]; norm_num
  set S := stairSet ca sa w I J₁ R₀ (fun i => (x i)) (fun i => (J i : ℝ)) n L l with hS
  set N : Fin I ⊕ Fin J₁ ⊕ Fin R₀ → ℕ := fun k => match k with
    | .inl i => n i | .inr (.inl j) => L j | .inr (.inr r) => l r with hN
  have hsum : ∑ k, N k = (∑ i ∈ Finset.range I, n i) + (∑ j ∈ Finset.range J₁, L j) + ∑ r ∈ Finset.range R₀, l r := by
    rw [Fintype.sum_sum_type, Fintype.sum_sum_type, ← add_assoc]
    simp only [hN]
    rw [Fin.sum_univ_eq_sum_range (fun i => n i), Fin.sum_univ_eq_sum_range (fun j => L j),
      Fin.sum_univ_eq_sum_range (fun r => l r)]
  -- each piece is a packing of its region
  have hp : ∀ k ∈ (Finset.univ : Finset (Fin I ⊕ Fin J₁ ⊕ Fin R₀)), PacksIn (N k) (S k) := by
    rintro (i | j | r) -
    · exact column hca hsa (x i) (J i) (n i)
    · simpa [hS, hN, stairSet] using (packsIn_block (L j) 1).toBox (w - L j) j
    · simpa [hS, hN, stairSet] using (packsIn_block (l r) 1).toBox 0 r
  have hc : ∀ k ∈ (Finset.univ : Finset (Fin I ⊕ Fin J₁ ⊕ Fin R₀)), IsClosed (S k) := by
    rintro (i | j | r) -
    · exact isClosed_colSet _ _ _ _ _
    · exact isClosed_box _ _ _ _
    · exact isClosed_box _ _ _ _
  -- separations
  have colU : ∀ i : Fin I, S (.inl i) ⊆ {p | ca * p.1 + sa * p.2 ≤ ca * x i + sa * J i + 1} ∩
      {p | ca * x i + sa * J i ≤ ca * p.1 + sa * p.2} ∩ {p | (J i : ℝ) ≤ p.2} := by
    intro i p hp
    obtain ⟨h1, h2, h3, -⟩ := hp
    exact ⟨⟨h2, h1⟩, h3⟩
  have layU : ∀ j : Fin J₁, S (.inr (.inl j)) ⊆ {p | ca * (w - L j) + sa * j ≤ ca * p.1 + sa * p.2} ∩
      {p | p.2 ≤ j + 1} ∩ {p | (j : ℝ) ≤ p.2} := by
    intro j p hp
    simp only [hS, stairSet, box, mem_prod, mem_Icc] at hp
    obtain ⟨⟨h1, -⟩, h3, h4⟩ := hp
    refine ⟨⟨?_, h4⟩, h3⟩
    simp only [mem_ofPred_eq]
    nlinarith [mul_le_mul_of_nonneg_left h1 hca, mul_le_mul_of_nonneg_left h3 hsa]
  have wedU : ∀ r : Fin R₀, S (.inr (.inr r)) ⊆ {p | ca * p.1 + sa * p.2 ≤ ca * l r + sa * (r + 1)} ∩
      {p | p.2 ≤ r + 1} ∩ {p | (r : ℝ) ≤ p.2} := by
    intro r p hp
    simp only [hS, stairSet, box, mem_prod, mem_Icc, zero_add] at hp
    obtain ⟨⟨-, h2⟩, h3, h4⟩ := hp
    refine ⟨⟨?_, h4⟩, h3⟩
    simp only [mem_ofPred_eq]
    nlinarith [mul_le_mul_of_nonneg_left h2 hca, mul_le_mul_of_nonneg_left h4 hsa]
  have dU : ∀ {A B : Set (ℝ × ℝ)} {c : ℝ}, A ⊆ {p | ca * p.1 + sa * p.2 ≤ c} → B ⊆ {p | c ≤ ca * p.1 + sa * p.2} →
      Disjoint (interior A) (interior B) := fun hA hB => halfplane_disjoint hαβ hA hB
  have dY : ∀ {A B : Set (ℝ × ℝ)} {c : ℝ}, A ⊆ {p | p.2 ≤ c} → B ⊆ {p | c ≤ p.2} →
      Disjoint (interior A) (interior B) := fun hA hB => halfplane_disjoint (α := 0) (β := 1) (by norm_num)
        (fun p hp => by simp only [mem_ofPred_eq, zero_mul, zero_add, one_mul]; exact hA hp)
        (fun p hp => by simp only [mem_ofPred_eq, zero_mul, zero_add, one_mul]; exact hB hp)
  have natsep : ∀ {a b : ℕ}, a < b → (a : ℝ) + 1 ≤ b := fun h => by exact_mod_cast h
  have hd : ∀ k ∈ (Finset.univ : Finset (Fin I ⊕ Fin J₁ ⊕ Fin R₀)), ∀ k' ∈ Finset.univ, k ≠ k' →
      Disjoint (interior (S k)) (interior (S k')) := by
    -- the separations, one ordered pair at a time
    have cc : ∀ i i' : Fin I, (i : ℕ) < i' → Disjoint (interior (S (.inl i))) (interior (S (.inl i'))) :=
      fun i i' h => dU (fun p hp => (colU i hp).1.1)
        (fun p hp => le_trans (hslab i i.isLt i' i'.isLt h) (colU i' hp).1.2)
    have cL : ∀ (i : Fin I) (j : Fin J₁), Disjoint (interior (S (.inl i))) (interior (S (.inr (.inl j)))) := by
      intro i j
      rcases hLC i i.isLt j j.isLt with h | h
      · exact dU (fun p hp => (colU i hp).1.1) (fun p hp => le_trans h (layU j hp).1.1)
      · have h' : ((j : ℕ) : ℝ) + 1 ≤ (J i : ℝ) := by exact_mod_cast h
        exact (dY (fun p hp => (layU j hp).1.2) (fun p hp => show ((j : ℕ) : ℝ) + 1 ≤ p.2 from
          le_trans h' (colU i hp).2)).symm
    have cW : ∀ (i : Fin I) (r : Fin R₀), Disjoint (interior (S (.inl i))) (interior (S (.inr (.inr r)))) :=
      fun i r => (dU (fun p hp => (wedU r hp).1.1)
        (fun p hp => le_trans (hlC i i.isLt r r.isLt) (colU i hp).1.2)).symm
    have LL : ∀ j j' : Fin J₁, (j : ℕ) < j' →
        Disjoint (interior (S (.inr (.inl j)))) (interior (S (.inr (.inl j')))) :=
      fun j j' h => dY (fun p hp => (layU j hp).1.2) (fun p hp => le_trans (natsep h) (layU j' hp).2)
    have LW : ∀ (j : Fin J₁) (r : Fin R₀),
        Disjoint (interior (S (.inr (.inl j)))) (interior (S (.inr (.inr r)))) :=
      fun j r => (dU (fun p hp => (wedU r hp).1.1)
        (fun p hp => le_trans (hlL j j.isLt r r.isLt) (layU j hp).1.1)).symm
    have WW : ∀ r r' : Fin R₀, (r : ℕ) < r' →
        Disjoint (interior (S (.inr (.inr r)))) (interior (S (.inr (.inr r')))) :=
      fun r r' h => dY (fun p hp => (wedU r hp).1.2) (fun p hp => le_trans (natsep h) (wedU r' hp).2)
    have ne_lt : ∀ {m : ℕ} (a b : Fin m), a ≠ b → (a : ℕ) < b ∨ (b : ℕ) < a := fun a b h =>
      Nat.lt_or_gt_of_ne (fun e => h (Fin.ext e))
    rintro (i | j | r) - (i' | j' | r') - hne
    · rcases ne_lt i i' (fun e => hne (by rw [e])) with h | h
      · exact cc i i' h
      · exact (cc i' i h).symm
    · exact cL i j'
    · exact cW i r'
    · exact (cL i' j).symm
    · rcases ne_lt j j' (fun e => hne (by rw [e])) with h | h
      · exact LL j j' h
      · exact (LL j' j h).symm
    · exact LW j r'
    · exact (cW i' r).symm
    · exact (LW j' r).symm
    · rcases ne_lt r r' (fun e => hne (by rw [e])) with h | h
      · exact WW r r' h
      · exact (WW r' r h).symm
  have P := PacksIn.biUnion Finset.univ S N hc hp hd
  rw [hsum] at P
  refine P.mono ?_
  -- every piece lies in the end region
  simp only [Finset.mem_univ, iUnion_true]
  refine iUnion_subset ?_
  rintro (i | j | r) p hp
  · obtain ⟨-, -, h3, h4, h5, h6⟩ := hp
    have a0 := hx0 i i.isLt
    have a1 := hxw i i.isLt
    have a2 := hcc i i.isLt
    refine ⟨by linarith, by linarith, by linarith [(Nat.cast_nonneg (J i) : (0 : ℝ) ≤ J i)], ?_⟩
    nlinarith [mul_le_mul_of_nonneg_left h4 hγ, mul_le_mul_of_nonneg_left h5 hσ]
  · simp only [hS, stairSet, box, mem_prod, mem_Icc] at hp
    obtain ⟨⟨h1, h2⟩, h3, h4⟩ := hp
    have a0 := hLw j j.isLt
    have a1 := hLc j j.isLt
    refine ⟨by linarith, by linarith, by linarith [(Nat.cast_nonneg (j : ℕ) : (0 : ℝ) ≤ (j : ℕ))], ?_⟩
    nlinarith [mul_le_mul_of_nonneg_left h4 hγ, mul_le_mul_of_nonneg_left h1 hσ]
  · simp only [hS, stairSet, box, mem_prod, mem_Icc, zero_add] at hp
    obtain ⟨⟨h1, h2⟩, h3, h4⟩ := hp
    have a0 := hlw r r.isLt
    have a1 := hlc r r.isLt
    refine ⟨h1, by linarith, by linarith [(Nat.cast_nonneg (r : ℕ) : (0 : ℝ) ≤ (r : ℕ))], ?_⟩
    nlinarith [mul_le_mul_of_nonneg_left h4 hγ, mul_le_mul_of_nonneg_left h1 hσ]

end UnitSquarePacking.Ryu
