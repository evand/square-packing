import Sqpack.StripRect

/-!
# The (0,8) strip leaf for `s(12) = 4`: statement and stitch

`notes/strip-split.md` §4–5.  Twelve closed unit squares in `box t`, `t ≤ 4`, pairwise disjoint
as closed sets.  Four of them (the *strip*, indices `s`) sit against the left wall.  The other eight
are to their right.  The leaf is closed by two statements, kept here as named `Prop`s so that
nothing is assumed silently:

* `SeamOK` (4-square side): every unit square in the container that is disjoint from the four strip
  squares lies in the slab `[x0, t] × [0, t]`.  Measured: holds with `x0 = 1` for classes A, B, F
  (`notes/strip-split.md` §4a–4c); it is the "optimistic seam".
* `RegionOK` (8-square side): no eight closed unit squares, pairwise disjoint, lie in that slab and
  avoid one given square (the strip's top square).  Measured: holds at `x0 = 1`, `t = 3.99` for the
  top squares of classes A, B, F.

`strip_leaf` is the stitch: the two statements and the packing are contradictory.  The proof is
bookkeeping; the point is that the hypotheses are exactly these and nothing else.  The class cover
("every strip of width `< w1` is in class A, B, F, or goes to another branch") and the proofs of
`SeamOK` / `RegionOK` per class are the remaining obligations.
-/

open Finset
open scoped Classical

namespace SquarePacking

/-- Every unit square in `box t` that is closed-disjoint from the four given squares lies in the
slab `[x0, t] × [0, t]`. -/
def SeamOK (t x0 : ℝ) (sc : Fin 4 → ℝ × ℝ) (sθ : Fin 4 → ℝ) : Prop :=
  ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ box t → (∀ k, Disjoint (sq c θ 1) (sq (sc k) (sθ k) 1)) →
    sq c θ 1 ⊆ vslab x0 t t

/-- No eight closed unit squares, pairwise disjoint, lie in `[x0, t] × [0, t]` and avoid the square
`sq tc tθ 1`. -/
def RegionOK (t x0 : ℝ) (tc : ℝ × ℝ) (tθ : ℝ) : Prop :=
  ∀ (c : Fin 8 → ℝ × ℝ) (θ : Fin 8 → ℝ),
    (∀ i, sq (c i) (θ i) 1 ⊆ vslab x0 t t) →
    (∀ i, Disjoint (sq (c i) (θ i) 1) (sq tc tθ 1)) →
    (∀ i j, i ≠ j → Disjoint (sq (c i) (θ i) 1) (sq (c j) (θ j) 1)) → False

/-- **The stitch.**  Twelve closed-disjoint unit squares in `box t`; an embedding `e` picks the four
strip squares and `k` names the strip's top square.  If the strip satisfies `SeamOK` and the region
to its right satisfies `RegionOK` for that top square, the packing does not exist. -/
theorem strip_leaf {t x0 : ℝ} (ctr : Fin 12 → ℝ × ℝ) (ang : Fin 12 → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) 1 ⊆ box t)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1))
    (e : Fin 4 ↪ Fin 12) (k : Fin 4)
    (hseam : SeamOK t x0 (fun m => ctr (e m)) (fun m => ang (e m)))
    (hregion : RegionOK t x0 (ctr (e k)) (ang (e k))) : False := by
  -- the eight squares outside the strip, in some order
  set R : Finset (Fin 12) := (univ.map e)ᶜ with hR
  have hcard : R.card = 8 := by
    rw [hR, card_compl, card_map, card_univ, Fintype.card_fin, Fintype.card_fin]
  let f : Fin 8 ↪ Fin 12 :=
    ⟨fun i => (R.equivFin.symm (Fin.cast hcard.symm i) : Fin 12), by
      intro a b hab
      have := Subtype.ext hab
      simpa using this⟩
  have hfR : ∀ i, f i ∈ R := fun i => (R.equivFin.symm (Fin.cast hcard.symm i)).2
  have hfnot : ∀ i m, f i ≠ e m := by
    intro i m h
    have := hfR i
    rw [hR, mem_compl, mem_map] at this
    exact this ⟨m, mem_univ _, h.symm⟩
  refine hregion (fun i => ctr (f i)) (fun i => ang (f i)) ?_ ?_ ?_
  · intro i
    exact hseam _ _ (hin (f i)) (fun m => hdisj _ _ (hfnot i m))
  · intro i
    exact hdisj _ _ (hfnot i k)
  · intro i j hij
    exact hdisj _ _ (fun h => hij (f.injective h))

end SquarePacking
