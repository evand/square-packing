import Sqpack.ExTree
import Sqpack.LebMass

/-!
# The pose-tree boxes of the Valid7 run and `LebMass`

Correspondences between the ℕ-encoded checker of the Valid7 kernel run (`ExTree.CovT` over the mixed
cover `mcoverP`, with its Lebesgue rectangles `rectMass`) and the real-valued `LebMass` (`LebMass.CovT`,
Lemma U `LebMass.leb_sound`).  No module imported by the generated root modules imports this file, so
adding it changes no root.

* `measure_mcoverP`, `one_le_measure_mcoverP_iff`: the measure of a square is its mass `/ W`.
* `CovT.toLebMass`, `CovT.ofLebMass`: an `ExTree.CovT` node is a `LebMass.CovT` box of the same measure
  (`ExTree.CovT` has an open lower end in `u`; `LebMass.CovT` a closed one and `u > 0`, so the first
  direction needs the lower end `0`).
* `volume_inter_le_mcoverP`: a rectangle of `mcoverP` with density `≥ 1` (`W ≤ w`) is `≥ λ` on any
  square `[a,b]²` inside it, which is the hypothesis of `LebMass.leb_sound`; `CovT.of_lebOK` turns the
  `LebMass` LEB test into an `ExTree.CovT` leaf of our cover.
-/

open MeasureTheory

namespace SquarePacking

namespace ZMTreeM

open BoxTree ZMTree

variable {D S Mq R W : ℕ} {pts : PTree} {segs : STree} {rects : List RectE}

lemma mcoverP_nonneg : (mcoverP D W pts segs rects).Nonneg :=
  ⟨fun _ _ => by simp only [mcoverP]; positivity,
   fun _ _ => by simp only [mcoverP]; positivity,
   fun _ _ => by simp only [mcoverP]; positivity⟩

/-- The measure of a measurable set is its mass `/ W` (as in `cover_of_tree`). -/
theorem measure_mcoverP {Qs : Set (ℝ × ℝ)} (hQs : MeasurableSet Qs) :
    (mcoverP D W pts segs rects).measure Qs = ENNReal.ofReal
      ((ptMass D pts.toList Qs + segMass D segs.toList Qs + rectMass D rects Qs) / W) := by
  rw [(mcoverP D W pts segs rects).measure_apply mcoverP_nonneg hQs]
  congr 1
  simp only [MixedCover.mass, mcoverP, ptMass, segMass, rectMass, gw_areaFrac]
  rw [add_div, add_div, Finset.sum_div, Finset.sum_div, Finset.sum_div]
  congr 1
  · congr 1
    refine Finset.sum_congr rfl fun e _ => ?_
    ring

theorem one_le_measure_mcoverP_iff (hW : 0 < W) {Qs : Set (ℝ × ℝ)} (hQs : MeasurableSet Qs) :
    1 ≤ (mcoverP D W pts segs rects).measure Qs ↔
      (W : ℝ) ≤ ptMass D pts.toList Qs + segMass D segs.toList Qs + rectMass D rects Qs := by
  have hWr : (0 : ℝ) < W := by exact_mod_cast hW
  rw [measure_mcoverP hQs, ENNReal.one_le_ofReal, one_le_div hWr]

/-- **`ExTree.CovT` → `LebMass.CovT`** (lower end `u = 0`). -/
theorem CovT.toLebMass (hW : 0 < W) {x0 x1 y0 y1 u1 : ℕ}
    (h : CovT D S Mq R W pts.toList segs.toList rects x0 x1 y0 y1 0 u1) :
    LebMass.CovT ((Mq : ℝ) / D) (mcoverP D W pts segs rects).measure
      ((x0 : ℝ) / (D * S)) ((x1 : ℝ) / (D * S)) ((y0 : ℝ) / (D * S)) ((y1 : ℝ) / (D * S))
      0 ((u1 : ℝ) / R) := by
  intro c u hx0 hx1 hy0 hy1 _ hu1 hu hu' hsub
  rw [one_le_measure_mcoverP_iff hW (measurableSet_sq _ _ _)]
  exact h c u hx0 hx1 hy0 hy1 (by simpa using hu) hu1 hu' hsub

/-- **`LebMass.CovT` → `ExTree.CovT`** (any lower end). -/
theorem CovT.ofLebMass (hW : 0 < W) {x0 x1 y0 y1 u0 u1 : ℕ}
    (h : LebMass.CovT ((Mq : ℝ) / D) (mcoverP D W pts segs rects).measure
      ((x0 : ℝ) / (D * S)) ((x1 : ℝ) / (D * S)) ((y0 : ℝ) / (D * S)) ((y1 : ℝ) / (D * S))
      ((u0 : ℝ) / R) ((u1 : ℝ) / R)) :
    CovT D S Mq R W pts.toList segs.toList rects x0 x1 y0 y1 u0 u1 := by
  intro c u hx0 hx1 hy0 hy1 hu0 hu1 hu' hsub
  have hu : 0 < u := lt_of_le_of_lt (by positivity) hu0
  exact (one_le_measure_mcoverP_iff hW (measurableSet_sq _ _ _)).1
    (h c u hx0 hx1 hy0 hy1 hu0.le hu1 hu hu' hsub)

/-- **Lemma U's hypothesis for `mcoverP`**: a rectangle of density `≥ 1` (`W ≤ w`) that contains
`[a,b]²` gives `λ(· ∩ [a,b]²) ≤ μ`. -/
theorem volume_inter_le_mcoverP (hW : 0 < W) {e : RectE} (he : e ∈ rects) (hw : W ≤ e.2.2.2.2)
    {a b : ℝ} (hab : LebMass.lsq a b ⊆ rectSet D e) (S' : Set (ℝ × ℝ)) (hS' : MeasurableSet S') :
    volume (S' ∩ LebMass.lsq a b) ≤ (mcoverP D W pts segs rects).measure S' := by
  classical
  have hWr : (0 : ℝ) < W := by exact_mod_cast hW
  have hfin : volume (S' ∩ rectSet D e) ≠ ⊤ :=
    ne_top_of_le_ne_top (volume_rectSet_ne_top D e) (measure_mono Set.inter_subset_right)
  set v := (volume (S' ∩ rectSet D e)).toReal with hv
  have hv0 : 0 ≤ v := ENNReal.toReal_nonneg
  have hterm : (e.2.2.2.2 : ℝ) * v ≤ rectMass D rects S' :=
    Finset.single_le_sum (f := fun e => (e.2.2.2.2 : ℝ) * (volume (S' ∩ rectSet D e)).toReal)
      (fun _ _ => mul_nonneg (Nat.cast_nonneg _) ENNReal.toReal_nonneg) (List.mem_toFinset.mpr he)
  have hwv : v ≤ (e.2.2.2.2 : ℝ) * v / W := by
    rw [le_div_iff₀ hWr]
    have : (W : ℝ) ≤ e.2.2.2.2 := by exact_mod_cast hw
    nlinarith
  have hp := ptMass_nonneg D pts.toList S'
  have hs := segMass_nonneg D segs.toList S'
  calc volume (S' ∩ LebMass.lsq a b) ≤ volume (S' ∩ rectSet D e) :=
        measure_mono (Set.inter_subset_inter_right _ hab)
    _ = ENNReal.ofReal v := (ENNReal.ofReal_toReal hfin).symm
    _ ≤ (mcoverP D W pts segs rects).measure S' := by
        rw [measure_mcoverP hS']
        apply ENNReal.ofReal_le_ofReal
        calc v ≤ (e.2.2.2.2 : ℝ) * v / W := hwv
          _ ≤ rectMass D rects S' / W := by gcongr
          _ ≤ (ptMass D pts.toList S' + segMass D segs.toList S' + rectMass D rects S') / W := by
              gcongr; linarith

/-- **The `LebMass` LEB test as a leaf of our cover**: `LebMass.lebOK` on a box, with a rectangle of
density `≥ 1` around `[a,b]²`, gives the `ExTree.CovT` node. -/
theorem CovT.of_lebOK (hW : 0 < W) {e : RectE} (he : e ∈ rects) (hw : W ≤ e.2.2.2.2) {a b : ℚ}
    (hab : LebMass.lsq (a : ℝ) b ⊆ rectSet D e) {x0 x1 y0 y1 u0 u1 : ℕ}
    (h : LebMass.lebOK a b ⟨(x0 : ℚ) / (D * S), (x1 : ℚ) / (D * S), (y0 : ℚ) / (D * S),
      (y1 : ℚ) / (D * S), (u0 : ℚ) / R, (u1 : ℚ) / R⟩ = true) :
    CovT D S Mq R W pts.toList segs.toList rects x0 x1 y0 y1 u0 u1 := by
  apply CovT.ofLebMass hW
  have := LebMass.leb_sound (m := (Mq : ℝ) / D)
    (fun S' hS' => volume_inter_le_mcoverP (pts := pts) (segs := segs) hW he hw hab S' hS') h
  simpa [LebMass.PBox.Cov] using this

end ZMTreeM

end SquarePacking
