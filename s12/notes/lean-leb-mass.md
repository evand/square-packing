# Lean: uniform-area (Lebesgue) mass primitives, Lemma U and Lemma K (2026-10-03)

Task `lean-leb-mass` (brief: `tasks/lean-leb-mass/README.md`).  **Status: done** — Lemma U, Lemma K (with the
width lemma and the tangent-cap refinement) proved sound, and all 689 LEB + 374 CAP leaves of run V3 checked in the
kernel (§6).  §1–§5 are the design as written before the code; the code follows it (deviations in §6).  Math: `search/QUADRANT_EXACT.md` §4.2
(LEB = Lemma U, CAP = Lemma K), §4.3 (chord ends); the tests as run: `search/qx2_zm.py` `cert_leb`, `cert_cap`
(~l. 939–1010).  Lean file: `lean/Sqpack/LebMass.lean` (generic), `lean/Sqpack/LebMass7.lean` (smoke test).

## 1. Where the leaves live: `CovT`, a box predicate on the measure

`ValidTilt m μ` (`ValidSplit.lean`) is stated on a measure, and `ValidTilt7 = ValidTilt 7 box7Cover.measure`, with
`box7_grid : box7Cover.measure = (gridCover 7 9 10¹² (wP box7Packed 35)).measure` already proved.  The existing tree
predicates do not fit as they are: `BoxTree.Cov` / `ZMTreeM.CovM` count point and segment mass over a list (no
polygon), on the `ℕ`-scaled box `[x0/(D·S), …] × [u0/R, u1/R]` with `u ∈ [0, 29/70]`.  So the leaf statements here
use a new predicate on **any measure**, in exactly `ValidTilt`'s pose region:

```lean
def CovT (m : ℝ) (μ : Measure (ℝ × ℝ)) (x0 x1 y0 y1 u0 u1 : ℝ) : Prop :=
  ∀ (c : ℝ × ℝ) (u : ℝ), x0 ≤ c.1 → c.1 ≤ x1 → y0 ≤ c.2 → c.2 ≤ y1 → u0 ≤ u → u ≤ u1 →
    0 < u → u ^ 2 + 2 * u ≤ 1 → sq c (2 * Real.arctan u) 1 ⊆ box m → 1 ≤ μ (sq c (2 * Real.arctan u) 1)
```

with `CovT.splitX/Y/U` (gluing) and `validTilt_of_covT : CovT m μ 0 (m/2) 0 (m/2) 0 (1/2) → ValidTilt m μ` (the
root: `u² + 2u ≤ 1 ⇒ u ≤ 1/2`).  A leaf box is a `PBox` of six rationals (the Python leaves are `Fraction`s; a future
`ℕ`-scaled tree passes `x/Q`, `u/R` as rationals at the leaf).  Points and segments are not lost: a later `CovM`-style
certificate lifts to `CovT` of `gridCover` by monotonicity (gridCover = segments + Lebesgue part ≥ the segments).

**Leaf tests in `ℚ`, by `decide +kernel`.**  The run has only 689 LEB and 374 CAP leaves (run V3), so per-leaf cost
does not matter (measured: ~4 ms per rational inequality in the kernel); `ℚ` keeps the tests readable and literally
the Python ones.  No `native_decide`.

**The angle bin and `ŵ`.**  `w(u) = cos θ + sin θ = widU u` (`ZeroMargin.lean`).  Bound used (zeromargin's
`bin_data` `whi`, specialised to what `ValidTilt` admits):
`whi u1 = if u1² + 2u1 ≤ 1 then widU u1 else 14143/10000`; sound for every `u ∈ (0, u1]` with `u² + 2u ≤ 1`: `widU`
is increasing on `[0, √2 − 1]` (`widU_strictMonoOn`, existing), and `widU ≤ √2 < 1.4143` always.  (`bin_data` takes
`max(w(u0), w(u1))` on a bin entirely past 45°; such bins are `SYM` leaves, never LEB/CAP of `ValidTilt`.)  Every
square then lies in `[cx − ŵ/2, cx + ŵ/2] × [cy − ŵ/2, cy + ŵ/2]` (`abs_sub_le_wid`, existing).

## 2. Lemma U (LEB)

```lean
def lebOK (a b : ℚ) (B : PBox) : Bool    -- a ≤ x0 − ŵ/2, a ≤ y0 − ŵ/2, x1 + ŵ/2 ≤ b, y1 + ŵ/2 ≤ b
theorem leb_sound {m : ℝ} {μ : Measure (ℝ × ℝ)} {a b : ℚ}
    (hμ : ∀ S, MeasurableSet S → volume (S ∩ Set.Icc (a:ℝ) b ×ˢ Set.Icc (a:ℝ) b) ≤ μ S)
    {B : PBox} (h : lebOK a b B = true) : B.Cov m μ
theorem gridCover_vol_le (hA : 2 * A ≤ 5 * K) (S) (hS : MeasurableSet S) :
    volume (S ∩ BentzFam.lebSq A K) ≤ (gridCover K A den w).measure S
```

New: **`volume (sq c θ 1) = 1`** (the rotated square is the preimage of `[−½,½]²` under `coord c θ`, an affine map
with linear part of determinant `cos² + sin² = 1`; via `MeasurableEquiv.finTwoArrow`, `Matrix.toLin'` and
`addHaar_preimage_linearMap`), and the polygon term of `gridCover` = `volume.restrict lebSq` (its weight is its area).

## 3. Lemma K (CAP)

Pose-level statement (θ = 2 arctan u, `0 < u`, `u² + 2u ≤ 1`, `Q = sq c θ 1 ⊆ box K`, `U = [a,b]²`, `a = A/5`,
`b = K − A/5`): if `Q ⊆ {x ≤ b, y ≤ b}` and, for the line `y = a`, either `Q ⊆ {y ≥ a}` or (`a ≤ cy`, `Q ⊆ {y ≥ a −
d}`, and every unit grid segment `(false, i, A)` with `i` in a cell range covering the chord has `d ≤ 5 w/den`), and
likewise for `x = a` (segments `(true, A, j)`), then `μ(Q) ≥ 1`.  Proof, as QUADRANT_EXACT §4.2:

1. `1 = λ(Q) ≤ λ(Q ∩ U) + λ(Q ∩ {y < a}) + λ(Q ∩ {x < a})`.
2. **Width lemma** (new, generic): `Q` convex, closed, bounded, centrally symmetric about `c`, `y ≤ a ≤ cy` ⇒
   `λ₁(Q_y) ≤ λ₁(Q_a)` for the horizontal slices `Q_y = {x | (x,y) ∈ Q}`.  Proof: `Q_y = [l, r]` (compact, convex),
   `Q_{2cy−y} = 2cx − Q_y`, and with `λ = (cy − y + cy − a)/(2(cy − y)) ∈ [½, 1]`, the two points
   `λ l + (1−λ)(2cx − r)`, `λ r + (1−λ)(2cx − l)` are in `Q_a` (convexity), at distance `r − l`.  Then by
   `Measure.prod_apply_symm` (Fubini on `volume = volume.prod volume`): `λ(Q ∩ {y < a}) ≤ d · λ₁(Q_a)` when
   `Q ⊆ {y ≥ a − d}`.  The `x = a` case is the same lemma for `Prod.swap ⁻¹' Q`.
3. **Segments on the line** (new): for `s = (false, i, A)`, `segFrac s Q = 5 λ₁(Q_a ∩ [i/5, (i+1)/5])` (an affine
   change of variables), so `Σ_i (w_i/den) segFrac ≥ ρ λ₁(Q_a)` when the cells cover `Q_a` and `ρ ≤ 5 w_i/den`.
4. `μ(Q) ≥ λ(Q ∩ U) + ρ_y λ₁(Q_a) + ρ_x λ₁(Q^x_a) ≥ 1 − d_y λ₁ + ρ_y λ₁ − … ≥ 1`.

Chord range, as `cert_cap`: crude `[x0 − ŵ/2, x1 + ŵ/2]`; **tangent cap**: `Q ⊆ BL + cone{(c,s), (−s,c)}` (BL the
lowest vertex, `(cx − (c−s)/2, cy − (c+s)/2)`), so the chord at depth `t = a − y_BL ≤ d` lies in
`[x_BL − (s/c)t, x_BL + (c/s)t]`, bounded over the box by `[x0 − (c0−s0)/2 − (s1/c1)d, x1 − (c1−s1)/2 + (c0/s0)d]`
(`c − s` decreasing, `s/c` increasing in `u` on `[0,1]`).  This holds for every `t ≥ 0`, so the kernel intersects
the two ranges whenever `u0 > 0` (`cert_cap` also asks `d ≤ min(s0, c1)`, which is only needed for tightness).  For
`x = a`: TL vertex, `[y0 + (c1−s1)/2 − (c0/s0)d, y1 + (c0−s0)/2 + (s1/c1)d]`.  Cell indices: `⌊5 t0⌋ … ⌈5 t1⌉ − 1`,
clipped to the grid (`Q ⊆ box K` bounds the chord anyway).

```lean
def capOK (K A den : ℕ) (w : SegIx → ℕ) (B : PBox) : Bool
theorem cap_sound (hA : 2 * A < 5 * K) (h : capOK K A den w B = true) :
    B.Cov K (gridCover K A den w).measure
```

## 4. Smoke test

`LebMass7.lean`: every LEB and CAP leaf of run V3 (`qx2_data/cert/runV3_6294052a_leaves.jsonl.gz`: 689 + 374
boxes, after `clip_bin`), as literal `PBox` lists, checked by `decide +kernel` against the packed `μ₇` weights
(`wP box7Packed 35`), giving `CovT 7 box7Cover.measure` for each.  Generator `lean/scripts/gen_lebmass_data.py`.

## 5. Size (estimate before writing; actual in §6)

| piece | new / reused | est. lines |
|---|---|---|
| `CovT`, splits, root | new (trivial) | 60 |
| `ŵ` bound | reuses `widU_strictMonoOn`, `wid_two_arctan`, `abs_sub_le_wid` | 50 |
| `λ(sq) = 1`, `gridCover ≥ λ|_U` | new | 80 |
| Lemma U leaf | new | 60 |
| width lemma + Fubini | new | 200 |
| segment change of variables, cells | new (reuses `segFrac`, `segMeasure_apply`) | 120 |
| tangent cap, ranges, `capOK` soundness | new | 200 |
| smoke test + generator | new | 100 + data |

About 900 lines: one to two sessions.  The segment lemmas of batch 1 (S/T/L/P) are not needed here: U and K only use
the Lebesgue part and the two boundary lines `x = a`, `y = a`, by exact slicing; the chord-end formulas of §4.3 are
replaced by the cone bound (the same numbers, with a containment proof instead of four options).

## 6. Result (what was built)

Files (default build): `lean/Sqpack/LebMass.lean` (951 lines), `lean/Sqpack/LebMass7Data.lean` (generated, the
leaf boxes), `lean/Sqpack/LebMass7.lean` (smoke test), `lean/scripts/gen_lebmass_data.py` (generator + exact Python
mirror of `lebOK` / `capOK`).  Every theorem below prints `[propext, Classical.choice, Quot.sound]`
(`lean/Axioms.lean`); no `sorry`, no `native_decide`.

```lean
-- LebMass.lean (namespace SquarePacking.LebMass)
def CovT (m : ℝ) (μ : Measure (ℝ × ℝ)) (x0 x1 y0 y1 u0 u1 : ℝ) : Prop      -- §1, + splitX/Y/U
theorem validTilt_of_covT (h : CovT m μ 0 (m / 2) 0 (m / 2) 0 (1 / 2)) : ValidTilt m μ
theorem volume_sq (c : ℝ × ℝ) (θ : ℝ) : volume (sq c θ 1) = 1
theorem leb_sound {a b : ℚ} (hμ : ∀ S, MeasurableSet S → volume (S ∩ lsq a b) ≤ μ S)
    {B : PBox} (h : lebOK a b B = true) : B.Cov m μ                                       -- Lemma U
theorem gridCover_vol_le (hA : 2 * A < 5 * K) (S) (hS : MeasurableSet S) :
    volume (S ∩ BentzFam.lebSq A K) ≤ (gridCover K A den w).measure S
theorem vol_hsl_le (hconv : Convex ℝ Q) (hcomp : IsCompact Q)
    (hsym : ∀ p ∈ Q, (2 * c.1 - p.1, 2 * c.2 - p.2) ∈ Q) (hya : y ≤ a) (hac : a ≤ c.2) :
    volume (hsl Q y) ≤ volume (hsl Q a)                                                  -- width lemma
theorem sq_below_le (c) (θ) (hac : a ≤ c.2) (hb : ∀ p ∈ sq c θ 1, a - d ≤ p.2) :
    volume (sq c θ 1 ∩ {p | p.2 < a}) ≤ ENNReal.ofReal d * volume (hsl (sq c θ 1) a)
theorem cap_sound (hA : 2 * A < 5 * K) (hden : 0 < den) {B : PBox}
    (h : capOK K A den w B = true) : B.Cov K (gridCover K A den w).measure                -- Lemma K
-- LebMass7.lean (namespace SquarePacking.Bentz)
theorem leb7_cov : ∀ B ∈ lebLeaves7, B.Cov 7 box7Cover.measure
theorem cap7_cov : ∀ B ∈ capLeaves7, B.Cov 7 box7Cover.measure
```

**Deviations from §1–§4.**  None of substance.  (i) `capOK` checks the cells one by one (`d·den ≤ 5 w_i` for each
cell meeting the chord range) instead of taking the minimum density; same test.  (ii) The tangent-cap range is
intersected with the crude one whenever `0 < u0`, `u1 < 1` (the cone bound holds at every depth), where `cert_cap` asks
`d ≤ min(sin θ₀, cos θ₁)` first; so `capOK` is at least as strong.  (iii) The chord-end table of §4.3 is not used:
the tangent cap is proved from the cone inclusions `C(p₁ − x_BL) + S(p₂ − y_BL) ≥ 0`, `S(p₁ − x_BL) − C(p₂ − y_BL)
≤ 0` (`cone_BL`; `cone_TL` for `x = a`), each an identity `= X + ½`, `= −Y − ½` in the rotated coordinates.
(iv) The `x = a` cap reuses the width lemma through `swapXY` (`sq_left_le`), the segment lemma is proved for both
orientations (`segMeasure_h`, `segMeasure_v`: a unit grid segment captures `5 λ₁(slice ∩ cell)`).  (v) The whole
final inequality is in `ℝ≥0∞`, with no subtraction: `1 = λ(Q) ≤ λ(Q∩U) + λ(Q∩{y<a}) + λ(Q∩{x<a}) ≤ λ(Q∩U) +
Σ_cells(line y = a) + Σ_cells(line x = a) ≤ μ(Q)`.

**Smoke test.**  `leb7_ok`, `cap7_ok` (`decide +kernel`) run the tests on all 689 LEB and 374 CAP leaves of run V3
(the full census of those kinds, not a sample) against the packed `μ₇` weights; `leb7_cov` / `cap7_cov` give `CovT`
of each box for `box7Cover.measure` via `box7_grid`.  The Python mirror rejects none (and the kernel agrees).
Mutation checks (`#eval`, not committed): with all weights 0, `capOK` rejects all 374 CAP leaves; without the tangent
refinement (`u0 := 0`) it rejects 124 of them (so the refinement is exercised); with every weight scaled by 0.99 it
still accepts all (the line densities are far above the cap depths at these leaves).  No CAP leaf passes `lebOK`.

**Timings** (`LEAN_NUM_THREADS=4`, no pinning, load ~3): `LebMass` 12 s, `LebMass7Data` 6 s, `LebMass7` 10 s (both
kernel checks together a few seconds; the file is dominated by the Mathlib import); full default `lake build` after
the change: 34 s wall (everything else cached).

## 7. What is left for `ValidTilt7` (estimate)

Run V3's leaves: PIECE 15,926, EXACT 9,763 + EXACT0 381 + EXACT45 759, CAP 374, LEB 689, SYM 1,171, AXIS 61, EMPTY 2,955
(no ADM/P1/MIX/SPLIT: the cover has no points).  Done in Lean: LEB, CAP (here), AXIS (Lemma Z, `validAxis7`),
SYM/EMPTY/clip (`ValidSplit`, `clip_bin_no_loss`, `sq_subset_box_iff`).  Missing:

1. **The tree layer** (`CovT` tree with `X/Y/U/XM/YM/UM` splits, clip and SYM nodes, leaf dispatch; generator from the
   leaf dump): ~400 lines + generator, mechanical; the dump already records every leaf box.
2. **PIECE with the polygon** (zm_mixed's `piece_bound`: Lemma S/T/L on the segments, which `ZMTreeM`/`LBlock` have,
   plus **Lemma S(b)** for the Lebesgue square: area of `U ∩ K`, `K` the intersection of certified half-planes, a
   convex-polygon clip in the kernel and its soundness), and the bridge from `CovM`-style segment mass to the grid
   cover's measure.  Also L′/V if the census needs them (`notes/lean-segments.md` §5).  Estimate 1–2 k lines, 3–5
   sessions; kernel ~16 k leaves × ~0.1–0.5 s.
3. **Lemma E** (11 k leaves, the bulk of the Python CPU): the E′/E″ concave minorants of `λ(Q ∩ U)` (cap areas
   `g(d)`, the corner3/xcut inclusion–exclusion, McCormick), the line-arrangement argument (piecewise concavity on
   cells, minimum at vertices, continuity), and the vertex certificates (rational vertices `X/Δ`, sign of `Δ`,
   general-degree Bernstein, S-procedure, every combination of alternatives).  The width lemma and slices here help
   with E′, but the arrangement/vertex part is new and large.  A Lean-friendly reformulation (the generator emits, per
   leaf, the explicit cells or a per-vertex certificate the kernel only checks) is advisable.  Estimate 3–6 k lines,
   **2–4 weeks of sessions**, and a kernel run likely tens of CPU-hours.

So `ValidTilt7` in Lean is **~5–8 k lines beyond this batch, several weeks**; `ValidTilt9` reuses everything (a
second run's leaf dump and build).  The order: tree layer + PIECE first (closes 17 k of 32 k leaves with what
exists), then Lemma E.
