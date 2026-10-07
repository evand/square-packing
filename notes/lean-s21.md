# `s(21) = 5` in Lean, from one computational hypothesis (`lean/Sqpack/{MixedMeasure,SegTree,S21Data,S21}.lean`)

**What this is.**  The mixed cover `certificates/s21/s21_mixed_cover_5.txt` (`certificates/s21/FORMAT.md` v1:
7,536 points and 1,872 segments, no polygons, total `2089474919732/10¹¹ = 20.89474919732 < 21`; until 2026-09-27
`runs/line-cover_m5_candidate_x1003.txt`, same bytes) is certified VERIFIED-D4 by `search/zm_mixed.py`
(`search/ZM_MIXED.md` §7.5, §8; the shipped run is `certificates/s21/zm_mixed_d4/`) and independently by `zmx2`
(`search/ZMX2.md`; `certificates/s21/zmx2_{d4,full}/`).  The Lean files state and prove `s(21) = 5` **from a single named
hypothesis, the checker's statement on the D4 region**, with everything else proved: the reduction from
*measures* (not just point sets) to packings, the D4 reduction for measures, the D4 invariance and the total of the
actual cover (as data), the scaling argument and the 5 × 5 packing.  Mathlib `v4.33.1`, `lake build` clean, no
`sorry`, standard axioms only, no `native_decide`.  Mirrors `notes/lean-s32.md`.

```lean
theorem s21_eq_five (h : S21RegionCover) : minSide 21 = 5
theorem s21_eq_five_of_checker (h : S21CheckerCover) : minSide 21 = 5
theorem s21_packs : Packs 21 5                                   -- upper bound, no hypothesis
```

`Packs`, `minSide` are those of `S32.lean` (`n` closed unit squares, any centres and angles, in `[0,s]²`, pairwise
disjoint interiors; `minSide n = sInf {s | Packs n s}`).

---

## 1.  The reduction for measures (`MixedMeasure.lean` §1–2)

```lean
theorem packing_le_measure (μ : Measure (ℝ × ℝ)) (C : Set (ℝ × ℝ))
    (hcover : ∀ c θ, sq c θ 1 ⊆ C → 1 ≤ μ (sq c θ 1))
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr ang) (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt ..) (sqInt ..)) : (n : ℝ≥0∞) ≤ μ C

theorem not_packs_of_measure (m : ℝ) (μ : Measure (ℝ × ℝ))
    (hcover : ∀ c θ, sq c θ 1 ⊆ box m → 1 ≤ μ (sq c θ 1))
    (n : ℕ) (htot : μ (box m) < n) {s : ℝ} (hs : s < m) : ¬ Packs n s
```

Any measure at all (values in `ℝ≥0∞`, no finiteness needed).  Proof as FORMAT.md says: the concentric closed unit
squares of the `L`-squares are closed (`isClosed_sq`), hence measurable, pairwise disjoint (inside the disjoint
interiors), and inside `C`, so `n ≤ Σ μ(Q_i) = μ(⋃ Q_i) ≤ μ(C)`.  `not_packs_of_measure` is `not_packs_of_cover`'s
scaling argument on top of it.

## 2.  D4 for measures (`MixedMeasure.lean` §3)

```lean
def D4InvM (m : ℝ) (μ : Measure (ℝ × ℝ)) : Prop :=
  ∀ s, MeasurableSet s → μ (reflX m ⁻¹' s) = μ s ∧ μ (swapXY ⁻¹' s) = μ s

theorem d4_reduction_measure (m) (μ) (hinv : D4InvM m μ)
    (hreg : ∀ c θ, c.1 ∈ Icc 0 (m/2) → c.2 ∈ Icc 0 (m/2) → θ ∈ Icc 0 (π/4) → sq c θ 1 ⊆ box m →
      1 ≤ μ (sq c θ 1)) :
    ∀ c θ, sq c θ 1 ⊆ box m → 1 ≤ μ (sq c θ 1)
theorem d4_reduction_measure_u ...   -- hreg over θ = 2 arctan u, u ∈ [0, 1/2]
```

Same proof as `d4_reduction` (`D4.lean`): the pose bookkeeping `d4_reduce` is reused verbatim; the only new input
is `μ(sq (g c) (−θ) 1) = μ(sq c θ 1)` for the two generators, from `g⁻¹'(sq (g c) (−θ) 1) = sq c θ 1`
(`mem_sq_reflX`, `mem_sq_swapXY`) and invariance.

## 3.  Mixed covers as measures (`MixedMeasure.lean` §4–5)

`MixedCover ι κ ν`: finite families of weighted points `pt i, pw i`, segments `[sa j, sb j]` with mass `sw j`, and
regions `poly k` (any sets; convex polygons in the files) with mass `gw k`.  Its measure is

```
M.measure = Σ_i ofReal(pw i) • δ_{pt i}
          + Σ_j ofReal(sw j) • segMeasure (sa j) (sb j)     -- push-forward of Lebesgue on [0,1] by t ↦ a + t(b − a)
          + Σ_k ofReal(gw k) • areaMeasure (poly k)         -- (volume P)⁻¹ • volume.restrict P
```

The push-forward of Lebesgue on `[0,1]` is exactly "mass spread uniformly by length" for a non-degenerate segment
(and it needs no Hausdorff measure).  Proved:

| lemma | statement |
|---|---|
| `MixedCover.measure_apply` | for measurable `S` (e.g. a closed square), `M.measure S = ofReal (M.mass S)`, with `M.mass S = Σ_{pt i ∈ S} pw i + Σ_j sw j · segFrac (sa j) (sb j) S + Σ_k gw k · areaFrac (poly k) S` |
| `segFrac a b S` | `(volume ((t ↦ a + t(b−a))⁻¹' S ∩ [0,1])).toReal`: the parametric fraction of FORMAT.md, `= |[a,b] ∩ S|/|[a,b]|` for convex `S` |
| `areaFrac P S` | `area(S ∩ P)/area(P)` (`areaFrac_eq`; `0` for a zero- or infinite-area region) |
| `MixedCover.measure_univ_le` | `M.measure univ ≤ ofReal M.total` (the file's total), always |
| `MixedCover.measure_univ_eq` | equality when every polygon has `0 < area < ∞` |
| `segMeasure_comm`, `segFrac_comm` | orientation of a segment does not matter |
| `MixedCover.measure_preimage_invol`, `MixedCover.d4InvM` | `D4InvM m M.measure` from entry-level invariance: involutions of the index sets mapping points to points (by `x ↦ m−x`, resp. `x ↔ y`, same mass), segments to segments (endpoints mapped, **either orientation**, same mass), polygons to polygons (`poly (g k) = g⁻¹' poly k`, same mass; uses that both generators preserve Lebesgue measure, `measurePreserving_reflX/swapXY`) |

Non-negativity of the masses is the only side condition (`MixedCover.Nonneg`).  The containment "all pieces in
`[0,s]²`" is not needed: the bound uses `μ(box m) ≤ μ(univ) ≤ total`.

## 4.  The `s(21)` cover and the hypothesis (`S21.lean`)

**Data** (`S21Data.lean`, 339 KB, generated by `lean/scripts/gen_s21_data.py` from the certificate, sha256
`8b415ceeb5f20b02c4e338ce2feb21bc206051ef6bc5c7e433bebae2ef39fc23` in the header): points as a `PTree` of
`(X, Y, w)` (`Cover.lean`), segments as an `STree` of `(X0, Y0, X1, Y1, w)` (`SegTree.lean`, same design: search
tree on the key, `mem_sound` holds for any tree, so ordering is untrusted).  The generator puts each segment's
endpoints in lexicographic order (harmless by `segMeasure_comm`), checks the header (`mixed 1`, `s = 5/1`,
`D = 1000`, `W = 10¹¹`, `npg = 0`) and refuses repeated points/segments and zero-length segments.
`python3 lean/scripts/gen_s21_data.py --check` (from `s12/`) regenerates and compares byte for byte.

```lean
def pt (e) := (X/1000, Y/1000)        def pw (e) := w/10¹¹           -- point entries
def sa (e) := (X0/1000, Y0/1000)      def sb (e) := (X1/1000, Y1/1000)      def sw (e) := w/10¹¹
noncomputable def cover : MixedCover (ℕ×ℕ×ℕ) SegE Empty     -- no polygons
noncomputable def μ := cover.measure
```

**Proved about the data** (`decide +kernel`, ~25 s for the file in all):

```lean
theorem S21Data.check_ok : check = true   -- both trees strictly increasing (no repeats); X ≤ 5000; segments stored
                                          -- normalised (so non-degenerate), X0, X1 ≤ 5000; the images of every
                                          -- entry under x ↦ 5 − x and x ↔ y (segments re-normalised) are entries
theorem S21Data.pwsum_tree : ptree.wsum = 350685505788
theorem S21Data.swsum_tree : stree.wsum = 1738789413944
theorem S21Data.card_pentries : pentries.card = 7536
theorem S21Data.card_sentries : sentries.card = 1872
theorem S21Data.total_eq : cover.total = 2089474919732 / 100000000000
theorem S21Data.mu_box_lt : μ (box 5) < 21
theorem S21Data.d4 : D4InvM 5 μ
theorem S21Data.mu_sq (c θ) : μ (sq c θ 1) = ofReal (Σ_{pt e ∈ sq c θ 1} pw e + Σ_e sw e * segFrac (sa e) (sb e) (sq c θ 1))
```

**The hypothesis** — the only input not proved in Lean:

```lean
def S21RegionCover : Prop :=
  ∀ (c : ℝ × ℝ) (θ : ℝ), c.1 ∈ Set.Icc 0 (5 / 2) → c.2 ∈ Set.Icc 0 (5 / 2) →
    θ ∈ Set.Icc 0 (Real.pi / 4) → sq c θ 1 ⊆ box 5 →
      1 ≤ ∑ e ∈ pentries.filter (fun e => pt e ∈ sq c θ 1), pw e
        + ∑ e ∈ sentries, sw e * segFrac (sa e) (sb e) (sq c θ 1)

def S21CheckerCover : Prop :=      -- the same over c ∈ [0,5/2]², θ = 2 arctan u, u ∈ [0, 1/2]
theorem S21CheckerCover.region (h : S21CheckerCover) : S21RegionCover
```

In words: every closed unit square inside `[0,5]²` with centre in `[0,2.5]²` and angle in `[0°,45°]` (resp.
`θ = 2 arctan u`, `u ∈ [0,½]`) gets mass `≥ 1`: the points in it (closed square) plus, for each segment, its mass
times the fraction of the segment inside it (a segment on the square's edge counts in full).  The hypothesis is
real-valued and measure-free apart from `segFrac`'s one-dimensional Lebesgue measure.

## 5.  Statement ↔ computation

| `S21CheckerCover` | `zm_mixed.py --d4 --cert-mode` run `certificates/s21/zm_mixed_d4/` (header: sha256 of checker, reader, `zeromargin.py`, cover; settings), and `zmx2 --d4` (`certificates/s21/zmx2_d4/`) |
|---|---|
| `c.1, c.2 ∈ [0, 5/2]`, `u ∈ [0, ½]` | the 40,000 closed root boxes (pitch `1/20`, 16 `u`-bins), checked against `d4_roots` |
| `θ = 2 arctan u` | the checker's parameter `u = tan(θ/2)` |
| `sq c θ 1 ⊆ box 5` | admissible (`Adm 5 c θ`, `sq_subset_box_iff`); non-admissible poses exempt in both |
| `pt e ∈ sq c θ 1` | `|X|, |Y| ≤ ½` in the rotated frame (`coord`) — closed square |
| `segFrac (sa e) (sb e) Q` | FORMAT.md's `|seg ∩ Q| / |seg|` (parametric fraction), closed `Q` |
| `pw`, `sw`, `pt`, `sa`, `sb` | the file's integers over `W = 10¹¹`, `D = 1000` |
| sums over `pentries`, `sentries` | sums over file lines: no repeats (`check_ok`), counts 7,536 / 1,872 match the header |
| `1 ≤ …` | every leaf EMPTY, or `PIECE`/`ADM`/`CHAIN`/`SPLIT` with a certified lower bound `≥ 1` at every admissible pose |

So the run's verdict **is** `S21CheckerCover`, provided `zm_mixed.py` is correct (its lemmas are proved on paper in
`ZM_MIXED.md` §2 and tested, not formalised; audit `ZM_MIXED_AUDIT.md`).  `zmx2`'s `--d4` verdict is the same
statement from an independent implementation (with the float caveat of `certificates/s21/README.md`), and its
`--full` run establishes the unreduced statement directly.

**What Lean now does that the checker's write-up argued:** the reduction from mixed measures to packings (FORMAT.md's
"`packing_le_weight` needs generalising" — done, for arbitrary measures), the symmetry argument for measures
(`ZM_MIXED.md` §2 "Symmetry" — `d4_reduction_measure` + `S21Data.d4`, which checks invariance of the actual entries
including segment orientation), and the total `< 21`.

**Not Lean-verified, plainly:** (1) `S21CheckerCover` (the checker run); (2) that `S21Data.lean` is the file the
checker ran on (generator + `--check`, sha256 in the header).  Everything else in `s21_eq_five_of_checker` is proved.

## 6.  Files, build, axioms

* `lean/Sqpack/MixedMeasure.lean` — 494 lines: measures → packings, D4 for measures, `MixedCover`.
* `lean/Sqpack/SegTree.lean` — 144 lines: `STree`, `segNorm`, soundness lemmas.
* `lean/Sqpack/S21Data.lean` — generated (339 KB): the cover.
* `lean/Sqpack/S21.lean` — 319 lines: data facts, hypothesis, theorem.
* `lean/scripts/gen_s21_data.py` — the generator.
* `lean/Sqpack.lean` imports the four modules; `lean/Axioms.lean` gains 28 `#print axioms` lines.  All show
  `[propext, Classical.choice, Quot.sound]` or a subset (`STree.mem_sound`, `S21Data.check_ok`: `[propext]`;
  `pwsum_tree`, `swsum_tree`: none).  No `Lean.ofReduceBool`, no `sorry`.  `S32.lean` and everything before it are
  untouched; `s32_eq_six_of_checker` still builds with the same axioms.

Build (2 cores, `taskset -c 12-13`, Mathlib from cache): `MixedMeasure` ~4 s, `S21Data` + `S21` ~30 s (kernel
evaluation of `check_ok` dominates).
