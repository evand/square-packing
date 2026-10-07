# Lean: `Valid7` / `Valid9` reduced to the tilted run's region; D4 and Lemma Z in Lean (2026-10-03)

Task `lean-valid-split`.  Files: `lean/Sqpack/ValidSplit.lean` (generic), `lean/Sqpack/ValidSplit7.lean`,
`lean/Sqpack/ValidSplit9.lean` (the two instances), `lean/Sqpack/ValidSplitData.lean` (generated: the box files'
weights packed into one natural number each), `lean/scripts/gen_validsplit_data.py` (the generator, with a Python
mirror of the checks).  All in the default build.  Predecessors: `notes/lean-bentz-reduction.md` (`Valid7`),
`notes/lean-k2m4-reduction.md` (`Valid9`).

`certificates/k2m3/README.md` "What is not machine-verified" listed (i) `Valid7` itself, (ii) the D4 reduction to the
fundamental domain, (iii) Lemma Z's reduction to finitely many corner limits.  (ii) and (iii) are now kernel-checked,
for `Valid7` and for `Valid9`; what remains of (i) is exactly the pose region of the `qx2_zm.py` box run.

## The theorems

```lean
-- generic (ValidSplit.lean)
def ValidTilt (m : ℝ) (μ : Measure (ℝ × ℝ)) : Prop :=
  ∀ (c : ℝ × ℝ) (u : ℝ), c.1 ∈ Set.Icc 0 (m / 2) → c.2 ∈ Set.Icc 0 (m / 2) → 0 < u →
    u ^ 2 + 2 * u ≤ 1 → sq c (2 * Real.arctan u) 1 ⊆ box m →
      1 ≤ μ (sq c (2 * Real.arctan u) 1)
def ValidAxis (m : ℝ) (μ : Measure (ℝ × ℝ)) : Prop :=
  ∀ c : ℝ × ℝ, sq c 0 1 ⊆ box m → 1 ≤ μ (sq c 0 1)
theorem valid_of_tilt_axis {m : ℝ} {μ : Measure (ℝ × ℝ)} (hinv : D4InvM m μ)
    (ht : ValidTilt m μ) (ha : ValidAxis m μ) :
    ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ box m → 1 ≤ μ (sq c θ 1)

-- k² − 3 (ValidSplit7.lean, namespace SquarePacking.Bentz)
def ValidTilt7 : Prop := ValidTilt 7 box7Cover.measure
def ValidAxis7 : Prop := ValidAxis 7 box7Cover.measure
theorem d4_box7 : D4InvM 7 box7Cover.measure
theorem valid7_of_tilt_axis (ht : ValidTilt7) (ha : ValidAxis7) : Valid7
theorem validAxis7 : ValidAxis7                                    -- Lemma Z, proved
theorem valid7_of_tilt (ht : ValidTilt7) : Valid7
theorem bentz_of_validTilt7 (ht : ValidTilt7) : ∀ k : ℕ, 6 ≤ k → minSide (k ^ 2 - 3) = k

-- k² − 4 (ValidSplit9.lean, namespace SquarePacking.Bentz4): the same with 9, box9Cover, Valid9,
theorem bentz4_of_validTilt9 (ht : ValidTilt9) : ∀ k : ℕ, 8 ≤ k → minSide (k ^ 2 - 4) = k
```

`box7Cover`, `box9Cover`, `Valid7`, `Valid9` are the unchanged definitions of `Bentz.lean` / `Bentz4.lean` (the box
files verbatim).  `#print axioms` → `[propext, Classical.choice, Quot.sound]` for every theorem above (all listed in
`lean/Axioms.lean`).  No `sorry`, no `native_decide`, no new axioms.

So `s(k² − 3) = k` (`k ≥ 6`) and `s(k² − 4) = k` (`k ≥ 8`) now rest, beyond Lean's kernel, on `ValidTilt7` /
`ValidTilt9` alone: the claims of the two `qx2_zm.py` box runs.

## `ValidTilt` is what the run certifies (checked against `qx2_zm.py`; no mismatch found)

`qx2_zm.py` (sha `6294052a…`, the file of run V3 and of `qx2_k4x_k008/qxzm_full`), read for this:

* **Roots.**  `zeromargin.d4_roots(m, 1/10, 8)`: centre cells of pitch `1/10` covering `[0, m/2]²` (closed), and
  `u = tan(θ/2)` in 8 bins covering `[0, 1/2]`.  Both runs of record are full (argv has no `--cx-lo` etc.; V3:
  9,800 roots; `qxzm_full`: 16,200 roots, header argv checked).
* **Leaves that certify nothing by themselves.**  `AXIS` (`u₁ = 0` after `clip_bin`: only `θ = 0` poses left) is
  delegated to Lemma Z; `SYM` (`u₀² + 2u₀ − 1 ≥ 0`, i.e. `θ ≥ 45°` on the whole box) is delegated to the diagonal
  symmetry; `EMPTY` and the slabs cut by `clip_bin` contain no admissible pose (admissible = `w(θ)/2 ≤ c ≤ m −
  w(θ)/2`, which is `sq c θ 1 ⊆ box m`: `sq_subset_box_iff`).  `EXACT45` certifies only the `θ ≤ 45°` part of its
  box.
* **So the run's own claim** is: every admissible pose with centre in `[0, m/2]²` and `0 < θ ≤ 45°` (`u > 0`,
  `u² + 2u ≤ 1`) has mass `≥ 1`.  That is `ValidTilt`.  The `θ = 45°` poses (`u = √2 − 1`, irrational) are never in
  a `SYM` leaf (a `SYM` leaf has rational `u₀ ≥ √2 − 1`, so `u₀ > √2 − 1`), so they are covered by certifying leaves.
  The roots' extra range `45° < θ ≤ 53.13°` is not needed and not assumed.
* **What the reduction needs** (`d4_reduction_measure`, `D4.lean`): centres in `[0, m/2]²`, `θ ∈ [0, π/4]`.  `θ = 0`
  is `ValidAxis`; `0 < θ ≤ π/4` is `ValidTilt` with `u = tan(θ/2)` (`exists_u_of_theta_pos`).
* **Angle convention.**  `zeromargin.in_rot_square` (`x = a c + b s`, `y = −a s + b c`) is Lean's `coord`; and the
  domain is in any case invariant under `θ ↦ −θ` up to the swap `x ↔ y` (which keeps `[0, m/2]²`), so a sign slip in
  a checker primitive could not move the run off the fundamental domain.

**Minor (documentation):** `search/K2M4_MARGIN.md` reports Lemma Z for `K4_k008_box9.txt` as "min 1, 348 tight".
`qx2_zm.py axis` on that file gives min 1 with **208** tight corners (1,600 corners on `[1/2, 9/2]²`); 348 is the LP
projection's count of tight `θ = 0` rows (`runs/qx2_k4x_k008/project.out`: 35,646 equalities − 35,297 germ rows − σ).
The Lean check agrees with `qx2_zm.py axis`: over the whole box 832 = 4·208 tight corners (and 752 = 4·188 for
`L4_k02_box7.txt`, whose `lemmaZ.out` says 188).

## How it is proved

1. **Grid covers** (`ValidSplit.gridCover K A den w`): mass `w s / den` on each unit segment `s` of the `1/5`-grid in
   `[0,K]²`, uniform by length, plus Lebesgue on `[A/5, K − A/5]²`.  Both box files are grid covers:
   `box7Cover_measure` (`Bentz.lean`) and `box9Cover_measure` (`Bentz4.lean`) already identify them with the
   families' `μ₇`, `μ₉`; `measure_eq_gridCover` / `famCover_measure_eq` (two layers added) turn those into grid
   covers.  For speed the weights are then replaced by a packed table (step 4), tied to the family by the kernel
   (`fitOK`: equal on the whole grid).
2. **D4** (`d4InvM_gridCover`): with the index maps `reflIx` (`x ↦ K − x`: horizontal `i ↦ 5K − 1 − i`, endpoints
   swapped; vertical `i ↦ 5K − i`) and `swapIx` (horizontal `(i, j)` ↔ vertical `(j, i)`), the existing
   `MixedCover.d4InvM` (`MixedMeasure.lean`) gives `D4InvM K` from invariance of the weights (`SymW`, decided by
   `symOKP`), the Lebesgue square being symmetric.
3. **Lemma Z** (`validAxis_gridCover`).  The axis square `[x₀, x₀+1] × [y₀, y₀+1]` (`x₀, y₀ ∈ [0, K−1]`); write `5x₀
   = p + σ`, `5y₀ = q + τ`, `p, q ∈ {0, …, 5K−6}`, `σ, τ ∈ [0,1]` (`exists_cell`; at the far wall `σ = 1`).  At the
   corner `(a, b) ∈ {0,1}²` of the cell the "inside" limit counts in full the horizontal unit segments `p+a ≤ i ≤
   p+a+4`, `q+1 ≤ j ≤ q+5` and the vertical ones `p+1 ≤ i ≤ p+5`, `q+b ≤ j ≤ q+b+4` (`cornerSet`, `cornerSum`: 50
   segments), and the Lebesgue overlaps `lov(p+a)·lov(q+b)/25` (`lov K A n = min(5K − A, n+5) − max(A, n)`, in
   fifths).  Soundness:
   * per segment (`seg_lb`), its length fraction in the closed square is at least the bilinear interpolation
     `Σ λ_ab [s ∈ cornerSet_ab]`, `λ_ab = (1−σ or σ)(1−τ or τ)`: the end segments `i = p` / `i = p+5` have
     fractions `≥ 1 − σ` / `≥ σ` (a sub-interval of the parameter range lies in the square, `segFrac_ge`), the
     middle ones `1`; a line `j = q` / `j = q+5` *on* the square's boundary is in the closed square, so the
     indicator is at least the cell's interior value — this is where closedness enters;
   * the Lebesgue overlap of `[A, 5K − A]` with `[n, n+5]` is at least linear on each cell (`ovl_interp`: four
     cases by integer comparisons, each `max(c + dσ, 0)` with integer `c`, `d ∈ {−1, 0, 1}`, `interp_max`), so the
     area, a product of two overlaps (`vol_inter`), is at least the bilinear interpolation;
   * so mass `≥ Σ λ_ab (cornerSum_ab/den + lov·lov/25) ≥ Σ λ_ab = 1` given `AxisW`: every corner value `≥ 1`.
   This is Lemma Z of QUADRANT_EXACT.md §4.1 in "cell + four closed corners" form; it covers the whole box (3,600
   corner limits for `K = 7`, 6,400 for `K = 9`), no symmetry used.
4. **Kernel checks** (`ValidSplit7.lean`, `ValidSplit9.lean`): `outOK` (family weights vanish outside `[0,K]²`),
   `fitOK` (packed table = family weights on the grid), `symOKP`, `axisOKP` (all four corners of all cells), by
   `decide +kernel`.

## Engineering: why packed weights

A first version stated `axisOK` with `Finset` sums over `range 5 ×ˢ range 5` of the family weight function: the
kernel took ~3 s per corner (the `Finset`/`Multiset` machinery), i.e. hours.  With plain `f 0 + … + f 4` sums it is
~60 ms per corner, dominated by the weight function itself: an evaluation of `Bentz.wN` / `BentzFam.wN` (`ℤ`
comparisons and `%` in the codes, a list lookup) costs ~1.2 ms in the kernel (`ℤ` arithmetic goes through the
`Int` constructors; only `ℕ` has GMP acceleration).  So the weights are packed into one natural-number literal per
box file (48 bits per grid segment, `ValidSplitData.lean`, ~95 KB), read with `>>>` and `%` (GMP, ~0.1 ms with all
`ℕ` index arithmetic), and tied to the family by one pass of `fitOK` (one family-weight evaluation per grid segment).
`outOK` is the only other check that evaluates the family weights.

## Timings

Kernel checks (`decide +kernel`, one `lake env lean` per check, cores 8–15):

| check | `K = 7` | `K = 9` |
|---|---|---|
| `outOK` | 0.15 s | 0.09 s |
| `fitOK` (2·36², 2·46² grid segments) | 4.3 s | 8.0 s |
| `symOKP` | 0.24 s | 0.46 s |
| `axisOKP` (3,600 / 6,400 corners) | 4.4 s | 9.0 s |

Build (`lake build`, `LEAN_NUM_THREADS=8`, cores 8–15): `ValidSplit` 10 s, `ValidSplitData` 2 s, `ValidSplit7` 12 s,
`ValidSplit9` 20 s.  `python3 lean/scripts/gen_validsplit_data.py` (from `s12/`, deterministic, ~5 s) regenerates
the data and runs its mirror (symmetry; Lemma Z: `K = 7`: 3,600 corners, min exactly 1, 752 tight; `K = 9`: 6,400
corners, min exactly 1, 832 tight).  Mutation sanity (`#eval`, not committed): lowering one positive weight of a
tight corner by one unit (`10⁻¹²` resp. `(2·10¹²)⁻¹`) makes `axisOKP` false for both files; lowering one weight makes
`symOKP` false.

## What is left for `Valid7` / `Valid9` in Lean

`ValidTilt7`, `ValidTilt9`: the `θ > 0` part, which the Python runs certify with the zero-margin leaf primitives
(PIECE = Lemmas S/T/L/R with the polygon Lemma S(b); LEB = Lemma U; CAP = Lemma K; EXACT = Lemma E with E′, E″).
The existing kernel verifiers (`ZMTree`, `ZMTreeM`, `CovM`) cover points and segments but not the Lebesgue polygon
and not Lemma E; that is the next gap (QUADRANT_EXACT.md §3–4, ZM_MIXED.md §2).
