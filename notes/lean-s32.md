# `s(32) = 6` in Lean, from one computational hypothesis (`lean/Sqpack/{D4,Cover,S32Data,S32}.lean`)

**What this is.**  The computer-assisted proof of `s(32) = 6` (`search/S32_EXACT.md` §10, §11) is a weighted
point cover of `[0,6]²` of total weight `31.7135… < 32` such that every closed unit square in `[0,6]²` captures
weight `≥ 1`; the exact checkers verify that only on the D4 fundamental region.  The Lean files state and prove
the theorem **from a single named hypothesis about that region**, with everything else proved: the D4 reduction,
the D4 invariance and the total weight of the actual cover (as data), the scaling argument, and the packing in
side 6.  Mathlib `v4.33.1`, `lake build` clean, no `sorry`, standard axioms only, no `native_decide`.

```lean
theorem s32_eq_six (h : S32RegionCover) : minSide 32 = 6
theorem s32_eq_six_of_checker (h : S32CheckerCover) : minSide 32 = 6
```

---

## 1.  The statements

**Packings and `s(n)`** (`S32.lean` §1–2):

```lean
def Packs (n : ℕ) (s : ℝ) : Prop :=
  ∃ (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ), (∀ i, sq (ctr i) (ang i) 1 ⊆ box s) ∧
    ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) 1) (sqInt (ctr j) (ang j) 1)

noncomputable def minSide (n : ℕ) : ℝ := sInf {s | Packs n s}

theorem s32_packs : Packs 32 6                                            -- no hypothesis
theorem s32_not_packs (h : S32RegionCover) {s : ℝ} (hs : s < 6) : ¬ Packs 32 s
theorem s32_isLeast (h : S32RegionCover) : IsLeast {s | Packs 32 s} 6
```

`Packs n s`: `n` closed unit squares, any centres and angles, inside the closed square `[0,s]²`, pairwise disjoint
interiors — the usual definition of `s(n)`, stated with `Basic.lean`'s `sq`/`sqInt` and `Chord.lean`'s `box`.
`s(n)` is defined as the infimum; `s32_isLeast` says more (the minimum is attained, at 6), and `s32_eq_six` is
`IsLeast.csInf_eq`.  The upper bound is `packs_grid` (the `n × n` axis-parallel grid holds any `k ≤ n²`
squares), used at `n = 6, k = 32`.  The lower bound is `not_packs_of_cover` (the classical scaling: a packing in
`box s`, `s < 6`, scaled by `6/s > 1` is a packing of larger squares in `box 6`, and `packing_le_weight` bounds
it by the total weight).

**The cover** (`S32.lean` §3), transcribed from `runs/s32-close_candidate.txt` into `S32Data.tree`:

```lean
def entries : Finset (ℕ × ℕ × ℕ) := tree.toList.toFinset                      -- the (x, y, w) lines
noncomputable def pt (e : ℕ × ℕ × ℕ) : ℝ × ℝ := ((e.1 : ℝ) / 1000, (e.2.1 : ℝ) / 1000)
noncomputable def wt (e : ℕ × ℕ × ℕ) : ℝ := (e.2.2 : ℝ) / 100000000000
```

**The hypothesis** (`S32.lean` §4) — the only input that is not proved in Lean:

```lean
def S32RegionCover : Prop :=
  ∀ (c : ℝ × ℝ) (θ : ℝ), c.1 ∈ Set.Icc 0 3 → c.2 ∈ Set.Icc 0 3 →
    θ ∈ Set.Icc 0 (Real.pi / 4) → sq c θ 1 ⊆ box 6 →
      1 ≤ ∑ e ∈ entries.filter (fun e => pt e ∈ sq c θ 1), wt e

def S32CheckerCover : Prop :=
  ∀ (c : ℝ × ℝ) (u : ℝ), c.1 ∈ Set.Icc 0 3 → c.2 ∈ Set.Icc 0 3 →
    u ∈ Set.Icc 0 (1 / 2) → sq c (2 * Real.arctan u) 1 ⊆ box 6 →
      1 ≤ ∑ e ∈ entries.filter (fun e => pt e ∈ sq c (2 * Real.arctan u) 1), wt e

theorem S32CheckerCover.region (h : S32CheckerCover) : S32RegionCover
```

In words: every closed unit square inside `[0,6]²` whose centre is in `[0,3]²` and whose angle is in `[0°,45°]`
(resp. `θ = 2 arctan u`, `u ∈ [0,½]`, i.e. `[0°, 53.13°]`) contains certificate points of total weight `≥ 1`.
`S32CheckerCover` is the stronger statement and the one the checker run literally establishes (§3).

**Proved about the data** (kernel evaluation by `decide +kernel`; ~30 s for `check_ok`, under 1 s for the rest):

```lean
theorem S32Data.check_ok : check = true          -- keys strictly increasing, x ≤ 6000, both reflections present
theorem S32Data.wsum_tree : tree.wsum = 3171350535386
theorem S32Data.card_entries : entries.card = 13085
theorem S32Data.total_lt : ∑ a ∈ A, w a < 32     -- 3171350535386 / 10¹¹ = 31.713505…
theorem S32Data.d4Inv : D4Inv 6 A w              -- invariant under x ↦ 6 − x and x ↔ y
```

`(A, w)` is the cover as `packing_le_weight` takes it (`coverA`, `coverW`: the points, and the total weight at each
point); `sum_filter_coverA` shows that the weight `(A, w)` captures in a square is the hypothesis' sum over the
entries in it.

---

## 2.  The D4 reduction (`D4.lean`)

```lean
def D4Inv (m : ℝ) (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) : Prop :=
  ∀ a ∈ A, (reflX m a ∈ A ∧ w (reflX m a) = w a) ∧ (swapXY a ∈ A ∧ w (swapXY a) = w a)

theorem d4_reduction (m : ℝ) (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hinv : D4Inv m A w)
    (hreg : ∀ (c : ℝ × ℝ) (θ : ℝ), c.1 ∈ Set.Icc 0 (m / 2) → c.2 ∈ Set.Icc 0 (m / 2) →
      θ ∈ Set.Icc 0 (Real.pi / 4) → sq c θ 1 ⊆ box m →
        1 ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a) :
    ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ box m →
      1 ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a

theorem d4_reduction_u ...   -- the same with hreg over θ = 2 arctan u, u ∈ [0, 1/2]
```

The conclusion is verbatim `packing_le_weight`'s `hcover` for `C = box m`.  The proof is §11.3(d) of
`S32_EXACT.md`, in four pieces:

| step | Lean |
|---|---|
| the axis square is invariant under the quarter turn: `sq c (θ + π/2) 1 = sq c θ 1` | `sq_add_pi_div_two` |
| `x ↦ m − x` and `x ↔ y` map `sq c θ 1` onto `sq (g c) (−θ) 1` | `mem_sq_reflX`, `mem_sq_swapXY` |
| they preserve `sq ⊆ box m` (via `Adm`, `sq_subset_box_iff`) and, for a `D4Inv` weight, the captured weight | `sq_subset_box_{reflX,swapXY}_iff`, `capt_reflX`, `capt_swapXY` |
| every pose `(c, θ)` with `c ∈ box m` is carried into `[0,m/2]² × [0, π/4]` | `d4_reduce` |

`d4_reduce` is abstract (any pose predicate invariant under the three moves): reduce `θ` mod `π/2` into
`[0, π/2)` (`toIcoMod`); if `θ > π/4`, apply `x ↦ m − x` (`θ ↦ −θ ≡ π/2 − θ`); then the quarter turn
`(x, y) ↦ (m − y, x)` = `reflX ∘ swapXY`, which keeps `θ`, takes the centre's closed quadrant to `[0,m/2]²`.
`D4Inv` mentions only the two generators, exactly what `Checker.symmetric_d4` and `zmcheck --d4` test.

`exists_u_of_theta`: for `θ ∈ [0, π/4]`, `u = tan(θ/2)` has `2 arctan u = θ` and `0 ≤ u ≤ ½` (from
`cos θ ≥ sin θ`, i.e. `1 − u² ≥ 2u`); so the checkers' domain `u ∈ [0, ½]` contains the fundamental region and
`d4_reduction_u`, `S32CheckerCover.region` follow.

---

## 3.  Statement ↔ computation

**What establishes `S32CheckerCover`.**  The `zeromargin.py` D4 run, `runs/s32py_d4/SUMMARY.txt`
(`S32_EXACT.md` §11): checker `runs/s32py_d4/checker/zeromargin.py`, sha256
`640fe453c1a32f4aa580ca2b1261c6406923a4d7c131f65604a432c7fc2086ab`, on `runs/s32-close_candidate.txt`, sha256
`a0d2d38fc9a585a166b9e06c5069fdae9bdce44ca9fb64dda75abcc99e3c2144`; all **7,200 roots** of
`[0,3]² × u ∈ [0,½]` (centre cells of pitch `1/10`, 8 `u`-bins of width `1/16`) present and certified,
`UNCERTIFIED 0` (7,199 at depth 24, one root, `x ∈ [1.5,1.6], y ∈ [0.5,0.6], u ∈ [0,1/16]`, by the
`--deepen 30` record).  Term by term:

| `S32CheckerCover` | `zeromargin.py` |
|---|---|
| `c.1, c.2 ∈ [0,3]`, `u ∈ [0,½]` | the union of the 7,200 *closed* root boxes, which is exactly `[0,3]² × [0,½]` (`d4_roots`; §11.3(d)(5)) |
| `θ = 2 arctan u` | the checker's parameter `u = tan(θ/2)` (`notes/lean-zeromargin.md`; `cos_two_arctan`, `sin_two_arctan`) |
| `sq c θ 1 ⊆ box 6` | *admissible*: `Adm 6 c θ`, equivalent by `sq_subset_box_iff`; poses of a root box that are not admissible are exempt in both |
| `pt e ∈ sq c θ 1` | `|X|, |Y| ≤ ½` with `(X, Y) = coord c θ p` (literally `Basic.lean`'s `coord`) |
| `pt e = (x/1000, y/1000)`, `wt e = w/10¹¹` | the certificate's `D = 1000`, `W = 100000000000`; the checker reads the same integers |
| sum over `entries` (duplicates would add) | the checker's weight of the captured points; the file has no repeated point (`check_ok`'s strict key order) |
| `1 ≤ …` | a root is certified only when each leaf box is EMPTY (no admissible pose) or has points of weight `≥ 1` proved captured at every admissible pose of the box (ADM / CHAIN leaves; exact arithmetic) |

So the run's verdict is the statement `S32CheckerCover`, provided the checker is correct; the soundness of its
exact primitives is what `lean/Sqpack/ZeroMargin.lean` formalises (`notes/lean-zeromargin.md`), but the program
(subdivision, bookkeeping, parsing) is not verified.  Independently, `zmcheck --d4` (§10) certifies the same
region (3,600 roots of `u`-width `1/8`; its 5 germ roots by `zeromargin.py`), a second witness for the same
hypothesis.

**What Lean now does that the checkers did.**  The D4 argument (§7 / §11.3(d)) and the invariance test
(`symmetric_d4`, `--d4`) are replaced by `d4_reduction` and `S32Data.d4Inv`, and the total `31.7135… < 32` by
`S32Data.total_lt`; the checker run is needed only for its 7,200 per-root verdicts.

**The transcription.**  `lean/Sqpack/S32Data.lean` (439 KB, committed because the build needs it) is generated
by `lean/scripts/gen_s32_data.py` from the certificate: one tree node `(x, y, w)` per certificate line, sorted by
`(x, y)`, the certificate's sha256 in the file header; the script checks the header (`m = 6`, `D = 1000`,
`W = 10¹¹`, `n = 13085`) and refuses repeated points.  `python3 lean/scripts/gen_s32_data.py
runs/s32-close_candidate.txt --check` (from `s12/`) regenerates it in memory and compares byte for byte.  In Lean,
`card_entries = 13085` and `wsum_tree = 3171350535386` match the file's count and total.  The search-tree
order is not trusted: `PTree.mem_sound` holds for any tree, so a mis-ordered tree can only fail `check_ok`.

**Not Lean-verified, plainly:** (1) `S32CheckerCover` (the checker run); (2) that `S32Data.lean` is the
certificate the checker ran on (the generator plus its `--check`; the sha256 in the header).  Everything else in
`s32_eq_six_of_checker` is proved.

---

## 4.  Files, build, axioms

* `lean/Sqpack/D4.lean` — 298 lines: the D4 reduction.
* `lean/Sqpack/Cover.lean` — 222 lines: integer data → `(A, w)`; `PTree` and its checks with soundness lemmas.
* `lean/Sqpack/S32Data.lean` — 214 lines (generated, 439 KB): the cover.
* `lean/Sqpack/S32.lean` — 319 lines: `Packs`, `minSide`, scaling, grid, the data facts, the theorem.
* `lean/scripts/gen_s32_data.py` — the generator.
* `lean/Axioms.lean` — now 63 `#print axioms`.  All show `[propext, Classical.choice, Quot.sound]` or a subset
  (`PTree.mem_sound`, `check_ok`: `[propext]`; `wsum_tree`: none).  No `Lean.ofReduceBool`: the data facts
  use `decide +kernel`, which the kernel checks by evaluation, not `native_decide`.

Build (8 cores, Mathlib from cache, own modules only): `S32Data` 10 s, `S32` 36–39 s (of which `check_ok` ~30 s
of kernel type checking), the rest ~25 s; whole `lake build` from a clean `.lake/build` 76 s.
