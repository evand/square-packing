# The Lean lower-bound ladder

Goal: kernel-checked lower bounds `s(n) ≥ t` from point certificates — **kernel reduction only**
(`decide +kernel`; no `native_decide`, no new axioms), with **one generic verifier proved sound
once**, and per result only data plus a one-line theorem.

Status (2026-09-27): rung 1 (pilot) done — `SquarePacking.s12_ge_35_9 : (35/9 : ℝ) ≤ minSide 12`,
`#print axioms` = `[propext, Classical.choice, Quot.sound]`, no `sorry`.  The same verifier, unchanged,
also proves rung-3 results from our weighted certificates: `s12_ge_3920_997 : (3920/997 : ℝ) ≤
minSide 12` (`S12WLower.lean`, 224 points, in the default build) and `s11_ge_3040_797 :
(3040/797 : ℝ) ≤ minSide 11` (`S11Lower.lean`, 680 points; **opt-in**, not imported by
`Sqpack.lean`: `lake build Sqpack.S11Lower`, ~26 min on 3 cores).  All print the same three
axioms.  `s(12) ≥ 3.968616` (1736 points) needs a 328k-leaf tree; not built (estimate below).

**Rung 2, the first zero-margin rung, done (2026-09-27):** `SquarePacking.s13_ge_4 : (4 : ℝ) ≤
minSide 13` and `SquarePacking.s13_eq_4 : minSide 13 = 4` (`S13Lower.lean`), from the case-free
3,621-point cover `certificates/rung2/s13_closed_cover_4.txt` (`notes/s13-casefree.md`), `#print
axioms` = `[propext, Classical.choice, Quot.sound]`, no `sorry`, no `native_decide`.  At side 4 the
bound is sharp (the 4×4 tiling), so monotone witness leaves provably cannot do it (`search/RUNG2.md`
§2); the tree uses the zero-margin leaf types of the new sibling verifier `ZMTree.lean` (below),
proved sound once (`ZMTree.sound`, in the default build).  The data is **opt-in**: `lake build
Sqpack.S13Lower`, 538 s wall / 33 CPU-min on 4 cores, ≤ 13.2 GB RSS per process.

## Files

| file | what |
|---|---|
| `Sqpack/BoxTree.lean` | the generic verifier `BoxTree.check` (natural-number arithmetic only), its soundness `BoxTree.sound`, the gluing lemmas `Cov.splitX/Y/U`, the tree decoder `BoxTree.dec`, and the end-to-end `BoxTree.le_minSide` |
| `scripts/gen_boxtree.py` | builds a box tree for a certificate (exact mirror of `check`) and writes the Lean data |
| `Sqpack/S12U/{Pts,Part0..3,Cov}.lean` | generated: the 81 points; 124 chunk theorems; `cov_root` |
| `Sqpack/S12Lower.lean` | `s12_ge_35_9` (hand-written, one `le_minSide` call) |
| `Sqpack/S12W/*`, `Sqpack/S12WLower.lean` | `s(12) ≥ 3920/997` from `certificates/s12_lower_3.931795_sparse.txt` |
| `Sqpack/S11/*`, `Sqpack/S11Lower.lean` | `s(11) ≥ 3040/797` from `certificates/s11_lower_3.8143.txt` (opt-in) |
| `Sqpack/ZMTree.lean` | the zero-margin verifier `ZMTree.check` (tree type `ZT`: leaves `Z`, `E`; nodes `X/Y/U`, `XM/YM/UM`, `F`, `C`), its soundness `ZMTree.sound` (for the same `BoxTree.Cov`, so `BoxTree.le_minSide` is reused), the decoder `ZMTree.dec` |
| `scripts/gen_zmtree.py` | zero-margin tree search (`search/zeromargin.py` as a read-only oracle) + exact integer mirror of `ZMTree.check` + leaf pruning + Lean emission |
| `Sqpack/S13/{Pts,Part0..3,Cov}.lean`, `Sqpack/S13Lower.lean` | generated: the 3,621 points (a `PTree` and the literal `ptsL`); 209 chunk theorems; `cov_root`.  `s13_ge_4`, `s13_eq_4` (opt-in) |

Regenerate: `python3 lean/scripts/gen_boxtree.py certificates/s12_uniform_7of81_3.888.txt --n 12
--name S12U --outdir lean/Sqpack/S12U` (25 s; deterministic).  The others:
`… s12_lower_3.931795_sparse.txt --n 12 --name S12W --outdir lean/Sqpack/S12W` (60 s) and
`… s11_lower_3.8143.txt --n 11 --look 0 --parts 24 --name S11 --outdir lean/Sqpack/S11` (4 min).
Rung 2: `python3 lean/scripts/gen_zmtree.py certificates/rung2/s13_closed_cover_4.txt --n 13
--name S13 --outdir lean/Sqpack/S13 --nproc 4` (156 s wall / 580 s CPU on 4 cores; deterministic —
a fresh run reproduces the files byte for byte; gitignored, `lean/scripts/gen_data.sh S13`).
Large data sets (`S11`, `S13` today, all future big ones) are **gitignored**: run `lean/scripts/gen_data.sh` before
`lake build Sqpack.S11Lower`.  Small ones (`S12U`, `S12W`) stay committed; both regenerate byte-identically.

## What is proved, and how

`BoxTree.le_minSide`: let `pts` be entries `(X, Y, w)` (point `(X/D, Y/D)`, weight `w/W`) with

* no repeated entry (`PTree.chainB`, kernel), invariant under `x ↦ Mq − x` and `x ↔ y`
  (`d4Check`, kernel), total `Σw < n·W` (kernel);
* `Cov … root`: every closed unit square inside `[0, Mq/D]²` with centre in `[0, Mq/(2D)]²` and
  angle `θ = 2 arctan u`, `u ∈ [0, 29/70]` (`29/70 > tan 22.5°`), contains entries of weight `≥ W`.

Then `Mq/D ≤ minSide n`.  The D4 reduction (`D4.lean`: `d4_reduction`), the cover-to-packing
argument (`S32.lean`: `not_packs_of_cover`, `packs_grid`, `Packs`, `minSide`) and the aggregation
of entries (`Cover.lean`: `coverA`, `coverW`, `D4Inv_cover`, `sum_filter_coverA`) are reused
unchanged.

`Cov` is established by `BoxTree.sound`: `check … t box c = true → Cov … box` for every tree `t`.
A tree splits the pose box `[x0,x1]×[y0,y1]×[u0,u1]` (centre over `Q = D·S`, `u` over `R`) at
midpoints (`X`, `Y`, `U` nodes), prunes the candidate list (`F` nodes: `near`, a pure speed
device), and ends in leaves `L sel`.  A leaf is accepted iff `u1 ≤ R` and

1. **walls**: `WL = ⌊Q(R² + 2U0R − U1²) / (2(R² + U1²))⌋` is a lower bound for `w(θ)/2 = (|cos θ| +
   |sin θ|)/2` on the bin (`wlo_le`; `w = (1 − u² + 2u)/(1 + u²)` by `wid_two_arctan`); the
   admissible centres of the box lie in the clipped rectangle `[max(x0,WL), min(x1, M−WL)] × …`
   (`sq_subset_box_iff`).  If it is empty the leaf holds vacuously;
2. **containment**: otherwise the selected candidates (`sel` = gaps into the candidate list) must
   each pass `ptOk` and weigh `≥ W` (`capSel`).  `ptOk` certifies `p ∈ sq c θ 1` for every
   admissible pose of the box: `p ∈ sq c (2 arctan u) 1 ⟺ G_k ≤ 0, k = 0..3`
   (`mem_sq_iff_gval`); `G_k` is affine in the centre with signs fixed on `u ∈ [0,1]`, so the
   worst centre is a corner of the clipped rectangle (Lemma A); and quadratic in `u`, so on
   `[U0, U1]` it is bounded by its Bernstein coefficients `G(U0)`, `G(U1)`, `G(U0,U1)` (polar
   form) — `lin3`, from the identity
   `(U1−U0)² G(v) = (U1−v)² G(U0) + 2(v−U0)(U1−v) G(U0,U1) + (v−U0)² G(U1)`.
   The 3 × 4 tests are linear in the point with all negative terms moved across (`triOk`), so
   everything is `Nat.add/mul/ble` — kernel-GMP operations, no `Int`, no `ℚ`.

**Trigonometry** is handled by the tan-half-angle parametrisation: `cos θ = (1−u²)/(1+u²)`,
`sin θ = 2u/(1+u²)` exactly, so no trig bounds, no approximation of `π`, and every test is a
polynomial inequality in rational `u`.  The only real-analytic facts used are `cos/sin (2 arctan
u)` (`ZeroMargin.lean`) and `θ ∈ [0, π/4] ⇒ tan(θ/2) ≤ 29/70` (`exists_u_of_theta'`).

**Nothing about the tree is trusted.**  `sound` holds for every tree, so the generator, the
chunking and the decoder `dec` (a tree ships as one hex numeral, a base-`B` digit stream) need no
proofs; a wrong tree only makes `decide` fail.  The certified statement is `Cov` of the root box,
glued from the chunk theorems by `Cov.splitX/Y/U` (a proof term the generator writes out).

## Engineering notes (what it took to make the kernel fast enough)

* **`Nat` only.**  The first version used `ℤ` with `decide (a ≤ b)`: 40 leaves/s and 25 GB.
  Rewriting every test as `Nat.ble (Σ positive terms) (Σ negative terms)` gave ~12×.
* **One declaration per chunk** (≤ 600 leaves).  The kernel's whnf cache lives per declaration;
  one 33k-leaf `decide` needed 32 GB and 167 s, the same tree in 144 declarations 7 GB.  Chunks
  are also the unit of parallelism: the generator splits them over `--parts` files that `lake`
  builds concurrently.
* **Leaves carry their selection.**  Testing every candidate at every leaf (and pruning the
  candidate list at every node) cost ~2×; now a leaf names its ≥ 7 points and only those are
  tested, and pruning happens only where it shrinks the list by ≥ 15 %.
* **Trees as numerals.**  Elaborating a 600-leaf nested constructor term cost about as much as
  its kernel check (~1 ms/leaf, mostly `OfNat` numerals in the leaf lists); one hex numeral per
  chunk, decoded in the kernel, costs ~13 % of the kernel time instead.
* Hand-writing the recursion with `BT.rec` instead of structural recursion: < 5 %, not kept.
  Pre-scaling the points (dropping `X·S`): 3 %, not kept.

## Measurements (2026-09-27, 16-core machine, pinned to physical cores 12–15)

Timings are CPU time on one core unless stated; "kernel" is the `decide +kernel` time, measured
as the difference against the same file with the decisions replaced by an axiom.  RSS includes
~6.6 GB for importing Mathlib.

| certificate | points | tree leaves (splits) | depth | points claimed / leaf | kernel throughput | build |
|---|---|---|---|---|---|---|
| `s12_uniform_7of81` (35/9) | 81, uniform | 29,529 (29,528) | 34 | 7 | **430 leaves/s** (7,594 leaves: 17.6 s) | clean `lake build` 41 s wall / 99 s CPU on 4 cores; 8.2 GB RSS per file |
| `s12_lower_3.931795_sparse` | 224, weighted | 25,127 (25,126) | 35 | ~21 | 170 leaves/s (6,459 leaves: 38 s) | 96 s wall / 251 s CPU on 4 cores; 10 GB RSS |
| `s11_lower_3.8143` | 680, weighted | 89,689 (89,688) | 37 | ~85 | 47 leaves/s (11,412 leaves: 240 s) | 26 min wall / 64 min CPU on 3 cores (oversubscribed); 17.6 GB RSS for an 11k-leaf file |
| `s12_lower_3.9686` (not built) | 1736, weighted | 328,275 | 39 | ~150 (est.) | ~25 leaves/s (est.) | ~4 h CPU (est.) |

Tree generation (Python, one core): 25 s, 60 s, 3.5 min, 29 min.  Elaboration of the data is now
negligible (one numeral per chunk); the kernel is > 85 % of a data file's build time.

**What dominates.**  The kernel does roughly 10⁶ `Nat` operations per second (GMP-backed; the
numbers here are 60–110 bits: `Q ≈ 3·10⁷`, `R ≈ 7·10⁷`, products of three or four of them).  Rational
sizes are not the issue (everything is scaled to integers once; no gcds, no `ℚ`), and there are no
trigonometric bounds at all (tan-half-angle, see above).  A leaf costs about

  `1.0 ms` (decode, walls, clipping, tree traversal) `+ 0.23 ms × (points claimed)`

— a claimed point is 3 Bernstein coefficients × 4 violation polynomials = 12 linear tests.  This
fits all three rows (2.6, 5.8, 20.6 ms predicted vs 2.3, 5.9, 21 ms measured).  So cost scales with
*leaves × points per leaf*, and points per leaf is `1 / (typical weight)`: uniform certificates
are cheap, dense weighted covers expensive.  Memory grows with the file (the 24-file split of
S11 kept each process under ~12 GB); keep ≲ 5k leaves of dense certificates per file.

**Extrapolation to the zero-margin rungs** (made before rung 2, with `BoxTree`'s leaf test; superseded
by the measurements in the rung-2 section below — s(13) took 33 CPU-min with the zero-margin leaves):

| target | boxes | points per leaf (≈ 1/avg weight) | est. kernel CPU | on 4 cores |
|---|---|---|---|---|
| s(13) = 4, 3621-point cover, avg weight 0.0036 | 16,872 | ~280 | 16,872 × 65 ms ≈ 18 min | ~5 min |
| s(32) = 6, 13,085 points, avg 0.0024 | 164k | ~410 | 164k × 95 ms ≈ 4.3 h | ~1.1 h |
| s(21) = 5, 7,536 points + 1,872 segments | 10⁵–10⁶ | ~250 + segment terms | 2–20 h | 0.5–5 h |

These are feasible, but the per-leaf point count is the lever: a leaf only needs *some* certified
subset of weight ≥ 1, and claiming the heaviest points first (as the generator now does) helps
little on these covers (6 % fewer digits on S12W).  Bigger wins, in order: (1) aggregate points
into clusters that a leaf can claim as one unit (a precomputed "sub-cover" whose points all lie in
a small disc: one containment test per cluster via the disc's bounding box); (2) a Mathlib-free
data layer, so each chunk file imports only the `Nat` checker (import time and the 6.6 GB baseline
disappear; soundness stays in the Mathlib file); (3) larger `F` pruning thresholds for dense
certificates (the candidate walk is ~5 % now).

## Rung 2: zero-margin leaves (`ZMTree.lean`), `s(13) = 4`

**Why new leaves.**  `BoxTree`'s only leaf is a *monotone witness set* (points in every admissible
square of the box).  At a sharp container side that is provably insufficient: any finite closed
subdivision certified that way forces total weight `≥ m²` (`search/RUNG2.md` §2, Theorem 2).  So
`ZMTree` adds the leaf types of `search/zeromargin.py`, each an exact integer test with its own
soundness lemma (reusing `ZeroMargin.lean`: `mem_sq_iff_gval`, `gval_eq`, `le_maxQuad`,
`condPoly_eq_gval`, `qeval4_le_maxBern`, `xnR_spec`/`xnW_spec`, `widU`), and states its result as the
same `BoxTree.Cov`, so the final step (`le_minSide`, D4, cover ⇒ packing) is unchanged.

| node / leaf | what the kernel checks | soundness |
|---|---|---|
| `E` | `c₁ < w(θ)/2` (in `x` or `y`) at both ends of the bin (`wgt`); `w` is quasi-concave on `u ≥ 0` | `E_cov` (`widU_ge_min`) |
| `C us` (`clip_bin`) | `w(us)/2 ≥ min(x₁,y₁)/Q` and `(R+us)(R+U₁) < 2R²` (`w` increasing on `[us,U₁]`): poses with `u > us` are inadmissible; child box `[U₀, us]` (often the degenerate bin `us = U₀`) | `C_cov` |
| `Z` | claimed entries (gaps into the candidate list), each with a tag; chains A, B of pivots `(X, Y, kind)`; an emptiness staircase | `Z_cov` |
| `X/Y/U`, `XM/YM/UM`, `F` | midpoint splits, splits at an explicit coordinate (the non-dyadic `1/10` root grid), candidate pruning | `Cov.split*` |

**Inside a `Z` leaf** (signed integers are pairs of naturals `(p, n)`; every test is `Nat.ble`):

* `ADM` (Lemma A with the walls, `_adm_cond_ok_int`): condition `k` of a point at its own centre
  corner, the lower centre bound being the box side (`'R'`: exact quadratic maximum `qOk`, i.e.
  `maxQuad`) or the wall `w(θ)/2` (`'W'`: a quartic `cpoly` = `Q R⁴ condPoly`, degree-4 Bernstein
  `bOk`).  An entry with tag `4` has all four conditions by `ADM` — the witness set `T`.
* `CHAIN` (Lemmas E–H): `pairOk l p q` = `λ_a G_p + λ_b G_q ≤ 0` on the whole box (four corners,
  exact quadratic maximum), `(λ_a, λ_b) = (1,−1)` (chain comparison / down reason), `(1,1), (2,1),
  (1,2)` (up reason, empty product region).  A chain's consecutive pivots are checked monotone; at
  any pose the pivots then cut the box into regions `r = 0..k` (`regions`).  An entry with swing kind
  `kp < 4` has the other three conditions by `ADM` and condition `kp` from a *down* reason
  `G_p ≤ G_{q_d}` (counts in regions `r ≥ d`) and/or an *up* reason (counts in `r < u`), in chain
  A and/or chain B.  Regions `(r, s)` of the product with `r < ka`, `s < e_r` are certified empty by
  `pairOk` between `A_{r+1}` and `B_{e_r}`; every other region must reach `W`.  One chain is the
  special case `kb = 0`, pure `ADM` the case `ka = kb = 0`.

The generator uses `zeromargin.py` (read-only) as the oracle for which leaf and which witnesses,
recomputing `T` without `P1` and without inherited points (the Lean leaf has neither), and checks
every leaf with an exact integer mirror of `ZMTree.check`; on s(13) the mirror accepts every leaf
zeromargin proposes, so the Lean tree is essentially zeromargin's own D4 tree (8,434 boxes against
its 8,452; `ADM` 1,416 and `EMPTY` 1,738 identical, `CHAIN` 2,663 against 2,672 — the clip value is
rounded to the `2⁻³²` grid).  It then prunes each chain leaf: pivots are dropped greedily while all regions still reach `W` (reasons remapped by
transitivity, `λ` unchanged), then unneeded entries and reasons.  That takes the chains from 94
pivots (one chain) / 76 × 65 (product) on average to about 3 / 3 × 5, and the claimed entries from
2.23 M to 1.33 M.

**The s(13) tree.**  400 root cells (`1/10` pitch, D4 region `[0,2]²`) × 8 `u`-bins of `[0, ½]`;
10,024 search boxes, max depth 13.  Leaves: **`E` 1,738; `Z` 4,079 = `ADM` only 1,416 + one chain
1,133 + two chains 1,530**; 1,590 `C` nodes, 2,617 midpoint splits in the cells, 374 `F` nodes.
Claimed entries 1,325,532 (1,122,165 `ADM` witnesses, 203,367 chain entries, 16,717 pivots), 325 per
`Z` leaf; 1.80 M base-2²⁰ digits in 209 chunks, 4 files.  Scale `Q = 1000·4096`, `R = 2³²`.

**Kernel cost** (cores 12–15): `lake build Sqpack.S13Lower` **538 s wall, 1,995 s CPU**, the four
part files 491–505 s each in parallel, **13.2 GB max RSS** per process (6.6 GB of it the Mathlib
import); `Pts` 23 s, `Cov` 3 s.  The kernel is ~98 % of a part file (elaboration alone: 7.8 s).
That is ~0.47 s CPU per `Z` leaf, or ~1.4 ms per claimed entry.
Measured on the heaviest chunk (32 `Z` leaves, 8,914 `ADM` entries, 825 chain entries): an `ADM`
entry costs ~0.5 ms (0.2 ms of it the test itself), a chain entry ~2.5 ms, the region count ~1.5 s
per chunk.  What it took (same chunk, 27.4 s → 7.6 s kernel):

* **Fast paths** that imply the exact tests, so every verdict is the exact one: `ADM` first by
  `ptOkK` (`BoxTree.ptOk` with one kind skipped: 12 linear tests on shared products; 0.2 ms against
  0.45 ms for four `qOk`), pair tests first by degree-2 Bernstein on per-box triples (`cornFast`),
  and a *single* corner when the combination does not depend on the centre (`cfree`: `G_p − G_q`
  of equal kinds, `G_p + G_q` of opposite kinds — most down/up reasons).
* **One numeral per `Z` leaf** instead of one per chunk: digit extraction from a 300-kbit numeral
  cost 85 µs per digit (2.1 s per chunk); now negligible.
* **The point list as a literal** (`ptsL`, `pts_toList : pts.toList = ptsL` proved once):
  evaluating `pts.toList` cost 0.75 s in *every* chunk declaration.
* The `ADM` weight is summed once per leaf; regions iterate the chain entries only.
* `lake` builds files concurrently only with `LEAN_NUM_THREADS > 1` (with `=1` it serialised the
  four parts).

**Extrapolation** (per-`Z`-leaf cost × leaf count × entries per leaf ∝ 1/average weight):

| target | `Z` leaves (zeromargin D4 census) | entries/leaf vs s(13) | est. kernel CPU | on 4 cores |
|---|---|---|---|---|
| s(13) = 4 (measured) | 4,079 (ADM 1,416, CHAIN 2,663) | 1 | 33 min | 9 min |
| s(32) = 6, 13,085 points | 82,208 (ADM 12,201, CHAIN 70,007) | ~1.5 | ~17 h | ~4–5 h (16 files) |
| s(21) = 5, 7,536 points + 1,872 segments | — | — | needs segment entries | — |

s(32) needs no new leaf type (same `ZMTree`, `m = 6`, root grid `[0,3]²` by `XM/YM` at pitch
`1/10`), only the generator run (zeromargin's own D4 sweep was 2.8 CPU-h; with the mirror and the
pruning expect ~3×) and the build; the remaining levers are the per-entry overhead in context
(0.5 ms against 0.2 ms for the bare test) and a Mathlib-free data layer.  Not yet run on s(32):
whether the mirror accepts every zeromargin leaf there without `P1` / inherited points (it did on
s(13)) is the first thing to measure; a rejected leaf is split further, never accepted.  s(21)
additionally needs segment entries (`MixedMeasure.lean`) in `Cov` and in the `Z` count.

## Remaining gaps / next steps

* No `sorry`; nothing is assumed about the tree, the generator or the decoder.
* `BoxTree`'s leaf test is sound everywhere but a tree can only terminate where the certificate
  has slack; the zero-margin leaves are in `ZMTree` (rung 2 above, `s(13) = 4` done).  Next:
  `s(32) = 6` with the same leaf types (estimate above); segments (`s(21)`) need `MixedMeasure` in
  `Cov` and a segment entry type in the `Z` leaf.
* D4 symmetry is assumed (`d4Check`); a certificate without symmetry needs a second root lemma
  (centre in `[0,m]²`, `u ∈ [0,1]`, via `θ ↦ θ + π/2`), about 30 lines on top of `d4_reduce`'s
  step D.  All certificates in `certificates/` checked so far are D4-invariant.
* Others' certificates (rung 2) only need a file in this format (`certificates/FORMAT.md`).

