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

## Files

| file | what |
|---|---|
| `Sqpack/BoxTree.lean` | the generic verifier `BoxTree.check` (natural-number arithmetic only), its soundness `BoxTree.sound`, the gluing lemmas `Cov.splitX/Y/U`, the tree decoder `BoxTree.dec`, and the end-to-end `BoxTree.le_minSide` |
| `scripts/gen_boxtree.py` | builds a box tree for a certificate (exact mirror of `check`) and writes the Lean data |
| `Sqpack/S12U/{Pts,Part0..3,Cov}.lean` | generated: the 81 points; 124 chunk theorems; `cov_root` |
| `Sqpack/S12Lower.lean` | `s12_ge_35_9` (hand-written, one `le_minSide` call) |
| `Sqpack/S12W/*`, `Sqpack/S12WLower.lean` | `s(12) ≥ 3920/997` from `certificates/s12_lower_3.931795_sparse.txt` |
| `Sqpack/S11/*`, `Sqpack/S11Lower.lean` | `s(11) ≥ 3040/797` from `certificates/s11_lower_3.8143.txt` (opt-in) |

Regenerate: `python3 lean/scripts/gen_boxtree.py certificates/s12_uniform_7of81_3.888.txt --n 12
--name S12U --outdir lean/Sqpack/S12U` (25 s; deterministic).  The others:
`… s12_lower_3.931795_sparse.txt --n 12 --name S12W --outdir lean/Sqpack/S12W` (60 s) and
`… s11_lower_3.8143.txt --n 11 --look 0 --parts 24 --name S11 --outdir lean/Sqpack/S11` (4 min).
Large data sets (`S11` today, all future big ones) are **gitignored**: run `lean/scripts/gen_data.sh` before
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

**Extrapolation to the zero-margin rungs** (leaf test as here; the zero-margin germ/chain leaves
of `ZeroMargin.lean` are *not* yet in this verifier and would add per-leaf cost):

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

## Remaining gaps / next steps

* No `sorry`; nothing is assumed about the tree, the generator or the decoder.
* The leaf test is sound everywhere but a tree can only terminate where the certificate has
  slack; zero-margin
  (`s(13)`, `s(32)`, `s(21)`) needs the germ/chain leaf types of `ZeroMargin.lean` as additional
  `BT` constructors, each with its soundness lemma — the gluing (`Cov.split*`) and the final
  theorem `le_minSide` are reusable as they are.  Segments (`s(21)`) need `MixedMeasure` in `Cov`.
* D4 symmetry is assumed (`d4Check`); a certificate without symmetry needs a second root lemma
  (centre in `[0,m]²`, `u ∈ [0,1]`, via `θ ↦ θ + π/2`), about 30 lines on top of `d4_reduce`'s
  step D.  All certificates in `certificates/` checked so far are D4-invariant.
* Others' certificates (rung 2) only need a file in this format (`certificates/FORMAT.md`).

