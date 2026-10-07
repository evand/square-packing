# Mixed covers (points + segments) in the zero-margin Lean verifier (2026-09-28)

Brief: `tasks/lean-segments/README.md` (batch 1: the mathematics — soundness lemmas and a generator
validated on toy covers; no big kernel builds).  Everything here is kernel reduction only
(`decide +kernel`), no `sorry`, no new axioms: every end theorem prints
`[propext, Classical.choice, Quot.sound]` (`lean/Axioms.lean`).

## 1. What is done

| target | status | where |
|---|---|---|
| 1. `CovM` + gluing + end-to-end + D4 for segments | **done** | `Sqpack/CovM.lean`: `CovM`, `CovM.split{X,Y,U}`, `CovM.of_cov`, `segD4Check`, **`le_minSide_mixed`** (via `MixedCover.measure_apply`, `MixedCover.d4InvM`, `d4_reduction_measure_u`, `not_packs_of_measure`) |
| 2. Lemma P, Lemma S, PIECE leaf, `sound` for the extended tree | **done** | `Sqpack/ZMTreeM.lean` (S-blocks, T-groups, `pc_sound`), `Sqpack/ZMTreeX.lean` (tree `ZTM`, `checkM`, **`soundM`**, decoder `decM`); `ZMTree.lean` untouched |
| 3a. Lemma T | **done** (matching form) | `SegParts.lean` `pair_capture` (+ `gval_pairV/H`); `ZMTreeM.lean` T-groups |
| 3b. Lemma L + "general ends" + Corollary L | **done** (box corners, no wall corners) | `LemmaL.lean`, `LBlock.lean`, `LBlockSound.lean` (**`lblk_sound`**) |
| 3c. Lemma L′ (short chord, hull form), Lemma V | not done | §5 |
| 4. Lemma R / `SPLIT` | not done | §5 |
| 5. Generator | **done** for S, T, L and points (Lemma P with `ZMTree` point leaves) | `lean/scripts/gen_zmmtree.py`, `lean/scripts/lblock.py`, toys `lean/scripts/mk_toys.py` → `lean/toys/` |

Toy end theorems (all `[propext, Classical.choice, Quot.sound]`):

* `SquarePacking.s3_ge_2 : (2 : ℝ) ≤ minSide 3` (`S3Lower.lean`, default build) — pure segments (`M2`:
  lines `x = 1`, `y = 1` of `[0,2]²` at density `5/8`), 178 `Z` leaves (81 S-only, 97 with an L-block).
* `SquarePacking.s3_ge_2_mixed` (`S3LowerP.lean`, default build) — points + segments (`M2p`, total 2.85):
  282 `Z` leaves, 222 of them `ADM` point leaves with the piece bound as phantom (Lemma P).
* `SquarePacking.s16_ge_4 : (4 : ℝ) ≤ minSide 16` (`S16Lower.lean`, **opt-in**, data gitignored) — the
  `ZM_MIXED.md` §4.3 grid cover T1 (six grid lines of `[0,4]²`, density 5/8, margin 3.55 %): 9,451 `Z`
  leaves (8,256 with an L-block, 1,742 with T-groups at the tile germs).  (Trivial bound; full exercise of
  T at germs and L off them.)

## 2. The Lean design, lemma by lemma (statement ↔ what the kernel checks)

**`CovM`** (`CovM.lean`).  `W ≤ Σ_{points in Q} w + Σ_segments w · segFrac(Q)`, `segFrac` the Lebesgue
fraction of the segment's parameter interval mapped into `Q` (`MixedMeasure.lean`).  `le_minSide_mixed`:
D4-invariant data (`d4Check` for points, `segD4Check` for segments: sorted keys, normalised, both
reflections present with equal masses), `pts.wsum + segs.wsum < n W`, `CovM` of the root box
`[0, Mq/(2D)]² × [0, Um/R]` with `R ≤ 2 Um` ⇒ `Mq/D ≤ minSide n`.

**Parts** (`SegParts.lean`).  A part `(e, l, h)` is a sub-interval of an axis-parallel segment whose
interior lies in `Q`.  `parts_le_segMass`: parts pairwise non-overlapping on each entry certify
`Σ w (h − l)/|e| ≤ segMass`.  Every block below produces parts at each pose; the claims (segment + tag)
make parts of different blocks live on different entries.

**Lemma P** (`ZMTreeX.checkM`).  The leaf computes the piece value `Lp` (units of `1/W`) and runs the
unchanged `ZMTree.zOk` with target `W − Lp` (`Nat.sub`; `Lp ≥ W` is the `PIECE` leaf).  So the phantom is
exactly zeromargin's: present in every region, never tested geometrically.

**Lemma S** (S-block `(dir, K, A, B)`).  Kernel: `ZMTree.admAll` (all four conditions) at the two ends
`(K, A)`, `(K, B)`; each tagged segment contributes `⌊w |seg ∩ [A,B]| / |seg|⌋`.  Soundness: at a pose
each condition is a half-plane in the point (`gval` affine), so the certified set of the line is convex
(`certified_interval`).  *Versus zm_mixed*: zm_mixed takes `J4 = ∩_k conv ∪_χ I_{k,χ}` (Bernstein
intervals per bound choice χ); the kernel needs both ends of `[A,B]` to pass some χ for every condition.
An end of `J4` can lie in a gap between two χ-intervals of another condition (audit note); the generator
rounds `J4` inward to the `1/Q` grid, checks with the mirror and bisects inward if needed.  On the toys
this never cost a leaf (M2: census identical to `zm_mixed.py --no-lin`, 376 PIECE / 664 EMPTY).  The
Lean `admK` uses the exact quadratic maximum for the box-side/box-side choice where zm_mixed's
`line_cond_iv` uses degree-4 Bernstein: the kernel is at least as strong there.

**Lemma T** (T-group).  Lines `x = ξ`, `x = ξ + 1` (vertical; singular kinds 1, 0) or `y = η + 1`,
`y = η` (kinds 2, 3).  The whole content of Lemma T used is the identity
`G_up(t) + G_down(t′) = 4u (t′ − t − u)` (`gval_pairV`, `gval_pairH`: the two singular violation
polynomials).  **The Lean form is the dual of Corollary T**: instead of `min_T f(T)` over breakpoints,
the certificate is a monotone *coupling* of the two lines' masses — pairs (up piece `[a,b]`, down piece
`[c,d]`, mass `m`) with `c ≤ a + u₀`, `d ≤ b + u₀`, `m ≤` each piece's mass, pieces sorted — and
`pair_capture` proves that at every pose each pair captures `≥ m` (a threshold argument on one pair;
no minimisation, no `±∞` cases).  Summing, the group captures `≥ Σ m`.  By max-flow/min-cut on the
line (nested neighbourhoods), the best coupling equals `inf_T f(T)`: the generator's greedy FIFO coupling
reproduces zm_mixed's decisions exactly on the T1 germ cell (`[1.4,1.6]²`, S + T only: 36,576 boxes,
18,304 PIECE leaves, depth 11 in both).  Plus: the J3 intervals (three conditions, `admAll` with the
singular kind skipped) and optional points `τ` certified for the singular condition (one point suffices,
by monotonicity along the line: zm_mixed instead requires the certified singular set to be a half-line,
weaker — audit (c)1).  Kernel cost per pair: `O(1)`.

**Lemma L + general ends + Corollary L** (L-block).  Per line: up end `a`, down end `b` (`b ≤ a`),
`Δ↑`, `Δ↓`, roles (for each condition: untyped / typed at the up end / typed at the down end; the kernel
checks the role against the condition's slope along the line, `kend`, and requires `u₀ > 0` for the
`±4u` slopes and `u₁ < 1` for `±2(1 − u²)`), pieces of the three zones, and affine minorants `s x + i`
(`s ≥ 0`, fine units `1/(WQ)`) of the two gains, checked piece by piece against a floor-rounded density
(`gainOk`, `gain_lb`).  Per box corner: an option per end (the cap, or a typed condition's threshold) and
a slack.  At a corner every option is `N(v)/σ̂(v)` with `N` quadratic in `v = Ru`, `σ̂ ∈ {4Rv,
2(R² − v²)}`, so all corner inequalities (slacks, and `Σ(chosen − σ) ≥ lg`) are quartics checked by
`ZMTree.bOk` (degree-4 Bernstein) — no general-degree Bernstein needed.  Soundness (`lblk_sound`):
`chord` (every point of the line between `b − Y` and `a + X` is in `Q`, `X = min(Δ↑, min_k x_k)`,
`x_k = −Q G_k(a)/σ_k`), the gain bound (`up_gain`, `dn_gain`), the per-end minimum over options
`emin`, and `Phi_box`: `Φ = Σ_lines (emin↑ + emin↓)` is concave in the centre (`ConcXY`), so `Φ ≥ lg` on
the box follows from the four corners (`concave_corners`).  **Difference from zm_mixed**:
`lemma_l_joint` uses a containing rectangle whose sides are the box sides *or the wall terms* `w(u)/2`
(float-chosen sub-bins, `corner_choices`); the kernel uses the box corners only (valid: any containing
rectangle; weaker near the walls, and it keeps every corner inequality of degree ≤ 4).  Slacks are per
corner (as in zm_mixed) but not per sub-bin.  On T1 (whole D4 region): zm_mixed without SPLIT 7,054 boxes /
3,514 PIECE leaves / depth 7; the Lean tree 21,519 boxes (incl. 1,811 clips) / 9,451 leaves / depth 11
(×2.7 leaves).  M2: 150 vs 178.

## 3. zm_mixed.py: statement vs code (these matter)

From a line-by-line audit (a forked agent, read-only, with small exact experiments) and from mirroring
the lemmas in the kernel.  **No reachable soundness gap was found**; every certifying decision is
`Fraction`/integer arithmetic, floats only steer or reject.

* **B1 — latent, would be unsound if reached.**  `region_phi` (Lemma R's region bound) returns `'EMPTY'`
  on a zero-width angle bin: with `u0 == u1` the sub-bin loop never runs and the initial `best = 'EMPTY'`
  is returned unchecked (`zm_mixed.py` ~l. 862–905; confirmed on `(1, 11/10, 1, 11/10, u = 1/4..1/4)`).
  `cert_split` treats `'EMPTY'` as a proven-empty region and skips it.  Zero-width bins with `u0 > 0`
  do occur (`clip_bin` returns `u0` at walls).  **Unreachable today** only because on a zero-width bin
  `lemma_l_joint` and `vertex_split` also return None, so `lparts = None` and `cert_split` is never called;
  the fallback `if lparts is None: lparts = ([], 0, Lbox)` in `cert_split` is dead code that would certify
  with nothing checked if reached.  A Lean mirror of `SPLIT` must not accept an empty region without a
  proof (e.g. require `U0 < U1` for emptiness).
* **Weaker than needed (incompleteness only).**  Typed Lemma L inequalities must have `I_k` of the exact
  form `(−∞, α]` / `[β, ∞)`, else the line gets no Lemma L; Lemma T's `S↑` must be `[t↑, ∞)`, else `τ = +∞`.
  Monotonicity along the line makes one certified point enough (the Lean blocks use one point).  On
  zero-width bins Lemma L, V and SPLIT are all unavailable (exactly the wall boxes `clip_bin` produces).
  `lemma_l_data` gives up when both hull slopes are `≤ 0`.  SPLIT and Lemma V refuse `u0 = 0`.
* **Sound but different from the note.**  Corollary L does not need slopes `≥ 0` (a minimum of affine
  functions is concave for any slopes; they are `≥ 0` anyway).  `Δ` in Lemma L is any positive number
  (the code's `move_bound` only tightens).  Lemma V: any `D` is sound; the code also drops the other
  short-chord lines.  Lemma R's `λ ∈ {1, ½, 2}` is `(a, b) ∈ {(1,1), (2,1), (1,2)}` (equivalent).  The
  float-chosen wall/side corner and split points of `corner_choices` are sound for any choice.  The
  reach filter never raises a bound (drops whole segments; all bounds monotone in the kept set; the
  Lemma L′ majorant is a bound for the kept sub-measure).  `Lpar` inheritance is valid (children are
  sub-boxes of the clipped parent); children inherit `L`, not the points-moved `L2` (correct).  The
  phantom never reaches `P1` (no `inh` there), and `cert_split` excludes it.
* **What the Lean leaf does differently, by design.**  Lemma T in matching form (equal value up to
  `1/Q` rounding); S-block ends tested by `admK` (exact quadratic maximum where zm_mixed uses Bernstein);
  per-segment flooring of S-block masses (`≤ 1/W` per segment) and of densities in L-blocks (fine units
  `1/(WQ)`); no `P1`, no inherited witnesses in the point part (as `gen_zmtree.py`).

## 4. Measurements (cores 12–15)

Kernel time per `Z` leaf (one-leaf chunks, difference against the same file with the decisions
replaced by `sorry`; T1 data): **S-only 23 ms, T-group 27 ms, L-block 125 ms** (avg 3 lines: ~40 ms per
L-line).  M2 part (178 leaves, 15 chunks): 20.8 s kernel.  Point entries cost as in `ZMTree` (s(32): 1.06 ms
per claimed entry).  Generator (Python, mirror included): T1 21,519 boxes in 283 CPU-s (13 ms per box).

**T1 built in full** (`s16_ge_4`, opt-in): 6 part files, ~420 s each (≈ 42 CPU-min, 15 min wall on
3 cores, ≤ 11 GB RSS per part).

**s(21) pilot (generator only, the interior germ cell `[1.4,1.6]²`, pitch 1/10, 8 u-bins, 32 roots):**
the mirror verifies it: 3,646 boxes, depth 15, 1,839 `Z` leaves (`ADM` 1,440, one chain 339, two
chains 60; 1,263 with T-groups, 576 with L-blocks), 230 claimed points per leaf, 292 CPU-s (80 ms per
box).  zm_mixed's shipped run needed 5,470 boxes (depth 17, pitch 1/20) for this cell.  Kernel, 30
one-leaf chunks each: **0.69 s per leaf (T + points), 0.89 s per leaf (L + points)** — the points
dominate.

## 5. Open obstacles

1. **SPLIT (Lemma R)** is not in the Lean leaf.  It needs the Lemma L bound minimised over a region
   polygon whose extra vertices (pivot line ∩ box edge) are rational functions of `u` with denominators
   `2(1 − u²)` or `4u`: at such a vertex the option values have degree up to ~6–7, so the kernel needs a
   general-degree Bernstein lemma (`P(U0 + H t) = Σ C(n,i) β_i tⁱ(1−t)ⁿ⁻ⁱ`, integer coefficients
   `C(n,i)β_i = Σ_j C(n−j, i−j) b_j`), plus the region logic of `ZMTree`'s chains with a per-region piece
   value.  The same general-degree lemma would allow the wall corners of Corollary L.
2. **Lemma L′ (short chord, hull form) and Lemma V** are not in the Lean leaf (vertex lines get their
   core only).  L′ is a variant of the L-line with a majorant `ℓ₂` at the down end (options of the form
   `ℓ₁(min(r↑, B) − p) − ℓ₂(max(r↓, A) − p)`), all still quartic at box corners; V needs a split
   constraint (degree as in 1).
3. **Kernel cost of L-blocks** (40 ms per line) is ~10× my operation count; not yet profiled
   (decoding was ruled out).  Candidates: repeated evaluation of `capU`/`optT` inside `cornerOk`,
   `List.getD` lookups of claims.
4. **s(21)** (not attempted, per the brief): the Lean tree will be larger than zm_mixed's without SPLIT,
   L′/V and wall corners; see the report for the extrapolation.
