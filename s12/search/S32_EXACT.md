# The exact check of the `s(32)` candidate: first pilots (2026-09-25/26)

Goal: `s(32) = 6` from `runs/s32-close_candidate.txt` (13,085 points of `[0,6]²`, `D = 1000`, `W = 10¹¹`, total
`3171350535386/10¹¹ = 31.713505354 < 32`; built in `S32_COVER.md` §9).  If `verify2/target/release/zmcheck cert`
returns `VERIFIED` on it, that is the whole proof: the reduction (cover ⇒ packing bound) is `notes/s13-casefree.md`
§1 / Lean `packing_le_weight`, and `s(32) ≤ 6` is the 6×6 tiling.  Then an independent re-check (`zeromargin.py`,
as for `s(13)`), shipping, and a literature check that `s(32) = 6` is open.

Labels: **certified** (exact, `zmcheck`, partial sweeps only so far), **heuristic** (float sampling).
Nothing below is a full sweep; no verdict on `s(32)` yet.

## 0. Answer so far

* **Everything except the interior tile germs certifies.**  A full column (`c_x ∈ [3.4, 3.5]`, all `c_y`, all angles,
  480 roots): 476 roots clean, **68.8 CPU-h** in total.  It crosses 6 tile centres, both wall bands and the
  grid-line families.
* **All remaining failures sit at interior tile germs**: poses within `0.003` of a tile centre at `|θ| ≤ 0.56°`.
  In the column: 60 uncertified boxes at depth 22, in the 4 roots touching `(3.5, 1.5)` from below and `(3.5, 4.5)`
  from above (mirror images under `y ↦ 6 − y`, as the cover's symmetry predicts), 13.6 CPU-h.  The germs
  `(3.5, 0.5)`, `(3.5, 2.5)`, `(3.5, 3.5)`, `(3.5, 5.5)` close from this side.
* **None of the failures is a coverage hole (heuristic).**  Every uncertified box so far was float-sampled
  (3–40k random poses + corners each): minima `1.010` (column germs), `1.011` (germ cell, branch-cap 640),
  `1.009 / 1.014` (off-germ, default caps), `1.033` (the `×1.0195` cover).  They are checker limits.
* **Two checker settings matter:** `--branch-cap 640` (default 160) and, at germs, `--node-cap 4000000 --depth 22`.
  Both only widen the search, so they are sound.
* **Estimate of a full sweep with these settings:** about 12 of the 60 columns lie on the germ lines `x = k + ½`
  and cost like this one; the rest are much cheaper (off-germ `≈ 700 CPU-s` per `0.1 × 0.1` cell).  Roughly
  **1,000–1,500 CPU-h (≈ 2 days on 28 threads)**, plus whatever the germs need.  A D4 symmetry reduction, as
  `zeromargin.py` already does, would give ≈ 4×.

## 1. Why it was slow: the branch cap (default 160)

`perf` (via `nix shell nixpkgs#perf`) on the first pilot: **94 % of the time in `Disj::search`**.  `ZM_DEBUG` on
an uncertified box shows the mechanism.  Near `θ ≈ 0` with the square's edge on a grid line, the monotone witness
set `T` weighs just under 1 (e.g. `0.99921`), and the forced two-chain order walks **136 pivots**: the points on
the grid line under the moving edge.  The branch list is truncated at `branch_cap = 160`, so the pivots that
would finish the proof (the opposite edge's line) never enter.  With `--branch-cap 320` or `640`, the same box
certifies in **3 nodes, 36 ms** (vs ~900 ms and failure).  The same held for tilted (`≈ 43°`) off-germ boxes.

Point structure (why `m = 4` never hit the cap):

| | `s(13)` cover | `s(32)` candidate |
|---|---|---|
| points / density | 3,621 / 226 per unit² | 13,085 / 363 per unit² (interior 461) |
| median point weight | 0.0021 | 0.0015 |
| points carrying 50 % of the weight | 554 | 1,880 |
| points on one grid-line unit segment next to a tile centre | 163 (weight 0.488) | 152 (weight **0.339**) |

The chain lengths are similar; the weight per chain step is lower, so the searches need more combinations.

## 2. The germs

At an interior tile centre with `θ → 0`, all four edges of the square lie on grid lines at once.  The
disjunctive search then has to combine two long chains of line points, and their product grows faster than the
node cap or depth can follow: `T` drops to `0.75–0.92` and one box took 447k regions with `--node-cap 4000000`.
The germ cell `x ∈ [3.4, 3.6], y ∈ [2.5, 2.7]` (32 roots):

| settings | wall (threads) | boxes | max depth | uncertified |
|---|---|---|---|---|
| default | > 2 h (28), killed | — | — | — |
| `--branch-cap 640` | 2,845 s (20) | 4,894 | 18 | 48, all in the 4 `θ ≈ 0` roots |
| `--branch-cap 640 --node-cap 4000000 --depth 22` | 3,963 s (28), 6.4 CPU-h | 3,468 | 18 | **0** |

The same settings leave the germs `(3.5, 1.5)` and `(3.5, 4.5)` open in the column run (§0).  Germ roots cost 2–4.4
CPU-h each on one thread.  With 16 interior germs (the 4 × 4 tile centres off the wall row; edge-row germs such as `(3.5, 0.5)` close cheaply) × up to 8 roots each, the germs are the dominant and least predictable cost.

## 3. Other pilots

| run | cell | settings | result |
|---|---|---|---|
| `s32zm_tilecell` | binding edge tile `x ∈ [3.4, 3.6], y ∈ [5.4, 5.6]` | default | 59 s (6), 678 boxes, depth 11, **0 uncertified** |
| `s32zm_offgerm` | `x ∈ [3.1, 3.3], y ∈ [2.1, 2.3]` | default | 353 s (8), 1,088 boxes, 2 uncertified (`≈ 43°`; both close with `--branch-cap 640`) |
| `s32zm_intcell_x10195` | germ cell | `r22 × 1.0195` (total 31.983), default | 5,331 s (20), 35,930 boxes, 2,008 uncertified: more margin did **not** help |
| `s32zm_intcell_nc20k` | germ cell | `--node-cap 20000` | no faster; killed |
| exact `pose` at the confirmed float dip `(3.4996, 5.4925, u = 0.000785)` | | | `1.0137` (the rounded pose misses the thin cell) |

## 4. Next (Evan decides)

1. **Point shifting at the germs** (Evan's suggestion, recommended first): fewer, heavier grid-line points next to
   the tile centres (merge, or move a little weight onto them), then re-check validity with the float scans
   (`s21_cover_eval.py hardscan` / `close_shifted.py confirm`).  Test case: the `(3.5, 1.5)` germ, whose roots are
   `--xlo 3.4 --xhi 3.4 --ylo 1.4 --yhi 1.4`.  The budget is `32 / 31.3716 = 1.020` over the `r22` weights.
2. A germ-specific primitive in `zmcheck` (one tilted edge crossing one grid line is a 1-D problem), with its
   soundness lemma.
3. Brute force on the germ roots only: `--depth 26`, a larger `--node-cap`.
4. D4 reduction in `zmcheck` (4×), keeping an unreduced run or `zeromargin.py` as the independent check.
   **Done (8× in roots):** `zmcheck cert --d4` and the runner `search/s32_sweep.py`, §7.

## 5. Code

`verify2/src/main.rs`: instrumentation only, gated by `ZM_ROOTLOG` (stderr): one `ROOT … boxes / maxdepth / leaf
census / uncert / seconds` line per root, and a `HEARTBEAT` every 60 s (roots done, boxes, depth histogram).
Verdicts and stdout are unchanged: the tile-cell census reproduces exactly (678 boxes, depth 11, ADM 99 / DISJ 184
/ EMPTY 72).  The pilots ran a build from a scratch `CARGO_TARGET_DIR` (same source); `verify2/target/release/zmcheck` was
rebuilt from it at commit time and reproduces the tile-cell census.

## 6. Files (`runs/`, gitignored)

`s32zm_{pilot,tilecell,intcell,intcell_x10195,intcell_nc20k,intcell_v2,intcell_bc,germ2,offgerm,col34}.{log,err,sh}`,
`s32zm_cells/` (the aborted per-cell sweep, 600 s caps), `s32_r22x10195.txt` (`r22 × 10195/10000`, total 31.983301).

```
Z=verify2/target/release/zmcheck      # rebuilt with this commit
ZM_ROOTLOG=1 $Z cert runs/s32-close_candidate.txt --threads 28 --branch-cap 640 --node-cap 4000000 --depth 22 \
    --xlo 3.4 --xhi 3.5 --ylo 2.5 --yhi 2.6                   # germ cell: 0 uncertified, 3963 s
ZM_ROOTLOG=1 $Z cert runs/s32-close_candidate.txt --threads 20 --branch-cap 640 --node-cap 4000000 --depth 22 \
    --xlo 3.4 --xhi 3.4                                       # column: 60 uncertified, all at germs, 68.8 CPU-h
```

## 7. The D4-reduced sweep (`zmcheck cert --d4`, `search/s32_sweep.py`)

The candidate is exactly D4-symmetric (all 13,085 points match under `x ↦ 6−x`, `y ↦ 6−y`, `x ↔ y`,
checked 2026-09-26 on the integer coordinates; `--d4` re-checks it in every run).  So one eighth of the root
boxes suffices.

**What `u` is.**  A pose is `(x, y, u)` with `u = tan(θ/2) ∈ [0, 1]`, `θ ∈ [0°, 90°]`, and the square is
`Q(c, θ) = c + R(θ)[−½, ½]²` (`zmcheck` tests `|X|, |Y| ≤ ½` with `(X, Y) = R(−θ)(p − c)`).  Root boxes:
`x, y` cells of pitch `1/10` over `[0, 6]` (cell `i` is `[i/10, (i+1)/10]`), `u` in 8 bins `[k/8, (k+1)/8]`;
`--xlo/--xhi/--ylo/--yhi` select cells by their lower-left corner.

**How D4 acts.**  Every `g ∈ D4` (about `(3, 3)`) maps the axis square `[−½, ½]²` to itself, so
`g Q(c, θ) = Q(gc, θ')` with `θ' ≡ θ` (mod 90°) for the rotations and `θ' ≡ −θ ≡ 90° − θ` for the four
reflections.  In `u`: rotations keep `u`, reflections send `u ↦ (1−u)/(1+u)`, which does not preserve the
`u`-bins (bin 0 ↦ `[7/9, 1]`); the centre grid *is* preserved (3 is a grid line).

**Fundamental region** (`--d4`): centres in the quadrant `[0, 3]²` (cells `i, j = 0..29`) and `u ∈ [0, ½]`
(bins 0..3): **3,600 of the 28,800 roots**.  A reflection brings any pose to `θ ≤ 45°`, i.e.
`u ≤ tan 22.5° = √2 − 1 < ½`; a rotation, which keeps `θ`, then brings the centre into the closed quadrant.
`u ≤ ½` is `θ ≤ 53.13°`, 18 % more angle than needed, but it keeps the bin boundaries (the same choice
`zeromargin.py` makes for its 4× reduction).  The alternative with the stock binary (centre triangle
`y ≤ x ≤ 3`, all 8 bins, per-column `--ylo`) needs 3,720 roots and has 20 interior-germ roots instead of 16,
so the flag was worth it.

**Soundness.**  Let `f(p)` be the total weight at `p` (duplicates summed), `W(Q) = Σ_{p∈Q} f(p)`.
(i) `f∘g = f` for all `g ∈ D4`: checked exactly by `--d4` on the generators `x ↦ 6−x` and `x ↔ y`
(otherwise `ERROR:`, exit 2).  (ii) `g` is an isometry with `g[0, 6]² = [0, 6]²`, so it maps the closed unit
squares contained in `[0, 6]²` onto the same family (admissibility is preserved), and
`W(gQ) = Σ_{q∈Q} f(gq) = W(Q)`.  (iii) Every admissible pose `(c, θ)` has a `g` with `g(c, θ) = (c', θ')`,
`c' ∈ [0, 3]²`, `θ' ∈ [0°, 45°]` (above).  (iv) The region's root boxes are *closed*, their union is
exactly `[0, 3]² × [0, ½] ⊇ [0, 3]² × [0, √2−1]`, and `zmcheck` certifies every pose of each closed box, so
box boundaries need nothing beyond the unreduced sweep.  Hence if all 3,600 roots are certified,
`W(Q) = W(gQ) ≥ 1` for every closed unit square `Q ⊆ [0, 6]²`.  One `zmcheck cert FILE --d4` without range
flags prints `VERIFIED-D4:`; the per-job runner concludes from the census of all 3,600 roots instead
(partial runs never say VERIFIED).

**End-to-end test** on the shipped `s(13)` cover (`certificates/rung2/s13_closed_cover_4.txt`, D4-invariant):
`zmcheck cert … --depth 18 --threads 4 --d4`: 1,600 roots, `done in 400s: boxes 4164, max depth 10`,
`ADM 552 DISJ 1424 EMPTY 906 UNCERTIFIED 0`, `VERIFIED-D4:` (unreduced: 12,800 roots, 2,015 s at 4 threads,
30,258 boxes).  The candidate with one weight changed by 1 is refused (`ERROR: --d4: certificate is not
D4-invariant …`).

**Code** (`verify2/src/main.rs`; the verdict logic is unchanged): the `--d4` flag; `check_d4` (exact,
a map of aggregated weights); the root loop takes cells `0..5m−1` and bins `0..3` under `--d4`; a
`VERIFIED-D4:` line for a clean unrestricted `--d4` run; and, instrumentation only, under `ZM_ROOTLOG` one
`UNCERT <box> IN <root>` stderr line per uncertified box (≤ 80 per root), so a job's uncertified boxes are
all collected, not just stdout's first 40.

**Runner.**  `search/s32_sweep.py`, run from `s12/`; needs a `zmcheck` built with `--d4` (rebuild
`verify2/target/release/zmcheck`, or pass `--zmcheck PATH`):

```
python3 search/s32_sweep.py plan                          # 180 jobs: 30 columns × 6 blocks of 5 y-cells, 20 roots each
python3 search/s32_sweep.py run --jobs 7 --threads 4      # -> runs/s32d4/, resumable, germ-heavy jobs first
python3 search/s32_sweep.py summary                       # -> runs/s32d4/SUMMARY.txt
#   --only REGEX (job names cII_yAA-BB), --zmcheck PATH, --out DIR, --yblock B (fixed per DIR at first use)
```

Each job is `ZM_ROOTLOG=1 zmcheck cert runs/s32-close_candidate.txt --d4 --threads T --branch-cap 640
--node-cap 4000000 --depth 22 --xlo a --xhi a --ylo b --yhi c`.  A job's `.log` appears (renamed from
`.log.part`) only when it finished, so a rerun skips finished jobs and restarts interrupted ones.  The summary
checks each log (root count, the `D4:` and `D4-REDUCED` lines, depth 22), sums boxes and CPU (from the `ROOT`
lines), records the binary's sha256 per `run`, lists every uncertified box with its root, and prints
`D4 SWEEP CLEAN` only when all 3,600 roots are in and none is uncertified.

**Sanity test** (2026-09-26, 4 threads).  The tile cell `x ∈ [3.4, 3.6], y ∈ [5.4, 5.6]` maps by the half-turn
to cells `i ∈ {24, 25}, j ∈ {4, 5}`: `s32_sweep.py run --out runs/s32d4_test --yblock 2 --only 'c2[45]_y04-05'`,
2 jobs, 16 roots, 111 s + 94 s wall, 208 CPU-s, 342 boxes, max depth 11, **0 uncertified** (the unreduced
pilot: 32 roots, 354 CPU-s, 678 boxes, depth 11).  A rerun skips both.

**Germ roots in the region**: `u`-bin 0 and the four cells around each tile centre.  (The `θ → 0⁻` side of a
germ, in bin 7, is the `θ → 0⁺` side of a mirror germ, also in the region.)

| germ | orbit | cells `(i, j)`, bin 0 | status so far |
|---|---|---|---|
| `(1.5, 1.5)` | inner corners (4) | `(14–15, 14–15)` | never run |
| `(1.5, 2.5)` | inner middles (8) | `(14–15, 24–25)` | `(14, 25)` is the image of the failing column-34 roots `y 1.4, bin 0` and `y 4.5, bin 7`; `(15, 25)` certified there |
| `(2.5, 1.5)` | inner middles | `(24–25, 14–15)` | `(25, 14)` is the image of the failing roots `y 4.5, bin 0` and `y 1.4, bin 7` |
| `(2.5, 2.5)` | centre (4) | `(24–25, 24–25)` | images of the germ cell around `(3.5, 2.5)`: certified (§2) |
| wall row `(0.5, k+½)`, `(k+½, 0.5)` | | cells 4–5 in one coordinate | cheap (the sanity test is `(2.5, 0.5)`) |

So at these settings expect uncertified boxes at least in `(14, 25)` and `(25, 14)`, bin 0 (`|θ| ≤ 0.56°`,
within `0.003` of the germ); `(1.5, 1.5)` is unknown.

**Cost estimate.**  Column 34's bins 0..3 (29.2 of its 68.8 CPU-h) are exactly the region's row 25 ∪ column 25
(by `(x, y) ↦ (y, 6−x)` and the half-turn), so one germ line costs ≈ 14.6 CPU-h.  The 8 germ lines
(`i` or `j ∈ {14, 15, 24, 25}`): ≈ 8 × 14.6 − 25 (the 16 germ cells counted twice) ≈ **90 CPU-h**.  The other
441 admissible cells (cells `< 5` are EMPTY) at 350–2,000 CPU-s each (off-germ ≈ 700 CPU-s per 8-bin cell;
more on the integer lines): ≈ **60–130 CPU-h**.  Total **≈ 150–250 CPU-h** (unreduced: 1,000–1,500), i.e.
≈ 6–9 h wall on 28 threads (`--jobs 7 --threads 4`), plus the tails of single germ roots (2–4.4 CPU-h each).

## 9. Performance (`zmcheck` profiling, 2026-09-26)

Branch `worktree-agent-a9f0f2716bab1a8db` (rebased on `0fa5f1c`, not merged): 5 code commits that change no verdict,
plus `ZM_UBINS` (benchmarking only: `ZM_UBINS=0,2` keeps those angle bins of each root cell and makes the sweep
partial) and `verify2/bench/bench.sh`.  Timings were taken with 20–44 other threads running on the 32-thread
machine, so CPU-seconds vary by ±15 % or more.  User-mode instruction counts (`perf stat`) do not depend on
load, so they are the figures to compare.

**Benchmarks.**  *tile*: the cell `x ∈ [3.4, 3.6], y ∈ [5.4, 5.6]`, default settings, 32 roots, at 1 and 4
threads.  *germB* and *germ0*: the single roots `x 3.4, y 2.3`, `u`-bins 2 and 0, with the germ settings
(`--branch-cap 640 --node-cap 4000000 --depth 22`).  germB is dominated by branch ranking; germ0 runs the
long forced two-chain searches.

**Where the time went.**  In the baseline, 96 % of the time was in `Disj::search`, and about 90 % of the total
was in `weight_of`.  That is the per-point loop that sums the weight of a cover: it walks one `Vec<u8>` per
swing point, and the ranking calls it twice per branch candidate (up to 640) at every search node.  The
checker is **CPU-bound**: IPC 4.9, L1d miss rate 23 % (almost all L2 hits), negligible LLC misses, RSS under
20 MB.  Threads do not contend: the tile cell takes 128 CPU-s on 1 thread and 123 CPU-s on 4.  The only
shared state is one atomic root counter.

| commit | change | measured effect |
|---|---|---|
| `3bbf237` | word-parallel `weight_of`: point `fi` is nibble `fi % 16` of word `fi / 16`, so the test is whether that nibble of `need & !cover` is zero | germB 70.5 → 32.3 CPU-s; tile 123 → 70 |
| `228780b` | ranking as `weight(cover) +` per-word deltas over the nonzero words of each hypothesis mask; per-word weights from 4-point subset-sum tables (`i64` when the swing points' total fits, else the old loop) | germB 32 → 15; tile 70 → 37 |
| `7456201` | no clones of the hypothesis masks in `cover_child`; a free list reuses `Cov` buffers | germ0 −13 % instructions |
| `5403218` | `Ĝ_q` computed once per hypothesis in `exact_mask`; no per-point `Vec` clones | tile −4 % instructions, −9 % cycles |
| `b4e33e2` | `Cov` carries its per-word weights and total, updated only where a mask changes a needed bit; no full `weight_of` per node | −2 … −6 % instructions on all three |

**Total (baseline `4d49633` vs `b4e33e2`; each pair ran at the same time under the same load):**

| benchmark | baseline | now | CPU | instructions |
|---|---|---|---|---|
| tile, 1 thread | 149–173 CPU-s, 3,092 G instr | 23–29 CPU-s, 362 G | **5–7.6×** | 8.5× |
| tile, 4 threads | 160–216 CPU-s, 61–68 s wall | 27–28 CPU-s, 8–9 s wall | **5.7–7.9×** | 8.5× |
| germB | 80–112 CPU-s, 1,640 G | 11.5–11.7 CPU-s, 152 G | **6.8–9.7×** | 10.8× |
| germ0 | 719 CPU-s, 2,738 G cycles, 8,622 G instr | 410–433 G cycles at the same load, 1,154 G instr | **≈ 6.5×** | 7.5× |

On an idle machine, expect roughly 6–8× less CPU for the sweep.  §7's D4 estimate of 150–250 CPU-h would then
be about 25–40 CPU-h.

**Why the verdicts are unchanged.**  Every change computes the same integers, only in a different order or
with fewer redundant steps.  The following checks all passed:

* The census and `ROOT` lines are identical on every benchmark.
* The sorted `--dump` leaf files are **identical**: every DISJ leaf's region count, node count and seed, and
  every ADM witness list.  This holds for the tile cell, germB and germ0, so the search visits exactly the
  same nodes.
* Builds with `debug-assertions` and `overflow-checks` assert every new computation against the old one on
  tile, germB and germ0: the per-point `weight_of` reference, delta against full sum, and the incremental
  `Cov` against the union and weight recomputed from scratch.
* The `i128` fallback path gives the identical census when all weights and `W` are scaled by 10¹⁰.
* The rung-2 rejection tests pass 23/23, with output byte-identical to the baseline's.
* `--d4` root selection is unchanged (`main` against the branch head on two cells).

**Tried and not kept.**

* `target-cpu=native`: slower with the lookup tables.  A branch-free masked sum only matches the tables with
  AVX-512, so it is not portable.
* `codegen-units = 1`: no change in instruction count, and cycles within noise.  Adding `panic = "abort"`
  gave −3 to −5 % cycles on the tile cell, which is inconclusive under this load.  `bench.sh` can recheck it
  on an idle machine.
* A "newly covered" bit loop for the deltas: 2× slower (branch mispredictions).
* Evaluating the Lemma-I tests on corner polynomials (5 coefficients instead of 15; the same integers, even
  modulo 2¹²⁸): −5 % instructions, no gain in cycles.
* Sparse masks stored as structure-of-arrays: slower.

**Where the time is now.**

* germ0-type roots: the cover bookkeeping, mostly `or_hyp` and the Lemma-K closure loop in `cover_child`
  (about 50 %), then per-word weights.
* germB-type roots: the ranking's `weight_delta`.
* The tile cell: about 40 % in the exact `i128` Lemma-I tests (`exact_mask`) and 35 % in ranking.

IPC is 2.7–3.4, the L1d miss rate 2–15 %, and LLC misses are negligible: still core-bound.

**Algorithmic ideas (not done).**

1. Schedule germ roots first, in `s32_sweep.py` or in the root order.  A 2–4 CPU-h germ root that starts
   last sets the wall time.
2. Rank fewer candidates.  Every node scores all ≤ 640 candidates, but only the top `rank_exact = 64`
   matter.  An exact upper bound per candidate (`pc.w` plus the weight of the points its mask touches) could
   skip most of them and still give the same top 64, with ties broken by index.
3. An `i64` fast path for the exact tests at shallow depths, backed by a per-box proven magnitude bound (the
   bits of `Dc·M²` and the like) or by checked arithmetic falling back to `i128`.  This targets the ~40 %
   off-germ share.
4. Inherit exact hypothesis masks from the parent box.  A Lemma-I certificate on a box holds on every
   sub-box, as the ADM `cmask` already uses.  This is sound but not bit-identical: it can change which leaves
   certify.
5. Make the forced-chain covers incremental along the chain instead of re-closing Lemma K over the whole
   `active` list at each node.

**Benchmark script.**  `verify2/bench/bench.sh` builds each commit from `git archive` into its own target
directory.  It runs tile (1 and 4 threads), germB and germ0 `-n` times each, one run at a time and never more
than 4 threads.  It records user+sys CPU time, wall time, max RSS and `perf stat` counters (cycles,
instructions, IPC, cache and L1d misses).  It diffs every build's census against the baseline's and prints
medians and speedups.  To run it on an idle machine:

```
cd s12 && CERT=$PWD/runs/s32-close_candidate.txt verify2/bench/bench.sh -n 3     # ≈ 1.5–2 h; -q (no germ0) ≈ 35–45 min
```
