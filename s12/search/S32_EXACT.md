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

## 8. `zeromargin.py` as the independent re-check (2026-09-26, 2 cores)

**Verdict: it works at `m = 6` and closes everything it was given, including a germ root that `zmcheck`
leaves open.  It is too slow for a full sweep as it stands: an estimated 1,500–10,000 CPU-h (most likely
3,000–5,000), against 25–40 CPU-h for `zmcheck --d4` after §9.**

**Code.**  Nothing depends on `m = 4`: `m` is read from the certificate, the constants (`0.7072`, `14143/10000`)
are geometric, and the weights (`W = 10¹¹`, sum `3.2·10¹²`) fit its `int64` exact sums.  There are no branch or
node caps.  CHAIN builds complete monotone chains and tries every product of two of them, with fixed
`λ ∈ {1, ½, 2}`.  The only knobs are `--depth` (the default of 16 is too low; roots reached 16–18, so use
22–24) and `--chain-from`.  Its symmetry reduction is **4×, not D4**: it uses `x ↦ 6−x` and `y ↦ 6−y` only,
with `u ∈ [0, ½]` and `c_y ≤ 3`, which is twice §7's region.  Small plumbing added: `--cy-lo/--cy-hi` and
`--progress N`.  The defaults are unchanged: Friedman-14 gives a census identical to the previous file.
Timing driver: `runs/s32py_drv.py`, which runs the roots one per task, logs each root, and applies a per-root
time cap.  It is a timing tool, not a verifier.

**`--full` is broken near walls.**  `clip_bin` does not clip bins above 45°.  So in unreduced mode, the
`u ∈ [7/8, 1]` roots at the top wall of the tile cell do not converge: after 540 s they had 1,318 boxes,
against 27–133 for `zmcheck`.  The reduced mode, which the `s(13)` check used, never has `u > ½`.  Use only the
reduced mode.

**Measured on 1–2 cores.**  Each `zmcheck` cell was mapped to its images in the reduced domain.  The tile cell
maps to `x ∈ [2.4,2.6] ∪ [3.4,3.6]`, `y ∈ [0.4,0.6]`.  The off-germ cell maps to `x ∈ [2.7,2.9] ∪ [3.1,3.3]`,
`y ∈ [2.1,2.3]`.  The column-34 cell `y ∈ [2.2,2.3]` maps to `x ∈ [2.5,2.6] ∪ [3.4,3.5]`.  `zmcheck` below is
the pre-§9 binary; the current one is about 6–8× faster.

| cell | `zeromargin.py` | `zmcheck` (pre-§9) | ratio |
|---|---|---|---|
| tile | 464 CPU-s, 356 boxes, depth 9, **0 uncert.** | 230 CPU-s, 678 boxes | 2× |
| col. 34, `y 2.2` (germ line, not a germ) | 5,717 CPU-s, 720 boxes, depth 12, **0** | 597 CPU-s (640/4M/22) | 9.6× |
| off-germ | ≥ 20,150 CPU-s, ≥ 4,174 boxes, depth 18, **0** in the 61 finished roots; 3 roots stopped at the 1 h cap (708–823 boxes each) | 573 CPU-s, 1,088 boxes, 2 uncert. at default caps | ≥ 35× |
| germ root `x 3.4, y 1.4`, `u ∈ [0, ⅛]` | 9,244 CPU-s, 798 boxes, depth 16, **0 uncert.** | 12,441 CPU-s, 735 boxes, depth 22, **13 uncert.** (640/4M/22) | 0.74× |

The off-germ ≈ 43° box that needs `--branch-cap 640` sits in the root `x[3.2,3.3] y[2.2,2.3] u[3/8,7/16]`.
There it closes at default settings: 315 boxes, depth 17, 1,366 s.  One of the depth-22 germ boxes that
`zmcheck` leaves open, `x[3.49922,3.5] y[1.49922,1.5] u[0,1/2048]`, gives `zeromargin.py diag` → CHAIN in
12 s.  The 2-cut product of 181 × 209 pivots has worst region 1.0099; `T = 0.692`.  **So the uncapped
two-chain product is exactly what the germs need.**  It could also serve as the model for `zmcheck`'s germ
handling: §4 item 2, or dropping the caps at the germs.

**Where the time goes.**  99.9 % is in `cert_chain`, at 2–15 s per call with 13,085 points.  About 700 points
are in reach and 780 are swing candidates.  Each candidate costs a binary search of `Fraction` `_gmax`
evaluations (4 corners × up to 3 `λ`).  ADM, P1 and MIX together take under 0.1 % of the time.

**Full-sweep estimate** (reduced domain, 8 bins of `u ∈ [0, ½]`, about 1,350 admissible cells).  Take §7's
pre-§9 estimate: 90 CPU-h on the germ lines plus 60–130 CPU-h elsewhere.  Double it for the 2× larger region.
Scale by the measured ratios: about 1–10× on germ lines and 10–35× elsewhere.  That gives **1,500–10,000
CPU-h**, most likely 3,000–5,000.  That is 4–7 days on 28 threads, or roughly 100× the current `zmcheck --d4`.

**What would make it practical** (none of this was done):
1. Exact integer arithmetic in `_gmax`/`cert_chain`.  All inputs have denominators `1000 · 2^k`, so common
   denominators and Python `int` would do; `gmpy2` is not installed.  Expected 3–10×.
2. A numpy float pre-screen of the chain comparisons (`_gmax ≤ 0`), confirming exactly only near ties, as
   `cert_adm` already does.  Probably another several ×.
3. The diagonal reflection, as in §7: `u ↦ (1−u)/(1+u)` does not keep the bins, so use §7's region
   `[0,3]² × [0,½]` and check `x ↔ y` exactly.  Gives 2×.
4. `--depth 22` or more, and a per-root runner like `s32_sweep.py`.  The roots are embarrassingly parallel.

With 1–3 the sweep would be roughly 100–300 CPU-h, a few hours on the machine, which is practical as the
independent check.  No cap needs raising.

Files (`runs/`): `s32py_drv.py`, `s32py_{tilecell_red,offgerm_red,col34y22,germ15,germ15b}.{log,err}`,
`s32py_germdiag22.log`, `s32py_tilecell_full.{log,err}` (the `--full` blow-up).

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

## 10. Result: the candidate is certified by two exact checkers (2026-09-26)

**Claim (certified, pending review and an independent re-check):** every closed unit square in `[0,6]²` captures weight
`≥ 1` from `runs/s32-close_candidate.txt` (sha256 `a0d2d38f…`, total `31.713505354 < 32`), hence `s(32) = 6`.

* **zmcheck, D4 sweep** (§7; `search/s32_sweep.py summary`): all 3,600 roots of `[0,3]² × u ∈ [0,½]`, `--d4`
  (symmetry checked exactly), `--branch-cap 640 --node-cap 4000000 --depth 22`.  53,660 boxes, 80.2 CPU-h (6 jobs on
  build `67e760bf`, 174 on the §9 build `2043cff8`; the two builds give identical ROOT lines on the 40 roots of
  `c14_y10-19`, `runs/s32d4_xcheck/`).  **3,595 roots clean; 5 roots with uncertified boxes** (52 boxes), all `u ∈ [0,⅛]`:
  `x∈[1.5,1.6] y∈[1.4,1.5]`, `x∈[2.4,2.5] y∈[1.4,1.5]`, `x∈[2.5,2.6] y∈[1.4,1.5]`, `x∈[1.4,1.5] y∈[2.4,2.5]`,
  `x∈[1.4,1.5] y∈[2.5,2.6]`: the interior germs `(1.5,1.5)`, `(2.5,1.5)`, `(1.5,2.5)`.
* **zeromargin.py on those 5 roots** (`runs/s32py_germfix.sh`, driver `runs/s32py_drv.py` → `Checker.run_box`,
  `use_chain=True`, `--depth 24`, same `u = tan(θ/2)`): each root as its two half-bins `u ∈ [0,1/16], [1/16,1/8]`;
  **10/10 clean, 0 uncertified**, max depth 18, 7.6 CPU-h (`runs/s32py_germfix_*.log`).
* **Gluing.** Each checker certifies closed pose boxes; the 3,595 + 5 roots are exactly the `--d4` fundamental region,
  and `--d4`'s soundness argument (§7) needs only that every box of that region is certified by *some* exact checker.

Open before announcing: (1) review of the gluing and the two primitives used at germs; (2) an independent re-check of
the 3,595 zmcheck roots (zeromargin.py needs the §8 speedups: ~100–300 CPU-h), or a Lean replay; (3) a clean
single-binary rerun and a `certificates/s32/` bundle with a verify script.  Alternative certificate: `S32_SHIFT.md`
(`s32_shift_v1.txt` + `ZM_MIXPAIR`), stock-zmcheck germs still running.

## 11. The independent re-check: `zeromargin.py` over the whole D4 region (2026-09-26)

**Verdict (certified, second exact checker): `D4 RECHECK CLEAN`.**  `zeromargin.py` (independent of `zmcheck`: its own
primitives ADM / P1 / MIX / CHAIN, its own subdivision, no caps) certifies all **7,200 roots** of `[0,3]² × u ∈ [0,½]`
(pitch `1/10`, 8 bins of `u` of width `1/16`) on `runs/s32-close_candidate.txt` (sha256 `a0d2d38f…2144`), with the
D4 invariance of the *weighted* point set checked exactly.  7,199 roots certify at depth limit 24; one root needs
depth 30 (below).  164,130 boxes, max depth 27, leaves ADM 12,201 / CHAIN 70,007 / EMPTY 3,457 / UNCERTIFIED 0,
**2.77 CPU-h, 12 min wall** (6, then 28 processes).  Checker `runs/s32py_d4/checker/zeromargin.py`, sha256
`640fe453c1a32f4aa580ca2b1261c6406923a4d7c131f65604a432c7fc2086ab` (= `search/zeromargin.py` at this writing).
With §10 this makes the whole region certified twice, by two independently written exact checkers.

**The one root that needed more depth.**  `x ∈ [1.5,1.6], y ∈ [0.5,0.6], u ∈ [0,1/16]` (the wall-row tile germ
`(1.5, 0.5)`): at depth 24 one box is left, `cx ∈ [1.50117, 1.50156], cy ∈ [0.50195, 0.50234]`, `θ ∈ [0.258°,
0.269°]` (a `clip_bin` sliver along the wall-tight poses `c_y = w(θ)/2`).  `diag`: ADM ∪ P1 weight 0.917, the best
chain region 0.99863 (short by 0.0014).  Float sampling (heuristic; `runs/s32py_regress/fsample.py`, 40,000 random
admissible poses + corners): minimum **1.01178**, so a checker limit, not a hole.  The frozen original gives the same
depth-24 result (133 boxes, the same box, 115 CPU-s against 4 s).  The whole root rerun at `--depth 30`: **133 → 139
boxes, max depth 27, certified**.  The runner's `--deepen 30` pass records that as a second, complete run of the root.
`zmcheck` certified this root (it is not among §10's five).

**Cross-check without the diagonal reflection** (`runs/s32py_c2half.{sh,log}`): the other half of the old
C2×C2-reduced domain, `c_x ∈ [3,6]`, `c_y ∈ [0,3]`, `u ∈ [0,½]` (7,200 roots, CLI, `--depth 30 --theta-bias 1`):
164,992 boxes, max depth 27, ADM 12,453 / CHAIN 70,206 / EMPTY 3,437 / **UNCERTIFIED 0**, 430 s wall on 24 processes.
Together with the D4 run (`c_x ∈ [0,3]`) this is the whole domain of the old reduction (`x ↦ 6−x`, `y ↦ 6−y` only),
so the `s(32)` claim also holds by `zeromargin.py` without the new diagonal-reflection argument.

### 11.1 Profile (cProfile, before any change)

One off-germ root that certifies in one box (`x ∈ [3.1,3.2], y ∈ [2.1,2.2], u ∈ [0,1/16]`, 9.9 CPU-s): **97 % in
`_gmax`** (50,117 calls, each 4 corners × a Fraction quadratic), i.e. CHAIN's pair tests.  After the CHAIN change,
the next item was the Fraction ADM condition tests (`_adm_cond_ok`, 54 % of the rest), called for every swing
candidate.  Everything else was < 10 %.

### 11.2 Changes and measured speedups

One process each, `--depth 24`, same settings as the runner, the machine loaded (±15 %).  "int CHAIN" is (a)+(b)
below; "+ int ADM" adds (c).

| root | boxes | frozen / `--ref` | int CHAIN | + int ADM (final) | total |
|---|---|---|---|---|---|
| `x[3.1,3.2] y[2.1,2.2] u[0,1/16]` (1 box, off-germ) | 1 | 9.85 s | 0.17 s | 0.08 s | 123× |
| `x[3.4,3.5] y[0.5,0.6] u[0,1/16]` (wall tile) | 23 | 39.3 s | 1.15 s | 0.74 s | 53× |
| `x[2.4,2.5] y[0.5,0.6] u[0,1/16]` (wall tile germ) | 65 | 81.3 s | 2.85 s | 1.79 s | 45× |
| `x[1.5,1.6] y[1.4,1.5] u[1/16,1/8]` (interior germ) | 39 | 395.9 s | 6.79 s | 3.79 s | 104× |
| regression sample, 29 roots (below) | 357 | 1,471 s | | 18.1 s | 81× |
| §10's 10 germ half-roots (frozen logs) | 2,470 | 27,535 s | | 245 s | 112× |
| `s(13)` cover, reduced domain, `--depth 18 --disj` | 16,872 | 6,138 s wall on 8 procs | | 103 s wall on 6 procs | ≈ 80× |

(a) **Exact integer arithmetic for CHAIN's pair tests** (`_box_ctx`, `_gle0`, `_quad_le0`).  Per box, one common
denominator `Dc` for the points' coordinates and the box corners, one `Du` for `u`: at a corner, `Dc·G` has integer
coefficients and the test "max over `[u0,u1]` of `c0 + c1 u + c2 u² ≤ 0`" becomes two endpoint signs (scaled by
`Du² > 0`) and, when `c2 < 0` and the vertex is strictly inside (`U0(−2c2) < c1 Du < U1(−2c2)`), the sign of
`4 c0 c2 − c1²`: the same rule as `_max_quad`, so the same booleans.  `λ ∈ {1, ½, 2}` become the integer pairs
`(1,1), (2,1), (1,2)` (a positive rescaling).  Python ints; `gmpy2` is not installed and was not needed (the numbers
stay below ~10⁴⁰, and after (b) the exact calls are no longer the bottleneck).
(b) **Float screen, exact acceptance** (`cert_chain_fast`).  numpy evaluates the same maxima in floats for whole
vectors of pairs.  A float value `> 10⁻⁹` rejects a pair (it is skipped: never unsound); `|v| ≤ 10⁻⁹` is decided
exactly; `< −10⁻⁹` only steers.  The chain is built on a float pair matrix per swing kind, but each pivot is appended
only after an exact test.  The three binary searches of the reference (down-sets, up-sets, the empty-region staircase)
run in lockstep over all rows with the reference's own probe sequence, and the final accepted probe of every row is
re-decided exactly (`fconfirm`).  If that confirmation failed, the row would be redone all-exact (`fmismatch`: never
happened in any run).  Region weights are sums over the swing candidates only (`inT` is disjoint from them), as one
int64 matrix product per chain row; same integers as the reference's `enough`.  The chain sort key (a float of an exact
rational) is computed as `int / int`, which is correctly rounded, exactly like `float(Fraction)`, so the sort is
identical.  A pair matrix for the chain build, the sort key and a per-box cache of exact results did not change the
time measurably on the germ root (3.83 → 3.96 s); they are kept because they are simple.
(c) **Integer ADM condition tests** (`_adm_cond_ok_int`, `_cond_poly_i`, `_bern_le0_i`).  Everything times a common
denominator `D` (points, `m`, the box's `R` bounds).  The Bernstein coefficients are times `12·Du⁴`, so
`12·Du⁴·β_i = Σ_j C(i,j)·(12/C(4,j))·(b_j Du⁴)` is an integer with the sign of `β_i`.  Same booleans as the Fraction
Lemma-A test.  Used everywhere ADM is (T, MIX, the swing screen).
(d) **The D4 reduction** (`--d4`, `d4_roots`, `Checker.symmetric_d4`): 7,200 roots instead of the old reduction's
14,400 (**2×**).
(e) Plumbing: the CLI's `Pool` no longer pickles the whole `Checker` (13,085 Fractions) with every root task (fork +
a module global).  Not measured separately.
Not done: PyPy (not installed; the fast path is numpy-heavy anyway); multiprocessing was already per root.

`--ref` (CLI) or `chk.fast = False` gives the old Fraction code path verbatim (`cert_chain_ref`,
`_adm_cond_ok_ref`); `--selfcheck` runs both at every CHAIN call and every ADM condition test and asserts that they
agree.  The unreduced (`--full`) and old reduced paths are unchanged.

**Why the cost fell so far below §8's estimate.**  §8 extrapolated 1,500–10,000 CPU-h for the old 2× region.  The
measured speedup is 45–120× (81× on the random sample), so the old code would have needed roughly
2.77 h × 81 × 2 ≈ 450 CPU-h for its region.  The §8 extrapolation from a few heavy cells was pessimistic: 5,327 of the
7,200 roots certify in their first box.

### 11.3 Soundness of each change

* (a) is the reference's test, multiplied through by positive integers.  At a corner, `a = p_x − c_x = A/Dc` exactly,
  so `Dc·G` is the reference's `G` times `Dc > 0`.  The maximum over the pose box of a positive combination of the `G`s
  (affine in the centre) is still attained at a corner of the centre rectangle, as in the reference.  The endpoint,
  vertex-location and vertex-value tests are exact integer identities for the reference's rational ones (derivation
  in `_quad_le0`'s docstring).
* (b): every *acceptance* on which a certificate rests is an exact integer test.  The certificate needs (i) each
  consecutive chain pair `G_{q_r} ≤ G_{q_{r+1}}` on the box: each is tested exactly before the append; (ii) for a
  down-set entry `r ≥ lo(p)`: `G_p ≤ G_{q_lo}` exactly (or `p = q_lo`), and then `G_p ≤ G_{q_r}` by (i) and
  transitivity; (iii) for an up-set entry `r ≤ lo(p)`: some `λ > 0` with `G_p + λ G_{q_lo} ≤ 0` exactly, and then
  `G_p + λ G_{q_r} ≤ 0` since `G_{q_r} ≤ G_{q_lo}`; (iv) for an empty product region `(r, s)` with `s < e(r)`:
  `G_{A_r} + λ G_{B_e} ≤ 0` exactly, hence `G_{A_r} + λ G_{B_s} ≤ 0`, so `G_{A_r} > 0` and `G_{B_s} > 0` never
  hold together; (v) the weight of each region's witness set is exact int64 arithmetic on `W·weight` (≤ 3.2·10¹²).
  A float that is wrong can therefore only lose a certificate (reject a good pair or stop a search early), never add
  one.  The regions and the per-region witness rule are the reference's (RUNG2.md §§6–7).
* (c) is (a)'s argument for the Lemma-A polynomials: `g·D` has integer coefficients, and the Bernstein coefficients
  `β_i` of `g` on `[u0,u1]` scale by `12·Du⁴·D > 0`.  The Bernstein upper bound itself is the reference's.
* (d) **D4, checked independently of §7.**  `f(p)` = total weight at `p` (duplicates summed); `W(Q) = Σ_{p∈Q} f(p)`.
  (1) `Checker.symmetric_d4` checks exactly (Fractions, on the aggregated map `f`) that `f∘g = f` for
  `g = (x ↦ m−x)` and `g = (x ↔ y)`.  These two reflections generate the dihedral group of the square (their product
  is a quarter turn about `(m/2, m/2)`), so `f∘g = f` for all eight.  (2) Each `g` is an isometry of `[0,m]²` onto
  itself, so it maps the closed unit squares inside `[0,m]²` onto the same family, and `W(gQ) = Σ_{q∈Q} f(gq) = W(Q)`.
  (3) With `Q(c,θ) = c + R(θ)[−½,½]²` and `θ` taken mod 90° (the axis square is invariant under `R(90°)`): a quarter
  turn `ρ` about the centre gives `ρQ(c,θ) = Q(ρc, θ)`.  A reflection `σ` gives `σQ(c,θ) = Q(σc, −θ)`, because
  `σR(θ)σ = R(−θ)` and `σ[−½,½]² = [−½,½]²`.  This holds for both generators: the diagonal reflection also sends `θ`
  to `−θ`.  (4) Given any pose, take `θ ∈ [0°,90°)`; if `θ > 45°`, apply `x ↦ m−x` to get `90°−θ ∈ [0°,45°)`.  Then
  apply the quarter turn (a power of `ρ`) that maps the closed quadrant holding the centre onto `[0,m/2]²`; the four
  quarter turns map `[0,m/2]²` onto the four closed quadrants, which cover `[0,m]²`.  Result: `c' ∈ [0,m/2]²`,
  `θ' ∈ [0°,45°]`, `u' = tan(θ'/2) ≤ tan 22.5° = √2−1 < ½`.  (5) `d4_roots` returns closed boxes whose union is
  exactly `[0,m/2]² × [0,½]` (cells `i, j < (m/2)/pitch`, bins `k/16`); the assertion `(m/2)/pitch ∈ ℤ` guards it.
  `zeromargin.py` certifies every admissible pose of each closed root box, as in the old reduced mode (which already
  used `u ≤ ½` and relies on the same `clip_bin` / EMPTY logic for `u ≤ ½`).  Hence all 7,200 roots certified ⇒
  `W(Q) ≥ 1` for every closed unit square `Q ⊆ [0,m]²`.  §7's argument for `zmcheck` is the same; I rederived it
  rather than reusing it.  The cross-check above does not use the diagonal reflection at all.
* **A gap in the old symmetry check, now closed.**  `Checker.symmetric()` compared only the *point sets* under
  `x ↦ m−x`, `y ↦ m−y`, not the weights.  A certificate with symmetric points but asymmetric weights would have been
  reduced unsoundly.  Demonstrated: the `s(13)` cover with one weight `+1` is reported `symmetric: True` by the frozen
  original, and refused by the new check.  Both shipped covers are symmetric in their weights too (the new checks
  pass on `s13_closed_cover_4.txt` and on the candidate; `zmcheck --d4` checks it independently), so no earlier result
  is affected.  The new `symmetric()` / `symmetric_d4()` compare the aggregated weight maps exactly.  Refusal tests:
  one weight `+1` → refused by the C2×C2 check; `+1` on a C2×C2 orbit of 4 off-diagonal points
  (`(0.8,1), (3.2,1), (0.8,3), (3.2,3)`) → C2×C2 passes, `--d4` refused (`ERROR: … not D4-invariant`, exit 2).

### 11.4 Regression (all before the big run)

* **Friedman-14** (`friedman14 --tri --depth 14`): new = frozen, byte-for-byte in the census: 7,164 boxes, depth 10,
  ADM 3,496 / TRI 74 / EMPTY 3,212 / UNCERTIFIED 0, `VERIFIED`.  (RUNG2.md §4's 7,220 is an older revision.)
* **`s(13)` cover** (`cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --disj --chain-from 0`), reduced
  domain: 16,872 boxes, depth 13, ADM 2,867 / CHAIN 5,320 / EMPTY 3,449 / UNCERTIFIED 0, `VERIFIED`: identical to the
  census recorded in `verify.sh`, in 103 s on 6 processes (was 6,138 s on 8).  With `--d4`: 3,200 roots, 8,452 boxes,
  ADM 1,416 / CHAIN 2,672 / EMPTY 1,738, `VERIFIED-D4`, 90 s on 3 processes.
* **`s(32)` roots against the frozen original** `runs/zeromargin_germfix_c1997272.py` (imported read-only by
  `runs/s32py_regress/regress.py`, the germfix driver's settings, depth 24).  Sample: 20 random admissible roots of the
  D4 region, 8 structural ones (wall corner, EMPTY cells, tile germs, the centre cell) and one interior germ root:
  **29/29 identical** census dicts (boxes, max depth, every leaf count) and uncertified-box lists; CPU 1,471 s →
  18.1 s.  Against §10's frozen logs of the 10 interior-germ half-roots: **10/10 identical**.  The uncertified root
  above: identical at depth 24 (frozen original run separately).  `runs/s32py_regress/compare.txt`.
* **Per-box self-check** (`fast = 'check'` / `--selfcheck`: fast and reference evaluated at every CHAIN call and
  every ADM condition test, asserting equal kind and equal witness list): the whole 29-root sample (357 boxes,
  `runs/s32py_regress/selfcheck.jsonl`), plus 6 further roots (300 boxes) during development: no mismatch.
  `fmismatch = 0` in every run.
* **Box counts**: no change is meant to alter the search except (d), which changes *which* roots are run, not how a
  root is searched: every per-root census above is identical.

### 11.5 Runner: `search/zm_d4_sweep.py`

```
runs/s32py_d4_recheck.sh [A:B]      # = zm_d4_sweep.py run --out runs/s32py_d4 --depth 24 --nproc-auto A:B (nice 10)
python3 search/zm_d4_sweep.py run --out runs/s32py_d4 --depth 24 --deepen 30 --nproc 8    # uncertified roots again
python3 search/zm_d4_sweep.py summary --out runs/s32py_d4      # -> runs/s32py_d4/SUMMARY.txt, exit 0 iff CLEAN
```

On first use it copies `search/zeromargin.py` into `OUT/checker/` and imports only that copy.  `manifest.json` holds
the certificate and checker sha256, the runner's sha256 and the settings; a resumed or `--deepen` run refuses any
mismatch.  Each worker re-hashes both files and re-checks D4 invariance exactly.  One task per root, no time cap;
finished roots are appended (fsync) to `roots.jsonl` with census, uncertified boxes, CPU, depth, checker sha and
certificate sha, and a restart skips them.  `--nproc-auto A:B` runs A roots at a time while the `runs/s32d4_v1`
sweep is running, and B after.  A root counts as certified if **any** of its records (each a complete, independent
run of the root, same checker and certificate, depth ≥ 24) has no uncertified box.  `SUMMARY.txt` lists all 7,200
roots (verdict, depth limit, number of records, census, CPU, checker and certificate sha256), then any uncertified
boxes and bad records.  It prints `D4 RECHECK CLEAN` only if every root is present and certified and all hashes match.

### 11.6 Files

`search/zeromargin.py` (fast path, `--d4`, `--ref`, `--selfcheck`, weighted symmetry checks), `search/zm_d4_sweep.py`;
`runs/s32py_d4/` (manifest, frozen checker, `roots.jsonl`, `SUMMARY.txt`), `runs/s32py_d4.log`,
`runs/s32py_d4_recheck.sh`; `runs/s32py_d4_aborted_d24/` (a first launch, stopped after 53 roots to add `--deepen`;
superseded); `runs/s32py_c2half.{sh,log}`; `runs/s32py_regress/` (`regress.py`, `compare.py`, `fsample.py`,
`sample.txt`, `frozen.jsonl`, `new.jsonl`, `germfix_new.jsonl`, `unc1_*.jsonl`, `compare.txt`).

**Open:** the review of §10's gluing (unchanged, and now not needed: this re-check alone covers the region), a Lean
replay, the `certificates/s32/` bundle.

## 12. Close-out (2026-09-26)

| certificate | zmcheck `--d4` | zeromargin.py `--d4` (§11) |
|---|---|---|
| `runs/s32-close_candidate.txt` (31.7135, sha `a0d2d38f…`) | 3,595 / 3,600 roots (§10) | **7,200 / 7,200, 2.77 CPU-h** |
| `runs/s32_shift_v1.txt` (31.6979, sha `7598d39a…`, `S32_SHIFT.md`) | **3,600 / 3,600 with `ZM_MIXPAIR=1`, 6.0 CPU-h** (`runs/s32d4_v1`); stock leaves 31 boxes at the `(1.5,2.5)` germs | **7,200 / 7,200, 2.42 CPU-h** (`runs/s32py_d4_v1`) |

`ZM_MIXPAIR` (committed, opt-in, enabled by any value other than empty or `0`): in `chain_order`, if the two largest
branch families share an orientation, promote the largest family of the other orientation to second place.  Branch
order only; verdict logic untouched; edge-tile census identical with and without it.

Lean (`notes/lean-s32.md`): `s32_eq_six_of_checker : S32CheckerCover → minSide 32 = 6`, standard axioms only; the
cover's D4 invariance and total `< 32` are kernel-checked from `Sqpack/S32Data.lean`.  `S32CheckerCover` is the
fundamental-region statement that the zeromargin.py run above establishes for the candidate.
