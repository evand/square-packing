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
