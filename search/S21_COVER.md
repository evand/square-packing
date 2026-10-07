# `COVER^closed(5)`: the cover side of `s(21) = 5` (task s21-cover; 2026-09-23)

Route: `s(21) = 5` from a weighted closed point cover of `[0,5]²` of total `< 21`, the `m = 5` rung of
the `s(13) = 4` machinery (`RUNG2.md`).  This note is the cover side only: the converged heuristic
LP, the honest cost of explicit covers, and where the weight sits.  No exact zero-margin check was
started.  Labels: **heuristic** (sampled-row float LP, float scan, local search; no bound in either
direction), **certified** (exact `Fraction`), **estimate** (extrapolation).  Nothing here is certified
except the `m = 4` comparators quoted from `RUNG2.md`.

Code: `search/closed4.py` (extended, backwards-compatible, §6), `search/rung2_close.py` (extended, §6),
`search/family_rows.py` (unchanged; already takes `--s`), `search/s21_cover_eval.py` (new: `hardscan`,
`family`, `regions`).  Cores 14–27, 2026-09-22 23:45 → 2026-09-23 04:47 (5 h wall, 14 cores).

## 0. Answer

> **Converged LP (heuristic): `COVER^closed(5) ≈ 20.72–20.74`.**  Five runs with different row sets
> (row pitch 0.02 / 0.01, 59 / 104 angles, with and without the item-3 and near-tile pools) agree
> to `±0.01`; the value drifts down by `≲ 0.002` per round, dual coverage `1.02–1.06` (the `m = 4`
> main run ended at `1.018`).  The 2026-08-30 run (`FAMILY.md` §3.1, `20.75–20.9`) was at the top
> of this band when stopped; `S21_OVERHEAD.md` §3's "`s5conv` runs at `20.9–21.0` at iteration 4" read
> the first rounds of `s5convA/B`, before column generation brought them down.  **So `21 − LP ≈ 0.27`, i.e. `1.3 %`.**
>
> **Honest cost (heuristic, total / min capture).**  With the `m = 4` protocol (`closed4.py stress`,
> as in `FAMILY.md` §2's `~12.46`): the closing loop plateaus at **`20.91–21.00`** (LP `20.78–20.80`,
> stress min `0.989–0.994`).  With the stricter protocol of `S21_OVERHEAD.md` §2 (`s21_cover_eval.py
> hardscan`, which finds the tilted/near-tile dips `stress` steps over): **`21.016`** for the best
> closing-loop cover (`s5convR3`, LP `20.801`, min `0.98978`) and `21.57` for the plain converged-LP cover.  At `m = 4` the same strict
> protocol gives `12.58` for `closed4_best.txt` and `12.70` for the LP weights of the shipped
> certificate (`12.339 / 0.9714`), against `12.956` certified.
>
> **Verdict: no demonstrated room below 21.**  The whole budget is a factor `21 / 20.73 = 1.013` over
> the converged LP, for validity *and* checker margin together.  At `m = 4` validity alone cost
> `1.3–2.9 %` on the strict protocol (and `closed4_best.txt` is exactly invalid by `5.8 %` on the item-3 family) (`0.5 %` after 67 lenient closing rounds) and the shipped
> checker margin was `×1.05`; here the best closing loop is at `×1.0138` over the converged LP
> (`21.016 / 20.73`) on the strict protocol after `1.8 h`, still above 21, and a stronger loop
> (`s5convR5`) re-opened a `3.6 %` hole it had closed earlier (§2).  `W < 21` needs the `m = 4` total
> overhead (`5 %`) cut by a factor `≈ 4`.  It is not
> excluded — the lenient plateau is below 21 and nothing here is a bound — but it is marginal at
> best, consistent with `S21_OVERHEAD.md` §5 (`Λ* = 1.009–1.019` for `nu_f(5) ∈ [20.5, 20.7]`).

## 1. Runs

All `taskset -c 14-20` or `21-27`, `--nproc 7`, `--s 5 --col-pitch 0.05 --no-literature` (the
Bentz/Nagamochi/Friedman columns are `s = 4` sets, `FAMILY.md` §3.1).  `seedrows` =
`python3 search/family_rows.py runs/s5conv_seedrows.txt --s 5 --pitch 0.004` (99,910 item-3 poses,
no `--tilted`: its named poses are `s = 4` poses).  `neartile` = `runs/s5conv_neartile.txt`: the
25 tile centres `± 0.012` (step 0.003) × `θ = 0.02..1.5°` (step 0.02), 151,875 poses.

| tag | what | rounds / wall | LP (last) | dual cov. | notes |
|---|---|---|---|---|---|
| `s5convB` | cold `linprog`, lattice columns, `seedrows` as hard rows, row pitch 0.02 | 12 / 52 min | 20.7458 | 1.044 | LP solve grew to 1441 s/round: stopped |
| `s5convA*` | same with `--warm-lp` | aborted | — | — | one-shot warm start slower than cold after column additions (bench below) |
| `s5convC` → `C2` | `--from-cert` B it11, `--warm-lp --lazy-seed seedrows`, pitch 0.02 | 10 + 8 / 66 min | 20.7337 → **20.7206** | 1.021 / 1.055 | |
| `s5convD` → `D2` → `D3` | same, **row pitch 0.01, 104 angles** (`--deg-step 0.5`) | 8 + 7 + 10 / 137 min | **20.7367** | 1.06–1.18 | `D3` = the "converged LP" cover, `runs/s5convD3_it9_last.txt` |
| `s5convR4` | from `R2`, pitch 0.02, lazy `seedrows` **+ `neartile`**, colgen | 11 / 46 min | 20.7324 | 1.044 | near-tile rows cost `≈ 0` LP |
| `s5convR2` | `rung2_close.py` (fixed columns = D2's 942 orbits / 7,460 atoms; stress-driven rows) | 14 / 50 min | 20.7808 | — | stress min 0.9927 → cost 20.934 (lenient) |
| `s5convR3` | `rung2_close.py --near-tile` from R2 | 33 / 111 min | 20.8011 | — | in-loop stress min 0.99424 → 20.922 (lenient); **the best cover**, `runs/s5convR3_final_last.txt` |
| `s5convR5` | `rung2_close.py --near-tile --fine-tilted --nworst 150 --nseed 150` from R3 r15 | 8 / 60 min | 20.8426 | — | in-loop stress min 0.99270 → 20.996; 560k rows by round 7 (polish-visited rows are not thinned) |

Convergence of the LP (heuristic; `closed4.py run` columns: LP / lattice min / violated lattice poses / dual coverage):

| run, round | LP | latmin | viol | dualcov |
|---|---|---|---|---|
| B it2 (item-3 rows in from the start) | 21.037 | 0.958 | 1413 | 1.43 |
| B it7 | 20.783 | 0.975 | 1711 | 1.058 |
| C it9 | 20.7337 | 0.985 | 191 | 1.021 |
| C2 it7 | 20.7206 | 0.992 | 129 | 1.055 |
| D3 it9 (fine rows) | 20.7367 | 0.991 | 618 | 1.177 |
| R4 it10 (+ near-tile pool) | 20.7324 | 0.996 | 104 | 1.044 |

Not converged in the strict sense (dual coverage never reached `1.00`, lattice `viol` never `0`), as
at `m = 4`.  The fine-row LP sits `+0.016` above the coarse one (at `m = 4` the fine-row refinement
added `+0.108`, `CLOSED4.md` §4, but that was over a *fixed* column set).  Reading: **`20.73 ± 0.02`**.

## 2. Honest cost (heuristic)

Three float checks (`s21_cover_eval.py`, calibrated on the `m = 4` files):

* `stress` = `closed4.py stress` (pitch 0.005, integer degrees + a few tiny angles + random), the
  protocol of every `m = 4` "honest" number in `CLOSED4.md` / `FAMILY.md`.  Inside the closing loop it
  runs at pitch 0.003 (`--near-tile` adds `0.02..3°`).
* `hardscan` = `S21_OVERHEAD.md` §2's protocol: pitch 0.004, `θ ∈ {0} ∪ [0.02°, 3°] step 0.02 ∪
  [3°, 45°] step 0.25`, polish of the 1,400 worst.  Calibration: `s13_closed_cover_4.txt` →
  `1.0200245 = 1.05 × 0.97145` (S21_OVERHEAD: `0.9715`); `closed4.py stress` on the same file says
  `1.0314`, i.e. **the lenient protocol overstates the minimum by `~1 %`**.
* `family` = the item-3 family at pitch 0.0005 (the lattice and both scans step over it).  Calibration:
  `closed4_best.txt` → `0.9420217`, the known exact violation (`RUNG2.md` §4.3).

| cover | total | `stress` min / cost | `hardscan` min / cost | `family` min | status |
|---|---|---|---|---|---|
| `m = 4` `closed4_best.txt` | 12.4175 | 0.99267 / 12.509 | 0.98695 / 12.582 | **0.94202** | heuristic; the family value is the exact `9420217/10⁷` of `RUNG2.md` §4.3 |
| `m = 4` r5 LP weights (= shipped / 1.05) | 12.3390 | 0.9823 / 12.562 | 0.97145 / 12.702 | 0.9838 | heuristic (shipped values ÷ 1.05) |
| `m = 4` shipped certificate | 12.9560 | 1.0314 | 1.0200 | 1.0330 | **certified valid** (`RUNG2.md`) |
| `m = 5` `s5convD3_it9` (converged LP) | 20.7372 | — | 0.96153 / **21.567** | 0.96584 | heuristic |
| `m = 5` `s5convR2` r13 | 20.7810 | 0.99267 / 20.934 | 0.96597 / 21.513 | 0.99128 | heuristic |
| `m = 5` `s5convR3` r15 | 20.7873 | 0.99182 / 20.959 | 0.98196 / 21.169 | 0.99971 | heuristic |
| `m = 5` `s5convR3` r29 | 20.8002 | 0.99259 / 20.955 (in-loop) | 0.98764 / 21.060 | 0.99843 | heuristic |
| `m = 5` **`s5convR3` final (r32)** | **20.8014** | 0.99115 / **20.987** | 0.98978 / **21.016** | 0.99947 | heuristic |
| `m = 5` `s5convR5` r5 | 20.8329 | 0.99149 / 21.012 (in-loop) | 0.98843 / 21.077 | 1.00000 | heuristic |
| `m = 5` `s5convR5` r7 | 20.8428 | 0.99270 / 20.996 (in-loop) | **0.96420** / 21.617 | 0.99817 | heuristic: the edge-tile `0.34°` hole is back |

Where the strict scan finds the deficits (`m = 5`): the **edge tile nudged by `0.15–0.7°`**
(`(0.502, 3.498, 0.30°)` for R2, `(3.494, 4.494, 0.72°)` for D3 — `S21_OVERHEAD.md`'s binding
wall-type pose, transplanted), then, once near-tile rows are in, the **edge tile tilted `4.5–5.3°`**
(`(0.544, 1.471, 5.29°)` for R3 r15, `(0.539, 1.474, 4.6°)` for R4) and interior tilted poses at
`25–45°` within `0.5 %` of it; for the R3 final cover the binding poses are all tilted,
`0.98978` at `(0.705, 1.484, 38.5°)` (wall band) and `0.99056` at `(0.588, 1.448, 11.25°)`.  Each is fixed at a small LP cost once it is a row (the near-tile pool
moved the LP by `≈ 0`; the closing loops by `+0.05–0.08`); the difficulty is finding them, not paying for
them.  The closing loop is whack-a-mole exactly as at `m = 4` (`FAMILY.md` §2), and not monotone:
`s5convR5` r5 had the near-tile holes closed (strict `21.077`), r7 has the edge tile at `(1.498, 4.485, 0.34°)`
back at `0.9642` (strict `21.62`) — the dip is narrow in `c_x` (`1.4982–1.4985` across the polished
hits) and the pitch-0.003 in-loop lattice steps over it; only the polish finds it.  Strict-protocol trajectory of R3:
r15 `21.169`, r29 `21.060`, r32 `21.016` (`−0.005`/round at the end).

## 3. Where the weight sits

`s21_cover_eval.py regions` (exact sums of the exported rational weights; a point on an interior grid line is split
equally between the two cells it bounds, as `s21_overhead.py regions`):

| cover | corner cell (mean) | edge cell (mean) | interior cell (mean) | total |
|---|---|---|---|---|
| `m = 4` `closed4_best.txt` | 0.4977 | 0.7775 | 1.0517 | 12.4175 |
| `m = 4` shipped certificate | 0.5288 | 0.8098 | 1.0907 | 12.9560 |
| `m = 5` `s5convD3_it9` (LP) | **0.4997** | **0.8048** | **1.0090** | 20.7372 |
| `m = 5` `s5convR3` final (closing loop) | 0.5013 | 0.8051 | **1.0150** | 20.8014 |

`s5convD3_it9`, per cell (rows = `y` cells, top first):

```
0.500  0.797  0.821  0.797  0.500
0.797  1.014  1.006  1.014  0.797
0.821  1.006  1.000  1.006  0.821
0.797  1.014  1.006  1.014  0.797
0.500  0.797  0.821  0.797  0.500
```

* **Corners: exactly `4 × 1.0000`** inside the four closed corner squares `[0,1]²` (all `m = 4` and
  `m = 5` covers): the axis-aligned corner square is tight and nothing else is spent there.  Split
  over the cells, `0.50` per corner cell at both `m`.
* **Wall bands** (`min(dx,dy) ≤ 1 < max`): `9.99` of `20.74` (48.2 %), vs `6.12` of `12.42` (49.3 %)
  at `m = 4`.  Edge cells *rise* from `0.78` to `0.80`.
* **Interior** (`min(dx,dy) > 1`): `6.75` (32.5 %); interior cells *fall* from `1.05` to `1.01`, the
  centre cell to exactly `1.000`.  The `m = 4` interior cells all touch the edge ring; at `m = 5` the
  one cell that doesn't is at the tiling value.
* By distance to the wall: `0<d<1` 44.7 %, `d = 1` 22.7 %, `1<d<2` 24.4 %, `d = 2` 6.7 %, `d > 2` 1.5 %.
* 82 % on the integer grid lines (`x = 1, 4`: 2.166 each, `x = 2, 3`: 2.132), 4 % within 0.05 of one,
  13 % elsewhere (`m = 4`: 91 / 4 / 5 %).  0 on the walls and the half-integer lines.
* The cellular extrapolation from `closed4_best` (`S21_OVERHEAD.md` §3: `20.79`) overshoots by
  `0.05`, all of it the interior-cell drop.  Cell model `4a + 12b + 9c` with the `m = 5` means:
  per extra ring the LP pays `≈ 0.80` per new edge cell and `≈ 1.00` per new interior cell, so
  `LP(m) − (m² − 4) = 4a + 4(m−2)(b − 1) + (m−2)²(c − 1) ≈ 2.0 − 0.78(m−2) + 0.009(m−2)²` (estimate, `a, b, c` frozen at their `m = 5` values):
  `−0.27` at `m = 5`, `−1.0` at `m = 6`.  The room below `m² − 4` grows linearly in `m` *if* the
  overhead factor stays bounded; the overhead is `O(m²)` (`S21_OVERHEAD.md` §5), so the balance at
  `m = 6` is `≈ 1.0` of room against `≈ 0.05 × 31 = 1.55` at the `m = 4` factor (`0.62` at a `2 %` factor).

## 4. Room below 21

| quantity | value | status |
|---|---|---|
| converged cover LP | `20.73 ± 0.02` | heuristic |
| room `21 − LP` | `0.27` (`1.3 %`) | heuristic |
| best honest cost, lenient (`stress`, `m = 4`-comparable) | `20.92` (in-loop) / `20.99` (`closed4.py stress`), `s5convR3` final | heuristic |
| best honest cost, strict (`hardscan`) | `21.016` (`s5convR3` final) | heuristic; the true minimum capture is at most the scanned one, so the true cost is `≥` this for this point set |
| `m = 4`: LP → certified | `12.339 → 12.956` (`×1.05`) | certified |
| `m = 5` at the `m = 4` factor | `20.73 × 1.05 = 21.77` | estimate |
| break-even factor | `21 / 20.73 = 1.013` | heuristic |

## 5. Performance notes (why the old run stalled)

* The LP solve dominates (`>95 %` of the main process), not the scans.  With the 50k item-3 rows as
  hard rows a cold `linprog` reached 1441 s per round by round 11.
* `--warm-lp` (persistent `highspy` model).  Bench (`s = 5`, same model each round, cold / one-shot warm /
  two-phase warm / IPM): it3 `29 / 105 / 28 / 65 s`, it4 `73 / 183 / 39 / 90 s`.  One-shot warm is *slower*
  after columns are added (the basis is neither primal nor dual feasible); two-phase (primal simplex after
  `addCols`, then dual simplex after `addRows`) is `~2×` faster than cold.  That is what `solve_warm` does.
* `--lazy-seed` keeps the 100k-pose item-3 pool out of the LP and adds only violated poses each round
  (`captured` on the pool, parallel, seconds); with `--prune-at` the LP stays at 20–35k rows.
* Column growth (`+250` orbits/round) is what still makes rounds slow (100 s → 500 s over 10 rounds).
  `--col-prune-at` drops unused orbits with reduced cost `> 2 %`, but only ~10 % qualify, and every prune
  costs a cold solve, so it is now episodic (next prune after 50 % growth).  The restarts `C → C2`,
  `D → D2 → D3` (`--from-cert` the checkpoint: support-only columns) were the effective column prune.
* `rung2_close.py` rounds: 180–260 s (`--near-tile`), 350–550 s (`--fine-tilted --nseed 150`, which
  adds 40–120k polish-visited rows per round).

## 6. Code changes (all backwards-compatible; defaults reproduce the old behaviour)

`search/closed4.py`:
* `Model.solve_warm()` + `run --warm-lp`: two-phase warm-started `highspy` re-solve; falls back to `solve()`.
* `run --lazy-seed`: `--seed-rows` as a separation pool (`pool_violations`).
* `Model.prune_cols()` + `run --col-prune-at N --col-prune-rc 0.02` (default off).
* **Bug fix in the row prune** (`--prune-at`, default 90000, never reached at `s = 4`): `keep` had one
  entry per row *at solve time*, so `Model.prune` silently dropped every row added in that round
  (`zip` truncation) and `keep[-prune_keep:]` pointed at old rows.  Now padded to the current row
  count; tight rows (`Ax ≤ 1 + 1e-7`) are kept as well as dual-support rows.
* `stress(..., near_tile=False, fine_tilted=False, nworst=10)`: optional extra angles, more worst poses returned.
* log line gains `lp=<s>` (LP solve time); JSON history gains `lp_s`.
* No `s = 4` hard-coding on the paths used: `sanity` / `nagamochi` / the literature columns remain
  `s = 4`-only (not used at `s = 5`, `--no-literature`).

`search/rung2_close.py`: `--warm-lp`, `--near-tile`, `--fine-tilted`, `--nworst`; `--cap` documented
(default `13.0` is the `s = 4` budget: pass `--cap 21.5` at `s = 5`, else it aborts at round 0).

`search/family_rows.py`: unchanged (`--s 5` works; `--tilted` names `s = 4` poses, not used).

## 7. Files

| file | what |
|---|---|
| `runs/s5conv_seedrows.txt`, `runs/s5conv_neartile.txt` | the item-3 and near-tile pose pools |
| `runs/closed4_s5conv{B,C,C2,D,D2,D3,R4}.log` / `.json`, `runs/s5conv*.nohup.log` | the LP runs (§1) |
| `runs/s5convD3_it9_last.txt` | the converged-LP cover, 9,548 points, total 20.737213 |
| `runs/s5convC_stage1_last.txt`, `runs/s5convC2_stage2_last.txt`, `runs/s5convD2_stage2_last.txt`, `runs/s5convR4_stop_last.txt` | other checkpoints |
| `runs/closed4_s5convR{2,3,5}.log` / `.json`, `runs/s5convR2_last_final.txt`, `runs/s5convR3_r15_last.txt` | closing loops |
| `runs/s5convR3_final_last.txt` (= `closed4_s5convR3_last.txt`) | the best closing-loop cover, 4,800 points, total 20.801370; `runs/s5convR3_r29_last.txt`, `runs/s5convR5_r{5,7}_last.txt` the others in §2 |
| `runs/closed4_s5convA_aborted.log`, `runs/closed4_s5convA_v2aborted.*`, `runs/closed4_s5convR_v1.*`, `runs/closed4_s5conv_stressB10.log` | aborted / diagnostic |

All `_last.txt` files are certificate format (`s = 5`, `D = 1000`, `W = 10⁷`, weights rounded up),
D4-symmetric, **not** valid covers (min capture `< 1`); scale by `1/min` (`scale_cover.py`) before any
exact check.

## 8. Reproduce

```
python3 search/family_rows.py runs/s5conv_seedrows.txt --s 5 --pitch 0.004
# converged LP (the D chain; C chain = same without --pitch/--deg-step/--block)
python3 search/closed4.py run --s 5 --tag s5convB --nproc 7 --col-pitch 0.05 --no-literature \
        --seed-rows runs/s5conv_seedrows.txt --prune-at 300000 --prune-keep 100000 --time 14400   # stopped at it11
python3 search/closed4.py run --s 5 --tag s5convD --from-cert runs/s5convB_it11.txt --pitch 0.01 --deg-step 0.5 \
        --block 0.05 --warm-lp --lazy-seed --seed-rows runs/s5conv_seedrows.txt --cg-want 300 \
        --prune-at 60000 --prune-keep 20000 --nproc 7            # then D2, D3: --from-cert the previous _last, + --col-prune-at
# closing loops
python3 search/rung2_close.py runs/s5convR_r1_last.txt s5convR2 --warm-lp --cap 21.5 --nproc 7 --nrand 100
python3 search/rung2_close.py runs/s5convR2_last_final.txt s5convR3 --warm-lp --near-tile --cap 21.5 --nproc 7 --nrand 60
python3 search/rung2_close.py runs/s5convR3_r15_last.txt s5convR5 --warm-lp --near-tile --fine-tilted \
        --nworst 150 --nseed 150 --cap 21.5 --nproc 7 --nrand 40
# evaluation
python3 search/s21_cover_eval.py hardscan runs/s5convD3_it9_last.txt       # 0.9615340 -> 21.567
python3 search/s21_cover_eval.py family   runs/s5convD3_it9_last.txt       # 0.9658376
python3 search/s21_cover_eval.py regions  runs/s5convD3_it9_last.txt
python3 search/s21_cover_eval.py hardscan certificates/rung2/s13_closed_cover_4.txt   # calibration 1.0200245
```
