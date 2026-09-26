# `COVER^closed(7)`: the cover side of `s(45) = 7` (task s45-cover; 2026-09-24/25)

Route: `s(45) = 7` from a weighted closed point cover of `[0,7]²` of total `< 45`, the `m = 7` rung of the
`s(13) = 4` machinery (`RUNG2.md`), measured with the `m = 5` recipe of `S21_COVER.md`.  This note is the cover
side only: the converged heuristic LP, the honest cost of explicit covers, and where the weight sits.  No exact
zero-margin check was started.  Labels as in `S21_COVER.md`: **heuristic** (sampled-row float LP, float scan,
local search; no bound in either direction), **certified** (exact `Fraction`), **estimate** (extrapolation).
Nothing here is certified.

Code: `search/closed4.py` (new `--lp-ipm`, §5), `search/s45_prep.py` (new: `transplant`, `neartile`),
`search/family_rows.py`, `search/s21_cover_eval.py`, `search/rung2_close.py` (all unchanged; they already take
any `s`; `rung2_close.py --cap` was passed explicitly as `45.5`).  Cores 8–14 + 24–30, 2026-09-24 17:45 →
2026-09-25 02:50 (≈ 9 h wall, 7 physical cores; peak RSS ≈ 6 GB per run).

## 0. Answer

> **Converged LP (heuristic): `COVER^closed(7) ≈ 43.74–43.83`.**  The coarse-row chain (pitch 0.02, 59 angles,
> item-3 + near-tile pools) ends at **`43.738`** (dual coverage `1.02`, drifting `−0.003`/round); the fine-row
> chain (pitch 0.01, 104 angles) at `43.827` (`D3`, flat at `43.82–43.83` over 5 rounds, §1.1).  At `m = 5` the
> fine rows sat `+0.016` over the coarse ones; here `+0.09`.  Reading: **`LP ≈ 43.78 ± 0.05`**.
> **Room `45 − LP ≈ 1.22` (`2.7 %`), break-even factor `45 / LP = 1.027–1.029`** (`m = 5`: `0.27`, `1.3 %`, `1.013`).
>
> **Honest cost (heuristic, strict protocol `s21_cover_eval.py hardscan` + `family`): best `44.698`**
> (`s7convR1` round 10, LP weights total `43.9135`, strict min `0.98245`, item-3 min `0.99197`); the closing
> loop then plateaus (`45.37 → 44.78 → 44.71 → 44.70 → 44.80 → 44.82` over rounds 5, 8, 9, 10, 11, final): the
> strict cost of the closing-loop covers is **`44.70–44.82`**, a factor `1.022–1.025` over the coarse LP.  The plain
> converged-LP cover (`C3`) costs `45.27` strict.  Lenient (`closed4.stress`): `44.43` (r10, in-loop), `44.58` (final).
>
> **Verdict: `W < 45` is plausible with the known overhead components — the first rung where it is.**  The best
> strict cover is already `0.30` (`0.67 %`) under 45 before any checker margin, and every closing-loop cover from
> round 8 on is `≥ 0.18` under.  The exact checker needs
> `λ_c ≤ 45 / 44.698 = 1.0068` over the strict-scan minimum; at `m = 4` it certified the whole container at
> `λ_c = 1.0025` over the true minimum (`M4_MARGIN.md`), which here would give `W ≈ 44.81` (estimate).  The risk
> is validity, not margin: the strict scan must be close to the true minimum (at `m = 4` it was:
> `1.0200245` scanned vs `1.0201154` exact), i.e. the exact checker must not find a dip more than `0.42 %` below
> the strict minimum (`44.698 × 1.0025 × (1 + δ) < 45` ⇔ `δ < 0.42 %`).  At `m = 5` the strict scan missed dips of
> `2.3 %` (`s5convR5` r7) between rounds — so the next task should run more closing rounds (the loop was not
> converged) and then the exact checker on a cover scaled to `1.0025–1.005` over its strict minimum.

## 1. Runs

All `taskset -c 8-14,24-30`, `--nproc 7`, `--s 7 --col-pitch 0.05 --no-literature`, lazy pools `seedrows` =
`runs/s7_seedrows.txt` (`family_rows.py --s 7 --pitch 0.004`, 209,874 item-3 poses) and `neartile` =
`runs/s7_neartile.txt` (`s45_prep.py neartile 7`: tile centres `± 0.012` step 0.003 × `θ = 0.02..1.5°` step
0.02, one cell per C4 orbit — C4 fixes `θ` — 13 cells, 78,975 poses).

Warm start (no cold start was run): **transplant** of the `m = 5` covers (`s45_prep.py transplant`): the centre
column/row of unit cells `[2,3]` is replicated three times, so the corner and wall-band cells are carried over
unchanged and the centre cell tiles the new `5 × 5` interior.  `runs/s7_tpD3.txt` (from `s5convD3_it9`, 24,188
points, total `43.3575` = the cell-model estimate) and `runs/s7_tpR3.txt` (from `s5convR3_final`, 11,724 points,
`43.4030`).  As covers they are poor (`0.02` lattice min `0.753`: the seams); they serve as column sets.

| tag | what | rounds / wall | LP (last) | dual cov. | notes |
|---|---|---|---|---|---|
| `s7convA` | `--from-cert s7_tpD3`, `--warm-lp` (simplex), prune 120k/40k | 4 / 87 min | 43.7084 | 1.56 | rows not converged (latmin 0.94); LP solve 3525 s at it3: stopped |
| `s7convB` | same from `s7_tpR3` | 6 / 94 min | 43.8976 | 1.20 | LP solve 2951 s at it5: stopped |
| `s7convC` | `--from-cert A it3`, **`--lp-ipm`**, pitch 0.02, prune 100k/30k | 13 / 62 min | 43.8047 | 1.050 | IPM: 34–600 s per solve |
| `s7convD` | `--from-cert B it5`, `--lp-ipm`, **pitch 0.01, 104 angles**, `--block 0.05` | 10 / 47 min | 43.9204 | 1.081 | |
| `s7convC2` | `--from-cert C it12`, prune **45k/12k**, `--col-prune-at 2500` | 12 / 99 min | 43.7617 | 1.03 (1.61 last) | |
| `s7convD2` | `--from-cert D it9`, fine rows, prune 60k/15k | 10 / 78 min | 43.8350 | 1.043 | `−0.01`/round; stopped for CPU |
| `s7convC3` | `--from-cert C2 it11`, prune 45k/12k | 11 / 137 min | **43.7382** | 1.020 | **the converged coarse LP**, `runs/s7convC3_it10_last.txt` |
| `s7convD3` | `--from-cert D2 it9`, fine rows, prune 60k/15k | 5 / 98 min | 43.8270 | 1.21 | the fine-row LP, `runs/s7convD3_it4_last.txt` |
| `s7convR1` | `rung2_close.py --near-tile --cap 45.5 --nrand 60` from `C2 it8` (2,638 orbits / 20,920 atoms), `CLOSED4_LP=ipm` | 12 / 4.7 h | 43.9167 (final) | — | 17–30 min/round; **the best cover**, `runs/s7convR1_r10_last.txt` |

Convergence of the coarse chain (heuristic; LP / lattice min / violated lattice poses / dual coverage):

| run, round | LP | latmin | viol | dualcov |
|---|---|---|---|---|
| C it2 | 44.1159 | 0.987 | 2427 | 1.457 |
| C it7 | 43.9054 | 0.983 | 2794 | 1.277 |
| C it12 | 43.8047 | 0.993 | 3316 | 1.050 |
| C2 it5 | 43.7944 | 0.994 | 3059 | 1.041 |
| C2 it11 | 43.7617 | 0.985 | 3711 | 1.608 |
| C3 it5 | 43.7439 | 0.994 | 4406 | 1.026 |
| C3 it10 | **43.7382** | 0.990 | 4496 | 1.020 |
| D2 it9 (fine rows) | 43.8350 | 0.993 | 6731 | 1.043 |

Not converged in the strict sense (dual coverage `1.02`, lattice `viol` never 0), exactly as at `m = 4, 5`.
The value came *down* from the transplant runs' `44.1` (early rounds, support-only columns) by column generation,
and sits **`+0.38` above the transplant / cell-model estimate `43.36`** (§3: the edge cells are heavier at `m = 7`).

### 1.1 Fine rows

`s7convD3` (from `D2 it9`, pitch 0.01, 104 angles): it0 `43.8168` (80k rows, latmin 0.79, LP solve 1907 s), it1
`43.8223`, it2 `43.8289`, it3 `43.8286`, it4 `43.8270` (latmin 0.983, dual coverage 1.21); stopped at its 5400 s
limit.  The D chain came from the sparser `s7_tpR3` column set and had fewer column-generation rounds (`D → D2 → D3`:
31 LP rounds from `B` vs 40 for `A → C → C2 → C3`), so part of the `+0.09` may be columns, not rows.  Cell means `0.5058 / 0.8308 /
1.0075` (§3): the difference is spread over the interior (`+0.0024` per interior cell, `25 × 0.0024 = 0.06`).

## 2. Honest cost (heuristic)

Protocols as in `S21_COVER.md` §2 (`stress` lenient; `hardscan` = pitch 0.004, `θ ∈ {0} ∪ [0.02°,3°] step 0.02 ∪
[3°,45°] step 0.25`, polish of the 1,400 worst; `family` = item-3 family at pitch 0.0005).  `hardscan` at `s = 7`:
300–440 s on 7 cores (with other jobs running).

| cover | total | in-loop `stress` min / cost | `hardscan` min / cost | `family` min | status |
|---|---|---|---|---|---|
| `m = 4` `closed4_best.txt` | 12.4175 | 0.99267 / 12.509 | 0.98695 / 12.582 | 0.94202 | from `S21_COVER.md` |
| `m = 4` shipped certificate | 12.9560 | 1.0314 | 1.0200 | 1.0330 | certified valid (`RUNG2.md`) |
| `m = 5` `s5convD3_it9` (converged LP) | 20.7372 | — | 0.96153 / 21.567 | 0.96584 | from `S21_COVER.md` |
| `m = 5` `s5convR3` final (best) | 20.8014 | 0.99115 / 20.987 | 0.98978 / **21.016** | 0.99947 | from `S21_COVER.md` |
| `m = 7` `s7convC` it12 | 43.8055 | — | 0.94698 / 46.258 | — | heuristic |
| `m = 7` `s7convC3` it10 (converged LP) | 43.7397 | — | 0.96626 / **45.267** | 0.98313 | heuristic |
| `m = 7` `s7convR1` r5 | 43.8994 | 0.97417 / 45.063 | 0.96762 / 45.369 | 0.98924 | heuristic |
| `m = 7` `s7convR1` r8 | 43.9069 | 0.98483 / 44.584 | 0.98040 / 44.785 | 0.99209 | heuristic |
| `m = 7` `s7convR1` r9 | 43.9110 | 0.98715 / 44.482 | 0.98213 / 44.710 | 0.99060 | heuristic |
| `m = 7` **`s7convR1` r10** | **43.9135** | 0.98830 / 44.433 | 0.98245 / **44.698** | 0.99197 | heuristic; **best** |
| `m = 7` `s7convR1` r11 | 43.9163 | 0.98880 / 44.414 | 0.98032 / 44.798 | 0.98985 | heuristic |
| `m = 7` `s7convR1` final (all rows, = r12) | 43.9175 | 0.98504 / 44.584 (`closed4.stress` near-tile) | 0.97986 / 44.820 | 0.98979 | heuristic |

Where the strict scan finds the deficits (`m = 7`): the same families as at `m = 5`.  For the plain LP cover the
**wall tile nudged by `0.26–0.28°`** (`(5.498, 6.498, 0.26°)`, `0.9663`; `S21_OVERHEAD.md`'s binding pose) with
**tilted wall-band / interior poses at `32–35°`** within `0.1 %` of it (`(5.695, 1.393, 35.1°)`, `0.9672`).
Once the closing loop has those as rows the binding poses are all tilted (`25–41°`, e.g. `(2.970, 1.530, 25.7°)`
for r8, `(5.620, 5.675, 32.9°)` for r9), and at r10 the **centre tile nudged by `0.11°`** (`(2.500, 2.500,
0.113°)`, `0.98245`) — an interior near-tile dip, which at `m = 5` never bound; at r11 / final the binding poses are
interior tilted at `38–40°` (`(2.281, 2.335, 39.6°)`, `(1.530, 5.787, 38.4°)`).  As at `m = 5`, the closing loop is
whack-a-mole: each dip costs `≈ 0` LP once it is a row (LP `43.84 → 43.92` over 12 rounds, `+0.18 %`), the strict
minimum oscillates in `0.980–0.982` from round 8 on, and the difficulty is finding the dips, not paying for them.

## 3. Where the weight sits

`s21_cover_eval.py regions` (exact sums of the exported weights; grid-line points split equally):

| cover | corner cell (mean) | edge cell (mean) | interior cell (mean) | total |
|---|---|---|---|---|
| `m = 5` `s5convD3_it9` (LP) | 0.4997 | 0.8048 | 1.0090 | 20.7372 |
| `m = 7` transplant `s7_tpD3` (= cell model) | 0.4997 | 0.8115 | 1.0052 | 43.3575 |
| `m = 7` `s7convC2` it11 (LP) | 0.5038 | 0.8295 | 1.0063 | 43.7627 |
| `m = 7` **`s7convC3` it10 (LP)** | **0.5024** | **0.8302** | **1.0051** | 43.7397 |
| `m = 7` `s7convD2` it9 (fine-row LP) | 0.5053 | 0.8292 | 1.0092 | 43.8360 |
| `m = 7` `s7convD3` it4 (fine-row LP) | 0.5058 | 0.8308 | 1.0075 | 43.8281 |

`s7convC3_it10`, per cell (rows = `y` cells, top first):

```
0.502  0.814  0.841  0.841  0.841  0.814  0.502
0.814  1.019  1.006  1.001  1.006  1.019  0.814
0.841  1.006  1.000  1.000  1.000  1.006  0.841
0.841  1.001  1.000  1.000  1.000  1.001  0.841
0.841  1.006  1.000  1.000  1.000  1.006  0.841
0.814  1.019  1.006  1.001  1.006  1.019  0.814
0.502  0.814  0.841  0.841  0.841  0.814  0.502
```

* **Corners** `4 × 1.0000` inside the four closed corner squares, as at `m = 4, 5`; `0.50` per corner cell.
* **Interior cells are at the tiling value**: the `3 × 3` block not touching the edge ring is `1.000` to three
  decimals, the ring of interior cells next to the wall band `1.001–1.019`.  Interior mean `1.005` (`m = 5`: `1.009`).
* **Edge cells are heavier than the cell model assumed:** `0.814` next to a corner (`m = 5`: `0.797`), **`0.841`**
  in the middle of a wall (`m = 5`: `0.821`, the one middle edge cell).  That is the whole miss of the cell model:
  `20 × (0.830 − 0.8115) = 0.37` of the `+0.38`.  So the per-cell model of `S21_COVER.md` §3 is right in form but its
  edge constant is `b ≈ 0.83`, not `0.80`: `LP(m) − (m² − 4) ≈ 2.0 − 0.68(m−2) + 0.005(m−2)²` (estimate), giving
  `−1.26` at `m = 7` (measured `−1.26`) instead of `−1.64`.  The room grows by `≈ 0.68` per rung, not `0.78`.
* Wall bands (`min(dx,dy) ≤ 1 < max`): `18.05` (41.3 %); interior (`min(dx,dy) > 1`): `21.69` (49.6 %).
  By distance to the wall: `0<d<1` 34.7 %, `d = 1` 15.8 %, `1<d<2` 24.8 %, `d = 2` 8.5 %, `d > 2` 16.3 %.
* 72.8 % on the integer grid lines (`x = 1, 6`: 2.712; `x = 2, 5`: 2.632; `x = 3, 4`: 2.631), 5.5 % within 0.05
  of one, 21.7 % elsewhere (`m = 5`: 82 / 4 / 13 %) — heuristic: the IPM solutions are not vertices, so some weight
  is spread thinly over near-duplicate columns (the exported point count is 26,620 vs `m = 5`'s 9,548).

## 4. Room below 45

| quantity | value | status |
|---|---|---|
| converged cover LP | `43.74` (coarse rows) – `43.83` (fine rows); reading `43.78 ± 0.05` | heuristic |
| room `45 − LP` | `1.22 ± 0.05` (`2.7 %`) | heuristic |
| break-even factor `45 / LP` | `1.027–1.029` (`m = 5`: `1.013`) | heuristic |
| best honest cost, lenient (in-loop `stress`) | `44.41–44.43` (`s7convR1` r10, r11) | heuristic |
| best honest cost, strict (`hardscan` + `family`) | **`44.698`** (`s7convR1` r10; min `0.98245`) | heuristic; the true minimum is at most the scanned one, so the true cost of this point set is `≥ 44.698` |
| strict cost of the closing-loop plateau | `44.70–44.82` (rounds 8–12) | heuristic |
| validity factor reached, strict cost / LP | `44.698 / 43.74 = 1.022` (`m = 5`: `1.0138`, `m = 4`: `1.029`) | heuristic |
| checker factor left | `45 / 44.698 = 1.0068` | heuristic |
| `m = 4` checker factor needed | `1.0025` over the true minimum (`M4_MARGIN.md`) | certified partial (`m = 4`) |
| estimated certifiable `W` | `44.698 × 1.0025 ≈ 44.81` | estimate (assumes the strict scan is within `0.4 %` of the true minimum) |

## 5. Performance notes

* **The LP is the bottleneck, and at `s = 7` simplex is the wrong algorithm.**  Warm simplex (`--warm-lp`)
  took 3525 s for one re-solve at 90k rows × 3.4k orbits (`s7convA` it3), and got *slower* each round.  Bench on a
  166k × 906 model (25M nonzeros): dual simplex 1190 s, primal simplex on the dual 1018 s, **HiGHS IPM without
  crossover 127 s** (value agrees to `1e-11`), IPM on the dual 160 s; HiGHS `hipo` (4 threads) did not finish in
  10 min.  Hence `closed4.py --lp-ipm` / `CLOSED4_LP=ipm` (cold IPM each round; §6).
* With a cold solve every round the row prune is free, so prune hard: `--prune-at 45000 --prune-keep 12000` keeps
  the LP at 25–38k rows.  IPM time then scales with the column count (70 s at 2.1k orbits → 1000 s at 4.8k), so the
  `--from-cert` restarts (support-only columns) are still the effective column prune, every ~10 rounds.
* Fine rows (pitch 0.01, 104 angles) cost `≈ 2×` per round in the LP and `4×` in the scans.
* `rung2_close.py` at `s = 7`: 17–30 min per round (in-loop near-tile stress + IPM on 73k → 130k rows, sharing the
  cores with an LP chain and the strict scans); `hardscan` 5–7 min, `family` 2 min.

## 6. Code changes (backwards-compatible; defaults reproduce the old behaviour)

`search/closed4.py`:
* `Model.solve_ipm()`: cold HiGHS IPM, `run_crossover = off`; primal values and duals `< 1e-9` set to 0.
* module flag `LPSOLVER` (default `'simplex'`, or env `CLOSED4_LP=ipm`); `Model.solve()` dispatches to
  `solve_ipm` when it is `'ipm'`; `run --lp-ipm` sets it and then ignores `--warm-lp`.  `rung2_close.py`
  (unchanged) picks it up through the environment variable when run without `--warm-lp`.

`search/s45_prep.py` (new): `transplant IN M OUT` (grow an odd-`m0` D4 cover to `[0,M]²` by replicating the centre
cell column/row; identity for `M = m0`, checked) and `neartile M OUT` (the near-tile pool, one cell per C4 orbit).

Checked for hard-coded `s`: `s21_cover_eval.py` (all modes take `s` from the file), `family_rows.py --s`,
`scale_cover.py` (exact, `s`-free), `rung2_close.py --cap` (passed as `45.5`; the sibling task has since made the
default `s² − 3.5`) — nothing else needed changing.

## 7. Files (`runs/`, gitignored)

| file | what |
|---|---|
| `s7_seedrows.txt`, `s7_neartile.txt` | the item-3 and near-tile pose pools |
| `s7_tpD3.txt`, `s7_tpR3.txt` | the transplanted `m = 5` covers (column sets) |
| `closed4_s7conv{A,B,C,D,C2,D2,C3,D3}.log` / `.json`, `s7conv*.nohup.log` | the LP runs (§1) |
| `s7convA_stop_last.txt`, `s7convB_stop_last.txt`, `s7convC_stop_last.txt`, `s7convD_stop_last.txt`, `s7convC2_it{8,11}_last.txt`, `s7convD2_it9_last.txt` | restart checkpoints |
| **`s7convC3_it10_last.txt`** | the converged coarse-LP cover, 26,620 points, total 43.739676 |
| `closed4_s7convR1.log` / `.json`, `s7convR1_r{5,7,8,9,10,11}_last.txt`, `s7convR1_final_best.txt` (= `closed4_s7convR1_best.txt`) | the closing loop and its checkpoints |
| `s7convD3_it4_last.txt` | the fine-row LP cover, total 43.828 |
| **`s7convR1_r10_last.txt`** | the best cover, total 43.913506, strict min 0.98245 → 44.698 |
| `s7_eval_*.log`, `s7_hardscan_Cstop.log`, `s7_lpbench.log` | strict evaluations (§2), the LP bench (§5) |

All `_last.txt` files are certificate format (`s = 7`, `D = 1000`, `W = 10⁷`, weights rounded up),
D4-symmetric, **not** valid covers (min capture `< 1`); scale by `1/min` (`scale_cover.py`) before any exact check.

## 8. Reproduce

```
python3 search/family_rows.py runs/s7_seedrows.txt --s 7 --pitch 0.004
python3 search/s45_prep.py neartile 7 runs/s7_neartile.txt
python3 search/s45_prep.py transplant runs/s5convD3_it9_last.txt 7 runs/s7_tpD3.txt
POOL="--lazy-seed --seed-rows runs/s7_seedrows.txt,runs/s7_neartile.txt --no-literature --nproc 7 --cg-want 300"
python3 search/closed4.py run --s 7 --tag s7convA --from-cert runs/s7_tpD3.txt --warm-lp $POOL --prune-at 120000 --prune-keep 40000   # stopped at it3
python3 search/closed4.py run --s 7 --tag s7convC  --from-cert runs/s7convA_stop_last.txt  --lp-ipm $POOL --prune-at 100000 --prune-keep 30000
python3 search/closed4.py run --s 7 --tag s7convC2 --from-cert runs/s7convC_stop_last.txt  --lp-ipm $POOL --prune-at 45000 --prune-keep 12000 --col-prune-at 2500
python3 search/closed4.py run --s 7 --tag s7convC3 --from-cert runs/s7convC2_it11_last.txt --lp-ipm $POOL --prune-at 45000 --prune-keep 12000 --time 7200
#   fine rows (D chain, from s7_tpR3 via s7convB): add --pitch 0.01 --deg-step 0.5 --block 0.05, prune 60000/15000
# closing loop
CLOSED4_LP=ipm python3 search/rung2_close.py runs/s7convC2_it8_last.txt s7convR1 --near-tile --cap 45.5 --nproc 7 --nrand 60 --time 14400
# evaluation
python3 search/s21_cover_eval.py hardscan runs/s7convR1_r10_last.txt     # 0.9824455 -> 44.698
python3 search/s21_cover_eval.py family   runs/s7convR1_r10_last.txt     # 0.9919715
python3 search/s21_cover_eval.py regions  runs/s7convC3_it10_last.txt
```

## 9. Closing to a stable strict minimum (task s45-close, 2026-09-25; written up 2026-09-26 from the logs)

Everything here is **heuristic**.  The agent that ran it stopped after the confirmation scan without writing this
section; the numbers below are read from `runs/s45closeA.log`, `runs/s45closeA_eval_r*.log` and
`runs/s45close_confirm_r6.log`.  Code: `search/close_strict.py` (new; the `m = 7` sibling of `close_shifted.py`:
every round measures the LP weights with the strict protocol itself at pitch 0.002 and keeps every dip as a
permanent row; see its docstring).

**Starting point.**  `hardscan` at pitch 0.002 on `s7convR1_r10` (the §0 best cover): polished min `0.98055` at
`(5.800, 2.470, 37.8°)`, cost **`44.785`** (at pitch 0.004 it was `44.698`: a `0.2 %` deeper dip at half the pitch).

**Loop** (`s45closeA`, from `s7convR1_r10`, columns + `s7convC2_it8`, 8 rounds, 3.0 h):

| round | total | strict min p = 0.002 | cost | `hardscan` p = 0.004 min / cost | where |
|---|---|---|---|---|---|
| 0 | 43.9054 | 0.85609 | 51.29 | 0.85611 / 51.28 | `(2.514, 0.554, 6.6°)` wall band |
| 1 | 43.9233 | 0.91502 | 48.00 | 0.91504 / 48.00 | `(2.481, 6.439, 7.4°)` |
| 2 | 43.9353 | 0.97666 | 44.99 | 0.97928 / 44.86 | edge tile `0.28°` |
| 3 | 43.9520 | 0.97862 | 44.91 | 0.98013 / 44.84 | wall band `15°` |
| 4 | 43.9600 | 0.96452 | 45.58 | 0.96352 / 45.62 | edge tile `0.19°` |
| 5 | 43.9622 | 0.99427 | 44.22 | 0.99520 / 44.17 | interior `29–37°` |
| **6** | 43.9671 | **0.99456** | **44.208** | 0.99506 / 44.19 | interior `35°` |
| 7 | 43.9698 | 0.99418 | 44.227 | 0.99389 / 44.24 | interior `40–45°` |

Declared STABLE after round 7 (`0.99427, 0.99456, 0.99418`).

**Confirmation scan on r6** (`close_strict.py confirm`, pitch 0.001, 150 dip boxes ±0.1, whole container at the 318
interleaved angles, polish of the 1,400 worst): dip boxes `0.99279`; full container `0.98751` at `(0.596, 5.483,
12.4°)`; **polished `0.98170` at `(1.502, 0.502, 0.26°)`, the edge tile turned `0.26°`**.  Cost **`44.787`**.
The confirmation minimum sits **`1.29 %`** below the loop's stable strict minimum (`m = 6`: `0.58 %`), so the loop's
stopping rule was fooled: its lattice never resolved the thin edge-tile cells.

**Verdict: no-go (heuristic).**  A candidate at `1.004 ×` over the confirmed minimum would total `≈ 44.966`, `0.08 %`
under 45: no room for further unfound dips, and (per `S32_EXACT.md`) the exact checker's real cost at `m ≥ 6` is the
interior tile germs, not margin.  No candidate file was written.  Revisit only after `s(32)` closes, with the
near-tile family handled exactly rather than by scans.

```
bash runs/s45closeA.sh          # close_strict.py loop ... --hs-pitch 0.002 --eval-004 (STABLE after 8 rounds)
bash runs/s45close_confirm.sh   # close_strict.py confirm runs/s45closeA_r6.txt runs/s45closeA_dips_at_r6.txt --pitch 0.001 ...
```
