# `COVER^closed(6)`: the cover side of `s(32) = 6` (task s32-cover; 2026-09-24/25)

Route: `s(32) = 6` from a weighted closed point cover of `[0,6]²` of total `< 32`. This is the `m = 6`
rung of `n = m² − 4` in the `s(13) = 4` machinery (`RUNG2.md`), and the next rung after `S21_COVER.md` (`m = 5`, which failed
narrowly).  This note covers only the cover side: the converged heuristic LP, the honest cost of explicit covers, and
where the weight sits.  No exact zero-margin check was started.  Labels: **heuristic** (sampled-row float LP, float
scan, local search; no bound in either direction), **certified** (exact `Fraction`), **estimate** (extrapolation).
Nothing here is certified.

Code: `search/closed4.py` (unchanged), `search/rung2_close.py` (extended, §6), `search/grow_cover.py` (new: the
`m → m+1` transplant and the near-tile pool), `search/s21_cover_eval.py` / `scale_cover.py` / `family_rows.py`
(unchanged; already `s`-generic).  Cores 0–6 + 16–22 (`--nproc 7`), 2026-09-24 17:45 → 2026-09-25 04:10
(10.5 h wall).  Peak RSS was about 4 GB per run.

## 0. Answer

> **Update 2026-09-26 (exact-check pilots, `S32_EXACT.md`):** everything except the interior tile germs certifies
> with `--branch-cap 640`; the germs (`|θ| ≤ 0.56°`, within 0.003 of a tile centre) remain open, with no coverage
> hole found (float min 1.010).  No full sweep yet.
>
> **Update 2026-09-25 (task s32-close, §9):** a closing loop measured by the strict protocol itself reached a
> stable strict minimum `0.99902` (pitch 0.002).  An off-grid pitch-0.001 confirmation then found `0.99319`.  The
> candidate `runs/s32-close_candidate.txt` totals **`31.7135`**: that is `1.004 ×` over the confirmed minimum, and
> `0.29` (`0.9 %`) below 32.  Verdict: **go** for the exact checker (heuristic; see §9 for the risk).

> **Converged LP (heuristic): `COVER^closed(6) ≈ 31.22 ± 0.03`.**  Three chains agree to within `0.04`, and each was
> still drifting down when stopped:
> * the transplant chain on coarse rows (`s32convT`, from the `m = 5` LP cover grown by one cell): `31.196`;
> * the cold chain (`s32convB`, lattice columns): `31.214`, falling `0.004`/round;
> * the fine-row chain (`s32convD`, row pitch 0.01, 104 angles): `31.238`, falling `0.005`/round, dual coverage `1.036`.
>
> **So `32 − LP ≈ 0.78`, i.e. `2.4 %`** (at `m = 5` it was `0.27`, `1.3 %`).  Break-even factor `32 / 31.22 = 1.025`
> (at `m = 5`: `1.013`).
>
> **Honest cost (heuristic, total / min capture).**  The best closing-loop cover is `s32convR4` round 11:
> 12,025 points, total `31.2817`, LP `31.2811`.
> * Strict protocol (`s21_cover_eval.py hardscan`): **`31.700`** at pitch 0.004 (min `0.98682`), and
>   **`31.801`** at pitch 0.002 (min `0.98368`, a `17.75°` wall-band dip that the 0.004 lattice steps over).
>   The best strict honest cost is **`31.80`**.
> * `family`: min `0.99365`.
> * Lenient `m = 4` protocol (`closed4.py stress`): `31.608`.
>
> For comparison, the `m = 5` best is `21.016` / `21.019` at pitch 0.004 / 0.002.  The plain converged-LP cover
> (`s32convD` it6) costs `32.364` on the strict protocol (`m = 5`: `21.567`).
>
> **Verdict: `W < 32` is plausible, unlike `W < 21`.**  The strict honest cost sits `0.20` (`0.63 %`) below 32.  That
> is `×1.0186` over the converged LP; at `m = 5` it was `×1.0139`, still above 21.  With the known overhead components:
> * validity (the closing loop) has cost `1.4–1.9 %` at `m = 5, 6`;
> * the exact checker certified the whole of `[0,4]²` at `0.25 %` margin (`M4_MARGIN.md`), with a floor of `~0.18 %`;
> * `31.80 × 1.0025 = 31.88 < 32`.
>
> What the exact checker would have to hit: a total factor **`≤ 32 / 31.2817 = 1.0230`** over the `R4 r11` weights,
> covering validity and margin together.  Of that, `≥ 1.0166` (`1/0.98368`) is already spent on the dips the scans
> have found, which leaves **`≤ 0.63 %`** for unfound dips plus checker margin.  This is tight but not hopeless: it is
> 2.5× the `m = 4` checker margin.  The risk is validity, not margin.  The pitch-0.002 rescan moved the minimum by
> `0.3 %` here (`0.98682 → 0.98368`), where at `m = 5` it moved it by `0.01 %`.  So the `R4` weights are less
> "settled" than the `m = 5` ones, and further closing-loop rounds whack-a-mole among `10–18°` wall-band and `0.1–0.7°`
> edge-tile dips (§2).  Recommendation: a longer closing loop with the pitch-0.002 scan (or the exact checker as the
> separation oracle) before the zero-margin run.

## 1. Runs

All `taskset -c 0-6,16-22`, `--nproc 7`, `--s 6 --no-literature`.
* `seedrows` = `runs/s32_seedrows.txt`: `family_rows.py --s 6 --pitch 0.004`, 149,892 item-3 poses.
* `neartile` = `runs/s32_neartile.txt`: `grow_cover.py neartile --s 6`, 218,700 poses (the 36 tile centres
  `± 0.012` in steps of 0.003, `θ = 0.02..1.5°` in steps of 0.02°, as at `m = 5`).
* Both pools were lazy (`--lazy-seed`) in every run.
* `transplantD3` = `runs/s32_transplantD3.txt`: `grow_cover.py transplant runs/s5convD3_it9_last.txt`.  It is the
  `m = 5` converged-LP cover with its centre cell column and row duplicated (`x → {x : x < 3} ∪ {x+1 : x ≥ 2}`):
  16,081 points, total `31.047` (the cell model), D4-symmetric, lattice min `0.78` at the seams.

| tag | what | rounds / wall | LP (last) | dual cov. | notes |
|---|---|---|---|---|---|
| `s32convT` | `--from-cert transplantD3 --warm-lp`, pitch 0.02, 59 angles, `--cg-want 400` | 9 / 2.8 h | **31.1964** | 1.040 | LP solve grew to 2600 s/round; stopped |
| `s32convB` | cold, lattice columns (pitch 0.05), same rows | 12 / 1.8 h | 31.2139 | 1.185 | falling `0.004`/round; stopped (T dominated) |
| `s32convD` | `--from-cert` T it7, **row pitch 0.01, 104 angles** (`--deg-step 0.5 --block 0.05`) | 7 / 4.7 h | **31.2383** | 1.036 | first cold solve on 75k rows took 3400 s; `D it6` = "the converged-LP cover" |
| `s32convR1` | `rung2_close.py --near-tile` from T it7 (no pruning) | 2 / 0.7 h | 31.2555 | — | 88k rows by round 1; replaced by R2 |
| `s32convR2` | same from T it8, `--prune-at 50000` | 9 / 2.2 h | 31.2741 | — | pruning dropped the narrow near-tile rows, and the `0.3–0.6°` holes came back (0.955) |
| `s32convR3` | from R2, + `--fine-tilted --lazy-seed` both pools, protected rows | 5 / 0.9 h | 31.2786 | — | column set was R2's support only (1,483 orbits): replaced by R4 |
| `s32convR4` | `rung2_close.py` from **D it6** (1,983 orbits / 14,397 atoms), `--near-tile --fine-tilted --lazy-seed neartile,seedrows --prune-at 60000 --prune-keep 25000 --prune-slack 0.003 --init-thr 0.005` | 17 / 3.7 h | 31.2870 | — | **the best cover is round 11**, `runs/s32convR4_r11_last.txt` |
| `s32convR5` | R4 without pruning | 2 / 1.1 h | — | — | identical trajectory to R4 (same seeds), 45 min/round; killed |

Convergence of the LP (heuristic; `closed4.py run` columns):

| run, round | LP | latmin | viol | dualcov |
|---|---|---|---|---|
| B it2 | 31.2782 | 0.873 | 76,163 | 1.20 |
| B it7 | 31.2620 | 0.973 | 2,520 | 1.075 |
| B it11 | 31.2139 | 0.983 | 776 | 1.185 |
| T it3 | 31.2020 | 0.953 | 14,032 | 1.22 |
| T it6 | 31.2022 | 0.994 | 704 | 1.106 |
| T it8 | 31.1964 | 0.987 | 587 | 1.040 |
| D it2 (fine rows) | 31.2625 | 0.988 | 1,632 | 1.127 |
| D it4 | 31.2484 | 0.993 | 2,272 | 1.081 |
| D it6 | **31.2383** | 0.991 | 1,827 | 1.036 |

The fine rows sit `+0.04` above the coarse ones (`m = 5`: `+0.016`).  None of the chains is converged in the strict
sense: dual coverage never reached 1.00, and the lattice `viol` never reached 0.  All three were still falling slowly.
Reading: **`31.22 ± 0.03`**.  The transplant is cheap and was the best warm start: the transplanted `m = 5` weights
are already at `31.05` (below the LP, because they are invalid at the seams), and the first LP round gave `30.81`.

## 2. Honest cost (heuristic)

The same three float checks as `S21_COVER.md` §2 (`stress` lenient; `hardscan` strict; `family` = the item-3
family at pitch 0.0005).  New here: `hardscan --pitch 0.002`, because at `m = 6` the 0.004 lattice missed a dip by `0.3 %`.

| cover | total | `stress` min / cost | `hardscan` 0.004 min / cost | `hardscan` 0.002 min / cost | `family` min | status |
|---|---|---|---|---|---|---|
| `m = 5` `s5convD3_it9` (converged LP) | 20.7372 | — | 0.96153 / 21.567 | — | 0.96584 | heuristic (S21_COVER) |
| `m = 5` **`s5convR3` final** | 20.8014 | 0.99115 / 20.987 | 0.98978 / **21.016** | 0.98967 / **21.019** | 0.99947 | heuristic |
| `m = 6` `s32convT_it8` (coarse LP) | 31.1971 | — | 0.95842 / 32.551 | — | 0.97449 | heuristic |
| `m = 6` `s32convD_it6` (converged LP) | 31.2390 | — | 0.96525 / 32.364 | — | 0.97740 | heuristic |
| `m = 6` `s32convR4` r6 | 31.2755 | 0.97964 / 31.926 (in-loop) | 0.98000 / 31.914 | — | 0.99298 | heuristic |
| `m = 6` `s32convR4` r8 | 31.2798 | 0.98598 / 31.724 (in-loop) | 0.97798 / 31.984 | — | — | heuristic |
| `m = 6` **`s32convR4` r11** | **31.2817** | 0.98966 / **31.608** (`closed4.py stress`); 0.98772 in-loop | 0.98682 / **31.700** | 0.98368 / **31.801** | 0.99365 | heuristic |
| `m = 6` `s32convR4` r16 (last) | 31.2876 | 0.97877 / 31.966 (in-loop) | 0.97600 / 32.057 | — | — | heuristic |

(The cost is total / min.  The true minimum capture is at most the scanned one, so the true cost of each point set
is `≥` the listed cost.)

Where the dips are, all near-wall:
* the edge tile nudged by `0.1–0.7°`, e.g. `(0.512, 3.498, 0.68°)`, `(2.501, 4.499, 0.12°)`, `(1.498, 4.499, 0.35°)`
  (the `S21_OVERHEAD.md` binding pose, as at `m = 4, 5`);
* wall-band squares tilted `9–18°` near a corner, e.g. `(0.579, 1.445, 9.95°)`, `(0.629, 1.457, 17.8°)`, and
  `(4.543, 0.629, 17.75°)`, the pitch-0.002 minimum of r11;
* interior tilted squares at `30–41°`, e.g. `(2.295, 4.633, 41.4°)`, `(3.692, 1.377, 35.5°)`.

The closing loop is whack-a-mole, exactly as at `m = 4, 5`.  The in-loop stress min oscillated between `0.950` and
`0.988` over rounds 4–16, and the LP rose monotonically by `+0.06` (from `31.221` to `31.287`).  Each dip costs
`≈ 0` LP once it is a row; the difficulty is finding them.  Row pruning (needed to keep the LP solve under ~10 min)
makes it worse: a dropped narrow row re-opens its hole (`R2`, `R4` rounds 4 and 9).  §6 describes the protection
that limits this.

## 3. Where the weight sits

`s21_cover_eval.py regions` (exact sums of the exported rational weights).

| cover | corner cell (mean) | edge cell (mean) | interior cell (mean) | total |
|---|---|---|---|---|
| `m = 4` `closed4_best.txt` | 0.4977 | 0.7775 | 1.0517 | 12.4175 |
| `m = 5` `s5convD3_it9` (LP) | 0.4997 | 0.8048 | 1.0090 | 20.7372 |
| `m = 6` `s32_transplantD3` (m = 5 grown, invalid) | 0.4997 | 0.8090 | 1.0066 | 31.0473 |
| `m = 6` `s32convT_it8` (coarse LP) | 0.4955 | 0.8206 | 1.0054 | 31.1971 |
| `m = 6` **`s32convD_it6` (LP)** | **0.5035** | **0.8239** | **1.0026** | 31.2390 |
| `m = 6` `s32convR4` r11 (closing loop) | 0.5037 | 0.8191 | 1.0101 | 31.2817 |

`s32convD_it6`, per cell (rows = `y` cells, top first):

```
0.504  0.814  0.834  0.834  0.814  0.504
0.814  1.003  1.004  1.004  1.003  0.814
0.834  1.004  1.000  1.000  1.004  0.834
0.834  1.004  1.000  1.000  1.004  0.834
0.814  1.003  1.004  1.004  1.003  0.814
0.504  0.814  0.834  0.834  0.814  0.504
```

* **Corners:** exactly `4 × 1.0000` inside the closed corner squares, as at `m = 4, 5`.
* **Interior:** the four cells that do not touch the edge ring are at exactly `1.000`, as the `m = 5` centre cell
  was.  The interior mean falls further, `1.009 → 1.003`.
* **Edge cells *rise* again, `0.805 → 0.824`.**  The cells next to the middle of each wall carry `0.834`.  This is the
  whole miss of the cell model: `S21_COVER.md` §3 froze `b = 0.80` and predicted `LP ≈ 31.0`, but the actual edge
  mean is `0.824`, i.e. `+0.02 × 16 = +0.3`.  Class sums: wall bands `14.05` (45.0 %), interior `13.19` (42.2 %).
* By distance to the wall: `0<d<1` 39.4 %, `d = 1` 18.3 %, `1<d<2` 25.2 %, `d = 2` 8.4 %, `d > 2` 8.6 %.
* On the integer grid lines: 76.6 % (`m = 5`: 82 %); within 0.05 of one: 6 %; elsewhere: 17 %.
* Updated cell model (estimate), with `a, b, c = 0.50, 0.82, 1.00`:
  `LP(m) − (m² − 4) ≈ 2.0 − 0.72(m − 2) + 0.003(m − 2)²`, i.e. `−0.27` at `m = 5` (measured), `−0.78` at `m = 6`
  (measured), and `≈ −1.5` at `m = 7` if the edge mean stops rising.  It rose by `0.02` per rung so far; if it keeps
  rising, `≈ −1.1`.

## 4. Room below 32

| quantity | value | status |
|---|---|---|
| converged cover LP | `31.22 ± 0.03` | heuristic |
| room `32 − LP` | `0.78` (`2.4 %`) | heuristic |
| best honest cost, lenient (`closed4.py stress`) | `31.608` (`R4` r11) | heuristic |
| best honest cost, strict (`hardscan` 0.004 / 0.002) | `31.700` / **`31.801`** (`R4` r11) | heuristic; the true cost is `≥` this |
| strict / LP | `31.80 / 31.22 = 1.0186` (`m = 5`: `21.019 / 20.73 = 1.0139`) | heuristic |
| room left for the checker | `32 / 31.80 = 1.0063` | heuristic |
| factor the exact checker must hit over `R4` r11 (validity × margin) | `≤ 1.0230` | heuristic (total `31.2817` is exact) |
| `m = 4` checker margin (whole square) | `1.0025` (`M4_MARGIN.md`) | certified sweep, not of record |
| `m = 6` at the `m = 4` shipped factor | `31.22 × 1.05 = 32.78` | estimate |

## 5. Performance notes

* The LP solve dominates, and it is single-threaded HiGHS.  Its time is erratic: `130 s` right after a row prune,
  then `1–2.6 ks` a few rounds later at 40k rows × 3–4k orbits.  The cold first solve of a `--from-cert` restart on
  70–75k rows takes `25–60 min`.  About 3.5× the `m = 5` times, roughly (area ratio)².
* `rung2_close.py` without pruning: 15 min/round at 88k rows, growing.  With `--prune-at 60000`: 8–17 min/round.
* Two runs sharing the 14 logical CPUs slow each other only during the scans; the LP phases barely interact.
* The transplant warm start plus lazy pools reached the LP plateau in 5 rounds (1 h).  The cold start took 10 rounds.

## 6. Code changes (backwards-compatible; defaults reproduce the old behaviour)

`search/grow_cover.py` (new):
* `transplant IN OUT`: `m → m+1` by duplicating the middle cell column and row; asserts D4 symmetry.
* `neartile OUT --s S`: the near-tile pose pool.  At `s = 5` the defaults give the 151,875-pose
  `runs/s5conv_neartile.txt` recipe; that file was written by an unsaved scratch script, so this is a
  reconstruction and is not re-verified here.

`search/rung2_close.py`:
* `--cap` defaults to `13` at `s = 4` (as before) and `s² − 3.5` otherwise.
* `--prune-at / --prune-keep / --prune-slack`: row pruning as in `closed4.py run`.  It keeps the dual support and
  the rows with `Ax ≤ 1 + slack`, plus the newest `--prune-keep`.
* **Protected rows**: stress seeds, polished worst poses and lazy-pool violations are never pruned.
* `--lazy-seed FILES`: pose pools as a separation oracle, as `closed4.py --lazy-seed`.
* `--init-thr` (default `0.02`, as before): threshold for the initial lattice rows.
* With pruning or a lazy pool, the seeds are added before the polish-visited poses so that the protection is exact.
  Without them the old order is kept.

`s21_cover_eval.py`, `scale_cover.py`, `family_rows.py`, `closed4.py`: checked; no `s = 4 / 5` hard-coding on the paths used.

## 7. Files (`runs/`, gitignored)

| file | what |
|---|---|
| `s32_seedrows.txt`, `s32_neartile.txt`, `s32_transplantD3.txt` | pools and the transplant |
| `closed4_s32conv{T,B,D,R1..R5}.{log,json}`, `s32conv*.nohup.log` | all runs |
| `s32convT_it{7,8}_last.txt`, `s32convB_it11_last.txt` | coarse LP covers |
| **`s32convD_it6_last.txt`** | the converged-LP cover: 14,397 points, total 31.239035 |
| **`s32convR4_r11_last.txt`** | the best closing-loop cover: 12,025 points, total 31.281679 |
| `s32convR4_r{6,8,16}_last.txt`, `s32convR2_stop_last.txt`, `s32convR3_r5_last.txt`, `s32convR1_r1_last.txt` | other checkpoints in §2 |
| `s32_eval_*.log` | the hardscan / family logs |

All `_last.txt` files are in certificate format (`s = 6`, `D = 1000`, `W = 10⁷`, weights rounded up) and D4-symmetric.
None is a valid cover (min capture `< 1`).  Scale by `1/min` or more (`scale_cover.py`) before any exact check;
e.g. `R4 r11 × 1.0170` gives `31.813`.

## 8. Reproduce

```
python3 search/grow_cover.py transplant runs/s5convD3_it9_last.txt runs/s32_transplantD3.txt
python3 search/grow_cover.py neartile runs/s32_neartile.txt --s 6
python3 search/family_rows.py runs/s32_seedrows.txt --s 6 --pitch 0.004
P=runs/s32_seedrows.txt,runs/s32_neartile.txt
python3 search/closed4.py run --s 6 --tag s32convT --from-cert runs/s32_transplantD3.txt --no-literature --warm-lp \
        --lazy-seed --seed-rows $P --cg-want 400 --prune-at 60000 --prune-keep 20000 --col-prune-at 3500 --nproc 7 --time 12600
python3 search/closed4.py run --s 6 --tag s32convB --col-pitch 0.05 --no-literature --warm-lp --lazy-seed --seed-rows $P \
        --cg-want 400 --prune-at 60000 --prune-keep 20000 --col-prune-at 5000 --nproc 7 --time 12600      # stopped at it11
python3 search/closed4.py run --s 6 --tag s32convD --from-cert runs/s32convT_it7_last.txt --pitch 0.01 --deg-step 0.5 \
        --block 0.05 --no-literature --warm-lp --lazy-seed --seed-rows $P --cg-want 400 --prune-at 60000 \
        --prune-keep 20000 --col-prune-at 3500 --nproc 7 --time 14400                                     # stopped after it6
python3 search/rung2_close.py runs/s32convD_it6_last.txt s32convR4 --warm-lp --near-tile --fine-tilted --nproc 7 \
        --nrand 40 --init-thr 0.005 --lazy-seed runs/s32_neartile.txt,runs/s32_seedrows.txt --prune-at 60000 \
        --prune-keep 25000 --prune-slack 0.003 --time 25200                                               # r11 copied out; stopped at r16
python3 search/s21_cover_eval.py hardscan runs/s32convR4_r11_last.txt                 # 0.9868166 -> 31.700
python3 search/s21_cover_eval.py hardscan runs/s32convR4_r11_last.txt --pitch 0.002   # 0.9836802 -> 31.801
python3 search/s21_cover_eval.py family   runs/s32convR4_r11_last.txt                 # 0.9936483
python3 search/closed4.py stress runs/s32convR4_r11_last.txt --nproc 7 --tag s32stress   # 0.9896608 -> 31.608
python3 search/s21_cover_eval.py regions  runs/s32convD_it6_last.txt
```

## 9. Closing to a stable strict minimum (task s32-close, 2026-09-25)

Everything in this section is **heuristic**: float LP over sampled rows, float scans, local search.  Code: `search/close_shifted.py`
(new; see its docstring).  `closed4.py`, `rung2_close.py` and `s21_cover_eval.py` are unchanged.  Cores 0–6 + 16–22.
The loop used 14 workers and ran 09:10 → 14:00 (4.8 h).  The confirmation scan used 10 workers and ran 14:02 → 14:28.

**Loop.**  `close_shifted.py loop`.  The columns are the R4 r11 support plus the `D it6` columns (1,985 orbits / 14,397
atoms), fixed, with no pricing.  Each round runs these steps:
1. **Measure** the exported cover with exactly the strict protocol, at pitch 0.002: the `hardscan` lattice plus the
   polish of the 1,400 worst poses (same chunking and seeds as `s21_cover_eval.py`), and `family` at pitch 0.0005.
   Round −1, the R4 r11 cover itself, reproduces the known `0.9836802` / `0.9936483` exactly.
2. **Find new rows** on a lattice of the same pitch with a random sub-pitch offset and angles jittered within their
   step, drawn fresh each round, so the rows never sit on the measurement lattice.
3. **Keep permanently** every dip below 0.995 and never prune it.  This covers the polished dips, the worst
   polish-visited poses, the family poses, and the 4 worst shifted-lattice poses per angle.  A 3×3 stencil at
   ±0.001 (±0.01°) around the 40 deepest dips is also kept.  Everything is logged to `runs/s32closeA_dips.txt`.
4. **Re-solve** with a cold HiGHS IPM; each solve took 75–206 s.  Prunable rows are dropped at 60k (`--prune-keep 25000`).

It stops once the strict minimum changes by less than 0.1 % over 3 consecutive rounds.

| round | total | strict min p = 0.002 (hardscan ∧ family) | cost | where the minimum sat | rows (permanent) |
|---|---|---|---|---|---|
| −1 (R4 r11) | 31.2817 | 0.98368 | 31.801 | wall band `17.75°` | 22,927 (0) |
| 0 | 31.2672 | 0.84936 | 36.813 | `(2.516, 0.558, 7.14°)` | 60,953 (2,590) |
| 1 | 31.2981 | 0.92134 | 33.970 | edge tile `0.17°` | 58,182 (5,614) |
| 2 | 31.3093 | 0.96578 | 32.419 | edge tile `0.08°` | 51,388 (8,656) |
| 3 | 31.3247 | 0.98464 | 31.813 | edge tile `0.16°` | 55,520 (11,448) |
| 4 | 31.3326 | 0.98576 | 31.785 | wall band `10.6°` | 53,032 (13,038) |
| 5 | 31.3376 | 0.98861 | 31.699 | edge tile `0.24°` | 59,194 (14,197) |
| 6 | 31.3409 | 0.99478 | 31.505 | centre tile `0.09°` | 65,347 (15,430) |
| 7 | 31.3434 | 0.98658 | 31.770 | edge tile `0.06°` | 67,881 (15,567) |
| 8 | 31.3446 | 0.99692 | 31.441 | interior `40.0°` | 73,529 (16,159) |
| 9 | 31.3459 | 0.99783 | 31.414 | interior `32.5°` | 80,107 (16,160) |
| 10 | 31.3464 | 0.99386 | 31.540 | edge tile `0.54°` | 85,399 (16,160) |
| 11 | 31.3469 | 0.99789 | 31.413 | `18.25°` | 90,753 (17,408) |
| 12 | 31.3476 | 0.99065 | 31.644 | edge tile `0.20°` | 94,757 (17,408) |
| 13 | 31.3479 | 0.99495 | 31.507 | family (item 3) | 102,724 (17,738) |
| 14 | 31.3481 | 0.98353 | 31.873 | family (item 3) | 105,141 (17,739) |
| 15 | 31.3632 | 0.97948 | 32.020 | `44.75°` near the wall | 110,232 (17,741) |
| 16 | 31.3673 | 0.98356 | 31.892 | `31.2°` / family | 115,934 (19,900) |
| 17 | 31.3687 | 0.96815 | 32.401 | edge tile `0.11°` (lattice = family) | 118,447 (22,061) |
| 18 | 31.3699 | 0.99717 | 31.459 | interior `32.0°` | 122,992 (23,501) |
| 19 | 31.3705 | 0.99834 | 31.423 | `2.36°` | 123,452 (23,501) |
| 20 | 31.3711 | 0.99866 | 31.413 | interior `37.0°` | 125,924 (23,501) |
| **21** | 31.3713 | **0.99912** | **31.399** | `4.51°` | 127,179 (23,501) |
| **22** (final) | **31.3716** | 0.99902 | 31.402 | interior `29.7°` | 130,732 (23,501) |

* **Stable** over rounds 19–22 (`0.99834, 0.99866, 0.99912, 0.99902`).  The LP rose by `+0.09` over the loop (`31.28 → 31.37`).
* Rounds 13–17 were the last whack-a-mole phase: family (item-3) poses and near-tile poses reopened, down to `0.968`.
  Once the permanent set had grown to 23.5k, no deep hole reopened.  The stop is driven by interior `30–40°` cells
  at `~0.999`.
* Pitch 0.004 was not run per round.  On the final cover r22, `s21_cover_eval.py hardscan` at pitch 0.004 gives
  `0.99932`, cost `31.393`.

**Confirmation scan (r22; `close_shifted.py confirm --tiles --interleaved --full --dips runs/s32closeA_dips.txt`).**
The scan is genuinely off-grid:

| scan | lattice min |
|---|---|
| hardscan at pitch 0.002 (re-measured) | 0.99902 |
| pitch-0.001 boxes (±0.08) around the 80 deepest distinct dips, at ±0.3° step 0.01° (tilted) or ±0.06° step 0.002° (near-tile) | 0.99734 at `(5.455, 2.515, 5.38°)` |
| every tile pose ±0.015 at pitch 0.0005, 0..1.5° step 0.01° | 0.99794 at `(4.499, 2.500, 0.07°)` |
| whole container at pitch 0.001, the 319 hardscan angles | 0.99773 at `(3.610, 4.684, 32.25°)` |
| whole container at pitch 0.001, the 318 **interleaved** angles (0.01..2.99 step 0.02, 3.125..44.875 step 0.25) | 0.99571 at `(0.516, 3.500, 0.09°)`, the edge tile |
| polish of the 1,400 worst of all the above | **0.99319** at `(3.4996, 5.4925, 0.09°)`, the top-edge tile turned 0.09° |

* **Confirmed minimum `0.99319`, cost `31.587`.**  The gap below the loop's strict minimum is **`0.58 %`**.  At
  `m = 7` the sibling measured `1.29 %`.  The binding pose is again the edge tile nudged by `~0.1°`, as at
  `m = 4, 5, 6`.  The deepest cell sits inside the ±0.015 / 0.0005 near-tile box but was missed by it, so it is
  thinner than 0.0005.  Thin cells like this are the remaining risk.
* **Candidate:** `runs/s32-close_candidate.txt` = r22 × `10109/10000`, i.e. `≥ 1.004 / 0.99319` (`1.01088`).
  Total **`31.713505`**, `0.2865` (`0.90 %`) below 32.  Its heuristic minimum capture is `1.00402`.

**Verdict: go for the exact checker (heuristic).**  Factor budget over the r22 weights:
* total allowed: `32 / 31.3716 = 1.0200`;
* already spent on the confirmed dips: `1.0069`;
* left for unfound dips plus checker margin: **`1.0131`**.  The candidate's `1.004` uses part of it, which leaves the
  candidate itself `32 / 31.7135 = 1.0090` of headroom.

If the checker finds the candidate short, rescale r22 up to `32/31.3716`.  That keeps the total `< 32` as long as the
true minimum is `≥ 0.98036`, i.e. within `1.3 %` of the confirmed minimum.  The m = 4 checker needed `0.25 %` margin
(floor `~0.18 %`).  So an unfound dip up to `~1 %` below `0.99319` still fits.  The pitch-0.002 → confirmation gap
seen here (`0.58 %`) is about half of that.

Files (`runs/`, gitignored): `s32closeA.{log,json}`, `s32closeA_r{0..22}.txt`, `s32closeA_dips.txt` (21,349 dip
poses, `cx cy θ value`), `s32closeA_confirm_r22.log`, `s32closeA_eval004_r22.log`, **`s32-close_candidate.txt`**.

```
python3 search/close_shifted.py loop runs/s32convR4_r11_last.txt s32closeA --cols-from runs/s32convD_it6_last.txt \
        --lazy-seed runs/s32_neartile.txt,runs/s32_seedrows.txt --hs-pitch 0.002 --prune-at 60000 --prune-keep 25000 \
        --nproc 14 --time 36000                                               # STABLE after round 22
python3 search/close_shifted.py confirm runs/s32closeA_r22.txt --dips runs/s32closeA_dips.txt --tiles --interleaved --full --nproc 10
python3 search/scale_cover.py runs/s32closeA_r22.txt 10109 10000 runs/s32-close_candidate.txt    # 31.713505
```
