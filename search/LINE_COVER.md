# Mixed covers (points + uniform line densities): the cover side (task line-cover A, 2026-09-27)

Question (`tasks/line-cover/A.md`): does moving cover weight from discrete points to uniform densities on
segments (the grid lines) shrink the validity overhead (converged LP → honest cost)?  Gate: an `m = 5` mixed cover
with honest cost `< 20.947` (`= 21 / 1.0025`) on the strict protocol.  Format: `tasks/line-cover/FORMAT.md`
(mixed v1), reader `search/mixed_cover.py`.  Labels as in `S21_COVER.md`: **heuristic** (float LP over sampled rows,
float scans; no bound either way), **certified** (exact).  **Nothing in this note is certified**; the exact word is
agent B's `search/zm_mixed.py`.

## 0. Answer

> **GO (heuristic).**  `runs/line-cover_m5_candidate.txt`: **total `20.832252` (`= 520806311/25000000`)**, 7,536
> points + 1,872 segments, exactly D4-symmetric, `1.0025 ×` over its confirmed minimum.  That is `0.168` (`0.8 %`)
> below 21 and `0.115` below the gate `20.947`.
>
> | `m` | cover | LP | strict honest cost (p 0.002) | confirmed (p 0.001) | overhead over LP |
> |---|---|---|---|---|---|
> | 4 | point covers (`S21_COVER.md` §2: `closed4_best`, r5 LP weights) | `12.34` | `12.58 – 12.70` | — | `1.9 – 2.9 %` |
> | 4 | point control, **same pipeline as the mixed runs** (`m4ptsA`, 3 closing rounds) | `12.344` | `12.654 – 12.683` | — | `2.5 %` |
> | 4 | **mixed** `m4germ` round 21 (stable) | `12.3296` | `12.331` | **`12.3346`** (min `0.999594`) | **`0.04 %`** |
> | 5 | point covers (`s5convR3` final, best closing loop) | `20.73 / 20.80` | `21.016` | — | `1.4 %` |
> | 5 | **mixed** `m5germC` round 12 (stable) | `20.7488` | `20.763` | **`20.7803`** (min `0.998489`) | **`0.15 %`** |
> | 5 | mixed `m5nolpA` (no germ rows, 12 closing rounds, stopped) | `20.7491` | `20.756 – 20.83` | — | `0.03 – 0.4 %` |
>
> * **The validity overhead collapses from `1.4–2.9 %` to `0.04–0.15 %`.**  The converged mixed LP is the same as the
>   point LP (`12.33` vs `12.34`, `20.749` vs `20.73–20.74`): the line densities cost nothing in LP value, they only
>   remove the gap between the sampled LP and the honest cost.
> * **The near-tile holes at `θ ≲ 1/D` disappear.**  On the final covers the near-tile scan (every tile `± 0.015` at
>   pitch 0.0005, `0–1.5°` step 0.01) and the tile-germ blow-up family (§2) are at `≥ 0.99976` (`m = 5`) and
>   `1.0000006` (`m = 4`); the point covers had `0.845–0.964` holes there.  The residual dips are tilted interior /
>   wall-band poses (`9.13°`, `31.9°`, `14.9°`), the kind the closing loop handles routinely.
> * The whole budget used: `21 / 20.7488 = 1.0121` over the LP; spent: validity `1.0015` (confirmed) + checker margin
>   `1.0025` (the `m = 4` figure of `M4_MARGIN.md`, assumed to carry over) → `20.832`.  Remaining risk: unfound dips
>   (the pitch-0.001 confirmation moved the minimum by `0.08 %` below the in-loop strict value, `0.9993 → 0.9985`) and
>   whether B's exact checker terminates on segments at the germs; the file leaves `0.8 %` for both.
> * Piece counts: `m = 4` 1,572 points (740 with weight `> 10⁻⁴`) + 1,136 segments; `m = 5` 7,536 points (5,712 with
>   weight `> 10⁻⁴`) + 1,872 segments, vs 3,621 (shipped `m = 4`) / 4,800–9,548 (`m = 5`) points.  The segments carry
>   `94 %` (`m = 4`) / `83 %` (`m = 5`) of the weight; the points are mostly an interior-point-LP residue of tiny
>   weights and could be thinned (not tried).

## 1. Method

Code: `search/line_cover.py` (new), `search/mixed_cover.py` (new; the reader shared with B).  Reused: `closed4.py`
(`captured`, `price`, `column_lattice`, `angle_list`), `close_shifted.py` (`scan_shift`, `stencil`, the loop recipe),
`family_rows.py` (`item3`).

* **Columns** (D4 orbits): points (lattice 0.05 + the supports of `closed4_best` / `s5convR3` / `s5convD3`, plus dual
  pricing each phase-A round), and segment pieces `[k/q, (k+1)/q]`, `q = 50`, on the interior grid lines
  `x, y ∈ {1, …, m−1}`, mass uniform by length.  Row coefficient of a piece = length fraction of the piece inside
  the **closed** square (a piece on an edge of `Q` counts in full, at exact `θ = 0` poses).
* **`--no-line-points`**: point columns on the interior grid lines are dropped, so the lines carry only densities.
  All reported runs use it (the variant with line points allowed, `m4mixA`, gave the same LP in phase A but its LP
  solves were 2× slower and it was stopped at round 6).
* **Float evaluation**: points by `closed4`'s cumulative-sum scan; each line by its cumulative mass function `F`
  (piecewise linear) at the two ends of the line's chord in `Q` (`line_cover.chord`) — so a scan costs the same as for
  points.  Cross-checked against `mixed_cover.mass_in_square_float` (generic clipping; agreement to the `10⁻⁸`
  export rounding) and the LP matrix (`A x` = scan value to `10⁻¹⁴`).  Calibration: on
  `certificates/rung2/s13_closed_cover_4.txt` the strict evaluator gives `1.0201154` (= `M4_MARGIN.md`'s `--nproc 20` value).
* **Loop**: phase A (cutting planes: 0.01 lattice at `angle_list(0.5)` with random sub-pitch shifts, polish, near-tile
  and item-3 pools, point pricing; until the LP moves `< 5·10⁻⁴`), then phase B = the `close_shifted` recipe: every
  round **measures** the exported cover by the strict protocol (`hardscan` at pitch 0.002 + `family`, identical angle
  lists and polish), plus a near-tile scan and the germ family (§2); permanent dip rows below 0.995; shifted-lattice
  separation; HiGHS IPM re-solve.  Stop: strict minimum stable `< 0.1 %` over 3 rounds.
* **Confirmation** (`line_cover.py confirm`): strict at 0.002, near-tile boxes at pitch 0.0005, germ family, boxes
  `± 0.08` at pitch 0.001 around the 80 deepest dips, the whole container at pitch 0.001 at the 319 hardscan angles
  and the 318 interleaved ones, polish of the 1,400 worst.

## 2. Why the tile holes go away, and the one family that still needs rows

At a tile `[i, i+1] × [j, j+1]` tilted by a small `θ` with centre offset `θ·(t_x, t_y)`, the line `x = i` meets the
square in `y ∈ [c_y + t_x, c_y + 1/2]` and `x = i + 1` in `[c_y − 1/2, c_y + t_x]` (to `O(θ)`), for `|t_x| ≤ 1/2`.
With **equal** uniform densities on the two lines the two chords always add to one full edge: the capture is
independent of the offset/angle ratio, i.e. the `θ → 0⁺` discontinuity is a single jump (down from the `θ = 0`
value, where all four edges count) and there is nothing for a hole to slip through.  With discrete line points the
same sum is a step function of `t_x / (point spacing)`, which is where every point cover's near-tile holes came from.

With piecewise-constant (not equal) densities the capture is **piecewise linear in `(t_x, t_y)`** as `θ → 0`, additive
(vertical lines depend on `t_x`, horizontal on `t_y`), with vertices where `c ± t` hits the `1/q` grid.  So a finite
row family samples the germ completely (up to `O(θ)` and the few off-line points): `line_cover.germ_poses` = the 2-D
blow-up `(i + ½ + θ t_x, j + ½ + θ t_y, θ)` for `t ∈ (1/q)ℤ ∩ [−½, ½]` at every tile, plus the sliding version (`c_y` on
the `1/q` grid, `t_x` on the grid).  Without these rows (`m4nolpA`, `m5nolpA`) the strict minimum kept oscillating by
`0.1–0.3 %` between rounds at germ poses `(1.5, 0.5, 0.002°)`, `(3.5, 4.5, 0.003°)` (the 1-D family alone did not
reproduce them — the two axes have to sit at vertices simultaneously); with them as permanent rows from the start
(`--germ`: 21k rows at `m = 5`) both loops became stable (`m4germ`: last four strict minima `0.99981–0.99992`;
`m5germC`: `0.99806–0.99930`) and the near-tile + germ scans measure `≥ 0.9998` from round 6 on (`1.0000` in the last two rounds).

## 3. Runs

All `taskset -c 0-9,16-25`; machine memory in use stayed `≤ 35 GB` with three runs concurrent (limit 60 GB).  LP: HiGHS IPM without crossover,
`200–1100 s` per solve at `45–100k` rows (the bottleneck; single-threaded).

| tag | `m` | columns | rounds A + B | wall | LP (last) | strict cost per B round |
|---|---|---|---|---|---|---|
| `m4ptsA` | 4 | points only (control) | 9 + 3 | 1.6 h (stopped) | 12.3441 | 12.672, 12.683, 12.654 |
| `m4nolpA` | 4 | lines + off-line points | 5 + 10 | 2.0 h (stopped) | 12.3299 | 12.462 → 12.385 → … → 12.336, 12.351 |
| `m4germ` | 4 | same + germ rows | 5 + 17 | 3.0 h, **STABLE** | 12.3296 | 12.453, 12.369, …, 12.957 (a point at `(2,2)`, fixed next round), …, 12.331, 12.331, 12.331 |
| `m5nolpA` | 5 | lines + off-line points (warm rows from `s5convR3`) | 4 + 12 | 5.5 h (stopped) | 20.7491 | 21.317, 20.809, 20.805, 20.888, 20.794, 20.809, 20.785, 20.794, 20.781, 20.756, 20.783, 20.829 |
| `m5germC` | 5 | same + germ rows | 4 + 9 | 4.2 h, **STABLE** | 20.7488 | 20.986, 20.837, 20.793, 20.801, 20.762, 20.789, 20.776, 20.781, 20.763 |

The control was stopped after 3 closing rounds; the historical point loops (`S21_COVER.md`, `M4_MARGIN.md`) ran
30–67 rounds without getting below `1.4 %` at `m = 5` or `1.3 %` (strict) at `m = 4`, and the control's first
closing round (`2.7 %`) vs the mixed run's first (`1.1 %`, `m = 4`) already shows the split.

**Confirmation of the finals** (pitch 0.001, §1):

| file | total | in-loop strict min | confirmed min | at | cost |
|---|---|---|---|---|---|
| `runs/lc_m4germ_r21.txt` | 12.329582 | 0.99985 | **0.9995941** | `(1.6116, 2.4089, 14.876°)` (interleaved angle) | 12.3346 |
| `runs/lc_m5germC_r12.txt` | 20.748862 | 0.99930 | **0.9984894** | `(3.5590, 0.5730, 9.132°)` (interleaved angle; D4 image `(4.427, 3.559)`) | 20.7803 |

Next dips at `m = 5`: `0.99874` at `(2.337, 1.689, 31.88°)`, `0.99906` at `(2.581, 3.601, 13.36°)`; near-tile boxes
`0.99976` at `(1.501, 1.4995, 0.10°)`; germ family and item-3 `≥ 1.0000`.  All residual dips are between the hardscan
angles (the interleaved set found the two deepest), i.e. ordinary tilted cells, not germs.

## 4. Where the weight sits (`m = 5` final)

Segments `17.27` (`83 %`): per interior line `2.186` (`x = 1, 4`) and `2.131` (`x = 2, 3`) — the same line masses as
the point cover's (`S21_COVER.md` §3: `2.166 / 2.132`).  Density along a line: median `0.42` per unit, peaks `≤ 4.6`
next to the corner squares.  Points `3.48`: `1.30` in the wall bands, `0.76` within `0.05` of a grid line (the old
curves), the rest interior off-grid.  `m = 4`: segments `11.59` (`94 %`), points `0.74`.

## 5. Files

| file | what |
|---|---|
| `runs/line-cover_m5_candidate.txt` | **the `m = 5` candidate**, mixed v1, `lc_m5germC_r12 × 10025000/9984894`, total 20.832252 |
| `runs/line-cover_m4_candidate.txt` | `m = 4` candidate, `lc_m4germ_r21 × 10025000/9995941`, total 12.365441 |
| `runs/lc_{m4germ,m5germC,m4nolpA,m5nolpA,m4ptsA,m4mixA}.log / .json`, `runs/lc_*_r{R}.txt`, `runs/lc_*_dips.txt` | loop logs, every round's cover (mixed v1, masses rounded up), permanent dip poses |
| `runs/lc_{m4germ_r21,m5germC_r12}_confirm.log` | the confirmations |

Candidates are heuristic: minimum `≈ 1.0025` on every float scan run here, not a proof.  Both are exactly
D4-invariant (orbit columns; masses rounded up per atom).

## 6. Reproduce

```
taskset -c 0-9,16-25 python3 search/line_cover.py loop m5germC --s 5 --q 50 --no-line-points --germ \
    --cols-from runs/s5convR3_final_last.txt,runs/s5convD3_it9_last.txt --warm runs/s5convR3_final_last.txt \
    --nproc 10 --lp-rounds 10 --lp-tol 5e-4 --prune-at 90000 --prune-keep 20000 --rounds 30 --seed 3
taskset -c 0-9,16-25 python3 search/line_cover.py loop m4germ --s 4 --q 50 --no-line-points --germ \
    --cols-from runs/closed4_best.txt --nproc 6 --lp-rounds 14 --prune-at 90000 --prune-keep 20000 --rounds 25
python3 search/line_cover.py confirm runs/lc_m5germC_r12.txt --dips runs/lc_m5germC_dips.txt --germ-q 50 --nproc 12
python3 search/line_cover.py scale runs/lc_m5germC_r12.txt runs/line-cover_m5_candidate.txt --factor 1.00401667   # 1.0025/0.9984894
python3 search/line_cover.py eval runs/line-cover_m5_candidate.txt --pitch 0.002 --tiles --germ-q 50   # strict 1.0031762
python3 search/mixed_cover.py runs/line-cover_m5_candidate.txt                                          # well-formedness + total
# control: the same loop with --q 0 (points only); m4nolpA / m5nolpA: without --germ
```
