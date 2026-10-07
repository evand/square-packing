# How much margin the exact checker needs at `m = 4`, and whether the validity gap is cheap to close (task m4margin; 2026-09-23)

Question: `S21_OVERHEAD.md` splits the `m = 4` overhead `Λ = 1.05` into validity `λ_v` (the LP weights miss
unsampled poses) and checker margin `λ_c` (the slack `zeromargin.py` needs to terminate).  At `s = 5` the whole
budget is `Λ ≤ 21/20.73 ≈ 1.013`.  This note measures the two factors separately.

Code: `search/zeromargin.py` gains `--chunksize` (root boxes per worker task, default 4 = old behaviour; it load-balances
column sweeps and cannot change the result).  Nothing else in `search/` changed; `certificates/` untouched.
Scratch: `$S = /tmp/claude-1000/-home-evand-math/f433c3ef-9956-443d-90fb-0c6d2ee917e2/scratchpad/m4margin`.
Labels: **certified** (exact `Fraction`: `zeromargin.py cert` / `pose`), **heuristic** (float scan / LP / polish),
**estimate**.

## 0. Answer

> **The checker factor is negligible; validity is the whole problem.**
>
> * **`λ_c ≤ 1.0025` over the true minimum, full container (certified).**  The shipped cover scaled by
>   `40100000/40804617` (true min exactly `1.0025` at the binding pose) verifies on all of `[0,4]²`:
>   **`W = 49479474611/3886154000 = 12.732248`, 0 uncertified, max depth 13, 29,372 boxes** (1.74× the shipped run's 16,872; 27.5 CPU-h) (§2).
>   That is a new `s(13) = 4` certificate, `0.224` lighter than the shipped `12.955972` (not shipped; `certificates/` untouched).
> * On the hardest band (`cx ∈ [1.4, 1.6]`, the binding pose) the checker terminates down to **`0.02 %`** margin at depth 13,
>   with boxes growing only `1,508 → 2,798` from `2 %` to `0.02 %` (≈ logarithmic).  At `0` margin and at `−0.02 %` it fails,
>   in both cases in ≤ 52 boxes pinned at the wall tile germ `(1.5, 0.5, θ→0)` (§1).
> * The residual `λ_c` is **structural CHAIN slack at tile germs**, box-size independent: `0.009 %` at the germ `(1.5, 0.5)`,
>   **`0.18 %`** at the germ `(0.5, 1.5)` (the full run at `0.02 %` fails there; `0.25 %` passes) (§1.3).
> * **Validity `λ_v` is not cheaply fixable with sampled rows (heuristic).**  The true min of the `r5` LP weights is
>   `≤ 0.9715385` (exact pose), i.e. `λ_v ≥ 1.0293`, honest cost `12.7005`.  Re-solving the LP with 129,600 near-tile rows
>   (`θ = 0.02–2°`) + item-3 rows + colgen gives LP `12.3369` and strict min `0.97189` (honest `12.694`): **no gain**.
>   The stress-driven closing loop (fixed columns, near-tile + fine-tilted stress) reaches in-loop stress min `0.993`, but
>   the strict scan of its checkpoint finds `0.845` at `(1.5, 3.478, 0.017°)`: every round opens a new near-tile hole at a
>   smaller angle (§3).
> * **Implied `m = 5` total:** `Λ = λ_v · λ_c` with `λ_c ≈ 1.002–1.0025`.  The best `s = 5` cover is at honest cost
>   `21.016` (`S21_COVER.md`, strict protocol) → certifiable at `≈ 21.07` (estimate).  **`W < 21` needs `λ_v ≤ 1.011`
>   over the converged LP `20.73`, against `1.0138` achieved at `m = 5` and `1.029` at `m = 4`.**
>
> **Recommendation:** stop paying for checker margin (use `×1.0025` over the true min, not `×1.05` over the LP) and spend
> everything on validity: exact-oracle separation (`rung2_oracle.sh`, scaled to `λ_c = 1.0025`) or analytic near-tile
> rows (critical angles `~1/D` rad and below at the wall tiles).  A float-row LP loop does not converge on the near-tile family (§3).

## 1. Checker factor on the hardest band (certified partial sweeps)

Cover: `certificates/rung2/s13_closed_cover_4.txt` scaled by `f = (1+ε)/μ₀`, where

    μ₀ = 40804617/40000000 = 1.020115425     (certified upper bound on its true min)

is the exact capture at the admissible rational pose `cx = 25199/50000, cy = 249609/100000, u = 2323208/581517279`
(`(0.50398, 2.49609, 0.4578°)`), found by `s21_cover_eval.py hardscan` (float min `1.0201154`, three distinct
poses on one plateau).  **Correction to the brief and to `S21_OVERHEAD.md` §1:** the shipped cover's true min is
`≤ 1.0201154`, not `1.0277`.  `1.02766` is one exact spot check; the float `0.9715 × 1.05` in `S21_OVERHEAD.md` §2 was
right.  So the `×1.035` band run there already had only `0.555 %` margin.

`zeromargin.py cert FILE --disj --chain-from 0 --cx-lo 1.4 --cx-hi 1.6` (320 roots; binding pose image `(1.5036, 0.5037,
0.427°)` and interior tile `(1.5, 1.5)` inside):

| file | `ε` = true-min margin | total | boxes | max depth | ADM / CHAIN / EMPTY | uncert. | wall (procs, chunk) | status |
|---|---|---|---|---|---|---|---|---|
| shipped `×1.05` | `2.01 %` | 12.955972 | 1,508 | 13 | 85 / 734 / 95 | 0 | 3067 s (2) | certified partial |
| `×1.035` LP (`S21_OVERHEAD`) | `0.555 %` | 12.770887 | 2,324 | 13 | 120 / 1,095 / 107 | 0 | 4715 s (2) | certified partial |
| `c_e50` | `0.5 %` | 12.763999 | 2,358 | 13 | 120 / 1,112 / 107 | 0 | 1435 s (20, 4) | certified partial |
| `c_e25` | `0.25 %` | 12.732248 | 2,544 | 13 | 122 / 1,201 / 109 | 0 | 1900 s (20, 4) | certified partial |
| `c_e10` | `0.1 %` | 12.713197 | 2,730 | 13 | 127 / 1,287 / 111 | 0 | 1766 s (20, 4) | certified partial |
| `c_e2` | `0.02 %` | 12.703036 | 2,798 | 13 | 132 / 1,315 / 112 | 0 | 1094 s (20, 1); 274 CPU-min | certified partial |
| `c_e0` | `0` | 12.700496 | 3,130 | **24** (limit) | 166 / 1,383 / 124 | **52** | 1390 s; 337 CPU-min | NOT VERIFIED |
| `c_em2` | `−0.02 %` | 12.697956 | 2,956 | **18** (limit) | 147 / 1,355 / 115 | **21** | 1275 s; 327 CPU-min | NOT VERIFIED (correct: the cover is short) |

(`ε` is exact at the binding pose; the other rows' `ε` are upper bounds on the true band margin.  Wall-clock is
dominated by a few heavy roots; boxes are the comparable cost.)

### 1.1 Cost vs margin

Boxes `≈ 2,358 + 190 · log₂(0.5 %/ε)` from `0.5 %` down to `0.1 %`, then flat (`+68` for another factor 5).  Depth never
moves off 13.  From the shipped `2 %` to `0.02 %` the band costs `1.86×` the boxes.  No blow-up as `ε → 0⁺`: the capture
function is piecewise constant, and shrinking `ε` only forces the few boxes that straddle the binding plateau's boundary one
or two levels deeper.

### 1.2 The checker is sharp (certified partial)

`+0.02 %` verifies and `−0.02 %` does not, so the band's true minimum is certified in `[μ₀/1.0002, μ₀]`:
`1.019911 ≤ min_{band} ≤ 1.020115` for the shipped cover.  The 21 uncertified boxes at `−0.02 %` all sit in
`cx ∈ [1.5, 1.5063], cy ∈ [0.5, 0.5063], θ ∈ [0, 0.53°]`: the plateau that runs from the wall tile germ `(1.5, 0.5, 0)` out to
the binding pose.

### 1.3 Why exactly zero margin fails: CHAIN slack at tile germs (`zeromargin.py diag`)

At `ε = 0` the 52 uncertified boxes sit at the same germ, `cx ∈ [1.5, 1.5023]`, `cy ∈ [0.5, 0.5023]`, `θ < 0.24°`.  `diag` on
boxes `[1.5, 1.5 + h]² × [0, θ]` with `h = 1/640, 1/1280, 3.9·10⁻⁴, 10⁻⁴` finds the sampled cover capture `≥ 1.025` in every box, while the best
CHAIN region is `0.999910863` in each.  **The shortfall is box-size independent**: a fixed region of the one-cut chain whose
provable witness weight is `8.9·10⁻⁵` below 1.  Any `ε > 0.0089 %` absorbs it.

The full-container run at `ε = 0.02 %` (`$S/full_e2.log`, killed at 2,500/6,400 roots) had 46 uncertified boxes in
`cx < 0.9`.  Their coordinates were never printed, because the run was killed.  `diag` on the candidate wall and corner germs finds exactly one
failing germ, **`(0.5, 1.5, θ → 0⁺)`** (box `cx ≥ 0.5, cy ≥ 1.5`).  The cover captures `≥ 1.0247` there,
but the best CHAIN region is `0.998376` at `ε = 0.02 %`, again box-size independent.  **That germ needs `ε ≥ 0.18 %`**
(`1.0002/0.998376 − 1`).  The mirror box `cy ≤ 1.5` passes with a region worth `1.0153`.  So `λ_c` at `m = 4` is set by
the one-cut/two-cut CHAIN disjunction at a wall tile germ, not by the cover.  A finer disjunction there would take
`λ_c` to `≈ 1.0001` (estimate; not attempted).

## 2. Full container at `ε = 0.25 %` (certified)

`c_e25.txt` = shipped cover × `40100000/40804617`, total `49479474611/3886154000 = 12.732248`.  It runs as 40
disjoint `cx` columns `[i/10, (i+1)/10]`, `i = 0..39`.  Their root boxes partition the 6,400 roots of one full run
(`roots()` skips exactly the columns outside `[cx-lo, cx-hi]`), and every column runs the same file with the same
options, so the union is a full verification:

```
for i in 0..39: zeromargin.py cert $S/c_e25.txt --depth 20 --nproc 4 --disj --chain-from 0 \
                    --cx-lo i/10 --cx-hi (i+1)/10 --chunksize 1          # $S/col.sh, $S/sched.py
```

| | value |
|---|---|
| columns verified | **40 / 40, every one `VERIFIED (PARTIAL …)`, 0 uncertified** |
| boxes | **29,372** over 6,400 roots (shipped `×1.05` run: 16,872, so `1.74×`) |
| leaves ADM / CHAIN / EMPTY | 4,113 / 9,786 / 3,987 (CORE, P1, MIX, TRI 0).  CHAIN carries 70 % of the non-empty leaves (shipped: 65 %) |
| max depth | 13 (limit 20 never reached) |
| CPU | 27.5 CPU-h (shipped: `≈ 13.6` CPU-h = 8 × 6138 s) |
| wall | 1 h 46 min (21:14 → 23:00) on 20 cores (4–5 columns in parallel; the heaviest column, `cx ∈ [1.5, 1.6]`, took 160 CPU-min) |

Columns `0–4` and `35–39` are 160 `EMPTY` roots each, since only `θ = 0` is admissible at `cx ≤ 1/2` or `cx ≥ 7/2`.  Most of the cost sits in
`cx ∈ [1.3, 2.7]`.  **Certified: `COVER^closed(4) ≤ 12.732248`**.  The weighted point set `$S/c_e25.txt`
(3,621 points) certifies `s(13) = 4` with total `12.732 < 12.956`.  It is not a certificate of record: it is not in
`certificates/`, it has no independent `zmcheck` run, and it is a union of 40 partial sweeps.  To promote it, run it as one
`zeromargin.py cert` and through `verify2/zmcheck`.

## 3. Validity factor (heuristic)

Measure: `s21_cover_eval.py hardscan` (319 angles, pitch 0.004, polish 1,400 seeds), honest cost = total / min.

| cover | LP / total | strict min | honest cost | notes |
|---|---|---|---|---|
| `r5` LP weights (shipped ÷ 1.05) | 12.339021 | **`≤ 0.9715385` (certified upper bound, exact pose)**; float `0.97154` | **12.7005** | `λ_v ≥ 1.0293` |
| `m4marginA` best (`closed4.py run --from-cert` r5 LP weights, lazy pool = item-3 61,656 + **near-tile 129,600**, `--deg-step 0.5`, colgen, 29 it / 2 h) | 12.336852 | 0.97189 at `(3.460, 2.536, 4.77°)` | 12.6937 | in-loop polish min wandered `0.96–0.99` every round |
| `m4marginA` it28 checkpoint | 12.337203 | 0.97073 at `(0.501, 2.500, 0.057°)` | 12.7092 | new hole at `θ = 0.001 rad` |
| `m4marginB` (`rung2_close.py` from it28, `--near-tile --fine-tilted --warm-lp`, 25 rounds / 1.5 h) | 12.344490 | **0.84542** at `(1.500, 3.478, 0.017°)` | 14.60 | in-loop stress min `0.986–0.993` |

Reading.  The near-tile rows cost nothing in LP (`12.339 → 12.337`), but the LP simply moves weight onto near-tile poses
that are still unsampled.  The holes sit at `θ ≈ k/D` rad (`0.057° = 10⁻³` rad) and below (`0.017°`), where one
grid point (`D = 1000`) enters or leaves a nudged wall tile.  A `0.02°` angle grid cannot see them, and `closed4.stress`
does not sample below `0.1°` except at `10⁻⁴, 10⁻³, 10⁻²` degrees.  The in-loop stress min (`0.993`) and the strict scan
(`0.845`) disagree by 15 %.  So **no LP-derived cover within `0.3 %` of valid was found at any LP total**, and the
`λ_v ≈ 1.029` of the shipped construction stands (heuristic).  The tilted interior deficits (`0.981–0.986` at
`30–40°`) remain as well, as in `S21_OVERHEAD.md` §2.

## 4. What this implies at `m = 5` (estimate)

    cert(5) = λ_c · (honest cost),   λ_c ≈ 1.0025   (m = 4, set by CHAIN at wall tile germs; per-germ, so m-independent)

| `s = 5` cover (`S21_COVER.md`) | honest cost (strict, heuristic) | × `λ_c = 1.0025` | vs 21 |
|---|---|---|---|
| `s5convR3` (best closing loop) | 21.016 | 21.07 | **above** |
| break-even honest cost | 20.948 | 21.000 | — |
| needed `λ_v` over converged LP 20.73 | `≤ 1.0105` | | `m = 4` achieved 1.029, `m = 5` 1.0138 |

Checker cost at `m = 5`, `λ_c = 1.0025`: `28.3k × (25/16) ≈ 45k` boxes, `≈ 60 CPU-h`, `≈ 3 h` on 20 cores (estimate, area scaling).
Affordable.  The obstacle is a cover with true min `≥ 0.9895` at LP `20.73`.  That is a validity question, and float
sampling does not answer it (§3).

## 5. Recommendation

1. Scale candidate covers by `1.0025 / (true min)`, not `1.05`.  Get the true min from `hardscan` + `zeromargin.py pose`, and
   confirm by the band-bisection of §1.2 (`±0.02 %` resolves it at 20 CPU-min a band).
2. Validity: drive the LP with the exact checker as oracle (`rung2_oracle.sh`, with the scaling of step 1).  Only it sees
   the `θ ≲ 1/D` near-tile holes.  Or add analytic near-tile rows: for each wall/interior tile, the poses where a grid
   point crosses an edge, `θ = (distance)/(lever arm)` for all points within `~0.01` of the tile boundary.  Float lattices at
   any fixed angle step reopen holes below it.
3. Optional: a finer CHAIN disjunction at wall tile germs would take `λ_c` from `1.0025` to `≈ 1.0001`.  Worth `0.05` at
   `m = 5`.  Only worth doing once validity is within `0.5 %`.
4. `m = 5` is **not reachable with current validity tooling** (heuristic best `21.07` certifiable).  It becomes marginally
   reachable if exact-oracle separation brings `λ_v` to `≤ 1.010` at LP `≤ 20.74`.

## 6. Reproduce

```
S=/tmp/claude-1000/-home-evand-math/f433c3ef-9956-443d-90fb-0c6d2ee917e2/scratchpad/m4margin
python3 search/s21_cover_eval.py hardscan certificates/rung2/s13_closed_cover_4.txt --nproc 20      # 1.0201154, 60 s
python3 search/zeromargin.py pose certificates/rung2/s13_closed_cover_4.txt \
        --cx 25199/50000 --cy 249609/100000 --u 2323208/581517279                                    # 40804617/40000000
python3 search/scale_cover.py certificates/rung2/s13_closed_cover_4.txt 40100000 40804617 $S/c_e25.txt # eps = 0.25 %
#   (e50: 40200000, e10: 40040000, e2: 40008000, e0: 40000000, em2: 39992000; all over 40804617)
taskset -c 0-19 python3 search/zeromargin.py cert $S/c_e2.txt --depth 24 --nproc 20 --disj --chain-from 0 \
        --cx-lo 1.4 --cx-hi 1.6 --chunksize 1                                                           # §1
python3 search/zeromargin.py diag $S/c_e2.txt --box 0.5,0.50078125,1.5,1.50078125,0,0.0000611          # §1.3: 0.998376
python3 $S/sched.py 25 4 <cols> ; python3 $S/sched_slot.py 25 4 <cols>                                  # §2 (col.sh)
python3 search/closed4.py run --s 4 --from-cert $S/lpw_r5.txt --tag m4marginA \
        --seed-rows $S/m4margin_seedrows.txt,$S/m4margin_neartile.txt --lazy-seed --deg-step 0.5 --nproc 8 --time 7200
python3 search/rung2_close.py $S/m4marginA_it28.txt m4marginB --near-tile --fine-tilted --warm-lp --nproc 4 --time 9000
python3 search/s21_cover_eval.py hardscan runs/closed4_m4marginA_best.txt runs/closed4_m4marginB_last.txt
```

Files (persistent copies, `runs/m4margin/`): `c_e25.txt` (the §2 cover, md5 `ecd72697c7fd8e70eb6e87188079a7de`), `c_e2.txt`, `c_e0.txt`,
`full25_logs/col_*.log` (the §2 run), `band_*.log` (§1), `full_e2.log`, `m4marginA.out`, `m4marginB.out`, `col.sh`, `sched*.py`,
`exactpose.py` (float pose → snapped admissible rational pose → exact capture).  Also `runs/closed4_m4marginA*`, `runs/closed4_m4marginB*`.
The near-tile pool (`$S/m4margin_neartile.txt`, 129,600 poses) is the 16 tile centres `± 0.012` (step 0.003) × `θ = 0.02..2°` (step 0.02),
clipped to admissible by `closed4.read_poses`.  The other `c_e*.txt` are `scale_cover.py` one-liners (§6).
