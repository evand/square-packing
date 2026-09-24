# Where the `m = 4` certificate's overhead sits, and what it predicts at `m = 5` (task s21-overhead; 2026-09-22)

> **Correction (2026-09-23, `M4_MARGIN.md`):** the shipped cover's true minimum is `<= 40804617/40000000 = 1.0201154` (exact pose (0.50398, 2.49609, 0.4578°)), not `1.02766` — that was one exact spot check; the float `0.9715 × 1.05` was right.  So the LP weights' shortfall is `>= 2.9 %`, and the checker needs only `~0.25 %` (not `1–2 %`).  Numbers below citing `1.0277` are superseded.

Code: `search/s21_overhead.py regions` (new, exact `Fraction` sums: per-cell and per-class weights
of the covers and of the `nu_f` measure).  Float scans used `closed4.scan_angle` / `closed4.polish`
(imported, unmodified) from scratch scripts.  The exact capture checks used `zeromargin.py pose`, and the
exact checker used `zeromargin.py cert --cx-lo/--cx-hi` bands.  Labels: **certified** (exact
`Fraction`, a theorem), **heuristic** (float LP / float scan / local search), **estimate** (a
model extrapolation).

## 0. Answer

> At `m = 4` the certificate costs `12.9560 − 12.2688 = 0.687` over `nu_f`.  That is
> `0.070` from the LP not being `nu_f`, about `0.3` to make the LP weights valid at all, and about
> `0.3` of margin for the exact checker.  The last two come from a single **uniform** factor
> `×1.05`, so they scale with the total weight.  They are **interior-type, `O(m²)`**: `0.043` per
> unit area over `nu_f`, `0.039` over the LP.  Nothing in RUNG2.md's list of "what it took" cost
> weight except the rows and the scaling.  `ADM`, `clip_bin`, `CHAIN` and the suffix down-set are
> checker-side and cost 0 weight.
>
> **With the `m = 4` recipe, `m = 5` costs about `1.05 × (nu_f(5) + 0.11)`, which is 21.4 to 22.1 for
> `nu_f(5) ∈ [20.3, 20.9]` (estimate). That is not below 21.**  A cover below 21 needs the combined
> factor (validity × checker margin) cut from `1.05` to at most `21/(nu_f(5)+0.11)`:
> `1.029` at `nu_f = 20.3`, `1.019` at `20.5`, `1.009` at `20.7`, and impossible at `≥ 20.89`.
> At `m = 4` the exact checker verifies the hardest band (`cx ∈ [1.4, 1.6]`, which holds the binding
> pose and an interior tile pose) at a total factor of **`1.035`** over the LP weights (§4, certified
> partial sweep, `1.5×` the boxes).  The LP weights themselves were invalid by `2.2–2.9 %` (§2).  So of the 5 %,
> 2.2–2.9 % was row under-sampling and only `0.6–1.3 %` was checker margin.  The rest was a speed choice.  None
> of these is a theorem.  Realistic floor after full separation: a factor of
> `1.01–1.02`.  **So `W < 21` is plausible only if `nu_f(5) ≲ 20.5`.  It is marginal at `20.5–20.7`
> and out of reach at `≳ 20.8` with this method (estimate).**

## 1. The `m = 4` ledger

| quantity | value | status |
|---|---|---|
| `nu_f^closed(4)` | `≥ 12.2688039` (`COVER4.md`); float LP converged at `12.268806` | certified / heuristic |
| LP of the run of record (`closed4.py run --tag r5`, sampled rows + colgen) | `12.338773`; the weights are exactly `certified × 20/21` (`12.339021`) | heuristic |
| min capture of those LP weights | **`≤ 0.978724`** exact (`1.02766 / 1.05`, admissible rational pose, §2); float polish finds `0.9715` at the wall edge of admissibility | certified upper bound / heuristic |
| smallest valid scaling `λ_v` of the LP weights | `1.0217 – 1.029` → honest cost `12.61 – 12.70` | heuristic |
| scaling actually shipped | `×21/20` → `2591194431/200000000 = 12.955972` | certified |
| min capture of the shipped cover | `≤ 1.02766` (exact pose below); `closed4.py stress` `1.0314` | certified upper bound / heuristic |
| checker slack as shipped, `1.05/λ_v` | `1.020 – 1.028` | heuristic |
| checker slack actually needed on the hardest band | `≤ 1.035/λ_v = 1.006 – 1.013` (`×1.035` band verifies, §4) | certified partial sweep + heuristic `λ_v` |

Decomposition of `0.687` over `nu_f`:

| piece | weight | what it is | type |
|---|---|---|---|
| LP − `nu_f` | `0.070` (`0.0044`/unit area) | column/row discretisation of the LP; includes `≤ 0.035` for the item-3 rows (`RUNG2.md` §10: `12.3035 → 12.3388` with colgen) | spread across the container, `O(m²)` |
| validity (`λ_v − 1`) · LP | `0.27 – 0.36` | row under-sampling: poses the LP never saw (§2) | multiplicative, `O(m²)` |
| checker slack (`1.05/λ_v − 1`) · `λ_v`·LP | `0.25 – 0.34` | speed margin for the exact checker (`×1.03` "only 0.5 % margin, verifies far more slowly", `RUNG2.md` §10) | multiplicative, `O(m²)` |

Exact spot check (certified): `zeromargin.py pose certificates/rung2/s13_closed_cover_4.txt --cx
1436496011/956220587 --cy 256600555/509421602 --u 3390122/910485467` gives `51382989/50000000 = 1.02765978`, admissible.
That pose is `(1.50226, 0.50371, 0.427°)`: the edge tile `(1,0)`, nudged.

## 2. Where the LP weights are short (heuristic, float scan + polish)

The LP weights are the shipped cover divided by 1.05.  Scan: pitch `0.004`, `θ` from `0.02°` to `53°` in
`0.02–0.25°` steps.  7,540 seeds `< 1` were polished and binned by the centre's cell class, canonical under D4:

| centre class | `θ < 3°` (near-tile) | tilted |
|---|---|---|
| corner cell | `≥ 1.00001` (no deficit) | `≥ 1.0064` |
| edge cell (wall band) | **`0.9715`** at `(1.503, 0.503, 0.32°)` | `0.9817` at `(0.703, 1.332, 37.5°)` |
| interior cell | `0.9926` at `(1.505, 1.508, 1.07°)` | `0.9828` at `(1.343, 1.690, 40.5°)` |

The binding deficit is **wall-type**: the edge tile nudged by a fraction of a degree, the transpose of the
`ZEROMARGIN.md` item-3 family at the wall.  But interior tilted poses are only `1 %` behind it.  A
repair confined to the walls would still leave `1.7–1.8 %` to pay everywhere else.  So the deficit is
effectively uniform, and a uniform factor is the right model for it.  The shortfall is a sampling artefact of that
run, not a property of the problem.  `FAMILY.md` §2's closing loop over fixed columns reached a stress min of
`0.9958` (`0.4 %`).

Tiling-like poses (`θ → 0`) are LP-tight: capture `1.00001` in all three classes.  These are the
zero-margin configurations that `CHAIN` has to certify, and after `×1.05` each carries exactly the 5 % pad.

## 3. By region (certified: exact sums, `s21_overhead.py regions`)

Per-cell weight.  A point on an interior grid line is split between its cells.  Classes: 4 corner, `4(m−2)` edge,
`(m−2)²` interior.

| object | corner (mean) | edge (mean) | interior (mean) | total |
|---|---|---|---|---|
| shipped cover (`×1.05`) | 0.5288 | 0.8098 | 1.0907 | 12.9560 |
| LP weights (`÷1.05`) | 0.5036 | 0.7712 | 1.0388 | 12.3390 |
| `runs/closed4_best.txt` (LP 12.4175) | 0.4977 | 0.7775 | 1.0517 | 12.4175 |
| `nu_f` measure, pose centres | 0.8132 | 0.6801 | 0.8937 | 12.2688 |

* The `×1.05` pad by class: corner `0.101`, edge `0.308`, interior `0.208` (of `0.617`).  At `m = 4` half
  of it sits in the wall band because 8 of the 16 cells are edge cells.  The pad follows the weight, so its
  share is `O(m²)` for large `m`.
* The two LP covers (`r5` and `closed4_best`) agree per class to within `0.013`, so the per-class structure is stable.
* 91 % of the shipped weight is on the interior grid lines (`11.77` of `12.96`).  The tiling structure *is*
  the cover.
* Cover and measure per cell are not comparable locally (duality matches totals only).  The measure puts
  its corner mass on the 4 corner squares (`0.813` each).  Its tilted mass sits in edge and interior cells
  (`2.56 + 3.17`).  70 % of its mass is at `θ < 5°` (`COVER4.md`).  So a margin required *only* on
  near-tile poses would cost, to first order by LP sensitivity, `≈ 0.7 δ · nu_f` rather than `δ · nu_f`:
  a 30 % saving at best.

**Cellular extrapolation to `m = 5`** (`4a + 12b + 9c`, estimate):

| object | `m = 5` extrapolation |
|---|---|
| `closed4_best` classes | **20.79** (matches the heuristic `s = 5` LP of `FAMILY.md` §3.1, `20.75–20.9`; the current `s5conv` runs are at `20.9–21.0` at iteration 4, unconverged) |
| shipped-cover classes | **21.65** |
| `nu_f` measure classes | 19.46.  **Not reliable**: at `m = 4` the "interior" cells all touch the wall ring, so the measure's tilted middle ring does not extrapolate by cells |

**Stretch test (heuristic).**  Insert a copy of a middle strip into the shipped cover.  The resulting
`[0,5]²` set weighs `21.525` (two variants).  Float scan: every `θ < 1°` pose, including the new interior cells,
still captures `≥ 1.03`, so the tiling-like part is **cellular** and transplants unchanged.  Tilted interior poses
drop to `0.850`, so the tilted part does not.  At `m = 5` the tilted weights must be re-solved, which is what the LP does.
The zero-margin part of the certificate, the part `CHAIN` has to handle, has the same per-cell
structure at every `m`.

## 4. Classification of what `RUNG2.md` needed

| item | weight cost at `m = 4` | why needed | scaling |
|---|---|---|---|
| `ADM` (off-centre wall witnesses, Lemma A) | 0 | the `O(u)` → `O(u²)` wall tolerance | checker; wall-type compute |
| `clip_bin` | 0 | wall boxes have lower-dimensional admissible bins | checker; wall-type compute |
| `CHAIN` (one chain at wall tiles, product of two at interior tiles) | 0 directly; needs slack (below) | Theorem 1: `m²` disjoint PIN germs, one per tile | slack paid at `m²` tile germs, `O(m²)` |
| suffix down-set | 0 | completeness | checker |
| item-3 rows (both edges on grid lines, sliding) | `≤ 0.035` LP | `closed4.py`'s row lattice steps over them | one family per grid-line pair: `O(m)` families; `≲ 0.06` at `m = 5` (estimate) |
| global `×1.05` | `0.617` | validity (2.2–2.9 %) + checker slack (2.0–2.8 %) | **`O(m²)`**: 5 % of the total |

No piece of the overhead is corner-type.  Corner cells are the only class with no LP deficit (§2).

**Checker slack, measured on one band (certified partial sweeps, `cx ∈ [1.4, 1.6]`, reduced domain, depth 18, `--disj --chain-from 0`).**
This band holds the binding wall-type deficit pose `(1.502, 0.504, 0.4°)` and the interior tile pose `(1.5, 1.5, 0)`.

| cover | total | exact min capture near the band's worst pose | boxes | max depth | leaves ADM / CHAIN / EMPTY | uncertified | wall (2 procs) |
|---|---|---|---|---|---|---|---|
| shipped `×1.05` | 12.955972 | `1.0277` | 1,508 | 13 | 85 / 734 / 95 | **0** | 3,067 s |
| `×1.035` over the LP weights (`scale_cover.py … 69 70`) | 12.770887 | `≈ 1.013` (`0.9787 × 1.035`) | 2,324 | 13 | 120 / 1,095 / 107 | **0** | 4,715 s |

Both runs are `VERIFIED (PARTIAL)`.  Cutting the pad from 5 % to 3.5 % costs `1.54×` boxes on this band and does not
increase the depth.  The checker's own slack is therefore at most `1.035/λ_v = 1.006–1.013`.  The `2 %` shipped
with the cover was a speed choice, not a requirement.  (Not a full verification of the `×1.035` file: the other 38 columns were not run.  As a by-product it
certifies that the LP weights capture `≥ 1/1.035 = 0.9662` at every admissible pose with `cx ∈ [1.4, 1.6]`.
That is consistent with both the exact `0.9787` and the float `0.9715`, so it does not decide `λ_v`.)

## 5. Scaling model and prediction

    cert(m) = Λ · LP(m),   LP(m) = nu_f(m) + g(m),   g(4) = 0.070,   Λ = λ_v · λ_c

All three terms are proportional to the area.  `g` is spread discretisation.  `λ_v` is uniform because
the deficits at walls and in the interior are within `1 %` of each other (§2).  `λ_c` is the margin at the `m²`
tile germs and the `O(m²)` measure of grid-line families (Theorem 1; §3 stretch test).  So the overhead is
`≈ (Λ − 1)·nu_f(m) + g(m)`, which is `0.043 m²` at the `m = 4` settings.  There is **no second data point**:
no zero-margin certificate exists at any other `m`.  The closest proxies are `CLOSED4.md` at `s = 3.99`
(LP `12.12`, stress cost `12.70`: a `4.8 %` validity gap from the same row under-sampling) and `FAMILY.md`
§2 (`0.6 %` after 67 separation rounds).  Both put `λ_v` in the same `0.5–5 %` range.

| `nu_f(5)` | `Λ = 1.05` (m=4 recipe) | `Λ = 1.03` | `Λ = 1.02` | `Λ = 1.01` | break-even `Λ*` |
|---|---|---|---|---|---|
| 20.3 | 21.43 | 21.02 | 20.82 | 20.62 | 1.029 |
| 20.5 | 21.64 | 21.23 | 21.02 | 20.82 | 1.019 |
| 20.7 | 21.85 | 21.44 | 21.23 | 21.02 | 1.009 |
| 20.9 | 22.06 | 21.64 | 21.43 | 21.22 | < 1 |

(estimate: `g(5) = 0.11` by area scaling; ±0.05.)

`W < 21` at `m = 5` needs `Λ ≤ 1.02` if `nu_f(5) = 20.5`.  That is two of the three `m = 4` overhead components
cut by a factor of 2.5, and it needs:
1. **`λ_v → 1.005`.**  Run the separation to convergence with colgen *and* the exact checker as the oracle
   (`rung2_oracle.sh`), with near-tile rows (`θ ∈ [0.05°, 1°]`, centre within `0.01` of a tile pose, at the
   wall bands) seeded explicitly.  That is where the `m = 4` LP was blind.
2. **`λ_c → 1.01`.**  §4 shows the checker already tolerates `1.006–1.013` on the hardest `m = 4` band,
   at `1.5×` the boxes and the same depth.  At `m = 5` expect `(25/16)` × that per-area cost, and more as the
   margin shrinks further.  The run of record took 1 h 42 min on 8 processes at `×1.05`, so the check itself
   is affordable.  Untested: whether `~0.5 %` still terminates.
3. Targeted rather than uniform padding saves `≤ 30 %` of the pad (§3), so it cannot replace (1)–(2).

## 6. Reproduce

```
python3 search/s21_overhead.py regions                                        # §3 tables, exact, 2 s
python3 search/zeromargin.py pose certificates/rung2/s13_closed_cover_4.txt \
   --cx 1436496011/956220587 --cy 256600555/509421602 --u 3390122/910485467  # 1.02765978, §1
python3 search/scale_cover.py certificates/rung2/s13_closed_cover_4.txt 69 70 OUT_x1035.txt
python3 search/zeromargin.py cert OUT_x1035.txt --depth 18 --nproc 2 --disj --chain-from 0 --cx-lo 1.4 --cx-hi 1.6
```
