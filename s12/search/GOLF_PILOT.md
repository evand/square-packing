# Golf before Lean?  A measurement pilot on the s(21) mixed cover (2026-09-30)

Question: should the mixed-cover certificates (`certificates/s21/`, and by extension s(45), s(60)) be "golfed"
(re-shaped for a cheaper Lean tree) before the `S21CheckerCover` hypothesis is discharged in the kernel?
Code and data: `search/golf/` (new; nothing in `certificates/`, `search/zm_mixed.py`, `search/zeromargin.py`,
`verify2/` or `lean/Sqpack/` was touched).  All runs on physical cores 6–11.  Labels: **[measured]** = a run reported
as-is, **[model]** = the cost model of §1 applied to measured trees, **[proved]** = an exact checker's verdict (with
that checker's trust caveat).

## 0. Answer

> **Golf, but only the cheap way: scale the cover up to total ≈ 20.99 before generating the Lean tree.  Nothing
> else tried pays.**
>
> | cover | total | slack over the exact threshold | Lean tree (whole D4 region, estimated) | kernel CPU [model] low / **central** / high | part files at ≤ 15 GB |
> |---|---|---|---|---|---|
> | shipped `s21_mixed_cover_5.txt` | 20.894749 | 0.49 % | 633 k `Z` leaves, 436 point entries / leaf, 1.78 M L-lines, 1.13 G digits | 85 / **97** / 188 CPU-h | ~2,000 |
> | **A** = shipped × 20.99/20.894749 (`search/golf/covers/cand_A_scale.txt`) | 20.990000 | 0.95 % | 341 k `Z` leaves, 415 / leaf, 0.93 M L-lines, 0.60 G digits | 44 / **51** / 99 CPU-h | ~1,080 |
> | B = drop the 1,816 points of mass < 10⁻⁴, then scale to 20.99 (`cand_B_drop.txt`) | 20.990000 | 0.49 % | heavy roots only: 0.67 × current (A: 0.51 ×) | — | — |
>
> * **A halves the build** (paired, root for root: × 0.34–0.83, × 0.53 overall; sampling + model uncertainty
>   about ± 30 %), costs nothing mathematically (same geometry; total still `< 21`), and is checked: **zmx2 `--d4`
>   VERIFIED-D4** over all 2,500 roots [proved; 7 CPU-s], **zm_mixed.py `--d4 --cert-mode` VERIFIED on 512 sample
>   roots** (tile-germ block, the two hottest tilted cells, the two tilted-dip cells) [proved on the sample only],
>   exactly D4-invariant (`zmx2 d4`).  The Lean generator `gen_zmmtree.py` closes every sampled cell and root of both
>   covers (0 UNCERT).
> * **What the cost is.**  ~85 % of it sits in **8 of the 625 root cells** — the tilted cells around
>   `(1.2–1.5, 1.2–1.5)` and `(2.4, 1.2)` at `θ ≈ 30–53°` (exactly where zm_mixed needed SPLIT) — and it is two equal
>   halves: per-leaf overhead (number of leaves) and claimed points (≈ 430 per leaf, ~0.2–0.35 of the unit mass
>   collected from points of weight ~5·10⁻⁴), plus ~20 % L-lines.  Slack attacks the first (leaves × 0.54); the
>   second barely moves (points/leaf × 0.95).
> * **Geometric golf fails the exact threshold.**  Coarser segments (1/25 or 1/10 pieces) and merging points
>   (D4-orbit clusters on a 1/20 or 1/10 grid) raise the exact threshold above 21 (§4) — the cover's slack (0.49 %)
>   is far smaller than what such moves cost.  Dropping the smallest points costs 0.45 % of slack for 0.35 % of mass,
>   and is dominated by A.  An LP re-solve with a sparsity objective (fewer, heavier points) was **not** run (the
>   line-cover loop is ~4 h wall per stable run); its ceiling is bounded: point entries are ~40 % of the central cost,
>   so even halving points/leaf saves ≤ 20 %.
> * **Recommendation:** formalise **A** (or the cover scaled to `20.999`, another −5 % in zmx2 boxes), i.e. "scale
>   to the limit, then formalise as is".  Do not spend LP time on golf.  The remaining big levers are on the Lean side:
>   SPLIT (Lemma R) in the leaf, for the 8 hot cells; a cheaper per-leaf context (c_Z); and fewer digits per L-block
>   (they dominate the file sizes: ~1,800 digits per leaf vs s(32)'s 616).
> * **s(45), s(60):** the same holds qualitatively (same construction, same hot-cell geometry), with less to gain:
>   their shipped slack over the exact threshold is 0.71 % / 0.77 % and the room below `n` allows ≈ 1.22 % / 1.0 %, so
>   scaling should buy ≈ × 0.6–0.7 (s(45)) and ≈ × 0.8 (s(60)) by the zmx2-box law of §3 [heuristic; not run].

## 1. The cost model and its fit

Kernel CPU of the part files (the `Pts`/`Cov`/top files are minutes and ignored):

```
CPU = c_Z · (Z leaves) + c_adm · (ADM witness entries) + c_ch · (chain entries) + c_L · (L-lines) + c_T · (T-groups)
```

`E` leaves, `C` nodes and splits are a few `Nat` comparisons and are ignored.  Entry and line costs are the
profiled ones: `c_adm = 0.5 ms`, `c_ch = 2.5 ms` (heaviest s(13) chunk, `LADDER.md`), `c_L = 40 ms` per L-line,
`c_T ≈ 4 ms` (T1 one-leaf chunks, `notes/lean-segments.md` §4).  `c_Z` (decode, walls, candidate walk, per-leaf
context, region counts, file imports amortised) is fitted on the three full builds (`search/golf/costmodel.py`):

| build | measured kernel CPU | entries + lines at the profiled costs | fitted `c_Z` |
|---|---|---|---|
| S13 (`s(13) = 4`; 4,079 `Z`, 1.12 M ADM, 0.20 M chain entries) | 1,995 s | 1,070 s | 0.227 s |
| S32Z (`s(32) = 6`; 89,473 `Z`, 43.8 M ADM, 3.32 M chain) | 49,766 s | 30,200 s | 0.219 s |
| S16/T1 (segments only; 9,451 `Z`, ~24.8 k L-lines, 1,742 T) | 2,520 s | 998 s | 0.161 s |

Three parameter sets: **low** `c_Z = 0.15` (fits all three within 0.84–0.96 ×), **central** `c_Z = 0.22`
(S13 0.99 ×, S32Z 1.00 ×, T1 1.22 ×), **high** `c_Z = 0.46`, `c_adm = 1.0 ms`, `c_L = 60 ms`, calibrated on the
s(21) pilot's one-leaf chunks (0.69 s per T + 230-point leaf, 0.89 s per L + points leaf; those chunks had no `F`
pruning, so they likely overstate — the high set over-predicts the three full builds by 1.8–2.3 ×).  The truth for a
mixed s(21) leaf is between central and high until a real s(21) chunk is built; ratios between covers are much less
sensitive than absolute numbers (every term moved together in §2).

**Memory.**  RSS per part process ≈ 6.6 GB (Mathlib) + 15 KB × base-2²⁰ digits in the file (S13: 14.7 KB, S32Z:
15.0 KB per digit).  So memory is a layout choice: ≤ 15 GB per process ⇔ ≤ ~560 k digits per file.  The s(21) trees
carry ~1,800 digits per `Z` leaf (s(32): 616), mostly L-block data, hence ~1,100–2,000 part files (and ~3–5 GB of
generated Lean, s(32)'s 55 M digits being 267 MB).  At ≤ 30 GB per process the file count halves.

**Cheap proxies.**  (a) The Lean generator itself, run per root cell / per zm root (`search/golf/lean_census.py`, the
generator's search unchanged plus a tree walk): exact tree, features and digits, ~50 ms CPU per box.  (b) **The zmx2
`--d4` box count (7–35 CPU-s for the whole region)** tracks it: current → A is × 0.53 in zmx2 boxes and × 0.53 in
central Lean CPU; current → B × 0.61 vs × 0.67 on the heavy roots.  (c) The zm_mixed census is a poor proxy for the
*distribution* (the Lean tree has no SPLIT, so it is 2–3 × zm_mixed's boxes in the hot cells and ~1 × elsewhere),
but it is what the sample was stratified on.

## 2. The census: current cover vs candidates

**Sampling** (`search/golf/sample.py`, seed fixed).  The 625 root cells (pitch 1/10, all 8 u-bins, as
`gen_zmmtree.py`) are stratified by the shipped zm_mixed census aggregated per cell (`zm_census.py`): 8 heavy cells
hold 71 % of zm_mixed's 461 k boxes.  Strata S1–S5 (617 light cells) are sampled as whole cells (51 cells: S5 all 3,
S4 10/31, S3 10/55, S2 8/57, S1 20/471).  The 8 heavy cells (S6) are sampled as zm roots (1/20 × 1/20 × 1/32 in `u`,
512 roots) in five strata by zm boxes (18 roots: R1 6/377, R2 5/52, R3 4/44, R4 2/26, R5 1/13).  Expansion
estimator; candidates on the same units (paired).  A first attempt at whole heavy cells did not finish in 25 min
each (6 cells: 27–34 k boxes, 12–17 k `Z` leaves each, still running; `data/lc_current.jsonl`), which is why S6 is
sampled at root level.

| | current | A (scaled to 20.99) | A / current |
|---|---|---|---|
| `Z` leaves | 633 k (heavy cells 524 k) | 341 k (274 k) | 0.54 |
| point entries (ADM / chain) | 276 M (275.6 / 0.54) | 141 M (140.9 / 0.42) | 0.51 |
| point entries per `Z` leaf | 436 | 415 | 0.95 |
| L-lines / T-groups | 1.78 M / 109 k | 0.93 M / 68 k | 0.52 |
| digits → part files at ≤ 15 GB | 1.13 G → ~2,020 | 0.60 G → ~1,080 | 0.53 |
| generator CPU (Python) | 17.4 h | 8.9 h | 0.51 |
| **kernel CPU low / central / high** | 85 / **97** / 188 h | 44 / **51** / 99 h | 0.53 |
| of which the 8 heavy cells (central) | 82 h | 42 h | |
| sampling s.e. (central) | 18 h | 11 h | |

Central breakdown, current: leaves 39 h, ADM entries 38 h, L-lines 20 h, chains/T < 1 h.  The s.e. understates:
R5 (13 roots, 23 % of the heavy zm mass) rests on one root and R4 on two; a ratio estimator (Lean cost per zm box per
stratum) gives the same heavy total within 5 %.  Honest range for the current cover: **~70–190 CPU-h**, central ~100;
for A **~40–100 CPU-h**, central ~50 — consistent with the `TODO.md` estimate of 50–150 CPU-h.

Paired per root (central model, `python3 search/golf/estimate.py paired current A`): heavy roots × 0.34–0.59, light
cells × 0.44–0.83 (cells already closed at the root bin unchanged).  B (heavy roots only) × 0.37–1.14, overall
× 0.67: dropping the tiny points lowers points/leaf ~12 % but raises the leaf count in the tightest roots (it spends
slack); B / A = 1.01–3.2 per root.

## 3. Slack is the lever (zmx2 box law)

`search/golf/zmx2_curve.sh`: the shipped cover scaled to several totals, `zmx2 --d4` over all 2,500 roots
[proved at each total; `data/zmx2_curve.txt`].  The exact threshold is `20.7935` (`f* ∈ (0.99500, 0.99516]`,
`zbis/` bisection: the cover × 0.99500 is refused, × 0.9951562 verified).

| total | 20.80 | 20.82 | 20.85 | **20.8947** (shipped) | 20.92 | 20.95 | 20.97 | **20.99** (A) | 20.999 |
|---|---|---|---|---|---|---|---|---|---|
| slack over threshold | 0.03 % | 0.13 % | 0.27 % | 0.49 % | 0.61 % | 0.75 % | 0.85 % | 0.95 % | 0.99 % |
| zmx2 boxes | 6.42 M | 4.37 M | 2.79 M | 1.83 M | 1.50 M | 1.22 M | 1.08 M | 0.97 M | 0.92 M |

Roughly boxes ∝ slack^−0.7…−1.0 in this range; the Lean tree follows (§1 (b)).

## 4. Geometric golf: the variants and their exact thresholds

`search/golf/golf.py` (each output exactly D4-invariant, masses rounded up), `make_covers.sh`, `bisect_all.sh`;
threshold = smallest factor at which `zmx2 --d4` verifies, × total [proved upper end, refusal at the lower end]:

| variant | points / segments | threshold total | usable (< 21)? |
|---|---|---|---|
| shipped | 7,536 / 1,872 | 20.7935 | yes (0.99 % room) |
| segments merged 2 → 1 (1/25 pieces) | 7,536 / 944 | 21.078 | **no** |
| segments merged 5 → 1 (1/10 pieces) | 7,536 / 376 | 21.430 | **no** |
| points merged by D4 orbit, 1/20 grid | 1,492 / 1,872 | 21.123 | **no** |
| points merged by D4 orbit, 1/10 grid | 652 / 1,872 | 21.535 | **no** |
| points of mass < 10⁻⁴ dropped (1,816 points, mass 0.073) | 5,720 / 1,872 | 20.887 | yes (0.49 % room) → B |

The line densities are not smooth at the 1/50 scale (peaks next to the corner squares) and the points are a dense
LP residue whose positions matter at the 0.1 % level, so any move larger than the slack fails.  A real reduction of
points per leaf needs a new LP (point columns on a coarser lattice, or an L0-style reweighting loop) — not attempted;
by §2 it is worth at most ~20 % (points are ~40 % of central cost) against a 4+ h LP loop and a new exact check.

## 5. Candidate A in full

`search/golf/covers/cand_A_scale.txt`: mixed v1, `s = 5`, `D = 1000`, `W = 10¹¹`, 7,536 points + 1,872 segments
(the shipped geometry), every mass × `20990000000/20894749197` rounded up, total `524750000723/25000000000 =
20.990000029 < 21`.

* `zmx2 d4`: exactly D4-invariant.  `zmx2 cert --d4`: **VERIFIED-D4**, 2,500 roots, 968,422 boxes (shipped:
  1,826,222), max depth 24, 0 uncertified, 7 CPU-s [proved].
* `zm_mixed.py cert --d4 --cert-mode` with the shipped settings on 512 roots [proved on these roots; `zm/`]:

  | region | roots | A: boxes / CPU | shipped cover (from its records): boxes / CPU |
  |---|---|---|---|
  | tile germs `[1.4,1.6]²` | 256 | 2,358 / 258 s | 3,614 / 438 s |
  | tilted `[2.4,2.5]×[1.2,1.3]` | 64 | 9,600 / 917 s | 19,502 / 1,979 s |
  | tilted `[1.35,1.45]×[1.2,1.3]` | 64 | 25,544 / 2,150 s | 45,824 / 3,937 s |
  | dip 9.1° `[0.5,0.6]×[1.4,1.5]` | 64 | 1,412 / 114 s | 2,350 / 218 s |
  | dip 31.9° `[2.3,2.4]×[1.6,1.7]` | 64 | 1,206 / 215 s | 1,816 / 401 s |

  All VERIFIED-D4 (PARTIAL), 0 uncertified; × 0.53 boxes, × 0.52 CPU.  A full `--d4` run (≈ 6–7 CPU-h by these
  ratios) is needed before A replaces the shipped file in a bundle.
* Lean generator: closes every sampled unit (0 UNCERT), §2.

## 6. Files

| file | what |
|---|---|
| `search/golf/golf.py` | transforms: `info`, `scale`, `merge`, `coarsen`, `drop` (exact D4 check on output) |
| `search/golf/lean_census.py` | per-cell / per-root census of the Lean tree (`gen_zmmtree.py`'s search, unchanged) |
| `search/golf/costmodel.py`, `estimate.py` | the model and its fit; stratified extrapolation, `paired` comparison |
| `search/golf/zm_census.py`, `sample.py` | zm_mixed census per 1/10 cell; the stratified sample (`data/strata*.json`, `data/sample_*.txt`) |
| `search/golf/zbisect.sh`, `bisect_all.sh`, `zmx2_curve.sh`, `make_covers.sh`, `run_census.sh`, `zm_sample.sh` | the runs |
| `search/golf/covers/cand_A_scale.txt`, `cand_B_drop.txt` | the candidates (the unscaled variants `v_*` are regenerated by `make_covers.sh`) |
| `search/golf/data/lc_{current,A,B}_{cells,roots}.jsonl` | census records (`lc_current.jsonl`: the six 25-min heavy-cell timeouts) |
| `search/golf/zm/`, `search/golf/zbis/*.bisect`, `data/zmx2_curve.txt` | zm_mixed sample records; bisections; the slack curve |

Pilot cost: ~7.8 CPU-h (Lean-generator census 6.7 h, zm_mixed sample 1.0 h, zmx2 minutes).

## 7. Reproduce (from `s12/`)

```
python3 search/golf/zm_census.py certificates/s21/zm_mixed_d4/roots.jsonl > search/golf/data/zm_cells_current.json
python3 search/golf/sample.py
python3 search/golf/golf.py scale certificates/s21/s21_mixed_cover_5.txt search/golf/covers/cand_A_scale.txt --total 20.99
bash search/golf/run_census.sh current certificates/s21/s21_mixed_cover_5.txt     # RT=3600 for root R4_1
bash search/golf/run_census.sh A search/golf/covers/cand_A_scale.txt
python3 search/golf/costmodel.py; python3 search/golf/estimate.py current A; python3 search/golf/estimate.py paired current A
bash search/golf/make_covers.sh; bash search/golf/bisect_all.sh v_seg2 v_seg5 v_merge05 v_merge10 v_drop1e4
bash search/golf/zmx2_curve.sh; bash search/golf/zm_sample.sh A search/golf/covers/cand_A_scale.txt
```
