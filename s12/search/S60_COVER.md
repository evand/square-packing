# `s(60) = 8`?  The `m = 8` line-density cover (task s60-cover, 2026-09-28)

The `k = 8` rung of the `s(k² − 4) = k` family, by the `s(45)` route (`S45_COVER.md` §10): a line-density (mixed v1)
cover of `[0,8]²` of total `< 60` in which every closed unit square captures mass `≥ 1`.  Labels: **[proved]** = an
exact checker run as named (with that checker's trust caveat), **[measured]** = float LP / float scan / exact
arithmetic on a finite object reported as-is, **[heuristic]** = reading or extrapolation.

Code (new): `search/s60_prep.py` (`transplant`: grow an odd-`m0` D4 mixed cover to `[0,M]²` by replicating the centre
cell; `ring` / `ring-table`: the exact savings profile of §4), `search/m8_zbisect.sh` (`m7_zbisect.sh` with the cores
as parameters).  Unchanged: `search/line_cover.py`, and the pinned checkers `search/zm_mixed.py` (`ee3e2915…`),
`mixed_cover.py` (`bb89de15…`), `zeromargin.py` (`640fe453…`), `verify2/target/release/zmx2`.  All runs on physical
cores 0–13; peak machine memory in use `≈ 50 GB`.

## 0. Answer

> **GO.**  `runs/s60_mixed_candidate_8.txt` (mixed v1, `s = 8`, 23,744 points + 5,216 grid-line segments, exactly D4
> invariant) has total **`748233441/12500000 = 59.85867528 < 60`**, and every closed unit square in `[0,8]²` captures
> mass `≥ 1`:
> * **`zmx2 --d4`: `VERIFIED-D4`** (6,400 roots, 2,617,534 boxes, 0 uncertified, 23 CPU-s) and **`zmx2 --full` (no
>   symmetry): `VERIFIED`** (51,200 roots, 21,036,120 boxes, 0 uncertified, 182 CPU-s) [proved, modulo zmx2's
>   outward-rounded binary64 chord-end enclosure, `ZMX2.md` §5];
> * **`zm_mixed.py --d4 --cert-mode`** (pinned sha `ee3e2915…`, the `s(21)`/`s(45)` shipped settings): **`VERIFIED-D4`**
>   — 102,400 roots, 500,134 boxes, max depth 18, 0 uncertified, **19.6 CPU-h** (107 min on 12 processes) [proved:
>   exact `Fraction`/integer arithmetic; the program is not formally verified].
>
> With the reduction of `certificates/s21/FORMAT.md` this gives `s(60) ≥ 8`, hence **`s(60) = 8`** (the grid gives
> `≤ 8`), at the same trust level as the `s(45)` bundle.  Not bundled, not committed (brief); no Lean.
>
> Route: the `m = 7` line-cover recipe at `s = 8`, seeded by the `s(45)` LP cover with its centre cell duplicated,
> **plus a coarse uniform row family** (without it the first LP collapses onto the lines, §1).  Chain B round 4: LP
> `58.3984`, exact threshold `f* ≤ 1.0171875` (zmx2 bisection) → exact cost `59.402`; candidate `= B4 × 41/40`
> (`0.77 %` over `f*`, `0.24 %` under 60).  ≈ 10 h wall from launch to verdict (LP 8 h, zm_mixed 1.8 h).
>
> **Savings (§4):** `64 − 58.398 = 5.60` at the LP, **`4.60` at the exact threshold**, `4.14` shipped.  Across
> `k = 5, 7, 8` the achievable saving (round cover × `f*`) is `4.21, 4.54, 4.60`: above 4 with room, not shrinking.
> It is entirely a wall-ring deficit (`+4.36 → +5.99` at the LP, `0.27 → 0.21` per unit of wall), paid back partly by
> the interior, whose excess at `f*` grows with the area (`−0.13, −0.70, −1.01`).

## 1. Runs

**Seed.**  `s60_prep.py transplant runs/lc_m7A_r1.txt 8 runs/s8_tp_m7A1.txt`: the `m = 7` A1 round cover (the one the
`s(45)` certificate is scaled from) with its centre cell column/row duplicated (`[3,4] → [3,4] ∪ [4,5]`), 30,268 points
+ 5,248 segment pieces, total `58.1657` (the "cell-model" value).  Exactly D4-invariant (checked); the identity
transplant `7 → 7` reproduces the input exactly.  As a cover it is poor: `zmx2 --d4` on the transplant of the *shipped*
`s(45)` cover (total `59.47`) finds dips of float value `0.924` at `(3.70, 3.70, θ ≈ 41°)`, at the new centre seam
[measured].  It serves as the column set and warm-row source.

**Chain A** (`runs/lc_m8A.sh`, the `m = 7` command verbatim at `s = 8`: `--q 50 --no-line-points --germ`, point
columns lattice 0.05 + transplant support = 6,209 point orbits / 49,332 atoms, 800 segment orbits / 5,600 pieces;
warm rows = 800 worst poses per angle of the transplant at pitch 0.01; germ family `q = 50`: 168,402 poses → 68,888
permanent rows).  **A0 collapsed onto the lines**: LP `56.16`, segments `55.4`, 1,160 points, lattice minimum
**`0.010`** (tilted squares around `(1.35, 1.40, 35–53°)` capture almost nothing).  The warm rows were all at the
transplant's seams, so most of the container had no rows at mid angles.  It recovered only slowly (A1–A4: `56.59 /
0.72`, `57.05 / 0.82`, `57.51 / 0.85`, `57.80 / 0.89`) and was stopped after A4.

**Chain B** (`runs/lc_m8B.sh`): as A plus a permanent **uniform row family** `runs/s8_uniform_rows.txt` (centres
pitch 0.1 in `[0,4]²` × `θ ∈ {2, 5, 10, …, 45}°`; 12k new rows after D4 dedup).  That was the fix: B0 was already
at lattice minimum `0.87`.

| round | LP [measured] | segment mass | rows | float lattice / polished / near-tile / item-3 min [measured] | IPM solve | `f*` (zmx2 bisection) | `total · f*` |
|---|---|---|---|---|---|---|---|
| B0 | 57.794301 | 44.35 | 181,881 | 0.87014 / 0.85888 / 0.98977 / 0.98855 | 5,515 s | — | — |
| B1 | 58.098389 | 42.72 | 119,564 | 0.94858 / 0.94690 / 0.99645 / 0.99056 | 3,742 s | — | — |
| B2 | 58.281521 | 40.99 | 134,617 | 0.96862 / 0.96576 / 0.99668 / 0.99496 | 4,590 s | `(1.035625, 1.036250]` | 60.395 |
| B3 | 58.363114 | 40.51 | 131,359 | 0.97686 / 0.97446 / 0.99858 / 0.99779 | 4,758 s | `(1.0265625, 1.0268750]` | **59.932** |
| **B4** | **58.398382** (exported 58.398548) | 40.37 | 134,906 | 0.98694 / 0.98436 / 0.99938 / 0.99899 | 5,439 s | **`(1.016875, 1.0171875]`** | **59.402** |

(B5 was stopped in its LP solve once the B4 candidate had verified.)

`f*` = the smallest factor at which `zmx2 --d4` verifies the scaled round cover (`search/m8_zbisect.sh`); the upper end
is [proved] (zmx2), the lower end a refusal [measured] (zmx2 loses `< 0.1 %`, `ZMX2.md` §0).  Each B round is
`≈ 75–95 min`, dominated by the single-threaded HiGHS IPM (`3.7–5.5 ks` at 120–180k rows × 7–7.7k orbit columns; `m = 7`:
`2.0–2.9 ks` at 88–114k rows).  B3 was already exactly below 60 (`59.932`, `0.11 %` room); B4 is the candidate round.

**Where the last dips are** (zmx2 UNCERT boxes, float value, in the scaled cover): B3 at `(3.152, 1.548, 46.1°)` and
its reflection `(1.548, 3.152, 43.9°)` — interior, tilted near 45°, one cell in from the wall ring.

## 2. The candidate

`runs/s60_mixed_candidate_8.txt` = `lc_m8B_r4 × 41/40` (masses rounded up; `line_cover.py scale`), mixed v1, `s = 8`,
23,744 points + 5,216 segments, total **`748233441/12500000 = 59.85867528 < 60`**; room `0.1413` (`0.24 %`), `0.77 %`
over B4's exact threshold `f* ≤ 1.0171875` (`s(45)`: `0.71 %`).  sha256 `2d0e456e…3b3b12a41`.
`zmx2 d4`: exact D4 invariance; `mixed_cover.py`: well formed.

| check | result |
|---|---|
| `zmx2 cert --d4` (7 threads) | **`VERIFIED-D4`** — 6,400 roots, 2,617,534 boxes, 0 uncertified, max depth 24, 23 CPU-s [proved, modulo zmx2's binary64 chord-end enclosure, `ZMX2.md` §5] |
| `zmx2 cert --full` (no symmetry) | **`VERIFIED`** — 51,200 roots, 21,036,120 boxes, 0 uncertified, max depth 24, 182 CPU-s [proved, same caveat] |
| `zm_mixed.py --d4 --cert-mode` (s(21)/s(45) settings: `--disj --chain-from 0 --depth 24 --pitch 1/20 --ubins 16`), 102,400 roots | **`VERIFIED-D4`** — 102,400 roots, 500,134 boxes, max depth 18, 0 uncertified; leaves ADM 132,043 / CHAIN 72,081 / SPLIT 59,067 / PIECE 8,666 / EMPTY 29,410 (Lemma T raised the piece bound at 51,244 leaves, Lemma L at 219,739); 70,609 CPU-s = 19.6 CPU-h, 6,435 s wall on 12 processes; manifest header shas = the pinned files, input `2d0e456e…` [proved] |

## 3. Reproduce

```
python3 search/s60_prep.py transplant runs/lc_m7A_r1.txt 8 runs/s8_tp_m7A1.txt
python3 -c "…"                     # runs/s8_uniform_rows.txt: see its header (pitch 0.1 in [0,4]^2 x 10 angles)
bash runs/lc_m8B.sh                # stop after "[m8B] A4" (~7 h; LP 4-5.5 ks per solve)
ZCORES=0-6 ZTHREADS=7 bash search/m8_zbisect.sh runs/lc_m8B_r4.txt 1.005 1.025 6    # f* in (1.016875, 1.0171875]
python3 search/line_cover.py scale runs/lc_m8B_r4.txt runs/s60_mixed_candidate_8.txt --factor 1.025
verify2/target/release/zmx2 cert runs/s60_mixed_candidate_8.txt --d4   --threads 7    # VERIFIED-D4
verify2/target/release/zmx2 cert runs/s60_mixed_candidate_8.txt --full --threads 7    # VERIFIED
bash runs/zm_mixed_s60.sh          # zm_mixed.py --d4 --cert-mode, 102,400 roots
python3 search/s60_prep.py ring-table <cover[:factor]> ...                             # §4
```

## 4. The savings structure, `k = 5..8`

**Convention.**  For the axis box `B_w = [w, k−w]²`, `μ_½(B_w)` counts mass strictly inside fully, mass on an edge
of `B_w` (points on it; grid-line segments lying along it) at `½`, a point at a corner of `B_w` at `¼`.  Then
`def(w) = area(B_w) − μ_½(B_w)` and the ring deficit `def(w) − def(w+1)` split the grid-line mass on `x, y = w+1`
half/half between the two rings, so the ring deficits add up exactly to the saving `k² − total` [exact `Fraction`,
`s60_prep.py`].  (With closed boxes, §2.3 of `CORNER_DEFICIT.md` gives `def ≤ 0` for every integer `w ≥ 1`; with
the half convention the interior boxes are also `≤ 0` in every cover here, and the whole saving sits in the wall
ring.)  A cover scaled by `F` has ring deficits `area − F·μ`.  Half-integer `w` (`ring` mode) is dominated by the
line mass it excludes/includes (e.g. `def(½) ≈ −8.9` at `k = 8`) and says nothing more; not tabulated.

Rows: the LP round cover (`F = 1`, what the LP "wants"), the same at its exact threshold `F = f*` (what is actually
achievable with that cover, zmx2), and the shipped/candidate file.  `k = 6` has no line-density cover; its row is the
shipped point cover (`s(32)`), which is not comparable in overhead.  `runs/s60_ring_table.txt`.

| k | cover | F | F·total | **saving** `k²−F·total` | wall ring `[0,1]` | per unit wall `/(4(k−1))` | interior rings `[1,2]`, `[2,3]`, … | centre | interior sum |
|---|---|---|---|---|---|---|---|---|---|
| 5 | `lc_m5germC_r12` (LP) | 1 | 20.7489 | 4.2511 | +4.3590 | 0.2724 | −0.1077 | −0.0001 | −0.108 |
| 5 | same at `f*` | 1.0020781 | 20.7920 | **4.2080** | +4.3348 | 0.2709 | −0.1246 | −0.0022 | −0.127 |
| 5 | shipped `s21` | 1 | 20.8947 | 4.1053 | +4.2771 | 0.2673 | −0.1647 | −0.0071 | −0.172 |
| 6 | shipped `s32` (point cover) | 1 | 31.7135 | 4.2865 | +4.7067 | 0.2353 | −0.3495 | −0.0707 | −0.420 |
| 7 | `lc_m7A_r1` (LP) | 1 | 43.7882 | 5.2118 | +5.5256 | 0.2302 | −0.2915, −0.0223 | +0.0000 | −0.314 |
| 7 | same at `f*` | 1.0152930 | 44.4578 | **4.5422** | +5.2431 | 0.2185 | −0.5406, −0.1450 | −0.0153 | −0.701 |
| 7 | shipped `s45` | 1 | 44.7735 | 4.2265 | +5.1099 | 0.2129 | −0.6581, −0.2029 | −0.0225 | −0.883 |
| 8 | `lc_m8B_r4` (LP) | 1 | 58.3985 | 5.6015 | +5.9858 | 0.2138 | −0.3169, −0.0624 | −0.0051 | −0.384 |
| 8 | same at `f*` | 1.0171875 | 59.4023 | **4.5977** | +5.6075 | 0.2003 | −0.6661, −0.2698 | −0.0739 | −1.010 |
| 8 | candidate | 1.025 | 59.8587 | 4.1413 | +5.4354 | 0.1941 | −0.8249, −0.3641 | −0.1052 | −1.294 |

(`k = 5`'s LP round is a 12-round phase-B-closed cover, `k = 7, 8` are phase-A rounds 1 and 4, so the `f*` overheads
are not like-for-like: `k = 5` was closed to `0.2 %`, `k = 7, 8` stop at `1.5–1.7 %`.)

Readings:
* **[measured] The whole saving is a wall-ring effect.**  At the LP the interior (everything at distance `≥ 1` from
  the wall) is within `0.1–0.4` of Lebesgue in total, and every interior box `[w, k−w]²`, `w ≥ 1` integer, has
  `μ_½ ≥ area` (deficit `≤ 0`) in every cover in the table.  The wall ring `[0,k]² \ (1,k−1)²` carries
  `+4.36, +4.71, +5.53, +5.99` (`k = 5..8`, LP / point cover for 6) — more than the saving; the interior pays it back.
* **[measured] The wall-ring deficit grows, but sub-linearly per unit wall**: `0.272 → 0.235 → 0.230 → 0.214` per unit
  of wall length (the ring has `4(k−1)` unit cells).  Increments of the LP saving: `+0.96` (5→7, two rungs), `+0.39` (7→8).
* **[measured] What the exact validity costs grows with the area**: the `f*` overhead `(f*−1)·total` is `0.04`
  (`k = 5`, closed), `0.67` (`k = 7`), `1.00` (`k = 8`), and it is paid in the interior rings (at `f*` the interior
  sum is `−0.13, −0.70, −1.01`) — i.e. the interior at `f*` sits `≈ 1.5–1.7 %` above Lebesgue, which over area
  `(k−2)²` is `0.5–0.6` at `k = 7–8`.
* **Savings at `f*` (the achievable number for that cover): `4.21, —, 4.54, 4.60`** (`k = 5, 7, 8`) [proved upper ends
  via zmx2 on the scaled files; the shipped `k = 5, 6, 7` values `4.11, 4.29, 4.23`].  So far the savings stay above 4
  with room, and the room did not shrink from `k = 7` to `k = 8` (`+0.06`).
* **[heuristic] Trend.**  The LP saving grows by `≈ 0.4–0.5` per rung while the phase-A validity overhead grows by
  `≈ 0.3` per rung (area × `≈ 1.6 %`), so the achievable saving is creeping up slowly (`+0.06` from 7 to 8), not
  falling.  If the per-unit-wall deficit keeps decaying (`0.23 → 0.21`, `−7 %` per rung) the LP saving gain per rung
  shrinks (`≈ 4 × 0.2 − decay`), while the overhead at a fixed relative validity gap grows like `(k−2)²`.  Unless
  the relative overhead is driven down with `k` (the phase-B closing loop reached `0.2 %` at `k = 5`; at `k = 7, 8` it was
  not needed and not run), the two cross at a `k` of the order of 10–12 (extrapolation of four points; nothing more).
  Closing loops (or the zmx2-UNCERT → rows loop) that bring `f*` to `≈ 1.003` would keep the achievable saving at
  `≈ LP saving − 0.2`, i.e. growing with `k`.

## 5. Performance notes, and what would make `k = 9` cheaper

* The LP, as at `m = 7`, is the bottleneck: single-threaded HiGHS IPM, 1–1.5 h per solve at `k = 8`.  Total wall time
  from launch to candidate: `≈ 8 h` (chain B: 5 rounds, `7.0 h`).
* **Initial rows matter.**  Warm rows from a transplant are concentrated at its seams; without a coarse uniform row
  family the first LP collapses onto the lines (chain A).  A permanent family of `≈ 12k` rows (pitch 0.1 × 10 angles)
  fixed it at no measurable LP cost.
* The germ family is `≈ 69k` permanent rows at `k = 8` (`52k` at `k = 7`); it dominates the row count after pruning.
  At `k = 9` it will be `≈ 90k`.  Options: generate it lazily (separate it on the current cover each round instead of
  keeping all of it permanent); prune `--prune-keep` harder; or drop the phase-A lattice density in the interior.
* `zmx2` stays a fast oracle at `k = 8` (`d4`: 3–4 s wall on 7 threads for a verifying run, up to a few minutes on
  refused ones with thousands of uncertified boxes), so bisection costs minutes per round.

## 6. Files (`runs/`, gitignored)

| file | what |
|---|---|
| `s8_tp_m7A1.txt`, `s8_tp_s45.txt` | transplants of `lc_m7A_r1` and of the shipped `s(45)` cover to `s = 8` |
| `s8_uniform_rows.txt` | the uniform row family of chain B |
| `lc_m8A.sh`, `lc_m8A.{log,nohup.log,json}`, `lc_m8A_r{0..4}.txt` | chain A (stopped after A4) |
| `lc_m8B.sh`, `lc_m8B.{log,nohup.log,json}`, `lc_m8B_r{0..}.txt` | chain B |
| `m8zmx/*.log` | zmx2 bisection runs (per-factor logs with UNCERT boxes) |
| **`s60_mixed_candidate_8.txt`** | the candidate, total `748233441/12500000` |
| `zmx2_s60/cand_{d4,full}.{log,out}` | zmx2 on the candidate |
| `zm_mixed_s60.sh`, `zm_mixed_s60/{roots.jsonl,manifest.json,run.log}` | the `zm_mixed.py` run |
| `s60_ring_k567.txt`, `s60_ring_k8.txt`, `s60_ring_table.txt` | §4 profiles |

Not done: no `certificates/s60/` bundle, no commit (brief); no Lean.
