# s(20) > 3 + 4√2/3 by a point cover at side 4.886 (task s20-lb, 2026-10-03/04)

> **Review (2026-10-04, `tasks/s20-review/REPORT.md`): no soundness defect; all runs reproduced box for box.**  Superseded as a bound by wand125's s(20) ≥ 1959/400 (jlevy); kept as an independent confirmation.

Goal (from `CEILINGS_17_20.md` §4): a closed cover of `[0,t]²`, `t > 3+4√2/3 ≈ 4.8856181`, with total `< 20`,
checked by `zmx2` and `zm_mixed`.  With Wainwright's packing of 19 squares in side `3+4√2/3` this gives
`s(19) ≤ 3+4√2/3 < s(20)`.  Labels as in `CEILINGS_17_20.md`: **[proved]** = exact checker verdict / exact
arithmetic, **[measured]** = float, **[heuristic]**.

## 0. Answer

> **`search/s20lb_cover_4886.txt`** (mixed v1, points only): side `t = 2443/500 = 4.886`, 12,864 points
> (all positions distinct), exactly D4-invariant, **total `249862891/12500000 = 19.98903128 < 20`**.
> sha256 `e0e25a06…6338af97` (full: `e0e25a0644fe8db1c2a346a1295eb8e966aedb2b36658d89d75813396338af97`).
>
> | checker | mode | verdict | roots / boxes / max depth | CPU |
> |---|---|---|---|---|
> | `zmx2` (commit `b8157df`) | `--d4 --sym-atoms` | **`VERIFIED-D4`**, 0 uncertified | 2,500 / 683,684 / 30 | 44 s |
> | `zmx2` | `--d4` (default atoms) | **`VERIFIED-D4`** | 2,500 / 689,800 / 30 | 30 s |
> | `zmx2` | `--full --sym-atoms` (no symmetry used) | **`VERIFIED`** | 19,208 / 4,525,820 / 30 | 280 s |
> | `zmx2` | `--full` | **`VERIFIED`** | 19,208 / 4,570,032 / 30 | 207 s |
> | `zm_mixed.py` (unmodified, HEAD sha `1fd20346…`) | `--d4 --cert-mode --disj --chain-from 0 --depth 24 --pitch 2443/50000 --ubins 16` | **`VERIFIED-D4`**, 0 uncertified | 40,000 / 67,766 / 17 | 1,453 s |
> | `zm_mixed.py` | same, `--full` | **`VERIFIED`** (no symmetry), 0 uncertified | 320,000 (160,000 × cover and its `x↔y` image) / 545,216 / 17 | 11,297 s |
>
> [proved, modulo the checkers]  So `s(20) ≥ 4.886`.  **Exact side check:** `(t − 3)·3/4 = 2829/2000` and
> `(2829/2000)² = 8003241/4000000 > 2`, so `(t−3)·3/4 > √2`, i.e. `t > 3 + 4√2/3` (by `0.000382`).  Hence
> **`s(20) > 3 + 4√2/3 ≥ s(19)`**, i.e. `s(19) < s(20)` strictly.
>
> * **LP.** Phase A (warm from `lc_ceil_t4886_r4`) converged at **`19.6368`** (A4; increments `−0.022, −0.012, −0.006`);
>   with column generation in phase B the LP of the final round is **`19.6162`**.  Exact ceiling `L(4.886) = 19.3708`
>   (`CEILINGS_17_20.md`), so `ν_f^closed(4.886) ∈ [19.371, 19.616]`.
> * **Closing overhead.** Final LP round `lc_s20lb_C_r5`: float strict minimum `0.981714`, **exact threshold
>   `f* ∈ (1.0185156, 1.0186328]`** (`zmx2 --d4 --sym-atoms` bisection: refused at the lower end, verified at the
>   upper) → certified cost `total·f* = 19.98183`; overhead `1.86 %` over the LP.  Shipped at `f = 1.019`
>   (total `19.98903`, `0.055 %` under 20).  The float strict minimum and the exact threshold agree here
>   (`1/0.981714 = 1.01863`): the zmx2-UNCERT rows (§2) removed the hidden dips that made `f*` exceed the float
>   estimate by `0.3 %` in the earlier loop.
> * **Mutations** (both checkers refuse): the file `× 0.998` (total `19.949`): `zmx2 --d4 --sym-atoms` 88 uncertified
>   boxes, `--full --sym-atoms` 346, exact `μ = 0.998382 < 1` at `(1.312061, 2.484766, θ = 27.28°)`
>   (`zmx2_tools.py uncert`); `zm_mixed --d4 --cert-mode` 438 uncertified.  The heaviest point `(0.6, 1.0)` (weight
>   `0.0596`) removed: `zmx2 --full --sym-atoms` 3,671 uncertified, float μ `0.9594` at the corner germ.
> * Total CPU ≈ **20 CPU-h** (§4).

## 1. Checker changes (zmx2 only; zm_mixed untouched)

`zmx2 cert` refused sides that are not a multiple of 1/10 (1/5 under `--d4`).  Commit `b8157df` (separate from the
certificate; the message has the argument):

1. **Roots.**  The 1/10 root cells now cover `[0, ⌈10s⌉/10]` (`--d4`: `[0, ⌈5s⌉/10]`), a superset of the required
   region.  Every pose of the required region still lies in a closed root; the extra poses are genuine poses of
   the container (redundant) or are rejected by Lemma E, which uses the exact side.  At a side that is a multiple
   of 1/10 (1/5) the root list is unchanged.
2. **Atom lines.**  The unit-grid atom lines were `x = k`, `1 ≤ k < s`.  They are now `x = k` and `x = s − k`,
   `1 ≤ k ≤ s − 1` (same for `y`): identical for integer `s`; at `s = 4.886` the set `{1, 1.886, 2, 2.886, 3, 3.886}`
   is invariant under `x → s − x`, and `x = 4` (`0.886` from the wall) is dropped.  Atom assignment only decides
   which lemma bounds a point's mass (the line lemmas hold for any axis-parallel line, as for segment lines), so it
   affects strength, not soundness.  Why: with the old set, points on `x = s − 1` (the mirror of the heavy line
   `x = 1`, mass `2.0` each) were ordinary points and `zmx2` could not certify the `θ → 0⁺` germs along them
   (153 boxes at `y ≈ 2.386` on the warm-start cover); with `x = 4` kept, its 11 atoms near the left wall lost
   their certain-inside credit at the corner germ.

**Regression [proved]:** integer sides unchanged: `s(21)` `--d4` 1,826,222 boxes, `--full` 14,709,448;
`s(13)` `--d4` 47,162; `s(60)` `--d4` 2,617,534 with `roots.log` identical root for root to
`certificates/s60/zmx2_d4/roots.log`, `--full` 21,036,120 (all = the shipped censuses).
`search/zmx2_tests.sh --quick`: 62/62 PASS.  **Non-integer side:** toy grid covers at `s = 4.886`
(`runs/s20lb_zmx/toy_{ok,rej}.txt`: points on `{0,…,2, 2.443, …, 4.886}²`, weights 1 / 0.49): `toy_ok`
verifies (`--d4`, `--full`, and `zm_mixed --d4`), `toy_rej` refused by `--d4` and `--full --pair-points
--sym-atoms`.  **Soundness harness** on the final cover (`zmx2_tools.py sound`, exact rational `μ` vs the
printed box bound): 300 random boxes (5,308 poses), `--refl` 120 (2,069), and 120 each near `(3.386, 1.2)`, near
`(3.386, 2.0, small u)` (right edge on the mirror line `x = 3.886`, `θ → 0⁺`), the wall germ `(0.5, 2.386, small u)` and the far corner `(4.3, 4.3)`: **0 FAIL**.

`zm_mixed.py` needed no change: `zeromargin.d4_roots` only needs the pitch to divide `s/2`, so the shipped settings
run with `--pitch 2443/50000` (`= (s/2)/50`, 2,500 centre cells × 16 `u`-bins = 40,000 roots, the same count as
`s(21)`).

One `zmx2` blind spot remains and is avoided on the cover side: at the corner germ `(0.5, 0.5, θ → 0⁺)` a point at
`(1, 1)` is inside `Q` only to second order (`X = ½ − (w−1)²/2` at the admissible boundary), and `zmx2` never
credits it (it bounds boxes that include inadmissible poses).  The LP put weight `0.056` there; the loop bans that
orbit (`--ban 1,1`, LP cost `≈ 0.004`).  `zm_mixed` certifies such boxes (its ADM leaves).

## 2. Method

`search/line_cover.py loop … --q 0` (point columns; no segment densities: at a non-integer side the unit grid lines
carry no tile germs).  Two runs:

* **`s20lb_B`** (`runs/s20lb_B.sh`): warm start `--warm/--cols-from runs/lc_ceil_t4886_r4.txt`, `--ban 1,1`,
  phase A 5 rounds then phase B (the close_shifted recipe).  `zmx2` on its round covers showed hidden dips: B7 had
  float strict `0.98686` but `f* = 1.01904` (exact min `0.98131`): the float scans and the tolerant (`≤ ½ + TOL`)
  LP rows count every point on ∂Q at an arrangement vertex, so a pose where several points touch the boundary
  looks fine while poses `10⁻⁶` away lose them.
* **`s20lb_C`** (`runs/s20lb_C.sh`): restarted from B7 (`--lp-rounds 0 --colgen-b`), plus the
  near-tight-row lever of `notes/jlevy-s17-techniques.md` §2 C in float form: `search/s20lb_zfeed.sh` runs `zmx2
  --d4 --sym-atoms` on every round cover scaled by `1.005` and `1.012`, `search/s20lb_zrows.py` samples float
  poses in shells of radius `2·10⁻⁶ … 10⁻³` around each UNCERT box and keeps the lowest ≤ 6 below 1, and the loop
  takes them as permanent rows (`runs/lc_s20lb_C_inject.txt`, read once per round).  The restart cost three
  rounds of row rebuilding (its warm rows are only the coarse separation), then:

| run / round | LP | exported total | points | rows | float strict min | float cost | `zmx2` `f*` (→ `total·f*`) |
|---|---|---|---|---|---|---|---|
| ceil `r4` (start) | 19.6205 | 19.6206 | 5,036 | 96k | 0.948 (pol) | — | 1.0595 (with (1,1)) |
| B A0–A4 | 19.6363, 19.6762, 19.6547, 19.6423, **19.6368** | | 5.6k–10.4k | 28k–52k | lat 0.933 → 0.985 | | A3: 1.0403 |
| B B5 | 19.6419 | 19.6420 | 9,108 | 52k | 0.96793 | 20.293 | |
| B B6 | 19.6700 | 19.6701 | 15,188 | 66k | 0.98440 | 19.982 | 1.01875 (20.039) |
| B B7 | 19.6829 | 19.6830 | 12,096 | 78k | 0.98686 | 19.945 | 1.01904 (20.058) |
| C B0–B4 | 19.663, 19.496, 19.591, 19.606, 19.612 | | | 18k–84k | 0.880 → 0.954 | | |
| **C B5** | **19.6162** | **19.6162** | 12,864 | 95k (64k perm.) | **0.98171** | **19.982** | **1.0186328 (19.9818)** |

[measured] except the `f*` column [proved upper ends].  LP solves: HiGHS IPM, single-threaded, 6–37 min.
Stopped after C B5 (first round below 20).

## 3. Reading

* The point LP at 4.886 converges at `≈ 19.62–19.64`, below the `19.65–19.7` guessed in `CEILINGS_17_20.md`;
  room to 20 is `≈ 1.9 %`.  The validity overhead of the closed point cover is `1.86 %` here (point covers at
  `m = 5`: `1.4 %`), so the margin is thin (`0.055 %` in the shipped file).
* The overhead is dominated by tilted interior dips near `θ ≈ 27–42°` (the last refused boxes of the bisection
  sit at `(1.31, 2.48, 27.3°)`), not by wall or tile germs.
* `zmx2` tracks the true minimum on this cover (its `f*` equals `1/(float min)` to `10⁻⁵` once the hidden vertex dips
  are rowed), i.e. the checker overhead is negligible; the remaining cost is genuine LP-row sampling.

## 4. CPU

LP loops ≈ 12 CPU-h (B 3.3 h wall, C 3.5 h wall; single-threaded LP plus a 6-process pool for scans); `zmx2`
(bisections, feed, finals, regressions, tests) ≈ 2 CPU-h; `s20lb_zrows.py` ≈ 1 CPU-h (one early run used all
cores through numpy threads before `OMP_NUM_THREADS=4`); `zm_mixed` `--d4` 0.4 CPU-h, `--full` 3.1 CPU-h,
mutation 0.3 CPU-h.  Total ≈ **20 CPU-h** (estimate from wall times).

## 5. Files

| file | what |
|---|---|
| **`search/s20lb_cover_4886.txt`** | the cover: `lc_s20lb_C_r5 × 1019/1000`, total `249862891/12500000` |
| `runs/s20lb_{B,C}.sh`, `runs/lc_s20lb_{B,C}.{log,json}`, `runs/lc_s20lb_{B,C}_r*.txt` | loops (gitignored `runs/`) |
| `runs/s20lb_zmx/cand_{d4,full}.log`, `cand_full.out` | `zmx2` finals |
| `runs/s20lb_zmm/{run.log,manifest.json,roots.jsonl}`, `runs/s20lb_zmm_full/…` | `zm_mixed` runs |
| `runs/s20lb_zmx/mut_*`, `zmm_mut_x0998.log` | mutation tests |
| `search/s20lb_zbisect.sh`, `s20lb_zrows.py`, `s20lb_zfeed.sh` | tools (commit `442d9bc`) |

Not done (brief): no `certificates/` bundle, no site/README edits, no Lean.

## 6. Reproduce

```
bash runs/s20lb_B.sh                       # warm from runs/lc_ceil_t4886_r4.txt; stopped after B7
bash runs/s20lb_C.sh & bash search/s20lb_zfeed.sh s20lb_C 0 &   # stopped after B5
ZTHREADS=6 ZFLAGS=--sym-atoms bash search/s20lb_zbisect.sh runs/lc_s20lb_C_r5.txt 1.0 1.03 8   # f* in (1.0185156, 1.0186328]
python3 search/line_cover.py scale runs/lc_s20lb_C_r5.txt search/s20lb_cover_4886.txt --factor 1.019
Z=verify2/target/release/zmx2
$Z cert search/s20lb_cover_4886.txt --d4 --sym-atoms --threads 6;  $Z cert search/s20lb_cover_4886.txt --full --sym-atoms --threads 8
bash runs/s20lb_zmm.sh; bash runs/s20lb_zmm_full.sh                # zm_mixed --d4 / --full, --cert-mode
python3 -c "from fractions import Fraction as F; a=(F(2443,500)-3)*F(3,4); print(a*a>2)"   # t > 3+4*sqrt(2)/3
```
