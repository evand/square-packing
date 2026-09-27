# Germ failures of the `s(32)` check: point shifting vs. the checker (2026-09-26)

Task: make `zmcheck cert` certify the interior tile germs (S32_EXACT §2) by modifying the cover near them,
keeping total `< 32` and exact D4 symmetry.  Labels: **certified** (exact `zmcheck`, partial sweeps),
**heuristic** (float), **diagnosis** (read off `ZM_DEBUG` / the `ZM_CHAINS` scratch build).

## 0. Answer

* **Cover `runs/s32_shift_v1.txt`: total `3169794004915/10¹¹ = 31.697940 < 32`** (0.0156 *below* the candidate), 10,565
  points, exactly D4-invariant (`--d4` check), heuristically valid by the same confirmation protocol as the candidate
  (confirmed min `0.99481` before the final scale, `1.00405` after).  Construction (§3): every point on an interior grid
  line (`x = k` or `y = k`, `k = 1..5`, other coordinate in `(1, 5)`) snapped to the nearest multiple of `0.012`
  ("fewer, heavier line points, at common positions"), then the LP re-closed on that support (`close_shifted.py loop`,
  15 rounds, STABLE), then `× 10093/10000`.
* **All 16 interior-germ roots of the D4 region (u-bin 0, the four cells around each of `(1.5,1.5)`, `(1.5,2.5)`,
  `(2.5,1.5)`, `(2.5,2.5)`) certify, 0 uncertified** (certified), with `zmcheck_fast` + the `ZM_MIXPAIR` pair rule (§1) at
  `--branch-cap 640 --node-cap 16000000 --depth 22`: 9,628 CPU-s for all 16, max depth 16.  Off-germ spot checks on v1
  with the **stock** checker and sweep settings are clean (column `x ∈ [1.2, 1.3]` all 120 roots, 723 CPU-s; edge-tile and
  tilted-dip cells).  **Stock checker at the germs:** `(1.5,1.5)` and `(2.5,2.5)` cells clean (1,165 / 1,984 CPU-s); in the
  `(1.5,2.5)` and `(2.5,1.5)` cells **23 and 8 uncertified boxes** (one root each, depth 22, 22.6k / 19.2k CPU-s; runs
  `s32shift_v1plain_g{1525,2515}`), which the patched checker closes in 2.2k / 2.3k CPU-s.  So v1 fixes defect 2 of §1 (the one no
  checker setting fixes) but not defect 1 (pair choice, a count heuristic the cover cannot control robustly).
* **Not local.**  The snap touches 5,240 line points everywhere in `[1, 5]²` and the re-LP changes every weight (only
  7,317 of the positions are shared with the candidate; `Σ|Δw| = 6.9` on them).  Nothing certified for the candidate
  carries over: **the D4 sweep must be rerun on v1.**  A purely local change could not work (§2): each germ needs its whole
  tile boundary aligned, and the interior tiles' boundaries are all of the interior grid lines.
* Two checker findings (heuristic ordering only; soundness untouched) — §1.  For the *unchanged* candidate a small
  `zmcheck` patch plus a larger node cap closes the `(1.5,1.5)` and `(2.5,1.5)` germ cells but not `(1.5, 2.5)`
  (10 boxes; 4 of them still fail at 5·10⁸ nodes).  So for zmcheck alone the cover change is the better fix; the cap-free
  two-chain product (S32_EXACT §8) would be the checker-side one.

## 1. Why the germs fail (diagnosis)

At a germ, a box touching the corner `(c, θ = 0)` has two free thresholds: `t` along the two horizontal edge lines (bottom
line in for `x < t`, top in for `x > t`) and `s` along the vertical ones.  A certificate is the product `t-chain × s-chain`;
the partner line of each chain is coupled by Lemma H (`+B_i ⇒ −T_j` for `x_i ≤ x_j`, `μ = 1`, exact).
Two separate defects:

1. **Pair choice.**  `chain_order` forces the interleave of the **two largest families by count**.  The families are not
   the half-lines (≈ 70 points): they include the other half of some lines and ≈ 110 "junk" points beyond the edge, cut by
   the 640-branch cap.  At `(3.5, 1.5)` from below-left: `k2 187, k3 185, k0 162, k1 106` — both horizontal, i.e. one
   parameter twice; the vertical one gets 4 greedy levels ⇒ `FORCED FAIL w = 0.9974`, 4M nodes, 35 s.  Scratch patch
   `ZM_MIXPAIR` (promote the largest family of the other orientation, `kind / 2`): same box **DISJ, 18,491 regions, 3.6 s**.
   Diff: `runs/s32shift_mixpair.patch` (39 lines, on current `verify2/src/main.rs`; binary `runs/zmcheck_fast_mix`,
   `ZM_CHAINS=1` prints family sizes).  Not committed.
2. **Partner points between pivots.**  In a product region `a_i in, a_{i+1} out` the partner line's points with
   `x ∈ (x_i, x_{i+1})` are undetermined and drop out of the witness.  With ~70 unaligned points per half-line that is
   ≈ one partner point (`≈ 0.003–0.005`) per region, about the germ margin, so regions fail at `w ≈ 0.998` and the product
   blows up (candidate + patch at `(1.5, 2.5)`: 17M–333M nodes per depth-22 box, some > 512M).  **Snapping the lines to
   common positions removes exactly these losses**: the partner points then sit *on* the pivots.  The same `(1.5,1.5)` box
   that fails at 64M nodes on the candidate closes with 10–91 regions on `s32_shift_s12`.

More centre weight does not help (`+0.02` at each `(1.5,2.5)`-class centre: identical `FORCED FAIL`); merging pairs on
alternate lines (`s32_shift_m2`) does not either.

## 2. Snap pitch vs. certifiability (box tests, candidate weights, 16 depth-22 germ boxes)

| file | pitch | max move | points | float (targeted diffscan / hardscan 0.004) | germ boxes certified (patch) |
|---|---|---|---|---|---|
| candidate | — | — | 13,085 | 1.0066 / 1.0102 | 0 (1 at 16M nodes) |
| `s32_shift_s5` | 0.005 | 0.002 | 12,541 | — | 2 / 16 |
| `s32_shift_s10` | 0.010 | 0.006 | 11,041 | — | 10 / 16 |
| `s32_shift_s12` | 0.012 | 0.006 | 10,457 | 1.0017 / **1.0004** | **16 / 16**, 10–91 regions |
| `s32_shift_s15` | 0.015 | 0.009 | 9,865 | — | 16 / 16 |
| `s32_shift_s20` | 0.020 | 0.016 | 9,185 | **0.9914** (violation) | 16 / 16 (15 stock) |

Snapping costs coverage (≈ 1 % at pitch 0.012, all at tilted 33–39° poses), hence the re-LP (§3).  `s12` itself on the stock
checker: `(1.5,1.5)` cell clean (1,280 CPU-s), but thin tilted boxes (`θ ≈ 5.5°`) fail elsewhere — a margin problem, fixed by v1.

## 3. Building v1 (heuristic)

```
python3 snap.py runs/s32closeA_r22.txt runs/s32_shift_s12_r22.txt 12          # r22 snapped: strict 0.98967 (was 0.99902)
python3 snap.py runs/s32convD_it6_last.txt runs/s32_shift_s12_colsD.txt 12    # extra columns, snapped the same way
python3 snap.py runs/s32convR4_r11_last.txt runs/s32_shift_s12_colsR4.txt 12
runs/s32shiftL12.sh     # close_shifted.py loop, rows-from s32closeA_dips.txt; STABLE after round 15, LP 31.4053
runs/s32shiftL12_confirm.sh   # confirm --tiles --interleaved --full: min 0.9948064 at (0.52206, 3.49991, 0.09°), edge tile
python3 search/scale_cover.py runs/s32shiftL12_r15.txt 10093 10000 runs/s32_shift_v1.txt   # 1.004/0.9948064 = 1.00924
```

`snap.py` (scratch, copied to `runs/s32shift_snap.py`; needs `runs/ld.py`): snaps `t` to the nearest multiple of `P/1000` that is not an integer
(`3000 ≡ 0 mod P`, so exactly D4-symmetric), sums coincident points.  Loop: rounds 0–15, strict `0.9769 → 0.99929`,
LP `31.3706 → 31.4053` (r22 was `31.3716`: the snap costs `+0.034` LP).  The confirmed minimum is *better* than r22's
(`0.99481` vs `0.99319`), so the final scale is smaller and the total lands below the candidate.

| confirm scan (r15, unscaled) | min |
|---|---|
| hardscan 0.002 | 0.99929 |
| 80 dip boxes, pitch 0.001 | 0.99632 at (3.48466, 5.45484, 5.44°) |
| tiles ±0.015, pitch 0.0005, 0–1.5° | 0.99901 |
| full, 319 angles, pitch 0.001 | 0.99736 |
| full, 318 interleaved angles | 0.99687 |
| polish of the 1,400 worst | **0.99481** at (0.52206, 3.49991, 0.09°) |

## 4. Exact runs on v1 (certified; `--d4`, `ZM_UBINS=0` for the germ cells)

| run | checker | cells | result | CPU-s | max depth |
|---|---|---|---|---|---|
| `s32shift_v1mix_g15` | patch, 16M | around (1.5,1.5) | 0 uncert | 1,247 | 13 |
| `s32shift_v1mix_g1525` | patch, 16M | around (1.5,2.5) | 0 uncert | 2,224 | 14 |
| `s32shift_v1mix_g2515` | patch, 16M | around (2.5,1.5) | 0 uncert | 2,329 | 12 |
| `s32shift_v1mix_g25` | patch, 16M | around (2.5,2.5) | 0 uncert | 3,828 | 16 |
| `s32shift_v1plain_g15` | **stock**, 4M | around (1.5,1.5) | 0 uncert | 1,165 | 14 |
| `s32shift_v1plain_g1525` | stock | around (1.5,2.5) | **23 uncert**, all in root (14,25) bin 0 | 22,576 | 22 |
| `s32shift_v1plain_g2515` | stock | around (2.5,1.5) | **8 uncert**, all in root (24,14) bin 0 | 19,171 | 22 |
| `s32shift_v1plain_g25` | stock | around (2.5,2.5) | 0 uncert | 1,984 | 16 |
| `s32shift_v1_offA` | stock | edge tile (0.4–0.5, 2.4–2.5), (1.7, 2.3), (1.5, 2.4), all 4 bins | 0 uncert | 510 | 13 |
| `s32shift_v1_col12` | stock | column x ∈ [1.2, 1.3], all y, 4 bins (120 roots) | 0 uncert | 723 | 9 |

## 5. Same germ cells on the unchanged candidate (certified)

| run | checker | result |
|---|---|---|
| `s32zm_germ15` (other agent) | stock, old build | (1.5,1.5): 15 uncert in 2 roots |
| `s32shift_mix_r3414` (old build) | patch, 4M | root x[3.4,3.5] y[1.4,1.5] u0: **0 uncert** (5,694 s, depth 16; stock: 13 uncert, 12,441 s); killed before u7 |
| `s32shift_mix_germ15` (old build) | patch, 4M | same 15 uncert (there the pair was already mixed) |
| `s32shift_fm_g15` | patch, 16M | (1.5,1.5) region cells: 0 uncert, 3,831 CPU-s, depth 22 |
| `s32shift_fm_g2515` | patch, 16M | 0 uncert, 6,143 CPU-s |
| `s32shift_fm_g1525` | patch, 16M | **10 uncert** in (14,24) bin 0; box mode: 4 close at 128M, 2 more at 512M (333M nodes), 4 not at 512M |

## 6. Files (`runs/`)

`s32_shift_v1.txt` (**the result**), `s32shiftL12{.sh,.log,.out,.err,_r*.txt,_dips.txt,_confirm.*}`, `s32_shift_s{5,10,12,15,20,25}.txt`
(+ `.moved`), `s32_shift_c{005,01,02}.txt`, `s32_shift_m2.txt`, `s32shift_*.{sh,log,err}`, `s32shift_mixpair.patch`,
`zmcheck_fast_shiftcopy` (frozen stock), `zmcheck_fast_mix` (patched), `zmcheck_shift*` (older builds).

## 7. Next (Evan decides)

1. Rerun the D4 sweep on `s32_shift_v1.txt`, with the `ZM_MIXPAIR` rule committed (at least for the germ roots), or
   the stock checker leaves 31 boxes at the `(1.5,2.5)`-class germs of v1.
2. Independent re-check with `zeromargin.py` on v1 as planned.
3. Commit `ZM_MIXPAIR`: it changes nothing where the top two families already mix (edge-tile cell census identical,
   ADM 99 / DISJ 184 / EMPTY 72, `s32shift_tile_{plain,mix}`).  A cap-free two-chain product (S32_EXACT §8) is the
   more general checker fix; on the candidate it is needed at `(1.5, 2.5)` (> 5·10⁸ nodes), on v1 it is not.
