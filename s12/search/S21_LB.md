# A lower bound for `s(21)`: `s(21) >= 5000/1001 = 4.995005` (task s21-lb; 2026-09-23)

Method: exactly that of `s(12) >= 15680/3951` (`TIGHTEN.md`): a weighted point set in `[0,t]^2` such
that every closed unit square inside captures weight `>= 1`, total weight `< 21`.  At `t < 5` the
margin is positive, so the angle-net verifier `verify/` and `xcheck.py` apply; no zero-margin
machinery.  Labels: **certified** (exact check by `verify/` *and* `xcheck.py`), **verify-only**
(exact `verify/` verdict, no `xcheck.py` run), **heuristic** (float LP / not a bound).

## 0. Answer

> **`s(21) >= 5000/1001 = 4.995004995` (certified).**  `certificates/s21/s21_lower_4.9950.txt`: 4604
> points (580 D4 orbits), coordinates over `D = 1001`, weights over `10^7`, total weight
> `260057/12500 = 20.8045600 < 21`.  `verify` min `10000083/10^7` at `N = 6000` and at `N = 12000`;
> `xcheck.py --all` at `N = 6000` (all 2486 bins, 87 min on 12 workers): same minimum `10000083/10^7`, same
> bin `k = 684`.  sha256 `c8e8f878…0994f2ef` (`certificates/s21/SHA256SUMS`).
>
> Literature (§1): best previous lower bound **`4.7438`** (Friedman, DS7 survey, 2009; proof not
> published there); best *published-with-proof* general bound `1 + sqrt(14) = 4.741657` (Nagamochi 2005).
> Upper bound: **5** (the `5 x 5` grid; no packing of 21 below side 5 is known, and `s(22) = 5`, Bentz).
> So the gap for `s(21)` goes from `[4.7438, 5]` to **`[4.995005, 5]`**, width `0.2562 -> 0.0050`.

## 1. Literature (`../site/data/lower_bounds.json`, `notes/proof-anatomy.md` §5)

| bound | value | source | status |
|---|---|---|---|
| area | `sqrt(21) = 4.582576` | trivial | proved |
| `sqrt(n - 2 floor(sqrt n) + 1) + 1 = 1 + sqrt(14)` | `4.741657` | Nagamochi 2005, EJC 12 R37, Thm 2(ii) | proved (published) |
| DS7 Table 2 "21: 4.7438 Friedman" | `4.7438` | Friedman, *Packing unit squares in squares*, EJC DS7 (2009) | claimed; unavoidable set not shown |
| **this note** | **`5000/1001 = 4.995005`** | `certificates/s21/` | certified (exact, two checkers) |
| upper bound | `5` | `5 x 5` grid; Ellsworth's catalogue lists no record for `n = 21` (so best known is 5) | — |

Neighbours for scale: `s(22) = s(23) = s(24) = 5` (Bentz 2016/2018, Nagamochi 2005, Friedman 1999);
`s(20)` best known 5, lower bound `1 + sqrt(13)`.  Whether `s(21) = 5` is open; `S21_KILL.md`
(`nu_f^closed(5) >= 20.6478`, certified) does not rule it out, and `S21_COVER.md` found no valid
closed cover of `[0,5]^2` below 21.  This bound is the positive-margin shadow of that route.

## 2. Scaling (which way)

A cover `P` of `[0,5]^2` (unit squares) scaled by `t/5` is a point set in `[0,t]^2`; a unit square in
`[0,t]^2` is the image of a square of side `5/t > 1` in `[0,5]^2`, which contains a concentric unit
square.  So *keep the integer coordinates, shrink the container*: header `5000 Dp / Dp`, container
`5000/Dp` (`tighten.py --Dp`).  A **valid** closed `s = 5` cover would give every `t < 5` for free;
none exists, so the question is how much the extra `5/t - 1` buys against the cover's deficit.

## 3. Results

| `t` | points | total | check | status |
|---|---|---|---|---|
| `5000/1010 = 4.950495` | 4796 (R3 as is) | `20.801370` | `verify` N=2000 min `>= 1` VERIFIED | verify-only |
| `5000/1002 = 4.990020` | 4796 (R3 as is) | `20.801370` | `verify` N=2000 min `0.994230` (so `x 1/0.99423` -> `20.922`) | verify-only (file not built) |
| `5000/1001 = 4.995005` | 4796 (R3 as is) | `20.801370` | `verify` N=2000 min `0.973355` (k=1: net shrink `~1e-3` eats the margin) | — |
| **`5000/1001 = 4.995005`** | **4604** | **`20.8045600`** | `verify` N=6000 & 12000 min `1.0000083`; `xcheck.py --all` N=6000 | **certified** |
| `25000/5003 = 4.997002` | 4604 (above, reweighting only by scaling) | — | `verify` N=6000 min `0.965789` | not a certificate |
| `50000/10008 = 4.996003` | same | — | `verify` N=6000 min `0.977366` | not a certificate |

R3 = `runs/s5convR3_final_last.txt` (`S21_COVER.md`, strict float min `0.98978` at `s = 5`).

**How the certificate was made.**  `tighten.py reopt` (exact-verifier cutting planes at fixed points)
on R3 at `Dp = 1001`, `--topk 20 --prune-at 600000`.  LP (heuristic): `20.687` (it0) -> `20.7833` (it14),
then flat at `20.783269` with degenerate re-solves (`probe_min 0.9990-0.9994`, `176-540` violated
placements per round, `~300 s` per round).  Stopped at it19; the probe file of that round (weights
`x/(1+1e-6)`, `W = 10^12`, `runs/s21lb/B1001_snap.txt`, total `20.783248`) has `verify` min
`998986818183/10^12`; `search/rescale_weights.py` multiplies every weight by the inverse and rounds up
to `W = 10^7` -> total `20.8045600`.  Cost of stopping early instead of converging: `0.021` of weight.

The first attempt (`s21lbA1001`, default `--topk 6`, prune at 120k rows) cycled: pruning every
~7 rounds threw away tight rows (`probe_min` back to `0.93`), LP `20.776` after 27 rounds, killed.

## 4. Where it stops (heuristic)

* Headroom at `t = 4.995`: `21 - 20.8046 = 0.195` of weight; the LP at fixed points (`20.783`) is
  still `0.22` below 21.
* The limit toward 5 is the angle net, not the weight: the shrunk `sigma_k`-square of bin `k` loses
  `~2/N` relative size, which must stay below the margin `1 - t/5` (`1e-3` at `t = 4.995`).  `t = 4.998`
  needs `N >= ~12000` in the loop and `24000` in `finalize`, i.e. `2-4x` the cost per round.  With the
  loop's LP at `~20.78-20.80` and R3's closed LP `20.80` at `t = 5`, a certificate at `t = 4.998-4.999`
  looks feasible (heuristic) at ~2 h loop + ~3 h `xcheck` (N = 12000) on 12 cores.  Not run.
* `xcheck.py` time is the practical bottleneck: `~200 core-s` per bin with 4604 points (before the
  speed-up below), i.e. ~11 h on 12 cores.  Sparsifying (`tighten.py sparsify`) would cut it
  roughly quadratically; not run (each reweighted round is a full cutting-plane loop here, ~1 h).

## 5. Code changes (backwards compatible)

* `xcheck.py`: a fast path for **plain** certificates (no clique boxes in the bin, no anchor pieces,
  no region trailer).  There every clique/anchor credit and every region requirement is identically 0,
  and the old loop built 8 `Fraction`s per arrangement cell to compute those zeros.  The fast path
  computes the same integer `pre[k1] - pre[k0]` and the same witness.  Speed-up 8.5x
  (25 bins of this certificate: 458 s -> 54 s on 12 workers).  Per-bin minima identical to the
  unmodified checker on bins `k = 0, 500, 1000, 1500, 2000` (`runs/s21lb/xcheck_orig_s500.log` vs
  `xcheck_fast_s100.log`), and full `-v` output byte-identical to it (N = 2000, stride 7) on
  `s12_56points_3.8.txt`, `s11_lower_3.8143.txt` and a rejected mutant (the 56 points at container
  `1520/397`: both say `4/5`, NOT VERIFIED).  Certificates with cliques, anchors or a region take the old path unchanged.
* `search/tighten.py`: `--prune-at N` (default 120000 = the old hard-coded value).
* `search/rescale_weights.py` (new): uniform weight rescale by an exact verifier minimum, ceil to `W`.

## 6. Reproduce

```sh
cd verify && cargo build --release && cd ..
taskset -c 20-31 python3 search/tighten.py reopt runs/s5convR3_final_last.txt s21lbB1001 --Dp 1001 --n 21 \
    --topk 20 --prune-at 600000          # stop once the LP value is flat (it19 here); ~1 h
cp runs/tight_s21lbB1001_probe.txt runs/s21lb/B1001_snap.txt
verify/target/release/verify runs/s21lb/B1001_snap.txt 21 6000 12 0            # min 998986818183/10^12
python3 search/rescale_weights.py runs/s21lb/B1001_snap.txt certificates/s21/s21_lower_4.9950.txt \
    998986818183 1000000000000                                                # total 260057/12500
# checks
verify/target/release/verify certificates/s21/s21_lower_4.9950.txt 21 6000  12 0   # ~2 min
verify/target/release/verify certificates/s21/s21_lower_4.9950.txt 21 12000 12 0   # ~6.5 min on 6 threads
python3 xcheck.py certificates/s21/s21_lower_4.9950.txt 6000 --all --n 21 -j 12    # 5201 s wall, 12 workers; min 10000083/10^7 at k=684
python3 search/export_points.py --roundtrip certificates/s21/s21_lower_4.9950.txt certificates/s21/s21_lower_4.9950.json
```

The loop's snapshot round is not fixed by the command line (as in `TIGHTEN.md`); the snapshot used is
kept in `runs/s21lb/B1001_snap.txt` and regenerates the certificate byte for byte.

## 7. Files

| file | what |
|---|---|
| `certificates/s21/s21_lower_4.9950.txt` (+ `.json`, `SHA256SUMS`) | **the certificate** |
| `runs/s21lb/B1001.log`, `runs/tight_s21lbB1001.log` | cutting-plane loop |
| `runs/s21lb/B1001_snap.txt` | loop snapshot (W = 10^12) the certificate was rescaled from |
| `runs/s21lb/xcheck_full_B1001.log` | exhaustive `xcheck.py` run, per-bin minima |
| `runs/s21lb/A1001.log` | first loop (cycled; killed) |
| `runs/s21lb/r3_D*.txt`, `Bsnap_D*.txt` | scaled covers of §3 |

Not wired into `verify.sh` / CI (a full `xcheck.py` run is 87 min on 12 workers).
