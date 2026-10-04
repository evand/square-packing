# Pure-cover ceilings for n = 17–20 (task ceilings-17-20, 2026-10-03)

Brief: `tasks/ceilings-17-20/README.md`.  For each container side `t`, bracket what a **pure closed cover**
can prove: a nonnegative measure `w` on `[0,t]²` (points, grid-line densities, anything) with every closed unit
square capturing `w(Q) ≥ 1` and total `< n` proves `s(n) ≥ t`.  Two ends:

* **ceiling** `L(t)`: an exact D4-symmetric fractional packing of closed unit squares in `[0,t]²` with max
  coverage `M ≤ 1`.  Weak duality (`CEILING.md`): every such cover has total `≥ L(t)`, for any `w`, and also at every
  side `≥ t` (`ν_f` is non-decreasing).  So `L(t) ≥ n` proves **no pure closed cover gives `s(n) ≥ t'` for any
  `t' ≥ t`**.
* **cover LP** `U(t)`: the LP value of `search/line_cover.py loop --s t` (point covers, sampled rows).  It is a
  float LP over finitely many rows, so it is not a bound either way; it estimates `COVER(t)` from below while
  phase A is still adding rows.

Labels: **[proved]** = exact `Fraction`/integer check by `search/dual_exact.py check` (from the support file, no
LP; with `--n N` the exit status asserts `L > N`); **[measured]** = float LP / float scan reported as-is;
**[heuristic]** = extrapolation.

## 0. Answer

| side `t` | ceiling `L(t)` [proved] | cover LP [measured] (round, lattice min) | targets decided |
|---|---|---|---|
| **4.660** | **`17044695833/10⁹ = 17.044696`** | not run (decided by the ceiling) | (17, 4.6604), (17, 4.6755): **killed** |
| 4.823 | `18470272437/10⁹ = 18.470272` | `18.879573` (A4, 0.931); est. `18.9–18.95` converged [heuristic] | (19, 4.8229): **open, marginal** |
| **4.8856** | **`77500070540/3999999999 = 19.375018`** | — | (19, 4.8856): **killed** |
| 4.886 | `154966592872/7999999999 = 19.370824` | `19.620544` (A4, 0.962); est. `19.65–19.7` converged [heuristic] | (20, 4.8856): **open, room ≈ 1.5–1.9 %** |
| 5 | `20.647781` (`S21_KILL.md`) | `20.75` (`LINE_COVER.md`) | (20, 5): killed (sanity, as known) |

Per target:

| `n` | `t` | ceiling | cover LP | verdict |
|---|---|---|---|---|
| 17 | 4.6604 | `≥ 17.0447` (measure at 4.660, monotonicity) | — | **No pure closed cover reaches it** [proved].  Answers TODO E6: `ν_f^closed(4.6604) ≥ 17`. |
| 17 | 4.6755 | `≥ 17.0447` (same) | — | killed [proved] |
| 18 | 4.8229 | `ν_f(4.823) ≥ 18.47 > 18` | — | not a target (`s(18) ≤ (7+√7)/2` by Hämäläinen); consistent |
| 19 | 4.8229 | `18.4703` (at 4.823) | `18.8796` at A4, rising | **not killed**; LP room `0.64 %` at A4, `≈ 0.3–0.5 %` at convergence [heuristic]: below the brief's 0.5 % threshold once converged; point-cover overheads at `m = 4, 5` were 1.4–2.9 %.  Stopped. |
| 19 | 4.8856 | **`19.3750`** (at 4.8856 `< 3+4√2/3`) | — | **killed** [proved]: `s(19) > 3+4√2/3` is out of reach of every pure closed cover |
| 20 | 4.8856 | `19.3708` (at 4.886) | `19.6205` at A4 | **open, room `≈ 1.9 %` at A4** (`20/19.6205`), `≈ 1.5 %` est. converged.  Passes the 0.5 % test; certificate **not started** (CPU budget, §3) |
| 20 | 5 | `20.6478` | `20.75` | killed (known) |

So of the milestone pair, **`s(20) > 3+4√2/3` by a point cover at `t = 4.886` is plausible** (LP room 1.5–1.9 %, no
tile germs at a non-integer side), and **`s(19) > (7+√7)/2` is marginal** (LP within ≈ 0.5 % of 19, the exact
ceiling 18.47 does not exclude it).  Pure covers are **dead for s(17) at R071's 4.66044 and for s(19) at 4.8856**.

`ν_f` is steep here: `ν_f^closed` brackets `[17.04, ?]` at 4.660, `[18.47, ~18.9]` at 4.823, `[19.375, ~19.7]` at
4.8856–4.886, `20.65` at 5 — about 10–14 units of mass per unit of side.  So a 0.5 % change of mass is ≈ 0.007–0.01
of side: the n = 19 milestone sits right at the edge.

## 1. The ceilings

Code unchanged: `search/packing_dual.py` (float column generation, D4, closed squares, certified in floats at
every arrangement vertex), `search/dual_exact.py build` (snap at `Q = 10⁷`, `Dc = 10⁸`, exact `M`), `check`.
Seeds: the `t = 5` `s5killE` support (`--warm`, centres scaled by `t/5` and clamped) plus the `t = 4` `E1`
support scaled by `t/4` (`cover4_cg.py pool --scale`, as `--seed-file`); coarse seeds `0.15 / 15°`, rows `0.05`.

| run | side | float result | exact (`build`) | CPU |
|---|---|---|---|---|
| `dual_ceil4823` | 4.823 | stopped at round 33 (LP mass 18.47, `M = 1.0025`, not converged) | `--Q 10⁷` with LP polish: **`L = 18.470272437`**, 76 poses | ≈ 2.6 h |
| `dual_ceil4886` | 4.886 | **converged** round 38: `M = 1`, 0 bad vertices, `L = 19.37082` | `--no-lp`: **`L = 19.3708241`**, 174 poses | ≈ 3.0 h |
| `dual_ceil48856` | 4.8856 | warm = the 4.886 support; converged round 22: `L = 19.37502` | `--no-lp`: **`L = 19.3750176`**, 167 poses | ≈ 1.0 h |
| `dual_ceil4660` | 4.660 | warm = the 4.823 support; `BEST L = 16.91697` | `L = 16.916972794` (< 17) | ≈ 1.4 h |
| `dual_ceil4660b` | 4.660 | warm = 4660's own support, seeds `0.1 / 10°`; converged round 16: `L = 17.04470` | `--no-lp`: **`L = 17.044695833`**, 68 poses | ≈ 1.0 h |

All four exact supports pass `check` in the D4 fundamental domain and `check --full --stream` on the whole
container (49k–939k vertices, seconds to a minute each); `--n` asserts as in the table [proved].  Shipped (small)
as `search/ceil_nuf_{4.660,4.823,4.8856,4.886}_exact_support.txt` (sha256 prefixes `3f6412…`, `8e1b14…`,
`1389e5…`, `0de2b6…`; same format as `s21_nuf5_exact_support.txt`).

Notes.
* Rescaling a converged support to a slightly different side destroys it: the 4.886 support scaled to 4.8856
  (factor 0.99992, or a wall-preserving shrink of the interior only) re-polished to only `18.74 / 18.82` — the measure
  lives on exact contacts.  A short warm-started float run (20 min) recovers it fully.  The same at 4.660: the
  4.823 support scaled by 0.966 gives `L = 13`.
* The 4.823 value is a lower bound from an unconverged run (`M` still 1.0025 at the kill); the true `ν_f(4.823)` is
  in `[18.47, ≈ 18.9]`.  More column generation would raise the ceiling but cannot reach 19 if the cover LP
  converges below 19 (`ν_f ≤ COVER`).
* `L(4.8856) > L(4.886)` only says the 4.886 run converged to a slightly worse column set.
* Scope of the s(17) kill: it is about closed covers (any nonnegative measure: points, lines, area).  It says nothing
  about capacity-one rule atoms or strict-core charges (`notes/jlevy-s17-techniques.md` §2 A), which is how
  R068/R071 got past 4.614.

## 2. The cover LPs

`runs/ceil_17_20/cover.sh`: `line_cover.py loop TAG --s t --q 0 --nproc 3 --lp-rounds 8 --rounds 1` (point
columns: lattice 0.05 + dual pricing; no warm rows, no line densities — at a non-integer side the D = 1000 segment
pieces do not tile the line and there are no tile germs, so the point LP is the right object, as at `m = 5` where the
line and point LPs agree, `LINE_COVER.md` §0).  Both stopped in phase A after A4 (budget):

| round | 4.823 LP | lattice / polished min | 4.886 LP | lattice / polished min | LP solve |
|---|---|---|---|---|---|
| A0 | 16.000000 | 0 / 0 | 16.000000 | 0 / 0 | 350–370 s |
| A1 | 18.335875 | 0.630 / 0.623 | 19.000000 | 0.711 / 0.711 | 810–840 s |
| A2 | 18.627860 | 0.773 / 0.773 | 19.313694 | 0.711 / 0.704 | 1130–1240 s |
| A3 | 18.817514 | 0.900 / 0.900 | 19.585866 | 0.812 / 0.809 | 1710–1820 s |
| A4 | **18.879573** | 0.931 / 0.931 | **19.620544** | 0.962 / 0.948 | 2130–2220 s |

[measured]  Increments `+0.29, +0.19, +0.06` (4.823) and `+0.31, +0.27, +0.035` (4.886): both are flattening,
so converged values `≈ 18.9–18.95` and `≈ 19.65–19.7` [heuristic].  The exported round covers are
`runs/lc_ceil_t{4823,4886}_r{0..4}.txt` (mixed v1, `s = 4823/1000`, `4886/1000`); `zmx2 info` / `zmx2 d4` read
them (non-integer side accepted; exact D4 invariance confirmed on `lc_ceil_t4823_r2`).  No phase-B round was
reached, so no honest-cost measurement.

Cover LPs at 4.660 / 4.675 were not run: the ceiling already decides those targets.

## 3. CPU, and what was not done

Total ≈ **13 CPU-h** (float dual runs ≈ 9 h — the pricing threads, not the LP, dominate; cover loops 2 × 1.9 h; exact
builds and checks ≈ 0.2 h).  No certificate was attempted: the (20, 4.8856) target passes the brief's room test, but
a closing loop + `zm_mixed` would exceed the 15 CPU-h line (estimate below).

## 4. Next steps

* **s(20) > 3+4√2/3 (recommended, the one live milestone).**  Continue `lc_ceil_t4886` to convergence and run the
  phase-B closing loop at `s = 4.886` (`LINE_COVER.md` §1 recipe; no `--germ`, no tile family needed).  Need
  `total / f* < 20`, i.e. validity × checker overhead `≤ 20/19.7 ≈ 1.5 %`; the mixed-cover closing loops reached
  `0.15–0.2 %` at `m = 5`, point loops `1.4 %` (but those were dominated by tile germs, absent here).  Then bisect
  `f*` with `zmx2 --d4` (`search/m8_zbisect.sh`), verify with `zmx2 --full` and `zm_mixed.py --d4 --cert-mode`.
  Estimate [heuristic]: 4–8 h wall, 15–30 CPU-h (LP 30–40 min per solve at ~100k rows; zm_mixed a few CPU-h at this
  size).  Warm-start with `--warm runs/lc_ceil_t4886_r4.txt` to skip the cold A0–A3 rounds.
* **s(19) > (7+√7)/2: marginal.**  Worth one cheap check before committing: converge the 4.823 cover LP (2–3 more
  A rounds, ≈ 2 CPU-h).  Go on only if the converged LP is `≤ 18.90` (room ≥ 0.5 %).  In the other direction, more
  column generation at 4.8229 could raise the ceiling, but it can only kill the target if `ν_f(4.8229) ≥ 19`, which
  the cover LP (18.88 and flattening) makes unlikely.  Also note wand125's `s(19) ≥ 4.815` (rectangle densities,
  unreplayed) is the nearest known bound.
* **s(17): stop pure-cover work.**  `ν_f^closed(4.660) ≥ 17.0447` [proved], so no closed cover of any kind proves
  `s(17) ≥ 4.660`, a fortiori not R071's 4.66044 or Bidwell's 4.6755.  The lever is capacity-one rule atoms (N15).

## 5. Reproduce

```
D=runs/ceil_17_20
for t in 4.823 4.886 4.660; do python3 search/cover4_cg.py pool $D/E1_x$t.txt --scale $(python3 -c "print($t/4)") runs/dual_E1_support.txt; done
bash $D/dual.sh 4.823 ceil4823          # stopped at round 33 (~40 min)
bash $D/dual.sh 4.886 ceil4886          # converged, 3100 s
python3 search/packing_dual.py 4.8856 ceil48856 --warm runs/dual_ceil4886_support.txt --seed-file $D/c4886_w4.8856.txt \
    --seed-pitch 0.15 --seed-dth 15 --row-pitch 0.05 --threads 3 --time 1200 --polish-time 300 --rounds 200 --sweep-k 10
python3 search/packing_dual.py 4.660 ceil4660 --warm runs/dual_ceil4823_support.txt --seed-file $D/E1_x4.660.txt \
    --seed-pitch 0.15 --seed-dth 15 --row-pitch 0.05 --threads 3 --time 1500 --polish-time 300 --rounds 200 --sweep-k 10
python3 search/packing_dual.py 4.660 ceil4660b --warm $D/dual_ceil4660_r1_support.txt --seed-file $D/E1_x4.660.txt \
    --seed-pitch 0.1 --seed-dth 10 --row-pitch 0.05 --threads 3 --time 1200 --polish-time 300 --rounds 200 --sweep-k 10
python3 search/dual_exact.py build --t 4823/1000 --tag ceil4823 --src runs/dual_ceil4823_support.txt --Q 10000000 --Dc 100000000 --n 19
python3 search/dual_exact.py build --t 4886/1000 --tag ceil4886 --src runs/dual_ceil4886_support.txt --Q 10000000 --Dc 100000000 --n 19 --no-lp
python3 search/dual_exact.py build --t 48856/10000 --tag ceil48856 --src runs/dual_ceil48856_support.txt --Q 10000000 --Dc 100000000 --n 19 --no-lp
python3 search/dual_exact.py build --t 466/100 --tag ceil4660b --src runs/dual_ceil4660b_support.txt --Q 10000000 --Dc 100000000 --n 17 --no-lp
python3 search/dual_exact.py check search/ceil_nuf_4.660_exact_support.txt --n 17             # PASS; likewise 4.8856/4.886 --n 19, 4.823 --n 18
python3 search/dual_exact.py check search/ceil_nuf_4.660_exact_support.txt --full --stream --n 17
bash $D/cover.sh 4.823 ceil_t4823; bash $D/cover.sh 4.886 ceil_t4886   # stopped after A4 (~1.9 h each)
```
