# The floor of `COVER(4, closed)`: an exact certificate `COVER^closed(4) >= 12.2688`

Code: `search/dual_exact.py` (now parametrised by the container side, `--t`), `search/cover4_cg.py`
(new: exact restricted-master column generation over poses), `search/packing_dual.py`,
`search/nu_f.py`.  Certificate: `search/cover4_exact_support.txt` (= `runs/dual_exact_4_support.txt`);
logs `runs/dual_exact_4.log`, `runs/dual_exact_4_check_F.log`, `runs/dual_exact_4_check_full.log`,
numbers `runs/dual_exact_4.json`; the float search that chose the poses is `runs/dual_E1.log` / `.json`.

## Statement proved

Let `C = [0, 4]^2`.  There is an explicit measure `mu` on placements of **closed** unit squares
inside `C` — 294 poses `(p, q, cx, cy, m)` with rational rotation `theta = 2 arctan(p/q)` and
rational centre, each pose carrying `m/8` on each of its 8 dihedral images (2352 closed squares) — with

    mass(mu) =  2453760771 / 200000000  =  12.268803855          (exactly)
    M        =  max_{x in C} mu({Q : x in Q})  =  1999999999 / 2000000000  =  0.9999999995   (exactly)
    L        =  mass / M  =  24537607710 / 1999999999  =  12.268803861134401...

Hence `mu / M` is a feasible fractional packing measure and

    nu_f^closed(4)  >=  L  =  24537607710 / 1999999999  =  12.2688038611 ,

and by weak LP duality (`nu_f <= COVER`, elementary, `search/CEILING.md`)

> **No weighted point set (LP cover) of total weight `< 24537607710/1999999999 = 12.2688038611`
> covers every closed unit square inside `[0,4]^2`** — whatever point set, cell decomposition,
> LP or verifier is used, and a fortiori at any container side `>= 4`.

`L - 12 = 537607722/1999999999 = +0.26880386`.  Every load-bearing quantity is a Python integer or
`Fraction`; floats only *choose* which poses carry mass.

**Baseline improved.**  The previous floor at `s = 4` was the `s = 3.99` certificate of
`search/DUAL_EXACT.md` read at `s = 4` — valid because `[0,3.99]^2 subset [0,4]^2`, so the same
measure is admissible and has the same maximum coverage:
`L(3.99) = 96065846032/7999999981 = 12.008230782520`.  This note raises the floor by
`+0.2605731` to `12.2688038611`.

## The bracket

| side | value | status |
|---|---|---|
| **`L` (this note)** | **`24537607710/1999999999 = 12.26880386`** | **exact, certified — a theorem** |
| `nu_f^closed(4)` heuristic (float LP mass, all rows) | `12.2688` (converged: the float polish reached `M = 1` at mass `12.268806`) | heuristic |
| `COVER^closed(4)` heuristic cover LP | `12.39 - 12.42` total (`CLOSED4.md`, `FAMILY.md` §2) | heuristic |
| **`U` (best honest heuristic cover cost on record)** | **`~12.46`** (`FAMILY.md` §1: the rung-2 closing loop plateaus at `12.3883` LP total / `0.9958` stress minimum); the only *exported* point set, `runs/closed4_best.txt`, is `12.417` total / **`12.509`** honest | heuristic, **not** a proved upper bound |

    12.2688038611  <=  COVER^closed(4)  <=  14          (both ends proved; 14 is Friedman / Nagamochi)
    12.2688038611  <=  COVER^closed(4)  <=  U ~ 12.46   (U is the best honest heuristic cover cost, NOT a proof)

So the working bracket is `[L, U] = [12.2688038611, ~12.46]`, of width `< 0.2` — but only its left
end is a theorem.  The LP proof of `s(13) = 4`
(rung 2, `tasks/rung2-s13`, `search/RUNG2.md`) provably costs **at least `0.2688` more than 12** —
but `12.2688 < 13`, so nothing here rules rung 2 out.  What it does rule out is any hope that the
cover LP at `s = 4` gets near 12: it cannot go below `12.2688`, and the heuristic upper side has
been stuck at `12.39 - 12.51` for two sessions, which is consistent with the true value being
`~12.3 - 12.45`.

## What is exact, and how

Identical machinery to `search/DUAL_EXACT.md`, now at `t = 4/1` (`dual_exact.py --t 4`); only the
differences are repeated here.

**Poses.**  Each float pose `(cx, cy, theta)` of the float search was re-snapped with
`Q = 10^7`, `Dc = 10^8`: `p/q = round(tan(theta/2) * 10^7)/10^7` (reduced), so
`cos = (q^2-p^2)/(q^2+p^2)`, `sin = 2pq/(q^2+p^2)` are rational, and the centre was rounded to
`10^-8` and then **clamped exactly** into `[w/2, 4 - w/2]` with `w/2 = (|cos| + |sin|)/2` rational,
so the closed square lies in the closed container.  Maximum displacement from the float pose:
`1.0e-7` rad in angle, `< 1e-8` in the centre.  The 8 dihedral images are formed exactly and
admissibility is checked exactly for each (the four corners over the common denominator, `0 <= X`
and `X <= 4 D`, same for `Y`).

*`Q = 10^7` matters.*  At `Q = 10^5` (the `3.99` setting) the same 86-pose `PC1` support certifies
only `11.8346`: the LP optimum is degenerate on tight vertices and a `10^-5` rad perturbation of
the geometry breaks the ties, exactly as `DUAL_EXACT.md` records at `3.99`.  At `Q = 10^7` the
snapped support reproduces the float value to 8 decimals (`12.163060631` vs float `12.163060670`),
and the final support reproduces `12.268803861` vs float `12.268805583`.  The angle snap is the
binding approximation, not the centre snap.

**Arrangement vertices, coverage, and why the vertex set suffices.**  Unchanged from
`DUAL_EXACT.md`: square corners plus the exact integer intersection point of every pair of
non-parallel closed edges from different squares, an exact grid prefilter, and the integer
containment test `2|a dx + b dy| <= r D Dk`, `2|a dy - b dx| <= r D Dk`.  `cov` is a finite sum of
indicators of *closed* squares, hence upper semicontinuous; the intersection `P` of the squares
containing a maximiser is a nonempty compact convex polygon on which `cov >= max`, and every vertex
of `P` is a square corner or the meeting point of two non-parallel closed edges.  Masses have the
common denominator `10^9`, so `cov(v) * 8 * 10^9` is an integer and `M` is an exact maximum.

**D4 reduction and the full check.**  The support is exactly D4-invariant, so it suffices to
evaluate the vertices of the fundamental domain `F = {0 <= x <= y <= 2}`; both were run and agree:

| check | vertices | `M` | `L` |
|---|---|---|---|
| `F` (D4 reduction) | 448,330 | `1999999999/2000000000` | `24537607710/1999999999` |
| whole container (`--full`, no reduction) | 3,580,256 | `1999999999/2000000000` | `24537607710/1999999999` |

**Masses (the only heuristic step, not relied upon).**  The exact incidence of every arrangement
vertex in `F` of the snapped support gives the LP rows; `max sum m_k s.t. sum_k (count_{v,k}/8) m_k <= 1`
(HiGHS) chooses the masses; they are rounded **down** to multiples of `10^-9` and re-checked
exactly (`M < 1` already, so no rescaling was needed).  The `check` mode — re-read the file,
re-enumerate, no LP — is what the statement rests on.

## How the measure was found

The float side is `search/packing_dual.py` at `t = 4` with closed semantics, D4-symmetrised,
column generation over poses (the exact Rust verifier's worst placements at `N = 1000` plus the
exact per-angle closed-square sweep pricer, refined by coordinate descent), rows = the violated
arrangement vertices of the current support.  Run `E1` (4 threads, 3436 s, 20 rounds, warm-started
from the `2026-08-26` `PC1` support of `DUAL.md`): LP mass `13.30 -> 12.27`, `M` `1.56 -> 1.000`,
columns `2233 -> 5273 -> 525` as the rows-only polish phase ages them out.  It **converged**: at
round 19 the LP has mass `12.268806`, `M = 1` and **no violated arrangement vertex at all**
(`bad = 0` over 486,112 candidates), and the measure was re-certified in floats over all
`3,867,052` arrangement vertices of the whole container.  Its 311-pose support is what the exact
build snaps.

The exact side is not a single shot: the snapped LP is only as good as its columns, so the poses
were accumulated by an **exact restricted-master column generation** (`search/cover4_cg.py`): start from
a support, add a batch of ~130 candidate poses from the pool (the `E1` per-round supports, the
`2026-08-26` `C1`/`PC1` supports, and the `3.99` exact support rescaled by `4/3.99`), run
`dual_exact.py build` on the union — whose LP sees *all* arrangement vertices of that union as rows —
and keep the poses that come out with positive mass.  Every step's `L` is an exact `Fraction`, and
adding columns can lower `L` (their squares' vertices are extra, genuine `cov <= 1` constraints), so
the best step is kept; each certified value is a theorem regardless of the ones after it.

| step | `L` (exact) |
|---|---|
| `3.99` certificate read at `s = 4` (baseline) | `12.008230783` |
| `PC1` (86 poses) snapped at `Q = 10^5` | `11.834552382` |
| `PC1` snapped at `Q = 10^7` | `12.163060631` |
| `C1` + `PC1` (265 poses), `Q = 10^7` | `12.187432785` |
| exact restricted-master CG over that pool (3 batches) | `12.203964875` |
| + the `3.99` exact support rescaled to 4 (3 batches) | `12.206881232` |
| + the `E1` supports of its first ~9 rounds (3 batches) | `12.236230778` |
| + the `E1` supports up to round 16 (8 batches) | `12.247647654` |
| **`E1` round 19 support (311 poses), `Q = 10^7`** | **`12.268803861`** |

## The measure (294 poses, mass 12.2688)

Concentration: 25% of the mass on one pose, 50% on 4, 75% on 40, 90% on 97, 99% on 226.

| pose (cx, cy, theta) | mass | what it is |
|---|---|---|
| `(0.5, 0.5, 0 deg)` | 3.2527 | the four corner squares, 0.813 each (orbit of size 4, images counted twice) |
| `(0.5, 1.5, 0 deg)` (two poses) | 2.1639 + 0.7138 | axis-aligned squares hugging a wall next to a corner |
| `(1.5, 1.5, 0 deg)` | 0.4084 | the axis-aligned squares of the inner `2x2` |
| `(0.704, 1.732, 39.7 deg)` | 0.2577 | tilted squares leaning on a wall with one corner on it |
| `(1.223, 2.545, 35.8 deg)`, `(1.234, 1.550, 35.5 deg)`, ... | 0.13-0.15 each | the tilted middle ring |

Mass by angle (5 deg bins): `[0,5)` **8.634 (70%)**, `[5,20)` 0.116, `[20,25)` 0.173, `[25,30)`
0.503, `[30,35)` 0.943, `[35,40)` 1.227, `[40,45]` 0.674.  By distance of the canonical centre to
the nearest wall: `< 0.52` 7.308 (60%, 35 poses — corner and wall squares), `0.52-0.8` 1.386 (11%,
56 poses — tilted squares with a corner on a wall), `1.2-1.6` 3.229 (26%, 159 poses — the middle
ring), `1.6-2` 0.346 (3%, 44 poses).  Structurally this is the `3.99` measure of `DUAL.md` with
more room: the rigid frame (corners at 0.81, wall squares) is slightly lighter and the fractional
tilted mixture in `[25, 45] deg` is heavier (3.35 vs 2.33), which is where the extra `0.26` of mass
comes from.

## Numbers and timings (4 processes on a 32-core box shared with other jobs)

| step | count | time |
|---|---|---|
| float column generation `E1` at `t = 4` | 20 rounds, up to 5.3k columns / 10.7k rows | 3436 s |
| float full-container certification of round 19 | 3,867,052 vertices | ~2 min |
| exact restricted-master CG (20 `build`s on 130-300 pose unions) | — | ~50 min |
| `build` on the 311-pose support: vertices in `F` | 482,635 | 30 s |
| `build`: exact incidences | 168,458,726 (vertex, square) pairs | 114 s |
| `build`: LP (floats, HiGHS) | 334,562 rows x 311 cols, 56.7 M nz | 402 s |
| `build`: exact `M`, write support | 294 poses, mass `12.268803855` | < 1 s (total 579 s) |
| `check` in `F` (from the file, no LP) | 448,330 vertices, 154.8 M pairs | 117 s |
| `check --full` (whole container, no D4 reduction) | 3,580,256 vertices, 1,235.8 M pairs | 982 s |

Regression: `python3 search/dual_exact.py build` with no `--t` still writes a
`runs/dual_exact_3.99_support.txt` **byte-identical** to the `2026-08-26` certificate, and
`check` on it reproduces `L = 96065846032/7999999981` exactly.

## What is heuristic

Everything that *chooses* the measure: the float column generation (LP, pricing, row/column
ageing, warm starts — `DUAL.md` describes it), the batching of the exact restricted master, and
the final HiGHS LP that picks the masses.  None of it is relied upon: the certificate is the
support file, and `check` re-derives `mass`, `M` and `L` from it with integers only.  The upper
end `U ~ 12.46` of the bracket is *not* a proof — it is the honest cost (total weight divided by
the dense-scan minimum captured weight) of heuristic covers that are known **not** to certify
(`FAMILY.md` §2 exhibits an exact violating pose for `closed4_best.txt`).  A proved upper bound on
`COVER^closed(4)` would need an exactly verified cover; the cheapest one on record is still
Friedman's / Nagamochi's weight 14.

## Reproduce

```sh
cd verify && cargo build --release && cd ..; mkdir -p runs
# float column generation at t = 4 (warm start from the 2026-08-26 PC1 support), ~60 min
python3 search/packing_dual.py 4 E1 --warm runs/dual_PC1_support.txt --seed-file runs/dual_C1_support.txt \
    --seed-pitch 0.08 --seed-dth 7.5 --row-pitch 0.04 --time 3000 --polish-time 1500 --rounds 80 \
    --N 1000 --cg-want 2500 --threads 4 --final-full
# exact snap + polish + certification of its support (~10 min)
python3 search/dual_exact.py build --t 4 --tag 4 --src runs/dual_E1_support.txt --Q 10000000 --Dc 100000000 --procs 4
# (optional) exact restricted-master column generation over a larger pose pool
python3 search/cover4_cg.py pool runs/pool4.txt runs/dual_E1_support.txt runs/dual_C1_support.txt runs/dual_PC1_support.txt
python3 search/cover4_cg.py pool runs/pool399.txt --scale 1.002506266 runs/dual_exact_3.99_support.txt
python3 search/cover4_cg.py cg runs/pool4.txt search/cover4_exact_support.txt G --batch 130 --rounds 10
# independent exact re-certification from the support file alone (no LP)
python3 search/dual_exact.py check search/cover4_exact_support.txt --t 4 --tag 4 --procs 4            # F,    ~2 min
python3 search/dual_exact.py check search/cover4_exact_support.txt --t 4 --tag 4 --full --procs 4     # full, ~16 min
```

Support file format: `pose p q cx cy mass` with `cx, cy, mass` as fractions and the header
`# t = 4/1` (which `check` reads).  `build` is deterministic given HiGHS; `check` depends only on
the support file and Python integers.
