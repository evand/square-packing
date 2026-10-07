# Direct packing column generation: a rigorous `L(3.99) >= 12`

Code: `search/packing_dual.py` (new; one script with its C kernels compiled at run time into the
gitignored `runs/`).  Verifier, certificates and every other file are untouched.

**Outcome.  `L(3.99) = 12.00820 > 12`**: an explicit D4-symmetrised measure on 172 poses of
closed unit squares inside `[0, 3.99]^2`, of total mass `12.00820`, whose coverage is `<= 1`
(`1 + 1.5e-12` with the containment tolerance) at every point of the container — certified by
evaluating the coverage at **all 1,007,508 vertices of the arrangement** of the 1376 support
squares.  Hence `nu_f(s) >= 12.008` for every `s >= 3.99`, so **no weighted unavoidable set of
total weight `< 12` exists at any container side `>= 3.99`**, whatever point set, LP, cell
decomposition or verifier is used.  With the shipped certificate (`COVER(15680/3951) < 12`,
`search/TIGHTEN.md`) the ceiling `s* = sup{s : COVER(s) < 12}` of the pure method satisfies

    3.968616  <=  s*  <  3.99 .

`E(3.99) = L(3.99) − 12 = +0.0082`.  At `3.98` and `3.97` the same procedure reached `11.918` and
`11.807` (below 12; these are lower bounds and say nothing about `s*` there).

## Results

`L` is the certified value (`sum(mu) / M`, `M` = exact max coverage at all arrangement vertices);
"LP mass" is the heuristic value of the column-generation LP.  Wall clock on 32 cores, several
runs sharing the machine; the certified measure of every row is `runs/dual_<TAG>_support.txt`.

| t | colgen run | rounds | wall | max LP mass | polished run (rows only) | rounds | wall | **L(t)** | support (poses / squares) | vertices checked (full container) |
|---|---|---|---|---|---|---|---|---|---|---|
| 3.97 | `D2` (warm from 3.99) | 10 | 1820 s | 11.823 | `PD2` | 9 | 9 s | **11.80666** | 246 / 1968 | 2,085,292 |
| 3.98 | `D1` (warm from 3.99) | 5 | 1786 s | 11.949 | `PD1` | 8 | 7 s | **11.91783** | 193 / 1544 | 1,367,860 |
| **3.99** | `A2` (cold, grid seed) | 16 | 1618 s | 12.014 | `PA2` | 8 | 9 s | **12.00820** | **172 / 1376** | **1,007,508** |
| 3.99 | `B1` (seeded from the 3.9686 certificate) | 15 | 2667 s | 12.031 | – (LP never finished, best saved support gave 11.918) | | | 11.918 | | |
| 4.00 | `C1` (closed convention) | 17 | 1331 s | 12.279 | `PC1` | 8 | 3 s | **12.16306** | 86 / 688 | 299,620 |

Earlier exploratory runs at 3.99: `A1` (12 rounds, no sweep pricer, no polish) `L = 11.83`.
Every `L` in the table was produced by the certification step of the script (`M = 1` to within
`1e-12` after the polish; the `polished` runs were re-checked **without** the D4 reduction over
the whole container, `--final-full`, and from the written support file).  `runs/dual_PA2.log`,
`runs/dual_PA2.json`, `runs/dual_PA2_support.txt` are the deliverables for the best `t`.

## What is rigorous, what is heuristic

**Rigorous (the certification, `Model.certify` / C `vertex_cov`).**  The measure is a finite list
of poses `(cx, cy, theta)` with masses `mu`; every pose is D4-symmetrised (mass `mu/8` on each of
the 8 dihedral images, exactly the `imgs` convention of `nu_f.max_coverage`), so the support is
the list of all images.  Every pose is checked to be admissible (closed unit square inside the
closed container; poses are created with margin `1e-9` and the file's 13-decimal rounding keeps
`>= 9.99e-10`).  For a finite set of **closed** squares the coverage `cov(p) = sum mu_r [p in Q_r]`
is piecewise constant on the arrangement of the squares' edges and upper semicontinuous, so its
maximum over the container is attained at a vertex of the arrangement: a square corner, an
intersection of two edges, or a container corner (a square inside the container meets the
boundary only at its own corners or along an edge whose endpoints are corners).  The script
enumerates every such candidate in floating point and evaluates `cov` with a containment tolerance
`tol = 1e-8` in the square's frame (`|R p − u|_inf <= 1/2 + tol`), which can only over-count.  Two
details make the float enumeration sound: pose angles are snapped to multiples of `1e-5` rad, so
any two edge lines are either exactly parallel (no vertex; their arrangement vertices are corners,
which are enumerated) or meet at an angle `>= 1e-5` rad, and then the 2x2 solve has error
`< 1e-9`, inside `tol`; the on-segment tests use a tolerance `1e-7` (extra candidates are harmless).
So `M = max cov_tol` over the candidates dominates `max_C cov`, `mu/M` is a feasible packing
measure and `L = sum(mu)/M <= nu_f(t)`.  The main loop uses the D4 reduction (candidates in the
fundamental domain `0 <= x <= y <= t/2`, valid because the support and hence the arrangement are
symmetric up to `1e-16`); the reported values were re-certified over the whole container without
it.  Weak duality `nu_f <= COVER` is elementary (`search/CEILING.md`); nothing uses strong duality.
Floating point enters only through the tolerances above and the sums (`sum(mu)` of 172 terms).

**Heuristic (everything that chooses the measure).**  The LP (`max sum mu` s.t. coverage `<= 1` at
a finite set of points; HiGHS via `scipy.optimize.linprog`), the choice of rows (a `0.04` grid of
the fundamental domain, the corners of the seed squares, and then the violated arrangement
vertices returned by the certification, worst first), the pricing (the LP dual `y` is a
D4-symmetric weighted point set; new columns are closed unit squares capturing `< 1` of it, found
by the exact Rust verifier at `N = 1000` with `topk = 8` (`runs/dual_<TAG>_price.txt`), by an
exact per-angle sweep of the closed-square capture (C `sweep_price`, 0.1° steps, `K = 6`
well-separated minima per angle, also against a Wentges-smoothed dual), by perturbations of the
support, all refined by coordinate descent on the exact capture), the column/row ageing, and the
warm starts.  The verifier's `U = total / min` is printed but is uninformative here (20–190): the
dual atoms sit exactly on the boundaries of the support squares (they are arrangement vertices),
so the verifier's `sigma`-shrunk squares miss them while closed unit squares capture them — this
closed-square artefact is also why the raster/verifier witnesses alone gave slow progress and the
exact closed-square sweep pricer was added.

## The measure at `t = 3.99` (`runs/dual_PA2_support.txt`, 172 poses, mass 12.00820)

Concentration: 25% of the mass on one pose, 50% on 4, 75% on 29, 90% on 73, 99% on 140.

| pose (cx, cy, theta) | mu | what it is |
|---|---|---|
| (0.5, 0.5, 0°) | 3.3996 | the four corner squares, 0.85 each (orbit of size 4, images counted twice) |
| (0.5, 1.50025, 0°) | 1.7691 | axis-aligned squares hugging a wall next to a corner, 0.221 on each of the 8 images |
| (0.50607, 2.48385, 0.70°) | 0.5479 | wall squares in the middle of each side, slightly tilted |
| (1.21243, 2.49955, 44.1°) | 0.4463 | 45°-diamonds in the middle ring, 0.056 each |
| (1.50407, 1.50574, 0.46°) | 0.2304 | near-axis squares at the inner corners of the 1.99×1.99 middle |
| (0.70106, 1.835, 37.5°), (0.65328, 1.835, 22.5°), (0.68301, 1.835, 30°) | 0.21, 0.16, 0.09 | tilted squares leaning on a wall with one corner touching it |

Mass by angle (5° bins): `[0,5°)` **6.984 (58%)**, `[5,10°)` 1.019, `[10,15°)` 0.698, `[15,20°)`
0.310, `[20,25°)` 0.672, `[25,30°)` 0.287, `[30,35°)` 0.557, `[35,40°)` 0.638, `[40,45°]` 0.843.
Axis-aligned (`|theta| < 1°`) squares carry 57% of the mass, the rest is spread over the whole
range of angles with a secondary peak at 40–45°.

Where the centres sit (distance of the canonical centre to the nearest wall): `< 0.52`: 6.22 (52%,
11 poses — corners and wall squares); `0.52–0.8`: 2.09 (17%, 86 poses — tilted squares with a
corner on the wall); `1.2–1.6`: 3.38 (28%, 67 poses — the middle ring); `1.6–2`: 0.31 (3%, 8
poses at the very centre).  Coverage by region (integral of `cov`): the four corner unit squares
3.41 (mean coverage 0.857), the four `2×1` edge strips 5.11 (mean 0.642), the inner `1.99²` 3.46
(mean 0.868); a `0.1`-resolution map (`#` = coverage in `[0.95, 1]`, digits = tenths):

```
8888888888344434444444444344438888888888
8888888889444455556556555544449888888888
8888888888444455666776665544448888888888
8888888888445556667777666555448888888888
8888888888455566677777766655548888888888
8888888888555566677777766655558888888888
8888888888455667777777777665548888888888
8888888888456778877777788776548888888888
8888888888467889987777899887648888888888
89888888896999#9###77###9#99969888888898
3444454446778888887777888888776444544443
444455556977778#98888889#877779655554444
444555567987888##889988##888789765555444
44455567898788999#9999#99988789876555444
355566678#888999###88###999888#876665553
45566678998##999##9889##999##89987666554
456666789#89#9####8888####9#98#987666654
456677778#888####778877####888#877776654
466777777#7889#9877887789#9887#777777664
4577777777789988888778888899877777777754
4577777777789988888778888899877777777754
466777777#7889#9877887789#9887#777777664
456677778#888####778877####888#877776654
456666789#89#9####8888####9#98#987666654
45566678998##999##9889##999##89987666554
355566678#888999###88###999888#876665553
44455567898788999#9999#99988789876555444
444555567987888##889988##888789765555444
444455556977778#98888889#877779655554444
3444454446778888887777888888776444544443
89888888896999#9###77###9#99969888888898
8888888888467889987777899887648888888888
8888888888456778877777788776548888888888
8888888888455667777777777665548888888888
8888888888555566677777766655558888888888
8888888888455566677777766655548888888888
8888888888445556667777666555448888888888
8888888888444455666776665544448888888888
8888888889444455556556555544449888888888
8888888888344434444444444344438888888888
```

So the measure is **neither a few configurations nor fully diffuse**: half of the mass is a rigid
frame (four corners at 0.85, eight wall squares at 0.22–0.27), and the other half is a genuinely
fractional, many-angle mixture — tilted squares with one corner on a wall (mass 2.1 over 86 poses)
interleaved with a ring of diamonds and near-axis squares in the middle (3.4 over 67 poses).  The
edge strips are the least covered region (0.64), and the interior diamonds reach down to `y ≈ 1`
exactly where the wall squares end — the strip squares cannot be heavier because the tilted
interior already saturates the coverage along their top edges.  This mixture is the object a case
analysis has to kill: any argument beyond the LP must exploit that a real packing cannot put
fractional mass on interleaved tilted squares.

## How the numbers were obtained, and what limits them

* Cold start `A2` (`seed 0.08 / 7.5°`, `2054` poses): the LP mass rose 11.76 → 12.01 in 16 rounds
  of ~1 min, then LPs of `5k × 13k` (25M nonzeros) took 5 min each.  `M` stayed at 1.01–1.05
  during column generation, so the in-loop `L` never exceeded 11.90.
* Polish `PA2`: the saved support (396 poses) alone, rows only — a small LP (`≤ 400` columns) and
  the finite vertex set converge in 8 rounds / 9 s to `M = 1` exactly with mass 12.0082.  This
  rows-only step on the support is what turns an LP mass of 12.01 into a certified 12.008; the
  vertex rows (instead of grid points or dilated cells) are what makes `M = 1` reachable at all —
  with dilated cells the number of leaves needed is `~ 0.7 · (total edge length of the support) / h`,
  i.e. `10^8` for a support like this, which is why `CEILING.md`'s `M` stayed at 1.02–1.14.
* Warm starts at 3.98 / 3.97 (`D1`, `D2`: the 3.99 support with centres rescaled, plus the
  certificate-derived seed): LP masses 11.95 / 11.82 after 5–10 rounds of 3–5 min (the warm start
  brings ~20k rows, 30–45M nonzeros), certified 11.918 / 11.807 after polishing.  These are limited
  by wall clock, not by certification slack (`M = 1`); the LP masses had not stabilised.
* Seeding from the shipped 3.9686 certificate (`B1`): the certificate scaled to 3.99 has minimum
  capture 0.22 (a cliff, cf. `TIGHTEN.md`), so its 48k tight/hole placements make a good column
  seed (mass 12.04 at round 0, 12.03 after the vertex rows are in), but the LPs grew to 35M
  nonzeros and the last one did not finish in the budget.
* At `t = 4` (closed convention) the same machinery certifies `nu_f^closed(4) >= 12.163`
  (`CEILING.md` had 10.23 rigorous, 11.9–14.8 heuristic).

Limiting factor for going lower than 3.99: LP size.  The closed-square convention makes the LP
"dodge" point constraints with slightly shifted copies of support squares, so the row set (the
arrangement vertices) and the column set both grow by thousands per round and HiGHS time grows
super-linearly; the polish step then fixes `M` but cannot add mass.  A better implementation
would keep the LP small (aggressive column pruning to the support plus a few hundred priced
columns, rows restricted to vertices of the current support) and polish every few rounds.

## Reproduce

```sh
cd verify && cargo build --release && cd ..; mkdir -p runs
# column generation at 3.99 (cold start), then rows-only polish of its support
python3 search/packing_dual.py 3.99 A2 --seed-pitch 0.08 --seed-dth 7.5 --row-pitch 0.04 --time 2400 --polish-time 600 --rounds 80 --N 1000 --cg-want 2500
python3 search/packing_dual.py 3.99 PA2 --seed-file runs/dual_A2_support.txt --seed-pitch 10 --seed-dth 45 --row-pitch 0.04 --polish-after 0 --time 0 --polish-time 1500 --rounds 300 --final-full
# 3.98 / 3.97 warm-started from the certified 3.99 measure (runs/seed_B1.txt = tight placements of the scaled 3.9686 certificate)
python3 search/packing_dual.py 3.98 D1 --seed-pitch 0.12 --seed-dth 15 --row-pitch 0.04 --seed-file runs/seed_B1.txt --warm runs/dual_PA2_support.txt --time 1500 --polish-time 300 --rounds 80 --N 1000 --cg-want 2000
python3 search/packing_dual.py 3.98 PD1 --seed-file runs/dual_D1_support.txt --seed-pitch 10 --seed-dth 45 --row-pitch 0.04 --polish-after 0 --time 0 --polish-time 600 --rounds 300 --final-full
python3 search/packing_dual.py analyse runs/dual_PA2_support.txt        # structure of a measure
```

HiGHS is deterministic for a fixed input, but the column-generation runs stop on the wall clock
and share the machine, so they are not bit-reproducible; what is asserted is checked from the
support file (`--final-full` re-certifies it over the whole container).
