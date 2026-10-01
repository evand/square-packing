# Green's s(17), s(18) ≥ (40√2+19)/17: what Figure 34 shows, and what survives (2026-09-30)

Inputs and scripts: `search/green_fig34/`.  Short version:

* **Figure 34 (n = 17) shows a 16-point set: four rows of four points.**  Setting the DS7 lemmas'
  constraints to equality reproduces Green's number *exactly*, so that is almost certainly the set
  Green meant (§2).  **But that set is not unavoidable at its side.**  We exhibit a closed unit square
  in `[0, (40√2+19)/17]²` at θ ≈ 30.36° that contains none of the 16 points, checked exactly in
  `Q(√2)` (`green_counterexample.py`).  The cause is a band triangle with one side
  `√((1−e)² + V²) = 1.0110 > 1`, so DS7 Lemma 3 does not apply to it.  The float scan puts the
  largest point-free square at side **1.01090**.
* **Same picture, corrected numbers.**  Re-optimising the same four-row topology gives the set H*,
  critical at `t* = 2√2 + 13/√65 = 2√2 + √65/5 ≈ 4.440879` (§3).  It is 0.0043 below Green's
  value.  The figure's pixels fit H* exactly as well as Green's set.
* **Exact verification (closed squares, closed container, unit weights):**
  `s(17) ≥ 111/25 = 4.44` (and so `s(18) ≥ 4.44`), from a rational 16-point rounding of H*,
  checked by `search/unavoid13_check.py --tri --seg --full`: **VERIFIED**, 0 uncertified boxes.
* **Sanity check of the method (n = 19–20):** Friedman's 18-point set at `6√2 − 4 ≈ 4.485281` is
  numerically critical as drawn.  `s(19) ≥ 1121/250 = 4.484` is **VERIFIED** the same way.

**Recommendation for the site.**  Keep Green 2000 at "claimed" for n = 17–18, and add a note that the
set drawn in DS7 Figure 34 does not support (40√2+19)/17: a unit square fits at that side.  The
same configuration, re-tuned, proves `s(17) ≥ 4.44` here (critical value `2√2 + √65/5 ≈ 4.44088`).
This is historical only, since n = 17 and 18 have much larger preprint bounds (4.66044, 4.679).
For n = 19–20, Friedman's `6√2 − 4` stays "proved"; we re-verified `4.484` by machine.
Green may have used a different set or a sharper triangle lemma, which is why we say "not supported"
rather than "false".  But the unavoidable-set method itself cannot give `s(17) ≥ (40√2+19)/17`
with this topology (§3).

## 1. Sources

* E. Friedman, *Packing Unit Squares in Squares: A Survey and New Results*, Electron. J. Combin.
  Dynamic Survey DS7; maintained HTML <https://erich-friedman.github.io/papers/squares/squares.html>.
  Fetched 2026-09-30, saved as `green_fig34/src/squares.html`
  (sha256 `185a9b460fea6da9b89f172c1c0d0979f1c1829db640e7f07ee06062d80514e8`).
  Theorems 9–10 ("Trevor Green has shown [8]"), then: "Unavoidable sets illustrating some of the
  lower bounds on s(n) are shown in Figure 34."  Table 2 lists `17-18 (40√2+19)/17 ≈ 4.4452,
  Figure 34, Green` and `19-20 6√2−4 ≈ 4.4852, Figure 34, Friedman`.
* Figure 34 = two GIFs, both 171×171, 3-colour (black points, grey proof lines, white):
  * `pic/L17.gif` (caption `s(17) ≥ (40√2+19)/17`), sha256
    `2b6b3af60fe9b629e00824d80c937c8f7de9f9fa354b4b96af03d9ca5a6e8c54`;
  * `pic/L19.gif` (caption `s(19) ≥ 6√2−4`), sha256
    `66891a9fc2c0e9d1094ffe5610be6111ce155f7a204303fcaa1e1292b6e3efae`.
  Saved in `green_fig34/src/`.  The container frame runs from pixel 0 to pixel 170.  Points are 3×3
  black blobs, and their centres were read by connected components (`family.py`: `L17_PX`,
  `L19_PX`).

The figure's captions say only `s(17)` and `s(19)`.  The table's "17-18" and "19-20" follow by
monotonicity, `s(n+1) ≥ s(n)`.  Point counts match: **16 points ⇒ s(17)**, **18 points ⇒ s(19)**.

## 2. The n = 17 picture and Green's number

**What is drawn** (pixels, y downward; scale 170/t ≈ 38.24 px per unit):

    row y=35 : x = 37, 56, 94, 132          (type A: gaps 19, 38, 38)
    row y=68 : x = 37, 75, 113, 132         (type B: gaps 38, 38, 19)
    row y=101: x = 37, 56, 94, 132          (A)
    row y=135: x = 37, 75, 113, 132         (B)

The grey lines draw the hull rectangle and a triangulation of it into near-equilateral triangles.
This is the DS7 Lemma 1/2/3 proof pattern: a corner lemma, a wall-strip lemma, and the lemma that a
triangle with sides ≤ 1 traps a vertex.  The set is invariant under the 180° rotation and not under
the reflections.

**Parametrisation** (`family.green17`): margins `a` (columns) and `b` (rows), in-row spacing `s`,
end gap `e = t − 2a − 2s`, row spacing `V`, middle gap `t − 2b − 2V`:

    A row: x = a, a+e, a+e+s, t−a      B row: x = a, a+s, a+2s, t−a
    rows: y = b (A), b+V (B), t−b−V (A), t−b (B)

**Deriving the bound.**  Make the lemma constraints equalities with `s = 1`:

* DS7 Lemma 2 on the bottom wall (gap 1, distance b): `1 + 2b = 2√2`, so `b = √2 − 1/2`;
* Lemma 2 on the side walls (gap V, distance a): `V + 2a = 2√2`;
* Lemma 3 for a band triangle with horizontal offset e: `e² + V² = 1`;
* both directions add up to t: `t = 2a + e + 2 = 2b + 3V` (equal middle gap).

The unique solution is

    V = (2√2 + 12)/17 ≈ 0.872214,  e = (8√2 − 3)/17 ≈ 0.489042,
    a = (16√2 − 6)/17 ≈ 0.978138,  b = √2 − 1/2 ≈ 0.914214,
    t = (40√2 + 19)/17 ≈ 4.445208     (check: e² + V² = (137 + 152)/289 = 1).

This is exactly Green's value.  The denominator 17 comes from eliminating V = (3+e)/4: `e² + (3+e)²/16 = 1`, i.e. `17e² + 6e − 7 = 0`.  At
38.24 px/unit it draws as a = 37.4, e = 18.7, b = 35.0, V = 33.4 px; the figure has 37, 19, 35 and
33/33/34.  We take this set ("H1") to be Green's.

**H1 fails.**  With in-row spacing 1, a band between an A row and a B row has diagonals with
horizontal offsets e and 1 − e, alternating.  Green's equations bound only the first:
`√(e² + V²) = 1`.  The other is `√((1−e)² + V²) = √(0.2611 + 0.7608) = 1.0110 > 1`, and Lemma 3 does
not apply to such a triangle.  This is a real gap, not just a gap in the proof:

* the float scan (`scan_sets.py H1`, an exhaustive vertex enumeration per angle, 360 angles plus
  refinement) finds a point-free square of side **1.010899** at θ = 30.36° (and 59.64°);
* `green_counterexample.py` checks one pose exactly in `Q(√2)` (sympy): centre `(2.7063, 2.2130)`,
  `tan(θ/2) = 0.27133`, a rational rotation.  The closed unit square lies in `[0, t]²` (wall slacks
  ≥ 1.05).  Every point has `max(|x'|, |y'|) ≥ 0.505403 > 1/2` in the square's frame.
  Output: `REFUTED`.  The square sits in the middle band, in the triangle
  `(a+1, b+V), (a+2, b+V), (a+1+e, t−b−V)`;
* it survives scaling: at side 4.445 (the rounded certificate `certs/fig34_17_green_4.445.txt`) the
  largest point-free square is 1.0108.

`e = 1/2` would make both diagonals equal.  Lemma 3 then forces `V ≤ √3/2`, and `t = 2b + 3V` gives
only `2√2 − 1 + 3√3/2 ≈ 4.42650`.  Moving e towards 1/2 cannot recover Green's number.

## 3. The best set with Figure 34's topology: H*

Let s vary as well.  Lemma 3 on both diagonals forces `e = s/2` and `V = √(1 − s²/4)`.  The walls
give `b = √2 − s/2` and `a = √2 − V/2`.  Balancing the two directions,
`2√2 − s + 3V = 2√2 − V + 5s/2`, gives `V = 7s/8`.  So

    s = 8/√65 ≈ 0.992278,  e = 4/√65 ≈ 0.496139,  V = 7/√65 ≈ 0.868243,
    b = √2 − 4/√65 ≈ 0.918074,  a = √2 − 7/(2√65) ≈ 0.980092,
    t* = 2√2 + 13/√65 = 2√2 + √65/5 ≈ 4.440879.

All triangles have sides ≤ 1, and the walls are exactly tight for Lemma 2.  At 38.28 px/unit this
draws as a = 37.5, e = 19.0, s = 38.0, b = 35.1, V = 33.2 px.  That fits the figure as well as H1:
the pixels cannot tell `s = 1` from `s = 0.992`.

Numerics (`scan_sets.py Hstar`): at `t*` the largest point-free square is **1.000000000**, tight at
θ ≈ 60.26° and 29.75° (the band triangles).  Scaled by 0.999 or 0.9999 it drops to exactly
0.999 or 0.9999, so H* is critical at t*.  An independent Nelder–Mead search of the
4-parameter family (`opt17.py`, true unavoidability, no lemmas) started from the pixel reading and
reached `t* = 4.43935` at `(a, b, s, v) = (0.98004, 0.91707, 0.99170, 0.86840)`, heading for H*.
It is local and not exhaustive, so "H* is the optimum of the 180°-symmetric family" is supported
numerically, not proved.

## 4. The n = 19 picture (Friedman, 6√2 − 4)

Pixels (scale 170/t ≈ 37.90 px/unit): rows at y = 38, 69, 100, 132, i.e. `1, 1+d, 1+2d, 1+3d` with
`d = 2√2 − 2 ≈ 0.8284` and `t = 2 + 3d = 6√2 − 4`.  Rows 0 and 3 are `x = 1, 1+d, 1+2d, 1+3d`.
Row 1 is `x = 1, 1+u, 1+u+d, 1+u+2d, 1+3d` (pixels 38, 59, 90, 121, 132), and row 2 is its 180°
image: 18 points.  Lemma 2 is tight on all four walls (`d + 2·1 = 2√2`).  Lemma 3 needs
`d − √(1−d²)/2 ≤ u ≤ √(1−d²)`, i.e. `0.5483 ≤ u ≤ 0.5601`.  The midpoint
`u = (d + √(1−d²)/2)/2 ≈ 0.55420` matches the drawn 21 px = 0.554, and we use it
(`family19.py`).  Numerically the largest point-free square is **1.000000000** at θ = 45° (the
wall-strip poses), so the set is critical exactly at `6√2 − 4`, as Friedman states.

## 5. Exact checks

The irrational set at side t is scaled by `λ = t'/t < 1` to a rational side t'.  Scaling is
monotone: a unit square in `λ·[0,t]²` is a square of side `1/λ ≥ 1` in the original, and it
contains a concentric unit square, which contains a point.  Coordinates are then rounded to the grid
`1/D` (`mkcert.py`).  The rounding is **not** assumed sound: each certificate (`certificates/FORMAT.md`
format, unit weights) is checked as written.  Checker: `search/unavoid13_check.py`, i.e.
`zeromargin.py` (ADM/CORE/P1/TRI) plus SEG, all exact `Fraction` arithmetic.  Options
`--tri --seg --full --nproc 1` (full domain θ ∈ [0°, 90°], no symmetry assumed), pinned to
core 12.  Statement checked: *every closed unit square in the closed `[0, t']²`, at every angle,
contains (boundary counts) at least one of the points.*  With 16 points this gives `s(17) ≥ t'`
(FORMAT.md's rescaling argument), and `s(18) ≥ s(17)`.

| certificate (`green_fig34/certs/`) | set | t' | D | depth limit | boxes | max depth | leaves ADM / TRI+SEG / EMPTY | uncertified | verdict | time |
|---|---|---|---|---|---|---|---|---|---|---|
| `fig34_17_hstar_4.43.txt` | H*, 16 pts | 443/100 | 10⁴ | 22 | 42,362 | 14 | 20,933 / 1,512 / 14,224 | 0 | **VERIFIED** | 22 s |
| `fig34_17_hstar_4.435.txt` | H*, 16 pts | 887/200 | 10⁴ | 24 | 49,030 | 17 | 23,681 / 1,538 / 14,784 | 0 | **VERIFIED** | 31 s |
| `fig34_17_hstar_4.438.txt` | H*, 16 pts | 2219/500 | 10⁴ | 26 | 59,598 | 20 | 28,483 / 1,539 / 15,265 | 0 | **VERIFIED** | 44 s |
| `fig34_17_hstar_4.44.txt` | H*, 16 pts | **111/25** | 10⁵ | 30 | 87,976 | 26 | 40,511 / 1,543 (410 SEG) / 17,422 | 0 | **VERIFIED** | 80 s |
| `fig34_19_4.48.txt` | Friedman, 18 pts | 112/25 | 10⁴ | 30 | 58,208 | 17 | 29,673 / 388 / 14,531 | 0 | **VERIFIED** | 49 s |
| `fig34_19_4.484.txt` | Friedman, 18 pts | **1121/250** | 5·10⁵ | 32 | 94,470 | 23 | 44,942 / 394 / 17,387 | 0 | **VERIFIED** | 104 s |
| `fig34_19_4.485.txt` | Friedman, 18 pts | 897/200 | 10⁵ | 30 | — | — | 80,868 / 393 / 21,990 | 786 | NOT VERIFIED (depth limit; slack 6·10⁻⁵) | 236 s |

An earlier 4.44 run with D = 10⁴ and depth 22 left 1,938 boxes uncertified, all Lemma-2 wall poses
near θ = 45°.  The finer grid and depth 30 cleared them.

Proved here, by machine: **`s(17) ≥ 111/25 = 4.44`** (hence `s(18) ≥ 4.44`) and
**`s(19) ≥ 1121/250 = 4.484`** (hence `s(20) ≥ 4.484`).

Cross-checks:

* `scan_cert.py` is an independent float scan of the rounded certificates, sharing no code with the
  checker.  Largest point-free square: 0.999809 (`hstar_4.44`, worst at 60.26°) and 0.999715
  (`19_4.484`, at 45°), both < 1.  `green_4.445` gives 1.010817, which FAILS as expected.
* `zmcheck` (Rust) cannot be used: it requires an integer container side, and these problems are
  not scale-invariant.  `zeromargin_stress.py` does not parse SEG witnesses, so it was not run on
  these dumps.

Reproduce (from `s12/`):

    python3 search/green_fig34/mkcert.py 17 111 25 search/green_fig34/certs/fig34_17_hstar_4.44.txt --set hstar -k 4000
    taskset -c 12 python3 search/unavoid13_check.py cert search/green_fig34/certs/fig34_17_hstar_4.44.txt --tri --seg --full --nproc 1 --depth 30
    python3 search/green_fig34/mkcert.py 19 1121 250 search/green_fig34/certs/fig34_19_4.484.txt -k 2000
    taskset -c 12 python3 search/unavoid13_check.py cert search/green_fig34/certs/fig34_19_4.484.txt --tri --seg --full --nproc 1 --depth 32
    python3 search/green_fig34/green_counterexample.py
    python3 search/green_fig34/scan_sets.py            # H1, Hstar, F19 float scans (~1 min)

## 6. Files

`green_fig34/src/` (survey HTML, L17.gif, L19.gif); `family.py` (n = 17 family, pixel data);
`family19.py`; `smax.py` (float largest-empty-square scanner, exploration only);
`scan_sets.py`; `scan_cert.py`; `opt17.py`; `mkcert.py`; `green_counterexample.py`; `certs/`.
