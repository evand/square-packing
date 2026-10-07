# The cover LP at the container `[0,4]^2` itself, closed semantics (heuristic)

Code: `search/closed4.py`.  Outputs: `runs/closed4_*` (logs, JSON histories, certificate-format
point sets, PNGs).  One ~2.5-hour session on 8 of 32 shared cores (load ≈ 40 throughout).

**Everything in this note is a heuristic LP with sampled rows.  Nothing here is a certificate,
and nothing here is a bound on `COVER(4)` in either direction** — see "What is heuristic" at the
end.  The one exception is the rigorous lower bound `L` in §7, which is the float-checked
`nu_f.py lower` machinery of `search/CEILING.md` applied to the point set found here.

## 1. The question

`COVER^closed(s)` = least total weight of a weighted point set in `[0,s]^2` such that every
**closed** unit square contained in `[0,s]^2`, at every angle, captures weight ≥ 1, a point on the
boundary of the square counting (the convention of `certificates/FORMAT.md` and `verify/`).
At `s = 4`:

* Friedman's 14 unit-weight points (DS7 Thm 4) and Nagamochi's weighted set (N05 §3, a = b = 4)
  both cost 14; nothing cheaper is in the literature (`notes/proof-anatomy.md` §7.1).
* A cover of weight `< 12` would prove `s(12) = 4` outright: a packing of 12 unit squares in a
  container of side `t < 4`, rescaled to `[0,4]^2`, gives 12 pairwise disjoint closed unit squares,
  which capture ≥ 12 with no point counted twice.
* `search/CEILING.md` recorded that the exact-mode LP at `s = 4` "stalled at a degenerate cover:
  its optimum wants points exactly on the grid lines, which leaf-centred column generation never
  proposes".  This note runs that LP with the right columns.

**Answer (heuristic): the sampled LP does not go below 12; the dual side agrees (packing mass
12.2 on a 0.01 grid, rigorous `ν_f ≥ 10.68`).  It sits at 12.30–12.42 and climbs
as rows are added; the cheapest point set found costs 12.417 (LP value) and 12.51 honestly
(total / min captured weight over a dense scan).  Its weight sits almost entirely on the six grid
lines `x, y ∈ {1,2,3}`, as CEILING.md predicted.**

## 2. Set-up

*Semantics.*  `w(Q)` is upper semicontinuous in the pose of the closed square `Q`, so the cover
condition for every closed unit square with centre in the **open** admissible box already implies
it for the wall-touching squares.  Rows are therefore placed in the open box, offset `1e-7` from
its boundary (the near-wall limit is sampled explicitly by "band" rows along each wall).  A
consequence: a point *on* a container wall is captured by no row, so it is useless in closed
semantics (a wall-touching square is dominated by its interior neighbours, which miss the wall).
Wall points were offered as columns anyway; every run gave them weight 0.

*Columns.*  D4 orbits (as `tighten.py`) of: the lattice of pitch 0.05 on `[0,4]^2`, which contains
the integer and half-integer lines exactly (863 orbits / 6573 atoms); the points of Friedman's
14, Bentz's 16 and Nagamochi's 12 (8 segment endpoints + 4 points); and points added by pricing
on the D4-symmetrised dual, evaluated on a 0.01 grid (which again contains the grid lines exactly)
and refined on the 0.001 grid (`price`, copied from `tighten.py`).  All coordinates live on the
`D = 1000` grid.

*Rows.*  Closed unit squares, containment test `|R(-θ)(p - c)|_∞ ≤ 1/2 + 1e-9` (the boundary
counts, with the margin in the direction that counts a point as captured).  Centres on a lattice
of pitch 0.02 in the rotated frame of each angle, restricted to the open admissible box, plus the
wall bands; angles (degrees) `0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5, 1, 2, …, 45` and
`45 − {0.5, 0.1, 0.03, 0.01, 0.003, 0.001}` (59 angles; `x ↔ y` maps `45 + δ` to `45 − δ`); rows
reduced by the C4 rotations that fix the angle.  Cutting planes each round: (a) the full lattice
scan, worst violated pose per 0.15-block per angle (≤ 60 per angle); (b) a local search
(`polish`: random perturbations of centre and angle with radii shrinking from 0.03 / 1° to
`0.03·0.3⁶`, plus jumps to the exact angles `0, 45°, 1e-6, 1e-5, 1e-4, 1e-3 rad` and their
`45°` mirrors) seeded from the 48 worst lattice poses, the 48 tightest LP rows and, every third
round, 48 random poses; every visited pose with captured weight `< 1 − 1e-9` is added (≤ 240
per round after thinning).  Pricing every round from round 2, ≤ 150 new orbits per round.

## 3. Sanity checks (all passed)

| check | result |
|---|---|
| (i) Friedman's 14 points, unit weights, **unsymmetrised**, rows over the full range `[0°, 90°]` (117 angles, pitch 0.02, wall bands) | min captured = **1.0000**, 0 violated poses; local search from 64 random poses: min 1.0000 |
| (i') Bentz's 16 points, unit weights, same rows | min captured = 1.0000, 0 violated; polished min 1.0000 (so Bentz's box-unavoidable set is unavoidable for closed unit squares too) |
| (ii) LP over Friedman's 14 points alone (columns = the 14 points, rows as in (i), cut loop to convergence) | **14.000000**, all weights exactly 1 |
| (i'') Nagamochi's weighted set with the segments and the area discretised as atoms | pitch 0.01: total 14.000, min captured 0.953 (cost 14.70); pitch 0.005: min 0.972 (cost 14.41).  The deficit halves with the pitch, i.e. it is the discretisation of the segments/area (`0.5·pitch` per cut segment end), consistent with N05 Lemma 1 (`σ(S) > 1` for the continuous set) |
| axis-aligned rows only (angle 0) | LP = **9.000000**, converged in 2 rounds, support = exactly the nine integer points `{1,2,3}²`, each of weight 1 (the 3×3 packing of disjoint closed axis-aligned unit squares shows 9 is optimal) |
| (iii) the same code at `s = 3.99` | LP **12.12** after 15 rounds (45 min, time-limited); stress min 0.955, cost 12.70.  Consistent with the shipped certificates and `CEILING.md` (12.19 @ 3.97, 12.48 @ 3.98 for the cover LP over a fixed point set; the column generation here reaches lower).  At 3.99 the weight is *not* on the grid lines (4.1 of 12.1 on integer lines, 7.9 elsewhere): the degenerate line structure is specific to `s = 4` |

The transcription of the literature sets is therefore correct and the row generator agrees with
the published unavoidability statements.

## 4. The LP at `s = 4`: value per round

Main run: `python3 search/closed4.py run --s 4 --tag s4 --nproc 5 --time 4200 --perang 60 --block 0.15`
(73 min, 25 rounds, stopped by the time limit; each LP solve took 3–6 min at the end).

| round | LP value | rows | orbits / atoms | support (orbits) | lattice min | lattice viol. | polish min | rows added | cols added | dual cov. max | t (s) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 9.2706 | 7320 | 863 / 6573 | 18 | 0.0000 | 423458 | 0.0000 | 3780 | 0 | – | 6 |
| 1 | 11.8410 | 10250 | 863 / 6573 | 37 | 0.5037 | 36029 | 0.5037 | 2930 | 0 | – | 20 |
| 2 | 12.4065 | 12421 | 972 / 7341 | 63 | 0.8789 | 9072 | 0.8789 | 2171 | 109 | 1.2314 | 34 |
| 3 | 12.3522 | 13587 | 1087 / 8157 | 82 | 0.9059 | 2227 | 0.8457 | 1166 | 115 | 1.1726 | 54 |
| 4 | 12.3513 | 14656 | 1199 / 8933 | 130 | 0.9583 | 1683 | 0.9095 | 1069 | 112 | 1.2266 | 92 |
| 5 | 12.3044 | 15600 | 1311 / 9737 | 150 | 0.9208 | 1455 | 0.8997 | 944 | 112 | 1.1840 | 122 |
| 6 | 12.3113 | 16368 | 1402 / 10385 | 187 | 0.9618 | 1055 | 0.9350 | 768 | 91 | 1.2151 | 165 |
| 7 | 12.2754 | 16845 | 1524 / 11225 | 227 | 0.9685 | 396 | 0.9604 | 477 | 122 | 1.0740 | 222 |
| 8 | 12.2646 | 17353 | 1626 / 11961 | 277 | 0.9742 | 445 | 0.9450 | 508 | 102 | 1.0404 | 307 |
| 9 | 12.2704 | 17792 | 1714 / 12601 | 316 | 0.9498 | 356 | 0.9429 | 439 | 88 | 1.0872 | 408 |
| 10 | 12.2770 | 18234 | 1773 / 13013 | 324 | 0.9737 | 327 | 0.9472 | 442 | 59 | 1.0386 | 546 |
| 11 | 12.2784 | 18506 | 1835 / 13445 | 344 | 0.9936 | 50 | 0.9540 | 272 | 62 | 1.0315 | 726 |
| 12 | 12.2853 | 18811 | 1892 / 13829 | 357 | 0.9731 | 99 | 0.9731 | 305 | 57 | 1.0249 | 914 |
| 13 | 12.2845 | 19116 | 1938 / 14157 | 400 | 0.9910 | 99 | 0.9708 | 305 | 46 | 1.0186 | 1076 |
| 14 | 12.2866 | 19393 | 1984 / 14489 | 403 | 0.9959 | 59 | 0.9576 | 277 | 46 | 1.0305 | 1255 |
| 15 | 12.2905 | 19683 | 2019 / 14717 | 427 | 0.9958 | 63 | 0.9739 | 290 | 35 | 1.0191 | 1486 |
| 16 | 12.2899 | 19947 | 2060 / 14993 | 436 | 0.9959 | 33 | 0.9565 | 264 | 41 | 1.0184 | 1775 |
| 17 | 12.2941 | 20212 | 2099 / 15273 | 458 | 0.9963 | 32 | 0.9666 | 265 | 39 | 1.0196 | 2116 |
| 18 | 12.2951 | 20469 | 2127 / 15477 | 481 | 0.9949 | 29 | 0.9621 | 257 | 28 | 1.0153 | 2444 |
| 19 | 12.2984 | 20738 | 2162 / 15729 | 484 | 0.9862 | 49 | 0.9618 | 269 | 35 | 1.0240 | 2783 |
| 20 | 12.2992 | 20987 | 2191 / 15937 | 475 | 0.9972 | 14 | 0.9862 | 249 | 29 | 1.0175 | 3229 |
| 21 | 12.3009 | 21240 | 2209 / 16073 | 495 | 0.9954 | 23 | 0.9913 | 253 | 18 | 1.0119 | 3559 |
| 22 | 12.3018 | 21492 | 2232 / 16233 | 480 | 0.9987 | 16 | 0.9759 | 252 | 23 | 1.0261 | 3928 |
| 23 | 12.3022 | 21755 | 2256 / 16409 | 500 | 0.9895 | 32 | 0.9829 | 263 | 24 | 1.0163 | 4155 |
| 24 | 12.3044 | 22013 | 2280 / 16581 | 485 | 0.9956 | 23 | 0.9659 | 258 | 24 | 1.0175 | 4380 |
| final | **12.3094** | 22013 | 2280 / 16581 | | | | | | | | |

("lattice min" = least captured weight on the 0.02 row lattice under that round's weights;
"polish min" = least found by the local search; "dual cov. max" = largest symmetrised closed
coverage of the dual measure at any 0.001-grid point, i.e. `1 −` the most negative reduced cost
per atom.  The value went *up* from round 8 on because the cuts outweighed the new columns; the
pricing gain shrank to 1.2–1.7 % per atom.)

Stress test of the round-24 set (3049 points, total 12.3096; scan at pitch 0.005, 254 angles
including 200 random ones in `[0°, 90°)`, then local polish): min captured **0.9514** at centre
`(0.5056, 2.4953)`, angle `0.64°` — a square almost touching the left wall, tilted by half a
degree, a region the 0.02 / 0.5° row lattice undersamples.  Honest cost **12.94**.

**Refinement** on the fixed support (`--from-cert runs/closed4_s4_best.txt --no-colgen --pitch 0.01
--deg-step 0.5`, 104 angles, 480 orbits / 3049 atoms, 69 rounds in 25 min, time-limited):

| round | 0 | 1 | 5 | 10 | 20 | 30 | 40 | 50 | 60 | 68 | final |
|---|---|---|---|---|---|---|---|---|---|---|---|
| LP | 12.3096 | 12.3135 | 12.3307 | 12.3839 | 12.3999 | 12.4111 | 12.4147 | 12.4159 | 12.4168 | 12.4172 | **12.4174** |
| rows | 14.4k | 15.0k | 17.2k | 23.0k | 27.0k | 30.1k | 30.9k | 31.4k | 32.5k | 33.1k | 33.1k |
| lattice min | 0.951 | 0.930 | 0.982 | 0.923 | 0.979 | 0.993 | 0.995 | 1.000 | 0.996 | 0.987 | |

Stress test of the refined set (`runs/closed4_best.txt`, 1972 points, total 12.4175): lattice min
0.99267, polished min **0.99267** at `(3.3965, 1.4235)`, angle `76.41°` (equivalently `13.6°`
by symmetry; the next five worst are all at 14–15° / 75–76° near a wall).  **Honest cost
12.4175 / 0.99267 = 12.509.**

## 5. Where the weight sits (`runs/closed4_best.txt`, total 12.4174; picture `runs/closed4_best.png`)

| location | weight | share |
|---|---|---|
| on the six grid lines `x ∈ {1,2,3}` or `y ∈ {1,2,3}` (exactly, on the 0.001 grid) | **11.342** | 91.3 % |
| — of which at the nine integer points `{1,2,3}²` | 0.135 | 1.1 % |
| — per line: `x = 1`: 1.958 over `y ∈ [0.41, 3.59]`; `x = 2`: 1.822 over `[0.47, 3.53]`; `x = 3`: 1.958; the `y`-lines the same by symmetry | | |
| within 0.05 of a grid line but not on it (0.001–0.01 offsets, the LP's way of splitting the two sides of a line) | 0.484 | 3.9 % |
| interior arcs (four rings of short arcs around each of `(1.5,1.5)`, `(1.5,2.5)`, …, at distance 0.1–0.5 from the lines) | 0.592 | 4.8 % |
| on the half-integer lines | **0** | |
| on the container walls | **0** | |

The heaviest single atoms (0.175 each) are the eight points `(0.73, 1), (1, 0.73)` and images —
on the grid lines, 0.27 from the wall-side end, where the line measure has to be concentrated
because a square touching the wall at a small angle cuts a line only near there.  The profile
along `x = 1` in 0.25-bins of `y` is
`0.03, 0.25, 0.21, 0.10, 0.15, 0.13, 0.12, 0.12, 0.13, 0.15, 0.07, 0.24, 0.25, 0.03`
(from `y = 0.25` to `3.75`): dense near both ends, roughly uniform at ≈ 0.5 per unit length in
the middle — Nagamochi's segments `L_i` carry exactly 0.5 per unit length.  So the optimum is a
"`#`" of line measure on all six lines (Nagamochi puts it on four, plus area on `[1,3]²`) plus a
small interior correction; 12.4 versus his 14.

By origin of the column: 10.7 of the 12.31 (main run) came from priced columns, 1.5 from the
0.05 lattice, 0.09 from Bentz's points, nothing from Friedman's or Nagamochi's points as such.
The pricing did exactly what CEILING.md said leaf-centred generation could not: propose points
*on* the lines at 0.001 resolution.

## 6. Reading the numbers

* A row-sampled LP is a **relaxation** of the cover problem over its own columns, so 12.417 is
  (up to `1e-9` float tolerance) a *lower* bound on the cheapest cover supported on those 1972
  points, and the LP value can only rise as rows are added — which is what the two runs show
  (12.30 → 12.42 with denser rows, still rising slowly at the time limit).
* Restricting columns is the opposite: the true `COVER(4)` could in principle be lower with points
  not offered.  The pricing gain at the end of the main run was 1.2–1.7 % per atom on the 0.001
  grid, and the value had been climbing, not falling, for 16 rounds, so a drop from 12.3 to below
  12 from further columns is not indicated by anything here — but it is not excluded by anything
  here either.
* The honest cost of the best explicit set is **12.51**; a cover of weight `< 12` at `s = 4` would
  need a set 4 % cheaper than anything the LP found and 0.4 units below the LP's own relaxed value.
* Consistency with the rest of the repository: `COVER^closed` is non-decreasing in `s`, the
  shipped certificate gives 11.974 at 3.9686, `TIGHTEN.md` finds no cover below 12 at 3.9696, this
  note finds 12.12 (heuristic) at 3.99 and 12.3–12.4 at 4.  The crossing of 12 stays where
  `CEILING.md` put it, in `3.94 ≲ s* ≲ 3.97`.

## 7. Rigorous lower bound from the dual

`python3 search/closed4.py lower runs/closed4_best.txt --tag lower` runs `nu_f.lower_from_cert`
(cover LP over the file's points with a broad row set → dual support + 3000 tightest placements →
packing LP against the 0.01 grid → `max_coverage` certification by subdivision, worst leaves fed
back for up to 7 rounds; see `CEILING.md` for why every round gives a valid `L ≤ ν_f(s)`).

Result (11 min, 4 processes):

| step | value |
|---|---|
| cover LP over the 1972 points with `nu_f`'s own row set (24541 rows: lattice `(η, δθ) = (0.01, 0.01 rad)` poses with `w(Q) < 1.05` plus the verifier's worst placements) | 12.251 (our finer row set gives 12.417 on the same points — fewer rows, lower value, as expected) |
| packing support | dual support 140 placements + 3000 tightest |
| packing mass on the 0.01 grid (heuristic `ν_f`) | 12.79 → 12.42 → 12.39 → 12.29 → 12.23 → 12.20 → 12.18 over the 7 rounds |
| certified max coverage `M` | 1.85 → 1.38 → 1.25 → 1.25 → **1.145** → 1.24 → 1.21 |
| **rigorous `L = mass / M`** | best **10.68** (round 4) |

So `ν_f^closed(4) ≥ 10.68` rigorously (float-checked; CEILING.md had 10.24 from a worse point
set), and the heuristic packing mass is 12.2 — the dual side agrees with the cover side that the
value is a little above 12.  The certification slack `M − 1 ≈ 15 %` is the same compute-limited
looseness documented in CEILING.md ("closing the rigorous interval from above needs a much richer
packing support"); it does **not** push `L` to 12, so the statement "no closed-semantics cover of
weight `< 12` exists at `s = 4`" remains unproved.  It is, however, what every number in this
note points to.

## 8. What is heuristic here

* Rows are a finite sample (pitch 0.02, then 0.01, in the rotated frame; 59 then 104 angles; local
  search).  The stress scan (pitch 0.005, 254 angles) found poses 0.7 % below 1 for the refined set
  and 5 % below for the un-refined one, so the *true* minimum captured weight of either set is at
  most what is reported and probably a little lower.  The exact verifier (`verify/`) would give
  the true minimum; it was not run because the totals are above 12 and nothing rides on them.
* Columns are a finite set (0.05 lattice, literature points, pricing on 0.01 → 0.001 grids).
* Floats throughout (`1e-9` margins).  The LP values are HiGHS optima of the sampled LPs.
* Only §7 is rigorous (modulo float margins), and only as a lower bound on `ν_f`.
* Two runs were stopped by wall-clock limits, not convergence; the refinement was within `1e-4`
  per round of its limit.

## 9. Files

| file | what |
|---|---|
| `search/closed4.py` | the code (`run`, `sanity`, `axis`, `stress`, `nagamochi`, `lower`; `--help`) |
| `runs/closed4_best.txt` / `.json` / `.png` | the refined set: certificate format, `s = 4`, `D = 1000`, `W = 1e7` (weights rounded up), 1972 points, total 12.4175; JSON with the point list, breakdown, history and stress result; picture |
| `runs/closed4_s4.log` / `.json` / `_best.*` | the main run (round table above; un-refined set, 3049 points) |
| `runs/closed4_refine.log` / `.json` / `_best.*` | the refinement (same content as `closed4_best.*`) |
| `runs/closed4_s399*` | the `s = 3.99` check |
| `runs/closed4_axis*` | the axis-aligned reference (9 points) |
| `runs/closed4_sanity.log`, `runs/closed4_nagamochi*.log` | sanity checks |
| `runs/closed4_lower.log` / `runs/closed4_lower_lower.json` | §7 |
