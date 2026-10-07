# Occupancy census at side 4: 12 unit squares (2026-10-06)

Everything here is [measured]: float multistart optimisation, so each margin is the best found, a lower bound on
the class's true sup.  None of it is a proof.  Negative numbers mean "best found", not "impossible".

## Setup

* Side `L = 4`, strip width `h ∈ {1.1213, 1.15, 1.18, 1.2071}`.  Centre regions are closed boxes: `K_i` is corner
  `i` (`[1/2, h]²` and its images), `E_w` is within `h` of wall `w` only, and `Z = [h, 4−h]²`.  The class is
  `(K0..K3, E_B,E_R,E_T,E_L, z)` up to D4, labelled `k e z K.... E.... S....`, where `S` is the per-wall strip count
  `E_w + K_w + K_{w+1}`.
* Enumerated: every class with `K_i <= 1` and every strip count `<= 3`, which gives 235 D4 classes.  Note that
  `ΣE + 2ΣK <= 12` already forces `z >= k` for `N = 12`, so every enumerated class has `z >= k`.
* Margin: the min over pairs of the SAT gap (edge normals of both squares) and over squares of the wall distance.
  Equivalently, the uniform gap.
* Optimiser (`opt.py`): the epigraph form.  It maximises `δ` over the poses plus one free separating line
  `(φ_ij, c_ij)` per pair, with the 8 vertex inequalities per pair, the wall constraints
  `x ± (cosθ+sinθ)/2` vs `δ`, and the class as bound constraints on the centres.  It uses SLSQP with analytic
  Jacobians (checked by finite differences to 8e-10), and gets up to 4 rounds, each re-seeding the lines from the
  best SAT axis.  The reported margin is recomputed independently by `sat_margin`.  `export.py` rechecks all 989
  witness files (margin and box membership): 0 mismatches.
* Seeds (`census.py`) are mixed:
  * uniform random poses in the class boxes;
  * jittered 12-subsets of the 4×4 tiling fitted to the boxes;
  * a central 45° quincunx or ~26° pinwheel with the other squares near axis-parallel;
  * near-axis random poses;
  * a library of 1,041 zero-margin free optima (stage A) mapped into the class by D4 plus a Hungarian assignment.

  `intensify.py` adds basin hops from the best point.
* Sanity check: on the 5-in-Z subproblem (box side `3 − √2`) it finds **+0.12132** (a 45° quincunx).  The earlier
  Nelder–Mead witness `scratch/zwit.py` had +0.0156.  So `Z <= 4` is false locally with a lot of room.
* Effort: 1,200 free solves, 143,854 class solves in stage B, and 44,600 intensification solves, about 190k SLSQP
  solves (~0.3 s each).
  * Stage B stopped a class early once it reached `δ >= −1e-9`, after at least 6 restarts.
  * Otherwise it ran 200 restarts for `z <= 7` and 60 for `z >= 8`.
  * Intensification added 600 restarts to every near-0 negative class (`z = 5` with best `> −0.08`, and
    `z <= 4` classes with some `E_w = 3` and best `> −0.07`).
  * The two best `z = 5` classes got 2,000 extra restarts (4,800 at `h = 1.1213`).

Files:
* `runs/` is gitignored (repo convention): data and witnesses are local only, on evand's machine; `census.py`/`intensify.py` regenerate them (~190k solves).
* Data: `runs/classes.jsonl` (stage B, one line per (h, class)), `runs/intensify.jsonl`, `runs/free.jsonl`.
* Witnesses: `runs/witnesses/h<h>_<class digits>[_int].json`, holding poses `(x, y, θ rad)`.
* Full per-class tables: `runs/table_kz.md`, `runs/table_z5.md`, `runs/table_z0to4.md`.  `summarize.py`
  regenerates them.

## 1. Main table: best margin over the classes with given (k, z)  [measured]

The number in parentheses is how many classes in that group reach 0 (`>= −1e-9`).  Rows are sorted roughly by
`z − k`.

| k | z | #classes | h=1.1213 | h=1.15 | h=1.18 | h=1.2071 |
|---|---|---|---|---|---|---|
| 0–4 | 8–12 | 54 | <= −0.156 | <= −0.177 | <= −0.193 | <= −0.207 |
| 0–4 | 6–7 | 87 | <= −0.0905 | <= −0.0955 | <= −0.104 | <= −0.114 |
| **4** | **5** | 1 | **−0.0363** | −0.0410 | −0.0410 | −0.0410 |
| **3** | **5** | 5 | **−0.0339** | −0.0411 | −0.0413 | −0.0418 |
| 2 | 5 | 13 | −0.0691 | −0.0731 | −0.0773 | −0.0816 |
| 1 | 5 | 14 | −0.0753 | −0.0901 | −0.0969 | −0.1044 |
| 0 | 5 | 7 | −0.0945 | −0.1061 | −0.1154 | −0.1055 |
| 4 | 4 | 1 | **0** (1) | 0 (1) | 0 (1) | 0 (1) |
| 3 | 4 | 2 | 0 (2) | 0 (2) | 0 (2) | 0 (2) |
| 2 | 4 | 10 | 0 (7) | 0 (7) | 0 (7) | 0 (7) |
| 1 | 4 | 9 | 0 (2) | 0 (2) | 0 (2) | 0 (2) |
| 0 | 4 | 7 | 0 (1) | 0 (1) | 0 (1) | 0 (1) |
| 3 | 3 | 1 | 0 (1) | 0 (1) | 0 (1) | 0 (1) |
| 2 | 3 | 4 | 0 (2) | 0 (2) | 0 (2) | 0 (2) |
| 1 | 3 | 6 | 0 (1) | 0 (1) | 0 (1) | 0 (1) |
| 0 | 3 | 4 | −0.0708 | −0.0951 | −0.1039 | −0.1039 |
| 2 | 2 | 2 | 0 (1) | 0 (1) | 0 (1) | 0 (1) |
| 1 | 2 | 2 | −0.0529 | −0.0526 | −0.0506 | −0.0503 |
| 0 | 2 | 3 | −0.0803 | −0.0960 | −0.1084 | −0.1107 |
| 1 | 1 | 1 | −0.0819 | −0.0776 | −0.0730 | −0.0688 |
| 0 | 1 | 1 | −0.1029 | −0.1237 | −0.1205 | −0.1325 |
| 0 | 0 | 1 | −0.0990 | −0.0923 | −0.0966 | −0.0790 |

The full per-(k, z) rows are in `runs/table_kz.md`.  The set of zero classes is the same 18 classes at all four h.
No class has a positive margin, as expected.

## 2. Priority 1: crowded centre, z >= 5  [measured]

* **Best margin with `z >= 5`: −0.0339 at `h = 1.1213`** (class `k3 e4 z5 K1110 E1111 S3322`; 4,800 restarts plus
  hops; the best value was hit twice).  It is −0.0363 for `k4 e3 z5`.  At `h >= 1.15`, both classes give
  −0.041 (2,200 restarts each).  Every other `z = 5` class is `<= −0.053`, and `z >= 6` is `<= −0.090`.
* The free (unconstrained) search never produced a zero-margin point with `z >= 5`: 0 of 1,041 zero optima.
* What the optima look like:
  * At `h = 1.1213` (`runs/witnesses/h1.1213_111011115_int.json`): five Z squares tilted 18° form a pinwheel, with
    two Z centres pinned on the box boundary at coordinate `h`.  The margin therefore depends on `h` here.
  * At `h = 1.2071` (`h1.2071_111111105_int.json`): five Z squares at 42.5° form a near-quincunx, three E squares
    at ~41–42°, and the corner squares are tilted 4.9°.
* Compare the isolated centre: five squares fit in `Z` with gap **+0.121** (45° quincunx).  So filling the boundary
  takes away ≥ 0.155 of that, and the result is still ≥ 0.034 short.  **A crowded centre does cost boundary
  squares, with room ≈ 0.034 in margin units** (≈ 0.041 for `h >= 1.15`).

## 3. Priority 2: the equality class z = k = 4  (`k4 e4 z4 K1111 E1111 S3333`)  [measured]

* It reaches 0 at every `h`, from the permutation-hole tiling seeds within ≤ 6 restarts.  The best negative found
  among the non-zero optima is about −1e-9 to −0.008.
* Plateau: 800 extra restarts per h; **274–301 of them end at 0 (to 1e-14)**.  Most are axis-parallel.
  * 19–25 per h have a *touching* square (own gap ≤ 1e-7, so not a rattler) tilted more than 1°, and 10–17 more
    than 10°.
  * The largest touching tilt is **44.74°** (`h = 1.15`).  At the other h it is 41.6° / 36.9° / 32.7°.
  * Witnesses: `runs/witnesses/h<h>_111111114_plateau_maxtilt.json`.
  * These match `n12-gap.md` §3.1–3.2: the zero set is held shut by an axis row or column jam, while tilted squares
    sit against it at 0.

## 4. Priority 3: k <= 3, and which classes violate `z <= k − 1` at margin 0  [measured]

Every enumerated class violates `z <= k − 1`, since `z >= k` is forced.  **18 classes reach 0** at all four h.
Besides the equality class, that leaves 17 classes:

| class | walls with strip = 3 (saturated) | witness (h = 1.15) |
|---|---|---|
| `k4 e4 z4 K1111 E1111 S3333` (equality) | 4 | `runs/witnesses/h1.15_111111114.json` |
| `k3 e5 z4 K1110 E1022 S3233` | 3 | `h1.15_111010224.json` |
| `k3 e5 z4 K1110 E1121 S3332` | 3 | `h1.15_111011214.json` |
| `k2 e6 z4 K1010 E2112 S3223` | 2 | `h1.15_101021124.json` |
| `k2 e6 z4 K1010 E2121 S3232` | 2 | `h1.15_101021214.json` |
| `k2 e6 z4 K1010 E2211 S3322` | 2 | `h1.15_101022114.json` |
| `k2 e6 z4 K1010 E2220 S3331` | 3 | `h1.15_101022204.json` |
| `k2 e6 z4 K1100 E0222 S2323` | 2 | `h1.15_110002224.json` |
| `k2 e6 z4 K1100 E1212 S3313` | 3 | `h1.15_110012124.json` |
| `k2 e6 z4 K1100 E1221 S3322` | 2 | `h1.15_110012214.json` |
| `k1 e7 z4 K1000 E2212 S3213` | 2 | `h1.15_100022124.json` |
| `k1 e7 z4 K1000 E2221 S3222` | 1 | `h1.15_100022214.json` |
| `k0 e8 z4 K0000 E2222 S2222` (tiling minus the 4 corners) | 0 | `h1.15_000022224.json` |
| `k3 e6 z3 K1110 E1122 S3333` | 4 | `h1.15_111011223.json` |
| `k2 e7 z3 K1010 E2221 S3332` | 3 | `h1.15_101022213.json` |
| `k2 e7 z3 K1100 E1222 S3323` | 3 | `h1.15_110012223.json` |
| `k1 e8 z3 K1000 E2222 S3223` | 2 | `h1.15_100022223.json` |
| `k2 e8 z2 K1010 E2222 S3333` | 4 | `h1.15_101022222.json` |

Strip saturation in the zero set goes from 0 walls (`S2222`) to 4.  So the per-wall count 3 is not what makes these
configurations tight.

**An empirical dichotomy** (all 4 h, 235 classes) [measured]:
* A class reaches margin 0 **iff `z <= 4` and every wall has `E_w <= 2`**.  All 72 class-h pairs with `z <= 4` and
  `max E_w <= 2` reach 0.
* The rest are `<= −0.0339` when `z >= 5`, and `<= −0.0503` when `z <= 4` with some wall holding 3 edge squares
  (e.g. `k2 e8 z2 E1232` and `k1 e9 z2 E2322`; 800 restarts each).

## 5. Side remark: the strip lemma itself at margin 0

The free search reaches margin 0 with strip count **4** on a wall (a full axis row along a wall), e.g.
`k4 e5 z3 S4333`, the most common free optimum.  So "≤ 3 centres within h of a wall" holds only for positive
margin, i.e. side < 4.  It fails at closed side 4 on the zero set.  This is consistent with the dilation transfer
(which needs slack), but it is worth stating: at `δ = 0` the strip inequality is not available, so the count cannot
be run on touching configurations.

## Verdict

1. **"z >= 5 costs" is plausible as a joint lemma, with room ≈ 0.034.**  [measured]  For every class with `z >= 5`,
   the best margin is `<= −0.0339` at `h = 1.1213` and `<= −0.041` at `h >= 1.15`.  The local `Z`-only problem has
   +0.121.  In counting terms: if 12 squares fit with margin `> −0.03` then `z <= 4`, equivalently
   `ΣE + 2ΣK >= 8 + ΣK`.  A second lemma with similar room also holds: three edge squares on one wall (with
   `z <= 4`) cost ≥ 0.050.
2. **But the count does not close, and `z <= k − 1` has no room in 17 classes besides the equality one.**
   [measured]  "z >= k forces margin < 0 except at the equality case" is **false at margin 0**.
   * The 17 other classes, with `(k, z) ∈ {(3,4),(2,4),(1,4),(0,4),(3,3),(2,3),(1,3),(2,2)}`, are realised at
     margin exactly 0, robustly across all h.  An example is the tiling minus its four corner squares
     (`k = 0, z = 4`).
   * So after the `z >= 5` lemma, what is left of the count is the statement that each of these 18 zero-margin
     classes has no positive-margin realisation.  Each is a touching-only, rank-8-type rigidity statement of the
     same kind as the original (§3.1 of `n12-gap.md`).
   * The occupancy split partitions the zero set rather than shrinking it.  The equality case is not special: 17
     other classes sit on the same plateau.
3. If a counting proof is pursued anyway, the useful slack lemmas are:
   * (A) `z >= 5` ⇒ margin `<= −0.034`;
   * (B) some `E_w = 3` ⇒ margin `<= −0.050`.

   Those two together reduce the census to the 18 zero classes (`z <= 4`, all `E_w <= 2`).  Those still need the
   rigidity argument, and it can only work at δ = 0 exactly, without slack.
