# bandcut-scan: the tilted-band family from `T = 12` down to `T = 4`

*2026-09-20 (night).  Task `tasks/bandcut-scan/README.md`.  New code: `search/bandcut_scan.py`
(exact `Q[sqrt 2]` margin evaluator + sparse disjunctive LP with wall anchors),
`search/bandcut_scan_svg.py` (read a published packing), `search/bandcut_scan_exact.py`
(rational certificate), `search/bandcut_scan_family.py`, `search/bandcut_scan_bar.py`,
`search/bandcut_scan_pT.py`, `search/bandcut_scan_down.py`.  Runs in gitignored
`runs/bandcut_scan_*`; every number below is quoted from one of them.*

**Convention throughout** (as `S6_LOCAL.md`): every configuration produced here is a **feasible
point** of the disjunctive LP of `S6_SKELETON.md` §3.1, i.e. a **lower bound** on `delta*_T`.
`delta > 0` means the `n` closed unit squares are pairwise disjoint inside `[0,T]^2` with room
to spare, i.e. `s(T^2-T) < T`.  `delta = 0` means a closed packing of side exactly `T` and is
**not** a counterexample.

---

## 0. Verdict

1. **The instrument verifies on the known positives, in exact rational arithmetic.**  The
   published record packings for `n = 132` (`T = 12`) and `n = 110` (`T = 11`) were read off
   `kingbird.myphotos.cc`, re-optimised in the centres by our own LP at fixed angles, rounded
   to rationals, and evaluated with our own pair/wall rows in `Fraction` arithmetic:

       T = 12 :  delta = 26555450193717668780862269255541 / 39020156737952405591711252900000000
                       = +6.805572405066e-04  > 0      ==>  s(132) < 12, exactly
       T = 11 :  delta = 57695914840601803503345378222961 / 209883723780129332204734621000000000
                       = +2.748946597738e-04  > 0      ==>  s(110) < 11, exactly

   These are **verifications of other people's packings**, not new results; but they are the
   first time this repo's `delta*` machinery has seen a positive case, and they are exact.
2. **The mechanism is `T`-free; only the budget is not.**  The whole family runs on one
   inequality, which the authors state in the same words we derived it in (§1): a stack of `p`
   unit squares at a common tilt `t` spans `p cos t + sin t`, which is `< p` as soon as
   `t > t_p = 2 arctan(1/p)`.  That stack exists in any container of side `>= p`.  What
   depends on `T` is the **area budget**: `n = T^2 - T` leaves exactly `T` free cells, and the
   band has to fit its waste inside them.
3. **What runs out is the budget, and it runs out at exactly `T = 12`** (§4).  In the authors'
   own assembly scheme the square is cut into four rectangles `A, B, C, D`, with `C, D` integer
   blocks and `A, B` diagonally opposite *squeezable* rectangles; the identity is
   `waste(A) + waste(B) = T` exactly.  The smallest-waste squeezable rectangle in the
   literature is `(4,8)` with `26` squares, **waste 6**; so the scheme needs `T >= 12`.  At
   `T <= 5` the ingredient does not even fit (`32` cells `> T^2`).  Neither cut-off is an
   angle-geometry statement -- confirming `ARCH_TLEDGER.md` §5.5.
4. **Clean negative picture at `T <= 10`.**  Over three independent instruments (§3) the best
   margin found in the band family is **exactly `0`** at every `T` from `4` to `10`, never
   positive, and the binding row at each optimum is a *wall row* or the *first link of the band
   itself*, never a slack interior pair.  At `T = 4` the best band configuration is a
   margin-`0` closed packing (§5): a wall-to-wall `(4,1)` stack at `28.59°` with `8`
   axis-parallel squares, `k = 8` near-axis -- inside the hole `T <= k <= (T-1)^2` of
   `notes/proof-architecture.md` §0a, and still containing a wall-to-wall chain of four
   axis-parallel squares, so it is *not* a Cleemann-type obstruction.
5. **Our parametric reconstruction of the family is strictly weaker than the published one**
   (§6), and we say so: the minimal published ingredient `(4,8)/26` has `delta > 0.01777`
   (their Table 2) and our best parametric rebuild of it reaches `-7.04e-03`.  So the `T = 12`
   positive above rests on the published coordinates, not on our reconstruction.  Nothing in
   the scan is a claim that no better band configuration exists at `T <= 10`.

---

## 1. The mechanism, exactly

> **[pub]** Arslanov-Mustafin-Shangitbayev, *Improved packings of `n(n-1)` unit squares in a
> square*, Electron. J. Combin. **28**(4) (2021) #P4.22, DOI 10.37236/8586, §2:
> "In Figure 2 we see one of the main ideas for packing unit squares: using of stacks `(4, 1)`
> tilted by an angle `α = arcsin(8/17)`.  The main idea for squeezing a packing follows from
> it: tilting stacks `(4, 1)` by an angle `α + ε` so that the stack `(4, 1)` is located in a
> vertical strip of width `4 − δ`, where `ε` and `δ` are sufficiently small."

A **stack `(p,1)`** is `p` unit squares in a row at a common tilt `t`; its bounding box is
`(p cos t + sin t) x (cos t + p sin t)`.  The paper's defining equations are literally
`4x_1 + y_1 = 4 - δ`, `5y_2 + x_2 = 5 - δ`, `6y_1 + x_1 = 6 - δ` with `(x,y)` the unit vector
along the square's side; at `δ = 0` this is

    p cos t + sin t = p    <=>    tan(t/2) = 1/p    <=>    t = t_p = 2 arctan(1/p)

which is exactly Part 3 of `ARCH_TLEDGER.md` §4 with the **stack length** in place of the
container side, and `arcsin(8/17) = 2 arctan(1/4)` is the case `p = 4`
(`bandcut_scan.py selftest` checks `p cos t_p + sin t_p = p` in exact `Fraction` arithmetic at
`p = 2,3,4,5,6,12`).  Published tilts, all of this form:

| `p` | 3 | 4 | 5 | 6 | `2k` |
|---|---|---|---|---|---|
| orientation `(cos, sin)` | `(4/5, 3/5)` | `(15/17, 8/17)` | `(12/13, 5/13)` | `(35/37, 12/37)` | `((4k²−1)/(4k²+1), 4k/(4k²+1))` |
| `t_p` | `36.8699°` | `28.0725°` | `22.6199°` | `18.9246°` | — |

(the `p = 4, 5, 6, 2k` rows are **[pub]**, from §2, Fig. 3, Fig. 4 and Lemma 1; the closed form
`t_p = 2 arctan(1/p)` is our reading, consistent with all four).  **There are no `45°` squares
anywhere in the paper** -- the three `45°` squares of Cleemann's `T = 17` packing
(`ARCH_TLEDGER.md` §5.2) are a different construction.

So a bar of `p` squares tilted past `t_p` and laid across `p` columns leaves **slack**
`p - (p cos t + sin t) > 0` in those columns; that is the entire source of `s(T^2-T) < T`.
Its price is the bar's height `cos t + p sin t`, which has to be paid out of the `T` free cells.

---

## 2. The instrument, and the two exact certificates

`search/bandcut_scan.py` has three pieces.

* **`Q2`** -- exact arithmetic in `Q[sqrt 2]` with an exact sign test, enough for every angle
  the family uses (Pythagorean tilts are rational; `45°` needs the `sqrt 2`).
* **`exact_delta`** -- `min(min pair gap, min wall slack)` over the rows of `S6_SKELETON.md`
  §3.1, evaluated exactly.  Cross-checked against `s6skel.value`: **max `|diff| = 3.33e-16`**
  over 12 random mixed-angle configurations (`bandcut_scan.py selftest`).
* **`BandLP`** -- the disjunctive LP in the centres at fixed angles, sparse, restricted to
  pairs within a cut-off (dropping rows only *raises* the LP value, and the point it returns is
  re-evaluated over **all** pairs, so the answer stays a feasible point), plus optional
  **wall-anchor equalities** `x_i - p_i = delta`, `T - x_j - p_j = delta` which force a bar to
  be wall-to-wall.  Without them the LP simply slides the band into a corner and reports the
  axis-parallel plateau.

`search/bandcut_scan_svg.py` reads the published files at
<https://kingbird.myphotos.cc/packing/square-N.svg> (blocks of axis-parallel squares as filled
rectilinear regions on the integer grid, everything else as `use ... translate(cx cy) rotate(a)`).
`search/bandcut_scan_exact.py` then (i) replaces each angle by the nearest angle with rational
half-tangent (denominator `1e7`), (ii) re-solves the LP in the centres, (iii) rounds the centres
to denominator `1e9`, (iv) evaluates exactly.

| `n` | `T` | published side `s` | squares parsed | `delta` at `T`, published angles | after our LP | **exact** |
|---|---|---|---|---|---|---|
| 110 | 11 | `10.99679327401957174` | 110 ✓ | `+1.458027763e-04` | `+2.748954098e-04` | **`+2.748946598e-04`** |
| 132 | 12 | `11.99137344423646923` | 132 ✓ | `+3.596984033e-04` | `+6.805580134e-04` | **`+6.805572405e-04`** |
| 156 | 13 | `12.98208376048414436` | 156 ✓ | `+6.900371253e-04` | — | — |
| 182 | 14 | `13.97419105332569700` | 180 (2 short) | `+9.234504729e-04` | — | — |
| 210 | 15 | `14.97413341886404581` | 208 (2 short) | `+8.637087841e-04` | — | — |
| 240 | 16 | `15.97556282833087771` | 236 (4 short) | `+7.648297569e-04` | — | — |
| 272 | 17 | `16.96971602419903036` | 268 (4 short) | `+8.922947136e-04` | — | — |

Rows where the parse is short are **indicative only** (a missing square leaves extra room);
`T = 11, 12, 13` parse exactly and are the ones we rely on.  Every file's bounding box comes out
as `[0,s]^2` to `1e-12` and `delta` at side `s` is `0` to `1e-15`, which is the check that the
parse is right.

Structure of the two verified packings (angle census mod `90°`, `bandcut_scan_svg.py --map`):

| | axis-parallel | tilt cluster 1 | tilt cluster 2 | shape |
|---|---|---|---|---|
| `T = 11`, `n = 110` | 64 | `23.9 … 28.3°` (21) | `57.8 … 65.1°` (25) | one bent band, bottom-left to top-right |
| `T = 12`, `n = 132` | 79 | `22.2 … 28.3°` (30) | `59.3 … 66.4°` (21) | same |

Both are annealed refinements (Ellsworth Jan 2026, Stead Jun 2026) of the analytic originals,
so the tilts are spread rather than sitting on `t_p` exactly; the `T = 12` original of
Arslanov et al. is the two-angle version, `82` axis-parallel + `50` tilted at `t_4` and `t_5`
only.

**For `T <= 10` the records page lists no packing at all**, which by its own convention
("For the `n <= 324` not pictured, the trivial packing (with no tilted squares) is the best
known packing") means `s(T^2-T) = T` is still the record for `T = 4 … 10`.  So the frontier is
`T = 11` and the ledger of `ARCH_TLEDGER.md` §5.1 is current.

---

## 3. The scan: `delta(T)` for `T = 12 … 4`

Three instruments, each producing feasible points; the table reports the best of the three.

**(a) Seam contraction** (`bandcut_scan_down.py`).  Start from the real `T = 12` packing; one
step `T -> T-1` deletes every square whose centre lies in a unit vertical seam and a unit
horizontal seam, closes the container up to `[0,T-1]^2`, restores the count to `(T-1)^2-(T-1)`
and re-optimises.  20 seam choices per step, best kept.

**(b) Anchored bar + fill** (`bandcut_scan_bar.py`).  A bar of `p` squares at tilt `t`,
anchored wall-to-wall by LP equalities, optionally crossed with a second bar at `90 - t`
anchored in the other direction; everything else axis-parallel, LP-optimised, with a relocation
loop that may not touch the band.  Swept over `t in [20°, 70°]`, `p in {2,…,8}`, offsets.

**(c) Full-width `(T,1)` stack** (`bandcut_scan_pT.py`) -- the sharp end: `p` squares at a
common tilt just above `t_p`, wall-to-wall, `p` running over `T, T-1, T-2, T-3`, tilt swept in
`0.5°` then `0.02°` steps.

| `T` | `n` | (a) contraction | (b) anchored bar | (c) full-width stack | **best** | status |
|---|---|---|---|---|---|---|
| 12 | 132 | `+3.596984e-04` (the seed) | — | `-2.220446e-16` (`p=12`, `t=15.55°`) | **`+6.805572405e-04`** | **positive, EXACT** (published packing, our rows) |
| 11 | 110 | `-4.535302e-03` | — | `-4.440892e-16` (`p=11`, `t=16.41°`) | **`+2.748946598e-04`** | **positive, EXACT** (published packing, our rows) |
| 10 | 90 | `-4.931367e-03` | `+0.000000000e+00` | `-1.110223e-16` (`p=10`, `t = 14.44°`) | **`0`** | zero, attained; nothing positive found |
| 9 | 72 | `-5.377948e-03` | `+0.000000000e+00` | `-1.110223e-16` (`p=9`, `t=18.70°`) | **`0`** | zero, attained |
| 8 | 56 | `-9.526935e-03` | `+0.000000000e+00` | `+0.000000000e+00` (`p=8`, `t = 14.27°`) | **`0`** | zero, attained |
| 7 | 42 | `-9.160453e-03` | `+0.000000000e+00` | `+0.000000000e+00` (`p=7`, `t = 16.2802°`) | **`0`** | zero, attained |
| 6 | 30 | `-3.606212e-03` | `+0.000000000e+00` | `+0.000000000e+00` (`p=6`, `t = 18.9446°`) | **`0`** | zero, attained |
| 5 | 20 | `-1.110223e-16` | `-5.356104e-03` | `-5.867512e-05` (`p=5`, `t = t_5` exactly) | **`0`** | zero, attained (the band dissolves) |
| 4 | 12 | `+0.000000000e+00` | `-5.697534e-03` | `+0.000000000e+00` (`p=4`, `t = 28.5925°`) | **`0`** | zero, attained — see §5 |

(Dashes are combinations still running or not attempted at that `T`; `runs/bandcut_scan_pT_T*.txt`
and `runs/bandcut_scan_bar_T*.txt` have the per-angle sweeps.)

**One full-width stack is never enough, at any `T`.**  Column (c) is `0` to `1e-16` at
`T = 11` and `T = 12` as well -- where a positive certainly exists.  A single `(T,1)` stack
threaded across the container buys slack in the rows it crosses and pays for it in the columns,
and the two cancel exactly; the published positives come from the *assembled* scheme of §4 (two
squeezable rectangles), not from one stack.  So the scan's `0`s at `T <= 10` are consistent with
the family being alive at `T = 11, 12` and say only that our generator did not find the
assembly at smaller `T`.

**Reading.**  `delta(T) = 0` for every `T` from `4` to `10` -- attained, never exceeded.  That is
the expected and useful outcome: the band family does not produce a counterexample below
`T = 11`, and at `T <= 10` it lands on `Z_T` rather than `P_T`.  The two positive rows are the
published packings, and they are the *only* positives anywhere in the scan.

**What binds.**  At every `T <= 10` optimum the binding row is one of

* a **wall row** (`('wall', i, 'L')`, reported by `bandcut_scan_exact.py` at `T = 4, 5`), i.e.
  the configuration is pressed flat against the container, or
* the **first link of the band itself** (`pair (0,1)`), i.e. the bar is edge-to-edge tight and
  cannot be spread any further without leaving the walls.

Never an interior pair with slack to give.  The `T = 5` case is the most informative single
number in the scan: a full-width `(5,1)` stack at **exactly** `t_5 = 2 arctan(1/5)` with 15
axis-parallel squares misses by `-5.867512e-05` (exact: `-76277651/1300000000000`), and the
deficit grows monotonically as the tilt is pushed above `t_5` (`-6.06e-03` at `+1.4°`,
`-1.02e-02` at `+15°`, `runs/bandcut_scan_bar_T5.txt`).  The band is *already* at its best angle
and still short: what is missing is not tilt, it is room.

---

## 4. What runs out: the waste budget, and why it is `12`

This is the part of the answer that is a theorem-shaped statement about the family rather than a
measurement.

**[pub]** (§3 of the paper, Figure 1 "Scheme of squeezable packing"): for even `n >= 14` the
square is assembled from rectangles `A = (12,6)`, `B = (n-10, n-6)`, `C = (10, n-6)`,
`D = (n-12, 6)`; for odd `n >= 13` from `A = (10,5)`, `B = (n-9, n-5)`, `C = (9, n-5)`,
`D = (n-10, 5)`.  `C` and `D` are integer blocks (trivially packed, unsqueezable); `A` and `B`
are diagonally opposite and carry all the waste, and **[pub]** `δ(S, ·) >= min(δ(A, ·), δ(B, ·))`.

**The budget identity** (elementary, ours).  If `[0,T]^2 = A ∪ B ∪ C ∪ D` with `C, D` integer
blocks holding exactly their area, then since `n = T^2 - T`,

        waste(A) + waste(B)  =  T        exactly.

**The ingredients, and their waste** (all **[pub]**, Table 2 and Lemma 1 of the paper):

| rectangle | squares | area | **waste** | `δ` |
|---|---|---|---|---|
| `(4,8)` | 26 | 32 | **6** | `> 0.01777021751` |
| `(5,9)` | 38 | 45 | 7 | `> 0.020403` |
| `(5,10)` | 43 | 50 | 7 | `> 0.0009652493` |
| `(6,11)` | 58 | 66 | 8 | `> 0.01681735886` |
| `(6,12)` | 64 | 72 | 8 | `> 0.004908231774819` |
| `(2k, 2k+4)`, `k >= 3` (Lemma 1) | `4k²+6k−2` | `4k²+8k` | `2k+2` | — |

> **[pub] Lemma 1.** "For any `k >= 3` there exists a squeezable packing of `4k² + 6k − 2` unit
> squares in a rectangle `(2k, 2k + 4)` (the waste is equal to `2k + 2`)."

**Minimum waste over every squeezable rectangle the literature knows is `6`.**  Hence

        T  =  waste(A) + waste(B)  >=  6 + 6  =  12 ,

i.e. **the scheme is empty for every `T <= 11`**, which is exactly where the theorem stops.  The
paper's own boundary is the same one read off the recipes: the odd recipe needs `n - 9 >= 4`
(`n >= 13`), the even one `n - 10 >= 4` (`n >= 14`), Lemma 1 needs `k >= 3`, and `n = 12` is a
one-off (Figure 6, from `(8,4)` and `(5,10)`).  At `n = 11` the odd recipe gives `B = (2,6)`,
which has no squeezable packing, so `min(δ(A), δ(B)) = 0` and the method gives nothing.
(Cantrell's `T = 11` packing is *not* from this scheme; it is numerically optimised.)

**Two different things run out, at two different places.**

* `T <= 5`: the *geometry*.  The smallest ingredient `(4,8)` needs a `4 x 8 = 32`-cell region and
  `T^2 <= 25`.  At `T = 4` there is no room for a squeezable sub-rectangle at all; the only
  wall-to-wall stack that fits is the `(4,1)` at `t_4`, which spans the container exactly and
  leaves `4 x 1.2353` above it (§5).
* `6 <= T <= 11`: the *budget*.  The ingredient fits but `T < 12 = 6 + 6`.

Neither is about the tilt.  `t_p = 2 arctan(1/p)` is available at every `T >= p`; the `(4,1)`
stack "exists in any container of side `>= 4`" (`ARCH_TLEDGER.md` §5.5), and our scan confirms
it -- at `T = 4` the stack is there, at the right angle, and the margin is still `0`.  **The
`T`-dependence of the whole problem sits in the area/counting side, i.e. in `A3`/`A5`, and not
in `A2c` or in any angle-geometry lemma.**

The one thing the budget argument does *not* do is close the question: "6 is the smallest waste
of a squeezable rectangle" is a statement about the literature, not a theorem.  A squeezable
rectangle of waste `<= 5`, or of waste `6` in a smaller box, would move `T*` down immediately.
That is the sharpest open sub-problem this scan produces, and it is one dimension smaller than
`s(12) = 4`.

---

## 5. `T = 4`: the best band-type configuration is a margin-`0` packing

`bandcut_scan_pT.py --T 4 --fine` (`runs/bandcut_scan_pT_T4.txt`).  Best over the family:

    delta = 0   (float, and -3.59e-13 after rounding the centres to denominator 1e12),
    p = 4,  t = 28.5925 deg  (t_4 = 28.0725 deg, excess +0.52 deg),  k = 8 near-axis.

The configuration (`runs/bandcut_scan_pT_T4.json`):

| squares | positions |
|---|---|
| 4 axis-parallel, bottom row | `(0.5, 0.5) (1.5, 0.5) (2.5, 0.5) (3.5, 0.5)` |
| 2 axis-parallel, lower right | `(2.5, 1.5) (3.5, 1.5)` |
| 2 axis-parallel, upper left | `(0.5, 3.5) (1.5, 3.5)` |
| 4 tilted at `28.5925°`, wall-to-wall | `(0.67831, 1.67831) (1.51815, 2.30682) (2.45456, 2.67831) (3.32169, 3.32169)` |

Variants at `T = 4`, all measured (`runs/bandcut_scan_pT_T4.txt`, `runs/bandcut_scan_bar_T4.txt`):

| variant | best `delta` |
|---|---|
| full-width `(4,1)` stack, `t = 28.5925°` | **`0`** |
| full-width `(3,1)` stack, `t = 36.39°` | `-1.202949e-02` |
| `(2,1)` stack | no wall-to-wall stack fits (`t_2 = 53.13°`, the rise leaves the container) |
| `(5,1)`, `(6,1)` stacks | no wall-to-wall stack fits (`5 cos t + sin t <= 4` forces `t >= 49.6°`, and then the stack's height `cos t + 5 sin t = 4.46 > 4`) |
| anchored bar, free tilt, `p in {2,…,6}` | `-5.697534e-03` (`p = 3`, `t = 22°`) |
| **two bands** (cross: a bar at `t` anchored in `x`, a bar at `90-t` anchored in `y`) | worse than every single bar at every `t` tried |
| bent band (two drift segments, the elbow) | `< 0` everywhere |
| bands of `3`- and `5`-chains at `2 arctan(1/j)`, profiles `2,3,4,4,3,2`, `3,4,4,4,3`, `4,4,4,4` | `< 0` everywhere |

So the only `T = 4` band configuration that reaches `0` is the single `(4,1)` stack, and nothing
reaches more.

**Three things worth recording about this point.**

1. It is a **point of `Z_4`, not of `P_4`** -- a legitimate closed packing of `[0,4]^2` with
   margin exactly `0`, consistent with `s(12) = 4` and with the measured `max delta*_4 = 0` of
   `RANK8.md` / `CENSUS.md`.  Nothing here challenges the conjecture.
2. `k = 8` near-axis squares puts it **inside the hole** `T <= k <= (T-1)^2`, i.e. `4 <= k <= 9`
   (`notes/proof-architecture.md` §0a item 2), and it is a *new* kind of inhabitant: the
   previously known generic stratum there is "a chain of four axis-parallel squares beside eight
   squares tilted up to `12.8°`"; this one has a wall-to-wall **tilted** stack of four at
   `28.6°` and only eight near-axis squares.
3. It nevertheless **still contains a wall-to-wall chain of four axis-parallel squares** (the
   bottom row), so `A4`'s conclusion holds on it and the chain bound `B_4 <= 0` of
   `ARCH_TLEDGER.md` §2.1 applies.  Our search never produced a `T <= 10` configuration in which
   the band cuts *every* row *and* every column -- the Cleemann/Arslanov obstruction -- with
   margin `>= 0`.  At `T = 4` there is simply no room: a `(4,1)` stack has height
   `cos t + 4 sin t = 47/17 = 2.7647` at `t_4`, so a second, transverse stack cutting the columns
   would need another `2.76` of the remaining `1.235`.

At the exact window angle, with the band as a rigid straight stack (centres
`(23/34 + 15k/17, y_0 + 8k/17)`, `cos t_4 = 15/17`, `sin t_4 = 8/17`) and the eight axis-parallel
squares on the integer grid, the exact margin is `-361/5780 = -0.06246`; the `0` above is
reached only after the LP relaxes the band out of a straight line.  Both numbers are in
`runs/`.

---

## 6. Honest limits of this scan

* **Our parametric rebuild of the family is weaker than the published construction.**  The
  decisive test is the minimal ingredient: the published `(4,8)` with `26` squares has
  `δ > 0.01777` **[pub]**; sweeping our band generator over tilts `28.1–34°`, bar profiles
  (`2,3,4,4,3,2`, `3,4,4,4,3`, `4,4,4,4`, …), one- and two-segment drifts (the elbow) and
  anchors gives at best **`-7.038e-03`** (`runs/bandcut_scan_rect48.txt`).  We did not
  reproduce it.  Consequently the `T = 12` positive in §2 comes from the published coordinates,
  and the `T <= 10` zeros in §3 are the best our family generator found, not a bound.
* **Seam contraction is a weak instrument.**  Contracting the real `T = 12` packing gives
  `-4.5e-03` at `T = 11` where the truth is `>= +2.75e-04`, i.e. it loses about `5e-03` per
  step.  Its value is only as a third opinion.
* **The `chain x/y` column printed by `bandcut_scan_bar.py` over-counts**: it is a longest path
  in the tight-pair digraph, which in an axis-parallel block runs through diagonal contacts as
  well.  It is not the wall-to-wall chain length of `A4` and should not be quoted as such; the
  binding-row report (`wall` vs `pair`) is the reliable diagnostic.
* Nothing here is a *proof* that `delta* <= 0` at any `T <= 10`, in the band family or outside
  it.  Every entry is a feasible point.

---

## 7. What this does to the architecture

| claim | verdict from this scan |
|---|---|
| `notes/proof-architecture.md` §0a item 8, "`arctan(8/15) = 2 arctan(1/4)`: the band is built from chains of four tilted squares at the `T = 4` window angle" | **confirmed and generalised**: the authors use `(p,1)` stacks at `2 arctan(1/p)` for `p = 3,4,5,6,2k`, and state the equation `p x + y = p - δ` themselves. |
| `ARCH_TLEDGER.md` §5.5, "the `T`-dependence is split between `A4` (the band cuts every row and column) and `A3`/`A5` (whether the container is big enough for the band to pay for itself)" | **the second half is the operative one.**  The budget identity `waste(A) + waste(B) = T` with minimum ingredient waste `6` reproduces the frontier `T >= 12` exactly, with no angle input at all. |
| `ARCH_TLEDGER.md` §6, "at `T = 4` a `(4,1)` stack already spans the container exactly, so `n = 12` cannot afford the elbow" | **measured**: the `(4,1)` stack at `T = 4` reaches margin `0` and no more, and there is no room (`1.235` of height left) for a transverse stack to cut the columns. |
| `notes/proof-architecture.md` §2, "the hole is `T <= k <= (T-1)^2`" | §5 adds a new inhabitant of the hole at `T = 4`, `k = 8`, with a wall-to-wall tilted stack -- but it still contains an axis-parallel chain of four, so `A4` survives on it. |

**The sharpest question this leaves**, and it is strictly smaller than `s(12) = 4`: *what is the
minimum waste of a squeezable rectangle?*  If it is `6` (attained by `(4,8)/26`) then the
Arslanov scheme cannot reach below `T = 12` and the `T = 11` counterexample must come, as
Cantrell's does, from outside the scheme.  A waste-`5` squeezable rectangle would immediately
give `T* <= 10` by the scheme (`5 + 5`), a waste-`5` paired with the known waste-`6` gives
`T* <= 11`, and a waste-`4` one would give `T* <= 8`.  This is a rectangle
question in at most `6 x 12`, it is decidable by exactly the machinery in this file (the LP now
takes a rectangular container), and it is the natural next run.

---

## 8. Reproduce

    python3 search/bandcut_scan.py selftest                                    #  2 s

    # read the published packings (the .svg files are downloaded from kingbird.myphotos.cc)
    python3 search/bandcut_scan_svg.py square-132.svg --T 12 --map --out runs/bandcut_scan_p132.json
    python3 search/bandcut_scan_svg.py square-110.svg --T 11 --map --out runs/bandcut_scan_p110.json

    # the two exact certificates                                               # 40 s each
    python3 search/bandcut_scan_exact.py --in runs/bandcut_scan_p132.json --T 12
    python3 search/bandcut_scan_exact.py --in runs/bandcut_scan_p110.json --T 11

    # the scan
    python3 search/bandcut_scan_pT.py --T 4 --fine --out runs/bandcut_scan_pT_T4.json     #  6 min
    python3 search/bandcut_scan_pT.py --T 5 --fine --out runs/bandcut_scan_pT_T5.json
    setsid nohup python3 -u search/bandcut_scan_bar.py --T 4 --deg 20:70:2 --p 2,3,4,5,6 \
        --cross --offstep 1.0 --topk 12 > runs/bandcut_scan_bar_T4.txt 2>&1 &            # 10 min
    setsid nohup python3 -u search/bandcut_scan_down.py --start runs/bandcut_scan_p132.json \
        --T 12 --steps 8 --tries 20 --seed 3 > runs/bandcut_scan_down12.txt 2>&1 &       # 25 min
