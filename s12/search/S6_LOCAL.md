# s6-local: the deficit off the zero set, re-measured and made exact on single leaves

*2026-09-20 (afternoon).  Follow-up to `S6_SKELETON.md`; corrects three of its numbers.  Code:
`search/s6local.py` (instrument with structured starts), `search/s6exact.py` (symbolic dual on one
leaf), `search/s6leaves.py` (tight-leaf count).  Runs in gitignored `runs/s6local_*.txt`,
`runs/s6exact_closedforms.txt`, `runs/s6leaves_*.txt`; every number is quoted here.*

Notation as `S6_SKELETON.md` §3.1: `delta*(theta)` is the value of the disjunctive LP in the centres
at a fixed angle vector; here with **hard containment** (the wall rows carry no `delta`; only pair
rows do), which is what `fixed_angle_value` measured too.  Everything from `s6local` is a **lower
bound** on `delta*` (a feasible point); everything from `s6exact` is an **exact upper bound for one
assignment** (one leaf), not for `delta*`.

## 0. What changed

| claim in `S6_SKELETON.md` | now |
|---|---|
| `n = 12`: `delta*(t,..,t) = -(2/3) t^2` | **wrong** (under-optimised): `>= -(1/3) t^2 + (2/9) t^3 + ...`, stable over `0.05°–6.4°`, optimum is a pinwheel (tiling minus a 4-cycle) |
| `n = 6`: cubic direction `~ -0.17 t^3` (two of three points below LP tolerance) | **`-(1/4) t^3 + (1/2) t^4 + ...`**, clean over `0.4°–25.6°` with tolerance `1e-10`, and exact on its leaf |
| `n = 12`: "cubic constant is the number to measure next" | **there is no cubic direction in the uniform family at `n = 12`**: `j >= 4` axis-parallel gives exactly `0`, `j <= 3` gives the *same* `-(1/3) t^2` to 10 digits |
| (not measured) leaving `Z_12` from a generic point | **first order**: `delta* ~ -0.31 s` in the tilt `s` of one of the four axis-parallel squares (expected `-s/3` for a straight chain), saturating at the pinwheel value |

The verdict of `S6_SKELETON.md` (an angle-box tree cannot close) never depended on the cubic: a box
meeting `Z` has `sup delta* = 0` and any positive relaxation error keeps it alive.

## 1. The instrument

`fixed_angle_value` starts from uniform random centres.  At `n = 12` that finds the good family in
about 1 run in 6 (3 seeds x 4 tilts: `-0.3326` at `0.2°` and `-0.3271` at `1.6°`, otherwise `-0.66`);
the `0.8°` row of `runs/s6_scaling.txt` (`-0.330`) was that family showing through.  A multistart
maximum can only err low, so an outlier *above* the trend convicts the trend.  `s6local` starts
from every `n`-subset of the `T^2` tiling cells (labels dealt over the angle classes), jittered by
`0.08` — the jitter breaks the tie at diagonal neighbours, which is the choice between the row and
the pinwheel assignment — then continues the best configurations up and down in `t`; HiGHS
tolerances `1e-10`.  It reproduces `n = 6` uniform tilt exactly (`-0.74217` at `0.4°`).

## 2. Deficit by number `j` of axis-parallel squares (the rest tilted by a common `t`)

*Uniform families only, all tilts of one sign.  §5 shows the general picture is not a function of `j`.*

`n = 6`, `T = 3` (`delta*/t^2`, at `0.8° / 1.6° / 3.2° / 6.4°`):

| `j` | value | order |
|---|---|---|
| 0 | `-0.734, -0.719, -0.689, —` | `-(3/4) t^2` |
| 1 | `-0.4896, -0.4793, -0.4592, -0.4205` | `-(1/2) t^2` |
| 2 | `-0.4830, -0.4667, -0.4365, -0.3837` | `-(1/2) t^2` |
| 3 | `delta*/t^3 = -0.2465, -0.2431, -0.2363, -0.2232, -0.1986, -0.1554, -0.0894` at `0.4° … 25.6°` | **`-(1/4) t^3`** |
| 4, 5, 6 | `0` to `1e-16` | in `Z` |

The `j = 3` optimum is a straight axis-parallel row of three on one wall, two tilted squares in the
far corners and the third tilted square between them (`runs/s6local_n6_core3.txt`): three tilted
squares in a `3 x 2` box, where the first-order terms cancel exactly.

`n = 12`, `T = 4` (`delta*/t^2` at `0.4° / 0.8° / 1.6° / 3.2° / 6.4°`):

| `j` | value |
|---|---|
| 0, 1, 2, 3 | `-0.33178, -0.33023, -0.32711, -0.32087, -0.30832` — **identical to 10 digits for all four `j`** |
| 4 (`0^4 t^8`) and 8 (`t^4 0^8`) | `0` to `1e-15`, up to `12.8°` |

So at `n = 12` up to three squares of the pinwheel are idle (its certificate, §3, touches a chain of
four plus a few transverse squares), and `Z_12` near the origin contains the same-sign part of "`>= 4` axis-parallel" (mixed signs untested at
`n = 12`; at `n = 6` they matter, §5), not a thin subset of it: eight squares tilted by up to `12.8°` fit beside a (not necessarily
straight) wall-to-wall chain of four with margin exactly `0`.  That is the fat plateau of `RANK8.md`
and `CENSUS.md` read in angle coordinates.

Leaving `Z_12` at a generic point, `theta = (s, 0, 0, 0, 5°^8)` (`runs/s6local_offZ_n12.txt`,
continuation from `s = 0`; hit rates are low, `3–8` of `36,100` starts, so treat as a lower bound):

| `s` | `0.1°` | `0.2°` | `0.4°` | `0.8°` | `1.6°` | `3.2°` |
|---|---|---|---|---|---|---|
| `delta*` | `-5.43e-4` | `-1.07e-3` | `-1.48e-3` | `-1.65e-3` | `-1.93e-3` | `-2.30e-3` |
| `/s` | `-0.311` | `-0.306` | `-0.213` | `-0.118` | | |

First order (`~ -s/3`, a chain of four with three `delta`-carrying links lengthened by `s`) until the
pinwheel is cheaper (`-(1/3)(5°)^2 = -2.5e-3`).  **Picture at `n = 12`:
`delta* ~ max( -(1/3) s_4 , -(1/3) |tilt|^2 )`** with `s_4` the fourth-smallest tilt — a row
mechanism that dies at first order and a pinwheel mechanism that dies at second order, and no third
one seen.  The hard regime is all twelve tilts small, where the pinwheel term is what has to be
proved.

## 3. Exact closed forms on single leaves (`s6exact.py`)

For a fixed assignment the margin is an LP with trigonometric data; any dual-feasible `w(t) >= 0`
gives `delta <= -sum w_r(t) b_r(t)` exactly.  Taking the support of the HiGHS dual at `t = 1.6°` and
solving the dual symbolically (sympy; each agrees with the LP at the sample to `1e-15`):

**`n = 6`, pinwheel, uniform tilt** (9 rows: 4 walls, 5 pairs).  With `u = cos t + sin t`:

    delta(t)  <=  -3 (u - 1)^2 / (u^2 + 3)   =   -(3/4) t^2 + (9/8) t^3 - (1/2) t^4 + ...

manifestly `<= 0`, zero only at `t = 0 mod 90°`.  Multipliers: one chain of three across `y`
(`lo-y, pair, pair, hi-y`) at `1/2 - t/4 + ...`, a transverse chain at `t/2`, two links at `t^2/2`.

**`n = 6`, `(0,0,0,t,t,t)`** (9 rows):

    delta(t)  <=  ( 6 sin 2t - 9u + sqrt2 cos(3t + pi/4) + 8 ) / ( 4 (u + 1) )   =   -(1/4) t^3 + (1/2) t^4 + ...

**`n = 12`, pinwheel, uniform tilt** (15 rows: 6 walls, 9 pairs): a trig-rational function
(`runs/s6exact_closedforms.txt`) with series `-(1/3) t^2 + (2/9) t^3 + (1/54) t^4 + ...`.  The
dominant part of the dual is **one wall-to-wall staircase chain of four** in the tilted frame
(`lo-x, pair, pair, pair, hi-x`, weight `1/3 - 2t/9` each), closed by a transverse chain at `t/3` and
links at `t^2/3` — `RANK8.md` §3's four-square chain certificate, continued off the tiling.

These are one leaf and one line each.  What they show is the *shape* of a local theorem: per leaf, an
explicit trig-rational bound with multipliers continuous in `theta`; the chain length enters as the
number of `delta`-carrying links (`1/2`-weights at `T = 3`, `1/3` at `T = 4`).

## 4. How many leaves are tight? (`s6leaves.py`)

The decision tree at the exact axis-parallel point `theta = 0`, `n = 6` (`kills = N`, `delta` free
below): **3,045 LP nodes, 2,284 leaves, all 2,284 tight** (`|value| <= 1e-9`), none strictly negative.
Two dual signatures, 1,142 each: *two walls + two pair rows on three squares*, in `x` or in `y` — the
wall-to-wall chain of three.  No leaf needs the pinwheel certificate at `theta = 0`; the pinwheel
only appears once the chains are tilted.

So the hoped-for split ("most leaves strictly negative, survive by continuity; a handful tight, done
by hand") is **false at the origin**: every leaf is tight, all for the same reason.  A local theorem
at `theta = 0` cannot be leaf-by-leaf continuity; it has to be one statement about chains — for every
assignment, some wall-to-wall chain of `T` squares, whose length
`P_a + m_ab + ... + P_z >= T` with equality iff all its squares are axis-parallel *and* its links are
straight, plus the second-order bookkeeping for staircase links (the `t/2`, `t^2/2` multipliers
above).  At generic points of `Z` (`(0,0,0,0,20°,33°)` and `(0,0,0,0,2°,3°)`) the tree did not
finish in 150,000 nodes; open.

## 5. The whole small-angle cube at `n = 6` (`s6cube.py`, evening) — **§2 above is the uniform family only**

116 directions `v` (max-norm 1, **signed** tilts: 60 dense Gaussian, 48 with 1–4 zeros, sign patterns, the
uniform families), `delta*(eps v)` at `eps = 0.5°, 1°, 2°, 4°`; 6,000 jittered tiling starts with random label
dealing at `2°`, continuation to the other scales (`runs/s6cube_n6.jsonl`).  Continuation alone is **not**
reliable at the small scales (four `sparse3` rows read `-0.4 eps^2` at `0.5°` where 12,000 fresh starts give
exactly `0`); the numbers below are from `2°` or were re-checked with fresh starts.

1. **Sign coherence decides the order.**  `delta*(eps v) ~ -a(v) eps - b(v) eps^2`.  With `k` = number of
   squares in the minority sign (|v_i| > 0.05): `k >= 3` always has `a > 0` (15 directions, `a` from `0.026`
   to `0.250`); `k = 2` usually (`a` up to `0.250`, median `0.008`); `k <= 1` has `a = 0` to fit accuracy and
   is second order, `b` from `~0` to `1.18` (uniform `+`: `3/4`; one square reversed: `1.24`).  The sign
   patterns `(-,-,+,+,+,+)` and `(-,-,-,+,+,+)` give `a = 0.252, 0.256`: **`-eps/4`, first order**, the same
   value at all four scales.  Reason: a link between squares of tilt `+t` and `-t` costs
   `m - 1 = |D|/2 = t`, and a staircase can only repay a chain whose squares turn the same way.
2. **`Z` is not a zero count.**  Four axis-parallel squares do not suffice: `(0,0,0,0,+1,-0.65)` has
   `delta* = -0.077 eps^2` and `(0,0,0,0,+1,-0.96)` `-0.229 eps^2` (fresh starts agree at `0.5°` and `2°`);
   all ten same-sign `sparse4` directions are `0`.  Three can suffice: `(0,0,0,1,b,c)` is exactly `0` at
   `(1, 0.88, 0.39)`, `(1, 0.64, 0.53)`, `(1, 0.61, 0.19)`, `(1, 1, 0.5)`, `(1, 0.5, 0.5)`, while the cubic
   deficit holds on a cone around the diagonal (`delta*/t^3` at `6.4°`: `-0.199` at `(1,1,1)`, `-0.199, -0.201,
   -0.201, -0.199` at `c = 0.99, 0.95, 0.9, 0.8`, `-0.219` at `(1, 0.95, 0.9)`), and mixed signs there are
   second order (`(1,1,-1)`: `-0.42 t^2`).  So `Z` near the origin is a union of pieces cut out by **sign
   coherence and by inequalities between the tilts**, inside "`>= 3` axis-parallel".
3. **The certificate is always a chain.**  In all 116 directions the dual of the optimal leaf has exactly
   one group of rows at weight `>= 0.2`: two walls and two pair rows along one axis — a wall-to-wall chain of
   three — (111 directions, at the freshly searched scale `2°`), or two such chains sharing the weight (5, among
   them both first-order sign patterns); total support 4–13 rows.  The rest of the support carries weight `O(eps)`.

**What this does to the plan.**  A "chain lemma" cannot be stated in the angles alone: the single-chain
inequality `2 delta <= 3 - P_a - P_c - m_ab - m_bc + tan(th_o1) Dy_ab + tan(th_o2) Dy_bc` still contains the
transverse offsets `Dy`, and eliminating them *is* the LP dual, leaf by leaf.  What the data supports is a
two-level statement: (i) a **first-order theorem** — `delta_1(v) = lim delta*(eps v)/eps <= 0`, a piecewise-linear
function computable from the axis-parallel plateau alone (linearise `P`, `m` and the normals; the unknowns are
a plateau configuration and first-order displacements), zero exactly on a "coherent" cone; (ii) a
**second-order theorem on that cone**, where the pinwheel closed form of §3 is the model case.  (i) is finite
and looks provable; (ii) is where `RANK8.md` §3.3 stopped.

## 6. Reproduce

    python3 search/s6local.py tilt --n 12 --T 4 --pattern uniform --tdeg 0.05,0.1,0.2,0.4,0.8,1.6,3.2 --show 0.4     # 35 s
    python3 search/s6local.py tilt --n 6 --T 3 --pattern 0,0,0,1,1,1 --tdeg 0.4,0.8,1.6,3.2,6.4,12.8,25.6 --jitters 6
    python3 search/s6local.py tilt --n 12 --T 4 --pattern 0,0,0,1,1,1,1,1,1,1,1,1 --limit 6000 --jitters 3       # j = 3; likewise j = 1, 2, 4
    python3 search/s6local.py tilt --n 12 --T 4 --pattern 1,0,0,0,0,0,0,0,0,0,0,0 --base 0,0,0,0,5,5,5,5,5,5,5,5 \
        --tdeg 0,0.1,0.2,0.4,0.8,1.6,3.2 --limit 12000 --jitters 3 --keep 100                                      # 5 min, detached
    python3 search/s6exact.py --n 6 --T 3 --pattern uniform ;  ... --pattern 0,0,0,1,1,1 ;  --n 12 --T 4 --pattern uniform
    python3 search/s6leaves.py --n 6 --T 3 --theta 0,0,0,0,0,0
    setsid nohup python3 -u search/s6cube.py run --n 6 --T 3 --ndense 60 --nsparse 12 --out runs/s6cube_n6.jsonl > runs/s6cube_n6.log 2>&1 &   # 6 min
    python3 search/s6cube.py report runs/s6cube_n6.jsonl
