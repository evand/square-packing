# `n = 11`: what the mass-11 fractional packing is, and a rigorous ceiling `ν_f(3.84) ≥ 11`

Task: `tasks/s11-anatomy/README.md`.  Read against `search/N11.md` (which found the crossing and
left this question open), `search/DUAL_EXACT.md` and `search/COVER4.md` (the same instrument at
`n = 12`), `search/RANKDIAG.md` and `search/CLIQUE.md` (the diagnostics used here).
Code: `search/n11_anatomy.py`, `search/n11_ceiling.py`, `search/n11_lift.py`,
`search/n11_coverdual.py` (all new; nothing in `verify/`, `verify2/`, `certificates/`, `lean/`,
`TODO.md` was touched).  Runs live in this worktree's gitignored `runs/`, so every number is
quoted here.

---

## 0.  Verdict, up front

**Neither (a) nor (b).  The mass-11 fractional packing is not a smear of anything — it is four
poses, and it is the `3×3` grid with its *middle* square counted three times.**

At `t = 77/20 = 3.85` the exact optimum found here is a D4-symmetrised measure on **four poses /
32 closed unit squares** with

    mass = 11 exactly,     M = max coverage = 1 exactly,     L = mass/M = 11 ,

certified in exact rational arithmetic at every vertex of the arrangement, both in the fundamental
domain (58 vertices) and over the whole container (400 vertices).  Its shape is:

| what | poses | squares | mass each | total |
|---|---|---|---|---|
| the four **corner** squares `[0,1]²`, axis-parallel | 1 | 4 | **1** | **4** |
| **wall** squares, axis-parallel, two overlapping per wall | 1 | 8 | 1/2 | **4** |
| **interior**, tilted `33.70°` | 1 | 8 | 1/4 | 2 |
| **interior**, tilted `4.99989°` | 1 | 8 | 1/8 | 1 |
| | | | | **11** |

So the *ring* — four corners plus one unit of mass per wall — carries exactly **8**, and the
interior carries exactly **3**.  The ring is integrally exact: those 12 distinct squares contain
8 pairwise-disjoint ones (four corners, one per wall, in a pinwheel) and no more.  The interior is
not: its 16 squares carry mass 3, and **once the ring's eight are pinned only two interior squares
survive**.  The independence number of the whole support is

    α(G) = 10  (complete branch and bound),      mass − α = 1.000000  exactly.

That is the answer to the three-way question:

* **(a) Trump's tilted 11-packing smeared — no.**  A smeared 11-packing would have `α = 11`; the
  support has `α = 10`, proved by a complete search, and the measure's heavy poses are
  *axis-parallel* (`72.7 %` of the mass sits at exactly `θ = 0`), which Trump's packing is not.
* **(b) the `3×3` grid plus two in disguise — half.**  The grid's **ring of eight** is there, at
  full mass and integrally realisable; the grid's *ninth* (middle) square is replaced by **three**
  units of interior mass, where only two fit.  It is `8 + 3`, not `9 + 2`.
* **(c) so: something else, and it is the `n = 12` phenomenon in its smallest form yet** — an
  integrality gap of exactly one square, localised to a named region, on 28 distinct squares
  instead of the 1248 of `DUAL_EXACT.md`.

Two results fall out of this that `N11.md` explicitly left unbuilt.

* **A rigorous ceiling for `n = 11`.**  `ν_f(96/25) ≥ 11` exactly (`runs/dual_exact_n11_3.84_support.txt`,
  the same four-pose shape at `t = 3.84`), hence by weak duality (`search/CEILING.md`)
  **no weighted unavoidable set of total weight `< 11` exists at any container side `≥ 3.84`** —
  whatever point set, cell decomposition, LP or verifier is used.  `N11.md` could only say this
  heuristically ("the LP settles at exactly 11"); with the certificate `s(11) ≥ 3040/797 = 3.814304`
  the ceiling `s*(11) = sup{s : COVER(s) < 11}` of the pure method is now bracketed
  **`3.814304 ≤ s*(11) ≤ 3.84`**, a window of `0.0257`.
* **A better upper bound on `COVER` at `3.85`.**  A converged cover LP with exact-verifier
  separation gives an exact-verified cover of weight **`11.1162388`** at `t = 3800/987 = 3.850051`
  (`runs/n11_cd_CD385_cover.txt`, 520 points, `min covered = 1.000004` at `N = 6000` **and**
  `N = 12000`), against `N11.md`'s `11.2608`.  Together with the packing above,
  `11 ≤ COVER(3.850051) ≤ 11.1162388`.

---

## 1.  The measure

`runs/n11_try385.txt`, format of `certificates`-style exact supports (`dual_exact.py`,
`leaf_ceiling.read_measure`): `pose p q cx cy mass`, `θ = 2 arctan(p/q)`, mass `m/8` on each of the
8 dihedral images, `t = 77/20`.

```
pose 0     1      1/2          1/2            4      θ = 0          (1/2, 1/2)
pose 0     1      1/2          7501/5000      4      θ = 0          (1/2, 1.5002)
pose 30287 100000 130193/100000 97989/40000   2      θ = 33.69997°  (1.30193, 2.449725)
pose 2183  50000  313/200      13/8           1      θ = 4.99989°   (1.565, 1.625)
```

`dual_exact.py check` and `check --full`:

| check | vertices | (vertex, square) pairs | mass | `M` | `L` |
|---|---|---|---|---|---|
| fundamental domain `0 ≤ x ≤ y ≤ t/2` | 58 | 245 | `11` | `1` | `11` |
| whole container, no D4 reduction | 400 | 1640 | `11` | `1` | `11` |

The maximum coverage `1` is attained at the container corner `(0,0)` (and, by symmetry, on all of
each closed corner unit box — see the map below).  The orbit of `(1/2,1/2)` has only 4 distinct
images, each counted twice, which is why mass `4` on that pose means mass **exactly 1 on each
corner square**; the other three orbits have 8 distinct images each, so the support is 32 images /
**28 distinct squares**.

The same shape, re-found independently at other containers by the same pipeline:

| `t` | file | poses with mass | mass | `M` | `L` | note |
|---|---|---|---|---|---|---|
| `96/25 = 3.84` | `runs/dual_exact_n11_3.84_support.txt` | 4 | `11` | `1` | **`11`** | `θ = 0, 0, 17.2032°, −0.93705°`; checked full (432 vertices) |
| `77/20 = 3.85` | `runs/n11_try385.txt` | 4 | `11` | `1` | **`11`** | the table above |
| `77/20 = 3.85` | `runs/dual_exact_n11_3.85_support.txt` | 6 | `11` | `1` | **`11`** | polish of a 148-pose union; `θ = 0, 0, 21.400°, 33.800°, 0.015°, 45.000°` |
| `387/100 = 3.87` | `runs/dual_exact_n11_3.87_support.txt` | 4 | `11 − 10⁻⁹` | `1` | `11 − 10⁻⁹` | the `3.85` shape shifted by `+0.01`; `ν_f` is non-decreasing, so `3.87` adds nothing |
| `191/50 = 3.82` | `runs/dual_exact_n11_3.82_support.txt` | 28 | `10.545455` | `1 − 6.25·10⁻¹⁰` | `10.545455` | best exact measure obtained at `3.82`; **not** 11 — see §6 |

The `3.84` and the two `3.85` measures agree on the ring (`4 + 4`) to the digit and disagree only
on how the interior's `3` is split between angles — the LP is degenerate there.  That the interior
composition is arbitrary and the ring is not is the whole content of the object.

### Coverage map (`t = 3.85`, `39×39` cell centres; `#` = coverage in `[0.95,1]`, digit = tenths)

```
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########5577#5555#5555#7755##########
##########577##7555#5557##775##########
##########5#####775#577#####5##########
555555555555555552303255555555555555555
555555557#5777777377737777775#755555555
555555577#57777777###77777775#775555555
55555557##57777777###77777775##75555555
5555555###57777777#5#77777775###5555555
555555557#57777778###87777775#755555555
555555555757777778757877777757555555555
555555555723777885676588777327555555555
555555555537####76###67####735555555555
##########07##5#57###75#5##70##########
555555555537####76###67####735555555555
555555555723777885676588777327555555555
555555555757777778757877777757555555555
555555557#57777778###87777775#755555555
5555555###57777777#5#77777775###5555555
55555557##57777777###77777775##75555555
555555577#57777777###77777775#775555555
555555557#5777777377737777775#755555555
555555555555555552303255555555555555555
##########5#####775#577#####5##########
##########577##7555#5557##775##########
##########5577#5555#5555#7755##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
##########555555555#555555555##########
```

The four `1.0` blocks are the corner squares; the `0.5` bands are the wall squares, one layer deep
everywhere except the short overlap where the two squares of a wall meet; the interior is a
`0.7–0.9` haze.  Compare the `t = 3.99`, `n = 12` map of `search/DUAL.md`, where the corners are at
`0.86` and nothing is integral.

---

## 2.  Anatomy (`search/n11_anatomy.py`)

**Mass per region** (the coverage integral `∫_R cov = Σ_i m_i·area(Q_i ∩ R)`, which sums to the
mass; areas by exact convex polygon clipping).  `t = 3.85`, 4-pose measure:

| region | count | area each | mass | share | mean coverage |
|---|---|---|---|---|---|
| corner box `[0,1]²` | 4 | `1` | **`4.0000`** | `36.4 %` | **`1.0000`** |
| wall strip `[1, 2.85] × [0,1]` | 4 | `1.85` | `4.5181` | `41.1 %` | `0.6106` |
| centre box `[1, 2.85]²` | 1 | `3.4225` | `2.4819` | `22.6 %` | `0.7252` |

(the 6-pose measure: `4.0000 / 4.4736 / 2.5264`.)  The corner boxes are **saturated**: coverage is
identically `1` on each closed `[0,1]²`.  Nothing at `n = 12` does this — at `t = 3.99` the four
corner boxes take `3.44` of `12.008` at mean coverage `0.86` (`DUAL_EXACT.md`, re-measured here
with the same script), which `search/BRANCH.md` singles out as "something no integral packing does".

By the pose's own region label (`leaf_ceiling.region_boxes`, `r = 1`), the split is exact:

    C0 = C1 = C2 = C3 = 1      W0 … W7 = 1/2      I = 3          total 4 + 4 + 3 = 11

**Angle histogram** (mass by `|θ|`; the 8 dihedral images of a pose carry `±θ` and `90° ∓ θ`):

| bin | mass | share | poses |
|---|---|---|---|
| `[0,5°)` | `9.0000` | `81.8 %` | 3 |
| `[30,35°)` | `2.0000` | `18.2 %` | 1 |
| `|θ| < 0.01°` | `8.0000` | `72.7 %` | — |
| `|θ − 45°| < 1°` | `0.0000` | `0 %` | — |

For the 6-pose measure: `8.5` within `0.02°`, `1` at `21.400°`, `1` at `33.800°`, `0.5` at
`45.000°`.  So between `77 %` and `82 %` of the mass is axis-parallel or within five degrees of it,
and **the `45°` diamond that carries `5.3 %` at `n = 12` is absent or marginal here**.  The
rotation that the LP is buying at `n = 11` is one tilted pair (`17°–34°`, container-dependent) plus
a sub-degree nudge — exactly the `HONEST.md` item-3 mechanism (a grid square nudged and tilted a
fraction of a degree to dodge the grid-line points), and nothing else.

**Concentration.**  `25 %` of the mass on 1 pose, `50 %` on 2, `75 %` on 3, `99 %` on 4.  At
`n = 12`, `t = 3.99`: `25 %` on 1, `50 %` on 3, `75 %` on 24, `99 %` on 129.

---

## 3.  The largest pairwise-disjoint subfamily (`search/rankdiag.py --step 1`)

Closed-intersection graph of the 32 images (exact integer separating-axis test, touching counts),
complete branch and bound:

    n = 32 squares, 112 edges (density 0.226),  mass = 11.000000
    α(G) = 10   (complete: True)                gap = mass − α = 1.000000

A maximum independent set — a genuine packing of **ten** unit squares in `[0, 3.85]²`:

    (0.500,0.500) (3.350,0.500) (0.500,3.350) (3.350,3.350)    four corners,        θ = 0
    (0.500,1.500) (3.350,2.350) (2.350,0.500) (2.350,3.350)    one per wall,        θ = 0   (a pinwheel)
    (1.400,2.548)                                              interior,            θ = 56.30°
    (2.225,1.565)                                              interior,            θ = 85.00°

**Where the missing unit lives.**  Restricting to each family and running the same complete search:

| family | squares | mass | `α` |
|---|---|---|---|
| corner squares | 4 distinct | `4` | **`4`** |
| wall squares | 8 | `4` | **`4`** |
| **ring = corners ∪ walls** | 12 distinct | **`8`** | **`8`** |
| interior, `33.70°` | 8 | `2` | `4` |
| interior, `4.99989°` | 8 | `1` | `1` |
| **interior = both** | 16 | **`3`** | `4` |
| **everything** | 28 distinct | **`11`** | **`10`** |

The ring is tight (`mass = α = 8`) and the interior alone is slack (`mass 3 ≤ α 4`): **the whole
gap is in the interaction**.  Pin a maximum independent set of the ring — the four corner squares
and one square per wall — and exactly 10 of the 16 interior squares remain disjoint from all eight;
their independence number is **2** (complete search), so the largest packing extending the ring is
`8 + 2 = 10`, and the LP puts `8 + 3`.

> **The statement a proof needs at `n = 11`.**  *With the four corner unit squares and one unit
> square per wall placed (eight squares, which fit), at most two further pairwise-disjoint unit
> squares fit in `[0,t]²` for `t` below `s(11)`.*
>
> This is a rank-10 inequality `μ(X) ≤ 10` on 28 explicit squares.  Like the rank-8 statement that
> `notes/review-2026-09-13.md` localises at `n = 12`, it is true integrally and inexpressible by any
> single-square (degree-1) inequality.  Unlike it, it is a pure *packing* statement — no incidence
> pattern, no labelled points — and it lives on 28 squares confined to one region.

**What the certifiable rank sub-families are worth here** (`rankdiag.py --step 2/3/4`,
`runs/n11_rankdiag_385.json`):

* chordless `C5` with mass `> 2`: **0** (heaviest anywhere `1.625`); `C7` `> 3`: **0** (heaviest
  `2.25`); 7-antiholes `> 2`: **0**; 5-wheels `> 2`: **0** — the same verdict as `RANKDIAG.md` §2.
* the five-anchor pentagon family that *did* bite at `n = 12` (`0.076–0.436` of the gap there):
  best `μ = 2.000000` against the rank bound `2`, **excess exactly `0`**.  It is worth nothing here.
* size-vs-excess front: the smallest violated rank subset has `|X| = 11` (excess `+0.125`, `12 %`
  of the gap); `|X| = 16` reaches `+0.5`; `|X| = 24` reaches the whole `+1.0`.  A recurring
  `|X| = 8` witness has mass `2.25` against `α = 2` (excess `0.25`) and lives in
  `{I 1.25, W6 0.5, W7 0.5}` — the interior plus two wall squares, angles `{0, 33.7, 56.3, 85}`.
* **cliques do bite** (`search/clique_check.py`): the maximum-mass clique of the overlap graph has
  **mass `1.250000` on 8 poses with no common point** — non-Helly, so the point-cover LP does not
  have it.  All eight are interior (`33.70°/56.30°` pair plus six `4.99989°/85.00°` squares).  At
  `n = 12`, `t = 3.99` the same script gives `1.327` on 221 poses (`CLIQUE.md`).

---

## 4.  The cover side

`search/n11_coverdual.py` runs the `n = 11` cover LP (`tighten.py`'s cutting-plane loop with the
exact verifier as separator, 12 column-generation rounds) and dumps its dual — the fractional
packing on placements that blocks the method.

| tag | `t` | iterations | LP value | converged? | dual |
|---|---|---|---|---|---|
| `CD382` | `15200/3979 = 3.820055` | 43+ (still running at the stop) | **`11.000022`** from `it19` on, unmoved | no (`probe_min ≈ 0.86`) | — |
| `CD385` | `3800/987 = 3.850051` | 40 | **`11.1162061`** | **yes** (`probe_min = 1.0000010`, 0 violated placements) | mass `11.1162`, symmetrised coverage `1.114` (not feasible: the point set is not fully priced out) |

`11.000022 = 11·(1 + 2·10⁻⁶)` is `11` times the LP's own coverage margin, i.e. the LP value is
exactly `11`.  `N11.md` saw this plateau at `D = 3983, 3982, 3980, 3975` (`3.8162–3.8239`); it holds
at `3.820055` too, over 25 consecutive iterations while 15k cuts per round were added.  At
`3.850051` the same loop passes through exactly `9.000018` (`it3`, `it4`) and exactly `11.000022`
(`it7`–`it13`) before settling at `11.1162`.

`CD385`'s cover was re-run through the verifier independently of the driver:

    verify runs/n11_cd_CD385_cover.txt 11 6000  2 0  ->  min covered = 10000039/10000000 = 1.000004
    verify runs/n11_cd_CD385_cover.txt 11 12000 2 0  ->  min covered = 10000039/10000000 = 1.000004

so it is an exact-verified cover of weight `111162388/10⁷ = 11.1162388`: **`COVER(3800/987) ≤ 11.1162388`**
(it is of course not a certificate — the weight is above 11).

**Where the cover's weight sits.**  Of the `11.108` carried by the 344 points of the `CD385` probe
at `it14`, `69.8 %` lies within `0.005` of one of the lines `x` or `y ∈ {1, t/2, t−1}`, `76.2 %`
within `0.02`, and `97.5 %` within `0.1`.  These are the wall-square corner lines and the mid-lines
— the `n = 11` analogue of `RANKDIAG.md`'s finding that every violated pentagon at `t = 4` is
anchored on `x ∈ {1,2,3}`.  The dual picture is the same object seen from the other side: the
packing's ring squares have their corners on exactly those lines.

---

## 5.  Comparison with `n = 12`

| | `n = 12`, `t = 3.99` | `n = 11`, `t = 3.85` |
|---|---|---|
| certified measure | mass `12.008231`, `M = 1 − 2.4·10⁻⁹` | mass `11`, `M = 1` exactly |
| support | 156 poses / 1248 squares | **4 poses / 28 distinct squares** |
| `α` of the support | `11` | `10` |
| integrality gap | `1.008231` | **`1.000000`** exactly |
| mass on a corner square | `0.85` | **`1`** (saturated) |
| corner boxes, mean coverage | `0.86` | **`1.00`** |
| axis-parallel share | `57 %` (`|θ| < 1°`) | `73 %` (`θ = 0` exactly) |
| `45°` diamonds | `5.3 %` | `0 %` |
| max-mass non-Helly clique | `1.327` (221 poses) | `1.250` (8 poses) |
| pentagon rank-2 family | `+0.076 … +0.436` | **`+0.000`** |
| where the missing unit lives | 8 labelled squares, a rank-8 *incidence* statement | the interior of the ring, a rank-10 *packing* statement |

The `n = 11` object is the same phenomenon, one square smaller, and enormously cleaner: the
corner and wall parts are integral, so the entire relaxation error is a single named sub-problem.

---

## 6.  How far down does it go?  (the ceiling window)

`ν_f` is non-decreasing, and `COVER(3040/797) ≤ 10.8146708` is certified (`runs/n11_H3985.txt`,
`N11.md`), so `ν_f(3.814304) ≤ 10.8147 < 11`.  The threshold
`t₁₁ = min{t : ν_f(t) ≥ 11}` therefore satisfies

    3.814304  <  t₁₁  ≤  96/25 = 3.84       (the right end is what §1 certifies).

`search/n11_ceiling.py` searches the pose parameters of this family directly (Nelder–Mead over the
`3k + 1` numbers, with the mass LP over the exact arrangement inside).  With the corner pose, the
wall pose and **two** free interior poses it reaches mass `11` down to a bracket
`[3.849063, 3.849531]` and no lower; with **four** free interior poses it reaches `11` at `3.8450`
and `3.8400` and falls to `10` at `3.8300` and `3.8200`.  The objective is flat (piecewise constant,
with a wide plateau at exactly `10`), so these are upper bounds on what the family can do, not
lower bounds on `t₁₁`; a downward continuation with more free poses was still running at the stop
(`runs/n11_descend_k4.out`, `runs/n11_descend_k5.out`).  At `3.82`, the best exact measure obtained
by any route here is `10.545455` (`runs/dual_exact_n11_3.82_support.txt`), so **the `0.026` window
`[3.8143, 3.84)` is not closed** and the cover LP's plateau at exactly `11` from `3.8162` up is
still only heuristic evidence that `t₁₁ ≈ 3.8153`.

---

## 7.  Which `n = 12` instrument to point at `n = 11` next, and what it should buy

**The box-clique certificate (`search/BOXCLIQUE.md`, `search/boxclique.py`, `search/cliquelever.py`),
after one `search/clique_ceiling.py` run to price it.**  Three reasons, all with numbers from above.

The corner branch — `search/branch.py`, `search/BRANCH.md`, `search/leaf_ceiling.py --corners 1111`
— is the *wrong* lever here, and for the opposite reason to the one the brief anticipated.  The
corner boxes are not empty; they are **saturated**: the unconstrained optimum already puts mass
exactly `1` on each corner square and coverage identically `1` on each corner box, so the leaf
`k = 4` is not a restriction at all and the branch is vacuous.  At `n = 12` it was a restriction
(corner mass `0.85`, and pinning it to `1` moved `T4SCREEN`'s leaf value), which is why it was worth
running there and is worth nothing here.  Nor is the pattern/incidence tree of `BENTZ.md` any use:
there is no point set to branch on, and the missing unit is a disjointness statement about three
interior squares.

Cliques, by contrast, bite immediately and in exactly the right place: the optimum violates a
non-Helly clique inequality by **`0.25`** (`μ(K) = 1.25` on 8 interior poses with no common point),
and all eight members are interior squares — the family where the whole gap lives.  A box clique is
already an object the verifier checks (`cores_meet`, exact `i128` separating axis), already in
`certificates/FORMAT.md`, and already Lean-proved (`packing_le_weight_cliques`), so nothing new has
to be trusted.  The hard ceiling on what *any* rank-style family can ever be worth on this measure
is `mass − α = 1.000000` — a full square, against `0.32–0.79` at `n = 12` — and, unlike at `n = 12`,
every fraction of it converts directly into container side: the pure method dies at
`COVER(t) = 11`, so a clique-strengthened cover of weight `< 11` at `t` *is* a certificate for
`s(11) ≥ t`, and the window between the record `3.814304` and the ceiling `3.84` is `0.026` wide
with the plateau sitting flat at exactly `11` across it.

Concretely: run `clique_ceiling.py` at `t = 3.82, 3.83, 3.84` for `n = 11` to measure the
clique-strengthened packing optimum (if it is `< 11` at some `t` in the window, a clique certificate
exists there), then `boxclique.py` to build one.  Expected purchase: the first `0.25` of the missing
`1.0` is visible now; whether that is enough to move `s(11)` from `3.8143` towards `3.84` is
precisely what the ceiling run measures, and it is a one-afternoon experiment because every
container in the window is `< 4` where the sweep verifier, `branch.py --n 11`, the cliques and the
polygons all already work.  The `0.75` the clique family cannot reach is the rank-10 statement of
§3, and that is mathematics, not another loop — but it is a much smaller and more concrete piece of
mathematics than the `n = 12` rank-8 statement, because the eight ring squares are *pinned by the
optimum itself*, not by a case split.

---

## 8.  Rigorous vs heuristic

**Rigorous** (exact integers and `Fraction`s end to end; the only floats *choose* the poses and
masses, and every chosen object is re-checked exactly):

* `ν_f(96/25) ≥ 11` and `ν_f(77/20) ≥ 11`, hence `COVER(t) ≥ 11` for every `t ≥ 96/25 = 3.84`
  (`dual_exact.py check` and `check --full` on `runs/dual_exact_n11_3.84_support.txt` and
  `runs/n11_try385.txt`; the argument is `DUAL_EXACT.md`'s, unchanged — upper semicontinuity of a
  finite closed-square coverage puts its maximum at an arrangement vertex, and every vertex is
  enumerated in integers).
* `COVER(3800/987) ≤ 11.1162388` (`runs/n11_cd_CD385_cover.txt`, exact verifier at `N = 6000` and
  `N = 12000`, `min covered = 10000039/10000000`).
* `α(G) = 10` and `mass − α = 1` for the `3.85` support, and every entry of the family table in §3:
  integer separating-axis tests and a complete branch and bound, masses compared over the common
  integer denominator.  Likewise the `0` counts for `C5`, `C7`, antiholes, wheels and the pentagon
  excess `0`.
* The max-mass clique `1.25`: the overlap graph is built with a strict separating-axis test and the
  maximum-mass clique found by complete branch and bound; the masses are exact dyadic rationals.

**Heuristic**: which poses are in the measure at all (Nelder–Mead, HiGHS, the column generation of
`packing_dual.py`/`tighten.py`); every "LP value" in §4, including the plateau at exactly `11` at
`3.820055` (rows incomplete — the value can only rise with more rows, but the column set was cut
off after 12 pricing rounds at symmetrised dual coverage `1.06`, so it is not a bound in either
direction); the claim that `t₁₁ ≈ 3.8153`; the claim that the four-pose shape is *the* optimum
rather than *an* optimum (the LP is degenerate on the interior — three different interior splittings
of the same `3` were found at `3.84`, `3.85` and `3.85`).  Nothing about `s(11)` itself is claimed:
`s(11) ≥ 3040/797` remains the record and the measures here are obstructions, not packings of 11
squares.

---

## 9.  Reproduce

```sh
cd verify && cargo build --release && cd ..

# the statement:  nu_f(3.85) >= 11  and  nu_f(3.84) >= 11
python3 search/dual_exact.py check runs/n11_try385.txt --t 77/20 --tag TRY385           # 58 vertices
python3 search/dual_exact.py check runs/n11_try385.txt --t 77/20 --tag TRY385 --full    # 400 vertices
python3 search/dual_exact.py build --t 96/25 --tag n11_3.84 --src runs/n11_cand_3.84.txt --procs 2
python3 search/dual_exact.py check runs/dual_exact_n11_3.84_support.txt --full --tag n11_3.84

# the anatomy
python3 search/n11_anatomy.py runs/n11_try385.txt --t 77/20 --map 39
python3 search/rankdiag.py runs/n11_try385.txt --step 4 --time 200 --cands 400 --pent 60
python3 search/clique_check.py runs/n11_try385.txt --time 300

# the cover side (25-40 min each at 2 threads; use taskset to cap)
python3 search/n11_coverdual.py runs/n11_seed_union2.txt CD385 --Dp 3948 --colgen 12 --cg-want 150 --threads 2 --n 11
python3 search/n11_coverdual.py runs/n11_seed_union2.txt CD382 --Dp 3979 --colgen 12 --cg-want 150 --threads 2 --n 11
verify/target/release/verify runs/n11_cd_CD385_cover.txt 11 6000  2 0
verify/target/release/verify runs/n11_cd_CD385_cover.txt 11 12000 2 0

# how the poses were found in the first place
python3 search/packing_dual.py 77/20 N385 --threads 2 --rounds 40 --time 4500 --polish-after 16 \
        --seed-pitch 0.06 --seed-dth 5 --row-pitch 0.025 --warm runs/dual_PA2_support.txt --final-full
python3 search/n11_lift.py --t 77/20 --out runs/n11_union_3.85.txt runs/n11_smear_3.85_support.txt \
        runs/n11_smear_3.82_support.txt runs/n11_try385.txt
python3 search/dual_exact.py build --t 77/20 --tag n11_3.85 --src runs/n11_union_3.85.txt --procs 2
python3 search/n11_ceiling.py --t 77/20 --restarts 6            # re-derives the four poses
python3 search/n11_ceiling.py --bisect 3.82 3.85 --steps 6 --restarts 10
```

`runs/n11_try385.txt` and `runs/dual_exact_n11_3.84_support.txt` are the two load-bearing files;
copies that survive the gitignored `runs/` are `search/n11_exact_3.84_support.txt` and
`search/n11_exact_3.85_support.txt`.  `dual_exact.py check` reads only the support file and Python
integers, so the statement can be re-checked from those two files alone.
