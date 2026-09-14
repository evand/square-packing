# rank8-margin: the margin of the rank-8 statement, and how small its core is (2026-09-13)

Task: `tasks/rank8-margin/README.md`.  Code: `search/rank8.py` (new; `search/bentz.py` and
`search/leaf_ceiling.py` are read, not touched).  Nothing in `verify*/`, `certificates/`, `lean/`
is touched: this is a measurement.  Runs live in this worktree's `runs/` (gitignored), so every
number is quoted here.

Read against `search/BENTZ.md` §0, §5 (iii), §6, `notes/status.md`, `notes/review-2026-09-13.md`,
`notes/s13-casefree.md` §1 (closed semantics), `notes/proof-anatomy.md` §7.2.

**Semantics.**  `t = 4`; closed unit squares, closed containment, closed container; a packing is a
family pairwise disjoint **as closed sets** (`notes/s13-casefree.md` §1).  **Leaf A** is the fully
pinned Bentz leaf `(AB)^4`: twelve *labelled* squares

| label | pattern (exactly) | what it is |
|---|---|---|
| `T_i`, `i = 0..3` | `{a_i, b_i}` | the four corner squares |
| `C_j`, `j = 0..3` | `{c_j}` | the four mid-wall singletons |
| `D_j`, `j = 0..3` | `{d_j}` | the four interior singletons |

each square confined to its **pattern region**: it contains its own points of `P0` and contains
**none** of the other sixteen (`python3 search/bentz.py points`).

---

## 0. The four answers, up front

1. **The margin is exactly `0`, it is attained, and the maximisers are a large plateau that does
   *not* collapse to the tiling.**  The sup over leaf A of `min_{i<j} gap(S_i, S_j)` (the
   separating-axis signed distance, `< 0` = overlap, `= 0` = touching, `> 0` = disjoint as closed
   sets) is `0`.  The lower bound `>= 0` is **exact**: the sixteen "tiling minus four" families
   have all twelve patterns right in integer arithmetic and exactly `18` of their `66` pairs
   meeting (`rank8.py selftest`).  The upper bound `<= 0` is **measured**: `3,016` multistart local
   optimisations, `2,889` of which ended inside the leaf, never returned a value above
   `+0.000000000000e+00`.  So leaf A is realisable under *open* semantics and fails under *closed*
   semantics by touching only — margin zero, exactly as `notes/status.md` guessed.  What was not
   expected: `2,188` starts reached within `1e-13` of `0`, `1,765` of them further than `0.05` from
   any tiling family (max-norm on centres and on angles mod 90°), the furthest at `0.914`, and
   `174` of them with a square tilted by more than `1°`, **the largest tilt being `32.5242°`**.
   The zero set is a plateau, not a point (§1).

2. **There is no small core: the *only* non-realisable sub-configuration is the whole of leaf A.**
   Over the eight singletons with the four corner squares always present, all `51` `D4` classes
   were measured; over all twelve labels, all `618` classes were resolved.  Every configuration of
   eleven or fewer squares is realisable as a genuine closed packing, and **50 of the 51 witnesses
   were re-verified exactly** (rational snapping, integer pattern tests, exact `sq_meets_sq`).  The
   three eleven-square classes have margins `4.426325e-02` (drop an interior singleton),
   `8.100185e-03` (drop a wall singleton) and — the surprise — **`3.763847e-06`** (drop a *corner*
   square), whose optimum is a genuine exact packing with twelve gaps all equal and a
   **three-square block tilted by `-3.2568°`**.  So the rank-8 statement is *tight*: no rank-5,
   rank-6, …, rank-11 sub-statement carries it, on half the container or anywhere else (§2).

3. **The first-order system around the tiling is infeasible, for all sixteen families, and the
   certificate is a set of wall-to-wall chains of four squares.**  Linearising containment,
   exclusion, admissibility and pairwise disjointness in the centre offsets and tilts
   `(dx, dy, phi)` gives a homogeneous disjunctive LP (40 conjunctive rows, 23 disjunctions over
   the 23 adjacent pairs, 68 branch options), solved as a MILP: value `0` for every one of the
   sixteen tiling-minus-4 families, i.e. **no first-order perturbation opens all twelve gaps**.
   An irreducible 18-row certificate reads as three chains — along `x`: `T2, D1, D2, C3`; along
   `y`: `C1, D0, D1, C2` and `T1, C3, D2, D3, T3` — each a run of unit squares pinned between two
   opposite walls of a side-4 box, which has exactly zero slack.  The linearisation keeps its
   combinatorial structure for `rho + 0.7072 alpha < 0.086000` and its second-order remainder is
   `<= 2 rho alpha + (1 + 2 rho) alpha^2/2 + alpha^2`, whose leading term is **positive** — which
   is exactly why item 1 had to be measured globally and cannot be read off this LP (§3).

4. **Sherali–Adams level 2 on the 162-pose support sees about a tenth of the missing unit.**  Pair
   marginals on disjoint pose pairs, with the twelve clique (pattern) rows lifted by `x_u` and
   `1 - x_u`, give **`11.893347407`** against the clique LP's `12` and the integer optimum
   `alpha = 11`; restricted to the eight singleton labels (the rank-8 statement itself) it gives
   **`7.894639816`** against `8` and `alpha = 7`.  The dual is carried entirely by the four
   *edge-adjacent* interior pairs and the eight wall–interior pairs, and is **exactly zero on the
   two diagonal interior pairs `(D0,D3)`, `(D1,D2)` and on every wall–wall pair** (§4).

**The reading.**  The margin being zero and the core being the whole twelve say the same thing
twice: the rank-8 statement is a *single indivisible* geometric fact with no slack.  A proof of it
cannot be a local perturbation argument (§3 gives only `O(rho alpha + alpha^2)`), cannot be
decomposed into smaller sub-families (§2), and cannot be separated by any pair-level relaxation on
its own (§4 leaves `0.89` of the unit).  What §3 does say is what shape the argument has: it is a
**chain** statement — four unit squares pinned between opposite walls of a side-4 container have
zero slack — and the only thing standing between that and a proof is the strictness that closed
semantics gives away and that Nagamochi / Bentz Corollary 7 recover by hand.

---

## 1. The margin of leaf A

### 1.1 What is being maximised

For two closed unit squares at `(c_i, th_i)` and `(c_j, th_j)`, write `Dc = c_j - c_i` and
`Delta = th_j - th_i`.  The support half-width of a unit square in a direction `d` is
`1/2 (|d.u| + |d.v|)`, so on either square's own axes the combined half-width is the same number
and

    gap(S_i, S_j) = max( |Dc.u_i|, |Dc.v_i|, |Dc.u_j|, |Dc.v_j| ) - W(Delta),
    W(Delta) = 1/2 + 1/2 (|cos Delta| + |sin Delta|).

Separating-axis on the four edge normals is exact for two convex quadrilaterals, so
`gap > 0` iff the squares are disjoint as closed sets, `gap = 0` iff they touch, `gap < 0` iff they
overlap.  The objective is `F = min_{i<j} gap(S_i, S_j)` over the twelve labelled squares, each in
its pattern region.  `rank8.py selftest` checks the sign of this float `gap` against
`leaf_ceiling.sq_meets_sq` (the exact four-axis integer test) on 4,000 random pairs: **4,000 / 4,000
agreements**.

Two caveats recorded once and not repeated.  (a) Non-containment is an *open* condition, so the
region `rank8.py` optimises over is the closed relaxation (exclusion margin `<= 0`); the reported
sup is therefore an upper bound for the true leaf, which is the direction a "`<= 0`" conclusion
needs, and the *witness* at `0` satisfies the exclusions strictly (the `1e-3` nudge).  (b)
Optimisation is `scipy` SLSQP on a smooth sub-problem with the `max`/`|.|` branches frozen,
iterated until the branch assignment stops changing — floating point throughout.

### 1.2 The lower bound, exactly

Nudging a tiling wall square towards the middle of its wall collapses its doubled pattern to the
mid-wall singleton (`BENTZ.md` §5 (i)).  Taking the four corner tiles, the four interior tiles and
**one** of the two tiles of each wall, nudged by `1e-3`, gives `2^4 = 16` configurations.  In exact
rational arithmetic (denominator `1000`, `theta = 0`), through `bentz.pattern_of` and
`leaf_ceiling.sq_meets_sq`:

> **every one of the sixteen tiling-minus-4 families has all twelve patterns exactly right, and
> exactly `18` of its `66` pairs meet as closed sets.**  Its min gap is exactly `0`.

The 18 contacts are `T_i–D_j` (one each), `C_j–D_k` (two each) and all six `D–D` pairs: the four
interior squares are the hub, each with six contacts.  So `sup >= 0` is a theorem, and leaf A fails
to be a packing *only* by touching.

### 1.3 The upper bound, measured

`python3 search/rank8.py margin --nudge 800 --nrandom 1600 --nsupport 600 --nsamp 3000000
--sigma 0.08 --sigma-th 0.12 --support runs/pgonly_corner_exact.txt`
(`runs/rank8_margin.log`, 67 s on 8 processes).  `3,016` starts:

| start family | how many | source |
|---|---|---|
| the sixteen tiling-minus-4 families | 16 | §1.2 |
| nudged tiling (`sigma = 0.08` centre, `0.12` rad angle) | 800 | |
| one random support pose per label | 600 | the 162 poses of `search/pgonly_corner_exact.txt` |
| one random pose per label from rejection-sampled pattern pools | 1,600 | `3e6` random admissible poses, `92k`–`226k` per label |

`2,889` of the `3,016` ended inside the (relaxed) leaf.

> **BEST min pairwise closed gap over leaf A: `+0.000000000000e+00`.**

Nothing positive, ever, from any start.  That is the measurement; it is not a certificate.

### 1.4 The maximisers do **not** converge to the tiling

`2,188` of the runs landed within `1e-3` of `0` — in fact within `1e-13`, the plateau.  Of those:

* `1,765` are further than `0.05` from *any* tiling-minus-4 family, in the max-norm over centres
  and over angles mod 90°; the furthest is at distance `0.914000`;
* the largest tilt of any square over the whole plateau is **`32.5242°`**;
* `174` have a square tilted by more than `1°` mod 90 — e.g. `30.1811°` at distance `0.526759`,
  `16.3302°` at `0.285016`, `12.1601°` at `0.244085`, all at min gap `-5e-16`, i.e. `0`.

So the answer to "does any tilted family get within `1e-3` of `0` far from the tiling" is
**yes, emphatically**: tilted families sit exactly *at* `0`, tens of degrees of tilt away and half
a unit of centre away.  The reason is structural and worth stating plainly: `min gap = 0` says the
twelve squares have **pairwise disjoint interiors**, and the set of interior-disjoint realisations
of leaf A is a large, positive-dimensional set — one square jammed into a corner can rotate freely
while its neighbours give way.  A typical example (`runs/rank8_margin_leafA.txt`, `rank8.py
analyse`) is eleven axis-parallel squares in tile positions with `T3` at `(3.309859, 3.309562)`
tilted `+32.2781°` and jammed into the top-right corner, its single contact being with `D3` pushed
to `(2.414, 2.5)`.

The most tilted one found (`support64` in `runs/rank8_margin_leafA.txt`, min gap `-3.3e-16`, i.e.
`0`; distance `0.526759` from the nearest tiling family) is:

```
  T0 ( 0.50000000,  0.50093490)  0     C0 ( 0.50000000,  2.50000000)  0
  T1 ( 3.50000000,  0.50000000)  0     C1 ( 1.85300142,  0.50000000)  0
  T2 ( 0.50000000,  3.50000000)  0     C2 ( 1.50000000,  3.50000000)  0
  T3 ( 3.50000000,  3.49999996)  0     C3 ( 3.50000000,  2.49651834)  0
  D0 ( 1.50000000,  1.50000000)  0     D2 ( 2.69737837,  1.40241234)  +30.18107 deg   <--
  D1 ( 1.50000000,  2.50000000)  0     D3 ( 2.50000000,  2.58600000)  0
```

— eleven axis-parallel squares and one interior square turned `30.18°` into the pocket the leaf
leaves free, with `16` of the `66` pairs touching.  Nothing about the plateau is near the tiling.

This is the closed/open distinction of `notes/s13-casefree.md` §2 localised to one leaf: under open
semantics leaf A has twelve squares and a whole plateau of them; under closed semantics it has
none, and the deficit is carried entirely by boundary contact.

Saved: `runs/rank8_margin_leafA.txt` (the best 120 configurations, float, one block per
configuration; `rank8.py analyse` prints the contact graph, `rank8.py verify` re-derives patterns
and disjointness exactly after snapping).

---

## 2. The minimal non-realisable sub-configurations

A sub-configuration is **realisable** iff the sup of its min pairwise closed gap is `> 0`, i.e. iff
its squares can be a genuine closed packing inside their pattern regions.  Realisability is
downward closed, so the answer is a pair: the maximal realisable classes and the minimal
non-realisable ones.

**Method.**  `rank8.py subsets` enumerates the `D4` classes, runs the §1 multistart on each with
**basin hopping** (rounds of perturb-and-reoptimise around the best incumbents), and — because a
false "not realisable" is the dangerous error — measures every class exhaustively rather than
stopping at the frontier.  Every *positive* answer is then hardened (re-optimised with a strict
containment and exclusion margin `1e-4`) and re-checked **exactly**: snapped to rationals
(`Q = 1e7` for the angle, `Dc = 1e9` for the centre), patterns re-derived by
`leaf_ceiling.sq_contains` in integers, disjointness by `leaf_ceiling.sq_meets_sq` in integers.

Local optimisation alone gives false negatives here, and did: a first pass with 8 random starts per
class called `{C0,C1,C2,D0,D1,D2,D3}` non-realisable at `+0.000000e+00`, and basin hopping later
found `+8.100185e-03` with an exact witness.  Every "not realisable" below therefore carries the
usual caveat and nothing more; every "realisable" is a theorem about an explicit rational
configuration.

### 2.1 With all four corner squares: the eight singletons

`python3 search/rank8.py subsets --exhaustive --nstart 300 --bh 8 --bhn 600 --confirm 0
--nsamp 1500000 --support runs/pgonly_corner_exact.txt --out runs/rank8_subsets_singletons.txt`
(`runs/rank8_subsets_singletons.log`).  All `51` `D4` classes of the `256` subsets measured.

> **MINIMAL NON-REALISABLE, up to `D4`:  `{C0,C1,C2,C3,D0,D1,D2,D3}` — all eight.  There is no
> other.**
>
> **MAXIMAL REALISABLE:  `{C0,C1,C2,C3,D0,D1,D2}` at `+4.426325e-02` and
> `{C0,C1,C2,D0,D1,D2,D3}` at `+8.100185e-03`.**

The margin table (the best value found for each class, merged and re-polished over two independent
runs — each is a *lower* bound on that class's sup, and each of the `50` proper classes is verified
exactly as a closed packing.  The only class `verify` rejects is the full eight, and it rejects it
for the right reason: correct patterns, `18` meeting pairs):

| `|T|` | class | margin |
|---|---|---|
| 8 | `{C0,C1,C2,C3,D0,D1,D2,D3}` | `+0.000000000e+00` — **not realisable** |
| 7 | `{C0,C1,C2,C3,D0,D1,D2}` | `+4.426324532e-02` |
| 7 | `{C0,C1,C2,D0,D1,D2,D3}` | `+8.100184990e-03` |
| 6 | `{C0,C1,C2,C3,D1,D2}` | `+1.213203436e-01` |
| 6 | `{C0,C1,C2,C3,D0,D1}` | `+7.643406839e-02` |
| 6 | `{C0,C1,C2,D0,D1,D2}`, `{C0,C2,C3,D0,D1,D2}` | `+6.872930441e-02` |
| 6 | `{C0,C1,D0,D1,D2,D3}` | `+3.857127499e-02` |
| 6 | `{C1,C2,D0,D1,D2,D3}` | `+1.813678996e-02` |
| 5 | `{C0,C1,C2,C3,D0}` | `+3.500000000e-01` |
| 5 | the other nine `|T| = 5` classes | `+5.058179e-02` … `+1.308005e-01` |
| 4 | `{C0,C1,C2,C3}` | `+5.000000000e-01` |
| 4 | `{D0,D1,D2,D3}` | `+1.566542267e-01` |
| 4 | the other eleven `|T| = 4` classes | `+8.061024e-02` … `+3.500000e-01` |
| 3 | all ten classes | `+1.566542e-01` … `+5.000000e-01` |
| 2 | all six classes | `+1.720000e-01` … `+5.000000e-01` |
| 1 | `{C0}`, `{D0}` | `+5.000000e-01`, `+4.983324e-01` |
| 0 | `{}` (the four corner squares alone) | `+2.000000e+00` |

The full table is in the log.  Two things it says.  **The four interior singletons on their own are
comfortable** (`{D0,D1,D2,D3}` with all four corners: `0.157`) and so are the four wall singletons
(`0.500`); it is only the *combination* that is tight, which is `BENTZ.md` §5 (iii) in the
continuum rather than on a 162-pose support.  And **the margin decays smoothly and dies exactly at
eight**: `0.157 -> 0.121 -> 0.044 -> 0.008 -> 0`, so there is no "nearly infeasible" seven-square
statement to prove instead.

### 2.2 Over all twelve labels, corner squares included

`python3 search/rank8.py subsets --all …` (`runs/rank8_subsets_all.log`): `618` `D4` classes of
the `4096` subsets, resolved by the downward-closure pruning after `4` tests.

| class | squares | margin |
|---|---|---|
| all twelve | 12 | `+0.000000e+00` — **not realisable** |
| drop an interior singleton `D3` | 11 | `+4.426325e-02` |
| drop a wall singleton `C3` | 11 | `+8.100185e-03` |
| **drop a corner square `T3`** | 11 | **`+3.763847e-06`** |

All three eleven-square witnesses are **exact closed packings** after hardening and snapping
(`runs/rank8_subsets_all_hard.txt`, `rank8.py verify`): the hardened drop-a-corner witness has min
gap `+6.670578e-07` — hardening costs it most of its margin, which is the point — and its snapped
rational poses have all eleven patterns right and no meeting pair, in integers.  Hence

> **the unique minimal non-realisable sub-configuration of leaf A is leaf A itself.**

The third row is the surprise of the task.  Removing an entire corner square — a quarter of the
"pinned" structure — buys a margin of `3.76e-6`, four orders of magnitude less than removing a wall
singleton.  Its optimum (`runs/rank8_drop_corner.txt`) has **twelve pairwise gaps all equal to
`3.763846784e-06`** and is bought by a *rotated block*:

```
  T0 ( 0.50066901,  0.50066497)  th =  -0.0759454 deg
  T1 ( 3.50000000,  0.50000000)  th =  -0.0000000
  T2 ( 0.50000000,  3.50000000)  th =  +0.0000000
  C0 ( 0.52759777,  1.52784234)  th =  -3.2567593 deg   <-- the tilted block
  C1 ( 2.49999624,  0.50000000)  th =  -0.0000000
  C2 ( 1.55967190,  3.50000000)  th =  +0.0000000
  C3 ( 3.49967776,  2.49973320)  th =  -0.0369376
  D0 ( 1.52598650,  1.47103153)  th =  -3.2567593 deg   <--
  D1 ( 1.47098333,  2.47578278)  th =  -3.2567593 deg   <--
  D2 ( 2.57400109,  1.50000376)  th =  -0.0000000
  D3 ( 2.49894690,  2.50035231)  th =  +0.0464645
```

`C0`, `D0`, `D1` rotate together by `-3.2568°`; every other square stays within `0.08°` of axis
parallel.  This is `HONEST.md` §0 item 3 and `BENTZ.md` §0 C1 again — the gain is bought by a
*small* tilt of a *group* of squares, not by one square turning.

### 2.3 Does a rank-5/6 statement on half the container carry the obstruction?

**No.**  Every sub-configuration of eleven or fewer of the twelve squares is realisable, with an
exact witness; in particular every five- or six-square sub-family is realisable with a margin of at
least `1.72e-01` (§2.1's table: with the four corner squares present, the tightest six-square class
is `{D0,D1}` at `+1.720000e-01` and the tightest five-square class is `{D0}` at `+4.983324e-01`;
even eight squares never get below `+8.061024e-02`).  A
half-container statement would therefore have to be a *different* statement — not a sub-family of
leaf A, but a statement about a half-container *region* with its own boundary conditions — and
nothing in the measurement suggests one exists: the certificate of §3 that is closest to a
"half-container" object is a chain of four squares pinned between two opposite walls, which spans
the *whole* container by construction.

---

## 3. The first-order (linearised) system around the tiling

### 3.1 The model

Write the perturbed pose of tile `k` as `c_k = C_k + (dx_k, dy_k)`, `theta_k = phi_k`, with `C_k`
the tile centre of the base family.  To first order, with `e = p - C_k`:

* a point `p` sits in square `k` iff both `coord_u = e_x - dx_k + e_y phi_k` and
  `coord_v = e_y - dy_k - e_x phi_k` have modulus `<= 1/2`.  Containment is the conjunction of the
  four `sg * coord <= 1/2`; **exclusion is their disjunction**, and at the tiling it is a genuine
  disjunction only where several coordinates are tight at once;
* two squares with nominal centre difference `dd` have
  `gap = max_branches [ |dd_ax| + sg (Ddelta_ax + s_o phi_o) ] - 1 - (1/2) |phi_j - phi_i|`,
  over `ax in {u, v}` with `dd_ax != 0`, `sg = sign(dd_ax)`, owner `o in {i, j}`, and
  `s_o = dd_y` for `ax = u`, `-dd_x` for `ax = v`.  **The tilt term does not drop out for a
  diagonal pair**: `dd = (1,1)` gives `max(Ddelta_x + max(phi_i,phi_j), Ddelta_y - min(phi_i,phi_j))`,
  which is where a rotation first buys separation.  The bounding-box cost `(1/2)|Dphi|` is the
  `O(|theta|)` of the brief;
* admissibility `c in [w/2, 4 - w/2]` with `w = |cos| + |sin| ~ 1 + |phi|` becomes
  `dx_k >= (1/2)|phi_k|` at a tile touching the left wall, and its three images.  **Every corner
  and wall tile touches the container**, so these rows are active and load-bearing.

Every `max` entering a `>=` is a disjunction (a binary); every `|.|` entering with a minus sign is
a conjunction (two rows).  The system is **homogeneous**, so it is normalised by
`|dx|, |dy|, |phi| <= 1` and the only question is whether `max g` is `0` or positive.

For the family `pick = 0000`: `40` conjunctive rows, `23` disjunctions over the `23` adjacent pairs
of `66`, `68` branch options; solved as a MILP (HiGHS through `scipy.optimize.milp`) with exactly
one branch selected per disjunction.

### 3.2 The result

    python3 search/rank8.py linear                 # one family, with the certificate
    python3 search/rank8.py linear --all-picks     # all sixteen

> **FIRST-ORDER MILP value `-0.000000000000e+00` for every one of the sixteen tiling-minus-4
> families.  No first-order perturbation of the centres and tilts opens all twelve pairwise gaps:
> the tiling-minus-4 family is first-order rigid inside leaf A.**

The LP is degenerate, so the Farkas dual is not unique; `rank8.py` extracts an **irreducible**
certificate by greedily dropping rows while the value stays `0` (18 rows for the branch the MILP
chooses), and reads it back as chains:

```
  along x: 4 squares {T2, D1, D2, C3}, anchored at T2 (left wall), C3 (right wall)
  along y: 4 squares {C1, D0, D1, C2}, anchored at C1 (bottom wall), C2 (top wall)
  along y: 5 squares {T1, C3, D2, D3, T3}, anchored at T1 (bottom wall), T3 (top wall)
```

Each chain is a run of unit squares pinned between two **opposite** container walls of a side-4
box: four unit widths in a width of four, so exactly zero slack, and a relative tilt only costs
(`(1/2)|Dphi|` per link, plus `(1/2)|phi|` per wall anchor).  Note the second chain is a full
*column* of the tiling — `C1, D0, D1, C2` are the four tiles of `[1,2] x [0,4]` — and the first is a
staircase that still spans wall to wall.  That is the whole first-order obstruction, and it is
precisely Bentz's counting made geometric.

### 3.3 Where the linearisation is valid, and what it does not prove

The model keeps a fixed combinatorial structure: the same pattern rows active, the same pairs
adjacent, the same branch in each disjunction.  With `|centre offsets| <= rho` and `|phi| <= alpha`,
a point of `P0` moves by at most `rho + (sqrt2/2) alpha` in a square's own frame, so the structure
survives while every inactive slack stays positive:

| what must not change sign | slack at the tiling |
|---|---|
| an inactive pattern inequality | `0.086000` (the `y`-coordinate of `C2` holding `c2`) |
| a tile clear of a wall staying clear | `1.000000` |
| a non-adjacent pair staying apart | `1.000000` |

> **the combinatorial structure is unchanged for `rho + 0.7072 alpha < 0.086000`** — e.g. for
> centre offsets under `0.04` and tilts under `2.6°`.

The second-order remainder is

    | gap - gap_lin |  <=  2 rho alpha + (1 + 2 rho) alpha^2 / 2 + alpha^2

(projection error `(1 + 2 rho)(1 - cos alpha) + 2 rho |sin alpha|`, plus the `W(Dphi)` error
`(1/4) Dphi^2 <= alpha^2`), giving

| `rho` | `alpha` | remainder `<=` |
|---|---|---|
| `0.01` | `0.0100` rad (`0.57°`) | `0.000351` |
| `0.02` | `0.0175` rad (`1.00°`) | `0.001166` |
| `0.01` | `0.0500` rad (`2.86°`) | `0.004775` |
| `0.05` | `0.0500` rad (`2.86°`) | `0.008875` |

**This is the honest limit of the first-order argument.**  Since the MILP value is `0` and not
negative, the certificate bounds the true min gap by `O(rho alpha + alpha^2)`, not by `0`, and the
sign of the leading second-order term is *positive* (`W(Dphi) = 1 + |Dphi|/2 - Dphi^2/4 + …`, so a
relative tilt costs less than the linear model says).  The remainder at `alpha = 1°` is `1.2e-3`,
three orders of magnitude larger than §2.2's `3.76e-6` margin on eleven squares — i.e. second-order
effects are exactly the size of the phenomena we are measuring.  A linearised rigidity proof would
therefore have to be done with the exact inequality, not its linearisation; what the LP settles is
only that the obstruction is *first order* (so it does not need a delicate second-order analysis to
*find*), and what shape it has (chains).

---

## 4. Sherali–Adams level 2 on the 162-pose support

`python3 search/rank8.py sa2 runs/pgonly_corner_exact.txt --alpha` (3 s; `runs/rank8_sa2.log`).

**Formulation.**  Take the `162` poses of the certified `PGCORN` measure
(`search/pgonly_corner_exact.txt`); each carries exactly one of the twelve leaf-A patterns
(`BENTZ.md` §0), giving classes of sizes `T0:1, T1:4, T2:2, T3:1, C0:14, C1:14, C2:10, C3:10,
D0:27, D1:24, D2:26, D3:29`.  `7,644` of the `13,041` pairs are disjoint as closed sets (exact
`sq_meets_sq`).  Variables: `y_s` per pose and `y_{st}` per **disjoint** pair (`y_{st} = 0` for
meeting pairs, so those variables are simply absent).  Rows:

* the twelve clique rows `sum_{s in class k} y_s <= 1` (every pose of a class contains the same
  point of `P0`, so a class is a clique — `(P2)` of `BENTZ.md` §1);
* `y_{st} <= y_s`, `y_{st} <= y_t`, `y_s + y_t - y_{st} <= 1`;
* **the level-2 lift**: each clique row multiplied by `x_u` and by `1 - x_u`, for every pose `u`
  outside the class —
  `sum_{s in class k} y_{su} <= y_u` and `sum_{s in class k} (y_s - y_{su}) <= 1 - y_u`.

Objective `max sum_s y_s`.  Every row is valid for 0/1 solutions, so the value is an upper bound on
the integer optimum.

| labels | poses | clique LP | **SA level 2** | integer optimum `alpha` |
|---|---|---|---|---|
| all twelve | 162 | `12` | **`11.893347407`** | `11` |
| the eight singletons `C0..C3, D0..D3` | 154 | `8` | **`7.894639816`** | `7` |

(`alpha = 11` and `alpha = 7` are recomputed here by a complete branch-and-bound on the exact
disjointness graph, independently reproducing `BENTZ.md` §5 (iii) and `rankdiag.py`.)

So **level 2 closes `0.1067` of the missing `1`** on the full leaf and `0.1054` of it on the rank-8
sub-statement — about a tenth.  It is strictly stronger than every degree-1 family in the repo (all
of which sit at exactly `12`), and it is nowhere near enough.

**Which pairs carry the dual.**  Aggregating the positive dual multipliers by label pair
(`404` of `26,508` rows carry a positive dual on the 12-label instance, `393` of `21,721` on the
8-label one):

| label pair | dual mass, 8 labels | dual mass, 12 labels |
|---|---|---|
| `(D0,D1)` | `1.518773` | `1.524602` |
| `(D1,D3)` | `1.475649` | `1.458723` |
| `(D0,D2)` | `1.457876` | `1.436008` |
| `(D2,D3)` | `1.225440` | `1.248316` |
| `(C3,D3)` | `1.057191` | `0.971945` |
| `(C2,D1)` | `1.018508` | `1.049367` |
| `(C0,D1)` | `0.998358` | `0.981580` |
| `(C3,D2)` | `0.943033` | `0.890075` |
| `(C1,D0)` | `0.903376` | `1.065608` |
| `(C1,D2)` | `0.886887` | `0.946198` |
| `(C0,D0)` | `0.846984` | `0.876661` |
| `(C2,D3)` | `0.768184` | `0.794547` |
| `(D0,D3)`, `(D1,D2)` (the **diagonal** interior pairs) | `0` | `0` |
| every `C_i–C_j` pair | `0` | `0` |
| `(T1,T2)`, `(C1,T1)` | — | `2.000000`, `0.666667` |

The structure is exact and worth stating: **the certificate lives on the contact graph of the eight
non-corner tiles** — the four-cycle `D0–D1–D3–D2–D0` of edge-adjacent interior tiles, plus the eight
edges joining each wall singleton to its two interior neighbours.  It puts **zero** weight on the
two interior pairs that meet only at a corner (`D0–D3`, `D1–D2`) and zero on wall–wall pairs.  That
is the same eight-square object §2 shows is indivisible and §3 shows is chained, seen a third way.
(The `(T1,T2)` mass of `2` in the 12-label instance is an artefact of the tiny corner classes —
`T0` and `T3` have a single support pose each, `T2` two — and carries no geometric content.)

---

## 5. What is exact and what is float

**Exact** (integers or `Fraction`; no float decides anything):

* the twelve pattern regions — `bentz.pattern_of` = `leaf_ceiling.sq_contains`, four integer
  half-plane tests in the square's own frame;
* the disjointness of two squares — `leaf_ceiling.sq_meets_sq`, four edge normals in integers;
* §1.2: the sixteen tiling-minus-4 families are in leaf A with `18` meeting pairs
  (`rank8.py selftest`, rational poses with denominator `1000`);
* §2: **every positive answer** — `50` of the `51` singleton classes and all three eleven-square
  classes over the full universe are exhibited as rational configurations whose twelve (or eleven)
  patterns are re-derived in integers and whose pairs are checked pairwise disjoint in integers
  (`rank8.py verify` on the hardened files);
* §4: the disjointness graph of the 162 poses, the class of each pose, and both integer optima
  `alpha = 11`, `alpha = 7` (complete branch-and-bound).

**Float**:

* **the whole of §1's upper bound.**  "The sup is `0`" is `3,016` local optimisations that never
  returned a positive value.  It is evidence, not a theorem;
* **every "not realisable" in §2.**  These are multistart measurements and are known to give false
  negatives — one is documented in §2 — so each is a statement about a search, not about geometry.
  The two conclusions that matter (`sup = 0` on twelve squares; every eleven-square class
  realisable) are, respectively, float and exact;
* **§3 entirely**: the MILP is HiGHS in double precision, and the `0.086` / second-order numbers are
  float evaluations of exact expressions.  Nothing in §3 is quoted as a bound on anything
  certifiable; it is a statement about a linearisation;
* **§4's LP value** `11.893347407` / `7.894639816` is HiGHS in double precision.  The graph it runs
  on, and both `alpha`s, are exact, so the *shape* of the result (`alpha < SA2 < clique LP`) is
  exact and only the third decimal of the value is float.

---

## 6. Reproduce

```sh
python3 search/rank8.py selftest                       # the geometry, against leaf_ceiling's exact tests

# (1) the margin
python3 search/rank8.py margin --nudge 800 --nrandom 1600 --nsupport 600 --nsamp 3000000 \
       --sigma 0.08 --sigma-th 0.12 --support runs/pgonly_corner_exact.txt \
       --out runs/rank8_margin_leafA.txt --nsave 120          # 67 s, 8 processes
python3 search/rank8.py analyse runs/rank8_margin_leafA.txt -n 6

# (2) the sub-configurations
python3 search/rank8.py subsets --exhaustive --nstart 300 --bh 8 --bhn 600 --confirm 0 \
       --nsamp 1500000 --support runs/pgonly_corner_exact.txt \
       --out runs/rank8_subsets_singletons.txt
python3 search/rank8.py subsets --all --nstart 200 --bh 5 --bhn 400 --confirm 2 \
       --nsamp 1000000 --support runs/pgonly_corner_exact.txt --out runs/rank8_subsets_all.txt
python3 search/rank8.py subsets --all --only T0,T1,T2,C0,C1,C2,C3,D0,D1,D2,D3 \
       --nstart 400 --bh 8 --bhn 800 --nsamp 1500000 --support runs/pgonly_corner_exact.txt
# runs/rank8_subsets_final.txt merges and re-polishes the per-class best over both runs
python3 search/rank8.py harden runs/rank8_subsets_final.txt runs/rank8_subsets_final_hard.txt
python3 search/rank8.py verify runs/rank8_subsets_final_hard.txt -n 60     # 50 of 51 EXACT
python3 search/rank8.py harden runs/rank8_subsets_all.txt runs/rank8_subsets_all_hard.txt
python3 search/rank8.py verify runs/rank8_subsets_all_hard.txt -n 6        # all three 11-square EXACT

# (3) the first-order system
python3 search/rank8.py linear                 # one family, irreducible certificate, chains
python3 search/rank8.py linear --all-picks     # all sixteen tiling-minus-4 families

# (4) Sherali-Adams level 2
python3 search/rank8.py sa2 runs/pgonly_corner_exact.txt --alpha
python3 search/rank8.py sa2 runs/pgonly_corner_exact.txt --alpha \
       --labels C0,C1,C2,C3,D0,D1,D2,D3
```

Inputs, read-only: `runs/pgonly_corner_exact.txt` (the certified `PGCORN` measure) and
`runs/inputs-2026-09-12/` (unused here beyond the pattern census, which `bentz.py` already does).
Logs quoted above: `runs/rank8_{selftest,margin,subsets_singletons,subsets_all,linear,sa2,verify}.log`;
configurations: `runs/rank8_margin_leafA.txt` (the 120 best leaf-A configurations, all at `0`),
`runs/rank8_subsets_final{,_hard}.txt` (one per `D4` class of the eight singletons; the hardened
file is the one that verifies exactly), `runs/rank8_subsets_all{,_hard}.txt` (the three
eleven-square classes over all twelve labels), `runs/rank8_drop_corner.txt` (the `3.76e-6`
configuration of §2.2).

**Merging note.**  §2.1's numbers are the best of two independent `subsets --exhaustive` runs plus
a polishing pass (`runs/rank8_subsets_final.txt`); individual runs differ in the fourth significant
figure on a few classes (e.g. `{C0,C1,D0,D1,D2}` at `8.288e-02` in one run and `9.317e-02` in the
other).  Each entry is a *lower* bound on that class's sup, and the qualitative content — which
class is realisable — never differed once basin hopping was on.

---

## 7. What this says about the next step

`notes/review-2026-09-13.md` listed two shapes for the rank-8 statement.  The measurement bears on
both.

* **Analytic / local (linearised rigidity).**  §3 says the obstruction *is* first order and gives
  its shape — chains of four unit squares pinned between opposite walls of a side-4 container, with
  each link paying `(1/2)|Dphi|` for a relative tilt and each wall anchor `(1/2)|phi|`.  That is a
  clean, human-sized statement and it is the right skeleton.  But §3.3 says the linearisation alone
  cannot finish: the second-order remainder at a `1°` tilt is `1.2e-3`, and §2.2 shows margins of
  `3.8e-6` matter.  The proof has to run the chain argument on the exact width
  `1/2 + 1/2(|cos Dphi| + |sin Dphi|)`, which is what Nagamochi's and Bentz's strictness lemmas do
  by hand.  The good news is that a chain argument on the exact width is a *one-dimensional*
  statement about four squares and two walls, which is the smallest object anything in this repo has
  ever reduced the `n = 12` gap to.
* **Degree 2 (Sherali–Adams).**  §4 says level 2 on the optimum's own support is worth about a tenth
  of the unit, on a certificate carried by exactly the tile-contact graph of the eight singletons.
  It is real progress over every degree-1 family (all at `12.000`), and it is not close.  Level 3 or
  a stronger pair object is the only way that direction goes further, and the continuum-crediting
  question (`ALLMEET.md` §0) is untouched by any of it.
* **What is now closed.**  There is no smaller statement to prove.  §2 rules out every
  sub-configuration of eleven or fewer squares by exhibiting an exact packing for each, so the
  rank-8 statement cannot be weakened, split, or tested on half the container first.  And §1 rules
  out any argument that works by finding a positive margin somewhere: the margin is exactly zero,
  it is attained, and it is attained on a plateau reaching `32°` of tilt and `0.9` of centre away
  from the tiling.  Whatever proves the rank-8 statement has to be an argument that survives
  touching — i.e. it must use closed semantics as a hypothesis and not as a convenience.
