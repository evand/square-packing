# Level-2 branching: wall-slot occupancy (task F design, 2026-08-30)

Code written for this note (all new, nothing in the certificate pipeline touched):

| file | what it does |
|---|---|
| `search/level2_regions.py` | splits a fractional packing measure (branch dual, support file, exact support) into corner-box / wall-slot / interior mass by the **centre** of the pose |
| `search/level2_capacity.py` | `alpha(R)`: maximises the separating-axis margin of `k` admissible poses with centres in given boxes, i.e. searches for `k` squares of a packing centred in a region |
| `search/level2_leaves.py` | leaf counts up to the symmetries of the container (Burnside over the eight slots) |
| `search/level2_lp.py` | packing-side leaf value: full-container LP (poses x points) with per-region occupancy equalities, row generation against the exact arrangement-vertex certification of `packing_dual.py`, and reduced-cost grid pricing |
| `search/level2_cover_lp.py` | cover-side leaf value on the `k = 4` run's own saved instance (its point columns and its dual's pose rows), with per-region multipliers — directly comparable with the published `12.000024` |

Runs are in this worktree's `runs/` (`l2_*.log`, `cap_*.txt`, `cover98*.txt`); `runs/` is
gitignored, so they are not in the commit.

---

## 0. Verdict, up front

1. **Branching on the wall-slot *total* is provably useless, and so is branching per wall.**  The
   `k = 4` leaf's own optimum at `t = 3.98` (`runs/branch_t398hk4_dual_it16.txt`) has corner mass
   exactly **4**, wall-slot mass exactly **4**, interior mass exactly **4**, and exactly **1** per
   wall strip.  A shared multiplier over the eight slots (9 D4-symmetric leaves) or one per wall
   (four regions) therefore has that same measure feasible in the leaf `m = 4` / `(1,1,1,1)`, and
   cuts nothing there.  Do not build either.
2. **What is fractional is the per-slot mass: exactly 0.5 in each of the eight slots.**  Only a
   per-slot, asymmetric branch moves the LP off its optimum.
3. **The tree is small — 15 leaves — because of a classical lemma, not because of the boxes.**
   Nagamochi's Lemma 7(i) / Stromquist's chord lemma says at most **3** squares of a packing have
   centres within 1 of a given wall (each cuts the line at distance 0.9 from that wall in length
   `> 1`, the chords are disjoint, the wall has length `t <= 4`).  In the `k = 4` leaf two of those
   three are the wall's corner squares, so **at most one of a wall's two slots is occupied**.  That
   leaves `3^4 = 81` occupancy patterns, **15 up to the symmetries of the square**
   (`search/level2_leaves.py`), all with `m = sum K_j <= 4`.  Numerically confirmed: four squares
   centred in `[0, 3.98] x [0, 1]` miss by margin `-0.0067` (`runs/cap_wallstrip_k4.txt`), which is
   exactly the `(4 - 3.98)/3` the chord count predicts.
4. **`packing_le_weight_regions` is already enough** for the leaf statement: it takes an arbitrary
   finite family of pose predicates, arbitrary real multipliers and arbitrary counts.  No new Lean
   theorem is needed for the reduction.  What *is* new is the chord lemma that makes the 15 leaves
   exhaustive — a strictness statement (`> 1`, true for side `L > 1`, false for unit squares), so
   it cannot be produced by any LP and has to be supplied analytically (`notes/proof-anatomy.md`
   §7.4 already lists it as such).  A branch that avoids it costs 43 leaves instead of 15 (§2.5).
5. **Numbers.**  §4.  On the `k = 4` leaf's own cover-side instance the branch is worth
   **0.556 on the hardest leaf** and up to 2.0 on the easier ones, against a control that reproduces
   the published `12.000024` exactly (against **0.00** for the wall-total and per-wall branches, which
   is a structural fact, not a measurement); on the packing side, an independent LP with the opposite
   error direction, the four `m = 4` leaves cut **0.49–0.82** off their own control, in the same
   order.  At **`t = 4.0`** (closed semantics) the same leaves cut **0.51–0.90** — the cut does not
   shrink at the container itself.  All of these are optimistic (§4.1), but the `k = 4` leaf at 3.98
   is open by `2.4e-5` and the cuts are four orders of magnitude larger.
6. **Can level-2 branching reach `t = 4`?  Not on the evidence here, and probably not alone.**  The
   cuts in item 5 are large, but they are measured on an instance whose pose rows are the 520 orbits
   the `k = 4` run happened to keep, and the LP buys most of the cut by setting an empty slot's
   multiplier to `-1` — which against the *full* pose space (a verifier as separation oracle) is far
   more expensive than against 131 sampled poses per slot.  How much survives is exactly what has
   not been measured, and the level-1 precedent is not encouraging: the leaves an integral packing
   resembles (corners full) gained nothing at level 1 (`BRANCH.md`).  Against that, `t = 4` needs
   the pure `12.163` brought below 12 in *every* leaf of a tree that, at `t = 4`, has to hang under
   every corner leaf and not just `k = 4`.  The clique columns are the part that scales there
   (`CLIQUE_CEILING.md`: the clique relaxation sits at 11.8–12.0 at `t = 4` against a pure 12.163),
   and §1 shows the two are complementary rather than redundant — both drive the *frame* to `4 + 4`
   and differ only in what they do to the interior.  **The realistic target for this tree is
   `s(12) >= 3.98`, where the `k = 4` leaf is open by `2.4e-5`; `t = 4` wants cliques as well.**

---

## 1. Input: where the mass sits

`level2_regions.py` splits a D4-symmetrised measure by the centre of the pose into four corner
boxes, eight wall slots and the interior (§2.1).  These partition the container, so the columns
add up to the total mass.

| measure | `t` | total | corners | wall slots | interior | per slot | per wall strip (2 corners + 2 slots) |
|---|---|---|---|---|---|---|---|
| `branch_t398hk4_dual_it16` (the `k = 4` leaf dual) | 3.98 | 12.000 | **4.000** | **4.000** | **4.000** | **0.500** | **3.000** |
| `dual_PD1_support` (pure, certified) | 3.98 | 11.918 | 3.436 | 4.804 | 3.678 | 0.600 | 2.918 |
| `dual_PA2_support` (pure, certified) | 3.99 | 12.008 | 3.400 | 4.808 | 3.800 | 0.601 | 2.902 |
| `dual_PC1_support` (pure, certified, closed) | 4.00 | 12.163 | 3.200 | 5.400 | 3.563 | 0.675 | 2.950 |
| `cc_C99c_exact` (clique-feasible, certified) | 3.99 | 11.364 | **4.000** | **4.000** | 3.364 | 0.500 | **3.000** |

Three readings.

* **The `k = 4` leaf optimum is integral on every aggregate and half-integral on every slot.**
  `4 / 4 / 4`, `1` per wall strip, `0.5` per slot.  The 0.5 is forced by D4 symmetry once the wall
  total is 4, but the wall total being 4 is not: the corner constraint alone drives it there.  And
  `3` per wall strip is exactly the chord-lemma bound (§2.3) — the leaf optimum is sitting on that
  constraint, which is why the per-wall aggregate has no slack to give.
* **The clique relaxation does the same thing to the frame and then eats the interior.**  At 3.99
  the certified clique-feasible measure is `4 / 4 / 3.364`: corners and wall slots integral in
  aggregate, all of the residual excess over 11 interior.  Branching and cliques push on the same
  frame; they differ in what they do to the interior.
* **Inside a slot the mass is not one object.**  Of the 0.5 per slot at 3.98, about 0.33 is the
  axis-parallel wall square (`(0.5, 1.5, 0 deg)` with orbit mass 1.94, plus `(0.5, 1.5238, 0 deg)`
  with 0.68) and about 0.17 is tilted poses leaning on the wall.  A visible sub-cluster of those
  (orbit masses 0.035–0.043 at `(0.577, 2.090, 9.7 deg)` and neighbours, about 0.15 of the 4.0 wall
  total, i.e. 4 %) sits at `cy ~ 2.0 = t/2` — the **mid-wall**, exactly where the slot boundary is.
  A design that leaves a gap there (§2.5, Design B) gives that mass somewhere to hide.

---

## 2. The regions

### 2.1 Definition — membership by the centre

Container `[0,t]^2`, corner size `r` (default `r = 1`), `h = t/2`:

```
C_0 = [0,r]^2              C_1 = [t-r,t] x [0,r]      C_2 = [0,r] x [t-r,t]   C_3 = [t-r,t]^2
W_0 = [r,h]   x [0,r]      W_1 = [h,t-r] x [0,r]                     (bottom wall)
W_2 = [t-r,t] x [r,h]      W_3 = [t-r,t] x [h,t-r]                   (right wall)
W_4 = [h,t-r] x [t-r,t]    W_5 = [r,h]   x [t-r,t]                   (top wall)
W_6 = [0,r]   x [h,t-r]    W_7 = [0,r]   x [r,h]                     (left wall)
I   = [r,t-r]^2
```

`W_0 .. W_7` are in cyclic order round the boundary; a 90-degree rotation is the shift by two, so
the container's D4 acts on the slots as the order-8 subgroup of the octagon's dihedral group
generated by that shift and the mirrors.  This is the Nagamochi/Bentz partition of `[0,4]^2`
(corner unit squares, wall strips of height 1, interior `[1,3]^2`; `notes/proof-anatomy.md` §3.4)
with each wall strip cut at its mid-line, and it is the partition the trivial `4 x 4` grid
respects: 4 + 8 + 4 = 16 centres, one per region.

### 2.2 Membership by the centre, or by "contains a wall point"?

| | centre box `R = {(c,th) : c in B}` | point region `R_p = {(c,th) : p in S(c,th)}` |
|---|---|---|
| capacity `<= 1` | needs a geometric bound (§2.3) | **free** — `R_p` is a clique, and `card_filter_clique_le_one` in `lean/Sqpack/Basic.lean` already proves the count is `<= 1` |
| verifier support | **exists**: the cell's bounding box against an axis-parallel box, exact for convex cells (`certificates/FORMAT.md`, "Branch certificates") | new: needs "every pose of this cell contains `p`" (`p` in the cell's exact bin core, `search/zeromargin.py` §2) and "some pose misses `p`" (`p` outside the cell's hull); both implementable, neither implemented |
| Lean | `packing_le_weight_regions` | `packing_le_weight_regions` (same) |
| catches the mass? | yes: a slot's 0.5 is a wall square *plus* a tilted mid-wall cluster, and the box holds both | only the poses through `p`; the tilted cluster escapes unless a second point is added |
| relation to cliques | none | `R_p` with a **negative** multiplier *is* a point-clique column |

**Recommendation: centre boxes now, point regions later.**  The verifier already does the centre-box
test, so the level-1 machinery generalises from four hard-wired corner boxes to a list of boxes
almost verbatim.  Point regions are the more elegant object — their `<= 1` is a theorem already in
the Lean file rather than new geometry — and they are the natural meeting point with the clique
work; they should follow once the exact-core test of `zeromargin.py` is ported into `verify/`.

### 2.3 Capacity `alpha(R)`

`search/level2_capacity.py` maximises the disjointness margin

```
margin(P1,P2) = max over the four edge normals n of ( |<c2-c1, n>| - (1 + w(delta))/2 ),
w(a) = |cos a| + |sin a|,  delta = th2 - th1
```

over admissible poses with centres in given boxes (`margin >= 0` iff the closed unit squares have
disjoint interiors — the separating-axis test for two rectangles).  200k random grid starts +
Nelder-Mead polish.  **One-sided and heuristic**: it exhibits configurations, it never proves their
absence; a comfortably negative maximum over a large search is evidence, not proof.  At
`t = 3.98`, `r = 1`:

| region | `k` | best min pairwise margin | reading |
|---|---|---|---|
| corner box `[0,1]^2` | 2 | — | `alpha = 1`, **proved and already checked by the verifier**: admissible part `[1/2,r]^2` has diameter `(2r-1)/sqrt2 < 1` |
| wall slot `[1,1.99] x [0,1]` | 2 | **+0.0568** (both at 15.98 deg, centres `(1, 0.618)` and `(1.99, 1.0)`) | `alpha >= 2` unconditionally |
| wall slot `[1,1.99] x [0,1]` | 3 | −0.4429 | `alpha = 2` |
| corner + 2 in the adjacent slot | 3 | **−0.2550** (all axis-parallel, centres 0.5 / 1.245 / 1.99, gaps 0.745) | with the corner occupied the slot holds one |
| wall strip `[1,2.98] x [0,1]` (no corners) | 3 | +0.0115 (at 80 deg — **tilting beats axis-parallel here**, the axis-parallel gaps are 0.99) | `alpha >= 3` |
| wall strip `[1,2.98] x [0,1]` | 4 | −0.3244 | `alpha = 3` |
| **full wall strip `[0,3.98] x [0,1]`** | 4 | **−0.0067** | `alpha = 3` — and `-0.0067 = (4 - 3.98)/3` exactly, the chord bound below |
| interior `[1,2.98]^2` | 5 | +0.3951 | |
| interior `[1,2.98]^2` | 6 | +0.0884 | `alpha(I) >= 6`; the interior gives no useful cap |

**The one that matters is a classical lemma.**  (**Proved and formalised, 2026-09-07**:
`notes/chord-lemma.md`, Lean `wall_strip_le_three_of_packing` in `lean/Sqpack/Chord.lean`.  The
sketch below is correct in outline; the write-up fixes the cut height to any
`a ∈ [(3-√2)/2, √2-1/2]` and the chord bound to `>= 1`, not `> 1`, which is why closed
disjointness — not interior disjointness — is the hypothesis at `t = 4`.)

> **Wall-strip capacity.**  In a packing of squares of side `L > 1` inside `[0,t]^2` with `t <= 4`,
> at most **3** squares have their centres within 1 of a given wall.

*Proof.*  A square of the packing whose centre is within 1 of the bottom wall straddles the line
`y = 0.9`: its top is at `cy + w/2 >= w >= L > 0.9` and its bottom at `cy - w/2 <= 1 - w/2 <= 0.5`.
Its intersection with that line has length `> 1` (Nagamochi 2005 Lemma 7(i); the minimum over
admissible poses is attained by the axis-parallel square standing on the wall, whose chord is
exactly `L`).  The squares are pairwise disjoint as closed sets — this is exactly what
`unit_subset_interior` gives for `L > 1` — so the chords are disjoint sub-intervals of a segment of
length `t <= 4`.  Four of them would need total length `> 4`.  ∎

Two remarks.  (i) The lemma is a **strictness** statement: for closed *unit* squares the minimum
chord is exactly 1 (computed: 1.000000 at `theta = 0`, `cy = 0.5`, for every cut height
`a in (0.5, sqrt2 - 1/2]`), and four unit squares in a row at `t = 4` touch.  It is therefore not
LP-certifiable in eroded semantics and has to be supplied analytically — `notes/proof-anatomy.md`
§7.4 lists precisely this statement as one of the "integral geometric lemmas about squares near the
container walls" that must come from mathematics rather than from the LP.  (ii) The margin the
numerical search reports at `t = 3.98`, `-0.0067`, is `(4 - 3.98)/3`: three gaps sharing the 0.02
of slack.  At `t = 4` the margin is 0 and the search would report "feasible" — the numerics cannot
see the lemma, which is the point of (i).

**Corollary (the leaf alphabet).**  In the corner leaf `k = 4` each wall strip already contains its
two corner squares, so it contains at most one further centre.  Hence each wall has **at most one
occupied slot**, and an occupied slot holds exactly one square:

```
K in {0,1}^8   with   K_0 + K_1 <= 1,  K_2 + K_3 <= 1,  K_4 + K_5 <= 1,  K_6 + K_7 <= 1 .
```

Also `m = sum K_j <= 4`, hence the interior holds `8 - m >= 4` squares (consistent with
`alpha(I) >= 6`).  Note that the conditional statement measured in the table ("corner + 2 in the
adjacent slot", margin `-0.255`) is now a *consequence* of the lemma, not an extra hypothesis; the
number is worth keeping only as an independent numerical check.

**Capacity as a column, not only as a branch.**  `alpha(strip) = 3` is also usable directly as a
valid inequality `mu(strip) <= 3` on the packing side, and dually as a *rank column* on the cover
side: crediting `v >= 0` to every pose centred in the strip costs `3v` in the budget.  This is the
brief's "capacity/rank column"; §6 shows `packing_le_weight_regions` already proves it.  On the
`k = 4` leaf dual it is **tight but not violated** (strip mass exactly 3), so it adds nothing there
by itself — its value is as the branching lemma, not as a cut.

### 2.4 Critical wall-box length, for completeness

For a wall box `[a, a+L] x [0, r]` that does **not** touch a side wall the two extreme centres are
`(a, r)` and `(a+L, w/2)` at a common angle, and the binding separating axis gives

```
alpha = 1   <=   max over phi in (0, pi/2) of  L cos phi + (r - w(phi)/2) sin phi  <  1 .
```

For `r = 1` the maximum is at `phi = 16.6095 deg` and the critical length is
`L* = 0.930804892...` (measured: `L = 0.93 -> -0.00077`, `0.87 -> -0.0582`, `0.99 -> +0.0568`).
`L*` is transcendental, so a verifier cannot check `L <= L*` in rationals; the crude exactly
checkable criterion is the one the verifier already applies to corner boxes ("the admissible part
has diameter `< 1`"):

```
4 L^2 + (2r - 1)^2 < 4          (integers, exact)     ->   L < sqrt3/2 = 0.866025 at r = 1.
```

### 2.5 Two designs

**Design A (recommended).**  Slots are the half-strips of §2.1.  The alphabet is `{0,1}` with the
per-wall constraint, by the chord lemma.  The regions partition the container, so a leaf also pins
the interior count at `8 - m` and nothing can hide.  **15 leaves.**  Price: the chord lemma has to
be supplied (cited, or proved in Lean — §6).

**Design B (fallback, no new mathematics).**  Slots shortened to `L = sqrt3/2` measured from the
corner box, `alpha = 1` each by the exact criterion of §2.4, no per-wall constraint used.  The
regions no longer partition — a gap of `t - 2r - 2L = 0.248` (0.268 at `t = 4`) sits at the middle
of each wall, and §1 says about 4 % of the wall mass is already sitting exactly there.  `2^8 = 256` patterns,
**43 leaves** up to symmetry.  Larger `r` shrinks the gap (`r = 1.1 -> 0.18`, `r = 1.2 -> 0.15`,
both still `(2r-1)^2 < 2`) at the cost of a weaker corner branch (`BRANCH.md`: `k = 4` at
`r = 1.2` gains only 0.1 over pure).

So the chord lemma is worth `43 -> 15` leaves, i.e. roughly a factor 3 in the whole cover-side
budget of §8.  It is the single highest-value piece of mathematics in this design.

---

## 3. The tree

```
                       packings of 12 unit squares of side L>1 in [0,t]^2
level 1: corner boxes C_0..C_3, per-box multipliers, K in {0,1}^4      6 leaves up to D4
    0000   1000   1100   1010   1110   1111
    11.24  11.56  11.83   -     11.97   >= 12.000024        (values at t = 3.98, BRANCH.md)
    cert.  cert.  cert.        rising    OPEN
                                            |
level 2 (under 1111 only): wall slots W_0..W_7, per-slot multipliers,
        K in {0,1}^8 with K_2i + K_2i+1 <= 1                          15 leaves up to D4
    m=0: 00000000
    m=1: 00000001
    m=2: 00000101  00000110  00001001  00010001  00010010
    m=3: 00010101  00010110  00011001  00100101
    m=4: 01010101  01010110  01011010  01100110
```

* **Multipliers.**  One `lam_j` per region, twelve regions (4 corner boxes + 8 slots); the leaf
  asserts `W - sum_j lam_j K_j < 12`.  Lean: `packing_le_weight_regions` with `m = 12`.  Unequal
  `lam_j` break the container symmetry, so the verifier sweeps `[0 deg, 90 deg)` and every point is
  its own column — which is already what the level-1 per-box form does, so this is not a new
  verifier mode.
* **Why not coarser.**  Shared multiplier over all eight slots: 9 D4-symmetric leaves, cheap, and
  the leaf `m = 4` contains the `k = 4` optimum unchanged (§1) — cuts nothing.  One multiplier per
  wall (alphabet `{0,1,2,3}` since `alpha(strip without corners) = 3`, or `{0,1}` with the chord
  lemma): the `k = 4` optimum has wall mass exactly 1 per strip, integral — cuts nothing either.
  Both are provably dead ends, which is the main negative result of this note.
* **Why not finer.**  Per-slot with the alphabet `{0,1,2}` (which is what is needed under a corner
  leaf with an *empty* corner box, where the chord lemma leaves room for two slot squares): 873
  classes, 512 with `m <= 8`.  Out of reach.
* **Under `k < 4`.**  The tree above hangs under `k = 4` only, because `k = 0,1,2` are certified
  closed at 3.98 and `k = 3` is a level-1 problem.  A push to `t = 4` needs a slot tree under every
  corner leaf; there the chord lemma gives `2 k_strip + (slots of that wall) <= 3` per wall rather
  than `<= 1`, and the leaf count grows accordingly (for `k = 0`: `{0,1,2,3}` per wall with at most
  3 per strip — 15 patterns per wall pair, hundreds of leaves).  This is a strong argument for
  finishing `t = 3.98` first and reaching `t = 4` by a different route (cliques).

---

## 4. Leaf values

### 4.1 Two LPs, and which way each errs

**Cover side** (`search/level2_cover_lp.py`).  `branch.py`'s LP is
`min W - sum_j lam_j K_j` subject to "every pose row captures `1 + sum_j lam_j [row in region j]`",
over a sampled set of pose rows and a point column set.  The `t398hk4` run saved both
(`runs/branch_t398hk4_cols.txt`, the point columns; `runs/branch_t398hk4_dual_it16.txt`, the pose
rows carrying dual mass).  This script rebuilds exactly that LP with twelve per-region multipliers
and re-solves it per pattern.  Semantics is the verifier's (rows are the concentric `sigma`-shrunk
squares with the saved half-side `h`), so the numbers are on the same footing as the published
`12.000024`.  *Error direction*: the pose rows are a sample, so the minimum is **too small** — an
optimistic leaf value.  A value `>= 12` is conclusive ("this leaf does not close"); a value `< 12`
is not.  The pattern-to-pattern **differences** on this fixed instance are the meaningful signal.

**Packing side** (`search/level2_lp.py`).  Over the whole container (no D4 reduction — the leaves
are asymmetric),

```
max sum mu_S  s.t.  coverage(p) <= 1 at every row point,  mu(C_i) = 1,  mu(W_j) = K_j,  mu >= 0
```

with `S` ranging over closed unit squares in `[0,t]^2`.  By LP duality this is the cover-side leaf
value, so `< 12` means the leaf closes.  Rows are generated from the exact arrangement-vertex
certification of `packing_dual.py` until the certified maximum coverage is `<= 1 + 1e-9`; columns
come from the certified pure supports, the `k = 4` branch dual, a coarse pose grid, and
reduced-cost pricing on a `0.03 / 2.5 deg` grid with coordinate-descent refinement.
*Error directions*: the pose set is finite, so the value is a **lower** bound on the leaf's true
value — again optimistic.  The row side is closed out (every reported value is post-certification).

**Why the packing-side absolute level is well below 12, and why that is expected.**  `branch.py`'s
`12.000024` is a cover LP over `sigma`-shrunk squares, whose dual is a fractional packing of shrunk
squares — a *relaxation* of the closed-unit-square packing that `level2_lp.py` solves.  That dual's
support has coverage `1.06` as a unit-square measure (`CLIQUE.md`), so its properly certified value
is about 11.3, not 12.  For reference the certified pure packing values are `L(3.98) = 11.918`,
`L(3.99) = 12.008`, `L(4.00) = 12.163` (`DUAL.md`), and `nu_f(3.98)` is somewhere in
`[11.918, 12.02]`.  **Read the gap column**, not the level.

### 4.2 Cover-side values on the `k = 4` instance (`t = 3.98`)

Instance: the `t398hk4` run's own point columns (6,155 orbits, subsampled by 4 to 1,539 orbits =
12,276 point variables) and its dual's pose rows (520 orbits = 4,156 rows in the full instance,
the heaviest 150 = 1,196 rows in the reduced one), twelve per-region multipliers, HiGHS IPM.
`runs/cover98e.txt` (full rows) and `runs/cover98d.txt` (reduced rows).

**Full-row instance** — the control reproduces the published `k = 4` value to six decimals, which is
the check that the instance is the right one:

| leaf | `m` | value `W - sum lam_j K_j` | `W` | cut vs the control |
|---|---|---|---|---|
| control (`k = 4` only) | – | **12.000000** | 12.00 | – (published: `12.000024`) |
| `01010101` | 4 | **11.444** | 19.33 | **0.556** |
| `01011010` | 4 | **11.300** | 14.80 | **0.700** |
| `01010110` | 4 | **11.286** | 14.57 | **0.714** |
| `01100110` | 4 | **11.000** | 19.00 | **1.000** |
| `00010101` | 3 | **10.848** | 14.76 | **1.152** |
| `00000101` | 2 | **10.000** | 13.00 | **2.000** |
| `00000001` | 1 | **9.000** | 9.00 | **3.000** |
| `00000000` | 0 | **8.000** | 8.00 | **4.000** |

**Reduced-row instance** (1,196 rows), which is what the earlier, cheaper run produced and which
shows the same shape one notch lower: control **11.9946**; `01011010` 11.235, `01010101` 11.200,
`01010110` 11.000, `01100110` 11.000, `00010101` 10.681, `00000101` 10.000, `00000001` 9.000,
`00000000` 8.000.

**The result to take away.**  Against a control that is exactly 12, every per-slot leaf drops, and
the *smallest* drop over the four hardest (`m = 4`) leaves is **0.556** — against a drop of
**0.000** for the wall-total and per-wall branches, which is not a measurement but the structural
observation of §1 that the `k = 4` optimum is already integral on those aggregates.  The `k = 4`
leaf is open by `2.4e-5`; a cut of half a unit is two orders of magnitude more than it needs.

How the LP uses the multipliers is worth reading off.  For an **empty** slot it drives `lam_j`
towards `-1`, i.e. it stops requiring anything of the poses centred there — that is the leaf's real
content, "no square of the packing is centred in slot `j`, so the cover need not handle those
poses".  For an **occupied** slot it pays `+lam_j` on those poses and buys it back once.  The four
corner multipliers come out at 0: once the empty slots are free the corner branch is not binding,
which is a second way of seeing that the level-2 constraint is doing the work.

**How much of the drop is real.**  These values are optimistic: the pose rows are the 520 orbits
that carried dual mass in the `k = 4` run, not the whole pose space, so an empty-slot `lam_j = -1`
is cheaper here than it would be against a verifier as separation oracle (a real cover would not
pay `W = 19.33` to buy back `4 x ~2`).  The point columns are subsampled by 4, which pushes the
other way but by very little — the control moves by `0.005` between the 12,276-column and the
49,104-column instances.  So read the table as **"the branch cuts, and the cut is large on the
`k = 4` leaf's own instance"**, not as "the leaves close at 11.0–11.4".  What settles it is a
cover-side run with the verifier as oracle (§8).

### 4.3 Packing-side values (`t = 3.98`)

Model: 2,331 poses (all eight images of the certified pure support `dual_PD1_support.txt` at 3.98,
plus a `0.25` / `15 deg` grid), rows generated from violated arrangement vertices until the certified
maximum coverage stops moving, HiGHS IPM, one core.  `runs/l2_h98.log`.

The **control** (corner masses fixed at 1 each, slots free) converges to `LP = 11.824056` at 15.7k
rows with certified `M = 1.005`, i.e. a certified `L = 11.76`.  That is 0.18 below `branch.py`'s
`12.000024` for the same leaf, for the reason in §4.1 (shrunk-square versus closed-unit-square
semantics), and it is the number the leaf values below have to be compared against — not 12.

| leaf | `m` | LP value | gap to the control | (cover-side gap, §4.2) |
|---|---|---|---|---|
| control (`k = 4` only) | – | **11.824056** | – | – |
| `01010101` | 4 | **11.335200** | **0.489** | 0.556 |
| `01011010` | 4 | **11.250000** | **0.574** | 0.700 |
| `01010110` | 4 | **11.200000** | **0.624** | 0.714 |
| `01100110` | 4 | **11.000000** | **0.824** | 1.000 |
| `00010101`, `00000101`, `00000001`, `00000000` | | *still running when this note was written* | | 1.152 / 2.000 / 3.000 / 4.000 |

**The two instruments agree.**  All four `m = 4` leaves come out in the same order on both sides,
with cuts of `0.49–0.82` (packing) against `0.56–1.00` (cover): the packing side is uniformly
12–20 % smaller, which is the direction its coarser pose set predicts.  They share no code beyond
`packing_dual.py`'s geometry kernels and err in opposite directions on the rows, so the size of the
cut is not an artefact of either instrument, even though neither establishes its exact value.

**Reading.**  The packing-side run is the weaker of the two instruments here: its absolute level is
set by how good the pose set is, and 2,331 poses is a small set (`packing_dual.py` needed column
generation over thousands of rounds to certify 11.918 for the *pure* problem).  Its value is as an
independent check that the leaf constraints bite in the same order as the cover side, on a
completely different LP with the opposite error direction on the rows.  For a full screening of all
15 leaves the run to make is the same script with `--warm-rounds 3` and a few pricing rounds — a
few hours on one core, which is the "screen before you commit cover-side compute" of §8.

### 4.4 `t = 4.0` (closed semantics)

Same script, same settings, warm-started from the certified closed-container support
`dual_PC1_support.txt` (86 orbits, mass 12.163) instead of `PD1`; 1,471 poses, 13.8k rows.
`runs/l2_c40.log`.

| leaf | `m` | LP value at `t = 4.0` | gap to the control | the same gap at `t = 3.98` |
|---|---|---|---|---|
| control (`k = 4` only) | – | **11.896104** | – | – (11.824056) |
| `01010101` | 4 | **11.385844** | **0.510** | 0.489 |
| `01100110` | 4 | **11.000000** | **0.896** | 0.824 |
| `00000101` | 2 | **10.000000** | **1.896** | *(cover side 2.000)* |
| `00000000` | 0 | **8.376812** | **3.519** | *(cover side 4.000)* |

**The cut is the same size at `t = 4` as at `t = 3.98`**, leaf for leaf (0.510 / 0.896 / 1.896
against 0.489 / 0.824 / 2.000 for the same three patterns on the same kind of pose set; the `m = 0`
leaf gives 3.519 against 4.000).  That is the useful fact: the level-2 constraint is not a `t`-dependent
trick, it cuts the same half unit at the container itself.

What it does *not* do is close anything at `t = 4`, and the arithmetic says why.  The absolute
levels above are set by the pose sets (§4.1) and are lower bounds; the honest references are the
certified pure values, `L(3.98) = 11.918` and `L(4.00) = 12.163`.  A cut of `0.5` applied to a leaf
that starts at `12.163 + (whatever the corner constraint adds)` does not obviously land below 12,
and unlike at `3.98` — where the leaf is open by `2.4e-5` and any positive cut suffices — at `t = 4`
the branch has to beat a `0.16+` deficit with a cut whose true (non-optimistic) size is unknown.
Add that at `t = 4` no corner leaf closes either, so the tree must hang under all six corner leaves
rather than one, and the leaf count and compute both multiply.

There is no cover-side `t = 4` analogue of §4.2 because there is no `t = 4` branch run to take an
instance from: `branch.py` has never been run at `t = 4` in the `k = 4` leaf (`BRANCH.md`: nothing
closes even at 3.99).  Producing one is the 13–27 h of §8, not a side computation.

Other reference points at `t = 4`, for the design:

| quantity | `t = 3.98` | `t = 4.00` |
|---|---|---|
| pure fractional packing value `L(t)` (certified, `DUAL.md`) | 11.918 | **12.163** |
| pure cover value `COVER(t)` (unconverged, rising) | ~12.02 | 12.3–12.5 (`CLOSED4.md`) |
| clique relaxation (restricted pose sets, drifting up) | ≲ 11.8 in the `k = 4` leaf | 11.8–12.0 |
| corner / wall / interior split of the pure measure | 3.44 / 4.80 / 3.68 | 3.20 / 5.40 / 3.56 |

The last row is what matters for this design: at `t = 4` the pure measure puts **more** mass in the
wall slots (5.40, i.e. 0.675 per slot) and less in the corners, so a per-slot branch has more
fractional mass to cut — which is consistent with the measured cut being no smaller there.

---

## 5. Which leaves close, which need cliques, which need `t = 4`

The 15 Design-A leaves fall into four kinds.  `m` is the number of occupied slots; the interior
then holds `8 - m` squares (the regions partition, so this is forced, not assumed).

| kind | leaves | cover-side value (§4.2) | what the packing must look like | measured capacity check |
|---|---|---|---|---|
| `m = 4` | `01010101`, `01010110`, `01011010`, `01100110` | **11.00 – 11.44** | 4 corner + 4 wall (one per wall) + 4 interior | 4 unit squares must fit in the interior `[1, t-1]^2` **plus** whatever room the gaps between the frame squares leave; `s(4) = 2 > 1.98`, so they cannot all sit inside `[1,2.98]^2` — they have to lean into the wall strips between the frame squares.  This is the family the `k = 4` optimum is an average of, and the hardest kind. |
| `m = 3` | `00010101`, `00010110`, `00011001`, `00100101` | **10.85** (`00010101`) | one wall empty, 5 interior | |
| `m = 2` | `00000101`, `00000110`, `00001001`, `00010001`, `00010010` | **10.00** (`00000101`) | two walls empty, 6 interior | `alpha(interior) >= 6` measured (margin `+0.0884`, six squares at 34–44 degrees leaning out of the interior into the wall gaps) — so these are right at the interior's capacity |
| `m <= 1` | `00000000`, `00000001` | **9.00 / 8.00** | 7 or 8 interior | `alpha(interior) = 6` in the searches here (7 squares miss by `-0.0229`, a limited search).  **If that bound is real these two leaves are refuted by capacity alone**, leaving 13 leaves.  Not established. |

Which leaves need what:

* **Cliques.**  §1: both the corner branch and the clique relaxation drive the *frame* to `4 + 4`
  and leave the residue in the interior.  Every level-2 leaf pins the frame exactly, so what is
  left in each leaf is an interior question — which is exactly what the clique constraints attack
  (`CLIQUE_CEILING.md`: the certified clique-feasible measure at 3.99 has `4 / 4 / 3.364`, the whole
  excess over 11 interior, and the non-Helly violation lives on wall-plus-interior clusters).  So
  the honest expectation is that the `m = 4` leaves need clique columns *in addition* to the branch,
  and the low-`m` leaves — which are capacity-tight in the interior — may close on the branch alone.
* **`t = 4`.**  **Not measured** — see §4.4.  What can be said without measuring: at `t = 4` the
  pure fractional packing value is `12.163` (certified) against `nu_f(3.98) in [11.918, 12.02]`, so
  every leaf starts 0.15–0.25 higher than the numbers above, and the tree has to hang under *every*
  corner leaf rather than only `k = 4`.  §9 for the verification side.
* **Zero-margin.**  Only a `t = 4` run needs it, and §9 says the tight poses are of the kinds
  `ZEROMARGIN.md` handles, with the region boundaries deliberately placed at integers so that the
  membership tests are not themselves zero-margin at the grid poses.

---

## 6. Lean

**The leaf statement needs no new theorem.**  `packing_le_weight_regions` takes
`R : Fin m -> R x R -> R -> Prop` (centre boxes and "contains `p`" are both instances),
`lam : Fin m -> R` of arbitrary sign, and concludes

```
(n : R) + sum_j lam_j * #{i | R j (ctr i) (ang i)} <= sum_a w a .
```

Substituting the leaf's counts, `W - sum_j lam_j K_j < n` refutes the leaf.  Nothing is special to
four regions or to counts in `{0,1}`; `m = 12` is the same theorem.

**Three things are new, in increasing order of work.**

1. *A rank/capacity corollary (about five lines).*  If `lam_j <= 0` and the count in region `j` is
   at most `alpha_j`, then `packing_le_weight_regions` gives `n <= W + sum_j alpha_j |lam_j|`.
   Box cliques are the case `alpha_j = 1` (`packing_le_weight_cliques`); `alpha_j > 1` is the rank
   column of §2.3 ("at most 3 centred in a wall strip").  A corollary of a theorem already in the
   file.
2. *The exact capacity criterion for Design B*, `4L^2 + (2r-1)^2 < 4` implies at most one centre in
   a wall box.  Same diameter argument as the corner boxes, so the natural home is the same one:
   an exact check in `verify/` and `xcheck.py`, as `(2r-1)^2 < 2` is today.
3. *The chord lemma of §2.3*, which is what makes Design A's 15 leaves exhaustive.  Two pieces: the
   chord bound (a one-variable trigonometric inequality — the minimum of the chord of a square of
   side `L` at angle `theta` with `cy in [w/2, r]` cut by `y = 0.9`, attained at `theta = 0` with
   value `L`) and the disjointness of the chords (immediate from `unit_subset_interior`, which is
   already in `Basic.lean`).  This is the one genuinely new development, and it is the same lemma
   the literature has used since 1984.

**What Lean still does not cover, at either level**: that the leaf set is *exhaustive*.  At level 1
this is the verifier's `(2r-1)^2 < 2`; at level 2 it is item 2 or item 3 above.  Exhaustiveness is
checked outside Lean today and this design does not change that.

---

## 7. Certificate format and verifier

The `region corner r_num r_den / lambda L1..L4 / k K1..K4` trailer generalises to a list of boxes:

```
regions Q g                     # coordinate denominator Q, number of regions g
L_j K_j b_j                     # region j: multiplier numerator over W, the leaf's count, #boxes
X0 X1 Y0 Y1                     # b_j box lines, integers over Q -- boxes of CENTRES
...
```

asserting: a closed unit square inside `[0,s]^2` whose centre lies in region `j` captures
`>= 1 + sum_j L_j / W` summed over the regions containing it, every other one `>= 1`, and
`sum_p w_p / W - sum_j (L_j / W) K_j < n`.  The existing corner trailer is the special case of four
boxes with a common `L` and `K = sum K_j`.

Verifier work, in order of size:

* **Generalise the region test from four hard-wired corner boxes to a list** — the cell-versus-box
  membership test is unchanged (bounding box of the cell against an axis-parallel box, exact for
  convex cells), only the bookkeeping is: a cell that may meet region `j` needs `1 + L_j`, a cell
  that may leave every region needs `1`, a straddling cell needs both.  This is the bulk of the
  change and it is mechanical.
* **Well-formedness checks**, by analogy with `(2r-1)^2 < 2` today: for Design B, each box with
  `K_j <= 1` must satisfy `4L^2 + (2r-1)^2 < 4`; for Design A, the verifier must additionally be
  told (or check) that the per-wall constraint holds — the chord lemma is not a numeric check, so
  the honest arrangement is a documented hypothesis in `FORMAT.md` plus a Lean proof.
* **Symmetry**: unequal `L_j` make the claim asymmetric, so the sweep is `[0 deg, 90 deg)` and every
  point is its own column.  Already the case for `lambda L1 L2 L3 L4`; no new mode.
* **`xcheck.py`** mirrors all of the above in exact rationals.

---

## 8. Cover-side cost per leaf

From `search/LPSPEED.md` (measured) and `search/BRANCH.md` (round counts):

| item | figure |
|---|---|
| LP per round, restricted master (`BRANCH_SOLVER=restricted`), ~3k-column master | 150–400 s at 58–68k rows; `~750 s` extrapolated at 186k rows |
| verifier per round, `N = 2000` | 100–190 s |
| pricing per round | 0.1–0.7 s (reduced cost) plus witness handling |
| **round** | **≈ 20 min measured** (`runs/restr1110b.out`, the asymmetric `1110` leaf); ~15 min at full size |
| rounds to convergence, asymmetric leaf | `1110` was at round 51 and unconverged; a 12-region leaf has three times the multipliers — estimate **40–80** |
| **per leaf** | **13–27 h** |
| **Design A, 15 leaves** | **200–400 core-hours** = 4–8 days on the two cores, ~1–2 days on eight |
| **Design B, 43 leaves** | **560–1160 core-hours** = 12–24 days on two cores |
| exact re-check (`xcheck.py`, `N = 2000`, ~10k atoms) | 8 h+ per leaf: another 120 h (A) or 350 h (B) |

Two consequences for how to spend it.

* **Screen on the packing side first.**  All 15 leaf values come out of one shared `level2_lp.py`
  model in a few hours; only the leaves whose optimistic value is near 12 need a cover-side run.
* **Do not start the tree until the `k = 4` leaf is closed by *something*.**  Even Design A is a
  4–8 day commitment on this machine whose payoff at `t = 3.98` is a bound within 0.005 of the pure
  ceiling.  The clique columns in the same leaf (`TODO.md` Phase 2) are a much smaller bet for the
  same `s(12) >= 3.98`, and §0.6 says they are also the part that scales to `t = 4`.

---

## 9. `t = 4` and zero-margin verification

Which leaves would have to be closed *at* `t = 4` (closed semantics, verified by the exact checker
of `search/ZEROMARGIN.md`) and whether their tight poses are of the kinds it handles:

* **All of them.**  There is no `t < 4` shortcut inside this tree: the tree exists to close the
  corner leaf `k = 4`, and at `t = 4` no corner leaf closes on its own (`BRANCH.md`: nothing closes
  even at 3.99).  So a `t = 4` proof runs the level-2 tree under every corner leaf, at margin zero.
* **The region boundaries are in the right place.**  With `r = 1` and `t = 4` the region boundaries
  in the centre plane are at `x, y in {0, 1, 2, 3, 4}` and the sixteen grid poses — the tight poses
  of any cover at `t = 4` — sit at `(0.5 + i, 0.5 + j)`, i.e. `0.5` from every boundary.  So the
  region membership tests are *not* zero-margin at the tight poses, and a cell straddling a region
  boundary there is not forced to satisfy two thresholds at once.  This is a design property worth
  keeping: **put region boundaries at integers.**  Design B's shortened slots have boundaries at
  `1 + sqrt3/2 = 1.866` and `3 - sqrt3/2 = 2.134`, still `0.366` from the grid centres — acceptable.
* **The tight poses are of the kinds the checker handles** (`ZEROMARGIN.md` §4): corner squares by
  the 2x2-box lemma **P1**; axis-parallel wall squares by P1 with the wall points plus the exact bin
  **CORE** beyond 11.5 degrees; axis-parallel interior grid poses by CORE, provided the root box
  lattice is commensurable with halves (pitch `1/10` is); tilted one-parameter families from
  unit-distance point pairs by the triangle lemma **TRI**.
* **The known gap is unchanged by this design.**  `ZEROMARGIN.md` §5: with a *weighted* cover the
  triangle lemma only delivers `min_v w_v`, so a weighted optimum that is tight on a tilted family
  would need the two-region certification that is not built.  Whether the optimum has such a family
  is an empirical question about the LP; the level-2 branch does not make it worse (it adds
  axis-parallel constraints, which push the optimum towards the grid, i.e. towards the families the
  checker *does* handle).
* **The chord lemma is itself a zero-margin statement at `t = 4`** (§2.3, remark (i)): four
  axis-parallel unit squares in a row across `[0,4]` touch.  It is not verified numerically at any
  `t`; it is proved once, for side `L > 1`, and used as a hypothesis.  That is the same status it
  has in Nagamochi and in Bentz.

---

## 10. What to do

1. ~~Prove (or formalise) the wall-strip chord lemma.~~  **Done** (2026-09-07):
   `notes/chord-lemma.md`, Lean `wall_strip_le_three_of_packing`.  `notes/branch-semantics.md` §4
   redoes the leaf enumeration with it, for every corner leaf and not only `k = 4`.
2. Do not build a symmetric or per-wall level-2 branch: §1 shows the `k = 4` optimum is integral on
   both aggregates.
3. Screen all 15 Design-A leaves with `search/level2_lp.py` on one shared model before committing
   any cover-side compute.
4. Generalise the `region` trailer to a list of boxes in `verify/`, `xcheck.py` and
   `certificates/FORMAT.md`; the cell membership test itself is unchanged.
5. Run the clique columns in the `k = 4` leaf first (`TODO.md` Phase 2).  If they close it, this
   tree is not needed for `s(12) >= 3.98`; if they do not, this tree is the fallback, and at `t = 4`
   the two have to be combined anyway.


---

## 11. Reproduce

```sh
R=/home/evand/math/square-packing/s12/runs        # or this worktree's runs/

# where the mass sits (table of section 1)
python3 search/level2_regions.py $R/branch_t398hk4_dual_it16.txt
python3 search/level2_regions.py $R/dual_PD1_support.txt
python3 search/level2_regions.py $R/cc_C99c_exact.txt

# capacities (section 2.3); one-sided and heuristic, minutes each
python3 search/level2_capacity.py --t 3.98 --box 1 1.99 0 1 --k 2          # a wall slot: alpha >= 2
python3 search/level2_capacity.py --t 3.98 --box 0 3.98 0 1 --k 4          # a wall strip: alpha = 3
python3 search/level2_capacity.py --t 3.98 --box 0 1 0 1 --box 1 1.99 0 1 --k 3   # corner + 2 in a slot
python3 search/level2_capacity.py --t 3.98 --box 1 2.98 1 2.98 --k 6       # the interior

# leaf counts (section 3)
python3 search/level2_leaves.py

# cover-side leaf values on the k = 4 instance (section 4.2), ~2 min per pattern
python3 search/level2_cover_lp.py --cols $R/branch_t398hk4_cols.txt \
    --rows $R/branch_t398hk4_dual_it16.txt --method highs-ipm --col-stride 4 \
    --patterns "........,01010101,01010110,01011010,01100110,00010101,00000101,00000001,00000000"

# packing-side leaf values (section 4.3)
python3 search/level2_lp.py 3.98 TAG --corners 1111 --warm-rounds 0 --rounds 1 --rowloops 12 \
    --seed-pitch 0.25 --seed-dth 15 --row-pitch 0.08 --method highs-ipm --prune-rows 0 \
    --warm $R/dual_PD1_support.txt --threads 1 \
    --patterns 01010101,01010110,01011010,01100110,00010101,00000101,00000001,00000000
```

`level2_lp.py` imports `packing_dual.py` from the same `search/` directory (it falls back to the
main tree's path when run from a worktree that does not have it) and compiles its C kernels into
`runs/`.  Do **not** prune rows in `level2_lp.py` (`--prune-rows 0`): pruning by dual value alone
makes the cutting-plane loop cycle, exactly as `BRANCH.md` records for the cover side.
