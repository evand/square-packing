# Anchor cliques: a certifiable clique family for rotated unit squares (task G, part 1)

Code: `search/clique_family.py` (geometry, search, separation), `search/clique_cert.py` (exact
rational certification of the clique property), `search/clique_continuum.py` (the LP that uses
them).  Results and the go/no-go verdict: `search/CLIQUE_CONTINUUM.md`.

## 0. Semantics

A certificate at container side `t` refutes packings at side `t' < t`; rescaled to `[0,t]²`
those are pairwise **disjoint closed** unit squares.  A *clique* is therefore a set of poses
that pairwise **closed-intersect** (touching counts), and every fractional packing measure `mu`
satisfies `mu(K) <= 1` for every clique `K`.  Throughout:

* `C = [0,t]²`, `Q = [-½,½]²`, pose `(c, θ)`, `S(c,θ) = c + R_θ Q` closed;
* `S` is **admissible** iff `S ⊆ C`, i.e. `c ∈ [w/2, t - w/2]²` with `w(θ) = |cos θ| + |sin θ|`;
* `P_p = { S admissible : p ∈ S }` — the **point clique** of `p`.  Its constraint
  `mu(P_p) <= 1` is exactly the coverage constraint the point-cover LP already imposes, so a
  clique is only *new* to the extent that it is not contained in a point clique.
* `U(p) = P_p` when we think of it as a family to be pierced.
* `K(p) = { S admissible : S ∩ S' ≠ ∅ for every S' ∈ P_p }` — the **maximal point-anchored
  family**: the largest set of poses that could possibly be added to the point clique of `p`.
  (`P_p ⊆ K(p)`, since two squares through `p` meet at `p`.)

## 1. The family

**Definition (anchor clique).**  For a point `p ∈ C` and a nonempty compact convex `A ⊆ ℝ²`,

    K(p, A)  =  { S admissible : p ∈ S and S ∩ A ≠ ∅ }  ∪  { S admissible : A ⊆ S }.

**Lemma 0.**  `K(p, A)` is a clique, for *every* `p` and *every* nonempty `A`.

*Proof.*  Two members of the first part both contain `p`.  Two members of the second part both
contain `A ≠ ∅`.  For a mixed pair `S` (first part), `S'` (second part): `S ∩ S' ⊇ S ∩ A ≠ ∅`
because `A ⊆ S'` and `S` meets `A`. ∎

That is the whole validity argument — no geometry, no non-avoidance lemma.  Three remarks:

* `A = {p}` gives back the point clique, so the family contains the point constraints.
* The general `m`-anchor version is `K(A_1,…,A_m) = ⋃_i { S : A_i ⊆ S and S ∩ A_j ≠ ∅ for
  every `j` in a chosen "direction set" }`; for each unordered pair `{i,j}` one has to choose
  which side carries the "meets" filter.  Task B's **box cliques** are the special case in which
  all the anchors pairwise *meet* (`A_i ∩ A_j ≠ ∅`), so no filter is needed and the family is
  `⋃_i {S : A_i ⊆ S}`; the new content here is anchors that do **not** pairwise meet, which is
  where the non-Helly mass lives.
* The filter `S ∩ A ≠ ∅` is what makes the family valid without a lemma; the price is that
  `K(p, A)` need not contain all of `P_p`.  §3 says exactly when it does.

**What is new relative to `search/BOXCLIQUE.md`.**  A box clique can never contain a point
clique: the core of "all poses through `p`" is `{p}` and the verifier's cores are inward-rounded
rectangles.  `K(p, A)` contains `P_p` exactly (when Lemma 2 applies), which is why it can be a
*strengthening* of a point row rather than an extra row that the LP ignores.

## 2. Where cliques can help at all: the wall

**Lemma 1 (no wall, no gain).**  If `1 ≤ p_x ≤ t-1` and `1 ≤ p_y ≤ t-1` then `K(p) = P_p`.

*Proof.*  Let `S'` be admissible with `p ∉ S'`.  `S'` is a closed square, so one of its four
outward edge normals `n` has `⟨p,n⟩ > h_{S'}(n)`, i.e. `S' ⊆ {x : ⟨x-p,n⟩ ≤ -d}` for some
`d > 0`.  Put `σ_i = sign(n_i)` (either sign if `n_i = 0`) and let `S` be the axis-parallel unit
square with corner `p` in that quadrant, `S = p + σ ⊙ [0,1]²`, centre `c = p + σ/2`.  For
`x ∈ S`, `⟨x-p, n⟩ = Σ_i |n_i| σ_i (x_i - p_i) ≥ 0`, so `S ∩ S' = ∅`.  `S` is admissible because
`w = 1` and `c ∈ [½, t-½]²` — for `σ_1 = +1` this needs `p_x ≤ t-1`, for `σ_1 = -1` it needs
`p_x ≥ 1`, and likewise in `y`.  And `p ∈ S`.  So `S' ∉ K(p)`. ∎

**Sharp.**  Measured (`clique_family.py kset`, exact-arithmetic-free but a rigorous *outer*
approximation of `K(p)`, so an over-estimate of the extra): the extra volume is `0` at wall
distance `1.000` and positive below it.  Volumes in pose space `(c_x, c_y, θ)`, `t = 3.99`:

| `p` | vol `P_p` | vol `K(p)` | vol `K(p) \ P_p` | ratio |
|---|---|---|---|---|
| `(1.001, 1.999)` | 1.38766 | 1.38766 | 0 | 0 |
| `(1.000, 1.999)` | 1.38628 | 1.38628 | 0 | 0 |
| `(0.999, 1.999)` | 1.38489 | 1.38495 | 5.70e-5 | 4.1e-5 |
| `(0.99, 1.995)` | 1.37215 | 1.37391 | 1.765e-3 | 1.29e-3 |
| `(0.95, 1.995)` | 1.31224 | 1.33038 | 1.813e-2 | 1.38e-2 |
| `(0.90, 1.995)` | 1.23236 | 1.28067 | 4.831e-2 | 3.92e-2 |
| `(0.80, 1.995)` | 1.06283 | 1.19026 | 1.274e-1 | 0.120 |
| `(0.60, 1.995)` | 0.71219 | 1.06167 | 3.495e-1 | 0.491 |
| `(0.50, 1.995)` | 0.54278 | 1.04327 | 5.005e-1 | 0.922 |
| `(0.99, 0.99)` (corner) | 1.18167 | 1.20317 | 2.149e-2 | 1.82e-2 |
| `(0.95, 0.95)` | 1.07988 | 1.19416 | 1.143e-1 | 0.106 |
| `(0.90, 0.90)` | 0.95470 | 1.19246 | 2.378e-1 | 0.249 |
| `(1.2, 1.995)`, `(1.5, 1.995)`, `(1.995, 1.995)`, `(1,1)`, `(1.05,1.05)` | — | — | 0 | 0 |

**This is the structural fact of the whole task.**  A clique can beat the point constraint at
`p` only if `p` is *strictly* within distance 1 of a wall, and the gain is `O(1 - dist)`.  The
tight points of the certified extremal measure at `t = 3.99` sit at `p_x ≈ 0.99` (the wall
squares reach `x = 1`), where the whole extra family is `0.13 %` of the point clique.  Corners
help by a factor ≈ 12 at the same wall distance but are still small.

Numerically (15,930 sampled pairs over 9 points, `clique_family.py kclique`) `K(p)` is itself
pairwise closed-intersecting, i.e. it is *the* maximal clique containing `P_p`.  Not proved.

## 3. When the anchor clique contains the whole point clique

`K(p, A) ⊇ P_p` iff every admissible `S ∋ p` meets `A`; call such an `A` a **transversal** of
`P_p`.  Then `mu(K(p,A)) ≤ 1` strictly dominates the coverage constraint at `p`.

**The complete test.**  `S(c,θ)` meets `A` iff `c ∈ A ⊕ R_θ Q` (`Q` is centrally symmetric), a
convex set; and `{ c : S(c,θ) admissible, p ∈ S(c,θ) }` is a convex polygon (at most 8 sides:
the container box and the rotated square `p + R_θ Q`).  So

> `A` is a transversal  ⟺  for every `θ`, every **vertex** of that polygon `v` has
> `S(v,θ) ∩ A ≠ ∅`.

`clique_family.transversal_margin` evaluates this vectorised over `θ` (and returns the SAT
margin, which is what the search and the certifier use).

**Lemma 2 (the transversal threshold at a wall).**  Let `0 < p_x < 1`, `0 < ε < 1 - p_x`, and
let `A_ρ = {p_x + ε} × [p_y - ρ, p_y + ρ]` be the vertical segment at horizontal offset `ε`.
Put `θ₀ = arccos p_x` and

    ρ*  =  ε · cot θ₀  =  ε · p_x / sqrt(1 - p_x²).

Then every admissible closed unit square containing `p` meets `A_ρ` **iff** `ρ ≥ ρ*`
(the "only if" needs `w(θ₀) = p_x + sqrt(1-p_x²) ≤ p_y ≤ t - w(θ₀)`; the "if" needs nothing
about `p_y`).

*Proof.*  Squares are `90°`-periodic in `θ`, so take `θ ∈ [0°, 90°)`; then `cos θ, sin θ ≥ 0`,
`e_1 = (cos θ, sin θ)`, `e_2 = (-sin θ, cos θ)`, `w = cos θ + sin θ`.  Write
`s_i = ⟨p - c, e_i⟩`, so `p ∈ S` ⟺ `|s_1|, |s_2| ≤ ½`, and `c = p - s_1 e_1 - s_2 e_2`.

A convex polygon and a segment are disjoint iff one of three axes separates them: the segment's
normal `e_x`, and the square's two edge normals.

*Axis `e_x`.*  Separation needs `max_S x < p_x + ε` or `min_S x > p_x + ε`.  The second fails
because `p ∈ S` gives `min_S x ≤ p_x < p_x + ε`.  The first fails because admissibility gives
`max_S x = c_x + w/2 ≥ w/2 + w/2 = w ≥ 1 > p_x + ε`.  So this axis never separates.

*Axes `e_1, e_2`.*  With `m = (p_x + ε, p_y)` the midpoint of `A_ρ`, "the whole segment lies
beyond the edge `(i, σ)`" reads `σ⟨m - c, e_i⟩ - ρ |⟨e_i, e_y⟩| > ½`, and
`⟨m - c, e_i⟩ = ε⟨e_x, e_i⟩ + s_i`.  The four cases:

    (1,+):  ε cos θ + s_1 - ½  ≤  ρ sin θ   ?
    (1,-):  -ε cos θ - s_1 - ½ ≤ -ε cos θ ≤ 0 ≤ ρ sin θ           (since s_1 ≥ -½)   ✓
    (2,+):  -ε sin θ + s_2 - ½ ≤ 0 ≤ ρ cos θ                      (since s_2 ≤ ½)    ✓
    (2,-):  ε sin θ - s_2 - ½  ≤  ρ cos θ   ?

Admissibility in `x` bounds `s_1` and `s_2`:
`c_x = p_x - s_1 cos θ + s_2 sin θ ≥ w/2` gives, using `s_2 ≤ ½`,
`s_1 ≤ p_x / cos θ - ½`, and using `s_1 ≥ -½`, `s_2 ≥ ½ - p_x / sin θ`.

**(1,+).**  Need `s_1 - ½ + ε cos θ ≤ ρ sin θ` with `s_1 ≤ min(½, p_x/cos θ - ½)`.
If `cos θ ≤ p_x` (i.e. `θ ≥ θ₀`) the min is `½` and the requirement is `ρ ≥ ε cot θ`, implied by
`ρ ≥ ε cot θ₀` since `cot` decreases.  If `cos θ > p_x` (`θ < θ₀`) the requirement is
`f(θ) := p_x/cos θ + ε cos θ - ρ sin θ - 1 ≤ 0`.  Now `f(0) = p_x + ε - 1 < 0` by hypothesis, and
`f(θ₀) = ε p_x - ρ sqrt(1-p_x²) ≤ 0` exactly when `ρ ≥ ρ*`.  Finally
`f'(θ) = sin θ (p_x/cos²θ - ε) - ρ cos θ` is increasing on `[0°, 90°)` (both terms are), so `f`
is **convex** and hence `f ≤ max(f(0), f(θ₀)) ≤ 0` on `[0, θ₀]`.

**(2,-).**  Need `ε sin θ - s_2 - ½ ≤ ρ cos θ`.  If `sin θ < p_x` use `s_2 ≥ -½`: the requirement
is `ρ ≥ ε tan θ`, and `tan θ < tan(90° - θ₀) = cot θ₀`, so `ρ ≥ ρ*` suffices.  If `sin θ ≥ p_x`
use `s_2 ≥ ½ - p_x/sin θ`: the requirement becomes
`p_x/sin θ + ε sin θ - ρ cos θ - 1 ≤ 0`, which is `f(90° - θ) ≤ 0` with `90° - θ ≤ θ₀` — the case
just proved.

*Sharpness.*  Take `θ = θ₀` and `s_1 = s_2 = ½`, i.e. `p` at the corner `c + ½(e_1 + e_2)`.  Then
`c_x = p_x - ½cos θ₀ + ½ sin θ₀ = w(θ₀)/2` exactly (using `cos θ₀ = p_x`) and
`c_y = p_y - w(θ₀)/2`, both admissible when `w(θ₀) ≤ p_y ≤ t`.  The `(1,+)` quantity is
`ε cos θ₀ + ½ - ½ = ε p_x`, so the segment is strictly beyond that edge — the square misses
`A_ρ` — exactly when `ρ sqrt(1-p_x²) < ε p_x`, i.e. `ρ < ρ*`. ∎

**Stress test** (`clique_cert.py stress`, floats, independent Liang–Barsky clipping, no shared
code with the proof): 300,000 random `(p_x ∈ [0.30,0.999], p_y, ε, ρ = (1+s)ρ*, admissible pose
through p)`: **0** failures of "`S` meets `A`"; and 300,000 checks that the sharpness witness at
`ρ = 0.999 ρ*` really is admissible, contains `p` and misses `A`: **0** failures.

**Consequences.**

* `ε` is capped twice: by `ε < 1 - p_x` (axis `e_x`) and by `ρ = ρ* < ½` (a unit square must be
  able to contain a vertical segment of length `2ρ`), i.e. `ε < sqrt(1-p_x²)/(2 p_x)`.  For
  `p_x = 0.99` the first is binding: `ε < 0.01`, `ρ < 0.0702`.
* The extra family `{S : A_ρ ⊆ S} \ P_p` — squares that contain the anchor but *not* `p` — has
  positive measure: at `θ = 0` it is `c_x ∈ (p_x + ½, p_x + ε + ½]` (width `ε`) times
  `c_y ∈ [p_y + ρ - ½, p_y - ρ + ½]` (width `1 - 2ρ`).  Measured volumes, with `ε` optimised and
  `ρ = 1.05 ρ*` (`clique_family.py anchorvol`):

| `p` | `ε` | `ρ` | vol `Q_A` | vol `Q_A \ P_p` | vol `K(p) \ P_p` | fraction of the maximum |
|---|---|---|---|---|---|---|
| `(0.99, 1.995)` | 0.0098 | 0.0722 | 1.1644 | 1.132e-3 | 1.823e-3 | 0.62 |
| `(0.95, 1.995)` | 0.0490 | 0.1565 | 0.9181 | 1.032e-2 | 1.842e-2 | 0.56 |
| `(0.90, 1.995)` | 0.0980 | 0.2125 | 0.7625 | 2.503e-2 | 4.881e-2 | 0.51 |
| `(0.80, 1.995)` | 0.1817 | 0.2544 | 0.6403 | 5.911e-2 | 1.281e-1 | 0.46 |
| `(0.60, 1.995)` | 0.3348 | 0.2636 | 0.5837 | 1.650e-1 | 3.487e-1 | 0.47 |
| `(0.50, 1.995)` | 0.4542 | 0.2754 | 0.5691 | 2.499e-1 | 5.008e-1 | 0.50 |
| `(0.99, 0.99)` | 0.0098 | 0.0722 | 1.0506 | 1.128e-3 | 2.155e-2 | 0.05 |
| `(0.90, 0.90)` | 0.0980 | 0.2125 | 0.6876 | 2.412e-2 | 2.384e-1 | 0.10 |

  So one vertical segment recovers half to two thirds of the theoretical maximum along a wall,
  and only 5–10 % of it at a corner (where the right anchor is not a wall-perpendicular
  segment; a two-anchor family would do better — not built).

## 4. Exact, terminating certification

`search/clique_cert.py` proves statement **(T)** — *every admissible closed unit square
containing `p` meets `A`* — in exact rational arithmetic, by adaptive subdivision of pose space
`(c_x, c_y, u)`, `u = tan(θ/2)`, exactly as `search/zeromargin.py` does.  Conventions and
primitives:

* `cos = (1-u²)/(1+u²)`, `sin = 2u/(1+u²)` — rational at box corners; `cos` decreasing and `sin`
  increasing on `u ∈ [0,1]` (`θ ∈ [0°, 90°]`), so the interval enclosures are exact endpoints.
* `w_lo` over a bin is `min(w(u_0), w(u_1))` (`w` is unimodal on `[0°,90°]`), which is `≥ 1`; the
  centre box is **clipped** to `[w_lo/2, t - w_lo/2]²`, which keeps every admissible pose of the
  box.  That clip is also what makes the `e_x`-axis test automatic (`c_x + w/2 ≥ w_lo ≥ 1 > x_0`).
* **EMPTY** — the clipped box is empty, or the interval enclosure of `⟨p-c, e_i⟩` lies outside
  `[-½, ½]` for some `i`, so no pose of the box has `p ∈ S`.
* **MEET** — all six (axis, side) separations are refuted over the box by interval arithmetic:
  the anchor's normal `d` (two sides) and, for each square normal `e_i`, "some endpoint of `A` is
  on this side of the edge" (two sides).
* Otherwise split the longest scaled side (`u` scaled by 2).

**Termination.**  With `ρ > ρ*` strictly the statement holds with a uniform positive margin, and
because the over-tested inadmissible sliver of a bin has width `(w_hi - w_lo)/2 → 0`, the
subdivision terminates at bin width `O(margin)`.  Run, `t = 399/100`, `p = (99/100, 2)`,
`ε = 1/250`, `ρ = 1.05 ρ*` (a rational `230557/7822046`):

    201,068 boxes, depth 30, EMPTY 60,412, MEET 46,522, uncertified 0     (44 s)
    VERDICT: CERTIFIED

**Where subdivision does not terminate, and the primitive that fixes it.**  At `ρ = ρ*`
exactly — the *maximal* transversal family — the pose `θ = θ₀ = arccos p_x`, `c_x = w(θ₀)/2`,
`c_y = p_y - w(θ₀)/2` touches `A` in a single point, and every box around it contains poses with
`c_x` slightly below `w(θ₀)/2` — inadmissible, but not excluded by any box test, because the
admissibility boundary `c_x = w(θ)/2` is a *curve* and a box always straddles it.  This is
exactly the situation of `ZEROMARGIN.md §4.1` (the corner square: quadratic margin against linear
erosion), and no erosion or core argument reaches it.  **Lemma 2 is itself the primitive**: it
decides the maximal case in closed form, in one test, for the whole one-parameter family.  In
practice one never wants `ρ = ρ*`: taking `ρ = (1+s) ρ*` with `s ≈ 0.05` costs `5 %` of the
anchor length and buys a subdivision proof, so the certifier is the operational tool and Lemma 2
is the statement one would formalise.

The other margin-zero place is **membership**, not validity: a pose whose boundary passes exactly
through an endpoint of `A` is on the boundary of `{S : A ⊆ S}`.  For a *cover* certificate this
is harmless — under-crediting a clique is sound (it only makes the check harder to pass), and the
verifier simply does not credit a swept cell that is not entirely inside the region.

## 5. What it would cost in `verify/`, `xcheck.py` and Lean

The format and the reduction already exist (`certificates/FORMAT.md` "Clique certificates",
`lean/Sqpack/Basic.lean` `packing_le_weight_cliques`), so this is a change of *clique object*,
not of the certificate semantics.

**Certificate format.**  The present block is

```
cliques N Q c
w_1 b_1
k U0LO U0HI U1LO U1HI      (b_1 box lines, rectangles in the frame of bin k)
```

The anchor version keeps `cliques N Q c` and the weight line and replaces the box lines by
*piece* lines, each a (contains-anchor, meets-anchor-list) pair:

    piece  a  f_1 … f_r        this piece is { S : A_a ⊆ S and S ∩ A_{f_i} ≠ ∅ for all i }
    anchorP  X Y D             anchor: the point (X/D, Y/D)
    anchorS  X0 Y0 X1 Y1 D     anchor: the segment with those rational endpoints

A pose belongs to the clique iff it belongs to one of the pieces.  Well-formedness the verifier
must check (`ERROR` otherwise): for every ordered pair of pieces `(i, j)`, either the anchors
`A_i, A_j` intersect (exact rational test) or `j` is in `i`'s filter list or `i` is in `j`'s —
that is precisely Lemma 0's hypothesis, and it is the analogue of the current "cores pairwise
meet" refusal.  Box cliques stay as they are (they are the case where all anchors pairwise
intersect and no filters are needed).  Weight and the counted-once rule are unchanged.

**Verifier (`verify/src/main.rs`).**  The sweep is per angle bin `k` and per cell of centres.
Two exact `i128` predicates are needed, both over the bin's rational `cos, sin` (the verifier
already carries `σ_k` and rational rotations):

* `contains(anchor, cell, k)` — every pose of the cell contains the anchor.  For a point this is
  the existing "point in the shrunk square" test with the **exact bin core** rather than the
  `σ`-square (`ZEROMARGIN.md §2`), which is what makes the point clique exact at edge midpoints;
  for a segment it is the same test on both endpoints (the square is convex).
* `meets(anchor, cell, k)` — every pose of the cell meets the anchor: refute the six
  (axis, side) separations by interval arithmetic in the bin, exactly as `clique_cert.py` does.
  Conservative failure (no credit) is sound.

A cell that only partially satisfies a predicate gets **no** credit, and — as with box cliques —
its LP witness must be placed outside the region so the cutting-plane loop cannot cycle.  Cliques
still have no representable images under the container symmetries, so a certificate carrying them
is swept over the full `[0°, 90°]`.

**`xcheck.py`.**  The same two predicates in `fractions.Fraction`, with anchors *inflated* rather
than eroded (so `xcheck` accepts everything the Rust accepts), plus the credit rule in the cell
sweep.  Roughly the same amount of code as the box-clique block it already has.

**Lean.**  `packing_le_weight_cliques` already takes an arbitrary pose predicate
`K : Fin m → ℝ × ℝ → ℝ → Prop` with the hypothesis

```lean
(hK : ∀ j c c' θ θ', K j c θ → K j c' θ' → (sq c θ 1 ∩ sq c' θ' 1).Nonempty)
```

so nothing in the reduction changes; what is needed is a sibling of `clique_of_cores` that
discharges `hK` for the anchor family.  It is Lemma 0 — no geometry, no shrink lemma:

```lean
/-- **Anchor cliques are cliques.**  Pieces indexed by `ι`; piece `i` is the set of poses that
contain `anc i` and meet `anc j` for every `j` in `filt i`; the well-formedness hypothesis
`hpair` says every ordered pair is covered either by intersecting anchors or by a filter. -/
lemma clique_of_anchors {ι : Type*} (P : ι → ℝ × ℝ → ℝ → Prop) (anc : ι → Set (ℝ × ℝ))
    (hin  : ∀ i c θ, P i c θ → anc i ⊆ sq c θ 1)
    (hmeet: ∀ i j c θ, P i c θ → Meets i j → (sq c θ 1 ∩ anc j).Nonempty)
    (hpair: ∀ i j, (anc i ∩ anc j).Nonempty ∨ Meets i j ∨ Meets j i) :
    ∀ c c' θ θ', (∃ i, P i c θ) → (∃ j, P j c' θ') →
      (sq c θ 1 ∩ sq c' θ' 1).Nonempty
```

with a three-way case split on `hpair i j`: intersecting anchors give
`(anc i ∩ anc j) ⊆ sq c θ 1 ∩ sq c' θ' 1` (exactly `clique_of_cores`'s argument);
`Meets i j` gives `sq c θ 1 ∩ anc j ⊆ sq c θ 1 ∩ sq c' θ' 1` since `anc j ⊆ sq c' θ' 1`;
`Meets j i` symmetrically.  `clique_of_cores` stays for box cliques.  The verifier's two
predicates (§ above) are what discharge `hin` and `hmeet` for a swept cell, exactly as `hcore` is
discharged today.

Lemma 2 (`K(p,A) ⊇ P_p`) is *not* needed for soundness and would not be formalised: it says the
clique constraint dominates the point row, which is a statement about the strength of the
certificate, not its validity.  Formalising it is real one-parameter trigonometry (the convexity
of `f(θ) = p_x/cos θ + ε cos θ - ρ sin θ - 1` on `[0, θ₀]`), and the exact certifier of §4 is the
practical substitute.

## 6. The measured facts (what the family is worth)

On the exact certified extremal measure at `t = 3.99` (`runs/dual_exact_3.99_support.txt`, mass
`12.008230754`, coverage `≤ 1` at every arrangement vertex):

| object | max mass | note |
|---|---|---|
| point clique (= coverage) | `1.000000` | certified exactly by `dual_exact.py` |
| max-mass clique of the support's closed-intersection graph | **1.2139** (B&B, 60–120 s cap; task A reports 1.327 with a longer run) | valid, but its members are individual poses |
| best **anchor clique** `K(p,A)` found (`clique_family.py separate`, 60 anchor points × 5 seeds × 20 greedy steps) | **1.0502** | `p = (3.01, 2.61)`, `A` = one whole support square: `{S : A ⊆ S} = {A}` — a *single pose*, so the cut is fragile under column generation |
| best anchor clique with `A` of positive-measure "reach" (diam `≤ 0.01`) | **0.998** | i.e. no violation at all |
| regularising task A's max clique (`clique_family.py regular`): `p` = its heaviest common point, `A` = intersection of the members that miss `p` | **0.9695** (best over the greedy path) | the rule "contains `p` **and** meets `A`" keeps only 54 of the 137 support poses through `p` |

**Reading.**  The certified 12.008 measure really does violate valid clique constraints (by
0.21–0.33) — that is not an artefact.  But every violated clique found is carried by *individual
poses*: as soon as the clique is required to be a **region** described by anchors, so that a
perturbed pose is still charged, the violation drops from 0.21 to at most 0.05, and to 0 once the
anchor is small enough for `{S : A ⊆ S}` to be more than one pose.  This is the same phenomenon
`CLIQUE_CEILING.md` measured as "grazing contacts", stated in a form that says *why*: the
grazing members are exactly the ones that a describable family cannot hold on to.

## 7. A near miss, recorded

The first version of the search dropped one separating axis.  A square and a *segment* are
disjoint iff one of **three** axes separates them — the square's two edge normals *and the
segment's own normal* — and only the first two were tested.  With that bug the scan reported an
anchor clique of mass **1.3308** on the certified measure at `p = (0.895, 0.995)` with a diagonal
anchor at offset `1/7`; the counterexample is the corner square `[0,1]²`, which contains `p` and
misses the anchor entirely (the anchor's endpoints are at `(1.0974, 0.9946)` and
`(0.8946, 1.1974)`, and the segment between them passes above the corner of `[0,1]²`: at `x = 1`
it is at `y = 1.0920`).  The missing axis is exactly the constraint `ε < 1 - p_x` of Lemma 2, and
without it the "family" is not a clique at all.  Two independent implementations of the
transversal test (a scalar polygon one and a vectorised one) and the exact certifier now agree on
`ρ*` to 12 digits, and the closed form of Lemma 2 is checked against both.
