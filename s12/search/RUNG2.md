# Rung 2 (`s(13) = 4` by a weighted closed cover of `[0,4]²`, `W < 13`) — task rung2-s13, 2026-09-11

**Verdict.**  Two results, one negative and one positive.

> **Theorem 1 (§2).**  If every leaf of a finite, closed, space-filling subdivision of the
> admissible pose space of `[0,m]²` (`m ≥ 4`) carries a *monotone witness certificate* — a fixed
> point set of weight `≥ 1`, every member of which lies in `Q` at every admissible pose of that
> leaf — then the total weight is `≥ m²`.
>
> `CORE`, `P1`, the `ADM` primitive of §3 and all their unions produce exactly such certificates.
> So **no cover of `[0,4]²` with `W < 16`, and in particular none with the rung-2 target `W < 13`,
> can be certified by them at any depth.**  The proof is an explicit dual: `m²` of the constraints
> with multiplier `1`, whose constraint sets are pairwise disjoint.  `TRI` is the one primitive in
> the checker that is not monotone (it is disjunctive), which is why rung 1 — Friedman's 14 points,
> `W = 14 < 16` — is certifiable at all and why it fails without `--tri`.

This supersedes `FAMILY.md` §2b's "engineering gap in `cert_p1`/`cert_core`" diagnosis and kills
the "scale the cover" and "solve the LP with a margin" routes outright (§5).  It also says what to
build, which is the second result:

> **`CHAIN` (§6)**, a weighted disjunctive primitive.  A monotone chain of pivot inequalities
> partitions a pose box into `k+1` regions, each with its own witness set; two chains of different
> kinds may be multiplied, with provably empty product regions.  Every test is an exact
> `Fraction` maximum of a polynomial that is affine in the centre and quadratic in `u`, so it is
> decided at four corners with no Bernstein slack and no recursion (Lemmas E–H).
>
> The disjunction it has to support is not small: at the worst rung-2 pose a single pivot reaches
> only `0.845828` of the needed `1`, while a chain cutting between every pair of consecutive
> coordinates (56 regions) reaches `1.013756`.  That measurement is what fixes the design.

Result on the target: see §4.  Everything load-bearing is `fractions.Fraction`; floats appear only
as pre-filters (documented supersets of the exact test, which can only lose certifications) and in
the independent stress tests.  Semantics unchanged (`ZEROMARGIN.md` §1: closed unit squares, closed
containment, a point on `∂Q` counts).

---

## 1. Notation

Container `[0,m]²`, `m = 4`.  Pose `(c, θ)`, `Q(c,θ) = c + R_θ[-½,½]²`, admissible iff
`c_x, c_y ∈ [w/2, m - w/2]` with `w(θ) = |cos θ| + |sin θ|`.  `u = tan(θ/2)`, so
`cos θ = (1-u²)/(1+u²)`, `sin θ = 2u/(1+u²)` are rational in `u`.  A *pose box* is
`[cx₀,cx₁] × [cy₀,cy₁] × [u₀,u₁]`.  For `p ∈ P` write `d = p - c` and

    X(c,θ) = d_x cos θ + d_y sin θ ,      Y(c,θ) = -d_x sin θ + d_y cos θ ,
    p ∈ Q(c,θ)   ⟺   |X| ≤ ½  and  |Y| ≤ ½ .

---

## 2. Theorem 1: the `m²` obstruction

Throughout, `m ≥ 4`, the leaves of the subdivision are closed boxes
`[a₀,a₁] × [b₀,b₁] × [u₀,u₁]` with `u₀ < u₁`, finitely many, covering the admissible pose space.
For a tile `(i,j) ∈ {0,…,m-1}²` write `m_{ij} = (i+½, j+½)` and `T_{ij} = [i,i+1] × [j,j+1]`, so
that `Q(m_{ij}, 0) = T_{ij}`, and let `clip(z,θ) = min(max(z, w/2), m - w/2)` with `w = w(θ)`.

### 2.1 The germ the checker must certify

**Lemma 0 (a leaf that sees one octant).**  Fix signs `σ_x, σ_y ∈ {±1}` with `σ_x = +1` when
`i = 0` and `σ_x = -1` when `i = m-1` (likewise `σ_y` in `j`), and put

    γ(t, s, θ) = ( clip(i+½ + σ_x t, θ),  clip(j+½ + σ_y s, θ),  θ ),      t, s, θ ≥ 0.     (γ)

Every `γ(t,s,θ)` is admissible by construction, and `γ(0,0,0) = (m_{ij}, 0)`.  Then some leaf `B`
contains `γ(t,s,θ)` for **all** `t, s, θ ∈ [0, η]`, for some `η > 0`.

*Proof.*  `γ` is continuous and the leaves are finitely many closed sets covering its image.  Write
`Λ_B = {(t,s,θ) ∈ [0,1]³ : γ(t,s,θ) ∈ B}`, a closed set; the `Λ_B` cover `[0,1]³`.  Suppose no `B`
contains a cube `[0,η]³`.  Then for every `B` and every `k` there is a point of `[0,1/k]³` outside
`Λ_B`; since there are finitely many leaves, some `B` and a subsequence... — more simply, for each
`k` pick `x_k ∈ [0,1/k]³` and a leaf `B_k ∋ γ(x_k)`; finitely many leaves means one leaf `B*`
occurs for infinitely many `k`, and `B*` is closed, so `γ(0,0,0) ∈ B*`.  Now `B*` is a **box** in
pose space: `a₀ ≤ i+½ ≤ a₁`, `b₀ ≤ j+½ ≤ b₁`, `u₀ = 0 < u₁`.  For the sign `σ_x = +1` the
infinitely many `x_k` with positive first coordinate force `a₁ > i+½`; when `i = 0` the same is
supplied by the clip, which pushes `c_x = w/2 > ½` to the right as soon as `θ > 0`; when
`i = m-1`, `σ_x = -1` and `a₀ < i+½` for the mirror reason.  Hence `γ([0,η]³) ⊆ B*` for `η` small
enough.  ∎

Consequently the witness set `S` of `B*` satisfies `w(S) ≥ 1` and

    S ⊆ PIN(i, j, σ) := ⋃_{η>0} ⋂_{0 ≤ t,s,θ ≤ η} Q(γ(t,s,θ)).                              (PIN)

`search/rung2_bound.py`'s `pin_mask` evaluates exactly this: the faces `t ∈ {0, η}`,
`s ∈ {0, η}` and four values of `θ` including `0`, with `η = 10⁻³` and `θ ≤ 10⁻³`.

### 2.2 What `PIN` is

**Lemma 1 (the cusp computation).**  Let `1 ≤ i, j ≤ m-2` (an interior tile), `d = p - m_{ij}`.
Then `p ∈ PIN(i,j,σ)` iff `p ∈ T_{ij}` and

  (a) `σ_x d_x < ½` **or** `d_x = σ_x ½` **and** the edge condition below holds — precisely:
      `p_x ∈ (i, i+1]` if `σ_x = +1`, `p_x ∈ [i, i+1)` if `σ_x = -1`; likewise in `y`;
  (b) on the vertical edge that survives (`x = i+1` when `σ_x = +1`, `x = i` when `σ_x = -1`):
      `σ_x d_y ≤ 0`, **strictly** when `σ_y = -σ_x`;
  (c) on the horizontal edge that survives (`y = j+1` when `σ_y = +1`, `y = j` when `σ_y = -1`):
      `σ_y d_x ≥ 0`, **strictly** when `σ_x = σ_y`;
  (d) the four corners of `T_{ij}` are excluded.

*Proof.*  On the interior of `T_{ij}` all four inequalities are strict at `(m_{ij},0)` and stay so
for small `t,s,θ`, so the open tile is in `PIN`.  For the rest write `c = γ(t,s,θ)`, so
`p - c = (d_x - σ_x t, d_y - σ_y s)` (`clip` is inactive for an interior tile and small `θ`), and
use `X = (d_x - σ_x t)\cos θ + (d_y - σ_y s)\sin θ`, `Y = -(d_x - σ_x t)\sin θ + (d_y - σ_y s)\cos θ`.

*(a)* Take `θ = s = 0`: `X = d_x - σ_x t ∈ [-½,½]` for all `t ∈ [0,η]` forces `σ_x d_x < ½` unless
`σ_x d_x = ½` is the *surviving* edge, i.e. `-σ_x d_x ≥ -½ + σ_x·(-σ_x t)`… concretely with
`σ_x = +1`: `d_x - t ≥ -½` for all small `t` gives `d_x > -½`, i.e. `p_x > i`, and `d_x - t ≤ ½` is
automatic.  Mirror for `σ_x = -1`.

*(b)* Surviving vertical edge, `σ_x = +1`, `d_x = ½`.  The binding inequality is `X ≤ ½`:
`(½ - t)\cos θ + (d_y - σ_y s)\sin θ ≤ ½`.  At `t = 0` this is
`d_y - σ_y s ≤ ½·(1-\cos θ)/\sin θ = ½\tan(θ/2)`.  If `σ_y = +1` the term `-σ_y s = -s` only helps,
so the binding case is `s = 0, θ ↓ 0`, giving `d_y ≤ 0` **closed**.  If `σ_y = -1` the condition is
`d_y + s ≤ ½\tan(θ/2)` for all `s, θ ≤ η`, and taking `s = η`, `θ ↓ 0` gives `d_y ≤ -η`, i.e.
`d_y < 0` **strict**.  Mirror for `σ_x = -1` (then the strict case is `σ_y = +1`).

*(c)* Surviving horizontal edge, `σ_y = +1`, `d_y = ½`; binding inequality `Y ≤ ½`:
`-(d_x - σ_x t)\sin θ + (½ - s)\cos θ ≤ ½`, i.e. `d_x - σ_x t ≥ -½\tan(θ/2) - s/\sin θ·…`; at
`s = 0` it is `d_x - σ_x t ≥ -½\tan(θ/2)`, so `σ_x = -1` (the `+t` helps) gives `d_x ≥ 0` closed and
`σ_x = +1` gives `d_x ≥ η`, strict.  Mirror for `σ_y = -1`.

*(d)* At a corner both `|d_x| = |d_y| = ½`.  Two of the four inequalities read
`±½\cos θ ± ½\sin θ ≤ ½` with the same sign pattern, i.e. `w(θ) ≤ 1`, false for `θ > 0`.  ∎

**Lemma 1' (wall tiles).**  For `i = 0` the clip is active, `c_x = w/2 = ½ + θ/2 + O(θ²)`.  Then
(i) the wall edge `x = 0` is entirely lost — `X ≥ -½` there reads `\cos θ·w(θ) ≤ 1`, false for
`θ > 0`; (ii) the interior edge `x = 1` is kept for the whole **open** edge — the exact condition is
`d_y ≤ ½ - u² + O(u³)` (`u = \tan(θ/2)`), the `O(u²)` loss instead of the `O(u)` of a fixed centre,
which is the same computation as Lemma A in §3.1 — with no restriction to a half.  (iii) The
horizontal edges obey (c) of Lemma 1 unchanged, with `i+½` read as `w/2 → ½`.  For the corner tile
`(0,0)` both clips are active and `PIN(0,0,+,+) = (0,1] × (0,1]` — in particular the tile corner
`(1,1)` **is** kept, which is exactly the `P1` lemma (`|c - p|_∞ = 1 - w/2`).  Mirror statements
hold at `i = m-1`, `j = 0`, `j = m-1`.

*(The "crude form suffices" sentence of the previous draft was wrong: the sets `I_i^{σ_x} × J_j^{σ_y}`
of (a) alone do **not** give 16.  Weight ½ at each of the 24 interior edge midpoints satisfies all
36 of those constraints with total 12 — the corner tile `(0,0)` gets `(1,½) + (½,1) = 1`, the edge
tile `(1,0)` with `σ = (+,+)` gets `(2,½) + (1½,1) = 1`, an interior tile gets two of its four
midpoints.  It is (b)–(d), the cusp conditions, that rule this out.)*

### 2.3 Disjointness, and the bound

**Lemma 2.**  Let `σ_x` depend only on the column index `i` and `σ_y` only on the row index `j`,
with `σ_x(0) = σ_y(0) = +1`, `σ_x(m-1) = σ_y(m-1) = -1`, and with the `+ → -` switch in each
direction occurring between two **interior** indices (possible iff `m ≥ 4`; e.g.
`σ = (+,+,-,-)` for `m = 4`).  Then the `m²` sets `PIN(i, j, σ_x(i), σ_y(j))` are pairwise disjoint.

*Proof.*  Let `p` lie in two of them, for tiles `(i,j) ≠ (i',j')`.  By Lemma 1(a) the `x`-ranges are
`I_i = (i,i+1]` where `σ_x(i) = +1` and `[i,i+1)` where `σ_x(i) = -1`; two such intervals are
disjoint unless `i' = i+1` with `σ_x(i) = +1` and `σ_x(i+1) = -1`, when they share exactly the line
`x = i+1`; likewise in `y`.  So either `i = i'` or the tiles are `x`-adjacent across the switch, and
similarly in `y`.  Four cases.

(i) `i = i'`, `j = j'`: excluded.
(ii) `i = i'`, `j' = j+1` across the `y`-switch, so `σ_y(j) = +1`, `σ_y(j+1) = -1` and
`p_y = j+1`.  By Lemma 1(c) tile `(i,j)` keeps its top edge with `d_x ≥ 0`, strict iff
`σ_x(i) = σ_y(j) = +1`; tile `(i,j+1)` keeps its bottom edge with `d_x ≤ 0`, strict iff
`σ_x(i) = σ_y(j+1) = -1`.  The two tiles are in the same column, so `σ_x(i)` is the same for both,
and exactly one of "`σ_x(i) = +1`" and "`σ_x(i) = -1`" holds: exactly one of the two conditions is
strict.  Their intersection is therefore `{d_x = 0}` minus one of them, i.e. empty.
(iii) `j = j'`, `i' = i+1` across the `x`-switch: the mirror of (ii), using Lemma 1(b) and the
common `σ_y(j)`.
(iv) both adjacencies, so `p = (i+1, j+1)`.  Then `p` is a corner of all four tiles, excluded by
Lemma 1(d).

Finally the switch must be interior: if it were at `i = 0` (i.e. `σ_x(1) = -1`) then tile `(0,j)`
keeps the *whole* open edge `x = 1` by Lemma 1'(ii) — no half — and tile `(1,j)` keeps its left
edge `x = 1` with `d_y ≥ 0`, and the two overlap.  (`rung2_bound.py dual --switch 0` prints the
overlap.)  ∎

**Theorem 1.**  Let `P ⊂ [0,m]²` (`m ≥ 4`) be a finite weighted point set, and suppose every leaf of
a finite, closed, space-filling subdivision of the admissible pose space carries a *monotone
witness certificate* — a set `S ⊆ P` with `w(S) ≥ 1` such that **every** `p ∈ S` lies in `Q(c,θ)`
at **every** admissible pose of the leaf.  Then `w(P) ≥ m²`.

*Proof.*  Fix the signs of Lemma 2.  By Lemma 0 each tile `(i,j)` has a leaf `B*` containing the
germ `(γ)`, and its witness set satisfies `w(PIN(i,j,σ)) ≥ w(S) ≥ 1` by (PIN).  The `m²` sets are
pairwise disjoint by Lemma 2, so `w(P) ≥ Σ_{ij} w(PIN(i,j,σ)) ≥ m²`.  ∎

Equivalently, in LP-duality language: `λ = 1` on those `m²` of the constraints
`w(PIN(i,j,σ)) ≥ 1` and `0` on the rest is a feasible dual of value `m²`, because the
multiplicity `Σ λ·1_{PIN}` is at most `1` at every point of `[0,m]²`.

**Why the check is finite and covers all covers.**  Membership in `PIN(i,j,σ)` depends only on the
comparisons of `p_x` with `i, i+½, i+1` and of `p_y` with `j, j+½, j+1` (Lemmas 1, 1'), so it is
constant on every cell of the arrangement of the half-integer grid.  Each cell contains a point of
the **quarter-integer lattice** (open 2-cells contain their centre `(k+½)/2`, relatively open
1-cells their midpoint, 0-cells are themselves quarter-integers).  So verifying
`multiplicity ≤ 1` at the `(4m+1)²` quarter-lattice points verifies it on all of `[0,m]²` — no
lattice restriction is imposed on the cover.

```
python3 search/rung2_bound.py dual --m 4
  sigma_x per column: [1, 1, -1, -1]     sigma_y per row: [1, 1, -1, -1]
  sum of lambda = 16
  max multiplicity over the 1/4-lattice (289 cell representatives) = 1
  VERDICT: the 16 PIN sets are pairwise disjoint => every monotone-witness-certifiable
           cover of [0,4]^2 has total weight >= 16
```
`--switch 0` (the switch at the wall adjacency) prints multiplicity `2`, as Lemma 2 predicts;
`--m 5 --switch 2`, `--m 6 --switch 2`, `--m 7 --switch 3` all print `1`, so the bound is `m²` for
`m = 4,…,7` and, by the same proof, for every `m ≥ 4`.

**Machine confirmation (the LP, over all 36 constraints).**  `rung2_bound.py bound --octants`
minimises the total weight subject to *all* the constraints `w(PIN(i,j,σ)) ≥ 1` (4 corner tiles
with one sign pair, 8 edge tiles with two, 4 interior tiles with four):

| lattice | candidates | constraints | LP value |
|---|---|---|---|
| 1/4 (a complete cell system) | 289 | 36 | **16.000000** |
| 1/20 | 6,561 | 36 | **16.000000** |
| 1/40 | 25,921 | 36 | **16.000000** |

and its optimal dual is exactly the integral certificate above.  Dropping the one-sidedness — using
only `S ⊆ Q(m_{ij},0) ∩ ⋂ Q(γ(t,t,0))` over *both* signs — the `16` sets overlap and the LP value
halves to `8` (`rung2_bound.py bound`, no `--octants`); it is Lemma 0, the fact that a leaf touching
the tile pose from one side is itself a leaf that must be certified, that gives `m²`.

**Machine confirmation 1 (the shipped scaled covers).**

```
python3 search/rung2_bound.py check runs/closed4_best_x103.txt
```
`container [0,4]^2, 1972 points, total weight 3197502027/250000000 = 12.790008108` —
**24 of the 36 PIN sets have weight below 1** (exact `Fraction` sums), the worst being
`tile (0,1) sx=+1 sy=-1: 5661189/10000000 = 0.566119`, i.e. the pose `(0.5, 1.5, 0)`.  The four
corner tiles pass at `1.03000` (Lemma 1', the `P1` situation).  `0.566119` is exactly the number a
brute-force float scan gives for "the weight captured simultaneously at every pose of an
arbitrarily small box at `(0.5, 1.5, 0⁺)`", and it does not move as the box shrinks from `10⁻²` to
`10⁻⁶` — the signature of the obstruction.  This is precisely where `FAMILY.md` §2b's depth ladder
was stuck (`cx ≈ 0.5`, `cy ≈ 1.5`, `θ → 0`): those boxes are not converging slowly, they are not
converging at all.

**Machine confirmation 2 (rung 1).**

```
python3 search/rung2_bound.py check runs/friedman14.txt
```
Friedman's 14 points, `W = 14 < 16`: **exactly 6 PIN sets are empty (weight 0)**, at the four
*interior* tile poses `(1.5,1.5,0)`, `(2.5,1.5,0)`, `(1.5,2.5,0)`, `(2.5,2.5,0)`.  Those six leaves
are exactly the ones `TRI` certifies, and they are why rung 1 **does not terminate without
`--tri`** (`ZEROMARGIN.md` §3: 87,200 uncertified at depth 21 with `CORE`+`P1`; with `ADM` and no
`--tri`, still 71,961 uncertified at depth 20).  Theorem 1 predicts the failure and names the poses;
the checker finds them.

**Corollary.**  The same proof gives `w(P) ≥ m²` at every rung of `search/FAMILY.md`: rung 2
(`m = 4`) needs `16` against a target of `13`, `m = 5` needs `25` against `21`, rung 3 (`m = 7`)
needs `49` against `45`.  *Every* rung of the `s(m²-3)` / `s(m²-4)` family is beyond the monotone
primitives, by exactly the `3` or `4` that makes the theorem interesting.  Closed-square semantics
is not an incidental convenience of these proofs — it is their whole content, and a machine proof
has to reproduce the case analysis (Bentz's 6-leaf tree, DS7's Lemma 3) in some form.  §6 builds
one.

---

## 3. `ADM`: a new exact primitive (subsumes `CORE` and `P1`)

Although it cannot beat Theorem 1, the right monotone primitive is still needed for every leaf that
*is* certifiable, and the old pair was measurably weak near the walls.  `ADM` replaces both.

### 3.1 Lemma A (monotone corners)

> **Lemma A.**  Let `B = [cx₀,cx₁] × [cy₀,cy₁] × [u₀,u₁]` with `u₀ ≥ 0`, `u₁ ≤ 1` (so
> `θ ∈ [0°, 90°]` and `cos θ, sin θ ≥ 0`).  For `θ` in the bin put
>
>     A_x(θ) = max(cx₀, w(θ)/2),   B_x(θ) = min(cx₁, m - w(θ)/2),   A_y, B_y likewise.
>
> If for **every** `θ` in the bin
>
>   (i)  `(p_x - A_x) cos θ + (p_y - A_y) sin θ ≤ ½`
>   (ii) `(p_x - B_x) cos θ + (p_y - B_y) sin θ ≥ -½`
>   (iii)`-(p_x - B_x) sin θ + (p_y - A_y) cos θ ≤ ½`
>   (iv) `-(p_x - A_x) sin θ + (p_y - B_y) cos θ ≥ -½`
>
> then `p ∈ Q(c,θ)` for **every admissible pose** `(c,θ) ∈ B`.

*Proof.*  At a fixed `θ`, the admissible centres of `B` are exactly
`c_x ∈ [A_x, B_x]`, `c_y ∈ [A_y, B_y]` (that is the definition of admissible, intersected with the
box).  `X(c,θ) = (p_x - c_x)cos θ + (p_y - c_y) sin θ` is non-increasing in `c_x` and in `c_y`
because `cos θ, sin θ ≥ 0`, so `X ≤ X(A_x,A_y) ≤ ½` by (i) and `X ≥ X(B_x,B_y) ≥ -½` by (ii).
`Y(c,θ) = -(p_x - c_x)sin θ + (p_y - c_y)cos θ` is non-decreasing in `c_x` and non-increasing in
`c_y`, so `Y ≤ Y(B_x,A_y) ≤ ½` by (iii) and `Y ≥ Y(A_x,B_y) ≥ -½` by (iv).  Hence `|X|,|Y| ≤ ½`. ∎

(The four inequalities are evaluated at four *different* corners; that is strictly stronger than
`CORE`, which tests all four corners against all four inequalities.  If `A_x > B_x` the admissible
set at that `θ` is empty and the conclusion is vacuous but the test is still sound, since each
inequality is checked at the correct extreme of the — possibly empty — interval.)

**This is what the brief's "off-centre wall witness" asks for, and more.**  Take a left-wall box,
`cx₀ ≤ w/2`, so `A_x = w/2`, and a witness at `p_x = 1` (the first interior grid line) at height
`D_y = p_y - A_y` above the box's bottom.  Inequality (i) reads, after clearing denominators,

    D_y ≤ u + (cos θ - sin θ)/2 = ½ - u² + O(u³).

The admissible tolerance in `D_y` is therefore `½` up to a **quadratic** error in the angle.
`P1`'s corresponding test is `|c - p|_∞ ≤ 1 - w/2 = ½ - u + O(u²)` — a **linear** loss.  That
linear loss is exactly the "slow convergence at wall boxes with off-centre witnesses" that
`FAMILY.md` §2b describes, and Lemma A removes it in closed form, for a witness pair at any
heights, with no symmetry assumption.  (`P1` is still kept: it certifies the corner region
`Q ⊆ p + [-1,1]²` in one box at any bin width, which Lemma A only matches after subdivision.)

### 3.2 Lemma B (the bounds are polynomials in `u`)

> **Lemma B.**  Write `c = (1-u²)/N`, `s = 2u/N`, `N = 1+u²`, and each centre bound as
> `X_n(u)/(2N)` with
>
>     'R' a  (a rational constant):  X_n = 2a(1+u²)
>     'W'    (`c_x ≥ w/2`):          X_n = (1-u²) + 2u
>     'M'    (`c_x ≤ m - w/2`):      X_n = 2m(1+u²) - (1-u²) - 2u
>
> and set `U = 2N p_x - X_n`, `V = 2N p_y - Y_n` (quadratics in `u`).  Then (i)–(iv) of Lemma A are
> equivalent to `G(u) ≤ 0` with, respectively,
>
>     G₁ = U·C + V·S - N²,   G₂ = -U·C - V·S - N²,
>     G₃ = -U·S + V·C - N²,  G₄ =  U·S - V·C - N²,        C = 1-u², S = 2u,
>
> polynomials of degree `≤ 4` with rational coefficients.  When both bounds are of kind `'R'`,
> `U = U₀N` and `V = V₀N` and `G` has the factor `N > 0`, leaving the quadratics
> `g₁ = (U₀-1, 2V₀, -U₀-1)`, `g₂ = (-U₀-1, -2V₀, U₀-1)`, `g₃ = (V₀-1, -2U₀, -V₀-1)`,
> `g₄ = (-V₀-1, 2U₀, V₀-1)`.

*Proof.*  `X = (U·C + V·S)/(2N²)` and `Y = (-U·S + V·C)/(2N²)` by direct substitution
(`p_x - c_x = (2N p_x - X_n)/(2N)` etc.), and `2N² > 0`. ∎

`max(cx₀, w/2)` is not a polynomial in `u`.  Both `cx₀` and `w(θ)/2` are separately *lower* bounds
for it, so verifying an inequality with either one is sound; the implementation offers `'R'`
always, adds `'W'` whenever `cx₀ < w_hi/2` (i.e. whenever the wall can bind on the bin), and
accepts an inequality if **any** offered combination certifies it.  That is why `ADM ⊇ CORE`: on a
box where no wall binds, every bound is `'R'`, the polynomial is a genuine quadratic, and Lemma C
below is exact — so `ADM` reduces to `CORE` with a centre rectangle clipped at least as tightly.

### 3.3 Lemma C (the exact max over the bin)

> **Lemma C.**  For `deg ≤ 2` the maximum of `G` on `[u₀,u₁]` is attained at an endpoint or, if
> the leading coefficient is negative, at the vertex `-a₁/(2a₂)` when it lies inside — an exact
> rational computation.  For `deg 3, 4`, write `G` in the Bernstein basis of degree 4 on `[u₀,u₁]`:
> `b_j = Σ_{k≥j} a_k C(k,j) u₀^{k-j} h^j` (`h = u₁-u₀`), `β_i = Σ_{j≤i} [C(i,j)/C(4,j)] b_j`.  Then
> `max_{[u₀,u₁]} G ≤ max_i β_i`, with equality at the endpoints (`β₀ = G(u₀)`, `β₄ = G(u₁)`), and
> the overestimate is `O(h²)`.

*Proof.*  A polynomial of degree `n` equals `Σ_i β_i B_{i,n}(t)` with `B_{i,n} ≥ 0`, `Σ B_{i,n} = 1`
(the convex-hull property of the Bernstein basis), so it is a convex combination of the `β_i`. ∎

Both branches are exact `Fraction` arithmetic; the Bernstein branch is a *sound* over-estimate, so
`ADM` never certifies a box it should not, and the `O(h²)` slack vanishes under subdivision.

### 3.4 Implementation and the θ-biased splitting rule

`search/zeromargin.py`:
* `cert_adm` — float pre-filter `_adm_mask` (numpy, the same formulas with a `1e-9` lenient
  tolerance, a superset of the exact test by construction) over the whole point set, then exact
  `Fraction` confirmation `_adm_exact` of the survivors in decreasing weight order.
* `cert_mix` — `ADM` and `P1` certify *different* points of the same box; their **union** is again
  a legitimate witness set, so a box neither carries alone can still be a leaf.  Reported as `MIX`.
* `--theta-bias K` (default 4): when the box touches a wall (`cx₀ < w_hi/2` or the three mirror
  conditions) and `θ` is small, the `u`-extent is weighted by `K` when choosing which dimension to
  halve.  This is the brief's item (a).  It is a small win on rung 1 (below) and does not change
  any verdict; the wall problem it was aimed at is solved analytically by Lemma A instead.
* `--no-adm` restores the old `CORE`-first path, byte-for-byte.

---

## 4. Measurements

**Rung 1 regression (`search/zeromargin.py friedman14 --tri --depth 14 --nproc 4`).**

| run | boxes | depth | ADM | CORE | P1 | MIX | CHAIN | TRI | EMPTY | uncert |
|---|---|---|---|---|---|---|---|---|---|---|
| `--no-adm --theta-bias 1` (the old code path) | 6,958 | 10 | – | 3,356 | 70 | – | – | 74 | 3,179 | **0** |
| `ADM`, `--theta-bias 1` | **6,810** | 10 | 3,367 | 0 | 0 | 0 | 0 | 74 | 3,164 | **0** |
| `ADM`, `--theta-bias 4` (default) | 7,220 | 10 | 3,511 | 0 | 0 | 0 | 0 | 74 | 3,225 | **0** |

The first row reproduces `ZEROMARGIN.md` §3 exactly (`3356/70/74/3179/0`, 6,958 boxes), so the old
path is untouched.  `ADM` alone replaces **both** `CORE` and `P1` (`P1` drops to 0 leaves: the 70
corner/wall boxes `CORE` could never do are now `ADM` leaves) and needs 148 fewer boxes.  Without
`--tri` it is still `NOT VERIFIED` (71,961 uncertified at depth 20) — as Theorem 1 requires, since
`14 < 16`.

**Rung 2, monotone primitives only** (`cert runs/closed4_best_x103.txt --depth 10 --nproc 8`,
`taskset -c 16-23`): 307,420 boxes, `ADM 117,805 / CORE 0 / P1 0 / MIX 0 / TRI 0 / EMPTY 4,563 /
UNCERTIFIED 34,542`, 321 s.  `FAMILY.md` §2b's depth-10 row for the same file was `300,294` boxes,
`114,471 / 337 / 0 / 4,327 / 34,212`.  `ADM` does strictly more per box and the uncertified count
is unchanged to 1 % — those boxes are the obstructed ones, and Theorem 1 says no amount of depth or
monotone cleverness removes them.  All 34,542 sit at the 24 violated PIN sets of §2.

**Rung 2 with `CHAIN`** (`--depth 14 --nproc 8 --disj --chain-from 5`, `taskset -c 8-15`).
CHAIN_RESULT_PLACEHOLDER

**Per-box behaviour of `CHAIN`** (direct calls, `runs/closed4_best_x103.txt`):

| box | `ADM` | `MIX` | `CHAIN` | time |
|---|---|---|---|---|
| `[0.5,0.6] × [1.4,1.5] × [0°,7.2°]` (the root box at the worst pose) | — | — | **CHAIN**, 197 points | 0.34 s |
| `[0.5,0.51] × [1.49,1.5] × [0°,0.11°]` | — | — | **CHAIN**, 196 points | 0.23 s |
| `[0.4875,0.5] × [1.4875,1.5] × [0°,1.8°]` | — | — | **CHAIN**, 196 points | 0.25 s |
| `[1.5,1.6] × [1.4,1.5] × [0°,1.8°]` (interior tile pose) | — | — | **CHAIN**, 402 points | 3.1 s |
| the four `(±,±)` octants at `(1.5,1.5,0)`, bins down to `u₁ = 2⁻¹²` | — | — | **CHAIN** (product of two chains) | 2.7–3.2 s |

The interior-tile boxes need the *product* of two chains and the empty-region test of Lemma H; a
single chain fails there, which is what forced the two-chain extension.

**Independent float stress (`search/zeromargin_stress.py`, no shared code path).**
`runs/zm_f14_adm.txt` (6,810 `ADM` + 74 `TRI` leaves), 40 sampled poses per leaf including the box
corners: **0 failures**.  Primitive tests: `P1` 200,000 instances, core lemma 68,569 (21 angles
each), triangle 72,358 — 0 failures.  `ADM` (Lemma A) on 40,000 random `(box, point)` pairs with
the four inequalities evaluated directly in floats on a 201-point `θ` grid: 2,713 certifying pairs,
43,859 pose samples, **0 failures**.  `CHAIN`: the violation polynomials against the geometry on
200,000 random `(point, pose)` pairs — **0 disagreements**; the `_gmax` enclosure (Lemma E) on
20,000 random `(box, nonnegative combination)` pairs × 30 interior poses — **0 violations**.

---

## 5. Two routes that are now closed, with the numbers

**(a) Scaling the cover (`FAMILY.md` §2b) cannot work.**  Captured weight is linear in the weights,
so `×λ` multiplies every PIN weight by `λ`.  The worst PIN set of `closed4_best.txt` is
`0.549...`; closing it needs `λ ≥ 1.82`, i.e. `W ≥ 22.6`.  `×1.02` and `×1.03` move the worst PIN
from `0.5496` to `0.5661` — nowhere near 1, which is why the depth ladder oscillated instead of
converging.  The exact `pose` spot-checks in `FAMILY.md` §2b (`1.0655`, `1.0759`) were taken at
`cy = 1.4969`, a *non-degenerate* pose; the degenerate pose is `cy = 1.5` exactly, where the base
cover captures `1.0000072` — a genuinely tight family that the spot check missed.

**(b) Solving the cover LP with a uniform geometric margin cannot work either.**  A cover with
margin `μ` (a point counts only if it is `≥ μ` inside every edge) is automatically monotone-witness
certifiable, so Theorem 1 forces `≥ 16` — and the LP confirms it directly.
`search/closed4.py` now takes `--margin MU` (sets the containment tolerance to `-MU`; the two
`tol=TOL` default arguments were rebound so the module global actually takes effect).

```
python3 search/closed4.py axis --s 4 --margin 0.01     ->  AXIS-ALIGNED LP = 16.000000 (converged)
python3 search/closed4.py run  --s 4 --margin 0.010    ->  LP = 16.000000 from it1 on
python3 search/closed4.py run  --s 4 --margin 0.005    ->  LP = 16.000000 from it1 on
```
(The axis-only run converges in 8 iterations with zero violated poses, so `16` there is the honest
sampled-LP optimum, not an unconverged value.  The two full-domain runs were stopped at iteration 6
after ~930 s, not run to convergence; their value is flat at exactly `16.000000` from iteration 1
on while the row count grows `24k -> 47k` and column generation keeps adding points, and the LP is
bounded above by `16` anyway — the all-ones weighting of the `4x4` integer lattice is feasible — so
`16` there is the optimum, not an unconverged upper estimate.)  Equivalently, by scaling, the cover LP for the
container `[0, 4/(1-2μ)]²` jumps to `≥ 16` as soon as `μ > 0`: `COVER(s)` is **discontinuous at
`s = 4`**, because `[0,4]²` is exactly tiled by 16 unit squares and the closed semantics is what
lets the 16 tiles share their boundary weight.  Any hope of buying robustness by working at
`s = 4 + δ` and scaling back dies here.

A bug found on the way: `closed4.py run`'s per-round checkpoint called `export(m, x, …)` with an
`x` from before that round's column generation and crashed with `IndexError` (`index 863 is out of
bounds`) the first time pricing added a column.  Fixed by zero-padding `x` to the widened model.

---

## 6. `CHAIN`: the disjunctive primitive

Theorem 1 says a leaf's certificate must sometimes be *disjunctive*: a partition of the leaf into
regions, each with its **own** witness set.  `TRI` is the case "three regions, one point each, each
of weight `≥ 1`" — useless when the weights are `≈ 0.01`.  `CHAIN` is the weighted version.

### 6.1 The shape of the disjunction, measured

At the worst rung-2 pose, `(0.5, 1.5, 0)` with `σ = (+,-)`, the square is `[0,1] × [1,2]` and its
capture structure over a small box is a **sliding cut**.  With `c_x = w/2 + …`, `c_y = 3/2 - η`,
`θ` small, the exact conditions are

    point (x, 1) is captured  ⟺  (x - c_x)\sin θ ≤ ½ - (c_y - 1)\cos θ  ⟺  x ≲ ξ,
    point (x, 2) is captured  ⟺  (c_y - 2)\cos θ + ½ ≥ -(x - c_x)\sin θ  ⟺  x ≳ ξ,
    ξ = c_x + η/θ + O(θ).

So the row `y = 1` is captured to the left of a cut and the row `y = 2` to the right of the *same*
cut, and `ξ` sweeps all of `[0,1]` as `(η, θ)` range over any neighbourhood of `(0,0)` — which is
why no box subdivision controls it, and why this is exactly the residue Theorem 1 predicts.

Measured on `runs/closed4_best_x103.txt` at that box (`runs/` scratch script, dense float scan;
the stable part weighs `w(T) = 0.566119`, the two rows `0.525992` and `0.498610`):

| disjunction | worst region |
|---|---|
| a single pivot (k = 2 regions), best cut | **0.845828** — fails |
| cuts between every pair of consecutive `x`-coordinates (56 regions) | **1.013756** — works |

A two- or three-region disjunction is therefore *not* enough; the primitive must support a chain of
arbitrary length.  It is cheap because the regions are nested: they are the sign pattern of one
monotone family of polynomials.

### 6.2 The primitive

For a point `p` and a pose `(c_x, c_y, u)`, with `a = p_x - c_x`, `b = p_y - c_y`, `C = 1-u²`,
`S = 2u`, `N = 1+u²`, define the four **violation polynomials**

    G_{p,0} = 2aC + 2bS - N,   G_{p,1} = -2aC - 2bS - N,
    G_{p,2} = -2aS + 2bC - N,  G_{p,3} = 2aS - 2bC - N,

so that `p ∈ Q(c,θ)` iff `G_{p,k} ≤ 0` for all four `k` (multiply `|X| ≤ ½`, `|Y| ≤ ½` by `2N > 0`).

> **Lemma E (exact maximum).**  Every nonnegative combination `Σ λ_r G_{p_r, k_r}` is **affine** in
> `(c_x, c_y)` and **quadratic** in `u`.  Hence its maximum over a pose box is
> `max` over the four corners of the centre rectangle of the exact quadratic maximum over
> `[u₀,u₁]` — a finite exact `Fraction` computation, with no Bernstein slack and no subdivision.

*Proof.*  `a` and `b` are affine in `c_x, c_y` and appear linearly; `C, S, N` are quadratic in `u`
with `a, b` not involving `u`.  An affine function of `(c_x,c_y)` attains its maximum over a
rectangle at a corner, for each fixed `u`, and `max_u max_{corners} = max_{corners} max_u`. ∎

> **Lemma F (chain).**  Let `q_1, …, q_k` be points with chosen kinds such that
> `max_B (G_{q_r} - G_{q_{r+1}}) ≤ 0` for `r = 1, …, k-1` (so `G_{q_1} ≤ … ≤ G_{q_k}` on `B`, by
> transitivity).  Then the `k+1` sets
>
>     R_0 = {G_{q_1} > 0},  R_r = {G_{q_r} ≤ 0 < G_{q_{r+1}}} (1 ≤ r < k),  R_k = {G_{q_k} ≤ 0}
>
> partition `B`, and on `R_r` every `q_j` with `j ≤ r` satisfies its inequality.

*Proof.*  Given a pose, let `r` be the largest index with `G_{q_r} ≤ 0` (`r = 0` if none): the pose
is in `R_r` and in no other.  On `R_r`, `G_{q_j} ≤ G_{q_r} ≤ 0` for `j ≤ r`. ∎

> **Lemma G (the opposite side).**  If `max_B (G_a + λ G_q) ≤ 0` for some `λ > 0`, then at every
> pose of `B` with `G_q > 0` one has `G_a < 0`.

*Proof.*  `G_a ≤ -λ G_q < 0`. ∎

> **Lemma H (empty product regions).**  If `max_B (G_{q} + λ G_{q'}) ≤ 0` for some `λ > 0`, then no
> pose of `B` has `G_q > 0` and `G_{q'} > 0` simultaneously.

*Proof.*  Both positive would give `G_q + λ G_{q'} > 0`. ∎

**The primitive.**  `cert_chain` (`search/zeromargin.py`):

1. `T` = everything `ADM`/`P1` certify for the whole box.  A *swing* point is a point of weight `>0`
   within reach of the box, not in `T`, exactly one of whose four inequalities fails the `ADM` test
   — its `G_p` is the swing polynomial.  (Screening: a condition whose *float* `ADM` bound already
   fails is exactly failed, so points with two or more float-failing conditions are discarded with
   no `Fraction` arithmetic at all.)
2. Group the swing points by kind; within a group, sort by a float proxy and build a chain greedily,
   keeping a point only when the exact test of Lemma F passes against the previous one.
3. For every candidate `a`, binary-search the largest `r` with "Lemma G certifies `a` from
   `q_r`" — the set of such `r` is a prefix because `G_{q_r}` is non-decreasing along the chain —
   using `λ ∈ {1, ½, 2}`.
4. Require `w(T) + w({q_1,…,q_r}) + w(U_{r+1}) ≥ 1` for every region `r`.

**Two chains.**  At a *wall* pose one cut slides and one chain suffices.  At an *interior* tile pose
two cuts slide independently — with `c = (3/2 + α, 3/2 + β)` the row `x = 2` is captured for
`y ≲ c_y + α/θ` and the row `y = 2` for `x ≳ c_x - β/θ`, and `α, β` are independent — so `CHAIN`
also tries the **product** of two chains of different kinds, with regions `R_r × R'_s`.  A product
region can be empty, and Lemma H certifies that: at the interior tile pose the two lowest pivots
cannot both be violated (the algebra reduces to
`G + G' = 2[(4 - c_x - c_y)C + (y_1 - x_2 + c_x - c_y)S - N] ≤ 0`), which is what makes the
`(0,0)` region — the one with no witness at all — provably empty.  The empty set is again a
staircase in `(r,s)` (both `G` families are monotone), so one binary search per `r` finds it.

**Everything load-bearing is `Fraction`**: Lemmas E–H are decided by `_gmax`, which is exact; the
only floats are the reach/screen pre-filters and the sort proxy, all of which can only *lose*
certifications, never create one.

### 6.3 Independent check

`search/zeromargin_stress.py` (floats, no shared code path) now also tests:

* the **violation polynomials**: 200,000 random `(point, pose)` pairs, comparing
  `max_k G_{p,k} ≤ 0` with the geometric `p ∈ Q(c,θ)` — **0 disagreements**;
* the **`_gmax` enclosure** (Lemma E): 20,000 random `(box, nonnegative combination)` pairs, the
  exact `Fraction` bound against 30 random interior poses each — **0 violations**;
* every `CHAIN` leaf of a dump: at 40 sampled admissible poses per leaf, the weight captured *from
  the leaf's recorded witness list* must reach `1` (a different subset does it in each region; the
  stress test knows nothing about the regions) — see §4 for the counts.

---

## 7. What is exact, what is float

| quantity | status |
|---|---|
| Lemmas A, B, C and every `ADM`/`CORE`/`P1`/`TRI` leaf | exact `Fraction`; floats only as a documented superset pre-filter |
| Lemmas E, F, G, H and every `CHAIN` leaf | exact `Fraction` (`_gmax`: four rectangle corners × the exact quadratic maximum in `u`).  Floats appear only in the reach filter, the per-condition screen (a *float* failure implies an exact failure, so the screen only discards non-candidates) and the sort proxy — all of which can only lose certifications |
| Lemmas 0, 1, 1', 2 and Theorem 1 | proved on paper (§2); the disjointness claim is additionally machine-verified at every cell of the half-integer arrangement via the quarter-lattice (`rung2_bound.py dual`), which is exhaustive over `[0,m]²`, not a lattice restriction |
| `rung2_bound.py check` | weights summed in `Fraction`; PIN *membership* is decided by evaluating the pose inequalities in floats at `θ ∈ {0, 10⁻⁵, 10⁻⁴, 10⁻³}` with tolerance `10⁻¹²` — the cover's coordinates lie on a `1/1000` lattice, so the `O(θ)` and `O(θ²)` gaps are resolved with `≥ 50×` margin |
| the LP values `16.000000` (`rung2_bound.py bound`, `closed4.py --margin`) | float LP (HiGHS); its optimal dual is the integral certificate proved by hand in §2.3, so the LP is confirmation, not evidence |
| `0.566119`, `0.845828`, `1.013756` | the first exact (`Fraction`); the two disjunction numbers are float scans of the cover (design measurements, not part of any proof) |
| leaf counts in §4 | exact (the checker's own accounting) |
| `zeromargin_stress.py` | floats by design — an *independent* confirmation, not a proof |

---

## 8. Reproduce

```bash
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-11/closed4_best*.txt runs/

# Theorem 1: the explicit dual certificate, exhaustively verified (seconds)
python3 search/rung2_bound.py dual --m 4                 # max multiplicity 1 => bound 16
python3 search/rung2_bound.py dual --m 4 --switch 0      # the wall adjacency: multiplicity 2
python3 search/rung2_bound.py dual --m 5 --switch 2      # 25 ;  --m 7 --switch 3 -> 49
# the LP behind it (1/4 lattice is a complete cell system; 1/20, 1/40 confirm)
python3 search/rung2_bound.py bound --m 4 --pitch 20 --octants   # 16.000000
python3 search/rung2_bound.py bound --m 4 --pitch 20             # 8.000000 (two-sided only)

# Theorem 1 applied to the covers (exact)
python3 search/rung2_bound.py check runs/closed4_best_x103.txt   # 24/36 PIN sets < 1
python3 search/rung2_bound.py check runs/friedman14.txt          # 6 empty PIN sets = the TRI poses

# rung 1 regression, the ADM primitive, the independent stress test
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4 --no-adm --theta-bias 1
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4 --dump runs/zm_f14_adm.txt
python3 search/zeromargin_stress.py runs/zm_f14_adm.txt 40       # 0 failures everywhere

# rung 2 with ADM alone (5 min on 8 processes) and with ADM + CHAIN
taskset -c 16-23 python3 search/zeromargin.py cert runs/closed4_best_x103.txt --depth 10 --nproc 8
taskset -c 8-15  python3 search/zeromargin.py cert runs/closed4_best_x103.txt --depth 14 \
                 --nproc 8 --disj --chain-from 5 --dump runs/x103_disj_leaves.txt
python3 search/zeromargin_stress.py runs/x103_disj_leaves.txt 40 --cert runs/closed4_best_x103.txt

# the margin LP (the route Theorem 1 closes)
python3 search/closed4.py axis --s 4 --margin 0.01 --tag ax010 --nproc 1    # = 16
```

Files: `search/zeromargin.py` (`ADM`, `MIX`, `CHAIN`, `--disj`, `--chain-from`, `--theta-bias`,
`--no-adm`), `search/zeromargin_stress.py` (the `ADM` and `CHAIN` random tests, `--cert` for
weighted covers), `search/rung2_bound.py` (new: `bound`, `check`, `dual`),
`search/closed4.py` (`--margin`, checkpoint-export fix).
