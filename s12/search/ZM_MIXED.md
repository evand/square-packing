# Exact checker for mixed covers (points + segments + polygons): `zm_mixed.py` (2026-09-27)

Brief: `tasks/line-cover/B.md`.  Code: `search/zm_mixed.py` (checker), `search/zm_mixed_test.py` (tests, toy covers,
leaf stress).  Format: `tasks/line-cover/FORMAT.md`; reader `search/mixed_cover.py` (agent A).  `zeromargin.py` is
imported, not edited (sha256 `640fe453…086ab`, printed and compared at every run).  Labels: **certified** (exact),
**heuristic** (float), **estimate**.

## 0. Results

* **`search/zm_mixed.py`**: exact checker for `FORMAT.md` files.  It imports `zeromargin.py` (sha pinned) and
  certifies each pose box by a lower bound `L` on the segment/polygon mass plus zeromargin's own point primitives,
  fed `L` as a phantom point.  Every lemma is proved in §2; everything that certifies is `Fraction` arithmetic.
* **Germs** (the hard place of the brief) are handled by **Lemma T**: the two lines of a germ pair share one
  threshold `T` at every pose (`t ≥ T` on one line, `t ≤ T + u₀` on the other, including `θ = 0` with `T = ±∞`), so
  the pair's mass is `≥ min_T` of an explicit piecewise-linear function, computed exactly.  A uniform density has no
  gap between the pivots, so a germ box closes with `O(box)` loss.  Toy grid cover: a germ cell's 4 germ roots close
  in **12 boxes**; with Lemma T off, 1,340 boxes stay open at depth 16 (Theorem-1-type obstruction).  Points on the
  germ lines can join the groups (Corollary T′), which is what closes the germs of agent A's cover.
* **Off the germs**, the brief's core bound (Lemma S) terminates but costs `ε^{−2.2}` boxes at margin `ε`.
  **Lemma L** (chord ends are affine in the centre at fixed `u` ⇒ the lines' mass bound is concave ⇒ exact at
  the 4 corners, Bernstein in `u`) makes it first order: the 2 % cell drops from 31,278 boxes to 320, and the count
  is logarithmic in the margin down to 0.25 %.  **T1m025** (grid lines, true minimum `1.0025`, total 14.52):
  **VERIFIED-D4, 42,300 boxes, 884 CPU-s**.
* **Tests** (§4): point-only files give the same boxes and leaf census as `zeromargin.py cert` (two `s(13)`
  columns); randomised exact tests of every lemma (≈62k checks, 0 violations); rejection tests (remove / lighten /
  shift a segment off a grid line / zero and negative margin) all fail exactly around their holes — the shift test's
  hole is a germ pivot pose at `θ = 0.57°` — and leaf stress tests of every certified dump find no pose below 1.
* **Performance**: comparable line vs point covers (§4.6): 17,626 boxes / 304 CPU-s (lines) vs 19,998 / 606 CPU-s
  (the same lines as 2,406 points, `zeromargin.py`); germ cell 1,040 boxes / 19 s vs 1,932 / 86 s.
* **Agent A's `m = 4` mixed covers** (scaled by me, checker tests only): `r5 × 1.05` (total 12.929)
  **VERIFIED-D4** (27,648 boxes, 1.6 CPU-h, final code); `r7 × 1.02` (12.566) NOT VERIFIED, 74 boxes left at depth
  18 where the cover has `≥ 1.5 %` float margin: points and segments anti-correlated inside a box, which the
  one-number coupling (Lemma P) cannot see (§5, the main open item for tight mixed covers).


## 1. Statement and structure

Cover: point masses, segments carrying mass `w` uniformly by length, convex polygons carrying `w` uniformly by area.
Claim checked: every closed unit square `Q ⊆ [0,s]²` (any centre, any angle) has `μ(Q) ≥ 1`, a segment on `∂Q`
counting in full.  Poses, boxes, admissibility, root grids, `clip_bin` and the θ-biased splitting are exactly those
of `zeromargin.py` (`ZEROMARGIN.md` §2, `RUNG2.md` §3–4).  For each box the checker computes

1. **`L`**, a certified lower bound on the segment + polygon mass captured at **every admissible pose** of the box
   (Lemmas S, T, T′, L, L′ below); `L ≥ 1` closes the box (leaf `PIECE`);
2. otherwise it runs zeromargin's point primitives `ADM → P1 → MIX → CHAIN` unchanged, with one extra "phantom"
   point of weight `L` declared to be in `Q` at every admissible pose (Lemma P).

Children inherit the parent's `L` (a sub-box has fewer poses) and zeromargin's per-subtree witness cache.

**Notation.**  Pose `(c, θ)`, `u = tan(θ/2) ∈ [u₀,u₁] ⊂ [0,1)`, `Q(c,θ) = c + R_θ[−½,½]²`.  For a point `p` put
`X = (p−c)·(cos θ, sin θ)`, `Y = (p−c)·(−sin θ, cos θ)`; `p ∈ Q` iff the four inequalities
`(0) X ≤ ½, (1) X ≥ −½, (2) Y ≤ ½, (3) Y ≥ −½` hold.

## 2. The lemmas

**Lemma A/B/C** (`RUNG2.md` §3, used as is).  For a bound choice `χ` (each centre coordinate bounded by the box
side `'R'` or by the wall `'W'/'M'`, as `zeromargin.Checker._adm_specs` offers) and inequality `k`, there is a
polynomial `G_{k,χ}(u; p)` of degree `≤ 4` in `u` such that `G ≤ 0` on `[u₀,u₁]` implies inequality `k` for `p` at
every admissible pose of the box; `max G ≤ max_i β_i`, the degree-4 Bernstein coefficients on `[u₀,u₁]`.

**Fact 0 (affinity).**  `G_{k,χ}(u; p)` is affine in `p = (p_x, p_y)` (`U`, `V` of Lemma B are affine in `p`, and
`G` is affine in `U`, `V`), and the Bernstein transform is linear.  So each `β_i` is an affine function of `p`.

**Fact 1 (convexity).**  For a fixed pose, inequality `k` is a closed half-plane in `p`.  Hence the set
`C_k = {p : inequality k holds at every admissible pose of the box}` is convex and closed (an intersection of
half-planes), and so is `C = C₀ ∩ C₁ ∩ C₂ ∩ C₃ = {p : p ∈ Q at every admissible pose}`.

> **Lemma S (certified core of a line / a polygon).**  (a) Let `ℓ = {P₀ + t d}` be a line.  For each choice `χ`,
> `I_{k,χ} = {t : β_i(P₀ + t d) ≤ 0, i = 0..4}` is an interval with rational endpoints (or empty or unbounded), and
> `I_k := conv ⋃_χ I_{k,χ}` satisfies `P₀ + I_k d ⊆ C_k`.  With `J₄ = I₀ ∩ I₁ ∩ I₂ ∩ I₃`, every segment `σ ⊂ ℓ`
> of mass `w` has `μ_σ(Q) ≥ w·|σ ∩ J₄|/|σ|` at every admissible pose of the box.
> (b) Likewise `K_{k,χ} = {p ∈ F : β_i(p) ≤ 0}` (`F` a square containing every `Q` of the box) is a convex polygon,
> `K_k = conv ⋃_χ K_{k,χ} ⊆ C_k`, `K = ⋂ K_k ⊆ C`, and a polygon `P` of mass `w` has `μ_P(Q) ≥ w·area(P∩K)/area(P)`.

*Proof.*  (a) By Fact 0 each `β_i(P₀+td) = b_i + t e_i` is affine in `t`, so `I_{k,χ}` is an intersection of five
half-lines (or all of `ℝ`, or empty when `e_i = 0 < b_i`).  For `t ∈ I_{k,χ}` Lemma C gives `G_{k,χ} ≤ 0` on the bin
and Lemma A/B gives inequality `k` at every admissible pose, i.e. `P₀ + td ∈ C_k`.  `C_k ∩ ℓ` is convex (Fact 1),
so it contains the hull of the union.  Points of `σ ∩ J₄` are in `C`, i.e. in `Q` at every admissible pose, and
`μ_σ` is `w/|σ|` times length.  (b) The same with `β_i` affine in `p` (half-planes), and `F ⊇ Q` so `K ∩ F`
loses nothing: `F = [cx₀ − r, cx₁ + r] × [cy₀ − r, cy₁ + r]`, `r = 0.7072 > √2/2`. ∎

Lemma S is continuous in the pose except where `Q` has an edge on the segment's line: at an axis-parallel square
with an edge on a grid line the whole unit piece of line counts, at `θ = 0+` only a fraction depending on
(edge offset)/θ.  A box touching such a pose gets from Lemma S only the part certified for all its poses.  At a
**tile germ** (all four edges on grid lines, `θ → 0`) this is about half of the four lines' mass (`RUNG2.md`
Theorem 1 is exactly this obstruction for points), so Lemma S alone cannot terminate there; §4 measures it.

> **Lemma T (threshold lemma for a germ pair).**  Let `ξ` be rational, let `ℓ↑ = {x = ξ}` and `ℓ↓ = {x = ξ + 1}`,
> both parametrised by `t = y`.  For every pose with `θ ∈ [0°, 90°)` and `u ≥ u₀ ≥ 0` there is `T ∈ [−∞, +∞]` with
>
>  * every point `(ξ, t)` with `t ≥ T` satisfies inequality (1), and
>  * every point `(ξ + 1, t)` with `t ≤ T + u₀` satisfies inequality (0).
>
> The same holds for `ℓ↑ = {y = η + 1}` with inequality (2), `ℓ↓ = {y = η}` with inequality (3), `t = x`.

*Proof.*  Vertical pair.  If `θ > 0` (`s = sin θ > 0`, `c = cos θ`): inequality (1) at `(ξ, t)` reads
`(ξ − c_x)c + (t − c_y)s ≥ −½ ⇔ t ≥ T := c_y + (−½ − (ξ − c_x)c)/s`.  Inequality (0) at `(ξ + 1, t)` reads
`(ξ + 1 − c_x)c + (t − c_y)s ≤ ½ ⇔ t ≤ c_y + (½ − (ξ + 1 − c_x)c)/s = T + (1 − c)/s = T + u`, and `u ≥ u₀`.
If `θ = 0`: (1) at `(ξ,t)` is `c_x ≤ ξ + ½` and (0) at `(ξ+1,t)` is `c_x ≥ ξ + ½`, both independent of `t`; take
`T = −∞` if `c_x ≤ ξ + ½` (the first claim is then "(1) holds for all t", true; the second is vacuous) and `T = +∞`
otherwise (first vacuous, second "(0) holds for all t", true).  Horizontal pair, `θ > 0`: (3) at `(t, η)` is
`−(t − c_x)s + (η − c_y)c ≥ −½ ⇔ t ≤ T' := c_x + ((η − c_y)c + ½)/s`; (2) at `(t, η+1)` is
`t ≥ c_x + ((η + 1 − c_y)c − ½)/s = T' − u`.  Put `T = T' − u₀`: then `t ≥ T ⇒ t ≥ T' − u` (since `u ≥ u₀`) and
`t ≤ T + u₀ = T'`.  `θ = 0`: (3) is `c_y ≤ η + ½`, (2) is `c_y ≥ η + ½`; `T = +∞` resp. `−∞` as before. ∎

**Corollary T (the group bound).**  For the box, let `J₃↑ = ⋂_{k≠1} I_k(ℓ↑)` (the three other inequalities, Lemma
S) and `S↑ = I₁(ℓ↑)`; if `S↑ = [t↑, ∞)` put `τ↑ = t↑` (`−∞` if `S↑ = ℝ`), otherwise `τ↑ = +∞`.  Symmetrically
`J₃↓ = ⋂_{k≠0} I_k(ℓ↓)`, `S↓ = I₀(ℓ↓)`, `τ↓ = t↓` if `S↓ = (−∞, t↓]`, else `−∞`.  Define, for the segments
`σ ⊂ ℓ↑` resp. `ℓ↓` with densities `ρ_σ` (mass per unit length),

    F↑(a) = Σ_σ ρ_σ |σ ∩ J₃↑ ∩ [a, ∞)|,      G↓(b) = Σ_σ ρ_σ |σ ∩ J₃↓ ∩ (−∞, b]|,
    f(T)  = F↑(min(T, τ↑)) + G↓(max(T + u₀, τ↓)).

Then at every admissible pose of the box the pair's mass is `≥ min_{T ∈ [−∞,∞]} f(T)`, and this minimum is attained
at one of the finitely many breakpoints or on a constant tail.

*Proof.*  Fix a pose and its `T` from Lemma T.  A point of `ℓ↑` with parameter `t ∈ J₃↑` satisfies (0), (2), (3) at
this pose (Lemma S(a)); if also `t ≥ T` it satisfies (1) (Lemma T), and if `t ≥ τ↑` it satisfies (1) by
`S↑ ⊆ C₁` — in both cases it lies in `Q`.  So `Q ∩ ℓ↑ ⊇ J₃↑ ∩ [min(T, τ↑), ∞)` (when `τ↑ = +∞` the second case is
empty and only the first is used).  Symmetrically `Q ∩ ℓ↓ ⊇ J₃↓ ∩ (−∞, max(T + u₀, τ↓)]`.  Hence the pair captures
`≥ f(T) ≥ min f`.  With the conventions `F↑(−∞) = F↑` total, `F↑(+∞) = 0` etc. this includes `T = ±∞`.  Each
`F↑`, `G↓` is continuous, piecewise linear in its argument (a sum of clamped linear functions of the segment
endpoints), and `min(·, τ)`, `max(·, τ)` preserve that; so `f` is continuous and piecewise linear with breakpoints
among `{e : e endpoint of σ ∩ J₃↑, e ≤ τ↑} ∪ {τ↑}` and `{e − u₀ : e endpoint of σ ∩ J₃↓, e ≥ τ↓} ∪ {τ↓ − u₀}`, and
constant left of the smallest and right of the largest.  A continuous piecewise-linear function attains its
infimum at a breakpoint or on a constant tail; the code evaluates `f` exactly at all breakpoints and one point on
each tail. ∎

**Corollary T′ (line points in the group).**  Let the points of the cover lying exactly on `ℓ↑` or `ℓ↓` with
parameter in `J₃` be added to `f` as steps: a point `(t_p, w)` of `ℓ↑` contributes `w·[min(T, τ↑) ≤ t_p]`, one of
`ℓ↓` contributes `w·[t_p ≤ max(T + u₀, τ↓)]`.  Then the pair's segment mass **plus the mass of these points** is
`≥ inf_T f`.  *Proof.*  Same as Corollary T: such a point satisfies the other three inequalities (`t_p ∈ J₃`) and
the singular one when `t_p ≥ T` or `t_p ≥ τ↑` (resp. `≤`).  `f` is now `PL + STEP` with both parts changing only
at the breakpoints (segment ends, `t_p`, `τ`, all shifted by `u₀` on `ℓ↓`); on an open interval between two
consecutive breakpoints STEP is constant and PL linear, so the infimum there is STEP(midpoint) + the smaller PL
value at the two ends; the code takes the minimum of that over all intervals (and tails) and of `f` at the
breakpoints. ∎  Used as a **second attempt** only (the first attempt leaves all points to zeromargin): the points
counted in a group are withheld from the point primitives for that box (their three weight arrays set to 0 and
restored afterwards) — otherwise they would be counted twice — and a point on two group lines (a grid crossing)
is given to one group only.  This is the piece that closes tile germs of covers that keep point mass on the
grid lines next to segments (§4.4): zeromargin's `CHAIN` has to multiply two pivot chains there, while here
the two orientations are separate one-dimensional minimisations.

The bound for a group is `max(min f, core)` where core is the Lemma S bound of the same two lines (both are valid);
the groups and all remaining lines and polygons involve **disjoint sets of pieces**, so their bounds add.  Any piece
may be dropped from any bound (masses are `≥ 0`), which is what makes the float "reach" filter (pieces farther than
`0.7072 +` half the box diagonal from the box centre are skipped) and the float pre-screens harmless.

Why Lemma T terminates where Lemma S cannot: at a germ box, `J₃↑`, `J₃↓` are the whole unit pieces of the lines up
to `O(box)` at the two corners, and for `T` inside the tile `f(T) ≈ ρ(1 − T + η) + ρ(T − η) = ρ` — the two lines
share one threshold and a uniform density has no gap between the up-line part above `T` and the down-line part
below `T + u₀`.  (For points the analogous sum has a gap between consecutive pivots: `S32_SHIFT.md` §1 defect 2.)
The loss is `O(box)`, so a germ box closes as soon as the cover has margin.

> **Lemma L (linear chord ends, axis lines).**  Let `ℓ` be a vertical line `x = ξ` (`t = y`, base `b = c_y`,
> `a = ξ − c_x`) or a horizontal line `y = η` (`t = x`, base `b = c_x`, `a = η − c_y`).  With `S = 2u`,
> `C = 1 − u²`, `N = 1 + u²`, inequality `k` at the point of parameter `t` is, at a pose with `u > 0`,
> `t ≤ t_k` ("up") or `t ≥ t_k` ("lo"), `t_k = b + num_k(u)/den_k(u)`:
>
> | line | `k` | type | `num_k` | `den_k` |
> |---|---|---|---|---|
> | V | 0 | up | `(½ − a) + (½ + a)u²` | `S` |
> | V | 1 | lo | `(−½ − a) + (a − ½)u²` | `S` |
> | V | 2 | up | `½ + 2a u + ½u²` | `C` |
> | V | 3 | lo | `−½ + 2a u − ½u²` | `C` |
> | H | 0 | up | `½ − 2a u + ½u²` | `C` |
> | H | 1 | lo | `−½ − 2a u − ½u²` | `C` |
> | H | 2 | lo | `(a − ½) + (−a − ½)u²` | `S` |
> | H | 3 | up | `(a + ½) + (½ − a)u²` | `S` |
>
> (the `C` rows hold for `u = 0` too).  Call `k` *typed* on the box if `den_k = C`, or `den_k = S` and `u₀ > 0`.
> Suppose the Lemma S intervals are `I_k = (−∞, α_k]` for the typed up `k` and `[β_j, ∞)` for the typed lo `j`;
> put `a↑ = min α_k`, `b↓ = max β_j`, and assume `b↓ < a↑`.  Let `Δ > 0`, and suppose every untyped `I_k`
> contains `[b↓ − Δ, a↑ + Δ]`.  Let `G(x)` be the mass of the line's segments in `(−∞, x]`, `ρ↑` the minimum of
> the (summed) density on `[a↑, a↑ + Δ]`, `ρ↓` that on `[b↓ − Δ, b↓]`.  Then at every admissible pose `P`:
>
>     μ_ℓ(Q) ≥ [G(a↑) − G(b↓)] + ρ↑ · min( min_{k up} t_k(P) − a↑ , Δ ) + ρ↓ · min( b↓ − max_{j lo} t_j(P) , Δ ).

*Proof.*  The table: e.g. V, `k = 0`: `X = (aC + (t−c_y)S)/N ≤ ½ ⇔ (t − c_y)S ≤ N/2 − aC`, and `S > 0`; the other
rows are the same computation (`X = (aC + bS)/N`, `Y = (−aS + bC)/N` with `a, b` the offsets of the point).
`α_k ∈ I_k ⊆ C_k` means the point `t = α_k` satisfies inequality `k` at every admissible pose, i.e.
`α_k ≤ t_k(P)`; so `r↑ := min_{up} t_k(P) ≥ a↑` and `r↓ := max_{lo} t_j(P) ≤ b↓`.  At `P` the typed inequalities
hold exactly on `[r↓, r↑]` and the untyped ones hold on `[b↓ − Δ, a↑ + Δ]`, so
`Q ∩ ℓ ⊇ Z := [max(r↓, b↓ − Δ), min(r↑, a↑ + Δ)] ⊇ [b↓, a↑] ≠ ∅`, and
`μ_ℓ(Q) ≥ G(min(r↑, a↑+Δ)) − G(max(r↓, b↓−Δ)) = [G(a↑) − G(b↓)] + [G(min(r↑,a↑+Δ)) − G(a↑)] + [G(b↓) −
G(max(r↓, b↓−Δ))]`; the second bracket is `≥ ρ↑·(min(r↑, a↑+Δ) − a↑)` and the third `≥ ρ↓·(b↓ − max(r↓, b↓−Δ))`
because `G` increases at rate `≥ ρ` on those intervals. ∎

> **Lemma L′ (short chord).**  Same setting, but without `b↓ < a↑` (e.g. a vertex of `Q` crossing the line, where
> the certified core `[b↓, a↑]` is empty).  Put `A = min(a↑, b↓) − Δ`, `B = max(a↑, b↓) + Δ`, suppose the untyped
> `I_k` contain `[A, B]`, and let `ρ` be the minimum density on `[A, B]`.  Then at every admissible pose
>
>     μ_ℓ(Q) ≥ ρ · min( min_{k up} t_k(P), B ) − ρ · max( max_{j lo} t_j(P), A ).

*Proof.*  The typed inequalities hold exactly on `[r↓, r↑]`, the untyped ones on `[A, B]`, so `Q ∩ ℓ ⊇
[max(r↓, A), min(r↑, B)]`.  If this interval is non-empty it lies in `[A, B]`, where `G` grows at rate `≥ ρ`, so
`μ_ℓ(Q) ≥ ρ(min(r↑,B) − max(r↓,A))`; otherwise the right-hand side is `< 0 ≤ μ_ℓ(Q)`. ∎  The right-hand side is again
a (min of affine) − (max of affine) in `c` at fixed `u`, so it enters Corollary L unchanged.  Since it can be
negative, the joint bound is computed with and without the Lemma L′ lines and the larger is used (without them,
those lines contribute their core, here 0).

> **Corollary L (joint bound).**  For a set of lines satisfying Lemma L or L′, let `Φ(c,u)` be the sum of their
> right-hand sides.  (i) For fixed `u`, every `t_k` is affine in `c`, so `Φ(·, u)` (a sum of minima of affine
> functions with coefficients `ρ ≥ 0`, plus constants) is concave in `c`.  At fixed `u` the admissible centres of
> the box form the rectangle `[max(cx₀, w/2), min(cx₁, m − w/2)] × [same in y]` (`w = (C + S)/N`); replacing
> `max(cx₀, w/2)` by **either** `cx₀` **or** `w/2` (and `min(cx₁, m − w/2)` by either term) gives a rectangle that
> *contains* it, so Φ's minimum over admissible centres is `≥` its minimum over the 4 corners of such a rectangle,
> for any choice.  The code splits the bin at (rational approximations of) the `u` where `w(u)/2` crosses a box
> side and on each sub-bin takes the tighter term; all corner coordinates are rational functions of `u`.
> (ii) At a corner each minimum is over rational functions of `u`; take its option `i*` smallest at mid-bin and an
> exact upper bound `σ ≥ max_u (A_{i*} − A_i)` over the others, so `min_i A_i ≥ A_{i*} − σ`.  (iii) The sum of the
> chosen options is `P(u)/Den(u)`, `Den = S^a C^b N^c` (`S` only when `u₀ > 0`).  If all degree-`n` Bernstein
> coefficients `β_i(Den)` on the sub-bin are `> 0` (`n = max(deg P, deg Den)`), then
> `min_u P/Den ≥ min_i β_i(P)/β_i(Den)` (and `max ≤ max_i`), because `P − λ Den = Σ_i B_i(u)(β_i(P) − λ β_i(Den)) ≥ 0`
> for `λ = min_i β_i(P)/β_i(Den)` (convex-hull property, `B_i ≥ 0`, `Σ B_i = 1`).  The minimum over sub-bins and
> corners of (iii) minus the slacks is a lower bound on the lines' mass at every admissible pose of the box; if any
> `β_i(Den) ≤ 0` the corollary is not used for that box. ∎

Lemma L is what makes the checker first-order exact off the germs: as the square moves, the chord ends of every
line move linearly and the gains at one end offset the losses at the other **inside one exact expression** (the
corner/Bernstein bound is exact at the corners and has `O(h²)` slack in `u`), instead of both ends being charged
`O(box)` as in Lemma S.  It is the continuous-density counterpart of zeromargin's disjunctive `CHAIN`.  Measured
effect in §4: the `2 %`-margin cell drops from 31,278 boxes to 320, and the box count becomes logarithmic in the
margin.

**Assembling `L`.**  Groups (Lemma T), and per line the core (Lemma S) and the Lemma L data.  Two assemblies,
both valid, the larger is used: (A) `Σ_groups max(T-bound, core) +` joint Lemma L bound (or cores, whichever is
larger) over the non-group lines; (B) the joint Lemma L bound over all eligible lines (group lines included) +
cores of the rest.  Plus polygons (Lemma S).  Distinct summands always involve disjoint sets of pieces.

> **Lemma P (points + pieces).**  Let `L ≥ 0` be a lower bound on the piece mass at every admissible pose of a box.
> Add to the point set a phantom point `φ` of weight `L'` (`L` rounded **down** to a multiple of zeromargin's common
> weight denominator), and pass it in zeromargin's `inh` mask.  If `ADM`, `P1`, `MIX` or `CHAIN` certifies the box,
> then `μ(Q) ≥ 1` at every admissible pose of the box.

*Proof.*  Every one of these primitives (`zeromargin.py`, `RUNG2.md` §3, §6–7) proves: for every admissible pose of
the box, the weight of the witness points lying in `Q` is `≥ 1`, where a point flagged in `inh` is only ever used as
"lies in `Q` at every admissible pose of the box" (the per-subtree cache; `cert_adm`, `cert_mix` and the `T` set of
`cert_chain` add it to the witness set with no geometric test).  Replace, in that statement, "`φ ∈ Q` with weight
`L'`" by "piece mass `≥ L ≥ L'`": for every admissible pose `μ_points(Q) + μ_pieces(Q) ≥ (witness weight without φ)
+ L' ≥ 1`.  `φ` sits at `(−1000, −1000)`, so no geometric test (`ADM`/`P1` masks, `CHAIN` reach and swing
candidates) ever selects it by itself; it is `inh`-only by construction.  The code keeps zeromargin's three weight
arrays (`W` Fractions, `Wf` floats for pre-screens, `Wnum` integers over `Wden` for `CHAIN`) consistent: `Wnum[φ]
= ⌊L·Wden⌋`, `W[φ] = Wnum[φ]/Wden`, `Wf[φ] = float(W[φ])` (floats only reject). ∎

**Symmetry.**  `--d4` (roots `[0,m/2]² × u ∈ [0,½]`) and the default (x → m−x and y → m−y; roots `c_y ≤ m/2`,
`u ∈ [0,½]`) use the arguments of `ZEROMARGIN.md` §2 / `S32_EXACT.md` §7, which need the **measure** invariant under
the generators.  Checked sufficiently: the multisets of (point, weight), (unordered segment endpoints, weight) and
(polygon vertex set, weight) are mapped onto themselves (integer coordinates).  `--full` (no symmetry): the cover on
`u ∈ [0,½]` plus its image under `σ:(x,y) ↦ (y,x)` on `u ∈ [0,½]`, all centres; this covers every pose because
`σQ(c,θ) = Q(σc, −θ) = Q(σc, 90° − θ)` and `μ(Q) = (σ_*μ)(σQ)`, so `θ ∈ [45°, 90°]` for the cover is `θ ∈ [0°,45°]`
for its image (this differs from zeromargin's `--full`, which uses `u ∈ [0,1]`: Lemma T is formulated at `θ = 0`).

**Validation.**  `mixed_cover.validate` (all integers, pieces in the container, `w ≥ 0`, no zero-length segment,
polygons with left turns and positive area) plus `check_simple_convex` here: a polygon passing `validate` could wind
twice; it is rejected unless its strictly convex vertices are exactly its convex hull in CCW order.

**What is exact, what is float.**  Everything that certifies is `Fraction`: Bernstein coefficients, intervals,
clipping, the breakpoint minimisation, `L`, and all of zeromargin's decisions.  Floats: the reach filter (drops
pieces only) and zeromargin's own pre-screens (reject only).

## 3. What the checker does at a tile germ

A germ box touches a pose `(ξ + ½, η + ½, θ = 0)` with all four edges of `Q` on grid lines.  There:

* the vertical pair `x = ξ, ξ + 1` and the horizontal pair `y = η, η + 1` each form a Lemma T group; `J₃` of each
  line is its unit piece minus `O(box)` at the two corners (Lemma S for the three non-singular inequalities), and the
  group bound is `≈ ρ(1 − O(box))` whatever the threshold does — the `θ → 0+` "pivot" poses, where one line's piece
  above `T` and the partner's piece below `T + u₀` are in, included;
* the other lines and polygons get Lemma S / Lemma L; the points go to zeromargin's primitives through the phantom;
* if that fails and the cover has **points on the germ lines**, a second attempt moves those points into the groups
  (Corollary T′), which replaces zeromargin's two-chain product by two one-dimensional minimisations.

Measured (§4): on the toy grid cover the four germ roots of a tile close in **12 boxes** (depth 1, < 1 s); with
Lemma T switched off the same roots leave **1,340 boxes uncertified at depth 16** (the `RUNG2.md` Theorem 1
obstruction: a fixed witness set at the germ holds about half of the four lines).  On agent A's `m = 4` mixed cover
(`×1.05`) the germ cell `[1.4,1.6]²` left 52 boxes at depth 18 with the first attempt alone (float minimum there
`1.037`: a checker limit, not a hole) and closes in **566 boxes, depth 6** with Corollary T′ (109 leaves by it).

## 4. Tests

All runs: `taskset -c 10-14,26-30`; files in `runs/zmm/` (not committed); toy covers written by
`python3 search/zm_mixed_test.py toys runs/zmm`.  CPU = wall × processes.

### 4.1 Point-only files agree with `zeromargin.py cert` (certified)

`certificates/rung2/s13_closed_cover_4.txt`, `--disj --chain-from 0 --depth 18`, default D2 roots:

| column | `zeromargin.py` | `zm_mixed.py` |
|---|---|---|
| `c_x ∈ [0.9, 1.0]` (160 roots) | 182 boxes, depth 2, ADM 86 / CHAIN 39 / EMPTY 46 | identical |
| `c_x ∈ [1.4, 1.5]` (160 roots, tile germs `(1.5, ½)`, `(1.5, 1.5)`) | 910 boxes, depth 13, ADM 75 / CHAIN 411 / EMPTY 49 | identical (also re-run with the final code) |

Same boxes, same leaf census: with no pieces the phantom has weight 0 and every decision is zeromargin's.  A full
run was not repeated (the shipped one is 16,872 boxes; `M4_MARGIN.md` §2 puts a comparable full run at ~27 CPU-h).

### 4.2 Randomised tests of the lemmas (`zm_mixed_test.py selftest`)

Every check compares a certified bound with the **exact** (Fraction) mass at a random admissible rational pose of
the box (faces, corners, tiny `u` and the Lemma T pivot poses `c_x = ξ + 1/(2 cos θ)` favoured):

```
$ python3 search/zm_mixed_test.py selftest --n 300 --seed 31
[1] exact vs mixed_cover float mass, 400 poses: max |diff| = 3.953e-09  ok
[2] Lemma S (lines): 48768 (box, line, t, pose) checks, 0 violations  ok
[3] Lemma S (polygons): 496 (box, pose) checks of bound <= mass, 0 violations  ok
[4] Lemma T at tile germs: 4225 (box, pose) checks of bound <= exact mass, 0 violations, Lemma T raised the bound in 149/300 boxes  ok
[4c] Lemma T with line points: 4235 checks, 0 violations (12751 points moved)  ok
[4b] Lemma L, dense axis covers: 1730 (box, pose) checks, 0 violations, Lemma L used in 35/300 boxes  ok
[4d] Lemma L on near-tight grid covers: 1088 checks, 0 violations, Lemma L used in 135/300 boxes  ok
[5] PIECE bound, random covers/boxes: 1200 checks, 0 violations  ok
selftest: PASS in 166s
```

[1] also checks `mixed_cover.mass_in_square_float` (the float evaluator agrees with the exact one to `4·10⁻⁹` at
generic poses).  [3] uses boxes on which the polygon bound is positive (a separate run: 780 admissible poses in 300
random boxes, bound `> 0` in all 300, 0 violations).  [4d] is the tight case: jittered T1 densities at margin
`≈ 0.25 %`, boxes down to `1/1280`, where Lemma L (and its short-chord form) decides most boxes.

### 4.3 Toy covers that must certify (certified)

`T1`: the 6 interior grid lines of `[0,4]²` at uniform density `ρ = 5/8` (240 segments of length `1/10`, total 15).
Its minimum is `2 · 0.8284 ρ = 1.0355` at `θ = 45°` (a flat family of diamonds crossing two parallel lines, including
the corner square), every germ has `2ρ = 1.25` one-sided.  `T1m*`: the same with `ρ` lowered so that the minimum is
`1 + ε`.  All `--d4` unless a cell is named.

| cover | margin `ε` | run | boxes | depth | CPU | result |
|---|---|---|---|---|---|---|
| T1 | 3.55 % | Lemmas S + T only (before Lemma L) | 194,850 | 11 | 626 s | VERIFIED-D4 |
| T1 | 3.55 % | final | **17,626** | 9 | 304 s | VERIFIED-D4 |
| T1m2 | 2.0 % | cell `[1.4,1.5]²`, S + T only | 31,278 | 13 | 195 s | VERIFIED (cell) |
| T1m2 / m1 / m05 / m025 | 2 / 1 / 0.5 / 0.25 % | cell `[1.4,1.5]²`, final | 320 / 350 / 376 / 394 | 6–7 | 5 s each | VERIFIED (cell) |
| **T1m025** | **0.2513 %** | final | **42,300** | 13 | 884 s | **VERIFIED-D4** |

With Lemma S alone the box count grows like `ε^{−2.2}` (3.5 % → 2 %: ×3.4); with Lemma L it is logarithmic in `ε`
down to 0.25 %.  The stress test of the T1m025 leaves (207,866 sampled poses) finds minimum `1.0025128`, the designed
minimum, and no exact violation.

### 4.4 Rejection tests (certified that nothing wrong is certified; heuristic hole location)

Hole = float scan minimum (`zm_mixed_test.py holes`), then **exact** mass at the snapped rational pose.  A test
passes when the hole pose lies in an uncertified box and in no certified leaf (`zm_mixed_test.py locate`), and the
leaf stress test (`zm_mixed.py stress`, exact re-check of every sampled pose with float mass `< 1 + 10⁻⁶`) finds no
certified pose below 1.

| cover | construction | exact hole | run | uncertified | hole in UNCERT | stress of certified leaves |
|---|---|---|---|---|---|---|
| R_lighten | T1 × 0.95 | `0.984108` at `(0.7067, 0.7067, 43°)` | cell `[0.7,0.8]²`, depth 11 | 267 | yes | min 1.000064, 0 violations |
| R_remove | T1 minus the D4 images of one `1/10` segment (`x = 2`, `y ∈ [1.2,1.3]`) | `0.973034` at `(1.4971, 1.3971, 45°)` | cell `[1.4,1.5]×[1.3,1.4]`, depth 11 | 7,614 | yes | min **1.000000** (exact equality at a box corner: Lemma L is tight there), 0 violations |
| R_pivot | `ρ = 0.81`, pieces of the tile `[1,2]²` edges removed (`x = 1, y ∈ [1.5,1.9]`; `x = 2, y ∈ [1.1,1.5]`, D4 images) | `0.564541` at `(1.4972, 1.5472, 7°)` | cell, depth 10 | 8,128 | yes | 0 violations |
| **R_shift** | the same pieces **shifted off the grid lines** by `1/100` (`x = 1 → 0.99`; `x = 2 → 2.01` and `1.99` at half weight) | `0.813128` at the **pivot** `(1.500025, 1.49, θ = 0.573°)` | full D4, depth 14 | 268 (all next to the shifted pieces) | yes | min 1.0019, 0 violations |
| T1m0 | margin `−6·10⁻⁶` | at 45° | cell `[1.4,1.5]²`, depth 14 | 2,048, all at `θ ∈ [44.60°, 45.37°]` | yes | — |
| T1mneg | margin `−0.09 %` | at 45° | cell, depth 12 | 4,608, all at `θ ∈ [41.5°, 48.4°]` | yes | — |

R_shift is the germ-specific test: its hole exists only at `θ > 0` on the pivot curve `c_x ≈ 1.5 + θ²/4`, where the
threshold `T` of the pair `x = 1 / x = 2` sits at `y = 1.5` and both displaced pieces are outside `Q`; at `θ = 0` the
same squares are fine.  Earlier runs of the same tests with the pre-Lemma-L code (full D4, depth 9) and with the
intermediate versions also rejected at the hole poses, with 0 stress violations.

### 4.5 Agent A's `m = 4` mixed covers (certified for the scaled files; the scalings are mine)

| file | points / segments | total | run | boxes | depth | CPU | result |
|---|---|---|---|---|---|---|---|
| `lc_m4mixA_r5.txt × 1.05` | 3,461 / 1,016 (on `x, y ∈ {1,2,3}`, `q = 50`) | 12.929248 | D4, `--disj`, depth 18, final code | **27,648** | 14 | 1.6 h | **VERIFIED-D4** |
| `lc_m4mixA_r7.txt × 1.02` | 3,001 / 912 | 12.566207 | same | 116,176 | 18 | 4.7 h | NOT VERIFIED: 74 boxes |

The r5 file is thus an exactly certified mixed cover of `[0,4]²` with total `12.929 < 13` (another `s(13) = 4`
certificate; heavier than the lightest point cover, `12.732`, `M4_MARGIN.md`); its leaf stress test (136,948 poses)
has minimum `1.044`, 0 violations.  **r7 × 1.02**: all 74 uncertified boxes (depth 18, side `1/640`, at
`θ ≈ 32°` near `(1.24, 1.56)` and `(1.39, 1.31)`) have sampled minimum `≥ 1.0154` (`zm_mixed_test.py unc`): the
cover has margin there, the checker is short.  Diagnosis of one of them: the pieces carry `0.651–0.676` and the
points `0.356–0.378` over the box, **anti-correlated** (total minimum `1.018`, but `min + min = 1.007`); the piece
bound is `0.645` and the point primitives cannot find the remaining `0.355`.  This is the decoupling of Lemma P
(one uniform `L` per box, while `CHAIN` splits the box into regions): see §5.

### 4.6 Performance against the point checker (comparable covers)

| cover | checker | boxes | depth | CPU | CPU / box |
|---|---|---|---|---|---|
| T1 (lines, total 15, margin 3.55 %) | `zm_mixed.py` | 17,626 | 9 | 304 s | 17 ms |
| P_T1pts (the same lines as 2,406 points at spacing `1/100`, margin 2.5 %) | `zeromargin.py --disj` | 19,998 | 12 | 606 s | 30 ms |
| germ cell `[1.4,1.6]²`, T1 | `zm_mixed.py` | 1,040 | 6 | 19 s | |
| germ cell `[1.4,1.6]²`, P_T1pts | `zeromargin.py --disj` | 1,932 | 7 | 86 s | |
| the 4 germ roots (u-bin 0) of that cell, T1 | `zm_mixed.py` | 12 | 1 | < 1 s | |
| T1m025 (margin 0.25 %) | `zm_mixed.py` | 42,300 | 13 | 884 s | 21 ms |

For comparison, the `s(32)` point checks spend 2–4.4 CPU-h **per germ root** (`S32_EXACT.md` §2).  For a mixed cover
the cost is dominated by the points (zeromargin's primitives); the pieces cost ~10–20 ms per box in `Fraction`
arithmetic.

## 5. Known limitations

* **Axis-parallel segments only get Lemmas T and L.**  Segments in other directions and polygons get Lemma S (loss
  `O(box)`), which is sound but slow at small margins (for T1, Lemma S alone needed `ε^{−2.2}` boxes).  Parallel
  non-axis segments at distance exactly 1 would form tilted germs (at `θ` = their direction) that nothing here closes.
* Lemma T is formulated at `θ = 0`; with `u ≤ ½` roots that is the only germ.  Hence `--full` = cover + swapped cover
  on `u ∈ [0,½]` (§2), not zeromargin's `u ∈ [0,1]`.
* **Pieces and points are coupled only through one number per box** (Lemma P).  When a tight pose has a discrete
  point entering `Q` while the pieces' mass falls (r7 × 1.02, §4.5), `min(points) + min(pieces)` over a box stays
  below the joint minimum until the box no longer contains the jump, which is slow.  The fix is to hand `CHAIN` a
  region-wise piece bound (Lemma L is an explicit concave function of the pose, so it can be evaluated per
  `CHAIN` region); not built.
* Corollary T′ (points in the groups) is a second attempt per box, only for points exactly on the germ lines.
* The float reach filter and zeromargin's float pre-screens can only drop certifications (sound, possibly
  incomplete).  `L` is rounded down to zeromargin's weight denominator (`≤ 1/W` loss).
* The checker is Python + `Fraction`; there is no second independent implementation (zeromargin's primitives are
  independent code, and the lemmas are covered by the randomised exact tests of §4.2).  The Lean reduction
  `packing_le_weight` still needs generalising from point sets to these measures (FORMAT.md; not part of this task).

## 6. Reproduce

```
S=search; R=runs/zmm
python3 $S/zm_mixed_test.py selftest --n 300 --seed 31
python3 $S/zm_mixed_test.py toys $R
python3 $S/zm_mixed.py cert $R/zmm_T1m025.txt --d4 --depth 22 --nproc 4 --dump $R/leaves_T1m025.txt
python3 $S/zm_mixed.py stress $R/zmm_T1m025.txt $R/leaves_T1m025.txt --per 10
python3 $S/zm_mixed.py cert $R/zmm_R_shift.txt --d4 --depth 14 --dump $R/leaves3_R_shift.txt
python3 $S/zm_mixed_test.py locate $R/leaves3_R_shift.txt 359991/239990 149/100 600/119999
python3 $S/zm_mixed.py cert certificates/rung2/s13_closed_cover_4.txt --disj --chain-from 0 --depth 18 --cx-lo 1.4 --cx-hi 1.5
python3 $S/zm_mixed.py pose FILE --cx C --cy C --u U          # exact mu(Q) at one rational pose
python3 $S/zm_mixed_test.py unc FILE LEAFDUMP                  # float diagnosis of uncertified boxes
```

