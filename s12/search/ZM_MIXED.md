# Exact checker for mixed covers (points + segments + polygons): `zm_mixed.py` (2026-09-27)

Brief: `tasks/line-cover/B.md`.  Code: `search/zm_mixed.py` (checker), `search/zm_mixed_test.py` (tests, toy covers,
leaf stress).  Format: `tasks/line-cover/FORMAT.md` (public copy: `certificates/s21/FORMAT.md`); reader `search/mixed_cover.py` (agent A).  `zeromargin.py` is
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
  hole is a germ pivot pose at `θ = 0.57°`.  The evidence that matters for the certificate is the audit's
  (`ZM_MIXED_AUDIT.md` §3): **component tests** (every certified component bound — piece bound, phantom, SPLIT
  region bound, points — against the exact mass of that component at adversarial rational poses, ≈ 115 M exact
  checks on the real `m = 5` cover, 0 violations, several at slack exactly 0) and **rejection runs on three holed
  versions of the real `m = 5` cover** (tilted dip, SPLIT-heavy cell, germ-pivot hole), each refused with its exact
  hole inside an uncertified box.  (The leaf stress tests below re-check only poses whose *total* float mass is
  `< 1 + 10⁻⁶`; on a cover with margin that set is empty, so they are no evidence about the checker there: §8.)
* **Performance**: comparable line vs point covers (§4.6): 17,626 boxes / 304 CPU-s (lines) vs 19,998 / 606 CPU-s
  (the same lines as 2,406 points, `zeromargin.py`); germ cell 1,040 boxes / 19 s vs 1,932 / 86 s.
* **`m = 5` (§7.5): agent A's candidate `× 1.003`, total `20.89474919732 < 21`, is VERIFIED-D4** (40,000 roots,
  461,204 boxes, 19.7 CPU-h, 0 uncertified).  Since then (§8): re-run in **certificate mode** with full provenance
  (the shipped run, `certificates/s21/zm_mixed_d4/`), re-verified by the independent checker `zmx2`
  (`ZMX2.md`), and the measure version of the FORMAT.md reduction is proved in Lean (`notes/lean-s21.md`):
  **`s(21) = 5`**, bundled in `certificates/s21/`.  The unscaled candidate (20.832) left one root (float margin `0.33 %` there).
* **Agent A's `m = 4` mixed covers** (scaled by me, checker tests only): `r5 × 1.05` (total 12.929)
  **VERIFIED-D4**.  `r7 × 1.02` (12.566) failed in round 1 (74 boxes); **round 2 (§7)** adds the region-wise
  primitive `SPLIT` (Lemma R) and affine end minorants in Lemma L, and it **verifies (D4, 43,970 boxes, 4.0 CPU-h)**.
  At the binding cells the checker verifies down to `0.05 %` over the true minimum and rejects at `−0.0014 %` (§7).


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

**Lemma L, general ends (what the code uses).**  In Lemma L replace `ρ↑·min(r↑ − a↑, Δ↑)` by
`min(ℓ↑(r↑ − a↑), g↑(Δ↑))`, where `g↑(x) = G(a↑ + x) − G(a↑)` is the gain of the up end, `Δ↑ > 0` any number (with
the untyped inequalities holding on `[b↓ − Δ↓, a↑ + Δ↑]`), and `ℓ↑(x) = s x + i` any affine function with `s ≥ 0`
and `ℓ↑ ≤ g↑` on `[0, Δ↑]` (the intercept `i` may be negative); likewise at the lo end.  *Proof.*  With `x =
min(r↑ − a↑, Δ↑) ∈ [0, Δ↑]` the up gain is exactly `g↑(x)` (the proof of Lemma L).  If `r↑ − a↑ ≤ Δ↑`,
`min(ℓ↑(r↑ − a↑), g↑(Δ↑)) ≤ ℓ↑(x) ≤ g↑(x)`; otherwise `x = Δ↑` and the min is `≤ g↑(Δ↑) = g↑(x)`. ∎  The term is a
minimum of affine functions of `r↑ = min_k t_k` with non-negative slope, so Corollary L applies unchanged.  The code
takes `ℓ` = the edge of the lower convex hull of `g`'s vertices on `[0, Δ]` (the greatest convex minorant of the
piecewise-linear `g`, so every hull edge is `≤ g`) that contains a float guess of where the chord end sits (any
choice is sound), and `Δ` = an exact upper bound of how far the end can move over the box (`max` over sub-bins and
corners of Bernstein-ratio bounds of `t_k − a↑`, which is affine in the centre), capped by the crude
`2(dx + dy) + 8 du`.  The earlier `ρ_min` form is the special case `ℓ = ρ_min x`.  This matters for LP covers,
whose segment densities jump by factors of 10 between neighbouring `1/q` pieces (e.g. `0.10 → 3.58` at `x = 0.98` on
`y = 1` in r7): the cone through the origin sees the low one, the hull edge the real slope.

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

**Lemma L′, hull form (what the code uses).**  Fix `p ∈ [A, B]` and let `g(x) = G(p + x) − G(p)` on
`[A − p, B − p]` (increasing, piecewise linear, `g(0) = 0`), `ℓ₁ ≤ g` an affine minorant and `ℓ₂ ≥ g` an affine
majorant there, both with slope `≥ 0` (edges of the lower / upper convex hull of `g`'s vertices).  Then at every
admissible pose `μ_ℓ(Q) ≥ ℓ₁(min(r↑, B) − p) − ℓ₂(max(r↓, A) − p)`.  *Proof.*  `x = min(r↑, B) ≥ a↑ ≥ A` and
`y = max(r↓, A) ≤ b↓ ≤ B`, so both lie in `[A, B]`.  If `x ≥ y`, `Q ∩ ℓ ⊇ [y, x]` (typed inequalities hold on
`[r↓, r↑]`, untyped ones on `[A, B]`) and `μ_ℓ(Q) ≥ G(x) − G(y) = g(x − p) − g(y − p) ≥ ℓ₁(x − p) − ℓ₂(y − p)`; if
`x < y` the right-hand side is `≤ g(x − p) − g(y − p) < 0 ≤ μ_ℓ(Q)`. ∎  Both terms are minima of affine functions of the
`t_k` with the right signs (`ℓ₁, ℓ₂` increasing), so Corollary L applies.  (`p` = the float chord centre, the hull
edges at the float chord ends; any choice is sound.)

> **Lemma V (splitting on a short chord).**  For a Lemma L′ line with active inequalities `k` (up) and `j` (lo),
> `D = t_k − t_j` is, at fixed `u`, an affine function of the centre (the base coordinate cancels; the coefficient of
> the other coordinate is `0` or `±N²/(SC)`).  The admissible poses split into `R_a = {D ≤ 0}` and `R_b = {D ≥ 0}`;
> on `R_a` the line is bounded below by `0`, on `R_b` by its Lemma L′ term; with the other lines' Lemma L terms, the
> region minima (computed as in Lemma R, with the constraint `D`) give `min(Φ_a, Φ_b)`, a valid lower bound on the
> lines' mass over the box.

*Proof.*  Every admissible pose lies in `R_a ∪ R_b`; on each region the claimed bound holds pointwise (mass `≥ 0`,
respectively Lemma L′), and region_phi's candidate set contains every vertex of the region polygon at every `u`
(Lemma R's argument, which only needs the constraint to be affine in the centre with coefficients that are
rational functions of `u`: the edge intersections are then rational functions too, computed exactly after putting
the coefficient over one denominator `S^a C^b N^e`). ∎

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

> **Lemma R (region-wise coupling of pieces and points; primitive `SPLIT`).**  Fix a box with `u₀ > 0`.  Let `T`
> be a set of points in `Q` at every admissible pose of the box (`ADM`/`P1`/inherited).  Let `q₁, …, q_k` be
> *swing* points of one kind `κ` (each satisfies its three other inequalities at every admissible pose of the box)
> with `G_{q₁} ≤ G_{q₂} ≤ … ≤ G_{q_k}` on the whole box, where `G_q = α_q(u) + β(u) c_x + γ(u) c_y` is zeromargin's
> violation polynomial of kind `κ` (`G ≤ 0` ⇔ inequality `κ`; `β, γ ∈ {±2C, ±2S}` depend on `κ` only).  For
> `r = 0..k` put `R_r = {admissible poses : G_{q_r} ≤ 0 (if r ≥ 1), G_{q_{r+1}} ≥ 0 (if r < k)}`,
> `D_r = {p swing : G_p ≤ G_{q_r} on the box}`, `U_r = {p swing : G_p + λ G_{q_{r+1}} ≤ 0 on the box, some λ > 0}`.
> Let the piece bound of the box be `rest + max(core, Φ)` as assembled in §2 (Φ the Corollary L term of a set of
> Lemma L/L′ lines, `rest` the part carried by the other pieces), and let `Φ_r` be a lower bound of the Lemma L
> right-hand side over `R_r`.  If for every `r`
>
>     w(T) + w(D_r ∪ U_r) + rest + max(core, Φ_r) ≥ 1        (or R_r = ∅),
>
> then `μ(Q) ≥ 1` at every admissible pose of the box.

*Proof.*  Take an admissible pose and let `r = max{j : G_{q_j} ≤ 0}` (0 if none).  Monotonicity gives
`G_{q_{r+1}} > 0` (if `r < k`), so the pose is in `R_r`.  `T ⊆ Q`.  For `p ∈ D_r`: `G_p ≤ G_{q_r} ≤ 0`, and `p`
satisfies its other three inequalities, so `p ∈ Q` (the chain members `q_j`, `j ≤ r`, are in `D_r`).  For `p ∈ U_r`:
`G_p ≤ −λ G_{q_{r+1}} < 0`, so `p ∈ Q`.  The pieces carry `≥ rest + core` (Lemma S/T parts, valid at every pose of the
box) and `≥ rest + Φ_r` (Lemma L at this pose, `≥` its minimum over `R_r`).  Points and pieces are disjoint, so
`μ(Q) ≥ w(T) + w(D_r ∪ U_r) + rest + max(core, Φ_r) ≥ 1`. ∎

**Computing `Φ_r`.**  At fixed `u`, the admissible centres of `R_r` lie in `rect(u) ∩ {G_{q_r} ≤ 0} ∩ {G_{q_{r+1}} ≥ 0}`
(`rect(u)` any rectangle containing the admissible centre rectangle, as in Corollary L), a convex polygon; its two
constraint lines are **parallel** (same `β, γ`), so every vertex is a rectangle corner or a constraint line ∩ an edge
line.  These candidates are rational functions of `u` with denominators `S^a C^b N^c` (the edge intersection
divides by `β` or `γ`, i.e. by `C` or `S`; `S` needs `u₀ > 0`).  The Lemma L bound is concave in the centre
(Corollary L (i)), so its minimum over the polygon is `≥` its minimum over **any finite set whose convex hull
contains the polygon**; the code keeps every candidate that is not *proven* (exact Bernstein-ratio bounds over the
sub-bin) to violate a constraint or to lie outside its edge — every true vertex at every `u` survives, since a vertex
satisfies all constraints at its own `u` — and bounds the Lemma L expression at each kept candidate exactly as in
Corollary L (ii)–(iii).  If no candidate survives on a sub-bin, the polygon is empty for every `u` of it (a non-empty
compact convex polygon has a vertex).  If a Bernstein bound is unavailable, `Φ_r` is not used (the region keeps
`max(core, ·) = core`).

`SPLIT` is tried after `ADM → P1 → MIX → CHAIN` (and Corollary T′) fail, one chain per swing kind (heaviest kind
first, chain built greedily in the order of `G` at the box centre, each link checked exactly with zeromargin's integer
`_gle0`; `D_r`, `U_r` by the same binary searches as zeromargin's `CHAIN`, with `λ ∈ {1, ½, 2}`).  Regions already
covered by the box bound (`w(T) + w(D_r ∪ U_r) + L_box ≥ 1`, valid by the same argument with `Φ_r` replaced by the
box's bound) skip `Φ_r`.  What it adds over the phantom (Lemma P): the pieces' bound is taken **on the region**, i.e.
on the side of the pivot surface where the point is out — exactly where the piece mass is larger when points and
pieces are anti-correlated (r7, §4.5).

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
bound is `0.645` and the point primitives cannot find the remaining `0.355`.  **Corrected in §7:** the loss was
mostly in the piece bound (Lemma L's `ρ_min` ends against LP densities that jump tenfold between neighbouring
pieces, plus lines anti-correlated with each other), with the point/piece coupling second; both are fixed there
and r7 × 1.02 verifies.

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
* **Pieces and points: region-wise only along one chain.**  `SPLIT` (Lemma R) couples them per region of a
  single-kind pivot chain; products of two chains (as `CHAIN` does for points) are not built, `u₀ = 0` boxes are
  excluded (the edge intersections divide by `S`), and only Lemma L/L′ lines get region bounds (Lemma T groups,
  Lemma S lines and polygons stay box-wide constants).  Sufficient for r7 × 1.02 and for 0.05 % margins at the
  binding cells (§7).
* Corollary T′ (points in the groups) is a second attempt per box, only for points exactly on the germ lines.
* The float reach filter and zeromargin's float pre-screens can only drop certifications (sound, possibly
  incomplete).  `L` is rounded down to zeromargin's weight denominator (`≤ 1/W` loss).
* The checker is Python + `Fraction`.  (When this was written there was no second implementation; there is now
  `zmx2`, `ZMX2.md`, written without reading this file, and the Lean reduction for measures,
  `notes/lean-s21.md`; §8.)

## 7. Round 2 (2026-09-27): region-wise coupling, sharper Lemma L, r7 × 1.02 verified

**What changed in the code** (all in `zm_mixed.py`; `zeromargin.py` untouched, sha unchanged):

1. **`SPLIT` (Lemma R, §2):** after `ADM/P1/MIX/CHAIN` and Corollary T′ fail, one pivot chain per swing kind is built
   and every region gets its own piece bound (the Lemma L expression minimised over the region's polygon, whose
   vertices are enumerated exactly), next to its own point witness set.
2. **Lemma L, general ends (§2):** the end gain is bounded by an edge of the lower convex hull of the actual
   mass function (negative intercept allowed) instead of `ρ_min·x`, and the cap `Δ` is an exact bound on how far the
   end can move over the box.
3. `--resume FILE` (one JSON line per finished root; restart skips them), per-root CPU in the census.

**Diagnosis that led to (2)** (one of r7's 74 boxes, `[1.2281,1.2297]×[1.5516,1.5531]`, `θ ≈ 34°`): the piece bound
was `0.6429` against a float minimum `0.6481`; per line, the float minima summed to `0.6426` — the lines are
anti-correlated with **each other** (the square slides mass from one line to the next), which only a joint bound sees;
and Lemma L's joint bound gained nothing because the relevant chord end of `y = 1` sits on a density jump
(`0.10 → 3.58` at `x = 0.98`), so `ρ_min = 0`.  The points (`0.356–0.380`) were a smaller part of the gap: a
single-kind chain leaves `≈ 0.352`.  After (1)+(2), the box bound in that area is `0.6466–0.6514` and the boxes
close at depth `≤ 20`.

**Tests added** (`zm_mixed_test.py selftest`, all exact at rational poses): `[6]` Lemma R — random mixed covers
(jittered grid lines + 20–80 points around the box), random boxes with `u₀ > 0`; for each sampled pose the region
`r` is computed from the exact `G` values and the test checks that **every point claimed for region `r` is in `Q`**
and that **region `r`'s piece bound ≤ the exact piece mass**, with half of the poses placed exactly on (or `10⁻⁹` off)
a pivot surface `G_q = 0`: 2,338 checks, 0 violations.  `[4e]` the piece bound on the real segments of A's r7 (880
checks, 0 violations; Lemma L decisive in 103/250 boxes).  All earlier tests still pass (§4.2, rerun).

### 7.1 r7 × 1.02 (certified)

| run | boxes | depth | CPU | leaves (ADM / CHAIN / SPLIT / other) | result |
|---|---|---|---|---|---|
| round 1 code, depth 18 | 116,176 | 18 | 4.7 h | 30,824 / 26,265 / — | 74 uncertified |
| + SPLIT only | 111,086 | 18 | 6.3 h | 27,483 / 25,043 / 2,074 | 55 uncertified |
| + SPLIT + general ends, depth 18 | 43,902 | 18 | 3.7 h | 9,410 / 8,706 / 3,221 | 18 uncertified (2 roots; both close at depth 20) |
| **same, depth 20** | **43,970** | 20 | **4.0 h** | 9,446 / 8,718 / 3,225 | **VERIFIED-D4** |

`lc_m4mixA_r7.txt × 1.02`, total `12.566207`.  Leaf stress test: 217,703 poses, minimum `1.0135`, 0 violations.

### 7.2 How tight: scale factor over the minimum (certified per cell)

r7's float minimum (A's strict protocol, `line_cover.py eval --pitch 0.002 --tiles`) is `μ_f = 0.9858395` at
`(1.49977, 3.48477, 0.028°)`; a targeted search found a lower **exact** value, `μ₀ ≤ 0.9853327085` at the admissible
rational pose `(45792065/91561246, 63046589/42031580, u = 17858/923433989)` (the wall germ `(½, 1.5)`, `θ = 0.002°`, its
D4 image).  Files `r7 × λ/μ_f` (exact rational factors), cells of 32 / 8 roots, `--depth 22`:

| cell | what | λ = 1.005 | 1.002 | 1.001 | 1.0005 | 1.0002 | 0.999 |
|---|---|---|---|---|---|---|---|
| `[0.4,0.6]×[1.4,1.6]` | the binding wall germ | — | ✓ 1,176 | ✓ 1,300 (d 20) | **✗ 12** | **✗ 13** | — |
| `[1.4,1.6]×[0.4,0.6]` | its θ ↔ 90° − θ image | ✓ 852 | ✓ 1,204 | ✓ 1,288 | ✓ 1,314 | ✓ 1,332 | ✗ 35 |
| `[1.2,1.3]×[1.5,1.6]` | the anti-correlated `θ ≈ 33°` cell | ✓ 3,472 | ✓ 4,680 | ✓ 5,162 | ✓ 5,530 | ✓ 5,744 | ✗ 421 at depth 13 (local min `≈ 1.009`: depth, not a hole) |
| `[1.4,1.6]²` | interior germ | ✓ 704 | ✓ 832 | ✓ 894 | ✓ 942 | ✓ 964 | ✓ 1,066 |

(✓ = VERIFIED, number = boxes.)  In terms of the exact minimum, `λ/μ_f · μ₀ = 1.000485` at `λ = 1.001` (**verified**)
and `0.9999857` at `λ = 1.0005` (the cover has a hole of depth `1.4·10⁻⁵`: **rejected**, and the hole pose lies in an
uncertified box and in no certified leaf; same at `1.0002`, hole `0.99969`).  So at the binding cells the checker
needs **≤ 0.05 % over the true minimum**, and it resolves a `1.4·10⁻⁵` hole.  (A's float minimum overestimates the
true one by `0.05 %` at this germ: the minimum sits at `θ = 0.002°`, below the scan's angles.)

### 7.3 Rejection tests for SPLIT (certified that nothing wrong is certified)

`r7 × 1.002` (raw factor), cell `[1.2,1.3]×[1.5,1.6]` (where SPLIT does most of the work), depth 16: hole found by
float search and confirmed exactly, `0.9983051` at `(2561819/2065521, 1032242/661481, u = 10502771/36579205)`
(`θ = 32.04°`).  Result: 14,394 boxes, 1,731 SPLIT leaves, 2,590 uncertified; the hole pose is in an uncertified
box and in no certified leaf; leaf stress test over the 4,611 certified leaves (184,443 poses): minimum `1.00019`,
0 violations.  The wall-germ rejections of §7.2 (λ = 1.0005, 1.0002) contain SPLIT leaves too (26–27 per run);
stress tests of those dumps: 0 violations.

### 7.4 Round-2 code, further changes found on the `m = 5` pilot

On the first `m = 5` pilot (depth 20) the tilted-dip cell `[0.5,0.6]×[1.4,1.5]` and the wall column left 81 boxes
(float minimum there `≥ 1.004`: checker limits), at wall poses with `θ = 2.5–7°` where the top vertex of `Q` crosses
`y = 2` next to the grid point `(1, 2)`: the chord on `y = 2` swings between `0.005` and `0.062` inside one box and the
Lemma L′ term (`ρ_min`) saw nothing.  Two more pieces (both in §2, with proofs, and in the tests `[4e]`, which now
also run on the `m = 5` candidate's segments and count Lemma V's uses): **Lemma L′ in hull form** (affine minorant /
majorant of the mass function instead of `ρ_min`), and **Lemma V** (split the box on whether the short chord exists).
With them, and `--depth 24`, the pilot is clean.  (The r7 runs of §7.1–7.3 used the code before these two changes.)

### 7.5 The `m = 5` candidate `runs/line-cover_m5_candidate.txt`

7,536 points + 1,872 segments, total `520806311/25000000 = 20.83225244 < 21`, D4-invariant (checked), agent A's
`lc_m5germC_r12` scaled by `1.004017` (`1.0025×` over its float-confirmed minimum, `search/LINE_COVER.md`).  Not scaled
by me.  `--d4 --disj --chain-from 0 --depth 24`, `zm_mixed.py` at commit `d4b97f9`
(sha256 in `runs/zmm/m5/full_sha.txt`).

Pilot (certified, partial):

| cell | roots | boxes | depth | CPU | result |
|---|---|---|---|---|---|
| interior germ `[1.4,1.6]²` | 32 | 5,470 | 17 | 755 s | VERIFIED |
| tilted dip `9.13°`: `[0.5,0.6]×[1.4,1.5]` | 8 | 4,472 | 24 | 341 s | VERIFIED |
| tilted dip `31.9°`: `[2.3,2.4]×[1.6,1.7]` | 8 | 2,722 | 18 | 499 s | VERIFIED |
| wall-band column `c_x ∈ [0.5,0.6]`, all `c_y ≤ 2.5` | 200 | 7,240 | 24 | 627 s | VERIFIED |

**Unscaled candidate, first full attempt (certified per root, stopped):** 2,946 of 5,000 root boxes done, 35,657
CPU-s, all clean except one root, `[1.2,1.3]×[1.5,1.6]`, `u ∈ [5/16, 3/8]` (`θ = 34.7–41.1°`): 241 boxes uncertified at
depth 24, where the cover's float minimum is `1.0033` (`zm_mixed_test.py unc`: a checker limit, not a hole).  The
same root, split into 8 sub-roots, **verifies** on the candidate `× 1.003` (16,074 boxes, 1,093 CPU-s) and `× 1.006`
(8,036 boxes, 629 CPU-s).  Heavy roots cluster around `c ≈ (1.2–1.4, 1.3–1.6)` at `θ ≥ 28°` (up to 13k boxes,
2,000 CPU-s each), the same anti-correlated geometry as r7's.

**Scaled file for the full sweep:** `runs/line-cover_m5_candidate_x1003.txt` = the candidate with every mass `× 1003/1000`
(exact), total **`522368729933/25000000000 = 20.89474919732 < 21`** (the brief allows up to `× 1.006`, total `20.957`).
Full D4 sweep with root pitch `1/20`, 16 `u`-bins (40,000 roots, for load balance), `--depth 24`, resumable
(`runs/zmm/m5/x1003_full.jsonl`), code `d4b97f9`, shas in `runs/zmm/m5/x1003_sha.txt`:

| run | roots | boxes | max depth | CPU | wall (10 procs) | leaves ADM / CHAIN / SPLIT / PIECE / EMPTY | result |
|---|---|---|---|---|---|---|---|
| `line-cover_m5_candidate_x1003.txt`, full D4 region | 40,000 (all, checked against `d4_roots`) | **461,204** | 21 | **70,917 s = 19.7 h** | 2.0 h | 146,017 / 23,518 / 53,951 / 9,058 / 18,058 | **VERIFIED-D4, 0 uncertified** |

**Certified (by `zm_mixed.py` alone):** every admissible closed unit square of `[0,5]²` with centre in `[0, 2.5]²` and
`θ = 2 arctan u`, `u ∈ [0, ½]` has `μ(Q) ≥ 1` for the mixed cover `runs/line-cover_m5_candidate_x1003.txt` of total
`20.89474919732 < 21`; with its (checked) D4 invariance this covers all poses (§2, Symmetry).  With the FORMAT.md
reduction that is `s(21) ≥ 5`, hence **`s(21) = 5`**, *provided* (i) the checker is correct — one implementation,
Python + `Fraction`, lemmas proved here and tested by randomised exact checks, no independent re-check yet — and
(ii) the reduction from point sets to these measures (Lean `packing_le_weight` generalisation, FORMAT.md) is done.
The unscaled candidate (20.832) is not certified: one root needs more margin than 0.3 % or a stronger checker (above).

Leaf dumps of three cells (a re-run with `--dump`, same settings): the hard cell `[1.2,1.3]×[1.5,1.6]` (42,338
boxes; 21,201 certified leaves), the germ `[1.45,1.55]²` (744 leaves) and the wall/dip area `[0.45,0.6]×[1.4,1.55]`
(1,018 leaves).  Their *stress* test (total float mass, min `1.0064–1.0078`) is uninformative at this margin (§8);
the audit's **component tests** of the same dumps (43.7 M + 0.22 M + 0.78 M exact checks, each certified component
bound against the exact mass of that component) found 0 violations, and its **rejection runs** on three holed
versions of this cover refused each exact hole (`ZM_MIXED_AUDIT.md` §3, rerun in certificate mode in §8).

Reproduce (resumable): `python3 search/zm_mixed.py cert runs/line-cover_m5_candidate_x1003.txt --d4 --disj
--chain-from 0 --depth 24 --pitch 1/20 --ubins 16 --nproc 10 --resume runs/zmm/m5/x1003_full.jsonl`.

**`m = 4` candidate** (`runs/line-cover_m4_candidate.txt`, 12.3654): a D4 run with the code before Lemma V / L′-hull
had 15,914 uncertified boxes at 3,000 / 3,200 roots when stopped for the `m = 5` work; not rerun.




## 8. Certificate mode and provenance (2026-09-27, after the audit)

`ZM_MIXED_AUDIT.md` found no soundness defect and five should-fixes; all are addressed without touching any lemma's
logic.

* **S2/S3 (provenance, resume).**  `cert` prints the sha256 of `zm_mixed.py`, `mixed_cover.py`, `zeromargin.py` and
  the input, the argv and every setting that changes what is checked (mode, depth, pitch, `u`-bins, CHAIN and
  `--chain-from`, `--theta-bias`, clipping, Lemmas T/L, SPLIT and its chain cap 400, the reach constant, cert mode,
  any partial region).  With `--resume` these go into a header, the first line of the jsonl; a restart refuses a
  file whose header (shas, settings, total) differs, or that has no header.  `--manifest FILE` writes header,
  argv, census, verdict and the records' sha256.
* **S5 (trusted surface).**  `--cert-mode`: Corollary T′ (points on germ lines) and the polygon code are
  unreachable (the polygon branch asserts), and a cover containing polygons is refused (exit 2).  Points lying on
  a segment line are then ordinary points (sound, possibly weaker; the `s(21)` cover has none).
* **S1.**  The cover, the run records and the exact checker files that ran are in `certificates/s21/`.
* **S4 (evidence).**  §0, §7.5 now cite the component and rejection tests; the leaf stress test stays as a
  heuristic tool (`zm_mixed.py stress`) but is not evidence about the checker on a cover with margin.
* **Nits.**  N1: dead `min_density`, `cone_slope` removed.  N4: comment at the `_gle0` call (weights `±1` are valid;
  zeromargin's docstring is not edited).  N2/N3 noted here: the code uses Lemma L′'s hull form (§2, proved), and
  Lemma V's split is valid for any `D` (the L′ term is valid at every pose, so soundness reduces to Lemma R).  N5
  (moved points withheld though unused: sound, wasteful) and N6 (the D4 check is sufficient, not necessary) are
  unchanged.  The section numbering (§7 before §6) is kept because §7.5 is cited elsewhere; Reproduce is §9.

**Tests rerun at the final code** (sha256 `ee3e2915…d760ac`, cores 0–9): `zm_mixed_test.py selftest --n 100 --seed 41`
PASS (all sections, 0 violations); audit `file` (independent parser: total, D4 and all four symmetries), `oracle`
(300 poses, 0 mismatches), `tprime` (7,403 checks, 0 violations), `boxes` (100 adversarial boxes, 4.59 M exact
checks, 0 violations, min slack 0); rejection runs **in cert mode**: R1 (`× 0.994`, tilted dip) 2,059 boxes, 36
uncertified, the exact hole `0.99972` in an UNCERT box; R3 (germ-pivot hole) 5,011 boxes, 1,703 uncertified, the
hole `0.98718` in an UNCERT box — both censuses identical to the audit's; component test of R1's 826 certified
leaves: 328,959 exact checks, 0 violations, min total in a certified leaf `1.00023`.

**The certified run** (`search/s21_cert_runs.sh zm_mixed`, `certificates/s21/zm_mixed_d4/`):
`--d4 --cert-mode --disj --chain-from 0 --depth 24 --pitch 1/20
--ubins 16`, code sha256 `ee3e2915…d760ac` (commit `fe2ac0b`), 10 processes on cores 0–9: 40,000 roots, **461,204 boxes, max
depth 21, leaves ADM 146,017 / CHAIN 23,518 / SPLIT 53,951 / PIECE 9,058 / EMPTY 18,058 / UNCERTIFIED 0, VERIFIED-D4**,
12.9 CPU-h, 1.4 h wall.  The census of every one of the 40,000 roots is identical to §7.5's run (cert mode changes
nothing here: the cover has no points on lines and no polygons).  `python3 search/s21_records.py zm_mixed
certificates/s21/zm_mixed_d4 certificates/s21/s21_mixed_cover_5.txt` re-checks the records from scratch.

## 9. Reproduce

```
bash search/s21_cert_runs.sh zm_mixed       # the shipped s(21) run -> certificates/s21/zm_mixed_d4/ (resumable)
python3 search/s21_records.py zm_mixed certificates/s21/zm_mixed_d4 certificates/s21/s21_mixed_cover_5.txt
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

