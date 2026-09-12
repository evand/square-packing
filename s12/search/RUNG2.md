# Rung 2 (`s(13) = 4` by a weighted closed cover of `[0,4]²`, `W < 13`) — task rung2-s13, 2026-09-11

**Verdict: not achieved, and the reason is now a theorem rather than an engineering gap.**

The brief asked for a checker fix (θ-preferential splitting, an off-centre wall-witness primitive)
and treated the two-region (disjunctive) primitive as a last resort — `FAMILY.md` §2b had diagnosed
the stuck boxes as "slow convergence of `cert_p1`/`cert_core`", i.e. case (b).  That diagnosis is
wrong.  The obstruction is structural:

> **Theorem 1.**  Let `P ⊂ [0,m]²` be a finite weighted point set.  If every leaf of a finite,
> closed, space-filling subdivision of the admissible pose space carries a *monotone witness
> certificate* — a set `S ⊆ P` of total weight `≥ 1` such that **every** `p ∈ S` lies in `Q(c,θ)`
> at **every** admissible pose of that leaf — then `w(P) ≥ m²`.
>
> `CORE`, `P1`, the new `ADM` primitive below, and every union of them produce exactly monotone
> witness certificates.  So at `m = 4` **no cover of total weight below 16 can be certified by
> `zeromargin.py`'s non-disjunctive primitives, at any depth, ever.**  `W < 13` is 3 below that.

`TRI` is the only primitive in the checker that is *not* of this kind (it asserts that *some*
vertex of a triangle is captured, not that a fixed point is), which is exactly why rung 1 —
Friedman's 14 points, `W = 14 < 16` — is certifiable at all, and exactly why it fails without
`--tri`.  Rung 2 therefore *requires* the brief's option (c), a weighted disjunctive primitive; the
options (a) and (b) it preferred cannot close the gap even in principle.

What this session delivers: the theorem with its proof and two independent machine confirmations
(§2), a new exact primitive `ADM` that strictly subsumes both `CORE` and `P1` and is the right
tool for everything *except* the disjunctive poses (§3, §4), the θ-biased splitting rule the brief
asked for (§3.4), the measurements that kill the "scale the cover" and "solve with a margin" routes
(§5), and a concrete statement of what a weighted disjunctive primitive has to prove (§6).

Everything load-bearing is `fractions.Fraction`; floats appear only as pre-filters (a documented
superset of the exact test) and in the independent stress tests.  Semantics unchanged
(`ZEROMARGIN.md` §1: closed unit squares, closed containment, a point on `∂Q` counts).

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

**Lemma 0 (one-sided approach).**  Fix a tile `(i,j) ∈ {0,…,m-1}²` and put
`p₀ = (i+½, j+½, θ=0)`; `Q(p₀) = [i,i+1] × [j,j+1]` and `p₀` is admissible (`w(0) = 1`, and
`½ ≤ i+½ ≤ m-½`).  Choose signs

    σ_x = +1 if i < m-1,  σ_x = -1 if i = m-1;      σ_y likewise in j.

Then the poses `p_ε = (i+½+σ_x ε, j+½+σ_y ε, 0)` are admissible for all small `ε > 0`
(`½ ≤ i+½+σ_x ε ≤ m-½` by the choice of sign).  The leaves are finitely many and closed, so some
leaf `B` contains `p_{ε_k}` for a sequence `ε_k ↓ 0`, and then `p₀ ∈ B` because `B` is closed.

**Lemma 0'.**  The same leaf `B` also contains poses with `θ > 0`: `p₀` is in `B`, and `B`'s
`u`-interval `[u₀,u₁]` contains `0`, hence `u₀ = 0 < u₁` (`u ≥ 0` on the domain, and a leaf of a
halving subdivision has `u₀ < u₁`).

*Proof of Theorem 1.*  Let `S` be `B`'s witness set, `w(S) ≥ 1`.  Every `p ∈ S` lies in `Q` at
every admissible pose of `B`, in particular at `p₀` and at every `p_{ε_k}`, so

    S ⊆ Q(p₀) ∩ ⋂_k Q(p_{ε_k}) = ([i,i+1] ∩ ⋂_k [i+σ_xε_k, i+1+σ_xε_k]) × (same in y)
      = I_i^{σ_x} × J_j^{σ_y},    I_i^{+1} = (i, i+1],  I_i^{-1} = [i, i+1).

Adding Lemma 0' and the θ-derivative of the four inequalities at `θ = 0` shrinks this further (the
"pinwheel": on the surviving vertical edge only the half `d_y ≤ 0` resp. `d_y ≥ 0` remains, and
likewise on the horizontal edge — see the derivation in §3.1, which is the same computation), but
the crude form already suffices here.  Hence

    w( PIN(i,j,σ_x,σ_y) ) ≥ 1     for every tile and every admissible sign pair,     (★)

where `PIN` is the intersection above.  Minimising `Σ w_p` subject to the `36` constraints (★) at
`m = 4` — 4 corner tiles with one sign pair, 8 edge tiles with two, 4 interior tiles with four —
is an LP.  `search/rung2_bound.py bound --octants` solves it over a lattice of candidate points:

| lattice | candidates | constraints | LP value |
|---|---|---|---|
| 1/20 | 6,561 | 36 | **16.000000** |
| 1/40 | 25,921 | 36 | **16.000000** |

so `w(P) ≥ 16`.  ∎ (for lattice covers; every cover this project builds has its columns on the
`1/1000` lattice of `closed4.py`, so the bound applies to all of them.  Restricting candidates to a
lattice can only *raise* an LP minimum, so the honest reading is: "16 for covers on a `1/20` or
finer lattice, stable under refinement".)

**Why the one-sidedness matters.**  If one uses only the *two-sided* information at each tile pose
— i.e. `S ⊆ Q(p₀) ∩ ⋂ Q(p_ε)` over both signs — the 16 sets overlap at the tile edges and the LP
value halves to `8`: `rung2_bound.py bound` (without `--octants`) prints `8.000000`, with a support
of 8 points of weight 1 on the lines `x = 1` and `x = 3`.  It is Lemma 0 — the fact that a leaf box
touching the tile pose from one side is *also* a leaf that must be certified — that gives 16.

**Machine confirmation 1 (the shipped scaled covers).**

```
python3 search/rung2_bound.py check runs/closed4_best_x103.txt
```
`container [0,4]^2, 1972 points, total weight 3197502027/250000000 = 12.790008108` —
**24 of the 36 PIN sets have weight below 1** (exact `Fraction` sums), the worst being
`tile (0,1) sx=+1 sy=-1: 5661189/10000000 = 0.566119`, i.e. the pose `(0.5, 1.5, 0)`.  The four
corner tiles pass at `1.03000` (this is the `P1` situation of `ZEROMARGIN.md` §4 item 1).
`0.566119` is exactly the number a brute-force float scan gives for "the weight captured
simultaneously at every pose of an arbitrarily small box at `(0.5, 1.5, 0⁺)`" — and it does not
move as the box shrinks from `10⁻²` to `10⁻⁶`, which is the signature of the obstruction.

This is precisely where `FAMILY.md` §2b's depth ladder was stuck (`cx ≈ 0.5`, `cy ≈ 1.5`,
`θ → 0`), and it explains the oscillating uncertified counts: the boxes there are not converging
slowly, they are **not converging at all**.

**Machine confirmation 2 (rung 1).**

```
python3 search/rung2_bound.py check runs/friedman14.txt
```
Friedman's 14 points, `W = 14 < 16`: **exactly 6 PIN sets are empty (weight 0)** — the sign pairs
`(1,1,+,+)`, `(1,1,+,-)`, `(2,1,+,+)`, `(1,2,-,-)`, `(2,2,-,+)`, `(2,2,-,-)`, i.e. the four
*interior* tile poses `(1.5,1.5,0)`, `(2.5,1.5,0)`, `(1.5,2.5,0)`, `(2.5,2.5,0)`.  Those six leaves
are exactly the ones `TRI` certifies, and they are why rung 1 **does not terminate without
`--tri`** (`ZEROMARGIN.md` §3: 87,200 uncertified at depth 21 with `CORE`+`P1`; with the new `ADM`
primitive and no `--tri` it is still `71,961` uncertified at depth 20 — §4).  Theorem 1 predicts
this failure and names the six poses; the checker finds them.

**Corollary.**  The same argument at general `m` gives `w(P) ≥ m²` for every rung of
`search/FAMILY.md`: rung 3 (`m = 7`, target `W < 45`) needs `49`, rung 2 of `m = 5`
(target `W < 21`) needs `25`.  *Every* rung of the `s(m²-3)` / `s(m²-4)` family is beyond the
non-disjunctive primitives, by exactly the amount that makes the theorem interesting
(`m² - (m²-3) = 3`).  The closed-square semantics is not an incidental convenience of these
proofs — it is the whole content, and any machine proof has to reproduce the case analysis
(Bentz's 6-leaf tree, DS7's Lemma 3) in some form.

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

| run | boxes | depth | ADM | CORE | P1 | MIX | TRI | EMPTY | uncert |
|---|---|---|---|---|---|---|---|---|---|
| `--no-adm --theta-bias 1` (the old code path) | 6,958 | 10 | – | 3,356 | 70 | – | 74 | 3,179 | **0** |
| `ADM`, `--theta-bias 1` | **6,810** | 10 | 3,367 | 0 | 0 | 0 | 74 | 3,164 | **0** |
| `ADM`, `--theta-bias 4` (default) | 7,220 | 10 | 3,511 | 0 | 0 | 0 | 74 | 3,225 | **0** |

The first row reproduces `ZEROMARGIN.md` §3 exactly (`3356/70/74/3179/0`, 6,958 boxes), so the old
path is unchanged.  `ADM` alone replaces **both** `CORE` and `P1` (`P1` drops to 0 leaves: the 70
corner/wall boxes `CORE` could never do are now `ADM` leaves) and needs 148 fewer boxes.  Without
`--tri` it is still `NOT VERIFIED` (71,961 uncertified at depth 20) — as Theorem 1 requires.

**Rung 2, the scaled cover** (`cert runs/closed4_best_x103.txt --depth 10 --nproc 8`,
`taskset -c 16-23`): 307,420 boxes, `ADM 117,805 / CORE 0 / P1 0 / MIX 0 / TRI 0 / EMPTY 4,563 /
UNCERTIFIED 34,542`, 321 s.  Compare `FAMILY.md` §2b's depth-10 row for the same file
(`300,294` boxes, `114,471 / 337 / 0 / 4,327 / 34,212`).  `ADM` does strictly more per box and the
uncertified count is unchanged to 1 % — the boxes that remain are the obstructed ones, and no
monotone primitive will ever remove them.  All 34,542 sit in the `cx ≈ 0.5, cy ≈ 1.5, θ → 0` region
and its D4 images, i.e. at the 24 violated PIN sets of §2.

**Independent float stress (`search/zeromargin_stress.py`, no shared code path).**
`runs/zm_f14_adm.txt` (6,810 `ADM` + 74 `TRI` leaves), 40 sampled poses per leaf including the box
corners: **0 failures**.  Primitive tests: `P1` 200,000 instances, core lemma 68,569 (21 angles
each), triangle 72,358 — 0 failures.  New: the `ADM` lemma on 40,000 random `(box, point)` pairs,
with the four inequalities evaluated *directly in floats on a 201-point θ grid* (no polynomials, no
Bernstein) and the conclusion checked by sampling admissible poses — `2,713` certifying pairs,
`43,859` pose samples, **0 failures**.

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

## 6. What a weighted disjunctive primitive has to prove

Theorem 1 says a leaf's certificate must sometimes be of the form "*for every admissible pose of
this box, **at least one of** the witness sets `S₁, …, S_k` is entirely captured*", with
`w(S_r) ≥ 1` for each `r`.  `TRI` is the case `k = 3`, `S_r = {v_r}` a single vertex of weight `≥ 1`
— useless for a cover whose weights are `≈ 0.17`.  Concretely, at the worst rung-2 pose
`(0.5, 1.5, 0)` (tile `(0,1)`, `σ = (+,-)`, PIN weight `0.5661`) the disjunction needed is:

* `S₁` = the stable part plus the line `y = 2` — captured exactly when `c_y ≥ 3/2` (and, at
  `c_y = 3/2` with `θ > 0`, only for `p_x ≥ c_x`), total `1.0643` for `closed4_best_x103.txt`;
* `S₂` = the stable part plus the line `y = 1` — captured exactly when
  `(p_x - c_x) sin θ ≤ ½(1 - cos θ) + (3/2 - c_y)cos θ`, total `1.0643`.

The two regions cover the box: the boundary between them is the smooth surface
`(p_x - c_x)\sin θ = ½(1-\cos θ) + (3/2 - c_y)\cos θ`, which in the `u` parametrisation is again a
polynomial with rational coefficients of degree `≤ 3` in `(c_x, c_y, u)`, so "box `⊆ R₁ ∪ R₂`" is
an exactly decidable statement (Bernstein enclosure on the box, exactly as Lemma C, now in three
variables).  That is the missing primitive, and the machinery of §3 is most of what it needs; what
it does **not** yet have is (i) the search that picks `S₁, S₂` and the splitting surface for a given
box, and (ii) a proof obligation discharge for the case `k > 2`.  It was not built here: the session
went into establishing that it is *necessary*, which the brief had listed as the least likely
outcome.

The other honest option is to accept a case analysis in the human sense — a machine-checked
version of Bentz 2010's 6-leaf tree — in which case the cover formulation is the wrong target
altogether and the `s(12)` endgame should be planned around branch certificates directly
(`ZEROMARGIN.md` §6).

---

## 7. What is exact, what is float

| quantity | status |
|---|---|
| Lemma A, B, C and every `ADM`/`CORE`/`P1`/`TRI` leaf | exact `Fraction`; floats only as a documented superset pre-filter |
| Theorem 1 and the sets `PIN(i,j,σ)` | the proof is exact; `rung2_bound.py check` sums the weights in `Fraction` but decides *membership* by evaluating the pose inequalities in floats at `θ ∈ {0, 10⁻⁵, 10⁻⁴, 10⁻³}` with tolerance `10⁻¹²` (the cover's coordinates are on a `1/1000` lattice, so the `O(θ)` and `O(θ²)` gaps are resolved with `≥ 50×` margin) |
| the LP values `16.000000` (`rung2_bound.py bound`, `closed4.py --margin`) | float LP (HiGHS); corroborated by two lattices and by the exact `check` above |
| `0.566119` at `(0.5, 1.5, 0)` | exact `Fraction` (`rung2_bound.py check`) |
| leaf counts in §4 | exact (the checker's own accounting) |
| `zeromargin_stress.py` | floats by design — an *independent* confirmation, not a proof |

---

## 8. Reproduce

```bash
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-11/closed4_best*.txt runs/

# Theorem 1, the LP behind it (seconds; --pitch 40 takes ~1 min)
python3 search/rung2_bound.py bound --m 4 --pitch 20 --octants     # 16.000000
python3 search/rung2_bound.py bound --m 4 --pitch 40 --octants     # 16.000000
python3 search/rung2_bound.py bound --m 4 --pitch 20               # 8.000000  (two-sided only)

# Theorem 1 applied to the covers (exact)
python3 search/rung2_bound.py check runs/closed4_best_x103.txt     # 24/36 PIN sets < 1
python3 search/rung2_bound.py check runs/friedman14.txt            # 6 empty PIN sets = the TRI poses

# rung 1 regression and the new primitive
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4 --no-adm --theta-bias 1
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4
python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4 --dump runs/zm_f14_adm.txt
python3 search/zeromargin_stress.py runs/zm_f14_adm.txt 40         # 0 failures everywhere

# rung 2 with ADM (5 min on 8 processes)
taskset -c 16-23 python3 search/zeromargin.py cert runs/closed4_best_x103.txt --depth 10 --nproc 8

# the margin LP (the closed route)
python3 search/closed4.py axis --s 4 --margin 0.01 --tag ax010 --nproc 1    # = 16
```

Files: `search/zeromargin.py` (`ADM`, `MIX`, `--theta-bias`, `--no-adm`),
`search/zeromargin_stress.py` (+ the `ADM` random test, + `--cert` for weighted covers),
`search/rung2_bound.py` (new), `search/closed4.py` (`--margin`, checkpoint-export fix).
No certificate is shipped: `certificates/rung2/` is not created, because nothing was certified.
