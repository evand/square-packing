# `zmx2`: an independent second exact checker for mixed covers (2026-09-27)

Brief: `tasks/s21-finish/xcheck.md`.  Re-verify `runs/line-cover_m5_candidate_x1003.txt` (the `s(21) = 5` mixed
cover, `tasks/line-cover/FORMAT.md` v1) with a checker that shares no code and no lemma write-up with
`search/zm_mixed.py`.  **Independence:** `zm_mixed.py`, `zm_mixed_test.py`, `ZM_MIXED.md` and every `ZM_MIXED_AUDIT*`
were never opened.  Read: `FORMAT.md` (both), `LINE_COVER.md`, the header and CLI of `verify2/src/main.rs`
(`zmcheck`), `S32_EXACT.md` §7 (the D4 argument).  Own parser, own evaluator; `search/mixed_cover.py` was used once,
only to eyeball the candidate's structure.

Code: `verify2/src/bin/zmx2.rs` (a separate binary; `zmcheck`'s `src/main.rs`, `Cargo.toml`, `Cargo.lock` untouched,
checked by test T0), `search/zmx2_tools.py` (exact rational evaluator, harness, toy and perturbed covers),
`search/zmx2_tests.sh` (tests, §8), `search/zmx2_run.sh` (the official runs, writes `search/zmx2_manifest.txt`).

## 0. Verdict

> **The candidate is certified: `VERIFIED-D4` (2,500 roots, 1,826,222 boxes, 0 uncertified, 11 CPU-s) and, without
> the D4 argument, `VERIFIED` (20,000 roots, 14,709,448 boxes, 0 uncertified, 87 CPU-s).**  Total
> `522368729933/25000000000 = 20.89474919732 < 21`, so with the (generalised) reduction of `FORMAT.md`, `s(21) ≥ 5`,
> hence `s(21) = 5`.  This agrees with `zm_mixed.py`'s `VERIFIED-D4`.
>
> * **Tight.**  The checker loses `< 0.1 %`: the candidate scaled by `0.9955` (true minimum `≈ 1.0009`) still
>   verifies; scaled by `0.994` it is refused, and the deepest refused boxes sit at the cover's known dip
>   `(0.5736, 1.4406, θ = 9.21°)` (the D4 image of `LINE_COVER.md`'s `(3.559, 0.573, 9.13°)`), where the exact
>   rational `μ = 0.99969 < 1`.
> * **Germs are no obstacle.**  A tile germ with a hole (line `x = 2` zeroed on `y ∈ [1.4, 1.6]`) is refused exactly
>   at `(1.5, 1.5, θ → 0⁺)`; exact `μ = 1.749` at `θ = 0` (edges count) and `0.9176` at `θ = 0.003°`.
> * **Cross-checks on point covers.**  The same binary re-verifies the shipped `s(13)` cover (`VERIFIED-D4`,
>   47,162 boxes, 3 CPU-s; agrees with `zmcheck` on a column) and the shipped `s(32)` cover
>   (`certificates/s32/s32_closed_cover_6.txt`, `--pair-points`: `VERIFIED-D4`, 1,405,342 boxes, 873 CPU-s;
>   `zmcheck`/`zeromargin.py` needed 150–250 CPU-h), and refuses both when scaled below their minima.
> * **Soundness harness:** exact rational `μ` at 27,301 random admissible rational poses inside 1,600 boxes (random
>   small boxes, and the tightest certified leaves) in the test suite, plus ≈ 80,000 more in exploratory runs, never
>   below the box's bound (§8).  Test suite: 42/42 pass.
>
> Arithmetic: points and all mass comparisons exact integer; chord-endpoint geometry by outward-rounded IEEE
> binary64 interval arithmetic, rounded outward onto an integer grid (§5, Lemma R).  So "exact" in the sense of a
> rigorous enclosure, not of rational arithmetic throughout; this is the one place it differs from `zmcheck`.

## 1. Statement and set-up

`μ = Σ w_p δ_p + Σ_segments (w/|seg|)·λ|seg` on `[0,s]²` (all masses `/W`, coordinates `/D`).  Claim: for every closed
unit square `Q ⊆ [0,s]²`, `μ(Q) ≥ 1`.  Scope of `zmx2`: points and **axis-parallel** segments (non-axis segments and
polygons are refused with `ERROR: unsupported`; the candidate has neither).

**Poses.**  `Q(c, θ) = c + R(θ)[−½, ½]²`; `Q(c, θ + 90°) = Q(c, θ)` so `θ ∈ [0°, 90°)` is all.  `u = tan(θ/2)`,
`C = cos θ = (1−u²)/(1+u²)`, `S = sin θ = 2u/(1+u²)`, `w(u) = C + S`.  Admissible (`Q ⊆ [0,s]²`) iff
`w/2 ≤ c_x, c_y ≤ s − w/2`.  A point `p` is in `Q` iff `|X|, |Y| ≤ ½`, `(X, Y) = R(−θ)(p − c)`.

**Boxes.**  `B = [x0,x1] × [y0,y1] × [u0,u1]`, rational, centre denominators `10·2^j`, `u` denominators `8·2^k`
(`j, k ≤ 28`).  Roots: centre cells of pitch `1/10`, `u` bins `[k/8, (k+1)/8]`, `k = 0..3` (`u ≤ ½`, `θ ≤ 53.13°`).
For each box a lower bound `L(B)` is computed with

> **(Main property)** `L(B) ≤ μ(Q(c, θ))` for every **admissible** pose `(c, u)` of `B`,

and `B` is a leaf when `L(B) ≥ 1` (certified) or when no pose of `B` is admissible (Lemma E, `EMPTY`); otherwise `B`
is halved (the dimension with the largest width, `u`-widths weighted by `κ = 1.5`), down to depth 40, below which the
box is reported `UNCERTIFIED` with its exact coordinates and a float estimate of its minimum.  Halves are closed and
cover `B`, so a root with no uncertified leaf is certified at every pose.

**Theorem.**  If every root of the region (§7) is certified, the claim holds.  (Main property + covering by closed
leaves + §7.)

`L(B) = P(B) + Σ_families Φ(B)`: `P` = the exact weight of the ordinary points certainly in `Q` (Lemma P); for each of
the two families of lines (vertical, horizontal) a lower bound on the mass of the lines (segment densities and
line atoms, §4.5) from chords (Lemmas C, M, S, Z, W, H), combined by Lemma DP; horizontal lines via Lemma T.  The
pieces of `μ` counted by different terms are disjoint (Lemma DP), so the sum is a lower bound.

## 2. Points (Lemma P)

With `a = p_x − c_x`, `b = p_y − c_y` and `N = 1+u²`, `2N·X − N` etc. give four *violation polynomials*
`G0 = (−2a−1)u² + 4bu + (2a−1)` (`X ≤ ½`), `G1 = (2a−1)u² − 4bu + (−2a−1)` (`X ≥ −½`), `G2 = (−2b−1)u² − 4au + (2b−1)`
(`Y ≤ ½`), `G3 = (2b−1)u² + 4au + (−2b−1)` (`Y ≥ −½`); `p ∈ Q ⟺ all Gk ≤ 0`.

**Lemma P.**  `max_B Gk ≤ 0` iff at each of the four corners `(x, y)` of the centre rectangle,
`max_{u ∈ [u0,u1]} Gk ≤ 0`, and the latter is decided exactly by integers.
*Proof.*  `Gk` is affine in `(c_x, c_y)` for fixed `u`, so its maximum over the rectangle is at a corner; at a corner it is
a quadratic `q(u) = p2u² + p1u + p0` whose maximum on `[u0, u1]` is `max(q(u0), q(u1))`, or, if `p2 < 0` and the vertex
`−p1/(2p2)` is interior, `q(vertex) ≤ 0 ⟺ 4p2p0 − p1² ≥ 0`.  With `a = A/E`, `b = B/E` (`E = D·10·2^j`) and
`u = U/(8·2^k)` all tests are integer comparisons of size `< 2^100` (`|A|, E < 2^45`, `U < 2^31`), in `i128`. ∎

`P(B)` sums the weights of points with all four `max_B Gk ≤ 0`.  A float pre-filter only decides which points are
*tested* or *dropped* (dropping only lowers `L`); a point certain for a box stays certain for its sub-boxes and is
inherited.  Exactly the points are counted that are in `Q` at every pose of `B` (admissible or not).

## 3. Chords of an axis-parallel line

**Lemma C.**  For `0 < θ < 90°` and the vertical line `x = ℓ`, put `d = c_x − ℓ`.  Then
`{y : (ℓ, y) ∈ Q} = [c_y + max(f1, g1), c_y + min(f2, g2)]` (empty if the ends cross), with
`f1 = (dC − ½)/S = (d−½)/(2u) − (d+½)u/2`, `f2 = (dC + ½)/S = (d+½)/(2u) + (½−d)u/2` (from `X`),
`g1 = (−dS − ½)/C = −dT − R/2`, `g2 = −dT + R/2` (from `Y`), `T = tan θ = 2u/(1−u²)`, `R = sec θ = (1+u²)/(1−u²)`.
*Proof.*  For `p = (ℓ, y)`: `X = −dC + (y − c_y)S`, `Y = dS + (y − c_y)C`; `|X| ≤ ½` and `|Y| ≤ ½` are intervals in
`y − c_y` since `S, C > 0`; substitute `C, S` in `u`. ∎

Write `L1 = c_y + f1`, `H1 = c_y + f2`, `L2 = c_y + g1`, `H2 = c_y + g2` (the crossings of the line with the four edge
lines of `Q`).  All four depend on `c_y` additively and on `(d, u)` only through `f`, `g`, and `c_y`, `d`, `u` are
independent box coordinates.

**Lemma M (enclosures).**  For fixed `u`, `f1, f2, g1, g2` are affine in `d`, so their range over `B` is the hull of
their ranges over `u ∈ [u0, u1]` at `d = d0` and `d = d1`.  For fixed `d`:
* `f = α/u + βu` (`α, β` constants): if `αβ ≤ 0` it is monotone, range = hull of the endpoint values; if `αβ > 0` it
  has one critical point `u* = √(α/β)` and the range is the hull of the endpoint values and `f(u*)`.  At `u0 = 0` the
  endpoint value is the limit `sign(α)·∞` (`0·` if `α = 0`: `f = βu`).
* `g = −dT ∓ R/2`: derivative `∝ −2d(1+u²) ± 2u`; for `d ≤ 0` (`+R/2`) resp. `d ≥ 0` (`−R/2`) monotone, otherwise one
  critical point in `(0, 1)`, `u* = 2|d| / (1 + √(1 − 4d²))` (none if `|d| > ½`); range = hull of endpoints and `g(u*)`.
  `T`, `R` are increasing on `[0, 1)`, so `g` over a small interval is enclosed by `−d·[T(lo), T(hi)] ∓ [R(lo), R(hi)]/2`.
*Proof.*  Elementary calculus; `u*` is enclosed by interval arithmetic (Lemma R), and the value over the enclosing
interval is enclosed naively (a superset).  When a sign is not certain (a coefficient interval straddles 0) the
code falls back to the naive enclosure over the whole `[u0, u1]` (a superset). ∎

From these: `ζ-range = [y0 + inf f1, y1 + sup f1]` (the range of `L1`), `L1hi = y1 + sup f1`, `H1lo = y0 + inf f2`,
`L2hi = y1 + sup g1`, `H2lo = y0 + inf g2`.

## 4. The line bounds

Each line `ℓ` carries the measure `ν_ℓ` = its segment density (piecewise constant, cumulative `F_ℓ`, piecewise linear
and continuous) plus its atoms (§4.5).

**4.1 Lemma S (single line).**  For every pose of `B` with `u > 0`, the chord contains `[lo*, hi*]`,
`lo* = max(L1hi, L2hi)`, `hi* = min(H1lo, H2lo)`, so the density mass is `≥ F(hi*) − F(lo*)` (or 0); atoms count when
all four of their conditions hold on `B` (Lemma P per condition).  *Proof.*  Monotonicity of `F`; each end is
pessimised separately. ∎

**4.2 Lemma Z (pair lemma: lines at distance exactly 1).**  For lines `a: x = ℓ` and `b: x = ℓ+1` and `u > 0`, let
`ζ = L1_a`.  Then `H1_b = ζ + u` **exactly**.  Hence, for every pose of `B`,

`ν_a(Q) + ν_b(Q) ≥ Φ(ζ)`, `Φ(z) = [F_a(h_a) − F_a(max(z, l_a))]⁺ + [F_b(min(max(z + u0, H1lo_b), H2lo_b)) − F_b(l_b)]⁺ + atoms(z)`,

with `l_a = L2hi_a`, `h_a = min(H1lo_a, H2lo_a)`, `l_b = max(L1hi_b, L2hi_b)`, and therefore
`ν_a(Q) + ν_b(Q) ≥ inf_{z ∈ ζ-range} Φ(z)`.
*Proof.*  `d_b = d_a − 1`, and `f2(d − 1, u) = (d − ½)/(2u) + (3/2 − d)u/2 = f1(d, u) + u`.  The chord of `a` is
`[max(ζ, L2_a), min(H1_a, H2_a)] ⊇ [max(ζ, l_a), h_a]`; the chord of `b` is `[max(L1_b, L2_b), min(H1_b, H2_b)]` with
`H1_b = ζ + u ≥ ζ + u0` and `H1_b ≥ H1lo_b`; `F` monotone; the actual `ζ` lies in the `ζ`-range. ∎

This is what makes tile germs finite.  Near `θ = 0` the crossing `ζ` of the left edge with `x = ℓ` runs over all of
`ℝ` (it is `≈ c_y − (½ − d)/θ`), which kills any bound that pessimises the two lines separately; but whatever `ζ` is,
`b` gets the complementary part of the edge.  `Φ` is exactly the `θ → 0⁺` germ-limit function of `LINE_COVER.md` §2
(up to `O(box)`), and `inf Φ` is computed exactly (below).

**4.3 Computing `inf Φ` exactly.**  After Lemma R all parameters are integers on the grid `1/(D·2^K)`, `K = 30`, and
the breakpoints of the lines are grid points.  `Φ` is continuous and piecewise linear between the *candidates*:
`zlo`, `zhi`, `l_a`, `h_a`, `H2lo_b − u0`, `H1lo_b − u0`, `l_b − u0`, the breakpoints of `a`, the breakpoints of `b`
minus `u0`, and the atom positions (below).  On each open interval between consecutive candidates `Φ` is affine, so
`inf Φ` is the minimum of the one-sided limits at the candidates; each is an exact integer (`F` at a grid point is
`F(b_k) + dens_k·(y − b_k)`, integer in units `1/(W·Lc·2^K)`, `Lc` = lcm of segment lengths, `20` here).

**4.4 Lemma DP.**  Lines pair only with a partner at distance exactly 1, so the lines within reach (`|c_x − ℓ| ≤ 0.75`;
a square meets a line only if `|c_x − ℓ| ≤ w/2 ≤ 0.7072`) split into chains `ℓ, ℓ+1, ℓ+2, …` (positions congruent
mod 1).  Along each chain a DP takes the best partition into singles and adjacent pairs.  Valid: distinct lines carry
disjoint parts of `μ`; a vertical and a horizontal line meet in one point, where densities put no mass, and an atom
belongs to exactly one line; so any partition gives `Σ` group bounds `≤ μ(Q)`.  (A first version ran the DP over
the sorted list of all lines, so with extra `--pair-points` lines the integer pairs were no longer adjacent; that
lost completeness, never soundness, and was found by the `s(13)` run.)

**4.5 Line atoms.**  A point on a line (`x ∈` {interior unit grid lines, segment lines}, then the same for `y`; with
`--pair-points` also every abscissa/ordinate with a partner point at distance 1) is an atom of that line instead of
an ordinary point.  Per box each atom gets its four conditions decided exactly (Lemma P per `Gk`; `G0 = H1`, `G1 = L1`,
`G2 = H2`, `G3 = L2` in the line's frame).  In a pair: atoms with all four certain count always; `a`-atoms with
`H1, H2, L2` certain count iff `z ≤ y` (the missing condition is `y ≥ L1_a = ζ`); `b`-atoms with `L1, L2, H2` certain
count iff `y ≤ z + u0` (the missing condition is `y ≤ H1_b = ζ + u`).  The `a`-atom sum is left-continuous in `z`, the
`b`-atom sum right-continuous, so the one-sided limits at a candidate `c` are `A(c) + B(c⁻)` and `A(c⁺) + B(c)`; the
code takes the minimum of both at every candidate (also at `zlo`, `zhi`, where one of them may be a limit from
outside the range; it is still `≤ Φ` there, so the result stays a lower bound).  The candidate has no points on lines;
atoms are what let `zmx2` check point covers (`s(13)`, `s(32)`) at all (`--no-atoms` fails at every one-cut germ).

**4.6 Lemma W (walls).**  If `ℓ = X_L + 1` (distance 1 from the left wall) then for every admissible pose with
`u > 0`: `H1_ℓ ≥ c_y + q(u)`, `q = (u(1−C) + C)/2 = (1 − u² + 2u³)/(2(1+u²))`.  If `ℓ = X_R − 1`:
`L1_ℓ ≤ c_y − q(u)`.
*Proof.*  `f2` and `f1` are increasing in `d` (`∂/∂d = C/S > 0`).  Admissibility `c_x − X_L ≥ w/2` gives
`d ≥ (C+S)/2 − 1`, and `f2((C+S)/2 − 1) = ((1−C)² + SC)/(2S) = u(1−C)/2 + C/2` (using `(1−C)/S = u`).  The right wall
is the mirror computation. ∎  On `[u0, u1]`, `q ≥ (1 − u1² + 2u0³)/(2(1 + u1²))`.  Used as `H1lo := max(H1lo, y0 + q_lo)`,
`L1hi := min(L1hi, y1 − q_lo)`, the `ζ`-range upper end likewise, and for atoms (`G0` resp. `G1` treated as certain
when `y ≤ y0 + q_lo` resp. `y ≥ y1 − q_lo`).  This is what closes the wall germs: at `θ → 0⁺` a square against the
wall has its far edge beyond `x = 1`, and without admissibility the box always contains the inadmissible poses where
it is not.

**4.7 Lemma T (horizontal lines).**  The rotation `ρ(x, y) = (−y, x)` satisfies `ρQ(c, θ) = Q(ρc, θ)` (rotations commute;
the unit square is `ρ`-invariant) and maps `[0,s]²` onto `[−s, 0] × [0, s]` (admissibility preserved) and the
horizontal line `y = a` onto the vertical line `x' = −a`, the point `(x, a)` to `(−a, x)`.  So the horizontal family is
the vertical machinery applied to the box `x' ∈ [−y1, −y0]`, `y' ∈ [x0, x1]`, walls `x' = −s` (left) and `0` (right),
line position `−a`, coordinate along the line = the original `x`.  One set of formulas serves both families.

**4.8 Lemma E (empty boxes).**  `w(u)` is increasing on `[0, √2−1]` and decreasing on `[√2−1, 1]` (`w' ∝ 2 − 4u − 2u²`),
so `w ≥ w_lo = min(w(u0), w(u1))` on `B`; if `x1 < w_lo/2` (or `y1 < w_lo/2`, or `x0 > s − w_lo/2`, or
`y0 > s − w_lo/2`) no pose of `B` is admissible.

## 5. Arithmetic (Lemma R)

(i) *Intervals.*  Every `f64` operation is IEEE-754 round-to-nearest (Rust on x86-64: SSE2, no FMA contraction, no
fast-math); the exact result lies in `[next_down(r), next_up(r)]` of the rounded `r` (also on overflow to `±∞`), so
widening every result by one ulp in the safe direction gives an enclosure; `sqrt` is correctly rounded too.  Rational
inputs `n/d` have `|n|, d < 2^53` (asserted), so they convert exactly and the quotient is one rounding (widened by 4 ulps).
An exact zero factor gives an exact zero, so an exactly-zero `α` (a box edge through a tile centre, `d = ½`) stays zero
(without this the `ζ`-range of such boxes was `ℝ`; found on the first pilot, a completeness loss only).
(ii) *Onto the grid.*  A lower (upper) bound `v` becomes `⌊next_down(fl(v·G))⌋` (`⌈next_up(·)⌉`), `G = D·2^30` exact;
since `v·G ∈ [next_down, next_up]` of its rounding, this is `≤ v·G` (`≥`).  `u0` enters as `⌊u0·G⌋` (exact integers).
(iii) *Clamp.*  Values beyond `±(4s + 10)` are clamped to it.  All `F` are supported in `[0, s]`, all shifts are
`≤ 1`, and `Φ` is constant for `z ≤ −s − 1` and for `z ≥ s + 1`, so no `F` value and no `inf Φ` changes.
(iv) Everything after (ii) is exact `i128`: `F` values `< 2^76`, sums `< 2^80`.
A float **never** certifies anything except through (i)–(ii).

## 6. `θ = 0` exactly (Lemma H)

Boxes with `u0 = 0` contain axis-parallel poses, where Lemma C does not apply (`S = 0`) and segments on edges count
in full.  At `θ = 0` the chord of `x = ℓ` is `[c_y − ½, c_y + ½]` if `|d| ≤ ½`, else empty; `g1(0) = −½`, `g2(0) = ½`,
and Lemma M evaluates `g` at `u = 0` directly.  **Lemma H.**  At every admissible pose `(c, 0)` of `B` the true line
mass is `≥` the single bound, and `≥ Φ(z)` for some `z` in the computed `ζ`-range.
*Proof.*  The box contains `(c, u)` for all small `u > 0`, where `f1(d, u) → sign(d − ½)·∞`, `f2(d, u) → sign(d + ½)·∞`
(or to `0`, `c_y`-relative, at `d = ½` resp. `−½`); so the sup/inf of Lemma M include these limits.  Single: if
`d > ½`, `L1hi = +∞` (empty); if `d < −½`, `H1lo = −∞` (empty); if `|d| ≤ ½`, `[lo*, hi*] ⊆ [L2hi, H2lo] ⊆ [c_y − ½, c_y + ½]`.
Atoms are decided by exact `Gk` tests at all poses of `B`, `u = 0` included.  Pair (`d = d_a`): for `d < ½` take
`z = zlo = −∞` (then `b` counts nothing: `H1lo_b = −∞` because `d_b + ½ < 0`, and `a` counts a sub-chord of its full chord
if `|d| ≤ ½`, nothing if `d < −½` since then `H1lo_a = −∞`); for `d > ½` take `z = zhi = +∞` (`a` counts nothing; `b`
counts a sub-chord of its full chord if `d_b ≤ ½`, nothing if `d_b > ½` since then `L1hi_b = +∞`); for `d = ½` both
lines are full and any `z` works.  Walls: at an admissible `(c, 0)` the lemma-W line has `|d| ≤ ½` or lies on the side
where the unmodified bound already gives 0, and the `W`-modified quantities only enter chords contained in the full
chord `[c_y − ½, c_y + ½]`. ∎

So the checker's value at a germ box is the `θ → 0⁺` value (the smaller one), and `θ = 0` itself is covered.

## 7. Regions: D4 and unreduced

**`--d4`** (as `S32_EXACT.md` §7, for measures).  (i) `check_d4` verifies exactly that `μ` is invariant under
`x ↦ s − x` and `x ↔ y` (generators of D4): the aggregated point weights (atoms included), and each line's density in
canonical form (maximal intervals of constant non-zero density), vertical `x = a` ↦ `x = s − a` with the same
`y`-density, horizontal `y = b` ↦ itself reflected, and vertical `x = a` ↔ horizontal `y = a`.  Otherwise `ERROR`,
exit 2.  (ii) A D4 element `g` is an isometry of `[0,s]²`, `g Q(c, θ) = Q(gc, θ')` with `θ' = θ` (rotations) or
`90° − θ` (reflections), and `μ(gQ) = μ(Q)`.  (iii) Every pose maps to `θ' ≤ 45°` (`u ≤ √2 − 1 < ½`) and `gc ∈ [0, s/2]²`.
(iv) The roots cover `[0, s/2]² × [0, ½]` by closed boxes: 25 × 25 × 4 = **2,500 roots** for `s = 5`.

**`--full`** (no symmetry used).  Pass 0: the cover, all centres, `u ∈ [0, ½]`.  Pass 1: the cover reflected by
`y ↦ s − y`, same boxes.  The reflection maps `Q(c, θ)` to `Q(ρc, 90° − θ)`, so pass 1 covers `θ ∈ [36.87°, 90°]`.
50 × 50 × 4 × 2 = **20,000 roots**.  This needs neither the D4 check nor §7 (ii)–(iii).

## 8. Tests (`search/zmx2_tests.sh`; all pass, `runs/zmx2/tests.out`)

| # | test | result |
|---|---|---|
| T0 | `zmcheck` sources (`src/main.rs`, `Cargo.toml`, `Cargo.lock`) identical to `HEAD` | pass (so its behaviour and outputs are unchanged) |
| T1 | malformed / unsupported files: negative weight, point outside, degenerate segment, diagonal segment, polygon, trailing tokens, short file | all `ERROR`, exit 2 |
| T2 | D4: candidate invariant; one weight `+1` | refused under `--d4` (exit 2) |
| T3 | **harness** (`zmx2_tools.py sound/tight`): exact rational `μ` at random admissible rational poses (plus corners, `u = 0` faces) of random small boxes biased to germs, walls, lines (box sides `5·10⁻⁵ … 0.1`), and of the certified leaves with bound `< 1.0005` | 0 fails: 3 × 400 random boxes (candidate, reflected candidate, `s(13)`), 6,617 + 6,755 + 6,369 poses; 2 × 200 tight leaves (candidate, `s(13)`), 3,828 + 3,732 poses; `min(μ) − bound` down to `5·10⁻⁶` (0 for `s(13)`: the harness does reach tight boxes) |
| T4 | toy mixed covers (uniform / non-uniform "wave" densities on the unit grid lines, `s = 2, 3, 4`, with/without tile-centre points), scaled to `1.02 ×` and `0.98 ×` their float minimum; `--d4` and `--full` | all 8 `ok` covers verify, all 8 `rej` covers are refused at the float minimum (`θ ≈ 41–45°`) |
| T5 | candidate `× 0.994`; `× 0.9955`; line `x = 2` zeroed on `y ∈ [1.4, 1.6]` (`--full`) | refused at `(0.5736, 1.4406, 9.21°)`, exact `μ = 0.99969`; **verifies**; refused at `(1.5, 1.5, θ → 0⁺)`, exact `μ = 1.749` at `θ = 0` vs `0.9176` at `θ = 0.003°` |
| T6 | `s(13)` shipped point cover: column `c_x ∈ [1.2, 1.3]`, all `c_y`, all angles (320 roots) with `zmx2 --full` and with `zmcheck` | both clean (`zmx2` 14,644 boxes / 1 CPU-s; `zmcheck` 542 boxes, ADM 186 DISJ 152 EMPTY 93 / 18 CPU-s); `zmx2 --no-atoms` refuses at the one-cut germs `(1.2, 0.5, θ = 0)` as designed; `zmx2 --d4` whole cover `VERIFIED-D4`; `× 0.975` refused at the wall germ `(0.5, 1.5, θ → 0⁺)` |
| T7 | the candidate, `--d4` | `VERIFIED-D4` |

Point-only agreement: `zmx2` and `zmcheck` agree on the column, and on the whole `s(13)` and `s(32)` covers
(`VERIFIED-D4` both), by entirely different germ mechanisms (`zmcheck`: DISJ sign-splits and Lemma K; `zmx2`: the pair
lemma on line atoms).

## 9. Runs

`search/zmx2_run.sh` (resumable: per-root logs `runs/zmx2/*.log`, a rerun skips finished roots; header records file hash
and settings, a changed header is refused).  Manifest: `search/zmx2_manifest.txt` (git HEAD, sha256 of source, binary,
input and logs, commands, census, verdict).  Defaults: depth 40, `κ = 1.5`, `K = 30`, node cap `2·10⁷` per root.
All on `taskset -c 0-9`, 10 threads.

| run | roots | boxes | certified | empty | uncert. | max depth | CPU | verdict |
|---|---|---|---|---|---|---|---|---|
| candidate `--d4` | 2,500 | 1,826,222 | 888,945 | 25,416 | **0** | 26 | 11 s | **`VERIFIED-D4`** |
| candidate `--full` | 20,000 | 14,709,448 | 7,164,300 | 200,424 | **0** | 27 | 87 s | **`VERIFIED`** |
| candidate `× 0.9955 --d4` | 2,500 | 6,311,982 | 3,081,361 | 75,880 | 0 | 33 | 30 s | `VERIFIED-D4` |
| candidate `× 0.994 --d4` | 2,500 | 11,890,121 | 5,752,894 | 190,221 | 3,326 | 40 | 52 s | refused (true hole) |
| `s(13)` shipped, `--d4` | 1,600 | 47,162 | 22,136 | 2,245 | 0 | 20 | 3 s | `VERIFIED-D4` |
| `s(32)` shipped, `--d4 --pair-points` | 3,600 | 1,405,342 | 686,886 | 17,585 | 0 | 29 | 873 s | `VERIFIED-D4` |
| `s(32)` shipped, `--d4` (grid atoms only) | 3,600 | 1,419,262 | 693,540 | 17,853 | 38 | 40 | 87 s | refused: all 38 at `(1.5006, 0.5016, 0.179°)`, float `μ ≥ 1.0115`; the loss is the off-grid pair `(0.999, 0.999)`, `(1.999, 0.999)` on opposite edges, which `--pair-points` couples |

(The box counts include the internal nodes; certified + empty + uncertified = leaves.)

## 10. Verdict and caveats

**`s(21) = 5` confirmed by a second, independent exact checker**: the candidate's closed-square covering property
holds at every pose, with the D4 argument (2,500 roots) and without it (20,000 roots), 0 uncertified boxes.

Caveats, for the bundling step:
* The geometry enclosures use directed-rounding binary64 intervals (Lemma R), not rational arithmetic.  The claim
  rests on IEEE-754 correct rounding of `+ − × ÷ √` and on the compiler not contracting or reordering them (Rust
  guarantees both).  The point tests and all mass sums are exact integers.
* The reduction "valid mixed cover of total `< n` ⇒ `s(n) ≥ s`" for measures is the `lean` task's; `zmx2` checks only
  the covering statement.
* Harness and code were written by the same agent: the harness is an independent *implementation* (Python,
  `Fraction`s, generic segment clipping, no Lemmas C–W), not an independent *author*.
* Not supported: polygons, non-axis-parallel segments (refused, never mis-checked).
