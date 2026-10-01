# `zmx2` with area densities: a second implementation for `s(k² − 3) = k` (2026-09-30)

Brief: give the k2m3 result its second, independent implementation by extending `zmx2` (`ZMX2.md`) to area density
(Lebesgue measure on an axis-parallel square) and checking `certificates/k2m3/L4_k02_box7.txt` (the measure `μ₇`
on `[0,7]²`).  **Claim checked (Valid7):** every closed unit square `Q ⊆ [0,7]²`, at every position and angle, has
`μ₇(Q) ≥ 1`, where `μ₇` = 800 uniform segment masses on the `1/5`-grid + area density 1 on `R = [9/5, 26/5]²`.

Code: `verify2/src/bin/zmx2.rs` (new: `Rect`, the polygon parser, `mass_bound`, `area_bound`, `side_contained`,
`cap_ub`, `quad_area_lb`, `coupled_side`, `family_bound_ex`, `frames`, the θ = 0 mode `bound0`/`cmd_cert0`/`boxes0`,
`--umin`, `float_rect_area`, `fscan --outside --skip-wall --umax`), `search/zmx2_tools.py` (polygons in the exact
evaluator: Sutherland–Hodgman clipping with `Fraction`s; `sound --area`, `sound0`, `toy --kind lebesgue`,
`perturb --op rect-scale|rect-delta|seg-delta`), `search/zmx2_area_tests.sh` (§10; also run as T10 of
`search/zmx2_tests.sh`).  Logs: `search/zmx2_area/`.

## 0. Status

| part of the pose space | method | result |
|---|---|---|
| `θ = 0` (axis-parallel squares) | `zmx2 cert0` (Lemma Z0, §6; exact integers, no floating point) | **`VERIFIED-D4`** (900 roots) and **`VERIFIED`** without symmetry (3,600 roots), 0 uncertified, < 0.01 CPU-s |
| `θ ∈ [θ_min, 90°)`, `θ_min = 2 atan(2⁻¹³) = 0.013988°` | `zmx2 cert --umin 10` (ZMX2.md + Lemmas 1–3, K) | **clean** over the whole D4 region (4,900 roots, 4,392,260 boxes, 0 uncertified, 65 CPU-s) and over the whole unreduced space (§9) |
| `θ ∈ (0, θ_min)` | — | **open**: needs a first-order lemma at the *zero-limit germs* (§11); the present box-sum bound provably cannot close them |

So **Valid7 is not yet certified by `zmx2`**: what is missing is exactly the thin slab `0 < θ < 0.014°` near a few
families of germs where `μ₇(Q) → 1` as `θ → 0⁺` (walls, and the columns `c_x ∈ {1.5, 5.5}`, `c_y ∈ {1.5, 5.5}` around
the area square), §11.  Everything else is certified, including all the exactly tight places that are not θ → 0⁺
limits: squares inside `R` (`μ = 1`), squares touching `∂R` from inside at any angle (`μ = 1`), axis-parallel squares
on a wall (`μ = 1`).  D4 invariance of `μ₇` is checked exactly (`zmx2 d4`, §7).

**Independence.**  Not opened: `search/qx2_zm.py`, any `search/qx2_*.py`, `search/QUADRANT_EXACT.md`,
`search/QUADRANT.md`, `search/zm_mixed.py`, `zm_mixed_test.py`, `ZM_MIXED.md`, `ZM_MIXED_AUDIT*`,
`notes/lean-bentz-reduction.md`, `lean/Sqpack/Bentz*.lean`, `tasks/k2m3-review/`, and in `certificates/k2m3/` anything
but `README.md` and `L4_k02_box7.txt`.  **Read:** `search/ZMX2.md` (all), `search/ZMX2_AUDIT.md` (lines 1–80),
`verify2/src/bin/zmx2.rs` (all), `search/zmx2_tools.py` (all), `search/zmx2_tests.sh` (all),
`search/zmx2_census_cmp.py` (first 30 lines = all), `certificates/s21/FORMAT.md` (the mixed format; it defines the
polygon record), `certificates/s21/verify.sh` (lines 1–80), the first line of `certificates/{s21,s60}/zmx2_d4/roots.log`,
`certificates/k2m3/L4_k02_box7.txt` (header and data; parsed by my own code), and `certificates/k2m3/README.md` from
the **Claim** paragraph to the next heading — that printed range also contained its "Earlier work" paragraph, the
reference list and the line "Computer-assisted, not peer reviewed …" (no checker details; noted for completeness).
Directory listings showed the *names* of `certificates/k2m3/` and `search/` files.  Not needed and not read:
`tasks/line-cover/FORMAT.md` (not in the public tree), `search/LINE_COVER.md`, `S32_EXACT.md`.

## 1. Set-up and format

The polygon record of the mixed format (`certificates/s21/FORMAT.md`) is `k w X1 Y1 … Xk Yk`: a convex polygon,
counter-clockwise, mass `w/W` spread uniformly by area.  `zmx2` accepts polygons that are **axis-parallel rectangles**
`[x0,x1] × [y0,y1]` given counter-clockwise (any starting corner) and refuses everything else (`ERROR: unsupported`,
exit 2; also degenerate, outside, negative or `≥ 2⁵⁰` weights).  A rectangle `r` with weight `w_r` and area
`Pn_r/D²` has density `ρ_r = (w_r/W)·D²/Pn_r`.  The box file has one rectangle, `[9/5, 26/5]²`,
`w = 11.56·10¹² = W·area`, so `ρ = 1`.  Covers without polygons take exactly the old code path (`mass_bound` calls
`segment_bound`), with the same log header; covers with polygons get ` area=N` appended to the header.

Notation as `ZMX2.md` §1: pose `(c, u)`, `u = tan(θ/2)`, `C = cos θ`, `S = sin θ`, `w = C + S`, box
`B = [x0,x1] × [y0,y1] × [u0,u1]`, admissible iff `w/2 ≤ c_x, c_y ≤ s − w/2`.  Additionally `p = min(C,S)`,
`P = max(C,S)`, `m = |C − S|/2`, `M = w/2` (so `M − m = p`, `M + m = P`).  For a rectangle side, `h` is the signed
distance from the centre of `Q` to the side's line, positive when the centre is on the rectangle's side.  The
**cap** of a side is `Q ∩ H`, `H` the open half-plane beyond the side.

**Main property** (as `ZMX2.md`): the bound `L(B)` satisfies `L(B) ≤ μ(Q)` at every admissible pose of `B`.  The
area measure is a separate summand of `μ` (disjoint from the points, lines and other rectangles), so
`L(B) = P(B) + [lines] + [rectangles]` stays a lower bound if each summand is bounded below; Lemma K bounds a line
*together with* one cap, and then that line is left out of the line bound.

## 2. Exact containment (Lemma 1)

**Lemma 1.**  (a) A side with `h ≥ h₀ = hn/hd` on `B` has an empty cap at every pose of `B` if
`(2h₀+1)u² − 2u + (2h₀−1) ≥ 0` on `[u0,u1]`, decided exactly by integers (`side_contained`, Lemma P's
`qmax_le0`).  (b) A side lying on or beyond a container wall has an empty cap at every **admissible** pose.  If all
four sides of `r` are contained, the area term is `ρ_r·1` **exactly** (`TR_r = ⌊w_r Lc 2^K D²/Pn_r⌋`, `= target` for
the box file).
*Proof.*  (a) The cap is empty iff the extent of `Q` towards the side, `M = w/2`, is `≤ h`; `h ≥ h₀` and
`2h₀(1+u²) ≥ 1 + 2u − u² ⟺` the quadratic `≥ 0`, with `w = (1+2u−u²)/(1+u²)`.  (b) `Q ⊆ [0,s]²`. ∎

This is what certifies the zero-margin squares inside `R` (and the toy "Lebesgue on the container" everywhere) with
bound exactly `1`.

## 3. The cap bound (Lemma 2)

**Lemma 2.**  For `0 < θ < 90°` the area of the part of `Q` beyond a line at signed distance `h` from its centre is
`a(h, θ) = 0` (`h ≥ M`), `(M − h)²/(2pP)` (`m ≤ h ≤ M`), `½ − h/P` (`|h| ≤ m`), `1 − (M + h)²/(2pP)` (`−M ≤ h ≤ −m`),
`1` (`h ≤ −M`); at `θ = 0`, `a = clamp(½ − h, 0, 1)`.  `a` is non-increasing in `h`.  `cap_ub(h₀, B)` returns
`max` over the regions possible on `[u0,u1]` of an interval upper bound of the region's formula, which bounds the cap
of a side with `h ≥ h₀` at every pose of `B`.
*Proof.*  In coordinates centred at `c` the vertical chords of `Q` have length `L(x) = 1/P` for `|x| ≤ m`,
`(M − |x|)/(pP)` for `m ≤ |x| ≤ M`, `0` beyond (the vertices of `Q` have abscissae `±m`, `±M`; between them the
chord is bounded by two sides of slopes `−C/S·` and `S/C·`, reaching `1/P` at `|x| = m`).  `a(h) = ∫_{−M}^{−h} L`:
for `h ≥ m` it is `∫_h^M (M−x)/(pP) = (M−h)²/(2pP)`; for `|h| ≤ m` it is `p²/(2pP) + (m−h)/P = ½ − h/P`; for
`h ≤ −m` it is `1 − a(−h)` by central symmetry.  Monotone in `h` since `L ≥ 0`.  Over `B` the actual `h` is `≥ h₀`,
so the cap is `≤ a(h₀, θ)`; for `θ ∈ [θ0,θ1]` the point `(h₀, θ)` lies in one of the regions, and a region is
possible only if `h₀ ≥ m_lo` (ii), `|h₀| ≤ m_hi` (iii), `h₀ ≤ −m_lo` (iv) with `m` enclosed through
`C − S = (1 − 2u − u²)/(1+u²)` (decreasing).  The bounds used: (ii) `s = M − h ∈ [0, p]`, `a ≤ s/(2P)` and
`a ≤ s²/(2pP)` with `s ≤ min(M_hi − h₀, p_hi)`, `P ≥ P_lo`, `p ≥ p_lo`; (iii) `½ − h₀/P_hi` (`h₀ ≥ 0`) or
`½ + |h₀|/P_lo`; (iv) `1 − q²/(2 p_hi P_hi)`, `q = max(0, M_lo + h₀) ≤ M + h` (covers (v) with `q = 0`).  `M_hi`
from `w_hi` (`w` is increasing on `[0, √2−1]`, decreasing after, max `√2`), `M_lo = w_lo/2`. ∎

## 4. Corners (Lemma 3) and the area bound

**Lemma 3 (inclusion–exclusion).**  `area(Q ∩ R) = 1 − Σ_{4 sides} area(cap) + Σ_{4 corners} area(Q ∩ K_j)`, `K_j`
the open quadrant beyond corner `j` (e.g. `{x < x0_R, y < y0_R}`).
*Proof.*  `Q \ R = Q ∩ ⋃ H_i`; of the pairwise intersections of the four half-planes the two opposite pairs are empty
and the four adjacent ones are the quadrants; every triple contains an opposite pair. ∎

Lower bounds for `area(Q ∩ K_j)` over `B` (the larger is used): **(a)** `Q ⊇ c + [−r, r]²`, `r = 1/(2w)`
(`|±rC ± rS| ≤ r w = ½`), so `area ≥ λ_x λ_y` with `λ_x = clamp(x0_R − x1 + r, 0, 2r)` etc. (the overlap of the small
square decreases as `c` moves away from `K_j`), `r ≥ 1/(2 w_hi)`.  **(b)** `area(Q ∩ K_j)` is non-increasing as
`c` moves away from `K_j` (`area((Q+v) ∩ K) = area(Q ∩ (K − v))` and `K − v ⊆ K`), so it is `≥` its value at the
farthest corner `c*` of the centre rectangle; at `c*` and the mid angle `u_m` it is computed by clipping `Q` with
outward-rounded interval arithmetic (`quad_area_lb`: every vertex must be classified with certainty, else no bound);
and `|d/dθ area(Q_θ ∩ K)| ≤ ∫_{∂Q} |v·n| ≤ (√2/2)·4 = 2√2` (Reynolds; `|v| ≤ √2/2` for rotation about `c`), with
`|θ − θ_m| ≤ u1 − u0`, so subtracting `2√2 (u1 − u0)` gives a bound for all of `B`.

**Area term.**  With `caps = Σ cap_ub` (rounded up) and `corners = Σ` corner bounds (rounded down), the rectangle
contributes `TR − ⌈caps·(TR+1)⌉ + ⌊corners·TR⌋` bound units (`≤ ρ·(1 − caps + corners)·target` also when negative,
because `TR ≤ ρ·target < TR + 1`), clipped at 0 in config 0 below; exactly `TR` when every side is contained.

## 5. A side coupled with its line (Lemma K) and the assembly

Inside `R` near a side the margin is zero: a rotated square inside `R` that touches the side at a corner has
`μ = 1` exactly, and the poses where the corner pokes out lose `O(ε²)` area while the side's segments gain `O(ε)`.
Separately pessimised (inf of the area + inf of the line) a box around such a pose never reaches 1; the two must be
coupled.

**Lemma K.**  Let the side lie on the line `λ` (a vertical line of the frame; horizontal sides via Lemma T) and let
the centre of `Q` be on the rectangle's side of `λ` at every pose of `B` (`h ≥ 0`).  Let `ℓ` be the length of the chord
`Q ∩ λ`, `t = M − h` the depth of the cap.  **(a)** `area(cap) = ℓ·φ(t, p)` with `φ = t − min(t, p)/2` (for `t > 0`),
hence `area(cap) ≤ ℓ·φ_max`, `φ_max = t_max − min(t_max, p_min)/2`, `t_max = w_hi/2 − h₀`,
`p_min = min(S(u0), C(u1))`.  **(b)** With `G(y) = F_λ(y) − κ y`, `κ = ⌈w_r Lc φn/(Pn_r 2^K)⌉` (bound units per
grid unit; `φn ≥ φ_max·G` on the grid), at every pose of `B`
`ν_λ(Q) − ρ_r·area(cap) ≥ G(hi) − G(lo) (+ certain atoms)`, where `[lo, hi]` is the actual chord, and the infimum of
`G(hi) − G(lo)` over `lo ∈ I_lo`, `hi ∈ I_hi`, `lo ≤ hi` is computed exactly.  **(c)** The value used is that infimum
if the chord is certainly non-empty on `B`, `min(0, ·)` otherwise, `0` if no `lo ≤ hi` is possible.
*Proof.*  (a) Case `h ≥ m` (`t ≤ p`): the cap is the triangle cut off at a vertex; `ℓ = L(h) = t/(pP)` and
`a = t²/(2pP) = ℓt/2`.  Case `0 ≤ h ≤ m` (`t ≥ p`): `ℓ = 1/P`, `a = ½ − h/P = ℓ(P/2 − h) = ℓ(t − p/2)` (`t = M − h`,
`M − P/2 = p/2`).  At `θ = 0` (`p = 0`) the cap is a `1 × t` rectangle, `φ = t`.  `φ` is non-decreasing in `t` and
non-increasing in `p`.  (b) `ν_λ(Q) ≥ F(hi) − F(lo)` (+ atoms in `Q`) and `ρ_r·ℓ φ_max ≤ κ (hi − lo)` in bound
units; `G` is piecewise linear with kinks at the breakpoints of `λ`.  The ranges: `lo = c + max(f1, g1)`,
`hi = c + min(f2, g2)` (Lemma C) are enclosed by Lemma M (`f, g` affine in `d`, limits at `u0 = 0` included) and
clipped to `[y0 − M_hi, y1 + M_hi]` (a non-empty chord lies in `Q`).  Exact infimum: scan `hi` over the sorted
candidates (range ends and breakpoints) keeping the running maximum of `G` over the `lo`-candidates `≤ hi` (`y = hi`
included); between consecutive candidates `G(y') − max(K, G(y'))` is the minimum of a linear function and `0`, hence
concave, so its minimum is at a candidate.  (c) Lemma C covers `u > 0`; the chord is certainly non-empty for those
poses if `lo_max < hi_min`.  At `θ = 0` (`u0 = 0`): if `|d| < ½` the chord `[c − ½, c + ½]` is the `u → 0⁺` limit
chord, which lies in the closed ranges, and the cap is `ℓ t` exactly; if `|d| = ½` the cap is empty and the chord
contains the limit chord, so the value is `≥ F(limit chord) ≥ G`-difference of it; if `|d| > ½` (possible only if
`2 hn_max > den`, flag `far0`) the chord is empty and the true value is `0`, whence `min(0, ·)`. ∎

**Assembly (`area_bound`).**  *Config 0:* `Σ_r max(0, area term)` + the line bound of `ZMX2.md` (`segment_bound`,
with `--sym-atoms` if given).  If that reaches the target it is used.  *Config 1:* every side that is not exactly
contained and has `h ≥ 0` on `B` is coupled (one coupling per line position), its cap leaves the area term
(no clipping), the coupled lines leave the line DP (`family_bound_ex`), and the Lemma K values are added; with
`--sym-atoms` both atom assignments are tried.  `L = max(config 0, config 1)`.  Every term is a lower bound of a
disjoint part of `μ` (Lemma 3 is an identity, Lemma K bounds `ν_λ` minus a part of the area deficit), so `L` has the
Main property.

**Arithmetic.**  As `ZMX2.md` §5: geometry by outward-rounded binary64 intervals, then onto the grid `1/(D 2^K)`;
masses exact `i128`.  New sizes: `TR`, `κ` and `G` use `w_r·Lc < 2^70` (checked at start-up, else `ERROR`),
`φn < 2^(K+1) D`; `κ·y` and `w_r Lc φn` are checked multiplications (a panic, never a wrap, on overflow).  For the
box file `w Lc = 1.156·10¹³ < 2^44`, `target = 10¹² 2^30 < 2^70`.

## 6. `θ = 0` exactly (Lemma Z0, `zmx2 cert0`)

At `θ = 0`, `Q(c) = [c_x ± ½] × [c_y ± ½]` (closed).  For a centre box `[x0,x1] × [y0,y1]` (no angle):

**Lemma Z0.**  `μ(Q(c)) ≥ P₀ + min_{(x,y) ∈ X × Y} [ H(x) + V(y) + Σ_r ρ_r ox_r(x) oy_r(y) ]` for every `c` in the
box, where `P₀` = the points (and line atoms) in `Q(c)` for every `c` of the box; `V(y) = Σ` over the vertical lines
with `|c_x − ℓ| ≤ ½` on the whole box of `F_ℓ(y + ½) − F_ℓ(y − ½)`; `H` likewise for horizontal lines;
`ox_r(x) = |[x ± ½] ∩ [x0_r, x1_r]|`; `X`, `Y` the candidate grids (box ends, breakpoints `± ½`, rectangle sides `± ½`).
*Proof.*  Lines in `Q` for all `c` of the box carry `F(c ± ½)` exactly (closed chord); the others contribute `≥ 0`.
`V`, `H`, `ox`, `oy` are piecewise linear with kinks in `Y` resp. `X`, so on each cell of the grid the bracket is
bilinear and attains its minimum at a vertex. ∎

Everything is exact integer arithmetic (box ends are grid points for `D = 5`; otherwise rounded outward; the area mass
is `⌊w Lc ox oy/(Pn 2^K)⌋`).  Unlike the per-line infima of the angled bound, the joint minimum is exact where several
lines' windows move in opposite directions (e.g. at `(½, 2.5⁻)`: the line `x = 0.8` loses while `x = 1` gains).
Roots: the `1/10` cells of `[½, s − ½]²` (`--full`) or `[½, s/2]²` (`--d4`: D4 maps axis-parallel squares to
axis-parallel squares).  `zmx2 boxes0` is the batch interface for the harness (`zmx2_tools.py sound0`).

## 7. D4 and the reflected cover

`check_d4` additionally requires the multiset of weighted rectangles to be invariant under `x ↦ s − x` and
`x ↔ y` (sufficient for invariance of the area measure); `reflect_y` (pass 1 of `--full`) reflects rectangles.  The
box file: `D4: measure invariant under x->s-x and x<->y (exact)`; a non-symmetric rectangle and a single segment
`+1` are refused (test A2).

## 8. `--umin m`

`zmx2 cert … --umin m` starts root bin 0 at `u = 2⁻ᵐ/8` instead of 0 (log header gets ` umin=2^-m/8`; the run never
prints a verdict, only `REGION CLEAN`).  It isolates the slab of §11.

## 9. Runs (cores 10–11, 2 threads; logs in `search/zmx2_area/`)

| run | roots | boxes | certified | empty | uncert. | max depth | CPU | result |
|---|---|---|---|---|---|---|---|---|
| `cert0 --d4` (θ = 0) | 900 | 900 | 900 | — | 0 | 0 | < 0.01 s | **`VERIFIED-D4 (theta = 0)`** |
| `cert0 --full` (θ = 0) | 3,600 | 3,600 | 3,600 | — | 0 | 0 | < 0.01 s | **`VERIFIED (theta = 0)`** |
| `cert --d4 --umin 10` (θ ≥ 0.014°) | 4,900 | 4,392,260 | 1,974,266 | 224,314 | **0** | 40 | 65 s | **clean** |
| `cert --full --umin 10` | see §9.1 | | | | | | | |
| `cert --d4 --umin 14` / `--umin 20` | 4,900 | (stopped at 3 per root) | | | 1,893 / 1,104 | 40 / 23 | 27 s / 23 s | refused only in the germ columns of §11 (depth / `ul ≤ 27` limits) |
| `cert --d4` (all θ, bin 0 from `u = 0`), region `y ≥ 1.8` | 2,380 | 464,992 | 231,310 | 1,682 | 1,282 | 40 | 10 s | refused only in the cells of §11, all at `u0 = 0` |

### 9.1 (filled in below by the full-region runs)

## 10. Tests (`search/zmx2_area_tests.sh`, 40/40; also T10 of `search/zmx2_tests.sh`)

| # | test | result |
|---|---|---|
| A1 | parser: two ccw rectangles accepted; clockwise, degenerate, triangle, non-rectangular quadrilateral, weight `2⁵⁰`, outside, negative weight refused | pass |
| A2 | D4 with a rectangle: box file invariant; one segment `+1` refused; a non-symmetric rectangle refused | pass |
| A3 | toy: area density 1 on the whole container, `s = 2, 3`: `μ = 1` at **every** pose; `cert --d4/--full` and `cert0` `VERIFIED` (all roots at depth 0: Lemma 1(b) makes the bound exactly 1); density `− 10⁻⁶`: refused by `cert` and `cert0`, exact `μ = 0.999999` at the uncertified boxes | pass |
| A4 | the box file at θ = 0: `cert0 --d4` and `--full` `VERIFIED` | pass |
| A5 | rejection on the box file: area density `− 10⁻⁶` (refused by `cert0` and by `cert` around `(3.1, 3.1)`, exact `μ = 0.999999`); the wall-band piece `x = 0.8, y ∈ [2.8, 3.0]` `− 10⁻⁴` (refused by `cert0`, exact `μ(½, 3, 0) = 0.9999`); the piece `x = 1.6, y ∈ [3.2, 3.4]` `− 10⁻⁴` (refused by `cert0` and by `cert --umin 10` near `(1.5, 3.3)`; exact `μ = 0.99994` at `(1.5 + θ/2, 3.3 + 0.28θ, θ = 0.014°)`) | pass |
| A6 | harness (exact rational `μ`, polygon clipping with `Fraction`s, at random admissible rational poses and box corners): `sound0` (θ = 0, boxes biased to walls, rectangle sides and half-integer line offsets) on the box file and a weakened one; `sound --area` (boxes at the area square's sides, corners, inside-touching poses; plain and reflected); generic `sound`; weakened cover; the tightest certified leaves (bound `< 1.0005`) of a run around the left side of `R` | 0 FAIL in all (≈ 30,000 poses); `min(μ) − bound` = 0 on many boxes (the harness reaches the exactly tight places) |
| A7 | no area densities: `s(21) --d4` census = shipped `certificates/s21/zmx2_d4/roots.log` root for root, header unchanged | pass |

## 11. What remains: the zero-limit germs (`0 < θ < θ_min`)

**The obstruction.**  `zmx2` bounds a box by a *sum of separately pessimised pieces* (points, each line or line
pair, each rectangle).  At a pose where `μ = 1` is attained only as a limit `θ → 0⁺`, and where some piece decreases
to first order along the way while others increase, every box containing the limit has
`Σ inf(piece) = 1 − (first-order loss) < 1`; refining never helps.  Lemma Z removes this for the cut variable `ζ` of
two lines at distance 1, Lemma K for a side and its cap, Lemma Z0 at `θ = 0` (joint minimum); none removes it for the
angle.

**Where it bites in `μ₇`** (D4 region; found by the runs of §9 and by exact germ scans `μ(c₀ + (a,b)θ, θ)` at
`θ = 2·10⁻⁶`):
* the **walls**: `c_x = ½`, `c_y ∈ [2.5, 4.5]` (and D4 images).  At `θ = 0`, `μ = 1` on the whole band (between
  grid lines; `cert0` certifies it); for `θ > 0` against the wall `μ ≈ 1 + 0.5θ·`… — the line `x = 1` loses its lower
  end at rate `≈ ρ·u` while the horizontal lines gain beyond `x = 1`, because admissibility pushes `c_x ≥ w/2`;
  germ minima `(μ − 1)/θ ≈ 0.21–0.27`;
* the **columns** `c_x = 1.5`, `c_y ∈ [2.2, 3.5]` (Q between `x = 1` and `x = 2`, `R` covering its right fifth), and
  `(1.5, 1.5)`: the infimum `1` is approached along `c_x = 1.5 + θ/2` (the corner of `Q` touching `x = 1`), germ
  minima `(μ − 1)/θ ≈ 0.16–0.22`.

In all of them the germ value is exactly 1 and the first-order term in `θ` is positive, so `μ₇ ≥ 1` there (numerically),
but `zmx2` cannot certify a neighbourhood of `θ = 0⁺`.  With `--umin m` the B&B certifies down to `θ_min`
(`m = 10`: 65 CPU-s for the D4 region); smaller `θ_min` costs more (`m = 14` needs depth > 40) and never reaches 0.

**What is needed (design, not implemented): Lemma U (first order at `θ = 0`).**  For a box with `u0 = 0`, bound every
piece by `A_i + u·S_i` valid on the box, with `A_i` the `θ → 0⁺` limit and `S_i` a lower bound of its `u`-derivative
along admissible paths (`c` fixed, or `c_x(t) = w(t)/2 + e` against a wall), and certify if
`Σ A_i ≥ 1` (jointly minimised over `c`, as in Lemma Z0) and `Σ S_i ≥ 0`:
* interior lines: chord ends `c ± ½ + (affine in t)` from `T ∈ [2t, 2t/(1−u1²)]`, `R/2 − ½ ≤ t u1/(1−u1²)`, `q ≥ ½ − u1 t`
  (Lemma W), `(w − 1)/2 ∈ [t(1−u1)/(1+u1²), t]`; slopes times density bounds on the swept ranges;
* line pairs at germs (Lemma Z): `inf_ζ` of the pair at `t` is `≥` its value at `t = 0` minus the losses of the
  four outer ends (monotonicity + Lipschitz; the cut moving with `t` only gains);
* rectangles: `d area/dθ = Σ_edges ∫_{edge ∩ R} (−s) ds` (`s` along the edge from its midpoint) — lower-bounded by
  counting possibly-inside `s > 0` and certainly-inside `s < 0`;
* `A = Σ A_i` exact (Lemma Z0-type joint minimum over the base rectangle; pairs by Lemma Z at `u0 = 0`).
The numbers above say the slope budget (`≥ 0.16` per unit `θ`) is ample for the slack of such bounds.  Cost: only
the `u0 = 0` boxes near the listed germs (a few hundred roots) would use it.

## 12. Caveats

* The new geometry (cap bounds, corner clipping, chord ranges in Lemma K) uses outward-rounded binary64 intervals
  like the rest of `zmx2`; `cert0` is pure integer arithmetic.
* The harness and the code were written by the same agent (an independent *implementation* of `μ`, not an
  independent author), as in `ZMX2.md`.
* Lemma K couples a line with only one cap; a line that is the side of two rectangles is coupled once (the other
  cap stays in the union bound).  Only axis-parallel rectangles are supported.
