# `zmx2` with area densities: a second implementation for `s(k² − 3) = k` (2026-09-30)

Brief: give the k2m3 result its second, independent implementation by extending `zmx2` (`ZMX2.md`) to area density
(Lebesgue measure on an axis-parallel square) and checking `certificates/k2m3/L4_k02_box7.txt` (the measure `μ₇`
on `[0,7]²`).  **Claim checked (Valid7):** every closed unit square `Q ⊆ [0,7]²`, at every position and angle, has
`μ₇(Q) ≥ 1`, where `μ₇` = 800 uniform segment masses on the `1/5`-grid + area density 1 on `R = [9/5, 26/5]²`.

Code: `verify2/src/bin/zmx2.rs` (new: `Rect`, the polygon parser, `mass_bound`, `area_bound`, `side_contained`,
`cap_ub`, `quad_area_lb`, `coupled_side`, `family_bound_ex`, `frames`, the θ = 0 mode `bound0`/`cmd_cert0`/`boxes0`,
`--umin`, the first-order bound `--first-order` (`first_order_bound`, `first_order_sig`, `u_line`, `u_slope`,
`end_rate`, `one_bp`, `dens_range`, `u_pair_cell`, `pair_sweep`, `u_area_slope`, `UTerm`, `UPair`, `UVar`, the cone
split of Lemma V: `AxisMode`, `AxisSpec`, `shift_rates`, `germ_splits`), `float_rect_area`,
`fscan --outside --skip-wall --umax`), `search/zmx2_tools.py` (polygons in the exact evaluator: Sutherland–Hodgman
clipping with `Fraction`s; `sound --area --small-u`, `sound --near X,Y,U --r R --ru RU --fine`, `sound0`, `toy --kind lebesgue`,
`perturb --op rect-scale|rect-delta|seg-delta`, `ZMX2` from the environment), `search/zmx2_area_tests.sh` (§10; also run as T10 of
`search/zmx2_tests.sh`).  Logs: `search/zmx2_area/`.

## 0. Status

| part of the pose space | method | result |
|---|---|---|
| `θ = 0` (axis-parallel squares) | `zmx2 cert0` (Lemma Z0, §6; exact integers, no floating point) | **`VERIFIED-D4`** (900 roots) and **`VERIFIED`** without symmetry (3,600 roots), 0 uncertified, < 0.01 CPU-s |
| `θ ∈ [θ_min, 90°)`, `θ_min = 2 atan(2⁻¹³) = 0.013988°` | `zmx2 cert --umin 10` (ZMX2.md + Lemmas 1–3, K) | **clean** over the D4 region (4,900 roots, 4,392,260 boxes, 0 uncertified, 65 CPU-s) **and over the whole unreduced pose space** (`--full`, 39,200 roots, 35,123,736 boxes, 0 uncertified, 440 CPU-s) |
| all `θ` (`0 ≤ θ < 90°`) | `zmx2 cert --first-order` (Lemma U per base cell + Lemma V, §11.5–11.7) | **`VERIFIED-D4`** (4,900 roots, 1,158,478 boxes, 0 uncertified, 53 CPU-s) **and `VERIFIED`** without symmetry (`--full`, 39,200 roots, 9,286,868 boxes, 0 uncertified, 399 CPU-s) |

So **`zmx2` certifies Valid7** (2026-10-01): the single run `cert --full --first-order` covers every pose, at every
angle including `θ = 0` and `0 < θ < θ_min`, without using symmetry; `cert0` and `cert --umin 10` certify the
`θ = 0` and `θ ≥ 0.014°` parts once more by other means.  Everything that is exactly tight is certified: squares
inside `R` (`μ = 1`, Lemma 1), squares touching `∂R` from inside at any angle (`μ = 1`, Lemma K; at the corner of `R`
with Lemma K (d), §11.7), axis-parallel squares on a wall (`μ = 1`, Lemma Z0), and the `θ → 0⁺` limits at the walls,
the column germs and the **double germs** (`μ → 1`, Lemma U per base cell and Lemma V; the double germs were the
last gap, §11.4).  D4 invariance of `μ₇` is checked exactly (`zmx2 d4`, §7).  Lemma U/V are the newest and least
reviewed part; §13 is a self-audit (one proof gap found in the old Lemma U and closed; no counterexample).

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
`tasks/line-cover/FORMAT.md` (not in the public tree), `search/LINE_COVER.md`, `S32_EXACT.md`.  The 2026-10-01
session (§11.4–11.8, §13; a fresh agent) read only this file, `verify2/src/bin/zmx2.rs`, `search/zmx2_tools.py`
(the parts it changed or calls), `search/zmx2_area_tests.sh`, `search/zmx2_tests.sh` and the logs in
`search/zmx2_area/`, and parsed `certificates/k2m3/L4_k02_box7.txt` with its own scripts; it opened none of the files
listed above and grepped only the files it read.

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
**(d)** (2026-10-01, §11.7) The ranges `I_lo`, `I_hi` are also clipped exactly at `y0 − ½` when `d ≤ −u1/(2(1 − u1²))` on the box, resp. at
`y1 + ½` when `d ≥ u1/(2(1 − u1²))`.

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

## 9. Runs (2026-09-30: cores 10–11, 2 threads; 2026-10-01 (`_v2`): cores 8–15, 8 threads, shared machine; logs in `search/zmx2_area/`)

| run | roots | boxes | certified | empty | uncert. | max depth | CPU | result |
|---|---|---|---|---|---|---|---|---|
| `cert0 --d4` (θ = 0) | 900 | 900 | 900 | — | 0 | 0 | < 0.01 s | **`VERIFIED-D4 (theta = 0)`** |
| `cert0 --full` (θ = 0) | 3,600 | 3,600 | 3,600 | — | 0 | 0 | < 0.01 s | **`VERIFIED (theta = 0)`** |
| `cert --d4 --umin 10` (θ ≥ 0.014°) | 4,900 | 4,392,260 | 1,974,266 | 224,314 | **0** | 40 | 65 s | **clean** |
| `cert --full --umin 10` (θ ≥ 0.014°, both passes, no symmetry) | 39,200 | 35,123,736 | 15,788,204 | 1,793,264 | **0** | 40 | 440 s | **clean** |
| `cert --d4 --first-order` (all θ), 2026-09-30 code | 4,900 | 1,179,242 | 587,463 | 4,488 | 120 | 40 | 96 s | refused: only the double germs of §11.4 |
| `cert --full --first-order` (all θ), 2026-09-30 code | 39,200 | 9,455,676 | 4,710,212 | 36,078 | 1,108 | 40 | 797 s | refused: only their D4 images |
| **`cert --d4 --first-order` (all θ), final code** (`k7_d4_first_order_v2`) | 4,900 | 1,158,478 | 577,204 | 4,485 | **0** | 23 | 53 s (37 s on an idle run) | **`VERIFIED-D4`** |
| **`cert --full --first-order` (all θ, both passes, no symmetry), final code** (`k7_full_first_order_v2`) | 39,200 | 9,286,868 | 4,627,012 | 36,022 | **0** | 23 | 399 s (332 s idle) | **`VERIFIED`** |
| `cert --d4 --umin 10`, final code (`k7_d4_umin10_v2`) | 4,900 | 4,392,260 | 1,974,266 | 224,314 | 0 | 40 | 71 s | clean (same census as 2026-09-30) |
| `cert --full --umin 10`, final code (`k7_full_umin10_v2`) | 39,200 | 35,123,736 | 15,788,204 | 1,793,264 | 0 | 40 | 604 s | clean (same census) |
| `cert --d4` (all θ, no Lemma U), region `c_y ≥ 1.8` | 2,380 | 464,992 | 231,310 | 1,682 | 1,282 | 40 | 10 s | refused at every zero-limit germ (walls, columns) |
| `cert --d4 --umin 14` / `--umin 20` | 4,900 | (3 per root) | | | 1,893 / 1,104 | 40 / 23 | 27 / 23 s | refused only in the germ columns (depth / `ul ≤ 27`) |

(Refused runs stop a root after 200 uncertified boxes, so their counts are lower bounds.  CPU times of the `_v2`
runs vary by ≈ 1.5× with the load of the shared machine; the box censuses are deterministic.)

## 10. Tests (`search/zmx2_area_tests.sh`, 80/80; also T10 of `search/zmx2_tests.sh`)

| # | test | result |
|---|---|---|
| A1 | parser: two ccw rectangles accepted; clockwise, degenerate, triangle, non-rectangular quadrilateral, weight `2⁵⁰`, outside, negative weight refused | pass |
| A2 | D4 with a rectangle: box file invariant; one segment `+1` refused; a non-symmetric rectangle refused | pass |
| A3 | toy: area density 1 on the whole container, `s = 2, 3`: `μ = 1` at **every** pose; `cert --d4/--full` and `cert0` `VERIFIED` (all roots at depth 0: Lemma 1(b) makes the bound exactly 1); density `− 10⁻⁶`: refused by `cert` and `cert0`, exact `μ = 0.999999` at the uncertified boxes | pass |
| A4 | the box file at θ = 0: `cert0 --d4` and `--full` `VERIFIED` | pass |
| A5 | rejection on the box file: area density `− 10⁻⁶` (refused by `cert0` and by `cert` around `(3.1, 3.1)`, exact `μ = 0.999999`); the wall-band piece `x = 0.8, y ∈ [2.8, 3.0]` `− 10⁻⁴` (refused by `cert0`, exact `μ(½, 3, 0) = 0.9999`); the piece `x = 1.6, y ∈ [3.2, 3.4]` `− 10⁻⁴` (refused by `cert0` and by `cert --umin 10` near `(1.5, 3.3)`; exact `μ = 0.99994` at `(1.5 + θ/2, 3.3 + 0.28θ, θ = 0.014°)`) | pass |
| A6 | harness (exact rational `μ`, polygon clipping with `Fraction`s, at random admissible rational poses and box corners): `sound0` (θ = 0, boxes biased to walls, rectangle sides and half-integer line offsets) on the box file and a weakened one; `sound --area` (boxes at the area square's sides, corners, inside-touching poses; plain and reflected); generic `sound`; weakened cover; the tightest certified leaves (bound `< 1.0005`) of a run around the left side of `R` | 0 FAIL in all (≈ 30,000 poses); `min(μ) − bound` = 0 on many boxes (the harness reaches the exactly tight places) |
| A8 | `--first-order`: wall band and a single-germ column cell clean down to θ = 0 (and refused without the flag); rejection: the wall's first-order gain removed (h-lines `y = 2.6…3.4` zeroed on `x ∈ [1, 1.2]`: `μ = 1` on the wall at θ = 0 but `0.99999973` at `(0.50001, 3, u = 10⁻⁵)`) refused; the column piece of A5 refused; harness with `--small-u` (poses at angles `u1 2⁻ᵏ`, `k ≤ 40`, to reach the germ variables) near 8 germs, on the gain-removed cover, area boxes, reflected cover | pass, 0 FAIL (≈ 25,000 poses; `min(μ) − bound` down to `8·10⁻¹²` on the weakened cover) |
| A9 | (2026-10-01) double germs: `--d4 --first-order` on the column `c_x ∈ [1.5, 1.6]`, `c_y ∈ [1.5, 3.5]`, on the cell `(2.3, 2.3)`, and `--full` on `(1.5, 5.5)` clean down to `θ = 0`; fine harness (boxes of side `1/(10·2^a)`, `a = 12..18`, within `3·10⁻⁵` of the germ, tiny angles) near all 8 double germs, plain and reflected cover; **rejection 1** (pair cut): `+2·10⁻⁵` on `y = 3.6, x ∈ [1, 1.2]` and on `y = 2.6, x ∈ [1.6, 1.8]`, `−10⁻⁵` on `x = 1.6, y ∈ [3.2, 3.4]`: `cert0 --full` `VERIFIED`, `--umin 10` clean around the germ, but exact `μ(1.5000020011, 3.1000004, u = 2·10⁻⁶) = 0.99999065`: refused by `--first-order` at `(1.5, 3.1)`; **rejection 2** (corner of `R`): `y = 1.8` zeroed on `x ∈ [1.8, 2]`: `cert0` `VERIFIED`, exact `μ(2.3001, 2.30008, u = 10⁻⁴) = 0.999999001`, refused at `(2.3, 2.3)`; fine harness on both weakened covers | pass, 0 FAIL (≈ 120,000 poses; `min(μ) − bound = 0` at many boxes) |
| A7 | no area densities: `s(21) --d4` census = shipped `certificates/s21/zmx2_d4/roots.log` root for root, header unchanged | pass |

Also, with the final source (cores 10–11, `NPROC=2`): `certificates/s21/verify.sh` and `certificates/s60/verify.sh`
(fast tier: they rebuild `zmx2` from source and compare a fresh `--full` census with the shipped one root for root):
**`s(21) bundle: OK`**, **`s(60) bundle: OK`**; `search/zmx2_tests.sh`: **63/63** (T0–T9 unchanged, T10 = this
suite).  Record: `search/zmx2_area/verify_and_tests.out`.  2026-10-01 (final code, cores 8–15, `TH=8`):
`search/zmx2_tests.sh` **64/64** (T10 = this suite, 80/80), 12.8 min wall; records
`search/zmx2_area/verify_and_tests_2026-10-01.out`, `area_tests_2026-10-01.out`.  The s(21)/s(60) bundle scripts were
not re-run (the changed code is reached only with area densities or `--first-order`; A7 re-checks the s(21) census).

## 11. The zero-limit germs and Lemma U (`--first-order`)

### 11.1 The obstruction

`zmx2` bounds a box by a *sum of separately pessimised pieces* (points, each line or line pair, each rectangle).  At
a pose where `μ = 1` is attained only as a limit `θ → 0⁺`, and where some piece decreases to first order along the
way while others increase, every box containing the limit has `Σ inf(piece) = 1 − (first-order loss) < 1`; refining
never helps.  Lemma Z removes this for the cut variable `ζ` of two lines at distance 1, Lemma K for a side and its
cap, Lemma Z0 at `θ = 0` (joint minimum); none of them removes it for the angle.  In `μ₇` (D4 region; runs of §9 and
exact germ scans `μ(c₀ + (a,b)θ, θ)`, `θ = 2·10⁻⁶`):

* the **walls**: `c_x = ½`, `c_y ∈ [2.5, 4.5]`.  At `θ = 0`, `μ = 1` on the whole band; for `θ > 0` against the wall
  the line `x = 1` loses its lower end at rate `≈ ρ u` while the horizontal lines gain beyond `x = 1` because
  admissibility pushes `c_x ≥ w/2`; germ minima `(μ − 1)/θ ≈ 0.21–0.27`;
* the **column** `c_x = 1.5`, `c_y ∈ [2.2, 3.5]` and `(1.5, 1.5)` (`Q` between `x = 1` and `x = 2`, `R` covering its
  right fifth): the infimum `1` is approached along `c_x = 1.5 + θ/2` (the corner of `Q` touching `x = 1`); germ
  minima `(μ − 1)/θ ≈ 0.16–0.22`.

### 11.2 Lemma U

(As written on 2026-09-30.  Since 2026-10-01 the code implements (b′), (c) and (f) per base cell, with the
corrections also on the pair ends and per-cell pair forms, §11.6; and adds the cone split of Lemma V, §11.5.)

Let `B` have `u1 ≤ 1/16` (θ ≤ 7.2°, so `p = S`, `P = C`); Lemma U bounds the box `B'` = `B` with `u0` replaced by 0
(its poses contain those of `B`).  Choose path signs `σ_x, σ_y ∈ {−1, 0, 1}` (`+1` forced if `B` reaches
`c_x < w_hi/2 + …` at the low wall, `−1` at the high wall, otherwise all three are tried) and write every pose as
`c = b + σ (w(u) − 1)/2` per axis with a **base** `b` (at a wall `b` is admissible at θ = 0 since `c ≥ w/2`).  Using
`w − 1 = 2u(1 − u)/(1 + u²) ∈ [2u(1 − u1)/(1 + u1²), 2u]`, `T = tan θ ∈ [2u, 2u/(1 − u1²)]`,
`R/2 − ½ = u²/(1 − u²) ≥ 0`:

**(a) Interior lines.**  If on `B'` the `f`-ends of Lemma C are dominated (`f1 ≤ −½ + min(0, α_g u1)`,
`f2 ≥ ½ + max(0, β_g u1)`, Lemma M enclosures with `u → 0` limits), or the line is the Lemma W line of a wall, then at
every admissible pose the chord `[lo, hi]` satisfies `lo ≤ b − ½ + α u`, `hi ≥ b + ½ + β u` with
`α = α_g + σ⁺`, `β = β_g + σ⁻` (`α_g = −2 d_min` or `2|d_min|/(1−u1²)`, `β_g = 2|d_max|` or `−2 d_max/(1−u1²)`, the
Lemma W variants `max(α_g, u1)`, `min(β_g, −u1)` from `q ≥ ½ − u1 u`, and `σ⁺ = 1`, `σ⁻ = (1−u1)/(1+u1²)` for
`σ = +1`, `σ⁺ = −(1−u1)/(1+u1²)`, `σ⁻ = −1` for `σ = −1`).  (The `σ` part only shifts the along coordinate.)

**(b) Window slope.**  Then `F(hi) − F(lo) ≥ W(b) + u s − corr(b)` with `W(b) = F(b + ½) − F(b − ½)`, `s` = `β` times
the smallest (gain) or largest (loss) density on the range the hi end sweeps for bases in the box, minus the same for
the lo end, and **(b')** if an end moves the adverse way across exactly one breakpoint `p`, the density beyond `p` in
`s` and `corr(b) = (ρ_before − ρ_beyond)⁺ · min(dist(end(b), p)⁺, δ_max)`, by
`∫_e^{e+δ} ρ ≤ ρ_beyond δ + (ρ_before − ρ_beyond)⁺ min((p − e)⁺, δ)` (both modes are tried per family).

**(c) Germ pairs** (Lemma Z, `a` at `ℓ` with its lo end cut by `ζ`, `b` at `ℓ + 1` with its hi end cut by `ζ + u`):
with base ends `H = b + ½`, `L = b − ½` and `X_a(ζ) = F_a(H) − F_a(max(ζ, L))`, `X_b(ζ) = F_b(min(ζ, H)) − F_b(L)`,
for every `ζ` the a-term at `u` is `≥ X_a + u r_a` and `≥ 0` (`r_a` = rate of `a`'s hi end `H + β_a u` and of its lo
end `max(ζ, L + α_a u) ≤ max(ζ, L) + α_a⁺ u`), the b-term `≥ X_b + u r_b` and `≥ 0` (hi end
`min(ζ + u, H + β_b u) ≥ min(ζ, H) + min(1, β_b) u`).  Using the affine form where the term is positive on the
`ζ`-cell and 0 elsewhere, the pair is `≥ Q(b, u) = inf_ζ [ψ₀(ζ) + u r(ζ)]`, an infimum of affine functions of `u`,
hence **concave in `u`**; `Q(b, 0) = P(b)` is Lemma Z's limit value, computed exactly as a function of the base.

**(d) Rectangles.**  `d/dt area(Q(c(t), θ(t)) ∩ R) = θ'·Σ_edges ∫_{edge ∩ R} (−s) ds + c'·N` (Reynolds; `s` along an
edge from its midpoint; `N` = `∫ n` over `∂Q ∩ R`; `θ' ∈ [2/(1+u1²), 2]`, `c' = σ w'/2 ∈ σ [w'(u1)/2, 1]`), enclosed
over the box extended towards the base (64 pieces per edge; `s > 0` counted where possibly inside, `s < 0` where
certainly inside).  So `area ≥ area_base(b) + u s_A`.

**(e) Edge lines.**  A line whose lo (left-type) or hi (right-type) end is not dominated has, if the centre stays on
one side of it, a chord of length `≥ min(1, τ/S(u1))` (triangle cut, Lemma 2's chord profile) with
`τ ≥ ½ + ℓ − b_⊥` (left) resp. `½ + b_⊥ − ℓ` (right) **for every `σ`** (since `w ≥ 1`), so its mass is
`≥ ρ_min · min(1, τ/S(u1))`, a function of the perpendicular base coordinate.  A germ pair contributes the larger of
its pair value and its two edge terms; the bound tries either choice per pair.

**(f) Assembly.**  For every pose, `μ ≥ (certain points and atoms) + C(b) + u S + Σ_pairs Q_p(b, u)` where `C` collects
windows, corrections, edge terms and base areas, `S` the window and area slopes.  For each `b` this is concave in
`u`, so its minimum over `u ∈ [0, u1]` is at `u = 0` or `u = u1`; for each `u` it is, on every cell of the candidate
grid (breakpoints `± ½`, rectangle sides `± ½`, edge and correction kinks), concave in each base coordinate
(windows linear, pair values minima of linear functions, area bilinear), so the minimum is at a grid vertex.
**`L_U = (points) + min over the grid of min(total(b, 0), total(b, u1))`**, exact integers except the slopes and
rates (outward-rounded) and the pair-rate products (rounded down).  `L(B) = max(L(B), L_U)` for boxes that the other
bounds do not certify (`--first-order`; off by default, log header `atoms=grid+fo`).

### 11.3 What Lemma U certifies

At the zero-limit germs the base value is exactly the germ limit (1) and the slope the true first-order term minus
`O(box)`, so the boxes containing the limit certify: the **wall band** (`σ_x = +1`), the **single column germs**
(`σ_x = +1` with the edge term of `x = 1`: in region `c_x < 1 + w/2` the corner chord of `x = 1` grows like
`τ/S(u1)` and dominates; elsewhere the shift of the horizontal lines and the area give slope `≈ +0.43` per unit `u`
against a true `0.44`), and the wall/pair germs at `c_y ∈ ½ + ℤ/5` on the wall.

### 11.4 The double germs (closed 2026-10-01)

Until 2026-10-01, `cert --d4 --first-order` left 120 boxes (`--full`: 1,108), all at `c_x ∈ [1.5, 1.50006]`,
`u ≤ 7·10⁻⁵`, `c_y` within `5·10⁻⁵` of `1.5, 2.3, 2.5, 2.9, 3.1, 3.5` and one at `(2.3, 2.3)`: the column germ of
`x = 1` **and** a horizontal pair germ at once (resp. the corner of `R`).  The exact germ minimum there is `1 + 0.16 θ`
(e.g. `(1.5, 3.1)`: `1.0000003229` at `θ = 2·10⁻⁶`), and the bound was a few `10⁻⁶` short.  Diagnosed with the
per-term debug output (`ZMX2_DEBUG=1`, `ZMX2_DEBUG_CELL=1`, `ZMX2_DEBUG_PAIR=1` with `zmx2 box`), the losses were
four different ones:

1. **(1.5, 3.1), (1.5, 2.9), (1.5, 2.5), (1.5, 3.5)** — a base-anchored pair end.  With `σ_x = +1` the base range is
   `[1.5 − ext, …]`, so the lo end `b − ½` of the bottom pair line `y = 2.6` sweeps across `x = 1`, where its density
   drops from `0.276` to `0.115`; the old pair rate took the larger density for the whole range (a loss of
   `2·0.161·u` per unit `u`, i.e. `2.4·10⁻⁶` at these boxes), although for bases `b ≥ 1.5` it never meets the `0.276`
   piece (and for `b < 1.5` the edge term of `x = 1` is huge).  Fix: the breakpoint correction (b′) per end, now
   also for the two base-anchored ends of a pair, evaluated per base cell (§11.6).
2. **(1.5, 2.3)** — an uncut pair taken as cut.  The base function contains the pair value `P = inf_ζ`, independent
   of the perpendicular base offset `δ = b_y − 2.3`; but for `δ > u` the pair `y = 1.8 / 2.8` is not cut at all (the
   top line lies inside `Q`), and the vertical windows lose `0.11 δ`: base value `1 − 1.3·10⁻⁶` at `δ = 1.2·10⁻⁵`,
   `u = 0`, at every depth.  Fix: **Lemma V** (§11.5), the cone `|η| ≤ K u` of the pair.
3. **(1.5, 1.5)** — the pair forms.  On a `ζ`-cell ending at `H` the a-term was used in its affine form
   `X_a + u r_a` even at `ζ = H` where `X_a = 0` and the true bound is `max(0, ·) = 0`, so its hi-end loss
   (`0.326 u`) was charged although the b-line has no density there; and the `ζ`-cell rates took the smallest density
   over the cell plus the sweep, which crosses the next breakpoint.  Fix: per `ζ`-cell choice of the form (affine or
   `0`) and `ζ`-candidates `bp ± sweep` (§11.6).
4. **(2.3, 2.3)** — the corner of `R` (`Q` inside `R` touching its left and bottom sides: `μ = 1` along the whole
   corner-touching curve).  Here the plain bound (config 1, Lemma K on both sides) was short by **one grid unit**
   (`1.5·10⁻¹⁵`): the chord range of `y = 1.8` was rounded outward below `x = 1.8`, where its density is `0`.  Fix:
   **Lemma K (d)** (§11.7), an exact clip of the chord ends.

With these, `cert --d4 --first-order` and `cert --full --first-order` leave nothing (§11.8), the maximum depth drops
from 40 to 23 and the runs are 2.5× faster.

### 11.5 Lemma V: the cone of a germ pair

Take a germ pair of a family (frame as in Lemma T): line `a` at `ℓ`, line `b` at `ℓ + 1`, `a` with its hi end and
`b` with its lo end dominated, on the box `B'` (`u ∈ [0, u1]`, `u1 ≤ 1/16`).  Let `η = d_a − ½ = c'_⊥ − g'` with
`g' = ℓ + ½` (frame perp coordinate), and assume `|η| < 1` on `B'`.  Let `α_g` be the g-end slope of `a` and `β_g`
that of `b` (§11.2 (a): `g1 ≤ −½ + α_g u`, `g2 ≥ ½ + β_g u` on `B'`), and
`K_a = 1 + 2 max(0, −α_g) u1`, `K_b = 1 + 2 max(0, β_g) u1` (rounded up).

**Lemma V.**  Every pose of `B'` with `u > 0` satisfies (L) `η ≤ −K_a u`, (R) `η ≥ K_b u`, or (C) `c'_⊥ = g' + κ u`
with `κ ∈ [−K_a, K_b]`.  In (L), `f1_a ≤ −½ + min(0, α_g u1)`: the lo end of `a` is dominated, `a` is an ordinary
window line (`lo ≤ b − ½ + α u`, `α` as in (a)).  In (R), `f2_b ≥ ½ + max(0, β_g u1)`: `b` is a window line.  At
`u = 0` (axis-parallel squares) `η ≤ 0` gives (L) and `η ≥ 0` gives (R) in the same sense.
*Proof.*  By Lemma C, `f1_a = (d_a − ½)/(2u) − (d_a + ½) u/2 = η/(2u) − (1 + η) u/2 ≤ η/(2u)` (as `1 + η > 0`); in (L)
this is `≤ −K_a/2 = −½ − max(0, −α_g) u1 = −½ + min(0, α_g u1)`.  For `b`, `d_b = d_a − 1 = η − ½` and
`f2_b = (d_b + ½)/(2u) + (½ − d_b) u/2 = η/(2u) + (1 − η) u/2 ≥ η/(2u)` (as `1 − η > 0`), `≥ K_b/2` in (R).  The
domination then gives the end bounds exactly as in (a): `lo = c + max(f1, g1) ≤ c − ½ + α_g u` because
`−½ + min(0, α_g u1) ≤ −½ + α_g u` for `u ≤ u1`, and symmetrically `hi_b ≥ c + ½ + β_g u`.  The three cases cover
all `η`.  At `θ = 0` the closed axis-parallel square contains `a` iff `|d_a| ≤ ½`, i.e. (as `d_a > −½`) iff `η ≤ 0`,
and then its chord is `[c − ½, c + ½]`, the window at the base `c`; likewise for `b`. ∎

**Use.**  On an axis whose path sign is free (no wall within reach), if the germ `g` of a pair cutting that axis
(vertical pairs cut `x`: `g = ℓ + ½`; horizontal pairs cut `y`: `g_y = y_a − ½`) lies in the box's centre range and
`|c − g| < 1` on the box, the poses are split into three **pieces** (`germ_splits`, `first_order_bound`):

* **(C) cone:** `c = g + κ u` (real coordinate; vertical pair `κ ∈ [−K_a, K_b]`, horizontal pair `c_y = g_y − η`,
  `κ ∈ [−K_b, K_a]`).  Lemma U with the **single base `g`** on this axis and the path slope interval `[k0, k1]` in
  place of `σ`: in (a) `α += k1`, `β += k0` (the along coordinate of the other family is `g + κ u`); in (d)
  `c'_⊥ = κ ∈ [k0, k1]`, the Reynolds path `c_⊥(t) = g + κ t` stays in the box; in (e) the depth of a vertex beyond a
  line with this perp axis is `τ = ℓ − c' + w/2 = ℓ − g' + ½ + ((w − 1)/2 − κ' u) ≥ ½ + ℓ − g' − (q1 − wlo)⁺ u1`
  (left-type; right-type `≥ ½ + g' − ℓ − (−q0 − wlo)⁺ u1`), with `[q0, q1]` the frame slope interval (`= [k0, k1]`
  for the vertical family, `[−k1, −k0]` for the horizontal one) and `(w − 1)/2 ≥ wlo u`, `wlo = (1 − u1)/(1 + u1²)`:
  the edge terms are shifted by this slack.  The pair itself keeps its value `inf_ζ`.
* **(L), (R):** Lemma U on the half of the centre range on the side of `a` (resp. `b`) of `g`, with that line a
  window line (`lo_ok` resp. `hi_ok` forced) and the other one an unpaired edge line; path signs as usual.

Every piece's formula is a lower bound of `μ` at the poses of that piece, and the corner-minimum argument of (f)
(§11.6) uses only the algebraic form of the formula on the base cells covering the piece's bases (the single base
`g`; the half range), so it is irrelevant that, e.g., the (L) formula is not a lower bound outside (L).  The bound of
the box is the larger of the plain Lemma U bound and, over the non-empty subsets of the splittable axes (at most one
pair per axis), the minimum over the pieces of the best over the path signs of each piece; each subset is a
partition of the poses.  It is tried only when the plain bound fails.

### 11.6 Lemma U per base cell (replaces (b′), (c), (f) of §11.2)

The base grid of each axis (kinks of every term, rectangle sides `± ½`, range ends, and now also the correction
kinks `bp ± ½ ∓ sweep` and the pair candidates `bp ± sweep ± ½`) cuts the base range into **cells**; all rates,
corrections and pair forms are computed per cell, uniformly for the whole cell.  `sweep(k) = ⌈|k| u1 G⌉ + 1 ≥ |k| u1`
(grid units).

**(b) One moving end** (`end_rate`).  Let an end be `e = b + off` with `b` in the cell `[c0, c1]`, moving to `e + k u`,
`u ∈ [0, u1]`.  For a hi end, `F(e + ku) − F(e) ≥ k ρ u − corr(b)`; for a lo end, `F(e + ku) − F(e) ≤ k ρ u + corr(b)`,
where `ρ` is the smallest (gain: hi end with `k ≥ 0`, lo end with `k ≤ 0`) or largest (loss) density on the swept
range (`(c0 + off, c1 + off + sweep)` for `k ≥ 0`, `(c0 + off − sweep, c1 + off)` for `k < 0`) and `corr = 0`; or, for
a loss when the swept range contains exactly one breakpoint `p` (densities `ρ_l`, `ρ_r` on its two sides), `ρ` = the
density beyond `p` and `corr(b) = (ρ_before − ρ_beyond)⁺ · min(dist(e, p)⁺, sweep)` with `dist = e − p` (moving left)
resp. `p − e` (moving right).
*Proof.*  `F` is the integral of the density.  Gains and plain losses: the integral over the swept interval, which
lies in the range.  Correction (moving left from `e` by `δ = |k| u ≤ sweep`): if `e ≤ p` the interval lies left of
`p`, `∫ = ρ_l δ`; otherwise `∫ = ρ_r min(e − p, δ) + ρ_l (δ − min(e − p, δ)) ≤ ρ_l δ + (ρ_r − ρ_l)⁺ min((e − p)⁺, δ)`;
moving right symmetrically. ∎  `corr` is `u`-independent, with kinks at `p − off` and `p − off ∓ sweep`, which are
grid points.  A window is `F(hi) − F(lo) ≥ W(b) + [hi-end term] − [lo-end term]` (also when the chord is empty, as
`F(hi) − F(lo) ≤ 0` then); the two base-anchored ends of a pair (`a`'s hi end `b + ½`, `b`'s lo end `b − ½`) use the
same lemma.  Per cell either both options (correction or not) are kept as *variants* (below).

**(c) A pair on a base cell** (`u_pair_cell`).  Let `C = [c0, c1]`, `H = b + ½`, `L = b − ½`.  The `ζ`-axis is cut at
the **candidates**: `L`, `H` (moving with `b`) and the constants `bp`, `bp ± sweep_p` for the breakpoints of both lines
near `C` (`sweep_p` = the larger sweep of `a`'s lo end (rate `α_a⁺`) and `b`'s hi end (rate `min(1, β_b)`)).  A moving
candidate meets a constant one only at `b = const ∓ ½`, a grid point, so **the order of the candidates is the same at
every `b` inside `C`** (and its closure at the corners).  Between consecutive candidates (plus the two half-infinite
sentinel cells) lie the `ζ`-cells.  For every `ζ` (§11.2 (c), Lemma Z: `H1_b = ζ + u`):
`a-term ≥ max(0, X_a + u r_a) − corr_a`, `b-term ≥ max(0, X_b + u r_b) − corr_b`, `X_a = F_a(H) − F_a(max(ζ, L))`,
`X_b = F_b(min(ζ, H)) − F_b(L)` (unclamped), with `r_a` = `a`'s hi-end rate (by (b)) `− α_a⁺ ρ_max` over the range of
`max(ζ, L) + [0, sweep]` for `ζ` in the cell and `b ∈ C`, and `r_b` = `b`'s lo-end rate `+ min(1, β_b) ρ` over the
range of `min(ζ, H) + [0, sweep]` (`min(ζ+u, H + β_b u) ≥ min(ζ, H) + min(1, β_b) u`).  On each `ζ`-cell a **form** is
fixed for the whole base cell and all `u`: for each term either `0` (always valid) or the affine form (allowed only if
the term is *active*: `X > 0` at one of the cell's ends at one of the corners — then the cell lies left of `H`
(resp. right of `L`), because `H` (`L`) is a candidate, so `X ≥ 0` is its unclamped value on the whole cell).
*Claim:* the pair mass at every pose with base `b ∈ C` is `≥ Π(b, u) − corr_a(b) − corr_b(b)`, where `Π(b, u)` is the
minimum over the `ζ`-cells and their finite ends `z(b)` of `V(b, u) = form(X_a(z(b), b)) + form(X_b(z(b), b)) +
u (rates of the affine forms)`.
*Proof.*  For fixed `(b, u)` the bound is, on each `ζ`-cell, linear in `ζ` (the kinks of `X_a`, `X_b` in `ζ` are
breakpoints, `L`, `H`: candidates), so its infimum is at a finite end (on a sentinel cell it is constant in `ζ`:
there `max(ζ, L) = L` resp. `min(ζ, H) = H` and the other term is inactive).  Along a candidate `z(b)` (constant or
`b ± ½`) `X_a(z(b), b)` is linear in `b` on `C` (`F_a(b + ½)` and `F_a(max(z, b − ½))` cross no breakpoint, and the
`max` switches only at grid points); the rates are constants.  So `V` is jointly affine in `(b, u)` on `C × [0, u1]`
and `Π` is jointly concave. ∎  The forms are chosen per `ζ`-cell: among those that keep the `u = 0` values at the two
corners `≥ P` (the exact limit, the minimum of `X_a + X_b` over all cells and ends; the all-active forms do), the one
with the best `u = u1` values.  So `Π(b, 0) = P(b)` exactly and only the `u`-dependence is optimised.

**(f) Assembly per 2-d cell.**  For every 2-d base cell (`x`-cell × `y`-cell) the bound is the best, over the
correction variants of the two 1-d cells (all off / all on / per line the better one by its own worst value) and the
pair choices (pair value, or the edge terms of its two lines on the perpendicular axis; all `2^n` for `n ≤ 4` pairs),
of the minimum over the 4 corners and `u ∈ {0, u1}` of `konst + Σ terms + Σ Π − corr + area(corner) + u S`.  The box
bound is the minimum over the 2-d cells.  *Proof.*  Every pose's base lies in some closed 2-d cell, and for every
fixed choice the formula is a lower bound of `μ` at the pose ((a)–(e), (b), (c)).  For fixed `u` it is concave in
`b_x` for fixed `b_y` on the `x`-cell (windows and corrections linear, edge terms concave (`ρ min(G, t/s1)`, kinks on
the grid), `Π` concave, `ρ ox(b_x) oy(b_y)` linear in `b_x`), likewise in `b_y`; for fixed `b` it is concave in `u`
(affine except `Π`).  Hence on the cell × `[0, u1]` it is `≥` its minimum over the corners and `u ∈ {0, u1}`.  All
rates are bounds for `u ∈ [0, u1_hi]` with `u1_hi ≥ u1` the outward enclosure of `u1`, and a concave function on
`[0, u1_hi]` is `≥ min(f(0), f(u1_hi))` on `[0, u1]`. ∎  Rounding: rates times `u1` are outward-rounded products
(`rate_u1`, rounded down), every `X` and window is exact, the pair values at `u1` are floored after a relative
`10⁻¹⁵` margin (values `< 2⁸⁰`, as before).

*The old version had a gap here* (§13, finding U1): the old (c′) recomputed the `ζ`-cell rates and activity flags at
each base vertex from the `ζ`-cell ends *at that vertex*, with the whole box as base range, so its vertex values were
not values of one concave function and "the minimum is at a vertex" was not proved.

### 11.7 Lemma K (d): exact clip of the chord ends

In Lemma K (and only there), the ranges of the chord ends are additionally clipped exactly: **if `d ≤ −u1/(2(1 − u1²))`
on the box, then `lo ≥ c − ½ ≥ y0 − ½`; if `d ≥ u1/(2(1 − u1²))`, then `hi ≤ y1 + ½`** (frame, `y0`/`y1` the box's
along range; `y0 − ½` is rounded down to the grid, `y1 + ½` up).
*Proof.*  `g1 = −dT − R/2 = −dT − ½ − u²/(1 − u²)` and `−dT ≥ 2|d| u` for `d ≤ 0`; `u²/(1 − u²) ≤ 2|d| u` iff
`|d| ≥ u/(2(1 − u²))`, which holds for all `u ≤ u1` (the right side increases with `u`).  So `g1 ≥ −½` and
`lo = c + max(f1, g1) ≥ c − ½`; at `u = 0` the limit chord starts at `c − ½`.  The other side is symmetric
(`g2 = −dT + ½ + u²/(1 − u²) ≤ ½`). ∎  At the corner of `R` this removes the one grid unit of outward rounding that put
the chord's lo end into the zero-density part of `y = 1.8` (§11.4, item 4).

### 11.8 What it certifies; runs

All of `μ₇`'s first-order germs, single and double, and the corner of `R`: `cert --d4 --first-order` **`VERIFIED-D4`**,
`cert --full --first-order` **`VERIFIED`** (§9).  The cost is lower than before (fewer, shallower boxes: max depth 23
instead of 40).  The rejection tests A9 (§10) show that the new machinery does not certify covers that are invalid
only at `0 < θ < 0.003°` near a double germ.

## 12. Caveats

* The new geometry (cap bounds, corner clipping, chord ranges in Lemma K) uses outward-rounded binary64 intervals
  like the rest of `zmx2`; `cert0` is pure integer arithmetic.
* The harness and the code were written by the same agent (an independent *implementation* of `μ`, not an
  independent author), as in `ZMX2.md`.
* Lemma K couples a line with only one cap; a line that is the side of two rectangles is coupled once (the other
  cap stays in the union bound).  Only axis-parallel rectangles are supported.
* **Lemma U (§11.2, §11.6) and Lemma V (§11.5) are the least mature part** and carry the only certificate of
  `0 < θ < 0.014°`: many moving pieces (path bases, cones, per-cell concavity, pair forms, f64 rate products, the
  64-piece Reynolds enclosure); written by two sessions of the same kind of agent; §13 is a self-audit, not an
  independent review.  The mutation tests of §13 show that several of its refinements (cone width, pair activity,
  sweep candidates, the `u = 0` form constraint) are *not* exercised by any test on `μ₇` — their correctness rests on
  the written proofs alone.  Nothing of §§2–8 depends on them (except Lemma K (d), §11.7), and they are off unless
  `--first-order` is given.
* Refused runs stop a root after 200 uncertified boxes (`--uncert-cap`), so their counts understate the uncertified
  volume; the cells involved are listed in §11.4.

## 13. Self-audit of Lemma U and Lemma V (2026-10-01)

Adversarial re-derivation of every step, by the agent that wrote §11.4–11.8 (so **not** an independent review).
For each step: the claim re-proved, how it was attacked, the finding.

| # | step | check | finding |
|---|---|---|---|
| U1 | old (c′): pair value as a function of the base | re-derived the concavity argument of (f) | **Proof gap (now closed).**  The 2026-09-30 code computed the `ζ`-cell rates (`max(ζ, L)`, `min(ζ, H)` ranges) and the activity flags of the pair terms *at each base vertex* from the `ζ`-cell ends at that vertex, with the whole box as base range.  The vertex values are then pointwise lower bounds, but not the values of one concave function on a cell (e.g. a rate with a smaller density range at the right vertex than in the interior), so "the minimum over the cell is at a vertex" was not proved.  The new `u_pair_cell` fixes rates, forms and the candidate order per base cell (§11.6 (c)).  No violation was ever observed (the old harness had 0 FAIL, and the old certificate claimed only `θ ≥ θ_min` plus the boxes Lemma U certified), so this is a gap in the argument, not a known false certificate; everything the old code certified is re-certified by the new code. |
| U2 | (a) end bounds with the path shift | `(w − 1)/2 = u(1 − u)/(1 + u²) ∈ [wlo u, u]` (the factor `(1 − u)/(1 + u²)` decreases on `[0, 1]`); `g1 ≤ −½ + α_g u` from `T ∈ [2u, 2u/(1 − u1²)]`, `R ≥ 1`; Lemma W: `½ − q = u²(1 − u)/(1 + u²) ≤ u1 u` | correct |
| U3 | (b) one moving end, corrections | §11.6 (b) proof; corrections need *exactly one* breakpoint strictly inside the open swept range (`one_bp`), densities `0` off the support; kinks `p − off`, `p − off ∓ sweep` are on the grid (also for the pair ends) | correct; mutation `corr` (correction term dropped, density beyond kept) is caught by the fine harness (`1.5,2.3`, reflected: `μ − bound = −3.1·10⁻⁷`) and by exact `μ` on 8 of 470 boxes where it raised the bound |
| U4 | (c) Lemma Z identity, `ζ`-cells, forms | `f2_b − f1_a = u` for lines at distance 1; `min(ζ + u, H + β u) ≥ min(ζ, H) + min(1, β) u`; affine form only on cells left of `H` (resp. right of `L`), where the unclamped `X ≥ 0` | correct; mutations `act` (affine form everywhere), `sweep` (no sweep in the `ζ`-cell rates) and `ustar` (forms chosen without the `u = 0` constraint, `u = 0` value kept at `P`) are unsound in principle but **not caught**: `act` changes no bound on 9,205 germ boxes, `sweep`/`ustar` raise 56/84 bounds but stay below the sampled exact `μ` (gaps `≥ 1.2·10⁻⁷`) |
| U5 | (d) Reynolds | `J(p − c)·n = −s` (with `Jn = τ`, `Jτ = −n`), `θ' = 2/(1 + u²)`, `c' = σ w'/2`, `w'/2 = (1 − 2u − u²)/(1 + u²)²` decreasing; pieces end at multiples of `1/64` so none straddles `s = 0`; "maybe" pieces only add negative `−s` and `hull(n ds, 0)` | correct |
| U6 | (e) edge terms | chord `= min(1/P, (M − d)/(pP)) ≥ min(1, τ/s1)` (`P ≤ 1`, `p ≤ s1`); `τ ≥ ½ + ℓ − b` for every `σ ∈ {−1, 0, 1}`; cone slack §11.5 | correct |
| U7 | (f) assembly | per 2-d cell, every fixed choice valid, concave in each coordinate and in `u` (§11.6 (f)); rates for `u ≤ u1_hi` and evaluation at `u1_hi ≥ u1` is safe; `inw` = points certain on the original box, valid because only its poses are claimed | correct |
| U8 | `θ = 0` poses inside Lemma U boxes | windows (lines with dominated ends have `|d| ≤ ½`, closed square: full chord), pairs (`inf_ζ` ≤ either full window), edges, area at `u = 0` | correct (and `θ = 0` is certified independently by `cert0`) |
| U9 | rounding | `X`, windows, areas exact integers; rates via outward intervals (`rate_u1`); pair values at `u1` from `i128 → f64` (values `< 2⁸⁰`, relative error `2⁻⁵³`) with a `10⁻¹⁵` relative margin before the floor | correct |
| V1 | Lemma V | re-derived (§11.5): needs `1 ± η > 0` (checked: `|c − g| < 1` on the box), `K` rounded up and used identically for the regions and the cone interval; frame signs for horizontal pairs (`η = g_y − c_y`, so (L) is the *upper* half and `κ ∈ [−K_b, K_a]`) checked against the code | correct; mutation `K` (`K_a = K_b = ½`, cone too narrow) changes **no** bound on 9,205 germ boxes: the cone width is not load-bearing for `μ₇` and is not exercised by any test |
| K1 | Lemma K (d) | §11.7 proof | correct; mutation `clip` (clip at `y0 + ½`) is caught: rejection cover 2 becomes `REGION CLEAN`, and its fine harness finds 20 FAIL (`μ − bound = −5.2·10⁻⁷`) |

**Targeted exact checks.**  Fine harness near all 8 double germs, plain and reflected cover (16 × 300 boxes,
110,400 exact poses at angles down to `u1·2⁻⁴⁰`): 0 FAIL, `min(μ) − bound = 0` at many boxes (the harness reaches the
exactly tight poses).  Same on both rejection covers: 0 FAIL.  The rejection covers show that the new machinery
refuses covers that are invalid **only** for `0 < θ < θ_min` (rejection 1: `cert0` and `--umin 10` accept it).

**Mutation runs** (`mutate.sh`-style: the mutant built in a scratch tree; checks: both rejection covers, the `μ₇`
double-germ column, 12 fine-harness runs): unmodified: all pass; `corr`, `clip`: caught; `K`, `act`, `sweep`,
`ustar`: survive (see U3, U4, V1).  **Residual risk:** these four refinements are justified only by the written
proofs in §11.5–11.6; an error in them would not be detected by the present tests.  A second reader should check
them first (and, ideally, write an independent harness: the present one shares its author with the code).
