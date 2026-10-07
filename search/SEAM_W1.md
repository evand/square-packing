# Seam capacity at w = 1 (exact), the dual as unfilled area, and what changes at w = 2 (2026-10-02)

Setting: the band LP of QUADRANT.md §4 / FRIEDMAN.md §8.  `μ ≥ 0` on `H = {y ≥ 0}`, 1-periodic in x, Lebesgue on
`y > w`, band mass `μ([0,1) × [0,w]) = w` (σ = 0), valid (`μ(S) ≥ 1` for every closed unit square `S ⊂ H`).
`M(w) := sup μ({0} × [0,w])`, the **seam capacity**.  Labels **[proved] [measured] [heuristic]** as elsewhere.
Code: `seam_w1_check.py` (§2 inequality), `seam_cand.py` (oracle on a fixed profile), `seam_dual_load.py` (§4),
`seam_xu.py` (§5).  Context and the conjecture ladder this feeds: FRIEDMAN.md §9.

## 0. Summary

* **`M(1) = 3/4`** [upper bound proved; matching profile proved modulo a written-out 2-variable casework, §2.3].
  Upper bound: one 45° square (§1).  Optimal profile: seam segments `{n} × [½, 1]` at density 3/2 plus the line
  `y = 1` at density ¼ (§2).  Every optimum puts all its non-seam mass in one triangle (§1.2), and the band LP
  optimum `runs/quad_band_w1` matches this to 5 digits.
* **Seam-avoiding dual certificates give exactly `M(w) ≤ w − ¼`, never better** [proved, §3].  Tight at w = 1;
  useless at w ≥ 2.  A certificate for w ≥ 2 must use squares that cross the seam.
* **Dual = unfilled area** [proved, §4]: a periodic fractional packing with load `ℓ ≤ Λ` on the band and `ℓ ≤ Λ − 1`
  on the seam certifies `M(w) ≤ ∫_band (Λ − ℓ)`.  Checking one is a finite arrangement computation.
* **Correction**: the saved w = 2 dual (`runs/seam_dual/w2.npz`, 44 poses) is a **lattice artifact**: its load reaches
  182 against Λ = 92 on slivers crossing the seam between atom heights (§4.1).  The band values `M(w)` = 1.433,
  1.965, 2.116, 2.231 (w = 2, 4, 5, 6) are therefore **lower bounds** on the continuum `M(w)` (restricted-column
  LP, oracle-validated), not caps.  FRIEDMAN.md §2.2 / §7 C5 used them as caps; corrected there.
* **The x-uniform + seam class** (`μ = λ(dy) ⊗ Leb_x + ρ(dy) ⊗ δ_ℤ(x)`) is optimal at w = 1 (0.7500) and loses
  ≈ 9 % at w = 2 (≈ 1.300 vs 1.433) [measured, §5].  So at w ≥ 2 the optimum needs structure across the period.
* **Corner coupling** at w = 1 costs exactly ¼ [measured: D = 0.5004 at R = 2, QUADRANT §4]; unexplained.  Across
  w = 1..5 the ratio `D/M ≈ 0.6–0.67` [measured, heuristic, FRIEDMAN §9.3].

## 1. Upper bound `M(1) ≤ 3/4`  [proved]

Let `S_δ` be the closed unit square at 45°, centre `(½, ½ + √2/2 + δ)`, `0 < δ < ½`.  Its part below `y = 1` is the
triangle `T_δ = {½ + δ ≤ y ≤ 1, |x − ½| ≤ y − ½ − δ}`, of height `h = ½ − δ` and area `h²`; at `y = 1` it spans
`[δ, 1 − δ]`, so `T_δ` lies in `(0,1) × [0,1]` and misses every seam line `x ∈ ℤ`.  Validity:
`1 ≤ μ(S_δ) = μ(T_δ) + Leb(S_δ ∩ {y > 1}) = μ(T_δ) + 1 − h²`, so `μ(T_δ) ≥ h²`.  The period cell `[0,1) × [0,1]`
has mass 1 and contains the disjoint sets `{0} × [0,1]` (mass `m_v`) and `T_δ`, so `m_v ≤ 1 − h²`.  Let `δ → 0`.  ∎

The same square placed at the top of a width-w band gives `M(w) ≤ w − ¼` for every w.

### 1.1 Reading
The oracle's last worst pose in the band LP run was `(0.499, 1.2072, 44.76°)`, i.e. `S_δ`.  The LP dual at w = 1 is
this one square (plus the σ multiplier).

### 1.2 Structure of every optimum  [proved]
If `m_v = 3/4` then `μ(T_δ) ≤ ¼` and `≥ (½ − δ)²` for all δ, and the band mass outside `seam ∪ T_δ` is
`≤ 1 − ¾ − (½ − δ)² → 0`.  So **every optimal w = 1 profile has mass ¾ on the seam, ¼ on
`T_0 = {½ < y ≤ 1, |x − ½| ≤ y − ½}`, and 0 elsewhere.**  Check [measured] on `runs/quad_band_w1/support.txt`:
seam 0.75003, `T_0` 0.24997, outside 0.000000.

## 2. A matching profile: `M(1) ≥ 3/4`

**Profile A**: Lebesgue on `y > 1`; the line `y = 1` at density ¼ (all phases); seam segments `{n} × [½, 1]` at line
density 3/2.  Band mass per period `¼ + ¾ = 1` (σ = 0), `m_v = ¾`.

### 2.1 Reduction to two variables  [proved]
For a closed unit square `S` with lowest point at height `y0 < 1` let `W(y)` be its horizontal slice width,
`A = ∫_{y0}^{1} W` the band area, `N(y)` the number of integers in the (closed) slice.  Then
`μ(S) = (1 − A) + ¼ W(1) + (3/2) ∫_{1/2}^{1} N(y) dy`.  A closed interval of length `≥ 1` contains an integer, so
`N(y) ≥ 1[W(y) ≥ 1]`, and validity follows from

> `F(θ, y0) := ¼ W(1) + (3/2) |{y ∈ [½, 1] : W(y) ≥ 1}| − A ≥ 0`,

which no longer depends on the horizontal position.  By the mirror symmetry `θ ∈ [0°, 45°]` suffices.  With
`s = sin θ ≤ c = cos θ`, `H = s + c`, depth `d = y − y0`: `W = d/(sc)` on `[0, s]`, `1/c` on `[s, c]`,
`(H − d)/(sc)` on `[c, H]`, and `W ≥ 1` exactly for `d ∈ [sc, H − sc]`.

### 2.2 Cases done by hand  [proved]
* `y0 ≥ 1`: `A = 0`.
* `θ = 0` (`W ≡ 1`): `F = ¼ + (3/2)(1 − max(½, y0)) − (1 − y0)`, which equals `y0` for `y0 ≤ ½` and
  `¼ + ½(1 − y0)` for `y0 ≥ ½`; `≥ 0`, with equality only for the wall-resting squares `[t, t+1] × [0, 1]`.
* `θ = 45°`, `y0 = ½` (the §1 square): `A = ¼`, `W(1) = 1`, no `W ≥ 1` slice in the band: `F = 0`.
* `θ = 45°`, `y0 = 0` (diamond standing on the wall): `A = 1 − (√2 − 1)²`, `W(1) = 2(√2 − 1)`, `W ≥ 1` on
  depths `[sc, H − sc] = [½, √2 − ½]`, which meets `[½, 1]` in length `√2 − 1`: `F = ½(√2−1) + (3/2)(√2−1) − 2(√2−1) = 0`.

### 2.3 General `(θ, y0)`  [measured; exact casework TODO]
`seam_w1_check.py` evaluates `F` in closed form on a 4501 × 10001 grid of `(θ, y0)`: **min F = −3·10⁻¹⁶** (rounding),
attained only at the tight families above (and trivially at `y0 = 1`).  The oracle on profile A (`seam_cand.py`:
lattice, 3·10⁵ random, near-lattice, germ poses, polish) finds min 1.000000.  The exact proof is a piecewise
calculus exercise with breakpoints `h = 1 − y0 ∈ {sc, s, c, H − sc, ½ + sc}` (the last two from the ends
of the `W ≥ 1` interval meeting `y = 1` and `y = ½`); or exact interval arithmetic over two
variables.  **Tight set: three families** — wall-resting axis squares, the raised 45° square of §1, the 45° square
standing on the wall.  (The x-uniform LP of §5 finds a different optimum, seam concentrated near `y ≈ 0.8–0.85`; the
optimal set is not a point.  Raised axis squares force seam mass to sit at height `≥ ¾` up to the line's ¼.)

## 3. Seam-avoiding certificates cap at `w − ¼`  [proved]

A "Λ = 1" dual certificate uses only squares whose band part avoids the seams (load 0 on the seam, `≤ 1` elsewhere);
it certifies `m_v ≤ w − Σ y_S area(S ∩ band)`.  **Claim: `Σ y_S area(S ∩ band) ≤ ¼`, so the best such bound is
`w − ¼`.**  Proof: the seams cut the band into open strips of width 1, and a unit square with a slice of width
`≥ 1` inside the band meets a seam, so a seam-avoiding square crosses `y = w` and its band part is the region of depth
`≤ sc` below its lowest vertex: a triangle with top chord `L ≤ 1` on `y = w` and area `L² sc/2 ≤ L²/4`.  Load `≤ 1`
just below `y = w` gives `Σ y_S L_S ≤ 1` per period, so `Σ y_S L_S²/4 ≤ ¼ Σ y_S L_S ≤ ¼`.  ∎
(The union of all such triangles on one chord is a half-disk — apexes on the Thales circle — of area π/8 > ¼, which
is why the chord-length argument, not an area argument, is needed.)

## 4. The dual as unfilled area  [proved]

Band LP (QUADRANT §2, band mode): minimise `−m_v` s.t. σ row (multiplier λ) and pose rows `A x ≥ 1 − Leb_S`.
For dual weights `y_S ≥ 0` define the periodic load `ℓ(p) = Σ_S y_S #{j : p + j e₁ ∈ S}`, `Λ := −λ`.  If
`ℓ ≤ Λ` on the band and `ℓ ≤ Λ − 1` on the seam `{x ∈ ℤ}`, then for every valid μ (σ = 0):

> `m_v ≤ ∫_{[0,1) × [0,w]} (Λ − ℓ) dA`.

Proof: `Σ_S y_S (1 − Leb(S ∩ {y > w})) ≤ Σ_S y_S μ(S ∩ band) = ∫ ℓ dμ_band ≤ Λ w − m_v`, and
`1 − Leb(S ∩ {y > w}) = area(S ∩ band)`.  ∎  (Λ is a free scale; §1 is Λ = 1 with one square.)
`ℓ` is upper semicontinuous (closed squares), so near the seam `ℓ ≤ Λ − 1` on every arrangement cell whose closure
meets it: the bound pays for a neighbourhood of the seam, never a sliver unless square edges run along it.
**A certificate is a finite set of squares; validity is a finite arrangement check** — the easy side to make exact.

### 4.1 The saved w = 2 dual is a lattice artifact  [measured]
`seam_dual_load.py 2`: Λ = 92.0, 44 poses, two near-axis squares of weight ≈ 91 and 87.  The formula reproduces
1.4328 (LP 1.4331), but the load reaches **182** in slivers where a square tilted 0.2° overlaps its own periodic copy,
and on the seam line itself.  The LP only has columns on the pitch-0.1 atom lattice, so its dual only constrains the
load there.  Consequence: `M(w)` values from `quad_band_w*` are lower bounds (on the continuum, given the oracle),
and no rigorous upper bound on `M(w)` exists yet for w ≥ 2 beyond `w − ¼`.  Next: a dual-side LP (square weights
as variables, pointwise load rows generated from the arrangement), whose answers can be checked exactly.

## 5. The x-uniform + seam class  [measured]

`seam_xu.py w 0.05`: columns = full-width horizontal lines and slabs at pitch 0.05 in y, seam segments and seam
points; same oracle.  w = 1: **0.75003** (optimal; top line ¼, seam mass mostly on `[0.80, 0.85]`).
w = 2: **≈ 1.3006** after 11 rounds (oracle min 0.9938, still adding rows; stopped for CPU), vs 1.433 unrestricted.
w = 3: stopped after 2 rounds (1.50, falling).  So the phase-free reduction of §2.1 cannot by itself carry w ≥ 2.
The unrestricted w = 2 optimum (QUADRANT §6) uses a near-seam zone (cells at phases 0–½) and a far zone (horizontal
pieces avoiding ≈ ±0.1 around the seam): next ansatz = two zones per row.

## 6. Open / next
1. Exact casework for §2.3 (or interval arithmetic) → `M(1) = 3/4` fully proved.  Lean candidate later (FRIEDMAN §9.6).
2. Corner at w = 1: why exactly ¼?  Read the corner-coupled LP dual (`R2 w1`) the way §1 reads the band.
3. w = 2 from both sides: dual-side LP (rigorous upper bound), two-zone primal (better lower bound).  ~10³ CPU-s each
   on the loaded machine; ask first.
4. Fine-pitch band primal at w = 4–6: is the apparent saturation of `M(w)` the pitch-0.1 grid (SEAM_1D: efficient
   seam structure at small amplitude needs pitch ≤ 0.005)?
