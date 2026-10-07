# Exact certificates for the 1D seam price (FRIEDMAN.md §8, item (a))

Code: `search/seam_1d_exact.py` (checker, certificate writer, tests); LP: `search/seam_1d.py` (now with `--save`).
Certificates: `search/seam_data/seam1d_N{N}_e{e}.json` (+ `.sha256`, + `.check.txt` = checker output), float LP
solutions they were rounded from: `search/seam_data/lp_N{N}_e{e}.npy`.  Labels **[proved] [measured] [heuristic]**.

## 1. Statement

ρ is a measure on the circle [0,1): atoms `w_k ≥ 0` at `x_k = k/N` (`k = 0` is the seam, `g := w_0`) plus a uniform
density `η ≥ 0`.  Total `Σ w_k + η = 1 + e`.  A closed unit square at angle θ with centre abscissa `cx` gets

    F(θ, cx) = η + Σ_k w_k G_θ(x_k − cx),   G_θ(u) = Σ_{j∈ℤ} V_θ(u + j),

`V_θ(d)` = length of the vertical chord at abscissa d of the closed unit square centred at 0, rotated by θ.
(The uniform part contributes `η ∫V_θ = η · area = η`.)  **Claim certified:** `F(θ, cx) ≥ 1` for all θ and all cx.
Then the y-invariant measure `ρ × Leb_y` is valid in the bulk, with seam density g and excess e: κ = e/g is an
upper bound for the bulk price of seam at that g.

## 2. Lemmas

Notation: `c = cos θ, s = sin θ, a = (c+s)/2, b = (c−s)/2`, `r = dist(u, ℤ) ∈ [0, ½]`.

**Lemma 1 (chord formula, symmetries)** [proved].  For `s ≠ 0`:
`V_θ(d) = max(0, min(1/|s|, 1/|c|, (a' − |d|)/|sc|))` with `a' = (|c|+|s|)/2`.  For `θ ∈ (0°, 45°]` this is
`max(0, min(1/c, (a − |d|)/(sc)))`, the formula of the brief.  For `θ = 0`: `V_0(d) = [|d| ≤ ½]` (closed square).
*Proof.*  The square is `{p : |p·n₁| ≤ ½, |p·n₂| ≤ ½}` with `n₁ = (c, s)`, `n₂ = (−s, c)`.  On the line `x = d` the two
slabs cut the y-intervals with centres `−dc/s`, `ds/c` and half-lengths `1/(2|s|)`, `1/(2|c|)`.  Two intervals with
half-lengths `h₁, h₂` and centre distance δ intersect in length `max(0, min(2h₁, 2h₂, h₁ + h₂ − δ))`; here
`δ = |d|(|c/s| + |s/c|) = |d|/|sc|` and `h₁ + h₂ = (|c|+|s|)/(2|sc|)`. ∎
The formula depends on θ only through `{|c|, |s|}`, i.e. through the square as a set (θ and θ + 90° give the same
square; x ↦ −x maps θ to −θ and V is even in d).  So **V_θ = V_{−θ} = V_{θ+90°} = V_{90°−θ} and θ ∈ [0°, 45°] suffices.**
The brief's claim is right.  Caveat (found while writing the checker): for θ slightly *above* 45° the brief's
formula `min(1/c, …)` is wrong near `d = 0` (it must be `min(1/s, …)`).  The checker's θ-range `t ≤ 107/256`
(θ ≤ 45.37°) overshoots 45°, but the bounds are only claimed, and only proved, for θ ≤ 45°; see Lemma 3.
Test [1] compares the formula (after reduction to [0, 45°]) with a polygon-clipping chord at 20000 random
`θ ∈ (−180°, 180°)`: max difference 1.1e-14.

**Lemma 2 (periodised chord in closed form)** [proved].  For `0 < θ ≤ 45°`,

    G_θ(u) = max(P, min(1/c, (a − r)/(sc))),   P := (c + s − 1)/(sc) = 1/c − s/(c(1+c)),

and G_θ is nonincreasing in `r ∈ [0, ½]`.  *Proof.*  `½ < a ≤ 1/√2 < 1`, so among the images `u + j` (distances `r,
1−r, 1+r, 2−r, …` from 0) only `r` and `1 − r` can lie in the support `|d| ≤ a`: `G = V(r) + V(1−r)`.  Also
`0 ≤ b < ½` and `a + b = c ≤ 1`, so `b ≤ 1 − a ≤ ½`.  Then on `[0, ½]`:
`r ≤ b`: `V(r) = 1/c`, `V(1−r) = 0` (as `1 − r ≥ 1 − b ≥ a`), `G = 1/c`;
`b ≤ r ≤ 1−a`: `V(r) = (a−r)/(sc)`, `V(1−r) = 0`, `G = (a−r)/(sc)`;
`1−a ≤ r ≤ ½`: both on the ramp, `G = (a−r)/(sc) + (a−1+r)/(sc) = (2a−1)/(sc) = P`.
The max/min expression reproduces the three pieces because `P ≤ (a−r)/(sc)` iff `r ≤ 1 − a`, and `P ≤ 1/c` iff
`c ≤ 1`.  The identity for P: `(c+s−1)(1+c) = s(1+c−s)` reduces to `c² + s² = 1`. ∎
This is the "joint bound for the pair" of the brief, and it is exact: an atom near `u = ±½` at small tilt sits on
two ramps whose values always sum to exactly P, independent of r.  Test [2]: closed form = periodised definition,
exactly (Fractions), at 3000 rational poses incl. poses at and 10⁻⁹ around every kink.

**Lemma 3 (box bound)** [proved].  Let `T = [θ_lo, θ_hi]`, `t = tan(θ/2)`, and `C = [cx_lo, cx_hi]`.  For atom k let
`r_k⁺ = max_{cx∈C} dist(x_k − cx, ℤ)` (exact: ½ if the u-interval contains a half-integer, else the larger endpoint
value).  Put `A = 1/c(θ_lo)`; `a_L = a(θ_lo)`; `σ_U = s c(θ_hi)` if `θ_hi ≤ 45°` (tested exactly as
`t_hi² + 2t_hi − 1 ≤ 0`) else ½; `P_L = A − s(θ_hi)/(c(θ_hi)(1 + c(θ_hi)))`.  Then for all `θ ∈ T ∩ (0°, 45°]`,
`cx ∈ C`:

    F(θ, cx) ≥ η + Σ_k w_k · max(P_L, min(A, (a_L − r_k⁺)/σ_U)).

*Proof.*  By Lemma 2 and `w_k ≥ 0` it suffices to bound each `G_θ(r) ≥ G_θ(r_k⁺)` (monotone in r) from below.
On [0°, 45°]: c is decreasing, so `1/c ≥ A`; `a = sin(θ + 45°)/√2` is increasing, so `a ≥ a_L`;
`sc = sin(2θ)/2` is increasing on [0°, 45°] and `≤ ½` always, so `sc ≤ σ_U`; `s/(c(1+c))` is increasing on [0°, 90°)
(numerator up, positive denominator down), so `P ≥ P_L`.  Finally `a_L ≥ ½ ≥ r` makes the ramp numerator
nonnegative, so `(a−r)/(sc) ≥ (a_L − r)/σ_U`.  max and min are monotone. ∎
At `θ = 0` the bound also holds for boxes with `θ_lo = 0`: there `A = 1`, every per-atom term is `≤ 1 ≤ G_0`.
In t-coordinates every bound is a rational function of `t_lo, t_hi` (`c = (1−t²)/(1+t²)`, `s = 2t/(1+t²)`).
Implementation: with `r₁ = a_L − σ_U A`, `r₂ = a_L − P_L σ_U` the per-atom term is A if `r ≤ r₁`, the ramp if
`r₁ < r ≤ r₂`, `P_L` otherwise (case check: `A ≥ P_L`; the ramp is `≥ A` iff `r ≤ r₁` and `≥ P_L` iff `r ≤ r₂`).
cx is kept on the integer grid `1/(N·2⁶⁴)`, so `r_k⁺` is an integer over that denominator and the thresholds are
compared after an exact floor; the sum is formed as `η + (A·S_plat + (a_L·S_ramp − Σ W r)/σ_U + P_L·S_P)/Q`.
Tests [4] (300 random boxes × 5 points) and [4b] (every 50th leaf of an actual certificate × 4 random points,
F evaluated from the *definition*): zero violations.

**Lemma 4 (θ = 0)** [proved].  `F(0, cx) = η + Σ_k w_k #{j : |x_k − cx + j| ≤ ½} ≥ η + Σ_k w_k = 1 + e`, because every
real number has an integer translate in `[−½, ½]` (two when it is ≡ ½).  So θ = 0 is trivial.  The checker still
evaluates `F(0, ·)` exactly at every breakpoint `x_k ± ½` and every midpoint between breakpoints (F(0, ·) is piecewise
constant and upper semicontinuous, so this is the exact min).

**Small tilt (the θ → 0 subtlety)** [proved].  Lemma 2 makes it disappear: `G_θ ≥ P(θ)` everywhere, and for
`0 < θ ≤ θ₁`, `P(θ) ≥ 1 − s(θ₁)/(c(θ₁)(1 + c(θ₁))) ≈ 1 − θ₁/2`.  So `F ≥ η + (1 + e − η)(1 − θ₁/2)`, which is ≥ 1 for
`θ₁ ≲ 2e` (≈ 0.5° at e = 0.005), with no condition on cx.  In the branch and bound this is just Lemma 3 on a box
with `θ_lo = 0`; no separate lemma or special-case code is needed (the shallowest leaves anywhere have θ-depth 4–5,
t-width `TMAX/16`–`TMAX/32`).  The per-image bound the brief warns about (ramps of the two images bounded
separately, giving ≈ 0 instead of ≈ 1) is never used.

## 3. The checker

Domain: `t ∈ [0, 107/256]` (`107/256 ≥ √2 − 1`, checked exactly: θ up to 45.37°).  The part above 45° is redundant
(Lemma 1) and not claimed: a box with `t_lo ≥ √2 − 1` (exact test `t_lo² + 2t_lo − 1 ≥ 0`) is skipped, a box
straddling 45° is bounded on its part `θ ≤ 45°` only (σ_U = ½, a_L = a(θ_lo)); the union of the other boxes still
covers `[0, √2−1]` since `√2 − 1` is irrational and never a box edge.  cx ranges over `[0, ½]` when the measure is
mirror-symmetric (`w_k = w_{−k}`, checked exactly; then `F(θ, −cx) = F(θ, cx)` since G is even), else over `[0, 1]`.  Depth-first branch and bound with exact Lemma 3 bounds; a
box is a leaf when its bound is ≥ 1.  Split rule (heuristic, does not affect soundness): halve the coordinate whose
collapse to the midpoint gains more bound.  Failing boxes at depth ≥ 24 have their centre (θ clipped to ≤ 45°)
evaluated exactly; F < 1 there stops the run with an exact counterexample.  The worse child is explored first.
Max depth 90 (then the min of F at 9 points of the box is reported).  Coverage is by construction (each box is replaced by its two halves).

**Symmetry of the LP solutions.**  The brief says the LP solutions are mirror-symmetric.  **They are not**: the stored
`runs/seam_dual/1d_N200_e005.npy` has `max |x_k − x_{−k}| = 0.073` (the LP optimum is not unique).  Symmetrising
`(x + x∘(−·))/2` is free: same g and total, and `F_sym(θ, cx) = (F(θ, cx) + F(θ, −cx))/2 ≥ min F`.  It even raised
the float min from 0.99975 (brief) to 0.99984.  `prepare` symmetrises, and `check` verifies the symmetry exactly
before halving cx.

**Preparing a measure** (`prepare`, float, untrusted): LP of `seam_1d.py` (θ on a 0.5° grid + tiny tilts, cx on the
`1/2N` grid) at fixed e; then cutting planes: a float scan finds, for θ on a 0.005° grid with local refinement, the
exact min over cx (F(θ, ·) is piecewise linear with kinks at `x_k ± b`, `x_k ± (1−a)` by Lemma 2, so the min is at a
kink) and adds the violated poses (and mirrors) plus the argmin pose for every θ on a 0.05° grid; up to 6 rounds,
keeping the round with the best predicted κ.  Then `W_k = ⌊10⁹ x_k⌋` (mirror-symmetric), and
`η = ⌈(1 − min F_float) + margin⌉` (denominator 10⁹).

## 4. Results

All [proved] (modulo the trust chain in §6), N = 200, Q = 10⁹, mirror-symmetric, θ = 0 check passed for each.
Exact values are in the JSON files; `e` and `g` below are rounded.  CPU: one core, Python Fractions, final code.

| N | e (certified) | g (certified) | κ = e/g | η | η/e | leaves | max depth | CPU | gridded LP g (FRIEDMAN §8) |
|---|---|---|---|---|---|---|---|---|---|
| 200 | 0.00101577 | 0.04589573 | **0.022132** | 1.59e-5 | 0.016 | 5 783 126 | 31 | 16 068 s | 0.0475 (N = 400) |
| 200 | 0.00251011 | 0.07990124 | **0.031415** | 1.02e-5 | 0.004 | 1 728 776 | 32 | 3393 s | 0.0808 |
| 200 | 0.00502130 | 0.12654188 | **0.039681** | 2.14e-5 | 0.004 | 505 029 | 30 | 649 s | 0.1273 |
| 200 | 0.01004007 | 0.20638961 | **0.048646** | 4.00e-5 | 0.004 | 295 664 | 29 | 210 s | 0.2069 |

(+ a few hundred boxes with θ ≥ 45° skipped per run; leaves-by-θ-depth histograms in the `.check.txt` files.)

* The headline: **κ = e/g is certified ≤ 0.0221 at g ≈ 0.046, ≤ 0.0315 at g ≈ 0.080, ≤ 0.0397 at g ≈ 0.127,
  ≤ 0.0486 at g ≈ 0.206**,
  i.e. the float "sublinear price" table of FRIEDMAN §8 survives exactly, losing only 0.2–3.4 % of g (the gridded LP
  is slightly infeasible in the continuum; the cutting planes fix that) plus η/e = 0.4–1.6 % of excess.
  Local exponent from the certified points: log(κ ratio)/log(g ratio) = 0.51 (e 0.0025→0.005) and 0.42 (0.005→0.01),
  and 0.63 (0.001→0.0025), consistent with `e ≈ C g^{1.5}` (κ ∝ g^{0.5}) [measured on 4 certified points].
* e = 0.0005 (N = 200): measure prepared (`seam1d_N200_e0.0005.json`: g = 0.028604, e = 0.000506, κ = 0.0177,
  η/e = 0.013; float-feasible only, [measured]); its exact check was still running at the time of writing
  (progress in `seam1d_N200_e0.0005.progress.log`, result in `.check.txt` when done).  Note the gridded LP says
  g = 0.0307 there: the continuum loses 7 %, the gap grows as e shrinks (the optimum has finer structure).
* The raw stored LP solution (`runs/seam_dual/1d_N200_e005.npy`, symmetrised, no cuts) was also certified
  (margin 1e-4, 66 543 leaves, 11 s, earlier code version; file not kept): g = 0.127296, e = 0.0052623, κ = 0.04134.
  Cutting planes lower g a little but cut the needed top-up from 1.6e-4 to 1.4e-6, which wins (κ 0.0397).
* Cost grows as e shrinks (smaller margin η − deficit at fixed η/e, finer atom structure).


## 5. Tests (`python3 search/seam_1d_exact.py test --cert …`)

Outputs: `seam_data/tests_N200_e*.txt`.  All passed (`TESTS PASSED 0`).
* [1] chord formula (reduced to [0°, 45°] by Lemma 1's symmetries) vs an independent polygon-clipping chord at
  20 000 random θ ∈ (−180°, 180°): max |diff| 1.1e-14.  Confirms formula, symmetries and the reduction.
* [2] Lemma 2 closed form vs the periodised definition `Σ_{|j|≤3} V(u+j)`, exact, 3000 rational poses (a fifth of them
  at or within 3·10⁻⁹ of a kink or of u = ½): 0 mismatches.
* [3] exact F (Fractions) vs `seam_1d.V` (floats) at 300 random poses: max |diff| 4.4e-16.
* [4] Lemma 3 bound ≤ F(definition) at 5 random points in each of 300 random boxes: 0 violations;
  [4b] (e = 0.005) every 50th certificate leaf (10 100 leaves) × 4 random points: F ≥ bound ≥ 1, 0 violations.
* [5] adversarial (all must be rejected; each was, with an exact witness pose where F < 1):

| cert | measure | rejected at (θ, cx) | exact F − 1 there |
|---|---|---|---|
| e = 0.005 | η = 0, every atom × (1 − 10⁻⁴) | 40.18°, 0.000122 | −4.1e-5 |
| e = 0.005 | atom x = 0.41 (+ mirror) lowered 1.14e-4 (F at the float-min pose 17.99°, 0.37 becomes 1 − 2.2e-4) | 44.998°, 0.1829 | −5.7e-6 |
| e = 0.005 | same, one-sided (symmetry check correctly fails → cx ∈ [0, 1]) | 40.18°, 0.99988 | −1.4e-5 |
| e = 0.01 | η = 0, atoms × (1 − 10⁻⁴) | 44.998°, 0.1015 | −3.6e-5 |
| e = 0.01 | atom x = 0.46 (+ mirror) lowered 1.29e-4 | 44.998°, 0.1015 | −1.1e-4 |
| e = 0.0025 | η = 0, atoms × (1 − 10⁻⁴) | 36.17°, 0.01135 | −8.1e-6 |
| e = 0.0025 | atom x = 0.935 (+ mirror) lowered 1.08e-4 | 41.49°, 0.000427 | −1.6e-5 |
| e = 0.0025 | same, one-sided | 36.17°, 0.98865 | −3.2e-5 |
| e = 0.001 | η = 0, atoms × (1 − 10⁻⁴) | 44.88°, 0.08197 | −9.5e-5 |
| e = 0.001 | atom x = 0.40 (+ mirror) lowered 1.01e-4 | 44.88°, 0.08197 | −1.5e-4 |
| e = 0.001 | same, one-sided | 44.88°, 0.91809 | −6.3e-5 |
| any | all weights halved | θ = 0 check: min F = 0.50 | |

  "At the right place": the perturbed measures *are* violated at the targeted float-minimum pose (exact F − 1 ≈ −2e-4
  there, printed by the test), but depth-first search stops at the *first* violated box it meets, which is often a
  different near-tight pose (near θ = 45° or the seam, cx ≈ 0).  Every reported witness is an exact rational pose with
  F < 1, which is what matters.  Caveat found while building the test: removing η alone (deficit 1.4e-6) is a
  borderline measure with whole 2D patches where F ≈ 1 ± 10⁻⁶ (LP-tight poses); rejecting it needs boxes of size
  ~10⁻⁷ and ran > 1 h without finishing, so the test uses a 10⁻⁴ violation instead.


## 6. What is NOT checked / trust chain

* Trusted: Python's `fractions.Fraction` and integer arithmetic; the ~120 lines `trig`, `rdist`, `BnB.tq`, `BnB.lb`,
  `BnB.run`, `check_theta0`, `Measure` in `seam_1d_exact.py`; and the hand proofs of Lemmas 1–4 above.
  Floats are used only to choose the measure (and in the tests); the split rule and search order are heuristics
  that cannot make a false claim pass.
* The 1D reduction itself (a y-invariant profile in the plane receives exactly `F(θ, cx)` per unit square:
  `ρ × Leb_y` of the square = `∫ρ(dx) V_θ(x − cx)`) is FRIEDMAN.md §8's [proved] remark, not re-checked here.
* Coverage of the (θ, cx) domain is by construction of the recursion, not re-verified from a stored leaf list
  (no leaf list is stored; re-running `check` regenerates the whole tree deterministically).
* Optimality is not certified: g is a certified *achievable* seam at excess e (an upper bound on the price κ(g)),
  not the LP optimum.  The continuum LP value is slightly below the gridded one (the certified g is 0.2–3.4 % below the gridded LP g, 7 % at e = 0.0005).
* Nothing here addresses items (b)–(d) of FRIEDMAN.md §8 (wall row, slow variation, top transition), nor the corner.
