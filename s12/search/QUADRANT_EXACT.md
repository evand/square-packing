# Exact certificate for an (R = 2, w = 2) fixed-profile family: s(k² − 3) = k for all k ≥ 6 (task quadrant-exact-w2, 2026-09-28/29)

Brief: `tasks/quadrant-exact-w2/README.md`.  Background: `search/QUADRANT.md` (the float quadrant LP).
Code (all new, nothing pinned was edited): `search/qx2_lp.py` (LP re-solve), `search/qx2_exact.py` (exact
projection, θ = 0 face, box cover writer, family file), `search/qx2_family_check.py` (independent family ↔ cover
check), `search/qx2_zm.py` (zero-margin checker: imports the pinned `zm_mixed.py`/`zeromargin.py`), `search/qx2_zm_test.py`
(tests), `search/qx2_tools.py` (float diagnostics), `search/qx2_reduction.py` (independent reduction check, written
by a second agent that did not read QUADRANT.md first).  Data: `search/qx2_data/` (the exact family, the box cover,
the exact solution).  Runs: `runs/qx2_*` (gitignored).

Labels: **[proved]** = proved here on paper and/or certified by exact (`Fraction`) computation, with the program and
sha named; **[measured]** = float; **[heuristic]** = interpretation.

## 0. Answer

RUN_STATUS_PLACEHOLDER

* **Exact family** [proved, exact arithmetic]: a profile π of 52 axis-parallel uniform pieces per period (pitch 1/5)
  plus Lebesgue on `y ∈ [9/5, 2]`, and a corner module ν of 36 pieces plus Lebesgue on `[9/5, 2]²`, all masses rational
  (denominator ≤ 10¹²), **σ = 0 exactly**, mirror/diagonal symmetric, and
  **`D = 423621306389/500000000000 = 0.847242612778 > 3/4`**, so the box measures have total
  `μ_k([0,k]²) = k² − 4D = k² − 3.388970451112 < k² − 3` for every `k ≥ 6`.
  Files (`search/qx2_data/`): `L4_k02_family.txt` (sha256 `3f4771d3…0768`), `L4_k02_box7.txt` (the box `k = 7`,
  FORMAT.md v1, sha256 `c0a67509…694b`), `L4_k02_exact.json` (`8931e71d…0fb5`).  `qx2_family_check.py` rebuilds μ₇ from
  the family file with its own code: the identical multiset of 800 segments + the Lebesgue square `[9/5, 26/5]²`,
  total `5701378693611/125000000000 = 49 − 4D`.
* **What the certificate claims**: every closed unit square `Q ⊂ [0,7]²` has `μ₇(Q) ≥ 1`.  It is split as
  (i) the θ = 0 face, **Lemma Z**, certified exactly (`qx2_zm.py axis`: 900 one-sided limit corners, minimum exactly 1,
  188 exactly tight; and on the quadrant model, `qx2_exact.py axis`: 907 corners, min exactly 1);
  (ii) θ ∈ (0°, 53°] on the D4 fundamental domain, by the box checker `qx2_zm.py` (sha256 `7f6554be…d96b`), whose leaves
  are zm_mixed's primitives (pinned `zm_mixed.py` `ee3e2915…60ac` and `zeromargin.py` `640fe453…86ab`, unedited) plus new
  exact primitives with proofs in §4: **LEB** (Lemma U), **CAP** (Lemma K), **EXACT** (Lemma E with E′/E″, the
  zero-margin one), **AXIS** (boxes whose admissible poses are θ = 0 only).
* **Reduction** [proved, §2; independently re-derived by a second agent]: μ₇ valid ⇒ the quadrant measure μ_Q
  valid ⇒ μ_k valid for all integers `k ≥ 6`; with the accounting identity and the dilation argument this gives
  `s(k² − 3) ≥ k`, hence `s(k² − 3) = k` for every `k ≥ 6`.  (The box `k = 2R + 2 = 6` would *not* suffice for the
  converse step; `k = 7` does.)
* **Zero margin, and where it sits** [proved + measured, §1.4]: σ = 0 forces the axis-parallel wall squares
  `[t,t+1]×[0,1]` to have mass exactly 1 (dilated-grid argument).  In the exact solution the θ = 0 face is tight on two
  one-parameter families per wall (`cy = ½⁺`, `cy = 3/2⁺`, `cx ≥ 5/2`) and inside U; the tile-germ limits (θ → 0 at a
  centre where lines lie on edges) are tight at several double germs in the corner (e.g. `Q → [1,2]²`).  Off θ = 0 the
  LP was made to grow linearly in θ (tilt margin 0.2·min(θ, 3°)·(1 − area(Q ∩ U))).  Lemma E is exact at the tight
  points (no box loss) and certifies zero margin with first-order growth; Lemma U/K handle U; zm_mixed does the rest.
* **A genuine invalidity was caught** [proved, §1.5]: the first exact solution (`L2_k02`, D = 0.847254…) is *not*
  valid: at the double germ `(3/2 + 10⁻¹², 3/2 + 2·10⁻¹³, θ = 2·10⁻¹²)` its exact mass is `0.99993599 < 1`.  The float LP
  missed it because quadrant_lp's containment tolerance 10⁻⁹ over-counts chords by ~10⁻⁹/θ at tiny tilts; the box
  checker refused exactly the two boxes containing the bad germs.  Fixed by strict tolerance and exact germ-limit rows
  (`L4_k02`); that cover is kept as `qx2_data/REJECTED_L2_k02_box7.txt`.  **The old float optimum of QUADRANT.md
  (D = 0.945) was measured with the same lenient tolerance**, so its "0.14 % residual" is optimistic at tiny tilts.
* **Not proved / trust**: everything certified is Python `Fraction` arithmetic in two programs (`qx2_zm.py` importing
  the pinned `zm_mixed.py`/`zeromargin.py`; `qx2_zm.py axis` / `qx2_exact.py axis` for θ = 0).  There is no second
  independent implementation of the θ > 0 part: `zmx2` (Rust) refuses polygons, and the Lebesgue square cannot be
  removed.  Lemma E is new code (~800 lines) tested by adversarial random covers, zero-margin slope tests, exact
  germ scans and rejection covers (§5), and it refused the genuinely invalid L2 cover exactly where it is invalid;
  but it has not been audited by a second person.  zm_mixed's polygon Lemma S(b) was never used in a shipped
  certificate before; it is audited in §4.5.

## 1. The exact-friendly family (steps 1–2)

### 1.1 LP changes (`qx2_lp.py`, imports `quadrant_lp.py`)

* **Lebesgue layer**: μ is exactly Lebesgue on `U = [a, ∞)²`, `a = w − 1/5 = 9/5` (profile Lebesgue on `y ∈ [a, w]`,
  corner module Lebesgue on `[a, R]²`, no atom in the interior of U).  With `R = w` the layers and L tile U, which is
  convex, so a square inside U has mass exactly 1 by a one-line lemma (Lemma U).
* **Lines only, pitch 1/5**: no points, no cells (points are usc step functions, cells have convex area pieces; both
  obstruct exact zero-margin work).  Profile: horizontal pieces at `y ∈ (1/5)ℤ ∩ (0, a]`, vertical pieces at phases
  `(1/5)ℤ`; corner module: horizontal/vertical pieces on the `1/5`-grid.
* **Exact θ = 0 rows**: every one-sided limit corner of the θ = 0 breakpoint grid (Lemma Z) is a row with EXACT
  coefficients (`qx2_exact.ExactModel.coeffs`); HiGHS primal/dual feasibility tolerance 10⁻¹⁰.
* **Exact germ-limit rows** (`--germs`): for every centre c0 on the 1/10-grid of the window (both orders, i.e. both
  tilt directions) and every `t` on the 1/5-grid of `[−6/5, 6/5]²` (at a wall only `t ≥ 1`), the EXACT limit as
  u → 0⁺ of the row at `(c0 + u t, u)` (`ExactModel.germ_coeffs`: each chord end `N(u)/F(u)` converges in the extended
  reals; the captured fraction is continuous in the ends).  At a germ the blow-up `t ↦ lim f(c0 + u t, u)` is a sum of a
  function of `t_x` and one of `t_y`, each piecewise linear with breakpoints on the 1/5-grid (ends `c0x ∓ t_y/2`,
  `c0y ± t_x/2` crossing breakpoints of the 1/5-lattice), constant beyond `|t| ≥ 1`; so these rows imply every germ limit.
  13,112 distinct rows; added to the LP only when violated.
* **Strict containment** in the float model: quadrant_lp's `TOL = 10⁻⁹` is replaced by 0 (see §1.5).
* **Tilt margin**: a row at angle θ asks for `1 + κ·min(θ′, θ₀)·(1 − area(Q ∩ U))`, `θ′` = distance of θ to 0 mod 90°,
  `κ = 0.2`, `θ₀ = 3°`; the factor `(1 − area)` is necessary (squares inside U have mass exactly 1).
* **Vertex solutions**: HiGHS dual simplex (no IPM), crossover on.
* Row generation as in QUADRANT.md §2 (lattice, random, near-lattice, germs, polish), warm-started from the 40k rows
  of `runs/quad_R2_w2_g`.

### 1.2 D at each stage [measured, float LP]

| run (`runs/qx2_*`) | model | rounds | D_LP | oracle min (last) | note |
|---|---|---|---|---|---|
| `quad_R2_w2_g` (QUADRANT.md) | pitch 0.1, points + cells, IPM | 18 | 0.945057 | 0.9997 (heavy 0.99856; box k=7 0.99884 per the second agent) | starting point; invalid by ≥ 0.1 % |
| `band_p1` | band only, layer a = 1.8, pitch 0.1 | 41 | m_v^band = 1.29909 | 0.999994 | (1.433 without the layer) |
| `p1_pts`, `p1_nopts` | pitch 0.1 + layer | 0 | — | — | first simplex solve did not finish in 1 h; abandoned |
| `p2_nopts` | pitch 0.2, layer, cells, no points | 22 | **0.910540** | 0.99999994 | |
| `p2_lines` | pitch 0.2, layer, lines only | 28 | **0.907834** | 0.99999992 | |
| `p2_nopts_k02` | + tilt margin κ = 0.2, θ₀ = 3° | 14 | 0.872090 | (margin rows) | |
| `p2_lines_k02` | lines only + tilt margin | 20 | 0.847263 | (margin rows) | |
| `L2_k02` | + exact θ = 0 rows, tolerances 10⁻¹⁰ | 3 | 0.8472544250 | 0.99999995 | **invalid** at two corner germs (§1.5) |
| `L2_k03` | κ = 0.3, θ₀ = 4° | 0 | infeasible | — | the margin demand has a ceiling |
| `L3_k02` | + strict containment | 6 | 0.8472435 | 0.99999998 | still invalid at a germ limit (−2·10⁻⁶): rows at θ = 10⁻⁵ do not pin θ → 0 |
| `L4_k02` | + exact germ-limit rows | 3 | **0.8472426128** | 0.99999995 | **the certified one** |

So the price of exactness is `0.945 → 0.847`: layer + pitch 0.2 cost 0.035, dropping cells 0.003, the tilt margin 0.06,
the germ rows 10⁻⁵.  D stays 0.097 above 3/4.

### 1.3 The exact solution (`qx2_exact.py project`)

Rather than re-solving the float basis exactly (tried first: the tilted basis rows are float poses whose exact
snapping moves the coefficients by up to 5·10⁻⁴ at near-tangent poses, and the resulting exact vertex was garbage),
the exact solution is the *projection* of the float one onto the exact identities that matter at zero margin:

* every θ = 0 limit corner and every germ-limit row whose float value is `< 1 + 10⁻⁷` becomes an exact equality
  `= 1` (326 + 34,503 rows), plus `σ = 0`; the support (50 variables) is kept; the system has **rank 5**;
* the free variables are rounded to `1/10¹²`, the five pivots solved exactly.

Result: all masses nonnegative, max deviation from the float LP `1.2·10⁻¹²`, `D = 423621306389/500000000000`,
`σ = 0`.  The tilted rows are not used: they carry the margin, and validity there is the checker's job.  The trade
"D for margin in the corner module" (brief step 2) was **not needed**: the tilt margin is bought in the LP instead,
spread over band and corner, with σ = 0 kept exactly.

### 1.4 Where the zero margin is

* **Forced** [proved]: with σ = 0 the band squares `A_t = [t,t+1]×[0,1]` (θ = 0) have mass exactly 1 for all
  non-integer t: disjoint closed squares `A_t` and `B_t = [t,t+1]×[1+ε, 2+ε]` tile a column of the strip up to ε, the
  band has mass w = 2 per period, so `μ(A_t) + μ(B_t) ≤ 2 + O(ε)` while both are ≥ 1.
* **θ = 0 face** [proved, exact enumeration]: tight exactly on `{cy = ½⁺, cx ≥ 5/2}` and `{cy = 3/2⁺, cx ≥ 5/2}` (and
  the diagonal images), plus squares inside U; every other one-sided limit corner is ≥ 1.001.
* **Germ limits** [proved, exact scan `qx2_germscan.py`, 1,922 centres × 625 offsets at u = 10⁻¹²]: minimum exactly
  `1.000000000000`, 0 poses below 1; tight at the double germs `(3/2, 3/2)`, `(3/2, 23/10)` and others (the LP rows).
* **Near θ = 0** [measured, `qx2_tools.py`]: `min (f − 1)/(θ(1 − area(Q∩U)))` over a 0.02-grid of centres at
  θ = ±10⁻⁵ … ±3·10⁻² is `0.20 – 0.23` (the LP's κ), never negative; along the wall family 0.4 – 0.6 per radian; at the
  corner germs ≈ 0.7.  So the tight set of the pose space is the two families at θ = 0, the germ columns on them, the
  corner germs (θ → 0 only) and U-inclusion — all of which Lemmas Z, E, U handle exactly.

### 1.5 The germ violation that the checker found (L2_k02)

The first full checker runs on the `L2_k02` cover left uncertified boxes only at
`[3/2, 3/2 + 1/640]² × [0°, 0.112°]` and `[3/2, 3/2 + 1/640] × [23/10, 23/10 + 1/640] × [0°, 0.112°]`, at depth 18,
with Lemma E failing at every depth.  The float minimum there (QUADRANT's evaluator) was `1.0005`.  An exact
evaluation (`zm_mixed.exact_mass` and, independently, the quadrant model `ExactModel.coeffs`) at
`(3/2 + 10⁻¹², 3/2 + 2·10⁻¹³, u = 10⁻¹²)` gives `0.999935986 < 1`, and the value is the same for u from 10⁻¹² to 10⁻⁵:
a **tile-germ limit** of the square `[1,2]²` (lines `x = 1, 2` and `y = 1, 2` on its edges, blow-up offset
`t = (1, 1/5)` in u-units) is below 1 by 6.4·10⁻⁵.  The LP never saw it: `quadrant_lp.chord` counts a chord end up to
`TOL = 10⁻⁹` outside the square, i.e. a length error `≈ 10⁻⁹/sin θ` per chord end, which at the oracle's germ angle
θ = 10⁻⁵ over-counts by ~2·10⁻⁴ (float model 1.0001 where the exact value is 0.99994).  Strict containment (`L3`)
exposed it; the LP then satisfied the germ rows at θ = 10⁻⁵ with the margin 2·10⁻⁶ but let the limit θ → 0 sit at
`1 − 2·10⁻⁶` (the first-order term pays for it); the exact germ-limit rows (`L4`) close this.  The rejected cover is
kept as `qx2_data/REJECTED_L2_k02_box7.txt`: a natural rejection test that the checker passed.

## 2. The reduction (step 4) [proved]

Setting (QUADRANT.md §1.1, with the layer): `μ_Q = ν + π|_{[R,∞)×[0,w]} + π′|_{[0,w]×[R,∞)} + λ|_{L}`,
`L = {x > w, y > w} \ [0,R]²`, and λ also on the layers; here `U = [a,∞)²` carries exactly λ plus boundary pieces.
Box: `μ_k` = ν at the four corners (images under the symmetries of the square), π on `[R, k−R]` along each wall,
λ on `[a, k−a]²`.

1. **Accounting.**  `k² − μ_k([0,k]²) = 4D + 4σ(k − 2R)` with `D = R² − ν(C) − m_v`: `[R, k−R]` is `k − 2R`
   half-open periods plus the closed end, whose cross-section has mass `m_v`; the Lebesgue part is
   `(k − 2a)² − 4(R − a)²` on top of the layers, and the pieces tile `[a, k−a]²` up to null sets.  With σ = 0 the
   saving is `4D` for every `k ≥ 2R + 2`.  (Exact check for k = 7: §0; float check k = 6..12: table below.)
2. **Localisation (μ_Q valid ⇒ μ_k valid, k ≥ 2R + 2).**  μ_k is D4-invariant: `x ↦ k − x` maps the bottom band to
   itself because π is mirror symmetric and k is an integer, and maps the BL corner module to the BR one because ν is
   diagonal symmetric (`s = r ∘ swap`).  On `H = [0, k−R)²`, μ_k = μ_Q (the right/top bands start at `k − w ≥ k − R`,
   the other corner modules at `≥ k − R`).  A unit square has x-extent `≤ √2`, so if `k − 2R > √2` every square in
   `[0,k]²` lies in `g(H)` for some `g ∈ D4`.  Smallest k: `2R + 2 = 6`.
3. **One box suffices (μ₇ valid ⇒ μ_Q valid)**, for `k = 2R + 3 = 7`.  Let `S ⊂ Q` be a closed unit square.
   (a) `S ⊂ H = [0,5)²`: `μ_Q(S) = μ₇(S) ≥ 1`.  (b) `x_max(S) ≥ 5 > y_max(S)`: then `x_min ≥ 5 − √2 > R`; on
   `{x > R}` μ_Q is `π + λ|_{y ≥ a}`, invariant under integer x-shifts that stay in `{x > R}`; the admissible shifts
   form an open interval of length `k − 2R − W ≥ 3 − √2 > 1` (`W` = x-extent), so some integer shift puts S inside
   `(R, k−R) × [0, k−R) ⊂ H` (the closed ends are avoided on purpose: the lines `x = R`, `x = k − R` carry extra
   module mass).  (c) symmetric.  (d) both ≥ 5: `S ⊂ {x, y > R + 1}` ⊂ U: `μ_Q(S) = 1`.
   With `k = 6` the interval has length `2 − W < 1` for tilted squares, and the argument fails (a shift crosses a
   corner module or lands on a seam line).
4. **Consequence.**  Upper bound `s(k² − 3) ≤ k` (tiling).  Lower bound: if `n = k² − 3` unit squares pack in
   `[0, s]²` with `s < k`, dilate by `k/s > 1`; the concentric closed unit squares in the dilated squares are pairwise
   disjoint (each lies in the interior of its dilated square) and lie in `[0,k]²`, so
   `n ≤ Σ μ_k(T_i) = μ_k(∪ T_i) ≤ μ_k([0,k]²) = k² − 4D < k² − 3 = n`, a contradiction.  Validity at the single
   side k is all that is used (the closed-square semantics alone would not do: squares of a packing may share
   boundary points).  This is FORMAT.md's argument; `lean/Sqpack/MixedMeasure.lean` proves it for arbitrary measures.

**Independent check** (second agent, derived before reading QUADRANT.md, then compared): (1)–(4) agree, including
the σ term and `k₀ = 2R + 2` for (2); it found that (3) needs `k = 2R + 3` (QUADRANT.md §8.1 says the same), and it
supplied the dilation step (4), which QUADRANT.md did not spell out.  Numerics (`qx2_reduction.py`, its own box
builder and evaluator, cross-checked against `quadrant_lp.Model.contrib` to `5·10⁻¹⁵` on 5000 poses), on the exact
masses of `L4_k02` (the second agent's script, run again on the final solution):

| k | total | k² − total | 4D | oracle min (≈ 0.44–0.61 M poses: random, walls, corners, near-lattice ±10⁻⁷, polish) |
|---|---|---|---|---|
| 6 | 32.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 7 | 45.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 8 | 60.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 9 | 77.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 10 | 96.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 11 | 117.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |
| 12 | 140.611029549 | 3.388970451 | 3.388970451 | 1.000000000 |

(The minima 1 are squares inside U.  This float oracle, like the LP's, would not have seen the L2 germ violation of
§1.5: it is 6·10⁻⁵ deep but lives on a set of poses of width ~10⁻¹² around θ = 0.)  On the *old* float solution the same
script found `0.99884` at k = 7, below the old run's own oracle (0.9997): the old solution was more invalid than its log
said.

## 3. The certificate

`qx2_data/L4_k02_box7.txt` (FORMAT.md v1: 0 points, 800 segments on 60 lines, 1 polygon = `[9/5, 26/5]²` with mass =
area).  Claim: every closed unit square `Q ⊂ [0,7]²` has `μ(Q) ≥ 1`.  By the D4 invariance (checked exactly by
`Cover.symmetric_d4`) it suffices to treat centres in `[0, 7/2]²` and `θ ∈ [0°, 45°]`; the roots are zm_mixed's D4 roots
`[0, 7/2]² × u ∈ [0, 1/2]` (u = tan(θ/2); 1/2 is 53.13°), pitch 1/10, 8 u-bins: 9,800 roots.

* θ = 0: Lemma Z (§4.1).
* each root box is subdivided (zm_mixed's rule, depth ≤ 18); a box is a leaf if one of:
  AXIS (after zm's `clip_bin` its admissible poses are θ = 0 only: Lemma Z), LEB (Lemma U), CAP (Lemma K),
  zm_mixed's PIECE / ADM / P1 / MIX / SPLIT (with the polygon handled by its Lemma S(b), audited in §4.5),
  EXACT (Lemma E; certifies `u ∈ (u₀, u₁]`, the θ = 0 part again by Lemma Z when `u₀ = 0`: leaf kind EXACT0; for a box
  straddling 45° it certifies the part `θ ≤ 45°`: EXACT45), **SYM** (`θ₀ ≥ 45°`, exact test `u₀² + 2u₀ − 1 ≥ 0`: the
  diagonal reflection maps the box's poses to poses with `θ′ = 90° − θ ≤ 45°` and centre in `[0, 7/2]²`, which the
  roots cover with `u ≤ tan 22.5° < 1/2`; the same argument covers the `θ > 45°` part of an EXACT45 box).

## 4. The new lemmas

Notation: pose `(c, u)`, `u = tan(θ/2)`, `C = 1 − u²`, `S = 2u`, `N = 1 + u²` (cos θ = C/N, sin θ = S/N);
`Q(c,u) = c + R_θ[−½,½]²`; a point p is in Q iff `|X|, |Y| ≤ ½`, `X = (p − c)·(C,S)/N`, `Y = (p − c)·(−S,C)/N`.
The measure is `μ = λ|_U + Σ_ℓ μ_ℓ`, U = `[a, b]²` (`a = 9/5`, `b = 26/5`), μ_ℓ a piecewise-uniform measure on an
axis-parallel line ℓ with cumulative mass `G_ℓ(t)` (non-decreasing, piecewise linear; breakpoints = segment ends and
density changes, including the two ends of the support).

### 4.1 Lemma Z (θ = 0 face)

*At θ = 0, `m(c) = μ([cx−½, cx+½]×[cy−½, cy+½])`.  Let Γ be the set of all `v ± ½`, v ranging over the breakpoints
of every `G_ℓ`, the positions of all lines, and a, b.  On every open cell of the grid Γ × Γ, m is a polynomial of degree
≤ 1 in each of cx, cy; and m is upper semicontinuous.  Hence `inf m` over centres in a grid rectangle is the minimum
of the one-sided limits of m at the corners of the open cells.*

Proof.  A horizontal line `y = η` contributes `1[η ∈ [cy − ½, cy + ½]]·(G(cx + ½) − G(cx − ½))`; on an open cell the
indicator is constant (η ± ½ ∈ Γ) and `G(cx ± ½)` is affine (its breakpoints ± ½ are in Γ).  Vertical lines likewise;
U contributes `|[cx ± ½] ∩ [a,b]|·|[cy ± ½] ∩ [a,b]|`, a product of affine functions on a cell.  A multilinear
function on a rectangle takes its extrema at the corners, so its infimum over the open cell is the least corner limit.
Upper semicontinuity: for `c_n → c`, `lim sup Q₀(c_n) ⊂ Q₀(c)` (closed squares) and μ is finite, so
`lim sup m(c_n) ≤ m(c)`; thus the value at any point of a cell's boundary is ≥ the limit from any adjacent cell. ∎

Implementation: `qx2_zm.py axis FILE` on the box cover, centres `[½, 7/2]²` (D4), exact `Fraction`s, limits
computed with one-sided inclusion of the lines lying on the moving edges; and `qx2_exact.py axis` on the quadrant
model (window `cx ≤ R + 11/4`, `cy ≤ min(cx, R + 7/4)`, a superset of QUADRANT.md's window).  Both: minimum exactly 1.
The *rejection* covers (§5.3) all fail here (0.999999, 0.999998, 0.99996).

### 4.2 Lemma U (inside the Lebesgue square) and Lemma K (caps)

**Lemma U.**  *If every admissible Q of a box lies in U, μ(Q) ≥ λ(Q) = 1.*  (μ ≥ λ|_U because the polygon's mass equals
its area; `u_square` checks this exactly.)  Test: `cx₀ − ŵ/2 ≥ a`, `cx₁ + ŵ/2 ≤ b` and the same in y, with ŵ =
zeromargin's upper bound of `cos θ + sin θ` on the bin.

**Lemma K.**  *Let Q ⊂ `{x ≤ b, y ≤ b}` be a closed unit square with `cy ≥ a` (resp. `cx ≥ a`) whenever Q meets
`{y < a}` (resp. `{x < a}`), and let `d_y = a − min_y Q`, `d_x = a − min_x Q`.  If `d_y ≤ ρ_y` and `d_x ≤ ρ_x`, where
`ρ_y` is a lower bound of the density of the cover's pieces on the line `y = a` over the x-range of Q (and ρ_x
likewise), then μ(Q) ≥ 1.*

Proof.  `λ(Q ∩ U) ≥ 1 − λ(Q ∩ {y < a}) − λ(Q ∩ {x < a})`.  Width lemma: the horizontal slices of a square, as a
function of height, increase from the lowest vertex to the second, are constant to the third, decrease to the top;
the third vertex is at height `cy + |cos θ − sin θ|/2 ≥ cy ≥ a`, so all slices below `a` are no longer than the slice
at `a`, of length `c_y = |Q ∩ {y = a}|`; hence `λ(Q ∩ {y < a}) ≤ c_y·d_y`.  The pieces on `y = a` carry
`≥ ρ_y c_y` inside Q; the lines `y = a`, `x = a` meet in one point.  So
`μ(Q) ≥ 1 − c_y d_y − c_x d_x + ρ_y c_y + ρ_x c_x ≥ 1`. ∎
Test: `d ≤ a − cy₀ + ŵ/2`; `ρ` = the exact minimum density of the line's profile over a range containing the chord
`Q ∩ {y = a}` (0 if not fully covered): `[cx₀ − ŵ/2, cx₁ + ŵ/2]`, or, when `u₀ > 0` and `d ≤ min(sin θ₀, cos θ₁)` (the cap
is the triangle at the lowest vertex BL, legs along the unit edges `(c, s)` and `(−s, c)`), the tighter
`[x_BL − (s/c)d, x_BL + (c/s)d]` with `x_BL = cx − (c − s)/2` bounded over the box (for `x = a`: the triangle at the
leftmost vertex TL, `[y_TL − (c/s)d, y_TL + (s/c)d]`, `y_TL = cy + (c − s)/2`).  This matters at squares touching the
boundary of U with a vertex, where the crude range reaches a gap of the line near `(a, a)`.  Test (`captangent`):
2,421 CAP-certified boxes at such tangencies, 24,210 exact masses, 0 below 1.

### 4.3 Chord ends are affine in the centre

For a horizontal line `y = η` (parameter t = x) and u ∈ (0, 1), the point `(t, η)` is in Q iff
`max(t₁, t₂) ≤ t ≤ min(t₀, t₃)` with

| option | type | `F·t =` (F ∈ {C, S}) |
|---|---|---|
| t₀ (X ≤ ½) | up | `C·t = C cx + S cy + N/2 − ηS` |
| t₁ (X ≥ −½) | lo | `C·t = C cx + S cy − N/2 − ηS` |
| t₂ (Y ≤ ½) | lo | `S·t = S cx − C cy + ηC − N/2` |
| t₃ (Y ≥ −½) | up | `S·t = S cx − C cy + ηC + N/2` |

and for a vertical line `x = ξ` (t = y): `S t = S cy + C cx + N/2 − ξC` (up), `S t = S cy + C cx − N/2 − ξC` (lo),
`C t = C cy − S cx + N/2 + ξS` (up), `C t = C cy − S cx − N/2 + ξS` (lo).  (Direct from `|X|, |Y| ≤ ½`; `qx2_zm_test.py
opts` compares the formulas with an independent chord computation at 3000 random rational poses × 2 lines: 0
mismatches.)  So `μ_ℓ(Q) = max(0, G(r↑) − G(r↓))`, `r↑ = min(up options)`, `r↓ = max(lo options)`, and
`μ_ℓ(Q) ≥ G(r↑) − G(r↓)` always.

### 4.4 Lemma E (exact minimisation on a pose box; zero margin allowed)

Box `B = [cx₀,cx₁]×[cy₀,cy₁]×[u₀,u₁]`, `0 ≤ u₀ < tan 22.5°`, `u₁ ≤ ½`.  Claim certified:
`μ(Q(c,u)) ≥ τ₀ + λu` for every admissible pose with `u ∈ (u₀, u_e]`, `u_e = min(u₁, tan 22.5°)` (all polynomial checks
run on the larger `[u₀, u₁]`; every bound that uses monotonicity of cos, sin, w in u is taken on `[u₀, u_e]`, where
cos decreases and sin, w increase, with the rational bounds `cos ≥ 7071/10⁴`, `w ≤ 14143/10⁴`, `cos − sin ≥ 0` when
`u₁` is past 45°) (τ₀ = 1, λ = 0 in the certificate; other values
only in tests).  The admissible centres at fixed u form the rectangle `R(u) = [x_lo(u), cx₁] × [y_lo(u), cy₁]`,
`x_lo = cx₀` or `w(u)/2` (whichever holds on the whole bin — else the lemma is not applied), likewise y; `cx₁ + w/2 ≤ 7`
is checked.

**E′ (the Lebesgue part).**  For the line `x = a` let `d = a − min_x Q = a − cx + w/2`, `s = sin θ`, `c = cos θ`.
`λ(Q ∩ {x < a}) = g(d)` with `g = 0` (d ≤ 0), `d²/(2sc)` (0 ≤ d ≤ s), `(d − s/2)/c` (s ≤ d ≤ c), `1 − (s+c−d)²/(2sc)`
(c ≤ d ≤ s + c), 1 beyond (the vertical slice width from the leftmost vertex is `t/(sc)`, then `1/c`, then decreasing).
Let `ĝ` = g on `(−∞, s]` continued linearly by `(d − s/2)/c` for `d ≥ s`.  Then `ĝ ≥ g` everywhere (on `[c, s+c]`
the continuation is the tangent at c of the convex function `1 − g`; beyond, `ĝ(s+c) = 1 + s/(2c) ≥ 1`), ĝ is convex
(C¹ at 0 and s), and the parabola `d²/(2sc)` is ≥ ĝ everywhere.  Also, if `d ≥ c` on the whole box,
`λ(Q ∩ {x ≥ a}) = (s + c − d)²/(2sc)` on `[c, s+c]` and 0 beyond, a convex function of d, ≥ its tangent `T_{d*}` at any
`d* ≤ 1 ≤ s + c`.  With the same for `y = a`, and `Q ⊂ {x, y ≤ b}` (checked):

    λ(Q ∩ U) ≥ Φ_x + Φ_y − 1,   Φ = 1 − ĝ(d)  ("std")  or  Φ = T_{d*}(d)  ("tan", when d ≥ c on the box),

and Φ = 1 for a line Q never crosses; if Q never meets U, λ(Q ∩ U) = 0.  Every Φ is a **concave** function of c at
fixed u (d is affine in c).

**E″ (near the corner of U: exact inclusion–exclusion).**  `λ(Q ∩ U) = 1 − g_x − g_y + A_LL` exactly, with
`A_LL = λ(Q ∩ {x < a, y < a})`; the bound `Φ_x + Φ_y − 1` drops `A_LL`, which near `(a, a)` is not small (0.64 at
`Q ≈ [1,2]²`) and, near the tight germs at `(3/2, 3/2)` and `(3/2, 23/10)`, is of the same order as the growth of the
mass in θ, so it must be kept exactly.  Two configurations, each certified on the whole box before use
(`corner3_ok`, `xcut_ok`), cover the corner germs; outside them the plain `Φ_x + Φ_y − 1` is used.  Notation: BL, BR,
TR, TL = lowest, rightmost, top, leftmost vertex (θ ∈ (0°, 45°)), `z = βC − αS` with `α = (a − x_BL)N`,
`β = (a − y_BL)N`, and `quad(α, β) = (2Cαβ + S(β² − α²))/(2CN²)` (`= αβ` at θ = 0).

* **corner3** (frame at BL: `ξ` along BL→BR, `η` along BL→TL, unit square, Jacobian 1).  Conditions:
  `x_BL ≤ a ≤ x_BR`, `y_TL ≥ a`, `(a, a)` not to the right of the right edge (`X(a,a) ≤ ½`).  Then `Q ∩ {x<a, y<a}` =
  `{ξ, η ∈ [0,1], Cξ − Sη ≤ α, Sξ + Cη ≤ β}` and, integrating over ξ,
  - if `z ≤ 0` (the point of the bottom edge above `x = a` is at height ≥ a): the whole y-cap lies left of `x = a`,
    `A_LL = g_y`, so `λ(Q ∩ U) = 1 − g_x`;
  - if `z ≥ 0`: `A_LL = quad(α, β)`, in both the triangle y-cap (BR above a) and the strip y-cap (BR below a) — the two
    integrals (`∫₀^{ξ*}` with `ξ* = (αC+βS)/N²`, resp. `∫₀^{α/C} + ∫_{α/C}^{ξ*}`) give the same quadratic; so
    `λ(Q ∩ U) = 1 − g_x − g_y + quad(α, β)`.
  (In the strip y-cap, `z > 0` automatically.)  Concave cell functions: on `z ≥ 0` cells `1 − ĝ_x − ĝ_y + Mc(quad)`,
  on `z ≤ 0` cells `1 − ĝ_x`; the line `z = 0` is added to 𝓛, and at a vertex whose side is not certified the minimum
  of the two is used (it is below both cell functions).
* **xcut** (frame at TL: `ξ′` along TL→TR, `η′` along TL→BL).  Conditions: `x_BL ≥ a` (the x-cap is the triangle at
  TL) and `y_TL ≥ a`.  With `α′ = (a − x_TL)N`, `β″ = (a − y_TL)N ≤ 0`, `z′ = α′C + β″S`: `A_LL = (z′⁺)²/(2SCN²)`
  (the part of the x-triangle below `y = a`, a triangle with base `z′/(SC)` on the left edge and height `z′/N²`), so on
  `z′ ≥ 0`: `−g_x + A_LL = quad(α′, β″)`; on `z′ ≤ 0`: `A_LL = 0`.  Cell functions `Φ_y + Mc(quad(α′, β″))` and
  `Φ_y − ĝ_x`, line `z′ = 0` added, minimum at uncertified vertices.
* **McCormick** `Mc`: `αβ ≥ Aβ + Bα − AB` with `(A, B)` the values at the box corner `(cx₁, cy₁)` (then
  `(α − A)(β − B) ≥ 0`) or at `(cx₀, cy₀)` (then `(A − α)(B − β) ≥ 0`) — both quantities decrease in cx and cy —, both
  anchors tried (the second is exact at a tight germ at the box's low corner); `Sβ² ≥ S(2β*β − β*²)` (tangent, β* at
  mid-`cy`); `−Sα²` kept.  Loss `≤ h_x h_y + u h_y²/4`, zero at the anchor.  No `1/u` factors: the combinations
  `−g_x + A_LL` and `−g_y + A_LL` in which the triangle areas `∝ 1/(SC)` cancel are what is bounded.

Tests (`qx2_zm_test.py corner`, `cuts`): at 9,472 (corner3) and 3,237 (xcut) random rational poses of random boxes
satisfying the conditions, the identities hold exactly against polygon clipping (`z ≤ 0`: area = 1 − g_x; `z ≥ 0`:
area = 1 − g_x − g_y + quad; `A_LL = (z′⁺)²/(2SCN²)`, `−g_x + A_LL = quad(α′, β″)` for `z′ ≥ 0`), and the code's terms
(cap bounds, McCormick with both anchors, the z-functions) are on the right side of the true values: 0 failures.

**E (the lemma).**  *Fix u.  Let 𝓛(u) consist of the four sides of R(u) and the following lines in the centre plane
(each affine in c, coefficients polynomial in u):*
* *for every line ℓ within reach, every option k and every breakpoint b of G_ℓ that is **bad** for k —
  density increasing across b if k is an up option, decreasing if k is a lo option —: `{t_k = b}`;*
* *for every pair (up option, lo option) of ℓ: `{t_up = t_lo}` (where the chord can become empty);*
* *for every std-mode cap: `{d = 0}` and `{d = s}`;*
*omitting only lines certified not to meet R(u) for any u of the bin.  Then `f(c) := Σ_ℓ μ_ℓ(Q(c,u)) + Φ_U(c)` is
concave on every open cell of the arrangement of 𝓛(u) in R(u), continuous on R(u), and `f ≤ μ(Q)`; hence
`min_{R(u)} μ(Q) ≥ min f` over the vertices of the arrangement, i.e. over the points of R(u) where two lines of 𝓛(u)
meet.*

Here `Φ_U` is the concave minorant of `λ(Q ∩ U)` of E′, or, in the corner3 / xcut regimes, the cell-wise concave
minorant of E″ (with the line `z = 0` resp. `z′ = 0` in 𝓛); the regimes are used only when their conditions are
certified on the whole box.

Proof.  `Φ_U` is concave on each cell (E′, E″).  Take a line ℓ and an open cell K.  (i) The chord `[r↓, r↑]` is either nonempty on all
of K or empty on all of K: `r↑ − r↓` can only change sign where some `t_up = t_lo`, a line of 𝓛.  If empty,
`μ_ℓ = 0` on K.  (ii) If nonempty, `μ_ℓ = G(r↑) − G(r↓)`.  On K, `r↑` never equals a bad up-breakpoint b: at such a
point the minimising option `t_{k*}` would equal b, i.e. the point would lie on `{t_{k*} = b} ∈ 𝓛`.  So `r↑(K)` lies in
an open interval between consecutive bad up-breakpoints, on which G has only good kinks (density non-increasing) and
is therefore concave: `G = min_p A_p` over its affine pieces there, and `G(r↑) = min_p A_p(min_k t_k) = min_{p,k}
A_p(t_k)` (A_p non-decreasing), a minimum of affine functions of c.  Likewise `G` is convex on the range of `r↓`, so
`−G(r↓) = min_{p,j}(−A_p(t_j))`.  So `μ_ℓ` is concave on K; the sum is concave.  Continuity: for u > 0 no axis line is
parallel to an edge of Q, so chord lengths are continuous in c; Φ is continuous.  A continuous function concave on
each open cell attains its minimum over the closure of a cell (a convex polygon) at a vertex; the cells cover R(u). ∎

**How it is certified** (`Exact.certify`).  For every pair of lines of `𝓛` (as polynomial forms `A cx + B cy + K`) the
vertex is `v(u) = (X/Δ, Y/Δ)`, Δ = `A₁B₂ − A₂B₁`, X, Y polynomials.  On a sub-bin (bisection of the bin, depth ≤ 32,
≤ 96 sub-bins per vertex):
1. the sign of Δ must be certified constant (Bernstein coefficients; at `u₀ = 0` after dividing out the power of u);
   parallel pairs (Δ ≡ 0) have no vertex;
2. if `v(u)` is certified outside R(u) on the whole sub-bin (one side-constraint strictly negative) it is skipped;
3. otherwise a lower bound `LB(u) ≤ f(v(u), u)` is built as a sum over lines of `min(up alternatives) −
   max(lo alternatives)`: an alternative is `A_p(t_k)` for an option k and a piece p that t_k can occupy on the
   sub-bin (options are dropped only if certified dominated — `t_k ≥ t_{k′}` for an up option kept, `≤` for lo —
   and the occupied pieces are those between breakpoints certified to bound `t_k`; `G(t) ≥ min_p A_p(t)` and
   `G(t) ≤ max_p A_p(t)` hold for t in the union of the pieces); lines whose float chord is empty at the sample points
   are dropped (bound 0, always valid); plus the Φ terms of E′ (for a cap, the regime `d ≤ 0`, `d ≥ s` is used only
   where certified, else the parabola, an upper bound of ĝ everywhere);
4. `LB(u) − τ₀ − λu`, written as one rational function `P(u)/(C^α S^β N^γ Δ^δ)`, is ≥ 0 on the sub-bin if all Bernstein
   coefficients of `±P` (the sign from Δ^δ) are ≥ 0 — **every combination** of alternatives (the minimum over
   alternatives of a sum is the minimum over combinations), after pruning alternatives certified dominated; or, if the
   vertex leaves R(u) inside the sub-bin, an **S-procedure** certificate `LB − τ − ν·g ≥ 0` with `ν ≥ 0` and g the violated
   side constraint (then `LB ≥ τ` wherever `g ≥ 0`, i.e. wherever v(u) is a vertex in R(u)).
At `u₀ = 0` the numerator test on the closed `[0, u₁]` covers `u ∈ (0, u₁]` (the denominators are > 0 there); a
numerator vanishing at 0 (the zero-margin families, `f → 1` as θ → 0) is fine as long as the Bernstein coefficients
are ≥ 0, which they are when f grows at first order.  Everything is `Fraction` arithmetic; floats only choose options,
pieces, orderings and the S-procedure multiplier (any choice is sound).

Why this is exact where zm_mixed is not: zm_mixed's bounds are box-wide (Lemma S cores, Lemma T minima, Bernstein
ratio bounds), losing `O(box)` or `O(Δu)`; at a pose where the true mass is exactly 1 on a face of the box they can
never reach 1.  Lemma E evaluates the mass *at the vertices, as functions of u*, so at a tight vertex the bound is the
exact value.  Measured: on the failing box `[3, 3 + 1/80] × [½, ½ + 1/80] × [0, 0.112°]` zm_mixed's piece bound is
0.99926 (the line `y = 1` lies on the top edge at θ = 0, so it only gets Lemma S), Lemma E certifies it in 0.2 s.

### 4.5 Audit of zm_mixed's polygon bound (Lemma S(b)), which the certificate uses

zm_mixed disables polygons in `--cert-mode` because they had not been audited.  `qx2_zm.py` runs zm_mixed's
`piece_bound` with the polygon branch on (Corollary T′ is unreachable: no points).  The code (`cond_region`,
`clip_halfplane`, `clip_convex`, `poly_area`, `piece_bound` lines 1104–1125) implements ZM_MIXED.md Lemma S(b):
* `cond_region`: for each Lemma-B bound choice χ, the Bernstein coefficients of `G_{k,χ}(u; p)` are affine in p
  (built from the values at `base`, `base + e_x`, `base + e_y`; Fact 0), each `β_i ≤ 0` is a half-plane
  `bX x + bY y + (b0 − bX bx − bY by) ≤ 0`, clipped exactly (Sutherland–Hodgman keeps `a x + b y + c ≤ 0`); a degenerate
  coefficient (bX = bY = 0) empties the region iff b0 > 0.  The hull over the choices is taken; it is ⊂ `C_k` because
  `C_k` (points satisfying inequality k at every admissible pose) is an intersection of half-planes, hence convex.
  The frame `box ± 0.7072` contains every Q, so clipping to it loses nothing.
* `piece_bound`: `K = frame ∩ K₀ ∩ … ∩ K₃` (each CCW from `convex_hull`), `area(P ∩ K)` by `clip_convex` (valid for
  convex CCW polygons; `check_simple_convex` validates P), times `w/area(P)`; `K ⊂ ∩ Q` over the box, so this is a lower
  bound of the polygon's mass at every admissible pose.  Polygons and lines are disjoint pieces, so the bounds add.
* Tests: zm_mixed's own selftest `[3]` (496 checks); here `qx2_zm_test.py capleb` checks LEB/CAP leaves against exact
  masses (7,380 poses, 0 below 1).  (No leaf dump of the certificate run was stress-tested; the independent float
  oracles of §2 and the exact germ scan of §1.4 are the evidence that the cover is valid beyond the checker's verdict.)
Nothing wrong was found.  The polygon bound is sound; it is simply lossy (`O(box)` on the part of ∂Q inside U), which
is why the straddle region is slow.

## 5. Tests (all exact unless stated; `qx2_zm_test.py`, final code unless stated)

### 5.1 Unit tests
* `opts`: chord-option formulas vs direct exact chords, 6000 checks, 0 mismatches.
* `poly`: Bernstein nonnegativity vs dense exact evaluation, 400 random polynomials, 0 false positives.
* `capleb`: 557 LEB + 58 CAP boxes on the real cover, 7,380 exact pose masses, 0 below 1.
* `corner`, `cuts`: the E″ identities (corner3, xcut) at 12,709 exact poses, 0 failures (§4.4).

### 5.2 Lemma E soundness and completeness
* `exact` (real cover, random boxes in the wall family, the B family, the corner germs `(3/2, 3/2)`, `(3/2, 23/10)`,
  `(23/10, 3/2)`, `u₀ ∈ {0, …}`): certifying `τ = (min of 60 exact sampled masses) + 10⁻⁹` would be a bug: 0 of 516
  boxes (earlier code versions and the final one); completeness at `τ = min − 10⁻⁴`: 57–87 %.
* `adv` (random D4-symmetric covers with 80–112 random piecewise-uniform lines and a Lebesgue square; boxes away from
  U, straddling its boundary, or near its corner; `u₀ ∈ {0, 1/20, 1/10}`): `τ = float minimum (3000 samples + polish)
  + 10⁻⁷` must not certify: **0 of ≈ 1,400** boxes over 17 seeds.
* `gap`: the largest certified τ vs the float minimum over the box: median gap 7·10⁻⁵, 80 % below 1.3·10⁻⁴, max 8·10⁻³
  (the float minimum is itself only an upper bound), 0 unsound (40 boxes, final code).
* `slope` (zero margin at u → 0): boxes `[x, x+h] × [½, ½+h] × [0, u₁]` and at `cy = 3/2`: the claim `mass ≥ 1 + λu`
  with λ = (float minimum growth over the box) + 0.02 must be refused: refused in 140/140; with λ = growth − 0.05
  certified in 94/140.
* `leafstress`: the checker re-run with a leaf dump on 5 roots at the tight places (wall family, corner germs `(3/2, 3/2)`, `(3/2, 23/10)`, `(23/10, 3/2)`, left B family), 40 EXACT/EXACT0/PIECE leaves, 1,200 exact masses at random admissible rational poses of the leaves (corners, faces, tiny u favoured): minimum `1.0000048`, 0 below 1.
* A natural rejection (§1.5): on the invalid `L2_k02` cover the checker left uncertified exactly the boxes containing
  the two bad corner germs; every other root verified.

### 5.3 Rejection covers (`qx2_zm_test.py holes`)
| cover (from the final L4 cover) | construction | Lemma Z | box checker (wall family, cx ∈ [2.9, 3.3], cy ∈ [½, 0.6], u ≤ 0.01) |
|---|---|---|---|
| R_scale | all segment masses × (1 − 10⁻⁶) | min 0.999999: VIOLATED | NOT VERIFIED: 32 uncertified boxes = every u₀ = 0 box on the family (no EXACT0 leaf) |
| R_piece | line y = 1, x ∈ [3, 3.2] (+ D4 images) lighter by 10⁻⁴ | 0.999998: VIOLATED | NOT VERIFIED: 32 uncertified boxes = every u₀ = 0 box on the family (no EXACT0 leaf) |
| R_seam | seam x = 3, y ∈ [0.8, 1] (+ images) lighter by 10⁻⁴ | 0.99996: VIOLATED | NOT VERIFIED: 32 uncertified boxes = every u₀ = 0 box on the family (no EXACT0 leaf) |
| REJECTED_L2_k02 | the first exact solution (germ violation 6.4·10⁻⁵) | min 1: passes (the violation is at θ → 0⁺ only) | NOT VERIFIED exactly at the two bad germs (runs `L2_k02_diag*`) |

(On the earlier L2 cover the same three holes were rejected by both Lemma Z and the box checker.)

## 6. The run

RUN_RESULTS_PLACEHOLDER

## 7. Open items

* A second, independent implementation of the θ > 0 part (Lemma E with the Lebesgue square) — zmx2 cannot take
  polygons; the Lean kernel verifier handles segments (P/S/T/L lemmas) but not polygons or Lemma E.
* An outside audit of Lemma E and of `qx2_zm.py` (≈ 900 lines), and of the zm_mixed polygon path.
* The same pipeline for `k² − 4` (R = 3, w = 3, QUADRANT.md: float D ≈ 1.15): needs D > 1, i.e. the tilt margin and
  the layer must cost less than ≈ 0.15 there; untried.
