# Hat (continuous piecewise-linear) line densities vs uniform pieces in the qx2 LP (2026-10-01/02)

Code: `search/qx2_pl.py` (model + comparison driver), `search/qx2_pl_test.py` (tests).  Nothing pinned was edited
(`quadrant_lp.py`, `qx2_lp.py`, `qx2_exact.py` and `qx2_kmax.py` are imported).  Runs: `runs/qpl_*`, `runs/qplk_*`
(gitignored; `log.txt`, `hist.json`, `sol.npz`, `support.txt`, `binding.txt`, `final.json`).
**Everything here is [measured] (float LP, float oracle) or [heuristic].  Nothing is proved.**

## 1. Hypothesis and method

Uniform line pieces make μ(pose) kink wherever a chord end crosses a density breakpoint.  The hypothesis: continuous
piecewise-linear line densities ("hats") remove those kinks, so they should (a) raise D at fixed tilt margin κ and/or
(b) raise the largest feasible κ.

* **Ramp element**: a segment [p,q] on an axis-parallel line with linear density, normalised to the element mass x.
  Its CDF is F(t) = ((t−p)/(q−p))² (up) or 1 − ((q−t)/(q−p))² (down).  The captured fraction is
  F(min(hi,q)) − F(max(lo,p)), where [lo,hi] = `Q.chord` with strict tolerance.
  **Hat at node t_i** = an up-ramp on [t_{i−1}, t_i] plus a down-ramp on [t_i, t_{i+1}], each piece of mass x
  (peak density 2x/h, so the density is continuous).  A half-hat sits at each end of a line's allowed range.
* The line families are the same as qx2_lp's lines-only model at hc = hp = 0.2:
  * corner lines y = c with their diagonal mirrors, over the same x-ranges as the uniform pieces;
  * profile horizontal lines with periodic hats in phase.  The orbit is {φ, 1−φ}, and the mirror maps the hat at φ
    to the hat at 1−φ (an up-ramp becomes a down-ramp).  Every ramp lies inside its own period cell [j, j+1],
    j ≥ R, so the band is still π restricted to x ≥ R.
  * profile vertical lines with hats in y on [0, a].
  The obj and σ coefficients count the pieces' masses exactly as the uniform code does.  The Lebesgue layer on U and
  the margin are XModel's, unchanged.
* **`both`** (an extra, beyond the brief): uniform pieces and hats together, a superset of either basis.  It
  separates "continuity" from "a different basis of the same size".
* **Equal terms**: both element types run with float rows only.  There are no germ rows.  The exact axis rows are
  replaced by float θ = 0 rows at every combination of {g ± 1e-7, cell midpoints} in x and y over the breakpoint
  grid g: one-sided corners, edge midpoints and cell midpoints, 1976 poses.  Everything else follows qx2_lp: initial
  lattice pitch 0.1 at 19 angles plus near-lattice poses, `Q.oracle` (10⁶ random poses, pitch-0.02 lattice,
  near-lattice and germ poses at θ = ±1e-5 and ±1e-3, polish of the worst 2000), up to 4000 rows per round, HiGHS
  simplex at 1e-10 tolerances.
* **Stopping rule**: oracle min ≥ 1 − 1e-6 in 2 consecutive rounds, or D flat to 1e-4 over 5 rounds, or 40 rounds.
  The oracle min is min(μ − margin).

## 2. Tests (`qx2_pl_test.py`, all pass, 9.4 min on one core)

| test | result |
|---|---|
| brute force: each ramp split into N = 2000 uniform pieces of mass x·ΔF and evaluated by the uniform code (an `XModel` with fixed hs/vs pieces); random x; random, tiny-θ (10⁻⁹…10⁻², and near 90°), wall- and corner-resting, near-lattice ±1e-7 and germ poses | R=w=2, h=0.5, a=1.5 (phase node 0.5 is self-mirror): 18,035 poses, max error 7.2e-8.  Same with a=1.7 (off the node grid): 18,035 poses, 5.9e-8.  Production R=w=3, h=0.2, κ=0.1: 4,000 poses, 1.3e-7.  Expected discretisation error is about x/(4N²) per cut ramp. |
| LP rows A·x + const vs the oracle path contrib(P, x) | ≤ 2.3e-14 |
| σ and D bookkeeping vs direct geometry for a random x: σ-row LHS = mass in [R+1, R+2) × [0,w]; m_v = mass on the line x = R+1; D = R² − μ([0,R]²) + m_v | agree to 1e-12 (3 configurations) |
| symmetry: diagonal (cx,cy,θ) ↔ (cy,cx,90°−θ), and band mirror x ↦ 2R+5−x | ≤ 2.5e-13 |
| all-equal hats c vs the uniform model with every piece 2c | identical everywhere (≤ 3e-14), same D and σ (the half-hats make the end cells uniform too) |
| `both` = uniform + hat − Lebesgue (masses and D) | 9e-15 / 0 |

## 3. Results [measured]

R = w = 3, δ = 0.2 (a = 2.8), lines only, hc = hp = 0.2, nvar: uniform 308, hat 326, both 634.

| κ (θ0) | uniform D | hat D | both D |
|---|---|---|---|
| 0 | **1.13158** (min 0.99989, 23 rnd, flat) | 1.09871 (0.99936, 13, flat) | **1.15495** (0.99964, 23, flat) |
| 0.05 (0.5°) | **1.12110** (0.999992, 19, flat) | 1.09007 (0.99958, 16, flat) | **1.14795** (0.99954, 21, flat) |
| 0.1 (3°) | **1.04176** (0.99996, 33, flat) | 0.94222 (0.99996, 25, flat) | **1.07555** (0.99969, 24, flat) |
| 0.2 (3°) | LP fails at round 5 (HiGHS kUnknown; D 1.03 → 0.75 at rounds 3–4) | **infeasible** at round 4 | **infeasible** at round 4 |

R = w = 2 (a = 1.8), same settings:

| κ (θ0) | uniform D | hat D |
|---|---|---|
| 0 | 0.90787 (0.9999999, 13, flat) | 0.87871 (0.99990, 11, flat) |
| 0.05 (0.5°) | 0.90645 (0.99999, 13) | 0.87699 (0.99985, 11) |
| 0.1 (3°) | 0.89036 (0.99998, 14) | 0.86286 (0.99859, 11) |
| 0.2 (3°) | 0.84746 (0.99971, 17) | 0.80809 (0.99822, 13) |

The uniform float-row values match today's exact-axis-row plus germ-row runs: 1.132 / 1.121 / 1.043 at R=w=3, and
κ = 0.2 infeasible.  At R = w = 3 the float-only comparison therefore loses nothing for the uniform model.

**Max κ** (`qx2_pl.py kmax`, qx2_kmax's LP with D free, θ0 = 3°, R = w = 3): uniform κ* ≤ 0.1135 (34 rounds,
oracle min 0.983 at stop); hat κ* ≤ 0.1126 (20 rounds, 0.967).  Both runs stopped on "κ flat" with the oracle
still finding violations, so these are upper bounds.  κ = 0.1 converged for both, so κ* ∈ [0.10, 0.113] for both.
The `both` kmax was stopped at round 4 (κ = 0.154 and falling) for time.

Run times: hat runs took 4–30 min each at R=3.  Another session's seam jobs were sharing cores 13–15, so the hat
oracle rounds ran 15–230 s against 14–50 s for uniform.  Hat contrib does about 2× the element work.  The kmax runs
took about 55 min each.

## 4. Reading [heuristic]

* **Hats alone do not help.**  At equal pitch, the hat basis loses 0.03 in D at κ = 0 and 0.05 (both R = 2 and 3),
  and loses 0.10 at κ = 0.1 (R = 3).  κ* is unchanged (≈ 0.11 for both bases), and κ = 0.2 is infeasible at R = 3
  for every basis.  The uniform optima use density jumps: line densities that switch on or off at square-edge
  positions, matched to the wall squares.  A continuous density at pitch 0.2 cannot copy a jump.
* **The union helps D but not the slope.**  `both` gains +0.023 (κ = 0), +0.027 (0.05) and +0.034 (κ = 0.1) over
  uniform, but κ = 0.2 is still infeasible.  The gain looks like extra degrees of freedom (comparable to refining
  the pitch) rather than kink removal.
* **Where the tight poses sit** (`binding.txt`, dual weight).  In every optimum, hat or uniform, 26–52 % of the dual
  weight is on θ = 0 rows and 0.3–0.7 is on near-θ = 0 rows.  The θ = 0 rows are the wall squares at
  cy = 0.5, 1.5, 2.5 with cx ≈ 3.9–4.5, plus the corner square (0.5, 0.5).  The near-θ = 0 rows are the germ poses
  around those wall squares, at θ = 1e-3 rad (0.057°) and θ ≈ 0.25–0.4°, with centres offset 1e-4…3e-3 from
  lattice + ½.  The θ = 0 binding poses sit at grid corners in all but a few cases (one interior pose in total).
  So hats remove the breakpoint kinks, but the binding structure does not change.  The cap comes from the σ = 0
  zero-margin wall families (QUADRANT_EXACT.md §1.4) and their first-order tilt growth, not from kinks.  With κ > 0
  the near-θ = 0 dual share rises (0.56–0.72), and those wall-germ rows are what cap κ for both bases.
* **Caveats.**
  * The hats are not continuous across x = R where the corner meets the profile: a jump is allowed there, as in the
    uniform model.
  * Hats were tried only at pitch 0.2.  A hat basis at pitch 0.1 (twice the variables) was not run.
  * Hat models cannot use the exact axis or germ rows (qx2_lp's exact machinery assumes uniform pieces).  An exact
    certificate for hats would need new exact code: the θ = 0 mass is piecewise quadratic, not multilinear.
  * For uniform at R=3 and κ=0.2 the LP status was kUnknown rather than a clean infeasibility certificate.
    Together with κ* ≤ 0.114, this reads as infeasible.
