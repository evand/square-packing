# SOS probe: can exact Positivstellensatz certificates pass the zero-margin wall at `s(6) = 3`?

*2026-10-03.  Task `tasks/sos-probe/README.md`.  Code: `search/sos_probe_enc.py` (encoding + numeric sanity),
`search/sos_probe_sdp.py` (SOS programs), `search/sos_probe_exact.py` (exact rational rounding + verification),
`search/sos_probe_facial.py` (facial-reduction attempt at `eps = 0`), `search/sos_probe_run.sh`,
`search/sos_probe_sweep.sh`, `search/sos_probe_batch.sh`, `search/sos_probe_report.py`.  Venv `runs/sos_venv/`
(cvxpy 1.9.3, Clarabel 0.11.1, SCS, sympy, numpy, scipy, python-flint 0.9.0).  Raw results are in the
gitignored `runs/sos_probe/` (`res/*.json` one line per SDP, `exact.jsonl`, `facial.jsonl`); every number
below is quoted from there.*

Labels: **[proved]** = exact rational identity checked in exact arithmetic (modulo the encoding, §1);
**[measured]** = floating-point SDP or local optimisation; **[heuristic]** = inference or extrapolation.

## 0. Verdict

1. **Degree.**  Degree 4 is enough at `T = 3 - eps` for `eps >= 0.05` on the 5-square cores (all 32 sign
   orthants at `eps = 0.1`), but fails by `eps = 0.02`.  Degree 6 reaches `eps = 0.001` on 2-parameter angle
   families, and on 3 parameters it holds at the one `eps` tried (`0.01`).  With 4 free angles it barely holds
   at `eps = 0.01`: the margin drops 75x.  The full 10-angle system does not fit in memory at degree 6 (over 30 GB).
   No degree reached `eps = 0`, the exact `T >= 3`, on any system with two or more free angles
   [measured].
2. **Exactness.**  Wherever the SDP has margin (`T = 3 - eps`), the float solution rounds to an exact
   rational identity at once (Peyrl–Parrilo projection, exact LDL).  This is done for the full 10-angle row
   core in one orthant at `T = 2.9`, and for the 2-parameter row/pinwheel families at `T = 2.99` and
   `2.999` [proved].  At `eps = 0` the SDP has no interior: one-parameter families reach `gamma = 3.0000000`
   in floats at degree 4.  One facial-reduction step (the zero-set point plus the plateau slacks) leaves no
   interior (`mu = -6e-12`, `-3e-10`), so no exact `T >= 3` was obtained [measured].
3. **Locality.**  The row-of-three core alone is **false** off the axis-parallel stratum (minimum `T` is
   `2.44`).  So is every 4-square subsystem tried.  The smallest true subsystems have **5 of the 6 squares**
   [measured].  Certificates do **not** transfer between nearby separation types: the intersection of two
   types that differ only in their diagonal pairs is false (`T = 2.47`, `2.00`) [measured, float witnesses].
   Sign orthants must be split (`2^n` pieces): without the split even two squares fail at degree 6 [measured].

**Recommendation: record as a dead end for `s(12)`.**  For `s(6)`, a full proof by this formulation needs a
different formulation, not more degree.  The obstruction is obstruction X (`notes/n12-gap.md` §4.11) moved
from box width into degree: the certifiable `eps` shrinks with degree, and the zero-margin limit needs
multi-level facial reduction that depends on the second-order contact.  The `n = 12` multiplier is
prohibitive (§6).

## 1. Encoding

Variables `T`, and per square `x_i, y_i, c_i = cos theta_i, s_i = sin theta_i` with `c_i^2 + s_i^2 = 1` and
`theta_i in [-45, 45 deg]` (`c_i >= |s_i|`).  This covers every orientation mod `90 deg` exactly once.  Edge
normal 0 (resp. 1) stays within `45 deg` of `x` (resp. `y`), so a typed direction is honestly "near x" or
"near y".  With `[0, 90)` instead, a square at `89 deg` turns an "x-separation" into a y-separation.  That
made sub-systems collapse to `T = 2`, so the first encoding was discarded.  `(c, s)` replaces the brief's
`u = tan(theta/2)`; it describes the same set and keeps every constraint **quadratic**, where the `u`-form has
degree 3 walls and degree 5 pair rows.

* walls: `x_i - (c_i +- s_i)/2 >= 0`, `T - x_i - (c_i +- s_i)/2 >= 0`, the same in `y` (half-width `(c+|s|)/2`;
  in a fixed sign orthant one sign suffices);
* pair `(i, j)` typed `(o, k, sg)` (normal `k` of square `o`, sign `sg`, `q` the other square):
  `n.(c_q - c_o) - 1/2 - cosD/2 -+ sinD/2 >= 0`, two rows encoding `|sin D|`, with `cosD = n_{o,0}.e_q`,
  `sinD = n_{o,0}.f_q` (separating-axis theorem; `sos_probe_enc.u_form_check`: `1.1e-16` against vertex
  projection).
* Every 3x3-cell pair is typed, far pairs included.  Dropping far pairs makes the 6-square systems false
  (`T = 2.83`, `2.69`).

By scaling, `s(6) = 3` is equivalent to "`T >= 3` on every type's closed feasible set".

Systems (cells `(row, col)` of the `3x3` tiling; "low" = the normal of the left/bottom square; diagonal
pairs typed `x` unless stated):

| name | cells | |
|---|---|---|
| `ROW_A` | row `(0,0),(0,1),(0,2)` + `(1,1),(2,0),(2,2)` | wall-to-wall row of three, `S6_LOCAL` `j = 3` shape |
| `PINWHEEL` | `3x3` minus the main diagonal | `S6_SKELETON` §5 |
| **R5** | `ROW_A` minus `(2,0)` | minimal true row subsystem (§4) |
| **P5** | `PINWHEEL` minus `(2,0)` | minimal true pinwheel subsystem (§4) |

**Certificate forms.**  `--mult angle`: every multiplier is a polynomial in the angle variables only.  For
fixed angles the system is an LP in the centres, so such a certificate is an LP dual `w(theta) >= 0` that is
SOS in `theta`.  Monomials are reduced mod `s_i^2 = 1 - c_i^2`, so there are no free multipliers.  The
generators are the rows, plus products of an angle-only generator with any row (`--prod`), whose multipliers
are scalar or SOS of degree 2 (`--pdeg 1`).  The **infeasibility form** at fixed `T = 3 - eps` is
`-lam = sum`, with `lam` maximised subject to trace `<= 1`; `lam > 0` is a certificate.  The **bound form**
is `T - gamma = sum`, with a `(3 - T)` generator for the angle-dependent `T`-weight.  "Order `d`" means
degree `2d` in the angles.  An early Lasserre/Putinar test with multipliers in all 25 variables gave `1.85`
instead of `2` on **two** squares at degree 4, and 30–50 GB per 5-square run.  The angle form is the one used
throughout.

## 2. Literature (short, from memory, not re-checked)

[heuristic] Lasserre/SOS hierarchies on packing appear in *density* bounds (de Laat, Oliveira, Vallentin;
Dostert, Guzmán, de Oliveira Filho, Vallentin for translative packings).  Those certificates are made rigorous
by interval or rational post-processing at a positive margin.  Optimality proofs for finite circle packings
in a square (Markót, Csendes) are interval branch-and-bound, not SOS.  Exact rounding of strictly feasible
SOS: Peyrl–Parrilo (2008); for degenerate (zero-minimum) cases, Kaltofen–Li–Yang–Zhi (rational-function SOS,
Gauss–Newton) and Magron–Safey El Din (RealCertify), with facial reduction per Permenter–Parrilo.  Degree
bounds for Putinar grow polynomially in `1/eps` (Nie–Schweighofer; Baldi–Mourrain).  That is the shape
measured in §3.  No SOS or Positivstellensatz treatment of *square packing in a square* is known to me.  A
web search found none.

## 3. Measurement 1 — degree vs `eps`  [measured]

`lam(eps)` at `T = 3 - eps`, orthant `+` on every square (the coherent-sign orthant; at `eps = 0.1` it is the
hardest or nearly the hardest of the 32, see below).  A dash means not run.

**Full systems** (all angles free), order 2, `--prod --pdeg 1`, `sigma_0` split over square triples
(`--csp0 3`; lossless at `eps = 0.1`: `1.697e-3` vs `1.708e-3` dense):

| system | angle vars | SDP (rows / Grams) | `eps = 0.1` | `0.05` | `0.02` | `0.01` | s / solve |
|---|---|---|---|---|---|---|---|
| R5 | 10 | 2,991 / 305 (max 25) | `1.33e-3` | `1.5e-4` | `1.9e-6` | — | 25–40 |
| P5 | 10 | 2,991 / 305 | `1.51e-3` | `1.8e-4` | `1.7e-6` | `4e-9` | 25–120 |
| R6 (`ROW_A`) | 12 | 5,813 / 440 | `3.12e-3` | `1.17e-3` | `3.5e-4` | `1.5e-4` | 120–170 (+50 build) |
| P6 (`PINWHEEL`) | 12 | 5,813 / 440 | `2.58e-3` | `4.8e-4` | `2.0e-5` | `1.1e-6` | 115–175 |

At degree 4, `lam` falls roughly like `eps^3` and is gone by `eps ~ 0.02` (R5, P5) or `0.01` (P6).  The sixth
square helps the row core (R6 is still at `1.5e-4` at `eps = 0.01`).  Order 3 (degree 6) on the full R5:
`M = 18,691` monomial rows and 60 Grams up to 63.  Clarabel's KKT factorisation asked for more than 30 GB
and died (two attempts); SCS at `1e-9` was not competitive (>10 min on the degree-4 instance).

**Angle families** (squares with equal labels share one angle, `0` = exactly axis-parallel).  These are
restrictions of the full system, so the degree they need is a **lower bound** on the full system's degree:

| family (R5 / P5) | order 2 | order 3 | order 4 |
|---|---|---|---|
| any 1-parameter family (`a,a,a,a,a`, `a,a,a,0,0`, `a,0,0,0,0`, `a,a,0,0,0`, `0,a,a,a,a`) | bound `gamma = 3.0000000` | 3.0 | 3.0 |
| R5 `a,a,a,b,b` (row / context) `eps = .1/.03/.01/.003/.001` | `3.0e-3 / 2.1e-4 / 3.3e-5 / 1.2e-6 / 4e-11` | `8.6e-3 / 2.0e-3 / 6.4e-4 / 2.1e-4 / 6.3e-5` | `1.2e-2 / 2.5e-3 / 8.5e-4 / 2.3e-4 / 1.3e-4` |
| P5 `a,a,a,b,b` | `7.7e-3 / 1.9e-3 / 5.2e-4 / 1.3e-4 / 4.0e-5` | `1.1e-2 / … / 8.2e-5` | `1.4e-2 / … / 1.5e-4` |
| R5 `a,a,b,c,c` (3 params), `eps = 0.01` | — | `5.4e-4` (74 s) | — |
| R5 `a,b,c,d,d` (4 params), `eps = 0.01` | — | `8.6e-6` (1,500 s, 9 GB, `M = 6,187`) | — |

So: 1 parameter → degree 4 is exact in floats.  2 parameters → degree 4 dies like `eps^3`; degree 6 is
linear in `eps` (`lam/eps ~ 0.06`), so it would reach the `eps -> 0` limit.  4 parameters → degree 6 is
already marginal at `eps = 0.01`.  The degree needed grows with the number of independent angles.

**Sign orthants.**  Two squares with free sign of `s` fail at degree 4 for every `eps` and at degree 6
below `eps = 0.1` (`lam = 9e-8` at `eps = 0.01`).  The kinks `|s_i|` in the walls force LP multipliers
`max(s, 0)`, which are not polynomial.  The same two squares in the `++` orthant certify at degree 4 with
`lam = 0.124 eps` down to `eps = 1e-4`.  All 32 orthants of R5 and of P5 certify at `eps = 0.1`, order 2,
`pdeg 1` (`lam` `1.1e-3`–`4.7e-3`).  The smallest margins are the coherent orthants `+++++` and `-----`,
matching the sign-coherence picture of `S6_LOCAL.md` §5.  Scalar product multipliers (`pdeg 0`) fail in all
32: SOS multipliers on the products are needed.

## 4. Measurement 3 — locality  [measured]

Multistart SLSQP minimum of `T` over the type's feasible set (`sos_probe_enc.numeric_minT`, 40–60 starts):

* row of three alone (x-typed): `T = 2.44` (a tilted staircase); every 4-square subsystem of `ROW_A`
  containing the row: `2.56–2.68`.  5-square: `{row, (1,1), (2,2)}` = R5 and `{row, (2,0), (2,2)}` give `3`;
  `{row, (1,1), (2,0)}` gives `2.83`.
* pinwheel: every 3- and 4-square subsystem `< 3`; of the 5-square ones only P5 gives `3`
  (`2.91`, `2.94`, `2.83`, `2.98` for the others).
* **The core needs context, always.**  The row core is false at any small tilt (first order: `T >= 3 - t`
  for a staircase), and the minimal true systems use 5 of the 6 squares.
* **Transfer.**  R5 typed with diagonals in `x` vs in `y` differ on 4 pairs.  The common system (those 4 pairs
  dropped) is false: `T = 2.47` (`2.83` inside `+++++`).  For P5 it is `2.00`.  Dropping a single differing
  pair already gives `2.56`.  So no certificate can serve two nearby types.  Every type and owner variant
  certifies on its own at `eps = 0.1` (R5: `1.3e-3`, `1.8e-3`, `4.0e-4`, `2.0e-3`; P5: `1.5e-3`, `6.5e-4`,
  `2.6e-4`, `1.1e-3`).  Interior-point certificates use every pair row (support: 280–310 of ~305
  blocks), so support gives no locality either.
* Correlative sparsity also loses.  Restricting each row's multiplier to its own squares' angles (`--csp 2/3`)
  kills the `eps = 0.1` certificate (`7e-9`).  The LP dual weight of a wall of square `j` is a function of the
  *other* chain squares' angles (two squares: `c_1 (T - x_2 - p_2)`).

## 5. Measurement 2 — exactness

**Infeasibility form, `eps > 0`: rounds at once [proved, modulo §1].**  The procedure: re-solve at
`lam = lam*/2` maximising a common interior margin `mu`; round to `2^-30`; project exactly with
`A_R^T (A_R A_R^T)^-1 r` (python-flint, `R` an independent row set); check the identity on all rows and
every Gram by exact LDL; then a float check with the unreduced generators at random points of the variety.

| system | `T` | `lam` (exact) | `mu` | rows / rank | max denominator | exact-stage time | result |
|---|---|---|---|---|---|---|---|
| 2 squares `++` | 1.99 | `387/745048` | `1.6e-5` | 125 / 101 | 343 bits | 0.1 s | exact |
| P5 `a,a,a,b,b`, order 2 | 2.999 | `3/150784` | `7.5e-8` | 275 / 215 | 679 bits | 0.2 s | exact |
| R5 `a,a,a,b,b`, order 3 | 2.99 | `158/489981` | `8.5e-7` | 683 / 585 | 2,028 bits | 9 s | exact |
| **R5 full 10 angles, `+++++`, order 2** | **2.9** | `287/433091` | `7.4e-7` | 2,991 / 2,991 | 8,738 bits | 360 s | **exact** |

The last row is an algebraic proof that R5, with that type and every tilt in `[0, 45 deg]`, has no packing
in `[0, 2.9]^2`.  It covers 1 of the 32 orthants of 1 of the types, at `eps = 0.1`.

**Bound form, `eps = 0`: needs facial reduction, which did not close.**  The 1-parameter families give
`gamma = 3.0000000` at degree 4 (statuses mostly `optimal_inaccurate`).  At the zero-margin point every term
vanishes, so the SDP is not strictly feasible.  `sos_probe_facial.py` did one reduction step: Gram matrices
of generators that are slack somewhere on the plateau (40 LP vertices at `theta = 0`, `T = 3`) and `sigma_0`
restricted to `v0-perp`, and scalar products positive on the plateau dropped.  It leaves `mu = -5.9e-12` (P5
uniform) and `-3.4e-10` (R5 uniform): still no interior.  The next level is forced by the order of contact.
`T* - 3` is second order in the tilt while `s_i` and some rows are first order, so terms of order `t` must
vanish too.  That needs the near-optimal curve `x(t)`, i.e. the second-order analysis of `S6_LOCAL.md` §3,
done per family and per stratum.  The strict Stengle form (`T < 3`) has the same degenerate face.  Not
built.

## 6. Extrapolation to `n = 12`  [heuristic]

Every multiplier below is measured at `n = 6` or counted directly.

* **Pieces.**  Each separation type needs its own certificate (§4), and each sign orthant its own (§3):
  `2^12 = 4,096` orthants per type at `n = 12`, against 32 at `n = 6`.  The number of types is the leaf count
  of the disjunction (`S6_LOCAL.md` §4: 2,284 tight leaves at `theta = 0` already at `n = 6`).
* **Core size.**  At `n = 6` the smallest true subsystem has 5 of 6 squares.  Expect 9–12 of the 12 at
  `n = 12` (chain of four plus context); the leaf-A / rank-8 sub-configurations have 8–12.  That is 18–24
  angle variables.
* **Degree.**  Degree 4 fails at `eps ~ 0.02` at `n = 6`.  Degree 6 is linear in `eps` for 2 parameters but
  marginal for 4.  Because the `n = 12` pinwheel deficit is quadratic everywhere (`-(1/3) t^2`, no cubic
  direction), degree may grow less steeply than at `n = 6`, but nothing measured supports degree `<= 6`.
* **Size.**  The full `n = 6` R5 at degree 6 is `M = 18,691` rows and more than 30 GB of KKT.  With 24 angle
  variables at degree 6 there are about `10^5` angle monomials times 25 centre factors, so `M ~ 2 x 10^6`.
  There are ~200 Gram blocks of size ~300 (svec `4.5 x 10^4` each).  That is far outside any interior-point
  solver, and per piece.  Sparsity does not rescue it: correlative splitting kills certificates already at
  `n = 6` (§4).
* **And the result would be `T > 4 - eps`, not `T >= 4`**, unless the facial-reduction step at `eps = 0` is
  solved.  At `n = 12` that is the second-order rigidity theorem on the coherent cone.  The margin there is
  `3.76e-6` on the rank-8 sub-family (`notes/n12-gap.md` §4.7), so `eps` must be below `10^-6`.  At
  `n = 6`, degree 4 already fails at `eps = 0.02`.

## 7. CPU and what failed

CPU: ~20,000 CPU-s total [measured from logs, plus estimates for unlogged runs].  That splits into
~6,700 s of logged single-threaded SDP time, ~6,900 s from one early run that BLAS/faer ran on 12–32 threads
by default (since forced to 1 thread), ~2,000 s in three dense Putinar runs that reached 30–50 GB and were
killed, and ~1,500 s in numeric sanity, exact rounding and failed solver attempts.  No core pinning; peak
concurrency 10 single-threaded jobs.

Failures and fixes: (i) angle range `[0, 90)`: types degenerate at `90 deg` → `[-45, 45]`.  (ii) Putinar with
multipliers in all variables: `1.85 < 2` on two squares and 30–50 GB at 5 squares → angle multipliers.
(iii) No orthant split: `|s|` kinks → split.  (iv) Clarabel `NumericalError` from linearly dependent monomial
rows → drop dependent rows (pivoted QR of `A A^T`), check consistency on all rows after solving (residuals
`1e-11`).  (v) Degree 6 on 10 angles → memory, open.  (vi) `eps = 0` → no interior, one facial-reduction step
insufficient.

## 8. Reproduce

    python3 -m venv runs/sos_venv && runs/sos_venv/bin/pip install cvxpy clarabel scs sympy numpy scipy python-flint
    # degree-4 certificate, full R5, one orthant (25 s)
    search/sos_probe_run.sh --cfg ROW_A --sub 0,1,2,3,5 --orthant=+++++ --prod --pdeg 1 --order 2 --T 2.9 --csp0 3
    # its exact rational version (7 min)
    runs/sos_venv/bin/python search/sos_probe_exact.py --cfg ROW_A --sub 0,1,2,3,5 --orthant=+++++ --prod --pdeg 1 --order 2 --T 2.9 --csp0 3
    # eps sweep on a family
    search/sos_probe_sweep.sh "2 3 4" "0.1 0.03 0.01 0.003 0.001" --cfg ROW_A --sub 0,1,2,3,5 --orthant=+++++ --prod --classes a,a,a,b,b
    # eps = 0 bound + facial reduction attempt
    runs/sos_venv/bin/python search/sos_probe_facial.py --cfg PINWHEEL --sub 0,1,2,3,5 --orthant=+++++ --prod --bound --classes a,a,a,a,a --order 2
    runs/sos_venv/bin/python search/sos_probe_report.py          # tables from runs/sos_probe/res/

Orthant strings beginning with `-` must be passed as `--orthant=-+...`.
