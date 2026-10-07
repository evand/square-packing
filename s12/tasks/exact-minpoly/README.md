# Exact forms for every best-known packing, n ≤ 324 (opened 2026-10-06)

**Goal.**  For every n ≤ 324, the best known packing's side length is specified exactly, by a mechanized and
re-checkable process, as *the unique root of an explicit integer polynomial in an explicit rational interval*.  It should
also come with an exact configuration that realizes it.  This is the "known exactly" half of the local-optimality
statement (the other half is the local-optimality checks).  It should be citable from the eventual s(11) explainer.

## Deliverable per n

1. `p ∈ ℤ[x]`, irreducible and primitive, plus rationals `a < b` such that `p` has exactly one real root in `[a, b]`.
   `S* :=` that root.
2. An exact configuration over a number field `K = ℚ[t]/f(t)`: `f` irreducible, with a real root `t*` in an isolating
   interval.  Each square's `(x_i, y_i, c_i, s_i)` and `S` are elements of `K`, with `c_i² + s_i² ≡ 1` (mod f), and
   `p(S(t)) ≡ 0` (mod f).  The contact list (which corner touches which side or wall) is part of the data.
3. Checks that don't depend on how the data was produced (msolve, LLL, anything else):
   * every contact equation is ≡ 0 (mod f): exact polynomial arithmetic;
   * every other pair of squares is separated, and every square lies inside the container: interval arithmetic at `t*`;
   * `S(t*)` lies in `[a, b]`.

   Together these give `s(n) ≤ S*` exactly, which is sharper than today's rational certificate (`S' = S*(1 + 1e-20)`).
4. Cross-checks: `S*` agrees with `S_exact` in `search/exact/batch/results.json` (80 digits), and the degree agrees with
   independent sources where they exist (below).

The polynomial for S depends only on the contact structure.  The canonical placement (below) affects only the field in
which the coordinates live, and so the size of the proof object.

## State on 10-07 (pipeline built; read this first)

Tools (all in `search/exact/`, run via `./env.sh`, which pins numpy 2.4.2 / scipy 1.17.0 / python-flint 0.9.0,
msolve from nixpkgs):
* `exactsolve.py` writes `NAME.contacts.json` (closed contacts, load-bearing set, corner-corner touches, free squares).
  Precision fix: residual scaling only below 1e-200 (applying it at 80 digits flipped a rank decision at n = 202).
* `minpoly.py`: number field, exact configuration, minimal polynomial of S.  Angle classes (parallel or mirror mod
  90°, t = tan(θ/2)); elimination of centres over Q(t₁..t_k); k = 1: gcd of consistency numerators; k ≥ 2: their
  irreducible factors through t*, then msolve's parametrization; force balance when the contacts leave one direction
  (dS/dt = 0; Lagrange determinant for a curve of angles); flat directions settled on near-contacts, then pinned at
  30-digit rationals.  Refuses invalid exactsolve points; drops listed equations that do not hold at the point.
* `verify_exact.py`: independent stdlib checker (exact identities mod f, rational intervals; lazy exact values).
* `minpoly/solve_all.sh`, `minpoly/run_minpoly.sh` (systemd MemoryMax caps; the machine is shared), `minpoly/summary.py`
  → `minpoly/results.md`.
* `localmin.py`: grade-A local-minimum certificates (below); `lean_cert.py`, `localmin_lean.py`: Lean certificates.

Results (10-07): **258+ of 323 register records verified exact** (fixes since then cleared 87, 102, 127, 150, 202;
final table in `minpoly/results.md`); 81 cross-checks (findpoly, Ellsworth degrees) agree, none differ.  Open: ~25
many-class records time out (needs block-triangular elimination), ~12 with more than one force-balance direction, ~10
flat rotations parallel to fixed squares, msolve-hard (29, 68, 71, 126, 228, 55, 182, 83, 235), 108/206 (large field).

Lean (`lean/Sqpack/`, no `native_decide` anywhere):
* `ExactPack`, `ExactCheck`: exact data ⇒ `Packs`; kernel-checked certificates `Exact/N5c`, `Exact/N11c`
  (n = 17 compiles: 402 s, 15.6 GB, not committed).  Frontier to aim for: field degree ≤ 8 (249 records incl. the
  176 axis-parallel ones via `packs_grid`); needs a cheap far-pair test and integer arithmetic, not compute.
* `LocalMin` (`IsLocalMinPacking`: a pose-space ball around the record with no packing in a smaller square; rigidity
  lemma; row expansion with w = sin δ, u = cos δ − 1, |u| ≤ w², no Taylor), `LocalMinRows` (`isLocalMin_of_rows`),
  `LocalMinCheck` (certificates over Q[T]/f, soundness `isLocalMin_of_lcert`).  Instance: `Exact/N11L` (n = 11).

### Local minima: grades (exactsolve's checks over the register)
* **A** (20): first-order rigid with all forces > 0: n = 1, 11, 28 and the perfect squares (those are globally optimal
  by area anyway).  Proof: rows "point P of closed square j lies beyond side k of i" (valid near the record when P
  projects strictly inside the side: P ∉ open square i, the other three side values stay < 1/2) and walls; aligned
  side-side pairs become one *midpoint* row (the average of the two corner incidences; branch-free, no SAT needed,
  and n = 11, 28 stay rigid with the same max-min λ).  Certificate: exact λ ∈ K (Σλ L = e_S, λ > 0), rational G with
  ‖GL − I‖ ≤ 1/2, margins.  `localmin.py`: n = 11 (λ_min 0.048, ‖G‖ 18.9, margin 0.021) and n = 28 valid.
* **B** (1): n = 5, strict only at second order.
* **C** (295): strict modulo flat motions.  Next: the paper argument.  Removing squares only relaxes the
  constraints, so it suffices to treat the load-bearing squares; the open part is load-bearing squares in exact flat
  families (often translations).
* Out of scope for now: 177, 211, 230, 261, 263 (no strictly positive forces), 272 (not a local min), 292.
* The definition fails on broad plateaus with a far downhill exit (s(7)'s L arrangements): accepted for step one.

## State on 10-06 (exploration in the main session; prototypes in `proto/`, data in `results/`)

The 323 certified register packings (`search/exact/batch/`):

| class | count | status |
|---|---|---|
| axis-parallel only | 176 | S is an integer: `x − k` |
| tilted, found by `mpmath.findpoly` from the 68-digit `S_exact` (degree ≤ 12, coefficients ≤ 10⁶) | 78 | mostly quadratics in √2 or √7; rationals 643/41 (n = 230), 684/41 (n = 261), 53/7 (n = 50); degree 4–8 at n = 11, 28, 39, 54, 70, 107, 153, 178, 267, 301 (`results/findpoly68_tilted.txt`) |
| tilted, not found | 69 | 2 now done (n = 17, 37); the rest open |

Plus n = 292 (no exact KKT point) and our unposted records 266, 270, 272 (`search/exact/results/sw2/`).

Two routes worked on 10-06 (the prototypes hard-code 10-06 scratch paths and a scratch venv: fix before reuse):

* **LLL at high precision** (`proto/lll2.py`, `proto/minp.py`, python-flint): n = 17 at 2000 digits gives an
  irreducible degree-18 polynomial with 12-digit coefficients, residual 1e-1980.  This needs the exactsolve precision
  fix (`proto/exactsolve-precision.patch`): stock `exactsolve.py` converts residuals to f64, so they underflow below
  1e-308, and the projection step is capped at 40 iterations.  In practice it tops out at about 320 digits; the
  README's "60–1000 digits" is wrong.
* **Symbolic elimination from the contact graph** (`proto/build2.py` → `proto/reduce.py` → msolve → `proto/param.py`;
  driver `proto/run1.sh`).  It runs in four steps:
  * one equation per incidence, at tolerance 1e-200 on the 2000-digit point or 1e-50 on the 68-digit files;
  * **angle classes**: squares parallel mod 90° share one (c, s) and axis-parallel squares get (1, 0);
  * variables with constant linear coefficients are substituted out;
  * msolve computes a rational parametrization over ℚ, then a resultant and factorization give the polynomial for S.

  | n | angle classes | system | msolve | result |
  |---|---|---|---|---|
  | 17 | 2 | 19 equations / 19 unknowns, 20 solutions | 0.2 s | degree 18, identical to LLL's, plus a spurious double linear factor |
  | 37 | 2 | 33 / 27 | ~10 min (2 threads) | degree 8, `36x⁸ − 2496x⁷ + 59768x⁶ − 733760x⁵ + 5289248x⁴ − 23462672x³ + 63458276x² − 96673872x + 64068561` |
  | 29 | 5 | 38 / 31 | > 10 min, killed | — |

msolve comes from nixpkgs: `nix --extra-experimental-features 'nix-command flakes' shell nixpkgs#msolve -c msolve …`.
python-flint was installed in a scratch venv.

**Independent cross-checks available.**
* The register's frontier files (`_untrusted-third-party/jlevy-squares/packing/frontier/n-NNN.md`) have
  `algebraic_degree` and `minimal_polynomial` for 21 n: 11:8, 17:18, 28:6, 37:8, 39:5, 41:42, 51:12, 70:4, 83:24,
  87:44, 88:20, 108:144, 128:40, 129:20, 146:16, 153:4, 179:158, 205:40, 235:83, 266:32, 300:40.  These may describe
  an *older* record: n = 83's degree 24 is Cantrell's 2024 packing.  17 and 37 match ours.
* **Ellsworth's page** (https://kingbird.myphotos.cc/packing/squares_in_squares.html; checked 10-06): for every
  record of degree ≥ 3 it has the polynomial in a hidden `<span class="frames">`, some with coefficients in ℤ[√2].
  The degree shown there is over ℚ.  n = 83's polynomial is a Mathematica `Root[…, 27]` in the comment of
  `square-83.svg`.  Raw downloads are in `private/downloads/kingbird-2026-10-06/` (not republished); our parse and
  verdict are in `results/kingbird_crosscheck.json` (`proto/parse2.py`, `proto/cmp.py`).  The test: the page
  polynomial, normed to ℚ and factored, changes sign within 1e-60 of our 68-digit S.
  * **Agree (22)**, n: degree over ℚ: 11:8, 17:18, 28:6, 37:8, 39:5, 41:42, 51:12, 69:38, 70:4, **83:672** (724-digit
    coefficients), 87:41, 88:20, 108:144, 128:40, 129:20, 146:16, 153:4, 179:158, 205:40, 235:83, 266:32 (the register's
    266, not our new one), 300:40.  This confirms our findpoly results for 11, 28, 39, 70, 153, and our msolve and LLL
    results for 17 and 37.  It settles 17 of our 69 unidentified packings by citation: 17, 37, 41, 51, 69, 83, 87, 88,
    108, 128, 129, 146, 179, 205, 235, 266, 300.  The other **52** have no published polynomial that we know of.
  * **Differ (12)**: Couzo's 2026-09/10 packings are newer than the page, so its polynomial describes the previous
    packing: 102, 106, 123, 130, 172, 177, 199, 206, 228, 259, 269, 302.  292 has no exact point of ours.
  * n = 83's equations (SVG comment) are 3 angles a, b, c plus a 4th equation `Det[Grad[{s, f1, f2, f3}, {s, a, b, c}]]`:
    S is fixed by force balance, not by the contacts alone.  So the Lagrange/determinant route (roadmap 4) is needed
    for the n = 83 test, not optional.
* The register's frontier `minimal_polynomial` fields (above) look like older snapshots of Ellsworth's page (they give
  83:24 and 87:44, where the page now has 672 and 41).  They add nothing independent.

## Canonical placement for non-rigid packings

Rattlers and squares in flat (sliding) families don't change S, but they need exact positions.  Convention (Evan,
10-06, the usual display convention):
* make as many squares as possible axis-parallel or parallel to an existing group (angle class);
* make as many side-to-side contacts as possible;
* a square that slides in 1-D between a side-to-side contact at one end and a corner contact at the other goes to the
  side-to-side end.

Rattlers: rotate to a neighbour's or the grid's angle if they fit, then push to two side-to-side contacts.  The order
and tie-break still need fixing so the result is deterministic.

Implementation sketch: after exactsolve, a "settle" pass.
1. For each flat mode or rattler, find the endpoints of its motion with an LP or a 1-D search.
2. Pick an endpoint by the rule above.
3. Add the new contacts and re-solve.
4. Check that S is unchanged and that the contact Jacobian has full rank.

Settling adds equations and usually lowers the field degree.  Fallbacks if this is hard: force balance (max-min λ,
which is Ellsworth's analytic convention), or pinning a coordinate at a nearby rational.  Pinning is valid but not
canonical, and the prototype's angle-class substitution already does something like it.

Separate issue: packings where the contacts alone don't fix S (S is set by force balance, the under-determined case in
the exactsolve README).  These need the Lagrange equations, which are linear in λ, and they become Ellsworth's
determinant conditions after eliminating λ.  We haven't counted these yet (compare rank J_A with the number of
unknowns in the batch JSON).  n = 17 and 37 were contact-determined.

## Lean demo: one small n, the whole process by hand before automating (Evan, 10-06)

Two levels:
* (a) an exact configuration, its contact graph and the polynomial give `Packs n S*`, with S* exact (a closed form or
  "the unique root of p in [a, b]");
* (b) a full local minimum.

What we know that bears on the choice (`search/exact/batch/results.json`, `second` field):

* **n = 5**, `S = 2 + 1/√2`: everything lives in ℚ(√2), so signs are decidable with rational arithmetic.  No rattlers.
  A *strict* local minimum, but only at second order (reduced Hessian positive definite, not first-order rigid), so (b)
  needs a second-order sufficiency argument.  Easiest arithmetic.
* **n = 11** (octic): the only non-square n that is **first-order rigid** (null J_A = 0) with strictly positive
  multipliers (max-min λ = 3.6e-2).  So (b) is a linear-growth argument: with Σλ_k ∇g_k = ∇S and J_A injective,
  S − S* ≥ c·|Δ| near the point, with explicit Taylor constants.  This is the cleanest "jamming certificate ⇒ local
  min" demo, and topical given the s(11) announcement (where global optimality makes local moot, but it shows the
  mechanism).  The cost is arithmetic in a degree-8 field: signs come from isolating intervals.
* **n = 4 or 9**: rigid axis grids.  A trivial dry run of the plumbing only.
* Every other n has rattlers or flat modes: "strict modulo flat motions" needs more theory, so it comes after the demo.

Pieces either way:
* a non-strict separating-axis soundness lemma (touching closed squares have disjoint interiors);
* exact field arithmetic by kernel `decide`;
* bridging `Spec.unitSq c θ` (it uses `Real.cos θ`) to algebraic (c, s), via θ = arg(c + is);
* for (b): a definition of local optimality in `Spec.lean`, and the **relaxation lemma** (near the configuration, a
  valid packing satisfies g_k ≥ 0 for every active incidence, since the corner of j stays outside square i).  Only
  that direction is needed, and it turns a non-smooth problem into a smooth one.

Existing Lean: `lean/Sqpack/Spec.lean` (`Packs`, `minSide`, `unitSq`), with grid upper bounds only (`packs_grid`).
Nothing about local optimality yet.

## Roadmap

0. **exactsolve precision fix** into `search/exact/` (scale the residual before the f64 solve; scale the polish
   iteration cap with dps; fix the README claim), with a 2000-digit n = 17 regression test.
1. **`search/exact/minpoly.py`**: build the contact graph from exactsolve's own contact list instead of re-detecting
   it; angle classes; linear substitution; msolve; resultant and factorization; choose the factor by evaluating at
   high precision, and require the choice to be unambiguous; isolating interval with flint.  Output `minpoly/n-N.json`.
2. **Elimination that scales with the number of angle classes.**  Given the angles, the positions and S satisfy a
   linear system M(c, s)·z = 0.  Consistency means its minors vanish, which leaves k equations in the k angle variables
   (use t = tan(θ/2) to keep things polynomial); then S comes from Cramer's rule.  Alternatively, take the
   block-triangular (Dulmage–Mendelsohn) decomposition of the Jacobian and solve rigid clusters in dependency order.
   n = 29 (5 classes) is the first benchmark.
3. **n = 83** (3 angle classes, degree 672, 724-digit coefficients; S fixed by force balance): compare with Ellsworth
   and Chang's polynomial (`private/downloads/kingbird-2026-10-06/n83_kingbird_poly.json`).  LLL is hopeless here (it
   would need about 5·10⁵ digits); this is the test of the symbolic route at scale.  It needs step 4's Lagrange
   equations first.
4. **Settling pass** (previous section), plus the Lagrange equations for force-determined S.
5. **`verify_exact.py`**: the independent checker for deliverable item 3 (arithmetic mod f, interval gaps).  It must
   not trust msolve: msolve over ℚ uses multimodular reconstruction.  A second implementation can come later.
6. **Batch** all n into a table: degree, height, field degree, time, cross-checks (LLL where degree × height digits is
   below ~5000; register and Ellsworth degrees).
7. **Lean**: kernel `decide` on the mod-f identities for small fields gives `Packs n S*` with S* exact.  This feeds
   the local-optimality line and the explainer.
8. **Leftovers**: n = 292 (no exact point yet); the bound-only cases 177, 211, 230, 261, 263 (the polynomial is still
   defined for their contact structure, but 211 and 263 are degenerate); the new 266, 270 and 272.

## Pitfalls seen on 10-06

* Gröbner on the raw system, without angle classes: msolve blew up at degree 5 (a 171k × 211k matrix).
* msolve needs integer coefficients.  An input with `c12/2` came back as "no solution".
* msolve picks a separating linear form when S takes the same value on several solutions, so S comes out as a
  rational function of another variable, and a resultant is needed.
* `findpoly`'s `maxcoeff = 10⁶` missed n = 37 (coefficients ~10⁸).
* The 68-digit files support contact detection at tolerance 1e-50, but not identification beyond degree ~12.
