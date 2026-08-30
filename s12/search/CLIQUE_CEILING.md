# The clique ceiling, packing side (task A, 2026-08-29/30)

Code: `search/clique_ceiling.py` (new; imports `packing_dual.py`, `clique_check.py`, `dual_exact.py`;
nothing in the certificate pipeline changed).  Runs: `runs/cc_<TAG>.log|.json|_support.txt|_exact.txt`.

## Question

The value of the clique-strengthened relaxation

    max sum mu_S   s.t.  coverage(p) <= 1 (all points p),   sum_{S in K} mu_S <= 1 (all cliques K),   mu >= 0

over closed unit squares in `[0,t]^2`, as a function of `t`.  By weak duality this is a lower bound on
the weight of **every** clique certificate (points, box cliques, any pairwise-closed-intersecting pose
family) at `t`, so a certified feasible measure of mass `>= 12` at `t` is a ceiling: no clique
certificate below 12 exists at any `t' >= t`.  Measured at `t = 3.99, 3.995, 4.0` (closed semantics)
and for the corner leaf `k = 4` at `3.98` (corner mass exactly 4).

## Semantics: closed intersection

A certificate at `t` refutes packings at side `t' < t`; rescaled to `[0,t]^2` those are pairwise
**disjoint closed** unit squares.  A valid clique is therefore a set of poses that pairwise
**closed-intersect** (touching counts), and a feasible measure must respect all of them.
`clique_check.py`'s separating-axis test was strict; here the SAT test is `<=` (floats with a `1e-9`
tolerance to *choose* cuts, exact integers to *certify*).  Pairwise touching alone is already imposed
by the point constraints (a touching pair shares a point), so the difference lives only in non-Helly
cliques with touching members.

## Method

*Rows.*  Point constraints at the violated arrangement vertices of the current support (the
float-certified enumeration of `packing_dual.py`, `DUAL.md`), plus a coarse seed grid.

*Clique cuts.*  On the D4-expanded support (every image with mass `mu/8`): the exact max-mass clique
(branch and bound), a greedy clique through each of the 300 heaviest images, and the B&B-best
clique through each of the 10 heaviest.  A cut is stored as a union of **pose boxes** around its
member images, sized so that the cores of any two boxes still intersect (so every pose in the union
pairwise meets every other — the same object task B certifies); membership of *new* columns is
geometric, and the reduced cost of a candidate pose includes the duals of the binding cuts it lies in.

*Columns.*  `packing_dual.py`'s pricers (exact verifier witnesses at `N = 1000`, the exact
closed-square sweep, perturbations of the support), reduced cost = point capture + clique duals
(+ `lambda` for corner poses in leaf mode).

*Refinement ladder.*  Explicit cuts against column generation lose: a new pose enters outside every
cut, and the LP re-forms the same clique structure a hair away.  So the loop is staged: **fix the
pose set, run rows + cuts to convergence** (max clique `<= 1 + 1e-5`, no violated vertex) — that
value is a certifiable feasible mass on the current poses — **then price**, add up to 300 improving
columns, and repeat.  The sequence of stage values shows whether the clique-LP value climbs back
towards the pure value as the pose set is refined, and the final pricing gap says how far the last
stage is from optimal over all poses.

*Exact certification* (`--exact`): poses snapped to rational rotations `2 arctan(p/q)` and rational
centres, exactly admissible; coverage at every arrangement vertex in exact integers
(`dual_exact.py`, generalised to any rational `t`); closed intersection of every pair of images in
integers (SAT on integer corner coordinates, `<=`); the maximum clique mass by branch and bound on
**integer** masses; masses rounded down and scaled so that both maxima are `<= 1` exactly.  In leaf
mode the rows carry a `1e-6` margin and the corner mass is adjusted to exactly `k` after rounding.

## A structural finding: the excess is carried by grazing contacts

On the certified pure measure at `3.99` (`PA2`, mass 12.008): the max closed clique has mass
**1.327** (221 images) — the wall point-clique of mass exactly 1.000 plus interior tilted poses of
mass 0.327 that meet every one of them.  Per-member SAT margins (how far a pair may shrink and still
meet): min `8.7e-6`, 10th percentile `5e-3`, median `0.038`.  Pruning the members whose margin is
below `1e-3` — six of them — drops the clique's mass from 1.327 to **1.0096**: essentially the whole
non-Helly excess sits on poses that barely touch the wall squares.  The same is true of the greedy
cliques.  Consequences:

* Box cliques of any real width capture almost none of the measured violation: with boxes the LP
  simply moves the grazing poses just outside the boxes.  A box-clique *certificate* (task B) will
  therefore need very thin boxes along the grazing directions, or a lemma-based family that includes
  the touching limit (the point-clique of the wall point is exact; what is missing is "all poses that
  meet every square through `p`", whose boundary is the grazing set).
* On the packing side the honest computation is the ladder above with explicit (zero-width) cuts on
  a fixed pose set, then refinement.

## Results

(filled in as the runs finish; `runs/cc_C99.log`, `cc_C995.log`, `cc_C40.log`, `cc_L98.log`)

| t | mode | stage values (fixed poses, cuts converged) | last pricing gap | certified exact mass | max clique / M (exact) |
|---|---|---|---|---|---|
| 3.99 | pure | | | | |
| 3.995 | pure | | | | |
| 4.00 | pure, closed | | | | |
| 3.98 | leaf `k = 4`, `r = 1` | | | | |

Reference (pure point-LP, certified): `L(3.99) = 12.0082`, `L(4.00) = 12.163`; corner leaf `k = 4`
at 3.98 `= 12.000 +- 0.001` (`DUAL_EXACT.md`, `DUAL.md`, `CLIQUE.md`).

## Reproduce

```sh
R=runs   # or the main tree's runs/ for the warm starts
python3 search/clique_ceiling.py 3.99  C99  --warm $R/dual_PA2_support.txt --seed-pitch 0 --stages 10 --inner 150 --N 1000 --threads 1 --dmin 0 --cg-want 300 --exact
python3 search/clique_ceiling.py 3.995 C995 --warm $R/dual_PA2_support.txt --warm $R/dual_PC1_support.txt  (same flags)
python3 search/clique_ceiling.py 4     C40  --warm $R/dual_PC1_support.txt --warm $R/dual_PA2_support.txt  (same flags)
python3 search/clique_ceiling.py 3.98  L98  --warm $R/dual_PD1_support.txt --branch-dual $R/branch_t398hk4_dual_it16.txt --kmass 4 --r 1 --margin 1e-6  (same flags)
python3 search/clique_ceiling.py 3.99 x --check runs/cc_C99_exact.txt      # independent exact re-check
```
