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

## Summary

* **No ceiling was found.**  At `t = 3.99, 3.995, 4.0` and in the corner leaf `k = 4` at `3.98`,
  every clique-feasible measure produced has mass **well below 12**; the pure point-LP measures
  (12.008, 12.16) are not clique-feasible (max clique 1.33–1.44) and the LP loses 0.6–0.8 of mass
  when made so.  So this run does **not** exclude the clique method at any `t <= 4`.
* **What is certified (exact):** a clique-feasible measure of mass `125/11 − 6e-9 = 11.3636` at
  `t = 3.99` (and 11.19 from an earlier stage); certified scaled measures of mass 9.8–10.1 at 3.995
  and 4.0 and 8.9 in the leaf (uninformative — their pre-scaling measures were not clique-converged).
  These are lower bounds on the clique-LP value; the interesting direction (upper bounds) needs the
  cover side (task B).
* **What is measured (heuristic, restricted pose sets):** the fixed-pose clique-LP value is
  ≈ **11.2–11.45 at 3.99** (converged at 314–912 poses; 11.37–11.45 with residual max clique 1.02–1.04
  at 1029–1400 poses), ≈ **11.36–11.48 at 3.995** (residual max clique 1.03–1.08), ≈ **11.8–11.85 at
  4.0** (residual 1.13–1.14; not converged after 400 iterations), and ≈ **11.65–11.8 in the leaf**
  (residual 1.12–1.24).  The values drift *up* slowly with pose refinement (3.99: 11.19 → 11.36 →
  11.45 over ~1400 poses), so the true clique-LP value is higher than these; how much higher is the
  open question, and the column-generation pricing cannot tell (the LP is degenerate: the minimum
  reduced cost is 0 at every stage while the value barely moves).
* **Answer to "do cliques alone reach 4?"** — *not decided, but the evidence points the right way*:
  on every pose set tried, including the pure method's own extremal poses plus five rounds of
  pricing, the clique relaxation sits 0.15–0.8 below 12 at `t = 3.99–4.0`, and the leaf that blocks
  `s(12) >= 3.98` drops from 12.000 to ≲ 11.8.  Nothing here is a bound in the useful direction.

## Results

Stage value = LP mass on the fixed pose set when the inner loop stopped; "conv" = no violated
vertex and max clique `<= 1 + 1e-5` (a feasible measure), otherwise the residual max clique is
given (more cuts would lower the value further).  One core per run; the first attempts
(`cc_C99`, `cc_C40`, `cc_C995`, `cc_L98`) ran with a 150-iteration inner cap and lost to a solver
failure or a bookkeeping bug at a stage boundary; the `b`/`c`/`d` runs restart from the previous
run's last support with the fixes and a 400-iteration cap.

| t | mode | poses → stage value ("mc" = residual max clique; each run restarts from the previous run's support, so pose counts restart) | certified exact mass | notes |
|---|---|---|---|---|
| 3.99 | pure | `C99`: 172 → 11.205 (mc 1.03), 468 → 11.211 (mc 1.07); `C99b`: 476 → **11.191 conv**; `C99c`: 314 → 11.191 conv, 613 → 11.191 conv, 912 → **11.364 conv**, 1029 → 11.374 (mc 1.045), 1170 → 11.374 (mc 1.018); `C99d`: 427 → 11.365 (mc 1.012), 622 → 11.448 (mc 1.044) | **11.363636358** (`runs/cc_C99c_exact.txt`, copied to `search/clique_ceiling_3.99_exact.txt`; re-checked by `--check`: 11 poses, 88 images, M = 1, max clique = 1, exact) and 11.191489 (`cc_C99b_exact.txt`) | pure `L(3.99) = 12.008` |
| 3.995 | pure | `C995`: 245 → 11.58 (mc 1.13); `C995b`: 432 → 11.366 (mc 1.034), 728 → 11.481 (mc 1.075); `C995c`: 495 → 11.357 (mc 1.047), 748 → 11.377 (mc 1.23, cut short), 908 → 11.88 (mc 1.19, cut short) | 9.91, 10.10 (scaled, uninformative) | pure `L(3.995)` in `[12.008, 12.163]` |
| 4.00 | pure, closed | `C40` (150-cap): 384 → 11.883 (mc 1.16), 542 → 11.950 (mc 1.20), 650 → 12.042 (mc 1.19), 726 → 12.132 (mc 1.20); `C40d` (400-cap): 405 → **11.795 (mc 1.13)**, 513 → 12.034 (mc 1.20); `C40e` (400-cap, 3× cuts/iteration, still running at it 319): 207 → 11.849 (mc 1.14), descending 12.05 → 11.96 → 11.87 → 11.85 at it 100/200/300 | 9.82 (scaled, uninformative) | pure `L(4.00) = 12.163`; the t = 4 inner loop does not converge in 400 iterations |
| 3.98 | leaf `k = 4`, `r = 1` | `L98`: 1206 → 11.744 (mc 1.18); `L98b`: 1302 → 11.651 (mc 1.12), 1449 → 11.766 (mc 1.23); `L98c`: 1245 → 11.681 (mc 1.12), 1377 → 11.815 (mc 1.24, cut short) | 8.90 (scaled, corner mass 3.02 after scaling — not a leaf statement) | leaf value without cliques `12.000 ± 0.001`; HiGHS needs the dual-simplex fallback on most leaf LPs |

Wall clock: 6 h per run per restart, one core each (the inner loop is 5–120 s per iteration, dominated
by clique separation once there are ~30k cuts); ~40 core-hours in total.

**The certified measure at 3.99** (mass 125/11): corners **4.000** (the pose `(0.5, 0.5, 0°)`, one
full square per corner), the eight wall slots **4.000** (3.273 on `(0.5, 1.5, 0°)`, 0.727 on
`(0.56, 2.08, 7.5°)`), and **3.364** in the interior on seven tilted poses (0.09–0.55 each, 8.4°–37°,
centred 1.2–1.7 from the walls).  Against the pure measure at the same `t` (corners 3.40, walls
~4.3, interior ~4.3) the clique constraints have made the frame integral and cut the interior by
about one unit; the residual excess over 11 is 0.36 and it is entirely interior.

Reference (pure point-LP, certified): `L(3.99) = 12.0082`, `L(4.00) = 12.163`; corner leaf `k = 4`
at 3.98 `= 12.000 +- 0.001` (`DUAL_EXACT.md`, `DUAL.md`, `CLIQUE.md`).

## Reading, and what to do differently

1. **The ceiling question is not answerable from the packing side alone with this loop.**  Every
   stage value is a restriction (fixed poses), hence a lower bound on the clique-LP value, and the
   values drift upward with refinement.  A ceiling needs a clique-feasible measure of mass `>= 12`,
   and none of ~40 core-hours produced one anywhere in `[3.99, 4.0]`; that is evidence that none
   exists, not a proof.  The decisive object is now on the cover side: a box-clique *certificate*
   (task B) of weight `< 12` at, say, 3.98 in the `k = 4` leaf, then at 3.99 pure.  The measures
   here say where its cliques must be (walls + interior straddlers) and how thin they must be
   (grazing contacts, §"structural finding").
2. **Where the mass goes when cliques bind.**  Certified 3.99 measure: corners 4.000 (integral),
   wall slots 4.000 (integral up to one 7.5° tilt), interior 3.364 on seven tilted poses.  The clique
   constraints make the frame integral; the residual fractional excess over 11 is interior and
   ≈ 0.36 on this pose set.  In the leaf dual the same happens with corner mass pinned at 4.
3. **Degeneracy is the practical obstacle**, in two forms.  (a) The LP has a large optimal face:
   thousands of cuts are needed before the value moves (mass migrates between near-identical poses),
   and at `t = 4` the 400-iteration inner loop still leaves a residual max clique of 1.13–1.2.
   (b) Column pricing is uninformative: the minimum reduced cost is 0 at every stage, yet adding
   300 such columns moves the value by 0.0–0.2.  A proper answer wants either the cover side
   (where cliques are columns and the verifier is the oracle) or a stabilised column-and-cut
   generation (dual smoothing on both point and clique duals, in-out separation).
4. **Costs.**  Exact certification (rational SAT + integer-weight B&B + arrangement vertices) is
   minutes; the search is hours.  The `--check` mode re-certifies an exact support file independently.

## Reproduce

```sh
R=runs   # or the main tree's runs/ for the warm starts
python3 search/clique_ceiling.py 3.99  C99  --warm $R/dual_PA2_support.txt --seed-pitch 0 --stages 10 --inner 150 --N 1000 --threads 1 --dmin 0 --cg-want 300 --exact
python3 search/clique_ceiling.py 3.995 C995 --warm $R/dual_PA2_support.txt --warm $R/dual_PC1_support.txt  (same flags)
python3 search/clique_ceiling.py 4     C40  --warm $R/dual_PC1_support.txt --warm $R/dual_PA2_support.txt  (same flags)
python3 search/clique_ceiling.py 3.98  L98  --warm $R/dual_PD1_support.txt --branch-dual $R/branch_t398hk4_dual_it16.txt --kmass 4 --r 1 --margin 1e-6  (same flags)
python3 search/clique_ceiling.py 3.99 x --check runs/cc_C99_exact.txt      # independent exact re-check
```
