# s(110) landscape explainer (10-08/09)

Data and figures for `site/www/landscape110.html`.  Inputs: the certified s(110) census (exactsolve reports in
`runs/*/exact/`), the corner-movement matrix between minima.

* `morph.py`: square matching (`match2`: assignment by squared corner movement), straight-line frames, and
  `walk_carrot(s, A, B, perm)`: spread A to fill a box of side s, follow the straight-line morph to B with feasible LP
  steps (min weighted L1 distance to the moving target s.t. linearised non-overlap, side fixed; corner-corner pairs on
  the line the step violates least), then shrink to B.  Not symmetric in A and B.
* `build.py reps | barrier | barrier-rev | cycle`: k-medoids representatives (both records forced in; `--keep` grows
  an existing set); barrier = smallest s (bisection, 2e-4) at which `walk_carrot` connects a pair, an upper bound on
  the expansion for this path family (`--near k`: only each rep's k nearest reps); failed pairs retried in the other
  direction; tour = nearest neighbour + 2-opt on sum + 3 max of the excess; leg frames at the barrier side.
* `forces.py`: contact multipliers (f64) at the exact point for the force figure.
* `settle.py`: slp2 squeeze trajectories from the record loosened 3 % and nudged (sigma 0.01 / 0.03 / 0.1).
* `export.py`: one JSON for the page (`site/www/data/landscape110.json`).

Outputs in `runs/landscape110/` (gitignored).  10-09 run: 48 reps, 507 pairs computed (10 nearest each, plus all
pairs of the first 24), ~35 min on 10-14 procs; every tour leg connects.
