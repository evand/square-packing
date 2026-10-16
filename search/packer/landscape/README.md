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

Connectivity map (10-11, exploratory): `mapdata.py` (system python, from packer/: certified minima from the store's
`minimum` table, weak / unresolved left out; junction-aware families; all-pairs fingerprint distance, ~15 s) then
`./landscape/mapenv.sh python3 landscape/maps.py` (nix python with umap-learn / scikit-learn): UMAP, t-SNE, MDS layouts,
family graph, disconnectivity graph of the tour minima from the stored barriers.  Output `runs/landscape110/map/`.
Hover-map tuning (10-11): `maps2.py` (two-level disc layout: rejected, loses within-family structure), `sweep.py`
(t-SNE / UMAP on fingerprint distance ** p; crowd / trustworthiness / global-Spearman scores in `map/sweep.json`),
`render.py NAME ...` (high-res render of a sweep layout, labels per family island).  Front-runner:
UMAP p = 0.3, k = 250, min_dist 0.5, spread 1 (crowd 9.5 %, trust 0.986, global 0.61).  t-SNE stays ~80-90 % crowded.
Fingerprint checks vs exact matching distance (10-11, `fpres.py`, `fpangle.py`, on runs/count_sep/k150.json): grid size
is irrelevant at fixed blur, sharper blur is worse (blur 0.275 box units = current is best of those tried); angle channels
(a*cos 4θ, a*sin 4θ, sin negated under reflections) add little (best: current + angle a = 0.5, Spearman 0.938 vs 0.928,
near pairs 0.756 vs 0.737, 10-NN overlap 85 % both); corners image: better nearest neighbour (83 % vs 77 %), worse globally.
Fill vs separation (10-11, `sweep_sep.py`; `sweep.py --fill`): fill = occupied share of a 12 x 12 square frame; sep =
median nearest other-family distance / median nearest same-family distance; 5 seeds each.  min_dist / spread sets the
trade (smaller min_dist: tighter islands, more sep, less fill, more crowding); p barely moves sep; p = 0.85 fills more
than 0.7 (41 vs 37 %).  Candidates (p = 0.5, k = 250): md 2 / spread 4 seed 4 (fill 58 %, sep ~3.3, Evan: families a bit
close), md 1 / spread 2 seed 3 (51 %, 3.47), md 1.5 / spread 4 seed 2 (44 %, 3.63).
Interactive mock-up (10-11): `mock.py NAME ...` + `mock_template.html` -> `map/mock.html` (data embedded, opens from file://;
layout switcher, hover preview, click to hold, family ring, island labels).
Count families vs islands (10-11, `splits.py`): most big count families split into 2 fingerprint clusters ~0.05 apart
(= the typical between-family distance) while each half's nearest other count family is 0.002-0.005 away: counts change
by cheap corner moves (one square grid <-> band), islands are the structure.  Clusters at fp distance 0.02: 36 (21 with
>= 5; same at 0.015; 0.03 chains into one of 629).  Record island: 91 minima, 6 count families, all with 67 axis squares.
In about half the split families the two islands place the axis squares differently (corner block vs strips, +-1).
Decision (Evan, 10-11): map layout = UMAP p = 0.5, k = 250, min_dist 1 / spread 2 or min_dist 1.5 / spread 4 (either);
min_dist 2 / spread 4 rejected (too little separation between groups).  Seeds matter (store the seed with the layout).
Two deep pockets (10-11): islands at fp cut 0.02 sorted by best side: record island 10.9967833966 (91 minima, all 67
axis) and 10.9968869868 (88, all 68 axis; 19/19/68/4 x 80) hold 83 of the 88 minima below 10.998; next island
10.9976461385.  Fp distance between them min 0.026 / median 0.040 (median between islands 0.070).  Smallest computed
tour expansion between them ~11.023 (upper bound, few pairs).  Tour has 10 reps in the record island, 3 in the second.
Baseline (Evan, 10-11): `umap_p=0.5_k=250_md=1_spread=2_seed=3` (in map/sweep.npz).
Orientation (10-11).  Store keeps arrival orientation (23 % of census minima are the diagonal transpose of the record's).
Display: `mapdata.py` orients by tree propagation (record as is; along the fp-distance MST each minimum takes the box
symmetry closest to its oriented parent): neighbour pairs (5-NN) drawn inconsistently 23.2 % -> 0.2 % (count rule
"larger band left": 1.4 %); larger band changes wall in ~3 % of neighbour pairs under either rule (real count changes).
Tour: `morph.match2` ignores box symmetries.  `orient_check.py` (5 largest-expansion legs where the transpose lowers
match2's cost): re-oriented expansions worse in 4 / 5 (fail, 0.167 vs 0.086, 0.119 vs 0.045, 0.060 vs 0.021), better in
1 (0.012 vs 0.028).  Matching cost does not predict the walk; the stored legs are not inflated by orientation in general.
New tour legs: run walk_carrot for both orientations of B (identity, transpose) and keep the smaller.
walk_carrot is not transpose-symmetric (10-11): pair 10.998021150 / 10.998391887 connects at s = 11.2312 as stored and fails
up to +0.005 transposed (same match2 permutation).  Not rounding (1e-12..1e-9 jitter: still connects), not square order
(3 random renumberings: same split), not LP column order (x/y columns swapped: same split); source not found (candidates:
inc.model_fast, _correct, slp2.worst).  Use: walk a pair in either frame (`build.py --frame 0|1`; frame 1 = transposed,
path mapped back by `tframe`), keep the smaller s; frames are always shown in the display orientation.
Tour v2 (`runs/landscape110/tour2/`, `tour2_prep.py`): 56 reps (47 v1 still certified + 9: uncovered islands >= 5, second
pocket to 6), display orientation; 334 v1 pairs seeded (302 frame 0, 32 frame 1).  Pass 1: frame 0, near-10 by fp
distance (196 new pairs, 4 procs).  Pass 2 (next): frame 1 for the tour legs and the largest expansions, then cycle.
Tour v2 passes (10-11): pass 1 (frame 0, near-10, 196 pairs, 50 min on 4 procs) -> 530 records; pass 2 (frame 1: tour
legs, the two unconnected stops, pairs touching the bad legs; 158 pairs, 35 min) -> 688.  Frames on 170 pairs walked both
ways: frame 0 better 54, frame 1 better 44, tie 72, median difference 0; frame 1 rescued 1 of 34 frame-0 failures.
Dropped (Evan 10-11: show, not complete; rule = no two connections <= 0.05 unless a pocket bottom): 10.9982332904,
10.9990062819, 10.9993587202, 10.9994555129, 10.9998468069 (only 40-tilted stop), 10.9998083737 (leg 0.079, small island).
Kept 10.9980211505 (only stop of the third-deepest island, 41 tilted) despite a 0.088 leg.  `build.py cycle --drop`.
Tour v2 result (10-11, `runs/landscape110/tour2/cycle.json`): 50 stops, 50 / 50 legs connect (14 walked in frame 1 and
mapped back), expansions <= 0.065 (three legs 0.053-0.065, the rest <= 0.04); consecutive legs join exactly after square
matching (position 0, angle 2e-16), display orientation throughout.  Compute ~7 proc-h (wall x 4 procs; not cpu_s).
Not yet on the page: export.py still reads runs/landscape110/{reps,cycle}.json (v1).
Page v2 (10-11): `site/www/landscape110.html` + `data/landscape110.json` (export.py v2: store snapshot after census v3,
775 certified minima, store sha 364894392d35; layout.py seed 0 by the fill rule; tour v2; forces for 43 v1 tour stops).
Sections: census rows by tilted count (gold ramp) -> map (side teal ramp / tilted toggle) sharing one packing panel ->
tour (square-root side scale above the record) -> settling -> why stuck.  Ramps validated (dataviz validator, ordinal,
light + dark).  Tour stops are drawn from the tour's own file everywhere (copies of a minimum differ by flat motions /
slack up to ~0.56).  Checked: headless Firefox render, light mode; hover / click / dark mode not checked headless.
Weak / unresolved re-solve (10-11, `search/packer/resolve110.py`, `resolve110_sign.py`; outputs runs/resolve110/; store
backup ~/store-backups/store-2026-10-11c-pre-resolve110.sqlite): 61 targets (44 weak minima, 17 unresolved sides) in 112 s
on 15 procs.  Recorded (record_verdict, via resolve110): 2 weak -> certified, 1 unresolved -> certified, 1 unresolved ->
not-min; store now 778 certified / 42 weak below 11 (page still on the 775 snapshot).  15 unresolved: "contact equations
inconsistent" (projection stalls; likely unfinished near-misses, not exact configurations).  36 PSD-weak: the zero modes
are not flat (|dS| 1e-9..1e-7 at t = 1e-3 vs < 1e-30 for flat).  Signed dS (probe patched): 12 grow both ways on every
mode (higher-order minimum likely, unproved); 24 shrink somewhere, but only 3 clearly (dS < -1e-9), 3 in between, 18 tiny
(> -1e-15; the projection lets S float, so the sign may be the projection's choice).  Next: for dS < 0 cases, check
feasibility of the projected point and descend from it (a lower minimum = not-min); the 6 "negative curvature blocked"
cases need the contact-sign case split.
Page refresh (10-15): store after the exactsolve multiplier fix (1,076 certified minima below 11, store sha 6e1306b66996);
mapdata.py -> layout.py (rule picked seed 4, fill 47 %) -> export.py; tour / settle / forces unchanged.  Text: second
pocket island now also holds the best 41-tilted minima (fp cut 0.02); methods note on the 11 inconclusive minima and
the fixed solver bug.  Old 775 snapshot kept in runs/landscape110/map/snap775/.
