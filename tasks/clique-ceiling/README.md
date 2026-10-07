# A. Clique ceiling (packing side)

**Goal.** The value of the clique-strengthened relaxation as a function of `t`, computed on the
packing side with pose column generation, certified exactly — i.e. the ceiling of any
clique-certificate method, or its absence below 4.  Go/no-go for "cliques alone reach 4".

**Semantics (important).** A certificate at container `t` refutes packings at side `t' < t`; rescaled
to `[0,t]²` those are 12 pairwise **disjoint closed** unit squares.  So a valid clique is a set of
poses that pairwise **closed-intersect** (touching counts), and a certified packing measure must
satisfy `sum_{S in K} y_S <= 1` for every such `K`.  `clique_check.py` uses a strict separating-axis
test; a max-clique computed that way can under-report (a touching pair inside a larger clique is
missed) and a "ceiling" built on it is not one.  Use closed intersection (`<=` in the SAT test, exact
rationals) for every clique used in certification.  Pairwise touching alone is already imposed by
the point constraints (a touching pair shares a point), so the difference lives only in non-Helly
cliques.

**Plan.**
1. Extend `search/packing_dual.py` (LP over poses, point constraints, pose pricing) with clique rows:
   separation = max-mass clique of the current measure under closed intersection (adapt
   `clique_check.max_weight_clique`), plus the best clique through each of the ~40 heaviest poses
   per round (`clique_lp.py --per-round` style; mass migrates to neighbouring cliques otherwise).
2. Keep the certification step of `dual_exact.py`: coverage `<= 1` at every arrangement vertex, and
   now also max closed-clique mass `<= 1` on the support, all exact.  A certified clique-feasible
   measure of mass `>= 12` at `t` is a ceiling: no clique certificate below 12 exists at any `t' >= t`.
3. Run at `t = 3.99, 3.995, 4.0` (closed) and the corner leaf `k = 4` at `3.98` (`--branch --kmass 4`:
   corner mass exactly 4).  Stop each when either the certified mass crosses 12 or the pricing gain
   is `< 0.1 %` with mass `< 12`.

**Inputs.** Certified measures `runs/dual_P{A2,C1,D1,D2}_support.txt`, the leaf dual
`runs/branch_t398hk4_dual_it16.txt`, code `search/packing_dual.py`, `search/nu_f.py`,
`search/dual_exact.py`, `search/clique_check.py`, `search/clique_lp.py`.  (`runs/` is gitignored:
read from the main tree, write to your own `runs/`.)

**Compute.** Small LPs (1–3k poses); a few cores.  Max-clique is seconds at 1,500 poses.

**Done when.** A table `t → certified clique-feasible packing mass (closed cliques), pricing gain`
including the `k = 4` leaf, written up as `search/CLIQUE_CEILING.md`, with the crossing of 12
located to 0.005 if there is one.  State clearly which numbers are certified and which are
heuristic (unconverged in which direction).

## Status (2026-08-30)

Done as far as the packing side goes; see `search/CLIQUE_CEILING.md`.  **No clique-feasible
measure of mass >= 12 was found at t = 3.99, 3.995, 4.0 or in the k = 4 leaf at 3.98**, so no ceiling
is established and the clique method is not excluded below 4.  Certified (exact, closed-intersection
cliques, `--check` re-verified): a clique-feasible measure of mass 125/11 = 11.3636 at 3.99.
Heuristic fixed-pose clique-LP values: ~11.2–11.45 (3.99), ~11.36–11.48 (3.995), ~11.8–11.85 (4.0,
unconverged), ~11.65–11.8 (leaf).  Structural finding: the non-Helly excess of the pure measures is
carried by grazing contacts (pair margins down to 1e-5), so box cliques must be very thin.  Next:
the cover side (task B) is where the answer can actually be certified.  `clique_check.py` still uses
the strict test; `clique_ceiling.py` has the closed one (float and exact).
