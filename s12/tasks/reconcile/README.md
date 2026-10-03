# J. Reconcile the two sides: is the anchor-clique gain real on the cover side?

**The discrepancy.**  Task G (`search/CLIQUE_CONTINUUM.md`) measured, on the packing side at
`t = 3.99`, a stable anchor-clique gain of **0.10–0.12** on a pose lattice that reproduces the pure
value (pure 12.01–12.03, clique 11.90 ± 0.006, flat over an 8.5× pose refinement).  Task I
(`search/ANCHOR.md` §4–5) built the cover-side object — verified, Lean-checked — and measured, on the
`k = 4` leaf at 3.98 with matched pairs (same rows, same columns, cliques on/off), a gain of
**0.0000, 0.0000, 0.0004**, with 1180 cliques offered.

**Why the explanation in `ANCHOR.md` §4 ("a clique column costs the same as the point column it
contains, so it only pays on a strip of width ε") cannot be the whole story.**  On one and the same
finite instance — a fixed finite set of poses (cover rows = packing variables) and a fixed finite
family of cliques (cover columns = packing constraints), with the point columns being the cliques
`K(p, {p})` — the cover LP and the packing LP are exact duals and have the **same value**.  So a
cover-side gain of 10⁻⁴ and a packing-side gain of 10⁻¹ cannot both be right *on the same instance*.
The two measurements used different pose sets (G: an exact-centre lattice plus priced poses, D4
orbits, ~4.7k orbits; I: verifier witnesses of the current cover, ~35–80k single poses over
`[0°, 90°)`), and different separation depth (G: clique cuts separated repeatedly against the evolving
measure until none is violated; I: one batch of priced cliques per round, interleaved with row
generation and restricted-master churn, and the LP is dual-degenerate enough to dodge a single
batch — exactly what `search/BOXCLIQUE.md` saw).  Which of these accounts for the gap decides the
whole clique line, and it can be decided cleanly.

Read first: `search/CLIQUE_CONTINUUM.md` (Verdict, §2, §3), `notes/clique-family.md` (Lemma 0–2, §5),
`search/ANCHOR.md` (all), `search/BOXCLIQUE.md` ("The negative result"), `search/TIGHTEN.md`,
`search/LPSPEED.md` (Recommendation), `search/DUAL.md`, `search/DUAL_EXACT.md`.  Code:
`search/clique_continuum.py`, `clique_family.py`, `clique_exact.py` (G), `search/anchorclique.py`,
`anchorsep.py`, `search/branch.py` (I; `--matched`, `--lam-lo 0 --lam-hi 0` for the pure cover),
`search/tighten.py`, `verify/`, `xcheck.py`.

## The experiment, at `t = 3.99` (pure cover, no multipliers)

**Step 1 — same instance, both sides.**  Take task G's final calibrated pose set at 3.99 (the last
stage support of the `cq_*` run in `/home/evand/math/square-packing/s12/runs/`, ~4.7k orbits; expand
the D4 images) as the **fixed row set** of the cover LP.  Point columns: price at the arrangement
vertices of those poses (the same "coverage ≤ 1 at every arrangement vertex" that the packing side
imposed), or start from G's certified support points plus the shipped certificate's points.  Clique
columns: price with G's own separator (`clique_family.py` max-anchor-clique on the cover's
D4-averaged interior-point dual, which is a packing measure on exactly these poses), and **iterate
pricing to convergence** — solve, price cliques, add, re-solve — until no clique has `ȳ(K) > 1 + 1e-6`.
Report the matched pair (same rows and columns, cliques on/off) at every iteration.  **Expected by
duality: the converged cover value equals G's packing-side clique value on the same poses (≈ 11.91)
and the pure value equals ≈ 12.02.**  If it does not, one of the two implementations is wrong —
find which (G's separator vs I's credit rule; compare `μ(K)` for the same clique computed by both).

**Step 2 — let the verifier speak.**  From the converged step-1 cover (points + cliques), run the
exact verifier at `N = 2000` as separation oracle (`tighten.py` / `branch.py --lam-lo 0 --lam-hi 0`
with the anchor block, restricted master), adding witness rows and re-pricing cliques to convergence
inside each round.  Track: LP value, probe minimum, matched-pair gain, number and location of the
new violated poses (are they poses G's pricer never generated?  where are they — wall band, interior,
which angles?).  Two outcomes, both decisive:

* the value climbs back to `≥ 12` and stays there as rows converge → the packing side's flatness
  was an artefact of incomplete pose pricing, `V(3.99) ≥ 12`, and anchor cliques do **not** beat the
  pure ceiling at 3.99.  Then measure how far below 3.99 they do work (3.98 pure, then the leaf) —
  that is still the number that matters for `s(12) ≥ 3.98`.
* the value stays `< 12` and the verifier runs out of violated rows → **verify** the certificate at
  `N = 6000` and `xcheck.py --all`, ship it as `certificates/s12_lower_3.99_anchor.txt`: that is
  `s(12) ≥ 3.99`, past the pure ceiling, the first bound from cliques.

**Step 3 (only if step 2 is negative).**  The same two steps on the `k = 4` leaf at 3.98 with task
G's leaf pose set (`cq_LA98*`, `cq_LP98b*`) — where G's pure value is itself unconverged (11.954 vs
12.000), so calibrate first.

## Compute

The machine is free: up to 12 threads and 60 GB; pin with `taskset` (cores `n`, `n+16` are
siblings).  One background job of task I's may still be running (`branch.py … t398ik4p`, one core);
leave it.  `runs/` is gitignored: read `/home/evand/math/square-packing/s12/runs/` by absolute path,
write your own `runs/`.  **First check your worktree base** (`git merge-base branch-corners HEAD`
vs `git rev-parse branch-corners`; `git merge branch-corners` if they differ).

## Done when

`search/RECONCILE.md`: the step-1 table (iteration, cliques, pure value, clique value, gain, max
`ȳ(K)`), the duality check against G's numbers with the discrepancy explained, the step-2 table
(round, rows, new violated poses and where, LP value, probe min, matched gain), and the verdict —
which of the two outcomes, with the certificate shipped and verified if the second.  State for every
number whether it is verified, heuristic, or a bound in which direction.  Commit in your worktree;
do not touch `TODO.md`, `README.md` or other tasks' files.
