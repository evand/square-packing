# s11-cliques: the clique-strengthened packing optimum at `n = 11` in the window `[3.8143, 3.83375]`  (drafted 2026-09-13, not yet launched)

`search/N11_ANATOMY.md` §7.  The pure cover LP is pinned at exactly `11` across the window (rigorous ceiling
`ν_f(3067/800) ≥ 11`); the optimum violates a non-Helly clique by `0.25` on eight interior poses; rank headroom
is `mass − α = 1.000` exactly; the corner branch is vacuous (corner boxes saturated at mass `1`).  Run
`search/clique_ceiling.py` with `--n 11` at `t = 3.82, 3.825, 3.83` to measure the clique-strengthened packing
optimum; wherever it is `< 11`, build a box-clique certificate with `search/boxclique.py` (already in
`certificates/FORMAT.md`, the verifier and Lean) and verify at `N = 6000, 12000`.  Deliverable: a new rigorous
`s(11) ≥ t`, or the number saying how much of the `1.0` cliques reach.  ≤ 4 threads, an afternoon.  The `0.75`
cliques cannot reach is the rank-10 statement (ring of eight pinned, only two interior squares fit where the LP
puts three) — the margin-positive rehearsal for `tasks/leafa-proof`.
