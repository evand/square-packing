# rank-diag: what non-clique structure carries the residual excess of the QSTAB measures? (2026-09-11)

**Why.**  `CLIQUELEVER.md` §6.3: at convergence the QSTAB measures are essentially Helly (24 of 27
tight cliques are point cliques) yet sit at 11.32 / 11.67 / 11.70 — above any integral packing on
their own finite supports (12 pairwise-disjoint closed unit squares in `[0,4]²` would disprove the
conjecture).  So on each support the gap `mass − α(support graph)` is realised by valid
inequalities that are *not* clique inequalities: rank inequalities `mu(X) <= α(G[X])`, of which
the first non-clique ones are odd holes (`C5`: `<= 2`), odd antiholes and odd wheels.  The
certifiable version exists and needs no geometry: pieces `X_i = {S ⊇ A_i, S meets A_j for j ~ i}`
with a "guaranteed-meet" graph `H` on pieces give `α(∪ X_i) <= α(H)` by Lemma 0 of
`notes/clique-family.md` (cliques are `H` complete).  Before building that, measure whether it
would bite.

**What to compute** (exact geometry via `search/leaf_ceiling.py`: `closed_graph`, integer SAT with
`<=`; `max_clique_int` for reference), on each of the converged measures in
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/`: `cl_A0101L1_exact.txt` (11.321,
70 poses), `cl_B40KL1_exact.txt` (11.674, 102), `cl_A0101L2_exact.txt` (11.435), `cl_B40KL2_exact.txt`
(11.786), `cl_PUREQ_exact.txt` (11.702, 152):

1. `α` of the support's closed-intersection graph, exactly (max independent set on ≤ 152 vertices;
   a MIP or a branch-and-bound with the complement-clique bound is fine).  Report `mass − α`.
2. Every chordless odd cycle of length 5 and 7 (and odd antiholes of size 5 = `C5`, 7) with mass
   `> (k−1)/2`; odd wheels with hub mass + rim mass `> 2`.  Report the heaviest ten of each kind:
   members, mass, excess, where they are (centres, angles, regions per the leaf's frame), whether
   they sit around a grid vertex such as `(2,2)`.
3. Minimal violated rank subsets more generally: search for `X` with `mu(X) − α(G[X])` maximal
   (heuristic search is fine — the exact enumeration is not needed; report what you find and how you
   searched), and whether they are covered by the odd-cycle family.
4. If odd cycles bite: for the heaviest ones, is there a choice of anchors `A_1..A_5` (points or
   short segments, as in `clique-family.md`) such that piece `i` contains `A_i` and meets `A_{i±1}`,
   with consecutive pieces pairwise meeting *by construction* — i.e. is the certifiable anchor-graph
   version of the cycle a set of positive pose volume that still carries the mass?  Use
   `clique_family.py`'s transversal/meets machinery.  This step is optional; do 1–3 first.

**The number that decides.**  If the odd-cycle family recovers a substantial part of `mass − α`
(say `>= 0.1` on the leaf), the anchor-graph rank family is the next lever and gets built.  If the
violated rank structure is diffuse (large `X`, no small certifiable witness), say so — that is the
"nothing certifiable beyond cliques" verdict and it matters just as much.

**Budget.**  Minutes of compute; 4 threads; nothing detached needed.  pids 921360–921362 are other
runs — do not touch; never `pkill`.

**Deliverables.**  `search/rankdiag.py`, `search/RANKDIAG.md` (verdict up front; the table per
measure; anatomy; exact-vs-float statement; reproduce), commit on your worktree branch.  Do not
edit `TODO.md`, `README.md`, or other tasks' files.
