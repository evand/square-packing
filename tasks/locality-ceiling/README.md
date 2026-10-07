# locality-ceiling: how much can ANY local family ever prove? (2026-09-11)

**Why.**  Every certifiable family in this project is *local*: a point row, an anchor clique, a
general clique, the pentagon rank-2 family (`search/RANKDIAG.md`) — each inequality involves only
poses whose centres lie within a window of diameter ≈ 2.  The only global constraints are the
region counts of the branch tree.  `search/RUNG2.md` Theorem 1 is an obstruction theorem for a
class of verifiers; this task is the obstruction measurement for the class of *families*: on a
fixed pose set, the best value that every local valid inequality together can reach, as a
function of the window diameter `D`.  If the leaf sits at `α = 11` already for `D ≈ 2–2.5`, local
reasoning suffices in principle and the continuum is the only question; if it stays near 12 until
`D ≈ 4`, Bentz-style global counting is unavoidable and the plan changes now.

**The relaxation (window consistency).**  Pose set `P` (fixed), closed-intersection graph `G`
(`leaf_ceiling.closed_graph`, exact).  Windows `W_j`: a family of closed boxes of side `D` (say on
a pitch `D/4`, clipped to the container; also try discs).  For each window, `P_j = {S ∈ P :
centre(S) ∈ W_j}` and the packing polytope `Π_j = conv{ 1_I : I independent in G[P_j] }`.

    LOC_D(P) = max mass(μ)  s.t.  cov(x) ≤ 1 at every arrangement vertex (as `cliquelever`),
                                   region counts = the leaf's counts (pinned regions, chord rows as
                                   in `cliquelever`), and  μ|_{P_j} ∈ Π_j  for every window j.

Dantzig–Wolfe: variables `λ_{j,I} ≥ 0` over independent sets `I` of `G[P_j]`, `Σ_I λ_{j,I} ≤ 1`,
`μ_S = Σ_{I ∋ S} λ_{j,I}` for every `j` with `S ∈ P_j` (one linking row per (j, S)).  Column
generation: max-weight independent set in `G[P_j]` under the linking duals (exact B&B; the
window graphs are small — report their sizes).  `μ|_{P_j} ∈ Π_j` implies every valid inequality
supported in `W_j`, so `LOC_D ≤ QSTAB ≤` (cliques + pentagons) `≤ …` and `LOC_4 = α(P) = 11`.
Restricting poses lowers the value, so every number is a lower bound on the continuum `LOC_D`;
intermediate LP values with a row subset are upper bounds on `LOC_D(P)`.  Use HiGHS via
`t4leaf.Hi` (`--lp-tlim 0` semantics: straight to IPM) or your own; a restricted master over
the λ columns is natural (`search/CLMASTER.md` for the pattern).

**What to measure.**  On the recorded pose sets (read-only,
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/`; copy into your worktree's
`runs/`): the leaf `cl_A0101L2_poses.txt` (`--corners 1111 --patterns 01010101 --chord`; value
11.435 under all cliques, α = 11) and corner `k = 4` `cl_B40KL2_poses.txt` (`--patterns ........`;
11.786, α = 11).  Table: `D ∈ {1.5, 2, 2.5, 3, 3.5, 4}` × {boxes, discs} → `LOC_D`, converged, with the
final measure exactly certified for coverage and regions (`leaf_ceiling.py check`) and, for each
window, an exact check that the restriction decomposes into independent sets (the λ's are the
certificate — verify pairwise disjointness of every support set exactly).  Then the anatomy of the
gluing measure at the smallest `D` where `LOC_D > 11`: where the windows disagree, what the
"locally a packing, globally not" structure is.  Optional: the same at `t = 3.98` for calibration
against the certified 3.98 leaves.

**The number that decides.**  `LOC_2 ≈ 11` on the leaf ⇒ local families can close the leaf on
this pose set and the pentagon result is near the ceiling; `LOC_3 ≥ 11.5` ⇒ no local family can,
and the endgame needs global counting inequalities — say which.

**Budget.**  Packing side only; no verifier, no Lean.  2 LP threads per run, ≤ 4 runs at once
(the machine carries six E1 runs and an 8-process checker; pids 921360–921362, 1109693–1109695,
1245763 and its children are others' — never `pkill`, never touch them).  Detach anything over 10
minutes.  Report within ~5 h.

**Deliverables.**  `search/locality.py`, `search/LOCALITY.md` (verdict up front; the table; the
anatomy; exact-vs-float; reproduce), commit on your worktree branch.  Do not edit `TODO.md`,
`README.md`, or other tasks' files.
