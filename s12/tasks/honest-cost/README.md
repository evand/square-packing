# honest-cost: what does the best pure `t = 4` dual actually cost as a cover? (2026-09-12)

**Why.**  Every `t = 4` number in the repo is packing-side: `QSTAB(P) + rank rows` on a finite
pose set `P`, a lower bound on the continuum relaxation.  The pure instance is bracketed
`[11.800417910, 11.865961]` on E2Pg's 956-pose support (`search/RANKDIAG.md` §15–17), 0.13 below
12 — but the rung-2 loop measured the tax between "LP on rows" and "valid cover at `t = 4`":
LP `12.339` → honest cost `12.56` (min capture `0.982`), i.e. +1.8 % before any checker slack
(`search/RUNG2.md` §0, §10), and the pure dual's lattice reduced cost floors at `0.05–0.07`
(§17), i.e. it is 5–7 % short somewhere on the 0.04 lattice.  `11.87 × 1.018 = 12.08`.  So E1 —
"is the continuum value of pure QSTAB + polygons below 12?" — is on a knife-edge, and the number
that decides whether V1 (the disjunctive `t = 4` verifier with cliques and polygons) is worth
building this week is the **honest cost of the best current dual as a cover**:

    honest(w) = Θ(w) / min over ALL admissible poses of capture_w(pose)

for the dual `w` = (point weights, clique weights, polygon weights) of the converged pure LP.

**Semantics** (do not change): `t = 4`, closed unit squares, closed containment, a point on `∂Q`
counts (`search/ZEROMARGIN.md` §1).  A pose is *admissible* iff its closed square lies in `[0,4]²`.

**The cover.**  A dual of the packing LP of `search/CLIQUELEVER.md` §1 (+ `RANKDIAG.md` §8) is a
cover with accounting `Θ = Σ_p y_p + Σ_K z_K + Σ_Π (k−1)/2 · ζ_Π`, and the capture of an
*arbitrary* pose `S` (not only poses of `P`) is

    capture(S) = Σ_{p ∈ S} y_p  +  Σ_K z_K [S ∈ K]  +  Σ_Π ζ_Π [S ∈ ∪ X_i(Π)]

where the membership rules are the *certifiable* ones, the ones a verifier would use:
  * a point `p` counts iff `p ∈ S` (closed);
  * a clique `K` counts iff `S` closed-meets every member of `K` — the rule `cliquelever.py`'s pricer
    already uses (its docstring around line 1224: "a candidate is charged a clique's dual iff it
    closed-meets every member");  *note this is a valid-but-generous rule*: also report the
    stricter box/anchor-clique rule if the run's cliques have recorded cores (`BOXCLIQUE.md`) —
    if the two differ materially, say so;
  * a polygon `Π` with anchors `A_0..A_{k−1}` counts iff `S` contains both endpoints of some side
    `[A_i, A_{i+1}]` (`search/rankfamily.py` docstring).

**Inputs** (read-only; copy into your worktree's `runs/`):
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-12/` — `cl_SUPEPGF_*` (the certified
956-pose support state, LP `11.800418`; `_poses.txt`, `_rows.txt`, `_cliques.txt`, `_pgons.json`
with exact rational anchors, `_measure.txt`, `.json`) and `cl_E2Pg_*` (the full 12k-column loaded
set, LP `11.865961`).  The duals are not written to disk as such: re-solve the final LP of each
state on its fixed pose set with `cliquelever.py` (`--price 0 --lattice-every 0`, resume rows /
cliques / `--resume-pgons`, `--iters` small) and take `y`, `z`, `pgz` from `Lever.read_duals`;
add a `--dump-dual` (or similar) that writes them.  If dual degeneracy makes the dual non-unique,
say so and take the one HiGHS returns; the number is an *upper* bound on the honest cost of the
best dual on that pose set either way.

**Approach.**
1. Rebuild each dual as a cover; check it reproduces `Θ = LP` and `capture ≥ 1 − 1e-9` on every
   pose of `P` (a self-test; if it fails, the crediting rule is wrong — fix that first).
2. Minimise `capture` over admissible `(c_x, c_y, θ)` in floats: (a) evaluate on the 0.04 / 2.5°
   lattice and on the `family_rows.py` families (item 3 at pitch 0.004, `--tilted`), (b) local
   descent (Nelder–Mead or coordinate search; capture is piecewise-constant in the centre for
   points, so use a small-step pattern search with restarts) from the 200 worst lattice poses and
   from every recorded arrangement vertex of `P`, (c) a dense final scan (≥ 10⁶ poses) to confirm.
   D4 symmetry of the container is NOT a symmetry of the dual unless the run enforced it — check
   `sym` in the pose file header and treat the full domain.
3. Report `min capture`, the minimising poses (location, angle, what they capture from each column
   type), `honest = Θ / min capture`, for both SUPEPGF and E2Pg duals, under both clique rules.
4. Optional if time permits: one round of *cover-side* repair — add the 400 worst poses as rows,
   re-solve, re-scan — and report how much the honest cost moves.  Do not iterate further.

**Reading.**  `honest ≥ 12.3`: the current family is short; V1 is not urgent; the next work is the
family.  `honest ≤ 12.1`: build the cover loop.  In between: report the minimisers' anatomy, which
is what decides the next family.

**Budget.**  Floats only (this is a packing-side measurement of a cover; nothing is certified).
Up to 8 processes; launch anything over 10 minutes detached (`setsid nohup … &`, `< /dev/null`);
never `pkill`; pids 1109693–1109695 are other runs — do not touch.  Report within ~4 h even if
not done, with the exact state.

**Deliverables.**  `search/HONEST.md` (verdict up front; the table; the minimisers; what is exact
vs float; reproduce section), `search/honestcost.py`, commit on your worktree branch.  Do not edit
`TODO.md`, `README.md`, `verify/`, `lean/`, `certificates/`, or other tasks' files.
