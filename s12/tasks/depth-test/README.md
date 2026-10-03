# Depth test: branch the interior squares of the hardest slot leaf at t = 4

**Why.**  `TODO.md` "Next session, first".  At `t = 4` the family (corner branch + slot branch +
anchor cliques + chord lemma) leaves leaf `01010101` under corner leaf `1111` at 11.745 (cliques on,
`M = 1`, one clique violated by 8.6e-4) / 11.920 (cliques off, certified), still rising when
stopped (`search/T4LEAF.md`).  Whether this family can ever be a proof turns on whether one more
level of branching — on the four interior squares — has pruning power.  This task measures it.

**What to measure.**  Extend `search/t4leaf.py` (subclass or flag; do not change existing
behaviour) with **interior sub-region equalities**: split the interior `[1,3]^2` into the four
quadrants `[1,2]x[1,2]`, `[2,3]x[1,2]`, `[1,2]x[2,3]`, `[2,3]x[2,3]` and pin the mass centred in
each to a count `q_i`.  Boundary semantics must be the verifier-compatible one of T4LEAF §1.1 (a
pose within `--bnd-delta` of a quadrant boundary — `c_x = 2` or `c_y = 2` — gets one column per
adjacent quadrant; a pose at `(2,2)` gets four).  The children of the leaf are the count vectors
`q` with `sum q = 4` and each `q_i <= 2` (a quadrant's admissible centres lie in a box of side
`<= 1`, diameter `sqrt 2`, so at most 2 centres more than 1 apart — check and state this bound; if
you can only justify `q_i <= 3`, use that).  Quotient by the symmetry of pattern `01010101`
(it is C4-invariant) and list the child count.

Run, in order of expected hardness: the balanced child `1111`; then `2110`-type and `2200`-type
(one representative each, or all if cheap).  Warm-start from the leaf-search checkpoints
(`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a8e11c0e17344cba4/runs/tl_A01010101_{poses,dual,cliques}.txt`
and `tl_C01010101P_*`, read-only; copy into your worktree's `runs/`).  Cliques on (`--cliques
--cq-wall --cq-interior`), chord on, boundary duplication on.  Take each child to `M <= 1` and
`kmax <= 1` if at all possible; report the trajectory per stage exactly as T4LEAF §2 does, and
report the matched no-clique value as well.  Then certify the best measure of each child exactly
with `search/leaf_ceiling.py snap --polish` (add the quadrant equalities to it too if needed —
it is self-contained by design, keep it so).

**The number that decides.**  Drop per level = 11.745 − max over children.  ~0.3–0.7 in one level
⇒ the Bentz-shaped branch-and-bound has pruning power and becomes a finite engineering project;
≤ 0.1 ⇒ this family needs new mathematics.  Also report where the mass of the hardest child's
optimum sits (grid vs tilted, boundary mass) — the anatomy is as informative as the value.

**Budget.**  4 threads total (`--threads 4`, one process at a time or two at 2).  Launch every run
detached (`setsid nohup ... &`) with checkpoints so it survives your session; never `pkill`; never
touch the running `leaf_ceiling.py` process (pid 740647).  Aim to have the balanced child's number
within ~6 hours of wall clock; if a stage takes > 30 min, cut `--pose-max` rather than waiting.

**Deliverables.**  `search/DEPTH.md` (verdict up front, table, anatomy, reproduce section, and an
honest "which numbers are bounds" paragraph in the style of T4LEAF §3.5); code changes with a
regression check that the new flag off reproduces `t4leaf.py` bit for bit on a small instance;
commit on your worktree branch.  Do not edit `TODO.md`, `README.md`, or other tasks' files.
