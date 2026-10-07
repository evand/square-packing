# Clique lever: is the full clique LP the missing 0.75 at t = 4?

**Why.**  `search/LEAF_CEILING.md` §5: the exact max-weight clique of the closed-intersection
graph of every certified `t = 4` leaf measure is **1.54–1.65** (complete branch and bound), while
the *certifiable* anchor-clique family `K(p, A)` reaches only 1.27–1.49 on the same measures.  So
the measures violate general clique inequalities `mu(K) <= 1` by ~0.6 and the certificate family
captures at most half of that.  The leaf LP with anchor cliques sits at 11.745 on leaf `01010101`
(T4LEAF), 0.255 short of closing.  Question: **if every clique inequality of the closed-intersection
graph were enforced, where would the leaf LP go?**  This is the value of the "QSTAB" relaxation
`max sum mu s.t. cov <= 1, mu(K) <= 1 for every clique K, region masses = leaf counts`, and it
bounds from above what *any* clique-based strengthening of the certificate can achieve on that pose
set — and bounds from below (by the finite pose set) the true QSTAB value.  Note the true QSTAB
value over all poses is `>= 11` always (11 pairwise disjoint closed unit squares fit in `[0,4]^2`,
e.g. the 4x4 grid minus a row and one more with spacing), so the interesting range is `[11, 11.75]`.

**How.**  Packing side only; no verifier, no Lean.  Build a row-generation loop, on a **fixed pose
set** (no pricing): start from `t4leaf.py`'s leaf LP state for `01010101` (poses, point rows,
anchor-clique rows; checkpoints at
`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a8e11c0e17344cba4/runs/tl_A01010101_{poses,dual,cliques}.txt`
and `tl_C01010101P_*`, read-only — copy into your `runs/`), solve with HiGHS (`t4leaf.Hi` or
your own), and separate general cliques with `search/leaf_ceiling.py`'s `max_clique_int` /
`max_clique_through` (integer weights, exact B&B, complete flag) on the closed-intersection graph
of the *support* of the current solution (positive-mass poses only — a clique row on the support
is a valid clique row for the whole pose set once you extend it to every pose that closed-intersects
all its members... careful: it is only a clique if you extend *maximally and consistently*; the
simplest sound choice is to add the row over the support clique only, which is a valid inequality
but weaker, then optionally extend it to a maximal clique greedily over all poses, checking
pairwise closed intersection exactly).  Iterate: solve → coverage rows (reuse `t4leaf`'s exact
arrangement-vertex certification) → max-weight clique on the support → add → repeat, until
`M <= 1` and max clique `<= 1 + 1e-6`.  The B&B is exponential in the worst case; use the time
budget and the `lb` incumbent; report whenever it is incomplete (then the row is still valid, only
the stopping test is not).

Do it in this order: (1) on the 234-pose certified support `lc_A0101.txt` (the
`leaf_ceiling` worktree, `/home/evand/math/square-packing/s12/.claude/worktrees/agent-abf4981871f0804cf/runs/lc_A0101.txt`),
which is tiny and gives the answer on that support in minutes; (2) on the ~8k-pose leaf set;
(3) if the value dropped a lot in (2), one pricing stage with `t4leaf.py`'s stratified pricer to
see whether new columns bring it back up (that is the T4SCREEN calibration failure mode).
Also run the matched **corner leaf `k = 4`** control (`tl_B40K_*`) at step (2) if time permits.

**The number that decides.**  QSTAB value on the 8k-pose leaf set vs 11.745.  ~11.0–11.3 ⇒ cliques
are the lever, and the job becomes finding a certifiable clique family rich enough (report the
*shape* of the violated cliques: how many members, do they share a point, a segment, are they
"three squares around a point" pinwheels, what do their anchor-family approximations miss).
≥ 11.6 ⇒ cliques are not the lever either.  Report the shapes either way.

**Budget.**  4 threads; runs detached (`setsid nohup`), checkpointed; never `pkill`; never touch
pid 740647.  Answer to (1) within the hour, (2) within ~6 hours.

**Deliverables.**  `search/CLIQUELEVER.md` (verdict up front, table with `M`/clique-max per stage,
the clique anatomy, exact statement of what is exact vs float, reproduce section);
`search/cliquelever.py`; commit on your worktree branch.  Do not edit `TODO.md`, `README.md`,
or other tasks' files.
