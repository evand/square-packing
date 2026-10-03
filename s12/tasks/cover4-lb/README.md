# cover4-lb: pin the lower end of COVER(4, closed) (2026-09-11)

**Why.**  Rung 2 (`tasks/rung2-s13`, `search/RUNG2.md`) is chasing an exactly verified closed
cover of `[0,4]²` with `W < 13`.  Its floor is only known as `≥ 12.008` — the `s = 3.99` measure of
`search/DUAL_EXACT.md` read at `s = 4`.  The heuristic cover LP sits at `12.39–12.46`.  Pin the floor:
an exact fractional packing of closed unit squares in `[0,4]²` (closed semantics: pairwise
closed-disjoint is the packing side; a *measure* only needs coverage `≤ 1`) of mass as large as
possible.  This is the pure LP dual (no cliques — cliques are not valid for the cover LP's dual;
only the point-coverage constraint is).  Conjecture: the true value is `≈ 12.1–12.3`, so the LP
proof of `s(13) = 4` costs about a quarter of a unit more than 12 — worth stating exactly.

**How.**  `search/dual_exact.py` (`build` from a float support, `check` exact) and the float
search it snaps (`search/packing_dual.py`, `search/nu_f.py`, `search/DUAL.md`).  Run the float
dual at `s = 4` with closed semantics (the container is `[0,4]²` closed; admissible iff the closed
square is inside), D4-symmetrised, column generation over poses (the pricer wants poses whose
capture under the current cover is smallest — the `closed4.py` rows are the right candidates);
then `dual_exact.py build` at `s = 4` (`Q = 10^5`, `Dc = 10^6`, clamping into `[w/2, 4 − w/2]`
exactly), `check` with all arrangement vertices, report `mass`, `M`, `L = mass / M` exactly.
Every `L` you certify is a theorem: no closed cover of `[0,4]²` has `W < L`.  Also record the
best *upper* value known (the heuristic covers' honest cost, `search/FAMILY.md` §1) so the note
states the bracket `[L, U]`.

**Inputs** (read-only, copy): `/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/`
has the cover side; the `3.99` support is `runs/dual_exact_3.99_support.txt` wherever
`DUAL_EXACT.md` says (search the worktrees under `.claude/worktrees/*/runs/` if it is not in the
main `runs/`; copy, never modify).

**Budget.**  4 processes; detach anything over 10 minutes (`setsid nohup`); never `pkill`; pids
921360–921362, 1109693–1109695, 1245763 and its children are others'.  Report within ~3 h.

**Deliverables.**  `search/COVER4.md` (statement proved up front, exact numbers, the bracket,
exact-vs-float, reproduce), the certificate `runs/dual_exact_4_support.txt` quoted in the note,
commit on your worktree branch.  Do not edit `TODO.md`, `README.md`, or other tasks' files.
