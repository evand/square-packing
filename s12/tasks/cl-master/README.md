# cl-master: restricted master + continuous pricing for `cliquelever.py` (2026-09-11)

**Why.**  Every `cliquelever` iteration is ~98 % LP: 120–195 s of HiGHS IPM on the full 8.5–15k-column
model, preceded by a 10 s warm dual-simplex attempt that makes zero pivots every single time
(`runs/cl_B40KL2.log`: `HiGHS simplex -> kTimeLimit after 10s; ipm+crossover`, `it0`).  The support
is 180–240 columns.  `search/LPSPEED.md` measured all of this on the cover LP and *built* the fix
for `branch.py` (`BRANCH_SOLVER=restricted`, `BModel.solve_restricted`): IPX is single-threaded with
cost `∝ ncols^1.5`, warm simplex is a degenerate crawl in every form, row sifting is parity — so
keep IPM and shrink the master.  `cliquelever.py` never got it.  Three E1 pricing resumes are
running right now on the old loop (`runs/launch_2026-09-11.sh`: pids 921360, 921361, 921362 —
**do not touch them, never `pkill`**); the goal is to replace them with a loop that is 3–8× faster
per iteration and prices the full lattice continuously instead of 550 poses per 3-hour stage.

**What to build** (in `search/cliquelever.py` + `search/t4leaf.py`'s `Hi`; every existing default
must reproduce the old behaviour, new behaviour behind flags):

1. `--lp-tlim 0` (or a `--no-warm` flag) skips the warm simplex and goes straight to IPM.
2. `--master`: a restricted master.  Active columns = current support ∪ columns entering this
   iteration ∪ up to `--master-add N` (default 1000) columns with the most positive reduced cost
   under the last duals, priced over *all* loaded columns each iteration (a sparse mat-vec; clique
   rows charge a column iff it is a member of the row — memberships are already stored).  Trim
   columns at zero mass for `--master-trim K` (default 3) consecutive solves.  The loop's stopping
   test becomes: no violated arrangement vertex on the support, max clique `<= 1 + ktol` with every
   B&B complete, **and** no loaded column with reduced cost `> price_tol`.  Follow
   `BModel.solve_restricted` in `search/branch.py` for the mechanics (infeasible master iff a pinned
   region has no active column — cover it greedily).
3. `--lattice-every M`: every `M` iterations (default 0 = off), price the `price_pitch`/`price_dth`
   lattice (the existing stage pricer, 169,538 candidates with clique charging) and add the
   `--cg-want` best columns to the *loaded* set, extending clique rows to them as the stage code
   does.  This turns "stages" into one continuously converging run; keep `--price` working as before.

**Semantics to preserve** (state them in the note): the certified number is always the final exact
measure (`Lever.finalize` → `leaf_ceiling.py check`), a rigorous lower bound on the continuum
QSTAB; intermediate LP values of a capped master are neither bounds; at termination the master
optimum equals the full-LP optimum on the loaded pose set (no positive reduced cost).

**Validation, exact, before any launch.**  On the recorded pose sets under `--master`:
`runs/inputs-2026-09-11/lc_A0101.txt` → `111/10` (27 s); the leaf step (2) of `CLIQUELEVER.md` §3
(`tl_A01010101_{poses,dual}.txt` + `lc_A0101.txt`, `--corners 1111 --patterns 01010101 --chord`)
→ LP `373/33 = 11.303030` to `1e-6` and an exactly certified measure within `1e-7` of it; the corner
control (`tl_B40K_*`, `--patterns ........`) → `537/46 = 11.673913`.  The trajectory will differ;
the converged value must not.  Record old vs new seconds per iteration and per run.
Inputs: `/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/` (read-only; copy into your
worktree's `runs/`).  `leaf_ceiling.py check --anchor clique` on every final measure.

**Then launch** (detached, `setsid nohup`, in your worktree's `runs/`, 1 core each, `--lattice-every`
on): `B40KM` from `cl_B40KL2_{poses,rows,cliques}.txt` (`--corners 1111 --patterns ........ --chord
--row-pitch 0`), `A0101M` from `cl_A0101L2_*` (`--patterns 01010101`), `PUREM` from
`tl_L2PURE_{poses,dual}.txt` (`--corners .... --patterns ........`, no chord).  These are the
successors of the three running resumes; the machine has 16 physical cores and ~13 are idle.

**Deliverables.**  `search/CLMASTER.md` (verdict up front: speed-up per iteration and per run,
the validation table, exact-vs-float statement, reproduce section), the code, the launch script
`runs/launch_master.sh`, commit on your worktree branch.  Do not edit `TODO.md`, `README.md`, or
other tasks' files.  Budget: validation within ~3 h; the launches after that.
