# cl-master: restricted master + continuous lattice pricing for `cliquelever.py` (2026-09-11)

Task: `tasks/cl-master/README.md`.  Code: `search/cliquelever.py` (`Lever.solve_master`,
`Lever.price_columns`, `Lever.lattice_price`, `Lever.rebuild`) and one condition in
`search/t4leaf.py` (`Hi.run`).  Read against `search/LPSPEED.md` — the measurement this ports to
the clique loop — and `search/CLIQUELEVER.md` §1/§3 for what the loop computes.  Runs in this
worktree's `runs/` (gitignored; every number is quoted here).

## 0. Verdict, up front

(pending)

## 1. What is and is not a bound (semantics, unchanged)

The **certified number is always the final exact measure**: `Lever.finalize` rounds the float
masses down to `/10^9` per column, tops each pinned region back up on *exact* coverage and clique
slack, recomputes `M`, the max-weight clique (complete B&B) and the region masses in integers, and
writes `runs/cl_<TAG>_exact.txt`, which `leaf_ceiling.py check --anchor clique` re-verifies through
a code path sharing nothing with the loop.  That number is a rigorous lower bound on the continuum
QSTAB value of the leaf.  Nothing in this task touches it: the master changes only which columns
the float LP is allowed to use.

What *does* change is the status of the intermediate LP values.

* Old loop: every column of the pose set is in the LP, so a solve is the exact optimum over a
  finite subset of the rows, hence an UPPER bound on `QSTAB(P)` at every iteration
  (`CLIQUELEVER.md` §1).
* `--master`: the LP additionally omits columns, which lowers the value.  **An intermediate capped
  master value is therefore neither an upper nor a lower bound on `QSTAB(P)`** — it is an upper
  bound on the QSTAB of the *active* columns and a lower bound on the row-relaxed LP.  It is a
  diagnostic, not a bound.
* At **termination** the two coincide again.  The loop can only stop on an iteration where no
  loaded column has reduced cost `> --price-tol`; on such an iteration the master's dual is
  feasible for the dual of the full LP over all loaded columns, so the master optimum IS the
  full-LP optimum on the loaded pose set (to the LP's own tolerance), and the old reading of the
  final LP value returns unchanged.  To make that reachable, the `--master-passes` cap is lifted
  automatically on any iteration whose row and clique separation found nothing (`force_full_master`
  in `run_stage`), so the last iteration always prices to convergence.

`--lattice-every` only enlarges the loaded pose set, which can only raise `QSTAB(P)` — the
direction of `T4SCREEN.md` §1.2: restricting the poses lowers the value, so any `QSTAB(P)` is a
lower bound on the continuum value, and a bigger `P` is a better one.  Rows generated from a master
solution are valid whatever the master was: a coverage row is an arrangement vertex of an actual
measure and a clique row is verified pairwise before it enters the LP, so nothing about soundness
depends on the LP being solved over all columns.

## 2. What was built

### 2.1 `--no-warm` / `--lp-tlim 0` — skip the warm simplex (`t4leaf.Hi.run`)

`Hi.run` tried a warm dual simplex on every incremental solve and fell back to ipm+crossover when
it hit `--lp-tlim`.  On these LPs the fallback is what runs *every single time*
(`runs/cl_B40KL2.log`: `HiGHS simplex -> kTimeLimit after 10s; ipm+crossover`, `it0` — zero pivots
in the 10 s), which is LPSPEED's finding for the cover LP restated: warm dual simplex on this
degenerate geometry is a 40–150k-pivot walk whatever basis it starts from.  One condition now sends
`tlim <= 0` straight to ipm+crossover:

```python
if not self.basis or self.tlim <= 0:      # was: if not self.basis:
```

`--no-warm` sets `--lp-tlim 0`.  Every default is unchanged (`Hi`'s own default is 300 s,
`cliquelever`'s `--lp-tlim` 60 s), so this is `--lp-tlim` × 1 second per solve saved, no more —
10 s an iteration on the recorded runs.

### 2.2 `--master` — the restricted master (`Lever.solve_master`, `Lever.price_columns`)

Follows `BModel.solve_restricted` in `search/branch.py`.  The LP is a maximisation with `<=` rows
and `>= 0` columns, so the reduced cost of a column under the row duals is

```
rc_j = 1 - sum_r y_r A[r,j] - sum_k z_k K[k,j] - lam(region of j) - sum_d chd_d D[d,j]
```

— exactly what `price()` already evaluates on new lattice candidates, but on the loaded columns,
with their own region label and their own clique memberships (a clique row charges a column iff the
column is a member, and the memberships are stored in `self.K`).  `rc_j > 0` is the improving
direction.  `price_columns` is two sparse mat-vecs plus a gather, `0.05–0.3 s` on the biggest
instance here.

Per iteration: active columns = the current support (columns at positive mass), whatever entered
from `--lattice-every`, and the `--master-add` (default 1000) loaded columns of most positive
reduced cost; all rows are kept; ipm+crossover on the small model.  After the solve every loaded
column is priced and the improving ones are added, for up to `--master-passes` (default 4) passes.
Columns at zero mass for `--master-trim` (default 3) consecutive solves **and** non-positive reduced
cost leave the master; they stay loaded and are priced again every solve.  The first master of a run
is the support of the loaded masses (a `t4screen` or `cl_*` checkpoint carries 400–1100 positive-mass
poses among its 8–9.5k), not the whole set.

Two things had to be got right, both measured:

* **The trim rule.**  `branch.py` trims on mass and age alone.  Here that cycles forever: the LP is
  massively dual-degenerate, so a zero-mass column can carry `rc = +0.4` under the master's own
  optimal dual; it gets trimmed, priced straight back in, trimmed again.  On the 235-column
  `lc_A0101.txt` instance the loop ran 330+ iterations at `LP = 11.100000` with
  `rc+` oscillating 4 → 42 → 4 and never terminated.  Requiring `rc <= price_tol` as well fixes it.
* **The pass cap.**  Pricing to convergence every iteration (`--master-passes 0`) reproduces the old
  trajectory exactly — the first iterations of the leaf and control validations came out at
  `LP = 11.952147` and `12.000000`, the recorded values to the last digit — but needs 7–11 ipm
  solves per iteration while the duals chase new rows, and the last of those solves is as wide as
  the full LP.  Capping at 4 and lifting the cap once the separation has nothing left keeps the
  cheap iterations cheap and still forces the exact full-LP optimum at the only place it matters,
  the stopping test.

Infeasibility: a restricted master is infeasible only through the pinned-region equalities
(`mu(R) = k_R` with no active column in `R`).  `cover_regions` covers such a region greedily with
its 64 best-priced columns before every solve, and a failed solve escalates — first every column of
every pinned region plus the 5000 best-priced, then the full column set — with a line in the log.
It did not fire in any run here.

The stopping test becomes: no violated arrangement vertex on the support, `M <= 1 + 1e-9`,
max clique `<= 1 + --ktol` with every B&B complete, **and** no loaded column with `rc > --price-tol`.

### 2.3 `--lattice-every M` — continuous pricing (`Lever.lattice_price`)

Every `M` iterations, and additionally on any iteration where neither a coverage row nor a clique
row was separated (so a converged-but-not-terminated loop can never spin), the existing stage pricer
runs inside the loop: the `--price-pitch`/`--price-dth` lattice (169,538 candidates at the defaults)
with the same clique charging, the `--cg-want` best poses snapped and added to the LOADED set, the
coverage rows extended to the new columns and every clique row extended to the new poses that
closed-meet all of its members — exactly what the stage boundary in `cmd_run` does, but without
ending the stage.  `--price` (stages) is untouched and still works as before; `--lattice-every 0`
(the default) is off.

Two supporting changes:

* `Poses.contains_rows` takes `pose_lo`/`col_lo` and returns only the new column block.  Poses are
  only ever appended and never move, so the old block of `A` is unchanged and a pricing pass costs
  `rows x new poses` instead of `rows x all poses` (seconds instead of minutes at 30k+ rows).  Both
  default to 0, which is the old whole-matrix behaviour the stage code still uses.
* `price()`'s `clique_cost` now carries a shrinking index array instead of a full-length boolean
  mask.  The member loop is up to 1000 deep on a wide clique row and the survivors are a handful, so
  the mask version did ~1000 numpy passes over all 169,538 candidates per clique row for up to 500
  rows.  Identical arithmetic on identical values, same result, minutes instead of tens of minutes.

## 3. Validation

(pending)

## 4. Reproduce

(pending)
