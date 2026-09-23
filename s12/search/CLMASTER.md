# cl-master: restricted master + continuous lattice pricing for `cliquelever.py` (2026-09-11)

Task: `tasks/cl-master/README.md`.  Code: `search/cliquelever.py` (`Lever.solve_master`,
`Lever.price_columns`, `Lever.lattice_price`, `Lever.rebuild`) and one condition in
`search/t4leaf.py` (`Hi.run`).  Read against `search/LPSPEED.md` — the measurement this ports to
the clique loop — and `search/CLIQUELEVER.md` §1/§3 for what the loop computes.  Runs in this
worktree's `runs/` (gitignored; every number is quoted here).

## 0. Verdict, up front

**`--master --no-warm` is 4–10x per iteration on the LPs the clique loop actually runs, and
`--lattice-every` is what turns the runs that had stopped moving back into converging ones.**
Measured on this machine, both paths at once, same inputs:

| instance | rows x columns | old loop, s/iteration | `--master`, s/iteration | ratio |
|---|---|---|---|---|
| leaf `01010101` (validation 2) | 14.7k x 8,087 | **46.7** (LP 45.2) | **12.2** (LP 11.0) | **3.8x** |
| corner control `k = 4` (validation 3) | 14.6k x 7,941 | **78.6** (LP 76.5) | **9.1** (LP 7.7) | **8.6x** |
| the `B40KM` launch instance | 37.1k x 9,669 | **259.9** (`cl_B40KL3.log` `s0.0`) | **25.3** (`cl_B40KM.log` `s0.0`) | **10.3x** |
| 234-pose support (validation 1) | 32.5k x 235 | 0.7 | 1.4 | 0.5x — too small to restrict |

Per run, on the pair that both converged: **1,635 s -> 440 s (3.7x)**, in three *fewer* iterations.
The master carries 1,500–3,800 of the 8k–10k loaded columns; the ipm is 0.5 s at 149 columns and
260 s at 9,669 (`t ~ n^1.5`, LPSPEED.md's law again), and the warm dual simplex that used to burn
`--lp-tlim` before every single ipm fallback is gone.

**All three recorded values reproduce exactly.**  `111/10` (`11.100000`), `373/33` (`11.303030`)
and `537/46` (`11.673913`) to every digit the loop prints; the exactly certified measures are
`11.099999992`, `11.303030287` and `11.673912994`, each within `2e-8` of its rational and within
`4e-9` of the recorded one, each re-verified by `leaf_ceiling.py check --anchor clique`
(`M = 1`, regions exact, max-weight clique `= 1` with a complete B&B, "(b) PROVED").  The
trajectories differ; the converged values do not.

**Float vs exact, unchanged (§1):** the certified number is always the final exact measure from
`Lever.finalize`, a rigorous lower bound on the continuum QSTAB; the intermediate LP value of a
*capped* master is neither an upper nor a lower bound on `QSTAB(P)`; and at termination the master
optimum equals the full-LP optimum on the loaded pose set, because the loop can only stop on an
iteration where no loaded column has reduced cost `> --price-tol`.

**And the runs this replaces had stopped converging at all.**  `cl_B40KL2.log` ends with 400
iterations at `LP = 11.785783` adding nothing, 0.4 s apiece; `cl_B40KL3.log` (pid 921360) was in
that state 5 minutes after it launched.  With `--lattice-every` the same resume prices the full
169,538-candidate lattice the moment separation stalls (63 s: 499 new poses, 72,959 new clique
memberships) and the next solve moves: `11.785783 -> 11.792824`, 12 fresh violated cliques to cut.

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
from `--lattice-every`, and the `--master-add` (default 2000) loaded columns of most positive
reduced cost; all rows are kept; ipm+crossover on the small model.  After the solve every loaded
column is priced and the improving ones are added, for up to `--master-passes` (default 2) passes.
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
  the full LP.  Capping at `--master-passes` and lifting the cap once the separation has nothing
  left keeps the cheap iterations cheap and still forces the exact full-LP optimum at the only
  place it matters, the stopping test.
* **Wide jumps, not many narrow ones.**  A single ipm solve on the B40KL2 resume (37,133 rows,
  9,669 loaded columns, `CLDBG=1`, one core) costs

  | master columns | 149 | 449 | 749 | 1,049 | 1,349 | 1,949 | 2,149 | 2,807 | 9,669 (full) |
  |---|---|---|---|---|---|---|---|---|---|
  | ipm+crossover | 0.5 s | 1.8 s | 3.0 s | 3.4 s | 7.2 s | 17.6 s | 23.3 s | ~60 s | **260 s** |
  | model rebuild | 0.4 s | 0.2 s | 0.2 s | 0.2 s | 0.3 s | 0.5 s | 0.5 s | 0.5 s | — |

  which is LPSPEED's `t ~ n^1.5` again (`260 x (2149/9669)^1.5 = 27 s` against 23.3 measured) — so
  the cost of reaching a 2,800-column master in four 1,000-column steps (134 s measured) is five
  times the cost of reaching a 2,149-column one in a single 2,000-column step (24.7 s measured), and
  all three produced the same LP value.  Hence the defaults `--master-add 2000 --master-passes 2`,
  and no add cap at all on the passes that close the master exactly.  Rebuilding the HiGHS model
  from scratch for each pass (`Lever.rebuild(cols)`) costs 0.2–0.5 s even at 37k rows, so there was
  no reason to complicate the code with incremental `addCols`.

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

Three supporting changes, all of which compute exactly what the old code computed and none of which
changes a default:

* `Poses.contains_rows` takes `pose_lo`/`col_lo` and returns only the new column block.  Poses are
  only ever appended and never move, so the old block of `A` is unchanged and a pricing pass costs
  `rows x new poses` instead of `rows x all poses`.  Both default to 0, which is the old
  whole-matrix behaviour the stage code still uses.
* `Poses.joins_all` extends a clique row to the new poses with one vectorised early-exit pass per
  member, instead of one `meets_all` call per candidate.  The closed-intersection margin is
  symmetric in the two squares (swapping them permutes `m1,m2 <-> m3,m4` and only flips signs inside
  `abs`), so it is the same float test on the same values in the same order; the greedy
  mutual-compatibility pass over the survivors is untouched.  `2,645 clique rows x ~1,200 new
  poses` was 3M python-level calls.
* `price()`'s `clique_cost` now carries a shrinking index array instead of a full-length boolean
  mask.  The member loop is up to 1000 deep on a wide clique row and the survivors are a handful, so
  the mask version did ~1000 numpy passes over all 169,538 candidates per clique row for up to 500
  rows.

Measured together on the B40KM instance (37,133 rows, 2,645 clique rows, 9,669 columns, one core):
one full in-loop lattice pass — 169,538 candidates priced with clique charging, 400 kept per region
plus 400 refined and the neighbours of the support, 499 new poses snapped in, the coverage rows
extended to them and **72,959 new clique memberships** — takes **63 s** (`cl_V5dry.log`: the
trailing seconds in a log line are CUMULATIVE stage time, `25s` at the end of `s0.0` and `88s` at
the end of its lattice pass), against **127 s** for the same pass before `joins_all`, and against
the ~3 hours a `--price` stage used to sit between two pricings.  In production on the leaf resume
(1,246 clique rows, 11,584 columns) a pass is **84 s** (`cl_A0101M.log` `s0.14`: `1884s -> 1968s`).

## 3. Validation

Both paths run from the same read-only inputs
(`/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/`, copied into this worktree's
`runs/`), on this machine, at the same time, two physical cores each (`--threads 2 --procs 2`),
load average 15–23 throughout from the three E1 resumes and the other agents — so the two columns
are contended equally and the ratios are the honest ones; the absolute seconds are 10–30 % above
what an idle machine would give.  `old` = today's default path (`--lp-tlim 10`), `new` =
`--master --no-warm` on the shipped defaults.

| pose set | run | LP value | exact certified measure | `M` | max clique | it | s/it (LP) | total |
|---|---|---|---|---|---|---|---|---|
| (1) `lc_A0101.txt`, 234 poses / 235 cols | V1old | `11.100000` | `2219999999/200000000 = 11.099999995` | 1 | 1 (complete) | 7 | 0.7 (0.5) | **28 s** |
| | V1new | `11.100000` | `1387499999/125000000 = 11.099999992` | 1 | 1 (complete) | 7 | 1.4 (1.0) | **32 s** |
| (2) leaf `01010101`, 7,845 poses / 8,087 cols | V2old | `11.303030` | `2260606057/200000000 = 11.303030285` | 1 | 1 (complete) | 34 | 46.7 (45.2) | **1,635 s** |
| | V2new | `11.303030` | `11303030287/1000000000 = 11.303030287` | 1 | 1 (complete) | 31 | **12.2 (11.0)** | **440 s** |
| (3) corner control `k = 4`, 7,829 poses / 7,941 cols | V3old | `11.673913` at it 38, still cutting a `1.0109` clique at it 40 | (not reached) | — | — | 41+ | 78.6 rising to 90.1 (76.5 -> 88.2) | **> 3,514 s** |
| | V3new | `11.673913` | `5836956497/500000000 = 11.673912994` | 1 | 1 (complete) | 53 | **9.1 (7.7)** | **560 s** |

Every one of the five certified measures was re-checked by `leaf_ceiling.py check --anchor clique`
(the independent code path): coverage `M = 1`, regions exact (corners `1,1,1,1`; slots
`0,1,0,1,0,1,0,1` in the leaf, free in the control), chord strips `3,3,3,3`, max-weight clique `= 1`
with a complete B&B, "(b) PROVED".

**The converged values reproduce exactly.**

| target | recorded (`CLIQUELEVER.md`) | old path today | `--master` today |
|---|---|---|---|
| `111/10 = 11.100000000` | `11.099999995` | `11.099999995` | `11.099999992` |
| `373/33 = 11.303030303` | `11.303030286` | `11.303030285` | `11.303030287` |
| `537/46 = 11.673913043` | `11.673912998` | (still running) | `11.673912994` |

Each float LP value printed `11.100000`, `11.303030` and `11.673913` — the three rationals to the
six digits the loop prints, i.e. `< 1e-6` — and each exactly certified measure is within `2e-8` of
its rational and within `4e-9` of the recorded one.  The spread between the three columns is the
`/10^9` rounding of `finalize` on three different (equally valid) optimal vertices of a degenerate
face, not an error: the old path and the master land on different 66–71-pose supports of the same
optimal face.  The trajectories differ, as the brief expected (V3new takes 53 iterations where the
recorded control took 47 + 13); the converged value does not.

**Speed.**  Per iteration, wall clock, same machine, same moment:

| instance | rows x columns | old | `--master` | ratio |
|---|---|---|---|---|
| (2) leaf | 14.7k x 8,087 | 46.7 s (LP 45.2) | 12.2 s (LP 11.0) | **3.8x** (LP 4.1x) |
| (3) control | 14.6k x 7,941 | 78.6 s (LP 76.5) | 9.1 s (LP 7.7) | **8.6x** (LP 9.9x) |
| B40KM resume (the launch instance) | 37.1k x 9,669 | **259.9 s** (`runs/cl_B40KL3.log`, `s0.0`) | **25.3 s** (`runs/cl_B40KM.log`, `s0.0`) | **10.3x** |
| (1) 235 columns | 32.5k x 235 | 0.7 s | 1.4 s | 0.5x |

and per run, on the only pair that both finished, **1,635 s -> 440 s = 3.7x** (with the master
taking three *fewer* iterations).  The control's old run is worse than that ratio suggests: its
per-iteration cost climbs with the clique rows (78.6 s over the first 21 iterations, 90.1 s over
41), and after 41 iterations and 3,514 s it was still cutting `1.0109`-mass cliques on the
degenerate `537/46` face that `--master` had finished with at 482 s.  (1) is the control on the other side: with 235 columns there is
nothing to restrict, the master is the whole LP, and the four rebuilds an iteration make it twice
as slow — `--master` is for the 8k–15k-column LPs it was built for.

**The other half of the win is not in that table.**  The three runs this replaces had already
stopped making progress: `runs/cl_B40KL2.log` ends with 400 iterations at `LP = 11.785783`,
`M = 1.000000003`, `kmax = 1.000000`, `+0` rows and `+0` cliques, 0.4 s apiece — the loop had
separated everything its pose set had and could only spin, because `M` sat `2e-9` above the
`1 + 1e-9` stopping threshold and `--price` only prices at a stage boundary (550 poses per ~3 h).
`runs/cl_B40KL3.log` (pid 921360) was in the same state 5 minutes after launch.  With
`--lattice-every`, the same resume prices the lattice as soon as the separation stalls, takes in
499 new poses and 72,959 new clique memberships in 88 s, and the next solve moves the value:
`11.785783 -> 11.792824`, with 12 new violated cliques to cut.

## 4. Reproduce

```sh
python3 search/cliquelever.py selftest      # pinwheel 1.5 -> 1, old path AND --master, exact
R=runs                                      # copies of the read-only inputs of the task brief
#   /home/evand/math/square-packing/s12/runs/inputs-2026-09-11/  (read-only; copy into runs/)

# (1) the certified 234-pose support: 111/10, old path then --master
python3 search/cliquelever.py run V1old --exact $R/lc_A0101.txt --corners 1111 \
    --patterns 01010101 --chord --rows0 30036 --threads 2 --procs 2 --ckpt 1
python3 search/cliquelever.py run V1new --exact $R/lc_A0101.txt --corners 1111 \
    --patterns 01010101 --chord --rows0 30036 --threads 2 --procs 2 --ckpt 1 --master --no-warm

# (2) the leaf pose set (CLIQUELEVER.md 3): 373/33
C="--threads 2 --procs 2 --clique-time 60 --cq-want 60 --cq-top 40 --cq-age 0"
python3 search/cliquelever.py run V2old --exact $R/lc_A0101.txt \
    --load-poses $R/tl_A01010101_poses.txt --load-rows $R/tl_A01010101_dual.txt \
    --corners 1111 --patterns 01010101 --chord --lp-tlim 10 $C
python3 search/cliquelever.py run V2new --exact $R/lc_A0101.txt \
    --load-poses $R/tl_A01010101_poses.txt --load-rows $R/tl_A01010101_dual.txt \
    --corners 1111 --patterns 01010101 --chord --master --no-warm $C

# (3) the corner control (CLIQUELEVER.md 0): 537/46
python3 search/cliquelever.py run V3old --load-poses $R/tl_B40K_poses.txt \
    --load-rows $R/tl_B40K_dual.txt --corners 1111 --patterns ........ --chord --lp-tlim 10 $C
python3 search/cliquelever.py run V3new --load-poses $R/tl_B40K_poses.txt \
    --load-rows $R/tl_B40K_dual.txt --corners 1111 --patterns ........ --chord --master --no-warm $C

# every final measure, through the independent checker
for t in V1old V1new V2old V2new V3old V3new; do
  python3 search/leaf_ceiling.py check $R/cl_${t}_exact.txt --corners 1111 \
      --patterns 01010101 --chord --anchor clique          # --patterns ........ for V3*
done

# the per-pass build/solve timings of the table in 2.2 (the B40KM instance, one core):
CLDBG=1 python3 search/cliquelever.py run V5dry --exact $R/cl_B40KL2_poses.txt \
    --resume-rows $R/cl_B40KL2_rows.txt --resume-cliques $R/cl_B40KL2_cliques.txt \
    --corners 1111 --patterns ........ --chord --row-pitch 0 --threads 1 --procs 1 \
    --cq-age 0 --cq-want 60 --cq-top 40 --clique-time 60 --master --no-warm \
    --lattice-every 5 --iters 12 --no-finalize
```

`runs/cmp.py TAG.log ...` prints the per-iteration LP and wall seconds of any of these logs.

The three successor runs are `runs/launch_master.sh`:

```sh
CL="python3 search/cliquelever.py run"
COMMON="--threads 1 --procs 1 --no-warm --master --lattice-every 5 --cg-want 400 \
        --clique-time 60 --cq-want 60 --cq-top 40 --cq-age 0 --iters 100000 --time 1000000 --ckpt 5"
launch() { tag=$1; shift; setsid nohup $CL "$tag" "$@" $COMMON > "runs/$tag.out" 2>&1 < /dev/null & }
launch B40KM  --exact runs/cl_B40KL2_poses.txt --resume-rows runs/cl_B40KL2_rows.txt \
              --resume-cliques runs/cl_B40KL2_cliques.txt --corners 1111 --patterns ........ \
              --chord --row-pitch 0
launch A0101M --exact runs/cl_A0101L2_poses.txt --resume-rows runs/cl_A0101L2_rows.txt \
              --resume-cliques runs/cl_A0101L2_cliques.txt --corners 1111 --patterns 01010101 \
              --chord --row-pitch 0
launch PUREM  --load-poses runs/tl_L2PURE_poses.txt --load-rows runs/tl_L2PURE_dual.txt \
              --corners .... --patterns ........
```

## 5. New flags, and what to reach for

| flag | default | what it does |
|---|---|---|
| `--no-warm` | off | `--lp-tlim 0`: straight to ipm+crossover, no warm dual simplex |
| `--master` | off | restricted master + reduced-cost pricing over every loaded column |
| `--master-add N` | 2000 | improving columns taken per pass (no cap on the closing passes) |
| `--master-passes P` | 2 | passes per iteration; `0` = always to convergence; the cap is lifted on any iteration whose separation found nothing |
| `--master-trim K` | 3 | zero-mass, non-positive-rc columns leave after `K` consecutive solves |
| `--lattice-every M` | 0 (off) | price the full lattice every `M` iterations and whenever separation stalls |

Every default reproduces the old behaviour exactly; `--master`, `--no-warm` and `--lattice-every`
have to be asked for.  `CLDBG=1` prints per-pass `cols / build s / solve s`.

Reach for `--master --no-warm` on anything with more than ~2,000 columns, and for
`--lattice-every 5` instead of `--price N` whenever the point is to keep converging rather than to
take one measured pricing step.  On a pose set small enough to solve whole (`lc_A0101.txt`'s 235
columns) leave both off.

## 6. What is left

* The **rows** are untouched: 37k coverage rows plus 2,645 dense clique rows are what a
  restricted master still carries, and they set the floor (`0.5 s` at 149 columns is not `0 s`).
  LPSPEED.md §"Row sifting" measured row sifting at parity on the cover LP and says that if rows
  must be cut it should be done *by ipm*; nobody has tried that here.
* The **exact-closure iteration is the expensive one**: after a lattice injection the master must
  price to convergence to decide termination, and on the B40KM instance that took 6 passes / 358 s
  (`cl_V5dry.log` `s0.1`) against 25 s for a capped iteration.  It is still cheaper than the 260 s
  the old loop paid *every* iteration, and it is the iteration that does the work, but a cheaper
  termination certificate (a single dual-feasibility repair rather than a re-solve) would pay.
* `M = 1.000000003` against a `1 + 1e-9` stopping threshold is what left `cl_B40KL2` spinning for
  400 iterations, and `--lattice-every` routes around it rather than fixing it.  The threshold
  itself (`run_stage`: `M <= 1 + 1e-9` with `bad` selected at `cov > DM * (1 + 1e-9)`) is a
  separate question for whoever owns the stopping test.
