# C. LP speed

**Why.** A leaf of the corner branch is 3–8 h because the cover LP is 150–190k rows × 10–13k
columns (~800 nonzeros per row) and HiGHS IPM takes 1–2 h per solve at that size (superlinear).
Everything cover-side that remains — certifying `1110`, `k = 4` with clique columns, level 2 —
pays this per round.  Target: 10× per round.

**Step 0: get an instance.** `search/branch.py` has no LP dump (rows live only in memory;
`lp_search.py --dump-lp` exists to copy).  Add `--dump-lp PATH` that writes, at the round boundary,
`c`, the sparse `A` (scipy `.npz`), `b`, bounds, and the row keys.  Then produce (a) a mid-size
instance (~30–60k rows) by restarting `1110` from `runs/branch_t398j1110_cols_it51.txt` /
`_probe_it51.txt` for 3–4 rounds, and (b) a large one (150k+) — either continue that restart or
generate rows by the lattice scan against the `k = 4` columns `runs/branch_t398hk4_cols.txt`.
Record the reference optimum (IPM) for each instance.

**Step 1: benchmark**, each pinned (`taskset`, HiGHS `threads`) and timed per solve, same tolerances:
- HiGHS IPM (baseline; the current default), HiGHS PDLP (`solver=pdlp` in recent HiGHS), HiGHS
  dual simplex cold.
- Incremental dual simplex via `highspy`: build the model once, `addRows` per round, `run()` with
  the carried basis.  The existing `warm` backend was 200–600 s where scipy took 0.4 s; find out why
  (rebuilding per call? dense conversion? basis size mismatch?) before concluding anything.
- **Row sifting** against the archived rows: solve on the near-tight subset, scan the archive for
  violated rows (a sparse mat-vec), add, repeat to convergence.  Finite and monotone, unlike dual
  pruning (which cycled).  Validate at 30k rows: same optimum as the full LP, less wall time.
- Combine the winners (sifting + warm simplex, or sifting + PDLP).

**Step 2: drop-in.** A `BRANCH_SOLVER=` backend in `branch.py` (and `tighten.py` if trivial) with the
best method; the LP value on the benchmark instance must agree with IPM to `1e-7` relative.

**Compute.** 8–16 pinned cores, up to ~50 GB.  Nothing else heavy should run at the same time;
`tasks/m2minus4-family` waits for this.

**Done when.** `search/LPSPEED.md` with the table (instance, method, seconds, objective, iterations)
and the backend merged.  Do not touch `TODO.md` or other tasks' files.

## Status (2026-08-30)

Done; results in `search/LPSPEED.md`.  Summary: on the real `1110` leaf LP (66k rows × 13.9k
columns) nothing that re-solves the same LP beats HiGHS IPM (1,658 s single-threaded — IPX uses one
core here): cold/warm/incremental dual simplex are 40–150k-pivot degenerate walks per solve
(the old `warm` backend was slow for that reason, not a wrapper bug — verified by a 0-iteration
zero-change test), row sifting reaches the same optimum at parity, PDLP stops four digits short.
The lever is the column count: IPX's per-iteration cost is `∝ cols^1.5` and ~90 % of the master's
columns are dead.  `BRANCH_SOLVER=restricted` (column-sifted master, `branch.py`) is the drop-in;
`--dump-lp`, `--cols-raw` (asymmetric restarts were 2–8× too wide), `search/lp_bench.py`,
`lp_dump_lattice.py`, `lp_subset.py`, `lp_toytest.py` are the tooling.  Loop-hygiene defect noted
(not fixed): the lattice warm start contains no corner-box rows, so round 0 of any restart is wasted.
