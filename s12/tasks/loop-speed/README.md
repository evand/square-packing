# K. Loop speed: stop paying for exactness every round

**Why.** Measured on the live runs (`runs/branch_J16.log`, `runs/branch_J12.log`, cumulative `t=`
and the per-phase brackets): J16 round 9 took **3.8 h — verifier 3.0 h, LP 40 min, clique
separation 47 s**; J12 round 9 took **4.8 h — verifier 4.7 h, LP 55 s**.  The exact verifier is
80–99 % of every round and its cost doubled between rounds 8 and 9 (more clique pieces, more atoms,
contention).  Task C already fixed the LP (`BRANCH_SOLVER=restricted`).  A leaf is now 40–80 rounds
× hours; the per-slot tree of `notes/level2-design.md` is 15 leaves.  The exact sweep is being used
as the *separation oracle* — to tell the LP which poses it under-covers — a job a float scan does in
seconds; exactness matters once, for the final file.  Target: **≥ 10× per round**, and fewer rounds.

Read first: `search/LPSPEED.md` (what was measured and the loop-hygiene notes), `search/BRANCH.md`
("Compute" bullet: the sweep is quadratic in the atoms), `search/ANCHOR.md` §"Two things worth
recording" (the anchor-credit prefilter work already done), `search/branch.py` (the loop: rows from
verifier witnesses, `--topk`, `--prune-at`, `--row-grid`, restricted master, `--cq-interior`,
`--matched`), `search/tighten.py`, `verify/src/main.rs` (the sweep: per angle bin, arrangement of
atom squares, cell credit incl. anchor pieces; witness mode `topk`), `xcheck.py`,
`tests/rejection_tests.sh`, `tests/bitid.sh` (bit-identity harness against a reference binary),
`search/closed4.py` and `search/rung2_close.py` (they already have a dense float scan + polish
separation — reuse), `search/zeromargin.py --oracle` (another separation mode), `search/RECONCILE.md`
§"State of the runs" and `runs/J12.sh`, `runs/J16.sh` (the exact commands and checkpoints).

## Deliverables, in order

1. **Instrument first.**  Per round: wall and CPU seconds per phase, number of verifier calls, atoms,
   cells swept, clique pieces credited, rows added, rows violated by the float scan vs by the exact
   sweep.  A one-line summary per round in the log (extend the existing `t=… [lp … ver …]` bracket).
   Benchmarks on this loaded box are fine for ≥ 10× claims — report load average and what else ran
   next to every timing; a clean re-run can come later.  Use the checkpoints on disk as instances:
   `/home/evand/math/square-packing/s12/runs/branch_J16_{cols,cliques,probe}.txt`,
   `branch_t398ik4n_*`, `branch_t398hk4_*`, `branch_t398j1110_*_it51.txt`; **do not touch or restart
   the running J12/J16 processes** (they run from another worktree; leave them alone).

2. **Float separation oracle for the loop** (`branch.py`, `tighten.py`).  Each round: (a) a dense
   pose scan (centre pitch ~0.003–0.005, angle pitch ~0.25–0.5°, over the full admissible domain, or
   the D4 fundamental domain for symmetric certificates) evaluating captured weight of points *and*
   clique pieces in numpy; (b) local descent (Nelder–Mead / coordinate) from the worst grid poses and
   from last round's witnesses; (c) emit the `topk` worst per region as rows, same row keys and
   dedup as today (`--row-grid`).  Clique-piece membership in floats must be **conservative in the
   same direction as the verifier** (a pose near a piece boundary counts as *outside*), so the float
   oracle never claims coverage the exact sweep would deny.  Then: exact verifier only every `E`
   rounds (`--exact-every`, default 5) and at convergence (`probe_min >= 1` on the float oracle for
   two consecutive rounds → exact verify at `N = 2000`; if it fails, its witnesses become rows and
   the loop continues).  Nothing is *claimed* until Rust + `xcheck.py` accept the file — say so in
   the write-up and keep the final step unchanged.  Validate: from the same checkpoint, the float
   loop and the exact loop reach the same LP value ±0.002 and the same verified file passes.

3. **Bucket the exact sweep** (`verify/src/main.rs`, then `xcheck.py` if it shares the structure).
   Two atoms interact only if their unit squares can meet, i.e. centres within distance √2; a
   spatial grid (cell ~1) makes the per-bin arrangement near-linear in the atoms instead of
   quadratic.  Same for anchor-piece credit: pre-bucket pieces by the cells their `contains`/`meets`
   predicates can possibly hold in.  **Guardrail: `tests/bitid.sh` must report ALL IDENTICAL**
   (stdout, witness files, `TIGHT_DUMP`) against the pre-change binary on every reference
   certificate, and the 136 rejection tests must pass; `xcheck.py --all` agreement on the demo
   certificates.  Report speed-up per atom count (3.4k, 10k, 13k atoms; with and without cliques;
   `N = 2000` and `6000`).

4. **Fewer rounds.**  (a) `--seed-from LEAF`: start a leaf from a sibling's converged rows,
   columns and cliques (the per-slot leaves differ only in multipliers; today each starts from the
   lattice warm start and wastes round 0 with λ at its cap and no corner rows — fix that too by
   seeding the corner-box rows into the warm start).  (b) In–out / proximal stabilisation of the
   row step (add rows at a convex combination of the LP point and the previous stabilised point;
   or a small quadratic proximal term on the weights) to damp the probe sawtooth (`0.29 → 0.87 →
   0.46 → 0.94 …`).  Measure rounds-to-converge on the `k = 4` leaf from `branch_t398hk4_*` with and
   without.

5. **Write-up** `search/LOOPSPEED.md`: the phase table before/after, the speed-ups with their
   instance sizes and load, the validation (same LP value, same verified file), and what remains
   slow.  Update `README.md`'s "Reproducing" section only if flags changed — otherwise leave docs
   to the coordinator.

## Compute

≤ 8 threads, ≤ 40 GB.  J12 and J16 (6 verifier threads each) must keep running undisturbed; check
`uptime` and pin with `taskset` (cores `n`, `n+16` are hyperthread siblings).  `runs/` is
gitignored: read `/home/evand/math/square-packing/s12/runs/` by absolute path, write your own.
**First check your worktree base** (`git merge-base branch-corners HEAD` vs `git rev-parse
branch-corners`; `git merge branch-corners` if they differ).

## Done when

2 and 3 merged behind flags with defaults that keep today's behaviour unless opted in (`--sep float`,
`--exact-every`), `tests/bitid.sh` ALL IDENTICAL, rejection tests green, `verify.sh` green, and
`search/LOOPSPEED.md` with the numbers.  4 as far as time allows, measured.  Commit in your
worktree; do not touch `TODO.md` or other `tasks/*` folders.
