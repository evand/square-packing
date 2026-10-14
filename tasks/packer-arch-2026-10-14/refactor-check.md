# Bit-exact refactor check: one fq call path + shared neighbour transforms (2026-10-14)

Change class: bit-exact (README "Change classes").  Base 8e895d7; old = `git archive 8e895d7 search/packer`, new = this
commit; both trees ran the same fq / harden binaries (main checkout `target/release`, fq behind an argv-logging wrapper).

## What changed
- `pk/engine.py`: `fq_quench(inp, QOpt, mode='subprocess'|'serve', out=, tmp=, binary=, timeout=, mod90=)` is the one
  fq quench call; `QOpt.args()` the one flag builder; `QOpt` gained the remaining fq quench flags (round_kick, round_p,
  round_loosen, soft_p, d0, dmin, vtol, tr_grow, tr_shrink, sqp, finish, ident, v, dump_lp; appended, so existing
  QOpts give the same args); `QOpt.from_args` / `with_args` parse flag lists (last wins, as fq's `arg`); unknown flags
  raise.  `FQ.request` = raw serve request (FQ.quench uses it, unchanged results).  `write_cfg` / `read_cfg`;
  preset `POLISH_DEEP` (= polish, flip-top 64, pit 600).
- Call sites on it: `hop.quench` (legacy `EXTRA` + `extra` parsed into a QOpt, or `opt=`), explore (screen
  `SCREEN.but(loosen)`, polish `POLISH.with_args(--polish-extra)`, kick retry, start / final polish), chain seed_job,
  compose2.polish, harden.job, relax.job (`OPS` are QOpts), rectlib (`_fq`, compose polish), fresh (`repolish`).
  All production calls stay `mode='subprocess'`.
- `pk/moves.py`: `remove_k`, `insert_k(holes='pose'|'uniform', dilate=)`, `harden_xform`; chain.seed_job and
  harden.make_input call them.  Each branch consumes its rng (random.Random / numpy Generator) exactly as before.

## Results (scripts in `refactor-check/`; inputs: store best packings at n = 40, 41, 110, 111, 237, 238)
- `drive.py` old vs new, `compare.py`: 47 call-site results (hop.quench x 4 flag sets incl. a duplicate `--loosen`;
  explore._screen_polish incl. kick_retry with sigma 1.234567e-4 and `--stag-tol 1e-7` polish-extra; compose2.polish;
  harden.job; 7 relax ops; rectlib._fq with an np.float64 d; fresh re-polish; rectlib compose polish; 6 chain seed_job
  (k = 1 at 41->40 x 2, 238->237; k = -1, -2 insertions), TRIES = 4) and 84 harden make_input cases (none / insert:1,3 /
  remove:all / remove:R x 4 seeds, n = 40, 110, 237): **0 of 103 output / input files differ (byte compare); 0 of 47
  results differ** (floats compared by bits; only fq's wall-clock fields `sec`, `t_alm`, trace step time masked).
- fq argv: 95 calls; 62 identical, 33 differ only in flag order / number spelling / a dropped overridden duplicate
  (e.g. `--seed 11 --loosen 1.0 --no-alm` -> `--loosen 1.0 --seed 11 --no-alm`; `--stag-tol 1e-7` -> `1e-07`;
  `--loosen 1.02 --no-alm --loosen 1.0` -> `--loosen 1.0 --no-alm`); all 33 equal after fq's own parse (rposition, last
  wins; values as f64 / int).  Argv is canonical QOpt order now, so raw argv identity across sites is not kept; fqreplay's
  `--pit 8` / `--no-alm` / `--kick` categorisation and its `--out`-offset parse still hold.
- `serve_check.py`: fq_quench subprocess vs serve, 60 (input, option) pairs (12 option sets incl. kick, round / soft
  pre-stages, shrink, rect; file and in-memory inputs with angles outside [0, 90)), serve run twice in one worker in
  two request orders: **0 differ** (output file bytes + JSON, timing masked).  No state leaks across serve requests.
- `engine_check.py`: FQ.quench old vs new, 9 cases: identical.
- `pk/tests/test_packing.py`, `pk/tests/test_watch.py`: ok.

Not covered by a run: fresh / rectlib `run` drivers end to end (their fq calls are tested through the same helpers),
the many analysis scripts that still call `fq` directly (polish_trace, polish_value, round_ctl, removal_ctl,
bench_quench, anneal_test, fqreplay, chainreplay: left as they are).
