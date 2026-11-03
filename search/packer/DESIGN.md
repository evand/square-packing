# packer architecture (rewritten 10-27 for the event loop; the goal / unit loops of 10-24..26 are in git history)

**A new packing is an event; every event becomes a parent; workers claim (parent, move) options one at a time, ranked
by a scorer; every move is a stage chain; every application goes to the store's move log and every useful output back
to the store as a new parent.**  A new idea is a new stage (or a parameter of one), used as a move.  A new campaign is
a new campaign file.  There are no campaign scripts.

## Layers

| layer | where | what |
|---|---|---|
| store | `pk/store.py`, `pk.py` | every packing (dedupe, provenance, parent links, status; every driver records parents: `parents.json`, `_p<id>_` file names, the move log, or `pk.py import --parent`; `pk.py lineage` reads them all).  Read with `pk.py get` / queries; write through `pk/sink.py`.  New near-record non-45-degree packings are re-quenched on entry (`pk/converge.py`). |
| engine | `pk/engine.py` | `fq_quench(inp, QOpt)`: the one fq call path; `QOpt` = every fq flag |
| registries | `pk/moves.py` (`@move`), `pk/fill.py` (`@placer`), `pk/schedules.py` (named QOpts) | building blocks with declared parameter spaces |
| stages, chains | `pk/stages.py` | one contract over all registries; a chain is `a(k=v) \| b \| ...` |
| events, scorers | `pk/events.py`, `<plugin>/goals/*.yaml` | the loop below |
| pool | `pk/tasks.py` | short single-threaded units by priority: the default way compute starts |
| jobs | `pk/jobs.py` | long runs that cannot be split; the pool counts them as outside load |
| follow-ups | `pk/react.py` | new store rows -> pool units (converge, near misses, plugin rules); `hold` for isolated experiments |
| reports | `python3 -m pk.events report NAME`, `pk.py yield`, the plugin's status report | numbers from the store |


## Compute (11-01)

Three ways work runs; choose by shape, never by core count:

| Work | How | Notes |
|---|---|---|
| Short, independent (10-60 s, a few minutes at most) | `pk/tasks.py add` / `repeat`, priority high, normal or low | No preemption: a slot frees about once a second under real units |
| A script keeping state between tries | ordinary process (`jobs.py start` for a log) with `tasks.Executor` | The script sleeps; each call is a unit; its memory is lost if it dies, saved results are not |
| Long and unsplittable, renders, anything else | `jobs.py start`, or just run it | Counts as outside load; the pool makes room within about 30 s |

Event-loop campaigns: `pk/tasks.py repeat NAME --budget CPU_H --prio low -- python3 -m pk.events work --campaign NAME --max-units 3`.
`--campaign` matters: without it a worker claims from every active campaign, and pool priority stops meaning anything.

Parallel scripts use `tasks.Executor` (submit / map like ProcessPoolExecutor; `initializer=`, per-unit `timeout=`,
`map(chunksize=)` for calls much shorter than 10 s; `jobpool.run_jobs` is the same with a hard timeout).  No private
pools (`--procs`, `--local`, `Pool(N)`): a fixed pool can't grow when neighbours end (11-02: hv853 at 3 procs and
tiltx `--local 4` left about 8 cores idle for hours).

Follow-ups: `pk/react.py` reads new store rows once a minute and queues units from its rules (`@rule`; plugins add theirs).
Results feed back as new rows.  `python3 -m pk.react hold` stops follow-ups for isolated experiments; rows arriving during
a hold never get them.  `pk/tasks.py limit 0` drains the pool for timing runs.  A paused campaign (`FAIL_STOP` failures in a row) fails its `tasks.Executor` calls, so the script stops and says why; a new `Executor` on that name clears the pause, as does `cancel`; plain `add` warns and its units wait for `resume`.  Tools that write result files instead of
store rows need `pk.py sink` afterwards: nothing imports them automatically now.

Measured 11-01 on event-loop runs of 1.5 CPU-h: 15.85 cores in units (long-lived workers 16.00), CPU overhead 0.41 s per
process; feedback from a new row to its converge child, median 1.0 min, p90 1.7 (the former job exit step: 32 min, 7 h).
Open: whether the faster loop finds better packings is untested, and out of scope of that work.

## Stages

A stage is `f(st: State, rng, ctx, **params) -> (State, info)` registered with a kind and a parameter space
(unspecified parameters are drawn; every value used is logged).  `State`: the packing, a per-square `frozen` mask, the
`region` the next recreate fills, generator facts in `meta`.

| kind | does | examples |
|---|---|---|
| `gen` | a start from nothing | `gen.file` |
| `ruin` | open a region | (plugin) |
| `recreate` | fill the region up to ctx n | `recreate.<placer>` for every `pk/fill.py` placer |
| `constrain` | freeze / release | `constrain.freeze`, `constrain.release` |
| `relax` | a T = 0 descent | `relax.<schedule>` for every `pk/schedules.py` entry; default after constructions: `round.15c.p5` (CLAIMS C50) |
| `perturb` | a stochastic move at fixed n | `perturb.<move>` for every `pk/moves.py` move |
| `grow` | n upward, unrelaxed | `grow.add`; plugin L constructions |
| `shrink` | n downward, unrelaxed | `shrink.drop`, `shrink.lcut` |

n-changing and perturb stages refuse frozen squares.  `python3 pk/tests/test_stages.py` checks every registered stage
(parameters logged, deterministic per seed, frozen squares untouched); a new stage is covered by registering it.

## The loop (`pk/events.py`, 10-27)

```
campaign YAML  targets (ns), w(n), moves (name -> chain, dn, once / relax / source), policies (A/B conditions),
               window, excluded parent families, unit CPU, budget, predictions (required); per move conds (A/B on
               move sets), admit; admit_window (admission width drawn per unit)
pool           per condition: the best distinct store packings within the window at every parent n, then every output
               of the condition's own units (a condition never sees another's outputs; frozen conditions add none)
               (`parents: [store ids]`: those packings are the start pools instead, e.g. outputs of one move, to test
               what follows them without the scorer's selection; they still need to be within the window)
claim          any worker: sources first (generator builds, refreshes), then the head of the condition's claim queue.
               A refill (one worker at a time; scoring outside the write lock) queues relax options first (a parent
               not yet relaxed offers only relax), else a batch of options (parent, move, target n) by sampled
               utility x w(n), round-robin over several Thompson draws, one per parent; a new best at n drops that
               n's entries and triggers a full refill (coalesced, >= batch_gap s apart; batch_age a backstop).  Targets must pass the grid-only rule.  (10-29: per-claim rescoring of every option cost
               2.8 s at 260k options; near-best options number in the thousands, so batches lose nothing.)
run            the chain on the parent until the unit's CPU is spent; each application -> a move-log row (`expansion`,
               child = the stored output); outputs below the parent (same n) or within the window (cross n) -> store
record         tries / wins per (parent, move), new parents, the condition's best per n (its score gain)
```

Unit structure (CLAIMS C54, untested): a chain is one structural step plus the minimum repair that makes its output
store-ready ([fill] + relax); everything else composes through the store, where the scorer allocates and selection acts.

Grid-only rule (Evan 10-27): in triangle row k (k = ceil sqrt n, c = k^2 - n) only the first grid-only entry is a
target, and only once row k - 1 has that c; deeper entries are never easier (monotonicity).  Scorers (`events.SCORERS`,
by the condition name before '-'): `uniform`, `rule` (near the best, few tries; regular families are never parents),
`tilted` (11-02: proportional sampling over (1 + A_same)^-0.5, A_same = discounted same-kind claims on the parent's
lineage; every queued entry logs its propensity `_w / _den`, so offline analyses can reweight; `draws: 1`),
and plugin models fit on the move log.  A scorer returns sampled utilities, so exploration lives in it (Thompson draws
for a fitted model).  `ctx['cpu']` gives the mean CPU per application of each move so far.

Equivalent choices are mixed, not settled (Evan 10-27): when options do not differ by enough to matter for the
question at hand (scorers, schedules, placers, parameter values), draw among them at random (ranges in chains, a mixed
condition 'a+b'); every application stays a logged experiment, so the knowledge keeps building while attention is on
the larger questions.

Run one: `python3 -m pk.events add NAME`, then the pool `repeat` line under "Compute"; read it with `python3 -m pk.events report NAME`.  Campaign files live in a plugin's `goals/` (they name
plugin stages and hold unposted results); this repo keeps only `goals/example.yaml`.  The campaign file is the
pre-registration (the loader refuses one without predictions); write the reading of the result under it.

## Adding an idea (11-02)

A new move is a lab stage plus `pk.py try`:

1. Register it as a stage of the right kind with `status='lab'` and an `evidence` pointer (a plugin's `lab/*.py` is
   imported automatically; or next to its family: `pk/fill.py`, `pk/moves.py`, `pk/schedules.py`, `pk/stages.py`).
   Declare every knob, return an info dict, draw randomness only from `rng`.  `python3 pk/tests/test_stages.py` covers it.
2. `pk.py try 'new.stage | relax.round.15c.p5' --ns 791,853 --tries 8`: every arm and the control (`relax.round.15c.p5`
   alone) on the same parents and seeds, every application a pool unit (prio high) and a move-log row
   (`run = try:NAME`), children in the store with parent links; `pk.py try --report NAME` gives the paired reading.
   `--admit window` keeps score-neutral outputs too (within the window of the best at n), for moves that pay off later.
   `--hold` keeps follow-ups (`pk/react.py`) off this experiment's outputs.
   An n-changing arm is `'dn=D: CHAIN'` (parent at n0, ctx n = n0 + D); `@NAME` is a named move (a plugin's `MOVES`).
   Each row is scored against the store best at its target n when the plan is made (`--asof ISO`: among rows added
   before then, to replay a past win); wins need `s < ref * (1 - events.WIN_RTOL)`.  The paired diff is in score units
   (a failed or worse output counts as keeping the parent), so arms of different dn compare on one scale; a
   deterministic control (`stages.deterministic`) runs once per parent.  `--requires F,...` keeps parents whose cached
   features (`pk/features.py`) qualify; campaign moves take `requires:` too.  `constrain.best` continues a chain from
   its best packing so far (keep-if-better cycles in one application).
3. If it earns a place, name its chain as a move in a campaign (`admit: window` per move there too), set
   `status='default'` and the evidence; steering results go into `CLAIMS.md`.  `pk.py stages [--md] [--status lab]`
   is the move registry: registry fields joined with the move log's use of every stage.

Scripts that are not stages log through one call, `Store.record(parent, move, params, child, cpu_s, run=)` (or
`pk.record`): a row for every output, kept or not.  The scratchpad is for plots and one-off analysis.  If you are writing
a loop over parents in a scratch script, stop: that is `pk try`.

## Plugins

`pk/stages.py` imports `pk_plugin.py` from each directory on `PK_PLUGINS` (default `../plugins/packer` next to this
repo, if present) before adapting the registries.  A plugin registers stages with the same decorators and may define
`TEST_STARTS` (`{stage: chain}` or `(chain, start n, ctx n)` for the contract test), features (`pk.features.feature`), event-loop scorers (`pk.events.SCORERS[name] = f`) and
campaign files of its own.  Plugin scripts find this directory through their own bootstrap; `runs/` and the store
live here.

This side stays generic (Evan 10-25): extension points, public stages, and examples built from public stages.  A
plugin's stage / move names, its scorers' coefficients and priors, its benchmarks and its goal files live in the plugin, never here;
the pre-push hook enforces the names (`CLAUDE.md`, Public / private).

## Not done yet (tasks/explore-loop/DESIGN.md §8)

- Event loop: lineage credit for sources (a build is paid when a later move descends from it); closure credit (a
  packing at (n, s) is also worth its drop closure at every n' < n); walks (`explore.py`) and screened constructions
  (`pk/screen.py`) as sources; store imports (public drops) as events; above-grid parents / children at grid-only
  targets (today pools and children stop at s < k; Evan 10-27: later); a per-campaign intercept learned online if
  decisions across moves / sources need calibrated totals (scorer totals were off x2-6 per run; ranking within a run
  is unaffected); claim cost ~40 ms per pop plus one ~4 s refill per batch or per new condition best (coalesced, >= 10 s apart) at ~260k options (10-29); partial refills (rescore only the n a new best touches, under the batch's stored draws so utilities stay comparable) when the pool grows ~5x;
  a near-best diversity reward when there is a cheap measure (first try: the child's family / fingerprint is new in
  its condition's pool at n).
- Closures as store hooks (drop-one and exact L-add derived packings queued when a new best enters).
- Paring: `mcmin.py` (file I/O and moves into `pk/`), then the slp-era engine it imports (`slp*`, `inc`, `rigid`,
  `movegen`, `ftmc`) and the Rust `packer` binary modes; `explore.py` into a `walk` stage; the evaluation modules
  (`pk/battery`, `corpus`, `eventbench`, `benchkit`, `crossing`, `refset`, ...) into one; `trial.py` folded into campaigns.
