# packer architecture (rewritten 10-24 for the explore loop; earlier versions in git history)

**A goal is data; one scheduler spends CPU on arms; every arm is a stage chain run on donors selected from the store;
every output goes back to the store; reward is measured from the store.**  A new idea is a new stage (or a parameter of
one).  A new campaign is a new goal file.  There are no campaign scripts.

## Layers

| layer | where | what |
|---|---|---|
| store | `pk/store.py`, `pk.py` | every packing (dedupe, provenance, parent links, status).  Read with `pk.py get` / queries; write through `pk/sink.py`.  New near-record non-45-degree packings are re-quenched on entry (`pk/converge.py`). |
| engine | `pk/engine.py` | `fq_quench(inp, QOpt)`: the one fq call path; `QOpt` = every fq flag |
| registries | `pk/moves.py` (`@move`), `pk/fill.py` (`@placer`), `pk/schedules.py` (named QOpts) | building blocks with declared parameter spaces |
| stages, chains | `pk/stages.py` | one contract over all registries; a chain is `a(k=v) \| b \| ...` |
| goals, arms, scheduler | `pk/goal.py`, `pk.py goal`, `<plugin>/goals/*.yaml` | the loop below |
| jobs | `pk/jobs.py` | the only way compute starts; on exit every job runs sink, then converge |
| reports | `pk.py goal report`, `pk.py yield`, the plugin's status report | numbers from the store |

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

## The loop (`pk/goal.py`)

```
goal YAML   targets, progress (wanted families, window), value weights, budget, quantum, predictions (required)
plan        arms per target n: walk (donors at n), add (n - 1), drop (n + 1, n + 2), lcut (n + 2k + 1 + e), and each
            plugin's ARM_TEMPLATES; an arm = (donor offset, chain, donor family filter)
select      donors = wanted-family packings within `donor_window` of the best wanted-family packing at the donor n
            (family-relative near-record: at 148 the best non-diamond is 0.06 above the diamond)
run(job)    chain arms build starts for `construct_share` of the quantum, then walk the best four; walk arms walk the
            donors.  `walk` = explore.py's local search with a fixed move mix
finish      inside the job, on its own CPU quantum (the walk leaves a measured share for this): ingest (starts with
            their donor as parent, the walk archive, walk roots linked to their donor), converge the job's packings,
            fq-polish its candidates best first until the quantum is spent; job.json has the CPU split
count       in the scheduler (cheap, serial): the job's polished packings as new region clusters within the window of
            the polished family best at goal start, plus the family-best and store-best gains (side units x weight).
            The reference minima are polished once per n before any job starts; their CPU counts in the budget
allocate    Thompson sampling on reward per CPU-h (Gamma-Poisson, pooled prior), pending jobs counted, `procs` at a time
```

Run one: `python3 pk/jobs.py start NAME --procs P -- python3 pk.py goal run NAME`; read it with `pk.py goal report
NAME`; `pk.py goal plan NAME` lists arms and donor pools without running anything.  Goal files live in a plugin's
`goals/` (they name plugin arms and hold unposted results); this repo keeps only `goals/example.yaml`.  The goal file is the
pre-registration (the scheduler refuses one without predictions); write the reading of the result under it.

## Adding an idea

1. Register it as a stage of the right kind, where its family lives: a placer in `pk/fill.py`, a move in `pk/moves.py`,
   an fq recipe in `pk/schedules.py`, anything else in `pk/stages.py` (or a plugin).  Declare every knob, return an info
   dict, draw randomness only from `rng`.
2. `python3 pk/tests/test_stages.py`.
3. If it changes n or builds starts, it becomes an arm automatically when a goal's planner (or a plugin's
   `ARM_TEMPLATES`) names its chain; compare it with the existing arms by running a goal and reading the per-arm
   reward per CPU-h.  Within-n refinements are judged the same way: an arm earns reward only through refined, converged
   packings in the store.
4. Steering results go into `CLAIMS.md` (one line, with evidence).

## Plugins

`pk/stages.py` imports `pk_plugin.py` from each directory on `PK_PLUGINS` (default `../plugins/packer` next to this
repo, if present) before adapting the registries.  A plugin registers stages with the same decorators and may define
`ARM_TEMPLATES` (`f(goal, n) -> [arm dict]`), `TEST_STARTS` (`{stage: chain}` for the contract test) and goal files of
its own.  Plugin scripts find this directory through their own bootstrap; `runs/` and the store live here.

This side stays generic (Evan 10-25): extension points, public stages, and examples built from public stages.  A
plugin's stage / arm names, its weights and priors, its benchmarks and its goal files live in the plugin, never here;
the pre-push hook enforces the names (`CLAUDE.md`, Public / private).

## Next: the unit loop (design 10-25, agreed with Evan; steps 1-2 built 10-26: `pk/families.py`, `pk/units.py`)

Why.  The goal loop schedules 15-minute single-arm jobs from one dispatcher process per goal.  Measured 10-25 on two
30 CPU-h runs: utilization 48-79 % of the declared processes (serial setup, a 24-minute drain of the last wave,
hand-off to the next queued run, a queue that would not start a 15-process run behind one stray thread); references
31-43 % of CPU; sequential runs as experiments (different targets, a moving store: +-35 % noise on gain per CPU-h, so
~120 CPU-h per condition to see a 2x effect).  Everything here is embarrassingly parallel except decisions and bandit
updates, and those are small database transactions.  The whole project is a value-of-information problem: where to
spend the next CPU-minute (which n, which move, polish or explore, which variant to learn about).  Overhead decides
how many such decisions we can afford, so the rule is: **make decisions cheap**, and make every unit of compute a
decision.

Principles.
- The unit is one single-threaded task with a fixed CPU allowance (~1-3 min) that reads the store, computes and writes
  its results.  No job, quantum, dispatcher or per-goal process.
- Decide late: a unit picks what to do when it starts, from everything finished so far.
- Decisions read precomputed state (indexes, summary rows): milliseconds, never a scan.
- Everything is a row: goals, units, results, CPU.  Accounting, reports and experiment readouts are SELECTs.

Unit kinds (each a stage chain or an existing entry point, run under a CPU allowance; fq `--deadline` keeps the tails):
- `construct`: donors -> build stages -> screen R builds -> fully relax the best 1 (+ an audit share, random).
- `walk`: the local search from a chosen store packing for the allowance.
- `polish` / `converge`: one packing (or the new near-best packings at one n) to a fixpoint.  Queued implicitly: they
  are options whose value is the chance of a gain, like any move ("relax and polish are moves").
- `reference`: the polished family best at n before the goal (only what gain accounting needs; cluster census off).

Claim (one transaction, any worker):
```
BEGIN IMMEDIATE
  pick goal  (priority, budget left, CPU share)
  pick option (n, kind, arm, variant) by Thompson sampling on expected gain per CPU-s, from the summary tables
  INSERT unit (goal, n, kind, arm, variant, worker, lease_until, started)
COMMIT
run -> INSERT result rows (packings via the store API); UPDATE unit (cpu, wall, gain, status); UPDATE summaries
```
A lease that expires (crashed or killed worker) returns the option; nothing is lost but that unit's CPU.

Workers.  `pk.py work --procs N`: N persistent processes (imports and caches once, then claim / run / record).  The
count follows the machine (free cores from `jobs.py`; quick tests for an active conversation still go first), so two
goals or another session's work share cores at unit granularity; a run "drains" in one unit.

Decision state (kept current on write, read on claim):
- store: `family` column on packing, a `best` table (n, family -> best side, id) maintained on insert, index (n, s);
  donor choice is a lookup.
- targets: live value per n = (current side - bar) / sqrt n from the newest bar list, x expected gain per CPU-s;
  n retire when they cross the bar.  Goals state an objective and constraints, not a hand-picked list of n.
- model: gain per CPU-s per option with pooled priors (kind x n band, kind, n) and discounting (gains deplete as a
  target is worked: never-worked targets paid 2-4x per CPU-h on 10-25); credit for a polish / converge gain goes to the
  unit that made its parent (lineage).

Experiments.  A goal may declare variants (`clusters: [on, off]`, `screen.k: [2, 8]`, a polish policy, ...); each claim
draws one at random.  Variants then share targets, time and store state, so the comparison is paired within one run
and the readout is a grouped SELECT with a bootstrap.  Pre-registration stays: the goal file holds the predictions.

Accounting.  Unit rows carry CPU (`getrusage` of the worker around the unit, children included), wall and kind;
worker rows carry idle time.  Utilization = unit CPU / (wall x workers); overhead = claim + record time, reported per
unit.  Targets: utilization >= 95 %, decision + record overhead <= 1 % of unit CPU.

Migration.  1. `best` table, `family` column, index (cheap decisions; the current loop benefits too).  2. Units and
workers next to the current loop, running the same arms; smoke goal.  3. One fresh-target run with the 10-25 lean
settings (clusters off) for a like-for-like readout of utilization and gain per CPU-h.  4. Variants, then live
targets.  5. Retire `goal run` / `goal job` / the quantum knobs.
Step 2 as built (10-26): `python3 -m pk.units add GOAL` / `work --procs N` / `report GOAL`; kinds arm, reference (only
when a candidate at n needs it), polish (candidates = an arm unit's new wanted-family packings within `polish_margin`
1e-3 of the family best: polish moved 97 % of 7815 10-25 goal-run candidates by < 1e-3); a polish unit's CPU and gain go to the
arm that made its candidates.  Smoke (2 targets, 0.8 CPU-h, 3 workers beside a 14-process run): claim + record 0.6 ms
median (0.01 % of unit CPU), units' CPU = the jobs.py measurement.  Workers swap in the cached-family donor lookup for
their own process until goal.py adopts it.

Open.  Unit allowance per kind (screened construct needs ~R x 0.5 s + one relax; polish at n ~ 1000 is 1.5-3 min);
SQLite write rate (16 workers x a few writes per unit per minute: WAL should do); the cluster census if it returns
(needs serializing); how much of the walk's value is in long chains (cut into units, it restarts from the best so far).

## Not done yet (tasks/explore-loop/DESIGN.md §8)

- Reward v2: the pool potential with lineage credit (a construction arm is paid when a later walk descends from it).
- Closures as store hooks (drop-one and exact L-add derived packings queued when a new best enters).
- Paring: `mcmin.py` (file I/O and moves into `pk/`), then the slp-era engine it imports (`slp*`, `inc`, `rigid`,
  `movegen`, `ftmc`) and the Rust `packer` binary modes; `explore.py` into a `walk` stage; the evaluation modules
  (`pk/battery`, `corpus`, `eventbench`, `benchkit`, `crossing`, `refset`, ...) into one; `trial.py` folded into goals.
