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

## Not done yet (tasks/explore-loop/DESIGN.md §8)

- Reward v2: the pool potential with lineage credit (a construction arm is paid when a later walk descends from it).
- Closures as store hooks (drop-one and exact L-add derived packings queued when a new best enters).
- Paring: `mcmin.py` (file I/O and moves into `pk/`), then the slp-era engine it imports (`slp*`, `inc`, `rigid`,
  `movegen`, `ftmc`) and the Rust `packer` binary modes; `explore.py` into a `walk` stage; the evaluation modules
  (`pk/battery`, `corpus`, `eventbench`, `benchkit`, `crossing`, `refset`, ...) into one; `trial.py` folded into goals.
