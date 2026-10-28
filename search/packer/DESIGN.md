# packer architecture (rewritten 10-27 for the event loop; the goal / unit loops of 10-24..26 are in git history)

**A new packing is an event; every event becomes a parent; workers claim (parent, move) options one at a time, ranked
by a scorer; every move is a stage chain; every application goes to the store's move log and every useful output back
to the store as a new parent.**  A new idea is a new stage (or a parameter of one), used as a move.  A new campaign is
a new campaign file.  There are no campaign scripts.

## Layers

| layer | where | what |
|---|---|---|
| store | `pk/store.py`, `pk.py` | every packing (dedupe, provenance, parent links, status).  Read with `pk.py get` / queries; write through `pk/sink.py`.  New near-record non-45-degree packings are re-quenched on entry (`pk/converge.py`). |
| engine | `pk/engine.py` | `fq_quench(inp, QOpt)`: the one fq call path; `QOpt` = every fq flag |
| registries | `pk/moves.py` (`@move`), `pk/fill.py` (`@placer`), `pk/schedules.py` (named QOpts) | building blocks with declared parameter spaces |
| stages, chains | `pk/stages.py` | one contract over all registries; a chain is `a(k=v) \| b \| ...` |
| events, scorers | `pk/events.py`, `<plugin>/goals/*.yaml` | the loop below |
| jobs | `pk/jobs.py` | the only way compute starts; on exit every job runs sink, then converge |
| reports | `python3 -m pk.events report NAME`, `pk.py yield`, the plugin's status report | numbers from the store |

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
               window, excluded parent families, unit CPU, budget, predictions (required)
pool           per condition: the best distinct store packings within the window at every parent n, then every output
               of the condition's own units (a condition never sees another's outputs; frozen conditions add none)
claim          one transaction, any worker: sources first (generator builds, refreshes), then relax-first (a parent
               not yet relaxed offers only relax), else the option (parent, move, target n) with the largest sampled
               utility x w(n) from the condition's scorer; targets must pass the grid-only rule
run            the chain on the parent until the unit's CPU is spent; each application -> a move-log row (`expansion`,
               child = the stored output); outputs below the parent (same n) or within the window (cross n) -> store
record         tries / wins per (parent, move), new parents, the condition's best per n (its score gain)
```

Grid-only rule (Evan 10-27): in triangle row k (k = ceil sqrt n, c = k^2 - n) only the first grid-only entry is a
target, and only once row k - 1 has that c; deeper entries are never easier (monotonicity).  Scorers (`events.SCORERS`,
by the condition name before '-'): `uniform`, `rule` (near the best, few tries; regular families are never parents),
and plugin models fit on the move log.  A scorer returns sampled utilities, so exploration lives in it (Thompson draws
for a fitted model).  `ctx['cpu']` gives the mean CPU per application of each move so far.

Equivalent choices are mixed, not settled (Evan 10-27): when options do not differ by enough to matter for the
question at hand (scorers, schedules, placers, parameter values), draw among them at random (ranges in chains, a mixed
condition 'a+b'); every application stays a logged experiment, so the knowledge keeps building while attention is on
the larger questions.

Run one: `python3 -m pk.events add NAME`, then `python3 pk/jobs.py start ev --procs N -- python3 -m pk.events work
--procs N`; read it with `python3 -m pk.events report NAME`.  Campaign files live in a plugin's `goals/` (they name
plugin stages and hold unposted results); this repo keeps only `goals/example.yaml`.  The campaign file is the
pre-registration (the loader refuses one without predictions); write the reading of the result under it.

## Adding an idea

1. Register it as a stage of the right kind, where its family lives: a placer in `pk/fill.py`, a move in `pk/moves.py`,
   an fq recipe in `pk/schedules.py`, anything else in `pk/stages.py` (or a plugin).  Declare every knob, return an info
   dict, draw randomness only from `rng`.
2. `python3 pk/tests/test_stages.py`.
3. Name its chain as a move in a campaign (a source if it builds from nothing, `dn` if it changes n); compare it with
   the existing moves by running the campaign and reading the per-move rows of the report and the move log.
4. Steering results go into `CLAIMS.md` (one line, with evidence).

## Plugins

`pk/stages.py` imports `pk_plugin.py` from each directory on `PK_PLUGINS` (default `../plugins/packer` next to this
repo, if present) before adapting the registries.  A plugin registers stages with the same decorators and may define
`TEST_STARTS` (`{stage: chain}` for the contract test), event-loop scorers (`pk.events.SCORERS[name] = f`) and
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
  is unaffected); claim cost ~130 ms at ~2k parents x moves x 6 settings (cache options per condition if it grows);
  a near-best diversity reward when there is a cheap measure (first try: the child's family / fingerprint is new in
  its condition's pool at n).
- Closures as store hooks (drop-one and exact L-add derived packings queued when a new best enters).
- Paring: `mcmin.py` (file I/O and moves into `pk/`), then the slp-era engine it imports (`slp*`, `inc`, `rigid`,
  `movegen`, `ftmc`) and the Rust `packer` binary modes; `explore.py` into a `walk` stage; the evaluation modules
  (`pk/battery`, `corpus`, `eventbench`, `benchkit`, `crossing`, `refset`, ...) into one; `trial.py` folded into campaigns.
