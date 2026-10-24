# packer architecture (rewritten 10-22; the 10-09 redesign plan it replaces is in this file's git history)

Every packing flows through the store; every fq call through one engine.  Every idea is a registered **stage** with a
declared parameter space, and ideas compose into **chains**.  An idea is graded by a **trial**: its chain head, then
a standard tail (release, anneal, polish), against a baseline that gets the same tail.

## Layers
| layer | where | what |
|---|---|---|
| store | `pk/store.py`, `pk.py` | every packing (dedupe, provenance, status, frontier); never copy coordinates by hand (`pk.py get`) |
| engine | `pk/engine.py` | `fq_quench(inp, QOpt)`: the one fq call path; `QOpt` = every fq flag (incl. `freeze`, 10-22) |
| registries | `pk/moves.py` (`@move`), `pk/fill.py` (`@placer`), `pk/schedules.py` (named QOpts), `pk/anneal.py` (`Schedule`) | building blocks, each with a declared parameter space: unspecified parameters are drawn, every value used is logged |
| stages, chains | `pk/stages.py`, `pk.py stages` | one contract over all registries: kinds `gen`, `ruin`, `recreate`, `constrain`, `relax`, `perturb`; a chain is a string `a(k=v) \| b \| ...` |
| evaluation | `trial.py` (any n, suites), `pk/battery.py` (n = 110 idea battery), `pk/benchkit.py` (event races), `pk/corpus.py` (replay) | verdicts on refined outcomes per measured CPU, never on raw output |

State between stages (`pk.stages.State`): the packing, a per-square `frozen` mask, the `region` the next recreate fills,
and generator facts in `meta`.  Freezing is state, so every later stage honours it.  Relax passes it to fq as
`--freeze`.  Perturb stages refuse frozen squares (release first).

## Adding an idea
1. **Pick the kind.** Kinds: a start from nothing (`gen`), opening a region (`ruin`), filling it (`recreate`), freezing or
   releasing (`constrain`), a T = 0 descent (`relax`), a stochastic move (`perturb`).
2. **Register it where its family lives.**  A placer goes in `pk/fill.py` (`@placer`, becomes `recreate.<name>`), a
   move in `pk/moves.py` (`@move`, becomes `perturb.<name>`), an fq recipe in `pk/schedules.py` (becomes
   `relax.<name>`), anything else in `pk/stages.py` (`@stage(name, kind, **space)`).  A private idea registers the same
   way from a plugin (below).
   - Declare every knob in the space.
   - Return an info dict.
   - Draw all randomness from the `rng` argument.
3. **Check the contract.**  `python3 pk/tests/test_stages.py` checks that every stage:
   - logs all its parameters,
   - is deterministic per seed,
   - leaves frozen squares alone.

   It covers a new stage with nothing to write.  `pk.py stages` lists it.
4. **Look at a few runs.**
   ```
   from pk import stages as S
   st, log = S.run_chain('gen.file(path="runs/x/n{n}.txt") | constrain.freeze | recreate.bestofK(K=16) | relax.round.15c', dict(n=548), seed=1)
   ```
5. **Pre-register.**  Add a `### Pre-registration: <title> (<NAME>, <date>)` block to PACKER.md (or the plugin's
   NOTEBOOK) with hypothesis, metric and rule.  `trial.py run` refuses a NAME without one.
6. **Run the trial.**
   ```
   trial.py run NAME baseline --suite SUITE
   trial.py run NAME ARM --head SPEC --suite SUITE
   trial.py report NAME
   ```
   Then write the result under the pre-registration, and in `CLAIMS.md` if it steers strategy.

Why the trial appends a tail: a generator's output is a distribution of starts for refinement (CLAUDE.md rule).  Closing
a line on raw fit rate is the mistake this prevents (it happened once, 10-22).

## Plugins (10-23)
Stages, placers, moves and trial suites can live in a separate (e.g. private) checkout.  `pk/stages.py` imports
`pk_plugin.py` from each directory on `PK_PLUGINS` (default: `../plugins/packer` next to this repo's root, if present;
`PK_PLUGINS=` disables) before it adapts the registries.  A plugin registers with the same decorators and may define:
- `TRIAL_SUITES`: suites for `trial.py` (same fields as `trial.SUITES`; `tail` None = the standard tail);
- `NOTEBOOK`: a markdown file `trial.py` also searches for pre-registrations;
- `TEST_STARTS`: `{stage: chain}` for stages that need a particular start in `pk/tests/test_stages.py`.
Scripts in a plugin find this directory through their own bootstrap (they keep `runs/` and the store here).

## Not done yet (TODO.md "Packings")
- Loops: `repeat(k, ...)` and `accept(metropolis, T)` combinators, so ruin-and-recreate and basin hopping become chains.
- Hand-wired scripts that predate the stages and still compose by hand: `explore.py`, `chain.py`, `hop.py`,
  `compose2.py`, `relax.py`, `fresh.py`.  Port each one when a trial shows the port matches it.
- `trial.py` statistics: sequential stopping and rate posteriors, as in `pk/battery.py` / `pk/benchkit.py`; more
  suites (a public one, e.g. n = 110).
