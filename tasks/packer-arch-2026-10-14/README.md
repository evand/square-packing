# Packer architecture review and TT2 plan (2026-10-14, Evan + Claude)

Map of the system (engines, store, proposals, s(110) apparatus): https://claude.ai/artifact/M6VSo1eq5twYqb8T6w3XzB (private).

## Findings (checked in code 10-14)
- chain launched explore with `--adapt` but no `--novel-ref` and the legacy reward: the bandit rewarded run-local novelty
  (~35-50 % of proposals "rewarded"), against the 10-08 (global novelty) and 10-09 (`--reward value`) decisions.  Fixed
  10-14 (chain `--reward value --novel store` defaults; `--reward legacy --novel run` reproduces the old behaviour).
- No explorer run after 10-12 18:46 was imported into the store (NH / NB / H5 missing); store-based novelty / diversity
  needs every production run imported.  Imported 10-14; chain should import its own rounds (TODO).
- `packing.fp` NULL in every row; fq's contact fp is computed and was thrown away.  Archive entries keep it as `cfp` from
  10-14 (label-dependent: same packing relabelled -> different fp; fine within a lineage).
- KindBandit half-life is 600 s of wall time: ~10k proposals of memory at n < 50, ~250 at n > 250 (NH).  Study: A1.
- Parent / start choice is per-mechanism (chain best(n) snapshot + n+-k seeds + carries + run-local library; explorer
  cell-uniform; composer cell bandit over pieces; harden / relax / vband by hand); no cross-run memory of what was
  expanded (explorer `expanded` lives in the run).  Ledger table `expansion` added 10-14; chain writes it per round.
- Several "best known" definitions: `store.best`, `store.frontier`, `eventbench.frontier_now` (careful: events +
  store + envelope + L-step), `mondrian.frontier_best` (dead code), chain's snapshot.  Chain's silent fallback to the
  10-05 register removed (fails unless `--allow-missing`).
- Duplicated implementations: fq call wrappers (hop.quench, pk/engine fq serve, compose2.polish, harden.py, relax,
  rectlib, fresh; + slp2 via mcmin), packing I/O (mcmin.load_deg, packing.read, import_ext.read_any, contacts.load),
  neighbour transforms (chain seed_job, harden xforms, arms n+-3, crop, insert110, excise), parent pickers (7), moves
  (mcmin, explore, pk/moves wrappers), similarity (side key, image fp, tilt signature, contact fp, matching distance),
  grid checks (layout.full_lines, packing.grid_lines, pipeline), pools (jobpool, bare ProcessPoolExecutor, pk/jobs),
  piece stores (atlas refs, rectlib files, store `piece`).

## Change classes (track per change)
- bit-exact: same proposals and outputs (refactors); check by replay (`fqreplay.py`, `chainreplay.py`).
- additive: writes new data, changes no decision (ledger, cfp, cluster column); check decisions unchanged.
- functional: changes what is proposed / picked; needs a functional test (TT2 or a slice).

## Work list (10-14 triage)
Now: N1 chain reward / novelty / as-of / loud src (functional, done); N2 ledger + cluster column + cfp (additive, done);
N3 docs.  Agents: A1 bandit timescale study (`search/packer/analysis/bandit_timescale.md`), A2 similarity / cluster
study at 110 (`analysis/similarity110.md`), A3 ingestion watcher (`pk/watch.py`).
Next: parent sampler `pk/parents.py` (pool query n, n+-k, cluster, gap band, as-of; policies; uniform cluster then
uniform within = default), run config (test / record modes; grid-rule knob), chain imports its rounds, cfp into store
on import; parent-distance model after A2 (fit from logs, then a randomized gap-band slice); near-miss queue + finisher;
one fq call path through QOpt (bit-exact); shared neighbour transforms (bit-exact if RNG kept); piece sources (free
W x H, crop anywhere, n > 324 / catalogue revisions, rectangle records).

## TT2 (next major test; Evan 10-14: "our best guess at the strongest collection of settings")
Time travel: frontier and starts as of a cutoff T (eventbench `starts(T)`, `frontier_at`, `score` incl. the recovered
fraction of later gains); parent pool / novelty ref limited to store rows added before T (store times are local).
Arm = strongest current settings (value reward + global novelty, cluster-uniform parents, A1's bandit timescale, ledger,
fresh ingestion as of T) vs the pre-10-14 chain (`--reward legacy --novel run`), concurrent, compute-matched
(~2 x 30 CPU-h).  Metric: later records recovered per CPU-h + recovered fraction of later gains.  Pre-register in
PACKER.md before launch.  Prior pair (time-travel 10-12, 2 x 33 CPU-h) found the 154 record.
