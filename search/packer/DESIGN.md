# packer tooling redesign (10-09, branch `tooling`)

Why: exploration so far is good but ad hoc (Evan 10-09).  Costs: coordinates copied by hand between tools; every test of
a change is a 5-30 min explorer run with ~10x replicate spread; 57-88 % of CPU is full polish, much of it on basins we
already have; the archive's "new" count is screening-level (110, c110_h3: 5,756 "new" entries, ~17 certified new), so
bandit rewards are mostly noise; knobs (noise, pressure, schedules) are scattered flags, so sweeping them is expensive.

Goal: support tooling first.  A store every packing flows through, an engine with explicit knobs, and a replay bench
that grades a change in minutes on frozen proposals.  Modularity is the means; context economy (no file dumps, no hand
copying, one name per packing) is the payoff.

## Pieces (`search/packer/pk/`, CLI `pk.py`)
1. **Store** (`pk/store.py`, sqlite at `runs/store.sqlite`, gitignored).  One row per packing: n, s, coords (f64 blob),
   source / finder / ref (commit, issue, run path), parent, status (raw < screened < polished < certified | not-min |
   infeasible), S_exact, fp, D4-canonical hash (dedupe), tags.  A `frontier` table: best known per n (register + pending
   issues + ours).  Importers: any file format (`import_ext.read_any`), register witnesses, pending-issue cert repos,
   explorer archives, our candidates.  Export in any format by id or "best of n".  Nothing gets copied by hand again.
2. **Engine** (`pk/engine.py`): `fq serve` (persistent process, text protocol on stdin/stdout) instead of a process and
   two temp files per quench; stages screen -> polish with identity checks (side + contact fingerprint) between them,
   kick-retry on non-jammed stalls.  Every knob (loosen = pressure, mu0 = stiffness, pit, flip-top, stag-tol) is a field
   of one options object.
3. **Moves** (`pk/moves.py`): registry; each move declares its parameter space (sigma range, region radius, melt
   pressure / sweeps / rounding, ...), so a bandit arm can be (move, parameter bin) and a sweep is a list of settings.
4. **Replay bench** (`pk/corpus.py`): freeze a corpus of proposals (parent id from the store, move, params, seed) at a
   few n; replay through any engine / options variant; funnel per item (screen side, class, polished side, fp, time,
   certified class via a shared exact cache).  Metrics: distinct certified new basins per CPU-h, near-frontier yield,
   false-new rate, stage times.  Same proposals for every variant -> paired, low variance, minutes.
5. **Explorer v2** on 1-4: archive = store rows tagged with the run; pluggable parent policy (gap-stratified quotas,
   lineage) and arm set; one job queue across n for chain (async seeds, no idle slots).
6. **Library** = store queries + tags: grid-beating (s < k), near-frontier (gap to best known < delta), above-grid
   robust (s in (k, k + 0.05), no full lines, low grid-collapse rate under standard kicks; measured, cached).

Order: 1, 2 (+ fq serve), 4 (needs 1, 2), 3, 5, 6.  Each step keeps the old scripts working until v2 matches them on
the bench.

## Data point behind the arm design (10-09, `proposals.jsonl` of c110_h2/h3, hunt1, hunt2)
Parent gap above the run's starting best vs child outcome.  CPU mostly goes to far parents (hunt1: 29 of 33 CPU-h on
gap >= 1e-2); parents >= 5e-2 above gave 0 new basins within 1e-3 of the best in ~41k proposals; parents within 1e-3
give 35-175 per CPU-h.  Big single-step drops are rare with our (mostly local) moves; lineages and big-move arms are
untested, and the corpus bench is the place to test them.

## Status (end of 10-09 session)
Done: 1 (store + importers + `pk.py`), 2 (`fq serve`; staged polish tried and dropped: worse), 3 (moves incl. `anneal`,
`cross`), 4 (`pk/corpus.py`) and the battery (`pk/battery.py`) as the idea-grading front end; plus `anneal sched`,
`pk/crossing.py`, `pk/umbrella.py`, `pk/rediscovery.py`, `pk/bridging.py`.  Not done: 5 (explorer v2 on pk: the explorer
still has its own moves / pipeline; recipe fixes were applied to `explore.py` directly), 6 (library tags / grid-collapse
robustness per packing).  Store size ~0.8 GB (gitignored); rebuild: `pk.py sync-register`, `sync-pending`,
`import-run runs`, seed pools (see PACKER.md).
