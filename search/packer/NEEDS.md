# What the packing program needs (2026-10-09; evening update first)

## 10-09 evening update (tooling session; numbers in `PACKER.md` § Tooling session, design in `DESIGN.md`)
* **Framing (Evan):** compute-efficient new basins come from new *ideas* (starts, proposals, schedules); compute-efficient
  evaluation matters mostly because it speeds up grading ideas.  Small-n statistics can't resolve tuning-sized (< 2x)
  effects; don't spend effort there.
* **Tools now:** store + `pk.py` (no hand-copied coordinates), `fq serve`, move registry with explicit parameters,
  replay bench (`pk/corpus.py`), idea battery (`pk/battery.py`: one command, P(candidate > baseline)), `anneal sched`
  (pressure / corner-radius / rotation / Q4-bias curves), crossing harness, umbrella / WHAM.
* **Learned:** (a) the screen-return shortcut lost half the new basins: polish everything; (b) single-move reach is
  ~0.05 matching distance for every move incl. anneals; recombination reaches no further than its nearest parent but uses
  near pool members ~35x better than the arm mix (two parents, random mate, is best) -> progress = pool expansion toward
  targets + recombination; (c) shallow whole-box anneal preserves sub-grid structure (positive control 92 % vs 62 %);
  (d) deep / global anneals and Q4 bias go to the grid; at P = 30 axis alignment wins both entropy and density; 0 grid
  crossings at 110 in ~3,000 trials (one long jump from a start 0.069 outside the funnel); (e) census 110: 683.
* **Needs now:** port crossover + ashallow into the explorer; lineage test (chained pool expansion); battery: tighter
  leakage filter, second n, throughput; crossing: region re-packing (block moves), line-load bias.  The older needs list
  below still stands (generator benchmark = the battery; fill / completion; chain throughput; polish at large n).


Written at the end of the 10-09 session (targeted hunt h2, generator tests, hunt2) for whoever picks this up next.
Log with numbers: `PACKER.md` § "Targeted hunt h2".  Tools: `README.md`.

## Where we are

Factor the search into three pieces (Evan, 10-09):

1. **seeds from scratch**: build an arrangement that sits in a good funnel;
2. **perturbation of seeds**: move from a decent packing to better nearby minima;
3. **luck**: a perturbation lands below the record.

**(2) works and is measurable; (1) mostly does not.**

* Perturbation: the explorer (`explore.py`; kicks, structured moves, melt, Thompson move choice with global-novelty
  reward) finds new minima near existing sub-grid packings at a steady rate (s(110): 493 certified sub-11 minima, discovery
  curve still near-linear; s(109): ~3,400 screening-level 110-style basins in a few CPU-hours).  It found the records
  132 and 155 (hunt1) and our 266 / 270 / 272 by perturbing other people's recent packings.  Held-out evaluation works:
  rediscovery pairs (`rediscover.py`, `runs/rd1`: detection falls with matching distance, ~70 % below 0.02, rare above
  0.06), polish replay, graded benchmark tiers (`bench_explore.py`).
* Seeds from scratch: every from-scratch constructor tested so far reaches the *neighbourhood* (non-grid, 0.01-0.2 above
  the record) but not the *funnel*.  Evidence (10-09 unless dated):
  - 110, line-free layout `hyb` / `ftmc` seeds (11.0000067-11.05; 10-04 genomes, no sub-11 input) → explorer 4 procs x 32
    min: 980 basins, ~960 new, **0 below 11**; 3,938 proposals screen-grid.
  - 110, diagonal strips (`strips.py`) → explorer: 0 below 11, 2,662 / 2,900 proposals screen-grid.
  - 126, staircase seeds with **Xu's measured column shifts** (11.90, 11.92) → explorer: **Xu's 11.7426406871 in < 4 min**.
    Same structure found **blind** by `stairgen.py` (CEM; 11.9497, 11.9701) → explorer 8 min: 11.9084, not Xu.  The gap at
    the start is similar; the arrangement decides.
  - Gadget transplant n → n + 1 (`transplant.py`), insertion into pose holes (`chain.py --seed-ks -1`), removal from
    non-grid lineages (90 → 89): all jam at the grid or never gain the extra square.
  - Older (10-04 to 10-08): layout genomes 11.012-11.03 at 110, `ftmc` 11.008 ("relatives, not the records"); disk anneal
    0.04-0.1 above records, all grid at 110; crops of 110-style records enter the known funnel's rim.
* Consequences: **records where the grid is the record (90, 111, ...) are not our strong suit**, even at large n where a
  better packing likely exists.  Regular (fixed-family) records are mixed: large 45-degree diamonds do not perturb
  (s(89): 300 sub-10 basins, 3 structure classes, all the 7 x 7 diamond with a rattler moved); 110-style packings perturb
  easily (s(109): 3,397 non-diamond sub-11 basins, 29-69 tilted, but the floor is 10.951872, 0.0021 above the 45-degree
  record); thin strips should perturb (untested).

## Where to spend compute (until (1) improves)

* n with a sub-grid packing that is recent / single-source / not yet perturbed by us: hunt passes over all open n
  (hunt1 found 2 records in one pass of 112 n), other people's new packings (pending registrations, n > 324: Couzo
  327-379), and the per-n sub-grid library (below).
* s(110) as a census, not a record hunt: per 40 CPU-h (~500 new certified minima) roughly **5 %** chance of a sub-record
  minimum we can find (10-09 estimate from the low tail: 0.5 % of new minima land within 1e-4 of the record, ~6 % within
  1e-3, nothing below it in 474 non-start discoveries; count within delta ~ delta^0.5-0.8).  Ignores everyone else's compute.
* 10-09 outcome of this policy: n > 324 (Couzo's single-source packings) gave the one new basin of the day (375,
  −4.3e-5) in ~6 CPU-h; a second pass over 84-324 (hunt2, ~25 CPU-h) gave only a precision refinement (270); s(110)
  ~12 CPU-h: 30 new certified minima, none near the record.  So: the frontier of less-worked n is where perturbation
  still pays; extend it (n 380-500 from Couzo / mirror / own constructions), and keep a census-level budget at 110.
* Not: re-exploring saturated records (h2r1: 8 hard targets x 100 rounds, nothing moved; 89 / 147 / 290 / 291 / 109 had
  hunt1 rounds too).

## Needs, in priority order

1. **A benchmark for seed generators** (cheap, do first).  Held-out: n with known sub-grid records (110, 126, 132, 155),
   records and their descendants withheld; metric = CPU time to the first sub-k certified minimum and to within 1e-3 of
   the record, from generator seeds + a fixed explorer budget.  Today's runs are the first data points (110 layout: no
   sub-11 in 128 CPU-min; 126 informed staircase: record in 8 CPU-min; 126 blind: 11.908 in 16 CPU-min).  Without this,
   generator work cannot be graded (a generator's value only shows after perturbation).
2. **A fill / completion that can express new-style records.**  The exact axis fill (`gen.mis_fill2`) puts axis squares
   on wall-anchored lattice lines plus short staircases; ry-xu's 126 needs axis squares in pockets against the block (both
   coordinates block-anchored) and against near-axis rotated squares: with Xu's exact 36 45-degree squares the fill gives
   125, residue-class candidates (171k) cover only 61 / 90 of his axis squares.  Candidates: a soft fill (fq with the
   tilted part frozen or weighted, axis squares from a loose lattice), or fill + near-axis rotations as MILP columns.
3. **Completion without jamming.**  `fq quench` cannot absorb one extra overlapping square: the ALM pushes whole rows out
   and lands on the grid (147 = 146 + 1 → 12.71-13.0; 89 = 88 + 1 → 10.0; Xu-like 125 + 1 → 11.90-11.92).  Needed for
   n → n + 1 seeds and for generator completion.  Ideas: insert at a dilated side and shrink with rows locked; local melt
   (`anneal melt`) around the insertion; a squeeze that forbids full lines (`layout.full_lines`) during descent.
4. **Score seeds by side, always with n squares** (Evan 10-09).  `stairgen.py` does this now (fq-quenched side of a
   complete packing, integer-jammed results rejected); `blockgen.py` / `strips.py` still count squares at fixed s (fine
   for reproducing 45-degree records, wrong as a search objective).
5. **Sub-grid libraries everywhere** (Evan 10-09: sub-grid packings are good seeds because they are not grid-jammed).
   Done in `chain.py` (`--lib-init`, `--lib-k`; per-n `<out>/lib/`).  Open: certify library entries (screening-level
   now; ~20-25 % of polished candidates are not minima); keep them across runs (one global per-n store instead of per run);
   a diversity measure per n (structure classes: tilted count, angle groups) to rank n for perturbation vs seed work.
6. **Generalising across sizes.**  Removal seeds (n + k → n) mostly copy the record or jam; transplant and insertion fail
   (above).  The working cross-size path so far is other people's packings at nearby n as starts, and lineage chaining.
7. **Chain throughput**: seed jobs run synchronously in the main loop (`chain.py`, `run_jobs(seed_job, ...)`), so at
   n ~ 260 slots sit idle for minutes while removal seeds quench (hunt2: 4 of 6 slots busy).  Run seeds in the background
   or prepare them one round ahead.
8. **Polish at large n**: fq's flip search (`--flip-top 8` default) stops early: at 375 it left a corner-corner slip
   that exactsolve's MILP found (dS −4.9e-6); `--flip-top 128 --pit 600` finished it.  Scale flip-top with n, or feed
   the MILP's descent branch back into fq.  Also: polish every imported start once (Couzo's 270 / 375 were unconverged;
   the chain then reports the polished start as a "candidate").  And fq can stall at non-jammed states (s(110) 10-09:
   states 5e-5 above the record that exactsolve calls not jammed; `--kick 1e-4` + re-quench sends them to known minima):
   the explorer's screening novelty overstates (249 "new" → 13 certified new).  Cheap fix: a kick-retry when the jam
   test fails, before a basin counts as new.
9. **Grid-record n (90, 111, 183, ...)**: needs (1)-(3) first.  90: best non-grid lineage 10.00905; 45-degree rectangles
   reach 89 / 90 squares at 9.9999.

## Tools added 10-09

`../exact/batch/regsync.py` (register refresh + diff), `import_ext.py` (others' formats → seed pool; f64 check),
`transplant.py`, `strips.py`, `blockgen.py` (reproduces 45-degree records 65, 66, 89, 150 to 1e-9), `stairgen.py`
(staircase blocks, CEM, scored by quenched side), `chain.py` (`--seed-ks` incl. insertion, `--seed-tries`, multiple
`--carry`, sub-grid library).

## Gotchas met 10-09

* numpy 2: `repr(np.float64(x))` is `np.float64(x)`; it got written into seed files (fq then fails to parse).  Cast to
  `float` before writing.
* NixOS has no `/bin/bash`: run scripts as `bash script.sh`, not via the shebang.
* `pkill -f PATTERN` / kill-by-pattern from a shell whose own command line contains the pattern kills that shell; match
  on `$1 ~ /python/` in `ps` output, or kill by saved PID.
* Register sides can be rationals (`n/d`): parse with `Fraction` (`import_ext.side`).
