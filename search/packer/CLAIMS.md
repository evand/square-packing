# Claims register (packer)

Started 2026-10-11 (review: `REVIEW-2026-10-11.md`).  One line per strategy-steering claim, so TODO / NEEDS / plans cite an
ID instead of restating a number.  Fields: claim | definition / scope | evidence (path, query) | events, replicates |
status | stale when.  Status: **solid** (re-derivable, enough events, concurrent control), **plausible**, **weak** (few
events, no control, or confounded), **stale** (superseded by a later measurement), **untraceable** (no saved derivation).
Everything is n = 110 unless stated (transfer to other n / families: out of scope for now, Evan 10-11, `CONFUSIONS.md` X16).
Compute units: before 10-11 "CPU-h" / "proc-h" = worker wall-clock (wall x procs, or summed per-item wall); from 10-11
`manifest.json` `cpu_s` = measured CPU.  "Minimum" = `pk/identity.py` (exact S to 1e-17); "certified" excludes weak.

| ID | claim | definition / scope | evidence | events / reps | status | stale when |
|---|---|---|---|---|---|---|
| C01 | Polish everything beats the 1e-6 screen-return shortcut | new distinct certified minima, frozen corpus | `pk/corpus.py report c110a` (10-11: full 33 vs base 20 minima; side keys 37 vs 21; the 10-09 log's 30 vs 16 read an older basin table) | 1 corpus, paired | solid | report re-read without a snapshot (X03); other n |
| C02 | Single-move reach ~0.05 matching distance | withheld minima rediscovered by parent distance band | `runs/corpus/c110r`, `pk/rediscovery.py c110r full --dist runs/dist110.json`; landing spread `move_spread.py` (median 0.003-0.009) | 9 hits, all from 0.02-0.05; 0 / 1,144 from 0.05-0.1 | plausible (±0.01: rattlers inflate d) | rattler-masked distances; other n |
| C03 | Two-parent crossover uses near pool members far better than the arm mix | battery reach suite, withheld minima | `pk/battery.py report cross2 --vs baseline`: 11 vs 2 minima | baseline 2 events | direction solid, size weak (~35x, CI ~10-300x) | more baseline events |
| C04 | Lineages cover 0.05-0.08 inside the sub-11 region in minutes | time to the record from vetted starts | `runs/lineage/L1`, `L2` (PACKER 10-10) | every cell reached the record; 1 replicate per cell | solid as stated | — |
| C05 | "110 below 11 is one funnel" | all known sub-11 minima within 0.1 of the record (2 / 537 beyond) | L1 / L2 + census | census built from record-family starts | weak (circular) | census v2 far / from-scratch arm |
| C06 | 0 grid crossings at 110 | anneal-type trials from above-grid / grid-like starts reaching < 11 | `runs/crossing/{x110a,rate110,depth110,bias110,lag110}` | 0 / 3,130; positive controls 208 | solid as stated | definition-dependent: rev1 had 31 / 720 crossings from above-11 *minima* (kicks); from-start anneals had no expansion phase |
| C07 | Census 110 = 720 certified + 40 weak distinct minima | store table `minimum`, n = 110, below 11 | `census_certify.py report --n 110` (10-11) | — | solid | any certify / descend run adds rows |
| C08 | Census discovery rate ~1 new certified minimum per worker-h | live explorer at 110 vs a frozen reference | `runs/insert110` baseline arm | 2 events (CI ~0.1-3.6) | weak | census v2 (measured CPU, more events) |
| C09 | Share of explorer outputs that are not minima | sub-11 outputs, exactsolve not-min | 20 % ex0; 27 % c110_h2/h3; ~35 % of all sub-11 archive entries in L1/L2/insert110 (`parent_retro.py`); 90 % of *first-seen* sides in the 10-10 lineage runs | several | stale as one number (TODO "Decisions 10-08" still says 20-25 %) | — (cite the population) |
| C10 | The explorer melt bug lost 22-40 % of melt proposals to the grid; fixed melt loses ~0 | paired old / new melt, same seeds | `runs/melt_ab`, `melt_ab.py --report runs/melt_ab` | 320 pairs (110, 130) | solid | — |
| C11 | Kick + full re-quench is a hop, not a finisher; kick + polish-only is safe | controls = certified minima | `runs/kr_eval` (18 / 60 controls moved lower), `runs/kr_eval_noalm` (1 / 60; 58 % of not-min to descend destination, 7.4 s CPU / 4 tries) | 120 points each, same seeds | solid for these points | other n (X16) |
| C12 | A not-min parent behaves about like its own minimum | children certified / moved per 1000 proposals | `parent_retro.py` matched pairs: 572 vs 607, 493 vs 493 | 308 pairs, observational | plausible | census v2 controlled arm |
| C13 | Never-descended (far?) not-min parents have ~2x the discovery rate | within-run first reach | `parent_retro.py` split: 153 vs 77-79 per 1000 | 985 proposals, selected | weak (hint for Evan's distance rule) | census v2 far-state arm |
| C14 | Mining the 110 not-min backlog with kick + polish is not worth it | new certified minima per measured CPU-h | `runs/census/n110/finish_backlog` | 0 new in 0.56 CPU-h, 310 points | solid for the backlog (descend's hard tail); says nothing about inline finishing | — |
| C15 | Thompson move choice with global-novelty reward 2.3x | new basins, ex5 vs ex6 | PACKER 10-08 | 1 pair; screening-level novelty | weak (still a TODO "Decision") | census v2 |
| C16 | Copy-rich minima (2-7 % of minima) take 20-36 % of parent picks | archive entries mapped to descended minima | PACKER 10-10; no saved script | 8 runs | untraceable (now recomputable from the `copy` table) | — |
| C17 | ~5 % chance per 40 CPU-h of a sub-record minimum at 110 | tail fit of new minima near the record | NEEDS.md (10-09) | 0 events below the record | weak (premise ~500 new / 40 h is 5-12x above measured rates) | — |
| C18 | Mondrian formula composer: 0 wins / 112 ties; corner overlap 0 / 164 | n <= 324 vs store frontier, tie <= 1e-9 | `runs/mondrian/compose_all.json`, MONDRIAN.md | — | solid for composer v0; not a verdict on the method (pieces pruned, squeeze is real search, v2 untested) | Mondrian v2 |
| C19 | Insertion from 109 donors adds no new 110 family | new minima vs frozen 758 | `runs/insert110` | donors 0, baseline 2 | plausible ("no new family"); rate comparison unresolved | — |
| C20 | Seeds from scratch do not enter the 110 funnel | sub-11 from generator seeds + explorer | PACKER 10-04..10-09 (layout, ftmc, strips, stairgen) | 20-128 CPU-min per generator | weak as a method verdict (budget, fill limit, completion) | seed benchmark (NEEDS 1) |
| C21 | Temperature is not the coordinate that opens the funnel | minima PT, anneal rate vs depth | pt1 / pt2 (slp2 era), rate110 / depth110 | 2 PT runs, 30-trial cells | plausible | PT with fq + explorer arms |
| C22 | Fresh single-source packings are the best hunt targets | records per hunt pass | hunt1, hunt3, big1 (~10 records) | — | plausible ("records resist exploration" contradicted by ch2, hunt1 s(132)) | — |
| C23 | Dedicated runs do not beat the chain per proc-h | hunt3 follow-ups | PACKER 10-10 | not compute-matched | weak | compute-matched test (TODO) |
