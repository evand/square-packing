# Mondrian composer (10-10 evening; plan in TODO.md § Packings)

Copy of itsnaka's SQUISH "Mondrian composition" (jlevy/squares#481, 15 records n = 131..307; certs and README in
`runs/ext/itsnaka-squish-certs-d45669b/squish-submission-2026-10-09/`).  Code: `mondrian.py`.

## Sources
* Kingbird catalogue as mirrored 09-30 (`site/data/packings.json`, 545 valid of 549, incl. 208 `_r*` revisions).  The
  "former" records #481 cites (152_r2, 202_r3, 241_r3) were the current files on 09-30, so no fetch was needed.
* Store: every non-run sighting (register, pending, outside pools), every certified minimum, best 30 distinct sides per
  n from our runs: 22,976 source packings in total.

## Atlas v0 (corner crops)
Tight corner crops beating the grid (m > floor W floor H), Pareto per (m, floor W, floor H): 77,670 pieces, 20 s on 2
procs (`runs/mondrian/atlas_all.npz`).

## Composer v0 (formula count) — not enough
T1 (piece + L of uprights) and one guillotine cut (two full-height columns, one piece each; covers T2/T3) with a
closed-form upright count; `best(s)` is monotone, bisect per n (~0.05 s per evaluation).  n <= 324, all n, same-n
pieces excluded: **0 below the frontier, 112 exact ties (records that are T1 compositions), ~15 within 5e-3.**
Frontier-only pieces (first check): same picture.  Beyond 324 (out of scope for now, Evan 10-10): Kingbird's s(416)
minus a row and column gives s(375) <= 19.899495 vs Couzo 19.907053 / ours 19.906982 (coverage gap, not new structure).

## Replication of #481 n = 259 (the method, measured)
* Formula count for their two pieces: at most 251 (they have 259).  Exact axis fill (`gen.mis_fill2`) keeping the
  crops' own uprights: 257-258 for every crop / orientation / corner choice, also at S + 0.01.  Composing at a slack
  side (16.8) then fq-squeezing lands in other basins (16.643-16.67).
* MILP fill around *their* 49 tilted squares as placed recovers all 210 uprights (259): the fill is adequate.
* Rigid match of their clusters to the sources: s(123) piece = its left 8.59 columns (orientation unchanged) at the
  upper-right corner; s(50)_r5 piece = transposed, **shifted down one unit** (bottom row dropped), tilted squares
  rotated -0.27 deg by their squeeze.
* Rebuilt: tilted squares of both pieces in place, **all uprights re-filled by the MILP** -> 259 at S = s(123) + 5,
  feasible as composed (pen 2e-10); fq -> **16.5913781455 = their S**.  Keeping the crops' uprights gives 258.
* So: a piece is a set of tilted squares + an integer offset from a box corner (rows / columns of the source dropped);
  the uprights are always re-filled exactly.  The closed-form count is a screen (it undercounts by up to ~8); the
  composer needs an exact-fill stage.

## Replication harness (`mondrian.py replicate`, 10-10 night)
Benchmark (Evan): from #481's listed pieces only, can our pipeline (compose -> fq -> short explorer) beat the register
value *before* #481 (status quo ante, `ANTE_481`)?  Tiers: REPLICATED (<= their S), BEAT-ANTE, no.
* Piece model: crop of a source (8 orientations), i columns / j rows dropped at the wall side (<= 3), tight inner cut
  holding the listed m squares; only its tilted squares are kept, placed at LL (A) and UR / LR (B), conflicts dropped.
  Uprights always re-filled by the MILP.
* Pair screen: the rectangle formula over all pairs (vectorised), slack 8, then joint MILP on the best 150 (+ both
  conflict orders), then fq + polish.  The additive single-piece-excess estimate fails below the grid bound (the band
  crossing every row / column needs both pieces: superadditive), so it was dropped.
* Sources: #481's "former" Kingbird records 152_r2, 202_r3, 241_r3 are NOT our 09-30 files (152_r2 = 12.8467, angles
  35/55 deg; ours = 12.8310, 28-30 deg): fetched 10-10 to `site/data/ellsworth/svg/new_*.svg` (gitignored mirror).
  **Atlas gap: the 09-30 mirror lacks older revisions too (152 had no `_r*` locally); the store has none of them.**
* Forensics (rigid match of source tilted squares onto their final packings): offsets are small integers (259: (0,-1)
  and +3 columns; 292, 302 similar); for 131, 153, 207, 209, 270, 305 only part of the source clusters match (squeeze /
  "moved" squares change them); 153's final S = ry-xu s(129) + 1 exactly, angles 39/51 vs the source's 35/55.  So
  their squeeze is a real search, and exact-composition replication is the wrong gate for those n.
* Pass 1 (compose + fq, sides in [S, S + 0.02]): 259 REPLICATED (seconds); 236 composes at s(88) + 6 = 15.8825 but fq
  does not squeeze (theirs 15.8640, ante 15.8678); 131 / 153 / 207 / 209 / 232 short by 1-5 squares in that window.

Job convention (after three self-kills via `pgrep -f` / `pkill -f` matching the controlling shell's own command line,
a `setsid` fork that made `$!` wrong, and output lost behind `| tail` pipes): use `pk/jobs.py start NAME -- CMD` /
`status` / `stop NAME` (process-group kill, pid recorded at launch, always a log file).  Never pattern-match command lines.
Positive control: `mondrian.py selftest` (rebuilds #481 n = 259, PASS in ~20 s).

## Kingbird catalogue in the store (10-10 night, Evan)
`pk.py sync-kingbird` (`pk/kingbird.py`).  The catalogue links 862 SVGs from 10 pages (we had only crawled 8; most
revisions are linked only from `compared2` / `compared3`).  Snapshots: `runs/kingbird/2026-09-30/` (the old site
mirror, 549 files, imported as history) and `runs/kingbird/2026-10-10/` (636 fetched, 189 revisions reused, 37 dead
links): 825 parsed, 640 new store rows in total.  Flags kept in the manifest: 13 labelled variants overlap
(`__invalid`, `_approx`, `_start`, `_asymmetric`; stored infeasible); `square-103_r0` overlaps 0.017 (old revision,
inexact as published?); `square-108_r0` fails in the site parser (`missing ref #sideA`).  Store backed up first
(`runs/store_backup_2026-10-10_pre-kingbird.sqlite`).  Atlas v2 from the catalogue: 23,240 sources, 82,954 pieces.
`mondrian.load_src('kingbird:NAME')` reads the newest snapshot.

## Trimmed replication pass (10-11)
First attempt (10-10 night) wedged 10 h: a scipy/HiGHS `milp` in the joint fills ignored its `time_limit`; 12 workers
at 100% CPU, no output (my watcher only waited for the summary line).  `exact_fill` now runs every MILP in a forked
child under a hard wall-clock limit (3 x time_limit + 30 s; overrun -> killed, input saved to `runs/mondrian/hung/`,
tilted squares returned unfilled = safe undercount); jobs relaunched under `pk/jobs.py` (trim1..3, ~45 min, no overrun).
Settings: `--procs 4 --top 150 --window 0.1 --extra 2 --fq-top 4 --relax --insert-max 4 --tries 2 --explore-min 10`.

| tier | n |
|---|---|
| REPLICATED | 259, 263, 292 |
| BEAT-ANTE | 232 (15.7687164 vs ante 15.7781746; theirs 15.7674283) |
| no | 131, 153, 207, 209, 236, 269, 270, 302, 303, 305, 307 |

Forensics (rigid match of each source's tilted squares onto their final packing, 8 orientations; tight = 0.05 / 1 deg):
the three REPLICATED n are exactly those with both pieces rigid in the final (259 / 292: s(123) 33/33 at UR shifted 3
columns, s(50)_r5 16/16 at LL; 263: s(202)_r3 78/78, one row dropped).  232 and 302 have one rigid piece (one row
dropped); 131, 209, 307 partly rigid; 153, 207, 236, 269, 270, 303, 305 not recoverable (relaxed descendants).
Wherever both pieces are visible they sit in opposite corners separated by an axis cut with a 1-2 unit band of
uprights: no interlocking.  #481's pieces are mostly OFF the atlas frontier (259: +1 / +2 pieces in cells whose best is
+5 / +9; 292 +1 vs +6; 263 +6 vs +14) and from former records: fit (fractional sides vs S) and relaxability matter,
not raw excess.  The atlas's Pareto prune would even drop their 259 s(50)_r5 piece (a same-m strictly smaller crop
exists).  Next: diversity atlas, rank by squeeze gain, blind test (TODO).

## s(126) record reads as a corner-overlap composition (10-11, Evan's reading; counts and side consistent, construction route unknown)
ry-xu's s(126) = 7.5 + 3 sqrt 2 = Kingbird s(85b) + 2.  s(85b) is a 6x6 diamond at 45 deg in uprights; T1 with an L of
2 rows + 2 columns gives 85 + 22 + 18 = 125.  The record's 126th comes from a 2-square tilted corner piece at UR
(53.3 / 64.9 deg) whose rectangle overlaps the upright-only corner of the s(85b) box; its diamond has the top diagonal
row slid by one lattice step.  Our pair screen requires a guillotine cut, so it cannot generate this (TODO).
