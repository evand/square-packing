# Fill B: infeasible fill + expand / relax (`fill_relax.py`; 10-16 night)

Evan: "we can also fill with infeasible packings (least-bad among random options?) and then expand/relax from there."
Design note (Evan, via coordinator): make "best of k random options" an explicit knob (k = 1 pure random, large k
near greedy); compare a few k on the same tasks and seeds.

## Mechanism

1. **Least-bad infeasible placement.**  Optionally the clean upright cells (`meta['cells']`, zero overlap) go in
   first (`SHOVE_CELLS=1`, default).  The remaining squares go in one at a time; for each, draw `SHOVE_K` random poses
   and keep the one with the least *total* penetration (sum of SAT depths against every square present, old and new,
   plus wall excess; vectorised in `pen_vec`).  Positions: uncovered raster points (h = 0.1, re-rastered every 4
   insertions), 75 % from the strip guide dilated by 0.6, 25 % anywhere uncovered.  Angles: 50 % upright, 35 % a
   neighbour's tilt (within 2.5) + N(0, 2 deg), 15 % uniform.
2. **Expand / relax**, one fq quench per schedule from the same infeasible state (`SHOVE_RELAX`; keys of `RELAX`):
   `round.15c` (rounded squares r0 0.15, compress with side free, r -> 0 at the held box, then ALM + polish),
   `round.25c` (r0 0.25, 40 steps), `holdP1.01` (rounded r0 0.15 relaxation at the box held at 1.01 S, then SLP polish
   only, no ALM).  Also defined, not in the default: `L1.1`, `L1.2` (plain ALM at big loosen), `round.15cL1.1`,
   `shrinkhop`.
3. Return the relaxed feasible states, best internal side first (the bench evaluates up to 3, quenching again).
   With `SHOVE_NPLACE` > 1 the internal side selects the 3 of NPLACE x schedules that the bench evaluates.

Env knobs: `SHOVE_K` (default 64 after the k comparison below), `SHOVE_NPLACE` (independent placements per call,
each relaxed by every schedule; default 1), `SHOVE_RELAX`, `SHOVE_CELLS`, `SHOVE_RAW` (also return the unrelaxed state),
`SHOVE_LOG` (one json line per fill call: penetration, internal side per schedule).

## Dev tests (load ~22-29 from EX1 / FB2; 1-2 processes)

All at d = 1 unless noted; tasks are `tasks.json` indices 0, 40, 80 (n 101 LL blk 19, 126 LL blk 30, 151 UL blk 55),
m = 2k + 1; "eval" = `harness.evaluate` of the relaxed state.  Scripts in the session scratchpad (t2-t6.py; one-offs).

**T1: relax schedule screen, all m squares least-bad (K = 16, no cells).**  vsL of eval:

| schedule | 101 | 126 | 151 |
|---|---|---|---|
| plain ALM loosen 1.0 / 1.1 / 1.2 | grid 12 (+.46) all | grid 13 (+.23) all | grid 14 (+.18) / 14 / 14.09 |
| round r0 .15 compress | **+.187** non-grid | **+.072** non-grid | grid |
| round r0 .10 compress, hold 1.01 | +.207 non-grid | +.098 non-grid | grid |
| soft a0 .3 compress | grid | grid 13.07 | grid |
| shrink + 3 hops, loosen 1.1 | **+.085** non-grid | grid | grid |

Total penetration of the placement: 4.2 / 5.7 / 5.7.

**T2: clean cells first (16 / 13 / 8 cells), least-bad rest (K = 16).**  Penetration 1.8 / 2.2 / 5.6.
round .15c: +.241 / +.191 / grid; round .25c (40 steps): +.172 (grid line) / grid / +.134 (grid line);
round .15c at loosen 1.1: +.241 / +.190 / grid; shrink+hops: grid everywhere.

**T3: fixed-box relaxation (cells first, K = 16).**  hold 1.0 / 1.02 (rounded relax at the held box, then ALM):
101 +.17 (grid line) / +.29; 126 and 151 grid.  plain `--shrink` (grow to feasibility, shrink in fixed-side
relaxations): grid everywhere (one run at 151 diverged to s 157 before the eval quench repaired it).
Hold + SLP polish only (no ALM): hold 1.0: +.17 grid / grid; **hold 1.01: +.300 / +.164 non-grid**; r0 .05: grid.

**T4: one square at a time** (cells first, polish; then each remaining square least-bad + a quench).  With the hold
+ polish quench: n 101 stays at exactly s + 1 = 11.535534 for 3 insertions, then jumps to 12.41 / 12.01 (final eval:
grid 12); n 126 drifts 12.708 -> 12.624 (fewer squares than m at that point), final grid 13.  With a plain SCREEN quench
per insertion: collapses to the k + 1 grid at the first or fourth insertion.  So incremental insertion does not help:
whenever the remaining squares do not fit at the current side, the quench escapes to the k + 1 grid.

Side note: the shove alone (no fill) quenched by the evaluator also does not return to s (smoke_none: 10.69 with a full
line, 12, 13): ALM from a state with a unit-wide gap loses the jumble's structure; any fill has to leave the state
nearly tight before the evaluator sees it.

## Bench runs (`bench.py run ... --limit 6 --tasks-every 40 --procs 2`, seed 0; tasks n 101, 126, 151, 181, 210, 237)

Default schedules (round.15c, round.25c, holdP1.01), cells first, d = 1.  54 rows each (6 tasks x 3 m x 3 cands).

Command (per row; `SHOVE_K` etc. as listed, `SHOVE_LOG=<out>/fill_log.jsonl`):

    SHOVE_K=64 python3 shove/bench.py run --bench runs/shove/b1005 --fill relax --out runs/shove/b1005/relax_k64 \
        --limit 6 --tasks-every 40 --procs 2

| run (`runs/shove/b1005/…`) | knobs | non-grid / 54 | beat L | best vsL | median vsL (non-grid) | CPU s | n with non-grid |
|---|---|---|---|---|---|---|---|
| relax_k4 | K 4 | 13 | 0 | +.0607 | +.171 | 635 | 101 126 151 |
| relax_k16 | K 16 | 12 | 0 | +.0649 | +.166 | 542 | 101 126 151 |
| relax_k64 | K 64 | **18** | **1** | **-.0005** | +.119 | 485 | 101 126 151 210 237 |
| relax_k256 | K 256 | 14 | 0 | +.0003 | +.139 | 479 | 101 126 151 210 237 |
| relax_k64_nocells | K 64, no clean cells | 17 | 0 | +.0119 | +.121 | 576 | 101 126 151 |
| relax_k64_d095 | K 64, `--d 0.95` | 19 | 0 | +.0132 | +.150 | 553 | 101 126 151 210 |
| relax_k4_s1 | K 4, `--seed 1` | 9 | 0 | +.1159 | +.144 | 474 | 101 126 151 |
| relax_k64_s1 | K 64, `--seed 1` | 14 | 0 | +.0693 | +.162 | 432 | 101 126 151 |
| relax_k64_np3 | K 64, 3 placements x (round.15c, holdP1.01) | **20** | 0 | +.0235 | +.119 | 522 | 101 126 151 181 |

Mean total penetration of the placement (fill_log): K 4: 9.0, K 16: 4.9, K 64: 3.3 (K 256 lower still); the k knob
does what it should for the placement.  Evaluated yield, two seeds each for K 4 and 64: K 4 13 + 9 = 22 / 108
non-grid, best vsL +0.061; K 64 18 + 14 = 32 / 108, best -0.0005 (K 16: 12 / 54, K 256: 14 / 54, one seed).  So
large K gives somewhat more non-grid results (weak evidence: 22 vs 32 of 108 on the same 6 tasks), and no k beats the
L-step at m >= 2k + 1.  CPU is dominated by fq quenches (~0.02 CPU-h per task), not by K: K = 256 costs no more than K = 4
in measured CPU (load varied 22-29 during these runs, so CPU-s comparisons are rough).

The one "beat": `relax_k64/hits/237_UL_267_5.txt`, n 237 UL (block 90), **m = 30 = 2k** (one square fewer than the
L-step) -> n' 267 at 16.9106940074, s + 1 = 16.9111917 (vsL -5.0e-4), but the cutoff record for 267 is 16.8467
(vsT +0.064): not a useful packing, just the L-step minus one square slightly relaxed.  The internal fill sides
(fill_log) show the same picture: best internal vsL -0.061 at n 101 m 20 (n' 121, which the evaluator then collapses
to the 11 grid, the record), otherwise >= 0.

Note: during the session `bench.py score` changed its L-step criterion to "beatL": vsL < 0 **and** m >= 2k + 1 (more
squares than the L-step below s + 1).  Re-scored with it, every dev run above has 0 beats (the relax_k64 row was
m = 2k).

## Final run (00:15 10-17, load ~6 before launch)

    python3 pk/jobs.py start relax_full --procs 4 -- env SHOVE_K=64 SHOVE_NPLACE=3 \
        SHOVE_LOG=runs/shove/b1005/relax_full/fill_log.jsonl bash -c "python3 shove/bench.py run --bench runs/shove/b1005 \
        --fill relax --out runs/shove/b1005/relax_full --tasks-every 3 --procs 4 --cpu-h 4 > runs/shove/b1005/relax_full/stdout.txt 2>&1"

Config: K 64, 3 placements per (task, m), each relaxed by round.15c, round.25c, holdP1.01 (9 fq quenches per m), the 3
best internal sides evaluated.  d = 1.  The CPU budget stopped result collection at 83 of 98 tasks (747 rows); the
pool kept running the queued tasks after "cpu budget reached" (see harness notes), so I stopped the job by hand
(`jobs.py stop relax_full`) and scored it:

    runs/shove/b1005/relax_full: 747 rows (0 fill errors), 747 feasible, 476 non-grid; CPU 4.01 h
       beat L-step from its source (more squares below s + 1, or 2 more at s + 1): 1 on 1 tasks   below record as of cutoff: 2   rediscovered later records: 0   below event-bench now: 0
       vsL quantiles (s' - (s + 1)): min -4.908e-02  10% +3.457e-12  median +1.128e-01
         n 228 UR blk 4 m 30 -> n' 258 s' 16.5641504424  vsL -4.05e-02  vsT -6.92e-03  vsNow +7.02e-04
         (x2), then n 145 UR / LL m 25 -> 170 at 13.5355339059 (= s + 1, vsT +2.8e-12), 173 LL -> 200 14.6568542495, ...

By m (rows / non-grid / vsL < 0 / vsL < 0.01 / best vsL):

| m | rows | non-grid | vsL < 0 | vsL < 0.01 | best vsL |
|---|---|---|---|---|---|
| 2k | 249 | 152 | 19 (6 tasks) | 44 | -0.0491 |
| 2k + 1 | 249 | 178 | **1** | 23 | -0.0116 |
| 2k + 2 | 249 | 146 | 0 | 0 | +0.0256 |

Non-grid fraction by n: 100-149 0.68, 150-199 0.61, 200-249 0.56, 250-300 0.74.  39 non-grid rows reproduce the
L-step side exactly (|vsL| < 1e-8).  Internal winners (fill_log): round.15c 110, round.25c 94, holdP1.01 49.

**Results worth a look (no LIVE, no rediscovery; nothing posted):**
* **Beats the L-step at its own count**: n 126 UL (block 9), m = 23 = 2k + 1 -> n' 149 at **12.76318065** vs
  s + 1 = 12.77473513 (vsL -0.0116), but 0.056 above the cutoff record for 149 (12.7071).  Hit:
  `runs/shove/b1005/relax_full/hits/126_UL_149_14.txt`.
* **Below the cutoff record**: n 228 UR (block 4), m = 30 = 2k -> n' 258 at **16.5641504424** vs frontier_at(258, reg1005)
  16.5710678 (vsT -6.9e-3); above today's event-bench frontier 16.5634480 by 7.0e-4 (so not live, not a rediscovery of
  the exact later record).  Two candidates of the same call reached it (and a grid-flagged third at 16.5668).  Hit:
  `runs/shove/b1005/relax_full/hits/228_UR_258_58.txt`.  With block 4, the shove here is nearly "move the whole
  packing to the lower-left and open both wall strips": an L-step with one fewer square, where the jumble relaxes
  into the strips.
* m = 2k below s + 1 (one fewer square than the L-step, below its side): 19 rows on 6 tasks (102 UR, 103 UL,
  105 UR, 123 UL, 126 UL, 150 LL/UR, 228 UR); only 228 is below the cutoff record.

CPU: dev tests ~0.3 CPU-h in-process scripts (unmeasured, estimate) + dev bench runs 1.03 CPU-h (sum of the table's
CPU s) + final 4.01 CPU-h recorded (plus a few CPU-minutes of discarded queued tasks before the manual stop).

## Conclusions

* Infeasible least-bad fill + relax does not beat the L-step in general: at m = 2k + 1 one task of 83 (0.4 % of rows)
  goes below s + 1, by 0.012, and that one is far above the record; at m = 2k + 2 never.  Most non-grid results sit
  0.05-0.2 above s + 1.
* The relax schedule matters more than the placement: plain ALM from the overlapping state (any loosen 1.0-1.2),
  `--shrink`, soft squares and one-at-a-time insertion all collapse to the k + 1 grid on the dev tasks; the rounded
  pre-stages (round_r0 0.15 / 0.25 with compress, or a 1.01-held relaxation + SLP polish) keep non-grid structure.
* k (best of k random poses): penetration falls with k (9.0 / 4.9 / 3.3 at k 4 / 16 / 64); evaluated non-grid yield
  rises weakly (k 4: 22 / 108, k 64: 32 / 108 over two seeds); no k beat the L-step at m >= 2k + 1 in the dev runs.
  CPU is quench-dominated, so large k is free.
* The useful by-product is the m = 2k case: one fewer than the L-step, the jumble relaxes into the strip enough to
  land 0.04-0.05 below s + 1, and once (n 228 -> 258) below the cutoff record.  This is "L-step minus one + relax",
  so it is worth checking against the other fills / the explorer before reading anything into it.

## Harness notes (not changed; for whoever owns bench.py)

1. Hit files are named `{n}_{corner}_{n2}_{seed}.txt`, so several hit candidates of one task at the same n' overwrite
   each other (the last written is kept).  Add the candidate index.
2. `--cpu-h`: the loop breaks, but `ProcessPoolExecutor.map` has already submitted every task, and leaving the `with`
   block waits for all of them (results discarded, `man.finish` and score delayed until they finish).  Use
   `ex.shutdown(cancel_futures=True)` or submit lazily.
3. `evaluate` quenches in the shove's box S, not the candidate's own side; fine for this fill (the relaxed states sit
   inside [0, S]), but a fill returning a state that needs a bigger box would get repaired by fq.

## Not tested

d other than 1 and 0.95 (and d = 0.95 only on the 6 dev tasks); kicks / round_kick; schedules mixing an ALM stage at
1.05-1.2 after the rounded stage; a polish (POLISH_DEEP) of the near-miss rows; anything beyond the 83 tasks of the
final run (the last 15 tasks, n ~ 280-300, were cut by the CPU budget); the n' = 258 packing has not been checked
by a verifier beyond fq's own feasibility check, and nothing was added to the store.
