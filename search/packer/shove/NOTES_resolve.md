# Fill C: local re-solve (`shove/fill_resolve.py`, 10-16)

## Design (as it evolved)

1. **refill:r** (first try): remove every non-block square within r of the strip guide, greedily re-place band + m
   squares (clean max-contact poses, angles = 0, 45 and the band's own), least-penetration for the rest.
2. **keep:r**: keep the band in place, add the m squares (clean block-lattice cells first, then greedy), then relax
   **only band + new** (rest fixed) by overlap-energy minimisation (L-BFGS on sum pen^2 + wall^2, finite-difference
   gradient of each square's own energy, vectorised SAT).  An earlier Jacobi "push" relaxer (`local_relax`) was
   kept for reference: it stalls at the same residuals.
3. **squeeze:r**: keep:r built at d = 1, then re-shove the source at d' < 1 with the same added squares carried in the
   block's frame (so they overlap the jumble by about 1 - d'), and relax band + new in the box s0 + d'; bisection on d'
   for the smallest d' whose local relaxation reaches max penetration <= 1e-6.  This is the "re-solve" proper: the
   band has to absorb 1 - d' of misfit locally.

## Region (important fix)

The strip guide (`meta['strip']`) is the block's bounding box grown by d, not its staircase.  With a staircase block
the real opened gap runs along the steps, often 2-3 units from the guide (n = 126 LL: the gap along the block's top
steps at y = 3-4 is 2.5 from the guide).  So the region is: points within r + 1 of the block (Chebyshev distance to
the union of block squares) or r + 0.8 of the guide; the band is non-block squares within r + 0.5 of the block or r of
the guide.  Lattice cells on the block's own lattice, nearest the block first, then fill the gap: on most tasks they
reach m = 2k + 1 cleanly (the staircase gap under a diagonal (1, 1) shift holds one cell per unit step).

## Test runs

All dev runs single process, `uptime` load ~22-25 (other experiments), d = 1 unless noted.  Scripts in the session
scratchpad (t1/t2/t3.py: one task, one variant, then `harness.evaluate`; t2 also tries `fq --shrink` + polish).

**T1 refill:2, task 0 (n 101 LL, the 45-degree diamond packing), m = 21.**  Band 38 removed; the greedy re-places only
33 of 59 cleanly (worse than the jumble it replaced); evaluate 11.8229, grid.  Removing the band and re-placing it from
scratch throws away the jumble's own good arrangement: dropped.

**T2 keep:r (guide-based region, before the region fix), m = 2k + 1, r = 1 / 2 / 3; Jacobi push relaxer 600 / 3000
iterations, then L-BFGS relaxer.**

| task | clean / m | residual pen after local relax | evaluate s' | vs s + 1 |
|---|---|---|---|---|
| 101 LL | 17-18 / 21 | 0.015-0.078 (more iterations: no change) | 11.5962-11.7711, non-grid | +0.06 to +0.24 |
| 126 LL | 17 / 23 | 0.038-0.058 | 13.0, grid | +0.23 |
| 103 UR | 15-16 / 21 | 0.054-0.19 | 11.786-12.0 (r 1, 3 grid) | +0.08 to +0.30 |
| 108 LR | 10-12 / 21 | 0.028-0.14 | 12.0-12.035, grid | +0.07 |

(`fq --shrink` + POLISH from the same states: equal or worse.)  The relaxers converge to jammed local minima; the
picture (n = 126 LL) showed why: the guide rectangles miss the staircase gap, so the extra squares were crammed into
the wrong place.  -> the region fix above.

**T3 keep:r with the region fix + block-lattice cells, m = 2k + 1.**

| task | lattice cells / m | residual | evaluate s' | vs s + 1 | fT(n + m) |
|---|---|---|---|---|---|
| 101 LL | 21 / 21 | 0 | 11.535534 | 0 (ties the L-step) | 11.5355 |
| 126 LL | 23 / 23 | 0 | 12.774917 | +1.8e-4 | 12.7071 |
| 103 UR | 11 / 21 | 0.074 (r 1), 0.006 (r 2) | 11.704439 / 11.709234 | +9e-4 / +5.7e-3 | 11.6569 |
| 108 LR | 21 / 21 | 0 | **11.916157** | **-9.8e-3** | 11.8813 |

Grid collapse gone (0 / 8).  108 LR beats the L-step from its source (not the record: fT is 0.035 lower).

**T4 squeeze:1 / squeeze:2 (bisection on d' in [0.9, 1], tol 1e-6) on 108 LR, 126 LL, 103 UR.**  No d' < 1 ever
relaxed to zero residual; results identical to keep.  Direct check on 108 LR, r = 2, L-BFGS 10 x 1000: residual
1.3e-4 at d' = 0.999, 6.8e-4 at 0.995, 2.5e-3 at 0.98, 1.0e-2 at 0.95 (roughly 0.13-0.2 x the misfit: the band does
not absorb it), and evaluate from each lands in the same 11.916157 basin.  The band (r <= 2, rest fixed) is locally
rigid: the jam runs through fixed squares, so the strict local re-solve cannot take misfit; the global fq quench is
what finds the 0.0098.  squeeze is kept in the file but out of the default set (5-10x the CPU, no gain).

**try1 (bench): `RESOLVE_VARIANTS=keep:1,kick:1,tilt:1` (default at the time), `python3 shove/bench.py run --bench
runs/shove/b1005 --fill resolve --out runs/shove/b1005/resolve_try1 --limit 6 --tasks-every 40 --procs 2`**
(log of per-variant clean counts / residuals: `runs/shove/b1005/resolve_try1.log`).  kick = keep, then gaussian kick
sigma 0.03 (angle 0.6 deg) on band + new and local relax again; tilt = keep without lattice cells.

    54 rows (0 fill errors), 54 feasible, 19 non-grid; CPU 0.14 h
    beat L-step from its source: 7   below record as of cutoff: 0   rediscovered: 0   below event-bench now: 0
    vsL quantiles: min -1.014e-02  10% -8.516e-03  median -4.4e-13

* By variant: tilt (no lattice cells) is grid or worse on 17 / 18 (one non-grid at +0.21): dropped.  keep and kick
  mostly land in the same basin; kick differed on 210 LL m 28 (-0.0101 vs -0.0076), 210 LL m 29 (-0.0085 vs -0.0064)
  and 237 UL m 30 (kick -0.0073 non-grid, keep grid).
* By m: of the 7 L-step beats, 2 are at m = 2k + 1 (n 210 LL -> 239: -0.0085 kick, -0.0064 keep; record as of
  cutoff 0.0107 lower), 5 at m = 2k (easy: the frontier at n + 2k is usually well below s + 1).  m = 2k + 2: never
  non-grid below s + 1 (the lattice holds only 2k + 1 cleanly; the extra one jams).
* Ties: 101 LL m 21, 181 UR m 27, 151 UL m 24 land at exactly s + 1 (the L-step's value through a different route).
* Lattice deficit: 237 UL (diagonal 45-degree band against a triangular block) gets only 17-18 clean of 30-32: the
  opened gap is narrow along the diagonal and the band's rows don't continue into it.  Variant ext (block-lattice
  cells, then cells continuing the band's own local lattices, +-1 along each band square's frame axes) adds only 1 on
  that task (19 / 31).

**try2 (bench): `RESOLVE_VARIANTS=keep:1,kick:2,ext:1`, same command, out `runs/shove/b1005/resolve_try2`.**

    54 rows (0 fill errors), 54 feasible, 24 non-grid; CPU 0.08 h
    beat L-step from its source: 8   below record as of cutoff: 0   rediscovered: 0   below event-bench now: 0
    vsL quantiles: min -8.698e-03  10% -7.564e-03  median -4.4e-13

Per (task, m), candidates keep / kick / third (G = grid):

| task, m | try1 keep / kick:1 / tilt:1 | try2 keep / kick:2 / ext:1 |
|---|---|---|
| 181 UR 26 (2k) | -0.0032 / -0.0032 / G | -0.0032 / -0.0032 / -0.0032 |
| 181 UR 27 | 0 / 0 / G | 0 / 0 / 0 |
| 210 LL 28 (2k) | -0.0076 / -0.0101 / G | -0.0076 / -0.0087 / -0.0076 |
| 210 LL 29 | -0.0064 / -0.0085 / G | -0.0065 / G / -0.0064 |
| 237 UL 30 (2k) | G / -0.0073 / G | G / G / G |
| 126 LL 22-23 | ~0 / +2e-4 / G | ~0 / +4e-4 / ~0 |
| 101, 151, all m = 2k + 2 | ties or G | ties or G |

ext lands in keep's basin whenever the block lattice already supplies the cells (no diversity); kick with r = 2 is a
bigger perturbation and goes grid as often as it helps.  Best set so far: keep:1 + kick:1.

**try3 (bench, d = 1.03): `RESOLVE_VARIANTS=keep:1,kick:1,ext:2 ... --out runs/shove/b1005/resolve_try3_d103 --d 1.03`**

    54 rows, 54 feasible, 20 non-grid; CPU 0.05 h; beat L-step 6 (all m = 2k: 181 UR -0.0035, 210 LL -0.0100); 0 below cutoff record

Slack in the box hurts: 210 LL m 29 (the only m = 2k + 1 beat at d = 1) goes grid on all three; 151 UL m 24 ties ->
grid.  d stays 1.

**try4 (bench): screened default `keep:1` + 4 x `kick:1`, one fq SCREEN (loosen 1.0) each, best 3 distinct
screened sides returned; out `runs/shove/b1005/resolve_try4`, same 6 tasks.**

    54 rows, 54 feasible, 23 non-grid; CPU 0.12 h
    beat L-step from its source: 8   below record as of cutoff: 0   rediscovered: 0   below event-bench now: 0
    vsL min -1.043e-02 (210 LL m 28)

Gains: 210 LL m 28 -0.0104 (best so far), 237 UL m 30 non-grid twice (-0.0008, -0.0073), 126 LL m 22 / 23 one tie
each instead of +2e-4.  Loss: 210 LL m 29 grid on all three.  The log shows why: m 29 needs one square beyond the 28
lattice cells, and this time every greedy placement of it relaxed to residual ~0.0028 (try1 / try2 keep: 5.7e-6, a
different hole: the rng stream differs).  So the m = 2k + 1 beat hinges on where the one extra square goes, not on
the kick.  -> v_keep now restarts the greedy + local relax up to RESOLVE_TRIES = 6 times when the clean placement
falls short (restarts pick each overlapping square at random among the 30 least-penetrating samples) and keeps the
lowest residual.

**Ablation (t3.py, one task each, m = 2k + 1): keep:1 with restarts vs norelax:1 (keep without the local relax).**

| task | keep:1 (restarts) | norelax:1 |
|---|---|---|
| 210 LL m 29 | resid 1e-5, s' 15.966514 (vsL -0.0065), 9.7 s | 16.0, grid |
| 103 UR m 21 | resid 0.0107, s' 11.724901 (+0.0214), 21.9 s | 12.0, grid |

The local relaxation of band + new (rest fixed) is what keeps the quench off the grid when the lattice cells fall
short.  Residual is not a perfect predictor (103 UR: residual 0.074 earlier gave 11.7044, 0.0107 now gives 11.7249).

**try5 (bench): try4's screened default + greedy restarts (RESOLVE_TRIES 6); out `runs/shove/b1005/resolve_try5`.**

    54 rows (0 fill errors), 54 feasible, 28 non-grid; CPU 0.36 h
    beat L-step from its source: 11   below record as of cutoff: 0   rediscovered: 0   below event-bench now: 0
    vsL quantiles: min -1.043e-02  10% -8.870e-03  median -4.5e-13

| task, m | try5 (3 returned candidates, vsL) |
|---|---|
| 210 LL 28 (2k) | -0.0101 / -0.0104 / -0.0076 |
| 210 LL 29 (2k + 1) | **-0.0089 / -0.0086 / -0.0085** (s' 15.964131; cutoff record 15.953819, 0.0103 lower) |
| 237 UL 30 (2k) | -0.0051 / -0.0009 / +0.0386 |
| 237 UL 31 (2k + 1) | +0.0561 / G / G |
| 181 UR 26 (2k) | -0.0032 x 3 |
| 181 UR 27, 101 LL 21, 151 UL 24 | exact ties with s + 1 |
| 126 LL 22 / 23 | ties / +2e-4 |
| every other m = 2k + 2 | G |

CPU tripled (0.36 h vs 0.12 h for 6 tasks) because the restarts also run where the clean placement can never reach
m (all of m = 2k + 2).  For the final run restarts drop to 2 when m > 2k + 1 (no m = 2k + 2 candidate has left the
grid in any dev run).

## Final run

Load 9.4 at 00:12, launched 00:13:

    python3 pk/jobs.py start resolve_full --procs 4 -- python3 shove/bench.py run --bench runs/shove/b1005 \
        --fill resolve --out runs/shove/b1005/resolve_full --tasks-every 4 --procs 4 --cpu-h 4

(--tasks-every 4, 74 tasks, rather than 3: at try5's ~0.05 CPU-h per task, 98 tasks would hit the 4 CPU-h cap and
drop the largest n, since the bench processes tasks in n order.)
Finished 00:59 (job log `runs/jobs/resolve_full/log`; no "** LIVE?" line).  `python3 shove/bench.py score
runs/shove/b1005/resolve_full`:

    runs/shove/b1005/resolve_full: 666 rows (0 fill errors), 666 feasible, 418 non-grid; CPU 3.08 h
       beat L-step from its source: 65   below record as of cutoff: 0   rediscovered later records: 0   below event-bench now: 0
       vsL quantiles (s' - (s + 1)): min -3.196e-02  10% -1.240e-03  median +3.594e-12
       (the 8 lowest-vsT rows are all exact ties with s + 1 where the cutoff record is the L-step: 101 LL -> 122,
        104 UL -> 125, 122 LL -> 145)

By m (74 tasks, 3 candidates per m):

| m | rows | non-grid | beat L-step (rows / tasks) | ties with s + 1 (|vsL| <= 1e-9) |
|---|---|---|---|---|
| 2k | 222 | 173 | 45 / 17 | 118 |
| **2k + 1** | 222 | 191 | **20 / 10** | 148 |
| 2k + 2 | 222 | 54 | 0 / 0 | 15 |

Best m = 2k + 1 per task that beats the L-step (hit files in `runs/shove/b1005/resolve_full/hits/`):

| source | n' | s' | vsL | vsT (cutoff record) | vsNow | file |
|---|---|---|---|---|---|---|
| 267 UL | 300 | 17.8294685 | **-0.01720** | +0.00535 | +0.00535 | 267_UL_300_61.txt |
| 108 LR | 129 | 11.9161570 | -0.00976 | +0.03485 | +0.04614 | 108_LR_129_5.txt |
| 210 LL | 239 | 15.9642267 | -0.00877 | +0.01041 | +0.01578 | 210_LL_239_40.txt |
| 273 UR | 306 | 17.9799828 | -0.00394 | +0.01654 | +0.01655 | 273_UR_306_64.txt |
| 235 UL | 266 | 16.8232233 | -0.00338 | **+0.00016** | +0.00019 | 235_UL_266_49.txt |
| 239 UR | 270 | 16.9524195 | -0.00140 | +0.01461 | +0.02422 | 239_UR_270_51.txt |
| 237 UL | 268 | 16.9099170 | -0.00127 | +0.03110 | +0.03110 | 237_UL_268_50.txt |
| 300 UL | 335 | 18.8229331 | -0.00119 | **+0.00006** | +0.00006 | 300_UL_335_73.txt |
| 208 LL | 237 | 15.9262159 | -0.00032 | +0.01502 | +0.02361 | 208_LL_237_39.txt |
| 228 LR | 259 | 16.6046015 | -9.8e-7 (marginal) | +0.00203 | +0.01322 | 228_LR_259_43.txt |

Closest to a record: 300 -> 334 (m = 2k) within 1e-6 of the cutoff record (tie), 300 UL -> 335 at +6e-5, 235 UL -> 266
at +1.6e-4.  Largest L-step beat overall: 108 LR m = 2k -> 128 at -0.032 (cutoff record 0.069 lower).

## Conclusions

* Against the baselines: fill_example collapses to the grid on 9 / 9; this fill keeps 418 / 666 non-grid (191 / 222
  at m = 2k + 1).  The L-step itself gives exactly s + 1; this fill ties it on 148 m = 2k + 1 rows and beats it on 20
  rows / 10 of 74 tasks (best -0.0172), never below the record as of the cutoff, no rediscovery, no live candidate.
* What did the work: (1) measuring the region from the block's staircase, not the guide rectangles, so the
  block-lattice cells fill the opened gap (usually all 2k + 1 cleanly); (2) for the squares beyond the clean cells,
  least-penetration insertion + L-BFGS relaxation of band + new with everything else fixed (ablation: without it,
  grid on both test tasks), with restarts; (3) kicks of the band + one fq screen to pick 3 distinct candidates.
* What did not: re-placing the band from scratch (refill) loses squares; local squeezing (d' < 1 with the band
  relaxed, rest fixed) never absorbs any misfit (the jam goes through the fixed squares); d = 1.03 slack is worse;
  no-lattice tilted fills go to the grid; extending the band's own lattices adds ~1 cell.
* m = 2k + 2 (one square beyond the L-step at the same side) never beat s + 1: 15 ties, 0 beats.

## Not tested

Other d values besides 1.0 and 1.03 (e.g. d slightly < 1 through the bench's own shove), r other than 1-3, more than
5 kick variants per m, the 3 / 4 tasks the final run skipped (--tasks-every 4), longer polish (POLISH_DEEP) of the
near-record hits (300 -> 335, 235 -> 266, 300 -> 334), and exact certification of anything.

## Harness notes

* `meta['strip']` is the block's bounding box grown by d; for staircase blocks the real gap is elsewhere.  A
  staircase-following strip (or the block-distance region used here) in the harness would help every fill.
* The bench shares one rng across the three m values of a task, so changing one variant's rng use changes the others'
  outcomes (seen on 210 LL m 29, try4 vs try1); per-m seeding would make A/B comparisons cleaner.
* bench.py `score` counts m = 2k beats together with m = 2k + 1 beats; a per-m breakdown would be more telling
  (most of the "beat L-step" count is m = 2k, which the frontier usually beats by more).
