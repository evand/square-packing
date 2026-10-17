# Fill A: noisy greedy best fit (`fill_greedy.py`), 10-16

## Design

Per placement (m times):

1. Raster the covered region at h = 0.04 (squares + everything outside the box, so walls count as covered).
2. For each angle in a small palette (0 plus up to 5 clustered tilts of the squares within 1.6 of the uncovered
   region; `SHOVE_TILT=0` gives uprights only; `SHOVE_RANDANG=r` adds r random angles per candidate), FFT-convolve
   the covered raster with the square's footprint (shrunk by 1.5h so touching squares give no overlap pixels) and
   with a ring of width 0.08 around it: `score = overlap_area - LAM * contact_area` (LAM 0.05).
3. **Best of K random options** (Evan's knob): draw K random poses (centre = uniform uncovered pixel, angle = uniform
   from the palette), slide each to the local minimum of its angle's score map (repeated argmin over a 0.3 window),
   keep the lowest score.  K = 1 is random + slide; large K is near-greedy.
4. Continuous clean-up of the pixel pose: push it along the SAT separation vectors of its current overlaps
   (neighbours, walls) for up to 30 steps, keep the least-penetrating pose.  Paint it and repeat.

Overlaps are allowed (least-bad pose is taken when nothing fits); the harness's fq quench relaxes them.  `fill`
returns `SHOVE_C` = 3 candidates; `SHOVE_K` is a comma list cycled over candidates (default `64,16,16`).  Fill cost is
1-3 s per candidate; the quench dominates.

Lessons from building it:
* The first version scored with LAM 0.5: contact (covered pixels in a ring outside the square) grows with overlap,
  so it rewarded overlapping; nearly every placed square penetrated 0.1-0.5.  LAM 0.05 + local slide fixed the
  ranking; the remaining pixel-quantisation overlaps (~2h) were fixed by the SAT push (max pen 0.003-0.04 at n = 101).
* Geometry fact (worth knowing for every fill): for an LL block B that is a lattice down-set, (B + [0,1]^2) \ B is
  free after a d = 1 shove and consists of exactly w + h + 1 lattice cells (w, h: the block's width/height in
  cells); the two wall strips hold k - h and k - w uprights.  Total 2k + 1: an exact upright fill along the
  block's staircase always exists and reproduces the L-step side.  `meta['cells']` undercounts it because
  `strip_rects` uses the block's bounding box rather than its staircase.  So at d = 1 and m = 2k + 1 "matching the
  L-step" is the easy baseline; beating it needs the jumble to relax (s' < s + 1) or an extra square (m = 2k + 2 at
  s + 1).

## Metric note

bench `vsL < 0` counts any m, including m = 2k (n + 2k squares below s + 1), which is not a win over the L-step
(it adds 2k + 1).  My counts ("TRUE beat-L") use: m >= 2k + 1 and s' < s + 1, or m >= 2k + 2 and s' <= s + 1.
The bench's vsL does not see the second kind (vsL = 0).  Analysis script: scratchpad `greedy/an.py` (per-m
breakdown: non-grid, < L, = L, min vsL, min vsT).

Harness issue: hit files are named `{n}_{corner}_{n2}_{seed}.txt` with no candidate index, so several hits of the
same task and count overwrite each other (seen: two 210_LL_238_4 hits).  Suggest adding the candidate index.

## Test runs (all `python3 shove/bench.py run --bench runs/shove/b1005 --fill greedy --procs 2 --limit 6`, d 1 unless
said; load was 22-33 from other experiments throughout)

Task set te40 = `--tasks-every 40` (n 101 LL, 126 LL, 151 UL, 181 UR, 210 LL, 237 UL); te47 = `--tasks-every 47`
(n 101 LL, 128 LR, 170 UR, 202 LR, 232 UR, 264 UL).  18 candidates per m per run.  Columns: non-grid count at
m = 2k / 2k+1 / 2k+2; at 2k+1 the counts below / equal to the L-step; TRUE beat-L candidates (tasks); CPU.

| run (out dir `runs/shove/b1005/…`) | env | ng 2k/2k+1/2k+2 | 2k+1: <L / =L | TRUE beat-L | CPU h |
|---|---|---|---|---|---|
| greedy_dev1 (`--limit 3`) | K=16 | 5/8/2 of 9 | 0 / 3 | 0 | 0.014 |
| greedy_k4_d1 | K=4 | 3/3/3 | 0 / 1 | 0 | 0.079 |
| greedy_k16_d1 | K=16 | 7/10/2 | 2 / 3 | 2 (210 LL) | 0.045 |
| greedy_k64_d1 | K=64 | 6/7/3 | 0 / 5 | 0 | 0.037 |
| greedy_k16_d1_s1 (`--seed 1`) | K=16 | 6/7/2 | 0 / 3 | 0 | 0.040 |
| greedy_k16_d0.95 (`--d 0.95`) | K=16 | 8/8/3 | 0 / 3 | 0 | 0.046 |
| greedy_k16_d0.9 (`--d 0.9`) | K=16 | 5/3/5 | 0 / 2 | 0 | 0.055 |
| greedy_k16_d1_notilt | K=16 TILT=0 | 4/7/3 | 0 / 4 | 0 | 0.025 |
| greedy_k16_d1_lam02 | K=16 LAM=0.2 | 6/8/2 | 2 / 4 | 2 (210 LL) | 0.050 |
| greedy_k16_d1_ra2 | K=16 RANDANG=2 | 8/6/4 | 0 / 0 | 0 | 0.058 |
| greedy_k256_d1_te40 | K=256 | 6/9/5 | 0 / 6 | 0 | 0.044 |
| greedy_mix_d1_te40 | K=64,16,16 | 7/8/4 | 0 / 4 | 0 | 0.056 |
| greedy_k4_d1_te47 (te47) | K=4 | 7/7/7 | 0 / 1 | 0 | 0.080 |
| greedy_k16_d1_te47 (te47) | K=16 | 12/16/11 | 0 / 9 | 2 (232 UR) | 0.036 |
| greedy_k64_d1_te47 (te47) | K=64 | 12/17/9 | 0 / 17 | 5 (232 UR, 264 UL) | 0.033 |
| greedy_k256_d1_te47 (te47) | K=256 | 11/16/11 | 0 / 16 | 6 (232 UR, 264 UL) | 0.035 |

Notable candidates (none below the record as of the cutoff; no LIVE lines):
* 210 LL, m = 2k+1 = 29: n' 239 at s' 15.96708 (lam02) and 15.96865 (k16): below the L-step's s + 1 = 15.97293 by
  up to 5.9e-3 (record at cutoff 15.9538, vsT +1.3e-2).  The only s' < s + 1 at 2k+1; seen in 2 of the 4 K=16 runs
  that include the task (not with seed 1, not with the K list), so it is luck of the draw.
* 232 UR, m = 2k+2 = 32: n' 264 at s + 1 = 16.778175 exactly (one more square than the L-step at its side); equals
  the record as of the cutoff (vsT +3.5e-12), vsNow +3.4e-3.  K = 16, 64, 256.
* 264 UL, m = 2k+2 = 34: n' 298 at s + 1 = 17.778175, also equal to the cutoff record.  K = 64, 256.
* 151 UL, m = 2k+2: n' 177 at 13.823223, 3.5e-4 above s + 1 and 2.4e-4 above the cutoff record (most K >= 16 runs).

## K comparison (Evan's knob)

Same tasks and seeds.  K = 4 is clearly worst (3-7 non-grid of 18 per m; it rarely even reproduces the L-step:
1 / 18 at 2k+1) and costs the most CPU (more overlap, longer quenches).  K = 64 / 256 reproduce the L-step almost
always on te47 (17 / 18, 16 / 18) and find the +1-square cases; K = 16 gave the only s' < s + 1 result (210 LL) but
is less reliable.  CPU per run is similar for K >= 16 (0.033-0.056 h per 6 tasks); the fill itself is not the cost.
Decision for the final run: K list 64,16,16 (one near-greedy candidate + two noisier ones).

d < 1 (0.95, 0.9) at the same m never beat the L-step: the quench just expands back to s + 1 or collapses.

## Final run (10-17 00:15, load ~7 at launch, rose to ~18 from other jobs)

    python3 pk/jobs.py start greedy_full --procs 4 -- python3 shove/bench.py run --bench runs/shove/b1005 --fill greedy \
        --out runs/shove/b1005/greedy_full --tasks-every 1 --procs 4 --cpu-h 4          # default SHOVE_K=64,16,16

It **crashed after 205 of 293 tasks** (1845 rows): `hop.quench` -> `engine.stats` raised JSONDecodeError because fq
printed `"min_gap": inf` (and `"s": 322.47`, `s_coarse` 3.3e4: the SCREEN quench at loosen 1.05 diverged; loosen 1.0 and
1.02 were fine on the same input).  The exception kills the whole bench (`ex.map`).  I resumed tasks 205-292 with the
same seeds using a scratch driver that wraps `harness.evaluate` in try/except (no shared files edited):

    python3 pk/jobs.py start greedy_full_rest --procs 4 -- python3 <scratchpad>/greedy/rest.py 205 runs/shove/b1005/greedy_full_rest

One evaluation failed the same way there (saved: `runs/shove/b1005/greedy_full_rest/evalerr/270_3736038_469437648.npy`,
row 0 = [S, 0, 0], then the 270 poses; reproduces: loosen 1.05 -> s 322, min_gap inf).
CPU (manifests, measured): 5066 s + 4583 s = **2.68 h** (bench's own row sum says 2.37 h).

`bench.py score` (each dir, and the concatenation):

    greedy_full:      1845 rows (0 fill errors), 1845 feasible, 995 non-grid; beat L-step from its source: 82; below record
                      as of cutoff: 0; rediscovered: 0; below event-bench now: 0; vsL min -3.806e-02, median +3.9e-12
    greedy_full_rest: 792 rows, 791 feasible, 466 non-grid; beat L: 23; below cutoff record: 0; rediscovered 0; live 0
    combined:         2637 rows, 2636 feasible, 1461 non-grid; beat L: 105; 0 / 0 / 0

The 105 "beat L" are mostly m = 2k (90 rows; not a win, see Metric note).  Per m (combined, 879 rows each):

| m | non-grid | s' < s + 1 | s' = s + 1 | min vsT |
|---|---|---|---|---|
| 2k | 566 | 90 | 317 | +3.9e-12 |
| 2k + 1 | 599 | **15** | 348 | +2.4e-12 |
| 2k + 2 | 296 | 0 | **28** | +2.8e-12 |

**TRUE beat-L: 43 candidates on 31 tasks** (15 at m = 2k+1 with s' < s + 1, 11 tasks; 28 at m = 2k+2 with s' = s + 1,
20 tasks).  None is below the record as of the cutoff, none rediscovers a later record, no LIVE line.
* m = 2k+1, s' < s + 1 (hits in `runs/shove/b1005/greedy_full{,_rest}/hits/`): 105 UR -> 126 (vsL -1.55e-2, vsT
  +1.6e-2), 131 LL -> 154, 132 LL -> 155, 156 LL -> 181 (-1.2e-2), 182 LL -> 209, 205 UL -> 234 (-1.4e-5, vsT +1.6e-3),
  210 LL -> 239 (-7.4e-3, again), 228 UL / UR -> 259 (-9.8e-7), 235 UL -> 266 (-3.4e-3, **vsT +1.6e-4**), 266 UL -> 299
  (-2.6e-5, **vsT +1.6e-4**), 272 LL -> 305.  The closest to the cutoff record are 235 UL -> 266 (s' 16.8232233047)
  and 266 UL -> 299 (17.8230363819).
* m = 2k+2 at exactly s + 1 (one square more than the L-step at its side): sources 147 (all 4 corners) -> 173 at
  13.65685, 201 LR -> 231, 232 (4 corners) -> 264, 264 (4) -> 298, 290 (4) -> 326, 295 LR / UL -> 331.  All equal the
  record as of the cutoff (so the record already has this structure).  bench does not write hit files for these
  (vsL = 0).

By candidate slot (Evan's K): K = 64 slot (879 rows): 511 non-grid, 25 TRUE beat-L, 187 reproduce the L-step at 2k+1.
K = 16 slots (1758 rows, two per task): 950 non-grid, 18 TRUE beat-L, 161 at L.  K = 64 found most of the +1-square
cases (18 of 28); K = 16 found most of the large s' < s + 1 margins (8 of 15 rows, incl. 131, 132, 156, 182, 210, 272),
K = 64 the tiny ones and 105 / 235.

## Conclusions

* Noisy greedy best fit reliably does at least what the L-step does when the block is a staircase (it rebuilds the
  upright staircase band; 348 / 879 at 2k+1 match s + 1 exactly), unlike the example fill (9 / 9 grid collapse).
* It beats the L-step by side on 11 / 293 tasks (margins up to 0.0155) and by count on 20 / 293 tasks, but never the
  record as of the cutoff: the +1 cases tie it exactly, the best side cases miss by 1.6e-4.
* d < 1 does not help (dev runs): the quench re-expands.
* Very small K (4) is bad; K 16-256 similar in cost; mixing a near-greedy and noisy candidates covers both kinds of win.

## Harness issues found (not fixed; shared files)

1. An fq divergence (`inf` in its JSON stats) raises out of `evaluate` and kills the whole `bench run`; suggest
   try/except in `evaluate` or `one` (row with s None + error) and/or fq emitting valid JSON for inf.
2. Hit file names `{n}_{corner}_{n2}_{seed}.txt` lack the candidate index: hits of the same task / count overwrite.
3. `vsL < 0` counts m = 2k rows; and m = 2k+2 at s' = s + 1 (a real win by count) is not counted or saved.
4. `cpu` per row undercounts vs the manifest (2.37 h vs 2.68 h).

## Not tested

K = 1; softmax/top-k selection instead of best-of-K; seeds beyond 0/1; d > 1; the final run used d = 1 only; letting
existing jumble squares move in the fill (only added squares are placed); a ranking that targets the slack (e.g. a
tilted square in the wall-strip end gaps) for m = 2k+2 at s' < s + 1.
