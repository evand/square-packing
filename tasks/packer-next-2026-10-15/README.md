# Packer: open ideas and next steps (end of session 10-15)

State at the pause (Evan 10-15: "pause and evaluate" after L2e / GX1).  Results: `search/packer/PACKER.md` entries
of 10-15 (L2d, L2e, EV1, EV2 gen 0 / 1, FP1, GX1, s(110) census 19:12); EV2 detail `search/packer/evolve2/LOG.md`.
Nothing is queued or running.  Ranked by my estimate of value; each needs Evan's go before compute.

## 1. Live slack-composition hunt (Mondrian v1 blind, n <= 324)
* Why: L2e (blind #481 time travel, exact ranking `l2_screen.py --rank exact`) beat the pre-#481 record at 4 / 15 n
  (153, 232, 236, 302) -- the pre-registered ">= 3 -> live hunt next" rule fired.  L1 (named pieces) was 4 / 5.
* Blocker: cost.  The compose stage (keep 1000 x 4 sides = up to 8,000 joint fills + fq per n) took ~6 CPU-h per n;
  the explorer after it 0.8 CPU-h per n.  Need a cheaper screen first: e.g. keep 200, 2 sides, or a fill-count
  pre-screen with `fillfast` (58 ms / fill) before fq; measure on the L2e set that the 4 hits survive the cut.
* Then: library from the live store (`analysis/croplib.py`), all n with a sub-grid frontier, explorer 45 min on the
  best compositions; pre-register record rate.

## 2. EV2 gen 2 (Claude Code subagent prompt evolution, J-style generator)
* Gen 1 (holdout row 234-240): blind ties at 234 / 235, best 236 +0.0445, 237 +0.101, 238-240 never below the grid.
* Changes for gen 2: focus 236-240; give agents a staircase-block / interlocking-block helper (by-eye block starts
  never interlocked: fill 15-20 short); a cheaper quench budget in the evaluator (seeds, loosens, candidate cap, procs
  per agent -- each a quality trade to test; gated quench rejected, FP1); carry gen-1 lessons (LOG.md "Gen 1 summary").
* Keep the 234-240 records unrendered / unseen (possible future holdout).

## 3. Explorer seeded from EV2 gen 1's 55 novel near-record minima
* ~5 CPU-h; P(record) ~3 % (my estimate).  Minima indexed in `runs/evolve2/minima.jsonl` (corner signatures).
  Note 68 / 79 near-record minima are at 234 (degenerate variants): seed by distinct corner cell, not by count.

## 4. s(110): census corner-cell holdouts
* The corner signature (`analysis/cornersig.py`) splits the census into cells; below the grid ~97 % saturated (Chao1).
  Test: withhold whole corner cells from a census snapshot, give the searcher the rest, measure rediscovery per CPU-h
  by method (structural reach, not value).  Census now 1,078 minima (19:12 snapshot).
* GX1 lesson for this: descent from foreign above-grid structure crosses easily (45 / 47) but collapses to the record
  basin; features survive only when descent stops early (melt arms), ~2e-3 above the record.  Protected melt = random.
  A "stop early / constrained descent" mode would be needed to keep a feature to a sub-grid minimum.

## 5. Smaller / parked
* J generator proper (constructive, beyond EV2's agent programs) -- see TODO "Constructors".
* numba for fillfast: asked, not approved; not needed now (fill 58 ms; the quench dominates EV2 evals).
* EV1 (`claude -p` loop): superseded by EV2; keep code, no further runs.
* L2e follow-up: 236 was beaten by only 2.9e-4 over the ante record -- check whether its composition is the posted
  structure family (cheap, analysis only).
