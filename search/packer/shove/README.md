# Corner shove (CS1, 10-16; Evan's idea)

**Operator.** Take a packing (n squares, side s).  Pick a corner's grid-aligned block (`harness.corner_blocks`, from
`analysis/cornersig.py`).  Grow the box by d (default 1) and move the block with its corner while everything else moves
with the opposite corner.  That opens a strip of width d between the block and the jumble, plus along the walls beyond
the block: Friedman's L-step strip, routed through the interior where relaxation can happen.  A **fill** adds m squares
(the L-step adds m = 2k + 1 at side s + 1, k = floor s).  **Evaluate** quenches with the explorer's fq stages (SCREEN at
loosen 1.0 / 1.02 / 1.05 from the overlapping state, then POLISH) and checks grid collapse.

**What we already know (read before designing a fill).**
* The opened region follows the block's staircase, not its bounding box: clean upright lattice cells in it
  (`meta['cells']`) are well short of 2k + 1 (smoke test: 16 / 13 / 8 against 21 / 23 / 25 at n = 101 / 126 / 151).
  The rest must come from tilted squares and from the jumble relaxing into the gap.
* The naive fill (`fill_example.py`: the clean cells + least-overlap random poses) collapses to the integer grid on
  9 / 9 smoke candidates.
* From the 10-05 count-pair move tests (PACKER.md, "Count-pair moves"): an exact upright refill around a slightly
  moved band fails (it needs the side to grow, then collapses to the grid when squeezed); keeping the squares,
  inserting extras at least-overlap positions and penalty-squeezing from the overlapping state worked (106 / 300 below
  11 at n = 110).  `transplant.py` (tilted part + exact upright refill in a bigger box) never gained a square.

**Files.**
* `harness.py`: `corner_blocks`, `shove`, `strip_rects`, `clean_cells`, `free_mask`, `pen_one`, `half_extent`,
  `evaluate`, `lstep`.  The fill interface is in its docstring.
* `bench.py`: `build` (tasks from the as-of-cutoff record set; built: `runs/shove/b1005`, cutoff reg1005, n 100-300,
  293 tasks over 105 n, blocks >= 4), `run --fill NAME` (each task at m = 2k, 2k + 1, 2k + 2; up to --max-cands
  candidates per m), `score`.  Built-in fills: `lstep`, `uprights`, `none`.
* 10-20: the engines of `fill_greedy.py` and `fill_relax.py` moved to `pk/fill.py` (placers `bestofK`, `leastbad`); the files here are adapters with the same env knobs and the same output for a seed (checked).  Relax schedules: `pk/schedules.py`.  `free_mask` / `half_extent` / `_inside` now come from `pk.fill`.
* `fill_<name>.py`: one module per mechanism.  `fill_example.py` is the interface check.

**Metrics** (per evaluated candidate, non-grid only): beats the L-step from its source (s' < s + 1 at the same
source), below the record as of the cutoff, rediscovers a later record, below the event-bench frontier now (live
candidate: live-check with `candidate.py --live`, then `candidate.py`).

**Commands.**

    python3 shove/bench.py run --bench runs/shove/b1005 --fill NAME --out runs/shove/b1005/NAME_try1 --limit 6 --tasks-every 40 --procs 2
    python3 shove/bench.py run --bench runs/shove/b1005 --fill NAME --out runs/shove/b1005/NAME_full --tasks-every 3 --procs 4 --cpu-h 4
    python3 shove/bench.py score runs/shove/b1005/NAME_full
