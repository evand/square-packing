# EV2 prompt lineage and results

Unit of evolution: the prompt (`prompts/P###.md`, on top of `RULES.md`).  Each generation: one Claude Code subagent per
prompt (writes constructors in `runs/evolve2/gen<g>/P###/`, evaluates with `analysis/evolve_eval.py`, reports).
Selection and mutation of prompts: the main session, from the reports.  Targets: live frontier.

| prompt | gen | parent | direction | agent | best normalised gap per target | notes |
|---|---|---|---|---|---|---|
| P001 | 0 | - | distorted strip (178, 207) | Opus | **178 0.000 (tie, reproduced)**, 207 +0.049, 206 +0.124, 177 +0.134 (178 minus one), 208 +1.0 | two staircase upright grids with "events" (jump / stall per row), band centre from the staircase corners, 2-wide rungs; 207 = sheared strip 49° then 41° + 2 notch squares.  Caveat: it measured the reference PNGs by segmentation (centres, angles), so the 178 tie is partly image-traced.  Suggests stating per n the tilted count, events and per-half tilts; drop 208 (different structure) |
| P002 | 0 | - | bent strip (179, 180) | Opus | running | |
| P003 | 0 | - | two bulbs (132, 182) | Opus | running | |
| P004 | 0 | - | block + secondary block (127) | Opus | 127 +0.063 (11.822876 = the 24.3° family, as EV1), 128 +0.062, 129 +0.113, 107 +0.044 | count-driven polyline bands + per-region corner-anchored uprights worked (0.3 -> 0.06); literal rectangular lattice blocks failed (no interlock -> grid); plateau prog_8-12; suggests describing the main block in staircase terms (row shifts, which upright row each edge touches) |
| P005 | 0 | - | staircase diamond (126) | Opus | running | |
| P006 | 0 | - | open design (104, 109, 110) | Opus | 104 0.000, 109 0.000 (ties; regular 45° strip / 7x7 diamond), 110 +1.000 | built its own greedy upright fill + bisection proxy (tolerance 0.025 overstated counts by 2); suggests describing 110 explicitly (two mirrored tilted regions of sliding brick rows, ~43 tilted) |
