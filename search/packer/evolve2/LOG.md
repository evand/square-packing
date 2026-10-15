# EV2 prompt lineage and results

Unit of evolution: the prompt (`prompts/P###.md`, on top of `RULES.md`).  Each generation: one Claude Code subagent per
prompt (writes constructors in `runs/evolve2/gen<g>/P###/`, evaluates with `analysis/evolve_eval.py`, reports).
Selection and mutation of prompts: the main session, from the reports.  Targets: live frontier.

| prompt | gen | parent | direction | agent | best normalised gap per target | notes |
|---|---|---|---|---|---|---|
| P001 | 0 | - | distorted strip (178, 207) | Opus | running | |
| P002 | 0 | - | bent strip (179, 180) | Opus | running | |
| P003 | 0 | - | two bulbs (132, 182) | Opus | running | |
| P004 | 0 | - | block + secondary block (127) | Opus | running | |
| P005 | 0 | - | staircase diamond (126) | Opus | running | |
| P006 | 0 | - | open design (104, 109, 110) | Opus | 104 0.000, 109 0.000 (ties; regular 45° strip / 7x7 diamond), 110 +1.000 | built its own greedy upright fill + bisection proxy (tolerance 0.025 overstated counts by 2); suggests describing 110 explicitly (two mirrored tilted regions of sliding brick rows, ~43 tilted) |
