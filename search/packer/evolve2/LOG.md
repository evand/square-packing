# EV2 prompt lineage and results

Unit of evolution: the prompt (`prompts/P###.md`, on top of `RULES.md`).  Each generation: one Claude Code subagent per
prompt (writes constructors in `runs/evolve2/gen<g>/P###/`, evaluates with `analysis/evolve_eval.py`, reports).
Selection and mutation of prompts: the main session, from the reports.  Targets: live frontier.

| prompt | gen | parent | direction | agent | best normalised gap per target | notes |
|---|---|---|---|---|---|---|
| P001 | 0 | - | distorted strip (178, 207) | Opus | **178 0.000 (tie, reproduced)**, 207 +0.049, 206 +0.124, 177 +0.134 (178 minus one), 208 +1.0 | two staircase upright grids with "events" (jump / stall per row), band centre from the staircase corners, 2-wide rungs; 207 = sheared strip 49° then 41° + 2 notch squares.  Caveat: it measured the reference PNGs by segmentation (centres, angles), so the 178 tie is partly image-traced.  Suggests stating per n the tilted count, events and per-half tilts; drop 208 (different structure) |
| P002 | 0 | - | bent strip (179, 180) | Opus | 179 +0.163 (straightened into a 178-type strip), 180 +1.0, 181 +1.0, 178 +0.152 | lattice-band polyline arms + greedy fill; outcomes a lottery across near-identical variants; loose clean starts (s ~14.2-14.4) and heavily overlapping tight ones both fail; rigid placement too sparse (166-176 of 180).  Arm chirality matters (lattice axis ~10° steeper than the arm).  Suggests building the staircases first and nesting the tilted squares (as P001 did), tilted count per n, starts within ~0.05 of target |
| P003 | 0 | - | two bulbs (132, 182) | Opus | 132 +0.495, 182 +0.121, 176 +0.445; 133 / 183 targets were the grid (my error) | first 10 evals returned the grid: layouts built in a box of side G hid full upright lines (count honesty again); bulbs = sheared strips on a wall, strip of k spans k - δ (k cos t + sin t = k - δ fixes the angle), bulb B = transposed copy of A; measured the PNGs (39 + 93 at 132, 59 + 123 at 182) |
| P004 | 0 | - | block + secondary block (127) | Opus | 127 +0.063 (11.822876 = the 24.3° family, as EV1), 128 +0.062, 129 +0.113, 107 +0.044 | count-driven polyline bands + per-region corner-anchored uprights worked (0.3 -> 0.06); literal rectangular lattice blocks failed (no interlock -> grid); plateau prog_8-12; suggests describing the main block in staircase terms (row shifts, which upright row each edge touches) |
| P005 | 0 | - | staircase diamond (126) | Opus | **126 0.000 (tie, reproduced; the J structure)**, 145 0.000 (tie), 123 +0.061, 125 +0.121 | 45° diamonds with slid rows + exact-ish upright fillers; 126 came from a "seam pocket": the gaps between left/right- and bottom/top-anchored upright grids cross near a corner, and a 2-square diagonal strip there absorbs both.  Own fit screen unreliable (125 winner ranked 16th); hole insertion (build for n-1 / n-2, add squares in holes) useful.  Caveat: the record was shown (gen 0).  Suggests stating the seam-gap mechanism; 125 / 123 are other families |
| P006 | 0 | - | open design (104, 109, 110) | Opus | 104 0.000, 109 0.000 (ties; regular 45° strip / 7x7 diamond), 110 +1.000 | built its own greedy upright fill + bisection proxy (tolerance 0.025 overstated counts by 2); suggests describing 110 explicitly (two mirrored tilted regions of sliding brick rows, ~43 tilted) |

## Gen 0 summary (12:10)
Reproduced (record shown in the refs): 178 (P001, distorted strip), 126 and 145 (P005, staircase diamond / seam
pocket), 104 and 109 (P006, regular).  Not reproduced: 207 (+0.049), 179-181, 127 (stuck on the 24.3° family), 132,
182, 110.  No records.  Pre-registered: (1) >= 1 of 9 reproduce targets tied, 45 %: yes (178, 126); (2) record 7 %: no.
Cross-agent lessons (inputs to gen 1):
1. Count honesty: every agent built its own count proxy and most were wrong (P006 tolerance, P003 box of side G,
   P005 screen ranked its winner 16th) -> gen 1 uses `--fill` (the exact MILP upright fill).
2. Build the upright staircases first, the tilted structure in the channel (P001 worked; P002 failed the other way).
3. Starts must be near the target side with small overlaps (P002, P004); loose clean starts relax to the grid.
4. Mechanisms worth stating: staircase "events" (jump / stall) and band centre from staircase corners (P001); strip
   span k cos t + sin t = k - δ (P003); arm chirality (P002); seam pocket where corner-anchored grids' gaps cross (P005);
   hole insertion (build for n - 1, add squares in holes; P005).
5. Agents segmented the reference PNGs (P001, P003): reproduction partly image-traced -> gen 1 shows only near-record
   packings (`evolve_eval.py near`), never the record, and scores on row segments with a soft-min (Evan 10-15).

## Gen 1 (holdout row 234-240; near-record refs only)
| prompt | parent | direction | best per n (normalised gap) | notes |
|---|---|---|---|---|
| P010 | P004 + P006 | blocks at two angles / composition | **234 0.000, 235 0.000 (hidden records tied blind)**, 236 +0.128, 237-240 1.0 / failed | verified 14:20: feasible (pen 0), 36 tilted each; nearest refs shown were 5e-3 / 7e-3 worse.  234 = 15 + (sqrt 7 - 1) / 2 (the 24.3° catalogue family); 235 = 15.8266056334 is not a family value: the first real holdout reproduction.  The winner is NOT the prompted motif: a 2-wide corner-to-corner band of 18 rungs at 24-28° (19 rungs ~40-45° at 236); the optimiser places the slips.  Every two-angle / composed structure failed (fill 15-20 short).  Method: read its own optimised packings rung by rung to fit rung count and end insets; overlaid candidates on the ref PNGs to check geometry without spending evals.  Found an evaluator bug: softmin drops failed n (a failed n raises the score); fix after gen 1.  Suggests: bands first; blocks only as rungs meeting upright rows at integer heights; per-n tilted counts (~50 at 237 ... ~85 at 240) |
