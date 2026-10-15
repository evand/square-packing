# EV2 rules, generation 1 (shared by every gen-1 prompt; 10-15)

You are one constructor-designer in a prompt-evolution experiment on packing n unit squares (rotations allowed) in the
smallest square box.  Your job: write a Python constructor program that builds a specific STRUCTURE, test it with the
evaluator, look at the renderings, and revise until it reproduces or beats the targets.  Then report.

## The program
A single file `prog_<k>.py` in your work directory defining

    def construct(n, rng):
        # return a list of up to 24 candidates (s, squares): s = box side (float), squares = n triples (x, y, theta_deg),
        # centres in [0, s]^2, lower-left origin.  Small overlaps / a slightly small s are fine.

It is a STARTING LAYOUT generator: the evaluator runs a strong local optimiser (augmented-Lagrangian squeeze + SLP
polish on the exact non-overlap constraints) on every candidate.  The optimiser repairs overlaps and shrinks the box,
but it essentially never changes the combinatorial structure (which squares are tilted, roughly at what angle, which
touch which, how rows / columns are offset).  Getting the structure right is the whole job; precise coordinates are
not needed.  Use the 24 candidates for genuinely different parameter choices (segment lengths, step positions, angles,
which corner, block sizes), not noise.  Typical good layouts: an axis-aligned grid everywhere except the tilted
structure, with the uprights next to it shifted by fractional amounts so rows interlock with the tilted squares.

Sandbox: only `math`, `numpy`, `random`, `itertools`, `functools` imports; no file / process / network access; a few
CPU seconds per call.  Programs are run in a sandbox and rejected by a static check otherwise.

## Tools
Run from `/home/evand/math/square-packing/public/search/packer`:

    python3 analysis/evolve_eval.py eval WORKDIR/prog_<k>.py --n 234,235,236,237,238,239,240 --fill --out WORKDIR/eval_<k> --tag <PROMPT_ID> --procs 2

prints per n: optimised side, target (the best known value today), gap, normalised gap (0 = ties the best known,
1 = no better than the grid, < 0 = new record), a description, and a SCORE line: softmin over the n of the
normalised gaps (one excellent n counts; being bad at some n is fine).  Optimise the softmin.  Writes
WORKDIR/eval_<k>/n<N>.png (look at them with the Read tool: grey = axis-aligned, colour = tilt) and n<N>.txt.

`--fill` (use it): your candidates may contain ONLY the tilted structure (any number of squares <= n, e.g. 40 tilted
squares); the evaluator adds axis-aligned squares by an exact MILP fill (the most uprights that fit around your
squares at your s, on any positions, not just a grid), retrying at s x 1.01 / 1.02 / 1.04 if short, and reports how
many squares the fill reached ("exact fill at the given s reached 233 / 236").  That count is honest: trust it over
any count proxy of your own (generation 0's designers all lost time to wrong home-made counts).  You may still return
full n-square layouts (uprights included); they are then used as given if they have n squares.

References: `runs/evolve2/near/near_n<N>_<i>.png` (N = 234-240, i = 0-3): VALID packings that are WORSE than the
best known (the title says by how much).  The best known packings are deliberately not shown: your score is against
them, so learn the structure from these and beat them.

Constraints: at most 12 eval calls; keep --procs 2 (other experiments share the machine); work only inside WORKDIR
(write programs there; do not edit anything else in the repo).  You may read: this file, your prompt, the parent
program named in your prompt, `runs/evolve2/near/`, and your WORKDIR.  Do NOT read anything else under runs/
(in particular not runs/evolve2/refs/, other designers' work, the store), do not run `evolve_eval.py ref`, and do not
read packing coordinates of any packing you did not construct.  The pictures in runs/evolve2/near/ are fine to study
closely.

## Report (your final message; keep it under 400 words)
1. Best SCORE (softmin) and best result per target n (side, normalised gap), with prog / eval dir.
2. What structure your best program builds, in 3-5 sentences.
3. What worked, what failed, and why you think so (the most useful part for the next generation).
4. A suggested change to the PROMPT (not the program) that would help the next designer.
Any normalised gap < 0 (CANDIDATE RECORD in the evaluator output): say so on the first line.
