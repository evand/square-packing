You write Python constructor programs for the square-packing problem: place n non-overlapping unit squares (rotations
allowed) inside the smallest possible square box of side s.  Lower s is better.  For most non-square n the best known
packings beat the trivial grid bound ceil(sqrt(n)) by placing a group of tilted squares cleverly among axis-aligned
ones.

Your program does NOT need to produce an exactly feasible packing.  It produces STARTING LAYOUTS; a strong numerical
optimiser (an augmented-Lagrangian squeeze followed by sequential linear programming on the exact non-overlap
constraints) then shrinks the box and removes overlaps, converging to a nearby local optimum.  That optimiser is good
at local repair and very bad at changing the combinatorial structure (which squares are tilted, their rough angles,
which squares touch which, how rows and columns are arranged and shifted).  The structure is your job.

Interface (the whole file is your program):

    def construct(n, rng):
        # n: number of unit squares (int).  rng: numpy.random.Generator (use it for any randomness).
        # return a list of up to 24 candidates; each candidate is (s, squares) with s the box side (float) and
        # squares a list of n triples (x, y, theta_deg): centre coordinates in [0, s] x [0, s] (lower-left origin)
        # and rotation in degrees.  Small overlaps / slightly too small s are fine: the optimiser repairs them.

Rules: only `math`, `numpy` (as np), `random`, `itertools`, `functools` may be imported; no file, network or process
access; at most a few seconds of CPU per call; it must work for every n from 30 to 330 (return fewer good candidates
rather than crashing).  Use the candidates to vary structure (different parameters, positions, angles), not just noise.

Background (public knowledge): known record packings (Erich Friedman's "Squares in Squares" catalogue and its
successors) mostly combine axis-aligned grid blocks with one or more bands or clusters of tilted squares.  Common tilt
angles are near 45 degrees (Göbel-style diagonal strips), around 40-42 degrees, and families with fixed angles such
as about 24.3 degrees (side k + (sqrt 7 - 1) / 2).  Squares adjacent to a tilted band are often
shifted by fractional amounts so rows / columns interlock with it.  Nobody knows the optimal structure for most n: new
designs are welcome, and a design may be specific to ranges of n.

Each evaluation runs your construct(n, rng) on a set of target n, optimises every candidate, and reports per n the
best side reached, its gap to the best known value, and a description of the best optimised layout.  Fitness is the
mean normalised gap (0 = ties the best known, 1 = no better than the grid bound; lower is better; negative = a new
record).

Reply format: first at most 4 sentences describing the structural idea and what you changed, then exactly one
```python block containing the complete program.
