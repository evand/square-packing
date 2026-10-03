# SOS probe: can exact algebraic certificates pass the zero-margin wall? (rehearsal on s(6) = 3)

Why: every s(12) route closed because interval/box methods lose at second order against a deficit that is zero
on a plateau (`notes/n12-gap.md` §4.7, §4.11 obstruction X).  Exact Positivstellensatz/SOS certificates on a fixed
separation type do not go through intervals.  Rule (n12-gap §4.11): no s(12) route that cannot do s(6) = 3.
This is a bounded probe: measure, don't build the full pipeline.

## Setting
Square i: centre (x_i, y_i), u_i = tan(θ_i/2), θ_i ∈ [0, π/2) (cos, sin rational in u_i).  Containment in [0,T]²
and "pair (i, j) separated by edge direction d of square i or j" are polynomial inequalities after clearing
denominators.  A **separation type** fixes one separating direction (and side) per interacting pair; far pairs
need no constraint.  Question per system: is `T ≥ 3` on its feasible set (equivalently: the system with T < 3 is
empty), and can SOS prove it with an exact rational certificate?

The binding configurations at n = 6, T = 3 (read `search/S6_SKELETON.md` §1, §4.1–4.2, `search/S6_LOCAL.md`):
1. **Row core**: grid-minus-3 plateau held shut by a wall-to-wall row of three (S6_SKELETON §1.2); margin exactly 0
   on a plateau.
2. **Pinwheel core** near axis-parallel, the hard direction (deficit `−(1/4)t³`, S6_SKELETON §4.2, S6_LOCAL §2–5).
Start with the core sub-systems (fewest squares that hold the point shut), then add context squares.

## Measurements (the deliverable)
1. **Degree**: the smallest SOS/Putinar (or Stengle infeasibility) degree that proves each system, or the largest
   tried without success.
2. **Exactness**: does the float SDP solution round to an exact rational certificate (check the identity in exact
   arithmetic, e.g. sympy/fractions)?  At a degenerate minimum the SDP is typically not strictly feasible;
   report whether facial reduction or the infeasibility form (`T < 3` strict, Stengle) is needed and works.
3. **Locality**: does a certificate for the core alone suffice, or does it need the other squares?  Does one core
   certificate transfer across nearby separation types (same core, different context)?
Also: SDP sizes, solve times, and an honest extrapolation to n = 12 (leaf A / rank-8: 8–12 squares; use the
sparsity: constraints only couple neighbouring pairs).

## Steps
0. Literature: SOS / Positivstellensatz / Lasserre applied to square packing or similar rigid-body packing
   lower bounds (who, what degree, exact or float).  Short.
1. Encode the two core systems; sanity-check by sampling feasible points at T = 3 (the plateau) and checking
   infeasibility numerically at T = 3 − ε.
2. SOS at increasing degree; record (1)–(3).  Stop at the first conclusive outcome either way per system.

## Rules
- Software: create a venv under `runs/sos_venv/` (ignored) and pip-install only well-known open-source packages
  (cvxpy, clarabel, scs, sympy, numpy, scipy, python-flint).  No system installs, no commercial solvers.
- CPU: no pinning; keep total busy processes ≲ 16; budget ~2·10⁴ CPU-s, ask before 3× that.
- Write-up `search/SOS_PROBE.md` ([proved]/[measured]/[heuristic]); code `search/sos_probe*.py`; runs `runs/sos_*`.
  Commit code + write-up to public main with only your own files (no `git add -A`); don't push.
