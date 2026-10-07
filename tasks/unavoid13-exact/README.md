# Brief: exact certificates for the two unavoidable-set integrality results (2026-09-22)

Working tree: `~/math/square-packing/public/s12` (the public repo's `s12/`; commit there, **never push**).
`runs/` is gitignored; `tasks/`, `TODO.md` are private symlinks — don't commit them.  Budget: a working day; stop and
report if the approach below hits a wall.  Launch anything over 10 minutes detached (`setsid nohup … &`, log to a file,
kill by PID, never `pkill`); leave no `until … sleep` / `tail -f` watchers behind.  Use at most 16 threads.

## The two claims

1. **T = 3.**  The minimum pure unavoidable set for `[0,3]²` (a finite point set that every closed unit square inside
   `[0,3]²` contains) has exactly 7 points.  `≤ 7`: `certificates/unavoid13/ks7_rational_3.txt`, certified by
   `search/unavoid13_check.py` (SEG primitive) — **out of scope, already exact.**  `≥ 7` is the claim here: the 575-square
   family `runs/unavoid13_t3a/final_family.txt`, shrunk to 87 squares (`runs/unavoid13_t3a/recheck_shrunk_family.txt`),
   admits no 6-point hitting set.
2. **C2.**  No 13-point unavoidable set for `[0,4]²` is invariant under the half-turn `p → (4,4) − p`.  Family
   `runs/unavoid13_symC2/F_round01.txt` (221 square orbits; orbit-candidate instance in `runs/unavoid13_symC2/*symcheck*.json`).
   D4, C4, V and the diagonal Klein group all contain the half-turn, so C2 alone carries all of them.

Read `notes/unavoid13.md` §1 (soundness: vertex-maximality, why arrangement vertices dominate every point), §2 (T = 3),
§4 (symmetric case: orbit candidates, points on the centre, why `F` must be G-closed) and `notes/publish-inventory.md`
footnote `*` (the adversarial check that motivates this task).

## What is wrong today

Geometry is exact (Fraction vertices, exact incidence, domination).  **Infeasibility is not**: it rests on `highspy`,
`scipy.milp` and a B&B with float LP bounds.  Reliable, not a certificate.

## What to build

1. **Promote the families** into `certificates/unavoid13/` (rational poses, the exact format they already use, plus a
   short `certificates/unavoid13/README.md` saying what each file certifies and how to check it).  For T = 3 ship the
   87-square shrunk family (confirm it alone gives `h ≥ 7` exactly; keep the 575 in `runs/`).
2. **An exact infeasibility certificate** for each instance: a branch-and-bound tree over the hitting-set IP
   (`min Σ c_j x_j` s.t. `A x ≥ 1`, `x ∈ {0,1}`; for C2 the variables are orbit candidates with cost = orbit size,
   rows = square orbits) whose every leaf is closed by either (a) an uncoverable row under the fixings, or (b) a
   **rational dual-feasible solution of the leaf LP** whose objective exceeds `k` (T = 3: `k = 6`; C2: `k = 13`, costs
   are 1 or 2 so the bound must exceed 13).  Generate duals in floats, then rationalise and repair to exact dual
   feasibility (safe bounding: take `y ≥ 0` on the covering rows, and cover any reduced-cost violation with the
   upper-bound duals — the bound stays valid for any `y ≥ 0`).  Use costs' integrality to prune at `> k`, i.e. `≥ k+1`
   after rounding up only where the argument allows it — state exactly which rounding you use.
3. **An independent checker**, `search/unavoid13_exactcheck.py`, in `fractions.Fraction`, that shares **no code** with
   the `unavoid13*` modules: from the family file alone it (i) enumerates the arrangement vertices exactly (square
   corners, pairwise edge intersections, and for C2 the orbit-specific points `notes/unavoid13.md` §4 requires),
   (ii) computes exact incidence, (iii) builds the candidate rows and checks every vertex's row is dominated by a
   candidate's, (iv) reads the tree certificate, checks that the branches partition the space and that every leaf
   certificate is valid exactly.  Output `VERIFIED` / `NOT VERIFIED`, like the other checkers.  Rejection tests
   (a mutated family that should become feasible, a tampered dual, a dropped subtree) in `tests/unavoid13/`.
4. Wire both checks into `verify.sh` (fast tier if they run in minutes at 4 threads, else `--full`), and update
   `notes/unavoid13.md` (the two **[proved]** statements, now with the certificate and checker) and the
   `s12/README.md` bullet that currently says "the infeasibility steps currently rest on floating-point MIP solves".

## Report

A short section appended to `notes/unavoid13.md` ("Exact certificates, 2026-09-2x"): sizes (nodes, leaves, candidates),
checker runtime, what exactly is verified and what is still trusted, anything that did not hold.  If either instance
cannot be certified within budget, say what the obstacle is and label the claim honestly ("exact instance;
infeasibility by two independent MIP solvers") in the README instead.  Commit in `public/` (not pushed) with a clear
message.  Final reply: ≤ 200 words.
