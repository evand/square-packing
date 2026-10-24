# Nakajima, soft-to-rigid packing: full paper read (2026-10-24)

Source: github.com/yoheinakajima/soft-to-rigid-packing @ efe7a92 (MIT; clone in `~/math/_untrusted-third-party/nakajima-soft-to-rigid-packing`,
read only, nothing run), paper `paper/main.tex` "Twelve Unit Cubes in a Cube of Side 2.9315, Found by Hardening Balls into
Cubes" (8 Oct 2026, Zenodo 10.5281/zenodo.23248095), `EXPERIMENT_LOG.md`, `src/squares/sim.js`.  Supersedes the README-only
summary in `search/packer/lit-2026-10-12.md` §6.  **[read]** = from the paper / code; **[inference]** = ours.

## Method [read]
- Shape: square of half-side 1/2 - r Minkowski-dilated by radius r; r = 1/2 disk, r = 0 unit square; every member inside the
  unit square, so hardening only grows shapes.  **Same family as fq `--round-r0` and `src/anneal.rs`.**
- Energy: sum over pairs of (2r - sd)_+^2 (sd = signed core distance: Euclidean when separated, minus SAT depth when
  overlapping) + squared wall violations + mu L.  A penalty: small overlaps are tolerated throughout.
- Dynamics: heavy-ball gradient descent (beta 0.85, lr 0.04, per-step velocity cap 0.03), box side L a dynamic variable under
  pressure mu, resized symmetrically.  No constrained solver until the end.
- Schedule (squares code default; paper's cube run is ~10x longer): compress disks at 3 mu (2,500 steps), morph r: 1/2 -> 0
  linearly at mu = 0.02 with the box still free (5,000 steps), rigid settle with mu decaying 1e-4x (1,500 steps).
  "Annealing" = Gaussian kicks to x, y, theta decaying linearly to 0 over compression (0.015) and morph (0.004).
- Finish: legalise by scaling centres about their mean, then SLSQP with corners-in-box and per-pair separating-line
  constraints (our SLP polish does the same job).

## Results [read]
- Cubes: s3(12) <= 2.9315185 (2 of 192 starts below the previous record); n = 13, 14 nothing.
- Squares, from random starts: recovers 11, 17, 18, 19, 26, 27 (hit rates 0.2-10 %); nothing at 28, 29; no square record in
  ~5,000 starts at n >= 17.  Hit rates fall with n at fixed budget; at 26 they rise with hardening time, at 27-29 not.
- Rigid control on the same schedule is worse (n = 11: 6/160 vs 0/120; n = 12 cubes: 13 vs 1 of 32 below the trivial side).
- Octagon intermediate: no help (0/320 vs 4/320).  Angles lock in only once r < 0.06; noise is ~0 by then; "a schedule that
  peaks there is untested" (their limitation v).
- Grid cases: rigid starts win where the grid has empty cells, blobs win for full grids.

## Against ours [inference]
- Shape: identical.  Protocol differs on the point Evan raised (10-24): theirs never holds the box below the hard side; the
  box breathes under constant pressure while r falls, with every overlap a small penalty.  Ours (`round.15c`) compresses at
  r0 = 0.15, then holds the box at that side for 20 hardening relaxations solved to 1e-7, and the final quench expands.
- Start: theirs from disks at random (a generator, n <= 27); ours from structured packings (a re-hardening move).  Our
  10-08..12 runs found deep rounding (r >= 0.2) destroys sub-grid structure, so their r0 = 1/2 is a generator setting, not
  a move setting.  Closest existing analogues here: `src/anneal.rs` (NPT hard-particle MC, disk -> square; 10-08, n = 37
  register hit, all-grid at 110) and the explorer's `melt` (local rounding up to 0.05-0.2 + MC, then re-square).
- Untested ideas from the paper worth a bench: (1) pressure-driven hardening with the box free at every step (vs our hold);
  (2) noise that peaks near the orientation lock-in (r ~ 0.06) rather than decaying to 0 before it; (3) penalty dynamics as
  a cheap low-precision relaxation for early steps (no ALM / SLP until the shape is rigid).
