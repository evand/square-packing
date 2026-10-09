# Tilings are never optimal: numerical test at k = 2 (2026-10-07/08)

Conjecture (Evan, 2026-10-07; problems.html §4, WISHLIST N17, Lean `Conjectures.TilingsNeverOptimal`): if s(n) is not an
integer, then s(k²n) < k·s(n) for every k ≥ 2.  The large-k case follows from the waste bounds (Lean, conditionally:
`tilingsEventuallyNotOptimal_of`); k = 2 is open.  This note tests k = 2 numerically on the smallest open cases.

**Result.**  For n = 241, 273, 307 (n² − n + 1 family, records by Francisco Couzo 2026 in jlevy's register), every 2×2
start improves, and the best by 0.115–0.127.  **Certified** (exact rational checks, `verify_cert.py` and the independent
`verify_cert2.py`, both VALID):

| 4n | 2 · record(n) | certified s(4n) ≤ | gain | grid |
|---|---|---|---|---|
| 964 | 31.976264879 | **31.859493882** | 0.117 | 32 |
| 1092 | 33.967851925 | **33.840529687** | 0.127 | 34 |
| 1228 | 35.962061097 | **35.847313543** | 0.115 | 36 |

Neither Ellsworth's catalogue nor jlevy's register lists n = 964, 1092, 1228 (checked 10-07: no entries; register case
pages 404), so these are the best packings we know of for those n.  **Evidence at the level of best-known values only**:
it shows s(4n) < 2·record(n), while the conjecture is about the true s(n), which may be below the record.

## Setup
`packer/tile_hop.py` (seeds, quench, hop; `--summary`), `packer/mkcert.py` (f64 packing → rational certificate: centres
spread by 1 + 10⁻⁶, t = tan(θ/2) rounded, side ×(1 + 10⁻⁶) + 10⁻¹²).  Seeds from `exact/batch/inputs/n-N.txt`:

* **T**: plain 2×2 copies.  **M**: copy (i, j) mirrored in x if i = 1, in y if j = 1 (seams meet their mirror images).
* **R**: pinwheel (copy q rotated by 90°·q).  **S**: scale ×2 and split every square into a 2×2 block (tilted blocks).

All seeds checked valid (no overlaps) before running.  Each start: `fq quench` at loosen 1.0 and 1.02, then `hop.search`
(T = 10⁻⁵, patience 40) for 10 minutes from the better one.  8 processes (half the machine; another session was running
the s(110) explorer), worktree build of `fq` at 7f6e9d6 (so a concurrent rebuild could not change the binary mid-run).
Runs: `packer/runs/tile0` (gitignored); table `tilings/summary.txt`; certificates `tilings/n{964,1092,1228}.cert`.

## Observations
* **S (split) gains most**: one quench takes off 0.08–0.10, the 10-minute hop another 0.02–0.04.  The tilted 2×2 blocks can
  shear along their internal cuts, which is the seam heuristic in its cleanest form.
* **T and R gain 0.008–0.038**: the seams relax a little.  **M barely moves** (964: 7 · 10⁻⁷, 1092: 2.7 · 10⁻³, 1228: 4.3 · 10⁻³):
  mirrored seams meet their own reflections, so wall contacts become matched contacts.  The case the heuristic says is hard.
* **The results are not local minima.**  `exactsolve.py` on the 964 T quench: "not jammed", corner-corner MILP descent
  dS ≈ −0.052 (fq's flip search, top 8 candidates, misses it at n ≈ 1000).  The first exactsolve attempt (default tol)
  crashed in `max_support` (linprog: empty contact set at tol 1e-9, presumably; f64 output converged only to ~1e-8); with
  `--tol 1e-7` it runs.  So the certified values are far from optimised: one quench costs 2–7 minutes at n ≈ 1000, and
  each hop got only 1–3 proposals.
* **Control n = 4 · 65 = 260** (record 16.657, below the doubled 17.071 and the grid 17): quenches do not move the tilings;
  hopping (25–150 proposals, 10 min) reaches 16.901 (S) / 17.03–17.04 (T, M, R).  At this budget the search rediscovers
  nothing near the record, as expected from the explorer work at 110.

## Split rigidity: when is the S start itself jammed? (10-08, later)
Question (Evan): can "splitting relieves shear" become a proof?  Template: if s(4n) = 2s(n), the split S(P) of every
optimal P is optimal, hence jammed in every corner-corner branch.  So **one optimal P whose split is not jammed gives
s(4n) < 2s(n)**; T, M, R and mixtures give further necessary conditions (M is generically jammed: seams match wall to
mirrored wall).  Jammed split ⇔ P has a stress that survives splitting: (a) no shear resultant across either midline
of a loaded square, no tension in a cut (frictionless dry joints); (b) no loaded contact of a differently oriented
square at an edge midpoint (it becomes a vertex-vertex hinge).  Lemma (easy): all loaded squares axis-parallel ⇒ s ∈ ℤ
(a horizontal force path crosses flush unit widths).

Test: split every non-integer register record n ≤ 100 (`tile_hop.py` variant S), `fq quench`, then `exactsolve` on
the splits that did not move (scratch, not kept; ~1 CPU-h):
* **s(5), s(10): not jammed** (s(5) split: dS = −0.47 per unit motion; quench goes to 5.0, all axis-parallel).  Mechanism
  (b): with the 4 centre quarters' rotations locked the MILP finds no descent; the quarters pinwheel-rotate, their
  vertices slide off the corner points.  The s(5) stress itself has no midline shear (a) — the naive check passes.
* **Jammed (local minimum modulo exact flat motions, jammed in every branch):** splits of 27, 40, 52, 65, 67, 82, 89
  (all 0°/45°), 84 (0°/45°, first pass only: re-solve after tiny overlaps timed out) and **86** (18 squares at 24.3°,
  S = (17+√7)/2, first pass only).  86 is a 2-wide staggered band, steps (1, ½) in the square frame, flush the whole
  length: split, it is a strip of the tilted unit lattice (steps (2, 1)), no new hinges.
* Not jammed: 11, 17–19, 26, 28, 37–39, 50, 51, 53, 55, 68, 83, 87 (quench drop), 29, 41, 54, 69–71, 88 (MILP descent).
  Inconclusive: 66 (KKT diverged), 85.
* **None of the jammed splits is competitive**: 2s(n) > ⌈2√n⌉ for all of them.
* 126 ph14 staircase band (`trio126/`): split not jammed, dS = −0.22; still −0.084 with all rotations locked
  (pure slip = mechanism (a)); locking all tilted quarters, or all axis quarters, kills it: the slip moves band and grid
  together.  `fq quench` on that split goes *up* to 24 (it cannot switch corner-corner branches).

Consequences.  "Non-grid P ⇒ split not jammed" is false, so a local proof must use competitiveness.  k = 2 is
nontrivial only for n ∈ [k² − k + 1, k² − 3] with s(n) < k (register n ≤ 324: 211, 241, 273, 307 only; 211 untested).
Those records carry many generic angles; flush same-orientation bands (86) are the structure such a proof must exclude.
Also: if the k-tiling is optimal so is every divisor tiling, so prime k suffice.

## Status
Paused (Evan, 10-08): the numerical evidence is in, and it fits the seam picture; optimising these packings waits for
the optimizer work in `packer/`.  Open: k = 3 (2169, 2457, 2763 squares) is out of reach of the current quench; a seam
lemma (a 2×2 tiling of a wall-touching packing is never jammed) would be the theoretical version.
