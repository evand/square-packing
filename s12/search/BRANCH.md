# Branching on corner occupancy: certificate + case split (B3, first level)

Code: `search/branch.py` (LP with one Lagrange multiplier), `verify/` (region trailer),
`xcheck.py` (exact re-check), `lean/Sqpack/Basic.lean` (`packing_le_weight_region`),
`tests/rejection_tests.sh` (nine trailer tests), `certificates/FORMAT.md` ("Branch certificates").
Status 2026-08-28: **method built and verified end to end (single and per-box multipliers); it
does not close `s = 3.99`, and at `s = 3.98` the all-corners leaf is exactly 12** — see the tables.  Everything numeric below is a cutting-plane LP
value that can only *rise* as the loop continues (rows are added; column generation is priced
out to 1–3 %), so an unconverged value above 12 is final and one below 12 is not.

## The reduction

Pure method: cover `w` with `w(S) >= 1` for every closed unit square `S` in `[0,s]^2`, total `W`;
`n` interior-disjoint squares give `n <= sum_i w(S_i) <= W`.  It is dead for `s >= 3.99`
(`DUAL_EXACT.md`: a fractional packing of mass 12.008 exists), and its dual optimum spreads mass
0.85 over each corner — something no integral packing does.

Branch.  Let `R` = the four closed corner boxes `[0,r]^2` and images, taken by the **centre** of a
pose.  For a leaf `k` (the number of squares of the packing centred in `R`), a cover `w` and a
real `λ` with

    w(S) >= 1 + λ   (centre in R),        w(S) >= 1   (otherwise)

give  `n + λ k <= sum_i w(S_i) <= W`, so `W − λ k < n` refutes every packing with exactly `k`
squares in `R`.  Two interior-disjoint unit squares have centres `>= 1` apart, and the admissible
part `[1/2, r]^2` of a box has diameter `(r − 1/2)√2`, so for `(2r−1)^2 < 2` each box holds at
most one centre, `k ∈ {0,…,4}`, and five leaves cover everything.  `λ` is a free LP variable:
`min W − kλ` s.t. `A x − flag_r λ >= 1`.  By LP duality the leaf value is exactly the maximum
mass of a fractional packing with total corner mass `k` — branch-and-bound with the cover LP
as node bound.  (Lean: `packing_le_weight_thresh` for a general per-pose threshold,
`packing_le_weight_region` for `1 + λ·[R]`; 0 sorries, standard axioms.)

Leaf `k = 0` drops the corner rows (`λ = −2`); `k = 4` wants `λ > 0`; `1 ≤ k ≤ 3` is free.
Because one `λ` is shared by the four boxes the leaves stay D4-symmetric and the verifier's
`[0°, 45°]` reduction survives; the price is that a mixed leaf cannot reward its occupied
corners and free its empty ones at once (see "What was learned").

**Per-box multipliers** (`--k 1100`, trailer `lambda L1 L2 L3 L4 / k K1 K2 K3 K4`): one `λ_j` per
box and an occupancy pattern `K ∈ {0,1}^4`; `n + Σ_j λ_j K_j ≤ W`, so the leaves are the sixteen
patterns, six up to symmetry (`0000, 1000, 1100, 1001, 1110, 1111`; `0000` and `1111` coincide
with the single-`λ` leaves `k = 0, 4`).  Unequal `λ_j` break the symmetry of the claim, so the
verifier sweeps `[0°, 90°)`, every point is its own LP column, and pricing runs over the whole
container.  (Lean: `packing_le_weight_regions`.)

## Verification

`verify/` reads the trailer `region corner r_num r_den / lambda L / k K`, checks `(2r−1)^2 < 2`
exactly, and in the arrangement sweep requires `1 + λ` of every cell that may meet a box, `1` of
every cell that may leave the boxes, both of a straddling cell — the membership tests are on
the cell's bounding box (exact for convex cells against axis-parallel boxes), in `f64` with a
padding that can only make them stricter.  Witnesses carry a flag (which threshold failed) and
are chosen inside the cell clipped to the admissible polygon on the required side of the box
boundary: the first version used the cell midpoint, which for the large empty cells of the
`k = 0` leaf lands outside the admissible box or inside the region and fed the LP constraints
the leaf does not need (the `k = 0` LP came out at 14.0 instead of 11.6 on the 224-point
certificate — a heuristic defect, not a soundness one, but worth recording).  With `λ = 0` the
output is bit-identical to the plain verifier.  `xcheck.py` classifies the same cells in exact
rational arithmetic with no padding.  Both agree on every test in `tests/rejection_tests.sh`
(57 checks, 9 for the trailer).

## Results (container side `t`, corner box `r`, LP value `W − kλ`; `< 12` closes the leaf)

Starting point in every run: the shipped 1736-point certificate rescaled to `t` (`--mul/--Dp`),
weights re-optimised with the verifier as separation oracle at `N = 6000`, 40 rounds of column
generation.  Iteration counts in brackets; `*` = still running when tabulated.

`t = 3.99`

| leaf | r = 0.8 | r = 1.0 | r = 1.2 |
|---|---|---|---|
| k = 0 | 11.83 [12*] | 11.35 [18*] | 10.71 [31*] |
| k = 1 | 12.00 [12*] | 11.65 [16*] | 11.40 [39*] |
| k = 2 | 12.07 [9] | 11.92 [10] | 11.89 [32*] |
| k = 3 | 12.17 [9] | 12.12 [9] | 12.17 [19] |
| k = 4 | 12.00 [37] | 12.00 [9] | 12.10 [17] |
| pure | | **12.20 [45*]** | |

`t = 3.98` and `3.975`, `r = 1`

| leaf | t = 3.975 | t = 3.98 |
|---|---|---|
| k = 0 | 11.17 [12*] | 11.20 [14*] |
| k = 1 | 11.45 [10*] | 11.51 [11*] |
| k = 2 | 11.71 [9*] | 11.79 [10*] |
| k = 3 | 11.90 [10*] | 11.99 [12*] |
| k = 4 | 11.94 [12*] | 11.97 [12*] |
| pure | 11.96 [25*] | **12.02 [26*]** |

The `k = 4` runs at `r ≤ 1` repeatedly settle at `W − 4λ = 12·(1 + margin)` exactly (e.g.
`W = 15.05, λ = 0.76`): the dual has found a fractional packing of mass 12 with mass exactly 1
in each corner box, so that leaf cannot close at 3.99.

`t = 3.98`, `r = 1`, restarted loops (IPM, `N = 2000`, column checkpoints; 2026-08-28)

| leaf | value | probe min | status |
|---|---|---|---|
| k = 0 | 11.243 | 0.996 | converging |
| k = 1 | 11.557 | 0.996 | converging |
| k = 2 | 11.832 | 0.996 | converging |
| k = 3 (shared λ) | **12.000** | – | open: crossed 12 at round 14 |
| k = 4 | **12.000024** | – | open: degenerate `12·(1+margin)` structure, 10 orbits |
| 1110 (per-box λ) | 11.93 [12*] | 0.87 | rising ~0.015/round |

The `k = 4` endpoint is worth recording.  The LP settles on an 80-point cover with dyadic weights
`1/2, 3/8, 1/4, 1/8` on eight orbits of points hugging the grid lines `x ≈ 1, y ≈ 1, 2` —
`(0.94, 0.99), (0.97, 1.00), (0.84, 0.99), (0.91, 1.99), (1.61, 1.77), (1.15, 1.52), …` — total
weight 18 with `λ = 1.5`, i.e. cost exactly `18 − 4·1.5 = 12`.  This is a Nagamochi/Friedman-style
grid cover; its dual is a fractional packing of mass exactly 12 with one unit of mass in each
corner box.  So in the all-corners case the corner relaxation is already tight: the packing puts a
whole square in each corner, and the entire LP gap sits with the *eight non-corner squares* in
the wall strips and the interior.  Per-box multipliers cannot help this leaf (it is symmetric),
and no `r` fixes it at 3.98: smaller boxes make `k = 4` easier and `k = 0` harder in step.

## Shipped leaf certificates (`certificates/branch/`, `verify_branch.sh`)

| file | leaf | points | `W` | `λ` | `W − λk` | Rust `N=2000 / 4000` | exact Python (`xcheck.py`, N=2000) |
|---|---|---|---|---|---|---|---|
| `s12_t3.98_corner_k0.txt` | k = 0 | 9984 | 11.2447760 | −2 (corner free) | 11.2447760 | 1.0000206 / 1.0000206 | 1.0000206 (8 h 15 min, agrees) |
| `s12_t3.98_corner_k1.txt` | k = 1 | 11892 | 11.2721100 | −0.2879786 | 11.5600886 | 1.0000107 / 1.0000107 | queued |
| `s12_t3.98_corner_k2.txt` | k = 2 | 11168 | 11.3399992 | −0.2475135 | 11.8350262 | 1.0000009 / 1.0000009 | queued |

Each refutes every packing of 12 unit squares with exactly `k` squares centred in the corner
boxes `[0,1]^2` (and images) of a container of side `< 3.98`.  They are not a bound: `k = 4` is
open.  They are the worked examples of the format at full scale.  The `1110` per-box run
(11.966, probe 0.94 at round 37) was still finalizing when this was written.

## What was learned

* **The pure cover at 3.99 costs ≈ 12.20**, not ≈ 12.01: `L(3.99) = 12.008` (`DUAL_EXACT.md`) is
  a lower bound on `COVER(3.99)` from a non-converged dual, and the primal is 0.2 above 12.
  Pure `COVER(3.98) ≈ 12.02`, `COVER(3.975) ≈ 11.96` (unconverged, rising): the pure ceiling
  `s*` is a little below 3.98, consistent with `TIGHTEN.md` (column generation at 3.9696 already
  failed there, but from a poorer starting point).  So the excess grows by roughly `0.1` per
  `0.01` of `t` near the ceiling, and a case analysis at `s = 4` (closed semantics) faces the
  `12.3–12.5` of `CLOSED4.md`, not `12.008`.
* **One level of corner branching recovers 0.03–0.1 on the worst leaf.**  The order of the
  leaves is what the dual predicts (emptier corners are cheaper: `k = 0` gains 0.4–1.5), but the
  mixed leaves are the hard ones, and `k = 3` is the worst for every `r`.
* **Box size is a genuine trade-off**: a large box lets the fractional "corner square" smear over
  many poses inside it (`k = 4` at `r = 1.2` gains only 0.1 on pure), a small box pins it but makes
  `k = 0` expensive.  No `r` closes all five leaves at 3.99.
* **The corner level is exhausted at about the pure ceiling.**  At 3.98 the leaves `k = 0, 1, 2`
  close with room (11.24 / 11.56 / 11.83) but `k = 4` is exactly 12 and `k = 3` (per-box `1110`)
  is at 11.93 and rising; at 3.99 nothing closes.  The gains where they exist are large (0.8 for
  `k = 0`), but the leaves an integral packing actually resembles — corners full — gain nothing,
  because there the fractional packing is already integral in the corners.  A new bound from this
  level alone would sit within ~0.005 of the pure ceiling (`3.975–3.978`).  The next level is
  wall-strip occupancy: the `k = 4` cover above is the natural seed, and its dual measure says where
  the eight units of mass go.
* **The shared `λ` is the visible weakness of mixed leaves.**  In the `k = 3` leaf the LP settles at `λ ≈ −0.13`:
  it gives the occupied corners' surplus away in order to exempt the empty one.  Per-corner
  multipliers `λ_1..λ_4` (leaves `1000, 1100, 1010, 1110` up to symmetry; verifier over
  `[0°, 90°)`, non-symmetric point sets) strictly improve every mixed leaf and are the next
  step at this level.  The level after that is wall-strip occupancy (eight boxes), Bentz-style.
* **Compute** (measured 2026-08-28).  A round = LP + verifier + pricing.  LP: HiGHS IPM is the best
  of the scipy methods (16 s vs 23 s dual simplex at 31k rows × 388 columns) but its cost is
  steeply superlinear — 15–20 min at 100k rows × 1000 orbit columns, hours for the dual simplex
  at 120k rows; a highspy warm-basis backend as first written was *slower* (200–600 s where scipy
  took 0.4 s) and is left disabled.  Verifier: quadratic in the atoms — 54 CPU-min per call at
  `N = 6000` on a 3400-atom probe (12 min on 5 threads), 2.4× less at `N = 2000`, which costs
  nothing against leaf margins of `10^-2`.  Loop hygiene that mattered: keep near-tight rows when
  pruning (pruning by dual alone cycles), never drop the rows added in the same round, checkpoint
  the full column set (a restart from the probe file loses the zero-weight columns and a 0.01
  pricing grid does not recover finely placed points), warm-start a restart with the whole
  lattice (`--warm-thr 3`), keep pricing on for asymmetric leaves (their objective froze the
  round pricing stopped), and merge near-duplicate rows (`--row-grid`).  Asymmetric leaves are
  ~8× wider (single points, not orbits) and need far more rows to pin down.  Even so a leaf is
  3–8 h; a second branching level with dozens of leaves needs a different LP strategy (the dual
  packing LP directly, or a column-and-row generation with a basis that survives).

## Reproduce

```sh
cd verify && cargo build --release && cd ..; mkdir -p runs
# leaf k at container 3.99 (= 15680*57 / 224000), corner box r = 1
python3 search/branch.py certificates/s12_lower_3.9686.txt L99k4 --k 4 --r 1 --Dp 224000 --mul 57 --colgen 40
# 3.98: --mul 199 --Dp 784000;   3.975: --mul 159 --Dp 627200;   pure control: search/tighten.py reopt with the same --Dp/--mul
./verify/target/release/verify runs/branch_L99k4.txt 12 6000 32 0      # "VERIFIED: (branch k=4) ..." or NOT VERIFIED
python3 xcheck.py runs/branch_L99k4.txt 6000 --all --n 12                # independent exact re-check
```
