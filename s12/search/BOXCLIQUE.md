# Box cliques: verifiable clique columns (task B, 2026-08-29)

Code: `verify/src/main.rs` (parser, exact core check, sweep credit), `xcheck.py` (the same in exact
rationals, independently), `tests/rejection_tests.sh` (§18–19, 25 checks), `lean/Sqpack/Basic.lean`
(`packing_le_weight_cliques`, `clique_of_cores`, `card_filter_clique_le_one`),
`search/boxclique.py` (LP driver), `certificates/FORMAT.md` ("Clique certificates").

## What was built

`search/CLIQUE.md` showed that the extremal fractional packings of the cover LP put mass 1.3–1.5 on
*cliques* of the pose-overlap graph that have no common point (rotated squares are not Helly),
and that clique cuts move the packing LP by 0.1–0.3 near `t = 4`.  For a certificate a clique
must be an object the verifier can check.  The family chosen here:

* A **box** is `(k, [U0LO,U0HI] × [U1LO,U1HI])`: all poses with angle in bin `k = [θ_k, θ_{k+1}]` of
  the verifier's net and rotated centre `R(−θ_k)·c` in the closed rectangle (integers over `Q`).
* A **box clique** is a finite union of boxes.  Its **core** per box is the axis-parallel rectangle
  `[U0HI − h_k, U0LO + h_k] × [U1HI − h_k, U1LO + h_k]` in the bin's frame, `h_k = σ_k/2` with the
  verifier's rounded-down `σ_k`: every unit square with a pose in the box contains it (the same
  shrink lemma the whole verifier rests on).  **The verifier refuses a clique unless every pair of
  its boxes (each with itself) has intersecting cores** — an exact separating-axis test between two
  rotated rectangles in `i128` (`cores_meet`); an empty core (rectangle wider than `2h_k`) is refused
  too.  Nothing else about a clique is trusted.  No new geometric lemma is needed, which is why this
  family was preferred to `K(p, A)`.
* In the sweep of bin `k` a cell is credited `w_K` iff it lies inside one of `K`'s boxes of that bin
  (exact integer comparison of the four sides over the common denominators); a cell that meets a box
  without lying inside gets nothing, and if it is violated its LP witness is placed *outside* the
  box (a corner of the clipped cell), so the cutting-plane loop cannot cycle on a straddling cell.
* Cliques have no representable images under the container's symmetries (the image of a bin is not
  a bin), so a certificate with cliques is always swept over `[0°, 90°]`, and the `N` in the block
  must be the `N` the verifier is run with.

Format (after the `m` point lines, before an optional `region` trailer):

```
cliques N Q c
w_1 b_1
k U0LO U0HI U1LO U1HI      (b_1 lines)
...
```

Total weight = points + cliques; the claim is `captured(S) = Σ_{p∈S} w_p + Σ_{K∋S} w_K ≥ 1` for every
closed unit square `S` inside the container, and `total < n`.

## Verification of the verifier

* **Bit-identity without cliques.**  The new binary was diffed against the pre-change binary on:
  the main 1736-point certificate at `N = 2000` (single-threaded, and in witness mode with
  `topk = 6`: stdout and the 3652-line witness file identical), the 224-point and 56-point
  certificates, and a per-box branch trailer mutant (`lambda -5 -5 0 0`) in witness mode.  All
  identical.  (Multi-threaded runs differ only in which of the first eight `FAIL` lines the racing
  threads print — that was already nondeterministic.)
* **Rejection tests** (`tests/rejection_tests.sh`, now 87 checks, was 62): a harmless clique
  verifies for `n = 13` and is reported as such with the full sweep; its weight counts in the total
  (`4/5` more → `WEIGHT NOT`); cores that do not meet, an empty core, the wrong `N`, a negative
  weight, a clique with no boxes, `k = N`, `LO > HI`, a short block, and a block *after* the region
  trailer are all `ERROR`s; a block *before* the trailer verifies as a branch certificate; a clique
  of weight 1 placed away from the failure of the one-point-deleted mutant leaves it rejected with
  every witness flagged 0 (no leak of clique weight).
* **`xcheck.py`** parses the block, checks the pairwise cores in rationals with the *exact*
  `σ_k` (its cores are supersets of the Rust ones, so anything the Rust accepts it accepts), credits
  cells exactly, and agrees with the Rust verifier on every smoke file and on the demonstration
  certificate below.
* **Lean** (`lean/Sqpack/Basic.lean`, 0 sorries, axioms `propext, Classical.choice, Quot.sound`):
  `packing_le_weight_cliques` — points `A, w ≥ 0` and cliques `K_j, v_j ≥ 0` with every closed unit
  square inside `C` covered by `≥ 1` ⟹ `n ≤ Σ w + Σ v` for any packing of `n` squares of side
  `L > 1`; `card_filter_clique_le_one` (a packing has at most one square in a clique);
  `clique_of_cores` (pieces with pairwise-meeting cores contained in every square of the piece form
  a clique — exactly the verifier's check, with the shrink lemma as the hypothesis `hcore`).

## The driver (`search/boxclique.py`)

`tighten.py`'s loop (D4-orbit point columns, rows from the verifier's witnesses, exact verifier as
separation oracle) plus **clique-orbit columns**: the eight images of a clique share one weight
(cost 8, coincident images merged), so the LP stays symmetric.  Pricing on the D4-averaged dual
`ȳ` — a genuine fractional packing at every atom, so `ȳ(K) > 1` is a real non-Helly gain over the
current atom set and the column's reduced cost is `8 (1 − ȳ(K))`.  The dual comes from an
interior-point solve without crossover (`highspy`): a vertex dual of this very degenerate LP is a
useless pricing signal — the first version priced cliques of "mass 5.9" on the simplex dual that
the LP then never used.  Candidates: max-mass clique of the closed-overlap graph of the support
images (squares shrunk by `eps`), then the best clique through each of the heaviest images; every
image clique is boxed (`± eps` around each pose, outward-rounded over `Q`, grown to `± 3 eps` if
the cores still meet) and pruned to a family that passes the exact pairwise test with the
verifier's own `σ_k`, so every clique written is one the verifier accepts.

## Demonstration, and what the LP said

All runs: the 224-point certificate `s12_lower_3.931795_sparse.txt` (29 orbits), `N = 2000`, full
`[0, 90)` sweep, 4 threads; a round is seconds.

| run | container | columns | LP total | clique weight | verdict |
|---|---|---|---|---|---|
| `ctrl224` | 3920/997 = 3.93180 | points | 11.985894 | – | VERIFIED |
| `cl224b` | 3920/997 | points + 43 clique orbits priced over 8 rounds (ȳ(K) up to 1.86 on 48–234 poses) | 11.985894 | 0 (none used) | – |
| `ctrl224_1992` | 980/249 = 3.93574 | points | 12.015331 | – | not a certificate (≥ 12) |
| `cl224_1992` | 980/249 | points + 17 clique orbits (ȳ(K) up to 2.13) | 12.015331 | 0 (none used) | not a certificate |
| `force224` | 3920/997 | points + first clique orbit held at 0.01 | 12.065894 | 0.080 | ≥ 12 |
| **`demo224`** | 3920/997 | points + first clique orbit held at 0.001 | **11.993903** | 0.008 | **VERIFIED** (Rust at N = 2000; `xcheck.py --all` agrees) |

**Shipped:** `certificates/s12_boxclique_demo_3.9318_N2000.txt` — 224 points and 8 cliques (the
D4 images of one 48-box clique on 33 angle bins, weight 10000/10^7 each), total 11.9939028 < 12,
`s(12) ≥ 3920/997`.  It is in `verify.sh`, `SHA256SUMS` and the rejection tests (§20: verifies;
one box shifted → `ERROR`; clique weights zeroed → the points still verify).  It is a
demonstration of the *format*: the bound is the one the points give alone, and the clique is not
load-bearing — the LP put the point weights at exactly their point-only optimum and the clique
weight on top.

**The negative result.**  On a *fixed* point set the LP never uses a box clique, at the set's own
container (11.9859) or 0.004 above it where the points alone give 12.015.  Reduced-cost pricing
is not the reason: the priced cliques have ȳ(K) = 1.7–2.1 under the dual returned (simplex or
interior), but with 29 orbit columns and 20–40k rows the optimal dual face is enormous and the
LP simply moves to another optimal dual with mass ≤ 1 on every offered clique — adding the
columns changes nothing.  Two readings, not distinguished here:

1. the non-Helly mass of `search/CLIQUE.md` (1.3–1.5 on the *converged* measures, all points
   available) is not present in the fixed-224-point sub-relaxation, whose extremal duals are far
   from unique and can dodge any finite family of cliques; or
2. box cliques as generated (small boxes around support poses, `eps = 0.005`) are too thin to
   catch it — a box clique can never represent a point clique exactly (the core of "all poses
   through `p`" is `{p}`, and boxes are rounded inward), so what they add is interior mass only.

The way to settle it is the one already planned: clique columns *inside* a point column-generation
loop (`tighten.py --colgen`, where the atom set is rich and the dual is pinned) or, cheaper, the
packing side of task A with the same box-clique constraints.  What this task establishes is the
infrastructure: a certificate may now carry cliques, both checkers verify them exactly, the
reduction is in Lean, and the LP driver can price, box, certify and write them.

## Reproduce

```sh
cd verify && cargo build --release && cd ..
./verify/target/release/verify certificates/s12_boxclique_demo_3.9318_N2000.txt 12 2000 8 0
python3 xcheck.py certificates/s12_boxclique_demo_3.9318_N2000.txt 2000 --all --n 12
./tests/rejection_tests.sh
python3 search/boxclique.py certificates/s12_lower_3.931795_sparse.txt demo224 --colgen 2 --per-round 2 --N 2000 --threads 4 --force 0.001
python3 search/boxclique.py certificates/s12_lower_3.931795_sparse.txt cl224_1992 --colgen 12 --per-round 4 --N 2000 --threads 4 --Dp 1992
```
