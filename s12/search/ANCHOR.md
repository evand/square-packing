# Anchor cliques on the cover side: the certificate object and the first leaf (task I)

Code: `certificates/FORMAT.md` ("Anchor cliques"), `verify/src/main.rs` (parser, Lemma-0 refusal,
the two exact predicates, sweep credit, witness placement), `xcheck.py` (the same in `Fraction`s),
`tests/rejection_tests.sh` (§21–27, 45 new checks), `lean/Sqpack/Basic.lean`
(`clique_of_anchors`, `packing_le_weight_anchor_cliques`), `search/anchorclique.py` (the family as
an LP column and as a certificate block), `search/anchorsep.py` (separation on the dual),
`search/anchorstress.py` (independent stress test), `search/branch.py` (`--cliques`).
Mathematics: `notes/clique-family.md`.  Packing-side measurements: `search/CLIQUE_CONTINUUM.md`.

*(status: work in progress — the leaf run of §5 is still going as this is written)*

## 1. What was built

**The object.**  A certificate may now carry an `anchors` block:

```
anchors A c                 # A anchors, c anchor cliques
anchorP X Y D               # anchor: the point (X/D, Y/D)
anchorS X0 Y0 X1 Y1 D       # anchor: the closed segment between two rational points
w_1 P_1                     # clique 1: weight numerator over W, number of pieces
piece a r f_1 … f_r         # { S : A_a subseteq S and S meets A_{f_i} for every i }
```

A clique is the union of its pieces, and the **only** thing trusted about it is Lemma 0's
hypothesis, which the verifier checks exactly: for every ordered pair of pieces, either the two
anchors intersect (exact rational segment intersection, closed) or one piece filters the other's
anchor.  Everything else — the weights, the counted-once rule, the reduction — is the clique
machinery that was already there for box cliques (`search/BOXCLIQUE.md`).

Two differences from box cliques are the point of the exercise:

* a one-piece clique with a **point** anchor and no filter *is* the point clique of `p`, i.e. the
  coverage constraint at `p`.  A box clique can never be that (`notes/clique-family.md` §1).  The
  anchor family therefore *strengthens* an existing row instead of adding one the LP ignores.
* the anchors do not refer to the verifier's angle net, so an `anchors` block carries no `N` and
  is meaningful at every `N` (a finer net credits more cells).  Box cliques are tied to their net.

**The two predicates** (both exact `i128`, both only ever *under*-crediting):

* `contains(anchor, cell, bin)` — every pose of the cell contains the anchor.  For a fixed centre,
  the intersection of the closed unit squares over the bin is the **exact bin core**
  `R_{θ_k}Q ∩ R_{θ_{k+1}}Q ∩ {x : dir(x) mod 90° ∈ [θ_k, θ_{k+1}] ⇒ |x| ≤ ½}`
  (`search/ZEROMARGIN.md` §2), which is convex, so it is enough to test the four corners of
  (anchor endpoint − cell): two rotated-square tests, a rational sector test and one squared norm.
  This is what makes a point anchor an *exact* point clique — an edge midpoint lies in every
  rotation of the square and the `σ_k`-square would lose it.
* `meets(anchor, cell, bin)` — every pose of the cell meets the anchor.  Every pose of the bin
  contains the concentric `σ_k`-square, so it is enough that the cell lie in `A ⊕ [−h, h]²`,
  `h = σ_k/2`: four half-planes from the anchor's bounding box and **two from the segment's own
  normal** — the axis whose omission is the recorded near miss of `notes/clique-family.md` §7.

A cell is credited a clique's weight iff **one piece holds every pose of the cell**; a cell that
only partly satisfies the predicates gets nothing and its LP witness is placed at a pose of the
cell outside the piece, so the cutting-plane loop cannot cycle on it.  Cell and anchors are rounded
**outward** to a `10⁻⁹` grid in the bin's frame before the exact tests, which keeps the `i128`
arithmetic far from overflow (and an overflow is an `ERROR`, never a verdict) and can only shrink
the credited region.

## 2. Verification of the verifier

**Bit-identity** (`tests/bitid.sh`, single thread, against the pre-change binary): identical output
on the main 1736-point certificate, the 764-, 224- and 56-point ones and the box-clique demo at
`N = 2000`; identical in witness mode (`topk = 6`) on the main certificate — stdout *and* the
3652-line witness file — and on the box-clique demo; identical on a per-box branch trailer in
witness mode; identical on the shipped `k = 1` branch certificate; identical in `TIGHT_DUMP` mode.

**Rejection tests**: `tests/rejection_tests.sh` is now **132 checks** (was 87).  The ones that
matter:

| test | expected | why |
|---|---|---|
| a deleted point replaced by its own point clique as a one-piece anchor clique | VERIFIED, min back to exactly `5/5` | the credit is *exact*, not approximate: the anchor piece covers precisely the poses the point covered |
| the same anchor displaced by `0.05` | REJECT | a verifier crediting cells that merely *meet* a piece would pass it |
| `K(p, A)` with the filter segment moved far away | REJECT | the `meets` predicate is really tested and is not vacuous |
| pieces whose anchors neither meet nor filter | ERROR | Lemma 0's hypothesis |
| filter naming its own anchor / a missing anchor, anchor outside the container, negative weight, no pieces, denominator 0, unknown keywords, block after the region trailer | ERROR | well-formedness, never a verdict |
| zeroed anchor weights | VERIFIED from the points alone | no weight leaks |
| box clique + anchor clique in one file | VERIFIED, weights add | the two blocks coexist |
| six of the files re-checked with `xcheck.py` | same minimum | independent implementation |

**`xcheck.py`** re-derives both predicates in `Fraction`s with the anchors inflated (exact `σ_k`,
exact anchor coordinates, no outward rounding), so it credits at least what the Rust credits.  On
the shipped demonstration certificate (§3a) an **exhaustive** run — all 6000 bins of the `N = 6000`
net, 1783 s — gives minimum covered weight `5000011/5000000 = 1.0000022`, exactly the Rust's
`10000022/10000000`; the two implementations agree cell for cell on a real anchor certificate.  One
real defect was found this way and fixed: the verifier's `y`-breakpoints contain `ylo` and `yhi`,
so its cells never stick out of the admissible range, while `xcheck`'s did — the anchor credit now
sees the clipped cell.  Before the fix the two disagreed (`4/5` against `5/5`) on the
point-clique-replacement test.

**Independent stress test** (`search/anchorstress.py`): random `(bin, cell, anchor)` cases; where
the exact predicates credit the cell, poses of the cell are sampled and tested directly from the
definition, by code sharing nothing with either checker (an explicit rotation for `contains`,
Liang–Barsky clipping for `meets`).  **0 failures.**

**Lean** (`lean/Sqpack/Basic.lean`, `lake build` clean, 0 sorries, axioms `propext`,
`Classical.choice`, `Quot.sound` — unchanged): `clique_of_anchors` is Lemma 0 (three-way case
split, no geometry), and `packing_le_weight_anchor_cliques` instantiates
`packing_le_weight_cliques` with it, so the reduction from an anchor-clique certificate to
`n ≤ Σ w + Σ v` is machine-checked end to end.  A gap left by the box-clique work is closed at the
same time: the leaf certificate is a *branch* certificate carrying a clique block, and that
combination had no theorem — `packing_le_weight_regions_cliques` now proves
`n + Σ_j λ_j·#{i : pose_i ∈ R_j} ≤ Σ w + Σ v`, which is exactly what such a file asserts.

## 3. The LP side

`search/anchorclique.py` turns `K(p, A)` into a column and into a certificate block;
`search/branch.py --cliques` prices them inside the cutting-plane loop.

* A candidate is a wall point `p` (Lemma 1: nowhere else can a clique beat the coverage row at
  `p`) with the perpendicular segment `A` of Lemma 2 at offset `ε = frac·(1 − dist)` and
  half-length `ρ = 1.1 ρ* + pad`, so `K(p, A) ⊇ P_p`: the clique column has the **same cost** as
  the point column of `p` and covers at least the same rows — it *dominates* it — and what it adds
  is the dual mass of `{S : A ⊆ S} \ P_p`.
* Membership of a row is decided with the same conservative predicates the verifier uses on a cell
  (`σ_k`-square for `contains`, hexagon for `meets`), so the LP never over-credits itself relative
  to the verifier.
* Pricing is on the interior-point dual, D4-averaged: the reduced cost of a clique orbit is
  `|orbit|·(1 − ȳ(K))`, so a column is worth adding iff `ȳ(K) > 1` — which is the same statement
  as "the dual packing violates the clique constraint `μ(K) ≤ 1`".

**The separation measurement that made the leaf run worth starting.**  On the *converged* dual of
the main tree's `k = 4` leaf (`runs/branch_t398hk4_dual_it16.txt`: 520 poses, mass `12.000000`,
maximum D4-averaged point coverage `1.003849`, i.e. the point loop is converged):

| object | value |
|---|---|
| best point row `ȳ(P_p)` over the wall grid | `1.003849` |
| best anchor clique `ȳ(K(p, A))`, `p = (0.99, 1.44)`, `ε = 0.0095`, `ρ = 0.0753` | **`1.054176`** |
| the extra, `ȳ(K) − ȳ(P_p)` | `+0.0503` |

So at the point where the pure leaf stalls at `12.000`, there is a clique the dual packing violates
by `5 %` while it violates no coverage row by more than `0.4 %`.  That is the cover-side
counterpart of task G's packing-side `+0.024`–`0.10`, and it is a **heuristic** number (a float
evaluation of an LP dual), not a bound in either direction.

## 3a. The demonstration certificate, and the one practical fact that makes the family work

`certificates/s12_anchorclique_demo_3.9318.txt` (`search/anchordemo.py`): the shipped 224-point
certificate with **one point deleted** and the anchor clique `K(p, A)` of Lemma 2 put at that point
with the same weight — `p = (1375, 3875)/1994`, wall distance `0.6896`, `ε = 0.2949`,
`ρ = 0.3109 > ρ* = 0.2808`.  Since `K(p, A) ⊇ P_p` the covering is at least as strong and the total
weight is unchanged, `11.9834372 < 12`.  It verifies at `N = 6000` and at `N = 12000` (minimum
covered weight `1.000002` at both — an anchor block refers to no angle net, unlike a box clique),
and with the clique's weight set to `0` the minimum drops to `0.880188`: **the clique is
load-bearing**, which the box-clique demonstration could not be.

Building it turned up the fact that decides whether the family is usable at all.  The first version
*failed*, at `0.958275` — exactly one clique weight short.  The reason: the sweep's cells are the
atoms' breakpoints, and a cell is credited only if it lies **wholly** inside a piece, so with `p`
gone from the atom list the cells straddle the boundary of `{S : p ∈ S}` and a whole band around it
loses the credit — precisely where the covering is tight.  The fix is to keep the anchors' own
points in the file as **zero-weight atoms** (allowed: weights need only be `≥ 0`), which puts the
cell boundaries back where the piece boundaries are.  `branch.py`'s export does this for every
anchor of every clique it writes, and `certificates/FORMAT.md` records it.

## 4. The leaf run

*(to be completed)*

## 5. Verdict

*(to be completed)*

## Reproduce

```sh
cd verify && cargo build --release && cd ..
./tests/rejection_tests.sh
sh tests/bitid.sh <pre-change binary>                      # bit-identity without anchors
python3 search/anchorstress.py --n 40000 --samples 25
( cd lean && lake build )
M=/home/evand/math/square-packing/s12/runs
python3 search/anchorsep.py $M/branch_t398hk4_dual_it16.txt --pitch 0.01 --top 400 --out runs/cq_seed_k4.txt
sh runs/leafk4.sh                                          # the k = 4 leaf at t = 3.98
```
