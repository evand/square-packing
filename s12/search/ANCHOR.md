# Anchor cliques on the cover side: the certificate object and the first leaf (task I)

Code: `certificates/FORMAT.md` ("Anchor cliques"), `verify/src/main.rs` (parser, Lemma-0 refusal,
the two exact predicates, sweep credit, witness placement), `xcheck.py` (the same in `Fraction`s),
`tests/rejection_tests.sh` (§21–27, 45 new checks), `lean/Sqpack/Basic.lean`
(`clique_of_anchors`, `packing_le_weight_anchor_cliques`), `search/anchorclique.py` (the family as
an LP column and as a certificate block), `search/anchorsep.py` (separation on the dual),
`search/anchorstress.py` (independent stress test), `search/branch.py` (`--cliques`).
Mathematics: `notes/clique-family.md`.  Packing-side measurements: `search/CLIQUE_CONTINUUM.md`.

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

## 4. The leaf run (`k = 4`, `r = 1`, `t = 3.98`)

Setup: a restart of the main tree's `t398hk4` run — same container `3120320/784000 = 3.98`, same
corner box `r = 1`, same net `N = 2000`, same margin `2·10⁻⁶`, same column checkpoint
(`runs/branch_t398hk4_cols.txt`, 6155 orbits) and probe, restricted master — **plus** anchor-clique
columns, seeded from the separation of §3 on that run's converged dual and re-priced every round.
The pure run it restarts converged to `12.000024 = 12·(1 + margin)`, i.e. exactly 12
(`search/BRANCH.md`); that is the number to beat.  Every LP value below is a **lower** bound on the
leaf's value that rises as rows are added and falls as columns are priced in, so an unconverged
value below 12 proves nothing and one above 12 is not final either.

The last and longest of the runs (`runs/branch_t398ik4g.log`; `runs/branch_t398ik4[d-h].log` are
the others):

| round | rows | point cols | clique cols (used) | clique weight | LP value | probe min | best ȳ(K) / ȳ(P_p) | time |
|---|---|---|---|---|---|---|---|---|
| `0` | 33,269 | 8322 | 460 (0) | 0 | `2.000016` | −1.500 | 1.0311 / 1.0311 | 0m23s |
| `1` | 39,269 | 8446 | 520 (7) | 1.0138 | `11.951010` | −0.853 | 1.4715 / 1.4667 | 0m53s |
| `2` | 45,268 | 8585 | 580 (1) | 2.0000 | `12.000024` | 0.375 | 1.1426 / 1.1364 | 0m38s |
| `3` | 51,268 | 8679 | 640 (9) | 0.2366 | `12.045093` | 0.884 | 1.1332 / 1.1332 | 1m16s |
| `4` | 57,265 | 8864 | 700 (11) | 1.0438 | `12.083954` | 0.923 | 1.2081 / 1.2081 | 1m46s |
| `5` | 62,934 | 9055 | 760 (1) | 4.0000 | `12.000024` | 0.000 | **1.0626 / 1.0082** | 1m24s |
| `6` | 68,934 | 9238 | 820 (13) | 3.8115 | `11.992229` | 0.498 | 1.0156 / 1.0132 | 3m36s |
| `7` | 74,934 | 9347 | 880 (14) | 2.8953 | `11.996570` | 0.938 | **1.3366 / 1.0000** | 3m58s |
| `8` | 75,289 | 9492 | 940 (1) | 4.0000 | `12.000024` | 0.000 | 1.0132 / 1.0132 | 3m12s |
| `9` | 81,079 | 9687 | 1000 (10) | 1.8342 | `12.342232` | 0.713 | 1.5361 / 1.5361 | 3m21s |

**What it shows.**  The mechanism runs end to end: the LP takes clique columns (up to 14 orbits at
once, carrying up to `3.8` of the total), the separator keeps finding violated ones, and the two
rows in bold are the interesting ones — at round 5 the best clique is violated by `+0.063` while the
worst coverage row is violated by only `+0.008`, and at round 7 by `+0.337` while **no** coverage
row is violated at all (`ȳ(P_p) = 1.0000`).  That is the non-Helly mass, on the cover side, in a
running LP.  What it does not show is convergence: the LP value oscillates between `11.95` and
`12.34` and the probe minimum between `0` and `0.94` as rows accumulate and columns are priced in.
**No verified certificate of value `< 12` was produced.**

The oscillation is the ordinary behaviour of this loop far from convergence (the pure reference took
18 rounds of 1–2 h each, its value climbing `11.996 → 12.000024` from below, before it settled); the
runs here are 10 rounds of 1–4 min from a restarted row set, i.e. an order of magnitude short — and
the LP value of a *restricted* master is an upper bound on the round's true optimum that moves as
the master's column set moves, so those excursions should not be read as "the leaf value moved"
either.

### The matched pair — the number that does not need convergence

The way to get a number out of a loop that has not converged is `search/CLIQUE_CONTINUUM.md`'s:
re-solve the **same** LP — same rows, same columns — with the clique columns removed, so the gain is
read off one row set instead of across runs.  `branch.py --matched` does that
(`runs/branch_t398ik4n.log`):

| round | rows | clique cols (used) | clique weight | clique LP | pure LP, same rows and columns | **gain** |
|---|---|---|---|---|---|---|
| `0` | 35,332 | 1060 (1) | 2.0000 | `2.000016` | `2.000016` | `0.000000` |
| `1` | 41,332 | 1120 (3) | 6.4000 | `11.600023` | `11.600023` | `0.000000` |
| `2` | 47,332 | 1180 (11) | 2.8298 | `11.811348` | `11.811739` | **`+0.000392`** |

So on the cover side, with 1180 anchor cliques offered and 11 of them carrying `2.83` of the total,
the LP is **0.0004 cheaper** than the same LP without them.  Against `≈ 0.10` on the packing side
(`search/CLIQUE_CONTINUUM.md` §7a, same leaf, same container).  Three rounds is a thin sample and
the row set is far from converged, but the sign and the order of magnitude are what they are, and
they say the cover side is not yet seeing what the packing side saw.  (Getting this right needed one
correction worth recording: the matched solve must run *before* the round's column generation — a
column added in between is marked active by the next solve, so the "pure" LP ran with 328 columns
against the clique LP's 168 and came out `0.05` *lower*, which is impossible for a minimisation over
a subset of columns.)

**Why the two sides can disagree.**  On the packing side a clique is a *cut*: `μ(K) ≤ 1` removes
mass from a fixed pose lattice, and 600 of them removed `0.10`.  On the cover side a clique is a
*column* of the same cost as the point column it contains (`K(p, A) ⊇ P_p`, Lemma 2), so it can only
help by the weight it saves on the poses of `{S : A ⊆ S}` outside `P_p` — a strip of width
`ε ≤ 1 − p_x`, which at the leaf's tight points (`p_x ≈ 0.99`) is `ε ≤ 0.01`.  The dual does see the
violation (`ȳ(K) = 1.054` on the converged pure dual, §3), but a violated dual constraint lowers the
primal only as far as the degeneracy allows, and this leaf is extremely degenerate — `search/BRANCH.md`
found its optimum sitting on an 80-point cover with dyadic weights and cost exactly
`18 − 4·1.5 = 12`.

### Two things worth recording from the engineering

* **The verifier became the bottleneck, and then did not.**  With 152 anchor cliques in the probe
  the sweep took `4300` CPU-s against `95` CPU-s for the same 1332-atom file without the block —
  520 s of a 976 s round.  Three fixes (a division-free per-cell prefilter in the sweep's own units,
  a fast-accept path for `contains` that tests corners only in a thin band around a piece's
  boundary, and hoisting the per-cell `Vec`s out of the closure — the last was the largest) brought
  it to `620` CPU-s with results unchanged wherever they can be compared exactly.
* **Row pruning cycles when the sweep doubles.**  A certificate with cliques is swept over
  `[0°, 90°)` instead of `[0°, 45°]`, so with `topk = 6` a round produces 12,000 witnesses instead
  of 5,000; the row set then hits `prune_at` every round, the prune drops rows whose slack is still
  large, and they come straight back — `new_rows = 12000` every round against `~1800` in the pure
  reference, and the LP value flat-lined at `11.94` with the probe minimum oscillating between
  `0.29` and `0.87` (`runs/branch_t398ik4d.log`).  `--topk 3 --prune-at 250000` restores the
  reference's row discipline.

## 5. Verdict

**The object is done and checked.  The dual violation is confirmed on the cover side — the family
does cut off the packing that pins the `k = 4` leaf at exactly 12.  The leaf is not closed, and the
matched pairs say that on the cover side the family is so far worth `10⁻⁴`, not the `10⁻¹` the
packing side measured.**

*Done and verified (not heuristic).*  A certificate may now carry anchor cliques; the verifier
refuses anything Lemma 0 does not certify, credits a cell only when one piece provably holds every
pose of it, and is bit-identical to the previous binary on every certificate without an anchor
block (10 checks, plain / witness / `TIGHT_DUMP` modes).  `xcheck.py` re-derives the whole thing in
`Fraction`s and agrees exactly — `5000011/5000000` against `10000022/10000000` — on an exhaustive
run over all 6000 bins of the shipped demonstration.  The Lean reduction is machine-checked
(`clique_of_anchors`, `packing_le_weight_anchor_cliques`, and — a gap the box-clique work left —
`packing_le_weight_regions_cliques` for a branch certificate that carries cliques), 0 sorries,
axioms unchanged.  136 rejection tests, 40,000 random cases of the two predicates against an
independent geometry implementation, 0 failures.  `certificates/s12_anchorclique_demo_3.9318.txt`
is a **load-bearing** demonstration: one point of the 224-point certificate replaced by its own
`K(p, A)`, verified at `N = 6000` and `N = 12000`, rejected the moment the clique's weight is
zeroed.

*Measured, heuristic (a float evaluation of an LP dual, no bound in either direction).*  On the
converged dual of the pure `k = 4` leaf — the packing of mass exactly `12.000000` that pins that
leaf — the best anchor clique has `ȳ(K) = 1.054176` while the worst coverage row has
`ȳ(P_p) = 1.003849`: **the family cuts off the packing that makes the leaf exactly 12, by `+0.05`,
and it is the only known describable family that does.**  That is the cover-side counterpart of task
G's packing-side `+0.024`–`0.10`.  Inside the running LP the same separation keeps firing (`+0.34`
at round 7 of the table above, with no coverage row violated at all).

*Measured, and the discouraging half of the answer.*  The **matched pairs** — the same LP with and
without the clique columns, so no convergence is needed — put the cover-side gain at `0.0000`,
`0.0000` and `+0.0004` over the three rounds measured, with 1180 cliques offered and 11 of them
carrying `2.83` of the total.  The packing side's `≈ 0.10` on this same leaf does **not** transfer,
at least not on these row sets.  §4 gives the reason to expect a gap: on the cover side the clique
is a column of the same cost as the point column it contains, so it can only pay for the weight
saved on `{S : A ⊆ S} \ P_p`, and at the leaf's tight points that strip is `0.01` wide.

*Not established.*  No verified leaf certificate of value `< 12`, and the leaf's value with cliques
was not pinned down: the loop's LP values (`11.95 … 12.34` over ten rounds) are restricted-master
upper bounds moving with the master's column set, not measurements of the leaf.  The loop is
compute-bound, not stuck — the pure reference needed 18 rounds of 1–2 h from the same start, and
the clique loop's row set was restarted several times while the verifier and the pricing were being
fixed.  A continuation is left running (`runs/leafk4q.sh`, from `runs/branch_t398ik4p_cols.txt` and
`runs/branch_t398ik4p_cliques.txt`: 10,395 point columns, 1,421 clique columns).  What it should be
asked next is not "does it close" but **"does the matched gain grow as the row set converges"** —
if it stays at `10⁻⁴` the family is not the lever on the cover side however violated the dual is,
and the design question moves to anchors that are *not* contained in a point clique (the two-anchor
corner families of `notes/clique-family.md` §3, which task G designed and did not build).

*The `1110` leaf was not run* (the brief makes it conditional on `k = 4` closing).  Its own pure
run stands at `11.967276` with probe minimum `0.94` and rising ~`0.0002`/round after 51 rounds
(`runs/branch_t398j1110.log`), i.e. it may well close without cliques.

*Deliverable 8 (the pure cover at `t = 3.99`) was not attempted.*  `branch.py --lam-lo 0 --lam-hi 0`
now pins the multipliers at zero, which is the pure cover LP with a vacuous trailer, so the same
driver can do it.

## Reproduce

```sh
cd verify && cargo build --release && cd ..
./tests/rejection_tests.sh
sh tests/bitid.sh <pre-change binary>                      # bit-identity without anchors
python3 search/anchorstress.py --n 40000 --samples 25
( cd lean && lake build )
M=/home/evand/math/square-packing/s12/runs
python3 search/anchorsep.py $M/branch_t398hk4_dual_it16.txt --pitch 0.01 --out runs/cq_seed_k4.txt
python3 search/anchordemo.py certificates/s12_lower_3.931795_sparse.txt out.txt --i 0 --D 1994
# the k = 4 leaf (`runs/` is gitignored, so the command in full; add --matched for the pair):
BRANCH_SOLVER=restricted BRANCH_RESTRICTED_PASSES=2 BRANCH_CQ_MAX=3000 \
python3 search/branch.py runs/branch_t398ik4p_probe.txt t398ik4q --k 4 --r 1 --N 2000 \
    --cols runs/branch_t398ik4p_cols.txt --warm-thr 3 --cliques 400 --colgen 400 --cq-want 60 \
    --cq-load runs/branch_t398ik4p_cliques.txt --prune-at 250000 --topk 3 --threads 8
```
