# A machine-checked lower bound for s(12)

[![verify](https://github.com/evand/square-packing-12/actions/workflows/verify.yml/badge.svg)](https://github.com/evand/square-packing-12/actions/workflows/verify.yml)

Readable write-up with the point diagram: **https://evand.github.io/square-packing-12/**

`s(n)` is the side of the smallest square into which `n` unit squares can be packed, with
rotations allowed.  `s(12) = 4` is conjectured but **open**: `s(13) = 4` is proved (Bentz 2010)
and `s(11) = 3.877083…` (Trump 1979), so the largest `n` with `s(n) < 4` is either 11 or 12,
and deciding which *is* the open problem.

This repo contains an exact certificate and its verifier for

> **s(12) ≥ 15680/3951 = 3.968616…**

improving the previous published bound `2 + 4/√5 = 3.788854…` (Stromquist 2003, inherited from
n = 11 by monotonicity; listed as the record for n = 11–12 in Friedman's Dynamic Survey DS7,
Table 2).  The upper bound remains the trivial `4`, so the gap narrows from
`[3.788854, 4]` to `[3.968616, 4]`.

It also contains a second, independent result produced by the same machinery: a **case-free
machine proof of `s(13) = 4`** — one weighted closed cover of `[0,4]²` of total weight
`12.955972 < 13`, verified exhaustively at margin zero by two checkers that share no code, where
Bentz's 2010 proof needs a six-leaf case analysis.  See
[`s(13) = 4` without case analysis](#s13--4-without-case-analysis) below.

Run everything with:

```sh
./verify.sh          # builds the verifier and re-checks every certificate
```

## The method

A **weighted unavoidable set** is a finite set of points with non-negative weights such that
*every* closed unit square inside the container captures total weight ≥ 1.  Twelve squares with
pairwise disjoint interiors capture ≥ 12 in total while no point is counted twice, so a
certificate of total weight < 12 rules them out.  This is the fractional (LP) relaxation of the
classical *unavoidable set of points* method — see credits below.

The main certificate has 1736 points and total weight `11.9738036`.  Two more are shipped:
764 points proving `s(12) ≥ 980/247 = 3.967611`, and a 224-point one at the earlier bound
`3920/997 = 3.931795` (the original 788-point certificate for that bound is kept for the record).

Weaker but purely combinatorial certificates are included too, because they are **uniform** —
every point has the same weight, so no fractions are involved and the statement is one sentence:

> **81 points** in `[0, 35/9]²` such that every closed unit square inside, at every angle,
> contains **at least 7** of them — twelve disjoint squares would need 84.

So `s(12) ≥ 35/9 = 3.888…`, already above Stromquist's `3.7889`, by a plain hitting-set
argument.  `search/uniform/UNIFORM.md` has the family for `k = 1 … 8` points-per-square
(`certificates/s12_uniform_<k>of<m>_<s>.txt`); they all share one describable skeleton — a `#`
of lines one unit from each wall plus a small ring at the centre — and the original 56-point,
5-per-square set at `19/5` is kept as `s12_56points_3.8.txt`.

## What is verified, and how

| | |
|---|---|
| `lean/` | Lean 4 + Mathlib formalisation of the reduction step (weighted set of total weight `W` ⟹ at most `W` squares), including the rescaling lemmas, the branch/clique variants, and the wall-strip chord lemma of `notes/chord-lemma.md` (at most 3 squares of a packing centred within 1 of a wall, for container side ≤ 4). **0 sorries**; axioms are only `propext`, `Classical.choice`, `Quot.sound`. |
| `verify/` | Exact `i128` verifier. Checks the covering property over the **entire continuum** of placements — no sampling. Angles are enumerated as rational rotations `θ_k = 2·arctan(k/N)` (so all trigonometry is rational); a unit square at any angle in `[θ_k, θ_{k+1}]` contains the concentric square of side `σ_k = 1/(cos δ + sin δ)` at angle `θ_k`, and for each such angle the minimum over all centres is computed exactly by an arrangement sweep. Verified at `N` = 6000 and 12000; at `N` ≤ 4000 the net's `σ`-shrink exceeds this certificate's slack and the verifier correctly refuses it (see `VERIFICATION.md`). |
| `tests/` | 42 rejection tests: mutated, malformed and adversarial certificates that the verifier must refuse.  Writing them found two soundness bugs in the verifier (neither affecting the shipped certificates); see `VERIFICATION.md`. |
| `xcheck.py` | Independent re-implementation in exact Python, re-derived from the definition rather than from the Rust; checks **every** angle bin (2486 bins at `N` = 6000 for the main certificate) and agrees with the Rust verifier bin by bin. |
| `verify2/` | Exact `i128` **zero-margin** checker `zmcheck`, for covers of the closed container `[0,4]²` where the margin is exactly 0 and the angle-net verifier above provably cannot work (`search/ZEROMARGIN.md` §1).  Adaptive subdivision of pose space `(c_x, c_y, u = tan(θ/2))` with one polynomial class, one exact-maximum routine and one inference rule; no sampling, no symmetry reduction.  This is what checks the `s(13)` cover; `search/RUNG2_XCHECK.md`.  23 rejection tests in `tests/rung2/`. |
| — | A third check, a dense float scan over 181 angles spanning the full 0–90° range (~500k centres each, not using the symmetry reduction), returns the same minimum. |

The `D4` symmetry of the point set — which is what reduces angles to `[0°, 45°]` — is itself
checked exactly by the verifier.

See `certificates/FORMAT.md` for the file format and the closed-square convention, and
`VERIFICATION.md` for the full log.

## Reproducing the certificate, not just checking it

`verify.sh` checks the shipped certificates.  To regenerate one from nothing:

```sh
# 1. LP with cutting planes over a rigorous cell decomposition of placement space.
#    Every LP iterate is already a valid certificate; this writes runs/cert_<tag>.txt.
#    (s, cell size `fine`, cell side `eta`, angle width `dt`, time limit, tag)
mkdir -p runs
python3 search/lp_search.py 3.92 0.005 0.005 0.005 7200 mytag

# 2. Scale the finished certificate up to its critical container size.
#    This step alone moved the first bound from 3.92 to 3.931795.
python3 search/scale_to_critical.py runs/cert_mytag.txt --n 12 --N 6000

# 3. Tighten: re-optimise the weights with the exact verifier as separation oracle,
#    sparsify, and grow the point set by column generation at a larger container.
#    This is what moved the bound from 3.931795 to 3.968616; see search/TIGHTEN.md.
python3 search/tighten.py --help
```

Step 2 is pure arithmetic: the integer coordinates never change, only the denominator `D`,
since scaling the whole picture by `λ` is exactly `D → D/λ`, `s → λs`.  A certificate
therefore proves a *family* of bounds and the best one is at the critical `D`; see
`certificates/FORMAT.md`.

Step 1 with a time limit is not reproducible bit for bit (it stops on the wall clock).  With
`--iters N --seed S` instead of a time limit it is a pure function of its arguments — two runs
give byte-identical certificates — and `--dump-lp` archives the exact final LP so the weights
can be regenerated from it with `--resolve`; see `search/REPRODUCIBILITY.md`.  The shipped file
predates those options.  The shipped file is pinned by `certificates/SHA256SUMS`, and what it asserts
is checked independently of how it was found — which is the entire point of the format.

The packing search that looked for a counterexample from the other side (L-BFGS + basin
hopping, validated against `s(5)`, `s(10)`, `s(11)`) is `search/pack_src/main.rs`.

## A bound for s(11) as well

The same pipeline run with `n = 11` gives **`s(11) ≥ 3040/797 = 3.814303…`**, improving
Stromquist's `2 + 4/√5 = 3.788854` (2003).  Certificate `certificates/s11_lower_3.8143.txt`
(680 points, total weight `10.8146708 < 11`), verified at `N` = 6000 and 12000 and by
`xcheck.py --all --n 11`; it is at its critical container.  It is a separate certificate:
`s(11) ≤ s(12)`, so the `n = 12` files say nothing about `n = 11`.  The LP crosses 11 at
about `3.815`, so for `n = 11` the pure method stops ~0.06 below the conjectured value
`3.877083` (Trump's packing); details and the crossing bracket in `search/N11.md`.

## Limits of the method

By LP duality the least possible certificate weight at container side `s` is at least the
**fractional packing number** `ν_f(s)`, so the method proves `s(12) ≥ s` only if `ν_f(s) < 12`.
**This ceiling is now pinned rigorously.**  An explicit fractional packing of mass `12.0282` (earlier `12.00823`)
at `s = 399/100`, certified in exact rational arithmetic (`search/DUAL_EXACT.md`,
`search/dual_exact.py`, `search/CLIQUE_CONTINUUM.md` §3), shows that no cover of weight `< 12` exists at any `s ≥ 3.99`; with
the shipped certificate, the ceiling `s*` of this entire family of arguments lies in
**`[3.968616, 3.99)`**.  It is probably a little below 3.99: non-converged runs give
`L(3.98) ≥ 11.918` and `L(3.97) ≥ 11.807` (`search/DUAL.md`), and column generation at
`3.9696` no longer gets below 12 (`search/TIGHTEN.md`).  The same measure shows that a
closed-semantics cover of `[0,4]²` — the object a limit argument at `s = 4` would need —
costs at least `12.008`; the best explicit one found costs `12.51` (heuristic,
`search/CLOSED4.md`), with 91 % of its weight on the grid lines `x, y ∈ {1,2,3}`.

Closing the remaining gap to 4 therefore needs case analysis layered on top of a certificate,
in the style of Bentz's `s(13)` proof (`notes/proof-anatomy.md` dissects those proofs).  What
the LP contributes to that is the size of the excess budget: between `0.008` and about `0.5`,
against the 3 units that unit-weight point sets leave for `n = 12`.  `TODO.md` separates what
is established from what is only suggested.

Separately, an extensive search for a packing of 12 unit squares into a square of side < 4
(L-BFGS + basin hopping, validated by reproducing `s(5)`, `s(10)`, `s(11)` to 5 decimals) found
nothing below 4; every run collapsed to the compressed 4×4 grid.  Code in `search/pack_src/`.

## Beyond the ceiling: branch certificates (2026-08-28)

The first step past the pure method is now built and measured (`search/BRANCH.md`).  A **branch
certificate** carries a region — the four corner boxes `[0,r]²` — a multiplier `λ` (or one per
box) and an occupancy count `k`, and asserts that squares centred in a box capture `≥ 1 + λ`,
all others `≥ 1`, with `W − λk < 12`; it then refutes every packing with exactly `k` squares
centred in the boxes, and the five (or, per box, sixteen) leaves together refute all of them.
The format is in `certificates/FORMAT.md`, the verifier and `xcheck.py` check it exactly, the
reduction is in Lean (`packing_le_weight_region`, `packing_le_weight_regions`), and
`search/branch.py` produces the certificates with the verifier as separation oracle.

What it gave at `s = 3.98`, where the pure cover costs ≈ 12.02: the leaves `k = 0, 1, 2` close
with room (`certificates/branch/`, verified; `verify_branch.sh`), the mixed leaf `1110` is at
≈ 11.75–11.97 and still rising, and the **all-corners leaf sits at `12.000 ± 0.001`**, straddled
by the row and column steps of the cutting-plane loop (an earlier claim that it was *exactly* 12 with
a dyadic grid cover was withdrawn — that run never produced a valid cover; see `search/CLIQUE.md`).
Its dual puts a full unit in each corner, 3.75 on the eight wall slots and 4 in a tilted interior
ring, and violates *clique* constraints by 0.5: rotated squares are not a Helly family, so the
point-cover LP misses valid inequalities `Σ_{S∈K} y_S ≤ 1` over pairwise-intersecting pose sets
`K` with no common point (`search/CLIQUE.md`, `search/CLIQUE_CEILING.md`).  Clique columns are
now a verifiable certificate object (`search/BOXCLIQUE.md`: exact Rust and Python checks, Lean
reduction, a demonstration certificate), and an exact zero-margin checker for closed containers
re-proves `s(15) = 4` from Friedman's 14 points (`search/ZEROMARGIN.md`).  Pure covers, for the
record: `COVER(3.99) = 12.2009` (converged), `COVER(3.98) ≈ 12.02`, `COVER(3.975) ≈ 11.96`.

**At the container itself (2026-09-07/08).**  Measured at `t = 4` in closed semantics, the whole
family — corner branch, per-slot branch, anchor cliques, the chord lemma (now proved in Lean,
`notes/chord-lemma.md`) and line-chord count cuts — leaves the hardest leaf at `11.75–11.78`
(certified `11.776`) and the corner leaf at `11.996`, every margin shrinking under refinement,
and no certified fractional packing of mass `>= 12` either; the extremal object is always the
`4×4` grid smeared and tilted.  The full corner × slot tree has 4213 leaves.  Notes:
`search/T4SCREEN.md`, `T4LEAF.md`, `LEAF_CEILING.md`, `LINECUTS.md`, `WITNESS.md`,
`notes/branch-semantics.md`, `notes/review-2026-09-07.md`; the plan is in `TODO.md`.

**`s(13) = 4` without case analysis (2026-09-12).**  As a milestone for the `t = 4` verifier, the
same machinery re-proves Bentz's `s(13) = 4` from a single object.  It is a separate result and
has its own section below.

**What the `t = 4` numbers for `n = 12` are worth (2026-09-12).**  Read as covers over the
continuum, the sub-12 packing-side duals cost 18–20, because the cliques carrying them are
non-Helly (grazing tangencies) and no sound positive-volume rule can credit them; without cliques
the certifiable family sits at exactly 12 on the corner leaf and above 12 on the pure instance
(`search/HONEST.md`, `notes/review-2026-09-12.md`).  `TODO.md` has the current critical path.

## s(13) = 4 without case analysis

`s(13) = 4` is Bentz's 2010 theorem.  His proof is a weighted unavoidable set plus a six-leaf case
analysis on where a square may sit (`notes/proof-anatomy.md` §2.3).  He needs the case tree because
the best *pure* unavoidable set of `[0,4]²` has 14 points (DS7 Theorem 4) against the 12 the
argument can afford — a deficit of 2 (`notes/proof-anatomy.md` §7.2).  The machinery built here for
`s(12)` re-proves the result from a **single object**, with no case tree at all.

> **Theorem.**  `certificates/rung2/s13_closed_cover_4.txt` is a set of **3,621 points** of
> `[0,4]²` with rational coordinates (denominator `D = 1000`) and rational weights (denominator
> `10⁹`), of total weight
>
>     2591194431/200000000  =  12.955972155  <  13
>
> such that **every closed unit square contained in `[0,4]²`, at every centre and every angle,
> captures total weight `≥ 1`** — closed containment, a point on the boundary of the square counts.
> Consequently no 13 unit squares fit in a square of side `< 4`; and 16 unit squares tile `[0,4]²`,
> so **`s(13) = 4`.**

**The arithmetic, in full.**  This is the reduction of `certificates/FORMAT.md` ("What the file
asserts"), formalised in `lean/Sqpack/Basic.lean`.  Suppose 13 unit squares pack into a container
of side `s' < 4`.  Rescale the picture by `4/s' > 1`: the 13 squares become squares of side
`4/s' > 1` with pairwise disjoint interiors inside `[0,4]²`, and each one *strictly* contains a
concentric closed unit square.  Those 13 closed unit squares are therefore pairwise **disjoint**,
not merely interior-disjoint, so no point of the cover is counted twice.  Each captures weight
`≥ 1` by the theorem, so `13 ≤ 12.955972155`, which is false.  Hence `s(13) ≥ 4`; the `4 × 4`
tiling gives `s(13) ≤ 4`; so `s(13) = 4`, and the bound is **sharp** — the 16 tiling squares are
what make `4` attainable, and they are also what makes the weight `13` (rather than `16`) the
thing to beat.

The closed semantics is not a convenience, it is the content.  Every published `s(m²−3)` proof
works with squares of side `1 + ε` for all `ε > 0`, which in the limit is exactly "closed unit
square, closed containment" (`notes/proof-anatomy.md` §1, §7.1; DS7 §5).  In *open* semantics the
16 slightly eroded squares of the tiling are disjoint and the cost of any cover of `[0,4]²` is
`≥ 16`; in closed semantics the tiling squares share their boundary grid points, that obstruction
does not exist, and `12.96 < 13` is available.

**Why the proof has to be disjunctive.**  `search/RUNG2.md` **Theorem 1**: if every leaf of a
finite, closed, space-filling subdivision of the admissible pose space of `[0,m]²` (`m ≥ 4`)
carries a *monotone witness certificate* — a fixed point set of weight `≥ 1`, every member of
which lies in the square at *every* admissible pose of that leaf — then the total weight is
`≥ m²`.  The proof is an explicit dual: for each of the `m²` tiles of the grid, a leaf must
contain the germ of poses approaching the tile pose from one chosen octant, its witness set is
forced into an explicit set `PIN(i,j,σ)`, and the `m²` such sets are pairwise disjoint
(`rung2_bound.py dual --m 4` verifies the disjointness exactly over a complete cell system of 289
quarter-lattice representatives).  So **no cover of `[0,4]²` with `W < 16` — the rung-2 target
`W < 13` included — can be certified by fixed-witness primitives at any depth**, and a cover of
weight `12.96` must be certified *disjunctively* almost everywhere: "at every pose of this box,
either this witness set is inside the square, or that one is".  Both checkers below measure
exactly that, independently: the disjunctive primitive carries 65 % of the non-empty leaves in
each of them.

**The two checkers.**  Two exhaustive checkers, written from the statement rather than from each
other, sharing no code, no subdivision rule and no primitive set:

| | `search/zeromargin.py` | `verify2/zmcheck` |
|---|---|---|
| what | the original checker, `search/ZEROMARGIN.md` + `search/RUNG2.md` | an independent re-implementation, `search/RUNG2_XCHECK.md` |
| arithmetic | Python `fractions.Fraction` (floats only as pre-filters, which can lose a certification but never create one) | Rust exact `i128` on the whole load-bearing path; refuses any box whose scale could overflow |
| primitives | `EMPTY`, `CORE`, `P1`, `ADM`, `MIX`, `TRI`, `CHAIN` | one polynomial class `P₄`, one exact-max routine (Lemma D/D′) and one inference rule (Lemma I); every primitive is a special case |
| domain | **reduced**: `c_x ∈ [0,4]`, `c_y ∈ [0,2]`, `u = tan(θ/2) ∈ [0,½]`, legitimate because the checker first verifies **exactly** that the point set is invariant under `x ↦ 4−x` and `y ↦ 4−y`, and refuses to run reduced otherwise | **full**: `c_x, c_y ∈ [0,4]`, `u ∈ [0,1]` (`θ ∈ [0°,90°]`); no symmetry is assumed or checked |
| census | `boxes 16872, max depth 13` (limit 18) — `ADM 2867  CORE 0  P1 0  MIX 0  CHAIN 5320  TRI 0  EMPTY 3449  UNCERTIFIED 0` | `boxes 30258, max depth 10` (limit 18) — `ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0` |
| disjunctive share | `CHAIN` 5,320 of 8,187 non-empty leaves = 65 % | `DISJ` 9,477 of 14,591 non-empty leaves = 65 % |
| time | 6,138 s (1 h 42 min) on 8 processes | 1,064 s (18 min) on 8 threads, 2,015 s (34 min) on 4 |
| determinism | re-run from the committed path reproduces the census to the last box (`5828 s`, same leaf counts) | identical census at every thread count tried |

Neither run reaches its depth limit, so in both cases the subdivision terminated on its own —
which is what "exhaustive" means here: every admissible pose of the container lies in some leaf
whose certificate is exact, and **0 boxes are uncertified**.

Two agreements worth more than the box censuses, because they are pointwise and exact:
`zmcheck pose` computes the captured weight of one rational pose with no boxes, no subdivision and
no primitives — just the four inequalities per point — and reproduces all four exact `Fraction`
values `RUNG2.md` publishes, to the last of nine digits, including a pose where a point sits
exactly on `∂Q` (`RUNG2_XCHECK.md` §0).  And a dense float scan of the cover by
`search/closed4.py stress`, which uses none of the box machinery, puts the minimum captured weight
at `1.0313820` — a 3.1 % margin.

**Rejection tests.**  `tests/rung2/rejection_tests.sh`, 23 checks, all passing: the unmodified
certificate over the hard band `c_x ∈ [0.5,0.6]`, `c_y ∈ [1,2]` (0 uncertified); five mutations
and their controls — weights cut 1 % (still valid), all weights cut 5 % (a violating pose is
exhibited, capture `0.997505942 < 1`), one point deleted (still valid), all weights halved
(violation), the set scaled by `0.995` (violation), one point moved by `0.01` (still valid), the
whole set moved by `0.01` (violation); the two invalid covers from the development history, at
their exact violating weights `0.970282351` and `0.9420217`; and nine malformed inputs — negative
weight, a point outside the container, a truncated list, `s_den ∤ s_num·D`, trailing data, `W = 0`,
a non-integer header, a missing file, `u ∉ [0,1]` — each of which must exit 2 with `ERROR:` and
**no verdict word**, because a checker that gives a verdict on garbage is not refusing.  No panics:
a panic is not a rejection.

**Lean.**  `lean/Sqpack/ZeroMargin.lean` (970 lines) formalises the **soundness of every primitive
both checkers rest on**: `sq_subset_box_iff` (the support-function characterisation of
admissibility, both directions and all four sides), Lemma A in eight-corner form, Lemma B (the
degree-4 coefficient rows), Lemma C (both branches, Bernstein by explicit degree-4 identities),
Lemma E (the maximum is attained, so the `_gmax ≤ 0` test is an equivalence), Lemmas F–H and the
chain covering, and `clip_bin_no_loss`.  **0 `sorry`**; `lean/Axioms.lean` prints the axioms of 36
theorems and every one is `[propext, Classical.choice, Quot.sound]`.  Mathlib `v4.33.1`;
`notes/lean-zeromargin.md` maps each theorem to the line of Python it covers.

**What is not machine-checked, plainly.**  The Lean covers the primitives, not the programs.  The
**subdivision and the exhaustiveness argument of each checker — the control flow that claims the
leaves cover the pose space, the weight bookkeeping, the certificate parser, the float
pre-filters — are not in Lean**, in either checker.  What stands behind them is that two programs
with different subdivisions, different primitive sets and different domains (one reduced by a
symmetry it checks exactly, one not reduced at all) reach `0 uncertified` on the same file, agree
to nine digits on individual poses, and refuse 23 things they should refuse.  That is a strong
claim, and it is not the same claim as a formal proof.  `notes/s13-casefree.md` is the
self-contained write-up, with a "what would make this wrong" paragraph.

**Where the cover came from, and how good it is.**  The LP that built it is in `search/RUNG2.md`
§10: a cover LP over `[0,4]²` with column generation, plus the explicit row families of
`search/family_rows.py` — the poses that the LP's own row lattice steps over, and whose omission
made every earlier candidate cover *invalid*.  `search/COVER4.md` proves
`COVER^closed(4) ≥ 24537607710/1999999999 = 12.2688038611` exactly, so this cover is within 5.6 %
of the best possible and **no cover argument at `[0,4]²` can go below `12.2688`** — comfortably
under 13, and comfortably over 12, which is why the same object says nothing about `n = 12`.

**Reproduce.**

```sh
./verify.sh          # includes the zmcheck sweep and the 23 rejection tests

# either checker on its own:
(cd verify2 && cargo build --release)
verify2/target/release/zmcheck cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --threads "$(nproc)"                      # ~17 min on 8 threads
./tests/rung2/rejection_tests.sh                             # ~4 min

# the slow path: the Python checker, 1 h 42 min on 8 processes
python3 search/zeromargin.py cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --nproc 8 --disj --chain-from 0
```

## Credits and prior art

The unavoidable-point-set method is due to Göbel, and was developed by Stromquist, Friedman,
Kearney–Shiu, Nagamochi and Bentz; Friedman's Dynamic Survey **DS7** is the standard reference.
The idea of replacing an unavoidable *set* by LP-optimised *weights*, and mechanising the
verification in exact rational arithmetic, is due to the August 2026 work on `s(17)` by
**Sam Burns** and **Gustavo Massaccesi** (unrefereed blog posts); this repo applies that idea to
`n = 12`, which as far as I know had not been tried.  `CREDITS.md` has the full lineage.

What is new here beyond that: cutting planes generated over a rigorous *cell* decomposition of
placement space (so every LP iterate is a valid certificate), column generation for the point
locations, the exact arrangement-sweep verifier over rational rotations, the Lean formalisation
of the reduction, and the observation that a finished certificate can be **scaled up to its
critical size**, and — the step that moved the bound from 3.93 to 3.97 — column generation
with exact reduced-cost pricing and the exact verifier as the separation oracle (`search/tighten.py`).

## References

**Results this work depends on**

* F. Göbel, *Geometrical packing and covering problems*, in *Packing and Covering in
  Combinatorics* (A. Schrijver, ed.), Math. Centrum Tracts **106** (1979), 179–199.
  The unavoidable-set argument.
* W. Stromquist, *Packing 10 or 11 unit squares in a square*, Electron. J. Combin. **10**
  (2003), #R8.  https://www.combinatorics.org/ojs/index.php/eljc/article/view/v10i1r8 —
  the previous bound, `s(11) ≥ 2 + 4/√5`, which transfers to `n = 12` by monotonicity.
* W. Bentz, *Optimal packings of 13 and 46 unit squares in a square*, Electron. J. Combin.
  **17** (2010), #R126.  https://www.combinatorics.org/ojs/index.php/eljc/article/view/v17i1r126 —
  `s(13) = 4`, which makes `n = 12` the boundary case.  It is also the result **re-proved
  case-free** here: Bentz's Theorem 9 is a weighted unavoidable set plus a six-leaf case analysis
  (dissected in `notes/proof-anatomy.md` §2.3), where the certificate of
  `certificates/rung2/s13_closed_cover_4.txt` is a single object with no case tree.  The theorem
  is his; only the proof is new.
* E. Friedman, *Packing unit squares in squares: a survey and new results*, Electron. J.
  Combin. Dynamic Survey **DS7** (last revised 2009).
  https://www.combinatorics.org/ojs/index.php/eljc/article/view/DS7 — Table 2 is the only
  tabulation of lower bounds for `s(n)`; a maintained copy is at
  https://erich-friedman.github.io/papers/squares/squares.html .
* D. Ellsworth, *Squares in Squares*, https://kingbird.myphotos.cc/packing/squares_in_squares.html —
  current record packings (including Trump's `s(11) = 3.877083…`, 1979) with exact constants;
  the "previous record" claim above was checked against this and DS7.

**The mechanised-certificate line this work belongs to (all August 2026, all unrefereed)**

* S. Burns, *Proposing a better lower bound for n=17 square packing* (6 Aug 2026),
  https://sam-burns.com/posts/proposing-better-lower-bound-for-n17-square-packing/ —
  weighted atoms, exact rational verification.  `s(17) ≥ 4.4811`.
* G. Massaccesi, *Another better lower bound for n=17 square packing* (21 Aug 2026),
  https://gus-massa.blogspot.com/2026/08/another-better-lower-bound-for-n17.html , and
  *Linear programming for square packing*,
  https://gus-massa.blogspot.com/2026/08/linear-programing-for-square-packing.html —
  LP-optimised weights over a rational angle net.  `s(17) ≥ 4.5058`.  The direct ancestor of
  the method here.
* S. Fort, https://github.com/stanislavfort/17squares — `s(17) > 4.456575`; unweighted points,
  exact subdivision of pose space, CI re-verification.
* Mira, https://github.com/Mira-acc/17squares — `s(17) > 4.468292`; 16 unweighted points,
  exact dyadic subdivision of pose space, a triangle-piercing lemma, three independent checkers,
  rejection tests, and a write-up at
  https://github.com/Mira-acc/17squares/blob/main/paper/17squares-lower-bound.pdf .  Several
  verification practices here (rejection tests, self-describing certificate files, shipping
  checkers rather than claims) are taken from this repo.

The two 2026 families use different conventions — Fort and Mira use *open* squares and
pigeonhole, Burns, Massaccesi and this repo use *closed* squares with a concentric-shrink
argument — so certificates are not interchangeable between their checkers and ours without
reworking the disjointness step.  As of this writing none of the four `n = 17` results, nor this
one, appears in DS7 Table 2 or on Wikipedia's *Square packing* page.

**Further reading**

* M. J. Kearney and P. Shiu, *Efficient packing of unit squares in a square*, Electron. J.
  Combin. **9** (2002), #R14.  https://www.combinatorics.org/ojs/index.php/eljc/article/view/v9i1r14
* H. Nagamochi, *Packing unit squares in a rectangle*, Electron. J. Combin. **12** (2005),
  #R37.  https://www.combinatorics.org/ojs/index.php/eljc/article/view/v12i1r37 — weighted
  points, segments and areas as "resources", the closest classical precedent for weights.
* W. Bentz, *Optimal packings of 22 and 33 unit squares in a square*, arXiv:1606.03746,
  https://arxiv.org/abs/1606.03746 .
* Wikipedia, *Square packing*, https://en.wikipedia.org/wiki/Square_packing .
* `google-deepmind/formal-conjectures` issue #646, *Seventeen square packing problem*,
  https://github.com/google-deepmind/formal-conjectures/issues/646 — a Lean formalisation
  venue for statements in this area.

## Status

Not peer reviewed.  The "previous record" claim rests on DS7 plus David Ellsworth's current
record tables; if a better published bound for `n = 12` exists, I did not find it.
Independent checking is welcome, and is the point of the format being this boring.
