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
  `s(13) = 4`, which makes `n = 12` the boundary case.
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
