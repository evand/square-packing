# `s(13) = 4` from one weighted cover, without case analysis

*2026-09-12.  Companion to `search/RUNG2.md` (the method and the obstruction theorem),
`search/RUNG2_XCHECK.md` (the second checker) and `search/ZEROMARGIN.md` (the semantics and the
original primitives).  Everything below is reproducible from this repository; every number names
its source.*

---

## 1. Statement

Let `s(n)` be the side of the smallest square into which `n` unit squares can be packed, rotations
allowed.  `s(13) = 4` is a theorem of Bentz [B10, Theorem 9]: a weighted unavoidable set plus a
six-leaf case analysis, needed because the best *pure* unavoidable set of `[0,4]²` has 14 points
[DS7, Theorem 4] against the 12 the argument can afford — a deficit of 2
(`notes/proof-anatomy.md` §7.2).  This note re-proves the theorem from a single certificate, with
no case analysis at all.

> **Theorem 1.**  Let `P` be the 3,621 weighted points of
> `certificates/rung2/s13_closed_cover_4.txt`: rational coordinates with denominator `D = 1000`,
> rational weights with denominator `10⁹`, all in `[0,4]²`, of total weight
>
>     w(P)  =  2591194431/200000000  =  12.955972155  <  13.
>
> Then every **closed** unit square `Q ⊆ [0,4]²`, at every centre and every angle, satisfies
> `w(P ∩ Q) ≥ 1`, where a point of `P` on `∂Q` counts as lying in `Q`.

> **Corollary.**  `s(13) = 4`.

*Proof of the corollary from Theorem 1.*  Suppose 13 unit squares are packed, with pairwise
disjoint interiors, in a container of side `s' < 4`.  Rescale the whole picture by `4/s' > 1`.
The 13 squares become squares of side `4/s' > 1` with pairwise disjoint interiors inside `[0,4]²`,
and each of them *strictly* contains the concentric closed unit square.  Those 13 closed unit
squares therefore lie in the pairwise disjoint **interiors** of the dilated ones, so they are
pairwise disjoint as sets: no point of the plane lies in two of them.  Each captures weight `≥ 1`
by Theorem 1, and no point of `P` is counted twice, so `13 ≤ w(P) = 12.955972155`, which is false.
Hence `s(13) ≥ 4`.  Sixteen unit squares tile `[0,4]²`, so `s(13) ≤ 4`.  ∎

The tiling is what makes the bound sharp, and it is also what fixes the target: the machine must
beat 13, and it has the four units between 13 and 16 to lose.  §4 explains why that number 16 is
not an accident.

## 2. Semantics: why *closed*, and why it is the whole content

The certificate is a statement about **closed** unit squares with **closed containment**: `Q` is a
closed square of side exactly 1, `Q ⊆ [0,4]²` as closed sets, and a point of `P` lying on `∂Q`
counts.  This is not a convenience.

Every published lower-bound proof in this literature obtains the *exact* bound `s(n) ≥ k`, rather
than `s(n) > k − ε`, by the same rescaling.  DS7 §5 states it directly: shrinking an unavoidable
set by `(1 − ε/k)` gives a set of `n−1` points in a square of side `k − ε` that every unit square
meets *in its interior*, so `s(n) > k − ε` for every `ε > 0`, hence `s(n) ≥ k`.  Stromquist puts
the same move the other way round: rather than shrink the container, enlarge the squares being
packed.  Either way every certificate in this literature is really a statement about squares of
side `1 + ε` for all `ε > 0` in the *closed* container, and the limit of that family is precisely
"closed unit square, closed containment" (`notes/proof-anatomy.md` §1, §7.1).

The two conventions differ by exactly one thing, and it is decisive at `t = 4`:

* **Open semantics** (squares of side `1 − ε`, boundary does not count).  The sixteen squares of
  the `4 × 4` tiling, eroded a hair, are pairwise disjoint.  Sixteen disjoint squares each needing
  weight `≥ 1` force `w(P) ≥ 16`.  No cover of `[0,4]²` below 16 exists, so nothing below 16 can
  ever be proved this way.
* **Closed semantics** (squares of side `1 + ε` for all `ε > 0`).  The sixteen tiling squares meet
  along the grid lines `x, y ∈ {1,2,3}`, and a point on a shared edge is a shared resource: it is
  credited to both squares.  The "16" obstruction does not exist, and the true floor for a cover of
  `[0,4]²` is lower.  How much lower is known exactly: `search/COVER4.md` proves
  `COVER^closed(4) ≥ 24537607710/1999999999 = 12.2688038611`.

So the shipped cover, at `12.955972`, is within 5.6 % of the best any cover of `[0,4]²` can do,
and the room between the floor `12.2688` and the target `13` is all the room this proof has.
Note also what the same floor says about the neighbouring problem: `12.2688 > 12`, so **no cover
of `[0,4]²` proves anything about `n = 12`**.  This certificate is about 13 and only about 13.

Formally, the pose of a square is `(c, θ)` with `Q(c,θ) = c + R_θ[−½,½]²`, and `Q(c,θ) ⊆ [0,m]²`
iff `c_x, c_y ∈ [w/2, m − w/2]` with `w(θ) = |cos θ| + |sin θ|` (*admissibility*).  Both checkers
parametrise the angle by `u = tan(θ/2)`, so that `cos θ = (1−u²)/(1+u²)` and `sin θ = 2u/(1+u²)`
are rational functions of `u` and the entire computation stays in `ℚ`.  The support-function
characterisation of admissibility is proved in Lean as `sq_subset_box_iff`
(`notes/lean-zeromargin.md` §0).

## 3. The certificate

`certificates/rung2/s13_closed_cover_4.txt`, in the plain format of `certificates/FORMAT.md`:
container side `4/1`, coordinate denominator `D = 1000`, weight denominator `10⁹`, then `3621`
lines `X Y w` meaning the point `(X/1000, Y/1000)` carrying weight `w/10⁹`.  Every token in the
file is an integer; there is no floating-point number anywhere in it.

    sha256  ea303acea08cc17a13cecc24d3714c2df409f91eba048cd5546050ed064b53ed
    total   2591194431/200000000 = 12955972155/1000000000 = 12.955972155

pinned in `certificates/SHA256SUMS`.

How it was found is irrelevant to whether it is true, which is the point of the format, but for the
record (`search/RUNG2.md` §10): a cover LP over `[0,4]²` with column generation, seeded from
`search/closed4.py`'s best earlier cover and scaled by the exact rational factor `21/20`.  The one
substantive lesson from the search is a *soundness* lesson: every candidate cover before this one
was **invalid**, and invalid in the same place — the LP's own row lattice stepped over a family of
poses near `(3/2, 1461/2000, θ → 0⁺)`, where the best earlier cover captures only
`9420217/10000000 = 0.9420217`.  `search/family_rows.py` emits that family explicitly; with those
rows the LP produces a cover that is valid with about 4 % margin.  A dense float scan of the
shipped cover (`python3 search/closed4.py stress …`, which uses none of the box machinery of §5)
puts its minimum captured weight at `1.0313820`.

## 4. The method, and why the proof must be disjunctive

A checker for a statement of this kind subdivides the three-dimensional pose space
`(c_x, c_y, u)` into boxes with rational corners and discharges each box by an exact argument; if
it cannot, it splits the box and retries, down to a depth limit.  `VERIFIED` means every box was
discharged and **none** was abandoned — i.e. every admissible pose of the container lies in some
leaf whose certificate is exact.

The natural argument for a box is a **monotone witness certificate**: exhibit a fixed set `S ⊆ P`
with `w(S) ≥ 1` such that *every* `p ∈ S` lies in `Q(c,θ)` at *every* admissible pose of that box.
This is what `ZEROMARGIN.md`'s `CORE` and `P1`, and `RUNG2.md` §3's `ADM`, produce.  It is also
provably not enough:

> **Theorem 2** (`search/RUNG2.md` §2).  Let `m ≥ 4` and let `P ⊂ [0,m]²` be a finite weighted
> point set.  Suppose every leaf of a finite, closed, space-filling subdivision of the admissible
> pose space carries a monotone witness certificate.  Then `w(P) ≥ m²`.

The proof is an explicit LP dual.  For each grid tile `T_{ij} = [i,i+1]×[j,j+1]` with centre
`m_{ij}`, so that `Q(m_{ij}, 0) = T_{ij}`, fix signs `σ_x(i), σ_y(j) ∈ {±1}` that are `+1` at
index `0`, `−1` at index `m−1`, and switch between two *interior* indices — possible exactly
because `m ≥ 4`.  Push the tile pose off in the chosen octant: `γ(t,s,θ)` is the pose with centre
`m_{ij} + (σ_x t, σ_y s)`, clipped to admissibility, at angle `θ`.  The leaves are finitely many
closed sets covering the image of `γ`, so *some* leaf contains `γ([0,η]³)` for some `η > 0`
(Lemma 0), and its witness set is forced into

    PIN(i,j,σ)  =  ⋃_{η>0} ⋂_{0 ≤ t,s,θ ≤ η} Q(γ(t,s,θ)) ,

an explicit subset of `T_{ij}`: the tile, minus the half of each boundary edge that the chosen
octant sweeps off, minus all four corners (Lemma 1).  The `m²` sets `PIN(i,j,σ)` are **pairwise
disjoint** (Lemma 2) — two can share only a boundary line, and there the one-sidedness conditions
are complementary, one of them strict — so `w(P) ≥ Σ_{ij} w(PIN(i,j,σ)) ≥ m²`.  ∎

Membership in `PIN` depends only on comparisons of `p_x` with `i, i+½, i+1` and `p_y` with
`j, j+½, j+1`, so it is constant on each cell of the half-integer arrangement, the quarter-integer
lattice has a representative of every cell, and `python3 search/rung2_bound.py dual --m 4` verifies
the pairwise disjointness exactly over all 289 of them (maximum multiplicity 1).  With
`--switch 0` — the switch at the wall adjacency, which Lemma 2 forbids — it prints multiplicity 2,
so the hypothesis `m ≥ 4` is being checked, not assumed.

**Consequence.**  At `m = 4`, no cover of `[0,4]²` with `w(P) < 16` can be certified by monotone
witness primitives at any depth, and that includes every cover with `w(P) < 13`.  A proof of
Theorem 1 therefore *must* be disjunctive almost everywhere: most boxes have to be discharged by an
argument of the form "at every pose of this box, either this witness set is inside the square, or
that one is", with different sets in different regions of the same box.

The shape of the required disjunction is also measured rather than guessed
(`RUNG2.md` §6.1).  At the container's worst monotone pose, `(½, 3/2, 0)`, the square is
`[0,1]×[1,2]` and its capture structure over a small box is a *sliding cut*: the row `y = 1` is
captured to the left of a threshold `ξ` and the row `y = 2` to the right of the *same* `ξ`, and
`ξ` sweeps all of `[0,1]` as the centre offset and the angle range over any neighbourhood of zero.
No subdivision controls that.  A single pivot (two regions) reaches only `0.845828` of the needed
`1`; cutting between every pair of consecutive `x`-coordinates (56 regions) reaches `1.013756`.
So the primitive has to support a chain of arbitrary length — which is cheap, because the regions
are nested: they are the sign pattern of one monotone family of polynomials.

Both checkers implement a disjunctive primitive of this shape, by different means, and both
*measure* Theorem 2 from the other side: the disjunctive rule carries **65 %** of the non-empty
leaves in each of them (§5).  `zmcheck --nodisj` leaves 3,473 boxes uncertified in exactly the band
the disjunction closes (`notes/review-2026-09-12.md`).

## 5. The two checkers

Theorem 1 is verified twice, by programs written from the statement rather than from each other,
sharing no code, no subdivision rule and no primitive set.

|  | `search/zeromargin.py` | `verify2/zmcheck` |
|---|---|---|
| language, arithmetic | Python, `fractions.Fraction` on the whole load-bearing path; floats only as pre-filters | Rust, exact `i128` on the whole load-bearing path; refuses any box whose scale could overflow (`4 d_u + d_xy > 74`) |
| primitives | `EMPTY`, `CORE`, `P1`, `ADM`, `MIX`, `TRI`, `CHAIN` — seven tests, implemented separately | one polynomial class `P₄` (affine in the centre, degree ≤ 4 in `u`), one exact-maximum routine, one inference rule (nonnegative combinations of a node's hypotheses); every primitive is a special case |
| disjunction | `CHAIN`: a monotone chain of pivot inequalities partitions the box into `k+1` nested regions, each with its own witness set; two chains multiplied where two cuts slide independently, with provably empty product regions | recursive sign splitting with a closure step |
| domain | **reduced**: `c_x ∈ [0,4]`, `c_y ∈ [0,2]`, `u ∈ [0,½]` — legitimate only because the checker first verifies **exactly** that `P` is invariant under `x ↦ 4−x` and `y ↦ 4−y`, and refuses to run reduced otherwise | **full**: `c_x, c_y ∈ [0,4]`, `u ∈ [0,1]`, i.e. `θ ∈ [0°,90°]`.  No symmetry is assumed or checked |
| census | `boxes 16872, max depth 13` (limit 18): `ADM 2867  CORE 0  P1 0  MIX 0  CHAIN 5320  TRI 0  EMPTY 3449  UNCERTIFIED 0` | `boxes 30258, max depth 10` (limit 18): `ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0` |
| disjunctive share | `CHAIN` 5,320 of 8,187 non-empty leaves = 65 % | `DISJ` 9,477 of 14,591 non-empty leaves = 65 % |
| wall time | 6,138 s (1 h 42 min), 8 processes | 1,065 s (17 min) at 8 threads; 4-thread figure in `search/S13_WRITEUP.md` |
| determinism | re-run from the committed path: `done in 5828s`, census identical to the last box | census identical at every thread count tried |

Neither run reaches its depth limit, so in both cases the subdivision terminated on its own rather
than being truncated.  `CORE`, `P1` and `MIX` are 0 in the Python census because `ADM` subsumes all
three, and `TRI` is 0 because a weighted cover has no unit-weight triangles to apply it to.

Two further agreements are worth more than the censuses, because they are pointwise and exact.
`zmcheck pose` computes the captured weight of one rational pose with no boxes, no subdivision and
no primitives at all — just the four inequalities per point — and reproduces every digit of all
four exact `Fraction` values that `RUNG2.md` publishes, including `970282351/1000000000` at
`(3/2, 1461/2000, u = 1/40000)` with the same count of 175 captured points, a pose at which one
point sits exactly on `∂Q` (`RUNG2_XCHECK.md` §0).  Two implementations of the *semantics* agreeing
to the last digit of a nine-digit rational on a boundary-touching pose is a stronger check than any
lemma below it.  And `search/zeromargin_stress.py` re-verifies each Python leaf's own recorded
witness in floats through no shared code path: 40 sampled admissible poses on each of the 11,636
leaves, 0 failures, plus 200,000 random `P1`, 68,545 core-lemma, 72,175 triangle, 2,713 `ADM` and
220,000 `CHAIN` instances, all 0 failures (`certificates/rung2/README.md`).

## 6. Lean

`lean/Sqpack/ZeroMargin.lean` (970 lines, Lean 4 + Mathlib `v4.33.1`) formalises the **soundness
of every primitive both checkers rest on**, stated over `ℝ` with hypotheses that are exactly what
the code tests:

* `sq_subset_box_iff` — `sq c θ 1 ⊆ box m ↔ Adm m c θ`, both directions, all four sides: the
  definition of admissibility the whole pipeline uses;
* Lemma A (monotone corners) in eight-corner form, matching `_adm_cond_ok`'s four conditions and
  their four distinct corners;
* Lemma B — the bounds are polynomials in `u`, with the degree-4 coefficient rows machine-checked;
* Lemma C — the exact maximum over a bin, both branches, with the Bernstein branch proved by
  explicit degree-4 identities;
* Lemma E — the maximum is attained, so the checker's `_gmax ≤ 0` test is an *equivalence*, not
  merely sufficient;
* Lemmas F, G, H and the chain covering — the disjunctive primitive of §4;
* `clip_bin_no_loss` — the completeness fix of `RUNG2.md` §4.6.

**0 `sorry`.**  `lean/Axioms.lean` runs `#print axioms` on 36 theorems and every one reports
`[propext, Classical.choice, Quot.sound]` — Lean's standard axioms, nothing else.  Build:
`Build completed successfully (8710 jobs)` (`runs/lean_build_2026-09-12.log`).
`notes/lean-zeromargin.md` maps each theorem to the Python function it covers, and records three
discrepancies the formalisation found — including a sign error in the *task brief's* statement of
Lemma G (the note and the code were right).

**What is not in Lean, plainly.**  The subdivision, the control flow that claims the leaves exhaust
the pose space, the weight bookkeeping, the certificate parser and the float pre-filters — in
*either* checker.  The Lean covers the primitives; it does not cover the programs.  See §9.

## 7. Rejection tests

`tests/rung2/rejection_tests.sh`, 23 checks, `23 passed, 0 failed, 0 panics`
(`runs/rung2_rejection_2026-09-12.log`).  A checker that never says no is worthless, so each input
below has a required verdict and the script fails if it is not produced.

* **Baseline.**  The unmodified certificate over the hard band `c_x ∈ [0.5,0.6]`, `c_y ∈ [1,2]` —
  which contains the pose `(½, 3/2, 0)` of §4, where the certificate has least room — gives
  0 uncertified boxes.
* **Five mutations and their controls.**  One point's weight cut 1 % (still valid); every weight
  cut 5 % (a violating pose is *exhibited*, exact capture `0.997505942 < 1`); one point deleted
  (still valid); all weights halved (violation); the set scaled by `0.995` (violation at
  `(7/2,7/2,0)`, where it then captures exactly nothing); one point moved by `0.01` (still valid);
  the whole set translated by `0.01` (violation).  A *violation* is exit 1 with `*** VIOLATION` —
  a disproof, not a refusal.
* **Two invalid covers from the project's own history**, at their exact violating weights
  `970282351/10⁹` and `9420217/10⁷`, both at the `(3/2, 1461/2000, u = 1/40000)` pose of §3.
* **Nine malformed inputs**: negative weight, a point outside the container, a truncated point
  list, `s_den ∤ s_num·D`, trailing data after the last point, `W = 0`, a non-integer header, a
  missing file, `u ∉ [0,1]` in pose mode.  Each must exit 2 with `ERROR:` and **no verdict word**:
  refusing to judge is the correct response to garbage, and a *panic is not a rejection* — any
  panic fails the whole script.

The mutation sweeps use the restricted band rather than the whole container so that the suite runs
in about four minutes; a restricted sweep prints `PARTIAL SWEEP` and can never print `VERIFIED`,
so those tests compare uncertified-box counts instead of verdicts.  The two historical-cover tests
need a development file kept outside the repository and print `skip` when it is absent (20 passed).

## 8. Reproduce

```sh
./verify.sh          # everything, including the sweep below and the 23 rejection tests

# the Rust checker alone: full domain, ~17 min on 8 threads
(cd verify2 && cargo build --release)
verify2/target/release/zmcheck cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --threads 8
./tests/rung2/rejection_tests.sh

# one exact rational pose, no boxes and no primitives at all
verify2/target/release/zmcheck pose certificates/rung2/s13_closed_cover_4.txt \
        --x 7/2 --y 7/2 --u 0

# the Python checker: reduced domain, 1 h 42 min on 8 processes
python3 search/zeromargin.py cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --nproc 8 --disj --chain-from 0 --dump runs/leaves.txt

# Theorem 2's dual, exhaustively verified over a complete cell system (seconds)
python3 search/rung2_bound.py dual --m 4

# the Lean
cd lean && lake build          # 0 sorries; lean/Axioms.lean prints the axioms of 36 theorems

sha256sum -c certificates/SHA256SUMS
```

## 9. What would make this wrong

The corollary in §1 is a three-line argument a referee can check by hand; everything else rests on
Theorem 1, i.e. on the two checkers.  In decreasing order of how much it would worry me:

1. **A shared error in the *statement*, not the code.**  Both checkers were written from the same
   note.  If the semantics encoded in both — closed square, closed containment, admissibility as
   `c ∈ [w/2, m − w/2]²`, `u = tan(θ/2)` covering `θ ∈ [0°,90°]` — is wrong or incomplete, two
   independent implementations agree and are both wrong.  This is the residual risk that
   independence does *not* address, and the one I would attack first.  Against it:
   `sq_subset_box_iff` proves the admissibility characterisation in Lean from the definition of the
   square rather than asserting it, and `zmcheck pose` evaluates the raw membership inequalities
   with no admissibility reasoning at all and agrees with the sweep.
2. **The exhaustiveness of a subdivision.**  Neither checker's control flow is formalised; a
   recursion that could drop a box — a split not covering its parent, a root grid not covering the
   domain, a worker's boxes lost — would report `VERIFIED` over a hole.  Mitigations, not proofs:
   the two subdivisions are completely different (16,872 boxes to depth 13 against 30,258 to depth
   10, different root grids, different split rules), one runs the full domain, and the Rust
   checker's census is invariant under thread count — which is what a lost-work bug breaks first.
3. **The symmetry reduction.**  If `zeromargin.py`'s exact invariance check were wrong, its result
   would cover only part of the pose space.  This is why the second checker runs the **full**
   domain with no symmetry at all; the Python run is corroboration, not the proof.
4. **The certificate file**, which is pinned by sha256 and was re-run from the committed path.  A
   file with a different hash proves nothing here.
5. **Overflow in the Rust.**  `i128` is finite; the checker carries an a-priori bound and
   **refuses**, with a message, any box deeper than it allows, so no verdict can depend on
   wrap-around.  The deepest box in any run reported here is depth 10.
6. **A wrong reading of the literature.**  `s(13) = 4` is Bentz's theorem [B10]: this is a new
   proof, not a new result.  `notes/proof-anatomy.md` §1 quotes the sources on the rescaling step
   the corollary uses.

What would *not* break it: a bug in the LP search that produced the cover — the file is checked
against the definition, not against the process that found it.  Nor would anything in the `n = 12`
work here: the two results share machinery and nothing else, and `COVER^closed(4) ≥ 12.2688 > 12`
says the object of §3 can never speak about twelve squares.

## References

* **[B10]** W. Bentz, *Optimal packings of 13 and 46 unit squares in a square*, Electron. J.
  Combin. **17** (2010), #R126.
  https://www.combinatorics.org/ojs/index.php/eljc/article/view/v17i1r126 — Theorem 9 is
  `s(13) = 4`, by a weighted unavoidable set plus a six-leaf case analysis, dissected in
  `notes/proof-anatomy.md` §2.3.
* **[DS7]** E. Friedman, *Packing unit squares in squares: a survey and new results*, Electron. J.
  Combin. Dynamic Survey **DS7** (last revised 2009).
  https://www.combinatorics.org/ojs/index.php/eljc/article/view/DS7 — §5 states the rescaling
  argument of §2; its pure unavoidable sets for `[0,4]²` cost 14.
* `search/RUNG2.md` — the method, Theorem 2 (§2), the `ADM` and `CHAIN` primitives (§3, §6), the
  measurements (§4) and how the cover was built (§10).
* `search/RUNG2_XCHECK.md` — the second checker: its single polynomial class, single exact-maximum
  routine and single inference rule, and every primitive derived as a case of them.
* `search/ZEROMARGIN.md` — the semantics (§1) and why the angle-net verifier of `verify/` provably
  cannot check a cover at the container itself.
* `search/COVER4.md` — the exact lower bound `12.2688038611` on any closed cover of `[0,4]²`.
* `notes/lean-zeromargin.md` — the Lean statements, one per primitive, and what is out of scope.
* `notes/proof-anatomy.md` — §1 (how "side ≥ m" rather than "> m − ε" is obtained), §2.3 (Bentz's
  `n = 13` case analysis), §7.1 (which semantics is provably blocked).
