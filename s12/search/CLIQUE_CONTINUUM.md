# Cliques in the continuum: is the non-Helly mass certifiable? (task G, 2026-08-30)

Code: `search/clique_family.py` (the family and its geometry, separation), `search/clique_cert.py`
(exact rational certification of the clique property), `search/clique_continuum.py` (the
clique-strengthened packing LP on a fixed pose lattice), `search/clique_exact.py` (exact
certification of a clique-feasible measure).  Mathematics, lemmas and proofs:
`notes/clique-family.md`.  Runs: `runs/cq_*`, `runs/cqx_*`.

## Verdict

**GO — the clique lever is real in the continuum, it is about 14× the pure method's own excess,
and unlike task A's the measurement is stable under pose refinement.  It is not a proof, and from
the packing side it cannot be one.**

The headline pair, on the lattice that passes calibration and on one and the same pose set:

    pure LP            12.0115 → 12.0258        (target L(3.99) = 12.008: PASS at every stage)
    anchor-clique LP   11.8965 → 11.9127        (gain 0.104 – 0.120)

over eleven stages in which the pose set grows from 558 to 4739 orbits, each stage ending with no
violated anchor clique the separator can find.  The clique value stays pinned at `11.90 ± 0.006`
across that 8.5× refinement while the pure value on the very same poses climbs — flat, where task
A's unrestricted-clique ladder drifted by `0.26` over a comparable refinement.  The same
measurement in the `k = 4` corner leaf at `t = 3.98` gives the same gain, `0.103–0.105` (§7a).

Three things are new relative to tasks A and B.

1. **A certifiable family that cannot be dodged.**  For a point `p` and a nonempty compact convex
   `A`,
   `K(p, A) = {S : p ∈ S and S ∩ A ≠ ∅} ∪ {S : A ⊆ S}`
   is a clique for **every** `p` and `A` — no geometric lemma, no non-avoidance argument, two
   lines (`notes/clique-family.md` Lemma 0).  Membership is a geometric test on the pose, so a
   newly priced column is charged correctly: unlike task A's explicit cuts, column generation
   cannot step around it.  The family contains the exact point clique (`A = {p}`), which box
   cliques (task B) never can.
2. **The non-Helly mass survives, but shrinks by a factor of ten.**  On the certified extremal
   measure at `t = 3.99` (mass `12.008230754`): the max clique of the support's
   closed-intersection graph is `1.214` (task A: `1.327`), the best *anchor* clique is `1.057`
   but its `{S : A ⊆ S}` is a single pose, and the best anchor clique whose `A` is a genuine
   segment — so that `{S : A ⊆ S}` is a positive-measure pose region — is **`1.024`**.  So the
   region-shaped violation is `+0.024`, not `+0.33`; but `0.024` is still three times the pure
   method's whole excess of `0.008` at `3.99`.
3. **Where it can live is pinned exactly.**  A clique can beat the coverage constraint at `p`
   only if `p` is *strictly* within distance 1 of a wall (Lemma 1, sharp and proved in three
   lines), and the size of the gain is governed by an exact threshold
   `ρ* = ε·p_x/sqrt(1-p_x²)` (Lemma 2, proved and stress-tested).  The violated region cliques
   found sit at wall distances `0.55, 0.75, 0.93` — inside that band, as predicted — while the
   *tight points* of the pure measures sit right against the line: `(0.98962, 0.99835)` for the
   old certified measure and `(0.999939, 0.999939)` for the new one, i.e. wall distance `0.99` and
   `0.99994`, where the whole extra family is `0.13 %` and `4e-6` of the point clique.  The clique
   lever is therefore **not** at the tight points of the pure optimum: the LP has to move mass
   inward to reach it, and that is exactly why the cuts cost it something.

**What is not settled, and why the packing side cannot settle it.**  Every number below is an LP
value on a *restricted* pose set, hence a **lower** bound on the continuum clique-LP value `V`;
they rise as the pose set is refined.  A "go" verdict — `V < 12` — would need an **upper** bound
on `V`, and upper bounds live on the cover side (a clique certificate).  So this task establishes
that the certifiable clique family is worth `≈ 0.1` at `t = 3.99` on every pose set rich enough to
reproduce the pure value, and that the mechanism is understood and exactly certifiable; it does
not prove that the continuum value is below 12.  **The clique line is not dead, and the next
decisive experiment is on the cover side** (clique columns in the `k = 4` leaf at `3.98`,
`TODO.md` Phase 2), which is where the `TODO` already points.

## 1. The certifiable family

Full statements and proofs: `notes/clique-family.md`.  Summary:

* **Lemma 0 (validity).**  `K(p, A)` is a clique for every `p` and every nonempty compact convex
  `A`.  Two members of the first part share `p`; two of the second share `A`; a mixed pair `S`,
  `S'` has `S ∩ S' ⊇ S ∩ A ≠ ∅`.  The general form is
  `K(A_0, …, A_{m-1}) = ⋃_i {S : A_i ⊆ S and S meets A_j for every j > i}`; task B's box cliques
  are the case where all anchors pairwise meet, so no filter is needed.
* **Lemma 1 (where cliques can help at all; sharp, proved).**  If `1 ≤ p_x ≤ t-1` and
  `1 ≤ p_y ≤ t-1` then the point clique `P_p` is a **maximal** clique: no admissible pose outside
  it meets every admissible pose inside it.  Proof: the outward edge normal `n` of `S'` that
  separates `p` gives the axis-parallel unit square with corner `p` in the quadrant of
  `sign(n)`, which is admissible exactly under that hypothesis and is disjoint from `S'`.
  Measured: the extra volume is `0` at wall distance `1.000` and `5.7e-5` at `0.999`.
* **Lemma 2 (the exact threshold).**  For `0 < p_x < 1`, `0 < ε < 1-p_x`, and `A_ρ` the vertical
  segment `{p_x+ε} × [p_y ± ρ]`, every admissible unit square containing `p` meets `A_ρ` **iff**
  `ρ ≥ ρ* = ε·p_x/sqrt(1-p_x²) = ε·cot(arccos p_x)`; the worst pose is `θ = arccos p_x`,
  `c_x = w(θ)/2` — a square jammed against the wall.  Proved (convexity of
  `f(θ) = p_x/cos θ + ε cos θ - ρ sin θ - 1` on `[0, θ₀]`), reproduced to 12 digits by two
  independent numerical implementations, and stress-tested on 300,000 random instances (0
  failures in both directions).  When `ρ ≥ ρ*` the clique constraint **dominates** the coverage
  constraint at `p`.

**Exact certification.**  `search/clique_cert.py` proves the transversal statement in exact
rational arithmetic by subdivision of pose space (`zeromargin.py` conventions: `u = tan(θ/2)`,
exact bin cosines, clip to `[w_lo/2, t-w_lo/2]²`, six (axis, side) separations refuted by interval
arithmetic).  `t = 399/100`, `p = (99/100, 2)`, `ε = 1/250`, `ρ = 1.05 ρ*`:

    201,068 boxes, depth 30, EMPTY 60,412, MEET 46,522, uncertified 0    (44 s)
    VERDICT: CERTIFIED

It **terminates for every `ρ > ρ*`** and **does not terminate at `ρ = ρ*`**: the critical pose sits
on the admissibility boundary `c_x = w(θ)/2`, a curve that every box straddles — the exact
analogue of `ZEROMARGIN.md §4.1`, quadratic margin against linear erosion.  Lemma 2 is itself the
closed-form primitive for that case.  The required depth scales like the reciprocal of the
relative slack: at `ρ/ρ* = 1.05` it finishes at depth 30, while at `ρ/ρ* = 1.0003`
(`ρ = 351/12500`) it is still reporting undecided boxes at depth 22, and they cluster at
`θ ≈ 81.9° = 90° − arccos p_x` with `c_x ≈ w(θ)/2` — the mirror image of Lemma 2's witness pose,
exactly where the theory says the margin vanishes.  In practice `ρ = 1.05 ρ*` costs 5 % of the
anchor length and buys a subdivision proof.

## 2. Calibration and values at `t = 3.99`

The LP is `max Σ mu_S` over poses on a lattice, subject to coverage `<= 1` at every arrangement
vertex of the current support (`packing_dual.py`'s exact vertex kernel) and to the clique cuts of
the mode.  Poses are D4 orbits, rows live in the fundamental domain, and column generation is over
the lattice: `h > 0` prices the whole centre lattice at every angle with the C kernel; `h = 0` uses
`packing_dual.py`'s **arrangement sweep pricer**, which returns the exact minimum captured weight
over *all* centres at each angle — so `h = 0` is the `h → 0` limit of the lattice and its pricing
gap is a real optimality certificate for the LP over `{all admissible centres} × {angle grid}`.

**Calibration rule.**  Target `L(3.99) = 12.008` (`DUAL_EXACT.md`, exact).

| lattice `(h, dθ)` | pure value (stage sequence) | calibration |
|---|---|---|
| `h = 0` (exact centres), `dθ = 1°` | `12.0083, 12.0095, 12.0166, 12.0227, 12.0282` — passes at **round 0** and keeps rising | **PASS** |
| `h = 0.01`, `dθ = 1°` | `11.878, 11.914, 11.956, 11.981, 12.006` (rising, run stopped there) | **marginal** (`0.0025` low after 5 stages, still climbing) |
| `h = 0.02`, `dθ = 2°` | `11.746, 11.913, 11.966`; a restart from that support gave `11.909, 11.963` | **FAIL** (`0.045` low) |

So the clique numbers below were computed **only** on the lattice that passes calibration, which
is the point of the rule.  It is worth saying why the coarse lattices are not merely weaker but
actively misleading in this context: on a centre lattice a great many pose pairs *touch exactly*,
so the closed-intersection graph has far more edges than the continuum's and its cliques are
correspondingly larger — a coarse lattice inflates the clique numbers and deflates the pure one at
the same time.

**Values on the calibrated lattice `(h = 0, dθ = 1°)`.**  A stage value is the LP mass when the
inner loop (rows + cuts on a fixed pose set) stopped; "conv" means no violated row and no violated
cut of that family were found.

| mode | family | stage values | residual max clique of the support graph |
|---|---|---|---|
| `pure` | — | `12.0083 → 12.0095 → 12.0166 → 12.0227 → 12.0282` (rising) | `1.21–1.33` (the pure measures are not clique-feasible — task A) |
| `anchor` | `K(p, A)`, geometric membership | see the matched table below | `1.16–1.37` (unrestricted cliques are still violated; anchor cliques are not) |
| `sat` | arbitrary finite pose sets (task A's family) | inner iterations `11.975 → 11.932 → 11.911 → 11.902` on a fixed 514-pose set (834 cuts, still falling when the run was stopped) | `1.07` and falling |

**Matched pairs — the number that matters.**  At every stage the run also solves the *same* LP on
the *same* pose set with the cuts removed, so the clique gain is read off one pose set and not
across lattices (`runs/cq_A99sw5.log`):

Nine of the eleven inner loops ended with `+rows 0 +cuts 0` and residual max **anchor** clique
`1.0000` — that is, with the separator unable to find any violated anchor clique at all.

| stage | poses (orbits) | anchor cuts | pure value | anchor-clique value | gain | residual max anchor clique |
|---|---|---|---|---|---|---|
| `r0` | 558 | 422 | `12.011523` | `11.896525` | `0.115` | `≈ 1.00` (converged) |
| `r1` | 965 | 458 | `12.011310` | `11.902459` | `0.109` | `≈ 1.00` (converged) |
| `r2` | 1340 | 490 | `12.011327` | `11.903163` | `0.108` | `≈ 1.00` (converged) |
| `r3` | 1787 | 545 | `12.011641` | `11.907901` | `0.104` | `≈ 1.00` (converged) |
| `r4` | 2204 | 709 | `12.015987` | `11.895710` | `0.120` | `≈ 1.00`; row loop hit its cap (`maxcov 1.0013`), so this stage is a little optimistic |
| `r5` | 2649 | 784 | `12.018611` | `11.901056` | `0.118` | `≈ 1.00` (converged) |
| `r6` | 3109 | 816 | `12.020539` | `11.903942` | `0.117` | `≈ 1.00` (converged) |
| `r7` | 3515 | 850 | `12.023412` | `11.909723` | `0.114` | `≈ 1.00` (converged) |
| `r8` | 3916 | 882 | `12.024525` | `11.909865` | `0.115` | `≈ 1.00` (converged) |
| `r9` | 4316 | 913 | `12.025325` | `11.912448` | `0.113` | `≈ 1.00` (converged) |
| `r10` | 4739 | 932 | `12.025821` | `11.912723` | `0.113` | `≈ 1.00` (converged) |

Every pure value clears the calibration target `12.008`, and the certifiable clique family costs
the LP `0.104–0.120` — **thirteen to fifteen times** the pure method's entire excess over 12 at
`t = 3.99`.  **The ladder is flat.**  Over these eleven stages the pose set grows by a factor of
**8.5** (558 → 4739 orbits) and the clique value stays pinned at `11.905 ± 0.008`, while the pure
value on the very same poses climbs from `12.0115` to `12.0258`.  The gap does not close; it
widens slightly.  That is the qualitative difference from task A, whose stage values drifted `11.19 → 11.36 →
11.45` over a comparable pose refinement: A's cuts were explicit lists of poses, so every new
column arrived outside every cut; these cuts are geometric regions, so a new column that lands
inside one is charged and the LP has to move it out of the region — which it evidently cannot do
cheaply.

The unrestricted (`sat`) family reaches `11.90` with 834 explicit cuts and is still descending, so
on a fixed pose set the certifiable family recovers essentially all of what the unrestricted one
achieves.  The difference between them is not the value: it is that the anchor cuts survive
pricing and the `sat` cuts do not.

**Caveat, stated once and meant throughout.**  Every stage value is an LP over a finite pose set,
so it is a **lower** bound on the value over the continuum of poses at those angles, and it rises
as poses are added.  This is exactly the ladder task A could not resolve.  What is different here
is that the ladder is now *meaningful*: with geometric cuts a newly priced pose that lands inside
`K(p, A)` is charged, so the LP cannot escape a cut by moving a pose an epsilon — it has to move
it out of the region, and the region has positive measure.

## 3. Exactly certified numbers

`search/clique_exact.py` snaps a float support to rational rotations `2 arctan(p/q)` and rational
centres clamped exactly into `[0, 399/100]²`, enumerates the arrangement vertices exactly
(`dual_exact.py` machinery), separates anchor cliques, re-solves, and rounds the masses down so
that the maximum coverage and every cut are `<= 1` in exact rationals.

### A better exact ceiling constant at `t = 3.99` (side result)

The pure run on the calibrated lattice produced a measure that certifies, in exact rationals,

    mass = 12028160771 / 1000000000 = 12.028160771
    M    = max coverage = 3999999987 / 4000000000 = 0.999999996750   (exactly)
    L    = mass / M = 48112643084 / 3999999987 = 12.028160810092  >  12

so **`ν_f(399/100) >= 12.0281608…`**, improving `DUAL_EXACT.md`'s `12.00823078`.  Support file
`runs/cqx_PURE99_support.txt` (240 poses with positive mass, 1920 images).  Everything
load-bearing is integer/`Fraction`: the poses are rational rotations `2 arctan(p/q)` with rational
centres clamped exactly into the closed container, the arrangement vertices are exact integer
triples, the coverage test is the integer one of `dual_exact.py`, and the masses were rounded
**down** to multiples of `10^-9`.  Floats only chose the masses.

**Independently re-checked** by the pre-existing checker, which re-reads the file and re-derives
everything with no shared code path with the search:

```
python3 search/dual_exact.py check runs/cqx_PURE99_support.txt
  EXACT: 255893 vertices checked; mass = 12028160771/1000000000;
  M = 3999999987/4000000000 (attained at ~(0.999939, 0.999939));
  L = 48112643084/3999999987 = 12.028160810092  (>= 12);  63.8 s
```

Consequence, by the duality of `search/CEILING.md`: no weighted unavoidable set of total weight
`< 12.0281` exists at any container side `>= 3.99` — a slightly larger margin than before on the
statement that pins the pure method's ceiling.

### An exactly certified anchor-clique-feasible measure

`clique_exact.py` was also run with the cuts on, end to end: snap the anchor run's support to
rationals, enumerate the arrangement vertices exactly, separate anchor cliques and re-solve until
none is violated, then round the masses down and re-check everything in integers, with the
**exact** rational membership test for `K(p, A)` (three-axis SAT in `Fraction`s) rather than the
LP's float one.  Result (`runs/cqx_A99.log`, `runs/cqx_A99_support.txt`):

    mass    = 595075949/50000000  = 11.901518980
    M       = 3999999999/4000000000 = 0.999999999750   <= 1
    37 anchor cuts, exact maximum value 0.999999998250  <= 1

So the float LP value `11.9097` is reproduced exactly at `11.9015` on the snapped poses, with every
cut verified in exact arithmetic — the pipeline is sound end to end.  As a *bound* this is a lower
bound on the clique-LP value and therefore decides nothing in the useful direction; it is here as
validation.

| other exact numbers | value | note |
|---|---|---|
| pure LP over the anchor solution's own 77 support poses | `11971698781/1000000000 = 11.971698781`, `M = 3999999999/4000000000` | certified; those poses alone cannot reach 12 |

## 4. What the measures say, directly

Independently of any LP, on the **certified** extremal measure at `t = 3.99`
(`runs/dual_exact_3.99_support.txt`, mass `12.008230754`, coverage `<= 1` exactly):

| object | max mass | status |
|---|---|---|
| point clique (= coverage) | `1.000000` | certified exact |
| max-mass clique of the support's closed-intersection graph | `1.2139` (60–120 s B&B cap; task A: `1.327`) | float, valid, a lower bound on the true maximum |
| best **anchor clique** `K(p, A)` found (greedy over support squares **and** the targeted `(d, ε, ρ)` scan of Lemma 2) | `1.0571` | its `A` is a whole support square, so `{S : A ⊆ S}` is the **single pose** `A` — valid, but a measure-zero cut a cover cannot use |
| best anchor clique whose `A` is a genuine segment, so `{S : A ⊆ S}` is a pose **region** | `1.0242` at `p = (0.55, 2.99)` (`A` a segment of length `0.735`, 8 support images contain it); then `1.0193, 1.0157, 1.0083, …` | **the usable number** |
| task A's max clique, regularised: `p` its heaviest common point, `A` the intersection of the members missing `p` | `0.9695` | the rule "contains `p` **and** meets `A`" keeps only 54 of the 137 support poses through `p` |

**The mechanism, and it is a factor of ten.**  The certified measure violates *some* valid clique
by `0.21–0.33`; it violates the best anchor clique by `0.057`, but that one is carried by a single
pose; and it violates the best *region* anchor clique by `0.024`.  The region cliques sit at wall
distances `0.55, 0.75, 0.93` — exactly the band Lemma 1 allows — and the violation shrinks as the
distance approaches 1.  This is `CLIQUE_CEILING.md`'s "the excess sits on grazing contacts" with
the mechanism named: the grazing members are the ones a describable family cannot hold on to, and
what survives is an order of magnitude smaller — but still an order of magnitude larger than the
pure method's own room.

## 5. How big the extra family can ever be

`clique_family.py kset` computes, for a point `p`, the maximal point-anchored family
`K(p) = {S : S meets every admissible S' ∋ p}` as an exact convex region per angle (an *outer*
approximation, so an over-estimate of the extra), and its volume in pose space `(c_x, c_y, θ)`:

| `p` (wall distance) | vol `P_p` | vol `K(p) \ P_p` | ratio |
|---|---|---|---|
| `>= 1.000` from all walls | — | **0** | 0 |
| `0.999` | 1.38489 | `5.7e-5` | `4.1e-5` |
| `0.99` | 1.37215 | `1.77e-3` | `1.3e-3` |
| `0.95` | 1.31224 | `1.81e-2` | `1.4e-2` |
| `0.90` | 1.23236 | `4.83e-2` | `3.9e-2` |
| `0.80` | 1.06283 | `1.27e-1` | `0.120` |
| `0.60` | 0.71219 | `3.49e-1` | `0.491` |
| `0.50` | 0.54278 | `5.01e-1` | `0.922` |
| corner `(0.99, 0.99)` | 1.18167 | `2.15e-2` | `1.8e-2` |
| corner `(0.90, 0.90)` | 0.95470 | `2.38e-1` | `0.249` |

A single certifiable segment anchor recovers 46–62 % of that maximum along a wall and 5–10 % at a
corner (`clique_family.py anchorvol`; at a corner the right anchor is not wall-perpendicular, and
a two-anchor family would do better — designed, not built).  `K(p)` is itself pairwise
closed-intersecting on 15,930 sampled pairs over 9 points, i.e. it appears to *be* the maximal
clique containing `P_p` (measured, not proved).

## 6. What would change in `verify/`, `xcheck.py` and Lean

Design; `notes/clique-family.md §5` has the detail.  The `cliques N Q c` block keeps its header and
weight line; the box lines become *piece* lines over point and segment anchors:

    piece  a  f_1 … f_r        { S : A_a ⊆ S and S ∩ A_{f_i} ≠ ∅ for all i }
    anchorP  X Y D
    anchorS  X0 Y0 X1 Y1 D

Well-formedness the verifier must refuse otherwise: for every ordered pair of pieces `(i, j)`,
either `A_i ∩ A_j ≠ ∅` (exact rational test) or `j ∈ filt(i)` or `i ∈ filt(j)` — the analogue of
today's "cores pairwise meet".  Two exact `i128` predicates per swept cell: *every pose of the
cell contains the anchor* (the existing test, but against the **exact bin core** of
`ZEROMARGIN.md §2` rather than the `σ`-square, which is what makes a point anchor exact at edge
midpoints) and *every pose of the cell meets the anchor* (refute six (axis, side) separations by
interval arithmetic).  A cell that only partially satisfies a predicate gets no credit and its LP
witness goes outside the region — the same discipline as box cliques.  `xcheck.py` mirrors both in
`Fraction`s with anchors inflated rather than eroded.

Lean: `packing_le_weight_cliques` is unchanged (it already takes an arbitrary pose predicate with
a pairwise-intersection hypothesis).  One new lemma, a sibling of `clique_of_cores`, discharges it:

```lean
lemma clique_of_anchors {ι} (P : ι → ℝ × ℝ → ℝ → Prop) (anc : ι → Set (ℝ × ℝ))
    (hin   : ∀ i c θ, P i c θ → anc i ⊆ sq c θ 1)
    (hmeet : ∀ i j c θ, P i c θ → Meets i j → (sq c θ 1 ∩ anc j).Nonempty)
    (hpair : ∀ i j, (anc i ∩ anc j).Nonempty ∨ Meets i j ∨ Meets j i) :
    ∀ c c' θ θ', (∃ i, P i c θ) → (∃ j, P j c' θ') → (sq c θ 1 ∩ sq c' θ' 1).Nonempty
```

a three-way case split, no shrink lemma.  Lemma 2 is *not* needed for soundness — it says the
clique constraint dominates the point row, a statement about strength, not validity — so it need
not be formalised; the exact certifier of §1 is the practical substitute.

## 7. A near miss, recorded

The first version of the search dropped one separating axis.  A square and a *segment* need three
SAT axes — the square's two edge normals **and the segment's own normal** — and only the first two
were tested.  With that bug the scan reported an anchor clique of mass **1.3308** on the certified
measure, matching task A's `1.327` and looking like a breakthrough.  It is refuted by the corner
square `[0,1]²`, which contains the anchor point `p = (0.895, 0.995)` and misses the anchor segment
entirely (the endpoints are `(1.0974, 0.9946)` and `(0.8946, 1.1974)`, and at `x = 1` the segment
is at `y = 1.0920 > 1`).  The missing axis is exactly the constraint `ε < 1 - p_x` of Lemma 2.
Two independent implementations of the transversal test and the exact certifier now agree on `ρ*`
to 12 digits.  Recorded also in `notes/clique-family.md §7`.

## 7a. `t = 3.98`, corner leaf `k = 4`

`clique_continuum.py --kmass 4 --r 1` adds the leaf's equality row (total mass of poses centred in
the four corner boxes `[0,1]²` equals 4) with its multiplier in the pricing, the packing-side dual
of the branch certificate of `BRANCH.md` (`n <= W - λk`).  Warm starts: the branch dual
`runs/branch_t398hk4_dual_it16.txt` (520 poses), `dual_PD1_support.txt`, and — rescaled by
`3.98/3.99` — the new certified 3.99 measure.

### A pricing bug, found and fixed

The first leaf runs stalled at `11.916`.  Two faults, both in the leaf multiplier:

* the reduced cost was computed as `1 - capture + λ` for corner poses.  Dual feasibility for
  `max Σμ  s.t.  Aμ <= 1, Fμ = k` is `Aᵀy + Fᵀλ >= c`, so it is `1 - capture - λ`;
* worse, on the exact-centre lattice (`h = 0`) the pricer is `packing_dual.py`'s arrangement
  sweep, which knows nothing about `λ` and filters candidates on `capture < 1`.  With `λ < 0`
  the improving corner columns are exactly those with `capture ∈ [1, 1 - λ)`, so **every one of
  them was invisible**.

Measured `λ = -0.917` at the first stage, and the new explicit corner pricer (`price_corner`,
which rasters the origin corner box — the other three are its D4 images) found **574 improving
corner poses at once** that the sweep could not see.  The fix moved the leaf pure value from
`11.916` to `11.949` and drove `|λ|` from `0.92` down to `0.29`.  The `t = 3.99` runs are
unaffected: they have no equality row, and their row counts stayed below the `maxatoms` pricing
cap that also had to be raised here (16k rows at 3.98 against 5.9k at 3.99).

### Matched pairs (`runs/cq_LA98.log`), same protocol as §2

| stage | poses | anchor cuts | pure value on the same poses | anchor-clique value | gain | corner mass | `λ` |
|---|---|---|---|---|---|---|---|
| `r0` | 1707 | 302 | `11.913024` | `11.815254` | `0.098` | `4.000000` | `-0.952` |
| `r1` | 2126 | 438 | `11.913105` | `11.813382` | `0.100` | `4.000000` | `-1.044` |
| `r2` | 2944 | 573 | `11.913980` | `11.812264` | `0.102` | `4.000000` | `-0.991` |
| `r3` | 3367 | 601 | `11.916297` | `11.812316` | `0.104` | `4.000000` | `-1.123` |

The gain is `0.098–0.104` and creeping **up** as the pose set doubles — the same size and the same
flatness as the `0.104–0.120` measured at `t = 3.99`.  The residual max clique of the support
graph is `1.11–1.33`, i.e. unrestricted cliques are still violated while anchor cliques are not.

### Calibration: NOT reached — `0.051` short, and why

The best pure leaf value obtained is **`11.9487`** (`runs/cq_LP98b.log` r4, 4880 poses,
`λ = -0.289`, still rising), against the target `12.000 ± 0.001`.  So the `≈ 0.10` clique gain
**cannot yet be quoted against 12**: on the pose sets reachable here the pure reference is
`11.91–11.95`, not `12.00`.

The shortfall is a column-generation deficit at `t = 3.98` in general, **not** something specific
to the leaf.  The control that shows this: the *unconstrained* pure LP at `t = 3.98`, same
machinery, same time, reaches only `11.934` (`runs/cq_P98u.log`) — *below* the leaf LP's
`11.949`, even though the leaf LP carries an extra equality constraint and so has the smaller
feasible set.  Both are simply short of their common ceiling.  The reason is the warm start: at
`3.99` the run begins from the certified `12.008` measure, which is already essentially optimal,
so calibration passes at round 0; at `3.98` no such measure exists.  The best available seed,
`dual_PD1`, is itself a non-converged `11.918` (`DUAL.md`), and the branch dual is a cover-side
object.  Rescaling the new certified `3.99` measure down to `3.98` was tried (`--warm-from 3.99`,
`runs/cq_LP98c.log`) and starts a different but not better path (`11.770 → 11.885 → 11.918` and
climbing).

One further caveat about the target itself: `12.000 ± 0.001` is **uncertified**.  `CLIQUE.md`
records that the `t398hk4` run reached it at a degenerate vertex with a 2 % pricing gain, with the
value "straddled by the row and column steps" (`11.9995` on the column step, `12.000024` on the
row step).  It is a soft target, and the true leaf value could be a little below it.

**Leaf verdict: undecided, because the pure reference is not yet converged.**  What is
established: the leaf's anchor-clique gain is `0.098–0.104`, stable, the same as at `3.99`.  What
is not: whether the leaf's pure value is really `12.000`.  If it is, the same gain would put the
anchor-clique leaf value at `≈ 11.90` — the leaf would close with room and `s(12) >= 3.98` would
follow — but that is an extrapolation from a `0.05` extrapolated reference, not a measurement, and
this write-up does not claim it.  What would settle it is not more clique work but a converged
pure LP at `3.98`: an exactly certified `t = 3.98` measure of the kind §3 produced at `3.99`,
which is a self-contained follow-up worth about a day of compute.

No exactly certified anchor-clique-feasible leaf measure is reported: the rule stated in §3 is
that such a measure is a lower bound on the clique-LP value and decides nothing in the useful
direction, and here it would additionally be quoted against an unconverged reference.

## 8. What was not done

* **A converged pure LP at `t = 3.98`.**  Without it the leaf's clique gain (measured, `0.10`) cannot be quoted against 12 — see §7a.  This, not more clique machinery, is the blocker on the `k = 4` leaf.
* **The cover side.**  Nothing here is a certificate.  The decisive experiment — a clique
  certificate of weight `< 12` — is on the cover side and needs the verifier changes of §6.
* **Multi-anchor cliques (`m ≥ 3`).**  Implemented (`clique_family.kmass_multi`,
  `anchor_multi`) and correct, but only lightly searched; at a corner they are the right shape and
  a single perpendicular segment recovers only 5–10 % of the available volume there.

## Reproduce

```sh
M=/home/evand/math/square-packing/s12/runs          # the main tree's runs/, for warm starts
python3 search/clique_family.py rho                                  # the rho* law
python3 search/clique_family.py kset --nth 180 --nthp 90             # vol K(p) \ P_p
python3 search/clique_family.py anchorvol --nth 120 --nthp 90        # what one segment recovers
python3 search/clique_family.py kclique --nth 120 --n 60             # is K(p) a clique?
python3 search/clique_family.py separate --pitch 0.02 --thr 0.97 --max-pts 30   # separation
python3 search/clique_family.py regular --tl 120                     # regularise A's max clique
python3 search/clique_cert.py cert --px 99/100 --py 2 --eps 1/250 --depth 30    # exact, 44 s
python3 search/clique_cert.py stress --n 300000
python3 search/clique_continuum.py 3.99 P99sw1 --mode pure   --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_continuum.py 3.99 A99sw5 --mode anchor --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_continuum.py 3.99 S99sw1 --mode sat    --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_continuum.py 3.99 P99h01 --mode pure   --h 0.01 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_exact.py runs/cq_P99sw1_support.txt PURE99 --no-cuts
python3 search/clique_exact.py runs/cq_A99sw5_support.txt A99 --rounds 8
```
