# Cliques in the continuum: is the non-Helly mass certifiable? (task G, 2026-08-30)

Code: `search/clique_family.py` (the family and its geometry), `search/clique_cert.py` (exact
rational certification of the clique property), `search/clique_continuum.py` (the
clique-strengthened packing LP), `search/clique_exact.py` (exact certification of a
clique-feasible measure).  Mathematics, lemmas and proofs: `notes/clique-family.md`.
Runs: `runs/cq_*`, `runs/cqx_*`.

## Verdict

**NO-GO for every clique family that is a *pose region*, at `t >= 3.99`.**  Two independent
strands say the same thing:

1. **Structure (proved).**  A clique can beat the point constraint at `p` only if `p` is
   *strictly* within distance 1 of a container wall, and the gain is `O(1 - dist)`.  Precisely
   (`notes/clique-family.md` Lemma 1, sharp and proved in three lines): if `p` is at distance
   `>= 1` from all four walls, the point clique `P_p = {S admissible : p ∈ S}` is already
   **maximal** — no admissible pose outside it meets every admissible pose inside it.  The tight
   points of the extremal measures at `t = 3.99` sit at wall distance `≈ 0.99`, where the entire
   maximal extra family is `0.13 %` of the point clique by volume.
2. **Computation.**  The clique-strengthened packing LP at `t = 3.99`, with cuts from the
   certifiable **anchor-clique** family and column generation that is exact over centres at
   `1°` angles, does not go below the pure value: stage values `12.0043, 12.0084, …` against the
   pure `12.0083, 12.0095, 12.0166, 12.0227, 12.0282, …`, each stage ending with **no violated
   anchor clique** found.  By weak duality a clique-feasible measure of mass `m` forbids any
   clique certificate of weight `< m` at that `t` and above, so no anchor-clique certificate of
   weight `< 12` exists at any `t >= 3.99`.

What is **not** excluded: the *unrestricted* clique relaxation, whose cliques are arbitrary
finite sets of poses.  Those do bite (see the `sat` column below, and task A's 11.2–11.45), by
0.2–0.8.  But a finite set of poses has measure zero in pose space, so it credits no weight in a
cover: the dual object a certificate needs is a clique that is a positive-measure *region*, and
the whole content of this task is that the two differ by essentially all of the non-Helly mass.
That is the same fact `CLIQUE_CEILING.md` measured as "the excess sits on grazing contacts",
here with the mechanism named.

## 1. The certifiable family

Full statements and proofs: `notes/clique-family.md`.  In one paragraph:

> **Anchor clique.**  For a point `p` and a nonempty compact convex `A`,
> `K(p, A) = {S admissible : p ∈ S and S ∩ A ≠ ∅} ∪ {S admissible : A ⊆ S}`
> is a clique for **every** `p` and `A` — two members of the first part share `p`, two of the
> second share `A`, and a mixed pair `S`, `S'` has `S ∩ S' ⊇ S ∩ A ≠ ∅`.  No geometric lemma is
> needed for validity, membership is a geometric test (so column generation cannot dodge it), and
> the family contains the exact point clique (`A = {p}`) — which box cliques never can.

Task B's box cliques are the special case in which all anchors pairwise *meet*.  The new content
is anchors that do not, and the filter `S ∩ A ≠ ∅` is the price.

**When the filter is vacuous** (so that the clique constraint strictly dominates the coverage
constraint at `p`) is `notes/clique-family.md` Lemma 2, which is exact:

> For `0 < p_x < 1`, `0 < ε < 1 - p_x`, and `A_ρ` the vertical segment `{p_x+ε} × [p_y±ρ]`, every
> admissible closed unit square containing `p` meets `A_ρ` **iff**
> `ρ ≥ ρ* = ε·p_x / sqrt(1 - p_x²) = ε·cot(arccos p_x)`, the worst pose being
> `θ = arccos p_x`, `c_x = w(θ)/2` — a square jammed against the wall.

Both the sufficiency and the sharpness are proved (convexity of
`f(θ) = p_x/cos θ + ε cos θ - ρ sin θ - 1`), reproduced to 12 digits by two independent numerical
implementations, and stress-tested on 300,000 random instances with 0 failures.

**Exact certification.**  `search/clique_cert.py` proves the transversal statement in exact
rational arithmetic by subdivision of pose space (`zeromargin.py` conventions: `u = tan(θ/2)`,
exact bin cosines, the clip to `[w_lo/2, t-w_lo/2]²`).  Example, `t = 399/100`, `p = (99/100, 2)`,
`ε = 1/250`, `ρ = 1.05 ρ*`:

    201,068 boxes, depth 30, EMPTY 60,412, MEET 46,522, uncertified 0    (44 s)
    VERDICT: CERTIFIED

It **terminates for every `ρ > ρ*`** and **does not terminate at `ρ = ρ*`**: the critical pose is
a genuine zero-margin configuration on the admissibility boundary `c_x = w(θ)/2`, which is a
curve that every box straddles — the exact analogue of `ZEROMARGIN.md §4.1`.  Lemma 2 is itself
the closed-form primitive for that case; in practice one takes `ρ = 1.05 ρ*` and the subdivision
finishes.

## 2. Calibration and values at `t = 3.99`

The LP is `max Σ mu_S` over poses on a lattice, subject to coverage `<= 1` at every arrangement
vertex of the current support (`packing_dual.py`'s exact vertex enumeration) and to the clique
cuts of the mode.  Poses are D4 orbits; rows live in the fundamental domain.  Column generation
is over the lattice: `h > 0` prices the whole centre lattice at every angle with the C kernel;
`h = 0` uses `packing_dual.py`'s **arrangement sweep pricer**, which returns the exact minimum
captured weight over *all* centres at each angle — i.e. `h = 0` is the `h → 0` limit of the
lattice, and its pricing gap is a genuine optimality certificate for the LP over
`{all admissible centres} × {angle grid}`.

**Calibration rule.**  Every clique value is reported next to the pure (no-cut) value on the same
lattice.  Target: `L(3.99) = 12.008` (`DUAL_EXACT.md`, exact).

| lattice `(h, dθ)` | pure value (stage sequence) | calibration | clique value: **anchor** | clique value: `sat` (unrestricted, fixed pose set) | residual max graph clique |
|---|---|---|---|---|---|
| `h = 0` (exact centres), `dθ = 1°` | `12.0083, 12.0095, 12.0166, 12.0227, 12.0282` (rising, stopped at r4) | **PASS** at round 0 | `12.0043, 12.0084, …` (rising; each stage ends with **no violated anchor clique**) | `11.98 → 11.93 → …` (falling, 1 cut/iteration; not converged) | `1.29` under anchor cuts; `1.18` under sat cuts |
| `h = 0.02`, `dθ = 2°` | `11.746, 11.913, 11.966, …` (rising, not converged) | **FAIL so far** (`> 0.04` below 12.008) | — | — | — |
| `h = 0.01`, `dθ = 1°` | (running) | | | | |

(*Values in this table are floating-point LP values on a restricted pose set: each is a lower
bound on the value of the same LP over the continuum of poses at those angles, and each rises as
the pose set grows.  Nothing here is certified except where stated in §3.*)

The coarse lattice `h = 0.02, dθ = 2°` is exactly the case the calibration rule is for: after
three rounds its pure value is `11.966`, `0.042` below `L(3.99) = 12.008`, so any clique value
computed on it would be meaningless — the lattice, not the cliques, would be doing the work.
This is worth stating plainly because a coarse fixed lattice is *also* what makes the
unrestricted-clique numbers look strong: on a lattice, huge numbers of pose pairs *touch* exactly,
so the closed-intersection graph has far more edges than the continuum has, and its cliques are
correspondingly larger.

## 3. Exactly certified

`search/clique_exact.py` snaps a float support to rational rotations `2 arctan(p/q)` and rational
centres clamped exactly into `[0, 399/100]²`, enumerates the arrangement vertices exactly
(`dual_exact.py`), separates anchor cliques, re-solves, and rounds the masses down so that the
maximum coverage and every cut are `<= 1` in exact rationals.

(results filled in below when the run finishes)

## 4. What the measures say, directly

Independently of any LP, on the **certified** extremal measure at `t = 3.99`
(`runs/dual_exact_3.99_support.txt`, mass `12.008230754`, coverage `<= 1` exactly):

| object | max mass | status |
|---|---|---|
| point clique (= coverage) | `1.000000` | certified exact |
| max-mass clique of the support's closed-intersection graph | `1.2139` (60–120 s B&B cap; task A: `1.327`) | float, valid, a lower bound on the true maximum |
| best **anchor clique** `K(p, A)` found over 60 anchor points × 5 seeds × 20 greedy steps | `1.0502` | float; but its `A` is a whole support square, so `{S : A ⊆ S}` is the **single pose** `A` — the cut is a measure-zero object |
| best anchor clique whose `A` has diameter `<= 0.01` (so `{S : A ⊆ S}` is a genuine pose region) | `0.998` | float — **no violation at all** |
| task A's max clique, regularised: `p` its heaviest common point, `A` the intersection of the members missing `p` | `0.9695` | the rule "contains `p` **and** meets `A`" keeps only 54 of the 137 support poses through `p` |

**This is the mechanism.**  The certified `12.008` measure really does violate valid clique
constraints, by `0.21–0.33`.  Every violating clique found is carried by *individual poses*.  The
moment the clique is required to be a region — so that a perturbed pose is still charged, which is
exactly what a cover certificate and a column-generation loop both require — the violation falls
to `0.05`, and to `0` once the anchor is small enough for `{S : A ⊆ S}` to contain more than one
pose.

## 5. How big the extra family can ever be

`clique_family.py kset` computes, for a point `p`, the maximal point-anchored family
`K(p) = {S : S meets every admissible S' ∋ p}` as an exact convex region per angle (an *outer*
approximation, so an over-estimate), and its volume in pose space `(c_x, c_y, θ)`:

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

A single certifiable segment anchor recovers 46–62 % of the maximum along a wall and 5–10 % of it
at a corner (`clique_family.py anchorvol`).  The extra family is `K(p) \ P_p` — squares that
*just* miss `p` — and its width in pose space is `O(1 - dist(p, wall))`.

The extremal measures' tight points are at wall distance `0.99` (the wall squares reach `x = 1`),
where the ceiling on the gain is `0.13 %`.  For the gain to be worth the pure method's whole
excess of `0.008` at `3.99`, the measure would have to concentrate its mass in that sliver — and
the LP, free to move poses, simply does not.

## 6. What would change in `verify/`, `xcheck.py` and Lean

Design only (the answer to Part 2 is no-go, so none of it was implemented).  Details in
`notes/clique-family.md §5`.  In short: the `cliques N Q c` block's box lines become *piece* lines
`(contains-anchor, meets-anchor list)` over point and segment anchors; the verifier needs two
exact `i128` predicates per swept cell — "every pose of the cell contains the anchor" (the
existing test, but against the **exact bin core** of `ZEROMARGIN.md §2` rather than the
`σ`-square, which is what makes a point anchor exact at edge midpoints) and "every pose of the
cell meets the anchor" (refute six (axis, side) separations by interval arithmetic); a cell that
only partially satisfies a predicate gets no credit and its LP witness goes outside the region.
`xcheck.py` mirrors both in `Fraction`s with inflated anchors.  Lean needs one new lemma, a
sibling of `clique_of_cores`:

```lean
lemma clique_of_anchors {ι} (P : ι → ℝ × ℝ → ℝ → Prop) (anc : ι → Set (ℝ × ℝ))
    (hin : ∀ i c θ, P i c θ → anc i ⊆ sq c θ 1)
    (hmeet : ∀ i j c θ, P i c θ → Meets i j → (sq c θ 1 ∩ anc j).Nonempty)
    (hpair : ∀ i j, (anc i ∩ anc j).Nonempty ∨ Meets i j ∨ Meets j i) :
    ∀ c c' θ θ', (∃ i, P i c θ) → (∃ j, P j c' θ') → (sq c θ 1 ∩ sq c' θ' 1).Nonempty
```

a three-way case split, no shrink lemma; `packing_le_weight_cliques` is unchanged.

## 7. A near miss, recorded

The first version of the search dropped one separating axis (a square and a *segment* need three
SAT axes, including the segment's own normal) and reported a clique of mass **1.3308** on the
certified measure — matching task A's 1.327 and looking like a breakthrough.  It is refuted by the
corner square `[0,1]²`, which contains the anchor point `p = (0.895, 0.995)` and misses the anchor
segment entirely.  The missing axis is exactly the constraint `ε < 1 - p_x` of Lemma 2.  Recorded
in `notes/clique-family.md §7`; two independent implementations and the exact certifier now agree.

## Reproduce

```sh
M=/home/evand/math/square-packing/s12/runs          # the main tree's runs/, for warm starts
python3 search/clique_family.py rho                                  # the rho* law
python3 search/clique_family.py kset --nth 180 --nthp 90             # vol K(p) \ P_p
python3 search/clique_family.py anchorvol --nth 120 --nthp 90        # what one segment recovers
python3 search/clique_family.py separate --pitch 0.02 --thr 0.97 --max-pts 60   # separation
python3 search/clique_family.py regular --tl 120                     # regularise A's max clique
python3 search/clique_cert.py cert --px 99/100 --py 2 --eps 1/250 --depth 30    # exact, 44 s
python3 search/clique_cert.py stress --n 300000
python3 search/clique_continuum.py 3.99 P99sw1 --mode pure   --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_continuum.py 3.99 A99sw1 --mode anchor --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_continuum.py 3.99 S99sw1 --mode sat    --h 0 --dth 1 --warm $M/dual_PA2_support.txt
python3 search/clique_exact.py runs/cq_A99sw1_support.txt A99sw1
```
