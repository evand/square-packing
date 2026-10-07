# Branch semantics at a region boundary

What a branch certificate means when a pose sits exactly on the boundary between two regions, which
of the possible meanings the shipped `verify/` + `lean/` pair implements, and what makes the leaf
tree exhaustive.

Written in answer to `search/T4SCREEN.md` §5 ("the leaf optimum sits on the slot boundary … no
cell-based verifier can implement a `1e-14` tie-break") and §7 item 1.
Companion: `notes/chord-lemma.md` (the lemma the enumeration needs), `notes/level2-design.md`
(the region design), `search/BRANCH.md` (level 1), `certificates/FORMAT.md` (the file format).

---

## 0. Verdict, up front

1. **The worry in T4SCREEN §5 is answered, and the answer is that no tie-break has to be
   implemented.**  The verifier's actual per-cell rule — *a cell straddling a boundary must meet
   every adjacent region's threshold* — is sound **simultaneously for every tie-break**.  A
   certificate that passes it refutes the leaf however a boundary pose is booked.  Its packing-side
   dual is that the fractional packing may split boundary mass between the adjacent regions in any
   proportion, which is exactly the freedom T4SCREEN's optima were exploiting.

2. **At level 1 (four corner boxes) the question is vacuous** and the shipped pair is sound as it
   stands: the corner boxes are pairwise disjoint, so a pose lies in at most one closed region,
   "closed" and "half-open" coincide, and `packing_le_weight_regions` covers the verifier exactly.

3. **At level 2 (thirteen regions sharing boundaries) `packing_le_weight_regions` does *not* cover
   the verifier's rule**: it demands the *sum* of the multipliers of all regions containing the pose,
   the verifier enforces the *max*.  The gap is closed by a new theorem,
   `packing_le_weight_regions_choice` (`lean/Sqpack/Basic.lean`, 0 sorries), stated and proved in §3.

4. **One latent soundness gap in `verify/`**, harmless for every shipped certificate but real: the
   region trailer checks `(2r-1)^2 < 2` but never `2r < s`, and with `2r >= s` the four corner boxes
   overlap and the per-box `required` rule is no longer the max.  §3.4 gives the one-line fix.

5. **Leaf counts** (`search/branch_leaves.py`, Burnside-checked), using only the chord lemma and the
   corner-box diameter bound, for the whole two-level tree up to `D4`:

   | corner leaf | `k` | slot leaves |
   |---|---|---|
   | `0000` | 0 | 1300 |
   | `0001` | 1 | 1830 |
   | `0011` | 2 | 546 |
   | `0101` | 2 | 351 |
   | `0111` | 3 | 171 |
   | `1111` | 4 | **15** |
   | | | **4213** |

   The `15` reproduces `search/level2_leaves.py`.  Capping a single slot at two centres (measured,
   not proved) brings the total to `2697`.

6. **The `m <= 1` leaves are not refutable by an interior capacity certificate at `t = 4`, and this
   is provable from work already in the repo.**  `search/T4SCREEN.md` §0 reports the interior LP
   (all frame regions pinned to zero) *converged* — coverage `M = 1`, clique max `1` — at
   `7.433677`, `7.508327`, `7.528038` and still climbing.  By LP duality a cover of the interior
   region of total weight `< 7` would force that value below `7`, so no such certificate exists at
   `t = 4`.  Separately, the *interior-disjoint* interior capacity at `t = 4` is at least **9**, by
   the explicit `3 x 3` unit grid; `search/level2_capacity.py --t 4.0 --box 1 3 1 3 --k 7` reports
   `-0.018354` and is therefore simply wrong on this instance (§6).  The shipped certificate format
   cannot state a region-restricted claim anyway; §6 says what the minimal extension is.

---

## 1. The reduction

Container `C = [0,s]^2`.  A **pose** is a pair `(c, θ)` whose closed unit square `S(c,θ)` lies in
`C`.  A **cover** is a finite weighted point set `(A, w)`, `w >= 0`, of total weight `W`; write
`cap(c,θ) = Σ_{a ∈ A ∩ S(c,θ)} w_a` for the weight a pose captures.

**Regions** `R_1, …, R_m` are arbitrary sets of poses, **multipliers** `λ_1, …, λ_m` arbitrary
reals (either sign).  For a packing — `n` squares of side `L > 1` with pairwise disjoint interiors
inside `C` — let `k_j` be the number whose pose lies in `R_j`.

> **Theorem (the reduction).**  If every pose satisfies
> `cap(c,θ) >= 1 + Σ_j [ (c,θ) ∈ R_j ] λ_j`, then
>
> ```
> n + Σ_j λ_j k_j  <=  W .
> ```
>
> Hence a certificate with `W - Σ_j λ_j k_j < n` refutes every packing with exactly that occupancy
> vector `k`.

This is `packing_le_weight_regions` in `lean/Sqpack/Basic.lean`, verbatim:

```lean
theorem packing_le_weight_regions
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (m : ℕ) (R : Fin m → ℝ × ℝ → ℝ → Prop) (lam : Fin m → ℝ)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C →
        1 + ∑ j : Fin m, (if R j c θ then lam j else 0)
          ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L)) :
    (n : ℝ) + ∑ j : Fin m, lam j * (((univ : Finset (Fin n)).filter
        (fun i => R j (ctr i) (ang i))).card : ℝ) ≤ ∑ a ∈ A, w a
```

The proof is one line of counting: the concentric closed *unit* squares of the packing are pairwise
disjoint (`unit_subset_interior`, needs `L > 1`), so no atom is captured twice, so
`Σ_i cap_i <= W`; and `Σ_i cap_i >= Σ_i (1 + Σ_j [pose_i ∈ R_j] λ_j) = n + Σ_j λ_j k_j`.
Nothing in it cares whether the regions overlap, partition, or are measurable.

---

## 2. Two ways to make the regions cover the container

The regions of interest are **centre boxes**: `R_j = {(c,θ) : c ∈ B_j}` for an axis-parallel box
`B_j`.  For the tree to be exhaustive the boxes must cover the container, and the natural boxes —
`notes/level2-design.md` §2.1's four corner boxes, eight wall slots and interior at `r = 1` — are
closed and share their boundary segments.  Two readings.

### (i) Closed regions, thresholds add

Take the `B_j` closed, exactly as drawn.  A centre on a shared boundary lies in **every** box
containing it, and the theorem's hypothesis at such a pose is

```
cap(c,θ)  >=  1 + Σ_{j : c ∈ B_j} λ_j                              (sum rule)
```

and the same pose is counted in **every** `k_j` for those `j`.  This is self-consistent and is
exactly what `packing_le_weight_regions` says; it needs no tie-break at all.  Its costs:

* the cover constraint at a boundary pose is the *strongest* of all — `1 + λ_i + λ_j` where the
  slots `i, j` meet, and `1 + λ_{C} + λ_{W} + λ_{I}` at a triple point.  With positive multipliers
  (the `k = 4` leaf wants `λ > 0`) this is a genuinely harder cover than either single constraint;
* the leaf indexing is by a vector `k` that may total more than `n`: the `3 x 3` grid of §6 has four
  centres on the boundary of the interior box and each is counted twice.  The enumeration has to
  allow for that, which roughly squares the leaf count near boundaries;
* the verifier does **not** implement it (§3).

### (ii) Half-open partition with a fixed tie-break

Replace `B_j` by `B_j^◦`, a half-open version — concretely, assign each centre to the region of
least index among the closed boxes containing it, which is what `search/level2_regions.classify`
does.  Then `R_1^◦, …, R_m^◦` partition the pose set, each pose lies in exactly one, the hypothesis
is

```
cap(c,θ)  >=  1 + λ_{j(c,θ)}                                       (partition rule)
```

and `Σ_j k_j = n` exactly.  Because `B_j^◦ ⊆ B_j`, every capacity bound proved for the closed box
(the corner-box diameter bound, the chord lemma) still holds for the half-open one, so the
enumeration of §4 is unaffected by the choice of tie-break.

**This is the design to use**, and §3 shows the verifier already implements it — for *every*
tie-break at once.

---

## 3. What a cell-based verifier can actually check, and why it is sound

### 3.1 The verifier's rule

`verify/`'s arrangement sweep works on **cells**: for each angle bin it partitions the admissible
centres into rectangles on which the captured weight is constant.  It cannot evaluate a pose; it
can only decide, conservatively, which boxes a *cell* meets.  The rule is (`verify/src/main.rs`,
the `required` closure in the bin sweep):

```rust
let required = |may: u8, inside: usize| -> i128 {
    if may == 0 { return 0; }                       // meets no box: threshold 1
    if inside < 4 { return lams[inside]; }          // wholly inside box `inside`: 1 + lam
    let mut req = 0;                                // straddling: the max, and 0 for "outside"
    for j in 0..4 { if may & (1 << j) != 0 && lams[j] > req { req = lams[j]; } }
    req
};
```

and the cell passes iff `weight >= W + required`.  `may` is a superset of the boxes the cell meets
and `inside` is `4` unless the cell is provably inside one box; the tests are on the cell's
bounding box with a padding that can only enlarge `may` and shrink `inside`
(`region_flags`, `RPAD = 1e-7`).  So, in words:

> **a cell that may meet several regions must satisfy every one of their thresholds** (and the
> threshold `1` of "outside every box", unless the cell is provably inside one).

No tie-break appears anywhere, and none is decided at `1e-14`.

### 3.2 Soundness — the lemma

> **Lemma (straddle soundness).**  Let `R_1, …, R_m` be arbitrary regions and `λ_1, …, λ_m` arbitrary
> reals.  Suppose the cover satisfies, for every pose,
>
> ```
> cap(c,θ) >= 1 + λ_j       for every j with (c,θ) ∈ R_j .              (max rule)
> ```
>
> Let `σ` assign to each square `i` of the packing some index `σ(i)` with `pose_i ∈ R_{σ(i)}`.  Then
>
> ```
> n + Σ_i λ_{σ(i)}  <=  W ,        i.e.       n + Σ_j λ_j · #{i : σ(i) = j}  <=  W .
> ```

*Proof.*  `Σ_i cap_i <= W` as in §1 (no atom is captured twice).  And `cap_i >= 1 + λ_{σ(i)}` by the
hypothesis applied at `pose_i` with `j = σ(i)`.  Sum. ∎

The hypothesis is exactly the verifier's rule: if a cell meets `R_j` then in particular every pose
of that cell that lies in `R_j` is required to reach `1 + λ_j`, and the verifier requires it of the
whole cell, which is stronger.  The "outside" region is the case `λ = 0`, and a cell provably inside
one box legitimately escapes it.

In Lean (`lean/Sqpack/Basic.lean`, new):

```lean
theorem packing_le_weight_regions_choice
    (A : Finset (ℝ × ℝ)) (w : ℝ × ℝ → ℝ) (hw : ∀ a ∈ A, 0 ≤ w a)
    (C : Set (ℝ × ℝ)) (m : ℕ) (R : Fin m → ℝ × ℝ → ℝ → Prop) (lam : Fin m → ℝ)
    (hcover : ∀ (c : ℝ × ℝ) (θ : ℝ), sq c θ 1 ⊆ C → ∀ j : Fin m, R j c θ →
        1 + lam j ≤ ∑ a ∈ A.filter (fun a => a ∈ sq c θ 1), w a)
    (n : ℕ) (L : ℝ) (hL : 1 < L) (ctr : Fin n → ℝ × ℝ) (ang : Fin n → ℝ)
    (hin : ∀ i, sq (ctr i) (ang i) L ⊆ C)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sqInt (ctr i) (ang i) L) (sqInt (ctr j) (ang j) L))
    (σ : Fin n → Fin m) (hσ : ∀ i, R (σ i) (ctr i) (ang i)) :
    (n : ℝ) + ∑ i : Fin n, lam (σ i) ≤ ∑ a ∈ A, w a
```

with `sum_assign_eq_sum_counts` rewriting `Σ_i lam (σ i)` as `Σ_j lam j · #{i : σ i = j}`.
Both are sorry-free; `#print axioms` gives `propext, Classical.choice, Quot.sound`.

**Consequence.**  A certificate passing the verifier refutes the leaf `k` under *every* tie-break
simultaneously: the conclusion is quantified over all assignments `σ`.  In particular it refutes
the leaf under the half-open partition of §2(ii), whatever the tie-break, and the verifier never
needs to know which one.  T4SCREEN §5's "printing the pose checkpoint to 13 decimals instead of 17
flipped the pose's region" is then a property of the *search*, not of the certificate: the
certificate is valid for both readings of that pose.

### 3.3 The packing-side dual

Dualise the leaf LP.  The max rule is the pair of linear constraints
`cap(p) >= 1 + λ_i`, `cap(p) >= 1 + λ_j` at a pose `p ∈ R_i ∩ R_j`, not one constraint with the
larger right-hand side; so in the dual the pose `p` appears as **two columns**, "p booked to `i`"
and "p booked to `j`", both with the same coverage footprint, and the region-mass equalities
`μ(R_j) = k_j` see them separately.  Hence:

> a boundary pose may be assigned to either adjacent region — at the packing's choice, and in the
> fractional relaxation with any split of its mass between them.

That is precisely what T4SCREEN §5 observed: `01011010` step 1 puts a whole unit of mass at
`(3.499547, 2.000000)`, on the `W_2 / W_3` boundary, and satisfies `μ(W_3) = 1, μ(W_2) = 0`.  Under
the correct semantics that measure is *feasible* for the leaf, and it is also feasible for the leaf
with the mass booked to `W_2` instead.  The consequence is not that the value is meaningless — it is
that the value is the max over the bookings, so the branch is weaker than a fixed tie-break would
suggest.  **The leaf values in T4SCREEN §3 are therefore values of an LP that a certificate *can*
dualise**, correcting §5 item 1 and §7 item 1; what §5 rightly shows is that the mid-wall line is
where the branch buys nothing, which is the argument for Design B's mid-wall gap on *value*
grounds, not on verifiability grounds.

### 3.4 What the shipped pair implements today, and the one gap

| | |
|---|---|
| `verify/` region trailer | four **closed** corner boxes `[0,r]^2` and images, membership by the centre; the `required` rule of §3.1 |
| `certificates/FORMAT.md` | "the four corner boxes … (closed)"; straddling cell "is required to reach both" |
| `xcheck.py` | the same classification in exact rationals, no padding |
| Lean | `packing_le_weight_regions` (sum rule); now also `packing_le_weight_regions_choice` (max rule) |

**Level 1 is sound and fully covered.**  For `2r < s` the four corner boxes are pairwise disjoint
closed sets, so a pose lies in at most one; the sum rule and the max rule coincide, semantics (i)
and (ii) coincide, and `packing_le_weight_regions` is an exact match for what `verify/` checks.
The shipped `certificates/branch/s12_t3.98_corner_k*.txt` have `r = 1`, `s = 3.98`, so `2r < s`
holds with room.

**The gap.**  `read_cert` checks `(2r-1)^2 < 2` (so each box holds at most one centre) but never
`2r < s`.  With `s <= 2r <= 1 + √2` the boxes overlap, and then

* the per-box form (`lambda L1 L2 L3 L4 / k K1 K2 K3 K4`) is **unsound as verified**: for a cell
  inside two boxes, `region_flags` sets `inside` to the larger index and `required` returns
  `lams[inside]`, ignoring the other — so a pose in both boxes is only required to reach
  `1 + λ_{inside}` where the max rule needs `1 + max(λ_i, λ_j)` (and the closed/sum reading needs
  `1 + λ_i + λ_j`);
* the single-`λ` form is unaffected: there the region is the *union* of the four boxes and overlap
  is invisible.

No shipped certificate is anywhere near this (`s ≈ 3.98`, `2r = 2`), and a container of side
`<= 1 + √2 = 2.414` cannot hold more than a handful of unit squares, so this is latent rather than
live.  Two ways to close it, both small and both outside the scope of this note (`verify/` is not
modified here):

* **(a), one line, exact integers.**  In the region trailer validation, additionally require
  `2 · r_num · s_den < s_num · r_den`.  Then the boxes are pairwise disjoint and everything above
  applies unchanged, with `packing_le_weight_regions` as the reduction.
* **(b), three lines, and strictly more general.**  Make `required` the max over *all* boxes the
  cell meets, i.e. delete the `if inside < 4` early return and instead compute
  `max { λ_j : j ∈ may }` when the cell is inside some box, `max({0} ∪ { λ_j : j ∈ may })` otherwise.
  Then the reduction is `packing_le_weight_regions_choice` and overlapping regions are legal.
  (b) is what a level-2 trailer with thirteen sharing boxes will need anyway.

**Level 2 is not implemented at all**: `region corner r_num r_den` is the only region kind the
parser accepts, and `notes/level2-design.md` §7's list-of-boxes trailer is unbuilt.  When it is
built, (b) is the rule to build, and `packing_le_weight_regions_choice` is the theorem it
discharges.

---

## 4. Exhaustiveness: what has to be enumerated

Fix the half-open partition of §2(ii) with any tie-break, so every packing has a well-defined
occupancy vector.  By §3 a certificate that passes the verifier refutes its leaf under *this*
tie-break among others, so:

> **The tree is exhaustive iff, for every packing of `12` unit squares in a container of side
> `s' < t`, the occupancy vector under the fixed tie-break is one of the enumerated leaves.**

Since the half-open regions are subsets of the closed ones, any *upper* bound on the count of
centres in a closed region is also a bound on the half-open count, so the enumeration may be done
entirely with closed-box capacity statements.  Three are used.

**(C1) A corner box holds at most one centre.**  The admissible part of `B = [0,r]^2` is
`[1/2, r]^2` (a unit square inside the container has `c_x, c_y >= 1/2`), of diameter
`(r - 1/2)√2`; two squares of a packing have centres more than `1` apart (a separating axis gives
`|⟨Δc, n⟩| > 1` for some edge normal).  So `(2r-1)^2 < 2` suffices — the exact check `verify/`
already performs.  At `r = 1`: `1 < 2` ✓.

**(C2) The chord lemma** (`notes/chord-lemma.md`, Lean `wall_strip_le_three_of_packing`): for
`t <= 4`, at most **3** squares of the packing have their centre within distance `1` of any one
wall.

**(C3) The regions partition the container**, so the interior count is `12 - Σ K - Σ n` and is
`>= 0`.

### 4.1 The corner level

Regions: the four corner boxes `A_0, A_1, A_2, A_3` in cyclic order (`A_0` bottom-left, `A_1`
bottom-right, `A_2` top-right, `A_3` top-left).  By (C1) each holds `0` or `1`, so the alphabet is
`{0,1}^4` and there are `16` patterns, **6 up to `D4`**: `0000, 1000, 1100, 1010, 1110, 1111`
(`search/BRANCH.md`).  The shared-`λ` form aggregates these into the five leaves `k = 0..4` of
`certificates/FORMAT.md`.  **No chord lemma is used at this level.**

### 4.2 The slot level under a corner leaf

Regions: additionally the eight wall slots `W_0..W_7` (cyclic; wall `w` owns `W_{2w}, W_{2w+1}` and
sits between corner boxes `A_w` and `A_{w+1 mod 4}`) and the interior `I`.  Let `n_j >= 0` be the
slot counts.  **This is where (C2) enters, and it is the only place it enters:**

```
for each wall w:      n_{2w} + n_{2w+1} + K_{A_w} + K_{A_{w+1}}  <=  3 .           (chord)
```

The bottom strip `[0,t] x [0,1]` is exactly `A_0 ∪ A_1 ∪ W_0 ∪ W_1` at `r = 1` (this is why `r = 1`
is the right corner size for the level-2 tree: it makes the corner boxes and the slots tile the
strip the chord lemma talks about).  Every centre in any of those four regions is within `1` of the
bottom wall, so (C2) bounds their total by `3`.  Note a corner box contributes to **two** walls.

Without (C2) the per-slot alphabet would be bounded only by (C1)-style diameter arguments, which do
not apply to a slot (a slot at `r = 1` has length `t/2 - 1 ≈ 0.99 > √3/2`, above the exactly
checkable threshold `4L^2 + (2r-1)^2 < 4` of `level2-design.md` §2.4), so there would be no finite
alphabet at all without either (C2) or Design B's shortened slots.

With (C2) the alphabet is `n_j ∈ {0,1,2,3}` and the wall constraint above; adding
`Σ K + Σ n <= 12` from (C3) gives the leaf set.

### 4.3 The counts

`search/branch_leaves.py` enumerates the pairs `(K, n)` satisfying (C1) + (C2) + (C3) and quotients
by `D4` acting on slots and corners together (the script checks that each group element carries
walls to walls, which is what makes the chord constraint `D4`-invariant); the orbit count is
cross-checked by Burnside.

```
$ python3 search/branch_leaves.py
slot cap 3, total <= 12
 corner leaf   k  slot leaves
        0000   0         1300
        0001   1         1830
        0011   2          546
        0101   2          351
        0111   3          171
        1111   4           15
       TOTAL             4213
```

* **`k = 4` gives 15**, reproducing `search/level2_leaves.py` and `notes/level2-design.md` §3: with
  both corner boxes of a wall occupied the chord constraint reads `n_{2w} + n_{2w+1} <= 1`, so the
  alphabet collapses to "at most one occupied slot per wall, holding one square", `3^4 = 81`
  patterns, `15` classes.  This is the count Design A quotes, and the chord lemma is exactly what
  earns it (Design B, which avoids the lemma by shortening the slots, has `43`).
* **`k = 1` has more leaves than `k = 0`** because the stabiliser of `0001` in `D4` has order `2`
  while that of `0000` has order `8`; the raw pattern counts are monotone the other way.
* Capping a slot at `2` centres — measured at `t = 3.98` (`level2-design.md` §2.3: three in one slot
  misses by `-0.44`) but **not proved** — gives `546 / 1176 / 438 / 351 / 171 / 15`, total `2697`.
  The rigorous number is `4213`.
* At `t = 3.98` only the `k = 4` subtree is needed (`k = 0,1,2` are closed by the shipped
  certificates and `k = 3` is a level-1 problem), i.e. **15 leaves**.  At `t = 4` the whole `4213`
  is needed, which at `level2-design.md` §8's `13–27 h` per cover-side leaf is not a computation
  anyone is going to run.  That is the real message of the count.

---

## 5. What is still checked outside Lean

`packing_le_weight_regions`, `packing_le_weight_regions_choice` and `wall_strip_le_three_of_packing`
are Lean theorems.  The following are not, at either level:

* **that the enumerated leaf set is exhaustive** — the combinatorics of §4.3, done in Python;
* **(C1)** — the corner-box diameter argument, done as an exact integer check inside `verify/`;
* **the identification of the verifier's `required` with the max rule** — i.e. that `region_flags`
  really returns a superset of the boxes a cell meets.  This is the padding argument of
  `certificates/FORMAT.md`, checked independently by `xcheck.py` in exact rationals, and it is the
  step (b) of §3.4 would change.

This is the same boundary as before; the new theorem moves the *reduction* for the boundary case
inside Lean, not the exhaustiveness.

---

## 6. The `m <= 1` leaves and a region-restricted capacity certificate

Under the corner leaf `k = 4`, a leaf with `m = Σ n` occupied slots leaves `8 - m` squares in the
interior `I = [r, t-r]^2`.  The two leaves `m = 0` (`00000000`) and `m = 1` (`00000001`) therefore
need `8` resp. `7` squares centred in `I`, and `notes/level2-design.md` §5 asks whether they are
refuted by interior capacity alone.  Three separate questions; the answers differ.

**(a) Interior-disjoint capacity at `t = 4` is at least 9, not 6.**  The nine axis-parallel unit
squares centred at `(i, j)`, `i, j ∈ {1,2,3}` occupy `[0.5, 3.5]^2 ⊆ [0,4]^2`, have all nine centres
in `I = [1,3]^2`, and have pairwise disjoint interiors (they tile, sharing edges).  Checked against
the repo's own margin function:

```
$ python3 -c "... from level2_capacity import margin_pair ..."
3x3 grid, centres in [1,3]^2, t=4: n = 9  min pairwise margin = 0.0
```

So `search/T4SCREEN.md`'s "six unit squares … with margin `+0.0567`" is a lower bound that badly
undersells the truth, and `level2_capacity.py --t 4.0 --box 1 3 1 3 --k 7`, which reports
`best min pairwise margin = -0.018354  (no k squares)`, is **wrong on this instance** — a
21-dimensional Nelder-Mead search missing the exact grid.  Any `alpha(I)` number from that script at
`t = 4` should be discarded.

**(b) The closed-disjoint capacity is a different question, and the LP route to it is dead at
`t = 4`.**  What a leaf actually needs is the *closed*-disjoint count, since the packing's
concentric unit squares are pairwise disjoint closed sets; the `3 x 3` grid does not qualify (its
squares touch).  A region-restricted cover certificate would settle it: a weighted point set of
total weight `< 7` such that **every closed unit square inside `[0,4]^2` whose centre lies in
`[1,3]^2` captures `>= 1`** proves "at most 6 such squares in any packing".  By LP duality such a
cover exists only if the corresponding fractional packing value is `< 7`.  `search/T4SCREEN.md` §0
measures exactly that value — the interior sub-problem, every frame region pinned to zero — and
reports it **converged** (coverage `M = 1`, clique max `1`, so the value is attained by an
admissible measure) at `7.433677`, `7.508327`, `7.528038`, still rising.  Therefore

> **no cover of total weight `< 7` for the interior region exists at `t = 4`**, and the `m <= 1`
> leaves cannot be closed by an interior capacity certificate there.

(Those LP numbers are float, not exact; the margin to `7` is `0.5`, far outside any plausible
float error, but this is not a machine-checked statement and is not shipped as a certificate.)

**(c) At `t = 3.98` it is open, and the format cannot say it.**  There the interior is
`[1, 2.98]^2` and the grid argument fails: three columns would sit at `1, 1.99, 2.98`, spacing
`0.99 < 1`, so the analogous nine squares overlap by `0.01`.  `level2-design.md` §5's numbers
(6 fit with `+0.0884`, 7 miss by `-0.0229`) are from the same heuristic that failed in (a) and
should be treated as weak evidence; the honest status of "at most 6 unit squares centred in
`[1, 2.98]^2` inside `[0, 3.98]^2`" is **open**.  The fractional value there is not measured; the
`t = 4` value of `7.53` and the `0.15–0.25` per `0.02` of `t` slope of `search/BRANCH.md` suggest it
is near `7`, i.e. right at the threshold, which is exactly the regime where only an exact
certificate would settle it.

**No certificate is shipped**, because the format cannot express the claim and the rules of this
task allow `certificates/` only for files `verify/` accepts.  What the claim needs:

* the claim is `n + λ·k_out <= W` with **one region** `R_out` = "centre outside `[1,3]^2`" and
  `λ = -1`: the threshold is `1` inside the interior box and `0` outside, and the conclusion
  `k_in = n - k_out <= W` gives `k_in <= 6` for `W < 7`;
* `R_out` is not a box.  The **minimal format extension** is `notes/level2-design.md` §7's
  list-of-boxes trailer, `regions Q g` followed by `L_j K_j b_j` and box lines, with `R_out`
  written as the four boxes `[0,s] x [0,1]`, `[0,s] x [3,s]`, `[0,1] x [1,3]`, `[3,s] x [1,3]` and a
  common `λ = -1`.  Those four boxes share boundary segments, so the trailer is *exactly* the case
  that needs §3.4(b)'s max rule and `packing_le_weight_regions_choice`: a centre on `y = 1` may be
  booked to either the strip or the interior, the certificate is valid for both, and `Σ_j λ_j k_j`
  is `-k_out` under either booking.  Nothing else in `verify/` changes — the cell-versus-box test is
  the one already there.
* Producing the certificate is then an ordinary cutting-plane run with the pose set restricted to
  centres in the interior box; `search/branch.py`'s machinery does everything except the box list,
  and `search/lp_search.py` would need the same restriction.  Neither is attempted here.

---

## 7. Reproduce

```sh
# leaf counts (section 4.3), seconds
python3 search/branch_leaves.py                 # rigorous: chord lemma only
python3 search/branch_leaves.py --slotcap 2     # with the measured (unproved) slot cap
python3 search/level2_leaves.py                 # the k = 4 subtree, 15 leaves, unchanged

# the 3x3 interior grid of section 6(a)
python3 -c "import sys; sys.path.insert(0,'search'); \
from level2_capacity import margin_pair; import itertools; \
P=[(i,j,0.0) for i in (1,2,3) for j in (1,2,3)]; \
print(min(margin_pair(p,q) for p,q in itertools.combinations(P,2)))"

# the Lean additions (warm Mathlib cache: ~10 s)
cd lean && lake build && lake env lean Axioms.lean
```
