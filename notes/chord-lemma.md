# The wall-strip chord lemma, proved

Lean: `lean/Sqpack/Chord.lean` (`wall_strip_le_three`, `wall_strip_le_three_of_packing`;
0 sorries, axioms `propext, Classical.choice, Quot.sound`).
Used by: `notes/branch-semantics.md` §4 (the leaf enumeration), `notes/level2-design.md` §2.3
and §10 item 1, `search/T4SCREEN.md` §7 item 6.
Sources: Stromquist 1984-I pp. 5–8; Nagamochi 2005 Lemma 7(i); `notes/proof-anatomy.md`
§3.3, §5.1, §7.4.

---

## 0. The statements

Throughout, a *unit square* is a closed square of side `1` at an arbitrary angle, `sq c θ 1` in the
Lean file; `[0,t]^2` is the closed container.

> **B1 (closed form).**  Let `t <= 4`.  Let `S_1, …, S_N` be unit squares contained in `[0,t]^2`,
> **pairwise disjoint as closed sets**, each with centre `c_i` satisfying `(c_i)_y <= 1`.
> Then `N <= 3`.

> **B2 (packing form).**  Let `t <= 4` and `L > 1`.  Let `S_1, …, S_N` be squares of side `L`
> contained in `[0,t]^2` with **pairwise disjoint interiors**, each with centre `c_i` satisfying
> `(c_i)_y <= 1`.  Then `N <= 3`.

B2 is B1 applied to the concentric unit squares: `unit_subset_interior` (already in
`lean/Sqpack/Basic.lean`) says `sq c θ 1 ⊆ sqInt c θ L` for `L > 1`, so disjoint interiors of the
`L`-squares give **disjoint closed** unit squares, and the centres are unchanged.  This is the
form the branch tree uses: a packing of unit squares in a container of side `s' < t` rescales by
`L = t/s' > 1`.

Both are stated for the bottom wall; the other three follow by the symmetries of `[0,t]^2`.
"Centre within distance 1 of the wall" is `(c_i)_y <= 1`, i.e. the centre lies in the closed strip
`[0,t] x [0,1]`.

**Why the hypotheses are exactly right.**

* *`t <= 4` and closed disjointness are both needed in B1.*  At `t = 4` the four axis-parallel
  unit squares centred at `(1/2 + i, 1/2)`, `i = 0..3`, lie in `[0,4]^2`, have centres at height
  `1/2`, and are pairwise disjoint **except that consecutive ones touch**.  So the interior-disjoint
  version of B1 is false at `t = 4`, and B2's `L > 1` is what rules this out.  For `t < 4` the
  interior-disjoint version is true (§4).
* *The lemma is a strictness statement*, in the sense of `notes/proof-anatomy.md` §7.4: it is true
  for boxes (open squares of side `> 1`) and false, at `t = 4`, for closed unit squares with
  margin `0`.  No eroded LP can see it, and `search/level2_capacity.py` cannot either — at
  `t = 3.98` it reports `-0.0067 = (4 - 3.98)/3` for four squares in the full wall strip, i.e. it
  measures the slack and reports "infeasible" only because `t < 4`; at `t = 4` it would report
  feasible.  That is why the lemma has to be supplied analytically.
* *The bound `3` is not about strips of height 2.*  Stromquist 1984-III p. 10 packs **four** unit
  squares into a `1.9 x 3.9475` rectangle (`notes/proof-anatomy.md` §5.1).  That configuration has
  an axis-parallel square at the top-left, whose centre is at height `1.4 > 1`, so it is not a
  counterexample: B1 constrains *centres* in a strip of height `1`, not squares in a strip of
  height `1.9`.

---

## 1. Notation and the chord

Write `C = cos θ`, `S = sin θ` for the angle of a square, and

```
w(θ) = |C| + |S| ∈ [1, √2],     p = w(θ)/2 ∈ [1/2, √2/2]
```

so `p` is the square's half-height (and half-width): the unit square with centre `c` and angle `θ`
has bounding box `[c_x - p, c_x + p] x [c_y - p, c_y + p]`, and the four bounding-box extremes are
attained at vertices.  Set `u = |C| + |S| = 2p`; then `u ∈ [1, √2]` and

```
u^2 = 1 + 2|C||S|,      i.e.   |C||S| = (u^2 - 1)/2 .
```

For a horizontal line `y = a`, the **chord** of the square is `S ∩ {y = a}`, an interval; write
`ℓ` for its length and `d = |a - c_y|` for the distance from the centre to the line.

---

## 2. The chord length

**Lemma 1 (chord length).**  For a unit square with `|C|, |S| > 0` and `d <= p`,

```
ℓ = min( 1/|C| ,  1/|S| ,  (p - d)/(|C||S|) ) .
```

For `|S| = 0` (axis-parallel) `ℓ = 1` when `d <= 1/2` and the chord is empty otherwise; likewise
for `|C| = 0`.

*Proof.*  Put the origin at the centre; a point at horizontal offset `ξ` and vertical offset `e`
from the centre lies in the square iff, in the square's own frame,

```
|ξ C + e S| <= 1/2      and      |e C - ξ S| <= 1/2 .                        (∗)
```

(This is the definition of `sq`: the two coordinates of `coord c θ p`; `sq_mem_iff_slab` in the
Lean file.)  Take `e = a - c_y`, so `|e| = d`.  With `C, S > 0` each condition is an interval in
`ξ`:

```
first:   ξ ∈ [ (-1/2 - eS)/C , (1/2 - eS)/C ]      of length 1/C
second:  ξ ∈ [ (eC - 1/2)/S , (eC + 1/2)/S ]       of length 1/S
```

so the chord is the intersection of two intervals and its length is the minimum of the four
differences (upper endpoint minus lower endpoint):

```
(1/2 - eS)/C - (-1/2 - eS)/C = 1/C
(eC + 1/2)/S - (eC - 1/2)/S  = 1/S
(1/2 - eS)/C - (eC - 1/2)/S  = [ S(1/2 - eS) - C(eC - 1/2) ] / (CS) = ( (C+S)/2 - e ) / (CS)
(eC + 1/2)/S - (-1/2 - eS)/C = [ C(eC + 1/2) + S(1/2 + eS) ] / (CS) = ( (C+S)/2 + e ) / (CS)
```

using `C^2 + S^2 = 1`.  The last two are `(p ∓ e)/(CS)`, whose minimum is `(p - |e|)/(CS)`.  The
sign cases reduce to `C, S > 0` because replacing `(C, S)` by `(-C, -S)` leaves both expressions
in (∗) unchanged up to sign, and replacing `(C, S, e)` by `(C, -S, -e)` does the same
(`slab_neg_neg`, `slab_neg_snd` in the Lean file).  The degenerate cases are immediate from (∗). ∎

Two checks against the literature.  At `θ = 45°` the chord is `≥ 1` exactly for
`d <= p - CS = √2/2 - 1/2 = (√2 - 1)/2`, which is Nagamochi's Lemma 2 threshold.  At `θ = 0` the
chord is exactly `1` for every `d <= 1/2`: the "`> 1`" of Nagamochi's Lemma 7(i) is a statement
about squares of side `λ > 1`, and for unit squares the inequality is not strict.

**Corollary 2 (when the chord is long enough).**  `ℓ >= 1` if and only if

```
d  <=  D(θ) := p - |C||S| = u/2 - (u^2 - 1)/2 = ( u - u^2 + 1 ) / 2 .
```

*Proof.*  `1/|C| >= 1` and `1/|S| >= 1` always, so by Lemma 1 the binding term is
`(p - d)/(|C||S|) >= 1`, i.e. `d <= p - |C||S|`.  (When `|C||S| = 0` the condition reads
`d <= 1/2 = D`, matching the degenerate case.)  Note `D <= p`, with equality only at `θ ≡ 0`. ∎

`D` decreases from `1/2` at `u = 1` to `(√2 - 1)/2 ≈ 0.2071` at `u = √2`.

---

## 3. Choosing the cut height

A unit square inside `[0,t]^2` has `c_y >= p` (its lowest vertex is at height `c_y - p >= 0`), and
the hypothesis gives `c_y <= 1`.  Note `p <= √2/2 < 1`, so the interval `[p, 1]` is non-empty.

**Lemma 3 (admissible cut heights).**  Let `a` satisfy

```
(3 - √2)/2  <=  a  <=  √2 - 1/2          i.e.   0.792893… <= a <= 0.914214…
```

Then every unit square contained in `[0,t]^2` whose centre has `c_y <= 1` meets the line `y = a`
in a chord of length `>= 1`, and the line meets the **interior** of the square.

*Proof.*  By Corollary 2 the chord condition is `|a - c_y| <= D`, i.e. `c_y ∈ [a - D, a + D]`.
Since this is an interval condition in `c_y` and `c_y` ranges over `[p, 1]`, it suffices to check
the two endpoints.

*Endpoint `c_y = p = u/2`.*  We need `|a - u/2| <= u/2 - (u^2-1)/2`.
 - If `a >= u/2`: `a - u/2 <= u/2 - (u^2-1)/2` ⟺ `a <= u - (u^2-1)/2 =: f(u)`.
   `f'(u) = 1 - u <= 0` on `[1, √2]`, so `min f = f(√2) = √2 - 1/2`.
 - If `a < u/2`: `u/2 - a <= u/2 - (u^2-1)/2` ⟺ `a >= (u^2-1)/2 = |C||S| <= 1/2`, which holds since
   `a >= 0.7929 > 1/2`.

*Endpoint `c_y = 1`.*  Since `a < 1`, we need `1 - a <= u/2 - (u^2-1)/2`, i.e.
`a >= 1 - u/2 + (u^2-1)/2 = (u^2 - u + 1)/2 =: g(u)`.  `g` is increasing on `[1, √2]`, so
`max g = g(√2) = (3 - √2)/2`.

So `g(√2) <= a <= f(√2)` is exactly the stated range, and it is non-empty because
`(3-√2)/2 < √2 - 1/2` ⟺ `3 - √2 < 2√2 - 1` ⟺ `4 < 3√2` ⟺ `√2 > 4/3`. ✓

*The line meets the interior.*  We need `d < p` strictly.  Upwards:
`c_y - a <= 1 - a <= 1 - (3-√2)/2 = (√2-1)/2 < 1/2 <= p`.  Downwards: `a - c_y <= a - p`, and `a - p < p` ⟺ `a < u`, which holds since
`a <= 0.9143 < 1 <= u`. ∎

**`a = 9/10` is admissible** (`0.7929 <= 0.9 <= 0.9142`), which is the classical choice — Stromquist's
and Nagamochi's `y = 0.9`.  The Lean proof uses `9/10`; the two endpoint inequalities become, after
clearing `u`,

```
u^2 - 2u + 4/5 <= 0        (from c_y = p)          ⟺  (u-1)^2 <= 1/5
u^2 -  u - 4/5 <= 0        (from c_y = 1)
```

both of which hold on `1 <= u <= √2` — the Lean proof discharges them from `1 <= u <= 1.415` and
`0 <= (u-1)(1.415 - u)`.

**Remark (how much room there is).**  The same argument with `c_y <= h` instead of `c_y <= 1`
works as long as `h - a <= min_u (u/2 - (u^2-1)/2) = (√2-1)/2`, so with `a = √2 - 1/2` the strip of
centres may be taken as high as `h = 3√2/2 - 1 = 1.12132…`.  The height `1` we use has `0.12` of
slack; the lemma is not tight at `h = 1`.

---

## 4. Proof of B1 and B2

Fix `a = 9/10` and let `I_i = S_i ∩ {y = a}` be the chords, `I_i = [α_i, β_i]`.

By Lemma 3, `β_i - α_i >= 1` for every `i`.  Since `S_i ⊆ [0,t]^2`, `I_i ⊆ [0,t] x {a}`, so
`0 <= α_i` and `β_i <= t`, hence

```
α_i ∈ [0, t - 1] ⊆ [0, 3]        (using t <= 4).
```

**B1.**  The `S_i` are pairwise disjoint, so the `I_i` are pairwise disjoint closed intervals.  If
`i ≠ j` and, say, `α_i <= α_j`, then `α_j ∉ I_i` (it is in `I_j`), and since `α_j >= α_i` this
forces `α_j > β_i >= α_i + 1`.  So the `α_i` are pairwise **more than 1 apart**.

Suppose `N >= 4`.  Four points of `[0,3]` pairwise more than `1` apart are impossible: order them
`α_(1) < α_(2) < α_(3) < α_(4)`; then `α_(4) > α_(3) + 1 > α_(2) + 2 > α_(1) + 3 >= 3`,
contradicting `α_(4) <= 3`.  Hence `N <= 3`. ∎

(The Lean proof avoids sorting: it takes any four indices and case-splits each of the six
`|α_i - α_j| > 1` into its two directions, `linarith` closing all 64 branches.)

**B2.**  `unit_subset_interior` turns disjoint interiors of the `L`-squares into disjoint closed
concentric unit squares, and containment `sq c θ L ⊆ [0,t]^2` gives `sq c θ 1 ⊆ [0,t]^2`; the
centres are unchanged.  Apply B1. ∎

**The interior-disjoint version at `t < 4`.**  If the `S_i` are unit squares with pairwise disjoint
*interiors* and `t < 4`, then the open chords `int(S_i) ∩ {y = a}` are pairwise disjoint open
intervals of length `>= 1` (Lemma 3 gives that the line meets each interior, and for a compact
convex set met by a line in its interior the open and closed chords have the same length),
contained in `[0,t]`.  Four of them would need total length `>= 4 > t`.  So `N <= 3`.  This is the
form that fails at `t = 4`, and it is the reason B1 is stated with closed disjointness.

---

## 5. The corollary the tree uses

**Corollary 4 (wall strip capacity).**  In a packing of `12` squares of side `L > 1` with pairwise
disjoint interiors in `[0,t]^2`, `t <= 4`, at most `3` squares have their centre in the closed
strip within distance `1` of any one wall.

With the level-2 regions of `notes/level2-design.md` §2.1 at `r = 1` — corner boxes `A_w` (cyclic)
and wall slots `W_0..W_7` — the bottom strip `[0,t] x [0,1]` is exactly `A_0 ∪ A_1 ∪ W_0 ∪ W_1`,
and similarly for the other three walls.  So for each wall `w` (slots `2w, 2w+1`, corner boxes
`A_w, A_{w+1 mod 4}`),

```
n_{2w} + n_{2w+1} + K_{A_w} + K_{A_{w+1}}  <=  3 .
```

In the corner leaf `K = 1111` this reads `n_{2w} + n_{2w+1} <= 1`: **at most one occupied slot per
wall, holding exactly one square** — which is what cuts the level-2 leaf count from `43` (Design B,
`2^8` slot patterns) to `15`.  `notes/branch-semantics.md` §4 does the enumeration and says exactly
where this is used; the counts for the other corner leaves are there too.

A corner box counts towards **two** strips (its centre is within `1` of two walls), which is why
the four constraints are not independent; summing them gives the global frame bound
`Σ n + 2 Σ K <= 12`.

---

## 6. What is formalised

`lean/Sqpack/Chord.lean`, 0 sorries:

| name | statement |
|---|---|
| `box t` | the container `[0,t]^2` |
| `slab P Q e` | the set of horizontal offsets `ξ` satisfying (∗); `sq_mem_iff_slab` ties it to `sq` |
| `slab_neg_neg`, `slab_neg_snd` | the two sign symmetries that reduce to `cos θ, sin θ >= 0` |
| `exists_unit_interval_subset_slab` | **Lemma 1 + Corollary 2**, in the only direction needed: if `|e| <= (P+Q)/2 - PQ` then the slab contains a closed interval of length `1` |
| `half_height_le_centre` | `(\|cos θ\| + \|sin θ\|)/2 <= c_y` for a square inside the box |
| `exists_chord` | **Lemma 3** at `a = 9/10`: an interval of length `1` of points of the square on the line `y = 9/10` |
| `no_four_spread` | four points of `[0,3]` cannot be pairwise more than `1` apart |
| `wall_strip_le_three` | **B1** |
| `wall_strip_le_three_of_packing` | **B2** |

The formalisation proves a *lower bound* on the chord (it exhibits a sub-interval of length exactly
`1`), which is all the counting needs; the exact chord length of Lemma 1 is not formalised.
`lake build` is ~8 s for this file on a warm Mathlib cache.

**What is still outside Lean.**  That the *leaf set* is exhaustive — the enumeration of §5's
corollary into a finite list of `(K, n)` patterns — is combinatorics done in Python
(`search/branch_leaves.py`), not in Lean; the same is true of the level-1 statement
`(2r-1)^2 < 2 ⟹ at most one centre per corner box`, which lives in `verify/`.  See
`notes/branch-semantics.md` §5.
