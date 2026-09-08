# `leaf_ceiling.py`: an exact no-go certifier for a level-2 leaf

**Verdict, up front.**  The instrument is built and validated; **no no-go was found**.  It
reproduces both published exact `nu_f(3.99)` certificates to the digit from a code path that shares
nothing with the one that produced them; it raises the best certifiable anchor clique on those two
measures from `1.220095` and `1.216331` to **`1.274316176`** and **`1.254313293`** — exact, and
above the *unrestricted*-clique number `1.2139` that `CLIQUE_CONTINUUM.md`'s own branch and bound
found; and at `t = 4` it turns `search/T4SCREEN.md`'s corner-leaf `LP = 12.000000, M = 1.009` — a
value on a row set that 55 generation loops never converged, hence a bound on nothing — into an
exact statement: the largest **exactly certified** measure at `t = 4` with corner masses exactly
`1, 1, 1, 1` and coverage `<= 1` everywhere has mass **`11.939983685`**, `0.060` short of 12, and it
violates the anchor-clique family by `0.29`.  So the corner branch alone at `t = 4` is *not* proved
unable to close — the evidence points the other way, and it now points there exactly rather than
heuristically.

Code: `search/leaf_ceiling.py` (self-contained: no import from `packing_dual.py`, `dual_exact.py`,
`clique_ceiling.py`, `t4screen.py`, so the certification shares no code path with the search that
produces the measures).  Tests: `python3 search/leaf_ceiling.py selftest`.
Runs and logs: this worktree's `runs/lc_*`.

## 0. What this is for

`search/T4SCREEN.md` §7 item 4 names the only decisive experiment available on the packing side at
`t = 4`:

> Exhibit, in one leaf, a measure exactly certified feasible (coverage `<= 1` at every exact
> arrangement vertex, `mu(K) <= 1` for every anchor clique, region masses exactly the leaf's
> counts) with mass `>= 12`.

Such a measure is a **theorem-grade no-go**: by weak LP duality (the argument of
`search/CEILING.md`, with the leaf's region multipliers and the anchor-clique columns of
`notes/clique-family.md`), the weight of *every* certificate built from point cliques, anchor
cliques `K(p, A)` and per-region multipliers is at least the mass of any measure feasible for
those constraints.  A certified mass `>= 12` in a leaf therefore proves that **that leaf does not
close at that `t` with that family**, whatever pose set, cell decomposition or verifier is used.

This file is the instrument.  It takes a finite measure on poses of closed unit squares in
`[0,t]^2` (`t` rational) and certifies, with every load-bearing comparison in Python integers or
`Fraction`s and floats used only to *choose* candidates:

| | what | exact? | complete? |
|---|---|---|---|
| (a) | `cov(x) <= 1` at every point `x` of the container | yes | yes |
| (b) | `mu(K(p, A)) <= 1` for every anchor clique of the certifiable family | yes | see §2.5 |
| (c) | the mass centred in each closed leaf region is exactly the leaf's count, boundary poses assignable either way; optionally the chord inequality | yes | yes |
| | the mass itself | yes | — |

## 1. Coverage (a): why the arrangement vertices suffice

Let `S_1, …, S_n` be closed unit squares with masses `m_i > 0`, and `cov(x) = sum_{i : x in S_i} m_i`.
Let `V` be the set of **arrangement vertices**: every corner of every `S_i`, and every intersection
point of two closed edges of two distinct squares.

> **Lemma 1.**  For every subcollection `D`, `max_x cov_D(x) = max_{v in V} cov_D(v)`.

*Proof.*  `cov_D` is a finite sum of indicators of *closed* convex sets, so it is upper
semicontinuous and takes finitely many values; let `x*` attain the maximum and `S = {i in D : x* in S_i}`.
Every point of `P = ∩_{i in S} S_i` lies in all of them, so `cov_D >= cov_D(x*)` on `P`, and `P` is a
nonempty compact convex polygon, hence has a vertex `v`.  At a vertex of an intersection of convex
polygons two distinct supporting lines meet, each carrying an edge (a closed segment) of some
`S_i`, `i in S`.  Either both edges belong to the same square — `v` is a corner of it — or to two
squares — `v` is the intersection of two non-parallel closed edges.  Both kinds are in `V`.  (If
`S` is empty the maximum is `0`.)  ∎

This is `search/DUAL_EXACT.md`'s argument; it is repeated here because §2 needs it for arbitrary
subcollections `D`, not only for the whole measure.

**Exactly.**  Each pose is `theta = 2 arctan(p/q)` with `cos = (q^2-p^2)/(q^2+p^2)`,
`sin = 2pq/(q^2+p^2)` and a rational centre, so the four corners are rational; every square is
checked to lie in the closed container in integers.  Edge crossings are computed as reduced integer
triples `(X, Y, D)` (a `t`-parameter cross-product solve, kept iff both parameters are in `[0,1]`
— closed segments).  Containment of `(X/D, Y/D)` in the square with centre `(CX/Dk, CY/Dk)` and
`cos = a/r`, `sin = b/r` is the integer test

    2 |a dx + b dy| <= r D Dk    and    2 |a dy - b dx| <= r D Dk,
    dx = X Dk - CX D,  dy = Y Dk - CY D.

Masses are integers over a common denominator `DM`, so `cov(v) * DM` is an exact integer and `M =
max_V cov` is an exact `Fraction`.  Both the edge-pair enumeration and the incidence evaluation are
prefiltered by an **exact** integer grid (`floor(coord * G)` by integer division), so the prefilter
cannot drop a real crossing.  Nothing is D4-reduced: the whole container is enumerated, because a
leaf's slot pattern is not symmetric.

## 2. Anchor cliques (b): the exact maximiser

### 2.1 The family and the objective

For a point `p` and a nonempty compact convex `A`,

    K(p, A) = { S : p in S and S ∩ A ≠ ∅ }  ∪  { S : A ⊆ S }

is a clique for every `p` and every `A` (`notes/clique-family.md` Lemma 0: two members of the first
part share `p`, two of the second share `A`, and a mixed pair `S`, `S'` has `S ∩ S' ⊇ S ∩ A ≠ ∅`).
A certificate can carry `A` a **point** (`anchorP`) or a **closed segment** (`anchorS`)
— `certificates/FORMAT.md` has nothing else — so that is the family certified here.

Write `Cont(A) = {i : A ⊆ S_i}`, `Meet(A) = {i : S_i ∩ A ≠ ∅}` (so `Cont ⊆ Meet`) and
`P(p) = {i : p in S_i}`.  Then

    K(p, A) = Cont(A) ∪ (Meet(A) ∩ P(p)),
    mu(K(p, A)) = mu(Cont(A)) + mu( (Meet(A) ∩ P(p)) \ Cont(A) ).                (∗)

Two degenerate cases are settled at once and are *not* part of the search:

* **`A` a point `a`.**  Then `{S : A ⊆ S} = P(a)` already contains the first part, so
  `K(p, {a}) = P(a)` and `mu = cov(a) <= M`.  **Point anchors are exactly the coverage constraint**
  — check (a) covers them.
* **`Cont(A) = ∅`.**  Then `K(p, A) ⊆ P(p)` and `mu <= cov(p) <= M`.

So a violation needs a genuine segment `A` contained in at least one square, and the whole question
is a maximisation over `(p, A)`.

### 2.2 Reduction to one dimension

Two reductions make the objective computable.

> **Lemma 2 (the point).**  For a fixed `A`, `max_p mu(K(p, A))` is attained at some `v in V`.

*Proof.*  By (∗) the `p`-dependent term is `cov_D(p)` for the subcollection `D = Meet(A) \ Cont(A)`;
apply Lemma 1.  The vertices of the arrangement of `D` are a subset of `V`.  ∎

> **Lemma 3 (the segment, on a fixed line).**  Fix a line `ℓ` and let `I_k = ℓ ∩ S_k = [lo_k, hi_k]`
> (empty for the squares `ℓ` misses).  For a segment `A = [s, u] ⊆ ℓ`,
>
>     Cont(A) = { k : lo_k <= s  and  u <= hi_k },      Meet(A) = { k : lo_k <= u  and  s <= hi_k }.
>
> The maximum of `mu(K(p, A))` over segments of `ℓ` is attained with `s` and `u` both among the
> `2n` breakpoints `{lo_k, hi_k}` — or with `A` a point, which §2.1 has already settled.

*Proof.*  Both descriptions are immediate from convexity (`A ⊆ S_k` iff both endpoints are in
`S_k`) and from `[s,u] ∩ [lo_k,hi_k] ≠ ∅ ⟺ lo_k <= u and s <= hi_k`.  Let `A* = [s*, u*]` be optimal
with the optimal `p`, and let `U = Meet(A*) ∩ P(p)`, `C = Cont(A*)`.  Put

    u0 = max_{i in U} lo_i,        s0 = min_{i in U} hi_i.

Because `A*` meets every `I_i`, `i in U`, we have `u0 <= u*` and `s0 >= s*`.  If `s0 <= u0`, the
segment `A' = [s0, u0] ⊆ A*` still meets every `I_i` (`s0 <= hi_i` and `lo_i <= u0` by
construction), and `A' ⊆ A*` gives `Cont(A') ⊇ C`; hence `K(p, A') ⊇ K(p, A*)` and `A'` is at least
as good, with **both endpoints breakpoints** (`u0 = lo_{i1}`, `s0 = hi_{i2}` for some `i1, i2 in U`).
If `s0 > u0` then any point `x in [u0, s0] ⊆ A*` lies in every `I_i`, `i in U`, so the *point* anchor
`{x}` gives `K(p, {x}) ⊇ K(p, A*)` and the value is `<= M` (§2.1).  ∎

Lemma 3 turns the segment search on a line into a finite two-index problem, and — this is what
makes it fast — the two masses `mu(Cont)` and `mu(Meet)` for **all** pairs of breakpoints at once
are a two-dimensional prefix sum over the matrix `A[a][b] = ` (mass of the squares whose interval
has `lo` of rank `a` and `hi` of rank `b`): with `i`, `j` the ranks of `s`, `u`,

    mu(Cont(i,j)) = sum_{a <= i} sum_{b >= j} A[a][b],        mu(Meet(i,j)) = mu(Cont(j,i)).

All of this is integer arithmetic on the mass numerators.  Only the `p`-maximisation of Lemma 2 is
expensive, and it runs only on the pairs whose upper bound
`mu(Cont) + min(M, mu(Meet) - mu(Cont))` beats the current incumbent.

### 2.3 The candidate lines — the completeness theorem

> **Theorem.**  Let `V` be the arrangement vertex set.  Then
>
>     max over all points p and all point-or-segment anchors A  of  mu(K(p, A))
>       =  max( M ,  max over  p in V,
>                              lines ℓ through two distinct points of V,
>                              pairs (s, u) of breakpoints of ℓ           of  mu(K(p, [s,u])) ).
>
> Moreover in the second term one may require that the two points of `V` spanning `ℓ` lie **in a
> common square `S_j` of the measure** — which is how `leaf_ceiling.exhaust` organises the search:
> an outer loop over the "primary container" `j`, an inner loop over pairs of arrangement vertices
> inside `S_j`.

*Proof.*  By §2.1 we may assume `A` is a genuine segment with `C = Cont(A) ≠ ∅`, and by Lemma 2 the
point `p` may be taken in `V`; fix it, and let `U = Meet(A) ∩ P(p)`.  If `U ⊆ C` the value is
`mu(C) <= cov(x) <= M` for any `x in A`, which is case (i); so assume `U \ C ≠ ∅`.

**Step 1 (grow the segment to a chord).**  Let `R = ∩_{k in C} S_k`, a nonempty compact convex
polygon containing `A`, let `ℓ` be the line of `A` and `A' = ℓ ∩ R ⊇ A`.  Then `A' ⊆ R ⊆ S_k` for
every `k in C`, so `Cont(A') ⊇ C`; and `A ⊆ A'` gives `Meet(A') ⊇ Meet(A)`.  Hence
`K(p, A') ⊇ K(p, A)` and `A'` is at least as good.  So we may assume `A = ℓ ∩ R`, and the value
depends only on `ℓ`.

**Step 2 (push the line onto two tangencies).**  Put `T_i = S_i ∩ R` for `i in U`; these are
nonempty (each `S_i` meets `A ⊆ R`) compact convex polygons, and `ℓ` meets every one of them.  Let
`N` be the set of lines meeting every `T_i`, `i in U`; `N` is closed in line space, and for **every**
`ℓ' in N` the chord `A'' = ℓ' ∩ R` satisfies `Cont(A'') ⊇ C` and `Meet(A'') ⊇ U`, so
`mu(K(p, A'')) >= mu(K(p, A))`.  We may therefore replace `ℓ` by any convenient member of `N`:

* Translate `ℓ` along its normal until it is about to leave `N`.  At that position `ℓ` **supports**
  some `T_{i0}`: `ℓ ∩ T_{i0}` is a face of `T_{i0}` — either an edge (then `ℓ` contains two vertices
  of `T_{i0}`, and we are done) or a single vertex `a`.
* In the latter case rotate `ℓ` about `a` until it is about to leave `N`.  At that position `ℓ`
  supports some `T_{j0}`, touching it along an edge (two vertices) or at a vertex `b`.  If `b ≠ a`
  we are done.  If a full rotation about `a` never leaves `N`, then every line through `a` meets
  every `T_j`, `j in U`; but a compact convex set not containing `a` misses some line through `a`,
  so `a in T_j` for all `j`, and then the line through `a` and any other vertex of any `T_j` is in
  `N` and we are done.  (If every `T_j` equals `{a}`, the point anchor `{a}` already gives the value
  and §2.1 applies.)

So `ℓ` may be taken through two distinct vertices of the polygons `{T_i}` and `R`.

**Step 3 (those vertices are arrangement vertices).**  `R = ∩_{k in C} S_k` and `T_i = S_i ∩ R` are
intersections of closed unit squares.  A vertex of an intersection of convex polygons is a point
where two distinct supporting lines meet, each carrying an edge of one of the squares involved —
that is, a corner of a square or a crossing of two closed square edges.  Both are in `V`.  And all
of them lie in `R ⊆ S_j` for every `j in C`, and `C ≠ ∅`, which is the "common square" claim.

**Step 4 (the endpoints).**  Apply Lemma 3 on the final line `ℓ`: the optimal segment may be taken
with both endpoints at breakpoints of `ℓ` (or the value is `<= M`).  ∎

Two remarks.

* Step 2 is the only place where the *convexity* of `A` and the two-piece shape of `K(p, A)` are
  used; it is the reason the candidate lines are `O(|V|^2)` and not `O(|V|^4)` (a naive
  four-parameter argument on the two endpoints would pin each endpoint separately).
* The theorem also says what a *heuristic* separator can miss.  `anchorsep.separate` and
  `clique_family.anchor_local` sample `(p, direction, offset, half-length)` on grids; the exact
  optimum sits on a line through two arrangement vertices with endpoints at square crossings, which
  no such grid contains.  §4 measures how much that is worth.

### 2.4 The sound over-approximation: the maximum-weight clique

Every `K(p, A)` is a clique of the **closed-intersection graph** of the support squares (two closed
squares that touch count as adjacent — the semantics of `search/CLIQUE_CEILING.md`, because a
certificate at `t` refutes packings at `t' < t`, which rescale to pairwise disjoint *closed* unit
squares).  Hence

    max over the anchor family  <=  max-weight clique of the closed-intersection graph,

and the right-hand side is computed exactly: the graph in integers (four-edge-normal separating
axis test with `<=`, over the integer corner coordinates), the maximum by weighted branch and bound
with a colouring bound on **integer** masses, reporting whether it finished.  **If that maximum is
`<= 1`, (b) is proved** — and more than proved, since the measure is then feasible for every clique,
not only the certifiable ones.  This is complete, exact, and cheap enough for every measure here.

### 2.5 The three instruments, and their reach

| mode | what it computes | direction |
|---|---|---|
| `--anchor clique` | exact max-weight clique of the closed-intersection graph (integer weights, branch and bound, with a completeness flag) | **upper** bound on the family max; `<= 1` PROVES (b) |
| `--anchor exhaust` | the complete candidate enumeration of §2.3, with a time budget and two sound prunes | the exact family max **when it reports `complete`**; a lower bound otherwise |
| `--anchor scan` | floats choose candidates (anchor point over the arrangement vertices of high coverage, 24 directions, offsets and half-lengths on a grid — `clique_scan.py`'s parameterisation with the offset cap lifted, `search/RECONCILE.md` §2), every candidate re-evaluated in exact rationals, then an **exact local ascent** over the candidate set of §2.3 (re-solve the 1-D problem on the current supporting line; move an endpoint to a nearby arrangement vertex) | **lower** bound with an exact witness `(p, A)` |

Left alone, `exhaust` is `O(sum_j |V ∩ S_j|^2)` calls of the 1-D solver, and that is hopeless on a
real measure: a `t = 4` leaf measure with 405 squares already has `1.1 x 10^5` arrangement vertices,
a median of `1.3 x 10^4` of them inside each square, and `4.3 x 10^10` candidate lines in total.  Two
**sound** prunes cut that down, and they are what make the mode usable:

* **the primary container.**  If `j in Cont(A)` then the whole clique `K(p, A)` lies in
  `{j} ∪ N[j]` and contains `j`, so `mu(K(p, A)) <= ` the max-weight clique **through `j`** — an
  exact branch and bound on the neighbourhood of `j`, seeded with the current incumbent so that it
  usually terminates on the first bound.  Squares failing it are skipped outright.  The loop visits
  the squares in decreasing neighbourhood mass so the incumbent rises first.
* **the clique ceiling.**  If the incumbent ever reaches the exact max-weight clique of the whole
  graph (§2.4), then — since that number upper-bounds every `K(p, A)` — the maximum is attained and
  the search stops with `complete: true`.  This is the exit on the small cases (the three-square
  pinwheel of §4 takes it, and prints so).

**Measured reach.**  On the hand-made measures of §4 `exhaust` finishes in seconds and is a proof.
On the `t = 4` leaf measures of §5 (142–234 squares, `1.2 x 10^4`–`3 x 10^4` arrangement vertices,
2,000–3,900 of them inside a single square) it processes **one** primary square in 15–20 minutes —
about `2 x 10^6` candidate lines at 10 ms each — and reports `complete: false`.  It is therefore a
much better *search* than the grid scan (on `lc_K40polish` it lifts the best certified anchor clique
from `1.1717` to `1.4885`, §5.2) and not yet a proof of feasibility at that scale.  When it does not
close, the decisive statement for (b) is the clique upper bound of §2.4; on the `t = 3.99` measures
(1248 and 1920 squares, `10^6` vertices) `exhaust` is out of reach altogether and only `scan` and
`clique` are run.  Closing that gap is §7 item 1.

## 3. Region masses (c) and the boundary tie-break

The thirteen closed regions of a level-2 leaf (four corner boxes `[0,r]^2`, eight wall slots, the
interior `[r, t-r]^2`) follow `search/level2_regions.classify` exactly, and a pose is placed by its
**centre**.  The regions are closed boxes that cover the container, so a centre on a shared boundary
lies in two or four of them; `regions_of` returns all of them and the certifier then asks whether
**some** assignment of the boundary poses makes every pinned region's mass exactly its leaf count
(exact integer subset-sum, depth-first with a deficit bound; the boundary poses are few).  This is
the right semantics: the level-2 branch is a partition only through a tie-break at the boundary, so
a packing is free to declare a boundary pose in either adjacent region, and a measure must be
refuted for *every* declaration.  `search/T4SCREEN.md` §5 records that `1–11 %` of the mass of the
`m = 4` leaf optima — a whole unit at a single pose in two leaves — sits exactly on a slot boundary,
so this is not a corner case: it is where those optima live.

`--chord` additionally checks `mu(strip) <= 3` for the four wall strips `[0,t] x [0,r]` and images.

`leaf_ceiling.py snap` builds an exact measure from a float one (`t4screen.py`'s pose checkpoint or
a `packing_dual.py` support file): rational rotations `2 arctan(p/q)`, rational centres clamped
**exactly** into `[w/2, t - w/2]^2`, masses rounded **down** to multiples of `1/DM`, and — if the
exact maximum coverage still exceeds 1 — scaled by `1/M` and rounded down again.  Rounding down
leaves each pinned region a hair short of its count, and adding mass raises coverage, so the deficit
is spent greedily on the poses that have **exact** coverage slack: for a pose `i`,
`slack_i = min over the arrangement vertices of S_i of (1 - cov(v))`, and adding `d <= slack_i` to
pose `i` provably keeps `cov <= 1`.  The result is re-certified from the file by `check`, which
shares no state with `snap`.

## 4. Tests

`python3 search/leaf_ceiling.py selftest` (a few seconds).  Every case is a measure whose answer is
known by hand:

| case | expected | why it is there |
|---|---|---|
| two disjoint unit squares, mass 1 each | mass 2, `M = 1`, max clique 1 | the trivial certify |
| two squares sharing an edge, mass 1 each | `M = 2` | touching counts as covered on both sides |
| two overlapping squares, mass 1 each | `M = 2 > 1` | must NOT certify |
| four disjoint corner squares, mass 1 each | `M = 1`, exhaustive anchor max `= 1`, complete | the complete enumeration terminates and finds nothing |
| **pinwheel**: three unit squares at `0, 30, 60 deg`, centres at distance `3/5` from the centre in the directions `90, 210, 330 deg`, mass `1/2` each | `M = 1` (coverage feasible), max clique `3/2`, **exhaustive anchor max `3/2`, complete** | the discriminating case: they pairwise closed-intersect but have empty triple intersection (3-fold symmetric and convex, so a common point would force the centroid in, which is outside every square), so the coverage never sees the violation and the anchor family does |
| three unit squares in a row at `x = 1/2, 1, 3/2`, mass `1/2` | `M = 3/2`; the anchor `[0.75,0.5]–[1.25,0.5]` has mass `3/2`; the direct definition of `mu(K(p,A))` at `p = (1, 1/2)` agrees | the anchor evaluator against the definition, on two code paths |
| segment/square SAT, 400 random instances | exact and float agree everywhere | the three-axis test (the missing-axis bug of `notes/clique-family.md` §7 is not present) |
| a mid-wall pose `(2, 1/2)` | in both `W0` and `W1`; mass 1 there can satisfy `W0 = 1, W1 = 0` **or** `W0 = 0, W1 = 1`, and cannot satisfy `W0 = W1 = 1` | the boundary tie-break of §3 |
| `(1,1)`, `(3/2,1/2)`, `(1/2,1/2)` | `{C0, W0, W7, I}`, `{W0}`, `{C0}` | `level2_regions.classify` agreement |

The pinwheel also makes a three-line measure file, which is the shortest documentation of the input
format and of the clique-ceiling exit of §2.5 (`check` on it prints *"incumbent 1.500000000 attains
the exact clique ceiling 1.500000000; the maximum is proved"* and `EXHAUST … complete: true`):

```
# t = 4/1 ; sym 1
pose 0 1 2 13/5 1/2
pose 26795 100000 1480385/1000000 17/10 1/2
pose 57735 100000 2519615/1000000 17/10 1/2
```

## 5. Validation and results

Every measure below is **exactly certified**: coverage `<= 1` at every arrangement vertex, in
integers, and the region masses exactly the leaf's counts where a leaf is named.  The last column is
the best anchor clique found (an exact `Fraction`, hence a certified lower bound on the family's
maximum: `> 1` means the measure is **not** anchor-clique feasible).

| measure | `t` | leaf | exactly certified mass | `M` | best `mu(K(p,A))` (exact) | exact max-weight clique |
|---|---|---|---|---|---|---|
| `dual_exact_3.99_support.txt` | `399/100` | — | `12.008230754` | `7999999981/8000000000` | `1.274316176` | `>= 1.2139` (`CLIQUE_CONTINUUM.md`) |
| `cqx_PURE99_support.txt` | `399/100` | — | `12.028160771` | `3999999987/4000000000` | `1.254313293` | — |
| `lc_K40polish` (support of `t4_K40`) | `4` | corner `k = 4` | `11.853337054` | `1` | `1.488450644` | `1.652009445` (complete) |
| `lc_UNIONk4` (union of three supports) | `4` | corner `k = 4` | **`11.939983685`** | `1` | `1.294097438` | `1.559218311` (complete) |
| `lc_A0101` (hardest `m = 4` leaf) | `4` | `1111 / 01010101` | `11.776361678` | `1` | `1.282857187` | `1.538759710` (complete) |

**No measure reached 12, so no no-go was produced.**  The two `t = 4` leaves come in `0.060` and
`0.224` short, and every one of them violates the certifiable clique family by `0.28`–`0.49`, so
adding the anchor cliques would push them further down, not up.

### 5.1 The certified `nu_f(3.99)` measures — coverage reproduced exactly

Both published exact measures were re-certified from their support files by this implementation,
which shares no code with `dual_exact.py` (it re-derives the poses, the arrangement, the incidence
test and the maximum).  Every digit agrees:

| measure | poses / squares | arrangement vertices (whole container, no D4 reduction) | mass | exact `M` | published |
|---|---|---|---|---|---|
| `runs/dual_exact_3.99_support.txt` | 156 / 1248 | 836,576 | `6004115377/500000000 = 12.008230754` | `7999999981/8000000000` at `~(0.989619, 0.998349)` | `DUAL_EXACT.md`: identical |
| `runs/cqx_PURE99_support.txt` | 240 / 1920 | 2,042,328 | `12028160771/1000000000 = 12.028160771` | `3999999987/4000000000` at `~(0.999939, 0.999939)` | `CLIQUE_CONTINUUM.md` §3: identical |

(`dual_exact.py check` enumerates the fundamental domain and reports 104,947 and 255,893 vertices;
eight times those, minus the vertices fixed by a symmetry, is what this run enumerates over the
whole container.)

**The anchor cliques of those measures.**  `search/RECONCILE.md` §2 is the reference: on
`dual_exact_3.99_support.txt` a segment-anchor scan with the offset cap lifted reaches `1.220095`,
above the `1.2139` that `CLIQUE_CONTINUUM.md`'s own unrestricted-clique branch and bound found.
The exact ascent of §2.5 goes further, and every number below is an exact `Fraction` with an
exactly re-checked witness:

| measure | published best (float scan) | **this tool (exact)** | witness |
|---|---|---|---|
| `dual_exact_3.99_support.txt` (12.008230754) | `1.220095` (`RECONCILE.md` §2) | **`79644761/62500000 = 1.274316176`** | `p ~ (2.820522, 2.422532)`, wall distance `1.1695`, `\|A\| = 0.76677`, `A ~ (2.99163, 2.09932)-(2.99000, 2.86609)`; **188 members**, 83 of them containing `A` |
| `cqx_PURE99_support.txt` (12.028160771) | `1.216331` (`RECONCILE.md` §2) | **`10034506341/8000000000 = 1.254313293`** | `p ~ (2.360503, 2.745881)`, wall distance `1.2441`, `\|A\| = 0.77009`, `A ~ (2.91662, 2.99009)-(2.14654, 2.99000)`; **295 members**, 133 containing `A` |

So the certifiable one-point-one-segment family reaches **`1.2743`** on the measure that pins the
`nu_f(3.99) >= 12.008` bound and **`1.2543`** on the one that pins `12.028` — `+0.054` and `+0.038`
above the best previously recorded for that family, and above the `1.2139` that
`CLIQUE_CONTINUUM.md` §4's *unrestricted*-clique branch and bound reported.  The mechanism is
`RECONCILE.md`'s: the anchor point is **interior** (wall distance `1.17` and `1.24`, well outside
the band Lemma 1 talks about) and the anchor is a long segment (`0.77`) pressed flat against the far
wall.

**Both witnesses were re-checked from the definition**, on a code path that shares nothing with the
maximiser: rebuild `K(p, A)` by testing every support square for `p in S` and `S ∩ A ≠ ∅` (or
`A ⊆ S`), sum the masses, and test **every unordered pair of members** with the exact four-axis
square-square SAT.  `188` members and `17,578` pairs for the first, `295` and `43,365` for the
second; **`0` non-intersecting pairs** in both, and the mass recomputed from the definition equals
the maximiser's `Fraction` exactly.  (This is the check that catches the missing-axis bug of
`notes/clique-family.md` §7, which turned a `1.024` into a spurious `1.3308`.)  The same check now
runs inside `check --anchor scan` on every witness it reports.

### 5.2 The corner leaf `k = 4` at `t = 4` — what is exactly certified

`search/T4SCREEN.md` §3 records the corner leaf's *pure* value landing on exactly `12.000000` in
three independent runs.  That number is an LP value on a row set that never converged: re-run here
(`runs/K40P.out`, `t4screen.py` on the `K40` 5516-pose checkpoint, cliques off, corner equality on,
row generation only, row ageing off, `row-cap 30000`, 91 minutes), the LP sits at `12.000000` at
**every one of 55 row-generation loops**, the row set grows from 7,300 to 24,538, and the certified
maximum coverage `M` of its own solution wanders between `1.009` and `1.49` without ever reaching 1
(final state `LP = 12.000000, M = 1.009021956`).  Per §1.2 of `T4SCREEN.md` that is an **upper**
bound on the value of a restricted pose set and a bound on nothing else.

The exact instrument settles what those poses actually carry.  The exact-row LP with the four
corner equalities (`snap --polish --corners 1111`, §5.4) converges in 7 rounds on the support of
that checkpoint (405 poses) and in 7 rounds on the union of the positive-mass supports of `K40`,
`INTP` (the refined interior) and `L01010101` (1256 poses), and certifies:

| support | poses in | poses out | exactly certified mass | `M` | corner masses | chord (four wall strips) |
|---|---|---|---|---|---|---|
| `t4_K40` | 405 | 142 | `11853337054/1000000000 = 11.853337054` | `1` exactly | `1, 1, 1, 1` exactly | `3, 3, 2.973314, 3` |
| union `K40 + INTP + L01010101` | 1256 | 182 | `11939983685/1000000000 = 11.939983685` | `1` exactly | `1, 1, 1, 1` exactly | `3, 3, 3, 3` exactly |

So the corner branch alone at `t = 4` is **not** proved unable to close: the largest exactly
certified corner-leaf measure obtained here is `0.060` short of 12.  What `12.000000` was, on the
evidence here, is the *first* LP value in the row-generation sequence — the union polish opens with
exactly `12.000000` on 27,413 rows and 14,515 violated vertices, and falls to `11.9400` as the rows
close in.

The anchor cliques of those measures are far from feasible.  On the two certified measures the
exact local ascent reaches

| measure | best `mu(K(p, A))` (exact, a lower bound) | witness | exact max-weight clique of the whole graph (B&B **complete**) |
|---|---|---|---|
| `lc_K40polish` (11.8533, 142 squares) | scan+ascent `585831419/500000000 = 1.171662838`; **`exhaust` 20 min: `372112661/250000000 = 1.488450644`** | `p ~ (1.338944, 1.765805)`, `A ~ (1.78408, 0.97700)-(1.94961, 0.99894)` (the ascent witness) | `330401889/200000000 = 1.652009445` (size 32) |
| `lc_UNIONk4` (11.9400, 182 squares) | `647048719/500000000 = 1.294097438` (12 contain `A`, 31 meet it); `exhaust` 15 min did not improve it | `p ~ (1.502795, 1.999999)` — **on the slot boundary `y = 2`** — `A ~ (1.33723, 2.11314)-(1.13790, 1.99239)` | `1559218311/1000000000 = 1.559218311` (size 48) |

So these measures violate the certifiable clique family by `0.29` to `0.49`: adding the anchor
cliques to the leaf LP moves its value **down** from `11.94`, not up towards 12.  Two things are
worth reading off.  First, on `lc_K40polish` the complete enumeration finds `1.4885` where the
`24 x 12 x 10` grid scan of `clique_scan.py`'s family finds `1.1717` — the exact candidate set of
§2.3 is worth `+0.32` over the grid, and it reaches `90 %` of the *unrestricted* clique number.
Second, the `lc_UNIONk4` witness sits exactly on the slot boundary `y = t/2`, which is
`search/T4SCREEN.md` §5's obstacle showing up again, now on the clique side.

### 5.3 The hardest leaf `01010101` at `t = 4`

The leaf-search agent's checkpoint for the hardest `m = 4` leaf (`tl_A01010101_poses.txt`, 12,828
columns, float `LP = 11.746690` with `M = 1.002` and `kmax = 1.05` — again not converged) has 417
poses with positive mass.  Snapped, polished on the exact rows with **all twelve** region
equalities, and re-certified:

    mass  = 11776361678 / 1000000000 = 11.776361678     (234 poses with positive mass)
    M     = max coverage = 1  exactly                    (over the whole arrangement)
    corners  = 1, 1, 1, 1  exactly
    slots    = 0, 1, 0, 1, 0, 1, 0, 1  exactly           (the leaf's own pattern)
    exact max-weight clique of the closed-intersection graph = 153875971/100000000 = 1.538759710

`0.224` short of 12.  The polish reaches it in 5 row-generation rounds (`11.777457959` at 22,019
rows down to `11.776361751` at 22,175 rows, 0 violated vertices), i.e. the exact value on that
support is `0.03` **above** the float screen's `11.7467` — the screen was column-starved on its own
support, not row-starved.  Its best certified anchor clique is `1282857187/1000000000 =
1.282857187` at `p ~ (1.490136, 2.878035)` (4 squares contain `A`, 44 meet it), so it too violates
the family, by `0.28`.

### 5.4 A by-product worth stating: the exact-row polish

The check `(a)` enumerates the **whole** arrangement of a support, and that vertex set is complete
for every sub-measure of it (Lemma 1).  So the LP

    max sum mu   s.t.  cov(v) <= 1 at every arrangement vertex v,  region masses = the leaf counts

solved on the poses of a support is not a screen: its value is a value those poses **really
attain**, with no row generation to converge and no `M > 1` caveat.  `snap --polish` solves exactly
that (by row generation over the exact vertex set, since the whole set is `10^5`–`10^6` rows and
`10^8` incidences), rounds the masses down and tops the regions back up on exact coverage slack.
That is how §5.2 turns `t4screen.py`'s `LP = 12.000000, M = 1.02` — an upper bound on a restricted
value, and a bound on nothing — into an exact number in one minute instead of an hour of row
generation.  The price is that it is confined to the poses of the input support: adding columns is
what the leaf search is for.

## 6. Reproduce

```sh
python3 search/leaf_ceiling.py selftest                        # ~20 s

# (1) the certified nu_f(3.99) measures: coverage reproduced exactly, best anchor clique
python3 search/leaf_ceiling.py check $R/dual_exact_3.99_support.txt --anchor scan --procs 4 \
        --scan-pts 200 --scan-verify 40 --ascent-time 900
python3 search/leaf_ceiling.py check $R/cqx_PURE99_support.txt      --anchor scan --procs 4 \
        --scan-pts 200 --scan-verify 40 --ascent-time 900

# (2) the corner leaf k = 4 at t = 4.  The float side (a screen; never converges its rows):
python3 search/t4screen.py 4.0 K40P --corners 1111 --load-poses runs/src_K40_poses.txt \
        --stages 1 --rowloops 60 --row-age 0 --row-cap 30000 --pose-max 0 --time 5400
#     The exact side: snap + exact-row polish with the corner equalities, then certify.
python3 search/leaf_ceiling.py snap runs/src_K40_poses.txt --t 4 --out runs/lc_K40polish.txt \
        --corners 1111 --polish --procs 4
python3 search/leaf_ceiling.py check runs/lc_K40polish.txt --corners 1111 --chord \
        --anchor all --procs 4 --clique-time 900 --ascent-time 300 --budget 1200
#     the same on the union of the positive-mass supports of K40, INTP and L01010101 (1256 poses)
python3 search/leaf_ceiling.py snap runs/src_UNION_poses.txt --t 4 --out runs/lc_UNIONk4.txt \
        --corners 1111 --polish --procs 4

# (3) the hardest leaf 01010101 (checkpoint of the leaf-search agent)
python3 search/leaf_ceiling.py snap runs/src_A01010101_poses.txt --t 4 --out runs/lc_A0101.txt \
        --corners 1111 --patterns 01010101 --polish --procs 4
python3 search/leaf_ceiling.py check runs/lc_A0101.txt --corners 1111 --patterns 01010101 \
        --chord --anchor all --procs 4
```

## 7. What a follow-up needs

1. **Close the gap in (b) on a real measure.**  §2.5: the complete enumeration is
   `O(sum_j |V ∩ S_j|^2)` lines, and on a `t = 4` leaf measure that is `10^9`–`10^{10}` at 10 ms
   each.  Both prunes already in it (the exact max-weight clique through the primary container; the
   clique ceiling) fire only when the incumbent is high, and on these measures the incumbent stops
   `0.16–0.36` below the clique number.  The untried attack is a branch and bound **in line space**:
   for a box `B` of lines, `mu(Cont) + min(M, mu(Meet))` computed with `Cont` and `Meet` taken over
   all lines of `B` is an exact upper bound and is monotone under subdivision.  A second, cheaper
   one: the two spanning vertices both lie in `R_C`, so a pair `(v, w)` can be skipped when
   `mu(P(v) ∩ P(w)) + M` is at or below the incumbent — one sparse mat-vec per `v` gives that for
   every `w` at once.
2. **The instrument does not produce the measure.**  Everything here certifies; the mass still has
   to come from the leaf search, and §5 records how far short of 12 the available measures are
   (`0.060` in the corner leaf, `0.224` in the hardest `m = 4` leaf).  The exact-row polish of §5.4
   is confined to the poses of the support it is given, so the productive loop is: leaf search
   prices columns → `snap --polish` certifies exactly what those columns carry → repeat.  Because
   the polish is complete in the rows, a leaf value it produces cannot be argued away, which is more
   than can be said for any `LP/M` pair in `search/T4SCREEN.md` §3.
3. **The anchor cliques are the reason the leaf values will not rise.**  Every certified measure in
   §5 violates the certifiable family by `0.28–0.49` and the unrestricted one by `0.54–0.71`.  A
   clique-strengthened leaf LP therefore pushes the value **down** from `11.94`, and a `>= 12`
   measure would have to be simultaneously coverage-feasible, region-exact and anchor-feasible — the
   three constraints this file checks, and the reason it exists.
