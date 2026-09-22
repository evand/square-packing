# counting-ladder: the chain-counting lemma proved, the three-level ladder in closed form, and what fails at `T = 11`

*2026-09-21.  Task `tasks/counting-ladder/README.md`.  New code, nothing in `search/` modified:
`search/ladder_hform.py` (imports `search/t3_chain_hform.py`; also imports `t3_chain.py` and
`t11_chains.py` for the `T = 11, 12` sections).  Runs: `runs/ladder_check.txt`,
`runs/ladder_cell89.txt`, `runs/ladder_pinwheel.txt`, `runs/ladder_tables.txt`,
`runs/ladder_t11.txt`, `runs/ladder_t11h.txt`, and the per-chain JSON
`runs/ladder_bandcut_scan_exact{110,132}.json`, `runs/ladder_h_bandcut_scan_exact{110,132}.json`.
Labels: **[proved]** = derivation done by hand here and re-checked in floating point by a script
here (exact cancellation to `< 1e-12`); **[measured]** = the value of an explicit configuration or
an LP on it; **[heuristic]**; **[guess]**.*

---

## 0. Verdict, up front

1. **The chain-counting lemma is true, and the proof is four lines — but the repo states it one
   notch too strongly.**  What `k > (T-1)^2` near-axis squares force (far squares arbitrary,
   `eps < 1/(T-1)`) is a chain of `T` in the **threshold order** `i <_x j  iff  x_j - x_i >= tau`,
   `tau = (1+delta)/W(eps)`, `W(eps) = cos eps + sin eps` — *not*, by this argument, a path in the
   x-DAG of `BANDCUT_K.md` §1.2.  The two coincide at `eps = 0` and differ at `O(eps)`.  The
   quantitative pay-off is
   `delta <= (T-1)(W(eps) - 1)/(T - 1 + 2 W(eps))  ~  (T-1) eps/(T+1)`, which is `<= 0` **only** at
   `eps = 0`.  §1.  **[proved]**
2. **Chain existence is not what `eps` buys.**  Run the same counting on the *x-DAG* itself: the
   only thing the tilt bound is needed for is **acyclicity** (so that "longest path ending here"
   is defined); the counting `k <= h_x h_y <= (T-1)^2` then needs nothing else.  So among any
   `k > (T-1)^2` squares whose two DAGs are acyclic a wall-to-wall chain of `T` exists, and with
   `n = T^2 - T > (T-1)^2` for every `T >= 2` **chain existence is free at every `T`** (the order
   form unconditionally, the DAG form wherever acyclicity holds, for which the same tilt bound is
   sufficient) — which is what `T11_CHAINS.md` §0(1) and `t3-chain.md` §0 measure at `T = 11, 12`
   and `T = 3` (`nochain = 0` at all 720 samples; five tight chains of eleven in Cantrell's
   packing, eight of twelve).  What
   `eps` buys is that the chain is **tight**: `h_x <= T`, so the `T` squares span the container.
   §1.4.  **[proved]**
3. **No counting statement forces a chain *through* tilted squares.**  Counting is insensitive to
   the tilts; what dies is the conclusion, not the hypothesis.  With the chain's members at tilt
   `<= eps` the forced bound is `(T-1)(W(eps)-1)/(T-1+2W(eps))`, which is `0.414` at `eps = 45°`
   for large `T` — a positive number, hence vacuous.  (And the `x`/`y` typing of normals, on which
   the *DAG* form rests, is degenerate at exactly `45°`: `t3-chain.md` §0's erratum.)  §1.5.
   **[proved]**
4. **The three-level ladder has a closed form and it reproduces every measured number.**  For the
   mixed ladder (crossbar of `T` with two adjacent tilted members, `L` transverse links, one rung)

       H1(T, A, B, l1, l2, eps)  =  -N/D ,
       N  =  (L + 2 - T) eps  +  (A - T + 3) eps^2  +  (2T - 2L - 5/2) eps^3/3  +  O(eps^4),
       D  =  (T + 1) + (L + 2) eps + O(eps^2),        L = l1 + l2,  A + B = T - 3,

   so **`A - T + 3 = -B`** and the whole lemma is again two integers: `L` against `T - 2`
   (first order) and `B` against `0` (second order).  At `L = T-2`, `B = 0` the first two orders
   vanish and

       H1(T, eps)  =  - eps^3 / (2(T+1))  +  O(eps^4) .

   At `T = 4` that is `-0.1 eps^3` against `runs/t4_cycles_cell89.txt`'s `-0.0989 eps^3` and
   `BANDCUT_K.md`'s measured chain-free maximum `-0.1004 eps^3`.  Without the rung the same ladder
   is `+eps^2/(T+1)` — `1/(T+1)` being exactly the far-field constant of `BANDCUT_K.md` §0(2), now
   proved rather than fitted.  §2.  **[proved]**
5. **The ladder needs the tilted block one square in from the wall.**  If the tilted pair sits
   *against* the main wall (`A = 0`, which is the `T = 11`/`T = 12` band shape) the second-order
   coefficient is `3 - T` instead of `0` and the ladder is `+(T-3) eps^2/(T+1) > 0` for every
   `T > 3`, at every `L`.  §2.4.  **[proved]**
6. **At `T = 11` and `T = 12` the clause that fails is leg existence — by a factor of two — and
   behind it the wall-block geometry.**  Chain existence holds (five tight chains of eleven, eight
   of twelve); the closed-form inequality holds (it is an identity); the common-normal hypothesis
   fails on the crossbar (kinks of `24-30°`) *and* on the legs.  Quantitatively: the band's kinks
   (`24-30°`) make the wall-block ladder need `L = 13-14` transverse links at `T = 11` and `14-16`
   at `T = 12`, against `T - 2 = 9, 10` for a clean ladder, and the packing offers `8-15`
   counted generously (both directions, any one normal) — and **exactly one link, in each
   direction, at every band kink**, which is the quantification of `T11_CHAINS.md` §3's
   unmeasured "unpinned".  Independently, the best H-shaped dual on the tight
   `x`-chain, *given every transverse pair row in the packing*, is `+6.9e-3` — 25 times `delta`
   itself.  §3.  **[measured]**
7. **`L*(T,0) = T - 2` is the same integer in both ladders**, and the number of levels needed to
   reach `<= 0` is **three at every `T`**; `T` enters only as (i) the number `T-2` of transverse
   links that must exist and (ii) the constant `1/(2(T+1))`.  §4.  **[proved]**

---

## 1. The chain-counting lemma

### 1.1 Setting, stated for Lean

Fix an integer `T >= 2`.  All angles are in radians; `W(x) = |cos x| + |sin x|`;
`P(th) = W(th)/2` (the half-width of a unit square of angle `th` in the `x` or `y` direction);
`m(i,j) = 1/2 + W(th_j - th_i)/2`; `tilt(th) = dist(th, (pi/2) Z) in [0, pi/4]`.

> **Hypotheses.**  Let `Q_1, ..., Q_n` be closed unit squares, `Q_i` with centre `c_i = (x_i, y_i)`
> and angle `th_i`, and let `delta in R`.
>
> * **(H1)** `delta > -1`.
> * **(H2)** *(containment, with margin)*  for all `i`:
>   `P(th_i) + delta <= x_i <= T - P(th_i) - delta` and the same for `y_i`.
> * **(H3)** *(packing, closed semantics, with margin)*  for all `i != j` there is an edge normal
>   `n` of `Q_i` or of `Q_j` (a unit vector at angle `th_i + k pi/2` or `th_j + k pi/2`) with
>   `|n . (c_j - c_i)| >= m(i,j) + delta`.
> * **(H4)** `K` is a set of `k` indices with `tilt(th_i) <= eps` for `i in K`; the indices outside
>   `K` are unconstrained.
> * **(H5)** `0 <= eps < pi/4`.
> * **(H6)** `k > (T-1)^2`.
>
> Write `tau = (1 + delta)/W(eps) > 0`, and for `i, j in K`
> `i <_x j  :<=>  x_j - x_i >= tau`, `i <_y j  :<=>  y_j - y_i >= tau`.

Two remarks on the hypotheses, both load-bearing.

* **(H3) at `delta = 0` is exactly "pairwise disjoint interiors, closed squares".**  Two convex
  polygons have disjoint interiors iff some edge normal of one of them separates them, and for
  unit squares the separation threshold along an edge normal of `Q_i` is
  `1/2 + W(th_j - th_i)/2 = m(i,j)`: the half-width of `Q_i` along its own normal is `1/2`, that of
  `Q_j` is `W(th_j - th_i)/2`.  Touching (`gap = 0`) is allowed and the inequality is non-strict.
  This is `notes/s13-casefree.md` §1 / `notes/branch-semantics.md` semantics, and it is why the
  lemma must be read **at the configuration's own margin** and not at gap `>= 0`
  (`BANDCUT_K.md` §0(4)).
* **(H2) is what "wall-to-wall" means.**  The two wall rows are available for *every* square, so a
  chain is "wall to wall" by fiat — it is closed by the wall rows of its two end squares, whether
  or not those squares touch the walls.  The resulting inequality is an equality only if they do.
  Nothing in the lemma requires the ends to touch, and nothing should: at `delta > 0` no square
  touches.

### 1.2 The lemma  **[proved]**

> **Lemma CC (order form).**  Under (H1)–(H6) there are `T` indices `i_1, ..., i_T in K` with
> either `i_1 <_x i_2 <_x ... <_x i_T` or `i_1 <_y ... <_y i_T`.
>
> **Corollary CC1 (the inequality it supports).**  With that chain, say in `x`,
>
>     (T-1) tau  <=  x_{i_T} - x_{i_1}  <=  T - P(th_{i_1}) - P(th_{i_T}) - 2 delta ,
>
> whence
>
>     delta  <=  [ T - P_1 - P_T - (T-1)/W(eps) ] / [ (T-1)/W(eps) + 2 ]
>            <=  (T-1) (W(eps) - 1) / (T - 1 + 2 W(eps)) .
>
> In particular **`eps = 0` gives `delta <= 0`**, with equality iff the `T` squares are
> axis-parallel, abut one another and touch both walls.
>
> **Corollary CC2 (heights).**  If moreover `delta >= 0` and `W(eps) < T/(T-1)` — for which
> `eps < 1/(T-1)` is sufficient, and `eps < eps_max(T) := arcsin(T/(sqrt2 (T-1))) - pi/4` is the
> sharp form of the same condition (`W(eps) = sqrt2 sin(eps + pi/4)`) — then every `<_x`-chain and every `<_y`-chain in `K` has at most `T`
> members, so the chain of Lemma CC is a *maximum* one and the `k` squares of `K` carry level maps
> into `{1..T}^2`.

*Proof of Lemma CC.*  Four steps.

1. **`<_x` and `<_y` are strict partial orders.**  Irreflexive and antisymmetric because
   `tau > 0` (H1); transitive because `x_j - x_i >= tau` and `x_l - x_j >= tau` give
   `x_l - x_i >= 2 tau >= tau`.  (This is the whole content of "acyclic and transitive enough for
   Mirsky": a threshold on a *global* coordinate is automatically both.  It is what the
   normal-classification route of `proof-architecture.md` §2 fails to give, as
   `proof-architecture-review.md` §5 says.)
2. **Comparability.**  Let `i != j in K` and let `n` be the normal of (H3), an edge normal of a
   square of tilt `<= eps`, so `|n_x| + |n_y| = W(tilt) <= W(eps)` (`W` is increasing on
   `[0, pi/4]`).  Then
   `1 + delta <= m(i,j) + delta <= |n . (c_j - c_i)| <= (|n_x| + |n_y|) max(|Dx|, |Dy|)`,
   using `m(i,j) >= 1` (because `W >= 1`).  Hence `max(|Dx|, |Dy|) >= (1+delta)/W(eps) = tau`, i.e.
   every pair of `K` is comparable in `<_x` or in `<_y`.
3. **Counting.**  Suppose no `<_x`-chain and no `<_y`-chain has `T` members.  Let `L(i)` be the
   number of members of a longest `<_x`-chain ending at `i`, and `M(i)` the same for `<_y`; by
   assumption `L, M : K -> {1, ..., T-1}`.  If `L(i) = L(j)` with `i != j` then `i, j` are
   `<_x`-incomparable (a relation would raise one of the values), so by step 2 they are
   `<_y`-comparable, so `M(i) != M(j)`.  Therefore `i |-> (L(i), M(i))` is injective and
   `k <= (T-1)^2`, contradicting (H6).  ∎
4. *(Corollary CC1)* is the display; *(CC2)* is `(h-1) tau <= T - 1 - 2 delta` with
   `P >= 1/2`, i.e. `h <= 1 + (T-1-2 delta) W(eps)/(1+delta) <= 1 + (T-1) W(eps) < T+1` when
   `W(eps) < T/(T-1)`, and `W(eps) = sqrt2 sin(eps + pi/4)`.  `W(eps) <= 1 + eps` gives the
   sufficient form `eps < 1/(T-1)`.  ∎

`eps_max(T)` in degrees: `45, 45, 25.529, 17.114, 13.052, 6.061, 5.479, 3.703` for
`T = 2, 3, 4, 5, 6, 11, 12, 17` (and `1/(T-1)` in degrees is `57.3, 28.6, 19.1, 14.3, 11.5, 5.73,
5.21, 3.58`); `eps_max(T) - 1/(T-1) = O(1/T^2)`, so the review's constant is right to first order
and slightly pessimistic.  **[proved]**

### 1.3 Sharpness: the threshold `(T-1)^2` cannot be lowered, at any `T`  **[proved]**

Put `(T-1)^2` axis-parallel unit squares at the centres `(s(a + 1/2), s(b + 1/2))`,
`0 <= a, b <= T-2`, `s = T/(T-1)`.  They are pairwise disjoint and inside `[0,T]^2` with margin

    delta  =  min( s/2 - 1/2 , (s-1)/2 )  =  1/(2(T-1))  >  0

(`1/4` at `T = 3`, `1/6` at `T = 4`, `1/20` at `T = 11`), they are exactly axis-parallel (so any `eps`), and they
occupy `T-1` distinct `x`-levels and `T-1` distinct `y`-levels, so the longest chain — in `<_x`, in
`<_y`, in the x-DAG and in the y-DAG alike — has `T-1` members.  So at `k = (T-1)^2` the
conclusion is false with room to spare, for every `T`.  (This is the `T = 4` witness of
`proof-architecture-review.md` §1 with a better constant: `{0.6, 2.0, 3.4}` has `delta = 0.1`,
`{2/3, 2, 10/3}` has `delta = 1/6`.)  A second, less artificial `T = 3` witness — a real packing, with a far
square present — is the repo's own filter case: the optimal five-square packing has
`k = 4 = (T-1)^2` near-axis squares, `h_x = h_y = 2`, one square at `45°`, and
`delta = +0.0540971 > 0` after the dilation of S0.

**What breaks for `k <= (T-1)^2` is therefore not slack in the argument but the statement.**  The
counting is an exact pigeonhole on a product of two level sets, and the product is full.  Any lemma
covering `T <= k <= (T-1)^2` must be geometric — it must use `n = T^2 - T`, i.e. the far squares —
and that is the hole of `proof-architecture.md` §0a item 2, unchanged by anything here.

### 1.4 The same counting on the DAG, and what `eps` is really for  **[proved]**

Define the x-DAG exactly as `BANDCUT_K.md` §1.2: `i -> j` iff some **x-type** edge normal `n` of
`Q_i` or `Q_j` (x-type = nearer `±e_x` than `±e_y`; unambiguous for `tilt < pi/4`) has
`n . (c_j - c_i) >= m(i,j) + delta`.  Then:

* **(a) Comparability is free.**  For `tilt < pi/4` every edge normal is x-type or y-type, so by
  (H3) every pair is an edge of the x-DAG or of the y-DAG.
* **(b) Acyclicity is exactly what the tilt bound buys.**  If `i -> j` in the x-DAG and both tilts
  are `<= eps`, then with `n = (cos ph, sin ph)`, `|ph| <= eps`,
  `cos(ph) Dx = n.D - sin(ph) Dy >= 1 + delta - eps |Dy| >= 1 + delta - eps (T-1-2 delta)`
  and `cos(ph) > 0`, so `Dx > 0` as soon as `eps < (1+delta)/(T-1-2 delta)` — in particular whenever
  `eps < 1/(T-1)` and `delta >= 0`.  Then `x` strictly increases along every edge and the x-DAG is
  acyclic, so `L(i) = ` longest path ending at `i` is well defined.
* **(c) Counting.**  With (a) and (b), the proof of §1.2 step 3 applies verbatim to the
  reachability orders: `k <= h_x h_y`, so `k > (T-1)^2` forces a path on `T` vertices in the x-DAG
  or the y-DAG.
* **(d) But the bound it supports is a factor `T-1` worse.**  Summing the `T-1` link rows and
  projecting on the global axis (`BANDCUT_K.md` §1.2) leaves the rise:
  `(T-1)(1+delta) <= sum_s n_s.D_s <= (T - P_1 - P_T - 2 delta) + eps sum_s |Dy_s|`, and the only
  a-priori bound on the non-telescoping `sum_s |Dy_s|` is `(T-1)(T-1-2 delta)`, giving
  `delta <= eps (T-1)^2/(T+1)` against the order form's `eps (T-1)/(T+1)`.

So the two forms of the conclusion are genuinely different objects: **the order chain is the one
that carries a sharp inequality, the DAG chain is the one that always exists.**  An order edge need
not be a DAG edge (a pair with `Dx = 1.5` and `Dy = -(T-1)` at common tilt `eps ~ 1/(T-1)` is
separated only in `y`), and a DAG edge need not be an order edge (`Dx` can be as small as
`1 + delta - eps(T-1)`).  Whether the `T` squares of Lemma CC can always be chosen to form a DAG
path is **open**; it is measured true at `T = 3, 4` (`BANDCUT_K.md` §0(1): every `delta >= 0`
configuration found carries a wall-to-wall chain of `T` among its near-axis squares).

The gap between the two is *exactly* the quantity the ladder of §2 exists to control: the rise
term `eps sum_s |Dy_s|`.  The counting lemma and the ladder lemma are the same statement at two
resolutions.  **[heuristic]**

### 1.5 A chain through tilted squares: what a counting argument can and cannot do  **[proved]**

Counting never looks at the tilts; it looks at `tau`.  Run Lemma CC with `K` = *all* `n` squares
and `eps -> pi/4`: comparability holds at `tau = (1+delta)/sqrt2`, transitivity is free, and
`n = T^2 - T > (T-1)^2` for every `T >= 2`, so **an order chain of `T` exists among all the
squares, tilted or not**, unconditionally.  (The corresponding *DAG* chain also exists whenever the
two DAGs are acyclic, §1.4 (b)–(c); that a wall-to-wall DAG chain of `T` is present is what the
measurements report at `T = 3, 4, 11, 12`.)  What fails is the conclusion: Corollary CC1 returns

    delta  <=  (T-1)(W(eps) - 1)/(T - 1 + 2 W(eps)) ,

which at `eps = 45°` is `0.1716, 0.2132, 0.2426, 0.3229, 0.3520, 0.4027` for `T = 3, 4, 5, 11, 17, 100`
— positive, increasing in `T`, and tending to `sqrt2 - 1 = 0.414`.  So:

> **There is no counting statement that forces a *tight* chain through tilted squares.**  The
> tilts enter only through `W`, and `W(eps) = 1` iff `eps = 0`.  A chain whose members have tilts
> `t_1 .. t_T` costs `sum_s (W(t_s) - 1)/2 ~ sum |t_s|/2` in the numerator — the same first-order
> price per mismatched link as `bandcut-cost.md` §1 — and counting has no way to pay it.

Two further riders.  (i) At exactly `45°` the x/y typing of normals is degenerate and the DAG is
not well defined (the erratum of `t3-chain.md` §0 about `bandcut_k.chain_certificate`); Lemma CC in
order form does not care (it never types a normal), which is one more reason to prefer it.
(ii) Acyclicity (§1.4 (b)) is proved only for `eps < (1+delta)/(T-1-2delta)`; for larger tilts I
have neither a proof nor a counterexample, so the DAG counting at large tilt is **open**, and the
tight chains of `T11_CHAINS.md` §2 are found by search, not forced by counting.

### 1.6 The two filters

| filter | what the lemma gives | verdict |
|---|---|---|
| **`T = 17`, Cleemann's 272 in side `< 17`** (`ARCH_TLEDGER.md` §5.3: 199 axis-parallel within `3°`, 35 at `+arctan(8/15)`, 35 at `-arctan(8/15)`, 3 at `45°`) | `199 <= 256 = (T-1)^2`, so **(H6) fails and the lemma says nothing about this packing.**  What it does say: since `delta > 0`, at most `(T-1)^2 = 256` squares can be exactly axis-parallel, i.e. `>= 272 - 256 = 16` must be tilted — and `73` are.  At `eps > 0` it says `>= 16` squares have tilt `> eps_0`, where `eps_0` solves `16(W(eps_0)-1)/(16+2W(eps_0)) = delta` (Cleemann's exact side is unpublished, so `delta` is unknown; `ARCH_TLEDGER.md` §5.2). | **passes, silently.**  The lemma is *true* at `T = 17` and carries no `T`-dependence, confirming `proof-architecture-review.md` §4(b) against `proof-architecture.md` §3.3's guess (β). |
| **`n = T^2 - 4`, `T = 3`** (`s(5) = 2 + 1/sqrt2 < 3`) | The optimal packing has `k = 4 = (T-1)^2` near-axis squares and one at `45°`, so again **(H6) fails by exactly one**.  Had all five been near-axis, Corollary CC1 would give `delta <= 2(W(eps)-1)/(2+2W(eps))`, i.e. `delta <= 0` at `eps = 0`: five axis-parallel unit squares in `[0,3]^2` indeed have `delta <= 0`.  The real packing's `delta = +0.0540971` would need `eps >= 6.998°`; the fifth square is at `45°`. | **passes, and by the tilt**, exactly as `t3-chain.md` §6.3 reports for Lemma H.  Note `6 > 4` does force a chain of three as soon as five of the six squares are near-axis — which is the `5 -> 6` step. |

---

## 2. The three-level ladder in closed form

### 2.1 What is being generalised

`t3-chain.md` §1.3's Lemma H is the two-level object: a crossbar of `T` at weight `1`, `L = a + c`
transverse leg links at weight `tan t`, six wall rows, of which the two **closing main-direction
wall rows** `hi-ax(b_1)`, `lo-ax(d_c)` carry weight `sin t tan t`.  The ladder replaces a closing
wall row by a **rung**: one further main-direction link at weight `tan^2 t` onto a fresh square,
plus that square's own wall rows.  `T4_CYCLES.md` §3.1 reads exactly this off the `eps^3` cell —
levels at `w_0`, `w_0 tan eps`, `w_0 tan^2 eps`, ratios `0.017450` and `0.017456` against
`tan 1° = 0.0174551` — and §3.2 names it `H1`.

Two ladders are needed, because the `eps^3` cell is **not** at uniform tilt.

### 2.2 The mixed ladder (the `eps^3` cell), square by square  **[proved]**

> **Definition (the mixed ladder `H1m(T, A, B, l1, l2, t)`).**  A crossbar `s_0 -> ... -> s_{T-1}`
> in the x-DAG whose members are axis-parallel except `s_A` and `s_{A+1}`, which have tilt `t`;
> `A >= 1`, `B = T - 3 - A >= 0`.  Link normals: links `0..A-1` and `A+2..T-2` are axial
> (`(1,0)`), links `A` and `A+1` are `(cos t, sin t)`.  One **entry leg** of `l1` links hangs
> below `s_A`, its first link along `(-sin t, cos t)` (the transverse normal owned by the tilted
> square `s_A`) and the rest axial; one **exit leg** of `l2` axial links hangs above `s_{A+2}`.
> `L = l1 + l2`.  Wall rows: `lo-x(s_0)`, `hi-x(s_{T-1})`, `lo-y(` foot of the entry leg `)`,
> `hi-y(` foot of the exit leg `)`, and at the entry leg's **kink** (its first square, where the
> normal turns from tilted to axial) either a closing `hi-x` wall row or a **rung**: one axial
> main-direction link to a fresh square plus that square's `hi-x` wall row.

**The weights, by the same bookkeeping as `t3-chain.md` §1.3.**  Write `C = cos t`, `S = sin t`,
`u = C + S`, `g = m(0, t) = (1+u)/2`.  Normalise the crossbar's first link to weight `1`.

| rows | weight | `b` |
|---|---|---|
| `lo-x(s_0)` | `1` | `1/2` |
| crossbar links `0 .. A-2` (axial, `0 <-> 0`) | `1` | `1` |
| crossbar link `A-1` (axial normal, `0 <-> t`) | `1` | `g` |
| crossbar link `A` (tilted normal, `t <-> t`) | `C` | `1` |
| crossbar link `A+1` (tilted normal, `t <-> 0`) | `C` | `g` |
| crossbar links `A+2 .. T-2` (axial) | `C^2` | `1` |
| `hi-x(s_{T-1})` | `C^2` | `1/2 - T` |
| entry leg, first link (tilted transverse normal) | `S` | `g` |
| entry leg, links `2 .. l1` (axial) | `S C` | `1` |
| `lo-y(` entry foot `)` | `S C` | `1/2` |
| the **rung** at the entry kink | `S^2` | `1` |
| `hi-x(` rung end `)`  (or `hi-x(` kink `)` if no rung) | `S^2` | `1/2 - T` |
| exit leg, `l2` axial links | `C S` | `1` |
| `hi-y(` exit foot `)` | `C S` | `1/2 - T` |

*Derivation (one square at a time, as in `t3-chain.md` §1.3).*  At `s_0` the wall row and the first
link cancel in `x` and nothing else appears.  At `s_A` (tilt `t`) the incoming axial link
contributes `(1,0)`, the outgoing tilted link `-p(C,S)` and the leg `+q(-S,C)`: the `y` equation
gives `q = p tan t` and the `x` equation `1 = p(C + S^2/C) = p/C`, so `p = C` and `q = S` — *this
is where the leg's weight `w_0 tan t` comes from, and why a leg must hang exactly at the square
where the normal turns.*  At `s_{A+1}` both links carry the same normal and telescope.  At
`s_{A+2}` (tilt `0`) the incoming tilted link leaves `(C^2, CS)`: the `x` part is taken by the
outgoing axial link (or by `hi-x`) at weight `C^2`, the `y` part by the exit leg at weight `CS`.
At the entry leg's kink the tilted-normal link leaves `(S^2, -SC)`: the `y` part is taken by the
next (axial) leg link at weight `SC`, and the `x` part `S^2` is the residual that the closing wall
— or the rung — must absorb.  At the rung's far square, `hi-x` at weight `S^2`.  Every centre
coordinate cancels; `search/ladder_hform.py check` re-derives the vector and reports the maximum
uncancelled coefficient as `< 1e-12` for `T in {4,5,6,11}`, every `A`, `l1, l2 in {1,2,T-3}`,
`t in {0.5, 1, 5, 20}°`, with and without the rung, and all weights `>= 0`.  ∎

> **Closed form.**  `delta <= H1m = -N/D` with
>
>     D  =  (A+1) + 2C + (B+1)C^2 + S + (l1+l2+1) S C + rho S^2 ,          rho = 2 with the rung, 1 without
>     N  =  1/2 + (A-1) + g + C + C g + B C^2 + C^2 (1/2 - T)
>           + S g + (l1-1) S C + S C/2 + [rho - 1] S^2 + S^2 (1/2 - T)
>           + l2 C S + C S (1/2 - T) .
>
> Series (`A + B = T - 3`, `L = l1 + l2`):
>
>     N  =  (L + 2 - T) t  -  B t^2  +  (2T - 2L - 5/2) t^3/3  +  (T/3 - A/3 - 5/4) t^4 + ... ,
>     D  =  (T + 1) + (L + 2) t - B t^2 + ... ,
>
> and **without** the rung `N`'s `t^2` coefficient is `-B - 1` instead of `-B`.

**Consequences, in order.**  `N_1 = L + 2 - T`: the ladder kills at first order iff `L > T - 2`,
is degenerate at `L = T - 2 = L*(T, 0)` — the same integer as Lemma H's threshold, and for the same
reason (each leg link is worth one unit of `m` at weight `~t`, each transverse wall costs `T` at
weight `~t`, and the two tilt-mismatched crossbar links pay `2 x t/2`).  `N_2 = -B`: at `L = T-2`
the ladder is `+B t^2/(T+1) > 0` unless `B = 0`, i.e. unless the tilted pair is adjacent to the
crossbar's far end.  With `L = T-2` and `B = 0`, `N = t^3/2 - t^4/4 + ...` and

    **H1m(T, T-3, 0, 1, T-3, t)  =  - t^3/(2(T+1))  +  O(t^4)** ,

and numerically `H1m <= 0` for every `t in (0°, 45°]` and every `T in {4, 5, 11, 12}` tested, not
only in the small-angle limit (`runs/ladder_tables.txt`).

### 2.3 Verification against `runs/t4_cycles_cell89.txt`  **[proved] + [measured]**

`T = 4`, `A = 1`, `B = 0`, `l1 = l2 = 1`, one rung — the exact shape `T4_CYCLES.md` §3.1 reads off
the dual (`lo-x(1), pair(1,4), pair(4,11), pair(7,11), hi-x(7)` at `0.19725`; `pair(3,4)`,
`lo-y(3)`, `pair(7,8)`, `hi-y(8)` at `0.00344`; `pair(3,6)`, `hi-x(6)` at `0.00006008`).  The
closed-form weights are `w_0 (1, 1, C, C, C^2; S, SC, CS, CS; S^2, S^2)`, i.e. at `eps = 1°`
`w_0 = 0.19724638` and

| row | closed form | `runs/t4_cycles_cell89.txt` (`k = 9`, `eps = 1°`) |
|---|---|---|
| `lo-x(s_0)`, crossbar link `0` | `w_0 = 0.19724638` | `0.19724638` |
| crossbar links `1, 2` | `w_0 C = 0.19721634` | `0.19721634` |
| `hi-x(s_3)` | `w_0 C^2 = 0.19718630` | `0.19718630` |
| entry leg link (tilted normal) | `w_0 S = 0.00344224` | `0.00344224` |
| `lo-y`, exit leg link, `hi-y` | `w_0 SC = 0.00344172 … 0.00344190` | `0.00344172, 0.00344190, 0.00344190` |
| rung, `hi-x(` rung end `)` | `w_0 S^2 = 0.00006008` | `0.00006008` |

and the values:

| `eps` | `CW` closed | `CW` LP | `H` closed | `H` LP | `H1` closed | `H1` LP | `delta*` |
|---|---|---|---|---|---|---|---|
| `1°` | `+6.873e-03` | `+6.872e-03` | `+5.956e-05` | `+5.955e-05` | `-5.197e-07` | `-5.285e-07` | `-5.260e-07` |
| `5°` | `+3.235e-02` | `+3.226e-02` | `+1.363e-03` | `+1.363e-03` | `-5.931e-05` | `-5.805e-05` | `-6.189e-05` |
| `10°` | `+6.027e-02` | `+5.956e-02` | `+4.905e-03` | `+4.903e-03` | `-4.233e-04` | `-4.062e-04` | `-4.596e-04` |

divided by the predicted power of `eps` (radians):

| `eps` | `CW/eps` closed / LP | `H/eps^2` closed / LP | `H1/eps^3` closed / LP | `delta*/eps^3` |
|---|---|---|---|---|
| `1°` | `0.3938 / 0.3938` | `0.1955 / 0.1955` | `-0.0978 / -0.0994` | `-0.0989` |
| `5°` | `0.3707 / 0.3697` | `0.1790 / 0.1790` | `-0.0892 / -0.0874` | `-0.0931` |
| `10°` | `0.3453 / 0.3413` | `0.1610 / 0.1610` | `-0.0796 / -0.0764` | `-0.0864` |
| `-> 0` | `(T-2)/(T+1) = 0.4` | `1/(T+1) = 0.2` | `-1/(2(T+1)) = -0.1` | (measured `-0.1004`, `BANDCUT_K.md`) |

So all three levels of `T4_CYCLES.md` §3.2's table are closed-form, and the three constants
`(T-2)/(T+1)`, `1/(T+1)`, `-1/(2(T+1))` are exact.  The middle one settles the observation of
`T4_CYCLES.md` §3.2(2) — "`0.1955 -> 1/5 = 1/(T+1)`, the far-field constant of `BANDCUT_K.md`
§0(2)" — as an identity, not a coincidence: the two-level `H` on a crossbar with two tilted members
is `+eps^2/(T+1)`, and `1/(T+1)` is the number of `delta`-carrying rows on a wall-to-wall chain.
`CW` is the leg-free ladder (both interface residuals taken by transverse wall rows); its closed
form is in `ladder_hform.CWm` and its limit is `(T-2)/(T+1)`.

*(Where the LP beats the closed form — `H1` at `1°`, `-5.285e-7` against `-5.197e-7` — it is the
same support with re-optimised weights, and where it beats `delta*` it is inside the recorded
tolerances of that run; the closed form is the valid bound in every row.)*

### 2.4 The wall-block variant, and why the ladder needs one square of clearance  **[proved]**

If the tilted pair is `s_0, s_1`, so that a tilted square sits against the main wall (`A = 0`; this
is the band shape of `T = 11`, `T = 12` and of Cleemann), the same bookkeeping gives `lo-x(s_0)` at
weight `1` with `b = u/2` (a *gain* of `(u-1)/2`, since a tilted square against a wall is wider)
but loses the mixed axial crossbar link `A-1` (a loss of the same `(u-1)/2`), and leaves `T-3`
axial links after the block.  Series:

    N  =  (L + 2 - T) t  +  (3 - T) t^2  +  (2T - 2L - 5/2) t^3/3 + ... ,   D = (T+1) + (L+2) t + ... ,

i.e. `N_2 = 3 - T` (the `B = T-3` case of §2.2), so at `L = T - 2`

    H1mw(T, eps)  =  + (T-3) eps^2/(T+1)  +  O(eps^3)   >  0   for every  T > 3 .

`H1mw` is in `ladder_hform.H1mw`, checked the same way.  At `T = 11` and a band kink of `29.5°` it
takes `L = 14` transverse links to bring `H1mw` to `0` (`+0.1277, +0.0960, +0.0108, -0.0147` at
`L = 9, 10, 13, 14`), against `L = T - 2 = 9` for the clearance-one ladder; at `T = 12` and `28°`,
`L = 15`.  **This is one of the two clauses `T = 11` violates** (§3.2), and the one that
raises its `L need`; it is invisible in the two-level Lemma H, whose statement has no notion of
*where* along the crossbar the tilt sits.

### 2.5 The uniform-tilt ladder, and the two pinwheel duals  **[proved]**

The same rung construction on Lemma H itself (`ladder_hform.uladder` / `H1u(T, L, t, r)`,
`r in {0,1,2}` rungs): at the leg's foot the residual `(S^2/C, -S)` is taken, instead of by the
closing wall `hi-ax` at weight `S tan t`, by a main-direction link at weight `nu = tan^2 t` plus
`hi-ax` and `hi-tr` at weights `nu C`, `nu S` at its far end, with the transverse wall weight
raised from `S` to `S + nu S`.  `r = 0` reproduces `t3_chain_hform.H` exactly.

The rungs help but do **not** move Lemma H's threshold off `L*(T, t)`: at `L = T-2` and `t = 1°`,
`H1u` is `+1.86e-4 / +1.13e-4 / +4.08e-5` for `r = 0, 1, 2` at `T = 3` and
`+2.38e-4 / +1.81e-4 / +1.24e-4` at `T = 4` — positive throughout.  **At uniform tilt the ladder
is `O(t)` positive at `L = T-2` and only the mixed ladder degenerates to `O(t^3)`**; the `eps^3`
cell is a mixed-tilt phenomenon, as `T4_CYCLES.md` §3 says.

What the uniform ladder *does* reproduce exactly is the **weight structure** of the two exact
pinwheel duals of `S6_LOCAL.md` §3.  In hard containment (the normalisation of `s6exact.py`: the
*pair* weights sum to `1`) the ladder's crossbar weight is

    w_0  =  1 / [ (T-1) + L tan t + 2 tan^2 t ]  =  1/(T-1)  -  L t/(T-1)^2  +  O(t^2) ,

legs `w_0 tan t`, rungs `w_0 tan^2 t`.  With `L = T - 2` and `r = 2`:

| | `S6_LOCAL.md` §3 | ladder, closed form | ladder at `t = 1°` |
|---|---|---|---|
| `n = 6`, `T = 3`, `L = 1` | `1/2 - t/4`, transverse `t/2`, two links `t^2/2` | `1/2 - t/4`, `w_0 tan t`, `w_0 tan^2 t` | `0.495524`, `0.017455 w_0`, `0.000305 w_0` |
| `n = 12`, `T = 4`, `L = 2` | `1/3 - 2t/9`, transverse `t/3`, links `t^2/3` | `1/3 - 2t/9`, `w_0 tan t`, `w_0 tan^2 t` | `0.329433`, `0.017455 w_0`, `0.000305 w_0` |

— both first-order corrections (`-1/4` and `-2/9`) come out of `-L/(T-1)^2` with `L = T-2`, which
is a non-trivial check that the pinwheel duals *are* three-level ladders with `L*(T,0)` legs.
Their **values** are not reproduced (`-3(u-1)^2/(u^2+3) = -0.75 t^2` at `n = 6` against the generic
ladder's `+0.26 t^2`): the pinwheel dual closes its two rungs with four wall rows in all, where the
generic ladder needs eight.  So the ladder family of §2.2 is one family among several with the same
level structure, and it is not the cheapest at uniform tilt.  **[measured]** for that gap.

### 2.6 The existence clause, in one sentence

> **Existence clause.**  In the separation graph at the configuration's own margin there is a
> wall-to-wall chain of `T` whose link normals take exactly two values — the global axis and one
> tilt `t` — changing exactly twice, at a tilted square `s_A` with `A >= 1` and at the square
> `s_{A+2} = s_{T-1}` next to the far wall (equivalently: `B = 0`, the tilted block is adjacent to
> the crossbar's far end and does not touch the near wall), and transverse chains with a common
> normal hang at those two squares, in the directions the two kinks demand, carrying
> `l1 + l2 >= T - 2` links in total, with one further main-direction link available at the kink of
> the tilted-normal leg.

---

## 3. `T = 11` and `T = 12`: which clause fails

*Input: `runs/bandcut_scan_exact110.json` (Cantrell, `delta = +2.748946598e-04`, exact rational)
and `runs/bandcut_scan_exact132.json` (`T = 12`, `delta = +6.81e-04`).  Instruments
`ladder_hform.t11` (the legs available at each kink of each tight chain) and `ladder_hform.t11h`
(the best restricted dual on each tight chain, the same instrument `t3_chain.py` uses at `T = 3`).
Every number here is **[measured]** on one configuration.*

### 3.1 The table

Per tight wall-to-wall chain (slack `< 1e-6`, deduplicated by its tilted sub-chain).  **kinks** =
squares where the link normal changes, counting the two walls as axial normals — each one needs a
leg.  **legs** = the longest transverse chain hanging at that kink whose links all carry a *single*
normal, that normal being an edge normal of the kink square itself (`up/down` = towards the hi/lo
transverse wall), with the crossbar's own squares excluded; in brackets, the same allowing any one
of the normal classes present on the chain (over-generous: it ignores which normal cancels the
kink and which direction the kink's sign demands).  **L need** = the least `L` making the
wall-block closed form `H1mw(T, 1, L-1, maxkink)` non-positive (§2.4).  **H1mw** = that closed form
at the `L` the packing actually offers.  **H(chain)** = the best restricted dual on that chain
given *every* wall row and *every* transverse pair row of the packing (`t3_chain.h_support` mode
`H`, the same instrument as `t3-chain.md` §2) — an upper bound on what any leg structure whatever
can achieve there.

| | chain (tilts along it) | kinks | legs at each kink, own normal (best single normal) | L avail | L need | `H1mw` at L avail | `CW` | `H`(chain) | `delta` |
|---|---|---|---|---|---|---|---|---|---|
| `T = 11` x | `-30 -29 -29 -25` then 7 x `0` | 3 (`29.5°` max) | `1/1` (5/3); `1/2` (5/3); `6/3` (6/3) | 9 | 14 | `+0.1275` | `+0.500` | `+6.91e-03` | `+2.75e-04` |
| `T = 11` x | `-25 -25 -25` then 8 x `0` | 2 (`24.9°`) | `1/1` (6/2); `8/1` (8/1) | 9 | 13 | `+0.0935` | `+0.333` | `+6.91e-03` | |
| `T = 11` y | 8 x `0` then `+24 +24 +24` | 2 (`23.9°`) | `1/8` (1/8); `1/1` (2/6) | 9 | 13 | `+0.0870` | `+0.285` | `+1.60e-03` | |
| `T = 11` y | 7 x `0` then `+27 +28 +28 +28` | 3 (`28.1°`) | `3/6` (3/6); `1/1` (3/5); `1/2` (4/4) | 9 | 14 | `+0.1169` | `+0.318` | `+2.68e-03` | |
| `T = 11` y | 9 x `0` then `+28 +28` | 2 (`28.3°`) | `0/9` (0/9); `1/1` (1/7) | 10 | 14 | `+0.0876` | `+0.200` | `+5.20e-03` | |
| `T = 12` x | `-24 -24 -24 -24` then 8 x `0` | 2 (`23.7°`) | `1/1` (6/3); `1/1` (7/2) | 2 | 14 | `+0.3343` | `+0.268` | `+6.09e-03` | `+6.81e-04` |
| `T = 12` x | `-24 -24 -24` then 9 x `0` | 2 (`23.7°`) | `1/1` (7/1); `9/1` (9/2) | 10 | 14 | `+0.0890` | `+0.300` | `+5.39e-03` | |
| `T = 12` x | `-31 -31` then 10 x `0` | 3 (`30.7°`) | `1/0` (8/1); `1/1` (9/1); `10/0` (10/0) | 12 | 16 | `+0.0842` | `+0.182` | `+6.09e-03` | |
| `T = 12` x (the elbow) | `0 0 -26 -25 +8 +24 +22 +7` then 4 x `0` | 7 (`33.5°`) | `4/6`; `1/2`; `1/2`; `1/1`; `1/1`; `1/1`; `1/1` | 14 | 16 | `+0.0474` | `+0.617` | `+2.06e-03` | |
| `T = 12` x | `-25 -25 -25 -25` then 8 x `0` | 2 (`25.4°`) | `1/1` (6/3); `7/3` (7/3) | 8 | 15 | `+0.1574` | `+0.444` | `+6.09e-03` | |
| `T = 12` y | 10 x `0` then `+28 +28` | 2 (`28.3°`) | `0/10` (0/10); `1/1` (1/8) | 11 | 15 | `+0.0952` | `+0.182` | `+3.06e-03` | |
| `T = 12` y | 9 x `0` then `+25 +25 +25` | 2 (`25.0°`) | `1/9` (1/9); `1/1` (2/7) | 10 | 15 | `+0.0988` | `+0.300` | `+1.84e-03` | |
| `T = 12` y | 8 x `0` then `+24 +24 +24 +24` | 4 (`24.0°`) | `3/7`; `1/1`; `1/4`; `2/3` | 15 | 14 | **`-0.0221`** | `+0.289` | `+2.24e-03` | |

The one entry where the proxy closed form goes negative — the last `T = 12` row — is where the
proxy is *invalid*: that chain has **four** kinks, so a ladder on it needs four legs and four
transverse wall rows, and `H1mw` charges for two.  Its honest reading is the `H` column, `+2.2e-3`.
The same caveat applies, less severely, to every row with `nchg > 2`.

*(`C`, the crossbar's own link rows plus its two end wall rows and nothing else, is `+inf` on every
one of these chains: with kinked normals that support admits no non-negative cancelling combination
at all.  It is `+inf` in the `T = 4` `eps^3` cell too, so this is not distinctive — it is what the
legs are for.)*

### 3.2 Which clause fails

| clause of the ladder lemma (§2.6) | at `T = 11`, `T = 12` | evidence |
|---|---|---|
| **the inequality** `H1 <= 0` at `L >= T-2` | **holds.**  It is an identity in `(T, A, B, l1, l2, t)`, verified to `1e-12`, and for the clearance-one ladder it is non-positive at *every* tilt up to `45°`. | §2.2, `runs/ladder_check.txt` |
| **chain existence** | **holds** — five tight wall-to-wall chains of eleven, eight of twelve.  (Measured, not forced: §1.4's counting licenses a DAG chain only where the DAG is acyclic, and its sufficient condition `eps < 1/(T-1) = 5.7°` is far below these band tilts.  The *order* chain of §1.2 is forced unconditionally, but it is the weak form.) | `T11_CHAINS.md` §2, §4 |
| **the common-normal hypothesis on the crossbar** | **fails**, but that is not fatal by itself: the ladder is *designed* for a crossbar with two normals and two kinks.  What is fatal is the *size* of the kink: `24-30°` instead of `eps`, which multiplies every leg weight by `tan(kink) ~ 0.5` and hence multiplies the cost `T x (transverse wall weight)` by the same. | `taus` column, `runs/ladder_*.json` |
| **the wall-block clause** (`B = 0`, block not against the wall) | **fails.**  In all thirteen tight chains the tilted block runs into the wall, so the second-order coefficient is `3 - T` and not `0`, which is what raises `L need` from `T-2 = 9, 10` to `13-16`. | §2.4 |
| **leg existence** | **fails, and it is the binding one.**  At every *band* kink the longest transverse chain with a common normal is **exactly one link** in each direction — which is the measurement `T11_CHAINS.md` §3 left undone ("transverse-chain search through the band squares to quantify 'unpinned'").  Total available `2-15` (counted generously, and only `8-12` on the nine chains with `nchg <= 3`) against `13-16` needed. | the table |

**And a check that does not depend on the ladder's shape at all:** the best dual supported on the
chain, every wall row and *every* transverse pair row of the packing is `+1.6e-3` to `+6.9e-3`
at `T = 11` and `+1.8e-3` to `+6.1e-3` at `T = 12`, i.e. `3` to `25` times `delta` and on the wrong
side of zero.  No transverse structure that exists
in Cantrell's packing closes any of its tight chains, whatever ladder one assembles from it.
(This is forced, of course — `delta > 0` means every valid dual is `> delta` — but the *factor*
`6-25` says the failure is not marginal, in contrast with `T = 4`, where the closed form lands within `1.2 %`
of `delta*` and the LP on the same support within `0.5 %`.)

**[measured]**, on two configurations.



---

## 4. `L*(T,0) = T-2`, and how the `eps`-orders depend on `T`

Collecting §2:

| object | `N_1` | `N_2` | leading value at `L = L*(T,0) = T-2` |
|---|---|---|---|
| bare crossbar + walls (`CW`) | `2 - T` | — | `+(T-2) eps/(T+1)` |
| two-level `H` (mixed, `B = 0`) | `L + 2 - T` | `-B - 1 = -1` | `+eps^2/(T+1)` |
| three-level `H1` (mixed, `B = 0`) | `L + 2 - T` | `-B = 0` | `-eps^3/(2(T+1))` |
| three-level `H1` (mixed, `B > 0`) | `L + 2 - T` | `-B` | `+B eps^2/(T+1)` |
| wall-block `H1mw` | `L + 2 - T` | `3 - T` | `+(T-3) eps^2/(T+1)` |
| uniform-tilt Lemma H / `H1u` | — | — | `> 0`; needs `L >= L*(T,t) > T-2` |

**Answer to the brief's question 4.**  In the closed form the number of levels needed to reach
`<= 0` is **three at every `T`** — one level per order of `eps`, `O(eps) -> O(eps^2) -> O(eps^3)`,
exactly as `T4_CYCLES.md` §3.2(3) observes at `T = 4` — and the depth does **not** grow with `T`.
What depends on `T` is entirely *which rows must exist*: `L >= T - 2` transverse links (so the
certificate spans `2T - 2` of the `T^2 - T` squares, `t3-chain.md` §6.2 item 3) and `B = 0` (the
tilted block adjacent to an end).  The only place `T` appears in the value is the denominator
`T + 1` — the number of `delta`-carrying rows on the crossbar.  This is the same conclusion as
`t3-chain.md` §6: **the inequality is `T`-uniform and cannot be the thing that fails at
`T = 11`; the existence clauses are.**  (`tasks/break-h1` measures the depth at `T = 5`; the
prediction here is that the depth is again three and that the `eps^3` constant is `-1/12`.)

`L*(T, 0) = T - 2` itself now has three independent derivations that agree: the removable
singularity of `t3-chain.md` §1.3's corollary; the first-order coefficient `L + 2 - T` of the mixed
ladder; and the ledger "each leg link gains `1` and each pair of transverse walls costs `T`, at the
same weight `~ w_0 tan t`, with the two tilt-mismatched crossbar links contributing `2 x 1/2`".

---

## 5. What in the brief's premises, and in the repo's notes, turned out to be wrong

1. **"a wall-to-wall chain of `T` among the near-axis squares in the x- or y-DAG" (the brief;
   `proof-architecture.md` §0a item 1; `proof-architecture-review.md` §5).**  The global-axis
   ordering proves the **order** form, not the DAG form; an order edge need not be a DAG edge at
   `eps > 0`.  The DAG form is provable by the same counting (§1.4) but then the tilt bound is used
   for *acyclicity*, not for the height, and the inequality the DAG chain supports is a factor
   `T-1` weaker.  Both forms coincide at `eps = 0`.  §1.4.
2. **"`eps < 1/(T-1)` is what forces the chain."**  It is not: chain existence is pure pigeonhole
   and survives arbitrary tilts (given acyclicity).  `eps < 1/(T-1)` is what forces the chain to be
   **maximum** (`h_x <= T`), i.e. to span the container.  The sharp constant is
   `eps_max(T) = arcsin(T/(sqrt2(T-1))) - pi/4`, which is `25.53°` at `T = 4` where `1/(T-1)` is
   `19.10°`, and `eps_max - 1/(T-1) = O(1/T^2)`.
3. **"`L >= L*(T, t)` is the whole lemma" (`t3-chain.md` §0(1), for the two-level H).**  For the
   three-level ladder it is *two* integers: `L` against `T-2` **and** `B` against `0`.  A ladder
   with the right number of legs in the wrong place (`B > 0`, or the tilted block against the wall)
   is `+O(eps^2)`, never `-O(eps^3)`.  §2.2, §2.4.  This clause is new and is the one `T = 11`
   violates structurally.
4. **`T11_CHAINS.md` §3's "the content of a chain certificate is entirely in bounding `C`"**, read
   as "the rise is the only obstruction", is too narrow: the ladder never mentions the rise (it is
   a Farkas certificate, and `t3-chain.md` §1.4's `omega` identity says the same), and what the
   `T = 11` band denies it is *transverse rows with a common normal*, not rise control.  The two
   are the same obstruction seen from the primal and the dual side; the dual side is the one with
   a closed form.
5. **`T11_CHAINS.md` §0(3)'s sign convention is easy to misread.**  "The two interfaces cost
   `-0.16` and `-0.18`" are contributions to `L`, i.e. they make the bound *stronger*; a tilted
   square against a wall is wider and helps the certificate at first order.  What hurts is the band
   links' rise credit `C`.  (Nothing in that file is wrong; the sign is stated once and then easy
   to invert.)
6. **Not wrong, but worth recording:** the uniform-tilt three-level ladder does not reproduce the
   *value* of the `n = 6` / `n = 12` pinwheel duals, only their weights.  The pinwheel certificate
   closes with four wall rows; the generic ladder needs eight.  A cheaper closing is the obvious
   next thing to look for, and it is what would make the uniform-tilt ladder `O(t^2)` negative.
   §2.5.

---

## 6. Reproduce

    python3 search/ladder_hform.py check     > runs/ladder_check.txt      # every closed form, residual < 1e-12
    python3 search/ladder_hform.py cell89    > runs/ladder_cell89.txt     # against runs/t4_cycles_cell89.txt
    python3 search/ladder_hform.py pinwheel  > runs/ladder_pinwheel.txt   # against S6_LOCAL.md sec 3
    python3 search/ladder_hform.py tables    > runs/ladder_tables.txt     # the eps-orders in T
    python3 search/ladder_hform.py t11  runs/bandcut_scan_exact110.json > runs/ladder_t11.txt
    python3 search/ladder_hform.py t11  runs/bandcut_scan_exact132.json >> runs/ladder_t11.txt
    python3 -W ignore search/ladder_hform.py t11h runs/bandcut_scan_exact110.json > runs/ladder_t11h.txt
    python3 -W ignore search/ladder_hform.py t11h runs/bandcut_scan_exact132.json >> runs/ladder_t11h.txt

Total compute: about four minutes, one thread.  `search/ladder_hform.py` imports
`t3_chain_hform` (the closed forms it extends), and, for §3 only, `t3_chain` (the restricted-dual
LP) and `t11_chains` (the tight-chain enumerator); it modifies nothing.
