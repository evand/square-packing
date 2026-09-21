# s6-skeleton: a separation/transversal branch-and-bound for `s(6) = 3`, and why it does not close

*2026-09-20.  Task `tasks/s6-skeleton/README.md`.  Code: `search/s6skel.py` (new; `rank8.py`,
`chains.py`, `leaf_ceiling.py` imported, not modified).  Runs are in the gitignored `runs/`
(`s6_margin.jsonl`, `s6_margin_shuf.jsonl`, `s6_cores_shuf.txt`, `s6_zeroset.txt`,
`s6_scaling.txt`, `s6_killradius.txt`, `s6_narrowbox.txt`, `s6_selftest.txt`, `s6_bb.jsonl`,
`s6_bb.log`, `s6_bb4.jsonl`, `s6_bb4.log`), so every number is quoted here.*

Semantics throughout as `notes/s13-casefree.md` §1: closed unit squares, closed containment,
*packing = pairwise disjoint as closed sets*.  The theorem being re-proved is

> **six closed unit squares, pairwise disjoint as closed sets, do not fit in `[0,3]^2`**
> (Kearney–Shiu 2002, `notes/proof-anatomy.md` §5.2).

---

> **Correction (2026-09-20, `S6_LOCAL.md`).**  Three numbers below are superseded: the `n = 12` uniform-tilt
> coefficient is `-(1/3) t^2`, not `-(2/3)` (§6; `fixed_angle_value` was under-optimised — the `0.8°` row already
> showed it); the `n = 6` cubic constant is `-(1/4) t^3`, not `-0.17` (§0.4, §4.2; two of the three points were
> below the LP tolerance); and at `n = 12` there is no cubic direction to measure (§6, last paragraph).  The
> verdict does not depend on any of them.

## 0. Verdict, up front

**NOT CLOSED, and the skeleton provably cannot close by refinement.**  What came out instead is
sharper than a leaf count, and all four items transfer to `n = 12`:

1. **The problem is `n`-dimensional, not `3n`-dimensional.**  With the *angles* fixed, containment
   is linear in the centres and pairwise disjointness is a disjunction of linear inequalities
   (separating-axis theorem).  So `delta*(theta) = max over centres of min slack` is an exact
   disjunctive LP, and

       s(6) = 3   <=>   delta*(theta) <= 0  for every theta in [0, 90 deg)^6 .

   A complete decision procedure for one `theta` is implemented and validated to 9 decimals
   against an independent instrument (§3).  At `n = 12` the same reduction turns `s(12) = 4` into
   a statement in **12** variables, not 36.

2. **The kill lemmas.**  The brief's K1 band `(cos t - sin t)/2` is the wrong constant — it
   vanishes at `45 deg` and kills nothing there.  The right band is `D(t) = p(t) - |cos t||sin t|`
   (`notes/chord-lemma.md` Cor. 2), which is `>= (sqrt2 - 1)/2 = 0.2071067811865476` at **every**
   tilt.  K1 then generalises to a one-parameter family **K2** (the *fractional chord lemma*,
   §2.2): *no `k` squares share a point of their `lambda`-bands whenever `k*lambda >= T`.*  Both
   are exact, both have open complements (so they are box-checkable), and the chain length enters
   as the single inequality `k*lambda >= T` (§2.3 — this is what makes them silent at `n = 17`).

3. **The residue, geometrically.**  `delta*(theta) = 0` is not generic: over 80 random angle
   vectors with `j` of the six angles forced to a multiple of `90 deg`, the value `0` was reached
   `0/80` times for `j = 0, 1, 2`, `3/80` for `j = 3`, `49/80` for `j = 4`, `80/80` for
   `j = 5, 6`.  Generic angle vectors are strictly negative (40 uniform samples: max
   `-0.029555`, median `-0.113101`, min `-0.155984`).  **The residue is a neighbourhood of the
   codimension-3 set `Z = {theta : at least three angles are multiples of 90 deg}`** — exactly the
   angles of the wall-to-wall row of three that the core analysis finds under 98.5 % of the zero
   set (§1.2).

4. **Why refinement diverges — the numbers that matter.**  A box *meeting* `Z` has
   `sup delta*_true = 0` on it, so a sound method kills it only if its relaxation is **exactly
   tight** there; every interval/cone relaxation errs strictly positively, and the error is
   measured on `[0, w]^6` (where the sup is exactly `0`, so the reported value *is* the error) to
   be `+0.4984 w^2 -> +w^2/2`.  A box *near* `Z` is killed only if that error is below the true
   deficit, and along the worst direction out of `Z` — tilting only the three squares of one
   wall-to-wall row — the deficit is measured to be a mere `-0.17 t^3`.  So a box of width `W` at
   distance `rho` from `Z` needs `W^2/2 < 0.17 rho^3`, i.e. `W/rho < 1.7 sqrt(rho)`: **the
   refinement factor needed grows like `rho^{-1/2}`, so the box count diverges as the tree
   approaches the zero set.**  Same phenomenon as `RANK8.md` §3.3 at `n = 12` ("the second-order
   remainder is positive"), now with both constants measured and the two exponents separated:
   error `2`, deficit `3`.

The honest one-line summary: *the skeleton is a correct and complete decision procedure at every
fixed angle vector, and margin zero defeats it on every angle box of positive width.*  §6 says
what would have to replace box arithmetic.

**Arithmetic disclaimer.**  Every "kill" below is a floating-point LP certificate (HiGHS through
`scipy.optimize.linprog`).  The *lemmas* K1 and K2 are exact and stated with proofs; the *box
counts* are not rigorous — making them rigorous means re-deriving each Farkas multiplier in
rational arithmetic, which is not done here.  No claim of a proof of `s(6) = 3` is made.

---

## 1. The instrument: the picture at `n = 6`, `T = 3`

    python3 search/s6skel.py margin --nstart 20000 --nudges 8 --kick-rounds 3 \
        --kick-seeds 64 --kicks 24 --nproc 14 --out runs/s6_margin.jsonl
    python3 search/s6skel.py cores runs/s6_margin_shuf.jsonl --limit 2000 --ntry 2

### 1.1 Margin exactly zero, on a plateau

Max-min of the closed pairwise gap (`rank8.pair_gap`, which is `T`-free) plus the wall slacks,
by assignment-fixed SLSQP from 84 tiling sub-configurations, their nudges, 20,000 random starts
and 3 rounds of kicks — **25,364 local optimisations, 13 s on 14 processes**:

| quantity | value |
|---|---|
| best value over all starts | `+0.000000000000e+00` |
| starts landing at `\|value\| <= 1e-9` | `20,730` |
| starts with value `> +1e-12` | `0` |
| largest tilt on the zero set | `45.0000 deg` |
| zero-margin points further than `0.05` (max-norm, centres and tilt) from every tiling sub-configuration | `16,852` of `20,730` |

This is the `n = 12` picture of `RANK8.md` §1 reproduced at `n = 6`: margin exactly `0`, attained,
on a plateau that does **not** collapse to the tiling and that reaches the full `45 deg` of tilt.

### 1.2 What holds a plateau point shut: a wall-to-wall row of **three**

`chains.py analyse` (its `T` rebound to `3`; the first-order analysis itself is `T`-free) on
2,000 zero-margin configurations drawn at random from the whole run:

| irreducible first-order core | count |
|---|---|
| one straight wall-to-wall row of three axis-parallel squares (`x3LR` / `y3BT`) | **1,930** (992 + 938) |
| signature containing such a row | 1,942 |
| corner-contact compound core (the "pinwheel") | 11 |
| core reported `TILTED` | 8, all with max tilt `<= 0.024 deg` (numerically axis-parallel) |
| flagged first-order feasible and skipped | 41, all at `<= 1.103e-07` (optimiser stopping short) |
| analysed | 1,959 |

**No core ever needed a tilted square** — identical to `CENSUS.md` §2 at `n = 12` (99 % straight
row of four, the rest corner-contact, none tilted).  The `n = 3` rehearsal has the same two core
types as `n = 4`, with `3` for `4` throughout.

---

## 2. The exact lemmas

### 2.1 K1 — the transversal chord lemma, with the right constant

For a unit square at angle `t` write `C = |cos t|`, `S = |sin t|`, `p = (C+S)/2` (its half-width
and half-height) and

    D(t) = p - C S       in  [(sqrt2 - 1)/2 , 1/2] ,   = 1/2 iff t = 0 mod 90 deg,
                                                      = (sqrt2-1)/2 iff t = 45 deg mod 90 deg.

> **K1.**  Let `Q_1, ..., Q_k` be closed unit squares in `[0,T]^2`, pairwise disjoint as closed
> sets.  If there is a height `a` with `|a - (c_i)_y| <= D(t_i)` for every `i`, then `k < T`.
> The same holds for vertical lines with `(c_i)_x`.

*Proof.*  `notes/chord-lemma.md` Lemma 1 gives the chord of `Q_i` on `y = a` as
`min(1/C_i, 1/S_i, (p_i - |a - (c_i)_y|)/(C_i S_i))`, and Corollary 2 says it is `>= 1` exactly
when `|a - (c_i)_y| <= D(t_i)`.  So the `k` chords are pairwise disjoint closed intervals
`[alpha_i, beta_i]` of length `>= 1` inside `[0,T]`.  Order them: disjointness of closed
intervals forces `alpha_(i+1) > beta_(i) >= alpha_(i) + 1`, so
`beta_(k) >= alpha_(k) + 1 > alpha_(1) + k >= k`, and `beta_(k) <= T`. ∎

`search/s6skel.py lemma` audits the constants: `D` over a grid of 200,001 angles has minimum
`0.2071067811865475` at `45.0000 deg`, against `(sqrt2-1)/2 = 0.2071067811865476`; the brief's
proposed band `(cos t - sin t)/2` is `0.5, 0.3536, 0.1830, 0.0616, 0.0123, 0` at
`0, 15, 30, 40, 44, 45 deg` — it is strictly smaller everywhere and **zero at `45 deg`**, so it
would kill nothing on the tilted part of the zero set.  `D` never drops below `0.207`.

*This is Nagamochi's Lemma 2 threshold and Stromquist's `y = 0.9` cut, restated for an arbitrary
line instead of a wall strip.*  In particular K1 **subsumes** `notes/chord-lemma.md` B1 (the
wall-strip lemma): a square with `p <= c_y <= 1` has `|0.9 - c_y| <= D(t)` for every `t`
(that file's Lemma 3), so B1 is K1 at `a = 9/10`.

### 2.2 K2 — the fractional chord lemma (the answer to "what kills the pinwheel")

The brief asks what replaces K1 when a square only *pokes a corner* across the line.  The answer
is to stop counting chords and start measuring them.  For `lambda` in `(0,1]` put

    D_lambda(t) = p(t) - lambda * C * S      ( >= D(t), with equality at lambda = 1 )

> **K2 (fractional chord lemma).**  Let `Q_1, ..., Q_k` be closed unit squares in `[0,T]^2`,
> pairwise disjoint as closed sets, and let `lambda` in `(0,1]`.  If there is a height `a` with
> `|a - (c_i)_y| <= D_lambda(t_i)` for every `i`, then `k * lambda < T`.  Equivalently:
>
>     *** no k squares of a closed packing in [0,T]^2 share a point of their lambda-bands,
>         for any (k, lambda) with k*lambda >= T ***

*Proof.*  As for K1: `|a - (c_i)_y| <= D_lambda(t_i)` makes the chord
`(p_i - |a - (c_i)_y|)/(C_i S_i) >= lambda`, and `1/C_i, 1/S_i >= 1 >= lambda`, so every chord has
length `>= lambda`; `k` pairwise disjoint closed intervals of length `>= lambda` in `[0,T]` force
`T >= beta_(k) > k*lambda`. ∎

At `T = 3` the family is `(k, lambda) = (3, 1), (4, 3/4), (5, 3/5), (6, 1/2)`.  Its content is
entirely about *tilted* squares, which is what `CENSUS.md` §3 asked for: at `t = 0`, `C*S = 0` and
every member has `D_lambda = 1/2`, so they coincide; at `45 deg` they separate —

| `t` | `D_1` | `D_{3/4}` | `D_{3/5}` | `D_{1/2}` |
|---|---|---|---|---|
| `0 deg` | 0.5000 | 0.5000 | 0.5000 | 0.5000 |
| `15 deg` | 0.3624 | 0.4249 | 0.4624 | 0.4874 |
| `30 deg` | 0.2500 | 0.3583 | 0.4232 | 0.4665 |
| `45 deg` | 0.2071 | 0.3321 | 0.4071 | 0.4571 |

**Two corollaries worth naming.**

* *(the axis-parallel theorem, with no point set)*  For axis-parallel squares `D = 1/2` exactly, so
  the `y`-bands are closed unit intervals; K1 says the interval family has depth `<= T-1`, so
  (interval graphs are perfect) it splits into `T-1` classes of pairwise band-disjoint squares;
  within a class every pair has `|Delta y| > 1`, and the same counting gives `<= T-1` per class.
  **At most `(T-1)^2` axis-parallel unit squares fit.**  `T = 3`: `4 < 6`.  `T = 4`: `9 < 12`
  (= `BENTZ.md` C1).  `T = 17`: `256 < 272`.
* *(a total-band inequality)*  Each of the `T-1` classes is a family of pairwise disjoint closed
  intervals in `[0,T]`, so its lengths sum to `< T`; adding up,
  **`sum_i D(t_i) < T(T-1)/2`**, strictly.  Since `D <= 1/2`, a packing of `n = T^2 - T` squares
  has `sum_i D(t_i) <= (T^2-T)/2 = T(T-1)/2` with equality iff every square is axis-parallel — so
  the inequality alone says *not all of them are*, at `T = 3`, `T = 4` and `T = 17` alike, and it
  says nothing more.  (It is exactly critical at `n = T^2 - T` for every `T`; that is one way to
  see why the `n^2-n` family is the hard one.)

### 2.3 Where the chain length `3` enters, and what the lemmas say at `n = 17`

K1 and K2 have exactly one quantitative hypothesis, `k * lambda >= T`:

| `T` | cheapest member | what it forbids |
|---|---|---|
| `3` | `(3, 1)` | 3 squares sharing a unit-band point |
| `4` | `(4, 1)` | 4 squares sharing a unit-band point |
| `17` | `(17, 1)` | 17 squares sharing a unit-band point |

At `T = 17` **no member of the family constrains fewer than 17 squares**, and the axis-parallel
corollary gives `256 < 272` — true, and silent about Cleemann's packing, which is tilted.  The
family therefore cannot prove `s(n^2-n) = n` uniformly: it says nothing at all about a packing
whose squares are spread over 16 band levels.  Sanity filter passed.

---

## 3. The skeleton: branch on the angles, decide the centres exactly

### 3.1 The reduction

> **Proposition.**  Fix angles `theta_1, ..., theta_n`.  Then `n` closed unit squares at those
> angles fit in `[0,T]^2`, pairwise disjoint as closed sets, **iff** `delta*(theta) > 0`, where
> `delta*(theta)` is the value of
>
>     maximise  delta  over the centres c_1..c_n, subject to
>       (walls)  p_i + delta <= (c_i)_x <= T - p_i - delta,  same in y
>       (pairs)  max over o in {i,j}, k in {0,1}, s in {+,-} of  s * n_{o,k}.(c_j - c_i)
>                    >=  1/2 + W(theta_j - theta_i)/2 + delta ,   W(D) = |cos D| + |sin D|
>
> and `n_{o,0} = (cos theta_o, sin theta_o)`, `n_{o,1} = (-sin theta_o, cos theta_o)`.

Containment is exact because a unit square's bounding box is `[c_x - p, c_x + p] x [c_y - p, c_y + p]`;
the pair rows are the separating-axis theorem for convex polygons, with `1/2 + W/2` the half-width
of the Minkowski difference along an edge normal of square `o`.  So the whole question is a
**disjunctive LP in `2n` variables**, and the angles are the only continuum left.

### 3.2 The decision procedure

Big-M MILP is useless here: on six squares at generic fixed angles HiGHS does not move the dual
bound off the cap in 30 s (measured).  What works is lazy branching (`Decider` in `s6skel.py`):
solve the LP with only the rows branched on so far, look at the centres it returns, find the pair
(or the K2 subset) it is cheating on, branch on that pair's disjunction.  Pruning at `delta <= 0`
is sound because dropping rows only increases `delta*`.

**Validation.**  `delta*` from the decider against an independent assignment-fixed LP ascent
(`fixed_angle_value`, multistart, shares no code path with the decider):

| configuration | decider | instrument |
|---|---|---|
| 4 squares, `T = 3`, all axis-parallel | `+0.333333333` | `+0.333333333` (`= 1/3`, exact) |
| 5 squares, `T = 3`, all at `45 deg` | `+0.050252532` | `+0.050252532` |
| 6 squares, `T = 3.2`, all axis-parallel | `+0.050000000` | `+0.050000000` (`= (3.2-3)/4`) |
| 5 squares, `T = 3`, all axis-parallel | `kill` (`delta* = 0`) | `+0.000000000` |
| 6 squares, `T = 3`, all axis-parallel | `kill` (3,213 LP nodes) | `+0.000000000` |
| 6 squares, `T = 3`, all at `45 deg` | `kill` (53,445 LP nodes) | `-0.102944` |

**Cost, and its limit.**  A kill is a proof by exhaustion of the disjunction, and it is not cheap:
`3,213` LP nodes at exact axis-parallel angles, `53,445` at exact `45 deg`, and at a *generic*
exact angle vector the search was still running at the `150,000`-node cap (`runs/s6_selftest.txt`)
— the theorem is true there (the instrument returns `-0.1`) but this decider did not finish
proving it inside the budget.  That is the practical ceiling of the method at `n = 6`.

All of this is re-checkable with `python3 search/s6skel.py selftest` (`runs/s6_selftest.txt`),
which also checks `chord_len` against a direct polygon–line intersection (max discrepancy
`4.44e-16` over 400 random `(theta, d)`), that `D >= (sqrt2-1)/2` and `chord(theta, D(theta)) = 1`
over 20,001 angles, and that the brief's band `(C-S)/2` is `<= D` everywhere and `0` at `45 deg`.

**Soundness.**  On the *full* angle cube `[0, 90 deg)^n` the decider must never kill, and does not:
4 in `[0,3]^2` survives at `delta = 0.1716`, 5 in `[0,3]^2` at `0.1314`, 6 in `[0,3.2]^2` at
`0.0191`, and 6 in `[0,3]^2` at `0.00155` (the cube contains the tiling, so it must survive).

### 3.3 The angle-box relaxation, and its error

Over a box of angles the normals sweep a cone and must be relaxed.  The naive bound
`|A| + sin(w)|B|` (`A`, `B` the coordinates of `c_j - c_i` in the mid-angle frame, `w` the
half-width) has error `O(w)` and is hopeless — at `11.25 deg` bins it killed `0` of the first 100
boxes.  The bound actually used splits interior from endpoint:

* the maximum over the cone is `sqrt(A^2+B^2)` if attained in the interior, and then
  `|A| >= sqrt(A^2+B^2) cos w`, so the row `|A| >= m cos w` is valid;
* otherwise it is attained at `mid ± w`, and the row at that angle is **exact**.

Error `m (1 - cos w) = O(w^2)`.  The three interval constants the box relaxation needs —
`p^lo`, `W^lo`, `D_lambda^lo` — must *under*-estimate on every interval, and `selftest`
brute-forces that against a 601-point grid on 2,000 random intervals and six values of `lambda`
(`0` violations).  One real trap was found and fixed there: `D_lambda` is concave in
`u = |cos| + |sin|` with `d/du = 1/2 - lambda u`, so its minimum sits at `45 deg` only when
`lambda >= 1/2`; the `T = 4` family reaches `lambda = 1/3`, where the minimum is at `0 deg`
instead, and the naive "minimum at `45 deg`" rule would have produced *unsound kills* at `n = 12`.

Measured on `[0, w]^6` (`runs/s6_scaling.txt`):

| box `[0,w]^6` | `delta_relaxed` | `/ w^2` |
|---|---|---|
| `w = 4 deg` | `+2.199e-03` | `0.4511` |
| `w = 2 deg` | `+5.784e-04` | `0.4747` |
| `w = 0.5 deg` | `+3.758e-05` | `0.4935` |
| `w = 0.25 deg` | `+9.457e-06` | `0.4967` |
| `w = 0.125 deg` | `+2.372e-06` | `0.4984` |

> **relaxation error `-> + w^2 / 2`.**

(These are the first witness the DFS returns, so they are lower bounds on the relaxed optimum;
the `w = 1 deg` row of the run, `0.2414`, is such an under-report.)

---

## 4. What the skeleton kills, and what is left

### 4.1 The zero set in angle space

`runs/s6_zeroset.txt`, `fixed_angle_value` with 60 multistarts, 80 samples per row: `j` of the six
angles forced to `0`, the other `6-j` uniform on `[0, 90 deg)`.

| `j` | `delta* = 0` (to `1e-9`) | largest `delta*` seen |
|---|---|---|
| 0 | 0/80 | `-0.021494158` |
| 1 | 0/80 | `-0.021070888` |
| 2 | 0/80 | `-0.001047479` |
| 3 | **3/80** | `-0.000000000` |
| 4 | 49/80 | `+0.000000000` |
| 5 | 80/80 | `+0.000000000` |
| 6 | 80/80 | `+0.000000000` |

> **`Z = {theta : delta*(theta) = 0}` is contained in the codimension-3 set
> `{theta : at least three of the six angles are multiples of 90 deg}`, and is a proper subset of
> it (with `j = 3` only `3/80` samples reached `0`: the three free squares still have to fit).
> The residue of the skeleton is a neighbourhood of `Z`.**

That is exactly the core analysis of §1.2 read in angle coordinates: a zero-margin configuration
carries a wall-to-wall row of three *axis-parallel* squares (98.5 %) or a corner-contact compound
core of four or five axis-parallel squares (0.6 %); the remaining squares are unconstrained and may
take any angle.  So `Z` is a union of three-dimensional strata, and the pinwheel sits on the
smaller (codimension 4 or 5) ones.

### 4.2 How fast the margin decays off `Z` — and why the pinwheel is the hard direction

`runs/s6_zeroset.txt`:

| direction out of `Z` | measured | order |
|---|---|---|
| all six squares tilted by `t` | `delta* = -0.749 t^2` at `t = 0.05 deg`, `-0.748` at `0.1 deg`, `-0.746` at `0.2 deg` | **second**: `-(3/4) t^2` |
| only three of the six tilted by `t` (the other three left axis-parallel) | `delta*/t^3 = -0.169, -0.164, -0.236` at `t = 0.1, 0.4, 1.6 deg` (`t` in radians) | **third**: `~ -0.17 t^3` |

The configuration attaining the first row **is the pinwheel**, checked directly: the optimum at a
uniform tilt of `0.5 deg` (200 LP ascents, `delta* = -5.6371e-05`) is

        y = 2.4957                   x = 1.49567    x = 2.49566
        y = 1.5000    x = 0.50434    x = 1.50440
        y = 0.5043    x = 0.51307                   x = 2.48708

— the `3 x 3` tiling minus a cyclic permutation: two squares in every row, two in every column,
**no straight wall-to-wall row of three at all**.

Why it beats the row: a `3 x 2` grid (two full rows of three) at a common tilt `t` has centres
forming, in the tilted frame, an axis-aligned `2(1+delta) x (1+delta)` rectangle that must fit in
a square of side `3 - cos t - sin t` tilted by `t`, i.e.
`(1+delta)(2 cos t + sin t) <= 3 - cos t - sin t`; differentiating at `t = 0` gives
`delta <= -t + O(t^2)` — **first order**.  The pinwheel reorganises and pays only `-(3/4) t^2`.
**That is the precise sense in which the pinwheel is the escape the brief asks about**, and it is
also why `RANK8.md` §3 finds the obstruction first order but cannot finish: on the pinwheel the
first-order term is zero, and the whole question lives in the second- and third-order terms that
§3.3 there bounds but does not sign.

The second row is the direction that defeats box refinement.  It keeps the maximum number of
axis-parallel squares that `Z` requires (three) and tilts the rest, so it leaves `Z` while staying
inside the codimension-3 set that contains it — and it costs only *third* order to do so.

### 4.3 The refinement verdict

A box of width `W` whose closest zero-margin angle vector is at distance `rho` has
`delta_relaxed ~ sup_box delta*_true + W^2/2`, and it is killed iff `sup_box delta*_true < -W^2/2`.
Along the quadratic direction (`-0.75 rho^2`) that needs `W < 1.22 rho`: a *constant* refinement
factor, scale-invariant, fine.  Along the cubic direction (`-0.17 rho^3`) it needs

    W^2 / 2  <  0.17 rho^3      i.e.      W / rho  <  1.7 sqrt(rho) ,

so the refinement factor must grow like `rho^{-1/2}` as `rho -> 0`.  **The box count diverges.**

The two ingredients of that inference are each measured cleanly and independently: the error
`W^2/2` on `[0,w]^6`, where `sup delta*_true` is exactly `0` so the reported number *is* the error
(§3.3), and the cubic deficit along the core-tilt direction (§4.2).  The combination is an
inference, not a measurement.  Corroboration (`runs/s6_killradius.txt`, node cap 400,000): boxes
of full width `11.25 deg` centred at uniform tilt `c` all survive, with first-witness values
`+1.3162e-02` (`c = 5 deg`), `+6.4035e-03` (`c = 10 deg`), `+9.0212e-03` (`c = 20 deg`), against
`delta*` at the box *centre* of `-0.00499`, `-0.01735`, `-0.05128`.  Two caveats keep this from
being a clean fit: `sup_box delta*_true` is not `delta*` at the centre (a box of mixed angles can
be much less negative — the best of 40 uniform random vectors is `-0.0296`), and the decider
returns the *first* surviving witness, so every `delta_relaxed` quoted here is a lower bound on
the relaxed optimum.  Order of magnitude agrees; the exponents are what the argument rests on.

### 4.4 Box counts

Two levels were run.  Both are quotiented by the label symmetry (the angle intervals may be taken
non-decreasing) and by the container reflection `t -> 90 deg - t`.

| bins | box width `W` | boxes | processed | **N** | **K1** | **K2** | survived | node budget exhausted |
|---|---|---|---|---|---|---|---|---|
| 4 (`--warm 4`) | `22.5 deg` | 44 | 31 | `0` | `0` | `0` | `31` | `0` (median 141 LP nodes) |
| 8 (`--warm 8`) | `11.25 deg` | 868 | 550 | `0` | `0` | `0` | `105` | `445` (cap `60,000`) |

**Zero boxes killed at either resolution**, and the reason is §4.3 arithmetic, not a defect of the
lemmas: `W^2/2` is `0.077` at `22.5 deg` bins and `0.0193` at `11.25 deg` bins, while `|delta*|` is
`0.113` in the median and only `0.0296` at the top of 40 uniform samples — so at `22.5 deg` the
relaxation error exceeds the deficit *everywhere*, and at `11.25 deg` it exceeds it on the large
part of the cube that lies within roughly `16 deg` of `Z`.  Where the error is small enough, the
node budget becomes binding instead: `445` of `550` boxes at `11.25 deg` hit the `60,000`-node cap
(a kill needs `3,213` nodes at exact axis-parallel angles and `53,445` at exact `45 deg`, and an
angle **box** raises the per-pair branching factor from 8 to 24).

**Not resolved, honestly: whether the skeleton kills any box of positive width.**
`runs/s6_narrowbox.txt` puts the question directly — boxes of width `0.1`, `0.5` and `1 deg`
around `45 deg`, and `0.2 deg` around `20 deg`, `0.1 deg` around `5 deg`, i.e. exactly the regime
where the §4.3 arithmetic says a kill should be possible (`W^2/2` between `1.5e-6` and `1.5e-4`
against `|delta*|` of `0.103`, `0.051`, `0.005`).  The `W = 0` control **is** killed, in `53,445`
LP nodes and `44 s`.  None of the positive-width boxes returned inside the budget that could be
spent on them (see the note below), so the question stands open: the obstruction there is search
cost, not the lemmas.

*A note on the compute actually available.*  On this machine no single process survived much more
than a minute or two of continuous work — several detached runs were reaped mid-flight (the
`--warm 8` level reached `550` of `868` boxes before its pool died, and its workers had to be
cleaned up by hand).  Every number quoted in this file comes from a run that finished inside that
window; the ones that would need more are named as open, not estimated.

**So the box counts here measure the budget, not the mathematics**; the mathematics is §4.3.

### 4.5 Why K2 never fires in the branch-and-bound

K2 enters the LP as "some pair of this `k`-subset is `lambda`-band separated", i.e. as a *clique*
constraint on the band graph.  For `lambda < 1` the bands are wider, but the constraint requires
**all** pairs of the subset to be band-close before it bites, and the surviving witnesses avoid
that: in the `[0, 1 deg]^6` witness (the pinwheel, §5) the two squares of one row are `0.9827`
apart in `y` while `D_1 + D_1 = 0.9824` — K1 is satisfied by `2.8e-4` — and the `4`-subset formed
by two adjacent rows has a pair exactly `1.0000` apart against `D_{3/4} + D_{3/4} = 0.9911`.  The
configuration threads every member of the family.  The *chord-sum* form of K2
(`sum_i chord_i(a) < T`) is strictly stronger than the clique form, but `chord_i(a)` is concave in
`(c_i)_y`, so putting it in an LP needs a lower bound that is linear in `(c_i)_y`, i.e. branching
on the **centres** as well — which is the `3n`-dimensional problem the reduction of §3.1 was
meant to avoid.  That trade-off is the real open question this task produced.

---

## 5. The residue, drawn

The witness the relaxation returns on `[0, 1 deg]^6` (`delta_relaxed = 1.9037e-05`) is

        S4 (0.5000, 2.4827)   S0 (1.4999, 2.4999)
        S5 (0.5000, 1.4827)                         S2 (2.4823, 1.5000)
                              S3 (1.4999, 0.5000)   S1 (2.5000, 0.5002)

i.e. the `3 x 3` tiling **minus its main diagonal** — a pinwheel: every row holds two squares,
every column holds two, and there is no straight wall-to-wall row of three at all.  Three links
sit at `L_inf` distance `0.9998` against the required `1`, and they are of both kinds: `S0/S4`
(the two squares of the top row, separated in `x`), `S1/S2` (the two of the right column,
separated in `y`) and `S3/S5` (a **corner** contact, `dx = -0.9998`, `dy = +0.9827`).  Each
shaves `~2e-4` off a chain that has to span exactly `2` between the wall-offset centres — e.g.
`S5 -> S3 -> S1` runs from `x = 0.5000` to `x = 2.5000` in two steps.

So the residue of this skeleton is, exactly:

> **angle boxes that contain a pinwheel or a row-of-three angle vector — i.e. boxes meeting `Z`,
> a codimension-3 union of strata — together with a `sqrt`-shrinking collar around them whose
> width is set by the relaxation error `W^2/2` against a cubic deficit.**

It is not a third mechanism; it is the same mechanism (margin zero) seen through box arithmetic.

---

## 6. Transfer to `n = 12`

**What carries verbatim.**

1. *The reduction.*  `s(12) = 4` iff `delta*(theta) <= 0` for every `theta` in `[0, 90 deg)^12`,
   with `delta*` the value of a disjunctive LP in the **24** centre coordinates, `66` pair
   disjunctions, `48` wall rows.  This is a real dimension reduction of the object `RANK8.md`
   studies (36 dimensions) and it is exact.  It also makes `RANK8.md` §2's sub-configuration
   census cheaper to redo.
2. *K1 and K2*, with `T = 4`: `(k, lambda) = (4,1), (5,4/5), (6,2/3), (7,4/7), (8,1/2), ...`.  The
   `(4,1)` member is `notes/chord-lemma.md` B1 for an arbitrary line rather than a wall strip, and
   its axis-parallel corollary is exactly `BENTZ.md` C1's `<= 9`.
3. *The core taxonomy*: `CENSUS.md` §2 already reports the `n = 12` analogue (99 % straight row of
   four, the rest corner-contact, none tilted), so `Z_12` should be
   `{theta : at least four angles are multiples of 90 deg}`, codimension **4** in 12 dimensions.

**What changes, measured.**  The same instrument at `n = 12`, `T = 4` (`runs/s6_scaling.txt`):

| | `n = 6`, `T = 3` | `n = 12`, `T = 4` |
|---|---|---|
| `delta*` at uniform tilt `t` | `-(3/4) t^2` (`-0.74902` at `0.05 deg`) | `-(2/3) t^2` (`-0.66599` at `0.05 deg`, `-0.66531` at `0.1 deg`, `-0.66280` at `0.2 deg`) |
| only the core row tilted | `~ -0.17 t^3` | third order, but the `n = 12` multistart is not reliable at this scale (`delta*/t^3` = `-0.332` at `0.1 deg`, `-47.5` at `0.4 deg`, `-11.6` at `1.6 deg` — local optima, not a fit) |

The second-order coefficient at the real target is therefore `2/3`, and the *shape* of the
obstruction is identical.  `RANK8.md` §2's `3.763847e-06` margin on eleven squares is three orders
below `W^2/2` unless the box width is under `0.157 deg` (`W^2/2 = 3.76e-6` at
`W = 2.74e-3 rad`), so an `n = 12` angle-box tree would have to reach bins finer than `0.16 deg`
— about `573` bins per angle, in `12` angles — before it could even *see* that sub-family.

**Box-count estimate, for the record.**  Quotienting by the label symmetry only, a uniform
subdivision into `B` bins per angle gives `C(B+11, 12)` boxes; the residue (boxes meeting a
codimension-4 stratum) is about `C(B+7, 8)`.  At `B = 8` that is `50,388` boxes with a `6,435`
residue; at `B = 16`, `17,383,860` with `490,314`.  A single `n = 6` kill already costs `3,213`
to `53,445` LP nodes; at `n = 12` the lazy DFS has `66` pair disjunctions instead of `15`, so the
per-box cost is worse by orders of magnitude, and the only regime where the relaxation could kill
anything (`W^2/2` below the local `|delta*|`) needs `B` in the thousands.  **The angle-box tree is
not the way to spend `n = 12` compute.**

**What to do instead** (the one thing this task actually argues for).  The two exponents are the
whole story: the box relaxation errs at order `2` in the box width, and the true margin deficit is
order `3` in the distance from the zero set.  Any method that reads the angles only through
interval bounds will lose.  The missing ingredient is a *first-variation* statement — an exact
inequality of the form

> tilting the squares of a zero-margin configuration off the axis-parallel core (centres free)
> strictly decreases the available margin, by at least a cubic form in the tilts

— the measured instance being `delta* <= -0.17 t^3` when three of the six squares are tilted by
`t`, and `delta* <= -(3/4) t^2` when all six are.  The `n = 12` analogue of the quadratic constant
is `2/3`, measured above; the cubic one is the number to measure next.  That is a statement about
`T` squares and two opposite walls — `RANK8.md` §7's "smallest object anything in this repo has
reduced the gap to" — and it is exactly what Nagamochi's and Bentz's strictness lemmas do by hand
on the exact width `1/2 + (|cos D| + |sin D|)/2`.  The contribution of this note is that the
required order is now known (**cubic, not linear and not quadratic**) and that the constant is
measurable with `fixed_angle_value` at `n = 12` directly.

---

## 7. Reproduce

    # the instrument (13 s, 14 processes)
    python3 search/s6skel.py margin --nstart 20000 --nudges 8 --kick-rounds 3 \
        --kick-seeds 64 --kicks 24 --nproc 14 --out runs/s6_margin.jsonl
    shuf runs/s6_margin.jsonl -o runs/s6_margin_shuf.jsonl
    python3 search/s6skel.py cores runs/s6_margin_shuf.jsonl --limit 2000 --ntry 2

    # the chord constants, and the validation of everything the kills rest on
    python3 search/s6skel.py lemma
    python3 search/s6skel.py selftest

    # the branch-and-bound over angle boxes  (hours; launch detached)
    setsid nohup python3 -u search/s6skel.py bb --warm 8 --depth 2 --kinds N,K1,K2 \
        --node-cap 60000 --nproc 12 --out runs/s6_bb.jsonl > runs/s6_bb.log 2>&1 &

`Decider(n, T, kills=...)` decides one angle box; `fixed_angle_value(theta, T, rng)` is the
independent instrument; `D_lam`, `chord_D`, `chord_len`, `lam_family` are the lemma constants.
Nothing in `verify*/`, `certificates/`, `lean/` is touched, and no existing `search/*.py` is
modified.
