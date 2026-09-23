# bandcut-cost: what a cut through every row and column costs, in length and in area  (2026-09-20)

Task `tasks/bandcut-cost/README.md`.  Labels used throughout: **[proved]** = one line of
arithmetic or a checked algebraic identity; **[measured]** = a number produced by a script here or
in `runs/`; **[pub]** = from Arslanov–Mustafin–Shangitbayev EJC 28(4) 2021 P4.22 as transcribed in
`search/BANDCUT_SCAN.md`; **[heuristic]** = a fit or an extrapolation over published data;
**[guess]** = neither.

Checks: `search/bandcut_cost.py` (no dependencies, runs in < 1 s).

## 0. Verdict

1. **The chain inequality generalises cleanly to mixed tilts (MT, §1) and it is `T`-free.**
   Mixed tilts cost `(W(D)-1)/2 ~ |D|/2` per link — *first* order in the tilt difference, against
   *second* order for a uniform tilt.  MT's numerator is identical to `ARCH_TLEDGER.md` §2.1's
   `B_T`; its denominator is larger by `2 cos t`, so MT is strictly sharper.  **[proved]**
2. **The chain route dies at mixed *normals*, not mixed tilts.**  The sum telescopes only for a
   common separation normal; a spread of `beta` in the normals adds `sin(beta)(k-1) tau`,
   `tau = (1+sqrt2)/2`, to the numerator, which at `beta = 1°`, `T = k = 12` is `+0.018` — enough to
   destroy any certificate.  That is the exact mechanism by which the Cleemann/Arslanov band
   escapes.  **[proved modulo the `tau` bound]**
3. **The length/chain accounting produces no threshold, near `11–12` or anywhere** (§2.3): gain and
   cost are both linear in `T` and the gain's coefficient is the free parameter `eps`.  Said plainly
   as the brief asks: **candidate (a) is not the `T`-dependent statement.**
4. **(a) and (b) are the same *structure* (a tilted/squeezable set meeting every row and column) in
   two currencies — length and area — but they are not the same statement, and only (b) carries
   `T`,** because only (b) is measured against a budget with a fixed coefficient (`waste = T`).
5. **The waste accounting produces exactly `12`, and it does so by integrality, not by a balance.**
   New here: every squeezable rectangle in the literature satisfies `waste = min(side) + 2` and
   `max >= min + 4` (10/10, `check E`), which makes the budget identity
   `waste(A) + waste(B) = T` an identity for *every* `T`; the only obstruction below `12` is
   `a_A, a_B >= 4`, i.e. `4 + 4 <= T - 4`.  **[proved from the fit; the fit is [measured]]**
6. **And that is one step too many.**  `T = 11` is a verified positive (`delta = +2.749e-04`,
   exact), so "`delta > 0 => T >= 12`" is **false**.  The clauses that fail at `T = 11` are exactly
   the two unproved structural ones (the rectangular decomposition and the exact integer
   complement).  **The waste route is a theorem about the Arslanov scheme, not about packings**
   (§6.1).
7. **A structure-free waste bound is refuted by `s(5)`**: `s(5) = 2 + 1/sqrt2` wastes `2.3284`, so
   no absolute constant `W_0 > 2.33` is true (§6.2).  The `n = T^2 - T` hypothesis is load-bearing.
8. **Net recommendation.**  Neither candidate is the missing statement as posed.  The one place the
   accounting points to a genuine, smaller, provable question is `(S3a)`: **no `(3, b)` rectangle is
   squeezable, for any `b`** — and its `T = 4` shadow, which MT can almost reach.  The rest of
   Lemma S (§5) is a correct derivation from two clauses that are false at `T = 11`.

## 1. The mixed-tilt chain inequality

Setting: the disjunctive LP of `search/S6_SKELETON.md` §3.1.  A pair `(i,j)` is separated when, for
some `o in {i,j}` and one of `o`'s two edge normals `e`, `|e.(c_j - c_i)| >= 1/2 + W(theta_j -
theta_i)/2 + delta`, `W(D) = |cos D| + |sin D|`.  Write `m(D) = 1/2 + W(D)/2` for the link
half-width; `m(0) = 1`, `m(45°) = (2+sqrt2)/4 = 1.2071`.

> **MT (mixed-tilt chain bound).**  **[proved]**  Let `Q_1, ..., Q_k` be squares of a packing of
> `[0,T]^2` with margin `delta`, with `Q_1` against the wall `x = 0` and `Q_k` against `x = T`, and
> suppose consecutive pairs are separated along a **common** unit normal `e = (cos phi, sin phi)`,
> `cos phi >= 0`.  Put `d_i = c_{i+1} - c_i`, `R = (c_k - c_1)_y` (the *rise*),
> `ubar = (u_1 + u_k)/2` with `u = cos theta + sin theta`, and
>
>     M  =  sum_{i=1}^{k-1} m(theta_{i+1} - theta_i)  =  (k-1)/2 + (1/2) sum_i W(theta_{i+1}-theta_i).
>
> Then
>
>     delta  <=  MT(k, phi, R, theta)  :=  [ cos(phi) (T - ubar) + sin(phi) R - M ]
>                                          / ( k - 1 + 2 cos(phi) ) .

*Proof.*  Sum the `k-1` link rows: `e.(c_k - c_1) = sum_i e.d_i >= M + (k-1) delta`.  The two wall
rows give `(c_1)_x >= u_1/2 + delta` and `(c_k)_x <= T - u_k/2 - delta`, so
`(c_k - c_1)_x <= T - ubar - 2 delta`; substitute, using `cos phi >= 0`, and solve for `delta`. ∎

**Checks (`bandcut_cost.py B`).**  At common tilt `t`, `phi = t`, `k = T` the numerator is
*identical* to that of `ARCH_TLEDGER.md` §2.1's `B_T = [(T-u)cos t + R sin t]/(T-1) - 1`, verified
at `T = 3,4,5,12` and four tilts and five rises; MT's denominator is larger by `2 cos t`, so **MT is
the sharper form** — it charges `delta` at the two walls as well as at the `k-1` links, which
`B_T` drops.  Signs agree everywhere.

**Three readings of MT.**

1. **The rise is still the only free quantity** (`dMT/dR = sin phi / (k-1+2cos phi) > 0`), so
   `ARCH_TLEDGER.md` §2's picture survives verbatim under mixed tilts.
2. **Mixed tilts cost at first order.**  `W(D) - 1 = |cos D| + |sin D| - 1 = |D| - D^2/2 + O(D^3)`
   (`check A`: `0.01730` at `1°`, `0.3529` at `t_4 = 28.07°`).  Each tilt-mismatched link subtracts
   `(W-1)/2 ~ |D|/2` from the numerator.  A *uniform* tilt, by contrast, costs only `O(t^2)`
   (`B_T(1,t) = -c_T t^2 + ...`).  **[proved]**  So along a chain the penalty for angular variation
   is one order worse than the penalty for angular magnitude: `sum_i |theta_{i+1} - theta_i| / 2`
   is the exact first-order price of a chain that changes tilt.
3. **Straight-stack specialisation.**  For a straight stack (`theta_i = t`, `R = (k-1) sin t`,
   wall squares also at `t`) the numerator collapses to `cos t [ T - u - (k-1) cos t ]`, so
   **[proved]**
   
       MT > 0   <=>   k cos t + sin t  <  T ,
   
   verified over `(T,k) in {(4,4),(4,3),(12,4),(12,12)}` and six tilts, 24/24 agreements.  This is
   the `p cos t + sin t = p` of `BANDCUT_SCAN.md` §1 with the container side in place of `p`.

**Why the chain route does not extend to mixed *normals*.**  If the links use different normals
`e_i` the sum does **not** telescope: `sum_i e_i.d_i` is not a function of `c_k - c_1`.  One can
still bound it by writing `a.d_i = cos(psi-phi_i)(e_i.d_i) - sin(psi-phi_i)(f_i.d_i)` and using
`|f_i.d_i| <= tau := (1+sqrt2)/2` for a link chosen as the max-surplus axis; with `beta` the
half-spread of the `phi_i` this gives **[proved, modulo the `tau` bound which is
[heuristic]]**

    delta [ cos(beta)(k-1) + 2 cos(phi) ]  <=  cos(phi)(T-ubar) + sin(phi) R
                                               - cos(beta) M  +  sin(beta) (k-1) tau .

The repair term `sin(beta)(k-1) tau` is **first order in `beta` and proportional to `k`**: at
`beta = 1°`, `T = k = 12`, `theta = 0`, `R = 0` it moves the bound from `0` to `+0.0180`
(`check B`), and at `T = 4` from `0` to `+0.0127`.  **So the chain certificate survives only on
chains whose link normals are near-common, i.e. on near-uniform-tilt chains.**  That is the precise
technical reason the Cleemann/Arslanov band escapes: its squares alternate between two tilts and its
links alternate normals, and there is no wall-to-wall near-common-normal chain at all.

## 2. What the inequality says about a band: gain, and the three costs

### 2.1 The gain

A `(p,1)` stack at tilt `t` spans `p cos t + sin t`; its **slack** is `s(p,t) = p - (p cos t +
sin t)`.  In exact `Fraction` arithmetic at `t_p = 2 arctan(1/p)` (`check D`, `p = 2 … 12`):

    cos t_p = (p^2-1)/(p^2+1),   sin t_p = 2p/(p^2+1),
    s(p, t_p) = 0   and   ds/dt |_{t_p}  =  p sin t_p - cos t_p  =  1   for EVERY p.    **[proved]**

So a band chain of any length, tilted `eps` past its window angle, frees exactly `eps` of column
width, to first order, **independently of `p` and of `T`**.  Turned into margin by MT (straight
stack, `k = p = T`-crossing): `d(delta)/dt = cos t_p / (p - 1 + 2 cos t_p)`; at `p = 4`, `T = 4`
this is `15/81 = 0.1852`, so the measured `eps = +0.52°` of `BANDCUT_SCAN.md` §5 buys
`delta ~ 1.7e-3` on the stack's own row.  **[proved]**  (The configuration's measured `delta` is
`0`: the rest of the packing eats it — see §2.3.)

### 2.2 The three costs, priced by MT

| cost | where | price, from MT | at `t_4 = 28.07°` |
|---|---|---|---|
| **interface** | each contact between a band square and a near-axis square | `(W(t)-1)/2` of length per contact | `3/17 = 0.17647` |
| **elbow** | the band must turn to cut both rows and columns; the turn is a tilt change `D` inside the band | `(W(D)-1)/2` per link, `~ |D|/2` | `0.2071` at a `45°` elbow square |
| **walls** | the band's two ends pay `u/2 = (cos t + sin t)/2` instead of `1/2` | `(u-1)/2` per end | `3/17` per end |

All three are **first order in the angle**, while the gain (§2.1) is first order in the *excess*
`eps` only.  Since the interface angle is `t_p ~ 2/p` and the excess is `eps << t_p`, **every
contact with the near-axis field costs `~1/p` and every band chain gains `~eps`**: the band is
viable only if the number of interface contacts per chain is `O(eps p)`, i.e. only if the band is
*thick* — many band squares per unit of interface.  That is the qualitative content of "a band,
not a single tilted square".  **[proved for the per-contact prices; [heuristic] for the conclusion.]**

### 2.3 The net, as a function of `T` — and the honest negative

Write the length ledger for a band that cuts every row and every column of `[0,T]^2`:

* number of band chains needed to cut `T` rows: `~ T / (cos t + p sin t) ~ T/3` (the stack height is
  `(3p^2-1)/(p^2+1) in [2.2, 3)`, `check D`) — **linear in `T`**;
* gain: `~ eps` per chain → total `~ eps T / 3` — **linear in `T`**;
* interface length: the band runs from wall to wall in both directions, so `L ~ c T` with
  `c in [1, 2 sqrt2]` — **linear in `T`**; interface cost `~ (W(t_p)-1) L / 2` — **linear in `T`**;
* elbow and wall costs: `O(1)` — **constant in `T`**.

**Both sides are linear in `T` with constant corrections.**  A ledger of the form
`gain = a T + b` vs `cost = c T + d` changes sign at `T = (d-b)/(a-c)`, which is a threshold only
if `a - c` is small and `d - b` is a dozen times smaller still.  Nothing here fixes `a` and `c`
independently — `a` is proportional to the free parameter `eps`, which the optimiser sets — so
**the length/chain accounting does not produce a threshold at all, near `11–12` or anywhere.**
Stated plainly, as the brief asks: *(a) is not the `T`-dependent statement.*  It is a correct and
sharp local certificate (§1) with no `T`-threshold in it, which is exactly the diagnosis
`ARCH_TLEDGER.md` §2.3 and `notes/proof-architecture.md` §0a already reached for `c_T`.

Confirmation on the one exact instrument we have (`check C`, the measured `T = 4` optimum of
`BANDCUT_SCAN.md` §5, `delta = 0`): **[measured]**

| chain | rise `R` | MT bound on `delta` |
|---|---|---|
| the four axis-parallel squares of the bottom row | `0` | `0.000000000` (tight) |
| the four tilted squares of the `(4,1)` stack | `+1.6434` | `+2.2600e-02` (slack) |

The *band's own* chain is far from binding; what certifies `delta <= 0` at `T = 4` is the surviving
**axis** chain.  And mixing tilts along a 4-chain of the same geometry costs exactly what §1(2)
predicts (same `R`, same normals, alternating tilts):

| tilt spread | `sum (W-1)/2` | MT bound |
|---|---|---|
| `0°` | `0` | `0.000000` |
| `1°` | `0.02595` | `-0.005198` |
| `5°` | `0.12503` | `-0.025189` |
| `28.07°` | `0.52941` | `-0.110854` |

## 3. The area side: waste, and the exact arithmetic of `T = 12`

**The budget.**  `waste([0,T]^2) = T^2 - n = T` exactly, for `n = T^2 - T`.  **[proved]**  This is
the only quantity in the problem that is *linear in `T` with coefficient 1 and no free parameter*,
and it is therefore the only place a clean threshold can live.

**A fit nobody seems to have written down.**  Every squeezable rectangle in the literature
(`BANDCUT_SCAN.md` §4, Table 2 and Lemma 1 of the paper) satisfies, with `a = min(side)`,
`b = max(side)`:

>     waste(a, b)  =  a + 2        and        b  >=  a + 4 .          **[measured, 10/10]**

| `(a,b)` | `(4,8)` | `(5,9)` | `(5,10)` | `(6,10)` | `(6,11)` | `(6,12)` | `(8,12)` | `(10,14)` | `(12,16)` | `(14,18)` |
|---|---|---|---|---|---|---|---|---|---|---|
| squares | 26 | 38 | 43 | 52 | 58 | 64 | 86 | 128 | 178 | 236 |
| waste | 6 | 7 | 7 | 8 | 8 | 8 | 10 | 12 | 14 | 16 |
| `a+2` | 6 | 7 | 7 | 8 | 8 | 8 | 10 | 12 | 14 | 16 |

(`check E`; the `(2k, 2k+4)` rows are Lemma 1, `4k^2 + 6k - 2` squares, so the fit is a *theorem*
on that infinite family and a 5-point coincidence on the sporadic ones.)  `a >= 4` is the
literature's floor: no squeezable `(3, b)` is known, and `BANDCUT_SCAN.md` §5 found none.

**The arithmetic of `T = 12`, exactly.**  In the scheme `[0,T]^2 = A u B u C u D` with `C, D`
integer blocks (waste `0`) and `A, B` the two squeezable rectangles, `A` covering the bottom rows
and the left columns and `B` the top rows and the right columns:

    T  =  waste(A) + waste(B)  =  (a_A + 2) + (a_B + 2)  =  a_A + a_B + 4 ,     a_A, a_B >= 4
      =>  T  >=  12 ,   with equality iff  a_A = a_B = 4.                       **[proved from the fit]**

The equality case is forced and it exists: `A = (4,8)` at `x in [0,4], y in [0,8]`, `B = (8,4)` at
`x in [4,12], y in [8,12]`, `C = (4,4)`, `D = (8,8)`; squares `26 + 26 + 16 + 64 = 132 = 12^2 - 12`
(`check E`).  Every row of `[0,12]^2` is cut by `A` or `B` and every column likewise.  The paper's
two recipes, read through the fit, reproduce the budget identity `waste(A) + waste(B) = T` for
**every** `n` (checked at `n = 13, 14, 20, 30`, both parities) — so the scheme is *exactly critical
at every `T`*, and the only thing that fails below `12` is the integer constraint `a >= 4`.

**This is the key structural point.**  The waste accounting is not a balance that tips at `T = 12`;
it is an identity that holds for all `T`, constrained by two integers.  `T >= 12` is
`4 + 4 <= T - 4`.  The threshold is **integrality, not a margin**.  That also explains why
`BANDCUT_SCAN.md` §4 finds *two different* failure modes: at `T <= 5` the geometry (a `(4,8)`
does not fit in `[0,T]^2`), at `6 <= T <= 11` the budget.  Both are the same inequality
`a_A + a_B + 4 <= T` read with `a >= 4` or with `b >= a+4 <= T`.

## 4. Are (a) and (b) the same statement?

**Same *structure*, different *currency*, and only one of them carries `T`.**

* Same structure.  (a)'s hypothesis "no wall-to-wall chain of `T` near-axis squares in either
  direction" is, by MT §1 (3), exactly "the tilted set meets every row and every column".  (b)'s
  hypothesis "`A u B` covers the bottom-and-left and the top-and-right" is, in the scheme of §3,
  exactly "the squeezable set meets every row and every column" — `C` and `D` are integer blocks
  precisely because they meet no complete row or column by themselves.  The two hypotheses are the
  same set-theoretic condition on the tilted/squeezable set.  **[proved, given A4's conclusion]**
* Same mechanism.  MT prices one contact between a band square and a near-axis square at
  `(W(t)-1)/2` of *length*.  Integrated along the band's interface, that length is a *gap*, and
  the gap's area is waste.  So (b)'s `waste` is, up to the geometry of the interface, the line
  integral of (a)'s penalty.  **[heuristic]**  (Quantitatively the match is only order-of-magnitude:
  `(W(t_4)-1)/2 = 0.1765` per unit length against the fit's `1` unit of waste per unit of short
  side.  I did not close that factor and deliberately did not tune it.)
* Different currency, and this is what decides.  (a) is an inequality in *length*, homogeneous of
  degree 1 in everything, with a free parameter `eps` on the gain side; §2.3 shows both sides grow
  linearly in `T`, so no threshold.  (b) is an inequality in *area* against a budget `T` that is
  fixed with coefficient 1 and has no free parameter.  **Only (b) can carry `T`, and it does so
  through integrality (`a >= 4`), not through a balance.**

So: **not the same statement.**  They are two prices for the same object; (a) is the sharp local
certificate, (b) is the global budget.  The brief's question "which is closer to provable" has an
uncomfortable answer: **(a) is essentially proved already and is `T`-free; (b) has the right shape
and is not provable in the form that gives `12`** (§6).

## 5. The candidate lemma, with constants

The strongest thing the accounting suggests, stated so that every constant is visible.  I state it
for the shape the brief asked for and then say at once which clause is the fiction.

> **Lemma S (squeezable-transversal waste).**  *Constants:* `a_0 = 4`, `w_0 = 2`, `eta = 1°`.
> Let `P` be a packing of `n = T^2 - T` closed unit squares in `[0,T]^2` with margin `delta > 0`.
> Let `X` be the union of the squares whose tilt is `>= eta` (mod `90°`), together with every cell
> of the integer grid that `X` meets.  Then
>
> * **(S1)** `X` meets every horizontal and every vertical line of `[0,T]^2`.
> * **(S2)** `X` is contained in the union of two axis-aligned rectangles `A`, `B` with disjoint
>   interiors, `A` meeting all rows below some height and all columns left of some abscissa, `B`
>   the complementary rows and columns.
> * **(S3)** each of `A`, `B` is *squeezable*, and a squeezable rectangle with short side `a` has
>   `a >= a_0 = 4` and wastes at least `a + w_0 = a + 2`.
> * **(S4)** `[0,T]^2 \ (A u B)` is a disjoint union of integer blocks packed exactly.
>
> Then `T = waste(P) = waste(A) + waste(B) >= (4+2) + (4+2) = 12`.

**Constants, and where each comes from.**

| constant | value | status |
|---|---|---|
| link half-width `m(D) = 1/2 + W(D)/2` | `W(D) = |cos D|+|sin D|` | **[proved]** (S6_SKELETON §3.1) |
| mixed-tilt penalty per link | `(W-1)/2 = |D|/2 + O(D^2)` | **[proved]** (§1, `check A`) |
| stack window angle | `t_p = 2 arctan(1/p)`, `ds/dt|_{t_p} = 1` | **[proved]** (`check D`) |
| stack height | `(3p^2-1)/(p^2+1) in [2.2, 3)` | **[proved]** (`check D`) |
| min short side of a squeezable rectangle | `a_0 = 4` | **[measured]** — literature floor, no proof |
| waste offset | `w_0 = 2` | **[measured, 10/10]** — fit, theorem on the `(2k,2k+4)` family only |
| budget | `waste = T` | **[proved]** |
| resulting threshold | `T >= 12` | **[proved from the above]** |

**What would have to be proved, exactly.**

1. **(S1)** — *A4 in the shape the review left it.*  "If `P` has `delta > 0` then the near-axis
   squares contain no wall-to-wall chain of `T`", i.e. the negation of A2a's hypothesis forces a
   transversal tilted set.  Status: A4 is proved only for `k > (T-1)^2` near-axis squares
   (`proof-architecture.md` §0a item 1); at `T = 4` the hole `4 <= k <= 9` is open and inhabited
   (`BANDCUT_SCAN.md` §5 sits in it with `k = 8`).  **This is the known open lemma; nothing here
   changes its status.**
2. **(S1) with a *minimum tilt*.**  MT gives `delta <= 0` for any wall-to-wall chain of `k` squares
   at near-common tilt `t` with `k cos t + sin t >= T`; to conclude "tilt `>= eta`" one needs that
   the near-axis squares cannot themselves form such a chain, with `eta` quantified.  MT's §1
   mixed-normal repair says the chain certificate dies at a normal spread of `beta ~ 0.1°`, so
   `eta` cannot be taken small: **an honest `eta` is `O(1/T)` at best, i.e. `2 arctan(1/T)`.**
   Unproved, and the `beta`-sensitivity is an obstruction, not a gap.
3. **(S2)** — *the rectangular decomposition.*  Purely a property of the published family.  **There
   is no reason a packing must decompose this way, and §6 shows it does not.**
4. **(S3a)** `a >= 4`: prove no `(3,b)` rectangle is squeezable, for any `b`.  This is a real,
   one-dimension-smaller problem (`BANDCUT_SCAN.md` §4 calls it "the sharpest open sub-problem").
   MT gives the `T = 4` half of it: in `[0,3]^2` a `(3,1)` stack at `t_3 = 36.87°` has height
   `2.6`, and `BANDCUT_SCAN.md` §5 measures `delta = -1.20e-02` for the best `(3,1)` band at
   `T = 4`.  Not a proof.
5. **(S3b)** `waste >= a + 2`: prove the fit.  Lemma 1 of the paper gives it on `(2k,2k+4)` as an
   upper bound (a construction); the *lower* bound is entirely open and is the mathematical content.
6. **(S4)** the complement is exactly packed: equivalent to "all the waste is in `A u B`".

## 6. Filters, and what is wrong with the strongest form

### 6.1 Filter "must be false at `T = 12`" — Lemma S **fails by exactly one step**

Lemma S concludes `T >= 12`.  At `T = 12` that is satisfied and tight (§3's two-`(4,8)` witness),
so the lemma is consistent there.  But `T = 11` is a **verified positive**: `n = 110`,
`delta = +2.748946598e-04` in exact arithmetic (`BANDCUT_SCAN.md` §2).  **Lemma S is therefore a
false statement.**  It over-proves by exactly one unit of `T`.

Which clause fails at `T = 11` is not a mystery: the `T = 11` packing is 64 axis-parallel squares
plus **46** tilted ones in one *bent* band with two tilt clusters (`23.9–28.3°`, `57.8–65.1°`), not
two axis-aligned squeezable rectangles.  **(S2) and (S4) are false at `T = 11`.**  And they are the
two clauses of Lemma S with no proof attached.  So:

> **The waste-budget route, as the literature realises it, is a theorem about the
> Arslanov et al. *scheme* (which really does stop at `T = 12`), not about packings.**  Making it a
> theorem about packings requires replacing `a_0 + w_0 = 6` by a bound that admits the `T = 11`
> bent band, i.e. by a constant `<= 11/2 = 5.5` per transversal piece — and then the same lemma no
> longer yields `12`.  **[proved, from the exact `T = 11` certificate]**

This is the sharpest negative result of this note, and it is what the brief's filter was for.

### 6.2 Filter "must not prove `s(5) = 3`"

`n = T^2 - 4 = 5` at `T = 3`.  A *structure-free* waste lemma — "a strictly squeezable packing
wastes at least `W_0`" — applied to `[0,3]^2` with `n = 5` sees `waste = 9 - 5 = 4`, and with any
`W_0 >= 5` concludes `s(5) = 3`.  That is **false**: `s(5) = 2 + 1/sqrt2 = 2.70711`, whose waste is
`s(5)^2 - 5 = 2.3284`.  Hence

> **No absolute-constant waste lower bound with `W_0 > 2.33` is true.**  **[proved]**

Lemma S survives only because (S2)/(S4) exclude the `s(5)` packing (one `45°` square, no integer
blocks) and because it is stated for `n = T^2 - T`.  **The `n = T^2 - T` restriction and the
structure hypothesis are load-bearing, not cosmetic** — consistent with `S6_SKELETON.md` §2.2's
observation that the total-band inequality `sum_i D(t_i) < T(T-1)/2` is *exactly critical* at
`n = T^2 - T` and slack at `n = T^2 - 4` (slack `(T-4)/2`; at `T = 3, n = 5` it is vacuous).

### 6.3 A ledger that comes out with the wrong sign — recorded so nobody re-derives it

"Each tilted square costs `c` of waste; the transversal needs `N(T)` tilted squares; budget is `T`;
threshold where `c N(T) = T`."  Measured **[measured, `BANDCUT_SCAN.md` §2]**: `N(11) = 46`,
`N(12) = 51`.  Available waste per tilted square is `T/N`: `4/11 = 0.364` (`T = 4`, §5 of that
note), `11/46 = 0.2391`, `12/51 = 0.2353`.  It **decreases** with `T`: this ledger says large `T` is
*harder*, the opposite of the truth.  Any lemma of the form "waste `>= c` per tilted square" is
therefore refuted as a source of the threshold.  What large `T` actually buys is not a bigger waste
budget per tilted square but **room for a larger integer block**, i.e. a smaller *fraction* of
tilted squares — an integrality effect again, not a density one.

## 7. Reproduce

    python3 search/bandcut_cost.py all      # < 1 s, no dependencies
    python3 search/bandcut_cost.py A        # W(D) - 1 is first order
    python3 search/bandcut_cost.py B        # MT vs B_T; straight-stack sign; beta sensitivity
    python3 search/bandcut_cost.py C        # MT on runs/bandcut_scan_pT_T4.json
    python3 search/bandcut_cost.py D        # s(p,t_p) = 0, ds/dt = 1, exact Fractions
    python3 search/bandcut_cost.py E        # waste = min+2, the budget identity, the filters

Inputs used: `runs/bandcut_scan_pT_T4.json` (the measured `T = 4` band optimum), Table 2 and
Lemma 1 of the paper as transcribed in `search/BANDCUT_SCAN.md` §4, and the exact `T = 11, 12`
certificates of `BANDCUT_SCAN.md` §2.  Nothing was fetched.

**Not done / honest limits.**

* The conversion factor between MT's per-contact length penalty `(W-1)/2` and the fit's `1` unit of
  waste per unit of short side is not closed (§4).  I stopped rather than tune it, per the brief.
* `tau = (1+sqrt2)/2` in the mixed-normal repair (§1) is argued from "the chain link is the
  max-surplus axis"; it is not a clean lemma.
* The fit `waste = a + 2` is a theorem only on the `(2k, 2k+4)` family (where it is Lemma 1 read
  backwards) and a 5-point coincidence on the sporadic rectangles.  No lower-bound proof exists for
  any of them.
* No new configuration was searched; every geometric number here is from `runs/`.
