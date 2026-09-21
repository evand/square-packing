# Adversarial review of `notes/proof-architecture.md`  (task `tasks/arch-review`, 2026-09-20)

*Brief: assume the sketch is wrong and find where.  Everything below is derivation or a quote from a
file in the repo; the four numeric checks are elementary and reproduced at the end.  I edited no
existing file and committed nothing.  Sibling tasks `arch-farfield` and `arch-tledger` were in
flight while this was written; where I cite their `runs/arch_*` output I say so and treat it as
unfinished.*

---

## 0. Verdict

**The spine survives; the composition does not, and it fails at one arithmetic step.**

The note places the A3/A4 handover at `k = T` near-axis squares.  The only counting argument
available — the one A4 itself uses — places it at `k = (T-1)^2`.  That single misplacement is what
makes A3 look cheap (cap `3` rather than cap `9` at `T = 4`) and what makes the admitted hole look
like a corner case.  Put the threshold where it belongs and three things happen at once:

1. A4 becomes **free** in the presence of far squares (the counting never touches them), so the
   note's proposed repair "prove A4 with the far squares as obstacles" costs nothing — and buys
   nothing beyond `k > (T-1)^2`;
2. the hole becomes exactly `T <= k <= (T-1)^2`, which contains **every configuration the repo's own
   sanity filter is built from** (`s(5) = 2 + 1/sqrt2` at `T = 3`; Cleemann's `272` at `T = 17`) and
   the fat stratum of `Z_4` the note itself describes in §0.2 (four axis-parallel, eight free);
3. A3's cap is either unjustified (`T-1`) or vacuous (`(T-1)^2 = 9`, against an optimum that puts
   only `7.60` at `theta = 0`).

So the hole is not "one gap visible only when written down".  **The hole is the theorem.**  And the
note's guess about where `T` enters (§3.3 (β), "A4 leaks") is wrong in both directions: the
statement §3.3 names is false at `T = 4` as well as `T = 17`, and the corrected counting is
*provably* true at `T = 17`.

What I could not break: **S0** (checked, including the dilation constant), S1, S2, the `(T-1)^2`
axis-parallel corollary, the `K1`/`K2` silence at `T = 17`, and `tan(t/2) >= 1/T`.

---

## 1. The objections

### 1. The A3/A4 handover threshold is `(T-1)^2`, not `T`.  **FATAL**

A4's own argument, run on the near-axis squares alone, gives:

> if the `k` near-axis squares contain no wall-to-wall chain of `T` in `x` **and** none in `y`, then
> `h_x <= T-1` and every Mirsky class (an antichain in `x`, hence a `y`-chain) has `<= T-1` members,
> so `k <= h_x * max class <= (T-1)^2`.

Contrapositive: **a chain of `T` is forced exactly when `k > (T-1)^2`**, which at `k = n = T(T-1)`
is the note's A4 and at `T = 4` is `k >= 10`, not `k >= 4`.

The threshold is tight, and the witness is the repo's own filter case.  The optimal `5`-square
packing (`s(5) = 2 + 1/sqrt2 = 2.7071`) is four axis-parallel squares at
`x, y in {1/2, 1 + 1/sqrt2}` plus one at `45°`.  Its four near-axis squares occupy **two** `x`-levels
and **two** `y`-levels, so `h_x = h_y = 2 = T-1` and there is no chain of `3`; and `k = 4 = (T-1)^2`
exactly.  Dilating its centres by `T/s = 1.1082` (S0) puts it in `[0,3]^2` with
`delta = 0.0541 > 0`, matching S0's predicted `(T/s - 1)/2 = 0.05409709` to 15 digits.  So

> **"`k >= T` near-axis squares contain a wall-to-wall chain of `T`" is false already at `T = 3`,**

and it is false at `T = 4` too: nine axis-parallel squares centred on `{0.6, 2.0, 3.4}^2` are
pairwise disjoint in `[0,4]^2` with `delta = 0.1` (min pair gap `0.4`, min wall slack `0.1`) and
have three `x`-levels — no chain of four.

**Consequence for A3.**  If the composition is to be closed at all, A3 must carry
`k <= (T-1)^2`, so its LP row must read *"mass on poses of tilt `< eps` is `<= (T-1)^2 = 9`"*.  The
evidence the note offers is evidence for a different row: "the optimum puts `7.60` of its mass at
exactly `theta = 0`" (`notes/review-2026-09-13b.md`) means the cap `9` is **slack at the known
optimum** and deletes nothing.  The capped value is then still `>= 12.173398` (pure, certified) and
`>= 11.999999926` (corner `k = 4`, certified).  A3 as priced in the note buys the composition
nothing.

*Repair.*  Either (i) state A3 with cap `(T-1)^2` and accept it is a different and much harder
lemma — see objection 3, which says a cover LP cannot prove it; or (ii) keep cap `T-1` and own the
fact that the lemma covering `T <= k <= (T-1)^2` is now the entire content of the proof (objection 2).

### 2. A4 extends to "far squares as obstacles" for free — and that is the bad news.  **FATAL (re-scoping)**

§2 says "no version [of A4] exists when `T <= k < n`, which is most of `Z`", and proposes "A4 must be
proved with the far squares as obstacles".  But **the counting never mentions the far squares.**  It
uses only that the `k` near-axis squares are pairwise disjoint, inside `[0,T]^2`, and have tilt
`< eps`.  So A4 holds verbatim as

> **A4′.**  `k > (T-1)^2` near-axis squares, the remaining `n - k` at arbitrary angles
> `=>` some `T` of the near-axis squares form a wall-to-wall chain.

and the hole shrinks from `[T, n)` to `[T, (T-1)^2]` — at `T = 4`, `k in {4,...,9}`; at `T = 3`,
`k in {3, 4}`.  No further combinatorics can shrink it: objection 1's witnesses sit at the top of
that range.  So the "obstacle" version the note asks for is either free (A4′) or it is a
**geometric** statement that must use `n = T(T-1)` quantitatively — i.e. it is exactly the lemma the
sanity filter demands and nobody has.

The hole is also not thin.  The note's own §0.2 says `Z_4` is, near the origin, "some four squares
axis-parallel and the other eight, at arbitrary angles, still fit beside them"; `S6_LOCAL.md` §2
measures that stratum (`j = 4`: `delta* = 0` to `1e-15`, eight squares tilted to `12.8°` beside a
wall-to-wall chain of four).  That is `k = 4 <= 9`: **the note's description of `Z_4` is a
description of the hole.**  And `S6_LOCAL.md` §2 says a chain of four *is* present there — so in the
hole the chain exists and no counting argument can find it.

### 3. A capped cover LP cannot prove A3 with the cap the composition needs.  **FATAL (for the candidate, not the statement)**

Every degree-1 family at `t = 4` is `>= 12` because the LP can smear a zero-margin configuration
(`BENTZ.md` §0; `notes/status.md`).  A mass cap kills that only if **no** zero-margin configuration
satisfies the cap.  Two repo measurements say one does:

* `T = 3`, `n = 6`: `S6_LOCAL.md` §5 (2) — "all ten same-sign `sparse4` directions are `0`", i.e.
  `theta = (0,0,0,0,a,b)` with `a, b` same-sign up to `4°` has `delta* = 0` exactly.  That is a
  touching packing of six squares in `[0,3]^2` with `k = 4 = (T-1)^2` near-axis (any `eps < 4°`).
* `T = 4`, `n = 12`: `S6_LOCAL.md` §2 — `j = 4`, `delta* = 0` to `1e-15` with eight squares tilted to
  `12.8°`.  `k = 4 <= 9`.

Both satisfy the cap `(T-1)^2`.  Whatever the LP does with them, the packing side of the capped LP is
`>= n` at every `eps` for which they survive, so **the far field cannot be closed by
`capped LP < n`** at the honest cap.  A3 as a *statement* may still be true; A3 as a *certificate
scheme* is refuted.

Note also that the sibling `arch-farfield` run uses `cap = T-1 = 2` at `T = 3`
(`runs/arch_farfield_T3e*.out`): unconverged iterates (`M > 1`, so not bounds) read
`mass = 6.02–6.05` at `eps = 1°` with the cap active, `5.90` at `5°`, `5.66` at `10°`.  Even a
clean negative at that cap would not license the architecture, because `2` is not the cap A4 can
justify.  **Whatever that task returns, it must be re-run at `cap = (T-1)^2` before it bears on the
composition.**

*Repair direction.*  The note's parenthetical alternative — "no wall-to-wall strip holds `T`
near-axis squares" — is the right shape, because it is the chain statement rather than a count.  But
it is a clique-type rule over poses, and `ALLMEET.md`'s "sound cliques of any positive-volume rule
have honest cost `>= 18`" is the standing obstacle to putting it in this LP.  Say so in the note.

### 4. §3.3's guess (β) is wrong twice, and its Cleemann test uses the wrong number.  **FATAL (as a guess); the repair is a better test**

(a) The statement §3.3 says dies at `T = 17` — *"`>= T` axis-parallel squares contain a wall-to-wall
chain of `T`"* — is false at `T = 4`, by the nine-square witness of objection 1.  A statement false
at the target `T` is not a lemma of a proof of Claim(4), so it cannot be the lemma that carries
`T < T*`.

(b) The corrected counting is **provably true at `T = 17`**.  A chain of `17` near-axis squares
cannot sit in `[0,t]^2`, `t <= 17`; so any packing of `272` squares in side `< 17` has
`k <= (T-1)^2 = 256` near-axis squares, i.e. **at least `272 - 256 = 16` of Cleemann's squares are
not near-axis**.  A4/A4′ therefore survive at `T = 17` with room to spare and carry no
`T`-dependence at all.

So in the repaired architecture the `T`-carrier is the lemma covering `T <= k <= (T-1)^2` (the hole),
or A3 if the split is pushed up to `(T-1)^2` — **never A4**.  The note has it exactly backwards.

*Better Cleemann test* (replaces the one in §3.3): count the squares that are **not** axis-parallel
and compare with `n - (T-1)^2 = 16`.  If the tilted band has `>= 16` squares, the picture is
consistent and the failure is in the far field / the hole.  If it has `< 16`, then A4′ is false at
`T = 17` and the counting is the carrier after all — which would be a genuine surprise and worth
knowing.  Also worth measuring: the slack `17 - t` of the packing, since a global rotation by
`alpha <= (17-t)/t` turns it into a witness with `k = 0`, which would refute A3 at `T = 17`
directly if `alpha > eps`.

### 5. A4's `eps` is not `~1/T` by the route the note gives; there is a cleaner route where it is.  **REPAIRABLE**

The note's route is "for tilts `< eps ~ 1/T` every separating normal is within `eps` of an axis;
'`x`-separated' is acyclic; Mirsky".  The first clause is true (all eight edge normals of two
near-axis squares are within `eps` of `±x`, `±y`), and "every pair is `x`- or `y`-separated" is true
(SAT for disjoint convex polygons).  The other two are where it leaks: the normal belongs to *one of
the two squares*, so different pairs use different normals, and

    Delta_x  =  (n.Delta) cos phi - (n_perp.Delta) sin phi  >=  cos eps - sqrt2 * T * sin eps ,

which is not even positive unless `eps < 1/(sqrt2 T)`; forcing `h_x <= T` through this estimate
needs `eps = O(1/T^2)`.  Acyclicity has the same problem (a cycle closes at cost
`k * eps * diam`).

*Repair (clean, and it gives the constant the brief asks for).*  Order by the **global** axes, not
by the separating normal.  For any disjoint pair of squares with tilts `<= eps`, the SAT row gives
`n.Delta >= 1/2 + W/2 >= 1` with the normal at angle `|phi| <= eps`, and

    n.Delta  <=  max(|Delta_x|, |Delta_y|) * (cos phi + sin phi)  <=  max(|Delta_x|, |Delta_y|) * (1 + eps)

(`cos phi + sin phi <= 1 + phi` is checked numerically below), so

    max(|Delta_x|, |Delta_y|)  >=  1/(1 + eps)  for every disjoint pair.

Define `i <_x j` iff `c_{j,x} - c_{i,x} >= 1/(1+eps)`.  Then

* **acyclicity is free** (`Delta_x > 0` along every edge);
* **every pair is comparable** in `<_x` or `<_y`, by the display above;
* a `<_x`-chain of `h` has `(h-1)/(1+eps) <= T - 2p <= T - 1` (as `p = (cos+sin)/2 >= 1/2`), so
  `h <= T` as soon as **`eps < 1/(T-1)`**;
* a Mirsky antichain is pairwise `y`-separated, so sorting by `c_y` bounds it by `T` the same way.

So A4/A4′ are correct with `eps < 1/(T-1)` — `19.10°` at `T = 4`, `3.58°` at `T = 17` — and the
note's "`eps ~ 1/T`" is right in order but should be stated as `1/(T-1)` with this proof, not the
normal-classification proof.  This is the one item in the note marked **[new, unchecked]** that I
believe survives, in this form.

### 6. A2b's inequality is a strict over-claim at every `t > 0`, and `T` enters it at second order.  **REPAIRABLE**

A2b reads "chain length `>= T + (1/2) sum |Dphi| + (1/2)|phi_wall| - (staircase credit)".  Take the
model case the note itself analyses in §3.2: `T` squares at a common tilt `t` in a *straight*
wall-to-wall row.  There `sum|Dphi| = 0` and `|phi_wall| = 2t`, so A2b claims `length >= T + t`.
The exact width is

    (T-1) cos t + 2p  =  T cos t + sin t  =  T + t - T t^2/2 + O(t^3)  <  T + t   for every t > 0.

The over-claim is `T(1 - cos t) ~ T t^2 / 2`: at `T = 4`, `t = 1°` it is `6.1e-4`, three orders
above the `3.76e-6` margins `RANK8.md` §2 says matter and the same size as the deficits
`S6_LOCAL.md` measures.  So A2b is not merely loose past the window — it is false as an inequality
from second order on, and the error term is the one that carries `T`.

*Repair.*  State the chain bound as the exact width, `extent = T cos t + sin t + (links)`, whose
sign flips at `tan(t/2) = 1/T`; the chain certificate is valid exactly on `t < 2 arctan(1/T)`, with
margin `sin t - T(1-cos t) ~ t - T t^2/2`.  Then §3.1's row "first-order rigidity … mechanism is
`T`-free" should read "`T`-free at first order; `T` enters at order 2 through the wall anchors with
coefficient `T/2`".

Good news attached to this: the two radii are compatible.  `1/(T-1) < 2 arctan(1/T)` for every
`T >= 3` (`28.6° < 36.9°`, `19.1° < 28.1°`, `14.3° < 22.6°`, `3.58° < 6.73°` at `T = 3,4,5,17`), so
A4's `eps` sits inside A2b's window and A5 is not obstructed *on this side*.  See objection 8 for
the other side.

### 7. A2c is false as literally stated.  **REPAIRABLE**

"On the coherent cone, `delta* <= -c_T |tilt|^2`."  `S6_LOCAL.md` §2 reports `delta* = 0` to `1e-15`
at `theta = (0^4, t^8)` for `t` up to `12.8°` at `n = 12` — a point with `|tilt| = sqrt8 * 12.8°`
and `delta*` exactly `0`.  So `|tilt|` must mean the **chain's** tilts, with the other `n - T`
angles as free parameters.  Once that is fixed, the constant in A2c is not the uniform-tilt `c_T`
tabulated in §3.2, and §3.3 (α)'s test ("`c_T` changes sign `=>` A2c fails") is testing a different
coefficient from the one A2c needs.  The `c_T` ledger is only a statement about A2c's regime for
tilts below `eps(T) ~ 1/(T-1)`, which at `T = 17` means below `3.58°`.

### 8. A5 is a two-sided window and the note states only one side.  **REPAIRABLE — and this is where I would put the `T`-dependence**

A5 says "`eps` in A3 `<=` radius in A2b/c".  But A3's cap only bites for `eps` **bounded away from
zero**: as `eps -> 0` the cap forbids a vanishing band of poses and the LP simply moves mass from
`theta = 0` to `theta = eps` at `O(eps)` cost in the coverage rows, so the capped value tends to the
uncapped one.  (The sibling `T = 3` iterates are consistent: `~6.03` at `eps = 1°`, below `6` only
from `eps >= 5°`.)  So A5 needs

    eps_min(T)   <=   eps   <   min( 1/(T-1),  2 arctan(1/T) )  ~  1/T ,

a **window that closes** once `eps_min(T)` stops shrinking as fast as `1/T`.  That is a fourth
birth mechanism for `P_T`, missing from §3.3:

> **(δ) the glue fails: the `eps`-window `[eps_min(T), 1/(T-1)]` becomes empty at some finite `T`.**

It is the mechanism most consistent with everything else in the file — the lemmas stay individually
true, the tube shrinks like `1/T`, and the far field needs a fixed angular resolution — and it makes
a cheap, decisive measurement: **`eps_min(T)`, the smallest `eps` at which the capped far-field LP
drops below `n`, at `T = 3` and `T = 4`, against `1/(T-1)`.**  I would put that in §4 in place of the
uncapped-`eps` sweep.

### 9. §3.1's `n = T^2 - 4` reasoning is a non-sequitur.  **COSMETIC, but it produces a wrong to-do**

The row reads: "`5 > 4` still forces a chain of 3, but `5 < 6` leaves the other class slack — so
existence of one chain cannot be the theorem; the proof must use that the complement of the chain is
also tight".  The filter case does not bite A4: the `s(5)` packing has a `45°` square, so "all
near-axis" fails, and its four near-axis squares sit exactly at `k = (T-1)^2`, where no chain is
forced.  Existence of **one** chain of `T` whose squares are axis-parallel *is* enough — A2a is two
lines (`T` centres with consecutive gaps `> 1` plus two wall offsets `>= 1/2` need `> T`).  What the
`n = 5` case actually shows is the tightness of the `(T-1)^2` threshold (objection 1) — the opposite
lesson from the one drawn, and the "complement is also tight" to-do it generates is unmotivated.

### 10. Smaller points.  **COSMETIC**

* §2 A2a is stated for chain squares "exactly axis-parallel", but A4 delivers tilt `< eps`.  Only the
  deepest stratum of `Z` is covered by A2a; everywhere else A2b/A2c are load-bearing.  The table
  should say so rather than "follows from A4 where A4 exists".
* A2b is stated for sign-coherent tilts and A2c on "the coherent cone", and the incoherent case is
  covered only by a measurement at `n = 6` (`S6_LOCAL.md` §5 (1), `-eps/4`).  No lemma in the table
  states it.  Add it or cite the measurement as an assumption.
* A3's conclusion `delta* < 0` (strict) is false on the nose: `S6_LOCAL.md` §5 (2) gives
  `delta* = 0` with `k = 4` at `T = 3` (see objection 3).  The composition only needs `<= 0`; weaken
  the statement.
* `tan(t/2) >= 1/T` is the **width** condition only; the row also needs `T sin t + cos t <= T`, which
  is slack at `t*` for every `T >= 3` (`3.0` against `17` at `T = 17`).  The stated window is right;
  the derivation should mention the second condition.  Verified: `T cos t* + sin t* = T` exactly, and
  `t* = 53.1301°, 36.8699°, 28.0725°, 22.6199°, 6.7329°` at `T = 2,3,4,5,17` — the note's numbers.
* A plateau observation worth adding to §3.2: `RANK8.md` §0's plateau reaches `32.52°` of tilt at
  `T = 4`, just **outside** the straight-chain window `28.07°`.  If that is not a coincidence it is a
  prediction — no chain certificate can explain zero margin past `2 arctan(1/T)` — and it is
  checkable on the existing plateau samples.
* A chain lemma cannot be proved by projecting onto the chain's own normal.  The admissible centre
  box `[p, T-p]^2` projects onto a direction at angle `t` with length `(T - 2p)u`, `u = cos t + sin t`,
  and `(T-2p)u - (T-1) = (u-1)(T-u-1) > 0` for every `t > 0` and `T >= 3`: `+0.0173` at `T = 4`,
  `t = 0.5°`, `+0.160` at `5°`, `+1.243` at `T = 17`, `5°` — two orders above the `-(1/3)t^2` being
  measured, and **growing linearly in `T`**.  Transverse containment is not optional in A2b, and the
  room it has to give grows with `T`.  This is the cleanest one-line statement of where `T` enters
  the chain side.

---

## 2. Direct answers to the brief

**(1) Does S0+S1+A2..A5 imply Claim(T)?  Which `theta` is uncovered?**  No.  With
`k = #{i : tilt_i < eps}` the table covers `k <= T-1` (A3) and `k = n` (A4 as written), or
`k > (T-1)^2` (A4′, objection 2).  **Uncovered: `T <= k <= (T-1)^2` with no wall-to-wall chain of `T`
among the near-axis squares** — at `T = 4`, `k in {4,...,9}`.  That region contains the `j = 4`
stratum of `Z_4` (eight free angles), the `s(5)`-type configurations at `T = 3`, and Cleemann at
`T = 17`.  Secondary gaps: A2a does not reach the near-axis (non-zero) chains A4 delivers (objection
10); the sign-incoherent case has no lemma (objection 10); A5 has an unstated lower bound on `eps`
(objection 8).

**(2) The `[new, unchecked]` items.**  S0's dilation is **correct** under closed semantics: with a
packing in side `t < T`, dilating centres by `T/t` about the origin gives wall slack
`(T/t - 1)p >= (T/t - 1)/2` and pair slack `(T/t - 1)(1/2 + W/2) >= T/t - 1`, so
`delta* >= (T/t - 1)/2 > 0`; reproduced numerically on the `s(5)` packing to 15 digits.  A4: the
acyclicity and "every pair `x`- or `y`-separated" claims are **true**, but the normal-classification
route needs `eps = O(1/T^2)`; the global-axis route (objection 5) gives the clean constant
**`eps < 1/(T-1)`**.  The `h_x = T-1` case is correct (`h_x in {T-1, T}` and with `h_x = T-1` every
Mirsky class is exactly `T`).  `tan(t/2) >= 1/T` is **correct** (width condition; the height
condition is slack for `T >= 3`).  A2b's inequality is **false** as stated (objection 6); A2c is
**false** as stated (objection 7).

**(3) The `T`-filter.**  `S0, S1, S2`, `K1/K2`, the `(T-1)^2` corollary, the total-band inequality,
A4/A4′, A2a: all uniformly true, all true at `T = 17` and at `n = T^2-4, T = 3`.  A2b/A2c are true
only inside a radius that shrinks like `2/T`.  **The only lemma that can be false at `T = 17` is the
one covering `T <= k <= (T-1)^2`** — the hole — **or A3 if the split is pushed to `(T-1)^2`.**  The
note's guess (β) is **not defensible**: the statement it names is false at `T = 4` (nine-square
witness), and A4′ is provably true at `T = 17` (it forces `>= 16` of Cleemann's squares to be
tilted, which the construction supplies).  My own guess is (δ) of objection 8 — the `eps`-window
closes — with (γ)/the hole as the visible symptom.

**(4) Is A3 plausible?  Does a near-axis mass cap cut the smeared tiling?**  *With cap `T-1 = 3`:
yes, it cuts it* — the certified optimum carries `7.60` at exactly `theta = 0`
(`notes/review-2026-09-13b.md`), every measured zero-margin configuration at `T = 4` is held shut by
an **axis-parallel** wall-to-wall row of four or an axis-parallel pinwheel core (`CENSUS.md`,
`S6_SKELETON.md` §1.2: "no core ever needed a tilted square"), and `review-2026-09-13b` lists tilt
bands among the few hypotheses that do cut the core.  *But cap `T-1` is not a cap the composition can
justify* (objection 1), and *with the honest cap `(T-1)^2 = 9` it cuts nothing* (`7.60 < 9`), and
*at that cap the scheme is refuted by explicit zero-margin witnesses* (objection 3).  The honest
summary for the note: **A3 is plausible only at a cap A4 cannot deliver.**  That is the pincer the
architecture is in, and it should be stated in §0 rather than discovered again later.

---

## 3. What the note should change

1. §2: replace `k <= T-1` / `k >= T` by `k <= (T-1)^2` / `k > (T-1)^2`, state A4′ (far squares
   arbitrary — it is free), and re-draw the hole as `T <= k <= (T-1)^2`, "the chain is there
   (measured) and no counting finds it".
2. §2 A2b: replace the first-order inequality by the exact width `T cos t + sin t + (links)` and
   name its validity window `t < 2 arctan(1/T)`.  §2 A2c: quantify over the **chain's** tilts.
3. §3.1: A4's row is `T`-uniform *and true at `T = 17`*; move the `T`-carrier to the hole/A3.
   Correct the `n = T^2-4` row (objection 9).
4. §3.3: add (δ) and replace the Cleemann test by "count tilted squares against `n - (T-1)^2 = 16`",
   plus the packing's slack `17 - t` versus `eps` (a global rotation of a slack packing is a free
   `k = 0` witness).
5. §4/§5: the decisive cheap experiments are now (a) the capped far-field LP at **cap `(T-1)^2`**,
   `T = 3` — I predict `>= n` at every `eps`, by objection 3; (b) `eps_min(T)` versus `1/(T-1)`
   (objection 8); (c) a falsification search inside the corrected hole: `T <= k <= (T-1)^2` near-axis
   squares with **no** chain of `T`, the rest tilted, maximise `delta*`.  (c) is the one that either
   finds the `T = 4` analogue of Cleemann or tells you what the missing lemma must assume.

---

## 4. Reproduce the numeric checks

    # straight-chain window, A2b over-claim, projection slack, cos+sin <= 1+phi, eps windows
    python3 - <<'PY'
    import math
    for T in (2,3,4,5,17):
        t=2*math.atan(1/T)
        print(T, math.degrees(t), T*math.cos(t)+math.sin(t))            # = T exactly
    for T,td in ((4,1),(4,28.07),(17,5)):
        t=math.radians(td); print(T,td, T*math.cos(t)+math.sin(t), T+t) # true < claimed
    for T,td in ((4,0.5),(4,5),(17,5)):
        t=math.radians(td); u=math.cos(t)+math.sin(t)
        print(T,td,(T-u)*u-(T-1))                                       # projection slack > 0
    for T in (3,4,5,17):
        print(T, math.degrees(1/(T-1)), math.degrees(2*math.atan(1/T))) # eps < window
    PY

The `s(5)` witness (four axis-parallel at `x,y in {1/2, 1+1/sqrt2}`, one at `45°` at the centre,
side `2+1/sqrt2`): min pair gap `0`, min wall slack `0`; dilated by `3/s` into `[0,3]^2` it has
`delta = 0.05409709377719363`, against S0's `(3/s - 1)/2 = 0.05409709377719396`.  The `T = 4`
nine-square witness (centres `{0.6, 2.0, 3.4}^2`): min pair gap `0.4`, min wall slack `0.1`, three
`x`-levels.
