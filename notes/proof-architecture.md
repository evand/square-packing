# Proof architecture for `s(12) = 4`: "if we had these lemmas, the theorem follows" — and where each lemma stops being true in `T`

*2026-09-20 (night).  A sketch, not a proof.  Purpose: check that the pieces the TODO is working on
compose into the target, find the lemmas nobody is working on, and track — per lemma — the container
side `T` at which it must fail.  Sources are quoted; anything marked **[new, unchecked]** was derived
while writing this note and has been verified by nobody.*

Notation.  `T` = container side, `n = T^2 - T` squares, `theta in Theta_T = [0, 90°)^n` the angles,
`delta*_T(theta)` the value of the disjunctive LP in the centres (`search/S6_SKELETON.md` §3.1),
`Z_T = {delta* = 0}`, `P_T = {delta* > 0}`.  "Tilt" of a square is its distance to `0 mod 90°`.
**Claim(T)** is `P_T = empty`, i.e. `s(T^2 - T) = T`.  Known: Claim(2) classical, Claim(3) true (Kearney–Shiu 2002);
Claim(4) is the target; Claim(5) conjectured (Stromquist: "very unlikely" to fail at 12 or 20,
`notes/proof-anatomy.md` §5.1); Claim(17) **false** (Cleemann, Friedman DS7 Fig. 8); asymptotics force
failure for all large `T`.  The first failing `T*` is in `5..11` (§0a items 9, 11: `T = 11`, `12` verified here exactly).

## 0a. Corrections after adversarial review (`notes/proof-architecture-review.md`, same night) — read first

1. **The handover threshold was wrong.**  Chain counting forces a chain of `T` among `k` near-axis squares only when
   `k > (T-1)^2`, not `k >= T` (nine squares on `{0.6, 2.0, 3.4}^2` have `delta = 0.1` and no chain of four; the optimal
   `s(5)` packing is the tight witness at `T = 3`).  In exchange A4 then holds **with the far squares arbitrary**, for every
   `T`, with `eps < 1/(T-1)` by ordering along the global axes.  A4 is a theorem-shaped, `T`-uniform, provable lemma.
2. **The hole is `T <= k <= (T-1)^2`** (`k in {4..9}` at `T = 4`), and it is not an artefact: the generic stratum of `Z_4`
   lives there (a chain of four axis-parallel squares beside eight squares tilted up to `12.8°`, margin exactly `0`).
   No capped cover LP can close it: for cap `>= 4` and `eps <= 12.8°` a real closed packing is feasible, so the value is
   `>= 12`.  A3 at cap `T - 1` is (measured) *true* and the cap that cuts the smeared tiling; A3 at cap `(T-1)^2` is what
   A4 can justify and is slack (`7.60` at `theta = 0`).  **That pincer is the state of the proof.**
3. **Guess (β) in §3.3 named the wrong lemma.**  "`>= T` axis-parallel squares contain a chain of `T`" is false already at
   `T = 4`; corrected A4 is true at `T = 17` (it forces `>= 272 - 256 = 16` of Cleemann's squares to be tilted).  The
   `T`-dependence is in the hole or in A3, never in A4.  Reviewer's proposal, **mechanism (δ)**: A3 needs
   `eps >= eps_min(T)`, A4 needs `eps < 1/(T-1)`; the window closes as `T` grows.  `eps_min(T)` is measurable.
4. A2b as written over-claims by `T(1 - cos t) ~ T t^2 / 2` (`6.1e-4` at `T = 4`, `1°`, against margins of `3.76e-6`):
   that term is where `T` enters the chain side, and it is the same `2/T` as the straight-chain window of §3.2.  A2c must
   quantify over the *chain's* tilts only.  A5 needs a lower bound on `eps` as well as an upper one.
5. Verified sound by the reviewer: S0 (dilation, numerically to 15 digits), S1, S2, the `(T-1)^2` corollary, `tan(t/2) >= 1/T`.

6. **`c_T` does not change sign** (`search/ARCH_TLEDGER.md`).  `c_5 = 3/8` exactly on its leaf
   (`-3t^2/8 + 5t^3/32 + ...`), and an exhaustive scan of all 6,814 dihedral classes of hole set finds nothing better, so
   tonight's `0.3745` was not under-optimised.  Conjecture (two exact points per branch; `T = 6` extrapolates to `0.39953` against `2/5`, from feasible points at `1.6°`, `3.2°` — an upper bound on `c_6`, two-point fit):
   `c_T = (T-2)/(2(T-1))` for `T >= 4`, `T/(2(T-1))` for `T = 2, 3`; increasing to `1/2`.  §3.2's "extrapolates through
   zero" read two branches as one curve.  The optimum is the tiling minus a *permutation* hole set; the chain dual weight is
   `1/(T-1)`.  Scenario (α) — birth at the tiling — is out: **the tiling is a strict second-order local maximum at every `T`**
   if the conjecture holds, so a local theorem at `theta = 0`, however sharp, carries none of the `T`-dependence.
7. **One exact chain line** (checked here by hand): a wall-to-wall chain of `T` squares at common tilt `t` whose centres
   rise by `R` (container frame) has `delta <= [(T - u) cos t + R sin t]/(T-1) - 1`, `u = cos t + sin t`, which is `<= 0`
   iff `R <= T tan(t/2) + cos t - sin t`.  The straight-chain window `tan(t/2) >= 1/T` is the case `R = (T-1) sin t` (equality checked at `T = 4`, `t = 28.07°`).
   Per `ARCH_TLEDGER.md` §2: `delta <= 0` as soon as some wall-to-wall chain rises by at most `1`, and the bound is lost
   once every chain rises by `2` — the `T`-dependence on the chain side is in the *rise*, an integer-like quantity.
8. **Cleemann lives in the hole.**  199 of 272 squares axis-parallel (`T = 17 <= 199 <= 256 = (T-1)^2`); a bent band at
   `± arctan(8/15)` with three `45°` squares at the elbow cuts every row and every column, so there is no wall-to-wall
   chain of 17.  `arctan(8/15) = 2 arctan(1/4)` exactly: the band is built from chains of **four** tilted squares at the
   `T = 4` window angle.  The ledger agent filed this under "(β), A4 carries `T`" against the uncorrected A4; with item 1
   it reads: **the `T`-dependence is in the region `T <= k <= (T-1)^2`**, which at `T = 4` is `k in {4..9}`.  Both agents
   and the sketch now agree on where the theorem is.
9. Smallest known counterexample is `T = 11` per the ledger agent's search (Cantrell 2025, unrefereed; refereed `T = 12`,
   Arslanov et al. 2021) — so `T*` is in `5..11`, not `5..17`.  Not independently checked here.

10. **A3 by a capped degree-1 LP is refuted at small `eps`** (`search/ARCH_FARFIELD.md`; packing side, exact certificates,
   points-only family).  `V(4, eps, cap 3) >= 12.098614247` at `eps = 0.5°, 1°` and `>= 12.0855` at `2°`, against
   `12.101126` uncapped on the same pool: the cap costs `0.0025`.  Mechanism: the LP moves the capped mass to tilt just
   above `eps` (heaviest pose a corner square at `1.2°`) — tilt-smearing is continuous, so a cap at `eps` is evaded at
   `eps+`.  Same at `T = 3`, where Claim is a theorem: `V(3, eps, cap 2) >= 6.08 … 6.017` for `eps = 0.5° … 5°`, and the
   cap only starts to bite past `10°`.  Cap `9` is never binding.  Not settled: `eps >= 5°` at `T = 4` (pool-limited: the
   pool's own uncapped value is `11.905`), and the strongest family (corner `k = 4` + polygons + regions + chord, which sits
   at exactly `12` uncapped) was **not** tested.  No cover-side number exists.
   Consequence: if A3 has an LP proof at all it is at `eps_min >~ 10°`, so the tube of A2 must be `>= 10°` wide — not a
   small-angle regime — and mechanism (δ)'s window is `[>~10°, 1/(T-1) rad = 19°]` at `T = 4`: narrow, possibly empty.
   *Erratum to the brief's addendum (mine):* "cap `>= 4` must give `>= 12` because a margin-0 packing exists" is not a
   valid argument — touching squares have coverage 2 at contacts and are infeasible for the closed-semantics LP.  The
   conclusion held numerically anyway.

11. **Band-cut scan** (`search/BANDCUT_SCAN.md`).  The record `T = 12` and `T = 11` packings verify with our rows in exact
   rational arithmetic (`delta = +6.81e-4`, `+2.75e-4`): the instrument has now seen a positive, and Claim(T) can only hold
   for `T <= 10`.  The published `T >= 12` scheme is two integer blocks + two *squeezable* rectangles with
   `waste(A) + waste(B) = T` exactly.  At `T = 4` the best band-type configuration is a margin-`0` point of `Z_4` of a new kind
   (wall-to-wall `(4,1)` stack at `28.59°` + eight axis-parallel, `k = 8`) which still contains an axis-parallel chain of four;
   no `T <= 10` configuration was found with a band cutting every row and column at margin `>= 0`.
12. **Band-cut accounting** (`notes/bandcut-cost.md`).  Mixed-tilt chain inequality MT (proved): for a wall-to-wall chain
   whose links share a normal `e = (cos phi, sin phi)`, `delta <= [cos phi (T - u_bar) + sin phi R - M]/(k - 1 + 2 cos phi)`,
   `M = sum (1/2 + W(Dtheta)/2)`; mixed tilts cost first order (`|D|/2` per link), and the certificate dies when the link
   *normals* spread (`+0.018` at `1°`, `T = k = 12`) — that is how a band escapes.  **The length accounting has no threshold
   in `T`**: gain, number of chains, interface length are all linear in `T`.  In the published scheme the threshold is
   integrality (`a_A + a_B + 4 <= T`, `a >= 4`), which over-proves by one: `T = 11` is positive and outside the scheme.
13. *Caveat on item 12's "provable target" (mine, checked by hand):* "no `(3,b)` rectangle is squeezable" is **false** under
   the natural reading (both sides can shrink): the `s(5)` packing (side `2 + 1/sqrt2 = 2.7071`) plus a column of two puts
   `7 = 12 - 5` squares in `3.7071 x 2.7071`, a squeezable `(3,4)` of waste `5 = a + 2`; Stromquist's `1.9 x 3.9475` with four
   squares is a squeezable `(2,4)` of waste `4 = a + 2` (`notes/proof-anatomy.md` §5.1).  So `a_0 = 4` and "minimum waste 6"
   are facts about one paper's table, not floors, and `waste = a + 2` is not a floor either: the `s(5)` packing itself is a
   squeezable `(3,3)` of waste `4 = a + 1`, and `s(11) < 4` a squeezable `(4,4)` of waste `5 = a + 1`.  Neither new rectangle
   moves `T*` by the two-rectangle scheme (the complement would need `(m, m+1)` with `m^2 + 1` squares).  In this language
   **Claim(T) says: a squeezable `(T,T)` square wastes at least `T + 1`** — true at `T = 2, 3`, the target at `T = 4`, false at
   `T = 11` (waste `11`).  The natural family to embed it in is `w_min(a, b)`, the least waste of a squeezable `(a,b)`
   rectangle, for small `a, b`: a table of mostly classical values, with `w_min(4,4) >= 5` the one we want.

14. **The hole, measured** (`search/BANDCUT_K.md`; feasible points).  **No chain-free configuration with margin `>= 0` was
   found at `T = 3` or `T = 4`**: every `delta* = 0` configuration carries a wall-to-wall chain of `T` among its near-axis
   squares.  Chain-free maximum `= -gamma eps^p` (`eps` in rad): far field `k <= T-1`: `eps^2/4` at `T = 3`, `eps^2/5` at `T = 4`
   (`1/(T+1)`?); in the hole `~ eps^2/8` (`T = 3`), `0.1–0.2 eps^2` (`T = 4`, `k = 4..7`, `k = 4` under-explored), and only
   **`0.1004 eps^3` at `T = 4`, `k = 8, 9`** (solid) — `5e-7` at `eps = 1°`.  Best chain-free configurations are grid-like, not
   band-like (tiling minus a permutation hole set, central `2 x 2` a pinwheel at tilt `eps`).  Two riders: "chain" must be
   read at the configuration's own margin (at gap `>= 0` the statement is vacuous); and what holds a chain-free optimum shut
   is a chain of `T` **through tilted squares**.  So the lemma the data supports is not "a near-axis chain exists" but
   **"every configuration has a tight wall-to-wall chain of `T` squares, tilted or not, on which MT (item 12) gives `<= 0`"** —
   which is exactly the statement a band defeats at `T >= 11` by spreading the link normals.

The tables below are left as first written so the review can be read against them; where they disagree with §0a, §0a wins.

## 0. Verdict, up front

1. **The spine is sound and uniform in `T`.**  The reduction to `delta* <= 0` and the Farkas form hold
   for every `T`, including `T = 17` where the conclusion is false.  So the spine carries no content;
   all of it is in three lemmas below (A2, A3, A4), of which the TODO is working on a corner of one.
2. **The "local theorem" is not local.**  `Z_4` is (near the origin) the set where some four squares
   are axis-parallel *and the other eight, at arbitrary angles, still fit beside them*.  A theorem on a
   neighbourhood of `Z` is local in 4 coordinates and **global in 8**.  The `theta = 0` theorem of
   `search/S6_LOCAL.md` §5 is its deepest stratum only.  Even *on* `Z` — four angles exactly `0` —
   `delta* <= 0` is a continuum statement in the other eight angles, and nothing on the TODO addresses it.
3. **The far field has a candidate the repo already knows how to certify** — the cover LP with a cap on
   near-axis mass (A3 below) **[new, unchecked, cheap to test]** — which would replace the 12-dimensional
   angle-box tree that `S6_SKELETON.md` §6 rightly calls infeasible.
4. **On "where does it break in `T`"**: every exact lemma the repo has is uniform in `T` and critical
   exactly at `n = T^2 - T` (§3).  The measured second-order constant at the tiling is
   `c_T = 1, 3/4, 1/3` at `T = 2, 3, 4` and `<= 0.375` at `T = 5` (one run tonight, not trusted: see
   §3.2).  It has not gone negative, so the break is not yet visible at the tiling; §3.3 lists the
   three places it can be born and the measurement that distinguishes them.

## 1. The spine (true for every `T`)

| | statement | status |
|---|---|---|
| **S0** | Claim(T) `<=>` `delta*_T(theta) <= 0` for all `theta`.  (`=>`: a point with `delta > 0` sits in `[delta, T - delta]^2`.  `<=`: dilate the centres of a packing in side `t < T` by `T/t`; containment and every separating-axis gap gain strictly.) | exact; `S6_SKELETON.md` §3.1 has the fixed-angle half.  Lean-able in a day. |
| **S1** | `delta* = max_sigma LP_sigma(theta)` over finitely many assignments `sigma` (one separating normal per pair); `LP_sigma(theta) <= 0` iff there is `w >= 0` with `w^T A_sigma(theta) = e_delta`, `w^T b_sigma(theta) <= 0`. | exact (LP duality). |
| **S2** | Closed semantics: "pairwise disjoint as closed sets in `[0,T]^2`" is what `delta* > 0` means; everything at `t = T` must survive touching (`search/RANK8.md` §7). | `notes/branch-semantics.md`. |

Nothing here knows what `T` is.  A Lean skeleton should consist of S0, S1 and the *statements* of §2 with
`sorry`, plus the one-line composition; do not formalise inside any lemma until the DAG below is stable.

## 2. Shape A: tube around `Z` + far field

Fix `eps > 0`.  Call a square *near-axis* if its tilt is `< eps`.  Split `Theta_T` by the number `k` of
near-axis squares.

| | lemma | what it needs | status |
|---|---|---|---|
| **A3** (far field) | `k <= T - 1`  `=>`  `delta* < 0`. | Candidate **[new]**: the `t = T` cover LP (points + polygons, the machinery of `search/BENTZ.md` §7) with the extra row "total mass on poses of tilt `< eps` is `<= T - 1`" has value `< n`.  Evidence it might: axis-parallel alone is worth exactly `9`, and `>= 2.6` of the missing `3` is bought under `1°` (`notes/status.md`); the optimum puts `7.60` of its mass at exactly `theta = 0`.  Capping that mass at `3` removes the smeared tiling, which is the only thing holding every degree-1 family at `12`. | **untested.**  One LP family in `eps`; certifiable cover-side with the existing exact pipeline.  If the value is `>= 12` for every `eps`, shape A's far field has no candidate and this note should say so. |
| **A4** (chain exists) | `k >= T` and all near-axis `=>` some `T` of the squares form a wall-to-wall chain along `x` or along `y`. | **[new, unchecked]** For tilts `< eps ~ 1/T` every separating normal is within `eps` of an axis; "x-separated" is acyclic; Mirsky on its transitive closure: longest `x`-chain `h_x <= T`, and the squares split into `h_x` classes of pairwise `y`-separated squares, each a `y`-chain of `<= T`.  With `n = T(T-1)`: either a chain of `T` exists in one direction, **or `h_x = T - 1` and every class is a full `y`-chain of `T`** — either way a chain of `T`, and in the second case `T - 1` disjoint ones covering everything.  For mixed tilts (some squares at `45°`) normals do not classify and this argument is silent. | combinatorics looks routine for all-near-axis; **no version exists when `T <= k < n`**, which is most of `Z`. |
| **A2a** (on `Z`) | The `T` squares of a chain exactly axis-parallel, the rest arbitrary `=>` `delta* <= 0`. | Trivial *given* the chain: `T` unit widths in width `T`.  The content is A4. | follows from A4 where A4 exists. |
| **A2b** (near `Z`, first order) | Chain squares tilted by `s_1..s_T`, small: chain length `>= T + (1/2) sum |Dphi| + (1/2)|phi_wall| - (staircase credit)`; the credit needs sign-coherent tilts and transverse room. | `search/S6_LOCAL.md` §5 (i): finite LP/MILP over the axis-parallel plateau.  `RANK8.md` §3.2 has it at the `T = 4` tiling. | **TODO item 1(i).**  Well-posed.  But stated at `theta = 0` only; A-shape needs it *at every point of `Z`*, with the eight free squares as parameters. |
| **A2c** (near `Z`, second order) | On the coherent cone, `delta* <= -c_T |tilt|^2` with an explicit radius. | Model case: closed forms of `S6_LOCAL.md` §3, which are exact in `t` on their leaf, **not asymptotic** — `-(u-1)^2` at `T = 2`, `-3(u-1)^2/(u^2+3)` at `T = 3`, both `<= 0` for every `t`. | **TODO item 1(ii).**  Radius must be explicit or A3 cannot meet it. |
| **A5** (glue) | `eps` in A3 `<=` radius in A2b/c, uniformly over `Z`. | needs constants from both sides. | nothing yet. |

**Composition.**  S0 + S1 + [A3 for `k < T`] + [A4 + A2a/b/c for `k >= T`] + A5 `=>` Claim(T).  The gap
in the composition, visible only when written down: **A4 is stated for "all near-axis" but A3 hands over
everything with `k >= T`**, including configurations with four near-axis squares and eight at `45°`.
Those are not near `Z` (the instrument says `delta* ~ -0.1` there) but no lemma in the table kills them.
Either A3 must be strengthened (cap the near-axis mass *per row/column strip*, a pose-region
hypothesis the LP can express: "no wall-to-wall strip holds four near-axis squares"), or A4 must be proved
with the far squares as obstacles.  **This is the hole.**  It is the `T = 4` face of `RANK8.md` §1's
"plateau reaches `32°`".

## 2′. Shape B: one chain theorem, no tube

Suggested by the closed forms being global in `t`: for every assignment `sigma` there is a chain-type
dual `w_sigma(theta)`, trig-rational in `theta`, valid on a whole semi-algebraic piece of `Theta`.
Then there is no `eps`, no glue, and the proof is "every assignment contains a chain; every chain
certificate is `<= 0`".  Against: `S6_LOCAL.md` §5 — the single-chain inequality still contains the
transverse offsets `Dy`, and eliminating them is the LP dual leaf by leaf; at `theta = 0`, `n = 6` there
are 2,284 tight leaves, and at `45°` a kill costs 53,445 nodes with no chain in sight.  Shape B is
plausible only on `k >= T`; it is a way to *do* A4 + A2, not a replacement for A3.  So the two shapes
are the same shape, and the far field needs its own idea either way.

## 3. The `T`-ledger: which lemma carries "`T < T*`"?

### 3.1 What is uniform

| ingredient | `T`-dependence | note |
|---|---|---|
| S0, S1, S2 | none | |
| K1, K2 (`S6_SKELETON.md` §2) | one hypothesis, `k lambda >= T`; silent below `k = T` | passes the filter by being silent |
| axis-parallel count `(T-1)^2` and total-band inequality `sum D(theta_i) < T(T-1)/2` | uniform, **exactly critical at `n = T^2 - T` for every `T`** | says "not all axis-parallel" and nothing more, at `T = 3, 4, 17` alike |
| A4's counting | uniform; uses `n > (T-1)^2` for existence and `n = T(T-1)` for "all chains tight" | at `n = T^2 - 4`, `T = 3` (`s(5) < 3`): `5 > 4` still forces a chain of 3, but `5 < 6` leaves the other class slack — so **existence of one chain cannot be the theorem**; the proof must use that the complement of the chain is also tight |
| first-order rigidity (links pay `|Dphi|/2`, anchors `|phi|/2`) | mechanism is `T`-free | measured at `T = 3, 4` only |

So every exact statement in the repo is `T`-uniform.  A proof assembled only from these is wrong.

### 3.2 What is not uniform (the places `T` enters as a number)

* **The straight-chain window.**  `T` squares at a common tilt `t` in a straight row need width
  `T cos t + sin t`, which is `<= T` iff `tan(t/2) >= 1/T` **[new; elementary]**: `t >= 53.1°, 36.9°,
  28.1°, 22.6°, 6.7°` at `T = 2, 3, 4, 5, 17`.  At `T = 2` the window is empty (`> 45°`).  From `T = 3`
  on, a tilted row fits *by itself*; what stops it is stacking.  The tube radius in A2 cannot exceed
  `~2/T`, and the tilt needed to beat a chain shrinks like `2/T`.
* **The second-order constant at the tiling**, uniform tilt, `delta* = -c_T t^2 + ...`:

  | `T` | 2 | 3 | 4 | 5 |
  |---|---|---|---|---|
  | `c_T` | `1` (exact, `-(u-1)^2`) | `3/4` (exact on its leaf) | `1/3` (exact on its leaf) | `<= 0.3745` (`runs/s6local_n20_T5_uniform.txt`, 6,000 of 53,130 subsets, hit by `91/18,004`; a feasible point, so the true `c_5` can only be smaller) |

  `1, 3/4, 1/3` extrapolates through zero before `T = 5`; the measured `T = 5` value is *above* `T = 4`,
  which is either a parity effect or the same under-optimisation that produced `-2/3` for `-1/3` at
  `n = 12` (`S6_LOCAL.md` §0).  A re-run with 60,000 starts is in `runs/s6local_n20_T5_uniform_seed7.txt`.
  **A general-`T` closed form for the pinwheel leaf (`s6exact.py`) is the single most informative
  computation for this section**: if `c_T` changes sign at some `T_loc`, the tiling stops being a local
  maximum there and A2c is the lemma that carries `T < T_loc`.
* **The far-field LP value** (A3) as a function of `T` — not measured at any `T`.

### 3.3 Three ways `P_T` can be born, and how to tell

`delta*` is continuous, so as `T` grows `P_T` appears either (α) **at the tiling**: `c_T` changes sign;
A2c fails, shape A's skeleton stays valid and the lemma carries the bound.  (β) **elsewhere on `Z`**:
a stratum with some near-axis chains and a tilted block; A4 fails — the chain lemma "leaks".  (γ) **away
from `Z`**: A3 fails.  Test: take Cleemann's 272-square packing (DS7 Fig. 8), read off its angle vector,
and classify: how many squares are axis-parallel (if `>> 17`, as in most Friedman-survey packings —
axis-parallel blocks cut by a tilted band — then it is (β), and the statement that dies is "`>= T`
axis-parallel squares contain a wall-to-wall chain of `T`": the tilted band cuts every row and column).
**My guess is (β)**, which would mean the `T`-dependence sits in A4 — the lemma with no owner — and not
in the local theorem at the tiling.

## 4. What this does to the TODO

| TODO item | serves | verdict |
|---|---|---|
| local theorem at `n = 6`, (i) first order | A2b at the deepest stratum | keep, but write it `(n, T)`-generic and **parametrised by a base point of `Z`**, not only `theta = 0` |
| (ii) second order on the cone | A2c | keep; add "general-`T` closed form, sign of `c_T`" |
| signed-tilt redo at `n = 12` | A2b | subsumed by the first-order instrument at `T = 4` |
| **(missing)** cover LP with near-axis mass cap, `T = 3` then `T = 4` | A3 | **add; cheapest decisive experiment in this note** |
| **(missing)** chain lemma with obstacles (`T <= k < n`) | A4 | **add; first as a falsification search** — four near-axis squares not in a chain, eight tilted, `delta*` as close to `0` as possible |
| **(missing)** Cleemann classification; `c_T` closed form | §3.3 | add |
| narrow-box kill test at `n = 6` | the angle-box far field | demote: only matters if A3 fails |
| corner triage, level-2 slot leaves, rank-10 rehearsal, 16-point certification | none of A2–A5 | drop (as proposed earlier today) |

## 5. How to falsify each lemma with what exists

* A3: `search/bentz.py` / `leaf_ceiling.py` with a tilt-band mass row; at `T = 3` first, where Claim is a theorem and the far field is known to close by point lattices (Kearney–Shiu).
* A4: `s6local.py` with `--base` putting `T` squares at `0` and the rest at large tilts, maximise `delta*` subject to "no wall-to-wall chain among the near-axis squares" (needs a chain detector: `search/chains.py`).
* A2b: the 116 signed directions of `runs/s6cube_n6.jsonl` (one-sided: exact `a(v)` `<=` measured).
* A2c / §3.2: `s6exact.py` at `T = 5, 6`; `s6local.py` at `T = 5` with all 53,130 subsets.
* Filter, binding on every proposed lemma: state it for general `T`; say whether it is true at `T = 17`
  and at `n = T^2 - 4`, `T = 3`.  Uniformly true lemmas are allowed — but the set of lemmas must contain
  at least one that is false at `T = 17`, and the note must say which.
