# `s(12) = 4`: the negative record  (2026-09-22)

What the repo has established about the `n = 12` case, for a reader deciding whether to park it.
Every number carries its source; every claim carries one of four labels:

* **[proved]** — an exact certificate re-checked in integer/rational arithmetic, a Lean theorem, or a
  hand identity re-verified by a script in the repo;
* **[measured]** — the value of an explicit configuration or a float LP: a one-sided bound (a feasible
  point can only under-report a maximum), never a theorem;
* **[heuristic]** — a fit, an extrapolation, or an argument with an unproved step;
* **[guess]** — neither.

Summaries this note was checked against: `notes/status.md` (2026-09-22) and
`notes/review-2026-09-22.md`; where a summary and its source disagree, the source wins and §7 says so.
Semantics throughout: `t = 4`, closed unit squares, closed containment, closed container `[0,4]^2`;
a *packing* is a family pairwise disjoint **as closed sets** (`notes/s13-casefree.md` §1).

---

## 1. The statement, and why closed semantics is the content

`s(n)` is the side of the smallest square holding `n` unit squares, rotations allowed.  Known:
`s(13) = 4` (Bentz 2010; re-proved case-free here by one weighted closed cover of `[0,4]^2` of weight
`2591194431/200000000 = 12.955972 < 13`, `notes/s13-casefree.md`) and `s(11) = 3.877083…` (Trump 1979).
Certified here: **`s(12) >= 15680/3951 = 3.968616`** (README; 1736-point weighted cover of weight
`11.9738036`, exact `i128` verifier, Lean reduction).  `s(12) = 4` is conjectured and open; deciding it
decides whether the largest `n` with `s(n) < 4` is `11` or `12`.

Three equivalent forms of the target, each **[proved]** equivalent to the others:

1. no 12 unit squares fit in a square of side `< 4`;
2. (dilate by `4/s'` about the centres, `notes/s13-casefree.md` §1, `notes/proof-anatomy.md` §7.1)
   **no 12 closed unit squares in `[0,4]^2` are pairwise disjoint as closed sets** — touching is
   forbidden, so this is a statement at margin zero;
3. (S0, `notes/proof-architecture.md` §1; the fixed-angle half in `search/S6_SKELETON.md` §3.1)
   `delta*(theta) <= 0` for every angle vector `theta in [0, 90°)^12`, where `delta*(theta)` is the value
   of a **disjunctive LP in the 24 centre coordinates** (containment linear; each of the 66 pairs a
   disjunction of separating-axis rows with half-width `m(D) = 1/2 + (|cos D| + |sin D|)/2`).

The `4 x 4` tiling has margin exactly `0`, and so does every tiling-minus-4 — twelve closed unit squares
in `[0,4]^2` with pairwise disjoint interiors and `18` of `66` pairs touching (`search/RANK8.md` §1.2,
exact).  So the object to be excluded sits at the boundary of the feasible set, and every method below
fails in the same place: at margin `0`, under an infinitesimal tilt.  A cover-type certificate of
`[0,4]^2` in *open* semantics costs `>= 16` (the 16 eroded tiling squares are disjoint); in closed
semantics the 16 tiling squares share boundary points and `12.956 < 13` is attainable
(`notes/proof-anatomy.md` §7.1; the case-free `s(13)` cover realises it).

---

## 2. What is proved about `n = 12` at `t = 4`

### 2.1 Every degree-1 certifiable family is `>= 12` at the container

All rows are lower bounds on the value of the named LP relaxation, given by an exact rational measure
whose maximum coverage `M` and every other row were re-derived from the measure file in integers.  A
lower bound `>= 12` proves the family **cannot** certify `s(12) = 4`; a lower bound `< 12` proves nothing
(`search/T4SCREEN.md` §0 calibration, `search/BENTZ.md` §8).

| family (`t = 4`, closed) | certified lower bound | source | label |
|---|---|---|---|
| pure, points only | `nu_f >= 24537607710/1999999999 = 12.2688038611` (294 poses, 2352 squares, `M = 1999999999/2000000000`) | `search/COVER4.md`, `search/cover4_exact_support.txt` | [proved] |
| points + odd polygons (2,002 rows) | `12173398461/10^9 = 12.173398` (`M = 62499993/62500000 < 1`) | `search/pgonly_pure_exact.txt`, `notes/review-2026-09-12.md` | [proved] for the relaxation *with these rows*; "no polygon goes below 12" is [heuristic] (separator is heuristic) |
| corner `k = 3` + polygons + regions + chord | `3006818329/250000000 = 12.027273` (LP `12.038` still rising; corner masses `0.9991` by a `1e-3` boundary tie-break) | `search/BENTZ.md` §7.2 | [proved] |
| corner `k = 4` + polygons + regions + chord | `5999999963/500000000 = 11.999999926`; improved to `2999999989/250000000 = 11.999999956`; the LP value is `12.000000000` on every pose set ever tried and has never exceeded it | `search/pgonly_corner_exact.txt`, `search/FARFIELD_STRONG.md` §3 | lower bound [proved]; "`= 12`" [measured] |
| fully pinned Bentz leaf `(AB)^4` | **`11.999999926 <= value <= 12` exactly** | `search/BENTZ.md` §0 | [proved] both ends (§2.2 item 6) |
| corner `k = 4`, `theta = 0 mod 90°` only | **exactly `9`** | `search/BENTZ.md` §0 C1 | [proved] (§2.2 item 1) |
| corner `k = 4`, `|theta mod 90°| <= 1.0000034°` | `11599390413/10^9 = 11.599390` (`M = 1`, regions OK) | `search/BENTZ.md` §7.2 | [proved] |
| corner `k = 4` family + cap `mu(tilt < eps) <= 3` | `11.999999870` at `eps = 1°`, `11.999999855` at `2°` (cap multiplier `0`) | `search/FARFIELD_STRONG.md` §4.1, §7b | [proved] |
| points only + cap `mu(tilt < eps) <= 3` | `12098614247/10^9 = 12.098614` at `eps = 0.5°, 1°`; `12.085541` at `2°` | `search/ARCH_FARFIELD.md` §0, §4 | [proved] |
| points only + cap `<= 9 = (T-1)^2` | `>= 12.101126` for `eps <= 15°` (row never binding) | `search/ARCH_FARFIELD.md` §6 | [proved] |
| level-2 slot leaf `01010101` under `k = 4`, no cliques | `11.920067` (float LP, `M = 1.000000000`), rising `+0.002–0.006` per stage when stopped | `search/T4LEAF.md` §0 | [measured] |
| any sound clique rule (positive-volume or not) | honest cover cost `>= 18.109` (`E2Pg`) / `>= 20.162` (`SUPEPGF`) — a floor for *every* sound rule | `search/HONEST.md` §0, `search/ALLMEET.md` §0 item 1 | [measured] |

Reading: pure `12.27` → polygons `12.17` → corner `k = 3` `12.03` → corner `k = 4` `12.000` → pinned leaf
`[11.999999926, 12]`.  No leaf of any count, slot or incidence tree at `t = 4` has been shown `< 12`, and
every one that converged is `>= 12` or pinned at `12`.  Chord rows have dual `0` in every `m = 4` leaf
(`search/T4LEAF.md` §3.1) and in every Bentz run (`search/BENTZ.md` §4); polygons are worth exactly `0` at
the corner leaf (`notes/review-2026-09-12.md`).

### 2.2 The theorems

1. **Axis-parallel is `<= 9`.**  [proved]  The nine points `{1,2,3}^2` are a closed cover for
   axis-parallel unit squares in `[0,4]^2` (`[a, a+1]` with `a in [0,3]` contains `1`, `2` or `3`), so any
   measure with coverage `<= 1` there has mass `<= 9`; the corner `k = 4` LP restricted to `theta = 0`
   converges to an integral nine-square packing of mass exactly `9` (`search/BENTZ.md` §0 C1, 71 s).
   Point-set-free proof: for axis-parallel squares K1's bands are unit intervals, so at most `(T-1)^2`
   fit (`search/S6_SKELETON.md` §2.2 corollary).  Consequence: all three of the missing units are
   rotation, and `11.599390 - 9 = 2.599` of them are certified to live under one degree of tilt.
2. **`COVER^closed(4) >= nu_f^closed(4) >= 12.2688038611`.**  [proved]  `search/COVER4.md`: the measure
   `mu/M` is a feasible fractional packing, and weak LP duality.  Consequence drawn in `notes/t3-chord.md`
   §6.2: **no pure unavoidable set of `[0,4]^2` has `<= 12` points**; the slack-`0` route (an `11`-point
   set) is dead.  Upper side: Friedman's `14` (DS7 Thm 4, certified by `zeromargin.py friedman14`); the
   fractional optimum is believed `~12.4` [heuristic, `search/COVER4.md` bracket].
3. **The wall-strip chord lemma.**  [proved, Lean]  `lean/Sqpack/Chord.lean` `wall_strip_le_three`:
   for `t <= 4`, unit squares in `[0,t]^2` pairwise disjoint as closed sets with centre `y <= 1` number
   at most `3` (`notes/chord-lemma.md` B1; packing form B2 for side `L > 1`).  So the row
   `mu({c_y <= 1}) <= 3` is valid for the closed-semantics relaxation.  (Several later notes still call
   this row a "flagged hypothesis", quoting `search/T4LEAF.md` §1.2 of 2026-09-07 — the day the Lean
   proof landed; see §7.)  The `1e-9` tolerance `T4LEAF.md` puts on `c_y` is argued there, not in Lean.
4. **The chain identity and its rise threshold.**  [proved]  For a wall-to-wall chain of `T` squares at
   common tilt `t` with consecutive pairs separated along `e = (cos t, sin t)` and rise `R`:
   `delta <= B_T(R,t) = [(T-u) cos t + R sin t]/(T-1) - 1`, `u = cos t + sin t`, and
   `B_T <= 0 iff R <= R*(T,t) = T tan(t/2) + cos t - sin t` (`search/ARCH_TLEDGER.md` §2).  In particular a
   chain of rise `<= 1` gives `delta <= 0` for every `t` and every `T >= 3`, equality only at `t = 0`.
   Mixed-tilt form MT (`notes/bandcut-cost.md` §1): a tilt mismatch costs `(W(D)-1)/2 ~ |D|/2` per link,
   first order; the sum telescopes only for a common normal, and a normal spread of `1°` adds `+0.0127` at
   `T = 4` (that repair term is [proved modulo a heuristic `tau` bound]).  Straight-stack window: `T`
   squares in a row at tilt `t` fit width `T` iff `tan(t/2) >= 1/T` (`28.07°` at `T = 4`); this is exactly
   Cleemann's `arctan(8/15) = 2 arctan(1/4)` (`search/ARCH_TLEDGER.md` §4–5).
5. **Chain counting (Lemma CC).**  [proved, Lean-ready]  `notes/counting-ladder.md` §1: among
   `k > (T-1)^2` squares of tilt `<= eps` in a configuration of margin `delta > -1`, some `T` form a chain
   in the threshold order `x_j - x_i >= (1+delta)/W(eps)` (or the same in `y`), with the far squares
   arbitrary; the chain supports `delta <= (T-1)(W(eps)-1)/(T-1+2W(eps))`, which is `<= 0` **only at
   `eps = 0`**.  Sharp: `(T-1)^2` axis-parallel squares at spacing `T/(T-1)` have margin `1/(2(T-1))` and
   no chain of `T`.  The tilt bound `eps < eps_max(T) = arcsin(T/(sqrt2 (T-1))) - pi/4` (`25.53°` at
   `T = 4`) is what makes the chain span the container, not what makes it exist.  With `n = 12 > 9` a
   chain of four always exists in the order form; whether it can be taken as a DAG path is open.
6. **The pinned leaf's upper bound.**  [proved]  On the Bentz leaf `(AB)^4` the twelve admitted patterns
   are distinct and each pattern region lies in `{S : p in S}` for a point `p` of its pattern, so
   `mu(R_pi) <= 1` by the coverage row at `p` and the objective is `<= 12` on any pose set
   (`search/BENTZ.md` §0).  The lower bound is the `PGCORN` measure, which lies inside the leaf.
7. **Farkas machinery in angle space.**  [proved]  Lemma H closed form
   `H(T,L,t) = [(T-u)M - (T-1) - L tan t]/[(T-1) + L tan t + 2M]`, `M = (1 + sin^2 t + sin t cos t)/cos t`,
   with `H <= 0 iff L >= L*(T,t)`, `L*(T,0) = T-2` (`notes/t3-chain.md` §1.3); Lemma W' (every
   normalised certificate has `delta <= Omega (T+1)/2 - 1 - B`, `Omega` = wall weight, `B >= 0`) and the
   wall-free master formula `delta <= [sum_i |r_i|_1 (T-u_i)/2 - sum_e w_e m_e]/[sum_e w_e + sum_i |r_i|_1]`
   which reproduces the full LP dual (`notes/t3-existence.md` §2.1–2.2); Lemma Z'' (a wall-to-wall chain
   of `T` with exactly axial link normals gives `delta <= [T - (u_1+u_T)/2 - sum m_s]/(T+1) <= 0`,
   `notes/t3-existence.md` §2.4); the corrected cycle lemma (turn cost is the `L1` jump `2 cos t`, not
   `2 sin(alpha/2)`, `notes/t3-existence.md` §2.3); the exact `45°` certificates
   `delta <= (93 - 66 sqrt2)/7` at `T = 3` and `(7 sqrt2 - 10)/2` at `T = 4` (`notes/t3-chain.md` §3.2,
   `search/T4_CYCLES.md` §4.3).
8. **The three-level ladder closed form.**  [proved]  `notes/counting-ladder.md` §2: for a crossbar of
   `T` with two adjacent tilted members at tilt `eps` placed one square in from the wall (`B = 0`),
   `L = T-2` transverse links and one rung, `H1(T, eps) = -eps^3/(2(T+1)) + O(eps^4)`; without the rung
   `+eps^2/(T+1)`; with the tilted block against the wall `+(T-3) eps^2/(T+1) > 0` at every `L`.  At
   `T = 4` the closed-form weights reproduce the LP dual of the thinnest cell digit for digit.
9. **Chord rows (Prop. CC3).**  [proved]  `notes/t3-chord.md` §3.1: three (at `T`: `T`) squares whose
   centres lie within `(cos t - sin t)/2` of one axis-parallel line have interior chords `sec t` each, so
   `delta <= -T (sec t - 1)/(T+1) = -(T/(2(T+1))) t^2 + O(t^4)` — exactly tight at `t = 0` and correct to
   second order off it, the only row in the repo with both properties.  It is a combination of existing
   pair, wall and centre-box rows, so it adds a *leaf* (centre branching), not a new valid inequality
   (`notes/t3-chord.md` §1.3).
10. **The `T = 3` far field, and its `T = 4` failure.**  [proved]  `notes/t3-existence.md` §3: centres
    of a packing are `n` points at pairwise distance `>= 1 + delta` in a square of side `T - u_min - 2 delta`,
    so `(T - u_min - 2delta) d_n >= 1 + delta`.  With Graham's `d_6 = sqrt13/6`, six unit squares all tilted
    `>= 25.8431°` do not fit in `[0,3]^2` at any margin; since Lemma H fails only above `28.7959°`, the
    existence clause has no far field left on the `T = 3` uniform-tilt line (below `25.84°` that line
    stays open, `notes/t3-existence.md` §8.1 row 4).  At `T = 4`: `4 - sqrt2 = 2.5858 > 1/d_12 = 2.5725`,
    so the mechanism never fires, at any angle, by `0.5 %`.
11. **The `T = 3` axis-parallel strata.**  [proved]  Theorem AP: six axis-parallel unit squares in
    `[0,3]^2` have `delta <= 0` and at `delta = 0` contain a wall-to-wall (staircase) chain of three; AP5:
    the same with five axis-parallel squares; at `k = 4` the count is exactly critical (`delta <= 1/3`,
    attained) (`notes/t3-existence.md` §4.2).  These, with item 10, are the only proved parts of the
    existence side of the angle-space route.
12. **Helly is false for squares at `t = 4`.**  [proved]  Three unit squares with edges on three half
    planes at `120°` pairwise meet with margin `+0.358`, thicken to positive-volume pose boxes that stay
    pairwise meeting, and admit no common point (exact Farkas certificate, `search/ALLMEET.md` §4).
13. **The counterexample frontier.**  [proved]  The record packings of `n = 110` in side `11`
    (Cantrell 2025, unrefereed) and `n = 132` in side `12` (Arslanov–Mustafin–Shangitbayev, EJC 2021,
    refereed) verify with the repo's own rows in exact arithmetic at margin `+2.748946598e-04` and
    `+6.805572405e-04` (`search/BANDCUT_SCAN.md` §2).  The tiling is a strict second-order local maximum
    at every `T`: `c_T = 1, 3/4, 1/3, 3/8` exact on the optimal leaf at `T = 2..5` (`search/ARCH_TLEDGER.md`
    §1–3), and `B_T(1, t) <= 0` for all `t`, `T >= 3` (item 4).
14. **`s(13) = 4` without case analysis.**  [proved]  One `3,621`-point closed cover of weight
    `12.955972155 < 13`, checked exhaustively at margin zero by two checkers sharing no code, `0` boxes
    uncertified in each, `23` rejection tests, Lean for every primitive (`notes/s13-casefree.md`).  It says
    nothing about `n = 12` because `12.2688 > 12` (item 2).

What is *conditionally* proved at `t = 4` about an actual 12-packing, in one line each: all
axis-parallel ⇒ at most `9` squares (item 1); at most `3` squares centred within `1` of any wall
(item 3); at `t <= 3.98`, every 12-packing has `>= 3` corner boxes occupied (verified `k = 0, 1, 2`
branch certificates, `notes/review-2026-09-13b.md`).  Nothing else.

---

## 3. What is measured

Everything here is a multistart or a float LP: a feasible point, hence a one-sided bound, never a
theorem.  Reliability is stated where the source states it.

### 3.1 The missing unit, localised

* **Integrality gap exactly one square, on the optimum's own support.**  [measured, exact on the
  support]  The 162 poses of the certified corner-`k = 4` measure carry mass `11.999999926`; their exact
  disjointness graph (`7,644` of `13,041` pairs disjoint as closed sets) has independence number
  **`alpha = 11`** by complete branch-and-bound, reproduced by `rankdiag.py`; no violated `C5`, `C7`,
  odd antihole or 5-wheel (`search/BENTZ.md` §5 (iii)).  The four corner pairs, the four `{c_j}` and the
  four `{d_j}` are each realisable alone (`4/4`); **the eight non-corner singletons together admit
  `7`**.  The statement a proof needs is therefore rank 8: with four corner squares holding `{a_i, b_i}`,
  at most seven further squares each contain one of `c_0..c_3, d_0..d_3`.  True integrally (one wall
  square per wall between pinned corners, then `s(4) = 2`), and not expressible by any single-square
  inequality [heuristic: no degree-1 row the repo has cuts it, `search/BENTZ.md` §6].
* **Its margin is zero, attained, on a plateau.**  [measured]  `search/RANK8.md` §1: over `3,016`
  starts the sup of the min pairwise closed gap on leaf A is `+0.000000000000e+00`; `2,188` starts reach
  `0` to `1e-13`, `1,765` of them further than `0.05` from any tiling-minus-4 family, the furthest at
  `0.914`, with tilts up to **`32.5242°`**.  The lower bound `>= 0` is exact (the 16 tiling-minus-4
  families, `18/66` touching pairs).  So leaf A is realisable under open semantics and fails under closed
  semantics by touching only, and the maximisers do not converge to the tiling.
* **No smaller core.**  [proved for every "realisable"; measured for the one "not"]  All `51` `D4`
  classes of the eight singletons and all `618` classes over twelve labels: every sub-configuration of
  eleven or fewer squares is a genuine closed packing with an exact rational witness (`50/51`, and all
  three eleven-square classes, re-verified in integers).  Margins: drop an interior singleton
  `4.426e-02`, a wall singleton `8.100e-03`, **a corner square `3.763847e-06`** (a three-square block
  rotated `-3.2568°`, twelve equal gaps).  The unique minimal non-realisable sub-configuration is leaf A
  itself (`search/RANK8.md` §2).
* **First-order rigid; second order unsigned.**  [measured]  The linearised system around every one of
  the sixteen tiling-minus-4 families has MILP value `0`; an irreducible 18-row certificate reads as
  three wall-to-wall chains of four (`T2, D1, D2, C3` along `x`; `C1, D0, D1, C2` along `y`; a five-square
  column).  The linearisation is combinatorially valid for `rho + 0.7072 alpha < 0.086`, with remainder
  `<= 2 rho alpha + (1 + 2 rho) alpha^2/2 + alpha^2` — `1.2e-3` at `1°`, three orders above the `3.76e-6`
  margin it would need to see, and of positive leading sign (`search/RANK8.md` §3).
* **Sherali–Adams level 2 sees a tenth.**  [measured]  On the 162-pose support, pair marginals on
  disjoint pairs with the twelve clique rows lifted: **`11.893347407`** against the clique LP's `12` and
  `alpha = 11`; `7.894639816` on the eight singletons against `8` and `7`.  The dual lives on the
  tile-contact graph of the eight non-corner tiles only (four edge-adjacent interior pairs, eight
  wall–interior pairs); zero on the diagonal pairs `(D0,D3)`, `(D1,D2)` and on wall–wall pairs
  (`search/RANK8.md` §4).

### 3.2 The census of pattern leaves

[measured, `search/CENSUS.md`]  Over all `10,945` `D4` classes of families of twelve pairwise-disjoint
realisable patterns on Bentz's 16 points (`86,403` leaves): **no class has positive margin**; leaf A is one
of **`>= 2,779`** zero-margin classes — `1,425` containing an exact tiling-minus-4 (`274,269` (tiling,
assignment) pairs enumerated, seeds re-checked), and `1,297 + 57` at `0` (to `1e-14`) with no tiling in
them, tilted up to `44°`, with patterns no tile holds.  Pass-1 negatives in `[-0.2, -0.02)` (1,601
classes) are not to be trusted (57 of 65 refined cells moved to exactly `0`).  Of `3,052` disjoint
pattern pairs exactly one `D4` orbit is pairwise infeasible (`-0.040244`): negative leaves die through
three or more squares.

**What holds a plateau point shut:** the first-order core is one straight wall-to-wall row of four
axis-parallel squares in `2,137/2,145` leaf-A plateau points and `2,677/2,722` census optima, the rest
corner-contact "pinwheel" cores at exactly axis-parallel configurations; **no core ever needed a tilted
square**; no first-order-feasible plateau point (`search/CENSUS.md` §2).  So the point-pattern tree is
blind to the obstruction, which lives in the separation structure.

### 3.3 The `s(6) = 3` rehearsal and the deficits off the zero set

[measured, `search/S6_SKELETON.md`, `search/S6_LOCAL.md`, `search/ARCH_TLEDGER.md`]

* `n = 6`, `T = 3`: margin exactly `0` on a plateau reaching `45°` of tilt (`25,364` local optimisations,
  `20,730` at `0`, none positive); first-order core a wall-to-wall row of three in `1,930/1,959`
  (`98.5 %`), the rest pinwheels; none tilted.  Zero set in angle space: `Z ⊂ {>= 3 angles ≡ 0 mod 90°}`
  (`0/80` at `j <= 2`, `3/80` at `j = 3`, `49/80` at `j = 4`, `80/80` at `j >= 5`); generic angle vectors
  are strictly negative (max `-0.0296` of 40).
* Deficit by number `j` of axis-parallel squares, the rest at a common tilt `t` of one sign:
  `n = 6`: `-(3/4) t^2` (`j = 0`, the pinwheel; on its leaf exactly `-3(u-1)^2/(u^2+3)`), `-(1/2) t^2`
  (`j = 1, 2`), **`-(1/4) t^3`** (`j = 3`; exact on its leaf), `0` (`j >= 4`).  `n = 12`:
  **`-(1/3) t^2`, identical to 10 digits for `j = 0..3`** (on its leaf `-(1/3)t^2 + (2/9)t^3 + (1/54)t^4`),
  exactly `0` for `j >= 4` up to `12.8°` — eight squares tilted to `12.8°` beside a wall-to-wall chain of
  four.  No cubic direction at `n = 12`.  Leaving `Z_12` at a generic point is first order (`~ -s/3` in
  the tilt `s` of one of the four axis-parallel squares) until the pinwheel is cheaper
  (`search/S6_LOCAL.md` §0–2).
* Mixed signs (`n = 6`, 116 signed directions): `k >= 3` squares in the minority sign always costs first
  order (`-eps/4` for the `(-,-,+,+,+,+)` patterns); four axis-parallel squares do **not** give `0` if the
  other two turn opposite ways (`-0.077 eps^2`); three can, for unequal tilts.  The optimal leaf's dual is
  in `111/116` directions dominated by one wall-to-wall chain of three (`search/S6_LOCAL.md` §5).
* At `theta = 0`, `n = 6`: all `2,284` leaves of the decision tree are tight, every one by a chain of
  three (`search/S6_LOCAL.md` §4) — no leaf-by-leaf continuity argument exists.
* `c_T` (uniform-tilt second-order constant at the tiling): `1, 3/4, 1/3, 3/8` exact on the leaf,
  `2/5` to seven digits at `T = 6` (permutation family only); conjecture `c_T = (T-2)/(2(T-1))` for
  `T >= 4`, increasing to `1/2` [heuristic].  The optimal hole set is always a permutation matrix; the
  binding chain's rise is `1 + O(t^2)` (`search/ARCH_TLEDGER.md` §1–3).

### 3.4 The angle-space cells at `T = 4`

[measured, `search/BANDCUT_K.md`, `search/T4_CYCLES.md`, `search/BREAK_H1.md`]

* Unrestricted, `k >= 4` near-axis squares at tilt `<= eps`: `delta* = 0` at every `eps in {1, 5, 10°}`,
  hit rates `41–695` of `240–696`.  Far field `k <= 3`: `-eps^2/(T+1) = -eps^2/5` (`gamma -> 0.1999`),
  and the same `1/(T+1)` law at `T = 3` (`0.2491`).
* **Every `delta* = 0` configuration found, at `T = 3` and `T = 4`, carries a wall-to-wall chain of `T`
  among its near-axis squares**, certifying `<= 0` through item 4 of §2.2; no chain-free configuration of
  margin `>= 0` was found (`search/BANDCUT_K.md` §0, §4).
* Chain-free (the chain among near-axis squares forbidden at the configuration's own margin):
  `-0.1..0.2 eps^2` for `k = 4..7`, and **`-0.1004 eps^3`** at `k = 8, 9` — `-5.26e-7` at `eps = 1°`,
  thinner than the `3.76e-6` of §3.1 and inside the `1.2e-3` linearisation remainder.  The configuration
  is the tiling minus a permutation hole set with the central `2 x 2` a pinwheel at tilt `eps`.  What holds
  it shut is a wall-to-wall chain of four **through two tilted squares** plus a transverse chain at
  weight `tan(eps)/5` and a rung at `tan^2(eps)/5` — a tree, not a cycle (best cycle bound `+3.39e-03`,
  four orders on the wrong side) (`search/T4_CYCLES.md` §3).
* Certificate shapes over `727` sampled optima with `delta* < 0`: `H <= 0` at `658`, best cycle `<= 0`
  at `365`, **chain-or-cycle at `690/727` — the dichotomy is false**; `H1` (H plus one further
  main-direction row) at **`727/727`**, plus `352/352` more under adversarial structured search
  (`search/BREAK_H1.md`).  Far field `k <= 3`: `124/126` duals cyclic but `104/126` have `mu >= 2`,
  median `13` pair rows on `11` squares — a thicket.  At `45°`, `T = 4`: `delta* = (7 sqrt2 - 10)/2` and
  `H(4, 5, 45°)` is exactly tight (`L* = 8 sqrt2 - 7 = 4.31 < 8`), so the `T = 3` large-tilt failure does
  not transfer; `H` fails on the uniform line at `35°` and `40°` for lack of a chain with a common normal.
* Ladder depth `2, 3, 2` at `T = 3, 4, 5`; the `eps^3` cell is a `T = 4` accident (a permutation hole set
  around a merged level cell needs `T - 2 <= 2`, a counting proof); the `T = 5` top of the hole is
  `-eps^2/6 = -eps^2/(T+1)`, closed by the bare H (`search/BREAK_H1.md` §2.1, §5).  An `eps^4` cell at
  `T = 5` is not found, not excluded.
* Best band-type configuration at `T = 4`: a wall-to-wall `(4,1)` stack at `28.5925°` beside eight
  axis-parallel squares, margin exactly `0` — a point of the zero set, still cut by an axis-parallel row
  of four; no `T <= 10` configuration with a band cutting every row and column at margin `>= 0`
  (`search/BANDCUT_SCAN.md` §5).

### 3.5 Duals that are not covers, and caps that do not bind

* The best sub-12 `t = 4` packing-side duals (`11.80` certified support, `11.86` loaded set) read as
  covers over the continuum cost `20.162` and `18.109`; `16.09` after one repair round; nothing below
  `15.8` under any variant.  The pricer's clique rule is unsound (`304/317` rows admit two disjoint
  credited squares); every dual-carrying clique row is non-Helly with an empty core (deficits
  `-0.0007..-0.10`); the `0.04 / 2.5°` lattice hides the worst pose by `6x` — the tiling corner square
  nudged `1e-7` captures `0.655` (`search/HONEST.md` §0).  Sound positive-volume boxes exist for every row
  (`317/317`, `638/638`, median volume `2 %` of pose space) and the honest cost is then `145` to `1543`
  (`search/ALLMEET.md` §2).
* A cap `mu(tilt < eps) <= 3` costs the points-only relaxation `0.002512` at `eps <= 1°` (the LP tilts
  the corner square to `1.1997°`), `0.0156` at `2°`; at `T = 3`, where `s(6) = 3` is a theorem, the cap
  `<= 2` leaves `V >= 6` for every `eps <= 5°` and bites only past `10°` (`search/ARCH_FARFIELD.md`
  §3–4).  In the strongest family the same cap is free at `1°` and `2°` (multiplier `0`) and costs
  `0.11–0.58` on the pools tried at `3°–30°`, numbers that recover by half when one family of columns is
  added — pool-limited, silent (`search/FARFIELD_STRONG.md` §4).  Cap `9` never binds anywhere.
* The verified `T = 11` packing has no wall-to-wall path of axis-parallel squares; its five tight chains
  of eleven run through the band with `2–4` tilted squares each; on each, the rise-free budget is
  `-0.62..-0.96` and the rise credit slightly larger; a band link nets `+0.10..0.13`, an interface
  `-0.16..-0.18`, and three band links break even — the whole `+2.7e-4` (`search/T11_CHAINS.md`).  The
  ladder's clause that fails there is leg existence (one link per band kink, `L` available `8–15`
  against `13–16` needed) and the wall-block placement (`notes/counting-ladder.md` §3).

---

## 4. The routes tried, and why each closed

In the order they were tried.  "Closed" means: measured or proved unable to certify `s(12) = 4` by
itself, with the reason named.

### 4.1 The pure weighted-cover LP (the method that gave `3.968616`)

The method proves `s(12) >= s` iff `nu_f(s) < 12`.  Exact fractional packings of mass `12.0282` at
`s = 3.99` (`search/DUAL_EXACT.md`, `search/CLIQUE_CONTINUUM.md` §3) and `12.2688` at `s = 4`
(`search/COVER4.md`) pin the ceiling of the pure method to `s* in [3.968616, 3.99)` [proved]; the
extremal measure is the tiling smeared: `70 %` of its mass under `5°` of tilt, four corner squares at
`0.81` each (`search/COVER4.md`).  Nothing below `t = 4` transfers to a proof of `s(12) = 4`, and at
`t = 4` the LP is `0.27` above `12`.  **Closed by §2.2 item 2.**

### 4.2 Clique rows (rotated squares are not a Helly family)

The corner leaf's dual violates clique constraints `sum_{S in K} mu_S <= 1` by `0.5` at `t = 3.98`
(`search/CLIQUE.md`), and clique columns pushed the `t = 4` packing-side numbers to `11.75–11.80`
(`notes/review-2026-09-12.md`).  Read as covers over the continuum those duals cost `18.1–20.2`; every
dual-carrying clique row is non-Helly with an empty core, the pricer's crediting rule is unsound, and no
sound rule can push the honest cost below the `meet` floor (`search/HONEST.md` §0, `search/ALLMEET.md`
§0 item 1) [measured].  Sound positive-volume boxes exist for every row and cost `145–1543`
(`search/ALLMEET.md` §2) [measured].  The Helly-type statement that would have made cliques worth points is
false at `t = 4` (§2.2 item 12) [proved].  The verifier `V1` for cliques was never built because no
certificate-shaped object below `12` ever existed.  **Closed: no sub-12 number carried by a clique row is a
bound on anything certifiable.**

### 4.3 Points + polygons + regions + chord, and the count trees

Without cliques the certifiable family sits at `12.17` pure, `12.03` at corner `k = 3`, exactly `12.000`
at corner `k = 4` (§2.1) [proved lower bounds]; polygons are worth `0` at the corner leaf and the chord
row's dual is `0` everywhere (`search/T4LEAF.md` §3.1) [measured].  The level-2 slot branch, once made
verifier-compatible, is worth `0.08` (no cliques) to `0.25` (with anchor cliques) on its hardest leaf,
not the `0.4–1.7` the screen measured on starved pose sets; that leaf reached `11.920` certified and was
still rising (`search/T4LEAF.md` §0, §3.3) [measured].  The full corner x slot tree has `4,213` leaves;
the tree under `k < 4` was never run.  **Closed as a strategy** (`notes/review-2026-09-13.md` decisions:
no more count or pattern trees at `t = 4`): every leaf measured to convergence is `>= 12`, the one not
converged is rising, and the extremal object is always the tiling smeared.

### 4.4 The Bentz template: incidence patterns on his 16 points

Bentz's `s(13)` proof is a 16-point set plus a six-leaf case tree; at 12 boxes his counting reads
`K + u = 4` (identity, `search/BENTZ.md` §2) [proved] and the four corner boxes can each hold a doubled
pair — exactly the tiling minus four (`notes/proof-anatomy.md` §7.2).  Measured at `t = 4`: a closed unit
square holds `93` of the `65,536` subsets of the 16 points, never 5 and never 0 (exact witnesses; "no
others" is a lattice scan) [measured]; `86,403` pattern leaves, `10,945` up to `D4`; the corner-pattern
level has `43` classes and its `(AB)^4` class is already fully pinned.  That leaf deletes `12 %` of the
poses and **none of the 162 that carry the optimum**: bracket `[11.999999926, 12]` (§2.1) [proved].  The
16 points are a closed cover of `[0,4]^2` on every test but at margin `0` (the tiling) and `2.2e-4`
(the `44°` band), so `zeromargin.py` cannot certify it (`3,516` boxes uncertified) [measured].  A pattern
region is a union of pose clouds by construction, and the leaf's mass is carried by `10–29` poses per
pattern differing by `1e-3` in centre and tens of degrees in angle (`search/BENTZ.md` §5 (ii)).
**Closed: incidence branching on any finite point set has no power at `t = 4`.**

### 4.5 The rank-8 statement and leaf A

The obstruction localised to eight labelled squares (§3.1).  A proof of the rank-8 statement cannot be a
local perturbation argument (the linearisation leaves `O(rho alpha + alpha^2)` of positive sign, three
orders above the margins that matter), cannot be split into smaller sub-families (every eleven-square
sub-configuration is an exact packing), cannot find a positive margin anywhere (it is `0`, attained, on a
plateau reaching `32.5°`), and is not seen by level-2 Sherali–Adams beyond a tenth (`search/RANK8.md`
§7).  Bentz's replacement lemma applies to the eight singleton boxes but is degree-1 and cannot kill by
double coverage (`notes/review-2026-09-13b.md`).  Degree 2 in the continuum has no crediting rule
(`search/ALLMEET.md` §0 untouched by it).  **Closed as a target for existing tools; it survives as the
statement any proof must contain.**

### 4.6 The census of leaves, and the separation structure

A per-leaf hand proof is not a strategy: `>= 2,779` zero-margin classes, a quarter of them without any
tiling in them (§3.2).  What the census pointed at — branch on separation type per pair rather than
point incidence, kill row-of-four leaves by an exact linear argument — became the angle-space route.
**Closed as a point-pattern tree; superseded by 4.7–4.10.**

### 4.7 The angle-space box tree (rehearsed on `s(6) = 3`)

With angles fixed the problem is a disjunctive LP in the centres (`search/S6_SKELETON.md` §3.1) [proved],
decidable at any single angle vector (validated to 9 decimals; `3,213` LP nodes at the axis-parallel
point, `53,445` at `45°`, not finished at a generic vector in `150,000`) [measured].  Over an angle box the
relaxation errs by `+W^2/2` in the box width, exactly (`0.4984 w^2` at `0.125°`), and any box meeting the
zero set has `sup = 0`, so no positive-width box meeting `Z` can be killed; the deficit off `Z` along the
worst direction is `-(1/4) t^3` at `n = 6`, so the refinement factor grows like `rho^{-1/2}` and the box
count diverges [measured + inference].  `0` boxes killed at `22.5°` and `11.25°` bins; at `n = 12` an
angle-box tree would need bins under `0.16°` in 12 dimensions to see the `3.76e-6` sub-family.  The exact
kill lemmas K1 (transversal chord, band `D(t) >= (sqrt2-1)/2`) and K2 (fractional chord, `k lambda >= T`)
are [proved] and never fire in the branch-and-bound (`search/S6_SKELETON.md` §2, §4.5).  **Closed: any
method reading the angles through interval bounds loses at order 2 against a deficit of order 3
(`n = 6`) or order 2 with no radius (`n = 12`).**

### 4.8 The architecture: tube around the zero set + far field

`notes/proof-architecture.md` composed the pieces (`S0 + S1 + A3 far field + A4 chain existence + A2
local theorem + A5 glue`) and its adversarial review found the hole: chain counting forces a chain of
`T` only for `k > (T-1)^2 = 9` near-axis squares (Lemma CC, sharp), the far-field cap must handle
`k <= 3`, and **`4 <= k <= 9` is covered by nothing** — and is inhabited (the `j = 4..8` strata of §3.3,
the `(4,1)` band stack of §3.4).  The local theorem at the tiling is `T`-uniform (`c_T > 0` at every `T`,
`B_T(1,t) <= 0` for all `T >= 3`) and is true at `T = 11` where the conclusion is false, so it carries
none of the theorem (`search/ARCH_TLEDGER.md` §0).  The far-field cap on the degree-1 LP is refuted at
small `eps` in the points-only family (cost `0.0025`; the LP tilts the corner square to `1.2°`) and in the
strongest family (cost `0`, multiplier `0` at `1°`, `2°`), and at `T = 3` — where the conclusion is a
theorem — it does not bite until past `10°`, far outside any tube the local lemmas can supply
(`search/ARCH_FARFIELD.md`, `search/FARFIELD_STRONG.md`) [proved refutations, measured elsewhere].  No
cover-side instrument exists for the strong family (its duals are not covers, `search/HONEST.md`), so no
region of the architecture is proved by machine either.  **Closed: the spine is sound and `T`-uniform;
the two lemmas that would carry `T` (`A3`, `A4` in the hole) have no candidate.**

### 4.9 Band cuts and the waste budget

The `T >= 12` counterexamples are two integer blocks plus two squeezable rectangles with
`waste(A) + waste(B) = T`; the fit `waste = min(side) + 2`, `a >= 4` gives `T >= 12` by integrality
(`notes/bandcut-cost.md` §3) [measured fit, proved arithmetic] — and over-proves by one: `T = 11` is a
verified positive whose band is bent, not two rectangles (§2.2 item 13).  The length/chain accounting has
no threshold in `T` at all (gain and cost both linear in `T`, the gain's coefficient a free parameter)
[proved]; a structure-free waste bound is refuted by `s(5)` (`waste 2.33`) [proved]; "no `(3,b)` rectangle
is squeezable" is false (`s(5)` plus a column of two) [proved].  At `T = 4` no band cuts every row and
column at margin `>= 0` (§3.4) [measured].  **Closed: the waste route is a theorem about one paper's
scheme, not about packings.**

### 4.10 Chains, cycles and the ladder (the `T = 3` rehearsal)

The candidate "every configuration has a tight wall-to-wall chain of `T` on which MT gives `<= 0`" is
false at `T = 11` (rise control fails, §3.5) and MT is a bound, not a Farkas certificate (it reads the
rise off the configuration; `notes/t3-chain.md` §2.7).  Making it a certificate is the H-lemma
(§2.2 item 7), which at `T = 3` fails above `28.8°` (only `L = 2` transverse links available against
`L* = 2.485`) and sporadically in the coherent small-tilt cone; the repair there is one extra
main-direction row, which at `T = 3` closes a cycle and at `T = 4` usually grows a branch
(`search/T4_CYCLES.md` §3).  Chain-or-cycle is false at `T = 4` (`37/727`); `H1` holds at every sampled
optimum (`727 + 352`), has a closed form, and its `T`-dependence is entirely in two integers of an
existence clause (`L >= T-2` transverse links; tilted block one square in from the wall)
(`notes/counting-ladder.md` §2.6, §4) [proved forms, measured existence].  Kearney–Shiu's `s(6) = 3` is
not a chain proof: unavoidable-set counting with the tilt quantified away once (`notes/t3-chain.md` §4).
**Reduced to existence clauses; all `T`-dependence in them; none proved.**

### 4.11 The closure: obstruction X and step S  (2026-09-22, Evan's rule)

Rule applied: no route for `s(12)` that cannot handle `s(6)`.

* **Scope.**  The existence clause quantifies over packings (`delta >= 0`); every recorded H-failure
  (`68/720` at `T = 3`) and dichotomy failure (`37/727` at `T = 4`) is at `delta < 0`, outside it.
  Inside scope, at all `749` margin-`0` packings found at `T = 3`, a **bare** wall-to-wall chain of three
  with exactly axial normals certifies `delta <= 0` (Lemma Z''); every such packing has `>= 3` exactly
  axis-parallel squares, never fewer (`notes/t3-existence.md` §1, §6) [measured].  The `45°` pentagon,
  the cycle and the far-field thicket all concern configurations that are not packings; at `T = 3` there
  are none above `25.8431°` at all (§2.2 item 10) [proved].
* **Obstruction X.**  [proved that each named mechanism fails; the obstruction itself is a diagnosis]
  The surviving clause — three exactly axis-parallel, consecutively separated squares — must separate
  packings from configurations of margin `-0.36 t^2` (uniform tilt) and `-t^3/4` (three exact zeros) as
  `t -> 0`.  Chain counting is exactly tight at `t = 0` and reverses at first order (`x`-step
  `1 - 2t` against box `2 - t`); the centre pigeonhole is `20 %` on the wrong side at `t = 0`; unavoidable
  points are angle-free by construction; Menger/Dilworth is chain counting.  Every mechanism stable under
  an `O(t)` perturbation of the angles is blind at `O(t^2)`.  So X is a second-order rigidity statement
  about the tiling-adjacent stratum — the "(ii) second-order theorem on the coherent cone" of
  `search/S6_LOCAL.md` §5 — and the angle-space route does not reduce it: it is it.  `T = 4` inherits X
  one order flatter (`-0.1 eps^3`) and loses the far-field pigeonhole (`notes/t3-existence.md` §4.4, §8).
* **Step S (the chord-row salvage).**  Chord rows give the first certificate that is tight at `t = 0` and
  second-order right (Prop. CC3, `-0.375 t^2` against the truth `-0.3637 t^2`, `3 %` apart) [proved].  But
  the leaf it needs — three centres within a window of length `cos t - sin t` — is produced by the sharp
  pigeonhole only if `(T-u)/2 <= C - S`, i.e. only at `t = 0`, with deficit `(T-2)t/(T-1)` at every `T`
  (Prop. PG) [proved]; the existence clause `(E3-chord)` is **false** at the pinwheel (`748/749`, the
  staircase is the gap) [proved]; and the truth at the uniform-tilt optimum misses the leaf by
  `0.63 t^2`, second order, so proving non-emptiness needs `O(t^2)` accuracy in the centres — X restated.
  Kearney–Shiu contains exactly one Farkas certificate (five rows, weights `(1,1,1,1,1)`,
  `delta <= (13/6 + 2 - 3 sqrt2)/2 = -0.037987`); everything else is unavoidable-set counting with no
  row-language form, and that counting is the existence step (`notes/t3-chord.md` §2, §4).
* **Verdict.**  `s(6) = 3` is not proved in the row language with or without chord rows; by the rule the
  route is closed for `s(12) = 4` (`notes/review-2026-09-22.md`).  Kept: Lemma CC, Prop. CC3, Corollary P3
  at `T = 3`, the ladder closed form, the master formula.

---

## 5. The sanity filter every proposed mechanism must pass

`s(n^2 - n) < n` for **all `n >= 12`** (Arslanov–Mustafin–Shangitbayev, Electron. J. Combin. 28(4) 2021,
P4.22, refereed) and for `n = 11` (Cantrell, Feb 2025, unrefereed, `s(110) <= 10.9967932740…`); both
records verified in this repo in exact rational arithmetic at margin `+6.81e-4` and `+2.75e-4`
(`search/BANDCUT_SCAN.md` §2) [proved].  So the conjecture `s(n^2 - n) = n` can hold at most for
`n <= 10`, `n = 4` being the target and `n = 5..10` open.  Separately `s(n^2 - 4) < n` at `n = 3`
(`s(5) = 2 + 1/sqrt2`).  Consequences, each checked against the lemma it filters:

* An argument uniform in `n` is wrong.  Every exact inequality in the repo — K1/K2, the `(T-1)^2`
  count, the total-band inequality `sum D(t_i) < T(T-1)/2` (exactly critical at `n = T^2 - T`), Lemma CC,
  Lemma H and its `L*`, Lemma W', the master formula, the cycle lemma, the ladder, MT, `B_T` — is
  `T`-uniform and true at `T = 11, 12, 17` (`notes/t3-chain.md` §6, `notes/counting-ladder.md` §1.6, §4).
  The `T`-dependence of the whole problem is quarantined in existence clauses (a chain with `L >= T-2`
  legs placed one square in; a chain through the near-axis squares in the hole `4 <= k <= 9`), and in
  those alone.
* A proof must use the chain length `4` quantitatively.  What `T = 11` has that `T = 4` cannot afford: a
  chain of eleven with room for three full-credit band links, seven axis-parallel squares and two
  interfaces; at `T = 4` a chain with three band links is the band itself (`search/T11_CHAINS.md` §0(4))
  [heuristic].  Cleemann's band cuts every row and column of `[0,17]^2` with `199` of `272` squares
  axis-parallel, so `k = 199 <= 256 = (T-1)^2` and every counting lemma is silent there, exactly as
  Lemma CC's sharpness predicts (`search/ARCH_TLEDGER.md` §5).
* A proof must not prove `s(5) = 3`: at `n = T^2 - 4`, `T = 3`, the optimal packing has
  `k = 4 = (T-1)^2` near-axis squares and one at `45°`; every chain of three through the tilted square
  pays `(sqrt2 - 1)/2` per link and Lemma H returns `+0.146` (`notes/t3-chain.md` §6.3); Lemma CC's
  (H6) fails by exactly one; Corollary P3 is silent (`S d_5 >= 1.12`).  All passed, by the tilt, never by
  counting.
* A proof must survive margin `0` as a hypothesis, not a convenience: every restricted-dual bound is
  `>= -4e-12` at every margin-`0` point tested (`30` points, `search/BREAK_H1.md` §6), and the one
  formula that went negative there (`-0.1399`, the cycle lemma without its axis-aligned-turn hypothesis)
  was thereby caught as false (`search/T4_CYCLES.md` §6).

---

## 6. What is open

1. **A 13-point pure unavoidable set for `[0,4]^2`** (`tasks/unavoid13`, running now; cite as **open**).
   `>= 13` is a theorem (§2.2 item 2); `<= 14` is Friedman's set.  A 13-point set is what the
   Kearney–Shiu slack-`1` accounting needs at `n = 12` (`12` squares against `13` points, the same
   `7 - 6 = 1` as at `n = 6`); the `90°` duality, Lemma 1 and Lemma 2 + (7) transfer, only Lemma 3's
   `5/3` is container-specific (`notes/t3-chord.md` §6.2).  No LP bound can decide it: the certified
   fractional floor is `12.2688`, the fractional optimum is believed `~12.4`, and a weighted cover of
   weight `12.956` exists, so both `13` and `14` are consistent with everything fractional.  It is an
   integrality question, decidable in principle by finite hitting-set computation in either direction.
   Even with the set in hand, the rest of a K–S-style argument at `T = 4` is unwritten.
2. **The existence clauses** of the angle-space route: at `T = 3`, (E3-row) on the strata "exactly `3`
   or `4` axis-parallel + tilted rest" and non-existence on the coherent cone (obstruction X); at `T = 4`,
   the `H1` conjecture (`727/727`, `352/352`) and "every `delta* = 0` configuration carries a chain of
   four among its near-axis squares" (`BANDCUT_K.md` §0(1)).  All measured, none proved, and X says why no
   `O(t)`-stable mechanism can prove them.
3. **Whether Lemma CC's order chain can be taken as a DAG path** at `eps > 0`
   (`notes/counting-ladder.md` §1.4); **DAG acyclicity at large tilt** (open beyond
   `eps < (1+delta)/(T-1-2 delta)`).
4. **The far field at `T = 4`** (`k <= 3` near-axis squares): the centre pigeonhole never fires
   (`0.5 %` short), the chord pigeonhole never fires, the capped LP does not close, and the true dual there
   is a thicket with `mu >= 2`.  No candidate mechanism.
5. **`V(4, eps, 3)` for `eps >= 5°`** in either family (pool-limited both ways); `eps^4` cells at `T = 5`
   (not found, not excluded); `c_T` beyond `T = 6`; the `T = 5..10` cases generally.
6. **Whether the 16 Bentz points are a closed cover of `[0,4]^2`** (evidence only; margin `0` at the
   tiling and `2.2e-4` in the `44°` band; not load-bearing for anything above).
7. Explicitly *not* worth pursuing, by decision (`notes/status.md`): interior/quadrant branching on
   symmetric leaves; `3.975–3.978` from the corner level; SDP / Lovász-theta in the continuum; pure
   `3.99` to convergence; line-chord cuts at `t = 4`; any further count or incidence tree at `t = 4`;
   testing new inequalities axis-parallel first; level-2 branching by wall total; the `s(12) >= 3.98`
   consolation prize.

The shape of what a proof must contain, as the record leaves it: a statement about a few squares at
once (rank 8 at least; the eight non-corner tiles), true at margin exactly `0`, accurate to second order
in the tilts near the tiling and to third order in the `k = 8, 9` cell, using the chain length `4` as a
number, silent at `T = 11` and at `n = 5`, and not reducible to any single-square inequality, any
point-incidence pattern, any interval relaxation of the angles, or any `O(t)`-stable counting.

---

## 7. Superseded claims and errors caught (pointers, so a reader knows which notes are stale)

Each line: the stale statement → where it is corrected.

* **"Corner leaf exactly `12` with a dyadic grid cover" (early September)** → withdrawn; that run never
  produced a valid cover (`search/CLIQUE.md`, README "Beyond the ceiling").
* **The heuristic closed covers of `search/CLOSED4.md` (best `12.51`)** → invalid at poses their row
  lattice stepped over (captures `0.942`, `0.970`, `search/RUNG2.md` §10); the exact floor is `12.2688`.
* **`T4SCREEN.md` §3's "`k = 4` pure reads exactly `12.000000`"** → an LP value on an unconverged row set
  with `M = 1.5186` and a violated clique; not a value (`search/T4LEAF.md` §4, `search/LEAF_CEILING.md`
  §5.2).  The slot branch's `0.4–1.7` → `0.08–0.25` (`search/T4LEAF.md` §3.3).
* **Sub-12 `t = 4` numbers `11.75–11.80` (2026-09-07/08) as evidence for a certificate** → carried by
  non-Helly cliques; honest cost `18–20` (`search/HONEST.md`); "E1 for the certifiable family: no"
  (`notes/review-2026-09-12.md`).
* **The Helly-type conjecture of the 2026-09-12 TODO** → false, exact triple (`search/ALLMEET.md` §4).
* **"Chord row is a flagged hypothesis" (`search/T4LEAF.md` §1.2, repeated in
  `search/FARFIELD_STRONG.md` §7 item 2, `notes/review-2026-09-22.md`, `notes/status.md`)** → the
  wall-strip lemma is proved in Lean for closed-disjoint unit squares at `t <= 4`
  (`notes/chord-lemma.md`, `lean/Sqpack/Chord.lean`, 2026-09-07); the row is sound in the repo's
  semantics.  The residual caveat is only the `1e-9` coordinate tolerance argued in `T4LEAF.md` §1.2.
* **`S6_SKELETON.md`: `n = 12` uniform-tilt coefficient `-(2/3) t^2`; `n = 6` cubic `-0.17 t^3`;
  "cubic constant at `n = 12` is the number to measure next"** → `-(1/3) t^2` (under-optimised);
  `-(1/4) t^3` exact on its leaf; no cubic direction at `n = 12` (`search/S6_LOCAL.md` §0).
* **`proof-architecture.md` §3.2 "`c_T` extrapolates through zero before `T = 5`"** → two branches read
  as one; `c_5 = 3/8` exact, `c_T` never changes sign (`search/ARCH_TLEDGER.md` §6).  §0a item 1
  "`k >= T` near-axis squares force a chain" → needs `k > (T-1)^2`; and "`eps < 1/(T-1)` forces the
  chain" → it only makes the chain span the container, and the lemma gives an order chain, not a DAG
  path (`notes/counting-ladder.md` §5).  §0a item 10's erratum: "cap `>= 4` must give `>= 12` because a
  margin-`0` packing exists" is invalid — touching squares have coverage `2` (`search/ARCH_FARFIELD.md`
  §5(a)).  "First failing `T*` in `5..17`" → `5..11` (`search/ARCH_TLEDGER.md` §5.1).
* **`bandcut-cost.md` / brief: "no `(3,b)` rectangle is squeezable"; "`waste >= a + 2` in every small
  case"** → false (`s(5)` + a column of two; `(3,3)/5` has `a + 1`) (`notes/proof-architecture.md` §0a
  item 13).
* **`BANDCUT_K.md` §4's `B_T` column and `bandcut_k.py chain_certificate`** → evaluates `B_T` at the
  mean tilt without a common-normal check and converts degrees twice; invalid at large tilt, approximate
  to `O(eps)` on near-axis chains (`notes/t3-chain.md` §0 erratum).  §5's "chain-free rows carry weight
  `0`" → some are exactly tight; the true statement is about `t3_chain.all_rows` (`search/BREAK_H1.md`
  §7 item 5).
* **`t3-chain.md`**: §5.3 cycle lemma's turn constant `2 sin(alpha/2)` → the `L1` jump `2 cos t`; the
  "every turn axis-aligned" hypothesis is load-bearing and nearly empty at `T = 4` (`search/T4_CYCLES.md`
  §4.1, §6; `notes/t3-existence.md` §2.3).  §0(4) "one extra edge is one cycle" → at `T = 4` it is usually
  a branch (`search/T4_CYCLES.md` §3.1).  §2.3 "H at `45°` is a five-square certificate" → it is a
  three-square one; the legs are worth nothing at `45°` (`notes/t3-existence.md` §7).  §0(1) "the whole
  lemma is one integer" → the ladder needs two (`notes/counting-ladder.md` §5).  §4.4 "`s(12) >= 4`
  needs an `11`-point unavoidable set" → that is the slack-`0` count, now excluded; the slack-`1` count is
  `13`, open (`notes/t3-chord.md` §6.2).  §3.4 "cycle of eight at `45°`, `T = 4`" → exists but a
  `4`-cycle with `q = 2` suffices (`search/T4_CYCLES.md` §4.3).
* **`t3-existence.md` §5.3's chord row** → false as written twice: the chord is not centred at the
  centre, and "meeting the line" must be interior-meeting under closed semantics; chord rows add centre
  branching, not a new valid inequality (`notes/t3-chord.md` §1.2, §7; `search/S6_SKELETON.md` §4.5 had
  it right).
* **`T11_CHAINS.md` §0(3) interface signs** → the review of 2026-09-22 says "inverted";
  `notes/counting-ladder.md` §5 item 5 says nothing there is wrong, only that the convention (interface
  terms are contributions to `L`, strengthening the bound) is easy to misread.  §3 "the content is
  entirely in bounding the rise" → too narrow; the dual side is transverse rows with a common normal
  (`notes/counting-ladder.md` §5 item 4).
* **`notes/status.md` (2026-09-22)**, small: quotes `11.99999987` for both `eps = 1°, 2°` cells of
  `FARFIELD_STRONG` (sources `11.999999870`, `11.999999855`); still carries the "five-square" reading
  above; labels the corner-`k = 4` family "`<= 12`" in the same breath as "certified" (the `<= 12` is a
  measured LP value for the family and a theorem only for the pinned leaf); "at least `2.6` of the `3`" is
  `2.599`.  The `sub-degree` figure in `search/BENTZ.md` §7.2 is also written "`>= 2.60`".
* **Brief-level premises found wrong along the way** (each in the note it belongs to): "the rise is
  bounded by transverse chains, so the lemma is H-shaped" (`notes/review-2026-09-21.md`); "test the
  candidate lemmas at non-optimal configurations with `delta >= 0`" (vacuous: every such configuration
  is an optimum, `notes/t3-existence.md` §1.1); "the chain-free `k = 8, 9` cell is held by a cycle" (a
  tree, `search/T4_CYCLES.md` §8); "ladder depth grows with `T`" (`2, 3, 2`, `search/BREAK_H1.md` §7);
  the `T = 3` analogue of the strongest family is degenerate at exactly `5` (`search/FARFIELD_STRONG.md`
  §7 item 3).

Older `notes/review-*.md` (09-07, 09-08, 09-11) predate every number in §2 and should be read as history
only; `notes/proof-architecture.md` is to be read with its §0a first, and `search/S6_SKELETON.md` with the
correction box at its top.
