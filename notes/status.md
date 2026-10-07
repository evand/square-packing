# Where the proof stands (2026-09-22)

> **Frozen at 2026-09-22.**  Current status of every result: the results table in [`README.md`](../README.md).

Moved out of `TODO.md` on 2026-09-13; update when a row changes.  Status log: `notes/review-*.md`.

## Since 2026-09-22 (added 2026-10-03)

The rest of this file is frozen at 2026-09-22 and is about `s(12)`, which is still open.  Exact results since, all
computer-assisted certificates with a `verify.sh`, none peer reviewed (per-result trust: `README.md` § Results):

* `s(32) = 6` (2026-09-26), `certificates/s32/`; kernel-checked in Lean with no hypothesis (`s32_eq_6`).
* `s(21) = 5`, `s(45) = 7` (2026-09-27), `certificates/{s21,s45}/`; two exact checkers each.  Lean: `s(21)`
  conditional on the checkers' region statement; `s(45)` not in Lean.
* `s(60) = 8`, hence `s(61) = 8` (2026-09-28), `certificates/s60/`; two exact checkers; not in Lean.
* `s(k² − 3) = k` for all `k ≥ 6` (2026-09-29), `certificates/k2m3/`; checker independently corroborated (wand125).
* `s(k² − 4) = k` for all `k ≥ 5` (2026-10-03), `certificates/k2m4/`: `k ≥ 8` by one family and one finite statement
  certified by one implementation; `k = 5…8` are the bundles above.  `s(77) = 9` was first proved by wand125
  (2026-10-01).

Lean state (2026-10-03): `bentz_of_validTilt7 : ValidTilt7 → ∀ k ≥ 6, minSide (k² − 3) = k` and
`bentz4_of_validTilt9 : ValidTilt9 → ∀ k ≥ 8, minSide (k² − 4) = k`, with the all-k reduction, the D4 reduction and
Lemma Z kernel-checked; the tilted runs `ValidTilt7`, `ValidTilt9` are the remaining hypotheses
(`notes/lean-valid-split.md`).  (Also: `s(11)` was settled by others on 2026-09-29, `README.md`.)

Certified: `s(12) >= 3.968616`; `s(11) >= 3.8143`, and for `n = 11` the pure method's ceiling is
bracketed `3.814304 <= s*(11) <= 3.83375` (`search/N11_ANATOMY.md`; `COVER(3.8191) <= 11.00004` verified).  **`n = 11` is now on
the same footing as `n = 12`** (`search/N11_LADDER.md`, 2026-09-14): every rung run once — clique-strengthened
relaxation exactly `32/3` on `[3.83375, 3.87]` (cliques take `1/3` of the missing unit, from one interior pose);
corner branch vacuous; box and anchor cliques worth exactly `0` on the cover side; SA-2 on the 28-square support
closes the whole unit (`= alpha = 10`, exact) where at `n = 12` it closes a tenth.  Pure gap to the conjecture:
`>= 0.043` at `n = 11` vs `<= 0.031` at `n = 12`, rigorous both ways.  **`s(13) = 4` case-free** (`certificates/rung2/`,
`search/RUNG2.md`): two exhaustive checkers sharing nothing, 23 rejection tests, Lean for every
soundness lemma, CI green, written up (`notes/s13-casefree.md`, README, `docs/`).

For `n = 12` at the container (`t = 4`, closed semantics), every degree-1 certifiable family is `>= 12`
(`search/BENTZ.md`, `notes/review-2026-09-13.md`):

| what | value | status |
|---|---|---|
| pure, points only | `>= 12.2688` exact | `COVER4.md` |
| pure, points + polygons | `12.173398` certified | `pgonly_pure_exact.txt` |
| corner `k = 3`, points + polygons + regions + chord | `12.027273` certified (LP `12.038` rising) | `BENTZ.md` §7.2 |
| corner `k = 4`, same family | `11.999999926` certified; LP pinned at `12.000000000` on every pose set (measured — `<= 12` is a theorem only for the pinned leaf, next row) | `pgonly_corner_exact.txt`, `LEAF_CEILING.md` §5.2 |
| corner `k = 4`, fully pinned Bentz leaf `(AB)^4` | `11.999999926 <= value <= 12` exactly | `BENTZ.md` §0 |
| corner `k = 4`, axis-parallel only | exactly `9` (theorem) | `BENTZ.md` §0, C1 |
| corner `k = 4`, `|θ mod 90°| <= 1°` | `11.599390` certified | `BENTZ.md` §7.2 |
| sound cliques of any positive-volume rule | honest cost `>= 18.109` (`>= 20.162` in the superset family) | `ALLMEET.md` |
| level-2 slot leaves under `k = 4`, no cliques | `11.92`, rising, not yet `>= 12` | `T4LEAF.md` |

The obstruction is rotation, and infinitesimal rotation: axis-parallel is `<= 9`, at least `2.599` of the
missing `3` is bought by tilts under one degree, and the LP averages over that continuum.  The
integrality gap is exactly one square and visible in integers: on the optimum's own 162-pose support
`α = 11`; the four corner pairs, four `{c_j}`, four `{d_j}` are each realisable alone, and **the eight
non-corner singletons together admit only seven**.  The theorem a proof needs is that rank-8 statement
— with four corner squares holding `{a_i, b_i}`, at most seven further squares each contain one of
`c_j, d_j` — true integrally (one wall square per wall, then `s(4) = 2`), inexpressible by any
single-square inequality, and (since tiling-minus-4 touches) almost certainly of margin zero.

Measured (`search/RANK8.md`, 2026-09-13): its margin **is** zero — the sup of the min pairwise closed
gap over the leaf is `0`, attained, and attained on a plateau reaching `32°` of tilt and `0.9` of
centre away from the tiling, so the maximisers do not converge to the tiling.  It has **no smaller
core**: every sub-configuration of eleven or fewer of the twelve squares is realisable, each with an
exactly verified rational witness, the tightest being "drop a corner square" at `3.76e-6`.  The
linearised system around the tiling is first-order rigid for all sixteen tiling-minus-4 families, with
an irreducible certificate made of four-square chains pinned between opposite walls; its second-order
remainder is `O(rho alpha + alpha^2)`, so the linearisation locates the obstruction but cannot close
it.  Sherali-Adams level 2 on the 162-pose support gives `11.893347` (and `7.894640` on the eight
singletons alone), i.e. about a tenth of the missing unit, on a dual carried by the interior-interior
and wall-interior tile contacts only.

Census (`search/CENSUS.md`, 2026-09-19; multistart measurements, not certificates): over all 10,945 `D4`
classes of pattern leaves no class has positive margin, and leaf A is one of `>= 2,779` zero-margin classes
(1,425 contain a tiling-minus-4; the rest sit at `0` on configurations tilted up to `44°`).  What holds a
plateau point shut is one axis-parallel wall-to-wall row of four (99 %) or a corner-contact pinwheel core
(axis-parallel configurations only); no core ever needed a tilted square.  So the point-pattern tree is
blind to the obstruction, a per-leaf hand proof is not a strategy, and the next step is a skeleton with
exact kills on open conditions, rehearsed on `s(6) = 3` (`tasks/s6-skeleton`).

`s(6) = 3` rehearsal (`search/S6_SKELETON.md`, 2026-09-20; measurements): same picture with 3 for 4 (98.5 % rows
of three, pinwheels, no tilted core).  With angles fixed everything is a disjunctive LP in the centres, so the
theorem lives in angle space (6 dims; 12 at `n = 12`), but a box tree there cannot close: the relaxation error is
`W^2/2` and any box meeting the zero set has `sup = 0`.  Rows die at first order; the pinwheel is the hard
direction and has no exact kill yet.

Deficit off the zero set, re-measured with structured starts and made exact on single leaves
(`search/S6_LOCAL.md`, 2026-09-20; corrects `S6_SKELETON.md`): `n = 6`: `-(3/4) t^2` uniform tilt (pinwheel; on its
leaf exactly `-3(u-1)^2/(u^2+3)`, `u = cos t + sin t`), `-(1/2) t^2` with one or two squares axis-parallel,
`-(1/4) t^3` with three (exact on its leaf; the earlier `-0.17` was below LP tolerance), `0` with four or more —
**for a common tilt of one sign only**: over 116 signed directions of the small-angle cube (`S6_LOCAL.md` §5) mixed signs
cost first order (`-eps/4` for two or three squares turned the other way), four axis-parallel squares do not give `0`
if the other two turn opposite ways, three can if the other tilts are unequal; the optimal leaf's dual is always
dominated by one wall-to-wall chain of three.
`n = 12`: **`-(1/3) t^2`, not `-(2/3)`** (the old multistart was under-optimised), the same to 10 digits with 0–3
squares axis-parallel, exactly `0` with `>= 4` (eight squares tilted to `12.8°` beside a chain of four): no cubic
direction, and leaving the zero set from a generic point is first order (`~ -s/3`) until the pinwheel is cheaper.
The `n = 12` pinwheel dual is a wall-to-wall staircase chain of four at weight `1/3` plus transverse links at
`t/3`, `t^2/3`.  At `theta = 0`, `n = 6`, all 2,284 leaves of the decision tree are tight, every one by a chain of
three: no leaf-by-leaf continuity argument; the local theorem is one statement about chains.

Conditional results actually proven: at `t = 4` only "all axis-parallel ⇒ `<= 9`"; at `t <= 3.98`, every 12-packing has `>= 3` corner boxes occupied.  Trivially `<= 9` orthogonal squares in any packing.

**Angle-space route closed at `T = 3` (2026-09-22, `notes/review-2026-09-22.md`).**  Tested by Evan's rule (no route for `s(12)` that
cannot handle `s(6)`): `notes/t3-existence.md` finished the certificate side hypothesis-free, proved existence on three strata
(axis-parallel, `>= 5` axis-parallel, the whole far field `>= 25.8431°` by centre pigeonhole vs Graham's `d_6`), and found that the
surviving clause — three exactly axis-parallel, consecutively separated squares — must separate packings from margin `-0.36 t^2`
configurations as `t -> 0`, which counting, pigeonhole, unavoidable points and Menger/Dilworth all lose at first order
(obstruction X).  `notes/t3-chord.md`: Kearney–Shiu contains exactly one Farkas certificate; its existence step is unavoidability;
chord rows give a second-order-correct certificate (Prop. CC3) but the pigeonhole substitute loses `(T-2)t/(T-1)` at every `T`.
`T = 4` inherits X one order flatter and loses the far-field pigeonhole (`4 - sqrt2 = 2.586 > 1/d_12 = 2.572`).  Kept: Lemma CC
(chain counting proved, `notes/counting-ladder.md`), the ladder closed form `H1 = -eps^3/(2(T+1))` (a description of the LP dual,
depth `2, 3, 2` at `T = 3, 4, 5`, `search/BREAK_H1.md`), and one open door: a **13-point** pure unavoidable set for `[0,4]^2`
(needed by the K–S slack-1 accounting; `<= 12` excluded by `COVER >= 12.2688`, 13 open, an integrality question).  The far-field mass cap is also
refuted in the strongest family (`search/FARFIELD_STRONG.md`: exact `11.999999870`, `11.999999855` at `eps = 1°, 2°` with the cap slack; no cover-side
instrument at any `eps`), so no region of the architecture is proved by machine either; note the chord row in that family's
table rows above was flagged as a hypothesis (`T4LEAF.md` §1.2) but is **proved** for closed-disjoint unit squares at `t <= 4`
(`lean/Sqpack/Chord.lean` `wall_strip_le_three`, `notes/chord-lemma.md`; only `T4LEAF`'s `1e-9` coordinate tolerance is argued outside Lean) — the flag in
`T4LEAF.md`, `FARFIELD_STRONG.md` §7 and `review-2026-09-22.md` is stale (caught by `notes/n12-gap.md` §7).

**13-point door tested (2026-09-22, `notes/unavoid13.md`, `notes/unavoid13-no.md`; briefs `tasks/unavoid13`, `tasks/unavoid13-no`).**  Hitting-set
cutting-plane loop: a finite family `F` of closed unit squares, the hitting-set IP over the arrangement cells (`h(F)` is a lower bound on the pure
unavoidable number when the candidates are all cells), exact verification of the IP's 13-sets, violated poses fed back.  `T = 3` **[proved]**: the
minimum pure unavoidable set for `[0,3]^2` is exactly **7** (575-square family with `h = 7`, exact incidence, three solvers; K–S set with `a = 9/10`
certified) against fractional `5.53–5.99` — integrality gap `>= 1.0`.  `T = 4` **[proved]**: no 13-point set is invariant under the 180° rotation
(hence none under `C4`, `V`, the diagonal Klein group, `D4`); single reflections open (`symRx` running at ~2 h/round, `Rd` never started).  `T = 4`
general: **undecided** — `h(F) = 13` at every solved round up to `|F| = 1,761` (34k exact cells, LP `12.25`; the finite cover LP over the `COVER4`
support is exactly `12.268804`), the largest family (2,809 squares) undecided after 2 h of MIP (the 1,617-square one later gave a 13-set, `unavoid13-no.md` §0 update), worst violation of polished 13-sets
down to `-0.016`, every 13-set found one skeleton up to `D4` (columns `3+3+4+3` at `x ≈ 0.90, 1.40, 2.25, 3.10`), maximin never above `-0.028`.
Evidence leans "no 13-set" (the `T = 3` gap transported says 14; asymmetry forced) with no proof; reaching `h(F) = 14` would need a prover that does
not pay per cell (row-pruned lazy IP, or geometric region branching with column generation), days of compute, outcome uncertain.  Side results:
`zeromargin.py` cannot certify K–S-type sets (the witness point switches along a curve inside every box) — fixed by a two-witness **SEG** primitive
(`search/unavoid13_check.py`, proved in `unavoid13.md` §2.2, unit weights only); a `1e-7` incidence tolerance silently weakened the first day's `T = 4`
LPs (`11.8` vs `12.27` exact).  The complete negative record for `n = 12` is `notes/n12-gap.md` (2026-09-22).

Architecture (2026-09-20 night, `notes/proof-architecture.md` §0a): the reduction, chain counting (`> 9` near-axis squares
force a wall-to-wall chain of four) and the local theorem at the tiling are all true at every `T` — the tiling is a strict
second-order maximum with `c_T = (T-2)/(2(T-1))` rising to `1/2` (`search/ARCH_TLEDGER.md`) — so none of them carries the
theorem.  What does: `4 <= k <= 9` near-axis squares among tilted ones (where Cleemann-type bands live at large `T`), and
the far field `k <= 3`, which a near-axis mass cap on the degree-1 LP does not close at small `eps` (`search/ARCH_FARFIELD.md`).

Band-cut work (2026-09-21, `notes/proof-architecture.md` §0a items 11–14, log `notes/review-2026-09-20c.md`): the record `T = 12`
and `T = 11` packings verify with our rows in exact arithmetic (`delta = +6.81e-4`, `+2.75e-4`).  Mixed-tilt chain inequality MT
proved (`notes/bandcut-cost.md` §1); the length accounting of a band has no threshold in `T`.  Measured (`search/BANDCUT_K.md`): no
chain-free configuration with `delta* >= 0` at `T = 3, 4`; chain-free maximum `-0.1 eps^3` at `T = 4` with 8 or 9 near-axis squares.
**Candidate theorem: every configuration has a tight wall-to-wall chain of `T` squares, tilted or not, on which MT gives `<= 0`** —
false at `T = 11`, where band chains spread their link normals.  In waste language Claim(T) is `w_min(T,T) >= T + 1`.

`T = 3` rehearsal of the chain architecture (2026-09-21, `notes/t3-chain.md`, Opus; log `notes/review-2026-09-21.md`): the H-lemma
(chain of `T` + `L` transverse links) has the closed form `H <= 0` iff `L >= L*(T,t)`, `L*(T,0) = T-2` [proved]; it fails at
`T = 3` above `28.8°` (only `L = 2` available; the H there is the three-square certificate `(10-7 sqrt2)/2 > 0` (`t3-existence.md` §7; legs are worth nothing at `45°`)).  What works in
the far field is a **cycle**: a closed loop of separations winding once, pinned at the four walls, no rise term — exact at
`45°`, `delta <= (93 - 66 sqrt2)/7` [proved]; the true dual is cyclic in 539/720 samples.  MT is a bound, not a certificate
(it reads the rise off the configuration); H ∪ MT covers 716/720.  Kearney–Shiu is not a chain proof.  `s(6) = 3` reduces to
two unproved existence clauses (chain with `L >= L*`, or cycle); all `T`-dependence is in them, none in the inequalities.
`search/T11_CHAINS.md`: at `T = 11` five tight chains of eleven, all through the band; existence holds, **rise control fails**
(band link nets `+0.12`, interface `-0.17`, three band links break even — the whole `+2.7e-4`).

`T = 4` cycle census (2026-09-21 night, `search/T4_CYCLES.md`, Opus, 727 + 11 optima): **chain-or-cycle is false at `T = 4`** (37 failures, 36 in
the chain-free cells).  The `-0.1 eps^3` cell (`k = 8, 9`) is held by a **tree** — a nested ladder: x-chain of four through two tilted
squares at `1/5`, y-legs at `tan(eps)/5`, an x-link on a leg's foot at `tan^2(eps)/5`; each level buys a power of `eps` and the ladder
cannot be truncated.  Best cycle bound there `+3.4e-3`, four orders wrong.  What holds 727/727 is `H1` = H + one further main-direction
row.  At `45°` H is exactly tight (`delta* = (7 sqrt2 - 10)/2`, `L* = 8 sqrt2 - 7 < 8`): the `T = 3` large-tilt failure does not transfer.
The `t3-chain.md` §5.3 cycle formula reproduces the LP at 1/17 cyclic duals and is **false without its axis-aligned-turn hypothesis**
(`-0.14` at a margin-`0` point); that hypothesis is nearly empty at `T = 4`.  Far field: 98 % cyclic, but `mu >= 2` at 104/126.

Sanity filter for any proposed mechanism (2026-09-20, corrected that night): `s(n^2 - n) < n` for **all `n >= 12`**
(Arslanov–Mustafin–Shangitbayev, Electron. J. Combin. 28(4) 2021, P4.22; refereed) and at `n = 11` (Cantrell, Feb 2025,
unrefereed; both verified here in exact arithmetic, `search/BANDCUT_SCAN.md` §2), so the conjecture can only hold for `n <= 10`;
and `s(n^2 - 4) < n` at `n = 3`.  An argument that
works uniformly in `n` is wrong; a proof must use the chain length `4` quantitatively.

Rules of thumb: explore on the packing side, certify on the cover side; nothing below `t = 4` transfers
to a proof of `s(12) = 4`; launch anything over 10 minutes detached (`setsid nohup`); never `pkill`;
16 threads is the machine's useful width.  Leave no `until … sleep` / `tail -f` watchers behind, and say so in agent
briefs; a `pgrep -f <name>` wait loop matches its own command line and never exits.

**`s(21)`, the `m = 5` rung (2026-09-23; `search/S21_{KILL,COVER,OVERHEAD,LB}.md`, `search/M4_MARGIN.md`).**  Certified: **`s(21) >= 5000/1001 = 4.995005`** (`certificates/s21/`, Rust verifier + `xcheck.py`; previous best `4.7438`, Friedman DS7) and `nu_f^closed(5) >= 20.647781` (exact; float `~20.65`), so unlike `n = 12` the fractional route is not killed (`L - n = -0.352` vs `+0.269`).  Heuristic: converged cover LP `20.73`, best honest cover `21.016` (strict scan).  Overhead anatomy at `m = 4`: the `×1.05` was almost all of it and scales with area; split into validity (LP weights `<= 0.97154` at unsampled near-tile poses, `θ ≲ 1/D`) and checker margin — the exact checker certifies the whole `[0,4]^2` at `0.25 %` margin (`W = 12.732248`, 40 column sweeps, not a certificate of record), its floor `~0.18 %` set by CHAIN at tile poses.  So `s(21) = 5` needs validity `<= 1.0105` over the LP; sampled rows and the closing loop do not get there (`m = 4`: `1.029`; `m = 5`: `1.0138`).  Not reachable as things stand.  Weight-per-cell extrapolation (estimate) gives room `~1.0` at `m = 6`, `~1.9` at `m = 7`: `s(32)`, `s(45)` may fit even the current overhead.

Decided not worth it: interior/quadrant branching on symmetric leaves; squeezing `3.975–3.978` from the
corner level; SDP / Lovász-ϑ in the continuum; pure `3.99` to convergence; line-chord cuts at `t = 4`;
any further count or incidence-pattern tree at `t = 4`; testing new inequalities axis-parallel first;
level-2 branching by wall total; the `s(12) >= 3.98` consolation prize (2026-09-13: an exact result is
what matters).
