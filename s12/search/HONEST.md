# honest-cost: what the best pure `t = 4` dual actually costs as a cover (2026-09-12)

Task: `tasks/honest-cost/README.md`.  Code: `search/honestcost.py` (+ a `--dump-dual` hook in
`search/cliquelever.py`; nothing in `verify/`, `lean/`, `certificates/` is touched).  Runs in this
worktree's `runs/hc_*` (gitignored; every number is quoted here).  Read against
`search/RANKDIAG.md` §8 and "Verdict after rounds 2-5", `search/CLIQUELEVER.md` §0-1,
`search/RUNG2.md` §0/§4.3/§10, `search/BOXCLIQUE.md`.

## 0. Verdict, up front

**The honest cost is `18.1` (E2Pg) and `20.2` (the certified support state) — not `12.1-12.6`.  The
`t = 4` rank family is nowhere near a cover, so V1 (the disjunctive `t = 4` verifier with cliques
and polygons) is NOT worth building on this dual: the next work is the family.**  The brief's
decision rule (`honest >= 12.3` ⇒ the family is short) fires with a factor of `1.5` to spare, and
it fires under every variant of the measurement I could construct.

| dual | `Theta` = LP | min capture, clique rule `meet` | **honest** | min capture, clique rule `core` | honest |
|---|---|---|---|---|---|
| `E2Pg` (14,611 poses, the full loaded set) | `11.864926892` | **`0.655199`** | **`18.109`** | `0.073029` | `162.47` |
| `E2Pg`, a second dual of the same optimum | `11.864926892` | `0.622941` | `19.047` | — | — |
| `SUPEPGF` (956-pose certified support) | `11.800418089` | **`0.585278`** | **`20.162`** | `0` | `inf` |
| `SUPEPGF` + one repair round (1,312 poses) | `11.835811134` | `0.735603` | `16.090` | — | — |

Robustness of the headline `18.109`: `16.79` if every candidate pose is forced `>= 1e-6` off the
grid lines (so no containment tolerance can be doing the work), `15.78` even with an
indefensibly generous `1e-6` slack on the clique-meet test.  The smallest number anywhere in this
measurement is `15.8`.

**Three findings, in order of how much they change the picture.**

1. **The generous clique rule is not merely generous — it is unsound.**  The brief's primary rule
   (`cliquelever.price`'s: a pose is credited `z_K` iff it closed-meets every member of `K`) does
   not define a valid weighting, because `{S : S meets every member of K}` need not be a clique.
   Measured (`honestcost.py sound`): on `SUPEPGF`, **304 of the 317 clique rows carrying dual admit
   two pairwise-DISJOINT admissible squares that both credit the row**, and those rows carry
   `6.782753` of the `7.352132` clique dual — `92 %`; on `E2Pg`, `330` of `638` rows, carrying
   `3.375913` of `7.610316` — `44 %` (and that is a lower bound: only 200 credited poses per row
   from a `0.04 / 2.5 deg` lattice are sampled).  Example (`SUPEPGF` row 37, `z = 0.206960`, 38
   members): `(0.77388, 2.46231, 67.5 deg)` and `(1.40880, 3.33271, 72.5 deg)` each closed-meet all
   38 members and are disjoint from each other.  So the `meet` numbers above are *optimistic*:
   `min capture` under `meet` is an upper bound on the min capture of any sound rule, hence
   `honest` under `meet` is a **lower** bound on the honest cost.  `18.1` and `20.2` are floors.
2. **Every clique row that carries dual is non-Helly, so the verifier's own primitive credits it
   nothing.**  `BOXCLIQUE.md`'s box clique is credited by `clique_of_cores`: a pose gets `w_K` iff
   it contains a core, cores pairwise meeting.  For a row whose "boxes" are single poses — all a
   member list gives you; the checkpoints record no boxes and no cores — the core is the common
   intersection of the member squares, and on both duals **all** of them are empty (`317/317` on
   `SUPEPGF` carrying `7.352132`; `638/638` on `E2Pg` carrying `7.610316`).  Re-derived
   independently by a max-margin LP over the member half-planes: the deficits are `-0.0007` to
   `-0.042`, a genuine non-Helly gap, not a rounding artefact.  Under the sound `core` rule the
   clique mass credits nobody at all off the pose set, `min capture` reaches exactly `0` (an
   axis-parallel square at `(1.5, 2.5)` captures *nothing* of `SUPEPGF`'s dual — its whole point
   mass is 67 atoms), `honest = inf` there and `162.47` on `E2Pg`, whose 403 atoms and 44 polygon
   rows leave every square at least `0.073029`.
   The two rules do not "differ materially"; they differ maximally, and nothing in
   between exists in the repo today.  Cliques being non-Helly is the whole reason they are a lever
   (`CLIQUELEVER.md` §0, `CLIQUE.md`) — and it is exactly what makes them uncreditable.
3. **The `0.04 / 2.5 deg` pose lattice hides the worst poses by a factor of six.**  On `E2Pg` the
   lattice's min capture is `0.945215` (honest `12.553`), which reproduces the pricer's
   reduced-cost floor of `0.05-0.07` (`RANKDIAG.md` §17) to three decimals — but
   `search/family_rows.py`'s item-3 family at pitch `0.004` finds `0.655199` (honest `18.109`), a
   reduced cost of `0.345`: five to seven times the floor the pricer has reported for eight
   consecutive injections.  The minimising pose is **the corner square of the trivial `4 x 4`
   tiling, nudged `1e-7` off the grid lines**: `(1.5000001, 3.4999999, 0 deg)`, the closed square
   `[1.0000001, 2.0000001] x [2.9999999, 3.9999999]`.  At `(1.5, 3.5, 0)` exactly, capture is
   `1.045`; `1e-7` to the right it is `0.655`, because the square then misses every atom on the
   line `x = 1` (the dual has 51 of them, mass `0.420`; `0.158` of that lies in this square's
   `y`-range) and, separately, clique rows worth `0.232`.  This is `RUNG2.md` §4.3 one level up,
   and it is the same fix: `family_rows.py` already emits these poses, and the pricer does not
   look at them.

## 1. What is computed

For a dual `w = (y, z, zeta)` of the packing LP of `CLIQUELEVER.md` §1 + `RANKDIAG.md` §8 on a
fixed pose set `P`:

```
Theta(w)     = sum_p y_p  +  sum_K z_K  +  sum_Pi (k-1)/2 * zeta_Pi          (= the LP value)
capture_w(S) = sum_{p in S} y_p  +  sum_K z_K [S credited K]
                                 +  sum_Pi zeta_Pi [S credited Pi]
honest(w)    = Theta(w) / min over ALL admissible poses S of capture_w(S)
```

Semantics, unchanged from the brief: `t = 4`, closed unit squares, closed containment, a point on
`dQ` counts (`ZEROMARGIN.md` §1).  A pose `(cx, cy, theta)` is admissible iff
`cx, cy in [w(theta)/2, 4 - w(theta)/2]`, `w(theta) = |cos| + |sin|`; `theta in [0, 90)` suffices
because a *unit square* is invariant under a quarter turn (that is not a container symmetry).  Both
pose files say `sym 1` and both runs pin nothing (`--corners .... --patterns ........`, so
`leaf_ceiling.region_targets` returns `None` and the LP has no region rows and no chord rows), so
the whole domain is scanned and no `D4` reduction is used anywhere.

**Crediting rules.**
* points: `p` counts iff `p in S` (closed).
* polygons: `Pi` with anchors `A_0..A_{k-1}` counts iff `S` contains both endpoints of some side
  `[A_i, A_{i+1}]` (`rankfamily.py`).  This rule is **sound** over the continuum — the lemma of
  `RANKDIAG.md` §8 is proved for arbitrary admissible squares.
* cliques, rule `meet` (the brief's primary, and `cliquelever.price`'s): `K` counts iff `S`
  closed-meets every member of `K`.  Generous, and unsound — §0 item 1.
* cliques, rule `core` (the brief's "stricter box/anchor-clique rule"): `K` counts iff `S` contains
  the core `C_K = intersection of the member squares`; rows with `C_K = {}` credit nobody.  Sound
  (two squares that both contain a non-empty `C_K` meet), and it is `BOXCLIQUE.md`'s rule with the
  boxes degenerated to single poses.

**The duals.**  They are not on disk, so each state's LP is re-solved on its **fixed** pose set
from its own checkpoint (`cl_TAG_poses.txt` as the only `--exact`, `_rows.txt`, `_cliques.txt`,
`_pgons.json` resumed, every clique row re-verified pairwise and every polygon row re-derived from
its exact rational anchors), with the restricted master OFF so the dual is feasible for **every**
loaded column:

| state | poses / columns | coverage rows | clique rows | polygon rows | LP re-solved | `Theta - LP` | max rc over the loaded columns |
|---|---|---|---|---|---|---|---|
| `SUPEPGF` | 956 / 962 | 14,881 | 1,006 | 746 | `11.800418089` (5 s) | `-3.6e-15` | `+3.1e-14` |
| `E2Pg` | 14,611 / 14,700 | 46,964 | 1,563 | 2,002 | `11.864926892` (1,110 s) | `+7.2e-13` | `+6.5e-13` |

`SUPEPGF`'s `11.800418089` is its certified measure to nine places (`RANKDIAG.md` §17's
`11.800417910 = 1180041791/100000000`).  `E2Pg`'s `11.864926892` is `0.000001` *below* the
`11.865961` quoted in `RANKDIAG.md` §17, which is expected rather than a discrepancy: the
checkpoint is written **after** the last iteration's separation, so it carries the 65 coverage rows
and 19 clique rows that iteration added but never solved with.  More rows, lower value.

**Dual degeneracy is real but does not matter here.**  These LPs are very degenerate (`CLMASTER.md`:
"a vertex dual of this very degenerate LP is a useless pricing signal").  Solving `E2Pg` the other
way — closing the restricted master exactly (`honestcost.py dual --master`, 24 min instead of
19 min) — lands on the **same** objective `11.864926892` with a visibly **different** dual (points
`3.808723` vs `3.709598`, cliques `7.511191` vs `7.610316`, 631 vs 638 rows with dual).  Its honest
cost is `19.047` (min capture `0.622941`), i.e. slightly worse, at the **same** minimising pose
`(1.5000001, 3.4999999, 0 deg)`, with the **same** clique credit there (`0.465143` from 43 rows, to
six places) and less point credit (`0.148462` against `0.180720`); its `0.04` lattice minimum is
`0.945215040`, identical to nine places to the first dual's.  So the number reported is the one
HiGHS returns, as the brief allows, and the non-uniqueness moves it by `1` in `18`, not by `6`.

**Self-test (step 1).**  `Theta` reproduces the LP to `1e-12` on both states, and under `meet`

    min capture over the 956 poses of SUPEPGF's P   = 1.000000000  (0 poses below 1 - 1e-9)
    min capture over the 14,611 poses of E2Pg's P   = 1.000000000  (0 poses below 1 - 1e-9)

which is the crediting rule being right: on the loaded poses the `meet` rule credits **exactly** the
LP row's members (checked row by row — for `SUPEPGF`'s clique 0, 32 members, 0 members not credited
and 0 non-members credited: the rows really are maximal over the loaded set), and the polygon rule
reproduces `rankfamily.row_members` on 40 of 40 rows tested.  Max capture over `P` is `1.42`, mean
`1.014`, as a degenerate optimum should look.

Under `core` the same self-test gives `min capture = 0.029972` on `SUPEPGF` and `0.082091` on
`E2Pg`, with **every** pose of `P` below 1 — the direct statement of §0 item 2: the sound rule does
not even certify the pose set the LP was solved on.

## 2. The numbers

### 2.1 `E2Pg` — the full loaded set (LP `11.864926892`)

`Theta` = points `3.709598` (403 atoms with `y > 0`) + cliques `7.610316` (638 rows with `z > 0`,
793 to 1,806 members, median 1,152) + polygons `0.545013` (44 rows).  **`64 %` of the total is
clique mass; `4.6 %` is polygon mass.**

| candidate set | poses | min capture | honest |
|---|---|---|---|
| `0.04 / 2.5 deg` lattice, full domain | 169,689 | `0.945215` | `12.553` |
| `family_rows.py` item 3 + tilted, pitch `0.004` | 61,656 | **`0.655199`** | **`18.109`** |
| the 14,611 poses of `P` (self-test) | 14,611 | `1.000000` | `11.865` |
| the 400 heaviest coverage-row points as centres, 36 angles | 14,400 | `0.968330` | `12.253` |
| pattern search, 400 seeds, steps `0.04 -> 1.5e-5` | 400 | `0.655199` | `18.109` |
| knife-edge refinement, `+-eps` at `1e-4 .. 1e-9` | 100 | `0.655199` | `18.109` |
| dense `0.02 / 1 deg` lattice | 1,689,197 | `0.801930` | `14.795` |
| pattern search + knife edge from the dense worst 400 | 400 | `0.781059` | `15.191` |

The descent does not move off `0.655199393` at any step size from `0.04` down to `1.5e-5`, nor
under `+-eps` probes at `1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9` in all three coordinates: a plateau,
not a knife edge (see the profile in §3).  The 1.69M-pose dense pass at `0.02 / 1 deg` gets to
`0.801930` — better than the coarse lattice's `0.945215` and still `0.15` above the family's
`0.655199`, which is the point of §3: refining a lattice does not reach a `1e-6`-wide dip.

The `core` rule on the same candidate sets: lattice `0.099474`, families `0.092799`, the poses of
`P` themselves `0.088385`, descent + knife edge `0.073029` (honest `162.47`), at
`(2.50023, 1.50093, 0.141 deg)` — `0.030124` from 7 point atoms, `0.042906` from 5 polygon rows,
and **`0` from all 638 clique rows**.

### 2.2 `SUPEPGF` — the certified 956-pose support (LP `11.800418089`)

`Theta` = points `3.276397` (only 67 atoms in the whole container) + cliques `7.352132` (317 rows,
27-131 members) + polygons `1.171889` (40 rows).  **`62 %` clique mass.**

| candidate set | poses | min capture | honest |
|---|---|---|---|
| `0.04 / 2.5 deg` lattice | 169,689 | `0.591071` | `19.964` |
| `family_rows.py`, pitch `0.004` | 61,656 | `0.861534` | `13.697` |
| the 956 poses of `P` | 956 | `1.000000` | `11.800` |
| row points as centres | 14,400 | `0.901772` | `13.086` |
| pattern search from the 400 worst | 400 | **`0.585278`** | **`20.162`** |
| knife-edge refinement `1e-4 .. 1e-9` | 200 | `0.585278` | `20.162` |
| dense `0.02 / 1 deg` lattice | 1,689,197 | `0.585278` | `20.162` |

Here the dense scan and the descent agree to nine places: on a 956-pose support the dual is short
almost everywhere, so a lattice has no trouble finding the bottom.  Under `core`, `min capture` is
exactly `0` at `(1.49999, 2.50001, 89.995 deg)` — an axis-parallel square in the middle of the left
half of the container that captures none of the 67 atoms and none of the 40 polygon rows.

### 2.3 The minimisers, and what they capture

**`E2Pg`** — the corner square of the trivial tiling, `1e-7` off the lines:

| pose | capture | points | cliques | polygons |
|---|---|---|---|---|
| `(1.5000001, 3.4999999, 0 deg)` | `0.655199` | `0.180720` (35 of 403 atoms) | `0.465143` (43 of 638 rows) | `0.009336` (4 of 44) |
| `(1.5000002, 3.488, 5.2e-5 deg)` | `0.655199` | `0.180720` (35) | `0.465143` (43) | `0.009336` (4) |
| `(3.4999999, 1.5000001, 0 deg)` | `0.667162` | `0.153877` (27) | `0.492585` (50) | `0.020701` (5) |
| `(1.5000001, 0.5000001, 0 deg)` | `0.669365` | `0.157444` (34) | `0.460436` (37) | `0.051485` (7) |
| `(1.5, 3.476, 5.7e-5 deg)` | `0.706515` | `0.180720` (35) | `0.516459` (50) | `0.009336` (4) |

The atoms it does see are all on `y = 3` and `x = 2` (`(2.0, 3.0486)`, `(1.5409, 3.0004)`,
`(1.6415, 3.0004)`, `(1.7886, 3.0)`, ...) — the dual's point mass sits on the wall-square corner
lines `x, y in {1, 2, 3}`, exactly where `RANKDIAG.md` §5 says the binding pentagons wrap — and a
square whose left edge has just left `x = 1` sees only the `y = 3` half of it.

The cleanest way to see how short this dual is: take the **16 unit squares of the trivial `4 x 4`
tiling of `[0,4]^2`** and evaluate `capture` at each, first exactly and then nudged `1e-7` (the best
of the admissible sign choices).  Rows are `cy = 3.5` down to `cy = 0.5`, columns `cx = 0.5 .. 3.5`:

| exactly on the tiling | | | | | nudged `1e-7` | | | |
|---|---|---|---|---|---|---|---|---|
| `1.0000` | `1.0451` | `1.1045` | `1.0000` | | **`1.0000`** | **`0.6552`** | **`0.7220`** | **`1.0000`** |
| `1.0337` | `2.0368` | `1.9334` | `1.0915` | | `0.6638` | `0.8335` | `0.8460` | `0.7247` |
| `1.1354` | `2.1108` | `2.0498` | `1.0527` | | `0.7396` | `0.9524` | `0.8352` | `0.6672` |
| `1.0000` | `1.0572` | `1.1089` | `1.0000` | | `1.0000` | `0.6694` | `0.7334` | `1.0000` |

Exactly on the tiling every square is captured `>= 1` — the four container corners at exactly
`1.0000`, the four interior cells at `1.93-2.11`, which is the dual paying twice for the poses it
knows about.  Nudge by `1e-7` and **12 of the 16 fall below 1**, as low as `0.6552`; only the four
corner squares survive, and they survive because there is nowhere for them to go (every admissible
nudge keeps the two container-boundary edges).  That is the whole story of this dual in one table.

**`SUPEPGF`** — a `45-55 deg` square pressed against the right wall, halfway up:

| pose | capture | points | cliques | polygons |
|---|---|---|---|---|
| `(3.29894, 1.42701, 52.50 deg)` | `0.585278` | `0.089295` (6 of 67) | `0.428359` (15 of 317) | `0.067624` (6 of 40) |
| `(3.30311, 1.41326, 54.75 deg)` | `0.598958` | `0.103103` (7) | `0.428359` (15) | `0.067496` (5) |
| `(3.29300, 1.48803, 46.00 deg)` | `0.604210` | `0.054595` (5) | `0.448672` (17) | `0.100943` (10) |

Its six atoms are again on the corner lines (`(3.0779, 1.0)`, `(3.1082, 1.0)`, `(3.0, 1.7671)`, ...)
and it sees `0.089` of the `3.276` point mass.  The *bulk* of the low-capture region, though, is
**interior**: of the 2,000 worst lattice poses, 292 have their centre within `0.25` of `(1.5, 1.5)`,
275 of `(2.5, 2.0)`, 258 of `(2.0, 1.5)`, spread evenly over all 36 angles.  A 956-pose support
simply does not constrain the interior of the container.

## 3. The lattice blind spot, measured

On `E2Pg` the `0.04 / 2.5 deg` lattice bottoms out at capture `0.945215`, i.e. reduced cost
`0.054785` — inside the `0.05-0.07` band `RANKDIAG.md` §17 reports the pricer plateauing in for
eight consecutive injections.  The measurement here computes the same quantity the pricer does
(`capture = 1 - rc` when nothing is pinned), on the same lattice, with the same arithmetic, so the
agreement cross-checks both.

The item-3 family finds `0.655199`: reduced cost `0.344801`, **6.3 times** the lattice floor.  The
profile along `cx`, at `cy = 3.5`, `theta = 0` (capture is piecewise constant, so these are the
cells, not samples):

| `cx - 1.5` | `0` | `1e-9 .. 9e-9` | `1.5e-8 .. 8e-7` | `1.3e-6 .. 5e-5` | `7e-5 .. 2e-4` | `3e-4` | `4e-4` | `>= 7e-4` |
|---|---|---|---|---|---|---|---|---|
| capture | `1.045` | `0.813` | **`0.655`** | `0.895` | `0.907` | `0.926` | `0.976` | `1.000` |

so the dip is about `1e-6` wide in the free coordinate and `0.345` deep, and the pricer's grid step
is `0.04` — four orders of magnitude too coarse.  The `1.5e-8` step is where the square's left edge
clears `packing_dual.capture`'s `1e-8` containment tolerance and the `x = 1` atoms drop out; the
`1.3e-6` step is where it starts closed-meeting a new cluster of clique members (pose centres are
snapped to a `1e-6` grid, `Dc = 1e6`, so grazing contacts are dense at exactly that scale).

**Is this a tolerance artefact?**  No.  Forcing every candidate at least `delta` off every grid
line, `delta = 1e-6 … 1e-3` (i.e. rebuilding the item-3 family with a bigger admissibility margin),
the minimum is `0.706515` — `honest = 16.794` — at `(1.5, 3.476, 5.7e-5 deg)`, and it does not move
between `delta = 1e-6` and `delta = 1e-3`.  And going the other way, allowing an
*indefensibly* generous `1e-6` slack on the clique-meet test (which would credit a square that
genuinely misses a member by up to `1e-6`; the verifier's `cores_meet` is exact in `i128`), the
family minimum is `0.751679` — `honest = 15.785`.  Every variant is far above `12.3`.

Two consequences worth writing down.

* **`E2Pg`'s descending LP sequence is not approaching the continuum value.**  `RANKDIAG.md` §17
  reads the rc floor as "the pricer never runs out of improving columns at this pitch".  The
  stronger statement is available now: the pricer is not looking where the improving columns are.
  Column generation on these instances should price `family_rows.py`'s families alongside the
  lattice, exactly as `closed4.py --seed-rows` does on the cover side.
* **`RANKDIAG.md` §16's E1 test is intact, but "the discretisation is not what binds" needs a
  footnote.**  Both lattices in that test returned the identical best reduced cost because both are
  lattices; the family that beats them both by `0.29` is on neither, at any pitch.

## 4. Step 4: one round of cover-side repair

`honestcost.py repair` adds the worst poses found by the scan as **columns** of the packing LP —
which is what "adding them as rows of the cover" means on the dual side: a new column forces
`capture >= 1` there — extends every clique row to the new poses by the loop's own greedy rule,
re-derives every polygon row from its exact anchors, and re-solves.  On `SUPEPGF`, the 400 worst
poses (356 distinct after snapping to `Q = 1e5`, `Dc = 1e6`):

| | poses / columns | LP | min capture | honest | minimiser |
|---|---|---|---|---|---|
| before | 956 / 962 | `11.800418089` | `0.585278` | `20.162` | `(3.29894, 1.42701, 52.5 deg)` |
| after | 1,312 / 1,318 | `11.835811134` | `0.735603` | **`16.090`** | `(0.50158, 2.50216, 89.8 deg)` |

The clique rows absorbed `14,768` new memberships and the polygon rows `51,353` — the `3.5x` ratio
of `RANKDIAG.md`'s "why the family keeps working" table, reproduced.  The LP rose `+0.035`, the
honest cost fell `20.16 -> 16.09`, the minimiser jumped to the opposite wall and became
axis-parallel, and `Theta`'s composition shifted from cliques toward points (cliques `7.35 -> 6.86`,
points `3.28 -> 4.26`).  One round is not a trend, but the direction and the size are informative:
the repair is doing real work and there is a very long way to go.  Per the brief the loop was not
iterated further.

## 5. What is exact and what is float

**Exact**: the pose set, the coverage-row points, the clique member lists and their pairwise
verification, and the polygon anchors and their re-derivation — all `cliquelever.py` /
`rankfamily.py`'s own exact machinery re-run here (`resume_rows`, `resume_cliques`,
`rankfamily.resume`, `leaf_ceiling.sq_contains` in integers).

**Float**: everything downstream.  The LP and its dual are HiGHS doubles; `capture`'s point term is
`packing_dual.capture`'s C kernel with its `1e-8` containment tolerance; the clique `meet` test is
`cliquelever.price.clique_cost`'s 4-axis separating-axis test with a `--geo-tol` slack (default
`1e-9`; `0` and `1e-9` give identical numbers everywhere here); the polygon test is
`rankfamily.pgon_cost`'s, same slack; the cores are a float Sutherland-Hodgman clip with a `1e-9`
inflation (needed, or a point core is lost to rounding after several hundred clips) and every
non-Helly deficit was re-derived independently by a max-margin LP over the member half-planes.
Every slack is in the direction that makes `capture` **larger** and `honest` **smaller**, so the
numbers err on the side of the cover looking good.

Nothing here is certified and nothing here bounds `s(12)`.  This is a packing-side measurement of
how far a known dual is from being a weighting.

**The accounting itself.**  `Theta` bounds `sum_i capture(S_i)` over a packing `S_1..S_n` of
pairwise disjoint closed unit squares only when each row's credited set meets the packing at most
as often as its right-hand side allows: at most once for a point (the squares are disjoint — the
boundary case at side exactly 4 is handled the way `ZEROMARGIN.md` §1 handles it), at most
`(k-1)/2` times for a polygon row (its lemma, valid for arbitrary admissible squares), and at most
once for a clique row **iff the credited set is a clique**.  That last hypothesis is what §0 item 1
measures and finds false for `meet` on `92 %` of `SUPEPGF`'s and `44 %` of `E2Pg`'s clique dual.

## 6. Reproduce

```bash
# inputs (read-only), copied into a gitignored runs/
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-12/* runs/
python3 search/family_rows.py runs/hc_families.txt --pitch 0.004 --tilted    # 61,656 poses

# the duals: re-solve each state's LP on its own fixed pose set (5 s / 1,110 s)
python3 search/honestcost.py dual SUPEPGF --threads 4
python3 search/honestcost.py dual E2Pg    --threads 6
python3 search/honestcost.py dual E2Pg    --threads 6 --master --out-tag E2PgM   # the other dual
#   -> runs/hc_TAG_dual.npz, the Theta = LP check, and the capture >= 1 self-test on P

# the scan: lattice + families + descent + knife edge + a 1.69M-pose dense confirmation
python3 search/honestcost.py scan SUPEPGF --rule meet --threads 4      # ~4 min
python3 search/honestcost.py scan SUPEPGF --rule core --threads 4      # seconds
python3 search/honestcost.py scan E2Pg    --rule meet --threads 8      # 3,363 s (the dense pass)
python3 search/honestcost.py scan E2Pg    --rule core --threads 8
#   -> runs/hc_TAG_RULE_scan.json  (every candidate set's min, the minimisers' anatomy)
#      runs/hc_TAG_RULE_worst.txt  (the 400 worst poses, for repair)

# the knife-edge stage alone, from a finished scan's worst poses (minutes)
python3 search/honestcost.py fine E2Pg  --rule meet --threads 6
python3 search/honestcost.py fine E2PgM --rule meet --threads 6

# is the `meet` rule a valid cover rule at all?  (100 s / 210 s)
python3 search/honestcost.py sound SUPEPGF --ncliques 400 --sample 200
python3 search/honestcost.py sound E2Pg    --ncliques 700 --sample 200

# step 4: one repair round
python3 search/honestcost.py repair SUPEPGF --rule meet --nadd 400
python3 search/honestcost.py scan   SUPEPGFrep --rule meet
```

`cliquelever.py` also grew `--dump-dual PATH` / `--dump-dual-stop`, which writes the same npz from
inside the loop after any solve, so the honest cost of a live run can be measured without
re-solving it.

## 7. What this says about the next step

1. **A clique row is not a coverable object as it stands.**  The rows carrying dual are non-Helly by
   `0.0007-0.042`, so no `clique_of_cores` block credits them, and the extension rule the pricer
   uses is not a weighting.  Making the family creditable means either (a) splitting each row into
   box cliques whose cores pairwise meet — and the boxes can only be as wide as the pairwise
   meeting slack allows, which on these rows is thousandths, so the block count explodes — or (b)
   using a family whose membership test *is* a containment (`S` contains a segment, or a point),
   which is exactly what the odd-polygon family already is.  The polygon rows carry only `0.545` of
   `11.86` on `E2Pg` and `1.17` of `11.80` on `SUPEPGF`.  **The measurement says: move the mass out
   of cliques and into polygons.**  Every polygon row is sound and creditable by a
   `cores_meet`-style test (`RANKDIAG.md` §9); no clique row here is.
2. **Price the families, not only the lattice.**  A `0.345` reduced cost sitting outside the
   pricer's candidate set for eight injections is the largest single correctable error in the loop,
   and `family_rows.py` already emits the poses.  Until that is fixed, `RANKDIAG.md` §17's
   descending sequence should not be read as an estimate of the continuum value.
3. **`honest` is a cheap instrument and belongs in the loop.**  One LP re-solve plus minutes of
   scanning turns any checkpoint into a go/no-go number, on the same arithmetic the pricer uses.
   Run it on a candidate family *before* the verifier work, not after.
