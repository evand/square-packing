# t4-cycles: which certificates the `T = 4` duals actually use, and what holds the thinnest cell shut

*2026-09-21.  Task `tasks/t4-cycles/README.md`.  Code (new; nothing in `search/` modified —
`t3_chain`, `t3_chain_shape`, `bandcut_k`, `s6skel`, `s6local` are imported):
`search/t4_cycles.py` (the dual census, the tight-separation graph, the pure-cycle bound, the
`sec 5.3` closed form with its hypotheses checked, the ladder, the filter),
`search/t4_cycles45.py` (the uniform-tilt line and the `45 deg` point in full).
Runs: `runs/t4_cycles_census1.jsonl` (516 recorded optima) + `runs/t4_cycles_census_report.txt`,
`runs/t4_cycles_fresh1.jsonl` (282 fresh) + `runs/t4_cycles_fresh_report.txt`,
`runs/t4_cycles_cell89.txt`, `runs/t4_cycles_45.txt`, `runs/t4_cycles_uniform.{txt,json}`,
`runs/t4_cycles_uniform_ladder.txt`, `runs/t4_cycles_cf.txt`, `runs/t4_cycles_dich.txt`,
`runs/t4_cycles_h1_cyconly.txt`, `runs/t4_cycles_filter_key.txt`, `runs/t4_cycles_filter.txt`.
Labels: **[proved]** = an identity checked by hand and re-checked in float arithmetic here (the
Farkas residual of the weight vector is printed); **[measured]** = produced by a multistart or by
reading one LP dual, hence a feasible point and a LOWER bound on the quantity it estimates;
**[heuristic]**; **[guess]**.*

---

## 0. Verdict, up front

1. **The chain-or-cycle dichotomy is FALSE at `T = 4`.**  At `37` of `727` sampled optima with
   `delta* < 0` neither the best H-restricted dual nor the best pure-cycle dual reaches `0`.
   `36` of the `37` are in the chain-free `(k, eps)` cells of `BANDCUT_K.md`, one is a coherent
   small-tilt (`coh1`) sample; the uniform-tilt line adds a `38`th at `t = 35 deg`, where
   `H = CYC = +3.2e-03` against `delta* = -4.2e-02`.  **What does hold, at `727 / 727`, is one
   notch further out: `H1` = the H plus ONE further main-direction pair row.**  Its support
   contains the H's, so `H1 <= H` whenever the H-optimal chain admits any extra main-direction
   row — which it did at every point measured here.  It repairs all `37` failures (`37/37`), and
   it also covers the `32` optima where only the cycle worked (`32/32`).  **[measured]**  (§7.)

2. **The chain-free `k = 8, 9` cell — the thinnest margin anywhere, `-0.100 eps^3` — is held shut
   by a TREE, not a cycle, and not by anything close to a cycle.**  At `eps = 1 deg`
   (`delta* = -5.26e-07`) the dual is a tree at every weight threshold and at row tolerances `1e-7`
   and `1e-10`; the best bound over all `123` simple cycles of the tight graph is **`+3.39e-03`**,
   i.e. `O(eps)` and four orders of magnitude on the wrong side.  The dual is a **three-level
   nested ladder** — a wall-to-wall x-chain of four through two tilted squares at weight
   `1/(T+1)`, a transverse y-chain at weight `tan(eps)/(T+1)`, an x-link on the leg's foot at
   `tan^2(eps)/(T+1)` — and it **cannot be truncated** (re-solving on the heavy rows alone gives
   `+inf`).  Each level buys one power of `eps`:
   `CW = +0.39 eps`, `H = +0.196 eps^2`, `H1 = -0.099 eps^3 = delta*`.
   **The chain route has been chasing the right object here.**  **[measured]**  (§3.)

3. **At `45 deg`, `T = 4`, the H is exactly tight and no cycle is needed.**
   `delta*((45 deg)^12, T = 4) = (7 sqrt2 - 10)/2 = -0.050252531694167` (15 digits), and
   `L*(4, 45 deg) = 8 sqrt2 - 7 = 4.3137 < L_max = 8`, so
   `H(4, 5, 45 deg) = (7 sqrt2 - 10)/2 = delta*` **exactly** — the `T = 3` large-tilt failure of
   Lemma H **does not transfer to `T = 4`**.  The optimum is degenerate: a `4`-cycle with `q = 2`
   turns, a straight diagonal chain of five with two walls at each end, and a `mu = 2` object are
   all optimal duals of the same value.  A **cycle of eight** does exist (`~1500` of them) and
   does give `-0.0502525`, so `t3-chain.md` §3.4's guess is realised — but it is **not
   necessary**: `k >= q sin(alpha/2)(T - u)` reads `k >= 3.66` for `q = 2`.  **[proved]** for the
   certificates, **[measured]** for `delta*`.  (§4.3, §4.4.)

4. **The `t3-chain.md` §5.3 cycle formula reproduces the LP value at exactly ONE of the sampled
   `T = 4` optima.**  Of `17` optima whose dual support is exactly one simple cycle, the formula
   (link weight `1`, turn weight `2 sin(alpha/2)`) matches the LP at `1` — the `45 deg` point —
   and is `+0.014` to `+0.172` at the other `16`, whose `delta*` runs from `0` to `-6.14e-02`
   (so at the two margin-`0` ones it is on the right side of the filter — just uselessly far).  The reason is the **uniform link weight**: a `T = 4` "cycle" dual is almost always
   an H whose crossbar and legs have closed into a loop, with the two arcs at weights in ratio
   `tan(eps)` (`10^2` apart at `eps = 1 deg`).  **[measured]**  (§4.2.)

5. **Filter: passed by every Farkas bound, and it catches a real bug in the cycle lemma.**  At the
   `(4,1)` band stack (`runs/bandcut_scan_pT_T4.json`, `delta = 0` exactly) and the `p = 4` bar,
   `CW = H = HP = CYC = -2.8e-17`, i.e. `0`, and the §5.3 closed form **with its hypotheses
   checked** is exactly `0` — tight, not negative.  But the same formula with the
   **"every turn is axis-aligned" hypothesis dropped** returns **`-0.1399`** there, on a `6`-cycle
   whose weight vector has Farkas residual `1.0` — it cancels nothing.  **That hypothesis is
   load-bearing and `t3-chain.md` §5.3 is false without it.**  Separately it is nearly *empty* at
   `T = 4`: for a square of tilt `t` the jump between its own two edge normals is axis-aligned
   only at `t = 45 deg`.  **[measured]**  (§6.)

6. **The shape of the dual, by regime.**  Where a chain of four exists and nothing forbids it
   (`51 / 58` unrestricted `k >= 4` optima are `mu = 0`, and `48` of those `51`, plus the band
   stack) the dual is exactly five rows — three links and two walls, each at weight
   `1/(T+1) = 0.2`.  In the far field
   (`k <= 3`) `124 / 126` duals are cyclic, i.e. **`98 %` of the sampled far field**, but `104` of
   those have `mu >= 2`, median `13` pair rows on `11` squares: a thicket, not a cycle.  In the
   chain-free cells the ladder of §3 dominates.  **[measured]**  (§2.)

7. **Caveat that governs all of the above: the optimal dual is not unique.**  At `45 deg` three
   local optima of the same value gave three structurally different optimal duals (`mu = 1` with
   2 walls, `mu = 0`-flavoured with 4, `mu = 2`).  `mu` is a property of the HiGHS vertex, not of
   the configuration.  Conclusions here are stated in terms of *restricted-dual bounds* (`H`,
   `CYC`, `H1`), which are vertex-independent, wherever possible.

---

## 1. What is measured, and what each object is

Semantics and row system exactly as `search/t3_chain.py` (= `BANDCUT_K.md` §1.1 =
`S6_SKELETON.md` §3.1): the **wall-carrying** LP in the `2n` centre coordinates at fixed angles,
every row of the form `a . c - delta >= b`, walls `lo/hi-x/y(i)` with `b = P_i`, `P_i - T`, pair
rows `(i,j,o,kind,s)` with `b = m_ij`.  A **Farkas certificate** is `w >= 0` supported on rows
the configuration satisfies, with `sum_r w_r a_r = 0` and `sum_r w_r = 1`; it proves
`delta <= -sum_r w_r b_r` for **every** configuration satisfying those rows.  `T = 4`, `n = 12`
throughout.

Five bounds are computed at each sampled optimum.  All but the last are restricted duals, i.e.
Farkas bounds, and therefore `>= delta*` by construction — **they cannot be below `0` at a
margin-`0` point, and that is why the filter of §6 is a test of the closed FORM, not of the LP.**

| name | support it is allowed |
|---|---|
| `CW` | one wall-to-wall chain of `T` (its `T-1` links) + **every** wall row |
| `H` | the H-lemma: that chain + every wall row + every **transverse**-type pair row |
| `HP` | `H` + the main-direction rows between squares *of the chain* (staircase closure) |
| `H1` | `H` + **one** further main-direction pair row (`t3-chain.md` §2.5's minimal repair), best over all choices |
| `CYC` | the pair rows of **one simple cycle** of the tight separation graph + every wall row |
| `CF` | the closed form of `t3-chain.md` §5.3 on an explicit cycle, **with its hypotheses checked** |
| `CFnaive` | the same formula with the hypotheses **not** checked (not a bound; see §6) |
| `F` | everything the configuration satisfies (`= delta*` at a global optimum) |

**The tight graph.**  By complementary slackness the optimal dual is supported on rows tight at
the optimum, so the cycles that can carry a certificate of value `delta*` live in the graph whose
edges are pairs carrying a **tight** pair row.  That graph has `25-50` edges on `12` vertices at a
typical `T = 4` optimum, and `10^2` to `>2 x 10^5` simple cycles of length `<= 9`; when there are
more than a few thousand a uniform random sample is taken (fixed seed), so every `CYC` and `CF`
number here is a **lower bound on what a cycle can do** — the honest direction, since a
certificate search can only fail to find one.  Two restrictions to keep in mind when reading a
*negative* cycle result: the cycles searched are those of the **tight** graph (complementary
slackness forces that only for a certificate of value exactly `delta*`; a certificate that merely
reaches `0` could in principle use a slack row), and the length is capped at `9`.  "No cycle
certifies here" therefore means "no cycle of length `<= 9` on tight rows, in this sample" — a
search result, not a theorem.

**Two caveats that bite at `T = 4` and did not at `T = 3`.**

1. **The optimal dual is not unique.**  At `theta = (45 deg)^12` three different local optima of
   the same value produced three structurally different optimal duals — a `4`-cycle with two wall
   rows, a `6`-edge object with four, and a `mu = 2` object (§4.2).  So "the shape of the dual"
   is the shape of *one vertex* of the optimal dual face; `mu` is not an invariant of the
   configuration.  Every shape statistic below is a statistic of the HiGHS vertex.
2. **The `eps = 1 deg` cells sit at `|delta*| ~ 5e-7`**, which is only five times the default
   row-inclusion tolerance `1e-7` of `t3_chain.all_rows`.  Every number for those cells was
   recomputed at `tol = 1e-10` (§3); the conclusions do not move, but the third digit does.

---

## 2. The dual census  [measured]

`516` recorded `T = 4` optima (`runs/bandcut_k_T4{a,free,free_pilot,_hunt}.jsonl`, the `(4,1)`
band stack of `BANDCUT_SCAN.md` §5, the two best-bar points) plus `282` fresh ones (47 angle
families x 6 seeds, `delta*` by multistart from `600` tiling starts + `200` random starts per
sample; `n = 12` multistart is weak, so every `delta*` is a lower bound).  `runs/t4_cycles_census1.jsonl`,
`runs/t4_cycles_fresh1.jsonl`; reports `runs/t4_cycles_census_report.txt`,
`runs/t4_cycles_fresh_report.txt`.

### 2.1 Three regimes, and they look completely different

| regime | `N` | `mu = 0` | `mu = 1` | `mu >= 2` | median `E`, `V` | `H <= 0` | `CYC <= 0` |
|---|---|---|---|---|---|---|---|
| **dense / hole, unrestricted** (`k >= 4`, no chain cut) | 58 | **51** | 6 | 1 | `3`, `4` | 58/58 | 39/58 |
| **far field** (`k <= 3`, plus `bigrand`/`generic`/`near20`) | 126 | **2** | 20 | 104 | `13`, `11` | **126/126** | 45/126 |
| **chain-free cells** (`k = 4..9`, the chain among near-axis squares forbidden) | 444 | 49 | 228 | 167 | `9`, `9` | 377/444 | 240/444 |
| the `(4,1)` band stack and the two bars | 3 | 3 | 0 | 0 | `3`, `4` | 3/3 | 3/3 |

Read off:

1. **Where a chain of four exists and nothing forbids it, the dual is the bare chain and nothing
   else.**  `51 / 58` unrestricted `k >= 4` optima have `mu = 0`, and `48` of those `51` — and the
   `(4,1)` band stack — have a dual of exactly **five rows**: `T - 1 = 3` pair rows and two wall
   rows on one axis, **each at weight exactly `1/(T+1) = 0.2`** (`E = 3, V = 4`; the other three
   are `(4,5)` and `(5,6)` staircases).  `54 / 58` of the `58` have `delta* = 0`.  This is
   `S6_SKELETON.md` §1.2's signature, unchanged at `T = 4`.
2. **In the far field the dual is not a cycle — it is a thicket.**  `124 / 126` far-field duals
   are cyclic (`mu >= 1`), so **`98 %` of the sampled far field is cyclic** — but the median is
   `13` pair rows on `11` squares with `mu >= 2` in `104 / 126`.  A single simple cycle with one
   wall row per turn is the *exception*, not the rule: it occurs in `20 / 126`.  At `T = 3`
   (`t3-chain.md` §3.1) the far-field median was `4-7` pair rows on `4-6` squares with
   `mu = 1` dominant; at `T = 4` the object has roughly doubled in size and gained loops.
3. **The chain-free cells are where the interesting duals live**, and they are neither bare chains
   nor single cycles: a *nested ladder* of alternating-axis chains of four whose weights fall by a
   factor `tan(eps)` per level (§3).

### 2.2 Which walls, and how many turns

Among the `444` chain-free-cell duals and all `126` far-field duals, **every** support uses at
least one row on each of the four walls (`4walls` column of the report is `N/N` in every such
group); the only supports that use fewer are the bare chains of §2.1(1), which use two wall rows
on one axis.  So "one wall row per wall" does not distinguish a cycle from an H at `T = 4`: an H
uses all four walls too (`lo-ax`, `hi-ax` for the crossbar, `lo-tr`, `hi-tr` for the legs).  The
statistic that does distinguish them is `mu`, and `mu` is **not an invariant** (§1, caveat 1).

**Turns.**  The brief asks for the number of turns.  "Turn" is only defined once the support has
been oriented as a closed loop, i.e. only for a dual whose pair rows form one simple cycle; there
are `17` such among the `798` (§4.2), with `q = 2` to `9`.  For the other `781` — bare chains,
ladders, thickets — the notion does not apply, and the substitute statistic is the number of wall
rows, tabulated above.

### 2.3 The uniform-tilt line

Every `unif` family sampled here (`1, 5, 10, 20, 28.07, 30, 45 deg`) is `CYC <= 0` at `6/6` and
`H <= 0` at `6/6` (except `unif20`, `5/6`).  **That is an artefact of which tilts were sampled**:
the finer sweep of §5 finds that both fail at `35 deg` and that `H` fails at `40 deg`.  §4.4 says
why the `T = 3` large-tilt failure of Lemma H nevertheless does **not** reappear at `45 deg`.

---

## 3. The `eps^3` cell: `k = 8, 9` chain-free is held shut by a TREE, not a cycle  [measured]

`BANDCUT_K.md` §3(b) reports the thinnest margin anywhere: at `k = 8` and `k = 9` the chain-free
maximum is `-0.100 eps^3` (`-5.26e-7` at `eps = 1 deg`).  `search/t4_cycles.py cell89`
(`runs/t4_cycles_cell89.txt`) reads the dual there and widens the support one geometric fact at a
time.  Every row is re-run at `tol = 1e-10` as well as `1e-7`; the numbers below are the
`tol = 1e-10` ones where they differ.

### 3.1 The dual is a three-level nested ladder

At `k = 9`, `eps = 1 deg` (`delta* = -5.259969e-07`; eight near-axis squares at `0.000 deg`, four
far squares at `1.000 deg`) the dual has exactly these weights:

| weight | rows | what it is |
|---|---|---|
| `0.19724638` | `lo-x(1)`, `pair(1,4)`, `pair(4,11)`, `pair(7,11)`, `hi-x(7)` | a wall-to-wall **x-chain of four**, `1 - 4 - 11 - 7`, whose two middle squares `4, 11` are FAR (tilt `eps`) |
| `0.00344` | `lo-y(3)`, `pair(3,4)`, `pair(7,8)`, `hi-y(8)` | two **y-legs**, one down from chain member `4`, one up from chain end `7` |
| `0.00006008` | `pair(3,6)`, `hi-x(6)` | one **x-link on the leg's foot**, `3 - 6`, closing on the far `x` wall |
| `1.8e-7` | `lo-y(0)`, `pair(0,1)` | a fourth level |

The weight ratio between consecutive levels is `tan(eps)` to three digits:
`0.0034419 / 0.19724638 = 0.017450` and `0.00006008 / 0.0034419 = 0.017456`, against
`tan(1 deg) = 0.0174551`.  Same at `eps = 5 deg` (`0.0868, 0.0875` vs `tan 5 = 0.08749`).  That is
Lemma H's leg weight `tan t` (`t3-chain.md` §1.3), applied **twice**.

**Is it a tree or a cycle?**  At `eps = 1 deg`, `mu = 0` at weight threshold `1e-6` **and** at
`1e-9`: a **tree**, `E = 7, V = 8` (`k = 9`) / `E = 6, V = 7` (`k = 8`).  At `eps = 5 deg` and
`10 deg` the vertex HiGHS returns has `mu = 1` and `mu = 3`, but those cells are `100x` and
`1000x` thicker in `|delta*|` and are not the cell the brief is about.  **The cell that carries
the `eps^3` is a tree.**  The heavy part alone does not certify: re-solving on the rows of weight
`> 1e-6` only returns `+inf` (no cancelling combination exists) — **the ladder cannot be
truncated**.

### 3.2 The best H bound, the best cycle bound, and what carries the `eps^3`

| `eps` | `CW` (chain + walls) | `H` | `H1` = H + **one** extra main-direction link | best cycle `CYC` | `delta*` |
|---|---|---|---|---|---|
| `1 deg` | `+6.872e-03` | `+5.955e-05` | **`-5.285e-07`** | `+3.394e-03` | `-5.260e-07` |
| `5 deg` | `+3.226e-02` | `+1.363e-03` | **`-5.805e-05`** | `+1.514e-02` | `-6.189e-05` |
| `10 deg` | `+5.956e-02` | `+4.903e-03` | **`-4.062e-04`** | `+2.605e-02` | `-4.596e-04` |

divided by the right power of `eps` (radians):

| `eps` | `CW / eps` | `H / eps^2` | `H1 / eps^3` | `CYC / eps` | `delta* / eps^3` |
|---|---|---|---|---|---|
| `1 deg` | `+0.394` | `+0.1955` | `-0.0989` | `+0.1945` | `-0.0989` |
| `5 deg` | `+0.370` | `+0.1790` | `-0.0874` | `+0.1735` | `-0.0931` |
| `10 deg` | `+0.341` | `+0.1610` | `-0.0764` | `+0.1493` | `-0.0864` |

Four things follow, and they answer the brief's second question.

1. **A cycle does not hold this cell shut, and is not close.**  The best pure-cycle bound is
   `O(eps)` — the same order as the bare chain-plus-walls — and is four orders of magnitude
   (`+3.4e-3` against `-5.3e-7`) on the wrong side at `eps = 1 deg`.  `123` simple cycles of length
   `<= 9` exist in the tight graph there; none of them helps.  **The chain route has not been
   chasing the wrong object here.**
2. **Nor does the H.**  `H = +0.1955 eps^2`; note that `0.1955 -> 1/5 = 1/(T+1)`, the far-field
   constant of `BANDCUT_K.md` §0(2).  The H is short by exactly one order of `eps`.
3. **What carries the `eps^3` is one extra main-direction link on top of the H** — the third rung
   of the ladder, `pair(3,6)` at weight `tan^2(eps)`.  `H1` recovers `delta*` to `0.5 %` at
   `eps = 1 deg` (`-5.285e-07` vs `-5.260e-07`) and to `6-12 %` at `5, 10 deg`.  Each level of the
   ladder buys exactly one power of `eps`: `O(eps) -> O(eps^2) -> O(eps^3)`.
4. **That extra link is NOT a rung and does NOT make a cycle.**  `t3-chain.md` §0(4) reads "one
   extra edge on a tree is one cycle"; at `T = 4` in this cell the extra edge arrives with its own
   wall row and a new square, so it is a **new branch**, and the support stays a tree.  The
   `T = 3` repair and the `T = 4` repair are the same move on the row list and different objects
   in the graph.

**The object that holds the chain-free `k = 8, 9` cell shut is therefore: a wall-to-wall chain of
`T` through two tilted squares, plus a transverse chain at weight `tan eps`, plus a second
transverse chain at weight `tan^2 eps` — Lemma H with its two closing wall rows replaced by two
more chain links.**  It is exactly the three-level object `S6_LOCAL.md` §3 reports for the
`n = 12` pinwheel (`1/3`, `t/3`, `t^2/3` there; `1/5`, `tan(eps)/5`, `tan^2(eps)/5` here, the
difference being hard containment vs. the wall-carrying LP).

---

## 4. The cycle lemma at `T = 4`

### 4.1 The arithmetic, and which of its hypotheses is load-bearing

Restating `t3-chain.md` §5.3: round a cycle `Q_1 ... Q_k` with link normals `n_1 ... n_k` at
weight `1`, the Farkas residual at vertex `s` is `n_{s-1} - n_s`, so the wall rows there must
supply `d_s = n_s - n_{s-1}`.

> **[proved] The lemma's "assume the `q` wall rows split evenly between `lo` and `hi`" is not an
> assumption — it is automatic.**  Summing the link coefficient vectors round a *closed* cycle
> telescopes to `0`, so the wall coefficients sum to `0` **axis by axis**; `lo-x` weight equals
> `hi-x` weight and likewise in `y`.  (Verified: `lo = hi = sqrt2` in the `45 deg` certificate
> below.)

> **[measured] The lemma's "every turn is axis-aligned, cancelled by ONE wall row" is a
> `45 deg`-only hypothesis, and it is where the lemma nearly dies at `T = 4`.**  For a square of
> tilt `t` the two edge normals are `n_0 = (cos t, sin t)` and `n_1 = (-sin t, cos t)`, and
> `n_0 - n_1 = (cos t + sin t, sin t - cos t)` is axis-aligned **only at `t = 45 deg`**.  So at a
> common tilt `t != 45 deg` no turn *within one square's own normal pair* can be cancelled by one
> wall row.  Of the `17` sampled `T = 4` optima whose dual support is exactly one simple cycle,
> the hypothesis is satisfiable at **`1`** — the `45 deg` point (`runs/t4_cycles_cf.txt`).

The repair is to allow **two** wall rows at a turn, one per axis; the turn weight is then the
**`L1`** norm `|d_x| + |d_y|` of the jump rather than its Euclidean norm `lambda_j`
(`cycle_closed_form2`).  That version has no hypothesis, applies to every cycle, and was checked
to be an exact Farkas certificate at all `17` (max residual `5.8e-12`).

### 4.2 Does the formula reproduce the LP value?  Once out of seventeen.  [measured]

`runs/t4_cycles_cf.txt`.  For each `T = 4` optimum whose dual support is exactly one simple cycle,
the generalised closed form (link weight `1`, turn weight `|d_x| + |d_y|`) against the LP value:

| point | `k` | turns `q` | wall rows | closed form | LP dual |
|---|---|---|---|---|---|
| `unif45` | 4 | 2 | 2 | **`-5.025253e-02`** | `-5.025253e-02` |
| `k=5 eps=1 chainfree` | 8 | 4 | 8 | `+7.184e-02` | `-1.178e-04` |
| `k=4 eps=1 chainfree` | 10 | 6 | 12 | `+1.099e-01` | `-8.961e-05` |
| `unif20` | 8 | 4 | 8 | `+1.429e-01` | `-2.558e-02` |
| `hole3_20` | 9 | 7 | 14 | `+1.449e-02` | `-6.138e-02` |
| `coh5` | 10 | 9 | 18 | `+9.898e-02` | `-2.249e-03` |
| (11 more) | 6-9 | 4-7 | 7-14 | `+4.5e-02 .. +1.72e-01` | `0 .. -3.36e-02` |

**`1 / 17`.**  The reason is not the turn hypothesis — it is the **weight `1` on every link**.
When the dual support of a `T = 4` optimum happens to be a single cycle, it is almost always an H
whose crossbar and legs have closed into a loop, and the LP puts weight `~1/(T+1)` on the
crossbar arc and `~tan(eps)/(T+1)` on the leg arc — a ratio of `10^2` at `eps = 1 deg`.  A cycle
with uniform link weights is then a **different and much weaker** certificate.

> **So the answer to the brief's question 3 is: no.  The `t3-chain.md` §5.3 formula reproduces the
> LP value at exactly one of the sampled `T = 4` optima, the `45 deg` point, and nowhere else.**

### 4.3 `45 deg` exactly  [measured `delta*`, proved certificate]

`search/t4_cycles45.py deep45` (`runs/t4_cycles_45.txt`), `theta = (45 deg)^12`, `T = 4`,
`2000` starts:

    delta*  =  -0.050252531694168   =   (7 sqrt2 - 10) / 2   (all 15 digits)

The optimum is degenerate and **three different local optima produced three different optimal
duals of the same value**:

1. a **`4`-cycle** `6 - 8 - 9 - 7` with `q = 2` turns, link weight `1`, two wall rows
   (`lo-x(6)`, `hi-x(8)`) of weight `sqrt2`.  The `t3-chain.md` §5.3 closed form applies
   (`axis_ok`, Farkas residual `1.1e-16`) and gives
   `[ (1/2)(2 sqrt2)(4 - sqrt2) - 4 ] / (4 + 2 sqrt2) = (4 sqrt2 - 6)/(4 + 2 sqrt2) = (7 sqrt2 - 10)/2`
   exactly.  **[proved]**
2. a **straight diagonal chain of five** squares, all four links carrying the same normal, with
   **two** wall rows at each end (weight `1/sqrt2` each) because the end residual `+-n` is not
   axis-aligned.  Four wall rows, total wall weight `2 sqrt2` — the same as (1).
3. a `mu = 2` object with `7` pair rows and four wall rows.

All three have `omega = (sqrt2 - 1)/2` and Lemma W returns the same number.

**Does a cycle of eight exist at `45 deg`?  Yes.**  The tight graph has `1472-1545` simple
`8`-cycles, and the best of them, with `q = 4` axis-aligned turns, gives **exactly**
`-0.050252532` — the same value again.  So the brief's `t3-chain.md` §3.4 guess is *realised*,
but it is **not necessary**: the general condition is `k >= q sin(alpha/2) (T - u)`, which at
`T = 4, 45 deg` reads `k >= 7.31` only for `q = 4`; for `q = 2` it reads
`k >= sqrt2 (4 - sqrt2) = 3.66`, i.e. **a cycle of four**, and that is what the `T = 4` optimum
actually uses.  *The "cycle of eight" guess assumed `q = 4` because `T = 3` needed one wall row
per wall; at `T = 4` the optimum pins itself against only two walls, both on the same axis.*

### 4.4 Why there is no `T = 4` analogue of the `T = 3` large-tilt failure  [proved + measured]

`t3-chain.md` §0(3): `H <= 0` needs `L >= L*(T, t)`, and at `T = 3, 45 deg` only `L = 2` links are
available against `L* = 2.485`, so the H fails and only a cycle works.  At `T = 4`,
`L*(4, 45 deg) = 8 sqrt2 - 7 = 4.3137` and `L_max = n - T = 8`, so `L = 5` suffices:

| `L` | 0 | 1 | 2 | 3 | 4 | **5** | 6 | 8 |
|---|---|---|---|---|---|---|---|---|
| `H(4, L, 45 deg)` | `+0.4983` | `+0.3431` | `+0.2171` | `+0.1127` | `+0.0248` | **`-0.0502525`** | `-0.1151` | `-0.2213` |

and `H(4, 5, 45 deg) = (7 sqrt2 - 10)/2 = delta*` **exactly**.  Measured directly at the `45 deg`
optimum: the H-restricted dual returns `-5.025253e-02` (`runs/t4_cycles_uniform_ladder.txt`).
*(Pretty coincidence: `H(3, 2, 45 deg) = (10 - 7 sqrt2)/2` is the exact negative of this.)*

> **Corollary (why "chain or cycle" is the wrong frame at uniform tilt).  [proved]**  By Lemma W
> (`t3-chain.md` §1.4) every normalised certificate at a common tilt has value
> `omega (T + 2 - u) - 1`, a strictly increasing function of the wall mass `omega` **alone**.  A
> cycle of `k` links with `q` turns of weight `lambda` has `omega = q lambda / (2(k + q lambda))`;
> an H with `L` legs has `omega = M / ((T-1) + L tan t + 2M)`.  At `T = 4, 45 deg` the `4`-cycle
> with `q = 2`, the `8`-cycle with `q = 4` and the H with `L = 5` all give
> `omega = (sqrt2 - 1)/2`.  **They are the same certificate in three costumes; there is no
> mechanism dichotomy here, only a competition for the smallest `omega`** — i.e. for the most pair
> rows per unit of wall weight.

---

## 5. The uniform-tilt line, and the ladder  [measured]

`search/t4_cycles45.py uniform` (`runs/t4_cycles_uniform.txt`, `.json`) plus the ladder
(`runs/t4_cycles_uniform_ladder.txt`).  `delta*` by multistart at `theta = t^12` — a lower bound.

| `t` | `delta*` | `CW` | `H` | `H1` | best cycle | dual shape |
|---|---|---|---|---|---|---|
| `1 deg` | `-6.0077e-05` | `+6.81e-03` | `-5.80e-05` | `-6.004e-05` | `-5.90e-05` | `mu=1`, `E=9,V=9`, 4 walls |
| `5` | `-1.4191e-03` | `+3.09e-02` | `-1.18e-03` | `-1.398e-03` | `-1.27e-03` | `mu=1`, `E=9,V=9` |
| `10` | `-5.2779e-03` | `+5.49e-02` | `-3.56e-03` | `-4.984e-03` | `-4.04e-03` | `mu=2` |
| `15` | `-1.1012e-02` | `+7.32e-02` | `-5.79e-03` | `-9.692e-03` | `-6.61e-03` | `mu=2` |
| `20` | `-1.8088e-02` | `+8.70e-02` | `-6.86e-03` | `-1.439e-02` | `-7.26e-03` | `mu=2` |
| `28.0725` (the band angle) | `-3.0966e-02` | `+1.02e-01` | `-4.38e-03` | `-2.266e-02` | `-4.38e-03` | `mu=2` |
| `30` | `-3.4072e-02` | `+1.04e-01` | `-2.83e-03` | `-2.500e-02` | `-2.83e-03` | `mu=1` |
| **`35`** | `-4.1655e-02` | `+1.09e-01` | **`+3.21e-03`** | `-3.167e-02` | **`+3.21e-03`** | `mu=1` |
| `40` | `-4.7644e-02` | `+1.12e-01` | **`+4.04e-03`** | `-3.976e-02` | `-4.26e-03` | `mu=1` |
| `44` | `-5.0129e-02` | `+1.13e-01` | `-5.68e-04` | `-4.790e-02` | `-4.07e-02` | `mu=2` |
| `45` | `-5.0253e-02` | `+1.13e-01` | `-5.0253e-02` | `-5.0253e-02` | `-5.0253e-02` | `mu=1`, `E=4,V=4`, 2 walls |

Three things.

* **`H` is not monotone in `t` and it does fail at `T = 4` — at `35 deg` and `40 deg`, not at
  `45 deg`.**  The failure is not the `T = 3` failure (there `L*` outran `L_max`); here `L* = 3.89`
  at `35 deg` with `L_max = 8`, so the closed-form inequality is satisfied with room.  What fails
  is **existence**: the `35 deg` optimum does not contain a wall-to-wall chain of four carrying
  four transverse links at a common tilt with a common normal — `26` chains of four exist there
  and the best H-restricted dual over all of them is `+3.2e-03`.
* **At `35 deg` the cycle fails too** (`+3.21e-03`), so this is a `uniform`-family counterexample
  to the chain-or-cycle dichotomy that the `(k, eps)` census does not contain.
* **`H1` — H plus one further main-direction pair row — is `<= 0` at every tilt on the line**,
  and tracks `delta*` to within a factor `1.3` above `20 deg`.

---

## 6. The filter  [measured]

The brief's rule: at a point with `delta* = 0` no chain or cycle bound may be `< 0`.
`search/t4_cycles.py filter --only 'band stack' 'bar'` (`runs/t4_cycles_filter_key.txt`) on the
`(4,1)` band stack of `BANDCUT_SCAN.md` §5 (`runs/bandcut_scan_pT_T4.json`, `delta = 0` exactly,
a genuine closed packing) and the `p = 4` anchored bar (also `delta = 0`).  The remaining
`53` recorded optima with `delta* >= -1e-12` are covered by the census scan of §2, which computes
`H` and `CYC` at every one of them, and by the tolerance re-check below;
`runs/t4_cycles_filter.txt` is the same sweep over all of them (slow: `~1` min per point, because
every cycle of the tight graph gets both an LP and a closed-form enumeration).

### 6.1 Everything that is a Farkas bound passes, as it must

| point | `delta*` | `F` | `CW` | `H` | `HP` | `CYC` | `CF` (hyps checked) | `CFnaive` |
|---|---|---|---|---|---|---|---|---|
| `(4,1)` band stack | `0` | `+5.6e-17` | `-2.8e-17` | `-2.8e-17` | `-2.8e-17` | `-2.8e-17` | `-0.0e+00` | **`-1.399e-01`** |
| `p = 4` bar | `0` | `+5.6e-17` | `-2.8e-17` | `-2.8e-17` | `-2.8e-17` | `-2.8e-17` | `-0.0e+00` | **`-1.377e-01`** |

`CW`, `H`, `HP` and `CYC` are restricted duals, hence Farkas bounds, hence `>= delta* = 0` by
construction; they come out at `-2.8e-17`, i.e. `0`.  `CF` — the `§5.3` closed form **with its
hypotheses checked and its Farkas residual verified** — also comes out at exactly `0`, tight.
**The cycle lemma as stated passes the filter at `T = 4`.**

Three "violations" at the `1e-9` level appeared in the census scan (`k=6 eps=10 free`,
`k=9 eps=1 free`, `k=11 eps=1 free`; `H` between `-2.2e-9` and `-7.8e-9`).  These are the
row-inclusion tolerance, not bounds: re-solved at `tol = 1e-10` and `1e-12` the same three give
`H = -6.9e-17, -5.6e-16, -1.7e-16` and `CYC = +inf, -5.6e-16, -1.1e-16`, i.e. `0` (verified,
this file).  **Nothing below `-1e-8` was produced by any Farkas bound at any margin-`0` point.**

### 6.2 What the filter does catch: the formula without its hypothesis

**`CFnaive = -0.1399` at the `(4,1)` band stack.**  That is the `§5.3` formula applied to the
`6`-cycle `[2, 9, 5, 7, 4, 11]` with `q = 4` turns, none of them axis-aligned, using the Euclidean
jump `lambda_j = |n_j - n_{j-1}|` as the single wall weight.  The resulting weight vector has
Farkas residual `1.0`, i.e. **it does not cancel the centre coordinates at all** — it is not a
certificate, and it "proves" `delta <= -0.14` at a packing of margin exactly `0`.  The same thing
happens at `45 deg`, where the naive form returns `-0.347` against a true `delta*` of `-0.0503`.

> **Conclusion.  The "every turn is axis-aligned" hypothesis of `t3-chain.md` §5.3 is
> load-bearing, not cosmetic: dropping it makes the lemma FALSE at `T = 4`, and the `(4,1)` band
> stack is an explicit counterexample.**  Any use of the cycle lemma must check, turn by turn,
> that the normal jump is axis-aligned — or use the `L1` generalisation of §4.1, which needs two
> wall rows per turn and is correspondingly weaker.

---

## 7. The dichotomy, tested  [measured]

`798` sampled `T = 4` optima (`516` recorded + `282` fresh), of which `727` have `delta* < 0`,
plus the `11` uniform-tilt points of §5.

| statement | holds at |
|---|---|
| `H <= 0` | `658 / 727` |
| best cycle `<= 0` | `365 / 727` (`CYC` is a *sampled* minimum over cycles of length `<= 9` on tight rows, so this under-counts) |
| **`H <= 0` or cycle `<= 0`** (the candidate lemma) | **`690 / 727`** |
| `H1 <= 0` (H + one extra main-direction link) | **`727 / 727`** |

The `37` failures of the dichotomy (`runs/t4_cycles_dich.txt`): `36` are in the chain-free
`(k, eps)` cells of `BANDCUT_K.md` (`k = 4..9`, every `eps`), `1` is a `coh1` sample (the coherent
small-tilt cone, which is also where `T = 3` failed).  The uniform-tilt line adds one more at
`t = 35 deg` (§5) that the `(k, eps)` census does not contain.  **At every one of the `37` the
minimal repair is the same and it works `37/37`: one further main-direction pair row**, and at
`19 / 37` `H1` recovers `delta*` to better than `1 %`.

At the `32` further optima where the cycle saves a failing `H`, `H1` also works (`32/32`,
`runs/t4_cycles_h1_cyconly.txt`) — **so nothing in this sample needs a cycle that the ladder
cannot do**.  `H1`'s support contains the H's, so `H1 <= H` whenever the H-optimal chain admits
an extra main-direction row; that held at every point measured here, and no point was found with
`H1 > H`.

### The lemma the data supports

> **Conjecture (measured, `T = 4`, `n = 12`, `798 + 11` sampled optima).**  Let `n = T(T-1)` unit
> squares sit in `[0,T]^2` with margin `delta < 0`.  Then there is a wall-to-wall chain of `T`,
> together with transverse pair rows and **one further main-direction pair row**, whose Farkas
> combination gives `delta <= 0`.  Measured `727 / 727`; no counterexample.  The `T = 3` version
> of the same statement is `t3-chain.md` §2.5 (`68/68` there), where the extra row closes a cycle;
> at `T = 4` it usually does not.

*What this does not say.*  It is an existence statement about the separation graph of a
hypothetical packing, exactly like `t3-chain.md`'s (O1)/(O2)/(O3), and nothing here proves it.
What the measurement adds is that **the third object the architecture needs is not a cycle**: the
two branches "chain with enough transverse links" and "cycle" do not cover `T = 4`, and one more
level of the transverse ladder does, in this sample.

### Two-line verdict

> **No.  The chain-or-cycle dichotomy is false at `T = 4`: at `37 / 727` sampled optima with
> `delta* < 0` — all but one of them in the chain-free `k = 4..9` cells, plus the uniform
> `35 deg` point — neither the best H nor the best cycle reaches `0`, and the thinnest cell of all
> (`k = 8, 9`, margin `-0.1 eps^3`) is held shut by a TREE whose best cycle bound is `+0.0034`,
> four orders of magnitude on the wrong side.**
>
> **What does hold at `727 / 727` is one notch further out: `H1`, the H plus one further
> main-direction pair row — which at `T = 3` closed a loop and at `T = 4` usually grows a new
> branch instead.  `98 %` of the sampled far field (`k <= 3`: `124 / 126`) is cyclic (`mu >= 1`),
> but it is a thicket, not a cycle: `mu >= 2` at `104 / 126`, median `13` pair rows on `11`
> squares.**

---

## 8. What in the brief's premises turned out to be wrong

* **"Is the chain-free `k = 8, 9` cell held by a cycle — if so the chain route has been chasing
  the wrong object there."**  It is **not**.  The dual there is a tree (`mu = 0` at both weight
  thresholds, at `tol = 1e-7` and `1e-10`), the best pure-cycle bound over all `123` cycles of the
  tight graph is `+3.4e-03` against `delta* = -5.3e-07`, and what carries the `eps^3` is a third
  rung of the transverse ladder at weight `tan^2 eps`.  The chain route is chasing the right
  object; it just needs one more level of it.  (§3)
* **"A cycle of eight at `45 deg`" (`t3-chain.md` §3.4).**  A cycle of eight does exist and does
  certify, but the guess derived `k >= 8` from `q = 4` turns, carried over from `T = 3` where the
  cycle is pinned by one row per wall.  The actual `T = 4` `45 deg` optimum uses a **cycle of
  four with `q = 2` turns**, pinned by two wall rows on the *same* axis, and an equally valid
  optimal dual is a straight diagonal **chain of five** with two walls at each end.  All three
  have the same `omega` and the same value `(7 sqrt2 - 10)/2`.  (§4.3)
* **"Chain H fails for tilts above `28.8 deg`, so the far-field certificate must be a cycle."**
  That is a `T = 3` fact and it does **not** transfer.  At `T = 4`, `L*(4, 45 deg) = 8 sqrt2 - 7 =
  4.31` against `L_max = 8`, and `H(4, 5, 45 deg) = (7 sqrt2 - 10)/2 = delta*` **exactly** — the
  H is *tight* at `45 deg`.  Where the H does fail on the uniform line is `35 deg` and `40 deg`,
  and for a different reason (chain **existence**, not the closed-form inequality).  (§4.4, §5)
* **Erratum to `t3-chain.md` §0(4), found here.**  "One extra edge on a tree is one cycle" is
  false as a general statement about this repair: at `T = 4` the extra main-direction row arrives
  with its own wall row on a new square, so the support stays a tree.  (§3.1)
* **A hole in `t3-chain.md` §5.3, found here.**  Its "each turn is axis-aligned" hypothesis is
  **load-bearing**, and the lemma is FALSE without it: at the `(4,1)` band stack (`delta = 0`
  exactly) the formula applied to a `6`-cycle with non-axis-aligned turns returns `-0.1399`.
  Separately, the hypothesis is nearly **empty** at `T = 4` — of the `17` sampled optima whose
  dual is exactly one simple cycle it is satisfiable at `1`, because for a square of tilt `t` the
  jump between its own two edge normals is axis-aligned only at `t = 45 deg`.  Its
  "assume the wall rows split evenly between `lo` and `hi`" clause, on the other hand, is a
  theorem, not an assumption.  (§4.1, §6.2)
* **Not wrong, but worth recording:** `search/bandcut_k.py chain_certificate` was not used here
  (it is the known-buggy function); nothing in `search/` was modified.

---

## 9. Reproduce

    # the dual census over the 516 recorded T = 4 optima          (~50 s on 12 processes)
    python3 search/t4_cycles.py census --out runs/t4_cycles_census1.jsonl --nproc 12
    python3 search/t4_cycles.py report runs/t4_cycles_census1.jsonl

    # 282 fresh samples, 47 angle families x 6 seeds             (~3 min on 12 processes)
    python3 search/t4_cycles.py fresh --reps 6 --limit 200 --extra 200 --nproc 12 \
        --out runs/t4_cycles_fresh1.jsonl
    python3 search/t4_cycles.py report runs/t4_cycles_fresh1.jsonl

    # the eps^3 cell, the ladder, both row tolerances             (~15 min)
    python3 search/t4_cycles.py cell89

    # the uniform-tilt line and 45 deg in full                    (~5 min / ~15 min)
    python3 search/t4_cycles45.py uniform --out runs/t4_cycles_uniform.json
    python3 search/t4_cycles45.py deep45 --keep 3

    # the dichotomy failures, re-tested with H1
    python3 search/t4_cycles.py dich runs/t4_cycles_census1.jsonl runs/t4_cycles_fresh1.jsonl

    # the filter                                                  (~40 min)
    python3 search/t4_cycles.py filter
    python3 search/t4_cycles.py filter --only 'band stack' 'bar'   # the two points that matter
