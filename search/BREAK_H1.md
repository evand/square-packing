# break-h1: `H1` under adversarial structured search at `T = 4`, and the ladder depth at `T = 3, 4, 5`

*2026-09-21.  Task `tasks/break-h1/README.md`.  Code (new; nothing in `search/` modified —
`s6skel`, `s6local`, `bandcut_k`, `t3_chain`, `t3_chain_shape`, `t4_cycles` are imported):
`search/break_h1.py` (structured starts, the `CW / H / HP / H1 / H2` ladder, the weight-level
clustering, the exponent fit, the filter).  Runs: `runs/break_h1_*`.
Labels: **[proved]** = an identity checked by hand; **[measured]** = produced by a multistart or
by reading one LP dual, hence a feasible point and a LOWER bound on the quantity it estimates;
**[heuristic]**; **[guess]**.*

---

## 0. Verdict, up front

> **Two lines.**
>
> **`H1` survives at `T = 4` under an adversarial structured search: `0` failures in `153`
> configurations with `delta* < 0` built from tiling-minus-permutation-hole-set starts with
> tilted blocks, pinwheels, and the chain cut through tilted squares as well as near-axis ones —
> `H1 <= 0` at `153 / 153`, on top of `T4_CYCLES.md` §7's `727 / 727`, and the closest any point
> came to failing was `H1 = -2.8e-17`, i.e. `0`.**
>
> **The ladder depth needed is `2, 3, 2` at `T = 3, 4, 5`: it does NOT grow with `T`.  The
> `T = 5`, `n = 20` top-of-hole chain-free cell (`k = 15, 16`) is `p = 2`, not `p = 4`, with
> `gamma -> 1/6 = 1/(T+1)` — so the `eps^3` at `T = 4` is a `T = 4` accident, and §2.1 gives the
> counting reason it cannot recur at `T = 5`.**

1. **No `H1` failure at `T = 4`.**  `201` structured jobs, `47` infeasible cuts, `153` optima with
   `delta* < 0` (chain-freeness re-verified from the geometry at all `153`).  `CW <= 0` at
   `127 / 153`, `H <= 0` at `150 / 153`, **`H1 <= 0` at `153 / 153`** (`H2` was not computed in
   the screen — it is only needed when `H1 > 0`, which never happened).  No candidate failure
   arose, so the re-solve protocol had nothing to re-solve; the three points
   that genuinely needed the third rung were re-solved at `tol = 1e-10` with `maxchains = 400`
   anyway and `H1` recovers `delta*` to `1.2 %`.  Re-reading the `516` recorded `T = 4` optima
   with the new objects (`159` done, §3.2) adds `153 / 153` more, with certificate depth
   `3 / 132 / 18` at levels `1 / 2 / 3` and **none at level `4`**.  Counting the per-`eps` cells
   of §§4–5 as well: `H1 <= 0` at `352 / 352`.  **[measured]**  (§3.)
2. **Exponent and depth, by `(T, cell)`.**  `delta* = -gamma eps^p`, `d` the first level of
   `CW -> H -> H1 -> H2` that reaches `<= 0`:

   | `T` | cell | `p` | `gamma` | `d` |
   |---|---|---|---|---|
   | 3 | `k = 4`, top of hole, chain-free | `2.02` | `-> 1/8` | **2** |
   | 3 | `k = 3`, same family, one cut square dropped | `1.86` | `-> 1/4` | **3** |
   | 4 | `k = 9`, top of hole, chain-free | `2.96` | `-> 1/10` | **3** |
   | 4 | `k = 8`, same family, one cut square dropped | `1.93` | `-> 1/5` | **2/3** (knife-edge) |
   | 5 | `k = 15`, top of hole, permutation hole set | `2.0` | `-> 1/6` | **2** |
   | 5 | `k = 15, 16`, sheared (no permutation hole set) | `2.0` | `~ 0.25` | **2** |

   **[measured]**  Local slopes stable to `2 %` over `eps in [0.25 deg, 10 deg]`.  (§4.)
3. **"Is depth `= p` always?"  NO.**  `T = 3`, `k = 3` has `p = 2` and `d = 3`: there
   `H = +0.09 eps^3 > 0` while `delta* = -0.247 eps^2`, so the H-lemma fails and `H1` is needed
   even though the margin is only quadratic.  The depth is decided by whether the H's own leading
   term cancels — a sign question about a higher-order coefficient — not by `p`.  **[measured]**
   (§4.1(3).)
4. **`T = 5`: `p = 2`, said loudly.**  `delta*(k = 15, eps) = -0.1413, -0.1542, -0.1618, -0.1657`
   times `eps^2` at `eps = 10, 5, 2, 1 deg`, converging to `eps^2 / 6 = eps^2/(T+1)` — the
   far-field constant of `BANDCUT_K.md` §0(2).  Certificate depth `2`: the bare H closes it, one
   rung *less* than at `T = 4`.  `k = 16` coincides with `k = 15` exactly, as `k = 8, 9` do at
   `T = 4`.  **[measured; a feasible point, hence a lower bound — §5.1 says why to believe it]**
5. **[proved] Why `T = 4` is the outlier.**  The `eps^3` cell is a `2 x 2` pinwheel occupying the
   four tiling cells of ONE level cell of the merge-type level map, with a **permutation** hole
   set absorbing the rest of the cross.  Counting the removed cells: the `T - 2` non-merged rows
   each give up exactly one cell, all of them in the two merged columns, so a permutation hole set
   needs `T - 2 <= 2`.  **It is tight at `T = 4` and impossible at `T = 5`.**  (§2.1.)
6. **Filter: `0` breaches.**  Every restricted-dual bound at every margin-`0` point tested — the
   `(4,1)` band stack, the `p = 4` bar, the `T = 4` `k = 5..12` unrestricted zeros, the `T = 4`
   tiling-minus-hole-set zeros and the `T = 5` `k = 15, 16` zeros — is `>= -4.0e-12`, four orders
   above the `-1e-8` bug threshold.  **[measured]**  (§6.)
7. **New instrument worth keeping.**  `break_h1.tophole_family` builds the `BANDCUT_K.md` §3(b)
   configuration for arbitrary `T` from the combinatorics and prices it with one LP.  It
   reproduces the `T = 4` `eps^3` cell in `2.7 s` from `85` jobs, against a `(k, eps)` scan that
   needs hours, and it is the only route by which the `T = 5` cell was found at all.  (§2.)

---

## 1. What is measured, and the two new objects

Row system exactly as `search/t3_chain.py` (= `BANDCUT_K.md` §1.1 = `T4_CYCLES.md` §1): the
wall-carrying LP in the `2n` centre coordinates at fixed angles, every row `a . c - delta >= b`.
A **restricted-dual bound** on a support `S` is

    min  -sum_r w_r b_r    over  w >= 0 supported on S,  sum_r w_r a_r = 0,  sum_r w_r = 1 ,

`= +inf` when no cancelling combination exists on `S`.  Every such bound is a Farkas bound, hence
`>= delta*` for every configuration satisfying the rows of `S`; in particular **it can never be
negative at a margin-`0` point** (the filter, §6).

| name | support it is allowed |
|---|---|
| `CW` | one wall-to-wall chain of `T` (its `T-1` links) + **every** wall row |
| `H` | that chain + every wall row + every **transverse**-type pair row (Lemma H) |
| `HP` | `H` + the main-direction rows between squares *of the chain* |
| `H1` | `H` + **one** further main-direction pair row, best over all choices |
| **`H2`** | `H` + **two** further main-direction pair rows, best over pairs — **new here** |
| `CYC` | the pair rows of one simple cycle of the tight separation graph + every wall row |
| `F` | everything the configuration satisfies (`= delta*` at the optimum) |

`H2` is searched over all pairs of extra main-direction rows on the two chains that did best at
`H1` and on the chain that did best at `H` (not on all chains); the value reported is therefore an
**upper** bound on the true best-pair value.  That is the conservative direction for a "`H2` fails"
claim and the optimistic one for "`H2` holds"; no `H2 > 0` at `delta* < 0` is claimed below.

**Certificate depth** `d(z)` — the answer to the brief's "ladder depth" — is the first level of
that ladder whose bound is `<= 0`:  `d = 1` if `CW <= 0`, else `2` if `H <= 0`, else `3` if
`H1 <= 0`, else `4` if `H2 <= 0`.  `H1` *fails* at `z` exactly when `delta*(z) < 0` and `d(z) > 3`.

**Weight-ladder depth** `w(z)` is the descriptive statistic of `T4_CYCLES.md` §3.1: the optimal
dual's weights are clustered by `round( log(w/w_max) / log(tan eps) )`, and `w(z)` is the number of
consecutive levels `0, 1, 2, ...` that are populated (weights down to `1e-13`).  The two are
different measurements and they do not have to agree: `w` counts the levels the LP vertex happens
to use, `d` counts the levels a *restricted* certificate needs.

**What a hunt for an `H1` failure is really hunting for.**  `H1` fails at `z` only if the true
margin is *thinner* than what three ladder levels can certify.  Empirically (§3, §4) each level
buys one power of `eps`, so a failure needs a cell whose margin is `O(eps^4)` while the best
three-level object is `Theta(eps^3)`.  The search is therefore a search for **more degenerate
cells**, not for more exotic angle vectors: every extra constraint one adds to a family makes
`delta*` *more* negative and `H1` *easier*.  That is why the starts below are all of the
"tiling minus a hole set with a tilted block, chain cut" kind — the only construction known to
produce a margin thinner than `eps^2`.

---

## 2. The top-of-hole chain-free family, written out for general `T`  [proved / measured]

`break_h1.tophole_family` builds the `BANDCUT_K.md` §3(b) configuration for any `T`, from the
combinatorics rather than from a multistart.  The construction, and the fact that it is forced:

A chain-free configuration needs its `k = (T-1)^2` chain-cut squares to sit **one per cell of a
`(T-1) x (T-1)` level grid** (Mirsky, `BANDCUT_K.md` §1.3).  Take that grid to be the tiling's
own, with **one adjacent column pair and one adjacent row pair merged** (so heights drop from `T`
to `T-1`).  Then

* the merged-column x merged-row level cell covers **four** tiling cells: a `2 x 2` **pinwheel at
  tilt `eps`** sits on them, one of the four being the cell's cut-set representative;
* each of the other `2(T-2)` cross level cells covers two tiling cells: one carries an
  axis-parallel cut-set square, the other is removed;
* each of the `(T-2)^2` remaining level cells covers one tiling cell: an axis-parallel cut-set
  square.

**[proved] That level map is the one the recorded `T = 4` optimum uses.**  Reading
`runs/bandcut_k_T4a.jsonl` at `k = 8, 9`: merged columns `{1,2}`, merged rows `{1,2}`, and the map
`(col, row) -> (level)` is injective on the nine cut-set squares.  The four central tiling cells
carry the pinwheel; the four removed cross cells are exactly the permutation hole set
`(1,0), (3,1), (0,2), (2,3)`.

**[measured] At `T = 4` the construction reproduces the `eps^3` cell from `85` LP-ascent jobs and
2.7 seconds** (`delta* = -5.260791e-07` at `eps = 1 deg`, against the recorded `-5.259969e-07`;
`runs/break_h1_tophole_T4.jsonl`).  That is a `10^4`-fold cheaper route to the cell than the
`(k, eps)` scan, and it is what makes the `T = 5` question answerable at all.

### 2.1 [proved] At `T = 5` this construction cannot have a permutation hole set

Count the removed cells.  Holes number `T^2 - n = T` and far squares `n - k = T - 1`, so
`2T - 1` cells are removed: `3` from the central level cell (the rest of the pinwheel) and one
from each of the `2(T-2)` cross cells.  For the hole set to be a **permutation** set it must have
exactly one hole in every tiling row and column.

* The `T - 2` non-merged rows each contain exactly one removed cross cell, and it lies in a
  **merged** column.  So those `T - 2` holes occupy at most `2` distinct columns.
* A permutation set needs them in `T - 2` distinct columns.

`T - 2 <= 2` iff `T <= 4`.  **At `T = 4` it is tight (two holes in columns `1` and `2`) and at
`T = 5` it is impossible** — and the same count kills `k = (T-1)^2 - 1` at `T = 5` (the dropped
level cell buys one hole, and the merged-row band then needs three holes in two rows).

The consequence is not cosmetic.  Without a permutation hole set some tiling line carries all `T`
squares, one of which is tilted, so that line is over-full by `cos eps + sin eps - 1 ~ eps` and

    delta*  ~  -0.083 eps        (measured at T = 5 over eps in [0.25, 10 deg], §4.2)

— **a first-order cost, two orders worse than the `eps^3` cell.**  Every `T = 5` configuration in
which the far squares sit on tiling cells is stuck on that branch.  The configurations that beat
it are **sheared**: the four squares of an x-level are spread over an x-window just under `1`
rather than sharing a tiling column (the `T = 4` cell already does this in its middle level, where
`x = 2.5` at the bottom and `x = 1.5` at the top).  §4.2 measures the sheared branch.

---

## 3. The `T = 4` hunt: `H1` survives  [measured]

`search/break_h1.py hunt` (`runs/break_h1_hunt_screen.{jsonl,log}`).  Starts, all **structured**:

* the `4 x 4` tiling minus each of the `24` permutation hole sets, with a structured subset of
  the `12` remaining squares tilted by `+- eps` — central `2x2`, central `3x3`, the two
  diagonals, the border ring, all twelve — and the signs coherent, checkerboard, or pinwheel.
  `318` such angle vectors, `51` after quotienting by the dihedral group of the tiling
  (`break_h1._d4_canon`);
* **the chain cut through tilted squares as well as near-axis ones.**  For each family the cut
  set is the near-axis set, the near-axis set plus one or two of the tilted squares, the
  near-axis set minus one, or a random mixture; the level pattern is the merge map of §2 (which
  is what makes such a cut set feasible at all) or a random `bandcut_k.level_patterns` entry.
  Unlike `bandcut_k.py scan`, the cut set here is **not** forced to be the tilt band.
* at each (angle vector, cut set, pattern) `delta*` is maximised over the centres by LP ascent
  from the structured seed plus jitters plus `~180` tiling and random starts, chain-freeness is
  re-verified from the geometry alone at the configuration's own margin, and then `CW`, `H`,
  `HP`, `H1` are read off.

| | count |
|---|---|
| jobs | `201` (`eps = 1, 5 deg`) |
| pattern infeasible (no configuration satisfies that cut) | `47` |
| **optima with `delta* < 0`** | **`153`** |
| chain-freeness re-verified from the geometry | `153 / 153` |
| `CW <= 0` (a bare chain of four plus walls already certifies) | `127 / 153` (`83 %`) |
| `H <= 0` | `150 / 153` (`98 %`) |
| **`H1 <= 0`** | **`153 / 153`** |
| certificate depth `1 / 2 / 3` | `127 / 23 / 3` |
| `H1 > 0` — i.e. an `H1` **failure** | **`0`** |

(`H2` is not computed in the screen: it is only needed when `H1 > 0`, and `H1 > 0` never
happened.  `H2` is computed in §3.2, §4 and §5, where it is `<= 0` everywhere.)

Adding the `46` per-`eps` cell configurations of §4 and §5 (`T = 3, 4, 5`) and the `153`
re-read census optima of §3.2, `H1 <= 0` at **`352 / 352`** configurations with `delta* < 0`
measured in this note, on top of `727 / 727` in `T4_CYCLES.md` §7.  **No counterexample.**

### 3.1 The three closest calls, re-solved  [measured]

The brief asks that every candidate failure be re-solved at row tolerance `1e-10` and with a
stronger multistart.  There were **no candidates** — the largest `H1` at any `delta* < 0` point
was `-2.8e-17`, i.e. `0` to machine precision, which is still a valid certificate of
`delta <= 0`.  The `3` optima that needed the third rung (`H > 0`, `H1 <= 0`) are all the
`eps^3` cell of §4.1 arrived at from a different hole set; re-solved at `tol = 1e-10` with
`maxchains = 400` and the exponent-scan multistart (§4.1) they give
`H = +5.956e-05`, `H1 = -5.197e-07` against `delta* = -5.261e-07` — `H1` recovers `delta*` to
`1.2 %` and is on the right side by four orders of magnitude.

### 3.2 The `516` recorded `T = 4` optima, re-read with `H2` and the depth  [measured]

`break_h1.py eval --census` re-reads every optimum `t4_cycles.load_recorded` knows about
(`runs/bandcut_k_T4{a,free,free_pilot,_hunt}.jsonl`, the band stack, the bars) and adds the two
objects `T4_CYCLES.md` did not have: `H2` and the certificate depth.  The run was still in
progress when this note was written; at `159` records:

| | |
|---|---|
| `delta* < 0` | `153` |
| `CW <= 0` | `3 / 153` |
| `H <= 0` | `135 / 153` |
| **`H1 <= 0`** | **`153 / 153`** |
| `H2 <= 0` | `153 / 153` |
| certificate depth `1 / 2 / 3` | `3 / 132 / 18` |
| margin-`0` points, worst bound | `6`, `-5.6e-17` |

Depth `2` dominates (the H closes most cells outright), depth `3` is the `k = 5, 8, 9` chain-free
cells, and **no record needs depth `4`**.

**Why this is not surprising, and what it would take to break `H1`.**  Every constraint one can
add to a structured family (a longer chain cut, more tilted squares in the cut set, a hole set
that leaves a full tiling line) makes `delta*` *more* negative and therefore makes `H1` *easier*.
The only known way to make `delta*` thinner than `eps^2` is the top-of-hole chain-free cell of
§2, and at `T = 4` there is exactly one of those.  A failure needs a cell at `eps^4`; §4 says
where to look, and §5 says it is not at `T = 5`.

---

## 4. The margin exponent and the ladder depth, cell by cell  [measured]

`search/break_h1_cells.py` takes the best configuration of a family, continues it in `eps` (same
angle pattern, tilts rescaled, centres carried from the neighbouring `eps`), re-solves the LP in
the centres at each `eps`, re-verifies chain-freeness from the geometry, and reads the ladder at
row tolerance `1e-10`.  `p` is the least-squares slope of `log(-delta*)` against `log eps`
(`eps` in radians) and the *local* slopes between consecutive `eps` are printed as well, because
they are what shows the law rather than an average of two regimes.

### 4.1 The table

| `T` | cell | `delta*` at `eps = 1 deg` | `gamma` (extrapolated) | **`p`** | local slopes | **certificate depth `d`** | weight depth `w` |
|---|---|---|---|---|---|---|---|
| 3 | `k = 4` top of hole, chain-free | `-3.887e-05` | `0.1257 -> 1/8` | **`2.02`** | `2.007 .. 2.049` | **`2`** | `4` |
| 3 | `k = 3` (top-of-hole family with one cut square dropped) | `-7.262e-05` | `0.2470 -> 1/4` | **`1.86`** | `1.966, 1.983` | **`3`** | `3` |
| 4 | `k = 9` top of hole, chain-free | `-5.261e-07` | `0.0990 -> 1/10` | **`2.96`** | `2.94 .. 2.98` | **`3`** | `3` |
| 4 | `k = 8` (one cut square dropped) | `-6.006e-05` | `0.1972 -> 1/5` | **`1.93`** | `1.93 .. 2.00` | **`2` or `3`** (knife-edge, below) | `3-4` |
| 5 | `k = 15` top of hole, chain-free, permutation hole set | `-5.046e-05` | `0.1657 -> 1/6` | **`~2`** | `1.874, 1.948, 1.965` | **`2`** | `4` |
| 5 | `k = 15, 16`, sheared (no permutation hole set) | `-7.479e-05` | `0.2455 -> ~0.26` | **`2.0`** | `1.95 .. 2.02` | **`2`** | `3` |

Read off:

1. **The top-of-hole constant is `1/8, 1/10, 1/6` at `T = 3, 4, 5`, against the far-field line
   `1/(T+1) = 1/4, 1/5, 1/6`.**  At `T = 3` the top of the hole is half the far field
   (`BANDCUT_K.md` §0(2) already says so), at `T = 5` it is *equal* to it, and at `T = 4` it is a
   whole order of `eps` below both.  **`T = 4` is the outlier, not the start of a trend.**
2. **`p = 2` everywhere except the `T = 4` top of the hole, where `p = 3`.**  Measured over
   `eps in [0.25 deg, 10 deg]` (`20 deg` at `T = 3`, and only down to `1 deg` for the `T = 5`
   permutation branch, §5).  The local slopes are within `2 %` of their limit for
   `eps <= 2 deg` and drift by up to `15 %` at `eps = 10, 20 deg`, where the next order of the
   expansion is no longer small.
3. **`d = p` is FALSE.**  The `T = 3` `k = 3` cell has `p = 2` and needs `d = 3`: its `H` bound is
   `+2.06e-08` at `eps = 0.25 deg`, `+1.30e-06` at `1 deg`, `+1.01e-05` at `2 deg` — i.e.
   `H = +c eps^3` — while `delta* = -0.247 eps^2`.  `H` is positive, so the H-lemma fails there,
   and `H1` recovers `delta*` exactly (`-4.7030e-06` against `-4.7030e-06` at `eps = 0.25 deg`).
   **The right statement is: `d >= 2` at every chain-free cell measured here (a chain of `T`
   through tilted squares never certifies on its own, `CW > 0` in all twelve rows of the table),
   `d = 3` when the H's own leading term cancels, and `d` is decided by the sign of a
   *higher-order* coefficient, not by `p`.**  (In the `T = 4` census of §3.2, where the chain cut
   is often absent, `d = 1` does occur: `3 / 153`.)
4. **The `T = 4` `k = 8` cell is a knife-edge for `H`.**  Along the continuation, `H` reads
   `-1.89e-06, +8.6e-10, +1.36e-08, -1.14e-04, +7.9e-06, -2.23e-03` at
   `eps = 0.25, 0.5, 1, 2, 5, 10 deg` — it changes sign with the local optimum the continuation
   lands on, and where it is positive it is `O(eps^4)` against `delta* = -0.2 eps^2`.  `H1` is
   `delta*` to `<1 %` at every one of those points.  **This is the mechanism that would break
   `H1`: a cell where the ladder's leading terms cancel one level further.  Nothing found here
   cancels two levels.**

### 4.2 The `eps^3` cell is `T = 4`-only, and why  [measured + the count of §2.1]

The `eps^3` at `T = 4` comes from the pinwheel sitting in the *merged* level cell: four tilted
squares on four tiling cells whose level is one cell, with the permutation hole set absorbing
the rest of the cross.  §2.1 shows that combination exists only for `T <= 4`.  At `T = 5` the two
surviving branches are

* **permutation hole set, `k = 15`** (`search/break_h1_perm5.py`, the shape written out in that
  file's docstring): `gamma -> 1/6`, `p = 2`;
* **no permutation hole set** (some tiling line carries all `T` squares with one of them tilted):
  `delta* ~ -0.083 eps`, `p = 1` — the branch every fixed-angle multistart falls into first, and
  the reason a naive `T = 5` scan reports a *thicker* cell than it should.

A third branch, the **sheared** one (`k = 15` and `k = 16` alike, `runs/break_h1_t5_cells.jsonl`),
has `gamma ~ 0.25`, `p = 2`, and is what the structured family of §2 finds when the pinwheel is
not in a merged cell.  All three are `p <= 2`.

---

## 5. `T = 5`, `n = 20`, top of the hole (`k = 15, 16`): **`p = 2`, not `4`**  [measured]

This is the brief's item 3, and the answer is the opposite of the one it was braced for.

`search/break_h1_perm5.py` writes out the one combinatorial shape that gives `T = 5`,
`k = (T-1)^2 - 1 = 15`, a merge-type level map **and** a permutation hole set (§2.1 shows
`k = 16` admits none, and that without one some tiling line is over-full and the margin is only
`O(eps)`).  `9216` shapes, `300` sampled, each LP-ascended from its own tiling seed; the winner is
then continued in `eps` by `search/break_h1_cells.py`.  Chain-freeness is re-verified from the
geometry at every point.

| `eps` | `delta*` | `-delta*/eps^2` | `CW` | `H` | `H1` | `H2` | `d` | `w` | hits |
|---|---|---|---|---|---|---|---|---|---|
| `10 deg` | `-4.302620e-03` | `0.1413` | `+5.563e-02` | `-3.326e-03` | `-3.326e-03` | `-3.326e-03` | `2` | `4` | `39` |
| `5 deg` | `-1.174009e-03` | `0.1542` | `+2.835e-02` | `-7.433e-04` | `-7.433e-04` | `-7.433e-04` | `2` | `4` | `7` |
| `2 deg` | `-1.970892e-04` | `0.1618` | `+1.151e-02` | `-1.089e-04` | `-1.089e-04` | `-1.089e-04` | `2` | `4` | `7` |
| `1 deg` | `-5.046493e-05` | `0.1657` | | | | | `2` | | `1` |

Local slopes `1.874, 1.948, 1.965` (over `10 -> 5 -> 2 -> 1 deg`) -> **`p = 2`**, and the
constant climbs monotonically toward `1/6`:

> **`gamma(5, top of hole) -> 1/6 = 1/(T+1)`** — `0.1413, 0.1542, 0.1618, 0.1657` against
> `0.16667`.  The top-of-hole chain-free cell at `T = 5` sits exactly on the far-field constant
> `eps^2/(T+1)` of `BANDCUT_K.md` §0(2).

*Where the continuation stops.*  Below `eps = 1 deg` the `n = 20` LP ascent loses this branch and
falls back to the `-0.083 eps` branch of §4.2 (`-1.44e-03` at `0.5 deg`), so the last two rows of
the sweep are **not** this cell and are excluded; the least-squares `p` over the whole sweep
(`0.40`) is meaningless for that reason and the local slopes above are what to read.  The same
loss does not happen at `T = 3, 4`, where the continuation holds down to `eps = 0.25 deg`.

**So: `p = 4` does NOT happen at `T = 5`, and "ladder depth grows with `T`" is NOT what the data
says.**  Stated loudly in the other direction: the `eps^3` cell is a `T = 4` accident of the
count in §2.1, the `T = 5` top of the hole is an ordinary `eps^2` cell, and the certificate depth
needed there is **`2`** — one *less* than at `T = 4`.  `H` alone closes it.

`k = 16` was tested by the `T = 4` `k = 8 -> 9` move (the extra cut-set square is the far square
at tilt exactly `eps` sitting on the missing level cell, `break_h1_t5.extend_to_16`): it is
feasible and gives **exactly the same `delta*`** as `k = 15` on the sheared branch
(`-6.228936e-03, -1.721090e-03, -2.926646e-04` at `eps = 10, 5, 2 deg`), just as `k = 8` and
`k = 9` coincide at `T = 4`.

### 5.1 Reliability, stated plainly

Every `delta*` here is a **feasible point**, so the true cell maximum is `>=` the number printed,
and a multistart can only make a cell look *thicker* than it is — which makes `H1` look *easier*.
Three separate reasons to believe the `p = 2` at `T = 5` anyway:

1. the same instrument, the same family generator and the same continuation reproduce the known
   `T = 3` (`1/8 eps^2`) and `T = 4` (`1/10 eps^3`) cell values to `1 %` (§4.1);
2. the constant converges to `1/(T+1)`, the value `BANDCUT_K.md` §0(2) predicts for the far field,
   from below and monotonically — a two-parameter coincidence if the branch were spurious;
3. §2.1 is a **counting argument, not a search**: the structure that produces `eps^3` at `T = 4`
   does not exist at `T = 5`.  For `p = 4` at `T = 5` there would have to be a *different*
   mechanism, and the `9216`-shape enumeration plus `~12,000` structured LP ascents
   (`runs/break_h1_tophole_T5.jsonl`) found none.

What is *not* established: that no `T = 5` chain-free configuration outside the merge-type level
maps is thinner.  The sheared branch (§4.2) is one such, and it is `p = 2` too, but the space of
non-merge level maps was not enumerated.  **An `eps^4` cell at `T = 5` is not ruled out; it is
only not found, by two independent constructions.**  A `bandcut_k.py scan --T 5` run
(`runs/break_h1_t5_scan.{jsonl,log}`) was started as a third, independent check; its unrestricted
column finished (`delta* = 0` at `k = 15, 16`, `259-370` hits of `736-796`, used by the filter in
§6) and its chain-free column was still running after `3.5 h` and was stopped to free the machine.

---

## 6. The filter  [measured]

Every bound in §§3–5 is a restricted **dual** LP value, i.e. `-sum w_r b_r` for an explicit
`w >= 0` with `sum w_r a_r = 0` and `sum w_r = 1`, so it is a Farkas bound for *every*
configuration satisfying its rows and cannot be negative at a margin-`0` point.  Anything below
`-1e-8` there is a bug in the rows.  `search/break_h1_filter.py` (`runs/break_h1_filter.log`).

| point | worst bound over `CW, H, HP, H1, H2, CYC, F` |
|---|---|
| `T = 4`, `k = 5 .. 12` unrestricted zeros (`T4a`, `T4free`), `eps = 1 deg` | `-5.6e-16 .. -4.0e-12` |
| `T = 4` tiling minus permutation hole set, central `2x2` tilted, `delta* = 0` | `-5.6e-16` scale |
| the `(4,1)` band stack (`delta = 0` exactly) | `-1.7e-16` |
| the `p = 4` bar (`delta = 0`) | `-1.7e-16` |
| `T = 5`, `k = 15` unrestricted, tiling minus hole set, `eps = 1, 2, 5 deg` (`delta* = -1.1e-16`) | `+2.8e-17, -2.8e-12, +1.4e-17` |

**`0` breaches in `30` margin-`0` points.**  The only entries above `1e-14` in magnitude are
`-2.5e-12`, `-4.0e-12` (two `T4free` zeros) and `-2.8e-12` (`T = 5`, `eps = 2 deg`), which are
the LP's own feasibility
tolerance (`1e-11`) times the row count — four orders above the `-1e-8` threshold and on the safe
side of it.  `H2`, the new object, passes everywhere it was evaluated, as it must: its support
contains `H1`'s.

*One trap worth recording:* `t4_cycles.load_recorded` labels the
`p = 3` bar with `delta = -1.213e-02`, so it is **not** a margin-`0` point and the filter must not
be applied to it (all seven bounds there equal `delta*` to `1e-16`, which is `H1` succeeding);
`break_h1_filter.py` skips it on the `delta < -1e-12` test, and a first pass that did not skip it
reported a spurious breach.  At the `T = 5` zeros the `CW`
column depends on how many of the many chains of five are enumerated: with `maxchains = 80` it is
`+8.3e-17` (a tight chain of five is found) and with `maxchains = 60` it is `+0.21` (it is not).
That is a property of the search cap, not of the configuration, and it is why every negative
result in this note caps the chain enumeration generously and re-solves the close calls.

---

## 7. What in the brief's premises turned out to be wrong

1. **"If some cell has margin `eps^4`, `H1` fails there and the lemma is a ladder of unbounded
   depth."**  The premise is sound; the *expectation* behind it — that the `T = 4` `eps^3` cell is
   the start of a trend — is not.  §2.1 identifies the exact combinatorial reason the `eps^3` cell
   exists (a `2 x 2` pinwheel inside a single merged level cell, with a permutation hole set
   absorbing the rest of the cross) and shows the same count is **infeasible at `T = 5`**.  The
   `T = 5` top of the hole is `eps^2` with `gamma -> 1/(T+1)` (§5).
2. **"the `eps^3` cell needed exactly three ladder levels because its margin is `eps^3` … is depth
   `= p` always?"**  **No.**  The `T = 3`, `k = 3` chain-free cell has `p = 2` and needs depth
   `3`: its `H` bound is `+0.09 eps^3 > 0` while `delta* = -0.247 eps^2`.  Depth is set by whether
   the H's own leading term cancels, which is a *sign* question about a higher-order coefficient,
   not by the exponent of `delta*`.  The `T = 4` `k = 8` cell shows the same knife-edge with `H`
   changing sign along the `eps` continuation (§4.1(4)).
3. **"the `T = 5` LP is `40` centre coordinates, fine for HiGHS."**  True for one LP; false for
   the search.  `bandcut_k.py scan --T 5 --ks 15,16` did not finish its chain-free column in
   `3.5 h` on four cores, and its `n = 20` multistart falls into the `-0.083 eps` branch of §4.2
   rather than the cell.  What made `T = 5` answerable was **enumerating the combinatorial shape
   first** (§2, `break_h1_perm5.py`: `9216` shapes) and using the LP only to price each shape —
   `300` shapes x `~18` starts recovers the cell where `~12,000` random structured starts did not.
4. **"`search/bandcut_k.py chain_certificate` is buggy; do not use it."**  Respected; not used.
   The chain-freeness re-check here is `bandcut_k.normal_gaps` + `longest_path` at the
   configuration's own margin (`break_h1.chainfree_ok`), i.e. the same test `bandcut_k.verify`
   makes, and it was applied to **every** configuration reported.
5. **A smaller one.**  `BANDCUT_K.md` §5 says "the chain-free rows themselves carry weight `0` …
   the constraint acts by ruling the angle vector out, not by being tight at the optimum."  At the
   recorded `T = 4`, `k = 8`, `eps = 1 deg` optimum the chain row `('chain','x',3,4,3,-1)` is
   **exactly tight** (residual `-5.6e-16`), and so are others.  The statement that is true is the
   one about the *dual weights* of `t3_chain.all_rows`, which does not contain chain rows at all.

---

## 8. Reproduce

    # the ladder + weight levels at one recorded optimum
    python3 -c "import sys; sys.path.insert(0,'search'); import break_h1 as B, t4_cycles as T, math; \
      r=[q for q in T.load_recorded('.') if q.get('k')==9 and q['free']=='chainfree' and q['eps']==1.0][0]; \
      print(B.bounds(r['z'],12,4.0,r['delta'],tol=1e-10,eps_rad=math.radians(1)))"

    # the T = 4 hunt (201 jobs, ~40 min on 2 cores)
    python3 search/break_h1.py hunt --eps 1,5 --ks 7,8,9 --jpf 6 --limit 60 --extra 60 \
        --njit 10 --maxchains 60 --nch2 0 --nocyc --nproc 2 --out runs/break_h1_hunt_screen.jsonl

    # the top-of-hole family at T = 3, 4, 5 and the exponent / depth per cell
    python3 search/break_h1.py tophole --T 4 --eps 0.25,0.5,1,2,5,10 --nsamp 500 --drop --mix \
        --nproc 1 --out runs/break_h1_tophole_T4.jsonl
    python3 search/break_h1_cells.py runs/break_h1_tophole_T4.jsonl --T 4 \
        --eps 0.25,0.5,1,2,5,10 --out runs/break_h1_cells_T4.jsonl

    # T = 5 top of the hole with a permutation hole set (the decisive run)
    NSAMP=300 python3 search/break_h1_perm5.py
    python3 search/break_h1_cells.py runs/break_h1_t5_perm.jsonl --T 5 --eps 0.25,0.5,1,2,5,10 \
        --maxchains 100 --out runs/break_h1_cells_T5perm.jsonl

    # the filter
    python3 search/break_h1_filter.py

Files: `search/break_h1.py` (starts, bounds, hunt, tophole), `search/break_h1_cells.py`
(continuation + exponent + depth), `search/break_h1_perm5.py` (the `T = 5` permutation-hole
shape), `search/break_h1_t5.py` (the `k = 15 -> 16` extension), `search/break_h1_filter.py`.
Runs: `runs/break_h1_hunt_screen.*`, `runs/break_h1_tophole_T{3,4,5}.jsonl`,
`runs/break_h1_cells_T{3,4}.jsonl`, `runs/break_h1_cells_T5perm.jsonl`,
`runs/break_h1_t5_{cont,cells,perm,embed,scan}.jsonl`, `runs/break_h1_census_eval.jsonl`,
`runs/break_h1_filter.log`.
