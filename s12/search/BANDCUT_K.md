# bandcut-k: the margin with `k` near-axis squares, with and without a chain of `T`

*2026-09-20/21 (night).  Task `tasks/bandcut-k`.  Instruments `search/bandcut_k.py` (the `(k, eps)`
scan, the chain detector, the dual reader) and `search/bandcut_khunt.py` (stage 2: push the
zero-margin configurations into the chain-free problem).  Runs `runs/bandcut_k_T3*`,
`runs/bandcut_k_T4*`.  This is a **measurement**: every number below is the value of an explicit
configuration, hence a LOWER bound on the cell's true maximum.  A lower bound can refute a trend
(an outlier above it) but never establish one; the reliability of each cell is reported.*

## 0. Verdict, up front

1. **No chain-free configuration of margin `>= 0` was found**, at `T = 3` or at `T = 4`, at any
   `eps` tested.  Every `delta* = 0` configuration found in the hole `T <= k <= (T-1)^2` carries a
   wall-to-wall chain of `T` **among its near-axis squares**, and that chain certifies
   `delta <= 0` through the exact identity of `ARCH_TLEDGER.md` §2.1.  (Section 4.)
2. **The far-field constant is `eps^2/(T+1)`, and at `T = 3` the hole's chain-free constant is
   exactly half of it** — but **at `T = 4` the top of the hole degenerates to third order.**
   Writing the chain-free maximum as `c(k, eps) = -gamma eps^p`, `eps` in radians, the table gives
   `gamma eps^p` extrapolated to `eps = 0`:

   | `T` | `k = 0` | `1 <= k <= T-1` (far field; no chain of `T` is possible) | hole, `k` below the top two | hole, `k = (T-1)^2 - 1, (T-1)^2` |
   |---|---|---|---|---|
   | 3 | `0.375 eps^2` | `0.2491 eps^2` (`~ 1/4`) | `0.1250 eps^2` (`k = 3`) | `0.1249 eps^2` (`k = 4`) |
   | 4 | `0.387 eps^2` (`1°` only) | `0.1999 eps^2` (`~ 1/5`) | `0.100-0.200 eps^2` (`k = 4..7`) | **`0.1004 eps^3`** (`k = 8, 9`) |

   `1/4` and `1/5` are `1/(T+1)`, the number of `delta`-carrying rows on a wall-to-wall chain of
   `T` (`T-1` links plus two walls).  **"No chain `=>` `delta <= -c eps^2`" is therefore FALSE at
   `T = 4`**: at `k = 8, 9`, `eps = 1°`, there is a chain-free configuration of margin `-5.26e-7`.
   The configuration is the `4x4` tiling minus a *permutation hole set* with the central `2x2`
   replaced by a pinwheel at tilt `eps` — `ARCH_TLEDGER.md` §2.2's optimal family, `S6_LOCAL.md`
   §2's `j = 8` row — whose *unconstrained* member has `delta* = 0` and **does** carry a chain of
   four among its eight axis-parallel squares (verified directly, §3).  Cutting that chain is what
   costs only `0.1 eps^3`.
3. **What holds a chain-free optimum shut is a wall-to-wall chain of `T` that uses a TILTED
   square.**  The dual of the LP at the optimum is carried by two walls and `T-1` pair rows along
   one axis — the familiar chain certificate — except that one, two or all of its squares are far
   squares at tilt exactly `eps`.  `CENSUS.md` §2, `S6_SKELETON.md` §1.2 and `chains.py` all report
   that no `delta = 0` core ever needed a tilted square; that is not contradicted (every zero found
   here is chained through its near-axis squares), but at `delta < 0` in the hole the certificate
   changes character, and a proof of the hole has to price a chain with tilted members.  (§5.)
4. **The chain must be defined at the configuration's own margin, not at gap `>= 0`.**  With the
   naive definition the lemma is false for a trivial reason: in a zero-margin configuration whose
   chain link is at gap exactly `0`, overlapping that one link by `1e-6` deletes the chain and
   costs exactly `1e-6` of margin, so `sup {delta : no chain} = 0` is attained in the limit by
   configurations that are chains.  Measured and then designed out (§1.3).
5. **`k` in the hole is not uniformly harmless even unrestricted.**  At `T = 3`, `k = 3`, the
   *unrestricted* maximum is `0` for `eps <= 5°` but `-4.42e-4` at `10°` and `-2.27e-3` at `20°`:
   with exactly `T` near-axis squares the chain must be all of them, and the `n - T` far squares
   cannot then be accommodated at margin `0` once `eps` is large.  `k >= T+1` is `0` at every
   `eps` tested.

## 1. What is measured

### 1.1 The quantity

`delta*(theta)` is the value of the disjunctive LP in the centres of `S6_SKELETON.md` §3.1 —
**with the walls carrying `delta`**, `p_i + delta <= (c_i)_x <= T - p_i - delta`.  (`s6local.py`
uses hard containment instead, which agrees in sign but not in value once `delta < 0`, and this
whole file lives at `delta < 0`; `search/bandcut_k.py` has both and uses the wall-carrying one.)
A square is **near-axis** if its tilt — its distance to `0 mod 90°` — is `<= eps`, **far** if the
tilt is `>= eps`.  A cell of the table is a pair `(k, eps)` and the search is over *both* the
angles and the centres, `k` squares confined to `|theta| <= eps` and `n - k` to
`eps <= |theta| <= 45°` with a free sign (sign coherence matters, `S6_LOCAL.md` §5.1, so the
signs are sampled per start and then fixed by the optimiser's constraints).

### 1.2 What a chain is

For a pair `(i, j)` and each of the four edge normals, the separating gap is
`|n.(c_j - c_i)| - m_ij` with `m_ij = 1/2 + (|cos D| + |sin D|)/2`.  A normal is **x-type** if it
is nearer the global `x` axis than the `y` axis (unambiguous for tilt `< 45°`).  Put `i -> j` in
the x-DAG when some x-type normal separates the pair with `j` on the `+` side.  A **chain of `T`**
is a path on `T` vertices of the x-DAG (or of the y-DAG), **not necessarily a straight row**:
staircases count, as they must (`RANK8.md` §3.2's certificate at `T = 4` contains one).

Why a chain of `T` is the obstruction, exactly: the `T-1` link rows give
`n_s.(c_{s+1} - c_s) >= m_s + delta` and the two end walls `delta` each, so summing along the
chain and projecting on the global axis,

    (T + 1) delta  <=  sum_s sin(tau_s) Dy_s  -  sum_s (m_s - 1)  -  (P_1 + P_T - 1)
                       -  sum_s (1 - cos tau_s) Dx_s ,

`tau_s` the tilt of the owner of link `s`.  Every term on the right except the first is `<= 0`, so
**a chain of `T` axis-parallel squares forces `delta <= 0` outright**, and a chain of near-axis
squares forces it unless the chain *rises*.  For a common tilt `t` this is exactly
`ARCH_TLEDGER.md` §2.1: `delta <= B_T(R,t) = [(T-u) cos t + R sin t]/(T-1) - 1`, `u = cos t + sin t`,
which is `<= 0` iff the rise `R <= T tan(t/2) + cos t - sin t`.  That threshold is

| | `eps = 1°` | `2°` | `5°` | `10°` | `20°` |
|---|---|---|---|---|---|
| `R*(3, eps)` | `1.0086` | `1.0169` | `1.0400` | `1.0736` | `1.1267` |
| `R*(4, eps)` | `1.0173` | `1.0343` | `1.0837` | `1.1611` | `1.3030` |

so at these `eps` a chain that climbs **one** cell still kills, and one that climbs two does not
(`proof-architecture.md` §0a item 7).

### 1.3 How chain-freeness is enforced — and at which level

**The level.**  The DAG edge test is `gap >= delta`, the configuration's own margin, not
`gap >= 0`.  With `gap >= 0` the measurement is vacuous: the optimiser breaks a chain by
overlapping one link by `eta` and pays exactly `eta` (measured: at `T = 3`, `k = 3`, `eps = 5°`
the "best chain-free" configuration was a wall-to-wall column of three near-axis squares with the
lower link at gap `-1e-6` and margin `-1e-6`).  At the level of `delta` that configuration is a
chain again, and the numbers below stop being an artefact of `eta`.

**Completeness.**  For tilts `< 45°` every separating normal is x-type or y-type, so any two
disjoint near-axis squares are comparable in the x-DAG or in the y-DAG.  By Mirsky, "x-height
`<= T-1` and y-height `<= T-1`" is equivalent to the existence of level maps
`L, M : near -> {0..T-2}` strictly increasing along the respective edges, and comparability makes
`(L, M)` **injective**.  So a chain-free configuration is exactly one admitting an injective
**pattern** `P : near -> {0..T-2}^2`, and

    L(i) < L(j)  =>  forbid the x-edge j -> i;     L(i) = L(j)  =>  forbid both,

each "forbid" being the row `m_ij - s n_{o,0}.(c_j - c_i) + delta >= eta`, linear in the centres.
The chain-free maximum is therefore a maximum over the `C((T-1)^2, k)` patterns (deduplicated by
the dihedral group of the level grid: `1, 3, 8, 16, 23, 23, 16, 8, 3, 1` classes for
`k = 0..9` at `T = 4`) of a smooth program.  **Nothing is assumed about the far squares** — they
are free obstacles at any angle with tilt `>= eps`, which is what the hole demands.

**Verification.**  The pattern rows are never trusted: the value reported for a cell comes from a
configuration whose tilt band and whose chain-freeness are re-checked from the geometry alone, by
walking the DAG at the level of that configuration's own margin (`tol 1e-7`, `eta = 1e-6`).

### 1.4 The instrument

Per start: an LP in the centres at fixed angles (HiGHS, exact for a fixed separating-axis
assignment, iterated), then SLSQP over centres *and* angles with the assignment fixed, then the LP
again, twice.  Starts: the near squares on the `(T-1)^2` level grid at three scales, all `n`
squares on distinct cells of the `TxT` tiling, a wall-to-wall row/column of near-axis squares
(free runs), uniform random, and — stage 2 — every configuration any other cell produced,
relabelled so its `k` least-tilted squares are the near-axis set, plus Cleemann-type starts (a
diagonal band of far squares at `2 arctan(1/T) = 28.07°` cutting the container with axis-parallel
blocks on either side).  `eps` is swept down and up with continuation, which is what makes the
`n = 12` multistart usable at all (`S6_LOCAL.md` §1).

## 2. `T = 3`, `n = 6`.  The hole is `k in {3, 4}`

`runs/bandcut_k_T3.jsonl` + `runs/bandcut_k_T3_hardcontain.jsonl` (the second file is misnamed:
it is an earlier, SLSQP-only run at the **same**, wall-carrying definition, kept because it is an
independent search; the two agree to 10 digits wherever both have a value).
`eps in {1, 2, 5, 10, 20}°`, `~600-850` verified local optima per cell.

**(a) unrestricted** — `max delta*` with exactly `k` near-axis squares:

| `k` | `1°` | `2°` | `5°` | `10°` | `20°` |
|---|---|---|---|---|---|
| 0 | `-1.1079e-04` | `-4.2982e-04` | `-2.4500e-03` | `-8.3920e-03` | `-2.4370e-02` |
| 1 | `-7.2932e-05` | `-2.7989e-04` | `-1.5986e-03` | `-5.4982e-03` | `-1.6228e-02` |
| 2 | `-7.2330e-05` | `-2.7505e-04` | `-1.4835e-03` | `-4.8140e-03` | `-1.4971e-02` |
| **3** | `0` | `0` | `0` | `-4.4202e-04` | `-2.2700e-03` |
| **4** | `0` | `0` | `0` | `0` | `0` |
| 5 | `0` | `0` | `0` | `0` | `0` |
| 6 | `0` | `0` | `0` | `0` | `0` |

**(b) chain-free** — `c(k, eps)`, the same with no chain of `3` among the near-axis squares
(`k >= 5 > (T-1)^2 = 4` is infeasible: chain counting forces a chain, which is the theorem):

| `k` | `1°` | `2°` | `5°` | `10°` | `20°` |
|---|---|---|---|---|---|
| 1 | `-7.2932e-05` | `-2.7989e-04` | `-1.5986e-03` | `-5.4982e-03` | `-1.6228e-02` |
| 2 | `-7.2330e-05` | `-2.7505e-04` | `-1.4835e-03` | `-4.8140e-03` | `-1.4971e-02` |
| **3** | `-3.7573e-05` | `-1.4823e-04` | `-8.8660e-04` | `-3.2712e-03` | `-1.0897e-02` |
| **4** | `-3.6932e-05` | `-1.4329e-04` | `-8.1732e-04` | `-2.8052e-03` | `-8.1899e-03` |

(`k <= 2 = T-1`: no chain of `3` can exist among `k` squares, so (b) `=` (a) there, and the
agreement of the two independent searches is a check on the instrument.)

`gamma = -c/eps^2` and its linear-in-`eps` extrapolation to `eps = 0`:

| `k` | `1°` | `2°` | `5°` | `10°` | `20°` | `eps -> 0` |
|---|---|---|---|---|---|---|
| 0 | `0.3637` | `0.3528` | `0.3217` | `0.2755` | `0.2000` | **`0.3747`** (`~ 3/8`) |
| 1 | `0.2394` | `0.2297` | `0.2099` | `0.1805` | `0.1332` | **`0.2491`** (`~ 1/4`) |
| 2 | `0.2374` | `0.2257` | `0.1948` | `0.1580` | `0.1229` | **`0.2492`** |
| 3 (cf) | `0.1233` | `0.1217` | `0.1164` | `0.1074` | `0.0894` | **`0.1250`** (`= 1/8`) |
| 4 (cf) | `0.1212` | `0.1176` | `0.1073` | `0.0921` | `0.0672` | **`0.1249`** |

The hole's two cells give the same constant to four digits, and it is exactly half the far field's.
The optimal pattern is the whole level grid in both cases (`k = 3`: the `L`-tromino of the `2x2`
grid; `k = 4`: all four cells), i.e. **the near-axis squares spread one per level cell** — they are
not band-like, they are grid-like, which is the `T = 3` shadow of the Cleemann picture.

Best chain-free configuration at `k = 4`, `eps = 5°` (`c = -8.1732e-04`; four near-axis squares at
the corners of a `2x2` level grid, two far squares at `-5°` wedged between them):

      sq 0  x=0.502460  y=0.499183  tilt= -0.0000   NEAR
      sq 2  x=2.500817  y=0.499183  tilt= +0.0000   NEAR
      sq 5  x=0.540858  y=1.540041  tilt= -5.0000
      sq 4  x=1.540129  y=1.459380  tilt= -5.0189
      sq 1  x=1.501635  y=2.500817  tilt= +0.0000   NEAR
      sq 3  x=2.500817  y=2.420028  tilt= +0.0000   NEAR

## 3. `T = 4`, `n = 12`.  The hole is `k in {4, ..., 9}`

`runs/bandcut_k_T4a.jsonl` (all cells, cheap), `runs/bandcut_k_T4free.jsonl` and
`runs/bandcut_k_T4free_pilot.jsonl` (the unrestricted column, deep), `runs/bandcut_k_T4_hunt.jsonl`
(stage 2).  `eps in {1, 5, 10}°`.  All `74` chain-free level patterns of `k = 4..9` were run; the
cell value is the best over them.

**(a) unrestricted.**  `delta* = 0` for **every** `k >= 4`, at every `eps`, with hit rates
`41-695` out of `240-696` verified local optima — solid.  For `k <= 3` (where no chain of `4` can
exist at all, so this row is also the far field and also the chain-free value):

| `k` | `1°` | `5°` | `10°` |
|---|---|---|---|
| 0 | `-1.1778e-04` | `-1.4191e-03` | `-5.2779e-03` |
| 1, 2, 3 | `-6.0077e-05` | `-1.4191e-03` | `-5.2779e-03` |

`gamma` for `k = 1, 2, 3` reads `0.1972, 0.1863, 0.1733` and extrapolates to **`0.1999`**, i.e.
`1/5 = 1/(T+1)` — the `T = 3` far-field constant was `0.2491 ~ 1/4 = 1/(T+1)` on the nose, so
**the far-field constant is `eps^2/(T+1)` at both `T`.**  (`k = 0` is under-explored: its `5°` and
`10°` entries have collapsed onto the `k <= 3` value, which means the search found a configuration
with a square at tilt exactly `eps`; the `1°` entry, `gamma = 0.387`, is the only honest one.)

**(b) chain-free, `c(k, eps)`:**

| `k` | `1°` | `5°` | `10°` | `-c/eps^2` at the three `eps` | law (`eps -> 0`) |
|---|---|---|---|---|---|
| 4 | `-5.9564e-05` | `-1.3638e-03` | `-4.9117e-03` | `0.1955, 0.1791, 0.1612` | `-0.200 eps^2` |
| 5 | `-3.0569e-05` | `-7.7469e-04` | `-3.1481e-03` | `0.1004, 0.1017, 0.1033` | `-0.100 eps^2` |
| 6 | `-4.5023e-05` | `-9.3278e-04` | `-3.9227e-03` | `0.1478, 0.1225, 0.1288` | `-0.154 eps^2` |
| 7 | `-4.0436e-05` | `-9.3278e-04` | `-2.8994e-03` | `0.1327, 0.1225, 0.0952` | `-0.135 eps^2` |
| **8** | `-5.2600e-07` | `-6.1889e-05` | `-4.5961e-04` | `0.0017, 0.0081, 0.0151` | **`-0.100 eps^3`** |
| **9** | `-5.2600e-07` | `-6.1889e-05` | `-4.5961e-04` | `0.0017, 0.0081, 0.0151` | **`-0.100 eps^3`** |

**The headline of this file is the last two rows.**  At `k = 8` and `k = 9` the chain-free maximum
is *cubic* in `eps`, not quadratic: `-c/eps^3` reads `0.0989, 0.0931, 0.0865` and extrapolates to
`0.1004`, while `-c/eps^2` reads `0.0017, 0.0081, 0.0151` and goes to `0`.  At `eps = 1°` the best
chain-free configuration has margin `-5.26e-7`; the whole of the hole's deficit has collapsed by
a factor of `eps`.  Any lemma of the form "no chain `=>` `delta <= -c eps^2`" is **false at
`T = 4`**, and `T = 3` gave no warning of it (there the top cell `k = (T-1)^2 = 4` is still
quadratic, `-1/8 eps^2`).

The configuration that does it is the repo's own deepest stratum with its chain cut.  At
`k = 9`, `eps = 10°` (`c = -4.5961e-04`, `43/142`):

      sq 0  x=0.647685  y=0.499851  tilt= +0.036   NEAR      sq 4  x=1.578409 y=1.432380 tilt=+10.000
      sq 3  x=2.500701  y=0.499851  tilt= +0.036   NEAR      sq11  x=2.567620 y=1.578409 tilt=+10.000
      sq 6  x=3.500149  y=0.647685  tilt= +0.036   NEAR      sq 9  x=1.432380 y=2.421591 tilt=+10.000
      sq 1  x=0.499851  y=1.499299  tilt= +0.036   NEAR      sq10  x=2.421591 y=2.567620 tilt=+10.000
      sq 7  x=3.500149  y=2.500701  tilt= +0.036   NEAR
      sq 2  x=0.499851  y=3.352315  tilt= +0.036   NEAR      near-axis heights x = 3, y = 3
      sq 5  x=1.499299  y=3.500149  tilt= +0.036   NEAR      all-squares heights x = 4, y = 4
      sq 8  x=3.352315  y=3.500149  tilt= +0.036   NEAR

(the four right-hand squares carry the full tilt `eps`; `sq 4` is the one the `k = 9` labelling
counts as near-axis at tilt exactly `eps`, which is why the `k = 8` and `k = 9` rows coincide).
The winning level pattern at `k = 8` is the `3x3` level grid **minus its centre cell**, and at
`k = 9` the whole grid: the near-axis squares ring the container and the centre is left to the
tilted ones, the same picture at every `eps`, with the `k = 4..7` winners sub-patterns of it.

The configuration is the `4x4` tiling **minus a permutation hole set** — holes at `(1.5,0.5)`, `(3.5,1.5)`,
`(0.5,2.5)`, `(2.5,3.5)` — with the central `2x2` replaced by a **pinwheel of four squares at tilt
exactly `eps`** and the eight border squares left (almost) axis-parallel.  That is
`ARCH_TLEDGER.md` §2.2's optimal family and `S6_LOCAL.md` §2's `j = 8` row, which is the
`theta = (0^8, t^4)` stratum of `Z_12`.  Checked directly here (`s6local` structured starts, LP in
the centres, soft walls):

    theta = (0^8, t^4),  t = 1°, 5°, 10°:   delta* = -1.1e-16, -3.3e-16, -4.4e-16   (i.e. 0)
    near-axis (the eight axis-parallel) heights:  x = 4 in all three cases

**so the zero-margin member of the family does contain a wall-to-wall chain of four among its
near-axis squares, and it is chain-freeness alone that costs the `0.1 eps^3`.**  The permutation
hole set leaves `x`-height `3` when the eight border squares are arranged the chain-free way, and
the price of that arrangement is third order.  No counterexample — but the margin is a factor
`eps` thinner than anywhere else in the table.

## 4. Are the `delta* = 0` configurations always chained?

`bandcut_k.py chains` walks, for the best configuration of every cell, **every** path on `T`
vertices of the near-axis DAG at that configuration's margin, and reports the chain with the
smallest rise together with the bound it certifies.  *Erratum (2026-09-21, `notes/t3-chain.md`): `chain_certificate` evaluates `B_T` at the chain's **mean** tilt without checking that the links share a normal, and converts the tilt to degrees twice; on near-axis chains (tilt `<= eps`) the `B_3`/`B_4` columns below are approximate to `O(eps)`, and at large tilt they are not bounds at all.  The chain counts and rises are unaffected.*  At `T = 3` (`runs/bandcut_k_T3*`):

| cell | `delta*` | chains of 3 among the near-axis squares | best rise | `B_3` (wall-to-wall) |
|---|---|---|---|---|
| `k=3`, `1°` | `-5.6e-17` | 1 | `+0.014` | `0` |
| `k=3`, `5°` | `-1.1e-16` | 1 | `0.000` | `0` |
| `k=3`, `10°` | `-4.4202e-04` | 1 | `-0.693` | `-6.96e-05` |
| `k=3`, `20°` | `-2.2700e-03` | 1 | `-0.001` | `-2.14e-04` |
| `k=4`, `1°` | `0` | 2 | `-0.000` | `0` |
| `k=4`, `10°` | `-1.1e-16` | 2 | `+0.369` | `-1.1e-16` |
| `k=5`, `5°` | `0` | 2 | `+1.915` | `0` |
| `k=6`, any | `0` | 8 | `~0` | `0` |
| `k<=2`, any | `<0` | 0 (impossible) | — | — |

At `T = 4` (`runs/bandcut_k_T4a`, `_T4free`, `_T4free_pilot`), every unrestricted cell
`k = 4..12`, every `eps`:

| cell | `delta*` | chains of 4 among the near-axis squares | best rise | `B_4` (wall-to-wall) |
|---|---|---|---|---|
| `k=4`, `1/5/10°` | `0` | 1 | `0.000 / 0.000 / -0.345` | `0 / 0 / 0` |
| `k=5`, `1/5/10°` | `0` | 2 / 1 / 1 | `+1.000 / -0.000 / +0.404` | `1.5e-10 / 0 / 0` |
| `k=6`, `1/5/10°` | `0` | 1 / 1 / 2 | `+0.000 / 0.000 / +1.419` | `0 / 0 / -1.1e-09` |
| `k=7`, `1/5/10°` | `0` | 3 / 3 / 1 | `-0.989 / -0.218 / +0.197` | `-1.1e-14 / 0 / 0` |
| `k=8`, `1/5/10°` | `0` | 6 / 4 / 3 | `+0.000 / -1.000 / +1.000` | `-4.1e-06 / 0 / -1.3e-10` |

**Every zero-margin configuration found, at `T = 3` and at `T = 4`, has a chain of `T` among its
near-axis squares, and the chain bound it certifies is `<= 0`.**  No zero-margin configuration was
found whose near-axis squares are chain-free, with or without the cross-feed of stage 2.

Two remarks on the rises.  (i) When the chain's squares are (numerically) axis-parallel the rise
is irrelevant — `sin t = 0` kills the `R sin t` term — which is why rows with rise `1.42` or `1.92`
still certify `0` even though `R*(4, 10°) = 1.161`.  The rise threshold only bites for a chain that
carries tilt.  (ii) The single sharpest test available is the `(0^8, t^4)` stratum of §3: its eight
axis-parallel squares have `x`-height exactly `4` at `t = 1°, 5°, 10°`, so the chain is there, and
forbidding it is what costs the `0.1 eps^3`.

## 5. What holds a chain-free optimum shut

`bandcut_k.py cert` reads the dual of the LP in the centres at the recorded optimum (the chain rows
included).  At `T = 3`, `k = 4`, `eps = 10°` (`c = -2.8052e-03`), support and weights:

    w=0.220051  ('lo','x',1)            near-axis
    w=0.223445  ('pair',1,4,...)        square 4 is FAR (10.000 deg)
    w=0.220051  ('pair',3,4,...)        square 4 is FAR
    w=0.220051  ('hi','x',3)            near-axis
    w=0.038801  ('lo','y',0) + ('pair',0,1) + ('hi','y',4)     a second, transverse chain
    total: wall 0.5177, pair 0.4823, chain rows 0.0000

That is a **wall-to-wall x-chain of three, `1 - 4 - 3`, whose middle square is tilted by `eps`** —
two walls and two pair rows at weight `~0.22` each, `4 x 0.22 + ...= 1`, the signature of
`S6_SKELETON.md` §1.2 with one square replaced by a tilted one.  At `k = 3`, `eps = 5°` the chain is `5 - 4 - 3`, **all three of them far squares at `-5°`**
(`lo-y(5)`, `pair(4,5)`, `pair(3,4)`, `hi-y(3)` at weights `0.2347, 0.2338, 0.2338, 0.2329`).  And
at `T = 4`, `k = 9`, `eps = 10°` (`c = -4.5961e-04`, the cubic cell):

    w=0.170565  ('lo','y',3)         w=0.173177  ('pair',3,4)   square 4 at tilt 10.000
    w=0.173177  ('pair',4,9)         w=0.175830  ('pair',5,9)   square 9 FAR at tilt 10.000
    w=0.175848  ('hi','y',5)
    w=0.029966  ('lo','x',1) + ('pair',1,9) + ('pair',3,6) + ('hi','x',6)   -- a transverse chain
    total: wall 0.4065, pair 0.5827, chain rows 0.0000

a wall-to-wall `y`-chain of **four**, `3 - 4 - 9 - 5`, of which the two middle squares carry the
full tilt `eps`, plus a transverse `x`-chain of four at a tenth of the weight.  The chain-free rows
themselves carry weight `0` in every case: the constraint acts by ruling the angle vector out, not
by being tight at the optimum.

**This is the answer to the second question in the brief.**  The censuses report that a zero-margin
*core* never needs a tilted square (`CENSUS.md` §2, `S6_SKELETON.md` §1.2); what is re-tested here
is the weaker, angle-space version of that — **every** zero-margin configuration found, in every
`(k, eps)` cell at both `T`, has a wall-to-wall chain of `T` among its *near-axis* squares (§4), so
nothing in this regime contradicts them.  What is new is that a **negative**-margin obstruction in
the hole *does* need a tilted square: once the near-axis squares are forbidden a chain of `T`, the
only chains of `T` left run through far squares, and the certificate is exactly such a chain.  (The
first-order core analysis of `chains.py` was not re-run on these `delta < 0` configurations — it is
built for `delta = 0` — so this is a statement about the LP dual, not about an irreducible core.)  A proof of the hole must therefore either (i) bound the chain
of `T` with `1` or `2` tilted members — for which `ARCH_TLEDGER.md` §2.1 already has the exact
line, the tilted links costing `|D|/2` each and the rise repaying `R sin t` — or (ii) never look at
chains through far squares at all, in which case it needs a mechanism this measurement does not see.

## 6. The lemma the data supports

> **Conjecture (measured; `T = 3` at `eps <= 20°`, `T = 4` at `eps <= 10°`).**  Let `n = T(T-1)`
> unit squares sit in `[0,T]^2` with margin `delta`, let `k` of them have tilt `<= eps` and
> `n - k` have tilt `>= eps`, `T <= k <= (T-1)^2`.  Build the x- and y-DAGs on the near-axis
> squares with the edge test `gap >= delta`.  Then:
>
> * if some DAG has a path on `T` vertices, `delta <= B_T(R, t) <= 0` whenever that path's rise
>   `R` is at most `T tan(t/2) + cos t - sin t`, `t` its common tilt (exact,
>   `ARCH_TLEDGER.md` §2.1; in particular always, when the path's squares are axis-parallel);
> * if neither DAG has one, `delta < 0` — strictly — but the size of the deficit **depends on `k`**:
>
>       delta  <=  -gamma_T(k) eps^2    for  T <= k <= (T-1)^2 - 2,
>       delta  <=  -gamma_T    eps^3    for  k = (T-1)^2 - 1  and  k = (T-1)^2,
>
>   with (extrapolated to `eps = 0`, the values accurate to `~1%`)
>
>   | | `T = 3` | `T = 4` |
>   |---|---|---|
>   | far field `k <= T-1` (no chain is possible) | `eps^2/4` | `eps^2/5` — i.e. `eps^2/(T+1)` |
>   | hole, chain-free, `k` below the top two | `eps^2/8` (`k = 3`) | `eps^2/10` to `eps^2/5` (`k = 4..7`) |
>   | hole, chain-free, top two `k` | `eps^2/8` (`k = 4`) | **`eps^3/10`** (`k = 8, 9`) |

Four riders the data attaches to it.

* **The exponent is not `2` everywhere.**  The clean form "no chain `=>` `delta <= -c eps^2`" is
  what `T = 3` suggests and what `T = 4` refutes, at exactly the cells `k = (T-1)^2 - 1, (T-1)^2`
  where chain counting is about to take over.  If a proof is going to use a quantitative version of
  the lemma it must survive `eps^3`, i.e. it must be as sharp as the third-order behaviour of the
  pinwheel-plus-permutation-hole family.  `eps = 1°` at `T = 4` means a margin of `5e-7`, which is
  smaller than the `3.76e-6` of `RANK8.md` §2.2 and inside the `1.2e-3` second-order remainder of
  `RANK8.md` §3.3: **no linearised argument can see it.**
* **The far-field constant is `eps^2/(T+1)` at both `T`,** `T+1` being the number of
  `delta`-carrying rows of a wall-to-wall chain (`T-1` links and two walls).  That is a cleaner and
  more useful number than anything measured for A3 so far, and it is what A5 would have to match.
* **`gamma` shrinks like `1/T`.**  Both the far field and the hole are `O(1/T)`; the hole does not
  become easier as `T` grows, and the top-of-hole cell degenerates by a further factor of `eps`
  between `T = 3` and `T = 4`.  If that degeneration continues (`eps^4` at `T = 5`?) then the
  lemma's *order*, not its constant, is what carries `T < T*`.  **Measuring `k = (T-1)^2` at
  `T = 5` is the single cheapest follow-up this file suggests.**
* **The conclusion the architecture needs is only `delta <= 0`,** which all of this implies with
  room to spare; the quantitative part matters only for A5's glue and for any numerical
  certificate, and it is there that `eps^3` hurts.

## 7. Reliability: which cells are under-explored

The asymmetry to keep in mind: a multistart can only **refute**.  Every `c(k, eps)` here is a
configuration that exists, so the true maximum is `>= ` the number printed; a cell reported at
`-0.15 eps^2` may really be at `-0.10 eps^2` or at `0`.  What the file can assert without
qualification is (i) the values as feasible points, (ii) that no chain-free configuration of
margin `>= 0` turned up, and (iii) the `eps`-scaling *within* a cell, which is read off three or
five points of the same family and is therefore much more robust than the constant.

| cell | starts / verified optima | hit rate at the reported value | verdict |
|---|---|---|---|
| `T = 3`, all `k`, all `eps` | `530-850` | `15-750` | **solid**; two independent runs agree to 10 digits |
| `T = 4`, `k >= 4`, unrestricted | `240-696` | `41-695` | **solid** (`delta* = 0` is easy to hit once `chainrow` starts exist) |
| `T = 4`, `k = 0`, `eps = 5, 10°` | `300` | `1, 9` | **broken**: the value has collapsed onto the `k <= 3` row via a square at tilt exactly `eps`.  Do not quote |
| `T = 4`, `k = 1, 2, 3` | `240-300` | `1-99` | fair; the `gamma` extrapolation is stable across the three `k`, which is the check |
| `T = 4`, `k = 4`, chain-free | `167-229` per pattern | `1-52` | **weakest cell of the table**: `k = 4` reads `-0.200 eps^2` against `-0.100` at `k = 5`, which is not credible as a shape.  Almost certainly an under-estimate; read `c(4,·)` as "`>= -0.200 eps^2`" only |
| `T = 4`, `k = 5, 6, 7`, chain-free | `80-223` per pattern | `3-90` | fair; the non-monotone `0.100, 0.154, 0.135` again suggests the true profile is flatter, probably `~0.100` throughout |
| `T = 4`, `k = 8, 9`, chain-free | `142-216` per pattern | `2-43` | **solid where it matters**: the cubic law is read off three `eps` and reproduces `0.0989, 0.0931, 0.0865`, and the family is identified (it is a known stratum), so the value is close to its maximum |

`eps` was pushed to `20°` only at `T = 3`; at `T = 4` nothing above `10°` was measured, and
`10°` is below the straight-chain window `2 arctan(1/4) = 28.07°`, so the whole `T = 4` table lives
in the regime where a level chain still kills.  The window itself is untested here.  `T = 5` is
untested entirely, and §6's last rider says that is the measurement to do next.

Two systematic caveats.

* **The band is closed on both sides** (`|theta| <= eps` for near, `>= eps` for far), so a
  configuration with a square at tilt exactly `eps` is feasible for cell `k` and for cell `k-1`.
  That is how `T = 4`, `k = 0` collapsed, and it is why the `k = 8` and `k = 9` rows are identical
  (they are one configuration with one square at `eps`, counted either way).  Under the strict
  reading "near-axis means tilt `< eps`" every such witness belongs to the smaller `k`, and
  chain-freeness is inherited by subsets, so the table is still a valid set of witnesses; only the
  labels shift by one.
* **`eta = 1e-6`** is the slack by which a forbidden separation must be beaten, and the
  verification tolerance is `1e-7`.  At `T = 4`, `k = 8`, `eps = 1°` the reported `-5.26e-7` is of
  the same order as `eta`; it is not an artefact (the same family reads `-6.19e-5` at `5°` and
  `-4.60e-4` at `10°`, both `>> eta`, and `0.1 eps^3` predicts `-5.3e-7` at `1°` from those), but
  that one entry should not be quoted to more than one digit.

## 8. Reproduce

    python3 search/bandcut_k.py scan --n 6 --T 3 --eps 1,2,5,10,20 --starts 200 --polish 1 \
        --patmink 3 --nproc 3 --out runs/bandcut_k_T3.jsonl                        # ~25 min
    python3 search/bandcut_k.py scan --n 12 --T 4 --eps 1,5,10 --starts 60 --pstarts 30 \
        --polish 0 --cycles 1 --rounds 6 --patmink 4 --nproc 8 --out runs/bandcut_k_T4a.jsonl
    python3 search/bandcut_k.py scan --n 12 --T 4 --eps 1,5,10 --ks 4,5,6,7,8,9,10,11,12 \
        --mode free --starts 150 --polish 2 --nproc 3 --out runs/bandcut_k_T4free.jsonl
    python3 search/bandcut_khunt.py --n 12 --T 4 --eps 1,5,10 --ks 4,5,6,7,8,9 \
        --seeds runs/bandcut_k_T4a.jsonl runs/bandcut_k_T4free.jsonl --maxseeds 40 \
        --nproc 11 --out runs/bandcut_k_T4_hunt.jsonl
    python3 search/bandcut_k.py report runs/bandcut_k_T3.jsonl runs/bandcut_k_T3_hardcontain.jsonl
    python3 search/bandcut_k.py report runs/bandcut_k_T4a.jsonl runs/bandcut_k_T4free.jsonl \
        runs/bandcut_k_T4free_pilot.jsonl runs/bandcut_k_T4_hunt.jsonl
    python3 search/bandcut_k.py chains --T 4 runs/bandcut_k_T4a.jsonl runs/bandcut_k_T4free.jsonl
    python3 search/bandcut_k.py cert   runs/bandcut_k_T4a.jsonl --T 4 --k 9 --eps 10 --wtol 0.01
    python3 search/bandcut_k.py show   runs/bandcut_k_T4a.jsonl --k 9 --eps 10

The `(0^8, t^4)` check of §3 (the one that decides whether the cubic cell is a counterexample):

    python3 -c "import sys,math,numpy as np; sys.path.insert(0,'search')
    import bandcut_k as B, s6local
    rng=np.random.default_rng(0)
    for td in (1.,5.,10.):
        th=[0.]*8+[math.radians(td)]*4; best=(-9,None)
        for (X,Y) in s6local.structured_starts(12,4.,[0]*8+[1]*4,rng,2500,3,0.08):
            z=np.zeros(37); z[1::3],z[2::3],z[3::3]=X,Y,th
            w,v=B.lp_polish(z,12,4.,8,None,1e-6)
            if v>best[0]: best=(v,w.copy())
        print(td, best[0], B.chain_report(best[1],12,range(8),best[0])['hx_near'])"
