# t3-chain: the H-lemma at `T = 3`, stated exactly, tested over the whole angle space, and what it cannot do

*2026-09-21.  Task `tasks/t3-chain/README.md`.  Code (new, nothing in `search/` modified):
`search/t3_chain.py` (the restricted-dual instrument and the angle-space scan),
`search/t3_chain_hform.py` (the closed forms and their verification),
`search/t3_chain_shape.py` (what shape the true dual has),
`search/t3_chain_exist.py` (does an H exist in the configuration at all),
`search/t3_chain_repair.py` (the minimal repair of a failing H),
`search/t3_chain_mt.py` (MT on chains with a genuinely common separating normal).
Runs: `runs/t3_chain_scan1.jsonl`, `runs/t3_chain_shape1.txt`, `runs/t3_chain_exist1.txt`,
`runs/t3_chain_hform.txt`, `runs/t3_chain_filter.txt`, `runs/t3_chain_repair1.txt`,
`runs/t3_chain_mt1.txt`, `runs/t3_chain_counterex.txt`, `runs/t3_chain_pinwheel.txt`.
Labels: **[proved]** = an identity checked by hand and re-checked in exact/float arithmetic by a
script here; **[measured]** = a number produced by a multistart, hence a feasible point and a
LOWER bound on the quantity it estimates; **[heuristic]**; **[guess]**.*

---

## 0. Verdict, up front

1. **The H-lemma can be stated exactly, and in closed form the whole of it is one integer.**
   At a common tilt `t`, a wall-to-wall chain of `T` carrying `L` transverse links has an explicit
   Farkas dual whose value is
   `H(T, L, t) = omega (T + 2 - u) - 1`, `u = cos t + sin t`,
   `omega = M / ((T-1) + L tan t + 2M)`, `M = (1 + sin^2 t + sin t cos t)/cos t`,
   and `H <= 0` **iff `L >= L*(T, t) = (M (T - u) - (T-1)) cot t`**, with `L*(T, 0) = T - 2`.
   **[proved]** (§1.3, `t3_chain_hform.py check`: the weights cancel every centre coordinate to
   `1e-12` for `T in {3,4,5,11}`, `a, c in {0,1,2}`, `t in {0.5, 5, 20, 45}` deg).

2. **MT is not a certificate, and making it one is where `T = 3` already hurts.**  The
   architecture's claim (`proof-architecture.md` §0a item 14) is *every configuration has a
   wall-to-wall chain of `T` on which MT gives `delta <= 0`*.  Tested correctly — enumerating
   chains **link normal by link normal**, so that MT's common-normal hypothesis actually holds —
   it is true at `602 / 720` sampled optima, and at the other `118` **no wall-to-wall chain of
   three with a common separating normal exists at all** (`search/t3_chain_mt.py`,
   `runs/t3_chain_mt1.txt`).  In the uniform families it is *exactly tight*: `MT = delta*` to
   `1e-9` at `unif1, 5, 10, 20, 30, 40, 44, 45`.  **But MT reads the rise `R` off the
   configuration**; it is a valid bound, not a Farkas certificate (it does not cancel the
   transverse coordinates, it substitutes their measured values).  Turning it into one is exactly
   the H, and **the H is strictly weaker**: at `theta = (45 deg)^6`,
   `MT = delta* = (93 - 66 sqrt2)/7 = -0.048299` while the best H gives
   `(10 - 7 sqrt2)/2 = +0.050253`.  **The price of eliminating the rise is `0.0986` of margin at
   `T = 3, 45 deg` — twice the margin itself.**  [proved] for both closed forms, [measured] for
   the rest.

3. **The H-lemma does not hold everywhere at `T = 3`.**  `H <= 0` at `652 / 720` sampled optima;
   it fails at `68`, systematically for uniform tilt `>= 28.8 deg` (`unif30/40/44/45`: `0/12`
   each) and sporadically in the coherent small-angle cone (`10 / 72`; `coh2` is `H <= 0` at only `7/12`).  The large-tilt
   failure is exact and has a one-line reason: `H(3, L, t) <= 0` needs `L >= L*(3, t)`, and
   `L*(3, t)` crosses `2` at `t = 28.795907 deg`, but a chain of three at `45 deg` admits only
   `L = 2` transverse links (`12/12` measured).  An H with `L = 2` spans **five** squares, and
   `(10 - 7 sqrt2)/2` is *also* `delta*` of five unit squares at `45 deg` in `[0,3]^2`
   (`S6_SKELETON.md` §3.2): **an H-shaped certificate at `T = 3` in the far field is literally a
   five-square certificate, and five squares do fit.**

4. **The far-field mechanism that is missing is a CYCLE.**  Taking the two mechanisms together,
   `H <= 0` **or** `MT <= 0` at `716 / 720`; the remaining `4` are far-field mixed-sign large
   tilts (`20-45 deg`) where no common-normal chain exists *and* the H is short.  At **all `68`**
   H-failures, adding back the pair rows of **one further main-direction link** makes the bound
   `<= 0`, and at `62` of them it recovers `delta*` exactly
   (`search/t3_chain_repair.py`, `runs/t3_chain_repair1.txt`).  One extra edge on a tree is one
   cycle.  At `theta = (45 deg)^6` the LP's own optimal dual is a **5-cycle** of pair rows at
   weight `1` on five of the six squares, anchored by four wall rows (one per wall) at weight
   `sqrt2`, giving `delta <= (93 - 66 sqrt2)/7` exactly — a genuine Farkas certificate of the same
   value as the rise-reading MT.  Over the whole scan the cyclomatic number of the true dual is
   `0` in `181/720`, `1` in `464`, `>= 2` in `75`, and `>= 1` at **all `48`** samples at exactly
   `45 deg`.  **[measured]**, §3.

5. **Kearney–Shiu is not a chain proof.**  Their `s(6) = 3` has no chain, no rise and no
   angle-space case split at all: it is a counting argument on `7 + 7 - 1 = 13` points of two
   unavoidable lattices, and the tilt is quantified away *once*, in their inequality (7), before
   any of the combinatorics starts.  In the dual language: their "chain" is the interval
   bookkeeping on the lines `x = 1`, `x = 2`, `y = 2`, their "rise bound" is the **joint**
   inequality `(1+t^2)/(1+t) + (1+2t-t^2)/2 >= 3/2` (the two summands have their minima at
   *different* angles, which is the whole content), and their far field is handled by the fact
   that unavoidability is angle-free.  Their `s(7) = 3` — seven squares, seven points, one
   bijection, one lemma — is the part with a clean `T = 4` analogue; their `s(6) = 3` is not,
   because the slack "`n+1` points, `n` squares" is `1` only at `T = 3`.  (§4.)

6. **Filter.**  The closed-form inequality `L >= L*(T, t)` is *`T`-uniform and gets easier as `T`
   grows* (`L*/L_max = 0.83, 0.54, 0.41, 0.17, 0.11` at `45 deg` for `T = 3, 4, 5, 11, 17`), so it
   **cannot** be the statement that fails at `T = 11`.  The statement that must fail there is the
   **existence** clause, and `T` enters it as the chain length: *every configuration of
   `T^2 - T` unit squares in `[0,T]^2` with margin `>= 0` contains a wall-to-wall chain of `T`
   carrying `L*(T, t)` transverse links.*  Cantrell's `T = 11` packing has
   `delta = +2.749e-04 > 0` (`BANDCUT_SCAN.md`, unrefereed), so some existence clause must fail
   there; **which one was not measured here** (§6.2).  The
   `s(5) = 3` filter is passed by a *different* clause, and this is the sharp point: five squares
   in `[0,3]^2` **do** carry wall-to-wall chains of three (measured: two of them at
   `theta = (0,0,0,0,45 deg)`, where `delta* >= +0.0858`), but every such chain runs **through the
   tilted square**, so each of its two links pays the mixed-tilt penalty
   `(W(45 deg) - 1)/2 = (sqrt2 - 1)/2 = 0.2071` and Lemma H returns `+0.1464 > 0`.  It is the
   *common-tilt* clause of Lemma H, not chain existence, that stops it proving `s(5) = 3` (§6).
   **[measured]** for both.

**What in the brief's premises turned out to be wrong.**

* **"The rise is bounded by transverse chains, so the lemma to prove is H-shaped."**  The first
  half is right and the second does not follow.  Eliminating the rise *by transverse chains* is
  one way to do it and it is lossy: `H - MT = 0.0986` at `T = 3`, `45 deg`.  The LP itself
  eliminates the rise a different way — by a **cycle**, which closes the transverse displacement
  around a loop instead of pinning its two ends against walls — and that way is exactly tight.
* **"Where H fails or no chain exists."**  A wall-to-wall chain of three exists at **every one**
  of the `720` sampled optima (`nochain = 0` in every family).  What can be missing is a chain
  with a **common separating normal**: that fails at `118 / 720`, in the mixed-sign small-tilt and
  the mixed-sign large-tilt families, never in a uniform one.
* **Erratum to an existing instrument (found here).**  `search/bandcut_k.py chain_certificate`
  builds its chain from the x/y-DAG and then evaluates `B_T` at the *mean* tilt of the chain,
  without checking that the links share a normal.  At `theta = (45 deg)^6` the x/y typing is
  degenerate (`|cos| = |sin|`), and the function returns `B_3 = -0.796` and `-1.034` for chains
  whose links differ by `180 deg` in normal direction, against a true margin of `-0.0483` — i.e.
  **invalid (too strong) bounds**.  It also reports `tilt_deg` through `math.degrees` twice
  (`-2578.31 deg` at `45 deg`).  The `B_T` columns of `BANDCUT_K.md` §4 are computed by this
  function; at the small tilts of that table (`<= 20 deg`) the typing is not degenerate and the
  numbers are believable, but nothing from it should be quoted at large tilt.  Nothing in
  `search/` was modified here.

---

## 1. The H-lemma, stated exactly

### 1.1 The LP and what a certificate is

Semantics as `notes/s13-casefree.md` §1 (closed squares, closed containment, packing = pairwise
disjoint as closed sets).  Fix the angles; by `S6_SKELETON.md` §3.1 the question is a disjunctive
LP in the `2n` centre coordinates, with the walls carrying `delta` (`BANDCUT_K.md` §1.1).  Every
row has the form `a . c - delta >= b`:

| row | `a` | `b` |
|---|---|---|
| `lo-x(i)` | `+e_{x_i}` | `P_i` |
| `hi-x(i)` | `-e_{x_i}` | `P_i - T` |
| `lo-y(i)`, `hi-y(i)` | `+e_{y_i}`, `-e_{y_i}` | `P_i`, `P_i - T` |
| `pair(i,j,o,k,s)` | `s n_{o,k}` at `j`, `-s n_{o,k}` at `i` | `m_ij` |

with `P_i = (|cos th_i| + |sin th_i|)/2`, `n_{o,0} = (cos th_o, sin th_o)`,
`n_{o,1} = (-sin th_o, cos th_o)`, `m_ij = 1/2 + (|cos D| + |sin D|)/2`, `D = th_j - th_i`.

> **Farkas form.**  **[proved]**  Let `w >= 0` be supported on rows that a configuration `z`
> satisfies, with `sum_r w_r a_r = 0` and `sum_r w_r = 1`.  Then `delta(z) <= -sum_r w_r b_r`.

(Sum the rows with weights `w`; the centre terms cancel and `-delta` is left with coefficient
`1`.)  Conversely the best such `w` over *all* satisfied rows reproduces `delta*` exactly when `z`
is the global optimum for its angle vector — verified in every sample of §2 (`fulldual` column,
agreement to `1e-13`).  So "restricting the support" is exactly "restricting which geometric facts
the proof is allowed to use", and the gap between the restricted and the full value is the price of
the restriction, in the units of the theorem.

### 1.2 The H support

> **Definition (H).**  Fix an axis `ax in {x, y}` and let `tr` be the other one.  An **H** consists
> of
>
> * a **crossbar**: a wall-to-wall chain `s_1 -> s_2 -> ... -> s_T` in the `ax`-DAG at the
>   configuration's own margin (`BANDCUT_K.md` §1.2: edge `i -> j` when some `ax`-type edge normal
>   separates the pair with `j` on the `+` side and gap `>= delta`; staircases allowed), together
>   with the two wall rows `lo-ax(s_1)`, `hi-ax(s_T)`;
> * two **legs**: a `tr`-chain `b_1 -> ... -> b_a -> s_1` running down to the `lo-tr` wall, and a
>   `tr`-chain `s_T -> d_1 -> ... -> d_c` running up to the `hi-tr` wall, with the wall rows
>   `lo-tr(b_1)`, `hi-tr(d_c)`, and the two closing wall rows `hi-ax(b_1)`, `lo-ax(d_c)`.
>
> Write `L = a + c` for the number of transverse links.  All `T + L` squares are distinct.

Two remarks the brief's phrasing does not anticipate.

* **The legs sit at the crossbar's two ends, not at its members.**  For a chain whose links all
  carry the same normal the Farkas residual after summing the crossbar is supported **only on
  `c_{s_1}` and `c_{s_T}`** (§1.3); the interior members cancel by themselves.  Legs through
  interior members are needed only for a staircase, i.e. when the link normals differ.
* **The four closing wall rows `hi-ax(b_1)`, `lo-ax(d_c)` are not optional.**  They carry weight
  `sin t tan t = O(t^2)`, and without them the support admits no cancelling combination at all.
  That is why a support restricted "to the chain plus transverse chains" with no wall rows is
  infeasible (mode `C` of `t3_chain.py`: `+inf` at every tilted sample).

### 1.3 The closed form  [proved]

> **Lemma H.**  Let the `T + L` squares of an H all have the same tilt `t in [0, 45 deg]`, let
> every crossbar link be separated by the normal `e = (cos t, sin t)` and every leg link by
> `f = (-sin t, cos t)`.  Put `C = cos t`, `S = sin t`, `u = C + S` and
>
>     M(t)  =  (1 + S^2 + S C) / C      ( = 1/C + S + S^2/C ,   M(0) = 1 ) .
>
> Then the weights
>
>     crossbar links (T-1 rows)            1
>     lo-ax(s_1), hi-ax(s_T)               1/C
>     leg links (L rows)                   tan t
>     lo-tr(b_1), hi-tr(d_c)               S
>     hi-ax(b_1), lo-ax(d_c)               S tan t
>
> cancel every centre coordinate exactly, and give
>
>     delta  <=  H(T, L, t)  :=  [ (T - u) M(t) - (T - 1) - L tan t ]
>                                / [ (T - 1) + L tan t + 2 M(t) ] .

*Proof.*  Coefficient bookkeeping, one square at a time.  The crossbar with all links at weight
`1` leaves `(1 - C, -S) . ` at `c_{s_1}` ... written out: the residual is `-e` at `c_{s_1}` and
`+e` at `c_{s_T}`, cancelled in `ax` by `1/C` on the two `ax`-walls up to the `S`-part, which the
legs supply; the leg at weight `tan t` contributes `+ (S/C) f` at `c_{s_1}`, whose `tr`-component
`S` cancels the crossbar's `-S` and whose `ax`-component `-S^2/C` is what pushes the `ax`-wall
weight from `C` to `1/C`.  At the leg's far end the residual is `(S^2/C, 0)` resp. `(-S^2/C, 0)`,
cancelled by the two closing walls, and `(0, -S)` resp. `(0, S)`, cancelled by the two `tr`-walls.
Interior squares of the crossbar and of the legs cancel because consecutive links share a normal
and a weight.  `sum w = (T-1) + L tan t + 2M`; `sum w b = (T-1) + L tan t + (u - T) M`. ∎

`search/t3_chain_hform.py check` re-derives the weight vector square by square and reports the
maximum residual coefficient and the agreement with the closed form: `OK` (residual `< 1e-12`) for
`T in {3, 4, 5, 11}`, `a, c in {0, 1, 2}`, `t in {0.5, 5, 20, 45} deg`.

> **Corollary (the whole lemma is one integer).**  **[proved]**
> `H(T, L, t) <= 0`  iff  `L >= L*(T, t) := ( M(t) (T - u) - (T - 1) ) cot t`, with
> `L*(T, 0) = T - 2` (a removable singularity; `L*` is increasing in `t`).

`runs/t3_chain_hform.txt`:

| `T` | `L_max = n - T` | `L*` at `0` | `1 deg` | `5 deg` | `10 deg` | `20 deg` | `30 deg` | `45 deg` |
|---|---|---|---|---|---|---|---|---|
| 3 | 3 | 1.000 | 1.043 | 1.207 | 1.396 | 1.731 | 2.036 | 2.485 |
| 4 | 8 | 2.000 | 2.069 | 2.334 | 2.642 | 3.189 | 3.670 | 4.314 |
| 5 | 15 | 3.000 | 3.095 | 3.461 | 3.888 | 4.647 | 5.304 | 6.142 |
| 11 | 99 | 9.000 | 9.251 | 10.223 | 11.363 | 13.396 | 15.108 | 17.113 |
| 17 | 255 | 15.000 | 15.408 | 16.986 | 18.839 | 22.144 | 24.912 | 28.083 |

and at `T = 3` the bound itself:

| `t` | `L = 0` | `L = 1` | `L = 2` | `L = 3` |
|---|---|---|---|---|
| `0` | `0` | `0` | `0` | `0` |
| `1 deg` | `+0.004512` | `+0.000186` | `-0.004103` | `-0.008355` |
| `5 deg` | `+0.025167` | `+0.004235` | `-0.015860` | `-0.035167` |
| `20 deg` | `+0.124503` | `+0.049063` | `-0.016891` | `-0.075043` |
| `28.8 deg` | `+0.190` | `+0.0868` | `~0` | `-0.0741` |
| `45 deg` | `+0.324583` | `+0.171573` | **`+0.050253`** | **`-0.048299`** |

**`L*(3, t) = 2` at `t = 28.795907 deg`** (`runs/t3_chain_filter.txt`).  Below it a chain of three
with two legs — five squares — suffices; above it the H must contain **all six**.

### 1.4 The same thing with the rise eliminated: the `omega` identity  [proved]

> **Lemma W (uniform tilt).**  If every square has tilt `t`, so that every `P_i = u/2` and every
> `m_ij = 1`, then **every** normalised Farkas certificate satisfies
>
>     delta  <=  omega (T + 2 - u)  -  1 ,
>
> `omega` = the total weight the certificate puts on `lo` wall rows (which equals the total on
> `hi` rows).  In particular `delta* <= 0` iff some certificate has `omega <= 1/(T + 2 - u)`.

*Proof.*  Summing the `x`-coefficient equations over all squares kills every pair row (its two
coefficients are opposite) and leaves `w(lo-x) - w(hi-x) = 0`; same in `y`.  So
`W_lo = W_hi = omega` and `W_pair = 1 - 2 omega`.  Then
`sum w b = omega (u/2) + omega (u/2 - T) + (1 - 2 omega) = 1 - omega (T + 2 - u)`. ∎

At `t = 0` the threshold is `omega <= 1/(T+1)`, which is exactly the chain certificate's
`1/(T+1)` per row (`BANDCUT_K.md` §5).  The H of Lemma H has
`omega = M / ((T-1) + L tan t + 2M)`, and **each extra link drives `omega` down** — the whole
content of "more transverse structure" is "more pair rows per unit of wall weight".  This is a
cleaner statement of the mechanism than the rise: the rise never appears.

### 1.5 Mixed tilts  [proved]

Lemma H survives verbatim if the link **normals** are all owned by squares of tilt `t` while the
*other* endpoints have different tilts: the cancellation is unchanged, only the `b`'s move.
Writing `m_s` for the crossbar links and `m_r` for the leg links and `P` for the four anchors,

    delta  <=  [ T M - (1/C)(P_{s_1} + P_{s_T}) - (S + S^2/C)(P_{b_1} + P_{d_c})
                 - sum_s m_s - tan t * sum_r m_r ]  /  [ (T-1) + L tan t + 2 M ] .

Each tilt-mismatched crossbar link costs `m_s - 1 = (W(D) - 1)/2 ~ |D|/2` — **first order**, as
`notes/bandcut-cost.md` §1 says — at weight `1`, and each mismatched leg link the same at weight
`tan t`.  At `t = 0` this degenerates to `delta <= (T - P_{s_1} - P_{s_T} - sum_s m_s)/(T+1)`, in
which `L` has disappeared: **an axis-parallel crossbar kills on its own and needs no legs**, and
a crossbar with one `45 deg` member does not kill at all (§6).  If the link normals themselves
spread the sum stops telescoping and the certificate dies (`bandcut-cost.md` §1, the `sin(beta)
(k-1) tau` repair term); Lemma H is the exact statement of the case where they do not.

---

## 2. The test over the whole angle space  [measured]

### 2.1 What was sampled and how

`search/t3_chain.py scan --reps 12 --limit 300 --extra 600 --nproc 14`:
**60 families x 12 samples = 720 angle vectors**, 321 s on 14 processes
(`runs/t3_chain_scan1.jsonl`).  Per angle vector: `delta*` by assignment-fixed LP ascent
(HiGHS, tolerances `1e-10`) from 900 jittered tiling starts with random label dealing plus 600
uniform random starts; then, at the best configuration `z` found, the full row system of §1.1
restricted to rows `z` satisfies at its own margin, and four dual LPs on nested supports:

| mode | support |
|---|---|
| `C` | the chain's links and its two end walls only |
| `CW` | `C` plus **every** wall row |
| `H` | `CW` plus **every transverse-type pair row** — this is the H-lemma, best over all chains |
| `F` | everything; equals `delta*` when `z` is the global optimum |

`F` agreed with `delta*` to `1e-13` at every sample, which is the instrument's self-check.  `H` is
maximised over every wall-to-wall chain of three present at the configuration's own margin (both
DAGs, staircases included).  **`H` is generous**: it is given *all* transverse pair rows and *all*
wall rows, so it is an upper bound on what any particular H of §1.2 can achieve; when `H > 0` no
H-shaped certificate exists at that configuration at all.

The families (all at `T = 3`, `n = 6`): `Z3..Z6` (`j` angles exactly `0`, the rest uniform on
`[0,90 deg)`); `near eps` (all six tilts in `[0, eps]`, random signs) and `coh eps` (same, one
sign) for `eps = 0.5, 1, 2, 5, 10, 20 deg`; `hole k_eps` for `k = 0,1,2` (far field) and
`k = 3, 4` (the hole) at `eps = 1, 2, 5, 10, 20 deg`; `unif t` for
`t = 1, 5, 10, 20, 30, 40, 44, 45 deg`; `split a_t` (`a` at `+t`, `6-a` at `-t`) for
`a = 1,2,3`, `t = 2, 5, 20, 45 deg`; `generic` (uniform on the whole cube); `bigrand eps` (all six
tilts `>= eps`, random signs).

*Reliability.*  Every `delta*` is a feasible point, so a **lower** bound; an under-found `delta*`
makes the H test **easier** (a lower margin admits more chain edges), so the failures below are
not an artefact of a weak search — but the successes may be.  The three families where `delta*`
is known independently reproduce: `unif1 = -1.1079e-04`, `unif5 = -2.4500e-03`,
`unif20 = -2.4370e-02` match `BANDCUT_K.md` §2(a) `k = 0` to all printed digits.

### 2.2 Result

**`H <= 0` at `652 / 720` sampled optima (90.6 %).  It fails at `68`.**  The failures are not
spread out; they are two identifiable regimes.

| regime | families | fails | character |
|---|---|---|---|
| **large coherent tilt** | `unif30`, `unif40`, `unif44`, `unif45` | `48 / 48` | systematic, `H - delta*` up to `+0.0986` |
| **the coherent small-angle cone** | all six `coh*` | `10 / 72` | `H > 0` by a *third*-order amount |
| far field, scattered | `hole0_5, hole0_20, hole1_20, hole3_20, bigrand20, bigrand30, generic` | `10 / 600` | mixed signs and large tilts |
| everything else | the other `43` families (`Z3..Z6`, all `near*`, all `split*`, `unif <= 20`, most `hole*`, `bigrand5/10`) | `0 / 516` | — |

Per-family hit rates (`runs/t3_chain_scan1.jsonl`, `t3_chain.py report`), abbreviated:

| family | `N` | max `delta*` | median `H` | max `H` | `H <= 0` | `CW <= 0` |
|---|---|---|---|---|---|---|
| `Z3` / `Z4` / `Z5` / `Z6` | 12 each | `0` | `-5.5e-3 / 0 / 0 / 0` | `0` | `12/12` each | `10,12,12,12` |
| `near0.5 .. near20` | 12 each | `-4.4e-6 .. -1.1e-2` | | `-4.4e-6 .. -9.5e-3` | `12/12` each | `0-3 / 12` |
| `coh2` | 12 | `-7.6e-5` | `-1.7e-5` | `+5.1e-5` | **`7/12`** | `0/12` |
| `coh5` | 12 | `-3.8e-4` | `-3.6e-4` | `+1.2e-4` | `10/12` | `0/12` |
| `hole3_*`, `hole4_*` (the hole) | 12 each | `-8.0e-4 .. -2.2e-2` | | | `11-12 / 12` | `0-2 / 12` |
| `hole0_*`, `hole1_*`, `hole2_*` (far field) | 12 each | `-1.3e-2 .. -5.4e-2` | | | `10-12 / 12` | `0-5 / 12` |
| `unif1 / 5 / 10 / 20` | 12 each | `-1.1e-4 .. -2.4e-2` | `-3.6e-5 .. -4.3e-4` | same | `12/12` | `0/12` |
| `unif30 / 40 / 44 / 45` | 12 each | `-3.9e-2 .. -4.8e-2` | `+2.9e-3 .. +5.0e-2` | same | **`0/12`** | `0/12` |
| `split1_45 / 2_45 / 3_45` | 12 each | `-4.8e-2` | `-4.8e-2` | same | `12/12` | `0/12` |
| `generic` | 12 | `-2.5e-2` | `-3.1e-2` | `+5.5e-2` | `10/12` | `2/12` |
| `bigrand5 / 10 / 20 / 30` | 12 each | `-2.6e-2 .. -5.6e-2` | | | `12,12,10,11` | `0-4 / 12` |

Two things to read off.

* **`CW` (chain + walls, no transverse pair rows) is useless.**  It is `<= 0` in only `191/720`
  samples, essentially only where the chain is axis-parallel (`Z4..Z6`, `split*` at small `t`).
  So the transverse chains really are the load-bearing part of the H, exactly as the brief says.
* **Mixed signs are easy and coherent signs are hard.**  `split a_t` (some tilts `+t`, some `-t`)
  is `12/12` at every `a` and every `t` up to `45 deg`, while `unif t` fails from `30 deg`.  The
  reason is `S6_LOCAL.md` §5.1's: a link between oppositely-tilted squares costs `|D|/2` at first
  order, so mixed-sign configurations are already far below `0` and any certificate reaches them.

### 2.3 The large-tilt failure, exactly

At `theta = (45 deg)^6` the measurement is not a measurement at all — both sides are algebraic.

    delta*( 45 deg, ..., 45 deg )  =  (93 - 66 sqrt2) / 7  =  -0.048299302374896...   [measured to 15 digits,
                                                                    and exact on its leaf, sec 3.2]
    best H-restricted dual bound   =  (10 - 7 sqrt2) / 2   =  +0.050252531694167...

and `(10 - 7 sqrt2)/2` is **three different things at once**:

1. `H(3, 2, 45 deg)` from the closed form of §1.3 — a chain of three with two transverse links;
2. the value of the best H-restricted dual over the six-square `45 deg` optimum (`t3_chain.py`);
3. `delta*` of **five** unit squares at `45 deg` in `[0,3]^2` (`S6_SKELETON.md` §3.2 reports
   `+0.050252532` from the independent `fixed_angle_value` instrument; re-measured here as
   `+0.050252532` with the chain-aware instrument, `runs/t3_chain_filter.txt`).

That triple identity is the content of the failure.  At `45 deg` a chain of three admits at most
`L = 2` transverse links (`t3_chain_exist.py`: `L = 2` at `12/12` `unif45` samples, and `L* = 2.485`),
so the H spans **five** squares; the sixth can only be attached by a link in the *main* direction,
which the H forbids.  **An H-shaped certificate at `T = 3` in the far field is literally a
five-square certificate, and five unit squares do fit in `[0,3]^2`.**

The crossover is exact: `L*(3, t) = 2` at `t = 28.795907 deg`, and the scan brackets it —
`unif20` is `12/12`, `unif30` is `0/12`.  (`L*(3, 28) = 1.976`, `L*(3, 29) = 2.006`.)

The decisive line, measured side by side (`runs/t3_chain_filter.txt`):

    5 squares at 45 deg in [0,3]^2:  delta* = +0.050252532 ,  best available L = 2  (chain 1-4-2, a = c = 1)
    6 squares at 45 deg in [0,3]^2:  delta* = -0.048299302 ,  best available L = 2  (chain 0-3-1, a = c = 1)

**Both admit exactly `L = 2`, and `H(3, 2, 45 deg) = +0.050252532` is attained by the five-square
one.**  So the H with `L = 2` is an exactly tight certificate for five squares and therefore
cannot say anything at all about six; and `L = 3` is not available.  This is not a numerical
accident and it is not fixable by sharpening the H's weights.

### 2.4 The small-angle failure

In the coherent cone (`coh2`: all six tilts of one sign in `[0, 2 deg]`) `H` fails at `5/12` with
`H` in `[+1.2e-5, +5.1e-5]` against `delta*` in `[-2.3e-4, -1.4e-4]` — both **third** order in the
tilt (`t^3 = 4.2e-5` at `2 deg`), which is `S6_LOCAL.md` §2's `-(1/4) t^3` direction.  Here `L = 2`
is available and `L*(3, 2 deg) = 1.086`, so Lemma H's *inequality* is satisfied with room; what
fails is Lemma H's *common-tilt hypothesis*.  The exact dual at such a point (read at
`theta = (5 deg)^6`, `t3_chain.py`) is

    0.2342 lo-x(0)   0.2342 hi-x(2)   0.2333 pair(0,5) [x-link]   0.2333 pair(2,5) [x-link]
    0.0205 lo-y(3)   0.0205 hi-y(2)   0.0204 pair(0,3) [y-link]   0.0018 pair(2,4) [x-link]

— a chain of three at `~1/4 - O(t)`, a transverse chain at `~(1/4) tan t = 0.0219`, **and one more
main-direction link at `(1/4) tan^2 t = 0.0019`**.  That last row is `S6_LOCAL.md` §3's "two links
at `t^2/2`"; it is a main-direction row off the chain, so the H support does not contain it, and
without it the equality system has no cancelling solution of comparable value.  The H therefore
loses an amount of order `t^3`, which is exactly the size of the margin in that cone.

### 2.5 The minimal repair: one rung  [measured, 68/68]

`search/t3_chain_repair.py` adds back, to each failing H, the main-direction pair rows of **one**
further pair, minimised over the fifteen pairs (`runs/t3_chain_repair1.txt`):

> **At all `68 / 68` failures, one extra main-direction link turns `H > 0` into `<= 0`**, and at
> `62` of them it recovers `delta*` **exactly**.  At `45 deg`: `H = +0.050253` becomes
> `-0.048299 = delta*`.  In the `coh2` cone: `+3.3e-5` becomes `-1.43e-4 = delta*`.

One extra edge on a tree makes exactly one cycle.  So the object that works everywhere at `T = 3`
is not an H but an **H with a rung** — a `Theta`, cyclomatic number `1`.

### 2.6 Under-explored cells

* **`eps > 20 deg` with mixed signs** is thin: `bigrand30` and `split*_45` are the only families
  above `30 deg` with sign freedom, `12` samples each.
* **The `Z3` stratum** (exactly three angles `0`, three free) got `12` samples; `S6_SKELETON.md`
  §4.1 reports only `3/80` of that stratum reaches `delta* = 0`, so the interesting part of `Z3`
  is a thin subset that random sampling will rarely hit.
* **Configurations that are not `delta*`-optimal** were not tested at all.  The H-lemma is a
  statement about *every* configuration with `delta >= 0`; testing it at optima only is the weakest
  possible reading, and it is the reading under which it already fails.
* **`n = 12`, `T = 4` was not run here.**  The closed form extends (§1.3 table) but nothing was
  measured.

### 2.7 MT tested correctly — and why it is not a certificate  [measured / proved]

`search/t3_chain_mt.py` enumerates, for each sampled optimum, every unit normal `e` along which
some ordered pair is separated at the configuration's own margin, builds the DAG **of that one
normal**, and evaluates MT (`notes/bandcut-cost.md` §1) on every path of three, for both wall
pairs.  Every value it prints is then a valid upper bound on `delta` (checked: the worst
`MT - delta*` over all 720 samples is `-2.0e-08`, i.e. zero to the `1e-7` edge tolerance).

    MT  <=  0  at  602 / 720 .   At the other 118 there is NO chain of three with a common normal.

| family group | `MT <= 0` | `H <= 0` |
|---|---|---|
| `Z3..Z6` | `12/12` each (`MT = 0`) | `12/12` each |
| `unif1 .. unif20` | `12/12` each, `MT = delta*` to `1e-9` | `12/12` each |
| `unif30, 40, 44, 45` | **`12/12` each, `MT = delta*` exactly** | **`0/12` each** |
| `split*` (mixed signs, any `t`) | `12/12` each | `12/12` each |
| `near0.5 .. near20` (mixed signs, small `t`) | `4-9 / 12` | `12/12` each |
| `coh1, coh2, coh5, coh10, coh20` | `6-10 / 12` | `7-12 / 12` |
| `bigrand20, bigrand30, generic, hole0_*` | `5-9 / 12` | `10-12 / 12` |

**The two mechanisms are complementary.**  `H <= 0` **or** `MT <= 0` holds at `716 / 720`; the
four residual samples are far-field mixed-sign configurations at `20-45 deg` with
`delta* in [-0.059, -0.050]`, and every one of them is repaired by the single rung of §2.5.

But MT is not a certificate.  Its derivation substitutes the configuration's own
`R = (c_k - c_1)_perp` into the chain sum; the transverse coordinates are never cancelled, so MT
proves `delta <= MT(R)` **for a configuration whose rise is `R`**, not for the angle vector.  To
turn it into a statement about `theta` one must bound `R`, and there are exactly two ways:

1. **pin the two ends against the transverse walls through transverse chains** — the H.  Cost:
   `H - MT`, which is `+0.0986` at `T = 3`, `45 deg`, `L = 2` (larger than `|delta*|` itself).
2. **close the transverse displacement around a loop** — the cycle.  Cost: `0` at `45 deg` (the
   pentagon of §3.2 reproduces `MT = delta*` exactly).

That is the single sharpest thing this task found.


---

## 3. The far-field mechanism: the dual is a cycle

### 3.1 Shape statistics  [measured]

`search/t3_chain_shape.py` builds, for each sampled optimum, the graph whose vertices are the
squares touched by the dual and whose edges are the pair rows in its support, and reports the
cyclomatic number `mu = E - V + components`.  **An H is a tree: `mu = 0`.**
(`runs/t3_chain_shape1.txt`.)

| max tilt of the configuration | `mu = 0` | `mu = 1` | `mu >= 2` |
|---|---|---|---|
| `0-5 deg` | 57 | 87 | 13 |
| `5-10` | 27 | 37 | 10 |
| `10-20` | 6 | 31 | 7 |
| `20-30` | 30 | 77 | 3 |
| `30-40` | 31 | 76 | 12 |
| `40-45` | 30 | 112 | 26 |
| exactly `45` | **0** | 44 | 4 |
| all | 181 | 464 | 75 |

Median support: `2` pair rows on `3` squares in `Z5`/`Z6` (the bare axis-parallel chain of three),
`4-7` pair rows on `4-6` squares everywhere else.  Every `unif` family from `5 deg` up is `12/12`
`mu >= 1`.  The `mu = 0` duals at large tilt are the ones where an *axis-parallel* chain still
exists (mixed-tilt families).

### 3.2 The `45 deg` certificate, exactly  [proved]

At `theta = (45 deg)^6` the optimum found is (up to the symmetry group)

      sq0 (0.65881, 1.33176)   sq1 (1.66824, 1.66824)   sq2 (0.99528, 2.34119)
      sq3 (1.33176, 0.65881)   sq4 (2.34119, 2.34119)   sq5 (2.34119, 0.99528)

and its dual has exactly nine rows:

    0.132705  lo-x(0)    0.132705  hi-x(4)    0.132705  lo-y(3)    0.132705  hi-y(2)
    0.093836  pair(0,1)  0.093836  pair(1,4)  0.093836  pair(2,3)  0.093836  pair(0,3)
    0.093836  pair(2,4)

The five pair rows form the **5-cycle** `3 - 0 - 1 - 4 - 2 - 3` on five of the six squares; `sq5`
is idle.  Each square has two edges of the cycle; the four that also carry a wall row are exactly
the four at which the cycle's link normal **switches** between the two edge normals
`n_0 = (1,1)/sqrt2` and `n_1 = (-1,1)/sqrt2`.  That is the whole mechanism, and it is one line:

> at a vertex with an incoming link of normal `a` and an outgoing link of normal `b`, both at
> weight `1`, the residual is `a - b`; if `a = b` it vanishes, and if `a, b` are the two edge
> normals of a `45 deg` square then `a - b = (+-sqrt2, 0)` or `(0, +-sqrt2)`, which one wall row
> at weight `sqrt2` cancels.

So the certificate has five pair rows at weight `1` and four wall rows at weight `sqrt2`, total
`5 + 4 sqrt2`, wall mass `omega = 2 sqrt2 / (5 + 4 sqrt2)`.  Lemma W (§1.4) then gives, with
`u = sqrt2`, `T = 3`:

    delta  <=  omega (T + 2 - u) - 1  =  2 sqrt2 (5 - sqrt2) / (5 + 4 sqrt2)  -  1
            =  (6 sqrt2 - 9) / (5 + 4 sqrt2)  =  **(93 - 66 sqrt2) / 7**  =  -0.0482993023748963

which matches the measured `delta*` to all 15 digits.  **[proved]** as a Farkas bound for every
configuration containing that cycle; **[measured]** that it is the value of `delta*`.

Contrast the H of the same configuration: `L = 2` transverse links are available, and
`H(3, 2, 45 deg) = (10 - 7 sqrt2)/2 = +0.0502525` — the five-square value.  **The cycle's fifth
pair row is worth `2 omega (T + 2 - u) = 0.0985` of margin**, because it enters the denominator of
`omega` at weight `1` while a transverse leg link enters at weight `tan t`.

### 3.3 What the cycle is, structurally

In every cyclic dual read here the same three features recur.  **[measured]**, `n = 6` only:

1. **Four wall rows, one per wall.**  An H uses `lo-ax`, `hi-ax`, `lo-tr`, `hi-tr` plus two
   closing rows on the same axis as the crossbar; a cycle uses exactly one of each, and no square
   ever carries two opposite walls.
2. **The cycle winds once.**  Following the link normals around the cycle, the normal direction
   advances by `2 pi` in total; the number of "turns" (vertices where the normal changes) equals
   the number of wall rows, i.e. four.  At `45 deg` that is a pentagon with four turns; at smaller
   tilt the same picture appears as a chain of three plus a rung (§2.5), which is a quadrilateral
   or pentagon with two of its sides nearly straight.
3. **It does not need all `n` squares.**  At `45 deg` it uses five of six.  This is the opposite
   of the H, which at `45 deg` would need all six and cannot get them.

So the far-field mechanism the repo is missing is: **a closed loop of separations winding once
around, pinned at all four walls** — the dual of the pinwheel, not of the row.  It is the same
object `S6_SKELETON.md` §4.2 and `S6_LOCAL.md` §3 call "the pinwheel certificate", now identified
as a cycle rather than as a chain-plus-corrections, and now known to be the *only* thing that
works above `28.8 deg`.

### 3.4 What a T = 4 analogue would be  [guess]

Nothing was measured at `T = 4`.  The uniform-tilt identity (Lemma W) is `T`-free, so the question
is entirely "which cycles exist and what is their `omega`".  For a cycle of `k` links with `q`
turns at a common tilt `t`, the same vertex bookkeeping gives link weight `1` and turn weight
`2 sin(alpha/2)` where `alpha` is the turn angle, so `omega = q sin(alpha/2) / (k + 2 q sin(alpha/2))`
and `delta <= 0` iff `k >= q sin(alpha/2) (T - u)`.  At `T = 3, t = 45 deg, q = 4, alpha = 90 deg`
that reads `k >= 2 sqrt2 (3 - sqrt2) = 4.49`, i.e. `k >= 5` — the pentagon, on the nose.  At
`T = 4` it would read `k >= 2 sqrt2 (4 - sqrt2) = 7.31`, i.e. a **cycle of eight** among twelve
squares.  Untested.

---

## 4. Kearney–Shiu 2002 in this language

Source: Electron. J. Combin. 9 (2002) #R14, doi 10.37236/1631, read in full for this note
(a transcript of the relevant passages is in `notes/proof-anatomy.md` §5.2 in summary; the
details below are from the paper itself).  Their §2 is `s(7) = 3`, their §3 is `s(6) = 3`.

### 4.1 What the proof actually is

**There is no chain, no rise, and no angle-space case split anywhere in it.**  The machinery is:

* an **unavoidable set** of 7 points in `[0,3]^2` (Friedman's, sharpened):
  `{(sqrt2-1/2, 1), (3/2, 1), (7/2-sqrt2, 1), (3/2, 3/2), (sqrt2-1/2, 2), (3/2, 2), (7/2-sqrt2, 2)}`
  — every unit square inside `[0,3]^2` contains one of them;
* its **90 deg rotation** about the centre `(3/2, 3/2)`, which is unavoidable too *because the
  container is a square*: the "green" and "red" lattices, sharing the centre `C = (3/2, 3/2)`, so
  `7 + 7 - 1 = 13` distinct points, classified as `1` C-point, `4` B-points at distance `1/2` from
  `C`, `8` A-points;
* **Lemma 1** (their only fully proved lemma): *any unit square covering `C` also covers a
  B-point*.  Proof: put the **square** in standard position `[0,1]^2` and let `C = (x_0, y_0)` with
  `0 <= x_0 <= y_0 <= 1/2`; the circle of radius `1/2` about `C` meets the square's boundary in an
  arc subtending `>= 90 deg`, and the four B-points are `90 deg` apart on that circle;
* **Lemma 2** (proof omitted, "elementary coordinate geometry"): a unit square with a corner on the
  `x`-axis at angle `theta`, with `(0,1)` on the opposite edge, has the points
  `((1+t^2)/(1+t), 1)` and `(1, (1+2t-t^2)/2)` on two of its edges, `t = tan(theta/2)`;
* their **inequality (7)**, for `0 <= t <= 1`:
  `(1+t^2)/(1+t) >= 2 sqrt2 - 2`,  `(1+2t-t^2)/2 >= 1/2`,  and the **sum `>= 3/2`**;
* **Lemma 3** (proof omitted): a unit square covering `(3/2,3/2)` with `(1,2), (2,2), (2,3/2)` on
  three of its edges meets `x = 1` at height `<= 5/3`.

`s(7) = 3` is then one line: 7 squares, 7 green points, each square covers `>= 1`, so a bijection;
same in red; the square covering `C` has `C` as its unique green *and* unique red point, but
Lemma 1 forces it to cover a B-point too, which is a second point of one lattice.  Contradiction.

`s(6) = 3` is the same counting with **one unit of slack** (6 squares, 7 points, so at most one
square covers two points of a lattice), spent in a two-case argument:
(a) `C` uncovered — then the matching is perfect, and an intercept computation on `x = 2`
(Lemma 2 applied twice plus the *sum* inequality of (7), giving `>= 3/2`) shows that covering one
A-point forces covering a B-point;
(b) `C` covered — Lemma 1 forces a B-point with it; four pairings are then forced, the square at
`C` is boxed away from `(1,3/2), (1,2), (2,3/2), (2,2)`, Lemma 3 caps its edge on `x = 1` at
`5/3`, so the free interval on `x = 1` has length `<= 5/3 - (sqrt2 - 1/2) = 13/6 - sqrt2 =
0.75245`, while Lemma 2 + the first inequality of (7) demands `>= 2 sqrt2 - 2 = 0.82842`.
Contradiction by `0.076`.

### 4.2 The dictionary

| our object | their object |
|---|---|
| a wall-to-wall **chain** of `T`, `sum_s m_s + P_1 + P_T >= T` | the **interval bookkeeping on a line**: the lengths that squares must occupy on `x = 1`, `x = 2`, `y = 2` add up to more than is there.  Their final contradiction `13/6 - sqrt2 < 2 sqrt2 - 2` is exactly a chain inequality on the line `x = 1`, with `5/3` playing the role of the wall row |
| the **rise bound** (what stops `sum_s sin tau_s Dy_s` from being large) | the **joint** inequality of (7), `(1+t^2)/(1+t) + (1+2t-t^2)/2 >= 3/2`.  The two summands are two intercepts of the *same* square at the *same* tilt; each alone is minimised at a different angle (`t = sqrt2 - 1`, i.e. `45 deg`, and `t = 0`), so the sum cannot be at both minima at once.  That is precisely our "the chain cannot both be short and rise" |
| the **transverse chains** | nothing, or rather: *the second lattice*.  Where we bound a transverse displacement by a chain of separations in the other direction, they use the `90 deg`-rotated point set to get a second, independent counting constraint in the other direction |
| **chain existence** | **unavoidability** — every unit square covers a point of the lattice.  This is where all their angle analysis is spent (Friedman's Lemmas 1, 2, 3, 5, 6, 7, each minimised at `theta = 45 deg`), and it is done once, globally, before any case split |
| the **far field** (all squares well tilted) | **does not exist as a case.**  Unavoidability is angle-free, and (7) has quantified `theta` away before the combinatorics starts |
| `omega <= 1/(T + 2 - u)` (Lemma W) | the counting `6 squares, 7 points, slack 1` |

### 4.3 How they dispose of all-tilted configurations

By never letting the tilt into the combinatorics.  The only place `theta` appears as a variable in
§§2–3 is the statement of Lemma 2, and inequality (7) replaces its two `theta`-dependent
coordinates by the three `theta`-free constants `2 sqrt2 - 2`, `1/2`, `3/2`, uniformly over
`0 <= theta <= 90 deg`.  Everything after that is incidence combinatorics on 13 fixed points.
**This is the exact opposite of the architecture in this repo**, which carries `theta` all the way
into the certificate and then has to handle the far field as a separate regime.  The price they pay
is that the constants `2 sqrt2 - 2`, `1/2`, `3/2`, `5/3`, `sqrt2 - 1/2` are all *worst-case over
`theta`*, i.e. attained only at `45 deg` or at `0`, and the argument has to be tight enough to
survive that slack — which at `T = 3` it is, by `0.076`.

### 4.4 `T = 4` analogues

| step | `T = 4` analogue? |
|---|---|
| unavoidable set of `n+1 = 7` points in `[0,3]^2` | **the step under strain.**  `s(15) = 4` uses Friedman's `14`-point set in `[0,4]^2` (`DS7` §5, quoted in `proof-anatomy.md` §5.3); `s(12) >= 4` would need an unavoidable set of `11` points in `[0,4]^2`, and **no such set is known** (`DS7` Table 1 has no pure point set below `15`, and `proof-anatomy.md` §7.1 is the repo's record of why).  `11` is `5` below `T^2 = 16`, where `6` is `1` below `3^2 = 9`.  Nothing here proves an `11`-point set cannot exist; it is simply the step with no instance |
| the `90 deg` rotation / duality | **survives verbatim**: the container is a square at every `T`, so any unavoidable set gives a second one |
| Lemma 1 (`C` `=>` B-point) | **survives**: it is a statement about one unit square and a circle of radius `1/2`, with no `3` in it.  At `T = 4` the analogue would need the right `B`-ring around the relevant point |
| Lemma 2 + inequality (7) | **survives**: purely local, one square, one corner, `[0,1]^2` |
| Lemma 3 (`y <= 5/3` on `x = 1`) | **container-specific**: `5/3` is the only constant in the paper that knows `T = 3` |
| the counting `6 squares vs 7 points, slack 1` | **dies**: at `T = 4` the slack would be `16 - 12 = 4` (or, against the largest known unavoidable set for `[0,4]^2`, still `>= 2`), and a slack of `2` already destroys the "at most one square covers two points" step that both their cases rest on |

So the transferable part is the *local* lemma toolbox (Lemma 1, Lemma 2, (7)) and the duality; the
part that carries `T = 3` is the counting slack of exactly one, which is the same "one unit of
slack" that `notes/proof-anatomy.md` §7.2 identifies as the whole difficulty at `n = 12`.

---

## 5. The proof attempt

### 5.1 The skeleton, and what is actually proved

> **S0 (reduction).**  `s(6) >= 3` iff `delta*(theta) <= 0` for every `theta in [0, 90 deg)^6`.
> **[proved]**, `S6_SKELETON.md` §3.1, reviewer-checked (`proof-architecture.md` §0a item 5).

> **S1 (Farkas).**  It suffices, for every `theta` and every configuration `z` with
> `delta(z) > 0`, to exhibit `w >= 0` on rows `z` satisfies with `sum w a = 0`, `sum w = 1` and
> `sum w b >= 0`.  **[proved]**, §1.1.

> **S2 (chain counting, A4 corrected).**  If more than `(T-1)^2 = 4` of the six squares have tilt
> `< eps` with `eps < 1/(T-1) = 1/2` rad, then the x- or y-DAG on those squares has a path of
> three.  **[proved]** elsewhere: `proof-architecture.md` §0a item 1.

> **Lemma H.**  §1.3 above.  **[proved]**, with the `L >= L*(T, t)` criterion.

> **Lemma W.**  §1.4 above.  **[proved]**.

> **Lemma Z (axis-parallel closure).**  If the crossbar's links all use the normals of
> axis-parallel squares, Lemma H at `t = 0` reads
> `delta <= (T - P_{s_1} - P_{s_T} - sum_s m_s)/(T+1) <= 0`, with equality iff the two end squares
> and all link partners are axis-parallel.  **[proved]** (§1.5; the legs contribute nothing at
> `t = 0`, weight `tan t = 0`).

Lemma Z alone closes the whole axis-parallel stratum and, more, every configuration containing a
staircase chain of three axis-parallel squares.  In particular it closes the pinwheel
(`S6_SKELETON.md` §5: the `3x3` tiling minus its main diagonal has **no straight row of three**,
but it has **16** staircase chains of three at level `0`, e.g.
`(0.5,2.5) -> (1.5,2.5) -> (2.5,1.5)`, and `C = CW = H = F = 0` on every one of them;
`runs/t3_chain_pinwheel.txt`).

### 5.2 What is NOT proved, and the exact shape of the hole

The proof **does not close**.  Three statements are needed and none is proved.

> **(O1) Chain existence.**  Every configuration of six unit squares in `[0,3]^2` with
> `delta >= 0` contains a wall-to-wall chain of three at level `delta`.
> **[measured]** `720/720` sampled optima, with no counterexample at any `delta`.  Nothing here
> proves it.  Chain *counting* (S2) gives it only when `>= 5` squares are near-axis; the hole
> `k in {3, 4}` and the far field `k <= 2` are untouched by counting.

> **(O2) Enough transverse structure.**  That chain carries `L >= L*(3, t)` transverse links at a
> common tilt `t`, or has a common separating normal and rise `R <= R*(3, t)`.
> **[measured false]** at `delta < 0` for uniform tilt `>= 28.8 deg` (the maximum `L` available is
> `2`, and `L* = 2.485` at `45 deg`).  At `delta >= 0` it is not refuted — but it is also not
> derivable from anything local, because **five** squares at `45 deg` in `[0,3]^2` have
> `delta = +0.0503 > 0` and do carry a chain of three with `L = 2`.  Any proof of (O2) must
> therefore use the sixth square, and an H with `L = 2` cannot see it.

> **(O3) The cycle.**  Every configuration with `delta >= 0` and no chain satisfying (O2)
> contains a closed cycle of `k` pair separations, with `q` "turns" (vertices where the link
> normal changes), each turn cancelled by one wall row, with `k >= q sin(alpha/2) (T - u)`
> (`alpha` the turn angle).  **[open]**.  §5.3.

**The region left open**, stated as precisely as the measurement allows: the set of angle vectors
at which *no* wall-to-wall chain of three has both a common separating normal and enough
transverse structure.  Sampled, that is

* the coherent large-tilt cone `{ all six tilts of one sign, all >= 28.8 deg }` — H fails, MT
  works (so the rise is small enough, but it cannot be *proved* small by an H);
* a thin far-field set of mixed-sign configurations with tilts in `[20 deg, 45 deg]` — both fail,
  `4 / 720` samples (`generic`, `bigrand20`, `bigrand30`, `hole0_5`);
* the coherent small-tilt cone `{ all six tilts of one sign, `<= 20 deg` }`, where H loses an
  amount of order `t^3`, the same order as `delta*` there — `10 / 72` samples.

In all three the repair is the same and is measured to work `68/68`: **one extra main-direction
link**.

### 5.3 The sub-lemma that would close it

> **Cycle lemma (the statement that is missing).**  Let `Q_1, ..., Q_k` be squares of a packing of
> `[0,T]^2` with margin `delta`, arranged in a cyclic sequence in which consecutive pairs are
> separated, `n_i . (c_{i+1} - c_i) >= m_i + delta` (indices mod `k`), by unit normals `n_i`.
> Suppose that going round the cycle the normal changes at exactly `q` vertices, and that at each
> such vertex `j` the jump `n_i - n_{i-1}` is `+- lambda_j e_x` or `+- lambda_j e_y`
> (`lambda_j = |n_i - n_{i-1}|`), i.e. every turn is axis-aligned and is cancelled by **one** wall
> row of weight `lambda_j`.  Assume the `q` wall rows split evenly between `lo` and `hi`.  Then,
> summing the `k` link rows at weight `1` and the `q` wall rows at weight `lambda_j`,
>
>     delta  <=  [ (1/2) sum_j lambda_j (T - u_j)  -  sum_i m_i ]  /  ( k + sum_j lambda_j )
>
> — a Farkas certificate with **no rise in it**, because the cycle closes.  Equivalently, by
> Lemma W at a common tilt, `omega = (1/2) sum_j lambda_j / (k + sum_j lambda_j)` and
> `delta <= 0` iff `k >= q sin(alpha/2) (T - u)` for a cycle with `q` turns of equal angle
> `alpha`.

*Status.*  The Farkas arithmetic is one line and is **[proved]** in the case verified in §3.2
(`T = 3`, `t = 45 deg`, `k = 5`, `q = 4`, `lambda = sqrt2`, giving `(93 - 66 sqrt2)/7`).  The
general statement above is **[heuristic]**: it is the obvious generalisation and it reproduces
`45 deg` exactly, but it was not checked at other tilts and the "each turn is axis-aligned"
hypothesis is restrictive (it is what forces exactly four turns at `45 deg`).

*What is genuinely open* is the **existence**: that every configuration with `delta >= 0` and no
good chain contains such a cycle.  That is the far-field statement the repo has been missing, and
it is the natural home for the "pinwheel" that both `S6_SKELETON.md` §4.2 and `S6_LOCAL.md` §3
identify as the hard direction — a pinwheel *is* a cycle.

### 5.4 Honest summary of the attempt

At `T = 3` the architecture reduces `s(6) = 3` to **(O1) + (O2 or O3)**, i.e. to two existence
statements about the separation graph of a hypothetical packing.  Both are combinatorial-geometric
statements about six squares with `delta >= 0`, and neither is easier than what Kearney–Shiu
prove; what the architecture buys is that the *inequalities* are now exact, closed-form and
`T`-uniform (§1.3, §1.4, §5.3), so the whole `T`-dependence is quarantined in the existence
clauses.  That is a real gain for `T = 4` planning and it is not a proof of `s(6) = 3`.

---

## 6. The filter

### 6.1 The lemmas, stated `T`-generically

* **Lemma W** (`delta <= omega (T + 2 - u) - 1` at uniform tilt): `T`-generic, **true at every
  `T`** including `T = 11, 12, 17`.  It is an identity, so it cannot fail; `T` enters only as the
  container side in `b = P - T`.
* **Lemma H** / `L >= L*(T, t)`: `T`-generic, **true at every `T`**.  `T` enters as the chain
  length (the `T - 1` links and the `T - u` wall span).  Its *hypothesis* gets **easier** as `T`
  grows: `L*(T, t)/(n - T)` at `45 deg` is `0.828, 0.539, 0.409, 0.332, 0.173, 0.158, 0.110` for
  `T = 3, 4, 5, 6, 11, 12, 17`.  So **this is not the lemma that fails at `T = 11`.**
* **Lemma Z**: `T`-generic, true at every `T`.
* **(O1) chain existence**: `T`-generic in statement; **at least one of (O1), (O2) is false at
  `T = 11`** (§6.2), but which was not measured here.
* **(O2) enough transverse structure / small rise**: `T`-generic; measured false at `T = 3`
  already at `delta < 0` (uniform tilt `>= 28.8 deg`).
* **(O3) cycle existence**: `T`-generic, status unknown at every `T`.

### 6.2 Which one is false at `T = 11`, and where `T` enters as a number

**(O1) or (O2) — and which one was not measured here.**  Cantrell's `11 x 11` packing of `110`
unit squares (2025, unrefereed) verifies at `delta = +2.749e-04` in exact rational arithmetic
(`BANDCUT_SCAN.md`, `proof-architecture.md` §0a items 9, 11).  Lemma H and Lemma W are identities,
so `delta > 0` forces one of the two *existence* clauses to be false there; **which one is exactly
the measurement `notes/review-2026-09-20c.md` asks for, and it was not done here.**  The published
`T >= 12` scheme is a band that cuts every row and column (`proof-architecture.md` §0a items 8,
11), which would make (O1) the failing one, but `T = 11` is *outside* that scheme, so attributing
the failure to (O1) is a **[guess]** until the wall-to-wall chains of that packing are enumerated.
Where `T` enters as a number:

1. **the chain length is `T`**, and a chain of `T` unit squares spanning `[0,T]` has
   `sum_s m_s + P_1 + P_T >= T` with equality only in the tiling — the inequality is critical at
   exactly `n = T^2 - T`;
2. **the band angle** `2 arctan(1/T)` (`28.07 deg` at `T = 4`, `10.30 deg` at `T = 11`): a stack of
   `T` squares at that tilt spans exactly `T`, so the band can be inserted without a chain
   surviving.  This is the `T`-dependence `proof-architecture.md` §0a items 7–8 identifies, and it
   is the one number that gets *small* as `T` grows, which is why the counterexamples start at
   moderate `T`;
3. **`L*(T, 0) = T - 2`**: the H needs `T - 2` transverse links even at zero tilt, i.e. the
   certificate must span `2T - 2` of the `T^2 - T` squares.

The quantity that does **not** carry `T` is the closed-form inequality, and that is the useful
negative result of this section: any amount of work sharpening `H(T, L, t)` is work that cannot
distinguish `T = 4` from `T = 11`.

### 6.3 They must not prove `s(5) = 3`  [measured]

`s(5) = 2 + 1/sqrt2 = 2.7071 < 3`, so five unit squares fit in `[0,3]^2` with room; any lemma that
applied to them would be false.  Measured (`runs/t3_chain_filter.txt`, `search/t3_chain.py`):

| five squares in `[0,3]^2` | `delta*` (measured) | chains of three at that level | best H-restricted dual |
|---|---|---|---|
| `theta = (0,0,0,0,45 deg)` | `+0.085786` | **2** | `+0.085786` |
| `theta = (0,0,0,0,30 deg)` | `+0.077350` | 2 | `+0.077350` |
| `theta = (45 deg)^5` | `+0.050253` | 2 | `+0.050253` |
| `theta = 0` | `+0.000000` | 4 | `+0.000000` |

So **chain existence does not fail for five squares** — the chains are there.  What fails is the
common-tilt clause: with only four axis-parallel squares and `(T-1)^2 = 4` level cells, chain
counting gives no axis-parallel chain of three (A4 is exactly critical at `k = (T-1)^2`), so
**every** chain of three must run through the tilted square, and its two links each pay
`m - 1 = (W(45 deg) - 1)/2 = (sqrt2 - 1)/2 = 0.2071`.  Lemma H at `t = 0` then returns

    (T - P_{s_1} - P_{s_T} - sum_s m_s) / (T+1)  =  (3 - 1/2 - 1/2 - 2 x 1.2071)/4  =  +0.1464  >  0 ,

which is exactly the mixed-tilt first-order penalty of `notes/bandcut-cost.md` §1 doing its job.
The filter is passed, and it is passed by the *`|D|/2` per link* term, not by any counting.
(Note that the sixth square is what makes `5 -> 6` work: with six squares, `6 > (T-1)^2 = 4`
forces an axis-parallel chain of three as soon as five of them are near-axis.)

### 6.4 One more filter: `T = 2`

`L*(2, t) = 0` at `t = 0` and `0.657` at `45 deg`, `L_max = T(T-2) = 0`.  So Lemma H proves
`s(2) = 2` only for `t = 0` and says nothing at `45 deg` — correct, since `s(2) = 2` needs the
`45 deg` case, which is Friedman's centre-point lemma.  The lemma set is silent where it should
be.

---

## 7. Reproduce

    # the angle-space scan (720 samples, 321 s on 14 processes)
    python3 search/t3_chain.py scan --reps 12 --limit 300 --extra 600 --nproc 14 \
        --out runs/t3_chain_scan1.jsonl
    python3 search/t3_chain.py report runs/t3_chain_scan1.jsonl

    # the closed forms and their verification (instant)
    python3 search/t3_chain_hform.py                      > runs/t3_chain_hform.txt

    # what shape the true dual has
    python3 search/t3_chain_shape.py  runs/t3_chain_scan1.jsonl > runs/t3_chain_shape1.txt

    # is a long enough H available in the configuration at all  (~2 min)
    python3 search/t3_chain_exist.py  runs/t3_chain_scan1.jsonl > runs/t3_chain_exist1.txt

    # the minimal repair of each H failure  (~3 min)
    python3 search/t3_chain_repair.py runs/t3_chain_scan1.jsonl > runs/t3_chain_repair1.txt

    # MT with a genuinely common separating normal  (~8 min)
    python3 search/t3_chain_mt.py     runs/t3_chain_scan1.jsonl > runs/t3_chain_mt1.txt

`runs/t3_chain_filter.txt` and `runs/t3_chain_counterex.txt` are the two ad-hoc checks quoted in
§2.3, §3.2 and §6.3 (the `L*` crossover, the five-square table, the `45 deg` configuration and its
link normals); the commands are at the top of each file's history in this note.

Library entry points: `t3_chain.all_rows(z, n, T, lvl)` (the row system a configuration
satisfies), `t3_chain.dual_bound(A, b, keep)` (the restricted dual LP),
`t3_chain.h_support(tags, TH, ax, path, mode)` (`mode in {C, CW, H, HP, F}`),
`t3_chain.best_h(...)`, `t3_chain.delta_star(theta, T, rng)`,
`t3_chain_hform.H(T, L, t)` / `Lstar(T, t)` / `omega_H(T, L, t)`,
`t3_chain_mt.best_mt(z, n, T, lvl)`.  `s6skel`, `s6local` and `bandcut_k` are imported, never
modified.
