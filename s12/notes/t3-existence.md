# t3-existence: the existence clause at `T = 3` — its real scope, what is now proved, and what is left

*2026-09-21.  Task `tasks/t3-existence/README.md`.  Code (new; nothing in `search/` modified —
`t3_chain`, `s6skel`, `s6local`, `bandcut_k` are imported):
`search/t3_exist_w.py` (Lemma W without the uniform-tilt hypothesis; the scope audit),
`search/t3_exist_master.py` (the wall-free form of the dual; the corrected cycle lemma),
`search/t3_exist_pigeon.py` (the centre-pigeonhole and its `T`-filter),
`search/t3_exist_scope.py` (how tilted a real packing can be),
`search/t3_exist_row.py` (the stress test of the surviving clause, many optima per angle vector),
`search/t3_exist_zero.py` (the axial-chain detector `zchain`, used by the two above; its
`walk` command over the zero set is superseded by `t3_exist_scope.py` and no number here comes
from it).
Runs: `runs/t3_exist_w1.txt`, `runs/t3_exist_pigeon.txt`, `runs/t3_exist_master1.txt`,
`runs/t3_exist_scope1.jsonl` (+`.log`, `runs/t3_exist_scope1_report.txt`),
`runs/t3_exist_row1.jsonl` (+`.txt`), `runs/t3_exist_hunt1.jsonl` (+`.txt`).
Labels: **[proved]** = proved here by hand and re-checked numerically by a script here;
**[measured]** = produced by a multistart or by reading an LP dual, hence one-sided;
**[heuristic]**; **[guess]**.*

---

## 0. Verdict, up front

> **(E3) is open, and the obstruction is not the one the repo has been working on.  Every
> certificate failure on record at `T = 3` and at `T = 4` — the `68/720` H-failures of
> `t3-chain.md` §2.2, the `37/727` dichotomy failures of `T4_CYCLES.md` §7 — is at a configuration
> with `delta < 0`, i.e. at a configuration that is not a packing and that (E3) does not quantify
> over.  Inside (E3)'s actual scope the H is not needed either: at every margin-`0` packing found
> here the *bare chain of three* already certifies `delta <= 0` exactly.**
>
> **The irreducible obstruction is that the surviving existence clause is a statement about exact
> incidences — "some three of the six squares are exactly axis-parallel and form a wall-to-wall
> chain" — which has to be separated from configurations of margin `-0.36 t^2` (uniform tilt) and
> `-0.12 t^3` (three exact zeros) as the tilts go to `0`.  Every mechanism in the brief's item 2
> loses that at first order in `t` (§4.4), so X is the rigidity statement `S6_LOCAL.md` §5 asks
> for, not a combinatorial one.  The `T = 4` route inherits it verbatim,
> and inherits it in a worse form: the one mechanism proved here that does dispose of the far
> field at `T = 3` (a centre-pigeonhole against Graham's 6-point constant) is arithmetically empty
> at `T = 4`, by `0.5 %`.**

Six things are new here, in decreasing order of how much they change the plan.

1. **Scope.**  (E3) quantifies over configurations with `delta >= 0`.  Every configuration with
   `delta >= 0` is automatically an optimum of its own angle vector, and the set of angle vectors
   it can have is exactly `Z = {theta : delta*(theta) = 0}`.  All `68` H-failures of
   `t3-chain.md` have `delta <= -8.8e-06`; the `34` in-scope samples of the same scan have
   `H <= 0`, `CW <= 0` **and** `C <= 0` — `34/34`, and `C` is the *bare chain*, four rows.
   **[measured]**, §1.  So the `45 deg` five-square identity, the cycle, the pentagon, the far-field
   thicket of `T4_CYCLES.md` §2.1 — all of it lives outside (E3).
2. **Lemma W holds without the uniform-tilt hypothesis** (§2.1): for *every* normalised Farkas
   certificate, `delta <= Omega (T+1)/2 - 1 - B`, where `Omega` is the total weight on wall rows
   and `B >= 0` is an explicit tilt bonus.  **[proved]**, verified on all `720` duals of
   `runs/t3_chain_scan1.jsonl` to `1.3e-07` (the LP tolerance).
3. **The walls can be eliminated from the architecture entirely** (§2.2).  Every wall row is
   satisfied by every packing, so a certificate is determined by its *pair-row* weights alone, and
   > `delta <= [ sum_i |r_i|_1 (T - u_i)/2 - sum_e w_e m_e ] / [ sum_e w_e + sum_i |r_i|_1 ]`,
   `r_i` = the residual of the link combination at square `i`.  Minimising the right-hand side
   reproduces the full LP dual **exactly** (checked at `250` optima, max discrepancy `1.0e-07`).
   The whole certificate question at `T = 3` is then one line: *is there a nonnegative weighting of
   the separating links with `sum_i |r_i|_1 (3 - u_i)/2 <= sum_e w_e m_e`?*  **[proved]**
4. **That formula repairs `t3-chain.md` §5.3's cycle lemma**, whose "every turn is axis-aligned"
   hypothesis `T4_CYCLES.md` §6 showed to be load-bearing and nearly empty: the correct turn cost
   is the **`L1`** norm of the normal jump, not `2 sin(alpha/2) = ` its `L2` norm, and with `L1`
   the hypothesis can be dropped.  At a common tilt `t` a `q`-turn `k`-cycle certifies
   `delta <= 0` iff `k m >= q C (T - u)`, `C = cos t`; at `t = 45 deg` this is the old formula
   (there `2C = sqrt2`), everywhere else it is strictly weaker and *valid*.  **[proved]**, §2.3.
   With the correct constant the **four-turn** cycle of `t3-chain.md` §3.3 — "the far-field
   mechanism the repo is missing", the pinwheel dual — needs `t >= 38.41 deg` for a pentagon and
   `t >= 25.73 deg` even for a hexagon on all six squares, while item 5 below shows there is no
   packing above `25.8431 deg`: **that branch of the architecture has an empty intersection with
   (E3)'s scope at `T = 3`, to within `0.12 deg`.**
5. **A proved far-field statement, and the `T`-filter that goes with it** (§3).  The centres of a
   packing are `n` points at pairwise distance `>= 1 + delta` in a square of side `T - u_min`.
   With Graham's `d_6 = sqrt(13)/6`:
   > **[proved]** If every one of six unit squares in `[0,3]^2` has tilt `>= 25.8431 deg`, no
   > packing exists — *whatever* the margin.
   Lemma H fails above `28.7959 deg`.  **The two overlap**, so the uniform-tilt line is completely
   covered and the "far-field mechanism is missing" problem of `t3-chain.md` §0(4) is void at
   `T = 3`.  At `T = 4` the same computation needs `4 - sqrt2 = 2.5858 >= 1/d_12 = 2.5725`, which
   is **true at every tilt**: the mechanism is empty there, by `0.5 %`.  It also passes both
   filters: silent at `n = 5` (`(3-u)/sqrt2 >= 1.1213 > 1` always) and correctly silent at `T = 2`,
   `t = 0`.
6. **The existence clause that the data actually supports is much smaller than an H** (§4): a
   wall-to-wall chain of three whose two link normals are *exactly* axial.  It has a two-line
   identity (Lemma Z'', §2.4), it is `T`-uniform, it needs no legs, no rung and no cycle, and it
   holds at every margin-`0` packing found here.  What it needs is three exactly-axis-parallel
   squares, and that is exactly the statement no combinatorial argument produces.

**What in the brief's premises turned out to be wrong** is collected in §7.

---

## 1. The scope of (E3), and why every recorded failure is outside it

### 1.1 Three elementary observations  [proved]

Write `delta(z)` for the margin of a configuration `z` (the largest `delta` for which every wall
row and one separating row per pair holds), and `delta*(theta) = max_z delta(z)`.

> **(S-a)** A configuration has `delta(z) >= 0` **iff** it is a packing of six closed unit squares
> in `[0,3]^2` (pairwise disjoint interiors, inside the container).
>
> *Proof.*  `delta(z) >= 0` says every wall constraint and, for each pair, some separating-axis
> constraint holds with gap `>= 0`; for convex polygons that is exactly disjointness of interiors
> (`S6_SKELETON.md` §3.1), and the wall rows are exactly containment. ∎

> **(S-b)** If `s(6) = 3` — i.e. `delta*(theta) <= 0` for all `theta` — then every configuration
> in (E3)'s scope has `delta(z) = 0` **and is a global optimum of its own angle vector**.
>
> *Proof.*  `0 <= delta(z) <= delta*(theta) <= 0`. ∎

> **(S-c)** At such a configuration every certificate of value `0` is supported on **tight** rows.
>
> *Proof.*  Complementary slackness: if `sum_r w_r (a_r . c - b_r - delta) = 0` with `w >= 0` and
> each bracket `>= 0`, every row in the support has slack `0`. ∎

(S-b) has a consequence the brief's item 3 does not anticipate: **"test the candidate lemmas at
configurations that are not LP optima" is vacuous.**  There are no non-optimal configurations with
`delta >= 0`.  What *is* under-tested is the opposite thing — the **variety** of optima at a single
`theta in Z`, since at `theta = 0` the optimal set is the whole (positive-dimensional) family of
packings of six axis-parallel unit squares in `[0,3]^2`.  §4 tests that instead.

### 1.2 The audit  [measured]

`python3 search/t3_exist_w.py scope` over `runs/t3_chain_scan1.jsonl` (`runs/t3_exist_w1.txt`):

| | |
|---|---|
| samples | `720` |
| in (E3)'s scope (`delta >= 0`) | `34` |
| samples with `H > 0` | `68` |
| samples with `H > 0` **and** `delta >= 0` | **`0`** |
| largest `delta` among the `68` H-failures | `-8.83e-06` |

and at the `34` in-scope samples (all with `delta = 0`, all in the `Z3..Z6` families, i.e. three to
six angles exactly `0` and the rest free up to `43.07 deg`):

| bound | support | holds at |
|---|---|---|
| `C` (the chain's two links + its two end wall rows — four rows) | minimal | `34/34`, value exactly `0` |
| `CW` (chain + every wall row) | | `34/34` |
| `H` (chain + walls + every transverse pair row) | | `34/34` |
| the axial chain of three of §4 (`t3_exist_zero.zchain`) | | exists at `34/34`, value `0` |

(`CW` and `H` are columns of `runs/t3_chain_scan1.jsonl`; `C` and the axial chain were recomputed
here, since the scan does not record them.)

At `T = 4` the same reading applies without a new measurement: `T4_CYCLES.md` §7's tables are
taken over the `727` sampled optima **with `delta* < 0`**, and the `37` dichotomy failures are
`37` of those `727`, so every one of them has `delta < 0`.  Hence: **`T4_CYCLES.md` §7's counterexamples to
"chain-or-cycle" are counterexamples to a statement about `delta < 0` configurations, which no
existence clause asserts.**

*What this does and does not buy.*  It does not make (E3) easier to prove; it makes the **target**
smaller and different.  The `45 deg` regime that `t3-chain.md` §0(3–4) identifies as the place
where the architecture needs a cycle is empty of packings (§3 proves it), so a proof of (E3) never
has to produce a cycle there — it has to produce the *emptiness*, which is a different job.

---

## 2. The identities, restated without their hypotheses

### 2.1 Lemma W without uniform tilt  [proved]

> **Lemma W'.**  Let `w >= 0` be supported on rows of the wall-carrying LP with
> `sum_r w_r a_r = 0` and `sum_r w_r = 1`.  Put
> `omega_x = w(lo-x rows)`, `omega_y = w(lo-y rows)`, `Omega = 2(omega_x + omega_y)`, and
> `B = sum_{walls} w_r (u_{i(r)} - 1)/2 + sum_{pairs} w_r (W(D_{ij}) - 1)/2 >= 0`.
> Then `w(lo-x) = w(hi-x)` and `w(lo-y) = w(hi-y)`, and
>
>     delta  <=  Omega (T+1)/2  -  1  -  B .
>
> In particular `delta <= 0` as soon as `Omega <= 2(1 + B)/(T+1)`.

*Proof.*  Summing the `x`-coefficient equations over all squares annihilates every pair row (its
two coefficients are opposite) and leaves `w(lo-x) - w(hi-x) = 0`; same in `y`.  Hence the wall
weight splits as `omega_x, omega_x, omega_y, omega_y` and
`sum_r w_r b_r = [Omega/2 + B_wall] - T Omega/2 + [(1 - Omega) + B_pair] = 1 - Omega(T+1)/2 + B`,
using `P_i = 1/2 + (u_i - 1)/2` and `m_{ij} = 1 + (W - 1)/2`.  The Farkas bound is
`delta <= -sum_r w_r b_r`. ∎

`t3-chain.md` §1.4 states this only at a common tilt, where `B = Omega(u-1)/2` and the bound reads
`omega(T + 2 - u) - 1` with `omega = Omega/2`; the two agree.  Verified on all `720` duals of the
`T = 3` scan: worst `|` identity `-` direct `sum w b| = 1.3e-07`, worst `|w(lo-ax) - w(hi-ax)| =
7.0e-08` (both at the LP's own tolerance), `runs/t3_exist_w1.txt`.

### 2.2 The wall-free form of the dual  [proved]

Every wall row is satisfied by every packing (`lo-x(i)` reads `x_i >= P_i + delta`, which is
containment).  So the wall rows are *free* — they are available at every square of every
configuration in (E3)'s scope — and a certificate is determined by its pair-row weights.

> **Master formula.**  Let each satisfied pair row `e` be read as a directed link
> `tail(e) -> head(e)` with unit normal `n_e` and rhs `m_e`, and give it weight `w_e >= 0`.  Put
>
>     r_i  =  sum_{head(e) = i} w_e n_e  -  sum_{tail(e) = i} w_e n_e        (the residual at i).
>
> Cancelling `r_i` with wall rows at `i` costs exactly `|r_i|_1` of weight (walls carry
> `+-e_x, +-e_y`), and this completion is optimal.  The resulting bound is
>
>     delta  <=  [ sum_i |r_i|_1 (T - u_i)/2  -  sum_e w_e m_e ]
>                / [ sum_e w_e  +  sum_i |r_i|_1 ]  ,
>
> and minimising over `w >= 0` reproduces the full dual at `z` whenever that value is below
> `(T - u_max)/2` (`>= 0.79` at `T = 3`), which is every case of interest.  In particular
>
>     delta <= 0     iff     sum_i |r_i|_1 (T - u_i)/2  <=  sum_e w_e m_e .

*Proof.*  The residual at `i` must be cancelled by rows whose coefficient vector is supported on
`c_i`, i.e. by the four wall rows at `i`, which carry `+e_x` (`lo-x`, `b = P_i`), `-e_x` (`hi-x`,
`b = P_i - T`) and the two `y` ones.  Cancelling a residual with components `(r_x, r_y)` costs
`|r_x| + |r_y|` of weight and no less.  Adding matched weight `lambda` to a `lo/hi` pair on the
same axis at square `i` cancels nothing and moves the bound from `N/D` to
`(N + lambda (T - u_i))/(D + 2 lambda)`, which is increasing in `lambda` exactly when
`N/D < (T - u_i)/2`.  So **the `L1` completion is optimal among certificates of value below
`(T - u_max)/2`**, which at `T = 3` is `>= (3 - sqrt2)/2 = 0.7929` — i.e. among all certificates
that say anything at all.  Writing
`L = sum_i |r_i|_1`: since `sum_i r_i = 0`, the `lo` and `hi` weights on each axis are equal, so
the wall part of `sum w b` is `sum_i |r_i|_1 P_i - T L/2 = sum_i |r_i|_1 (u_i - T)/2`, and the
pair part is `sum_e w_e m_e`.  Divide by `sum w = sum_e w_e + L`. ∎

Checked against the LP: at `250` sampled optima the master formula evaluated on the *pair-row part
of the LP's own dual* reproduces the LP value, worst discrepancy `1.0e-07`
(`python3 search/t3_exist_master.py verify`, `runs/t3_exist_master1.txt`).

Two things this buys.

* **The wall rows disappear from the architecture.**  `C`, `CW`, `H`, `HP`, `H1`, `CYC` differ only
  in which *links* they are allowed; the walls are never a modelling choice.  `t3-chain.md` §1.2's
  "the four closing wall rows are not optional" is the statement that `|r_i|_1 != 0` at a leg's
  foot.
* **The criterion is one inequality between two sums over the separation graph,** with `T`
  appearing once, linearly.  A chain of `k` links with a common normal `n` has `L = 2|n|_1 = 2u(n)`
  and certifies iff `k m >= u(n)(T - u)`, which at `u = 1` is `k >= T - 1`: *exactly* the
  wall-to-wall chain, at every `T`.  That is the cleanest statement of why `T = 3` is critical —
  the chain of three is an equality case of the criterion, and nothing weaker survives.

### 2.3 The cycle lemma, corrected  [proved]

`t3-chain.md` §5.3 assumes every turn of the cycle is axis-aligned and then charges `2 sin(alpha/2)`
per turn; `T4_CYCLES.md` §6 showed the hypothesis is load-bearing (dropping it gives `-0.1399` at a
`delta = 0` configuration, i.e. an invalid bound) and nearly empty (at a common tilt `t` the jump
between a square's own two edge normals is axis-aligned only at `t = 45 deg`).  The master formula
supplies the repair with no hypothesis at all: the turn at a vertex with incoming normal `a` and
outgoing normal `b`, both at weight `1`, has residual `a - b` and costs `|a - b|_1`, not
`|a - b|_2 = 2 sin(alpha/2)`.

> **Cycle lemma (corrected).**  A cycle of `k` links at weight `1` whose normal changes at `q`
> vertices, with jumps `j_1, ..., j_q`, certifies
> `delta <= [ sum_v |j_v|_1 (T - u_v)/2 - sum_e m_e ] / (k + sum_v |j_v|_1)`, hence `delta <= 0`
> iff `sum_v |j_v|_1 (T - u_v)/2 <= sum_e m_e`.  At a common tilt `t <= 45 deg` with all turns
> between the square's own two edge normals, `|j|_1 = 2 cos t` and the criterion is
> `k m >= q cos t (T - u)`.

`|a - b|_1 >= |a - b|_2` always, so the corrected lemma is never stronger than the old one, and it
coincides with it exactly when the jump is axis-aligned — `2 cos 45 deg = sqrt2`.  Numbers:
`python3 search/t3_exist_master.py cycle`.  At `T = 3`, `q = 4`, `t = 45 deg` it still reads
`k >= 4.485`, i.e. the pentagon; at `t = 30 deg` it reads `k >= 5.66`, i.e. `k >= 6` — a
Hamiltonian cycle on all six squares — where the old formula said `4.62`, i.e. `k >= 5`.

The `q`-dependence is what this costs the architecture at `T = 3`.  At a common tilt the criterion
is `k >= q cos t (3 - u)`, so with six squares (`k <= 6`):

| turns `q` | pentagon `k = 5` needs | hexagon `k = 6` needs |
|---|---|---|
| `4` (the `t3-chain.md` §3.3 picture, one wall row per wall) | `t >= 38.4059 deg` | `t >= 25.7266 deg` |
| `2` (`T4_CYCLES.md` §0(3)'s shape at `T = 4`) | every `t` | every `t` |

So the **four-turn** cycle — the one `t3-chain.md` §3.3 identifies as "the far-field mechanism the
repo is missing", the pinwheel dual — exists only above `25.73 deg`, and §3 shows there is no
packing above `25.8431 deg`: **the four-turn cycle branch has an empty intersection with (E3)'s
scope, to within `0.12 deg`.**  Two-turn cycles survive the arithmetic at every tilt; none was
needed at any packing measured in §6.  Since §3 below shows there is no packing at all above `25.8431 deg`, this costs nothing.

### 2.4 Lemma Z, corrected and strengthened  [proved]

`t3-chain.md` §5.1's Lemma Z asks that "the crossbar's links all use the normals of axis-parallel
squares" and concludes `delta <= (T - P_1 - P_T - sum m_s)/(T+1) <= 0` *with equality iff the two
end squares and all link partners are axis-parallel*.  The master formula shows that the tilts of
the chain's own squares are not a hypothesis at all — they only help.

> **Lemma Z''.**  Let `s_1 -> ... -> s_T` be squares of a packing with, for each `s`, an exactly
> axial separation `x_{s+1} - x_s >= m_s + delta` (equivalently: a pair row with normal exactly
> `e_x`; it is available whenever one endpoint of the link is exactly axis-parallel).  Then
>
>     delta  <=  [ T - (u_{s_1} + u_{s_T})/2 - sum_s m_s ] / (T + 1)  <=  0 ,
>
> with equality iff `u_{s_1} = u_{s_T} = 1` and every `m_s = 1`, i.e. iff all `T` squares are
> axis-parallel.  Tilting any square of the chain makes it strict: an end square at tilt `t` gains
> `(u - 1)/2`, an interior tilt gains `(W - 1)/2` on each of its two links.

*Proof.*  In the master formula `r_{s_1} = -e_x`, `r_{s_T} = +e_x`, all interior residuals vanish,
so `sum_i |r_i|_1 (T - u_i)/2 = (T - u_{s_1})/2 + (T - u_{s_T})/2` and `sum_e w_e m_e = sum_s m_s`;
`sum w = (T-1) + 2`.  Since every `m_s >= 1` and every `u <= sqrt2`, the numerator is at most
`T - 1 - (T - 1) = 0`. ∎

At `T = 3` that is `delta <= [3 - (u_1 + u_3)/2 - m_{12} - m_{23}]/4`, and e.g. a chain whose middle
square is at `45 deg` and whose ends are axis-parallel gives `(3 - 1 - 1.2071 - 1.2071)/4 =
-0.1036`.  **This is the whole of the certificate side of (E3) at `T = 3`** — see §4.

---

## 3. A proved far-field statement: the pigeonhole on centres

### 3.1 The mechanism  [proved]

> **Proposition P.**  Let `n` unit squares be packed in `[0,T]^2` with margin `delta >= 0`, and let
> `u_min = min_i (|cos th_i| + |sin th_i|)`.  Let `d_n` be the largest possible minimum pairwise
> distance of `n` points in the unit square.  Then
>
>     (T - u_min - 2 delta) d_n  >=  1 + delta .
>
> *Proof.*  Containment gives `c_i in [u_i/2 + delta, T - u_i/2 - delta]^2`, and every one of those
> boxes is contained in the largest of them, so all `n` centres lie in the closed square
> `[u_min/2 + delta, T - u_min/2 - delta]^2` of side `S = T - u_min - 2delta`.
> Disjointness of squares `i, j` gives a unit normal `n` with
> `n . (c_j - c_i) >= m_ij + delta >= 1 + delta`, hence `|c_j - c_i|_2 >= 1 + delta`.  Rescaling the
> square of side `S` to the unit square, `n` points at pairwise distance `>= (1+delta)/S` exist, so
> `(1+delta)/S <= d_n`. ∎

At `T = 3`, `n = 6`, `delta = 0`: `d_6 = sqrt(13)/6 = 0.600925213` (Graham 1963; the optimum for
`n <= 9` is proved), so a packing forces `3 - u_min >= 6/sqrt(13) = 1.6641006`, i.e.
`u_min <= 1.3358994`, i.e. **some** square has `u <= 1.3358994`:

> **Corollary P3.**  **[proved]**  If every one of six unit squares has tilt `>= 25.8431077 deg`,
> they do not fit in `[0,3]^2` at any margin `>= 0`.  (`u = sqrt2 sin(t + 45 deg)`, so the
> condition `u_min <= 1.3359` is `tilt <= 25.8431 deg` for at least one square.)

`runs/t3_exist_pigeon.txt` tabulates the slack `1 - S d_6` along the uniform-tilt line: `-0.2019`
at `0 deg`, `-0.0326` at `20 deg`, `-0.0042` at `25 deg`, `+0.0001` at `25.8431 deg`, `+0.0133` at
`28.7959 deg`, `+0.0471` at `45 deg`.

### 3.2 What it settles, and what it costs

* **Lemma H fails exactly above `28.7959 deg` (`t3-chain.md` §2.3, `L*(3,t) = 2` there).  Corollary
  P3 kills everything above `25.8431 deg`.  The two regimes overlap**, so along the uniform-tilt
  line the existence clause has no far field left: either some square has tilt `< 25.8431 deg`, or
  there is nothing to prove.  The `(93 - 66 sqrt2)/7` pentagon, the "five-square certificate"
  identity and `t3-chain.md` §0(4)'s "the far-field mechanism that is missing is a CYCLE" all
  concern configurations that do not exist.
* It is a statement about `u_min`, not about the average: **one** near-axis square switches it off.
  So it does *not* cover the mixed regime (say four squares near-axis and two at `45 deg`), which
  §4 shows is where the real difficulty sits.

### 3.3 Filters  [proved]

| filter | what P gives | verdict |
|---|---|---|
| `s(5) = 2 + 1/sqrt2 < 3` | `S d_5 = (3 - u)/sqrt2 >= (3 - sqrt2)/sqrt2 = 1.1213 > 1` at every tilt | never fires at `n = 5`: **passed** |
| `T = 2`, `n = 2` (`s(2) = 2`) | `S d_2 = (2-u) sqrt2 >= 1` iff tilt `<= 21.0943 deg` | fires only for two squares both tilted `>= 21.09 deg`, which indeed do not fit in `[0,2]^2`; silent at tilt `0`, where they do: **passed** |
| `T = 4`, `n = 12` | needs `4 - u_min >= 1/d_12 = 2.572479`; `4 - sqrt2 = 2.585786 > 2.572479` | **never fires, at any angle — by `0.5 %`** |
| `T = 11`, `n = 110` | `11 - sqrt2 = 9.5858` against `1/d_110 ~ 9.734` (best known packing, so this is *not* a proof) | would fire above `18.5 deg` *if* `d_110` were known optimal; unproved |

The `T = 4` row is the sharpest thing in this section.  The one mechanism that disposes of the
`T = 3` far field outright is arithmetically empty at `T = 4`, and it is empty by half a percent:
twelve points at pairwise distance `1` fit in a square of side `2.5858` with `0.5 %` to spare.
**Whatever proves the `T = 4` far field, it is not this.**

---

## 4. The existence clause: what it really is at `T = 3`

### 4.1 The clause the data supports

At **every** one of the `749` margin-`0` packings examined here (§6.1) the following holds, and
nothing weaker is needed:

> **(E3-row)**  Every packing of six unit squares in `[0,3]^2` with margin `delta >= 0` contains
> squares `s_1, s_2, s_3` and an axis `ax` with both consecutive pairs separated in the **exact**
> `ax` direction at level `delta`.

By Lemma Z'' (§2.4) that gives `delta <= [3 - (u_{s_1} + u_{s_3})/2 - m_{12} - m_{23}]/4 <= 0`.
(E3-row) is clause (a) of the brief with `t = 0` and `L = 0` — **no legs, no rung, no cycle, and
no `L >= L*(3,t)` condition at all.**  It is `T`-uniform (a wall-to-wall chain of `T` with exact
axial normals certifies at every `T`) and it implies `s(6) = 3`.

An exactly axial pair row is available whenever one endpoint of the link is exactly axis-parallel;
so (E3-row) says, in the worst case, *three of the six squares are exactly axis-parallel and
consecutively separated along one axis*.

### 4.2 What is proved

> **Theorem AP.**  **[proved]**  Let six **axis-parallel** unit squares be packed in `[0,3]^2` with
> margin `delta`.  Then `delta <= 0`, and at `delta = 0` the separation graph contains a
> wall-to-wall chain of three (staircases allowed), so (E3-row) holds.
>
> *Proof.*  Centres lie in `[1/2 + delta, 5/2 - delta]^2` and every pair has
> `max(|dx|, |dy|) >= 1 + delta`.  Let `i ->_x j` mean `x_j - x_i >= 1 + delta`, and likewise in
> `y`; every pair is comparable in at least one, and both relations are subrelations of the strict
> coordinate orders, hence acyclic.  Let `f(i)`, `g(i)` be the lengths of the longest `x`- and
> `y`-paths ending at `i`.  An edge raises `f` (resp. `g`) strictly, so `(f, g)` is **injective**,
> giving `6 <= L_x L_y`.  A path of `L` squares needs `(L-1)(1 + delta) <= 2 - 2 delta`, so
> `L <= 3`, and `L <= 2` as soon as `delta > 0`.  If `delta > 0` then `6 <= 4`, absurd; so
> `delta <= 0`.  At `delta = 0`, `L_x L_y >= 6` with `L_x, L_y <= 3` forces `L_x = 3` or
> `L_y = 3`. ∎

> **Corollary AP5.**  **[proved]**  If **five** of the six squares are exactly axis-parallel,
> (E3-row) holds (apply the injectivity to those five: `5 <= L_x L_y` with `L <= 3` forces some
> `L >= 3`), and `delta <= 0` follows from Lemma Z''.  The same argument applied to `k`
> axis-parallel squares alone gives `delta <= 0` for `k >= 5` and is **exactly tight at
> `k = 4`**: four points pairwise `L∞`-separated by `1 + delta` in `[1/2+delta, 5/2-delta]^2`
> need `2 - 2delta >= 1 + delta`, i.e. `delta <= 1/3`, which is `S6_SKELETON.md` §3.2's measured
> `delta*(4 axis-parallel) = +0.333333333`.  So the counting is sharp at every `k` that matters:
> `+1/3` at `k = 4`, `0` at `k = 5, 6`.

> **Corollary P3** (§3).  **[proved]**  If all six tilts are `>= 25.8431 deg` there is no packing,
> so (E3) is vacuously true there.

> **Proposition NA.**  **[proved]**  If all six tilts are `< arcsin(1/(T-1)) = 30 deg`, both DAGs
> (typed by which normal separates) are still acyclic — an `x`-type normal at angle `a` with
> `n.(c_j - c_i) >= 1` and `|dy| <= 2` forces `dx >= (1 - 2 sin a)/cos a > 0` — and every pair is
> still comparable, so `(f, g)` is injective and a wall-to-wall chain of three **exists**.
> (This is `proof-architecture.md` §0a item 1's A4 with the sharper constant
> `arcsin(1/(T-1))` in place of `1/(T-1)` rad: `30 deg` instead of `28.65 deg` at `T = 3`.)

**Proposition NA does not close anything**, and this is the answer to the brief's item 1(i).  The
chain it produces has link normals at the *tilts* of their owners, and by the master formula the
bare chain's value is

    [ |n_1|_1 (3-u_1)/2 + |n_1 - n_2|_1 (3-u_2)/2 + |n_2|_1 (3-u_3)/2 - m_12 - m_23 ]
    / ( 2 + |n_1|_1 + |n_1-n_2|_1 + |n_2|_1 )

whose first-order expansion in the tilts (`a`, `b` the signed tilts of the two normal owners,
`t_i` the tilts of the chain squares, `D` the tilt differences across the links) is

    [ |a| + |b| + |a - b|  -  (|t_1| + |t_3|)/2  -  (|D_12| + |D_23|)/2 ]  /  4 .

At a **common** tilt `t > 0` this is `+t/4 > 0`, and exactly
`[u(3-u) - 2] / (2 + 2u) = +0.004214` at `t = 1 deg` against `t/4 = 0.004363`
(`python3 search/t3_exist_master.py near`).  So chain counting plus Lemma Z closes the exactly-axis-parallel case
and nothing else; **Lemma Z is discontinuous at tilt `0` and the discontinuity is first order**.

Adding a transverse leg at the chain's end `s_1`, at weight `|a|`, changes the residual at `s_1`
from `-(1, a)` to `(-1 - a t_p, 0)` and creates a residual `|a|(1 + |t_p|)` at the leg's foot: the
`|a|` on the left-hand side moves but does not grow, while the right-hand side gains `|a| m`.  So
**each leg buys exactly the tilt of the normal it cancels, to first order**, and with two legs the
criterion becomes `|a - b| <= (|t_1|+|t_3|)/2 + (sum |D|)/2 + O(t^2)` — at a common tilt,
`0 <= 0`.  That is Lemma H being *exactly critical to first order*, which is why its sign is
decided at second order in the uniform family (`delta*(unif t) = -0.36 t^2 + ...`,
`runs/t3_chain_scan1.jsonl`) and at third order in the coherent cone with three exact zeros
(`S6_LOCAL.md` §2, §5: `-(1/4) t^3`).

**The pinwheel stratum** (brief item 1(iii)) is closed by Theorem AP and nothing else is needed.
The `3 x 3` tiling minus its main diagonal is axis-parallel with `delta = 0`; it has no straight
row of three but `L_x = L_y = 3` in the *staircase* sense, and `t3-chain.md` §5.1 already records
`16` staircase chains of three at level `0` with `C = CW = H = F = 0` on every one.  At a uniform
tilt `t > 0` the pinwheel has `delta = -3(u-1)^2/(u^2+3) < 0` (`S6_LOCAL.md` §3), so it is not a
packing and leaves (E3)'s scope immediately.  What makes the stratum interesting is not the
certificate but that it is the `L_x = L_y = 3` corner of Theorem AP's counting, i.e. the case where
`6 <= L_x L_y` is *least* tight.

### 4.3 The sub-region table

| sub-region of the angle space | status | why |
|---|---|---|
| all six exactly axis-parallel (mod `90 deg`) — includes the pinwheel stratum | **proved** | Theorem AP; `delta <= 0` outright, and the chain exists |
| `>= 5` exactly axis-parallel | **proved** | Corollary AP5 + Lemma Z'' |
| all six tilts `>= 25.8431 deg` | **proved** (vacuous) | Corollary P3: no packing |
| `>= 3` squares with tilt `>= 10 deg` | **open**; measured empty | `0` packings in `768` samples of `runs/t3_exist_scope1.jsonl` (`k = 3,4,5,6`); max `delta* = -3.1e-03` |
| exactly `4` axis-parallel + `2` tilted | **open**; (E3-row) measured to hold at every packing found | counting is *exactly critical*: `4 = (T-1)^2`, and the witness of criticality is the `s(5)` packing (four corner squares, `L_x = L_y = 2`, no chain) |
| six tilts in `(0, 25.84 deg)`, no exact zero (the coherent cone) | **open**; measured empty | `delta*(unif t) ~ -0.36 t^2`; `delta*` on the three-zero cone `~ -t^3/4`.  This is the rigidity case |
| mixed signs, any tilts | **open**; measured easy | every `split*`/`+-` family here has `delta* < 0` by a first-order amount (`|D|/2` per link), so no packing exists and the clause is vacuous |

### 4.4 Why the last two rows cannot be closed by any of the mechanisms in the brief

A configuration of six unit squares all at tilt `t = 0.001 deg`, placed at the tiling positions
minus three cells with the centres re-optimised, is geometrically indistinguishable from a packing:
its margin is `-0.36 t^2 = -1.1e-10`.  It contains **no** exactly axial row, no common-normal chain
of three with legs at the required weights, and no cycle, so (E3) has to exclude it — and the only
thing that excludes it is a quantity of size `1e-10`.

Make that quantitative.  Along the uniform-tilt line the statement to be proved is
`delta*(t) <= 0`, and `delta*(0) = 0` exactly with `delta*(t) = -0.36 t^2 + O(t^3)`: **tight at
`t = 0`, and only quadratically true just off it.**  So an argument must be (i) *exactly* tight in
the axis-parallel case and (ii) accurate to second order in `t`.  The mechanisms of the brief's
item 2 fail one or the other, and it is easy to see where:

* the **pigeonhole** (§3) compares `S d_6` with `1` and its slack at `t = 0` is `-0.2019`, on the
  wrong side by `20 %` — it fails (i) by a wide margin and only recovers above `25.84 deg`;
* **chain counting** (Theorem AP) *is* exactly tight at `t = 0`, but its inequality reverses at
  **first** order in `t`: the engine is "a chain of three needs `2(1 + delta)` of span against
  `2 - 2 delta` available", and at tilt `t` the `x`-step of an `x`-type separation is only
  `(1 - 2 sin t)/cos t = 1 - 2t + O(t^2)` while the box is `3 - u = 2 - t + O(t^2)`, so
  `2(1 - 2t) <= 2 - t` holds with room — the contradiction evaporates at rate `Theta(t)`, ten
  orders of magnitude faster than `delta*` turns negative;
* **unavoidable points** are angle-free by construction (`t3-chain.md` §4.3: the tilt is quantified
  away once, before the combinatorics), so they cannot distinguish `t` from `0` at all;
* **Menger/Dilworth** is chain counting.

> **Obstruction X.**  The existence clause has to separate configurations whose margin is
> `-O(t^2)` (at `T = 3`; `-O(eps^3)` at `T = 4`) from packings, as the tilts go to `0`.  Every
> mechanism that is *stable* under an `O(t)` perturbation of the angles — counting, pigeonhole,
> unavoidable points, Menger/Dilworth — is blind at that scale by construction.  The clause is
> therefore a second-order (at `T = 4`: third-order) **rigidity** statement about the
> tiling-adjacent stratum, exactly the object `S6_LOCAL.md` §5 calls "(ii) a second-order theorem
> on the coherent cone", and the angle-space route does not reduce it: it *is* it.

---

## 5. Three mechanisms, tried

The brief asks for at least three genuinely different attacks on the region that counting does not
reach.  Here they are, with what each proves and where each stops.

### 5.1 Pigeonhole / strip counting on the centres  — **proves a sub-region** [proved]

§3.  The centres are `n` points at pairwise distance `>= 1 + delta` in a square of side
`T - u_min - 2delta`; against Graham's `d_6 = sqrt(13)/6` this kills every configuration whose
*minimum* tilt is `>= 25.8431 deg`.  It is the only mechanism here that closes a positive-measure
region of the angle space outright.

*Where it stops.*  It is governed by `u_min`, so a single near-axis square switches it off, and the
`4 + 2` regime of §4.3 is invisible to it.  Sharpening it would mean replacing the Euclidean
exclusion by the true octagonal one — but for two squares at a **common** tilt the exclusion is the
`L∞` ball of the rotated frame, and then the statement *is* the uniform-tilt case of `s(6) = 3`:

> **The clean sub-problem the uniform-tilt case reduces to** (equivalent, no LP in it):
> *prove that six points with pairwise `L∞` distance `>= 1` cannot lie in a square of side
> `3 - cos t - sin t` rotated by `t`, for any `t in (0 deg, 45 deg]`.*
> At `t = 0` the statement is false with equality (side `2`, the `3 x 3` grid minus three), which
> is exactly why it is delicate.  [measured true; the deficit is `-0.36 t^2`.]

### 5.2 Menger / Dilworth: no chain ⇒ a transverse cut — **is chain counting, and it is exactly critical** [proved]

The cut is not an extra object: it is the other axis.  If the `x`-DAG has no path of three, its
longest path is `2`, so by Mirsky the squares split into two `x`-antichains; an `x`-antichain is a
set pairwise separated in `y`, i.e. a `y`-**chain**.  If the `y`-DAG also has no path of three,
each antichain has `<= 2` elements and `n <= 4`.  This is the `(f, g)`-injectivity proof of
Theorem AP, i.e. `proof-architecture.md` §0a item 1's A4, and it is `T`-uniform: `n <= (T-1)^2`.

*Where it stops.*  Two places, both sharp.

1. It needs the DAGs acyclic, which holds iff every tilt is `< arcsin(1/(T-1))` (`30 deg` at
   `T = 3`, `19.47 deg` at `T = 4`).  Beyond that a "chain" can run backwards and the count is
   void; the mechanism and Corollary P3 leave the band `[30 deg, ...]`-with-a-near-axis-square
   uncovered.
2. `6 > (T-1)^2 = 4` is a slack of **two squares**, and the `4 + 2` regime spends both: four
   axis-parallel squares in a `2 x 2` pattern (the `s(5)` packing's four corner squares) realise
   `L_x = L_y = 2` exactly.  So the counting bound is tight at the only place it is needed, and
   the two tilted squares are what must be priced.  This is the same "one unit of slack" as
   `proof-anatomy.md` §7.2 and as Kearney–Shiu's `6` squares vs `7` points.

### 5.3 Kearney–Shiu unavoidable points as a *source* of incidences — **the one with room left** [heuristic]

Can an unavoidable set *produce* the chain rather than replace it?  Directly, no: "square `i`
contains point `p`" localises `c_i` to within the square's half-diagonal (`<= 0.7071`) of `p`,
while the spacing of their lattice is `0.5858`, so incidences do not imply separations and no chain
falls out.

What their proof actually uses is something the repo's row system **does not contain**, and this is
worth recording as a concrete extension:

> **Line-transversal (chord) rows.**  Fix a line `x = c`.  The squares that meet it cut out
> pairwise disjoint intervals on it, so for any two of them
> `|y_i - y_j| >= (h_i(c) + h_j(c))/2`, where `h_i(c)` is the chord `Q_i ∩ {x = c}`; and the total
> of the chords is `<= T`.  These are **linear** inequalities in the centres at fixed angles, they
> are satisfied by every packing, and they are *not* implied by any single separating-axis row —
> their validity uses the disjunction, not one leaf of it.  `s6skel.chord_len` already computes
> `h`.  A Farkas certificate is allowed to use them (it only needs rows the configuration
> satisfies), so the support-shape list can be extended by them at no cost to soundness.

That is exactly what Kearney–Shiu's contradiction is: their `13/6 - sqrt2 < 2 sqrt2 - 2` is a chord
count on the line `x = 1` (`t3-chain.md` §4.2 reads it as a chain inequality, which is the right
dictionary but the wrong row system — the rows are chord rows, not separation rows).  At tilt `0` a
chord row degenerates to the ordinary `y`-separation row (`h = 1`), so this extension is **invisible
in the axis-parallel stratum and first order in the tilt off it** — i.e. it is a candidate for
exactly the region §4.4 says nothing stable can reach.  Untried here; it is the single most
promising concrete next step, and it is what `notes/proof-anatomy.md` §5.2's reading of the paper
has been pointing at all along.

*Its `T`-filter, however, is bad.*  The unavoidable-set side of Kearney–Shiu dies at `T = 4`
(`t3-chain.md` §4.4: no `11`-point unavoidable set for `[0,4]^2` is known, and `DS7` Table 1 has
none below `15`), and the chord rows themselves are `T`-uniform but their *accounting* — "`n`
squares, `n+1` points, slack `1`" — is not.

---

## 6. The measurements behind (E3-row), and the filters

### 6.1 What was tested  [measured]

Four independent sets of margin-`0` packings, all of them genuine feasible points:

| source | samples | packings found | (E3-row) holds | max value of the Lemma Z'' bound |
|---|---|---|---|---|
| `runs/t3_chain_scan1.jsonl`, the in-scope samples | `720` | `34` | `34/34` (`C = CW = H = 0`) | `0` |
| `runs/t3_exist_scope1.jsonl`, `k` squares at tilt `>= tau`, the rest exactly `0` | `960` | `106`, **all with `k = 2`** | `106/106` | `0` |
| `runs/t3_exist_row1.jsonl`, up to `24` distinct optima at each of `40` angle vectors | — | `387` | `387/387` | `0` |
| `runs/t3_exist_hunt1.jsonl`, the same at `38` angle vectors with only `2` or `3` exact zeros | — | `222`, **all with exactly `3` axis-parallel squares** | `222/222` | `0` |

**`749` margin-`0` packings, `749` axial chains of three, no failure.**  In every one of the
`609` packings of the last two rows the certifying chain consists of three *exactly axis-parallel*
squares, and the number of exactly-axis-parallel squares in a packing was `3` (251 packings), `4`
(246), `5` (97), `6` (15) — **never fewer than three, and never a packing with only two**.

Two structural facts fall out, both one-sided in the reliable direction (every packing reported is
an explicit configuration; "none was found" is not a proof):

* **[measured] `delta* < 0` as soon as three squares carry a common nonzero tilt.**
  `three at t deg` (three angles `t`, three exactly `0`): `delta* = -6.4e-07, -6.8e-05, -4.4e-04,
  -2.3e-03, -6.6e-03` at `t = 1, 5, 10, 20, 45 deg` — third order in `t`, and never `0`.
  Likewise every `k >= 3` cell of `runs/t3_exist_scope1.jsonl` is empty
  (`runs/t3_exist_scope1_report.txt`; `k` squares drawn uniformly from `[tau, 90 - tau]`, the rest
  exactly `0`, `24` samples per cell, max `delta*` shown):

  | `k` \ `tau` | `10` | `15` | `20` | `25` | `30` | `35` | `40` | `45` |
  |---|---|---|---|---|---|---|---|---|
  | `2` | `0` (12/24) | `0` (9) | `0` (8) | `0` (6) | `0` (13) | `0` (13) | `0` (21) | `0` (24) |
  | `3` | `-3.1e-03` | `-1.0e-02` | `-8.6e-03` | `-7.2e-03` | `-6.1e-03` | `-8.4e-03` | `-7.1e-03` | `-6.6e-03` |
  | `4` | `-2.0e-02` | `-2.8e-02` | `-2.5e-02` | `-2.5e-02` | `-2.9e-02` | `-3.0e-02` | `-3.2e-02` | `-3.2e-02` |
  | `5` | `-2.8e-02` | `-3.0e-02` | `-3.2e-02` | `-3.7e-02` | `-4.4e-02` | `-4.2e-02` | `-4.4e-02` | `-4.8e-02` |
  | `6` | `-3.7e-02` | `-3.8e-02` | `-5.0e-02` | `-5.0e-02` | `-5.3e-02` | `-5.1e-02` | `-5.2e-02` | `-4.8e-02` |

  `0` packings in the `768` samples with `k >= 3`; all `106` packings the scan found have `k = 2`,
  and every one of them has `C = CW = H = 0` and an axial chain of three.  The `k = 6`,
  `tau >= 30 deg` cells are the ones Corollary P3 **proves** empty.
* **[measured] but with unequal tilts, three nonzero tilts do occur**: `theta = (0,0,0, 43.07,
  6.70, 0.65) deg` is a packing (`delta* = 0`), as is `(0,0,0, 0.31, 30.91, 33.43) deg`.  So the
  invariant is not "at most two tilted squares"; it is "at least three *exactly* axis-parallel
  squares, and they form a chain".
* **[measured] opposite signs never admit a packing.**  Every `two at +-t deg` family has
  `delta* <= -8.2e-04` — first order in `t`, the `|D|/2`-per-link penalty of `S6_LOCAL.md` §5.1.

### 6.2 Filters  [proved]

> **`s(5) = 2 + 1/sqrt2 < 3`.**  (E3-row) is *false* for five unit squares in `[0,3]^2`, and for
> the sharpest possible reason: any wall-to-wall chain of three forces
> `x_{s_3} - x_{s_1} >= m_{12} + m_{23} + 2 delta >= 2 + 2 delta`, while centres lie in
> `[1/2 + delta, 5/2 - delta]`, giving span `<= 2 - 2 delta`.  So **(E3-row) implies `delta <= 0`
> immediately**, and at the `s(5)` optimum (`delta = +0.0858`) no such chain can exist.  Indeed
> the four axis-parallel squares of the `s(5)` packing sit in a `2 x 2` pattern with
> `L_x = L_y = 2`.  **Passed.**
> **`T = 2`.**  The `T`-generic clause is "a wall-to-wall chain of `T` with exactly axial
> normals".  At `T = 2`, `n = 2`, two axis-parallel squares give the chain and `delta <= 0`
> follows (span `1 - 2delta >= 1 + delta` is false for `delta > 0`); at `45 deg` no axial row
> exists and the clause is **silent**, which is correct — `s(2) = 2` at `45 deg` is Friedman's
> centre-point lemma, not a chain.  **Passed**, and it is silent in exactly the same place
> Lemma H is (`t3-chain.md` §6.4).

### 6.3 The `Z3` stratum, which `t3-chain.md` §2.6 flags as under-sampled

It is the stratum that matters, and it is now sampled: `three at t` (empty at every `t > 0`),
`mixed 43/6.7/0.65` and `mixed 30.9/33.4/0.3` (both packings, `15` and `14` distinct optima each,
(E3-row) at all of them), plus the `k = 3` cells of `runs/t3_exist_scope1.jsonl`.  Nothing in it
threatens (E3-row); what it shows is that the stratum's packings have **exactly three**
axis-parallel squares, so the chain has no slack at all — all three are forced onto it.

---

## 7. What in the brief's premises, and in the repo, turned out to be wrong

* **"Every existence clause is unproved" — right; "so the object to guarantee is an H, an H1 or a
  cycle" — wrong.**  Inside the scope of (E3) the object to guarantee at `T = 3` is a bare
  wall-to-wall chain of three with axial normals: four rows, value exactly `0`, no legs, no rung,
  no cycle, no `L >= L*(3,t)`.  The H's transverse legs, the rung of `t3-chain.md` §2.5 and the
  cycle of §3.2 are all responses to failures at `delta < 0`, where no existence clause applies.
* **Brief item 3 ("test at `delta >= 0` and not only at LP optima; perturb optima; sample margins
  slightly negative") is partly vacuous and partly out of scope.**  Every configuration with
  `delta >= 0` *is* an optimum of its angle vector (§1.1 (S-b)), so there is nothing to perturb
  towards; and a configuration of margin "slightly negative" is not in (E3)'s hypothesis, so a
  lemma that fails there is not refuted.  What the brief should have asked for — and what §6 does
  — is *many distinct optima at one angle vector*, since at `theta in Z` the optimal set is
  positive-dimensional.
* **Brief clause (c) quotes a false lemma.**  `t3-chain.md` §5.3's `k >= q sin(alpha/2)(T - u)` is
  invalid without its "every turn is axis-aligned" hypothesis (`T4_CYCLES.md` §6) and the
  hypothesis is nearly empty.  §2.3 replaces `2 sin(alpha/2) = |jump|_2` by `|jump|_1` and the
  lemma becomes unconditional.
* **Brief item 1(i): Lemma Z does *not* close the `>= 5` near-axis case.**  Lemma Z needs the link
  normals to be *exactly* axial; at a common tilt `t > 0` the bare chain's value is `+t/4 + O(t^2)`
  (§4.2), and the H that repairs it is exactly critical to first order.  What chain counting gives
  is existence of a chain, never its value.  (Chain counting's own constant can be sharpened:
  `tilt < arcsin(1/(T-1)) = 30 deg` at `T = 3`, not `1/(T-1) = 28.65 deg`;
  `proof-architecture.md` §0a item 1.)
* **Brief item 1(ii): "is `omega <= 1/(5-u)` provable from `delta >= 0` alone?"**  At a *uniform*
  tilt `t > 0` the hypothesis `delta >= 0` is contradictory — proved here for `t >= 25.8431 deg`
  (Corollary P3), measured (`delta* = -0.36 t^2`) below it — so the clause is vacuous and
  "provable from `delta >= 0` alone" means "provable *because* the hypothesis is empty", i.e. it is
  the uniform-tilt case of `s(6) = 3` itself.  §5.1 gives it as a clean `L∞`-geometry problem with
  no LP in it.  Separately, Lemma W's uniform-tilt hypothesis is unnecessary (§2.1).
* **`t3-chain.md` §5.1's Lemma Z can be weakened and strengthened at the same time** (§2.4): only
  the *normal owners* need be axis-parallel, and the tilts of the chain's own squares only make the
  bound more negative.
* **Lemma H's closed form is not optimal on its own support at `L = 0`, and its two legs are
  worth exactly nothing at `45 deg`.**  The weight vector of `t3-chain.md` §1.3 puts `1/C` on
  `lo-ax(s_1)` *and* `S tan t` on `hi-ax(s_1)` — opposite walls on one square, which §2.2 shows is
  never optimal (`1/C - S tan t = C`).  Completing the bare chain optimally instead gives, at a
  common tilt,
  > `delta <= [ u (T - u) - (T - 1) ] / [ (T - 1) + 2u ]`,   **[proved]**
  which is `t/4 + O(t^2)` (`0.004214` at `1 deg` against `H(3,0,1 deg) = 0.004512`), and at
  `t = 45 deg`, `T = 3` equals `(3 sqrt2 - 4)/(2 + 2 sqrt2) = (10 - 7 sqrt2)/2 = +0.050252532`
  **exactly** — i.e. the same number as `H(3, 2, 45 deg)`, whose weight vector turns out to be
  twice this one.  So `t3-chain.md` §2.3's "an H-shaped certificate at `T = 3` in the far field is
  literally a five-square certificate" is reading a *three*-square certificate: the chain of three
  plus its four optimal wall rows already gives `+0.050252532`, and the two transverse legs at
  `45 deg` contribute nothing at all.  (At `20 deg` they do: `+0.044334` bare against
  `H(3,2,20 deg) = -0.016891`.)  `python3 search/t3_exist_master.py near`.
* **`t3-chain.md` §0(3)'s headline — "an H-shaped certificate at `T = 3` in the far field is
  literally a five-square certificate, and five squares do fit" — is true and irrelevant to
  existence**: at `theta = (45 deg)^6` there is no packing, by Corollary P3 with `0.047` to spare.
* **Corroboration at `T = 4`, already in the repo and not read this way:** `BANDCUT_K.md` §0(1)
  reports that *every* `delta* = 0` configuration found at `T = 3` **and at `T = 4`** carries a
  wall-to-wall chain of `T` **among its near-axis squares**, certifying through the exact identity.
  That is precisely (E3-row) at `T = 4`, measured and already believed.  The `T = 4` "far-field
  thicket" of `T4_CYCLES.md` §2.1 (`mu >= 2`, median `13` pair rows on `11` squares) is entirely a
  `delta* < 0` phenomenon.

---

## 8. Verdict

> **(E3) is open at `T = 3`.  What is proved here is the whole *certificate* side — Lemma W'
> without its uniform-tilt hypothesis, the wall-free master formula that reproduces the full dual,
> Lemma Z'' with only the normal-owners axis-parallel, and the corrected cycle lemma — together
> with three strata of the *existence* side: the axis-parallel stratum (including the pinwheel),
> the `>= 5`-axis-parallel stratum, and the entire far field, where a centre-pigeonhole against
> Graham's `d_6 = sqrt(13)/6` shows that six unit squares with every tilt `>= 25.8431 deg` do not
> fit in `[0,3]^2` at all.  Inside (E3)'s real scope no H, no rung and no cycle is ever needed: a
> bare wall-to-wall chain of three with exactly axial normals certifies `delta <= 0` at every one
> of the `749` margin-`0` packings tested here.**
>
> **The irreducible obstruction X is that the surviving clause — "three of the six squares are
> exactly axis-parallel and consecutively separated along one axis" — has to separate packings
> from configurations whose margin is `-O(t^2)` as the tilts `t` go to `0` (uniform tilt:
> `-0.36 t^2`; three exact zeros: `-0.12 t^3`).  Every mechanism that survives an `O(t)`
> perturbation of the angles — counting, pigeonhole, unavoidable points, Menger/Dilworth — is
> blind at that scale, so X is a rigidity statement about the tiling-adjacent stratum and not a
> combinatorial one.  `T = 4` inherits X verbatim and one order flatter (`-O(eps^3)`,
> `T4_CYCLES.md` §3, `BANDCUT_K.md` §0(2)), and additionally loses the only far-field mechanism
> that works at `T = 3`: twelve points at pairwise distance `1` fit in a square of side
> `4 - sqrt2 = 2.5858` with `0.5 %` to spare, so the pigeonhole never fires there.**

### 8.1 Sub-regions, proved / open

| # | sub-region | status | mechanism / obstruction |
|---|---|---|---|
| 1 | all six exactly axis-parallel (mod `90 deg`), incl. the pinwheel | **proved** | Theorem AP: `(f,g)`-injectivity, `6 > (T-1)^2`, then Lemma Z'' |
| 2 | `>= 5` exactly axis-parallel | **proved** | Corollary AP5 (`5 > 4`) + Lemma Z'' |
| 3 | all six tilts `>= 25.8431 deg` | **proved**, vacuous | Corollary P3 (Graham `d_6`) |
| 4 | uniform tilt `t in (0, 25.84 deg)` | **open**, measured vacuous (`-0.36 t^2`) | reduces to a clean `L∞` problem, §5.1; nothing proved |
| 5 | `>= 3` squares at tilt `>= 10 deg` | **open**, measured vacuous (`0` packings in `768` samples) | no mechanism; the `k = 3,4,5,6` cells of `runs/t3_exist_scope1.jsonl` |
| 6 | exactly `3` or `4` exactly axis-parallel + tilted rest | **open**, (E3-row) measured at `749/749` | counting is *exactly* critical (`4 = (T-1)^2`, realised by the `s(5)` packing's four corner squares); the two tilted squares must be priced geometrically |
| 7 | six tilts in `(0, 25.84 deg)` with no exact zero (coherent cone) | **open**, measured vacuous | obstruction X: decided at `O(t^2)`/`O(t^3)`; no stable argument reaches it |
| 8 | mixed signs, any tilts | **open**, measured vacuous by a *first-order* amount | `|D|/2` per link; the easiest region, and the only one where a stable argument plausibly works |

*(One caveat on the measurements: they establish (E3-row) at the packings found, all of which have
`delta = 0`.  If `s(6) = 3` were false there would be packings with `delta > 0` too — but the same
multistarts that produced these packings are the ones maximising `delta`, and none ever returned a
positive value.  The three **proved** rows do not have this caveat: Theorem AP and Corollary P3
both handle `delta > 0` directly and derive a contradiction.)*

Rows 1–3 are a proof of (E3) on their union.  Rows 4, 5, 7, 8 are vacuous (no packing exists) and
therefore need a *non-existence* proof, not a support.  Row 6 is the only region where (E3) has
real content and real configurations, and there the clause holds at everything measured.

### 8.2 Can any finite support shape have a provable existence clause?

**Yes for the shape; no for the proof, and the shape was never the problem.**  Inside the scope,
one shape suffices at `T = 3` (the axial chain of three) and the same is already measured at
`T = 4` (`BANDCUT_K.md` §0(1): every `delta* = 0` configuration found carries a wall-to-wall chain
of `T` among its near-axis squares).  The supports available at a configuration are drawn from a
finite list (`24` wall rows and `120` pair rows on six squares, so at most `C(144,13)` vertex
supports), so "a finite family of shapes covers `Z`" is not in doubt; what is in doubt — and what
X says is hard — is the implication `delta >= 0 => the shape is present`.  Sharpening the
identities cannot help: §2.2 shows the criterion is a single inequality
`sum_i |r_i|_1 (T - u_i)/2 <= sum_e w_e m_e` whose equality case is the axial chain of `T` at every
`T`, so there is no slack anywhere in it to be recovered.

The one concrete thing that might: **extend the row system** by the line-transversal (chord) rows
of §5.3.  They are valid, linear, satisfied by every packing, absent from the current system,
first-order in the tilt off the axis-parallel stratum, and they are the rows Kearney–Shiu's proof
actually uses.

### 8.3 Reproduce

    python3 search/t3_exist_w.py verify              > runs/t3_exist_w1.txt     # Lemma W', 720 duals
    python3 search/t3_exist_w.py scope                                          # the scope audit
    python3 search/t3_exist_master.py verify runs/t3_chain_scan1.jsonl 250 \
                                                      > runs/t3_exist_master1.txt
    python3 search/t3_exist_master.py cycle                                     # corrected cycle lemma
    python3 search/t3_exist_master.py near                                      # chain + legs at small tilt
    python3 search/t3_exist_pigeon.py                > runs/t3_exist_pigeon.txt # Corollary P3 + filters
    python3 search/t3_exist_scope.py scan --reps 24 --nproc 2 \
        --out runs/t3_exist_scope1.jsonl             # ~2 h; how tilted a packing can be
    python3 search/t3_exist_scope.py report runs/t3_exist_scope1.jsonl
    python3 search/t3_exist_row.py test --reps 24    > runs/t3_exist_row1.txt   # ~25 min, 387 packings
    python3 search/t3_exist_row.py hunt --reps 40 --nrand 14 \
                                                     > runs/t3_exist_hunt1.txt  # the thin strata

Library entry points: `t3_exist_w.decompose(w, tags, TH, T)` (the `Omega`/`B` split of a dual),
`t3_exist_master.master(wpair, tags, TH)` (the master formula on a link weighting),
`t3_exist_pigeon.tilt_threshold(T, n)`, `t3_exist_zero.zchain(z, tags, TH, lvl)` (the axial chain
of three and its Lemma Z'' value), `t3_exist_row.test_one(theta, rng, reps)`.
