# arch-tledger: where does `s(T^2 - T) = T` break?

*2026-09-20 (night).  Answers §3 of `notes/proof-architecture.md`.  New code: `search/arch_ct.py`
(hole-set search for the uniform-tilt optimum), `search/arch_chain.py` (the wall-to-wall chain
identity, symbolic + configuration analyser).  Runs in gitignored `runs/arch_tledger_*.txt`; every
number quoted here.  Convention of `S6_LOCAL.md`: **`s6local`/`arch_ct` values are feasible points,
i.e. lower bounds on `delta*` and therefore UPPER bounds on `c_T`**; `s6exact` is exact for one
leaf.*

## 0. Verdict

1. **`c_5 = 3/8`, `c_6 = 2/5`, and `c_T` never changes sign.**  Three independent searches agree at
   `T = 5`, and the whole table is explained by one exact inequality (§2), which gives

       c_T = (T - 2) / (2 (T - 1))     for T >= 4,        c_2 = 1,  c_3 = 3/4 .

   That is **increasing** in `T`, with limit `1/2`.  It is `> 0` for every `T >= 3`, at `T = 17` as
   much as at `T = 4`.  **Scenario (α) of §3.3 is refuted**: the tiling does not stop being a strict
   local maximum, so `A2c` is *not* the lemma that carries `T < T*`.
   (`T = 2, 3, 4, 5` are exact on their leaves; `T = 6` extrapolates to `0.3999999` against
   `2/5 = 0.4`, §1.4.)
2. **The mechanism is a single elementary inequality** — the wall-to-wall chain bound (§2) —
   and its `T`-dependence is *entirely* in one integer: the **rise** of the chain.  `delta <= 0`
   holds for every `t` and every `T >= 3` as soon as some wall-to-wall chain of `T` squares rises by
   at most `1`; it fails as soon as every such chain rises by `2`.  Chain rise is an `A4` quantity,
   not an `A2c` quantity.
3. **Part 3 (straight-chain window) verified** and folded into the same identity: the window
   `tan(t/2) >= 1/T` is the special case `R = (T-1) sin t` of the chain bound (§4).
4. The same identity says what *does* die with `T`: the tube in which the local theorem lives has
   radius `~ 2/T`, and the deficit at its edge is `~ 2/T^2`.  The local theorem stays true and
   becomes vacuously shallow.  That is `A5` (the glue), not `A2c`.
5. **Cleemann's packing is scenario (β)** (§5): 199 of its 272 squares are axis-parallel, and its
   bent tilted band **cuts every one of the 17 rows and every one of the 17 columns**, so there is no
   wall-to-wall chain of 17 — `A4` is what fails, exactly as §3.3 guessed.  Its tilt angle
   `arctan(8/15)` is **exactly `2 arctan(1/4)`**, i.e. the `k = 4` straight-chain window of Part 3
   at equality (`4 cos t + sin t = 4` exactly).
6. **The `T`-ledger in §3 is two steps out of date.**  `T = 17` is no longer the smallest known
   counterexample: `s(12^2-12) < 12` (Arslanov–Mustafin–Shangitbayev, EJC 2021, refereed) and
   `s(11^2-11) <= 10.99679327…` (Cantrell, Feb 2025, not refereed).  **The first failing `T*` is in
   `5..11`, not `5..17`.**

---

## 1. `c_T` measured

### 1.1 The instrument: hole sets

At `n = T^2 - T` on the `T x T` tiling, `T` of the `T^2` cells are **empty**; the configuration is
its hole set.  `arch_ct.py holes` enumerates hole sets (all `C(T^2,T)` of them up to the dihedral
group of the grid, or just the `T!` permutation matrices), jitters each by `rho = 0.08`, runs
`s6local.Fixed.ascent`, and — this is the point — classifies the results by the **final** hole set
(round each relaxed centre back to a cell), because the ascent routinely migrates from its start.

| `T` | `n` | hole sets | dihedral classes | family used |
|---|---|---|---|---|
| 3 | 6 | `C(9,3) = 84` | 16 | all (exhaustive) |
| 4 | 12 | `C(16,4) = 1,820` | 252 | all (exhaustive) |
| 5 | 20 | `C(25,5) = 53,130` | 6,814 | all (exhaustive) |
| 6 | 30 | `C(36,6) = 1,947,792` | — | permutations only (`720 -> 115` classes) |

**In every case the optimal final hole set is a permutation matrix** (one empty cell per row and per
column), and the top of the table is a short list of permutation classes.

### 1.2 The numbers

`delta*/t^2` from `runs/arch_tledger_T{3,4,5,6}_*.txt`, and the `t -> 0` extrapolation
(`arch_chain.py fit`, quadratic in `t`):

| `T` | `delta*/t^2` measured | `c_T` extrapolated to `t = 0` | closed form | status |
|---|---|---|---|---|
| 2 | — | `1` | `-(u-1)^2`, exact | **exact** (`S6_LOCAL.md` §3) |
| 3 | `-0.734390` @ `0.8°` | `0.75004` | `-3(u-1)^2/(u^2+3)`, exact on its leaf | **exact on leaf**, `= 3/4` |
| 4 | `-0.33178 … -0.32087` @ `0.4…3.2°` | `0.333337` | trig-rational, exact on its leaf, `-(1/3)t^2 + (2/9)t^3 + (1/54)t^4` | **exact on leaf**, `= 1/3` |
| 5 | `-0.374727 … -0.365505` @ `0.1…3.2°` | **`0.3749983`** | trig-rational, exact on its leaf, `-(3/8)t^2 + (5/32)t^3 + (33/128)t^4` | **exact on leaf**, `= 3/8`; exhaustive over all 6,814 hole classes |
| 6 | `-0.396385 … -0.371908` @ `0.4…3.2°` | **`0.3999999`** | — | `= 2/5` to 7 digits; permutation family only |

**`T = 5` is not under-optimised.**  §3.2 of `notes/proof-architecture.md` flagged `c_5 <= 0.3745`
as suspicious (non-monotone).  Three independent searches now agree to every digit printed:

| `t` | `0.1°` | `0.2°` | `0.4°` | `0.8°` | `1.6°` | `3.2°` | `6.4°` |
|---|---|---|---|---|---|---|---|
| `runs/s6local_n20_T5_uniform.txt` (6,000 subsets, 18,004 starts) | | `-0.37445` | `-0.37390` | `-0.37277` | `-0.37044` | `-0.36551` | `-0.35460` |
| `runs/s6local_n20_T5_uniform_seed7.txt` (20,000 subsets, 60,000 starts) | | | `-0.37390` | | `-0.37044` | | `-0.35460` |
| `runs/arch_tledger_T5_holes_all.txt` (**all** 6,814 hole classes x 10 jitters) | `-0.374727` | `-0.374451` | `-0.373897` | `-0.372769` | `-0.370440` | `-0.365505` | |

The non-monotonicity `3/4, 1/3, 3/8` is real, and §2 explains it: `T = 2, 3` sit on a different
branch of the same formula.

### 1.3 The optimal configuration at `T = 5`

Winning final hole set (`runs/arch_tledger_T5_holes_all.txt`, identical at every `t` tried),
holes `(0,1), (1,3), (2,4), (3,0), (4,2)` up to the dihedral group — a **permutation of cycle type
`3 + 2`**, not a 5-cycle:

    # . # # #        row r is empty at column sigma(r),  sigma = (1, 3, 4, 0, 2)
    # # # . #        cycles:  0 -> 1 -> 3 -> 0   and   2 -> 4 -> 2
    # # # # .
    . # # # #
    # # . # #

Two more permutation classes tie to `2e-4` relative (`sigma = (1,2,4,0,3)`, a 5-cycle, and
`sigma = (1,4,2,0,3)`, type `4+1`); everything else is on the `-5/8 t^2` branch.  For the record the
`T = 4` optimum is `sigma = (1,3,0,2)`, a **4-cycle** ("tiling minus a 4-cycle", `S6_LOCAL.md` §0),
and the `T = 3` optimum is `sigma = (1,0,2)`, a transposition.

**The whole spectrum has exactly two leading values.**  At `T = 4` the final hole sets split into
`-0.330227 t^2` (one class) and `-0.651…-0.661 t^2` (all the others); at `T = 5` into
`-0.3705 t^2` (three classes) and `-0.595…-0.609 t^2`.  Extrapolated these are

    (T-2)/(2(T-1))   and   T/(2(T-1)) ,     difference exactly 1/(T-1).

`1/(T-1)` is the dual weight of a chain link (`1/2` at `T = 3`, `1/3` at `T = 4`,
`S6_LOCAL.md` §3) — see §2, where the gap is derived as `(1 unit of rise) x (chain weight)`.

### 1.4 `T = 6`: `c_6 = 2/5`, confirmed to seven digits

`n = 30`, `C(36,6) = 1,947,792` hole sets — too many to enumerate, so only the `720 -> 115`
permutation classes were searched (`runs/arch_tledger_T6_perm*.txt`; note a feasible point is an
**upper** bound on `-delta*(t)/t^2`, so a *better* search reports a *smaller* number, and at `T = 4`
6 jitters per permutation class report `0.6605` where 400 report the true `0.3302`):

| `t` | 60 jitters/class | **200 jitters/class** | `2/5 - 0.52 t - …` |
|---|---|---|---|
| `0.4°` | — | **`0.396385`** | `0.396370` |
| `0.8°` | — | **`0.392800`** | `0.392739` |
| `1.6°` | `0.391029` | **`0.385719`** | `0.385479` |
| `3.2°` | `0.381996` | **`0.371908`** | `0.370958` |

Extrapolating the 200-jitter column to `t = 0` (`arch_chain.py fit --T 6`; cubic fit on all four
points, quadratic on the three smallest — they agree):

    c_6  ->  0.3999999     against   (T-2)/(2(T-1)) = 2/5 = 0.4000000     (difference 1e-07)

**so the conjecture of §3.1 is confirmed at `T = 6` to seven digits**, with `d_6 ~ 0.520` — *not*
the `3/25 = 0.12` of the two-point third-order fit of §3.1.

**The mechanism transfers cleanly.**  `arch_chain.py analyse` on a `T = 6` optimum finds the same
certificate: a wall-to-wall `x`-chain of **six** (and a `y`-chain of six by the pinwheel symmetry),
span `4.972468 = T - u` to 7 digits, **rise `R = 1.000865`** at `t = 1.6°`, bound
`B_6 = -3.04931e-04` against `delta* = -3.04933e-04`.  Link sequence `+t, +t, -0.94t, 0, (1-1.02t)`
— the `T = 4, 5` bookkeeping with one more spare link.  Optimal hole sets are permutations, as at
every smaller `T` (best at `t = 3.2°`: cycle type `4+1+1`, `sigma = (2,1,5,0,4,3)`).

**The spectrum is a ladder with spacing `t/(T-1)`, not two branches.**  The top four final hole
classes at `t = 3.2°` are

    0.371908,  0.381996,  0.392085,  0.392682          gaps  0.010088,  0.010089,  0.000597

and `t/(T-1) = 0.011170`.  A gap of `t/(T-1)` in `-delta*/t^2` is exactly **one unit of `r_2` in the
rise `R = 1 + r_2 t^2`** (since `dB/dR = sin t/(T-1)`).  So these classes all share the *same*
leading constant and differ only in the `t^3` term: reading them against `c_6 = 2/5` gives
`r_2 = 2.02, 1.11, 0.21, 0.16` — consecutive integers plus a small offset.  That is consistent with
`c_6 = (T-2)/(2(T-1)) = 2/5`, with the caveat that the `t^3` coefficient is `d_6 ~ 0.50` for the best
class, not the `3/25 = 0.12` of the two-point fit in §3.1.

**The one caveat at `T = 6`:** only the permutation family was tried (`C(36,6)` is too large to
enumerate).  At `T <= 5` the exhaustive scans justify that restriction; at `T = 6` the `r_2` ladder
means the optimum is a class the `T <= 5` optima do not predict (`r_2 ~ 2.1` against `1/6` and
`1/8`), so a non-permutation optimum cannot be excluded there the way it can at `T <= 5`.
Reproduce with

    python3 search/arch_ct.py holes --T 6 --family perm --tdeg 0.4,0.8,1.6,3.2 -j 200 --nproc 8 --show

**The sign is untouched either way**: `c_6 ~ 0.4` is positive and larger than `c_4 = 1/3`, as the
formula requires.

---

## 2. Why: one exact inequality

### 2.1 The chain bound

> **Chain lemma.**  Let `Q_1, ..., Q_k` be unit squares at a **common** tilt `t`, centres `c_i`, such
> that consecutive pairs are separated along the common edge normal `e = (cos t, sin t)` with gap
> `>= 1 + delta` (this is the pair row of `S6_SKELETON.md` §3.1: for equal tilts the width
> `m = 1/2 + (|cos 0| + |sin 0|)/2 = 1`).  Summing the `k-1` rows telescopes:
>
>     (c_k - c_1) . e  >=  (k - 1)(1 + delta) .
>
> If in addition `Q_1` touches the wall `x = 0` and `Q_k` the wall `x = T` (so `(c_1)_x = P`,
> `(c_k)_x = T - P` with `P = u/2`, `u = cos t + sin t`), then with `R := (c_k - c_1)_y` and `k = T`
>
>     delta  <=  B_T(R, t)  :=  [ (T - u) cos t + R sin t ] / (T - 1)  -  1 .

Nothing is approximated; this is one line of arithmetic from the LP rows.  `R` — the chain's
**rise** — is the only free quantity, and

    B_T(R, t)  =  [ (R - 1) t + (1 - T/2) t^2 + O(t^3) ] / (T - 1) .

So with `R = 1 + r1 t + O(t^2)`,

    B_T  =  - c(r1, T) t^2 + O(t^3) ,      **c(r1, T) = ( T/2 - 1 - r1 ) / (T - 1)** .

`dB/dR = sin t / (T-1)`: the sensitivity of the bound to the rise is exactly the chain-link dual
weight `1/(T-1)` of `S6_LOCAL.md` §3, which is why the two branches in §1.3 differ by `1/(T-1)`.

### 2.2 The two branches

`arch_chain.py analyse` reads an optimal configuration, finds the tight links, splits them by which
normal they use, walks every wall-to-wall chain and reports the smallest rise.  Results:

| `T` | `t` | binding chain | span (`= T-u`) | rise `R` | `r1 = (R-1)/t` | `B_T` | `delta*` |
|---|---|---|---|---|---|---|---|
| 3 | `0.8°` | `y`-chain of 3 | `1.986136` | `+0.986376` | `-0.9757` | `-1.431726e-04` | `-1.431731e-04` |
| 4 | `0.8°` | `x`- and `y`-chain of 4 | `2.986136` | `+1.000033` | `+0.0024` | `-6.437992e-05` | `-6.437945e-05` |
| 5 | `3.2°` | `x`- and `y`-chain of 5 | `3.945738` | `+1.000585` | `+0.0105` | `-1.140111e-03` | `-1.140113e-03` |

The chain bound **is** `delta*` to 6 digits: the optimum's certificate is this one chain.  And the
rise is quantised:

* `r1 = 0`  (`R = 1 + O(t^2)`): the chain climbs **exactly one cell**, `c_T = (T-2)/(2(T-1))`.
  Observed at `T = 4, 5` (and conjectured for all `T >= 4`).
* `r1 = -1` (`R = 1 - t + O(t^2)`): the chain climbs one cell and loses a `t` zigzag,
  `c_T = T/(2(T-1))`.  Forced at `T = 2` (there `R <= T - u = 1 - t + t^2/2`, the container itself)
  and at `T = 3`; at `T >= 4` it is the *second*-best family.

Anatomy of the rise (read off the `T = 4` optimum, `runs/arch_tledger_T4_show.txt`).  In the
relaxed tiling a chain link is one of

| link | `(Dx, Dy)` | contribution to `R` |
|---|---|---|
| level | `(~1, 0)` | `0` |
| level, zigzag up | `(1 - t^2/2 + delta, +t)` | `+t` |
| level, zigzag down | `(1 + t^2/2 + delta, -t)` | `-t` |
| diagonal up-right | `(1 - t, 1 - t)` | `1 - t` |

**An `x`-chain can never descend by a cell.**  A *down*-diagonal step `(1,-1)` has
`e . (1,-1) = cos t - sin t < 1`, so it is not a link of an `x`-chain at all: it is separated along
`e'`, and is a link of a `y`-chain — which for the same reason runs up-and-**left**.  That
handedness is the pinwheel swirl, and it is why `R >= 1 - O(t)` whenever no full row of `T` cells is
occupied (which is forced: `n = T(T-1)` with a permutation hole set leaves every row one short).

Observed link sequences (measured, `arch_chain.py analyse`):

| `T` | links of the binding chain | `R` |
|---|---|---|
| 3 | diagonal `(+1)`, level zigzag `(-t)` | `1 - t` |
| 4 | level zigzag `(+t)`, diagonal `(1-t)`, level `(0)` | `1` |
| 5 | zigzag `(+t)`, zigzag `(+t)`, diagonal `(1-t)`, zigzag `(-t)` | `1` |

**At `T >= 4` the chain has `>= 3` links, and the configuration can spend one of the spare level
links zigzagging upward; the `T = 3` chain has only two links and its single level link zigzags the
wrong way.**  That is the whole difference between the two branches, and why `c_2, c_3` are off the
line.  *(Observation, not proof: whether the spare-link argument is forced at every `T >= 4` is the
combinatorial gap in the conjecture of §3.1.)*

### 2.3 The global statement, and the sign

Solving `B_T(R,t) <= 0` for `R` (`arch_chain.py threshold`, verified symbolically and on a
`6 x 8999` grid with 0 mismatches):

>     B_T(R, t) <= 0    <=>    R  <=  R*(T, t)  =  T tan(t/2) + cos t - sin t
>
>     R*(T, t) = 1 + (T/2 - 1) t - t^2/2 + O(t^3) .

Three specialisations, all exact and all elementary:

* **`R = 1`** (the tiling).  `R*(T,t) >= 1  <=>  T >= 1 + cos t + sin t = 1 + u`.  Since
  `max u = sqrt2`, this holds **for every `t` as soon as `T >= 1 + sqrt2 = 2.414…`, i.e. `T >= 3`**,
  with equality only at `t = 0 mod 90°`.  At `T = 2` it fails for every `t` in `(0°, 90°)` — which
  is why `s(2) = 2` needs the container argument and not a chain.
* **`R = (T-1) sin t`** (a straight row of `T`: consecutive centres differ by the edge vector).
  `B_T > 0  <=>  T cos t + sin t < T  <=>  tan(t/2) > 1/T` — the straight-chain window of
  `notes/proof-architecture.md` §3.2, i.e. **Part 3**; see §4 below.
* **`R = 0`**: `B_T = (T-u) cos t/(T-1) - 1 < 0` for every `t` in `(0°, 90°)`.

**Consequence for §3.3.**  `c_T = (T-2)/(2(T-1))` is positive and increasing for all `T >= 3`, so
the second-order constant at the tiling **never changes sign**, at `T = 17` no less than at `T = 4`.
Scenario (α) — "`P_T` is born at the tiling" — is refuted for the uniform-tilt direction, and, by
the global form above, refuted for *every* `t` on any configuration that keeps one chain of rise
`<= 1`.  What a counterexample at `T = 17` must do instead is destroy the chain: either there is no
wall-to-wall chain of `T` squares at all (the tilted band of scenario (β) cuts every row and
column), or every such chain rises by `2` or more.  Both are statements about **A4**, the lemma with
no owner.

*Caveat, stated plainly.*  `arch_ct` produces feasible points, so §1 gives `c_T <= (T-2)/(2(T-1))`;
the exhaustive hole-set scan at `T <= 5` says nothing better exists **on the tiling**, and
`s6exact` (§3) makes it exact on one leaf, but `delta*` is a max over all configurations and the
search is a search.  The sign conclusion is robust for a different reason: `B_T(1,t) <= 0` for all
`t` and all `T >= 3` is a *theorem*, not a measurement, so the tiling stratum cannot be where a
counterexample is born as long as a rise-`1` chain survives.

---

## 3. Exact closed form on the `T = 5` leaf

`python3 search/s6exact.py --n 20 --T 5 --pattern uniform --tdeg 1.6`
(`runs/arch_tledger_s6exact_T5.txt`): the optimal leaf has **27 dual rows of 270**, and solving the
dual symbolically on that support gives a trig-rational `delta(t)` (the full expression is in the
run file; it agrees with the LP at the sample to `3e-13`) with series

    delta(t)  <=  - (3/8) t^2  +  (5/32) t^3  +  (33/128) t^4  +  O(t^5)        [n = 20, T = 5]

**`c_5 = 3/8` exactly on the leaf.**  The dual is the same shape as at `T = 3, 4`: a wall-to-wall
staircase chain of **five** squares (`lo-x 1, pair(1,6), pair(6,10), pair(10,14), pair(14,18),
hi-x 18`) at weight `1/4 - (3/16) t + …` each, closed by a transverse chain at `t/4` and links at
`t^2/4`.  So the brief's guess is right: **the chain weight is `1/(T-1)`**, and in fact

| `T` | chain weight | transverse | links | `c_T` | `d_T` (`t^3`) |
|---|---|---|---|---|---|
| 3 | `1/2 - t/4` | `t/2` | `t^2/2` | `3/4` | `9/8` |
| 4 | `1/3 - 2t/9` | `t/3` | `t^2/3` | `1/3` | `2/9` |
| 5 | `1/4 - 3t/16` | `t/4` | `t^2/4` | `3/8` | `5/32` |

i.e. chain weight `= 1/(T-1) - (T-2) t/(T-1)^2 + … = [1 - 2 c_T t]/(T-1)`, transverse `t/(T-1)`,
links `t^2/(T-1)`.

### 3.1 The conjectured formula

> **Conjecture (leading order — this is the answer to §3.2's question).**  At `n = T^2 - T` in
> `[0,T]^2`, for `T >= 4`,
>
>     delta*(t,…,t)  =  - c_T t^2 + O(t^3) ,       **c_T = (T - 2) / (2 (T - 1))** ,
>
> attained by the tiling with a permutation hole set; equivalently, the rise of the binding
> wall-to-wall chain of `T` is `R = 1 + O(t^2)`.
> For `T = 2, 3` the other branch applies: `c_T = T/(2(T-1))` (`c_2 = 1`, `c_3 = 3/4`).
>
> `c_T` is `> 0` for every `T >= 3` and increases strictly to `1/2`.  **It never changes sign.**

*Evidence:* exact on the leaf at `T = 2, 3, 4, 5`; at `T = 5` confirmed by an exhaustive scan of all
6,814 dihedral classes of hole set; at `T = 6` the permutation-family search extrapolates to
`0.3999999` against `2/5` (§1.4).

**A weaker, third-order guess, which `T = 6` does not support.**  Writing `R = 1 + r_2 t^2` and
feeding it into §2.1 gives `d_T := [t^3] delta* = (1/2 + r_2)/(T-1)`; `d_4 = 2/9` and `d_5 = 5/32`
both give `r_2 = 1/(2(T-1))`, i.e. `d_T = T/(2(T-1)^2)`.  That is a one-parameter fit through two
points, and `T = 6` breaks it: the best class there needs `d_6 ~ 0.52`, not `3/25 = 0.12`, and the
hole classes form a ladder of `r_2` values spaced by `1` (§1.4).  **Treat the `t^3` coefficient as
unknown for `T >= 6`;** only the `t^2` coefficient is claimed.

Checks (`3/8 - (5/32) t - (33/128) t^2` against the exhaustive search at `T = 5`):

| `t` | `0.1°` | `0.2°` | `0.4°` | `0.8°` | `1.6°` | `3.2°` |
|---|---|---|---|---|---|---|
| leaf series | `0.3747265` | `0.3744514` | `0.3738966` | `0.3727681` | `0.3704356` | `0.3654692` |
| exhaustive search `-delta*/t^2` | `0.374727` | `0.374451` | `0.373897` | `0.372769` | `0.370440` | `0.365505` |

(the search row is `runs/arch_tledger_T5_holes_all.txt`, all 6,814 hole classes; it agrees with the
exact leaf series to **6-7 digits at every `t` from `0.1°` to `3.2°`**, which is as strong a confirmation of `c_5 = 3/8` as
a feasible-point method can give).

**Sign.**  `c_T = (T-2)/(2(T-1))` is `0` at `T = 2`, `1/3` at `T = 4`, and increases strictly to
`1/2`.  **`c_T` does not change sign at any `T`**, and neither branch does (`T/(2(T-1))` also
increases to `1/2`).  The sequence
`1, 3/4, 1/3, 3/8, 2/5, 5/12, …` is non-monotone only because of the branch change at `T = 4`;
§3.2 of `notes/proof-architecture.md` read `1, 3/4, 1/3` as three points of one curve and
extrapolated through zero.  They are two points of one branch and one of the other.

### 3.1a The leaf bounds need no radius on `(0°, 45°]`

`A2c` in `notes/proof-architecture.md` asks for an *explicit radius*.  On these leaves there isn't
one to give: `delta*(t,…,t) = delta*(90° - t)` (reflect the whole configuration in `y = x`), so
`(0°, 45°]` is the whole story, and `arch_chain.py leaf` (`runs/arch_tledger_leafcheck.txt`)
evaluates the four exact forms on a `0.01°` grid:

| `T` | `5°` | `15°` | `28°` | `45°` | max on `(0°,45°]` | formula's first sign change |
|---|---|---|---|---|---|---|
| 2 | `-0.006947` | `-0.050510` | `-0.124199` | `-0.171573` | `-0.000000` (at `t -> 0`) | none below `90°` |
| 3 | `-0.004994` | `-0.033674` | `-0.077158` | `-0.102944` | `-0.000000` | none below `90°` |
| 4 | `-0.002390` | `-0.018799` | `-0.052917` | `-0.085786` | `-0.000000` | `58.38°` |
| 5 | `-0.002738` | `-0.021898` | `-0.060388` | `-0.077988` | `-0.000000` | `53.62°` |

**Strictly negative on the whole of `(0°, 45°]`, at every `T` for which an exact leaf form exists.**
(The `T = 4, 5` formulas turn positive above `~54°`; that is the fixed assignment ceasing to be
optimal, not `delta*` turning positive, and it is outside the range anyway.)

### 3.2 The certificate in the `T = 5` configuration

The binding chain (`arch_chain.py analyse`, §2.2's table) is
`(0.527,3.417) (1.524,3.472) (2.522,3.528) (3.469,4.473) (4.473,4.417)` at `t = 3.2°`: two level
links zigzagging `+t` each, one diagonal up-right (`~1-t`), one level link zigzagging `-t`; net
rise `1.000585 = 1 + O(t^2)`, against the predicted `1 + t^2/8 = 1.00039`.  Its span is
`3.945738 = T - u` to 7 digits — wall-to-wall — and its bound `B_5 = -1.140111e-03` equals
`delta* = -1.140113e-03`.

---

## 4. Part 3: the straight-chain window `tan(t/2) >= 1/T`

`arch_ct.py window`.  **Verified, with one correction to the statement in
`notes/proof-architecture.md` §3.2.**

*Proof.*  A straight row of `T` unit squares at a common tilt `t in (0°, 90°)` has bounding box
`(T cos t + sin t) x (cos t + T sin t)`.  For the first factor,

    T cos t + sin t <= T
      <=>  sin t <= T (1 - cos t)
      <=>  2 sin(t/2) cos(t/2) <= 2 T sin^2(t/2)
      <=>  cos(t/2) <= T sin(t/2)                       (sin(t/2) > 0 on (0°, 180°))
      <=>  tan(t/2) >= 1/T .                                                          QED

**Correction.**  The row must fit *in both directions*; the transverse extent
`cos t + T sin t <= T` is the same condition with `t -> 90° - t`.  So the window is
`t in [t_T, 90° - t_T]` with `t_T = 2 arctan(1/T)`, not a half-line, and it is **empty iff
`t_T > 45°`, i.e. iff `1/T > tan 22.5° = sqrt2 - 1`, i.e. iff `T < 1 + sqrt2 = 2.4142…` — only
`T = 2`.**  The thresholds quoted in `notes/proof-architecture.md` §3.2 are right:

| `T` | 2 | 3 | 4 | 5 | 6 | 10 | 17 | 100 |
|---|---|---|---|---|---|---|---|---|
| `t_T = 2 arctan(1/T)` | `53.1301°` | `36.8699°` | `28.0725°` | `22.6199°` | `18.9246°` | `11.4212°` | `6.7329°` | `1.1459°` |
| window | empty | `[36.87, 53.13]` | `[28.07, 61.93]` | `[22.62, 67.38]` | `[18.92, 71.08]` | `[11.42, 78.58]` | `[6.73, 83.27]` | `[1.15, 88.85]` |
| `2/T` in degrees | `57.30°` | `38.20°` | `28.65°` | `22.92°` | `19.10°` | `11.46°` | `6.74°` | `1.15°` |

`T cos t_T + sin t_T = T` to machine precision at every `T` in the table (0 mismatches of the
equivalence over `6 x 8999` grid samples, `arch_ct.py window`).

**What it costs the architecture.**  `t_T = 2 arctan(1/T) = 2/T - O(1/T^3)`, so the tube radius in
`A2` cannot exceed `~2/T`, and at the edge of the tube the deficit the local theorem can deliver is

| `T` | 3 | 4 | 5 | 6 | 10 | 17 | 100 |
|---|---|---|---|---|---|---|---|
| `B_T(1, t_T)` | `-6.00e-2` | `-6.46e-2` | `-5.18e-2` | `-4.03e-2` | `-1.72e-2` | `-6.42e-3` | `-1.98e-4` |
| `-c_T t_T^2` | `-1.04e-1` | `-8.00e-2` | `-5.84e-2` | `-4.36e-2` | `-1.77e-2` | `-6.47e-3` | `-1.98e-4` |

`c_T -> 1/2` but `t_T^2 -> 4/T^2`, so the handover to `A3` is worth `~2/T^2` and shrinks
quadratically.  **This is the `T`-dependence that a "uniform in `T`" reading of `A2c` hides**: the
lemma stays true with a *better* constant and a *worse* domain.

---

## 5. Part 2: Cleemann's packing and the smallest known counterexample

*Literature sweep (web).  Verbatim quotes marked **[pub]**; measurements taken off the DS7 figure
by a subagent running connected-component analysis on the `343 x 343` GIF are marked **[fig]** and
are **not published data**; inferences are **[inf]**.*

### 5.1 `T = 17` is no longer the frontier — the ledger in §3 is out of date by two steps

| `T` | `n = T^2-T` | best known `s(n)` | who, when |
|---|---|---|---|
| 4 … 10 | 12 … 90 | `= T` (trivial packing best known) | **open** |
| **11** | 110 | **`10.99679327401957…`** | **David W. Cantrell, February 2025** (exact analytic solution; not refereed) |
| 12 | 132 | `< 11.99790201731` | Arslanov–Mustafin–Shangitbayev, March 2019 (**refereed**, EJC 2021) |
| 13, 14, 15 | 156, 182, 210 | `< T` | same |
| 16 | 240 | `15.9755…` | Károly Hajba, September 2015 |
| 17 | 272 | `< 17` | **Lars Cleemann, between 1991 and 1998** |

**The smallest `T` for which `s(T^2-T) < T` is known is `T = 11`** (Cantrell, Feb 2025,
`s(110) <= 10.99679327401957223263170298166840268704638386123877`, full coordinates and angles to
~50 digits in <https://kingbird.myphotos.cc/packing/square-110.svg>).  The smallest *refereed* one
is `T = 12`:

> **[pub]** M. Z. Arslanov, S. A. Mustafin, Z. K. Shangitbayev, *Improved packings of `n(n-1)` unit
> squares in a square*, Electron. J. Combin. **28**(4) (2021) #P4.22, DOI 10.37236/8586, abstract:
> "Let `s(n)` be the side of the smallest square into which we can pack `n` unit squares.  The
> purpose of this paper is to prove that `s(n²−n) < n` for all `n ⩾ 12`.  Besides, we show that
> `s(18²−17) < 18`, `s(17²−16) < 17`, and `s(16²−15) < 16`."

> **[pub]** same, p. 1: "An important question is to find the minimum `n` for which `s(n²−n) < n`.
> For small `n`, only `s(2) = 2` and `s(6) = 3` have been proved, but we dont even know the proof of
> `s(12) = 4`."

So **`notes/proof-architecture.md` §3's "the first failing `T*` is somewhere in `5..17`" should read
`5..11`** (or `5..12` if only refereed results count).  `T = 4` (our target) and `T = 5..10` are
open; `T = 5` is *nine* steps from a known counterexample, not twelve.

Current-records page (maintained by David Ellsworth, successor to Friedman's Packing Center):
<https://kingbird.myphotos.cc/packing/squares_in_squares.html>; its `n = 110` entry reads **[pub,
non-refereed]** "Found by David W. Cantrell in February 2025.  Bounds the s(n^2-n)=n conjecture to
n < 11."  For `n = 12` and `n = 20` the page lists nothing, which by its own convention means the
trivial packing is still the best known.  Best known bounds: `2 + 4/sqrt5 = 3.7888 <= s(12) <= 4`
(Stromquist, DS7 Table 2) and `6 sqrt2 - 4 = 4.4852 <= s(20) <= 5` (Friedman, DS7 Table 2).

### 5.2 Cleemann's 272-square packing: what is published

> **[pub]** E. Friedman, *Packing Unit Squares in Squares: A Survey and New Results*, Electron. J.
> Combin. Dynamic Survey **DS7**, DOI 10.37236/28 (v5, 14 Aug 2009;
> <https://www.combinatorics.org/files/Surveys/ds7/ds7v5-2009/ds7-2009.html>):
>
> "It was conjectured that s(n2-n)=n whenever n is small.  The smallest known counterexample of this
> conjecture, due to Lars Cleemann, is s(172-17)<17.  That is, 272 squares can be packed into a
> square of side 17 in such a way that the the square can be squeezed together slightly (see
> Figure 8).  **Three squares are tilted by an angle of 45o, and the other tilted squares are tilted
> by an angle of arctan(8/15).**"
>
> Caption: "s(272)<17", **Figure 8** (v5 2009); the same passage is **Figure 10** in v1 (1998) and
> **Figure 11** in v2 (2000).  Image: `pic/17x17.gif`.

**Published quantitative content: the strict inequality `< 17`, the two tilt angles, and "three
squares at 45°".  No coordinates, no exact side, no `epsilon`.**  Arslanov et al. confirm this
**[pub]**: "in [2] a sporadic squeezable packing of 272 unit squares in a square (17,17) is given,
proving that `s(17²−17) < 17`, but from this it does not follow that `s(18²−18) < 18` etc."  Date:
DS7 carries no citation; Ellsworth's page says **[pub, non-refereed]** "Originally found by Lars
Cleemann between 1991 and 1998."  (The `s(272) = 16.96971…` now listed there is *not* Cleemann's
packing — it is converted from Arslanov et al.'s `s(210)`.)

### 5.3 The structure, read off the figure

**[fig]** Connected-component analysis of the `343 x 343` GIF finds **exactly 272 components** and

| class | count |
|---|---|
| axis-parallel (within 3°) | **199** (73 %) |
| `+arctan(8/15) = +28.0725°` | 35 |
| `-arctan(8/15)` (mirror, `61.93°`) | 35 |
| exactly `45°` | **3** |

— the "three at 45°" matches DS7's text exactly, which validates the extraction.  Occupancy map
(`#` = axis-parallel square on the integer grid, `/` = cell cut by the tilted band, row 0 on top):

       01234567890123456
     0 ##//#############
     1 ////#############
     2 ////#############
     3 #////############
     4 #////############
     5 #////############
     6 #////############
     7 ##////###########
     8 ##////###########
     9 ##////###########
    10 ##////###########
    11 ###///////#######
    12 ####//////////###
    13 #####////////////
    14 ######///////////
    15 ##########//////#
    16 ##############//#

**The band is bent, not straight**: a steep `+28°` branch down the upper left, a shallow `-28°`
branch across the lower right, with the three `45°` squares at the elbow as the hinge.  The 199
axis-parallel squares sit on the exact integer grid, in four blocks.

**[fig] The decisive observation for §3.3: no row and no column of the container carries a
wall-to-wall run of 17 axis-parallel squares — the band cuts every one of the 17 rows and every one
of the 17 columns.**  (Read the map: every row has at least one `/`, and so does every column.)

### 5.4 Classification: this is **(β)**

§3.3's three scenarios: (α) `c_T` changes sign at the tiling; (β) elsewhere on `Z` — the chain
lemma leaks; (γ) away from `Z` — the far field.

* (α) is refuted outright by §2: `c_T = (T-2)/(2(T-1)) > 0` and increasing, and `B_T(1,t) <= 0` for
  every `t` once `T >= 3`.
* Cleemann's packing has **199 of 272 squares axis-parallel** — so `k = 199 >> T = 17`, i.e. it is
  deep inside the region `A3` hands to `A4`, not in the far field.  And it has **no wall-to-wall
  chain of 17 axis-parallel squares**.  That is exactly the hypothesis of `A2a`/`A4` failing.
* Therefore: **(β)**, as §3.3 guessed.  The statement that dies is `A4`:
  *"`k >= T` near-axis squares contain a wall-to-wall chain of `T`"*.  Cleemann's band is a
  codimension-one obstruction threaded through the whole container: it cuts every row and column, so
  Mirsky/Dilworth on the near-axis squares alone cannot produce a chain — the tilted squares are not
  near-axis and break the acyclicity argument that `A4` relies on.

`A4`'s counting also fails for a more elementary reason here: `199 <= (17-1)^2 = 256`, so the
axis-parallel corollary of K2 (`S6_SKELETON.md` §2.2) is satisfied with room to spare.  The
`(T-1)^2` bound and the total-band inequality are *silent*, exactly as §2.3 of `S6_SKELETON.md`
predicted.

### 5.5 The tilt angle is our Part-3 threshold, exactly

    arctan(8/15) = 2 arctan(1/4) = 28.0724869358529…°,     tan(t/2) = 1/4 exactly,
    4 cos t + sin t = 4·(15/17) + 8/17 = 68/17 = 4          EXACTLY.

So **Cleemann's tilt is the `k = 4` straight-chain window at equality** (§4): a straight stack of
four unit squares at that angle spans a width-4 strip with *no slack at all* while occupying only
`cos t + 4 sin t = 47/17 = 2.7647` in the other direction.  Arslanov et al. use the same device
explicitly **[pub]**: "In Figure 2 we see one of the main ideas for packing unit squares: using of
stacks (4,1) tilted by an angle `α = arcsin(8/17)`."  (`arcsin(8/17) = arctan(8/15)`, same angle.)

**[inf] What that means for the ledger.**  The `(4,1)` stack exists in *any* container of side
`>= 4`; there is nothing 17-ish about it.  The window inequality of §4 is `T`-free once you read the
`4` as the *stack length*, not the container.  What is 17-ish is having enough room around the band
to pay for the elbow: the band wastes area locally (four squares of area 4 occupying a
`4 x 2.7647` box) and has to be repaid by the `199` axis-parallel squares elsewhere.  That is an
**area/counting** statement — `A3`'s cover LP — and *not* an angle-geometry statement.  So the
`T`-dependence of the whole problem is split between:

* `A4`: the band cuts every row and column, so no chain (this is the qualitative break, **(β)**);
* `A3`/`A5`: whether the container is big enough for the band to pay for itself (this is the
  quantitative threshold, and it is what moved `T* <= 17` down to `T* <= 11` between 1998 and 2025).

Neither is `A2c`.

### 5.6 The one artefact worth fetching next

Cleemann's packing has **no published coordinates** — everything in §5.3 is measured off a
`343 x 343` raster.  The `T = 11` counterexample does have them, to ~50 significant digits, with an
exact analytic solution:

    https://kingbird.myphotos.cc/packing/square-110.svg        s(110) = 10.99679327401957223263…

Per-square centres and rotation angles are in the SVG.  That is the right input for the `A4`
falsification search of `notes/proof-architecture.md` §4: run `search/chains.py` / `arch_chain.py`
on a *real* `s(T^2-T) < T` packing and see exactly which chain statement it violates, at the
smallest `T` where one exists.  Caveat: unlike Cleemann's, it is a numerically optimised packing
with many distinct tilt angles (`23.92°, 26.41°, 27.33°, 28.12°, 28.28°, 60.52°, …`), so it is less
diagnostic of a *mechanism* and more of an existence proof.  Cleemann's two-angle packing is the
better model; the `s(110)` file is the better data.

---

## 6. What this does to §3 of `notes/proof-architecture.md`

| §3 claim | verdict |
|---|---|
| "`c_T = 1, 3/4, 1/3` at `T = 2,3,4` and `<= 0.375` at `T = 5` … extrapolates through zero before `T = 5`" | **wrong extrapolation.**  `c_5 = 3/8` exactly; the sequence is `1, 3/4, 1/3, 3/8, 2/5, …` with a *minimum* at `T = 4` and limit `1/2`.  `T = 2, 3` are a different branch of the same formula (§2.2). |
| "the measured `T = 5` value is … either a parity effect or the same under-optimisation" | **neither.**  Exhaustive over all 6,814 dihedral classes of hole set; three searches agree. |
| "if `c_T` changes sign at some `T_loc`, … `A2c` is the lemma that carries `T < T_loc`" | **the antecedent is false.**  `c_T` does not change sign at any `T`. |
| "**My guess is (β)**, which would mean the `T`-dependence sits in `A4`" | **confirmed** by the Cleemann figure (§5.4). |
| "the first failing `T*` is somewhere in `5..17`" | **`5..11`** (Cantrell 2025), or `5..12` refereed (Arslanov et al. 2021). |
| "`tan(t/2) >= 1/T` … `t >= 53.1°, 36.9°, 28.1°, 22.6°, 6.7°`" | **verified**; add that the window is two-sided, `[t_T, 90° - t_T]`, and empty only at `T = 2` (§4). |
| "The tube radius in `A2` cannot exceed `~2/T`" | **verified and sharpened**: the deficit available at the edge of the tube is `~2/T^2` (§4), which is what `A3` must beat. |

**Revised §3.3 answer.**  `P_T` is born at (β), on `Z` but away from the tiling, and the lemma that
carries `T < T*` is **A4** — "`k >= T` near-axis squares contain a wall-to-wall chain of `T`" —
with `A3`/`A5` carrying the quantitative threshold.  `A2c` is uniform in `T` and should be stated
that way; the honest version is not "`delta* <= -c_T |tilt|^2`" but

> **A2c′ (chain form, no `eps`, no constant to tune).**  If a configuration contains a wall-to-wall
> chain of `T` unit squares at a common tilt `t` whose rise is at most `1`, then
> `delta <= B_T(1,t) <= 0` for every `t`, with equality only at `t = 0 mod 90°`, for every
> `T >= 1 + sqrt2`.

which is §2.3, is exact, is two lines, and is *true at `T = 17`*.  The falsifier that `A4` now needs
is concrete: **Cleemann's band, which cuts every row and column of a `17 x 17` container while
leaving 199 axis-parallel squares.**  The `T = 4` version of the falsification search
(`notes/proof-architecture.md` §4, "chain lemma with obstacles") should be run with that shape in
mind: near-axis squares in blocks, a bent band of tilted squares at `2 arctan(1/4) = 28.07°` with a
`45°` hinge, and the question is how few squares such a band needs.  At `T = 4` a `(4,1)` stack
already spans the container exactly, so `n = 12` cannot afford the elbow — that, and not the angle
geometry, is what has to be proved.

## 7. Reproduce

    python3 search/arch_chain.py table        # c(r1,T), the two branches, exact B_T on each
    python3 search/arch_chain.py threshold    # R*(T,t) = T tan(t/2) + cos t - sin t, and (a)(b)(c)
    python3 search/arch_chain.py leaf         # the exact leaf bounds, sign over (0, 45] deg
    python3 search/arch_ct.py window          # Part 3
    python3 search/arch_ct.py holes --T 3 --family all  --tdeg 0.8 -j 40 --nproc 2 --show      #  2 s
    python3 search/arch_ct.py holes --T 4 --family all  --tdeg 0.8 -j  6 --nproc 8 --show      #  2 s
    python3 search/arch_ct.py holes --T 4 --family perm --tdeg 0.8 -j 400 --nproc 2            # 12 s
    setsid nohup python3 -u search/arch_ct.py holes --T 5 --family all \
        --tdeg 0.1,0.2,0.4,0.8,1.6,3.2 -j 10 --nproc 8 --topk 15 --seed 3 --show \
        > runs/arch_tledger_T5_holes_all.txt 2>&1 &                                            # 25 min
    setsid nohup python3 -u search/arch_ct.py holes --T 6 --family perm --tdeg 1.6,3.2 \
        -j 60 --nproc 2 --show > runs/arch_tledger_T6_perm.txt 2>&1 &                          # 10 min
    setsid nohup python3 -u search/arch_ct.py holes --T 6 --family perm \
        --tdeg 0.4,0.8,1.6,3.2 -j 200 --nproc 8 --seed 23 --show \
        > runs/arch_tledger_T6_perm_big.txt 2>&1 &                                             # 20 min
    python3 search/arch_chain.py analyse --T 4 --file runs/arch_tledger_T4_show.txt --tol 3e-6
    # (for T = 5, 6: cut the `winning configuration at t=...` block of the run into a file first)
    python3 search/arch_chain.py fit --T 5 --tdeg 0.2,0.4,0.8,1.6 --vals=-0.37445,-0.37390,-0.37277,-0.370440
    setsid nohup python3 -u search/s6exact.py --n 20 --T 5 --pattern uniform --tdeg 1.6 \
        --limit 3000 --jitters 4 --nproc 2 > runs/arch_tledger_s6exact_T5.txt 2>&1 &
