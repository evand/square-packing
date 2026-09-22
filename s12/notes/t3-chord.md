# t3-chord: `s(6) = 3` in the row language with line-transversal (chord) rows — the salvage attempt, and its result

*2026-09-22.  Task `tasks/t3-chord/README.md`.  Code (new; nothing in `search/` modified —
`s6skel`, `t3_chain`, `rank8` are imported):
`search/t3_chord_geom.py` (exact chord geometry, the leaf-linear endpoint forms, the two
premise checks), `search/t3_chord_ks.py` (Kearney–Shiu §3 in this language: the 13 points,
unavoidability, Lemma 1, Lemma 2 + (7), Lemma 3, and the configuration-(b) Farkas
certificate), `search/t3_chord_scan.py` (the chord functional on the recorded packings, the
coherent cone, the band-vs-pigeonhole table), `search/t3_chord_oblique.py` (chord rows on
lines of arbitrary direction).
Runs: `runs/t3_chord_geom.txt`, `runs/t3_chord_ks.txt`, `runs/t3_chord_scan.txt`,
`runs/t3_chord_oblique.txt`.
Source read for §2: Kearney–Shiu, *Efficient packing of unit squares in a square*, Electron. J.
Combin. 9 (2002) #R14, doi 10.37236/1631, §§2–3 in full (not the repo's summary).
Labels: **[proved]** = proved here by hand and re-checked numerically by a script here;
**[measured]** = produced by an exact scan or a multistart, hence one-sided; **[heuristic]**;
**[guess]**.*

---

## 0. Verdict, up front

> **K–S reproduced — but only its last three lines are a Farkas certificate; the pigeonhole
> substitute fails at step S = the existence step, and it fails at *first* order in the tilt
> while the truth is second order.  The angle-space route is closed.**
>
> In the brief's three forms the answer is the **second and the third at once, and they name
> the same step**: *K–S does not fit the row system at step `S` = the unavoidable-set
> counting*, which is its existence step and has no centre-variable expression; and *the
> pigeonhole substitute for that same step fails*, by `(T-2)t/(T-1)`, first order in the tilt.
> `s(6) = 3` is **not** proved in the row language.
>
> Concretely, three things.
>
> 1. **Chord rows are real, and on the certificate side they are exactly what obstruction X
>    asks for.**  A *banded column of three* — three squares whose centres lie within
>    `(C - S)/2` of one vertical (or horizontal) line — certifies
>    `delta <= -3(sec t - 1)/4 = -0.375 t^2 + O(t^4)` at a common tilt `t`, with equality `0`
>    at `t = 0`.  That is exactly tight in the axis-parallel stratum **and** correct to second
>    order off it, which no row in the existing system is (`t3-existence.md` §4.2: the bare
>    chain of three gives `+t/4`, first order and the wrong sign).  **[proved]**
> 2. **The existence side is strictly worse than before, not better.**  The clause the
>    certificate needs — "some axis-parallel line is met by three squares in band" — is
>    **false** at an explicit `delta = 0` packing (the `3 x 3` tiling minus a permutation: the
>    pinwheel, `max_c F = 2` on every axis-parallel line, and no line of any direction carrying
>    more chord than its container chord allows), whereas
>    `t3-existence.md`'s (E3-row) holds at all `749` recorded packings.  And the pigeonhole
>    that would produce it is **exactly critical at `t = 0` and fails at first order for every
>    `t > 0`**: six centres in an interval of length `T - u` give three inside a window of
>    length `(T-u)/2`, the certificate needs a window of length `C - S`, and
>    `(T-u)/2 <= C - S` iff `3 cos t - sin t >= 3`, i.e. iff `t = 0`.  The deficit is
>    `t/2 + O(t^2)` at `T = 3` and `(T-2)t/(T-1)` at general `T`.  **[proved]**
> 3. **Obstruction X is restated, not reduced.**  At the uniform-tilt optimum the leaf
>    hypothesis fails by `1.9e-04 ~ 0.63 t^2` at `t = 1 deg`, against `delta* = -1.1e-04
>    = -0.3637 t^2` — the *same* order.  (It must fail: the certificate would otherwise give
>    `-0.375 t^2`, which is **below** `delta*`; the two second-order constants agree to `3 %`.)
>    So the chord existence clause is a second-order statement about the centres, exactly like
>    X, and the only tool for it — pigeonhole — is accurate to `O(t)`.  **[measured]**

And on Kearney–Shiu itself: **their argument is not a Farkas certificate except in its final
three lines.**  Configuration (b) ends with a genuine certificate — five rows, weights
`(1,1,1,1,1)`, exact cancellation of every centre variable, value
`delta <= (13/6 + 2 - 3 sqrt2)/2 = -0.037987010226`.  Configuration (a), and every
intermediate step of (b), end instead in a **counting** contradiction ("this square would cover
two points of one lattice"), which is a statement about the unavoidable-set bijection and has
no expression in the centre-variable row language at all.  That is the step that resists, and
it is the same step — existence — that kills the pigeonhole substitute.

The `T = 4` transfer is in §6.  Short form: (a) a **13**-point pure unavoidable set for
`[0,4]^2` is *not* excluded by anything in the repo, and cannot be excluded by the repo's
method, because the fractional bound that would have to exceed `13` is itself believed to be
`~12.4`; what the repo *does* prove is that **no pure set of `12` or fewer points exists**, so
the slack-`0` route (`11` points, `t3-chain.md` §4.4) is dead and `13` is the smallest
conceivable size.  (b) The chord pigeonhole gives **nothing** at `T = 4` that the centre
pigeonhole does not: it is exactly critical at `t = 0` for every `T` and fails at first order
immediately off it, and it needs a *straight* column of `T`, which the `T x T`-tiling-minus-a-
permutation denies at every `T`.

What in the brief's premises, and in the repo's notes, turned out to be wrong is collected in
§7 — including two premises of this brief that are false as stated.

---

## 1. The chord row, exactly

### 1.1 The chord, and where it sits  [proved]

Let `Q` be a closed unit square with centre `c = (x, y)` and angle `th`, `C = cos th`,
`S = sin th`, `p = (|C| + |S|)/2`.  For the vertical line `L : x = a` put `e = a - x`.  A point
`(a, y + eta)` lies in `Q` iff

    |e C + eta S| <= 1/2      and      |eta C - e S| <= 1/2 ,

so `eta` ranges over the intersection of

    branch 0 :  [ (-1/2 - eC)/S , (1/2 - eC)/S ]        (length 1/|S|)
    branch 1 :  [ ( eS - 1/2)/C , ( eS + 1/2)/C ]       (length 1/|C|)

and therefore

> **(C0)**  `alpha(L) = y + max(lo_0, lo_1)`, `beta(L) = y + min(hi_0, hi_1)`,
> `h = beta - alpha = min( 1/|C|, 1/|S|, (p - |e|)/(|C||S|) )`  (`= s6skel.chord_len(th, |e|)`),
> and **on a leaf that fixes which branch is active at each end, `alpha` and `beta` are exactly
> affine in `(x, y)`** — with coefficient `(C/S, 1)` on branch `0` and `(-S/C, 1)` on branch `1`.

The leaf that fixes the branches is a **single-variable box row on the centre**: branch `0` is
active at the top iff `hi_0 <= hi_1` iff `e >= (C - S)/2`, i.e. iff `x <= a - (C-S)/2`.  The
horizontal-line case is the same with `x, y` and `C, S` exchanged.  Verified against an exact
polygon–line clip (max error `9.8e-15` over `4000` random cases) and against `s6skel.chord_len`
(max error `3.3e-15`): `python3 search/t3_chord_geom.py selftest`.

Two constants used throughout:

    band_full(th) = (C - S)/2   :  |e| <= band_full  <=>  h = 1/max(C,S) = sec(tilt)   (its maximum)
    D(th) = p - CS              :  |e| <= D          <=>  h >= 1                       (chord-lemma.md Cor. 2)

`band_full` decreases from `1/2` at `t = 0` to `0` at `45 deg`; `D` decreases from `1/2` to
`(sqrt2 - 1)/2`.

### 1.2 Two premises of the brief are false  [proved]

The brief states the row as: *"the chords `Q_i ∩ {x = c}`, `Q_j ∩ {x = c}` are disjoint closed
intervals, so `|y_i - y_j| >= (h_i + h_j)/2`."*  **Both halves are wrong.**

* **The chord is not centred at the centre.**  The two branch intervals have midpoints
  `-eC/S` and `+eS/C`, which agree only at `e = 0`, so the chord's midpoint is
  `y + mu(th, e)` with `mu` a nonzero piecewise-linear function of `e`.  Measured maximum
  offset `0.4728` (at `th = 87.9 deg`, `e = -0.5175`), i.e. the midpoint can be two thirds of a
  half-diagonal away from `y`.
* **Hence `|y_i - y_j| >= (h_i + h_j)/2` is false for disjoint squares.**  Explicit
  counterexample found by search and re-checked with `s6skel.pair_gap >= 0`:

  | | `th` | centre | chord on `x = 0.911` |
  |---|---|---|---|
  | `Q_1` | `73.003 deg` | `(1.5, 1.5)` | length `0.1263` |
  | `Q_2` | `73.718 deg` | `(0.4, 1.5084)` | length `0.4057` |

  `|y_1 - y_2| = 0.0084` while `(h_1 + h_2)/2 = 0.2660`: the inequality is violated by
  `0.2576`.  The two squares are disjoint; their chords are disjoint; it is the *midpoints*,
  not the centres, that are `(h_1 + h_2)/2` apart.

The correct pair statement is the pair of **endpoint** rows of (C0), and that is what is used
below.  (This is not a cosmetic repair: the false form would have made a chord row look like a
strengthening of the ordinary `y`-separation row `|y_i - y_j| >= m_ij`, which it is not.)

### 1.3 The three chord rows, with their leaf hypotheses  [proved]

Throughout, `delta >= 0` and the configuration satisfies the wall rows (so every square lies in
`[delta, T - delta]^2`) and, for each pair, one separating-axis row at level `delta`.

> **(C1) endpoint rows.**  On a leaf fixing the active branches, `alpha_i(L)` and `beta_i(L)`
> are affine in `c_i`.  *No inequality; a change of variables.*
>
> **(C2) chord-length rows.**  On a leaf fixing the active branch and the sign of `e_i`, the
> statement `h_i(L) >= H` is the single-variable box row
> `|e_i| <= p_i - H |C_i| |S_i|`, i.e. `x_i` in an interval of length `2(p_i - H C_i S_i)`
> around `a`.  In particular `h_i >= sec(tilt_i)` iff `|e_i| <= (C_i - S_i)/2`, and
> `h_i >= 1` iff `|e_i| <= D(th_i)`.
>
> **(C3) chord-chain rows.**  Let `L` be *any* line and `K` a set of `k` squares whose
> **interiors** meet `L`.  Then
>
>     sum_{i in K} h_i(L)  <=  len( L ∩ [delta, T-delta]^2 )  -  (k-1) delta ,
>
> and for `L` vertical or horizontal this is
>
>     sum_{i in K} h_i(L)  <=  T  -  (k+1) delta .

*Proof of (C3).*  Each `L ∩ int(Q_i)` is a relatively open subinterval of `L` of length `h_i`,
and they are pairwise disjoint because the interiors are.  All of them lie in
`L ∩ [delta, T-delta]^2` because each `Q_i ⊆ [delta, T-delta]^2` (wall rows).  For two of them,
with `Q_i` before `Q_j` along `L`, the separating row gives a unit normal `n` with
`Q_i ⊆ {n.(z - c_i) <= 1/2}` and `Q_j ⊆ {n.(z - c_i) >= 1/2 + delta}`; writing `n_L` for the
component of `n` along `L`'s direction, `n_L != 0` (if `n_L = 0` the slab is parallel to `L` and
at most one of the two interiors can meet `L`), and moving along `L` from the end of the first
chord to the start of the second costs `delta/|n_L| >= delta` since `|n_L| <= 1` and
`delta >= 0`.  So the `k` chords plus the `k-1` internal gaps fit inside
`L ∩ [delta, T-delta]^2`.  For `L` vertical, `len(L ∩ [delta, T-delta]^2) = T - 2 delta`. ∎

Three remarks, all load-bearing.

* **The interiors matter.**  Under the repo's closed semantics (`s13-casefree.md` §1) two
  squares may share an edge, and then their *closed* chords on that edge's line **coincide**.
  Counting closed chords makes (C3) false: on the axis-parallel `3 x 3` tiling the line `x = 1`
  carries four closed chords of length `1` in a container of length `3`.  With the interior
  convention the same line carries **none**.  This is the reason `t3-existence.md` §5.3's
  informal statement "`sum_i h_i <= T` over all squares meeting the line" has to be read with
  strict incidence, and it is what makes the chord clause fail at the pinwheel (§4).
* **`delta >= 0` is a hypothesis of (C3), not a conclusion.**  At `delta < 0` the gap bound is
  `delta/|n_L|`, which is *weaker*, and the row can be violated by up to `|delta|/(C S)` — an
  amplification of order `1/t` at small tilt.  Measured: at the `coh0.5` optimum
  (`delta* = -4.7e-06`) the axis-parallel chord sum overshoots `T` by `1.3e-03`, i.e. by `280x`
  the margin.  Any use of (C3) as a *detector* at `delta < 0` configurations must carry this
  factor; the tables in §4 report `max_c F` rather than a spurious "bound".
* **(C3) adds no new valid inequality to the centre LP.**  Its proof is a nonnegative
  combination of (i) the leaf's pair rows, (ii) the wall rows, and (iii) the single-variable
  box rows of (C2)/(C0) that fix the branches.  So the chord extension does not enlarge the
  cone of valid rows; what it enlarges is the **leaf system** — it forces branching on the
  *centres* as well as on the angles.  `S6_SKELETON.md` §4.5 already says exactly this
  ("`chord_i(a)` is concave in `(c_i)_y`, so putting it in an LP needs a lower bound that is
  linear in `(c_i)_y`, i.e. branching on the centres as well — which is the `3n`-dimensional
  problem the reduction of §3.1 was meant to avoid"), and this note confirms it: **the
  usefulness of chord rows is entirely a question of whether the centre-leaves they need are
  provably exhaustive.**  §4 says they are not.

---

## 2. Kearney–Shiu §3 as Farkas certificates over leaves

Read from the paper, not from the repo's summary; `runs/t3_chord_ks.txt` is the verification.

### 2.1 The objects, verified  [measured]

| object | statement | check |
|---|---|---|
| the green lattice (6) | `{(sqrt2-1/2,1), (3/2,1), (7/2-sqrt2,1), (3/2,3/2), (sqrt2-1/2,2), (3/2,2), (7/2-sqrt2,2)}` | `13` distinct points in green ∪ red, `1` C-point, `4` B-points at distance `1/2`, `8` A-points at distance `0.7702` — as the paper says |
| unavoidability | every unit square in `[0,3]^2` covers a green point | `min over squares of max_p inside(p) = +0.000000`, argmin the axis-parallel tiling square; **tight, zero slack** (`721` angles x `361^2` centres) |
| Lemma 1 | a unit square covering `C` covers a B-point | `min over squares covering C of max_B inside(B) = +0.000000`, argmin again axis-parallel; **tight** (`901 x 361^2`) |
| Lemma 2 | the two stated points lie on edges of the standard square | `max |inside| = 2.2e-12` over `2e5` values of `t` |
| (7) | `(1+t^2)/(1+t) >= 2 sqrt2 - 2`; `(1+2t-t^2)/2 >= 1/2`; sum `>= 3/2` | minima `0.828427125` at `t = sqrt2 - 1`, `0.5` at `t = 0`, `1.5` at `t = 0` — exact |
| Lemma 3, in the form the proof uses | `V` covers `C` and `(3/2,2)` and avoids `(1,3/2),(1,2),(2,3/2),(2,2)` ⇒ `V` meets `x = 1` and the bottom of its chord is `<= 5/3` | `100 %` of the admissible grid meets `x = 1`; max bottom `1.662235 <= 5/3 = 1.666667` (`1801` angles x `501^2` centres) |

Two things worth recording.  First, **the first inequality of (7) is the repo's own wall-strip
chord lemma**: `min_th (C + S - 1)/(C S) = 2 sqrt2 - 2` at `45 deg` is exactly `chord_lemma.md`
Lemma 1 evaluated at `d = 1 - p`, i.e. Stromquist's `f(th, 1)`, which is proved in Lean
(`lean/Sqpack/Chord.lean`).  So one of K–S's two omitted lemmas is already formalised in this
repo.  Second, **the Lemma-3 hypothesis is load-bearing**: without "`V` covers the green
B-point `(3/2,2)`", `12.5 %` of the admissible squares do not meet `x = 1` at all and the
`5/3` cap says nothing.  The paper states Lemma 3 only in its extremal form (three points *on*
edges) and applies it in the general form; the general form is what is measured above.

### 2.2 Configuration (b): a genuine Farkas certificate  [proved]

The last three lines of the paper.  In configuration (b) the pairings are forced, and on the
line `L : x = 1` there are three relevant squares:

* `W` covers the red A-point `(1, sqrt2 - 1/2)` — so `beta_W >= sqrt2 - 1/2`;
* `X` covers the green A-point `(sqrt2 - 1/2, 1)` and the red B-point `(1, 3/2)` — so by
  Lemma 2 + the first inequality of (7), `h_X = beta_X - alpha_X >= 2 sqrt2 - 2`;
* `V` covers `C` and the green B-point in region 2 — so by Lemma 3, `alpha_V <= 5/3`;

and `W`, `X`, `V` meet `L` in that vertical order.  The rows, in the centre variables
`(x_W, y_W, x_X, y_X, x_V, y_V)` with the endpoint forms of (C0) on the leaf:

| | row | `delta`-coefficient | kind |
|---|---|---|---|
| `r1` | `alpha_X - beta_W >= delta` | `1` | (C3), `k = 2` |
| `r2` | `beta_W >= sqrt2 - 1/2` | `0` | incidence leaf |
| `r3` | `alpha_V - beta_X >= delta` | `1` | (C3), `k = 2` |
| `r4` | `-alpha_V >= -5/3` | `0` | Lemma 3 leaf |
| `r5` | `beta_X - alpha_X >= 2 sqrt2 - 2` | `0` | (C2) on the Lemma-2 leaf |

> **Weights `w = (1, 1, 1, 1, 1)`.**  Every centre variable cancels **identically**, for *any*
> branch assignment — `a_X - b_X + (b_X - a_X) = 0` at `X`, `-b_W + b_W = 0` at `W`,
> `a_V - a_V = 0` at `V` — and
>
>     0  >=  2 delta + (sqrt2 - 1/2) - 5/3 + (2 sqrt2 - 2) ,
>     delta  <=  (13/6 + 2 - 3 sqrt2)/2  =  -0.037987010226 .

Verified numerically: residual `max |w.A| = 0.000e+00` exactly, value agreeing with the closed
form to `2.2e-16` (`python3 search/t3_chord_ks.py certb`).  K–S's own display
`13/6 - sqrt2 = 0.752453 < 2 sqrt2 - 2 = 0.828427` is the same statement with the `2 delta`
suppressed; the gap `0.075974` is twice the certificate value.

So `t3-chain.md` §4.2's dictionary entry — "their final contradiction is exactly a chain
inequality on the line `x = 1`, with `5/3` playing the role of the wall row" — is **right in
shape and wrong in detail**: it is a chord chain, not a separation chain, and neither of the
two end rows is a wall row.  Both ends are chord-disjointness rows against *other squares*
(`W` below, `V` above); the container never enters.  That matters for `T`-transfer (§6): the certificate has no wall row and no `T` in its
*structure*, and two of its three constants (`sqrt2 - 1/2`, `2 sqrt2 - 2`) are `T`-free
one-square constants; only `5/3` (Lemma 3) knows `T = 3`.

### 2.3 What resists: everything else  [proved]

The rest of §3 of the paper is not a certificate and cannot be made into one.

* **Unavoidability of (6)** is an *existence clause* about a single square — "every unit square
  in `[0,3]^2` covers one of these seven points".  It plays exactly the role (E3) plays in this
  repo, and it is not a row: it is a disjunction over incidences, i.e. a leaf-generation
  device.  It is also **exactly tight** (§2.1), so it has no slack to spend.
* **The counting** — `6` squares, `7` points, each square covering `>= 1` of each lattice, so
  at most one square covers two of one lattice — is combinatorics on incidences.  Nothing in
  the centre-variable row language expresses "these two squares cover different points".
* **Configuration (a) ends in that counting, not in `delta <= 0`.**  Its numeric input is the
  *third* inequality of (7) (`sum >= 3/2`); its output is "the square covering the red A-point
  in region 9 must have upper intercept on `x = 2` at least `3/2`", i.e. it must also cover the
  red B-point in region 6 — which contradicts the bijection.  As a row statement this is
  `beta >= 3/2` and `alpha <= sqrt2 - 1/2`, hence `(2, 3/2) in` the chord: a perfectly good
  *linear consequence*, but the contradiction it feeds is the counting, which has no Farkas
  form.
* **Every intermediate step of (b)** (the two "if `W` covers `(1,1)`" branches) is of the same
  type.
* **Two steps of the paper are asserted, not proved**: the proofs of Lemma 2 and Lemma 3 are
  omitted as "elementary coordinate geometry", and the sentence *"Elementary geometric
  considerations show that any other positioning of the second square will increase the
  intercept `3/2`"* is the load-bearing generalisation of Lemma 2 out of its extremal position.
  §2.1 measures both; neither is in doubt, but neither is in the paper.

> **Step S, named.**  The only part of Kearney–Shiu that is a certificate in the repo's row
> language is the five-row combination of §2.2.  **Everything that produces the leaves — the
> unavoidable set, its `90 deg` dual, and the `6`-squares-vs-`7`-points bijection — is the
> existence side, and it has no row-language substitute.**  This is precisely the calibration
> the brief asked for, and it says the same thing `t3-existence.md` §4.4 says: the certificate
> side was never the problem.

---

## 3. What chord rows buy on the certificate side

### 3.1 The banded column of three  [proved]

> **Proposition CC3.**  Let `L` be a vertical (or horizontal) line and let `i, j, k` be three
> squares whose centres satisfy `|x_i - a| <= (C_i - S_i)/2` (tilts in `[0, 45 deg]`, so the
> band is nonempty).  Then each chord is `h = 1/max(C, S) = sec(tilt)`, the interiors all meet
> `L`, and (C3) gives
>
>     sec(t_i) + sec(t_j) + sec(t_k)  <=  T - 4 delta ,
>
> so at `T = 3` and a common tilt `t`,
>
>     delta  <=  -3 (sec t - 1)/4  =  -(3/8) t^2 + O(t^4)  =  -0.375 t^2 + O(t^4) ,
>
> with equality `0` exactly at `t = 0`.

This is the first row in the repo that is **simultaneously** exactly tight in the axis-parallel
stratum and correct to second order off it.  Compare (`t3-existence.md` §4.2, §7): the bare
chain of three plus its optimal wall completion gives `[u(T-u) - (T-1)]/[(T-1) + 2u]
= +t/4 + O(t^2)`, first order and the *wrong sign*; Lemma H is exactly critical at first order;
Lemma Z'' is discontinuous at `t = 0`.  Chord rows are the first mechanism that is not.

They are also strictly stronger than Lemma Z'' where they apply: Lemma Z'' needs a link normal
that is **exactly** axial, i.e. an exactly axis-parallel square; CC3 needs no exact incidence at
all, only a band condition on three centres.

### 3.2 But they are not new inequalities

By §1.3 remark 3, (C1)–(C3) live in the cone generated by the existing wall and pair rows
together with single-variable centre-box rows.  So the gain of §3.1 is bought entirely by the
new leaf — "three centres in a band of width `C - S` around one line" — and the whole question
is whether that leaf is forced.  That is §4.

---

## 4. The pigeonhole substitute, and where it fails

### 4.1 The two pigeonholes that do work  [proved]

> **(P-count)** `∫_0^T #{i : int(Q_i) meets the line x = c} dc = sum_i u_i >= n`, with equality
> iff every square is axis-parallel.  At `n = 6`, `T = 3` this gives average count `>= 2`, and
> **strictly `> 2` as soon as one square is tilted**, so some vertical line meets `>= 3`
> squares.  (The brief's item 2 is right about this.)
>
> **(P-area)** `∫_0^T F(c) dc = sum_i area(Q_i) = n = 6 = 2T`, where `F(c) = sum_i h_i(c)`.  So
> `mean_c F = 2` and the certificate's target `F(c) >= T = 3` is `50 %` above the mean: the
> gain has to come entirely from *choosing* `c`, exactly as the brief anticipates.

Both are `T`-uniform and both are proved.  Neither is enough.

### 4.2 The clause is false at a `delta = 0` packing  [proved]

> **Counterexample.**  Take the **pinwheel**: six axis-parallel unit squares at the cells of the
> `3 x 3` grid minus a permutation, e.g. centres
> `(0.5,1.5), (0.5,2.5), (1.5,0.5), (1.5,2.5), (2.5,0.5), (2.5,1.5)`.  This is a packing with
> `delta = 0`.  Then **for every line of every direction**, the chord functional is at most the
> container chord, with `max_c F = 2` on every axis-parallel line (exact scan over the
> breakpoints, `search/t3_chord_scan.py`), and `max_L [sum h - len(L ∩ [0,3]^2)] = 0`, the
> maximum being attained at a corner-clipping line with `k = 2` (`search/t3_chord_oblique.py`,
> `721` directions x `3000` offsets).  So no line of any direction carries a chord sum
> exceeding what the container allows, and no chord certificate of positive strength exists at
> this packing.

So the chord existence clause

> **(E3-chord)** every packing of six unit squares in `[0,3]^2` with `delta >= 0` has an
> axis-parallel line met by three squares whose chords sum to `>= T`

is **false**.  Its content at `t = 0` is exactly "three of the squares share a common
transversal", i.e. three `x`-coordinates (or three `y`-coordinates) inside an **open** window of
length `1`; the pinwheel's three columns hold two squares each and its `x`-coordinates
`{0.5, 0.5, 1.5, 1.5, 2.5, 2.5}` have no three inside an open window of length `1` (the closest,
`{0.5, 0.5, 1.5}`, spans exactly `1`), and the same in `y`.  What it *does* have is a
**staircase** chain of three — precisely the case Theorem AP of `t3-existence.md` §4.2 has to
allow (`L_x = L_y = 3` in the staircase sense), and precisely the case Lemma Z'' certifies and
the chord row does not.  **(E3-chord) is strictly weaker than (E3-row), and the gap is the
staircase.**

Measured over the recorded packings (`search/t3_chord_scan.py packings`, exact over the
breakpoints, `runs/t3_chord_scan.txt`):

| `#` exactly axis-parallel | packings | (E3-chord) holds | `max_c F` |
|---|---|---|---|
| `3` | `255` | `255` | `3.000000` |
| `4` | `358` | `358` | `3.000000` |
| `5` | `109` | `109` | `3.000000` |
| `6` | `27` | **`26`** | `3.000000`, one at `2.000000` |
| total | `749` | `748` | median `3`, min `2` |

`748/749`, and the single failure is the pinwheel.  Every firing case fires at value **exactly
`0`** — the chord chain is exactly tight, never strict, at a margin-`0` packing.  The
corner-clipping oblique lines are tight too, but their tightness is the exact-tiling
coincidence "two squares abut one another and the wall so as to tile the triangle `L` cuts
off"; at `delta = 1/3` (four axis-parallel squares) the same scan returns `-1.000000`, so
nothing there survives a positive margin, and the hypothesis is an exact incidence of the
Kearney–Shiu kind, not a pigeonhole consequence.

### 4.3 Step S: the pigeonhole cannot produce the band  [proved]

The certificate of §3.1 needs three centres within `(C - S)/2` of one line, i.e. three of the
six `x`-coordinates (or `y`-coordinates) inside a **window of length `C - S`**.  The centres
lie in `[p, T-p]`, an interval of length `T - u`, `u = C + S`.  Six points in an interval of
length `l` always have three inside a window of length `l/2`, and `l/2` is **sharp** (the
`2-2-2` placement realises it).  So:

> **Proposition PG.**  The chord pigeonhole produces the leaf of Proposition CC3 iff
>
>     (T - u)/2  <=  C - S       <=>       3 cos t - sin t  >=  3       (at T = 3),
>
> which holds **only at `t = 0`**, where it holds with equality.  For `t > 0` the deficit is
>
>     (T - u)/2 - (C - S)  =  (T - 2)/(T - 1) · t  +  O(t^2)      (general T)
>                         =  t/2 + O(t^2)                         (T = 3) .

At general `T` the window needed is still `C - S` (the chord is `sec t` there and `T` of them
overflow `T`), the interval is `T - u` long, and pigeonhole splits it into `T - 1` windows of
length `(T-u)/(T-1)`, one of which holds `>= T` of the `T^2 - T` centres; the criterion
`(T-u)/(T-1) <= C - S` is again an equality at `t = 0` and fails at rate `(T-2)/(T-1)` per
radian.  **The pigeonhole substitute is exactly critical at the tiling and first-order wrong
everywhere else, at every `T`.**  This is step **S**.

### 4.4 The numbers  [measured]

`python3 search/t3_chord_scan.py band` at the uniform-tilt optima of
`runs/t3_chain_scan1.jsonl`:

| `t` (deg) | needed `C - S` | pigeonhole `(T-u)/2` | pigeonhole deficit | smallest realised `3`-window | its excess | `delta*` |
|---|---|---|---|---|---|---|
| `0.5` | `0.991235` | `0.995656` | `0.004420` | — | — | — |
| `1` | `0.982395` | `0.991350` | `0.008955` | `0.982586` | `1.91e-04` | `-1.11e-04` |
| `5` | `0.909039` | `0.958325` | `0.049286` | `0.913859` | `4.82e-03` | `-2.45e-03` |
| `10` | `0.811160` | `0.920772` | `0.109612` | `0.834714` | `2.36e-02` | `-8.39e-03` |
| `20` | `0.597672` | `0.859144` | `0.261471` | `0.704559` | `1.07e-01` | `-2.44e-02` |
| `30` | `0.366025` | `0.816987` | `0.450962` | `0.571649` | `2.06e-01` | `-3.89e-02` |
| `45` | `0.000000` | `0.792893` | `0.792893` | `0.672954` | `6.73e-01` | `-4.83e-02` |

The two right-hand columns are the whole story: the **truth** (the realised window is too long
by `1.91e-04 ~ 0.63 t^2` at `1 deg`) is second order and commensurate with `delta*`; the
**argument** (the pigeonhole is too weak by `8.96e-03 ~ t/2`) is first order and `47x` larger.

And the chord functional itself, `python3 search/t3_chord_scan.py unif`:

| `t` (deg) | `delta*` | `T - max_c F` | `(T - max_c F)/t` |
|---|---|---|---|
| `1` | `-1.108e-04` | `1.047e-02` | `0.600` |
| `5` | `-2.450e-03` | `4.406e-02` | `0.505` |
| `10` | `-8.392e-03` | `8.294e-02` | `0.475` |
| `20` | `-2.437e-02` | `1.400e-01` | `0.401` |

so `T - max_c F = 0.6 t + O(t^2)`: **the chord sum at the uniform-tilt optimum falls short of
`T` by a first-order amount**, one order flatter than `delta*`.  The mechanism of the
amplification is visible in (C0): placing the line at the middle of the best `3`-window leaves
the two outer centres `9.6e-05` outside the band at `t = 1 deg`, and each loses
`(|e| - (C-S)/2)/(C S) = 9.6e-05/0.017449 = 5.5e-03` of chord — `O(t^2)/O(t) = O(t)` — for a
total of `1.10e-02` against the measured `1.05e-02`.  **A second-order displacement of the
centres becomes a first-order loss of chord, because the chord's slope in the centre is
`1/(C S) ~ 1/t`.**  That is the whole reason the chord route fails where its rows succeed.

---

## 5. Obstruction X, with chord rows

X (`t3-existence.md` §4.4) says: the existence clause must separate packings from
configurations of margin `-0.36 t^2` as `t -> 0`, and every mechanism stable under an `O(t)`
perturbation of the angles is blind at that scale.  Chord rows change the picture in one
direction and make it worse in another.

**Better:** the *certificate* is no longer blind.  Proposition CC3 gives `-0.375 t^2`, which is
second order and of the right sign; the counting/pigeonhole/unavoidable-point mechanisms all
gave `O(t)` of the wrong sign.  So the brief's hope — that chord rows are "first order in the
tilt off the axis-parallel stratum" and therefore reach where nothing stable reaches — is
correct **about the rows**.

**Worse:** the *leaf* is now a second-order statement with almost no room in it.

> **The `3 %` observation.**  `delta*(t) = -0.36371 t^2 + O(t^3)` (`runs/t3_chain_scan1.jsonl`,
> `-1.107934e-04` at `t = 1 deg`), while CC3 would certify `-0.375 t^2`.  A valid certificate
> can never give a bound **below** `delta*`, so the leaf of CC3 **must** be empty at the
> uniform-tilt optimum — and it is, by `1.9e-04 = 0.63 t^2` of window length.  The two
> second-order constants differ by only `0.0113 t^2`, i.e. by `3.1 %`: the row is *just* too
> strong to be available, and the configuration sits *just* outside its leaf.  **[measured]**

That is the honest statement of where the route ends.  The certificate constant and the truth
agree to `3 %`, which is the good news; the price is that the leaf is missed by `0.63 t^2`, a
second-order amount, so **any proof that the leaf is non-empty at a packing must be accurate to
`O(t^2)` in the centre positions** — which is obstruction X restated, not reduced.  The
pigeonhole available is accurate only to `O(t)` (§4.3), one whole order too coarse, and that is
the same order by which chain counting misses (`t3-existence.md` §4.4).  Chord rows move the
*certificate* from `O(t)`-wrong to `O(t^2)`-right and leave the *existence* side exactly where
it was.

**No configuration family had to be constructed for the brief's item 3**, because the one the
repo already has does the job: the uniform-tilt family `unif t` of `runs/t3_chain_scan1.jsonl`.
At `t = 1 deg` its optimum is (centres, `theta = 1 deg` throughout)

    (2.491167, 1.508576)  (1.491430, 1.491125)  (1.508276, 2.491461)
    (2.491461, 0.508539)  (0.508539, 0.508539)  (0.508539, 2.474010)

— the `3 x 3` tiling minus the transposition `{(0,1), (1,0), (2,2)}`, re-optimised.  Every
chord certificate's leaf fails there: the smallest `3`-window in `x` is `0.982586` and in `y`
is `0.982586`, both above `C - S = 0.982395`; `max_c F = 2.989533 < 3`; and
`delta* = -1.107934e-04`, i.e. `1.1e-04` from being a packing.  `delta*` here is the
`t3_chain` instrument's value (multistart, `12` samples, all agreeing); the `s6skel` `Decider`
does not terminate at a generic exact angle vector inside its node cap (`S6_SKELETON.md` §3.2),
so no `Decider` kill is available at these angles and none is claimed.

---

## 6. Transfer to `T = 4`

### 6.1 What is `T`-uniform and what is not

| ingredient | `T`-uniform? |
|---|---|
| the chord formula (C0), the bands `(C-S)/2` and `D(th)` | **yes**, no `T` in them |
| the chord-chain row (C3) | **yes**: `sum h <= T - (k+1) delta` |
| Proposition CC3 (a banded column of `T`) | **yes**: `T sec t <= T - (T+1) delta`, i.e. `delta <= -T(sec t - 1)/(T+1) = -(T/(2(T+1))) t^2` |
| the K–S certificate of §2.2 | **yes and empty**: it contains no `T` at all, but its leaf is produced by the counting, which is not `T`-uniform |
| (P-count), (P-area) | **yes**: `sum u_i >= n`, `∫F = n`, so `mean_c F = n/T`; at `n = T^2 - T` this is `T - 1`, and the target is `T` — the **same** `1`-unit deficit at every `T` |
| the pigeonhole of §4.3 | **yes in form, empty in content**: exactly critical at `t = 0`, deficit `(T-2)t/(T-1)` |
| the counting `n` squares vs `n+1` points | **no** — §6.2 |

The `mean_c F = n/T = T - 1` line is worth isolating: **at every `T`, the chord functional has
to be pushed from its mean `T - 1` to `T`, i.e. by exactly one unit, and the mean is fixed by
area alone.**  That is the same "one unit of slack" as `proof-anatomy.md` §7.2 and K–S's
`6` squares vs `7` points, arrived at from a completely different direction.

### 6.2 (a) Is a `13`-point pure unavoidable set for `[0,4]^2` excluded?

**No, and it cannot be excluded by the repo's method.**  The relevant facts, all already in the
repo:

* `search/COVER4.md` proves, exactly, `COVER^closed(4) >= nu_f^closed(4) >= 12.2688038611`.  A
  pure point set is a cover of integer weight, so **every pure unavoidable set of `[0,4]^2` has
  at least `13` points.**  In particular the `11`-point set that `t3-chain.md` §4.4 names as
  what "`s(12) >= 4` would need" (the slack-`0` route: an unavoidable set of `N-1 = 11` points
  gives `s(12) >= 4` in one line) is **excluded outright** — as is `12`.  That is a proved
  statement the repo has not drawn.
* So `13` is the *smallest conceivable* size, and it is exactly the size the K–S slack-`1`
  argument needs at `n = 12` (`13 - 12 = 1`, the same slack as `7 - 6 = 1` at `n = 6`).  The
  brief's count of `13` is right and `t3-chain.md` §4.4's `11` is the count for a different
  (slack-`0`) argument.
* Nothing in the repo excludes `13`.  The best known pure set is Friedman's `14` (DS7 Thm 4,
  `notes/s13-casefree.md` §1), and the disjoint-square lower bound is useless here: `s(13) = 4`
  forces at most `12` pairwise disjoint *closed* unit squares in `[0,4]^2`, and `s(12) = 4`
  would force at most `11` (`proof-anatomy.md` §7.1).
* **And the LP route can never decide it.**  The repo's own estimate of the fractional optimum
  is `COVER^closed(4) ~ 12.39-12.46` (`CLOSED4.md`, `FAMILY.md`), and the exact certified
  weighted cover of `s13-casefree.md` has weight `12.955972155 < 13`.  A fractional lower bound
  can exclude `13` only by exceeding `13`, and the fractional optimum is believed to be below
  `12.5`.  **The question "is there a `13`-point pure unavoidable set for `[0,4]^2`?" is an
  integrality question — the certified fractional floor is `0.73` below `13` and the fractional
  optimum itself is believed to be `~0.5` below it — and no bound of the kind the repo produces
  can answer it.**

> **What it would take to decide.**  (i) *Yes*: exhibit `13` points and verify unavoidability
> exactly.  This is a search in `26` real dimensions with an exact verifier already in the repo
> (`dual_exact.py` / `verify2/`, the same machinery that certifies the `3621`-point weighted
> cover); the natural attack is to take the `12.96`-weight certified weighted cover of
> `certificates/rung2/s13_closed_cover_4.txt` (`3621` points) and look for an integral
> `13`-point cover on or near its support — an integer program with a `0.04`-unit budget,
> tight but not obviously infeasible.  **[guess]**  (ii) *No*: an integrality obstruction, e.g. a
> family of `13` unit squares in `[0,4]^2`, pairwise disjoint **as closed sets except at a
> controlled set of contacts**, forcing two of them to share a point of any `13`-point cover.
> Nothing of this shape exists in the literature.  Status: **open, decidable in principle by
> (i), not decidable by any LP bound.**

For completeness: even a `13`-point set would only start the K–S argument at `T = 4`.  The
`90 deg` duality survives (`t3-chain.md` §4.4), Lemma 1 and Lemma 2 + (7) survive (they are
one-square statements with no `T`), and Lemma 3's `5/3` is the one constant that is
container-specific and would have to be redone.  The certificate of §2.2 would transfer in
form, with `5/3` replaced by the `T = 4` analogue of Lemma 3.  **The whole of the `T = 4` K–S
route hangs on the `13` points, and on nothing else that is currently unknown.**

### 6.3 (b) Does the chord pigeonhole give anything at `T = 4`?

**No — nothing at all, and in particular nothing where the centre pigeonhole of
`t3-existence.md` §3 is silent (which is everywhere).**

* *Small tilt.*  §4.3 at `T = 4`: `12` centres in an interval of length `4 - u`, split into `3`
  windows of length `(4-u)/3`, one holding `>= 4`; the certificate needs a window of length
  `C - S`; `(4-u)/3 <= C - S` iff `4 cos t - 2 sin t >= 4`, i.e. iff `t = 0`, with deficit
  `2t/3 + O(t^2)`.  Exactly the `T = 3` failure, one third worse.
* *At `t = 0` exactly.*  The clause needs a **straight** column of `4`, and the `4 x 4` tiling
  minus a permutation matrix (`12` squares, `3` per row and `3` per column) has none — the same
  pinwheel obstruction as at `T = 3`, and at `T = 4` there are `4! = 24` such patterns rather
  than `3! = 6`.  So (E3-chord) is false at `T = 4` too, and more often.
* *Far field.*  At `45 deg` the band `(C-S)/2` is empty, so axis-parallel chord rows say
  nothing.  Taking the line along an edge normal instead makes every chord `1` but stretches
  the container chord to `T sec t`, so the row reads `k <= T sec t = 5.657` and needs
  `k >= 6` squares inside a **unit** window of the projection; the `12` projections spread over
  an interval of length `(4 - u)u = 4 sqrt2 - 2 = 3.657`, which splits into `4` windows of
  `0.914 < 1`, so pigeonhole guarantees only `k >= 3`.  Nothing.
* The centre pigeonhole of `t3-existence.md` §3 never fires at `T = 4` (`4 - sqrt2 = 2.5858 >
  1/d_12 = 2.5725`, by `0.5 %`).  The chord pigeonhole never fires at `T = 4` either, and for a
  different reason: it is exactly critical at `t = 0` and the tiling realises the critical case.
  **The two mechanisms are empty in disjoint ways, and their union is still empty.**

---

## 7. What in the brief's premises, and in the repo, turned out to be wrong

* **The brief's row is wrong twice** (§1.2).  The chord of a tilted square is *not* centred at
  the square's centre — the offset reaches `0.4728` — and consequently
  `|y_i - y_j| >= (h_i + h_j)/2` is **false** for disjoint squares (explicit counterexample,
  violation `0.2576`).  The correct object is the pair of endpoint rows `alpha_i, beta_i`,
  which are affine in the centre on a branch leaf.  The brief's remark that `h_i` is concave
  piecewise-linear in `|x_i - c|` is right; what it misses is that the *position* of the chord
  is piecewise-linear too, and not equal to `y_i`.
* **`t3-existence.md` §5.3's statement of the chord row is wrong under the repo's own
  semantics.**  "`sum_i h_i <= T` over all squares meeting the line" is **false** for closed
  semantics: on the axis-parallel `3 x 3` tiling the line `x = 1` has four squares meeting it
  with closed chords of length `1` each, total `4 > 3`.  The row is true only for squares whose
  **interiors** meet the line — and with that correction the same line carries *no* chord at
  all.  This is not a technicality: it is exactly what makes (E3-chord) fail at the pinwheel.
* **`t3-existence.md` §5.3's "they are *not* implied by any single separating-axis row — their
  validity uses the disjunction, not one leaf of it" is wrong in the direction that matters.**
  On a leaf, a chord row *is* a nonnegative combination of that leaf's pair rows, the wall rows,
  and single-variable centre-box rows (§1.3).  Chord rows add no valid inequality; they add
  **centre branching**.  `S6_SKELETON.md` §4.5 had this right and `t3-existence.md` §5.3 lost
  it.
* **`t3-chain.md` §4.2's dictionary entry for the final contradiction is right in shape, wrong
  in detail** (§2.2): `5/3` is not "playing the role of the wall row" — both ends of the chain
  are chord-disjointness rows against other squares, and the container does not appear.  The
  certificate is `T`-free.
* **`t3-chain.md` §4.4's "`s(12) >= 4` would need an unavoidable set of `11` points" is the
  count for the *slack-`0`* argument, and that count is now excluded by the repo's own
  `COVER^closed(4) >= 12.2688`** (§6.2).  The K–S argument proper needs `13`, which is the
  smallest number the repo does not exclude.  Also, `t3-chain.md` §4.4's "the counting `6`
  squares vs `7` points, slack `1`: **dies** at `T = 4`" is conditional, not absolute — it dies
  *against the largest known* set (`14`, slack `2`); against a `13`-point set the slack is
  exactly `1` and the structure transfers.
* **The brief's "any tilt makes `sum_i (|cos| + |sin|) > n = 6` ... so some vertical or
  horizontal line meets `>= 3` squares with total chord `<= 3`" is right about the count and
  useless about the chord.**  The count is proved ((P-count), §4.1); the chord sum at such a
  line is `2` at the tiling and `3 - 0.6 t` at the uniform-tilt optimum, never `>= 3`.
* **`t3-existence.md` §5.3's forecast — "this extension is invisible in the axis-parallel
  stratum and first order in the tilt off it — i.e. it is a candidate for exactly the region
  §4.4 says nothing stable can reach" — is half right.**  It is invisible in the axis-parallel
  stratum (worse than that: it *fails* there, at the pinwheel).  Off it, it is **second** order,
  not first, on the certificate side — better than forecast — and the existence side is first
  order, which is what kills it.
* **Two of Kearney–Shiu's own steps are unproved in the paper** (Lemma 2, Lemma 3, and the
  "elementary geometric considerations" generalising Lemma 2).  All three are confirmed
  numerically here, and the first inequality of (7) turns out to be the repo's already-formalised
  wall-strip chord lemma.

---

## 8. Verdict, restated

> **`s(6) = 3` is not proved in the row language, with or without chord rows.  Kearney–Shiu is
> reproduced only as far as it is a certificate at all — its final five-row combination, weights
> `(1,1,1,1,1)`, value `delta <= (13/6 + 2 - 3 sqrt2)/2 = -0.037987010226`, residual `0` exactly
> — and the rest of their §3 is unavoidable-set counting, which has no row-language form.  The
> pigeonhole substitute fails at step S = producing the certificate's leaf: the leaf is "three
> centres in a window of length `C - S`", the sharp pigeonhole gives a window of length
> `(T - u)/2`, and `(T-u)/2 <= C - S` holds only at `t = 0`, with deficit `(T-2)t/(T-1)` — first
> order in the tilt, at every `T`.  Obstruction X is not reduced; it is sharpened, from
> `O(t^2)` to `0.03 · O(t^2)`.  The chord extension is a genuine improvement on the certificate
> side (`-0.375 t^2`, exactly tight at `t = 0`) and a regression on the existence side (false at
> the pinwheel, where (E3-row) holds).  Since the route cannot close `s(6) = 3`, by Evan's rule
> it is closed for `s(12) = 4`.**

### 8.1 The one thing that survives, for the record

Proposition CC3 is a valid, `T`-uniform, second-order-correct certificate that the repo did not
have, and it is cheap: three centres in a band around one line.  If some *other* route ever
supplies the incidence "three squares of a packing share a transversal in band" — a rigidity
statement, in the sense of `S6_LOCAL.md` §5, not a combinatorial one — then CC3 closes the
tilted case in one line at every `T`.  The row is not the missing piece; the incidence is, and
it was the missing piece before chord rows too.

### 8.2 Reproduce

    python3 search/t3_chord_geom.py selftest      > runs/t3_chord_geom.txt     # ~2 min
    python3 search/t3_chord_ks.py all             > runs/t3_chord_ks.txt       # ~25 min
    python3 search/t3_chord_scan.py packings                                   # the 749
    python3 search/t3_chord_scan.py cone                                       # coherent cone
    python3 search/t3_chord_scan.py unif                                       # the X table
    python3 search/t3_chord_scan.py band                                       # step S
    python3 search/t3_chord_oblique.py pinwheel   > runs/t3_chord_oblique.txt  # ~8 min
    ( the four t3_chord_scan commands, concatenated, are runs/t3_chord_scan.txt )

Library entry points: `t3_chord_geom.chord(x, y, th, a, axis)` (exact endpoints and the active
branches), `t3_chord_geom._affine_endpoints(th, a, axis, bl, bh)` (the leaf-linear forms),
`t3_chord_scan.max_F(z, n)` and `t3_chord_scan.best_line_bound(z, n)` (the chord functional,
exact over the breakpoints), `t3_chord_ks.cmd_certb()` (the K–S certificate),
`t3_chord_oblique.scan(z, n)` (all directions).
