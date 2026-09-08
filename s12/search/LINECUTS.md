# The chord lemma on every line, as a family of packing cuts

Code: `search/linecuts.py` (the family, the predicate, the screens), `search/t4leaf.py`
(`--lines ...`, a flag: with it off `t4leaf.py` is bit-for-bit what `search/T4LEAF.md` ran).
Runs live in this worktree's `runs/` (gitignored), so every number is quoted here.
Nothing in `verify/`, `lean/`, `certificates/`, `xcheck.py` or `branch.py` is touched.

Read against `notes/chord-lemma.md` (the wall-strip lemma this generalises), `search/T4LEAF.md`
(the instrument and the three instances), `notes/branch-semantics.md` §1–2 (what a region
inequality means) and `search/T4SCREEN.md` §1 (the LP).

---

## 0. Verdict, up front

**The general-line chord inequality is a valid cut family, it is worth exactly `0.000000` on all
three instances, and there is a structural reason it cannot be worth much: the whole axis-parallel
family, aggregated optimally, says only `mu(axis-parallel poses) <= 9` at `t = 4`, and the optima
already sit at `6.7 – 8.0` on that quantity.**

| instance | with line cuts | matched, without them | line cuts worth | certified? |
|---|---|---|---|---|
| (a) `k = 4` corner leaf, cliques (`L1B40K` st0) | `11.992612` | `11.992612` | **`0.000000`** | `M = 1.0370` no |
| (a) the same, best stage (`L1B40K` st1) | **`11.995221`** | not valued | — | `M = 1.000000000`, `kmax = 1.000000` **yes** |
| (a') `k = 4`, **no cliques** — the smeared grid (st0) | `12.000000` | `12.000000` | **`0.000000`** | no |
| (b) **pure relaxation, no branch**, cliques (`L2PURE` st1) | `12.153620` | `12.153620` | **`0.000000`** | `M = 1.0377` no |
| (c) leaf `01010101`, cliques, best stage (`L3A01010101` st2) | `11.773816` | `11.773816` | **`0.000000`** | `M = 1.000000000`, `kmax = 1.0097` no |
| (c) the same, best *certified* stage (`L3A01010101` st0) | **`11.742058`** | `11.742058` | **`0.000000`** | `M = 1.000000000`, `kmax = 1.000000` **yes** |

(§5 has the trajectories, the matched `2 x 2`, and the convergence flags.  The two bold rows are
genuine lower bounds on their instance's value over its pose set; `L3A01010101` st0 is a slightly
*better*-certified leaf bound than `T4LEAF.md`'s `11.745193`, which still had an anchor clique
violated by `8.6e-4`, and `L1B40K` st1 improves on `T4LEAF.md`'s certified `11.988925`.)
**Nothing goes below 12 that was not already below 12, and (b) — "pure + line cuts `< 12` would
mean no branching tree is needed at all" — fails decisively: the certified `dual_PC1` measure of
mass `12.163060670` satisfies every one of `8 002` line cuts at pitch `0.001` with `0.05` to
spare.**

Four things are established, in decreasing order of how much they matter:

1. **The general-line lemma is true and its proof is the wall proof verbatim** (§1).  The
   short-line version is `N <= ceil(Lam) - 1` where `Lam` is the length of the line's intersection
   with the container — so a `45 deg` line near a corner gives `<= 1` or `<= 2`.  §1.3.
2. **The whole axis-parallel family has an aggregate ceiling of exactly 12** and, in its sharpest
   weighting, says `mu(axis-parallel) <= 9` (§3).  The three lines `y = 1, 2, 3` realise that
   sharpest weighting — this is exactly the observation the task started from, and it is the *best*
   the family can do, not a first step.  §5.1's positive control shows the LP's own dual finding
   that weighting unaided and landing on `9.000000`, so the bound is attained and not merely an
   upper estimate; and `9 + mu(tilted)` is all the family can ever say about the total.
3. **Every existing optimum already satisfies every line cut**, on a screen of `8 002` lines at
   pitch `0.001` (§4).  The largest excess anywhere is `+2.5e-8`, i.e. LP tolerance.  So the LP
   value with the cuts equals the LP value without them *exactly*, by the trivial argument: the old
   optimum is feasible for the new model.
4. **The reason is tilt, and it is quantitative.**  The extremal measures carry `4.6 – 7.2` of
   *non*-axis-parallel mass out of `12`, and a `45 deg` square's catch band has half-width
   `D = 0.2071` against `0.5` for an axis-parallel one.  The `4 x 4` grid *would* be cut hard — the
   sanity check of §2.3 shows `mu(B_L) = 4` on `y = 0.9` and `8` on `y = 1` against a cap of `3` —
   but the LP optimum is not the grid: it is the grid **smeared and tilted**, and the smearing is
   exactly what buys it room under the line cuts (`T4LEAF.md` §4).

---

## 1. The general-line chord inequality

### 1.1 The statement

> **Theorem L (line cut).**  Let `K` be a convex subset of the plane, `L` a line, and
> `Lam = length(L ∩ K)`.  Let `S_1, …, S_N ⊆ K` be **pairwise disjoint closed** sets, each of which
> meets `L` in a segment of length `>= 1`.  Then
>
> ```
> N < Lam ,        equivalently        N <= ceil(Lam) - 1 .
> ```

> **Theorem L' (packing form).**  Let `L > 1` and let `S_1, …, S_N` be squares of side `L` inside
> `K` with pairwise disjoint **interiors**.  Then at most `ceil(Lam) - 1` of their concentric
> **unit** squares meet a given line in a chord of length `>= 1`.

For the container `K = [0,t]^2` with `t <= 4` and `L` a horizontal or vertical line, `Lam = t <= 4`
and the bound is `3`.  That is `notes/chord-lemma.md`'s `wall_strip_le_three`, with the line's
height free instead of pinned at `9/10`.

### 1.2 The proof

*Proof of L.*  `L ∩ K` is a segment (an intersection of a convex set with a line); identify `L`
with `ℝ` by arclength so that `L ∩ K = [0, Lam]`.  Each `I_i = S_i ∩ L` is a closed interval
`[a_i, b_i] ⊆ [0, Lam]` with `b_i - a_i >= 1`, and the `I_i` are pairwise disjoint because the
`S_i` are.

The left endpoints are distinct (two intervals with the same left endpoint meet there), so relabel
`a_1 < a_2 < … < a_N`.  For each `i`, `a_{i+1} ∈ I_{i+1}` so `a_{i+1} ∉ I_i`; since
`a_{i+1} > a_i`, this forces `a_{i+1} > b_i >= a_i + 1`.  Chaining,

```
Lam  >=  b_N  >=  a_N + 1  >  a_{N-1} + 2  >  …  >  a_1 + N  >=  N ,
```

so `N < Lam`. ∎

This is `notes/chord-lemma.md` §4's argument with `t` replaced by `Lam`; there the last step reads
"four points of `[0,3]` pairwise more than `1` apart are impossible", which is `N < 4` for
`Lam = 4`.  **Nothing in it uses that the line is near a wall, is horizontal, or that `Lam = t`.**

*Proof of L'.*  `unit_subset_interior` (`lean/Sqpack/Basic.lean`): for `L > 1` the concentric
closed unit square is contained in the open `L`-square, so disjoint interiors of the `L`-squares
give **pairwise disjoint closed** unit squares, all inside `K`.  Apply L to them. ∎

Note that L is about *closed* disjointness and is false for interior-disjointness at `Lam = 4` —
the four axis-parallel unit squares at `(1/2 + i, 1/2)` in `[0,4]^2` meet `y = 1/2` in four chords
of length `1` that share endpoints.  This is the same strictness `notes/chord-lemma.md` §0 flags,
and it is what makes §2.2's boundary convention worth stating carefully.

### 1.3 Short lines: the exact bound

`N <= ceil(Lam) - 1` is exact as an integer bound (`N` is the largest integer strictly below
`Lam`), and it does not degenerate at integer `Lam`: at `Lam = 4` it gives `3`, not `4`.  For the
container `[0,4]^2`:

| line | `Lam` | bound |
|---|---|---|
| any horizontal `y = c`, `0 <= c <= 4`, or vertical `x = c` | `4` | `3` |
| `x + y = c`, `0 < c <= 2^{1/2}` (i.e. `Lam <= 2`) | `c sqrt2` | `1` |
| `x + y = c`, `sqrt2 < c <= 3/sqrt2 = 2.1213` | `c sqrt2` | `2` |
| `x + y = c`, `2.1213 < c <= 2 sqrt2 = 2.8284` | `c sqrt2` | `3` |
| `x + y = c`, `c > 2 sqrt2` | `> 4` | `>= 4` (kept out of the LP: weak) |

and the images under the container's `D4` symmetry.  So the diagonal lines through a corner region
are the strong ones — a line cutting off a corner triangle of leg `2` can carry at most **two**
long chords.  `linecuts.build(..., diag_pitch=h)` adds both diagonal families and
`--line-maxcap 3` drops the ones whose bound is `4` or more.  In the runs below the family is

```
570 lines: 56 with bound 1, 56 with bound 2, 458 with bound 3
```

(the axis-parallel grid at pitch `0.02` plus `y, x = 1, 2, 3` exactly, and the two diagonal
families at pitch `0.05`).

---

## 2. The predicate, and the boundary convention

### 2.1 The exact criterion

`notes/chord-lemma.md` Lemma 1 and Corollary 2 are stated for a horizontal line; a rotation is an
isometry taking unit squares to unit squares and the square of angle `theta` to the square of angle
`theta - phi`, so for a line `L` of direction angle `phi`, writing

```
psi = theta - phi ,        d = dist(centre, L) ,        p = (|cos psi| + |sin psi|)/2
```

the chord `S ∩ L` has length `min( 1/|cos psi| , 1/|sin psi| , (p - d)/|cos psi sin psi| )` when
`d <= p` and is empty otherwise, and

```
|chord| >= 1     <=>     d  <=  D(psi) := p - |cos psi sin psi|
                              = (|cos psi| + |sin psi|)/2 - |cos psi sin psi| .
```

`D` runs from `1/2` at `psi ≡ 0 (mod 90 deg)` down to `(sqrt2 - 1)/2 = 0.207107` at `psi = 45 deg`.
This is `linecuts.Dfun` / `linecuts.member`, evaluated in double precision from the pose's
`(c_x, c_y, theta)` and the line's `(n, c_0, phi)` — no lattice, no rounding, no region labels.

### 2.2 The boundary convention: subsets are safe, supersets are not

`mu` is monotone, so `mu(B') <= mu(B_L) <= ceil(Lam) - 1` for **every** `B' ⊆ B_L`: a subset of the
true "chord `>= 1`" set gives a valid (possibly weaker) row, a **superset does not**.  So the safe
direction is to *exclude* poses on the boundary `d = D(psi)`.

`linecuts.member(F, P, eps)` tests `d <= D(psi) + eps`:

* `eps = 0` — the **exact closed criterion**.  A chord of length exactly `1` is a chord of length
  `>= 1`, so `B_L` really is closed and this is the true set, not a superset.  In exact arithmetic
  it is the strongest valid row.
* `eps < 0` — a strict subset of `B_L`, valid *whatever* the floating-point error in `d` and `D`
  (which is `~1e-16`).  **This is the default (`--line-eps -1e-12`) and it is what every run below
  used.**
* `eps > 0` — a superset, **not valid**.  Never used.

The two useful conventions are not equivalent in principle — §3 shows the sharpest aggregate bound
is `9` under `eps = 0` and `12` under `eps < 0` — but they are **numerically identical on every
measure in this repo**, because `packing_dual.ADM = 1e-9` holds every pose `1e-9` inside the
container: a corner grid pose is at `c_y = 0.500000001`, not `0.5`, so its distance to `y = 1` is
`0.499999999 < D = 0.5` and it is strictly inside `B_{y=1}` under both conventions.  Measured mass
sitting exactly on the knife edge of one of `y, x = 1, 2, 3`:

```
tl_B40K_measure.txt               0.0000        tl_C01010101P_measure.txt    0.0000
tl_B40K_measure_pure_nochord.txt  0.0000        tl_A01010101_measure.txt     0.0000
tl_B40K_measure_pure.txt          0.0000        dual_PC1_support.txt         0.0000
```

and every `mu(B_L)` in §4 is identical to twelve digits under the two conventions.

### 2.3 Sanity: the `4 x 4` grid

`python3 search/linecuts.py --check`, on the sixteen axis-parallel unit squares at the grid centres
`(1/2 + i, 1/2 + j)` with mass `1` each (every one of these lines has `Lam = 4`, cap `3`):

```
  closed  d<=D        y=0.5:4  y=0.9:4  y=1:8  y=1+1e-7:4  y=1.5:4  y=2:8  y=2.5:4  y=3:8
  strict  d<D         y=0.5:4  y=0.9:4  y=1:0  y=1+1e-7:4  y=1.5:4  y=2:0  y=2.5:4  y=3:0
```

Reading it:

* **`y = 1` closed gives `8`**: the four row-`0` squares (top edges on the line, `d = D = 1/2`) and
  the four row-`1` squares (bottom edges).  This is the touching the coverage-`<= 1` relaxation
  cannot see, and it is `8` against a cap of `3`.
* **`y = 1` strict gives `0`**: all eight are exactly on the knife edge.  This is the one place the
  convention bites, and §3 turns it into the difference between an aggregate ceiling of `9` and one
  of `12`.
* **every other line gives `4`**, under either convention, against a cap of `3`.  So the grid is
  cut by *any* generic horizontal line, not just by `y = 1, 2, 3`.

That the grid violates the cut is the consistency check, not a discovery: the `4 x 4` grid is not a
packing of closed disjoint squares (consecutive squares share an edge), and mass `1` on each of the
sixteen is already infeasible for the coverage rows (`T4LEAF.md` §4).

---

## 3. What the whole family can possibly say

This is the part that decides the question, and it is independent of any LP.

**Proposition A (the aggregate identity).**  For a pose `(c_x, c_y, theta)` the set of heights `c`
at which the horizontal line `y = c` cuts a chord of length `>= 1` is the interval
`[c_y - D(theta), c_y + D(theta)]`, of length `2 D(theta)`; and it is contained in `[0, 4]`, because
a unit square inside the container has `p <= c_y <= 4 - p` and `D <= p`.  Hence for **every**
measure `mu` on poses,

```
integral_0^4  mu(B_{y=c})  dc   =   integral  w(theta)  dmu ,
      w(theta) := 2 D(theta) = |cos| + |sin| - 2|cos sin|  ∈  [ sqrt2 - 1 , 1 ] ,
```

with `w = 1` exactly on axis-parallel poses and `w = sqrt2 - 1 = 0.414214` at `45 deg`.

**Corollary B (the uniform ceiling).**  The line cuts give `mu(B_{y=c}) <= 3` for every `c`, so
`integral w dmu <= 12`.  On an axis-parallel measure `w ≡ 1`, so the **whole horizontal family
bounds the mass of an axis-parallel measure by exactly `12`, never below**; in general it bounds
the mass by `12/(sqrt2 - 1) = 12(sqrt2 + 1) = 28.97`.

**Corollary C (the sharpest weighting, and where `9` comes from).**  A weighting `omega_c >= 0` of
the horizontal lines yields `mass <= 3 (sum omega) / m`, where `m` is the least total weight any
admissible pose collects on its catch interval.  Under the **closed** criterion, the three unit
spikes at `c = 1, 2, 3` collect `>= 1` from every axis-parallel pose — such a pose has
`c_y ∈ [1/2, 7/2]`, so the closed interval `[c_y - 1/2, c_y + 1/2]` always contains one of
`1, 2, 3` — giving

```
mu(B_{y=1}) + mu(B_{y=2}) + mu(B_{y=3})  <=  9        =>        mu(axis-parallel)  <=  9 ,
```

and `3` is optimal (`[0,1]`, `[3/2, 5/2]`, `[3,4]` are three disjoint windows each needing weight
`1`).  So **`9` is the best the horizontal family can do about axis-parallel mass, not a first
step.**  Under the **strict** criterion the windows are open, `(0,1), (1,2), (2,3), (3,4)` are four
disjoint windows each needing weight `1`, and the best is `3 x 4 = 12` — Corollary B again.

**What this means.**  Write `mass = mu(axis-parallel) + mu(tilted)`.  The horizontal family alone
gives `mass <= 9 + mu(tilted)`, so it can only certify `< 12` if the measure carries less than `3`
of tilted mass.  It does not:

```
measure                            mass       non-axis-parallel mass    integral w dmu (ceiling 12)
tl_B40K_measure.txt               11.9957            7.2195                     10.9255
tl_B40K_measure_pure_nochord.txt  12.0000            4.6227                     10.7915
tl_C01010101P_measure.txt         11.9200            6.5132                     10.4726
tl_A01010101_measure.txt          11.7452            4.8708                     10.4466
dual_PC1_support.txt              12.1631            6.3043                     10.1889
```

Every one of them has `1.1 – 1.8` of slack in the aggregate.  **The tilted mass is precisely what
pays for the line cuts**, and `T4LEAF.md` §4's anatomy of the `k = 4` optimum already said why it is
there: the grid cannot be used at mass `1` per cell (two closed squares sharing an edge give
coverage `2`), so the LP smears and tilts, and the tilt that buys it coverage room buys it line-cut
room at the same time.

The same computation for the diagonal families is worse for the cuts: a diagonal line's catch band
has half-width `D(theta + 45 deg)`, which is `0.207` for the axis-parallel mass that dominates
these measures, so the diagonal cuts see almost none of it.

### 3.1 The obvious strengthening, also worth zero

Theorem L's proof never uses `1`: `N` chords of length `>= tau`, pairwise disjoint and closed, in a
segment of length `Lam` satisfy `N < Lam/tau`, so

```
mu( { S : |chord of S on L| >= tau } )  <=  ceil(Lam/tau) - 1                      (the TAU CUT)
```

for every `tau > 0`, with `tau = 1` the line cut.  This is a genuinely different family:
`tau = 4/3` on a horizontal line has bound `2`, and its set contains only squares with
`psi ∈ [41.41 deg, 48.59 deg]` — within `3.6 deg` of `45 deg`, since a chord of length `>= 4/3`
needs both `1/|cos psi| >= 4/3` and `1/|sin psi| >= 4/3` — i.e. exactly the tilted mass the
`tau = 1` cut cannot see.  It is complementary in principle.

Screened over `1 064` lines (axis-parallel at pitch `0.01`, both diagonals at pitch `0.05`) and ten
thresholds `tau ∈ {1, 1.05, 1.1, 1.2, 1.25, 4/3, 1.35, 1.4, 1.41, 2}` against all five measures:
**`0` violated rows, every time.**  The tightest anything gets is `mu = 1.000000` against a bound
of `1` on the short diagonal lines that cut off a corner (`x - y = -3.1`, `Lam = 1.273`, bound `1`)
and `mu = 3.000000` against `3` on the axis-parallel lines through the wall bands.  The largest
excess over the whole `5 x 10 x 1064` screen is `+0.000000` (`< 1e-7`).  `linecuts.chord_length` /
`linecuts.cap_tau` and `--tau` implement it; it was not put in the LP because there is nothing for
it to cut.

---

## 4. The screen: every optimum already satisfies every cut

`python3 search/linecuts.py --measure FILE ...` and the finer sweep in §6's reproduce block.  On
`8 002` horizontal and vertical lines at pitch `0.001` (plus `1, 2, 3` exactly), cap `3` each:

```
measure                            mass       max mu(B_L)      on         excess    #lines over cap by >1e-7
tl_B40K_measure.txt               11.995681   3.000000025   x = 0.501   +2.52e-08            0
tl_B40K_measure_pure_nochord.txt  12.000000   3.000000000   y = 0.651   +9.70e-13            0
tl_C01010101P_measure.txt         11.920049   3.000000000   y = 3.067   -1.47e-13            0
tl_A01010101_measure.txt          11.745193   3.000000000   x = 0.488   +2.99e-13            0
dual_PC1_support.txt              12.163061   2.950016332   y = 0.497   -5.00e-02            0
```

identical under both boundary conventions.  And the three-line quantity of Corollary C:

```
measure                            mu(B_{y=1}) + mu(B_{y=2}) + mu(B_{y=3})   (cap 9)      the x images
tl_B40K_measure.txt                    2.3197 + 2.6335 + 1.7014 = 6.6547                    6.9902
tl_B40K_measure_pure_nochord.txt       2.6358 + 2.5943 + 2.8141 = 8.0441                    8.1368
tl_B40K_measure_pure.txt               2.6527 + 2.4667 + 2.8055 = 7.9250                    8.2018
tl_C01010101P_measure.txt              2.3536 + 2.2878 + 2.3540 = 6.9953                    7.1381
tl_A01010101_measure.txt               2.6891 + 2.2105 + 2.7152 = 7.6148                    7.6468
dual_PC1_support.txt                   2.3593 + 1.9401 + 2.3593 = 6.6587                    6.6587
```

`0.96 – 2.35` of slack against `9` on every one.

Two consequences, both exact:

* **Adding the line cuts to any of these models cannot change its value.**  The old optimum is
  feasible for the new model and optimal for the old one, so the two optima agree.  The `0.000000`
  column of §0 is not a measurement of a small effect; it is this argument, confirmed numerically.
* **The certified packing-side value `12.163060670` survives the line cuts.**  The `dual_PC1`
  measure is `M`-certified at
  `1.000000000` (`runs/dual_PC1_support.txt`'s own header) and satisfies every line cut with
  `0.05` to spare, so the pure relaxation with the whole line family added is still `>= 12.163` —
  and in particular **`> 12`.  There is no branching-free proof here.**

---

## 5. The LP runs

Three detached processes, `t = 4.0`, closed semantics, `r = 1`, HiGHS backend, boundary
duplication on, chord rows on, anchor cliques on, the `570`-line family of §1.3 with
`--line-eps -1e-12`; warm-started from `T4LEAF.md`'s checkpoints (copied read-only from
`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a8e11c0e17344cba4/runs/`).
`runs/launchL{1,2,3}.sh`, logs `runs/tl_L*.log`.

| tag | instance | corners | pattern | warm start |
|---|---|---|---|---|
| `L1B40K` | (a) `k = 4` corner leaf, slots free | `1111` | `........` | `T4LEAF.md`'s `B40K` at `11.9957` |
| `L2PURE` | (b) **pure relaxation, no region equalities at all** | `....` | `........` | `dual_PC1` (mass `12.163061`) + `B40K`'s poses |
| `L3A01010101` | (c) the hardest level-2 leaf | `1111` | `01010101` | `T4LEAF.md`'s `A01010101` at `11.7452` |

`5 100 s` of wall clock each, `3 + 3 + 2 = 8` threads in total, run concurrently.

### 5.1 Positive control: the cuts are wired in, and they can be worth 9

Before the null results, a control that the rows are really in the LP and that the pricer sees
their duals.  Take the *pure* relaxation (no region equalities), an **axis-parallel-only** pose set
(`--seed-dth 90`), and a deliberately starved coverage row set (`--row-pitch 0.6`, one row loop, so
`61` rows):

```
python3 search/t4leaf.py 4.0 XPOS --corners .... --patterns "........" --lines --line-pitch 0.05 \
    --seed-pitch 0.25 --seed-dth 90 --row-pitch 0.6 --corner-rows 0 --rowloops 1 --stages 1 \
    --pose-max 0 --threads 1 --variants

  [s0..........final] cols=144 rows=61 LP=9.000000 M=1.000000000 int=1.000000 conv=True
  [s0..........lines] 3 of 162 lines carry dual, sum 3.000000, 105 tight;
                      top  x=0.85 v=1.0000 mu=3.0000,  x=1.95 v=1.0000 mu=3.0000,  x=3 v=1.0000 mu=3.0000
  RESULT ........ = LP 9.000000 / nolines 18.000000
```

**The line cuts are worth `9.000000` there** (`18 -> 9`), the LP lands on `9` exactly, and the LP's
own dual picks out **three vertical lines carrying multiplier `1` each** — which is precisely the
optimal weighting of Corollary C, found by the simplex rather than by hand.  So the machinery
works, the pricer is line-aware, and Corollary C's `9` is attained.

The same control also shows what the cuts cannot do: `9` is the *ceiling* the family imposes on
axis-parallel mass, and
`9 < 12`, so on a pose set with *only* axis-parallel poses the line cuts would close `t = 4`
outright.  Every run below fails to reproduce that only because the LP has tilted poses available.

### 5.2 The trajectory

`python3 search/t4leaf_table.py --table runs/tl_L*.json`.  A stage is one pricing round taken after
the row and clique generation loop has run to its limit and the model has been re-solved (the
settle solve of `T4LEAF.md` §1.3); `conv` = `M <= 1 + 1e-9` **and** `kmax <= 1 + 1e-6` at that
settle solve, i.e. the value is then a genuine lower bound on the instance over its own pose set.
`ldual` is the sum of the line duals; `lines` is `LP_nolines - LP`, the matched worth of the whole
`570`-line family.

```
tag           st pattern            LP         M     kmax        pure     nolines  pure_nolin    lines      int   cols   rows    cq   ldual conv
L1B40K         0 ........    11.992612  1.037041  1.00718   12.000000   11.992612   12.000000  0.00000  3.99261   8670  17749   247  0.0045  no
L1B40K         1 ........    11.995221  1.000000  1.00000      (not valued: --variant-every 2)  3.99522   8024  24037   208  0.7649 YES
L2PURE         0 ........    12.129881  1.003563  1.00307   12.204111   12.129881   12.204111  0.00000  3.62481   9115  24583   459  0.0005  no
L2PURE         1 ........    12.153620  1.037718  1.00207         ---    12.153620         ---  0.00000  3.62974   8046  28871   511  0.0000  no
L3A01010101    0 01010101    11.742058  1.000000  1.00000   11.877321   11.742058   11.877321  0.00000  3.74206   8572  14634   165  5.1015 YES
L3A01010101    1 01010101    11.759752  1.104573  1.01388      (not valued: --variant-every 2)  3.75975   8090  19160   206  6.6343  no
L3A01010101    2 01010101    11.773816  1.000000  1.00972   11.873041   11.773816   11.873041  0.00000  3.77382   8075  21704   131 10.1972  no
```

Each run stopped when its `5 100 s` budget ran out inside the following stage's row loop, after
settling the stage above.  The trajectories, all still rising:

```
L1 (k = 4 control, cliques + lines)   11.992612 -> 11.995221   (+0.0026; st1 fully certified)
L2 (pure relaxation, cliques + lines) 12.129881 -> 12.153620   (+0.0237; towards the certified 12.163061 from below)
L3 (leaf 01010101, cliques + lines)   11.742058 -> 11.759752 -> 11.773816   (+0.0177, +0.0141)
```

`L2PURE` climbing to `12.1544` with `M = 1.0025` against the certified `12.163061` of
`runs/dual_PC1_support.txt` is the consistency check that the pure instance is set up right: the
line cuts and the anchor cliques together have not moved `nu_f(4)` at all.

**Every `lines` entry is `0.000000`.**  In `L1B40K` stage 0 and `L2PURE` stage 0 the `nolines`
value was obtained by a real matched solve on the same rows and columns (the line duals were
nonzero, so the exact-by-duality shortcut of `t4leaf.py` did not apply); in `L3A01010101` stage 0
it was too, with ten lines carrying duals summing to `5.101465`.  **The line rows are in the
optimal dual basis and still change nothing** — they are tight and degenerate, which is exactly the
situation §3 predicts.

Which lines bind, when any do (`L3A01010101` stage 0, the largest line duals of any record):

```
L3A01010101 st0 (10 lines, dual sum 5.101465)   L3A01010101 st2 (11 lines, dual sum 10.197214)
  x=0.76 v=1.0278 mu=3.0000                      x=3.16 v=4.3263 mu=3.0000
  y=3.12 v=0.8823 mu=3.0000                      x=3.12 v=1.6566 mu=3.0000
  x=0.78 v=0.8216 mu=3.0000                      y=3.12 v=1.0225 mu=3.0000
  x=3.22 v=0.7658 mu=3.0000                      x=3.2  v=0.8104 mu=3.0000
  x=3.1  v=0.7307 mu=3.0000                      x=0.84 v=0.5934 mu=3.0000
  y=3.22 v=0.3948 mu=3.0000                      y=3.1  v=0.4324 mu=3.0000
  y=3.1  v=0.3100 mu=3.0000                      x=0.82 v=0.4146 mu=3.0000
  y=0.8  v=0.0756 mu=3.0000                      y=3.26 v=0.2980 mu=3.0000
```

— all of them axis-parallel lines *through the wall bands*, at `mu(B_L) = 3.000000` exactly, i.e.
they are re-expressing the wall-strip content the `chord` rows and the leaf's own region
equalities already carry (`T4LEAF.md` §3.1).  **No line through the interior `[1,3]^2` ever carried
a dual in any record, and no diagonal line ever did.**  That is the whole story: the cut family
speaks about the frame, where the LP has no slack to give, and is silent about the interior, which
is where the leaf's remaining `3.74` of mass lives.

### 5.3 The matched 2 x 2: cliques against lines

Task item (d): the four cells are computed on **one** pose set and **one** row set, at the settle
solve of each stage, by `--variants` (`LP`, `LP_pure`, `LP_nolines`, `LP_pure_nolines`).

| instance | cliques + lines | cliques only | lines only | neither | cliques worth | lines worth |
|---|---|---|---|---|---|---|
| (a) `L1B40K` st0 | `11.992612` | `11.992612` | `12.000000` | `12.000000` | `+0.007388` | `0.000000` |
| (b) `L2PURE` st0 | `12.129881` | `12.129881` | `12.204111` | `12.204111` | `+0.074230` | `0.000000` |
| (b) `L2PURE` st1 | `12.153620` | `12.153620` | (not valued) | (not valued) | — | `0.000000` |
| (c) `L3A01010101` st0 | `11.742058` | `11.742058` | `11.877321` | `11.877321` | `+0.135263` | `0.000000` |
| (c) `L3A01010101` st2 | `11.773816` | `11.773816` | `11.873041` | `11.873041` | `+0.099225` | `0.000000` |

The two levers are **not** complementary; the line lever is simply zero.  The clique lever is worth
`0.007 – 0.135`, and the chord rows (`LP_nochord` in the trajectory) are worth `0.000000`
everywhere too, reproducing `T4LEAF.md` §3.1 on a different pose set.

### 5.4 Cross-check of the implementation

Same model, HiGHS backend against `scipy.optimize.linprog` (`--backend scipy`), with the chord
rows, the line rows and the boundary duplication all on:

```
seed-pitch 0.4  / line pitch 0.1  / diag 0.1     highs LP 9.000000  nolines 9.000000
                                                 scipy LP 9.000000  nolines 9.000000
seed-pitch 0.2  / line pitch 0.05 / diag 0.05    highs LP 11.000000 nolines 11.000000
                                                 scipy LP 11.000000 nolines 11.000000
```

---

## 6. Reproduce

```sh
# the sanity check of 2.3, the screen of 4, and the tau screen of 3.1
python3 search/linecuts.py --check
python3 search/linecuts.py --pitch 0.05 --diag-pitch 0.05 --top 8 --tau 1.0,1.1,1.3333333 \
    --measure runs/tl_B40K_measure.txt --measure runs/tl_C01010101P_measure.txt

# the positive control of 5.1
python3 search/t4leaf.py 4.0 XPOS --corners .... --patterns "........" --lines --line-pitch 0.05 \
    --seed-pitch 0.25 --seed-dth 90 --row-pitch 0.6 --corner-rows 0 --rowloops 1 --stages 1 \
    --pose-max 0 --threads 1 --variants

# the three LP runs, and the trajectory table
sh runs/launchL1.sh      # (a) k = 4 control, cliques + chord + lines
sh runs/launchL2.sh      # (b) pure relaxation, cliques + chord + lines
sh runs/launchL3.sh      # (c) leaf 01010101, cliques + chord + lines
python3 search/t4leaf_table.py --table runs/tl_L*.json
```

Warm starts are copied read-only from `T4LEAF.md`'s worktree,
`/home/evand/math/square-packing/s12/.claude/worktrees/agent-a8e11c0e17344cba4/runs/`
(`tl_B40K_{poses,dual,cliques}.txt`, `tl_A01010101_{poses,dual,cliques}.txt`), and `dual_PC1` from
`/home/evand/math/square-packing/s12/runs/dual_PC1_support.txt`.

`--lines` is off by default, so `t4leaf.py` without it is exactly the program `T4LEAF.md` ran; the
only other change to it is that `clique_cost` now also returns the line part of a pose's reduced
cost, which is `0` when there are no lines.  That one override is what makes the pricer, `pose_rc`
and `rc_check` all see the line duals (`t4screen.py` takes the whole non-coverage part of a
reduced cost from `clique_cost`).

---

## 7. Verdict

**The idea is correct, the cut family is valid and now implemented, and it is worth exactly
`0.000000` at `t = 4` on every instance measured — including the smeared grid at `12.000000`
itself.  Nothing goes below 12.  In particular (b) fails: the pure relaxation with the whole line
family is still `12.13`, and the certified `nu_f(4)` measure of mass `12.163060670` satisfies
every one of `8 002` line cuts with `0.05` to spare, so no branching-free proof comes out of this.**

Answering the deliverable's questions in order.

*How much do line cuts remove from the smeared grid at `t = 4`?*  **`0.000000`.**  `L1B40K` stage 0
values the `k = 4` corner leaf with no cliques at exactly `12.000000` with the line cuts on and
exactly `12.000000` with them off, and the measure attaining it — `T4LEAF.md` §4's `4 x 4` grid
smeared — has `mu(B_L) <= 3.000000000` on every one of `8 002` lines at pitch `0.001`, with a
largest excess of `+9.7e-13`.

*Does anything go below 12?*  Only what already was: the `k = 4` control (a certified `11.995221`)
and the leaf (a certified `11.742058`, rising to `11.773816` with the coverage rows certified and
one clique violated by `9.7e-3`), by exactly the amounts the anchor cliques buy.  Nothing new — and
the pure relaxation is at `12.153620` and climbing towards `12.163061`.

*Which lines bind?*  Only axis-parallel lines through the wall bands (`x, y ∈ [0.76, 0.84]` and
`[3.10, 3.26]`), at `mu(B_L) = 3.000000` exactly, and only degenerately — they carry dual (up to
`4.33` on one line, `10.2` summed) but do not move the objective by `1e-6`.  **No line through the
interior `[1,3]^2` and no diagonal line ever carried a dual, in any record of any run.**

*Why, structurally?*  §3.  The whole horizontal family, aggregated with any weights, can say no
more than `mu(axis-parallel poses) <= 9`, and `9` is attained by the three lines `y = 1, 2, 3` —
so the observation the task started from is not the first line of an argument, it is the last.  A
measure of mass `12` therefore needs only `3` of tilted mass to be immune, and the optima carry
`4.6 – 7.2`.  The `4 x 4` grid *is* cut hard (`8` against `3` on `y = 1`, `4` against `3` on every
generic line), but the LP does not use the grid: it uses the grid smeared and tilted, because the
coverage rows forbid the grid outright (two closed unit squares sharing an edge give coverage `2`),
and the tilt that buys it coverage room buys it line-cut room for free.  The obvious strengthening
— the `tau`-threshold family of §3.1, which catches tilted squares instead of axis-parallel ones —
is also never violated, over `1 064` lines and ten thresholds.

*What is worth keeping.*  The general-line lemma itself (§1) is a strictly stronger statement than
`notes/chord-lemma.md`'s `wall_strip_le_three` at no extra proof cost — it is the same proof with
`t` replaced by `Lam` — and its short-line corollary (a line cutting off a corner triangle of leg
`2` carries at most **two** long chords) is a real geometric fact about corners that no wall strip
expresses.  It costs nothing to state in Lean alongside the existing lemma.  But as an **LP cut
family at `t = 4` it is worth zero**, for the same structural reason `T4LEAF.md` §3.1 found the
wall-strip rows worth zero: at `t = 4` the closed-square coverage constraint already implies
everything the chord counting can say about the frame, and the frame is all the chord counting can
see.

*What would be needed instead.*  Anything that constrains the **interior** `[1,3]^2` and sees
**tilt**.  The line family sees the frame and prefers axis-parallel poses; the anchor cliques
(`search/CLIQUE_CONTINUUM.md`) are the only lever in this repo that does the opposite, and they are
the only one that moves any of these numbers.  A line-cut analogue with a chance would have to be
posed on the *tilted* geometry — e.g. a family of `tau`-cuts on lines whose direction is chosen per
candidate pose rather than from a fixed grid, separated rather than enumerated — but §3.1's screen
says there is nothing there to separate at these optima.
