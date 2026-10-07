# farfield-strong: the near-axis mass cap in the strongest certifiable `t = 4` family (2026-09-22)

*Task: `tasks/farfield-strong/README.md`.  Code: `search/farfield_strong.py` (new).  It **imports**
`search/cliquelever.py` (the family's LP: the exact pose set, the exact arrangement-vertex row
generation, the exact finalisation), `search/leaf_ceiling.py` (arrangement, regions),
`search/rankfamily.py` (the odd-polygon rows) and `search/arch_farfield.py` (the exact near-axis
test) and modifies none of them.  The only new object is one LP row,
`sum_{tilt(S) < eps} mu_S <= K`, plus the guards that keep the exact finalisation inside it.  Runs
in the gitignored `runs/farfield_strong_*`, certificates in `search/farfield_strong_*_exact.txt`,
so every number is quoted here.  Semantics as everywhere in this repo (`notes/s13-casefree.md` §1):
closed unit squares, closed containment, closed container.*

*All cells finished; nothing was left running, and no watcher processes were left behind.*

---

## 0. Verdict, up front

**There is no `eps_ff`, no cell is certified `< 12` cover-side, and no region of the `s(12) = 4`
architecture is proved.**  Two findings, in order of how much they change the picture.

> **1.  The cap is *refuted* in the strongest family at `eps = 1°` and `eps = 2°`.**
> `V_F(4, 1°, 3) >= 11999999870/10^9` and `V_F(4, 2°, 3) >= 11999999855/10^9`, each by an exact
> certificate — a rational measure with maximum coverage `M = 1` exactly over the complete
> arrangement of its support, corner masses `1, 1, 1, 1` exactly, all four chord rows `<= 3`, all
> odd-polygon row of the `1,911` the run carries satisfied, and near-axis mass exactly
> `2.999999977 <= 3`.  The family's
> own **uncapped** certificate on the same pool is `11999999938/10^9`; the two differ by `6.8e-8`,
> which is the `10^-9` round-down of `finalize` and nothing else, and the LP value is
> `12.000000000` in both cases with the cap row's multiplier `lam = 0` — *the row is not even
> binding*.  So `search/ARCH_FARFIELD.md`'s verdict ("A3's candidate fails, and it fails at the
> small `eps` the composition actually needs") now holds in the family that sits at exactly `12`,
> not only in the points-only family that sits at `12.2688`.
>
> **2.  At `eps >= 3°` the packing side cannot decide the question at all, and this is structural.**
> The uncapped value of this family *is* `12`.  A packing-side certificate is a lower bound, and the
> exact finalisation rounds down, so the best it can ever print is `12 - 4.4e-8`: **the rule "a cell
> `>= n` refutes the cap" is unattainable here by construction**, and only a cover-side certificate
> could show a cell `< 12`.  §6: there is no instrument for one — in the points-only family a cover
> below `12` is impossible at small `eps` (the packing side is already `>= 12.0855` there), and in
> this family the LP duals are not covers at all (`search/HONEST.md`: honest cost `18.1`, floor
> `15.8`).  The pattern argument that bounds the pinned Bentz leaf by `12` gives nothing under a cap,
> because no pattern of `P0` constrains the angle (`|A - B| = 0.1216`, against the diameter `1.414`).

What the runs do deliver is the first measurement of **the cost of the cap in the strongest family**,
on a pool that reproduces that family's optimum exactly, uncapped:

| `eps` | `1°` | `2°` | `3°` | `4°` | `5°` | `10°` | `15°` | `20°` | `30°` |
|---|---|---|---|---|---|---|---|---|---|
| cost of `mu(tilt < eps) <= 3` | **`0`** | **`0`** | `0.2357` | `0.2408` | `0.1096` | `0.2290` | `0.3489` | `0.4318` | `0.5811` |

(the `5°` and `10°` entries are `0.2520` and `0.3155` on the same pools as their neighbours and
`0.1096` and `0.2290` once the pool is given one extra family of columns — §4.1; that
non-monotonicity in `eps` *is* the pool caveat, in two numbers)

and, one container down where `Claim(3)` is a theorem, the same shape with `3` for `4`: the cap
`<= 2` is **refuted at `5°`** (`6.008472573 >= 6`) and costs `0.103 / 0.341 / 0.438 / 0.504` at
`10° / 15° / 20° / 30°` — `ARCH_FARFIELD.md` §3.4 reproduced on a pool whose uncapped value
(`6.042662562`) is above `n`.

**Reading the `>= 3°` numbers.**  They are lower bounds on a pool, not values.  The move that makes
the `1°` and `2°` cells free is the one `ARCH_FARFIELD.md` §4 identified — *tilt the corner square to
just past `eps`* — and at those two `eps` the pose it needs (`(0.51036, 0.51036)` at `1.1997°`,
`(0.51715, 0.51715)` at `2.0002°`) is already a **native pose of the family optimum's own 162-pose
support**.  That support has `0.137` of mass in `[1°, 2°)` and only `0.010` in the whole of
`[2°, 5°)`, so from `3°` on the LP has to be *given* the column, and the two substitutes tried here
(rigid rotations of the whole container, and "band-edge" columns that re-angle each pose in place)
are cruder than a real column generation.  The trend is the finding; the individual values are not.

---

## 1. The question, and which side of it every number is on

For a pose `S` write `tilt(S)` for the distance of its angle to `0 mod 90°`.  Let `F` be the
strongest certifiable family the repo has at `t = 4` (`notes/status.md`, `search/BENTZ.md` §7.2,
`search/LEAF_CEILING.md` §5.2): the LP

    V_F(T, eps, K) = sup mu(Pi)
      s.t.  mu({S : p in S}) <= 1            for every point p of [0,T]^2      (coverage)
            mu(X_0 u ... u X_{k-1}) <= (k-1)/2   for every odd-polygon row     (polygons)
            mu({S : centre(S) in C_i}) = 1   i = 0..3                          (corner k = 4)
            mu({S : c_y(S) <= 1}) <= 3       and its three images              (chord)
            mu({S : tilt(S) < eps}) <= K                                       (the cap)
            mu >= 0

whose **uncapped** value is exactly `12.000000000` as an LP and `11.999999926` as the exactly
certified measure published in `notes/status.md` (`search/pgonly_corner_exact.txt`, the `PGCORN`
configuration of `runs/launch_2026-09-12_pgonly.sh`); §3 improves that certificate to
`11.999999956`.

**Direction, which is the whole discipline of this note** (`search/ARCH_FARFIELD.md` §1,
`notes/review-2026-09-13b.md` "Clarifications"):

* **Packing side.**  A measure supported on a *finite* pose set that satisfies every row above is
  feasible for the continuum problem, so its mass is a **lower bound** on `V_F`.  Restricting the
  pose set can only lower it, and the near-axis set is over-estimated (below), so the cap is if
  anything too strong.  **A cell `>= n = 12` refutes the cap at that `eps`; a cell `< 12` proves
  nothing at all** unless the same pool's *uncapped* value is itself `>= 12`, and even then it is
  only evidence about that pool.
* **Cover side.**  An upper bound needs a dual — point weights `y`, polygon multipliers `zeta`,
  four region multipliers `pi_i` (free sign, the rows are equalities), four chord multipliers and
  the cap multiplier `lam >= 0` — whose *capture* is `>= 1` at **every** admissible pose of tilt
  `>= eps` of the continuum and `>= 1 - lam` at every near-axis one.  The bound is then
  `|y| + sum_i pi_i + 3 sum chord + sum (k-1)/2 zeta + K lam`.  **A cell certified `< 12` on that
  side would be a theorem.**  §6 reports what happened to this side.

**A structural point the brief's step 1 cannot be given literally, and it changes the reading of
every cell.**  The brief asks for "a pool whose uncapped value is `>= 12`, otherwise the cell
measures the pool".  In *this* family that is unattainable: the family's own uncapped value is
`12.000000000` (the LP value, reproduced here and in `PGCORN`, `BZA`, `BZAC`, `T4SCREEN` §3 and the
`UNION2` polish of `LEAF_CEILING.md` §5.2 — it has never gone above `12` on any pose set), and the
exact finalisation rounds every mass **down** to a multiple of `10^-9`, so the best certified number
the family can produce is `12 - O(10^-8)` (`11.999999926` published, `11.999999956` here).
Consequently:

> **In the strongest family the packing side can never reach `12`, so it can never refute the cap by
> the rule "a cell `>= n`".  The refutation test here is the *cost* of the cap: a cell whose exact
> value equals the uncapped cell's, to the last digit of the round-down, is a cell where the cap
> costs the family nothing.  And a cell strictly below `12` is, as always, silent — only a
> cover-side certificate could turn it into a theorem.**

That is not a defect of the pool: the pool below reproduces the family optimum exactly, uncapped,
which is the strongest calibration `ARCH_FARFIELD.md` §5's rule can ask for (its `P4c` was `0.36`
short of the truth *before* any cap).

**The near-axis test is exact and errs towards a stronger cap.**
`sin(tilt) = min(|cos|,|sin|) = min(|q^2-p^2|, |2pq|)/(q^2+p^2)` for `theta = 2 arctan(p/q)`, and
the test is `min(|q^2-p^2|,|2pq|) * 10^12 <= ceil(sin(eps) 10^12) (q^2+p^2)` — integers only, with
`sin(eps)` rounded **up**, so the set called near-axis is a **superset** of the true one
(`arch_farfield.near_axis_exact`, unchanged).

---

## 2. The instrument

`search/farfield_strong.py` has four modes.

**`pool OUT --exact F... --float8 F... --rot A,B,...`** builds a `sym 1` exact pose pool: exact
measure files (the 8 dihedral images expanded when the file says `D4-symmetrised`), float support
files, and **rigid rotations** of either.  A rotation of the whole container about its centre
preserves coverage exactly and shifts every angle by `alpha`, so it converts a known optimum into
columns whose originally-axis-parallel mass sits at tilt exactly `alpha` — the far-field columns of
`arch_farfield.py rotseed`, which is what a capped LP needs and what
`search/ARCH_FARFIELD.md` §5(b) says the `eps >= 5°` pools were short of.

**`cap TAG --pool P --eps E... --cap K... [--uncapped]`** is the certificate.  It builds
`cliquelever.Poses` on the pool, resumes the `1,919` odd-polygon rows of
`runs/inputs-2026-09-12/cl_E2Bg_pgons.json` (membership **re-derived from the exact rational
anchors** over this pose set, so no stale member list is inherited) and the `47,296` exact coverage
rows of `cl_E2Bg_rows.txt`, pins the four corner boxes at `1`, adds the four chord rows and the cap
row, and then row-generates on the **exact arrangement of its own support** until no arrangement
vertex is violated.  At that point the LP is the exact optimum of the capped family **restricted to
that pose set** (`LEAF_CEILING.md` §5.4 Lemma 1: the arrangement of a support is a complete row set
for every sub-measure of it).  The masses are then rounded **down** to multiples of `10^-9`; if the
exact maximum coverage still exceeds `1`, or the cap is exceeded, the measure is scaled down
exactly; the four corner equalities are topped back up **only on exact slack of every row of the
family** (coverage, polygons, chord) *and of the cap*.  Every load-bearing comparison is a Python
integer or `Fraction`.

All cells of one run share the pose set and the row set (rows are only ever added and every row is
valid for every cell), so the exact arrangement is enumerated once and each `(eps, K)` is a fresh
LP plus its own row generation.

**`check FILE --eps E --cap K --pgons J`** re-derives mass, `M`, near-axis mass, corner masses,
chord masses and the polygon rows from the measure file alone, on a code path that shares nothing
with the LP.

**`table TAG...`** collects the JSON records.

**Float** (chooses *which* measure gets certified, never quoted as a bound): the HiGHS solve, the
pose snapping in `pool`, the float prefilter of the incidence tests (decided exactly inside a
`10^-7` band).  **Exact**: the poses, the arrangement vertices, every incidence, the mass, `M`, the
near-axis mass, the region and chord masses, and every polygon row.

---

## 3. Calibration: the pool, and its uncapped value

Two pools were built (`runs/farfield_strong_pool4*.txt`).  Both start from the family's own
certified optimum and add the columns a capped LP needs — tilted ones.

| pool | poses | built from |
|---|---|---|
| `pool4` | `9,642` | `pgonly_corner_exact` (162) + `pgonly_pure_exact` (1,423) + `cover4_exact_support` D4-expanded (2,344) + the five `P4c` and two `P4a` capped exact supports of `ARCH_FARFIELD.md` (1,544) + `arch_farfield_pool4` D4-expanded (1,184) + `pgonly_corner_exact` rigidly rotated by `±1, ±2, ±3, ±5, ±8, ±12, ±16, ±20, ±25, ±30, ±35, ±40°` (2,985) |
| `pool4b` | `5,039` | the same without `pgonly_pure_exact`, `cover4_exact_support` and `arch_farfield_pool4`, with rotations at `±2, ±4, ±6, ±8, ±11, ±14, ±17, ±21, ±25, ±30, ±35, ±40°` |

`pool4b` is the one the cells below run on: on `pool4` (`9,642` columns against `60,000`+ exact
rows) a single LP solve costs `4–6` minutes, which is not affordable across eleven cells, and the
uncapped value is the same.  **The rotated columns are the point of the pool**: a rigid rotation of
the container by `alpha` preserves coverage exactly and moves every angle by `alpha`, so it turns
the family optimum — whose mass is overwhelmingly near-axis — into an equally good configuration
all of whose axis-parallel mass now sits at tilt exactly `alpha`.  That is precisely the move
`ARCH_FARFIELD.md` §4 found the LP making by hand (the corner square at `1.2°`), supplied in
advance at every `alpha` the `eps` grid needs, and it is what §5(b) of that note says the
`eps >= 5°` pools were missing.

**How far the cap is from the uncapped optimum.**  The family optimum is even more near-axis than
the points-only one: `search/pgonly_corner_exact.txt` puts `7.597978` of its `12` at tilt **exactly**
`0` (41 poses) and `9.150436` below `0.5°`, so at every `eps` of the grid the cap `3` has to remove
two thirds of the measure and the LP has to rebuild it from different columns
(`arch_farfield.py hist`):

| `eps` | `1°` | `2°` | `5°` | `10°` | `15°` | `20°` | `30°` |
|---|---|---|---|---|---|---|---|
| near-axis mass of the uncapped optimum | `9.196589` | `9.333330` | `9.343817` | `9.697728` | `10.237279` | `10.327624` | `10.736481` |
| restrict-and-refill floor `far + min(near, 3)` | `5.803411` | `5.666670` | `5.656183` | `5.302272` | `4.762721` | `4.672376` | `4.263519` |

So nothing is inherited: every capped cell below is a different fractional packing, not a
restriction of the known one.

**The uncapped value on the pool, certified** (`runs/farfield_strong_S4A.log`, cell `_uncapped`):

    LP = 12.000000000          exactly, at every iteration
    EXACT mass = 2999999989/250000000 = 11.999999956
    M = max coverage = 1 exactly;  corner masses 1, 1, 1, 1 exactly;
    chord 2.999999996 … 2.999999998 <= 3;  1,696 polygon rows all satisfied
    -> search/farfield_strong_S4A_uncapped_exact.txt

`11.999999956` is `4.4e-8` below `12` — the round-down of `finalize`, and nothing else.  It is the
best certified measure this family has produced (`search/pgonly_corner_exact.txt` is
`11.999999926`), and it settles the calibration question in the only form it can take here: **the
pool is not what is short.**

---

## 4. `T = 4`, `n = 12` — the table

**How to read it.**  `packing side` is the exactly certified mass of a measure that satisfies every
row of the family *and* the cap — a rigorous **lower** bound on `V_F(4, eps, K)`.  `cover side` is a
certified **upper** bound; there is none (§6), and the column says so rather than being left blank.
`status` is one of

* **cap free** — the certified value equals the uncapped one to the last digit of the `10^-9`
  round-down, and the LP value is `12.000000000`: the cap costs the family nothing, and A3 fails at
  that `eps` in this family;
* **cap costs `c`** — the LP value is `12 - c` with `c > 0` on a pool that reproduces the family
  optimum uncapped.  **This is not a proof that `V_F < 12`** (a restricted pose set can only lower
  the LP), and it is not a proof that `V_F >= 12` either.  It is the size of the bite;
* **pool-limited** — the cell did not converge inside its clock, so its number is a scaled-down
  iterate and is quoted only as a floor.

### 4.1 The table (`cap = 3 = T - 1`)

Every entry is an exactly certified measure: `M = 1` exactly over the complete arrangement of its
own support, corner masses `1, 1, 1, 1` exactly, all four chord rows `<= 3`, every polygon row
satisfied, near-axis mass `<= 3` exactly — the `certified=True` of the `RESULT` lines.

| `eps` | coarse pool `4b` | fine pool (`4c`/`4d`/`4e`: rotations at `eps + 0.05°`) | band-edge pool `4f` | **packing side, best** | cover side | status |
|---|---|---|---|---|---|---|
| — (uncapped) | `11.999999956` | `11.999999938` / `…936` / `…948` | `11.999999941` | `11.999999956` | none | the family value: `12` up to the `10^-9` round-down |
| `1°` | `11.779080538` | `11.999999870` | — | **`11.999999870`** | none | **cap free** — A3 fails at `1°` |
| `2°` | `11.768448016` | `11.999999855` | — | **`11.999999855`** | none | **cap free** — A3 fails at `2°` |
| `3°` | — | `11.764338358` | — | `11.764338358` | none | cap costs `0.2357` on this pool |
| `4°` | — | `11.759237923` | — | `11.759237923` | none | cap costs `0.2408` on this pool |
| `5°` | `11.746694150` | `11.747961588` | **`11.890363723`** | `11.890363723` | none | cap costs `0.1096` on the best pool |
| `10°` | `11.684456589` | `11.682792169` | **`11.770986233`** | `11.770986233` | none | cap costs `0.2290` on the best pool |
| `15°` | `11.651051049` | `11.648871483` | — | `11.651051049` | none | cap costs `0.3489` on this pool |
| `20°` | `11.568239631` | `11.558462035` | — | `11.568239631` | none | cap costs `0.4318` on this pool |
| `30°` | `11.413752361` | `11.418864512` | — | `11.418864512` | none | cap costs `0.5811` on this pool |

Cell files: `search/farfield_strong_S4{A,B,C,D,E,F,G,H}<cell>_exact.txt`, runs
`runs/farfield_strong_S4*.{log,json}`.

**The `5°` column is the warning label for the whole right-hand side of this table.**  `pool4f`
differs from `pool4c` by one family of columns — every pose of the optimum's support kept where it
is and **re-angled** to tilt `5.05°, 5.3°, 6°, 7°, 10.05°, 10.5°, 12°` (the `--tilt` option; the
whole-container rotation does the same thing but drops poses out of the box).  That one family takes
`V_F(4, 5°, 3)` from `11.7480` to `11.8904` and `V_F(4, 10°, 3)` from `11.6845` to `11.7710` — it
recovers `56 %` and `27 %` of the apparent cost.  Nothing says the next family of columns would not
recover the rest.

### 4.2 The mechanism, and why the `1°` and `2°` rows are the important ones

`ARCH_FARFIELD.md` §4 diagnosed the evasion in the points-only family: *told it may keep only `3`
units of near-axis mass, the LP rotates the corner square by `1.2°` and leaves everything else where
it was.*  The two pools above turn that diagnosis into a controlled experiment, because the only
difference between them is **whether the pose set contains poses at tilt just above `eps`**:

* `pool4b` rotates the family optimum by `±2, ±4, ±6, ±8, ±11, ±14, ±17, ±21, ±25, ±30, ±35, ±40°`.
  At `eps = 1°` its nearest "just outside the band" column is at `2°`, and the cap costs `0.221`.
* `pool4c` adds `±1.05, ±1.2, ±1.5, ±2.05, ±2.3, ±3, ±5.05, ±5.5, ±7, ±10.05, ±11, ±14°`.  At
  `eps = 1°` the LP now has a column at `1.05°`, and the cap costs **nothing**: `LP = 12.000000000`
  with the cap row's multiplier `lam = 0` — the row is not even binding — and the certified measure
  is `11.999999870` with near-axis mass exactly `2.999999977 <= 3`.  Same at `2°`.

So the `0.221` of `pool4b` was an artefact of the angle lattice, exactly as `ARCH_FARFIELD.md` §5(b)
warned, and **in the strongest certifiable family the cap `mu(tilt < eps) <= 3` is worth nothing at
`eps = 1°` and `eps = 2°`.**  That is the one conclusion of this note that has the strength of the
`ARCH_FARFIELD.md` verdict, and it extends that verdict from the points-only family (`12.0986`, cap
cost `0.0025`) to the family that actually sits at `12`.

**Where the LP puts the mass** (`heaviest` in `runs/farfield_strong_S4D.json`; the measure is not
symmetrised, so each dihedral image is its own pose).  The uncapped optimum is the smeared tiling:
the four corner squares carry mass exactly `1` each at `(0.5, 0.5)`-type centres and tilt `0`, and
`8.587` of the `12` is near-axis.  At `eps = 1°, cap = 3` the LP does one thing:

```
  uncapped optimum                     eps = 1 deg, cap = 3            eps = 10 deg, cap = 3
  0.50000 0.50000   0.0000 deg 1.000   0.51036 0.51036 +1.1997  0.443   3.50000 3.50000 0.0 0.345 near
  0.50000 3.50000   0.0000 deg 1.000   3.48964 0.51036 +1.1997  0.422   3.50000 0.50000 0.0 0.340 near
  3.50000 0.50000   0.0000 deg 1.000   0.51036 3.48964 -1.1997  0.422   0.50000 3.50000 0.0 0.325 near
  3.50000 3.50000   0.0000 deg 1.000   3.48964 3.48964 -1.1997  0.420   0.50000 0.50000 0.0 0.312 near
  0.50000 2.50000   0.0000 deg 0.415   3.48964 3.48964 +1.1997  0.404   2.50000 3.50000 0.0 0.189 near
```

`(0.51036, 0.51036)` at `1.1997°` is **the same pose** `ARCH_FARFIELD.md` §4 reported as the single
heaviest atom of the points-only capped optimum ("the corner square, tilted", carrying `1.81075`).
Here it is supplied to the LP in advance as a rotation column, and the LP takes it, splits it over
the two chiralities, and loses nothing.  At `eps = 10°` that move is no longer available at a price
the corner equalities can pay, and the optimum goes back to four axis-parallel corner squares
sharing the cap's `3` between them.

From `eps = 3°` on the same trick stops paying for itself *on these pools*: the fine rotations are
present (`3.05°`, `4.05°`, `5.05°`, `10.05°`, `15.05°`, `20.05°`, `30.05°`) and the cells still read
below `12`, because rotating the whole configuration by more than a couple of degrees drops too many
poses out of the container for the LP to refill.  Adding the per-pose band-edge columns recovers
more than half of the apparent cost at `5°` (§4.1).  **But a cell below `12` is silent** (§1), and
nothing here says whether `V_F(4, 5°, 3)` or `V_F(4, 10°, 3)` is above or below `12`.

There is one concrete reason to suspect the `>= 3°` cells are still about the pool.  Look at which
pose does the evading: at `eps = 1°` it is the corner square at `1.1997°`, at `eps = 2°` the corner
square at `2.0002°` — and **both are native poses of the family optimum's own 162-pose support**
(its tilt histogram has `0.137` of mass in `[1°, 2°)` and only `0.010` in all of `[2°, 5°)`).  The
support simply has nothing at `3°–5°` to tilt onto, and the substitutes this note supplies —
whole-container rotations, and "band-edge" columns that re-angle each pose in place — are both
cruder than what a real column generation would find.  So the honest reading of the `>= 3°` rows is
**"the cap costs *this pool* that much"**, and the number to watch is not the value but the trend.

### 4.3 The cap ladder at `eps = 10°`, and the `cap >= 4` filter

| `cap` | `1` | `2` | `3` | `4` | `9` |
|---|---|---|---|---|---|
| packing side at `eps = 10°` | `11.306579704` | `11.558727880` | `11.684456589` | `11.773057383` | **`11.999999931`** |

Monotone in the cap, as it must be.  `cap = 9 = (T-1)^2` — the cap the composition actually needs —
is **not binding at all**: the LP sits at `12.000000000`, the certified measure is `11.999999931`
with near-axis mass `8.999999975`, and the cap is free.  Exactly as `ARCH_FARFIELD.md` §6 caveat 3
found in the points-only family, now in the strongest one.

**The brief's filter — "a cell with `cap >= 4` at
`eps <= 12.8°` reading `< 12` is a bug, because a genuine margin-`0` packing with `>= 4` near-axis
squares exists there" — does not apply**, and this note does not treat those cells as bugs:
`ARCH_FARFIELD.md` §5(a) and `notes/proof-architecture.md` §0a item 10's own erratum settle it.  A
margin-`0` configuration has touching squares; in closed semantics two touching closed unit squares
give coverage `2` at the contact, so the configuration is *infeasible* for this relaxation, and with
a wall-to-wall chain of four there is no slack to smear the contacts away — the chain of four can
carry at most `2`, so the family is worth at most `10` to the LP, not `12`.

---

## 5. `T = 3`, `n = 6` — the control

`Claim(3)` is a theorem (Kearney–Shiu), so if the cap mechanism is any good at all it must close
one container down; `ARCH_FARFIELD.md` §3.4 found that in the points-only family it only starts to
bite past `eps = 10°`, and even at `30°` takes `V` only `0.49` below `6`.  The question here is
whether the polygon rows change that.

Family: points + odd polygons, **no** region equalities and **no** chord rows — §7 item 3 explains
why the literal `T = 3` analogue of the `t = 4` family is degenerate (it reads exactly `5` uncapped).
Pool `runs/farfield_strong_pool3.txt`, `8,105` poses: all thirty-two `P3a` and twenty `P3b` capped
exact supports of `ARCH_FARFIELD.md` §3.3 with their eight dihedral images, plus the uncapped and
the `30°`-capped optima rotated by `±1 … ±40°`.  The polygon rows are separated here (`--pent
5,7,9`), since no `t = 3` polygon file existed.

**The polygon rows had to be dropped, and the reason is the calibration rule itself.**  With the
separated odd-polygon rows in, the pool's own **uncapped** value is `5.965582893` — *below* `n = 6`
with no cap at all (`runs/farfield_strong_S3B.log`), which by `ARCH_FARFIELD.md` §5 makes every
capped cell under it meaningless.  Without them it is `6.042662562 > 6`
(`runs/farfield_strong_S3C.log`), so the control is run points-only and the polygon variant is
recorded here as a pool failure, not as a result.

| `eps` | packing side, `cap = 2` | cover side | status |
|---|---|---|---|
| — (uncapped) | `6.042662562` | none | calibration passes: the pool is above `n = 6` |
| `5°` | **`6.008472573`** | none | **`>= 6`: the cap is refuted at `5°`** |
| `10°` | `5.939504605` | none | cap costs `0.103` |
| `15°` | `5.701550603` | none | cap costs `0.341` |
| `20°` | `5.604719146` | none | cap costs `0.438` |
| `30°` | `5.538869201` | none | cap costs `0.504` |

This reproduces `ARCH_FARFIELD.md` §3.3–§3.4 on a better pool and confirms its reading: **at `T = 3`,
where `Claim(3)` is a theorem (Kearney–Shiu), the cap `mu(tilt < eps) <= T - 1` does not close the
relaxation until somewhere past `10°`, and even at `20°` it takes `V` only `0.44` below `n`.**  The
mechanism therefore fails to produce the conclusion in a case where the conclusion holds, at every
`eps` a tube of radius `~1/T` could reach.

---

## 6. The cover side — why no certificate exists here, and what would be needed

**No cover-side number was obtained, and the reason is structural, not a shortage of compute.**

1. **In the points-only family a cover-side bound below `12` is impossible at any `eps`, `K >= 0`.**
   The capped dual is `min |y| + K lam` subject to `capture_y(S) >= 1` for every pose of tilt
   `>= eps` and `capture_y(S) >= 1 - lam` for every near-axis pose.  Dropping the cap row entirely
   only lowers the objective, and what is left is a cover of the *tilted* poses; but the capped
   value is bounded below by the mass of any feasible capped measure, and `ARCH_FARFIELD.md` §3.2
   certifies `V_points(4, eps, 3) >= 12.0986` at `eps <= 1°` and `>= 12.0855` at `2°`.  So at small
   `eps` the points-only cover side cannot go below `12` — it is already refuted from below.
2. **In the strong family the dual is not a cover.**  `search/HONEST.md` §0 measured exactly this
   object: the honest cost of the best `t = 4` points + polygons (+ clique) dual, read as a
   weighted cover of `[0,4]^2` over the **continuum** of poses, is `18.1` and `20.2`, with `15.8`
   the smallest number anywhere in that measurement — against the `12` it would have to beat.  The
   LP duals of the runs below are feasible only for the poses that are loaded; the `0.04 / 2.5°`
   lattice pricer under-reports the worst pose by a factor of six (`HONEST.md` §0 item 3), so
   "no improving column" is not a convergence claim over the continuum, and nothing in this
   pipeline produces a dual that is feasible off its own pose set.
3. **The one structural upper bound the repo has does not survive the cap.**  `BENTZ.md` §0 bounds
   the fully pinned leaf `(AB)^4` by `12` with a two-line argument: the twelve admitted patterns
   partition the admissible poses, and `mu(R_pi) <= 1` by a coverage row at any point of `pi`.  To
   push that below `12` with the cap one would need `>= 4` of the twelve patterns to be realisable
   **only** by near-axis poses (then `sum_pi mu(R_pi) <= 8 + K = 11`).  The patterns are the four
   corner pairs `{a_i, b_i}` and the eight singletons `{c_j}`, `{d_j}`; the pair `A(1, 457/500)`,
   `B(457/500, 1)` is at distance `0.086 sqrt 2 = 0.1216`, far from the unit square's diameter
   `sqrt 2`, so a square containing both is under no angle restriction at all, and a singleton
   pattern is under none by construction.  **The pattern argument gives nothing at any `eps`.**

So the cover side of this question is open in the same sense it was open before this task, and
`ARCH_FARFIELD.md` §6 caveat 4 stands unchanged.  What would settle it is a dual certified over the
continuum — the `verify/` exact verifier does that for *point-weight* covers, but it has no
primitive for polygon rows, region equalities, chord rows or a near-axis multiplier, and the
points-only route is closed by item 1.

---

## 7. Caveats, and what turned out to be wrong in the brief's premises

1. **"Build a pool whose uncapped value is `>= 12`" cannot be done in this family, and the rule it
   comes from does not apply here.**  §1: the family's own uncapped value *is* `12.000000000`, and
   the exact finalisation rounds masses **down** to `10^-9`, so `11.999999928` is the ceiling of
   what a certificate can print.  `ARCH_FARFIELD.md` §5's calibration rule ("a cell below `n` whose
   pool is below `n` uncapped measures the pool") was written for the *points-only* family, where
   the truth is `12.2688` and a pool can sit well above `12`.  Here the right calibration is
   "does the pool reproduce the family optimum uncapped?", and it does, exactly.  **Consequently no
   cell of this note can be a refutation by the rule "packing side `>= n`", and none is claimed.**
2. **The chord rows are a flagged hypothesis, not a theorem.**  *(Superseded 2026-09-22: they are a theorem, `lean/Sqpack/Chord.lean` `wall_strip_le_three`; see the correction in `search/T4LEAF.md` §1.2. The paragraph is kept as written.)*  `search/T4LEAF.md` §1.2:
   `mu({c_y <= 1}) <= 3` is Nagamochi 2005 Lemma 7(i) / Stromquist's chord lemma for squares of side
   `L > 1`, "**still a hypothesis in this repo** — it is a strictness statement, false for closed
   unit squares, and is being proved separately".  `notes/status.md`'s row "corner `k = 4`, points +
   polygons + regions + chord — `11.999999926` certified" therefore rests on it, and so does every
   cell here that carries `--chord`.  This is harmless for the packing side — an extra row only
   *restricts* the measure, so a certified value stays a valid lower bound on the chord-free
   value — but it would invalidate any cover-side bound, and it is worth saying plainly that the
   "strongest **certifiable** family" of the brief is not certifiable in one of its four
   ingredients.  §4 reports the `--chord`-free variant as well.
3. **The `T = 3` analogue of the family is degenerate, and the control had to be weakened.**  With
   `corners = 1111` **and** the chord rows at their `T = 3` value `T - 1 = 2`, the LP sits at
   exactly `5.000000000` from the first solve, for a structural reason: each wall strip contains
   two corner boxes, so `mu(strip) <= 2` with the corner masses pinned at `1` forces all eight wall
   slots to `0`, and the interior `[1,2]^2` carries at most `1`.  `4 + 0 + 1 = 5 < 6` **with no cap
   at all**, so that family closes trivially and measures nothing about the cap.  The `T = 3`
   control in §5 is therefore run in points + odd polygons, with no region equalities and no chord
   rows — the strongest family that is not degenerate one container down.
4. **The addendum's "cap `>= 4` at `eps <= 12.8°` must be `>= 12`" filter is not a theorem**, and
   the note does not treat a cell below `12` there as a bug: `ARCH_FARFIELD.md` §5(a) and
   `notes/proof-architecture.md` §0a item 10's erratum — a margin-`0` configuration has touching
   squares, whose closed-square coverage is `2` at the contact, so it is *infeasible* for this
   relaxation and its `t`-average is worth at most `10`, not `12`.

---

## 7b. The certificates, and the two independent `check` commands

The two cells that carry §0's finding, re-derived from the measure file alone with no LP and
integers only, on two code paths that share nothing with each other or with the run:

```sh
# (i) this task's checker: mass, M, near-axis mass, corner masses, chord, polygon rows
python3 search/farfield_strong.py check search/farfield_strong_S4D_e1p0_c3p0_exact.txt \
        --t 4 --eps 1 --cap 3 --corners 1111 --chord-rhs 3 \
        --pgons runs/inputs-2026-09-12/cl_E2Bg_pgons.json --procs 2
#   266 poses, 44272 arrangement vertices (whole container)
#   mass      = 1199999987/100000000 = 11.999999870
#   M         = 1                    OK <= 1
#   near-axis = 2999999977/1000000000 = 2.999999977   OK <= cap 3
#   corner masses ['1','1','1','1'];  chord 2.99999999… <= 3 OK
#   1749 polygon rows re-derived from the exact anchors, OK, worst violation +0.000000000

python3 search/farfield_strong.py check search/farfield_strong_S4D_e2p0_c3p0_exact.txt \
        --t 4 --eps 2 --cap 3 --corners 1111 --chord-rhs 3 \
        --pgons runs/inputs-2026-09-12/cl_E2Bg_pgons.json --procs 2
#   314 poses, 70654 vertices; mass = 2399999971/200000000 = 11.999999855; M = 1;
#   near-axis = 1499999991/500000000 = 2.999999982 <= 3;  corners 1,1,1,1;  1756 polygon rows OK

# (ii) the repo's own checker (search/leaf_ceiling.py, unchanged), coverage + regions + chord
python3 search/leaf_ceiling.py check search/farfield_strong_S4D_e1p0_c3p0_exact.txt \
        --t 4 --corners 1111 --chord --anchor none --procs 2
#   MASS = 11.999999870;  M = max cov = 1 over 44272 vertices  OK;
#   C0..C3 = 1.000000 each, I = 4.000000;  chord S0..S3 = 3.000000 each <= 3 OK
```

Every cell of §4 and §5 has its own file of the same form,
`search/farfield_strong_<TAG><cell>_exact.txt`, with `<TAG>` and `<cell>` as in the tables
(`S4A_uncapped`, `S4B_e30p0_c3p0`, `S4F_e10p0_c9p0`, `S3C_e5p0_c2p0`, …), and each is checked by
the same two commands with its own `--eps`, `--cap` (and `--t 3 --corners '' ` for the `T = 3`
files, which carry no region or chord rows).

---

## 8. Reproduce

```sh
# ---- the pools (all four are `sym 1` exact pose files; `pool4b` is the coarse one) -------------
COMMON="--exact search/pgonly_corner_exact.txt --rot-src search/pgonly_corner_exact.txt"
python3 search/farfield_strong.py pool runs/farfield_strong_pool4b.txt --t 4 $COMMON \
  --exact runs/arch_farfield_P4c_e5p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e10p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e15p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e20p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e30p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4a_e30p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4a_e5p0_c3.0_exact_support.txt \
  --rot 2,-2,4,-4,6,-6,8,-8,11,-11,14,-14,17,-17,21,-21,25,-25,30,-30,35,-35,40,-40      # 5,039

python3 search/farfield_strong.py pool runs/farfield_strong_pool4c.txt --t 4 $COMMON \
  --exact runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4a_e1p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4a_e2p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4a_e5p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e5p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e10p0_c3.0_exact_support.txt \
  --rot 1.05,-1.05,1.2,-1.2,1.5,-1.5,2.05,-2.05,2.3,-2.3,3,-3,5.05,-5.05,5.5,-5.5,7,-7,\
10.05,-10.05,11,-11,14,-14                                                               # 4,851

# the band-edge pool: every pose of the optimum KEPT WHERE IT IS and re-angled to tilt A (--tilt)
python3 search/farfield_strong.py pool runs/farfield_strong_pool4f.txt --t 4 $COMMON \
  --exact runs/arch_farfield_P4a_e5p0_c3.0_exact_support.txt \
  --exact runs/arch_farfield_P4c_e10p0_c3.0_exact_support.txt \
  --rot 5.05,-5.05,6,-6,8,-8,10.05,-10.05,12,-12 \
  --tilt 5.05,-5.05,5.3,-5.3,6,-6,7,-7,10.05,-10.05,10.5,-10.5,12,-12                     # 4,871

# ---- THE CERTIFICATES: the capped LP of the strong family, one cell per (eps, cap) ------------
# phase 1 row-generates on the pool; the polish then restricts to the support and closes the rows
# exactly (LEAF_CEILING.md 5.4) -- on a pool this size the unpolished loop sits in a degenerate
# tail with M oscillating just above 1 and never converges.  Detached, ~20 min per cell.
C4="--t 4 --corners 1111 --chord --resume-pgons runs/inputs-2026-09-12/cl_E2Bg_pgons.json \
    --seed-rows search/pgonly_corner_exact.txt --rows0 16000 --row-cap 10000 \
    --iters 8 --polish-iters 30 --time 900 --threads 2 --procs 1 --lp-tlim 0"
setsid nohup python3 search/farfield_strong.py cap S4A --pool runs/farfield_strong_pool4b.txt \
    $C4 --uncapped --eps 5 10 --cap 3 --recheck > runs/farfield_strong_S4A.out 2>&1 &
setsid nohup python3 search/farfield_strong.py cap S4B --pool runs/farfield_strong_pool4b.txt \
    $C4 --eps 15 20 30 --cap 3            > runs/farfield_strong_S4B.out 2>&1 &
setsid nohup python3 search/farfield_strong.py cap S4D --pool runs/farfield_strong_pool4c.txt \
    $C4 --uncapped --eps 1 2 5 10 --cap 3 > runs/farfield_strong_S4D.out 2>&1 &   # the headline
setsid nohup python3 search/farfield_strong.py cap S4F --pool runs/farfield_strong_pool4c.txt \
    $C4 --eps 10 --cap 1 2 4 9            > runs/farfield_strong_S4F.out 2>&1 &   # the ladder
setsid nohup python3 search/farfield_strong.py cap S4H --pool runs/farfield_strong_pool4f.txt \
    $C4 --uncapped --eps 5 10 --cap 3     > runs/farfield_strong_S4H.out 2>&1 &   # band edge

# ---- the T = 3 control (points + odd polygons is pool-limited; points only is the control) -----
python3 search/farfield_strong.py pool runs/farfield_strong_pool3b.txt --t 3 \
  --exact runs/arch_farfield_P3a_e0p5_c1000000.0_exact_support.txt \
  --exact runs/arch_farfield_P3a_e5p0_c2.0_exact_support.txt \
  --exact runs/arch_farfield_P3a_e10p0_c2.0_exact_support.txt \
  --exact runs/arch_farfield_P3a_e20p0_c2.0_exact_support.txt \
  --exact runs/arch_farfield_P3a_e30p0_c2.0_exact_support.txt \
  --exact runs/arch_farfield_P3b_T3e5_e5p0_c2.0_exact_support.txt \
  --exact runs/arch_farfield_P3b_T3e10_e10p0_c2.0_exact_support.txt \
  --rot-src runs/arch_farfield_P3a_e0p5_c1000000.0_exact_support.txt \
  --rot 3,-3,7,-7,12,-12,18,-18,25,-25,35,-35                                             # 4,132
setsid nohup python3 search/farfield_strong.py cap S3C \
  --pool runs/farfield_strong_pool3b.txt --t 3 --uncapped --eps 5 10 15 20 30 --cap 2 --recheck \
  --corners '' --seed-rows runs/arch_farfield_P3a_e0p5_c1000000.0_exact_support.txt \
  --rows0 12000 --row-cap 6000 --iters 8 --polish-iters 30 --time 600 --threads 2 --procs 1 \
  --lp-tlim 0 > runs/farfield_strong_S3C.out 2>&1 &
# with `--pent 5,7,9` (odd polygons separated at t = 3) the pool's UNCAPPED value is 5.965582893
# < 6, i.e. pool-limited: that variant is `S3B` and is not read as a result.

# ---- the checks (7b), the tilt profiles, and the table ----------------------------------------
python3 search/arch_farfield.py hist search/farfield_strong_S4A_uncapped_exact.txt \
        search/farfield_strong_S4D_e1p0_c3p0_exact.txt \
        search/farfield_strong_S4D_e10p0_c3p0_exact.txt --eps 1 2 5 10 15 20 30 --cap 3
python3 search/farfield_strong.py table S4A S4B S4C S4D S4E S4F S4G S4H S3C
```

All nine runs (`S4A`–`S4H`, `S3C`) finished inside their clocks; every cell reports
`certified=True`.  No detached job and no watcher process was left running at hand-off.
