# arch-farfield: does a near-axis mass cap close the far field?  Lemma **A3** measured at `T = 3` and `T = 4` (2026-09-20)

*Task: `tasks/arch-farfield/README.md`, plus the coordinator's addendum of the same evening (sweep the
cap over `{T-1, T, (T-1)^2}`, not only `T-1`).  Code: `search/arch_farfield.py` (new).  It **imports**
`search/packing_dual.py` (C kernels, `Model`, the exact sweep pricer, the refiner, the
arrangement-vertex certification) and `search/dual_exact.py` (rational pose snapping, exact
arrangement vertices, exact coverage) and modifies nothing.  All runs are in the gitignored
`runs/arch_farfield_*`, so every number is quoted here.  Semantics as everywhere in this repo
(`notes/s13-casefree.md` §1): closed unit squares, closed containment, closed container.*

---

## 0. Verdict, up front

**A3's candidate fails, and it fails at the small `eps` the composition actually needs.**  Capping the
near-axis mass at `T-1` costs the `t = T` fractional-packing / cover relaxation essentially nothing:

> **`V(4, eps, 3) >= 12098614247/1000000000 = 12.098614247 > 12`** for `eps = 0.5°` and `eps = 1°`,
> and `>= 12.085541166` for `eps = 2°` — each by an **exact certificate** (an explicit rational
> D4-symmetrised measure whose near-axis mass is exactly `3` and whose maximum coverage `M` is exactly
> `3999999999/4000000000 < 1`, re-checked independently over all `303,860` arrangement vertices of the
> **whole** container by `search/dual_exact.py check --full`).
>
> **`V(3, eps, 2) >= 6.080714, 6.073422, 6.010225, 6.017318` at `eps = 0.5°, 1°, 2°, 5°`** — the same
> one container down, where `Claim(3)` is a theorem (Kearney–Shiu).  So the mechanism does not produce
> the conclusion even in a case where the conclusion is known to hold.
>
> At `T = 3` the cap only starts to bite past `eps = 10°` (§3.4), which is far outside any tube radius
> the local lemmas A2b/A2c could supply, so it is of no use to the composition even there.

On the *same* pose pool the **uncapped** value is `12.101126115`.  So at `T = 4`:

| `eps` | `0.5°` | `1°` | `2°` | `5°`, `10°` (pool `P4c`) |
|---|---|---|---|---|
| cost of `mu(tilt < eps) <= 3` | **`0.002512`** | **`0.002512`** | `0.015585` | `0.165895` |

**Two thousandths of a unit.**  The missing `0.2688` between the pure relaxation and `12` is not
bought by near-axis mass at all, and the picture in `notes/status.md` ("at least `2.6` of the missing
`3` is bought by tilts under one degree") is about the *complementary* restriction — poses **confined**
to the band — not about the *cap*, which is what A3 proposes.  The LP's answer to "you may keep only
`3` units of near-axis mass" is one move: **tilt the corner square by `1.2°`**.  In the certified
measure the single heaviest pose is the corner square at `(0.51036, 0.51036)` at `+1.1997°` carrying
`1.81075`, and the rest of the configuration is unchanged (§4).

The result is robust to every *other* sound row the repo has: the certified measure satisfies the K1
band rows (`max_a mu(band_a) = 2.948652 <= 3`), the corner-box rows (`0.845955 <= 1` on each of the
four corners) and the wall-strip rows (`2.948652 <= 3`).  It was **not** tested against the odd-polygon
rows (§6, caveat 2).

**What is *not* settled at `T = 4`: the cells `eps >= 5°`, and the reason is provably the pose pool.**
The best pool there (`P4c`, 248 poses) gives `V(4, 5°, 3) >= 11.739389` and `V(4, 10°, 3) >= 11.739389`
— but **that same pool's own uncapped value is `11.905285`**, already below `12` and `0.364` below the
certified truth `12.268804`.  The cap costs it `0.166`; the pool costs it `0.364`.  A cell below `12`
whose pool is below `12` *with no cap at all* measures the pool, so nothing is concluded there in
either direction.  The `T = 3` rehearsal shows what happens when the pool is fixed: the same `eps = 5°`
cell read `5.972535` on a thin pool and `6.017318` on the pool generated for it — it crossed `n` purely
by getting better columns.  §5 also corrects the addendum's sanity cell: a margin-`0` configuration is
*not* a closed packing, and neither it nor its tilt-average is worth `12` to this relaxation, so `< 12`
at `cap = 4, eps = 5°` is not by itself a bug.

**The cap the composition actually needs — `(T-1)^2 = 9` — is dead everywhere tested**, including
`eps > 12.8°`: the uncapped optimum never uses more than `9.02` of near-axis mass, so the row
`mu(tilt < eps) <= 9` is essentially never binding (`V(4, eps, 9) >= 12.101126` for `eps <= 15°`,
`>= 12.070847` at `eps = 30°`).

---

## 1. The question, and which side of it every number is on

For a pose `S` write `tilt(S)` for the distance of its angle to `0 mod 90°`, and let

    V(T, eps, K)  =  sup  mu(Pi)   s.t.   mu({Q : p in Q}) <= 1  for every p in [0,T]^2
                                          mu({Q : tilt(Q) < eps}) <= K
                                          mu >= 0

A3 (`notes/proof-architecture.md` §2) proposes `V(T, eps, T-1) < n = T^2 - T`.

**Direction, which is the whole discipline of this note** (`notes/review-2026-09-13b.md`
"Clarifications", `search/T4SCREEN.md` §0):

* A measure supported on a **finite** pose set that satisfies the coverage rows and the cap is feasible
  for the continuum problem.  Its mass is therefore a **lower bound on `V`** — the packing side.
  Restricting the pose set can only lower the LP value.  **A value `>= n` refutes A3 at that `eps`;
  a value `< n` proves nothing at all.**
* An **upper** bound needs a weighted point set `y >= 0` and a scalar `lam >= 0` with `y(S) >= 1` for
  every admissible pose of tilt `>= eps` and `y(S) >= 1 - lam` for every near-axis pose; the bound is
  then `|y| + K lam`.  That has to hold for every pose in the continuum, not for a sample.  Nothing in
  this session produced one (§6).

Every headline number below is on the packing side, and the verdict rests only on the cells where that
side is `>= n`.

---

## 2. The instrument

`search/arch_farfield.py` has six modes; the two that carry the result are `cg` and `exact`.

**`cg T TAG --eps E --cap K`** — column generation for the capped LP.  It subclasses
`packing_dual.Model` (`CapModel`), whose only change is one extra row in `solve()`:
`sum_{tilt(pose) < eps} mu <= K`, read back as the multiplier `lam`.  Columns are priced by
`packing_dual`'s exact per-angle closed-square sweep pricer plus its refiner and neighbour moves, with
the reduced cost corrected to `1 - capture - lam*[near-axis]`; rows are the violated arrangement
vertices of the current support, from the same C kernel.  Seeds: `packing_dual.seed_poses`, extra
angles hugging the band edge, and the **rotated seeds** of the `rotseed` mode (below).  This is a
*proposal* engine: nothing it prints is quoted as a bound except through `exact`.

**`exact --src POOL --t T --eps E... --cap K...`** — the certificate.  It snaps every pose of the pool
to a rational `(p, q, cx, cy)` with `theta = 2 arctan(p/q)` (`dual_exact.snap_pose`, `Q = 10^6`),
forms the 8 dihedral images exactly, enumerates **every** vertex of the arrangement of the images in
the fundamental domain as an exact integer triple, computes the exact incidence of every vertex, and
solves

    max sum mu   s.t.  (coverage at every arrangement vertex) <= 1,  (near-axis mass) <= K

with HiGHS.  The masses are then rounded **down** to multiples of `10^-9`, the maximum coverage `M` and
the near-axis mass are recomputed as exact `Fraction`s, and if either exceeds its bound the measure is
scaled down exactly.  The printed value is the mass of a measure that is **exactly** feasible.  Because
the row set is the complete arrangement of the pool (not a point grid), this LP is the *exact optimum
of the capped problem restricted to that pose set* — which is why it beats the float `cg` value by a
long way (at `T = 3`: `5.53` float `-> 6.02` exact on the same poses).

Why the vertex set suffices is `dual_exact.py`'s own argument, unchanged: `cov` is a finite sum of
indicators of closed squares, hence upper semicontinuous; the intersection of the squares containing a
maximiser is a nonempty compact convex polygon on which `cov` is at least the maximum, and every vertex
of that polygon is a square corner or the crossing of two non-parallel closed edges.

**The near-axis test is exact.**  `sin(tilt) = min(|cos|, |sin|) = min(|q^2-p^2|, |2pq|)/(q^2+p^2)`,
and the test is `min(|q^2-p^2|, |2pq|) * 10^12 <= ceil(sin(eps) * 10^12) * (q^2+p^2)` — integers only.
`sin(eps)` is rounded **up**, so the set called "near-axis" is a **superset** of the true one and the
cap is, if anything, too strong: the certified value stays a valid lower bound on `V(T, eps, K)`.

**`rotseed --src S --alpha A`** — the far-field candidate that makes the whole thing work.  A rigid
rotation of the container about its centre preserves coverage **exactly** and shifts every angle by
`alpha`, so the sub-family of a known optimal measure that still fits in `[0,T]^2` after rotating by
`alpha` is a feasible measure all of whose originally-axis-parallel mass now sits at tilt exactly
`alpha`.  At `t = 4`, rotating `runs/dual_E1_support.txt` (the `12.2688` support of `COVER4.md`) keeps
`60.8 %` of its mass at `alpha = +0.6°` rising to `63.8 %` at `+11°` — not enough on its own, but the
right *columns*, and the LP then refills the corners it freed.

**`hist` / `bands`** — diagnostics: the tilt histogram and the "restrict-and-refill" bound
`mass(tilt >= eps) + min(near, K)`; and, for a measure, `max_a mu(band_a)` over the K1 bands
`D(t) = (|cos t| + |sin t|)/2 - |cos t sin t|` (`search/S6_SKELETON.md` §2.1: at most `T-1` squares of
a packing share a point of their `D`-bands, so `mu(band_a) <= T-1` is a **sound** row), plus the
corner-box and wall-strip masses.

**`dualcost`** — the cover-side cost of a run's dual, with an angle-**sampled** minimum capture.  A
heuristic, never a certificate.

---

## 3. Results

### 3.1 Why the existing measures say nothing by themselves

The optimum of the uncapped problem is overwhelmingly near-axis, so simply restricting it is hopeless
(`arch_farfield.py hist`):

| measure | mass | `tilt = 0` exactly | `tilt < 0.5°` | restrict-and-refill at `eps = 0.5°`, `K = 3` |
|---|---|---|---|---|
| `search/cover4_exact_support.txt` (`COVER4.md`) | `12.268804` | `6.538698` | `7.558070` | `7.710733` |
| `runs/dual_E1_support.txt` (its float parent) | `12.268806` | `5.532910` | `7.256146` | `8.012660` |
| `search/pgonly_corner_exact.txt` (`PGCORN`) | `12.000000` | `7.597978` | `9.150436` | `5.849564` |

So the `12.0986` below is **not** a restriction of a known measure; the LP had to find a different
fractional packing, and it did.

### 3.2 `T = 4`, `n = 12` — the target

Two pools, each a support of one round of a capped column generation whose seeds included the rotated
`E1` poses: `P4a` = `runs/arch_farfield_T4e1_support.txt` (191 poses, `166,897` arrangement vertices in
the fundamental domain, `108,818` distinct incidence rows) and `P4c` = `runs/arch_farfield_T4e10_support.txt`
(248 poses, `148,943` rows).  Every entry is an exact certificate, hence a **lower bound** on
`V(4, eps, cap)`; the table shows the better of the two pools per cell.

| `eps` | cap `3` | cap `4` | cap `9` | cap `none` | pool |
|---|---|---|---|---|---|
| `0.5°` | **12.098614** | **12.101126** | **12.101126** | **12.101126** | P4a |
| `1.0°` | **12.098614** | **12.101126** | **12.101126** | **12.101126** | P4a |
| `2.0°` | **12.085541** | **12.095866** | **12.101126** | **12.101126** | P4a |
| `5.0°` | 11.739389 | 11.826920 | **12.101126** | **12.101126** | P4a, P4c |
| `10.0°` | 11.739389 | 11.826920 | **12.101126** | **12.101126** | P4a, P4c |
| `15.0°` | 11.701882 | 11.797305 | **12.101126** | **12.101126** | P4a, P4c |
| `20.0°` | 11.593434 | 11.703149 | **12.101101** | **12.101126** | P4a, P4c |
| `30.0°` | 11.473717 | 11.609875 | **12.070847** | **12.101126** | P4a, P4c |

**Calibration — read this before reading any cell below `12`.**  Each pool's own **uncapped** value
(the last column of the table, and here per pool) says how good that pool is; the certified truth is
`nu_f^closed(4) >= 12.268804` (`COVER4.md`).

| pool | poses | exact rows | its **uncapped** value |
|---|---|---|---|
| `P4a` | 191 | 108,818 | `12.101126` |
| `P4c` | 248 | 148,943 | `11.905285` |

So `P4a` is `0.168` short of the truth before any cap is applied and `P4c` is `0.364` short.  A cell
below `12` whose own pool is already below `12` uncapped (every `eps >= 5°` cell of `P4c`) measures the
pool.  The cells that decide the question are the ones **above** `12`, and they are exact certificates.

**The cost of the cap, on a fixed pool** (`cap = 3` against the same pool's uncapped value):

| `eps` | `0.5°` | `1°` | `2°` | `5°` (`P4c`) | `10°` (`P4c`) |
|---|---|---|---|---|---|
| cost of `mu(tilt < eps) <= 3` | `0.002512` | `0.002512` | `0.015585` | `0.165895` | `0.165895` |

Even at `eps = 10°` the cap costs `0.166` of a relaxation whose distance to `12` is `0.2688`.  The
shadow price of the cap row in the `eps = 1°` column generation at round 1 was `lam = 0.0247`.

### 3.3 `T = 3`, `n = 6` — the rehearsal, where `Claim(3)` is a theorem

Pools: one per `eps`, each the support of that `eps`'s column generation, plus the joint pool `P3a`.

| `eps` | cap `2` | cap `3` | cap `4` | cap `none` | pool |
|---|---|---|---|---|---|
| `0.5°` | **6.080714** | **6.081236** | **6.081236** | **6.081236** | P3b_T3e0p5 |
| `1.0°` | **6.073422** | **6.074460** | **6.074460** | **6.074460** | P3b_T3e1 |
| `2.0°` | **6.010225** | **6.026468** | **6.026854** | **6.026854** | P3a |
| `5.0°` | **6.017318** | **6.039755** | **6.047619** | **6.047619** | P3b_T3e5 |
| `10.0°` | 5.951334 | 5.986395 | 5.993703 | **6.026854** | P3a, P3b_T3e10 |
| `15.0°` | 5.617647 | 5.835732 | 5.966679 | **6.026854** | P3a |
| `20.0°` | 5.603352 | 5.812104 | 5.962794 | **6.026854** | P3a |
| `30.0°` | 5.540000 | 5.796016 | 5.959502 | **6.026854** | P3a |

Calibration, same rule (the truth at `t = 3` is not certified in this repo; the best uncapped value any
pool here reached is `6.081236`):

| pool | poses | exact rows | its **uncapped** value |
|---|---|---|---|
| `P3a` | 151 | 121,649 | `6.026854` |
| `P3b_T3e0p5` | 127 | 114,806 | `6.081236` |
| `P3b_T3e1` | 127 | 123,155 | `6.074460` |
| `P3b_T3e10` | 162 | 176,392 | `5.986395` |
| `P3b_T3e2` | 33 | 5,617 | `5.980769` |
| `P3b_T3e5` | 158 | 150,503 | `6.047619` |

`P3b_T3e2` has only 33 poses and is `0.10` below its siblings *uncapped* (`5.980769 < 6`), so on its own
it would have put the `eps = 2°` row below `6`; the joint pool `P3a` settles that cell at `6.010225`.
This is the calibration rule in action, and it is why the sub-`n` cells at `T = 4` are not read as
evidence.  From `eps = 10°` on, by contrast, `P3a` is above `6` uncapped while every capped cell is
below — there the cap really is doing the work (§3.4).

**At `T = 3`, `cap = T-1 = 2`, every `eps` from `0.5°` to `5°` is certified `>= 6`.**  `Claim(3)` is a
theorem (Kearney–Shiu), so A3 is not merely unproven there — the mechanism it proposes demonstrably
does not produce the conclusion in a case where the conclusion holds.

The column generation reached the same verdict independently in floats, converging to `M = 1.000000`:
`runs/arch_farfield_T3e{0p5,1,5}.out` end with `BEST ... L = 6.080786, 6.073867, 6.017318`, all `>= 6`.

### 3.4 Where the cap *does* start to bite, and why it is the wrong place

At `T = 3` one pool (`P3a`) covers every `eps` and is good enough to exceed `n` **uncapped** at all of
them (`6.026854 > 6`).  On that fixed pool the cap's cost is measurable across the whole range:

| `eps` | `0.5°` | `1°` | `2°` | `5°` | `10°` | `15°` | `20°` | `30°` |
|---|---|---|---|---|---|---|---|---|
| `V(3, eps, 2)` on `P3a` | `6.024370` | `6.020355` | `6.010225` | `5.972535` | `5.830856` | `5.617647` | `5.603352` | `5.540000` |
| cost of `mu(tilt < eps) <= 2` | `0.002484` | `0.006499` | `0.016629` | `0.054319` | `0.195998` | `0.409207` | `0.423502` | `0.486854` |

(The main table's `eps = 10°` entry, `5.951334`, comes from a richer pool built for that `eps`; the row
above is a single fixed pool so that the *cost of the cap* is comparable across `eps`.)

So the cap is worth thousandths below `2°`, hundredths at `5°`, and only becomes a real constraint
somewhere past `10°` — and even at `30°` it takes `V` only `0.49` below `6`, i.e. it would have to be
combined with something else.  **This is the wrong place for A3 to work.**  A5 (the glue,
`notes/proof-architecture.md` §2) needs `eps` no larger than the radius of the A2b/A2c tube, and the
relevant scales there are the `~2/T` straight-chain window and a second-order deficit `-(1/3) t^2`; an
`eps` of `10°–30°` is far outside any tube those lemmas can supply, and A4 is stated for tilts
`< eps ~ 1/T`.  A3 has to hold at small `eps` or it is of no use to the composition, and at small `eps`
it is refuted above.

At `T = 4` the same transition is **not resolved**: the only pools that reach past `eps = 5°` are
`P4c` (uncapped `11.905285`, below `n`) and `P4a` (built for `eps = 1°`).  Nothing here says whether
`V(4, 10°, 3)` is above or below `12`.

---

## 4. The fractional packing that does it

`runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt`: **89 poses, mass `12098614247/10^9`, near-axis
mass `2999999999/10^9 <= 3`, `M = 3999999999/4000000000 < 1`.**  Re-derived from the file alone, with
no LP and integers only, by the repo's own checker over the **whole** container:

```
python3 search/dual_exact.py check runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt \
        --t 4 --tag archff_P4a_e0p5_c3full --full --procs 2
  303860 vertices checked; mass = 12098614247/1000000000;  M = 3999999999/4000000000
```

Tilt profile (`arch_farfield.py hist`), and the K1 / region rows (`arch_farfield.py bands`):

| tilt band | mass | | sound row | value | allowed |
|---|---|---|---|---|---|
| `= 0` exactly | `2.794606` (3 poses) | | `max_a mu(band_a)` (K1) | `2.948652` | `3` |
| `(0, 0.5°)` | `0.205394` (1 pose) | | corner box `[0,1]^2` (x4) | `0.845955` | `1` |
| **near-axis total** | **`3.000000`** (cap binding) | | wall strip `[0,1] x [0,4]` (x4) | `2.948652` | `3` |
| `[1°, 2°)` | `3.167289` (12) | | centre in the open interior | `3.687827` | — |
| `[2°, 5°)` | `1.560598` (13) | | | | |
| `[5°, 10°)` | `1.007819` (6) | | | | |
| `[10°, 20°)` | `0.261032` (5) | | | | |
| `[20°, 30°)` | `1.289637` (17) | | | | |
| `[30°, 40°)` | `1.086645` (24) | | | | |
| `[40°, 45°]` | `0.725593` (8) | | | | |

Heaviest poses (`cx cy theta mass`; the measure is D4-symmetrised, so each carries `mass/8` on each of
its 8 dihedral images).  50 % of the mass is on 9 poses, 90 % on 50:

```
 0.51036  0.51036   +1.1997 deg   1.81075     <- the corner square, tilted
 0.50000  0.50000    0.0000 deg   1.57307     <- the corner square, axis-parallel   ]
 0.50000  1.50000    0.0000 deg   1.00146     <- a wall square, axis-parallel       ] near-axis, total 3
 1.50000  1.50000    0.0000 deg   0.22007     <- an interior square, axis-parallel  ]
 0.50175  2.49838   +0.2000 deg   0.20539     <- the 4th near-axis pose             ]
 0.52006  1.48027   +2.2999 deg   0.36481
 0.51715  1.52000   +2.0002 deg   0.30988
 1.22511  2.53466  +40.5001 deg   0.26508
 1.50395  2.48815   +1.7000 deg   0.26508
```

**Read it as one sentence.**  The uncapped optimum of `COVER4.md` puts `6.54` of its `12.2688` at tilt
exactly `0` — four poses, the corner square being the heaviest.  Told it may keep only `3`, the LP
rotates the corner square by `1.2°` and leaves everything else where it was.  A tilt of `1.2°` costs a
wall-to-wall chain of four `4 cos t + sin t - 4 = +0.020060` of width — a chain of four at that tilt
does **not** fit (it only fits again from `28.1°` on, `notes/proof-architecture.md` §3.2) — and the
relaxation pays for that with `0.0025` of mass, because it is a *measure* and can smear the discrepancy over a cloud of
poses — the item-3 family of `search/HONEST.md` §0 again, one degree up.  That is why the mechanism
fails: the near-axis mass in the optimum is an artefact of where the pose lattice happens to put it,
not a structural commitment of the relaxation.

---

## 5. The coordinator's sanity cell, and what it proves about the pools

The addendum supplies a built-in sanity check: *for any cap `>= 4` and `eps <= 12.8°` the value MUST be
`>= 12`, because a real closed packing exists there — four axis-parallel squares in a wall-to-wall chain
beside eight squares tilted up to `12.8°`, margin exactly `0` (`search/S6_LOCAL.md` §2, `n = 12` table,
`j = 4`); if the LP reports `< 12` at `cap = 4, eps = 5°` that is a bug or a pose-set restriction.*

The best pool reports `V(4, 5°, 4) >= 11.826920`, i.e. it fails the cell by `0.173`.  **It is a
pose-set restriction, and in addition the cell is not the theorem the addendum takes it for.**  Both
halves matter, so both are spelled out.

**(a) The cell is not a theorem, because margin `0` is not a closed packing.**  `delta* = 0` means the
best configuration at those angles has *some pair touching*: the squares have disjoint interiors but
are not pairwise disjoint as closed sets.  The relaxation's rows are coverage rows for **closed**
squares, so at a contact point of two squares the indicator measure of that configuration has
coverage `2`, and it is **infeasible**.  Nor does smearing repair it in general: with four axis-parallel
unit squares forming a *wall-to-wall* chain in width `4` there is no slack, so the three contacts
persist for every tilt `t` of the other eight, and the three coverage rows at those contacts force
`m_1 + m_2 <= 1`, `m_2 + m_3 <= 1`, `m_3 + m_4 <= 1` — the chain of four can carry at most **`2`**, not
`4`, however the family is averaged over `t`.  The indicator (or `t`-average) of the `j = 4` family is
therefore worth at most `12 - 2 = 10` to this LP, not `12`.  So "a real closed packing exists there"
does not transfer to the `t = T` closed-coverage relaxation, and `< 12` at `cap = 4, eps = 5°` is not
*per se* a bug.  (At `t < T`, where a packing has positive margin, the implication would be sound.)

**(b) The pools are nevertheless what is short at `eps >= 5°`, and it is visible in one number.**  The
pool that covers `eps >= 5°` (`P4c`) has uncapped value **`11.905285`** — *below `12` with no cap at
all*, against the certified truth `12.268804`.  Its shortfall, `0.364`, is more than twice what the cap
costs it (`0.166`).  Every pool in this session was column-generated around a *small*-tilt optimum (its
seeds are `packing_dual`'s lattice plus the `COVER4` support rotated by `0.6°…11°`), and at `eps >= 5°`
the cap forces the measure into a regime those seeds never explored.  The `T = 3` rehearsal shows the
repair working: the `eps = 5°` cell read `5.972535` on a pool built for other `eps` and `6.017318` on
the pool generated for it, crossing `n` purely by acquiring better columns.

**Consequently:** the cells `eps >= 5°` are reported as lower bounds and **no conclusion is drawn from
them in either direction.**  The verdict of §0 rests only on `eps <= 2°`, where the certificates are
above `n`.

---

## 6. Caveats, honestly

1. **Only the cells `>= n` are conclusive.**  `V(4, 5°, 3) >= 11.739389` does *not* say A3 might hold
   at `eps = 5°`; §5(b) says why the number is about the pool (that pool is below `12` uncapped).
   Deciding `eps >= 5°` at `T = 4` needs either a column generation seeded in the `5°–13°` regime and
   run to `M = 1` (a few hours — the three `cg` runs left going are a start), or a cover-side
   certificate.  At `T = 3`, where the pools are good enough, the cap does begin to bite past `eps = 10°`
   (§3.4) — but at an `eps` no tube can reach.
2. **A3 is stated for "points + polygons", and this note measures points only.**  The LP here has the
   *complete* coverage-row family (every arrangement vertex of the support) but no odd-polygon rows and
   no clique rows, so it is the **pure** relaxation whose value is `nu_f^closed(4) >= 12.2688`
   (`COVER4.md`); the points+polygons family is certified at `12.173398` and the corner `k = 4` family at
   `11.999999926 <= 12` (`notes/status.md`).  What is shown here is that the *cap* costs `0.0025` at
   `eps <= 1°`, two orders of magnitude less than the `0.2688` a polygon/region family would have to
   find; and the certified measure satisfies every region and K1 row that was checked (§4).  It was not
   checked against the `1,919` odd-polygon rows of `runs/inputs-2026-09-12/cl_E2Bg_pgons.json`
   (`rankdiag.py --pgons` reads a `leaf_ceiling` measure file, and the format conversion was not done).
   **If someone wants to keep A3 alive, that is the one remaining crack** — and note that it must find
   `0.2688` of slack that the cap is not providing.
3. **`cap = 9 = (T-1)^2` — the cell the composition actually needs — is dead by a wide margin at every
   `eps` tested, including `eps > 12.8°`.**  The uncapped optimum of the pool never uses more than
   `9.02` of near-axis mass even at `eps = 20°`, so the row `mu(tilt < eps) <= 9` is essentially never
   binding: `V(4, eps, 9) >= 12.101126` for `eps <= 15°`, `>= 12.101101` at `20°` and `>= 12.070847` at
   `30°`.  This closes the addendum's "open cells `cap in {4..9}` at `eps > 12.8°`" for `cap = 9`.
4. **No cover-side number was obtained.**  §6 of the brief asks for one separately; the only dual in
   hand is from an unconverged round-0 LP and costs `16.94` honestly (`dualcost`, angle-sampled), which
   certifies nothing.  Since the packing side is already `>= 12` at the cells that matter, an upper
   bound there could only confirm the negative.
5. **The variant "no wall-to-wall unit strip holds `T` near-axis squares" is expressible but redundant.**
   It is the row `mu({near-axis poses with |a - c_y| <= 1/2}) <= T - 1`, which is implied by K1
   (`mu(band_a) <= T-1` for *all* poses, `search/S6_SKELETON.md` §2.1), a row that is already sound
   without any near-axis hypothesis.  The certified measure of §4 satisfies the stronger row with room
   to spare (`2.948652 <= 3`), so the variant does not cut it either.

---

## 7. What is exact and what is float

**Exact** (Python `int` / `Fraction`; no float decides anything):

* the rational poses `(p, q, cx, cy)`, their 8 dihedral images, and the admissibility of every image
  (four rational corners inside the closed container);
* every arrangement vertex, as a reduced integer triple `(X, Y, D)`, and every incidence
  (`2|a dx + b dy| <= r D Dk`);
* the mass, the maximum coverage `M` and the near-axis mass of the certified measure, all as
  `Fraction`s over the common denominator `10^9`; the LP solution is rounded **down** before any of
  them is computed, and the measure is scaled down exactly if `M > 1` or the cap is exceeded;
* the near-axis test `min(|q^2-p^2|, |2pq|) * 10^12 <= ceil(sin(eps) 10^12) (q^2+p^2)`, with `sin(eps)`
  rounded **up** so the near-axis set is a superset (the cap is if anything too strong);
* the independent re-check by `search/dual_exact.py check --full`, which re-reads the support file and
  re-derives `mass` and `M` from scratch.

**Float** (chooses *which* measure gets certified, and is never quoted as a bound):

* the whole `cg` mode — HiGHS, the sweep pricer, the refiner, row/column ageing.  Its `L` column is a
  crude scaling of an unconverged iterate (`L = 10.505` where the exact LP on the *same poses* gives
  `12.0986`), so it is quoted only as trajectory;
* the HiGHS solve inside `exact` that picks the masses on the exact rows;
* the pose snapping itself (`Q = 10^6` for the angle, `10^8` for the centre) — the snapped pose is a
  different, but exactly admissible, pose, so nothing is lost;
* `dualcost`'s angle-sampled minimum capture.

**Direction discipline.**  A restricted pose set can only lower the LP value, and the near-axis set is
over-estimated, so every certified number is a **lower bound** on `V(T, eps, K)`.  A cell `>= n`
refutes A3 there; a cell `< n` is silent.

---

## 8. Reproduce

```sh
# tilt histograms of the measures already in the repo, and the restrict-and-refill bound
python3 search/arch_farfield.py hist search/cover4_exact_support.txt runs/dual_E1_support.txt \
                                     search/pgonly_corner_exact.txt --eps 0.5 1 2 5 10 --cap 3

# the far-field columns: rigid rotations of the 12.2688 support (coverage is preserved exactly)
python3 search/arch_farfield.py rotseed --src runs/dual_E1_support.txt --t 4 \
        --alpha 0.6 1.1 2.2 5.5 11 --out runs/arch_farfield_rot4.txt

# capped column generation (proposal engine; detached, ~2 h)
setsid nohup python3 search/arch_farfield.py cg 4 T4e1 --eps 1 --rounds 80 --time 5400 \
        --polish-time 1800 --threads 1 --seed-pitch 0.08 --row-pitch 0.03 --cg-want 1500 \
        --seed-file runs/arch_farfield_rot4.txt runs/dual_E1_support.txt \
        > runs/arch_farfield_T4e1.out 2>&1 &

# THE CERTIFICATES: exact capped LP over the complete arrangement of a pose pool
python3 search/arch_farfield.py exact --src runs/arch_farfield_pool4.txt --t 4 \
        --eps 0.5 1 2 5 10 15 20 30 --cap 3 4 9 1000000 --tag P4a --Q 1000000 --procs 2
python3 search/arch_farfield.py exact --src runs/arch_farfield_pool3b.txt --t 3 \
        --eps 0.5 1 2 5 10 15 20 30 --cap 2 3 4 1000000 --tag P3b --Q 1000000 --procs 2
python3 search/arch_farfield.py table P4a P3b

# independent exact re-certification of the headline measure (no LP, integers only)
python3 search/dual_exact.py check runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt \
        --t 4 --tag archff_P4a_e0p5_c3full --full --procs 2

# the sound rows the certified measure is tested against
python3 search/arch_farfield.py bands runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt --t 4
python3 search/arch_farfield.py hist  runs/arch_farfield_P4a_e0p5_c3.0_exact_support.txt --eps 0.5 --cap 3

# the cover side (heuristic, angle-sampled: NOT a certificate)
python3 search/arch_farfield.py dualcost runs/arch_farfield_T4e1_support.txt --eps 1 --cap 3 --threads 2
```

Detached runs left behind at hand-off: the `cg` runs `T4e0p5`, `T4e1`, `T4e10` (each on a `5400 s + 1800 s`
clock, checkpointing `runs/arch_farfield_<TAG>_{last_,}support.txt` every round).  They can only add
columns; nothing in §0 depends on them.  No watcher processes were left running.
