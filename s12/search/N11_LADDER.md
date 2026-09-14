# `n = 11` on the same footing as `n = 12`: every rung of the machinery, run once (2026-09-14)

Task: `tasks/s11-ladder/README.md`.  Read against `search/N11.md` (the pure certificate),
`search/N11_ANATOMY.md` (the mass-11 packing and the rigorous pure ceiling), and, for each rung's
`n = 12` counterpart, `search/CEILING.md`, `search/CLIQUE_CEILING.md`, `search/CLIQUELEVER.md`,
`search/BRANCH.md` / `search/BENTZ.md`, `search/BOXCLIQUE.md`, `search/ANCHOR.md`, `search/RANK8.md` §4.
New code: `search/n11_sa2.py`.  New exact files: `search/n11_qstab_{382,383375,385,386,387}_exact.txt`.
Runs in the main tree's gitignored `runs/` (`cc_N11Q*`, `cc_N11L385`, `branch_N11AQ3980*`,
`boxclique_N11BC3980*`, `n11_sa2.log`); every number is quoted here.

---

## 0.  Why, and the verdict

Every comparison of the two problems so far compared the *pure* method at `n = 11` with the *fully
instrumented* method at `n = 12`.  At `n = 12` the certified record `3.968616` is pure-LP, but the
ceiling was then measured with every stronger family — cliques, the corner branch, polygons /
regions / chord, anchor cliques on the cover side, Sherali–Adams level 2 — all as packing-side
measurements at `t = 4`.  At `n = 11` only the pure rung had been run, plus one clique diagnostic.
So "the required gap at `n = 11` is twice the one at `n = 12`" (`N11.md`) was a pure-vs-pure
statement, and "the clique family cannot reach the `0.75`" (`N11_ANATOMY.md` §7) had not been
measured.  This note runs the same rungs at `n = 11`.

**Verdict.**  No new lower bound: `s(11) ≥ 3040/797 = 3.814304` stands, and it stands for the same
reason `s(12) ≥ 3.968616` stands — every certifiable clique family is worth **exactly zero on the
cover side** at `n = 11` too (anchor cliques: 400 offered, 0 used, matched gain `0.000000`; box
cliques: 0 used).  On the packing side the clique lever is real and clean: the clique-strengthened
relaxation is **exactly certified at `32/3` on the whole window `3.83375 ≤ t ≤ 3.87`** (and `10` at
`3.82`), i.e. cliques remove one third of the missing unit, by cutting one third of the mass of one
interior pose and nothing else; no clique-feasible measure of mass `≥ 11` was found below Trump's
`3.877083`, so from the packing side the clique method's ceiling is not located — the same
situation as `n = 12` at `t = 4`.  The one rung where the two problems differ in kind is
**Sherali–Adams level 2**: on the 28-square extremal support it closes the whole unit at `n = 11`
(exact rational certificate, `SA-2 = SA-3 = α = 10`), where at `n = 12` it closes a tenth.

The baselined gap statements are in §8.

---

## 1.  The ladder

Conjectured values `s(12) = 4` (Bentz's `s(13) = 4` plus monotonicity is the only upper bound) and
`s(11) = 3.877083` (Trump 1979).  "Exact" = integers and `Fraction`s end to end; "float" = an LP
value with no bound in either direction unless stated.

| rung | `n = 12` | `n = 11` |
|---|---|---|
| 0. literature | Stromquist `2 + 4/√5 = 3.788854` (inherited by monotonicity) | Stromquist `3.788854` |
| 1. pure certificate | **`3.968616`** (1736 points, `11.9738`) | **`3.814304`** (680 points, `10.8147`) |
| 2. pure ceiling, rigorous | `< 3.99` (`ν_f(3.99) ≥ 12.008`, `DUAL_EXACT.md`) | `≤ 3.83375` (`ν_f = 11` exactly, `N11_ANATOMY.md`); and `COVER(3.8191) ≤ 11.00004` (§5, new) |
| 2′. pure window | `[3.968616, 3.99)`, width `0.021`; heuristic crossing `3.965–3.98` | `[3.814304, 3.83375]`, width `0.0195`; heuristic crossing `3.8153–3.8162` |
| 3. clique-strengthened (QSTAB) measure, exact, packing side | `11.364` at `3.99` (`cc_C99c`); at `t = 4` nothing converged in 40 core-h; fixed-pose leaf QSTAB `11.303–11.674` (`CLIQUELEVER.md`) | **`10` at `3.82`; `32/3` at `3.83375, 3.85, 3.86, 3.87`**, four poses each (§2) |
| 4. corner branch | `k = 4` leaf at `t = 4`: `11.9999999` with all degree-1 families; `11.674` under all cliques (fixed poses) | vacuous: the pure optimum already has corner mass exactly `4`; under cliques the leaf reproduces `32/3` with corner mass `4` (§3) |
| 5. box-clique certificate, cover side | never load-bearing (`BOXCLIQUE.md`) | round 0 at `3.8191`: 8 cliques priced, 0 used; separation `O(n²)` on a 28k-row dual, 1.8 h/round, not completed (§4) |
| 6. anchor-clique certificate, cover side | matched gain `10⁻⁴` on the `k = 4` leaf at `3.98` (`ANCHOR.md`) | matched gain **`0.000000`** over 40 rounds at `3.8191`; 400 cliques offered, 0 used (§5) |
| 7. SA-2 on the extremal support | `11.893` on 162 poses: `0.107` of the missing unit | **`10`** on 28 squares: the whole unit, exact certificate; SA-3 `= 10` (§6) |

---

## 2.  Rung 3: the clique-strengthened relaxation, exactly certified at `32/3`

`search/clique_ceiling.py --exact` (the `CLIQUE_CEILING.md` loop: points + every closed-intersection
clique, column generation on poses, staged; each converged stage's measure is snapped to rationals
and certified — coverage `≤ 1` at every arrangement vertex, max clique mass `≤ 1` by a complete
integer branch and bound), warm-started from the four-pose mass-11 measures of `N11_ANATOMY.md`,
flags as in `CLIQUE_CEILING.md` "Reproduce" with `--inner 400 --time 14400`, 2 threads each.

| `t` | stages | converged stage values (all `M = 1`, max clique `= 1`, B&B complete) | unconverged stages (LP at the cap; residual max clique) | **certified exact** |
|---|---|---|---|---|
| `191/50 = 3.82` | 9 | `9, 9, 10, 10, 10` | `10.06 (1.22), 10.15 (1.61), 10.23 (1.58), 10.35 (1.52)` — 22k → 86k cuts | **`10`** (`n11_qstab_382_exact.txt`) |
| `3067/800 = 3.83375` | 3 | `32/3, 32/3` | `10.72 (1.13)` | **`32/3`** (`n11_qstab_383375_exact.txt`) |
| `77/20 = 3.85` | 5 | `32/3, 32/3, 32/3` | `10.84 (1.19), 10.81 (1.18)` | **`32/3`** |
| `193/50 = 3.86` | 6 | `32/3, 32/3, 32/3, 32/3` | `10.89 (1.22), 10.88 (1.21)` | **`32/3`** |
| `387/100 = 3.87` | 5 | `32/3, 32/3, 32/3, 32/3` | `10.81 (1.17)` | **`32/3`** (`n11_qstab_387_exact.txt`) |

(The certified files carry `5333333333/500000000 = 10.666666666`, the exact rounding-down of `32/3`
at denominator `5·10⁸`.)

**What the measure is.**  At every `t` from `3.83375` to `3.87` the certified `32/3` measure is the
`N11_ANATOMY.md` four-pose object with one change: the mass of the near-axis interior pose is cut
from `1` to `2/3`.

    corner   (1/2, 1/2, 0°)            mass 4       ring: 8, unchanged
    wall     (1/2, ~1.51, 0°)          mass 4
    interior tilted  (20.5°–33.7°)     mass 2       unchanged
    interior near-axis (2.8°–5.0°)     mass 2/3     was 1

So the clique constraints — the eight-pose non-Helly clique of mass `1.25` that `N11_ANATOMY.md` §3
found — remove exactly one third of the missing unit, from exactly the pose family the rank-10
statement is about, and touch nothing else.  At `3.82` the certified measure is `8 + 1 + 1` on two
tilted interior poses at `18–19°`, and the pure LP there is exactly `11` with a 28-pose dual
(`N11_ANATOMY.md` §6), so `3.82` is a container where the clique relaxation is at least `0.65`
below the pure one on every pose set the loop built (the stalled LP values are upper bounds on QSTAB
of *their* pose sets).

**What is and is not established.**  Every certified value is a rigorous *lower* bound on the
continuum QSTAB value at its `t`: `QSTAB(t) ≥ 32/3` for `t ≥ 3.83375`, `≥ 10` at `3.82`.  Nothing
here is an upper bound on QSTAB in the continuum: the unconverged stages are restricted to their
pose sets, and pricing reports `min reduced cost = 0` at every converged stage (degenerate — no
strictly improving column, but 300 columns of zero reduced cost each time, after which the LP
climbs to `10.8–10.9` and the inner loop cannot cut the residual `1.13–1.22` clique in 400
iterations; the same degeneracy `CLIQUE_CEILING.md` §"Reading" describes).  A clique-feasible
measure of mass `≥ 11` at some `t < 3.877` would locate the clique method's ceiling below the
conjecture; none was found, exactly as none of mass `≥ 12` was found at `t = 4` for `n = 12`.  The
recurring exact `32/3` after every convergence, on pose sets of 300 to 1160 columns, is evidence —
not proof — that the continuum QSTAB is `32/3` throughout the window.

---

## 3.  Rung 4: the corner branch is vacuous under cliques as well

`clique_ceiling.py --kmass 4 --r 1 --margin 0` at `77/20` (the `k = 4` corner leaf: corner-box
mass pinned to exactly `4`).  `--margin 1e-6`, the `n = 12` setting, is *infeasible* here: with the
corner boxes saturated, corner mass `4` forces coverage exactly `1` at the corner squares' common
point, which a `1 − 10⁻⁶` row forbids — itself a statement of the vacuity.  Stages 0 and 1: value
`32/3`, `M = 1`, max clique `1`, corner mass `4.000000` — identical to the unconstrained run.  Stage
2 was still cutting at hand-off (LP `10.74`, residual clique `1.17` on cliques of 450+ poses, ~100 s
per separation).  At `n = 12` the `k = 4` leaf was a real restriction (corner mass `0.85` in the pure
optimum) and the branch's other leaves die; at `n = 11` the pure and the clique-strengthened optima
both sit *in* the `k = 4` leaf, so a corner-count tree has one leaf, the whole problem.

---

## 4.  Rung 5: box cliques on the cover side

`boxclique.py runs/n11_H3985.txt --Dp 3980 --n 11 --N 2000` (the `3.8143` certificate rescaled to
`760/199 = 3.8190955`, 88 orbits / 680 atoms, 24k lattice warm-start rows).  Round 0: points
`10.715469`, probe minimum `0.64`, 8 clique orbits priced, **0 used** — and 6423 s for the round,
against "seconds" at `n = 12`: `price_cliques` builds a dense closed-overlap graph over the D4 images
of every row with positive interior-point dual, ~10⁵ images on this row set, `O(n²)`.  Killed after
4 h in round 1.  A retry without the lattice warm start (`--no-warm --N 1000 --per-round 4`) was
running at hand-off (`runs/boxclique_N11BC3980b*`).  The measurement that matters — cliques used —
is already `0`, as at `n = 12`, and §5 makes the same measurement properly with point column
generation in the loop.

---

## 5.  Rung 6: anchor cliques on the cover side, with column generation — gain exactly `0`

`branch.py runs/n11_H3985.txt N11AQ3980 --k 4 --r 1 --lam-lo 0 --lam-hi 0 --n 11 --Dp 3980 --N 2000
--colgen 40 --cliques 40 --cq-want 40 --matched --topk 3` (`BRANCH_SOLVER=restricted`): the pure
cover LP with a vacuous trailer (`ANCHOR.md` "Deliverable 8"), point column generation *and* anchor
clique columns `K(p, A)` priced every round, and every round also solved without the clique columns
(the matched pair).  40 rounds, 904 s.

| round | total | pure (matched) | gain | cliques offered / used | probe min |
|---|---|---|---|---|---|
| 0 | `10.715469` | — | — | 40 / 0 (`ȳ(K) = 1.4085` = the point value) | `0.654` |
| 1 | `11.012070` | `11.012070` | `−0.000000` | 80 / 0 | `0.762` |
| 2 | `11.009962` | `11.009962` | `+0.000000` | 120 / 1 (weight `0.318`) | `0.902` |
| 3 | `11.000022` | `11.000022` | `+0.000000` | 160 / 0 | `0.803` |
| 10 | `11.000022` | `11.000022` | `+0.000000` | 276 / 0; separator finds `ȳ(K) = 1.000`, nothing violated | `0.890` |
| 20 | `11.000022` | `11.000022` | `+0.000000` | 400 / 0 | `0.971` |
| 30 | `11.000022` | `11.000022` | `+0.000000` | 400 / 0 | `0.997` |
| 40 | `11.000022` | `11.000022` | `−0.000000` | 400 / 0 | **`1.0000010`, 0 violated** |

Finalised: 304 points, total `11.0000400`, verified `min = 2500007/2500000` at `N = 2000` **and**
`N = 4000` — a valid cover, not a certificate (`≥ 11`).  Hence, rigorously, **`COVER(760/199) ≤
11.00004`**; with `COVER(3040/797) ≤ 10.8147` (the certificate) the pure cover value rises from
`10.81` to `11.00004` over `0.0048` of container, i.e. the pure method's death at `n = 11` is a
cliff, as `TIGHTEN.md` found at `n = 12`.  The clique columns are offered the same non-Helly mass
the packing side exploits (§2: the pure dual here is the mass-11 object whose interior violates an
eight-pose clique by `0.25`), and use none of it: `ANCHOR.md` §4's reason applies unchanged — an
anchor clique costs what its point costs and can only pay for `{S : A ⊆ S} \ P_p`, a strip of
width `~0.01` at the tight points.

---

## 6.  Rung 7: Sherali–Adams level 2 closes the whole unit on the support

`search/n11_sa2.py` (agent run, 11 s per support, `runs/n11_sa2.log`), the formulation of
`RANK8.md` §4 (`rank8.py sa2`) on the 28 distinct closed squares of the mass-11 measure at
`3067/800` and at `77/20` — which have the **same** closed-intersection graph, 108 edges, 270
disjoint pairs, 29 maximal cliques (four of them non-Helly).  Variables `y_s` and `y_{st}` on
disjoint pairs; base rows either the coverage rows at the 400 arrangement vertices (deduplicated) or
the 29 maximal cliques; each base row lifted by `x_u` and `1 − x_u` for every pose `u` outside it.

| support | `α` (exact B&B) | pure LP | clique LP | SA-2, point rows | SA-2, clique rows | SA-3, clique rows |
|---|---|---|---|---|---|---|
| `3.83375` | `10` | `11` | `32/3` | **`10`** | **`10`** | `10` |
| `3.85` | `10` | `11` | `32/3` | **`10`** | **`10`** | `10` |

The SA-2 and SA-3 values are **exact**: the positive duals were rounded to rationals and the
combination re-verified in `Fraction`s to dominate `Σ y_s` with right-hand side exactly `10`.  The
certificate is elementary and worth stating (`3.85`, point rows, all multipliers `1`): corners `≤ 4`;
one wall pair `≤ 1`; one non-Helly triple `{wall, wall, interior} ≤ 1`, lifted by a case split on the
interior square; and "the left half `≤ 2`" and "the right half `≤ 2`", each obtained by a single
case split on one wall square, under which the point-clique rows become rank rows.  Its dual lives
entirely on wall–interior contacts (`100 %` of the lifted mass for point rows, `96.4 %` for clique
rows), zero on corners, zero on wall–wall — the `n = 11` analogue of `RANK8.md` §4's contact-graph
structure.  The clique LP alone reaches `32/3`, the number §2 certifies in the continuum.

Against `n = 12`: SA-2 on the 162-pose leaf-A support closes `0.107` of the unit and the certificate
needs the interior four-cycle; here it closes all of it with two single-variable splits.  Caveats
for the comparison: both are fixed-support measurements, and the `n = 11` support is six times
smaller; on a refined pose set SA-2 will be weaker.  But it is the same instrument, and it reaches
the truth at `n = 11` and a tenth of it at `n = 12`.

---

## 7.  What is exact and what is float

**Exact.**  The five QSTAB measures (`search/n11_qstab_*_exact.txt`; rational poses, coverage at every
arrangement vertex, closed intersection of every image pair in integers, max clique by a complete
integer B&B — `clique_ceiling.py --check FILE` re-verifies each independently); the anchor-run
cover (`runs/branch_N11AQ3980.txt`, exact `i128` verifier at `N = 2000, 4000`), hence `COVER(760/199)
≤ 11.00004`; `α = 10`, the intersection graph, the maximal cliques, and the SA-2 / SA-3 values `10`
on both supports (rational dual certificates).

**Float / heuristic.**  Every unconverged stage value in §2; the claim that QSTAB `= 32/3` in the
continuum; the anchor run's matched gain (a difference of two LP values, `0.000000` to the printed
precision); the box-clique round-0 numbers; the pure LP `11` and clique LP `32/3` on the SA-2
supports (known exactly from `N11_ANATOMY.md` and §2 respectively, but the script's are HiGHS
floats).  Nothing about `s(11)` itself changes.

---

## 8.  The baselined comparison

Gaps to the conjectured value, side units, both `n`:

| | `n = 12` (to `4`) | `n = 11` (to `3.877083`) |
|---|---|---|
| pure certificate | `0.0314` | `0.0628` |
| pure ceiling, rigorous | gap `> 0.010` (`ν_f(3.99) ≥ 12`) | gap `≥ 0.0433` (`ν_f(3.83375) ≥ 11`) |
| pure crossing, heuristic | `0.02–0.035` | `0.061–0.062` |
| as a fraction of the side | `0.5–0.9 %` | `1.1–1.6 %` |

So the pure-vs-pure statement survives and is now rigorous in both directions: the pure method's
gap at `n = 11` is at least `0.043`, and at `n = 12` at most `0.031`.  "Twice" is fair.

Beyond the pure rung the two problems look alike, rung for rung:

* **Cover side (what can be certified):** at both `n`, box and anchor cliques buy `0`; the corner
  branch buys nothing at `n = 11` (vacuous) and, at `n = 12`, cannot close its `k = 4` leaf.  The
  certified records at both `n` are pure-LP and will stay so with these families.
* **Packing side (what the relaxations are worth):** at `n = 11` cliques certifiably remove `1/3`
  of the unit throughout the window and heuristically exactly that; at `n = 12` they remove
  `0.3–0.7` of the unit on the leaf pose sets at `t = 4` and a certified `0.64` at `3.99`
  (`12.008 → 11.364`).  Neither problem yields a clique-feasible measure at the conjectured side, so
  neither clique ceiling is located.
* **The one asymmetry:** SA-2 on the extremal support closes the unit at `n = 11` and a tenth at
  `n = 12`.  In the language of `notes/status.md`, the `n = 11` obstruction is a rank-10 *packing*
  statement with an elementary two-split proof on its support, where the `n = 12` obstruction is a
  rank-8 *incidence* statement of margin zero.  This is the reason the `n = 11` statement was
  proposed as the rehearsal for `tasks/leafa-proof`, and the measurement supports the proposal:
  the rehearsal is easier than the performance in exactly the place that matters.

What this does **not** support is the earlier phrasing that cliques "cannot reach the `0.75`": what
they cannot reach is `2/3` of the unit, and what no certifiable family reaches, at either `n`, is any
of it.

---

## 9.  Reproduce

```sh
cd /home/evand/math/square-packing/s12
S=/path/to/scratch
for t in 3.83375 3.85; do python3 search/n11_exact2float.py search/n11_exact_${t}_support.txt $S/n11_warm_${t}.txt $t; done
CCF="--seed-pitch 0 --stages 10 --inner 400 --N 1000 --threads 2 --dmin 0 --cg-want 300 --exact --time 14400"
python3 search/clique_ceiling.py 191/50   N11Q382    --warm $S/n11_warm_3.83375.txt $CCF
python3 search/clique_ceiling.py 3067/800 N11Q383375 --warm $S/n11_warm_3.83375.txt $CCF
python3 search/clique_ceiling.py 77/20    N11Q385    --warm $S/n11_warm_3.85.txt $CCF
python3 search/clique_ceiling.py 193/50   N11Q386    --warm $S/n11_warm_3.85.txt $CCF
python3 search/clique_ceiling.py 387/100  N11Q387    --warm $S/n11_warm_3.85.txt $CCF
python3 search/clique_ceiling.py 77/20    N11L385    --warm $S/n11_warm_3.85.txt --kmass 4 --r 1 --margin 0 $CCF
python3 search/clique_ceiling.py 191/50   chk --check search/n11_qstab_382_exact.txt      # and 3067/800, 77/20, 193/50, 387/100 for the others
python3 search/clique_ceiling.py 3067/800 chk --check search/n11_qstab_383375_exact.txt
BRANCH_SOLVER=restricted BRANCH_RESTRICTED_PASSES=2 python3 search/branch.py runs/n11_H3985.txt N11AQ3980 \
    --k 4 --r 1 --lam-lo 0 --lam-hi 0 --n 11 --Dp 3980 --N 2000 --colgen 40 --cliques 40 --cq-want 40 --matched --topk 3 --threads 4
verify/target/release/verify runs/branch_N11AQ3980.txt 11 4000 3 0
python3 search/boxclique.py runs/n11_H3985.txt N11BC3980 --Dp 3980 --n 11 --N 2000 --threads 3 --rounds 60 --colgen 30
OMP_NUM_THREADS=1 python3 search/n11_sa2.py search/n11_exact_3.83375_support.txt search/n11_exact_3.85_support.txt
```
