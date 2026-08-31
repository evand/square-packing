# The m^2-4 family and rung 2 (task H; 2026-08-30)

Code: `search/zeromargin.py` (exact checker, vectorised, task H step 1), `search/closed4.py`
(heuristic cover LP, container side `s` a free parameter), `search/rung2_close.py` (new: a
stress-driven separation loop trying to close the gap between a converged LP cover and the exact
checker's zero-margin requirement), `search/nu_f.py` (rigorous, float-checked lower bound on
`nu_f(s)`; fully parametrised by the container side read from a certificate, or passed directly
to `run`). Every number below is labelled **certified** (exact `Fraction` arithmetic, or
float-checked with an explicit `~1e-9`-scale margin the way `CLOSED4.md` §7 already treats as
rigorous), **heuristic** (a sampled LP or a local search value, not a bound), or a **bound**
(direction stated).

## 1. Table: `m -> cover LP value, certified cover, certified packing mass, n-1`

| `m` | `n = m^2-4` | cover LP value (heuristic) | certified closed cover | certified packing mass (rigorous, float-checked lower bound on `nu_f`) | `n-1` |
|---|---|---|---|---|---|
| 4 | 12 | 12.4175 total / **12.509 honest** (`runs/closed4_best.txt`, `CLOSED4.md`); closing loop plateaus at **12.388 total / ~12.46 honest** (§2); scaled to `x1.02` / `x1.03` (**12.666** / **12.790**, §2b) removes the known exact violation but leaves other boxes uncertified at depth 8-17 (positive margin, slow convergence) | **no** (§2, §2b) | **10.68** (`CLOSED4.md` §7) | 11 |
| 5 | 21 | **20.75-20.9**, still oscillating/rising when stopped (§3.1) | not attempted (bracket only, per the brief) | **16.05** (§3.2, this task) | 20 |

(`m = 6, 7` not started — gated on `m = 4, 5` being finished or the time budget being spent,
per the brief; neither happened, see §4.)

## 2. Rung 2, `m = 4`, `W < 13` — **not achieved**

**Step 1 (oracle) — done.** `search/zeromargin.py`'s float pre-filter is vectorised with `numpy`
over the whole point set (`ZEROMARGIN.md` §7). Rung 1 is byte-identical before/after (`6,958`
boxes / `3356/70/74/3179/0` reduced domain, `27,396` / `13556/298/280/12364/0` full domain,
`zeromargin_stress.py` `0` failures). A `--depth 6` exact check of the 1,972-point
`runs/closed4_best.txt` cover now takes **74 s on 8 processes** (was 575 s on 4) for identical
leaf counts (`111,020` boxes; `24905` CORE / `350` P1 / `0` TRI / `4459` EMPTY / **`28996`
UNCERTIFIED**). New CLI modes: `--oracle FILE` (candidate cutting-plane poses from uncertified
boxes) and `pose CERT --cx --cy --u` (exact captured weight at one rational pose, no subdivision).

**The baseline cover is not valid, exactly, not just by a float estimate.** `CLOSED4.md` reported
a stress-test estimate of `0.99267` captured weight at `(3.3965, 1.4235, 76.41°)`. Using the new
exact `pose` mode with a nearby rational angle (`u = 7873211803/10000000000`, `θ ≈ 76.4282°`,
admissible):

```
python3 search/zeromargin.py pose runs/closed4_best.txt --cx 3.3965 --cy 1.4235 --u 0.7873211803
EXACT captured weight = 1240843/1250000 = 0.9926744  ***VIOLATION: < 1***
```

This is a **certified** counterexample (`Fraction` arithmetic, 157 points checked exactly): the
2026-08-26 `closed4_best.txt` cover of `[0,4]^2` genuinely fails at this closed pose, not merely
"unconverged in the checker's box subdivision." Rung 2 needed either a repaired cover or a proof
this specific one certifies; it does neither.

**Step 2 (separation loop) — plateaus below 1, does not close.** `search/rung2_close.py` ran the
separation loop the brief asks for over the fixed 1,972-point column set: solve the LP, stress-
test at pitch `0.003` (`search/closed4.py`'s own dense scan) to find the actual worst poses, add
those (+ a local polish) as hard rows, resolve, repeat — each round targets exactly where the
previous round's row set was blind, a genuine (float) cutting-plane loop, not another generic
refinement. Run: `taskset -c 0-4`, `nproc=5`, 67 rounds, 5,210 s (~87 min); rows grew
`14,046 -> 239,723`.

| rounds | LP value | stress min (pitch 0.003 + polish) | honest cost |
|---|---|---|---|
| 0 (initial, thr 1.02 lattice only) | 12.252 | 0.799 | 15.33 |
| 1-9 | 12.25 - 12.37 | 0.81 - 0.98 (recovering from the fresh row set) | - |
| 10-30 | 12.36 - 12.38 | 0.87 - 0.99 (whack-a-mole across wall/tilted poses) | - |
| 47-67 (plateau) | **12.3868 - 12.3883** (flat) | **0.9928 - 0.9958** (flat band, best seen 0.99580) | **~12.46** |

(Full 67-round table: `runs/closed4_rung2c.log` / `.json`.) The plateau is stable for 20+
consecutive rounds with no further improvement despite continued row generation (each round adds
2,000-4,500 new hard cutting-plane rows and the LP value moves by `< 0.001`): **this is strong
evidence, not proof, that the fixed 1,972-point column set has an intrinsic residual violation of
about half a percent that no reweighting can remove** — closing it would need new columns (points)
placed at the persistent violation loci (the worst poses recur near `cx ≈ 0.5-0.7` or `3.3-3.5`,
`cy ≈ 1.3-2.6`, `θ ≈ 20-90°`, i.e. the same wall/tilted-wall region `CLOSED4.md` already flagged),
which is a further column-generation step this run did not take. No certificate was exported from
this run: `rung2_close.py` (like `closed4.py run`) only writes `_best.txt` after its loop exits
normally, and the process was stopped at the plateau (continuing showed no further gain within the
observed 20+ rounds) rather than run to its full 3-hour budget for a marginal, uncertain gain — the
in-memory LP weights at the plateau were therefore not preserved. The shipped, checkable artefact
for rung 2 remains `runs/closed4_best.txt`, and it is exactly and provably **not** a valid rung-2
cover (see the counterexample above).

*No tilted tight family of the `ZEROMARGIN.md` §4 item 4 kind was identified* — the residual
violations are ordinary unconverged tilted-wall poses (small but nonzero margin, not exactly-1
tight families), consistent with `CLOSED4.md`'s own diagnosis. The two-region primitive
`ZEROMARGIN.md` flags as possibly-needed was therefore not needed for what was found, but it may
still be needed if a genuinely tight family appears once the column set is fixed by adding points.

**Verdict: rung 2 (`W < 13`, exactly zero-margin certified at `m = 4`) is open.** `s(13) = 4` is
not re-proved by this work. What is established: (a) the existing best heuristic cover is provably
invalid at an exact rational pose; (b) 67 rounds of a genuine stress-driven cutting-plane loop over
its column set push the honest cost from `12.509` to about `12.46` but plateau roughly `0.5%`
short of a valid cover; (c) the next step, not taken here for lack of time, is column generation
(new points) targeted at the plateau's residual violation region, or accepting a small increase in
compute for the same loop with priced columns enabled (`rung2_close.py --allow-colgen` was
stubbed as a CLI flag but the column-adding logic itself was not implemented in this session).

## 2b. Coordinator follow-up (2026-08-30): scale the cover instead of reweighting it

**Idea.** Captured weight scales linearly with the point weights, so if a cover's true minimum
captured weight is `mu > 0` everywhere, multiplying every weight by any `lambda >= 1/mu` gives a
valid closed cover of total `lambda * (original total)`. Since the exact worst pose found for
`runs/closed4_best.txt` is `0.9926744` (`mu >= 0.9804` needed for `lambda = 1.02`), scaling by
`1.02` or `1.03` should fix ordinary small-margin poses everywhere at once (P1 is unaffected by
scaling: `|c-p|_inf <= 1-w/2` doesn't involve the weight at all), leaving only genuinely tight
(margin-zero) or actually-worse-than-found poses as open questions.

**Scaled certificates.** `search/scale_cover.py IN.txt P Q OUT.txt` multiplies every integer
weight numerator by `P` and the weight denominator by `Q` -- exact, no floats, no rounding:

| file | scale | total |
|---|---|---|
| `runs/closed4_best_x102.txt` | `x 51/50 = 1.02` | `6332916636/500000000 = 12.665833` |
| `runs/closed4_best_x103.txt` | `x 103/100 = 1.03` | `12790008108/1000000000 = 12.790008` |

Both `< 13`.

**Exact checker, increasing depth, `--nproc 8`, `taskset -c 0-7`:**

| cert | depth | boxes | CORE | P1 | TRI | EMPTY | UNCERTIFIED | wall time |
|---|---|---|---|---|---|---|---|---|
| x1.02 | 8  | 214,956   | 63,449  | 338 | 0 | 4,481 | 42,410 | 133 s |
| x1.02 | 10 | 370,588   | 128,351 | 363 | 0 | 4,692 | 55,088 | 275 s |
| x1.02 | 12 | 579,610   | 255,238 | 395 | 0 | 5,090 | 32,282 | 535 s |
| x1.02 | 14 | 742,658   | 317,471 | 412 | 0 | 5,561 | 51,085 | 622 s |
| x1.02 | 17 | 1,092,112 | 446,398 | 441 | 0 | 6,250 | 96,167 | 840 s |
| x1.03 | 8  | 189,814   | 61,004  | 314 | 0 | 4,149 | 32,640 | 115 s |
| x1.03 | 10 | 300,294   | 114,471 | 337 | 0 | 4,327 | 34,212 | 207 s |
| x1.03 | 12 | 427,784   | 191,623 | 357 | 0 | 4,590 | 20,522 | 377 s |
| x1.03 | 14 | 534,220   | 227,550 | 372 | 0 | 4,894 | 37,494 | 745 s |

**Neither scaling reaches 0 uncertified, and the count does not shrink monotonically** (it dips
then rises again for both: `42k -> 55k -> 32k -> 51k -> 96k` for x1.02, `33k -> 34k -> 21k -> 37k`
for x1.03) -- this is the "boxes stop shrinking" stopping condition, not a crash or a runaway
count, so per the brief's instructions this is reported honestly rather than pushed further on
faith. `--tri` was not used (`C(1972,3)` is intractable; no candidate triangle was identified
either, see below).

**Diagnosis (step 3): case (b), not (a) or (c).** The uncertified boxes cluster at the same
location for both scalings: a wall-touching square with its left edge on the container wall
`x = 0` and its right edge on the interior grid line `x = 1`, centred near `cy approx 1.5`
(and its three D4-type images by the `x -> 4-x`, `y -> 4-y` symmetry), angle `theta -> 0`.

*Not case (a).* The exact `pose` mode confirms strictly positive margin there, for both scalings:

```
python3 search/zeromargin.py pose runs/closed4_best_x102.txt --cx 0.5015625 --cy 1.4969 --u 0.000872665
EXACT captured weight = 266364789/250000000 = 1.065459156  (OK: >= 1)
python3 search/zeromargin.py pose runs/closed4_best_x103.txt --cx 0.5015625 --cy 1.4969 --u 0.000872665
EXACT captured weight = 537952417/500000000 = 1.075904834  (OK: >= 1)
```
(and a second point at the reported uncertified box's corner, x1.02: `106655127/100000000 =
1.06655127`, also `>= 1`). So the true minimum here is **not** below `1/1.02` or `1/1.03` -- the
scaling did what it was supposed to at this pose.

*Not case (c).* The margin is a healthy `6.5-7.6%`, nowhere near zero, so this is not an
unhandled zero-margin tight family (no triangle or two-region primitive is implicated).

*Case (b): why CORE/P1 converge slowly here, specifically.* Direct inspection of `cert_core` /
`cert_p1` on the reported box (clipped to its admissible slice, `[0.5,0.5] x [1.49688,1.5] x
[0,~0.0004]deg`): **P1 sums to only `0.533`** of the needed `1` (`search/zeromargin.py`'s P1 wall
shortcut only admits points with `p_x` within `[0,1]` for a box this close to the LEFT wall --
points near the mirror wall `p_x approx 3` fail the *other*, non-wall-shortcut inequality because
they are geometrically far from this box, not because the near-wall shortcut doesn't apply to
them); **the exact CORE sum is only `0.306`** over the same clipped rectangle (CORE must hold for
*every* pose in the box simultaneously, a strictly harder requirement than the `1.065` captured at
any single pose in it). Root cause: the natural witness for a wall-square whose far edge sits on
`x = 1` would be a cover point at exactly `(1, cy)` for the box's centre `cy`, but this cover
(`CLOSED4.md`) deliberately splits its weight on the `x = 1` line to points at `y = 1.45` and
`y = 1.55` rather than placing one at `y = 1.5` (the same degenerate-line "cusp" structure
`ZEROMARGIN.md` §4 item 3 describes for *interior* squares) -- so there is no single coincident
witness, only two off-centre ones, each of which needs the box to shrink further before it alone
carries enough weight. Both primitives provably converge to the true value (`>1.06`) as the box
shrinks to a point (CORE and P1 are monotone in box size), but the convergence here is unusually
slow because the admissible-width margin `w(theta)/2 - 1/2` grows *linearly* in `theta` near a
wall-touching pose (not quadratically, as in the pure corner case `ZEROMARGIN.md` §4 item 1 where
P1 closes in one box) -- so precision in `cx` requires proportionally finer `theta`, and the
subdivider splits whichever of `(dx, dy, 2 d(theta))` is currently longest, which does not
preferentially resolve the actual bottleneck. This is an **engineering gap in the current
primitives**, not a mathematical obstruction: candidates for closing it are (i) a wall-square
lemma for *off-centre* witnesses (generalising `ZEROMARGIN.md` item 2 to two points straddling the
symmetric line, analogous to DS7 Lemma 2's `(1,y)`-or-`(1+x,y)` disjunction), (ii) biasing the
subdivision to split `theta` preferentially near `theta = 0` boxes, or (iii) enough further depth
-- we did not locate where (if anywhere within reach) the oscillation turns into clean shrinkage.

**Verdict on this experiment: rung 2 still not achieved.** No certificate is shipped (the
coordinator's step 4 was conditional on certification). The scaling idea is validated in the sense
that it removes the previously-identified *exact* violation (the `0.9926744` pose now reads
`>1.06` at both scalings) and narrows the open question to a specific, well-characterised, and
apparently fixable convergence weakness in `cert_p1`/`cert_core` for wall-adjacent squares with
off-centre witnesses, rather than a structural deficiency of the point set itself.

**Process improvement (coordinator note, implemented).** `search/rung2_close.py` and
`search/closed4.py run` now export the current LP weights to `runs/<tag>_last.txt` **every
round** (`C.export` / `export`), not only at loop exit -- a plateau's weights are no longer lost
if the process is stopped or killed. (`search/nu_f.py run` already did this per-iteration via
`write_cert_scaled`, no change needed there.)

## 3. `m = 5`, `n = 21`, `W < 21` — bracketed, gap not closed

### 3.1 Cover LP (heuristic)

`python3 search/closed4.py run --s 5 --tag s5 --nproc 3 --col-pitch 0.05 --no-literature --time 10800`
(the `s = 4`-specific Bentz/Nagamochi literature columns do not generalise — their `d4` mirror
images are computed at a hard-coded container side `4.0` inside `closed4.py`, so they were
excluded rather than silently mis-mirrored at `s = 5`; the grid-lattice and dual-pricing columns
are fully general in `s` already and need no code change). Stopped after 10 iterations (~29 min on
3 cores) once the trend was clear and cores were reallocated to the certified bound (§3.2):

| iteration | LP value | rows | support | lattice min | dual coverage max | t (s) |
|---|---|---|---|---|---|---|
| 0 | 16.000 | 18300 | 5 | 0.000 | - | 11 |
| 1 | 19.909 | 26862 | 138 | 0.537 | - | 54 |
| 2 | 20.829 | 32841 | 144 | 0.901 | 1.328 | 120 |
| 3 | 20.908 | 34531 | 201 | 0.935 | 1.144 | 196 |
| 4 | 20.892 | 35567 | 245 | 0.979 | 1.345 | 323 |
| 5 | 20.824 | 37022 | 275 | 0.972 | 1.101 | 424 |
| 6 | 20.781 | 38611 | 330 | 0.964 | 1.094 | 646 |
| 7 | 20.785 | 39685 | 355 | 0.986 | 1.064 | 888 |
| 8 | 20.772 | 40877 | 426 | 0.947 | 1.063 | 1323 |
| 9 | 20.753 | 41756 | 497 | 0.984 | 1.073 | 1705 |

(Full table: `runs/closed4_s5.log` / `.json`.) Dual coverage was still above `1.06` (pricing keeps
finding columns worth adding) and `viol` still in the hundreds when stopped, i.e. **not
converged** — as at `m = 4`, this value is a relaxation over its own sampled rows/columns and can
only move (mostly rise, occasionally dip when new columns relax it) with more compute. It was
already oscillating in a `20.75-20.9` band, i.e. **below `n = 21`** but not decisively so and with
no exactness behind it at all (sampled rows, `1e-9` float containment). A second, independent
cover-LP run inside `search/nu_f.py run 5 ...` (built from scratch with its own smaller/coarser
column set, see §3.2) reached a comparable but less refined `21.5` after 28 iterations —
consistent with, not contradicting, the `closed4.py`-specific run's `20.75-20.9`.

### 3.2 Certified packing-side lower bound

`search/nu_f.py`'s `ExactLP`, `max_coverage`, and `lower_from_cert` are already fully
parametrised by the container side (read from a certificate file, or passed as a `Fraction`
directly to `run`) — no `m = 4`-specific code was found, so **no generalisation work was needed**,
only running it:

```
python3 search/nu_f.py run 5 10800 s5nuf 0.01 0.01
```

(`taskset -c 5-7`, `nproc=3`; the optional exact-verifier cross-check step fails harmlessly every
round — `verify/target/release/verify` is not built in this worktree — the `try/except` around it
in `lower_from_cert`/`run` means this only forgoes a few extra cutting-plane rows, not correctness.)
28 iterations, ~63 min, stopped at a plateau (`L` oscillating `13.0-16.1` for the last 10+
iterations with no new best):

| it | t (s) | cover LP (this run's own, from scratch) | `M` (max certified coverage) | `L = mass/M` |
|---|---|---|---|---|
| 0 | 1 | 22.00 | 2.500 | 8.800 |
| 6 | 83 | 25.00 | 1.900 | 13.158 |
| 11 | 531 | 24.32 | 1.583 | **15.371** |
| 19 | 1295 | 22.71 | 1.415 | **16.053** (best) |
| 24 | 2384 | 21.50 | 1.370 | 15.693 |
| 27 (last) | 3769 | 21.53 | 1.403 | 15.347 |

**Rigorous (float-checked, `1e-9`-scale margins, same standard as `nu_f(4) >= 10.68`):
`nu_f(5) >= 16.053`.**

### 3.3 Reading the bracket

`16.053 <= nu_f(5) <= COVER(5) <=` (heuristic, unconverged) `~20.8`, against `n - 1 = 20` and
`n = 21`. Neither side closes: the certified lower bound (16.05) is well short of 20, and the
heuristic cover (20.8) is close to but not decisively below 21, and has no certification behind it
at all. **The fractional gap looks structurally similar to `m = 4`** (there, `nu_f(4) >= 10.68`
against `n - 1 = 11`, cover heuristic `12.4-12.5` against `n = 12` — the certified side reaches
roughly `89-97%` of `n-1` at both `m = 4` and `m = 5` with the compute spent here, the heuristic
cover sits within `1-4` units above `n`), so this task's data point is consistent with, but does
not prove, task D's hypothesis that the deficit is universal across the `m^2-4` family rather than
a small-`m` artefact. A firm answer needs both sides converged (more rounds/columns), which the
compute budget here did not reach.

**No certified cover, no attempted exact zero-margin check, at `m = 5`** — out of scope for this
session per the brief ("bracket... only `m = 6` if everything else is done").

## 4. What is certified, heuristic, or a bound, explicitly

| statement | status |
|---|---|
| Rung 1 (Friedman's 14 points certify `s(15) = 4`) | **certified**, exact `Fraction`, re-verified byte-identical after vectorisation |
| `zeromargin.py`'s vectorised pre-filter is sound (superset of the exact test) | **certified** by construction (same lenient `1e-9` inequalities, now over `numpy` arrays) and by the byte-identical regression |
| `runs/closed4_best.txt` fails at `(3.3965, 1.4235, θ ≈ 76.428°)`, captured weight `1240843/1250000` | **certified**, exact `Fraction`, 157 points checked |
| `runs/closed4_best.txt` costs `12.4175` (LP) / `12.509` (honest, stress-checked) | **heuristic** |
| The `rung2_close.py` plateau (`12.388` total, stress min `~0.994-0.996`) | **heuristic** (a float local/lattice search, not a bound; no exact check was run on its point set, which was not exported) |
| `nu_f(4) >= 10.68` | **bound**, rigorous, float-checked `1e-9` margins (`CLOSED4.md` §7; not re-run here) |
| Rung 2 (`W < 13`, exactly zero-margin certified at `m = 4`) | **not established**; the specific counterexample above shows the obvious candidate fails |
| `m = 5` heuristic cover LP value | **heuristic**, unconverged, `20.75-20.9` (§3.1) at the point this session stopped it |
| `nu_f(5) >= 16.053` | **bound**, rigorous, float-checked `1e-9` margins, same standard as `nu_f(4)` (§3.2, new this session) |

## 5. Files

| file | what |
|---|---|
| `search/zeromargin.py` | exact checker (vectorised), `friedman14` / `cert` / `pose` modes, `--oracle` |
| `search/rung2_close.py` | rung-2 separation loop (stress -> hard rows -> resolve) over a fixed column set; now checkpoints `runs/<tag>_last.txt` every round |
| `search/scale_cover.py` | exact rational weight scaling of a certificate (`IN.txt P Q OUT.txt` scales by `P/Q`) |
| `runs/closed4_rung2c.log` / `.json` | the `m = 4` closing-loop run (67 rounds; no `_best.txt`, see §2) |
| `runs/closed4_rung2c_attempt1.log` | the first (discarded) attempt, killed at round 6 for regressing below the baseline before being judged prematurely -- kept for the record |
| `runs/closed4_s5.log` / `.json` | the `m = 5` heuristic cover LP run (10 iterations; no `_best.txt`, stopped deliberately, §3.1) |
| `runs/nuf_s5nuf.log` (`= s5nuf.nohup.log`) / `.json` / `_cur.txt` | the `m = 5` certified packing-mass run (§3.2); `_cur.txt` is the last iteration's own (unconverged, unverified) cover-format export |
| `runs/closed4_best.txt` | the `m = 4` reference cover (2026-08-26, `CLOSED4.md`); **exact counterexample above shows it is not a valid rung-2 cover** |
| `runs/closed4_best_x102.txt` / `_x103.txt` | the `x1.02` / `x1.03` exactly-scaled covers (§2b); totals `12.665833` / `12.790008`, both `< 13`; **not certified** (uncertified boxes remain at depth 8-17, diagnosed as slow-converging positive-margin boxes, not violations) |

To reproduce the exact counterexample: `python3 search/zeromargin.py pose runs/closed4_best.txt
--cx 3.3965 --cy 1.4235 --u 0.7873211803`. To reproduce the rung-1 certification:
`python3 search/zeromargin.py friedman14 --tri --depth 14 --nproc 4`. To reproduce the scaled-cover
positive-margin check: `python3 search/zeromargin.py pose runs/closed4_best_x102.txt --cx 0.5015625
--cy 1.4969 --u 0.000872665`. To reproduce a depth-check on a scaled cover: `python3
search/zeromargin.py cert runs/closed4_best_x102.txt --depth 14 --nproc 8`.
