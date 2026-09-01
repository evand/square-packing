# Loop speed: the exact verifier stops being the separation oracle (task K)

Code: `search/floatsep.py` (new), `search/branch.py` and `search/tighten.py` (`--sep`,
`--exact-every`, `--sep-pitch`, `--sep-dtheta`, plus `--stab`, `--seed-from`, `--corner-rows` in
`branch.py`), `verify/src/main.rs` (the sweep, and the `VERIFY_STATS` counters).
**Every default reproduces the old behaviour** — nothing changes unless one of the new flags is
given, so the reproduce commands in `README.md` and the run notes still mean what they say.

**Soundness is untouched.**  The float oracle only decides *which rows the LP gets*.  Every row it
emits is a valid constraint by construction (`floatsep.py`, "Rows are valid": for any centre in a
bin's admissible box there is an angle in the bin at which the unit square is inside the container,
and that square contains the concentric `sigma_k`-square the row is written for), whatever the scan
computed.  The cutting-plane loop only ever *stops* on an exact verdict, and the file that is
finally claimed still goes through the unchanged Rust verifier and `xcheck.py`.  The sweep changes
in `verify/src/main.rs` are guarded by bit-identity against the pre-change binary.

Machine for every timing below: 32 logical / 16 physical cores, pinned with `taskset` to 8
distinct physical cores (`8-15`) unless stated; the load average at the time is quoted with each
table.  `runs/verify_ref` is the pre-change binary, built from the source at `de40c27`.

## 1. Instrumentation, and where the time went

`VERIFY_STATS=1 verify ...` prints one line to **stderr** after the verdict (nothing at all without
the variable, so stdout, the witness file and `TIGHT_DUMP` are untouched):

```
VERIFY_STATS atoms= bins= threads= wall= | strips= cells= active= (avg /strip)
  | anchor: pieces= calls= tests= skip= ok= core= | cliqueboxes: cands= tests= | witnesses= | ns/cell=
```

`search/branch.py`'s per-round bracket gains an `fsep` phase and, in float mode, a
`[sep float_min=… fviol=… clean=… exact]` group; `tighten.py` marks a float round's `probe_min`
with `(float)`.

### Baseline (`runs/verify_ref`, N = 2000, 8 threads, load 20–31, `runs/baseline_ref.txt`)

| instance | atoms | anchor pieces | bins swept | topk | wall |
|---|---|---|---|---|---|
| `branch_t398hk4_probe` | 180 | 0 | 830 (`[0,45°]`) | 3 | 0.64 s |
| `branch_J16_probe` | 144 | 48 | 2001 | 3 | 1.39 s |
| `branch_t398ik4n_probe` | 1084 | 176 | 2001 | 3 | 134.3 s |
| `branch_J6_probe` | 2268 | 240 | 2001 | 6 | 259.1 s |
| `branch_J12_probe` | 1308 | 400 | 2001 | 6 | 545.1 s |
| `branch_t398j1110_probe_it51` | 3285 | 0 | 2001 | 3 | 546.2 s |

The counters say what that time is (`runs/ab_out.txt`, 4 threads):

| run | cells | active/strip | anchor tests | of which succeed | witnesses | ns/cell |
|---|---|---|---|---|---|---|
| `t398ik4n` (176 pieces), N = 400 | 470 M | 351 | 2.94 G | 95.4 % | 41 M | 172 |
| `J12` (400 pieces), N = 300 | 467 M | 417 | 11.57 G | 94.8 % | 24 M | 646 |
| `t398ik4n` with the `anchors` block stripped, N = 2000 | 977 M | 352 | — | — | 294 M | 125 |

Three things, in order of size:

1. **The anchor-clique credit.**  6 (176 pieces) to 25 (400 pieces) `(cell, piece)` tests per cell
   at ~25 ns each — about 70 % of a clique run.  And **95 % of those tests succeed**: the sweep
   was not searching for the pieces that hold a cell, it was re-proving the same ones cell after
   cell, through a `Vec` indirection per test, with a linear "already credited?" scan.
2. **The witness.**  The Sutherland–Hodgman clip that places the LP witness ran for *every*
   violated cell — a third of the arrangement in an unconverged probe — and made six heap
   allocations per call.  On a certificate without cliques that was ~3/4 of the whole sweep.
3. **Two sorts per strip** (the active atoms by `u1`, then their `2A` breakpoints) plus two binary
   searches per cell — a `log A` factor on top of an arrangement that is already `Θ(m·A)`.

## 2. The float separation oracle (`search/floatsep.py`)

Per angle bin of the verifier's own `N`-net:

* the bin's rational `cos`/`sin`, its `sigma_k` (rounded down to 1e-6) and its admissible centre
  box `[L, s-L]^2` with `L = w_min(k)/2` come from the verifier's own integer formulas
  (`bin_geometry`), so a row this module emits has exactly the shape of a verifier witness and
  lands in the same `--row-grid` bucket;
* the captured atom weight is evaluated on a regular grid of centres in the bin's rotated frame
  (`--sep-pitch`, default 0.004) by a prefix-sum sliding window (`lp_search.scan_angle`'s scan);
* every anchor-clique **image** is credited its weight on the grid poses it provably holds.  In the
  bin's frame `contains(A)` is an axis-parallel rectangle and `meets(F)` is a rectangle intersected
  with the segment's own normal slab, so a clique costs one rasterised rectangle plus a couple of
  half-plane tests, not a scan of the grid;
* for a branch certificate the region threshold `1 + lambda_j` is subtracted where the centre lies
  in corner box `j`, and the row carries that box's flag;
* the `topk` worst cells per bin are taken one per `block × block` tile, so the rows are spread out
  the way the verifier's per-bin `topk` spreads them.

**Conservative in the verifier's direction.**  Clique membership uses the `sigma_k`-square for
`contains` (the verifier credits by the *larger* exact bin core) and the hexagon `A + [-h,h]^2`
including the segment's own normal for `meets`, with `TOL = 1e-9` shaved off every inequality in
the safe direction — the identical predicates `search/anchorclique.py` uses to build the LP column.
So the float oracle never claims coverage the exact sweep would deny.  It can of course miss a
violation because the grid is finite, which is why the loop still runs the exact sweep periodically
and never stops without it.

Cross-check (`runs/check_floatsep.py`): the grid credit against `anchorclique.member`, pose for
pose, on `certificates/s12_anchorclique_demo_3.9318.txt` and `runs/inst/branch_J12_probe.txt`
(200 image cliques), 7 bins spread across the net — **1.50 M poses, 744 k of them credited,
0 mismatches, max |diff| = 0**.

### Cost, against the exact sweep it replaces

8 threads, `--sep-pitch 0.004 --sep-dtheta 0.4 --topk 6`, load 9–16
(`runs/fsep_out.txt`, `runs/pair_out.txt`):

| instance | atoms / pieces | exact sweep, old | exact sweep, new | **float scan** | float vs old |
|---|---|---|---|---|---|
| `branch_t398hk4_probe` | 180 / 0 | 0.63 s | 0.26 s | 0.6 s | 1.0× |
| `branch_J16_probe` | 144 / 48 | 1.45 s | 0.93 s | 2.3 s | 0.6× |
| `branch_t398ik4n_probe` | 1084 / 176 | 140.4 s | 73.8 s | 5.2 s | **27×** |
| `branch_J6_probe` | 2268 / 240 | 236.4 s | 85.8 s | 6.6 s | **36×** |
| `branch_J12_probe` | 1308 / 400 | 552.9 s | 339.0 s | 11.0 s | **50×** |
| `branch_t398j1110_probe_it51` | 3285 / 0 | 548.9 s | 165.3 s | 1.4 s | **392×** |
| `s12_t3.98_corner_k1` | 11892 / 0 | — | — | 1.5 s | — |
| `s12_t3.98_corner_k0` | 9984 / 0 | — | — | 1.2 s | — |

On the two tiny leaves the exact sweep is already sub-second and the float scan is slower; the
float oracle is for the instances where the loop actually spends its life.

How close is the float minimum to the exact one, on the same file?

| instance | exact min | float min |
|---|---|---|
| `branch_J12_probe` | 0.891139 | 0.891139 |
| `branch_J16_probe` | 0.416667 | 0.416667 |
| `branch_J6_probe` | 0.974619 | 0.983162 |
| `branch_t398ik4n_probe` | 0.523590 | 0.582004 |
| `branch_t398hk4_probe` | 0.891305 | 0.913044 |
| `branch_t398j1110_probe_it51` | 0.939586 | 0.981502 |

The float minimum is an over-estimate wherever the grid misses the worst pose, which is exactly why
`probe_min >= 1` on the float oracle is treated as a hint (run the exact sweep now) and never as a
verdict.

### In the loop

`branch.py --sep float` (default `exact`):

* every round the float oracle runs and its worst poses become the rows;
* the exact verifier runs when `it % --exact-every == 0` (default 5, so round 0 is exact) **and**
  as soon as the float oracle has reported `probe_min >= 1` in two consecutive rounds;
* the loop's stopping test insists on an exact verdict (`mv is not None and mv >= 1`), so a clean
  float round can never end the loop; if the exact check fails, its witnesses become rows and the
  loop continues exactly as before;
* the probe checkpoint is still written every round, so a restart loses nothing.

`tighten.py` has the same two flags with the same defaults.

**Amortised separation cost per round** at `--exact-every 5` — four float rounds and one exact
sweep (new binary) in five, against today's exact sweep every round:

| leaf | today | `--sep float`, per round | speed-up |
|---|---|---|---|
| `branch_t398ik4n_probe` | 140.4 s | (4×5.2 + 73.8)/5 = 19.0 s | 7.4× |
| `branch_J6_probe` | 236.4 s | (4×6.6 + 85.8)/5 = 22.4 s | 10.6× |
| `branch_J12_probe` | 552.9 s | (4×11.0 + 339.0)/5 = 76.6 s | 7.2× |
| `branch_t398j1110_probe_it51` | 548.9 s | (4×1.4 + 165.3)/5 = 34.1 s | 16.1× |

At `--exact-every 10` the same four leaves give 13×, 20×, 12.6× and 30×.  In a *float* round alone
the separation is 27–392× cheaper (table above).  The live rounds these came from were 80–99 %
verifier (J16 round 9: 3.8 h of which 3.0 h verifier; J12 round 9: 4.8 h of which 4.7 h), so on
those the round-level speed-up is the separation speed-up: **7–16× at the default, 12–30× at
`--exact-every 10`**, after which the LP is the round again (task C's territory: J16's LP alone was
40 min a round, J12's 55 s).

### Validation: the same LP value, and the same verified file

* **Same LP value.**  The `k = 4` leaf of `branch_t398hk4_*` run five ways from the same
  checkpoint (`runs/rounds.sh`, `runs/rounds_out.txt`): exact oracle, exact + `--corner-rows`,
  float + `--corner-rows`, float + `--corner-rows --stab 0.3`, and exact + `--seed-from`.  All five
  converge to **obj = 12.139560** — identical to six decimals, against a tolerance of ±0.002 — and
  each stops on an exact `probe_min = 1.0000010`.
* **Same verified file.**  `runs/e2e.sh` runs `tighten.py reopt` on
  `certificates/s12_lower_3.931795_sparse.txt` both ways to completion.  Exact oracle: 6 rounds,
  224 points, total **11.9859020**.  Float oracle (`--sep float --exact-every 5`): 8 rounds, 224
  points, total **11.9859028** (a difference of 8e-7).  Both files are then re-read by the
  unchanged verifier — `VERIFIED … s(12) >= 3.931795386`, min `400001/400000` at `N` and at `2N` —
  and by `xcheck.py`, the independent `Fraction` implementation, which agrees on the minimum.

## 3. The exact sweep (`verify/src/main.rs`)

Three changes, all bit-identical by construction:

1. **The anchor credit stops re-deriving what the strip already fixed.**  The filter lists are
   flattened into one array per bin (they were a `Vec` per piece — a pointer chase per test);
   the `u0` half of each `meets` half-plane and the piece's `u0` window verdicts are constants of
   the strip and are computed once per strip per piece (~30 pieces) instead of once per
   `(cell, piece)`; and "has this clique already been credited in this cell" is a stamp array
   instead of a linear scan of the cliques credited so far.
2. **The witness clip runs in reusable scratch buffers.**  `clip_rect_poly_buf` is the same
   Sutherland–Hodgman, same arithmetic, same order — without the six `Vec` allocations per call.
3. **No sorting inside a strip.**  Each bin builds its atom list in `u1` order once (`gy`); a
   strip's active list is a linear filter of the `[ylo - h, yhi + h]` band of that list (an atom
   outside it is in no window of any cell of the strip *and* contributes no breakpoint inside
   `[ylo, yhi]`), the `u1` breakpoints are the two sorted sequences `yv ± h` merged with
   `(ylo, yhi)` in one pass, and the two `partition_point`s per cell are two sweeping pointers,
   since the cell bounds only increase.

The spatial restriction the brief asks for is in 1 and 3: the strip already keeps only the atoms
whose squares can meet it, 3 restricts them again to the `u1` band the strip's admissible centres
can reach, and 1 does the same for the anchor pieces by hoisting everything the strip decides out
of the cell loop.

### Measured, ref and new back to back on the same instance (`runs/pair.sh`, `runs/pair_out.txt`)

8 threads, cores 8-15, load 14–22 throughout (each pair ran within a minute of the other, so load
drift affects both equally).

| instance | atoms | anchor pieces | N | topk | ref | new | speed-up |
|---|---|---|---|---|---|---|---|
| `branch_t398hk4_probe` | 180 | 0 | 2000 | 3 | 0.63 s | 0.26 s | **2.4×** |
| `branch_J16_probe` | 144 | 48 | 2000 | 3 | 1.45 s | 0.93 s | **1.6×** |
| `branch_t398ik4n_probe` | 1084 | 176 | 2000 | 3 | 140.4 s | 73.8 s | **1.9×** |
| `branch_J6_probe` | 2268 | 240 | 2000 | 6 | 236.4 s | 85.8 s | **2.8×** |
| `branch_J12_probe` | 1308 | 400 | 2000 | 6 | 552.9 s | 339.0 s | **1.6×** |
| `branch_t398j1110_probe_it51` | 3285 | 0 | 2000 | 3 | 548.9 s | 165.3 s | **3.3×** |
| `s12_t3.98_corner_k0` | 9984 | 0 | 500 | 3 | 422.1 s | 88.9 s | **4.7×** |
| `s12_t3.98_corner_k1` | 11892 | 0 | 500 | 3 | 664.5 s | 140.6 s | **4.7×** |
| `s12_lower_3.9686` (shipped) | 1736 | 0 | 6000 | 0 | 89.1 s | 22.3 s | **4.0×** |

Read down the atom-only rows: 2.4× at 180 atoms, 3.3× at 3285, 4.7× at 10–12 k.  The win grows
with the atom count because what was removed (the per-strip sorts, the per-cell binary searches)
is exactly the `log A` factor.  The anchor-clique rows gain less — 1.6× on the 400-piece J12 —
because the piece tests are still the dominant term there; see "What is still slow".

### Guardrails

* `tests/bitid.sh` against `runs/verify_ref`: **ALL IDENTICAL** (stdout, witness files and
  `TIGHT_DUMP`), and a wider set in `runs/quickid.sh` — the anchor-clique demo, the box-clique
  demo, a per-box branch trailer, and the five real leaves (`t398ik4n`, `J12`, `J6`, `j1110`,
  `corner_k1`) in witness mode — **QUICK ALL IDENTICAL** at every stage of the change.
* `tests/rejection_tests.sh`: **136 passed, 0 failed, 0 panics.**
* `verify.sh`: **exit 0**, 24 `VERIFIED:` verdicts, no `REJECTED` and no `NOT VERIFIED`
  (`runs/verify_sh_full.txt`) — every shipped certificate re-verified from scratch, including the
  main one at `N = 6000` and `N = 12000`, plus the rejection suite and the JSON round-trips.

One thing bit-identity does *not* cover, and did bite once: the per-bin top-`k` selection is
`sort_unstable_by_key` on a `Vec` whose element type is `(i128, f64, f64, u8)`, and Rust's unstable
sort permutes equal keys differently for a different element type.  An earlier version of change 2
deferred the witness by putting a bigger struct on that heap, and the witness file changed —
same values, different cells among the ties.  The witness is therefore still computed eagerly; only
its allocations are gone.

## 4. Fewer rounds (`branch.py`)

* `--corner-rows` seeds the warm start with poses whose centre lies in an **active** corner box, on
  a lattice, whatever they capture.  `lattice_rows` keeps only the poses that capture least under
  the input weights, and on the `k = 4` leaves that is never a corner pose, so round 0 had *no* row
  charged to any multiplier: `lambda` went straight to its cap and the first LP value was
  meaningless.  For a symmetric model the D4 images of a pose have the same coefficient vector, so
  only box 0 over `[0, 45°]` is emitted (1322 rows on the `t398hk4` leaf).
* `--seed-from LEAF` starts from a sibling's checkpoint: its columns, its clique columns and the
  **binding rows of its dual** (`runs/branch_LEAF_dual.txt`, written every round, one line per row
  with `y > 0`).  The per-slot leaves differ only in the multipliers, so the sibling's binding row
  set is almost the right one; rows charged to a box this leaf has switched off are dropped by
  `row_active` on the way in.
* `--stab F` is in–out stabilisation: in the float rounds the oracle separates at
  `(1-F)*x + F*(previous separation point)` instead of at the LP vertex.  The exact rounds always
  separate at the true LP point, so the stopping rule is unchanged.

Measured on the `k = 4` leaf of `branch_t398hk4_*` (`runs/rounds.sh`, `runs/rounds_out.txt`), all
five reaching **obj = 12.139560**:

| run | rounds to converge | `probe_min` at round 0 | LP value at round 0 |
|---|---|---|---|
| S0 exact, as today | 7 | −0.1666638 | 10.000032 |
| S1 `--corner-rows` | **6** | +0.2750010 | 12.000024 |
| S2 `--corner-rows --sep float` | 11 | +0.2750010 | 12.000024 |
| S3 S2 `--stab 0.3` | 18 | +0.2750010 | 12.000024 |
| S4 `--corner-rows --seed-from S1` | **4** | +0.4011656 | **12.139560** (already final) |

* `--corner-rows` is a clear small win and fixes the round-0 pathology: `lambda` is paid for from
  the start, so the first LP value is the real one instead of `10.0` with `lambda` free at its cap.
* `--seed-from` is the big one here: 7 rounds → 4, and the LP is at its final value in round 0.
  The remaining rounds are spent re-proving the covering, not moving the value.
* The float oracle needs *more* rounds (11 vs 6) because it emits far fewer rows per round
  (99 against 1044 on this leaf: 213 scanned bins × topk against 2001 exact bins × topk).  On a
  leaf whose verifier costs 0 s that is a loss; on J12, where a round costs 553 s of verifier, 11
  cheap rounds beat 6 expensive ones by a wide margin.  A finer `--sep-dtheta` buys rows back at
  a proportional cost in scan time.
* **`--stab 0.3` did not help** on this leaf: 18 rounds against 11.  Damping the separation point
  also damps the rate at which genuinely violated poses are found.  It is left in, off by default,
  and should not be turned on without measuring.

## 5. What is still slow

* **The anchor-piece tests.**  On the 400-piece J12 leaf the sweep is still 25 `(cell, piece)`
  tests per cell and 339 s a sweep.  95 % of those tests *succeed*, and within a strip a piece's
  rectangle test is monotone in the cell index (both cell bounds only increase), so the set of
  cells a piece credits is an interval and the whole credit could be an event sweep — `O(pieces)`
  a strip instead of `O(cells × pieces)`, worth an estimated further 3–6× on clique leaves.  It is
  not done here because the cells a piece *partly* meets feed the witness, and getting that set
  wrong changes the witness file, i.e. breaks the bit-identity guarantee.  Doing it safely needs
  the witness's `partial`/`apart` sets recomputed on demand for the ~5 % of cells that produce a
  witness.
* **The LP.**  Once the separation is off the critical path the LP is the round again.  On the
  `t398ik4n` checkpoint *with* its 10 174-column checkpoint loaded, a single round is 117 s → 239 s
  → 409 s of `highs-ipm` and the separation is invisible; that is task C's problem, not this one.
* **The float oracle's row yield.**  At `--sep-dtheta 0.4` it looks at 213 of 2001 bins, so it
  emits ~1/10 the rows of an exact round.  The right setting is probably per-leaf and should be
  chosen by watching `fviol` and the round count in the log.
* **`--stab`** needs a better rule (a real in–out with the stability centre updated only on
  progress, or a proximal term in the LP) before it earns its place.

## Reproducing

```sh
sh runs/pair.sh                     # ref vs new exact sweep, back to back per instance
sh runs/fsep_bench.sh 0.004 0.4     # the float oracle on the same instances
sh runs/rounds.sh                   # rounds-to-converge: corner rows / float / stab / seed-from
sh runs/e2e.sh                      # tighten.py reopt both ways, then verify + xcheck both files
python3 runs/check_floatsep.py certificates/s12_anchorclique_demo_3.9318.txt
sh tests/bitid.sh runs/verify_ref   # the guardrail
sh tests/rejection_tests.sh
```
