# Witness credit: the LP's rows must carry the verifier's credit, not their own (2026-09-07)

Code: `verify/src/main.rs` (`record_cred`, the witness writer), `search/branch.py`
(`export`/`cqmap`, `read_witnesses`, `add_rows`, `add_clique`, `prune`), `search/tighten.py`
(`read_witnesses`), `search/floatsep.py` (`descend`'s seed unpacking),
`tests/rejection_tests.sh` (§25), `tests/bitid.sh`, `runs/diag_rowvalue.py`, `runs/cutcount.py`.
Read against `search/ANCHOR.md` (the anchor-clique object) and `search/LOOPSPEED.md` §6 (the
previous, unrelated, oracle bug).

**Soundness is untouched.**  Nothing about the verdict path changes: the same cells, the same
exact predicates, the same arithmetic, the same exit codes, the same stdout.  What changes is one
*output* — the witness file the LP reads — and only for a certificate that carries an `anchors`
block.  `tests/bitid.sh` against the pre-change binary is the guardrail, and the certificate
format is not touched at all, so `xcheck.py` and every shipped file are unaffected.

## 1. The defect

`search/branch.py`'s cutting-plane loop feeds the exact verifier's witnesses back as LP rows.  A
witness is one **pose** `(cx, cy, theta_k, sigma_k/2, flag)`, and `BModel.add_rows` used to build
its anchor-clique coefficients with `anchorclique.coeff`, a **pose** predicate.  The verifier
credits a clique to a **cell** — and only when one piece provably holds *every* pose of it
(`certificates/FORMAT.md`, "Anchor cliques"; a cell that straddles a piece's boundary gets
nothing, and `witness_of` then deliberately places the witness *outside* the pieces the cell only
partly meets, but it can only avoid the pieces the sweep recorded as partial, and it is a pose in
a cell that straddles the boundary, so it may perfectly well land inside one).

Two consequences, both measured on the `1110` leaf checkpoint
(`runs/branch_L1110f_probe.txt` in the main tree: `s = 3.98`, `r = 1`,
`lambda = [+1.5, +1.5, +1.5, -2]`, 3531 atoms, 219 anchor cliques, `N = 2000`, `topk 6`):

**(1) Over-credit.**  At the worst witness — `v = 0.969246`, `theta = 73.30 deg`, centre
`(3.3576, 2.4833)`, flag 0 — the LP credited clique 54 (weight `0.2266`) because the pose sits
`0.00325` inside that piece's `meets` hexagon, while the verifier refused it for the cell.  The
LP therefore read the row at `1.195814` and believed it satisfied; the cell was never fixed, and
the loop's four exact rounds it60, it65, it70, it75 read `probe_min = 0.9641, 0.9655, 0.9633,
0.9625` with the objective flat at `11.8695 ± 0.0006` (`runs/branch_L1110f.log`; the checkpoint's
own sweep at `topk 6` is the `0.969246` above).  `finalize()` cannot converge against that: it exports,
the verifier fails, the witnesses come back, the LP says they are already satisfied, and nothing
moves.

Over the whole witness file (`runs/diag_rowvalue.py`, 11,661 witnesses):

| rule | max `LP − verifier` | over-credited (`> 1e-6`) | rows the LP calls satisfied while the verifier fails |
|---|---|---|---|
| `anchorclique.coeff` at the pose (old) | **`+0.226569`** | 2,133 of 11,661 | **2,047** |

**(2) Dropped rows.**  The same 11,661 witnesses collapse to **4,196** dedup keys under the run's
`--row-grid 0.002 0.005` (`add_rows`: `key in self.rkey -> skipped`), 2,461 keys being shared by
more than one witness.  A witness at `0.986` whose key matches a row added rounds earlier at a
*different* pose inside the same grid cell is silently dropped — and for most witnesses the LP row
rule reproduces the verifier's value exactly, i.e. those rows were simply never in the LP.  A grid
of `0.002` in position and `0.005 rad` (`0.29 deg`) in angle is far too coarse to identify a row
on a converged leaf, where the deficits are `10^-2`.

## 2. The fix

**A. The verifier says what it credited.**  Every witness line keeps its five columns
(`value theta_k cx cy flag`) and, **for a certificate with an `anchors` block only**, gains two
more: the number of anchor cliques the sweep credited to the cell the witness came from, and their
file indices (0-based, in the order of the `anchors` block).  The count is always written for such
a file, even when it is `0`, so a reader can tell "credited nothing" from "an old binary".  Box
cliques are not listed — no producer here puts them in a cutting-plane run, and omitting a credit
only makes the LP's row *stricter* than the verifier's, which is the safe direction.

The credited set is deliberately kept **out** of the top-`k` heap.  That heap is selected with
`sort_unstable_by_key` on a `Vec<(i128,f64,f64,u8)>`, and Rust's unstable sort permutes equal keys
differently for a different element type — the trap that bit task K (`search/LOOPSPEED.md` §3,
"Guardrails") and would have changed the witness file of certificates with no anchors at all.  So
the sets live in a side map keyed by the heap entry's own bits, pruned whenever the heap is
(bounded by `4*topk` entries), written only when an `anchors` block is present, and never
allocated otherwise.  On the rare key collision the map keeps the **intersection**, which is
credit the verifier gave either way.

**B. The LP builds the row from that.**  `export` writes `cqmap`, the map from file clique (one
per D4 image) back to the LP's clique **orbit**, next to the probe as
`runs/branch_TAG_probe.txt.cqmap` and on the model as `m._cqmap`.  `read_witnesses` maps the
credited file indices through it and counts, giving `{orbit: number of credited images}` — exactly
the shape `anchorclique.coeff` produced, read off the sweep instead of a pose predicate.
`add_rows` uses those coefficients verbatim for an exact-verifier row.

Float-oracle rows, the lattice warm start, `--corner-rows` and a sibling's dual rows keep
`AC.coeff`: they are heuristic *choices of which row to add*, every one of them is a valid
constraint however it was scored, and the exact rounds are the safety net — the loop only ever
stops on an exact verdict, and an exact round re-emits any pose whose cell is still short with the
verifier's own numbers.  The code says so in a comment at both places.

`add_clique` gives a newly separated clique coefficient `0` on the existing **exact** rows rather
than `AC.coeff`.  Those rows carry what the sweep gave their cell for the cliques that were in the
certificate at the time; a clique separated afterwards was not in that file and gave those cells
nothing, so `0` is the truthful entry, and it keeps the invariant that an exact row's clique
coefficients are the verifier's own.  The new column earns its credit on those poses at the next
exact round, which re-emits them.

**C. Exact witnesses are never deduplicated on the row grid.**  There are at most `topk × bins` of
them a round (12,006 at `topk 6`, `N = 2000`, full range), so nothing is saved by collapsing them,
and §1(2) is what collapsing them costs.  Their key still carries the credited set, and they stay
out of `rkey` altogether, so a float row can never be swallowed by an exact one or vice versa.
The `--row-grid` dedup is unchanged for float rows.

**Every parser of the witness file.**  `branch.read_witnesses` (uses the new columns),
`tighten.read_witnesses` (ignores everything past the first four fields; its model has no clique
columns), `floatsep.descend`'s seed unpacking (was `for (cx, cy, th, _h, _fl) in seeds`, now
positional), `tests/rejection_tests.sh` §16/§17/§19/§25, and `search/boxclique.py` through
`tighten.read_witnesses`.  One more was found and deliberately **not** changed:
`search/nu_f.py` filters its witnesses with `len(q) == 4`, which has matched nothing since the
region-flag column was added long before this work — its `pts` list is always empty.  Relaxing it
would feed the `nu_f` model rows none of its published numbers were measured with, so it carries a
comment saying so instead, to be fixed together with a re-measurement.

**D. Nothing else moves.**  Region flags and the per-box `lambda` are untouched: the flag column is
where it was, `row_active` and `flagmat` see the same values, and a witness whose verifier value is
already `>= 1` after rounding is impossible — the sweep pushes a witness for threshold `j` only
when `sum < W + lambda_j`, which is exactly `v < 1` for the threshold that names it (checked in the
sweep and in the loop's stopping test, which still requires `added == 0` **and** an exact
`mv >= 1`).

## 3. Validation

### Bit-identity and the rejection suite

`sh tests/bitid.sh runs/verify_ref` (the pre-change binary built from `51b93f3`), single thread:
**ALL IDENTICAL** — the 1736-, 764-, 224- and 56-point certificates and the box-clique demo at
`N = 2000`; witness mode (`topk 6`) on the main certificate (3652 witnesses, stdout *and* the
witness file) and on the box-clique demo; a per-box branch trailer in witness mode; the shipped
`k = 1` branch certificate; `TIGHT_DUMP` mode (19,754 cells).  A new case pins the one place the
output is *meant* to change: on an anchor-clique certificate in witness mode stdout is identical
and the witness file is the old one with the credited-clique columns appended — `cut -f1-5`
reproduces it byte for byte.

`./verify.sh`: exit 0, **24 `VERIFIED:` verdicts**, no `REJECTED` and no `NOT VERIFIED` — every
shipped certificate re-verified from scratch, including the main one at `N = 6000` and `N = 12000`
and the anchor-clique demonstration at both, plus `xcheck.py` over every bin of the 56-point set
and 15 JSON round-trips.

`sh tests/rejection_tests.sh`: **138 passed, 0 failed, 0 panics** (was 136).  §25 now
checks the new columns both ways: on a file whose single clique sits where nothing fails every
witness carries an empty credit list, and on one where the clique is credited exactly where the
covering breaks (a point dropped to weight 0 plus its own point clique) most witnesses name it.

### The 1110 leaf: no over-credit anywhere in the file

Same checkpoint, same run as §1 (`N = 2000`, `topk 6`, 8 threads, about 640 s wall), verdict and
minimum unchanged — `min = 969245966245/1000000000000 = 0.969246` at bin
1488, `NOT VERIFIED`, 11,661 witnesses — and `runs/diag_rowvalue.py` over **every** one of them:

| rule | max `LP − verifier` | min | over-credited (`> 1e-6`) | rows the LP calls satisfied while the verifier fails |
|---|---|---|---|---|
| `LP_pose` — `anchorclique.coeff` (old) | `+0.226569` | `-0.010791` | 2,133 / 11,661 | 2,047 |
| `LP_cell` — the verifier's credit (new) | **`+0.000000`** | `-0.001418` | **0** / 11,661 | **0** |

At the worst witness the sweep names **21** credited cliques; the pose predicate credits those 21
**and clique 54**, of weight `0.2266` — nothing else, exactly the diagnosis above.  The clique
credit therefore drops from `0.3917` (pose) to `0.1651` (cell) and the row value from `1.195814`
to **`0.969246`**, which is the verifier's own number to six decimals.  Over
the whole file only 6 witnesses of 11,661 move at all under the new rule, all of them *downward*
(by at most `0.0014`, atoms sitting exactly on a cell boundary): the LP's row is now never weaker
than the cell the verifier failed, and on 11,655 of 11,661 it is exactly it.

The verdict is unchanged, which is the point: `cut -d' ' -f1-5` on the new witness file and the old
one are the same 11,661 lines (differing only in the order of ties, which follows the thread count),
and the `min covered weight ... = 0.969246 (at angle k=1488)` line is identical.

### How much of an exact round reached the LP at all

The sharpest measurement is on the checkpoint itself, and it needs no LP solve
(`runs/cutcount.py`): take the 11,661 witnesses of the sweep above, turn each into the row the
loop would build, and ask whether that row is violated by the very weights the witness was
separated from — i.e. whether it is a **cut** at all.

| path | rows added | of them, cuts | lost |
|---|---|---|---|
| old — pose credit, `--row-grid 0.002 0.005` dedup | 4,196 | **3,307** | 7,465 to the dedup, 889 more that are no cut |
| new — verifier credit, no dedup on exact rows | **11,661** | **11,661** | none |

So an exact round of the stalled leaf was delivering `3,307` usable cuts out of `11,661`
witnesses.  That is the whole of §1 in one line, and it is why the objective could sit at
`11.8695 ± 0.0006` for four exact rounds while the sweep kept failing the same cells.

### Two matched runs of the loop

Three consecutive **exact** rounds from the `L1110f` checkpoint, restarted with `--seed-from` (its
columns, its 800 clique orbits and the 2,883 binding rows of its dual), no column or clique
generation, no float oracle, `topk 3`, `N = 2000`, 8 threads — run once with each credit rule, so
the only difference is what the witnesses become:

| round | `probe_min` (new) | rows added (new) | `probe_min` (old) | rows added (old) |
|---|---|---|---|---|
| `it0` | `-1.4999985` | 6,000 | `-1.4999985` | 2,439 |
| `it1` | `+0.1123757` | 6,000 | `+0.1120711` | 2,339 |
| `it2` | `+0.1622294` | 6,000 | `+0.1625915` | 2,529 |

Both reach `obj = 11.8668276` (the checkpoint's own `11.870423`, to four decimals), and the exact
minimum **rises monotonically over the three rounds** instead of standing still.  The two rules do
not separate here, and that is expected and worth stating plainly: the row set of a `--seed-from`
restart is thin, the deficits in these rounds are `0.85`–`1.5`, and a credit disagreement of at
most `0.227` on 18 % of the witnesses does not decide which constraint binds at that distance from
convergence.  What does separate, in every round, is how much of the round survives: 6,000 rows
against ~2,400.  The regime where the *credit* itself decides the round is the converged one — the
one the leaf was actually stuck in — and there the measurement is the table above: `3,307` cuts out
of `11,661`, with the worst cell of all among the 889 the LP was told it had already fixed.

## 4. Reproduce

`runs/` is gitignored, so the commands are given in full.  `runs/verify_ref` is the pre-change
binary (build `verify/` at `51b93f3` into a separate `CARGO_TARGET_DIR`); `$M` is the main tree's
`runs/`, which holds the `L1110f` checkpoint read-only.

```sh
cd verify && cargo build --release && cd ..
sh tests/bitid.sh runs/verify_ref
sh tests/rejection_tests.sh
M=/home/evand/math/square-packing/s12/runs
# the exact sweep on the leaf checkpoint (about 10 min at 8 threads)
./verify/target/release/verify $M/branch_L1110f_probe.txt 12 2000 8 6 runs/L1110f_sep_new.txt
python3 runs/diag_rowvalue.py $M/branch_L1110f_probe.txt runs/L1110f_sep_new.txt 2000
python3 runs/diag_rowvalue.py $M/branch_L1110f_probe.txt $M/diag_L1110f_sep.txt 2000   # the old file
python3 runs/cutcount.py $M/branch_L1110f_probe.txt runs/L1110f_sep_new.txt 2000
# three consecutive exact rounds from the checkpoint, with the new credit rule and with the old
# one (BRANCH_WITNESS_CREDIT=pose); ~20 min each.
export BRANCH_SOLVER=restricted BRANCH_CQ_MAX=800 OMP_NUM_THREADS=8
for mode in cell pose; do
  BRANCH_WITNESS_CREDIT=$mode python3 search/branch.py $M/branch_L1110f_probe.txt L1110$mode \
      --k 1110 --r 1 --cols-raw --no-warm --seed-from $M/branch_L1110f \
      --prune-at 250000 --row-grid 0.002 0.005 --N 2000 --sep exact --threads 8 --max-iters 3
done
```

`BRANCH_WITNESS_CREDIT=pose` restores the pre-2026-09-07 rule (`anchorclique.coeff` at the pose,
dedup as for any other row).  It exists for that A/B and nothing else; `cell` is the default and
the only rule the loop should be run with.

## 5. What this does and does not settle

*Fixed and measured.*  The LP's exact rows now carry the verifier's own clique credit — nowhere in
an 11,661-witness file does a row read higher than the cell it came from — and an exact round
delivers every witness as a cut instead of a third of them.  The mechanism that made `finalize()`
unable to converge on a clique leaf is gone: a cell the verifier fails now always produces a row
the LP has to move for, so the loop's `added == 0 and mv >= 1` stopping test means what it says.

*Not settled by this.*  Whether either leaf **closes** is a different question.  `finalize()`
terminating means the covering can be made to hold at the leaf's multipliers; the leaf closes only
if `total - sum_j lambda_j k_j < 12` at that point, and `search/ANCHOR.md` puts the `k = 4` leaf at
`12.000` with the anchor family worth `10^-4` on the cover side, while the `1110` leaf's last
converged LP value was `11.8704` with the covering still `0.03` short.  Nothing here moves those
numbers; it removes the reason the loop could not get to them.
