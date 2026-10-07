# LP speed: backends for the cover cutting-plane LP (task C, 2026-08-29)

Code: `search/lp_bench.py` (backends + row sifting), `search/lp_dump_lattice.py` (build an instance
from a leaf's column checkpoint + lattice rows), `search/lp_subset.py` (cut sub-instances and a
`it0 ⊂ it1` round pair), `branch.py --dump-lp PREFIX` (dump the real loop's LP every round).

## Instances

All from the `1110` leaf at `t = 3.98` (asymmetric, per-box multipliers), columns = the 13,043 points
of its round-51 checkpoint (`runs/branch_t398j1110_cols_it51.txt`, the run's actual column set —
`branch.load_cols` would have symmetrised them to 69k, which is what the first attempt at a restart
did and why its round-0 LP did not finish in 55 min), rows = lattice poses (`eta 0.005`, `dt 0.01`
over `[0, 90°)`, 1000 lowest-capture poses per angle under the round-51 weights, merged on the run's
`--row-grid 0.002 0.005`).

| dump | rows | cols | nnz | how |
|---|---|---|---|---|
| `lp1110_pool` | 138,348 | 13,047 | 237.4M | full lattice pool |
| `lp1110_mid_it0` | 40,000 | 12,747 | 66.0M | first 40k pool rows, last 300 point columns dropped |
| `lp1110_mid_it1` | 43,000 | 13,047 | | `it0` + the 3,000 pool rows most violated by `it0`'s optimum + the 300 columns |

~1,650–1,700 nonzeros per row (a unit square captures ~1/8 of 13k points at `t = 3.98`).  The real
round-51 LP was 186k rows × 13k columns, i.e. ~1.35× the pool.

## Scheduling discipline (coordinator correction, 21:25)

Cores `n` and `n+16` are hyperthread siblings.  From 21:25 on: one multi-threaded solve at a time on
`taskset -c 8-15,24-31` with `threads=8`; single-threaded jobs on distinct physical cores with the
sibling idle.  Every timing row records `load` (1/5-min averages), affinity and what else ran.

**Contaminated timings (kept only as objective-agreement checks, re-run alone before entering the
table):**
- `mid_it0` scipy `highs-ipm` + `highspy` IPM, 8 threads on `24-31`, started 21:26 — ran alongside
  PDLP on `9-15` (siblings of `25-31`) and the single-threaded dual simplex of `lp_subset.py` on `8-15`.
- `mid_it0` PDLP `kkt_tolerance 1e-6`, 7 threads on `9-15`, started 21:27 — same neighbours.

## Faithful instances: the real loop, restarted (`lp1110_itN`)

`branch.py … --cols … --cols-raw --max-iters 4 --dump-lp runs/lp1110` restarted the `1110` leaf from
its round-51 probe + checkpoint with the run's own 13,043 columns (a new `--cols-raw`: without it
`build_model`/`load_cols` symmetrise every point of an asymmetric leaf and the restart's LP was
28k–69k columns wide).  Cores `10-15,26-31`, verifier `--threads 12`.

| round | rows | cols | LP (scipy `highs-ipm`) | obj | note |
|---|---|---|---|---|---|
| it0 | 56,983 | 13,043 | 1,812 s | 4.3297 | lattice warm start only, **0 box rows**, λ at cap, probe −1.5 |
| it1 | 58,613 | 13,339 | 1,377 s | 11.4033 | + 1,630 witness rows (5,453 in boxes) + 296 columns |
| it2 | 62,436 | 13,639 | 1,546 s | 11.7025 | + 3,823 rows + 300 columns |
| it3 | 66,053 | 13,912 | 1,658 s | 11.7590 | + 3,617 rows + 273 columns; run stopped (`--max-iters 4`) |

The IPM process sits at 101 % CPU: HiGHS's IPX is effectively **single-threaded** on this LP, so an
IPM solve is a one-core job whatever `threads` says, and the "8 threads" of the earlier contaminated
runs bought nothing.  Load average 8.5–9 throughout (other agents on cores 0–7).

Structure of the optimum (it2): 1,108 columns in the support, **963 rows with positive dual**, but
**13,495 rows within 0.02 of tight** (slack quantiles 10 % = 0.006, 25 % = 0.023, median 0.055).
A flat, massively degenerate LP: the cover is nearly tight on a fifth of the poses, and the
simplex has to walk through that degeneracy whatever basis it starts from.

Loop-hygiene defect seen on the way (not fixed here, outside the task): the lattice warm start keeps
the `perang` lowest-capture poses per angle, which never includes the corner-box poses (they need
`1 + λ`), so round 0 of every restart is a wasted ~30-min solve with `λ` at its cap and probe `−1.5`.

## Results

**Killed (lower bounds only; single-threaded, each alone on its physical core from 22:17):**
- `mid_it0` (40k × 12.7k) HiGHS dual simplex cold: **> 7,900 s**, not finished (killed 23:27).
  ≥ 8× the IPM on the same instance.  Cold dual simplex is not an option at this size.
- `mid_it0` HiGHS PDLP (CPU cuPDLP-C, `kkt_tolerance 1e-6`, crossover off): **> 7,250 s**, not
  finished (killed 23:27).  Not an option either.

### Warm-start sanity checks (coordinator's list, 2026-08-30 01:50)

`search/lp_toytest.py`.  **Toy zero-change test** (same instance, no new rows/columns, carried
basis; synthetic 3,000 × 404 pair): re-solve takes **0 simplex iterations** on the same `Highs`
object, through `run_inc` (addRows/addCols path) and through `run_warm` (`setBasis` path), whether
the pre-solve was IPM + crossover or dual simplex; basis lengths equal the model's.  The wrapper is
not the problem.

**Alignment `it2 → it3`** (exact checks on the dumps): `it2`'s 62,436 row keys are a prefix of
`it3`'s 66,053 (no reordering, dedup or pruning happened — `prune_at 250000`); the old-rows ×
old-columns block, the λ block, `c/lo/hi` on old columns, λ bounds, `b` and the box flags are
identical; new rows are appended `kBasic`, new columns `kLower`, λ statuses carried.  Under `it2`'s
recorded dual `y`: **0** old columns with reduced cost `< −1e−7`, **272 of the 273 new columns**
with negative reduced cost (min −0.39) — as column generation guarantees.  So the mathematically
expected start for `inc` is ≤ 273 dual infeasibilities at objective 11.7025.

What HiGHS reported instead: "Solving LP with useful basis so presolve not used", then dual
phase 1 with `Du: 1128 (sum 24.6)` — the printed objective −24.6 is the phase-1 objective (= minus
that sum), not the LP value.  ~850 dual infeasibilities beyond the new columns, i.e. the basis HiGHS
holds after IPM → crossover → postsolve → "solving the original LP" is not cleanly dual feasible on
the original model.  `inc2` measures this directly (zero-change re-run before any modification) and
then avoids phase 1 altogether: rows first with dual simplex, columns second with primal simplex.

### Row sifting on `it3` (66,053 × 13,912), single thread, core 11 alone

Start = the 13,495 rows of `it2` that were tight (`y > 0` or slack `< 0.02`); each pass appends the
≤ 5,000 most violated rows to the same `Highs` object (dual simplex, basis kept) and rescans all
rows (a sparse mat-vec, 0.07 s).

| pass | active rows | obj | violated after | simplex it | cumulative s |
|---|---|---|---|---|---|
| 1 (cold) | 13,495 | 11.6849360 | 5,775 | 44,287 | 492 |
| 2 | 18,495 | 11.7447577 | 1,838 | 37,867 | 1,054 |
| 3 | 20,333 | 11.7587000 | 247 | 34,769 | 1,664 |
| 4 | 20,580 | 11.7590347 | 47 | 4,990 | 1,756 |
| 5–7 | 20,644 | **11.7590443** | 0 | 946 | **1,794** |

Objective identical to IPM's (`dev 3e−16`), 122,859 iterations in total, 20,644 rows active at the
end (31 %).  **1,797 s vs 1,658 s for IPM: parity, not a win.**  The cost is not the number of
rows but the degenerate walk: passes 1–3 each take 35–45k iterations at ~13 ms whether they add
5,000 rows or 1,800; only once the basis is near-optimal do passes get cheap (4,990 it for 247 rows,
796 for 47).  Sifting does reproduce the full-LP optimum exactly, so it is a valid but not a faster
substitute; in the loop its pass 1 would be warm (previous round's basis) rather than cold, which
saves the 492 s of pass 1 but not the two big passes after it.

Variant `--sift-add 15000` (append every violated row at once): pass 2 appended all 5,775 and took
**149k iterations / 3,063 s** (vs 37.9k / 562 s when capped at 5,000); killed at 3,558 s.  Row
batches must stay small.

### Column sifting (restricted master) on `it3`, single thread, core 11 alone

Where the IPM's time goes: on `it2` (62k × 13.6k) 61 interior iterations at ~22 s each = 1,352 s,
crossover ~100 s.  IPX dualises the model, so each iteration factors the normal matrix on the
**column** side (13.6k, effectively dense: 13.6k³/3 ≈ 8·10¹¹ flops, single-threaded).  Only ~1,100
–1,500 of the 13.9k columns ever carry weight.  `csift`: IPM over an active column set (= `it2`'s
support, 1,108, + the round's 273 new columns = 1,385; all 66k rows), price every column by its
reduced cost `c_j − y·A_j` under the master's dual (0.1 s), add the 1,000 most negative, repeat.
A restricted master is infeasible iff some row has no active column — none here; the code covers
such rows greedily.

| pass | master cols | IPM s (it) | obj | cols with rc < −1e−7 (min) |
|---|---|---|---|---|
| 1 | 1,385 | **51** (54) | 12.2382 | 2,024 (−1.00) |
| 2 | 2,385 | 148 (56) | 11.8187 | 3,273 (−0.38) |
| 3 | 3,385 | 240 (54) | 11.7796 | 2,289 (−0.17) |
| 4 | 4,385 | 377 (55) | 11.7636 | 1,070 (−0.05) |
| 5 | 5,385 | 534 (59) | 11.7594645 | 250 (−0.017) |
| 6 | 5,635 | 578 (57) | 11.7590749 | 47 (−0.009) |
| 7 | 5,682 | 604 (60) | 11.7590461 | 13 (−0.010) |
| 8 | 5,695 | 601 (61) | **11.7590443** | 12 (−3.8e−5) — stopped |

The objective agrees with the full IPM (11.7590443) to 1e−8 after 8 passes / 3,133 s: on a *single*
instance to exact optimality, column sifting is no faster than the full IPM either (the tail of
slightly-negative columns costs a 600-s pass each).  What it measures is the **cost curve of the
master**: 51 s at 1.4k columns, 148 s at 2.4k, 240 s at 3.4k, 377 s at 4.4k, 534 s at 5.4k, 1,658 s
at 13.9k — `t ∝ n^1.5`, ~1 s per interior iteration at 1.4k columns against 27 s at 13.9k.  Every
solve reproduces the same 54–61 IPM iterations; the iteration count does not depend on the width.

### Warm-started dual simplex across a real round (`it2 → it3`: +3,617 rows, +273 columns)

All three variants start from the IPM + crossover optimal basis of `it2` (clean: the zero-change
re-run takes 0 iterations), single thread, one physical core each, alone:

| variant | what HiGHS does | state after ~30 min of simplex |
|---|---|---|
| `inc` (same object, addCols + addRows, dual simplex) | dual phase 1: 1,128 dual infeasibilities (≤ 273 expected from the new columns) | phase 1 sum 24.6 → 0.95 after 12,400 it / 2,840 s; not finished |
| `warm` (new object, `setBasis` with rows by key) | identical: phase 1, 1,547 dual infeasibilities | same trajectory |
| `inc2` stage 1 (rows only, dual simplex) | **`Du: 0` — dual feasible, no phase 1**, as the algebra says | phase 2 crawl: 28,400 it / ~1,800 s, obj 11.7025 → 11.7239 of 11.7590, ~31k primal-infeasible rows along the way; not finished |

So the extra ~850 dual infeasibilities of `inc`/`warm` come from HiGHS's handling of a simultaneous
row + column append (it rebuilds/repairs the basis), not from the carried basis; separating the
two (`inc2`) removes phase 1 exactly as expected — and it does not matter.  The dual simplex on
this LP needs tens of thousands of degenerate pivots per re-solve (13 ms each on 66k rows) however
small the change: 3,617 new rows, of which only the violated ones matter, cost more than a full
IPM.  That is also the explanation of the original observation "the `warm` backend took 200–600 s
where scipy took 0.4 s" (at 31k rows × 388 columns): not a wrapper bug, the same walk at a smaller
size.  (The wrapper *was* also 8× wider than it needed to be on asymmetric leaves — `--cols` without
`--cols-raw` — but that is a separate defect.)

## Comparison table

Single thread throughout (IPX is single-threaded on this LP; the dual simplex and the CPU PDLP
are single-threaded), each job alone on its physical core with the sibling idle, load average
7–12 from other agents on cores 0–7 in every row.  `mid` = `lp1110_mid_it0` (40,000 × 12,747,
66.0M nnz, lattice rows only); `it1`/`it3` = the real loop's LP at rounds 1/3 (58,613 × 13,339,
93M nnz / 66,053 × 13,912).  Reference objectives: `mid` 4.326687287, `it1` 11.403312252,
`it3` 11.759044336.

| instance | method | seconds | objective (rel. dev.) | iterations | note |
|---|---|---|---|---|---|
| mid | scipy `highs-ipm` (baseline) | **1,053** | 4.326687287 | 55 | clean rerun, core 10 |
| mid | scipy `highs-ipm`, 8 threads | 947 | same | 55 | contaminated (shared cores with PDLP) — threads buy nothing |
| mid | `highspy` IPM, 8 threads | 1,011 | same (4e−16) | 56 | contaminated, same neighbours |
| mid | HiGHS PDLP, kkt 1e−4 | 990 | 4.325835 (**−2.0e−4**) | 6,840 | clean; four digits short of what the 2e−6 probe margin needs |
| mid | HiGHS PDLP, kkt 1e−6 | > 7,250 | — | — | killed, unconverged |
| mid | HiGHS dual simplex, cold | > 7,900 | — | — | killed, unconverged |
| it1 | scipy `highs-ipm` (harness) | 1,355 | 11.403312252 | 61 | clean, core 12; the loop itself took 1,377 s |
| it3 | scipy `highs-ipm` (the loop) | **1,658** | 11.759044336 | — | reference for the rows below |
| it3 | row sifting, +5k rows/pass, from `it2`'s tight rows | 1,797 | 11.759044336 (3e−16) | 122,859 | 7 passes, 20,644 rows active |
| it3 | row sifting, +15k rows/pass | > 3,558 | — | > 193k | killed after pass 2 |
| it3 | column sifting, +1k cols/pass, from `it2`'s support | 3,133 (to 1e−8) | 11.7590443 | 8 IPM solves | pass 1 (1,385 cols) **51 s**; pass 3 (3,385) 240 s |
| it2→it3 | `inc` warm dual simplex (crossover basis) | > 5,139 | — | > 32,000 | phase 1 unfinished, killed |
| it2→it3 | `warm` (`setBasis`) | > 5,100 | — | — | same state, killed |
| it2→it3 | `inc2` rows-then-columns | > 1,900 (stage 1) | 11.7239 → … | > 28,400 | no phase 1, still a crawl; process ended without result |

## Recommendation

**Keep IPM; shrink the master.**  Nothing that re-solves the *same* 14k-column LP beats HiGHS's IPM:
the dual simplex (cold, warm, incremental, sifted) is a 40–150k-pivot degenerate walk per solve,
PDLP does not reach certificate accuracy, and IPX's cost per iteration is set by the column count
(`t ∝ n^1.5`: 1 s at 1.4k columns, 27 s at 13.9k) while its iteration count (54–61) is not.  The
loop already prices columns every round (`T.price` / `price_asym`, 300 per round) but leaves every
column it ever generated inside the LP; ~90 % carry no weight.  `BRANCH_SOLVER=restricted`
(`BModel.solve_restricted`) keeps a master of the current support + the round's new columns + what
reduced-cost pricing pulls in (`BRANCH_RESTRICTED_ADD`, default 1,000 per pass; `BRANCH_RESTRICTED_PASSES`,
default 2; `0` = to convergence, i.e. the exact full-LP optimum), trims columns that stay at zero
for three rounds, and prices all columns of the model by reduced cost each round (0.1–0.7 s).  The
master's dual is what the loop's pricing and verifier need; the LP value of a capped master is an
upper bound on the full LP's that converges over rounds — which is how column generation converges
anyway.  With a ~3k-column master a round's LP is 150–400 s instead of 1,400–1,700 s.

Measured in the loop (`runs/restr1110b.out`, same restart as the IPM reference; rows below):

| round | rows | LP value | note |
|---|---|---|---|
| it0 | 58,622 | 4.361 | same wasted round 0 as the IPM reference (no corner rows in the lattice warm start) |
| it1 | 61,804 | 9.252 | |
| it2 | 65,087 | 11.674 | |
| it3 | 67,934 | 11.749 | run stopped (`--max-iters 4`) |

`restr1110c`: four rounds in ≈ 82 min wall clock *including* verifier and pricing, i.e. ≈ 20 min per
round, against ≈ 30–45 min per round for the IPM reference loop above (whose LP alone took
1,377–1,812 s per round).  The values are not directly comparable to the reference's (the witness
rows differ between runs), but they track it (11.749 vs 11.759 at round 3).  The full-size (186k-row)
speed-up is the extrapolation below, not a measurement.

## Large-instance extrapolation

The real round-51 LP was 186k rows × 13k columns and took 6,800 s of IPM.  Rows enter IPX's cost
through the (dualised) normal matrix assembly `A Aᵀ`-side products and the crossover, roughly
linearly (1,053 s at 40k rows / 12.7k cols; 1,355 s at 58.6k / 13.3k; 1,658 s at 66k / 13.9k;
6,800 s at 186k / 13k — the last with a busier machine and 1,600–1,700 nnz per row), so at 186k
rows a 3k-column master should sit at roughly `6,800 × (3/13)^1.5 ≈ 750 s`, ~9× — with the
verifier (100–190 s at `N = 2000`) and pricing unchanged, a round of the `1110` leaf goes from ~2 h
to ~15 min.  Row sifting would cut the rows too, but its passes are dual-simplex walks; if rows
must be cut, do it *by IPM* (solve on the near-tight rows, scan, add) — the one combination not
measured here, because the column side alone already gives the target.

Not measured: multi-threaded dense Cholesky (IPX has none; a solver with a parallel factorization —
none available in this environment — would multiply the gain), and PDLP on a GPU.


