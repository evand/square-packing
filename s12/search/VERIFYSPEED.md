# Verifier speed on anchor-clique certificates (task: verify-speed, 2026-09-08)

Code: `verify/src/main.rs` (`anchors_d4_closed`, the event-driven anchor credit in `min_cover_k`,
the scheduler in `main`), `tests/rejection_tests.sh` (§29–31, 34 new checks), `tests/bitid.sh`,
`certificates/FORMAT.md` ("Anchor cliques", the symmetry paragraph), `VERIFICATION.md`.
Machine: 32 logical / 16 physical cores, load 4–5 throughout (one unrelated 4-process job);
`runs/verify_ref` is the pre-change binary built from `008ef89`, `runs/verify_new` the new one.
Instances: the `k = 4` anchor-clique probe `runs/branch_J16i_probe.txt` (7844 atoms, 3092 of them
zero-weight anchor atoms, 1208 anchor cliques of two pieces over 2416 anchors, corner trailer
`lambda 1.5 / k 4`, `t = 3.98`) and the `1110` probe `runs/branch_L1110f_probe.txt` (3531 atoms,
219 cliques, per-box lambdas `[1.5, 1.5, 1.5, −2]`, not symmetric), both at `N = 2000`, the net
the cutting-plane loop uses.

**Every change keeps the verdict logic and the reported numbers identical.**  The three levers
were validated as follows, all on the final binary: `tests/bitid.sh runs/verify_ref` is ALL
IDENTICAL (single thread; plain, witness, `TIGHT_DUMP`, branch, box-clique and anchor-clique
modes); every verdict of `verify.sh` is unchanged; `tests/rejection_tests.sh` is 172 checks, 0
failures, 0 panics (was 138); on the two probes the per-bin stdout *and* the witness file are
byte-identical to the old binary on every bin compared (§3), and the incremental credit agrees
with the retained per-cell test at every cell of those bins (`VERIFY_ANCHOR_XCHECK=1`).

## 1. Profile: where the 20–120 s per bin went

`VERIFY_STATS=1 VERIFY_BINS=k:k`, one thread, old binary (`runs/base_*.out`):

| instance, bin | cells | (cell, piece) anchor tests | succeed | corner tests | wall | ns/cell |
|---|---|---|---|---|---|---|
| `k = 4`, bin 0 | 12.1 M | 1.375 G (233 per credited cell) | 97.7 % | 13.2 M (1 %) | 16.2 s | 1341 |
| `k = 4`, bin 400 | 63.8 M | 8.58 G | 98.1 % | 30.3 M | 101.5 s | 1591 |
| `k = 4`, bin 829 | 65.5 M | 9.59 G | 97.1 % | 42.7 M | 115.6 s | 1765 |
| `k = 4`, bin 1000 | 60.5 M | 8.32 G | 97.4 % | 26.6 M | 100.8 s | 1665 |
| `k = 4`, bin 1500 | 65.1 M | 8.74 G | 98.5 % | 24.7 M | 103.6 s | 1593 |
| `k = 4`, bin 1999 | 69.2 M | 10.32 G | 98.5 % | 89.3 M | 118.0 s | 1706 |
| `k = 4`, bin 0, **anchor block stripped** | 12.1 M | — | — | — | 0.90 s | 74 |
| `1110`, bin 0 | 6.1 M | 49.9 M (35 per cell) | 94 % | 0.36 M | 0.62 s | 102 |
| `1110`, bin 1000 | 13.4 M | 210 M | 95 % | 1.0 M | 2.38 s | 178 |
| `1110`, bin 1488 | 13.7 M | 213 M | 95 % | 0.9 M | 2.48 s | 181 |

(The stripped run's 74 ns/cell includes a witness placement for half of its cells, since
without the clique weight they all fail; the pure sweep is ~50 ns/cell.)  So on the `k = 4`
probe **94 % of the time is the anchor credit**, and the anchor credit is not geometry: only
0.3 % of the (cell, piece) tests reach the exact bin-core corner test, 97–98 % of them succeed,
and what costs 1.5 µs per cell is walking ~230 pieces per cell through the window checks, the
filter half-planes and the "already credited" stamp — the same pieces at every cell of a strip,
with the same answer.  That is the shape the review's "incremental credit" estimate assumed, and
the estimate (~5×) was, if anything, low: the whole anchor cost is removable, and the remainder
is the sweep itself.  The full sweep is 2000 bins (`k = 0..1999`; the banner prints `k=0..2000`) at ~100 s ≈ 55 CPU-hours; the loop's
log (`branch_J16i.log`, round 55) records 64,362 s of exact verification on its 3 cores.

Levers, in order of the measured payoff: **2** (incremental credit, 5–19× per bin), **1** (the
D4 reduction, 2.41× on the bin count, 2.18× of wall), **3** (scheduling, 2 % here).  The brief's order
(1, 2, 3) is by the review's estimates; the profile put 2 first, and it was done first.

## 2. The three levers

### 2.1 Incremental anchor credit (`min_cover_k`)

**What the per-cell test decides.**  In the sweep of one bin, a strip is a u0-interval `[a, b]`
of centres (over DEN) and its cells are the intervals `[by[t], by[t+1]]` of the strip's
u1-arrangement.  For a piece `pi` the old code decided, per cell, with the cell rounded outward
to `[a0s, b0s] × [c0s, c1s]` on the 1e-9 grid (`c0s = floor(c0·ASC/DEN)`, `c1s = ceil(c1·ASC/DEN)`):

    credited(pi, cell)  ⟺  (W1) the cell is in the piece's u1-act window (awind over DEN)
                        ∧  (W3) sg_u0[pi]  ∧  c0s ≥ w1l  ∧  c1s ≤ w1h
                        ∧  (C)  fast path  ∨  every corner (v0, v1) of (endpoint − cell) is in the bin core
                        ∧  (F)  for every filter f:  G_f(c0s − fp1h_f)  ∧  G_f(c1s − fp1l_f)

with `G_f(r) := flo_f + n1_f·r ≥ −bnd_f ∧ fhi_f + n1_f·r ≤ bnd_f` (the segment-normal
half-planes; `flo/fhi` are constants of the strip), the fast path being
`sg_u0cm[pi] ∧ c0s ≥ w1l + m ∧ c1s ≤ w1h − m` with `m = core_margin`.  (The old
`lo = flo + min(n1·r_l, n1·r_h) ≥ −bnd ∧ hi = fhi + max(…) ≤ bnd` is exactly `G(r_l) ∧ G(r_h)`.)

**Every one of these is an interval condition on `c0s` and on `c1s` separately**, given the
strip: (W1) and (W3) are windows, each `G_f` is a linear inequality in one variable, so
`{r : G_f(r)}` is an integer interval `[rlo, rhi]` (two exact divisions per filter per strip),
and the fast path is a window.  Since `c0s` and `c1s` are monotone in the cell index `t`, each
condition holds on a contiguous range of cells, found by a search in `by` (the conversions
`c0s ≥ A ⟺ c0 ≥ ⌈A·DEN/ASC⌉`, `c0s ≤ B ⟺ c0 ≤ ⌈(B+1)·DEN/ASC⌉ − 1`,
`c1s ≥ A ⟺ c1 ≥ ⌊(A−1)·DEN/ASC⌋ + 1`, `c1s ≤ B ⟺ c1 ≤ ⌊B·DEN/ASC⌋` are exact; the
strip-independent ones are computed once per bin, `athr`).  So per (strip, piece) the code now
computes three cell ranges — the window run (W1 ∧ ¬skip, needed for the witness's `part` list),
the candidate run (W1 ∧ W3 ∧ F) and the sure run (candidate ∧ fast path) — and runs the exact
corner test only on the candidate cells outside the sure run: the two bands of width `m` at the
ends of the window, one or two cells each.  The credited cells of the piece are emitted as
start/end **events** on cell indices; the cell loop keeps a counter per clique and a running
credit, so a cell costs O(events at that cell) instead of O(pieces).  `done` (credited cliques)
and `part` (pieces the cell meets but that do not hold it) for a witness cell are read off the
active sets; they are the same sets as before (`record_cred` sorts, `witness_of` uses `all`).

**Exactness, not just soundness.**  The decision per (cell, piece) is the old one by construction
— the conditions are the same, only evaluated per range — and the old per-cell closure is kept
in the binary: with `VERIFY_ANCHOR_XCHECK=1` it is re-run at every cell and the credit, the
`done` set and the `part` set are compared; any difference is exit 3 with no verdict.  Rejection
tests §30 run six anchor certificates (point and segment anchors, filters, straddling cells, a
per-box trailer, witness mode) under that flag; the `k = 4` bin 0 (12.1 M cells, 5.9 M credited)
and the `1110` bins 0/1000/1488 were run under it as well, all clean, and their stdout and witness
files are byte-identical to the old binary's.  Under-crediting would be sound but would change the
reported minima; there is none.

Three further exact reductions, each measured (bin 400 of `k = 4`, one thread): the per-cell
rounding `c0s, c1s` by an f64 estimate corrected with exact multiplications instead of a 128-bit
division (13.5 → 12.3 s); the searches in `by` by exponential search from the previous piece's
positions (pieces are visited in u1 order); and skipping a corner's `in_core` call when the
fast-path lemma already proves it inside the core — `|v0| ≤ h2 − m` (the strip is `m` inside the
u0 window) and `|v1| ≤ h2 − m1(v0)`, `m1 = ⌈(|v0| sin δ − h2(1 − cos δ))/cos δ⌉⁺ + ⌈h2(1 − cos δ)⌉ + 2`
(12.3 → 8.05 s).  The lemma: with `|v0| ≤ h2 − m`, `R_0Q` is immediate; `R_δQ`'s first axis needs
`|v0| cos δ + |v1| sin δ ≤ h2`, which `m ≥ h2 (sin δ + 1 − cos δ)` gives for every `|v1| ≤ h2`
(`(1 − c)(1 + c − s) ≥ 0`); its second axis needs `|v0| sin δ + |v1| cos δ ≤ h2`, i.e. the first
term of `m1`; the sector condition applies within `δ` of an axis, and there `|v|² ≤ v0²/cos²δ`
(near the v0-axis, covered by `m ≥ h2(1 − cos δ)`) or `|v|² ≤ v1²/cos²δ` (near the v1-axis, the
second term of `m1`).  A skipped call is one whose answer is `true`; the decision is unchanged.

**Per-bin result** (`runs/new3_*.out`, one thread, output byte-identical in every case):

| instance, bin | old | new | speed-up | new ns/cell |
|---|---|---|---|---|
| `k = 4`, bin 0 | 16.23 s | 2.84 s | 5.7× | 234 |
| `k = 4`, bin 400 | 101.5 s | 8.05 s | 12.6× | 126 |
| `k = 4`, bin 1000 | 100.8 s | 7.66 s | 13.2× | 127 |
| `1110`, bin 0 | 0.62 s | 0.35 s | 1.8× | 57 |
| `1110`, bin 1000 | 2.38 s | 0.59 s | 4.0× | 44 |
| `1110`, bin 1488 | 2.48 s | 0.59 s | 4.2× | 43 |

The `1110` bins are at 43–57 ns/cell, i.e. at the cost of the point sweep itself; the `k = 4`
bins keep ~75 ns/cell of per-(strip, piece) work (10 M evaluations per bin, ~5 band cells each).

### 2.2 The `[0°, 45°]` reduction with an anchor block (`anchors_d4_closed`)

**Soundness.**  Let `g` be a symmetry of the container `[0, s]²`.  It maps closed unit squares to
closed unit squares, preserves inclusion and intersection, so the image of the piece
`{S : A_a ⊆ S, S ∩ A_f ≠ ∅ ∀f}` is the piece with anchors `gA_a`, `gA_f`, and the image of a
clique (a union of pieces) is the clique with every anchor mapped.  If for each generator `g`
of D4 the multiset `{(w_K, gK)}` equals `{(w_K, K)}` — the family is *D4-closed* — then
`Σ_{K ∋ gS} w_K = Σ_{K ∋ S} w_K` for every pose `S`; with D4-symmetric atoms (`check_symmetry`,
unchanged) and equal region thresholds (the corner boxes are D4-symmetric) the covered weight
minus the threshold is D4-invariant.  Every pose at an angle in `[45°, 90°)` is the diagonal
image of one at an angle in `(0°, 45°]`, so the bins `k = 0..829` of the `N = 2000` net, which
cover `[0°, 45.03°]`, bound every pose.  Nothing else changes: the per-bin lower bounds are the
same conservative ones.  The check is exact: anchors are canonicalised as unordered pairs of
reduced rational endpoints (the image of `x/d` under `x ↦ s − x` is `(s_num d − x s_den)/(s_den d)`,
reduced), pieces as (anchor, sorted filter anchors), cliques as sorted piece lists, and the sorted
list of (weight, clique) must be identical after each generator.  A closed family written
differently (a duplicated piece, say) fails the check and gets the full sweep, which is always
sound; box cliques still get it always (a box is tied to a bin of the net, and the image of a bin
is not a bin).  `VERIFY_FULLSWEEP=1` forces the full range.

**The trap it must not fall into** is rejection test §29: test 8's symmetric family, which fails
only at 39.8° and is repaired by the D4 orbit of one point, with the eight repair points turned
into anchor point cliques.  All eight → VERIFIED over `[0°, 45°]`, bins `0..829`.  Seven of them,
missing the one mirror image → the family is not closed, the sweep is full, and the certificate is
REJECTED with every violation at 50.2°: a verifier that reduced because the *atoms* are symmetric
would have said VERIFIED.  Eight anchors with one weight zeroed → not closed, REJECTED.  A
two-piece `K(p, A)` at the centre with a vertical `A` → not closed under the diagonal (`[0°, 90°)`);
with its horizontal image → closed (`[0°, 45°]`); with a box-clique block as well → full.

**Where it applies.**  `branch.py`'s single-lambda leaves keep D4-symmetric atom sets and add
anchor cliques as whole orbits with one weight (`anchorclique.images`), so the `k = 4` probe is
closed (1208 distinct cliques, invariant under all three generators; checked independently in
Python as well): 2000 bins → 830, **2.41×** on the bin count.  The witness file then covers
`[0°, 45°]` only, which is what the loop's row set wants: a symmetric model has D4 orbits as
columns, so a pose and its image are the same row, and `ANCHOR.md` §4 records that the doubled
sweep was what made the row pruning cycle.  The `1110` leaf is a per-box pattern with raw point
columns (`branch.py`: "every point is its own column"), so its atoms are not even swap-invariant
and no reduction is available to it — the lever is worth exactly nothing there, as the profile
says, and its gain comes from lever 2 alone.

### 2.3 Dynamic bin scheduling (`main`)

Bins are handed out one at a time from an atomic counter instead of contiguous blocks of
`nbins/threads`; bin costs vary 7× across the net (16 s at bin 0 against 100–118 s elsewhere on
`k = 4`, 0.6–2.5 s on `1110`), so a static block can leave threads idle for a large fraction of
the last block.  Nothing semantic depends on the schedule, and two things that *used* to depend on
it no longer do: the minimum is taken lexicographically in (value, bin), and every witness carries
(bin, index) so that the file is the single-threaded order — bins in increasing order, the stable
sort by value on top — whatever the thread count.  Rejection test §31 checks the witness file
(11,635 lines) and the minimum line are identical at 1 and 4 threads and under `VERIFY_STATIC=1`
(the old split, kept for measurement).

## 3. Whole-probe validation and speed-ups

All runs by `runs/val.sh` (detached; `runs/verify_ref` = old binary, `runs/verify_new` = new),
3 threads for the whole sweeps, one thread per bin for the per-bin pairs; load 4–5 throughout.

### 3.1 Per-bin identity and speed on the `k = 4` probe (`runs/samp_*`, one thread each)

Every 100th bin, old binary against new binary with the full sweep forced.  For all 20 bins
(`k = 2000` is outside the net) the stdout minus the `D4-symmetric` and `VERIFY_STATS` lines
and the witness file are identical:

```sh
for k in $(seq 0 100 1900); do
  diff <(grep -v -E 'VERIFY_STATS|^D4' runs/samp_old_$k.out) <(grep -v -E 'VERIFY_STATS|^D4' runs/samp_new_$k.out) \
  && cmp runs/samp_old_$k.sep runs/samp_new_$k.sep; done      # 20 of 20 identical
```

| bin | old | new | speed-up | minimum (both) | stdout + witness file |
|---|---|---|---|---|---|
| 0 | 15.8 s | 2.92 s | 5.4× | 0.999747 (999747269846/10^12) | identical |
| 100 | 110.3 s | 9.19 s | 12.0× | 0.999550 | identical |
| 200 | 100.9 s | 8.16 s | 12.4× | 0.999894 | identical |
| 300 | 103.7 s | 8.16 s | 12.7× | 0.999019 | identical |
| 400 | 101.3 s | 7.84 s | 12.9× | 1.000000 (1000000070482/10^12) | identical |
| 500 | 101.1 s | 8.00 s | 12.6× | 1.000000 (1000000366643/10^12) | identical |
| 600 | 104.8 s | 7.98 s | 13.1× | 0.999875 | identical |
| 700 | 109.1 s | 8.02 s | 13.6× | 1.000000 (1000000003469/10^12) | identical |
| 800 | 112.5 s | 7.71 s | 14.6× | 0.999199 | identical |
| 900 | 111.5 s | 7.84 s | 14.2× | 1.000001 (1000000742139/10^12) | identical |
| 1000 | 99.4 s | 7.21 s | 13.8× | 1.000001 (1000000656614/10^12) | identical |
| 1100 | 102.1 s | 7.33 s | 13.9× | 1.000000 (1000000014967/10^12) | identical |
| 1200 | 99.8 s | 7.29 s | 13.7× | 1.000000 (1000000366643/10^12) | identical |
| 1300 | 102.3 s | 7.16 s | 14.3× | 0.997514 | identical |
| 1400 | 101.8 s | 7.20 s | 14.1× | 1.000000 (1000000151139/10^12) | identical |
| 1500 | 99.8 s | 7.04 s | 14.2× | 0.998809 | identical |
| 1600 | 97.2 s | 5.62 s | 17.3× | 1.000000 (1000000233342/10^12) | identical |
| 1700 | 100.9 s | 7.05 s | 14.3× | 0.999686 | identical |
| 1800 | 107.8 s | 7.02 s | 15.4× | 0.999812 | identical |
| 1900 | 116.9 s | 6.15 s | 19.0× | 1.000001 (1000000749920/10^12) | identical |
| **sum, 20 bins** | **1999 s** | **144.9 s** | **13.8×** | | |

The old binary's full sweep was not repeated (55 CPU-hours); its cost is the sample's mean of
99.95 s per bin over the 2000 bins `k = 0..1999`: **199,900 CPU-s = 55.5 CPU-h, i.e. ~66,600 s
wall at 3 threads** — which is what the loop's log shows for the same probe (`branch_J16i.log`,
round 55: `ver 64362s` on 3 cores, the round whose exact minimum `0.9970929` the runs below
reproduce).

### 3.2 Whole-probe runs, 3 threads (`runs/val_*`)

| run | bins | wall | minimum | witness file | verdict |
|---|---|---|---|---|---|
| `1110`, old binary, `topk 6` | 2000 | **1699.5 s** | 0.969246 (969245966245/10^12) at k=1488 | 11,661 lines | NOT VERIFIED |
| `1110`, new binary, `topk 6` | 2000 | **416.5 s** | 0.969246 (969245966245/10^12) at k=1488 | 11,661 lines | NOT VERIFIED |
| `k = 4`, old binary, full sweep (extrapolated, §3.1; log: 64,362 s) | 2000 | **~66,600 s** | 0.9970929 (log) | | NOT VERIFIED |
| `k = 4`, new binary, full sweep forced (`VERIFY_FULLSWEEP=1`), `topk 3` | 2000 | **5158.3 s** | 0.997093 (997092898038/10^12) at k=411 | 2,131 lines | NOT VERIFIED |
| `k = 4`, new binary, D4 reduction, `topk 3` | 830 | **2367.2 s** | 0.997093 (997092898038/10^12) at k=411 | 1,021 lines | NOT VERIFIED |

*Decisions and minima.*  `1110`: stdout identical line for line (`diff` of the two `.out` files
minus the `VERIFY_STATS` lines is empty) and the witness files are identical as multisets
(`sort` both, then `cmp`) — the old 3-thread file's order among equal values was
schedule-dependent, the new one is the sequential order; the sweep placed 144,659 witnesses in
both, of which `topk 6` per bin keeps 11,661.  `k = 4`: the full sweep by the new binary and the
reduced sweep report the **same minimum at the same bin**, `997092898038/10^12` at `k = 411`,
which is the old binary's `0.9970929` from the loop's log; the argmin lies below 45°, so the
reduction changes neither the decision nor the number here.  (In general the reduced sweep's
minimum can differ from the full sweep's when the argmin lies above 45°, since the net is not
mirror-symmetric about 45°; both are valid lower bounds and the decision is the same.)  The
witness file of the reduced sweep is the `[0°, 45°]` half, 1,021 of the 2,131 lines.

*Speed-ups at 3 threads.*

| lever | instance | before | after | measured |
|---|---|---|---|---|
| 2 (+3): incremental credit | `1110`, whole sweep | 1699.5 s | 416.5 s | **4.08×** |
| 2 (+3): incremental credit | `k = 4`, whole full sweep | ~66,600 s (log 64,362 s) | 5158 s | **12.9×** (12.5× against the log) |
| 2: incremental credit | `k = 4`, per bin, 20 sampled bins, 1 thread | 1999 s | 144.9 s | **13.8×** (5.4–19.0× per bin) |
| 1: D4 reduction | `k = 4`, whole sweep | 5158 s (2000 bins) | 2367 s (830 bins) | **2.18×** (2.41× on bins; the `[0°, 45°]` bins are the dearer ones, 2.85 s against 2.58 s per bin-thread) |
| 3: dynamic scheduling | `k = 4`, bins 400..439, 4 threads | 81.2 s (static split) | 79.5 s | **1.02×** |
| **all three** | **`k = 4`, whole sweep** | **~66,600 s = 18.5 h** | **2367 s = 39.5 min** | **28.1×** (27.2× against the log) |
| all | `1110` (no symmetry available) | 1699.5 s | 416.5 s | 4.08× |

The review's estimates were 2.4× / ~5× / ~1.3×.  The profile said the anchor credit was 94 % of
the sweep and removable in full, and lever 2 came out at 13–14× rather than 5×; lever 1 is 2.18×
of wall rather than 2.41× because the low bins are the expensive ones; and lever 3 is worth 2 %
rather than 30 %: on these probes the bins are near-uniform in cost (97–117 s, only bin 0 is
cheap), so the static split was already nearly balanced.  It stays because it is free and makes
the output schedule-independent.

### 3.3 Threads (bins 400..439 of `k = 4`, new binary, full sweep, `runs/scal_*`)

| threads | wall | per-thread (bins / wall) | efficiency |
|---|---|---|---|
| 1 | 312.9 s | 40 / 312.8 s | 100 % |
| 2 | 157.0 s | 20 / 156.8 s; 20 / 157.0 s | 100 % |
| 4 (dynamic) | 79.5 s | 10 / 78.8; 10 / 79.2; 10 / 79.3; 10 / 79.5 s | 98 % |
| 4 (static split, `VERIFY_STATIC=1`) | 81.2 s | 10 / 79.9; 10 / 80.1; 10 / 80.2; 10 / 81.2 s | 96 % |

**Projection to 16 threads.**  The sweep is embarrassingly parallel over bins and scales at 98 %
to 4 threads here; the review measured 1.8× from 8 to 16 threads on the old binary (16 physical
cores; 16 → 32 gains nothing), so 16 threads are taken as 14.4 single-thread equivalents
(90 %).  On that basis, with the measured CPU-seconds (3 threads × wall × 0.98):

| sweep | CPU-s | at 4 threads (measured/derived) | **at 16 threads (projected)** |
|---|---|---|---|
| `k = 4`, old binary, full | 199,900 | ~51,000 s = 14.2 h | ~13,900 s = 3.9 h |
| `k = 4`, new, full sweep forced | 15,165 | ~3,870 s = 64 min | ~1,050 s = 17.5 min |
| `k = 4`, new, D4 reduction | 6,960 | ~1,780 s = 30 min | **~480 s = 8 min** |
| `1110`, old binary | 4,997 | ~1,275 s | ~350 s |
| `1110`, new | 1,225 | ~312 s | **~85 s** |

So a `k = 4` round's exact separation goes from the 12–18 h the loop was paying to about 8
minutes at 16 threads (30 minutes at 4), and the `1110` leaf's from ~6 to ~1.5 minutes.

### 3.4 `xcheck.py`, bin by bin

`certificates/s12_anchorclique_demo_3.9318.txt` at `N = 2000` (the demonstration needs
`N = 6000` to verify, so this is a rejected file, which is the harder case for agreement):
`python3 xcheck.py … 2000 --all --n 12 -v -j 1` (exhaustive, 2000 bins, `runs/xc_py.txt`)
against the Rust per-bin minima from 2000 `VERIFY_BINS=k:k` runs (`runs/xc_rust.txt`), compared
as exact fractions (`runs/xccmp.py`): **2000 of 2000 bins agree, 0 differ**; minimum
`9801147/10^7 = 0.9801147` at bin 111 in both.  The rejection suite's §27 (six anchor files,
`N = 400`, both exhaustive) and §29–30 add the D4-closed families and the cross-checked ones.

### 3.5 The shipped certificates

`./verify.sh` from the repository root with `$(nproc)` pinned to 4 (`runs/verify_sh.out`):
exit 0, the same 24 `VERIFIED:` verdicts as before (every certificate in `certificates/`, the
`N = 12000` re-runs, the box-clique and anchor-clique demonstrations, `n = 11`), the two
`xcheck.py` exhaustive runs, the 172 rejection checks and every JSON round-trip; 4 min 11 s
wall.  None of the shipped certificates is a D4-closed anchor family (the anchor demonstration
carries one clique, not an orbit), so none of them changes path, and `tests/bitid.sh` pins
their output byte for byte.

## 4. Reproduce

```sh
cd verify && cargo build --release && cd ..
./tests/rejection_tests.sh                       # 172 checks; §29-31 are the levers' traps
sh tests/bitid.sh runs/verify_ref                # ALL IDENTICAL against the pre-change binary
VERIFY_STATS=1 VERIFY_BINS=400:400 verify/target/release/verify runs/branch_J16i_probe.txt 12 2000 1 3 out.sep
VERIFY_ANCHOR_XCHECK=1 VERIFY_FULLSWEEP=1 VERIFY_BINS=0:0 verify/target/release/verify runs/branch_J16i_probe.txt 12 2000 1 3 out.sep
verify/target/release/verify runs/branch_J16i_probe.txt 12 2000 4 3 probe.sep      # the reduced sweep, 830 bins
sh runs/val.sh                                   # the whole-probe runs of section 3
```
