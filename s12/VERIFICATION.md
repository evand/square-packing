# Verification log

## Result
`s(12) >= 15680/3951 = 3.968616`  — 12 unit squares cannot be packed into any square of side < 15680/3951.

Second result, same machinery, separate certificate and separate checker:
`s(13) = 4`, case-free (Bentz 2010 re-proved without its six-leaf case analysis).  See
**s(13) = 4 (rung 2) (2026-09-12)** below.

## Certificates
| file | s | points | total weight | verify N=6000 min | N=12000 | xcheck --all N=6000 |
|---|---|---|---|---|---|---|
| `s12_lower_3.9686.txt` | 15680/3951 = 3.968616 | 1736 | 11.9738036 | 10000056/10^7 VERIFIED | VERIFIED | 1250007/1250000 VERIFIED (5 min, 32 cores) |
| `s12_lower_3.9676.txt` | 980/247 = 3.967611 | 764 | 11.9962288 | 10000033/10^7 VERIFIED | VERIFIED | VERIFIED |
| `s12_lower_3.931795_sparse.txt` | 3920/997 | 224 | 11.9834372 | 10000022/10^7 VERIFIED | VERIFIED | VERIFIED |
| `s12_lower_3.931795.txt` (original) | 3920/997 | 788 | 14916233/1250000 = 11.9329864 | 10000023/10^7 VERIFIED | VERIFIED | VERIFIED |

All coordinates have denominator 2·s_den (1994 for the 3920/997 files, 494 and 3951 for the
others); weights have denominator 10^7.  How the three newer files were produced (weights
re-optimised at fixed points with the exact verifier as separation oracle; reweighted-L1
sparsification; column generation with exact reduced-cost pricing) is in `search/TIGHTEN.md`.
The checks were re-run independently of the search session on the merged tree.

## s(11) certificate (2026-08-26)
`certificates/s11_lower_3.8143.txt` — `s(11) >= 3040/797 = 3.814304`; 680 points, coordinates
over 3985, weights over 10^7, total 27036677/2500000 = 10.8146708 < 11.  Checks: `verify` with
n = 11 at N=6000 (min 10000042/10^7) and N=12000: VERIFIED; `xcheck.py --all --n 11` at N=6000:
identical minimum; `scale_to_critical.py --n 11 --N 6000`: already critical (D = 3985).  The
N=6000 run was repeated by hand on the merged tree.  Search: `tighten.py --n 11` from the
n = 12 certificates rescaled, with column generation (`search/N11.md`).  Now in `verify.sh`.

## Verifier diagnostic mode (2026-08-26)
`verify/src/main.rs` gained an opt-in mode `TIGHT_DUMP=<path> TIGHT_THRESH=<int>` that lists the
arrangement cells below a captured-weight threshold (`search/TIGHTSET.md`).  With the variable
unset the code path and output are unchanged: the default output for the main certificate at
N=6000 was diffed against the previous binary, and `verify.sh` plus the 42 rejection tests pass.

## Witness file: credited cliques (2026-09-07)
`verify/src/main.rs` now records, on each line of the optional witness file (`topk > 0`), the
anchor cliques its sweep credited to the cell that witness came from — two extra columns, a count
and the file indices, and **only** for a certificate carrying an `anchors` block.  The verdict
path, the exit codes and stdout are untouched, and the witness file of every certificate without
such a block is byte-identical: `tests/bitid.sh` against the previous binary is ALL IDENTICAL on
ten cases (plain, witness, `TIGHT_DUMP`, branch trailer, box cliques), and on an anchor-clique file
the new witness lines cut back to five fields reproduce the old ones exactly.  The rejection suite
is 138 checks, 0 failures, 0 panics, and `verify.sh` is exit 0 with 24 VERIFIED verdicts.  Why: the cutting-plane loop was building its LP rows with a
*pose* predicate while the verifier credits a *cell*, so the LP believed rows the verifier keeps
failing — `search/WITNESS.md` has the defect, the fix and the numbers.

## Verifier speed-up (2026-09-08)
`verify/src/main.rs`: three changes for the anchor-clique certificates of the cutting-plane loop
(`search/VERIFYSPEED.md`): the `[0°, 45°]` reduction is applied again when the atoms AND the
anchor-clique family are D4-invariant (a new exact closure check, `anchors_d4_closed`); the anchor
credit of a strip's cells is computed once per (strip, piece) as a cell range instead of once per
(cell, piece), with the per-cell test kept as an optional cross-check (`VERIFY_ANCHOR_XCHECK=1`);
bins are scheduled dynamically over the threads, and the witness file and the reported minimum
are now independent of the thread count.  Checks: `tests/bitid.sh` against the previous binary is
ALL IDENTICAL (11 cases, single thread, plain / witness / `TIGHT_DUMP` / branch / box-clique /
anchor-clique modes); `./verify.sh` at 4 threads is exit 0 with the same 24 VERIFIED verdicts; the rejection suite is 172
checks (was 138; the 34 new ones are the levers' over-credit traps), 0 failures, 0 panics; on the
`k = 4` and `1110` probes the per-bin output (stdout and witness file) is byte-identical to the
old binary on every bin compared, and the incremental credit agrees with the per-cell test at
every cell of those bins.  Whole-probe before/after numbers are in `search/VERIFYSPEED.md`.

## s(13) = 4 (rung 2) (2026-09-12)

A second, separate result: one weighted closed cover of `[0,4]^2` of total weight < 13, so no 13
unit squares fit below side 4, and 16 of them tile `[0,4]^2`, so `s(13) = 4` (Bentz 2010) with no
case analysis.  Note: `verify/` **cannot** check this file — at the container itself the margin is
exactly 0 and its angle-net erosion is linear against a quadratic margin (`search/ZEROMARGIN.md`
§1), which is why `verify2/` exists.

| file | container | points | total weight | sha256 |
|---|---|---|---|---|
| `certificates/rung2/s13_closed_cover_4.txt` | `[0,4]^2` | 3621 | `2591194431/200000000 = 12.955972155` | `ea303acea08cc17a13cecc24d3714c2df409f91eba048cd5546050ed064b53ed` |

Coordinates have denominator `D = 1000`, weights `10^9`; every number in the file is an integer.
Now pinned in `certificates/SHA256SUMS`.

**Checker 1 — `search/zeromargin.py`** (Python, `fractions.Fraction` throughout; floats only as
pre-filters, which can lose a certification but never create one).  Symmetry-reduced domain
(`c_x in [0,4]`, `c_y in [0,2]`, `u = tan(th/2) in [0,1/2]`), legitimate because the checker first
verifies **exactly** that the point set is invariant under `x -> 4-x` and `y -> 4-y`
(`Checker.symmetric`) and refuses to run reduced otherwise.

```
python3 search/zeromargin.py cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --nproc 8 --disj --chain-from 0 --dump runs/leaves.txt

done in 6138s: boxes 16872, max depth 13
  leaves: ADM 2867  CORE 0  P1 0  MIX 0  CHAIN 5320  TRI 0  EMPTY 3449  UNCERTIFIED 0
VERIFIED
```

1 h 42 min on 8 processes; depth limit 18 never reached (max 13), so the subdivision terminated on
its own.  Re-run from the committed path it reproduces the census to the last box
(`done in 5828s: boxes 16872 ... ADM 2867 / CHAIN 5320 / EMPTY 3449 / UNCERTIFIED 0, VERIFIED`),
so the check is deterministic and the shipped file is the one that verifies (`search/RUNG2.md` §0).

**Checker 2 — `verify2/zmcheck`** (Rust, exact `i128` on the whole load-bearing path; its own
subdivision, its own primitive set, no shared code, written from the lemmas rather than from the
Python).  **Full** domain: `c_x, c_y in [0,4]`, `u in [0,1]` (`th in [0,90]` deg), no symmetry
assumed or checked.  `runs/zmcheck_main_2026-09-12.log`, main checkout:

```
verify2/target/release/zmcheck cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --threads 8

certificate certificates/rung2/s13_closed_cover_4.txt: container [0,4]^2, 3621 points, D=1000 W=1000000000
total weight = 12955972155/1000000000 = 12.955972155
12800 root boxes (pitch 1/10 in x,y; 8 bins of u in [0,1]); depth limit 18; disj true; 8 threads
done in 1065s: boxes 30258, max depth 10
  leaves: ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0
VERIFIED: every closed unit square in [0,4]^2 captures weight >= 1; total weight 12955972155/1000000000 = 12.955972155
```

17 min on 8 threads; depth limit never reached (max 10 of 18).  Timed again for this write-up on
the same machine, same binary, same file: `done in 1064s` at 8 threads and `done in 2015s` (33 min
35 s) at 4, both `boxes 30258, max depth 10` with
`ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0` — the census is identical in all three runs, so
neither the verdict nor the box count depends on the thread count
(`runs/zmcheck_t8_2026-09-12.log`, `runs/zmcheck_t4_2026-09-12.log`).  4 threads is the GitHub
runner's core count, and 33.6 min is inside the budget, so `verify.sh` runs the full sweep on
every push; see `search/S13_WRITEUP.md`.

Two checkers, two subdivisions, two primitive sets, two domains, `0 uncertified` on the same file.
The disjunctive primitive carries 65 % of the non-empty leaves in each (`CHAIN` 5320 of 8187;
`DISJ` 9477 of 14591) — which `search/RUNG2.md` Theorem 1 predicts, since no cover of `[0,4]^2`
below weight 16 can be certified by fixed-witness primitives at any depth.

**Rejection tests.**  `runs/rung2_rejection_2026-09-12.log`:

```
./tests/rung2/rejection_tests.sh

23 passed, 0 failed, 0 panics
REJ_EXIT=0
```

Five mutations of the certificate and their controls (1 % weight cut: still valid; 5 % cut, halved
weights, set scaled by 0.995, set translated by 0.01: a violating pose is **exhibited** each time,
exit 1 with `*** VIOLATION`; one point deleted, one point moved by 0.01: still valid), the two
invalid covers from the development history at their exact violating weights (`0.970282351`,
`0.9420217`), and nine malformed inputs that must exit 2 with `ERROR:` and **no verdict word**.
The mutation sweeps run over the restricted band `c_x in [0.5,0.6]`, `c_y in [1,2]` — which
contains the pose `(1/2, 3/2, 0)` that `RUNG2.md` §2 identifies as the container's worst
monotone-witness pose — and are compared on uncertified-box count, because a restricted sweep
prints `PARTIAL SWEEP` and can never print `VERIFIED`.  One test, the `closed4_best_x103`
violation, needs a development-history file kept outside the repo and prints `skip` (20 passed)
when it is absent; that is the case in CI.

**Lean.**  `lean/Sqpack/ZeroMargin.lean`, 970 lines, Mathlib `v4.33.1`
(`runs/lean_build_2026-09-12.log`):

```
cd lean && lake build
Build completed successfully (8710 jobs).
BUILD_EXIT=0
```

0 `sorry`; `lean/Axioms.lean` prints `#print axioms` for **36** theorems and every one is
`[propext, Classical.choice, Quot.sound]`.  Formalised: `sq_subset_box_iff` (admissibility, both
directions, all four sides), Lemma A (eight-corner form), Lemma B (the degree-4 coefficient rows),
Lemma C (both branches), Lemma E (the maximum is attained, so the `_gmax <= 0` test is an
equivalence), Lemmas F–H and the chain covering, `clip_bin_no_loss`.  **Not** formalised, in
either checker: the subdivision, the exhaustiveness/control flow, the weight bookkeeping, the
certificate parser, the float pre-filters (`notes/lean-zeromargin.md`).

**In `verify.sh` / CI.**  `verify.sh` now builds `verify2` and runs the `zmcheck` sweep at
`--depth 18 --threads "$(nproc)"`, failing on `NOT VERIFIED` or `PARTIAL SWEEP` exactly as `chk`
does for `verify`, then runs `tests/rung2/rejection_tests.sh`.  The `zeromargin.py` sweep is a
commented slow path next to it.  `search/S13_WRITEUP.md` has the timing table and the CI decision.

`./verify.sh` run end to end at 16 cores after these changes
(`runs/verify_full_2026-09-12.log`, 1054 s, exit 0): **25 `VERIFIED` verdicts** (was 24; the new
one is the rung-2 sweep, `done in 626s`, same census), no `NOT VERIFIED`, no `PARTIAL SWEEP`, the
rung-2 suite `23 passed, 0 failed, 0 panics`, the main suite `172 passed, 0 failed, 0 panics`, and
every `points.json` round-trip byte-identical.  `sha256sum -c certificates/SHA256SUMS`: 36 of 36 OK
before and after.

## Checks performed on the original 788-point certificate (2026-08-23)
## Companion certificate (weaker but human-readable)
`certificates/s12_56points_3.8.txt` — 56 points in [0,19/5]^2 (scales up to [0,1520/397]^2), weight 1/5 each (total 56/5 = 11.2).
Statement: every closed unit square inside [0,3.8]^2, at any angle, contains at least 5 of the
56 points; 12 disjoint squares would need 60.  Checks: `verify` at N=2000/4000/8000, and the
independent exact Python re-check `xcheck.py` over every bin at N=2000 (829 bins, identical
to the Rust per-bin minima) and N=8000 (3314 bins): minimum exactly 1.

## Uniform certificates (2026-08-25)
`certificates/s12_uniform_<k>of<m>_<s>.txt`, k = 1..8: m points, weight 1/k each, m < 12k,
every closed unit square contains >= k of them.  Headline: 81 points in [0,35/9]^2 with
k = 7, i.e. s(12) >= 35/9 = 3.888889.  Each file is written at its critical container.
Checks: `verify` at N=2000 and N=8000 (all nine files), `xcheck.py --all` at N=2000 (all
nine) and at N=4000 for the 81-point set, JSON round-trip for all.  Search method (ILP over
D4 orbits + verifier-driven polish) in `search/uniform/UNIFORM.md`.

## Formalisation
`lean/Sqpack/Basic.lean` (Lean 4 + Mathlib) proves the reduction:
a weighted set whose closed unit squares all carry weight >= 1 bounds the number of squares of
side L>1 packable with disjoint interiors by the total weight; plus the scaling lemmas.
0 sorries; axioms: propext, Classical.choice, Quot.sound.

## Verifier defects found by the rejection tests (2026-08-25), and fixed

Recorded because a verifier's history matters as much as its current state.  None affects
the shipped certificates or the bound; each was found by writing a test that the verifier
should fail, and watching it not fail.

1. **Negative weights were not rejected.**  The reduction `n <= sum w` requires `w >= 0`
   (the Lean proof's hypothesis `hw`).  The old verifier accepted the 56-point set plus one
   point of weight `-2` placed *outside* the container, claimed for `n = 11`: the total
   dropped to 10.8 while the covering was untouched, and it printed VERIFIED.  Now an
   `ERROR`.  The shipped certificates have only positive weights.
2. **Centre box too small above 45°.**  For each angle bin the admissible-centre box used
   the bounding-box width `w = cos θ + sin θ` at `θ_k` only.  Below 45° that is the bin's
   minimum (correct); above 45° `w` decreases, so it was the bin's maximum and a strip of
   width `(w(θ_k) − w(θ_{k+1}))/2` along the container edges was never checked.  Found by
   the exhaustive Python checker disagreeing with the Rust one on a non-symmetric mutant
   (bin k=1530 of N=2000: Rust 2/5, exact 1/5).  Only the `[0,90°)` path was affected; both
   shipped certificates are D4-symmetric and use `[0,45°]` only.  Fixed to
   `min(w(θ_k), w(θ_{k+1}))`.  The tests now pin the per-bin minimum on that mutant.
3. **Panics on malformed input** (header-only file, non-integer token, short point line,
   point count mismatch, empty or missing file, `s_den ∤ s_num·D`).  A panic is not a
   rejection.  All now exit 2 with `ERROR:` and no verdict word.
4. Points outside the container were not rejected; now `ERROR`.

5. **`verify.sh` (hence CI) did not fail on NOT VERIFIED.**  `verify` exits 0 on a verdict
   either way, and `set -e` only sees exit codes, so a rejected certificate would have left
   the badge green.  Every verdict is now checked explicitly; tested by breaking a
   certificate's header and watching `verify.sh` exit 1.

`tests/rejection_tests.sh` now has 42 checks (was 7); against the pre-fix binary 15 fail and
8 panic.
