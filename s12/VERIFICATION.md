# Verification log

## Result
`s(12) >= 15680/3951 = 3.968616`  — 12 unit squares cannot be packed into any square of side < 15680/3951.

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
