# s(45) = 7

**Claim.**  No 45 unit squares fit in a square of side less than 7.  Since 49 unit squares tile the 7 × 7
square, **s(45) = 7**.  With `s(21) = 5` ([`../s21/`](../s21/README.md)) and `s(32) = 6` ([`../s32/`](../s32/README.md))
this is the third exact value of `s(k² − 4)` for `k ≥ 4`.  Previous lower bounds: `1389/200 = 6.945` (wand125,
2026-09-25/26, rectangle-density certificate); `1 + √34 = 6.8310` (Nagamochi 2005, a general bound; its published proof is incomplete, Karakuş, arXiv:2609.37410).

Computer-assisted, not peer reviewed, **not (yet) in Lean**.  The proof is designed to be re-checked: `./verify.sh`
here.

## The proof in one paragraph

`s45_mixed_cover_7.txt` is a **mixed cover** (format: [`../s21/FORMAT.md`](../s21/FORMAT.md), the same as for
`s(21)`): 19,989 weighted points plus mass spread uniformly along 3,912 segments of the interior grid lines
`x, y ∈ {1, …, 6}`, total `2238676387 / 5·10⁷ = 44.77352774 < 45`, such that **every closed unit square contained
in `[0,7]²`, at any position and any angle, captures mass at least 1** (a point on the square's boundary counts, and
a segment lying along one of its edges counts in full).  If 45 unit squares with disjoint interiors fit in a
square of side `s < 7`, scale their centres by `7/s > 1`: they become 45 pairwise *disjoint* closed unit squares in
`[0,7]²`, so together they would capture mass `≥ 45 > 44.77`.  Contradiction.  As for `s(21)`, the mass on the grid
lines (`32.8` of the `44.77`) is what makes this possible at margin zero: a square sitting exactly on a grid cell
gets its whole boundary's line mass, and a slightly tilted or shifted one loses part of an edge on one line and
gains the complementary part on the parallel line one unit away.  The cover is invariant under the symmetries of
the square, so it suffices to check squares with centre in `[0, 7/2]²` and angle in `[0°, 45°]`; that check is an
exhaustive, exact subdivision of pose space by a computer, done by two separately written checkers (shared point-test lineage: below), and the
second one also checks the whole pose space without using the symmetry.

## What is checked by what

| step | how | trust |
|---|---|---|
| cover ⇒ `s(45) ≥ 7` (scaling, disjointness), and the D4 reduction for measures | on paper: [`../s21/FORMAT.md`](../s21/FORMAT.md) ("Why that bounds `s(n)`").  The same statements are proved in Lean for every side `m` (`lean/Sqpack/MixedMeasure.lean`: `not_packs_of_measure`, `d4_reduction_measure_u`), but **no `s(45)` data file or top theorem exists**: nothing about *this* cover is checked in Lean | paper; generic Lean lemmas |
| `s(45) ≤ 7` | the `7 × 7` grid | trivial |
| every admissible closed unit square with centre in `[0,7/2]²`, `θ = 2 arctan u`, `u ∈ [0, ½]`, captures mass `≥ 1` | **`search/zm_mixed.py --d4 --cert-mode`**: all 78,400 root boxes certified, 0 uncertified (`zm_mixed_d4/`); the checker files are byte-identical to those of the `s(21)` bundle | exact `Fraction`/integer arithmetic; the program is not formally verified |
| the same, separately written (shared point-test lineage, below), and the whole pose space without the symmetry argument | **`verify2/src/bin/zmx2.rs` (`zmx2`, Rust)**, the binary of the `s(21)` bundle (same source and binary sha256): `--d4` 4,900 / 4,900 roots (`zmx2_d4/`), and `--full` (no symmetry used) 39,200 / 39,200 roots (`zmx2_full/`); 0 uncertified | exact integers for points and masses; chord-end geometry in outward-rounded binary64 intervals (caveat below) |
| total `< 45`, D4 invariance, well-formedness | `search/mixed_records.py cover` (its own parser) and `zmx2 d4` in `verify.sh`; `zm_mixed.py` and `zmx2` also check invariance before a D4 run | scripts |

**Note (2026-09-30).**  `zmx2.rs` has since gained the opt-in flag `--sym-atoms` (`search/ZMX2.md` §4.9, §12; source sha256 `92a4cfe8…`, binary `ed31d3ee…`).  Without the flag its decisions are unchanged: `verify.sh` rebuilds `zmx2` from the current source and its fresh run reproduces the shipped census root for root.  The shipped `zmx2` records were made with the earlier source and binary named in their manifests.

## The float caveat of `zmx2`, plainly

`zm_mixed.py` certifies with rational arithmetic only (floats only choose which exact test to try).  `zmx2`
decides point containment and adds up all masses in exact integers, but it encloses the *positions of chord
ends* (where a grid line crosses the rotated square's edges) with IEEE-754 binary64 interval arithmetic, every
operation widened outward by one ulp and the result rounded outward onto an integer grid of pitch `1/(1000·2³⁰)`
(`search/ZMX2.md` §5, Lemma R).  So its verdict is a rigorous enclosure *provided* the hardware and compiler give
correctly rounded `+ − × ÷ √` without contraction or reordering (IEEE-754 requires the former; Rust guarantees the
latter; x86-64 SSE2).  It is a different kind of trust from `zm_mixed.py`'s, which is the point of having both.

## Files

| file | |
|---|---|
| `s45_mixed_cover_7.txt` | the cover, mixed format v1: `s = 7`, `D = 1000`, `W = 10⁸`, 19,989 points, 3,912 axis-parallel segments on the 12 interior grid lines, no polygons, all integers, exactly D4-invariant, no point on a segment line.  (It is the second cutting-plane round of the line-cover LP, `lc_m7A_r1`, scaled by `409/400` with masses rounded up; its first comment line says so.  `search/S45_COVER.md` §10.) |
| `zm_mixed_d4/` | the `zm_mixed.py` run: `roots.jsonl` (a header with the sha256 of the checker, reader, `zeromargin.py` and cover and all settings, then one record per root: census, uncertified boxes, CPU), `manifest.json` (header, command line, result, records sha256), `run.log`, and `checker/`, the exact files that ran.  How the records came to be here: below |
| `zmx2_d4/`, `zmx2_full/` | the `zmx2` runs on the bundled file: `roots.log` (one line per root), `run.out`, `manifest.txt` (git HEAD, source/binary/input sha256, command, census, verdict) |
| `verify.sh` | re-check; below |
| `SHA256SUMS` | |

**Provenance of `zm_mixed_d4/`.**  The sweep was run once (2026-09-27, 16.5 CPU-h, 12 processes) on the same bytes
under their working name `runs/s45_mixed_candidate_7.txt` (sha256 `5da180d1…7c6c`), with `--resume` records.  For
the bundle those records were not recomputed: `search/s45_cert_runs.sh import` copied them verbatim, changed only
the `input` *path* in the header line to `certificates/s45/s45_mixed_cover_7.txt`, and ran `zm_mixed.py` itself
with `--resume` on them, which recomputes its header from the files (all four sha256, settings, total), refuses a
file whose header differs, found all 78,400 roots done and wrote `manifest.json` (hence its `wall_s` of 0 and the
`--nproc 1` in its argv).  `run.log` is the original run's log followed by that resume's log.  `verify.sh --full`
recomputes everything from scratch.  The two `zmx2` runs were made fresh from the bundled file (`search/s45_cert_runs.sh
zmx2`); their census equals that of the earlier runs on the working copy, root for root.

## Re-checking

```
certificates/s45/verify.sh          # ~30 s on 6 cores: hashes; exact total and D4 invariance; shipped zm_mixed.py
                                    #   and zmx2 records re-summarised from scratch; checker files = the s(21)
                                    #   bundle's; a FRESH zmx2 run over all poses, no symmetry (39,200 roots,
                                    #   ~150 CPU-s) compared root for root with the shipped one
certificates/s45/verify.sh --full   # also re-runs zm_mixed.py --d4 --cert-mode over the whole D4 region
                                    #   (~16.5 CPU-h; about 1.5 h on 12 processes)
```

`zm_mixed.py` settings (those of the `s(21)` bundle): centre cells of pitch `1/20` over `[0, 7/2]²`, 16 bins of
`u = tan(θ/2)` of width `1/32` over `[0, ½]` (`θ ≤ 53.13°`, more than the `45°` needed), subdivision to depth 24,
zeromargin's CHAIN from depth 0, `--cert-mode` (Corollary T′ and the polygon code unreachable; the cover has neither
points on its lines nor polygons).  Totals: 78,400 roots, 437,510 boxes, max depth 19, leaves ADM 122,892 /
CHAIN 45,866 / SPLIT 54,276 / PIECE 9,458 / EMPTY 25,463 / UNCERTIFIED 0, 16.5 CPU-h (83 min on 12 processes).
Every acceptance is a `Fraction` or integer test; the lemmas are proved in `search/ZM_MIXED.md` §2; `zeromargin.py`
(its point primitives, imported unchanged) is the file pinned by the `s(32)` bundle.

`zmx2` settings: centre cells of pitch `1/10`, 4 bins of `u` of width `1/8`, depth ≤ 40.  `--d4`: 4,900 roots,
2,071,984 boxes, max depth 24, 19 CPU-s.  `--full`: the cover on `u ∈ [0, ½]` plus its mirror image `y ↦ 7 − y`
(which covers `θ ∈ [36.87°, 90°]`), 39,200 roots, 16,648,752 boxes, max depth 24, 144 CPU-s.  Lemmas and tests:
`search/ZMX2.md`.

How the cover was found (line-cover LP with germ rows, exact threshold by `zmx2` bisection, the `409/400` scaling):
`search/S45_COVER.md` §10.  By that bisection the unscaled round cover verifies at factor `1.0152930` and is
refused at `1.0152344`, so this cover sits about `0.7 %` above its exact threshold (the lower end is a refusal,
not a proof of a dip, since `zmx2` may lose up to `≈ 0.1 %`).

## What is not machine-verified

* **No Lean for `s(45)`.**  Unlike `s(21)` and `s(32)`, there is no `S45Data.lean` and no top theorem
  `minSide 45 = 7`; the reduction from the cover to the packing bound is on paper here (it is the `m = 7` case of
  lemmas proved in Lean for general `m`).
* The checker programs (subdivision, bookkeeping, parsing, the interval arithmetic).  Their lemmas are proved on
  paper (`ZM_MIXED.md` §2, `ZMX2.md` §2–7), not in Lean.  The two checkers share no code: `zmx2` was written
  without opening `zm_mixed.py` or its write-up, with its own parser, its own lemmas and its own arithmetic.
  They are not independent in every respect: `zmx2`'s point test (Lemma P, `ZMX2.md` §2) uses the violation
  polynomials `G0`–`G3` of `zeromargin.py`'s write-up (`RUNG2.md` §6.2), and `zm_mixed.py` tests points by calling
  `zeromargin.py` itself, so the two share a point-test lineage and an error in that formulation would be common to
  both.  Their treatments of the line masses were derived separately.
* The run records: `verify.sh` re-summarises them; `--full` regenerates the `zm_mixed.py` one, and the default
  tier regenerates the `zmx2 --full` one and compares it root for root.

## Review record

Both checkers are the exact files audited for `s(21)` (same sha256: `zm_mixed.py 1fd20346…`, `mixed_cover.py
bb89de15…`, `zeromargin.py 640fe453…`; `zmx2.rs 6b7f0f79…`, binary `0247012e…`); see the review record in
[`../s21/README.md`](../s21/README.md) (`search/ZM_MIXED_AUDIT.md`, `search/ZMX2_AUDIT.md`).  Nothing specific to this
cover has been audited separately: no adversarial component tests or holed-cover rejection runs were made on it,
and there is no independent review of this bundle yet.  The two checkers reach the same verdict on it.

**Re-run 2026-09-29 (checker `zm_mixed.py` `1fd20346…`, was `ee3e2915…`).**  One guard added after the audit: `region_phi` (Lemma R) now refuses a zero-width angle bin instead of returning an unproved `EMPTY` (latent bug B1 of the lean-segments audit, `notes/lean-segments.md` §3; unreachable in the earlier run).  The whole `--d4 --cert-mode` run was repeated with the new file (records and `checker/` replaced); its census is identical to the previous run's root for root (every root, every leaf count).
