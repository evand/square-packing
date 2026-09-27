# s(21) = 5

**Claim.**  No 21 unit squares fit in a square of side less than 5.  Since 25 unit squares tile the 5 × 5
square, **s(21) = 5**.  With `s(32) = 6` ([`../s32/`](../s32/README.md)) this is the second exact value of
`s(k² − 4)` for `k ≥ 4`.  Previous lower bounds: `s(21) ≥ 5000/1001 = 4.995005` (this project, 2026-09-23;
`s21_lower_4.9950.txt` here); `399/80 = 4.9875` and before it `249/50 = 4.98` (wand125, 2026); `122/25 = 4.88`
(jlevy, 2026); `4.7438` (Friedman, DS7 survey); `1 + √14 = 4.7417` (Nagamochi 2005, a general bound).

Computer-assisted, not peer reviewed.  The proof is designed to be re-checked: `./verify.sh` here.

## The proof in one paragraph

`s21_mixed_cover_5.txt` is a **mixed cover** ([`FORMAT.md`](FORMAT.md)): 7,536 weighted points plus mass spread
uniformly along 1,872 segments of the interior grid lines `x, y ∈ {1, 2, 3, 4}`, total
`522368729933 / 25·10⁹ = 20.894749197 < 21`, such that **every closed unit square contained in `[0,5]²`, at any
position and any angle, captures mass at least 1** (a point on the square's boundary counts, and a segment
lying along one of its edges counts in full).  If 21 unit squares with disjoint interiors fit in a square of side
`s < 5`, scale their centres by `5/s > 1`: they become 21 pairwise *disjoint* closed unit squares in `[0,5]²`, so
together they would capture mass `≥ 21 > 20.89`.  Contradiction.  The mass on the grid lines is what makes this
possible at margin zero: a square sitting exactly on a grid cell gets its whole boundary's line mass, a slightly
tilted one loses part of an edge on one line and gains the complementary part on the parallel line one unit
away, so nothing can slip between discrete points.  The cover is invariant under the symmetries of the square,
so it suffices to check squares with centre in `[0, 5/2]²` and angle in `[0°, 45°]`; that check is an
exhaustive, exact subdivision of pose space by a computer, done by two independently written checkers.

## What is checked by what

| step | how | trust |
|---|---|---|
| cover ⇒ `s(21) ≥ 5` for any measure (scaling, disjointness), `s(21) ≤ 5` (grid), D4 reduction for measures, invariance and total of *this* cover | **Lean 4 + Mathlib**, `lean/Sqpack/{MixedMeasure,SegTree,S21Data,S21}.lean`: `s21_eq_five_of_checker : S21CheckerCover → minSide 21 = 5`; standard axioms only, no `sorry`, no `native_decide` (`notes/lean-s21.md`) | kernel |
| `S21CheckerCover`: every admissible closed unit square with centre in `[0,5/2]²`, `θ = 2 arctan u`, `u ∈ [0, ½]`, captures mass `≥ 1` | **`search/zm_mixed.py --d4 --cert-mode`**: all 40,000 root boxes certified, 0 uncertified (`zm_mixed_d4/`) | exact `Fraction`/integer arithmetic; the program is not formally verified |
| the same, independently, and the whole pose space without the symmetry argument | **`verify2/src/bin/zmx2.rs` (`zmx2`, Rust)**, written without reading `zm_mixed.py` or its write-up: `--d4` 2,500 / 2,500 roots (`zmx2_d4/`), and `--full` (no symmetry used) 20,000 / 20,000 roots (`zmx2_full/`); 0 uncertified | exact integers for points and masses; chord-end geometry in outward-rounded binary64 intervals (caveat below) |
| that the Lean data is this file | `lean/scripts/gen_s21_data.py --check` (byte for byte), sha256 in the Lean file header | script |
| total `< 21`, D4 invariance, well-formedness | Lean (`S21Data.total_eq`, `S21Data.d4`) and again in `verify.sh` by `search/s21_records.py cover` (its own parser) and `zmx2 d4` | kernel; scripts |

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
| `s21_mixed_cover_5.txt` | the cover, mixed format v1 ([`FORMAT.md`](FORMAT.md)): `s = 5`, `D = 1000`, `W = 10¹¹`, 7,536 points, 1,872 segments, no polygons, all integers.  (It is the cover side's candidate scaled by `1003/1000`; its first comment line says so.) |
| `FORMAT.md` | the mixed format and why a valid file bounds `s(n)` |
| `zm_mixed_d4/` | the `zm_mixed.py` run: `roots.jsonl` (a header with the sha256 of the checker, reader, `zeromargin.py` and cover and all settings, then one record per root: census, uncertified boxes, CPU), `manifest.json` (header, command line, result, records sha256), `run.log`, and `checker/`, the exact files that ran |
| `zmx2_d4/`, `zmx2_full/` | the `zmx2` runs: `roots.log` (one line per root), `run.out`, `manifest.txt` (git HEAD, source/binary/input sha256, command, census, verdict) |
| `verify.sh` | re-check; below |
| `SHA256SUMS` | |
| `s21_lower_4.9950.txt`, `.json` | the earlier `s(21) ≥ 5000/1001` point certificate (positive margin, `search/S21_LB.md`), superseded |

## Re-checking

```
certificates/s21/verify.sh          # ~1 min: hashes; exact total and D4 invariance; Lean data = this file;
                                    #   shipped zm_mixed.py and zmx2 records re-summarised from scratch;
                                    #   a FRESH zmx2 run over all poses, no symmetry (20,000 roots, ~100 CPU-s)
certificates/s21/verify.sh --full   # also re-runs zm_mixed.py --d4 --cert-mode over the whole D4 region
                                    #   (~20 CPU-h; about 2 h on 10 processes)
cd lean && lake build && lake env lean Axioms.lean                         # the Lean side
```

`zm_mixed.py` settings: centre cells of pitch `1/20` over `[0, 5/2]²`, 16 bins of `u = tan(θ/2)` of width `1/32`
over `[0, ½]` (`θ ≤ 53.13°`, more than the `45°` needed), subdivision to depth 24, zeromargin's CHAIN from depth
0, `--cert-mode` (Corollary T′ and the polygon code unreachable; the cover has neither points on its lines nor
polygons).  Totals: 40,000 roots, 461,204 boxes, max depth 21, leaves ADM 146,017 /
CHAIN 23,518 / SPLIT 53,951 / PIECE 9,058 / EMPTY 18,058 / UNCERTIFIED 0, 12.9 CPU-h (1.4 h on 10 processes).  Root for
root the census equals that of the first run of the same cover (`ZM_MIXED.md` §7.5, before the provenance changes).  Every acceptance is a `Fraction` or integer test; the lemmas are proved in
`search/ZM_MIXED.md` §2; `zeromargin.py` (its point primitives, imported unchanged) is the file pinned by the
`s(32)` bundle.

`zmx2` settings: centre cells of pitch `1/10`, 4 bins of `u` of width `1/8`, depth ≤ 40.  `--d4`: 2,500 roots,
1,826,222 boxes, max depth 26, 12 CPU-s (27 s in the shipped run, on a loaded machine).  `--full`: the cover on `u ∈ [0, ½]` plus its mirror image `y ↦ 5 − y`
(which covers `θ ∈ [36.87°, 90°]`), 20,000 roots, 14,709,448 boxes, max depth 27, 96 CPU-s (183 s shipped, loaded).  Both runs are by the
binary patched after the audit (`ZMX2.md` §11); root for root they have the same census as the runs before the patch.  Lemmas and tests:
`search/ZMX2.md`.

## What is not machine-verified

* The checker programs (subdivision, bookkeeping, parsing, the interval arithmetic).  Their lemmas are proved on
  paper (`ZM_MIXED.md` §2, `ZMX2.md` §2–7), not in Lean.  The two checkers share no code: `zmx2` was written
  from the format statement alone, with its own parser, its own lemmas (including its own
  derivation of the germ mechanism, a pair lemma for lines at distance 1) and its own arithmetic.
* The run records: `verify.sh` re-summarises them; `--full` regenerates the `zm_mixed.py` one, and the default
  tier regenerates the `zmx2 --full` one and compares it root for root.

## Review record (2026-09-27)

* **`zm_mixed.py`**: adversarial audit, `search/ZM_MIXED_AUDIT.md`: every lemma re-derived against its code, no
  soundness defect; ≈ 115 M exact checks of every certified *component* bound against exact masses at
  adversarial rational poses on this cover, 0 violations (several at slack exactly 0); three deliberately holed
  versions of this cover (a tilted dip, a SPLIT-heavy cell, a germ-pivot hole that exists only at `θ > 0`) each
  refused with its exact hole inside an uncertified box.  Its provenance should-fixes (shas and settings in the
  records, a header-checked `--resume`, a certificate mode with a smaller trusted surface) are done
  (`ZM_MIXED.md` §8), and the component and rejection tests were rerun on the final code in certificate mode.
* **`zmx2`**: self-tests (`search/zmx2_tests.sh`, 45/45: malformed and overflowing input, holes refused
  at the right pose, an exact-rational harness at 27,301 random poses, agreement with `zmcheck` on point covers, the
  audit's differential test) and a separate adversarial audit, `search/ZMX2_AUDIT.md`: no defect touching this
  certificate (≈ 3.5 M adversarial boxes, 117 sharp differential runs, all sound); **one must-fix, in the parser
  only** (unchecked `i128` overflow on crafted input integers could print a false side or total next to a
  verdict), with no effect on this cover.  Fixed by bounding the inputs (`ZMX2.md` §11), and both shipped `zmx2`
  runs were repeated with the patched binary.  It also re-verifies the shipped `s(13)` and `s(32)` point covers.
* Agreement: the two checkers reach the same verdict on this cover, and both refuse it scaled by `0.994`, at the
  same tilted dip `(0.5736, 1.4406, θ = 9.21°)` where the exact mass is `0.99969 < 1`.
* Lean: `notes/lean-s21.md`; the one hypothesis is the checkers' region statement, stated with the closed-square,
  parametric-fraction semantics of `FORMAT.md`.
