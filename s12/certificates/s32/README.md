# s(32) = 6

**Claim.**  No 32 unit squares fit in a square of side less than 6.  Since 36 unit squares tile the 6 × 6
square, **s(32) = 6**.  As far as we know this is the first exact value of `s(k² − 4)` for any `k ≥ 4`
(s(12), s(21), s(32), s(45), … were all open; `notes/literature-s32.md`).  The previous best lower bound was
`s(32) ≥ 119/20 = 5.95` (wand125, rectangle-density certificates, 2026-09-26); before 2026, Nagamochi's
general bound `1 + √23 = 5.7958` (its published proof is incomplete: Karakuş, arXiv:2609.37410).

Computer-assisted, not peer reviewed; checked by two exact checkers and, since 2026-09-28, kernel-checked in full in
Lean 4 with no hypothesis (`s32_eq_6`, `lean/Sqpack/S32Lower.lean`).  The proof is designed to be re-checked: `./verify.sh` here.

## The proof in one paragraph

`s32_closed_cover_6.txt` is a weighted point set in `[0,6]²` of total weight
`3171350535386 / 10¹¹ = 31.713505354 < 32` such that **every closed unit square contained in `[0,6]²`, at any
position and any angle, contains points of total weight at least 1** (a point on the square's boundary counts).
If 32 unit squares with disjoint interiors fit in a square of side `s < 6`, scale their centres by `6/s > 1`:
they become 32 pairwise *disjoint* closed unit squares in `[0,6]²`, so no point is in two of them, and together
they would capture weight `≥ 32 > 31.71`.  Contradiction.  The cover is invariant under the symmetries of the
square, so it suffices to check squares whose centre lies in `[0,3]²` and whose angle is in `[0°, 45°]`; that
check is an exhaustive, exact subdivision of pose space by a computer.

## What is checked by what

| step | how | trust |
|---|---|---|
| cover ⇒ `s(32) ≥ 6` (scaling, disjointness), `s(32) ≤ 6` (grid), D4 reduction, invariance and total of *this* cover | **Lean 4 + Mathlib**, `lean/Sqpack/{S32,D4,Cover,S32Data}.lean`: `s32_eq_six_of_checker : S32CheckerCover → minSide 32 = 6`; standard axioms only, no `sorry`, no `native_decide` (`notes/lean-s32.md`) | kernel |
| `S32CheckerCover`: every admissible closed unit square with centre in `[0,3]²`, `θ = 2 arctan u`, `u ∈ [0, ½]` captures weight `≥ 1` | **`search/zeromargin.py`**: all 7,200 root boxes certified, 0 uncertified (`zeromargin_d4/`) | exact rational/integer arithmetic; the program is not formally verified |
| the same, separately | **`verify2/` (`zmcheck`, Rust)**, a separately written checker (no shared code; written from `zeromargin.py`'s lemmas, so the point-test formulation is common): 3,595 / 3,600 roots of the same region (`zmcheck_d4/SUMMARY_candidate.txt`); its 5 open roots, at interior tile germs, are closed by `zeromargin.py` | exact integer arithmetic |
| the same, a third time (added 2026-09-27) | **`zmx2`** (`verify2/src/bin/zmx2.rs`, Rust), the mixed-cover checker written for `s(21)` (`../s21/`), which handles points as atoms of the grid lines: `zmx2 cert s32_closed_cover_6.txt --d4 --pair-points`, all 3,600 roots certified, 1,405,342 boxes, 0 uncertified, 873 CPU-s (`search/ZMX2.md` §9; repeated with the binary patched after its audit, `ZMX2.md` §11, same census; not shipped as a record here) | exact integers for points and masses; interval geometry (`../s21/README.md`) |
| the whole pose space, **no symmetry argument** (added 2026-09-30) | **`zmx2`** with the opt-in `--sym-atoms` (`search/ZMX2.md` §4.9 Lemma A, §12; source `zmx2.rs` sha256 `92a4cfe8…`, binary `ed31d3ee…`): `zmx2 cert s32_closed_cover_6.txt --full --pair-points --sym-atoms`, all 28,800 roots of `[0,6]²` × `u ∈ [0,½]` × both passes certified, 11,268,760 boxes, 0 uncertified, 11,592 CPU-s (`zmx2_full_sym/`).  Without `--sym-atoms` the full sweep leaves 152 boxes at the ±90° images of one tile germ (found by jlevy/squares): a point on two candidate lines was always made an atom of the vertical one; Lemma A bounds each box with that assignment and its mirror and keeps the larger.  So this checker no longer relies on the D4 fold | as the row above |
| `S32CheckerCover` itself, inside Lean's kernel (added 2026-09-28) | **`lean/Sqpack/S32Lower.lean`**: a generated tree of pose boxes (5,990 chunk theorems) decided by `decide +kernel` against the zero-margin verifier `ZMTree.check`, proved sound once (`ZMTree.sound`); gives `s32_checkerCover : S32CheckerCover` and `s32_eq_6 : minSide 32 = 6` with **no hypothesis**.  The tree generator is untrusted.  Opt-in: `lean/scripts/gen_data.sh S32Z`, then `lean/scripts/build_parts.sh S32Z Sqpack.S32Lower 4` (13.8 CPU-h; `lean/LADDER.md`) | kernel |
| that the Lean data is this file | `lean/scripts/gen_s32_data.py --check` (byte-for-byte), sha256 in the Lean file header | script |

A second cover, **`s32_shift_v1.txt`** (total `31.6979…`, `search/S32_SHIFT.md`), is certified in full by *both*
checkers: `zeromargin.py` 7,200 / 7,200 (`zeromargin_d4_shift_v1/`) and `zmcheck --d4` 3,600 / 3,600 with the
opt-in branch order `ZM_MIXPAIR=1` (`zmcheck_d4/SUMMARY_shift_v1.txt`).  Either cover proves the theorem; the
first is the one transcribed into Lean.

## Files

| file | |
|---|---|
| `s32_closed_cover_6.txt` | the cover: `certificates/FORMAT.md` plain format, `m = 6`, `D = 1000`, `W = 10¹¹`, 13,085 points, all integers |
| `s32_shift_v1.txt` | the second cover (same format) |
| `zeromargin_d4/` | the `zeromargin.py` run on the first cover: `manifest.json` (hashes, settings), `roots.jsonl` (one record per root: census, uncertified boxes, CPU, checker and certificate sha256), `SUMMARY.txt`, and `checker/zeromargin.py`, the exact checker file that ran |
| `zeromargin_d4_shift_v1/` | the same for `s32_shift_v1.txt` |
| `zmx2_full_sym/` | the `zmx2 --full --pair-points --sym-atoms` run: `roots.log.xz` (one line per root), `run.out`, `manifest.txt` (git HEAD, source/binary/input sha256, command, census, verdict) |
| `zmcheck_d4/` | `zmcheck --d4` sweep summaries (per-root logs are not shipped; rerun with `search/s32_sweep.py`) |
| `verify.sh` | re-check; below |
| `SHA256SUMS` | |

## Re-checking

```
certificates/s32/verify.sh          # seconds: hashes; Lean data = this file; both shipped runs re-summarised
                                    #   (7,200 roots each, hashes, exact D4 invariance); two germ roots re-run
certificates/s32/verify.sh --full   # re-runs zeromargin.py over the whole region for both covers
                                    #   (~2.8 + 2.4 CPU-h; about 12 min on 28 processes)
cd lean && lake build && lake env lean Axioms.lean                         # the Lean side
```

`zeromargin.py` settings: centre cells of pitch `1/10` over `[0,3]²`, 8 bins of `u = tan(θ/2)` of width `1/16`
over `[0, ½]` (i.e. `θ ≤ 53.13°`, more than the `45°` needed), subdivision to depth 24; one root (the wall tile
germ `x ∈ [1.5,1.6], y ∈ [0.5,0.6], u ∈ [0,1/16]`) needs depth 30 and is certified by a complete second run at that
depth.  Totals: 164,130 boxes, max depth 27, leaves ADM 12,201 / CHAIN 70,007 / EMPTY 3,457 / UNCERTIFIED 0,
2.77 CPU-h.  The checker's floats only pre-screen: every acceptance a certificate rests on is an exact integer or
`Fraction` test (`search/S32_EXACT.md` §11.3), and `--selfcheck` re-runs the reference `Fraction` path beside the
fast one at every call.

## What is not machine-verified

* The checker programs (subdivision, bookkeeping, parsing).  Their mathematical primitives are proved in Lean
  (`lean/Sqpack/ZeroMargin.lean`, `notes/lean-zeromargin.md`), the programs are not.  The two checkers share no
  code.  Since 2026-09-28 the result no longer rests on them: the covering statement is also checked inside Lean's
  kernel (above), so what remains trusted for `s32_eq_6` is Lean's kernel, Mathlib and the definitions in the
  statement.
* The run records: `verify.sh` re-summarises them, but the per-root censuses are the checker's own reports.
  `--full` regenerates them.

## Review record (2026-09-26)

Fresh-eyes audits before publication, none finding anything that affects the result: the Lean statement is the
standard `s(32)` (closed unit squares, arbitrary angles, disjoint open interiors, infimum attained) and its one
hypothesis matches the checker's conventions; every acceptance in `zeromargin.py`'s fast path is exact;
angles above 45° are handled soundly; `--selfcheck` and `--ref` reproduce the shipped censuses; and three
deliberately broken D4-symmetric covers (a corner orbit removed, a germ orbit removed, orbits removed so that
a pose at 52° falls below 1) are each refused.  History and all pilots: `search/S32_COVER.md`,
`search/S32_EXACT.md`, `search/S32_SHIFT.md`.
