# s(k² − 4) = k for every k ≥ 5

**Status.**  Working in public.  For `k ≥ 8`: an exact certificate of the finite statement
`Valid9`, made by the same checker code as the `k² − 3` certificate ([`../k2m3/`](../k2m3/README.md): six agents
reviewed that code adversarially, and wand125's independent checker corroborated its claim `Valid7`); in Lean, the
all-k reduction, the D4 reduction and the axis-parallel face (Lemma Z) are kernel-checked, so `s(k² − 4) = k` for all
`k ≥ 8` rests on one finite statement, `ValidTilt9` (the tilted squares of the run's region), certified by our run and,
since 2026-10-06, re-checked by wand125's separately written exact checker
([`valid7-independent-check`](https://github.com/wand125/valid7-independent-check/tree/c561dbb3fcea9178599db017f51f56b42ad9f91f) `c561dbb3`, release `records-tilt9-v1`; see below).  jlevy/squares has reviewed that
check and records it as reported evidence; T-081 is not yet register-verified.  Not yet externally reviewed or fully
formalised.  The fast `verify.sh` does not recompute the positive-tilt
part of the proof; `--full` does (below).  `k = 5, 6, 7` are separate bundles (below).

**Claim.**  For every integer `k ≥ 5`, no `k² − 4` unit squares fit in a square of side less than `k`.  Since `k²`
unit squares tile the `k × k` square, **`s(k² − 4) = k` for all `k ≥ 5`**: `s(21) = 5`, `s(32) = 6`, `s(45) = 7`,
`s(60) = 8`, `s(77) = 9`, … .  So the deficit-4 case of Friedman's Conjecture 1 [DS7] holds ("if `s(n² − 4) = n`
then `s((n+1)² − 4) = n + 1`"; we write `F₄`, [`../../search/FRIEDMAN.md`](../../search/FRIEDMAN.md) §0), whatever the
open `s(12)` is: at `n = 3` the hypothesis fails (`s(5) = 2 + 1/√2 < 3`), and at `n = 4` the conclusion `s(21) = 5`
holds.  The strength of the evidence differs by `k`; it is stated per `k` below.

* **`k ≥ 8`: this bundle** (the family below; Lean `bentz4_of_valid9` + `Valid9` by the certificate here).
* **`k = 5`: `s(21) = 5`**, [`../s21/`](../s21/README.md): two separately written exact checkers (shared point-test
  lineage); kernel-checked in Lean conditional on the checkers' region statement (`s21_eq_five_of_checker`).
* **`k = 6`: `s(32) = 6`**, [`../s32/`](../s32/README.md): two exact checkers, and kernel-checked in Lean in full,
  with no hypothesis (`s32_eq_6`).
* **`k = 7`: `s(45) = 7`**, [`../s45/`](../s45/README.md): two separately written exact checkers; not in Lean; no
  independent review yet.
* (`k = 8` also has its own bundle, `s(60) = 8`, [`../s60/`](../s60/README.md).)

**Earlier work.**  Friedman [DS7] conjectured that `s(n² − c) = n` implies `s((n+1)² − c) = n + 1`.  As far as we
know, no exact value of `s(k² − 4)` for `k ≥ 4` was known before 2026 (`s(12)` is still open; Nagamochi's general bound
is weaker, and its published proof is incomplete: Karakuş, arXiv:2609.37410; Kearney–Shiu and Bentz treat `k² − 3`).
Ours: `s(32) = 6` (2026-09-26), `s(21) = 5` and `s(45) = 7` (09-27), `s(60) = 8` (09-28).  wand125 proved `s(77) = 9`
first (2026-10-01, by inserting a band into our `s(60)` cover and checking it with our checkers; jlevy/squares#279).  The
cases `k ≥ 10`, and `k = 9` by a second route, are new here.  (`s(k² − 4) = k` implies `s(k² − 3) = k`, so this also re-derives the `k² − 3`
family for `k ≥ 5`.)

Computer-assisted, not peer reviewed.  The proof is designed to be re-checked: `./verify.sh` here.

## The proof in one paragraph

For each `k ≥ 7` there is a measure `μ_k` on `[0, k]²` (the family, `K4_k008_family.txt`: a corner module in each
corner, a periodic profile of period 1 along each wall in a band of width `w = 3`, Lebesgue measure on the middle
square `[14/5, k − 14/5]²`; all the non-Lebesgue mass sits on uniform pieces of the lines of the `1/5`-grid) with
total mass `k² − 4D`, `D = 214770225571/200000000000 = 1.073851… > 1`, such that (for `k ≥ 8`) **every closed unit
square in `[0, k]²`, at any position and angle, has mass at least 1**.  If `k² − 4` unit squares packed a square of
side `s < k`, scaling their centres by `k/s > 1` would give `k² − 4` pairwise disjoint closed unit squares in
`[0, k]²`, of total mass `≥ k² − 4 > k² − 4D`: a contradiction.  For `k ≥ 8 = 2R + 2` (edge-zone width `R = 3`) the
validity of `μ_k` reduces to that of `μ₉`: a unit square in `[0, k]²` has, under `μ_k`, the same mass as some
integer translate of it (one shift per axis) has under `μ₉`, because the profile is periodic and a unit square is
less than 2 wide.  So everything rests on one finite statement, **`Valid9`**: every closed unit square in `[0, 9]²`
has `μ₉`-mass `≥ 1`, where `μ₉` is the box cover `K4_k008_box9.txt` (2,076 segments on 80 grid lines plus the
Lebesgue square `[14/5, 31/5]²`, total `3835229774429/50000000000 = 76.7046… = 81 − 4D < 77`).  `Valid9` is the
certificate's claim: at `θ = 0` by an exact enumeration (Lemma Z), and for `θ > 0` by an exhaustive exact
subdivision of pose space (`qx2_zm.py`, run `qx2_k4x_k008`) over centres in `[0, 9/2]²` and `θ ∈ (0°, 45°]`, which
suffices because the cover is invariant under the symmetries of the square.  The margin is zero (the wall squares
and many tile-germ limits have mass exactly 1), so the checker's leaves include the exact zero-margin primitive
Lemma E, as for `k² − 3`.  The family was found by the `k² − 3` pipeline at `R = w = 3` with a tilt margin
(`κ = 0.08`, `θ₀ = 3°`) and exact projection ([`../../search/K2M4_MARGIN.md`](../../search/K2M4_MARGIN.md) §4).

## What is checked by what

| step | how | trust |
|---|---|---|
| `Valid9 ⇒ s(k² − 4) = k` for all `k ≥ 8`: the family `μ_k` for every `k`, `μ₉` = the box file, total `k² − 4D < k² − 4`, localisation to `μ₉`, dilation, the `k × k` grid | **Lean 4 + Mathlib**, `lean/Sqpack/{BentzFam,Bentz4,Bentz4Data,MixedMeasure}.lean`: `bentz4_of_valid9 : Valid9 → ∀ k, 8 ≤ k → minSide (k ^ 2 - 4) = k`, with `Valid9` stated on the box file verbatim (`box9Cover`); `famCover fam4 9 = box9Cover` by kernel evaluation (`box9Cover_measure`); `#print axioms`: `propext, Classical.choice, Quot.sound`; no `sorry`, no `native_decide` ([`../../notes/lean-k2m4-reduction.md`](../../notes/lean-k2m4-reduction.md)) | kernel |
| that the Lean data is these files | `lean/scripts/gen_bentzfam_data.py` regenerated from the bundled cover and family, compared byte for byte in `verify.sh` | script |
| the family: `σ = 0`, symmetries, `D`, it rebuilds the box file exactly, total `81 − 4D` | `search/qx2_family_check.py --k 9 --R 3 --w 3` (own code, reads only the two text files; its printed thresholds `4D > 3`, `< k² − 3` are the `k² − 3` ones); in Lean (`box9Cover_measure`, `famCover_total4`) | script; kernel |
| the box cover is well-formed, total `< 77` (i.e. `D > 1`), D4-invariant | `search/qx2_records.py cover` (own parser); `qx2_zm.py` also refuses a non-invariant cover | script |
| `Valid9` at `θ = 0` | **Lemma Z** (`search/QUADRANT_EXACT.md` §4.1): `qx2_zm.py axis`, exact enumeration of the 1,600 one-sided limit corners over centres `[½, 9/2]²` (D4): minimum exactly 1, 208 corners exactly tight (`qx2_zm/lemmaZ.out`) | exact; **and kernel-checked in Lean** (`validAxis9`, whole box, no symmetry) |
| `Valid9` for `θ > 0` on the D4 domain | **`qx2_zm.py` run `qx2_k4x_k008`** (`qx2_zm/`): 16,200 / 16,200 root boxes certified, 0 uncertified; leaves are exact primitives with paper proofs (`QUADRANT_EXACT.md` §3–4, `ZM_MIXED.md` §2) | exact `Fraction` arithmetic (floats only choose which exact test to try, and omit lines by a float distance test with ≈ 10⁻⁴ slack); the program is not formally verified |
| the D4 reduction (centres in `[0, 9/2]²`, `0 < θ ≤ 45°`, plus `θ = 0`, suffice) | **Lean**: `valid9_of_tilt_axis` (`lean/Sqpack/ValidSplit*.lean`), with the cover's D4 invariance checked by the kernel (`d4_box9`); the remaining hypothesis `ValidTilt9` is exactly the run's region (`notes/lean-valid-split.md`) | kernel-checked |
| the run record: settings, roots, every leaf, coverage | `search/qx2_records.py record` (below), independent of `qx2_zm.py`'s control flow.  Structure and coverage only; it does **not** recompute any leaf's mass bound | script |
| tiny-tilt germs (supplementary; not a step of the proof) | `search/qx2_germscan.py --hi 9/2`: 3,362 centres × 625 offsets at `u = 10⁻¹²`, min exact mass 1, none below (`qx2_zm/germscan.out`, shipped, not re-run) | exact, sampled |

**`qx2_records.py record`** (the same program as for `k² − 3`, generic in the box side) re-checks the record from
scratch: the header's sha256 are the files in `qx2_zm/checker/` and the box cover; the argv has the certificate
settings (`--depth 18 --exact-umax 1/2 --exact-from 3 --dump-leaves`, no region restriction); the `.out` repeats the
shas, the argv, the container `[0,9]²` and the root count; the roots are exactly the D4 grid `[0, 9/2]² × u ∈ [0, ½]`
(pitch `1/10`, 8 u-bins: 16,200), each once; per root, UNCERT 0, every leaf kind a certifying one, the census equal to
the leaf list, `boxes = 2·leaves − 1`, labels consistent (AXIS: `u₁ = 0`; SYM: `θ₀ ≥ 45°`; EMPTY: no admissible pose,
exactly); and **coverage**: the leaves have pairwise disjoint interiors, and every part of the root they leave
uncovered (the 3,043 slabs cut away by `clip_bin`) contains no admissible pose, decided exactly from
`w(u) = cos θ + sin θ = (1 + 2u − u²)/(1 + u²)` against the centre range.  (The record lists leaves only: the slabs, and the gaps left by zero-thickness EMPTY leaves — in 10 roots mid-way through the u-range — are not recorded, so the leaf tree cannot be rebuilt by midpoint splits alone; the record check and the `leaves` reviewer each reconstructed them and found every one inadmissible.)  So leaves + slabs tile each root (volume
`81/8` in total, exact).  Finally the census totals equal the `.out`.  Mutation tests on this record (a leaf dropped,
a leaf's bin shortened by 1/1000, a leaf relabelled EMPTY, a leaf duplicated, a region restriction in argv) are all
refused.  Like for `k² − 3`, a clean record is a complete, well-formed proof skeleton made by the shipped checker
files, not a fresh geometric check of `θ > 0` and not a second checker; only `verify.sh --full` recomputes the leaves,
with the same program.

## Files

| file | |
|---|---|
| `K4_k008_box9.txt` | the box cover `μ₉` (`k = 9`), mixed format v1 ([`../s21/FORMAT.md`](../s21/FORMAT.md)): `s = 9`, `D = 5`, `W = 2·10¹²`, 0 points, 2,076 axis-parallel segments on 80 lines of the `1/5`-grid (16 unit segments listed twice, once per piece: corner module and profile on the band's end lines), 1 polygon `[14/5, 31/5]²` with mass = area (Lebesgue).  Copy of `search/qx2_data/K4_k008_box9.txt` (sha `4151d7c4…`) |
| `K4_k008_family.txt` | the family (`R = w = 3`, pitch `1/5`, `σ = 0`): the profile and the corner module, exact rational masses.  Copy of `search/qx2_data/K4_k008_family.txt` |
| `K4_k008_exact.json` | the exact LP solution the two files were written from (`search/qx2_exact.py project`: 35,646 tight equalities, rank 13, 0 inconsistent, 0 negative); copy of `runs/qx2_k4x_k008/sol_exact.json` |
| `qx2_zm/checker/` | the exact files that ran the certificate (= the header's shas = `../k2m3/qx2_zm/checker/`): `qx2_zm.py` `6294052a…e737`, `zm_mixed.py` `1fd20346…ba95`, `zeromargin.py` `640fe453…86ab` (the file pinned by the `s(32)` bundle), `mixed_cover.py` `bb89de15…aae5` |
| `qx2_zm/run_k4x_k008_leaves.jsonl.gz` | the run's record (`runs/qx2_k4x_k008/qxzm_full.jsonl`, gzipped): a header (sha256 of the four checker files and the cover, argv), then one line per root: census, uncertified boxes (none), CPU, and **every leaf** (box after `clip_bin`, kind) |
| `qx2_zm/run_k4x_k008_leaves.out` | the run's log (`runs/qx2_k4x_k008/qxzm_full.out`) |
| `qx2_zm/lemmaZ.out` | the output of `qx2_zm.py axis K4_k008_box9.txt` (Lemma Z) |
| `qx2_zm/germscan.out` | the exact germ scan (`runs/qx2_k4x_k008/germscan.out`; supplementary) |
| `verify.sh`, `SHA256SUMS` | re-check; below |

## Re-checking

```
certificates/k2m4/verify.sh          # ~20 s: hashes; cover total < 77 and D4 invariance (own parser); family = box
                                     #   cover; Lemma Z re-run; the record re-checked (settings, roots, every leaf,
                                     #   coverage, census) -- structure only, no leaf's mass bound recomputed;
                                     #   Lean data regenerated and compared
certificates/k2m4/verify.sh --full   # also re-runs qx2_zm.py (the shipped checker/) with the run's settings
                                     #   (~815,000 CPU-s ≈ 226 CPU-h; about 16 h on 15 processes) and compares
                                     #   root for root
cd lean && lake build && lake env lean Axioms.lean                         # the Lean side
```

Run `qx2_k4x_k008` (2026-10-02/03): `qx2_zm.py sol_exact_box9.txt --depth 18 --nproc 15 --exact-umax 1/2
--exact-from 3 --dump-leaves`: root boxes of centre pitch `1/10` over `[0, 9/2]²` × 8 bins of `u = tan(θ/2)` of width
`1/16` over `[0, ½]` (`θ ≤ 53.13°`), 16,200 roots; subdivision to depth ≤ 18, Lemma E tried from depth 3.  Totals:
214,336 boxes, max depth 17, 115,268 leaves: PIECE 84,152 / EXACT 20,577 / EXACT0 1,424 / EXACT45 1,392 / LEB 683 /
CAP 393 / SYM 2,676 / AXIS 81 / EMPTY 3,890 / UNCERTIFIED 0; 815,343 CPU-s (56,808 s on 15 processes).  The leaf
kinds are those of `../k2m3/README.md` (PIECE: Lemmas S, T, L, R; LEB: U; CAP: K; EXACT*: E with E′, E″; SYM, AXIS,
EMPTY).  Before it, the Rust checker `zmx2` (first-order lemmas only) certified this cover except 4,844 boxes, all at
`θ < 0.05°` (`K2M4_MARGIN.md` §4); those tiny-tilt germ boxes are what needed Lemma E.

## What is not machine-verified

* `ValidTilt9` (the tilted part of `Valid9`).  Lean proves `ValidTilt9 → s(k² − 4) = k` for all `k ≥ 8`; `ValidTilt9` is the Python certificate's
  claim, made by our implementation (`qx2_zm.py` with the zm_mixed/zeromargin primitives), whose lemmas are proved on
  paper, not in Lean.  The code is byte-identical to the reviewed `k² − 3` checker, but this run's cover (box 9,
  `R = 3`, Lebesgue from `14/5`, corner pieces on the band's end lines) has not been reviewed.
  **A second implementation has checked `ValidTilt9`:** wand125's exact checker (the one that corroborated `Valid7`),
  [`valid7-independent-check`](https://github.com/wand125/valid7-independent-check/tree/c561dbb3fcea9178599db017f51f56b42ad9f91f) commit `c561dbb3`, records in release `records-tilt9-v1` (2026-10-06): centres
  `[0, 9/2]²`, `u ∈ [0, 7/16]`, no symmetry used, 28,350 roots, none uncertified.  jlevy/squares'
  [review](https://github.com/jlevy/squares/blob/main/docs/project/reviews/review-2026-10-06-wand125-validtilt9-independent-check.md) found that its region contains `ValidTilt9`'s and no defect, and records it as reported evidence on
  T-081, which stays V0/C1 (not register-verified).  We have not re-run its record check here.  (`zmx2` leaves 4,844
  tiny-tilt boxes open.)
* Kernel-checked in Lean since 2026-10-03 (`lean/Sqpack/ValidSplit*.lean`, `notes/lean-valid-split.md`): the D4 reduction (`valid9_of_tilt_axis`: the box cover is D4-invariant, and the tilted region the run covers — centres in `[0, m/2]²`, `0 < u`, `u² + 2u ≤ 1`, i.e. `0 < θ ≤ 45°` — plus the axis face gives every pose) and Lemma Z itself (`validAxis9`: every axis-parallel unit square in the box, by exact corner limits over the whole box).  So the only unformalised step is the tilted run (`ValidTilt9`).  Lean: `bentz4_of_validTilt9 : ValidTilt9 → ∀ k ≥ 8, minSide (k² − 4) = k`.
* Floats in the checker only choose which exact test to try, with the one exception noted for `k² − 3`
  (`lines_in_reach`, slack `≈ 9·10⁻⁵`).
* `k = 5, 6, 7`: as stated in their bundles (above); only `s(32) = 6` is hypothesis-free in Lean.

## Review record

Three adversarial review agents (2026-10-03), each briefed as a hostile referee on what is new relative to the reviewed
`k² − 3` certificate, each with its own code; reports kept in the project's private notes.  **No BREAKS, no GAP in the
mathematics; minor items fixed in this README and in `search/K2M4_MARGIN.md`.**

* **break-it** (own exact and float evaluators; none of the shipped checker code run).  Exact enumeration of all
  111,392 tiny-tilt one-sided limits on the 1/5-grid (both tilt signs, wall constraint): minimum exactly 1, none below;
  2.55 M float + 307 k exact poses at small tilts around the 10,824 limits equal to 1, aimed at the doubled end lines
  `3, 6`, the corner/band junction and the Lebesgue edges `14/5, 31/5`; 5.6 M multi-scale poses near `θ = 0`; ~6 M
  random/grid poses and ~1,500 local minimisations.  No pose below 1 (every float value below 1 was re-checked exactly).
  Notes: squares inside the Lebesgue block touching no segment have mass exactly 1 at every angle (zero margin by
  design); the smallest local minimum off that set is ≈ 1.0000448 (near `(4.904, 3.407)`, 44.5°); doubled segments are
  summed correctly by the checker.  Evidence, not proof, at general angles.
* **leaves** (the run record).  Own reconstruction of coverage: the roots are exactly the 45 × 45 × 8 grid, each root
  exactly covered by leaves plus gaps with no admissible pose; D4 + `θ ↦ π/2 − θ` reach every pose; Lemma E's hypotheses
  hold at box 9.  All 683 LEB and 393 CAP leaves re-proved with own code; own branch-and-bound proved most of a sample
  of PIECE / EXACT / EXACT45 leaves (the rest ran out of budget, nothing suspicious); the checker's primitives re-run on
  277 recorded leaves (bounds shown tight, not vacuous); ~0.9 M exact masses over every leaf kind: none below 1.  Own
  Lemma Z over the whole face: minimum exactly 1, 832 tight.  Minor: the record omits the clip_bin slabs and EMPTY gaps
  (documented above); the CPU-location wording in K2M4_MARGIN.md (fixed).  Not checked: Lemma E internals beyond the
  `k² − 3` review, independent proof of most EXACT leaves, a full root re-run.
* **claim** (statement chain and literature).  Lean `Valid9`/`bentz4_of_valid9` match the box file and FORMAT.md
  semantics exactly (closed squares, parametric segment fraction, density-1 polygon; 2076 segments in order, 16 doubled,
  D4-invariant); `k = 5, 6, 7` match the s21/s32/s45 bundles; no exact `s(k² − 4)`, `k ≥ 4`, was known before 2026.
  Gap (formalisation, not truth): the D4 reduction to the run's region was argued only in a docstring — **since closed
  in Lean** (`valid9_of_tilt_axis`, `validAxis9`).  Minor (fixed): "F₄" is our notation for the deficit-4 case of
  Friedman's Conjecture 1; state trust per `k`; credit wand125's earlier `s(77) = 9`.

The checker code itself (`qx2_zm.py` 6294052a…, which differs from the reviewed cdade4b6… only by an assert and a dump
flag) was reviewed by six agents for the `k² − 3` certificate ([`../k2m3/`](../k2m3/README.md)).  Not yet externally
reviewed.
