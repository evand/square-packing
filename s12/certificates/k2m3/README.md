# s(k² − 3) = k for every k ≥ 6

**Status.**  Working in public: a single-implementation exact certificate, adversarially reviewed by six independent
agents with no errors found; the all-k reduction is kernel-checked in Lean, conditional on the finite statement
`Valid7` (not proved in Lean).  Not yet independently re-implemented (a second implementation, `zmx2` with area
density, is in progress and does not yet cover the smallest tilts: `search/ZMX2_AREA.md`), externally reviewed, or
fully formalised.  The fast `verify.sh` does not recompute the positive-tilt part of the proof; `--full` does (below).

**Claim.**  For every integer `k ≥ 6`, no `k² − 3` unit squares fit in a square of side less than `k`.  Since `k²`
unit squares tile the `k × k` square, **`s(k² − 3) = k` for all `k ≥ 6`**: `s(33) = 6`, `s(46) = 7`, `s(61) = 8`,
`s(78) = 9`, … .

**Earlier work.**  Friedman [DS7, 1998] conjectured that `s(n² − k) = n` implies `s((n+1)² − k) = n + 1`, which
together with `s(6) = 3` [Kearney–Shiu 2002] would give `s(k² − 3) = k` for all `k ≥ 3`.  Bentz, having proved
`s(13) = 4` and `s(46) = 7` [Bentz 2010], stated that `s(m² − 3) = m` should hold for all `m ≥ 3` (the statement is
often called Bentz's conjecture) and proved the cases `m = 5, 6` in an arXiv preprint [Bentz 2016].  The cases
`k ≥ 8` were open (apart from `k = 8`, which this project obtained on 2026-09-28 as a corollary of `s(60) = 8`,
[`../s60/`](../s60/README.md)); the result here covers all `k ≥ 6`.  Since `s` is non-decreasing, `s(k² − 4) = k`
implies `s(k² − 3) = k`, so this project's `s(21) = 5`, `s(32) = 6`, `s(45) = 7` ([`../s21/`](../s21/README.md),
[`../s32/`](../s32/README.md), [`../s45/`](../s45/README.md)) also give `k = 5, 6, 7`, and the case-free `s(13) = 4`
([`../rung2/`](../rung2/README.md)) gives `k = 4`: no case `k ≥ 4` depends on Bentz's proofs, though his came first.
Only `k = 3` is taken from the literature [Kearney–Shiu 2002].  (For `k = 2` the statement is false: `s(1) = 1`.)

References: Kearney & Shiu, Electron. J. Combin. 9 (2002) #R14; Bentz, "Optimal packings of 13 and 46 unit squares in
a square", Electron. J. Combin. 17 (2010) #R126; Bentz, "Optimal packings of 22 and 33 unit squares in a square",
arXiv:1606.03746 (2016, preprint); Friedman, "Packing unit squares in squares: a survey and new results", Electron. J.
Combin. Dynamic Survey DS7 (1998, latest 2009); Nagamochi, Electron. J. Combin. 12 (2005) #R37 (states `s(k² − 2) =
s(k² − 1) = k`; its Lemma 1 is false, so the published proof is incomplete: chelokot's Lean counterexample and
replacement proof of `s(k² − 2) = k`, https://github.com/chelokot/square-packing-archive ; Karakuş, arXiv:2609.37410,
which re-proves `s(k² − 1) = k`).

**Corollary: `s(k² − 2) = k` for every `k ≥ 2`, without Nagamochi's Lemma 1.**  `s` is non-decreasing, so
`s(k² − 3) ≤ s(k² − 2) ≤ s(k² − 1) ≤ k`.  Hence `k ≥ 6` from this family (for `k ≥ 9` resting on this single-checker
certificate alone among our results; `k = 6, 7, 8` also from `s(32)`, `s(45)`, `s(60)`), `k = 5` from `s(21) = 5`,
`k = 4` from `s(13) = 4`, `k = 3` from `s(6) = 3` [Kearney–Shiu 2002] (and El Moumni 1999 directly), `k = 2` from
`s(2) = 2` [Göbel 1979].  None of these certificates uses Nagamochi's Lemma 1 or rectangle bound.  chelokot's
kernel-checked proof is an earlier, separate route for all `k ≥ 2`.

Computer-assisted, not peer reviewed.  The proof is designed to be re-checked: `./verify.sh` here.

## The proof in one paragraph

For each `k ≥ 6` there is a measure `μ_k` on `[0, k]²` (the family, `L4_k02_family.txt`: a corner module in each
corner, a periodic profile of period 1 along each wall, Lebesgue measure on the middle square `[9/5, k − 9/5]²`; all
the non-Lebesgue mass sits on uniform pieces of the lines of the `1/5`-grid) with total mass `k² − 4D`,
`D = 423621306389/500000000000 = 0.84724… > 3/4`, such that **every closed unit square in `[0, k]²`, at any position
and angle, has mass at least 1**.  If `k² − 3` unit squares packed a square of side `s < k`, scaling their centres by
`k/s > 1` would give `k² − 3` pairwise disjoint closed unit squares in `[0, k]²`, of total mass
`≥ k² − 3 > k² − 4D`: a contradiction.  The validity of every `μ_k` reduces to that of one of them: a unit square in
`[0, k]²` has, under `μ_k`, the same mass as some integer translate of it (one shift per axis) has under `μ₇`,
because the profile is periodic and a unit square is less than 2 wide.  So everything rests on one finite statement,
**`Valid7`**: every closed unit square in `[0, 7]²` has `μ₇`-mass `≥ 1`, where `μ₇` is the box cover
`L4_k02_box7.txt` (800 segments on 60 grid lines plus the Lebesgue square `[9/5, 26/5]²`, total
`5701378693611/125000000000 = 45.611… = 49 − 4D < 46`).  `Valid7` is the certificate's claim: at `θ = 0` by an exact
enumeration (Lemma Z), and for `θ > 0` by an exhaustive exact subdivision of pose space (`qx2_zm.py`, run V3) over
centres in `[0, 7/2]²` and `θ ∈ (0°, 45°]`, which suffices because the cover is invariant under the symmetries of
the square.  The margin is zero: the wall squares `[t, t+1] × [0, 1]` and many tile-germ limits have mass exactly 1,
so the checker's leaves include an exact zero-margin primitive (Lemma E) that bounds the mass at the vertices of an
arrangement as rational functions of `u = tan(θ/2)`.

## What is checked by what

| step | how | trust |
|---|---|---|
| `Valid7 ⇒ s(k² − 3) = k` for all `k ≥ 6`: the family `μ_k` for every `k`, `μ₇` = the box file, total `k² − 4D < k² − 3`, localisation to `μ₇`, dilation, the `k × k` grid | **Lean 4 + Mathlib**, `lean/Sqpack/{Bentz,BentzData,MixedMeasure}.lean`: `bentz_of_valid7 : Valid7 → ∀ k, 6 ≤ k → minSide (k ^ 2 - 3) = k`, with `Valid7` stated on the box file verbatim (`box7Cover`); `famCover 7 = box7Cover` checked by kernel evaluation (`box7Cover_measure`); `#print axioms`: `propext, Classical.choice, Quot.sound`; no `sorry`, no `native_decide` ([`../../notes/lean-bentz-reduction.md`](../../notes/lean-bentz-reduction.md)) | kernel |
| that the Lean data is these files | `lean/scripts/gen_bentz_data.py` regenerated from the bundled cover and family, compared byte for byte in `verify.sh` | script |
| the family: `σ = 0`, symmetries, `D`, it rebuilds the box file exactly, total `49 − 4D` | `search/qx2_family_check.py` (own code, reads only the two text files); and in Lean (`box7Cover_measure`, `famCover_total`) | script; kernel |
| the box cover is well-formed, total `< 46`, D4-invariant | `search/qx2_records.py cover` (own parser); `qx2_zm.py` also refuses a non-invariant cover | script |
| `Valid7` at `θ = 0` | **Lemma Z** (`search/QUADRANT_EXACT.md` §4.1): `qx2_zm.py axis`, exact enumeration of the 900 one-sided limit corners over centres `[½, 7/2]²` (D4): minimum exactly 1, 188 corners exactly tight (`qx2_zm/lemmaZ.out`) | exact `Fraction` arithmetic; paper lemma |
| `Valid7` for `θ > 0` on the D4 domain | **`qx2_zm.py` run V3** (`qx2_zm/`): 9,800 / 9,800 root boxes certified, 0 uncertified; leaves are exact primitives with paper proofs (`QUADRANT_EXACT.md` §3–4, `ZM_MIXED.md` §2) | exact `Fraction` arithmetic (floats only choose which exact test to try, and omit lines by a float distance test with ≈ 10⁻⁴ slack, below); the program is not formally verified |
| the D4 reduction (centres in `[0, 7/2]²`, `θ ∈ [0°, 45°]` suffice) | the cover's exact D4 invariance (above); the reduction itself on paper ([`../s21/FORMAT.md`](../s21/FORMAT.md)); the SYM leaves use it for `θ > 45°` | paper |
| the run record: settings, roots, every leaf, coverage | `search/qx2_records.py record` (below), independent of `qx2_zm.py`'s control flow.  It checks the record's structure and that the leaves cover every root; it does **not** recompute any leaf's mass bound | script |

**`qx2_records.py record`** re-checks the V3 record from scratch: the header's sha256 are the files in
`qx2_zm/checker/` and the box cover; the argv has the certificate settings; the roots are exactly the D4 grid
`[0, 7/2]² × u ∈ [0, ½]` (pitch `1/10`, 8 u-bins), each once; per root, UNCERT 0, every leaf kind a certifying one, the
census equal to the leaf list, `boxes = 2·leaves − 1`, labels consistent (AXIS: `u₁ = 0`; SYM: `θ₀ ≥ 45°`; EMPTY:
no admissible pose, exactly); and **coverage**: the leaves have pairwise disjoint interiors, and every part of the
root they leave uncovered (the slabs cut away by `clip_bin`) contains no admissible pose, decided exactly from
`w(u) = cos θ + sin θ = (1 + 2u − u²)/(1 + u²)` against the centre range.  So leaves + slabs tile each root (volume
`49/8` in total, exact).  Finally the census totals equal the `.out`.  Mutation tests (a leaf dropped, a leaf's bin
shortened by 1/1000, a leaf relabelled EMPTY, a leaf duplicated, a region restriction in argv) are all refused.

**What `qx2_records.py` does not do.**  It does not re-prove any leaf.  A `PIECE`, `EXACT`, `LEB` or `CAP` leaf is
accepted on its label, its interval and the coverage checks above; the inequality the label stands for (mass `≥ 1` on
that box) is not recomputed.  (EMPTY, AXIS and SYM labels are checked exactly, as above: no admissible pose; only `θ = 0`, which Lemma Z
covers; `θ ≥ 45°`.)  So the
fast `verify.sh` establishes that the shipped record is a complete, well-formed proof skeleton made by the shipped
checker files, plus Lemma Z (re-run), but it is not a fresh geometric check of `θ > 0` and not a second checker.  Only
`verify.sh --full`, which re-runs `qx2_zm.py`, recomputes the leaves; it uses the same program.

## Files

| file | |
|---|---|
| `L4_k02_box7.txt` | the box cover `μ₇` (`k = 7`), mixed format v1 ([`../s21/FORMAT.md`](../s21/FORMAT.md)): `s = 7`, `D = 5`, `W = 10¹²`, 0 points, 800 axis-parallel segments on 60 lines of the `1/5`-grid, 1 polygon `[9/5, 26/5]²` with mass = area (Lebesgue).  Copy of `search/qx2_data/L4_k02_box7.txt` |
| `L4_k02_family.txt` | the family (`R = w = 2`, pitch `1/5`, `σ = 0`): the profile (52 pieces per period + Lebesgue on `y ∈ [9/5, 2]`) and the corner module (36 pieces + Lebesgue on `[9/5, 2]²`), exact rational masses |
| `L4_k02_exact.json` | the exact LP solution the two files were written from (`search/qx2_exact.py`) |
| `qx2_zm/checker/` | the exact files that ran V3: `qx2_zm.py` `6294052a…e737`, `zm_mixed.py` `1fd20346…ba95`, `zeromargin.py` `640fe453…86ab` (the file pinned by the `s(32)` bundle), `mixed_cover.py` `bb89de15…aae5` |
| `qx2_zm/runV3_6294052a_leaves.jsonl.gz` | run V3's record: a header (sha256 of the four checker files and the cover, argv), then one line per root: census, uncertified boxes (none), CPU, and **every leaf** (box after `clip_bin`, kind) |
| `qx2_zm/runV3_6294052a_leaves.out` | run V3's log |
| `qx2_zm/lemmaZ.out` | the output of `qx2_zm.py axis L4_k02_box7.txt` (Lemma Z) |
| `verify.sh`, `SHA256SUMS` | re-check; below |

## Re-checking

```
certificates/k2m3/verify.sh          # ~3 s:  hashes; cover total and D4 invariance (own parser); family = box cover;
                                     #   Lemma Z re-run; the V3 record re-checked (settings, roots, every leaf,
                                     #   coverage, census) -- structure only, no leaf's mass bound recomputed;
                                     #   Lean data regenerated and compared
certificates/k2m3/verify.sh --full   # also re-runs qx2_zm.py (the shipped checker/) with the V3 settings
                                     #   (~81,000 CPU-s; about 3 h on 8 processes) and compares root for root
cd lean && lake build && lake env lean Axioms.lean                         # the Lean side
```

Run V3 (2026-09-29): `qx2_zm.py L4_k02_box7.txt --depth 18 --nproc 8 --exact-umax 1/2 --exact-from 3 --dump-leaves`:
root boxes of centre pitch `1/10` over `[0, 7/2]²` × 8 bins of `u = tan(θ/2)` of width `1/16` over `[0, ½]`
(`θ ≤ 53.13°`, more than the `45°` needed), 9,800 roots; subdivision to depth ≤ 18, Lemma E tried from depth 3.
Totals: 54,358 boxes, max depth 17, 32,079 leaves: PIECE 15,926 / EXACT 9,763 / EXACT0 381 / EXACT45 759 / LEB 689 /
CAP 374 / SYM 1,171 / AXIS 61 / EMPTY 2,955 / UNCERTIFIED 0; 81,377 CPU-s (10,232 s on 8 processes).  The leaf kinds:
PIECE = zm_mixed's piece bound (Lemmas S, T, L, R; the Lebesgue square through its polygon Lemma S(b)); LEB (Lemma U:
every admissible square inside the Lebesgue square); CAP (Lemma K); EXACT (Lemma E with E′, E″; EXACT0 at `u₀ = 0`,
EXACT45 certifies the `θ ≤ 45°` part of a box straddling 45°); SYM (`θ ≥ 45°` on the whole box: the diagonal image
is a pose of the domain with `θ ≤ 45°`); AXIS (only `θ = 0` poses left: Lemma Z); EMPTY (no admissible pose).
Proofs: `search/QUADRANT_EXACT.md` §3–4, `search/ZM_MIXED.md` §2.

## What is not machine-verified

* `Valid7` itself.  Lean proves `Valid7 → s(k² − 3) = k` for all `k ≥ 6`; `Valid7` is the Python certificate's claim.
  The certifying programs (`qx2_zm.py`, ~1,300 lines, with the zm_mixed/zeromargin primitives it imports) are one
  implementation; their lemmas are proved on paper, not in Lean.  There is no complete second, independent
  implementation of the `θ > 0` part yet.  The Rust checker `zmx2` of the other bundles has been extended to area
  density (`search/ZMX2_AREA.md`, 2026-09-30, written without opening `qx2_zm.py` or its write-ups): it certifies
  `θ = 0` and every `θ ≥ 0.014°` over the whole pose space with no symmetry assumed, but at `0 < θ < 0.01°` it leaves
  120 boxes of the D4 region (1,108 of the full space) uncertified at the double germs (ZMX2_AREA.md §11.4), so it
  does not yet certify `Valid7`.
* The D4 reduction to the fundamental domain and Lemma Z's reduction to finitely many corner limits (paper).
* Floats in the checker only choose which exact test to try, with one exception (review lemmaE-code, minor):
  `lines_in_reach` omits lines from Lemma E's arrangement by a float distance test (radius `0.7072 + ½·diagonal`; a
  unit square lies within `√2/2 = 0.70711` of its centre, a slack of `≈ 9·10⁻⁵`, far above float error).

## Review record (2026-09-29)

Six agents reviewed the certificate adversarially and independently (hostile-referee brief: derive the statements
before reading the proofs, design own exact tests; `Fraction` arithmetic for any claimed violation).  **No BREAKS
and no soundness GAP in any of the six**; one auditability GAP, now closed.

* **break-it** (empirical attack): own exact evaluator from FORMAT.md; exact scans of the `θ = 0` face (4,500 limit
  corners over the whole `[½, 13/2]²`, no symmetry used), of tile germs at `u = 10⁻¹²…10⁻¹⁵` and of the wall families:
  nothing below 1; the tight germ limits grow like `1 + 0.319|u|`.  The tightest pose off the known tight set is a
  diamond near `45°` (mass `1 + 1.03·10⁻⁴`).  Six D4-symmetric mutants: the one run through the box checker (a
  tile germ at `1 − 2·10⁻⁴`, invisible at `θ = 0`) was refused at exactly that germ.
* **lemmaE-code** (every line of `Exact` read): 449 boxes (random, at tight corners and germs, sub-boxes of the
  run-B root) on real and mutated covers, each asked to certify `τ` just above its exact minimum: never certified.
  Minor: an empty list of alternatives would have certified vacuously (unreachable) — **guard added** (`qx2_zm.py`
  `6294052a`: assert); the float `lines_in_reach` (above).
* **lemmaE-math** (Lemma E, E′, E″ derived independently; 20,000 exact poses for E″, 12,000 end-to-end trials,
  mutations of the lemma refused): minor wording — Lemma E's `f` is continuous on each closed cell of the arrangement
  (not on all of `R(u)` in the corner regimes, where the code takes the minimum of the two cell functions at vertices
  on `z = 0`, which is what the proof needs); the McCormick loss is `≤ h_x h_y + u h_y²/(4(1 − u²))`; neither affects
  soundness.
* **reduction-coverage**: the reduction `μ₇ ⇒ μ_k (k ≥ 6)` and the dilation derived independently and checked
  exactly for `k = 6..14`; roots tile the D4 domain; 128 roots re-run with leaf dumps: leaves + `clip_bin` slabs tile
  them.  **GAP (auditability): run V2's record had per-root totals only**, so coverage could not be audited without a
  re-run — **closed by run V3**, which records every leaf, and by `qx2_records.py`, which checks the coverage of all
  9,800 roots.  Minor: the small cases (attribution above; `k = 2` is false); the format file is
  `certificates/s21/FORMAT.md`.
* **zm-polygon** (zm_mixed's polygon Lemma S(b), unaudited before this certificate): re-derived; 188,408 exact poses
  on 1,686 boxes of this cover and 25,872 poses on 462 certificate PIECE leaves: the polygon bound never exceeds the
  true area.  Minor wording only.
* **z-u-k-leaves** (Lemmas Z, U, K, `clip_bin`, EMPTY/AXIS/SYM): own Lemma Z over the whole `[½, 13/2]²` (3,600
  limits, min exactly 1); 2.6 M exact `clip_bin`/EMPTY/AXIS tests; 4,763 roots re-run with leaf dumps, every AXIS and
  EMPTY leaf and most LEB/CAP leaves checked against exact masses; about 1.1 M exact tiny-tilt poses: none below 1.
  Minor: `qx2_zm.py axis` covers `[½, 7/2]²` only, relying on the D4 invariance it does not itself check (checked
  here by `qx2_records.py cover` and by `qx2_zm.py`'s main run).

**Guards and re-run.**  After the reviews: `zm_mixed.py` `1fd20346…` (was `ee3e2915…`: `region_phi` refuses a
zero-width angle bin instead of returning an unproved EMPTY; latent bug B1 of the lean-segments audit, unreachable
here) and `qx2_zm.py` `6294052a…` (was `cdade4b6…`: the empty-alternatives assert, `--dump-leaves`, and the final
message naming `qx2_zm.py axis`).  **Run V3** repeated the whole certificate with these files and every leaf
recorded: VERIFIED-D4, 0 uncertified, **census identical to run V2's root for root** (all 9,800 roots, every
count).  Earlier runs A, B, V2: `search/qx2_data/cert/`, `search/QUADRANT_EXACT.md` §6.  Reviews:
`private/s12/tasks/k2m3-review/` (not public).
