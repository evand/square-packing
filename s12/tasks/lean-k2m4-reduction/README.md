# Lean: s(k² − 4) = k for all k ≥ 8, conditional on Valid9 (2026-10-03)

Goal: the k² − 4 analogue of `SquarePacking.Bentz.bentz_of_valid7` (`lean/Sqpack/Bentz.lean`, notes
`notes/lean-bentz-reduction.md`):

    theorem bentz4_of_valid9 (h : Valid9) : ∀ k : ℕ, 8 ≤ k → minSide (k ^ 2 - 4) = k

(check the threshold: the family's box saving holds for k ≥ 2R + 2 = 8; if the localisation argument needs more, say so and
state what you proved).  `Valid9` = every closed unit square in [0,9]² has mass ≥ 1 under the box-9 cover taken verbatim
as a `MixedCover`, exactly as `Valid7` is for k² − 3.

Inputs (certificate of record, run 2026-10-02/03, VERIFIED-D4 by qx2_zm, 0 uncertified; axis + germ scan exact):
* `runs/qx2_k4x_k008/sol_exact_box9.txt`  (box k = 9, 2076 segments on the 1/5 grid + Lebesgue [14/5, 31/5]², total 81 − 4D)
* `runs/qx2_k4x_k008/sol_exact_family.txt` (R = w = 3, pitch 1/5, D = 214770225571/200000000000 ≈ 1.07385, σ = 0)
`runs/` is gitignored: copy both verbatim into `search/qx2_data/` as `K4_k008_box9.txt`, `K4_k008_family.txt` (record
sha256s; same role as `L4_k02_*`).

Approach: generalise `lean/scripts/gen_bentz_data.py` and `Bentz.lean` (codes: offsets 0..5R−1 from the nearer wall, band
phases, outside; Lebesgue corner a = 14/5; localisation shifts for box 2R + 3).  Prefer a version parametric in R if that's
cheap and keeps `bentz_of_valid7` intact (k² − 5 at R = 5 is next); otherwise a sibling `Bentz4.lean` + `Bentz4Data.lean` is
fine.  The kernel must cross-check the family table against the box file (as `box7Cover_measure` does) and prove the
accounting total = k² − 4D for every k, and k² − 4D < k² − 4.

Rules: kernel reduction only (no native_decide, no new axioms, no sorry); `#print axioms` = [propext, Classical.choice,
Quot.sound].  `bentz_of_valid7` and the default build must still pass.  Pin builds to physical cores (taskset -c 8-15 or
fewer; lake ignores taskset for module parallelism, use LEAN_NUM_THREADS / build_parts.sh as LADDER.md says); ≤ 8 cores.
Commit on your branch (commit message with the Co-Authored-By line from the session); don't push, don't merge.
Deliverables: the theorem in the default build if it checks in a few minutes (else opt-in, documented); a section in
`lean/LADDER.md` and a note `notes/lean-k2m4-reduction.md` (what's proved, the exact hypothesis, the data sha256s, timings);
report anything about the family that differs from k² − 3 (e.g. threshold, codes, band-zone ends).
