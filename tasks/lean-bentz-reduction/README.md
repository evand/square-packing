# lean-bentz-reduction: the conditional all-k theorem for s(k²−3) = k (2026-09-29)

**Goal.**  A kernel-checked (no `sorry`, no new axioms; end theorems print `[propext, Classical.choice, Quot.sound]`)
Lean theorem of the form

    theorem bentz_of_valid7 (h : Valid7) : ∀ k : ℕ, 6 ≤ k → minSide (k^2 - 3) = k

where `Valid7` says: every closed unit square `sq c θ 1 ⊆ [0,7]²` has mass ≥ 1 under **the exact measure μ₇ of the box
cover `search/qx2_data/L4_k02_box7.txt`** (800 segments + the Lebesgue polygon `[9/5, 26/5]²` with mass = area), built
as a `MixedCover` (`Sqpack/MixedMeasure.lean`).  So the whole claim reduces to one finite statement that the existing
Python certificate asserts and a future Lean certificate can discharge.  This is the hand-proved part of the argument:
`search/QUADRANT_EXACT.md` §2 (accounting, D4 invariance, localisation, one box suffices with k = 2R + 3 = 7, dilation)
and `search/QUADRANT.md` §1, §8.1.  Background on the family: QUADRANT_EXACT.md §1, the family file
`qx2_data/L4_k02_family.txt`, `qx2_exact.py` / `qx2_family_check.py` (how μ_k is built from the family).

**Suggested shape** (your call; say what you chose and why):
1. A generic structure for fixed-profile families (period-1 profile π on the band `[0,w]`, mirror-symmetric; corner
   module ν on `[0,R]²`, diagonal-symmetric; Lebesgue on the rest; σ = 0) and the measure `μ_k` on `[0,k]²`.
   Generic in the data where it is cheap: the same argument is wanted later for k²−4 (R = w = 3).
2. Accounting: `μ_k([0,k]²) = k² − 4D` for all k ≥ 2R + 2 (a sum over k − 2R periods; exact rationals).
3. D4 invariance of μ_k; localisation; the integer-shift argument (μ₇ valid ⇒ μ_k valid for every k ≥ 6), with
   care at the closed ends `x = R`, `x = k − R`, the layer, the seams.  Closed squares throughout.
4. Instantiate with the data: μ₇ from the family equals the box cover's MixedCover measure (or state Valid7 on the
   family's μ₇ and prove equality with the box-file cover — the hypothesis must be the thing the certificate proves);
   `D = 423621306389/500000000000 > 3/4`, σ = 0, symmetries, nonnegativity by `decide`/`norm_num`.  Generate Lean data
   from the files with a script in `lean/scripts/` (exact rationals; record input sha256s).
5. Dilation + upper bound ⇒ `minSide (k²−3) = k`: reuse `not_packs_of_measure`, `le_minSide_*`, and whatever the s(13)
   / s(32) end theorems use for the upper bound (tiling).  k = 3..5 are known results, not needed in Lean (optional).

If a step of the paper argument turns out wrong or incomplete, **stop and report it** (that is a finding, not an
obstacle to work around).  If the full generic version is slow, a version specialised to this family is fine.

**Rules.**  Work in your git worktree of `public/` (you were given one); new files only under `lean/` plus a note
`notes/lean-bentz-reduction.md`; don't edit existing Lean files unless a small additive lemma is needed (say so);
keep the default build green; commit on your worktree branch (no push, no merge to main).  Don't touch `search/`,
certificates, site.  **Compute:** physical cores 12–15 only (lake ignores taskset: use `LEAN_NUM_THREADS` /
`lake env lean` under `taskset -c 12-15`, see `lean/LADDER.md`, `scripts/build_parts.sh`), ≤ 32 GB.  Mathlib is in the
lake cache; don't rebuild it.  Report: the end theorem's exact statement, axioms printout, what is generic vs
specialised, build time, and any deviation from the paper argument.
