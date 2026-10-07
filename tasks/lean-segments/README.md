# lean-segments, batch 1: mixed covers (points + grid-line segments) in the zero-margin Lean verifier (2026-09-28)

**Why.**  s(21) = 5, s(45) = 7, s(60) = 8 are proved by `search/zm_mixed.py` (Python, exact but unverified).
`lean/LADDER.md` "Next: segments" scopes the hypothesis-free Lean rung.  This batch is the **mathematics**: the
soundness side and a generator validated on toys.  Not the big kernel builds (15–40 CPU-h each) — those come later.

**Read first.**  `lean/LADDER.md` (all; esp. rung 2 / `ZMTree` and "Next: segments"), `search/ZM_MIXED.md` §1–2, §7–8
(the lemmas: P, S, T, L/L′, R/SPLIT, PIECE, `rf_bound`), `lean/Sqpack/{ZMTree,MixedMeasure,SegTree,S21}.lean`,
`lean/scripts/gen_zmtree.py`.

**Targets, in order** (each independently useful; stop where the budget ends):
1. **`CovM`**: `Cov` for the mixed mass (points inside + Σ segment weight × inside fraction), `splitX/Y/U` gluing,
   end-to-end `le_minSide_mixed` (via `not_packs_of_measure`, `d4_reduction_measure_u`), D4 check for segments
   (`SegTree` exists).  A toy end-to-end theorem that exercises it.
2. **Lemma P + Lemma S** in the `Z` leaf: phantom point of weight `L` from a certified core of a segment
   (Bernstein coefficients affine in the point → rational interval on the line).  Plus the `PIECE` leaf.  Exact
   `Nat`-only tests, one soundness lemma each, `ZMTree.sound`-style theorem for the extended tree type (new file
   `ZMTreeM.lean` or equivalent; do not change `ZMTree.lean`'s existing theorems — s(13)/s(32) must still build).
3. **Lemma L / L′ and Lemma T** (chord-length lower bounds on an axis line; the germ threshold lemma).  This is the
   genuinely new real analysis; formalise **what `zm_mixed.py` tests**, and where the note and the code differ,
   record it (`notes/lean-segments.md`).  If a lemma resists, state it precisely in Lean, prove what you can, keep
   it out of the default build, and explain the obstacle.
4. **Lemma R / `SPLIT`.**
5. **Generator**: `gen_zmtree.py`'s analogue with `zm_mixed.py` as read-only oracle and an exact integer mirror of the
   new `check`.  Validate on the toy covers of `ZM_MIXED.md` §4.3 and, if cheap, a rejection test (a cover that must
   fail).  Report the leaf census and a per-leaf kernel cost measured on a small chunk, and extrapolate s(21).

**Rules.**  Kernel reduction only (`decide +kernel`; no `native_decide`, no new axioms); `#print axioms` on every new
end theorem = `[propext, Classical.choice, Quot.sound]` (`lean/Axioms.lean`).  No `sorry` in anything imported by
`Sqpack.lean`.  Keep existing builds green (`lake build` default target).  Nothing about trees/generator is trusted.

**Compute.**  Physical cores 12–15 only, ≤ 4 Lean processes, ≤ 35 GB (lake ignores taskset: use
`lean/scripts/build_parts.sh` or `taskset` on `lake env lean`).  Do **not** build S21/S45/S60 data sets or any
multi-hour kernel run.  Mathlib: symlink `/home/evand/math/square-packing/public/s12/lean/.lake` into the worktree.

**Deliverables.**  Lean files + `Sqpack.lean`/`Axioms.lean` updates, generator script, `notes/lean-segments.md`
(statement ↔ implementation correspondence, what is proved, what isn't), a LADDER.md section update; commits on
your worktree branch.  Don't edit `TODO.md`, `README.md`, `search/*.py`, `verify*/`, `certificates/`, `site/`.
**Report** every ~4–6 h of work at most: which targets compile, axioms, line counts, open obstacles.
