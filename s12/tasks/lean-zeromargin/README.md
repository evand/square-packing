# lean-zeromargin: Lean statements and proofs for the zero-margin checker's lemmas (2026-09-12)

**Why.**  `search/zeromargin.py` proves `s(13) = 4` from `certificates/rung2/s13_closed_cover_4.txt`
(`search/RUNG2.md` §0) and is the base of V1, the `t = 4` verifier the `s(12)` endgame needs.  Its
correctness rests on a handful of paper lemmas in `RUNG2.md`.  `lean/Sqpack/Basic.lean` already
formalises the reduction (weighted cover of total `W` ⟹ at most `W` squares) and `Chord.lean` the
chord lemma, 0 sorries, Mathlib `v4.33.1`.  This task adds the checker's lemmas in the same style.

**Targets, in order** (stop at whichever the budget allows; each is independently useful):
1. **Lemma A** (`RUNG2.md` §3.1, monotone corners): for a box `[cx₀,cx₁]×[cy₀,cy₁]×[u₀,u₁]` with
   `θ ∈ [0, π/2]`, if the four inequalities (i)–(iv) hold at the four *specified* corners for every
   `θ` in the bin, then `p ∈ sq c θ 1` for every admissible pose of the box.  State it with
   `SquarePacking.sq` / `coord` from `Basic.lean`.  Proof is monotonicity of two affine forms in
   `(c_x, c_y)` with `cos θ, sin θ ≥ 0`.
2. **Lemma E** (§6, the four-corner maximum): a function of the form
   `g(c_x, c_y, u) = a(u) c_x + b(u) c_y + d(u)` with `a, b, d` polynomials of degree ≤ 2 in `u`
   attains its maximum over a box at one of the four centre corners for each `u`, and the max over
   `u ∈ [u₀,u₁]` of a quadratic is at an endpoint or the interior vertex (Lemma C's `deg ≤ 2`
   branch).  State exactly what the checker uses.
3. **Lemma G** (§6, the region test): if `max_B (g_i + λ G) ≤ 0` for some `λ ≥ 0` then `g_i ≤ 0` on
   `B ∩ {G ≤ 0}`.  (One line; but it is what makes `CHAIN` sound, so state it.)
4. **`clip_bin`** (§4.6): the admissible set of a wall-touching box at angle `θ` is
   `c_x ∈ [max(cx₀, w(θ)/2), min(cx₁, 4 − w(θ)/2)]`, with `w(θ) = cos θ + sin θ`, i.e. the
   admissibility of a pose is equivalent to those bounds — `Basic.lean` may already have the
   support-function fact for `sq ⊆ [0,m]²`; reuse it.
5. **Lemma B/C** (§3.2–3.3, the polynomial forms in `u = tan(θ/2)` and the Bernstein bound) only
   if 1–4 are done; the Bernstein convex-hull bound may exist in Mathlib (`Polynomial.bernstein`).

Do **not** attempt to formalise the subdivision or the checker's control flow; the goal is that
each primitive's *soundness* is a Lean theorem whose hypotheses are exactly what the Python tests.
Where the note's statement and the Python's test differ, formalise what the Python tests and
record the difference in the write-up.

**Build.**  Mathlib is 7.5 GB and lives at `/home/evand/math/square-packing/s12/lean/.lake`; a
fresh worktree has no `.lake`.  Symlink it: `ln -s /home/evand/math/square-packing/s12/lean/.lake
<worktree>/lean/.lake` (nothing else builds Lean right now; the `Sqpack` build outputs are small).
`cd lean && lake build` must end with 0 errors, and `#print axioms` on each new theorem must show
only `propext`, `Classical.choice`, `Quot.sound` (see `lean/Axioms.lean`).  No `sorry` in
committed files; if a target is unfinished, leave it out of `Sqpack.lean` and say so.

**Budget.**  Lean builds are single-threaded per file; this task uses no other compute.  Report
within ~4 h even if not done, with which targets compile.

**Deliverables.**  `lean/Sqpack/ZeroMargin.lean` (or one file per lemma), `Sqpack.lean` updated,
`Axioms.lean` updated, `lean/README` or `search/RUNG2.md` §-pointer note (`notes/lean-zeromargin.md`
with the statement/implementation correspondence), commit on your worktree branch.  Do not edit
`TODO.md`, `README.md`, `search/*.py`, `verify/`, `certificates/`, or other tasks' files.
