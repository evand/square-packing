# I. Anchor cliques on the cover side: certificate object, checkers, Lean, and the first leaf

**Why.** Task G (`search/CLIQUE_CONTINUUM.md`, `notes/clique-family.md`) found a clique family that
is valid by a two-line argument and cannot be dodged by column generation:

    K(p, A) = { S : p ∈ S and S ∩ A ≠ ∅ } ∪ { S : A ⊆ S }        (Lemma 0: a clique for every p, A)

and measured, on the packing side, that it is worth **≈ 0.10–0.12** at `t = 3.99` and in the
`k = 4` corner leaf at `3.98`, stable across an 8.5× pose refinement, against a pure excess of
0.008–0.026 at 3.99 and a leaf value of `12.000 ± 0.001` at 3.98.  Packing-side numbers are lower
bounds on the clique-LP value; the proof that the value is `< 12` is a **certificate**, and that is
this task.  The cover-side design is already written: `notes/clique-family.md` §5 (format, the two
exact predicates `contains` / `meets`, the credit rule, the Lean lemma `clique_of_anchors`), and
`search/CLIQUE_CONTINUUM.md` §6.  Box cliques (`search/BOXCLIQUE.md`) are the precedent for every
piece of plumbing: parser, well-formedness refusal, sweep credit, witness placement outside the
region, `xcheck.py` mirror, rejection tests, `verify.sh`, `SHA256SUMS`.

Read first: `notes/clique-family.md` (all), `search/CLIQUE_CONTINUUM.md` (Verdict, §1, §5, §6, §7a),
`search/BOXCLIQUE.md`, `certificates/FORMAT.md`, `verify/src/main.rs` (clique block, sweep credit,
witness placement), `xcheck.py`, `tests/rejection_tests.sh`, `lean/Sqpack/Basic.lean`
(`packing_le_weight_cliques`, `clique_of_cores`), `search/boxclique.py`, `search/branch.py`
(restricted master `BRANCH_SOLVER=restricted`, `--dump-lp`, `--cols-raw`), `search/LPSPEED.md`
(Recommendation), `search/BRANCH.md` (the `k = 4` leaf, its dual `runs/branch_t398hk4_dual_it16.txt`),
`search/ZEROMARGIN.md` §2 (the exact bin core, which replaces the σ-square for `contains`).

## Deliverables, in order

1. **Format.** `certificates/FORMAT.md`: the anchor-clique block per `clique-family.md` §5 (`piece`,
   `anchorP`, `anchorS` lines; weight over `10^7`; counted once).  Keep box cliques as they are.
2. **Verifier** (`verify/src/main.rs`).  Well-formedness: for every ordered pair of pieces the
   anchors intersect (exact) or one filters the other — refuse otherwise (`ERROR`).  Two exact `i128`
   predicates per (cell, bin): `contains(anchor)` — every pose of the cell contains the anchor,
   using the exact bin core (point: core test; segment: both endpoints) — and `meets(anchor)` — every
   pose of the cell meets the anchor, by refuting the separating axes over the bin (conservative).
   A cell gets `w_K` iff it lies wholly in some piece; partial cells get nothing and their LP witness
   is placed outside the piece.  With no anchor blocks the output must be **bit-identical** to today
   (diff against the pre-change binary on the main certificate and the box-clique demo, single
   thread, witness mode included).  Full `[0°, 90°]` sweep whenever anchor cliques are present.
3. **`xcheck.py`**: the same in `Fraction`, anchors inflated (accepts everything the Rust accepts).
4. **Rejection tests**: anchors that neither meet nor filter; a piece whose filter list names itself
   or a missing piece; a segment anchor outside the container; weight credited to a straddling cell
   (must not be); zeroed anchor weights (points must still verify); a valid anchor clique reported as
   such for `n = 13`; malformed lines.  Existing tests must stay green.
5. **Lean**: `clique_of_anchors` as stated in `clique-family.md` §5 (three-way case split; no
   geometry), 0 sorries, axioms unchanged; connect it to `packing_le_weight_cliques` the way
   `clique_of_cores` is.  `lake build` must pass.
6. **LP driver**: anchor-clique columns inside the cutting-plane loop — extend `search/branch.py`
   (so the `k = 4` leaf and the mixed leaves can use them) and, if cheap, `search/tighten.py`.
   Pricing on the D4-averaged interior-point dual, as `boxclique.py` learned: a candidate is a wall
   point `p` (strictly within distance 1 of a wall — Lemma 1 says nowhere else helps) plus a segment
   `A` chosen by Lemma 2 so that `K(p, A) ⊇ P_p` (then the clique column *dominates* the point
   column and can replace it), priced by the dual mass of `{S : A ⊆ S}` beyond `P_p`; also the
   max-anchor-clique heuristics of `search/clique_family.py` / `clique_continuum.py` on the current
   dual's support.  Every clique written must pass the verifier's own well-formedness check.
7. **The first leaf.** `k = 4`, `r = 1`, `t = 3.98` with anchor-clique columns, restricted master,
   `N = 2000`, seeded from `runs/branch_t398hk4_cols.txt` / `_probe.txt` (main tree).  The pure
   loop sits at `12.000 ± 0.001`; the packing side says cliques are worth ≈ 0.10 here.  Goal: a
   **verified** leaf certificate `certificates/branch/s12_t3.98_corner_k4.txt` with value `< 12`.
   Then, if it closes, the remaining per-box leaf `1110` (restart from
   `runs/branch_t398j1110_cols_it51.txt` / `_probe_it51.txt`, `--cols-raw`) with cliques allowed —
   together with `k = 0, 1, 2` (shipped) that is **`s(12) ≥ 3.98`**.
8. If time remains: the pure cover at `t = 3.99` with anchor cliques (`tighten.py` path).  The
   pure value is 12.20 there; the packing side says the clique value on rich pose sets is ≈ 11.91.
   A verified certificate `< 12` at 3.99 would be `s(12) ≥ 3.99`, past the pure ceiling.

## Compute

Up to 8 threads and 50 GB for the LP loops (a leaf round ≈ 20 min with the restricted master; a
leaf is 40–80 rounds — run it in the background while you do 3–6).  Another unrelated job
(`statarb`) may use ~10 cores; check `uptime`, pin with `taskset` (cores `n`, `n+16` are siblings).
`runs/` is gitignored: read `/home/evand/math/square-packing/s12/runs/` by absolute path, write your
own `runs/`.  **Before anything else, check your worktree's base** (`git merge-base branch-corners
HEAD` must be the current tip of `branch-corners`; if not, `git merge branch-corners` first).

## Done when

1–7 with `verify.sh`, `verify_branch.sh`, the rejection tests and `lake build` green, and a note
`search/ANCHOR.md`: what was built, the bit-identity check, the leaf run table (round, rows,
columns, clique columns used, LP value, probe min, wall time), and the verdict — leaf closed or
not, and if not, where the LP value stalled and what the dual says.  State for every number whether
it is verified, heuristic, or a bound in which direction.  Commit in your worktree; do not touch
`TODO.md`, `README.md` or other tasks' files.
