# G. Cliques in the continuum: is the non-Helly mass real, and what family certifies it?

**Why.** Phase 1 left the clique lever in doubt.  `search/CLIQUE_CEILING.md` (task A) found that the
non-Helly excess of every extremal measure sits on **grazing contacts** — dropping the six members
with pairwise contact margin `< 1e-3` takes the max clique at 3.99 from 1.327 to 1.0096 — and
`search/BOXCLIQUE.md` (task B) found that on a fixed point set the cover LP never puts weight on a
box clique.  A grazing contact is not stable under perturbation: nudge the wall square and the
interior square no longer touches it.  So the measured clique mass may be partly an artefact of
finite pose sets, and the *continuum* clique-LP value at 3.99 may be `>= 12` after all.  A's stage
values drift up with pose refinement (11.19 → 11.36 → 11.45).  **Nobody has determined whether cliques
buy anything in the continuum.**  This is the go/no-go for the whole clique line, and half of it is
mathematics.

Read first: `search/CLIQUE.md`, `search/CLIQUE_CEILING.md` (esp. "structural finding"),
`search/BOXCLIQUE.md` ("The negative result"), `search/ZEROMARGIN.md` §2 (exact bin cores, the
2×2-box lemma, the triangle lemma — the primitives available for zero-margin certification),
`search/DUAL_EXACT.md` (how a packing measure is certified), `notes/proof-anatomy.md` (the
non-avoidance lemmas of the literature), `certificates/FORMAT.md`.

## Part 1 — the family (mathematics first, then code)

Semantics: a certificate at container `t` refutes packings at side `< t`; rescaled, those are pairwise
**disjoint closed** unit squares, so a clique is a set of poses that pairwise **closed-intersect**
(touching counts), and a clique's weight counts once.

1. Characterise the cliques the measures actually use (A's 3.99 measure: 163 poses through a wall
   point of mass 1.000 + 58 interior tilted poses of mass 0.327 touching them all).  Which of these
   survive as cliques of a *continuum* of poses?  Specifically: for the wall point `p`, the family
   `{S admissible : p ∈ S}` is an exact clique; which interior pose *regions* `R` have the property
   that every `S' ∈ R` closed-meets every admissible `S ∋ p`?  Is that set empty, a lower-dimensional
   set, or a region of positive measure?  (Guess: near a wall, admissibility prunes the squares
   through `p` enough that a positive-measure region exists; find it or refute it, with proof.)
2. Define a certifiable family that includes the touching limit.  Candidates: **point-anchored pose
   boxes** `{S ∈ B : q ∈ S}` for a pose box `B` and point `q` (point cliques are `B` = everything);
   `K(p, A) = {S ∋ p} ∪ {S ⊇ A}` (valid iff every admissible `S ∋ p` meets `A`); unions.  For each pair
   type, state the pairwise closed-intersection condition and an *exact, terminating* certification:
   exact-core subdivision (`search/zeromargin.py` machinery) where margins are positive, and a named
   lemma primitive where they are zero — write the lemma, prove it (one-parameter trigonometry is
   fine), and stress-test it numerically the way `zeromargin_stress.py` does.  Say explicitly where
   subdivision fails to terminate without the lemma.
3. Say what changes in `verify/` and `xcheck.py` (credit rule in the sweep) and in the Lean statement
   (`packing_le_weight_cliques` already takes any pairwise-intersecting family; what remains is the
   clique property of the new family as a hypothesis or lemma).  Design only; do not implement the
   Rust unless the answer to Part 2 is "go" and time remains.

## Part 2 — the continuum value at 3.99 (compute, packing side)

A's column-generation loop was degenerate and uninformative.  Replace it by a **fixed fine lattice**
of poses: D4-reduced centres at pitch `h` (0.02, 0.01, 0.005) × angles at pitch `dθ` (2°, 1°, 0.5°),
plus the poses of the certified measures `runs/dual_PA2_support.txt` etc.  Packing LP over those
poses with (a) point rows at all arrangement vertices (or a fine point lattice, then certify coverage
`<= 1` at *all* arrangement vertices with `dual_exact.py` as before) and (b) lazily separated
closed-intersection cliques (`clique_ceiling.py`'s exact SAT with `<=`; greedy + B&B separation).
Every certified measure is a valid **lower bound on the continuum clique-LP value** at that `t`.

**Calibration rule (mandatory):** on the same lattice, the LP *without* clique cuts must reproduce the
pure value `L(3.99) = 12.008` to within ~0.01; if it does not, the lattice is too coarse and its
clique value says nothing.  Report both numbers per lattice.

Deliver, per `(h, dθ)`: pure value, clique value, residual max clique, certified mass.  Read the
trend: if the clique value climbs to `>= 12` as the lattice refines (or a certified clique-feasible
measure `>= 12` appears) the clique method is dead at `t >= 3.99` — say so, that is a result.  If it
stays clearly below 12 on lattices fine enough to pass calibration, estimate the limit and say which
cliques carry the constraint (support poses, pairwise margins: are they grazing, or do they have
positive-measure overlap?).  Then the same at `t = 3.98` in the corner leaf `k = 4` (`--kmass 4`,
`r = 1`), whose no-clique value is `12.000 ± 0.001`.

## Compute

`<= 4` cores, `<= 30 GB`.  Another job (`statarb`, not ours) uses ~10 cores; check `uptime` before a
long run.  `runs/` is gitignored: read the main tree's `/home/evand/math/square-packing/s12/runs/`
by absolute path, write to your own `runs/`.  Build the verifier in your worktree if you need it
(`cd verify && cargo build --release`).

## Done when

`search/CLIQUE_CONTINUUM.md`: the calibration/value table with a clear verdict (go / no-go / undecided
and why), and `notes/clique-family.md`: the family, the lemmas with proofs, the certification
procedure, and the verifier/Lean changes it would need.  Code under `search/`.  State for every
number whether it is certified, heuristic, or a bound in which direction.  Commit in your worktree.
Do not touch `TODO.md` or other tasks' files.
