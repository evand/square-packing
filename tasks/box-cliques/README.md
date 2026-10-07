# B. Box cliques: verifiable clique columns

**Goal.** A certificate may consist of weighted *cliques* of the pose-overlap graph, points being
the special case; the reduction is one line (a packing has at most one square per clique).  Make
that verifiable end to end.

**Family.** A **box clique** is a finite union of pose-space boxes `B_i = [cx-,cx+] × [cy-,cy+] ×
[θ_k, θ_{k+1}]` (centre rectangle in container coordinates × one angle bin of the verifier's net).
Its **core** `core(B_i)` = the intersection of all closed unit squares with pose in `B_i` — the
concentric σ-square at `θ_k` that the verifier already uses for the angle bin, shrunk further by the
centre rectangle (exact for rational rotations).  `K = ∪ B_i` is a clique iff `core(B_i) ∩ core(B_j)
≠ ∅` for all `i, j` (closed intersection suffices — see `tasks/clique-ceiling` for why touching
counts).  No geometric lemma beyond the σ-lemma is needed.

Why this and not `K(p, A) = {S ∋ p} ∪ {S ⊇ A}`: `K(p,A)` needs a non-avoidance lemma per `(p, A)`, and
the measured cliques (`search/CLIQUE.md`: 163 poses through a wall point + 58 interior tilted
poses) only need the *support's* poses through `p`, which a union of boxes covers directly.  Keep
`K(p,A)` as a fallback if the point part must stay whole.

**Deliverables.**
1. `certificates/FORMAT.md`: a `clique` block — list of boxes (rational endpoints, angle bin index
   or `p/q` endpoints), weight over `10^7`.  Total weight counts each clique once.
2. Rust `verify/`: in the sweep, a cell gets `w_K` iff the cell's angle bin is `K`'s and its centre
   cell is inside one of `K`'s rectangles (rectangle edges become arrangement segments; straddling
   cells get nothing — conservative).  Before the sweep: check pairwise core intersection exactly
   and refuse otherwise.  With no clique blocks the output must be bit-identical to today.
3. `xcheck.py`: the same in exact rationals, independently.
4. `tests/rejection_tests.sh`: boxes whose cores do not meet; a clique straddling angle bins; a
   rectangle outside the admissible box; weight applied to a cell not inside; malformed blocks.
5. `lean/Sqpack/Basic.lean`: `packing_le_weight_cliques` — for a family of pose sets `K_j`, each
   pairwise closed-intersecting, weights `w_j ≥ 0` with `∑_{j : S ∈ K_j} w_j ≥ 1` for every admissible
   `S`, `n ≤ ∑ w_j`.  State the box-core clique property as a hypothesis; prove it if cheap.
6. `search/tighten.py` (or a small new driver): clique columns priced by max-weight clique on the
   dual measure, boxed (one small box per member pose, members dropped if their cores stop meeting).
7. **One verified certificate** using at least one non-point clique, at a `t` with slack (3.975 is
   fine; the point is the format, not the bound).

**Compute.** Light until the final verification (`N = 2000` is enough for a demonstration).

**Done when.** 1–7, with `verify.sh` and the rejection tests green, and a note `search/BOXCLIQUE.md`.
Do not touch `TODO.md` or other tasks' files.

## Status (2026-08-29)

Done, 1–7; write-up `search/BOXCLIQUE.md`.

1. `certificates/FORMAT.md` "Clique certificates": `cliques N Q c` block, boxes `(k, rect)` in the
   frame of bin `k`, closed rectangles over `Q`, weights over `W`; the core and the pairwise
   condition are spelled out.  Tied to `N`.
2. `verify/`: parser, `check_cliques` (exact `i128` separating-axis test on the cores with the
   verifier's own `σ_k`; empty core refused), per-cell credit in the sweep (exact integer
   containment), straddling cells get nothing and their witness is placed outside the box, full
   `[0,90)` sweep whenever cliques are present.  **Bit-identical** to the previous binary without a
   clique block (main certificate single-threaded and in witness mode, 224- and 56-point files,
   a per-box region trailer in witness mode).
3. `xcheck.py`: independent, exact rationals, exact `σ_k`; agrees with the Rust on every file.
4. `tests/rejection_tests.sh`: 62 → 91 checks (§18–20).
5. Lean: `packing_le_weight_cliques`, `card_filter_clique_le_one`, `clique_of_cores` — 0 sorries,
   axioms `propext, Classical.choice, Quot.sound` (checked with `#print axioms`).
6. `search/boxclique.py`: clique-orbit columns priced on the D4-averaged interior-point dual,
   boxed, grown, pruned to a pairwise-certified family, exported.
7. `certificates/s12_boxclique_demo_3.9318_N2000.txt`: 224 points + 8 cliques (48 boxes each),
   total 11.9939 < 12, verified by Rust and `xcheck.py --all` at `N = 2000`; in `verify.sh`,
   `SHA256SUMS`, rejection tests.  **The clique is not load-bearing**: on the fixed 224-point set
   the LP never uses a clique, at its own container or at 3.9357 where the points alone give
   12.015 — the priced cliques have ȳ(K) up to 2.1 under the returned dual, but the optimal dual
   face of a 29-column LP is large enough to dodge every offered clique.  Whether box cliques
   carry the 0.1–0.3 of `CLIQUE.md` must be tested inside point column generation (rich atom
   set, pinned dual) or on the packing side (task A); see `search/BOXCLIQUE.md`.

Not done: no clique+region combination in Lean (the verifier accepts both blocks together; the
Lean statement for that combination is a routine merge of `packing_le_weight_regions` and
`packing_le_weight_cliques`, not written); `TIGHT_DUMP` refuses clique certificates.
