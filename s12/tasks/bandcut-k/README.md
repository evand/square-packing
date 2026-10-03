# bandcut-k: what does the zero set look like with `4 <= k <= 9` near-axis squares, and is there always a chain?  (2026-09-20)

**Why.**  `notes/proof-architecture.md` §0a items 1–2, 8: chain counting only forces a wall-to-wall chain of `T` near-axis
squares when `k > (T-1)^2`; the generic stratum of `Z_4` (chain of four + eight squares tilted up to `12.8°`) has `k = 4`;
Cleemann-type packings at large `T` have `T <= k <= (T-1)^2` and **no** chain.  The lemma the proof needs is some form of
"no wall-to-wall chain among the near-axis squares `=>` `delta* <= -c < 0`" at `T = 4`.  Measure it before anyone tries to prove it.

**Do.**  `T = 3, n = 6` first (hole is `k in {3, 4}`), then `T = 4, n = 12`.  "Chain" = a path of `T` squares in the DAG of
x-type (or y-type) separations, staircases included (`search/RANK8.md` §3.2, `search/arch_chain.py`, `search/chains.py`).
(1) For each `k` and `eps in {1°, 5°, 10°}`: max `delta*` over angle vectors with exactly `k` squares of tilt `< eps` and
the rest `>= eps`, (a) unrestricted, (b) restricted to configurations whose near-axis squares contain no chain of `T`.
Signed tilts (`search/S6_LOCAL.md` §5: sign coherence matters).  (2) Classify the `delta* = 0` configurations found with
`k` in the hole: is a chain of `T` near-axis squares always present?  Are there zero-margin configurations whose
obstruction involves a tilted square (none seen in `search/CENSUS.md` / `S6_SKELETON.md` §1.2 — re-test in this regime)?
(3) The constant: `c(k, eps)` = best (least negative) `delta*` without a chain; how it scales with `eps`; which
configurations attain it (are they band-like?).  Everything is a feasible point (lower bound on `delta*`); the `n = 12`
multistart is weak (`S6_LOCAL.md` §1) — use structured starts, report hit rates, and say which cells are under-explored.
**Deliverable.**  `search/BANDCUT_K.md`, scripts `search/bandcut_k*.py`, runs `runs/bandcut_k_*`.
