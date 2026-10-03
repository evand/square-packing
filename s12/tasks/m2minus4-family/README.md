# D. The m²−4 family: is the fractional gap universal?

**Question.** For `n = m² − 4`, `m = 5, 6, 7` (n = 21, 32, 45), does the closed-semantics cover LP
at `[0,m]²` go below `n`?  The literature's best pure point sets have deficit 2, 2, 1 against
`n − 1` (`notes/proof-anatomy.md §7.2`); Bentz's 45-point set at `[0,7]²` already gives `W = 45`,
so for `n = 45` *any* weighted improvement plus exact zero-margin verification (task E) is a new
theorem `s(45) = 7`.  And if the LP is `≥ n` for all three, the gap at `n = 12` is not a small-`m`
wall artefact.

**Plan.** Generalise `search/closed4.py` (container side `m`, D4 orbits, grid-line columns, lattice
+ polish rows) and `search/nu_f.py lower` / `dual_exact.py` for the packing-side lower bound.  Run
`[0,5]²` first (cheapest), then 6, then 7.  Report cover value (heuristic, sampled rows — rising)
and certified packing mass (rigorous lower bound on the cover).

**Compute.** Heavy: `[0,7]²` is ~3× the area of `[0,4]²`; expect hours on 8+ cores.  Serialize after
`tasks/lp-speed`'s benchmark.

**Done when.** A table `m → cover LP value, certified packing mass, n − 1` in `search/FAMILY.md`.
