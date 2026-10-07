# bandcut-scan: the band construction from `T = 12` down to `T = 4`  (2026-09-20)

**Why.**  `notes/proof-architecture.md` §0a: every lemma we can state is true at every `T`; the `T`-dependence is in
configurations with `T <= k <= (T-1)^2` near-axis squares cut by a tilted band.  `s(T^2-T) < T` is a theorem for all
`T >= 12` (Arslanov–Mustafin–Shangitbayev, EJC 28(4) 2021, P4.22, https://www.combinatorics.org/ojs/index.php/eljc/article/download/v28i4p22/pdf/)
and reported at `T = 11` (Cantrell 2025; `search/ARCH_TLEDGER.md` §5).  Our `delta*` machinery has never seen a positive case.

**Do.**  (1) Rebuild the `T = 12` packing (and `T = 11` if a description can be found) as explicit centres + angles; verify
`delta > 0` with the exact pair/wall rows of `search/S6_SKELETON.md` §3.1 — this validates the instrument on a known positive.
(2) Parametrise the family (band position, angle(s), elbow, which cells are removed) so it is defined for every `T`, and for
`T = 12, 11, 10, ..., 4` maximise `delta` over the family's parameters and over the centres (fixed-assignment LP in the centres,
outer search on angles; `search/s6local.py` has the LP, `search/arch_chain.py analyse` reads chains and rises).  Report best
`delta(T)`, the binding constraints (dual), and what runs out as `T` falls — length for the band? the elbow? the rise budget?
(3) At `T = 4`: the best band-type configuration, its `k`, and how far below `0` it sits.  Try variants (straight band,
bent band, two bands, bands of 3-chains / 5-chains at `2 arctan(1/j)`).
Everything found is a feasible point (lower bound on `delta*`).  A positive `delta` at any `T <= 10` would be a new result:
verify it exactly (rational arithmetic) before saying so.
**Deliverable.**  `search/BANDCUT_SCAN.md`, scripts `search/bandcut_scan*.py`, runs `runs/bandcut_scan_*`.
