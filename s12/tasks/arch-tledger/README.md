# arch-tledger: where does `s(T^2-T) = T` break?  (§3 of `notes/proof-architecture.md`)  (2026-09-20)

(1) **`c_T`.**  Uniform tilt, `delta* = -c_T t^2`: `c_2 = 1`, `c_3 = 3/4`, `c_4 = 1/3` exact on their leaves
(`search/S6_LOCAL.md` §3); `c_5 <= 0.3745` measured (`runs/s6local_n20_T5_uniform.txt`, re-run in
`runs/s6local_n20_T5_uniform_seed7.txt`) and suspicious (non-monotone; `n = 12` was under-optimised the same way).
Get the optimal `T = 5` leaf (`s6local.py --show`), its exact closed form (`s6exact.py --n 20 --T 5`), look at the
configuration (which cells are empty: a permutation? a 5-cycle?), and try structured starts from every permutation-type
hole set (120) — and at `T = 6` if affordable.  Conjecture a formula for `c_T` and its dual (chain weight `1/(T-1)`?).
Does it change sign?  (2) **Cleemann.**  Find the 272-square packing in side `< 17` (Friedman, "Packing unit squares in
squares", Electronic J. Combin. DS7, Fig. 8 or nearby; web search allowed) and any other known `s(T^2-T) < T`.  Classify:
how many squares axis-parallel, tilt angles, do axis-parallel squares form wall-to-wall chains of 17, which of (α)(β)(γ)
in §3.3 it is.  If coordinates are not published, reason from the figure's structure and say so.  What is the smallest
`T` for which a packing beating `T` is known?  (3) Straight-chain window `tan(t/2) >= 1/T`: verify.
Compute: `<= 8` threads until ~23:00; detached launches over 10 min; never `pkill`; no watchers.  Deliverable:
`search/ARCH_TLEDGER.md`, new scripts `search/arch_*.py`.  `delta*` from `s6local` is a feasible point (lower bound on
`delta*`, upper bound on `c_T`); `s6exact` is exact for one leaf only.  Do not edit existing files.  Do not commit.
