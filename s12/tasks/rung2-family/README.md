# H. Rung 2 and the m²−4 family: exact closed covers at the container itself

**Why.** `search/ZEROMARGIN.md` (task E) built an exact checker for closed covers of the closed
container `[0,m]²` — bin cores + 2×2-box lemma + triangle lemma — and re-proved `s(15) = 4` from
Friedman's 14 points in under a second.  Its verdict on rung 2 (`W < 13` at `[0,4]²`, i.e. `s(13) = 4`
by machine) was "compute, not new mathematics": the `closed4.py` LP loop must be run with the exact
checker as its separation oracle, and the weighted core test must be vectorised first.  Then the
same at `m = 5` (`n = 21`, `W < 21`) answers task D's question — is the fractional gap at `n = m²−4`
universal? — and, if the LP goes below 21, is a **new theorem** `s(21) = 5`.

Read first: `search/ZEROMARGIN.md` (all of it), `search/CLOSED4.md`, `search/closed4.py`,
`search/zeromargin.py`, `search/zeromargin_fan.py`, `tasks/zero-margin/README.md`,
`tasks/m2minus4-family/README.md`, `notes/proof-anatomy.md` §5, §7.  Semantics: closed unit squares,
closed containment in `[0,m]²`, points on the boundary count (`certificates/FORMAT.md`); a closed
cover of `[0,m]²` of weight `W < n` proves `s(n) >= m`.

## Steps

1. **Make the checker an oracle.**  Vectorise the float pre-filter of `zeromargin.py` (numpy over the
   point set) so that a 2,000-point weighted cover runs to depth 10–12 in minutes, not hours; add an
   output mode that emits the uncertified boxes' worst poses (box corners and centre, exact) as LP
   rows.  Keep the exact `Fraction` path as the final word; re-run `zeromargin_stress.py`-style checks
   after any change.  Confirm rung 1 (Friedman's 14 points) still certifies, byte-identical leaf
   count or explain the difference.
2. **Rung 2, `m = 4`, `W < 13`.**  Run the `closed4.py` cover LP (D4 orbits, grid-line columns) with
   the checker as separation oracle from `runs/closed4_best.txt` (12.4175, fails by 0.7 % at tilted
   wall poses).  Goal: a cover the checker certifies with **0 uncertified boxes** and `W < 13`.
   Ship it as `certificates/s13_closed_4.txt` with the checker's leaf summary, and record whether
   any tilted tight family (ZEROMARGIN §4 item 4) appeared and how it was handled.  If the loop
   stalls, report the residual uncertified boxes (where, what margin) honestly.
3. **`m = 5`, `n = 21`, `W < 21`.**  Generalise `closed4.py` to container side `m` (orbits, grid lines
   `x, y ∈ {1, …, m−1}`, lattice + polish rows) and run it with the checker as oracle.  Report the LP
   value trajectory (heuristic, sampled rows — it can only rise), and in parallel a **certified
   packing-side lower bound** (generalise `search/nu_f.py lower` / `search/dual_exact.py` to `m = 5`,
   exact rationals) so the answer is bracketed: if the packing measure reaches `>= 21` the pure method
   is dead at `m = 5` too and that is the result; if the cover goes `< 21` and certifies, that is
   `s(21) = 5`.  Only then `m = 6` if time remains.

## Compute

Up to `8` threads and `50 GB`; another job (`statarb`, not ours) uses ~10 cores — check `uptime`, and
pin long runs with `taskset` to distinct physical cores (cores `n` and `n+16` are hyperthread
siblings).  `runs/` is gitignored: read `/home/evand/math/square-packing/s12/runs/` by absolute
path, write to your own `runs/`.

## Done when

`search/FAMILY.md`: table `m → cover LP value (heuristic), certified cover (if any), certified packing
mass, n − 1`, plus the rung-2 outcome; `search/ZEROMARGIN.md` updated with the oracle mode and the
rung-2 result.  Certificates under `certificates/` with the exact checker command that verifies them.
State for every number whether it is certified, heuristic, or a bound in which direction.  Commit in
your worktree.  Do not touch `TODO.md` or other tasks' files.
