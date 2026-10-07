# unavoid13-no: one bounded push for the rigorous "no 13-point unavoidable set for `[0,4]^2`"  (2026-09-22)

**Where it stands.**  `tasks/unavoid13` (`notes/unavoid13.md`, read all of it first) decided `T = 3` (minimum 7, proved),
excluded every 13-set with a 180° symmetry at `T = 4` (proved; `Rx`, `Rd` still running in `runs/unavoid13_symRx`, `_symRd`),
and stalled on the general question: `h(F) = 13` at every round, the IP the bottleneck (10,673 binaries, presolve removes
nothing, 10+ min/round), the LP over every family only `11.3–11.6` although `ν_f^closed(4) >= 12.2688` is a theorem
(`search/COVER4.md`).  The evidence (`T = 3` gap `>= 1.0`; asymmetry forced; maximin margins never above `−0.028`) leans
"no".  **This task's goal is a finite family `F` of closed unit squares with `h(F) >= 14`, exact and re-solved, i.e. a
theorem that no 13-point set exists** — or a precise statement of why that cannot be reached with this machinery.  A "yes"
would show up as a by-product (an IP 13-set that certifies); verify it with the SEG-augmented checker
(`search/unavoid13_check.py`) if it does, but do not search for one.

**Tools, all in place** (import, never modify): `search/unavoid13_lib.py` (arrangement cells, exact incidence, dominance),
`_loop2.py` (lazy-row feasibility loop; `--support search/cover4_exact_support.txt`, `--init tiles,grid2,support:N`),
`_recheck.py` (exact re-solve, three solvers, `--shrink k`), `_check.py` (checker with the SEG primitive), `_sym.py` /
`_symcheck.py`.  Formats: `certificates/FORMAT.md`; `search/cover4_exact_support.txt` (294 poses `p q cx cy mass`,
`theta = 2 arctan(p/q)`, each standing for its 8 `D4` images — 2,352 squares; this is the *fractional* packing measure of
mass `12.2688`).  Running detached from the previous task — **do not kill or restart them**; read their logs when useful:
`runs/unavoid13_t4L/round_log.txt` (feasibility loop), `runs/unavoid13_symRx/`, `_symRd/` (reflection cases; when both are
infeasible and exactly re-checked, "no 13-set has any symmetry of the container" is a theorem — record it if it lands during
this task), `runs/unavoid13_maximin*`.

**Do, in this order; stop as soon as `h(F) >= 14` is exact.**

1. **Built-in completeness test, then the right `F`.**  With `F ⊇` the full 2,352-square `COVER4` support and candidates =
   *all* cells of the arrangement (or a dominance-equivalent set), the hitting-set LP over `F` must be `>= 12.2688038611`
   (its dual contains the certified measure).  If your LP is below that, the candidate set is incomplete — fix that before
   anything else; this is the soundness check for every lower bound that follows.  Report the LP value.  Then add the
   previous task's families (`runs/unavoid13_t4b/F_round07.txt` or later, `t4L`) and the tiles/shifted tilings.
2. **Reduce before solving.**  Cells: keep containment-maximal ones only (sound for hitting sets).  Squares: if
   `cells(Q1) ⊆ cells(Q2)` drop `Q2` (hitting `Q1` hits `Q2`).  Report sizes before/after.
3. **Structure first, brute force second.**  Compute `h` exactly on the *axis-parallel* subfamily alone (tiles, every shifted
   `4×3`/`3×4` tiling at a fine shift grid, the `3×3` shifted tilings) and describe its 13-sets (do they all put 9 points near the
   interior lattice points?  where do the 4 extra go?).  Then add the `45°` squares only, then the near-axis (`< 5°`) ones,
   then the rest of the support — tracking `h` and the *shape* of the 13-sets at each stage.  The aim is to understand what
   a 13-set must look like so that the eventual `F` with `h = 14` is small enough to read, and so a hand argument
   ("9 near-lattice points are forced, then the 4 free points cannot serve both the `45°` wall squares and …") becomes
   visible.  If the 13-sets over rich `F` share a rigid skeleton, say so precisely and branch on it (skeleton / no
   skeleton) — each branch is a smaller IP.
4. **A purpose-built branch-and-bound**, not a generic MIP, for the feasibility question "13 points hitting `F`?": branch on
   the square with the fewest candidate cells (every hitting set hits it), depth `<= 13`, LP bound (or a cheap matching
   bound) at each node, incremental incidence.  Exploit `D4`: `F` is `D4`-closed, so if `P` hits `F` so does `gP`; use that to
   cut the root branching (e.g. the point hitting a `D4`-invariant square orbit representative may be taken in a fundamental
   domain of its stabiliser — write the rule down and argue its soundness before using it).  Compare against
   `highspy` on the same instance for calibration; if HiGHS with a long time limit (hours, detached) proves infeasibility
   on a family where your B&B is slow, that is an acceptable proof too **provided** the instance is dumped and re-solved
   by a second method (scipy `milp`, or your B&B, or HiGHS with different settings/seed).
5. **The family loop, restarted from the richer `F`.**  Only after 1–4: lazy rows (`_loop2.py`) or your own loop; add the
   worst violated poses (exact `pose` confirmation) plus their `D4` images; the goal is infeasibility, so favour adding
   *many* diverse violated squares per round over polishing.  If `h(F)` stays 13 with LP `>= 12.27` and the B&B is
   exhausting the tree in minutes, report the sequence of 13-sets and their violation depths — that is the "cannot reach"
   statement.

**Sanity checks you must run.**  Step 1's LP `>= 12.2688`.  Friedman's 14 (`certificates/unavoid13/friedman14_4.txt`) hits every
`F` you build (exact).  Any `h(F) >= 14` claim: dump `F` (rational poses), cells, incidence; re-solve with two independent
methods; then shrink `F` (drop squares while `h` stays 14) and re-solve the shrunk family — the shrunk family is the
certificate.  Any 13-set that the checker certifies: report it first, before anything else.

**Deliverable.**  `notes/unavoid13-no.md`: verdict up front ("13 excluded: `F` of N squares, `h(F) = 14`, certificate at …",
or "not reached: LP/IP numbers, the 13-set skeleton, what would be needed"), then steps 1–5 with numbers, then what in this
brief was wrong.  Scripts `search/unavoid13no_*.py`, runs `runs/unavoid13no_*`, certificate `certificates/unavoid13/no13_*`.
Label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  Do not edit `TODO.md`, `notes/status.md`, existing
`notes/*.md`, existing `search/*`; do not commit.

**Budget and style.**  One agent-day.  `<= 10` threads (the detached runs hold ~8).  Anything over 10 minutes: `setsid nohup`,
never `pkill`, **no `until … sleep` / `tail -f` / `tail -F` watchers** (the previous agent left one; a `pgrep -f <name>` loop
matches itself).  Write the note section by section as results land.  Be adversarial with this brief and with
`notes/unavoid13.md`: check its soundness argument for the symmetric exclusions (§4.2 (a)/(b)) and for the SEG primitive
(§2.2) as you use them, and record anything wrong.
