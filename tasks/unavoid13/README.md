# unavoid13: is there a 13-point pure unavoidable set for `[0,4]^2`?  Decide it by a hitting-set loop  (2026-09-22)

**Why.**  The angle-space route is closed (`notes/review-2026-09-22.md`).  The one open door is the Kearney–Shiu
slack-1 accounting at `T = 4`: 12 squares against a **13-point** set unavoidable for closed unit squares in `[0,4]^2`
(`notes/t3-chord.md` §6.2).  Known: `>= 13` is a theorem (`search/COVER4.md`: `COVER^closed(4) >= 12.2688`, integrality);
`<= 14` (Friedman, DS7 Thm 4; certified exactly by `search/zeromargin.py friedman14 --tri`); the fractional optimum is
believed `~12.4` and a certified weighted cover of weight `12.956 < 13` exists (`certificates/rung2/s13_closed_cover_4.txt`).
**No LP bound can decide 13 vs 14.**  This task decides it by finite computation — in either direction.

**Read first.**  `search/ZEROMARGIN.md` §1–2 (semantics: closed unit square, closed containment, pose `(cx, cy, u = tan(θ/2))`,
admissible iff `cx, cy ∈ [w/2, m − w/2]`, `w = |cos| + |sin|`; a point on `∂Q` counts), `certificates/FORMAT.md`,
`notes/t3-chord.md` §2.1 (the K–S 7-point set and how tight it is) and §6.2, `search/COVER4.md` (the 294-pose / 2352-square
packing measure — its support file `search/cover4_exact_support.txt` is the natural starting family), `notes/proof-anatomy.md`
§7.2 (the deficit table).  Tools: `search/zeromargin.py cert FILE [--tri] [--oracle OUT]` (exact checker; **an uncertified box is
not a violation** — confirm a suspected violation with `zeromargin.py pose FILE --cx --cy --u`, exact), `highspy` and
`scipy.optimize.milp` are installed (no pulp/ortools).

**The loop (the whole method).**  Let `F` be a finite family of closed unit squares in `[0,m]^2`.  Any set hitting every
closed unit square hits `F`, so `h(F) := min |P| with P ∩ Q ≠ ∅ for all Q ∈ F` is a **lower bound** on the pure unavoidable
number, valid **only if the candidate points range over every cell of the arrangement of `F`** (the point set is
otherwise unconstrained; a cell is characterised by which squares of `F` contain it, so one candidate per cell, or per
maximal cell, loses nothing).  Iterate:
1. Build candidates from the arrangement of `F` (cells: pairwise edge intersections, square corners, plus one interior
   sample per face — or, simpler and still exact for the lower bound, every vertex of the arrangement plus the
   containment-maximal faces; say how you did it and why it loses nothing).  Solve the hitting-set IP (`highspy`), exactly
   (integrality gap closed, not a heuristic).
2. If `h(F) >= 14`: **no 13-point set exists** — `F` plus the IP's optimality proof is the certificate.  Make it
   reproducible: dump `F` (rational poses), the candidate cells, the incidence matrix, and re-solve from the dump.  Try to
   shrink `F` (drop squares while `h` stays 14) so the certificate is small enough to reason about.
3. If `h(F) = 13`: take the 13-point solution, snap to rationals, check with `zeromargin.py cert` (unit weights, `W = 1`).
   If certified: **done, a 13-point set exists** — publish it as a certificate file.  If not: find an actual violated
   pose (float search for `min_Q max_p inside(p, Q)` seeded at the uncertified boxes / `--oracle` output, then `pose` mode
   to confirm exactly), add it (and its `D4` images) to `F`, go to 1.  Also add a few *nearby* violated poses per round,
   not one, to keep the round count down.
4. Report `h(F)` per round.  If it stalls at 13 with the verifier finding ever-smaller violations, that is the fractional
   gap talking: switch to structured candidate families (below) and to *symmetry-reduced* IPs before giving up.

**Rehearsal first, `T = 3` (Evan's rule).**  Same loop on `[0,3]^2`.  It must rediscover a 7-point set (K–S's is
`{(√2−½,1),(3/2,1),(7/2−√2,1),(3/2,3/2),(√2−½,2),(3/2,2),(7/2−√2,2)}`; the loop may find a different one — fine, certify it).
Then push: **is 7 the minimum at `T = 3`?**  If `h(F) = 7` for some finite `F` at `T = 3` that is a theorem nobody has stated;
record it.  Calibration numbers to report alongside: the `T = 3` fractional floor (heuristic is fine — a float LP over a
dense pose grid with the cover-side value and the packing-side measure both reported; label **[heuristic]**; only if it is
cheap with `search/nu_f.py` / `search/packing_dual.py`, otherwise skip) so that the `T = 3` integrality gap
(fractional → 7) can be set beside `T = 4` (`~12.4` → 13 or 14).  Do not spend more than ~1 h on `T = 3`.

**Structured families to try at `T = 4` if the raw loop stalls.**  (i) K–S's pattern transferred: 4 points on each of
`y = 1, 2, 3` plus `(2,2)`; my back-of-envelope says 45° squares centred on `y = 1.5, 2.5` kill it (line points have only
`±0.207` of `x`-tolerance for them) — check, and record the killing pose.  (ii) `C4`- or `D4`-symmetric sets: 13 odd forces
the centre `(2,2)` plus three `4`-orbits (or one `4`-orbit and one `8`-orbit) — 6 or 4 real parameters, cheap to
optimise directly on `max_Q min_p dist(p, Q)` by multistart.  (iii) Friedman's 14 minus one point plus one moved point.
(iv) The support of the certified `12.956` cover (`3621` points) as the candidate set, IP restricted to it (an upper-bound
search, not a lower bound).

**Sanity checks you must run.**  The loop at `T = 4` with Friedman's 14 loaded as a warm start must accept it (checker
passes, as `ZEROMARGIN.md` reports).  Every `h(F)` claim: re-solve the IP from the dumped incidence matrix with a second
method (e.g. `scipy.optimize.milp` vs `highspy` directly, or a hand branch-and-bound on the 13-vs-14 question) before
writing "theorem".  Every "unavoidable" claim: `zeromargin.py cert` with `0` uncertified boxes, or every uncertified box
individually confirmed by `pose`/`diag`.

**Deliverable.**  `notes/unavoid13.md`: verdict up front — one of "**13 exists** (certificate `certificates/unavoid13/…`)",
"**13 excluded** (finite family `F`, `h(F) = 14`, reproducible)", or "undecided: `h(F) = 13` stalled at round N, smallest
violation ε, what was tried"; then the `T = 3` rehearsal result (is 7 minimal?), the calibration numbers, the round log,
what in this brief turned out wrong.  Scripts `search/unavoid13_*.py` (import, never modify, existing ones), runs
`runs/unavoid13_*`.  Label **[proved]** / **[measured]** / **[heuristic]** / **[guess]**.  Do not edit `TODO.md`,
`notes/status.md`, existing `notes/*.md`, existing `search/*`; do not commit.

**Working style.**  `<= 6` threads (another agent is writing a note; the machine's useful width is 16).  Anything over
10 minutes: `setsid nohup`, never `pkill`, no `until … sleep` / `tail -f` watchers (a `pgrep -f <name>` loop matches
itself).  Write the note section by section as results land.  Stop when the verdict is clear; a clean "13 excluded, here
is `F`" is a full deliverable and the best possible outcome after "13 exists".
