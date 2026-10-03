# s11-ladder: put `n = 11` on the same footing as `n = 12`  (launched 2026-09-14, supersedes `tasks/s11-cliques`)

**Why.**  Every comparison of the two problems so far ("the required gap at `n = 11` is twice the one at
`n = 12`", `N11.md`; "the clique family cannot reach the `0.75`", `N11_ANATOMY.md` §7) compares the *pure*
method at `n = 11` with the *fully instrumented* method at `n = 12`.  At `n = 12` the certified record
(`3.968616`) is pure-LP, but the ceiling was then measured with every stronger family — cliques
(`CLIQUE_CEILING.md`, `CLIQUELEVER.md`), corner branch (`BRANCH.md`, `BENTZ.md`), polygons / regions / chord
(`LEAF_CEILING.md`), anchor cliques on the cover side (`ANCHOR.md`), SA-2 (`RANK8.md` §4).  At `n = 11`
only the pure rung has been run (`N11.md`, `N11_ANATOMY.md`), plus one clique diagnostic (`+0.25` on 8 poses).
The gap statements are therefore not baselined.  This task runs the same rungs at `n = 11`, and a new
rigorous `s(11) ≥ t` is a by-product if any cover-side rung buys one.  Not the goal: an exact `s(11)`.

**The ladder** (n = 12 value → n = 11 value to be filled in; conjectured values `4` and `3.877083`):

| rung | `n = 12` | `n = 11` |
|---|---|---|
| 0 literature | Stromquist `3.788854` (inherited) | Stromquist `3.788854` |
| 1 pure certificate | `3.968616` | `3.814304` (`certificates/s11_lower_3.8143.txt`) |
| 2 pure ceiling (exact measure of mass `≥ n`) | `< 3.99` | `≤ 3.83375`; window `[3.8143, 3.83375]` |
| 3 QSTAB (all cliques) exact measure, packing side | at `t = 4`: no measure `≥ 12` found in 40 core-h; fixed-pose leaf values `11.30–11.67` | `clique_ceiling.py --exact` at `3.82, 3.83375, 3.85, 3.86, 3.87` |
| 4 corner branch | `k = 4` leaf at `t = 4`: `11.9999999`; `k ≤ 2` dead at `t ≤ 3.98` | vacuous for the pure optimum (corner boxes saturated); re-measure under cliques (`--kmass 4`) |
| 5 box-clique certificate, cover side | never load-bearing on a fixed point set | `boxclique.py --n 11` on `n11_H3985.txt` at `Dp = 3980, 3975` (`t = 3.8191, 3.8239`) |
| 6 anchor-clique certificate, cover side | matched gain `10⁻⁴` | `branch.py --cliques --matched --n 11 --lam-lo 0 --lam-hi 0` (after 5) |
| 7 SA-2 on the extremal support | `11.893` on 162 poses, `0.107` of the missing unit | `search/n11_sa2.py` on the 28-square support |

**Runs** (all detached from the main tree, `setsid nohup`, ≤ 16 threads total; logs `runs/cc_N11Q*.log`,
`runs/boxclique_N11BC*.txt`, stdout in the session scratchpad).  Warm starts: float conversions of
`search/n11_exact_{3.83375,3.85}_support.txt`.  Flags as in `CLIQUE_CEILING.md` "Reproduce" with
`--inner 400 --time 14400`.

**Deliverable.**  `search/N11_LADDER.md`: the table above filled in, each entry marked exact / float, with
the same "gap to the conjecture" column for both `n`; any new certificate into `certificates/` with
`verify.sh`, `xcheck.py --n 11`, README and site updated; `notes/status.md` line for `n = 11` rewritten.

**Done 2026-09-14:** `search/N11_LADDER.md`.  All runs finished or stopped; §3 (leaf, certified `32/3` with corner mass `4`) and §4 (box-clique retry, `0` cliques used) are final.
