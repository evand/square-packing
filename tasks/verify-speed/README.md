# Verifier speed: ~15x on anchor-clique certificates

**Why.**  `notes/review-2026-09-07.md` Diagnostics: the exact verifier costs 20–120 s per angle
bin at `N = 2000` on the `k = 4` anchor-clique probe (`1650 ns/cell`, ~134 (cell, piece) anchor
tests per cell, 99 % succeeding) against `210 ns/cell` without anchors, and sweeps all 2001 bins
because the D4 check skips certificates with a clique block.  It is the separation oracle of
`search/branch.py`, so it sets the cost of every cover-side leaf (12 h/round for `k = 4`).  Three
levers were identified and not taken: **D4 symmetry with an anchor block (2.4x)**, **incremental
anchor credit (~5x)**, **dynamic bin scheduling (~1.3x)**.  Take them.

**Rules.**  `verify/` is the trusted component.  Every change must keep: (a) the 42+138 rejection
tests passing (`./verify.sh`, `tests/`); (b) bit-identical VERIFIED/REFUSED decisions and identical
reported minimum on every shipped certificate (`certificates/`, `certificates/branch/`) and on the
`k = 4` / `1110` probes (`runs/branch_J16i_probe.txt`, `runs/diag_L1110f_*`) — put the before/after
numbers in the report; (c) soundness arguments written down for each lever (the D4 reduction with
anchors needs the anchor set to be D4-closed, or the orbit images added; the incremental credit must
be conservative at every cell; scheduling changes nothing semantically).  Use
`VERIFY_BINS=lo:hi` for timing before you touch anything, so you have a baseline per bin.  Profile
first (`perf` or `cargo flamegraph` if available; otherwise instrumented timers) and report the
actual breakdown — if the 5x estimate for incremental credit is wrong, say so and take whatever the
profile says is real.  `xcheck.py` must still agree bin by bin on at least one anchor certificate.

**Budget.**  ≤ 4 threads for validation runs (the verifier takes a thread count; measure at 4 and
report projected 16); detached for anything > 10 min; never `pkill`; never touch pid 740647.

**Deliverables.**  Code in `verify/`, `search/VERIFYSPEED.md` (baseline table, per-lever gain,
soundness notes, before/after decisions on every certificate), tests for any new code path
(rejection tests where a lever could over-credit); commit on your worktree branch.  Do not edit
`TODO.md`, `README.md`, or other tasks' files.
