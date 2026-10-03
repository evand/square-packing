# line-cover: s(21) = 5 via mixed (point + line + region) covers, 2026-09-27

Two agents in parallel, same working tree `~/math/square-packing/public/s12`, disjoint files.  Format: `FORMAT.md` here.

- **A (cover side)**: `A.md`.  Go/no-go gate for the whole project.
- **B (checker side)**: `B.md`.  Exact checker for mixed covers.

Gate (A, heuristic): at `m = 5` a mixed cover with honest cost (strict protocol: total / strict-scan minimum) `< 20.947`
(= 21 / 1.0025, leaving the checker margin measured at `m = 4`).  If A fails the gate, B's work still serves lighter
`s(32)`/`s(45)` covers.

Common rules: commit your own files only (`git add <paths>`), **never push**; `runs/`, `tasks/`, `TODO.md` are not
committed.  **Do not modify `search/zeromargin.py`, `verify2/`, or `certificates/`** (checker sha is pinned in the shipped
bundles); import, don't edit.  Detached launches (`setsid nohup … &`, kill by PID, never `pkill`); no leftover polling loops.
Don't hand back while your own jobs are running.  Labels: **certified** / **heuristic** / **estimate**, as in `S21_COVER.md`.
