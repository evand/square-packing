# Brief: close the `s(32)` cover to a stable strict minimum (m = 6), 2026-09-25

Working tree: `~/math/square-packing/public/s12`.  Commit there (`git add` your own files only; a sibling agent works
on `s(45)` in the same repo; re-read shared files like `rung2_close.py`, `closed4.py`, `s21_cover_eval.py` right before
editing, keep changes backwards-compatible), **never push**.  `runs/`, `tasks/`, `TODO.md` are not committed.
Budget: one working day wall-clock.  **Resources:** `taskset -c 0-6,16-22` (7 physical cores + siblings), RSS `<= 45 GB`.
Detached launches (`setsid nohup … &`, kill by PID, never `pkill`); no `until … sleep` / `tail -f` / `pgrep -f` loops left behind.
**Don't hand back while your own jobs are running**: wait on them (one long Monitor/wait) and report only when done.

## Where it stands

`search/S32_COVER.md` (read §0 first).  Best cover: `runs/s32convR4_r11_last.txt`.  Strict protocol = `s21_cover_eval.py hardscan` + `family`.
At m = 6 halving the hardscan pitch (0.004 → 0.002) found a `17.75°` wall-band dip `0.3 %` deeper, and the closing
loop swaps holes round to round (at m = 5 a closed hole reopened `2.3 %` deep).  The exact zero-margin checker needs
`~0.25 %` margin over the *true* minimum, so the gap between strict-scan and true minimum must be `<~ 0.4 %`.

## What to do

1. Strict protocol at **pitch 0.002** (and family) on the current best cover.  Record the honest cost.
2. Closing loop (`rung2_close.py --near-tile`, fine tilted angles) whose separation includes the pitch-0.002 hardscan
   dips (and their polished neighbourhoods), with every dip pose kept as a **permanent** row so holes can't reopen.
   Add what's needed to the tooling (e.g. a `--rows-from` pool, or hardscan output fed back each round).
3. Stop when the strict minimum at pitch 0.002 is stable (changes `< 0.1 %`) over 3 consecutive rounds, or at budget.
   Then a **denser confirmation scan** on the final cover (pitch 0.001 in the wall bands / near-tile / tilted bands where
   dips have appeared) to estimate the remaining scan-vs-true gap.
4. Write the candidate for the exact check: the final cover scaled to `1.004 ×` over its strict minimum
   (`scale_cover.py`), as `runs/s32-close_candidate.txt`, with its total.  **Do not run the exact checker** (next task, the user decides).

## Report

Append a dated section to `search/S32_COVER.md`: trajectory table (round, total, strict min p=0.004 / 0.002, cost), where the
dips were, confirmation-scan result, candidate total vs `32` and the margin left.  Verdict: go / no-go for the exact
checker, and its estimated factor budget.  Commit (not pushed).  Final reply <= 200 words, numbers first.
