# Brief: the cover side of `s(45) = 7` (the `m = 7` rung of `n = m² − 4`), 2026-09-24

Working tree: `~/math/square-packing/public/s12`.  Commit there (`git add` your own files only; a sibling agent is
committing to the same repo on the `m = 6` rung), **never push**.  `runs/` is gitignored; `tasks/`, `TODO.md` are private
symlinks, don't commit them.  Budget: one working day wall-clock; stop and report if the approach hits a wall.
**Resources:** pin everything with `taskset -c 8-14,24-30` (7 physical cores + SMT siblings; `--nproc 7` or so), keep RSS
`<= 45 GB` total (a sibling run shares the 125 GB box; at `m = 5` an LP polish hit 81 GB and had to be killed).  Launch
anything over 10 minutes detached (`setsid nohup … &`, log to a file, kill by PID, never `pkill`).  Leave no
`until … sleep` / `tail -f` watchers behind (a `pgrep -f <name>` loop matches itself and never exits).

## Why

`s(45)` is open; the best known packing is the trivial `7`, and its neighbour `s(46) = 7` is proved (Bentz).
A weighted closed point cover of `[0,7]²` of total `< 45`, checked exactly at zero margin, proves `s(45) = 7` (the
`s(13) = 4` machinery, `search/RUNG2.md`).  At `m = 5` this failed narrowly: converged cover LP `20.73`, room
`21/20.73 = 1.013`, best honest cover `21.016` on the strict protocol (`search/S21_COVER.md`,
`S21_OVERHEAD.md`, `M4_MARGIN.md`; overview `notes/status.md` last paragraph).  The cell model of
`S21_COVER.md` §3 (corner ≈ 0.50, edge ≈ 0.80, interior ≈ 1.00 per cell) *estimates* the LP at `m²−4 − 1.7`,
i.e. room `~4 %`, against total overhead `~5 %` shipped at `m = 4` and a floor
`~1.0–1.4 %` observed so far (validity `~1–3 %`, checker margin `~0.25 %`).  This task measures whether the room is real.

## What to do

1. **Converged cover LP at `s = 7`** with `search/closed4.py run --s 7 --no-literature` and the `m = 5` recipe
   (`S21_COVER.md` §1, §5, §8): item-3 seed rows (`family_rows.py --s 7`) as a lazy pool, a near-tile pool (the
   `m²` tile centres, small tilts), `--warm-lp`, column/row pruning, restarts `--from-cert` as the effective column
   prune.  LP size grows ~area²; if the cold start is too slow, warm-start columns by transplanting the `m = 5`
   cover (`runs/s5convR3_final_last.txt`, `runs/s5convD3_it9_last.txt`): its corner/wall bands onto the new
   walls, its centre cell (weight ≈ 1.000) replicated across the interior.  Your call.  Check every tool for
   hard-coded `s = 4`/`5` (`s21_cover_eval.py`, `rung2_close.py --cap`, `scale_cover.py`) and fix backwards-compatibly.
2. **Honest cost**: closing loop (`rung2_close.py --near-tile`, `--cap` set for `s = 7`) and the **strict**
   protocol (`s21_cover_eval.py hardscan` + `family`) on the best covers.  Report total / min capture.
3. **Where the weight sits** (`s21_cover_eval.py regions`): per-cell table, checks the cell model.

Optional, only if the LP is clearly below `45`: the certified `nu_f(7)` lower bound (`S21_KILL.md`) is *not*
needed. The converged cover LP is heuristically `>= nu_f`, so an LP `< 45` already means the kill test cannot fire.
Don't start the exact zero-margin checker; that is the next task and depends on the answer.

## Report

`search/S45_COVER.md` in the style of `S21_COVER.md` (labels **certified / heuristic / estimate**; §0 answer
first; runs table; honest-cost table with the `m = 5` rows for comparison; cell table; reproduce block).  Headline:
converged LP, room `45 − LP` in absolute and %, best strict honest cost, and a verdict: is `W < 45` plausible with
the known overhead components, and what factor would the exact checker have to hit?  Commit (not pushed).  Final
reply: <= 200 words, numbers first.
