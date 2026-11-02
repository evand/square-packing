# square-packing: working notes for agents

Start at `TODO.md` (open items; conventions in `~/math/TODO.md`), then `README.md` (every result, its bundle, how to re-check it).
`Completed.md` is the done log (one line each, newest first).

## Layout
- `certificates/<result>/`: published bundles (cover, run records, `SHA256SUMS`, `verify.sh`, README).  Files listed in a
  `SHA256SUMS` are frozen: never edit them; add a new file instead.  READMEs and `verify.sh` are not hashed.
- `verify/` (Rust `verify`), `verify2/` (`zmcheck`, `zmx2`), `xcheck.py`, `search/zeromargin.py`, `search/zm_mixed.py`: checkers.
  Exit status of `verify`/`zmcheck`: 0 VERIFIED, 1 NOT VERIFIED, 2 bad input, 3 internal, 4 partial; `zmx2` still exits 0 (read its verdict line).
- `lean/`: Lean 4 + Mathlib (`cd lean && lake build`; `lake env lean Axioms.lean` must show no `sorryAx`).  Big generated data dirs are gitignored (`lean/scripts/gen_data.sh`).
- **Public / private (10-23).**  This repo is public.  New search methods (moves, generators, constructors, bandits,
  priors, acceptance rules), their configs, notes, TODO items and results for unposted packings go in the private plugin
  repo `../plugins` (`search/packer/DESIGN.md` "Plugins"; run its scripts as `search/packer/priv/<script>`), never here.
  Framework changes go here as generic extension points.  Bug fixes and minor updates to existing public moves /
  generators go public too; major new methods are a judgment call (Evan 10-23).  `.git/hooks/pre-push` refuses private paths and added lines
  that mention Kaggle (review them, then `PUBLIC_OK=1 git push`).  Unsure which side: ask Evan.
  Public code stays generic (Evan 10-25): never name a plugin's stages, arms or modules here (code, docstrings, TODO,
  examples), and keep per-plugin specifics (weights, priors, benchmarks, configs) in the plugin; public examples use
  public stages.  The pre-push hook refuses added lines naming a plugin stage (`@stage('x.y')`) or a word in
  `../plugins/PUBLIC_DENY.txt` (add new bare arm names there); `PLUGIN_OK=1` only for a reviewed exception.  The hook
  checks the push's net diff but history keeps every commit: fix a slip before it is pushed, ideally before committing.
- `search/`: research code + one `UPPERCASE.md` log per investigation.  `search/packer/` (packing search; `README.md` map, `DESIGN.md` architecture, `CLAIMS.md` findings; lab log to 10-24 in `archive/`), `search/exact/` (exact forms; `EXACT_FORMS.md`).
- Packings: use the store (`search/packer/pk.py`: `frontier`, `ls N`, `get best:N --fmt ...`, `lineage ID` (credit start, every parent link kind), `sync-register`, `watch --apply`, `sync-issues`; `pk/README.md`), never copy coordinates by hand.  Add verified packings to the store (ours or others') without asking: a packing is checkable.  Provenance and "who was first" are best effort (socially verified): record what we know, fix it when we learn better (Evan, 10-12).  There is no separate records session, only notes on how to do it.  Packing search is the event loop + stages (`search/packer/DESIGN.md`): a new idea is a stage used as a move, a campaign is a campaign file run by `python3 -m pk.events`; the file is the pre-registration and holds the result, so it lives in the plugin repo (`../plugins/packer/goals/`), not here.
  Explorer archive files (`b#####.txt`) may be absent (deduplicated against the store, 10-20): read them with `pk.packing.read` / `mcmin.load_deg`, or call
  `pk.archive.ensure(path)` before handing a path to fq; `pk.py archive-restore RUN` rehydrates a run.  Store dedupe (`pk.py dedupe`): never `--apply` without
  discussing the rule with Evan in that session.
- Short work (units of 10-60 s, a few minutes at most) goes to the pool: `python3 search/packer/pk/tasks.py add CAMPAIGN --prio high|normal|low -- CMD` or `repeat CAMPAIGN --budget CPU_H -- CMD`; no core counts (`tasks.py -h`).  Renders and other ordinary scripts run outside it; it makes room.
- Launch long runs that cannot be split with `python3 search/packer/pk/jobs.py start NAME --procs N -- CMD` (`jobs.py free` suggests N; it caps math-library threads at 1 per process, since an in-process numpy / exactsolve pool otherwise runs ~47 threads per worker: load ~300 on 10-14).  Check `jobs.py status` after launch: it flags oversubscription.
- Strategy-steering numbers and open confusions live in `search/packer/CLAIMS.md` (cite IDs, update status there).
- `docs/`: HTML write-ups (Pages serves them at /s12/, /s13/, /s21/, ... via `.github/workflows/pages.yml`); `site/`: the Atlas (`site/www/` is the site root).
- `literature/`: literature scans and reading lists; `literature/README.md` indexes every scan, including older ones elsewhere.  Save new lit-search results here and add a row to the index.
- `notes/`, `tasks/<task>/`: working notes and per-task material (history; dated files are not rewritten).
- `runs/`: big local outputs, gitignored.  `outreach` and `notes/publish-inventory.md` are symlinks into `../private/` (gitignored; never commit them).
- `s12/README.md` is only a signpost: everything lived under `s12/` until 2026-10-07.

## Checking
- Tiers: `./verify.sh --spot` (~25 s; what CI runs on push), `./verify.sh` (fast: + fresh sweeps, rejection suites; ~4 min on 16 cores; run it after touching a checker or certificate), `./verify.sh --full` (slow sweeps).
- Lean: CI builds only the spot target (`lean/AxiomsSpot.lean`); `cd lean && lake build && lake env lean Axioms.lean` locally (`Sqpack.Bentz` peaks at 23.5 GB).  CI is a courtesy spot check, not a reproduction; don't engineer it for hosted-runner limits.
- Rejection suites: `tests/rejection_tests.sh`, `tests/rung2/rejection_tests.sh`.  Exact batch: `search/exact/batch/verify_all.sh`.
- Run scripts from the repo root unless they `cd` themselves.

## Rules
- External posts (GitHub issues/comments on other repos, X, email): draft in `outreach/`, ask Evan explicitly per post.
  Record submissions: follow `search/packer/SUBMITTING.md` (vocabulary: *refinement* only for the same local minimum,
  else *a new packing found from X*; credit; template; checklist); draft with `search/packer/submission.py NAME`.
- Beyond the register (Evan 10-23): public packings not in jlevy's register get posted with credit to their finder; where
  only a side length is public (Kaggle per-puzzle list), we may claim "first publicly known packing" with details, naming
  the side-only claim.  `search/packer/regcheck.py` lists these cases; policy in `search/packer/SUBMITTING.md`.
- Before any bound-chasing or record claim, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`, `packing/frontier/RESULTS.md`); the current register is the store's sparse clone `search/packer/runs/jlevy` (`packing/frontier/n-<n>.md`; refresh with `pk.py sync-register`), not the older read-only copy in `~/math/_untrusted-third-party/jlevy-squares`.
- `site/www/problems.html` is the reader-facing copy of `search/WISHLIST.md`'s object-level items: edit both together.
- `~/math/_untrusted-third-party/`: read only, never execute.
- Licensing (Evan, 10-17): never copy others' code, but packings / certificates are mathematical facts, not
  copyrightable: use them freely (store, Lean, gluing), like a result from a paper, with careful, courteous citation.
- Work on `main` (single checkout, no long-lived worktrees); branch for big restructures.  Commit messages end with the attribution lines the session gives.
- Commits and pushes (Evan 10-22): commit locally whenever it makes sense, in good logical groupings, without asking; always ask before pushing to any remote.
