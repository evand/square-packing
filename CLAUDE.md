# square-packing: working notes for agents

Start at `TODO.md` (open items; conventions in `~/math/TODO.md`), then `README.md` (every result, its bundle, how to re-check it).
`Completed.md` is the done log (one line each, newest first).

## Layout
- `certificates/<result>/`: published bundles (cover, run records, `SHA256SUMS`, `verify.sh`, README).  Files listed in a
  `SHA256SUMS` are frozen: never edit them; add a new file instead.  READMEs and `verify.sh` are not hashed.
- `verify/` (Rust `verify`), `verify2/` (`zmcheck`, `zmx2`), `xcheck.py`, `search/zeromargin.py`, `search/zm_mixed.py`: checkers.
  Exit status of `verify`/`zmcheck`: 0 VERIFIED, 1 NOT VERIFIED, 2 bad input, 3 internal, 4 partial; `zmx2` still exits 0 (read its verdict line).
- `lean/`: Lean 4 + Mathlib (`cd lean && lake build`; `lake env lean Axioms.lean` must show no `sorryAx`).  Big generated data dirs are gitignored (`lean/scripts/gen_data.sh`).
- `search/`: research code + one `UPPERCASE.md` log per investigation.  `search/packer/` (packing search; `PACKER.md`, `README.md`, `NEEDS.md`, `DESIGN.md`), `search/exact/` (exact forms; `EXACT_FORMS.md`).
- Packings: use the store (`search/packer/pk.py`: `frontier`, `ls N`, `get best:N --fmt ...`, `sync-register`, `watch --apply`, `sync-issues`; `pk/README.md`), never copy coordinates by hand.  Add verified packings to the store (ours or others') without asking: a packing is checkable.  Provenance and "who was first" are best effort (socially verified): record what we know, fix it when we learn better (Evan, 10-12).  There is no separate records session, only notes on how to do it.  Grade new search ideas with `pk/battery.py` against the frozen baseline.
- Launch background compute with `python3 search/packer/pk/jobs.py start NAME --procs N -- CMD` (`jobs.py free` suggests N; it caps math-library threads at 1 per process, since an in-process numpy / exactsolve pool otherwise runs ~47 threads per worker: load ~300 on 10-14).  Check `jobs.py status` after launch: it flags oversubscription.
- New packer launchers: call `pk.manifest.start(out, name, ...)` / `finish()` (provenance + measured CPU) and `scratch.scope(name)` (temp dirs); report compute as measured CPU (`cpu_s`), not wall x procs.  Review + pre-run fix list: `search/packer/REVIEW-2026-10-11.md`.  Strategy-steering numbers live in `search/packer/CLAIMS.md` (cite IDs, update status there); open anomalies in `search/packer/CONFUSIONS.md`; pre-register any test run (hypothesis, metric, rule) in PACKER.md before launch.
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
- Before any bound-chasing or record claim, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`, `packing/frontier/RESULTS.md`); the local clone's `frontier/n-<n>.md` is in `~/math/_untrusted-third-party/jlevy-squares` (refresh per its PROVENANCE.txt).
- `site/www/problems.html` is the reader-facing copy of `search/WISHLIST.md`'s object-level items: edit both together.
- `~/math/_untrusted-third-party/`: read only, never execute.
- Licensing (Evan, 10-17): never copy others' code, but packings / certificates are mathematical facts, not
  copyrightable: use them freely (store, Lean, gluing), like a result from a paper, with careful, courteous citation.
- Work on `main` (single checkout, no long-lived worktrees); branch for big restructures.  Commit messages end with the attribution lines the session gives.
