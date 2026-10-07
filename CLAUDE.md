# square-packing: working notes for agents

Start at `TODO.md` (open items; conventions in `~/math/TODO.md`), then `README.md` (every result, its bundle, how to re-check it).
`Completed.md` is the done log (one line each, newest first).

## Layout
- `certificates/<result>/`: published bundles (cover, run records, `SHA256SUMS`, `verify.sh`, README).  Files listed in a
  `SHA256SUMS` are frozen: never edit them; add a new file instead.  READMEs and `verify.sh` are not hashed.
- `verify/` (Rust `verify`), `verify2/` (`zmcheck`, `zmx2`), `xcheck.py`, `search/zeromargin.py`, `search/zm_mixed.py`: checkers.
  Exit status of `verify`/`zmcheck`: 0 VERIFIED, 1 NOT VERIFIED, 2 bad input, 3 internal, 4 partial; `zmx2` still exits 0 (read its verdict line).
- `lean/`: Lean 4 + Mathlib (`cd lean && lake build`; `lake env lean Axioms.lean` must show no `sorryAx`).  Big generated data dirs are gitignored (`lean/scripts/gen_data.sh`).
- `search/`: research code + one `UPPERCASE.md` log per investigation.  `search/packer/` (packing search; `PACKER.md`, `README.md`), `search/exact/` (exact forms; `EXACT_FORMS.md`).
- `docs/`: HTML write-ups (Pages serves them at /s12/, /s13/, /s21/, ... via `.github/workflows/pages.yml`); `site/`: the Atlas (`site/www/` is the site root).
- `notes/`, `tasks/<task>/`: working notes and per-task material (history; dated files are not rewritten).
- `runs/`: big local outputs, gitignored.  `outreach` and `notes/publish-inventory.md` are symlinks into `../private/` (gitignored; never commit them).
- `s12/README.md` is only a signpost: everything lived under `s12/` until 2026-10-07.

## Checking
- `./verify.sh` (fast tier, ~15 min on 16 cores; CI runs it on every push touching code/certificates) and `./verify.sh --full` (slow sweeps).
- Rejection suites: `tests/rejection_tests.sh`, `tests/rung2/rejection_tests.sh`.  Exact batch: `search/exact/batch/verify_all.sh`.
- Run scripts from the repo root unless they `cd` themselves.

## Rules
- External posts (GitHub issues/comments on other repos, X, email): draft in `outreach/`, ask Evan explicitly per post.
- Before any bound-chasing or record claim, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`, `packing/frontier/RESULTS.md`).
- `~/math/_untrusted-third-party/`: read only, never execute.
- Work on `main` (single checkout, no long-lived worktrees); branch for big restructures.  Commit messages end with the attribution lines the session gives.
