# s60-bundle: ship s(60) = 8 as certificates/s60/, wire it in, commit (2026-09-28)

**What.**  `s(60) = 8` is certified (`search/S60_COVER.md`, commit 86cb2d8): cover `runs/s60_mixed_candidate_8.txt`
(sha256 `2d0e456e…`), zmx2 --d4/--full VERIFIED, zm_mixed.py --d4 --cert-mode VERIFIED-D4 (records in
`runs/zm_mixed_s60/`).  Bundle and wire it in exactly as s(45) was: commits d200158 (bundle) and 5d8bcd0 (wiring)
are the template — read both diffs first and mirror them.

**Steps.**
1. `search/s60_cert_runs.sh` = `search/s45_cert_runs.sh` for s60 (paths, root counts 6,400 / 51,200 / 102,400).
2. `certificates/s60/`: `s60_mixed_cover_8.txt` (byte copy of the candidate; check sha), `README.md` in the s45
   style (numbers from S60_COVER.md §0/§2; same trust caveats; "no Lean"), `SHA256SUMS`, `verify.sh`.
3. `zm_mixed_d4/` by **`import runs/zm_mixed_s60`** (no re-sweep: the 19.6 CPU-h run is not to be repeated).
   `zmx2_d4/`, `zmx2_full/`: fresh runs on the bundled file (≈ 3–4 CPU-min).  Use `CORES=0-6`, at most 7 processes.
4. Wiring as in 5d8bcd0: `s12/verify.sh` (+ `S60_SEPARATE` for the zm_mixed re-sweep), `.github/workflows/verify.yml`,
   `site/data/lower_bounds.json` and `site/www/data/lower_bounds.json` (n = 60 best = 8, history: find the prior
   published lower bound for s(60) as the s45 commit did for wand125's s(45) — Friedman's page / the existing json;
   don't invent one), `site/www/sources.html`, `README.md`, `s12/README.md`.  If the `/s45/` write-up page presents
   s(45) as a rung of the k²−4 family, add s(60) there in a sentence or a table row; no new page.
5. Run `certificates/s60/verify.sh` (default tier) and `s12/verify.sh` fast tier; both must pass.

**Don't.**  Push; edit `TODO.md`, `lean/`, `search/*.py` checkers, or other certificate bundles; re-run zm_mixed.
**Done when.**  Two commits (bundle; wiring) on main in `/home/evand/math/square-packing/public`, message style as
the s45 ones, ending `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.  Report: shas, verify output, what
you changed on the site, anything that didn't mirror cleanly.
