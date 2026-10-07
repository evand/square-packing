# s13-writeup: publish `s(13) = 4` case-free (2026-09-12)

**Why.**  Rung 2 is finished and trust-complete: one weighted closed cover of `[0,4]²`
(`certificates/rung2/s13_closed_cover_4.txt`, 3,621 points, `W = 2591194431/200000000 = 12.955972 < 13`),
verified exhaustively at margin zero by two checkers that share nothing — `search/zeromargin.py`
(Python `Fraction`, 1.7 h, `search/RUNG2.md`) and `verify2/zmcheck` (`i128` Rust, full domain, 17 min
on 8 threads, `search/RUNG2_XCHECK.md`) — with 23 rejection tests (`tests/rung2/rejection_tests.sh`)
and every soundness lemma in Lean (`lean/Sqpack/ZeroMargin.lean`, 36 theorems, standard axioms only,
`notes/lean-zeromargin.md`).  `TODO.md` item 1 of the critical path is to write it up so it is
checkable by a stranger the way the `s(12)` bound is.  Nothing mathematical is left to do here; this
task is documentation, `verify.sh`, CI and the Pages site.  Read first: `README.md` (the two
2026-09-12 paragraphs near the end of the "Beyond the ceiling" section are the seed), `search/RUNG2.md`
§0–§2 and §9, `search/RUNG2_XCHECK.md` §0–§1, `notes/review-2026-09-12.md`, `notes/lean-zeromargin.md`,
`VERIFICATION.md`, `docs/index.html`.

**Deliverables** (commit on your worktree branch; do not edit `TODO.md` or other tasks' files; do not
touch `certificates/rung2/*`, `verify2/src`, `search/zeromargin.py`, `lean/` — if you believe one of
those needs a change, write it down in the note instead):

1. **`README.md`**: a proper section `## s(13) = 4 without case analysis` placed after the `s(12)`
   material and before "Credits" — theorem statement (the cover, closed semantics, the `< 13`
   arithmetic, why `16` tiling squares make the bound sharp), the two checkers as a two-row table
   (what each is, language, domain reduced vs full, time, leaf census, determinism), the rejection
   tests, the Lean (which lemmas, 0 sorries, axioms), Theorem 1 of `RUNG2.md` as the one-paragraph
   explanation of why a proof of a `< 16` cover has to be disjunctive, and what is *not* machine-checked
   (the subdivision/exhaustiveness code of each checker is not in Lean — say it plainly).  Fold the
   existing 2026-09-12 paragraphs into it; keep the README's voice (flat, specific, numbers with
   sources).  Update the one-line summary near the top and the `What is verified, and how` table with
   a `verify2/` row.  Add Bentz 2010 as the result being re-proved in "References" if not already
   sufficient.
2. **`VERIFICATION.md`**: a dated `## s(13) = 4 (rung 2)` entry in the existing log style — files,
   hashes, both invocations verbatim, both censuses, the rejection suite count, the Lean build line
   (`runs/lean_build_2026-09-12.log`, `runs/zmcheck_main_2026-09-12.log`,
   `runs/rung2_rejection_2026-09-12.log` have the numbers).
3. **`verify.sh`**: build `verify2` and run `zmcheck cert certificates/rung2/s13_closed_cover_4.txt
   --depth 18 --threads "$(nproc)"`, checking for the `VERIFIED:` line exactly as `chk` does for the
   Rust verifier (a `PARTIAL SWEEP` or `NOT VERIFIED` must fail the script); then
   `./tests/rung2/rejection_tests.sh`.  Add the rung-2 certificate to `certificates/SHA256SUMS`.  Leave
   the 1.7 h `zeromargin.py` sweep as a commented "slow path" command next to it, the way the
   `xcheck.py` N=6000 runs are.  **Measure** the `zmcheck` wall time at 4 threads (GitHub runners have
   4 cores) — if it is over ~45 min, say so in the note and put the full sweep behind
   `workflow_dispatch`/the monthly schedule in CI while every push still runs the rejection suite plus
   a restricted-band sweep; if it is under, run it on every push.  Extend `.github/workflows/verify.yml`
   accordingly (cache `verify2/target` too).  Run `./verify.sh` end to end once, detached
   (`setsid nohup ./verify.sh > runs/verify_full_$(date +%F).log 2>&1 &`), and quote its tail.
4. **`docs/index.html`** (the GitHub Pages write-up): a section for `s(13) = 4` mirroring the README
   section at reader level — what a closed cover is, why weight `< 13` suffices, what "margin zero"
   and "exhaustive" mean, the two-checker table, a link to `RUNG2.md`.  Match the page's existing
   style and structure; no new JS unless the page already has a pattern for it.  Do not draw the
   3,621 points unless the page already renders certificates from JSON and it is a one-line addition.
5. **`notes/s13-casefree.md`**: a self-contained 2–4 page note in the shape of a short paper
   section: statement, semantics (closed squares, closed containment, `t = 4`, the dilation argument
   for why closed is the right semantics — `notes/proof-anatomy.md` §1 and §7), the certificate, the
   checking method in one page (the primitives, the disjunction, Theorem 1), the checkers, the Lean,
   the rejection tests, reproduce commands, and a "what would make this wrong" paragraph.  Cite
   Bentz 2010 and DS7 as in the README.  This is the piece that could be sent to someone.

**Do not**: quote any `n = 12` `t = 4` packing-side number as a bound (see `TODO.md`, `HONEST.md`);
touch pids 1109693–1109695 or anything under `.claude/worktrees/`; run anything longer than 10 min
in the foreground; use more than 16 threads.  Run `sha256sum -c certificates/SHA256SUMS` before and
after.  Every number you write must come from a file in the repo or a run you made — name it.

**Finish** with `search/S13_WRITEUP.md`: what was changed (file by file), the `verify.sh` timing
table (8 threads here, 4 threads), the CI decision and why, and anything you found inconsistent
between `RUNG2.md`, `RUNG2_XCHECK.md`, `ZEROMARGIN.md`, the Lean note and the README (list, don't
silently fix outside the files above).  Budget: ~3 h plus the detached `verify.sh` run.
