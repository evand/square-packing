# Writing up `s(13) = 4` (task `s13-writeup`, 2026-09-12)

Documentation, `verify.sh`, CI and the Pages site for the rung-2 result.  No mathematics was done
in this task and no mathematical claim is new here: every number below comes from a file already in
the repository or from a run made in this task and named as such.  Nothing under
`certificates/rung2/`, `verify2/src`, `search/zeromargin.py` or `lean/` was touched.

## 1. What changed, file by file

| file | change |
|---|---|
| `README.md` | New section **`## s(13) = 4 without case analysis`**, placed after the `s(12)` material and before "Credits and prior art": theorem statement, the rescaling arithmetic in full, why closed semantics is the content rather than a convenience, `RUNG2.md` Theorem 1 as the one-paragraph reason the proof has to be disjunctive, the two checkers as a table, the 23 rejection tests, the Lean, what is **not** machine-checked, provenance and the `COVER4` floor, and reproduce commands.  The 2026-09-12 paragraph in "Beyond the ceiling" is folded into it and replaced by a one-line pointer.  Top-of-file summary gained a paragraph on the second result with a link to the section.  "What is verified, and how" gained a `verify2/` row.  The Bentz 2010 reference now says it is the result being re-proved, and that the theorem is his and only the proof is new. |
| `VERIFICATION.md` | New dated entry **`## s(13) = 4 (rung 2) (2026-09-12)`** in the existing log style: file, container, point count, total weight, sha256; both checkers' invocations verbatim with their censuses; the note that `verify/` provably *cannot* check this file; the rejection suite; the Lean build line; the re-timing done for this task; what changed in `verify.sh`.  The "Result" section at the top names the second result and points at the entry. |
| `verify.sh` | Builds `verify2` alongside `verify`.  New `chkzm()` with the same contract as `chk()` — `zmcheck` exits 0 on both `NOT VERIFIED` and `PARTIAL SWEEP`, so both must fail the script.  New block running `zmcheck cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads "$(nproc)"`, then `./tests/rung2/rejection_tests.sh`.  The 1 h 42 min `zeromargin.py` sweep sits next to it as a commented slow path with its expected census, the way the `xcheck.py` N=6000 runs already did. |
| `certificates/SHA256SUMS` | `certificates/rung2/s13_closed_cover_4.txt` added: `ea303acea08cc17a13cecc24d3714c2df409f91eba048cd5546050ed064b53ed`.  `sha256sum -c` is clean before and after (36 files). |
| `.github/workflows/verify.yml` | `verify2/target` added to the cargo cache, and `verify2/Cargo.lock` to the cache key.  `timeout-minutes: 300` so a hang fails loudly instead of burning the 6 h default.  A second `sha256sum -c` after `./verify.sh`, so a certificate silently changing under the run cannot pass.  A comment records the timing measurement and names the exact fallback if a GitHub runner turns out to be too slow. |
| `docs/index.html` | New section **`07 · A second result — s(13) = 4, without case analysis`**, between "Limits" and "Reproduce" (which becomes `08`), at reader level: the theorem in a `card thm`, what a closed cover is and why weight under 13 suffices, what "margin zero" and "exhaustive" mean (including why the section-05 verifier provably cannot do this job), the two-checker table, the `zmcheck` transcript, and links to `RUNG2.md`, `RUNG2_XCHECK.md` and `notes/s13-casefree.md`.  Existing style and structure only; no new JS, no new figure.  The artifacts table gained the rung-2 rows and two dead paths were corrected (see §4). |
| `notes/s13-casefree.md` | **New.**  The self-contained note, in the shape of a short paper section: statement and the corollary proved in three lines; semantics, with the dilation argument and the open/closed contrast; the certificate and its provenance; the method and Theorem 1 in one page, with the measured shape of the required disjunction; the two checkers; the Lean; the rejection tests; reproduce commands; and a "what would make this wrong" section listing six failure modes in decreasing order of concern.  Cites Bentz 2010 and DS7. |
| `search/S13_WRITEUP.md` | **New.**  This file. |

## 2. Timing

All on the development machine (32 cores, load from three unrelated `python3` jobs at 100 %
throughout), the same binary and the same file each time:
`verify2/target/release/zmcheck cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads T`.

| threads | wall time | boxes | max depth | leaf census | log |
|---|---|---|---|---|---|
| **8** | **1064 s** (17 min 44 s) | 30258 | 10 | `ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0` | `runs/zmcheck_t8_2026-09-12.log` |
| **4** | **2015 s** (33 min 35 s) | 30258 | 10 | `ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0` | `runs/zmcheck_t4_2026-09-12.log` |
| 8 (earlier, main checkout) | 1065 s | 30258 | 10 | same | `runs/zmcheck_main_2026-09-12.log` |
| 8 (earlier, worktree; `RUNG2_XCHECK.md` §0) | 1017 s | 30258 | 10 | same | — |

The two rows measured here were run **concurrently** (8 + 4 = 12 threads on a 32-core machine that
also had three unrelated `python3` jobs pinned at 100 % throughout), so both numbers are if
anything pessimistic; the 8-thread figure came out at `1064 s` against `1065 s` for the earlier
solo run in the main checkout, which says the contention cost was negligible.  Speed-up from 4 to 8
threads: `2015/1064 = 1.89×`.

The census is identical in all four runs, which is the point: the verdict, the box count, the max
depth and the leaf breakdown do not depend on the thread count.

Full `./verify.sh`, end to end, detached and capped at 16 cores
(`setsid nohup taskset -c 0-15 ./verify.sh > runs/verify_full_2026-09-12.log 2>&1 &`; the `taskset`
is because `verify.sh` passes `$(nproc)` and this machine has 32).  **1054 s (17 min 34 s)**,
21:18:28 → 21:36:02, exit 0, `runs/verify_full_2026-09-12.log`:

* **25 `VERIFIED` verdicts** (was 24 before this task; the new one is the rung-2 sweep), no
  `NOT VERIFIED`, no `PARTIAL SWEEP`, no `REJECTED`;
* the rung-2 sweep inside it: `done in 626s: boxes 30258, max depth 10 / ADM 5114  DISJ 9477
  EMPTY 6938  UNCERTIFIED 0 / VERIFIED` — the same census again, at 16 threads;
* `23 passed, 0 failed, 0 panics` (rung-2 suite, with the `closed4_best_x103` file present) and
  `172 passed, 0 failed, 0 panics` (the main suite);
* every `points.json` round-trip byte-identical.

`sha256sum -c certificates/SHA256SUMS` is clean before and after: 36 of 36 OK.

## 3. The CI decision

**The full sweep runs on every push.**  The measurement is `2015 s = 33.6 min` at 4 threads, the
GitHub runner's core count — inside the `~45 min` budget the brief set, with about 25 % of headroom.
So `verify.sh` runs it unconditionally, `.github/workflows/verify.yml` keeps its existing triggers
(push, pull_request, workflow_dispatch, monthly schedule), and nothing is moved behind
`workflow_dispatch`.  The only workflow changes are mechanical: `verify2/target` added to the cargo
cache and `verify2/Cargo.lock` to its key, `timeout-minutes: 300`, and a second `sha256sum -c`
after the run.

Three things a future reader should know about that decision:

* **A GitHub core is slower than this machine's.**  The 25 % headroom is measured here, not there,
  and the runner also has to do everything else `verify.sh` does (the `N = 6000` and `N = 12000`
  sweeps of the `s(12)` certificates, the uniform family, `xcheck.py`, 172 rejection checks) on the
  same 4 cores.  If the job starts hitting the timeout, the fix is written down in a comment at the
  top of the job and is the brief's fallback: move the full sweep to `workflow_dispatch` plus the
  monthly `cron`, and give every push `tests/rung2/rejection_tests.sh` — which already contains a
  restricted-band sweep over `c_x ∈ [0.5,0.6]`, `c_y ∈ [1,2]`, the container's hardest band — in
  its place.  `zmcheck` already supports that: `--xlo/--xhi/--ylo/--yhi`.  Nothing needs to be
  built to take that route.
* **The restricted-band fallback is not a weak substitute.**  It is the band containing
  `(½, 3/2, 0)`, which `RUNG2.md` §2 identifies as the worst monotone-witness pose of the
  container, and it is where every mutation test is discriminated.  What it does *not* do is print
  `VERIFIED` — a restricted sweep always prints `PARTIAL SWEEP` and then `NOT VERIFIED` — so under
  that route the green badge would mean "the checker still refuses everything it should, and the
  hard band is still clean", not "the theorem is re-verified".  That is the reason to prefer the
  full sweep while it fits.
* **One rung-2 rejection test skips in CI.**  `tests/rung2/rejection_tests.sh` case (e) needs
  `runs/inputs-2026-09-11/closed4_best_x103.txt`, a development-history file kept outside the
  repository (`runs/` is `.gitignore`d).  The script prints `skip` and reports `20 passed` instead
  of 23 when it is absent, and does not fail.  That is deliberate but undocumented outside the
  script; it means the CI log will not match the 23 quoted in `README.md` and `VERIFICATION.md`.

## 4. Inconsistencies found between the documents

Listed, not silently fixed, except where the file is one of this task's deliverables (`README.md`,
`VERIFICATION.md`, `docs/index.html`) — those cases are marked **[fixed here]**.

1. **`docs/index.html` §06 "Limits" contradicts `README.md` "Limits of the method" on the ceiling.**
   The page says `ν_f` "crosses 12 at roughly `s ≈ 3.95`", that "the ceiling of this entire family
   of arguments for n = 12 is about `3.95–3.96` — not 4", and that "the bound proved here is within
   about `0.02` of that ceiling".  `README.md` says the ceiling is now **pinned**: an exact
   fractional packing of mass `12.0282` at `s = 399/100` puts `s*` in `[3.968616, 3.99)`.  The
   shipped bound `3.968616` is the *lower* end of that bracket, so "within about 0.02" is not
   right either.  The page predates the pinning (`search/DUAL_EXACT.md`) and needs rewriting
   against it.  **Not fixed** — it is a substantive rewrite of a section this task's brief does not
   cover, and it needs the owner's reading of `search/DUAL.md` and `search/CLIQUE_CONTINUUM.md` §3.
2. **`README.md` "What is verified, and how" says `tests/` has "42 rejection tests".**
   `VERIFICATION.md` (2026-09-08) says the suite is **172 checks** (was 138, was 42); the script
   itself has 115 `check` invocations.  `docs/index.html` also said 42.  **[fixed here]** in the
   Pages artifacts table by dropping the stale count; the `README.md` table row is left for the
   owner, since the correct number depends on whether "checks" or "cases" is meant.
3. **`docs/index.html` artifacts table listed two paths that do not exist**: `cert5.py` (now
   `search/lp_search.py`) and `pack/` (now `search/pack_src/`).  **[fixed here]**, since the same
   table needed the rung-2 rows.
4. **`search/ZEROMARGIN.md` §2 says "a box is a leaf when one of *four* tests succeeds"** and lists
   `EMPTY`, `CORE`, `P1`, `TRI`.  `search/zeromargin.py` now also has `ADM`, `MIX` and `CHAIN`
   (`RUNG2.md` §3, §6), and its census prints all seven.  §2 is the reader's first description of
   the checker and is stale by three primitives.  `RUNG2.md` §3 does say `ADM` "replaces both"
   `CORE` and `P1`, but nothing in `ZEROMARGIN.md` points forward to it.
5. **`search/ZEROMARGIN.md` §2 describes the root grid as "eight `u`-bins on `[0, 1/2]`"**
   unconditionally.  That is the *symmetry-reduced* domain; `RUNG2.md` §0 records that the checker
   verifies the invariance exactly and refuses to run reduced otherwise, and `ZEROMARGIN.md` §3
   itself mentions a `--full` mode.  §2 never says the reduction is conditional or that it is
   checked, which is the single most load-bearing caveat about that checker's result.
6. **The `zmcheck` wall time is quoted as three different numbers.**  `RUNG2_XCHECK.md` §0 and
   `RUNG2.md` §0 both say `1017 s` ("17 minutes"); `notes/review-2026-09-12.md` records the
   main-checkout reproduction at `1065 s`; `tests/rung2/rejection_tests.sh`'s header comment says
   the full sweep "takes ~20 min on 8 threads".  All three are real measurements of the same run on
   the same machine, so nothing is wrong — but no single document says which one a reader should
   expect, and the `~20 min` in the test header is the only one with no log behind it.  This task's
   own measurements (§2) are a fourth set.
7. **`certificates/rung2/README.md` does not mention the second checker at all.**  Its "How it was
   verified" section gives only the `zeromargin.py` invocation and census, and its "Independent
   checks" list is `zeromargin_stress.py` and `closed4.py stress`.  `verify2/zmcheck` — the full
   domain run, the one thing that makes the result trust-complete — is absent, as are the 23
   rejection tests and the Lean.  A reader who opens the certificate directory gets the weaker
   story.  **Not fixed**: `certificates/rung2/*` is out of bounds for this task.
8. **`notes/lean-zeromargin.md` §6 item 2 is already resolved upstream.**  It records that
   `RUNG2.md` §6.2's Lemma E should say "real combination" rather than "nonnegative"; `RUNG2.md`
   §6.2 now reads "Every real combination".  The discrepancy list does not say it was fixed, so a
   reader is left checking.  Items 1, 3, 4 and 5 of that section are live and correctly describe
   the current files.
9. **`README.md`'s `verify/` row and the new `verify2/` row overlap in a way worth watching.**
   Nothing in the repository states, in one place outside `ZEROMARGIN.md` §1, that `verify/`
   *cannot* check a cover at the container itself — the angle-net erosion is linear against a
   quadratic margin, so no `N` works.  Added to the `VERIFICATION.md` entry and the Pages section
   by this task; the `README.md` table row for `verify/` still reads as though it were the general
   tool.

None of the above touches the correctness of the `s(13) = 4` result or of any number in it.
