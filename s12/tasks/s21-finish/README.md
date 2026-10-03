# s21-finish: from "VERIFIED-D4 by zm_mixed.py" to a shipped s(21) = 5, 2026-09-27

Candidate: `runs/line-cover_m5_candidate_x1003.txt` (mixed v1, `tasks/line-cover/FORMAT.md`), total
`522368729933/25000000000 = 20.89474919732 < 21`.  VERIFIED-D4 by `search/zm_mixed.py` (commit e72f317; `search/ZM_MIXED.md`
§7: 40,000 roots, 461,204 boxes, 0 uncertified).  Cover side: `search/LINE_COVER.md`.  Standard to meet: what `s(32)` got
(`search/S32_EXACT.md` §10–12, `certificates/s32/`, `notes/lean-s32.md`).

Three parallel agents, then a bundling step:
- **audit** (`audit.md`): adversarial review of zm_mixed.py's lemmas, code-vs-proof, the run's coverage.  Cores `10-11`.
- **xcheck** (`xcheck.md`): an independent second checker, written without reading zm_mixed.  Cores `0-9`.
- **lean** (`lean.md`): the reduction for mixed measures + D4 for measures + `s21_eq_five_of_checker`.  Cores `12-13`.
- then bundle `certificates/s21/`, verify.sh/CI, write-up (after audit + xcheck agree; Evan approves pushing).

**Machine:** 16 physical cores (logical `i` and `i+16` are siblings); hyperthreading does not help these workloads, so
pin to physical cores only, as listed.  Common rules: commit own files only, **never push**; don't modify
`search/zeromargin.py`, `certificates/s32/`, or zmcheck's existing behaviour; detached launches, kill by PID, wait for
your own jobs before replying.  Final replies `<= 200` words, numbers first.
