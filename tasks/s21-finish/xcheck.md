# Brief: an independent second exact checker for mixed covers, 2026-09-27

Goal: re-verify `runs/line-cover_m5_candidate_x1003.txt` (the `s(21) = 5` cover) with a checker that shares no code and no
lemma write-up with `search/zm_mixed.py`, as `zmcheck` (Rust) and `zeromargin.py` did for `s(32)`.
Cores `taskset -c 0-9` (physical only).  Budget: two working days.

**Independence rules.**  Do **not** read `search/zm_mixed.py`, `search/zm_mixed_test.py`, `search/ZM_MIXED.md`, or any
`ZM_MIXED_AUDIT*`.  You may read: `tasks/line-cover/FORMAT.md` (the statement), `search/LINE_COVER.md` (the cover side),
`search/mixed_cover.py` (reader; or write your own parser, better), `certificates/FORMAT.md`, and the point checkers
(`verify2/src/main.rs` = `zmcheck` + `search/RUNG2.md`, `search/ZEROMARGIN.md`, `search/zeromargin.py`) for pose-space
subdivision and point primitives.

Suggested route (your call): extend `zmcheck` as a **new** binary or mode in `verify2/` (keep existing `zmcheck` behaviour and
its outputs byte-identical; `certificates/s32/` pins it), exact integer/rational arithmetic, with your own sound lower bounds
for segment and polygon mass over a pose box.  The known hard places: tile germs (an axis-parallel square with edges on grid
lines counts those segments in full, a slightly tilted one only partly; the mass is discontinuous there) and tight tilted cells
(float margin `~0.3 %` at `θ ≈ 35–41°`).  Write every lemma with a proof in `search/ZMX2.md`.

Tests: point-only agreement with `zmcheck` on the shipped `s(13)` cover (a column); toy mixed covers that must certify; rejection
tests (holes must fail at the right pose).  Then the full D4 run on the candidate (`--d4` region `[0,5/2]² × u ∈ [0,1/2]`, D4
argument as in `S32_EXACT.md` §7, check D4 invariance of the measure exactly), resumable, with a manifest (checker sha, settings,
roots, result).  Report: `search/ZMX2.md` (lemmas, tests, run table, verdict).  Commit own files, never push.
