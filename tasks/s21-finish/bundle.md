# Brief: bundle and ship-ready s(21) = 5 (do NOT push), 2026-09-27

Read `tasks/s21-finish/README.md`, then: `certificates/s32/` (README, verify.sh, layout: the model to follow),
`search/ZM_MIXED_AUDIT.md` (should-fixes 1–4 and nits), `search/ZMX2.md` + `search/zmx2_manifest.txt`, `notes/lean-s21.md`,
`search/LINE_COVER.md`, `s12/verify.sh`, `.github/workflows/{verify,pages}.yml`, `docs/s32.html`, `site/` (how s(32) is wired
into the site: `site/build.sh`, `site/data`, `site/www`), top-level `s12/README.md`.  Cores `taskset -c 0-9` for runs.
Budget: half a day.  A sibling agent is auditing `zmx2` (read-only); if it reports a must-fix, the lead will tell you.

1. **Audit should-fixes in `zm_mixed.py`** (no change to any lemma's logic): record shas of `zm_mixed.py`, `mixed_cover.py`,
   `zeromargin.py` and the input file plus every setting (SPLIT, θ-bias, chain-from, depth, caps) in the run log/manifest;
   `--resume` refuses on any mismatch; a `--cert-mode` (or similar) that disables Corollary T′ and polygon code paths, and
   refuses covers that contain polygons.  Update `ZM_MIXED.md` (item 3: cite component + rejection tests, not leaf stress).
   Rerun `zm_mixed_test.py selftest`, the audit's `zm_mixed_audit.py` quick tests, and the rejection tests.
2. **Fresh certified runs with the final code**: `zm_mixed.py --d4 --cert-mode` on the cover (≈ 2 h on 10 cores) and `zmx2`
   D4 + no-symmetry, each writing a manifest with shas + settings + result, as `certificates/s32/zeromargin_d4/` does.
3. **`certificates/s21/`**: the cover (`s21_mixed_cover_5.txt`), `FORMAT.md` for the mixed format (from `tasks/line-cover/FORMAT.md`,
   public wording), README in the s(32) style (claim, one-paragraph proof, "what is checked by what" table incl. Lean
   `s21_eq_five_of_checker`, the two checkers, the float-interval caveat of zmx2 stated plainly, audits), run manifests/summaries,
   `SHA256SUMS`, `verify.sh` (default tier: sha check, exact total, D4 invariance, Lean data `--check`, **zmx2 full no-symmetry run**
   (~90 s); `--full`: also the zm_mixed D4 run).
4. **Wire in**: `s12/verify.sh`, `.github/workflows/verify.yml` (zmx2 run in CI; zm_mixed full run local only, like S32_SEPARATE),
   `pages.yml` (serve `/s21/`), write-up `docs/s21.html` (style and depth of `docs/s32.html`: the line-density idea is the story;
   credit prior bounds: our 5000/1001, wand125 249/50 → 399/80, jlevy 122/25, Nagamochi, Friedman 4.7438), site data/records, README.
   Also mention s(32) is now re-verified by a third checker (zmx2) if you add it to s32's README table.
5. Run `./verify.sh` (default tier) and the site build locally; everything green.  Commit (own files), **never push**; leave a
   summary of all commits since 167d842 for Evan's review.  Final reply `<= 250` words: what's in the bundle, run results, commits.
