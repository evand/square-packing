# Commits since 167d842 for review (s(21) = 5), 2026-09-27 -- NOT pushed

* `b37ecf2` mixed_cover.py: reader/validator/float evaluator for mixed (point+segment+polygon) covers, format v1
* `443eac0` zm_mixed.py: exact checker for mixed (point + segment + polygon) covers
* `78b75a1` line-cover A: mixed (point + line-density) covers; m=5 candidate 20.832 at 1.0025x confirmed min (heuristic GO)
* `d4b97f9` zm_mixed.py: region-wise piece/point coupling (SPLIT, Lemma R), hull-edge Lemma L ends, Lemma V, --resume
* `e72f317` ZM_MIXED.md: round 2 results; m = 5 candidate x 1.003 (total 20.8947) VERIFIED-D4 by zm_mixed.py
* `0462eff` Lean: s(21) = 5 from the zm_mixed D4 checker statement; packing bound and D4 reduction for measures, mixed covers
* `5d91695` zmx2: independent second exact checker for mixed covers; s(21) candidate VERIFIED-D4 and VERIFIED (unreduced)
* `b9d1e07` zmx2_manifest.txt: official runs at 5d91695 (VERIFIED-D4 2,500 roots; VERIFIED unreduced 20,000 roots; 0 uncertified)
* `14eb9bd` ZM_MIXED_AUDIT.md: adversarial audit of zm_mixed.py and the s(21) run (no soundness defect; provenance should-fixes); zm_mixed_audit.py component/rejection tests
* `457f876` ZMX2_AUDIT.md: adversarial audit of zmx2 (no claim-affecting defect; must-fix: unchecked i128 input overflow gives VERIFIED for false statements); zmx2_audit/ tests
* `dd6f63a` zmx2: bound every input integer in the parser (ZMX2_AUDIT.md F1: i128 wrap); tests: overflow inputs, differential cert (A7)
* `e49dc67` ZMX2.md sec 11: the post-audit parser fix, tests and reruns; s(32) README: zmx2 as a third checker
* `fe2ac0b` zm_mixed.py: provenance and certificate mode (ZM_MIXED_AUDIT.md S2, S3, S5, N1, N4); no lemma changed
* `c7e2606` Lean s(21): point at the bundled cover certificates/s21/s21_mixed_cover_5.txt (same bytes, same sha256); comments only
* `086a129` certificates/s21: s(21) = 5 bundle (mixed cover, zm_mixed.py --d4 --cert-mode run, zmx2 --d4 and --full runs, verify.sh)
* `5ae38b9` s(21) = 5: wire in (verify.sh + CI with the zmx2 full run, Pages /s21/), write-up docs/s21.html, site bounds/sources/nav, READMEs

Bundle-step notes: zm_mixed cert-mode run 40,000/40,000 roots, 461,204 boxes, 0 uncertified (census identical root
for root to e72f317's run); zmx2 (after the F1 parser fix dd6f63a) --d4 2,500 and --full 20,000 roots, 0 uncertified,
census identical to the pre-fix runs; s(32) --pair-points re-verified with the patched binary.  ./verify.sh fast tier
green (incl. certificates/s21/verify.sh), site/build.sh clean (only lower_bounds.json changed), lake build clean.
To check before pushing: wand125's 399/80 date (no date recorded in lower_bounds.json / docs); site/www/sources.html
now lists jlevy/squares; search/zmx2_manifest.txt (sibling's, runs/ paths, pre-fix binary) is superseded by
certificates/s21/zmx2_*/manifest.txt.
