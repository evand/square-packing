# k2m4-site: project docs + site for s(k² − 4) = k (2026-10-03)

Repo `~/math/square-packing/public` (main checkout, branch main).  Source of truth: `s12/certificates/k2m4/README.md` (read
in full: claim, per-k trust, earlier work, review record, what is not machine-verified) and `s12/notes/lean-valid-split.md`.
Model everything on how the k² − 3 result (`certificates/k2m3/`, 2026-09-29) appears in each place — match tone, length
and structure; don't overclaim (state trust per k; credit wand125's s(77) = 9 of 2026-10-01; "the deficit-4 case of
Friedman's Conjecture 1", not "F₄", in public prose).

1. `s12/README.md`: a result entry for s(k² − 4) = k, k ≥ 5, next to the k² − 3 one; update the k² − 3 entry's Lean
   sentence (now `bentz_of_validTilt7 : ValidTilt7 → …`; D4 reduction and Lemma Z kernel-checked) and its "known before /
   ours also give" sentence if k² − 4 changes it.
2. `s12/notes/status.md` (frozen at 2026-09-22): add a short dated section listing the exact results since (s(32), s(21),
   s(45), s(60)/s(61), k² − 3, k² − 4) and the Lean state; keep it short.
3. `s12/search/FRIEDMAN.md`: only the F₄ status lines (§0 line ~20, §4, §5, §6 open items, §9 line ~381): F₄ now holds,
   with k ≥ 8 resting on ValidTilt9.  This file has someone else's uncommitted edits in other sections — edit only those
   lines and tell me exactly which lines you changed.
4. Site (`site/`; read `site/README.md`, `site/build.sh`): 
   * `site/data/lower_bounds.json` (+ the identical copy `site/www/data/lower_bounds.json`; it's `cp`'d by build.sh): for
     every n = k² − 4 present in the file, k ≥ 5, add/adjust entries in the style of the k² − 3 ones (source
     "EvanDaniel2026", status "preprint", url the k2m4 bundle, note with per-k basis); where an n already has our own
     bundle entry (21, 32, 45, 60) add only a note mention if the k2m3 pattern did that; check what's already there for
     77 (wand125).  Validate JSON; keep the two files identical.
   * `site/www/sources.html` (the line near the k² − 3 item), `site/www/proofs.html` (an entry next to "Three short of any
     square"), `site/www/index.html` (a card like the k2m3 one), and a page `site/www/k2m4/` modelled on `site/www/k2m3/`
     (shorter is fine; reuse its CSS/structure; no new external scripts).  Update the subnav lists that list k2m3.
   * If the k2m3 page or others state "only the all-k reduction is in Lean", update for ValidTilt7.
5. Run whatever site checks exist (build.sh steps that don't fetch from the network; any link/HTML checker in the repo).

Rules: don't commit, don't push.  Touch only the files above (plus new `site/www/k2m4/`).  Cores 0–3 (taskset).  Report:
files changed, the exact FRIEDMAN.md lines, checks run, anything you were unsure about (≤ 300 words).
