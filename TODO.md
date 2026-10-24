# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths from the repo root.  Done: `Completed.md`.  Conventions: `~/math/TODO.md` (one line per entry, details in `tasks/`).
Pre-tidy text (10-15, has details dropped here): `tasks/todo-tidy-2026-10-15/TODO-before.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`).**

## Next session (Evan 10-23)
- [ ] Public / private split (10-23): new search methods (moves, generators, bandits, priors, Kaggle tooling) and their notes / TODO live in the private plugin repo `../plugins` (`search/packer/DESIGN.md` "Plugins"); this repo keeps the framework, store, verification, exact work and posted records.  `.git/hooks/pre-push` refuses private paths and unreviewed Kaggle mentions.  Posting default (Evan 10-23): post what others already have (e.g. sides that match public ones); revisit later.
- [ ] Then posts beyond the register (`search/packer/regcheck.py`, policy `search/packer/SUBMITTING.md` "Beyond the register"): **268, 303 FIRST** (side-only on the Kaggle per-puzzle list 10-22, we hold the packings); **268, 303, 304 PUBLIC** (Kaggle notebook packings below the register: credit their authors).  Every one: exact check, same-day live check incl. a fresh side list, Evan OK per post.

## Waiting on Evan
- [ ] Review the Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY) → site page.
- [ ] Compute OKs: s(29) grade-C (`search/exact/homotopy/N29.md` §Pickup); 17 grade-C timeouts uncapped; regularize the 16 raw records; #279 replay (`tasks/s59-replay/`).
- [ ] Decide: n = 88 strip-family test set (~16 CPU-h; TT2 is done, results not 110-shaped → drop?).
- [ ] Posts: VibeMathed s(21), s(32), k²−3, k²−4 (`outreach/drafts-s21-s32.md`); variant lists to jlevy (after the write-up below).
- [ ] Time-lapse records video (`site/timelapse/`): Evan to return; store sighting times are local, not UTC.

## Watch
- [ ] jlevy #702 (our s(343) ≤ 18.977682482405, beyond the corpus, posted 10-20, cert 88385ab), #703 (our s(270) ≤ 16.928162252197, posted 10-20, cert 88385ab), replies 10-20 on #664 (271: a different minimum, finder us), #669 ("by" bug), #674 (229 questions); #680 (our s(183) ≤ 13.999580178737, first 183 below the grid, posted 10-18, cert 3a9948c; lineage Couzo 182 -> insertion -> explorer), #674 (our s(229) ≤ 15.650163254859, posted 10-17, dilation cert e0338c7; updated same day to ≤ 15.650161513896, cert 00e68df; start: #672 / #668), #673 (credit correction: found_by at 146 and 308 to us, per the same-minimum rule; posted 10-17), #669 (our s(155) ≤ 12.952260944651, posted 10-17, cert 59326b8), #664 (our s(271) ≤ 16.948947733714, posted 10-16, cert 757b40b), #659 (our s(236) ≤ 15.861934390608, posted 10-16, cert 5ad0069), #639 (exact forms, 22 current records), #606 (our 146: T-187 V0/C0 10-15, awaiting their replay review), #545 (s29 / grade-C retained 10-15; their replays wait on budget; they note exact-point non-overlap for s(29) open), #596, #584, #552, #489, #465, #419 / #420 (exact forms, no reply), #375, #3 (preprint, left open); chelokot#18 (repo unknown: not under bend / bend-math).

## Packings (`search/packer/`: `HANDOFF-2026-10-14.md`, `tasks/packer-arch-2026-10-14/`)
Thesis (Evan 10-14): diversity pays; local search at current budgets has mined the frontier.  No compute above n = 324 without a specific reason.
- [ ] Sink after every job (10-23: 4,291 explorer packings from two runs sat outside the store until a manual `pk.py import-run`): `jobs.py start` runs `pk.py sink` on exit (one sink at a time).
- [ ] familymap.classify: D/S needs 80% of tilted within 1 deg of 45, so slightly distorted diamonds (our 372 best: 44-46 deg) come out "2"; widen to ~3 deg or add a distorted-diamond class, and store the code as a packing column/tag so "irregular" is a SELECT.
- [ ] Kaggle per-puzzle side list is a 10-22 snapshot (`runs/kaggle/lb_perpuzzle_2026-10-22.json`, staff notebook output): refresh it (re-download the notebook output, or its successor) before any live check; `candidate.py --live` should run `regcheck.py` on its n (it checks register / issues / kingbird / store, not side-only lists), and `pk.py ingest` could print regcheck's table.
- [ ] Live check misses franciscouzo's claims (10-23): `pk.py watch` can't fetch his repo ("files over 1073741824 bytes"), so #763 (436 bounds, exact certificates in `upper_bounds.json` / `certificates/` at 964ee2b) never reached the store; we nearly posted 268 / 303 / 304 duplicates.  Fetch his `upper_bounds.json` (15 MB) + the certificates it names instead of the whole repo, and have `candidate.py --live` read it.
- [ ] Kaggle CSV imports have no finder (sighting finder NULL; the notebook author is only in the source path): record the author as finder in `kaggle_import.py` so PUBLIC posts credit the right person.
- [ ] Fill module (10-20, PACKER "Fill module"): `pk/fill.py` + `pk/schedules.py` + `fillbench.py` built; next a pre-registered refill bench (placer x relax x stop x K, incl. K = 1 / softmax), then defaults.
- [ ] Stages / chains / trials refactor (Evan 10-22: only new work moved so far; `search/packer/DESIGN.md` "Not done yet"): loop combinators (`repeat`, `accept`), **first an accept-loop tail for `trial.py`** (its 10-22 tail is one perturb, which only worsens good starts), then port the hand-wired scripts (explore, chain, hop, compose2, relax, fresh) one at a time with a trial showing each port matches, `trial.py` sequential statistics + a public suite (110).
- [ ] Archive pass (Evan 10-20, own session): move unused-but-interesting code (old generators on `gen.mis_fill*`, legacy fills, one-off scripts) to an archive dir.
- [ ] s(110) corner-cell holdouts (census 1,078 minima, 19:12 10-15) + a stop-early descent mode (GX1: crossings collapse to the record basin).  §4.
- [ ] Detour, parked (Evan 10-15): our own register to n = 1024-4096 (every published regular packing, L groups, mechanical constructions), then hunt / submit records > 324.
- [ ] Dedupe (`pk/ident.py`, PACKER "Conservative dedupe"): applied once 10-20 at 110 / 272-274 (-91k rows); rule must be re-discussed with Evan before any further `--apply`.  Next ideas: contact-graph refusal from exactsolve reports; run-side cluster weighting for seed selection.
- [ ] (#607 added 10-15 by `pk.py watch --apply`.)  scoring rule for our own posted results inside a scored run (146 in TT2).
- [ ] Run config with test / record modes; chain imports its own rounds; cfp into the store on import.
- [ ] Perf follow-ups (each ≤ ~10 %; PACKER "Performance session 10-13"): block LDL, SIMD `Prob::eval`, `FQ_SOC_ROUNDS=16`, exactsolve `probe_flat` reuse.
- [ ] `fresh.py`: fix its finishing gap (364).  (It runs after every sync via `pk.py ingest`.)
- [ ] exactsolve should snap near-axis corner-loaded squares itself (n = 292, `search/packer/s292.md`).
- [ ] s(110) weak residue: 12 side keys / 10 minima "negative curvature blocked" (cone-QP check, then branch-aware second order); exactsolve fallback fixed 10-15 (PACKER "exactsolve fix").
- [ ] Shallow-valley descender X01: `search/packer/nlpdesc.py` (second-order fixed-line NLP, IPOPT) resolved 34 / 40 unresolved s(110) keys 10-15 (PACKER "unresolved / weak: second-order descent"); next: the 225 stuck points, the 6 left (MILP-branch kick).
- [ ] s(110) explainer refresh + write-up (`tasks/landscape110-explainer/`; rewrite of `s110-landscape.md`'s census section).

- [ ] Next jlevy submission: bundle s(260)'s exact form (S* 16.645778786061741…, certified 10-17, `search/packer/runs/candidates/s260_cs3/`, store kind exact-form): the analytic optimum of de Winter's #668 packing, 1.6e-8 below the posted floats; credit de Winter.  Not its own issue (Evan 10-17).
- [ ] Packer bugs (10-22): `dgen.py` sweep frontier() ignores our pending claims (SW3 reported a fake 343 "record" vs the register, 4.2e-4 above our #702); `explore.py` final refine admits grid-locked results (17.000000 at 274) into the archive.


## Major: exact results
- [ ] s(12) strip leaf (paused 10-11): 4-strip width ≤ 1.90 + seam at x = 1.02 ⇒ no 12-packing, margin ~0.02; region lemma needs branching. `tasks/strip-leaf/`
- [ ] k²−5 for all large k: conditional GO, margin ≈ 0.05 (`K2M4_MARGIN.md`); make it cheaper first (cert golf `tasks/cert-slack/`, finer-pitch D at w = 4).
- [ ] Seam capacity M(∞) (FRIEDMAN §8–9, `SEAM_W1.md`, `SEAM_1D.md`): dual-side band LP w = 2, 3; wall-row + slow-variation lemmas.
- [ ] Corner coupling (R3): does D(w) ≥ κ M(w)?  Start with the corner dual at w = 1 (D = ½, corner cost ¼?).
- [ ] k² − c package (Fang, Ryu, new bounds), Evan 10-17/18: **Lean done 10-18** (`certificates/k2-minus-c/`, write-up `site/www/k2c/`): Fang's theorems; Ryu's 8⌈√(k−4)⌉ − 1 and Theorem 1.2 (k ≥ 8·10⁶); new c*(k) ≤ 4⌈√(2k−7)⌉ − 1 (k ≥ 7) and < 18.9 k^{2/5} + 35.7 k^{1/5} + 122 (k ≥ 3·10⁵).  Packings page `site/www/k2c/packings.html` (triangle of k and c, every construction's packing drawn, viewer; checked by `search/ryu/explorer_check.mjs`) built 10-18, Evan reviewing before Pages.  Open: push + posts (drafts in `outreach/`, each needs Evan's OK); better linear pieces (`tasks/fang-gluing/` step 3); beyond-324 records from the filled L-shape (k ≳ 45).  Junction flip (10-19, `search/ryu/junction.py`): c lower by 1-2 at 374 / 386 k ≤ 400, exact checks; Lean not started (`tasks/fang-gluing/` "Merged junction corner").  `tasks/fang-gluing/`
- [ ] Exponent of c*(k) from above (10-19): exponent calculator `search/ryu/exponents.py` reproduces 1/2, 2/5, 3/8 and shows deeper nesting stalls at 3/8; the blocker is the tilt-change (block transition) cost in Ryu's wall filler; free transitions → 1/3.  Next: a cheaper tilt-change primitive or a local lower bound.  `tasks/k2c-exponent/`
- [ ] Lower bounds on c*(k) with sane thresholds (Evan 10-18): limits of covers (is L(k) → ∞?), Ryu's constants + hybrid, new cover geometries; c = 5 open for 7 ≤ k < 2.5·10⁹.  `tasks/k2c-lower/`
- [ ] Proof-piece wishlist `search/WISHLIST.md` (10-03 picks: Lean `E(n) ≤ 0` + seam `D ≤ m_v`, then `M(1) ≤ ¾`, then `M(1) = ¾`).

- [ ] Certification pipeline (Evan 10-17): benchmark the stuck points with every existing tool, then refactor into one `certify.py`; first step done-now: candidate.py falls back to dilate_cert.  `tasks/certify-pipeline/`

## Exact forms (`tasks/exact-minpoly/`, `search/exact/EXACT_FORMS.md`)
- [ ] Target (10-07): exact + Lean `Packs` for all n ≤ 82 (open 29, 55, 68, 71; 83 deg 672), Lean local min ~70.  10-11: 269/324 exact, `Packs` 258, local min 178.
- [ ] Lean local minima: next n = 5 (grade B), one grade C (10 or 19).
- [ ] Exact forms for current records: 272 / 324 (10-15, `search/exact/cur/`: +23); 52 open by failure in its README (22 OOM, 18 timeout → homotopy with a path cap, 7 multi-direction force balance, 3 flat rotation).  Paused 10-15.  Posted as jlevy#639 (10-15).  Superseded certificates stay (Evan 10-10).
- [ ] Coverage of the 55 open, by cost: flat rotation, force balance, block-triangular elimination, msolve-hard 29 / 68 / 71 / 126 / 228 (task doc §Numerics).
- [ ] Grade-C certificates (`search/exact/gradec.py`, `grade-C.md`): 233/250 VALID at degree ≤ 4; left 17 timeouts, degree > 4.
- [ ] s(29) exact form (degree 7213, `search/exact/homotopy/N29.md`): grade-C cert (mod-q + arb), rigidity, EXACT_FORMS row; data `runs/n29-backup/`.
- [ ] Lean `Packs` for field degree > 20 (41, 51, 69, 87, ...): needs a cheaper field representation; big certs: split into modules.
- [ ] Explainer: local optimality of the n ≤ 324 records (`notes/local-optimality-explainer-draft.md`).

## Display regularization (`search/regularize/REGULARIZE.md`, `site/TODO.md`)
- [ ] 16 records drawn as posted (146 153 155 156 176 181 207 211 238 239 261 270 271 272 304 305, `raw: true`): regularize, then `site/build.sh provenance && lists`.
- [ ] Variants page + write-up with 3-5 cases (rotation 268, slide 38, straightening 301, 175's kept column, 26 / 267 / 302).
- [ ] Open: should a level-3 alternate ever be best; 261 still on its old input; novelty check of our 5 alternates.

## s(12) / s(17)
- [ ] s(12): exact elimination along wall-to-wall chains (RANK8 chains of four) with small exact QE pieces (`search/WALLCHAIN.md`).
- [ ] s(12) group-chat lemmas (low priority): `tasks/n12-chat-lemmas/`.
- [ ] s(17) wall-chain: width-1 ring DP misses m851903's obstruction (jlevy#494); width-2 band DP unpriced (`search/WALLCHAIN.md`).
- [ ] s(17) Idea B (unstarted): stationarity / monotonicity-test pruning for n17 capture (X-052 §2.8); theory + cheap pilot.
- [ ] s(17) lower bound: closed covers at their limit (`CEILINGS_17_20.md`); next needs rule atoms (WISHLIST N15, `notes/jlevy-s17-techniques.md` §2).
- [ ] s(17) exact via the s(11) method: on hold (CPU). `tasks/s17-core-isolation/`

## Correctness / review
- [ ] Nagamochi Lemma 1 with positions free (`search/NAGQE.md`): hand LP says irreparable for e > 0 (margin ~2e); make exact.
- [ ] `lean/Sqpack/Bentz.lean` peaks at 23.5 GB / 126 s: find the heavy term / tactic before the families grow.
- [ ] Format bridge: read wand125 / tokoharu rectangle-density certs in `zmx2`; first use s(19) ≥ 193/40 (T-103) ⇒ s(18) < s(19).
- [ ] `zmx2` exits 0 on NOT VERIFIED / INCOMPLETE; ~30 callers (incl. bisect scripts) to adapt first.
- [ ] k2m3 2nd implementation (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py rest (pitch assert + int64 guard done 10-15): `zm_d4_sweep.py summary` sha / settings gate; stale docstrings; re-pin s(32) bundle.
- [ ] Issue #3 preprint: Lemma 4.10 boxes, Lemmas 5.3 / 5.4 unchecked (`tasks/issue3-k2-minus-c/`); FRIEDMAN.md Roth–Vaughan constant (10⁻¹⁰⁰ vs 10⁻¹⁰).

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` open statement in `Spec.lean` (offer to formal-conjectures); axis-parallel ≤ 9; Lemma CC.
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares; bridges to Queuingtheorydotcom (chelokot moving to Bend: re-scope).  Draft → Evan.
- [ ] `ValidTilt7/9` (`notes/lean-leb-mass.md` §7): `CovT` tree layer + generator, then PIECE with polygon clip, then Lemma E (weeks).
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Unit squares in rectangles, one n at a time (Evan 10-20; `search/rectfront/RECTFRONT.md`): **All n, strips (10-22)**: `strip_low` (H < 2√2 − 1 ⇒ L ≥ n, every n); **`strip_odd_iff`: odd n exactly for 1.87 ≤ H < 2** (`StripC.lean`, slope walk + parity potentials).  Open: odd n on 2√2 − 1 ≤ H < 1.87 (3-step prefixes, length condition), even n (Morandi family); wand125 even n ≥ 12 entries beaten by the extended Morandi family (tell them? ask Evan).  R₂, R₃ done in Lean (`rect2_iff`, `rect3_iff`, chain proof, 10-21; boundary credit: Friedman's Math Magic, June 2010).  Golfed 10-21 (reusable steps in `RectChain.lean`).  R₄ done 10-21 (`rect4_iff`).  R₅ done 10-21 (`rect5_iff`, `Rect5Thm.lean`; wide part by point covers incl. 10-point 2/3–1/3 cover C, `search/rectfront/r5pts/R5PTS.md`; low strips by midline points; high strips `final5_high` by generated certificates, `search/rectfront/r5strip/FINAL5.md`).  Open: a better-style strip proof (Evan 10-21: "prove it first, better style later"; ideas in `search/rectfront/LESSONS.md`); web view of cover C (like the s(12) page) suggested by Evan 10-21.  Before any public mention: credit Friedman's page and its solvers; wand125's `problems/rectangles/` has the same edges as certified constructions.
- [ ] Arslanov rectangle decomposition in Lean (s(k²−k) < k for k ≥ 12?); read the paper first (E-JC doi:10.37236/8586).

## Outreach / site
- [ ] Site: show exact sides (minimal polynomial + isolating interval) for every n that has one (`search/exact/cur/RESULTS.md`, `search/exact/minpoly/`); none shown today (Evan 10-17: later, math first).
- [ ] Per-n history viewer: what changed between consecutive records (matched-square diff, like the s(146) vs #505 page); design first.
- [ ] Draw every packing we have numbers for (`tasks/site-all-packings/`; replaces the `NEWER` stopgap in viewer.js).
- [ ] Site build speed (when cores are idle): `site/tools/import_lists.py` -> `analysis.analyze()` takes 160-1,510 s per changed drawing on one core (10-19: step-229-1 1,508 s, step-229-0 759 s); profile it (cProfile on step-229-1: LPs vs `free_region` vs mpmath).  Also cap BLAS threads in import_lists itself (uncapped, 8 workers drove load to 108 on 10-19).
- [ ] Repo unification: outreach → a small private GitHub repo (the plugin repo `../plugins` is separate: teammates see it, outreach stays out of it).
