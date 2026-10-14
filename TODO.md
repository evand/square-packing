# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths relative to the repo root (`public/`; flattened from `public/s12/` on 10-07).  Public since 10-03 (outreach stays in `private/s12/outreach`, symlinked as `outreach/`).  Done: `Completed.md`.  Conventions: `~/math/TODO.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`) for current values.**

## Next (Evan, 10-07)
Cleanup session 10-07 done (Completed.md): flatten s12/ → root, exit codes, upstream sync, site refresh.  Next: optimizer work; exact-forms route paused after jlevy#419/#420.  Git: single checkout on `main`; exact-minpoly raw solve outputs (gitignored) in `search/exact/minpoly/solve/`.

Proof-piece wishlist (10-03; picks at the end): `search/WISHLIST.md`.  Picked then: Lean `E(n) ≤ 0` + seam bound `D ≤ m_v`, then `M(1) ≤ ¾` (FRIEDMAN §9.6 item 5); `M(1) = ¾` fully (SEAM_W1 §2.3).

## Exact forms (10-07; `tasks/exact-minpoly/`)
First deliverable: a table of exact forms of S + a handful of Lean local-optimality proofs.  n = 83 (deg 672) optional.
- [ ] Lean local minima: s(11), s(28) done (`Exact/N11L`, `N28L`; sparse G check 10-07); next n = 5 (grade B), one grade C (10 or 19).
- Exact certificates of superseded packings are kept on purpose (Evan, 10-10): the web viewer lists old packings too, and exact forms are a piece it may use later.
- [ ] Exact forms for current records without one (72 as of 10-10: the 55 open + 17 whose record moved on, e.g. SQUISH / Ryan Xu / Mishapolk packings): compute, after the optimizer run.
- [ ] Lean `Packs` for field degree > 20 (41, 51, 69, 87, 106, 128, 152, 177, 205, 266, 300): field arithmetic too heavy (n = 41 > 35 GB); needs a cheaper field representation.
- [ ] Target (Evan, 10-07): exact forms + Lean `Packs n S*` for all n ≤ 82 (open: 29, 55, 68, 71 = msolve-hard; 83 = deg-672 wall), Lean local min for ~70 of them (not 17/41/51/69).  State: 269/324 exact, `Packs` 258, local min 178 (`search/exact/EXACT_FORMS.md`).  Large n Lean: split cert into modules (B294/B310 die on per-file memory).
- [ ] Grade-C exact certificate (`search/exact/gradec.py`, 10-11; `tasks/exact-minpoly/grade-C.md` end): VALID for 233/250 records of field degree <= 4 (230, 261 with `GRADEC_DROPZERO=1`).  Left: 17 timeouts at 30 min (large flat families; rerun uncapped on 8 cores, several hours, Evan's OK), field degree > 4 (incl. 177, 211, 263), Lean grade C for n = 10 or 19.
- [ ] **s(29) exact form found (10-11, `search/exact/homotopy/N29.md`)**: degree 7213 minimal polynomial by homotopy continuation + lattice recovery of 4 missed roots; irreducible, linked to the record by Krawczyk (rigorous), parametrization identities modulo 50 random primes (Monte Carlo).  Data gitignored; backup `runs/n29-backup/` (Evan: Google Drive).  Next (design in N29.md §Pickup; needs Evan's OK for compute): grade-C certificate with a mod-q + arb backend, rigidity of the 19 rigid squares, exact corners in K, EXACT_FORMS row.
- [ ] Coverage, 55 open by cost: flat rotation (10), multi-direction force balance (11), block-triangular elimination (26 + msolve-hard 9), field 108/206.  msolve-hard: 29, 68, 71, 126, 228 are k multivariate equations in k class angles (n = 29: 5 eqs, degrees 12–26; field degree likely > 120), not a factor-choice problem (that fix landed 10-07).  Needs better elimination (block-triangular / triangular decomposition / other solver).  55: force-balance determinant deg 85, 5.7M terms; try multipliers as msolve variables.  Numerics question: task doc §Numerics.
- [ ] Explainer: local optimality of the n ≤ 324 records (for the s(11) explainer to cite). `notes/local-optimality-explainer-draft.md`.

## Display regularization (10-09; `search/regularize/REGULARIZE.md`)
- [ ] **Records drawn as posted, awaiting regularization** (10-11; Evan: hold the compute, track them): 146, 153, 155, 156, 176, 181, 207, 211, 238, 239, 261, 270, 271, 272, 304, 305 (jlevy#500, #501, #503, #504, #505; `raw: true` in `site/data/records2026.json`; poster's certificate, both our checkers, in `lists/certs/raw/`, `search/regularize/posted_cert.py`).  When cores are free: regularize them into `lists/best.json`, then `site/build.sh provenance && site/build.sh lists`.  Before any site refresh, check jlevy/squares issues newer than the last `pk/pending.yaml` entry (10-11: none after #505).
- [ ] **Next (handoff):** site viewer done 10-10 (Explore/Compare/Overview/Bounds draw `lists/best.json`; `site/TODO.md`).  Still to do: the variants page + a short write-up with 3-5 example cases (rotation 268, slide 38, straightening 301, kept 0.8 deg column 175, coincident faces 26 / alternates 267, 302).  Then: offer the lists to jlevy for the poster (draft in `outreach/`, ask Evan).  Open: should a level-3 alternate ever be the best choice; 261 still on its old input (ry-xu's registered packing does not certify: exactsolve −4.5e-14 among redundant contacts; the site draws it as posted); novelty check of our 5 alternates against publications.

## Major: exact results
- [ ] **s(12) strip leaf (10-11, depth-first; paused 10-11)**: class-A 4-strip of width <= 1.90 + seam at x = 1.02 ⇒ no 12-packing; margin ~0.02.  Lean: guillotine + stitch done.  Region lemma: a point cover fails (LP ~7.97), needs branching + count rows.  Resume options in `tasks/strip-leaf/` (Status).
- [ ] k²−5 for all large k: R = w = 5 κ = 0.02 D ≈ 1.30, conditional GO, margin ≈ 0.05 (`K2M4_MARGIN.md`); box 13 pitch 0.1. Make it cheaper first: cert golf (`tasks/cert-slack/`), finer-pitch D at w = 4 (insertable shortcut dead: price ≈ 1.8, `S2_INSERTABLE.md`).
- [ ] Seam capacity `M(∞)` (C5/R2; FRIEDMAN §8–9, `SEAM_W1.md`, `SEAM_1D.md`): `M(w) → ∞` ⇒ families give s(k²−c) = k for every c. Next: dual-side band LP w = 2, 3 (rigorous caps); wall-row + slow-variation lemmas.
- [ ] Corner coupling (R3): does D(w) ≥ κ M(w)? Start with the corner dual at w = 1 (D = ½, corner cost ¼?).

## s(12)
- [ ] Exact elimination along wall-to-wall chains (RANK8 chains of four) with small exact QE pieces: the s(12) home for the wallchain idea (10-10, `search/WALLCHAIN.md` reading).
- [ ] Group-chat lemmas (10-05, low priority, `tasks/n12-chat-lemmas/` lists them): replay their evidence tarball with our checkers; Lean strip demo h = 1.13 and the corner–edge lemma `U+V ≥ 3.91c − 1.91`.

## Correctness / review
- [ ] Nagamochi Lemma 1 with positions free (10-10, `search/NAGQE.md`): hand LP says irreparable for e>0 (margin ~2e); make exact. Whole-problem CAD (z3/cvc5) timed out.
- [ ] Owed on jlevy (drafts: `outreach/draft-jlevy-replies-2026-10-07.md`; #238, #375, #256 posted 10-07): #279 s(59) exact replay **paused 10-07** (Evan: poor RoI on CPU for now; resume when cores are idle): one `zm_mixed` sweep at depth 34 over the whole D4 region (wand125 needed two runs whose joint coverage jlevy couldn't check); 33,725/102,400 roots done, ~127 CPU-h total, the rest is the expensive part.  Resume: `python3 search/zm_mixed.py cert runs/s59_wand125/n59_mixed_cover_8.txt --d4 --cert-mode --disj --chain-from 0 --depth 34 --pitch 1/20 --ubins 16 --nproc 15 --progress 5000 --resume runs/s59_wand125/zm_mixed_d4_depth34.jsonl --manifest runs/s59_wand125/zm_mixed_d4_depth34_manifest.json`; then reply (draft (D)).  Posting needs Evan's OK.
- [ ] Watch jlevy#419: Lean custody + closed forms posted 10-10 evening (all 258 `Packs` at 6bbae37; asks exact `verified_upper_bound` for 41 counts; 15 replays; 232 caveat).
- [ ] (perf, scales) `lean/Sqpack/Bentz.lean` peaks at 23.5 GB / 126 s alone (next: S32 11.8 GB, ValidSplit9 9.8 GB, rest ≤ 7 GB): find the heavy term/tactic before the families grow.
- [ ] Format bridge: read wand125/tokoharu rectangle-density certificates in `zmx2` (converter or reader); first use: replay s(19) ≥ 193/40 (T-103), which alone settles s(18) < s(19) (jlevy may replay it too).
- [ ] `zmx2` exit status: still 0 on NOT VERIFIED / INCOMPLETE (verify, zmcheck changed 10-07); ~30 callers incl. bisect scripts to adapt first.
- [ ] k2m3 2nd impl (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py should-fixes (audit 09-26; + s20 review 10-04: divisibility assert in `roots()` like `d4_roots`, else a non-dividing `--pitch` silently skips a strip): int64 guard on `sum(W·weight)`; `zm_d4_sweep.py summary` sha/settings gate; stale docstrings. Re-pin s(32) bundle.
- [ ] Watch: #465 (new s(132) ≤ 11.9870993322, s(155) ≤ 12.9524989440; posted 10-08; #399 266/270/272 confirmed T-119, closed), #375 (n = 292; n = 17: replied 10-10 that the #419 exact form supersedes S', pinned 8c47162), #419/#420 (exact forms, local minima; no reply yet), chelokot#18 (quiet since 10-02).
- [ ] Issue #3 (squarepacker preprint, c*(k) → ∞, v1.2: κ = 0.0353 log k): marked **claimed (unrefereed)** 10-10 in FRIEDMAN.md, WISHLIST A3, problems.html, Conjectures.lean (Evan: no reimplementation for now; watch what others do).  Partial read: `tasks/issue3-k2-minus-c/review-2026-10-10.md` (no gap found in 4.5/4.9/4.10; Lemma 4.10 boxes and Lemmas 5.3/5.4 unchecked).  Replied on #3 10-10 (left open).  Also check FRIEDMAN.md's Roth–Vaughan constant (10⁻¹⁰⁰ vs the preprint's 10⁻¹⁰).

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` as an open statement in `Spec.lean` (formal-conjectures has only s(11), s(17): offer it upstream too); axis-parallel ≤ 9 (item 1); Lemma CC chain counting (item 5, "Lean-ready").
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares as shared spec; offer bridges to Queuingtheorydotcom (chelokot is moving proofs Lean → Bend, 10-07: re-scope). Draft → Evan approves.
- [ ] **PR #7 (wand125, Valid7 follow-up: ValidTilt7 without Face, LebMass bridge)**: spot-checked 10-12, all OK (`tasks/pr7-valid7-followup/spotcheck-2026-10-12.md`).  Merge on GitHub was blocked for Claude; Evan to merge, then post `outreach/draft-pr7-thanks-2026-10-12.md`.
- [x] **PR #6 (wand125) merged 10-10 (060c1a1)**: `Valid7` / `ValidTilt7` in the kernel, s(k²−3) = k with no hypothesis.  Opt-in `Sqpack.Valid7Checker`; spot check (checker + Bridge built; all 9,800 roots regenerated in a sandbox and scanned, 0 problems; roots not rebuilt, Evan's call): `tasks/pr6-valid7-kernel/`.  Sample root builds done 10-12 with the PR #7 check: 6 roots incl. the heaviest, 939 modules, 12.1 CPU-h, all OK (`tasks/pr7-valid7-followup/`).  Thank-you posted 10-10 (asked for the qualified `example` in Main; watch).
- [ ] `ValidTilt7/9` in Lean (`notes/lean-leb-mass.md` §7): LEB + CAP done 10-03 (`LebMass.lean`). Next: `CovT` tree layer + generator (~400 lines), then PIECE with polygon clip (1–2k lines; with the tree, ~17k of V3's 32k leaves), then Lemma E (3–6k lines, weeks). k2m4 needs a leaf dump.
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Arslanov rectangle decomposition in Lean (squeezable rectangles glue to side < m; s(k²−k) < k for k ≥ 12?).  Read the paper first (2019 E-JC doi:10.37236/8586; 2021 secondhand).  None in Lean anywhere (10-05).
- [ ] Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY): Evan reviews → site page.

## s(17)
- [ ] Wall-chain elimination (10-10, `search/WALLCHAIN.md`): width-1 ring DP misses m851903's obstruction (float). Negative posted as jlevy#494; width-2 band DP unpriced.
- [ ] Idea B (10-10, unstarted): stationarity / monotonicity-test pruning as the "non-margin argument" for n17 capture (X-052 §2.8 soft directions); theory + cheap pilot.
- [ ] s(17): pure closed covers are at their limit (exact ceiling ν_f(4.660) ≥ 17.0447, `CEILINGS_17_20.md`); further lower-bound gains need rule atoms (WISHLIST N15).
- [ ] Techniques from jlevy/Kleddamag/Guzhou (`notes/jlevy-s17-techniques.md` §2): rule atoms as LP columns, near-tight cells as exact LP rows gmpy2 in exact checkers.  s(17) exact: their PR #307 is far ahead and wand125 built its independent check (jlevy PR #410): not ours to do.
- [ ] s(17) exact via the s(11) method: on hold (CPU). `tasks/s17-core-isolation/`

## Packings (10-14 triage with Evan: `tasks/packer-arch-2026-10-14/`; handoff `search/packer/HANDOFF-2026-10-14.md`)
Thesis (Evan 10-14): diversity pays (new mechanisms, new sources); local search at current budgets has mined the frontier.  No compute above n = 324 without a specific reason.
- [ ] Agents running 10-14: A1 bandit timescale (`search/packer/analysis/bandit_timescale.md`), A2 similarity / clusters at 110 (`analysis/similarity110.md`), A3 ingestion watcher (`pk/watch.py`).
- [ ] Parent sampler `pk/parents.py`: pool query (n, n±k, cluster, gap band, as-of), named policies; default uniform over clusters then within; reads the `expansion` ledger.
- [ ] Run config with test / record modes (sources, mechanisms, reference, reward, grid-rule knob, control); chain imports its own rounds; cfp into the store on import.
- [ ] **TT2** time-travel test, strongest settings vs pre-10-14 chain, ~2 × 30 CPU-h; plan in `tasks/packer-arch-2026-10-14/`; pre-register in PACKER.md.
- [ ] Parent-distance model (after A2): fit P(gain | parent gap, distance, n band, mechanism) from logs, then a randomized gap-band slice to de-confound.
- [ ] Near-miss queue (store tag) + a finisher that drains it (harden, long polish).
- [ ] Constructors (Evan 10-14): (1) off-diagonal piece library (free W × H, fq `--rect-dx/--rect-dy`) + offset-converter screen; (2) distorted-strip / V-band (fix `vband.py` fit: 41 / 80 tilted kept at 241); (3) two-bulb, lopsided bulb.
- [ ] Piece sources: crop anywhere (edge strips, interior), pieces from n > 324 and catalogue revisions, published rectangle-packing records.
- [ ] Refactors (bit-exact, check by replay): one fq call path through `QOpt`; shared neighbour transforms (remove / insert k, crop) in `pk/moves.py`.
- [ ] Perf follow-ups (each ≤ ~10 %, PACKER "Performance session 10-13"): block LDL, SIMD `Prob::eval`, `FQ_SOC_ROUNDS=16`, exactsolve `probe_flat` factor reuse; measure with `fqreplay.py` / `chainreplay.py`.
- [ ] Run `fresh.py` after every sync (A3 automates); fix its finishing gap (364).
- [ ] exactsolve should snap near-axis corner-loaded squares itself (n = 292, `search/packer/s292.md`).
- [ ] Shallow-valley descender X01 (Evan: eventually): follow the corner-corner MILP direction with exact contacts until jammed; test on the 225 stuck n = 110 points (PACKER "X01").
- [ ] Store backup (last 10-11, `~/store-backups/`; store now ~728k rows after the 10-14 run import); automate if repeated.
- [ ] s(110) explainer refresh, iterating (`tasks/landscape110-explainer/`); s(110) write-up = rewrite of `search/packer/s110-landscape.md`'s census section.
- [ ] Time-lapse records video (`site/timelapse/`; Evan to return): store sighting times are local, not UTC.
- [ ] Watch open jlevy issues: #596, #584, #552, #489.
- [ ] Evan to decide: n = 88 strip-family test set (~16 CPU-h): drop unless TT2 results look 110-shaped.

## Outreach / hygiene
- [ ] Site: draw every packing we have numbers for (Couzo's 49, itsnaka's 23 SQUISH (local copies + pins: `tasks/site-all-packings/README.md`), de Winter's 211, ours 266/270/272, n > 324; replaces the `NEWER` stopgap in viewer.js); coords with credit, no code/text. `tasks/site-all-packings/`
- [ ] Announce: VibeMathed entries for s(21), s(32) (and now k²−3, k²−4 / F₄); X drafts `outreach/drafts-s21-s32.md`.
- [ ] Repo unification, rest: outreach → a small private GitHub repo (versioned, backed up); `wip/*` branches + pre-push hook if needed.
