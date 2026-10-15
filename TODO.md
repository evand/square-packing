# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths from the repo root.  Done: `Completed.md`.  Conventions: `~/math/TODO.md` (one line per entry, details in `tasks/`).
Pre-tidy text (10-15, has details dropped here): `tasks/todo-tidy-2026-10-15/TODO-before.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`).**

## Waiting on Evan
- [ ] PR #7 (wand125, ValidTilt7 without Face): spot check OK 10-12 (`tasks/pr7-valid7-followup/`); merge (blocked for Claude), then post `outreach/draft-pr7-thanks-2026-10-12.md`.
- [ ] Review the Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY) → site page.
- [ ] Compute OKs: s(29) grade-C (`search/exact/homotopy/N29.md` §Pickup); 17 grade-C timeouts uncapped; regularize the 16 raw records; #279 replay (`tasks/s59-replay/`).
- [ ] Decide: n = 88 strip-family test set (~16 CPU-h; TT2 is done, results not 110-shaped → drop?).
- [ ] Posts: VibeMathed s(21), s(32), k²−3, k²−4 (`outreach/drafts-s21-s32.md`); variant lists to jlevy (after the write-up below).
- [ ] Time-lapse records video (`site/timelapse/`): Evan to return; store sighting times are local, not UTC.

## Watch
- [ ] jlevy #606 (our 146), #596, #584, #552, #489, #465, #419 / #420 (exact forms, no reply), #375, #3 (preprint, left open), PR #6 `example` ask; chelokot#18.

## Packings (`search/packer/`: `HANDOFF-2026-10-14.md`, `tasks/packer-arch-2026-10-14/`)
Thesis (Evan 10-14): diversity pays; local search at current budgets has mined the frontier.  No compute above n = 324 without a specific reason.
- [ ] TT3: gating PASSED (C42, 2.08x); melt / GLS v0 neutral (C43); parent value model PASSED on the 12-round subset (C44, ~1.7-2x over G) and on the full 49 rounds (PF49 10-15: 2.81x [2.22, 3.45]); no transfer at 110 (absolute gap bins).  Next: a relative-gap model (110), then the all-n hunt-style test (new minima + regions + records).
- [ ] Baseline on the major tools (Evan 10-15): per family (A composition, B soft-to-rigid, F/G removal, J constructors) level 1 = reproduce events from named parents, level 2 = blind time travel; improvements only after.  `tasks/composition-baseline/`, PACKER "Baseline on the major tools".
  - Status 10-15 04:35: rectangle blocks (Ostashev #528/#562) 5 / 6 blind (C45; 303 needs a (62,6) block at his 5.90996, ours 5.91536 after 10 CPU-h); Mondrian v1 (#481) L1 4 / 5 informative via slack composition + explorer (C46; harden not needed); #500 (v2 narrowed) 0: pieces / restructure missing; removal 2 / 12; graft 3 / 6; J (Xu #432) 0 / 10 by nearby search (needs a generator).  Next: L2-481 blind, then the live compose(slack) + explore hunt over n <= 324.
  - J generator: EV1 (automated `claude -p` loop) 0 / 18, superseded; EV2 (Claude Code subagents, prompts as the unit of evolution, `search/packer/evolve2/`): gen 0 reproduced 178 / 126 / 145 with the record shown; gen 1 running on the holdout row 234-240 (near-record refs only, `--fill`, softmin).  Assess after gen 1 (Evan 10-15).
- [ ] **Fill performance (Evan 10-15: "important performance work"; next after EV2 gen 1).**  `--fill` = `mondrian.exact_fill` per candidate (MIS MILP over ~1,700 candidates, 0.4-0.6 s, forked child, retry ladder x4) dominates EV2 eval time (median 15.6 min / call).  Plan: (1) safe reductions: corner rule (a candidate whose conflicting neighbours all lie weakly right and above is simplicial -> take it; 4 orientations) + domination + components; kernel = lattice seams + the band at the tilted squares; (2) compare on a corpus (gen-1 eval candidates, L2e composition fills): full MILP vs reduce + kernel MILP (must agree) vs reduce + greedy kernel vs tolerant fill (overlap tau); metrics: count error, time, and downstream quench outcome / side.  Use-case tiers: exact where the count is the claim; very good + fast for screens / evaluator; tolerant only where the quench follows and only at a tau the downstream test supports (P006: tau 0.025 overstated counts by 2).  Cheap wins: one retry only when short by <= 3, 1 seed / 12 candidates in fill mode, dedupe tilted sets, in-process MILP.  Latent bug: `exact_fill` forked after HiGHS has run in the same process can deadlock (60 s wait, returns unfilled); guard it.
- [ ] Detour, parked (Evan 10-15): our own register to n = 1024-4096 (every published regular packing, L groups, mechanical constructions), then hunt / submit records > 324.
- [ ] Store tiering before cheaper proposals raise the row rate: coords near the frontier + cluster representatives, lineage-only (replayable) rows elsewhere.  Design first.
- [ ] Near-miss drain `nmdrain1014b`: stopped 10-15 (0 gains in 42 items, 100 queued); low value, resume only if cores idle.
- [ ] (#607 added 10-15 by `pk.py watch --apply`.)  scoring rule for our own posted results inside a scored run (146 in TT2).
- [ ] Run config with test / record modes; chain imports its own rounds; cfp into the store on import.
- [ ] Parent-distance model: fit P(gain | parent gap, distance, n band, mechanism) from logs, then a randomized gap-band slice.
- [ ] Constructors (Evan 10-14): off-diagonal pieces (fq `--rect-dx/--rect-dy`) + offset converter; distorted strip / V-band (fix `vband.py` fit); two-bulb, lopsided bulb.
- [ ] Piece sources: crop anywhere (edge strips, interior), pieces from n > 324 and catalogue revisions, published rectangle-packing records.
- [ ] Perf follow-ups (each ≤ ~10 %; PACKER "Performance session 10-13"): block LDL, SIMD `Prob::eval`, `FQ_SOC_ROUNDS=16`, exactsolve `probe_flat` reuse.
- [ ] `fresh.py`: fix its finishing gap (364).  (It runs after every sync via `pk.py ingest`.)
- [ ] exactsolve should snap near-axis corner-loaded squares itself (n = 292, `search/packer/s292.md`).
- [ ] Shallow-valley descender X01 (eventually): corner-corner MILP direction with exact contacts until jammed; test on the 225 stuck n = 110 points.
- [ ] s(110) explainer refresh + write-up (`tasks/landscape110-explainer/`; rewrite of `s110-landscape.md`'s census section).

## Major: exact results
- [ ] s(12) strip leaf (paused 10-11): 4-strip width ≤ 1.90 + seam at x = 1.02 ⇒ no 12-packing, margin ~0.02; region lemma needs branching. `tasks/strip-leaf/`
- [ ] k²−5 for all large k: conditional GO, margin ≈ 0.05 (`K2M4_MARGIN.md`); make it cheaper first (cert golf `tasks/cert-slack/`, finer-pitch D at w = 4).
- [ ] Seam capacity M(∞) (FRIEDMAN §8–9, `SEAM_W1.md`, `SEAM_1D.md`): dual-side band LP w = 2, 3; wall-row + slow-variation lemmas.
- [ ] Corner coupling (R3): does D(w) ≥ κ M(w)?  Start with the corner dual at w = 1 (D = ½, corner cost ¼?).
- [ ] Proof-piece wishlist `search/WISHLIST.md` (10-03 picks: Lean `E(n) ≤ 0` + seam `D ≤ m_v`, then `M(1) ≤ ¾`, then `M(1) = ¾`).

## Exact forms (`tasks/exact-minpoly/`, `search/exact/EXACT_FORMS.md`)
- [ ] Target (10-07): exact + Lean `Packs` for all n ≤ 82 (open 29, 55, 68, 71; 83 deg 672), Lean local min ~70.  10-11: 269/324 exact, `Packs` 258, local min 178.
- [ ] Lean local minima: next n = 5 (grade B), one grade C (10 or 19).
- [ ] Exact forms for current records without one (72 on 10-10).  Superseded certificates stay (Evan 10-10: the viewer lists old packings).
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
- [ ] Arslanov rectangle decomposition in Lean (s(k²−k) < k for k ≥ 12?); read the paper first (E-JC doi:10.37236/8586).

## Outreach / site
- [ ] Per-n history viewer: what changed between consecutive records (matched-square diff, like the s(146) vs #505 page); design first.
- [ ] Draw every packing we have numbers for (`tasks/site-all-packings/`; replaces the `NEWER` stopgap in viewer.js).
- [ ] Repo unification: outreach → a small private GitHub repo; `wip/*` branches + pre-push hook if needed.
