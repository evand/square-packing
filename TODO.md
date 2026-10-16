# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths from the repo root.  Done: `Completed.md`.  Conventions: `~/math/TODO.md` (one line per entry, details in `tasks/`).
Pre-tidy text (10-15, has details dropped here): `tasks/todo-tidy-2026-10-15/TODO-before.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`).**

## Waiting on Evan
- [ ] Review the Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY) → site page.
- [ ] Compute OKs: s(29) grade-C (`search/exact/homotopy/N29.md` §Pickup); 17 grade-C timeouts uncapped; regularize the 16 raw records; #279 replay (`tasks/s59-replay/`).
- [ ] Decide: n = 88 strip-family test set (~16 CPU-h; TT2 is done, results not 110-shaped → drop?).
- [ ] Posts: VibeMathed s(21), s(32), k²−3, k²−4 (`outreach/drafts-s21-s32.md`); variant lists to jlevy (after the write-up below).
- [ ] Time-lapse records video (`site/timelapse/`): Evan to return; store sighting times are local, not UTC.

## Watch
- [ ] jlevy #659 (our s(236) ≤ 15.861934390608, posted 10-16, cert 5ad0069), #639 (exact forms, 22 current records), #606 (our 146: T-187 V0/C0 10-15, awaiting their replay review), #545 (s29 / grade-C retained 10-15; their replays wait on budget; they note exact-point non-overlap for s(29) open), #596, #584, #552, #489, #465, #419 / #420 (exact forms, no reply), #375, #3 (preprint, left open); chelokot#18 (repo unknown: not under bend / bend-math).

## Packings (`search/packer/`: `HANDOFF-2026-10-14.md`, `tasks/packer-arch-2026-10-14/`)
Thesis (Evan 10-14): diversity pays; local search at current budgets has mined the frontier.  No compute above n = 324 without a specific reason.
- [ ] TT3: gating PASSED (C42, 2.08x); melt / GLS v0 neutral (C43); parent value model PASSED on the 12-round subset (C44, ~1.7-2x over G) and on the full 49 rounds (PF49 10-15: 2.81x [2.22, 3.45]); no transfer at 110 (absolute gap bins).  Next: a relative-gap model (110), then the all-n hunt-style test (new minima + regions + records).
- [ ] Baseline on the major tools (Evan 10-15): per family (A composition, B soft-to-rigid, F/G removal, J constructors) level 1 = reproduce events from named parents, level 2 = blind time travel; improvements only after.  `tasks/composition-baseline/`, PACKER "Baseline on the major tools".
  - Status 10-15 19:00: rectangle blocks 5 / 6 blind (C45); Mondrian v1 (#481) L1 4 / 5 named (C46), **L2e blind 4 / 15 BEAT-ANTE** (exact ranking) -> live hunt next; #500 0; removal 2 / 12; graft 3 / 6; J 0 / 10 by nearby search.  **Next steps ranked: `tasks/packer-next-2026-10-15/`** (paused for Evan's review 10-15).
  - J generator: EV1 0 / 18 (superseded); EV2 (`search/packer/evolve2/`, LOG.md): gen 0 reproduced 178 / 126 / 145 / 104 / 109 (record shown); gen 1 on the hidden row 234-240 tied 234 + 235 blind, 236 +0.0445, 237 +0.101, 238-240 at the grid.  Gen 2 plan in `tasks/packer-next-2026-10-15/` §2.
- [ ] Live slack-composition hunt over n <= 324 (L2e rule fired); first a cheaper compose screen (~6 CPU-h / n now).  `tasks/packer-next-2026-10-15/` §1.
- [ ] EV2 quench budget: the evaluator quench (320-1,093 s / eval) is now the bottleneck (fill solved, FP1); test seeds / loosens / candidate cap / procs.
- [ ] Explorer seeded from EV2 gen 1's 55 novel near-record minima (~5 CPU-h).  `tasks/packer-next-2026-10-15/` §3.
- [ ] s(110) corner-cell holdouts (census 1,078 minima, 19:12 10-15) + a stop-early descent mode (GX1: crossings collapse to the record basin).  §4.
- [ ] numba for `fillfast`: asked 10-15, not approved, low value now.
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
- [ ] s(110) weak residue: 12 side keys / 10 minima "negative curvature blocked" (cone-QP check, then branch-aware second order); exactsolve fallback fixed 10-15 (PACKER "exactsolve fix").
- [ ] Shallow-valley descender X01: `search/packer/nlpdesc.py` (second-order fixed-line NLP, IPOPT) resolved 34 / 40 unresolved s(110) keys 10-15 (PACKER "unresolved / weak: second-order descent"); next: the 225 stuck points, the 6 left (MILP-branch kick).
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
- [ ] Arslanov rectangle decomposition in Lean (s(k²−k) < k for k ≥ 12?); read the paper first (E-JC doi:10.37236/8586).

## Outreach / site
- [ ] Per-n history viewer: what changed between consecutive records (matched-square diff, like the s(146) vs #505 page); design first.
- [ ] Draw every packing we have numbers for (`tasks/site-all-packings/`; replaces the `NEWER` stopgap in viewer.js).
- [ ] Repo unification: outreach → a small private GitHub repo; `wip/*` branches + pre-push hook if needed.
