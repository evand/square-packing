# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths relative to the repo root (`public/`; flattened from `public/s12/` on 10-07).  Public since 10-03 (outreach stays in `private/s12/outreach`, symlinked as `outreach/`).  Done: `Completed.md`.  Conventions: `~/math/TODO.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`) for current values.**

## Next (Evan, 10-07)
Cleanup session 10-07 done (Completed.md): flatten s12/ → root, exit codes, upstream sync, site refresh.  Next: optimizer work; exact-forms route paused after jlevy#419/#420.  Git: single checkout on `main`; exact-minpoly raw solve outputs (gitignored) in `search/exact/minpoly/solve/`.

## Next (picked 2026-10-03)
Proof-piece wishlist (new candidates N1–N15 + 65 open items, picks at the end): `search/WISHLIST.md`.
1. Lean, small and interesting: `E(n) ≤ 0` + seam bound `D ≤ m_v`, then `M(1) ≤ ¾` (FRIEDMAN §9.6 item 5).
2. `M(1) = ¾` fully proved: SEAM_W1 §2.3 casework (or interval arithmetic).

## Exact forms (10-07; `tasks/exact-minpoly/`)
First deliverable: a table of exact forms of S + a handful of Lean local-optimality proofs.  n = 83 (deg 672) optional.
- [ ] Exact forms batch (10-07): 269/324 exact (4 new: 102, 106, 152, 177), Lean `Packs` 258, Lean local min 178; `search/exact/EXACT_FORMS.md`.  Posted 10-07: jlevy#419 (exact forms), #420 (local minima), #375 comment (s(292)): watch for import.
- [ ] Lean local minima: s(11), s(28) done (`Exact/N11L`, `N28L`; sparse G check 10-07); next n = 5 (grade B), one grade C (10 or 19).
- [ ] Lean `Packs` for field degree > 20 (41, 51, 69, 87, 106, 128, 152, 177, 205, 266, 300): field arithmetic too heavy (n = 41 > 35 GB); needs a cheaper field representation.
- [ ] Target (Evan, 10-07): exact forms + Lean `Packs n S*` for all n ≤ 82 (open: 29, 55, 68, 71 = msolve-hard; 83 = deg-672 wall), Lean local min for ~70 of them (not 17/41/51/69).  Large n Lean: split cert into modules (B294/B310 die on per-file memory accumulation).
- [ ] (low, curiosity) n = 10 second order: Göbel's 10a record has dim null(J_A) = 4 on the two 45° squares (2 exact flat, 2 PD eig 0.339): what are the coupled motions?  10b (one 45° square, contacts at all 4 edge midpoints) is the n = 5 motif.  Not needed: s(10) global via Stromquist / chelokot bridge.
- [ ] Grade-C exact certificate (`grade-C.md`: paper proof done, numerics pass for 7 n).  Needs the second-order lemma too (hypothesis 3); build it on n = 5 first (grade B = no flat motions).
- [ ] Coverage, 59 open by cost: flat rotation (10), multi-direction force balance (11), block-triangular elimination (26 + msolve-hard 9), field 108/206; 292 (packer agent 10-07).
- [ ] msolve-hard: 29, 68, 71, 126, 228 are k multivariate equations in k class angles (n = 29: 5 eqs, degrees 12–26; field degree likely > 120), not a factor-choice problem (that fix landed 10-07).  Needs better elimination (block-triangular / triangular decomposition / other solver).  55: force-balance determinant deg 85, 5.7M terms; try multipliers as msolve variables.  Numerics question: task doc §Numerics.
- [ ] Explainer: local optimality of the n ≤ 324 records (for the s(11) explainer to cite). `notes/local-optimality-explainer-draft.md`.

## Major: exact results
- [ ] k²−5 for all large k: R = w = 5 κ = 0.02 D ≈ 1.30, conditional GO, margin ≈ 0.05 (`K2M4_MARGIN.md`); box 13 pitch 0.1. Make it cheaper first: cert golf (`tasks/cert-slack/`), finer-pitch D at w = 4 (insertable shortcut dead: price ≈ 1.8, `S2_INSERTABLE.md`).
- [ ] Seam capacity `M(∞)` (C5/R2; FRIEDMAN §8–9, `SEAM_W1.md`, `SEAM_1D.md`): `M(w) → ∞` ⇒ families give s(k²−c) = k for every c. Next: dual-side band LP w = 2, 3 (rigorous caps); wall-row + slow-variation lemmas.
- [ ] S1 measurement (optional): full insertable LP class (b) at k = 6, 7; is the price ≈ 1.8 k-independent? ~1–2e4 CPU-s (`S2_INSERTABLE.md` §5).
- [ ] Corner coupling (R3): does D(w) ≥ κ M(w)? Start with the corner dual at w = 1 (D = ½, corner cost ¼?).

## s(12): group-chat lemmas (10-05, `tasks/n12-chat-lemmas/`; low priority, not a route change)
- [ ] Ask the group for `n12_progress_evidence_2026-10-04.tar.gz`; replay wall-strip h = 1.2071 (8.2M leaves), p37_b24, E2 point certs with our checkers.
- [ ] Lean: strip demo h = 1.13 (`strip-demo/PROOF.md` §8: 26 boxes, 234 rational leaves + ordering/chain lemmas).
- [ ] Lean: their corner–edge lemma `U+V ≥ 3.91c − 1.91` (side-free, holds at t = 4; analytic + one Bernstein quartic).
- [ ] Counting route (README Outcome): partitions the zero set into 18 margin-0 classes, doesn't shrink it. Only (A) z ≥ 5, (B) 3 edge/wall have room.

## Correctness / review
- [ ] Owed on jlevy (drafts: `outreach/draft-jlevy-replies-2026-10-07.md`; #238, #375 posted 10-07): #256 confirm s60 re-run reading (README fixed in cca7bf1); #279 drop or run wand125's s59/s77 replay; #419 superseded-records note.  Each needs Evan's OK.
- [ ] Upstream data to pull (10-07 sync): itsnaka SQUISH certs (18 better packings, `itsnaka/squish-certs`, jlevy#401); Couzo HEAD 6042c56b (n = 208, 209, 228, 263, 272, 303, 306); re-pin `_untrusted-third-party` clones (all behind; wand125's repos merged into `wand125/square-packing`).
- [ ] k2m3 wording: `docs/k2m3.html` ("one exact checker" tag, status) and `docs/s32.html:147` still say single-checker; wand125's Valid7 check exists: check the register status, then reword.
- [ ] Format bridge: read wand125/tokoharu rectangle-density certificates in `zmx2` (converter or reader); first use: replay s(19) ≥ 193/40 (T-103), which alone settles s(18) < s(19) (jlevy may replay it too).
- [ ] `zmx2` exit status: still 0 on NOT VERIFIED / INCOMPLETE (verify, zmcheck changed 10-07); ~30 callers incl. bisect scripts to adapt first.
- [ ] Second implementation of `Valid9`: wand125 did one (ValidTilt9, `wand125/valid7-independent-check` c561dbb3, 10-06; register: reported, T-081 stays V0/C1).  Read it, then update site/READMEs that say "no second implementation".
- [ ] k2m3 2nd impl (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py should-fixes (audit 09-26; + s20 review 10-04: divisibility assert in `roots()` like `d4_roots`, else a non-dividing `--pitch` silently skips a strip): int64 guard on `sum(W·weight)`; `zm_d4_sweep.py summary` sha/settings gate; stale docstrings. Re-pin s(32) bundle.
- [ ] Watch: #399 (266/270/272 not registered yet), #375 (n = 292 not registered), #419/#420 (no reply yet), chelokot#18 (quiet since 10-02).

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` as an open statement in `Spec.lean` (formal-conjectures has only s(11), s(17): offer it upstream too); axis-parallel ≤ 9 (item 1); Lemma CC chain counting (item 5, "Lean-ready").
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares as shared spec; offer bridges to Queuingtheorydotcom (chelokot is moving proofs Lean → Bend, 10-07: re-scope). Draft → Evan approves.
- [ ] `ValidTilt7/9` in Lean (`notes/lean-leb-mass.md` §7): LEB + CAP done 10-03 (`LebMass.lean`). Next: `CovT` tree layer + generator (~400 lines), then PIECE with polygon clip (1–2k lines; with the tree, ~17k of V3's 32k leaves), then Lemma E (3–6k lines, weeks). k2m4 needs a leaf dump.
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Arslanov rectangle decomposition in Lean (squeezable rectangles glue to side < m; s(k²−k) < k for k ≥ 12?).  Read the paper first (2019 E-JC doi:10.37236/8586; 2021 secondhand).  None in Lean anywhere (10-05).
- [ ] Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY): Evan reviews → site page.

## s(17)
- [ ] s(17): pure closed covers are at their limit (exact ceiling ν_f(4.660) ≥ 17.0447, `CEILINGS_17_20.md`); further lower-bound gains need rule atoms (WISHLIST N15).
- [ ] Techniques from jlevy/Kleddamag/Guzhou (`notes/jlevy-s17-techniques.md` §2): rule atoms as LP columns, near-tight cells as exact LP rows (their re-solve of our s(12) points: ≥ 3.9702002), gmpy2 in exact checkers. s(17) exact: their PR #307 is far ahead; wand125 already built an independent check (jlevy PR #410, `wand125/square-packing` n17_kernel_verifier), so low priority for us.
- [ ] s(17) exact via the s(11) method: on hold (CPU). `tasks/s17-core-isolation/`

## Packings (active, 10-06)
Engine `search/packer/` (map + plan: `README.md` "Plan (10-06)"; log: `PACKER.md`; merged to main 10-06).  s(110) calibration
done enough: sub-11 funnel small and well sampled, plateau above 11 rough; the needle is entering the funnel.  Now: samplers that go uphill.
- [ ] Minima PT (`mcmin.py`): runs A and B at 110 both missed (A: hot replicas sink into the grid; B: 11.0076 funnel deep, staggered-chain trap at 11).  Temperature alone doesn't connect funnels: next, coordinated moves (vacancy → hole chain shift), bias in side.
- [ ] New moves (vacancy → hole chain shift first): judge on `abtest.py` at 110 and on the 20 rediscovery pairs at distance > 0.06 that kicks never re-find (PACKER "Rediscovery result").
- [ ] Record hunting, cheap: re-run the wide sweep n = 30–323 at slp2 budget 180 s (13 % of 60 s quenches at n ≥ 250 unconverged); kick every older packing we have (mirror, register histories, our off-record minima).
- [ ] Is the plateau → funnel barrier entropic (first-order)?  If so temperature can't open it: bias / multicanonical in side, or gap-spanning moves (lessons from an earlier MH/PT project).
- [ ] New records 10-06: s(266) ≤ 16.8230287508, s(270) ≤ 16.9378072284, s(272) ≤ 16.9681101458 (certified, `search/exact/results/sw2/`); posted jlevy#399 (10-06): watch for registration.
- [ ] (parked, low) p(δ): "a δ-better packing would have been found with probability p" from mixing/ESS on the no-grid measure, basin distance, δ-vs-weight; motivates mixing across the grid side (PACKER 10-06).
- [ ] Landscape measure: nested sampling / splitting on side from random starts, log X(s) at 110 (content-agnostic "how hard is this n"; `lit-statmech.md` A2).
- [ ] s(90): needle estimate from the 110 calibration (old step 4) before more compute; best so far 10.0095668 (certified local min).
- [ ] slp2 speed: contact model vectorised 10-06 (bitwise identical, ~1.6x); now ~500 LP solves per quench dominate (40–110 s at n ≈ 270): fewer LPs per descent, warm starts, or Rust.
- [ ] Exact batch (`search/exact/`, jlevy#375): 323/324 + n = 292 certified 10-07 (S' = 17.5972493911564650406, −4.31e-11, Couzo's packing; `search/packer/s292.md`): add to `batch/` (→ 324/324); exactsolve should snap near-axis corner-loaded squares itself.  Leftovers: 12/866 cen7 outputs unresolved.
- [ ] Packer at large n: default 2 % loosen kicks the unkicked s(292) record out of its basin (loosen 1.0 fine): check other large n, re-run sw2/sw3 there if needed; sweep 291 (now has a 292 input).
- [ ] Case study `search/packer/s110-landscape.md`: rewrite the census section from cen7 + rev1 + minima PT.
- [ ] Side bet: rectangle containers in slp2; Arslanov 26-in-(4−δ)×8 control, then jlevy H-049 (20 in (4−δ)×6).  Refresh clone first.
- [ ] Later targets for the same instrument: s(183/242/274/308), 299, 234, s(147) (`search/WISHLIST.md` §P).
- [ ] Reference hygiene: before claiming a packing, check jlevy's register and grep the local clone `frontier/n-<n>.md` (refresh per `_untrusted-third-party/PROVENANCE.txt`).
- [ ] (parked) contact-graph counts (force-bearing network), pose → closed form beyond `exact/` (`lit-provenance.md`).

## Outreach / hygiene
- [ ] Site `problems.html` is the reader-facing copy of WISHLIST's object-level items (§P questions, A1/A3, N16): edit both together.
- [ ] Site: draw every packing we have numbers for (Couzo's 49 for n = 68–307, de Winter's n = 211, n > 324); coords are fair game, credit + link, copy no code/text. `tasks/site-all-packings/`
- [ ] Announce: VibeMathed entries for s(21), s(32) (and now k²−3, k²−4 / F₄); X drafts `outreach/drafts-s21-s32.md`.
- [ ] Repo unification, rest: outreach → a small private GitHub repo (versioned, backed up); `wip/*` branches + pre-push hook if needed.
- [ ] Prune ~30 `worktree-agent-*` branches on the private `square-packing-research` remote.
