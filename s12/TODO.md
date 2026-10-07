# TODO — square packing

Goals (Evan): **major** an exact result for something open; **secondary** s(17) lower bound, Lean, packings (low expectation).
Priorities (10-01): correctness first, then interesting Lean, then cleverness over compute; others building on our tools is a win.
Paths relative to `public/s12`.  Public since 10-03 (outreach stays in `private/s12/outreach`).  Done: `Completed.md`.  Conventions: `~/math/TODO.md`.
External posts: ask Evan explicitly per post.  **Before any bound-chasing run, check jlevy's register (`jlevy.github.io/squares/cases/<n>.html`) for current values.**

## Next (picked 2026-10-03)
Proof-piece wishlist (new candidates N1–N15 + 65 open items, picks at the end): `search/WISHLIST.md`.
1. Lean, small and interesting: `E(n) ≤ 0` + seam bound `D ≤ m_v`, then `M(1) ≤ ¾` (FRIEDMAN §9.6 item 5).
2. `M(1) = ¾` fully proved: SEAM_W1 §2.3 casework (or interval arithmetic).

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
- [ ] Owed on jlevy (`notes/jlevy-s17-techniques.md` §4): #238 `s12/verify` exits 0 on NOT VERIFIED / partial runs, `cargo build | tail -1` hides failures; #256 stale s60 README paragraph, pin/print zmx2 source digest; #279 replay wand125's s59/s77 (he said yes).
- [ ] Second implementation of `Valid9` (F₄ k ≥ 8 rests on qx2 alone): deferred 10-03, Lean/simplification first; others may do it (k2m3 route: `ZMX2_AREA.md`).
- [ ] k2m3 2nd impl (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py should-fixes (audit 09-26; + s20 review 10-04: divisibility assert in `roots()` like `d4_roots`, else a non-dividing `--pitch` silently skips a strip): int64 guard on `sum(W·weight)`; `zm_d4_sweep.py summary` sha/settings gate; stale docstrings. Re-pin s(32) bundle.
- [ ] Watch: jlevy#281 (s(19) target offered to wand125, 10-04), jlevy#316 (k2m4 registration; when T-081 is verified: Overview `FAMILIES` += {4: 5}, `register: verified` on n = 96), #311 (is our s(61) replay independent of zmx2?), PR #290, chelokot#18, #375 (analytic batch).

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` as an open statement in `Spec.lean` (formal-conjectures has only s(11), s(17): offer it upstream too); axis-parallel ≤ 9 (item 1); Lemma CC chain counting (item 5, "Lean-ready").
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares as shared spec; offer bridges to chelokot (Frame), Queuingtheorydotcom. Draft → Evan approves.
- [ ] `ValidTilt7/9` in Lean (`notes/lean-leb-mass.md` §7): LEB + CAP done 10-03 (`LebMass.lean`). Next: `CovT` tree layer + generator (~400 lines), then PIECE with polygon clip (1–2k lines; with the tree, ~17k of V3's 32k leaves), then Lemma E (3–6k lines, weeks). k2m4 needs a leaf dump.
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Arslanov rectangle decomposition in Lean (squeezable rectangles glue to side < m; s(k²−k) < k for k ≥ 12?).  Read the paper first (2019 E-JC doi:10.37236/8586; 2021 secondhand).  None in Lean anywhere (10-05).
- [ ] Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY): Evan reviews → site page.

## s(17)
- [ ] s(17): pure closed covers are at their limit (exact ceiling ν_f(4.660) ≥ 17.0447, `CEILINGS_17_20.md`); further lower-bound gains need rule atoms (WISHLIST N15).
- [ ] Techniques from jlevy/Kleddamag/Guzhou (`notes/jlevy-s17-techniques.md` §2): rule atoms as LP columns, near-tight cells as exact LP rows (their re-solve of our s(12) points: ≥ 3.9702002), gmpy2 in exact checkers. s(17) exact: their PR #307 is far ahead; help with an independent check of its local theorem, not a restart.
- [ ] s(19) > (7+√7)/2 ⇒ s(18) < s(19): **on hold** (wand125 at 1927/400, 0.11 % short; our pure-cover LP has ≈ 0.8 % room, CEILINGS §6). Offer the target on jlevy instead of racing.
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
- [ ] Exact batch (`search/exact/`, jlevy#375): 323/324 certified (n = 105, 130 added 10-05, pushed 10-06; comment on #375 owed: ask Evan); open n = 292.  Leftovers: 12/866 cen7 outputs unresolved.
- [ ] Case study `search/packer/s110-landscape.md`: rewrite the census section from cen7 + rev1 + minima PT.
- [ ] Side bet: rectangle containers in slp2; Arslanov 26-in-(4−δ)×8 control, then jlevy H-049 (20 in (4−δ)×6).  Refresh clone first.
- [ ] Later targets for the same instrument: s(183/242/274/308), 299, 234, s(147) (`search/WISHLIST.md` §P).
- [ ] Reference hygiene: before claiming a packing, check jlevy's register and grep the local clone `frontier/n-<n>.md` (refresh per `_untrusted-third-party/PROVENANCE.txt`).
- [ ] (parked) contact-graph counts (force-bearing network), pose → closed form beyond `exact/` (`lit-provenance.md`).

## Outreach / hygiene
- [ ] Site `problems.html` is the reader-facing copy of WISHLIST's object-level items (§P questions, A1/A3, N16): edit both together.
- [ ] Site: draw every packing we have numbers for (Couzo's 49 for n = 68–307, de Winter's n = 211, n > 324); coords are fair game, credit + link, copy no code/text. `tasks/site-all-packings/`
- [ ] Announce: VibeMathed entries for s(21), s(32) (and now k²−3, k²−4 / F₄); X drafts `outreach/drafts-s21-s32.md`.
- [ ] `notes/status.md` frozen at 09-22: add k²−3, k²−4, s(21), s(32), s(45), s(60/61), Lean state.
- [ ] Repo unification, rest: outreach → a small private GitHub repo (versioned, backed up); `wip/*` branches + pre-push hook if needed.
- [ ] Prune ~30 `worktree-agent-*` branches on the private `square-packing-research` remote.
