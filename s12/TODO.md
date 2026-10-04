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

## Correctness / review
- [ ] Owed on jlevy (`notes/jlevy-s17-techniques.md` §4): #238 `s12/verify` exits 0 on NOT VERIFIED / partial runs, `cargo build | tail -1` hides failures; #256 stale s60 README paragraph, pin/print zmx2 source digest; #279 replay wand125's s59/s77 (he said yes); site s(17) floor R067 → R068.
- [ ] Second implementation of `Valid9` (F₄ k ≥ 8 rests on qx2 alone): deferred 10-03, Lean/simplification first; others may do it (k2m3 route: `ZMX2_AREA.md`).
- [ ] k2m3 2nd impl (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py should-fixes (audit 09-26): int64 guard on `sum(W·weight)`; `zm_d4_sweep.py summary` sha/settings gate; stale docstrings. Re-pin s(32) bundle.
- [ ] Watch: jlevy#316 (k2m4 registration), #311 (is our s(61) replay independent of zmx2?), PR #290, chelokot#18.

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` as an open statement in `Spec.lean` (check formal-conjectures first); axis-parallel ≤ 9 (item 1); Lemma CC chain counting (item 5, "Lean-ready").
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares as shared spec; offer bridges to chelokot (Frame), Queuingtheorydotcom. Draft → Evan approves.
- [ ] `ValidTilt7/9` in Lean (`notes/lean-leb-mass.md` §7): LEB + CAP done 10-03 (`LebMass.lean`). Next: `CovT` tree layer + generator (~400 lines), then PIECE with polygon clip (1–2k lines; with the tree, ~17k of V3's 32k leaves), then Lemma E (3–6k lines, weeks). k2m4 needs a leaf dump.
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY): Evan reviews → site page.

## s(17)
- [ ] s(17): pure closed covers are at their limit (exact ceiling ν_f(4.660) ≥ 17.0447, `CEILINGS_17_20.md`); further lower-bound gains need rule atoms (WISHLIST N15).
- [ ] Techniques from jlevy/Kleddamag/Guzhou (`notes/jlevy-s17-techniques.md` §2): rule atoms as LP columns, near-tight cells as exact LP rows (their re-solve of our s(12) points: ≥ 3.9702002), gmpy2 in exact checkers. s(17) exact: their PR #307 is far ahead; help with an independent check of its local theorem, not a restart.
- [ ] s(20) > 3+4√2/3 verified (`search/S20_LB.md`) but **superseded**: wand125 s(20) ≥ 1959/400 verified on jlevy (4.9 reported) already gives s(19) < s(20). Keep as independent confirmation (point cover, zmx2 non-1/10 sides); review running (`tasks/s20-review/`).
- [ ] s(19) > (7+√7)/2 ⇒ s(18) < s(19): **on hold** (wand125 at 1927/400, 0.11 % short; our pure-cover LP has ≈ 0.8 % room, CEILINGS §6). Offer the target on jlevy instead of racing.
- [ ] s(17) exact via the s(11) method: on hold (CPU). `tasks/s17-core-isolation/`

## Outreach / hygiene
- [ ] Announce: VibeMathed entries for s(21), s(32) (and now k²−3, k²−4 / F₄); X drafts `outreach/drafts-s21-s32.md`.
- [ ] `notes/status.md` frozen at 09-22: add k²−3, k²−4, s(21), s(32), s(45), s(60/61), Lean state.
- [ ] Repo unification, rest: outreach → a small private GitHub repo (versioned, backed up); `wip/*` branches + pre-push hook if needed.
- [ ] Prune ~30 `worktree-agent-*` branches on the private `square-packing-research` remote.
