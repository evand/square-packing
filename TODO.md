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
- [ ] Lean `Packs` for field degree > 20 (41, 51, 69, 87, 106, 128, 152, 177, 205, 266, 300): field arithmetic too heavy (n = 41 > 35 GB); needs a cheaper field representation.
- [ ] Target (Evan, 10-07): exact forms + Lean `Packs n S*` for all n ≤ 82 (open: 29, 55, 68, 71 = msolve-hard; 83 = deg-672 wall), Lean local min for ~70 of them (not 17/41/51/69).  State: 269/324 exact, `Packs` 258, local min 178 (`search/exact/EXACT_FORMS.md`).  Large n Lean: split cert into modules (B294/B310 die on per-file memory).
- [ ] Grade-C exact certificate (`grade-C.md`: paper proof done, numerics pass for 7 n).  Needs the second-order lemma too (hypothesis 3); build it on n = 5 first (grade B = no flat motions).
- [ ] Coverage, 55 open by cost: flat rotation (10), multi-direction force balance (11), block-triangular elimination (26 + msolve-hard 9), field 108/206.  msolve-hard: 29, 68, 71, 126, 228 are k multivariate equations in k class angles (n = 29: 5 eqs, degrees 12–26; field degree likely > 120), not a factor-choice problem (that fix landed 10-07).  Needs better elimination (block-triangular / triangular decomposition / other solver).  55: force-balance determinant deg 85, 5.7M terms; try multipliers as msolve variables.  Numerics question: task doc §Numerics.
- [ ] Explainer: local optimality of the n ≤ 324 records (for the s(11) explainer to cite). `notes/local-optimality-explainer-draft.md`.

## Major: exact results
- [ ] k²−5 for all large k: R = w = 5 κ = 0.02 D ≈ 1.30, conditional GO, margin ≈ 0.05 (`K2M4_MARGIN.md`); box 13 pitch 0.1. Make it cheaper first: cert golf (`tasks/cert-slack/`), finer-pitch D at w = 4 (insertable shortcut dead: price ≈ 1.8, `S2_INSERTABLE.md`).
- [ ] Seam capacity `M(∞)` (C5/R2; FRIEDMAN §8–9, `SEAM_W1.md`, `SEAM_1D.md`): `M(w) → ∞` ⇒ families give s(k²−c) = k for every c. Next: dual-side band LP w = 2, 3 (rigorous caps); wall-row + slow-variation lemmas.
- [ ] Corner coupling (R3): does D(w) ≥ κ M(w)? Start with the corner dual at w = 1 (D = ½, corner cost ¼?).

## s(12)
- [ ] Group-chat lemmas (10-05, low priority, `tasks/n12-chat-lemmas/` lists them): replay their evidence tarball with our checkers; Lean strip demo h = 1.13 and the corner–edge lemma `U+V ≥ 3.91c − 1.91`.

## Correctness / review
- [ ] Owed on jlevy (drafts: `outreach/draft-jlevy-replies-2026-10-07.md`; #238, #375, #256 posted 10-07): #279 s(59) exact replay **paused 10-07** (Evan: poor RoI on CPU for now; resume when cores are idle): one `zm_mixed` sweep at depth 34 over the whole D4 region (wand125 needed two runs whose joint coverage jlevy couldn't check); 33,725/102,400 roots done, ~127 CPU-h total, the rest is the expensive part.  Resume: `python3 search/zm_mixed.py cert runs/s59_wand125/n59_mixed_cover_8.txt --d4 --cert-mode --disj --chain-from 0 --depth 34 --pitch 1/20 --ubins 16 --nproc 15 --progress 5000 --resume runs/s59_wand125/zm_mixed_d4_depth34.jsonl --manifest runs/s59_wand125/zm_mixed_d4_depth34_manifest.json`; then reply (draft (D)).  Posting needs Evan's OK.
- [ ] (perf, scales) `lean/Sqpack/Bentz.lean` peaks at 23.5 GB / 126 s alone (next: S32 11.8 GB, ValidSplit9 9.8 GB, rest ≤ 7 GB): find the heavy term/tactic before the families grow.
- [ ] Format bridge: read wand125/tokoharu rectangle-density certificates in `zmx2` (converter or reader); first use: replay s(19) ≥ 193/40 (T-103), which alone settles s(18) < s(19) (jlevy may replay it too).
- [ ] `zmx2` exit status: still 0 on NOT VERIFIED / INCOMPLETE (verify, zmcheck changed 10-07); ~30 callers incl. bisect scripts to adapt first.
- [ ] k2m3 2nd impl (zmx2 area density): needs a 2nd reader for the 4 mutation-blind refinements (`ZMX2_AREA.md` §13).
- [ ] zeromargin.py should-fixes (audit 09-26; + s20 review 10-04: divisibility assert in `roots()` like `d4_roots`, else a non-dividing `--pitch` silently skips a strip): int64 guard on `sum(W·weight)`; `zm_d4_sweep.py summary` sha/settings gate; stale docstrings. Re-pin s(32) bundle.
- [ ] Watch: #399 (our s(266) ≤ 16.8230287508, s(270) ≤ 16.9378072285, s(272) ≤ 16.9681101458; not registered yet), #375 (n = 292), #419/#420 (exact forms, local minima; no reply yet), chelokot#18 (quiet since 10-02).

## Lean
- [ ] n = 12 batch (`notes/n12-gap.md` §2.2): `s(12) = 4` as an open statement in `Spec.lean` (formal-conjectures has only s(11), s(17): offer it upstream too); axis-parallel ≤ 9 (item 1); Lemma CC chain counting (item 5, "Lean-ready").
- [ ] Common spec (`lean/Sqpack/Spec.lean`): propose on jlevy/squares as shared spec; offer bridges to Queuingtheorydotcom (chelokot is moving proofs Lean → Bend, 10-07: re-scope). Draft → Evan approves.
- [ ] `ValidTilt7/9` in Lean (`notes/lean-leb-mass.md` §7): LEB + CAP done 10-03 (`LebMass.lean`). Next: `CovT` tree layer + generator (~400 lines), then PIECE with polygon clip (1–2k lines; with the tree, ~17k of V3's 32k leaves), then Lemma E (3–6k lines, weeks). k2m4 needs a leaf dump.
- [ ] s(21) Lean on cand A (`search/golf/candA_verified/`, ~51 CPU-h); levers: SPLIT leaf for 8 hot cells, per-leaf overhead.
- [ ] Arslanov rectangle decomposition in Lean (squeezable rectangles glue to side < m; s(k²−k) < k for k ≥ 12?).  Read the paper first (2019 E-JC doi:10.37236/8586; 2021 secondhand).  None in Lean anywhere (10-05).
- [ ] Lean explainer (`notes/lean-explainer-draft.html`, artifact FJPhwWoQLhdkcerwh4VBwY): Evan reviews → site page.

## s(17)
- [ ] s(17): pure closed covers are at their limit (exact ceiling ν_f(4.660) ≥ 17.0447, `CEILINGS_17_20.md`); further lower-bound gains need rule atoms (WISHLIST N15).
- [ ] Techniques from jlevy/Kleddamag/Guzhou (`notes/jlevy-s17-techniques.md` §2): rule atoms as LP columns, near-tight cells as exact LP rows gmpy2 in exact checkers.  s(17) exact: their PR #307 is far ahead and wand125 built its independent check (jlevy PR #410): not ours to do.
- [ ] s(17) exact via the s(11) method: on hold (CPU). `tasks/s17-core-isolation/`

## Packings (active, 10-08)
Engine `search/packer/`: `fq` (Rust quench: lifted-separator ALM + face-branch SLP polish) + `explore.py` (quality-diversity basin
explorer, Thompson move selection with global-novelty reward, elite / frontier budgets) + `explore_exact.py` / `arm_yield.py` /
`lineage.py`.  Log + night summary: `search/packer/PACKER.md` (end).  s(110) census: 496 certified sub-11 minima (`runs/known110_all.json`).
- [ ] SQP polish (rotation curvature, Clarabel QP): large-n stage 2 is ~56 SLP iterations × 2 LPs, 91–94 % of explorer CPU at n = 237.
- [ ] s(110) survey write-up (census 496, funnel structure, discovery curves per method, per-arm yields); a different sub-11 funnel needs a new constructor.  Cleanup: delete untracked `search/packer/seeds110/` (3 useless files; rm was blocked).
Decisions (10-08): compute serves testing (catalogue gains incidental); every comparison gets a concurrent control and a frozen novelty ref; certify before counting (~20–25 % of polished sub-k candidates are not minima); bandit reward = global novelty; elite share on for record hunting, off for census; resume a census run only on a new record.
- [ ] Start generation (next category step; each new method so far ~doubled the reachable set): SQUISH-style chaining (each new best seeds n ± 1, n ± 2), grafts of small-n records into patches of large ones, `layout.py` channel constructions at other angles finished with fq, crops of other non-grid records, fixed-side shrink from random starts.  Judge by new certified groups/basins per CPU vs a frozen ref and on the SQUISH benchmark seeds (`bench/seeds/`).
- [ ] Large-n throughput: LP solves (~0.07 s each at n = 237, ~92 % of explorer CPU); try a warm-startable solver (dual simplex / active set) for the SLP polish.  Replicates before further policy changes (tx237n vs tx237s differ by 3e-3).
- [ ] SQUISH gap: explorer + neighbour seeds (s108 seed: 98 % in 20 min; 237 from 238 − 1 still descending at 30 min); record starts don't move.
- [ ] Record hunting, cheap: re-run the wide sweep n = 30–323 at slp2 budget 180 s (13 % of 60 s quenches at n ≥ 250 unconverged); kick every older packing we have (mirror, register histories, our off-record minima).
- [ ] Is the plateau → funnel barrier entropic (first-order)?  If so temperature can't open it: bias / multicanonical in side, or gap-spanning moves (lessons from an earlier MH/PT project).
- [ ] Landscape measure: nested sampling / splitting on side from random starts, log X(s) at 110 (content-agnostic "how hard is this n"; `lit-statmech.md` A2).
- [ ] s(90): needle estimate from the 110 calibration (old step 4) before more compute; best so far 10.0095668 (certified local min).
- [ ] slp2 speed: contact model vectorised 10-06 (bitwise identical, ~1.6x); now ~500 LP solves per quench dominate (40–110 s at n ≈ 270): fewer LPs per descent, warm starts, or Rust.
- [ ] exactsolve should snap near-axis corner-loaded squares itself (the n = 292 failure, `search/packer/s292.md`); 12/866 cen7 outputs unresolved.
- [ ] Packer at large n: default 2 % loosen kicks the unkicked s(292) record out of its basin (loosen 1.0 fine): check other large n, re-run sw2/sw3 there if needed; sweep 291 (now has a 292 input).
- [ ] Case study `search/packer/s110-landscape.md`: rewrite the census section from cen7 + rev1 + minima PT.
- [ ] Side bet: rectangle containers in slp2; Arslanov 26-in-(4−δ)×8 control, then jlevy H-049 (20 in (4−δ)×6).  Refresh clone first.
- [ ] Later targets for the same instrument: s(183/242/274/308), 299, 234, s(147) (`search/WISHLIST.md` §P).
- [ ] Tiling starts: **paused (Evan, 10-08)** after the first batch (`search/TILINGS.md`): all three doublings beaten, certified s(964) ≤ 31.8595, s(1092) ≤ 33.8406, s(1228) ≤ 35.8474.  The packings are not local minima (exactsolve: corner-corner descent ≈ 0.05 at 964); optimise them only after the optimizer work.  `exactsolve.py` crashes at default tol on fq output at n ≈ 1000 (`max_support` linprog, empty contact set?); `--tol 1e-7` works.

## Outreach / hygiene
- [ ] Site: draw every packing we have numbers for (Couzo's 49, itsnaka's 23 SQUISH (local copies + pins: `tasks/site-all-packings/README.md`), de Winter's 211, ours 266/270/272, n > 324; replaces the `NEWER` stopgap in viewer.js); coords with credit, no code/text. `tasks/site-all-packings/`
- [ ] Announce: VibeMathed entries for s(21), s(32) (and now k²−3, k²−4 / F₄); X drafts `outreach/drafts-s21-s32.md`.
- [ ] Repo unification, rest: outreach → a small private GitHub repo (versioned, backed up); `wip/*` branches + pre-push hook if needed.
