# Site review, 2026-09-28: proofs.html vs /s12/, stale facts, Lean status

Not published: `pages.yml` copies only `s12/docs/{index,s21,s32}.html` from this directory.
Builds on `site/notes/site-review-2026-09-22.md` (its §3 "division of labour" still stands; this
memo updates it for s(21), s(32), s(45) and the Lean ladder).

## 1. What each page is for (as written today)

| page | audience | job |
|---|---|---|
| `proofs.html` | general, curious | how a lower bound is proved: the unavoidable-points trick, the 81-point playground, why 12 is open, the Bentz s(13) walkthrough, the case-free s(13) proof |
| `/s12/` | technical | write-up of s(12) ≥ 3.968616, s(11) ≥ 3.8143 and case-free s(13) = 4: result, context, certificate, verification, limits, the s(13) checkers, artifacts |
| `/s21/`, `/s32/` | technical | one result each, fixed 7-section template (result, context, certificate/idea, verification, Lean, caveats, reproduce) |

`/s12/` predates the template and carries two results (s(12)/s(11), and s(13) in §07).

## 2. Overlap, section by section

| content | proofs.html | /s12/ | should own it |
|---|---|---|---|
| weighted unavoidable set → counting argument | "The trick" | §01 ¶2 | proofs (idea); /s12/ one sentence + link |
| 81 points, 7 each, s(12) ≥ 35/9 | playground (interactive) | §01 card, §04 static SVG | proofs (playground); /s12/ keeps the SVG as a no-JS fallback, cut its caption |
| why 12 is open / LP-duality ceiling / 12.2688 | "Why twelve" (1 ¶) | §06 (full) | /s12/; proofs keeps its 1 ¶ |
| Bentz six-leaf walkthrough | yes | — | proofs |
| case-free s(13): theorem, closed convention, 12.2688, margin zero | "Thirteen again" (4 ¶) | §07 (theorem, 3 subsections, table, transcript) | a dedicated s(13) page (see §5); proofs keeps 1 ¶ |
| two checkers, rejection tests, Lean | "Thirteen again" ¶4 | §07 | /s12/ §07 (or s(13) page) only |
| s(21)/s(32)/s(45) | not mentioned | not mentioned (except "none of this moves s(12)") | own pages; proofs should get 1 ¶ each |

## 3. Factual inconsistencies / stale content

Fixed now (see §6): all Lean-status statements.  **Still open (need Evan):**

1. **s(11) is no longer our best claim.**  `notes/literature-s32.md` §4–5 and `s12/TODO.md`: Kleddamag
   s(11) > 31/8 = 3.875 (Sept 2026, replayed by jlevy; a *lower* bound, just under Trump's 3.877084
   packing), jlevy s(11) ≥ 3.827.  Our 3.8143 is below both.  Stale on: `/s12/` lede ("a new floor for
   eleven squares"), §02 table ("bound raised here to 3.8143"), `sources.html` §5, and
   `lower_bounds.json` n = 11 `best` (ours).  Needs a verified replay before the Atlas adopts 3.875.
2. **s(12) "previous published bound 3.788854" (`/s12/` §01, §02, lede "had stood since 2003").**
   jlevy T-017 s(12) ≥ 99/25 = 3.96 (2026-09-04, unrefereed).  Our 3.9686 is dated 2026-08-26 in
   `lower_bounds.json`, so it may not "predate" ours — but it should be cited; the `/s12/` caveat "if a
   better published bound exists, we did not find it" is now misleading.
3. **s(17).**  `sources.html` §5 and `lower_bounds.json` n = 17: Mira 4.613028635886 (2026-09-08).
   `s12/TODO.md`: Mira 4.6141535 (2026-09-20, 3,280 atoms).  `notes/literature-s32.md`: Guzhou0806 and
   Kleddamag s(17) > 4.6200.  Bounds chart and Sources are both behind.
4. **Uncited 2026 work.**  `sources.html` does not cite chelokot/square-packing-archive (Lean,
   kernel-checked s(n²−2) = n, s(13), s(22), s(33), s(46) …), tokoharu, Guzhou0806, Kleddamag.  chelokot
   matters for our wording: s(13) = 4 was already kernel-checked in Lean (via a different proof); we do
   not claim first, but Sources should say so.
5. **s(45) = 7 has no page and no mention on proofs/`/s12/`.**  Only `sources.html`, the bounds data and
   the repo README.  Site index has no card for it.
6. Repo READMEs (outside my edit scope) still say s(32) is "from that computational hypothesis" and
   that the Lean covers only the primitives for s(13): `README.md` (repo root) row s32/s12,
   `s12/README.md` lines ~40, 115, 358–370, `s12/certificates/s32/README.md` lines 26, 69.
7. Voice: `/s12/` and `/s21/`/`/s32/` say "we"; OK now.  (The 09-22 review's "I" issue is fixed.)
8. Known and unfixed from 09-22: bounds.js prints the year twice in the floor card; nav overflow on
   phones (9 items now); `/s12/` context table overflows at 360 px.

No broken internal links or anchors (checked every `href`/`src` in the assembled site, and every
`github.com/evand/square-packing/{tree,blob}/main/…` link against `git ls-files`).

## 4. Lean status (for reference; now reflected on the pages)

| result | Lean |
|---|---|
| s(13) = 4 | kernel-checked, no hypothesis (`S13Lower.lean`, opt-in) |
| s(32) = 6 | kernel-checked, no hypothesis (`S32Lower.lean`, opt-in, 13.8 CPU-h) |
| s(21) = 5 | from the checker hypothesis (`S21.lean`); segments not yet in the verifier |
| s(45) = 7 | none |
| s(11) ≥ 3040/797; s(12) ≥ 35/9, ≥ 3920/997 | kernel-checked |
| s(12) ≥ 15680/3951 | reduction only; kernel check in progress (marked `<!-- LEAN-3.968616 … -->` in `/s12/` ×2, `proofs.html`, `index.html`, `sources.html`) |

## 5. Recommended restructure (not implemented)

1. **proofs.html = the cross-result methods page.**  Keep the trick, the playground and the Bentz
   walkthrough (the site's best explainers).  Then a short "What the method has proved since" section:
   one paragraph + link per result — s(12)/s(11) floors, s(13) case-free, s(21) and s(45) (line
   densities), s(32) — plus one paragraph on checking (exact checkers, margin zero, what Lean's kernel
   now checks).  Cut "Thirteen again" to ~90 words (the 09-22 rewrite) and link out.
2. **`/s12/` = s(12) and s(11) only**, in the s21/s32 template: result, context (with jlevy/Kleddamag),
   certificate, verification (merge Lean into its own §), limits, reproduce.  Move §07 out.
3. **Give s(13) its own page, `/s13/`** (from `s12/docs/s13.html`, one more `cp` in `pages.yml`).  It is
   now the cleanest result on the site — hypothesis-free in Lean, 3,621 points, two checkers — and
   the method page for s(21)/s(32) (both link to `../s12/#s13` for "the method").  Put a cover plot
   there (the 09-22 review suggested one).  Keep a redirect-free anchor: leave a one-line stub at
   `/s12/#s13` linking to `/s13/`.
4. **s(45):** either a `/s45/` page on the s(21) template, or (cheaper) a section on `/s21/` "the same
   method at k = 7", plus an index card.  README currently is its only documentation.
5. Nav: 9 items already overflow on phones; with s(13)/s(45) group results under one "Results" entry
   (or a results strip on proofs.html) rather than adding more top-level links.

## 6. Edits made in this pass (uncommitted)

Lean status (job A):
- `s12/docs/s32.html`: tag "Lean 4 top theorem" → "kernel-checked in Lean 4"; §01 last ¶; §05
  rewritten (`s32_checkerCover`, `s32_eq_6`, `ZMTree.sound`, trusted base, gen_data/build_parts
  commands); §06 first caveat; §07 artifact rows for `S32Lower.lean`, `ZMTree.lean`, `LADDER.md`.
- `s12/docs/index.html` (`/s12/`): lede (drops "machine-checked from end to end", which overstated the
  3.9686 Lean status; adds s(13) Lean); §01 note that 35/9, 3920/997, s(11) are kernel-checked; §02
  s(13) row; §05 intro + Lean card (what is and isn't in the kernel); §07 Lean ¶ (ZMTree, `s13_eq_4`,
  trusted base); §08 Lean reproduce commands + four artifact rows.  Two `LEAN-3.968616` markers.
- `s12/docs/s21.html`: §05 new ¶ (s(13)/s(32) hypothesis now proved; segments are future work).
- `site/www/proofs.html`: playground aside (35/9 is a hypothesis-free Lean theorem); "Why twelve"
  marker; "Thirteen again" ¶4 (kernel check).
- `site/www/index.html`: s(32) and s(12)/s(13) cards; marker.
- `site/www/sources.html` §5: s(32), s(13), s(11)/s(12) small bounds Lean status; s(21) wording made
  precise ("from the checkers' covering statement"); marker.
- `site/data/lower_bounds.json` and `site/www/data/lower_bounds.json` (kept identical): s(32) notes
  "Lean top theorem from the checker hypothesis" → "kernel-checked in Lean with no hypothesis".

Small factual fixes (job B):
- `proofs.html`: "both [checkers] reject 23 deliberately broken certificates" → only the Rust one
  does (`tests/rung2/rejection_tests.sh` runs `zmcheck` only; `/s12/` §07 already said so).
- `s32.html` §06 and `s21.html` §06: "Nothing here settles … s(45)" was stale (s(45) = 7 since
  2026-09-27); now links the s(45) certificate, noting no Lean.

## 7. Open questions for Evan

1. Adopt Kleddamag's s(11) > 3.875 (and jlevy 3.96 for s(12)) after a replay?  How to phrase our
   s(11) result then ("an independent, Lean-checked 3.8143")?
2. s(17): which value goes in the Atlas (Mira 4.6141535 vs Guzhou0806/Kleddamag 4.6200)?
3. `/s13/` page and s(45) page/section: yes/no?
4. When the S12H kernel check lands: grep `LEAN-3.968616` (5 places) and upgrade the wording;
   `/s12/` lede can then say "checked in full inside Lean's kernel".
