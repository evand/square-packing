# Site review, 2026-09-22

Scope: the Square Packing Atlas (`site/www`, and its older published copy `public/site/www`) plus the
s(12) write-up (`s12/docs/index.html`, served at `/s12/`).  Read-only review; nothing edited.
How it was checked: source read in full.  All seven pages were served locally and screenshotted at
1400 px with headless Chromium.  Every local `href`/`src`/`#anchor` was checked by script, and every
JS `fetch` target was listed.  Horizontal overflow was measured at a 360 px viewport (pages loaded in
375 px iframes).  Headless Chromium will not render narrower than 500 px, so the "mobile" screenshots
are not trustworthy; the iframe measurement is.

---

## 0. Merge inventory: `site/` vs `public/site/`

Evan's decision: `public/` is the one repo, and `site/` folds into `public/site/`.  The two trees
have **diverged in both directions**.  `public/` carries the publication-day fixes of 2026-09-10:
the rename to "Square Packing Atlas", links to the new repo, and a batch of viewer and data fixes.
`site/` carries the 2026-09-12 content: the s(13) section and the s(11) floor.  Neither side can
simply be copied over the other.  There are no files present on one side only (apart from
`__pycache__`).

| file | newer side | what each side has | merge action |
|---|---|---|---|
| `www/index.html` | public | public: title/brand "Square Packing Atlas" | take public |
| `www/explore.html` | public | brand only | take public |
| `www/compare.html` | public | brand only | take public |
| `www/bounds.html` | public | brand only | take public |
| `www/js/viewer.js` | public | public: `TABLE_MAX`, stale-load guard (`S.gen`), junk-`n` guard, sub-page variants, no ellipsis on exact `s`, the record's prose/closed form not lent to alternatives, unencoded `%27` source links | take public |
| `www/js/bounds.js` | public | public: undated records kept (fixes n = 5 history) | take public |
| `www/js/overview.js` | public | public: skips optimiser starts and invalid files | take public |
| `www/proofs.html` | **both** | public: brand, fixed repo link, "fit in no box smaller than this" wording, rescaling explanation in the aside, "Why twelve" rewritten (credits Burns/Massaccesi, links `s12/`, drops the "ceiling" claim). site: the **"Thirteen again, with no cases"** section (lines 85–95) and the 12.27 sentence | start from public, paste site's s(13) section, then fix its links (§3); reconcile the two "Why twelve is still open" paragraphs by hand (public's wording, site's 12.27 figure if the negative record goes public) |
| `www/sources.html` | public | public: brand, Kearney–Shiu date, Bidwell/Hämäläinen wording, Wainwright + joint-credit row, "five items", Mira's later bound reworded, repo link, **"Code: …" link** (site still says "to be published") | take public; nothing site-only |
| `www/data/lower_bounds.json` | **both** | public: repo URLs → `square-packing/tree/main/s12`, reworded n = 82–85 typo note. site: **n = 11 floor `3040/797 = 3.814304`** (source `EvanDaniel2026`, URL to the archived repo) | take public; re-add site's n = 11 entry **only once `s12/certificates/s11_lower_3.8143.*` is in the public repo** (it isn't yet — public deliberately kept Stromquist), with the URL pointed at the new repo |
| `data/lower_bounds.json` | both | same as above (the build copies it to `www/data/`) | same |
| `notes/lower-bounds-notes.md` | both | public: relative data path; drops the "could not verify s(11) ≥ 3040/797" section. site: keeps that section, now stale in the other direction (the certificate exists) | take public; add one line recording the n = 11 decision |
| `build.sh`, `data/ellsworth/fetch.sh` | public | public: fetch list / user-agent pointing at the repo, not an email address | take public |
| `README.md` | public | public: renamed, explains that the fetched data is not committed and how Pages deploys | take public |
| `TODO.md` | site | site: the "Done 2026-09-12" block | take site, then retitle as needed |

**The s(12) write-up has the same problem.**  `s12/docs/index.html` (2026-09-13, private) and
`public/s12/docs/index.html` (2026-09-10) both changed since their common base, s12 commit
`1eb1363`.
- **Public has corrections the private file lacks.**  It says the uniform certificates use weight 1/7
  each; the private file still says "every point weight exactly 1/5", which is wrong for the
  1736-point certificate.  It also has the corrected resolution footnote: N = 6000/12000 for the main
  certificate, and the verifier refuses coarser nets.  And it adds a link back to the Atlas.
- **Private has content public lacks:** §07 (the case-free s(13) proof), and a §06 rewritten around
  the 12.0282-at-3.99 ceiling and `COVER4`.
- **The merge is not clean.**  A trial `git merge-file` (base `1eb1363`) gives **4 conflict hunks**
  in `docs/index.html` and 2 in `README.md`; `VERIFICATION.md` merges clean.
- **Links break on the move.**  The private §06–§08 link to `search/*.md` and `notes/*.md` on the
  archived `square-packing-12` repo, `blob/master`; public is `tree/main/s12`.  Whether those links
  resolve depends on piece 1: which notes and `search/` files are pushed.

---

## 1. Site map

| page | words | purpose | links out (beyond the top nav) | reached from |
|---|---|---|---|---|
| `index.html` Overview | 200 | hook + 324-tile grid of all n, coloured by packing kind / year / … | tiles → `explore.html?n=`; cards → Explore, Bounds, Proofs, Sources | nav; `/s12/` "Atlas" link (public only) |
| `explore.html` | 72 + panel | one packing: drawing, angles, freedom, gaps, shape | Ellsworth's SVG + page | overview tiles, Compare "explore →" |
| `compare.html` | 102 | all packings for one n, A/B morph | `explore.html?p=` | nav only |
| `bounds.html` | 221 | gap chart n ≤ 100, one n's history, table | none | nav, overview card |
| `proofs.html` | 1,493 | unavoidable points, 81-point playground, Bentz walkthrough, (site only) case-free s(13) | repo link (archived repo in `site/`; `s12/` in public) | nav, overview card |
| `sources.html` | 1,789 | bibliography, packers, finder tally, 2026 preprints, explainers | ~90 external | nav, overview card, Bounds footnote (text only, no link) |
| `s12/` write-up | 3,115 | the s(12) ≥ 3.9686 result (and, in private, s(13)) | GitHub; public adds "Atlas →" | Proofs (public only). **Not in the nav, not on the Overview** |

**Link check.**  No dead local links and no bad anchors in either tree; every JS `fetch` target
exists.  External links point to the archived repo in `site/`: 2 in proofs, 1 in sources, 4 in the
data, and 10 in the private `s12/docs`.  These are fixed in public, but will return if site's files
win the merge.

**The reader's likely path, and where it breaks:**
1. Overview grid → click 17 → Explore.  Works well.  Explore then leads nowhere:
   - there is no "compare the 6 packings of 17" link, although the variant dropdown knows they
     exist;
   - there is no "bounds for n = 17" link, although the floor is the obvious next question;
   - it doesn't link to Proofs even for proved n.
2. Bounds: the table rows and chart points don't link to Explore for that n.  The chart is
   click-to-history only.
3. Proofs → s(12): in `site/`, nothing links to `/s12/` at all.  In public there is one inline link.
   **The site's own original results (s(12) ≥ 3.9686, the case-free s(13)) are invisible from the
   Overview**, whose lede still reads as a pure catalogue ("This site collects the record packings…").
4. The `/s12/` page is a different design (no topbar, its own fonts and layout), so arriving there
   feels like leaving the site.  Public added one link back; the private version has none.
5. Compare is only in the nav.  The Overview cards list Explore/Bounds/Proofs/Sources but not
   Compare.

**Proposed map (minimal):**
- **Nav:** Overview · Explore · Compare · Bounds · Proofs · **s(12)** · Sources.  Either add `/s12/`
  to the nav, or make it a sub-page of Proofs with the Atlas topbar.  Give `/s12/` the Atlas topbar
  (a copy of the `<div class="topbar">`; the fonts already match).
- **Overview:** one sentence and one card for original work: "New here in 2026: a machine-checked
  s(12) ≥ 3.9686 and a proof of s(13) = 4 with no case analysis →".  Add a Compare card, or drop the
  cards entirely: the nav is right there and the cards repeat it.
- **Cross-links:**
  - Explore side panel: "all N packings of n → Compare", "floor for n → Bounds?n=", and "how it's
    proved → Proofs" where one exists;
  - Bounds table: n → Explore;
  - Bounds n = 12/13 rows → `/s12/` or `proofs.html#s13`;
  - give the Proofs `<h2>`s `id`s (none exist now) so they can be linked.
- **Sources:** leave the page alone; see §3 for trims.  Keep it last.

**Small screens (360 px, measured).**
- Real page overflow on two pages:
  - Overview: the fixed `width:260px` year slider pushes the page to 402 px;
  - `/s12/`: the context table is 397 px.
- The nav scrolls sideways inside the topbar, so on a phone "Sources" is cut off with no hint.
  With a 7th item this gets worse, so wrap it.
- Everything else stays within 360 px (Bounds table and Compare strip scroll in their own boxes).
  The site TODO's accessibility list stands; none of it was fixed in either tree.

---

## 2. Visualizer / explainer inventory

| piece | page | shows | state | earns its place? |
|---|---|---|---|---|
| **Overview grid** (324 tiles, 5 colourings, year slider) | index | every n at a glance; the slider replays 1979→2026 | works; clean and striking | yes, and a strong **announcement image** (the "kind of packing" colouring at "today") |
| **Packing viewer** | explore | the drawing with free-square slide arrows, contacts, gaps, symmetry, numbers; side panel of rotations, freedom, gaps, shape | works; n = 17 looks great. TODO bugs 1 and 3 (`catFill(-1)`, hover while panning) are still in both trees. The tag chip "NONE" (symmetry) is cryptic out of context | yes; the core of the Atlas. **Feature n = 17** with free squares on |
| **Compare morph** | compare | all packings of n on a strip; A/B slider morph | works. On load at t = 0, 9 of 17 squares are outlined red, which reads as an error; shift-click for B is impossible on touch (known) | yes, but secondary; it needs an entry point from Explore |
| **Gap chart** (n ≤ 100) | bounds | floor vs best, per n | works. The heading says "up to 100" but the default is 60; the `[√n]` baseline label is cryptic | yes |
| **History chart + story cards** | bounds | one n over time | works. **Bug:** the floor card prints the year twice, "(Mira 2026, 2026, unrefereed)", because `bounds.js:114` appends `ll2.y` to a `src` that already contains it. At n = 17 the 2026 floor dots stack on the label "floor 4.61303" | yes |
| **Bounds table** | bounds | all n ≤ 100 | works; no links out | yes; add n → Explore |
| **81-point playground** (drag/rotate square, live count, per-angle heatmap, "check all angles", 56-point alternative) | proofs | the unavoidable-set idea, touchable | works. The "show the count everywhere" checkbox sits centred on its own line above its label (layout glitch). TODO bug 2 is open: "check all angles" blocks the main thread, and arrow keys are captured page-wide | **yes: the best explainer on the site**, and the one to feature for s(12) |
| **Bentz walkthrough** (clickable leaf tree + figure) | proofs | Bentz's six-leaf s(13) proof, each leaf's count of 12 | works, and good.  But the intro paragraph above it and the figure caption for "The 16 starting points" say the same thing (see §3) | yes |
| **s(13) case-free section** | proofs (site only) | text only, **no figure** | not yet published | a plot of the 3,621 weighted points (dot size = weight) would earn it a place next to the Bentz figure. Optional |
| **81-point certificate figure** (three sample squares with capture counts) | `/s12/` §04 | static SVG | fine, but it duplicates what the playground does interactively | keep one sentence and link to the playground, or keep it as the non-JS fallback |
| **Terminal transcripts** (`verify`, `zmcheck` output) | `/s12/` §05, §07 | proof of checking | fine for this audience | yes, but one is enough on the page |
| **Rigidity check line** | sources §7 | our rigidity vs Ellsworth's "Rigid" marks | works | move it to Explore's Freedom section, or to a footnote; on Sources it's orphaned |

**For an announcement:** (1) the Overview grid, (2) Explore n = 17 with free squares showing,
(3) the 81-point playground, (4) the Bentz leaf tree.  All are real screenshots already.

---

## 3. Writing review

**Tone overall.**
- Consistent and good: plain, precise, first-person-plural, honest about status.
- The recurring weakness is **stacked qualification**: sentences that state a fact, then its
  provenance, then its caveat, then why the caveat matters, all in one run.
- The two long pages (Sources 1,789 words, `/s12/` 3,115) could lose 25–35 % without losing a fact.
- Proofs is about right apart from the s(13) section, which is the longest thing on the page and
  explains closed semantics at length.
- The `/s12/` page switches to first-person singular ("if a better published bound … I did not find
  it"), while the Atlas says "our".  Pick one.

### Duplicated explanations (feeds the "proof deduplication and consolidation" TODO)

The same four explanations appear in 3–4 places each:

| explanation | proofs.html | `/s12/` | `s12/README.md` | `notes/s13-casefree.md` |
|---|---|---|---|---|
| weighted unavoidable set → counting | §"The trick" + "Why twelve" | §01 | §The method | §1 |
| closed semantics: why 13 not 16 | s(13) ¶3 (≈150 words) | §07 "What a closed cover is" (≈220) | yes | yes |
| margin zero → why sampling can't work, adaptive boxes | s(13) ¶4 | §07 "Margin zero" | yes | yes |
| two checkers, 23 rejection tests, Lean covers the lemmas and not the subdivision | s(13) ¶4 | §07 table + ¶ | yes | yes |

**Suggested division of labour.**  proofs.html gives the *idea* at popular level: one paragraph
each, no checker details, linking to `/s12/`.  `/s12/` is the *technical* write-up, with the checker
table and transcripts.  The README is a *repo guide*: two lines per result, then how to run
`verify.sh`, with the rest linked to `/s12/`.  The notes are the record, linked from `/s12/` only.
Concretely:
- cut proofs' s(13) ¶3–¶4 to about 90 words (rewrite below);
- cut the README's result prose to the theorem statements.

### Per page: the passages to tighten

**index.html** (200 words)
1. Lede.  It lists what the site does but not that it contains new results.
   > *Rewrite:* "Call the answer s(n). For a perfect square it's obvious; for most other n nobody knows — only the best packing found and a proven floor beneath it. This atlas draws every record packing from Friedman's and Ellsworth's catalogue, lets you take them apart, shows how the exact results are proved — and includes two new computer-checked results of our own, for n = 12 and 13."
2. The foot line, "As of today: 182 of 324 … dates are read from its wording and rounded to the
   year."  Fine, but the attribution repeats the lede.  Drop the second sentence.
3. The hint "hover a tile · click one to explore, where the catalogue has an entry · dot = proved
   optimal".  The dot legend belongs in the legend row, not the hint.

**bounds.html** (221 words)
1. The lede's first two sentences are in the wrong order and repeat "known".
   > *Rewrite:* "For each n two numbers are known: the side of the best packing anyone has found, and a floor that a proof says no packing can beat. Where they meet, s(n) is settled; elsewhere the truth lies in the gap — and for most n the gap has barely moved in forty years. Floors are collected for n ≤ 100."
2. "The gap, for n up to 100" while the default shows 60: say "The gap" and let the control speak.
3. The footnote under the table is fine.  Link "Erich Friedman's survey" to Sources §1.

**explore.html / compare.html** (72 / 102 words)
- Compare's intro is one long instruction sentence.
  > *Rewrite:* "Every packing known for this n, oldest first. Click one for the left side (blue), shift-click for the right (gold), then drag the slider to morph. Outlined squares change place."
- The side note "Matching is by position and angle; when two packings are 'alternatives' …" →
  "The morph matches nearest squares; it is a visual aid, not a physical motion."
- Explore: the "Smallest gaps" note "Closest approach between two squares that do not touch, or a
  square and a wall, in units of one square's side. nothing above 0.25 is measured." has a lowercase
  sentence start.
  > *Rewrite:* "Closest non-touching pairs (square–square or square–wall), in side lengths; only gaps below 0.25 are listed. Click to zoom."
- There is no glossary for free / slides / moves-with-neighbours / wedged.  It's a known TODO item;
  a 4-line `<details>` under "Freedom" would do.

**proofs.html** (1,493 words)
1. **The Bentz intro paragraph** (≈190 words) repeats the figure caption for "The 16 starting
   points" almost line for line (the points A–D, the mirroring, the slack of 3).
   > *Rewrite of the paragraph:* "Bentz (2010) proves s(13) = 4 from 16 unavoidable points in [0,4]². Thirteen boxes against sixteen points leave three points of slack — enough to force two boxes to sit alone with a corner point, where a replacement lemma pins them further. The argument then splits six ways. Click through the leaves:"

   and leave the coordinates to the caption.
2. **s(13) ¶3, closed semantics** (≈150 words).
   > *Rewrite:* "The closed convention — a point on a square's edge counts — is what makes 13 the target. In the 4 × 4 tiling a point on a shared edge is credited to both squares; shrink the squares a hair instead and the sixteen become disjoint, so no cover could weigh under 16. Nor can any cover weigh under 12.2688 (proved exactly), which is also why no cover of this box can ever speak about twelve squares."
3. **s(13) ¶4, checking** (≈150 words).
   > *Rewrite:* "Here the margin is zero — at the corner placement the only captured point lies on the square's edge — so nothing can be sampled. Two independent exact checkers (Python rationals, Rust integers) split the space of positions and angles into boxes until every box is settled; both finish with none left over, and both reject 23 deliberately broken certificates. Their lemmas are proved in Lean; the subdivision code is not."
4. **s(13) closing aside.**  It names "the square-packing-12 repository" and shows `notes/s13-casefree.md`
   and `search/COVER4.md` as unlinked monospace paths.  Point it at `/s12/#…` and the new repo, and
   drop the paths.  It also says "only the proof is new", then "it is a computation, not something a
   reader checks by hand"; the second clause repeats ¶4.
5. The playground aside ("The map samples 120 × 120 positions …") is good; public's version of it is
   better.  The last sentence, "The skeleton is the same in every such set found…", is repeated on
   `/s12/` §04; keep it here only.

**sources.html** (1,789 words)
1. The intro, "Nothing here is original except the interactive parts and the contact/rigidity
   analysis", is **now false**: the s(12) bound and the case-free s(13) proof are original.
   > *Rewrite:* "The packings are Ellsworth's and the classical proofs are their authors'. Original here: the interactive tools, the contact and rigidity analysis, and the 2026 computer-checked bounds in §5 marked 'this project'. Preprints are unrefereed; this page says which is which."
2. §2's closing paragraph (Hämäläinen, Wainwright, Cottingham, Stenlund, Gustafsson) is ~90 words.
   Turn it into a 5-item list like the rest of the page, or cut it to one line: "Early records by
   Hämäläinen, Wainwright, Cottingham, Stenlund and Gustafsson reached Gardner by mail; attributions
   follow Friedman."
3. §3's asymptotic-waste entries (Erdős–Graham through McClenagan, 6 items, ~200 words) aren't about
   anything else on the site.  Put them under a `<details>` "Asymptotics (the large-n question)".
4. §7 mixes provenance, precision policy and a live JS check.
   > *Rewrite:* "Built in 2026 from Ellsworth's 30-digit constants. Gaps and free squares are computed in double precision; 'exact contact' flags are checked at 50 digits. Code: github.com/evand/square-packing."

   Move the rigidity-check line to Explore (§2).
5. §5, our own entry: "evand, 2026, square-packing-12" becomes the public wording, and add the
   s(13) re-proof and the s(11) floor (once published).  Label it as ours: "this project".

**`/s12/`** (3,115 words)
1. **Lede.**  "…closing about 85% of the remaining gap, and an honest account of how far it does not
   go."  Keep it.  The four tag chips repeat the lede ("Lower bound improved", "Main conjecture still
   open", …); drop all but the GitHub link.
2. **§03 "No packing below 4 was found"** (≈130 words) is a search report on a proof page.  Cut to
   two sentences:
   > "A basin-hopping search over all 36 coordinates, validated on n = 5, 10, 11, finds nothing below side 4 — every run collapses to the compressed 4 × 4 grid. Schadt's and Ellsworth's 2025–26 annealing campaigns, which improved many records, never improved 12."
3. **§04's second and third paragraphs** (the history of 3.92 → 3.9318 → 3.9486 → 3.9686, and "The
   points were found by linear programming…").  Useful to specialists, ~250 words.  Compress to one
   paragraph of ~80 words, or move to the README.  The private version still has the "1/5" error
   here (§0).
4. **§06 "What this does not do"** is the longest section and now partly superseded by
   `notes/n12-gap.md`: every degree-1 family ≥ 12 at t = 4, rank-8, obstruction X.  Rewrite it once
   piece 1 decides what's public: three sentences of conclusion plus a link to the negative record.
   Until then, public's hedged version is the safer text.  The private version's "It is not 4" and
   "real research problem, not a longer run" are good lines; keep them.
5. **§07** duplicates proofs.html (table above).  Since `/s12/` is the technical home, keep the
   checker table and one transcript here, and trim proofs.  Also: the page is titled "Packing Twelve
   Squares" but now carries a thirteen result.  Retitle it ("Twelve and thirteen squares"), or split
   §07 into `/s13/`.

**Jargon used before it's defined:**
- "free squares" (Overview colour option, Explore);
- "wedged", "moves with neighbours" (Explore);
- "rigid" (Sources §7);
- "boxes" (proofs, Bentz: it's defined, but after "boxes" appears in the case list);
- "closed semantics" (`/s12/` §06, which is before §07 defines it);
- "degree-1", "ν_f" (`/s12/` §06: ν_f is defined, degree-1 isn't);
- "fractional packing of mass 12.27" (proofs "Why twelve": undefined for a popular reader.  Public's
  rewrite dropped it; keep it dropped, or add "a spread-out 'fractional' arrangement").

**Factual staleness:**

| item | where | fix |
|---|---|---|
| "Code: to be published" | `site/` sources §7 | public already fixed |
| `square-packing-12` links | `site/` proofs ×2, sources ×1, data ×4; private `s12/docs` ×10, `s12/README.md` ×2 | → `github.com/evand/square-packing/tree/main/s12` (and `blob/main/s12/...` for files) |
| s(11) ≥ 3.8143 | only in `site/` data; not on `/s12/`'s context table (row 11 shows just Trump's packing), not in sources §5 | show it everywhere once the certificate is public |
| s(13) case-free | nowhere in public; the site's Bounds n = 13 row and Sources §5 don't mention it | add "re-proved case-free, 2026 (unrefereed)" to the n = 13 row note and to §5 |
| "Nothing here is original" | sources intro | rewrite (above) |
| weight "1/5" | private `s12/docs` §04 | take public's 1/7 wording |
| verification-resolution footnote | private `s12/docs` §05 | take public's |
| brand "Squares in Squares" | every `site/` page | that is Ellsworth's page title; public's "Square Packing Atlas" is right |
| "Mira 2026, 2026" | bounds story card | `bounds.js:114`: don't append the year when `src` already ends with it |

---

## 4. Top-10 changes, in priority order

1. **(M) Do the merge in §0 first.**  Take public for 12 of the 16 differing files.  Merge
   proofs.html and lower_bounds.json by hand, and three-way merge `s12/docs/index.html` with base
   `1eb1363`.  Everything below assumes the merged tree.
2. **(S) Make the original work findable.**  Put s(12) in the nav (or a sub-item of Proofs), add an
   Overview lede sentence plus a card, and give `/s12/` the Atlas topbar.
3. **(S) Fix the staleness table:** repo links, the "Nothing here is original" intro, s(13) in
   Bounds and Sources §5, and s(11) once the certificate is pushed.
4. **(S) Cross-links for flow:** Explore → Compare / Bounds / Proofs for the current n; Bounds
   rows → Explore; `id`s on Proofs `<h2>`s.
5. **(M) Deduplicate the proof explanations** by the division of labour in §3: proofs.html gives
   the idea, `/s12/` the technical write-up, the README is a repo guide.  About 400 words come out
   of proofs, and the README shrinks.
6. **(M) Trim `/s12/`:** §03 to two sentences, §04 history to one paragraph, §06 rewritten against
   the negative record (after piece 1), the tag chips dropped.  Retitle it or split out s(13).
7. **(S) Trim Sources:** the intro, §2 closing paragraph, §3 asymptotics under `<details>`, §7.
   Move the rigidity check to Explore.
8. **(S) Small bugs seen in this pass:**
   - the doubled year on the Bounds card;
   - the misplaced playground checkbox;
   - Overview slider `width:260px` → `max-width`;
   - `/s12/` table overflow → `overflow-x:auto` wrapper;
   - wrap the nav instead of scrolling it.
9. **(M) The site TODO's open bugs 1–3** (viewer `catFill(-1)`, playground main-thread
   "check all angles" and page-wide arrow keys, hover while panning).  Before announcing, fix the
   playground ones at least: it is the feature most likely to be tried from a tweet.
10. **(L, optional) A figure for the case-free s(13):** the 3,621 weighted points, dot area ∝ weight,
    next to the Bentz figure.  It turns the most novel result into something visible.  Then the
    accessibility list from the site TODO, if a public push wants it.
