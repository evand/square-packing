# TODO — the explainer site

Living list.  Items marked ✓ are done.  Findings from an audit on 2026-09-09; items marked
**(reproduced)** have been seen happening in a browser, the rest are still only the audit's reading
of the code and should be re-checked before they are fixed.  **(ship)** marks what should not go
public unfixed — see the Publishing section for why each one is on that list.

## Done 2026-09-22 (site review, `private/site/notes/site-review-2026-09-22.md`)

* Nav gains **s(12)**, on every page including `/s12/`, which now wears the Atlas topbar; the nav wraps
  instead of scrolling (on a phone "Sources" was cut off).  Explore's stage no longer assumes a 49 px
  topbar.
* Overview: lede says the site has original results; a "New: s(12) and s(13)" card and a Compare
  card; the year slider no longer widens the page on phones; three-column cards.
* Cross-links: Explore → Compare (when there are alternatives), Bounds (n ≤ 100), `/s12/` (n = 12),
  the Bentz walkthrough (n = 13); Bounds table n → Explore; Bounds footnote → Sources §1; ids on
  every Proofs and Sources heading.
* Writing: Bounds lede; Proofs "Why twelve", Bentz intro and the s(13) section tightened (the
  technical account now lives on `/s12/#s13`); Sources lede no longer says nothing here is
  original, §5 lists s(11) and the case-free s(13) as "this project", §3's asymptotics folded
  under `<details>`, §7 shortened; Explore gap note.  `/s12/` 29 % shorter, retitled "Twelve and
  Thirteen Squares", s(11) added, §06 rewritten against `s12/notes/n12-gap.md`.
* Bugs fixed: `catFill(-1)` black squares; playground calls before its data loads, page-wide arrow
  keys (now only when the diagram has focus), stale-load race, legend stuck on "7", "check all
  angles" gives feedback; Compare's play timer after a failed load; the doubled year on the Bounds
  card; the misplaced "show the count everywhere" checkbox.

## Done 2026-09-09

* ✓ **Bounds no longer calls s(5) and s(10) open.**  `lower_bounds.json` stored every floor rounded
  to 7 dp, which put the n = 5 floor at 2.7071068 — *above* the true 2 + 1/√2 = 2.707106781… — so it
  compared as greater than the matching upper bound.  Both rows read "open" with a gap of "-0.0000".
  Every entry already carried an `exact_form`, so nothing was lost: `value` is now the full
  double-precision evaluation of it (256 fields; largest correction 4.9e-8, invisible at the 6 dp the
  table shows).  Checked before writing that this flips exactly two verdicts and that no still-open
  case sits within 1e-6 of settled, so the 1e-9 comparison stays safe.  s(5) and s(10) now read
  "settled" with gap "—", no negative gap survives anywhere in the table, and the chart's summary
  went from "0 of the 36 shown are settled" to "2 of the 36".  17 entries have no parseable
  `exact_form` (n = 13, 17, 21 and the four flagged UNRELIABLE) and keep their curated decimal.
  The audit's knock-on claim did **not** hold: `overview.js` adds any n whose prose says "Proved"
  to the proved set, and both n = 5 and n = 10 do, so their tiles were never affected.
* ✓ **The tooltip no longer survives a packing change.**  It was cleared only on `pointerleave` and
  when the pointer left a square, so stepping n, picking a variant or clicking a gap row while the
  pointer rested over the stage left a tooltip describing a square from the packing you had just
  left — "square 5363 … centre (28.3891, 63.7444)" over a packing 4.68 wide.  Cleared at the top of
  `load()`, which covers the no-analysis path too: that one returns before `render()`.
* ✓ **Annotations no longer grow with the zoom.**  Gap markers and square numbers were sized in
  world units (`lw = s/700`, `0.28`), so zooming to a gap scaled them with the drawing — a two-line
  gap label stood taller than the stage.  `annScale()` pins them to the size they had when the
  packing was fitted and holds it from there in; below fit they stay tied to the squares, which is
  what reads correctly when a square is a few pixels across.  The gap circle keeps its true radius
  (`gap * 2`) — only the visibility floor counter-scales, so zoomed in the circle finally shows the
  gap at size.  Verified constant on screen through 8.26x.
* ✓ **Gap labels route around the square numbers.**  The numbers are pinned to square centres and
  cannot yield, so the gap labels take the detour.  Straight above the circle is still the default,
  but vertical-only nudging is not enough: a gap inside a column of squares has their numbers both
  above and below it, so nothing on the vertical is ever free and the search fell back to a position
  a full square away with its leader drawn down through a number (s(17), 6–16, through the 6).  The
  candidates are now eight directions at four radii, and one is refused if the label box *or* the
  leader segment would touch anything already placed — `segBox` is a Liang–Barsky segment/box test.
  If nothing at all is clear the label stays by its own circle rather than being flung across the
  drawing.  Only the neighbourhood can collide, so the 9465-number case is filtered to it first.
  Verified over 18 packings at fit and at four zoom levels: no label overlaps, no leader crossing a
  number or another label.
* ✓ **Every overlay counter-scales, not just the labels.**  The free/wedged outlines, the slide
  ghost and its direction line, the contact marks, the symmetry axes and the selection outline were
  all still sized in world units, so at gap-level zoom the dashes were slabs.  They hold their
  fitted size now too.  Stroke and marker sizes carry a `[constant, multiple of lw]` spec on the
  node (`elS`), so a zoom retunes them in place: the free-region layer paints and serialises a
  canvas per square — 4613 of them on the n = 9465 entry — and must not be rebuilt per frame.  The
  cheap layers (gaps, symmetry, selection) are redrawn instead, because their text offsets move
  with the scale.  ~11 ms per wheel event at n = 17/55/100 with every layer switched on.
* ✓ **Gap labels carry a halo** in the background colour (`paint-order: stroke fill`), so a square
  edge passing behind the digits no longer reads as a strikethrough.
* ✓ Square numbers are built once at unit size and scaled through the group transform, so a zoom
  retunes `transform` rather than rebuilding the layer: 29 ms rather than 280 ms at n = 9465.  The
  redraw is gated on `annScale()`, which is clamped at 1, so zooming out past fit and panning cost
  nothing.


* ✓ **The sidebar gap list was silently empty.**  `explore.html` names the sidebar list
  `<div class="gaplist" id="gaps">` and `viewer.js` names the SVG marker layer `<g id="gaps">`; the
  `<g>` is created first in document order, so `$('gaps')` resolved to it everywhere.  `renderSide`
  appended the rows as HTML `<div>`s *inside the SVG drawing*, where nothing paints them, and the
  sidebar div was never touched — no throw, so no console error to find.  The same collision let
  `renderSide`'s `G.innerHTML = ''` wipe the marker circles `drawGaps` had just drawn.  The SVG
  layer is `gapmarks` now; the page has no duplicate ids left.

* ✓ **Explore: click-to-select did nothing.**  `setPointerCapture` on `pointerdown` retargets the
  `pointerup` to the `<svg>`, so `closest('.sq')` always missed.  The square is now recorded on the
  down event; the selection is also outlined in the drawing.
* ✓ **Squares that slide along a zero-width channel were classed "moves with neighbours".**  The
  61x61 sampling grid can only see a 2-D region, and such a square's feasible translations are a
  *segment* (s(17) square 12, which slides 0.0712 into square 5).  `free_region` now line-searches
  every edge direction of the square and its neighbours as well; `region.kind` is `region` / `slide`
  / `turn`, and Explore draws a ghost outline at the far end of a slide.
* ✓ **Gap threshold raised** from `NEAR_MAX = 0.05` to `0.25`, and Explore states it rather than
  cutting the list off silently.  (s(17) had four gaps hidden just above the old ceiling.)
* ✓ `tools/patch_regions.py` — refresh gaps + regions + freedom on the exported JSONs without
  redoing the ~11 h rigidity search.
* ✓ **`?n=2` / `?n=3` no longer blank.**  Ellsworth's "2, 3" box links `square-2.svg` but embeds
  `square-3.svg`, which does not exist; `build_index` now falls back to the link and takes
  `pictured_n` from the file that is actually drawn.  `gotoN` says so instead of failing silently
  when there is no data at all, and Explore notes when the picture holds a different count than the
  n you asked for.
* ✓ **Optimiser start configurations dropped from the variant list** (`files[].start` in the index;
  a `?p=` link to one still loads).
* ✓ **Rigidity cross-check surfaced.**  Ellsworth's rigid page shows 18 packings he calls rigid
  beside 12 he calls not rigid or only semi-rigid; our independent analysis agrees on all 30.  The
  per-packing verdict is on Explore, the totals on Sources — both computed from the shipped data, so
  they move if the answer does.
* ✓ **Seven packings never loaded at all.**  `square-26_r1%27.svg` and six others are named with a
  literal `%27`; dropped into a fetch URL unescaped the server decodes it back to an apostrophe and
  the request 404s.  Names are encoded now, in `viewer.js` and `compare.js`.
* ✓ **A sidebar that cannot be built now says why.**  The "no exported analysis" path used to return
  early and leave the previous packing's panel on screen, and a throw inside `renderSide` left it
  half-built with no clue; both report into a red box at the top of the panel.  Optional annotation
  slots write through `setText`, so a cached older `explore.html` no longer takes the panel down.
* ✓ **Gap markers were drawn at the closest vertex**, which for two squares facing each other put the
  circle at the end of the facing edge and read as measuring the wrong pair.  They are centred on
  the overlap of the facing edges now, labelled with the pair as well as the distance, and nudged
  apart when several gaps pass through the same hole.
* ✓ **Slides truncated at the scan radius.**  The grid doubled its window when the region ran off the
  edge; the line search did not, so 302 slides were reported as exactly `R`.  Both grow now, and the
  exported direction and length carry enough digits to be replayed (a zero-width channel has no
  sideways clearance, so 6 significant figures could not be).
* ✓ **Claims trimmed to what the site does.**  "Slivers" dropped (jargon, and nothing measured it);
  Bounds says "n up to 100" — the range `lower_bounds.json` actually covers — and its two controls
  are capped there rather than at 150 / 324.
* ✓ **(ship) `compare.html?n=3` no longer goes blank** — and ten more n with it.  The filter dropped
  every file without `n_parsed`, so n = 3 (whose only drawing is `square-2.svg`, because Ellsworth's
  "2, 3" box embeds a `square-3.svg` that does not exist) emptied the list and `setup()` returned,
  leaving "0 packings on record" over an empty stage.  Compare now falls back to what the catalogue
  actually pictures, the way Explore already did, and says which n the drawing holds.  Three cases
  had been conflated and now read differently: 11 n borrow another n's picture (3, 147, 194, 223,
  232, 254, 264, 287, 290, 295, 322); 8 have a recorded s but no drawing at all (585, 852, 973,
  1169, 1310, 1536, 1697, 1953 — Ellsworth's own prose, quoted in the panel, says "SVG not made
  yet"); and anything else snaps to the nearest n *and says so*.  "0 on record" was itself untrue —
  the record exists, the drawing does not — so the count distinguishes them.
  Ten of those eleven were worse than blank: they are absent from `timeline.json`, so
  `loadN(TL[q.n] ? +q.n : 17)` sent `?n=147` to n = 17 and rewrote the URL to match, showing a
  different packing with no indication.  Navigation now runs over `timeline ∪ records`.
* ✓ **(ship) The 135 mislinked overview tiles no longer lie.**  Tiles for n with no catalogue entry
  linked `explore.html?n=N`, where `gotoN` silently snapped to the nearest n that had one — so
  clicking 12 showed 11 as though it were 12.  Those 135 are no longer links (`.nolink`), and the
  hover says "not in the catalogue — nothing to open"; the other 189 are unchanged.  Verified in the
  browser: 189 links, 135 not, and **zero** remaining links resolve to an n that would snap.  The
  same silent snap was still reachable by hand-typed URL, so `gotoN` now records what was asked for
  and Explore reports it, distinguishing that from the separate "one entry covers two n" case it
  already handled.  Checked that the claim it makes is true: n = 12 has no entry because s(12) = 4
  really is the grid, not because the entry is missing.
* ✓ **(ship) Compare can no longer hang the tab.**  Measured on the shipped data rather than
  estimated: the 8-symmetry Hungarian match costs 0.1 s at n = 201, 0.7 s at n = 1037, 7.7 s at
  n = 2135 and **29 s with a 603 MB cost matrix at n = 9465** — and the n box accepts anything in
  `timeline.json`.  Above `MATCH_MAX = 1200` the match is skipped, the morph controls are disabled
  rather than left to slide the walls around squares that cannot follow, and the panel says why.
  n = 9465 now loads in ~40 ms with a 6.5 ms event-loop turnaround; n = 1037 keeps its morph.
  Two things found while measuring and *not* done, deliberately: the SVG thumbnails looked like a
  second hang, but that was `requestAnimationFrame` being throttled in a background tab — measured
  properly they are 249 ms at n = 9465, so the canvas rewrite they seemed to need is unnecessary.
* ✓ **`loadN` was re-entrant.**  It awaits a fetch per packing, so two calls in flight let the older
  finish last and publish `A` from its own list, against which the newer render had already run —
  `A.rec` unset, `setup()` throwing on `A.rec.s`.  Found by stress, not by reading.  Newest call
  wins now; packings whose analysis fails to load are dropped *before* the strip is built, since the
  thumbnails carry their index into that list.  195 overlapping loads plus 120 arrow-mashes: clean.
* ✓ **A tie no longer reads as a change.**  n = 9465's two packings have identical s, which rendered
  as "+0.000000" — a change of nothing rather than no change.  It says which it is, and only claims
  an exact tie when the exported values are exactly equal.


## Done 2026-09-12 (merged into this repo 2026-09-22 from the untracked `site/` copy)

* ✓ **"Why twelve is still open" no longer says a fractional packing of mass 12 exists at side 4** —
  the proved figure is 12.27 (exact certificate, `s12/search/COVER4.md`); the sentence now says so.

* ✓ **`proofs.html` has a section on the case-free proof of s(13) = 4.**  Added after the Bentz
  walkthrough, as "Thirteen again, with no cases": the single weighted set of 3,621 points in
  `[0,4]²` of total weight 12.955972 < 13 that every closed unit square catches weight 1 from, the
  rescaling that makes disjointness (and so the count) work, why the *closed* convention is what
  puts the target at 13 rather than 16, the exact 12.2688 floor on any cover of that box — which
  also says why no such cover can ever speak about n = 12 — and how it is checked at margin zero:
  two exhaustive exact programs sharing no code (Python rationals, Rust integers), 23 rejection
  tests, soundness lemmas in Lean but not the subdivision code.  Says explicitly that the theorem
  is Bentz's (2010) and only the proof is new.  Facts and numbers from `s12/notes/s13-casefree.md`
  and `s12/README.md` "s(13) = 4 without case analysis"; links to the repository in the same form
  the playground aside already uses.  Hand-written HTML, so nothing to rebuild — `build.sh` only
  regenerates `www/data/`, never the pages.

## Done 2026-09-10 (pre-publication review)

* ✓ **Explore said "nobody has beaten the trivial packing" for 48 n that have a better one**
  (331 … 9465).  `gotoN` consulted only `index.records` — the main table, every n ≤ 324 plus five
  larger — so an n drawn only on a sub-page snapped to the nearest record and got the trivial-packing
  sentence.  It now opens the best of that n's own files first.  The sentence survives only for
  n ≤ `TABLE_MAX` = 324, where the data bears it out (no n ≤ 324 without an entry has any file below
  ⌈√n⌉); above that the page says the catalogue covers only selected n.
* ✓ **Explore lent the record's closed form and prose to alternatives** (39 and 94 files):
  `square-147_r2` showed s = 12.678779… over 7 + 4√2 = 12.657.  `s_tex`, `polys` and `prose` fall
  back to the record's only when the file *is* the record.
* ✓ Explore, smaller: junk `n` (`abc`, `0`, `-5`, `17.5`, an empty box) no longer prints "n = NaN";
  overlapping loads keep the newest (the race Compare already guarded); a failed load updates `S.n`;
  exact sides lose the ellipsis ("s = 2"); links to the seven `%27` files are no longer
  double-encoded (Ellsworth's server returned 404).
* ✓ **Bounds called s(5) open by 0.29289** in the one-n chart while the table said settled: the
  chart dropped undated packings, and the n = 5 record is undated.  Undated records now sit at the
  left edge, labelled "date not recorded".
* ✓ Overview no longer offers `square-51_r2__invalid.svg` (overlapping squares) as the 2023–25 record
  for tile 51: optimiser starts and files with non-`warn:` errors are skipped.
* ✓ **Data and attribution.**  n = 11's floor is back to Stromquist 2003 (the entry it replaced cited
  work that is not public yet).  `meta.caveats` no longer carries a local path.  The n = 82–85 note no
  longer calls Friedman's closed form a typo.  The finder tally — hand-written, not generated — had
  dropped Peter Shiu (Ellsworth writes "Kearney`<br>`and Peter Shiu") and Robert Wainwright ("Found
  *first* by").  Kearney–Shiu are credited with s(6) only; Bidwell's n = 17 is based on Hämäläinen's
  n = 17; Mira's 4.613 is attributed to the later weighted certificate; "All four" → five; Burns 6 Aug.
* ✓ Proofs names the 3.9686 floor as ours and unrefereed, credits Burns and Massaccesi, and no longer
  asserts a fractional packing at side 4 (that rests on unpublished work).  The 35/9 statement and
  the closed-convention aside now say the conclusion covers every smaller box.
* ✓ Links point at the new repository; the fetch User-Agent no longer carries a personal address.

## Bugs

1. `proofs.html`: the theorem box is written for the 81-point set; after switching to `cert_56` it
   still says 35/9 and 7 (the colour legend now follows the chosen set).
2. `overview.js:92` counts tilt angles as `n_angles − 1`, wrong for packings with no unrotated square
   (`50_r0r`, `200_r0r`); not seen at any slider year tried.
3. Explore's "fit" throws on a file with no analysis (`S.rec` is null).
4. Compare at t = 0 outlines most squares red, which reads as an error rather than "moves".

Refuted in the review, so not a bug: `hover()` does not trail the pointer while panning —
`setPointerCapture` retargets `pointermove` to the `<svg>`, and the tooltip stays hidden.

## Payload

* Bounds and Sources fetch the 589 KB `index.json` and use almost none of it (`IDX` is never read on
  Bounds; Sources needs only the 84-byte `checks`).  Overview reads ~40 % of it: `polys`, `s_tex` and
  `pages` could move to a separate file.
* Compare at n = 9465 fetches two 5 MB packings eagerly, though the match is skipped at that size.

## Data that is computed but never shown

* `lower_bounds.json` carries `sources[].url/title/venue`, `history[].url/date/note` and per-entry
  `note`s; Bounds cites "Mira 2026" with no link and never shows any `note` — e.g. n = 17's
  "external peer review still pending", Fort's not vouching for correctness, Burns's ChatGPT
  provenance.
* Bounds shows the n = 17 floor (4.613, Mira) above the n = 18 floor (4.445, Green).  Not
  propagating it by monotonicity was deliberate (`notes/lower-bounds-notes.md`), but the page never
  says so.
* Unused elsewhere: `analysis.n_corner_corner`, `contacts[].normal`, `contacts[].overlap`,
  `angle_groups[].count`, `region.area`; `index.summary.{waste,wedged,mobile,n_rotated,angles,min_gap}`;
  `timeline` entry `.prose/.rigid/.s_tex`.

## Promises the site does not keep

* No glossary.  free / slides / moves-with-neighbours / wedged / jammed, and the 1e-8 and 1e-20
  thresholds, exist only as tooltip fragments.  Which packings are drawn with truncated constants is
  never listed, though the "numeric contacts" tag hints at it.
* `proofs.html` never mentions the 2026 n = 17 certificates that `sources.html` §5 lists.

## Accessibility and small screens

* No `:focus`/`:focus-visible`, `aria-*`, `role`, `tabindex`, `<main>` or skip link anywhere.
* Keyboard-unreachable click targets: the viewer's legend rows, gap rows and freedom rows, and
  compare's thumbnails, are all `<div onclick>`.  Compare needs shift-click to set B — impossible on
  touch.
* Explore zoom is wheel-only with `touch-action:none`: on a phone you can pan but never zoom, and
  hover tooltips fire on tap.
* `site.css` hard-codes the topbar at 49 px (`height:calc(100% - 49px)`) but `.topbar` wraps, so the
  stage overflows between ~720 and ~900 px.
* Bounds charts use a fixed 1000x360 viewBox: at 360 px wide the 11 px axis labels render at ~4 px,
  and tooltips are `mouseenter`-only.
* `overview.js` never clamps its tooltip to the viewport, so it clips at the right and bottom edges.

## Later — grounding the epsilons

Three stages, cheapest first.  The motivating observation is what bug 2 actually looked like: the
comparison tolerance (`isSettled` at 1e-9) was *tighter than the precision of the data it compared*
(floors stored to 7 dp = 1e-7).  An epsilon is only meaningful between two quantities — the error in
the inputs, and the smallest difference worth resolving — and the code recorded neither, so there was
no way to notice the invariant was violated.  Making both explicit is the whole of the work.

### Stage 1 — stop discarding precision that already exists (cheap, no analysis)

`data/packings.json` and `www/data/p/*.json` keep `s` as a 28-digit **string**
(`"2.7071067811865475244008443621"`), and `viewer.js` then does `S.rec.s = +S.rec.s`.  Worse,
`timeline.json` — which is what Bounds actually compares — stores `s` as a JSON *number*, so it is a
double before it ever reaches the browser.  The square coordinates are likewise already doubles
(17 digits) in the export though the parser worked in higher precision.  So the first step is not
numerical analysis at all: keep the strings, and compare with a decimal type.  For the settled/open
question that alone is enough, because the values in question are exactly representable as decimals
or as short surds.

### Stage 2 — provenance-derived intervals instead of a fixed epsilon

Give every number an error bound derived from where it came from, and make the comparisons interval
predicates:

* scraped `s_dec`: 14 dp, and truncated rather than rounded in both cases checkable against a known
  exact value (n = 5 stores `2.70710678118654`, true value `…6547…`) — so `[d, d + 1e-14]`, or
  `d ± 1e-14` if Ellsworth's convention is not assumed.
* `exact_form` evaluated in floating point: IEEE 754 gives correctly-rounded `+ - * / sqrt`, so a
  k-operation expression is within about k/2 ulp; widen by k ulp and it is rigorous.
* certificate bounds (`15680/3951`, and Mira's rational) are exact — no interval needed.

All derived, none chosen, and all ~10^5 tighter than the 1e-9 in use today.

Worth being clear about what this buys: **intervals can prove "open" but never "settled"**.  Disjoint
intervals give a genuine, quotable proven gap.  Overlapping intervals only ever justify
"indistinguishable at the precision available" — which is a third verdict the table does not have and
probably should, instead of today's forced binary.

### Stage 3 — exact comparison, only where the verdict needs it

"Settled" is a claim about two reals being equal, so it needs symbolic evidence.  But the cheap
version is not an algebraic-number library: both sides already carry a closed form — `exact_form` on
the bound, `s_tex` and the minimal `polys` on the packing — and for n = 5 both are `2+1/sqrt(2)`.
Comparing normalised closed forms settles it without any algebra.  Only genuinely different-looking
forms, or the degree-18 roots, would need real algebraic-number arithmetic, and those can fall back
to stage 2's honest "indistinguishable".

### Not every epsilon is a precision epsilon

There are ~29 in the JS and ~51 in `tools/*.py`, and at a glance they all look like `1e-9`.  They are
three different things and only the first should be derived:

* **precision claims** — `isSettled`, the `rec.s >= t - 1e-9` trivial test, `overview.js`'s proved
  test.  These are the ones stage 2 is for.
* **modelling tolerances** — "tilt below 1e-9 means the author meant this square unrotated", the
  angle-merge threshold.  These encode intent, not measurement error; tightening one because it
  resembles a precision epsilon would be a regression.  The merge slider is the pattern worth
  copying: make the tolerance explicit and let the reader move it.
* **algorithmic guards** — divide-by-zero protection and the like, including the `|ux| > 1e-9` in the
  gap-label placement.  Leave them alone.

The first job is therefore to label the existing epsilons by kind; the count of genuine precision
claims is small.

## Publishing

The site lives in `site/` of github.com/evand/square-packing, next to the published s(12) proof in
`s12/`.  Pages deploys `site/www/` as the site root and `s12/docs/` as `/s12/`
(`.github/workflows/pages.yml`).  Active research stays in a private repository until it is ready.

Before announcing it anywhere, tell David Ellsworth: his prose for all 194 records ships verbatim in
`index.json`, and every packing is derived from his SVGs.  Offer to change or remove anything.
