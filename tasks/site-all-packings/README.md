# Site: draw every packing we have numbers for

Status 2026-10-04: not started.  Raised in the 10-04 site session.

## Goal

Draw (Explore, Compare, Overview, Bounds upper bounds) every packing whose coordinates are public,
not just Ellsworth's SVGs: including n > 324 and packings outside the catalogue.

## Policy (Evan, 2026-10-04)

- Pure coordinates of a packing are a mathematical object, not art; copyright doesn't really apply,
  and we act accordingly: we may parse and redraw them.
- Do our honest best to credit: easy links, the finder's name, verification status.
- Do not copy write-ups, papers or software without an explicit license that allows it
  (Couzo's repo has no LICENSE: take the numbers, not the code or README).
- Re-running searches ourselves to re-derive or improve packings is possible, but would muddy credit
  and make Compare/Bounds weird; not part of this task.

## What we know (from jlevy's register, live 2026-10-04)

`packing/atlas/known-best/bound-citations.json` and `packing/frontier/n-NNN.md` in jlevy/squares
(exact forms there, e.g. n = 68: 137481175581603/15625000000000).

- Couzo 2026 (github.com/franciscouzo/square-packing, f3c5a529, 2026-09-27): 49 packings, n = 68..307,
  each below the catalogue.  jlevy verified (exact + interval replay, T-056) all but 206, 259, 305
  (reported).  Gains 3e-8 (n = 172, 199) to 0.022 (n = 301).  In n <= 100 only n = 68 (9e-7).
- de Winter 2026: n = 211 at 14.99796 (< grid 15; no catalogue entry), verified (T-057).
- jlevy's "Friedman & Ellsworth" upper bounds for 147, 232, 264, 290, 295 are not new: shared
  catalogue entries (s(147) = s(148) entry).  The Overview mis-showed them; fixed in d9eb91b.

## Pieces

1. Data: an "outside the catalogue" packing source next to `site/data/ellsworth/` (coords, finder,
   link, date, verification status), parsed and run through `tools/analysis.py` like the SVGs.
2. Upper bounds from it on Bounds, Overview tiles, Explore's large-n notice (`upperFor`).
3. Explore/Compare draw them, with credit lines; Sources §1 sentence ("parsed from Ellsworth's SVG
   files") must change.
4. Ellsworth's own SVGs for n > 324 are already parsed (up to 9465); the Overview stops at 324.

## Local copies (2026-10-07)

Read-only clones in `~/math/_untrusted-third-party/` (no LICENSE in either: coordinates with credit only):
- `franciscouzo-square-packing/` @ 6042c56b (2026-10-03; supersedes f3c5a529 above, which changed
  n = 208, 209, 228, 263, 272, 303, 306): `nNNN.txt` = header `# s = ...` then `x y theta(rad)` per
  square, 49 n (68..307); README table has side values.
- `itsnaka-squish-certs/` @ e63e4e52 (2026-10-07; SQUISH, jlevy/squares#401): 23 n (88..303) below
  the catalogue/Couzo; `squish-submission-*/nNNN/nNNN.cert.json` = exact rationals (centre x, y,
  t = tan(θ/2), `s_exact`); README table says which folder is current per n.
- jlevy-squares @ be8172aa has rational witnesses for both (`packing/witnesses/squish-401-2026/`,
  `franciscouzo-2026-10-03/`, `kingbird-2026/`).
