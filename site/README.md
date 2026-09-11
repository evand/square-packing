# Squares in Squares — an explainer site

A static site about the problem *how small a square can hold n unit squares?* (`s(n)`): the record
packings, what can be said about them (angles, contacts, free squares, gaps, symmetry), how the
records changed over time, the proven lower bounds, and how the few exact results are proved.

Everything geometric is derived from the exact constants in David Ellsworth's SVG catalogue
(https://kingbird.myphotos.cc/packing/squares_in_squares.html, continuing Erich Friedman's survey).
Attributions and dates are quoted from his wording. See `www/sources.html`.

## Layout

| path | what |
|---|---|
| `www/` | the site: `index.html` (overview grid), `explore.html` (one packing, analysed), `compare.html` (history + morph), `bounds.html`, `proofs.html`, `sources.html`; static JSON under `www/data/` |
| `tools/parse_svg.py` | independent parser for Ellsworth's SVG dialect → exact squares (mpmath, 50 digits) |
| `tools/analysis.py` | contacts (exact-checked), first-order rigidity LP + nonlinear verification, sampled free regions, symmetry, angle groups |
| `tools/parse_all.py`, `tools/export.py`, `tools/timeline.py`, `tools/scrape_pages.py` | the pipeline, see `build.sh` |
| `data/` | `lower_bounds.json` (curated, cited).  The fetched pages/SVGs and intermediate JSON are not committed; `build.sh` regenerates them |
| `notes/` | research notes: prior work, lower-bound sourcing |

## Build

```sh
./build.sh all        # fetch (polite, ~5 min), scrape, parse (~20 min), export (analysis), timeline, bounds
python3 -m http.server -d www 8765
```

`build.sh export` is incremental (skips packings already exported); `tools/export.py --force` redoes all.

A fresh clone has only the exported `www/data/`, which is enough to serve the site; rebuilding it
needs `./build.sh fetch` first.  GitHub Pages deploys from `.github/workflows/pages.yml` at the
repository root: `site/www/` becomes the site root and `s12/docs/` becomes `/s12/`.

## Conventions

* Frames are y-up with the container at `[0, s]²`; angles are `θ mod 90` in degrees.
* Contacts closer than 1e-8 are contacts; those with a 50-digit gap above 1e-20 are flagged as
  numeric (the packing was not analytically optimised).
* Gaps (`analysis.near_misses`) are reported up to `NEAR_MAX = 0.25` — a quarter of a square's side.
  Explore lists the 12 smallest.
* "Free" = the square can move on its own; "moves with neighbours" = a verified finite collective
  motion; "wedged" = the linearised contacts allow a motion but the nonlinear projection found none;
  everything else is jammed (rigorous: the LP with corner–corner contacts dropped is a relaxation).
* A free square's `region` records *how* it is free, as `kind`:
  - `region` — a genuine 2-D patch of translations, drawn from the sampled mask;
  - `slide` — the feasible translations are a **segment**: the square is wedged between parallel
    faces and slides along a channel of exactly zero width. Such a set almost never meets a node of
    the sampling grid, so it is found by line searches along every edge direction of the square and
    its neighbours (`_slide_dirs` / `_slide_scan`), not by the grid. `region.slide` gives the unit
    direction and the length; the site draws the far end as a ghost outline.
  - `turn` — it cannot translate at all but can rotate on the spot.
  `dx` / `dy` are exact line-search reaches, not grid-quantised bounds.

## Checks

* Ellsworth's rigid page lists 18 packings he calls rigid and 12 he calls not rigid or only
  semi-rigid; our independent first-order analysis agrees on all 30 (`index.json` → `checks`, shown
  on Sources and per packing on Explore).
* Parsed square counts match the filename for all 530 SVGs; two files that are not valid packings
  (an optimiser start configuration and a file named `__invalid`) are flagged and skipped.
* Contacts: `exact` is decided at 50 digits from the SVG's own constants; packings drawn with
  truncated constants (e.g. `square-41_r4`, ~8 digits) show as "numeric contacts".
