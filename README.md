# square-packing

How small a square can hold *n* unit squares?  That side, `s(n)`, is known exactly for only a
handful of *n*.  This repository has two parts:

| | |
|---|---|
| [`site/`](site/) | An explorer for the record packings in David Ellsworth's catalogue: every packing drawn and analysed (angles, contacts, free squares, gaps, symmetry, rigidity), how the records changed over time, the proven lower bounds, and how the exact results are proved.  **https://evand.github.io/square-packing/** |
| [`s12/`](s12/) | A machine-checked lower bound **s(12) ≥ 15680/3951 = 3.968616…**: an exact certificate, a Rust verifier over the full continuum of placements, an independent Python re-check, and a Lean formalisation of the reduction.  Write-up: **https://evand.github.io/square-packing/s12/** |

[![verify](https://github.com/evand/square-packing/actions/workflows/verify.yml/badge.svg)](https://github.com/evand/square-packing/actions/workflows/verify.yml)
re-checks every s(12) certificate whenever `s12/` changes, and monthly.

## Status

Nothing here has been peer reviewed.  The s(12) bound is computer-assisted and meant to be
checked: `s12/verify.sh` rebuilds the verifier and re-checks every certificate.

## Data and licence

The code is MIT-licensed (`LICENSE`).  The packings themselves are David Ellsworth's
(https://kingbird.myphotos.cc/packing/squares_in_squares.html, continuing Erich Friedman's
survey).  `site/www/data/` is derived from his SVGs and quotes his attribution text; that
material is his and is not covered by the MIT licence.  Full credits are on the site's Sources
page (`site/www/sources.html`).
