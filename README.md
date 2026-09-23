# square-packing

How small a square can hold *n* unit squares?  That side, `s(n)`, is known exactly for only a
handful of *n*.  This repository has two parts:

| | |
|---|---|
| [`site/`](site/) | **Square Packing Atlas**, an explorer for the record packings in David Ellsworth's catalogue: every packing drawn and analysed (angles, contacts, free squares, gaps, symmetry, rigidity), how the records changed over time, the proven lower bounds, and how the exact results are proved.  **https://evand.github.io/square-packing/** |
| [`s12/`](s12/) | Machine-checked results on unit squares in a square: **s(12) ≥ 15680/3951 = 3.968616…** and **s(11) ≥ 3040/797 = 3.814304…** (exact weighted certificates, a Rust verifier over the full continuum of placements, an independent Python re-check), and a **case-free proof of s(13) = 4** (one weighted closed cover, checked at margin zero by two checkers sharing no code).  Lean for the reductions and the checkers' soundness lemmas.  Also the research log, including a detailed record of why these methods stop short of s(12) = 4.  Write-up: **https://evand.github.io/square-packing/s12/** |

[![verify](https://github.com/evand/square-packing/actions/workflows/verify.yml/badge.svg)](https://github.com/evand/square-packing/actions/workflows/verify.yml)
re-checks the `s12/` certificates (fast tier) and builds the Lean whenever `s12/` changes; the slow
sweeps run on demand (`s12/verify.sh --full`).

## Status

Nothing here has been peer reviewed.  The results are computer-assisted and meant to be
checked: `s12/verify.sh` rebuilds the verifiers and re-checks every certificate.

## Data and licence

The code is MIT-licensed (`LICENSE`).  The packings themselves are David Ellsworth's
(https://kingbird.myphotos.cc/packing/squares_in_squares.html, continuing Erich Friedman's
survey).  `site/www/data/` is derived from his SVGs and quotes his attribution text; that
material is his and is not covered by the MIT licence.  Full credits are on the site's Sources
page (`site/www/sources.html`).
