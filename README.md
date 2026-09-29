# square-packing

How small a square can hold *n* unit squares?  That side, `s(n)`, is known exactly for only a
handful of *n*.  This repository has two parts:

| | |
|---|---|
| [`site/`](site/) | **Square Packing Atlas**, an explorer for the record packings in David Ellsworth's catalogue: every packing drawn and analysed (angles, contacts, free squares, gaps, symmetry, rigidity), how the records changed over time, the proven lower bounds, and how the exact results are proved.  **https://evand.github.io/square-packing/** |
| [`s12/certificates/s60/`](s12/certificates/s60/) | **s(60) = 8**: a mixed cover of `[0,8]²` (weighted points plus mass spread uniformly along the grid lines) of total `59.8587 < 60`, certified at margin zero by the same two independently written exact checkers as s(21) = 5 and s(45) = 7, hence also **s(61) = 8** (s is monotone); not in Lean yet.  No separate write-up: the certificate directory's README is the documentation (summarised on the s(45) page). |
| [`s12/certificates/s45/`](s12/certificates/s45/) | **s(45) = 7**: a mixed cover of `[0,7]²` (weighted points plus mass spread uniformly along the grid lines) of total `44.7735 < 45`, certified at margin zero by the same two independently written exact checkers as s(21) = 5; not in Lean yet.  Write-up: **https://evand.github.io/square-packing/s45/** |
| [`s12/certificates/s21/`](s12/certificates/s21/) | **s(21) = 5**: a mixed cover of `[0,5]²` (weighted points plus mass spread uniformly along the grid lines) of total `20.8947 < 21`, certified at margin zero by two independently written exact checkers, with a Lean 4 top theorem from that computational hypothesis (the checkers' covering statement; segments are not yet in the kernel verifier).  Write-up: **https://evand.github.io/square-packing/s21/** |
| [`s12/certificates/s32/`](s12/certificates/s32/) | **s(32) = 6**: a weighted closed cover of `[0,6]²` of total weight `31.7135 < 32`, certified at margin zero by two independently written exact checkers, and kernel-checked in full in Lean 4 with no hypothesis (`s32_eq_6`).  Write-up: **https://evand.github.io/square-packing/s32/** |
| [`s12/`](s12/) | Machine-checked results on unit squares in a square: **s(12) ≥ 15680/3951 = 3.968616…** (still the best lower bound for s(12) we know of), **s(11) ≥ 3040/797 = 3.814304…** (since superseded by jlevy's 3.827 and Kleddamag's 3.875), and a **case-free proof of s(13) = 4** (one weighted closed cover, checked at margin zero by two checkers sharing no code, and kernel-checked in Lean with no hypothesis; chelokot's Lean archive kernel-checked Bentz's proof earlier).  Exact weighted certificates, a Rust verifier over the full continuum of placements, an independent Python re-check.  In Lean's kernel: s(11) ≥ 3040/797, s(12) ≥ 35/9 and ≥ 3920/997, s(13) = 4; for 3.968616 only the reduction so far.  Also the research log, including a detailed record of why these methods stop short of s(12) = 4.  Write-ups: **https://evand.github.io/square-packing/s12/** (s(11), s(12)) and **https://evand.github.io/square-packing/s13/** |

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
page (`site/www/sources.html`), including the 2026 work this builds on or runs beside: Burns, Massaccesi,
Fort, Mira, jlevy (Squares Project), Kleddamag, Guzhou0806, tokoharu, wand125 and chelokot.
