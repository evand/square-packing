# k2m3-site: Bentz page + site catch-up (2026-09-30)

Repo `~/math/square-packing/public`.  Site root `site/www/`; per-result write-ups in `s12/docs/*.html`, served by
`.github/workflows/pages.yml` (copies; add the new page and its cover file there).  Read `s12/docs/REVIEW-2026-09-28.md`
(site review memo: page roles, the 7-section template, stale facts) and the existing `s12/docs/s45.html`, `s21.html`,
`s32.html` as the template.  Also `site/TODO.md`.

1. **New page `s12/docs/k2m3.html`** served at `/k2m3/` (+ `cover.txt` = the box cover): s(k²−3) = k for all k ≥ 6
   (Bentz's conjecture, with the small cases).  7 sections as the template: result, context, certificate/idea (the
   fixed-profile family: corner modules + periodic wall band + Lebesgue interior, saving 4D per box for every k; one box
   k = 7 certifies all k ≥ 6; zero margin along walls/germs; a figure of the family if cheap — SVG from the family file),
   verification, Lean, caveats, reproduce.  Content source: `s12/search/QUADRANT_EXACT.md`, `s12/search/QUADRANT.md`,
   `s12/notes/lean-bentz-reduction.md`, the bundle `s12/certificates/k2m3/` (being written in parallel by another agent:
   link to `certificates/k2m3/README.md` and `verify.sh`; don't edit it), reviews in
   `~/math/square-packing/private/s12/tasks/k2m3-review/*/REPORT.md` (private: summarise, don't link).
   **Status line, verbatim, prominent:** "Working in public: a single-implementation exact certificate, adversarially
   reviewed by six independent agents with no errors found; the all-k reduction is kernel-checked in Lean. Not yet
   independently re-implemented, externally reviewed, or fully formalised."  Tone: promising result, not an announcement.
   Small-case attributions: mark with `<!-- LIT -->` (a literature check is running); current belief k = 3
   Kearney–Shiu 2002; k = 4, 7 Bentz 2010; k = 5, 6 Bentz arXiv:1606.03746 preprint.
2. **Catch-up** (stale per `s12/TODO.md` "Site: s(60)/s(61) not yet on index cards, Results subnav, proofs.html"): add
   s(60) = 8 and s(61) = 8 (from `s12/certificates/s60/README.md`) and the k²−3 family to index cards, the Results subnav
   on every page, `proofs.html` (one paragraph each), `bounds.html`/`sources.html` as appropriate, and
   `site/data/lower_bounds.json` + `site/www/data/lower_bounds.json` (check how they relate; an infinite family: add the
   concrete n = k²−3 values the data model can hold, e.g. up to the table's n range, with a note; ask nothing, decide
   and report).  s(60) needs no page of its own unless the pattern demands one (s45 has one: follow the pattern,
   a short s60 page is fine).
3. Build/check: `site/build.sh` if it exists; validate HTML (no broken internal links; check with a small script), dark
   mode and phone width per the existing CSS.  Also update `README.md` (repo root) and `s12/README.md` results lists.

**Rules.**  Work in the main checkout on main; touch only `site/`, `s12/docs/`, `.github/workflows/pages.yml`,
`README.md`, `s12/README.md`; don't commit; don't touch `s12/certificates/`, `s12/search/`, `s12/lean/`.  Cores 8–11.
Report ≤ 400 words: files changed, decisions made, anything left.
