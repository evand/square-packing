# Literature

This is the index of our literature searches.  Put new scans and reading lists here
(`literature/<topic>-<date>.md`).  The older scans below stay where they are: other files link to them, and dated notes are
not rewritten.  When you add a scan, add a row here.

Labels used in the scans: **[verified]** / [P] means we read the source or its bibliographic record; [S] means secondhand;
**[unverified]** / [I] means from memory or our own inference.

## Squares-in-squares: history, exact values, conjectures

| File | Date | What |
|---|---|---|
| [`site/notes/research-prior-work.md`](../site/notes/research-prior-work.md) | 09-08 | Broad prior-work survey: primary sources, explainers and tools, a timeline of lower bounds and exact results, rigidity, 2025–26 news. |
| [`notes/literature-s32.md`](../notes/literature-s32.md) | 09-26 (upd. 10-01) | Status of s(32), the k²−4 family and exact values before our results. |
| [`tasks/k2m3-lit/REPORT.md`](../tasks/k2m3-lit/REPORT.md) | 09-30 | Prior status of s(k²−3) = k: Friedman's conjecture, Kearney–Shiu, Bentz 2010/2016, Nagamochi; suggested attribution sentence. |
| [`tasks/friedman-lit/REPORT.md`](../tasks/friedman-lit/REPORT.md) | 09-30 | Roth–Vaughan and the Friedman survey (Roth–Vaughan itself not reached). |
| [`notes/wishlist-literature-2026-10-03.md`](../notes/wishlist-literature-2026-10-03.md) | 10-03 | What is known for each WISHLIST item. |
| [`notes/conjectures-literature-2026-10-04.md`](../notes/conjectures-literature-2026-10-04.md) | 10-04 | Double-check for the conjectures page: Chung–Graham prizes, Erdős #106, Ellsworth's {s(n²+1)} ≥ ½. |
| [`waste-bounds-2026-10-12.md`](waste-bounds-2026-10-12.md) | 10-12 | Primary-source read of Erdős–Graham 1975, Roth–Vaughan 1978 (exact theorem, c₀ = 10⁻¹⁰⁰), Stromquist's 1984 memoranda I–III, Bui 2504.09489, Bašić–Slivková 2018 (s(61) ≥ 7.8906, superseded); claim-by-claim check of our files. |

## Packing search: records, optimizers, methods

| File | Date | What |
|---|---|---|
| [`search/packer/lit-provenance.md`](../search/packer/lit-provenance.md) | 10-04 | Who found each best-known packing for n ≈ 50–135 and how (from Ellsworth's site and SVG comments). |
| [`search/packer/lit-optimizers.md`](../search/packer/lit-optimizers.md) | 10-04 | Numerical optimizers: record finders, circle packers, phi-functions, nesting, Torquato–Jiao ASC, Berthold et al. 2026. |
| [`search/packer/lit-statmech.md`](../search/packer/lit-statmech.md) | 10-04 | Stat-mech / global optimization and basin census, aimed at the k-row jamming failure. |
| [`search/packer/lit-2026-10-12.md`](../search/packer/lit-2026-10-12.md) | 10-12 | Pallet-loading DP, symmetry-constrained starts, grain boundaries, LNS, cluster tricks, continuation; ranked top 8. |
| [`group-chat-2026-10-12.md`](group-chat-2026-10-12.md) | 10-12 | Reading list from the group chat (irregular/polygon packing, phi-functions, MINLP); read in the next row. |
| [`irregular-packing-2026-10-12.md`](irregular-packing-2026-10-12.md) | 10-12 | The group-chat list read against our packer: separating-line / phi-function NLPs, sentinels, MIQCP, soft polygons.  Verdict: nothing beyond `fq`; one speculative idea (soft-square aspect homotopy for k-row jams). |
| [`soft-squares-2026-10-12.md`](soft-squares-2026-10-12.md) | 10-12 | Full read of Melashenko 2025 and Bennell 2024 (soft rectangles / stretched polygons: no hardening, no rigid baseline, no equal squares); design for log-aspect soft squares in `fq`; why in-place relief absorbs k-row jams; pre-registrable crossing test at 54 / 71 / 110 with a rigid-schedule control. |
| [`gensane-ryckelynck-2026-10-12.md`](gensane-ryckelynck-2026-10-12.md) | 10-12 | Full read of Gensane–Ryckelynck 2005 (DCG 34), the method that Schadt's SA reimplemented: hard-constraint maximin inflation ω, a feasible-only random walk, and deflate/re-inflate MBH with 1/1 step adaptation.  They report the same one-angle jam and fix it only with restarts.  Status of 11/17/18/29/37 (none a record now; their s(37) was already beaten by Cantrell's 2002 packing); comparison table; 4 ranked ideas (none targets the k-row jam); citing papers. |

## Method: computer-assisted proofs

| File | Date | What |
|---|---|---|
| [`computer-assisted-proofs-2026-10-12.md`](computer-assisted-proofs-2026-10-12.md) | 10-12 | Closest prior work on method: interval B&B for circles in a square (Markót–Csendes n = 28–30, Markót n = 31–33), Montanher et al. (3 squares in a circle), Fekete et al. (squares into a disk), Flyspeck, Hales–Kusner pentagons, Tammes 13/14, the 2026 s(11) proof (Lean with `native_decide`), the hand-proof lineage of unavoidable sets; comparison table, positioning, techniques worth borrowing (Flyspeck's sharp-corner lemma, face reduction, minimal-precision kernel replay), paywalled DOIs. |

Private (in `../private/`, not linked): `s12/notes/lit-statarb.md`, lessons from Evan's statarb sampler for `mcmin.py`.
