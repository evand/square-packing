# Provenance of best-known squares-in-squares packings, n ≈ 50–135 (and the n²−n frontier)

Compiled 2026-10-04 from David Ellsworth's live site plus the SVG source comments of individual records.
Pages were downloaded as HTML with curl and read as text. No downloaded code was executed.

Sources (all fetched 2026-10-04):
- Main list: https://kingbird.myphotos.cc/packing/squares_in_squares.html
- Older/alternative packings ("compared"), split into three pages:
  https://kingbird.myphotos.cc/packing/squares_in_squares__compared.html (n ≤ 100),
  …__compared2.html (101–196), …__compared3.html (197+; not read)
- s(n²−n−1) pattern page: https://kingbird.myphotos.cc/packing/squares_in_squares__n^2-n-1.html
- Analytic minimization note: https://kingbird.myphotos.cc/packing/squares_in_squares__analytic_minimization.html
- Ellsworth's edit of Friedman's survey: https://kingbird.myphotos.cc/packing/squares.html
- Göbel strips/squares and rigid pages (skimmed; nothing relevant to 50–135 beyond what the main list says)
- SA run statistics: https://kingbird.myphotos.cc/packing/square-51.stats.txt, …/square-55.stats.txt
- Per-record SVG source comments, which hold the detailed history and dates: https://kingbird.myphotos.cc/packing/square-NN.svg
  for NN = 50, 51, 53, 55, 68, 69, 71, 83, 87, 88, 103, 105, 110, 126, 129, 131, 132, 152, 154, 155, 156
- Tools repo: https://github.com/Davidebyzero/packing_squares_in_squares__tools
- Secondary source: H. Gill essay, https://tomrocksmaths.com/wp-content/uploads/2026/08/square-packing-trm-essay-harry-gill.pdf

Abbreviations used below:
- "SA-rand": started from randomness.
- "mod-SA": Ellsworth's modified version of Schadt's GPU simulated annealing program. Versions "#2" and "#3" are named on the page.
- "opt": analytic optimization, meaning an exact polynomial-root solution by Ellsworth in Mathematica.

---

## 1. Records in 50–135 that are not trivial or purely classical

Trivial (best known s = ⌈√n⌉, nothing pictured): 56–61, 72–78, 90–97, 111–118, 133–135.
Proven trivial (Nagamochi 2005): 62, 63, 79, 80, 98, 99, 119, 120.

Purely classical records with no modern activity (Göbel 1979 / Stenlund 1980 / Friedman 1997 constructions):
- 52, 65 (Göbel)
- 66 (Stenlund)
- 67, 84, 104, 125 (extend Göbel's s(52))
- 82, 101, 122 (add "L"s to Göbel's s(65))
- 85 (Friedman; the 26/85 series)
- 86 (Friedman; extends Gustafsson's alternative s(18))
- 127 (extends s(86))
- 89, 109, 124 (continue Göbel's pattern)

The genealogy of every other record:

| n | s (live) | Genealogy and method (dates from SVG comments where given) |
|---|---|---|
| 50 | 53/7 = 7.571428… | Previous record was Cantrell Sep 2002 (adds an "L" to his s(37)). **Schadt, 13 Dec 2025, SA-rand** (found s=7.57143688). Ellsworth optimized it 14 Dec 2025. Rational side length from the {3,4,5} tilt; resembles Ellsworth's non-record s(104) of Dec 2024. |
| 51 | 7.700799… (deg 12) | Previous record was Hajba Jul 2009 (7.70435). **Schadt, ≤15 Jan 2026, SA-rand.** Ellsworth refound it with mod-SA #2 SA-rand on 31 Jan 2026 and optimized it 1 Feb 2026. Statistics: 4 hits in 3004 sub-8.0 finds, about 4.9 GPU-hours per hit (see §3). |
| 53 | 7.822876… = 13/2+√7/2 | Cantrell Sep 2002, improved by Cantrell Dec 2024. **Ellsworth, 7 Feb 2026, mod-SA #3 SA-rand**, optimized the same day. |
| 54 | 7.846667… (closed form) | Cantrell Oct 2005, improved by DeVincentis Apr 2014. No activity since. |
| 55 | 7.945771… | Previous record was the DeVincentis Apr 2014 s(n²−n−1) pattern, improved by Ellsworth/Cantrell in 2023–24 and by Ellsworth's mod-SA in Dec 2025 (7.9542). **Schadt, 3 Jan 2026, GPU SA starting from a "cherry-picked" state (s=7.9951) produced by his reimplementation of the Gensane–Ryckelynck (G&R) algorithm, run SA-rand.** Ellsworth refound it with mod-SA #3 SA-rand on 5 Feb 2026, refined it, and optimized it 8 Feb 2026 (about 40 GPU-minutes per hit). |
| 68 | 8.798796… (not yet opt) | Brendberg Jun 2023 (his own program plus manual polishing). **Schadt, 19 Dec 2025, SA-rand.** Ellsworth opt 26 Dec 2025. **Joost de Winter, 14 Aug 2026**, "unspecified AI" with an evolutionary beam search seeded from the existing record with its symmetry relaxed. **Allen Chang, 10 Sep 2026**, GPT‑6 Astra (+"TheMagicAnimals"). **Ruilin Wang, 10 Sep 2026**, ChatGPT‑6 Astra. **Jake Loyd, 19 Sep 2026**, his own L‑BFGS program, from randomness. Five improvements in five weeks: the most contested n in the range. |
| 69 | 8.827195… (deg 38) | Morandi Jun 2010, improved by Cantrell Aug 2023. **hmbelvedere.com, 29 Jul 2026**, "undisclosed iterative method (probably with AI)". Ellsworth refound it with mod-SA (starting from Morandi's s(69)) and optimized it 7 Sep 2026. |
| 70 | 8.881667… (deg 4) | DeVincentis Apr 2014. No activity since. |
| 71 | 8.944072… | Cantrell Oct 2005, matched by the DeVincentis 2014 n²−n−1 pattern. **Schadt, 29 Dec 2025, GPU SA from an s=9.00000194 approximation of Cantrell's s(71).** Ellsworth improved it with mod-SA from Cantrell's s(71) on 25 Jan 2026 and optimized it the next day. |
| 83 | 9.634758… (deg 672) | Stenlund 1980 (adds an "L" to s(66)). Hajba Sep 2024, Cantrell Nov 2024 (now extends Bidwell's s(17)). **Allen Chang, 4–10 Sep 2026, GPT‑5.6 Sol / GPT‑6 Astra**, who also found the polynomial-root form. Ellsworth optimized it and verified it to 10⁶ digits. |
| 87 | 9.838815… (deg 41) | Ellsworth Dec 2024 (from Hajba's s(107) and DeVincentis's s(54)), improved by Cantrell 8 Jan 2025. **Ellsworth, 2 Feb 2026, mod-SA starting from his own Dec 2024 s(87).** **Allen Chang, 10 Sep 2026, GPT‑6 Astra.** Ellsworth reconstructed the full polynomial 16 Sep 2026. |
| 88 | 9.888153… (deg 20) | See §2. Friedman 1997 (generalizes Cottingham's s(41)), improved by Cantrell Aug 2002. Ellsworth 25 Nov 2024 (adapts Cantrell's s(102) improvement) and 29 Nov 2024 (adapts Cantrell's 2002 s(37) technique). Then Cantrell and Ellsworth independently improved it in Jan 2025 (Ellsworth's SVG: 16 Jan 2025) with the "new technique". **"Improvement by Thomas Schadt pending."** |
| 102 | 10.611388… (deg 8) | Hajba Sep 2024 (extends Stenlund's s(37)). Cantrell and Ellsworth Nov 2024 (Cantrell's s(37) technique), Ellsworth Dec 2024, Cantrell Jan 2025 (the "new technique"). No SA activity. |
| 103 | 10.703782… (not yet opt) | Before Dec 2025 nothing beat s(104), a Göbel s(52) extension. Ellsworth's best hand attempt (Dec 2024–Jan 2025) was **invalid** (two overlapping pairs). **Schadt, 19 Dec 2025, SA-rand.** Ellsworth refound it with mod-SA SA-rand on 1 Jan 2026, plus manual nudging of squares fed back into the program. **Tej Stead, 16 Jun 2026, "working with Claude Fable 5".** |
| 105 | 10.807619… (not yet opt) | Previous record added an "L" to Friedman's s(86). **Schadt, 30 Dec 2025, SA-rand.** Ellsworth refined it the same day, then refound and refined it with mod-SA #3 on 4–6 Mar 2026. |
| 106 | 10.822980… (deg 32) | Ellsworth Nov 2024 (based on Cantrell's s(53)), improved by Cantrell Dec 2024. |
| 107 | 10.846667… (closed form) | Hajba Nov 2024: an alternative extension of DeVincentis's s(54). |
| 108 | 10.925919… (deg 144) | Hajba Oct 2024, improved by Ellsworth Nov 2024 (now extends Trump's s(11)). |
| 110 | 10.996793… (not yet opt) | No nontrivial packing was known before. **Cantrell, 16 Feb 2025**, with an exact analytic solution. **Ellsworth, 17 Jan 2026, GPU mod-SA.** **Tej Stead, 16 Jun 2026, Claude Fable 5.** |
| 123 | 11.601400… (deg 12) | Ellsworth Dec 2024 (extends Hajba's s(102) plus Cantrell's s(37) technique), improved by Cantrell Jan 2025 (the new technique). |
| 126 | 11.774735… (not yet opt) | Friedman (2 "L"s on s(86)). Ellsworth 14 Dec 2024 (from Cantrell's s(39)), improved by Cantrell 16 Dec 2024. **Schadt, 16 Dec 2025, SA from a simplified previous best.** Ellsworth re-improved it with mod-SA and optimized it 8 Jan 2026. **Joost de Winter, 14 Aug 2026, AI-assisted evolutionary beam search with SA, seeded from the Dec 2025 s(105).** |
| 128 | 11.825092… (deg 40) | Ellsworth Nov 2024 (from Morandi's s(69)), improved Dec 2024 with Cantrell's s(53) technique. Cantrell and Ellsworth improved it jointly in Jan 2025. |
| 129 | 11.881306… (deg 20) | Ellsworth Dec 2024 (Cantrell's s(37) technique). **Ellsworth, 31 Dec 2025, mod-SA.** Optimized 10 Jan 2026. The SA "introduces a new technique … of putting two squares of different angles against each other in a row". |
| 130 | 11.911191… (deg 8) | Ellsworth, Cantrell, and Ellsworth again, all Nov 2024. Extends Friedman's s(88) and Cantrell's s(37) technique. **"Further improvement pending."** |
| 131 | 11.956525… (not yet opt) | Hajba Nov 2024 (DeVincentis n²−n−1 pattern). Schadt and Ellsworth tried in Dec 2025 from a constructed state (not records). **Ellsworth, 23 Jan 2026, GPU mod-SA from Hajba's s(131)**, which "resulted in significant reorganization". **Tej Stead, 16 Jun 2026, Claude Fable 5.** |
| 132 | 11.991373… (not yet opt) | No nontrivial packing was known before. **Arslanov, Mustafin & Shangitbayev, Mar 2019** (E-JC, doi:10.37236/8586). **Cantrell, 21 Mar 2025**, exact. **Ellsworth, 17–18 Jan 2026, GPU mod-SA.** **Tej Stead, 16 Jun 2026, Claude Fable 5.** |

The same pattern holds just past the range:
- s(152): Ellsworth/Cantrell Jan 2025, then **Grigoriy Dyachkov, 10 Sep 2026, Claude Opus 5 + Fable 5**, applying the Jul 2026 s(69) trick.
- s(154): Ellsworth Dec 2024, Cantrell Feb 2025, Ellsworth mod-SA Jan 2026, Stead Jun 2026.
- s(155): mod-SA from s(182) with 27 squares removed, Jan 2026, then Stead.
- s(156): AMS 2019, Ellsworth Dec 2024, Cantrell Mar 2025, Ellsworth mod-SA Jan 2026, Stead Jun 2026.

---

## 2. The s(88) note

Quoted exactly from the main list (the same text appears on the compared page):

> **88** — $s = 9.88815305375857$ (degree-20 root)
> Found by David Ellsworth in November 2024, by adapting and extending the $s(37)$ improvement found by David W. Cantrell in September 2002.
> Improves upon the $s(88)$ found by Erich Friedman in 1997, which extended the $s(41)$ found by Charles F. Cottingham in 1979.
> Improved independently by both David W. Cantrell and David Ellsworth in January 2025.
> This new technique, which will need to be applied to about 13 additional packings previously thought to be finished, independently found by both David W. Cantrell and David Ellsworth.
> **Improvement by Thomas Schadt pending.**

The s(88) SVG source comment (square-88.svg) does not mention Schadt. Its last entry is "Improved by David Ellsworth on January 16, 2025."

What this implies:
- The current s(88) = 9.8881530537… is known not to be the best available. Thomas Schadt has a better packing, presumably from his SA program, that Ellsworth has not yet posted or optimized.
- I could not find the value, date, or method of the pending improvement on any page or SVG. Mark: **unverified / not public**.
- The "about 13 additional packings previously thought to be finished" remark appears at s(88), s(102), and s(123). It means Cantrell and Ellsworth know of a still-unapplied construction technique from Jan 2025 that should improve about 13 more records. The page does not list which ones.
- Two other records in or near the range carry "Further improvement pending": s(130), which extends s(88), and s(172)/s(199)/s(228).
- So the s(88) family (88, 130, 207, …) is actively being improved.

---

## 3. What the site says about the search software

### Schadt's program
None of the pages describes the moves, objective, or cooling schedule. Everything known comes from the attribution lines, SVG comments, and stats files:
- It is a simulated-annealing program Schadt wrote himself, and it runs on a GPU: SVG comments for 55, 71, 110, 131, 132, 155, and 156 say "GPU simulated annealing program".
- His first record appeared 13 Dec 2025 (s(50)). Through Dec 2025–Jan 2026 he set records at 39, 41, 50, 51, 55, 68, 71, 103, 105, and 126, and improved non-record points at 131, 155, 182, 240, 272, …
- The usual mode is "starting from randomness". Other starting points were used:
  - an approximation of a known record: s(71) from an s=9.00000194 approximation;
  - a simplified previous best: s(126);
  - a constructed n²−n−1 "mock-up": s(182);
  - a "cherry-picked state" from his earlier, separate program that "reimplements the Gensane & Ryckelynck algorithm": s(55).
- **Key difficulty statement**, from square-55.svg:
  > "The cherry-picked state has s=7.9951486563, which was important because **without special modifications, the simulated annealing algorithm almost always gets stuck just above the trivial size (8.0 in this case) due to stacked rows and/or columns.** When left to continue being optimized by the G&R reimplementation, it went in a different direction and got stuck at s=7.9577055242."
- The s(n²−n−1) page says DeVincentis's 2014 pattern "held the record for over a decade on s(29), s(41), and s(71) … but thanks to Thomas Schadt's simulated annealing program has been shown to be inoptimal for all n > 3". Simulated annealing also showed that the odd-n pattern "becomes more chaotic when optimized (although it retains symmetry)".
- Code availability: **not public as far as I can find**. It is not in Ellsworth's GitHub. A web search found only a secondary essay (Gill 2026), which gives no technical detail.

### Ellsworth's modified versions
- He refers to "his modified version", and also to "version #2" (Jan 2026, used for s(51)) and "version #3" (Feb–Mar 2026, used for s(53), s(55), and s(105)).
- s(51) notes say Schadt's later refinement was "possibly incorporating the refinement technique from David Ellsworth's version of the program". So Ellsworth's main change appears to be a better final-refinement stage. The changes are otherwise **undocumented**.
- Hardware and parameters, from square-51.stats.txt: "NVIDIA RTX 3080 Ti with the simulated annealing program set to use 65536 threads". Runs were tuned with "the median average temperature and cooldown period taken from the previous sessions' finds". So the knobs are a temperature and a cooldown period.
- **s(51) statistics** (9 sessions, 31 Jan–1 Feb 2026), counting runs that finished below 8.0:
  - 2816 refine to Hajba's 2009 s=7.70435;
  - 184 refine to s≤7.70375;
  - **4** refine to the record 7.70080 ("exceedingly rare").
  - Each find took about 23.6 s, so one record hit takes about **4.9 h** of GPU time.
- **s(55) statistics** (5+1 sessions, 4–5 Feb 2026):
  - 1893 runs got below 8.0, and 24% of those got below 7.96;
  - 816 refined to the old 7.95417 record;
  - 5 refined toward the new record 7.94577, some only after "more delicate refinement";
  - about **40 GPU-minutes per record hit**.
- Workflow notes from SVG comments:
  - s(103): "manually moving some squares in the lower-right area, and feeding this back into the program, to coax it into finding a particular improvement";
  - s(53): "first found the pattern using version #3, then started refining it; started optimizing it after the refinement passed the mark of beating the previous record";
  - s(131): significant reorganization; he suggests seeding from s(156) minus 25 squares or s(132) minus 1 square, which "hasn't been tried yet";
  - s(155), s(181), s(208), s(209): seeded from larger n²−n records with 27–33 squares deleted.
- The SA output is then **analytically optimized** in Mathematica. The analytic_minimization page explains the method: it solves the contact equations plus nested Jacobian-determinant conditions when there are more free variables than constraints. This was first tested on Schadt's s(39).

### Public code: Ellsworth's tools repo
https://github.com/Davidebyzero/packing_squares_in_squares__tools contains:
- `refine_packing.cpp`: a GMP high-precision overlap fixer/refiner with SAT overlap tests, `--snap-grid`, `--shift` (moves an "I" strip), and `--check`. Commit messages say much of it was written with Claude Sonnet 4.5.
- `parse_svg_packing.py`: SVG to coordinates at full precision (Claude Sonnet 4.5 / Opus 4.5).
- `check_packing.py`: a Decimal overlap checker.
- `enum-45deg-rect-packings.cpp`: enumerates Göbel strip/square-type constructions.

**The SA program itself is not in the repo.** I read these files' text only and did not run them.

### Other methods named
- Brendberg: his own program plus manual polish (s(68), 2023). Program link: https://hohmiyazawa.github.io/squares_in_rects/squaredrawer.html
- Gensane & Ryckelynck 2004 program (s(29)).
- Joost de Winter (s(68), s(126)), in his own words:
  > "I essentially vibe coded an annealing program and tuned its meta-parameters. Then I fed it known-good packings from previous waves, which made it an evolutionary algorithm, and added constraints such as a fixed number of unrotated squares. Finally I fed it the existing record, with the symmetry assumptions of David Ellsworth's analytic solution relaxed …, and that worked … in minutes on my CPU."

  His notes describe the method as "structured perturbation of known good packings, followed by numerical refinement (sequential linear programming over separating-axis constraints, margin homotopy, clearance ladder)".
- Jake Loyd (s(68)), in his own words:
  > "L-BFGS (slightly modified LBFGSpp) … each square's x/y/rotation at a fixed container scale S. The objective combines squared overlap penalties based on SAT tests with squared wall protrusion penalties. Analytic gradients … An outer loop reduces S and reoptimizes … basin hopping, targeted perturbations and ruin and recreate … final polish … per square coordinate descent."

  He worked from random starts with no record input, on CPU.
- Cantrell and Hajba: hand or analytic constructions. The site says nothing about software for either. Cantrell sends "attached mathematical data" or exact analytic solutions by e-mail.
- AI-assisted entries:
  - Tej Stead with "Claude Fable 5", 16 Jun 2026: refined about 10 SA records (103, 110, 131, 132, 154, 155, 156, 180, 181, 182, …). These are small refinements; method not described.
  - Allen Chang with GPT‑5.6 Sol / GPT‑6 Astra: 68, 83, 87.
  - Ruilin Wang with ChatGPT‑6 Astra: 68.
  - hmbelvedere.com, probably AI: 69.
  - Dyachkov with Claude Opus 5 + Fable 5: 152.

---

## 4. n = 90, s(n²−n) = n, s(110), s(132)

- n = 72 (9²−9) and n = 90 (10²−10) are both **trivial** on the live page (s = 9 and s = 10). No nontrivial packing or attempt is shown on either the main or the compared page.
- s(110) (Cantrell, Feb 2025): "**Bounds the s(n²−n)=n conjecture to n < 11.**" So the conjecture s(n²−n)=n is still open exactly for n ≤ 10, and **s(90) is the largest open case**.
- s(132) (Arslanov–Mustafin–Shangitbayev, Mar 2019, E-JC doi:10.37236/8586): "Bounded the s(n²−n)=n conjecture to n < 12." Cantrell's s(110) then lowered the bound.
- s(156), s(182), and s(210) "Show s(n²−n)<n" for n = 13, 14, 15.
- The survey page (squares.html) keeps Friedman's historical text: "conjectured that s(n²−n)=n whenever n is small", with Cleemann's s(272)<17 as the first counterexample.
- s(110) progression:
  - 10.99679327… (Cantrell 2025, deficit 0.0032);
  - improved by Ellsworth mod-SA, Jan 2026;
  - refined by Stead/Claude Fable 5, Jun 2026;
  - not yet analytically optimized.
- s(132) progression:
  - AMS 2019: 11.99790;
  - Cantrell, Mar 2025;
  - Ellsworth mod-SA, Jan 2026;
  - Stead, Jun 2026;
  - live value 11.99137 (deficit 0.0086).
- Deficits 11 − s(110) = 0.0032 and 12 − s(132) = 0.0086 shrink quickly toward smaller n. This is consistent with s(90) being very hard or false to beat.
- **Has anyone tried n = 90?** The site says nothing. There is no "attempt", "invalid", or "closest attempt" entry for 90 like the one for s(103), where Ellsworth posted an invalid near-miss in Dec 2024–Jan 2025.
- The s(55) quote in §3 indicates that plain SA from randomness stalls just above the trivial size because of stacked rows/columns. That is exactly the regime of n = 90.
- Mark: **no evidence either way that someone has run SA at n = 90**.

---

## 5. Patterns

- **Construction genealogies that dominate 50–135:**
  1. Göbel's s(52)/s(65)/s(40)/s(18) families, extended by "L"s and pattern continuation: 52, 65, 67, 82, 84, 89, 101, 104, 109, 122, 124, 125, plus 85, 86, 127 from Friedman's 26/85 and Gustafsson's s(18).
  2. The **Cantrell 2002 s(37)-improvement technique**, adapted by Ellsworth in Nov–Dec 2024 and extended by the shared Jan 2025 "new technique": 88, 102, 123, 129, 130, 146, and many above 135.
  3. **Morandi s(69) / Cantrell s(53)** lineage: 106, 128, 152, 177, 205.
  4. **DeVincentis s(54)** lineage: 87, 107, 178.
  5. **n²−n−1 pattern (DeVincentis 2014)** and the **AMS 2019 n(n−1) packings**: 55, 71, 110, 131, 132, 154–156. These have since been reshaped by SA.
- **Found from scratch (random start) by search:** 50, 51, 103, and 105 (Schadt SA); 53 (Ellsworth mod-SA #3); 55 (Schadt, via a G&R-reimplementation seed from randomness). Jake Loyd's s(68) is the only from-scratch record from a non-SA (L-BFGS) method. Outside the range: 39 and 41 (Schadt).
- **Found by seeding search with a known record:** 71, 87, 110, 126, 129, 131, 132, plus the AI-era refinements of 68, 69, 83, 87, 126.
- **Untouched since 2014 or earlier:** 54, 70, and all the classical constructions.
- **Activity rate:** activity is bursty but large. Counting record changes in 50–135:
  - **Nov 2024–Mar 2025** (the Ellsworth/Cantrell/Hajba construction wave): about 25 changes.
  - **Dec 2025–Mar 2026** (Schadt SA and Ellsworth mod-SA): about 20 changes, at 50, 51, 53, 55, 68, 71, 87, 103, 105, 110, 126, 129, 131, 132.
  - **Jun 2026**: Stead/Claude Fable 5 refinements at 103, 110, 131, 132.
  - **Jul–Sep 2026** (AI-assisted outsiders): 9 changes at 68 (×4), 69, 83, 87, 126, and 152.

  Several records are flagged "Not yet analytically optimized". Pending items include s(88) (Schadt) and s(130), plus about 13 records awaiting the Jan 2025 technique. More changes in this range should be expected.

---

## Cross-check against our mirror (site/www/data/index.json, mtime 2026-09-30)

- All 194 n in the mirror exist on the live page, and the live page has no n that the mirror lacks.
- s values for 45–160 match the live page to the printed digits.
- Prose matches for every n, apart from parsing artifacts on my side: "Explore group" link text, and line-wrapping in s(106), which is "Found by David Ellsworth in November 2024, based on the s(53) found by David W. Cantrell in September 2002" in both.
- Mismatches also showed at n = 1453 and 2043, outside this range; not investigated.
- **Conclusion:** the mirror was current as of 2026-10-04 for this range, including the s(88) "Improvement by Thomas Schadt pending" line and Jake Loyd's 19 Sep 2026 s(68).
- The mirror holds only main-list prose. The richer history above (exact dates, seeds, SA statistics, method notes from de Winter and Loyd) comes from the SVG comments and the compared pages, which the mirror does not capture.

## Unverified / not found
- What Schadt's pending s(88) improvement is (value, date).
- What Schadt's SA program does: moves, objective, temperature schedule. It is not published anywhere I could find.
- What Ellsworth's modifications #1–#3 consist of.
- Which "about 13 additional packings" the Jan 2025 technique applies to.
- Whether anyone has attempted n = 72 or n = 90.
- What Tej Stead's Claude Fable 5 refinements did.
- compared3.html (n ≥ 197) was not read.
