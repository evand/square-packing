# Squares-in-squares: prior work, sources, and existing explainers

Research notes compiled 2026-09-08. Format: author, year, title, URL, one-line note. Items marked **[?]** are unverified or uncertain. Items marked **[in hand]** were already known; only new details are given.

## 1. Historical / primary sources

### Gardner and Göbel
- Gardner, M. (Oct 1979). "Some packing problems that cannot be solved by sitting on the suitcase," *Mathematical Games*, Scientific American 241(4). Original column posing the problem and Göbel's packings; includes "Gardner's conjecture" that n=11 is the first n needing non-45° tilts. Follow-ups (reader responses/addenda) appeared in the Nov 1979, Mar 1980, and Nov 1980 columns (titles of those issues are about other topics: "The random number omega...", "Graphs that can help cannibals...", "Taxicab geometry..."). Column list: https://en.wikipedia.org/wiki/List_of_Martin_Gardner_Mathematical_Games_columns
- Gardner, M. (1992). *Fractal Music, Hypercards and More*, W. H. Freeman, ch. "Packing Squares" pp. 289–306. Book reprint of the 1979 column with addenda (the source MathWorld and Friedman cite). Stanford's Gardner papers finding aid lists a "Packing of Squares" correspondence folder: https://archives.stanford.edu/download/sc0647.pdf **[?]** exact folder contents unverified.
- Göbel, F. (Frits) (1979). "Geometrical packing and covering problems," in A. Schrijver (ed.), *Packing and Covering in Combinatorics*, Math. Centre Tracts 106, pp. 179–199. Amsterdam. Source of the s(5)=2+1/√2 and s(10) packings, the "Göbel strip" family 2a²+2a+b² squares in side a+1+b/√2 (28, 40, 65, 89), and the first table of packings. First name "Frits" per Ellsworth's rigid-packings page; no online copy found.
- Gardner private communication to Friedman (1998) is ref [6] in the DS7 survey.

### Trump n=11 provenance
- Trump, W. (1979). n=11 packing, s ≈ 3.877084 (tilt ≈ 40.182°). Friedman's survey: "Many people have independently discovered this packing. The original discovery has been incorrectly attributed to Gustafson and Thule." Ellsworth's "compared" page says Evert Stenlund was "one of the 11 independent rediscoverers" before March 1980 (i.e., it arrived via Gardner's readers). Wikipedia: Walter Trump is a Bavarian retired teacher known for magic squares. https://en.wikipedia.org/wiki/Walter_Trump
- Trump, W. (Mar 2023). "Packing of 11 unit squares in a square with minimum size," ResearchGate preprint. https://www.researchgate.net/publication/368988287_Packing_of_11_unit_squares_in_a_square_with_minimum_size — Argues the 1979 packing is rigid and can't be improved within its combinatorial type; reports Thierry Gensane confirmed by email (Feb 2023) that their program couldn't beat it. (ResearchGate blocks fetch; summary from search snippets.)
- hohMiyazawa (20 Mar 2023). GitHub issue #3 on erich-friedman.github.io: "Is Gensane & Ryckelynck's 11-square packing identical to the 1979 packing of Walter Trump?" https://github.com/erich-friedman/erich-friedman.github.io/issues/3 — Asks whether G–R's 3.87708359 is a genuine improvement or a rounding artifact of Trump's; unresolved in the thread (Trump's 2023 note effectively answers: same packing).
- Gensane, T. & Ryckelynck, P. (2005, online 2004). "Improved dense packings of congruent squares in a square," Discrete Comput. Geom. 34, 97–109. https://doi.org/10.1007/s00454-004-1129-z **[in hand]** — "inflation" algorithm; new packings for n=11 (rediscovery), 29, 37; alternative optimal n=18. Authors at LMPA, Université du Littoral, Calais. No public code or coordinate site found.

### The Gardner-era contributors (all via Gardner correspondence; no independent web presence found)
- Cottingham, Charles F. (late 1979). Introduced diagonal width-2 strip packings for n≤49; his n=41 packing was reported in Gardner's Oct 1979 column (per Ellsworth "compared" page) and improved by Cantrell in 2005.
- Hämäläinen, Pertti (1980; letter dated 20 Apr 1980 per Stromquist 2003). Improved n=17 and n=18 (45° arrangements).
- Wainwright, Robert (1980). Improved n=19 (on Cottingham's idea).
- Stenlund, Evert (1980). Packings for all n≤100, incl. n=41, 66 (width-3 diagonal strip), 89; one of the 11 rediscoverers of Trump's n=11.
- Gustafsson, Mats (1981). Alternative optimal n=18 with tilt arcsin((√7−1)/4) ≈ 24.295° (spelled "Gustafson" in Friedman).
- Bidwell, John (1998, University of Hawaii undergraduate; some meme posts say 1997 **[?]**). n=17 (s≈4.6755, three angles) and n=29. Ref [1] of DS7 = private communication.
- Green, Trevor (2000, private communication). Lower-bound theorems 9–10 in DS7 and bounds for n=17–18, 22, 26–27, 28–30.
- Cantrell, David W. (2002–2005, and 2024–). n=19 (2002), 26, 37, 39, 41, 53, 54, 68, 70, 87, 88 (2005); first n=71 packing with s<9; many 2024–25 improvements on Ellsworth's page. Active on sci.math (circle packing threads) c. 2005.
- DeVincentis, Joe (Apr 2014). n=41 packing and an alternative n=54 extension (per Ellsworth's page).
- Hajba, Károly (2009: n=51, first s(51)<7+1/√2; Sep–Nov 2024: n=83, 102, 107, 108, 131, 172, 240).
- Morandi, Maurizio (2010). n=69.
- Brendberg, Sigvart (Jun 2023). A packing with s=8.80345993651653 (own program + manual optimization), later beaten by Schadt (Dec 2025). **[?]** which n — check Ellsworth's page.
- Schadt, Thomas (Dec 2025–2026). Simulated-annealing program "starting from randomness"; records for n=28, 29, 39, 41, 50, 51, 55, 68, 71, 103, 105, 126 etc.; a rational-side record 7+4/7 with {3,4,5} tilt. Program is **not public**; no web presence found beyond Ellsworth's page.
- Ellsworth, David (2023–2026). Maintains the record page; refined Schadt's finds; own records n=87, 88, 102, 106, 123, 126–136, 146 (first "doubly semi-primitive" record), 152, 172, 177, 205, 228, 234–236, 259, 266, 268. GitHub: https://github.com/Davidebyzero (no dedicated square-packing repo found there).

### Stromquist
- Stromquist, W. (1984). Memoranda, Daniel H. Wagner Associates: "Packing unit squares inside squares, I (six unit squares)" (11 Sep 1984), "II (ten unit squares)", "III" (Gardner's conjecture for n=11; 15 Nov 1984). PDFs exist: http://walterstromquist.com/papers/squares1.pdf, squares2.pdf, squares3.pdf (all return 200; not listed on his site index).
- Stromquist, W. (2003). "Packing 10 or 11 unit squares in a square," Electron. J. Combin. 10, #R8. https://www.combinatorics.org/ojs/index.php/eljc/article/view/v10i1r8 **[in hand]** — proves s(10)=3+1/√2 and s(11)≥2+4/√5≈3.789 (settles Gardner's conjecture).

### Surveys / reference pages
- Friedman, E. (1998→2009, maintained). "Packing unit squares in squares: a survey and new results," EJC DS7. https://doi.org/10.37236/28 ; versions: 1998 https://www.combinatorics.org/files/Surveys/ds7/ds7v1-1998.pdf , 2000 .../ds7v2-2000/ds7-2000.html , 2009 .../ds7v5-2009/ds7-2009.html ; maintained copy https://erich-friedman.github.io/papers/squares/squares.html **[in hand]**. Note: theorems 1–8 (unavoidable / almost-unavoidable sets) are Friedman's; 9–10 Green's; 11 Stromquist's. Version diffs are useful for a timeline of when bounds/packings appeared.
- Friedman, E. "Erich's Packing Center — squares in squares" https://erich-friedman.github.io/packing/squinsqu/ — now a one-line pointer: "The squares in squares page is now being maintained by David Ellsworth."
- Weisstein, E. "Square Packing," MathWorld. https://mathworld.wolfram.com/SquarePacking.html — table of best-known s(n), asterisks for proven; lists proven set {2,3,5,6,7,8,10,13,14,15,22,23,24,33,34,35} + squares; refs: Ellsworth, Erdős–Graham, Friedman DS7, Gardner 1992, Göbel 1979, Hoffman *The Man Who Loved Only Numbers* p.174, Roth–Vaughan. Last updated 2 Sep 2026. **[?]** its proven-list omits 46–48 (Bentz 2010 / Nagamochi).
- Wikipedia, "Square packing." https://en.wikipedia.org/wiki/Square_packing — refs: Brass–Moser–Pach 2005, Nagamochi 2005, Bentz 2010, Friedman 2009, Stromquist 2003, Gensane–Ryckelynck 2005, MathWorld, Erdős–Graham 1975, Chung–Graham 2009, Wang–Dong–Li 2016, McClenagan 2024, Roth–Vaughan 1978, Abrahamsen–Stade 2024. States proven n = squares, 1–10, 13, 14, 15, 24, 34, 35, 46, 47, 48 (omits 22, 23, 33 from Bentz 2016 **[?]** — check current revision), n=11 lower bound 2+4/√5, Bidwell n=17. No rattler/symmetry discussion.
- Croft, Falconer & Guy (1991). *Unsolved Problems in Geometry*, Springer, section D5 (packing squares). Ref [2] in DS7.
- Brass, Moser & Pach (2005). *Research Problems in Discrete Geometry*, Springer. Notes the asymptotic growth rate of wasted area is open.
- OEIS: **no sequence for s(n) itself** (values are irrational). Related: A360610 (T(n,k) = number of side-k squares fitting in side-n square; cites DS7); A374528 (record n for squares in circles, Pfoertner 2024); A124484 / A374505 / A373008 (axis-parallel unit squares in a circle of radius/diameter n, cite Erdős–Graham 1975); A084616 (circles in square of area n). OEIS search: https://oeis.org/search?q=%22unit+squares%22+packing+square

## 2. Existing visualizations, explainers, and tools

### Record databases / interactive
- Ellsworth, D. "Squares in Squares" https://kingbird.myphotos.cc/packing/squares_in_squares.html — SVG of every best-known packing for n up to ~1800 (?) with high-precision s(n), exact algebraic forms where known, dated attributions, and an "SVG Edit Mode" where squares can be dragged. Sub-pages: triangular table view (`__triangular_table.html`), older/alternative packings (`__compared.html`), rigid packings (`__rigid.html`), Göbel strips (`__Göbel_strips.html`), Göbel squares (`__Göbel_squares.html`), s(n²−n−1) page (`__n^2-n-1.html`), and a mirror of DS7 (`squares.html`). Does NOT show: proof/lower-bound history, method explanations, timeline of when records fell, animations of proofs.
- Quilez, I. (2023). "Squares" https://iquilezles.org/maths/squares/ — static page with coordinates/angles for Bidwell's 17 (s=4.6755389909…) and a two-angle alternative (s=4.67765236…); no interactivity, no proofs.
- Crane, K. (17 Feb 2023). X thread https://x.com/keenanisalive/status/1626707731339173888 (Thread Reader: https://threadreaderapp.com/thread/1626707731339173888.html) — reproduces the n=17 packing in <5 min in the Penrose constraint-graphics editor (penrose.cs.cmu.edu/try/); a demo of Penrose, not a search tool.
- Observable notebooks (Esperança "Square Packing", "Hierarchical Square Packing"; Muyskens "Square Packing") — **different problem** (greedy packing of unequal squares / treemaps); not relevant.
- Wolfram Demonstrations "Packing Squares with Side 1/n" — different problem (harmonic squares); not relevant.
- Hobby GitHub repos: smartycope/SquarePacking (RL for n=11), leove4/Square-packing-simulation (matplotlib GUI simulation), simonpasi96/Square-packing, bmetenko/SquarePacking — none find records; all toy-level.
- 3D-print models of n=11 (MakerWorld "Square Packing n = 11", Printables "11 Squares Packing Puzzle") — physical puzzles.

### Videos
- Deckard (YouTube @deckardv), 27 Sep 2025, "Packing Squares Inside The Smallest Square Possible" https://www.youtube.com/watch?v=uL5wuiy34rs — the original of the two-video series; credits Friedman and Ellsworth. HN thread: https://news.ycombinator.com/item?id=45405476 .
- Deckard, ~mid-Aug 2026, "NEW SQUARE PACKING SOLUTIONS" https://www.youtube.com/watch?v=beFC3qRBG4E **[in hand]** — sequel covering the Schadt/Ellsworth 2025–26 records; points to Ellsworth's page. **[?]** exact upload date (search said "3 weeks ago" on ~8 Sep 2026).
- EddAardvark, "Packing 17 squares into a larger square" https://www.youtube.com/watch?v=MGHNIj6qOXA — short explainer of the n=17 meme. **[?]** date.
- Andy Math, "Square Packing" https://www.youtube.com/watch?v=jToq8C89r0I — short explainer. **[?]** date.
- OneMinuteThings, "Square packing is weird." https://www.youtube.com/watch?v=81DCjd5DPMs — one-minute short. **[?]** date.
- Combo Class, "The Insane World of Polygon Packings" https://www.youtube.com/watch?v=jWT08JVb-fk — broader polygon-packing survey, includes squares. **[?]** date.
- YouTube channel @17squaresinasquare "The Optimal Way To Pack 17 Squares Into A Square" — meme channel.
- No Numberphile / Stand-up Maths / 3Blue1Brown video on this problem was found.

### Social / meme / press
- Piker, D. (@KangarooPhysics), 10 Dec 2021 thread (origin) and 14 Feb 2023 repost https://x.com/KangarooPhysics/status/1625423951156375553 (5.1M views); follow-up on rigidity https://x.com/KangarooPhysics/status/1625627324296175617 ("not rigid – 3 of them can slide, and one even has room to wiggle"). Know Your Meme entry with timeline: https://knowyourmeme.com/memes/17-squares-in-a-larger-square .
- Munroe, R. (20 Feb 2023). xkcd #2740 "Square Packing" https://xkcd.com/2740/ ; explainxkcd https://www.explainxkcd.com/wiki/index.php/2740:_Square_Packing — joke about n=11 and a hydraulic press; explainer links G–R 2005, kingbird, Friedman.
- Hacker News: "Squares in Squares" (kingbird), 17 Feb 2023 https://news.ycombinator.com/item?id=34809023 ; Deckard video Sep 2025 https://news.ycombinator.com/item?id=45405476 ; Massaccesi bound Aug 2026 https://news.ycombinator.com/item?id=49390775 (30 pts, 8 comments; discussion of why weighted-point certificates are less "visual" than packings).
- Manifold markets: "Does a better packing of 17 squares into a square exist?" (Aspen) https://manifold.markets/aspentropy/does-a-better-packing-of-17-squares (41% YES); "Will proof of an optimal packing of 17 or fewer squares ... before 2024?" (Jim Hays) https://manifold.markets/JimHays/will-proof-of-an-optimal-packing-of — resolved YES because Friedman's page was updated to mark n=13 as proved (Bentz 2010) after community contact; a nice anecdote for the site.
- Plevris, V. (Medium, 2025?) "Eleven Squares, One Tiny Gap, and a Problem Still Unsolved" https://vplevris.medium.com/eleven-squares-one-tiny-gap-and-a-problem-still-unsolved-c6f47b447cfb — popular n=11 piece (403 on fetch; content unverified).
- "West Indian Archie" (Medium) "Deeply unsettling asymmetric patterns in mathematics: optimal packing of 17 squares" https://medium.com/@Mazequizzing/deeply-unsettling-asymmetric-patterns-in-mathematics-optimal-packing-of-17-squares-0171a792f165 — popular piece (unverified).
- Quora, TikTok "17 squares in a square" discoverable pages — meme-level only.
- Grokipedia / HandWiki / Brainbound.blog copies of Wikipedia — ignore.

### Gap analysis (what nobody has built)
- No page shows a **timeline** of record-breaking per n or an animated "how the record fell" view.
- No page shows the **lower-bound side**: neither Friedman's Table 2 nor Ellsworth's page visualizes unavoidable-point sets, weighted-point certificates, or the gap s_lower(n) vs s_upper(n) as a chart.
- No interactive explainer of proof methods (unavoidable sets, almost-unavoidable sets, Bentz's continuously varying families, Nagamochi's rectangle lemma, weighted/LP certificates).
- Rigidity/rattler info exists only as Ellsworth's list and Piker's tweet; no visualization of the free squares.

## 3. Proven lower bounds and exact results (timeline)

### Small n (exact values / specific bounds)
- Göbel 1979: s(2)=s(3)=2, s(5)=2+1/√2 (per DS7 figs 23–24; Friedman's DS7 restates as Theorems 1–2).
- Stromquist 1984 memos (unpublished): s(6)=3, s(10)=3+1/√2; claimed s(14)=s(15)=4, s(24)=5.
- Friedman 1998 (DS7 v1): unavoidable-set proofs of s(8)=3, s(15)=4, s(24)=5, s(35)=6; almost-unavoidable-set proofs of s(7)=3, s(14)=4; s(13)≥3.8437; s(19–20)≥6√2−4≈4.4852; s(21)≥4.7438.
- El Moumni, S. 1999. "Optimal packings of unit squares in a square," Studia Sci. Math. Hungar. 35(3–4), 281–290 — independent proofs of s(7)=s(8)=3, s(15)=4.
- Green, T. 2000 (private comm.): s(17)=s(18)≥(40√2+19)/17≈4.4452; s(22)≥2√2+2; s(26–27)≥2√2+(27+2√10)/13≈5.3918; s(28–30)≥2√2+6/√5≈5.5117; general Theorems 9–10 for n²+1 and n²+⌊n/2⌋+1.
- Kearney, M. J. & Shiu, P. 2002. "Efficient packing of unit squares in a square," Electron. J. Combin. 9, #R14 https://www.combinatorics.org/ojs/index.php/eljc/article/view/v9i1r14 — first published proof s(6)=s(7)=3 by a "duality" method; also n_r bounds for s(n²+1)≤n+1/r (n_2≤43).
- Stromquist 2003 (EJC 10 #R8): s(10)=3+1/√2; **s(11) ≥ 2+4/√5 ≈ 3.7888** (hence s(12) too, since s is nondecreasing). This is still the best n=11/12 lower bound as of 2026 (nothing newer found).
- Nagamochi, H. 2005. "Packing unit squares in a rectangle," Electron. J. Combin. 12, #R37 https://www.combinatorics.org/ojs/index.php/eljc/article/view/v12i1r37 — s(n²−2)=s(n²−1)=n for all n (so 23, 34, 47, 48, 62, 63, 79, 80, ...); general bound s(N)≥min{⌈√N⌉, √(N−2⌊√N⌋+1)+1}.
- Bentz, W. 2010. "Optimal packings of 13 and 46 unit squares in a square," Electron. J. Combin. 17, #R126 https://www.combinatorics.org/ojs/index.php/eljc/article/view/v17i1r126 — s(13)=4, s(46)=7 (i.e., s(m²−3)=m for m=4,7).
- Bentz, W. 2016 (arXiv 12 Jun 2016; v2 dated 20 Oct 2018). "Optimal packings of 22 and 33 unit squares in a square," arXiv:1606.03746 https://arxiv.org/abs/1606.03746 — s(22)=5, s(33)=6 via continuously varying families of unavoidable sets. **[?]** No journal reference found on arXiv or dblp; may be unpublished — cite as arXiv.
- Green/Burns/Massaccesi/Mira 2026 on n=17: see §5. Summary chain: 4.4452 (Green 2000) → 4.4811 (Burns + ChatGPT, 6 Aug 2026) → 4.5058 (Massaccesi LP, 21 Aug 2026) → 4.607 → **4.6130286** (Mira-acc/17squares, exact rational certificate, 1,620 atoms; date **[?]**, Aug–Sep 2026). Upper bound unchanged: Bidwell 4.6755. None peer-reviewed; all self-verifying scripts.

### Asymptotics (wasted area W(x) in a square of side x)
- Erdős, P. & Graham, R. L. 1975. "On packing squares with equal squares," J. Combin. Theory Ser. A 19, 119–123 https://www.sciencedirect.com/science/article/pii/0097316575900990 (PDF: https://www.math.ucsd.edu/~fan/ron/papers/75_06_squares.pdf) — W(x)=O(x^{7/11}) via tilted strips; asked whether W(x)=o(x^{1/2}) possible... (note: Erdős offered $ prizes for the exponent).
- Roth, K. F. & Vaughan, R. C. 1978. "Inefficiency in packing squares with unit squares," JCTA 24(2), 170–186 https://doi.org/10.1016/0097-3165(78)90005-5 — lower bound W(x) ≫ (x·|x−round(x)|)^{1/2}; in particular for x=n+1/2, W ≫ √x. So exponent lies in [1/2, 3/5].
- Chung, F. & Graham, R. 2009. "Packing equal squares into a large square," JCTA 116(6), 1167–1175 https://doi.org/10.1016/j.jcta.2009.02.005 (PDF https://mathweb.ucsd.edu/~ronspubs/09_03_square_packing.pdf) — W(x)=O(x^{(3+√2)/7} log x) ≈ O(x^{0.6306} log x).
- Wang, S., Dong, T. & Li, J. 2016. "A new result on packing unit squares into a large square," arXiv:1603.02368 — W(x)=O(x^{5/8}); also a covering corollary.
- Chung, F. & Graham, R. 2020. "Efficient packings of unit squares in a large square," Discrete Comput. Geom. 64, 690–699 https://doi.org/10.1007/s00454-019-00088-9 — claimed W(x)=O(x^{3/5}).
- Arslanov, M. Z. & Bui, H. D. 2025. "Note on 'Efficient packings of unit squares in a large square'," Discrete Comput. Geom. https://doi.org/10.1007/s00454-025-00767-w — identifies an angle-computation error invalidating the 2020 O(x^{3/5}) claim.
- Bui, H. D. 2025 (arXiv 13 Apr 2025). "Square packing with asymptotically smallest waste only needs good squares," arXiv:2504.09489 — reduction: only near-axis-aligned squares matter for lower-bound analysis.
- Bui, H. D. 2025 (arXiv 6 Aug 2025, rev. Mar 2026). "Square packing with O(x^{0.6}) wasted area," arXiv:2508.04603 — repairs the 3/5 exponent by a new construction.
- McClenagan, R. 2024. "Asymptotic square packing problems," MSc thesis, Univ. of Northern British Columbia, DOI 10.24124/2024/59553; and 2026 (arXiv 1 Feb 2026) "Optimally packing a large square by unit squares," arXiv:2602.01484 — independent O(x^{3/5}) proof. Current status: 1/2 ≤ exponent ≤ 3/5.
- Tangential (different "Erdős square packing" problem: max sum of side lengths of n squares in a unit square, f(k²+1)=k conjecture): Praton 2005 arXiv:math/0504341; "A note on the Erdős conjecture about square packing" arXiv:2411.07274 (2024); "An equivalence between Erdős's square packing conjecture and the convergence of an infinite series" arXiv:2506.23284 (2025). Don't confuse with s(n).
- Tangential: Abrahamsen & Stade 2024 (FOCS) hardness of packing polygons with unit squares; Montanher, Neumaier, Markót et al. 2019 "Rigorous packing of unit squares into a circle," J. Global Optim. 73, 547–565 (interval-arithmetic proof for 3 squares in a circle — the only computer-rigorous optimality method in the unit-square literature); Mok & Wong 2023 "Least optimal square packing in a square," Hang Lung Math Awards 10, 39–94 (the *dual* worst-case problem); "Covering a square by congruent squares" arXiv:2601.16535 (covering analogue).

## 4. Rattlers, rigidity, symmetry

- Ellsworth, D. "Squares in Squares: Rigid packings" https://kingbird.myphotos.cc/packing/squares_in_squares__rigid.html — defines rigid (no continuous motion without enlarging the container), semi-rigid (e.g., "carousel" group in n=28), alternative packing vs rearrangement; lists rigid record packings for n = 5, 11, 18, 28, 40, 50, 52, 149, 200, 203, 233, 265, 296, 299, 335, 373, 493, 740, 1037, 1044, 1384, 1781. Also defines "doubly semi-primitive" (non-rotated squares only on left/right edges, rotated squares protruding into top/bottom) — a structural, not symmetry-group, classification. The term "rattler" is **not** used anywhere in the squares literature; it's borrowed from sphere/disk-packing (Connelly; Torquato). Ellsworth's main page labels n=5, 11, 40 "[Rigid]".
- Piker 2023 tweet (above): n=17 has 3 sliding squares and one with wiggle room — the only explicit "free squares" statement for n=17 found.
- Trump 2023 (above): rigidity argument for n=11.
- Connelly, R. "Packings of circles and spheres" lecture notes https://pi.math.cornell.edu/~connelly/PackingsIII.IV.pdf — rigidity theory background (for circles), useful for terminology.
- No paper classifying square packings by symmetry group (C1/C2/C4/D1/D2/D4) was found; this would be new content. (Bentz's proofs and Friedman's unavoidable sets do use D4 symmetry of the container; Massaccesi's LP uses D4 symmetry reduction of weights.)

## 5. 2025–2026 news

- Ellsworth page updates (Sep 2024 → 2026): Hajba (Sep–Nov 2024), Cantrell (Nov 2024→), Schadt (Dec 2025→; ~a dozen records incl. rational side 7+4/7), Ellsworth (Dec 2025–Feb 2026, modified copies of Schadt's SA program; many records incl. n=146 first doubly-semi-primitive record). Ellsworth's page notes the s(n²−n)=n conjecture is now false for n=11 (a 110-square packing with s<11) — bounding the conjecture to n<11. **[?]** exact statement/date; verify on the page.
- Deckard video 2 (~Aug 2026) covers these.
- Burns, S. 5 Aug 2026. "The n=17 Square Packing Problem – Introduction" https://sam-burns.com/posts/n17-square-packing-problem-intro/ ; 5 Aug 2026 "A Near-Record Arrangement" https://sam-burns.com/posts/n17-square-packing-near-record-arrangement/ (s≈4.6776483, 11 axis-aligned + 6 tilted at 39.63°, found by numerical search; cf. Quilez's two-angle 4.6776524); 6 Aug 2026 "Proposing a Better Lower Bound for n=17" https://sam-burns.com/posts/proposing-better-lower-bound-for-n17-square-packing/ — s(17)≥4.4811 via 268 weighted points (total 16.9476, every unit square captures ≥1.0003); certificate "developed by ChatGPT 5.6 Pro" with Burns as facilitator; Python verifier; code/model at https://github.com/sam-bee/squarl .
- Massaccesi, G. 21 Aug 2026. "Another Better Lower Bound for n=17 Square Packing" https://gus-massa.blogspot.com/2026/08/another-better-lower-bound-for-n17.html — s(17)≥4.5058 via 168 weighted points on a 29×29 grid, weights from LP (scipy linprog) with D4 reduction; exact-rational Python verifier + Racket visualizer. Companion: "Linear Programing for Square Packing" (21 Aug 2026) https://gus-massa.blogspot.com/2026/08/linear-programing-for-square-packing.html . HN https://news.ycombinator.com/item?id=49390775 .
- Mira-acc/17squares (GitHub) https://github.com/Mira-acc/17squares — s(17) > 4.613028635886 (= sqrt of an explicit rational), 1,620 rational-weighted atoms, total mass 16.99798 < 17, swept over 2,881 rational directions with a 0.99985-side core, min mass 1.000002103; exact rational arithmetic (Python + C++17/Boost), no Lean; earlier packages 4.468292 and 4.607028598640 retained. README credits "Joshua Levy's squares project (CC BY 4.0)" for an initial 4.59 measure and a "weighted principle and sharp dilation lemma" from prior work; **no mention of Fort, Burns, or Massaccesi by name in the README and no AI mention** — the attribution of this repo to Burns/Massaccesi/Fort/"Mira" in the task brief is **[?]** unverified; I could not find any Stanislav Fort post on the topic (his blog https://stanislavfort.com/blog/ and X account returned nothing relevant in search). "Joshua Levy" also not found elsewhere.
- BoltonBailey, 2 Sep 2025. formal-conjectures issue #646 "Seventeen square packing problem" https://github.com/google-deepmind/formal-conjectures/issues/646 — proposes Lean formalization of s(17); status closed; **[?]** whether a Lean file was merged (no file found by search).
- AI-assisted packing results: AlphaEvolve (2025) improved circles-in-square and 11 hexagons-in-hexagon, **not** squares-in-squares; Berthold et al. 2026 arXiv:2605.04850 (global MINLP solvers) covers circles/polygons/Platonic solids, not unit squares; ThetaEvolve, ShinkaEvolve etc. target circle packing. So the only AI-assisted squares-in-squares result found is Burns's ChatGPT-derived n=17 certificate (Aug 2026). Sam Burns's "squarl" also trained an ML model to nudge random arrangements but records were found by numerical optimization.
- Wikipedia/Grokipedia mention nothing of the 2026 bounds yet **[?]** (fetched summary didn't show them).

## Open verification items
1. Exact upload dates for Deckard's second video and the EddAardvark/Andy Math/OneMinuteThings/Combo Class videos (YouTube pages don't render in fetch; oembed gives no date).
2. Whether Bentz's 22/33 paper was ever journal-published.
3. Provenance of Mira-acc/17squares and the roles of Fort/"Mira"/Joshua Levy; date of the 4.6130 certificate.
4. Which n Brendberg's June 2023 packing (s≈8.8035) is, and the exact n=110 conjecture statement on Ellsworth's page.
5. Contents of Stromquist's three 1984 memoranda PDFs (exist at walterstromquist.com/papers/squares{1,2,3}.pdf; not read).
6. Whether Trump's 1979 packing reached Gardner directly or via the March 1980 addendum (Ellsworth's compared page implies the latter via "11 independent rediscoverers before March 1980").
