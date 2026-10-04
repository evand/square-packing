# Conjectures-page literature double-check (2026-10-04)

Labels: [P] read in the primary source, [S] secondhand, [I] my inference.  Builds on `wishlist-literature-2026-10-03.md` and
`tasks/friedman-lit/REPORT.md`; nothing below contradicts them.  There are three additions the page should use: Chung–Graham's
prize challenges, the status of Erdős #106, and Ellsworth's implicit {s(n²+1)} ≥ ½ claim.  See "Flags" at the end.

## 1. Number of tilt angles in waste constructions — **nothing found on bounded-angle packings**
* **Erdős–Graham 1975** (scan text-extracted, pp.120–122) [P]: in the two border rectangles, every stack R(1,n+1) is "tilted at
  the appropriate angle θ" (one angle).  The leftover trapezoids T are cut into r = Ω(α^{4/11}) sub-trapezoids T_k, each with its own
  height s_k.  Their stacks touch the sloped top edge, and "the sum of the angles at the top vertices is O(α^{-1/11})", so adjacent
  stacks differ in angle.  So the number of distinct angles grows with α, at least polynomially [I].
* **Chung–Graham 2020** (preprint https://fanchung.ucsd.edu/wp/spacking.pdf) [P]: stacks in each B_i have different tilts.  The waste
  bound counts "the difference in the tilts of the two stacks".  **Bui 2508.04603** [P, grep of the PDF]: the angles θ_i, σ_1 and σ_2
  are recomputed at each step (θ_i ∈ Θ(x^{-ω})).  **McClenagan 2602.01484** [P]: it alternates angles φ, ψ, θ and θ′, and a second
  algorithm "undoes" the angle change at each step.  All of these use unboundedly many angles [I].
* None of these papers, nor Bui 2504.09489 ("only good squares matter", tilt ≤ 10⁻¹⁰), says anything about packings with a bounded
  number of angles.  Web searches found nothing either.  The question "is the waste linear when only finitely many angles are
  used?" looks unstated [I].

## 2. Algebraic nature of s(n) — **nothing in the literature; empirical data only**
* Ellsworth's page https://kingbird.myphotos.cc/packing/squares_in_squares.html [P] gives exact forms and minimal polynomials for
  best-known values.  Its degrees go up to **672 (s(83))**, then 198 (s(206)), 158 (s(179)) and 144 (s(108)).  37 entries are
  marked "Not yet analytically optimized".  A companion page, `squares_in_squares__analytic_minimization.html` [P], describes how to
  get exact roots: the side s and each distinct non-orthogonal angle are the variables, and contact equations are the constraints.
  It says nothing about algebraicity in general.
* **Rational non-integer best-known values exist** [P, Ellsworth]: s(293) ≤ 17+26/41 (Ellsworth Dec 2024, improved by Cantrell Jan
  2025, "the first record-setting packing found with a rational side length, thanks to the Pythagorean triple {20,21,29}"); s(230)
  ≤ 15+28/41; s(261) ≤ 16+28/41 (s(230) plus an "L"); and s(50) ≤ 7+4/7 (Schadt Dec 2025, triple {3,4,5}).  Non-record rational
  packings: s(104) and s(67) ≤ 8+5/7.  **No proved non-integer s(n) is rational**: s(5), s(10) = k+½√2 [P].
* DS7 (1998–2009) and Wikipedia have no discussion of algebraicity, degree, rationality or radicals [P, full-text grep].
* [I] Algebraicity of every s(n) follows from Tarski–Seidenberg: it is the minimum of a semialgebraic problem, with angles as
  (c,s), c²+s²=1.  No source states this.  No source discusses constructibility or radicals; s(11)'s octic is not solved in radicals.

## 3. Non-integer plateaus — **data only; no discussion found**
* Ellsworth [P] groups equal best-known values together: 147,148 (7+4√2); 232,233 (8+11√2/2); 264,265 (9+11√2/2); 290,291 (14+5√2/2,
  "the first s(n²+2) that has the same side length as the best known s(n²+1)"); 295,296 (17+√2/2).  The Göbel-square page [P]
  lists n = 148, 233, 265, so the largest member of each pair is a Göbel square or strip.  Ellsworth explains in general that a
  grouped label means "each smaller is represented by removing any square".
* DS7 (all versions) and Friedman's paper have no plateau statement [P grep].  Gardner's columns (Sci. Am. Oct/Nov 1979, Mar/Nov
  1980) were not accessible.
* Analogue: Erdős #106 (https://www.erdosproblems.com/106) [P]: "Erdős also asks for which n is it true that f(n+1)=f(n)", where f
  is the maximum sum of sides.  Erdős–Graham 1975 [P]: "We do not know as f(l) increases from 4k to 4k+4 how large the jumps are and
  where they occur."

## 4. Symmetry — **descriptive only; no theorems or conjectures**
* Ellsworth's "compared" page https://kingbird.myphotos.cc/packing/squares_in_squares__compared.html [P]: for s(17), a "Symmetric
  version found by David W. Cantrell in September 2023 … Didn't set an overall record, but is the best known symmetric packing"
  (4.68013 vs 4.67553).  It also lists several "alternative … with rotational symmetry" and "diagonally symmetric background"
  variants (n = 51-range, 66, 83, 86).  Main page [P]: Cantrell's "rotationally symmetric form" (s(146)), the "rotational symmetry
  technique" (s(230), s(261)), and "doubly semi-primitive" packings.  The s(n²−n−1) page [P] says the optimised odd-n pattern
  "becomes more chaotic … (although it retains symmetry)".
* No paper found on symmetry or asymmetry of optimal square packings.

## 5. Rigidity / local optimality / jamming — **informal definitions only**
* Ellsworth `squares_in_squares__rigid.html` [P]: "A packing is rigid when it cannot be continuously transformed into any other
  valid packing without changing the size of its enclosing square."  "The property of rigidity is rare among nontrivial best known
  packings.  Even the known recurring patterns only extend finitely before becoming inoptimal."  The page lists rigid packings: 5,
  11, 28 (Dec 2025), 40, the s(52) family (rigid alternatives for a ≡ 1 mod 5), and others.  It also uses "semi-rigid" (s(28)
  carousel) and "alternative" vs "rearrangement" (not reachable vs reachable by continuous motion).  DS7 [P]: s(40) and s(11) are
  called "rigid", with no definition.
* Gensane–Ryckelynck, "Improved dense packings of congruent squares in a square", DCG 34 (2005) 97–109 [S]: an inflation-based
  local optimiser.  No rigidity theory.
* Connelly-style infinitesimal rigidity, tensegrity and jamming results are all for disks and spheres (e.g. Connelly, EJC 2008
  "Rigidity of packings") [S].  **No square-in-square rigidity or jamming theory found.**

## 6. Problem collections
* **erdosproblems.com**: I fetched all 1221 problem pages and grepped them [P].  **Only #106** concerns squares packed in a square.
  "Draw n squares inside the unit square with no common interior point.  Let f(n) be the maximum possible sum of the side-lengths.
  Is f(k²+1)=k?"  Status: **DISPROVED (Lean)**.  "Claude Opus 5 (prompted by Silverstein) … proved that f(17)>4", hence
  f(k²+1) ≥ k + c/k for k ≥ 4.  f(10)=3 is open (page edited 28 Aug 2026).  **The waste problem W(x) is not on erdosproblems.com.**
* **formal-conjectures** (clone of 2026-10-03) [P]: `FormalConjectures/Wikipedia/SquarePacking.lean`.  Its open items are
  `least_eleven_square_packing_in_square` and `least_seventeen_square_packing_in_square` (IsLeast, answer(sorry)), plus the
  bounds 3.877084 and 4.6756.  Three squares in a circle is marked solved (Vu-Le, Kitamura).  The repo has no waste, Friedman or
  k²−k statements and no file for Erdős #106.
* **Brass–Moser–Pach** (2005) p.45 [S via Wikipedia, which cites it]: says the asymptotic waste is open "even for half-integer
  sides".  Not read.  **Croft–Falconer–Guy** problem **D4 "Packing equal squares in a square"** (p.111; DS7 cites pp.108–110) [P TOC
  only, https://www.gbv.de/dms/ilmenau/toc/188348913.PDF].
* **Chung–Graham 2020, concluding remarks** [P]: "one might be tempted to guess W(x)=O(x^{1/2}).  (The authors don't believe
  this.)"  **Challenge 1 ($100):** W(x)=o(x^{3/5}).  **Challenge 2 ($250):** O(x^{3/5−c}).  **Challenge 3 ($500):** W(N+½) ≫ N^{1/2+c}.
  Rewards are also offered for disproofs.  Bui 2508.04603 [P]: Conjectures 1–3 (trapezoid-packing exponents) and **Question 1: is
  o(x^{3/5}) possible when x−⌊x⌋ ∈ Θ(x^{−2/5})?**
* MathOverflow: no question on Friedman's conjecture, plateaus or W(x) found (search only) [S].

## 7. 2026 arXiv — **nothing new**
Searches turned up only items already listed: McClenagan, Bui v2, Karakuş, Sriswasdi, Dósa–Lángi–Tuza, and "Three squares in a
rectangle".  The arXiv API was rate-limited (HTTP 429), so coverage is limited to web search.  Non-arXiv items from 2026: Massaccesi's
blog has s(17) > 4.5058 (already known), and Ellsworth credits **Grigoriy Dyachkov (Sept 2026, "working with Claude Opus 5 and
Claude Fable 5")** with improved large-n packings, including s(2043).

## 8. Göbel strips and squares [P, Ellsworth sub-pages]
* **Strips:** a ∈ Z⁺, b = 1+⌊(a−1)√2⌋, n = (a+1)a+2+b, s = a+1+½√2.  "This yields the best known packing for all a < 44 except for
  a = 3" (n = 17).  "Göbel strips are inoptimal starting at s(2043) [a = 44] or possibly earlier."  The s(2043), s(2135) and s(3047)
  strips are beaten; s(1953) and s(2937) are still listed.  A second type starts at s(18): n = (a+1)a+3+b, s = 2+(b+1)½√2.
  Note that DS7's "a²+a+3+(a−1)√2" omits the floor.
* **Squares:** a−1 < ½b√2 < a+1, n = 2(a+1)a+b², s = a+1+½b√2 (better than trivial iff s < ⌈√n⌉).  Marked "Not optimal" at
  **n = 28** (beaten Dec 2025 by Ellsworth's rigid 5.82444), **1765** (Hajba/Ellsworth Nov 2024), **4009** and **9465**.  Listed
  without that mark: 40, 65, 89, 109, 124, 148, 233, 265, 376, 416, 445, 1544 and 1624.
* "Must eventually be beaten": the only statement is Ellsworth's empirical remark that recurring patterns "only extend finitely".
  [I] A proof is easy.  Both families have waste s² − n = Θ(s), e.g. ≈ a for strips.  Since W(x) = O(x^{7/11}) (Erdős–Graham), a
  square of side s−c already holds n squares for large s.  So every fixed-shape family with linear waste is eventually suboptimal.
  This is exactly Erdős–Graham's opening point, but no source applies it to Göbel families.

## Flags vs. the summary
1. **Chung–Graham 2020's $100/$250/$500 challenges** and their explicit *disbelief* in W = O(x^{1/2}) are missing from our notes.
   The authors' disbelief bears directly on Friedman's Conjecture 2 and on Erdős–Graham's "perhaps O(α^{1/2})".
2. **Erdős #106 is now disproved** (f(17) > 4, Lean-verified).  Our notes call it "a different problem", which is true, but give no
   status.
3. Ellsworth states an implicit pattern claim: **"Bounds {s(n²+1)} ≥ ½ to n < 42"** (s(1765), Hajba 2024; s(1850) bounds it to
   n < 43).  So someone once observed that the fractional part of best-known s(n²+1) is ≥ ½, and it fails at n = 42.  If the page
   has any s(n²+1) conjecture, cite this.
4. No contradictions with Friedman Conj. 1/2, Roth–Vaughan, Nagamochi, Arslanov et al., or Cantrell s(110) < 11 (Ellsworth:
   "Bounds the s(n²−n)=n conjecture to n < 11").

### URLs
Erdős–Graham https://static.renyi.hu/~p_erdos/1975-22.pdf · Chung–Graham 2020 https://fanchung.ucsd.edu/wp/spacking.pdf ·
Bui https://arxiv.org/abs/2508.04603 · McClenagan https://arxiv.org/abs/2602.01484 · Ellsworth main/rigid/strips/squares/compared:
https://kingbird.myphotos.cc/packing/squares_in_squares.html (+ `__rigid`, `__Göbel_strips`, `__Göbel_squares`, `__compared`,
`__n^2-n-1`, `__analytic_minimization`) · Erdős #106 https://www.erdosproblems.com/106 · formal-conjectures
https://github.com/google-deepmind/formal-conjectures/blob/main/FormalConjectures/Wikipedia/SquarePacking.lean · DS7 2009
https://www.combinatorics.org/files/Surveys/ds7/ds7v5-2009/ds7-2009.html · Wikipedia https://en.wikipedia.org/wiki/Square_packing
