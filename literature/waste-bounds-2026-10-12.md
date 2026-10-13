# Classical waste bounds and Stromquist's memoranda, read in the primary sources (2026-10-12)

Until now our notes described Erdős–Graham, Roth–Vaughan and Stromquist's 1984 memoranda secondhand.  This note
reads them, plus Bui arXiv:2504.09489 and Bašić–Slivková 2018.  It records the exact statements, the proof ideas,
what they give for our exact s(n) work, and a list of fixes to our files (§6).  No dated note or reader-facing page was
edited.

Labels: **[read]** = read in the primary source (page or theorem cited); [S] = secondhand; [I] = my inference.
Notation: s(n) = side of the smallest square holding n unit squares.  W(α) = α² − (max number of unit squares in a
square of side α), as defined by Erdős–Graham and by Roth–Vaughan.  ‖α‖ = distance from α to the nearest integer,
and [α] = ⌊α⌋.

## 1. Sources

| Source | Status | Where |
|---|---|---|
| P. Erdős, R. L. Graham, "On packing squares with equal squares", J. Combin. Theory Ser. A 19 (1975) 119–123 | **[read]**, all 5 pages | PDF supplied by Evan (`../private/literature/1-s2.0-0097316575900990-main.pdf`, not copied) |
| K. F. Roth, R. C. Vaughan, "Inefficiency in packing squares with unit squares", J. Combin. Theory Ser. A 24 (1978) 170–186, doi 10.1016/0097-3165(78)90005-5, received 1 Nov 1976 | **[read]**, all 17 pages | PDF supplied by Evan (`../private/literature/1-s2.0-0097316578900055-main.pdf`) |
| W. R. Stromquist, "Packing Unit Squares Inside Squares, I (Six Unit Squares)", internal memorandum to "450 File (Professional Leave)", 11 Sep 1984, 19 pp. | **[read]** (scanned page images, no text layer) | http://walterstromquist.com/papers/squares1.pdf |
| — "II (Ten Unit Squares)", Daniel H. Wagner Associates, Paoli PA, internal memorandum, 15 Oct 1984, 15 pp. | **[read]** | …/squares2.pdf |
| — "III (Cases with n ≤ 65 and Martin Gardner's Conjecture for n = 11)", Daniel H. Wagner Associates, 15 Nov 1984, 13 pp. | **[read]** | …/squares3.pdf |
| H. D. Bui, "Square Packing with Asymptotically Smallest Waste Only Needs Good Squares", arXiv:2504.09489 (v1 only, 13 Apr 2025; "submitted to CCCG 2025"; no journal reference) | **[read]**, all 7 pages | arXiv |
| B. Bašić, A. Slivková, "On optimal piercing of a square", Discrete Applied Mathematics 247 (2018) 242–251, doi 10.1016/j.dam.2018.03.048 (University of Novi Sad) | bibliographic record **[read]** (Crossref); Theorems 7, 9, 10 and Proposition 8 **[read]** in a PDF-to-text copy (jlevy's archive, `~/math/_untrusted-third-party/jlevy-squares/packing/resources/papers/basic-slivkova-2018-optimal-piercing-square.raw.md`), not the publisher's page | §5 |

## 2. Exact statements

### Erdős–Graham 1975 [read]
* Definition (p.120): for real α > 0, W(α) = α² − sup|𝒫| over packings 𝒫 of unit squares into a square of side α.
  The squares are "nonoverlapping open unit squares" (p.119).
* **Theorem (p.120, eq. (1)): W(α) = O(α^{7/11}).**  No explicit constant and no explicit threshold on α.  In their
  proof, "f = Ω(g)" is defined as two-sided (c·g < f < c′·g), i.e. it means Θ.
* Corollary (p.120, p.122): f(k² + [c·k^{7/11}]) > 4k, where f(l) is the maximal sum of *perimeters* of l squares packed
  in a unit square.  This is Erdős's sum-of-perimeters problem (now Erdős #106), not our s(n).
* Remarks (p.122): "It is rather annoying that we do not at present have any nontrivial lower estimate for W(α).  Indeed
  we cannot even rule out the possibility that W(α) = O(1).  Perhaps the correct bound is O(α^{1/2})."  They also say
  (p.120) that they do not know how large the jumps of f(l) are as f increases from 4k to 4k + 4, or where they occur.
* The paper offers no prizes, poses no o(α^{1/2}) question and makes no higher-dimensional statement.  Its closing
  questions are about convex curves 𝒞 and f(𝒞, k).

### Roth–Vaughan 1978 [read]
* Waste definition (p.170, eq. (1)): w(α) = α² − sup|𝒜|, for real α ≥ 1, exactly as in Erdős–Graham.  Packings are
  closed unit squares with disjoint interiors inside the closed square S(α) (p.171).
* **Theorem (p.170, unnumbered).  Suppose that α(α − [α]) > 1/6.  Then w(α) ≫ (‖α‖·α)^{1/2}.**
* **Explicit constant.**  The proof (p.173, eqs. (5)–(6), conclusion on p.186) shows that every packing has waste
  ≥ c₀(‖α‖α)^{1/2} with c₀ = c^{10} = 10^{−100}, where c = 10^{−10} (eq. (3)) is the tilt threshold for "good" squares.
  This holds when ‖α‖α > 1/6, assumption (2).  In the remaining case of the hypothesis (α just below an integer, with
  α(1 + [α] − α) ≤ 1/6), p.171 gives w(α) ≥ ½ ≥ (α(1 + [α] − α))^{1/2}.  So **α(α − [α]) > 1/6 ⇒ W(α) ≥
  10^{−100}·√(α‖α‖)** is a correct reading.  Stromquist's memo III, p.12, states it in exactly this form.
* The hypothesis matters.  If α(α − [α]) ≤ 1/6 (α just *above* an integer k), RV note (p.170) that sup|𝒜| = k², so
  W(α) = α² − k² ≈ 2α‖α‖.  This is "essentially best possible", and it is *smaller* than 10^{−100}√(α‖α‖) once
  α‖α‖ < ~10^{−200}.  The unconditional form "W(α) ≥ c√(α‖α‖) for all α" is therefore false as stated, although only
  absurdly close to an integer.
* Special case (abstract): side n + ½ ⇒ waste ≫ n^{1/2}.
* Context on p.170: Montgomery has an unpublished argument improving 7/11 to (3 − √3)/2 + ε = 0.633974… + ε.  RV are
  "dubious" of Erdős–Graham's O(α^{1/2}) but say that if it is true, then it is essentially best possible.  They do not
  call their bound "matching" the 7/11 upper bound.

### Bui, arXiv:2504.09489 [read]
* Definitions (§2): W(x) = min over packings A of x² − |A|; good square = inclination ≤ c = 10^{−10} (RV's c);
  W*(x) = min waste over packings using good squares only.
* **Proposition 22: there is a constant c₂ > 1 such that for every packing A there is a good packing A′ with
  W(A′) ≤ c₂·W(A).  Theorem 1: W(x) = Θ(W*(x)).**  The constants are not explicit.
* Bui states RV's result only as "W(x) ∉ o(x^{1/2})" (abstract, §1).  That is weaker than RV's pointwise statement.
* Open problems (§6): can the tilt threshold be taken as O(1/log x) or O(x^{−ε}); Conjecture 1 (the radius r in
  Prop. 21 can be arbitrarily small, with waste ≥ c′·min(1, r²)·θ); algorithmic questions about the flow construction.
* It is not a new lower bound on W.  It simplifies RV's Lemma 4 (the bound on the number of bad squares) away.

### Stromquist memoranda, 1984 [read]
* **I (11 Sep 1984).  Theorem (p.3): six pairwise non-intersecting "blocks" (open squares of side 1 + ε, 0 < ε < 10^{−4})
  cannot lie in the interior of a square of side 3.**  So s(6) = 3, and with it s(7) = s(8) = s(9) = 3 (p.2: "This will
  settle the case of n = 6, and the cases of n = 7, n = 8, and n = 9 as well").  §2 also gives self-contained solutions
  for n ≤ 5; for n = 5 the proof is credited to Göbel.  Section 3 is a key-point argument: points at spacing ½ are
  allocated to the blocks, and the only surviving allocation (Fig. 14) is refuted in Fig. 15.
* **II (15 Oct 1984).  Theorem (p.2): ten pairwise non-intersecting blocks cannot lie in the interior of a square of
  side 3 + ½√2.**  So s(10) = 3 + 1/√2 ≈ 3.7071.  Four "nonavoidance lemmas" (§1), e.g. Lemma 2: a block whose centre
  lies in a triangle with sides ≤ 1 contains a vertex of the triangle.  The final step: every block must contain one of
  13 points (Fig. 14), all of which are denied to one block.  Fig. 1 shows three different optimal packings.
* **III (15 Nov 1984), in three parts.**
  1. Table 1 (p.2) lists the best known s for n ≤ 65.  New packings: n = 18, s = 7/2 + ½√7 ≈ 4.823 (θ ≈ 24.295°), and
     n = 26, s ≈ 5.650629 (θ ≈ 27.583°).  Also a 19-square packing, s ≈ 4.888 (Fig. 2d, θ ≈ 23.944°), which is not a
     record.  The 3.877 packing for n = 11 is credited to Gustafsson & Thulin (via Gardner, Nov 1980).
  2. **Theorem (p.7): with s = 2 + (4/3)√2 ≈ 3.886, eleven blocks of orientation 0° or 45° cannot lie inside the
     square.**  This settles Gardner's conjecture that n = 11 is the first case needing a tilt other than 45°.  The proof
     uses ten marked points plus a further 12-point set (Figs. 5–8).
  3. Parenthetical, p.10: "By reducing the value of s, essentially the same argument can be made to work for general
     packings.  The result is that 11 unit squares cannot be packed inside a square with s < 2 + (4/5)√5 ≈ 3.789."
     This is a 1984 sketch of the bound s(11) ≥ 2 + 4/√5 that Stromquist published in EJC 10 (2003) #R8.
  * Lower bounds (p.6): Göbel proved n ≤ 5; memos I–II extend to n ≤ 10.  "Certain other cases (n = 14, n = 15, n = 24)
    are fairly easy to settle by elementary means" (no proof given).  Göbel attributes s(7) = 3 to G. Bajmóczy.  "It
    seems not to be known, even, whether there is a nontrivial packing for some n of the form n = k² − 1."
  * §3 (pp.10–12) only restates Erdős–Graham (W(s) = O(s^{7/11})), Montgomery (exponent 3/2 − ½√3 + ε) and RV
    (W(s) ≥ 10^{−100}√(‖s‖s) provided s(s − [s]) > 1/6), and remarks that RV caps any improvement of the exponent at ½.
    It also notes (p.3) that Erdős–Graham "compel the existence of nontrivial packings with n = k² + k" for all large k,
    and that Gardner quotes Graham as having found one at n = 40² − 40 = 1560.
* The memoranda prove nothing asymptotic.  Their content is s(n) for n ≤ 10, the 45° theorem for 11, the sketched
  s(11) ≥ 3.789, and the table.

## 3. Proof ideas

* **Erdős–Graham (pp.120–122).**  Pack an N × N grid with N = [α − α^{8/11}].  The two border strips have width
  β = α − N = Θ(α^{8/11}) and length γ = Θ(α).  Fill each with parallel stacks R(1, n + 1), n = [β], tilted at a single
  angle θ = Θ(α^{−4/11}) so that each stack touches both long sides.  Each stack wastes two border triangles of area
  Θ(α^{−4/11}), so the border waste is Θ(α^{7/11}).  The end trapezoids T (vertical sides D = Θ(α^{2/11}),
  D′ = Θ(α^{4/11})) are cut into r = Θ(α^{4/11}) sub-trapezoids of base m = [α^{4/11}].  Each gets a grid block plus its
  own tilted stacks, with waste Θ(α^{3/11}) each.  The tilt scale is (strip width)^{−1/2}, not 1/α.
* **Roth–Vaughan.**
  * Fundamental lemma (p.172).  Let S₁ be a square, and S₂ a square or the boundary, at distance ≤ 1 with relative
    inclination θ.  Then some open disc of radius 2 contains S₁ and waste of area ≥ cθ.  Proof: a strip of width
    1/(4m + 1), m ≤ 20, parallel to the line of centres, crosses no square vertex, and some trapezoid in it is pure waste.
  * Main argument (pp.173–185), by contradiction from waste < c₀√(α‖α‖).  From each point Y on the left side, run a path
    to the right, following the nearly horizontal sides of good squares.  Show (Lemmas 2–5) that most such paths are
    "chains" that reach the right side with waste < c‖α‖ along them.  Bad squares are few (Lemma 4, by projections).
    Pick ≫ α chains spaced ≥ 20 apart (Lemma 6 controls how close chains can come).
  * Key step (p.185).  A chain carries k ≥ [α] + 1 squares whose horizontal projections cos θᵢ are disjoint inside
    width α.  So Σ(1 − cos θᵢ) ≥ k − α, which forces max θᵢ ≫ (‖α‖/α)^{1/2}.  Summing the fundamental lemma along the
    chain gives waste ≫ (‖α‖/α)^{1/2} near each chain, and ≫ α chains give ≫ (α‖α‖)^{1/2}.
* **Bui.**  Mark the ⌈x⌉ × ⌈x⌉ grid cells that meet bad squares.  Max-flow/min-cut (with vertex splitting) gives f
  disjoint cell paths from marked cells to the boundary, plus rectangles of total perimeter ≤ 4f covering every marked
  cell.  Along each path, a generalized fundamental lemma (Prop. 21) charges waste ≥ c′c, so W(A) = Ω(f).  Delete
  everything inside the rectangles and regrid; this costs O(perimeter) = O(f).
* **Stromquist.**  Replace unit squares by open (1 + ε)-blocks, so that every boundary case becomes strict.  Use
  nonavoidance lemmas (a block whose centre lies in a region must meet given points or edges).  Show that every block
  must contain one of a finite point set, so more blocks than points is impossible, or the allocation is refuted case
  by case.

## 4. Does any of this matter for exact s(n)?

* **RV gives nothing for side k − ε with n = k² − c [I, now checked against the paper].**  For α = k − ε we have
  ‖α‖ = ε and [α] + 1 = k.  The chain argument forces tilt only if a chain carries k squares.  Its own output is
  waste ≥ 10^{−100}√(kε), which tends to 0 as ε → 0.  Then n ≤ (k − ε)² − 10^{−100}√(kε) excludes nothing at the
  thresholds we care about.  The obstruction is structural, not just the constant: the argument never forces a line to
  meet k squares.  The issue #3 preprint (Ryu 2026, `tasks/issue3-k2-minus-c/review-2026-10-10.md`) adds exactly the
  missing ingredient.  On any line at most k − 1 chords have length ≥ 1 (integrality), and this is combined with a
  sharpened RV lemma (waste ≥ θ(2 − θ)Σλ²/4).  That preprint, not RV, is the RV-style route to c*(k) → ∞.  It is
  unrefereed.
* RV's explicit constant 10^{−100} and Bui's non-explicit c₂ rule out any use at small n.
* **Erdős–Graham's upper bound** gives s(k² − C·k^{7/11}) < k for large k [I; Arslanov et al. 2021 state it].
  Explicit constants are not in the paper.  Making them explicit is still WISHLIST A13's second half.
* **Stromquist's method is the ancestor of our certificates.**  The (1 + ε)-block reformulation and the "every block
  contains one of these points" argument are the classical form of our point-set / piercing certificates.  Memo III is
  also a citable proof that the best 45°-only packing of 11 needs s ≥ 2 + (4/3)√2.  That is relevant whenever we
  compare against restricted-angle families.
* **Bašić–Slivková's Theorem 7** is a uniform formula for piercing sets (§5).  It is a cheap source of starting point
  sets, and a sanity floor, for any n where we build an LP or point certificate.

## 5. Bašić–Slivková: what the n = 61 bound is

* **Citation:** B. Bašić and A. Slivková, "On optimal piercing of a square", *Discrete Applied Mathematics* 247 (2018)
  242–251, doi 10.1016/j.dam.2018.03.048 (Crossref record checked).  No arXiv version was found.
* **Statement (Theorem 10):** s(61) ≥ 7√3/2 + 2√2 − 1 ≈ 7.8906.  It is a **lower bound**.  The paper says that before
  it "no nontrivial lower bound on s(61), that is, better than s(61) ≥ √61 ≈ 7.8102, was known".
* **Method.**  π(U_x) is the piercing number: the fewest points meeting the interior of every unit square inside a
  square of side x.  Proposition 8 says at most π(U_x) unit squares fit.  Theorem 7 bounds π(U_x) using rows of an
  equilateral-triangle lattice:

      π(U_x) ≤ ⌊x⌋(m + 2) [+ ⌊(m + 2)/2⌋ when {x} ≥ ½],  where m = ⌊(2/√3)(x + 1 − 2√2)⌋.

  Just below x = 7√3/2 + 2√2 − 1 this gives π ≤ 7·8 + 4 = 60.  Theorem 9 reproves s(n) for n = 8, 15, 23, 24, 34, 35, 46,
  47, 48 in the same framework.
* **Relation to our result:** **superseded** by s(60) = s(61) = 8.  jlevy's register (entry T-087, secondhand) also
  lists it as superseded, by a floor of 1 + 4√3 ≈ 7.9282 (T-063/T-069, not checked here).
* **Suggested cite:** "The best previously published lower bound was s(61) ≥ 7√3/2 + 2√2 − 1 ≈ 7.8906 (Bašić and
  Slivková, Discrete Appl. Math. 247 (2018), Thm 10)."  If the 7.928 floor turns out to be published, say "best
  published before 2026" with care.

## 6. Claims in our files

Verdicts are against the primary sources above.  Dated notes are not rewritten: their "fix" is for whoever next cites
them, or for a dated addendum.

| File:line | Our claim | Verdict | Fix |
|---|---|---|---|
| `site/www/sources.html:85` (reader-facing) | RV: "the matching lower bound, W(x) ≫ √x at half-integer sides" | **wrong** ("matching") | RV's √x does not match Erdős–Graham's 7/11; it would match only the *conjectured* O(√x).  Use: "the only lower bound: W(x) ≫ √(x‖x‖) when x(x − ⌊x⌋) > 1/6, so W ≫ √x at half-integer sides." |
| `site/www/problems.html:85` (reader-facing) | "Lower bound: W(s) ≥ c·√(s·‖s‖) … (Roth–Vaughan 1978)" | imprecise | Add the hypothesis "for s(s − ⌊s⌋) > 1/6 (c = 10⁻¹⁰⁰ works)".  Without it the statement is false for s extremely close above an integer.  The rest of the item (EG quote, exponent between 1/2 and 3/5) is fine. |
| `site/www/problems.html:158` (reader-facing) | "a small-tilt perturbation (tilts near 2/k, as in Erdős–Graham)" | **wrong** attribution | Erdős–Graham's tilts are Θ(α^{−4/11}) (≈ (strip width)^{−1/2}, p.121), not ≈ 2/k.  Drop "as in Erdős–Graham", or say "(small tilts; Erdős–Graham's are ≈ (strip width)^{−1/2})". |
| `site/www/sources.html:71` (reader-facing) | Memoranda "settled s(6) and s(10) and first argued Gardner's n = 11 conjecture" | correct (incomplete) | Optional: "settled s(6) (hence 7–9) and s(10), proved Gardner's n = 11 conjecture for 45° packings, and sketched s(11) ≥ 2 + 4/√5." |
| `site/www/sources.html:84` (reader-facing) | EG: W = O(x^{7/11}) via tilted strips; opened the asymptotic question | correct | — |
| `site/www/data/timeline.json`, `index.json`, n = 11 (3.8856 entry; prose copied from Ellsworth) | "proved by Walter R. Stromquist in 2002 to be the optimal 45° packing" | imprecise | It was already proved in memo III (15 Nov 1984), Theorem p.7.  Fix in the generator's source or override, not in the JSON. |
| `search/WISHLIST.md:86` (reader-facing via problems.html) | ≥ 5-group packings sit "just below k², where the Erdős–Graham small-tilt packings live" | imprecise | EG's construction is for any large non-integer α (their motivating example is α = N + 1/100, just *above* an integer).  Use "where small-tilt packings live". |
| `search/WISHLIST.md:128` (A13) | Needs "the exact form of Roth–Vaughan Thm 1 (primary source)" | now resolved | Theorem (p.170): α(α − [α]) > 1/6 ⇒ W ≫ √(α‖α‖), explicit c₀ = 10⁻¹⁰⁰ (pp.173, 186).  Keep only the second half of A13 (an effective O(x^{0.6}) constant). |
| `search/WISHLIST.md:17` | "the Roth–Vaughan route to A" is refuted | correct | — |
| `search/WISHLIST.md:33` (N9) | EG: "cannot even rule out W(α) = O(1)" | correct | — |
| `search/FRIEDMAN.md:15` | RV does not give A "[proved, given the bound as secondhand sources state it]" | correct | Drop the qualifier: the bound is now confirmed in the primary source. |
| `search/FRIEDMAN.md:183` | To-do: confirm RV Theorem 1 if ever in hand | resolved | Point to this note. |
| `search/QUADRANT.md:295` | "W(x) ≥ c·√(x(x − ⌊x⌋)), c ≈ 10⁻¹⁰" | wrong (already corrected in place, line 298) | The correction there is right: the factor is ‖x‖ and the constant is 10⁻¹⁰⁰. |
| `site/notes/research-prior-work.md:101` | EG "asked whether W(x)=o(x^{1/2}) possible… (note: Erdős offered $ prizes for the exponent)" | **wrong** | EG say they "cannot even rule out W(α) = O(1)" and "perhaps the correct bound is O(α^{1/2})".  There is no o(x^{1/2}) question and no prize in the paper.  The prizes are Chung–Graham 2020's $100/$250/$500 challenges. |
| `site/notes/research-prior-work.md:102` | RV: W ≫ (x·\|x − round x\|)^{1/2}; W ≫ √x at half-integers; "so exponent lies in [1/2, 3/5]" | imprecise | Add the hypothesis x(x − ⌊x⌋) > 1/6.  3/5 rests on Bui 2508.04603 / McClenagan 2602.01484 (preprints), since the 2020 claim is in error. |
| `site/notes/research-prior-work.md:36` | Memo III = "(Gardner's conjecture for n=11; 15 Nov 1984)" | imprecise | Full title: "III (Cases with n ≤ 65 and Martin Gardner's Conjecture for n = 11)".  Memo II is dated 15 Oct 1984.  II and III are on Daniel H. Wagner Associates letterhead; I is to "450 File (Professional Leave)". |
| `site/notes/research-prior-work.md:89` | "Stromquist 1984 memos: s(6)=3, s(10)=3+1/√2; claimed s(14)=s(15)=4, s(24)=5" | imprecise | Add s(7) = s(8) = s(9) = 3 (memo I, p.2), the 45° theorem for n = 11 and the sketched s(11) ≥ 2 + 4/√5 (memo III, p.10).  For 14, 15, 24 the memo says only "fairly easy to settle by elementary means", with no proof. |
| `site/notes/research-prior-work.md:107` | Bui 2504.09489: only near-axis-aligned squares matter for lower-bound analysis | correct | Could add: W = Θ(W*), non-explicit constants, v1 only (CCCG 2025 submission). |
| `site/notes/research-prior-work.md:137` | Memoranda "not read" | outdated | Read: see §2 here. |
| `tasks/friedman-lit/REPORT.md:8–9` | W(α) ≥ 10⁻¹⁰⁰√(α‖α‖) under s(s − ⌊s⌋) > 1/6 (via DS7) | correct | Now [read]: RV p.170, pp.173/186. |
| `tasks/friedman-lit/REPORT.md:7` | Abstract: side n + ½ ⇒ waste ≫ n^{1/2} | correct | — |
| `tasks/friedman-lit/REPORT.md:11–13` | RV does not give s(k² − c) = k; "check Theorem 1 of the paper if ever in hand" | correct | Checked: the factor is ‖α‖ (nearest integer).  The paper notes that for α(α − [α]) ≤ 1/6 only the trivial bound holds. |
| `tasks/friedman-lit/REPORT.md:14–16` | Method: good/bad squares, fundamental lemma (disc radius 2, waste ≥ cθ), "counting along vertical lines/strips" | imprecise | The main structure is *horizontal chains* from the left side to the right side; the tilt is forced by the projection sum Σ(1 − cos θᵢ) ≥ [α] + 1 − α (p.185).  Vertical lines appear only in Lemmas 4 and 6. |
| `tasks/friedman-lit/REPORT.md:17` | Bui Thm 1: W = Θ(W*) | correct | — |
| `tasks/friedman-lit/REPORT.md:37–40` | EG 7/11; Montgomery (3 − √3)/2 + ε; "Lower: only Roth–Vaughan" | correct | (Bui 2504.09489 is a reduction, not a new lower bound.) |
| `notes/wishlist-literature-2026-10-03.md:81–84` | EG quote p.122; Montgomery exponent | correct | — |
| `notes/wishlist-literature-2026-10-03.md:90` | "only Roth–Vaughan: W(α) ≥ 10⁻¹⁰⁰√(α‖α‖) [S]" | imprecise | Add "when α(α − ⌊α⌋) > 1/6"; now [read]. |
| `notes/wishlist-literature-2026-10-03.md:90–91` | Bui: only good squares matter; W = Θ(W*); "its open problems are about tilt bounds" | correct (roughly) | Its §6 also asks about the radius in Prop. 21 and about algorithms. |
| `notes/wishlist-literature-2026-10-03.md:94–95` | RV is geometric rigidity, not a weighting argument | correct | — |
| `notes/wishlist-literature-2026-10-03.md:115` | EG stacks between parallel lines at distance x waste O(x^{−1/2}) per unit length | correct | Matches p.121: θ = Θ(α^{−4/11}) for β = Θ(α^{8/11}). |
| `notes/wishlist-literature-2026-10-03.md:130–131` | EG: no higher-dimensional statement; open questions about convex curves | correct | — |
| `notes/conjectures-literature-2026-10-04.md:8–11` | EG: one angle θ in the border strips; r = Θ(α^{4/11}) sub-trapezoids; the angle sum is "O(α^{−1/11})" | correct | The paper writes Ω, which it defines as Θ. |
| `notes/conjectures-literature-2026-10-04.md:16` | Bui 2504.09489 says nothing about a bounded number of angles | correct | — |
| `notes/conjectures-literature-2026-10-04.md:44–45` | EG quote on jumps of f(l) from 4k to 4k + 4 | correct | f is the sum of *perimeters* in EG (#106 uses sides; equivalent up to a factor of 4). |
| `notes/conjectures-literature-2026-10-04.md:102–104` | EG ⇒ every fixed-shape family with linear waste is eventually suboptimal; "exactly EG's opening point" | correct | EG p.119: the grid at α = N + 1/100 wastes about α/50. |
| `notes/conjectures-literature-2026-10-04.md:108` | "Erdős–Graham's 'perhaps O(α^{1/2})'" | correct | — |
| `tasks/k2m3-lit/REPORT.md:14` | "n=61 bound ≈ 7.8906 (Bašić–Slivková, unverified)" | correct, now verified | Lower bound: DAM 247 (2018) 242–251, Thm 10; superseded by s(61) = 8 (§5). |

Not checked here: Friedman DS7's own statements, Chung–Graham, Wang–Dong–Li, McClenagan, Bui 2508.04603, Stromquist
2003 (EJC), Nagamochi.

## 7. Provenance notes for the store and the site

* Memo III dates Stromquist's n = 26 packing (s ≈ 5.650629) and his n = 19 packing (s ≈ 4.888) to Nov 1984.  This
  agrees with the site.
* Memo III calls the n = 18 packing 7/2 + ½√7 "New" in Nov 1984.  The site and Stromquist 2003 credit Hämäläinen (1980),
  so Stromquist was unaware of it in 1984.  No change is needed.
* Memo III credits the 3.877 packing of 11 to Gustafsson & Thulin (Gardner, Nov 1980).  The site credits Walter Trump
  (1979), following Ellsworth.  Possibly both are right (independent finds, or first publication versus first find).
  Worth a line in the provenance notes, and not a change without checking Gardner's column.
