# Wishlist: proof pieces, conjectures, and partial results (2026-10-03)

Reader-facing version (object-level problems and conjectures only): the site's `problems.html` ("Open problems").

A working list of statements whose proof (or refutation) would move this project forward, for us and for anyone who
wants to help.  Each item: statement · status · what it unlocks · source.  Labels as elsewhere in the repo:
**[proved]** (exact or Lean), **[measured]** (float LP / scan), **[heuristic]**, **open**.  **(Lean-ready)**: a full
paper proof exists and no new kernel machinery is needed.  Paths are relative to `s12/`.

Notation: `c*(k) = max{c : s(k²−c) = k}`.  `L(k)` = best saving `k² − (closed-cover LP value on [0,k]²)`; `L(k) > c ⇒
s(k²−c) = k`, so `c*(k) ≥ ⌈L(k)⌉ − 1`.  `S_ins(k)` = best saving of an insertable certificate (FRIEDMAN §2.1).
`M(w)` = seam capacity of a σ = 0 band of width w (SEAM_W1).  `D(w)` = per-corner saving of a fixed-profile family.

Sections: **N** new candidates (literature-checked 2026-10-03, citations in `notes/wishlist-literature-2026-10-03.md`), **P** packing targets (2026-10-04), then **A–F** items already stated in our notes (A
Friedman structure, B seam/wall/corner, C threshold certificates, D n = 12, E n = 11 and 17, F Lean), **G** algebraic
degree of s(n) (notes only), then picks.
Superseded or refuted items are omitted on purpose: n12-gap §7's list, the Roth–Vaughan route to A, the slide
conjecture (n = 89), corner-local savings (`D*(R) = 0`), S2 at c = 3, 4 (`S2_INSERTABLE.md`), exact SOS certificates for n = 12 (`SOS_PROBE.md`, n12-gap §4.12), band-LP values as caps
on M(w), clique certificates below 12, the s(12) ≥ 3.98 consolation prize.

## N. New candidates

| # | Statement | Status / literature | Would unlock | Source |
|---|---|---|---|---|
| N1 | **Quantitative rigidity.** Λ-periodic μ on a region where every unit square at every angle slides through a full period, density `1 + ε`, every square of mass ≥ 1 ⇒ `|μ̂(ξ)| ≤ 2ε / max_θ |1̂_{S_θ}(ξ)|` for ξ ∈ Λ* \ 0; at ε = 0, μ is Lebesgue. | [proved] (sketch: `f_θ = μ * 1_{S_θ} ≥ 1` with mean `1+ε` has `‖f_θ − mean‖₁ ≤ 2ε`).  ε = 0 is a case of the **Pompeiu property** of the square (Pompeiu 1929; Brown–Schreiber–Taylor, Ann. Inst. Fourier 1973).  No quantitative/stability version for covering measures found. | a price for every non-uniform periodic structure; explains why thin (high-frequency) seams are cheap; input to N5, N6, B8 | S2_INSERTABLE §2; this list |
| N2 | **Axis blindness.** `1̂` of an axis-parallel unit square vanishes at every ξ ∈ ℤ² \ 0, so axis-only covers carry free period-1 structure (grid lines, `{1,2,3}²`); rotation removes it. | [proved] (one line). | the conceptual statement "rotation is the content" of n = 12 (add to n12-gap §1) | — |
| N3 | **One-direction rigidity for bands.**  x-Fourier modes `μ̂_n(y)` of an x-periodic σ = 0 band, tested by height kernels `K_{θ,n}`; generalises the row lemma. | open. | rigorous caps on M(w) (only `w − ¼` proved); pairs with B8 | FRIEDMAN §8 |
| N4 | **Other dimensions and shapes.**  The insertion lemma (slab width 1 + max projection) and crossing-cube rigidity hold for unit cubes in ℝ^d and any convex body allowed to rotate.  Cube analogue of F_c: `s₃(k³ − c) = k` for all large k. | proofs carry over [proved, routine].  Literature: **nothing** on lower bounds, `k³ − 1` analogues or waste for cubes (Friedman's "Cubes in Cubes" lists packings only). | a new problem family where our methods apply; 3D checkers are far costlier, so start with LP go/no-go at small k | — |
| N5 | **Interior-Lebesgue conjecture.**  ∃ w₀: for all k, ε-optimal closed covers of `[0,k]²` can be taken Lebesgue on `[w₀, k − w₀]²` with loss → 0. | [measured] hints (S60 §4 ring table: interior within 0.1–0.4 of Lebesgue).  **Cheap test**: LP at k = 5..8 with interior forced Lebesgue beyond w; loss ε(k, w) vs w (~10³–10⁴ CPU-s). | ⇒ S1 (A7) with δ = wall-crossing cost ⇒ almost-monotone fractional Friedman | S60_COVER §4 |
| N6 | **Interior near-uniformity, quantitative.**  Interior Fourier coefficients of optimal covers → 0 (rate in k), via N1 on local averages. | open; data cheap (Fourier-analyse the k = 5..8 LP covers). | step toward N5 | — |
| N7 | **Price formula for insertability.**  price(k) ≈ wall saving lost on the four crossing segments of length 1 + √2. | [heuristic]; naive count `4 · 2.414 · 0.27 ≈ 2.6` vs measured 1.87 / 1.73: partial. | predicts `k_ins(c)` (C4) | S2_INSERTABLE §5 |
| N8 | **Strong duality, closed model.**  No gap between the closed-cover LP and ν_f on `[0,k]²`; attainment. | **Partly known**: jlevy/squares (2026-09-10 note, unrefereed) proves no gap for the *interior* model with packing attainment.  Closed model: attainment of the cover optimum is easy (μ ↦ μ(S) upper semicontinuous for closed S); a gap can only sit at discontinuities of `L ↦ τ_cl(L)`.  Aharoni–Holzman 1992: gaps happen for general infinite hypergraphs. | makes L(k) a clean object both ways; licenses dual arguments (FRIEDMAN §7.1); Lean candidate later | — |
| N9 | **Growth of L(k).**  Is `L(k) → ∞`?  `L(k) = Θ(k^β)`? | open; **no LP/fractional version in the literature.**  Note `L(k) → ∞ ⇒ c*(k) → ∞` (A), so it is a route to A, not a weaker statement.  Upper `L ≤ c* + 1 = O(k^{0.6})` (waste 3/5: Bui arXiv:2508.04603, McClenagan arXiv:2602.01484; Chung–Graham 2020 erroneous).  Erdős–Graham 1975: "cannot even rule out W(α) = O(1)". | ties R1–R3 to one exponent | FRIEDMAN §0 |
| N10 | **LP completeness / integrality gap ≤ 1.**  `⌈L(k)⌉ − 1 = c*(k)` for k ≥ k₀. | open; **no published integrality-gap bounds** (jlevy defines the gap; a core-model gap only).  s(12) (k = 4) is the candidate exception. | says the dual route loses nothing asymptotically | — |
| N11 | **Monotonicity / Fekete for seam capacity.**  `M(w+1) ≥ M(w)`, or subadditivity giving `lim M(w)/w^β`. | open (not obvious: shifting a band up exposes its floor). | makes R2 (B2) a statement about one exponent | — |
| N12 | **One-direction families: fixed-height rectangles.**  For fixed height b, one x-periodic slab (period mass b) propagates a certificate on `[0,a]×[0,b]` to every a′ ≥ a: "ab − c unit squares need width a at height b" for all a ≥ a₀ from one box.  No crossing square, so N1 does not force Lebesgue: **points + segments certificates qualify**, and Lean checks those today. | open; **literature: nothing** on fixed height with growing width, nor on rectangle monotonicity (C8 = A12).  Nagamochi's rectangle bound `ab − (a+1−⌈a⌉) − (b+1−⌈b⌉)` has a proof gap (Karakuş, arXiv:2609.37410); Arslanov–Mustafin–Shangitbayev 2021 give squeezable rectangle packings (upper side). | first Friedman-type "all sizes" theorems that are **hypothesis-free in Lean**; half of C8 for small b | FRIEDMAN §2.1 |
| N13 | **Finite ladders.**  A one-step insertable cover (q > 1) proves its value on a window k..k+J, J ≈ slack/(q − 1). | [proved] bookkeeping (S2 §4).  Example: one k = 4 box gives `s(k²−3) = k` for k = 4..8 [measured]. | cheap base-case generator for per-k results | S2_INSERTABLE §4 |
| N14 | **Cheap s(12) gain: re-solve our 1,736 s(12) points with near-tight-cell rows.**  jlevy got `s(12) ≥ 3.9702002` this way (ours 3.9686155), with LP slack left. | [reported] on jlevy; not replayed here. | better Lean s(12) bound with existing machinery; the technique (cert-slack lever 3) generalises | notes/jlevy-s17-techniques.md §2 |
| N15 | **Rule atoms as LP columns** (capacity-one k-of-m and winning-subset rules, as in the s(17) floors R068/R071). | [measured] elsewhere: points-only stalled at ≈ 4.614 for s(17); rules reached 4.66044.  Lean: an instance of our clique lemma. | strict bounds (s(17), s(11)); probably not zero-margin exact proofs | notes/jlevy-s17-techniques.md §2 |
| N16 | **Non-integer plateaus.**  Least n with `s(n) = s(n+1) ∉ ℤ` (equivalently `s(n) = s(n+1) ≠ ⌈√n⌉`).  Provable frontier: the least n not yet ruled out. | **Conjecturally 147**: Friedman's table gives n = 147, 232, 264, 290, 295 as the (n+1) packing minus a square, at sides 7+4√2, 8+11√2/2, 9+11√2/2, 14+5√2/2, 17+√2/2 (site data).  **Provably the frontier is n = 18**: strictness proved for every non-integer pair below ((4,5), (5,6), (9,10), (10,11), (11,12), (16,17), (17,18)) and for (19,20) (wand125 s(20) ≥ 1959/400 > 3+4√2/3); (12,13), (20,21) are integer-plateau questions (s(12), c = 5).  (18,19) is open: s(18) ≤ (7+√7)/2 vs s(19) ≥ 1927/400. | a clean public frontier; each strict pair is a lower-bound milestone.  Weaker forms (2026-10-07, problems.html §3, Lean `Conjectures` with implications proved): finitely many ⇒ bounded length, ⇒ finitely many per fractional part; per-fract + finitely many plateau fracts ⇒ finitely many; all-n per-fract ⇒ plateau per-fract (not ⇒ finitely many) (next: s(19) > (7+√7)/2, offered on jlevy#281) | CEILINGS_17_20 §6 |
| N17 | **Tilings are never optimal.**  `s(n) ∉ ℤ ⇒ s(k²n) < k·s(n)` for all k ≥ 2 (Evan, 2026-10-07).  Tiling a k×k array of copies (or scaling by k and splitting) gives `≤`; the claim is strictness. | **Stated here**; literature check 2026-10-07 (DS7 2009, Ellsworth's catalogue, Kearney–Shiu, Erdős–Graham, Chung–Graham, arXiv 2508.04603 / 2602.01484) found nothing.  **[proved] for k ≥ K(n)**: tiling wastes `k²(s² − n)`, optimal side-ks packings waste `O(k^{3/5})` (7/11 suffices), giving `s(k²n) ≤ ks − kδ/(4s)`, δ = s² − n.  k = 2 open.  Table check: all 28 non-integer records with 4n ≤ 323 beat doubling by ≥ 0.41; k = 3, 4 margins ≥ 1.04, 1.59. | packing targets P6; a structural statement ("optimal packings are not self-similar") | this list; problems.html §4 |

---

## P. Packing targets: where conjectures predict a better packing  (2026-10-04)

The other sections ask for proofs.  This one asks for **packings**: n where a conjecture or a trend in the best-known
table predicts that the current record can be beaten.  A new packing at any of these n is a record, and for a family
member it also removes n from the family.  Data: `search/record_structure.py` over the site's best-known table
(Ellsworth's, as mirrored by the site build: complete for n ≤ 323, plus a few larger n).  These are best-known values,
not optima.  **Before searching, check jlevy's register for current values** (TODO rule).

| # | Target | Why | What it would show |
|---|---|---|---|
| P1 | **s(90) < 10** | The k²−k family at k = 10.  In the margin table below, the c = 10 column reads .153 (k = 8), .056 (k = 9); a straight extension crosses 0 at k ≈ 10.  The c = 11 column ends at .0032 (Cantrell's s(110)). | The last k with s(k²−k) = k is ≤ 9 (now: in [3, 10]; true at 2, 3; false at 11 and ≥ 12).  s(72) < 9 looks less likely (c = 9 column: .054 at k = 8). |
| P2 | **Frontier borderlines:** s(183) < 14 (c = 13), s(242) < 16 (c = 14), s(274) < 17 (c = 15), s(308) < 18 (c = 16) | Each column extends to ≈ 0 at the first "= k" entry.  The columns slow down (c = 16: .065, .041, .024, .012), so straight extension overstates the chance. | Lowers the record upper bound on c\*(k) at that k |
| P3 | **The √7 family**, side k + (√7−1)/2, angles 0° and 24.295°: records at n = 18, 53, 86, 127, 151, 176, 204, 234, 299 (k = 4, 7, 9, 11–15, 17).  Best bets: the largest k (299, 234, 204). | Waste grows linearly (s² − n: 5.3 at k = 4, 18.7 at k = 17), so O(s^0.6) packings eventually beat every member.  The family already has no record at k = 5, 6, 8, 10, 16.  None is proved optimal (s(18) ≥ 4.695, wand125, vs 4.8229).  At n = 53 and 151 packings with 4 and 3 rotation groups tie the family value. | Removes a member; evidence on when fixed-shape families die |
| P4 | **Other fixed-fractional-part families** (45°-based): 0.7071 = √2/2 (Göbel strips, 16 records with n ≤ 296; per Ellsworth best known for every a < 44 except a = 3, first beaten at a = 44, n = 2043); 0.5355 (11 records, n = 65 … 291); 0.6569 (9, n = 66 … 294); 0.7782 (8, n = 150 … 298); 0.8467 (4: 54, 107, 178, 267). | Same linear-waste argument as P3.  The largest-k members go first. | As P3 |
| P5 | **Conjectured non-integer plateaus** n = 147, 232, 264, 290, 295 (N16) | Each is drawn as the (n+1) packing minus one square, and **each lies in a P4 family**: 147 at 0.6569, 232 and 264 at 0.7782, 290 at 0.5355, 295 at 0.7071. | A better packing at n removes that plateau; one at 147 moves the conjectured first plateau |
| P6 | **Doubled Arslanov–Mustafin–Shangitbayev packings:** s(964) < 31.98161, s(1092) < 33.97641, s(1228) < 35.96544 (2× the records for 241, 273, 307) | N17 predicts each 2×2 tiling is beatable.  Neither Ellsworth's catalogue nor jlevy's register lists these n (checked 2026-10-07), so the tilings are the best known, and they beat the grid by 0.018–0.035.  Look at the seams: wall-held squares now press against the neighbouring copy. | Evidence for N17 at k = 2 (and records at those n) |

Margin table, `k − s(k²−c)` for the best-known packings ("= k": best known is the grid; blank: not in the k-interval):

| k \ c | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 |
|---|---|---|---|---|---|---|---|---|
| 8 | .0542 | .1533 | .1771 | .2929 | .2992 | .4286 | | |
| 9 | = k | .0559 | .1183 | .1728 | .2012 | .2929 | .3431 | .4645 |
| 10 | | = k | .0503 | .1118 | .1612 | .1771 | .2574 | .2929 |
| 11 | | = k | .0032 | .0503 | .0741 | .1533 | .1770 | .1924 |
| 12 | | | = k | .0086 | .0435 | .0888 | .1187 | .1749 |
| 13 | | | | = k | .0179 | .0416 | .0683 | .1183 |
| 14 | | | | | = k | .0258 | .0431 | .0649 |
| 15 | | | | | | = k | .0259 | .0414 |
| 16 | | | | | | = k | .0092 | .0244 |
| 17 | | | | | | | = k | .0118 |
| 18 | | | | | | | | = k |

Record upper bounds on c\*(k), k = 4..18: 4, 5, 6, 7, 8, 9, 10, 10, 11, 12, 13, 14, 14, 15, 16.  They are
nondecreasing, so Friedman's C1 (A1) points at no specific n today; it would if a new record broke the pattern.

**Structure of the records** (context for the questions below).  "Rotation groups" = distinct angles mod 90°.
* s(11) (proved) uses **one** tilt angle (49.818°) plus axis squares; its minimal polynomial is an irreducible
  octic (constructibility not decided: Galois group not computed).  s(28) (best known, Ellsworth 2025, rigid) uses one
  tilt angle (59.709°); its polynomial is an irreducible sextic with Galois group **S6** (sympy), so the best-known
  s(28) is not expressible in radicals.
* Packings with ≥ 5 rotation groups sit mostly at fractional parts ≳ 0.87, just below k², where the Erdős–Graham
  small-tilt packings live.  Exceptions: n = 68, 103, 123, 297.  Packings with ≤ 3 groups hold about half of each
  k-interval through k = 17 (9 of 18 at k = 17), with no visible decline.

**Questions these targets bear on** (object-order; not yet stated elsewhere in this list):
* Sharpened ascent `c*(k+1) ≤ c*(k) + 1` (proved: + 2, FRIEDMAN §1.2).  It would make {k : s(k²−k) = k} an
  initial segment; today nothing rules out the family resuming after a failure (s(12) < 4 would not decide s(20)).
* Is each fractional part `{s(n)}` attained by only finitely many n?  (Fixed-part families with linear waste must
  die; sublinear waste does not obviously forbid a repeat.)
* Does the degree of s(n) over ℚ tend to ∞ along non-integer optima?  Is any optimal side a non-integer rational?
  Rational *records* exist: s(50) ≤ 7+4/7 (3-4-5 tilt), s(230) ≤ 15+28/41, s(261) ≤ 16+28/41, s(293) ≤ 17+26/41
  (20-21-29), per `notes/conjectures-literature-2026-10-04.md` §2.
* `W_m(s)`: least waste using ≤ m distinct angles.  W₁ = Θ(s) (axis only).  Is W_m = Θ(s) for each fixed m (then
  optimal packings need ever more angles)?  Checked 2026-10-04: Erdős–Graham and every later construction vary the tilt per stack; nothing published on bounded angle counts.
* Best waste coefficient κ (waste ≈ κ s) of a fixed-shape family at a given fractional part.  Göbel strips: κ ≈ 1;
  the grid at fractional part 0.707: κ ≈ 1.41.  An explicit family with κ < 1 beating the strips from modest n (below
  the observed first fall, n = 2043)?  Are the strips optimal among two-angle families (ties to W_2)?  Existence of
  some explicit beating family is folklore (Erdős–Graham with explicit constants, A13); the crossover is bookkeeping.
* Does the many-angle regime spread down from the top of each k-interval as k grows?  Related: are optimal packings
  eventually asymmetric?
* (Evan, 10-04) How many locally optimal (rigid) packings exist at n = 90, and how are their sides distributed?  Guess: exponentially
  many jammed at exactly 10 (any full axis line of 10 pins the side), sub-10 basins rare; that would make "hard for annealing"
  precise.  Measurable with `search/packer` (quench, squeeze, dedupe by side + contact graph).  Cf. jlevy H-012.

---

## A. Friedman / c*(k) structure  (13 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| A1 | **C1 (Friedman's Conjecture 1):** `c*(k+1) ≥ c*(k)` for all k. | open.  Upper row of known c* is nondecreasing (atlas records).  Translation descent is refuted (n = 89 is short by 1–2 even with 25 translation classes or diamond moves). | the full conjecture | FRIEDMAN §0, §3, §7 C1 |
| A2 | **Block-descent structure lemma** (the replacement for the slide conjecture): every tilted region of a near-optimal packing can be shrunk in its own frame, compatibly with its neighbours. | open, no candidate statement yet.  Descent works on every atlas record with n ≤ 130 except n = 89 (7×7 diamond at 45°); the good descents re-pack the block (7×7 → 6×6 plus axis squares). | the primal route to C1 and S5 | FRIEDMAN §3, §4.1, §6 |
| A3 | **C7 / "A":** `c*(k) → ∞` (for every c, `s(k²−c) = k` for all large k). | open in the literature (friedman-lit REPORT).  Known: `3 ≤ c*(k)` for k ≥ 3, `c* ≤ k−1` for k ≥ 12, `c* = O(k^0.6)`. | first unbounded lower result on c* | FRIEDMAN §0, §7 C7 |
| A4 | **R1:** `c*(k) ≥ k^β` for some β > 0 (β ≤ 0.6 forced). | open.  Follows from `D(w) ≳ w^β`, i.e. from R2 ∧ R3 (B2, B6). | quantitative A | FRIEDMAN §9.2 |
| A5 | **C2 (fractional Friedman):** `L(k+1) ≥ L(k)`. | open.  Data `L(4) ≤ 4, L(5) ≈ 4.25, L(7) ≈ 5.21, L(8) ≈ 5.60` is consistent [measured].  The four-sub-box averaging identity is proved and reduces C2 to a wall-layer bound (A6). | each LP certificate gives `A_c` from its k; "no gaps" in LP-provable values | FRIEDMAN §7, §7.1, §9.4 S3 |
| A6 | **Wall-layer lemma for C2:** some optimal fractional packing at k+1 has a boundary frame (averaged over a shift law π on [0,1]²) cut by weight ≤ 2k+1. | open.  The naive bound fails: wall layers can carry ≈ √2 per unit length.  The dual side meets the same gap (repairing a smoothed cover costs ≥ √2 per unit wall naively). | C2 (A5) | FRIEDMAN §7.1 |
| A7 | **S1 (bounded splice cost):** there is δ < ∞ such that every valid μ on `[0,k]²` (k ≥ k₀) has an insertable valid μ′ with `μ′ ≤ μ + 2δ`.  This gives `L(k′) ≥ L(k) − 2δ` for k′ ≥ k. | open.  The measured insertability price is 1.87 (k = 4) and ≈ 1.73 (k = 5), i.e. δ ≈ 0.9 per direction [measured, float]. | almost-monotone fractional Friedman | FRIEDMAN §9.4; S2_INSERTABLE §5 |
| A8 | **Splice-cost LP / δ(k):** the exact cost of one wall crossing in isolation (an arbitrary wall on each side, a periodic patch of length ≈ 3.4), and whether the price stays ≈ 1.8 at k = 6, 7. | open.  It is the next measurement, not yet run (~1–2·10⁴ CPU-s). | decides S1's δ and whether S2-type shortcuts reopen for c = 5 | FRIEDMAN §9.6 item 2; S2_INSERTABLE §5; TODO |
| A9 | **C4:** the insertable class is asymptotically complete: `L(k) − S_ins(k)` → 0 or stays bounded. | open.  `lim S_ins = 4 sup D ≤ 4 M(∞)`.  If C5 holds (M(∞) < ∞) and L grows, C4 fails. | C7 via families | FRIEDMAN §7 C4 |
| A10 | **S4:** bounded family delay, `k_fam(c) ≤ 2 k_LP(c)`. | open.  Measured delay ≈ 2, 1.6, 1.3–1.7 for c = 3, 4, 5 [heuristic].  False if M saturates. | quantitative bridge from isolated certificates to families | FRIEDMAN §9.3, §9.4 |
| A11 | **S5:** Friedman up to bounded loss, `c*(k′) ≥ c*(k) − O(1)` for k′ ≥ k. | open.  Needs control of the integrality gap, so LP tools alone cannot reach it. | weak C1 | FRIEDMAN §9.4 |
| A12 | **C8 (rectangle Friedman):** `c*(a,b)` is nondecreasing in a and in b separately. | open.  It is the half-step of C1; each step is one periodic-slab insertion. | C1 one direction at a time | FRIEDMAN §7 C8 |
| A13 | **Literature checks behind the "≳ c^{5/3}" floor:** the exact form of Roth–Vaughan Thm 1 (primary source), and whether the O(x^0.6) waste constant is effective. | secondhand only.  The claim that no threshold is linear in c rests on the second check. | soundness of FRIEDMAN §0 and §9.0 as published | FRIEDMAN §6, §9.0 |

Proved here and reusable (not wishlist items, listed so they are not re-proved): the ascent `c*(k+1) ≤ c*(k)+2`
(§1.2); the insertion lemma (§2.1); C3 (`S_ins` is nondecreasing); the averaging lemma (§9.1.1, in Lean); the
four-sub-box identity (§7.1); the crossing-square Lemma 1 (S2_INSERTABLE §2).

## B. Seam, wall and corner  (13 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| B1 | **`M(1) = 3/4` fully proved:** profile A (seam `{n}×[½,1]` at density 3/2 plus the line y = 1 at density ¼) is valid, i.e. `F(θ, y0) ≥ 0` on `[0°,45°] × [0,1]`. | the upper bound ¾ is [proved].  The hand cases θ = 0, 45° are [proved].  The general (θ, y0) case is [measured]: min F = −3·10⁻¹⁶ on a 4501×10001 grid, with tight set = 3 families.  Missing: piecewise calculus with breakpoints `h ∈ {sc, s, c, H−sc, ½+sc}`, or 2-variable interval arithmetic. | first exact seam capacity; Lean target after B12 | SEAM_W1 §2.3, §6.1; TODO Next 2 |
| B2 | **C5 / R2 (seam capacity):** is `M(∞) < ∞`?  Conjecturally `M(w) ≳ w^{1/3}`. | open.  Lower bounds M ≥ 1.433, ≈1.755, 1.965, 2.116, 2.231 (w = 2..6) are [measured] and are not caps.  The only proved cap is `w − ¼`.  The heuristic `M(w) ≳ C w^{1/3}` comes from a budget of ½ and a price of `e ≈ 0.1 g^{1.5}` per row, but is partly superseded by the drift correction (B5). | if ∞: families can reach every c (with R3).  If finite: the dual route caps at c < 4M(∞) | FRIEDMAN §7 C5, §8, §9.2; SEAM_W1 §0 |
| B3 | **§8 missing lemma (b): wall row.** A valid wall row that spends the excess budget (≤ ½) while feeding the bulk seam profile. | open ("none looks hard"). | with (c), (d), B5: `M(∞) = ∞` | FRIEDMAN §8 |
| B4 | **§8 missing lemmas (c) and (d).** (c) Slow variation: a profile whose amplitude g(y) varies on scale ≫ 1 costs `O(|g′|)` extra excess (total O(max g)).  (d) The top transition to Lebesgue. | open ("none looks hard").  Item (a), exact 1D price certificates, is **done**: κ ≤ 0.0221 / 0.0315 / 0.0397 / 0.0486 at g ≈ 0.046 / 0.080 / 0.127 / 0.206 [proved, exact B&B]. | `M(∞) = ∞` | FRIEDMAN §8; SEAM_1D §4 |
| B5 | **Phase-constrained seam price:** in a real band the excess is a drift (≤ r per residue class, row lemma), not a uniform e.  Is the seam price still sublinear when the seam must sit where squares see excess? | open.  The next model is the 1D problem with a height phase.  The w = 6 optimum thins the seam just above integers [measured]. | rescues the w^{1/3} heuristic, i.e. B2 | FRIEDMAN §8 "Correction" |
| B6 | **R3 (corner costs a bounded fraction):** `D(w) ≥ κ M(w)`. | open.  Measured D/M ≈ 0.67, 0.66, 0.66, 0.64, 0.61 for w = 1..5 [measured; M is a lower bound for w ≥ 2]. | R2 ∧ R3 ⇒ R1 ⇒ A | FRIEDMAN §9.2–9.3 |
| B7 | **Corner cost at w = 1 is exactly ¼** (D = ½ against M = ¾): explain it with a few squares, as §1 of SEAM_W1 explains M(1) ≤ ¾ with one. | [measured] D = 0.5004 (R = 2).  Unexplained. | first exact case of R3; template for corner duals | SEAM_W1 §0, §6.2; FRIEDMAN §9.6 item 4; TODO |
| B8 | **Rigorous upper bounds on M(2), M(3)** from a dual-side band LP (square weights as variables, pointwise load rows from the arrangement), checked by the unfilled-area duality. | open.  The duality theorem is [proved] (SEAM_W1 §4).  The saved w = 2 dual is a lattice artifact (load 182 against Λ = 92). | first caps below `w − ¼`; calibrates B2 | SEAM_W1 §4.1, §6.3; FRIEDMAN §9.6 item 3 |
| B9 | **Is the apparent saturation of M(w) a pitch artifact?**  A fine-pitch band primal at w = 4–6 (efficient small-amplitude seam needs pitch ≤ 0.005). | open (measurement). | B2 evidence | SEAM_W1 §6.4; FRIEDMAN §9.3 |
| B10 | **Two-zone ansatz at w = 2:** the x-uniform + seam class loses ≈ 9 % at w = 2 (1.30 vs 1.433), so an analytic w = 2 profile needs structure across the period. | [measured]; no analytic profile yet. | an analytic witness for M(2), towards B2 | SEAM_W1 §5 |
| B11 | **Slope-ceiling bound:** σ = 0 forces every band row to be tight at θ = 0, so the first-order tilt margin satisfies κ*(w) ≲ 1/(2(w−δ)).  Measured κ* ≈ w^{−1.3} at pitch 0.1. | [heuristic] Fubini estimate, right order at w = 2, 3; the data decay faster. | the cost scaling of certificates for larger c | K2M4_MARGIN §2–3 |
| B12 | **Seam bound** `E(n) ≤ 0` and `D ≤ m_v ≤ w` for valid quadrant families. | [proved] on paper (dilated grid).  **(Lean-ready)**: nearly a corollary of `packing_le_measure`. | Lean base for B1 and families | QUADRANT §1.4; FRIEDMAN §9.6 item 5; TODO Next 1 |
| B13 | **Localisation lemma for a window checker:** a family is valid if the corner window plus one band period is (far-band poses reduce mod 1, the interior is Lebesgue). | stated as needed, not written (variant of QUADRANT §1.3 / §8.1). | much cheaper certificates for k²−5 and beyond | FRIEDMAN §5 golf 1 |

## C. Threshold certificates (`s(k²−c) = k`)  (9 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| C1 | **k²−5 family:** `s(k²−5) = k` for all k ≥ 2R+2, from one box (R = w = 5, box 13). | float GO: D ≈ 1.302 > 5/4 at κ = 0.02, margin ≈ 0.05 (still drifting 10⁻⁴/round) [measured].  w = 4 caps just above 1.25 (no-go). | C6 at c = 5 (new theorem); a step towards F₅ | K2M4_MARGIN §3; FRIEDMAN §4.2; TODO Major |
| C2 | **F₅ base cases:** `s(k²−5) = k` for k = 6, 7 (s(31), s(44)) and from k = 9 up to the family start.  The threshold `k₁(5)` itself involves s(20) = 5. | open.  The shipped s(32) cover saves only 4.29 at k = 6 (< 5); k = 8 (s(59)) is proved. | F₅ | FRIEDMAN §0, §4.2, §9.3 |
| C3 | **C6 via insertable boxes:** `S_ins(8 or 9) > 5` would give c = 5 from one box. | likely dead.  The insertability price ≈ 1.8 against L(8) ≈ 5.6 [measured] gives S_ins(8) ≈ 3.8. | C6 at c = 5 | FRIEDMAN §7 C6; S2_INSERTABLE §5 |
| C4 | **k_ins(4) ∈ [6, 8]:** pin it down; the prediction is ≈ 8–9. | [measured] no-go at k = 5 (price ≈ 1.73). | S1/S2 evidence | FRIEDMAN §9.1.3; S2_INSERTABLE §5 |
| C5 | **One-step ladder certificate:** the q > 1 insertable k = 4 cover proves `s(k²−3) = k` for k = 4..8 from one box (smallest margin 0.33). | [measured, float]; not certified.  Of methodological interest only (F₃ is a theorem). | a case-free "one box, five values" certificate | S2_INSERTABLE §4 |
| C6 | **Second-order / near-tangent primitive:** certify `μ = 1 + 0·θ + cθ²` at the σ = 0 tight rows (no margin), plus exact germ-family lemmas in place of deep subdivision. | open.  Lemma Z covers θ = 0 and Lemma E covers zero margin with first-order growth. | uniformity in w; cuts checker depth (k²−4 ran 815k CPU-s, 53 % in u ∈ [0, 1/16]) | K2M4_MARGIN §3; FRIEDMAN §5 golf 3 |
| C7 | **Independent verification of `ValidTilt9`** (F₄ for k ≥ 8 rests on one implementation, qx2_zm + Lemma E) and of `ValidTilt7` (k²−3; zmx2 refuses polygons). | open.  The reviews found no break (k2m4-review: 0 gaps, 2 minor).  Lemma E has not been audited by a second person. | removes the single-implementation caveat on F₃'s tail and F₄ | TODO Correctness; QUADRANT_EXACT §7; FRIEDMAN §6 |
| C8 | **Second reader for the 4 mutation-blind refinements** in the zmx2 area-density second implementation of k2m3. | open. | completes k2m3's second implementation | TODO; ZMX2_AREA §13 |
| C9 | **F_c decidability:** each `F_c` is decidable (Tarski) once `A_c` holds with an effective k₀. | [proved] remark.  Practical form: k₀(c) explicit for c ≤ 4 only. | — (context for C1, C2) | FRIEDMAN §1.1 |

## D. n = 12  (12 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| D1 | **s(12) = 4.**  Equivalently, no 12 closed unit squares in `[0,4]²` are pairwise disjoint as closed sets. | open.  Certified `s(12) ≥ 15680/3951 = 3.968616` (Lean).  Every degree-1 family is ≥ 12 at t = 4 [proved lower bounds]. | decides whether 11 or 12 is the largest n with s(n) < 4 | n12-gap §1–2 |
| D2 | **13-point pure unavoidable set for `[0,4]²`** (Kearney–Shiu slack-1 accounting at n = 12). | open.  `≥ 13` is [proved] (ν_f ≥ 12.2688); `≤ 14` is Friedman's set.  No family with h ≥ 14 has been found up to 1,761 squares.  The 2,809-square union family is undecided (HiGHS 2 h, no verdict).  It is a pure integrality question. | the K–S route to D1 (the rest of that argument is also unwritten) | n12-gap §6.1; unavoid13-no §0, §7 |
| D3 | **Leaf A admits no closed packing:** four corner squares holding `{a_i, b_i}` plus eight singleton squares.  Margin 0 on a plateau up to 32°; no sub-core. | open.  Every degree-1 family pins it at `[11.999999926, 12]` [proved].  Planned route: replacement hulls (step 1), then a combinatorial chain proof (step 2). | kills one leaf of the corner tree; template for D1 | tasks/leafa-proof; n12-gap §2.2 item 6; RANK8 |
| D4 | **Existence clauses at T = 4:** (i) the H1 conjecture; (ii) every `δ* = 0` configuration carries a chain of four among its near-axis squares. | [measured] 727/727 and 352/352 sampled optima.  None proved, and obstruction X says no O(t)-stable mechanism can prove them. | the angle-space route | n12-gap §4.10, §6.2; counting-ladder §2.6 |
| D5 | **Obstruction X (second-order rigidity):** packings vs configurations of margin −0.36t² near the tiling stratum.  At T = 3: (E3-row) on the "exactly 3 or 4 axis-parallel" strata, and non-existence on the coherent cone. | open.  It is a diagnosis: each named mechanism provably fails.  Rule: no s(12) route that cannot do s(6) = 3. | an s(6) = 3 proof in the row language, then T = 4 | n12-gap §4.11, §6.2; S6_LOCAL §5 |
| D6 | **Lemma CC as a DAG path:** can the order chain be taken as a DAG path at ε > 0?  And DAG acyclicity at large tilt, beyond `ε < (1+δ)/(T−1−2δ)`. | open. | stronger chain counting | n12-gap §6.3; counting-ladder §1.4 |
| D7 | **Far field at T = 4** (≤ 3 near-axis squares). | open, with no candidate mechanism.  The centre pigeonhole misses by 0.5 % (`4−√2 = 2.5858` vs `1/d₁₂ = 2.5725`) and the capped LP does not close. | the other half of the angle-space route | n12-gap §2.2 item 10, §6.4 |
| D8 | **V(4, ε, 3) for ε ≥ 5°; ε⁴ cells at T = 5; c_T beyond T = 6; the T = 5..10 cases.** | open (pool-limited).  Known: `c_T = 1, 3/4, 1/3, 3/8` at T = 2..5 [proved]. | ladder data | n12-gap §6.5 |
| D9 | **Are Bentz's 16 points a closed cover of `[0,4]²`?** | evidence only (margin 0 at the tiling, 2.2·10⁻⁴ in the 44° band).  Not load-bearing. | — | n12-gap §6.6 |
| D10 | **Value of `ν_f^closed(4)`** (bracket `[12.2688, 12.956]`, believed ≈ 12.4). | [proved] lower bound; the ~12.4 is [heuristic]. | calibrates every LP route at t = 4 | n12-gap §2.2 item 2; COVER4 |
| D11 | **`s(n²−n) = n` for n = 5..10** (the conjecture that survives the counterexamples at n = 11 and n ≥ 12). | open (literature).  The repo verifies the counterexamples exactly. | context: sanity filter for D1 mechanisms | n12-gap §5 |
| D12 | **Lean-ready n = 12 pieces:** axis-parallel ⇒ ≤ 9 squares (`{1,2,3}²` is a closed cover for axis squares); Lemma CC (chain counting, order form). | [proved] on paper.  **(Lean-ready)**.  See F6. | — | n12-gap §2.2 items 1, 5; TODO Lean |

## E. n = 11 and n = 17  (7 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| E1 | **Continuum version of the n = 11 SA-2 certificate:** on the 28-square extremal support SA-2 = SA-3 = α = 10 [proved, exact], via two single-variable splits.  Open: does the rank-10 statement survive on refined pose sets? | open (fixed support only). | an elementary rehearsal for leaf A (D3); a route to s(11) independent of the case split | N11_LADDER §6, §8 |
| E2 | **QSTAB = 32/3 in the continuum** at n = 11 for `3.83375 ≤ t ≤ 3.87`. | [proved] on four fixed poses; continuum [heuristic]. | locates the clique ceiling at n = 11 | N11_LADDER §2, §7 |
| E3 | **s(11) is settled** (Queuingtheorydotcom's computer-assisted proof, accepted on jlevy/squares as T-060).  Optional only: an independent check with our checkers, or finishing its Lean (6 `sorry`). | done elsewhere; optional cross-check | an independent confirmation | N11_LADDER §1; jlevy README |
| E4 | **s(17): local isolation of Bidwell's 13-square core** in `[0,T]²` (first-order rigidity of the tied rows after dropping the movable squares 4, 5, 10, 12). | open; task on hold (courtesy to the s(11) authors, and CPU).  Site analysis [measured, float]: squares 5 and 12 move alone, 4 and 10 move with them. | the endgame of an s(11)-style proof of s(17) | tasks/s17-core-isolation; TODO s(17) |
| E5 | **s(17) global sizing:** the minimum number of closed cells of diameter < 1/(U−1) covering [0,1]² (guess 26–30, vs 16 at n = 11), and whether the case analysis is feasible. | open (estimate only). | feasibility of E4's global part | tasks/s17-core-isolation §5 |
| E6 | **Pure-cover ceiling at t = 4.6604:** is `ν_f ≥ 17`?  If yes, pure covers are at their limit for s(17). | open; nobody has computed it (~1 CPU-h). | decides whether to stop pure-cover work on s(17) | TODO s(17) |
| E7 | **Green's `s(17) ≥ (40√2+19)/17`:** the set drawn in DS7 Fig. 34 is not unavoidable (counterexample exact in Q(√2)).  Is there another 16-point set or a sharper triangle lemma that supports the bound? | "not supported".  We proved `s(17) ≥ 4.44` from the corrected set H* (critical at `2√2+√65/5`).  The four-row topology cannot reach Green's value. | historical correction only (preprint bounds are larger) | GREEN_FIG34 §0, §3 |

## F. Lean  (11 items)

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| F1 | **`ValidTilt7` / `ValidTilt9`, part 1: the `CovT` tree layer** (X/Y/U/XM/YM/UM splits, clip and SYM nodes, leaf dispatch, generator from the leaf dump). | not started.  ~400 lines, mechanical.  k2m4 needs a self-contained leaf dump (the review found the clip_bin slabs unrecorded). | prerequisite for F2, F3 | lean-leb-mass §7.1; k2m4-review/leaves |
| F2 | **Part 2: PIECE with the polygon:** Lemma S(b), the area of `U ∩ K` for K an intersection of certified half-planes (a convex-polygon clip and its soundness), plus a bridge from CovM segment mass. | not started.  Est. 1–2k lines.  Closes ~17k of V3's 32k leaves. | most of `ValidTilt7` | lean-leb-mass §7.2 |
| F3 | **Part 3: Lemma E in Lean** (the E′/E″ concave minorants of `λ(Q∩U)`, the line-arrangement concavity argument, vertex certificates incl. general-degree Bernstein and the S-procedure). | not started.  Est. 3–6k lines, weeks.  A per-leaf certificate format (the kernel only checks) is advised. | **hypothesis-free** `s(k²−3) = k` (k ≥ 6) and `s(k²−4) = k` (k ≥ 8) in Lean | lean-leb-mass §7.3; LADDER |
| F4 | **Remaining segment lemmas:** Lemma L′/V, SPLIT (Lemma R), wall corners in Corollary L (the last two need a general-degree Bernstein lemma). | not done. | discharging `S21CheckerCover` / s(45) hypotheses without workarounds | LADDER "Segments"; lean-segments §5 |
| F5 | **s(21) = 5 hypothesis-free** (`S21CheckerCover` in the kernel, cand A, ~51 CPU-h). | open; it is now an engineering plus kernel-cost item. | s(21) = 5 fully in Lean | TODO Lean; lean-s21 |
| F6 | **n = 12 batch:** ~~state s(12) = 4 as an open statement~~ (done 2026-10-07: `Conjectures.S12` in `lean/Sqpack/Conjectures.lean`, with the other problems-page items); axis-parallel ≤ 9; Lemma CC. | **(Lean-ready)** paper proofs exist (D12). | — | TODO Lean; n12-gap §2.2 |
| F7 | **s(6) = 3 (Kearney–Shiu 2002) in Lean**, the last missing k of k²−3. | not done (check external repos first). | complete Lean coverage of `s(k²−3) = k` for k ≥ 3 (modulo F1–F3) | TODO Lean |
| F8 | **Seam bound** `E(n) ≤ 0`, `D ≤ m_v`; then **`M(1) ≤ ¾`** (one 45° square, triangle area, δ → 0). | **(Lean-ready)** (B12 and SEAM_W1 §1). | first Lean results on the Friedman side | FRIEDMAN §9.6 item 5; TODO Next 1 |
| F9 | **Validity of the w = 1 profile** (M(1) ≥ ¾), via the reusable slice lemma `N(y) ≥ 1[W(y) ≥ 1]`. | blocked on B1 (Cavalieri plus trig casework). | `M(1) = ¾` in Lean | FRIEDMAN §9.6 item 5 |
| F10 | **Insertion lemma in Lean:** parametrise `bentz_of_valid7` / `BentzFam` over c and the box predicate. | **(Lean-ready)**: the paper proof is short and `BentzFam` is already generic in R. | certificate-level Friedman (C3) in Lean | FRIEDMAN §2.1, §6 |
| F11 | **Ascent `c*(k+1) ≤ c*(k) + 2`**, the budget lemma, the row lemma (FRIEDMAN §8). | **Ascent [proved] in Lean 2026-10-07** (`cStarStepsOfTwo`, `lean/Sqpack/ConjecturesProofs.lean`; also `tilingBound`, s(k²n) ≤ k·s(n)).  Budget and row lemmas: **(Lean-ready)** elementary, with packing-side and measure-side proofs on paper.  Low payoff each. | groundwork for C2/C5 statements in Lean | FRIEDMAN §1.2, §8 |

## G. Algebraic degree of s(n)  (4 items; notes only, not for `problems.html`: Evan 2026-10-04, "neat but extremely niche")

Every s(n) is algebraic (routine: Tarski–Seidenberg, angles as (c, s) with c² + s² = 1).  Data: Ellsworth's
exact forms (`notes/conjectures-literature-2026-10-04.md` §2).  **s(83) ≤ 9.634757648631…** has degree **672**: we factored the
`Root[…, 27]` polynomial from `square-83.svg` with python-flint (2026-10-04).  It is irreducible over ℚ, has coefficients up to
724 digits and 52 real roots, and s is the 27th.  Structure (from the SVG's Mathematica model): 55 axis-parallel squares, and 28
tilted squares in four groups at three angles (a, b ≈ 43°, c ≈ 36°).  There are 3 contact equations in (s, a, b, c), so the
contacts leave a one-parameter motion.  The fourth equation is that s is stationary along it,
`det ∂(f1, f2, f3)/∂(a, b, c) = 0`, and that condition is where the degree comes from.  Rigid records stay small (s(11), s(37):
degree 8).  Next largest on the page: 198 (s(206)), 158 (s(179)), 144 (s(108)), 83 (s(235)).

| # | Statement | Status / evidence | Would unlock | Source |
|---|---|---|---|---|
| G1 | **Is deg s(n) unbounded?** | open.  Every proved non-integer value has degree ≤ 2 until s(11) (below).  Best-known values reach 672 (n = 83), but those are not optima. | — | this list |
| G2 | **Least n with s(n) provably of degree ≥ 3.** | **n = 11 once E3 is confirmed**: s(11) is Trump's octic root (degree 8).  It is least because s(n) is known for every n ≤ 10, and those values all lie in ℚ(√2) (s(5) = 2+½√2, s(10) = 3+½√2, the rest integers).  So s(11) is also the first s(n) outside ℚ(√2). | a one-line corollary of T-060 | E3; README |
| G3 | **Growth of the maximum degree** of a locally optimal packing of n squares (or of s(n)). | Upper bound **2^{O(n)}**: routine via quantifier-elimination / critical-point degree bounds (Basu–Pollack–Roy; cf. Nie–Ranestad, algebraic degree of polynomial optimization).  Probably really 2^{O(#tilted squares)}.  Not stated for packings anywhere we found.  Lower bound for local optima: probably easy by chaining mechanisms like s(83)'s; for s(n) itself it is G1. | — | conjectures-literature §2 |
| G4 | **Degree from the contact structure:** bound deg s in terms of the number of distinct angles and how far the contacts are from rigid (number of stationarity conditions). | open; data only (Ellsworth's table).  This is the same pipeline as TODO "pose → certified local optimum → closed-form s". | predicts which records hide large degrees (37 entries are "not yet analytically optimized") | TODO Packings |

---

## Picks (2026-10-03): payoff for the difficulty

1. **N12, fixed-height rectangles.**  A new theorem class nobody has looked at, provable with point/segment
   certificates and checkable in Lean today.  First step: LP go/no-go at b = 3, 4 (cheap).
2. **B1, `M(1) = 3/4` fully proved.**  A two-variable inequality with known breakpoints, numerically tight to 10⁻¹⁶.
3. **N1 + N3 → B8, rigidity write-up and rigorous M(w) caps.**  Turns seam capacity from one-sided data into brackets.
4. **N5, interior-Lebesgue LP test.**  ~10³–10⁴ CPU-s; tells us whether S1 (bounded splice cost) is the right target.
5. **F8 + F10 (Lean seam bound, `M(1) ≤ ¾`, parametrised insertion lemma)** and **N14** (better Lean s(12) bound).
6. **B7, the corner costs exactly ¼ at w = 1.**  First exact case of R3.

Hard, high payoff: B3–B5 (`M(∞) = ∞`), N9 (`L(k) → ∞`), D2 (13-point unavoidable set), D3 (leaf A), F1–F3
(hypothesis-free k²−3 / k²−4 in Lean), C1 (k²−5 family, after cost cuts).
