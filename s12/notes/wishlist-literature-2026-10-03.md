# Wishlist literature check (2026-10-03)

Labels: [P] read in the primary source, [S] secondhand (abstract, search snippet, another paper's restatement),
[I] my inference, not stated by any source.  "Known" below means a source states it.
Third-party repos were shallow-cloned to the scratchpad for reading only (nothing executed).

---

## 1. Quantitative rigidity of periodic covering measures

**Verdict: the ε = 0 case is known (a special case of the Pompeiu property of the square).  The quantitative bound
|μ̂(ξ)| ≤ 2ε / max_θ |1̂_{S_θ}(ξ)| appears to be new as a stated result, but it is a one-line Fourier estimate.  Present
it as an easy lemma with the Pompeiu citation, not as a theorem.**

* **Pompeiu problem.**  A region D has the Pompeiu property if the only f with ∫_{σD} f = 0 for every rigid motion σ
  is f = 0.  Pompeiu (1929) proved it for the square under a decay hypothesis.  Christov removed the hypothesis and
  added triangles and parallelograms.  Brown–Schreiber–Taylor, "Spectral synthesis and the Pompeiu problem", Ann.
  Inst. Fourier 23(3) (1973) 125–154, https://www.numdam.org/article/AIF_1973__23_3_125_0.pdf [P intro and §6]:
  D has the property iff the Fourier–Laplace transform of 1_D does not vanish on a whole variety
  {z₁²+z₂² = α}.  Every polygonal region has it (Thm 5.9), and so does every convex set with a corner (Cor. 5.12).
  Their §6 notes that for **bounded** f "one needs only to consider the Fourier–Stieltjes transforms of measures in the
  usual sense, and the rotations enter in only a superficial way".
  Later work: Machado–Robins, "The null set of a polytope, and the Pompeiu property for polytopes", arXiv:2104.01957
  [S abstract]; Kolountzakis et al., "Curves in the Fourier zeros of polytopal regions and the Pompeiu problem",
  Anal. Math. (2024), https://eigen-space.org/ps/pompeiu.pdf [S].
* **How the ε = 0 case follows [I].**  Let μ be Λ-periodic with μ(S) ≥ 1 for every unit square S and density 1.  Then
  v ↦ μ(S_θ+v) − 1 is ≥ 0 and has mean 0, so it is 0 a.e.  So ν = μ − Leb is a bounded periodic signed measure that
  integrates to 0 over a.e. congruent copy of the square.  Mollify radially and apply Pompeiu, or note directly that
  ν̂(ξ)·1̂_{S_θ}(ξ) = 0 for every θ, and for each ξ ≠ 0 some θ has 1̂_{S_θ}(ξ) ≠ 0.  The periodic case needs only that
  last fact, which is much weaker than Pompeiu (no complex variety, just the real circle |ξ| = r).
* **Quantitative form.**  No source states a stability or Fourier-decay bound of this kind for covering measures.  The
  only "stability of Pompeiu" literature found concerns perturbing the domain D, not the measure (Ebenfelt, "Some
  results on the Pompeiu problem", Ann. Acad. Sci. Fenn. 18 (1993), https://www.acadsci.fi/mathematica/Vol18/ebenfelt.pdf
  [S]).  The estimate itself is the standard fact that a nonnegative function's Fourier coefficients are bounded by its
  mean, applied to f_θ(v) = μ(S_θ+v) − 1.  The factor 2 depends on normalisation.
* **Nearby (same Fourier mechanism).**  Steinhaus problem: Kolountzakis–Wolff, "On the Steinhaus tiling problem",
  Mathematika 46 (1999), http://www.math.caltech.edu/papers/kowfinal.pdf [S]: if f ∈ L¹ tiles with every rotation of
  Z^d, then f̂ vanishes on every sphere through a nonzero lattice point.  No measurable Steinhaus sets exist for d ≥ 3;
  d = 2 is open.  Kolountzakis–Papadimitrakis, arXiv:math/0009207 [S].  Survey: Kolountzakis, "The study of translational
  tiling with Fourier analysis", arXiv:math/0304005.  Kolountzakis–Lev, "Tiling by translates of a function: results and
  open problems", arXiv:2009.09410.  These concern exact tilings (level sets), not one-sided covering inequalities.
  The step "a covering of density exactly 1 is a tiling" is folklore.
* Suggested framing: "the rigidity at ε = 0 is the (periodic, measure) Pompeiu property of the square [Pompeiu; Christov;
  BST 1973]; the quantitative version is a routine Fourier estimate we have not seen stated."

## 2. Strong duality for the continuum closed-cover LP vs fractional packing; attainment

**Verdict: partially known.  Strong duality and attainment of the packing optimum are stated for the *interior*
model in jlevy/squares (an unrefereed AI-assisted repo).  The *closed*-cover version (Evan's LP) is not stated
anywhere.  Abstract LP duality on compact spaces is standard, and so are counterexamples for general infinite hypergraphs.**

* **jlevy/squares**, `docs/project/research/research-2026-09-10-x027-fractional-duality.md` (2026-09-10, "GPT-6 Astra"
  lanes, "independently reviewed" in-repo; no novelty claimed) [P].  Pose space P_L compact, packing λ ∈ M₊(P_L) with
  **interior** depth ≤ 1 everywhere, cover μ with μ(int S_p) ≥ 1 for every legal pose.  Its "Analytical theorem" (§3):
  ν_∘(L) = τ_∘(L) = τ_fin(L) = τ_ac(L) (finite and absolutely continuous covers give the same value).  The packing
  supremum is attained.  It explicitly does **not** assert an attained covering optimum.  Proof: compactness gives a
  finite hitting set, finite LP duality, and weak-* compactness of packing measures.  Also: ν_ae (closed depth ≤ 1 a.e.)
  = ν_∘.
  The same note gives two **boundary counterexamples** relevant to the closed model.  (a) Closed coverage cannot be paired
  with interior capacities: one atom at a shared edge point covers two touching squares.  (b) The closed-depth feasible
  set of packings is not weakly closed: translates converge to a touching pair with closed depth 2.
* **General theory.**  Aharoni–Holzman (1992), https://holzman.net.technion.ac.il/files/2012/09/gcoptimal.pdf (cited in
  the jlevy note): value equality can fail for arbitrary infinite hypergraphs.  Infinite LP duality on compact spaces
  (Anderson–Nash style) and Sion's minimax are the standard tools.
* **Closed-cover attainment [I, easy; not found stated].**  For closed S, μ ↦ μ(S) is upper semicontinuous under weak
  convergence (Portmanteau: limsup μ_n(F) ≤ μ(F)).  So {μ ≥ 0 : μ(S) ≥ 1 for every closed S ⊂ [0,k]²} is weak-* closed,
  and with bounded mass it is compact.  Hence **the closed-cover LP minimum is attained**, unlike the interior cover.
  The no-gap question then reduces to τ_cl(k) =? τ_∘(k) = ν_∘(k).  τ_cl ≤ τ_∘ is trivial.  A smearing/rescaling argument
  gives τ_∘(L) ≤ lim_{L'↓L} τ_cl(L'), so a gap can only occur where L ↦ τ_cl is discontinuous [I].
* Kearney–Shiu's "duality method" (EJC 9 (2002) #R14) [P] is red/green-lattice duality of unavoidable sets, not LP
  duality.  Massaccesi's "Linear Programing for Square Packing" (blog, 2026-08-21) uses finite LPs only.
* Suggested framing: cite jlevy for the interior-model value equality (noting that it is unrefereed), and prove the
  closed version, plus attainment, ourselves.

## 3. Asymptotics of the LP saving L(k) and of the waste W(s)

**Verdict: the waste bounds are as in our notes.  Nothing is published on any fractional/LP version of the waste or of
the saving L(k).  "L(k) → ∞" (and any rate) appears new.**

* Upper bounds on W(s) (all constructions):
  * Erdős–Graham, JCTA 19 (1975) 119–123: O(s^{7/11}).  Read p.122 [P]: "It is rather annoying that we do not at
    present have any nontrivial lower estimate for W(α).  Indeed we cannot even rule out the possibility that
    W(α) = O(1).  Perhaps the correct bound is O(α^{1/2})."
  * Montgomery (unpublished): (3−√3)/2+ε.
  * Chung–Graham, JCTA 116 (2009): O(s^{(3+√2)/7} log s).
  * Wang–Dong–Li, arXiv:1603.02368: O(s^{5/8}).
  * Chung–Graham, DCG 64 (2020): O(s^{3/5}) claimed.  Error found by Arslanov–Bui, DCG 2025, doi 10.1007/s00454-025-00767-w.
  * **Bui, arXiv:2508.04603** (v1 2025-08-06, v2 2026-03-15) [S abstract]: O(x^{0.6}).
  * **McClenagan, arXiv:2602.01484** (2026-02-01): O(x^{3/5}).
* Lower bound: only Roth–Vaughan, JCTA 24 (1978) 170–186: W(α) ≥ 10⁻¹⁰⁰ √(α‖α‖) [S].  Bui, arXiv:2504.09489, shows that
  only "good" (tilt ≤ 10⁻¹⁰) squares matter, W = Θ(W*) [S].  Its open problems are about tilt bounds; nothing fractional.
* **Fractional/LP version: none found** in any paper, survey (Friedman DS7, Brass–Moser–Pach) or repo (jlevy, wand125,
  tokoharu, chelokot, DRMacIver, Guzhou0806).  The jlevy fractional-duality work is all at n = 11 (fixed small side).
* Consequences [I]: weak duality plus the constructions give L(k) ≤ c*(k)+1 = O(k^{0.6}).  Roth–Vaughan is a geometric
  rigidity argument, not a weighting argument, so it does not obviously give any lower bound on L(k).
* Adjacent covering problem: Soifer (2006) conjectured that a square of side > n cannot be covered by n² + O(1) unit
  squares.  Sriswasdi, arXiv:2609.15876 (2026-09-14) [S]: proves it for n = 4 (17 squares).  Wang–Dong–Li 2016 has a
  covering corollary x² + O(x^{5/8}).  Dósa–Lángi–Tuza, arXiv:2601.16535 (2026): covering, small n.

## 4. Unit squares in a fixed-height strip / a×b rectangle

**Verdict: partially known for rectangles shrinking in both directions.  Nothing found on the fixed-height /
growing-width version or on a "rectangle Friedman" monotonicity (C8).  Erdős–Graham's strip estimate is a device inside
their square construction, not a stated strip theorem.**

* **Nagamochi**, EJC 12 (2005) #R37 [S abstract]: for real a, b ≥ 2, at most ab − (a+1−⌈a⌉) − (b+1−⌈b⌉) unit squares fit
  in any a′×b′ with a′ < a, b′ < b.  For integers that is ab − 2, so ab − 1 squares need a′ = a or b′ = b.
  **Caveat: Karakuş, arXiv:2609.37410** (2026-09-29) [S abstract]: counterexamples to a scoring assertion used in the
  proof.  The published proof of the rectangle bound "is incomplete".  An independent strip-measure proof of a weaker
  rectangle bound recovers s(k²−1) = k but not s(k²−2) = k.  It also gives s(N) ≥ ½ + √(N − ⌊√N⌋ + ¼) for nonsquare N ≥ 8.
* **Arslanov–Mustafin–Shangitbayev**, EJC 28(4) (2021) P4.22 [P]: squeezability δ((Rx,Ry),m), the strip inequality
  δ((Rx,Ry),m) ≤ δ((Rx+1,Ry), m+Ry−1), and explicit squeezable rectangle packings.  Examples: 26 in (4,8) [= 32 − 6];
  64 in (6,12) [= 72 − 8]; 43 in (5,10) [= 50 − 7]; 58 in (6, 11−2/35); and a family 4k²+6k−2 in (2k, 2k+4), waste 2k+2.
  These are *upper* bounds on a rectangle c*(a,b) (both sides shrink).  No rectangle lower bounds beyond Nagamochi.
* **Strips.**  Erdős–Graham's tilted vertical stacks between two parallel lines at distance x waste O(x^{−1/2}) per unit
  length (restated by Bui 2508.04603 and Wikipedia [S]).  Here x is non-integer and large.  Nothing found for an integer
  height h with the width shrinking.  The online "strip packing" literature (e.g. arXiv:1010.4502) is a different
  problem.  "Three squares in a rectangle", arXiv:2608.13595, is tangential.
* Friedman's Packing Center has no squares-in-rectangles page (checked the index).

## 5. Cubes (d = 3)

**Verdict: apparently nothing.  No lower bound for unit cubes in a cube beyond volume was found.  s_3(k³−1) = k
(the analogue of Nagamochi) is not stated anywhere.  No waste asymptotics for cubes were found.**

* Friedman, "Cubes in Cubes", https://erich-friedman.github.io/packing/cubincub/ [S fetch]: best known upper bounds only.
  Non-trivial packings for n = 9–14 and 28–33 (Friedman 1998: 13 → 2.956, 14 → 2.989, 28–33 → 3.707 = 3+1/√2;
  **Haowei Lin, July 2026**: 11 → 2.88295, 12 → 2.93277).  Its text: "For all other values of n < 34, the trivial
  packing is the best known."  No proofs or lower bounds are mentioned.
* Erdős–Graham 1975 (read in full, pp.122–123) [P]: no higher-dimensional statement.  Its open questions are about convex
  curves C and f(C,k).  Roth–Vaughan, Chung–Graham and Bui are all 2-D.  Bui 2504.09489 says only that the method
  generalises to near-rectangular quadrilaterals.
* Searches for "n³−1 unit cubes", cube waste and tilted-cube packing found only unrelated work (cube tilings, Keller,
  online cube packing, perfect packings of harmonic cubes).

## 6. Friedman's Conjecture 1 and c*(k) → ∞: status as of 2026-10-03

**Verdict: Conjecture 1 is open in general.  c*(k) → ∞ is open, and no published source claims it or attacks it.  The
recent progress on specific c (c = 3, c = 4, and k = 8 for c = 5) is our own work plus wand125's.**

* Statement: Friedman DS7 (1998 v1 through the 14 Aug 2009 version) [P, per tasks/friedman-lit].  Evidence there is "true
  of all the best known packings".
* Known before 2026:
  * Nagamochi 2005: c*(k) ≥ 2.
  * Kearney–Shiu 2002: s(6) = 3.
  * Bentz 2010: s(13), s(46).
  * Bentz arXiv:1606.03746: s(22), s(33); "strongly suggest s(m²−3) = m".
  * Arslanov et al. 2021: s(n²−n) < n for n ≥ 12, so c*(k) ≤ k−1.
* **2026 activity** (from jlevy/squares `packing/frontier/RESULTS.md` and wand125/square-packing-bounds README, read [P]):
  * T-064: Evan Daniel, s(k²−3) = k ∀k ≥ 6.  Replayed in full by jlevy (Valid7, 9,800 roots).  wand125 wrote an
    independent checker on 2026-10-02.
  * T-081: Evan Daniel, s(k²−4) = k ∀k ≥ 5 (2026-10-03), resting on ValidTilt9.  jlevy rates it V0: not replayed yet.
  * wand125: s(59) = 8 (the k = 8 case of k²−5) and s(77) = 9 (k²−4 at k = 9, 2026-10-01).
  * jlevy issue #316 is the k²−4 registration request.
* **AI-assisted repos checked** (GitHub, updated ≤ 2026-10-03): jlevy/squares, wand125/square-packing-bounds and -tools,
  tokoharu/square-packing-density-bounds, chelokot/square-packing-archive (Lean-checked archive),
  Queuingtheorydotcom/11SquaresOptimal and 11SquaresFormalized (claimed computer-assisted proof that s(11) = 3.8770835…,
  the Trump value; does not touch s(12)), Guzhou0806/n17-square-packing, DRMacIver/square-packing-research,
  hsthanb4/jsp-000118-all-k (Erdős #106, the sum-of-sides problem f(k²+1) > k; a different problem).
  **None of them states or attacks c*(k) → ∞ or Conjecture 1 in general.**  The only "c*(k)" text in jlevy is its mirror
  of our FRIEDMAN.md.
* arXiv 2026 checked: McClenagan 2602.01484, Bui 2508.04603 v2, Karakuş 2609.37410, Sriswasdi 2609.15876.  None
  addresses Conjecture 1.  The Erdős f(k²+1) = k papers (math/0504341, 2411.07274, 2506.23284) are a different problem;
  don't conflate them.
* Not checked: X/Twitter (not searchable), and Ellsworth's page beyond what tasks/friedman-lit recorded.

## 7. Integrality gap (fractional vs integral unit-square packing in a square)

**Verdict: no published bounds.  Only informal or numerical discussion exists, all from 2026 and all in repos.  Any
proved gap, or a proof that there is no gap at some n, would be new.**

* jlevy/squares fractional-duality note §5–6 [P]: defines the integrality gap m(L) < ν_∘(L).  It has a fractional family
  of mass 11 at L* = 38200/9977 ≈ 3.8288, but s(11) ≥ 3.8264 is known, so no physical gap is established at n = 11.
  There is a gap in a *core* model: 88 cores of side 0.9977 at L = 191/50, excluded by a threshold-charge certificate
  (T-025).  It lists "is there a physical fractional gap for six squares below side three?" as open.  `ideas.md`
  H-034/H-064 (n = 11 fractional piercing and packing floors) are numerical and registered, not resolved.
* DRMacIver/square-packing-research `PUBLICATION-CANDIDATES.md` (audit dated 2026-08-08) [P] claims a "witness-measure
  barrier": no covering-measure argument can show s(17) > 4.5, because the LP optimum at 4.5 "already exceeds 17" (float,
  LP subset).  [I] This is **apparently contradicted** by the later exact certificates s(17) > 4.5058 (Massaccesi) and
  > 4.6130 (Mira-acc/17squares, a measure of mass 16.998 after rescaling cores).  Restricting a valid cover to a smaller
  box keeps it valid, so COVER(4.5) < 17.  Cite it only as a refuted float claim, if at all.
* Our own §9.1.2 notion (k_LP − k₁, with s(12) as the c = 4 candidate) has no literature counterpart.  Bašić–Slivková's
  integral piercing is the precedent for τ vs τ* of unavoidable sets (a different gap; via jlevy ideas.md).
* The CS "integrality gap" literature (stabbing unit squares, arXiv:2106.12385; geometric knapsack) concerns
  axis-parallel combinatorial LPs and is not relevant.

---

### Main URLs
* BST 1973: https://www.numdam.org/article/AIF_1973__23_3_125_0.pdf
* Kolountzakis–Wolff: http://www.math.caltech.edu/papers/kowfinal.pdf
* arXiv:2104.01957 (null set of a polytope / Pompeiu)
* jlevy duality note: https://github.com/jlevy/squares/blob/main/docs/project/research/research-2026-09-10-x027-fractional-duality.md
* jlevy RESULTS: https://github.com/jlevy/squares/blob/main/packing/frontier/RESULTS.md
* Aharoni–Holzman: https://holzman.net.technion.ac.il/files/2012/09/gcoptimal.pdf
* Erdős–Graham 1975 scan: https://static.renyi.hu/~p_erdos/1975-22.pdf
* Arslanov et al. 2021: https://www.combinatorics.org/ojs/index.php/eljc/article/view/v28i4p22
* Karakuş: https://arxiv.org/abs/2609.37410 ; Sriswasdi: https://arxiv.org/abs/2609.15876
* Bui: https://arxiv.org/abs/2508.04603 , https://arxiv.org/abs/2504.09489 ; McClenagan: https://arxiv.org/abs/2602.01484
* Cubes in Cubes: https://erich-friedman.github.io/packing/cubincub/
* wand125: https://github.com/wand125/square-packing-bounds ; DRMacIver: https://github.com/DRMacIver/square-packing-research
