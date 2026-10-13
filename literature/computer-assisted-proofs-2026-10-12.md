# Computer-assisted optimality proofs for packing problems: closest prior work on method (2026-10-12)

Our earlier scans covered the history of squares in squares and numerical optimizers.  This one covers the
closest prior work on **method**: computer-assisted (and formally verified) optimality proofs for packing
problems, read against our certificates (weighted closed covers of `[0, m]²` whose validity, "every closed unit
square captures weight ≥ 1", is checked by exhaustive exact subdivision of the 3-dimensional pose space of one
square, by two or three exact checkers, and for s(13), s(32) and s(k² − 3) also inside Lean's kernel).

Labels: **[read]**: full-text section read this session; **[verified]**: bibliographic record, abstract or
README fetched this session; **[unverified]**: from memory or inference; **[S]**: secondhand, reported by
another source read this session (named).  PDFs were downloaded to the session scratchpad only; no code was run.

---

## 1. Interval branch and bound for circles in a square (Markót, Csendes; predecessors)

**Problem.** Pack n equal circles in the unit square, maximizing the radius; equivalently, place n points in
`[0,1]²` maximizing the minimum pairwise (squared) distance ("point packing") [read, Markót 2021 §1].

**Lineage** (as summarized by Markót 2021 and Markót–Csendes 2006, both [read]):
- n ≤ 20: computer-aided proofs by de Groot, Peikert, Würtz (n = 10, 1990) and Peikert et al. (1992), in
  floating point [S, Markót 2021].
- 21 ≤ n ≤ 27: Nurmela & Östergård, DCG 22 (1999) 439–457, floating point, active regions approximated by
  polygons; Markót calls it "the most effective approximation of the non-interval methods" [S, Markót 2021;
  record verified; full text not reached, Springer paywall in practice].
- Markót–Csendes first built a validated (interval) method that re-verified the known optima for n ≤ 27 "with
  some exceptional problem cases (where the number of subintervals turned out to be prohibiting)" [read,
  Markót–Csendes 2006 §1].
- **n = 28, 29, 30**: Markót 2004 (Numer. Algorithms 37, n = 28) and Markót & Csendes, SIAM J. Optim. 16 (2005)
  193–219, the first fully interval-arithmetic proofs [verified records; method read in the 2006 companion].
- **n = 31, 32, 33**: Markót, J. Global Optim. 81 (2021) 773–803, open access [read].

**Method** [read, Markót 2021; Markót–Csendes 2006 preprint]:
- *Search space*: the full 2n-dimensional configuration space (all n point coordinates).  Interval B&B with
  directed outward rounding (C-XSC; the 2005 work also PROFIL/BIAS).
- *Global phase, "tiling"*: the square is cut into a k × l grid of tiles (7 × 6 for n = 28–30; 6 × 8 for
  n = 31–33, scaled to `[0,24]²` so tile corners are integers); each subproblem fixes which tiles hold a point.
  The raw counts are 48-choose-n ≈ 4.24·10¹², 2.25·10¹², 1.09·10¹² for n = 31, 32, 33.  They are never enumerated:
  half-square (6 × 4) tile sets are solved first and combined, with density filters ("no optimal packing of
  31–33 points has ≥ 10 points in a 0.5 × 0.5 region", from the known f₁₀*), pattern lookup in stored
  lower-dimensional results, and "cutoff monotonicity" (regions eliminated at a lower cutoff stay eliminated at
  a higher one, so n = 33 is solved first and reused for 32 and 31).
- *Symmetry*: only the rectangle symmetries (two reflections, 180° rotation): a tile combination is kept only if
  its 48-bit bitset is lexicographically smallest among its 4 images.  Markót notes this does not remove all
  square symmetries.
- *Key accelerating test, "method of active areas"*: each point's feasible region is an interval polygon; a part
  of a region within the cutoff distance of **all vertices** of another point's region is deleted (Lemma 1: a
  point within d of all vertices of a convex polygon is within d of all of it).  Correctness of the polygon
  update (Algorithm 2) is proved by hand case analysis; the 2006 paper exists precisely because "a detailed
  algorithmic description of this tool together with a proof of correctness is essentially required in order to
  allow checking for correctness and reproducing the results".
- *Local phase*: continued B&B to ε = 10⁻¹², clustering, round-robin bisection; "free points" (points that can
  move without affecting optimality, i.e. rattlers) are detected and not bisected.  No interval Newton/Krawczyk
  uniqueness step in the 2021 paper.
- *Output*: an **enclosure** of the optimum, e.g. f₃₂* ∈ [0.2131745625898764, 0.2131745625898779] (width
  ≈ 10⁻¹⁵), and enclosures of all optimal configurations; the best-known packings lie inside.  Not an exact
  algebraic value.

**Size**: n = 28, 29, 30: about 53, 50, 21 h on a Pentium IV, ~10⁶ stored subintervals [read, 2006].  n = 31, 32, 33:
26, 61, 13 CPU-h on one core of a 1.73 GHz laptop; 3.8·10⁷, 8.5·10⁷, 1.6·10⁷ tile combinations survive the filters
and go to B&B [read, 2021 Tables 1–3].  The 2021 method is "roughly about 40–100 times faster" than the 2005 one.

**Verification and publication**: peer-reviewed journal papers; single implementation; correctness argued by
proofs of the lemmas and of the polygon algorithm.  No code or data release is mentioned in the 2021 paper ("Data
sharing not applicable") [read].  No independent re-implementation and no formal proof found [verified: none
found in searches; absence not proved].  The 2007 Springer book *New Approaches to Circle Packing in a Square*
(Szabó, Markót, Csendes, Specht, Casado, García) is advertised "with program codes" [verified record].

## 2. Montanher, Neumaier, Markót, Domes, Schichl: unit squares in a circle (JOGO 2019)

"Rigorous packing of unit squares into a circle", J. Global Optim. 73 (2019) 547–565 [read: Optimization Online
preprint, 2018, authors Montanher–Neumaier–Markót; Domes and Schichl added in the published version, record verified].

- *Problem*: smallest circle holding n rotatable unit squares.  Proved: **n = 3** (n = 1, 2 were known),
  r₃ ∈ [1.2884705074, 1.2884705080], all optimal arrangements enclosed in 8 boxes of width ≤ 5.7·10⁻⁸.
- *Model*: a constraint satisfaction problem over (r, c_i, θ_i), the full 3n-dimensional configuration space.
  Containment = four vertex constraints (convexity of the disk).  Non-overlap via **sentinels**: two squares are
  disjoint iff neither contains any of the other's 9 sentinel points (4 vertices, 4 edge midpoints, centre),
  i.e. 9 non-smooth max-abs constraints per ordered pair.
- *Solver*: GloptLab (Java), interval B&B with forward–backward constraint propagation on interval unions, linear
  relaxation filtering, and feasibility verification (Ipopt point, then a small box checked rigorously);
  depth-first, bisect the widest coordinate.
- *Cutting the space*: centres' domain `[-r, r]²` is cut into 4p² isosceles triangles (each grid cell split by its
  diagonals) small enough that no triangle holds two centres (centres of disjoint unit squares are ≥ 1 apart);
  subproblems = choice of one triangle per square.  Iterate on the number of squares: infeasible 1- and
  2-square subproblems prune 3-square ones.  Symmetric triangle combinations (8 symmetries of the grid) are
  discarded without processing.  Rounding: the grid is **scaled to integer vertices** and the square side becomes
  a tiny interval s, so the tiling constraints are exact.
- *Size*: 36 triangles; 6 + 43 + 12 subproblems out of C(36,3) = 7140; phase 3 about 104 h, of which 90 h
  (≈ 4.5 days) on the one subproblem containing the optimum (i7 laptop).
- *Verification/publication*: peer-reviewed; AMPL-style models, statistics and logs published as supplementary
  zip [read: URL given in paper; not fetched].  Single solver; no formal proof.
- *Their own framing* (2018): for squares in a square, "none of these contributions rely on computer-assisted
  proofs" [read].

## 3. Squares in squares: earlier proofs (hand) and 2026 computer-assisted / formal work

- **The resource ("unavoidable points") method is classical and hand-checked.** Göbel 1979 (unit points),
  El Moumni (line segments: mass = intersection length), Nagamochi 2005 ("a combined system of (weighted)
  points, line segments, and a rectangular area"), Stromquist 2003 (nonavoidance lemmas: a box with centre in a
  region must meet certain points or boundary parts), Bentz 2010/2016 (continuously varying families of
  unavoidable point sets) [read: Bentz 2016 §2, which summarizes this lineage; Stromquist 2003 §1].  Stromquist's
  s(10) proof (EJC 10 (2003) #R8, 10 pages) shows no sign of computer use: the only "computations" are
  closed-form checks [read; grep of the full text].  Nagamochi 2005 likewise [read: no mention of computation].
  Our certificates are exactly this kind of argument (weighted points, segments, area), with the validity of the
  resource system checked by machine instead of by lemmas.
- **s(11) = 3.8770835… (Trump's 1979 packing optimal), Queuingtheorydotcom, 2026** [verified: READMEs of
  `github.com/Queuingtheorydotcom/11SquaresOptimal` and `11SquaresFormalized`, fetched 10-12].  Method (as the
  README states it): a "closed center cover" reduces any smaller packing to one of 2,184 canonical cell patterns;
  exact certificates (rational coordinates, angle intervals, algebraic data, dual weights) exclude 2,180; a
  square-symmetry argument reduces the last 4 to one case, closed by a complete branch cover that either
  contradicts or forces the packing into a certified isolation neighbourhood of Trump's packing; side
  `T = (6u+4)/(1+2u−u²)`, u a root of a degree-8 polynomial.  Verification: a 23-stage verifier (`VERIFY.py`);
  the Lean repository accepts 7,920 modules with zero `sorry`, but "selected expensive, exact numerical
  certificate checks use `native_decide`", so "the final theorem trusts Lean's kernel and native compiler; this
  is not a kernel-only verification claim".  It incorporates upstream certificate checkers under `Sqpack/` and
  `Interop/Wand125/` (our Lean namespace and wand125's work; exact extent not checked).  Unrefereed.  This is the
  only other computer-assisted *global optimum with a non-trivial (tilted) packing* for squares in squares that we
  know of.
- **chelokot, square-packing-archive** (Lean; kernel proof of s(6) = 3 and a counterexample to Nagamochi's
  Lemma 1 with a replacement proof of s(k² − 2) = k) and **jlevy/squares** (register; replays of our checkers)
  are already recorded in `certificates/k2m3/README.md` and `notes/literature-s32.md`; not re-read here.
- No 2020–2026 peer-reviewed computer-assisted optimality proof for squares in a square turned up [verified:
  searches only].

## 4. Fekete et al.: squares into a disk, worst-case density (SoCG 2021)

Fekete, Gurunathan, Juneja, Keldenich, Kleist, Scheffer, "Packing squares into a disk with optimal worst-case
density", SoCG 2021 (LIPIcs 189, 36), arXiv:2103.07258 [read: §1 related work, §3].
- *Problem*: critical density for packing arbitrary sets of squares into a disk: 8/(5π) ≈ 0.509.  Not equal
  squares, not a fixed n; the proof analyses a packing algorithm.
- *Method*: manual analysis plus interval-arithmetic lemmas: subdivide a k-dimensional box (k ≤ 9) into cuboids
  "in a tree-like fashion according to a fixed ordering of the k variables"; if the first j variables already
  exclude a region, the later ones are not split; adaptive refinement.  They chose this "simple technique ...
  arguably less error-prone than relying on deeper, more advanced methods".
- *Size*: all automatic proofs < 1.5 h on 4 laptop cores.
- *Verification*: code public (`github.com/phillip-keldenich/squares-in-disk`), peer-reviewed conference paper;
  no formal proof.  Their related-work paragraph also points to Bánhelyi et al. (B&B for covering polygons by
  disks) [read, not followed].

## 5. Flyspeck: formal proof of the Kepler conjecture (reference point)

Hales 2005 (Annals 162, 1065–1185) and Hales et al., "A formal proof of the Kepler conjecture", Forum Math. Pi 5
(2017) e2, arXiv:1501.02155 [read: §1–§5].
- *Structure*: the main statement is reduced, in HOL Light, to three computational claims: ~1,000 nonlinear
  inequalities, the classification of tame plane graphs (enumeration, formalized in Isabelle by Nipkow and Bauer),
  and a large family of linear programs.
- *Nonlinear inequalities*: interval arithmetic with **Taylor interval approximations**, plus monotonicity: a
  partial derivative of fixed sign reduces a k-dimensional box to a (k−1)-dimensional face.  The domains were
  partitioned into > 23,000 inequalities; "an informal search procedure ... finds an appropriate subdivision ...
  and other information", then "a formal verification procedure ... takes the result of the informal search
  procedure as input and produces a final HOL Light theorem".  All formal arithmetic is on HOL natural numerals
  (no hardware floats).  ~5,000 processor-hours on Azure (32 cores, < 1 week), rechecked at Radboud.
- *Sharp inequalities* ("f ≤ 0 with f(x₀) = 0", x₀ always a corner of the domain): (1) f(x₀) = 0 by exact
  arithmetic; (2) f < 0 off a small neighbourhood U of x₀ by Taylor intervals; (3) every partial derivative on U
  has the sign that puts the maximum at x₀ [read].  This is the closest precedent to our zero-margin leaves.
- *Main statement* checks in ~5 h (40 min from recorded proofs); peer-reviewed; the formal proof is the
  verification.

## 6. Other computer-assisted packing proofs, close in method

- **Hales & Kusner, densest packing of regular pentagons** (arXiv:1602.07220 v2, 2016) [read: abstract, §12]:
  reduces the pentagonal ice-ray conjecture to area minimization over clusters of ≤ 4 Delaunay triangles,
  solved by computer: interval arithmetic (OCaml, Alliot–Gotteland), automatic differentiation (1-jets) for local
  optimality near the optimum, a meet-in-the-middle search, cube subdivision.  ~5,000 lines OCaml on GitHub, ~60 h
  on a quad-core.  Journal publication status [unverified].  Congruent polygons, but a density problem in the
  plane, not a finite container.
- **Musin & Tarasov, Tammes problem N = 13 (DCG 48, 2012) and N = 14 (Exp. Math. 24, 2015)** [read: N = 14 abstract
  and §3 sketch]: enumerate candidate contact graphs with `plantri` (~1.5·10⁹ graphs for N = 14), eliminate them
  by linear approximations and LP, then handle the survivors with interval bounds.  Combinatorial-type enumeration
  rather than a geometric cover; not formalized.
- **Subercaseaux, Nawrocki, Gallicchio, Codel, Carneiro, Heule, "Formal verification of the empty hexagon number"**
  (ITP 2024, LIPIcs 309, 35; arXiv:2403.17370) [verified: abstract/records]: Lean verification of a SAT-based
  geometry result (h(6) = 30); verified encoding + checked unsatisfiability proof; found and fixed a flaw in the
  original symmetry-breaking argument.  Not packing, but the closest recent example of a large computation made
  kernel-trustworthy in Lean, and of formalization catching a symmetry-reduction bug.

---

## Comparison

| Work | Problem | n | Method | Space cut | Size | Checked by | Formal |
|---|---|---|---|---|---|---|---|
| Nurmela–Östergård 1999 | circles in square | 21–27 | polygon active regions, floating point [S] | 2n-dim configs | ? | one code | no |
| Markót–Csendes 2005; Markót 2004 | circles in square | 28–30 | interval B&B (C-XSC, PROFIL/BIAS), active areas, tiling | 2n-dim configs, 7×6 tiles | 53 / 50 / 21 h | one code; correctness proofs on paper | no |
| Markót 2021 | circles in square | 31–33 | interval B&B, interval polygons, tile patterns | 2n-dim, 6×8 tiles, ~10¹² combos → ~10⁷ | 26 / 61 / 13 CPU-h | one code, not released | no |
| Montanher et al. 2019 | unit squares in circle | 3 | interval CSP (GloptLab), sentinels, triangle tiling | 3n-dim configs | ~104 h (one case 90 h) | one solver; logs published | no |
| Fekete et al. 2021 | squares into disk, worst-case density | all sets | manual + interval subdivision lemmas | ≤ 9-dim parameter boxes | < 1.5 h | code public | no |
| Hales–Kusner 2016 | regular pentagons, plane density | — | interval + autodiff + meet-in-middle | Delaunay-cluster params | ~60 h | code public | no |
| Musin–Tarasov 2012/15 | Tammes on sphere | 13, 14 | contact-graph enumeration + LP + intervals | combinatorial types | ~1.5·10⁹ graphs | one code | no |
| Flyspeck 2014/17 | Kepler (spheres) | — | Taylor intervals, tame graphs, LPs | > 23,000 inequality boxes | ~5,000 CPU-h | the formal proof | **HOL Light + Isabelle, kernel** |
| Queuingtheorydotcom 2026 | squares in square | 11 | cell-pattern case split + exact certificates + branch cover | 2,184 patterns, then poses | ? | `VERIFY.py`; Lean | Lean, **with `native_decide`** |
| **ours** | squares in square | 13, 21, 32, 45, 60, 61; k²−3 (k ≥ 6); k²−4 (k ≥ 5) | weighted point/segment/area cover (LP-found) + exact zero-margin subdivision | **3-dim pose space of one square**, D4 fold or full | s32: 2.8 CPU-h (`zeromargin.py`), 11.6k CPU-s (`zmx2 --full`); k2m4 ~226 CPU-h | 2–3 separately written exact checkers, mutation suites, external replay | Lean kernel, no `native_decide`: s(13), s(32), k²−3 (wand125's 9,800-module build); k²−4 and s(21) conditional |

---

## How to position our results

- **The certificate is classical; the scale and the checking are new for this problem.**  Weighted unavoidable
  sets of points, segments and area (Göbel, El Moumni, Nagamochi, Stromquist, Bentz) are the standard lower-bound
  tool for squares in a square.  What we add is finding such measures by LP / column generation (thousands of
  atoms) and verifying them by exhaustive exact computation, not by hand lemmas.  We should say this plainly and
  cite the lineage, not present the cover idea as ours.
- **Structurally different from the interval B&B line.**  Markót and Montanher et al. search the 2n- or
  3n-dimensional configuration space of all n objects and must fight its symmetries (tilings, 10¹² tile
  combinations, subproblem trees).  Our check lives in the 3-dimensional pose space of a single square, and n
  enters only through the total weight, so the cost scales with the container, not with n.  That is why n = 60, 61
  are reachable and why a periodic family plus one finite check (`Valid7`) gives **all k** at once; we found no
  computer-assisted packing proof covering an infinite family this way [inference from the works read; not an
  exhaustive search].
- **Exact versus enclosure.**  The interval works prove an enclosure of the optimum (width ~10⁻¹⁵) and of the
  optimal configurations; our statements are exact equalities (s(32) = 6), because the bound is attained with zero
  margin and the checkers use exact integer/rational primitives for the tight cases.  Handling sharp cases exactly
  is not unprecedented: Flyspeck proves its sharp inequalities by an exact value at a corner plus derivative signs
  nearby.
- **But the method only reaches the "easy" optima.**  All our exact values are cases where the trivial grid is
  optimal (s(n) = ⌈√n⌉).  For non-trivial optima the cover method gives lower bounds that do not meet the upper
  bound (our s(11) ≥ 3.8143 versus 3.8771); the 2026 s(11) proof needed a configuration case split (closer to the
  Markót/Montanher style).  We should not suggest the method settles tilted optima.
- **Formal verification: in the tradition, with a stronger trust base than the nearest neighbour.**  Flyspeck is
  the precedent for kernel-checked packing computations (and used the same split: untrusted search emits a
  subdivision, a trusted checker replays it).  For squares in squares, the s(11) Lean proof trusts the native
  compiler; our s(13), s(32) and k² − 3 theorems are kernel-only (`decide +kernel`, standard axioms).  The k² − 3
  kernel build is wand125's and was not rebuilt here; k² − 4 and s(21) are still conditional in Lean.
- **Independent checking is above the norm of the published interval proofs.**  Those rest on one implementation
  (Markót 2021 released no code); ours have two or three separately written exact checkers, rejection/mutation
  suites, shipped run records, and external replays (jlevy, wand125).  This is a reproducibility claim, not a
  correctness proof of the programs.
- **What is standard**: exhaustive box subdivision with adaptive refinement, integer-scaled grids to avoid rounding
  (Montanher's scaled tiling, Markót's `[0,24]²`), symmetry folding, LP duals as certificates, informal search
  plus trusted replay.  None of these is a contribution.
- **Not peer reviewed**, unlike every interval-method paper in this list except the arXiv-only ones.  Any write-up
  should cite Markót 2021, Montanher et al. 2019 and Flyspeck as the methodological neighbours.

## Possibly useful techniques

1. **Flyspeck's sharp-inequality recipe for tile germs.**  Exact value at the germ, interval bound away from a
   small neighbourhood U, and a sign condition on U that puts the minimum at the germ.  Our deepest roots (the s32
   wall germ needing depth 30, the 152 `--full` boxes at tile germs) are near-tight neighbourhoods of a corner; a
   local monotonicity lemma could replace deep subdivision there.  Captured mass is piecewise constant in pose
   for point atoms, so this applies only to the segment/area parts of a measure, or to the vertex functions of
   Lemma E (rational in u), which are differentiable.
2. **Monotonicity face reduction** (Flyspeck; Hales–Kusner autodiff): where the Lebesgue/segment mass is smooth
   on a box and a partial derivative has fixed sign, check one face instead of the box.  Relevant to the k² − 3 /
   k² − 4 area covers, where `PIECE`/`LEB` leaves dominate.
3. **Fixed variable order with early stop** (Fekete et al.): split (x, y) first and u only when needed, and never
   split u on parts that x, y already settle.  Worth checking whether our bisection rule already does this; a
   cheap experiment on the s32 root census.
4. **Coarsest-precision replay for the kernel** (Flyspeck: the informal search "tries to find minimal precision
   necessary for the formal verification").  Rounding each box's endpoints and the certificate's data to the
   smallest denominators that still certify could cut `decide +kernel` cost (s32: 13.8 CPU-h; k² − 3: ~7,000
   process-hours on wand125's side).
5. **Cutoff monotonicity / reuse across certificates** (Markót: regions eliminated for n = 33 reused for 32, 31).
   Our analogue: boxes certified for one cover stay certified for any cover that dominates it pointwise; a cache
   keyed on local atom sets could share work between s(45), s(60) and the family bundles.
6. **Canonical-bitset symmetry filtering** (Markót) is weaker than what we do (exact D4 invariance, or the full
   space with no fold); nothing to take.  Montanher's sentinel lemma (two squares are disjoint iff neither contains
   the other's 9 sentinel points) is a primal-side tool; possibly useful for the packer's overlap test or for an
   independent primal check of small cases, not for the covers.
7. **Formal-proof bug finding**: the empty-hexagon formalization found an error in a symmetry-breaking argument.
   Our D4 fold is already in Lean for s(32) and k² − 3; the k² − 4 and s(21) folds should get the same treatment
   before publication.

## Paywalled or not reached (for Evan to fetch)

| Paper | DOI | Why |
|---|---|---|
| Markót & Csendes, "A new verified optimization technique for the 'packing circles in a unit square' problems", SIAM J. Optim. 16 (2005) 193–219 | 10.1137/S1052623403425617 | the n = 28–30 method paper; tiling phases, free circles, fixed-point uniqueness |
| Markót, "Optimal packing of 28 equal circles in a unit square – the first reliable solution", Numer. Algorithms 37 (2004) 253–261 | 10.1023/B:NUMA.0000049472.75023.0a | first interval result |
| Nurmela & Östergård, "More optimal packings of equal circles in a square", DCG 22 (1999) 439–457 | 10.1007/PL00009472 | n = 21–27 proofs, floating-point method |
| Peikert, Würtz, Monagan, de Groot, "Packing circles in a square: a review and new results" (1992), LNCIS | 10.1007/BFb0113271 | n ≤ 20 computer proofs |
| Szabó, Markót, Csendes, Specht, Casado, García, *New Approaches to Circle Packing in a Square* (Springer 2007) | 10.1007/978-0-387-45676-8 | book with program codes; the full interval method |
| Montanher et al., JOGO 73 (2019) 547–565 (published version; preprint read) | 10.1007/s10898-018-0711-5 | check changes since the preprint (Domes, Schichl added) |
| Solovyev & Hales, "Formal verification of nonlinear inequalities with Taylor interval approximations", NFM 2013, LNCS 7871 | 10.1007/978-3-642-38088-4_26 | how Flyspeck's subdivision replay is made cheap in a kernel |

Open and read or available: Markót 2021 (PMC8550790), Markót–Csendes 2006 preprint (`inf.u-szeged.hu/~csendes/intpoly.pdf`),
Montanher et al. preprint (Optimization Online 2018/04/6556), Fekete et al. (arXiv:2103.07258), Flyspeck (arXiv:1501.02155),
Hales–Kusner (arXiv:1602.07220), Musin–Tarasov (arXiv:1410.2536, 1002.1439), Stromquist 2003 (10.37236/1701), Bentz 2016
(arXiv:1606.03746), empty hexagon (arXiv:2403.17370).
