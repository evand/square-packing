# Literature scan, 2026-10-12: composition DP, symmetry, grain boundaries, LNS, cluster tricks, continuation

This scan adds to `lit-optimizers.md`, `lit-statmech.md` and REVIEW-2026-10-11 §4, and does not repeat them.  Markers:
**[verified]** means I fetched the paper or its abstract/bibliographic record this session.  **[read]** means I also read the
relevant section of the full text.  **[unverified]** means the item is from memory or a secondary citation.  No downloaded
code was executed.  Two PDFs were converted to text with `pdftotext`, in the scratchpad.

Two zero-cost **measurements on our own store** were made while writing this note (`runs/store.sqlite`, opened read-only,
best frontier packing per n, 30 ≤ n ≤ 420, 305 packings).  Scripts are in the scratchpad (`sym/sym.py`, `sym/tilt2.py`),
not in the repo.  They answer two of the falsification questions below directly (§2 and §3).

---

## 1. Pallet loading: recursive partitioning, a DP over sub-rectangles

**Birgin, Lobato, Morabito**, "An effective recursive partitioning approach for the packing of identical rectangles in a
rectangle", *J. Oper. Res. Soc.* 61 (2010) 306–320, doi:10.1057/jors.2008.141.  PDF:
https://www.ime.usp.br/~egbirgin/publications/blm.pdf **[verified, read]**.  Its building blocks:

- **Morabito & Morales**, *JORS* 49 (1998) 819–828: the recursive five-block heuristic **[verified via the BLM reference list]**.
- **Lins, Lins & Morabito**, "An L-approach for packing (ℓ,w)-rectangles into rectangular and L-shaped pieces", *JORS* 54
  (2003) 777–789 **[verified abstract]**.
- **Birgin, Morabito & Nishihara**, the data-structure note for the L-approach (*JORS* 2005) **[unverified volume]**.
- **Scheithauer & Terno**, "The G4-heuristic for the pallet loading problem", *JORS* 47 (1996) 511–522 **[verified abstract]**.
- Extension to non-identical pieces: Birgin, Lobato, Morabito, *JORS* 63 (2012) 183–200 **[verified record]**.

**Key tricks.**
1. *Raster points* (normalisation).  Every packing can be shifted so that all coordinates are integer conic combinations
   r·l + s·w, and then reflected to the reduced raster set R = {L − ⟨L − x⟩ : x ∈ S}.  Only cuts at raster points need to
   be tried.  This is the whole complexity story: there are |R_L|·|R_W| rectangle subproblems and |R_L|²|R_W|² L-shaped
   subproblems.
2. *First-order non-guillotine (five-block) cut*.  A quadruple (x1, x2, y1, y2) splits the rectangle into four blocks
   arranged in a pinwheel around a central fifth block.  Guillotine cuts are degenerate cases.  The DP is memoised and
   uses upper bounds (Barnes) to stop recursing when lower bound = upper bound.  Time:
   O(|R_L||R_W|(|R_L|²|R_W|² + |S_L| + |S_W|)).
3. *L-approach*.  The state is an L-shaped piece (X, Y, x, y).  There are nine ways (B1…B9) to divide a rectangle or an
   L into two L's.  Seven apply to a non-degenerate L, each indexed by one raster pair (x′, y′).  BLM add two new
   divisions.  Memory is O(|R|⁴), the binding constraint.  Phase 1 (five-block) supplies lower bounds, and phase 2
   (L-approach) runs only where a gap remains.  Together they solved all of Cover I/II (over 50k instances) to proven
   optimality.  G4 recursively applies a four-block structure in the same spirit.

**Transfer.**  Our route (a) is depth 1 of this DP: two library pieces in opposite corners, with the remainder (an L, or a
rectangle minus two corners) filled axis-aligned.  The generalisation:
- **Leaf value.**  f(w,h) = max(⌊w⌋⌊h⌋, max{m : library piece (m, W, H) with W ≤ w, H ≤ h, either orientation}).
- **Recursion.**  Guillotine, five-block and L-divisions, with each sub-rectangle or sub-L recursively filled.
- **Raster set.**  Our pieces have non-integer sides, so the raster set becomes {s − Σ(piece sides)}.  It stays small if
  the number of non-integer pieces crossed by any horizontal or vertical line is capped (2–3 is enough: records seldom
  stack more).
- **Output.**  For fixed s the DP returns N(s) = the most squares placeable as a disjoint composition, which gives
  s_DP(n) = min{s : N(s) ≥ n} by bisection.  The DP output is a seed for our exact polish, which closes seams (pieces
  in records interlock, and the DP does not model that).
- **Exhaustion bound (a by-product).**  If the DP over the library cannot beat the register at n, then **no
  rectangular-boundary composition of library pieces can**, so route (a) is exhausted at that n.

**Falsification test (about 1 day).**  Implement guillotine plus five-block over the library at the register s of every
n ≤ 400.
- *Pass:* it reproduces every route-(a) record already registered (the Mondrian/SQUISH ones), and it finds n where
  N(s_reg − 1e-6) ≥ n before polish.
- *Fail:* if it finds none, and the L-approach adds nothing over five-block, the composition space is saturated by
  others and we should stop investing in (a).

## 2. Crystal structure prediction: symmetry-constrained random starts

**Sources.**
- **AIRSS**: Pickard & Needs, *PRL* 97 (2006) 045504.  Review: *J. Phys. Condens. Matter* 23 (2011) 053201,
  doi:10.1088/0953-8984/23/5/053201 **[verified bibliographic]**.  "Hot-AIRSS": Pickard, arXiv:2407.06294 **[verified listing]**.
- **CALYPSO**: Wang, Lv, Zhu, Ma, arXiv:1205.2264 **[verified bibliographic]**.
- **USPEX fingerprints**: Oganov & Valle, *J. Chem. Phys.* 130 (2009) 104504.  Valle & Oganov, *Acta Cryst. A* 66 (2010)
  507 **[verified bibliographic]**.
- **Why symmetry and low energy co-occur**: Wales, "Symmetry, near-symmetry and energetics" (*CPL* 1998, **[unverified
  citation]**).  A 2026 formalisation: Schneider, arXiv:2603.00011 **[verified, abstract]**.  It is an argument, not a
  theorem, and does not discuss packings.

**Key trick.**  Pick a random space group (for us, a subgroup of D4), place a random half of the atoms, and generate the
rest as images.  Relaxation may then land in a *supergroup*.  Sensible structures add minimum separations and building
units.  Diversity comes from removing duplicates by fingerprint (USPEX/CALYPSO).  We already have this: the image
fingerprint in PACKER.md.

**Evidence for our problem (measured today, 305 frontier packings, tolerance 1e-6 in position and angle).**
- 84 / 305 best-known packings are **exactly invariant** under a non-trivial element of D4.  By element: diagonal mirror
  Md 53, C2 17, anti-diagonal mirror 7, axis mirrors 7.
- **53 of these contain tilted squares**, so they are not just grids.  Examples: n = 50, 52, 65, 67, 82, 101, 104, 109,
  122–125, 145–150, 170–174, 197–203, 226–233, 257–265, 290–298.
- The exactly symmetric share is flat across sizes (27/70, 28/100 and 27/100 for n in 30–99, 100–199 and 200–299), then
  2/35 for n ≥ 300.  The drop at n ≥ 300 is probably because that part of the frontier is young.
- Another 120 of the 221 asymmetric packings have at least 95% of squares matched by one mirror.  Part of this is the
  trivial matching of corner-anchored axis blocks, so it is weaker evidence.
- lit-provenance.md records that de Winter's and Schadt's improvements came from *relaxing* the symmetry of
  Ellsworth's analytic packings.  So symmetric funnels lead to records, and breaking the symmetry gives the last step.

**Transfer.**  Generate and quench in the Md-invariant (or C2-invariant) subspace:
- half the variables;
- squares on the diagonal constrained to θ ∈ {0, 45°};
- an even/odd n parity rule fixing how many squares sit on the axis.

Then release the symmetry and polish.  This also addresses "no constructor reaches new funnels": the subspace landscape
has far fewer minima.

**Falsification test (S, hours).**  Take battery n whose record is Md-symmetric (for example 101, 125, 148, 200).  Run
equal-budget quench campaigns, symmetric versus free, from the same generator, and compare hit rates of the record
basin.
- *Kill condition:* symmetric hit rate is no better than free.
- *Then:* run symmetric campaigns at n whose record is asymmetric and check that relaxing the symmetry reaches the
  record value or better.

## 3. Grain boundaries, CSL and commensurate tilt angles

**Sources.**
- Coincidence and similarity rotations of Z²: Glied & Baake, arXiv:0808.0109.  Baake & Zeiner, arXiv:1709.07317
  **[verified bibliographic]**.  Exact coincidence rotations of the square lattice are the rotations with rational
  cos and sin, i.e. **Pythagorean angles**: 36.87° (Σ5), 22.62° (Σ13), 28.07° (Σ17), 16.26° (Σ25), and so on.
  **atan(1/2) is *not* a CSL angle.**  It is a *similarity* rotation: rotation by atan(b/a) maps Z² onto
  (a²+b²)^{-1/2}·(a sublattice of index a²+b²).
- 2D square-lattice tilt boundaries show no deep energy cusps between sampled angles (arXiv:1806.00700 **[verified
  listing; details from the search summary]**).

**Mechanism for packings (my derivation, not from a paper).**  A tilted block at tan θ = b/a, cut by an axis-parallel
line, presents a boundary of period √(a²+b²).  Block widths between two parallel walls lie in (1/√(a²+b²))·Z.  Short
periods (√5, √10) mean short repeat units of waste along a wall or a seam.  Pythagorean angles additionally make the
period an integer, so they commensurate with an axis block.

**Measured today (tilted squares folded to [0°, 45°], ±0.15°; 8,100 tilted squares in 149 packings).**

| Angle | Tilted squares within ±0.15° | Records containing it |
|---|---|---|
| 45° | 2,408 | many |
| atan(1/3) = 18.43° | 231 | 13 |
| atan(1/2) = 26.57° | 192 | 21 |
| ≈ 24.4° | 178 | 13 |
| atan(3/4) = 36.87° (the Σ5 angle) | 95 | 8 |
| atan(8/15) = 28.07° | 86 | 14 |
| atan(5/12) = 22.62° | 19 | 4 |
| atan(7/24) = 16.26° | 7 | 3 |
| atan(2/3) = 33.69° | 57 | 4 |

The ≈ 24.4° peak sits at atan(5/11), but that ratio is probably not its cause.

**Verdict.**  The CSL/Pythagorean hypothesis is **falsified as an organising principle**: two of the four Σ ≤ 25 angles
are depleted.  Simple rational slopes (short boundary period) are enriched, and atan(2/3) is weak.

**Transfer.**  Give generators a tilt alphabet of {45°, atan 1/3, atan 1/2, ≈24.4°, atan 3/4, 28.1°} instead of 38–45°
only.  Before adding it, find the ≈ 24.4° family: which n contain it, and is it a known family?

**Remaining test (zero compute).**  Repeat the histogram restricted to records registered after 2026-09 (the newest
methods), to see whether the alphabet is drifting.

## 4. Large neighbourhood search, ruin and recreate, window MIP

**Sources.**
- **Schrimpf, Schneider, Stamm-Wilbrandt, Dueck**, "Record breaking optimization results using the ruin and recreate
  principle", *J. Comput. Phys.* 159 (2000) 139–171 **[verified record]**.  Shaw (CP 1998) introduced LNS **[unverified]**.
  Key point: large, structured destruction (spatial or radial ruin) followed by greedy rebuild, accepted by threshold.
  Big ruins beat small kicks on constrained landscapes.
- **Gomes & Oliveira**, "Solving irregular strip packing problems by hybridising simulated annealing and linear
  programming", *EJOR* 171 (2006) 811–829, doi:10.1016/j.ejor.2004.09.008 **[verified record]**.  SA chooses a
  *combinatorial* move (swap or relocate), and an LP compaction or separation model, with relative positions fixed,
  produces the neighbour exactly.  This is structurally what our SLP does.  The difference is that the move is
  combinatorial (which side of which neighbour) and the LP is the evaluation.
- **Kallrath & Rebennack**, "Cutting ellipses from area-minimizing rectangles", *J. Glob. Optim.* 59 (2014) 405–437,
  doi:10.1007/s10898-013-0125-3 **[verified record]**.  Its *polylithic* scheme: exact NLP (separating hyperplanes) is
  infeasible beyond about 14 objects, so objects are added sequentially, re-solving the small NLP each time in random
  order (GRASP-like).
- **Cherri, Carravilla, Toledo**, three-phase matheuristic with exact MIP subproblems, *Pesquisa Operacional* (2016)
  **[verified listing]**.
- Recent exact work, useful only for small windows: Lastra-Díaz & Ortuño, parallel branch-and-bound for nesting,
  arXiv:2503.21009 **[verified listing]**.  The combinatorial Benders decomposition for convex polygons, arXiv:2609.23461,
  is **[listing only]**.

No 2020–2026 paper does window LNS for *continuous-rotation* congruent squares.

**Transfer.**
- **Ruin.** Remove the k (5–15) squares inside a disc or strip around the largest hole, a seam, or the most-stressed
  contact.
- **Recreate.** Choose one of three, in increasing order of strength:
  - a greedy reinsertion at dynamic-lattice sites (§5);
  - an exact small global NLP over the k removed squares, in the Berthold et al. separating-hyperplane formulation
    (already in lit-optimizers; their code is public), with everything else fixed and objective = maximise the minimum
    clearance;
  - an MIP over "which side of which neighbour" disjunctions, with angles fixed to the alphabet of §3.
- **Accept.** Use our exact polish on the full packing as the acceptance step.

**Falsification test (M).**  On x110a and five battery n, compare 200 ruin-recreate steps against 200 basin-hopping kicks
at equal wall time.  Count new exact minima below the incumbent and new contact topologies.

## 5. Atomic-cluster global optimisation: lattice attractors and surface operators

**Sources.**
- **Dynamic lattice searching (DLS)**: Shao, Cheng, Cai, *J. Comput. Chem.* 25 (2004) 1693–1698 **[verified via
  reference list]**.  Interior-operation two-phase variant (DLS-TPIO): *Sci. China Chem.* (2011),
  doi:10.1007/s11426-011-4280-4 **[verified listing]**.
- **Current best description**: Lai, Hao, Guo, Wen, Fu, "Structural optimization of atomic clusters using iterated
  dynamic lattice search", *Comput. Phys. Commun.* (2025), doi:10.1016/j.cpc.2025.109655 **[verified, read]**.  It
  improved 47 of 300 Ag cluster records.

**Key trick (Lai et al. §2.4).**
- *Freeze the core.*  Fix the low-energy (interior) atoms.
- *Build the dynamic lattice.*  Run 2N "probe atom" relaxations: drop one extra atom at random on the surface and
  locally minimise *its* energy only.  The resulting stable sites, together with the sites of the N_mov highest-energy
  atoms, form the **dynamic lattice** (rebuilt every iteration).
- *Search the lattice.*  Assign the N_mov active atoms to lattice sites by multi-start greedy: move the worst atom
  (per-atom energy E(i) = E(X) − E(X∖i)) to the best free site.  Then L-BFGS-relax the m best lattice solutions.
- *Escape.*  The outer loop is iterated local search with a surface perturbation and Metropolis acceptance.  The
  "interior operator" in TPIO moves a high-energy surface atom inside.

**Related.**
- **Northby** lattice search, *J. Chem. Phys.* 87 (1987) 6166, doi:10.1063/1.453492 **[verified]**: occupancy search
  on a fixed icosahedral lattice, then relax.
- **Romero, Barrón, Gómez**, *CPC* 123 (1999) 87: lattice-based GA plus growth for n = 148–309 **[verified via
  reference list]**.
- **Doye, Miller, Wales**, *J. Chem. Phys.* 110 (1999) 6896, arXiv:cond-mat/9808265 **[verified]**: the LJ38 double
  funnel.  The global minimum sits in a narrow, deep funnel, and the wide funnel captures the dynamics.  That is our
  picture of k-row attractors versus record funnels.  The cluster community's answer was *biased* lattice searches
  (Northby, DLS with constructed core) *plus* unbiased ones, run separately.

**Transfer.**  This directly targets two of our weaknesses.
- *Fill.*  Our fill only uses wall-anchored lattice lines.  A probe-square relaxation, one free unit square minimising
  its own overlap against the fixed packing and walls, finds sites anchored to *tilted* neighbours too, including tilted
  sites in holes.
- *Absorbing an extra square.*  To add square n+1 without jamming to the grid:
  1. build the dynamic lattice of probe sites;
  2. pick the N_mov most-stressed squares (or the squares of the seam);
  3. greedily reassign;
  4. polish.

  This is exactly the growth step that the LJ codes use to go from N to N+1.  The "interior operator" corresponds to
  moving a square out of a loose rattler region into the tight interior.

**Falsification test (S–M, about 1 day).**  Remove one square from 20 battery records (n+1 → n), add it back with
dynamic-lattice insertion, and compare with our current quench.
- *Recovery:* how often is the record recovered?
- *Growth:* take records at n and attempt n+1 (Fang/Couzo style).
- *Kill condition:* recovery is no better than quench.

## 6. Soft-to-rigid, shape continuation, inflation

**What "Nakajima" refers to (resolved).**  Yohei Nakajima, "Soft-to-rigid packing", 2026:
https://yoheinakajima.github.io/soft-to-rigid-packing/, code at github.com/yoheinakajima/soft-to-rigid-packing (MIT;
not run) **[verified, page and README fetched]**.
- *Shape.*  Each object is "a cube of half-side ½ − r dilated by radius r" (a rounded square).
- *Protocol.*  Compress balls (r = ½) under wall pressure, lower r linearly to 0 while pressure continues, legalise by
  scaling centres apart, then run SLSQP with hard corner-in-box and per-pair separating-plane constraints.
- *Results.*  It rediscovers s(11), s(17), s(18), s(19), s(26), s(27) from random starts.  Longer hardening raised the
  n = 26 hit rate from 2.5% to 10.4%.  An octagon intermediate did not help.  It improved the n = 12 *cube* record
  (8 Oct 2026) but no square record.
- *Use by others.*  dostashev (#501) uses it on a GPU as a re-hardening step for a constructed V-band, not as a
  generator.

This matches REVIEW §4's "≈ our corner-radius anneal".

**Literature for the same idea.**
- **Avendaño & Escobedo**, "Phase behavior of rounded hard-squares", *Soft Matter* 8 (2012) 4675 **[verified
  bibliographic]**.  Depending on the rounding, compressed rounded squares go to a hexagonal rotator phase, a
  tetratic-to-square phase, or (in a narrow crossover) polycrystals.  Experiments (Zhao, Bruinsma, Mason, *PNAS* 2011,
  **[unverified]**) show a **rhombic** crystal.  Implication: *the rounding schedule decides which intermediate order
  forms.*  Mid-range r gives sheared or rhombic rows, i.e. tilted bands.  Small r gives square order, the k-row
  attractor.  A schedule that dwells at the crossover r is the physics version of "merge two tilted bands into a V".
- **Doye, Leary, Locatelli, Schoen**, "Global optimization of Morse clusters by potential energy transformations",
  *INFORMS J. Comput.* 16 (2004) 371 **[verified via reference list]**.  Optimise on a *transformed* (smoother,
  longer-range) potential, then map minima back.  This is the cluster-community form of shape homotopy, and it worked
  for the hard short-range Morse cases.
- Digital alchemy (shape as a thermodynamic variable; van Anders et al., *ACS Nano* 2015) **[unverified]** is the same
  idea pushed to a per-particle shape variable.

**Transfer.**
1. *Rounding schedule.*  Our corner-radius anneal (tried) should be re-tested as a **schedule** experiment: dwell at the
   rounded-square crossover r rather than ramping linearly.
2. *Partial hardening.*  Apply soft-to-rigid only to an *added* square or a seam region (rounded insert squeezes in;
   neighbours stay rigid).  That attacks the "extra square jams to grid" weakness.

**Falsification test (S).**  Measure the n = 26–29 and n = 50 hit rates for linear versus dwell schedules (3 dwell
radii).  Kill if no schedule beats linear by 2×.

## 7. Circle and equal-object records, nesting, AI-driven packing, 2020–2026

**New since the earlier notes:**
- **ImprovEvolve**: Kravatskiy, Khrulkov, Oseledets, arXiv:2602.10233 (ICML 2026 AI4Math workshop) **[verified]**.  The
  LLM evolves only three operators (init, local improvement, perturbation) inside a fixed basin-hopping loop.  It set new
  hexagon-in-hexagon records (n = 11, 12, 15, 16; 14, 17, 23 after light expert tuning).  This is the closest published
  analogue to "LLM evolves the operators of *our* engine" rather than a whole optimiser.  It is more sample-efficient
  than whole-program evolution because the expensive numerics stay ours.
- **Geometry-aware RL for 2D irregular nesting with free continuous rotation**, arXiv:2606.10611 (June 2026, preprint)
  **[listing only]**.  It claims parity with Sparrow, not better.  No transfer beyond confirming that Sparrow is still
  the state of the art.
- **Packomania** updates in 2026 (Specht; Zhou, He, He 2024 for unequal circles in a circle) **[listing]**.  No new
  method papers for equal circles in a square after Lai et al. 2022 and Zhou–He–He 2025 (already noted).
- **López & Beasley**, formulation space search for unequal rectangles and squares *with rotation* in a circle,
  arXiv:1802.07519 **[listing]**.  It is the only rotated-square FSS application found.
- **Equal rectangles or squares in rectangles with free rotation**: nothing new 2020–2026.  The rotation literature
  (Korf; Simonis & O'Sullivan; Caprara et al.) is about 90° rotation only.
- **Evolution Fine-Tuning (371 optimisation tasks)**, arXiv:2606.29082 **[listing only]**.  I did not confirm whether
  it includes squares-in-square.

**Honest gaps.**  I found no paper on "rectangle library + composition" for rotated congruent squares, and none on
window MIP for continuous-rotation squares.  Grain-boundary theory gives no quantitative prediction for finite packings
beyond the angle alphabet in §3.

---

## Ranked top 8 by expected value per engineering-day

1. **Symmetric-subspace generation, then relax** (§2).  S.  Evidence measured today: 84/305 records are exactly D4-invariant,
   53 with tilts, 53 via the diagonal mirror.  Test: symmetric versus free hit rates on Md-symmetric battery n.
2. **Dynamic-lattice insertion (probe-square sites plus greedy reassignment of the worst squares)** (§5).  S–M.  It fixes
   fill's wall-lattice limitation and is the N → N+1 growth step.  Test: remove-one/re-add recovery on 20 records, then
   n → n+1.
3. **Tilt alphabet from data**: 45°, atan(1/3), atan(1/2), ≈24.4°, atan(3/4), 28.1° in generators (§3).  S.  The
   histogram is done, and the CSL hypothesis is falsified.  Remaining work: identify the 24.4° family, then A/B the
   generators on the battery.
4. **Recursive-partitioning DP over the rectangle library** (five-block plus L-approach, with raster sets of piece-side
   sums) (§1).  M.  It generalises route (a) and certifies when (a) is exhausted at an n.
5. **Ruin-and-recreate windows with exact small subproblems** (Gomes–Oliveira LP-as-neighbour; Berthold-style global NLP
   on k ≤ 10 squares) (§4).  M.
6. **Rounding-schedule dwell / partial soft-to-rigid on inserted squares** (§6).  S, but we already tried a linear
   anneal, so the prior is lower.
7. **ImprovEvolve-style operator evolution inside our engine** (§7).  L.  Its value rises once items 1–3 define a good
   operator interface.
8. **Fingerprint-niched population** (USPEX/CALYPSO, §2).  S, but low marginal value: our fingerprint and
   lineage machinery already do most of it.
