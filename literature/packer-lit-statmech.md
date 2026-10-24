# Literature: stat-mech / global-optimization methods and basin census, for the square packer

Written 2026-10-04 (literature search, no code run).  The context is our packer (`PACKER.md`):
fixed side s, squared-SAT overlap energy, L-BFGS quench, basin hopping, and a descending target side.
The failure mode: many basins jam with a full axis-aligned row of k squares, so the basin minimum is
exactly side k.  At s = k−ε their energy (~ε²) looks the same as a genuine near-miss, and record basins
appear exponentially rare.

Markers: **[V]** means I checked the bibliographic data (authors, venue, year) by web search this session.
**[R]** means it is from memory and was not re-checked, so verify it before citing externally.
"Code" lists public repositories where they are known.  None were run.

---

## Topic A: global and stochastic optimization

### A1. Replica exchange with side length (pressure) as the exchange variable
- **Odriozola, "Replica exchange Monte Carlo applied to hard spheres", J. Chem. Phys. 131, 144107 (2009)**, arXiv:1010.2923 [V title/arXiv; R volume/page].  The replicas live in an extended isobaric ensemble and are exchanged in *pressure*, not temperature: "expanding volume instead of increasing temperature" unblocks crowded hard systems.  Follow-ups apply it to hard ellipsoids (arXiv:1105.2789) and to massive REMC at high pressure.
- **Okabe, Kawata, Okamoto, Mikami, Chem. Phys. Lett. 335, 435 (2001)**, isobaric REMC [R].
- **Mechanism for us.** Run a ladder of fixed sides s₁ < s₂ < … < s_M (for example from s_target up to k+0.3), with soft-energy MC/basin hopping at each rung, and swap configurations with the Metropolis rule on ΔE (rescaled into each rung's box).  A row-jammed configuration that reaches a looser rung can un-jam there: tilts open up and the row can break.  It returns to the tight rung only if it re-compresses well.  This is the cheapest way to give "decompress, then re-compress differently" detailed-balance-like bookkeeping, and it replaces ad-hoc shake/restart.
- **Cost.** M× the compute, which parallelizes trivially across our 16 cores.  Hard spheres need a pressure ladder about √N rungs dense; we likely need about 8–16 rungs.

### A2. Population annealing, nested sampling, and multilevel splitting in the side length
- **Population annealing:** Hukushima & Iba, AIP Conf. Proc. 690, 200 (2003) [R]; Machta, PRE 82, 026704 (2010) [R]; applied to hard-sphere mixtures by Callaham & Machta, PRE 95, 063315 (2017) [R].  Mechanism: anneal a population of R configurations through a schedule of decreasing s, reweight/resample at each step, and keep the families diverse.  The output includes a free-energy-like estimate of how much configuration space survives at each s.
- **Nested sampling:** Skilling, Bayesian Analysis 1, 833 (2006) [R]; Pártay, Bartók, Csányi, J. Phys. Chem. B 114, 10502 (2010) [R]; hard spheres: Pártay et al., PRE 89, 022302 (2014) [R].  **Superposition-enhanced NS:** Martiniani, Stevenson, Wales, Frenkel, *PRX* 4, 031034 (2014), arXiv:1402.6306 [V], which mixes known minima into NS to fix broken ergodicity.  Applied to us: keep K live configurations.  Each has a "feasible side" (the smallest s at which it quenches to E < tol).  Repeatedly discard the worst, clone a survivor, and decorrelate it by MC *constrained to side ≤ the discarded value*.  The volume compresses geometrically, X_i ≈ e^{−i/K}, so NS gives a direct estimate of the fraction of configuration space that can reach side s.  That is a census quantity (see Topic B).
- **Adaptive multilevel splitting:** Cérou & Guyader, Stoch. Anal. Appl. 25, 417 (2007) [R].  **Forward flux sampling:** Allen, Warren, ten Wolde, PRL 94, 018104 (2005) [R].  These estimate the probability that a trajectory reaches a rare region (here, side < s*) by cloning the trajectories that get furthest along an order parameter.  The algorithm is nearly identical to "keep the best K, re-branch", and it gives an *unbiased probability estimate* as a by-product.
- **Why it helps the attractor.** Selection by the *final side achieved*, rather than by energy at one fixed s, means row-jammed families stall at k and are culled.  Clones of rare families are amplified, which deals with exponential rarity: splitting turns a probability p into about log(1/p) stages.

### A3. Flat-histogram methods (Wang–Landau, multicanonical, landscape paving, metadynamics)
- **Wang & Landau, PRL 86, 2050 (2001)** [R]; **Berg & Neuhaus, PRL 68, 9 (1992)** multicanonical [R].
- **Landscape paving:** Hansmann & Wille, PRL 88, 068105 (2002) [R].  It is tabu search turned into a sampler: the energy is replaced by E + f(H(q)), where H is the visit histogram of a collective variable q.  Visited regions become uphill.
- **Mechanism for us.** Pick a collective variable that *identifies* the attractor.  One choice is R = the maximum, over horizontal and vertical lines, of the number of axis-aligned (|θ| < δ) squares whose interiors the line crosses; R = k means a full row.  Another is the number of near-axis-aligned squares.  Penalize visited R, or simply add a bias that is zero for R < k and positive for R = k.  This attacks the degeneracy at its root: the energy cannot tell ε² row-jams from near-misses, but a structural collective variable can.
- **Cost.** Little, for a bias or rejection test.  True Wang–Landau convergence in about 600 dimensions is expensive, and a flat histogram is not the goal anyway.

### A4. Swap / non-local moves, cut-and-splice, population basin hopping
- **Swap MC:** Ninarello, Berthier, Coslovich, *PRX* 7, 021039 (2017), doi:10.1103/PhysRevX.7.021039 [V].  Swapping particles of different sizes sped equilibration by more than 10¹⁰ in polydisperse glasses.  Our squares are identical, so a literal swap is a no-op.  The useful analogues are **remove–reinsert (vacancy hopping) moves**: take a square out of the jammed row and reinsert it in the largest hole, then quench.  Polydispersity helps swap MC because size moves are cheap.  Our analogue is an **orientation-transfer move**: rotate a whole column or block to a common tilt.
- **Cut-and-splice GA:** Deaven & Ho, PRL 75, 288 (1995) [R].  Cut two parents with a random plane, join the halves, and relax.  For us: cut along a random line (or diagonal), take one side from each of two good packings, repair the square count, and quench.  Records for n ≈ 50–300 are piecewise-structured (axis-aligned blocks plus tilted strips), so crossover recombines good substructures.
- **Packing-specific evidence:** Addis, Locatelli, Schoen, "Disk packing in a square: a new global optimization approach", *INFORMS J. Comput.* 20, 516 (2008), doi:10.1287/ijoc.1080.0263 [V].  They argued that the landscape is funnel-like and used monotonic basin hopping, plus a population variant, to improve **32 putative optima for n ≤ 130** (smallest n = 53).  This is the closest published analogue to our engine at our n.  Grosso, Jamali, Locatelli, Schoen, J. Glob. Optim. 47, 63 (2010) [R] applied population basin hopping with a dissimilarity measure to circles in a circle.
- **Squares-in-square specifically:** Gensane & Ryckelynck, "Improved dense packings of congruent squares in a square", *Discrete Comput. Geom.* 34, 97 (2005) [R details; V that their program found new square packings, e.g. an alternate n=18].  They used a stochastic perturbation/shrink algorithm, the closest prior work to ours.

### A5. Minima hopping (Goedecker)
- **Goedecker, J. Chem. Phys. 120, 9911 (2004)** [R].  Escape from each minimum with a short MD run at kinetic energy E_kin.  Accept or reject the new minimum by an energy threshold, and keep a *history of visited minima*: revisiting a known minimum raises E_kin, and finding a new one lowers it.
- **Why relevant.** The feedback is an automatic way to stop falling back into the same jammed row.  The MD escape direction (softened along low-curvature modes) favours collective moves such as row tilting over random single-square kicks.  The method needs a minimum fingerprint (see B4).  Cost: modest.

### A6. Compression protocols: Lubachevsky–Stillinger and the Torquato–Jiao adaptive shrinking cell
- **Lubachevsky & Stillinger, J. Stat. Phys. 60, 561 (1990)** [R].  Event-driven MD while the particles grow.  Slow growth crystallizes; fast growth jams in disordered states.  Lubachevsky and Graham used it to find many disk-in-square/triangle records [R].
- **Torquato & Jiao, Nature 460, 876 (2009)** and PRE 80, 041104 (2009), arXiv:0909.0940 [V].  The adaptive shrinking cell treats the cell's deformation and the particle moves as optimization variables, solved by MC and/or sequential LP.
- **Assessment.** These protocols *generate* jammed states.  They do not escape them.  For identical squares, slow compression tends to *produce* axis-aligned rows, which is our attractor.  The useful idea from the adaptive shrinking cell is to give the container extra degrees of freedom.  For example, optimize in a rectangle w×h with w·h decreasing and penalize |w−h| gradually: a row that is jammed in w can still let h shrink.  Cheap to try.

### A7. Event-chain and hard-particle MC (HPMC)
- **Event-chain MC:** Bernard, Krauth, Wilson, PRE 80, 056704 (2009) [R].  For polyhedra: Klement & Engel, J. Chem. Phys. (2021), "Newtonian event-chain MC and collision prediction with polyhedral particles" [V title].  Hard regular polygons, squares included: Anderson, Antonaglia, Millan, Engel, Glotzer, PRX 7, 021001 (2017) [R].  **Code: HOOMD-blue HPMC** (glotzerlab/hoomd-blue) supports convex polygons with rotation.
- **Assessment.** These are excellent for *equilibrium* sampling of hard squares at fixed density.  Combined with pressure replica exchange (A1), HPMC NPT is the most principled "anneal hard squares into a box" method.  But equilibrium is not the goal: the record is a single zero-entropy state, and in this regime equilibrium MC is slow.  It is more useful as a decorrelation kernel inside nested sampling or splitting (A2) than as a search on its own.

### A8. CMA-ES
- **Hansen & Ostermeier, Evol. Comput. 9, 159 (2001)** [R].  It learns a full covariance in d dims, costing O(d²) memory and needing about 10·d evaluations to adapt.  At d = 3n ≈ 300–900 with hard constraints it is poor.  **Use it only on low-dimensional parametrized families**: tilt angles, strip widths and offsets of a structured seed template (about 5–20 parameters), each evaluated by a quench.

### A9. LLM-guided program evolution
- **AlphaEvolve:** Novikov et al., Google DeepMind white paper (May 2025) [R]; Georgiev, Gómez-Serrano, Tao, Wagner, "Mathematical exploration and discovery at scale", arXiv:2511.02864 (2025) [V].  It improved circle-packing sum-of-radii for n=26 and n=32, among about 67 problems.  The LLM evolves *programs* (seed constructors plus local optimizers) that run under a time budget; the evolved heuristic, not the LLM, does the numerics.
- **ShinkaEvolve:** Lange, Imajuku, Cetin (Sakana AI), arXiv:2509.19349, ICLR 2026 [V].  Open source (Apache-2.0).  It beat AlphaEvolve's n=26 result in about 150 proposals.
- **OpenEvolve:** open-source AlphaEvolve clone (github.com/algorithmicsuperintelligence/openevolve) [V].
- **"Discovery Loop":** W. Sander, arXiv:2609.05093 (Sept 2026) [V].  It improved **10 Packomania sum-of-radii records for N = 101–114** for about $28.  Code: github.com/ucsandman/discovery-loop.
- **FrontierCS:** arXiv:2512.15699 (Dec 2025) [V paper].  It reportedly includes "n unit squares in a square, rotations allowed" as a benchmark task.  **[unverified: seen only in a search snippet]**
- **Other non-LLM learned starting points:** Cassioli, Di Lorenzo, Locatelli, Schoen, Sciandrone, "Machine learning for global optimization", Comput. Optim. Appl. 51, 279 (2012) [R].  They learned which starting points lead to good basins, for circle packing.
- **Assessment.** The best fit is to *evolve the seed generator*: a program that emits structured initial packings (blocks, tilted strips, staircase cuts) which our quench then finishes.  Records at n ≈ 50–300 are highly structured, and LLMs write structured constructors well.  The evidence is real but limited to circle problems so far.  Cost: low dollars, plus our engine as the evaluator.

---

## Topic B: counting jammed states and basin census

### B1. Core papers
- **Stillinger & Weber, PRA 25, 978 (1982)** (inherent structures) [R].  **Stillinger, "Exponential multiplicity of inherent structures", PRE 59, 48 (1999)** [V].  Bounds show that the number of distinct inherent structures grows as e^{αN} at fixed density.
- **Gao, Blawzdziewicz, O'Hern, PRE 74, 061304 (2006)**, arXiv:cond-mat/0606224 [V].  For bidisperse disks with N ≤ 14 they enumerated nearly all mechanically stable packings by compress/decompress plus minimization.  The packings are **not equiprobable**: frequencies span orders of magnitude.
- **Xu, Frenkel, Liu, PRL 106, 245502 (2011)**, arXiv:1101.5879 [V].  A free-energy MC method measures basin volumes and so estimates the *number* of minima without enumerating them.  It was validated against enumeration for small N.  Packing entropy is extensive.
- **Asenjo, Paillusson, Frenkel, "Numerical calculation of granular entropy", PRL 112, 098002 (2014)**, arXiv:1312.0907 [V].  This version outperforms direct enumeration by more than 200 orders of magnitude, for up to 128 polydisperse disks.  The entropy needs a 1/N! factor and is then extensive.
- **Martiniani, Schrenk, Stevenson, Wales, Frenkel, "Turning intractable counting into sampling…", PRE 93, 012906 (2016)**, arXiv:1509.03964 [V title/venue; R article number].  This is the 3D jammed-sphere count.  It found a **strong power-law correlation between a packing's pressure and its basin volume**.
- **Martiniani, Schrenk, Stevenson, Wales, Frenkel, "Structural analysis of high-dimensional basins of attraction", PRE 94, 031301(R) (2016)**, doi:10.1103/PhysRevE.94.031301 [V].  The basin-volume method uses a ladder of harmonic restraints around the minimum, a "is it still in the basin?" oracle (quench and compare), and multistate Bennett acceptance ratio (MBAR) to chain the free-energy differences.  It also gives the radial profile of the basin.
- **Martiniani, Schrenk, Ramola, Chakraborty, Frenkel, "Numerical test of the Edwards conjecture shows that all packings are equally probable at jamming", Nat. Phys. 13, 848 (2017)**, arXiv:1610.06328 [V].  After correcting for basin volume, packings are equiprobable only at unjamming.  Above unjamming they are not.
- **Suryadevara, Casiulis, Martiniani, "The basins of attraction of soft sphere packings are not fractal"**, arXiv:2409.12113 (2024; rev. 2026) [V].  True gradient-flow basins are smooth.  **FIRE and L-BFGS distort basin geometry badly**, and the "basin" is minimizer-dependent.  For us the relevant basin *is* the one our L-BFGS protocol defines, but measurements are only valid for that protocol.
- **Landscape tools:** Wales, *Energy Landscapes* (CUP 2003) [R]; GMIN/OPTIM/PATHSAMPLE (Wales group, Fortran); **pele** (github.com/pele-python/pele) with **mcpele** (MC/parallel tempering, used for the basin-volume work) [V].  Disconnectivity graphs, from Becker & Karplus, J. Chem. Phys. 106, 1495 (1997) [R], are the standard way to see whether a landscape is a funnel or a "palm tree".
- **Hard-particle IS counting by tiling:** Ashwin & Bowles, J. Non-Cryst. Solids (2009), arXiv:0903.1874 [V].  It builds jammed packings of confined hard discs from compatible local "tiles".  Conceptually close to our block/strip structure.
- **Counting local optima in combinatorial landscapes:** Garnier & Kallel, SIAM J. Discrete Math. 15, 122 (2002) [R], estimate the number of optima from repeat-hit statistics of random restarts, assuming multinomial hits weighted by basin size.  Reeves & Eremeev (2004) [R] and later work report errors of about 3–10% on test landscapes.  Ecology estimators: **Good, Biometrika 40, 237 (1953)** (Good–Turing) [R]; **Chao, Scand. J. Stat. 11, 265 (1984)** (Chao1) [R]; Chao & Lee, JASA 87, 210 (1992) (ACE) [R].

### B2. Scaling with n and the shape of the basin-volume distribution
1. **The number of jammed states is exponential in N** (Stillinger 1999; Xu et al. 2011; Asenjo et al. 2014).  The entropy is extensive once the 1/N! factor is included.  For us, also divide out the 8 square symmetries.
2. **Basin volumes are extremely broad**, spanning many orders of magnitude even for N ≈ 10–30 (Gao et al. 2006; Martiniani et al. 2016).  A random quench therefore samples minima in proportion to basin volume, not uniformly.  **Hit frequency equals basin volume divided by total volume**, for the starting distribution actually used.
3. **Basin volume anti-correlates with "constrainedness"**: higher-pressure packings have smaller basins, as a power law (Martiniani 2016).  The analogy is mine **[speculation]**: our record packings are the most compressed states at a given s, so we should expect their basins to be among the smallest.  That is consistent with "record basins are exponentially rare", and with the expectation that this gets worse with n.
4. Because of (2), species-richness estimators (Chao1, ACE) **underestimate badly** under heavy-tailed abundances.  They give lower bounds.  What remains well-estimated is the **missing probability mass** (Good–Turing F₁/M), which is the quantity that matters for search.

### B3. Estimating the probability of a rare basin (directly answers "how rare is the record?")
Take a known record configuration X* (from Friedman's page or our own) at side s* + δ:
- **Basin volume with our own quench** (Martiniani 2016 recipe): run MC with restraint U_λ = λ|X−X*|² for a ladder of λ, under the hard constraint "a quench from X returns to X* (up to symmetry)".  Combine with MBAR to get log V(X*).  Do the same for a few typical row-jammed minima.  **The ratio of volumes is the ratio of hit probabilities** under uniform starts.
- **Cheaper proxy:** the return probability p(r) under isotropic perturbations of radius r, and r₁/₂ where p = ½.  In d ≈ 3n dimensions, log V ≈ d·log r_eff, so even a modest ratio of radii means an exponential ratio of volumes.  This also suggests basin hopping kick sizes: kicks much larger than r₁/₂ of the record basin can never land in it, except by chance.
- **For our actual (non-uniform) start protocol,** estimate P(protocol reaches side < s*) by **adaptive multilevel splitting** on the order parameter "best feasible side so far".  This needs about log(1/p)/log(2) levels instead of 1/p runs.

---

## The 5 most promising ideas for the integer-side attractor

1. **Detect the attractor structurally, not energetically.**  Compute a row collective variable R (the maximum number of near-axis-aligned squares crossed by one horizontal or vertical line).  Once R = k is reached, either reject, apply a landscape-paving penalty, or apply a targeted escape move (tilt the row, or remove–reinsert one of its squares).  This costs almost nothing and removes the ε² degeneracy directly.  (A3, A4)
2. **Replica exchange in side length** (Odriozola-style) over 8–16 rungs from s_target to about k+0.3.  Rows are allowed to break in loose rungs, and swaps bring the reorganized configurations back down.  It uses the cores we have.  (A1)
3. **Population search with splitting, or nested sampling, on "achieved side".**  Keep K configurations, cull the worst by feasible side, and clone and decorrelate survivors.  This amplifies rare families instead of hoping a single trajectory hits them, and gives P(side < s) as a free census.  (A2)
4. **Structure-aware recombination and seeding.**  Use cut-and-splice crossover between good packings plus population basin hopping with a dissimilarity term (the Addis/Locatelli/Schoen evidence at n ≤ 130).  Optionally **LLM-evolve the seed constructor** (ShinkaEvolve/OpenEvolve, with our engine as the evaluator).  Records are block/strip-structured.  (A4, A9)
5. **Minima-hopping feedback.**  Keep a fingerprint history; revisits raise the escape energy, and escape uses low-curvature MD directions (collective row tilts) rather than random kicks.  (A5)

Honourable mention: a rectangle-relaxation container (w×h, with |w−h| penalized gradually; A6).

## Concrete approximate-census recipe

1. **Fingerprint.**  For each quenched minimum at fixed s: canonicalize over the 8 symmetries of the square (choose the lexicographically smallest).  Sort the squares.  Hash rounded (x, y, θ mod 90°) at about 10⁻³.  Also record a **coarse class**: the row variable R, the number of tilted squares, and the multiset of contact-graph degrees.  Count both fine minima and coarse classes.
2. **Sample.**  Run M quenches (for example 10⁴–10⁵ per n) from the production start protocol at 2–3 values of s (k−0.05, s_rec+0.02, …).  Record the abundance table: F_j = the number of minima seen exactly j times, and S_obs = the number of distinct minima.
3. **Estimate:**
   - The **missing mass** m₀ ≈ F₁/M is the probability that the next quench is new.  If m₀ ≈ 1, the census is hopeless at the fine level, so move to coarse classes.
   - **Chao1:** S ≳ S_obs + F₁²/(2F₂).  Report it as a lower bound, and **ACE** alongside it.
   - **Fraction attracted to rows:** P̂(R = k), with a binomial confidence interval.  This is the headline number for the failure mode.
   - **Record-class probability:** if it is never hit, the 95% upper bound is 3/M (rule of three).  Refine it by splitting (B3) or by the ratio of basin volumes from a known record.
4. **Basin-size check on 3–5 representatives** (record, best near-miss, typical row-jam).  Use the return-probability profile p(r), or MBAR with harmonic restraints if it is affordable using the pele/mcpele approach (reimplemented in our engine, not run from their code).  log(V_record/V_row) predicts how many quenches each record costs.
5. **Scaling.**  Repeat for several n (for example 50, 80, 120, 200).  Fit log P(row) and log P(record-class) against n.  Exponential decay of P(record) with a measured rate tells us how much compute each n deserves, and which method (2 or 3 above) must provide the amplification.
