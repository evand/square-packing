# Gensane & Ryckelynck 2005: the billiard/inflation method for squares in a square (2026-10-12)

T. Gensane, Ph. Ryckelynck, "Improved dense packings of congruent squares in a square", *Discrete Comput. Geom.* 34
(2005) 97–109, doi:10.1007/s00454-004-1129-z (received April 2004, online October 2004).  Read in full from the
user-supplied PDF (`../private/literature/gensane2004.pdf`, converted with `pdftotext -layout` in the session
scratchpad; not copied here).  Before this we only knew the paper secondhand (`search/packer/lit-optimizers.md` §1,
`lit-provenance.md`).

**Labels.**  **[read]**: from the full text, with section / proposition / procedure numbers.  **[verified]**: we
checked it ourselves (store, catalogue snapshot, a script, or a bibliographic record).  **[unverified]** /
**[inference]**: from memory or our own reasoning.

**Bottom line.**
- The method is a hard-constraint, maximin, stochastic "billiard".  Squares only ever grow: the size is always the
  largest inflation the current centres and angles allow.  Moves are random single-square moves, accepted only if
  feasible at the current size.  A jitter-everything, deflate, re-inflate, keep-if-better outer loop is
  **monotonic basin hopping** with an adaptive amplitude.  There is no overlap energy, no gradient, and no jam
  detection.
- They report **our failure mode**: the walk is "attracted" to configurations in which all squares share one angle,
  and these are "rarely good" (§6).  Their remedy is volume: a few thousand random starts, then cherry-picking.
- Every ingredient is something we already have, mostly in a stronger form.  Nothing in the paper targets the k-row
  jam.  Its value for us is (a) the lineage: Schadt's 2025 SA, which found s(29) and other records, began as a
  reimplementation of this method; (b) a closed-form pair "inflation distance" ψ; (c) a calibration target (n = 29, 37).
- Records: s(29) = 5.934342 stood for 21 years and fell in Dec 2025 (Schadt).  The s(37) improvement was already beaten
  by Cantrell's unpublished Sep 2002 packing.  None of their values is a record today.

---

## 1. Method [read]

### 1a. Objective: maximal inflation ω(C)
- Container [−L, L]² (L = 1 in the code).  The square q_{a,b,θ,c} has centre (a, b), angle θ, half-diagonal c, and
  vertices (a, b) + c(cos θᵢ, sin θᵢ), θᵢ = θ + iπ/2 (§2).
- **Angle convention.**  θ = 0 is a *diamond* relative to the container (|x| + |y| ≤ c).  So the paper's angles
  are ours + 45° (mod 90°).  [verified]: the 14 squares of their s(29) that they list at "45°" are the 14 axis-aligned
  squares of the catalogue copy (store id 211958; tilts 64.01°, 68.38°, 69.74°×3, 69.75°×3, 71.82°×7 = their
  19.01°, 23.37°, 24.67/24.74°, 26.82° plus 45°).  Read §6–7 with this in mind.
- **Pair inflation ψ(q′, q″)** (§3, Prop. 1).  This is the supremum k at which the two squares, scaled about their
  own centres to half-diagonal k, are still disjoint.  It is a growth distance: invariant under rigid motions and
  under common scaling (Lemma 1).  Closed form: in the frame of one square, ψ₀(a, b, θ) = (|a| + |b|) · min over the
  positive denominators of 1 / (1 − (sgn a cos θᵢ + sgn b sin θᵢ)), and ψ = max of ψ₀ evaluated in each square's frame.
  Parallel squares give (|a| + |b|)/2.  The non-overlap condition is ψ > c (Lemma 1(b)).
  - [verified] We checked this formula against bisection with an SAT overlap test on 3000 random pairs: max relative
    error 3e-15 (scratch script).
  - Small erratum: Remark 3 claims ψ(q₁,q₂) + ψ(q₂,q₃) < 5.08, but the sum is 0.5858 + 4.5 = 5.086.  The intended
    point, that the triangle inequality fails against 5.858, stands.
- **Wall inflation ϕ** (§4, Prop. 2): ϕ(a, b, θ) = (L − max(|a|, |b|)) / max(|cos θ|, |sin θ|).
- **ω(C) = min(minᵢ<ⱼ ψ, minᵢ ϕ)** (Prop. 3).  Rescaling all squares to half-diagonal ω(C) is admissible, and the
  global max of ω equals √2 L / s_n.  It is a nonsmooth maximin, the square analogue of the min-distance in
  disk "maximal separation" (§1).
- **Contact graph** (§5).  A bipartite graph of centres and contact vertices.  A square is "isolated" (a rattler)
  iff it has no edge (Prop. 4).  The graph is used for visualisation only, not in the search.

### 1b. Algorithm (§6, Procedures 1–4)
- **P1 RandomWalking(C, Nₐ, ε, c).**  Repeat Nₐ times:
  - pick a square (random or sequential);
  - propose a centre uniform in the ε-ball (clipped to the box) and θ uniform in (θᵢ − ε, θᵢ + ε);
  - accept iff min over j of ψ(z, qⱼ) ≥ c and ϕ(z) ≥ c, i.e. still feasible at the *current* size.
  Neutral moves are accepted, so non-bottleneck squares and rattlers diffuse freely.  This is hard-particle MC at
  infinite pressure.  The same ε is used for position (box half-width 1) and angle (radians).
- **P2 BilliardOfSquares(C, ε₁, ε₂, Nₐ).**  ε = ε₁, c = ω(C).  While ε > ε₂: run P1; if ω(C) > c, set
  c := ω(C) (inflate) and ε := 2ε, else ε := ε/2.  So inflation is implicit, after every walk.  There is no
  schedule, and the amplitude follows a 1/1 success rule.
- **P3 Perturbation(C, ε).**  Move *every* square uniformly within ε (centre and angle), then rescale all squares to
  ω(C).  This makes any overlap feasible by **deflating**, not by penalty resolution.
- **P4 WithPerturbations(C, ε₁, ε₂, factor, Nₐ).**  While ε > ε₂:
  - perturb (P3), then run a short billiard P2(C, ε, ε/factor, Nₐ);
  - if ω improved, keep it and set ε := 2ε; otherwise restore the incumbent and set ε := ε/2.
  This is monotonic basin hopping (accept only improvements) whose kick size adapts by the same 1/1 rule.
- **Stopping** is purely by amplitude: stop when ε < ε₂.  Nothing is time- or stall-based.
- **Parameters used.**
  - Screening: thousands of P2(C, 0.1, 1e−8, 1000) runs from random configurations.  They cherry-pick good C* from
    these.
  - The new results all come from P4(C*, 0.1, 1e−12, 1.5, 1000).  With factor 1.5, each inner billiard is a single
    walk, which repeats only while it keeps succeeding.
- **Implementation.**  C on an 800 MHz CPU.  **No run times, no success rates, no counts of restarts** beyond "some
  thousand", and no code.  The only success statements are qualitative: n = 11 "several times"; for multi-angle
  optima such as n = 17 the chance is "weak".
- **Finish.**  For n = 11 and 17 they derive exact polynomial systems (Maple / 20-digit solve).  For 29 and 37 the
  systems were "too complex"; they checked the bounds "with sufficient accuracy" (§7).

### 1c. Jams and how they escape (§6) [read]
- Quote: the billiard "often … gives jammed packings with ω(C) close to √2/(√n + 1) when the configuration obtained has
  just the angle θ = 0".  The procedures "seem to be 'attracted' by configurations with angle θ = 0 which are rarely
  good".  √2/(√n + 1) means s ≈ √n + 1.
- **[inference]** on what this attractor is.  In their convention θ = 0 is literally the 45° diamond, but parallel
  squares have ψ = (|a| + |b|)/2 whatever their common angle.  The sentence most plausibly means "all squares share
  one angle", i.e. a lattice-like jam.  For axis-aligned squares that is our grid / k-row jam (s = ⌈√n⌉ ≲ √n + 1).
  Thomas Schadt's s(55) notes (`search/packer/lit-provenance.md`) support this: his G&R reimplementation stalled at
  7.9577 against the trivial 8, and his SA "almost always gets stuck just above the trivial size … due to stacked
  rows and/or columns".
- **Escape mechanism.**  None is specific to jams: they rely on restart volume, cherry-picking, and global jitter plus
  deflation (P3–P4).  There is no move that breaks rows, no relocation, and no topology change beyond what jitter
  produces.

### 1d. Structural observations (§7) [read]
- n = 11: five central squares share one tilt.  cos θ is a root of a degree-8 polynomial over Q(√2), and a relation
  ties s to (cos θ, sin θ).  [verified]: our store's Trump s(11) (id 3) gives P(cos 4.818°) = −3.6e-11 (their
  θ = ours − 45°), and their s relation reproduces 3.8770836.
- n = 17 (Bidwell): two tilt angles; polynomial system of degree 7.
- n = 18: alternative optimal packing at (7 + √7)/2 with six isolated squares and six more movable by small
  translations.  Optimal packings are highly degenerate.
- n = 29: six distinct angles (19.01°, 23.37°, 24.67°×3, 24.74°×3, 26.82°×7, axis×14).  "Impossible … without a
  computer-aided method."  A second packing (Fig. 7) is at s ≈ 5.958702.
- n = 37: five distinct angles, six isolated squares.
- **Pattern** [inference]: their wins are multi-angle, rattler-rich, near-degenerate configurations, made of a block
  of axis-aligned squares and a band of nearly but not exactly parallel tilted squares.  This is the "chaotic"
  structure that SA later found for the n² − n − 1 family (lit-provenance).

## 2. Results and current status

| n | before (as cited) | G&R 2004 | now (store `pk.py ls`, 10-12) | status [verified] |
|---|---|---|---|---|
| 11 | ≤ 3.8772 (Friedman's rounded Trump bound) | 3.87708359 (same Trump packing, exact degree-8 form) | 3.877083590023 Trump | Not a new packing, only sharper digits.  The abstract's "improve n = 11" means the printed bound. |
| 17 | 4.6755 (Bidwell, rounded) | "4.6755300960455" | 4.67553009360455 Bidwell | Same packing.  Their printed value drops a digit ("…0093604…" → "…009604…"); it is 2.4e-9 too high. |
| 18 | (7+√7)/2 (Hämäläinen, Gustafsson) | alternative packing, same side | catalogue `square-18d` credits G&R | An alternative optimum, not an improvement. |
| 29 | 5.9648 (Bidwell) | **5.934342** (−3.05e-2); second packing 5.958702 | 5.933833462677 (Schadt 8 Dec 2025, analytic by Ellsworth 10 Dec 2025) | Stood 21 years.  Now 5.08e-4 above the record.  Schadt's repo README calls G&R's the previous record. |
| 37 | 6.6213 (Friedman) | **6.603236** (−1.8e-2) | 6.598619609244 (Cantrell, **Sep 2002**) | Already beaten when published: Cantrell's packing predates it but first appeared in Friedman's survey v4 (Nov 2005).  G&R's packing is 4.61e-3 above the record (catalogue `square-37_r3`, re-solved by Ellsworth 2023 to 6.6032338). |

Sources: the store's kingbird snapshot `search/packer/runs/kingbird/2026-10-10/manifest.json`, the register rows
via `pk.py frontier`, and Schadt's repo github.com/BalthasarStrauss/Squares-packing_S-29-_New-Record.  That repo holds
coordinates, a checker and a README only.  Its "Methodology" says only "C++ code with data type float", relaxed to
1e-100.  There is no search code; we did not execute anything from it.

## 3. Comparison with our engine

"Ours" = `packer` (penalty L-BFGS), `fq` (lifted separating-line ALM + face-branch SLP polish), `hop` / explorer,
`anneal.rs` (hard rounded-square NPT MC), `melt`, `hhop` (event chains), `fq --shrink` (see `search/packer/PACKER.md`).

| Item | G&R 2005 | Ours | Verdict |
|---|---|---|---|
| Objective | ω(C), nonsmooth maximin (inflation in a fixed box) | `fq`: min s with smooth lifted separators, then SLP to an exact KKT point; `packer`: fixed-side penalty | Same problem; our local solve is far sharper (exact contacts versus a random walk to ε = 1e−12). |
| Feasibility | Always feasible; overlaps never exist (P3 deflates) | Penalty / ALM phases pass through infeasible states; `repair` / `fq --shrink` grow to feasible | Deflate = uniform box growth, which we have (`--shrink`, repair).  Not new. |
| Local moves | Single-square uniform ε-ball translate+rotate, accepted iff feasible at current size (neutral moves accepted) | `anneal.rs` hard-square NPT MC (finite pressure, adaptive steps); `melt` (local hard MC); event chains (`hhop`) | We have it.  Our "squares-only MC" control (= G&R's walk at finite pressure) crystallised to the grid 17–24/24 at k² + m sizes: the same attraction G&R report. |
| Rotations | Same ε for angle (rad) and position; θ mod π/2; no special handling | Continuous θ in all engines; `fq` ties separator to owner face | Nothing to take. |
| Inflation schedule | Implicit: c := ω after every walk; no continuation | Squeeze / soft squeeze (μ ladder), `fq` min s, shrink-hopping | Ours has more control. |
| Step adaptation | 1/1 rule: ε×2 on improvement, ε/2 on failure, stop at ε₂ | hop kick σ drawn from fixed ranges (0.003–0.03, local 0.05–0.2); bandit over arms in explorer | **We lack an amplitude-adaptation rule inside an arm.**  Cheap to add. |
| Perturbation | Jitter *all* squares by ε, deflate to feasible, short re-inflate | Global and local kicks, melt, rowslide, relocate, cross, chain seeds | Ours is a superset. |
| Acceptance | Monotone (improve-only) MBH | Metropolis T 1e−5, reset to best; delayed acceptance; archive | Ours is a superset. |
| Restarts | Thousands of short random billiards → cherry-pick → long MBH | Census, explorer archive, chain / neighbour seeds, register seeds | Ours is a superset.  Seeding from neighbours is what wins at n ≥ 50 (SQUISH, Schadt). |
| Jam escape | None specific; admitted (§6) | Grid-obstructed screens, row-jam detector, rowslide, soft squares (pending test) | Neither solves the k-row jam; G&R add nothing here. |
| Rattlers / degeneracy | Notes isolated squares (18, 37) and flat motions; the walk lets rattlers diffuse | exactsolve rattler check, rattler-masked distances; rattlers are frozen during polish | Free rattler diffusion *between* polishes is the one dynamic we don't run on polished packings (see W2). |
| Exactness | Maple / 20-digit polynomial systems for 11, 17 | `search/exact/` exact forms, certified bundles | Ours is stronger. |
| Pair distance | Closed-form ψ (growth distance, scale-invariant) | SAT penetration depth; separator LP gaps | ψ is a cheap exact "how much can this pair grow" (W4). |

## 4. Worth trying (ranked; ≤ 4)

Nothing here targets the k-row jam directly.  G&R suffered from it and offered only restart volume.  If anything
breaks the jam it will be the soft-square / transversal-topology work, not this paper.

1. **Calibration at n = 29 and 37 (and 18) from random**.  Our calibration stops at 11 and 17.  Run `fq` hop or
   explorer from random starts at 29 and 37 (which have 5–6 distinct tilt angles and 6 rattlers) and record time to
   reach G&R's level (5.934342 / 6.603236) and the record (5.933833 / 6.598620).
   - *Failure mode:* multi-angle optima, which G&R say random billiards find only weakly (n = 17), and the
     one-angle attractor at small n.  It also gives a direct "do we beat 2004" baseline, and 29 is a record that
     Schadt's SA reached from randomness.
   - *Cost:* about 1–2 core-hours (20 starts × 2 n × a few minutes), with a script outside `search/packer`.
     Pre-register in PACKER.md.
2. **Rattler diffusion between polishes**: G&R's P1 on a polished packing.  At fixed s, run hard single-square moves
   accepted iff overlap-free (neutral moves allowed, about 10³·n proposals, ε adapted by the 1/1 rule), then `fq`
   polish.
   - *Failure mode:* degenerate, rattler-rich minima (our sub-10 n = 90 bests have up to 8 rattlers and 20–45 flat
     motions; G&R's 18 / 37 likewise) and vacancies that never propagate (hop1: "the vacancy is never absorbed").
     Free volume moves where gradient polish sees zero force.
   - *Cost:* about 1 day.  The overlap test exists in `anneal.rs` (r = 0), so this is a new explorer arm graded with
     `pk/battery.py`.
   - Expectation is modest [inference]: `melt` already does a local finite-pressure version.  The difference is
     that this one is global, at infinite pressure, on the polished state, and does not crystallise.
3. **1/1 adaptive kick amplitude inside hop / explorer arms**: σ ×2 after an accepted improvement, ÷2 after a
   failure, reset to the arm's range at a floor.  This replaces fixed σ ranges.
   - *Failure mode:* proposals wasted at the wrong scale.  hop0 found nothing from records at fixed σ ranges; polish
     dominates cost, so each wasted kick costs a full quench.
   - *Cost:* a few hours plus one battery run (about 2–4 core-hours).  It overlaps with the bandit, so measure
     against it.
4. **Closed-form ψ as a cheap feasibility / growth screen**: ω(C) in O(pairs) with no LP.
   - *Failure mode:* screening cost.  It gives exact feasibility after kicks or in the P1 walk (needed for item 2),
     and a scale-free gap that could rank which pairs bind ("pain" for relocate moves).
   - *Cost:* about 2 hours, since the formula has been checked to 3e-15 against bisection.  Low priority unless item
     2 is built.

Not worth it: G&R as a whole-pipeline start generator.  Our `anneal.rs` squares-only control is essentially that and
crystallises, and G&R report the same.

## 5. Later papers citing it (OpenAlex 13 / Semantic Scholar 18 citations, 2026-10-12) [verified bibliographic records]

**None improves their packings in print.**  The improvements came through Friedman's survey and Ellsworth's
catalogue (Cantrell 2002, Schadt 2025).  The five most relevant:

1. E. Friedman, "Packing unit squares in squares: a survey and new results", *EJC* Dynamic Survey DS7 (2009 version),
   doi:10.37236/28.  This records the G&R entries (n = 29, 37, 18), and Cantrell's s(37), which supersedes theirs.
2. T. Gensane, P. Ryckelynck, "Producing dense packings of cubes", *Discrete Math.* 308 (2008) 5230–5245,
   doi:10.1016/j.disc.2007.09.047.  The 3-D sequel (abstract not retrieved; whether the algorithm changed is
   [unverified]).
3. T. Gensane, P. Ryckelynck, "On the maximal inflation of two squares" (2005; CiteSeerX 10.1.1.88.4966, no DOI).
   A companion note on ψ (not read).
4. M. Abrahamsen, T. Miltzow, N. Seiferth, "Framework for ∃R-completeness of two-dimensional packing problems",
   FOCS 2020, doi:10.1109/FOCS46700.2020.00098; journal version *TheoretiCS* (2024), doi:10.46298/theoretics.24.11.
   Packing with rotations is ∃R-complete, which is the theoretical reason exact forms get hard.  [unverified]: what
   they cite G&R for.
5. M. Abrahamsen, J. Stade, "Hardness of packing, covering and partitioning simple polygons with unit squares", FOCS
   2024, doi:10.1109/FOCS61266.2024.00087; *SIAM J. Comput.* (2025), doi:10.1137/24M171824X.

Also citing it, but less relevant: Kamali & Nikbakht, cutting stock with rotation (COCOA 2020,
doi:10.1007/978-3-030-64843-5_36); Imamichi's NLP-for-packing thesis (Kyoto 2009, doi:10.14989/DOCTOR.K14773); an
arXiv preprint "A counterexample to Nagamochi's scoring lemma and a new rectangle packing bound" (2026, not read);
the CG:SHOP 2024 report (arXiv:2403.16203); and ∃R neural-network papers.  Secondary: H. Gill's 2026 essay (in
lit-optimizers §1) calls Schadt's program a "continuation" of the G&R algorithm.
