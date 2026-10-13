# Irregular / polygon packing papers from the group-chat list, read against our packer (2026-10-12)

Follow-up to [`group-chat-2026-10-12.md`](group-chat-2026-10-12.md).  It adds to `search/packer/lit-optimizers.md` and
`search/packer/lit-2026-10-12.md` and does not repeat them.

**Method.**  Bibliographic records came from the Crossref, Semantic Scholar, OpenAlex, arXiv and DOAJ APIs.  Abstracts
came from those APIs, publisher pages or web-search summaries of them.  Five open-access PDFs were converted with
`pdftotext` in the scratchpad and their method and results sections were read: Peralta et al. 2017/2018 (two papers),
Birgin et al. 2006, Rajagopalan et al. 2018, Callisaya & Kari 2025, plus the Guo et al. review.  All the Springer
(J. Glob. Optim., Wireless Netw., Mobile Netw. Appl.) and Elsevier full texts were paywalled or blocked by a bot check,
so for those items only the abstract was read.  No downloaded code was executed.

*Updated 2026-10-12 after reading full text:* Evan later supplied four of the paywalled papers (§6, §8, §12a, §12b).
Their method and results sections were read and those sections are now labelled [read], with section and equation
references.  The verdicts did not change, but several details did: Cherri's binaries pick an *edge* (not a free line),
LOFRT is a per-item trust box plus pair pruning, and the global-solver size limit is lower than the scan implied
(see §8 and Ranked item 2).

**Sources read in full** (method and results sections; DOIs only):
- Stoyan, Pankratov, Romanova, quasi-phi-functions and ellipses, J. Glob. Optim. 65 (2016), doi:10.1007/s10898-015-0331-2.
- Stoyan, Pankratov, Romanova, irregular objects with continuous rotations, J. Oper. Res. Soc. 67(5) (2016),
  doi:10.1057/jors.2015.94.
- Cherri, Cherri, Soler, MIQCP with continuous rotations, J. Glob. Optim. 72 (2018), doi:10.1007/s10898-018-0638-x.
- Kallrath, Romanova, Pankratov, Litvinchev, Infante, minimum-perimeter convex hulls, J. Glob. Optim. 85 (2023),
  doi:10.1007/s10898-022-01194-4.
- Earlier in this scan (open access): Peralta et al. (arXiv:1707.07177; doi:10.5220/0006602700710077), Birgin et al.
  (doi:10.1057/palgrave.jors.2602067, author preprint), Rajagopalan et al. (doi:10.56801/rebicte.v4i.77), Callisaya &
  Kari (doi:10.5540/tcam.2025.026.e01823), Guo et al. (doi:10.3389/fmech.2022.966691).

*Most useful:* Cherri 2018, for its hard numbers on how small a global solve must be and its symmetry-breaking
constraints for identical pieces.  Second, the quasi-phi paper, for LOFRT, the published form of our neighbour list.
The other two add nothing.

**Labels.**
- **[verified]**: the record or abstract was fetched this session.
- **[read]**: I also read the method or results section in the full text.
- **[unverified]**: from memory or a secondhand summary.

**Our baseline, for comparison.**
- *Engine:* fixed side; SAT penetration² plus wall energy; L-BFGS; basin hopping; squeeze.
- *`fq`:* one separating line (φ, d) per near pair, written as smooth constraints and solved with a PHR augmented
  Lagrangian.  Separators are snapped to faces, ambiguous pairs are flipped, and a trust-region SLP polish follows.
  The objective is s, so the container is already a variable.
- *Exact layer:* the incidence MILP jam test and KKT Newton in mpmath.

**One-paragraph verdict.**  *(Updated 2026-10-12 after reading full text: Cherri moved to the edge-selection group;
the item range was widened.)*  Every one of the ten "likely relevant" papers uses a non-overlap formulation we already
have:
- a per-pair free separating line (Kallrath 2009, Peralta 2017/18, Callisaya 2025, the quasi-phi model of Kallrath
  2023 §4.2; the same object as Berthold et al.);
- an edge-normal max–min phi-function, or a binary choice of separating edge (Stoyan 1996, Rajagopalan 2018, the
  phi model of Kallrath 2023 §4.1, Cherri 2018).  This is SAT, and it is our incidence model;
- a smooth point-in-polygon penalty that is weaker than SAT (Birgin 2006).

They solve these with IPOPT, GENCAN, BARON or SCIP, from multistart or bottom-left starts.  Global solves cover 2 to 6
items; local solves go up to 300.  None
mentions jamming, lattice attractors or anything like the k-row problem.  Their "variable container" is what `fq`
already does, since the side is the objective.  Only one item reports squares-in-a-square: Birgin et al. 2006 reproduce
the then-known packings for n ≤ 37 at fixed sides.

The only idea we have not tried is the **soft-item homotopy** (Litvinchev 2024; Bennell 2024; Melashenko 2025): items
change aspect ratio at fixed area and are then made rigid again.  It is an untested speculative continuation for the
k-row jam.  Those papers show only toy examples.

---

## 1. Berthold, Kamp, Mexi, Pokutta, Polik: "Out-of-the-box global optimization for packing problems"

**Already covered.**
- *Record:* arXiv:2605.04850 (2026), noted in `lit-optimizers.md` §1 and REVIEW-2026-10-11.
- *What we use:* `fq`'s lifted separator comes from this paper.
- *Not re-read.*

## 2. Stoyan, Novozhilova, Kartashov: "Mathematical model and method of searching for a local extremum for the non-convex oriented polygons allocation problem"

- *Record:* *Eur. J. Oper. Res.* 92(1) (1996) 193–210, doi:10.1016/0377-2217(95)00038-0 [verified].
- *Abstract:* from the IAOR summary only; there is no abstract in Crossref or OpenAlex, and the ScienceDirect page is
  paywalled.  Irregular objects are placed on a strip of given width that has defects, in fixed orientation
  ("oriented": no rotation).  The model uses only continuous variables, built on Φ-functions and systems of linear
  inequalities.  A local minimum is found from any feasible start using active-set ideas [verified].
- *Method:* an early phi-function active-set local search for *translations only*.  The later continuous-rotation
  version is Chernov–Stoyan–Romanova 2010 (`lit-optimizers.md` §2).
- *For us:* **nothing new.**  Our slp2 and `fq` SLP polish is an active-set method on face incidences with rotations.

## 3. Rajagopalan (Narendran R.), Mala C., Nobert: "An Optimization of Phi-function for Convex Polygons"

- *Record:* *Research Briefs on Information & Communication Technology Evolution* (ReBICTE) 4 (2018) art. 17,
  pp. 158–165, doi:10.56801/rebicte.v4i.77 [verified, read].
- *Method:*
  - The standard convex-polygon phi-function is max( max_i min_j u_ij, max_j min_i v_ji ), where u_ij is the signed
    distance from edge i of one polygon to vertex j of the other.  This is the SAT separation over edge normals.
  - The paper shows by a two-triangle example that this value is not the true Euclidean separation when the polygons
    are apart: different gaps can give the same value.  It proposes a corrected formula.
  - There is no packing algorithm and no experiments.
- *For us:* **nothing new.**  We only need the sign and the overlap side, where SAT penetration is already exact for
  squares.  The "anomaly" concerns separated pairs, where our energy is zero anyway.

## 4. Bennell, Litvinchev, Pankratov, Romanova: "Packing stretched convex polygons in an optimized rectangle"

- *Record:* *Wireless Networks* 30 (2024) 7369–7376, doi:10.1007/s11276-023-03642-9 [verified].
- *Abstract:* from the White Rose eprint page [verified].
  - Polygons may be compressed or stretched along their principal axes within limits, at constant area, under free
    translations and rotations.
  - The objective is to minimise the height of a rectangular container.
  - Non-overlap and containment are written with phi-functions as a nonconvex NLP, and solved by multistart.
  - The paper gives "a numerical example".  The motivation is porous media under pressure.
- *For us:*
  - The formulation and solver are nothing new (phi-functions plus a local NLP, with the container height as a
    variable, as in `fq`).
  - The **new idea is the shape degree of freedom**: aspect ratio as a variable at fixed area.  See "Ranked" item 1.
  - The paper gives no evidence that this helps optimisation: there is one example, and the stretch is part of the
    problem, not a continuation.

## 5. Litvinchev, Infante, Romanova, Martinez-Noa, Gutierrez: "Packing Soft Convex Polygons in an Optimized Convex Container"

- *Record:* *Mobile Networks and Applications* 29 (2024) 211–220, doi:10.1007/s11036-023-02286-5 [verified].
  A conference version is "Optimized Packing Soft Convex Polygons", COMPSE 2022 (Springer 2024),
  doi:10.1007/978-3-031-34750-4_7 [verified listing].
- *Abstract:* via search summary of the Springer page [verified].
  - Soft convex polygons change shape within limits, staying convex and keeping their area.
  - The NLP has non-overlap, containment, convexity and area-conservation constraints.
  - The optimal container is computed for fixed elasticity.  An inverse problem asks for the smallest elasticity
    that still gives the same optimal container.
  - Experiments pack soft triangles and pentagons into optimized circular and square containers.
- *For us:* the same as §4.  The inverse problem ("least softness that reaches the same container") is the one
  conceptually interesting piece.  Read as a continuation, it asks how much shape freedom is needed to escape a given
  jam.  There are no rigid-limit results and no equal squares.  **No direct use**; it feeds "Ranked" item 1.

## 6. Kallrath, Romanova, Pankratov, Litvinchev, Infante: "Packing convex polygons in minimum-perimeter convex hulls"

- *Record:* *J. Glob. Optim.* 85(1) (2023) 39–59, doi:10.1007/s10898-022-01194-4 [verified].
*Updated 2026-10-12 after reading full text.*
- *Problem (§2):* convex polygons, free translation and rotation, are packed into a convex m-gon of minimum perimeter.
  The container is a closed chain of m segments, each with a start vertex, angle and length (eqs. 1–2).  With m ≥ the
  total vertex count, this is the minimum-perimeter convex hull [read].
- *Two formulations (§3–4) [read]:*
  - the classical phi-function, max(max_i min_j, max_j min_i) of edge–vertex terms (eq. 5): SAT over edge normals.
    The nonsmooth feasible set is a union of η smooth subregions, and global optimality means solving all of them
    (eqs. 10–11; η = m_A + m_B for two polygons);
  - the quasi-phi-function (eq. 6): one free separating line (φ, γ) per pair, with the normal written as
    (cos φ, sin φ).  It gives σ = 4m + 3n + n(n−1) variables and O(n²) smooth constraints, all pairs, no neighbour list.
- *Results (§6) [read]:*
  - Examples 1–6 have 2, 2, 2, 3, 4 and **6** polygons.  They are "solved to optimality" by BARON on the phi-function
    model under a 48 h limit.  **No run times, gaps or node counts are reported.**  The 6-item case is six identical
    right triangles (0,0),(0,4),(2,4).
  - Examples 7–8 (10 pentagons; 13 trapezoids plus 12 triangles = 25) are IPOPT local optima on the quasi-phi model,
    with a 12 min limit and a single reported value each.
  - The appendix only conjectures the two-triangle optimum.
- *For us:* **nothing new; confirmed.**  "Up to 6" is right, but it is weak evidence: six copies of one small triangle,
  with no timing.  Cherri 2018 (§8) gives better-calibrated numbers.  A perimeter or hull objective is not our problem.

## 7. Kallrath: "Cutting circles and polygons from area-minimizing rectangles"

- *Record:* *J. Glob. Optim.* 43 (2009) 299–328, doi:10.1007/s10898-007-9274-6.  Online 2008, issue 2009 [verified].
- *Abstract:* via search summary of the Springer page [verified].
  - Circles, rectangles and convex polygons, at any orientation, are cut from one or several rectangles whose area
    is minimised; the rectangles have bounded width and length.
  - It is a nonconvex NLP, or an MINLP when there are several plates.
  - **Separating hyperplanes** give the non-overlap.  Symmetry breaking is used.  Lower bounds come from replacing
    each polygon by its largest inscribed circle.
  - It finds globally optimal solutions "for small numbers of objects".
  - The test instances are c1p1 … c6p3, i.e. at most about 9 objects; Peralta 2018 (§9) reuses them [read, via §9].
- *For us:* **nothing new.**  It is an ancestor of Berthold's separator.  The inscribed-circle lower bound is far too
  weak for s(n).  The polylithic follow-up (Kallrath–Rebennack 2014) is already in `lit-2026-10-12.md` §4.

## 8. Cherri, Cherri, Soler: "Mixed integer quadratically-constrained programming model to solve the irregular strip packing problem with continuous rotations"

- *Record:* *J. Glob. Optim.* 72(1) (2018) 89–107, doi:10.1007/s10898-018-0638-x [verified].
*Updated 2026-10-12 after reading full text.*
- *Model (TRPM, §3) [read]:*
  - Non-convex pieces are split into convex parts (Greene's decomposition).  Piece i has a reference point (x_i, y_i)
    and rotation variables (s_i, c_i) with s_i² + c_i² = 1 (eq. 3); vertices are bilinear in them (eqs. 4–5).
  - Non-overlap uses the D-function (eq. 2): the sign of a vertex relative to an oriented edge.  For each pair of
    convex parts there is one binary γ per edge **of either part**.  γ = 1 forces every vertex of the other part to the
    outer side of that edge by a big-M constraint (eq. 9).  **Exactly one** γ per part pair is 1 (eq. 10).
  - So the binary really picks a separating edge: SAT restricted to edge normals.  There is no free line variable, so
    this is the edge-branch model, not Berthold's lifted line.
  - It is quadratic, not trigonometric: the D-function of rotated points is bilinear in (s, c, x, y).
- *Two variants (§4.2, §5.1):* they differ only in big-M, M = d_ll'(H + L), with L either the variable length
  (TRPM_L) or a constant upper bound (TRPM_E).  The scan's guess that they "differ in how non-overlap is written" was
  wrong.
- *Symmetry breaking (§4.3):*
  - a piece with rotational symmetry angle φ gets c_i ∈ [cos(φ/2), 1], i.e. θ ∈ [−π/4, π/4] for a square;
  - identical pieces are ordered, x_i ≤ x_j;
  - otherwise one piece's y is restricted to [0, H/2] (reflection).
  - The exactly-one constraint (10) also removes duplicate edge choices.
- *Results (§5; GAMS, i7-4790, 7 h limit; COUENNE 0.5, SCIP 3.2, BARON 16.5):*
  - **Proved optimal:** "Three" (3 pieces) by SCIP in 34 s.  Jones' four-part (4 convex pieces that tile the board
    exactly) by SCIP in 17–64 s; Jones' QP-Nest needed 20 807 s.  COUENNE and BARON did not close four-part in 7 h.
  - **Not closed:** Jones' six-part (6 convex pieces, perfect fit).  Best is UB 12.000 (the given start) against LB
    7.903, and SCIP ran out of memory.  Poly1A subsets of 5 pieces: SCIP OOM after about 3.4 h with gaps of 50–80 %.
    Blaze2[1] and Shapes4 timed out.
  - Warm start: three back-to-back 7 h SCIP runs on Shapes4, each seeded with the previous result, reduced the length
    from 23.02 to 19.92.
- *For us:* **nothing new for the engine; a useful calibration.**
  - The edge binary is exactly our incidence MILP (`inc.py`) and `fq`'s face-branch flips.
  - A *whole-instance* global MIQCP closes at 3–4 pieces and fails at 5–6, even with exact-fit instances, warm
    starts and symmetry breaking.  That is below the "≤ 6–10" assumed for window solves (Ranked item 2).
  - Its symmetry constraints (square θ ∈ [−π/4, π/4], identical-piece ordering, exactly-one edge) are the obvious
    ones to put into any window MINLP of ours.

## 9. Peralta, Andretta, Oliveira: separation lines (two papers)

**(a) "Solving irregular strip packing problems with free rotations using separation lines"**
- *Record:* *Pesquisa Operacional* 38(2) (2018) 195–214, doi:10.1590/0101-7438.2018.038.02.0195, arXiv:1707.07177
  [verified, read].

**(b) "Packing Circles and Irregular Polygons using Separation Lines"**
- *Record:* ICORES 2018, pp. 71–77, doi:10.5220/0006602700710077 [verified, read in part].

**Method.**
- Each pair of convex parts has a separating line written as y = c·x + d.
- The line passes through a reference point (x̄, ȳ) at angle α, so there are 3 variables per line.  This is fewer than
  Kallrath's vector form (about 87 → 10 variables for two polygons).
- All vertices of one part lie on one side and all of the other's on the other side.
- Lines are created for **all pairs**: there is no neighbour list, so there are O(n²) lines.
- Rotations are free (vertex = R(θ)·local + t).  The strip length (in b, the rectangle area) is a variable.
- The solver is IPOPT with an L-BFGS Hessian approximation and a 1 h cap.  Starts come from bottom-left placement with
  no-fit raster at 4 orientations, best of 1000.

**Results.**
- In (a): ESICUP instances with 22–440 convex parts; poly20a has 290k variables and 717k constraints.  It is slightly
  worse than Stoyan–Pankratov–Romanova 2016 (§12a) and better than the "rubber band" method.
- In (b): Kallrath's 2009 circle-plus-polygon instances (at most about 9 objects), equal or better areas.

**For us: nothing new.**
- This is Berthold's or `fq`'s separator in a worse chart: the slope form is singular for vertical lines, while ours
  uses an angle.
- All-pairs lines are the cost wall that our neighbour list avoids.
- Results are local optima from heuristic starts, with no squares.

## 10. Birgin, Martínez, Mascarenhas, Ronconi: "Method of sentinels for packing items within arbitrary convex regions"

- *Record:* *J. Oper. Res. Soc.* 57 (2006) 735–746, doi:10.1057/palgrave.jors.2602067 [verified].
- *Full text:* the author preprint (June 2005), https://www.ime.usp.br/~egbirgin/publications/bmmr.pdf [read].

**Method.**
- *Sentinels.*  Each item gets a finite point set (sentinels) such that two copies overlap only if some sentinel of one
  lies strictly inside the other.  Finite sets exist when all interior angles are ≥ π/2.
- *Penalty.*  It is Σ_sentinels Π_edges max(0, b − ⟨c, ·⟩)², which is smooth.  Container violations are added as
  squared hinges on the vertices.  The whole thing is minimised unconstrained with GENCAN, a box-constrained
  active-set solver.
- *Minimal set for a unit square:* the 4 corners, the midpoints of two opposite sides, the centre, and two points on
  the other two sides.  That is 9 points.
- *Outer loop:* maximise the number of items m in a fixed region.  Add one item at a time, with random restarts until
  the penalty is zero.
- *Start trick.*  The authors note that the product penalty has "deep valleys" and many spurious minima.  So they
  first solve a disk-packing problem, with each rectangle approximately covered by a few circles, and use the result
  as the start.

**Squares-in-a-square (Table 2, problems 16–23).**
- Unit squares in squares of the then-best-known sides (Friedman 2000), rounded up: 2 + √2/2, 3 + √2/2, 3.878, 4.676,
  (7 + √7)/2, 3 + 4√2/3, 7/2 + 3√2/2, 6.621.
- It packs 5, 10, 11, 17, 18, 19, 26 and 37 squares respectively, i.e. it reproduces the known packings for n ≤ 37.
- Time: 0.02–2841 s on a 2005 Opteron (n = 37 is the slowest, 1253 restarts).

**For us: nothing new.**
- The sentinel penalty is strictly weaker than SAT penetration (it is zero only at the end, and has product valleys).
- Reproducing n ≤ 37 at a given side is what our engine already does from random starts (n = 11, 17 in seconds).
- The disk-cover pre-solve is the r = ½ end of the soft-to-rigid continuation we already have (`lit-2026-10-12.md` §6).

## 11. Other list items read for completeness

- **Guo, Zhang, Hu, Li, Wu, Peng, Zhang, "Two-dimensional irregular packing problems: A review."**
  - *Record:* *Frontiers in Mechanical Engineering* 8 (2022) 966691, doi:10.3389/fmech.2022.966691, open access
    [verified, read].
  - *Content:* an industrial nesting survey.  It covers the envelope-polygon and no-fit-polygon geometry, sequence
    metaheuristics (GA, SA, ACO, PSO), and a little on NLP models (it cites Peralta, Cherri and Stoyan 2016).
  - *Squares:* nothing on equal squares or free-rotation optimisation beyond those references.  The only square
    entry is Martello & Monaci 2015, "rectangles into the smallest square" (*Comput. Oper. Res.* 63), which allows
    only 90° rotation.
  - **Nothing for us**, and it is not even a useful index for continuous-rotation methods.
- **Caoduro & Sebő, "Packing, Hitting, and Coloring Squares."**
  - *Record:* arXiv:2206.02185 (2022) [verified; no journal version found].
  - *Content:* for a given finite family of squares in the plane, it bounds τ/ν (hitting number over packing number)
    by 6 for unit squares and 10 for general squares, not necessarily axis-parallel, and χ/ω by 6 and 9.
  - **Combinatorial geometry, not s(n).  Irrelevant.**
- **Callisaya & Kari, "Packing Polyhedra within Convex Sets."**  This is misfiled under 3D in the list: the experiments
  are 2D.
  - *Record:* *Trends Comput. Appl. Math.* 26 (2025) e01823, doi:10.5540/tcam.2025.026.e01823 [verified, read].
  - *Method:* a separating hyperplane (a, b) per pair with ‖a‖ = 1, or a = (cos η, sin η) in the "trigonometric
    model".  It is solved with IPOPT through JuMP, as a feasibility problem in an m → m + 1 loop with warm starts from
    the m − 1 solution.
  - *Result:* regular polygons in a circle of radius 4, including **38 squares** (density 0.74); there are also some
    3D icosahedra.
  - **Nothing new**: the same formulation as Berthold or Peralta, at toy size.

## 12. Related papers found while scanning (step 3)

**(a) Stoyan, Pankratov, Romanova, "Cutting and packing problems for irregular objects with continuous rotations:
mathematical modelling and non-linear optimization"**
*Updated 2026-10-12 after reading full text ((a) and (b)).*
- *Record:* *J. Oper. Res. Soc.* 67(5) (2016) 786–800, doi:10.1057/jors.2015.94 (advance online 16 Dec 2015)
  [verified record, read].
- *Model (§2) [read]:*
  - objects bounded by arcs and segments, decomposed into basic objects, with continuous rotations, in a rectangle
    (strip length or area) or a circle (radius);
  - **ordinary radical-free phi-functions, not quasi-phi**: no auxiliary variables, max–min of smooth functions
    (eq. 8);
  - the feasible set is a union of smooth subregions W_s (eqs. 9–13), each picking one max-branch per basic
    phi-function.
- *Feasible starts (§3.1) [read]:*
  - SPA1: a greedy ε-strip raster placement with discrete rotations.
  - SPA2: random centres for the circumscribed circles in an oversized container, then max λ (a common homothety
    0 → 1) subject to circle phi-functions, then random angles.  Since every object sits in its circle at λ = 1, the
    result is feasible.
- *Subregion procedure (§3.2–3.3, "IPLOA") [read]:*
  - At the current feasible point, each max-branch is fixed to the one with the **largest value** (§3.2).  For SAT
    pairs this is "snap the separator to the face with the largest gap".
  - Each vertex of the item's bounding rectangle may move only inside a δ-square around its current position (eqs.
    14–15, 16 nonlinear inequalities per item), with δ = the mean item size.
  - Pairs, or pairs of basic parts, whose δ-expanded rectangles cannot meet are dropped (Step 2, hierarchically).
  - IPOPT solves the smooth subproblem.  The branch choice and boxes are rebuilt at the new point, until the
    objective stops changing.
  - **There is no rule for ties or flips between branches** (cf. `fq`'s dual-ranked ambiguous flips).
- *Instances (§4) [read]:*
  - Burke/Whitwell Profile1–8 (better lengths than published, 1–220 s per local minimum);
  - Rocha's Poly10a and Poly20a, 150 and 300 objects (126.29 and 251.04 against 162.35 and 305.70; 619 s and 1209 s);
  - Leung's Swim (48 objects), Shapes2 and Dagli;
  - arc-approximated ellipses TC30/50/100 (slightly worse than (b), much faster);
  - new 52 to 160-object arc and ellipse instances.
  - **No rectangles-only or equal-square instances.**
- *Standing:* the reference for the phi-function school.  Peralta (§9) reports being slightly worse than it
  (not re-checked here).
- *For us:* **nothing new.**  IPLOA is our incidence-model polish (snap to the largest-gap face, re-snap each step,
  near pairs only) without the flip step.

**(b) Stoyan, Pankratov, Romanova, "Quasi-phi-functions and optimal packing of ellipses"**
- *Record:* *J. Glob. Optim.* 65 (2016) 283–307, doi:10.1007/s10898-015-0331-2 (accepted June 2015) [verified record,
  read].
- *Quasi-phi-function (Def. 4, Thm. 2, eqs. 1, 5) [read]:*
  - a function with extra variables u′ whose max over u′ is a phi-function.
  - **For convex polytopes the extra variables are exactly a separating line:** u_P = (θ_P, μ_P), with half-plane
    cos θ_P·x + sin θ_P·y + μ_P ≤ 0.  The constraint is min_i ψ_P(vertex of A) ≥ 0 and min_j −ψ_P(vertex of B) ≥ 0
    (eq. 5).  The normal is an angle, as in `fq`.
  - Appendix 1 notes that this needs one NLP where the max–min phi-function needs m₁ + m₂ branch NLPs.
  - **For ellipses (eq. 9) the extra variables are not a line:** they are two boundary parameters (t₁, t₂).  The
    constraint is "the tangent at t₂ on E₂ is anti-parallel to that at t₁ on E₁, and two points on E₂'s tangent lie
    outside E₁'s tangent".  A line (θ_P, μ_P) plus t's is used only with distance constraints (eq. 11).
- *Rotations:* a plain angle θ_i, entering through cos/sin of rotation matrices (§2, §5.1).  Homothety λ is a variable
  only in the start procedure.
- *Starts (SPA, §5.2.1):*
  - put tiny copies (λ₀ ≈ 0.01·min b / max a) at random non-overlapping centres in a container of side 2Σa_i, with
    random angles;
  - fit the auxiliary variables per pair by a small max-μ NLP;
  - maximise λ at the fixed container until λ = 1.
- *LOFRT (§5.2.2, App. 2–3) [read]:*
  - each ellipse gets a fixed axis-aligned box of half-side a_i + ε around its current centre, and only its centre
    is confined to that box.  ε = mean semi-minor axis was chosen from a time-vs-ε sweep on TC50; the time is flat for
    ε in [0.6, 1.6].
  - Pairs whose boxes do not overlap are dropped together with their auxiliary variables, and IPOPT solves the rest.
  - If the optimum hits a box face, the boxes are re-centred there, the dropped auxiliaries are re-fitted, and the
    solve is repeated.  If the optimum is interior, it stops.
  - This keeps O(n) pairs per subproblem instead of O(n²).
- *Compared with ours:*
  - The box is a hard per-item trust region.  It makes the pruning *rigorous* inside a subproblem: dropped pairs
    cannot collide.
  - `fq` uses a Verlet skin instead (cutoff √2 + skin, rebuild when max displacement > skin/2, checked between ALM
    outer iterations, `src/fq.rs` `relax_fixed`), then a trust-region SLP.
  - The two are the same idea.  LOFRT restarts only on box contact, while we rebuild on displacement.
- *Results (§6):*
  - non-identical ellipses, n = 5–120, with time limits of 1–48 h for "at least 10 local minima";
  - beats Kallrath–Rebennack on TC05–TC100 by 0.1–8.5 % (Table 1);
  - one 3D example with 30 polytopes.
  - No squares.
- *For us:* **confirmed, nothing new.**  The lifted separator and its angle chart are what `fq` has.  LOFRT is our
  neighbour list plus trust region.

Neither (a) nor (b) reports equal squares, jamming, or any escape move beyond multistart.  The only symmetry handling
is in Kallrath 2023 (one item fixed) and Cherri (§8).  **Nothing new for us.**  If we describe `fq` as a phi-function
method, cite (b) for the lifted separator and (a) for the active-branch polish.

**(c) Melashenko, Romanova, Litvinchev, Martínez Gómez, Yang, Sun, "A model-based heuristic for packing soft rotated
rectangles in an optimized convex container with prohibited zones"**
- *Record:* *Mathematics* 13(3) (2025) 493, doi:10.3390/math13030493, open access [verified abstract via DOAJ; the
  MDPI full text was blocked by a 403].
- *Method:* rectangles with fixed area and variable dimensions within limits, freely rotated, in a convex container
  whose size metric is minimised.  Phi-functions are combined with a model-based heuristic.  Tests use triangular,
  circular and polygonal containers.
- *For us:* the closest to "soft squares → rigid squares", but it still treats softness as part of the problem, not
  as a continuation.

I found no 2020–2026 phi-function paper on *equal squares in a square* (a targeted search turned up nothing).
`lit-2026-10-12.md` §7 reached the same conclusion for rotated equal rectangles in general.

## 13. 2D heuristics and 3D items (abstracts or titles only)

**2D heuristics**
- **Theodoracatos & Grimsley** (*Comput. Methods Appl. Mech. Eng.* 125 (1995) 53–70,
  doi:10.1016/0045-7825(95)00795-3) [verified record]: SA nesting of arbitrary polygons.  Superseded by everything in
  `lit-optimizers.md` §3.  Irrelevant.
- **Liu, Si, Hua, Jia**, GA + LP (*Appl. Sci.* 13 (2023) 12474, doi:10.3390/app132212474) [verified abstract]:
  NFP/bottom-left GA, then LP compaction, on ESICUP.  This is Gomes–Oliveira 2006 again (`lit-2026-10-12.md` §4).
  Irrelevant.
- **Kurpisz & Suter** (ESA 2023, LIPIcs 274:76; arXiv:2308.08523) [verified abstract]: approximation ratios for
  *translation-only* convex polygon packing.  Irrelevant.
- **Fang, Chen, Rao, Peng, Yan** (*Sci. Rep.* 14 (2024), doi:10.1038/s41598-024-81749-5) [verified abstract]: 2D
  irregular nesting sequences reused across tasks via a twin neural network, plus an improved bottom-left heuristic.
  Irrelevant.

**3D items**
- **Romanova, Bennell, Stoyan, Pankratov** (*EJOR* 268 (2018) 37–53, doi:10.1016/j.ejor.2018.01.025) [verified
  record]: concave polyhedra in a cuboid with quasi-phi-functions.  This is the 3D form of §12b.  Irrelevant.
- **Pankratov, Romanova, Litvinchev, "Packing Oblique 3D Objects"** (*Mathematics* 8 (2020) 1130) [verified
  abstract]: oblique cylinders and cones, phi-functions, cuboid of minimum volume.  Irrelevant.
- **Litvinchev, Pankratov, Romanova, "3D Irregular Packing in an Optimized Cuboid Container"** (*IFAC-PapersOnLine*
  52 (2019) 2014–2019) [verified record]: 3D.  Irrelevant.
- **Liu et al., HAPE3D** (*Front. Inf. Technol. Electron. Eng.* 16 (2015) 380–390) [verified record]: 3D constructive
  heuristic.  Irrelevant.
- **Chekanin & Chekanin** (LNME, 2023, doi:10.1007/978-3-031-48851-1_14) [verified record]: greedy placement of
  orthogonal polyhedra.  Irrelevant.
- **Ahn et al., "The minimum convex container of two convex polytopes under translations"** (*Comput. Geom.* 77 (2019)
  40–50) [verified record]: exact geometry for two objects.  Irrelevant.
- **Torquato & Jiao**: already covered (ASC).

---

## Ranked: worth trying

The list yields little.  Only item 1 is a new idea; items 2 and 3 are reuse notes for things already planned.

1. **Soft-square homotopy: aspect ratio at fixed area, then harden.**
   - *Sources:* §4, §5, §12c.
   - *Targets:* the **k-row jamming** failure.  A full row of k unit squares is rigid only because the widths sum to
     exactly k.  If each square may become w × 1/w (|log w| ≤ ε), a jammed row can shorten, and the slack moves into
     the cross direction.
   - *Schedule:* ε → 0 drives the state back to rigid squares, and either re-jams it or leaves a tilted defect.
   - *How it differs from the rounding schedule we tried:* it breaks the integer-width coincidence directly instead
     of smoothing corners.
   - *Cost:* S–M (1–2 days).  Add a per-square log-aspect variable to `fq` (corners scale by (w, 1/w) in the body
     frame, so the separator constraints stay smooth) and an ε schedule.
   - *Kill test:* at n = 54 and n = 71 from states jammed at s = k, the hardened outputs never go below k in about
     50 trials each.
   - *Prior:* low to moderate.  The papers give no optimisation evidence, and the hardening step may simply restore
     the grid.
*Items 2–3 updated 2026-10-12 after reading full text.*

2. **Exact window solves with a global MINLP (BARON/SCIP) on ≤ 4–6 free squares.**
   - *Sources:* §6, §7, §8.
   - *Targets:* **exact polishing** and local structure changes.
   - *What is new:* nothing beyond `lit-2026-10-12.md` rank 5, apart from a lower size calibration.
     - Cherri 2018 (§8): SCIP proves 3–4 pieces in under a minute.  It fails to close 5 or 6 pieces in 7 h (OOM, gaps
       of 35–80 %), even on exact-fit instances with warm starts.
     - Kallrath et al. 2023 (§6): "6" is six identical triangles with no timing reported.
     - Squares are simpler (4 vertices, θ ∈ [−π/4, π/4]) and a window has a fixed frame, so a few more may close.
       Plan around ≤ 4–6 free squares, not 10.
   - *What to borrow:* Cherri's symmetry constraints (§8).
   - *Cost:* M.  Use the Berthold code rather than any of these models.
3. **Cite, don't adopt.**  If we write up `fq`, describe it as a quasi-phi-function / separating-line NLP (§12b eq. 5;
   Kallrath 2009; Peralta 2018).  It is solved with a neighbour list (cf. LOFRT, §12b) and an active-branch SLP polish
   (cf. IPLOA, §12a), plus dual-ranked branch flips, which none of them have.
   - *Cost:* zero.
   - *Target:* documentation only; it does not address a failure mode.

**Not worth trying:**
- the sentinel penalty (weaker than SAT);
- the Rajagopalan phi-function correction (it concerns separated pairs only);
- all-pairs separation lines (O(n²));
- MIQCP for whole instances.
