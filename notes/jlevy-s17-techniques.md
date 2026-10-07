# jlevy/squares and s(17): methods, transferable techniques, bearing on our open items (2026-10-03)

Source studied read-only: `~/math/_untrusted-third-party/jlevy-squares` (github.com/jlevy/squares, main at
`79419cdf`, 2026-10-03; see `PROVENANCE.txt`), plus the unmerged PR #307 branch `claude/busy-goldberg-rouoli`
@ `234a07f4`, read through `gh api`.  Nothing of theirs was executed and no code was copied.  Paths below are
relative to their repo root unless they start with `s12/`.  Their code is MIT; their docs and data are CC BY 4.0
("Joshua Levy, the squares project").  The s(17) lower-bound certificates themselves are by Kleddamag and
Guzhou0806.  jlevy/squares reviews and replays those certificates; it did not produce them.

**Summary.**
(1) All s(17) bounds above 4.614 come from one architecture.  A per-angle-interval *strict core* is swept exactly
over the parent centres, and the charge is not just point mass: it includes capacity-one Boolean *rule atoms*
(k-of-m thresholds and pairwise-intersecting winning-subset families).  The rule atoms are what moved the bound
from about 4.614 (points only) to 4.66044.
(2) The rule atoms are, in effect, our all-meet anchor cliques, with the anchors taken to be small finite site
sets.  They can be checked inside the same rectangle sweep by Möbius inversion, and this is the biggest
transferable idea.
(3) R068 does **not** decide our ν_f ceiling test.  It bounds the *clique*-weighted value, and ν_f is larger
than that.
(4) In open PR #307, jlevy is already well into the s(11)-style exact route for s(17): the local half is proved,
and 63% of the global census is excluded.  This largely overtakes `s12/tasks/s17-core-isolation/`.

## 1. How the s(17) lower bounds were produced

### 1.1 Lineage (values as registered in `SYNOPSIS.md` and `packing/frontier/RESULTS.md`)

| Bound | Who | Charge | Where reviewed |
|---|---|---|---|
| 4.6141535 | Mira-acc | points only (2,881-direction net, "Condition 5") | `docs/project/reviews/review-2026-09-20-n17-r012-and-mira-4613-proof-review.md` |
| 4.61304 (461300/99999) | Guzhou0806 R012 | points only; 2,925-entry parent-angle catalogue | same review |
| 4.61979 (461300/99853) | Kleddamag v1.0.0 | R012 + 253 orbits of **2-of-3** rules (jlevy's T-025 threshold atom) | `review-2026-09-21-n17-kleddamag-461300-99853.md` F3 |
| 4.62002 | Guzhou0806 R052 | point/threshold, 15,721 intervals | `review-2026-09-25-n17-guzhou-r052.md` |
| 4.640020 | Kleddamag v1.1.0 | + weighted thresholds, + pairwise-intersecting winning-subset rules on 7 sites (Fano included), sites continuously optimised; 8,876 sites, 546 rule orbits, **coarser** 2,048 intervals | `review-2026-09-27-n17-kleddamag-4640020.md` §3, §8 |
| 4.66001 | Kleddamag | 20,856 sites, 889 rule orbits, 86 rule patterns on up to 10 sites, 2,168 intervals; points are only 29% of the budget M | `review-2026-09-27-n17-kleddamag-466001.md` §3 |
| 4.66018 | Guzhou0806 R067 | 4.66001's charge unchanged, finer catalogue (2,808 rows), rebuilt cores | `review-2026-09-28-n17-guzhou-r067-r068.md` |
| **4.66044** (116511/25000) | Guzhou0806 R068 | R067 + one site moved by 8e-5 + one weighted 4-site point orbit at (1.34, 1.34); 4,991 rows | same |

R070, R071 and Kleddamag's 4.6601 are **not** in jlevy's repo.  The frontier file `packing/frontier/n-017.md`
still reports R068 as best.  Our `site/notes/lower-bounds-notes.md` is the only place that tracks them.

### 1.2 Certificate type ("parent-core charge")

- **Scaling.** The container is fixed at L = 4613/1000 and the parent side is A = L/T < 1.
- **Catalogue.** The angle is parametrised by half-angle u = tan(θ/2).  The catalogue is a gapless chain of
  rational intervals [a, b] from 0 past tan(π/8); θ ∈ (π/4, π/2) is handled by reflecting the single parent
  across the diagonal, so no symmetry of the packing is assumed.  Each row fixes a *core*: a concentric closed
  square of side B at the midpoint orientation, with A − B·max_{u∈{a,b}}(cos δ + |sin δ|) > 0.  That condition
  puts the core strictly inside every parent in the interval; the endpoint check suffices, proved in the R067/R068
  review §3.  The centre envelope [r, L−r]², with r = A·min_{u∈{a,b}}(c+s)/2, contains every legal centre.
- **Charge.** F(Q) = Σ_R w_R·[R fires on core Q], with integer weights at scale 1e9 and complete D4 orbits
  carrying one weight each.  A rule R is a monotone Boolean predicate on the set of its sites captured by the core:
  - a point (1-of-1);
  - a weighted threshold Σ a_i·[site i captured] ≥ k with ⌊Σa_i/k⌋ = 1;
  - a list of winning masks that pairwise intersect.

  The budget is M = Σ w·cap·images.  Every rule in every certificate since v1.1.0 has **exact capacity 1**,
  checked by a DP over all 2^n capture patterns.
- **Theorem.** Γ := min F over all rows and all centres in the envelope.  Disjoint cores give 17Γ ≤ M, so
  17Γ > M excludes side T, and compactness gives the strict s(17) > T.  For R068, Γ = 1,000,026,844,
  M = 17,000,448,944 and 17Γ − M = 7,404, a relative 4.4e-7.

### 1.3 Checker

For each row, the sites are rotated into the core frame.  The set of centres at which a core captures site p is
then an axis-parallel B×B square.  Each rule is expanded by integer Möbius inversion into signed "all of U
captured" terms, and each term is again a rectangle (an intersection of squares).

The checker sweeps the exact arrangement of these rectangles: slabs plus a lazy range-add/range-min segment tree,
int64 accumulators under an a-priori bound Σ|w| < 2^50, and exact rational coordinates (BigInt, Boost `cpp_int`).
It takes the minimum over the open cells that meet the envelope.  Boundary and event-line centres are covered by
monotonicity, because captures are closed and rules are up-closed.  Cells are never enumerated: R068's
2.2·10^12 cells are implicit in the range-min queries.

Costs:
- **R068**: the C++ + Node paired run took 1,265 s on two workers (R067/R068 review §4).
- **4.66001**: Python took 163 s and Node about 3×120 s for 2,168 rows (466001 review F6).
- **jlevy's independent Python row sweep**: about 3.3 s per row.

### 1.4 Small tilts, zero margin, margins

There is **no special small-tilt handling and no zero-margin problem.**
- The chain starts at u = 0 like any other row, and the core shrink absorbs the angular spread.  B/A is about
  0.9998 at the base width (2.0e-4 in u), and 0.9999993 at R068's narrowest row (1/256 of the base width).
- Containment margins A − B·h are 1e-10 to 1e-12, exact rationals.
- Because each core lies strictly inside its *open* parent, the cores of touching parents are disjoint, and the
  target T sits strictly below s(17).
- The ledger is flat.  At 4.66001, R067 and R068, **every** row has the same minimum Γ, attained at the corner
  cell (r, r) (the parent tucked into the container corner).  At the five 4.66001 rows profiled, the next cell
  value is within 1e-7 of Γ.
- Every row's minimum sits at its corner.  At the rows they sampled, the R068 patch (the new 4-site point orbit)
  is not load-bearing.

### 1.5 How sites and weights were chosen

The producers are not published.  The sources mention:
- "numerical search", separate Codex tasks, "reoptimization", and "numerically proposed general intersecting
  rules" (Kleddamag `bounds/4.66001/METHOD.md`, retained under
  `packing/resources/web/n17-kleddamag-466001-2026-09-27/`);
- weights from an LP with "charge ≥ 1", rounded up at 1e9 and then re-audited exactly;
- exact "actual-parent probes" (a pose with 17F < M) used as counterexamples to drive refinement.

The sites are optimised continuously and end up on core boundaries.  As a result, the fixed-weight charge falls
in *steps* as the side grows: a 6.5% drop on row 0 within ΔS = 0.0005 (`review-2026-09-27-plan-4640020-lemma-check.md`
§4).  That is why each continuation (R067, R068) had to re-catalogue, move sites, or add mass.

Guzhou's next dictionary (FIT02, aimed at 4.6605) failed.  An exact counterexample sits at θ ≈ 12.6°, with a
deficit of 48,152 units per core (R067/R068 review F8).  jlevy's reading (same review F5) is that the
"fixed charge + finer catalogue + patches" recipe is spent at 4.66044, and that the next step needs a new
dictionary.

### 1.6 Compute

Verification takes minutes to about 20 minutes of wall time.  Production compute is not reported.
jlevy's own independent re-derivations took under 6 CPU-minutes (R067/R068 review) and 91 CPU-minutes (4.66001).

## 2. Techniques that could carry over to our pipeline

| # | Technique (theirs) | Where it plugs into ours | Expected gain | Confidence |
|---|---|---|---|---|
| A | **Capacity-one rule atoms** as extra columns: k-of-m with ⌊m/k⌋ = 1, weighted thresholds, pairwise-intersecting winning families on ≤ 10 sites. | The LP (`cover4_cg.py`, `quadrant_lp`) gets rule columns.  A rule's row coefficient on a pose is 1 iff some winning subset is inside the closed square.  In `zmx2`/`zm_mixed` a leaf credits w_R iff a winning subset ⊆ the surely-captured set, which is cheap because we already compute that set per leaf.  In `verify/` (the sweep) use Möbius terms = rectangle intersections, as they do.  Lean: an instance of our `clique_of_cores` / `card_filter_clique_le_one` (`lean/Sqpack/Basic.lean`), since the poses firing one rule pairwise share a site. | At n = 17 this is the difference between ~4.614 (points) and 4.66044.  At n = 12 our anchor-clique LP already measured about 0.10–0.12 of mass at 3.99 (`search/CLIQUE_CONTINUUM.md`).  Their rules are a cheaper, sweep-checkable member of the same family. | **High** for strict lower bounds (s(17), s(12)).  **Low** for zero-margin exact proofs: touching closed squares can share boundary sites, so capacity one needs an extra argument (sites off contact seams, or strict cores). |
| B | **Exact per-row centre sweep with a strict core** instead of a box tree over centres; **adaptive, non-uniform angle catalogue** (only failing intervals bisected: 825/994/349 pieces per base interval in R068); B chosen maximal per row rather than rounded down. | For any *strict* s(17) attempt, use the `verify/` design (it already sweeps per angle bin with a σ_k shrink), not `zmx2`'s pose box trees, which exist for zero margin.  Make the bin list adaptive instead of θ_k = 2arctan(k/N).  jlevy's s(12) run needed N = 96,000 (39,765 bins, about 800 s on 4 threads) to spend a 1.3e-4 scale; R068 needed 4,991 rows at n = 17. | Fewer rows by roughly 5–10× at equal loss (estimate).  It also removes germ boxes and the tiny-tilt problem for strict bounds. | Medium |
| C | **Tight-cell dump as LP rows** (jlevy's s(12) "Route B", `docs/project/research/research-2026-10-02-s12-beyond-rescaling.md`).  Our `verify` `TIGHT_DUMP` lists arrangement cells with charge ≤ threshold.  Membership is decided exactly in integers at the cell midpoint, the row is the captured *set*, and HiGHS weights are rounded up at 1e7.  16 rounds, 8.7k rows, about 2.3 CPU-h. | `s12/tasks/cert-slack/` lever 3.  Have `zmx2` and `zm_mixed` emit the captured set of each uncertified or near-tight leaf, not a rounded float pose.  This removes the "printed to 1e-6, useless near tangency" problem, because the row is combinatorial.  Also replaces part of the random hill-climb oracle. | Shown on *our* s(12) certificate: 1,736 points unchanged, reweighted, gave **s(12) ≥ 15680000/3949423 = 3.9702002** (ours was 3.9686155), with LP total still 0.024 below 12, so more is available. | High |
| D | **Certificate golf by rescaling + finer net** (squarepacker #309; their Route A).  Our s(12) file, scaled by 7902/7901, verifies at N = 24,000.  The binding square (bin 0, corner, 58 points) does not move under scaling.  The loss is net erosion, 1 − σ_k ≈ 2/N. | `cert-slack`: before re-optimising, try scale + finer net on any shipped certificate. | Small (+0.0005 on s(12)) but nearly free | High |
| E | **Python exact-arithmetic speedups** measured in PR #307 (`review-2026-10-02-n17-cost-reduction-performance.md` on that branch): gmpy2 `mpq` in place of `Fraction` 2.8–10.9×; homogeneous-integer predicates 17.6–23× over `Fraction`; Fraction dispatch/gcd is 55–65% of their Python provers' time; verifiers are 70% of total cost. | `zeromargin.py`, `zm_mixed.py`, `xcheck.py`, `dual_exact.py`, the qx2 exact paths: rebinding to gmpy2 is a few lines.  Rust `zmx2` is unaffected. | 3–10× on Python exact checkers | High (their measurements, similar workload) |
| F | Rust rationals: skip the final reduction after cross-cancelled multiplication (`research-2026-09-30-exact-arithmetic-verifier-performance.md`). | Anywhere we use `num-rational` in a hot loop. | 13–18% | Medium; minor |
| G | **Int64 accumulators with an a-priori mass bound** (Σ|w| < 2^50, integer weights at fixed scale) and exact geometry only for event coordinates. | `verify/` already uses i128.  Relevant if we write a strict-core rule sweep (B). | Constant factor | Medium |
| H | **D4 fold per parent** (reflect one parent across the diagonal; the charge is D4-invariant with one weight per orbit). | We already do this (zmx2 sym atoms). | None new | n/a |

Not useful: their LP tooling.  They also use HiGHS (single-threaded), and they publish no column generation for
n = 17.  Their native parent-core interval branch-and-bound (`devtools.verify_evand_angle_net_native`, about
0.35 s per row on our s(12) file) is a checker, not a producer.

## 3. Bearing on our open items

### 3.1 Ceiling test at t = 4.6604 (`s12/TODO.md` §s(17))

jlevy has a **clique-weighted ceiling lemma**: `review-2026-09-27-n17-kleddamag-4640020.md` §8, re-proved in
`review-2026-09-27-plan-4640020-lemma-check.md` §3.  If unit squares in [0,S]² carry y ≥ 0 with y(K) ≤ 1 on every
clique of the interior-overlap graph, and Σy ≥ 17, then no capacity-one strict-core certificate exists at S.
- **R068 already proves** that no such clique-weighted family exists at any S ≤ 4.66044.
- **That does not settle our test.**  Our ν_f constrains only point depth (Helly-type: μ{Q ∋ p} ≤ 1), so
  ν_f ≥ the clique value, and the lemma-check §3 shows by an explicit three-square example that the two differ.
  R068 says nothing about ν_f(4.6604).
- "Nobody has computed it" remains true.  Their H-248 (`packing/campaign/hypotheses/H-248-...md`) targets the
  clique value at 4.675 and 4.67; it is `blocked` and the instrument is not built (`packing/campaign/ledger.md`).
  Their only data point is a 2-minute triangle anneal at 4.65 that reached α* = 16.
- Our `search/CLIQUE_CEILING.md` (2026-08-29/30) computed the same kind of clique ceiling for n = 12 four weeks
  before their lemma.

Updated reading of our test:
- If ν_f(4.6604) ≥ 17, pure point/line/area measures cannot reach R068/R071, while capacity-one rules can.  The
  lever is then technique A, not more points or lines.  The prior strongly favours this: the points-only
  certificates stopped at 4.614, and points are 29% of M at 4.66001.
- A cheap, more decisive companion run: the clique-weighted value at 4.675 (their H-248) with our
  `clique_ceiling.py`.  A family there would cap the whole rule architecture below Bidwell.

### 3.2 s(17) exact via the s(11) method (`s12/tasks/s17-core-isolation/`, on hold)

jlevy is far along, in **open PR #307** (Sessions 167–168, unmerged; its body says "The n17 proof is
incomplete").

**Endpoint.** S* is proved to be the root of the degree-18 polynomial (H-255, H-265), the endpoint is exactly
feasible with all 68 walls and 136 pairs certified (H-256), and the rational ceiling is T-065.

**Local half (proved, as a capture-target theorem).**
- A packing of side ≤ S*, in the family's occupancy state, with its 45 non-slider coordinates within 1/5000 of
  the family's, lies on the family.
- A 58-row exact stress gives first-order stationarity in both corner branches (H-258).
- Exact duals exist for 90 signed directions; the worst ratio is 0.926 at r = 1/5000 (H-261).
- Slide coverage is H-268.
- Per-coordinate radii close down to 1/1216 (`review-2026-10-03-n17-local-radius.md`).  A uniform radius is
  limited to about 1/4630 by −ω₁₁.
- Their sliders: square 6 loose in a hole; 5, 11 and 13 slide.  This matches our 4 moving squares, which we
  label 0-based as {5, 12} alone plus {4, 10} together; the label mapping is not yet checked.

**Global half (63% done).**
- A 24-cell capacity-one cover at U = 1169/250 (4 corner, 12 wall, 8 interior cells, with a "depth-width wall
  lemma").  Our brief guessed 26–30 cells.
- 346,104 states in 43,593 D4 orbits; 15,953 orbits remain after sub-pattern exclusions W7, A, SW9 and N1.
- Their branch-and-bound stops at 7 cells.  The ownership-induction kernel (lifted from n = 11) needs adaptive
  angle rows on thin classes.
- The rest is priced at 4·10³–4·10⁴ CPU-h (`review-2026-10-02-n17-residue-process.md` on that branch).

**Capture.** Undecided: no contraction at 128 rows, and a 256-row run is in progress.

**Charge floors from R068** (H-262) were rejected.  At the cap, R068's charge collapses.  Any D4-symmetric linear
per-cell floor leaves ≥ 30,966 orbits, though an asymmetric floor reached 631 in a short search
(`review-2026-10-02-n17-charge-floor-pilot.md` on the branch).

**Implication.** Tasks 0–3 of our brief (exact packing, core system, first-order rigidity, radii) are done there
with exact checkers.  Restarting them would duplicate published work.  Options:
1. Leave it.
2. Write an independent re-check of their local theorem (a referee-style contribution).
3. Offer a fast pose-space excluder for their global residue.  Their branch-and-bound stalls at 8 cells; a
   zmx2-style float-enclosure box tree is the kind of tool that could help there.

Priority note: the s(11) method is Queuingtheorydotcom's, and jlevy is applying it to n = 17 openly.

### 3.3 k²−5 certificate cost; Lean

- Nothing in their n = 17 material prices k²−5.  X-048 ("New leads from the October Evand review") quotes our own
  k²−4 estimate back to us (5–10 CPU-h LP, 50–80 CPU-h check).
- Techniques C and E above apply directly to the k²−m loops: exact tight-cell rows, and gmpy2 in the Python exact
  paths.
- Technique A does not help zero-margin k²−m (see the A row).
- Lean: the rule-atom counting lemma is a few lines on top of our clique lemma.  jlevy has no Lean for n = 17.
  Their Lean work concerns our s(13) and the k²−4 reduction (§4).

## 4. Where their repo reviews, cites or relies on our work, and open asks

Sources: their reviews `docs/project/reviews/review-2026-09-27-evand-s32-s12.md`, `-09-28-evand-s21-s45-mixed-covers.md`,
`-10-01-evand-mathematical-transfer.md`, `-10-01-evand-source-coverage.md`, `-10-02-evand-s32-no-fold-run.md`,
`-10-02-evand-s60-geometric-premises.md`, `-10-02-valid7-independent-checker.md` and `-10-03-evand-k2m4-claim-chain.md`;
issues #238, #256, #279 and #316; PRs #298, #311, #320 and #322.

**Status.** None of their reviews found a mathematical defect in any of our claims.

| Rating | Results |
|---|---|
| V3/C3 | T-050, T-051 (s32), T-052 (s21), T-053 (s45), T-062 (s60), T-063 (s61), T-064 (k²−3, with `bentz_of_valid7` built on standard axioms) |
| Confirmed | T-006 (case-free s13, kernel-checked by them as `s13_eq_4`) |
| V0/C1, "reviewed" | T-081, k²−4 (PR #322 built `bentz4_of_validTilt9`; #316 is replaying `ValidTilt9`, priced at 226 CPU-h) |

**Our s(12) bound has been superseded using our own certificate.**  The new bounds are T-078 (squarepacker,
rescaling, 3.969118, #309) and T-079 (Levy, reweighting, 3.9702002); see §2 C and D.  T-049 (our 15680/3951) is
now marked superseded.

**Requests and open items addressed to us:**
- **#238** (closed, but we did not answer the questions):
  - `s12/verify` exits 0 on NOT VERIFIED and on partial runs.  They want a nonzero exit; squarepacker's replay
    instructions rely on it.
  - `cargo build | tail -1` masks failed builds in `search/zmx2_run.sh` and `certificates/s60/verify.sh`.
  - Smaller points: a latent `i as u32` cast in `zmx2.rs`; the verify header says "interior" but the code uses
    closed squares; they would like to cite `tasks/s21-finish/xcheck.md` if we publish it.
  - We offered a native zero-margin checker spec there and have not delivered it.
- **#256** (closed; question unanswered):
  - Confirm their reading of the stale "Provenance of zm_mixed_d4/" paragraph in the s60 README.
  - Pin or print the `zmx2` source digest, because `verify.sh` builds `1cd4dcbd` while the records came from
    `6b7f0f79`.
- **#279:** we offered to replay wand125's s59/s77 covers on current `zm_mixed` / `zmx2 --full --sym-atoms`, and
  wand125 said yes.  This is owed unless already done.  wand125 also reports s(76) (k²−5, k = 9) stalling at about
  6k uncertified boxes.
- **Review notes not yet answered:**
  - k²−4 R-4: nothing machine-links `qx2_zm`'s (c, u) boxes to the Lean `sq c (2 arctan u) 1`.
  - k²−4 §5: our three adversarial review reports are private, so they cannot be cited.  This is an implicit
    request to publish them.
  - s21/s45 F3: `verify2` has no `overflow-checks`.
  - s21/s45 F6, s60 F3: no Lean data for s45 or s60.
  - s60 F4: they offer their cover-mutation controls to us.
  - Source coverage: our site's sources §5 still shows R067 rather than R068 for s(17), and lacks k²−3 and
    Queuingtheorydotcom.
- **Closed:** the side-3.99 supports (published 2 Oct and acknowledged), and the anchor-clique Lemma 2 convexity
  repair (fixed 2 Oct, `AnchorLemma2.lean`).

**n = 17.**
- Neither `packing/frontier/n-017.md` nor PR #307 relies on our work.
- X-048 R7 proposes trying our endpoint mixed-measure idea at T17.  The "New leads from the October Evand
  review" section lists anchored cliques and exact limiting-pose constraints as transferable.
- jlevy has not registered R070/R071 or Kleddamag's 4.6601.

**Verifier.** jlevy's clean-room `sqverify_fast` (PR #311, Rust, directed-rounding interval branch-and-bound) does
not yet read our continuous-angle closed covers.  Its Milestone C plan
(`docs/project/specs/active/plan-2026-10-03-measure-verifier-milestone-c.md`) targets our points (s13, s32) and
then lines with Lemmas Z/H/DP (s21, s45, s60, n59, n77), timed black-box against our checkers.  It would be a
genuinely independent second route for C4 ratings.
