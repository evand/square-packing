# packer — packing search (started 2026-10-04)


**Naming (10-05):** "candidate" / `cand110` / run label `cand` throughout = **Couzo's 2026-09-27 s(110) record**
(`../exact/results/couzo110.*`; start file `candidates/s110_10.996783403.txt`), found by us independently a week later; "record"
/ `rec110` = Ellsworth's earlier 10.99679327….  Not ours; we add only its exact optimum (2.8e-12 below his printed bound).

Goal: better packings at WISHLIST §P targets, first **s(90) < 10** (P1).  jlevy register 10-04: s(90) ∈ [9.6, 10], grid best known
(wand125 claims ≥ 9.725, unreplayed).  Screening only (f64); certification would be a separate exact step.

## Tools
* `packer` (Rust, `cargo build --release`): fixed-side overlap energy (SAT penetration², walls), L-BFGS, Verlet list.
  `relax [--squeeze]` (squeeze = shrink side while feasible → the basin's minimal side); `hop` = basin hopping with descending
  target (moves: disk shake, relocate, re-angle, cluster rotate/translate, swap, strip slide, best-hole relocate).
  Text format: `n s` then `x y deg` per line.  `site2txt.py` / `txt2site.py` convert from/to the site JSON.
* `gen.py`: tilted set + axis fill; `mis_fill` = exact max independent set of corner-lattice axis candidates (HiGHS MILP).
  Reproduces the full count of records 54, 87, 88, 132 from their tilted squares alone; 71 and 110 need slid axis squares (−3, −1).
* `channel.py`: Cantrell-type gadget = left-wall + bottom-wall channels of tilted rows (w per row, period 1/cos θ).
* `surgery.py`: delete a horizontal + vertical unit band from a record, close the gaps (k → k−1).

## Calibration (from random starts)
* n = 11: 4/4 seeds reach Trump's 3.877084 in ≤ 8 s.  n = 17: 3/4 reach 4.67554 (record 4.67553) in ≤ 6 s.
  (jlevy's sqsearch: grid 5.0 at n = 17, 5/5 seeds.)
* n = 54, 71: fail from random, at s = k (grid basin) and at s = k − 0.02 (stalls at E ~ 1e-4: jammed axis rows).

## n = 90 so far (all negative)
* 110 → 90 surgery: 821 distinct seeds; pure relax+squeeze best 10.00015 (a single bottom 62° channel); hopping near-misses at
  9.9999 converges to **exactly 10** every time.
* **Positive control fails**: 132 → 110 surgery (1020 seeds) never relaxes below 11, though Cantrell's s(110) < 11 exists.
  Hopping the 12 best of them (10 min each at s = 10.999, all 8 moves) also fails: E stalls at 1.7–6e-7 (jammed at 11).
  So the 90 negatives say nothing yet: the instrument must pass 132 → 110 before a 90 result means anything.

## Diagnosis
Side k is a giant attractor: any state containing a jammed line of k axis squares has basin minimum exactly k, and at target
k − ε its energy (~ε²) is indistinguishable from a genuine near-miss.  Fixed-side energy therefore cannot rank candidates near k.
Ideas: target well below k (k − 0.01…0.05) and track E_min(s); a penalty or move that breaks full axis lines; search on gadget
topology (channel widths, row counts, junction squares) with the axis region solved exactly (MIS / LP), not by annealing.

## Open question (Evan, 10-04)
How many locally optimal (rigid) packings exist at n = 90, and how are their sides distributed?  Measurable: quench many starts,
squeeze, dedupe by side + contact graph.  Expect exponentially many basins jammed at exactly 10 and rare sub-10 ones; that would make
"hard for annealing" precise.  Cf. jlevy H-012 (record vs modal basin attraction ratio, n = 10, 11; not run).

## Literature (10-04; full reports `lit-optimizers.md`, `lit-statmech.md`, `lit-provenance.md`)
* Ellsworth on s(55): "without special modifications, the simulated annealing algorithm almost always gets stuck just above the
  trivial size … due to stacked rows and/or columns": our failure mode, seen by the record holders too; the modifications are unpublished.
* Records in 50–135 come mostly from constructions extended across n and from annealing seeded with known records; from-random annealing
  (Schadt, GPU) found 50, 51, 53, 55, 103, 105 at about GPU-hours per hit.  No public annealer.  Ellsworth's public repo has an exhaustive
  45°-gadget enumerator and a high-precision refiner (read-only reference; reimplement, don't copy).
* s(88): "Improvement by Thomas Schadt pending"; a Cantrell/Ellsworth Jan 2025 technique "will need to be applied to about 13 additional
  packings previously thought to be finished".  So best-known values at n < 100 are not settled.
* No squares-specific method published since Gensane–Ryckelynck 2005.  Transferable: LP shrink/rigidity step (Donev–Torquato–Stillinger–
  Connelly 2004); GLS overlap weights + global relocation (Sparrow, Gardeyn et al. 2025, Rust); population basin hopping (Addis–Locatelli–
  Schoen 2008: 32 circle records, n ≤ 130); contact-graph hashing to reject seen structures; replica exchange in side length; census via
  repeat-hit statistics (Good–Turing/Chao), cf. Frenkel/Martiniani basin-volume work.

## Plan (10-04)
Pass/fail test: reproduce s(110) < 11 without s(110)-specific hints (from family structure or s(132)-type seeds, not s(110) itself).
1. LP shrink/rigidity step: replaces squeeze; certifies jams; answers "rigid?" and feeds the census.
2. Row-jam detector + row-breaking moves (pull a square from a full line, reinsert tilted into the slack) + GLS overlap weights.
3. Gadget enumeration: Cantrell-type channel families, axis region solved exactly (MIS/LP), only θ and offsets continuous.
4. Approximate basin census at a few n (cheap once 1 exists) → "why good packings are hard to find" writeup.

## Step 1–2 results (10-04)
* `rigid.py`: first-order LP shrink/jam test + per-coordinate rigidity (FD-checked gradients).  Records are first-order jammed and
  nearly rigid (s(110): 104 rigid / 4 partly pinned / 2 free; s(89): 66/15/8; s(11): 11/0/0).  Our n = 90 jams at exactly 10 have
  0 rigid squares, 32–40 partly pinned (the jam lines), 50–58 free, and 77–173 contacts (records: 380–560).
* `slp.py` (LP step + packer relax): 10.00015 → 10.000093 on one near-miss, then stalls while the exact LP still finds a first-order shrink
  (second-order blocked).  Marginal as a squeeze; useful as a flag.
* Jam-line detector + row-break move + line penalty in `hop` (`--linepen`): no gain at n = 71 from random (also stuck with 0 lines,
  E ≈ 2e-4); compression from s = 10 reaches 9.000 in seconds and then sits on the grid attractor (12/12 runs, 5 min).
  Breaking a full grid needs many lines changed at once; single-square row breaks can't.  Overlap-weighting (GLS) untried.

## Step 3 + event chains (10-04)
* `gen.mis_fill2` (corner lattices + staircase candidates touching tilted squares, exact MILP) recovers the full count of records
  54, 71, 87, 88, 110, 132 from their tilted squares alone (~0.1 s).
* `gadget_mc.py` (state = tilted squares, score = |T| + exact fill): from s(132)'s tilted set in a 10.999 box, f = 109 of 110 at
  once, then 8/8 runs × 15 min never move: every Gaussian move or tilted-only relax destroys 3–10 squares of fill (needle landscape).
* `packer hhop` = loosen by γ, event chains (slide until contact, pass the rest to the square hit; rotations stop at contact),
  squeeze; Metropolis on side.  Unlock test n = 71 from states jammed at 9.0: 6/6 stay at 9.000000000x (15 min).  Pass/fail 132→110:
  best 11.0000000000 (8 runs).  Local moves, even long collective slides, don't change the topology.
* Topological reading: s < k means every one of the k rows and k columns must lose an axis square, so the tilted part must be a
  transversal of all rows and all columns (diagonal band; Cantrell's L of two channels + corner block).  Jammed-at-k states have
  full lines; getting from there to a transversal is a global change.  So: generate transversal topologies constructively, and use
  only contact-preserving moves inside a topology (rigid channel/row moves, event chains within the gadget, exact axis fill).

## Constructive layouts (10-04, `layout.py`)
* Genome = rows per channel (w, θ, set-back k); construction = each row on its wall, slid down/left until contact (transpose for the
  bottom channel), + exact axis fill.  Records' channels are tapered lattice blocks (s(132) left rows 2,3,4,4,3,2; lower rows set
  back 1–2 squares from the wall), hence the set-back gene.
* Integer count is the wrong objective (piecewise constant, jumps when a channel edge crosses an integer width): MC on count reached
  128/132 at s(132) and 106/110 at 10.999.  Count vs s for a record-like genome: 132 only at s ≈ 12.4 (construction is imprecise).
* Hybrid objective s_hyb = build at first side with count ≥ n, then `packer relax --squeeze`: record-like genome 12.45 → 12.068–12.094;
  hhop 5 min → 12.057 (record 11.991).  Genome MC on s_hyb running (~6 s/eval).
* Genome MC on s_hyb drifts to the grid (14/14 runs end at exactly k: channels shrink until the fill is a full grid).  With layouts
  containing a full line of k axis squares rejected (`layout.full_lines`): 132: 12.068 → 12.022 (record 11.991); 110: 11.030 →
  11.016 (need < 11; record 10.9968).  ~150 evals per 20-min run.
* hhop polish of those: without constraint 10/12 runs fall to k + 1e-7; with `--forbid-lines k`: 11.0149 → 11.0148, others stuck.
  The basins found are tight local optima ≈ 0.015–0.03 above the records.
* Record channels are not clean "wall-anchored lattice segments slid into contact": s(110)'s left rows have set-backs 0, 0, 0,
  ≈0.85, ≈1.0, ≈1.74 (not multiples of cos θ) and mixed angles; so the construction family probably doesn't contain the record.

## Throughput + squeeze dynamics (10-04)
* s_hyb cost was 72% the bisection squeeze (4.9 of 6.8 s).  New `relax --squeeze-pen`: L-BFGS on s + μ·E(config; s) over (config, s),
  μ = μ0 … 1e9 ×10 per stage, then minimal growth to strict feasibility.  Eval now ≈ 1.5 s.
* μ0 acts like a temperature.  Low μ0 (10) explores (construction bench: 11.037 vs 11.052 at μ0 = 1e5) but drifts off exact optima
  (s(110) record → 10.99746, s(71) → 8.94459); high μ0 (≥ 1e5) keeps records to ~1e-7.  The squeeze path picks the basin.
* `hhop --soft p`: soften-and-reharden move (penalty squeeze from μ0 ∈ [10, 1e4]); `--forbid-lines k` rejects states with a full line.
* Layout genome: continuous row set-backs (records have 0.85, 1.0, 1.74).
* 110 batch with fast eval + continuous set-backs (14 × 20 min, ~700 evals each): best 11.0123 (found at 151 s, then plateau).
* **Records are not rigid lattice rows.**  Chaining squares exactly one unit apart along u splits s(110)'s 46 tilted squares into
  27 chains (16 singletons); neighbours are staggered/sheared, the left channel mixes 24.9–32.1°, and many chains rest on axis squares
  or walls (invisible to the construction).  A genome read off the record by hand scores 11.45.  So the row-genome family has a floor
  above the records; it finds relatives, not the records.  The records look like continuous optimisation of a constructed start.

## ftmc: free tilted squares + exact fill + penalty squeeze (10-04, `ftmc.py`)
* State = squeezed packing; step = perturb tilted squares (jiggle, rotate cluster, delete, add neighbour, shift angle group,
  re-angle, move) → relax tilted alone → exact axis fill (axis–tilted overlap ≤ 0.03 allowed, squeeze resolves; axis–axis exact)
  → `--squeeze-pen` (μ0 = 1e3) → reject full lines → Metropolis on side.  ~0.7 it/s.
* Exact (1e-6) fill made every non-delete move infeasible (tiny tilted motions kill whole lattice rows): hence the 0.03 tolerance.
* 110 test, 14 × 20 min from layout results: best 11.0123 → 11.0091 (b_5, T = 0.004, 139/916 accepted).  Several runs accept ~0
  (most moves create full lines).  Our branch differs from the record's (short left channel, 5-tall bottom channel vs the record's
  long left channel, 4-tall bottom): branch changes need many squares to move between channels → diversity of starts matters.

## Basin entry (10-04, `basin.py`: loosen 2 %, Gaussian σ on all coords (angles 20σ deg), `--squeeze-pen` μ0 = 1e3; 20 trials/σ)
* s(110) record: never returns to 10.996793 (±1e-5), even unperturbed (→ 10.99813); but σ ≤ 0.01 → 60/60 below 11 (median ≈ 10.9983),
  σ = 0.03 → median 10.9982 (max 11.001), σ = 0.1 → median 11.0006.  Exact optimum sharp; the sub-11 funnel is broad (total
  displacement ~0.45 still lands < 11).
* Our 11.0091: different broad, rough funnel (11.003–11.008; perturbation sometimes improves).  ftmc 40-min continuation plateaus ≈ 11.008.
* Reading: the target is not a continuous needle; it's a discrete choice of funnel (channel counts, junction, branch), separated
  by many-square moves.  Next: funnel-level search/census over structure classes, with descent only inside a funnel.
* ftmc 40-min continuation (14 runs): best 11.008067; plateau.
* **Our local descent is broken near jams.**  Landings from the record (σ = 1e-4) differ (40/40 distinct sides, spread 6e-4), polish
  stably (stiff penalty squeeze, then bisection: < 1e-6 change), but fail the LP jam test (ds* ≈ −1e-4, 0 rigid, 365–408 contacts vs
  the record's 563).  slp.py stalls too (10.997922 → 10.997921, ds* still −9e-5).  Likely cause: nonsmooth SAT energy at face-to-face
  contacts (tied axes).  So the "distinct minima" were stall points, search comparisons were noisy at ~1e-3 (the record's margin is
  3.2e-3), and basin-return failures say more about the optimiser than about basins.
* Next: a correct local minimiser (min s s.t. separation; active set / SQP, right axis per contact, face–face handling), passing
  (1) record returns to 10.996793 from σ = 1e-4…1e-2, (2) every output passes the LP jam test.

## Correct local minimiser (10-04, `slp2.py`)
* Trust-region SLP on min s s.t. linearised SAT separation, no penalty relax; step error repaired by scaling config + side by (1+p).
  Two bugs found: (a) kinks of |cos α|,|sin α| (nearly aligned squares, α ~ 1e-6) were split only for |·| < 1e-9 → first-order
  step errors (≈ 1.9 R); now split when within 3R; (b) HiGHS's 1e-7 feasibility tolerance swamped small steps → LP solved in units
  of R with tolerance 1e-10.  Now converges to first-order jammed states (exact-contact ds* = 0), ~20 s per descent (most of it a
  long tail; within 3e-5 after 3 s).
* Stalled penalty landing u3: 10.99792 → 10.9968134 (jammed; distinct minimum 2e-5 above the record).
* Census around the record with slp2 (`census2.py`, 8 trials/σ): σ = 1e-4 → record ×4, 10.9968029, 10.9968652 ×3 (all jammed);
  σ = 1e-3 → record ×1 + 10.99684, .99685, .99689, .99812, .99834, 11.0028, 11.0071; σ = 1e-2 → 11.0077–11.081 (none < 11);
  σ = 3e-2 → 11.007–11.44.  ≥ 7 distinct jammed minima within 2e-3 of the record from 16 trials; repeats now occur.
* Two descent maps: strict local (slp2) sees a small sub-11 basin (σ = 1e-2 never returns < 11); soft penalty squeeze (loosen 2 %,
  μ0 = 1e3) lands < 11 at σ ≤ 0.03 → it acts as a funnel-level compaction.  Pipeline: soft squeeze, then slp2 to finish.
* Re-finish of all ftmc outputs with slp2 (`refin.py`): everything improves.  Best without full lines: 11.00765 (b_5 family,
  was 11.0082).  a_2 states 11.016 → 11.00012 and b_1 11.0098 → 11.0003, but with 3–8 full lines of 11 (heading to exactly 11).
  Some b_5 states jam at 11.0047–11.0063 but acquire a full line during the finish.  Stalled scores mis-ranked the states.
* HiGHS prints bound-shift warnings on some LPs (badly scaled rows); results unaffected so far, worth a look.

## Census with soft squeeze + slp2 finish (10-04, `census3.py`; 30 trials/σ; summary rebuilt from saved trials)
Escape curve P(s ≥ 11 | σ) from the s(110) record: σ = 1e-3, 3e-3, 1e-2: 0/83; σ = 0.03: 2/28 (1 full-line); σ = 0.1: 12/30
(7 full-line), median 10.99841.  So σ50 ≈ 0.1–0.15 under this descent ("hot enough" for the record's funnel).
From our 11.00765 state: never below 11 at any σ; σ ≤ 0.03 mostly back to 11.00756 (robust own funnel); σ = 0.1: 22/30 jam on full lines.
**Candidate improvement of s(110):** 32/141 trials end below the posted 10.99679327401957: two repeated minima, the record
(10.9967933, ~13 hits) and **10.9967834** (~19 hits).  Best: `candidates/s110_10.996783403.txt`, s = 10.996783403149346 (Δ = −9.87e-6);
40-digit check of the f64 coordinates: min pair separation +1.5e-15, min wall clearance +1.9e-16; no full lines.  Same structure
as the record (rms displacement 0.003, 4 squares > 0.01, max angle change 1.1°): a better-finished neighbour, consistent with the
record's "not yet analytically optimized" note.  To do before any claim: exact/rational certificate, check the live record and
jlevy register, ask Evan before any external post.

## Candidate vs record: structure (10-04; `contacts.py`, `stress.py`, `homotopy.py`)
* Difference is at the junction: record's near-axis squares #62, #63 (1.12°) and #65 (0.16°) become exactly axis-aligned in the
  candidate, snapping onto the top-right axis lattice (~0.01 shift); #53/#54 shift 0.014, #41/#42 0.005; #67 26.1° → 25.2°.
* Contact role inversion: record has corner 1 of #62 on side 3 of #76 (also #62's corners on #42, #48); candidate has corner 3 of
  #76 on side 2 of #62.  #65 loses its corner contacts on #54 and #63.  Going between them passes a corner–corner configuration:
  a different contact graph / stratum, not a point on the record's plateau.
* Angle homotopy (#62, #63, #65 rotated record → 0, rest re-minimised by slp2 at each step): s rises 10.996793 → peak 10.997468
  at ~0.5° (the non-jammed step = the flip), then 10.996846 at 0°; released 10.996845 (≠ candidate 10.996783: other moves also
  needed).  Barrier ≤ +6.7e-4 on this path; both in the same sub-11 funnel.  Straight-line interpolation is a bad path (overlap
  6.6e-3).
* Distinct contact graphs among 128 sub-11 census minima: raw near-contact graphs all distinct (floppy axis squares carry tiny
  rotations and incidental contacts); channel-square contact graphs 114 distinct (still dominated by floppy neighbours).  Needs the
  force-bearing (max-support dual) network: part of the exact-contact solver task.  Junction-state classes (squares tilted
  0.05–2°): 10 classes; none (candidate class) 42 configs 10.9967834–10.99886; {62, 63, 65} (record class) 75 configs
  10.9967933–10.99864; 8 rarer classes 10.99688–10.99886.
* Census push (`runs/cen4`, 480 trials from record + candidate, σ 0.01–0.1): below 11: 60/60, 57–58/60, 50–53/60, 34–36/60.
  330 jammed sub-11 minima w/o lines, 51 distinct sides; top basins 10.9971888 (94), candidate 10.9967834 (75), 10.9972187 (35),
  10.9968134 (18), 10.9968259 (16), record (15).  Nothing below the candidate.  Chao1 33 (σ ≤ 0.03) … ~700 pooled (unstable).

## Exact-contact solver (10-04, `../exact/`, see `../exact/README.md`)
* Input packing → load-bearing contacts (max-support jam dual) → KKT Newton in mpmath (80+ digits) → λ, second order,
  flat-mode probe, corner–corner MILP → rational certificate checked exactly (`../exact/verify_cert.py`, stdlib only).
* **Candidate:** exact KKT point S = 10.99678339663159163950…; certificate `../exact/results/cand110.cert` proves
  s(110) ≤ 10.996783396631591639612746267957 exactly (9.877e-6 below the posted record).  λ > 0, PSD modulo 17 exact
  flat slides, jammed in all corner–corner branches.  Record's own exact point: 10.99679327395374924222… (6.6e-11 below
  the posted value).  Not posted anywhere.
* **Census false jams:** slp2's fixed separating axis jams at corner–corner touches.  5 of 7 tested cen3 minima are not
  local minima (MILP descent); after a kick 2 go to the candidate, 2 to the record, 1 to 10.996866449; `../exact/pipeline.py`.

**Correction (10-04, late): not new.**  jlevy's register (https://jlevy.github.io/squares/cases/110.html) already lists Francisco Couzo, 2026-09-27: s(110) ≤ 10.996783396634359 (verified there by exact rational packing and interval arithmetic).  Our candidate is an independent rediscovery of that packing; our exact KKT point 10.99678339663159163950… and certificate S' = 10.996783396631591639612746267957 are 2.8e-12 below Couzo's printed bound (analytic optimisation of the same basin).  Ellsworth's page (and our mirror) still show 10.99679327401957.  Lesson: check the register before calling anything a candidate.
* **Register check (10-04):** jlevy's register is ahead of Ellsworth's page and our mirror at the frontier: Couzo 2026-09-27
  improved s(110) (10.996783396634359), s(130) (11.911187706548756, was 11.91119052015898), s(132) (11.991327887694469, was
  11.99137344423647; Casson 09-23: 11.99134529315214).  s(88) unchanged (Ellsworth 9.88815305375857, degree-20 algebraic).
  Always compare against the register.  Probe (perturb → soft → slp2, 20 trials × σ 0.003, 0.01): s(88) returns to 9.888153 every
  time (LP says not jammed: corner–corner false jam).
* Probe results (10-05): s(130), 39 trials (σ 0.003, 0.01): best 11.911191084, nothing below the mirror's 11.91119052 or Couzo's
  11.911187707; most land at 11.9124–11.9127.  s(132) not run.  **Bug:** one s(130) trial ran 9 h 19 min at 100 % CPU without
  finishing (slp2 / repair loop with no wall-clock cap), hanging the pool; killed.  Add a per-trial timeout (pool imap + timeout,
  or a time budget inside slp2 and a cap on repair growth iterations).

## slp2 rebuilt: incidence model, false jams gone (10-05; `inc.py`, `jobpool.py`)
* **Cause of the false jams was wider than corner–corner.**  rigid.pair_rows linearised a pair on one separating axis: "both
  near corners of t outside o's side line".  At an offset side–side contact one of those corners lies past the end of o's side,
  so feasible relative rotations were forbidden.  The fix is the exact solver's incidence model (each segment end: corner on the
  other's side), with true corner–corner touches as disjunctions (~60–80 at n = 110, vs ~210 tied axes in the old model,
  mostly grid diagonals).  For a corner–corner alternative "line L separates", every corner of the other square within range
  must be outside L, including ones past the side end (a 9e-5-rad tilted neighbour's far corner dips below L).
* Disjunctions in the loop: LP relax (drop corner–corner pairs) → round each to its least-violated line → restricted LP.  The
  exhaustive MILP over all axes was useless (210 binaries: 60 s limit hit without finding the LP's own point); over incidence
  lines (~75 binaries) it takes 0.1–2 s and is the final jam test.
* Validation, `../exact/inputs/cen3_*` (the 7 census minima the exact solver checked): the false jams now descend to the record /
  candidate; 0.003_10 stays (genuine); 0.1_28, unresolved before, converges to a strict local min (exact: PD, λ ≥ 2.3e-4,
  10.99756726622101311…); 0.01_24 reached 10.99686519433409872… (certified, λ ≥ 3e-5; below the pipeline's earlier
  10.996866449 whose min λ was 2.7e-7) and with the final model goes to the record.  5/5 slp2 outputs certify directly in
  exactsolve (no descend.py rounds).
* Speed: 20–50 iterations per descent (was 300–650: the over-restrictive model zigzagged), ~3–5 s at n = 110.  Second-order
  correction (side fixed, box ≈ 10x overlap; an unboxed correction LP returns an arbitrary vertex and moves everything by R)
  lets long descents progress (11.1884 → 11.1554 in 480 iterations where the old one stalled).
* Mini census `runs/cen6` (rec110, σ 0.01/0.03/0.1 × 28): 84/84 jammed, 0 timeouts; distinct sides 4 / 8 / 16 (old cen5 run with
  the axis model had many more), best = candidate.  Older census counts of "distinct minima" were inflated by false jams.
* **The 9 h hang** (s(130) σ 0.003 #11, reproduced with the old code): stuck inside one HiGHS `linprog` call.  Now: HiGHS
  `time_limit` on every LP/MILP, slp2 `budget`, and `jobpool.run_jobs` kills overdue trials (process per job).
* Unpinned OpenBLAS spins ~32 threads on import (4 s CPU per `import scipy.optimize`): rigid.py pins BLAS threads to 1.
* Open: exactsolve refuses two long-descent outputs (r_0.03_12 at 10.9984139, r_0.1_0 at 11.1554): "no equilibrium with smooth
  contacts", its MILP slope −1.2e-7 (noise).  Likely f64 convergence of near-parallel incidences.

## Clean s(110) census (10-05 evening; `runs/cen7`, commit 26436c6; `census3.py` + `census_exact.py`)
* Protocol: refs `seeds/rec110.txt` and `candidates/s110_10.996783403.txt`; σ ∈ {0.01, 0.03, 0.05, 0.1} × 125 trials each
  (1000); census3 (perturb, loosen 2 %, soft squeeze μ0 1e3, slp2 R 1e-3, budget 600 s); every sub-11 output through exactsolve,
  classes keyed by S_exact to 1e-20.  19.5 min wall, 3.6 CPU-h on 15 processes (≈ 6.5 CPU-s per trial + ≈ 4.5 per exact solve).
  0 failed / killed trials.  Pilot `runs/pil7` (64 trials) agrees.
* **19 distinct certified sub-11 minima** (812 certified outputs; Chao1 23; Good–Turing unseen mass f1/N = 4/812 ≈ 0.005 pooled).
  The top four hold 91 %: candidate .99678339663 (219), record .99679327395 (119), .99718876014 (236), .99719765403 (166);
  then .99841389 (19), .99841239 (12), .99756727 (7), …; nothing below the candidate.  rec and cand starts give the same
  classes in similar proportions.
* Per σ (both refs): distinct classes 4 / 6 / 8–9 / 15–16; unseen mass 0 / ≈ 0.01 / 0.02–0.03 / ≈ 0.09.  P(end ≥ 11) 0 / 2.4 %
  / 8.8 % / 42 %, so σ50 ≈ 0.11–0.12 (angles perturbed by 20σ deg).  **The old cen4 count (51 distinct sides, Chao1 to ~700)
  was false-jam inflation; don't cite it.**
* exactsolve gaps: 53 sub-11 outputs unresolved, 48 of them "solve failed" on the two .99841 classes (the known near-parallel
  incidence refusal: same f64 sides as the certified members, so almost surely the same minima); 1 TypeError ('NoneType' not
  subscriptable, a `runs/cen7/exact/*.json` with status error); 3 certified with a valid certificate but cls unresolved
  (rigid at .997568, two "local minimum" at .998412: certificate check failed?).  6 certified with second order inconclusive.
* ≥ 11 (134 trials): 71 at exactly 11 with full lines (tilted 40–49); ~60 jammed minima without full lines at 11.001–11.08
  (repeats: 11.0042943 ×7, 11.0009674/85, 11.0197992 ×2 …); 45° structures numerically at 9 + 3/√2 = 11.1213203 (×8) and
  7 + 3√2 = 11.2426407 (×4, 80–88 tilted).  f64 keys only.
* **Calibration reading.**  Inside the record's funnel the landscape is small and well sampled: ~20 sub-11 minima, the best is
  the 2nd most common basin, and ~100 trials at σ ≤ 0.05 find it with near certainty.  The needle is *entering* the funnel,
  not descending in it.  Next measurement: P(enter a sub-11 funnel) from defined outside starts (count-pair seeds, then PT).
* **Roles (`census_roles.py`).**  Every certified sub-11 minimum keeps the start's assignment, same indices: left 22, bottom 21,
  axis 67.  Exceptions are junction deformations, all ≥ 10.9978: #106–108 (left channel, at the junction) rotate 32° → ~46°
  (a third angle, not a transfer to the 62–66° channel; 3 classes + 1 with #106 alone), or 1–2 axis squares near the junction
  tilt 3–5°.  No channel ↔ channel or grid ↔ channel exchange below 11.  Role changes at all: 0 / 6 / 22 / 115 of 250 trials
  (σ 0.01 … 0.1); at σ = 0.1 mostly channel → axis (L→A 49, B→A 35), and those end ≥ 11, typically (L, B, axis) = (20, 21, 69)
  or (22, 19, 69): two channel squares go axial and complete a full line.  The census sampled one count pair; per-coordinate
  Gaussian heat breaks structure but can't make the coordinated move a transfer needs, so the count-pair table must be built.

## Count-pair moves (10-05 evening; `movegen.py`, `runs/mv0`–`mv2`)
* Moves on each band's own lattice (add / del k sites, transfer k at the junction, whole lattice row add / del, rigid shift,
  rotate, drop the junction squares), 1–2 per proposal, then evaluation, then slp2.
* Evaluation 'fill' (ftmc.evaluate: exact axis refill around the edited bands) fails: nearly every real edit needs the fill side
  to grow to 11.37–11.53 before n squares fit, and the squeeze back collapses into a full line (18/30 no result, rest ≥ 11).
  The record's axis squares are staircased exactly around the bands; an exact refill around a slightly moved band can't
  reproduce that.  Evaluation 'keep' (keep the start's axis squares, drop the most-overlapped or insert at least-overlap
  positions, penalty squeeze from the overlapping state) works: 224/300 evaluated, 215 without full lines, 106 below 11.
* (Bug fixed on the way: ftmc.load is degrees, slp2 / rigid.load radians; the first smoke run fed degrees to slp2.)
* **New count pairs below 11** (`runs/mv2`, 300 proposals from rec110, 5.5 min wall):  (L, B, axis, other) =
  (23, 21, 66, 0) at **10.998202** (×13; `add L 1`: one more tilted square in the left band, one fewer axis square);
  (22, 20, 68, 0) at 10.99994 (×1, `del B 1 far`); (22, 21, 66, 1) at 10.99782 (an axis square tilts).  (19, 24) and (21, 22)
  rows are the junction 46° artifact (squares at 40–50°), not transfers.  Above 11 but close: (24, 21, 65) 11.00132,
  (21, 21, 68) 11.00136, (22, 22, 66) 11.00165, (19, 22, 69) 11.00055.
* Exact status: (23, 21, 66) is first-order jammed in the corner–corner MILP (dS = 0) and pipeline.py finds no descent, but
  exactsolve can't form the KKT system ("no equilibrium with smooth contacts": the near-parallel gap, as for the .99841
  census classes).  Not certified yet; f64 valid with 1.8e-3 below 11.  (22, 20, 68): MILP dS −5.6e-11 (noise), same gap.
* So the record funnel's count is not rigid: a 23-square left band also packs below 11.  Next: iterate (moves from the best
  of each count pair), and fix exactsolve's near-parallel gap, which now blocks certifying new shapes.

## exactsolve fixes, re-solve, beam round 1 (10-05 late)
* exactsolve (`250ad0f`): (1) KKT Jacobian singular near convergence → re-choose the independent contacts at the current
  point (Newton had reached 3.7e-24 then diverged; cond 1.6e5 → 9e16 along the iteration); (2) jammed only with
  corner–corner disjunctions → one branch's corner–corner incidences as equations (every branch is jammed, so each has an
  equilibrium) before the near side–side ladder.  Test set identical.  Register n = 105 (fix 2) and n = 130 (fix 1) now
  certified KKT local minima, 2.09e-11 and 1.30e-11 below Couzo's values; batch 323/324, open n = 292.
* (23, 21, 66) s(110) packing: certified, S' = 10.998202452206362520354…, jammed in every branch; λ ≥ 0 not established in
  the chosen branch (bound + first-order jam only).
* cen7 re-solved: unresolved 53 → 12; certified 853, **21 distinct** sub-11 minima (new singletons 10.99756849…,
  10.99879066…), Chao1 30, f1/N 0.007; the .99841239 / .99841389 classes have 29 / 40 hits.  `data/cen7` updated.
* Beam round 1 (`beam.py`, `runs/mv3`: best of 12 count keys × 40 moves): 327/480 evaluated, 16 count keys (L, B, other;
  role windows 12–40 / 50–78, ~45° junction squares = other) below 11.  New: (19, 19, 4) 10.997163 (better than all but the
  top three census minima), (19, 21, 3) 10.997841, (19, 20, 4) 10.997964, (20, 20, 4) 10.998383, (18, 20, 4) 10.998384,
  (22, 20, 0) 10.999920.  Nothing below the candidate.  Most new keys restructure the junction (squares at ~45°).

## Beam round 2, and reverse connectivity from above 11 (10-05 late; `runs/mv4`, `runs/rev1`)
* Beam round 2 (16 keys × 40): 403/640 evaluated, 24 count keys below 11, nothing below the candidate.  Exact (both certified,
  jammed in every branch): **(L, B, axis, other) = (19, 19, 68, 4) at 10.99688704566278565…** (strict modulo 11 flat motions;
  a 4-square ~45° junction, 1.0e-4 above the candidate) and **(24, 22, 64, 0) at 10.99929042749055536…** (46 band squares).
* Above 11 so far (census + mv2/mv3): 696 trials ending ≥ 11 (349 with a full line); **207 distinct jammed no-line minima**
  (f64 side to 1e-7), 160 singletons, spread 11.0003–11.25; census and moves share 1 of them.  Badly undersampled.
* Reverse (`runs/rev1`, 34.8 min wall): 40 starts = the 30 above-11 minima nearest 11 + 10 in 11.005–11.1.
  Perturbation (σ 0.01/0.03/0.05 × 6, census3 multi-ref): of 720, 197 same minimum, **454 a different ≥ 11 minimum**, 38 full
  line, **31 below 11 (4.3 %)**, from only 8 of 40 starts (one, 11.0010975, returns 17/18: a shallow shoulder of the funnel).
  Even σ = 0.01 leaves the start 63 % of the time (from the record: 0 %): above-11 basins are small, the web between them
  dense.  Moves (10/start): 19/400 below 11, 8 starts.  Distance to 11 does not predict return (11.00034: 0/28; 11.0154: 2/28).
* Where returns land: 16 exact classes, all in the upper sub-11 ladder (10.99756–10.99969); **never the top four** (candidate,
  record, .99718876, .99719765) in one step.  6 classes not in cen7 (10.99884139…, .99895628…, .99919713…, .99930379…,
  .99952510…, .99968978…): sub-11 exact minima seen so far: 27 (+ the move-found count pairs).
* Reading: below 11 a smooth funnel (few minima, large basins, best basins most common); above 11 a rough plateau (many small
  basins, mostly connected to each other); a few gateway minima near 11 lead into the funnel's upper rim.  For the s(90)
  needle: the search has to find a gateway; once in the rim, cen7 says descent to the bottom is easy.
* **Caveat on "saturated" (Evan, 10-05):** Good–Turing f1/N bounds the total *measure* of unseen basins under this sampling
  (~0.7 %), not their number or depth.  A narrow optimal basin of measure ~1e-4 is consistent with all our data; claims about
  the extreme member of the set carry large error bars.  "Saturated" means only: basins with appreciable measure here.

## Junction sub-funnel, crop pass test, first s(90) runs (10-05 night; `runs/jx1`, `crop1*`, `pass1`, `c90*`)
* Junction family (`runs/jx1`): census from (19, 19, 68, 4) (σ 0.003/0.01/0.03 × 40): 120/120 below 11, 119 certified,
  **20 distinct minima** (Chao1 25), none in cen7; best **10.99688698684887…**; 10.99689–10.99692 cluster 1e-4 above Couzo's;
  commonest 10.99752 (40), 10.997157 (30).  Moves from 12 family starts: nothing below 10.996887.  A second sub-funnel one
  junction edit from the record's, invisible to the record census: the sub-11 set is ≥ ~47 exact minima and growing.
* **Crop pass test (132 → 110).**  `crop.py` on Couzo's s(132) record (11.9913; bands 18 left / 31 bottom, the bottom 4 rows
  thick): 600 crops, 214 jammed no-line, best 11.00438; the best cuts thin the bottom band to 3 rows (strip at y 3.6–4.4).
  Beam round 1 (480): best 11.00195; round 2 (480): **one sub-11, 10.998566** ((L, B, other) = (20, 21, 2)).  ~2 CPU-h, no
  s(110) data used.  Census from it (60 trials, exact): 42 stay, 6 → 10.99919, 12 ≥ 11; none descend toward Couzo's.  So the
  constructor + moves reach the funnel's upper rim (rare: 1/960 moves), and the rim minimum is a shallow isolated basin.
* **s(90) crops** (k = 10, n = 90; 120 each from rec110, couzo110, (19,19,4), (24,22,0), (23,21,66)): 231 jammed no-line, 170
  distinct, none below 10; best 10.018024 (from (23,21,66)), 10.019334 (from (19,19,4)); rec110 10.0305, couzo110 10.0347.
  Crop-stage gap to k: 0.018 at 90 vs 0.0044 at 110.  Beam rounds running (`runs/c90b–d`).
* Ops: beam.py had no --procs passthrough (one launch oversubscribed to 23 processes for ~2 min; killed, relaunched).

## s(90) first campaign: stop and regroup (10-05 night; `runs/c90*`)
* Pipeline = the 132 → 110 pass test's: crop (600, 5 sources) → beam ×3 (480 each, top 16 count keys) → one left-band-focused
  round (640: `--focus L`, moves uniform / rebuild / generic).  ~3700 trials, ~8 CPU-h.
* Best per stage: crop 10.018024 → 10.012239 → 10.012133 → **10.009566752689652885…** ((L, B, other) = (10, 18, 0); exact:
  strict local min modulo 21 flat motions, jammed in every branch, certificate valid; `data/s90_best_10.0095668.*`) → focused
  round no change (uniform best 10.01038, rebuild 10.01312: clean-lattice left bands lose to the improvised ones).
* Shapes: the right architecture (bottom band 60–65° ~3 rows along the bottom wall to the right wall; left band 20–35° up the
  left wall; transversal; axis block top right), certified local minima at 10.018–10.028 (`figs/s90_first_minima.png`).  Waste
  sits in the left band (shorter, multi-angle, stair-step gaps) and the junction.  The beam shortened the left band to ~10.
* Reading: no evidence either way on s(90) < 10.  At 110 the same pipeline needed 2 beam rounds to close a 0.0044 crop gap
  and reached only the rim (1/960); at 90 the crop gap is 0.018 and 3 rounds closed half of it, decelerating.  A greedy
  best-per-key beam is the wrong instrument for a 0.01 gap: it cannot accept uphill steps.
* **Regroup topics** (Evan: stop here, think about mixing and basin hopping):
  1. Basin hopping with Metropolis over minima (moves = movegen kinds + σ-kicks), then PT over temperature.  Scale from data:
     neighbouring minima differ by 1e-4–1e-3 inside a funnel; the above-11 plateau's basins are small (σ = 0.01 leaves 63 %).
  2. Calibrate rim → bottom at 110: from the crop-derived rim minimum (`runs/crop1c/m_1200238.txt`, 10.998566) can moves
     (not σ) reach Couzo's?  Cheap, and it is exactly the step s(90) would need after a gateway.
  3. Gateway fraction at 110: what share of above-11 minima return below 11, and a cheap predictor.
  4. Constructors for k = 10 other than cropping 110 (other two-band records, layout genomes seeded with the s(110) band
     shapes), since the crop gap grew from 0.0044 to 0.018.
  5. Content-agnostic swaps (pair exchanges, cluster rotations) alongside the band-aware moves.

## Minima PT (10-06; `mcmin.py`, `sweep.py`)
* **Caveat (10-07, with Evan):** neither minima PT nor kick+quench hopping samples exp(−s/T): polish is a many-to-one descent
  (zero reverse density, so not a valid MH component), and kick+quench is not symmetric between minima (P(B|A) != P(A|B) when
  basins differ in size).  Fine for the goal (find the tail, not integrate over it); acceptance/swap rates, censuses and
  detection rates are properties of the walk (protocol-relative), used to understand and tune it.  What they cannot show is
  where the walk never goes: only external known-better answers (SQUISH set, rediscovery pairs) test for blind spots.
* `mcmin.py`: Metropolis over quenched minima (E = side), replica exchange across a geometric T ladder (side units).  Moves:
  kick, local kick, **reinsert** (remove a random square, re-add at a hole picked with weight clearance^α, α = 2), cluster
  rotation, angle swap, band (movegen); 1 + Poisson(0.25) moves per proposal (Evan).  Quench = loosen 2 % + soft squeeze +
  slp2 (budget 60 s).  Epochs are time-based (60 s per replica), then R swap sweeps (swaps are free; lesson from an earlier MH/PT project).
* Smoke `runs/pt0` (12 replicas from the rim minimum, 1 epoch): 4–8 s per proposal, some 20–40 s; reinserts mostly land
  at 11.0 with full lines or 11.04–11.06; kicks ≤ 0.04 return to the rim minimum.
* **Pre-registered calibration (written before the runs):**
  - Run A `runs/pt1`: 12 replicas, all from the crop-derived rim minimum `runs/crop1c/m_1200238.txt` (10.998566, no s(110)
    data in its construction), T 3e-5 … 1e-2 geometric, 2 h wall on 12 processes (24 CPU-h).  **Success:** a line-free
    minimum ≤ 10.99680 (record or Couzo class).  **Partial:** ≤ 10.99720 (top four of cen7).  Prior (mine): partial ~50 %,
    success ~30 %.  Census from this rim minimum (60 trials) never went below 10.99919.
  - Run B `runs/pt2`: same settings from our own funnel 11.0075618 (`runs/cen3/our_0.03_4.txt`).  **Success:** any line-free
    sub-11 minimum; full success ≤ 10.99680.  Prior: sub-11 ~25 %, full ~5 %.  Census from it never went below 11.
  - Also report: acceptance per move kind and per T, swap acceptance per pair (and per walker), excursions to sub-11 /
    ≤ 10.99720 with N_eff = N / (1 + CV²), distinct minima visited.  Round-trip counts are diagnostic only (grazes).
* **Neighbour sweep n = 30–200** (`runs/sw1`, 76 open non-integer n, 8 × (n+1 minus a random square) + 4 kicks (σ 0.01, 0.03)
  each, 912 quenches, ~1.5 h on 4 processes): **nothing**.  Kicks return to the record's own minimum (the "below register"
  values, −3e-12 … −1.4e-11 at n = 68, 102, 103, 106, 110, 132, 154, 177, 182, are analytic polish of the same packings, as in
  the exact batch); n+1 minus a square never reaches s(n).  Records here are robust to cheap probes.  n = 201–323: `runs/sw2`.
* **Run A result (`runs/pt1`, 2 h, 12 replicas, 9230 proposals, 16.7 CPU-h): not met (neither partial nor success).**  Best
  10.9975673 (cen7 class .99756727), from 10.998566.  Hot half (T ≥ 2.5e-4) sat on full-line states at E = 11.0 exactly
  (acceptance 0.92: a flat, entropy-rich grid plateau); swap acceptance 0.40 at pairs 0–1 and 3–4; excursions below 11 N_eff
  5.6.  The grid is an entropy sink for hot replicas: "temperature is the wrong coordinate" (lesson from an earlier MH/PT project).
  Per kind (all T): reinsert 92 % full lines, kicks 85 % (mostly from replicas already on the grid).
* **Reinsert A/B (`runs/ab1`, 4 starts × 11 variants × 30):** every reinsert variant is destructive (37–97 % full lines, the
  rest median +0.02 … +0.3); pose-space best-pose holes are *worse* on lines than clearance holes (the best pose is often
  axis-aligned at a wall); load-weighted removal worse; loosen 0.5 % / 1 % no better; 1/1080 reached a better sub-11
  minimum.  kicksym (equal corner metric) = kick.  Reading: inserting a square anywhere in a jammed packing costs ~0.4
  penetration, and a global quench pays for it with a large rearrangement that snaps rows.  Fixes would need a local
  rearrangement (vacancy → hole chain shift) rather than a better hole.
* **Run B deviation (decided after A and the A/B, before B started):** full-line states rejected (line_pen 1e9), reinsert
  weight 0; otherwise as pre-registered (from 11.0075618, 12 replicas, T 3e-5 … 1e-2, 2 h).  Criteria unchanged.
* **Neighbour sweep n = 201–323 (`runs/sw2`, 6 rm + 4 kicks per n, 600 quenches): two new best-known packings.**
  Live register checked 10-06 (jlevy.github.io/squares/cases/<n>.html):
  - **s(266) ≤ 16.82303560648544271826283801102** (cert S'; KKT point 16.8230356064854427180946…), register 16.82306208283780
    (Ellsworth 2024): **−2.65e-5**.  Strict local min modulo 52 flat motions, jammed in every corner–corner branch.
  - **s(270) ≤ 16.937807228446029110709567408515** (KKT 16.9378072284460291105…), register 16.9378103293409541102 (Couzo
    09-27, our 10-05 exact refinement): **−3.10e-6**.  Strict local min modulo 28 flat motions, jammed in every branch.
  Both from a single kick of the record (σ 0.01 / 0.03) + quench; certificates VALID under verify_cert.py and verify_cert2.py
  (`../exact/results/sw2/`, copies in `candidates/`).  Not posted anywhere (ask Evan per post).  2 of 60 n improved with 4
  kicks each: large-n records are lightly polished.  Deeper kick census (18 per n, σ 0.003/0.01/0.03): `runs/sw3`.
* **Run B result (`runs/pt2`, 2 h, 12 replicas, 6525 proposals, 16.3 CPU-h, full lines rejected): not met.**  No sub-11 minimum
  (606 distinct line-free minima visited).  Best below the start: 11.0025675 (line-free), and a **line-free trap at exactly
  11.0** (68 visits; the coldest walker sat there for the last ~30 epochs): 78 axis squares, 32 at 19–22°, no full line by
  `full_lines`, but side exactly 11: a staggered vertical/horizontal chain of axis squares (consecutive x-intervals overlap)
  forces height ≥ k just like a full line.  `full_lines` under-detects the grid attractor; the rigorous obstruction is "a chain
  of k axis squares with consecutively overlapping projections".  Mid-T replicas (1e-4 … 2e-3) stayed in the 11.0076 funnel:
  its exits land at ≥ 11.03, rarely accepted.  Swap acceptance 0.25 at pair 0–1 (the 11.0 trap).
* Grid escape from run A's hot replicas: 7237 proposals from E = 11.0 states, ~2 % land line-free but all far above 11 (never
  accepted, even at T = 1e-2), 0 below 11.  Under our moves the grid is effectively absorbing; the useful measure is the one
  conditioned on "no grid state", plus per-basin absorption rates (Evan's rejection-counting idea; flux balance
  π(G) e_G = π(¬G) r̄ gives only a lower bound on π(G)/π(¬G) when e_G ≈ 0).
* **Kick census n = 201–323 (`runs/sw3`, 18 kicks per n, σ 0.003/0.01/0.03, 1080 quenches): nothing new.**  266 again
  (16.8230358, a neighbour of our 16.8230356); 270's improvement *not* re-found in 18 kicks (detection is far from certain at
  δ = 3e-6, distance 0.03).
* **Staggered axis chains** (`layout.axis_chain`, `83fca46`): `full_lines` now also flags a chain of ≥ k axis squares with
  consecutively overlapping projections (forces side ≥ k).  Run B's 11.0 trap: chains of 11 both ways, no straight line.
  451 saved states: no sub-11 state has a chain ≥ 11, no chain exceeds the side.
* **Rediscovery calibration (`rediscover.py`, `runs/rd1`):** 52 (older → newer) record pairs = every n where the live register
  beats the 09-30 site mirror (50, mostly Couzo Sep 2026) + register → ours at 266, 270.  Matching distance (best of 8
  container symmetries, Hungarian, corner metric): 0.003 … 0.62; δ 3e-8 … 2.2e-2; larger δ tends to sit farther.  Protocol:
  10 kicks per σ ∈ {0.003, 0.01, 0.03, 0.1} from the old packing; detected = side ≤ new + 1e-7, no grid obstruction.
* **Rediscovery result (`runs/rd1`, 51 pairs × 40 kicks, 2040 quenches, ~1.5 h on 15 processes):**
  - Detection (new record reached from the old one) falls fast with matching distance: dist < 0.02: 7/10 pairs found (often
    at most σ); 0.02–0.06: 4/21 (207, 177, 272, 236); > 0.06: 1/20 (106, at σ = 0.1).  δ matters less than distance: tiny-δ
    pairs at small distance (68, 199, 269) are found, large-δ far pairs never.
  - Improvement over the old packing (not necessarily to the new) is common even when detection fails: 31/51 pairs improve in
    ≥ 1 kick; kicks find intermediate minima between old and new records.
  - Never detected at any distance: 259, 156, 130 (small dist, δ ≤ 3e-6: tight neighbours that kicks of σ ≥ 0.003 overshoot?),
    and our own 266/270 from the register (σ grid differs from the sweep's lucky kick: detection at δ ~ 1e-5, dist 0.03 is rare).
  - **New best-known packing on the way: s(272) ≤ 16.968110145769600505005359235592** (cert S'; KKT 16.9681101457696005048…),
    register 16.968165867852158 (Couzo 10-03, our 10-05 exact polish): **−5.57e-5**.  One kick (σ 0.03, trial 8) of the 09-30
    mirror's older packing (16.9697160, Ellsworth/Casson era).  Strict local min modulo 40 flat motions, jammed in every
    branch; certificate VALID under both verifiers (`../exact/results/sw2/rd1_n272.*`, copy in `candidates/`).  Live register
    checked 10-06.  Lesson: kicking *older* records (other basins of the same family) finds things kicking the current one
    doesn't.
* **Deep census 266/270/272 (`deepcensus.py`, `runs/dc1`, 1000 quenches, slp2 budget 180 s, ~1.3 h on 15 processes):**
  starts = 09-30 mirror (270, 272), register, ours; σ 0.003 … 0.05 × 25; every output saved with distances.
  - **s(266) ≤ 16.823028750756475970708783711254** (KKT 16.8230287507564759705…): −6.86e-6 below this morning's, −3.33e-5 below
    the register; 9/125 kicks from ours reach it at every σ (d ≈ 0.008): a robust neighbouring basin, not a lucky hit.  Strict
    local min modulo 94 flat motions, jammed in every branch, certificate VALID under both verifiers (`dc1_n266.*`).
  - With the 180 s budget, kicks from the 266 register record reach this morning's 16.8230356 at σ ≥ 0.01 (rd1 at 60 s: 0/40):
    the descent budget changes detection at large n (13 % of 60 s quenches at n ≥ 250 were unconverged).
  - 270, 272: nothing below ours.  Landscape near the records is rich: f1/N 0.4–0.8 at σ ≥ 0.01 (most minima seen once),
    median travel 0.015–0.03 from the start.  At 270 the mirror's family stays ≥ 1.6e-3 above (a different funnel).
* **Cold chains from our bests (`runs/ch_266/270/272`, mcmin 4 replicas each, T 1e-7 … 3e-5, kicks/lkick/crot/aswap, slp2
  budget 180 s, 1.5 h, ~14 CPU-h): nothing below.**  Only 158–292 proposals per n (quench 40–110 s at n ≈ 270), 47–136 distinct
  minima each.  Total extra compute on the three: ~1000 census quenches + ~700 chain proposals (≈ 10× the original per-n
  effort).  Submission list: s(266) ≤ 16.8230287508, s(270) ≤ 16.9378072284, s(272) ≤ 16.9681101458.

## fq: smooth lifted quench + face-branch SLP polish (10-07; `src/fq.rs`, `bench_quench.py`, `runs/bq0`, `runs/bq1`)
Motivation: SQUISH (itsnaka, jlevy#401/#422: 23 n in 88–303, gains up to 1.6e-2) wins by volume (neighbour/graft seeds +
basin hopping + polish); our quench (soft squeeze + slp2) costs 4–8 s at n = 110, 16–40 s at n = 270 (~500 HiGHS LPs).
* **Formulation.**  Each nearby pair gets its own separating line (direction phi, offset d from the centres' midpoint);
  constraints "4 corners of i on one side, 4 of j on the other" + corners inside the walls, objective s.  Each constraint is
  smooth, so the PHR augmented Lagrangian is C^1 (no SAT max-over-axes kink).  Gradient checked (rel err 3e-8).
* **ALM alone is not enough:** linear multiplier convergence; at mu >= 1e5 L-BFGS is ill-conditioned and stops feasible but
  not stationary (|grad| 1e-2; e.g. 8e-5 above the slp2 minimum at n = 270).
* **Lifting does not remove the corner-corner disjunction** (my first guess): at a corner-corner touch, rotating the line
  changes both corners' rows equally and oppositely, so phi enters at second order; the LP sees one half-plane (a strict
  subset of the true union of two face branches) and crawls.  Fix = incidence model inside the lifted one: before each LP,
  snap every separator to a face (largest gap), phi tied to the owner's theta (column substitution + tie at trial points);
  on a stall, flip ambiguous pairs (>= 2 non-parallel faces within 1e-7 of contact; real in axis blocks: diagonal
  neighbours touch at corners, ~260 at n = 270), candidates ranked by LP row duals (top 8, first improving flip).
* **Polish** = trust-region SLP (Clarabel 0.11.1 interior point, pure Rust; variables scaled to the box, one elastic
  slack, l_inf merit s + 3s·viol, second-order correction when rho < 0.75), from the ALM point at mu <= 1e2.  ~10–25 SLP
  iterations, ~40 ms per LP at n = 270 (only columns that occur in near rows).
* **Head-to-head (`runs/bq1`, same kicked states, n = 110 / 270, sigma 0.01/0.03/0.1 × 5):** fq 0.3–5 s at 110, median
  ~3 s at 270 (outliers 12–45 s at sigma 0.1) vs slp2 pipeline 1–10 s / 16–186 s.  slp2 run on fq's output changes nothing
  (< 1e-10) in 22/30; the rest gain 3e-7 … 6e-4 (missed flips).  Basins differ from the slp2 pipeline's (different descent
  path); neither is better on average; at 270 sigma 0.01 fq lands in the 16.93797x family (near our record) 4/5 vs slp2 0/5.
  At 110 the result is independent of mu cap / skin (same 6 minima for every setting).
* Next: hopping driver on the SQUISH set (benchmark: from each pre-SQUISH register packing, CPU budget per n, reach their
  side?); speed (ALM evals dominate at mu <= 1e2 only for big kicks; LPs ~40 ms each).

## SQUISH benchmark v0, and how much polish per step (10-07; `hop.py`, `polish_value.py`, `runs/hop0`, `runs/pv0`, `runs/pv1`)
* **hop0:** fq basin hopping (global kicks 0.003–0.03 / local kicks r 1.5–3.5, sigma 0.05–0.2; Metropolis T 1e-5; reset to
  best after 40) from each pre-SQUISH register packing (batch inputs, 10-05), 10 min per n, 23 n of SQUISH's table: **nothing**
  except n = 126 (−1.9e-4, 13 % of SQUISH's gain).  Matches SQUISH's lineage (their big gains are new seeds: (n+k) − k from
  neighbours, grafts of small-n records, carving), not local search from the record.  Median quench 0.3 s (108) … 28 s (263,
  303).  Loosen 1.02 throws some records out of their basin before the first step (154: +1.0e-3, 238: +6e-5).
* **Polish value (Evan: are we over-polishing?).**  Per hop proposal: ALM-only side, side after k accepted SLP iterations,
  full; 240 proposals at n = 108–263.  Polish = 88–99 % of quench time.  |s_k − s_final| median 7e-5 (k = 0), 2e-6 (4),
  1e-8 (8), 5e-11 (12); 90 % tail ~1e-5 until the flip search (up to 3e-2 in rare cases).  Cost k = 4: 0.28, k = 8: 0.40.
  Decision discrepancy D(k, T) = mean |a_T(s_k) − a_T(s_final)| over same-n pairs: T = 1e-3: k = 3 → 0.02; T = 1e-4: k = 8 →
  0.015; T <= 1e-5: floor 0.03–0.06 until the full polish.  Lens (Evan): polish depth = energy-evaluation fidelity, so effort
  should scale with T; for greedy/cold chains use delayed acceptance (Christen & Fox 2005): screen at k ~ 4–8, full polish
  only within a few × the k-error of the current state.  A lazy rule on the ALM-only side alone saves ~2× but misses a
  quarter of the good basins (130: 17/24 at margin 1e-4).
* **Neighbour seeds (`runs/hop1_108`, `hop1_179`; 7 chains × 8 min, loosen 1.0):** 108 from record(110) − 2: seeds
  10.960–10.987, hopped to 10.9435 (register 10.9259, SQUISH 10.9099).  179 from Couzo's 180 − 1: seeds 13.922–13.924, hopped
  to 13.9167 (register 13.8953, SQUISH 13.8838, i.e. 4.4e-2 below the 180 it came from).  Our quench after a removal shrinks
  ~6e-3; the vacancy is never absorbed by a global rearrangement.  Guess at SQUISH's "shrink-hopping": fixed-side feasibility
  search (target just below current, overlap-only energy, hop until zero overlap, shrink again); overlap pressure spread over
  the packing drives vacancy propagation, which min-s descent does not.  Our old `packer hop` (descending target) was only ever
  run from random starts.  Next test.

## Search benchmark, delayed acceptance, identification (10-07; `bench_search.py`, `bench/seeds/`, `runs/bs_*`, `runs/pv2`)
* **Harness:** fixed seeds s108 = record(110) − 2 (10.995103), s179 = Couzo 180 − 1 (13.923576), r126/r237 = register
  records (fq `--no-alm` polish; fq's ALM from a feasible jammed start can leave the basin: 126 +1.8e-4).  5 replicates ×
  480 s; metric = median gain / SQUISH's gap at 60/120/240/480 s.
* **Waste in the walk (hop0/hop1 logs):** 18–82 % of full quenches return to the current minimum (sigma 0.003 kicks: 72–82 %);
  most seeded chains were still improving at the end (sample-limited).
* **Variants (frac of SQUISH gap at 480 s, medians; s108 / s179 / r126 / r237):** base (loosen 1.0, full quench every
  proposal) 0.01 / 0.03 / 0 / 0; l102 (loosen 1.02) 0.08 / 0.18 / 0.48 / 0; da (delayed acceptance: screen = ALM + <= 8 SLP
  iterations, no flips; reject if > 1e-4 above current, "same basin" if within 3e-7, else full polish + Metropolis) 0.05 /
  0.20 / 0.13 / 0 with 3–4× the proposals; da102 0.08 / 0.16 / 0 / 0.  Both changes ~5–10× over base, not additive (loosen
  1.02 acts as a useful global move).  Replicate spread ~10× (5 reps cannot rank l102 vs da).  r237 never moves.  One l102
  replicate reached 11.77373 at 126 (70 % of SQUISH's gap).
* **Identification (Evan's corner-distance idea; `fq --ident`, `runs/pv2`, 60 traced quenches at 130/180):** load-bearing
  contact network = support of the interior-point (max-support) LP duals, at contact level (pair + canonical face, corners
  dropped: which corners of a face contact carry load drifts along flat motions); ~2.1–2.4 contacts per square at n = 130;
  dual magnitudes have a clean gap (0–5 vs 10–11 decades below max).  Network settles by SLP iteration ~7–9 (s within 1e-7).
  Stop rule "network unchanged for 3–4 accepted iterations": cost 0.37–0.41 of full, |err| median 4e-9–5e-10, 90 % 3e-6, max
  4e-5 (fixed k = 8: 0.40, 1e-8, 1e-5, 3e-2): better tail at equal cost; ~10 % still gain > 1e-6 later via branch flips
  (stable under the current branch), which only a stage-2 polish catches.  Late drift: a few contacts drop out of the
  network after s has converged, so basin matching needs set similarity + side, not an exact hash.

## Breadth: basin explorer at s(110) (10-07; `explore.py`, `explore_exact.py`, `runs/ex_smoke`, `runs/ex0`)
Evan: more explore, less exploit; the sub-11 set at 110 was 48 certified minima (cen7, jx1, rev1, pass1, pil7;
`runs/known110.json`) in ~2 families.  Metric: distinct certified sub-11 minima and structural groups, not best side.
* **Explorer:** quality-diversity archive of polished basins (MAP-Elites style), descriptor (ntilted/2, tilt-angle classes,
  below k), parents uniform over cells (top 8 per cell, weight 1/(1+expanded)^2), stepping stones up to k + 0.05, grid
  states discarded; moves = mcmin's kick/kicksym/lkick/crot/aswap/band/reinsert + bigkick (sigma 0.05–0.2); quench = fq screen
  (<= 8 SLP iterations, no flips, loosen 1.0/1.02/1.05) → return if within 1e-6 of an archived side, else full polish
  (--no-alm --loosen 1.0).  Starts: rec110, couzo110, junction (19,19,68,4), crop rim 10.998566, our 11.0076 funnel.
* **Bugs fixed on the way:** fq defaults to loosen 1.02, so `--no-alm` polishes without an explicit `--loosen 1.0` polished a
  2 %-loosened copy (record → .99718876); fq flags now take the last occurrence (`hop.quench` prepends EXTRA).  **The
  `da102` benchmark variant's stage 2 was affected** (polished from a loosened copy, not the screened state): read da102 as
  "screen + loosened re-quench", not as delayed acceptance at 1.02.
* **Smoke (3 min, 15 procs):** 81 sub-11 basins; exactsolve: 68 certified distinct minima (1 archive entry each), 12 not
  minima, 1 unresolved; **40 new → the known sub-11 set goes 48 → 88**; 14 role-count groups (L, B, axis, other), e.g.
  (19,24,67), (22,22,66), (23,19,68), (25,19,66) beside (22,21,67) and (19,23,68).  New near-bottom trio 10.99682777 /
  .99682787 / .99682933 (all certified, distinct; ~4.5e-5 above Couzo's).  Nothing below the record.  Still two-channel
  structures: groups are variations within the family, not new funnel types.
* **ex0 (60 min, 15 procs, same starts; 22.6k proposals: 7.6k returns, 8.7k discarded ≥ 11.05 or grid, 5.5k new basins):**
  268 sub-11 basins → exactsolve: 205 distinct certified minima (54 not minima, 8 unresolved; 1 f64 split), 172 not in the
  known 48.  **Union of certified sub-11 minima at 110 (known + smoke + ex0): 229** (was 48).  28 role-count groups.
  Discovery curve (distinct certified by t): 68 / 88 / 106 / 126 / 149 / 170 / 193 / 205 at 5 / 10 / 15 / 20 / 30 / 40 / 50
  / 60 min: ~2.3 per minute and linear after the first 15 min (no sign of saturation).  Best: still Couzo's
  10.99678339663; closest new 10.9967871155 (+3.7e-6; (22,21,67,0)), then .99682787, .99682933, .99683259, .99686519.  New
  minima by move: lkick 48, kick 47, bigkick 41, aswap 41, kicksym 32, crot 16, band 9.  ~20 % of sub-11 archive entries are
  not minima (screen + polish lets saddles/unjammed through): certification stays in the loop.  Not resumed (no record;
  Evan: resume only on a record).
* **ex1 = ex0 + screen-stage discards (screened side > k + 0.05 + 1e-3 or grid-obstructed → no full polish), 8 procs ×
  30 min (240 proc-min):** 10.9k proposals (2.7× ex0's per proc-min; 2174 screen-grid, 1311 screen-discard), 218 sub-11
  basins → 156 certified distinct (ex0 at equal proc-min: ~112): **~1.4× certified discoveries per CPU**.  **Union of
  certified sub-11 minima (known + smoke + ex0 + ex1): 275** (`runs/known110_all.json`).  No record.
* **Shrink-hopping descent (`fq --shrink`, fixed-side overlap-only relaxation, grow to feasible then shrink by d, halving on
  failure; then polish):** on 40 explorer proposals from ex0 basins, vs the ALM quench: same basin 20/40; where they differ
  ALM lower 14, shrink lower 6; new vs ex0 archive: ALM 6, shrink 5; shrink 1.6× slower.  A different map (keeps as a
  diversity arm) but not more efficient per CPU; not the SQUISH secret by itself.

## Per-arm yields, new moves, adaptive selection, new starts (10-08 night; `arm_yield.py`, `runs/ex2`–`ex4`)
* **Instrumentation:** `proposals.jsonl` per run (kind, outcome, CPU s, new archive index); `arm_yield.py` = new certified
  sub-k minima (not in a reference set, credited to the kind that first produced the archive entry) per CPU-hour.
* **New coordinated moves:** rowslide (a contiguous run of axis squares slides 0.2–0.6 along its row/column, sometimes one
  end only), chainshift (remove a square, shift the chain on the segment to a clearance hole one step toward the vacancy,
  re-insert at the hole), mirror (reflect a local cluster across h/v/diagonal through a random point).
* **ex2 (8 procs × 30 min, all arms, fixed weights, ref = 275 known):** 13.7k proposals, 3.98 CPU-h, 140 distinct certified
  sub-11 (52 new vs 275).  New per CPU-h: kick 33.5, bigkick 21.7, kicksym 17.9, aswap 17.4, **rowslide 12.1**, band 7.7,
  lkick 7.7, **chainshift 2.8, mirror 0**, crot 0, reinsert 0.  Discard share (screen-grid / ≥ ceiling) tracks the
  failure: chainshift 78 %, mirror 79 %, reinsert 95 %, crot 65 %.  Reading: inside the funnel small perturbations (random
  or coordinated) make new minima; destructive moves leave the funnel (the quench lands ≥ 11 or on the grid).
* **Seeds from neighbours:** records 111–114 minus 1–4 squares → all exactly 11 (grid; those records are grid-like).
  Records 107–109 plus 1–3 squares at clearance holes, loosen 1.02–1.08 → 11.004–11.4; 6 distinct ≤ 11.05 kept
  (`seeds110b/`): new structures (21,0,89,0), (1,10,99,0), (4,1,89,16), (7,4,98,1), (4,14,83,9) + one pure grid at 11.0.
  (`seeds110/` holds 3 useless 111/114-derived files; left in place.)
* ex3 = ex2 starts + `--adapt` (Thompson sampling over kinds, reward = new below-k archive entry, rate per CPU-s,
  Gamma(1 + r, 60 s + c), half-life 10 min, 3 % floor per kind); ex4 = ex2 weights + the 6 seeds110b starts.  Both 8 procs ×
  30 min, concurrently.
* **ex3 (adaptive kinds, reward = new below-k entry *in this run*):** 148 certified sub-11, only 38 new vs 275 (ex2: 52).
  The bandit piled onto kicksym (2327 s of decayed cost vs ~500 for others): kicksym makes many basins new to the run but
  mostly re-finds globally known ones (88 certified, 20 new).  **Lesson: the reward must be global novelty.**  Fixed:
  `--novel-ref` (sorted sides; reward only basins not within 2e-9 of any).
* **ex4 (fixed weights + 6 seeds110b starts):** 146 certified, 38 new vs 275; **the new families produced no archive child at
  all**: ~130 proposals each, ~90 % screen-grid (these 107/108+k structures are 89–107 axis squares; any kick makes an axis
  chain of 11).  Near-grid families are not stepping stones; useful starts need rich tilted structure (cf. the 132 → 110
  crop that reached the rim).  The 11.0076 funnel start: 20 proposals, 13 grid, 8 basins, none below 11.
* Certified union of sub-11 minima at 110 (known + all explorer runs through ex4): **372** (`runs/known110_all.json`).
* ex5 (adaptive kinds, global-novelty reward vs `runs/ref110_ex5.json` = 528 sub-11 sides, certified + archive-level) and ex6
  (fixed weights, same ref, concurrent control): 8 procs × 30 min each.
* **ex5 (adaptive kinds, global-novelty reward) vs ex6 (fixed weights), concurrent, same ref (528 sides):** certified
  distinct sub-11 166 vs 107; **new vs ref 30 vs 13 (~2.3× per CPU)**.  Bandit mass (decayed): bigkick, aswap, band,
  kicksym; reinsert, mirror, rowslide starved.  Control rate 13 new / 4 CPU-h (ex2: 52): the current move set is saturating
  its reachable region (diminishing returns): time for a category step.  Census (certified union): **414**.
* **Crop seeds (`cropseed.py`, fq instead of slp2):** 600 crops of records 129–132 (one unit row + column strip removed,
  count fixed, loosen 1.02/1.05): 517 grid, 0 below 11, **34 distinct stepping stones 11.002–11.047 with new structures**:
  (18,23,69,0), (12,28,70,0), (14,34,61,1), (9,5,70,26), (18,11,81,0), (8,24,68,10), (4,27,75,4), (13,15,70,12), …
  (`seeds110c/`).
* **Frontier budget (`--frontier PREFIX --frontier-share f`):** archive entries carry a lineage root; with probability f the
  parent comes from a frontier lineage (uniform over live roots, then the usual cell pick inside it), else the normal pick.
  Cell-uniform selection gave new families ~1 % each (ex4); a cell bandit would starve them further.
* ex7: standard starts + 34 crop seeds, --adapt (global novelty vs `runs/ref110_ex7.json`), frontier share 0.5, 16 procs ×
  40 min.
* **ex7 (crop-seed frontier, 50 % budget, 16 procs × 40 min, 10.6 CPU-h):** 2 of 34 crop lineages entered below 11:
  c132_74 (seed 11.031 → 6 sub-11 basins, best 10.998351; groups (21,24,65), (21,23,66), (20,24,66), (23,22,65)) and c132_14
  (11.019 → 10.999359, (18,24,68)); the other 32 stalled at 11.0001–11.012.  All 7 certified; **5 already in the census**:
  both outside entries land on the upper rim of the same two-channel funnel (as the 10-05 crop pass test).  Either that
  funnel is the dominant sub-11 structure at 110 or our quench/moves funnel everything into it; not separable yet.
  Whole run: 113 certified distinct sub-11, 13 new vs ref 599 (low per CPU: half the budget on lineages that stall).
  Census: **429**.
* ex8 (adaptive + elite share 0.3: parent from the global top-20 sub-11 basins, weight 1/(1+expanded)) vs ex9 (adaptive),
  concurrent, ref 620; metric adds deep novelty (new certified < 10.9975).
* **ex8 (elite 0.3) vs ex9 (adaptive only), concurrent, ref 620:** certified sub-11 158 vs 186; new 18 vs 23; deep
  (< 10.9975) basins visited 25 vs 15, **deep and new 4 vs 0**.  The elite budget converts breadth into depth (~20 % of
  total novelty for the only deep novelty so far): census work wants elite ~0, record hunting wants it on.  Census: **473**.
* Transfer tests (record hunting on the SQUISH benchmark seeds with the explorer: --adapt, elite 0.3, 8 procs × 30 min):
  `runs/tx179` (s179 = Couzo 180 − 1, 13.923576; SQUISH 13.883795), `runs/tx108` (s108 = 110 − 2, 10.995103; SQUISH
  10.909940).  Hop variants this morning: median gains ≤ 8e-3 (s179) / 7e-3 (s108) per 8 single-proc minutes.
* **Transfer: the explorer as a record hunter (`runs/tx108`, `tx179`; adaptive, elite 0.3, 8 procs × 30 min).**  Best vs
  wall time (1 / 4 / 10 / 20 / 30 min): s108 10.98398 / 10.95506 / 10.91677 / **10.912010** / same (98 % of SQUISH's gap; below
  the 10-05 register 10.925919 and SQUISH's registered 10.920659, 2.1e-3 above their pending #422 10.909940: no claim);
  s179 13.91715 / 13.91651 / 13.90948 / 13.90856 / **13.908471** (38 %).  At equal CPU (first minute ≈ one 8-min hop
  replicate) similar to the hop variants; afterwards the explorer keeps making large gains (4e-2 – 8e-2) where hop chains
  plateaued.  Both flatten after 10–20 min (elite converged).  Next: record starts r126, r237 (`runs/tx126`, `tx237`).
* **Record starts (`runs/tx126`, `tx237`):** r126 → 11.774546701 at 4 min (13 %; the same basin hop0 found), then
  nothing; r237: 759 basins, **all ≥ the record**.  Exploration around a record stays in the record's region, whose bottom
  the record is; SQUISH's record-beaters come from other seeds (their lineage notes).  Test of the recipe (neighbour seeds +
  explorer): `runs/tx237n` from 238 − 1 / 239 − 2 seeds (15.9306–15.9498; record 15.911192, SQUISH 15.903676), 16 procs ×
  30 min.
* **tx237n (neighbour seeds 238 − 1 / 239 − 2 + explorer, 16 procs × 30 min):** 15.930596 → 15.929763 / 15.929178 /
  15.924934 / 15.921655 at 4 / 10 / 20 / 30 min (31 min: 15.921096, (28,38,170,1), via rowslide); still descending ~3e-3 per
  10 min, 1.0e-2 short of the record 15.911192.  ~6 s per proposal at n = 237 vs ~1 s at 110.
* **Profile at n = 237:** screen 0.6–0.9 s (ALM ~0.2 s); full polish 1.6–5 s normally but **300-iteration crawls (18–20 s)**:
  s settled by iteration ~12, then < 1e-8 total over 288 iterations (pred ~1e-10 > the 1e-14·s stop).  **Stagnation stop**
  (`--stag-w 15 --stag-tol 2e-9`, default): s decreased < 2e-9 over the last 15 iterations → stop.  24 proposals at 237:
  53 s vs 137 s (2.6×); |Δs| median 1e-10, max 1.2e-8 (4/24 > 2e-9: occasional archive duplicates; certification dedupes).
  Deployed (release binary rebuilt between runs).  Before/after: `runs/tx237s` = tx237n settings.
* **tx237s (stagnation stop on) vs tx237n:** 3100 vs 2518 proposals (+23 %, not 2.6×); best curves within noise
  (15.924719 vs 15.921655 at 30 min).  91–94 % of CPU is the stage-2 polish of "new?" candidates (median 10–13 s each).
* **Stage 2 at n = 237 is partly a long descent, not only a crawl** (30 proposals from tx237s lineages, screen then
  `--no-alm` polish): default stop 197 s, median 56 it; 1e-8/10 + cap 60: 90 s; no stop 400 it: 315 s, and it ends up to
  **1.3e-3 lower than the default stop and 2.2e-2 lower than the capped one**.  Truncation misjudges basin depth at large
  n, so polish speed is a search-quality problem there, not only throughput.  Suspect: trust-region cap 4e-3 per
  coordinate (= 4 r0).  `--rmax` flag added; test 2e-2 / 5e-2 vs a 400-iteration reference.
* **Trust-region cap is not the bottleneck:** 30 proposals at 237, rmax 4e-3 / 2e-2 / 5e-2: same iterations (median 56)
  and time (225–230 s); larger caps slightly worse (4 vs 1 proposal > 1e-6 above the best).  Default stop matches a
  400-iteration no-stop reference in 29/30 (one genuine long descent, 1.3e-3).  Stage-2 cost at large n = ~56 SLP
  iterations × 2 LPs (~0.07 s each): the lever is fewer iterations (SQP with rotation curvature; Clarabel takes QPs) or
  cheaper LPs.
* **ex10 (adaptive, global novelty, 14 procs × 30 min, ref 697):** 214 certified distinct sub-11, 23 new (3.3 / CPU-h vs
  ex5's 7.5 at a smaller ref): in-family saturation continues.  **Census: 496 certified sub-11 minima at 110** (48 on
  10-07 morning); best still Couzo's 10.99678339663.

### Night summary (10-08, for the next session)
* What worked (measured, concurrent controls where it mattered): quality-diversity explorer (48 → 229 in one hour);
  screen-stage discards (1.4× certified per CPU); Thompson move selection **with a global-novelty reward** (2.3×; with a
  run-local reward it backfired); elite share (only arm with deep novelty, costs ~20 % breadth); explorer as record hunter
  (s108 seed: 98 % of SQUISH's gap in 20 min; s179: 38 %); polish stagnation stop (2.6× on quenches at 237, smaller in the
  explorer).
* What didn't: chainshift / mirror / crot / reinsert moves (destructive: leave the funnel); shrink descent (different map,
  less efficient); 107–109 + k and 111–114 − k seeds (grid-trapped); crop frontier (enters the known funnel's rim, mostly
  known basins); larger trust region.
* Open: (1) **large-n polish = ~56 SLP iterations × 2 LPs** (stage 2 is 91–94 % of explorer CPU at 237): SQP polish with
  rotation curvature is the next engineering step; (2) record starts are robust to exploration (r126 +13 %, r237 nothing):
  record-beating needs other seeds (neighbours ± k worked partly at 237: 15.9306 → 15.9211 in 30 min, still descending);
  (3) the sub-11 census at 110 lives in one two-channel funnel (all 2 outside entries landed on its rim); a different
  funnel, if it exists, needs a different constructor.
* **SQP polish (`--sqp`): no gain.**  QP subproblem with the Lagrangian Hessian (rotation curvature of active rows, previous
  LP duals, phi tied to owner theta, convexified by diagonal dominance): 30 proposals at 237, median 52 vs 56 iterations,
  same time, slightly worse outliers.  A typical trace converges in ~22 iterations (s settled by ~8); the cost is per LP
  (~0.07 s at 237) × ~2 per iteration + flip-search LPs, not curvature.  Kept off by default.
* **Policy fix instead: `--polish-margin M`.**  At n = 237 the stepping-stone ceiling (k + 0.05 = 16.05) admits every
  proposal (~15.92), so every one paid the full polish.  Now only screened sides < best + M are polished; the rest enter
  the archive at the screened value (flag `unpolished`).  Test: `runs/tx237m` (tx237n settings + M = 3e-3).
* **tx237m (polish margin 3e-3):** 10.3k proposals (3.3× tx237s) but best 15.929641 / 15.929609 / 15.929609 / 15.926733 at
  4 / 10 / 20 / 30 min: no better than tx237n (15.921655) or tx237s (15.924719).  Likely: unpolished entries (sides off by
  up to ~1e-3) flood the archive as parents, so selection and acceptance get noisier; at large n per-proposal quality
  seems to matter more than count.  Single runs; tx237n vs tx237s already differ by 3e-3, so inconclusive.  Left off by
  default.  Next for large n: replicates before more policy changes; cheaper LPs (warm-started active-set / dual simplex
  instead of interior point) are the remaining throughput lever.

## Polish LP timing: Clarabel vs HiGHS, warm starts (10-08; `fq --dump-lp`, `lp_bench.py`)
Question: would a warm-started simplex (HiGHS dual simplex reusing the previous basis) cut the polish's LP cost at n = 237
by >= 2x?  Dumped every polish LP of 8 polishes (screen then `--no-alm --loosen 1.0`; r237 register kicked 0.01 / 0.05,
four tx237n archive states kicked / local-kicked, rec110 kicked 0.01 / 0.05): 1428 LPs, ~2500-2950 rows x ~1350 columns at
237 (~1100-1240 x 600 at 110).  Timed offline one core each (highspy 1.15.1, threads 1).
* **Cold, HiGHS is slower than Clarabel:** per-LP median at 237 Clarabel 32-46 ms, HiGHS dual simplex 44-93 ms (~1900
  simplex iterations), HiGHS IPM 50-114 ms; totals over the 6 polishes at 237: Clarabel 52 s, dual 118 s, IPM 113 s.
  Objectives agree to <= 1e-7 (abs, typical 1e-9).
* **Structure persists:** consecutive LPs keep 100 % of columns and (median) 100 % of rows (min 62 % after a pair-list
  rebuild); only coefficients and right-hand sides change.
* **Second-order-correction LPs (half of all LPs: SOC runs and is accepted on every iteration) warm-start for free:** same
  rows, shifted rhs; 1-2 simplex iterations, ~1-2 ms vs Clarabel ~37 ms; all optimal.  Total at 237: 25.3 s -> 1.2 s.
* **Step LPs do not:** after an accepted step (new point, new coefficients, new rhs) the old basis is primal and dual
  infeasible; warm dual simplex takes median ~1500-2700 iterations (tail 6000), total 67.6 s vs Clarabel 25.5 s; ~15 % of
  mapped bases make HiGHS return kError at once (not diagnosed).  After a rejected step (same point, smaller box) warm is
  0.5-2x Clarabel: no consistent gain.  Primal simplex warm start is worse.  Flip-search LPs: warm ~0.6x Clarabel (small).
* **Hybrid (Clarabel for step LPs, warm simplex / active-set solve for SOC and flips): LP time 52 s -> ~27.5 s (1.9x)**;
  LPs are ~73 % of polish wall time at 237 (71 s total), so ~1.5x on the polish.  Not the 2x+ hoped for.
* **SOC is worth its LP:** with SOC off (experiment build) iterations per polish hit the 300 cap in 4/6 at 237 and end
  1e-8 ... 1.5e-6 higher; per-iteration cost halves, net time similar or worse.
* **Where the iterations go (spotted on the way):** 25-55 % of iterations are rejected (rho < 0.1: 2 LPs wasted each, then
  R/4); R stays 6e-5 ... 1e-3, far below the 4e-3 cap (hence rmax didn't matter); before the SOC rho is hugely negative
  (median -246 on t237c): near stationarity the second-order violation (~kappa R^2 with kappa = 3s ~ 48) swamps the
  first-order prediction.  Long polishes are slow crawls along flat motions (r237big: 8e-8 over iterations 50-300, ~3e-10
  per iteration: just above the 2e-9 / 15-iteration stagnation stop, so it runs to the cap).
* Recommendation: no HiGHS replacement.  If anything, an active-set SOC (reuse the step LP's active rows / bound columns,
  one sparse solve, verify, fall back to Clarabel; `faer`, pure Rust, ~1 day) for ~1.5x on polish; `highs` crates need
  cmake + C++ (cmake not on PATH; nix-shell), ~1 day, same gain.  Probably larger and cheaper: the step-acceptance / stop
  dynamics (rejection rate, R never growing, crawl tail), e.g. a stop on gain per LP relative to s, or a filter /
  nonmonotone acceptance; measure against the benchmark first.

## Explorer benchmark pilot (10-08; `bench_explore.py`, `runs/be_pilot`)
16 frozen tasks (`bench/tasks/`: SQUISH seeds s108 s179 r126 r237 + 12 rd1 pairs, mirror packing fq-polished -> register
side), replicate = explorer 2 procs x T min, paired seeds across variants, score = mean clipped gap fraction over checkpoints.
* Pilot (quick tier 8 tasks, base = adapt + elite 0.3 vs elite0, 3 reps x 8 min, 12 procs, 65 min): elite0 - base = -0.15,
  95 % bootstrap [-0.20, +0.11]: not separable.  Tasks are mostly binary at this budget: rd106/rd130/rd177 saturate (hit in
  most reps), rd210 never moves, rd236 hit-or-nothing (2/3 vs 1/3); graded: s108 (0.23-0.52 spread), s179, rd263.  Next:
  quick tier = graded tasks only, more reps (>= 5).
* rd130 base r0 reached 11.9104141 (from the 09-30 mirror in 7 min on 2 procs): below Couzo's register 11.9111877 but above
  SQUISH's registered 11.904483 (10-06).  No claim; the rd targets are "a known better answer", not the current record.

## Polish-length census (10-08; `polish_trace.py`, `polish_trace_report.py`, `runs/pt_110.jsonl`, `pt_237.jsonl`)
Evan: polish until the contact graph is identified, dedupe on graphs, finish later?  Replayed the explorer pipeline
(move -> screen -> full polish, `fq --ident` with a per-accepted-step load-network fingerprint `fps`) on parents from
`runs/ex10` (110, 500 proposals -> 286 polished) and `runs/tx237n` (237, 150 -> 101).
* **Iterations are not short:** full polish median 20 (110) / 75 (237); p90 114 / 300; cap 300 hit by 3 % / 16 %.
  Polishes >= 80 iterations = 17 % / 48 % of polishes but **59 % / 87 % of polish time**; cap-hitters alone 21 % / 42 %.
* **Network settles early, the side settles with it:** final load network reached and stable from median iteration 8 / 38
  (p90 69 / 188); side error at that point median 2e-10 / 1e-9 (p90 7e-7 / 1.4e-6); **53 % / 32 % of polish time comes
  after settling**.  Summed max-corner displacement after settling median 6e-4 / 8e-4, p90 1.4e-2 / 2.4e-2.
* **Long polishes = crawls, mostly on a fixed network:** of the >= 150-iteration polishes, after iteration 100 most gain
  1e-8 ... 3e-5 total with 1-8 distinct networks over ~150 steps; a minority swap ~1 contact per step (13.8 at most).
* **Exact fingerprints never repeat across proposals:** 3019 / 2903 distinct fps, 0 shared between trajectories; 282 / 100
  distinct final basins from 286 / 101 polishes; 40 / 19 final networks differ from the last step's.  A transition table on
  exact graphs has nothing to match; it needs a coarser key.
* **Duplication vs the archive (by side, 2e-9):** polished proposals ending in an archived basin 28 % (110) / 12 % (237);
  and the explorer's screen-return rule (screened side within 1e-6 of an archive side -> no polish) is wrong 26 % (36/137)
  / 3 of 9 times (a new basin skipped).

## Newton finish: prototype and first Rust attempt (10-08; `polish_replay.py`, `polish_finish_report.py`, `bench/polish110`, `polish237`)
Evan: polish only until the contact graph is identified, finish exactly, no fixed thresholds.  Tier-0 harness: 200 (110) / 80
(237) frozen screened states (explorer proposals from ex10 / tx237n), replayed through polish variants; deterministic
(base vs base: 200/200 identical sides, wall time +-4 %).
* **Prototype = stop the polish at k iterations (`--pit k`), then `exactsolve` (unlifted KKT Newton + all its checks).**
  Default exactsolve (with its f64 jam test first): success 30 / 48 / 72 % at k = 8 / 16 / 64, and **every success equals
  the full polish's side (143/143)**.  Failures: "not jammed" + empty-contact-set crashes (fixed: a rung with no contacts
  now counts as not jammed).
* **`exactsolve --force-kkt`** (new, opt-in: skip the jam test, Newton on the max-support set of the best rung, tol up to
  1e-3): same side 81 / 105 / 145 of 200 at k = 8 / 16 / 64, but also "higher" certified points 10 / 10 / 5 (+3e-6 ...
  +1.9e-2).  The full check set rejects every higher one: corner-corner branch MILP finds a descent (vv_milp_dS < 0) or the
  second-order test finds negative curvature (saddle).  So the finish is self-checking once branch + second-order tests are
  in the acceptance rule.
* **But the savings are small: ideal policy (try at 8/16/32/64, fall back to the full polish) = 0.88x polish time at 110**
  (exactsolve's own Python cost excluded).  The long polishes, where the time is, are not near any stationary point until
  late: at k = 16 the 300-iteration state p0004 has no equilibrium on any contact set up to tol 1e-3 (residual 1.3e-4).
  The 10-08 census's "network settles at iteration 8, 53 % of time after" was contact-level stability, not convergence:
  the point keeps moving along {g_A = 0} with a first-order descent of 1e-6 ... 1e-4 left.
* **What the crawl is (fq Newton-finish debugging, `fq --finish`, opt-in, experimental):** on a fixed network the
  remaining descent is along a curved valley (rotations): a linear step of size R gains ~3e-4 R but costs ~R^2
  curvature error, so rho < 0.1, R shrinks to ~1e-4, ~3e-8 per step: the observed crawl.  Rows close one at a time
  (support 1062 -> 1073 rows over 40 iterations).
* **Rust finish in the lifted formulation does not work:** Gauss-Newton projection onto the dual-support rows converges
  (after a fix: use the post-step vector incl. separators), but equilibrium fails (residual 1e-6 ... 1e-4: the point is not
  stationary, as above); the KKT step with curvature needs multipliers, and in the lifted model the min-norm multipliers
  are meaningless (max 157, 339 / 1057 negative: separator rows are massively redundant), so H_L is garbage; QDLDL also
  needs an augmented-Lagrangian top-left (H + rho J^T J) to stay quasi-definite.  A real finish needs the unlifted
  corner-on-side contacts (exactsolve's form) in Rust.  Not pursued: the prototype bounds its value at ~12 % at 110.

## Cheap polish changes: trust-region rule, precision stop, search test (10-08; `fq --tr-grow/--tr-shrink`, `runs/be_st7`)
* **Crawl anatomy (`FQ_TR=1` per-iteration trace, p0004 at 110):** a 4-step trust-region cycle (R doubles at rho > 0.5,
  the larger step is rejected, R / 4), accepted steps gain ~half their 1e-9 prediction, slope pred/R ~ 1e-5 persists for
  hundreds of iterations: a long, gently sloped curved valley.
* **Replay (`bench/polish110`, `polish237`, vs the default polish):**
  - TR rule grow at rho > 0.75, shrink /2: 0.94x time at both n; at 237 22 lower / 11 higher (up to 5.6e-4): a different
    path, not systematically better.  Dropped.
  - Precision stop `--stag-tol 1e-7` (s gain < 1e-7 over 15 iterations; default 2e-9): 110 0.79x time (19/200 higher, p90
    4.5e-10); **237 0.63x time** (32/80 higher; p90 1e-6, max 3.5e-5; cap-hitters 12 -> 2).  Hard cap 60 (reference): 0.39x,
    max error 8.7e-3.
* **Search test (explorer `--polish-extra=--stag-tol 1e-7 --final-polish 5` vs base; graded tier s108 s179 rd263 rd237,
  5 paired reps x 1 proc x 5 min):** +30-50 % basins per run (s108 median 103 vs 77, rd263 28 vs 18), score difference
  -0.004, 95 % [-0.15, +0.17]: no measurable search gain.  More quenches per CPU do not convert into better packings at this
  budget: the limit is the proposals, not polish throughput.  Benchmark note: scores are quantised by basins (rd237 sits at
  0.43 = one basin in 9/10 reps); the graded tier still needs harder / finer tasks.

## Recombination of same-type basins (10-08, Evan; `cross.py`, `runs/cx1`, `cx2`, `cx3_P*`)
Parents = P nearest (matching distance, best of 8 symmetries, Hungarian) of 30 sampled sub-11 basins of `runs/ex10`, aligned
to a random reference's square IDs; start pose per ID from one parent; explorer screen + polish; novelty vs every side seen
before (certified census + all explorer archives, ~40k sides).  Variants: uniform, bestfit:K (Evan: random ID order, K of P
parents sampled per ID, least overlap with placed squares), field / softfield (smooth random plane-wave field per parent:
hard / soft patches), half (cut-and-splice control), kick control.  Paired trials (same reference + candidates).
* **Not destructive:** 0 grid, 0 discards in 3300 proposals; 87-98 % land below 11; fewer parent returns than kicks
  (26-49 vs 67 of 300).  Half-plane cuts return to a parent most (59/150).
* **Novelty ~ kicks:** distinct new sub-11 per 300 trials: kick 16; P = 5: K2 16, K3 28, K5 26; P = 10: uniform 18, K2 26,
  K3 14, K5 21, K8 22; P = 20: K3 20, K5 30, K8 21 (Poisson +-5).  No K / P trend; best cells ~1.7x kicks at the noise
  edge; deep (< 10.9975) new 0-1 per cell, kicks included.  The 110 funnel is saturated for every move (census 496,
  ~1 new per 10-20 proposals), so it cannot separate moves.  Literature for the K / P choice: gene-pool / multi-parent
  recombination (Muehlenbein & Voigt 1996; Eiben et al. 1994), tournament size (Goldberg & Deb 1991), greedy crossover
  (Grefenstette 1985, TSP), parameterized uniform crossover (Syswerda 1989; Spears & De Jong 1991), cut-and-splice (Deaven
  & Ho 1995).
* Explorer move `--cross W` (bestfit 3 of 5 nearest of 30 sampled archive basins); test on the graded benchmark: `runs/be_cross`.

## Disk -> square shape annealing as a constructor (10-08; `src/anneal.rs`, `anneal_test.py`, `runs/an1`-`an3_*`, `an_ex*`)
Evan: start as a hot disk packing, anneal to squares as it cools (reported to work, also in 3D).  Hard rounded squares
(inner square half-side 0.5 - r, Minkowski radius r; r = 0.5 disk, r = 0 square; exact inner-square distance, unit-tested
vs brute force), NPT compression MC (translate / rotate / ln-area box moves, adaptive steps), schedule over progress t:
bP = bp0 (bp1/bp0)^t (1 -> 3000), r = 0.5 (1 - clamp((t - t0)/(t1 - t0)))^p; shape steps accepted if overlap-free, else the
smallest box growth <= 5 %.  Finish: `fq quench` (loosen 1.02) + grid test.  ~3 s per 20k sweeps at n = 110.
* **n = 110: all grid.**  Every variant (morph early / mid / late, 4k ... 100k sweeps) and both controls (same MC with squares
  throughout; random placement) end at >= 11.0, grid-obstructed; the anneal itself crystallises (~100 axis squares, side
  11.03-11.13 before fq).  110 = 11^2 - 11: the square-lattice crystal is the equilibrium phase and the trivial packing is
  the natural attractor; sub-11 needs a delicate tilted structure the anneal never forms.
* **k^2 + m sizes (37 52 68 87 106 126; 24 runs per variant), gap to the register:** shape annealing beats both controls
  clearly: grid-obstructed 0-4 / 24 (slow, 100k sweeps) vs 17-24 / 24 (squares-only MC, random); best gaps 0.04-0.11 vs
  0.05-0.24 (controls' medians = the trivial grid value).  **n = 37: one slow run reached the register value (+1.4e-12).**
* **Schedule (n = 68 / 106, 16 runs):** slower is better, sub-linearly: median gap 20k sweeps +0.16 / +0.16, 100k +0.11 /
  +0.16, 300k +0.10 / +0.13, 1M +0.07 / +0.09 (47 / 84 s per anneal); best gaps barely move (0.04 / 0.055).  Morph timing
  matters: early morph (square while still fluid, t in 0-0.5) falls back to the grid (12-16 / 16); late (disk until dense,
  0.5-0.95) jams worse (+0.18 median); mid (0.2-0.8) best.  End pressure 3e4 vs 3e3: no gain.
* **As explorer starts (3 best 1M-sweep anneals, explorer adaptive + elite 0.3, 4 procs x 5 min):** 68: 8.840 -> 8.834
  (register 8.799), 106: 10.878 -> 10.870 (register 10.823): the anneal lands in other, worse funnels and local search does
  not bridge them.
* Reading: a real constructor of non-grid starts (the shape morph is what avoids the grid; squares-only MC crystallises),
  but at these budgets 0.04-0.1 above records except at small n (37 hit).  Not competitive at n = 110 or as explorer seeds
  yet.  Untried: much longer anneals / replica exchange in bP at fixed r schedule, r tied to packing fraction instead of
  progress, anneal from a record's neighbourhood (partial melt: r up to ~0.1 locally, then re-cool) as a large-move proposal.
* **Explorer test (`runs/be_cross`, graded tier, base vs base + `--cross 2`, 5 paired reps x 1 proc x 5 min):** -0.038, 95 %
  [-0.067, +0.135]; s108 / s179 worse (s179 basins per run 24 vs 47: the per-proposal alignment, 30 candidates x 8
  symmetries x Hungarian in Python, costs seconds at n = 179), rd263 / rd237 unchanged.  Parked as an available move (cheap
  fix if revisited: cache alignments per archive entry, fewer candidates).

## Benchmark tier v2 (10-08; `bench_explore.py --make-squish-seeds`, `runs/be_newtasks`)
New SQUISH-lineage tasks record(n + 1) - 1 (10-05 register, best of 4 removals), target SQUISH: base 3 reps x 1 proc x 5 min,
gap fraction: s130 0.03-0.45, s154 0.04-0.27, s208 0.11-0.14, s238 0.09-0.10 (graded); s88, s263, s303 0 (stuck at this
budget).  `graded` tier = s108 s179 s130 s154 s208 s238 (old 4-task tier kept as `graded1`).

## Size chaining (10-08; `chain.py`, `runs/ch2`; v1 `runs/ch1` discarded: it dropped seed-lineage progress every round)
Window 123-131 from the frozen 10-05 register (pre-SQUISH; seeds for the top from 132, 133), rounds = explorer 2 procs x 4 min
from best(n) + carried lineage basins (top 3) + fresh seeds best(n+1)-1, best(n+2)-2; frontier share 0.7 on the seed
lineages; 8 slots, 1 h (16 cores), 112 rounds, 11-14 per n.
* **126: 11.7747351 -> 11.7735853 (-1.15e-3, 80 % of SQUISH's gain; SQUISH 11.7733036 registered 10-06: no claim).**
  From the record's own lineage (bigkick in round 87, rowslide in 103), not the seeds: repeated fresh 4-min explorations of
  the carried best move a record that one run (bench r126: 13 %) does not.
* Seed lineages descend steadily but stay above the records: 129 11.9085 -> 11.8823 (record 11.8813), 130 11.9505 ->
  11.9140 (record 11.9112); others slow.  Many removal seeds are copies (n+1 record minus a square, polished back to the
  same side: 126's lineage sits at record(127) all hour); 133 = 12 x 12 grid, useless as a seed source.
* Next: reject non-shrinking removals / more removals per seed; live register + a window SQUISH has not covered.

## Regional melt move, seed fix, live register (10-08; `anneal melt`, explorer `--melt`, `runs/be_melt`)
* **Melt** (explorer kind `melt`: squares within R (1.5-3.5) of a random square rounded up to rmax 0.05-0.2, hard-particle
  MC at bP 2000 for 3000 sweeps around them, squared again; ~3.3 s per proposal at n 108-238): graded tier v2, 5 paired
  reps x 1 proc x 5 min, base vs base + `--melt 2`: +0.055, 95 % [-0.03, +0.12]; **s130 0.03 -> 0.55 (all 5 reps >= 0.43),
  s154 0.04 -> 0.18, s108 0.15 -> 0.27; s179 0.17 -> 0.00 (basins 26 vs 47)**; s208 / s238 unchanged.  Size-dependent:
  kept in the move set with the bandit deciding (adaptive, per-size statistics carried across chain rounds:
  `--bandit-state`, discount 0.5).  165 melt proposals: 57 new archive entries, 47 screen-grid.
* **Seeds v2** (chain): removals that do not shrink are rejected (copies of the n + k record), up to 12 tries; n = 123-131:
  copies 15/36 -> 0/36, no usable seed 0 -> 13/36 (rigid neighbours), best seed gap 2.2e-3 -> 2.1e-4.
* **Live register** (jlevy/squares main tarball 10-08 -> `../exact/batch/regnow/`, `inputs_live/`, both gitignored; table
  `runs/register_live_2026-10-08.json`): 247 open (102 trivial-grid bounds), 54 changed since 10-05 (Couzo 31, SQUISH 23);
  our 266 / 270 / 272 still below it and unregistered.

## Record hunt hunt1 (10-08; `chain.py` live mode, 112 open non-trivial n in 100-324, live register, `runs/hunt1`)
* **s(132) ≤ 11.987099332245063227179742877435** (cert S'; KKT 11.98709933224506322705987…), register 11.991327887691501
  (Couzo 09-27; re-fetched raw n-132.md 10-08, no pending PR / issue with a 132 packing): **−4.23e-3**.  Strict local min
  modulo 33 flat motions, jammed in every branch (92 corner-corner touches); certificate VALID (exactsolve exact check +
  verify_cert2).  Roles (L, B, axis, other) (18, 23, 91, 0).  First round at 132 (2 procs x 4 min), from the record itself:
  bigkick -> aswap -> kick -> aswap (12.0368, above k) -> crot -> chainshift -> kick -> aswap -> kick -> aswap, 227 s: an
  uphill excursion through stepping stones above 12 into another funnel.  `candidates/hunt1_n132.*`.
* **s(155) ≤ 12.95249894401400738196585056694** (cert S'), register 12.9525032026045129 (SQUISH 10-06; raw n-155.md
  re-fetched 10-08): **−4.26e-6**.  Strict local min modulo 2 flat motions, jammed in every branch (146 corner-corner
  touches); VALID under both verifiers.  `candidates/hunt1_n155.*`.
* **hunt1 final (3.5 h, 16 cores; 241 rounds; + hunt1b n = 87-90, 4 cores, 100 rounds):** records 132 (−4.23e-3) and 155
  (−4.26e-6), nothing else below the live register.  Closest seed lineages: 299 (+7.6e-9), 204 (+1.6e-7).  87-90 none;
  90 non-grid lineage 10.0095668 -> 10.0090596 (> 10).  CPU by kind (bandit): bigkick 25 %, aswap 18 %, chainshift 8 %,
  kicksym / kick / lkick / crot / band / mirror / rowslide 5-7 %, melt 4.5 %, reinsert 4 %.
* Scheduling flaw: after the first pass, the lineage-gain credit concentrated follow-ups (rounds per n: 94 n x 1, a few
  n x 11-17).  Next: a full second pass before any third, then priority.

## s(110) census re-check: slack-separated pairs (10-08)
* Found while explaining the census: two "distinct" certified minima, 10.99756849434 (`runs/ex0/b00567`) and
  10.99756849478 (`runs/ex1/b00454`), are 7.6e-7 apart in total corner movement.  Same touching pairs; at squares
  104/105 (37/38 in b00454; 0.011° apart, nearly parallel) the higher one rests a corner of 105 on 104's side, while
  105's side also clears every corner of 104 by 1.6e-8.  Non-overlap is a disjunction over separating lines, so that
  incidence is not a constraint there.  exactsolve fixed one line per pair and certified a branch artifact.
* Fix (`search/exact/exactsolve.py`, `slack_separated`; README "Limitations"): such pairs are released before the
  multiplier, second-order and MILP analyses.  The verdict is now also NOT a local minimum when the equilibrium
  residual is > 1e-8 or the corner-corner MILP finds a descent.  Before, `census_exact.classify` read only
  `second_order.status`, so two census entries with MILP dS < 0 were counted.  Regression input `inputs/s110_slackpair`.
* Re-certified every census minimum (493 with files, plus the 7 older jx1/rev1 ones; ~4 min on 15 procs):
  **3 not minima** (10.99756849478 → slp2 lands on …434; 10.99884293949 → 10.99884139849, known; 10.99957319574,
  MILP descent), 1 a minimum after release (10.99986693230, 2 pairs released), the rest unchanged.
  `runs/known110_all.json` 496 → **493** (9-digit sides; **497 distinct exact minima**, since 4 pairs of genuinely
  different packings agree to 9 digits: corner distance 0.017–13).  Run caches (`runs/*/exact.json`) reclassified.
* Exact register batch (`search/exact/minpoly/solve`, 323 n): same scan, only n = 150 has slack pairs (4, clearance
  2–3e-4); re-solved, still a local minimum (released pairs carry no load).  Certificates (upper bounds) are unaffected
  throughout.  Lean `IsLocalMinPacking` proofs (N11L, N28L) are sound: `RowOK` requires the other three side values at
  each contact point ≤ 1/2 − μ (`lean/Sqpack/LocalMinRows.lean`), which excludes this case.

## Targeted hunt h2 (10-09; 88, 89, 90, 109, 147, 148, 290, 291; `runs/h2r1`, `runs/seeds1009`)
* **Register refresh** (`../exact/batch/regsync.py PACKING_DIR`: tarball → `regnow/` (old kept as `regnow_<date>`),
  `inputs_live/`, `runs/register_live_<date>.json`, diff): 6 changed since 10-08 (68, 105, 292 Couzo/refinements; our
  266 / 270 / 272 registered); none at the targets.  Pending (open issues/PRs): ry-xu #432 s(146) ≤ 12.583782 (register
  12.600908), Mishapolk #470 (84, 86, 105, 108, ...), Gupta #438 refinements; nothing pending at a target n.
* **Seed pool** (`import_ext.py scan/best`: Ellsworth text, Couzo text, Witness/v2 yaml, ours; f64 SAT check): register,
  register git history (jlevy `known-best` per-n commits), Ellsworth mirror incl. `_r*` revisions, ry-xu, Mishapolk,
  Couzo, SQUISH, our hunt archives: ~15k packings over n in 84-93, 104-113, 142-152, 284-296 (`runs/seeds1009/pool`,
  per-n best in `best/`).  Structure of the targets: 89, 109, 147, 290 are pure 45-degree constructions (89 = 7 x 7
  diamond, 49 tilted); the new 88 (SQUISH) and 146 (ry-xu) use ~31-37 squares at 41-50 degrees.
* **chain.py**: `--seed-ks` (k < 0 = insertion into best(n + k) at pose holes after a small dilation), `--seed-tries`,
  multiple `--carry` per n (was: last one won).  Insertion seeds mostly quench to the grid (10.0, 11.0, 13.0): inserting
  into a jammed record pushes a whole row out.  Removal from the regular 89 / 148 / 291 records = copies (rejected).
* **Gadget transplant** (`transplant.py`: tilted set of a source at n + 1, exact axis fill, scan s): 146 → 147, 88 → 89,
  108 → 109 never gain the extra square below the target (the gadget itself has to change).
* **Strips** (`strips.py`: K parallel strips at theta, offsets, per-strip shifts, exact fill; `--cut` = no unit band free
  of tilted squares, `--deficit` = add missing squares at holes): at 110 the best fill is 104-105 at s = 10.999; quenches
  land at 11.0 (grid) or 11.00-11.03; **explorer from 19 such stepping stones (2 procs, 10 min): 2662 / 2900 proposals
  screen-grid, nothing below 11 (Evan's test: strips do not reach the 110 funnel).**  Stopped early.
* **Blocks** (`blockgen.py`: a x b lattice rectangle at theta, optional staircase rows (length, shift), placement as
  fractions of the free range or centre + quarter-lattice offset `c+dx,dy`, exact fill, fq): reproduces 89 (7 x 7, 45,
  9.949747468), 65 (5 x 5), 66 (8 x 3 at c + (1/4, -1/4)), 150 (11 x 6 at c + (-1/4, 1/4)) to 1e-9.  Not 69 (two angles,
  Cantrell-type) or ry-xu 126 (staircase diamond: rows of different lengths and shifts).
* **Round 1 (h2r1, 8 targets x 2 procs, 100 rounds, 1 h 25 min): nothing.**  No target moved; lineages: 290 17.5809 →
  17.5628, 291 17.5906 → 17.5869, others unchanged.  Stopped early (Evan: retarget compute to things with a reasonable
  chance; add features that open new targets).
* **Block sweep** (`runs/bgsw`, 600 + 600 evals per n, 45-degree rectangles with quarter-lattice centres, then 38-45 degree
  staircases): just below each record side the fill is 1-4 squares short at 147, 148, 89, 109, 90 (89 of 90 at
  9.9999), 290, 291; near-45 angles always lose more than they gain to the exact fill.
* **Staircase positive control (126):** seeds built from Xu's measured column shifts (+ fill + 1 square at a hole,
  quenched: 11.900, 11.923) → explorer 2 procs: **11.7426406871 = 15/2 + 3 sqrt 2 (Xu) in < 4 min**.  So arrangement-level
  seeds well above the record work.  Blind (`stairgen.py`, CEM over placement and column offsets): fills top out at
  123/126 at 11.76; the exact fill cannot express Xu (with Xu's exact tilted set and only the 36 45-degree squares it
  gives 125; the other axis squares sit in pockets against the block and against near-axis rotated squares; residue-class
  candidates: 171k, still 61/90 covered).  stairgen now scores by the fq-quenched side of a complete n-square packing
  (Evan: always n squares, score by side; deficit-insertion was a hack).  Bug fixed on the way: numpy-2 `repr` of
  np.float64 wrote `np.float64(...)` into seed files.
* **Sub-grid library** (`chain.py --lib-init J --lib-k K`; Evan 10-09: sub-grid packings are good seeds because they are by
  definition not grid-jammed): per n, every distinct (1e-9) sub-k basin the explorer finds is kept (`<out>/lib/n<n>/`,
  300 best); each round starts K random library entries as frontier lineages.  Pre-loaded for hunt2 from hunt1 / hunt1b /
  h2r1 round archives + seed pools: ~63k files over 118 n.
* **hunt2** (`runs/hunt2`): second pass over the 118 open non-trivial n in 84-324 (`--first-pass`), starts = best of
  register + pending (ry-xu, Mishapolk, Couzo, SQUISH) + ours (22 n start below the register), removal seeds k = 1, 2,
  library K = 3, 6 slots x 2 procs x 4 min, 2.5 h.
* **110 from-scratch test** (`runs/t110_ex1`; Evan: can we find sub-11 packings without any sub-11 seed?): 40 line-free
  seeds 11.0000067-11.05 from the 10-04 layout `hyb` / `ftmc` runs (hand-made two-channel genomes, no sub-11 input) →
  explorer 4 procs x 30 min, novelty vs the 493-minimum census.
* **hunt2 n = 270 "candidate" = refinement, not a new basin.**  Start = Couzo's pending 270 (16.936723155037; jlevy#469
  T-130, also in Mishapolk's repo): exactsolve says *not jammed* (residual 1.6e-5); `fq quench --loosen 1.0` on his file
  alone gives 16.9367230228797, = the chain's best (lineage: the start entry).  Exact: **s(270) ≤
  16.936723022876183407…** (cert S', VALID under verify_cert + verify_cert2; strict local min modulo 25 flat motions, MILP
  dS 0), 1.32e-7 below his reported side.  Construction Couzo's; ours is the exact optimum (`candidates/hunt2_n270.*`).
  Lesson: polish every imported start once before searching (pending packings may be unconverged); chain counts a
  polished start as a candidate.
* **Start-polish survey** (every start in `seeds1009/best2` + `best3`, 288 packings, `fq quench --loosen 1.0 --stag-tol
  1e-12`): everything registered or pending in 84-324 (SQUISH, ry-xu, Gupta, Mishapolk, register) is converged except
  Couzo's 270 (above).  n > 324: Couzo's 375 19.907024692 → 19.907024297, but exactsolve: not jammed, corner-corner MILP
  descent −4.9e-6; `--no-alm --flip-top 128 --pit 600` → **s(375) ≤ 19.907024202136836…** (cert S', VALID both
  verifiers, strict local min modulo 67 flat motions, MILP dS 0; `candidates/polish_n375.*`), 4.9e-7 below Couzo's file
  (`2d32a6e`).  Not in the register (n > 324).  So fq's default flip search (top 8) stops early at large n.
  (Mishapolk 343-360 drop by up to 0.2 but stay at / above the 19-grid.)
* **s(110) census continuation** (`runs/c110_h2`: record + Ellsworth + 20 random census minima, 4 procs x 85 min, elite
  0.3, global novelty): 283 sub-11 archive entries, 249 "new" at screening; exact: 203 certified (202 distinct), 77 not
  minima, 3 unresolved; **13 new certified minima, all 4.2e-4 … 3.2e-3 above the record** → `known110_all.json` 493 →
  **506** (old list kept as `known110_all.2026-10-08.json`).  The three screening "new" basins within 1e-4 of the record
  (10.9968326, .9968339, .9968812) were *not jammed* (exactsolve) and fq `--no-alm --flip-top 128` did not move them;
  `--kick 1e-4` + quench (6 seeds) sends every such state into a known minimum (record x3, Ellsworth's, 10.9968870 x4,
  10.9971899).  **Polish stall:** fq can stop at non-jammed states; the screening "new" count overstates (here 249 vs 13).
* **big1 (n 327-379, 22 n with sub-grid packings: Couzo `2d32a6e` + mirror; 2 slots x 2 procs x 6 min):
  s(375) ≤ 19.906981767290601595…** (cert S', VALID under verify_cert + verify_cert2; strict local min modulo 64 flat
  motions, MILP dS 0; `candidates/big1_n375.*`): **4.29e-5 below Couzo's 375** (19.907024692, his repo `2d32a6e`,
  unchanged since 10-08) and 4.24e-5 below our exact polish of it.  One `kicksym` from the polished start (round 25),
  i.e. a different basin next to his.  n > 324 is outside jlevy's register; no other 375 packing found.
* **c110_h3** (8 procs x 60 min from 13 near-record screening states + records + 10 census minima): 314 sub-11 entries,
  225 certified distinct, 85 not minima; **17 new certified, all ≥ 6.3e-4 above the record** → census **523**.  Both
  runs together (~12 CPU-h): 30 new certified minima, closest 4.2e-4 above the record, none within 1e-4 (consistent
  with the ~5 % / 40 CPU-h estimate in NEEDS.md; the funnel keeps producing minima away from its floor).
* **big1 final** (43 rounds, 2.9 h x 4 procs): only 375.  **big2** (from big1's bests, 4 slots x 3 procs x 8 min, 38
  rounds): nothing; closest seed lineages 335 (+3e-5), 372 (+1.4e-4), 332 (+5.8e-4): start there next time.
* **hunt2 final** (123 rounds = one pass over 118 n + a few): only the 270 refinement; lineages tie the record at 176
  and 299, 263 within 1e-6.  A second pass over 84-324 finds much less than hunt1's first (132, 155).

### Session summary (10-09)
Results: s(375) ≤ 19.906981767290601… (new basin, 4.29e-5 below Couzo's); exact refinements of Couzo's 270 (1.3e-7) and
375 (superseded); s(110) census 493 → 523 (30 new certified minima, none within 4.2e-4 of the record).  Draft post list:
`outreach/draft-results-2026-10-09.md` (not posted).  What the program needs: `NEEDS.md`.

## Frontier refresh (10-09 evening; jlevy main 657cc48 = #478 merge of our #465)
* Register sync: 37 n changed since 10-08 (mostly ry-xu #432 and SQUISH/SidG2k1 refinements imported; our 266/270/272 in).
* Pending (open issues): itsnaka **#481** 15 new packings 131-307 at 1.7e-4 … 1.1e-2 below the register (131, 153, 207,
  209, 232, 236, 259, 263, 269, 270, 292, 302, 303, 305, 307; + 154, 237), certs itsnaka/squish-certs `d45669b`; Couzo
  **#476** (132, 237, 263, 267, 270, 303; exact optima via our exactsolve, certs `02f9690`); Mishapolk #470 (28 n, 84-306,
  certs `dd3da5c`); Couzo #451; SidG2k1 #438; ry-xu #432.
* Ours: **132 beaten** (Mishapolk 11.986956, Couzo 11.986954 from it, vs ours 11.987099); **270 beaten** (itsnaka
  16.929780, new structure); 155 ties Couzo #451 (same basin); 266, 272 still best.  **375**: Couzo's repo `02f9690` now
  has 19.906981767294 = our big1 basin (ours 4e-12 lower, exact optimum): nothing to post.  10-09 draft superseded.
* Implication (NEEDS.md policy): ~25 fresh single-source sub-grid packings (#481, #470, #476) = the best perturbation targets now.

## Tooling session (10-09 evening, branch `tooling`; DESIGN.md, pk/README.md)
* `pk/` toolkit: store (`runs/store.sqlite`: register, pending issues, 1,001 run archives, seed pools; 169k packings),
  `fq serve`, move registry with explicit parameters, pipeline, replay bench (`pk/corpus.py`).
* **Bench c110a** (1,760 frozen proposals, 40 certified parents x 11 arms x 4): polish-everything finds 30 new certified
  basins (66 / CPU-h) vs 16 (48 / CPU-h) for the 10-08 recipe; its screen-return test (1e-6) misassigns 292 / 710 returns,
  24 of which polish to new certified basins.  Staged polish worse.  Per arm: aswap 177 / CPU-h, kick 157, kicksym 137,
  lkick 134, rowslide 131, bigkick 101, mirror 52, crot 49, band 48, chainshift 0, reinsert 0.  -> explore.py: `--same`
  default 0, chainshift / reinsert weight 0.
* `anneal sched` (curves for pressure, corner radius, rotation share; region; snapshots) + `pk/crossing.py`.
* **Crossing x110a** (n = 110; 1,120 trials + controls; 9.5 min x 16): **0 crossings in every cell** (95 % upper bound
  8.8 % per 40-trial cell).  Controls: pos (sub-11 packing dilated to 11.01 + kick 0.01) crosses back under quench 10/10,
  melt 9/10, partial global anneal (rmax 0.2) 4/10, full global anneal (rmax 0.5) 0/10; neg (n = 117) never below 11.
  Global anneals send above-grid starts to the grid 38-40/40 (local quench / melt / small kick keep them non-grid, at
  11.015-11.02).  So the hard-square NPT anneal is a grid-maker: it destroys sub-grid structure it is given.  Lock-in
  (quenched snapshots) already by t = 0.9; from scratch the line occupancy reaches k at t ~ 0.8-0.85 (corner radius ~0.1).
  Caveats: from-start anneals begin at bp1/10 (no expansion phase), so grid starts never melt; the occupancy proxy counts
  squares within 1 degree of axis (near-axis squares in sub-11 records trip it).
* **Rate vs depth (rate110, depth110; 30 trials per cell):** anneal length 30-30,000 sweeps changes nothing; melt depth
  (corner-radius amplitude) decides.  Positive control (sub-11 pushed to 11.01) crosses back at rmax 0.03-0.1: 28-29/30
  (record 10.9967834 reached several times), 0.15: 23/30, 0.2: 16-17/30 (same at 30 and 100 sweeps), 0.3: 4/30, 0.5: 0/30.
  Above-grid starts end in the grid 6/30 (0.03), 8 (0.1), 14 (0.15), 16-19 (0.2), 27 (0.3), 30 (0.5), at every length.
  From scratch: grid at every depth and length.  So "too much annealing" = too deep, not too slow: once shapes are round
  enough (r ~0.2+) the orientational structure is lost and squaring up re-forms rows.  Shallow global anneals (rmax <= 0.1)
  are a candidate explorer move (bench them).
* statarb review (research agent, `~/predictions/statarb`): their two-basin problem had the same signature (nested branches:
  heating connects continuously to the high-entropy basin; no temperature rung helps, `tasks/level1_exact_flip_oracle/
  findings_20260710.md`); fixed by block independence proposals with Metropolis correction, not by tempering.  Suggested
  here: measure F(q; P) by umbrella sampling in q = full-line count / line occupancy; bias ladder in q; wholesale region
  re-packing with sub-grid motifs + quench.  Code pointers: `scripts/respace_ladder.py`, `tasks/level1_exact_flip_oracle/
  barrier_profile.py`, `statarb/utils/excursions.py` (N_eff = N / (1 + CV^2)).
* **Q4 order parameter** (Q4 = mean cos 4 theta; the 110 record and all certified sub-11 minima: 0.472).  `anneal sched
  --bias lin|harm --bias-lam CURVE --bias-q0 CURVE`, `--log-every`; `pk/umbrella.py` (windows + WHAM).  u110a (P = 30,
  100, 300, 1000; 20 windows x 2 seeds): P = 30 is fluid (window sd 0.02-0.03), F rises ~1 kT per 0.01 of Q4, F(0.47) -
  F(1) ~ 55 kT, and the aligned states are also denser (<s> 11.38 vs 11.5-11.8): no entropy-vs-density trade in the
  fluid, the grid wins both.  P >= 100 is jammed (sd 0.003-0.006, seeds disagree): profiles non-equilibrium.
* **Biased-anneal crossing** (bias110, lag110; 30 trials per cell): 0 crossings.  Artifact found: from-start anneals that
  melt deep end with a forced uniform box growth of 19-30 % (shapes cannot square up in place at high pressure), and the
  quench from s ~ 13 lands on the grid.  `--lag-tol` (expand while shapes lag) removes it; grid share at full melt 30 ->
  25 / 30, and the Q4 harmonic bias at rmax 0.3 lowers it further (21 -> 6 / 30 above-grid), but the results are tilted
  glass at median 11.02, not sub-11.  The lag expansion hurts the positive control (rmax 0.15: 7 / 30 back below 11 vs
  23 / 30 with the final uniform growth, which preserves the arrangement), so it is off by default.
  Reading: Q4 can keep packings off the grid but the sub-11 funnel is a needle inside the Q4 = 0.47 slice; Q4 is not the
  order parameter that separates the funnel.
* **Shallow-anneal arm + reduced-start rediscovery** (`ashallow` = anneal sched diagonal, rmax 0.03-0.15, 300-2000
  sweeps; `pk/rediscovery.py`).  c110b (40 certified parents x 11 arms x 8, 3,520 proposals, 4 min): new certified basins
  per arm: kicksym 11, ashallow 6, rowslide 6, melt 5, aswap 5, kick 4, lkick 4, mirror 4, band 3, bigkick 1, crot 0.
  Per-arm counts are Poisson-small; ranks moved between c110a and c110b (aswap 9 -> 5, kicksym 6 -> 11), so only
  differences >= 2-3x are real: crot, bigkick, band weak; ashallow in the top group (59 / CPU-h).
  c110r (record + 31 basins within matching distance 0.02 withheld; 28 parents at 0.02-0.05, 0.05-0.1, >= 0.1; 2,464
  proposals): record hit once (rowslide), 9 withheld basins rediscovered, **all from parents at 0.02-0.05; none from
  0.05-0.1 (1,144 proposals) or beyond**.  Single-step reach ~0.05 for every arm incl. the shallow anneal (0 withheld).
  Reaching further needs lineages (multi-step), the next test.
* **Value-graded bandit reward** (`explore.py --reward value --reward-g0 1e-3 --reward-alpha 1`; Evan 10-09: novelty is
  novelty, grade it by distance from the frontier, no above-grid special case).  reward = novelty x min(1, (g0 / gap)^alpha)
  x parent robustness (posterior mean of its children's non-grid rate).  Calibration from lineages (c110_h2/h3, hunt1,
  hunt2; expanded states by gap at creation): P(subtree finds a new basin within 1e-3 of the start's best) 0.4-0.9 at gap
  < 1e-3, 0.01-0.07 at 1e-3..1e-2, <= 0.007 at 1e-2..3e-2, ~0 beyond; near-best finds per 100 descendants 2-35, 0.1-2.3,
  <= 0.25, ~0.  About one decade of value per decade of gap.  Learned under a policy that rarely chains far states
  (subtrees 5-15 descendants): likely understates far value; alpha is the knob, the lineage test measures it.
  `--reward legacy` (default) reproduces the old reward; `--reward-above` removed.
* **Idea battery** (`pk/battery.py`: frozen suite at 110 = 40 bench parents, 30 reach parents with the record's 0.02
  neighbourhood (31 basins) withheld, 60 clean above-grid starts, 20 positive controls; Gamma / Beta posteriors;
  P(candidate > baseline)).  Baseline = explorer arm mix, 15 min (7,920 proposals); candidates 10 min each.
  | metric | baseline | ashallow | cross2 | crossP5near | crossP5far |
  | new certified / CPU-h (bench) | 27 | 34 | 35 | 19 | 22 |
  | withheld rediscovered (distinct) | 2 | 3 | 11 | 7 | 9 |
  | children moving toward the record | 0.07 | 0.16 | 0.34 | 0.29 | 0.41 |
  | above-grid starts -> grid | 0.42 | 0.23 | 0.71 | 0.33 | 0.86 |
  | positive control returns below 11 | 0.62 | 0.92 | 0.82 | 0.93 | 0.81 |
  | CPU s / proposal | 1.9 | 1.7 | 1.2 | 1.5 | 1.4 |
  Reach caveat: every withheld hit, for every arm, involved a participant (parent or mate) within 0.05 of the record (the
  six nearest reach parents, 0.040-0.049).  Per proposal involving such a participant: baseline 0.002, ashallow 0.039,
  cross2 0.074, crossP5near 0.028, crossP5far 0.022 (hits / proposals).  So recombination does not extend reach beyond
  its nearest participant; two-parent crossover is the most efficient use of near pool members (~2x the shallow anneal,
  ~35x the arm mix); more parents dilute.  Bridging (`pk/bridging.py`): two-parent children with mates 0.05-0.1 apart land
  ~midway (median 0.047 from both); five-near children collapse onto the nearest mate.
  One long jump into the funnel: crossP5near child of an above-grid start (11.0085) + 4 above-grid mates -> certified
  10.998018 (known basin), 0.108 from its parent; but the parent is 0.069 from the nearest sub-11 basin (the suite's clean
  filter compared against the 30 best sub-11 references only: tighten to all basins at the next freeze).
* **s(110) census 523 -> 683** (`runs/known110_all.json`; old list `known110_all.2026-10-09.json`, gitignored copy):
  160 new certified sub-11 minima from the 10-09 evening bench / battery / crossing runs (every proposal fully polished,
  every distinct sub-11 side certified), deduplicated at 5e-9.  Closest to the record: 10.997517877 (+7.3e-4), then
  1.1-1.3e-3 above.  None within 1e-4.

## Night session 10-10 (lineage reach; explorer port; Hunt 3)
* **Explorer port** (`39dde93`, `e79fa6e`): `cross2` (two-parent crossover, variant drawn from field / softfield / half /
  bestfit2; random sub-k mate for sub-k parents, nearest of 8 for above-k parents) and `ashallow` (shallow whole-box anneal),
  weight 2 each by default (`--cross2 0 --ashallow 0` = old mix); mates recorded in the archive (`mate`) and proposal log.
  Battery split by crossover variant: no variant resolvable (1-11 new certified per cell); `bestfit2` lowest bench yield in all
  three crossover arms (4 / 2 / 2), suggestive only.  The mate rule, not the variant, carries the effect.
* **Per-move landing spread** (`move_spread.py`, battery trials, bench + reach suites; d = child -> parent matching distance of
  polished sub-k non-grid children): single moves land close (median 0.003-0.009, q90 0.012-0.044, <= 6 % beyond 0.05) and
  23-69 % return to the parent's own basin; big pushes are pulled back or lost (bigkick: 60 % lost, survivors median 0.009;
  band / crot 60 %, mirror 78 % lost).  Crossover is the only move that lands far while staying sub-grid: cross2 median 0.041,
  34 % beyond 0.05, 14 % lost; crossP5far median 0.051.  So "single-move reach 0.05" = the tail of small steps.
* **Rattler check** (exactsolve free squares, distance recomputed on pairs where neither square is free): for most packings free
  pairs carry <= 7 % of d^2, but some starts are inflated: b70455 0.056 -> **0.046** (38 % of d^2 from free squares), b59218
  0.075 -> 0.068, b61927 0.071 -> 0.068; group candidates b54624 0.055 -> 0.044, b59788 0.073 -> 0.065, b74548 0.072 -> 0.062.
  Band edges near 0.05 need the masked distance.  Same basin re-found reads d ~ 0.001-0.002 (noise floor).
* **L1 lineage test** (`lineage_reach.py`, `runs/lineage/L1`; 12 certified starts at d 0.050-0.083 from the s(110) record,
  closer basins unknown to the run (novelty ref = basins >= 0.05); new vs old mix, 4 procs x 60 min each, concurrent; hunt
  recipe `--adapt --elite-share 0.3 --reward value`): both reach the record.  new: < 0.02 at 3.5 min, record 4.0 min (from
  b70455, depth 7, chain starts with cross2), 108 basins < 0.05 / 26 < 0.02.  old: < 0.02 at 6.9 min, record **24.4 min from
  b5457 (masked d 0.063, the d ~ 0.062 family 1.4e-4 above the record in side), depth 14, local moves only** (kicks, bigkick,
  band, rowslide, aswap): a lineage crossing that cannot have borrowed.  new also reached the record's neighbourhood from b5457
  (11.2 min, depth 14, 4 cross2 steps; mates unrecorded in L1).  Caveats: b70455 is inside single-move reach once rattlers are
  masked; the elite share / value reward concentrate parents on the lowest-side starts (b5457 175-235 sub-11 descendants;
  most far-in-side starts 1-56), so far starts were barely tested; one replicate, wall-clock under ~22 busy processes.
* Census side note: the 15 lowest basins >= 0.05 from the record all sit at d 0.062-0.064 (gap 1.0e-4 ... 4.2e-4): a second
  structural family nearly as good as the record (lowest: 40601 = 10.9968870, the often-revisited census minimum; Ellsworth's
  10.99679327 is not in it); the 0.02-0.05 ring has nothing below 1.1e-3.
* **L2 lineage test** (`runs/lineage/L2`; one explorer per start group x recipe, 2 procs x 60 min each, 6 concurrent, so no
  cross-group starving or borrowing; 4 starts per group, all vetted by rattler-masked d; mates recorded):
  | group (masked d of starts; gap above record) | new: t(< 0.02) / t(record) | old: t(< 0.02) / t(record) |
  | far2 (0.071-0.082; 1.8-3.2e-3) | 3.2 / 4.4 min | 1.5 / 1.5 min |
  | far1 (0.051-0.065; 2.6-3.1e-3) | 6.0 / 6.7 min | 0.9 / 1.8 min |
  | near (0.063-0.064, the second family; 1.4-4.2e-4) | 4.6 / 24.6 min | 14.8 / 54.3 min |
  Every run reaches the record.  From the far-in-side starts it takes ~7 steps and minutes even with local moves only
  (far2_old0: b62620, d 0.071, 7 moves, 1.5 min).  Crossing from the second family is slower (25-54 min, depth 7-16).
  New vs old: old faster on far groups, new faster on the near group; one replicate, so no recipe verdict.  Borrowed mates
  appear in new-recipe chains (far1_new0 record chain: 2 of 2 cross2 steps borrowed from another start's lineage).
* **Reading (L1 + L2):** single-move reach (~0.05) is not the limit for lineages: within the s(110) sub-11 region chains
  cover 0.08 in minutes.  Nearly every certified sub-11 basin is within 0.1 of the record (2 of 537 beyond), and lineages
  from all three groups find it, so the whole sub-11 region at 110 behaves as one funnel with the record at its floor; the
  second family (d ~ 0.062, 1e-4 above) is a side valley, slower to leave.  L1's far starts failed only because the elite
  share / value reward starved them (1-56 descendants).  Consequences: (a) the hard step at 110 is entering the sub-grid
  region (0 grid crossings), not moving within it; (b) parent allocation matters more than the move mix once below k: a
  picker that starves high-side starts loses whole lineages; per-start (or per-cluster) budgets fix it; (c) the reach
  question for hunting is whether other n have several sub-k funnels; if so, per-start explorers are the way to sample them.
* **Hunt 3** (`runs/hunt3`; 22 n whose best known packing is pending and unperturbed by us: itsnaka #481, Couzo 132 / 267,
  TheSnakeFang #484 (308 / 343 / 344, new 10-10); `chain.py --first-pass`, 8 slots x 2 procs x 6 min, new recipe +
  `--reward value`; started 02:06).  First hour (41 rounds), certified (exactsolve exact optimum + verify_cert + verify_cert2):
  - **s(132) ≤ 11.98568019884580881140…** (strict local min mod 31 flat motions, MILP dS 0): **−1.27e-3** below Couzo's
    pending 11.986954193641.  New basin: depth 8 from a library seed (11.98736) via bigkick, rowslide, lkick, kicksym x2,
    cross2, kick.
  - **s(305) ≤ 17.95118356772045354314…** (mod 47): −1.26e-5 below itsnaka's 17.951196144147 (converged); one rowslide.
  - **s(308) ≤ 17.99846539275931059669…** (mod 52): −8.44e-4 below Kevin Fang's 17.999309855163 (register: the 18-grid).
    His file polished in place (fq --no-alm): his construction was unconverged; ours is its exact optimum.
  - **s(344) ≤ 18.99445051427740152092…** (mod 71): −1.04e-3 below Kevin Fang's 18.995489275431.  His file polished in
    place gives 18.9944595 (unconverged, like 308); that was not jammed (MILP dS −3.1e-4); fq `--flip-top 128 --pit 600` did not move it; `exact/descend.py` (MILP kick + slp2) did.
  - **s(343) ≤ 18.989552952452492** (dilation certificate, new `dilate_cert.py`: centres and box x (1 + 1e-12), rational
    rounding; VALID both checkers): −5.35e-3 below Kevin Fang's 18.994903529220, but not a local minimum (MILP dS −2.0e-4;
    descend.py's kick fails at every size tried).  Lineage: removal seed from best(344) (one square out) -> cross2, crot,
    ashallow; his 343 polished in place gives only 18.99117.  Round 43 already has 18.98650673 (re-certify at the end).
  Morning (Hunt 3 ended 07:43; dedicated explorers x343 4 procs 120 min, x343b 3 procs 75 min, x308 4 procs 115 min, all
  hunt recipe; store 169392-169394; candidate files replaced, superseded values stay in git):
  - **s(343) ≤ 18.97823252363561127954…** (exactsolve exact optimum, strict mod 19 flat motions, jammed in every branch;
    VALID both checkers): **−1.667e-2** below Kevin Fang's.  x343 (from the chain's 18.98009) reached 18.97824020 (not a
    local min, MILP dS −1.8e-6; descend.py moved it only 8e-10 then stalled at dS −2.9e-8); x343b from that reached a
    jammed 18.97823252.  Matching distance 0.24 from Kevin's 343, 0.11 from our first 343: a different basin, reached via
    his 344 minus one square.  Chain alone ended at 18.98009.
  - **s(308) ≤ 17.99826987952625587520…** (mod 52, jammed): −1.040e-3 below Kevin Fang's (the polish was −8.44e-4).  x308,
    7 steps from the chain best: cross2 <- cross2 <- aswap <- kick <- rowslide <- band <- cross2; d 0.033 from his file.
  - **s(305) ≤ 17.95113977220699754576…** (mod 44, jammed): −5.64e-5 below itsnaka's (first-hour value −1.26e-5); chain,
    d 0.018 from itsnaka's.
  - 209: chain's 14.9462232264 is itsnaka's raw file converged (−4.3e-7, exact form computed, `runs/hunt3/exact/c209`); not
    claimed.  132, 344 unchanged.  The other 16 n: no improvement in 5.6 h.
  - Reading (corrected 10-10 morning, Evan asked whether it was compute-corrected: it was not).  Proc-hours: chain gave 343
    2.4 (12 rounds: 18.9912 -> 18.9801, the two biggest jumps from removal seeds of best(344)), the dedicated x343 8 (-1.9e-3)
    + x343b 3.75 (-7.7e-6); chain gave 308 1.2 (-5.6e-6), x308 7.7, but its whole gain (-1.9e-4) came at 9 min (0.6
    proc-h).  Follow-up dedicated runs (343 / 344 / 305 / 132, 4 procs x 65 min, ~17 proc-h): zero gain.  So no evidence that
    dedicated runs beat the chain per proc-hour; on 343 the chain was better.  Candidate mechanism (untested): archive memory.
    x308's winning lineage is 7 steps in 9 min (~3 chain rounds of compute) entirely through basins *above* its start
    (17.9997 .. 17.9990 vs 17.99846); x343's is 12 steps over 92 min with a 32-min wait on one parent; chain rounds reach
    lineage depth median 5 / p90 8 and keep only best + 3 carries + 3 library entries.  Test: same starts and proc-hours,
    resetting 6-min rounds vs continuous vs rounds with the archive carried over, several seeds and n.
  Files `candidates/hunt3_n*`.  **Posted 10-10 ~08:55 as jlevy/squares#489** (132, 308, 343, 344; pinned c013f43).  305
  dropped: Couzo's #488 (08:08, before our push) has the identical exact optimum (same 31-digit rational side; both from
  itsnaka's packing via fq), so it is his.  #488 added to `pk/pending.yaml`.
* **n=110 census brought up to date** (10-10 morning, new `census_certify.py`: solves every side below k with packing rows
  but no `basin` row, keeps exactsolve output incl. contacts in `runs/census/n110/<key>/`).  453 solved (+7 certified rows
  that had no stored S): `basin` now 2183 side keys below 11 = 775 certified (**718 distinct exact sides**), 1399 not-min,
  9 unresolved.  Distinct certified within 1e-6 / 1e-5 / 1e-4 / 1e-3 / 3e-3 of the record: 1 / 3 / 9 / 59 / 622.  The
  "second side" 10.996783403 (6.5e-9 above) solves to the record's exact S: an unconverged copy, not a second basin.
  **Caveat:** of the 437 sides first seen in tonight's lineage runs, 394 (90%) are *not jammed* (corner-corner MILP finds a
  first-order descent; fq's fixed-axis linearisation stalls there) and only 41 certified.  So explorer "new basin" counts
  overstate distinct local minima by up to ~10x at 110, and the archive's parents are mostly unconverged points.  Open:
  where do the not-min points go when descended properly (new minima, or known ones)?  Same-basin is a contact-graph
  question; only the 456 newly solved keys have kept contacts so far.
* **Descending the n=110 census's not-min points** (10-10 afternoon, new `descend_census.py`: exact/descend.py's
  descend() with thr 1e-9, then exactsolve; random order; `runs/census/n110/descend.jsonl`, deduplicated
  `descend_dedup.json`).  912 points (all not-min sides below 11 with a stored packing): **605 land on a known certified
  minimum, 62 on a new one (40 distinct), 168 still not jammed, 77 unresolved (mostly 900 s timeouts)**.  Side drop median
  3.4e-7, p90 5e-5, max 8.5e-4; 184 s per point median.  Census now **754 distinct certified minima** below 11; new
  minima are all >= 10.99692 (two at 1.4e-4 above the record, the rest above 10.997).
  - Duplication: the 667 certified landings hit only 155 minima (93 once); clusters of 127, 51, 43, 36, 34 points; the top
    5 minima take 44% of landings.  A copy can sit up to 8.5e-4 above its minimum, so 1e-9 side keys (basin table, chain
    library, explorer archive) do not dedupe them.
  - Bias in starts (L1 / L2 lineage archives, every sub-11 entry mapped to its minimum): parent picks spread over ~85-120
    effective entries but only ~23-62 effective minima.  Minima with >= 3 archive copies (2-7% of minima) take 20-36% of
    picks.  The most-picked minimum in 7 of 8 runs is 10.9968870 (1.0e-4 above the record; 4-13 copies per archive,
    6-18% of picks); the record gets 3-5%.  Elite top-20 holds only 7-15 distinct minima.  Mechanism: elite weights and
    cell top-8 windows are per entry, so a minimum's weight scales with its copy count.
  - Matching distance does not separate copies from distinct minima: copies median 0.004, p90 0.020, max 0.046;
    side-adjacent distinct minima min 0.0001, p10 0.0018.  A 0.01 threshold catches 78% of copies but merges 35% of distinct
    pairs.  "Same minimum" needs the descent (or the contact graph); rattler masking may tighten the copy distances.
  - Not measured: whether the bias hurts.  At 110 every L2 run still reached the record, and 10.9968870 is a near-record
    basin, so part of the concentration is useful exploitation.
* **Insertion test at 110** (10-10 evening, Evan; `insert110.py`, `runs/insert110`; 12 threads, ~1 h).  Question: does
  inserting one square into non-regular sub-11 109 packings add a pile of new sub-11 110 minima, beyond the old technique's
  rate, and any outside the record's family (census: all certified minima within d 0.104 of the record)?
  - Baseline expectation, fixed before the runs: L1 / L2 (old technique, explorer from census starts, 20 proc-h) reached
    113 certified minima nothing else has = **5.7 new minima per proc-h** (9% of their points unmapped, so a bit more;
    expected lower now, since the census grew).  Live baseline arm = same explorer flags from the control's 17 source minima.
  - Seeds (171 in 164 s, 4 hole tries each): control (17 census minima, one square out, back in) 19 / 51 quench below 11
    directly (best 10.99787); donors (40 distinct sub-11 109 packings, not best(109) = 6 + 7/sqrt 2, stratified by side and
    tilted count 29-69) 4 / 120 below 11, 36 in 11-11.02, 27 to the grid.
  - Explorers (6 runs x 2 procs x 30 min, novelty reward vs every known sub-11 side, no elite share, starts below 11.02),
    then certify every new side and descend a random 15 not-min points per arm; new = exact S not in the frozen 758:
    | arm | proc-h | sub-11 sides | certified | new | not-min sampled -> new / known / still | new per proc-h |
    | baseline | 2 | 286 | 207 | 2 (d 0.044, 0.060) | 15 -> 0 / 8 / 7 | 1.0 |
    | control | 1 | 155 | 115 | 0 | 15 -> 0 / 9 / 6 | 0 |
    | donor | 3 | 186 | 136 | 0 | 15 -> 0 / 7 / 8 | 0 |
  - Which donors reached below 11: record-like ones (45-49 tilted: three of them back to the record itself, 10.996783) and
    two low-tilt ones (31, 33 tilted); every sub-11 descendant of the low-tilt donors had rearranged into the record family
    (42-46 tilted, d 0.03-0.063).  The other 35 donors' lineages stayed above 11 (best 11.0004-11.019).
  - Reading: no new family; insertion lands back in the record's funnel and rediscovers known minima.  The live baseline
    is 1.0 new / proc-h, not 5.7: the census is close to saturated for this explorer (or the 30-min novelty runs differ).
    Donors 0 in 3 proc-h vs baseline 1.0 / proc-h: below the baseline, but small counts (P(0 | rate 1) ~ 0.05).  Counter-
    example to "every new technique adds a pile of sub-11 packings" at 110, at this scale.

## Review + first fixes (10-11; `REVIEW-2026-10-11.md`)
* Review (architecture, methods ledger, claims audit, external scan) saved as `REVIEW-2026-10-11.md`; its section 5 is the
  pre-run fix list.  Store backed up (TODO).
* **Melt arm was broken in the explorer since 10-08**: `anneal melt` is NPT and returns a new box side (e.g. 10.997 -> 11.100
  in a smoke test), `explore.mv_melt` dropped it and the screen quenched the melted coordinates in the parent's box.  Every
  explorer `melt` result before 10-11 (be_melt +0.055, hunt recipes with `--melt`) measured that version; `pk.moves.melt`
  (battery / corpus) was right.  Fixed: `explore.melt_side` keeps the side; `mv_melt` stays as the lossy 2-tuple wrapper.
* `pk.pipeline.Policy.same` default 1e-6 -> 0 (polish everything, as explore.py since 10-09); `pk/corpus.py` variants pinned
  to an explicit `R1008 = Policy(same=1e-6)` so old bench names keep their meaning.
* Anneal binary: explorer and `pk/moves.py` used `target-dev/release/anneal` (10-08 working-tree build; it does have `melt`);
  now `target/release/anneal` (current src).
* Temp dirs: explore.py / chain.py now run under `scratch.scope()` (TMPDIR = one `/tmp/pk-<name>-*` dir, removed at exit and on
  SIGTERM / SIGHUP); `explore.work` uses a TemporaryDirectory.  Cleared 769k day-old leftover `tmp*` dirs from /tmp (+24 GB free).
* `pk.engine.FQ.quench(timeout=)` now enforced (watchdog kills the child; error 'timeout'; default 1 h hang guard).
* `explore.py --report` reads n from the archive (was 110 hard-coded).  0-byte decoy `runs/store.sqlite` at the repo root removed.
* **Melt fix A/B** (`melt_ab.py`, `runs/melt_ab`, job record `runs/jobs/melt_ab`; 4 procs, 6.9 min wall = 27.6 worker-min, ~4.2
  cores measured).  Paired: same parent + seed = identical melt, arms differ only in the screen side (old = parent's box,
  the bug; new = anneal's side); explorer screen + polish.  20 parents per n (110: L2 lineage archive entries; 130: be_melt
  base s130 entries) x 8 seeds.  Metrics fixed in the script before launch (no written hypothesis / threshold: do that next time).
  | n | arm | screen-grid | kept non-grid | parent-side return | other side (below parent) | d median (q90) |
  | 110 | old | 64/160 | 96 | 21 | 75 (44) | 0.0057 (0.0095) |
  | 110 | new | 0/160 | 160 | 31 | 129 (71) | 0.0054 (0.0086) |
  | 130 | old | 35/160 | 125 | 43 | 82 (57) | 0.022 (0.061) |
  | 130 | new | 1/160 | 158 | 53 | 105 (79) | 0.021 (0.064) |
  Reading: the bug was a pure loss (22-40 % of melt proposals to the grid), no extra reach (same d; both-kept pairs mostly end
  on the same side).  Fixed melt: 1 grid loss in 320, the most grid-safe move measured so far (10-10 landing spread: others
  lose 23-78 %).  Cheap: the anneal is ~0.1 s at these n (screen + polish sums 26.4 of the 27.6 worker-min), so the 10-08
  "3.3 s per proposal" was the whole proposal.  "Other side" is side-level, not certified: at 110, d ~0.005 is in the copy
  range (10-10 census), so many are the parent's minimum again; at 130 d ~0.02 is real movement.  Census tail (not the
  test's purpose): 110 finals 206 on certified basins, 46 on known not-min sides, 4 unseen sides (~2 distinct, 10.998342 /
  10.998392, +1.6e-3), uncertified; ~9-10 worker-min at 110 (expect ~0.15 new at the live 1 / proc-h).
  Compute units: "CPU-h" / "proc-h" in this log = worker wall-clock (wall x procs, or summed per-item wall in workers),
  never measured CPU; at load > 16 a worker-hour < a core-hour, so cross-run rate comparisons at different load are biased
  (concurrent controls are fine).  Fix with the run-manifest item (getrusage incl. children, load at start / end).
  Open: be_melt s130 / s179 cell with the fixed explorer (~50 worker-min) decides melt on by default in hunt recipes.
* **Run manifests + explorer provenance** (`pk/manifest.py`): every run dir of explore / chain / battery / corpus replay /
  crossing / melt_ab gets `manifest.json` (census_certify / descend_census: `manifest_<kind>_<time>.json` in
  `runs/census/n<n>/`): argv, args, git commit + dirty files, fq / anneal / packer sha256, load at start and end, wall,
  measured CPU incl. children (`cpu_s`), wall x procs (`worker_wall_s`), status ok / incomplete; battery adds the suite
  sha256, corpus the items sha256.  Explorer `proposals.jsonl` adds `seed`, `desc` (move + params), `loosen`, `s_screen`,
  `hit` (archive entry landed on: new, dup or return), `cpu` (measured, per proposal).  `explore.py --replay OUT:LINE`
  re-runs one proposal from parent + seed + mate: 16/16 reproduced exactly (one per kind x status, 1-min test run).
  Test run: cpu_s 120.0 vs worker_wall_s 121.9 at load ~9.6.  Not covered yet: mondrian.py, lineage_reach.py,
  insert110.py top level (their explore.py children have manifests), fq-internal knobs beyond the binary hash.
* **Identity / verdicts unified** (`pk/identity.py`, store tables `minimum` + `copy`, `identity_backfill.py`; tests
  `pk/tests/test_identity.py`).  A minimum = exact S to 1e-17 (neighbour keys merged); verdicts certified / **weak**
  (PSD zero modes or inconclusive second order: counted as certified before) / not-min / unresolved; the f64 1e-9 side key
  stays as the `basin` lookup cache.  All certifiers (corpus, battery, census_certify, descend_census) write through
  `Store.record_verdict` / `add_minimum`; descend_census no longer writes basin rows keyed by float(S) (could overwrite
  another side's row); explorer / chain same-side tolerance = `identity.SAME_SIDE` (2e-9; chain library was 1e-9);
  battery suite `best` precedence bug fixed (ignored the certified known set when a frontier row existed; frozen suites
  unaffected at 110, best = the register record).  Backfill (store snapshot first:
  `~/store-backups/store-2026-10-11b-pre-identity.sqlite`): every certified basin row's second-order status recovered from
  saved verdict files (519 rows had status 'import'); 42 basin rows relabelled weak.  **n = 110 census: 760 distinct minima
  below 11 = 720 certified + 40 weak** (the "754" of 10-10 counted weak as certified and merged 4 distinct minima under
  shared side keys); certified within 1e-4 / 1e-3 of the record: 9 / 57; 1,457 side keys linked as copies to the 760 minima
  (most-copied: 128, 52, 44, 37, 35).  Not yet switched: battery / corpus / rediscovery "new certified" metrics still
  count basin side keys (next: count minima), insert110 (historical).
* **Kick-retry, pre-registered check** (written before launch; `kr_eval.py`, `runs/kr_eval`; 4 procs, expected < 10 min):
  60 not-min census points with a known descend destination + 60 certified minimum points (controls), 4 kicked
  (sigma 1e-4) re-quenches each.  Hypothesis: kick-retry is a cheap substitute for descend before counting / parenting.
  Adopt as the census-v2 arm if (a) >= 50 % of not-min points reach their descend destination (or a lower known minimum)
  within 4 tries at <= 5 s CPU median, and (b) <= 2 / 60 controls move by more than 2e-9.  Otherwise report what the tries
  do and keep descend as the reference.
  **Result** (4.3 min wall, 953 s measured CPU): **fails both criteria; not adopted.**  not-min: after 1 / 2 / 4 tries
  12 / 15 / 20 of 60 at their descend destination, 15 below it after 4, 26 on a known minimum (43 %), 53 moved (median drop
  1.3e-6, max 8.5e-4); CPU median 8.4 s for 4 tries.  Controls (certified minima): **18 / 60 moved, all to lower sides**
  (all on the first try; drop up to 6.0e-4), 54 / 60 end on a known minimum.  Reading: kick + full re-quench (ALM repair
  stage included) is a short-range basin hop, not a finisher: ~30 % of certified 110 minima sit one sigma-1e-4 hop above a
  lower minimum.  As pipeline policy it would erase shallow genuine minima before certification (census undercount).
  `--kick-retry` stays, default off.  Next: kick + polish only (`--no-alm`) as a finisher candidate (same points); the hop
  as an explorer move (judged as a move); census: many shallow minima, so a "minima at kick resolution" count would be
  much smaller than the strict 760 (definition question, Evan 10-11).
* **Decisions (Evan 10-11).**  Census success metric = strict: certified exact local minima (exactsolve: no first- or
  second-order descent direction, to its precision) = table `minimum` cls certified (110: 720; the 40 weak reported beside,
  not counted).  Open (Evan's guess: no): are partially descended (not-min) states useful parents?  Reconstructing them may
  be the job of higher-temperature / annealing proposers.  Candidate rule to test: sub-grid packings farther than some d
  from every packing on record are worth keeping as parents, minimum or not.
* **Not-min states as parents, retrospective** (`parent_retro.py`, `runs/parent_retro.json`; 14 logged runs at 110:
  L1, L2, insert110; ~76.8k proposals; parent verdict by side key, minimum via the copy table; observational, parents were
  chosen by elite / value rules).  Per 1000 proposals, certified children / children on a different minimum than the
  parent's / within-run discoveries: certified parents 635 / 490 / 78; not-min 574 / 478 / 34; above 11: 49 / - / 9.
  Matched (308 pairs: a not-min copy and its own certified minimum, both expanded in one run): 572 vs 607 / 493 vs 493 /
  29 vs 45; but the minimum was expanded first in 236 / 308 pairs and discovery credit goes to the first reach, so the
  discovery gap is mostly order; the order-free rates are equal.  Reading: a not-min copy is about its minimum as a parent
  (no special value, no loss), so duplicate-aware (per-minimum) parent picks look safe.  Hint for Evan's distance rule:
  never-descended not-min parents (no known minimum; 985 proposals) have ~2x the discovery rate of other not-min parents
  (153 vs 77-79 per 1000): small, selected (timeouts, late additions).  Census v2: a controlled arm "keep sub-grid states
  farther than d from every known minimum as parents" with random parent choice inside the arm.
* **Kick + polish only, pre-registered** (written before launch; `kr_eval.py --no-alm`, `runs/kr_eval_noalm`; same seed =
  same 120 points and kick seeds as kr_eval): kick sigma 1e-4, uniform repair, SLP polish, no ALM stage.  Adopt as the
  finisher (before counting / parenting) if (a) >= 50 % of not-min points reach their descend destination or a lower known
  minimum within 4 tries at <= 5 s CPU median and (b) <= 2 / 60 controls move by more than 2e-9.  Else: no cheap finisher at
  this sigma; descend.py stays the reference.
* **Reports count minima** (`Store.side_minima`: side key -> minimum via the copy table): battery `new_cert`,
  `near_new`, `reach.withheld`; corpus report `cert` / `newcert` / `recall` / `<1e-3`; rediscovery.py all counts = distinct
  certified minima (weak excluded); the old side-key counts are printed beside (`*_sidekeys`, `keys` / `newkeys`).
  Battery cross2 vs baseline: bench new 33 vs 32 minima (35 vs 35 side keys), reach withheld 11 vs 2 (unchanged).
  c110a re-reported today: polish everything (`full`) 33 new minima vs 20 for `base` (side keys 37 vs 21); the 10-09 log
  said 30 vs 16: the reports read the live basin table, which has grown since (report drift, review §3): the direction
  holds, the ratio moved 1.9x -> 1.65x.  To do: snapshot the reference tables per report.
  **Result** (2.6 min wall, 565 s CPU): not-min: destination after 1 / 2 / 3 / 4 tries 14 / 29 / 35 / 35 of 60 (none below
  it), on any known minimum 27 / 43 / 51 / 49; CPU median 7.4 s (4 tries).  Controls: 1 / 60 moved (by 1.5e-6); CPU 1.4 s.
  Against the pre-registered rule: safety (b) passes, accuracy (a) passes (58 %), cost (a) fails (7.4 s > 5 s; ~5.5 s at 3
  tries, same accuracy) -> not adopted as written; adoption with a revised cost bound is Evan's call.  The ALM stage was
  what hopped out of genuine minima in the first check.  descend.py: ~184 s median, so kick + polish is ~25x cheaper at
  ~60-85 % of its effect.  `kick_retry(alm=False)` available; explorer flag still uses alm=True (default off anyway).
* **Finisher on the backlog, pre-registered** (written before launch; `finish_backlog.py`, `runs/census/n110/finish_backlog`;
  4 procs, ~1.1 CPU-h).  Question (Evan 10-11): is kick + polish-only (4 tries) + exactsolve a better use of CPU than more
  explorer time?  Break-even framing: new certified minima per measured CPU-h; explorer marginal at 110 ~1 / worker-h;
  descend + solve on 1,012 points: 40 distinct new / 56 CPU-h (0.7 / CPU-h).  Points: the 312 not-min / unresolved census
  side keys descend never settled (168 still not-min after descend, from its output; 78 descend unresolved; 59 never
  descended; 7 unresolved).  Prediction: 2-8 new certified minima per CPU-h (lower than the ~6 estimated from the random
  kr_eval sample, since these are descend's leftovers).  Reading: >= 2 / CPU-h = finisher beats explorer time at 110 (adopt as
  the finisher for novelty / parents, cost scoped to archive-new sub-k candidates); < 1 / CPU-h = no.  (Smoke: 2 points,
  both still not-min.)
* **Out of scope for now (Evan 10-11), recorded to return to:** whether finisher settings, the break-even and the
  measured rates transfer to other n and to other packing families (45-degree diamonds, strips, large n); everything
  above is n = 110.
  **Result** (8.7 min wall, 2,010 s measured CPU = 0.56 CPU-h, median 3.2 s / point): **0 new certified minima -> below
  the pre-registered 1 / CPU-h: not adopted.**  310 points: 57 certified-known + 1 weak-known, 245 still not-min (225 'not
  jammed'), 7 unresolved; by group: still-not-min-after-descend 25+1 / 141 / 1, descend-unresolved 5 / 69 / 2, never-
  descended 25 / 34 / 0, unresolved 2 / 1 / 4.  Expected new at the random-sample rates: ~58 conversions x P(new) ~0.04 =
  ~2.3, P(0) ~0.1.  Caveat (mine, after the fact): the population is descend's hard tail (19 % convert vs 83 % in the
  random kr_eval sample), so this settles "mine the backlog" (no, at 110) and not "finish fresh explorer output inline";
  that goes to census v2 as an arm (finisher on / off, end to end, certified minima per CPU-h).  Open confusion: 225
  points stay 'not jammed' (exactsolve MILP dS < -1e-6) after both descend.py and the finisher: tolerance-degenerate?
* **Frozen references + holdouts** (`pk/refset.py`, tracked dir `refsets/`).  Finals first: `pk.py import-finals --n 110`
  put 18,000 battery / corpus final packings with a certified / weak verdict into the packing table (13,463 + 4,537 new;
  they never reached the store before), so every one of the 760 minima has a representative packing (102 had none).
  Snapshot `census110_1011` (760 minima: 720 certified + 40 weak; sha db691fdafece).  Holdouts (40 explorer starts each,
  novel_ref = visible minima + their copies): `h110_rand50` (random half: 390 visible / 370 withheld, record withheld),
  `h110_near1e3` (side band: the 62 minima within 1e-3 of the record withheld, record included), `h110_asof1008` (visible =
  minima matched by the 10-08 census list: 509 / 251; record visible).  `refset.py score` checked on old lineage runs
  (L2 far1+far2 new: 172 / 370 withheld of h110_rand50 "rediscovered", 0 new, as expected for runs inside the snapshot).
* **Census v2, pre-registered** (written before launch; `census_v2.py`, `runs/census_v2`; Evan 10-11: random-half holdout,
  8 cores).  Holdout `refsets/h110_rand50` (snapshot census110_1011 sha db691fdafece; 390 visible / 370 withheld minima,
  record withheld).  Census recipe (insert110): `--adapt --reward legacy`, novelty vs the visible minima only, no elite
  share, 1 proc per run.  Arms x 2 replicates, all 8 runs concurrent, 60 min each (8 worker-h): control (40 random visible
  starts); melt (+ `--melt 2`, C10); finish (+ `--finish 3`: kick + polish-only on archive-new sub-k results, C11 / C14);
  far (control recipe; starts = the 40 visible minima farthest from the best visible minimum, matching distance 0.078-0.104,
  C05 / X13).  Then import + certify every new sub-11 side; score against the frozen holdout (`pk/refset.py score`).
  Primary metric: withheld minima rediscovered (certified) per measured CPU-h, pooled over replicates; rate ratio vs control
  with a 95 % CI (conditional binomial).  Secondary: genuinely new certified minima (not in the snapshot); not-min share of
  the runs' sub-11 sides; whether / when the withheld record is rediscovered.  Smoke test (1 min x 8): all arms run,
  minute-one rates ~200-600 / CPU-h (meaningless, saturating).
  Predictions: control 40-120 withheld per CPU-h (L2 runs re-scored against this holdout: ~43 / worker-h).  melt ratio
  0.8-1.25 (safe, but one arm of ~12, mid-pack in c110b).  finish ratio 0.7-1.1 (fewer proposals per CPU-h; fewer not-min
  copies) and a lower not-min share of sides than control.  far: if "one funnel" (C05), ratio 0.7-1.3 and the record
  rediscovered in both replicates (C04: lineages cross 0.08 in minutes); if far starts sit in other funnels, more
  genuinely new minima and fewer near-record withheld ones.  New minima: ~1 / worker-h (C08) -> ~2 per arm, too few to
  compare arms (reported, not decided on).
  Decision rule: an arm differs from control if its ratio's 95 % CI excludes 1 and both replicates fall on the same side
  of the control's per-replicate range; otherwise "no difference resolved".
  Not this test (Evan 10-11): "can you find the other funnels" (a holdout by funnel / kind; date holdouts may need to go
  back to ~40 known points); other n (X16).

* **Census v2, result** (10-11, 17:51; evaluated against the pre-registration above, unchanged).  8 runs x 1.00 measured
  CPU-h, all `ok`.  `census_v2.py certify`: 7,779 new entries, 183 unlabelled sub-11 sides (24 certified, 156 not-min, 3
  unresolved; 50 s, no descend).  `census_v2.py report` (`runs/census_v2/report.json`):

  | arm | withheld / CPU-h | ratio vs control (95 % CI) | per rep | within 1e-3 of the record | record found (s) | new | not-min share of sides |
  |---|---|---|---|---|---|---|---|
  | control | 110 | — | 96, 124 | 12, 14 | 81, 542 | 0, 2 | 14 %, 14 % |
  | melt | 101 | 0.91 (0.75-1.11) | 107, 94 | 16, 10 | 253, 425 | 2, 3 | 19 %, 15 % |
  | finish | 84 | 0.76 (0.62-0.93) | 70, 97 | 11, 11 | 3360, 592 | 0, 2 | 4 %, 10 % |
  | far | 82 | 0.75 (0.61-0.92) | 88, 76 | 3, 2 | no, no (best 10.99679327 = best visible; 10.99688984) | 8, 6 | 27 %, 14 % |

  Against the rule: **control** 110 / CPU-h, inside the predicted 40-120 (top end).  **melt**: no difference resolved
  (prediction 0.8-1.25 held).  **finish**: CI excludes 1 but replicate 97 lies inside control's 96-124, so no difference
  resolved; leaning lower, as predicted (0.7-1.1).  Not-min share lower (4-10 % vs 14 %), as predicted.  It also archived
  fewer distinct sides (146, 230 vs 248, 319).  Not adopted: it cleans up output but finds no more minima.  **far**:
  resolved lower (CI excludes 1; both replicates below 96).  The "one funnel" prediction **failed**: no replicate reached the
  withheld record in 1 h, and of the 34 withheld minima within 1e-3 of the record far re-found 3 and 2 vs control's 12 and 14.  The
  alternative prediction held: 14 genuinely new minima vs control's 2.  But all 14 sit at 10.9994-10.99999, the shallow
  edge just under 11.  These are new shallow minima, not a new deep family.
  Coverage: control r0 / r1 overlap 70 → Lincoln-Petersen ~170 withheld minima reachable from random visible starts.
  All 8 runs together found 224 / 370, and far found 31 that none of the 6 other runs found.  So catchability depends on
  the start region: the population is heterogeneous and that estimate is too low.  146 withheld minima were never re-found
  in 8 CPU-h.  New certified minima: control 1 / CPU-h (C08 confirmed in measured units).  The explorer re-finds withheld
  minima ~100x faster than it finds new ones, so the low live rate (X02) is mostly saturation of what these starts reach,
  not slowness.  Next ideas (not run): starts drawn from under-covered regions (far-type parents, Evan's distance rule) for
  new minima; "can you find the other funnels" holdout.

* **Role transitions in census v2 (10-11, analysis only).**  Sub-11 parent → sub-11 result over all 8 runs (5,388
  proposals): 58 % keep the role counts.  20 % change the axis count, so squares move between the tilted bands and the
  axis region (mostly ±1-2, up to ±5).  22 % keep the axis count but change the L/B split.  Crossover changes the axis
  count in 50 % of its outcomes, `band` (movegen add / del / transfer, already an explorer move: 4.4 % of proposals,
  110 new sub-11 / CPU-h, a weak move) 16 %, `aswap` 11 %, kicks / melt / rowslide / mirror 0-3.5 %.  The junction
  family (4 squares at ~45°, (19,19,68,4) in the 10-05 beam keys; 10.99688699 / ...705 / ...984) is reached in all 8 runs,
  control included.  55 % of sub-11 entries have ≥ 3 squares at 40-50°.  **Classifier artifact:** explore.roles cuts at
  45°, so those junction squares split at random into L or B ((19,23,68,0) vs (23,19,68,0) for one structure).  Much of the
  "split change" share is this artifact, and role counts can't see the family.  Junction class needed before any
  region-count holdout.  X13 now leans "side valley": every run enters it, and far_r0 passed through it on its way down
  (10.996887 at 15 min → 10.9967933 at 60 min).

* **Count separation + fingerprint distance (10-11, analysis; `count_sep.py`, `pk/fingerprint.py`).**  Junction-aware
  roles (`explore.roles_j`: L 10-40°, B 50-80°, axis < 3° off, J 40-50°, other).  150 random certified minima of
  census110_1011, all 11,175 exact matching distances (`runs/count_sep/k150.json`, 0.5 measured CPU-h).  Groups by
  (L, B, axis, J, other): within-group distance median 0.017 (p10 0.004, p90 0.066), between groups 0.068.  The nearest
  neighbour is in the same group 66 % of the time (chance 4 %); by axis count alone 95 % (chance 21 %), by J count 90 %
  (25 %).  **Counts separate configurations well**, so a region-count holdout is a real "unseen family" test.  Top
  groups: (19,19,68,4,0) 17 / 150, (19,21,67,3,0) 12, (19,20,67,4,0) 9, (24,22,64,0,0) 8.
  Fingerprint (2-channel smoothed 24×24 centre images, axis / tilted, min over 8 symmetries, L1 / n; ~30 µs per pair vs
  0.2 s exact): Spearman 0.93 vs exact on all pairs.  On per-minimum nearest-neighbour distance: Spearman 0.90, same
  nearest neighbour 77 %, top-quartile overlap 89 %.  Good enough to rank "far from everything known" inside the explorer.
* **External refresh (10-11).**  jlevy#506 (Couzo, franciscouzo@337078ca): beats the store at 175, 304, 305, 338, 339,
  340, 344, 375, 378, 379, plus 380 and 381 (below 20, none in the store).  Not yet verified or ingested; ingestion is
  owned by the records session (TODO "drawn as posted").  Ostashev (dostashev/square-packing@c3d17ab0, new source; clone
  in `_untrusted-third-party/dostashev-square-packing`): files only for 155 / 156 / 211, all already in the store
  (pending:#501).  Its other claims have no files, and its page says equal-or-better packings were already public.
  Method: explicit V-shaped band + integer program for the corner blocks, GPU soft-to-rigid, LP polish, exact KKT.
  Kingbird pages unchanged since 10-10.  Register frontier unchanged.  Full report:
  `notes/refresh-2026-10-11.md`.

* **Census v3, pre-registered** (written before launch; `census_v3.py`, `runs/census_v3`; Evan 10-11: near-record
  holdout, longer far runs, distance-rule parents; 8 cores).  Snapshot census110_1011b (780 minima, sha 5467fe3c04b4),
  holdout `refsets/h110_near_b`: the 62 minima within 1e-3 of the record withheld (record included), 718 visible (best
  visible 10.9978387).  Census recipe as v2, 1 proc per run, 8 runs concurrent, **180 min each** (24 CPU-h).  2 × 2 arms
  × 2 reps: starts far (the 40 visible minima farthest from the best visible, matching distance 0.076-0.096) or rand (40
  random visible) × parent rule plain or sparse (`--sparse-share 0.5`: half the picks weighted by (fingerprint distance
  to the nearest known configuration / median)² / (1 + expanded); bank = the 718 visible minima + the run's archive).
  Questions and predictions:
  (1) Is 110 below 11 one funnel (C05, X13)?  Record = the withheld 10.9967833966.  rand: re-found in both replicates,
  within 30 min (smoke: 10.996793 in < 2 min from 4 random starts).  far: if one funnel with slow mixing, re-found in
  ≥ 1 of 2 far replicates within 3 h (v2 far_r0 was 1e-5 away at 60 min).  Neither far replicate in 3 h (nor far_sparse)
  = evidence of a real barrier at this scale.
  (2) Does the distance rule find more new minima (C13: ~2× hint)?  Primary for sparse: genuinely new certified minima
  (not in census110_1011b) per measured CPU-h, sparse arms vs plain arms (pooled over start type), rate ratio with 95 % CI;
  prediction 1.0-2.0.  Secondary: withheld rediscovered (ratio 0.7-1.2), median nnd of sub-11 archive entries higher.
  Decision rule as v2: CI excludes 1 and, for both start types, both sparse replicates on the same side of the plain
  replicates' range; otherwise "no difference resolved".
  Not this test: the region-count holdout (withhold whole count groups: count separation above says it is meaningful);
  next after this.

* **Census v3, result** (10-11, 21:20; against the pre-registration above, unchanged).  8 runs × 2.99 measured CPU-h, all
  `ok`.  Certify: 664 unlabelled sub-11 sides (36 certified, 4 weak, 618 not-min, 6 unresolved; 127 s).
  `census_v3.py report` (`runs/census_v3/report.json`):

  | arm | withheld re-found (of 62) per rep | per CPU-h | record found at (s) | new minima per rep |
  |---|---|---|---|---|
  | far | 14, 17 | 5.2 | 1952, 2634 | 4, 8 |
  | far_sparse | 22, 29 | 8.5 (1.64×, CI 1.03-2.66) | 3341, 1354 | 11, 9 |
  | rand | 24, 30 | 9.0 | 108, 739 | 3, 2 |
  | rand_sparse | 29, 27 | 9.4 (1.04×, CI 0.70-1.53) | 382, 239 | 6, 8 |

  (1) **Funnel: the far prediction held.**  All 4 far-start runs reached the withheld record, in 23-56 min, vs 2-12 min
  from random starts.  So the configurations 0.076-0.096 out are in the same funnel, slower to drain.  C05 now has
  direct support at this scale, and X13 (second family) is a side valley that every run passes through.  v2's far runs
  missing it in 60 min was the tail of this 23-56 min spread.
  (2) **Distance rule: resolved, 2.0× new minima** (34 vs 17 in 12.0 CPU-h each, 95 % CI 1.09-3.81).  Both sparse
  replicates lie above the plain range for both start types: far 11, 9 vs 4, 8; rand 6, 8 vs 2, 3.  This is at the top
  of the predicted 1.0-2.0 and matches C13's 2× hint.  Secondary: withheld re-found far_sparse 1.64× (both replicates
  above, resolved), rand_sparse 1.04× (no difference).  The nnd comparison is not available: plain arms do not compute
  fingerprints.  To do post hoc if needed.  All 51 new minima are shallow (10.9979-10.9997): the distance rule fills in
  the upper rim, no new deep family.
  Rates are lower than v2's (withheld 5-9 / CPU-h vs 110) because the withheld set is 62 hard near-record minima, not a
  random half: not comparable.
  Decision: adopt `--sparse-share 0.5` in the census recipe.  Next: region-count holdout ("find an unseen family").

* **Night plan 10-11 (Evan: "the whole cpu is yours for the night"), pre-registered before launch.**  16 cores, ~9 h.
  **Census v4, region-count holdouts** (`census_v4.py`, `runs/census_v4`; 8 procs).  `refset.holdout --rule groups`
  (groups = explore.roles_j counts, L / B mirror pairs merged; 46 groups over census110_1011b's 780 minima):
  `h110_J0` withholds every group with no ~45° junction squares (18 groups, 186 minima, record included).  `h110_grand`
  withholds random groups (seed 3) to ≥ 24 % (16 groups, 217 minima, record visible).  `h110_recnb` (record's group +
  every group within 1 in each count: 9 groups, 156 minima) is batch b, after batch a.  Recipe: census v3 + distance rule
  (`--sparse-share 0.5`, bank = visible minima), 40 random visible starts, 1 proc per run; batch a = J0 × 4 reps + grand
  × 4 reps, 240 min each (32 CPU-h).
  Metric: fraction of the withheld minima re-found per measured CPU-h (Poisson rate per withheld minimum), J0 vs grand,
  rate ratio with 95 % CI.  Secondary: J0 record re-found (time), time to the first withheld minimum, new minima.
  Predictions: grand re-finds 25-50 % of its withheld minima per 4-h replicate.  J0 is harder, ratio J0 / grand 0.3-0.8
  (count changes are common: 20 % of sub-11 transitions in v2, but J 4 → 0 needs four junction squares re-assigned).  The
  record is re-found in ≥ 2 of 4 J0 replicates (one funnel, C05).
  Decision rule: J0 harder if the CI excludes 1 and all 4 J0 replicate fractions are below grand's lowest; easier
  symmetric; otherwise no difference resolved.  A J0 ratio near 1 means count groups are no barrier for the census search
  (unseen families are found as easily as hidden ones).  A J0 with no record and ratio < 0.3 means families are real
  barriers at this budget.
  **Time-travel record test** (`chain.py`, `runs/timetravel/{plain,sparse}`; 8 procs).  29 n in 84-211 whose best improved
  after the frozen 10-05 register (`../exact/batch/inputs`), gaps 1.3e-3 .. 3.2e-2.  Frozen truth
  `runs/timetravel/truth.json` (10-05 side, today's store best).  Chain from the 10-05 packings (seeds from 10-05
  neighbours), `--first-pass`, 2 slots × 2 procs × 4-min rounds, 9 h, two arms: plain (chain recipe as hunt1) and sparse
  (+ `--explore-extra "--sparse-share 0.5"`; no bank: distance to the round's own archive only; in chain rounds the
  sparse pick also overrides half of the frontier-lineage picks).  Expect ~9 rounds (~1.2 CPU-h) per n per arm.
  Metric per n: recovered gap fraction (s1005 − best) / (s1005 − s_now), capped at 1; full recovery = best ≤ s_now + 1e-9.
  Predictions (ch2: 126 reached 80 % in 11-14 rounds; bench: most n stuck at 5 min): median fraction 0.1-0.3 per arm,
  full recovery in 0-3 of 29 n.  Sparse vs plain paired by n: no resolvable difference expected (median ratio 0.8-1.5).
  Decision: the distance rule transfers to record hunting if the Wilcoxon signed-rank p < 0.05 and the median paired
  difference > 0.  Anything below s_now is checked against the live register before any word "record".

* **Census v4 batch a, result** (10-12 02:50; against the pre-registration, unchanged).  8 runs × 3.99 measured CPU-h,
  all `ok`.  Certify: 787 unlabelled sub-11 sides (76 certified, 3 weak, 703 not-min, 5 unresolved; 210 s).
  J0 (186 withheld, record included): re-found 89, 90, 82, 59 = 48, 48, 44, 32 %.  First withheld minimum at 32-541 s.
  **Record re-found 4 / 4, at 833-1235 s.**  grand (217 withheld): 107, 126, 112, 125 = 49, 58, 52, 58 %.
  Rate per withheld minimum per CPU-h: J0 / grand = **0.79 (95 % CI 0.69-0.92)**, and all 4 J0 replicates are below
  grand's lowest.  By the rule, **J0 is resolved harder, but mildly**: at the top of the predicted 0.3-0.8.  The record
  prediction (≥ 2 of 4) held, 4 / 4.  New minima: J0 61, grand 74, all shallow (≥ 10.99790).
  Reading: for the census search, a whole unseen count family (no junction squares, incl. the record) is a ~20 % slower
  find than random hidden families, not a barrier.  It is entered within minutes from junction-family starts.  Count
  groups separate configurations (count_sep) but the moves cross between them freely.
* **Census v4 batch b, predictions** (written 02:55, after launch at 02:44 and before any results were read: `recnb`,
  record's group + count neighbours withheld, 156 minima; 8 reps × 240 min, recipe as batch a).  The record is re-found
  in ≥ 6 of 8 replicates.  Per-replicate re-found fraction 30-50 %.  Comparison with batch a's grand is descriptive only
  (not concurrent: different load, and the time-travel chains run alongside both).
* **Census v4 batch b, result** (10-12 06:50).  8 runs × 3.99 measured CPU-h, all `ok`.  Certify (all 16 v4 runs, the
  new sides are batch b's): 441 unlabelled sub-11 sides (46 certified, 4 weak, 387 not-min, 4 unresolved; 96 s).
  recnb (156 withheld, the record's group + count neighbours): re-found 83, 81, 77, 81, 63, 53, 72, 71 = 53, 52, 49, 52,
  40, 34, 46, 46 % (mean 47 %).  First withheld minimum at 61-889 s.  **Record re-found 8 / 8, at 286-3689 s.**
  Against the predictions: record ≥ 6 / 8 held (8 / 8).  Fraction 30-50 %: 7 of 8 replicates inside, r0 just above
  (53 %), the mean inside.  Descriptive only: per withheld minimum per CPU-h, recnb 0.117 vs batch a's grand 0.136 and
  J0 0.108.  152 new minima.
  Reading: withholding the record's own count family and every adjacent family does not hide the record.  It is
  reached in every replicate from starts that differ from it by ≥ 2 in some count.  Together with batch a: at n = 110 count families are
  no barrier to the census search; the cost of an unseen family is ≤ ~20 %.
* **Time-travel interim: n = 154 below the best known** (10-12 06:50; tt_plain, round 26 of the chain from the 10-05
  register).  `runs/timetravel/plain/best_154.txt` side 12.9225664019, which `exactsolve.py` takes to an exact KKT
  point.  Certificate `runs/timetravel/cert154/n154_tt.cert`: **s(154) ≤ 12.922566401934059503377…**, VALID under both
  `verify_cert.py` (separating axes) and `verify_cert2.py` (polygon intersection).  Best known before it: itsnaka
  (SQUISH) 12.9230702023011 (#481, registered as T-136 on main; the witness file still shows the older
  12.9265622459).  Δ = 5.04e-4.  Checked live 06:45: jlevy/squares main @0d6010ff, open issues mentioning 154 (#481,
  #470, #500), Couzo #506 (no 154).  Not imported into the store and not posted (Evan's call; the records session
  owns ingestion).  Note: the chain started from the 10-05 packings and never saw itsnaka's packing, so this is an
  independent route that ends below it.
* **Time-travel record test, result** (10-12 07:47; against the pre-registration, unchanged).  Both arms `done` after
  9 h (543.9 min); plain 218 rounds, sparse 213 (measured CPU 33.4 / 33.2 h, ~1.1 CPU-h per n per arm; the chain's scheduler
  spends rounds unevenly: 153 / 154 took 28-30 in plain, 12 n got ≤ 2).
  Recovered gap fraction: plain median **0.00** (mean 0.076), sparse median **0.00** (mean 0.064); 22 / 29 n not moved
  at all in either arm.  Moved: 108 (plain 0.66, sparse 0.12), 126 (0.04, 0.17), 130 (0.05, 0.12), 154 (1.00, 0.43),
  155 (0.46, 0.61), 180 (0, 0.39), 208 (0, 0.01).  Full recovery: plain 1 / 29 (154, below s_now: see the interim
  entry: certified s(154) ≤ 12.9225664019341, 5.0e-4 below the best known), sparse 0 / 29.
  Against the predictions: median 0.1-0.3, **failed** (0.00: most n do not move from their 10-05 packing at this budget,
  as the bench's "most n stuck at 5 min" suggested more strongly than ch2's 126).  Full recovery 0-3 of 29, held (1).
  Sparse vs plain: better in 5, worse in 2, tied 22; median difference 0, Wilcoxon p = 0.94: **no difference
  resolved**, as predicted.  The distance rule (C13) does not transfer to record hunting at this budget (X16).
  Reading: 9 h of chaining from a 6-day-old register recovers about 7 % of the later community gains, concentrated in
  a few n (108, 154, 155), and once overshoots (154).  A per-n budget of ~1-2 CPU-h does not reproduce what other
  groups found in the following week; where it moves, it moves a lot.
* **n = 154 submission prep** (10-12 08:00).  Re-checked before posting: jlevy/squares main @4a491052 (witness and case
  page 12.9265622459; #481's 12.9230702023 registered as T-136), issues through #550, Ellsworth's catalogue
  (12.9316671296): nothing at or below 12.92256640.  Matching distance from itsnaka's 12.92307 packing 0.319, from the
  10-05 start (Couzo) 0.396: a different local minimum.  Exactsolve: KKT point, min λ 1.5e-4, reduced Hessian PSD
  modulo 45 exact flat motions (realised), 98 corner-corner touches jammed (MILP dS 0), 8 rattlers.
  Regularized (`search/regularize`, coincide style, conservative and fewest policies: identical; input via a scratch
  REFRESH dir so the shared `runs/regularize/exact/n-154` was not touched): same S*, rotation groups 24 → 16, level 2
  (same packing, angle moves path-checked), 80-digit verified; certificate VALID under both checkers.
  Files: `candidates/tt_n154{.cert,.exact.txt,.input.txt,.json}`, `candidates/tt_n154_regularized.{cert,json}`.
  **Posted 10-12 ~08:10 as jlevy/squares#552** (certificate pinned at 597e660).
* **Submission tooling** (10-12, after #552 was assembled by hand).  `candidate.py FILE --name N [--vs ..] [--chain RUNDIR]
  [--stage] [--import --ref ..]`: exactsolve, both checkers, read-only live check (register witness + case page, issues
  naming n with their nearby decimals and authors, Ellsworth's catalogue, store frontier), matching distances, optimality
  facts, chain lineage, regularized view (private REFRESH dir), staging into `candidates/`.  Regression: it reproduces
  every number of #552 in one command (`runs/candidates/test_n154/`; the live check lists #552 itself as the only tie).
  `chain_lineage.py RUNDIR N`: lineage through chain rounds.  **154 descends from Couzo's 10-05 155** (12.9556195921,
  one square removed at r0141), then 17 productive rounds / 36 moves at 154; #552's credit paragraph said "154, 155
  and 156" (correct but loose; comment drafted).  Other traces: plain 108 from 110's start (two squares removed), sparse
  155 and plain 155 from 156, sparse 180 from 181.  `chain.py` now logs each seed's source side (`seeds_src`).
  Store import of 154 not done (shared-store write blocked; TODO decision).
* **154 in the store** (10-12).  Evan: packings are verifiable, so add them without anyone's approval; provenance is
  best effort; there is no separate records session (the "owned by the records session" notes above were wrong).
  Imported: the winning round `r0192_n154` (archive, parent links) and the certified exact point (id 291398,
  `candidate:tt_n154`, ref #552); `pk/pending.yaml` has #552 (pinned 597e660, claim 12.922566401934060), and
  `sync-pending --issue 552` deduplicated the posted certificate onto it.  Frontier 154: ours (pending #552).

## Grid-border and diamond hunt (10-12; Evan: goal = number of new packings n ≤ 324; most interesting: the grid border, where a new packing is the first of its kind, and the diamonds 89 / 109)
Store refreshed first: #506 (Couzo, 15 new packings) via pending.yaml; register @25de5a4 (no new witnesses; n = 51
witness fails to parse: a list-valued field, look later).  Live starts: `pk.py export 2..326 --dir runs/live1012`
(`mcmin.load_deg` now reads n/d sides).
Border targets = the first "= k" entry of each column c of the margin table k − s(k² − c) (WISHLIST P1/P2; columns
shrink with k, so a later k in the same column is harder: 309 behind 274, 111 behind 90, 212 behind 183):
**90 (c = 10), 133 (11), 157 (12), 183 (13), 242 (14), 274 (15)**; 308 (c = 16) fell 10-10.  Extrapolating each
column one step: 183 ≈ −0.004, 274 ≈ −0.006, 242 ≈ −0.015, 157 ≈ −0.03, 133 ≈ −0.04: all long shots.  Route:
Fang's 308 grew the 307 two-wedge packing by one square; the chain's insertion seed (seed_k −1: scale up 0.2-1.2 %,
add a square at a pose hole, quench) is the closest operator we have.
* **H1, positive control: does the chain's insertion find 308 < 18 from the 10-05 register?**  `chain.py --ns-file
  [307, 308] --seed-ks 1,-1`, 10-05 inputs (308 = grid, 307 = 17.98103), 2 slots × 2 procs × 4-min rounds, 2 h
  (~8 CPU-h).  Success: any 308 below 18 (certified).  Prediction: 30 % (insertion must pull a 307 + 1 packing from
  above 18 to below in one quench; Fang's wedge cut is a targeted operator).  If it fails: the border needs a
  wedge-cut / stack-cut operator (new machinery) before more border compute.
* **H2, border + diamonds, live frontier**: `chain.py --ns-file [89, 90, 109, 133, 157, 183, 242, 274] --seed-ks
  -1,1,2 --src runs/live1012`, 2 slots × 2 procs × 6-min rounds, 2 h now (8 cores while another session runs),
  extended when cores free up.  Not a test: a hunt.  Success = any certified packing below the live frontier.
  Prediction: border 0 of 6 in this budget; diamonds 89 / 109 (sources 88 at 9.8825 and 108 at 10.9048 leave 0.067 /
  0.045 of slack for one inserted square): 25 % for at least one.
Chain vs dedicated (X15 / C23) is not run now; the chain logs (state.json events with t, per-round manifests with
measured CPU, seeds_src) suffice to compare against a later dedicated arm.
* **X01 investigated** (10-12; the 225 n = 110 points 'not jammed' after descend.py and the finisher;
  `runs/census/n110/finish_backlog/results.jsonl`; scratch analysis, 30-point sample of the 103 with a store packing).
  - Where: 144 / 225 in one side band 10.998342-10.998343 (1e-6 linkage), the rest in ~26 smaller clusters.  All 225
    have certified minima within ±1e-5 in side (median 12; median 9 around a certified minimum: the region is dense
    everywhere, no special comb).  Matching distance to the closest of the ~16 side-nearest certified minima
    (census110_1011b reps): 0.001-0.024 (single-move reach ~0.05; within-count-group median 0.017).  Side above that
    minimum: −2.7e-6 … +2.3e-5, mostly +1e-7 … +9e-7 (8 / 30 below it).
  - Not a near-contact artefact: the corner-corner MILP dS is unchanged when contacts up to 1e-3 (not 1e-6) enter the
    model (exactsolve find_contacts + vv_branch_test at 1e-6 / 1e-5 / 1e-4 / 1e-3).  Real first-order descents, but
    shallow: dS −4e-7 … −2.3e-5, equilibrium residual 1e-7 … 6e-6 ("almost KKT").
  - Realisable gain is tiny: |dS| × distance to the nearby minimum 3e-9 … 2e-7, the finisher's moves ~1e-8 (median),
    43 / 225 not moved at all.
  - `exactsolve --force-kkt` (contact set up to 1e-3, Newton decides): 22 / 30 solve failed; 8 converged to KKT points
    1e-7 … 5e-6 from every known minimum side, uncertified (saddles / degenerate, not checked).
  Reading: slow valleys, not tolerance.  The points sit 1e-7 … 1e-6 above, and 0.001-0.02 away from, known minima, on
  gentle real slopes our step-based descenders give up on.  For census counting they are near-copies (the finisher found
  0 new among them).  What they drain into (a known minimum or an uncertified one in the dense band) needs a descender
  that follows shallow MILP directions with exact contacts (active-set line search); not worth building for the count.
* **H3, s(109) from the 110 family** (10-12 10:00; Evan: the diamond 6 + 7/√2 looks out of place, and the non-diamond
  109 packings are real and close: use them as seeds, certification not a priority).  The store holds 2,871 distinct
  non-diamond sub-10.99 sides at 109 (29-45 tilted squares; roles L / B / axis ≈ 22 / 14 / 73: the s(110) two-band
  structure with a shorter second band; best 10.9518715, +2.1e-3 over the diamond; visual: claude.ai artifact
  TsFtMvxHD3bsy3PY2b811g).  `chain.py --ns-file [109] --seed-ks 1,-1 --src runs/live1012 --carry` the 3 best
  non-diamond packings, `--lib-init` 60 more (best 40 + 20 spread), 1 slot × 2 procs × 6-min rounds, 6 h (more slots
  when cores free).  Not a test.  Success: a certified 109 below 10.949747468.  Prediction: 20 % in 6 h (12 CPU-h);
  the family's best is 2.1e-3 away and at 110 the record sat ~1e-3 below most of the family's mass.
* **s(90) near the grid** (10-12; Evan: are there rigid non-grid packings just above 10, in the two-band family, or do
  they collapse to the grid?).  Store: 4,216 distinct non-grid sides in [10, 10.05), all uncertified polish output;
  none in (10, 10.005), 68 in [10.005, 10.01), 1,131 in [10.01, 10.02).  Best classes (roles_j L / B / axis / J / O):
  one dominant band (B 17-20) plus 0-10 L and 2-6 J, e.g. (4, 20, 64, 2, 0) at 10.0090456; the balanced two-band class
  (10, 18, 62, 0, 0) at 10.0095668.  exactsolve on the best of the 15 best classes: 10 certified local minima
  (10.0090596 … 10.0101264, strict modulo 20-36 flat motions, 1-8 rattlers), 5 not jammed with real descents (dS
  −0.0015 … −0.018).  descend.py + exactsolve on those 5: 4 land on certified minima at 10.0090596 (same S* as the best
  certified one), 10.0090720, 10.0092177, 10.0092301; the best side 10.0090456 descends to 10.0089808, still not
  jammed (dS −2.3e-7, an X01-type slow valley).  Nothing collapses toward the grid: rigid non-grid minima sit 0.009
  above it.  Scratch: runs in the session scratchpad (n90/).  Many rattlers (up to 8) in the best ones.
* **H1 / H2 result** (10-12 10:38, 2 h each, by the pre-registration).  **H1 failed**: insertion into 307 brought the
  308 lineage to 18.0005167 in 13 rounds (18.0027 after one), never below 18; plain insertion gets within 5e-4, the
  crossing needs a targeted operator (Fang's wedge cut).  **H2: no candidate.**  Closest border lineages: 274
  17.0084153, 242 16.0133959, 133 12.0164375; none at 90, 157, 183 (2-3 rounds each).  Diamonds: 89 lineage
  9.9522574, 109 10.9518715 (the stored family best; not improved).
* **16 cores from 10:40** (the other session finished).  Pre-registered before launch:
  - **H3 restarted** (state kept): 1 slot × 6 procs per round, `--explore-extra "--sparse-share 0.5"`, 6 h.  (A
    single-n chain runs one round at a time, so more slots do not help.)  Prediction unchanged (20 % for s(109)).
  - **H4, s(90)**: `--ns-file [90] --seed-ks 1,-1`, carry the three best (10.0089808 descended / not jammed, 10.0090596
    and 10.0090720 certified), library = 12 more rigid or descended n = 90 packings from the near-grid check, 1 slot
    × 4 procs × 6-min rounds, 6 h.  Success: any 90 < 10.  Prediction 2 %.  Secondary (descriptive): the best non-grid
    side reached (now 10.0089808).
  - **H5, count sweep**: the 107 non-grid neighbours (n ± 1, 2) of the 56 n improved since 10-05 (`runs/h5_ns.json`),
    excluding 89 / 90 / 109; live frontier, `--seed-ks 1,2,-1 --first-pass`, 3 slots × 2 procs × 4-min rounds, 8 h.
    Success = new packings below the live frontier (certified via candidate.py).  Prediction: 1-4 (time-travel: 1 full
    + 6 partial recoveries in 29 n at ~1.1 CPU-h each; here ~0.45 CPU-h per n, but the seeds are fresh).
* **12:15 notes.**  (a) **H3b crashed at launch** (10:39): the stopped H3 had left a half-written `r0007_n109`, and
  the resumed chain's `os.makedirs` failed on it, so ~1.5 h of 6 procs sat idle.  chain.py now renames such a dir to
  `*.partial<t>` (lineage skips it).  Relaunched 12:14 for 4.5 h (ends ~16:45, same as H4); prediction unchanged.
  (b) **First H5 "candidate" was a false positive**: n = 207 at 14.880960371559 vs register 14.880960373044 is
  BuddenPing's #503 packing re-polished (candidate.py: both checkers VALID at S* = 14.8809603715556, matching distance
  0.000 to the store best, case page shows 14.88096037).  Not a new packing; not imported (the store keeps #503).
  chain.py's CANDIDATE flag now needs a gain over `CAND_TOL = 1e-8` (it was 1e-9, inside polish noise); the
  running H5 still uses the old rule, so the watcher filters to > 1e-8.
* **H3b / H4 result** (17:15).  **H3b (s(109), 48 rounds): no candidate**; the lineage stayed at 10.9518715 (the
  stored family best, 2.1e-3 above the diamond), so nothing moved in 4.5 h.  Prediction was 20 %: a miss.  **H4 (s(90),
  59 rounds): no candidate** (prediction 2 %).  Secondary: the best non-grid side went 10.0089808 -> **10.0024260487**
  (`runs/hunt1012/h4/carry_90_0.txt`).  candidate.py (`runs/candidates/h4_n90_nongrid`): both checkers VALID,
  exactsolve says a strict local minimum modulo 7 flat motions (min lambda 8e-7, small), matching distance 0.684 to
  the grid.  So rigid non-grid 90s exist 0.0024 above the grid, 4x closer than the 10.009 seen this morning.  Live-check
  bug: the register value for n = 90 is the string "'10'" and fails to parse; fix in candidate.py.
* **H5b pre-registration** (17:21; Evan: give the idle cores to the sweep, wrap up at 18:40).  H5 cannot change slots
  live, so a second chain, `runs/hunt1012/h5b`: the 90 H5 targets that have had only their first-pass round
  (`runs/h5b_ns.json`; H5 keeps working the ones with active lineages, 52 / 83 / 206-208 / 68 / 53 / 84-87), same
  settings (live frontier, `--seed-ks 1,2,-1`, 4-min rounds) with 5 slots × 2 procs, 1.2 h so it drains by ~18:40.
  Success: a candidate with gain > 1e-8 that candidate.py certifies as a new packing.  Prediction: 0-1 (about one
  more round per n, ~0.13 CPU-h each).  Evaluated together with H5 at 18:40.
* **H5 / H5b result** (18:44).  **No new packings** (prediction H5 1-4: a miss; H5b 0-1).  H5: all 107 targets, 279
  rounds, 43.4 CPU-h measured; H5b: 62 of its 90 targets, 62 rounds, 10.0 CPU-h.  The only flag was the n = 207
  polish false positive.  Closest lineages to the register (lineage best, not the register packing): 71 +1.9e-12
  (re-found the record, presumably the same packing), 53 +7.3e-7, 208 +1.7e-6, 68 +1.2e-5, 211 +7.3e-5, 110 / 151
  +1.0e-4, 234 +1.6e-4.  Reading: neighbours of freshly improved n are already near their own basin floors; seeds
  from n ± 1, 2 re-find the records, they don't beat them.  The day's compute (H3/H3b 26.8, H4 23.6, H5 43.4, H5b
  10.0 CPU-h; plus H1/H2) gave 0 new packings.  The useful by-products are the rigid non-grid s(90) minimum at
  10.0024260 and the 109 family picture.  More plain chain compute on these targets looks unproductive.  The next lever
  is new machinery (the wedge-cut insertion operator for the border, the shallow-valley descender).
* **s(110) census brought up to date** (19:00).  Today's hunts ran three n = 110 rounds (H5 r0021 / r0270, H5b
  r0008).  `pk.py import-run --below 11`: 155 sub-11 entries, 150 new packing rows, but only 1 new 1e-9 side key,
  and that one is not-min.  So the census is unchanged: 4,312 side keys (958 certified / 3,277 not-min / 51 weak / 26
  unresolved, 0 unsolved), 890 certified + 49 weak distinct minima below 11; best 10.9967833966 (the record), 9 within
  1e-4 of it.  Every other chain run today was at n != 110.
* **Richness / "does a better s(110) exist", incidence statistics** (10-12 evening, analysis only; `census_incidence.py`).
  32 census replicates (v2-v4, 1-4 CPU-h each) as sampling units; a minimum's incidence = runs that reached it (starts
  excluded).  827 of the 890 certified minima reached.  Chao2 976, ICE 928 overall (lower bounds under heterogeneity),
  so ~40-90 unseen, mostly in the 2e-3 … 1e-2 band (Chao2 by band: <1e-5 3 = known, <1e-4 9 = known, 1e-4 … 1e-3 52
  vs 48 known, 1e-3 … 3e-3 766 vs 697, 3e-3 … 1e-2 153 vs 136).  Sample coverage 98.5 % (unseen minima ~1.5 % of
  incidences).  Incidence falls with depth band (median 22 / 13 / 16 / 6 / 4 of 32), but not strictly: two minima at
  +4.46e-5 / +4.58e-5 were reached by 1 / 2 runs.  So deep-and-rare exists.  Detection bound: a better minimum missed
  by all 32 runs has per-run detection p < 0.09 (95 %), i.e. rarer than 7 of the 9 minima within 1e-4.  No richness
  estimator speaks to the extreme member; that needs a tail model (side distribution of the low-incidence class), not done.
* **Baseline: compute per new record** (10-12 evening, ledger from this log; pre-10-11 figures are wall × procs, see
  line 1244, some estimated from logged settings; H3 / H4 / H5 / tt are measured manifests).  Record events n ≤ 324
  (below the live frontier when found): 8 (266, 270, 272 on 10-06; 132, 155 on 10-08; 132, 308 on 10-10; 154 on 10-12);
  155 was a tie with Couzo #451, 270 and the first 132 were later beaten.  Record-hunting compute ~540 h in all.
  | window | compute (h) | events n ≤ 324 | h per event (95 % Poisson CI) |
  |---|---|---|---|
  | 10-06 sweeps / kick censuses (lightly polished large n) | ~76 | 3 | ~25 (9-120) |
  | 10-08 hunt1 + hunt1b | ~69 | 2 (1 tie) | ~35 (10-300) |
  | 10-09 h2r1, hunt2 (+ big1/2, n > 324 only) | ~53 | 0 | – |
  | 10-10 Hunt 3 + dedicated (+ 343, 344 above 324) | ~126 | 2 | ~63 (17-520) |
  | 10-11 / 10-12 time travel, H2-H5 | ~178 | 1 | ~180 (32-7,000) |
  | all | ~540 | 8 | ~68 (34-157) |
  Falling with time (25 → 35 → 63 → 180 h per event) while the frontier is also moving under other people's work.  A
  single-rate Poisson fit is optimistic for the next one; the 10-11/12 window, the only one on current methods and
  frontier, says ~1 per 100-200 h at best, with a CI that cannot exclude ~1 per several thousand.

## Plan 10-12 evening (PLAN-2026-10-12.md; Evan: preliminary evaluations in ~32 core-hours)
* **fq: rectangle container + rounded squares** (`--rect-dx/--rect-dy`: container [0, s+dx] x [0, s+dy], s minimised;
  `--round-r0 R [--round-compress] [--hold-f F] [--round-steps K] [--round-p P] [--round-kick SIG]`: soft-to-rigid
  pre-stage = inner half-side 0.5 - r with margin r on every separator / wall (still smooth in the lifted variables),
  optional compression at r0 with the side free, then the box held at F x that side while r goes to 0 through K
  fixed-side overlap relaxations, then the usual quench + polish).  Checks: gradcheck 1.8e-8..4.6e-8 (square, r 0.2,
  rect dx 1.5 + r 0.1); s(110) record --no-alm polish unchanged (10.99678339663382); 110 in an (s+1) x s box 10.7058488
  (both walls checked); `repair` grows s so both rectangle walls scale (it diverged for dx > half-width).
* **Soft-to-rigid positive control, pre-registered** (`round_ctl.py`, `runs/round_ctl/r1`; written before launch).
  Seeds = #575's: 241 Couzo 15.98813244 (store 347; #575 -> 15.987134249), 211 de Winter 14.99796070 (store 94; #575 ->
  14.993694187).  Grid 552 runs per n: r0 {0.03 .. 0.3} x compress {y, n} x hold {1, 0.999, 0.997} x steps {10, 40}
  x p {1, 2} x kick {0, 0.01 x 3 seeds}; ~9 s each, ~2.8 CPU-h total.  Success (per n): any output > 1e-8 below the
  seed side.  Secondary: below the #575 value; which knobs.  Prediction: 241 yes (30-60 %: their move started from this
  exact packing), 211 yes (40 %), below #575 in neither.  If neither n moves: our rounding stage is not their move
  (look at their description again: "every square soft, box held") before any wider use.
* **Event bench harness** (`pk/eventbench.py`, `bench/events.json`; PLAN item 0): frontier_at(n, T), as-of-T start sets
  for chain.py `--src`, `score` with `--certify`.  Validation on the 10-05 time-travel runs reproduces PACKER's numbers
  except where later events moved the truth (Couzo #570 now leads 154 by building on our #552: plain 154 = 0.99).
  Diagnosis from the same runs: **179 had 1 round of 218** (153 / 154 got 28-30); its 180 - 1 removal seed quenched to
  13.9224 vs the cutoff frontier 13.8953.  Budget allocation, not (yet) an operator verdict.
* **Removal positive control, pre-registered** (`removal_ctl.py`, `runs/removal_ctl/r179`; before launch).  Every single
  removal from the 10-05 register's s(180) (Couzo 13.9278881) x completion modes plain / shrink (`fq --shrink`, our
  guess at SQUISH's shrink-hopping, never tested on its motivating case) / soft (round-r0 0.15 compress) / kick 0.02:
  720 fq runs, ~2 CPU-h.  Success: any 179 below the cutoff frontier 13.8953410700 by > 1e-8 (SQUISH #401 got
  13.8913127 this way).  Prediction: no mode alone gets below (one probe: plain 13.9199, shrink 13.9195, soft 13.9142);
  soft best of the modes.  Follow-up if the best removals are within ~0.02: explorer 15 min x 2 procs from the best 4.
* **Results 20:05 (controls, against the pre-registrations above).**  Soft-to-rigid r1 (17.8 k CPU-s): **0 / 552 below
  the seed at 241, 0 / 552 at 211** (best = the seed itself / +1.8e-4): our rounding pre-stage is not #575's move, as
  written.  Removal r179 (5.8 k CPU-s): no mode below 13.8953411; best soft 13.9142096, kick / shrink 13.9194628, plain
  13.9199427 (prediction held: soft best, none below).  Rectangle library lib1 (4.0 k CPU-s): 1,559 crops re-optimised,
  median gain over the crop's start side 0.05, 955 pieces beat the grid in their own box after Pareto pruning.
  **Found on the way:** composing the 10-10 crop atlas (formula composer) and *polishing the ties* gives s(129) =
  11.8700138504, Ostashev's #562 value (posted 10-12): composer v0 counted it as a tie with the 11.8720298 frontier
  and never polished it.
* **Follow-ups, pre-registered (20:08).**  (a) Explorer 15 min x 4 procs from the 4 best distinct removal outputs
  (`runs/removal_ctl/r179_ex`): success = 179 below 13.8953411; prediction 25 %.  (b) Soft-to-rigid grid g2 (Ostashev's
  README: "soften slightly, hold the box at its side while they re-harden, then release and compress"): no compression,
  hold 1.0, r0 {0.05, 0.1, 0.2} x steps {10, 40} x kick {0.03, 0.05, 0.1} x 10 seeds, 241 + 211 (`runs/round_ctl/r2`):
  success as r1; prediction: still 0 (if so, park the move; it is not reproducible from the description).
* **Results 20:45.**  r179_ex (explorer 15 min x 4 procs from the 4 best removals): best 13.9102726, not below
  13.8953411 (prediction 25 %: miss).  Soft-to-rigid g2 (Ostashev's wording, larger kicks): **0 / 180 at 241, 0 / 180
  at 211**.  Soft-to-rigid as implemented is parked; #575's move is not reproducible from its description with our tools.
  Composer on the 10-10 crop atlas, every layout constructed + polished (`rectlib.py construct`), scored with
  eventbench (now with the monotone envelope: s(n) <= s(m), m > n; candidate.py's live check too): benchmark successes
  at cutoff 10-11 05:51 = 129 (Ostashev #562's value exactly), 375, 380 (s(m) + k layouts below Couzo's later posts);
  live = 327 (certified: S* = s(123) + 7 = 18.5913781454977, both checkers VALID, strict local min mod 188 flat motions,
  83 rattlers, 5.9e-3 below Couzo; not posted), 364 (= 327 + 1 row/col), 417-420 (our frontier above 400 is thin:
  unverified).  The library adds nothing at n <= 324 at the formula level except 129.
* **Squeeze test, pre-registered (20:46)**: the 15 merged-atlas compositions within 5e-3 above frontier(10-11 05:51)
  (53, 86, 110, 151, 175, 199, 204, 211, 227, 234, 261, 266, 295, 299, 308), constructed + polished, then explorer
  10 min x 2 procs each (~5 CPU-h; `runs/rectlib/lib1/squeeze`).  This is the step SQUISH / Ostashev run after
  composing.  Success = below frontier(T) (benchmark) / frontier_now (live).  Prediction: 1-2 of 15.
* **Squeeze test, result (21:06)**: 15 runs x 1 proc x 20 min in parallel (Evan: use the idle cores).  **0 / 15 below
  frontier(T)** (prediction 1-2: miss).  11 / 15 did not move at all: several starts sit at integer sides (110 at 11,
  211 at 15, 308 at 18: grid-jammed L layouts); movers 151 (to the frontier, +2.6e-12), 175 (+3.6e-4), 204 (+3.3e-5).
* **Where the composition gap is: piece quality.**  Library best for Ostashev's decisive cell, 38 squares in
  (h + 2) x h: ours 5.8712 (re-optimised crop), his 5.7642 from random starts in the rectangle (his best crop 5.8203).
  Our pieces are crops (+ relaxation); his winning pieces were searched from scratch in their own rectangle.
* **Session tally (10-12 evening, ~19 measured CPU-h of the 32 budget).**  Benchmark successes: 129 (= #562), 375, 380
  (n > 324, s(m) + k).  Live: s(327) <= 18.5913781454977 certified (not posted).  Negatives: soft-to-rigid (0 / 1,464
  runs), removal 180 -> 179 (best 13.9103 vs 13.8953), formula composer with the library at n <= 324 (only 129), squeeze
  of 15 near compositions (0 / 15).  Lessons: polish every composition (129 was a tie in composer v0); the frontier
  needs the monotone envelope; the explorer cannot leave integer-side L layouts; next lever = from-scratch search in
  the rectangle for small (m, d) pieces (fq rect mode exists; the explorer needs rect plumbing).

## Rectangle pieces from scratch (10-12 night; `rectsearch.py`)
Multistart + MBH in the (h + d) x h box via `fq serve` (QOpt rect_dx, round_r0 / round_compress added to pk.engine).
Smoke (38 squares, d = 2, 2 proc-min each): hard random starts best 6.0 (grid), soft-to-rigid starts (r0 0.1-0.45,
compress) 5.8211, + kicks 5.7996.  ~0.09-0.2 s per quench.
* **Positive control, pre-registered**: cell (38, 2), 14 procs x 10 min (2.3 CPU-h), soft share 0.8.  Success = h <=
  5.76417431749 (Ostashev #528, from random starts on a GPU).  Prediction 30 %.  If it fails: report the best and the
  start-mode / kick yields, then decide between more budget and a better start model.
* **Control result (21:19)**: **passed.**  5.7641743175 (= #528's piece, to 1e-9) first reached at 264 s in one of 14
  workers (~1 CPU-h of search at the first hit), 769 hits in all; every hit came through MBH kicks from soft-to-rigid
  starts (soft starts alone 5.7888, hard random starts 5.8717).  So our CPU stack reproduces the GPU piece search.
* **Cell sweep, pre-registered (21:20)**: m = 10..60 x d = 1..4 (204 cells), 1 proc x 3 min each (~10 CPU-h,
  `runs/rectsearch/sw1`), soft share 0.8; then export -> merged atlas (crops + lib1 + sweep) -> composer + construct +
  polish -> eventbench at cutoff 10-11 05:51 and now.  Metrics: cells where the search beats the best crop piece;
  composer results below frontier_now at n <= 324.  Predictions: search beats crops in > 70 % of the cells both cover;
  composer wins at n <= 324: 0-3.

## Soft squares (log-aspect homotopy; 10-12 night; design `literature/soft-squares-2026-10-12.md` §2, Evan: pursue at
modest odds)
* fq `--soft-a0 E [--soft-compress] [--hold-f F] [--soft-steps K] [--soft-p P] [--soft-kick SIG]`: per-square
  b_i, a_i = E tanh(b_i), w = e^{a/2}, h = e^{-a/2} (area 1) as extra ALM variables in a pre-stage only (polish
  untouched); gradcheck 2.2e-8 .. 3.0e-8 with soft on (also with rect + rounding); 110 record polish unchanged.
  pk.engine QOpt soft_* fields; pk/crossing.py methods `soft:EPS:HOLD:STEPS:KICK` (EPS 1e-9 = matched rigid
  control) and `rectwalk:DX:J` (rigid squares in an (s + DX) x s box, DX walked back to 0 in J quenches).
* Pilot n = 54 grid starts, 6 trials per arm: 0 crossings in every arm; soft compress reaches 7.86-7.98 (< k = 8),
  hardening at the held box leaves overlap 0.003-0.027, the rigid finish returns to the 8.0 grid: the failure mode
  §2.4 predicts (reversible deformation, no topology change at kick 0.02).  The `above` class is empty at 54.
* **Test, pre-registered (22:00; runs after the cell sweep)**: n = 54, 71 (grid + above where available), 100 trials
  per arm: quench; soft:1e-9:1.0:20:0.1 (control); soft:{0.05, 0.2}:1.0:20:{0.02, 0.1}; soft:0.2:0.98:20:0.1 (box
  held 2 % tighter); rectwalk:1:10, rectwalk:2:10 (`runs/crossing/soft1`).  Metric: crossed (side < k, no full line).
  Kill rule (the literature session's): 0 crossings in all soft arms at 400 trials and line count unchanged in
  >= 90 %; at 100 trials, 0 crossings in every arm = stop here unless the trajectories show topology changes.
  Prediction: 0-2 crossings in total (the design's own estimate is 10 % at 54 / 71, mine is lower after the pilot).
* **Soft-square test, result (22:35)** (`runs/crossing/soft1`, 1,800 trials, ~2 CPU-h): **grid starts, n = 54 and 71:
  0 crossings in every arm** (incl. rectwalk 1 and 2); 99-100 / 100 end on the grid in every arm, so the line count
  did not change.  `above` class at 71 is degenerate (plain quench takes 100 / 100 below 9, to 8.9632675): there the
  soft arms are *worse* than the matched rigid control with the same kicks (rigid 72 / 100 below 9; soft 0.05: 64,
  0.2: 46, 0.2 held 2 % tighter: 36).  Against the rule: kill (0 crossings, lines unchanged).  Prediction (0-2) held.
  Soft squares parked; the container-axis variant (literature session) was conditional on line changes, which did
  not occur.

## Rectangle pieces: sweep and #528 rebuilt (10-12 night)
* **Cell sweep result** (`runs/rectsearch/sw1`, 204 cells x 3 proc-min, 36.6 k CPU-s): the search beats the best
  re-optimised crop in 90 / 136 cells both cover, ties 23, loses 23 (66 %; prediction > 70 %: just missed).
  3 min per cell is not enough for quality: (38, 2) reached 5.7996 here vs 5.7642 with ~1 CPU-h in the control.
  Formula composer on crops + lib1 + sweep: nothing new (129, 327 again).
* **#528's s(175) rebuilt.**  Two copies of our searched 38-square piece (control, h = 5.76417431749) in opposite
  corners (B transposed) + the exact MILP upright fill (`mondrian.exact_fill`) give **175 squares at S = h + 8 =
  13.7641743175 = Ostashev's s(175)** (#528, 10-12 07:26; frontier at the 10-11 05:51 cutoff 13.7671190: a benchmark
  success).  The formula composer cannot see it (guillotine columns: 170 by hand).  So the #528 method is reproduced
  end to end on our side; the weak link was the composer, not the pieces.
* Bug fixed in rectsearch.py (elite sorted only on some events; no effect on any best.txt: all equal their log minima);
  adaptive kick (`--adaptive`, Gensane-Ryckelynck) added, untested.
* **Composer v2** (`compose2.py`: two pieces in opposite corners + exact MILP upright fill at a target side; pieces
  tight for S: S - dim within 0.03 of an integer; area bound; tilted-overlap check).  Positive control at the 10-11
  05:51 cutoff: **175 HIT** (37 tight pieces, 35 fills, polished 13.764174317490 = #528).
* **Live run, pre-registered (22:45)**: pieces = lib1 + sw1 + c38 (1,044), every sub-grid n in 20..340, S =
  frontier_now(n) - 1e-9 (eventbench, envelope), 14 procs (`runs/compose2/live1`).  Success = count >= n, then
  candidate.py (both checkers + live check).  Prediction: 0-2 hits (the frontier at n <= 324 has had this method
  applied by Ostashev / SQUISH with their libraries; ours differs mainly in the re-optimised crop pieces).

## Night plan 10-13 (budgets; each run evaluated against its budget afterwards)
Principle (Evan 10-12): new progress comes from new technique, more compute only loosely; measurement only to steer
compute; every run has a CPU budget and is graded on what it got within it.  Calibration from EVENTS: Ostashev's few
searched pieces were each reused for several records (38-sq (h+2)xh: 175, 232, 264; 44-sq, 62-sq: 176, 204, 303), so
the leverage is in which cells get searched and in reusing every good piece.
* **Code (no CPU)**: compose2 pulls arms within a time budget per n through `pk/bandit.py` (KindBandit / CellBandit
  moved from explore.py, now shared), every fill logged; pieces live in the store (`piece` table, `Packing.dx`,
  `pk/pieces.py`; 1,762 imported) and rectsearch keeps every distinct local minimum (h to 1e-7) there; compose2
  `--probe D`: fill at frontier + D and bisect each fitting arm to its smallest side (near-miss = how much piece
  improvement a record needs).  Control: 175 at the pre-#528 frontier found in 2-24 s (3 seeds) vs 89 s exhaustive.
* **P1 store sync** (agent; ~0 CPU) + crops of the 14 new frontier packings -> pieces (budget 2 CPU-h).  Graded on:
  new pieces that beat the cell's best.
* **P2 bandit composer live**, n 20..340 sub-grid, 120 s per n (budget 6 CPU-h).  Success = hit -> candidate.py.
  Prediction 0-2 hits.  Kills the old exhaustive c2live when it starts (P2 supersedes it; its 113 n are logged).
* **P3 probe scan** at frontier + 0.02 (budget 6 CPU-h): near-miss list per n.  Its only value is choosing P4's cells
  and finding polish candidates (need < 2e-3: fq polish the composition; success = below frontier after candidate.py).
* **P4 targeted piece search** on the cells P3 ranks (budget 15 CPU-h, <= 1.5 CPU-h per cell; a cell stops early once
  its h meets the requirement, then compose at once).  Graded on hits per CPU-h and the fraction of cells whose
  requirement was met.
* **P5 (contingent)** removal from several near-best n + 1 parents (the 179 miss), budget 3 CPU-h.
* **c2live stopped (00:50)** at 113 / 159 n (20..262), 0 hits, ~25 CPU-h measured-ish (14 cores x 115 min): the
  exhaustive pair order stalled on n = 263 for 40 min.  Superseded by the bandit composer.
* **P1 result**: store synced (#528-#575; only #570's s(154) file missing, never pushed).  New frontier beyond the
  events table: 304 (Ostashev 10-13, 17.9294277581), 207 (Mishapolk #573, 1.5e-9).  frontier_now now includes the
  store envelope; our s(327) <= 18.5913781455 is in the store (certified, ours:candidates).  Crops of the 50 new
  packings re-optimised (`runs/rectlib/lib3`, 1,007 re-optimisations, 3.3 min x 14 procs ~ 0.8 CPU-h of 2 budgeted):
  **290 beat the cell's best piece, 96 new cells**; the big gains are at d = 6..8 (e.g. (50, 8) by 1.74, (28, 6) by
  1.44), cells our sweep (d <= 4) never covered.  Lesson: the library had a hole at large d.
* **P2 + P3 combined, pre-registered (00:55)**: compose2 `--probe 0.02 --budget 150`, n 20..340 sub-grid, 14 procs,
  all store pieces (`runs/compose2/probe1`).  Hits (need < 0) are polished and go to candidate.py.  Budget 7 CPU-h.
  Prediction: 0-2 hits (new pieces at d 6-8 make 1-2 a bit likelier than before); the near list feeds P4.
  *Restarted 01:05* after 5 min: + one-piece arms (B = none) and 104 frontier packings (n 10..250) as d = 0 pieces
  (control: 327 = s(123) + piece found on the first fill at the pre-ours frontier).  Same budget and prediction.
* **P2 + P3 result (probe1, 25.8 min x ~14 cores ~ 6 CPU-h of 7)**: 159 n.  3 hits, all n > 324: **325, 326 ->
  18.5714285714 (frontier = envelope of our s(327) 18.5913781455; nothing posted at 325 / 326 in the store), 339 ->
  18.9268989004 (frontier 18.9294338838, gain 2.5e-3)**; candidate.py running.  At n <= 324 no hit; near arms are
  mostly exact ties (need 0: the composer rebuilds how the frontier at 125, 149, 171, 174, ... was made); strictly
  positive needs >= 3e-3, mostly small well-searched cells.  Prediction 0-2: 3 (all above 324).
* **Seam polish of tight arms, pre-registered (01:50)**: `compose2.py polish runs/compose2/probe1 --max-need 2e-3
  --per-n 3` (fq polish, all squares free).  Hypothesis: a tie built from other pieces than the frontier's is a
  different packing and may relax lower.  Success: s_pol < frontier_now - 1e-9 -> candidate.py.  Budget 1 CPU-h.
  Prediction: 0-1 (the squeeze of 15 near compositions on 10-12 gave 0).
* **Certified (candidate.py, both checkers VALID, live check clean, in the store as ours:candidates)**:
  - **s(325) <= 130/7 = 18.571428571428571...** (exact; `runs/candidates/comp2_n325`): s(50) (53/7, frontier piece) at
    one corner + a searched (36, 2) piece (sw1, 3 min) at the other + exact upright fill (side = s(50) + 11).
    Nothing posted at 325 (register stops at 324; no issue through #575); previous bound = envelope of our s(327)
    18.5913781455 (gain 2.0e-2).  Local min, 54 rattlers.  326 (same side, (41, 3) piece) in exactsolve.
  - **s(339) <= 18.926898900356742** (`runs/candidates/comp2_n339`), Couzo #506 18.9294338838: gain 2.53e-3.  Two
    (h + 7) x h crops of Ostashev's 10-13 s(304) update (lib3, the crops made 40 min earlier) in opposite corners.
  Not posted (Evan's call; drafts later).
* **Seam polish result** (5.1 min x 12 procs ~ 1 CPU-h): 71 polishes; below the frontier only at 325, 326, 339 (the
  probe hits again); every tie stayed a tie (gains within +-1e-10).  Prediction 0-1 new: 0.  Ties don't relax.
* **P4 not run as designed**: the probe's positive needs at n <= 324 sit on rigid structural sides (227 at 109/7 vs
  15.5711, 308 at exactly 18): shrinking a piece would not move them, so the targeted piece search has no target.
* **Relax-a-composition control, pre-registered (02:20)** (`relax.py`).  Ostashev's README: the composer's 304 layout
  (17.939518) relaxed to 17.929428 "as its aligned blocks shifted" when every square was softened.  Our 304
  compositions, generated without pieces from #575's s(304) (crops of the answer), start at 17.947-17.959 (5 arms).
  Operators: rigid polish, fq shrink (SQUISH-style), kick + shrink, rounded soft r0 0.05 / 0.15 with the box held,
  soft + compress, loosen + kick; 4 seeds.  Success: any result < 17.929566 (pre-#575 frontier).  Secondary: drop
  > 5e-3 from the start (the mechanism exists in our tools).  Budget 1 CPU-h.  Prediction: no success (our starts
  are 1e-2 worse than his), secondary 30 %.
* **Relax control result (r304, 125 runs, 0.48 CPU-h)**: no success (best 17.94156 vs 17.929566; prediction held).
  **Secondary met**: soft operators drop compositions by up to 9.9e-3 (soft15c; soft05 8.6e-3, soft15 8.4e-3, loose
  6.7e-3; 18 / 80 soft runs > 5e-3), rigid polish / shrink by 3e-5.  Ostashev's 304 drop was 1.0e-2.  So the B-family
  step exists in our tools; our 10-12 control applied it to frontier packings, not to compositions.
* **relax1, pre-registered (02:30)**: every probe1 near arm at n <= 324 with 1e-6 < need < 0.015, up to 3 distinct
  pairs per n (95 inputs, 45 n), ops soft05 / soft15 / soft15c / loose x 6 seeds (`runs/relax/relax1`).  Success:
  s < frontier_now(n) - 1e-9 -> candidate.py.  Budget 7 CPU-h.  Prediction: 1-3 successes (the drop distribution
  reaches 1e-2; needs below 5e-3 exist at ~15 n, e.g. 53, 106, 227, 299 at 3.5e-4 .. 7e-4).
* **relax1 result**: 0 successes (prediction 1-3: miss).  74 of 95 inputs did not move under any soft op (jammed
  against grid-like upright blocks); the drops that exist are rare (104: +1.3e-2 relaxed to an exact frontier tie;
  205 5.6e-3 -> 1.6e-3; 53 7.0e-4 -> 3.5e-4).  Closest positive: 299 (2.8e-4), 106 (3.2e-4), 53, 227 (3.6e-4).
* **relax2, pre-registered (02:50)**: soft-relax the probe's exact ties at n <= 324 (need <= 1e-6: compositions that
  rebuild the frontier side from other pieces; 28 n: 38 67 82 101 104 122 125 145 147 149 170 171 173-175 197 198
  200 201 226 229 231 257 260 262 290 291 296), up to 3 pairs per n, same 4 ops x 6 seeds.  A tie starts at the
  record, so any soft drop > 1e-9 is a success -> candidate.py.  Budget 4 CPU-h.  Prediction: 0-1 (most of these
  frontier packings are themselves rigid compositions, s(m) + k type; relax1 showed 78 % of inputs don't move).
* **relax1 cost** 11.4 min x ~14 cores ~ 2.7 CPU-h (of 7).  **relax2 result** (8.6 min, ~2 CPU-h of 4): 1,656 runs,
  **no tie moved by even 1e-9** under any soft op.  0 successes (prediction 0-1 held).  The tie compositions (and so
  the frontier packings they rebuild) are rigid: s(m) + k / grid-row constructions with no slack for softening.
* **Flat frontier steps**: only 2 n <= 324 have s(n) = s(n + 1): 147 (12.6568542495) and 290 (17.5355339059).
  **Removal 148 -> 147, pre-registered (03:05)**: every single-square removal of s(148) (store id 59) x {shrink,
  soft05, soft15c, loose}, seed 1 (`runs/removal_ctl/r147`).  Success: < 12.656854249 - 1e-9.  Budget 0.5 CPU-h.
  Prediction: 0 (a classical construction; Couzo's and our chains have done removal at most n).
* **relax3, pre-registered (03:15)**: iterated soft relax (`--rounds 4`: each input restarts from its best result;
  monotonic basin hopping with soft kicks) on the only n <= 324 whose compositions moved in relax1 and still sit
  near: 53, 106, 151, 205, 227, 299 (up to 2 pairs each), ops soft05 / soft15 / soft15c / loose x 6 seeds per round.
  Success: < frontier_now - 1e-9 -> candidate.py.  Budget 3 CPU-h.  Prediction: 0-1 (closest gaps 2.8e-4 .. 3.6e-4
  after one round; drops shrink as a composition settles).
* **Removal 148 -> 147 result**: 592 runs, all exactly 12.6568542495; 0 (prediction held).  326: exactsolve failed
  (SVD did not converge) on the probe's layout; the seam-polish copies (same side) are being certified instead.
* **P6 probe above 340, pre-registered (03:25)**: compose2 `--probe 0.02 --budget 150 --grid`, n 341..480 (48 below
  grid, 92 at the grid: beating the grid is a record too), all pieces (`runs/compose2/probe2`), then seam polish of
  hits.  Success = candidate.py-certified packing below frontier_now.  Budget 6 CPU-h.  Prediction: 3-10 certified
  (325-340 gave 3 / 16 targets; the grid n add easy ones near (k - 1)^2).  Runs after relax3 (cores).
* **relax3 result** (4 rounds, 1,056 runs, ~1.5 CPU-h): 0 (prediction 0-1 held); 299 relaxed to an exact frontier tie
  (+6e-10), the rest stalled at their round-1 values.  Soft relax at n <= 324 is done for this library.
* **s(326) <= 130/7 certified** (`runs/candidates/comp2_n326_1`, the seam-polish copy; the probe's own layout failed
  exactsolve).  In the store.
* **P6 result** (probe2, 26.6 min x ~15 cores ~ 6.6 CPU-h of 6): 140 targets, **37 hits** (n 362-446; gains 3e-4 to
  8.6e-2 against frontier_now, which above 340 is mostly the grid or generic constructions).  Prediction 3-10:
  exceeded.  All 37 going through candidate.py (`cert2`).  Caveat before any claim: frontier_now above 340 knows only
  the store (issues, catalogue, register); general constructions in the literature (Göbel-type formulas, Friedman's
  survey) may beat some of these and must be checked n by n.
  (Times in this section after 01:00 are approximate; the machine clock read 02:19 at the cert2 check.)
* **cert2**: 17 / 37 certified on the probe's own layouts; 20 failed exactsolve ("no valid certificate") or are slow.
  Seam-polish copies (`runs/compose2/polish2`, 100 below the frontier) are being certified for those 20 (`cert3`).
* **Counterfactual probe, pre-registered**: compose2 `--cf-eps E` scales every rectangle piece (d > 0) to h (1 - E)
  (square frontier pieces untouched: a better square piece is a better s(m)), then asks the exact fill whether n fits
  below the frontier.  A cf hit names the cells whose improvement by E would give a record at n <= 324: the targets
  P4 lacked.  Run: n 20..324, E = 0.005 (what ~1 CPU-h of search bought at (38, 2) over 3 min), 60 s per n.
  Budget 3 CPU-h.  Then rectsearch bandit on the named cells (budget <= 8 CPU-h, targets = the cell's h (1 - E)
  or the bisected requirement).  Success (of the pair) = a certified record at n <= 324.  Prediction: cf names
  5-20 cells; the pair yields 0-1 records.
* **cf1 result** (22 min x 6 procs ~ 2.2 CPU-h of 3): **no counterfactual hit at n <= 324** (two tie artifacts with
  square frontier pieces, 200 = s(148) + 2 and 229 = s(147) + 2, within the fill's 1e-6 tolerance).  Every rectangle
  piece 0.5 % better still composes no record below 325.  So targeted piece search (P4) stays off: two-piece
  composition at n <= 324 is saturated for pieces at our quality + 0.5 %.  Prediction (5-20 cells): miss.
* **Refill test, pre-registered**: for each n in 20..340, the store's best packing's tilted squares + `exact_fill` of
  the uprights at the same side.  count > n at side f(n) with f(n + k) > f(n) = a record at n + k (insertion by
  refilling).  Budget 0.3 CPU-h.  Prediction: 0 (finders' upright blocks are normally maximal).
* **Refill test result** (1.1 min x 10): 321 n; count > n at 151 n but only where s(n + k) already equals s(n); 0
  records (prediction held).  Frontier upright blocks are maximal.
* **Three-piece compositions (built, not run)**: `compose2 --three` (LL, UR, UL corners).  Area screen kills it at
  n <= 324: the waste budget at the frontier is f^2 - n (~15 at n = 175) and each tilted piece wastes ~6-7 against its
  box (38 in 7.76 x 5.76), so 0 / 20,000 random triples pass at 175.  A composition below 325 carries at most ~2
  tilted blocks; this is also why two-piece composition is saturated there.
* **Bug: store_frontier dropped catalogue rows with pen NULL** (`pen<=1e-9` is false for NULL): Ellsworth's s(368) =
  19.7071067812 was missing from the envelope.  Fixed (`pen IS NULL OR`).  Consequence: **365, 366, 367 are ties
  with s(368), 368 is worse; not records.**  Every other hit re-checked against the fixed envelope (excluding our own
  new candidates): gains unchanged.
* **P6 certified so far: 29 / 37; records (excl. 365-368): 362 363 364 370 375 401 402 403 404 409 410 411 412 414
  417 418 419 443 444 446 454 456 457 459 460** (+ 413 442 450-453 455 458 in exactsolve).  In the store.  Many are
  s(m) + k arithmetic in untouched territory (443 = s(123) + 10; 409 / 411 / 413 = s(129) + 9 with Ostashev's piece).
  Before any claim: literature general constructions above 340 (Friedman DS7, Göbel) must be checked n by n.
* **Caveat on "saturated"**: probe1 tried a median 340 of ~14,000 arms per n (2.4 %), cf1 fewer.  At frontier + 0.02
  the best sampled count is n - 1 or n at every n <= 324, so the claim rests on sampling.
* **Depth test, pre-registered**: compose2 probe 0.02 at 10x the budget (1,500 s per n, 1 proc each) at 53 86 106
  129 151 205 227 234 299 308 (closest positive needs / most arms).  Success: a hit (need < -1e-9) -> candidate.py.
  Secondary (the decision this buys): min positive need improves > 2x at >= 3 of the 10 n => sampling depth matters
  and deep composition sampling deserves real compute; otherwise two-piece composition is done at n <= 324.
  Budget 4.5 CPU-h.  Prediction: no hit; improvement at <= 2 n.
* **Literature check above 324 (agent; Ellsworth's pages, Friedman DS7, AMS 2021): the L-step closure kills most
  of the n > 324 "records".**  Friedman's L-step (n squares in side s, s not an integer -> n + 2 floor(s) + 1 in
  s + 1) applied to the published n <= 324 values beats or ties 26 of 36: e.g. s(290) + L = 18.5355 at 325 / 326 (vs
  our 130/7), s(298) + 2L at 370.  Our frontier never applied it, and the composer missed it because frontier
  pieces stopped at n = 250.  **Fixed: `eventbench.frontier_now` = events + store, closed under the L-step, then the
  monotone envelope** (`_l_closure`, n <= 1200).  Re-checked against it (our own candidates excluded):
  **records standing: s(339) <= 18.9268989004 (2.5e-3), s(364) <= 19.5837822040 (7.6e-3), s(375) <=
  19.8824510305 (1.7e-2 vs Göbel s(376)), s(414) <= 20.8984034740 (1.1e-3).**  Not records: 325, 326, 362, 363,
  365-370, 401-404, 409-413, 417-419, 442-444, 450-460; 446 ties s(294) + 4L.  **Our 10-12 s(327) 18.5913781455 is now
  a tie**: s(292) = 17.5913781455 was posted since (synced tonight), and s(292) + L = s(327) exactly.
  P6 graded against the corrected frontier: 4 records, not 25+; prediction (3-10) held only by luck.  cert2 / cert3
  stopped (the remaining n are not records).
  Lessons: (1) the frontier must include every cheap closure we know (L-step, envelope, s(m) + k is the composer's
  one-piece case) before calling anything a record; (2) the probe2 / cert spend above 340 (~13 CPU-h) was mostly on
  n the L-step already covered; (3) the earlier "s(325) / s(326) certified record" claims in this log are withdrawn.
* **Depth test result** (deep1, 25 min x ~9 cores ~ 3.8 CPU-h of 4.5): 10x the arms (e.g. 106: 3,172 vs 207 fills)
  improved the closest miss at 1 / 10 n (151: 2.2e-3 -> 3.9e-5), and that one is a crop of Ellsworth's own s(151)
  refilled to the same side (polish and 2 rounds of soft relax: tie within 3e-12).  No hit.  Prediction held.
  **Decision: two-piece composition (our pieces, + 0.5 %, 10x sampling) is done at n <= 324.**
* **Night tally (10-12 22:25 .. 10-13 ~04:30)**: measured ~29 CPU-h with manifests + c2live ~25 (stopped) +
  certification ~10 => ~64 CPU-h.  Records: s(339), s(364), s(375), s(414) (n > 324; certified, in the store, not
  posted).  None at n <= 324.  Built: bandit composer (shared pk/bandit.py), pieces in the store, probe / rank /
  polish / counterfactual / three-piece modes, relax.py (soft relax works on compositions: 304 control), refill
  test, rectsearch bandit over cells (built, not run: no targets survived), frontier closed under the L-step.
* **Pre-post live check (10-13 08:20 MDT)**: synced #582 (Miles Gray s(381), claims only) and #583 (Couzo; pin 4596670
  pushed ~14:00Z: 154 211 241 302 339 379-382) + Couzo HEAD's unannounced n375 19.8939065342.  Public-only L-closed
  frontier vs ours: 339 +4.85e-4 (vs #583 18.9273836385), 364 +7.60e-3, 375 +1.15e-2 (vs Couzo HEAD), **414 now beaten**
  (Couzo 375 + L = 20.8939 < our 20.8984; our 375 + L = 20.8825 is the corollary).  Certificates committed (231d41d,
  `candidates/comp2_n{339,364,375}.*`); posted as jlevy/squares#584 (Evan OK, 10-13).
  364 = s(146) [BuddenPing #566] + 7; 375 = s(88) [SQUISH] + 10 (crops of s(88), s(129) [Xu]); 339 = two crops of
  Ostashev's s(304).  `pk/sources.py`: pending.yaml entries with `repo: null` are claims only.

## Day plan 10-13 (Evan: "run it in order"; track CPU estimate vs actual per run)
Order: (1) fresh-packing pass, (2) band-merge relax control on 241, (3) V-band constructor (plan to Evan before the
sweep), (4) effort map as a pk command; compute: estimator check on logs, 2 % + shrink what-if, border-n screen.

| run | estimate CPU-h | actual CPU-h | result |
|---|---|---|---|
| F0 fresh-pass control | 3 | ~3.4 (crops 0.56, compose 2.53, certify ~0.3 serial; 444 still certifying) | 3 / 3 recovered (pass) |
| B0 band-merge relax control | 2 | 0.39 | 0 / 144 runs moved (drop; duplicate of 10-12 g2) |
| W2 what-if 2 % + shrink | 2.5 | 2.71 | 0 cells (stop rectangle search at n <= 324) |
| G1 grid-boundary screen | 0.8 | 1.08 | 0 / 9 below the grid (composition can't reach the border) |
| A1 arm composer, cross-n control | 1.5 | 1.61 | miss: nothing below the grid side (anchoring bug + hard control) |
| A2 arm operators, +1 square from n-3..n+3 | 1.1 | 1.31 | miss: best count n - 1 in 4,780 evals |
| H1 holdout rebuild, n = 110 (160 of 890 certified minima) | 2.7 | 3.4 (+0.1 aborted first launch) | 0 / 160 rebuilt |
| E0 basin-entry baseline (rand + above, 8 + 8 x 60 min) | 16 | 16.04 | 16 / 16 entered; median ~11 CPU-min |
| R1 race: control / finish3 / nomirror / notop3 | 16 | 18.6 (overshoot: budget ignored running replicates; fixed) | all open; notop3 0.53 [0.25, 1.16] |
| V-band fit / rebuild checks (interactive) | - | ~0.1 | ceiling test: own tilted set + fill = n at 6 / 6 band records; idealised arms 17-23 short |

* **F0 fresh-pass control, pre-registered (09:45)**.  `fresh.py run --ids 291615 291590 291569` (Ostashev's s(304)
  #575, BuddenPing's s(146) #566, Xu's s(129) on the register: the packings 339 / 364 / 375 used) `--target public
  --exclude-src pending:#584` (as if we had not posted), budget 30 s per n (last night 150 s), 16 procs, n 20..480.
  Hypothesis: the pass finds what the broad sweep found, from the new packings alone.  Metric: certified hits below
  the public frontier at 339, 364, 375.  Rule: >= 2 of 3 = adopt as the standing post-sync step; <= 1 = look at why
  (budget, piece selection) before using it.  Estimate 3 CPU-h (crops ~0.2, composer ~2.5, certify ~0.3).
* **B0 band-merge relax control, pre-registered (09:45)**.  `relax.py run --ids 347 94` = the pre-#575 s(241)
  (15.9881324395, Couzo / our exact optimum; `pk.bands`: bands 48 + 32) and de Winter's s(211) (14.9979607050; one
  connected group of 74, two angle families by eye: 18-20 and 24-29 deg), ops shrink, kshrink, soft05, soft15, soft15c,
  loose; 6 seeds, 2 rounds (MBH from the best), 16 procs.  Ostashev (README, GPU soft-to-rigid): 241 -> 15.987134
  (one V of 85), 211 -> 14.993694 (one V of 76).  Hypothesis: our fq soft/shrink operators on CPU do the same merge.
  Metric: best side per input.  Rule: 241 <= 15.9876 (half his gain) or 211 <= 14.9959 = the operator works ->
  sweep every frontier packing at n <= 324 with >= 2 bands / families (estimate first); neither moves > 2e-4 = drop
  (our soft relax is not his soft-to-rigid; note what differs).  Estimate 2 CPU-h (144 fq runs x ~50 s).
  (`pk/bands.py` new: bands = connected tilted groups (link 1.6), families = folded-angle clusters; it does not
  separate de Winter's 211, whose bands touch.)
* **F0 result (10:40)**: 3 / 3 recovered from the three new packings alone (prediction >= 2: pass; adopt).  Hits below
  the public-minus-#584 frontier: 339 18.9268989004 (= ours), 364 19.5837822733 (polish 2.7e-8 above our certified
  S*; exactsolve found no certificate from it, a re-polish moves only to 2.6e-8: finishing gap), 375 19.8857007500 (a
  different arm than ours, 3.2e-3 above our 19.8824510305, still 8.2e-3 below the public frontier), plus 403 / 444 (L-step
  corollaries reached by composition).  Certified: 339, 375, 403 (live check: tied / beaten by our own #584, as
  expected).  CPU ~3.4 of 3 estimated, wall 15 min + serial certification (fixed: candidate.py runs in parallel now;
  summary written before the import check; exactsolve failure -> re-polish + retry once).
* **B0 result (10:20)**: 0 / 144 relax runs moved 241 or 211 by more than 1e-9 (soft15c / soft05 only got worse).
  Rule: drop.  **This repeated 10-12 round_ctl g2 (0 / 180 at each)** -- I had not checked PACKER before
  pre-registering; 0.4 CPU-h wasted.  Three tests now: our fq rounding soft-to-rigid does not do #575's band merge.
* **V-band constructor, first build (`vband.py`, `pk/bands.py`)**.  Findings so far (no search run yet):
  - Record arms (241 pre-/post-#575) are staircase blocks: 5-7 parallel tilted columns of 7-10 squares along the arm,
    column starts sheared (taper).  The two arms are separate bands at 241 old (48 + 32) and one connected band after the
    merge (85); de Winter's 211 is one connected group already (its bands touch).
  - **Fill ceiling**: the record's own tilted squares + `mondrian.exact_fill` give exactly n at 241 (x2), 211 (x2),
    156, 176.  The fill is expressive for this family (unlike Xu's 126).
  - **Idealised arms lose ~20 squares**: fitted arms (unit spacing, one angle per arm) keep 72-76 of 80-85 tilted and
    the fill then gets 218-224 of 241 at S + 0.01..0.1 (239 at S + 0.6): small deviations in the band boundary cost
    whole upright rows.  Same failure as layout.py (10-04) and stairgen (10-09).  Ostashev bridges it with GPU
    soft-to-rigid, which we cannot reproduce (3 negatives).
* **W2 what-if: 2 % better pieces + shrink tail, pre-registered (10:50)**.  `compose2.py run --cf-eps 0.02 --probe 0.03
  --budget 60 --ns 20-324 --procs 16` (`runs/compose2/cf2`; every rectangle piece's h scaled by 0.98, d = 0 pieces
  unchanged).  A cell "matters" if some arm with it has need < shrink allowance = 7e-4 x S (p90 of the measured
  post-construction shrink, 10-13 table).  Hypothesis (from cf1 at 0.5 %: 0): a short list of cells, mostly small m.
  Rule: empty list -> stop rectangle search at n <= 324; non-empty -> those cells (and only those) get the cell bandit
  (<= 10 CPU-h, separately pre-registered).  Estimate 2.5 CPU-h.
* **G1 grid-boundary screen, pre-registered (10:50)**.  `compose2.py run --grid --only 90,111,133,157,183,212,242,274,309
  --probe 0.05 --budget 300` (`runs/compose2/g1`): the first grid-best n per k.  Hypothesis: area-infeasible
  everywhere (no arm fits n even at k + 0.05).  Metric: best need per n.  Rule: all misses > 0.02 = composition is the
  wrong tool at the border (band construction only); any need < 0.005 = look at it.  Estimate 0.8 CPU-h.
* **Estimator check (logs, 0 CPU)**: 40 composer arm pairs at the same n differing only in one piece's h: dS_min / dh is
  1 (25 %: that piece binds) or 0 (75 %: the fill / other piece sets the side), nothing between.  A piece improvement
  moves the side one-for-one only where it binds, so all-pieces what-ifs are upper bounds; credit binding cells only.
* **W2 result (11:15)**: with every rectangle piece 2 % smaller, 36 arms come within the shrink allowance (7e-4 S) at 7 n
  (38, 67, 82, 101, 125, 211, 308): **all of them are whole square frontier packings + uprights** (s(m) + k ties, d = 0,
  not scaled); 0 cells with a binding rectangle piece.  Prediction (short list): miss.  Rule: **stop rectangle-piece
  search at n <= 324**.  The cell bandit is not run.  2.71 CPU-h of 2.5.
* **G1 result (11:25)**: at 90, 111, 133, 157, 183, 212, 242, 274, 309 (5 min each, 112-162 arms filled per n) the
  smallest fitting side is exactly the grid's k (S_min = k - 1e-9, the fill tolerance): every arm that holds n needs the
  whole box; none fits below.  Prediction (composition is area-infeasible below the grid at the border): held.  The
  border needs band construction (or nothing reaches it).  1.08 CPU-h of 0.8.
* **Day tally so far (11:25)**: F0 3.4 + B0 0.4 + W2 2.7 + G1 1.1 + interactive ~0.2 = **~7.8 CPU-h measured** (est.
  8.3).  Adopted: fresh pass.  Closed: band-merge relax with our operators (3rd negative), rectangle-piece search at
  n <= 324 (W2), composition at the grid border (G1).  Open: V-band constructor finisher (see V-band notes above).
* **Arm composer (`arms.py`, 10-13 after Evan: "real arms and variations of them")**.  Library `runs/arms/lib.pkl`:
  6,760 arms (8 symmetries) from 479 band packings, n 80..324, <= 4 distinct per n.  Operators on real squares: column
  end extend / trim, edge column add / drop, arm shift; ends near the other arm or a wall weighted 4x.  Own-arm rebuilds
  (no operators) give exactly n at the packing's side for 241 (x3), 211, 156, 176 (layout tolerance 1e-6, as exact_fill).
* **A1 cross-n control, pre-registered (12:05)**.  `arms.py run --only 241,211 --exclude-n` (no arm from a packing of the
  target n), dn 25, ds 0.6, max-ops 3, 8 reps x 300 s per n, 16 procs, finish = fq polish of the best 3 per n
  (`runs/arms/a1`).  Hypothesis: arms from neighbouring n, varied locally, rebuild a record-quality packing.  Metric:
  best finished side.  Rule: <= the pre-#575 record (241: 15.9881324395; 211: 14.9979607050) at either n = pass ->
  sweep the k(k+1) clusters (plan to Evan first); both above by > 2e-3 = look at what the operators miss before more
  compute.  Prediction: pass at one of the two (50 %); below today's frontier at neither (90 %).  Estimate 1.5 CPU-h.
* **A1 result (12:15)**: miss.  0 layouts below the grid at either n: the best that holds n needs side 16.0000077 (241) /
  15.0000042 (211), i.e. the fill gives up the band.  1.61 CPU-h of 1.5.  Diagnosis (0 CPU): (1) corner anchoring
  breaks the elbow junction when the side changes: 240's own arms in 241's box give 239, 239's in 240's give 238, 210's
  in 211's give 209; (2) consecutive records are not +-1 edits (240 -> 241: same 85 tilted, columns redistributed
  between the arms; 155 -> 156: 46 -> 67 tilted).  Fix: an arm pair from one packing is placed in its source frame
  (box grows at the top / right walls only): 240 -> 241's box keeps 240 in 8 / 8 symmetries, 210 -> 211's box 210 in 8 / 8.
* **A2, pre-registered (12:20)**: can the operators add the one square?  `--only 241,211 --exclude-n --dn 3` (sources
  n - 3 .. n + 3 minus n), max-ops 3, 16 reps x 120 s per n (`runs/arms/a2`).  Pass = any layout holding n below the
  grid side k (need < k - f); stronger = finished side <= the pre-#575 record.  Prediction: pass 25 %.  Estimate 1.1 CPU-h.
* **A2 result (12:30)**: miss (prediction pass 25 %).  2,541 evals at 211 and 2,239 at 241 from n - 3 .. n + 3 arms with
  the fixed-frame fix and 0-3 local operators: the best count at f + 0.03 is n - 1 (210: 196 times, 240: 173), never n.
  Local +-1 edits of real arms do not find the extra square; consecutive records differ by collective rearrangement
  (columns redistributed between arms), which is a continuous move.  Arms + exact fill reproduce records but do not
  improve them.  Parked with the V-band constructor until a continuous re-hardening step exists (Evan 10-13: a
  dedicated performance session, GPU-style soft-to-rigid at scale, is planned; no CPU brute force now).
* **Day tally (12:35)**: ~10.7 CPU-h measured (F0 3.4, B0 0.4, W2 2.7, G1 1.1, A1 1.6, A2 1.3, interactive ~0.2), est. 10.9.
* **Population check (Evan 10-13: is 44,000 the right baseline?)**: no.  At n = 110 the store has 43,807 sub-11 rows
  (copies), 3,610 distinct f64 sides (unconverged copies spread up to ~1e-3 above their minimum, and not-minima), and
  **890 certified exact local minima** (+ 49 weak, 3,277 screened sides shown not-min).  The coverage yardstick uses the
  890; each has a store packing within 1e-8 of its exact side.
* **H1 holdout rebuild, pre-registered (12:55)**.  `arms.py holdout --n 110 --k 160 --budget 60` (`runs/arms/h1`): hold out
  160 of the 890 certified minima, arm library from the other 730 (8 symmetries), 60 s per held-out minimum on one core
  (pairs + 0-3 operators + exact fill, stop at the first layout holding 110 at side <= s_h + 1e-6, then fq polish).
  Metrics: rebuilt fraction vs time; same minimum (polished side within 1e-8 of s_h) vs a different one at or below s_h.
  Purpose: a reusable generative-coverage yardstick (arms now; explorer / GPU soft-to-rigid later on the same split).
  Prediction: rebuilt 30 % (wide: 10-70 %), mostly different minima at or below s_h rather than the same one.
  Estimate 2.7 CPU-h.  (First launch stopped after ~2 min: its "rebuilt" counted lower training minima; now a hit =
  a layout fitting in (s_h - 2e-8, s_h + 1e-3] whose polish lands on s_h within 1e-8.)
* **H1 result (13:20)**: **0 / 160** held-out certified minima rebuilt (prediction 30 %: miss).  1,693 candidates fitted
  within (s_h - 2e-8, s_h + 1e-3] and were polished; every one landed on a different minimum (the 890 minima sit ~3.6e-6
  apart on average in [10.9968, 11]).  Recombining real arms from other minima + exact fill + polish reaches minima
  already in the training set, not specific new ones.  As a yardstick "rebuild this exact minimum" is very strict; the
  record-relevant version is "hold out the best K, train on the rest, does the mechanism find anything below the best
  training side" (not run; to Evan).  Day tally ~14.2 CPU-h.
* **Is depth regional? (Evan 10-13; `landscape_depth.py`, analysis only, ~0 CPU)**.  n = 110, 890 certified minima,
  fingerprint distance.  (1) Landscape: Spearman(side, mean side of the k nearest minima) = 0.975 / 0.956 / 0.847 at
  k = 1 / 5 / 20 (shuffled ~0, 97.5 % < 0.09).  The top 50 minima are concentrated: mean pair distance 0.022 vs
  0.061 for random 50 (p < 0.002), and every top-50 minimum's nearest neighbour is top-50 (random: 6 %).  **Corrected
  (Evan, from the t-SNE map): they are two separate basins**, A = 11 minima 10.9967834-10.9971977 (record) and B = 39
  minima 10.9968870-10.9975201, closest pair 0.049 apart (25th percentile of all pairs; nearest neighbours inside a basin
  ~0.002); the split is stable for linkage radius 0.01-0.03.  (2) Chain
  (220,664 census v2-v4 proposals with parent links): Spearman(parent side, child side) = 0.80.  First-in-run reaches
  of top-10 / top-50 minima per 1000 proposals by parent depth decile: deepest 11.3 / 44.0; deciles 2-8 0-1.1 /
  0.2-2.0; shallowest decile 1.9 / 6.7 (U-shape: far / near-grid parents are the second-best source, as C13).  Revisits
  dominate raw counts (deepest decile: 559 / 1000 revisits of top-100).  (3) Other n (41 n with >= 1000 sub-k
  proposals in hunt runs): median Spearman(parent, child) 0.84 (IQR 0.76-0.88); the pooled decile curve is confounded
  by tied parent sides and revisits (not used).
  Reading: draws are not random; deep minima sit in one region and depth is inherited along the chain.  Re-finding
  "the best 10" is easy whenever their neighbours are visible; the hard step is entering the deep region (C05 far starts:
  record not reached in 2 x 1 CPU-h; s110-landscape: "the hard part is choosing the funnel").
* **Basin holdout sizes (for a region benchmark)**: minima within r of basin A: 11 / 31 / 69 / 160 at r = 0.01 / 0.02 /
  0.03 / 0.04, none of them in B; around B: 40 / 89 / 160 / 279, none in A.  Withholding A + its r = 0.03 neighbourhood
  leaves B's 10.9968870 as the best visible: "find anything below the best visible" = find a new basin (the record-hunt
  analog at a fixed, well-mapped n).
* F0 leftover stopped (11:50): the n = 444 certification (a composition equal to our posted 364 + L) ran 98 min on one
  core without finishing; not needed.  (Stopping it with `pgrep -f` also killed my own shell: the pattern matched the
  shell's command line, the failure pk/jobs.py was written to avoid.  Use pk/jobs.py stop / explicit pids.)
* **E0 basin-entry benchmark, pre-registered (11:55)**.  `bench110.py run --arms rand,above --reps 8 --minutes 60`
  (`runs/bench110/e0`).  Holdout `refsets/h110_basinA` (snapshot census110_1013: 939 minima; withheld = basin A, the
  record's single-linkage cluster at 0.02 among the top 50, + every minimum within fingerprint distance 0.03: 78;
  visible 861, best visible = basin B's 10.9968869868).  Recipe = census v4's explorer (adapt, legacy reward,
  sparse-share 0.5, fp bank = visible, novel-ref = visible), 1 proc per replicate, 16 concurrent.  Arms: rand = 40
  random visible minima; above = 40 random n = 110 packings with 11 < s < 11.05 (grid crossing).
  Metric: entry = first proposal with a sub-11 side below 10.9968869868 (the new-basin event); CPU to entry, entering
  proposal kind; secondary: per-kind sub-11 / new-minimum rates.  Predictions: rand enters in 4 / 8 within 60 min
  (median ~30 CPU-min among entries; basin B visible may hold the chain); above enters in 1 / 8.  Baseline only: no
  adoption rule; it calibrates the benchmark (rand 8 / 8 in < 10 min = too easy -> widen the radius to 0.04; 0 / 8 =
  lengthen).  Estimate 16 CPU-h.
* **Proposal kinds, observational (census v2-v4 logs, 220 k proposals; 0 CPU)**: first-in-run reaches of top-50 minima
  per proposal CPU-h: kicksym 10.8, bigkick 9.5, cross2 9.4, kick 9.2, ashallow 6.6, band 5.8, aswap 5.8, lkick 5.4,
  rowslide 5.0, crot 4.4, mirror 1.0 (melt 17.5 on only 306).  New minima per CPU-h 380 (mirror) to 1,220 (cross2).
  Confounded by the adaptive kind choice (bandit per parent); ablations on the E-benchmark decide.
* **E0 result (12:57)**: **16 / 16 replicates entered basin A** (predictions 4 / 8 rand, 1 / 8 above: both missed, the
  benchmark is easier than I thought).  CPU s to entry: rand 11, 114, 541, 588, 695, 699, 1375, 1443 (median 640);
  above 108, 113, 508, 700, 707, 713, 1267, 1925 (median 704).  Every replicate also reached the record itself
  (10.9967833966) within its hour.  Entries per CPU-h 1.01 in both arms: starting above the grid costs nothing at 110.
  Entering kinds spread over 8 of 11 (cross2 3, kicksym 3, bigkick 3, rowslide 2, aswap / kick / lkick / ashallow /
  band 1 each; crot, mirror 0).  Per-kind new minima per CPU-h in E0 match the census table (mirror lowest, 382).
  Calibration (rule: 8 / 8 in < 10 min = too easy): not all < 10 min (max 24-32 CPU-min), so radius 0.03 is usable:
  time-to-entry spreads over 2 orders of magnitude, a mechanism comparison resolves rate ratios of ~1.5-2x with 16
  replicates per arm.  Reading: at 110 the explorer's distribution reaches a withheld deep basin from random visible
  minima or from above the grid in ~10 CPU-min; no single proposal kind carries it.  It does not test other families
  (Evan: overfitting to the double-band V family, not to s(110), is the hope; a test set from another family is next).
  Day tally ~30.2 CPU-h.
* **Race benchmarks (`pk/benchkit.py`, reusable; `bench110.py race`)**: slots run replicates back to back; a replicate
  ends at its first event (entry below the best visible) or a 30 CPU-min cap and a fresh one (new seed, new random 40
  starts) takes the slot.  Rate = events / CPU (exponential, censored); arm vs control = rate ratio, 95 % interval from
  sqrt(1/e1 + 1/e0) (+0.5 each); sequential rule: at every check (each arm another --check CPU-h) stop when every arm's
  interval excludes 1 (faster / slower) or lies inside [1/1.5, 1.5] (same); else at the budget.  CPU per replicate from
  /proc (process group: utime + stime + reaped children).  Smoke (0.1 CPU-h) clean.
* **R1, pre-registered (13:35)**: `bench110.py race --arms control:rand nomirror:rand:--weights=mirror=0
  notop3:rand:--weights=cross2=0,kicksym=0,bigkick=0 finish3:rand:--finish=3 --budget 16 --check 1 --cap-min 30`
  (`runs/bench110/r1`; holdout h110_basinA; E0's recipe).  Questions: (1) is mirror (the lowest new-per-CPU kind)
  dead weight; (2) do the three kinds that entered most in E0 carry entry; (3) does the inline kick + polish-only finish
  (C11: cleaner, not more productive in census v2) change the entry rate.  Predictions: nomirror same (ratio 1.0-1.3);
  notop3 slower (~0.6); finish3 same.  Each arm ~4 CPU-h -> ~20 events: resolves ~2x.  Rule: as the sequential rule;
  an arm resolved "faster" is adopted into the census / hunt recipe only after a second race confirms it.
* **R1 result (14:05)**: 18.6 CPU-h, 97 replicates, 72 entries.  Rates (entries / CPU-h): control 3.88 (19 / 22 reps),
  finish3 5.13 (24 / 27), nomirror 4.41 (20 / 25), notop3 2.01 (9 / 23).  Ratios vs control: finish3 1.32 [0.73, 2.39],
  nomirror 1.14 [0.61, 2.11], notop3 0.53 [0.25, 1.16]: **all open** (predictions: nomirror same, notop3 ~0.6 slower,
  finish3 same: consistent, unresolved).  Entering kinds spread over 8-10 kinds per arm; with cross2 / kicksym / bigkick
  removed, rowslide / ashallow / aswap take over.  Reading: no single kind is necessary or dead weight at a 2x resolution;
  the top-3 drop probably slows entry (~0.5x) but did not resolve; the inline polish-only finish is not worse.
  **Two flaws found and fixed in benchkit**: (1) interim checks at 95 % gave false "faster" (finish3 3.47x, nomirror 3.27x
  at 1 CPU-h, both ~1.2x at the end): interim checks now use 99.5 %; the run continued only because notop3 was open;
  (2) the budget ignored running replicates (16 -> 18.6 CPU-h): it now counts them.  Day tally ~49 CPU-h.
* **Session closed (14:30, Evan: pause, write up, performance session next)**: summary and next steps in
  `HANDOFF-2026-10-13.md`; claims C25-C36 in `CLAIMS.md`; TODO "Packings" points at the handoff.  Nothing running.

## Performance session 10-13 (Evan: bit-perfect speedups, then minor-output-effect ones; no algorithm changes)
* **Profile (n = 110, bench110 recipe, 16 x 6 min, 4,499 proposals, ~1.34 CPU-s each):** fq ~89 % of CPU (Clarabel LPs
  ~60 % of the total, ALM screen ~13 %, other fq ~8 %), Python + per-call fq spawn + parent ~11 %.  `perf stat`: IPC ~4.0,
  L1d miss 5 %, nearly no misses beyond L2: fq is compute-bound at 110 (the KKT system fits in L2), not memory-bound.
* **target-cpu=native**: 2000 / 2000 replayed calls bit-identical, ~1.06x.  **Clarabel iterative refinement off**: polish
  1.64x, screen 1.29x on the same 2000 calls, not bit-identical (97 % of polished sides within 1e-9, symmetric tails
  +-1e-4).  exactsolve n = 237: 45 % in probe_flat's per-step dense lstsq; gmpy2 backend 1.11x only (and changes digits
  near 1e-73).
* **Q0, pre-registered (fq replay equivalence; a cheaper acceptance test than the 16 CPU-h race for fq-local changes):**
  Corpus C110 = the 8,892 fq calls (screen + polish, inputs + args) logged from the 16 x 6 min explorer run above
  (`runs/fqcorpus/c110`); fixed seeded sample of 4,000.  Arms replayed on the same calls, interleaved in one 16-proc
  pool: `base` (current source, native build), `j1`, `j2` (base with input jitter x, y +-1e-13, theta +-1e-11 deg, two
  seeds: the null = "same search, re-rolled"), `noir` (iterative refinement off).  Per call: side s, feasibility, child
  CPU (rusage), class (screen: discard s > k + 0.051, grid = full_lines, else kept; polish: grid / sub-k / below
  BEST_VISIBLE 10.9968869868).  Second corpus C237 (16 x 5 min explorer from n = 237 store packings, 400 calls; arms
  base, j1, noir) because polish dominates there.  Certification: 60 sub-k line-free polish outputs per arm (same calls)
  through exactsolve; certified fraction and |S_exact - s_fq|.
  Metrics vs base: (M1) divergence = fraction |ds| > SAME_SIDE (2e-9); (M2) bias = sign test on diverging calls + class
  frequencies with 95 % intervals; (M3) failures; (M4) certified fraction, precision; (M5) CPU ratio per stage.
  **Rule: adopt noir without a race if** M2 shows no class frequency shift > 1 pp (95 % CI inside +-1 pp for sub-k and
  below-best; sign test p > 0.01), M3 failures <= base + 0.2 pp, M4 certified fraction within 10 pp and median precision
  not 10x worse, M5 >= 1.3x fq CPU; any failure -> investigate, then the race.  Divergence alone (M1 > jitter's) is not a
  failure as long as it is unbiased.  Predictions: j1 diverges from base on ~30 % of screens, ~3 % of polishes; noir
  similar or larger, unbiased; speed polish ~1.5x, screen ~1.25x (C237: polish ~1.6x).  Budget ~6 CPU-h.
* **Q0 result (16:05):** ~5 CPU-h total (corpora + replays + certification).  **Adopted: refinement off is the fq default
  (`FQ_LP_IR=1` restores it); native build via `.cargo/config.toml`.**  All rule conditions hold on both corpora:
  | corpus / stage | diverge (j1 / j2 / noir) | class shifts (noir) | noir lower / higher (sign p) | noir speed |
  |---|---|---|---|---|
  | C110 screen (2010) | 0.76 / 0.76 / 0.29 | 0 of 2010 calls | 302 / 287 (0.56) | 1.21x |
  | C110 polish (1990) | 0.032 / 0.030 / 0.030 | 0 | 28 / 32 (0.70) | 1.54x |
  | C237 screen (234) | 0.82 / - / 0.51 | 0 | 67 / 53 (0.24) | 1.14x |
  | C237 polish (166) | 0.120 / - / 0.114 | 0 | 11 / 8 (0.65) | 1.55x |
  Failures 0 in every arm.  noir's |ds| quantiles are below jitter's everywhere (C110 polish q99 2.4e-8 vs 6e-8).
  Certification (M4, C110, 60 sub-k polish outputs, 45 of them diverging): 19 / 60 both arms (paired: 17 both, 39
  neither, 2 / 2 discordant); |S_exact - s_fq| median 2.4e-12 both.  Deployed binary checked: refinement-off default ==
  the tested noir arm, `FQ_LP_IR=1` == base, 300 / 300 calls exact (side + fingerprint).  Predictions: speed right
  (polish 1.5x / screen 1.25x predicted); noir diverged less than jitter on screens (predicted "similar or larger").
  Expected explorer gain: ~1.35x per proposal at 110 (fq ~89 % of CPU), ~1.45x at 237 (polish-dominated).  Side
  observation (CONFUSIONS X18): a 1e-13 jitter lands polishes lower more often than higher (C110 42 / 21, p = 0.011).
  Benchmark for future fq-local changes: `fqreplay.py` (record / run / compare / certify), corpora `runs/fqcorpus/c110`,
  `c237` (gitignored), ~2 CPU-h per arm set; the 16 CPU-h race stays for search-level changes.
* **Perf score (16:30; `perfbench.py`, new)**: explorer proposals per measured CPU-h (Python included), 8 single-proc
  replicates each at n = 110 (bench110 recipe and start pool) and n = 237 (4 of the 200 best store packings per
  replicate), fixed seeds and frozen starts (`runs/perfbench/starts*`), 6 min; score = geometric mean of the two
  speedups vs `runs/perfbench/orig` (commit 6adb4f3 binaries: generic build, refinement on).  Counts only with an
  fqreplay equivalence pass.  **orig 2,792 / 764 -> current (efc88ff) 3,983 / 997 proposals per CPU-h: 1.43x / 1.31x,
  perf score 1.36.**
* **Q1 faer LDL backend (Clarabel `faer-sparse`, supernodal): negative.**  fqreplay C110 1000 + C237 200 calls: speed
  0.98-1.04x, outputs re-roll-equivalent.  The KKT factor is too sparse for supernodes.
* **LP census after Q0 (FQ_LPT with tags; 60 C110 + 20 C237 calls):** LPs are 66 % of fq at 110 (step 28 %, SOC 20 %,
  flip 16 %, final 2 %), 87 % at 237 (step 43 %, SOC 39 %, flip 4 %); 18-25 IPM iterations per LP; numeric LDL ~42 %
  of fq (qdldl).
* **Q2, pre-registered (active-set SOC; a minor algorithm change, so not bit-identical by design):** `FQ_SOC_AS=1`:
  the SOC LP (same matrix and objective as the step LP, shifted rhs) is replaced by the min-norm correction on the step
  LP's active set (rows with slack < 1e-7, box-active columns and tau = 0 fixed), one quasidefinite LDL solve (Clarabel's
  qdldl); if the corrected step violates any row / bound, fall back to the Clarabel SOC LP; the merit test is unchanged.
  Test = fqreplay on C110 (2000) + C237 (400), arms cur, j1, socas; Q0's rule (no class shift > 1 pp, unbiased sign
  test p > 0.01, failures <= base + 0.2 pp) plus certification of 60 polish outputs (within 10 pp) and median SLP
  iterations per polish not up by > 20 %.  Predictions: fallback < 20 %, fq speed C237 polish ~1.5x, C110 polish
  ~1.2x, screens ~1.05x; diverges more than Q0's noir (different SOC points) but unbiased.  Budget ~3 CPU-h.
* **Q2 result (17:40): faster, and not a re-roll: biased toward lower sides, i.e. better-converged polishes.**  fqreplay
  (`runs/fqreplay/q2_socas_c110`, `_c237`): speed screen 1.13x / 1.08x, **polish 1.47x / 1.62x** (C110 / C237); 0
  failures, 0 class shifts (all CIs contain 0); lower / higher C110 screen 233 / 133 (p < 0.001), polish 34 / 15 (p =
  0.009); C237 66 / 41, 23 / 9 (p = 0.02 each).  **Rule outcome: fails "sign test p > 0.01"** (favourable direction).
  Investigation (certification, 60 polish outputs per corpus, diverging calls first): certified C110 21 -> 29 (paired
  9 only-socas / 1 only-cur), C237 21 -> 28 (10 / 3); |S_exact - s_fq| identical where both certify (median 2.3e-12 /
  3.3e-12).  Reading: the active-set correction lets polishes go past stalls that the IPM SOC leaves (the kick_retry
  problem, X01): outputs are more often true jammed minima.  Fallback to the LP SOC: ~1 % (110), ~10 % (237) in a
  smoke test.  Per the rule, adoption needs the race (Evan's call: race vs adopt on this evidence).  Not deployed
  (`FQ_SOC_AS=1` to use it).
* **Q3, pre-registered (FQ_FASTDOT=1, 8-accumulator dot in the ALM L-BFGS; rounding-only change):** fqreplay C110
  1000 + C237 200, arms cur, j1, fastdot (+ the experiment binary in default mode, must be bit-identical to cur).  Q0's
  rule.  Prediction: re-roll-equivalent; screen speed ~1.05-1.1x, polish ~1.0x.
* **Q3 result (18:10): re-roll-equivalent, small.**  Screens C110 1.12x (507 calls; j1 0.98x), C237 0.94x (111 calls;
  j1 0.97x, and j1 read 1.28x in Q1: screen CPU per call at 237 is too noisy at this sample size); unbiased (p 0.39 /
  0.60), 0 class shifts; polishes bit-identical (no L-BFGS in polish); the experiment binary in default mode is
  bit-identical to deployed.  ~4 % of explorer CPU at 110.  Misses the 1.3x speed bar; kept as a knob
  (`FQ_FASTDOT=1`), not default.
* **Q2 race, pre-registered (18:20; Evan "keep going", the race was the pre-registered path for Q2):**
  `bench110.py race --arms control:rand socas:rand:+FQ_SOC_AS=1 --budget 16 --check 2 --slots 10
  --out runs/bench110/q2race` (bench110 recipe, entry events below BEST_VISIBLE; arms now take `+KEY=VAL` environment
  tokens, so fq knobs race without a separate binary).  10 slots, leaving ~6 cores for development.  Expected: ~4 entries
  / CPU-h control, so ~30 events per arm: resolves ~1.7x.  Prediction: socas 1.2-1.5x (fq ~1.3x faster per proposal at
  110, plus better-converged polishes), most likely "open".  **Rule: adopt FQ_SOC_AS=1 as the default unless the race
  resolves "slower"** (the speed gain is measured separately; the race guards against a search-level loss the replay
  can't see).
* **Where polish time goes with Q2 (18:15, dev measurements, Q2 on):** perf, 12 C237 polishes: qdldl numeric factor 49 %,
  triangular solves 16 %, Clarabel IPM other 8 %, per-LP setup (qdldl symbolic, triplet sorts) ~12 %.  LP rows ~2,800
  (+ 2 x ~1,400 box rows) x ~1,400 columns, ~18-24 IPM iterations; only ~0.1-1 % of rows are provably redundant (jammed:
  nearly every near row can activate).  Step LPs: ~40 % rejected, mostly the R shrink from r0 = 1e-3 at near-converged
  inputs.  Flip search (`FQ_FLIPDBG=1` census): C110 200 polishes: flip LPs 33 % of polish LP time, 65 % of them in
  searches that find nothing (114 of 332 stalls); winning pair ranks spread over 0-7 (16 % at rank >= 4, so a smaller
  --flip-top would lose real flips).  C237: flips 4 %, SOC-fallback LPs 16 % of LP time.  CPU split per fq stage with
  Q2: polish 59 % / screen 41 % (C110), 67 / 33 (C237).
* **Q4 (LP tolerances, dev measurement, negative):** `tol_gap_abs/rel` 1e-7 / 1e-6 (default 1e-8) cut IPM iterations by
  < 1 per LP and raised SOC fallbacks 2x / 5x (the active set is read from the LP solution).  Knob removed.
* **Q5 (flip LPs as decisions, dev measurement, negative):** a Clarabel termination callback stopping a flip LP once the
  primal (found) or dual (none) cost settles "pred > 1e-12 s": 2,761 / 2,761 decisions agree with the full solve, but it
  saves ~1 of 18 IPM iterations (failing flips sit at the threshold to within the final gap; residual bound 1e-7 or 1e-3
  makes no difference).  Removed.
* **Python share with Q2 (cProfile, in-process worker, 3 min each):** Python is 3.6 % (n = 110) / 4.9 % (n = 237) of
  explorer CPU; process spawning ~0.1 %, so `fq serve` is not worth wiring.  Half to most of the Python is
  `layout._straight_lines` (O(m^2) loop).  **Vectorised (bit-identical by construction, same float64 ops): 3,000 packings
  x 3 k, 0 mismatches, 16x faster** (~2-3 % of explorer CPU).  Old loop kept as `_straight_lines_ref`.
* **Q6 (adaptive trust-region shrink, dev measurement, negative):** on rejection shrink R by clamp(0.5 / (1 - rho), 4,
  cap) instead of 4: cap 16 / 64 / 256 all *raise* step LPs (C237 dev 40 polishes: 3,167 -> 3,269-3,313).  Long polishes
  carry the step LPs (n = 237: 72 % of SLP iterations are in calls with >= 100 iterations); only 13 % (110) / 4 % (237)
  of iterations come after s is within 1e-12 of its final value, so neither R0 nor termination is a lever.  Removed.
* **Dev hygiene from here on:** tuning uses the already-used slice of fqreplay's seed-1 shuffle (C110 calls 0-3999, C237
  0-399); confirmation runs use the untouched rest (`fqreplay.py run --skip`), so each knob is confirmed once on data it
  was not tuned on.
* **Q7, pre-registered (19:10; own LP interior point, `FQ_LP_IPM=1`; minor algorithm change: same LPs, different solver):**
  `src/lpipm.rs`: Mehrotra predictor-corrector, bounds native (no box rows), normal equations A'(Y/S)A + bound diagonal
  (~1,450 columns at 237, vs Clarabel's ~7,300-row KKT), qdldl with AMD per LP, static primal 1e-8 and dual 1e-8
  regularisation (Clarabel's values, in normal-equation form), tolerances gap 1e-8 abs / feasibility 1e-8, start primal
  feasible (w = 0, tau above every row), Clarabel fallback on breakdown (accepts the iterate if within 10x of
  tolerance).  Dev (dev slices): per LP 13 ms vs 23 ms at 237 (1.75x; factor 55 % of the rest), objectives agree within
  the 1e-8 gap tolerance (optimal faces are large: w differs O(1) between the two solvers' stopping points, as between
  Clarabel runs at different settings), fallback ~1 in 3,700 LPs; polish wall 1.37x (110) / 1.9x (237) on dev batches.
  Built on Q2 (all arms FQ_SOC_AS=1).  Frozen binary `runs/fqreplay/q7bin/fq` (sha256 d1abd440...); with FQ_LP_IPM unset it
  is bit-identical to deployed (100 / 100 calls).
  **Test** (held-out): fqreplay C110 `--skip 4000 --sample 2000`, C237 `--skip 400 --sample 400`; arms cur
  (FQ_SOC_AS=1), j1 (cur + jitter, the null), ipm (FQ_SOC_AS=1 FQ_LP_IPM=1); certification 60 sub-k polish outputs per
  corpus (cur vs ipm, diverging first).  **Rule: adopt FQ_LP_IPM=1 as default** (together with Q2, which it is built on)
  if: every class-frequency difference vs cur has a 95 % CI inside +-1 pp; failures <= cur + 0.2 pp; sign test p > 0.01
  in each of the 4 corpus x stage cells (4 tests: family-wise ~4 %) -- a favourable-direction failure goes to the
  certification as for Q2; certified fraction within 10 pp of cur and median |S_exact - s_fq| not 10x worse; polish fq
  CPU >= 1.2x on both corpora.  Predictions: polish 1.35x (110) / 1.7x (237), screens ~1.1x; divergence at least
  jitter's; unbiased.  Budget ~1.5 CPU-h replay + ~1 CPU-h certification, at 6 procs beside the Q2 race (arms are
  interleaved per call, so the ratios stay paired).
* **Q8 (Gondzio multiple-centrality correctors in lpipm, `FQ_IPM_GONDZIO=K`, dev measurement, neutral):** IPM
  iterations 17.8 -> 15.9 (K = 1) -> 14.9 (K = 2) at 110, 19.8 -> 16.4 at 237, but LP time at 110 rises 3-10 % (the
  factor is cheap there, so the extra solve + vector passes cost about an iteration); ~15-20 % faster at 237 on 2 calls
  (loaded machine).  Not worth a test slot; kept as a knob, off.
* **Bit-identical micro-opts (dev, 100 / 100 calls identical with IPM on and off):** eval reads the container / rounding
  atomics once per evaluation; L-BFGS reuses its history buffers.  Speed 1.00-1.02x (noise): free, negligible.
* **Q7 result, C110 held-out (2000 calls):** screens 1.15x, **polish 1.66x** (pred. 1.35x); 0 failures; class shifts all
  within +-0.3 pp; sign tests p = 0.59 (screen, 165 / 176) and 1.0 (polish, 26 / 25); divergence below jitter's (screen
  0.33 vs 0.75, polish 0.053 vs 0.040, |ds| q99 4e-5 vs 9e-1 / 3e-5 vs 1e-7).
* **Q2 race result (20:05): open, so per the rule Q2 is adopted (active-set SOC is the fq default; `FQ_SOC_AS=0`
  restores the SOC LP).**  16.0 CPU-h, 172 replicates, 159 entries: control 68 / 7.96 CPU-h = 8.55 per CPU-h, socas
  91 / 8.05 = 11.31; ratio **1.32 [0.97, 1.81]**.  Entering kinds spread as in R1.  Predictions: 1.2-1.5x, "open" --
  right; the event rate was 2.5x the pre-registered guess (9.9 vs ~4 per CPU-h), so the run resolved ~1.45x instead of
  ~1.7x.
* **Was 16 CPU-h the right amount? (Evan asked.)**  The race had one job: catch a search-level loss that the replay
  cannot see.  Power for that (rate ratio, one-sided 5 %, 80 % power, ~10 events / CPU-h): a true 0.5x is caught with
  ~5 CPU-h, 0.67x with ~15 CPU-h, 0.8x needs ~50 CPU-h, 0.9x ~220 CPU-h.  This run could catch 0.67x and worse.  What
  we needed: Q2 had already passed fqreplay (0 class shifts on 2,400 paired calls) and certification (more outputs
  certified), and the mechanism (fewer polish stalls) points toward better convergence, so the plausible failures were
  either nothing or something gross (e.g. the explorer starving some move kind); ~5 CPU-h (catches 2x) would have
  served the decision equally well.  Moderate harms (0.8-0.9x) are out of reach of any race we would run: for
  fq-local changes the paired replay is the strong instrument (most calls are identical or differ by 1e-10, so the
  noise is far below the race's per-replicate exponential noise).  **Policy from here:** fq-local changes that pass
  fqreplay + certification get a ~5 CPU-h guard race (stop early only on "slower"), not 16; search-level changes keep
  races sized from a stated minimal effect of concern (~50 CPU-h for 0.8x), and if that is too expensive we say the
  test can only catch large effects.
* **Q7 result (20:05): adopted (`FQ_LP_IPM=0` restores Clarabel).**  C237 held-out (400 calls): screens 1.14x, **polish
  1.67x**; 0 failures; 0 class shifts; sign p = 0.069 (screen 50 / 71) and 0.20 (polish 19 / 11).  Certification: C110
  cur 30 / 60 vs ipm 35 / 60, C237 32 vs 28 (both within 10 pp), median |S_exact - s_fq| identical (2.3e-12 / 3.3e-12).
  Every rule condition holds.  Deployed default == the tested arm (150 / 150 calls, side + fingerprint).
* **Perf score after Q2 + Q7 (20:15, quiet machine): 2.77x** (`runs/perfbench/q2q7`): 7,366 / 2,222 proposals per CPU-h
  at 110 / 237 vs 2,792 / 764 pre-session (2.64x / 2.91x); was 1.36x after Q0.  Consistent with the stage speedups
  (polish ~2.4-2.7x, screens ~1.3x).  CLAIMS C38-C40.  Next costs at 110: screens are now ~75 % ALM (L-BFGS two-loop
  dot latency 37 %, eval 25 %, sincos 11 %), so FQ_FASTDOT (Q3, equivalence already shown) is worth re-measuring.
* **Q9, pre-registered (fast dot on the Q2 + Q7 base; rounding-only change):** `FQ_FASTDOT=1` (8-accumulator dot in the
  ALM's L-BFGS, Q3).  Q3 showed it re-roll-equivalent on the old base; it missed Q3's 1.3x bar when LPs dominated, and
  screens at 110 are now ~75 % ALM (L-BFGS dot latency ~37 %).  No race: summation order is exactly what the jitter
  null models, and Q3 already passed equivalence.  Test (fresh held-out slices): fqreplay C110 `--skip 6000 --sample
  1500`, C237 `--skip 800 --sample 300`; arms cur (deployed), j1 (jitter), fd.  **Rule: adopt if** Q0's equivalence
  conditions hold (class CIs inside +-1 pp, failures <= cur + 0.2 pp, sign p > 0.01 per corpus x stage), polishes are
  bit-identical (no L-BFGS there), and screen fq CPU >= 1.05x at C110.  Predictions: screens 1.25x (110) / 1.1x (237),
  polish 1.00x, unbiased.  Budget ~0.5 CPU-h.
* **Q9 result, part 1 (20:40):** C110 (1,500 calls): screens **1.16x**, polishes bit-identical (730 / 730), 0 failures,
  unbiased (p 0.32), all class CIs inside +-1 pp.  C237 (300 calls): screens 1.15x, polishes bit-identical, unbiased
  (p 0.66), but **above-k -2.6 pp [-5.1, -0.1]: the +-1 pp condition fails** -- and the jitter null fails it too (its
  CIs span +-2.5-4 pp at 153 screens): the criterion was unattainable at that sample size.  Not adopted; Evan chose
  the higher-power test.
* **Q9b, pre-registered (rule fixed before the data):** new corpus **C237b**: same recipe and starts as C237
  (`--adapt --reward legacy --sparse-share 0.5`, 16 procs), seed 2371, 10 min, deployed fq (Q2 + Q7); expected ~4,000+
  calls, ~2,000-2,500 screens.  Replay **all screens** (`--stage screen`), arms cur, j1, fd.  Power: per-screen class
  variance for a re-roll, from the 153-screen run, gives 95 % CI half-widths of ~0.6-1.0 pp at 2,500 screens (grid the
  widest); a fixed +-1 pp margin would need ~5,000.  **Rule: adopt FQ_FASTDOT=1 if** every class difference fd vs cur
  has a 95 % CI inside **+-1.5 pp**, sign p > 0.01, failures <= cur + 0.2 pp, screen CPU >= 1.05x; **calibration:** j1
  must pass the same class condition, else the test is underpowered (inconclusive, not a failure).  Prediction: pass,
  screens ~1.15x.  Budget ~2.7 CPU-h recording + ~1.6 CPU-h replay.
  **Methodology note for future tests:** class-shift margins are set per corpus from the re-roll variance (and checked
  by the jitter arm passing), not a fixed +-1 pp.
* **Q9b result (21:20): passes; FQ_FASTDOT is the fq default (`FQ_FASTDOT=0` restores the sequential dot).**  C237b
  (`runs/fqcorpus/c237b`, 7,117 calls, 3,888 screens, all replayed): screens **1.21x**; 0 failures; sign p = 0.92
  (1,648 / 1,655); class differences fd vs cur: above-k -0.13 pp [-0.59, +0.34], discard +0.41 [-0.13, +0.95], grid
  -0.21 [-0.92, +0.51], sub-k -0.08 [-0.50, +0.34] -- all inside +-1 pp, well inside the pre-registered +-1.5 pp;
  calibration: j1 passes too (its CIs <= 0.7 pp).  Deployed default == the tested fd arm (150 / 150).  The Q9 C237
  "above-k -2.6 pp" was noise at 153 screens, as the power calculation said.
* **Perf score after Q9: 2.92x** (`runs/perfbench/q2q7q9`): 7,896 / 2,298 proposals per CPU-h at 110 / 237 (2.83x /
  3.01x).  Fast dot added 1.07x / 1.03x (predicted ~1.06-1.08 at 110).
* **Next-round profile (21:00, deployed Q2 + Q7 + Q9):** n = 110: qdldl factor 19 %, eval 18 %, L-BFGS 17 %, rest of
  the IPM ~22 %, sincos 9 %; n = 237: IPM ~80 % (factor 36 %).  Factor: 0.6-0.9 M multiply-adds in 0.21-0.29 ms
  (~3 GFlop/s; L ~27 k, column counts 5-70); faer simplicial / supernodal (relaxed amalgamation up to 256) is never
  faster than qdldl on a quiet machine (0.23-0.26 vs 0.21 ms) -- supernodes too small; a custom block (3x3 per square)
  factor is the remaining big lever, deferred.  IPM sections (FQ_IPM_PROF=1, 237): factor 47 %, setup 10 %, assembly
  10 %, step lengths 7 %, solves 7 %, residuals 5 %, A'y 5 %, Ax 4 %.
* **Dev negatives (no test needed):** gather-form A'y and A'DA assembly (bit-identical, 0.99-1.00x: reverted); IPM
  start mu0 0.1 - 1e-4 (iterations 20 -> 21-37: reverted); faer relaxed supernodal (above).
* **Bit-identical ALM speedups (deployed; 800 / 800 calls identical vs the previous deploy, both corpora, both stages):**
  (1) `term()` skips its two divisions when lam = 0 (inactive constraints; the general formula reduces exactly,
  signed zeros included); (2) the L-BFGS two-loop fuses each update of d with the next fast dot (same 8-lane order).
  Screens 1.17x (110) / 1.22x (237), polishes 1.00x.  **Perf score 3.06x** (`runs/perfbench/bitid1`: 8,377 / 2,384
  proposals per CPU-h; load ~7 from another session at the start, so possibly a little low).
* **Q10, schedule-replay LDL (`src/ldl.rs`; bit-identical, deployed):** qdldl's numeric factor re-walks the elimination
  tree for every row on every refactor (~40 % of its time is that bookkeeping).  Within one LP the pattern is fixed, so
  the symbolic pass records qdldl's schedule (per row: eliminated columns in qdldl's visiting order, target slot in L)
  and the numeric pass replays it -- the same flops in the same order; same AMD call, qdldl's permuted-column entry order
  (input order -- sorted rows gave 1e-9 differences), qdldl's regularisation and solve.  Bench (n = 237 step LPs):
  bitwise-equal D and L, 1.28-1.37x per refactor.  In fq: 800 / 800 replayed calls identical (both corpora, both
  stages); polish 1.14x (110) / 1.18x (237), screens 1.04x / 1.06x.  **Perf score 3.38x** (`runs/perfbench/bitid2`:
  8,885 / 2,740 proposals per CPU-h; load ~9 from other sessions at the start, so possibly low).  Also measured: AMD
  is ~0.25 ms of the ~0.4 ms symbolic setup per LP (~2 % of fq) -- reusing orderings is not worth a re-roll.
* **Q11, pre-registered (vectorised sincos in the ALM eval, `FQ_FASTSC=1`; rounding-only change):** every square angle
  and separator angle per eval in two batches through a branch-free Cody-Waite + Cephes sincos (<= 2 ulp from libm
  over 6 M random angles; 8.1x faster than libm when vectorised).  Dev (C110 dev slice, 300 screens): screens 1.11x;
  knob off is bit-identical to deployed.  Polishes do not evaluate the ALM objective, so they should be bit-identical.
  **Fresh held-out corpora** (also for the next re-roll tests): C110c = explorer, bench110 recipe, perfbench starts110
  rep0 + rep1 (80 starts), seed 1101, 16 procs x 6 min; C237c = C237b's recipe, seed 2372, 6 min; deployed binary.
  Test: all screens of each, arms cur, j1, fsc (+ 200 polishes per corpus for the bit-identity check).  **Rule (Q9b
  form):** adopt if every class difference has a 95 % CI inside +-1.5 pp, sign p > 0.01, failures <= cur + 0.2 pp,
  polishes identical, screen CPU >= 1.05x at C110; calibration: j1 passes the same class condition (else
  inconclusive).  Prediction: pass; screens ~1.1x (110) / ~1.04x (237).  Budget ~3.2 CPU-h recording + ~1.5 replay.
* **Q11 result (22:05): adopted (vectorised sincos is the default; `FQ_FASTSC=0` = libm).**  C110c (3,000 screens):
  1.09x, class CIs within +-0.75 pp, sign p 0.92; C237c (all 2,862 screens): **1.17x**, class CIs within +-1.4 pp
  (grid +0.63 [-0.11, +1.37]), sign p 0.65; j1 passes the same margin on both; 0 failures; 400 / 400 polishes
  identical.  Corpora: C110c 25,863 calls (13,737 screens), C237c 5,287 (2,862 screens).
* **AMD reuse (dev, dropped):** consecutive LPs have the same column set only 23 % (237) / 51 % (110) of the time;
  saving <= ~2-4 % of fq.  Not worth a re-roll.
* **Q12, pre-registered (active-set SOC repair rounds 4 -> 16, `FQ_SOC_ROUNDS=16`; minor change inside Q2's SOC):**
  every SOC fallback in a 30-polish C237 dev sample was "rounds exhausted" (650 fallback LPs, 21 % of LP time,
  concentrated in long polishes); dev: 8 rounds -> 362 fallbacks, 16 -> 245, 32 -> 245; LP time 39.3 -> 34.6 s.
  Test (held-out polish slices, disjoint from Q11's): C110c polishes 0-1999, C237c polishes 0-1199; arms cur
  (deployed incl. Q11), j1, r16 (FQ_SOC_ROUNDS=16); certification 60 sub-k outputs per corpus (cur vs r16).
  **Rule (Q7 form, Q9b margins):** adopt if class CIs inside +-1.5 pp (j1 must pass too), failures <= cur + 0.2 pp,
  sign p > 0.01 or a favourable bias confirmed by certification (as Q2), certified fraction within 10 pp and median
  precision not 10x worse, polish CPU >= 1.05x at C237c.  Predictions: most polishes identical (only fallback cases
  change), unbiased or slightly favourable; polish 1.05-1.1x at 237, ~1.0x at 110.  Budget ~2 CPU-h + certification.
* **Q12 result (22:25): not adopted; kept as a knob (`FQ_SOC_ROUNDS`, default 4).**  C237c (1,200 polishes): 1.06x,
  0 class shifts, 0 failures, but **biased lower** (76 / 37, p < 0.001; j1 78 / 82); C110c (2,000): 0.99x, 0.8 %
  diverge, unbiased.  Certification: C110c 46 vs 46 / 60; C237c 10 (cur) vs 8 (r16) of 60 diverging outputs (paired 5
  only-cur / 3 only-r16), median precision equal.  The rule needs a favourable bias *confirmed* by certification; it is
  neither confirmed nor contradicted, and the gain is small (~4 % of fq at 237, ~0 at 110).  Parked.
* **Policy clarification:** rounding-only changes (Q3/Q9 fast dot, Q11 sincos) are pure re-rolls by construction: replay
  only.  Changes to the algorithm's path (Q2, Q7, Q12) get the ~5 CPU-h guard race from the Q2 analysis -- Q7 did not
  get one (it was adopted before that policy); the cumulative guard race below covers it.
* **G1, pre-registered (cumulative guard race for this session's fq changes):** `bench110.py race --arms
  pre:rand:+FQ_LP_IPM=0;+FQ_SOC_AS=0;+FQ_FASTDOT=0;+FQ_FASTSC=0;+FQ_LP_IR=1 cur:rand --budget 5 --check 1 --slots 16
  --out runs/bench110/g1` (bench110 recipe; entry events below BEST_VISIBLE; "pre" = the session's fq math switched off
  -- Clarabel LPs with refinement, SOC LP, sequential dot, libm sincos; native build and the bit-identical changes stay).
  Two questions: (1) guard: does any search-level loss offset the speed (power: ~5 CPU-h catches a true ratio <= ~0.5
  relative to the speed gain); (2) does the ~3x proposal throughput convert into ~3x entries per CPU-h?  **Rule:** if cur
  resolves slower than pre, investigate and consider reverting Q2 / Q7; otherwise the adoptions stand.  Prediction:
  ratio 2.5-3.5x (CI roughly x/1.6).
* **G1 result (22:55): cur / pre = 3.04 [1.45, 6.35] entries per CPU-h ("faster"); the adoptions stand.**  5 CPU-h:
  pre 9 entries / 2.56 CPU-h (3.5 per CPU-h), cur 27 / 2.44 (11.1).  (1) Guard: no search-level loss; (2) the ~3x
  proposal throughput converts to ~3x entries per CPU-h (prediction 2.5-3.5x: right).  Event rates were ~3x below my
  planning guess, so the interval is wide (x/2.1), still far from "slower".  Cross-check: the cur arm uses the same
  seeds and start draws per replicate index as the Q2 race's socas arm; median CPU to an entry 91 s vs 180 s there.
* **perfbench noise, and the Q11 "regression" at 237 (23:15):** after Q11 the perf score read 3.16 / 3.22 (two runs) vs
  3.38 before; n = 237 fell 2,740 -> 2,431 / 2,528 while `FQ_FASTSC=0` on the same binary reproduced 2,720.  Paired by
  replicate (same starts): 0.895 [0.82, 0.98] and 0.927 [0.84, 1.03] (the two "on" runs share trajectories, so they are
  not independent); n = 110: 0.99.  Per-status CPU: screens cheaper, polished proposals ~20 % dearer.  **Chain replay**
  (`chain.py`: screen with arm X, then polish its output with deployed defaults; C237c screens, + jitter null):
  1,200 chains: polish CPU on / off 1.017 [0.96, 1.08], jitter / off 1.050 [0.99, 1.11]; whole chain (screen + polish)
  on - off -0.9 % [-6.5, +4.8], jitter +5.3 % [-0.6, +11].  **No chain effect: any re-roll moves the following polish's
  cost by ~+-5 % (heavy-tailed per chain); Q11's expected ~+4-5 % at 237 is below that floor and below perfbench's
  resolution.**  Q11 stands.  Lessons: (1) perfbench resolves ~+-10 % at 237 per pair of runs (replicate SE 3-9 % per
  run, and a re-roll redraws all 8 trajectories at once); quote perf scores as ~3.2-3.4x, not to 0.01; (2) stage speed
  should come from paired replays, end-to-end gains below ~10 % from chain replays with >= 1,000 chains.

## Soft-to-rigid engine (10-14; Evan: build new features first, background compute yields to tests)
* **Why**: C29's three negatives used fq's rounding pre-stage, a quasi-static ALM quench at 10-40 radius steps.  Nakajima's
  protocol (paper §2, `sim.js` step rule; read as data, scratch copy) is ~10^4 momentum steps on a penalty energy:
  E = sum (2r - sd)_+^2 + walls^2 + mu L, beta 0.85, lr 0.04, |v| cap 0.03, Gaussian noise.  Overlaps tolerated
  transiently, inertia; a different operator.  Ostashev #575: "every square soft, box held, then release and compress".
* **Built**: `src/harden.rs` (`harden fresh | reharden | test`; Verlet lists, legalisation by centre scaling) and
  `harden.py` (inputs x transform none / insert:K / remove:all|R x parameter grid x seeds -> harden -> fq polish ->
  vs target).  `harden test`: pair gradient vs FD 0 / 19,999 off by > 1e-5 (worst 2.8e-9), energy gradient 0 / 90.
* **Validation (Nakajima's n = 11 numbers, standard 5k / 45k / 3k schedule, 160 seeds each, raw legalised side)**:
  anneal 9 / 160 reach Trump 3.8770845 (paper 6 / 160); no anneal 0 / 160; rigid control 0 / 160, best 3.88710
  (paper 0 / 120, best 3.8877).  0.27 s per run.  Implementation reproduces the paper.
* **Dev (241 from store 347, 12 cells x 4 seeds, `runs/harden/dev241a`)**: 2-10 s per run incl. fq polish.  Noise 0.001
  returns to the seed basin (+3e-12); noise 0.004 lands 1.3-1.8e-3 above.  0 / 48 below the seed.
* **T1 band-merge positive control, pre-registered (before launch)**: inputs store 347 (241, 15.9881324395; #575 ->
  15.987134249) and store 94 (211, 14.9979607050; #575 -> 14.993694187).  Grid r0 {0.05, 0.1, 0.2} x noise {0.002,
  0.004, 0.008} x hold {8000, 30000}, 24 seeds per cell per n (864 runs, est. ~1.7 CPU-h).  Success (per n): polished
  side <= 15.9876 (241) / < 14.9979607 - 1e-8 (211).  Prediction: 20 % for at least one n (dev: the move is not
  common; Ostashev's batch sizes are unknown).  If 0: look at the trajectories (do bands re-tilt at all?) before more
  volume; volume x10 is cheap (~17 CPU-h) only if the outcome distribution has a tail toward the target.
* **T2 insertion (border growth), pre-registered**: the 10-05 register's 307 (`../exact/batch/inputs/n-307.txt`,
  17.981030549) + 1 square at a uniform random pose per seed (`--xform insert:1`), box held at hold_f x 17.98103.
  Grid hold_f {1.0, 1.001} x r0 {0.1, 0.2} x noise {0.002, 0.004}, 24 seeds (192 runs, ~0.5 CPU-h).  Baseline: H1
  (10-12, chain insertion, ~8 CPU-h): best 18.0005167, never below 18; Fang's #484 reached 17.99827 from 307.
  Success: any polished 308 < 18 - 1e-8 (then candidate.py).  Prediction: 15 %.
* **T3 plateau rebuild, pre-registered**: Göbel's 233 (store 102, 15.778174593) and 296 (store 128, 17.707106781) with
  every single square removed (`--xform remove:all`), 1 seed.  Cells: control (r0 0.001, noise 0, hold 2000 = polish
  of the removal) and r0 {0.05, 0.15} x noise 0.004 x hold_f 0.999 (box squeezed 0.1 %).  ~1,600 runs, ~3 CPU-h.
  The pre-break values are the plateaus (s(232) = s(233), s(295) = s(296)); itsnaka #481 broke 232 (15.7674283,
  later 15.7641743 Ostashev / Mishapolk), ry-xu #432 broke 295 (17.7042328).  Success (per n): polished side below the
  plateau by > 1e-8; strong: at or below the first break.  Prediction: any break 30 % (the control is the C31 /
  removal_ctl kind of move, expect 0), strong 5 %.
* **T1 result: pass at 211, fail at 241** (864 runs, 0 failures, ~25 min wall on 16 procs).  211: 11 / 432 below the
  seed, best polished **14.9936458090** (seed 14.9979607050, -4.3e-3; Ostashev #575 14.993694187; today's frontier
  14.9935589140 = #583 / #593 refinements of the #575 packing), all hits at r0 0.2 (cells noise 0.004 / 0.008, hold 30k
  and noise 0.008 hold 8k).  241: 0 / 432 (best = the seed itself; r0 0.2 cells all +5e-3 or worse).  285 / 269
  distinct polished sides per n: the operator does move between basins.  Prediction (20 % for one n): hit.
  Reading: the engine does #575's kind of move at 211; softness ~0.2 is where it happens.
* **Amendment to T2 / T3 (before their launch, after T1)**: T1's hits were all at r0 0.2, so both get a softer cell.
  T2 grid: hold_f {1.0, 1.001} x r0 {0.1, 0.2, 0.3} x noise {0.002, 0.004} (288 runs).  T3 cells: control, r0 0.15
  noise 0.004 hold_f 0.999, r0 0.25 noise 0.004 hold_f 0.999 (default hold 8000; launched so -- the noise 0.006 / hold 30k I first
  wrote here was not what ran) (separate dirs: t3ctl, t3).  Predictions unchanged.
* **T1 structure check**: seed 94 has 74 tilted squares (one connected group), the hit 14.9936458 has 81 in one band,
  the current record (#575 refined) 82 in one band, with nearly the same square indices: the hit is in the record's
  V family, reached from the pre-#575 seed.
* **T1b second round (Ostashev: "a second round was seeded from that V"), pre-registered**: the 5 best distinct T1 211
  outputs (14.99365 .. 14.99464) re-hardened, cells r0 {0.05, 0.1, 0.2} x noise {0.002, 0.004}, 16 seeds (480 runs,
  ~1 CPU-h).  Success (record): polished < 14.9935589140 - 1e-8 (frontier_now), then candidate.py.  Secondary: <
  14.99369418 (#575's first value).  Prediction: record 10 %, secondary 50 %.
* **T2 result: fail** (288 runs, 0 failures): no 308 below 18.  8 / 288 end on the grid (18 exactly); best non-grid
  18.00000023 (near-grid, not looked at further), then 18.0000010, 18.0014.  Better than H1's chain (18.0005 in ~8
  CPU-h) at ~0.5 CPU-h, but the crossing did not happen.  Prediction 15 %: miss.
* **T3 result: fail** (control 529 runs; soft 1,058): 0 removals below either plateau.  At 232 the median outcome IS
  the plateau (+3.6e-12, every arm): Göbel's strips relax straight back after a removal; 295 soft arms median +9e-4
  (they move, but up).  Both real breaks were new layouts (itsnaka composition, ry-xu), not local re-hardenings.
  Prediction 30 % any: miss.  Plateaus need a layout generator, not this finisher.
* **T1b result**: best 14.9936138 (+5.5e-5 vs frontier 14.9935589; below #575's first value 14.99369418): record miss
  (prediction 10 %), secondary hit (50 %).  Rounds keep improving: harden.py now has `--rounds R --keep K`
  (monotone; each source restarts from its K best distinct outputs).
* **T1c, pre-registered**: continue from T1b's 4 best distinct, 4 rounds, keep 4, r0 {0.05, 0.1, 0.2} x noise 0.002,
  8 seeds (~400 runs, ~0.6 CPU-h).  Success: < 14.9935589140 - 1e-8 (then candidate.py).  Prediction 15 %.
* **D241 (dev, not a test)**: 241 from 347 with a looser held box: hold_f {1.002, 1.005} x r0 {0.15, 0.25} x noise
  {0.004, 0.008}, hold 30000, 12 seeds (96 runs).
* **Explorer kind `harden`** (explore.py `mv_harden`; off unless `--weights harden=W`): re-harden the whole packing,
  r0 {0.05, 0.1, 0.2} x noise {0.002, 0.004} drawn per proposal, then the usual screen / polish.  Smoke test (110 record
  start, 2 procs x 1.5 min, weight 50): proposals run, harden children among the new sub-11 basins.
* **HR race, pre-registered**: `bench110.py race --arms control:rand harden:rand:--weights=harden=4 --budget 8
  --check 2 --slots 10 --out runs/bench110/hr` (10 slots: shares the machine with T1c / fresh) (weight 4 ~ 1 in 6 proposals).  Question: does harden raise basin entries
  per CPU-h at 110 (the explorer's home ground; harden proposals cost ~2-3x a kick)?  Rule: adopt into the default mix
  if the ratio's lower bound > 1; drop if the upper bound < 1; otherwise keep as an opt-in kind (n = 110 is not a band
  family, so a null here does not speak to band n).  Prediction: open (ratio 0.8-1.5).
* **D241 (dev)**: looser held box is worse: 0 / 96 below the seed, best +7.9e-3; many runs end at +1.19e-2 = the 16 grid.
  At 241 the dynamics fall to the grid rather than merge the bands.  241 stays open (Ostashev's "soft" may not be
  corner rounding; his batch size unknown).
* **MB multi-band hunt, pre-registered**: the store-best packing of each n in 100-324 whose tilted squares form >= 2
  bands (`pk.bands --scan`; 18 n: 126 132 153 171 175 182 198 204 206 232 259 264 272 273 290 291 292 303; ids in
  `runs/harden/mb_ids.json`), harden.py 3 rounds, keep 3, r0 {0.1, 0.2} x noise {0.004, 0.008}, 6 seeds (~170 runs per
  n, ~5 CPU-h).  Target frontier_now(n).  Success: a certified packing below the live frontier (candidate.py).
  Prediction: 1 of 18 (Ostashev: no effect on two identical pieces 175 / 204 / 303; the merge worked 1 / 2 in T1).
* **T1c result: miss** (prediction 15 %): stopped in round 3 (load 24; rounds 1-2 added nothing).  The chain sits at
  14.9936138080 (+5.5e-5 vs the frontier) from round 0 on: a fixpoint of the operator.  The record (#575 refined,
  14.9935589) is a different minimum in the same V family; re-hardening does not hop between neighbouring V minima.
* **Fresh pass f1014b** (#586 / #593 packings, 45 new; first try f1014 stopped for T1): 0 hits at n <= 324; 3
  certified hits above (candidate.py: both checkers, cert True): **365 19.6213203436** (store best 19.7071068;
  composition frontier 199 + 26), **370 19.7677307463** (19.8086632; Couzo 175 crop + rectsearch 38-piece), **460
  21.9599696830** (21.9838629; SQUISH #481 crop + Couzo #586 crop).  Composer target was frontier_now (L-step closed,
  envelope).  Not posted: check the beyond-horizon lists (Couzo, Mishapolk, Fang, register T-152) first; ask Evan.
* **FH (background)**: harden 1 round on the new public packings n <= 324 (r0 {0.1, 0.2}, noise 0.004, 4 seeds); success = below frontier_now.
* **HR result: open, leaning worse** (8.0 CPU-h): control 47 entries / 3.99 CPU-h (11.8 per CPU-h, matching G1's 11.1),
  harden 31 / 4.02 (7.7); ratio **0.66 [0.42, 1.03]**.  1 of 31 entries came from a harden proposal.  Rule: keep
  `harden` opt-in (off by default).  At 110 (not a band family) harden proposals cost CPU without adding entries.
* **MB result: 0 / 18** (3 rounds, ~5 CPU-h; prediction 1): 13 n come straight back to their seed (+3e-12); closest
  non-seed outputs 206 +1.4e-11 (a re-found copy), 303 +7.2e-5, 132 +2.5e-4, 273 +7.5e-4.  Re-hardening the
  frontier's multi-band packings gives nothing at these settings.
* **S1 stiffness variant, pre-registered (dev-sized)**: `harden --stiff0 k0` = rigid or rounded shapes with the pair
  penalty scaled by k0 during pre, ramped k0^(1-u) -> 1 while the box is held (soft = finite-stiffness overlap, my
  second reading of Ostashev's "every square soft").  241 (347) and 211 (94): r0 {0, 0.1} x stiff0 {0.01, 0.1}, noise
  0.004, hold 30000, 16 seeds (128 runs).  Question: does 241 merge (<= 15.9876)?  Prediction 10 %.
* **S1 result: miss** (prediction 10 %): 241 best 15.9879770 (one run, -1.6e-4 vs seed, r0 0.1 stiff0 0.1), not <= 15.9876;
  most runs again at the 16 grid.  211 worse than rounding alone (best +3e-4 vs seed).  Stiffness ramps are not the
  missing piece for 241.  Harden as built: a band-merge finisher that works at 211, nothing else so far.
* **Family map (Evan 10-14: "squint at the triangle: which packings are out of place?")**: `familymap.py` (k = ceil
  sqrt n, c = k^2 - n; codes G / D 45-deg block / S 45-deg strip / V one low-angle band / 2 bands / J jumble;
  `runs/familymap.json`).  The 45-degree cells sit at a fixed offset m = n - (k-1)^2 <~ 9 on every row (89, 104, 109,
  125, 149-150, 172-174, 200-203, 229-233, 260-265, 294-298) while the V / band region grows with k.  Neighbour bump
  s(n) - (s(n-1) + s(n+1))/2: positive = weak vs both neighbours.  45-deg cells with positive bumps: 89 +0.009, 125,
  150 +0.013, 201 +0.015, 203, 229 +0.013, 260 +0.022, 295 +0.022 (147 +0.037 is the plateau).  109 is a dip (-0.001).
* **FH result**: harden 1 round over 109 new public packings n <= 324: 0 below frontier_now.
* **Policy (Evan 10-14)**: no compute above n = 324 without an extremely specific test reason.  Posted jlevy#596 (s(365),
  s(370), s(460), fresh-pass compositions).  Time-lapse video: background agent, <= 3 procs.
* **M0 Mondrian relax positive controls, pre-registered** (Evan: where a Mondrian is the record it is a not-yet-relaxed
  variant; the trick is relaxing after construction).  relax1 (10-12): 74 / 95 compositions did not move under fq soft
  ops ("jammed against grid-like upright blocks"); Ostashev's 304 relaxed "as its aligned blocks shifted".  Inputs:
  (a) Ostashev's own 176 composition (store 291603, #562, 13.808663; his re-hardening -> 13.7958723), (b) our 5 composer
  304 layouts (`runs/relax/comp304.jsonl`, 17.9474-17.9593; fq soft best 17.94156).  Grid r0 {0.05, 0.1, 0.2, 0.3} x
  noise {0.002, 0.004, 0.008}; 6 seeds per input (launched so; I first wrote 16 for 176); 432 runs.  Success: 176 < 13.807637
  (pre-#575 frontier, Couzo) [strong: <= 13.7959]; 304 < 17.929566 (pre-#575 frontier).  Secondary: 304 drop > 1e-2
  (beats fq soft's 9.9e-3).  Prediction: 176 success 40 % (it is T1's kind of move: two half-bands through a junction),
  304 success 10 %, secondary 50 %.
* **M0 result**: 176 **success (weak)**: best 13.8027534 from Ostashev's seed (< Couzo 13.807637; >> his 13.7958723);
  many runs end exactly on Couzo's 13.80764 (same seed).  304: miss, best 17.94289 (drop 8.5e-3; fq soft got 9.9e-3);
  secondary miss.  Prediction (176 40 %, 304 10 %, secondary 50 %): 176 hit weak, 304 miss.
* **M0b, pre-registered**: 176 rounds from M0's 4 best distinct, 4 rounds, keep 4, r0 {0.1, 0.2, 0.3} x noise {0.002,
  0.004}, 8 seeds.  Success (strong control): <= 13.7959.  Prediction 25 %.
* **`--hold-move 1`** (harden reharden): the box is free under pressure mu during pre / hold, i.e. the squares soften and
  the box shrinks while they are soft, then harden under pressure (Nakajima's morph phase started from a packing).
  Hypothesis: the held box is why aligned blocks do not shift (relax1: 74 / 95 jammed; M0 304 drop 8.5e-3).
* **M1, pre-registered**: same inputs as M0 (176 seed + 5 layouts of 304), `--hold-move 1`, r0 {0.1, 0.2, 0.3} x noise
  {0.002, 0.004} x mu {0.02, 0.005}, 6 seeds (432 runs).  Success as M0 (304 < 17.929566; 176 <= 13.7959 strong).
  Secondary: 304 drop > 1e-2.  Prediction: 304 success 10 %, secondary 40 %; 176 strong 20 %.
* **M0b result: strong success** (prediction 25 %): round 0 13.7959023, round 1 **13.795872318439 = Ostashev's #575
  s(176) exactly** (frontier tie, no gain), then fixpoint.  So harden + rounds rebuilds #575's 176 and 211 moves from
  their seeds; 241 still not.  Mondrian-type composition -> relaxed V works when the pieces carry half-bands (176);
  pure upright-block relaxation (304) does not yet.
* **M1 result: miss** (prediction 304 10 %, secondary 40 %, 176 strong 20 %): soft compression with the box free is
  worse everywhere (304 best 17.94434, often the 18 grid; 176 best = Couzo's 13.80764).  The held box is not why
  blocks stay put.  Reading: harden's 304 drop (8.5e-3) is close to Ostashev's (1.0e-2 from 17.9395); our 304
  layouts start 8e-3 worse than his.  So test harden on compositions already near the frontier.
* **M2 Mondrian hunt, pre-registered**: relax1's 95 composer near-arms (45 n <= 324, need < 0.015 vs the 10-12
  frontier; fq soft moved only 21 / 95), harden 2 rounds, keep 2, r0 {0.05, 0.1, 0.2} x noise {0.002, 0.004}, 3
  seeds (~1,700 + 3,400 runs, ~7 CPU-h).  Target frontier_now(n).  Success: certified below the live frontier.
  Secondary: inputs that drop > 5e-3 (relax1: 3 / 95 with fq).  Prediction: success 1 of 45 n at most (5-10 %);
  secondary 20+ inputs.
* **M2 result: 0 successes** (prediction <= 1 n, 5-10 %; 5,130 runs, 0 failures, ~8 CPU-h): harden moves 85 / 95
  composer layouts (fq soft: 21 / 95), but only 5 drop > 5e-3 (prediction 20+: miss).  Exact frontier ties at 26, 27,
  38, 50, 85, 104, 174 (rigid compositions that are the record); closest non-tie 234 +1.1e-4, 106 +3.2e-4, 227
  +3.6e-4, 53 +3.5e-4.  Reading: the jam is broken (blocks do move) but the moves do not find the relaxed variant.

## Night 10-14/15 (Evan: "the cpu is yours for the night. Go record hunting"; ~128 CPU-h)
Split ~60 % hunt / ~40 % tests + features.  Nothing above n = 324.
* **NH night hunt, pre-registered**: chain.py over the 139 non-grid n in 30-324 (`runs/night1014_ns.json`, from
  familymap), live frontier starts (`runs/live1014`, register b8c9303 + #586 / #593 synced), `--seed-ks 1,2,-1
  --first-pass`, 6 slots x 2 procs x 4-min rounds, 7 h (~84 CPU-h), post-perf engine (3x, C40 / G1).  Success: a
  candidate.py-certified packing below the live frontier.  Prediction: 0-2, point 0.7 (10-11/12: ~1 per 100-200 CPU-h
  on the old engine; /3 for the engine; the frontier got a big refresh from #593 today, which hurts).
* **Converter screen deferred**: our piece library is (h + d) x h only (rectlib re-optimises with an integer
  difference), so every piece is on the frac(W) = frac(H) diagonal by construction; an off-diagonal what-if needs pieces
  re-optimised in free W x H boxes first (a build, not a screen).  Noted for Evan.
* **NM near-miss follow-up, pre-registered**: M2's 2 best distinct outputs at 53, 106, 227, 234, 258, 299 (+1.1e-4 ..
  +7e-4), harden 4 rounds keep 3, r0 {0.05, 0.1, 0.2} x noise {0.001, 0.002, 0.004}, 6 seeds, 4 procs.  Success:
  below frontier_now.  Prediction 5 %.
* **CP composer probe (feeds Mondrian hunt v2)**: compose2 n 20-324 vs today's frontier, probe 0.015, 60 s per n, 4
  procs; near arms (need < 0.015) -> harden 2 rounds (M2 recipe).  Success: below frontier_now.  Prediction 5-10 %.
* **CP result**: 425 near arms at 59 n <= 324 (need < 0.015), many exact ties (need 0: 38 67 82 101 104 122 125 145 147
  170 171 173 174 197 198 200 201 226 229 231 ...).  0 compositions below the frontier.
* **MH2 Mondrian hunt v2, pre-registered**: up to 3 distinct arms per n (inputs `runs/night1014/mh2in`), harden 2 rounds
  keep 2, r0 {0.05, 0.1, 0.2} x noise {0.002, 0.004}, 2 seeds.  Success: below frontier_now.  Prediction 5 % (M2 on
  stale arms: 0 / 45 n).
* **NM result: 0** (prediction 5 %): after 4 rounds 53 +3.5e-4, 106 +3.2e-4, 227 +3.6e-4, 258 +5.2e-4, 299 +1.6e-4
  (no movement after round 1); **234 still descending** (+1.06e-4 -> +4.1e-5 in round 3).  234 gets more rounds (bandit
  arm with momentum).  MH2 launched.
* **MH2 result: 0** (prediction 5 %): 128 arms / 59 n, 2 rounds; ties stay ties (no drop > 1e-8 from any record-tie
  composition), non-ties best 182 +1.8e-3, 129 +2.0e-3.  Harden does not relax our compositions below the frontier.
* **N234**: 234 continuation from NM's 4 best distinct, 6 rounds keep 4, same grid, 4 procs; success below frontier_now.
* **N234 result: 0** (62 min, 4 procs, ~4 CPU-h, 336 runs/round): stalls at +4.12e-5 from round 3 on (fixpoint, like
  211's T1c).  Harden's descent at 234 was a slow tail into a local fixpoint, not a path to the record.
* **NB bump-target chain, pre-registered** (freed N234 cores, ~2 h): chain.py on the familymap neighbour-bump / 45-degree
  cells 89 109 147 150 201 229 260 295, 12-min rounds (3x NH's), 2 slots x 2 procs.  Hypothesis: longer rounds reach
  what NH's 4-min rounds miss on the out-of-place cells.  Success: certified below live frontier.  Prediction 2 %.
* **NH result: 0** (prediction 0-2, point 0.7; 593 rounds over 139 n, 424 min x 12 procs, ~85 CPU-h): no candidate below
  the live frontier; every lineage best stays at or above the register.  First miss of the post-perf engine on a full
  non-grid sweep: the #593 refresh plus our own earlier passes have mined the chain arm at this budget.
* **NB result: 0** (prediction 2 %; 20 rounds, 233 min x 4 procs, ~15 CPU-h): lineage bests stay above the record at all
  8 cells (89 +2.5e-3, 109 +2.1e-3, 147 +2.9e-2, 201 +0.18); chain restarts from the 45-degree records never improve them.
  Reading: the out-of-place cells need constructors (Evan's thesis), not longer local search.
* **Night total**: ~128 CPU-h, 0 records at n <= 324.
* **Strip excision (dev, ~0 CPU; `excise.py`)**: carry a family one row down, n' = n + 2k - 1 (k = ceil sqrt n; same
  column m), cut a vertical + a horizontal unit strip from best(n'), side s(n') - 1.  Best cuts keep too few squares:
  240 -> 211 keeps 200 (allowance 29, short 11), 130 -> 109 keeps 103 (short 6); bands crossing the strips lose extra
  squares.  Short 6 at side 10.9 ~ 0.27 in side vs 0.045 slack (s(130) - 1 = 10.9045 vs diamond 10.9497): dead as a
  constructor.  (The row-above family beating the L-step by 0.045 at 109 is still the "pattern boundary" signal.)

### s(109) region census (10-15, analysis, ~0 CPU; uncertified)
Store + 87 chain-round archives at 109.  Sub-11 distinct sides ~31,700 (1e-8; 3,383 before 10-09), but side bins
saturate (4,600 / 5,000 at 1e-5).  Role groups (L/B/axis/J, L-B swap merged) by best side: < 11 716 (123 before 10-09),
< 10.98 445 (88), < 10.96 112 (23), < 10.955 27 (10).  New groups at every depth since 10-09 (H3's 48 rounds the bulk):
broad, not saturated; best non-diamond 10.9518715 unchanged since 10-12.  No exactsolve census at 109.

## Architecture review + TT2 (10-14 day; Evan: "TT2 should be our best guess at the strongest collection of settings")
Review, findings, triage: `../../tasks/packer-arch-2026-10-14/`.  Built today (classes: functional / additive):
* chain (functional): `--reward value` + store-wide novelty ref (`--novel store`), store library starts via `pk/parents.py`
  (`--lib store --lib-policy banded`), pooled move-kind bandit per n band with no decay and full carry (`--bandit band`;
  A1 `analysis/bandit_timescale.md`: no drift at round scale, per-n learning not viable, current setting keeps 60-95 %
  of the kind information), `--as-of` time travel, no silent 10-05 fallback, imports each round, ledger rows.
  Pre-10-14 behaviour = `--reward legacy --novel run --lib run --bandit n`.
* store (additive): `expansion` ledger, `packing.cluster` column, explorer `cfp` (fq contact fp; A2: not an identity key).
* `pk/cluster.py`: leader balls in matching distance, r = 0.42 / sqrt n, deepest first, image prefilter (A2
  `analysis/similarity110.md`: at 110 40 clusters, A / B apart, 72 % of explorer steps inside a cluster vs 6 % null);
  same-packing key = side at 2e-9.  At other n clusters are mostly 1-3 minima, so `pk/parents` draws a gap decade
  first (`banded`), then a cluster, then a minimum.
* `pk/watch.py` (A3): jlevy issues + 10 watched repos -> pending entries + every parseable file; first run +163 packings.
* Store import of NH / NB / H5 (+436k rows; nothing since 10-12 18:46 had been imported).

### TT2 pre-registration (time travel, cutoff 2026-10-10T06:00Z = 10-10 00:00 local)
* Why this cutoff: store `added` = import time, store built 10-09 20:15 local, so earlier as-of pools are empty;
  43 later events over 28 n at n <= 324, all inside the n set.
* n set: NH's 139 non-grid n in 37..308 (`runs/night1014_ns.json`, familymap; fixed before looking at events).
  Starts: `eventbench starts --cutoff 2026-10-10T06:00` (`runs/tt2/starts`).  Seeds 1, 2, -1; `--first-pass`.
* Arms, concurrent, compute-matched: 4 slots x 2 procs x 4-min rounds each, ~4 h wall (~2 x 32 CPU-h).
  - S (strong, today's defaults): `--reward value --novel store --lib store --lib-policy banded --bandit band --as-of T`.
  - C (control, pre-10-14 chain): `--reward legacy --novel run --lib run --bandit n` (same starts / n / seeds / budget).
  As-of filters keep both arms' own imported rounds out of S's library, novelty ref and ledger weights.
* Metrics (`eventbench score --cutoff T`): primary = n where the arm's best beats frontier_at(n, T) by >= 1e-8 (a
  record at T), per measured CPU-h; secondary = sum of recovered fractions of later public gains over the 28 event n;
  tertiary = certified (candidate.py) successes.
* Rule: keep S's settings as chain defaults unless S < C on both primary and secondary; a primary difference counts only
  with a Poisson rate-ratio interval excluding 1, else report as open.  Either way the per-component credit is unknown
  (bundle test).
* Prediction: primary 0-3 per arm (the 10-12 time-travel pair, 2 x 33 CPU-h, found 1: s(154)); P(S > C on primary) 0.55;
  secondary S ~1.3x C.  Not tested here: the composer / fresh ingestion (not chain settings).
* **TT2 census analysis (secondary; registered 10-14 ~12:30 local, before any scoring; only the 2-round health check
  had been seen: no gain at n = 39, 40 in either arm).**  Reference per n: store rows added before T; each arm scored
  alone.  New = a polished, grid-free sub-k row from that arm's rounds > 2e-9 (side key, A2) from every reference side.
  Per measured CPU-h per arm: (1) distinct new sub-k minima; (2) new regions = new minima that start their own cluster
  (`pk/cluster` leaders vs the as-of reference clusters, r = 0.42 / sqrt n); (3) near-frontier new minima, within 1e-3
  of frontier_at(n, T); (4) S / C overlap on (1), (2).  Certification sample: ~40 random new minima per arm (sub-k, within
  1e-2 of the frontier) through exactsolve for the false-new rate (historically 20-25 % not-min), ~1-2 CPU-h after TT2.
  Prediction: S / C on (2) >= 1.5x, (1) ~1.2x, (3) ~1.3x.  Caveat: S's reward and parent draw partly optimise these.
* **How the TT2 metrics combine (Evan 10-14):** the record metrics are the strongest evidence (actual later records
  recovered, recovered fraction of later gains, certified new records).  If the arms look similar there and one is
  clearly better on new-basin diversity ((1)-(3)), read that as a hint in its favour, not a decision.  No fixed blend.

### TT2 results, record metrics (scored 10-14 ~15:00; census analysis still to run)
* Both arms ran to completion (~57.5 measured CPU-h each).  Scoring fix first: chain imports every round into the store,
  so `frontier_now` scored an arm against its own rounds (n = 146) and f64 re-import noise (~1e-14 .. 2e-10) counted as
  "later events" (118 n instead of 28).  `eventbench score` now leaves out the scored dirs (`--runs-out`, give both arms)
  and needs a later gain > 1e-9 (`EV_THR`).  Corrected: 28 event n, as pre-registered.
* Primary (certified success, beat frontier_at(T) by >= 1e-8): **S 3 (146, 305, 308), C 5 (132, 209, 270, 305, 308)**
  (C's 146 failed exactsolve as raw output; after `exact/descend.py` it lands exactly on #505, a rediscovery).
  Secondary (sum of recovered fractions over the 28 event n): **S 1.91, C 2.84**; C rediscovered the later 209 record (Couzo, #488 / #506) exactly.
  Rule: the primary difference is far inside Poisson noise (open), so S's settings stay; point estimates favour C, and
  my prediction (P(S > C) 0.55, secondary S ~1.3x C) did not hold.
* **New record from S: s(146) <= 12.583782190152, posted as jlevy/squares#606** (Evan: "Let's submit it."), 1.38e-8 below
  Barnett's #505.  From Xu's #432: two bigkicks at round 180 -> 12.5837821905 (not jammed: corner-corner MILP dS -1.3e-5)
  -> `exact/descend.py` -> certified strict local minimum mod 59 flat motions.  #505 also certified (mod 61); regularized,
  the two agree to 1e-7 except a 10-square cluster near 42 deg (moves <= 0.056; no shared angle); contacts 556 vs 577,
  corner-corner 69 vs 63.  Comparison `candidates/tt2_n146_vs505.json`, page
  https://claude.ai/artifact/T675THXVLJSM48o7o9ptP2.  No closed form: 4 coupled angle classes, factors of degree 15-32.
* Lesson: 2 of 9 raw successes were not jammed and needed the corner-slip descent before certification (one of them is
  the record).  Candidate: run `descend.py` inside `candidate.py` / the finisher whenever exactsolve says not jammed.

### TT2 census analysis (pre-registered secondary; `analysis/tt2_census.py`, `runs/tt2/census.json`, `census_cert.json`)
Measured CPU-h S 57.6, C 57.5; reference = as-of store rows (sub-k, grid-free); 139 n.
| metric | S | C | S / C per CPU-h | prediction |
|---|---:|---:|---:|---|
| (1) distinct new sub-k minima | 67231 | 57036 | 1.18 | ~1.2 |
| (2) new regions (own leader ball) | 46502 | 39337 | 1.18 | >= 1.5 |
| (3) new minima within 1e-3 of frontier_at(T) | 901 | 604 | 1.49 | ~1.3 |
* (1) and (2) are driven by small n, where the as-of reference is nearly empty (one minimum at 37 and 39; 20 n with
  < 20 reference keys): S / C by band 37-71 1.49 / 1.46, 72-150 1.07 / 1.08, 151-240 0.88 / 0.85, 241-308 0.85 / 0.81.
  Per n, S > C at 60 n, S < C at 79 (sign test p = 0.13; regions p = 0.31).  No consistent diversity gain from S.
* (3) is consistent: S > C at 47 n, S < C at 20, 72 ties (p = 0.001); every band S / C 1.35-2.14; 304 alone gives 223 vs
  66, without it 678 vs 538 (1.26).
* (4) overlap: 2891 new side keys in both arms; 850 / 849 new regions inside the other arm's new regions (~2 % of each).
* Certification sample (40 per arm, within 1e-2): S 25 certified, 3 weak, 11 not-min, 1 unresolved; C 27 / 1 / 11 / 1.
  Same false-new rate (28 %, mostly 'not jammed'), so the count ratios stand as they are.
* **Reading (Evan's rule: record metrics first, diversity a hint):** records slightly favour C (inside noise); diversity
  favours S only near the frontier, and that part is consistent across n.  Keep S's settings (pre-registered rule); hint
  in S's favour for near-frontier exploration (plausibly the banded parents + value reward: bundle, no per-component
  credit).  Predictions: (3) held, (1) held only through small n, (2) failed.
* Ops note: the first certification launch ran in-process exactsolve without thread caps (~47 threads per worker, load
  ~300); `pk/jobs.py` now caps math-library threads at 1 by default (`5e9f662`).

* **Correction (10-15): TT2's CPU-h were double-counted.**  `analysis/tt2_census.py` summed the chain manifest (which
  includes its children, i.e. every round) and the round manifests: true measured CPU is **S 30.0, C 30.1 CPU-h** (rounds'
  own manifests 27.6 / 27.4; proposal-level CPU 27.3 / 27.2), not 57.6 / 57.5.  Every per-CPU-h rate above is 1.9x too
  low; ratios between arms are unchanged.  `cpu_h` fixed (chain manifest only).

### TT2 at n = 110 (one 4-min round per arm; `analysis/tt2_cert110.py`, `runs/tt2/cert110{,_desc}.json`)
Every new sub-11 minimum of both arms (vs the as-of store: S 6, C 12) through exactsolve; the 12 'not jammed' ones
through `exact/descend.py` and exactsolve again.  Result: **no new exact minimum, no new region, no record.**
* Directly certified: 6, all among the 939 known exact minima (10.996923243, 10.997025578, 10.998243315, 10.998443828,
  10.999196864, 10.999651054).
* After the descent, 11 / 12 land on known minima: S's queued near-miss 801876 (+4.45e-6) is the **record minimum
  itself** (10.99678339663159164; the finisher's polish got there too); C's five 10.99688699 rows are two known minima
  (…868489, …876103); the rest are known minima 10.99801257 (weak), 10.99834170, 10.99841239, 10.99999749.  C 801753
  (10.99692215) is still not jammed after 2 rounds (dS* -8e-7).
* All sit in known leader regions (matching distance to the nearest known minimum 0.001-0.004, one 0.021); 9 of 18 in
  the depth-rank-1 region (leader 10.99688705).  Reading: at 110 one round per arm only re-finds the mapped landscape;
  "new" against the as-of store here means "not yet imported", not new.  18 / 18 'new' rows were known minima or
  unconverged copies of them (the census's 28 % false-new rate understates it at heavily mined n).

## TT3: feature tests on a cheap replay design (10-15; Evan: one feature, one arm, keep power at lower compute)
Review 10-14/15 (Evan: are we doing too much exact descent?): in TT2 81 % of proposals got a full polish, only 1-4 % of
those landed within 1e-4 of the round's start best (2-11 % within 1e-3), and polish is ~60-67 % of fq time.  First feature:
**polish gating**; the movement-heavy arms (GLS feasibility search, bigger melts) come next, on top of it if it holds.

### Design (power analysis: `analysis/tt3_score.py`, numbers from TT2's per-round data)
* TT2 gave most n one 4-min round (n <= 240: 94 of 102 n one round).  Metric (3) of TT2 (near-frontier new minima,
  within 1e-3 of frontier_at(T)) has overdispersion ~2.4 vs Poisson; the bootstrap over n gives x/÷ ~1.35 at n <= 240.
* n set: the **39 n <= 240 where TT2 arm C found >= 1 near-frontier new minimum** (selected on C, so S and any replay of
  S are not selected on their own noise): 82 % of S's events (297 / 363) at 38 % of its CPU (6.5 / 17.2 CPU-h).
  Caveat: a feature that only helps at n with no near-frontier activity is invisible here; the final config test covers it.
* Rounds: `chain.py --schedule runs/tt3/schedule.json` replays S's 49 rounds on those n in S's order with S's round
  numbers (same explore seeds, removal / insertion seeds, store library picks as of T, starts `runs/tt2/starts`).  None of
  those 49 rounds used a neighbour seed source that S had already improved (checked), so the subset is self-contained.
  Not replayed: the pooled move-kind bandit state (S's came from all 205 rounds) and wall-clock timing.  Hence a
  **concurrent replay control A** (S's settings, no feature); later features reuse A if A ~ S.
* Arms (concurrent, 4 slots x 2 procs each, 4-min rounds, ~1 h wall, ~6.5 CPU-h each), all with TT2 S's chain settings
  (`--reward value --novel store --lib store --lib-policy banded --bandit band --as-of 2026-10-10T06:00 --seed-ks 1,2,-1`):
  - **A**: replay, no change.
  - **G**: replay + `explore.py --polish-gate 2e-3`: full polish only if the screened side < min(frontier_at(T), run best)
    + 2e-3 (above the 1e-3 metric window + the screen's error, so every candidate the metric can count is polished).  The
    rest enter the archive at the screened side as stepping stones (`unpolished`; store status `screened`): never elite
    parents, never the run best, identity and novelty at UNP_TOL = 1e-5 (screen error), bandit reward as usual (value
    reward, gap > 2e-3, so <= 0.5 x robustness).
* Primary: near-frontier new minima per measured CPU-h (round manifests), **G / A**, rate ratio with a 95 % bootstrap over
  the 39 n.  Rule: G passes (base for the next feature) if the interval excludes 1 from above; fails if it excludes 1 from
  below; otherwise open (then: decide on throughput + secondary, or one more replicate pair).
* Secondary (not decisive): records vs frontier_at(T) (eventbench score, G vs A); proposals per CPU-s and polished share;
  best gap per n (sign test G vs A); **A / S-historical** = replay fidelity (expect ~1; if far from 1, later tests use A).
* Predictions (written before launch): proposals per CPU-s G / A 1.6x (n < 72) .. 2.5x (n >= 150); primary G / A ~1.4,
  P(interval excludes 1 from above) ~0.5, P(excludes from below) ~0.1; A / S within [0.75, 1.33].
* Smoke 10-15 (`runs/tt3/smoke`, rounds 2 and 3 at 0.5 min, gate on): replay numbering works, 4 % of proposals polished,
  22.7 proposals per CPU-s vs S's 14.4 at n = 38-39 (screens dominate at small n).

### Result (10-14 evening, scored `analysis/tt3_score.py runs/tt3/G runs/tt3/A --json runs/tt3/score.json`)
Both arms ran on the tree later committed as `d8843de` (the run manifests say "dirty" at 1038892: launched minutes
before the commit, no edits in between).  Measured CPU (chain manifests) A 6.86, G 6.87 CPU-h; round manifests 6.54 / 6.53.
* **Primary: G / A = 2.08, 95 % bootstrap over n [1.75, 2.56]: PASS** (predicted ~1.4, P(pass) ~0.5).  Near-frontier
  new minima: A 282 (43.1 / CPU-h), G 590 (90.3 / CPU-h), S-historical 297 (45.4 / CPU-h).
* Replay fidelity: **A / S = 0.95 [0.79, 1.10]**, inside [0.75, 1.33]: the replay design holds, so later features reuse it
  with G (already run) as the baseline instead of a fresh control.
* Throughput (proposals per CPU-s, G / A): n 37-71 1.81x, 72-150 2.29x, 151-240 2.21x (predicted 1.6x .. 2.5x).  Polished
  share 81 % -> 6 %.  The gain is mostly more proposals, not better ones: near-frontier events per proposal are ~equal
  (A 282 / 79.0k, G 590 / 154.2k).
* Robust to the duplicate-minimum question (f64 polish noise at large n): recounting distinct sides at tolerance 1e-7 /
  1e-6 / 1e-5 instead of 2e-9 gives G / A 2.03 / 1.93 / 1.95.
* Spread over n: G higher at 29 n, lower at 6, tied at 4.  Best gap per n (sign test, no event = worse): G better 17,
  worse 7, tied 15, p = 0.064.
* Records vs frontier_at(T) (eventbench, underpowered as expected): A 2 (146, 209), G 3 (132, 146, 209); G's 132 is the
  current frontier (posted after T) re-found to 4e-12.  No new records.
* Not measured, so not claimed: far (non-near-frontier) minima and new regions.  G polishes far fewer far states, so the
  store gets far fewer polished far minima and many more `screened` rows (sightings A 28.7k, G 48.5k).  Evan's open
  question (are looser far packings useful add / remove donors?) is not tested by this design; the final config test
  on all n should include TT2's region metric.
* Decision: gating at 2e-3 becomes the base for the next feature (movement-heavy arms), scored as feature+G vs G.

### TT3 round 2: melt and GLS arms on top of G (pre-registration, 10-14 ~21:00 local, before launch)
Evan 10-14: try the other features on top of gating; no new baseline.  **Baseline = G** (TT3 run above): same replay
schedule, same chain settings, the replay design passed its fidelity check (A / S 0.95).  Load matched: two arms
concurrent, 4 slots x 2 procs each, as G / A were.  Bandit state starts empty in every arm (as in A and G).
* **M** = G + `--melt 1 --bigmelt 1`: the regional melt (R 1.5-3.5, rmax 0.05-0.2) plus a new `bigmelt` kind (R 3.5-7,
  rmax 0.1-0.25; 0.1-0.5 s per melt at n = 211, ~2-4 % box growth, the screen re-compresses).  Melt was off in TT2 / G;
  the fixed melt was neutral for the census rate at 110 (C10, 0.91 [0.75, 1.11]).
* **L** = G + `--gls 1`: new kind `gls` (`anneal gls`, src/anneal.rs; `explore.gls_side`), Sparrow-style guided local
  search at a fixed target side min(parent, gate line) x (1 - d), d log-uniform [1e-4, 3e-3]: per axis, squares beyond
  a random cut slide in by the shrink (overlap on two seams); then alternate gradient relaxation of the weighted squared
  SAT penetration depth (20 steps), relocation of each colliding square to the best of 32 global + 16 local samples
  (+ coordinate refinement), and GLS pair-weight updates (x1.2-2 on overlapping pairs after 8 passes without 1 % gain,
  decay 0.95 to 1).  Budget 2 ms x n wall per call.  Output (the least-overlap state) goes to the usual screen.
  **v0 caveat:** on bench cases it resolves slack spread through the box (a record scaled up 5 % then shrunk 4.5 %,
  n = 55) but not a single hole far from the seams (record minus one square, shrink 1e-3, n = 55) or loose n = 211:
  collective motion is weak.  A fail here speaks to this v0, not to fixed-side feasibility search in general.
* Both kinds are chosen by the move-kind bandit (reward per CPU-s, floor 3 % per kind), so a useless kind costs only
  about its floor share of CPU; a pass needs the kind to earn a large share.  Per-kind CPU share and new-entry rates
  are reported as diagnostics.
* Primary per arm: near-frontier new minima per measured CPU-h, **X / G**, 95 % bootstrap over the 39 n
  (`analysis/tt3_score.py runs/tt3/M runs/tt3/L runs/tt3/G --base runs/tt2/S`).  Same rule as G: pass if the interval
  excludes 1 from above, fail if from below, else open.  Two arms, no multiplicity correction (each is its own decision;
  a pass at the margin gets a replicate before adoption).
* Secondary: records (eventbench), best gap per n (sign test vs G), proposals per CPU-s, per-kind CPU share.
* Predictions: M / G ~1.0 (P pass 0.15, P fail 0.1); L / G ~0.95 (P pass 0.1, P fail 0.15); bandit CPU share of the new
  kinds: melt + bigmelt 8-15 %, gls 4-8 %.
* Smoke 10-14 (`runs/tt3/smoke2{M,L}`, rounds 2 (n 38) and 89 (n 211), 1 min): kinds fire, per-proposal CPU includes the
  anneal child (gls 0.1 s at n 38, 0.85 s at 211; bigmelt 0.08 / 0.7 s).

### Result, round 2 (10-14 ~22:00; `analysis/tt3_score.py runs/tt3/M runs/tt3/L runs/tt3/G --json runs/tt3/score2.json`)
Ran at `b90e5fe` (clean).  Round-manifest CPU 6.53 CPU-h per arm, as G.
* **M / G = 1.03 [0.85, 1.22]: open.  L / G = 0.99 [0.83, 1.18]: open.**  Near-frontier new minima M 610, L 586, G 590.
  Both inside the predictions (M ~1.0, L ~0.95); neither harms, neither helps at the ~1.2x this design resolves.
* Bandit CPU share of the new kinds (above the predicted 8-15 % / 4-8 %): melt 6.9 % + bigmelt 10.4 % in M; gls 10.2 %
  in L.  Near-frontier polished proposals per CPU-h by kind (repeats included, a diagnostic only): bigmelt 1,080, melt
  1,230, gls 880, vs the arm's kick / aswap / rowslide 1,450-1,990 and mirror ~300.  So the new kinds are mid-pack
  (melt) to below average (gls) at the near-frontier metric; the bandit funds them on the novelty reward.
* Bandit allocation is itself noisy between arms with identical kinds (cross2 13.8 % of CPU in G, 9.5 % in M, 5.9 % in
  L; kicksym 10.6 / 11.8 / 24.7 %), so per-kind shares from one run are rough.
* Secondary: records M 2 (146, 209), L 3 (132, 146, 209), G 3.  Best gap per n vs G: M better 13 / worse 7 / tied 19
  (p = 0.26); L 11 / 9 / 19 (p = 0.82).  Proposals per CPU-s: M 5.71, L 6.14, G 6.56 (the new kinds cost more per
  proposal).
* Not measured: reach / new regions (melt's earlier positive was reach at s130, PACKER 10-08).  GLS v0's weakness at
  collective motion (pre-registration) is consistent with its below-average yield; a v1 would need stronger continuous
  relaxation (e.g. L-BFGS on the weighted penalty), which is close to what fq's screen already does.
* Decision: config stays **G**.  Melt / bigmelt are cheap and harmless: candidates for an arm in the all-n config test
  with the region metric, where their claimed benefit (reach) can show.  GLS v0 parked.

### Parent / move value model, offline (10-14 late; `analysis/parent_model.py build && report [--label yl]`)
Evan: is parent choice a modeling problem that generalizes?  Today: kind = context-free bandit (rate per CPU-s per n
band); parent = fixed diversity rule (70 %: a seed lineage uniformly, then a uniform descriptor cell, then best-8 of
the cell by expansions; ~9 % elite; `--adapt-cells` exists but chain.py never sets it).  No model shares across parents.
Data: 433k proposals of TT3 G / M / L (19.3 CPU-h); label y = distinct new near-frontier minimum (first in the arm),
yl = y + 0.5 x descendants' yl within the round (lineage credit).  Poisson GLM, offset log CPU, 5-fold CV over n.
* Parent gap dominates: y per CPU-h 320-480 for parents within 1e-3 of frontier_T (17 % of CPU), 210 at 1e-3..3e-3, 50 at
  3e-3..1e-2, 23 / 8 / 3 at 1e-2 / 3e-2 / 1e-1+ (**74 % of CPU at >= 1e-2 for 7 % of y**).  Screened stepping stones
  as parents: 75 % of CPU, 12.6 per CPU-h (polished 409, store seeds 134).
* Lineage credit hardly rescues far parents (yl per CPU-h 11-22 at >= 3e-2 vs 610-810 near): within a round, far
  parents do not lead to near-frontier descendants.  (Not measured: far rows reused as library seeds in later rounds.)
* Held-out deviance: const 13,713; gap 9,550; + saturation (log expansions, log archive-new children) 9,147; + kind and
  n 9,078; + near x kind 9,098 (no gain: kinds rank the same near and far).  Saturation: rate x (1+exp)^-0.78
  x (1+kids)^+0.57 (productive parents stay productive).  Kind spread near the frontier: kicksym 563, kick 513, aswap
  499 .. mirror 211 per CPU-h.
* Reallocation bound (held-out ranking, logged CPU): the top 25 % of CPU by predicted rate yields 3.7x the average per
  CPU (gap alone 3.7x; top 10 % 6x).  Upper bound only: concentrating CPU raises expansions and saturation (the
  expansion term is the only guard), and the archive's future differs under another policy.  Only an arm measures it.

### TT3 round 3: parent value model arms (pre-registration, 10-14 ~23:00 local, before launch)
Evan 10-14: try the simple model (it shifts as compute is applied through the saturation terms, and it scales to many
parents).  `explore.py --pm-share q --pm-beta b`: with probability q the parent is drawn from every archive entry below
smax with weight rate^b, rate = the 'full' Poisson fit above without kind / n terms (gap bin vs frontier_s, screened,
start, (1 + expanded)^-0.78, (1 + archive-new children)^0.57; `explore.PM_*`); otherwise the old rule.  Proposals log
`pick: pm` and `pp` (the pick probability under the model), for later off-policy scoring.
* **Cheap design** (Evan: fast iteration): 12 of TT3's 49 rounds (`runs/tt3/sub12_schedule.json`, n 55 68 132 172 176
  181 182 204 208 211 234 240), chosen by the event rate of **S and A** (not G, so the baseline is not picked on its own
  luck): 25 % of the CPU (G 1.61 CPU-h) and 59 % of S + A's events; G has 332 events there.  Baseline = G on the same
  rounds (no new control).  Expected resolution ~1.3x; the bootstrap has only 12 clusters (rougher intervals).
* Arms (concurrent, 4 slots x 2 procs each, as before): **PH** = G + `--pm-share 0.5`; **PF** = G + `--pm-share 1`;
  both b = 1.  ~1.6 CPU-h each, ~15 min wall.
* Primary: X / G near-frontier new minima per CPU-h on the 12 rounds, 95 % bootstrap over n (`analysis/tt3_score.py
  runs/tt3/PF runs/tt3/PH runs/tt3/G --schedule runs/tt3/sub12_schedule.json`); rule as before (pass / fail / open).
  A pass here is confirmed on the full 49 rounds before adoption.
* Secondary: CPU share by parent gap bin, records, best gap per n, the model's predicted vs realised rate by bin
  (calibration under the new policy: saturation is the risk).
* Predictions: PF / G ~1.5 (P pass 0.6, P fail 0.15: herding on a few parents, repeats); PH / G ~1.3 (P pass 0.4).
  The logged-CPU bound (3.7x for the top 25 % of CPU) is an upper bound.
* Smoke (`runs/tt3/smoke3`, round 89 n 211, 1.5 min, PF): 96 % of CPU on parents within 1e-3 of frontier_T (G there:
  52 %), pick probabilities up to 0.9 early (one seed dominates until expansions accumulate).

### Result, round 3 (10-14 ~23:30; `analysis/tt3_score.py runs/tt3/PF runs/tt3/PH runs/tt3/G --schedule runs/tt3/sub12_schedule.json --json runs/tt3/score3.json`)
Ran at `9af0ef3` (clean).  1.60 CPU-h per arm (round manifests), 12 rounds.
* **PF / G = 1.98 [1.52, 2.63]: PASS.  PH / G = 1.67 [1.25, 2.22]: PASS.**  PF / PH = 1.18 [1.00, 1.41].  Near-frontier
  new minima PF 654, PH 554, G 332 (S-historical 171).  Predicted PF ~1.5, PH ~1.3: both above.
* Identity check (PF polishes 38 % of proposals vs G 9 %, so repeat polishes could add side keys by f64 noise): counting
  a row as new only if no reference / earlier row is within side 1e-5 **and** matching distance (pk.cluster.md) 1e-3:
  PF / G 1.90, PH / G 1.63; at md 1e-2 (merges rattler variants) 1.73 / 1.60; side 1e-4 + md 1e-2: 1.66 / 1.58.  Side
  alone at 1e-5 / 1e-4 gives 1.40 / 1.13, but that merges distinct configurations (so it understates).
* Where the CPU went (share / y per CPU-h by parent gap bin): parents within 1e-3: G 29 %, PH 59 %, PF 77 %; >= 1e-2:
  G 59 %, PH 27 %, PF 6 %.  Saturation as expected: y per CPU-h in the 1e-4..1e-3 bin G 782, PH 603, PF 511; near
  parents' mean expansions 9.5 / 18.3 / 25.7.  The shift outweighs it.
* Best gap per n vs G: PF better 6 / worse 1 / tied 5; PH 5 / 2 / 5.  Near-count higher at 11 / 10 of 12 n.
  Records on the subset: PF 1 (132), PH 1 (132), G 1 (132).
* Not measured, and the cost to watch: PF leaves far parents almost no CPU (6 % at >= 1e-2), so far-region discovery
  (TT2's region metric; Evan's donor question) will drop.  The all-n confirmation should carry the region metric.
* Next: confirm on the full 49 rounds (PF, maybe a sharper beta), then the all-n hunt-style test with regions + records.

### H110-PM: n = 110 holdout test of the parent value model (pre-registration, 10-14 ~23:55 local, before launch)
Evan 10-14: a short, high-power test that also covers far reach (the replay cannot).  `analysis/h110_pm.py`.
* Holdout `refsets/h110_rand50` (census110_1011, 390 visible / 370 withheld certified-or-weak minima, record withheld;
  the 40 visible random starts).  Census v2's control re-found ~110 withheld / CPU-h (old recipe), census v2 smoke and
  this smoke: ~1,000 / CPU-h in the first minutes (saturating), so runs are short.
* Arms (explore.py directly, current explorer config as in TT3 rounds: `--adapt --reward value --elite-share 0.3
  --polish-gate 2e-3`, novelty vs the visible minima, `--frontier-s` = best visible 10.9967871155): **G** (old parent rule)
  vs **PF** (+ `--pm-share 1`).  4 replicates per arm x 2 procs x 10 min, all 8 runs concurrent (16 procs), ~2.7 CPU-h.
  No chain: no library / seed / frontier-lineage picks (a difference from TT3's G, same in both arms).
* Primary: certified withheld minima re-found per measured CPU-h (side key -> store minima, as refset.score; polished
  archive entries only), **PF / G** ratio of mean log per-replicate rates, 95 % Welch t interval.  Pass / fail / open as
  in TT3.
* Strata (the far-reach question): withheld by side gap to the record (< 1e-3: 34; 1e-3..3e-3: 286; >= 3e-3: 50) and by
  matching distance of their representative to the record's (< 0.05: 71; >= 0.05: 299; max 0.106).  PF / G per stratum.
* Secondary: withheld found by 150 s / 300 s (saturation check); record re-found (reps); union over reps.
* Predictions: PF / G all 1.3 (P pass 0.5); gap < 1e-3 stratum 1.8; gap >= 3e-3 stratum 0.6; md >= 0.05 stratum 1.1
  (at 110 nearly every minimum is within 1e-2 of the record, so the model does not starve these as it does the far
  parents at other n; P(PF / G < 1 at md >= 0.05) 0.35).
* Smoke (1 min each, `runs/h110_pm_smoke`): G 38 withheld, PF 49 in ~125 CPU-s each; every PF pick logged `pm`.

### H110-PM, result (10-15 ~00:15; `analysis/h110_pm.py score --json runs/h110_pm/score.json`)
Ran at `3185510` (clean); 8 runs x 10 min x 2 procs, ~2.7 CPU-h.  Bookkeeping note (Evan): at 110 the grid is 3.2e-3
above the record, and the census tracks only sub-11 minima, so every withheld target is within 3.2e-3; above-grid
states are parents here, never targets.  This test covers sub-11 rediscovery and structural (matching-distance) reach only.
* **Primary: PF / G = 1.05 [0.85, 1.30]: open** (predicted 1.3).  Withheld per CPU-h G 211, PF 223.
* Strata: gap < 1e-3 (34): **0.70 [0.46, 1.07]** (predicted 1.8: wrong direction); 1e-3..3e-3 (286): 1.20 [0.98, 1.47];
  >= 3e-3 (50, the shallow edge 10.9998-11): 0 in both arms; md < 0.05: 0.90 [0.66, 1.23]; md >= 0.05: 1.11 [0.89, 1.38]
  (predicted 1.1).  Record re-found 4 / 4 reps in both.
* Saturation: by 150 s PF / G 1.19, by 300 s 1.06, at 600 s 1.05; union over reps G 102, PF 108 of 370: both arms hit
  the reachable pool from these 40 starts within minutes, which compresses any difference.
* Why PF did not concentrate: parent CPU by side: G 44 % below 10.998 (elite share 0.3 = the 20 deepest), 14 %
  10.998-11, 41 % above 11; PF 13 % / 85 % / 2 %.  The fitted model uses *absolute* gap bins (vs best visible
  10.99679): almost the whole sub-11 census sits in the 1e-3..3e-3 bin, the near bin holds few entries, and the
  saturation term spreads mass over the big bin, so PF put less CPU near the record than G's elite picks.  PF polished
  6,236 proposals (89 % of CPU) vs G 7,853 (57 %): fewer proposals, more expensive ones, no more minima.
* Reading: the absolute-gap model does not transfer to n = 110's compressed band (every sub-11 minimum within 3.2e-3);
  the TT3 subset gain (C44) is real there but not general.  Fixes to try: gap relative to the run's own distribution
  (rank / quantile among archive entries, or gap / (grid side - frontier)), and keeping an elite share beside the model.

## Baseline on the major tools (10-15; Evan: "get to baseline on the major tools of the modern trade", in order)
Diagnosis (EVENTS tables): at n <= 324, families A (composition) 32 % -> 61 % of the gain (10-05..12 -> 10-12..15),
B (soft-to-rigid) 0 -> 29 %, J (evolved constructors) 40 % -> 0, F/G (removal / insertion) 21 % -> 1 %, H (perturbation /
polish, ours) 5 % / 6 %.  TT3 (gating, parent model) tuned H.  Order now: per family, level 1 (method fidelity: given the
named parents, does our pipeline reproduce the event?), then level 2 (blind time travel from everything public at T),
then improvements.  Evaluation on n <= 324 + time travel; the n > 324 register rebuild is a detour (TODO).
* **F level 1, pre-registered (before launch)**: `removal_bench.py run --minutes 15 --procs 2 --par 7 --out
  runs/removal_bench/b1`.  The 13 F events with a named parent "s(m) - k" (SQUISH #401 / #422, ours #552): parent = the
  frontier of m at T or a packing posted in the same update (all 13 in the store); every single removal (k = 2: 400
  random pairs) quenched (chain seed quench), 8 best distinct -> explorer 15 min x 2 procs (gate 2e-3, value reward,
  elite 0.3).  Tiers: REPLICATED (<= posted + 1e-9), BEAT-ANTE (< frontier_at(n, T) - 1e-9), no.  ~8 CPU-h.
  Note: several "removals" are large restructurings (108: s(110) record 10.99678 - 2 -> 10.90994; 179: s(180) 13.92364
  - 1 -> 13.89131), so the search after removal is the method.  Prior: 10-12 time travel did not move 179 in 1.1 CPU-h.
  Prediction: BEAT-ANTE or better 3 / 13, REPLICATED 1 / 13.  If low, the missing component is SQUISH's post-removal
  search ("iterative shrink-hopping": fixed-side basin hopping, then shrink; README wording only).
