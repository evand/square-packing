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
