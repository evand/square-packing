# Case study: the landscape of 110 unit squares below side 11 (DRAFT, 2026-10-04)


**Naming (10-05):** "candidate" / `cand110` / run label `cand` throughout = **Couzo's 2026-09-27 s(110) record**
(`../exact/results/couzo110.*`; start file `candidates/s110_10.996783403.txt`), found by us independently a week later; "record"
/ `rec110` = Ellsworth's earlier 10.99679327….  Not ours; we add only its exact optimum (2.8e-12 below his printed bound).

**Status: draft.  Numbers numerical unless marked.  Not for external posting.**

**Correction (10-04, late): not new.**  jlevy's register (https://jlevy.github.io/squares/cases/110.html) already lists Francisco Couzo, 2026-09-27: s(110) ≤ 10.996783396634359 (verified there by exact rational packing and interval arithmetic).  Our candidate is an independent rediscovery of that packing; our exact KKT point 10.99678339663159163950… and certificate S' = 10.996783396631591639612746267957 are 2.8e-12 below Couzo's printed bound (analytic optimisation of the same basin).  Ellsworth's page (and our mirror) still show 10.99679327401957.  Lesson: check the register before calling anything a candidate.

## Why this packing
s(n² − n) = n was conjectured and is now known false for n ≥ 11: Cantrell (Feb 2025) packed 110 unit squares in a square of
side < 11; refined by Ellsworth (Jan 2026, modified Schadt annealer) and Stead (Jun 2026, with Claude Fable 5) to the posted
10.99679327401957, "not yet analytically optimized".  It is the smallest n where the conjecture fails, its margin below 11 is
only 3.2e-3, and it is the obvious model for the open case s(90) < 10.  We used it to understand what a sub-k packing looks
like as a point in a landscape: how isolated it is, what its neighbours are, and what a search has to do to find it.

## The structure
Two tilted channels plus an axis frame (figure: `figs/s110_record_vs_candidate.png`).  A left channel of 22 squares at
24.9–32.1°, a bottom channel of 21 at 61.7–66.1° (i.e. 23.9–28.3° the other way), meeting at a junction near (4.5, 4.5); the
record has three near-axis junction squares (#62, #63 at 1.12°, #65 at 0.16°).  The tilted part crosses every row and every
column of the 11 × 11 grid: any packing with side < k must have no full line of k axis squares, so its tilted part is a
transversal.  The channels are not rigid lattice rows: squares are staggered, sheared and mix angles, and many rest on axis
squares or walls.

## A slightly better neighbour (= Couzo 2026-09-27, rediscovered)
Perturbing the record (Gaussian σ on all coordinates) and descending (soft penalty squeeze, then an exact-contact-aware LP
squeeze) lands in two repeated minima: the record (10.9967933) and **10.9967834** (≈ 19/141 trials), better by 9.9e-6.
In the candidate the three junction squares are exactly axis-aligned and sit on the top-right lattice; four other squares
move by more than 0.01.  The contact graph changes at the junction by a role inversion: in the record a corner of #62 rests
on a side of the tilted square #76; in the candidate a corner of #76 rests on a side of #62.  Passing between them requires a
corner–corner configuration, so they are distinct strata.  An angle homotopy (rotate #62, #63, #65 to 0 while re-minimising
everything else) rises from 10.996793 to a ridge at 10.997468 (≈ 0.5° tilt, the flip) and descends to 10.996846: a barrier of
at most 6.7e-4, far below 11.  Exact (`exact/`): KKT point S = 10.99678339663159163950…, λ > 0, PSD modulo 17 exact flat slides,
certificate S' = 10.996783396631591639612746267957 (two independent exact verifiers).  The record's own exact point is
10.99679327395374924222… (6.6e-11 below its posted value).

## The sub-11 landscape (sampled)
**Correction (10-05): the counts below come from the false-jam-era descent and are inflated.**  The clean census (`runs/cen7`,
1000 trials, every sub-11 output solved exactly, PACKER.md) finds **19 distinct sub-11 minima** (Chao1 23), the top four holding
91 % of hits: candidate (219), record (119), 10.9971887601 (236), 10.9971976540 (166).  Rewrite this section from cen7.
* **Near-degenerate ladder.**  128 sub-11 minima from 141 perturbations: ≥ 10 junction classes (which junction squares are
  tilted 0.05–2°); the record's class {62, 63, 65} (75 configs, 10.9967933–10.99864) and the all-axis class (42 configs,
  10.9967834–10.99886) dominate; within a class, many distinct minima over ~2e-3.
* **Caveat (exact solver):** slp2 can false-jam at corner–corner touches (fixed separating side); of 7 census minima tested
  exactly, 5 were not local minima (after a kick: 2 → candidate, 2 → record, 1 → 10.996866449), 1 genuine distinct minimum
  (10.99719765…), 1 unresolved.  So distinct-minimum counts below are upper bounds (inflated).
* **Census push** (480 trials: from record and candidate, σ = 0.01, 0.03, 0.05, 0.1, 60 each; soft squeeze + LP finish): 330 jammed
  sub-11 minima without full lines, **51 distinct sides** (10.9967834–10.9988563).  Basin frequencies: 10.9971888 (94 hits), the
  candidate 10.9967834 (75), 10.9972187 (35), 10.9968134 (18), 10.9968259 (16), the record 10.9967933 (15).  So the candidate's
  basin is ~5× the record's under this descent, and neither is the largest.  Nothing below the candidate; just above it rare minima
  at 10.9967838, …839, …864, …912.  36 singletons → richness estimates (Chao1) ≈ 33 at σ ≤ 0.03, 70–90 at σ = 0.05–0.1, ~700 pooled
  (unstable: one doubleton): dozens to hundreds of near-degenerate jammed minima in this funnel, more reached as σ grows.
* **Funnel size.**  Escape curve P(s ≥ 11 | σ) from the record: 0/83 for σ ≤ 0.01, 2/28 at 0.03, 12/30 at 0.1; census push
  (record / candidate): 0/60, 2–3/60 at 0.03, 7–10/60 at 0.05, 24–26/60 at 0.1; σ50 ≈ 0.12–0.15 under soft squeeze + LP finish.  Under strict local descent alone the sub-11 basin is far smaller (σ = 0.01: 0/8 below 11).
* **Other funnels.**  From family structure alone (two-channel layouts, free-tilted search) our best is a different funnel at
  11.00765: shorter left channel, 5-tall bottom channel.  It never descends below 11 under perturbation; at σ = 0.1 most trials
  collapse to full lines (exactly 11).

## What this says about search
1. The trivial-side attractor (full lines) dominates every level; rejecting full lines is essential.
2. The target funnel is not a continuous needle (σ ≈ 0.1 per coordinate still lands < 11), but it is one of many funnels near
   11 separated by many-square rearrangements: the hard part is choosing the funnel (structure), not descending in it.
3. Local descent quality matters: a penalty-gradient squeeze stalls at kinks (face-to-face contacts) ~1e-3 above true minima,
   enough to mis-rank candidates and fake a "census" of distinct minima.
4. Inside a funnel, discrete junction choices create a ladder of near-degenerate optima; small "junction flip" moves explore it.

## Open
Exact certification; counting distinct force-bearing contact graphs; finding the record's funnel from outside (the pass test);
the same analysis for s(132) and other "not yet analytically optimized" records; s(90).
