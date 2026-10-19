# Symmetric-subspace search (10-18; idea 1 of the 10-18 review: `lit-2026-10-12.md` §2, ranked first there, never run)

Code: `symscan.py` (D4 invariance of store bests), `symstart.py` (free / Md-symmetric random starts), `symharden.py`
(stage-0 harness), `symfunnel.py` (funnel width), `harden reharden --sym-md 1` (src/harden.rs: Md partners matched once
from a symmetric input; positions + velocities projected onto the Md subspace after every step; diagonal squares on
x = y at 0 / 45 degrees).  Runs in `runs/symsub/`.  These were exploratory pilots, mostly not pre-registered (said here,
not hidden); the numbers below are what ran.

## Scan
Store public bests, tol 1e-6: 38 non-integer n in 30-324 with tilts are exactly Md-invariant (38 40 52 65 67 82 89 101
104 109 122 124 125 145 147-149 170 171 173-175 197 198 200-203 226 231 233 262 265 290 291 294 296 298); also 5, 10,
27 below 30.  Their diagonal squares always read 0-degree uprights at both ends + one contiguous 45-degree run (the
Göbel strip crossing the diagonal).

## Stage 0: from random starts (no)
Same soft-to-rigid run (harden reharden, r0 0.2 / 0.35, hold 30k, box free while soft) + fq polish, free vs symmetric.
* Symmetric starts with free dynamics: symmetry lost in every run, noise or not (smoke, smoke0): not a test.
* Enforced symmetric dynamics, diagonal angles iid (smoke2, small): sym worse everywhere; n 10 / 27 / 38 / 40: record
  hits free 10 / 128, sym 4 / 128; sub-grid outputs free 84, sym 35.
* Diagonal = contiguous 45 run, m from ceil(S) - 1 .. + 2 (run1, n 10 27 38 40 52 65, 24 seeds, r0 0.35): hits free 6
  (all n = 10), sym 3 (n = 10); sub-grid free 110 / 144, sym 29 / 144; sym medians +0.2 .. +0.7 worse.
  Reading: from scratch the symmetric subspace does not find funnels; random symmetric starts are worse starts.

## Funnel width (yes)
Positive control (symctl, scratch): record loosened x 1.02 / 1.08, harden r0 0.1, 3 seeds, n 27 52 65: sym returns to
the record 17 / 18, free 4 / 18 (free drifts to +0.02 .. +0.12 alternates).
**FW1** (`symfunnel.py --ns 27,52,65,82,101,125 --k 12`, ~0.3 CPU-h): record loosened x 1.05, Gaussian kick sigma on
x, y (20 sigma deg on angles; sym arm's kick symmetrised, so its effective sigma is ~sigma / sqrt 2), harden r0 0.1 +
polish.  Returns to the record (of 72 per sigma):
| sigma | free | sym |
|---|---|---|
| 0.05 | 18 | 60 |
| 0.1 | 15 | 58 |
| 0.2 | 8 | 53 |
| 0.4 | 0 | 8 |
Total 41 / 288 vs 179 / 288.  At matched effective kick (sym sigma vs free sigma / sqrt 2) still ~4x up to 0.2,
narrowing at 0.4.  Per n: 27 52 65 82 101 strong; 125 weak (sym 9 / 48; its record is Md only, not D4).  Free fails
even at 0.05 because the soft stage itself walks away from the record; the subspace pins it.
Reading: in the Md subspace the record funnels are several times wider; the from-scratch failure is the starts, not
the subspace.  Next: (a) symmetric basin hopping at Md-record n (the landscape has fewer minima; wider funnels mean
longer hops stay productive), released by the free polish; (b) time travel: symmetrise the pre-record packing of an
Md record and search the subspace -> the later record?

## TT-S1: time travel in the subspace, pre-registered (10-18 ~21:00; Evan: "Let's try it!")
`symtt.py run --out runs/symsub/tt1 --hops 60 --reps 4 --procs 15`.  Starts: the (up to) 2 best distinct public
packings above today's record at each Md-record n (store: kingbird revisions, older issues): 20 starts over 15 n (67 89
104 109 147 171 174 175 198 201 202 203 233 265), 0.003-0.05 above the record.  No record information: each start's own
best diagonal mirror (Hungarian mismatch 0-0.34; Ma rotated to Md), symmetrised by averaging matched pairs.  Arms (4
reps x 60 hops each, same seeds): sym = symmetric kick + harden --sym-md + free polish, monotone on the released side;
free = the same from the unsymmetrised start with free kicks / dynamics.  sigma cycles 0.05 / 0.1 / 0.2.  ~10 CPU-h.
Metric: per start, any rep reaching today's record (released side within 1e-7, or below); best gap; starts improved.
Rule: sym reaches the record on >= 2x as many starts as free AND on >= 4 starts -> stage (1), symmetric basin hopping
from current records at all 38 Md n (pre-registered separately).  Predictions: sym reaches the record on >= 4 starts
45 %; free on >= 4 starts 20 %; any start below today's record 3 %.
* **TT-S1 result (10-18 ~21:45; ~37 min wall on 15 procs; the per-run cpu field counts Python only, see manifest)**:
  0 / 20 starts reach today's record in either arm (0 / 80 runs each).  Free improved on its start at 3 starts, sym at 1
  (198: +0.0201 -> +0.0199).  Sym is usually far worse than its start: symmetrising a start with 0.12-0.34 mirror
  mismatch lands +0.02 .. +0.83 above it, and hopping doesn't recover.  On the starts that were already exactly Md
  (171 / 212229, 203, 233, 265) sym returns to the start every time (203 / 233 / 265: +0.0503, all 4 reps) and never
  goes lower; free drifts above them (+0.11 .. +0.13).  Predictions: sym on >= 4 starts 45 %: miss; free on >= 4 20 %:
  miss; below the record 3 %: none.  Rule not met: stage (1) (symmetric hopping at all 38 Md n) does not fire.
  Reading: the wide funnels of FW1 cut both ways.  The subspace holds whatever symmetric basin it is in, the record's
  or not, so it is a conservative tool (refine within a design), not an explorer across designs; and blind
  symmetrisation of an asymmetric packing is a big jump, not a projection.  Both sym and free hop operators (harden
  from a 1.03 loosen + kick) are too destructive to hold even their own start; a gentler kick ladder would be the next
  knob if this line is reopened.
* **DUP1, blend vs duplicate-half (10-18 ~22:10; Evan: "blend them, or duplicate one half and see if that relaxes?";
  `symdup.py`, ~0.3 CPU-h, not pre-registered)**: per TT-S1 start: blend (average), dup-lo / dup-hi (one half +
  its mirror, near-diagonal squares from the kept side snapped onto the diagonal, count fixed to n), one symmetric
  soft-to-rigid relax (r0 0.05 / 0.1 / 0.2 x 2 seeds) + free polish.  Best per start: blend 7, dup-lo 4, dup-hi 6
  (ties counted for each); typical released side +0.1 .. +0.6 above the record for both; 0 records; never better than
  the start except where the start was already symmetric (kept exactly).  Dup rescues a few blends (109: +0.050 vs
  +0.855; 198 / 212265: +0.020 vs +0.257).  The relax itself preserves symmetric inputs (171 / 212229 exact), so the
  damage is the seam + count fix: these older packings' halves are not mirror-compatible.
* **DUP2, pre-registered (10-18 ~22:20; Evan: "try both, in that order")**: `symdup2.py --out runs/symsub/dup2 --top
  60`: duplicate-half of the CURRENT record at the 60 asymmetric non-integer records (tilted, not Md / Ma invariant)
  closest to a diagonal mirror (Hungarian mismatch), both halves, symmetric relax r0 0.05 / 0.1 / 0.2 x 2 seeds + free
  polish (720 runs, ~1 CPU-h).  Metric: released side vs the frontier; below -> candidate.py (nothing posted).
  Predictions: a record 2 %; a tie (|d| < 1e-7, i.e. the record has a symmetric relative) at >= 1 n 25 %; median best
  gap > +0.05 70 %.
* **DUP2 result (10-18 ~22:35; 720 runs, ~1 CPU-h)**: 0 records.  Ties at the record side at 51 (mismatch 0.000),
  153 (0.018), 130 (0.006): the three most nearly symmetric records, so most likely the record with rattlers moved (not
  checked for distinctness).  Nearest non-ties: 172 +7.2e-4, 297 +8.3e-4, 207 +1.5e-3, 110 +1.6e-3, 308 +1.7e-3; gaps
  grow with mismatch (> 0.1 mismatch: +0.06 .. +0.6).  Predictions: record 2 %: none; tie at >= 1 n 25 %: hit (likely
  trivial); median best gap > +0.05 70 %: miss (median ~+0.048).  Line closed for now: the subspace refines within a
  design but every way of entering it from a non-symmetric packing (blend, duplicate-half, from scratch) costs more
  than it gives back at these n.

## DUP3: duplicate-half as a start distribution, pre-registered (10-18 ~23:00; Evan: "reduce the n, reduce compute,
see if this framing gets any interesting results at all")
`dup3.py gen / refine / report` (runs/symsub/dup3).  n = 51 110 153 172 207 297 (DUP2's ties / near misses).  Gen: 64
variants per n (written as 120 before launch; the grid is 64: seed 2 only at w = 0) (mirror Md / Ma x half x seam band w 0 / 0.7 / 1.5 (`harden --sym-band`: band deleted, refilled at
random, left free) x r0 0.05 / 0.15 x box loose / tight (held at the record side, overlaps persist) (+ seed 2 at w 0)),
symmetric soft stage + free polish.  Filter: within 5e-3 of the record, matching distance >= 0.005, distinct sides,
<= 6 per n.  Refine: explore.py 5 min each; baseline: the same number of 5-min explorer runs from the record itself.
~5 CPU-h.  Metrics: per arm, distinct polished minima within 1e-4 of the record (sides 1e-9 apart), how many at or
below the record, best refined side.  Predictions: dup finds more distinct near-record minima than base 55 %; any
minimum below a record 3 %; dup's kept starts number >= 3 at >= 4 of 6 n 60 %.
* **DUP3 result (10-18 ~23:30; gen 384 variants 95 s wall ~0.4 CPU-h; refine 26 x 5 min ~2.2 CPU-h)**: kept new starts
  (within 5e-3, md >= 0.005): 51: 1, 110: 5, 153: 2, 172: 0, 207: 3, 297: 2; every kept start came from w = 0 (seam
  bands never within 5e-3); best non-tie start 110 tight-box +1.1e-4.  After 5 min refinement: no record.  Dup starts
  reach the exact record side at 51 / 153 / 207 as DISTINCT packings (md 0.007-0.017; base from the record stays md
  0.000), at 110 collapse onto the record itself (md 0.001, 3 / 5), and 297 refines +8.8e-4 -> +1.8e-4 (md 0.014, new).
  Distinct minima within 1e-4: dup 14, base 18.  Predictions: dup more near-record minima 55 %: miss; below a record
  3 %: none; >= 3 kept at >= 4 n 60 %: miss (2 n).  Reading: as a start distribution, duplicate-half lands in the
  record's own funnel or on its exact-side plateau (like the shove ties), with one new near miss (297, +1.8e-4).
