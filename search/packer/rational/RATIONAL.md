# Rational constructions (10-18; Evan: "a family of exact constructions, compare side lengths, look for records or
near records, then optimize with those as start points")

## What the rational records are (census, ~0 CPU)
Store public bests n <= 2000 whose side is rational with denominator <= 200 (to 2e-9): only 6.
| n | s | tilt | tilted part |
|---|---|---|---|
| 50 | 53/7 | 3-4-5 (36.87) | one sheared 4x4-ish island (16) |
| 171, 198 | 95/7, 102/7 | 3-4-5 | two 16-square islands |
| 230 | 643/41 | 20-21-29 (46.40) | a wide diamond block (104) |
| 293, 328 | 723/41, 764/41 | 20-21-29 (43.60) | a 3 x 20 strip (rows sheared); 328 = L-step of 293 |
Denominators are a + b (7 = 3 + 4, 41 = 20 + 21).  fillfast reproduces every upright count from the tilted squares
alone at the exact side.

## Tool: fixed-angle LP (`fixlp.py`)
Fixed angles => non-overlap is linear per separating-axis choice => min s is an LP; iterate separators.  Exact rational
optimum per structure at Pythagorean angles.  Check: loosened (x 1.03 + jitter) copies of 50 / 171 / 293 return to
53/7, 95/7, 723/41 to 1e-15 (ms each); 230 stopped at 644/41 (one separator basin up).

## RG1 pilot, pre-registered (10-18 ~17:20)
`ratgen.py run --out runs/rational/rg1 --K 12-15 --triples 3,4,5:20,21,29:8,15,17:5,12,13 --grid 4 --procs 15`:
221,184 designs (angle and its complement, C = 1..6 rows x L = C..1.45(K+1), unsheared, placement 4 x 4, start
side K + {0.3, 0.55, 0.8, 0.97}; fill -> LP -> refill loop).  ~3 CPU-h.  Metric: per n, best s* vs frontier (store
public best / claims); rows below, ties, within 1e-3.  Predictions: a record < 3 % (unsheared single blocks are the
simplest family, Ellsworth / Göbel-type ground is well trodden); any non-integer n where the best design is within
1e-3 of the frontier: 50 %; rediscovers 723/41 at 293 exactly: 30 % (needs the right shear).
* **RG1 result (10-18 ~17:20; ~14 min wall on 15 procs, ~3 CPU-h)**: 0 records, 0 ties at non-integer n.  Best per n
  over 122-232 (`runs/rational/rg1/report.txt`): closest 174 at 562/41 (+2.1e-4; one 1 x 18 strip at 43.60 vs the
  Göbel 13 + 1/sqrt 2), then 145 / 197 / 226 at +1.05e-3 (5 x 5 / 5 x 12 blocks at 20-21-29 vs Göbel-type 45-degree
  records): near-45 Pythagorean angles imitate 45-degree designs at a small loss.  Everything else +6e-3 .. +0.3;
  3-4-5 / 8-15-17 / 5-12-13 single blocks never get below the integer side except as above.  Predictions: record
  (< 3 %): none; within 1e-3 at a non-integer n (50 %): hit (174); 723/41 at 293 (30 %): not testable, my error (K
  12-15 excludes 293).  Reading: unsheared single blocks from a loose fill are too weak a family; positive control
  next (sheared 3 x 20 rows at 293).
* **Positive control PC293 (10-18 ~17:45; ~0.7 CPU-h)**: ratgen gained per-row shears.  20-21-29, 3 x 20, row shears
  (0, a, b) with a, b in {-1 .. 1 step 1/4}, placement 6 x 6, s0 = 17 + {0.7, 0.85, 0.97}: 8,748 designs.  **293 at
  723/41 exactly** (shear (0, 3/4, 1/4)), i.e. the record; 289-292 also at 723/41 (record family's subsets), 294 at
  17.893.  So shear is the missing ingredient for RG1's family; the generator + fill + LP reach a known rational record
  from scratch.
* **Store (10-18 ~18:00)**: `ratgen.py import DIR` puts designs within `--gap` (0.03) of the frontier into the store:
  source `ours:rational` (private), kind `rational-design`, tags `rational` + `tri:a-b-c`, `s_exact` when the side is a
  fraction with denominator <= 1e5.  RG1's per-n bests imported (only the best per n was kept there).

## RG2, pre-registered (10-18 ~18:05; Evan: "go ahead with RG2 or a related plan")
`ratgen.py run2 --out runs/rational/rg2 --K 10-17 --samples 120 --procs 15`: triples 20-21-29, 48-55-73, 65-72-97,
119-120-169, 3-4-5, 28-45-53 (one orientation: the mirror covers the other; RG1 wasted half on it), C = 2, 3, 4 rows,
L = max(C, K - 4) .. 1.45 (K + 1), per (K, triple, C, L) 120 random (row shears on the 1/4 grid in [-1, 1), placement,
start side K + {0.6, 0.8, 0.97}): 207,360 designs, ~10-15 CPU-h.  Keeps every non-integer design within 0.03 of the
frontier (<= 10 sides x 3 per n) in keep.jsonl, imported to the store after.  Metric: per n, best s* - frontier; count
of n within 1e-3; records.  Then (pre-registered here, run if RG2 gives >= 5 n within 1e-3): explorer from the best
rational design per such n (harden_explore, 30 min x 1 proc).  Predictions: a record 4 %; >= 5 n within 1e-3 of the
frontier 45 %; 723/41 at 293 and 643/41 at 230 rediscovered (both in range) 60 % / 25 %.
* **Angle release (`release.py`, 10-18 ~18:20)**: fq POLISH_DEEP (free angles) from each kept design; released packing
  imported (source `ours:rational-release`, parent = the design, tag `rational-release`).  Test on RG1's four near
  misses (145, 174, 197, 226; 2 procs, ~5 min): every one polishes exactly onto the current 45-degree record (d
  +3e-12 .. +3e-11).  So RG1's near misses are rational shadows of the Göbel-type records, not new basins; the release
  step is the test that separates the two (a released design that lands on a *different* minimum is the interesting
  kind).
* **RG2 result (10-18 ~19:00; 5.1 CPU-h measured, 207,360 designs)**: 0 records, 0 n within 1e-3 (best +2.5e-3 at 93,
  +3.2e-3 at 110, +5.6e-3 at 299 / 199).  **Positive control failed**: 293 only 4217/239 (+1.0e-2), not 723/41; 230
  +0.146.  Predictions: record 4 %: miss; >= 5 n within 1e-3 45 %: miss (0); 723/41 at 293 60 %: miss; 643/41 at 230
  25 %: miss.  The explorer follow-up (rule: >= 5 n within 1e-3) does not fire.  57 designs imported.
* **Why (PC293b / c / d, ~1 CPU-h)**: the 293 cell (20-21-29, 3 x 20) on the PC293 grid hits 723/41 at n >= 293 in
  **1 of 8,748** designs (shear (0, 3/4, 1/4), block flush on the left wall).  RG2's 120 random samples per cell were
  hopeless.  Stress-removal ladder (`ratgen.ladder`: drop the square with the largest LP contact force (fixlp duals),
  polish, refill; best per n along the way): 0 / 600 random designs (no help).  Wall-biased placements (fx, fy in
  {0, 1} half the time): 1 / 600, ~15x better, still far too sparse (~1.3 CPU-s per design with the ladder; a full scan
  ~400 CPU-h).  Reading: the loose-side fill decides the count and almost never matches the tight structure; the
  family needs a count-aware construction (fill at the target side, or contact-built staircases) rather than sampling.
* **REL1 (release of all 118 kept designs, ~0.5 CPU-h)**: 0 below; released minima tie the record at 145, 150, 174, 197,
  226, 293 (the 45-degree / 723/41 records themselves); every other non-integer n stays >= 2.1e-3 above (199, 110,
  205, ...).  118 released packings in the store (`ours:rational-release`, parents linked).  Nit: RG1's import also
  took integer-side (grid) designs; `import` should skip those.
