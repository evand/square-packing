# Packing similarity and clusters at n = 110 (10-14)

`analysis/similarity110.py` (stages `features pairs reps edges analyze othern`; ~15 CPU-min at 8 procs).  Outputs:
`runs/similarity110/` (`report.txt` = every table, `summary.json`, `othern_s*.json`).  Store opened read-only.
Ground truth: 939 exact minima at 110 (`minimum`), 36,317 store rows labelled through `copy`.  Sample: 939
representatives + 4,112 copies (half of them unconverged, side > S* + 1e-9, chosen on purpose as the hard "same" cases).
Pairs: 13,643 same-minimum, 15,132 "near" (each packing's 3 nearest-in-side packings of another minimum), 6,000 random.
md = matching distance (`rediscover.distance`: RMS centre + angle displacement, best of D4 x assignment); mv = fraction of
squares moved > 0.1 under that matching; img = `pk/fingerprint` distance.

## Numbers
**Same minimum (task 2).**  Rates: false-split = same pairs called different; false-merge = different pairs called same.

| test | false-split | false-merge (near side) | false-merge (random) |
|---|---|---|---|
| \|ds\| < 2e-9 (stored side) | 0.113 | 0.022 | 0 |
| \|ds\| < 1e-7 | 0.069 | 0.136 | 0 |
| \|dsp\| < 2e-9 (after fq polish) | 0.092 | 0.022 | 0 |
| \|dsp\| < 1e-8 and img < 0.02 | 0.067 | 0.046 | 0 |
| img < 0.01 | 0.013 | 0.40 | 0.020 |
| md < 0.05 | 0.0001 | 0.52 | 0.16 |
| census3 tilt signature equal | 0.14 | 0.25 | 0.004 |
| fq contact fp equal (stored labels) | 0.999 | 0.0001 | 0 |
| fq contact fp equal (copy relabelled onto rep by the matching) | 0.999 | 0.0001 | - |
| label-free contact-graph WL hash (gap < 1e-8, rattlers stripped) | 0.91 | 0.008 | 0 |

- **Contact fingerprints are not identity keys at 110.**  fq's fp is deterministic, but a 1e-12 coordinate jitter
  changes it in 6 / 8 representatives (load-bearing set 263-318 contacts, 5-11 contacts change): the minima are
  degenerate (flat motions, weak contacts).  Labels are not the problem (relabelled: same 0.999).
- **Distinct minima with near-equal sides are geometric twins.**  12 consecutive exact minima are < 1e-9 apart, 43 < 1e-8:
  median md 0.002, img 0.001, no square moves > 0.1.  No geometric measure separates them and search need not.  The
  side-key false merges (2.2 % of near pairs) are these 28 twin pairs.
- Side false-splits are unconverged copies (s up to 1e-3 above S*; md to their minimum p50 0.003, p99 0.035).  fq
  polish converges about a third of them (copies > 1e-9 above S*: 10.7 % -> 7.5 %; split 0.113 -> 0.092).  Within-minimum md p99 0.035, max 0.059; different-minimum md
  nearest-in-side p50 0.047, random p50 0.069: geometry overlaps, side does not.

**Regions (task 3), 939 representatives.**  Basins reconstructed as in `landscape_depth.py` (A = 11, B = 39).
Transitions = 10,590 explorer parent -> child links between two different known minima (46 archives: ex*, census v2-v4,
c110_h*).  "within" = share of those links inside one cluster; null = same with shuffled endpoints.

| clustering | #cl | n_eff | largest | A one cluster, A / B apart | within | null | lift | Spearman(S, cluster mean S) |
|---|---|---|---|---|---|---|---|---|
| single-link img 0.02 (landscape/ islands) | 39 | 17 | 141 | yes / yes | 0.83 | 0.12 | 7.2 | 0.81 |
| single-link img 0.025 | 17 | 3.6 | **642** | yes / **no** | 0.93 | 0.65 | 1.4 | 0.57 |
| single-link md 0.035 | 30 | 18 | 132 | yes / yes | 0.85 | 0.10 | 8.5 | 0.84 |
| single-link md 0.04 | 15 | 2.7 | **717** | yes / **no** | 0.96 | 0.75 | 1.3 | 0.50 |
| **leader md 0.04** (deepest first) | 40 | 27 | 97 | yes / yes | 0.72 | 0.06 | 11.8 | 0.88 |
| leader md 0.03 | 53 | 36 | 93 | yes / yes | 0.67 | 0.05 | 13.0 | 0.90 |
| leader md 0.05 | 17 | 10 | 252 | yes / yes | 0.80 | 0.17 | 4.8 | 0.53 |
| leader md + 10 \|dS\| 0.05 | 37 | 24 | 90 | yes / yes | 0.70 | 0.07 | 10.8 | 0.91 |
| leader img 0.03 | 32 | 19 | 115 | yes / yes | 0.74 | 0.08 | 9.0 | 0.82 |

- Single linkage percolates: one step past the island cut (img 0.02 -> 0.025, md 0.035 -> 0.04) a giant cluster takes
  70-75 % and merges A with B.  It will get worse as the store fills in.  Leader balls do not percolate, and they keep A whole and
  A / B apart at every radius 0.01-0.05.
- Reachability: P(a pair of minima is ever linked by one explorer step) = 0.47 at md <= 0.003, 0.20 at <= 0.017, 0.03
  at <= 0.044, 0.003 at <= 0.068.  Transition md p50 / p90 / p99 = 0.010 / 0.047 / 0.077 (all pairs p50 0.068).  So
  one step reaches ~0.04-0.05 (C02), and a leader ball of 0.04 holds 72 % of steps (11.8x chance).
- Side term (Evan's "md + side gap"): md + 10 |dS| at 0.05 is about the same as md at 0.04 (lift 10.8 vs 11.8, a bit more
  depth-sorted).  At 110 the side term adds nothing; at larger n it would dominate (sub-k sides span ~0.1), so leave it
  out of the cluster and use depth inside the cluster instead.
- Leader order: random visit orders give 35-37 clusters vs 40, ARI 0.69-0.82 against deepest-first.  Use deepest-first.

**Other n (task 4; side key as the minimum proxy, 400 random sub-k side keys each; 88 also at 1,600).**

| n | sub-k side keys | parent -> child md p10/50/90 | NN md p50 (in sample) | all-pairs md p50 | leader md 0.04: clusters / sample | img per md (at md ~0.04) |
|---|---|---|---|---|---|---|
| 88 | 9,414 | 0.018 / 0.042 / 0.095 | 0.061 (0.040 at 1,600) | 0.21 | 328 / 400; 1,045 / 1,600 | 1.25 |
| 110 | 4,176 | 0.003 / 0.008 / 0.047 | 0.002 | 0.061 | 28 / 400 | 0.75 |
| 237 | 7,113 | 0.010 / 0.031 / 0.108 | 0.061 | 0.26 | 319 / 400 | 0.65 |

- 110 is a dense funnel.  The stored sub-k rows at 88 and 237 are spread out, so clusters of radius 0.04 there are mostly 1-3 minima.  The
  img : md ratio changes 2x across n (Spearman 0.89-0.97), so an img threshold cannot be carried across n.  Used as a
  prefilter, img < 2r still catches 99.7-100 % of md < r pairs at all three n.
- Step size falls with n (parent -> child p50 0.042 at 88, 0.031 at 237; ratio 1.4, sqrt(237/88) = 1.6), as expected
  for RMS over n when a step moves a fixed number of squares.  Hence the radius below is on the total displacement
  sqrt(n) * md.  Evidence is two n; treat the exponent as a default, not a measurement.

## Recommendation
- **Same packing (cheap key):** `side_key` to 2e-9 on a polished row, plus a geometric check that only guards against
  collisions: same iff |ds| < 2e-9 and md < 0.01 (img < 0.01 as the cheap stand-in).  False-merge 1.7-1.9 % of
  near-side pairs, all of them geometric twins (md ~0.002).  Do not use any contact fp as a key.  Unconverged copies
  stay split until polished (exactsolve is the only full fix); they fall into their minimum's cluster anyway.
- **Cluster:** leader balls on matching distance, deepest first, radius on total displacement
  sqrt(n) * md <= 0.42 (= md 0.040 at 110, 0.045 at 88, 0.027 at 237).  Prefilter leaders with img < 2r, then exact
  md on the survivors (~2.7 ms per pair at 110, ~10 ms at 237).  Sample a cluster uniformly, then a minimum (a side key,
  not a row) uniformly inside it.
