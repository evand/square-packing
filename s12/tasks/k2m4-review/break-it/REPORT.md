# break-it report: K4_k008_box9 (Valid9), 2026-10-03

Reviewer: adversarial agent, cores 4-5, about 1.5 CPU-h.  Everything here is in this directory, and all of it is my own code: an
independent parser and evaluator (`cover.py`: exact `Fraction`; `fastf.py`: vectorised float, agrees with exact to 1e-14 on 300
random poses).  I read none of qx2_zm/germscan's logic except the line-profile construction (see Minor).

## Findings

**BREAKS: none.**  No pose with exact mass < 1 was found.

**GAP: none found in my area.**

**MINOR**
1. **Large exactly-tight set (by design, not a bug).**  Density-1 Lebesgue on [14/5,31/5]^2 means every unit square inside it that
   meets no segment has mass exactly 1, at every angle.  Examples: the corner square [3,4]^2 with tiny tilt, and [14/5,19/5]^2 tilted
   inside the box.  Any perturbation of D or of the Lebesgue weight downward breaks Valid9 on an open set, so the margin is exactly 0
   there.  It is harmless for the claim, but `81-4D` cannot be "rounded".
2. Smallest non-flat generic-angle local minimum found: 1.0000448 (float), near (4.904, 3.407, 44.5 deg), a square poking about 0.1
   below the Lebesgue edge y = 14/5.  The (m-1)/protrusion ratio near the 14/5 edge is >= 0.005 at protrusion about 0.1, and it is
   linear from 1 for tiny protrusions: line 14/5 carries segments for y in [0.6,8.4] except gaps [1.2,1.4] and [7.6,7.8].
3. Duplicated segments: the 16 doubles on x,y in {3,6} add, by FORMAT and in my evaluator.  `qx2_zm.py` `Profile.__init__`
   (`search/qx2_zm.py:269-281`) sums the densities of all covering segments, and `zm_mixed.py:79-104` sorts endpoints and appends the
   duplicates, so the doubles add there too.  Reversed endpoints (for example `38 4 37 4`) are normalised (`zm_mixed.py:83,86`).

## Checked, OK (how)
* Structure: 2076 unit segments (1038 H, 1038 V) on the 1/5-grid, plus one square polygon of density 1.  The measure is invariant
  under all 8 elements of D4 (exact multiset check).  Total 3835229774429/50000000000 = 76.7046 < 77.
* **Own exact germ scan** (`germ.py`), with my derivation: for theta -> 0+, bottom keeps x <= X and top keeps x >= X; left keeps
  y >= Y and right keeps y <= Y (mirrored for 0-).  The wall forces the pivot corner.  I scanned all 41^2 grid x0,y0, all X,Y on the
  1/5-grid, and both signs: 111,392 germs, **min exactly 1**, 10,824 tight, 0 below.  Bilinearity per cell means this also covers
  the whole theta = 0 face and its one-sided limits.
* Tiny tilts at every tight germ (`tight_tilt.py`): 2.55M float poses, theta in {1e-4..3e-2}, cut positions perturbed by +/-0.25,
  order-1 slides.  0 violations.
* Exact rational tilts (t = 1e-6, 1e-4, 1/300) on a 5x5 lattice of cut offsets, at all 4,120 tight germs touching lines
  {3, 6, 14/5, 31/5} (the doubled end lines, the corner/band junction and the Lebesgue edge) (`exact_tight.py`): 307k exact poses,
  0 below 1, min (m-1)/t = 0 (flat set only).
* Multi-scale near-axis search (`nearaxis.py`): 5.6M poses, offsets and tilt from 1e-8 to 1e-1.  The 68 float "violations" at
  theta about 1e-8 are float noise; all are > 1 exactly.
* Global search (`rand2.py`, `gridsearch.py`, `nmstarts.py`, `band*.py`, `side45.py`): about 6M random and grid poses plus about
  1,500 Nelder-Mead descents.  Every float value < 1 was rechecked exactly: each was either > 1 or outside [0,9]^2 (theta =
  pi/2 float artefacts at the wall).

## Not checked
* That qx2_zm's leaves/Lemma E are sound (area `leaves`), and the Lean side (`claim`).  My searches are sampling plus local descent:
  they are evidence, not proof, at generic angles.
* Second-order germs (offsets about theta^2) were only sampled, not enumerated exactly.
