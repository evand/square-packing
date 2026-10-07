# Alternatives for the rigid s(149) Göbel-type packing (2026-10-04)

Exploratory; cited from the Open problems page (§6, "How many alternatives?").

**Model** (`model.py`).  Side S = K + r, r = √2/2.  Axis-parallel unit squares sit at coordinates in
ℤ or ℤ + r (index t ↦ t//2 + r·(t%2)); in index space two squares conflict iff both index gaps are ≤ 1
(king graph).  Rotated squares appear only inside s(5) blocks (4 axis squares + one 45° square, side 2 + r),
treated as solid boxes worth 5.  Exact ILP via HiGHS.

**Results**
* K = 12: model max = 149, and Ellsworth's rigid layout (blocks at (0,0),(3,2),(5,5),(8,7),(10,10)) reaches it.
* `complete.py`: every 149-layout has exactly 5 blocks at integer corners whose x-corners form a chain
  0 → 10 with steps {2,2,3,3} (6 chains), same for y (by symmetry).  Cases ≤4 blocks: 148; ≥6: 145;
  odd corner: 147; shared x: 147; other x-sets: 148.
* `enum.py`: candidates = 6 x-chains × 6 y-chains × 120 pairings (overlapping ones skipped).
  Stopped about halfway (`partial_149.jsonl`, 2277 checked): 473 reach 149, and in all 473 the
  axis squares are uniquely determined in the model (suggesting rigidity; not proved).  No pairwise
  neighbour-offset rule separates 149 from 148 layouts.
* K = 6 (`run.py 6 38`): max 38, 54 discrete layouts with 2 blocks that have room to slide, so many
  of them are one continuous component.

**Not shown**: that distinct layouts are distinct components (plausible: moving a block needs slack ~1),
rigidity in the continuous sense, the full count, or symmetry classes (divide by up to 8).
