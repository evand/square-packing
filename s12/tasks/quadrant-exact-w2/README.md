# quadrant-exact-w2: exact certificate for the (R = 2, w = 2) fixed-profile family → s(k²−3) = k for all k ≥ 6 (2026-09-28)

**Context.**  Read `search/QUADRANT.md` in full (esp. §1 accounting/localisation/seam bound, §5 stability, §6 structure,
§8 checker scope).  Float LP: D = 0.945 ± 0.002 for R = 2, w = 2 (runs in `runs/quad_*`); D > 3/4 exactly would give
s(k²−3) = k for every k ≥ 6 (Bentz's conjecture, with the known k ≤ 7).  The float solution still has captures down to
≈ 0.9986; the job is to turn it into an **exactly valid** family with D > 3/4, certified.

**Hard constraint.**  No ×f scaling of the band or interior (σ < 0 kills the all-k claim).  σ = 0 exactly; Lebesgue
interior exactly.  Slack may be bought only in the corner module (and D has 0.195 of room over 3/4 to spend).

**Steps.**
1. **Re-solve for an exact-friendly solution**: add the "Lebesgue layer" constraint (§8.3b: module and band equal to
   Lebesgue in a thin layer along their boundary with L, so slivers there are trivially exact), refine pitch near the
   seam, get a **vertex** (crossover) solution, keep closing rows until the float oracle's min capture is ≥ 1 − 1e-9 or
   the remaining violators are identified continuum-tight families.  Report D at each stage.
2. **Exact rational solution**: re-solve on the active rows in exact arithmetic (`search/dual_exact.py` style or an
   exact simplex on the basis), with σ = 0 kept as an equality.  Masses as fractions.  Then, if needed, trade D for
   margin in the corner module only (e.g. add small corner mass) and record the exact D.
3. **Exact check of one box.**  By §8.1 the box k = 2R + 3 = 7 certifies every k ≥ 6.  Build `μ_7` element by element
   (`quadrant_tools.box_saving`), write it as a mixed cover; total must be exactly 49 − 4D < 46 (needs D > 3/4).
   Polygon (Lebesgue) masses: `zm_mixed.py` has polygon code disabled in `--cert-mode`.  Do **not** edit the pinned
   `search/zm_mixed.py`; copy to `search/zm_poly.py` (or similar), enable polygons, audit the polygon lower bound
   (area of square ∩ polygon over a pose box) and write down its lemma; rejection tests (a cover that must fail).
   Zero-margin continuum-tight families (squares sliding along the wall; Lebesgue-boundary slivers; germs = Lemma T)
   need exact lemmas; state each, prove it on paper, implement, test.  zmx2 (Rust, `verify2/`) as the fast float-ish
   second opinion if it can take polygons; otherwise note it.
4. **Double-check the reduction** (quadrant ⇒ all k ≥ 2R+2; box 7 ⇒ quadrant, §8.1) independently of QUADRANT.md:
   re-derive it, and test it numerically by building μ_k for k = 6..12 and running the float oracle on each.

**Labels** [proved]/[measured]/[heuristic].  An exact result must say exactly which program, which sha, which lemmas.
**Compute:** physical cores 0–11 (≤ 12 processes), ≤ 60 GB.  **Don't** edit TODO.md, README.md, pinned checkers
(`zm_mixed.py`, `zeromargin.py`, `mixed_cover.py`, verify2 sources — add new files instead), certificates, site, lean.
**Deliverables:** `search/QUADRANT_EXACT.md`, code in `search/`, runs `runs/qx2_*`; commit locally on main in
`/home/evand/math/square-packing/public` (no push).  Report every ~6 h at most: stage reached, exact D, open lemmas.
