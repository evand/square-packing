# E. Zero-margin verification at `t = m`

**Why.** Every argument that ends at `s(12) = 4` (or `s(45) = 7`) is a certificate at the container
`[0,m]²` itself, closed semantics, with tight poses of coverage exactly 1 — the axis-parallel
corner and wall squares and their tilted neighbours.  Erosion verifiers (ours, Mira's subdivision)
cannot certify any closed cover of `[0,4]²` with `W < 16`: they check eroded squares, and 16 eroded
unit squares fit.  The literature handles these poses with non-avoidance lemmas (Bentz 2010
Lemmas 1–7, Nagamochi 2005 Lemmas 2–7, DS7 Lemmas 1–7; all in `notes/proof-anatomy.md`).  This task
mechanises that.

**Design to evaluate.**
- Subdivide pose space `(cx, cy, θ)` adaptively (Mira-style, exact rationals); a box is certified
  when its core captures weight `≥ 1` (points) — this handles everything except neighbourhoods of
  tight poses.
- **Lemma primitives** certify the rest: a box whose centres all lie in the corner box `[0,1]²` is
  certified by "contains `(1,1)`" (DS7 Lemma 1); wall-strip boxes by DS7 Lemma 2 / Bentz Cor. 3
  ("contains `(1,y)` or `(1+x,y)`", so `min(w_A, w_B)` is captured); triangles by Lemma 3.  Each
  primitive is a container-independent, rotation-general statement with a one-parameter
  trigonometric proof; list exactly which are needed for the calibration cover, state them, and
  check each numerically over a fine grid as a sanity test (formal proofs later).
- Where the cover's tightness comes from *interior* grid-line points (the `closed4.py` optimum puts
  91 % of its weight on `x, y ∈ {1,2,3}`), the capture region of a boundary point near a tight pose
  is a cusp (width `~ε²` in the centre for tilt `ε`), and no lemma helps; note this, and consider
  **segment resources** (weight = density × chord length; rational for rational rotations,
  piecewise linear in the centre, so exact on arrangement cells) as the robust replacement.

**Calibration ladder.**
1. Friedman's 14 points at `[0,4]²`, unit weights (DS7 Thm 4; coordinates in `proof-anatomy.md §5.3`):
   certify exactly that every closed unit square in `[0,4]²` contains one — `s(15) = 4`.  Slack 1.
2. A closed cover of `[0,4]²` with `W < 13` (`closed4.py` finds ~12.5 heuristically, but it may lean
   on grid lines; a robust one may need Nagamochi-style offsets): `s(13) = 4` re-proved by machine.
3. (Later, with task D) `[0,7]²`, `W < 45`.

**Deliverables.** A design note `search/ZEROMARGIN.md` (what is certified how, which primitives,
where the cusps are); a prototype exact checker (Python, rationals) that certifies rung 1; an
honest statement of what rung 2 needs.  Light compute.  Do not touch `TODO.md` or other tasks' files.

## Status (2026-08-29)

**Rung 1 done.**  `search/zeromargin.py friedman14 --tri` certifies Friedman's 14 points exactly for
closed unit squares in the closed `[0,4]²` at every angle (`s(15) = 4`): 6,958 boxes, depth 10,
3,356 leaves by exact bin cores, 70 by the 2×2-box lemma (walls/corners), 74 by the triangle lemma,
0 uncertified; also on the unreduced pose space.  Leaf dump `runs/zeromargin_friedman14_leaves.txt`,
independent float stress test `search/zeromargin_stress.py` (0 failures on 6,679 leaves × 40 poses
and on 340k random primitive instances).  Design note: `search/ZEROMARGIN.md`.

Findings: (a) the exact bin core (`R_θ₀Q ∩ R_θ₁Q ∩` the sector-arc condition) is what makes
axis-parallel tight families certifiable by boxes — the σ-square is not; (b) walls and corners need
the one-line 2×2-box lemma (quadratic margin); (c) Friedman's set has a *tilted* tight family at
`θ = arctan(3/4)` from the unit-distance pair `(1.6,1)–(1,1.8)`, which no box argument handles and
which needs the triangle lemma — the checker needs exactly Friedman's Lemma 3 and nothing else.

**Rung 2 assessed, not done.**  The best closed cover at `s = 4` (`runs/closed4_best.txt`, 12.42)
is robust at the grid-pose fans (1.04–1.11) and exactly P1 at the corner; its failures are
unconverged tilted wall poses (0.9927), so rung 2 is a separation loop at `s = 4` with this
checker as the oracle, plus the open question whether the converged cover has tight tilted
families (then a two-region primitive is needed, not built).  Segments buy nothing.
