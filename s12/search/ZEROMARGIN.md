# Zero-margin verification at the container itself (task E; 2026-08-29)

Code: `search/zeromargin.py` (exact checker), `search/zeromargin_stress.py` (independent float
stress test of a leaf dump and of the three primitives), `search/zeromargin_fan.py` (rung-2
diagnostic).  Everything load-bearing is `fractions.Fraction`; floats only pre-filter.

**Result (rung 1).**  Friedman's 14-point set is certified, exactly, to be unavoidable for closed
unit squares in the closed container `[0,4]²` at every angle — i.e. `s(15) = 4` re-proved by
machine — with 3,356 leaves certified by exact cores, 70 by the 2×2-box lemma at the walls and
corners, 74 by the triangle lemma, 3,179 empty (no admissible pose), 0 uncertified; max
subdivision depth 10, 6,958 boxes, under a second.  Without the triangle lemma the checker does
not terminate (87,200 uncertified boxes at depth 21), and the reason is a genuine one-parameter
family of tight poses at `θ = arctan(3/4)` that no erosion or core argument can handle (§4).

## 1. The statement and the semantics

`P ⊂ [0,m]²` finite, weights `w_p ≥ 0`.  Claim: for every closed unit square `Q` with
`Q ⊆ [0,m]²` (closed containment), at every angle, `∑_{p ∈ Q} w_p ≥ 1`, a point on `∂Q` counting.
Pose `(cx, cy, θ)`; `Q = c + R_θ [−½,½]²`; admissible iff `cx, cy ∈ [w/2, m − w/2]` with
`w = |cos θ| + |sin θ|`.  This is exactly the convention of `certificates/FORMAT.md` and of every
`s(m²−3)` proof in the literature (dilated boxes ⇔ closed unit squares; `notes/proof-anatomy.md
§1, §7.1`).

Why the existing verifier cannot do it: it checks the concentric σ-shrunk square over an angle
bin, and at the corner pose `([0,1]², θ = 0)` the only point captured, `(1,1)`, is on the boundary;
the tilted neighbours `(w/2, w/2, ε)` contain `(1,1)` with margin `(1−cos ε)(1−sin ε) ≈ ε²/2`,
while any bin of width `ε` erodes by `≈ ε`.  Quadratic margin against linear erosion: no `N`
works.  That is the "16 eroded squares fit" obstruction of `notes/TODO-archive-2026-08-28.md`.

## 2. The checker

Adaptive subdivision of pose space `(cx, cy, u)`, `u = tan(θ/2)`, into boxes with rational
endpoints (root grid of pitch `1/10` in `cx, cy`, eight `u`-bins on `[0, 1/2]`; the longest
scaled side is halved).  A box is a leaf when one of four tests succeeds; otherwise it is split,
down to a depth limit, below which it is reported as uncertified with its location.  Only the
*admissible* poses of a box matter.

**EMPTY.**  `cx₁ < w_lo/2`, or `cx₀ > m − w_lo/2`, or the same in `y`, with `w_lo = min(w(θ₀), w(θ₁))`
(`w` is unimodal on `[0, 90°]`, so this is the smallest `w` on the bin): no admissible pose.

**CORE (exact bin core).**  For a fixed centre, the intersection of the closed unit squares over
an angle bin is

    core(θ₀, θ₁) = R_θ₀ Q  ∩  R_θ₁ Q  ∩  { x : dir(x) mod 90° ∈ [θ₀, θ₁]  ⇒  |x| ≤ ½ }.

*Proof.*  `x = r e_φ ∈ R_θ Q` iff `r·max(|cos(θ−φ)|, |sin(θ−φ)|) ≤ ½`.  The function
`g(ψ) = max(|cos ψ|, |sin ψ|)` has period 90°, maxima 1 at `ψ ≡ 0`, minima at `ψ ≡ 45°`; on an
interval of length `< 90°` it is maximised at an endpoint unless the interval contains a point
`≡ 0 (mod 90°)`, where the maximum is 1.  ∎  This is **larger than the σ-square** exactly where
it matters: the midpoint of an edge, `(½, 0)`, lies in the core of every bin `[0, h]` (it has
`|x| = ½` and direction 0), whereas the σ-square loses it.  The core is convex, so
"`p ∈ Q(c, θ)` for all `c` in the centre rectangle and all `θ` in the bin" ⇔ the four corners
of `p − rect` lie in `core(θ₀, θ₁)`: two rotated-square tests (rational `cos, sin` from `u`) and
one sector test (fold `R_{−θ₀}(p − c)` into the first quadrant mod 90°; angle `≤ θ₁ − θ₀` ⇔
`q cos Δ ≤ p sin Δ` with `cos Δ, sin Δ` rational) plus `|v|² ≤ ¼`.  Weighted: sum `w_p` over the
points passing; leaf if `≥ 1`.  The centre rectangle is first clipped to the widest admissible
range on the bin (`cx, cy ∈ [w_lo/2, m − w_lo/2]`); the poses that this still over-tests are
inadmissible ones in a sliver of width `≈ (w(θ₁) − w(θ₀))/2`, which is sound and, at positive
margin, harmless.

**P1 (2×2-box lemma).**  *If `|c − p|_∞ ≤ 1 − w(θ)/2`, i.e. the unit square `Q(c, θ)` lies in
`p + [−1,1]²`, then `p ∈ Q(c, θ)`.*  Proof: with `d = p − c`, `|d|_∞ ≤ t := 1 − w/2`, the rotated
coordinates satisfy `|d_x cos θ + d_y sin θ| ≤ t w = w − w²/2 ≤ ½` because `(w − 1)² ≥ 0`. ∎
(Stromquist 1984-I Lemma 1 / DS7 Lemma 1 "a unit square inside a square of side 2 contains its
centre" is the case `p = (1,1)`, container `[0,2]²`.)  Applied to the admissible poses of a box:
each of the four inequalities is checked against the box with `w_hi = max w` on the bin
(`1.4143 > √2` if the bin contains 45°), *except* that a side on which `p` is within 1 of the
container wall is automatic (`p_x − 1 ≤ 0 ⇒ cx ≥ w/2 ≥ p_x − 1 + w/2`).  This is what makes P1
exact at the corner and wall poses: at `(w/2, w/2, ε)` equality holds for every `ε`, so the
quadratic-margin family is certified in one box.

**TRI (triangle lemma).**  *`T` a closed triangle with sides `≤ 1`; a closed unit square whose
centre lies in `T` contains a vertex of `T`* (Friedman / Stromquist, DS7 Lemma 3, Bentz 2010
Lemma 2, there stated for open boxes; the closed version follows by taking concentric boxes of
side `1 + ε → 0`, or directly:).  *Proof.*  Suppose no vertex is in `Q`; each vertex `v` then has
`|⟨v − c, e₁⟩| > ½` or `|⟨v − c, e₂⟩| > ½`.  Two vertices on opposite sides of the same slab are
`> 1` apart, so all vertices lie in `H₁⁺ ∪ H₂⁺` (say) and none in `H₁⁻ ∪ H₂⁻`.  Assign each vertex
to a type with `⟨v − c, e_i⟩ > ½`; write `c = ∑ λ_v v`.  Projecting on `e₁`:
`0 = ∑ λ_v ⟨v − c, e₁⟩ > λ_A/2 − λ_B/2`, so `λ_B > λ_A`; on `e₂`, `λ_A > λ_B`.  (If a type is empty
the other projection gives `0 > ½` directly.) ∎  Used as: centre rectangle `⊆ T` for a triangle of
points of weight `≥ 1` each (rung 1 has 10 such triangles) ⇒ leaf, any angle.

**Symmetry.**  If `P` is invariant under `x ↦ m − x` and `y ↦ m − y` (checked exactly; Friedman's
set is, but is *not* invariant under `x ↔ y`), the `x`-reflection sends `θ ↦ −θ ≡ 90° − θ`, so
`θ ∈ [0, 45°]` suffices, and the 180° rotation then gives `cy ≤ m/2`.  Root domain
`u ∈ [0, ½]` (θ up to 53°, an over-cover to keep endpoints rational), `cx ∈ [0, m]`, `cy ∈ [0, m/2]`.
`--full` runs the unreduced domain (`θ ∈ [0, 90°]`, `cy ∈ [0, m]`): also VERIFIED, 27,396 boxes.

## 3. Results on Friedman's 14 points (rung 1)

| run | boxes | depth | CORE | P1 | TRI | EMPTY | uncertified |
|---|---|---|---|---|---|---|---|
| reduced domain, `--tri`, depth ≤ 14 | 6,958 | 10 | 3,356 | 70 | 74 | 3,179 | **0** |
| full domain, `--tri` | 27,396 | 10 | 13,556 | 298 | 280 | 12,364 | **0** |
| root pitch 1/5, 5 bins, `--tri` | 1,560 | 10 | 714 | 22 | 44 | 500 | **0** |
| reduced, **no** `--tri`, depth ≤ 21 | 504,498 | 21 | 165,000 | 70 | 0 | 3,179 | 87,200 |

Leaf dump: `runs/zeromargin_friedman14_leaves.txt` (box, kind, witness).  Independent check
(`zeromargin_stress.py`, floats, no shared code): 40 random poses per leaf including box corners,
each witness verified directly — 0 failures over 6,679 leaves; the three primitives on 200,000
random instances (P1), 68,576 (core lemma, 21 angles per instance), 72,394 (triangle) — 0 failures.

## 4. Where the tight poses are, and what each needs

*Tight* = coverage exactly 1 (for unit weights: every captured point on `∂Q`).

1. **Corner squares** `(w/2, w/2, θ)`, all `θ`: the only point is `(1,1)`, margin `≈ θ²/2` —
   quadratic.  CORE never terminates (the admissibility boundary `cx = w(θ)/2` is a curve, and a
   bin of width `h` always contains the inadmissible corner `(w(θ₀)/2, ·, θ₁)`).  **P1** certifies
   the whole region `Q ⊆ [0,2]²` in one test.
2. **Axis-parallel wall squares** `(cx, ½, 0)`: a point of the row `y = 1` on the top edge.  The
   tilted neighbours keep a point with linear margin, but the inadmissible sliver of a bin kills
   CORE when the surviving point is to the *left* of the centre (`d = p_x − cx < −h/4`), which is
   again the "point on the boundary of the limit square" situation.  **P1** with the wall points
   `(1.6,1), (2.4,1), (3,1), (1,1)`: their `P1` intervals cover every `cx` for `w ≤ 1.2`
   (`θ ≤ 11.5°`); beyond that the wall poses near `cx = 2` have positive margin (`≥ 0.01`, the
   `x + 2y = 2.8 < 2√2` slack of DS7 Lemma 2) and CORE finishes them.
3. **Axis-parallel interior squares with points on their edges**, e.g. `(1.5, cy, 0)`,
   `cy ∈ [1.8, 2.2]` (four points on the vertical edges) or `(1.5, 1.5, 0)` (four on three edges):
   one-parameter families, but all at `θ = 0` with centre coordinates in `(1/10)ℤ`, hence on box
   boundaries.  Each adjacent box has a single point in its exact core (e.g. the box
   `[1.5,1.6] × [1.9,2.0] × [0,h]` is certified by `(2, 1.8)`, the box to its left by `(1, 2.2)`)
   because the exact core keeps edge-midpoint-like points that the σ-square drops.  No lemma.
   This depends on the tight pose lying on the box lattice: with a root pitch not commensurable
   with the coordinates the boxes straddle the family and CORE alone would not terminate.
4. **Tilted families from unit-distance pairs.**  `(1.6,1)–(1,1.8)` are exactly 1 apart in the
   direction `(−0.6, 0.8)`; a square at `θ = arctan(3/4) = 36.87°` (`u = 1/3`, rational) has them
   on two opposite edges and can slide along the edge direction `(0.8, 0.6)`: centres
   `(1.3, 1.4) + t (0.8, 0.6)`, and for `t ∈ [0.02, 0.3]` (between `(1,1)` leaving and `(2,1.8)`
   entering) no other point is inside.  A tight *curve in general position* — not on any box
   lattice, not at a bin endpoint — and the minimum-margin search (`/tmp`-style Nelder–Mead from
   3,000 random starts) finds exactly this: margin `0` at `(1.3227, 1.4170, 36.87°)`.  Mirror
   images: `(2.4,1)–(3,1.8)`, `(1,2.2)–(1.6,3)`, `(2.4,3)–(3,2.2)`.  The horizontal unit pairs
   `(1,1.8)–(2,1.8)` etc. give the `θ = 0` families of item 3.  **TRI** with `(1.6,1), (2,1.8),
   (1,1.8)` (sides 0.894, 1, 1) contains the whole tight segment in its interior and certifies
   every box around it for every angle.  This is the one place where a genuine non-avoidance
   lemma is unavoidable, and it is exactly Friedman's Lemma 3 — the checker rediscovers the
   proof's structure.

**Cusps.**  For a *point* `p` on the edge of a tight axis-parallel square, the set of nearby poses
capturing `p` is a cusp: `(1.5, 1.5 + δ, ε)` captures `(1.5, 1)` iff `δ ≤ ½(1/cos ε − 1) ≈ ε²/4`.
The exact core handles it because the cusp's boundary is the arc `|v| = ½` of the core, which is
in the core (closed), and the box on the other side (`δ ≥ 0`) is certified by a different point.
What no box argument handles is a tight set not aligned with the boxes (item 4) or with only
quadratic margin (item 1); those need a primitive.

## 5. Rung 2: a weighted closed cover of `[0,4]²` with `W < 13`

What the checker needs from a weighted cover: at every leaf, a *fixed* set of points, each in the
leaf's core (or P1-certified), with total weight `≥ 1`.  Near a tight pose the boxes on each
side converge to the one-sided limit poses (e.g. `(1.5+0, 1.5+0, 0)`, the square shifted by an
infinitesimal amount up and right, which captures the `y = 2` and `x = 2` edges but not `y = 1`,
`x = 1`), so the checker terminates iff the cover has weight `≥ 1` on every such limit pose *with
the boundary points of the lost edges excluded*, plus positive margin elsewhere.  These limit
poses are legitimate constraints of the closed cover (they are limits of admissible poses, and
coverage is upper semicontinuous), so a valid cover satisfies them, but a cover that is tight at a
limit pose cannot be certified by boxes and would need a union-of-regions argument.

Measured on `runs/closed4_best.txt` (1,972 points, 12.4175, the best heuristic LP cover of
`CLOSED4.md`, 91 % of its weight on the grid lines; `zeromargin_fan.py`):

| pose family | captured weight |
|---|---|
| corner square `(w/2 + d, w/2 + d, ε)`, all `ε ≤ 0.03`, `d ≤ 0.01` | exactly **1.0000** — from the points on `{1}×[0.414,1]` and `[0.414,1]×{1}`, all P1-certified |
| axis-parallel grid poses `(cx, cy, 0)`, `cx, cy ∈ {½,1,1½,2}` | 1.00 – 1.97 |
| fan `(1.5, 1.5 + δ, ε)`, `δ ∈ [10⁻⁶, 0.03]`, `ε ∈ [10⁻⁵, 0.03]`, all ratios | **1.039 – 1.11** |
| fan at the wall square `(1.5, ½ + δ, ε)` | **1.006 – 1.06** |
| worst over log-grid fans around all 28 grid poses, both signs | 1.0000 (the corner) |

So the LP cover is *not* fragile at the grid poses: the sampled LP (rows at angles down to
0.001° and wall bands at `10⁻⁷`) has already pushed 4–10 % of slack into the fans, and the corner
is exactly the P1 situation.  Its known failures (`CLOSED4.md`: 0.9927 at `(3.40, 1.42, 76.4°)`)
are ordinary unconverged rows at 14°/76° near a wall, not zero-margin poses.

The exact checker on this cover (`runs/zeromargin_closed4.log`, depth limit 6, 575 s on 4 cores;
`--tri` is useless for weights `≈ 0.17`): 111,020 boxes, 24,905 leaves by weighted CORE, 350 by
P1, 4,459 empty, 28,996 not certified *at that depth* (boxes of side 0.025 × 1.8°, along the
walls — `cx ≈ 0.5`, `cy ≈ 1.4–1.6`, `θ ≤ 7°` first — where the cover's margin is a few per
cent and the weight is spread over many points).  A depth-10 run was stopped after an hour: with
2,000 points the weighted core test costs ~8 ms per box, so this cover needs the float
pre-filter vectorised (numpy over the point set) before it is a practical oracle — an
engineering item, not a semantic one.  "Not certified at depth 6" is not a failure of the cover;
the known failures are the 0.9927 poses above.

**Verdict.**  Rung 2 needs (i) an LP cover at `s = 4` that is a valid cover — the current best
fails by 0.7 % at tilted wall poses, i.e. the `closed4.py` loop must be run with the exact
checker as its separation oracle (its uncertified boxes are the rows to add), which is the same
loop `tighten.py` runs with `verify/`; (ii) no tight tilted family of the item-4 kind in the
converged cover, because the weighted triangle lemma gives only `min_v w_v`.  Whether (ii) holds
is an empirical question about the LP optimum; the item-4 families arise from unit-distance
pairs with nothing else nearby, which a weighted optimum has no reason to produce, but there is
no proof.  If it does happen the tool needed is a two-region certification (box `⊆ R_p ∪ R_q` for
two point regions, exact polygon clipping at the bin endpoints plus the sector arcs), which is
implementable but not built.  Segment resources are *not* needed: their chord length is
discontinuous at the same poses as points (an edge on a grid line), so they buy nothing here;
Nagamochi's 0.9-offsets buy robustness only in the dilated-box regime he works in.  What rung 2
therefore costs is compute (the separation loop at `s = 4`, with 2,000-point covers the weighted
CORE test is ~100× slower than rung 1) rather than new mathematics — unless (ii) fails.

## 6. What transfers to the `s(12)` endgame

The certificate at `t = 4` for `n = 12` will be a clique/branch certificate, not a point cover,
but its verification has exactly this shape: boxes certified by cores, walls and corners by P1
(the region `Q ⊆ p + [−1,1]²` is a pose box, so it is also the natural shape for a box clique),
and tight tilted families by triangle-type lemmas.  The checker's structure (subdivide; certify
by core / lemma; report the rest with locations) is the verification half of the "certificate at
the limit" of `TODO.md` Phase 3, and its uncertified-box output is the separation oracle the
cover side needs at `s = 4`.  Rung 3 (`[0,7]²`, `n = 45`) is the same code with `m = 7` once
task D produces a candidate.

## 7. Task H (2026-08-30): the checker as a fast oracle, and the rung-2 closing loop

**Vectorisation (step 1 of the brief).**  `cert_core` and `cert_p1` (`search/zeromargin.py`) now
pre-filter the *entire* point set for a box with `numpy` (four vectorised corner tests over all
`N` points at once, superset of the exact test by the same `1e-9` lenient tolerance the old
per-point float pre-filter used) before doing any exact `Fraction` work; if the pre-filter's own
weight sum can't reach `1`, the box is rejected with zero exact arithmetic.  Points that pass are
confirmed exactly in decreasing-weight order so the accumulator reaches `≥ 1` in as few `Fraction`
operations as possible.  Everything load-bearing (`in_core`, `in_rot_square`, the exact confirm
loop) is untouched.

*Regression.*  Rung 1 is byte-identical before/after: reduced domain `6,958` boxes /
`3356/70/74/3179/0`, full domain `27,396` / `13556/298/280/12364/0`, `zeromargin_stress.py` `0`
failures on `6,679` leaves and `340k` primitive instances (all re-run and matching ZEROMARGIN.md
§3 exactly).  On the 1,972-point `runs/closed4_best.txt` cover (the rung-2 target), depth 6 now
takes **74 s on 8 processes** (was 575 s on 4) for the *identical* leaf counts (`111,020` boxes,
`24905` CORE / `350` P1 / `0` TRI / `4459` EMPTY / `28996` UNCERTIFIED) — about a 30× per-core
speed-up, matching the brief's "minutes, not hours" target and making a 2,000-point cover a
practical separation oracle.

**New CLI modes.**
- `--oracle FILE`: for every uncertified box, write its centre and the poses at its 4 corners ×
  its 2 angle-bin endpoints (`≤ 27` poses/box after dedup) as `cx cy theta_rad` rows — candidate
  cutting planes for a cover LP (`search/closed4.py` / `search/rung2_close.py` consume these
  indirectly via their own local search seeded at these locations).
- `pose CERT --cx C --cy C --u U`: compute the **exact** (`Fraction`) captured weight of one
  rational pose (`u = tan(θ/2)` rational ⇒ `cos θ, sin θ` rational, `trig(u)`) with *no*
  subdivision or depth limit — a way to turn a float-found candidate violation (from
  `closed4.py stress`/`polish`) into a rigorous yes/no on whether the current point set really
  fails there, independent of the box checker's adaptive-subdivision machinery.

**Rung 2, `m = 4`, `W < 13` (step 2 of the brief).**  `search/rung2_close.py` runs the separation
loop the brief asks for: solve the cover LP over `runs/closed4_best.txt`'s fixed column set (no
new points), stress-test at pitch `0.003` (`search/closed4.py`'s own dense scan) to find the
*actual* worst poses, add those (and a local polish around them) as hard LP rows, resolve, repeat
— each round is blind exactly where the previous round's row set was blind, so it is a genuine
(float) cutting-plane loop against the residual gap CLOSED4.md left open (LP converged at
`12.4174` on its own row lattice, but a finer stress scan found `0.99267` at `(3.40, 1.42,
76.4°)`).  Status and numbers: see `search/FAMILY.md` (filled in once the loop finishes or is
reported stalled, per the brief's instruction to report honestly rather than declare
near-success).  Whether the converged cover develops a tilted tight family needing a two-region
primitive (ZEROMARGIN §4 item 4 / §5(ii)) is recorded there too.
