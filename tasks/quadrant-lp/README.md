# quadrant-lp: can a Nagamochi-type family prove s(k²−4) = k (or s(k²−3) = k) for all k?  (2026-09-28)

**Setting.**  A *fixed-profile family*: Lebesgue interior (density exactly 1), a wall band `[0,w]` along each wall
carrying a 1-periodic profile, four identical corner modules on `[0,R]²`.  A unit square sees only ~√2 around it, so
validity reduces to the **quadrant** `Q = [0,∞)²`: if every closed unit square in `Q` captures mass ≥ 1, the family is
valid in `[0,k]²` for every integer `k ≥ k₀ ≈ 2R + 2`, with saving `k² − μ([0,k]²) = 4D`.
* `σ` (wall deficit per period) must be **exactly 0**: `σ > 0` is impossible (Erdős–Graham / Chung–Graham waste
  bound: savings ≤ C k^0.63), `σ < 0` makes the saving fall linearly in k.  Impose per-period mass = w as an equality.
* Needed: **D > 1** gives s(k²−4) = k for all k ≥ k₀; **D > 3/4** gives s(k²−3) = k for all k ≥ k₀ (Bentz's
  conjecture; proved only k ≤ 7, and s(61) = 8 from our s(60)).  Nagamochi 2005 is this form with D = 1/2
  (`notes/proof-anatomy.md` §3; note his Lemma 1 is false, chelokot counterexample, `notes/literature-s32.md`).
  `search/CORNER_DEFICIT.md`: with Lebesgue walls, D = 0 (proved).
* **Seam bound (derived in-session, unchecked):** corner boxes `[0,n]²`, n integer, have deficit ≤ 0 (dilated grid,
  CORNER_DEFICIT §2.3); the rest of the box splits into half-open wall cells (mass = area by σ = 0) and Lebesgue; the
  only unaccounted mass is on the seams, so **D ≤ m_v** := profile mass on one integer vertical cross-section
  `{i} × [0,w]` of the band.  Nagamochi: m_v = 1/2 = D (tight).  So a unit band (w = 1) gives D ≤ 1.

**Tasks, in order.**
1. **Write the seam bound up carefully** (`search/QUADRANT.md` §1): exact statement, conventions (closed squares,
   line mass at seams, D4 / mirror of the profile), proof, and check it on Nagamochi.  If it is wrong, say how, fix it.
2. **Quadrant LP** (`search/quadrant_lp.py`, reuse `search/line_cover.py` / closed-cover LP machinery and HiGHS):
   columns = corner-module atoms (points + grid-line segments on `[0,R]²`, diagonal symmetry) + periodic-profile
   atoms (points + segments in one period cell of the band, repeated along both walls by the diagonal reflection),
   Lebesgue elsewhere (fixed, area terms).  Rows = unit-square poses in a finite window that sees the corner plus
   one-period wall window (the far wall and interior are periodic / Lebesgue).  Objective: max D, with D defined so
   that the box saving is exactly 4D; write the definition down and check it on hand examples.  σ = 0 as equality.
3. **Sanity:** Lebesgue walls → D ≈ 0 (up to the lattice artifact); Nagamochi's structure fixed → D ≈ 1/2 or a
   measured violation (compare chelokot's counterexample).
4. **Grid:** w ∈ {1, 2, 3}, R ∈ {2, 3, 4}.  Report D_LP, m_v, the seam-bound gap, where the saving sits.
5. **Honesty about artifacts.**  Sampled-row LPs show fake savings (CORNER_DEFICIT §0: LP value halved when q doubled).
   Use row generation with a float separation oracle over continuous poses (x, y, θ) in the window, and report D as
   a function of lattice pitch / rounds until it stabilises.  A number you cannot show is stable is [heuristic].
6. **Only if** a stabilised D clearly exceeds 3/4 (or 1): rescale-free validity is required on the band and interior
   (no ×f trick — σ < 0 kills it), so scope what the exact checker needs (zm_mixed with a periodic window + Lebesgue
   polygon mass; continuum-tight pose families along the wall).  Scope only; don't build it.

**Labels** [proved]/[measured]/[heuristic] as in the other notes.  **Compute:** physical cores 7–11 only (≤ 5
processes; Lean agent is on 12–15), ≤ 40 GB.  **Don't** edit TODO.md, README.md, checkers, certificates, site, lean.
**Deliverables:** `search/QUADRANT.md`, code in `search/`, runs in `runs/quad_*` (gitignored); commit locally on main
in `/home/evand/math/square-packing/public` (no push).  Report within ~6 h even if unfinished: D table, stability
evidence, seam-bound verdict, and your go/no-go for k²−4 and for k²−3.
