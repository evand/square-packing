# Quadrant LP: fixed-profile families for s(k²−4) and s(k²−3) (task quadrant-lp, 2026-09-28)

Brief: `tasks/quadrant-lp/README.md`.  Code: `search/quadrant_lp.py` (model, row-generation LP, float oracle),
`search/quadrant_tools.py` (pose breakdown, deep local search, independent box accounting).  Runs:
`runs/quad_*/` (gitignored; `log.txt`, `hist.json`, `sol.npz`, `support.txt`).

Labels: **[proved]** (exact argument given here), **[measured]** (float LP / float oracle numbers, reported as-is),
**[heuristic]** (interpretation).

## 0. Answer

> * **Seam bound: correct as stated** [proved, §1.4]: for a valid fixed-profile family, `D = E(R) + m_v` with
>   `E(R) = R² − μ_Q([0,R]²) ≤ 0` (dilated grid), so **`D ≤ m_v ≤ w`**.  Tight for Nagamochi (`E = 0`, `D = m_v = ½`)
>   and, to within 1 %, for every LP optimum below.  Conventions that make it exact: closed corner box, band closed
>   at `x = R`, `m_v` = profile mass on one integer cross-section, profile mirror-symmetric (§1.1–1.2).
> * **Sanity** [measured]: Lebesgue walls give `D = 0` once the oracle has near-lattice (`ε`-gap) poses — without them
>   random separation "certified" `D = 0.0059`, a pure artifact (§3.1).  **Nagamochi's structure is invalid** in the
>   quadrant: a wall-resting square `≈ [1,2] × [0,1]` tilted −0.62° captures **0.955** (→ 0.95 as the tilt → 0),
>   because the corner Q-point weighs 0.45 (§3.2).
> * **`D` table** (§4): `w = 1`: **0.500** (`R = 2`); `w = 2`: **0.945** (`R = 2`, stable), 0.951 (`R = 3`, falling);
>   `w = 3`: **1.154** (`R = 3`, pitch 0.2, stable), 1.21 (pitch 0.1, falling), 1.17 (`R = 4`, pitch 0.2, falling).
>   Band-only bound `m_v ≤ 0.750 / 1.433 / ≤ 1.755` for `w = 1, 2, 3`.
> * **Where the saving sits** (§6): on the seams.  `E(R) ≈ 0` in every optimum; the band's integer cross-section
>   (a vertical line segment `{i} × [≈0.5, w]`) carries `m_v`, the corner module is essentially the lines
>   `x = 1`, `y = 1` continuing that profile.  `w = 2` reads as two stacked Nagamochi rows sharing the seam line.
> * **Go/no-go** (§7): `w = 1` no-go for both (proved `D ≤ 1`; measured `D ≤ m_v^band = ¾`, `D ≈ ½`).
>   **`k² − 3`: GO** to exact verification with `R = 2, w = 2` (`D = 0.945 ± 0.002`, margin 0.195 over ¾, residual
>   float violations ≤ 0.14 %).  **`k² − 4`: conditional GO** with `R = 3, w = 3` (`D ≈ 1.154` at pitch 0.2, margin
>   0.15 over 1, same residual level, younger runs).  Nothing is proved: all `D` values are float LP + float oracle.
> * **Exact checker scope** (§8): one box `k = 2R + 3` suffices; needs polygon masses in cert mode, an exact rational
>   re-solve with `σ = 0` exactly, and zero-margin handling of continuum-tight wall families (no ×f allowed).


## 1. The seam bound

### 1.1 Setting and conventions

* **Closed semantics.**  A measure `μ ≥ 0` on a closed region `X` is *valid* if `μ(S) ≥ 1` for every closed unit
  square `S ⊂ X` (any angle), mass on `∂S` counting.  For `X = [0,k]²` (k integer) the *saving* is `k² − μ([0,k]²)`.
* **Family** (parameters: integer `R ≥ 1`, band width `0 < w ≤ R`):
  * a *profile* `π`: a locally finite measure on the strip `ℝ × [0,w]`, invariant under `x ↦ x+1` and under the
    mirror `x ↦ −x`;
  * a *corner module* `ν`: a finite measure on `C = [0,R]²`, invariant under the diagonal reflection `(x,y) ↦ (y,x)`;
  * the *quadrant measure* on `Q = [0,∞)²`:
    `μ_Q = ν + π|_{[R,∞)×[0,w]} + π'|_{[0,w]×[R,∞)} + Leb|_L`, `L = {x > w, y > w} \ [0,R]²`,
    where `π'` is the diagonal reflection of `π`.  The band pieces are **closed at `x = R`** (resp. `y = R`); any other
    split between `ν` and the band is a relabelling.
* **Per-period deficit** `σ := w − π([0,1) × [0,w])`; **seam mass** `m_v := π({0} × [0,w])` (the mass of one
  integer cross-section of the band: points on it, and vertical segments along it; horizontal segments and area
  contribute 0).
* **Mirror symmetry of `π`** is what makes the D4 box consistent: the bottom band seen from the right-hand corner is
  `x ↦ k − x` of the band seen from the left, and for integer `k` this is the same measure iff `π` is mirror
  symmetric (on the phase-0 lattice).  Without it one needs two corner modules (a C4 box); not considered here.
* **Box measure** `μ_k` (integer `k ≥ 2R + 2`): `ν` and its images at the four corners, `π` on `[R, k−R] × [0,w]`
  along each wall (rotated), Lebesgue on `[w, k−w]²` minus the four corner squares.

### 1.2 Box accounting and the definition of `D`  [proved]

> `k² − μ_k([0,k]²) = 4D + 4σ(k − 2R)`,  with  **`D := R² − ν(C) − m_v`**.

Proof.  `μ_k([0,k]²) = 4ν(C) + 4π([R, k−R] × [0,w]) + area(Lebesgue part)`.  The interval `[R, k−R]` is `k − 2R`
half-open periods plus the closed end `{k−R}`, so `π([R,k−R]×[0,w]) = (k−2R)(w − σ) + m_v`.  The Lebesgue part has
area `(k−2w)² − 4(R−w)² = k² − 4R² − 4w(k − 2R)`.  Summing: `k² − 4(R² − ν(C) − m_v) − 4σ(k−2R)`.  ∎

So with `σ = 0` the saving is `4D` for every `k ≥ 2R+2`; `σ > 0` would make it grow linearly (impossible, waste bound),
`σ < 0` makes it fall linearly.  `σ = 0` is imposed as an equality.

Equivalently, with `E(n) := n² − μ_Q([0,n]²)` (the closed corner-box deficit in the quadrant):
`μ_Q([0,R]²) = ν(C) + 2 m_v` (the two band cross-sections at `x = R` and `y = R`; `[0,R]² ∩ L = ∅`), so
**`D = E(R) + m_v`**.  With `σ = 0`, `E(n) = E(R)` for every integer `n ≥ R` (the L-shell `[0,n+1]² \ [0,n]²` has mass
exactly `2n+1`: two half-open band cells of mass `w` each plus Lebesgue area `2n + 1 − 2w`).

Hand checks (`quadrant_tools.box_saving`, which builds `μ_k` element by element, independent of the formula):
* Lebesgue corner + Lebesgue band (`R = 3`, `w = 2`): saving 0, `4D = 0`, for `k = 8, 11`.
* Nagamochi's structure (§3.2): saving 2 = `k² − (k² − 2)` = his total score `ab − (a+1−⌈a⌉) − (b+1−⌈b⌉)` at
  `a = b = k`; `4D = 2`; `k = 6, 7, 10`.

### 1.3 Localisation: quadrant validity implies box validity  [proved]

If `μ_Q` is valid in `Q` then `μ_k` is valid in `[0,k]²` for every integer `k > 2R + √2` (so `k₀ = 2R + 2`).
The bottom-left chart (identity) agrees with `μ_k` on `[0, k−R)²`: there `μ_k` has the BL corner module, the bottom
and left bands (the bottom band agrees with `π|_{[R,∞)}` by periodicity, and near `x = k−R` by the mirror symmetry and
`k ∈ ℤ`) and Lebesgue; the right and top bands start at `x = k − w ≥ k − R`.  A closed unit square `S ⊂ [0,k]²` whose
centre lies in `[0, k/2]²` lies in `[0, k/2 + √2/2]² ⊂ [0, k−R)²`, so `μ_k(S) = μ_Q(S) ≥ 1`; other centres use the
other three charts (D4).  ∎

### 1.4 The seam bound  [proved; the brief's statement is correct]

> **If `μ_Q` is valid in `Q`, then `E(n) ≤ 0` for every integer `n ≥ 1`; hence `D ≤ m_v`.  Since the band cell has
> mass `w` and contains the cross-section, `m_v ≤ w`; so `w = 1` gives `D ≤ 1`, and `D > 1` needs `w > 1`.**

Proof (the dilated grid of `CORNER_DEFICIT.md` §1/§2.3, in `Q`).  For `ε > 0` the `n²` closed squares
`S_ij = [i(1+ε), i(1+ε)+1] × [j(1+ε), j(1+ε)+1]`, `0 ≤ i, j < n`, lie in `Q`, are pairwise disjoint, and their union
contains `[0,n]²` up to strips of width `ε`; all of them lie in `[0, n + (n−1)ε]²`.  Validity gives
`μ_Q([0, n+(n−1)ε]²) ≥ Σ μ_Q(S_ij) ≥ n²`; letting `ε ↓ 0` (continuity from above, `μ_Q` locally finite) gives
`μ_Q([0,n]²) ≥ n²`.  With `n = R`: `ν(C) + 2 m_v ≥ R²`, i.e. `D = R² − ν(C) − m_v ≤ m_v`.  ∎

Remarks.
* The closed box `[0,R]²` contains both seams `{R}×[0,w]` and `[0,w]×{R}`; in the box accounting each wall's open
  middle strip `(R, k−R) × [0,w]` misses one cross-section (mass `m_v`), which is exactly the slack `D − E(R)`.  The
  bound is the statement that the only saving is the "missing seam" of each wall: `k` closed unit squares cannot sit
  disjointly in a row of length `k`.
* It uses axis-parallel squares only, and only the corner box; nothing about the profile beyond `σ = 0`.
* Tightness: `D = m_v` iff the corner box has zero closed deficit, `E(R) = 0`.  Nagamochi: `E(2) = 4 − 3 − 2·½ = 0`,
  `m_v = ½ = D` (tight, as the brief says).  But Nagamochi's structure is *not valid* in closed semantics (§3.2).
* The bound says nothing about how large `m_v` can be; that is what the LP measures (§2–4).

## 2. The quadrant LP (`search/quadrant_lp.py`)

**Columns** (all nonnegative; one variable per symmetry orbit, each element of the orbit carrying the full mass):
* corner module on `[0,R]²`, diagonal orbits: points on the `hc`-lattice (wall points omitted: in closed semantics a
  point on `∂Q` is never needed, CLOSED4.md), axis-parallel segment pieces `[t, t+hc]` on the lattice lines
  (uniform by length), uniform cells of side `hca`;
* profile, phase cell `[0,1) × [0,w]`, mirror orbits `φ ↔ −φ`: points on the `hp`-lattice, horizontal and vertical
  pieces, cells of side `hpa`; each element is copied to `x = j + φ`, `j = R … R+4` along the bottom wall and
  reflected to the left wall.
* Lebesgue on `L` enters each row as a constant (exact polygon clipping).

**Objective** `max D = R² − ν(C) − m_v` (§1.2; `m_v` = phase-0 points and vertical pieces).  **Equality**
`π(cell) = w` (`σ = 0`).  **Rows**: closed unit squares `(cx, cy, θ)`, `θ ∈ [0, 90°)`, `cx ≥ cy` (diagonal symmetry),
`cy ≤ R + 0.75`, `cx ≤ R + 1.75`.  This window contains every pose that touches the corner module plus one full
period of the band beyond the corner's reach (`R + √2/2 + 1`); every other pose is a translate of one of these
or lies in Lebesgue.  Poses are kept `10⁻⁷` inside `Q`.

**Row generation / oracle** (float, `TOL = 10⁻⁹` closed containment in the safe direction).  Initial rows: a
pitch-0.1 lattice at 19 angles + the near-lattice family.  Each round the oracle evaluates the current solution on
* a pitch-0.02 lattice at 24 angles, and a shifted copy;
* 10⁶ random poses (30 % resting on the bottom wall, 8 % in the corner, 25 % at special angles);
* **near-lattice poses**: axis-parallel and `10⁻⁶`-tilted squares whose edges sit `±10⁻⁷, ±10⁻⁵` off the atom
  lattice (the `ε`-gap poses of CORNER_DEFICIT §0–1), plus 2·10⁵ random poses within `2·10⁻³` of them;
* **germ poses** (LINE_COVER §2): centres `lattice + ½ + θ(t_x, t_y)`, `t` on the atom grid of `[−½, ½]`,
  `θ = ±10⁻⁵, ±10⁻³` (tiles blown up at infinitesimal tilt, where line densities are discontinuous);
* a random-descent polish of the worst 2000 (diverse) poses, 8 rounds × 24 trials, with wall/corner projection;

and adds up to 4000 violated poses (`< 1 − 10⁻⁷`, one per `0.01 × 0.01 × 0.5°` bucket).  Stop: oracle min
`≥ 1 − 10⁻⁶` in three independent passes (two with doubled random sets and pitch 0.0154), or the round limit.
LP: HiGHS IPM (no crossover), one persistent model, rows appended.  One process per run, pinned to one physical core.

**Independent checks** (`search/quadrant_tools.py`):
* `check`: a heavier, differently seeded oracle (≈7.8·10⁶ poses: 4·10⁶ random, pitch-0.013 lattice at 43 angles,
  near-lattice at `d ∈ {10⁻⁸, 10⁻⁶, 10⁻⁴, 10⁻³}`, germs at 6 angles and at half pitch; polish of 6000 seeds × 12 rounds).
* `box`: builds the D4 box measure `μ_k` element by element (no quadrant charts, no periodicity), checks
  `k² − μ_k([0,k]²) = 4D` by direct summation, and runs a whole-box oracle (2·10⁶ random + near-lattice + polish).
* Point and segment masses were cross-checked against `line_cover.Cover` (independent chord code) on 2·10⁴ random
  poses of the `R = 2, w = 2` box: max difference `3·10⁻¹⁵`; cell/Lebesgue areas against Monte Carlo (within noise).

What the LP value means.  With finitely many rows, `D_LP` is an *upper* bound on the best `D` over the given atom
lattice (a relaxation), and it decreases as rows are added; it is a *lower* bound on nothing.  A value is only
trusted when (a) the oracle, and then the heavy check, find no poses below `1 − O(10⁻³)`, and (b) it has stopped
moving under further rounds.  Nested atom lattices (`0.2 ⊂ 0.1 ⊂ 0.05`) have nested feasible sets, so the *true*
values are monotone in the pitch; an LP value that *rose* when the atoms were coarsened would be an artifact.

## 3. Sanity checks

### 3.1 Lebesgue walls (band = Lebesgue, corner free, `R = 2`): `D = 0`  [measured; proved in CORNER_DEFICIT]

| rows | `D_LP` | oracle min |
|---|---|---|
| sampled lattice only (pitch 0.05, 19 angles) | **0.1926** | 0.807 |
| + 2 rounds of random/polish rows | 0.0273 → **0.0059** | 1.0000 (oracle blind!) |
| + near-lattice rows (from round 0) | **0.000000** | — |

The 0.1926 is the lattice artifact of CORNER_DEFICIT §0 (their `q = 20`: 0.198).  The instructive line is the
second: two rounds of random + polish separation brought `D_LP` to 0.0059 and then reported **no violation at
all**, although `D = 0` is a theorem.  The missed poses are the dilated-grid squares `[1+ε, 2+ε]²` (mass 0.9942
at `ε = 10⁻⁵…10⁻⁷`): an `ε`-gap that random poses essentially never hit.  This is why the near-lattice and germ
families are in the oracle; with them the LP returns `D = 0` at round 0.

### 3.2 Nagamochi's structure fixed (`R = 2`, `w = 1`): `D = ½`, **invalid**  [measured]

Corner: Lebesgue `[1,2]²`, segments `y = 1`, `x = 1` on `[0.9, 2]` at 0.5/length, points `(0.9,1)`, `(1,0.9)` at
0.45; profile: line `y = 1` at 0.5/length, point `(i, 0.9)` at 0.5.  `σ = 0`, `m_v = ½`, `E(2) = 0`, `D = ½`
(box saving 2 at `k = 6, 7, 10`, matching his total score `k² − 2`).  The oracle finds closed unit squares with mass
**0.9554** (deep search; pose `(1.49570, 0.50538, −0.620°)`): the square `≈ [1,2] × [0,1]`, resting on the wall,
tilted by a fraction of a degree so that its right edge passes left of the P-point `(2, 0.9)`.  It keeps the Q-point
`(1, 0.9)` (0.45) and the chord of `y = 1` (≈ 0.497) and nothing else; as the tilt → 0 the mass → `0.45 + 0.5 =
0.95`.  At mid-wall the same pose captures `0.5 + 0.5 = 1` (tight); at the corner the Q-point weight 0.45 < 0.5 is
what breaks it.  This is the corner analogue of chelokot's counterexample to Lemma 1 (side 1.0001, score 0.9775);
the whole-box oracle at `k = 6` finds the same failure near every corner (0.970 with a light polish), and the LP
  oracle with its germ family hits it directly: 0.9507 at `(1.4996, 0.5005, −0.057°)`.

### 3.3 Nagamochi's profile fixed, corner free (`R = 2`, pitch 0.1)  [measured]

`D_LP` 0.500 (round 0) → 0.4522 → **0.451826** (rounds 4–9, oracle min 0.9993).  With his profile the best
pitch-0.1 corner is ≈ 0.048 short of ½; a corner at `D = ½` for that profile would need atoms off this lattice
or does not exist.  The free-profile runs below reach ½ with a different (smeared) profile.

## 4. Results: the `D` table  [measured]

`hc = hp` = point/segment pitch (corner / profile), cells `0.25` (pitch 0.1) or `0.5` (pitch 0.2).  "rounds" = row-generation
rounds with the full oracle (near-lattice + germs), after any warm start.  "oracle" = min over the last round's oracle;
"heavy" = `quadrant_tools.py check` on the final solution; "box" = whole-box oracle on the element-built `μ_k`.
`E(R) = D − m_v` is the closed corner-box deficit (the seam bound says `≤ 0`).

| `w` | `R` | pitch | rounds | first → last `D_LP` | `m_v` | `E(R)` | oracle (last) | heavy | box | status |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 2 | 0.1 | 7+10 | 0.830 → **0.50042** | 0.50042 | 0.0000 | 0.9998 | — | — | converged (flat to 3·10⁻⁵ over 4 rounds) |
| 1 | 3 | 0.1 | 4+3 | 0.830 → 0.5021 | 0.5021 | 0.0000 | 0.951 | — | — | stopped, falling |
| 1 | 4 | — | — | — | — | — | — | — | — | not run (w = 1 is closed by the band bound) |
| 2 | 2 | 0.2 | 12 | 1.003 → **0.93675** | 0.93675 | 0.0000 | 0.9999 | — | — | converged |
| 2 | 2 | 0.1 | 9+18 | 1.059 → **0.945057** | 0.94506 | 0.0000 | 0.9997 | **0.99856** | **0.99955** (`k=6`), **0.99970** (`k=7`) | converged (flat to 4·10⁻⁵ over 8 rounds) |
| 2 | 2 | 0.05 | 1 | 0.9555 (warm: the 0.1 rows) | 0.9555 | 0.0000 | 0.972 | — | — | one round only (LP 25 min); not informative |
| 2 | 3 | 0.1 | 11 | 1.096 → 0.9507 | 0.9507 | 0.0000 | 0.991 | — | — | stopped, falling ~1·10⁻⁴/round |
| 2 | 4 | — | — | — | — | — | — | — | — | not run |
| 3 | 3 | 0.2 | 20 | 1.333 → **1.15429** | 1.1616 | −0.0073 | 0.9998 | **0.99853** (rd 13) | **0.99974** (`k=8`), **0.99858** (`k=9`) | converged (flat to 2·10⁻⁴ over 7 rounds) |
| 3 | 3 | 0.1 | 10 | 1.401 → 1.2092 | 1.2093 | −0.0001 | 0.9975 | — | — | stopped, falling (−0.0004 last round) |
| 3 | 4 | 0.2 | 8 | 1.356 → 1.1693 | 1.1788 | −0.0095 | 0.992 | — | — | stopped, falling (−0.003 last round) |
| 3 | 2 | — | — | — | — | — | — | — | — | excluded: needs `R ≥ w` (the two bands would overlap on `[R,w]²`) |
| ½-prof. | 2 | 0.1 | 10 | 0.500 → 0.45183 | 0.5 (fixed) | −0.048 | 0.9993 | — | — | Nagamochi's profile fixed (§3.3) |

`R3 w3` pitch 0.1 above pitch 0.2 and `R4 w3` above `R3 w3` are the expected directions (nested lattices, nested
corners); both were still falling when stopped, so only the pitch-0.2 `R = 3` value is a stable number.

Band-only LP (no corner: the far band alone, `max m_v` subject to `σ = 0` and validity of every square in the
half-plane `y ≥ 0` with Lebesgue above `y = w`).  By the seam bound this is an upper bound for `D` at that `w` and
that atom lattice: **`D ≤ m_v ≤ m_v^band(w)`**.

| `w` | `m_v^band` | rounds | oracle (last) | status |
|---|---|---|---|---|
| 1 | **0.750033** | 6 | 1.000000 in 3 passes | converged (stop criterion met) |
| 2 | **1.43308** | 12+23 | 0.99998–1.0 | converged (flat to 10⁻⁵ over 20 rounds) |
| 3 | ≤ 1.7554 | 15 | 0.9997 | stopped, falling slowly |

## 5. Stability: what is stable and what is an artifact

* **Artifacts seen and removed.**  (i) Lattice rows only: every run starts 0.3–0.5 too high (`R2 w1`: 0.830,
  `R2 w2`: 1.059, `R3 w3`: 1.401) with oracle minima 0.3–0.7; these are the `ε`-gap / germ artifacts, gone after 1–3
  rounds.  (ii) Random + polish separation alone is blind to `ε`-gaps (§3.1: it certified 0.0059 for a problem whose
  value is 0).  With near-lattice + germ families in the oracle the Lebesgue-wall LP returns exactly 0 at round 0.
  (iii) Germ poses: adding the germ family to the `R2 w2` run (restart with its 29,709 rows, round 8 → round 0 of
  `quad_R2_w2_g`) exposed 209 violated germ poses at once (min 0.977), but they cost only `4·10⁻⁵` in `D`.
* **`R = 2, w = 2` is the stable number.**  `D_LP`: 1.059 → 0.971 → 0.960 → 0.952 → 0.948 → 0.9468 → 0.9462 →
  0.9457 → 0.94554 (first run, 9 rounds) → 0.94551 → … → **0.945057** (restart with germs, 18 rounds: the last 8
  rounds move it by `4·10⁻⁵`).  Oracle min in the last rounds 0.9991–0.9997; heavy check min **0.99856** (39 of
  6000 polished poses below `1 − 10⁻⁶`, all near `(1.6, 1.6)` at 75–80°, inside the corner module); box oracle
  **0.99955** (`k = 6`) and **0.99970** (`k = 7`); box saving by direct summation `3.780229 = 4 D` at `k = 6` and `7`.
  The residual violations are ≤ 0.14 %; the last 19 rounds (violations 2.3 % → 0.1 %) cost `6·10⁻⁴` in `D`.
  **Estimate: `D(R=2, w=2, pitch 0.1) = 0.945 ± 0.002`** [measured].
* **Pitch.**  `R2 w2`: pitch 0.2 → **0.93675** (12 rounds, oracle 0.99988, flat for 6 rounds); pitch 0.1 → 0.94506.
  The coarse lattice is a subset of the fine one, and the value went *down* when coarsened, by 0.008: the direction a
  real value must move, not the collapse an artifact shows (CORNER_DEFICIT: halving with `q`).  Pitch 0.05 did not get past
  one round (25-min LP on 58k rows × 3752 columns); its round-0 value 0.9555 is a warm-started relaxation, not a measurement.
  `R3 w3`: pitch 0.2 → 1.1543 (converged); pitch 0.1 → 1.2092 (falling, not converged): same direction.
* **`w = 1`** (`R = 2`): 0.830 → 0.5096 → 0.5014 → … → **0.50042** (7 + 10 rounds, oracle 0.9998).  `R = 3`: 0.5021
  after 7 rounds, still falling (stopped).  Both are consistent with `D(w=1) = ½` and nothing above it.
* **`R = 3`** runs are 5–10× slower per round (LP 150 s + oracle 250 s) and were stopped before full convergence;
  their last values are upper-ish estimates (they only fall with more rounds), see the table.
* Nothing stabilised was checked in exact arithmetic.  Every number here is a float LP value with a float oracle.

## 6. Where the saving sits

* **Entirely on the seams**: in every run `E(R) = D − m_v ∈ [−0.01, 0]`, i.e. the closed corner box carries its
  full area (the seam bound is tight to within 1 %), and the box saving is `4 m_v` = the four "missing"
  cross-sections of §1.4.  Nothing is saved inside the corner module, as CORNER_DEFICIT requires.
* **`w = 1` optimum** (`R2 w1`, `D = 0.5004`): profile = a vertical segment `{i} × [0.5, 1]` of mass 0.4995 (a
  smeared Nagamochi P-point, density rising to 1.4–1.7 on `[0.8, 1]`), the top line `y = 1` on phases
  `[0.1, 0.5] ∪ [0.5, 0.9]` (0.32 per period), a density-1 cell `[0, ¼] × [¾, 1]` either side of the seam (0.125).
  Corner: the lines `x = 1`, `y = 1` on `[0.5, 2]` (0.83 each, the same density profile as the band's seam line near
  the wall) plus density-1 cells.  This is Nagamochi's architecture with the P-point smeared into a segment; it reaches
  his `½` where his own weights (§3.2) fail.
* **`w = 2` optimum** (`R2 w2`, `D = 0.945`): the seam line `{i} × [0.4, 2]` carries `m_v = 0.945` (density 0.6–0.9 on
  `[0.5, 0.9]`, a spike 2.3 on `[0.9, 1]`, 0.35–0.67 on `[1.1, 2]`); the rest of the period (1.055) is horizontal
  pieces at `y ≈ 0.7–1.0` and `y ≈ 1.5–2.0` on phases `0.1–0.5` (mirrored) and cells in `[0, ½] × [¾, 1]` and
  `[0, ½] × [1¾, 2]`.  Read as two stacked Nagamochi-type rows sharing one vertical seam line: each row contributes
  ≈ 0.47.  Corner: again just the lines `x = 1`, `y = 1` (0.99 each), continuing the seam-line profile, plus traces.
* **Band bound vs corner**: the band alone allows `m_v = ¾` (`w = 1`), `1.433` (`w = 2`), `≤ 1.755` (`w = 3`,
  not converged); the corner costs `≈ ¼` at `w = 1`, `≈ 0.49` at `w = 2`.  So the corner, not the band, is the
  binding constraint; this is why larger `R` might help (it does, slightly, at `w = 2`: see the table).
* Tight poses in the `R2 w2` solution (rows with slack `< 10⁻⁶`, IPM interior point, so a lower count): 60 %
  axis-parallel, 62 % resting on the wall, 53 % in the far band (e.g. every lattice phase
  of `[t, t+1] × [0,1]`: a continuum-tight family along the wall), the rest tilted 0–45° near the corner.

## 7. Go / no-go

* **`w = 1`: no-go for both** [proved + measured].  `D ≤ m_v ≤ w = 1` (seam bound, proved) already rules out
  `D > 1`.  Measured: the band alone allows `m_v ≤ 0.7500` (stable, 3 clean oracle passes), so `D ≤ ¾` and `D > ¾`
  is out too (up to the band LP's float accuracy, and on the pitch-0.1 atom lattice); the corner-coupled LP gives
  `D = 0.5004` (`R = 2`) and `≈ 0.502` (`R = 3`, falling).  Nagamochi's `½` is the ceiling of this width.
* **`k² − 3` (`D > ¾`): GO, for the exact-verification step** [measured, not proved].  `R = 2, w = 2` gives
  `D = 0.945 ± 0.002` at pitch 0.1 (0.937 at pitch 0.2): stable over 27 rounds, heavy-check violations `≤ 0.14 %`,
  whole-box check at `k = 6, 7` consistent, direct box sum `= 4D`.  The margin over ¾ is 0.195, a hundred times the
  residual violation depth, and it does not shrink under refinement.  If it survives an exact check it gives
  `s(k² − 3) = k` for every `k ≥ 6` (`k₀ = 2R + 2`), i.e. Bentz's conjecture (proved only for `k ≤ 7`).
* **`k² − 4` (`D > 1`): conditional GO at `w = 3`** [measured, less mature].  `w = 2` stays below 1 (`R = 2`: 0.945;
  `R = 3`: 0.951 and still falling slowly).  `w = 3, R = 3`: pitch 0.2 converges to **1.1543** (20 rounds; heavy check 0.99853,
  whole-box 0.99974 at `k = 8` and 0.99858 at `k = 9`, box sum `= 4D = 4.6175`); pitch 0.1 was still falling (1.21 after 7 rounds) when stopped.  Since the pitch-0.2
  lattice is contained in the pitch-0.1 one, the pitch-0.2 value is the one to bet on; it clears 1 by ≈ 0.15.  If it
  survives an exact check it gives `s(k² − 4) = k` for every `k ≥ 8`, which with the repo's `k = 5–8` covers
  (`s(21)`, `s(32)`, `s(45)`, `s(60)`; S60_COVER.md) would leave only `k = 4` (`s(12)`).
* **Why this is not absurd** [heuristic].  (i) The repo's own exact-checked box covers already save `> 4` at
  `k = 5, 7, 8` (S60_COVER §4: 4.21, 4.54, 4.60 at their exact thresholds, LP values up to 5.60), with the saving
  sitting in the wall ring, which is exactly the structure found here; a fixed-profile family is a restriction of
  those LPs, and its values (3.78 at `w = 2`, ≈ 4.6 at `w = 3`) sit below theirs.  (ii) The savings of a fixed family
  must stay `O(k^0.63)` only as `w` grows with `k`; for fixed `w` any constant is allowed.  (iii) If I recall
  Roth–Vaughan (1978) correctly (waste `W(x) ≥ c·√(x(x − ⌊x⌋))`, `c ≈ 10⁻¹⁰`), `s(k² − c′) = k` is already known
  for astronomically large `k`, so constant savings > 4 are the expected truth, and the open question is only
  small `k` — which a fixed family with `k₀ = 8` would close.  (Not re-checked; do check before quoting.)
* **What would make me say no-go instead.**  An exact checker finding a violation that the LP cannot repair
  cheaply, or `D` collapsing under a finer oracle (it did not: germs, near-lattice at four offsets, 7.8·10⁶ heavy
  poses and whole-box checks moved `R2 w2` by `< 10⁻³`).

## 8. Scope of the exact checker (task 6; scope only, nothing built)

1. **One box suffices.**  Quadrant validity ⇒ box validity for all `k ≥ 2R + 2` (§1.3).  Conversely every quadrant
   pose is congruent to a pose in the faithful region `[0, k−R)²` of the box `k = 2R + 3` (poses touching the corner
   module have `x, y ≤ R + 1.75 + √2/2 < k − R`; far-band poses reduce mod 1 to `cx ∈ [R + 0.71, R + 1.71]`).  So an
   exact check of the single box `k = 2R + 3` (`k = 7` for `R2 w2`, `k = 9` for `R3 w3`) with the existing mixed-cover
   machinery certifies the family for every `k ≥ 2R + 2`.  No periodic-window checker is needed.
2. **Polygons.**  The Lebesgue part (`[w, k−w]²` minus the corner squares: for `R = w` a single square) and the
   uniform cells need polygon masses.  `zm_mixed.py` has polygon code, but it is disabled in `--cert-mode`; it would
   have to be enabled and audited (or the cells replaced by dense grid lines, which changes the LP).
3. **Zero margin, no ×f.**  The shipped covers are certified at a scale `f* > 1`; here scaling the band is fatal
   (`σ < 0`, the saving falls linearly in `k`).  The float solution is tight on continuum families (§6: squares
   `[t, t+1] × [0, 1]` along the wall capture exactly 1 for every `t`; germs; squares barely overlapping the
   Lebesgue region capture `area` in the limit).  Needed: (a) an **exact rational re-solve** of the final LP on its
   active rows, with `σ = 0` exactly (masses as fractions, the equality row kept); (b) a checker that closes
   **zero-margin continuum-tight families**: Lemma T already does germs; the wall-sliding family and the
   Lebesgue-boundary slivers need their own exact lemmas (the second is trivial if the profile and module are
   exactly Lebesgue in a layer along their boundary with `L`, a constraint worth adding to the LP).  (c) Margin can
   be *bought* only in the corner module (finite mass, cost `O(δ R²)` in `D`); the band has to be exact.
4. **Size.**  `R2 w2` pitch 0.1: 1002 variables, of which ≈ 170 in the support; `R3 w3` pitch 0.2: 507 variables, ≈ 195 in the support.
   The box `k = 7` then has a few thousand points/segments — small next to the `s(60)` cover (23,744 points).
5. **Before any of that**: re-run the LP with the "Lebesgue layer" constraint, pitch 0.05 near the seam, and a
   crossover (vertex) solution, and re-measure; the numbers above are from IPM interior points.
