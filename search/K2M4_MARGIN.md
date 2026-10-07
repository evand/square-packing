# k² − 4 go/no-go and the tilt-margin ceiling (2026-10-01)

Labels: **[measured]** = float LP + float oracle (heuristic), **[heuristic]** = interpretation.  Code: `qx2_lp.py` (unchanged),
new `qx2_kmax.py` (maximise the margin slope κ, optionally with D ≥ d or, in `--band` mode, m_v ≥ d; writes the binding
rows by dual weight).  Runs: `runs/qx2_k4*`, `runs/qx2_k4d_*`, `runs/qx2_k5_*`, `runs/qxk_*` (gitignored).

## 1. The go/no-go (R = w = 3, lines only, pitch 0.2, strict tolerance, exact axis + germ rows)  [measured]

The k2m3 pipeline (QUADRANT_EXACT §1) at R = w = 3.  The tilt margin asks every row at angle θ for
`1 + κ·min(θ', θ₀)·(1 − area(Q∩U))`.

| κ | θ₀ | D (settled) | note |
|---|---|---|---|
| 0 | — | 1.132 | (lenient-tolerance QUADRANT §4: 1.154) |
| 0.02 | 3° | 1.122 | |
| 0.05 | 0.5° | **1.121** | R = 4: **1.130** |
| 0.05 | 3° | 1.103 | |
| 0.1 | 3° | 1.043 | |
| 0.2 | 3° | infeasible (round 4) | the k2m3 setting |
| 0.2 | 1° | 0.899 at round 5, collapsing | |
| 0.2 | 0.5° | 1.081 at round 4, falling | |

**GO for k² − 4 at κ ≈ 0.05, θ₀ = 0.5°**: float margin ≈ 0.12 over 1 (the exact re-solve cost ~10⁻⁵ at w = 2).  Price: the
checker's first-order slope near the θ = 0 tight families is ¼ of k2m3's (κ = 0.2), so ~4× finer boxes there.

## 2. The slope ceiling κ*  [measured]

`qx2_kmax.py`, κ maximised (D unconstrained unless stated).

| model | w | κ* |
|---|---|---|
| full (R = w) | 2 | 0.2134 |
| full, D ≥ 1 | 3 | 0.1120 |
| band only | 2 | 0.2132 |
| band only | 3 | 0.1129 |
| band only, m_v ≥ 1 | 3 | 0.1123 |
| band only | 4 | 0.0683 |
| band only, m_v ≥ 1.25 | 4 | 0.0641 |
| band only | 5 | 0.0389 |
| band only, m_v ≥ 1.5 | 4 | 0.0473 |
| band only, m_v ≥ 1.5 | 5 | 0.0226 |
| band only | 6 | 0.0073 (pitch 0.2 artifact, see below) |
| band only, **pitch 0.1** | 2 | 0.2457 |
| band only, pitch 0.1 | 3 | 0.1454 |
| band only, pitch 0.1 | 4 | 0.0993 |
| band only, pitch 0.1 | 5 | 0.0734 |
| band only, pitch 0.1 | 6 | 0.0563 |

* **The ceiling is a band phenomenon**: band-only = full model at w = 2 and 3; D ≥ 1 / m_v ≥ 1 barely moves it.
* **Binding rows** (dual weights, w = 3 full D ≥ 1 and band w = 2, 3): all in the far band, axis-near squares
  `[t, t+1] × [j, j+1]`, j = 0..w−1, tilted 0.06°–0.2°, phases spread over the period; **equal dual weight 1/w per row**,
  plus the σ = 0 row (weight 1/3 at w = 2, 1/4 at w = 3).  The certificate is "average a full column of w stacked unit
  squares over one period, plus the band-mass constraint".
* [heuristic] Why: a phase-averaged column of w axis-parallel unit squares captures exactly the band mass per period (= w
  by σ = 0), so **σ = 0 forces every row to be tight at θ = 0** (not only the wall row).  Tilting gains only at the
  column's boundaries (wall, Lebesgue layer), which must pay the margin in all w rows.  A first-order Fubini estimate
  (stacking at spacing cos θ makes adjacent tapers sum to 1/cos θ; the wall leaves a triangle of height sin θ) gives
  κ* ≲ 1/(2(w − δ)) — right order at w = 2, 3 but too high, and it predicts c/w decay, while the data decay faster.
* **Decay** (pitch 0.2): ratios 0.53, 0.61, 0.57, then **0.19** at w = 6 — but that is the lattice: at pitch 0.1 the
  w = 6 ceiling is ≈ 0.057 (8×) and w = 5 ≈ 0.074 (ratio ≈ 0.77).  The binding structure at w = 6 is unchanged (98 %
  near-axis, 1/6 per row, σ weight 1/7 = 1/(w+1) as at every w).  So the element lattice must refine with w; the
  pitch-0.1 series 0.246, 0.145, 0.099, 0.073, 0.056 (ratios 0.59, 0.68, 0.74, 0.77, rising) is a **power law,
  κ* ≈ w^−1.3**, not geometric.  Finer pitch might raise it further (pitch 0.05 not run).

## 3. Consequences  [heuristic]

* k² − 4 (w = 3): κ = 0.05 is half the ceiling at D ≥ 1; first-order lemmas (U/E) suffice, with deeper refinement.
* k² − c in general: with κ* ~ w^−1.3 (pitch 0.1), first-order checker cost near the tight rows grows only polynomially in
  w (~1/κ*), provided the element pitch refines with w.  A **second-order lemma** (certify μ = 1 + 0·θ + cθ² at the tight
  rows, no margin) is a constant-factor improvement, not a necessity, for each fixed c; it matters for uniformity in w.
* Element smoothness (piecewise-linear / quadratic densities): the averaging argument uses only that elements are measures
  with fixed band mass, so it should not raise κ*; its value would be on the checker side (smaller box loss per κ).
  **Measured (`QX2_PL.md`, R = w = 3, pitch 0.2, float rows): confirmed.**  Hats alone lower D at every κ (1.099 / 1.090 /
  0.942 vs uniform 1.132 / 1.121 / 1.042 at κ = 0 / 0.05 / 0.1); κ* unchanged (≤ 0.1126 vs ≤ 0.1135); κ = 0.2 infeasible for
  uniform, hat and uniform+hat.  Uniform+hat gains 0.02–0.03 in D (extra freedom, not kink removal).  The binding rows are
  the σ = 0 wall-row families and their germs in every basis.
* w = 4 (k² − 5), R = 4: κ = 0.05 → D = 1.218, κ = 0.02 → 1.248, κ = 0 → ≈ 1.257 (round 10, then a HiGHS kUnknown; warm
  restart running).  **No-go at R = 4, w = 4.**  R = 5 w = 4 κ = 0.02 settles at 1.2586 (R = 4 κ = 0: 1.256): the corner gains ~0.01 per
  unit R, so **w = 4 caps just above 1.25; k² − 5 needs w ≥ 5**; **R = w = 5 κ = 0.02: D ≈ 1.302 at round 31, still drifting ~10⁻⁴/round** (m_v 1.394, corner
  deficit ≈ −0.09): float margin ≈ 0.05 over 5/4 — a conditional GO for k² − 5, thinner than k2m3's 0.10 over ¾ after
  exactness.  Next levers: R = 6 (corner), pitch 0.1 (κ* and D), κ ≤ 0.02 (κ*(5) ≈ 0.07 at pitch 0.1).
* w = 4 earlier note: band allows κ = 0.064 at m_v = 1.25, but the full R = 4 model at κ = 0.05 sits at D ≈ 1.22, m_v ≈ 1.27
  (round 15: 1.219, settling < 1.25 → no-go at this setting): corner deficit E ≈ −0.05 and an m_v far below the band capacity 1.97.  The corner, not the slope, binds at
  w = 4 (QUADRANT §6 saw the same).  Runs: R = 5 w = 4, and R = 4 w = 4 κ = 0 to separate margin from corner.

## 4. Exact loop toward a k² − 4 certificate (2026-10-02)  [exact unless stated]

Pipeline per candidate (`qx2_xloop.sh RUNDIR 9`): `qx2_exact.py project` (tight θ = 0 + germ rows exact, support fixed,
free vars rounded to 10⁻¹²) → `axis` (Lemma Z, exact) → family + box k = 9 (= 2R + 3, R = 3) → `qx2_germscan.py --hi 9/2`
(exact, 3362 centres × 625 offsets) → zmx2 d4 / cert0 / `cert --d4 --first-order` → rows from uncertified boxes.

| candidate | exact D | zmx2 boxes | zmx2 uncertified | of which θ ≥ 0.05° | float-min < 1 |
|---|---|---|---|---|---|
| κ = 0.05, θ₀ = 0.5° (first projection) | 1.121226946718 | 140.6 M | 10,709 | — | 94 (2 confirmed exactly: 0.9999629 at (3.4875, 1.5122, 1.45°), 0.9999650 at (1.6, 2.6, 13°)) |
| κ = 0.05 + filtered zmx2 rows (`k4x_it2b`) | 1.120958235386 | 141.8 M | 6,562 | 227 | 1 (0.999999) |
| **κ = 0.08, θ₀ = 3°** (`k4x_k008`) | **1.073851127855** | **15.6 M** (772 CPU-s) | 4,844 | **0** | 0 |

All three: projection consistent (0 inconsistent, 0 negative, σ = 0), θ = 0 face min exactly 1, germ scan min exactly 1.
* **Lessons.**  (i) The float LP's oracle (random + lattice + random hill climb) misses genuine violations at generic
  angles; zmx2 finds them in minutes.  (ii) zmx2 prints poses to ~10⁻⁶: rows rebuilt from its log are wrong at
  near-tangent tiny-θ poses, and diagonal images at θ = 90° are outside the quadrant model's domain — together they
  dragged D to 1.080 then 0.896 before being filtered out (keep θ ∈ [0.05°, 89.95°], admissible poses, no images).  The
  genuine zmx2 rows cost only 2·10⁻⁴ in D.  (iii) Spending surplus D on margin (κ 0.05 → 0.08, θ₀ 0.5° → 3°) cuts the
  zmx2 work 9× and removes every generic-angle uncertified box; rank of the exact identities 21–23 → 13.
* What is left is tiny-θ germ boxes, which zmx2's first-order lemmas cannot close at these slopes but qx2_zm's Lemma E
  (exact, zero margin allowed) can: pilot on the κ = 0.05 box, cx ∈ [1.4, 1.6], cy ∈ [3.5, 4.1]: 50/96 roots, 0
  uncertified.  **Certificate run of record**: `qx2_zm.py` on the κ = 0.08 box, 16,200 roots, depth 18, Lemma E for
  u ≤ 1/2 from depth 3, `--dump-leaves` (`runs/qx2_k4x_k008/qxzm_full.*`).

**Result (2026-10-03 02:26): VERIFIED-D4, 0 uncertified.**  `qx2_zm.py` (sha `6294052a…`, zeromargin `640fe453…` =
pinned) on `sol_exact_box9.txt` (input sha `4151d7c4…`): 16,200 roots, 214,336 boxes, max depth 17, **815,343 CPU-s
(≈ 226 CPU-h; 56,808 s wall on 15 procs)** — 4.3× the ~190k CPU-s estimate, nearly all at cx ∈ [2.4, 2.6], cy ∈ [2.5, 3.3] (squares straddling the Lebesgue edge 14/5 and the doubled line 3) and in the u-bin [0, 1/16] (53 % of CPU); the slowest roots' true minimum mass is ≈ 1.0036, so the cost is bound looseness, not a tight pose (review `leaves`, 2026-10-03).  Leaves:
PIECE 84,152, EXACT 20,577 (+ EXACT0 1,424, EXACT45 1,392), SYM 2,676, EMPTY 3,890, LEB 683, CAP 393, AXIS 81, UNCERT 0.
With the θ = 0 face (Lemma Z, exact: `qx2_zm.py axis` min 1, 208 tight in the D4 region; the 348 in §4's xloop output is `qx2_exact.py axis`'s count of tight θ = 0 LP rows, 1,397 corners; Lean `validAxis9` re-proves it over the whole box, 6,400 corners, 832 tight) and the exact germ scan (min 1), `Valid9` holds for this box by
our single implementation: total `81 − 4D = 76.704595…< 77`.  Record: `runs/qx2_k4x_k008/qxzm_full.{out,jsonl}` (not
yet bundled).  Still to do for all k: the all-k reduction (Lean, analogue of `bentz_of_valid7`; task
`lean-k2m4-reduction`), k = 5..7 from s(21), s(32), s(45), bundle, review.
