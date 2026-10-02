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
