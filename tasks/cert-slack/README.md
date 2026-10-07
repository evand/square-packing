# Certificate slack and the pose oracle (2026-10-02)

Context: k²−4 exact loop (`search/K2M4_MARGIN.md`, `search/qx2_xloop.sh`).  The float LP's oracle
(`quadrant_lp.oracle` + qx2 additions) = initial pitch-0.1 lattice at 19 angles + near-breakpoint poses; exact θ = 0
limit corners and exact germ-limit rows; per round 0.02 lattice at 24 angles, 10⁶ random, near-lattice + germ poses,
random-perturbation hill climb from the 2000 worst seeds, ≤ 4000 rows/round (one per 0.01 × 0.5° bucket).  No gradients,
no exact minimisation.  It missed genuine violations at 1.45° and 13° (κ = 0.05 family) that zmx2 found at once.

Certificate size does not grow with rows (support 155 → 142 → 157 lines, 1976–2072 box segments: elements are fixed);
check cost is driven by TIGHTNESS (the LP pushes every seen constraint to its bound).  σ = 0 tightness is forced and
handled by exact lemmas; the rest is optional.

Levers (biggest first):
1. Spend surplus D on slack: second LP with D ≥ ~1.04, maximise margin (larger κ / θ₀); optionally reweighted-L1
   pruning of support.  First try: κ = 0.08, θ₀ = 3° (`runs/qx2_k4x_k008`), compared with κ = 0.05 on zmx2 counts.
2. Structured oracle: enumerate critical poses (≈ 3-contact configurations: square vertices/edges on lines and
   breakpoints — Lemma E's arrangement), or at least gradient polish instead of the random hill climb.
3. zmx2 as oracle with exact poses: a flag printing uncertified boxes' float-min poses at full precision (today's
   output is rounded to 1e-6, useless at near-tangent tiny-θ poses; we filter to θ ∈ [0.05°, 89.95°]).  Changes the zmx2
   sha → re-pin bundles.
Timing (horizon: first valid k²−4 certificate): 1 now (on the critical path); 2, 3 after, unless the loop stalls.
