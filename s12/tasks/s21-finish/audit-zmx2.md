# Brief: adversarial audit of zmx2 (the Rust mixed-cover checker), 2026-09-27

Read `search/ZMX2.md`, `verify2/src/bin/zmx2.rs`, `search/zmx2_tools.py`, `search/zmx2_tests.sh`, `search/zmx2_manifest.txt`,
`tasks/line-cover/FORMAT.md`.  Cores `taskset -c 10-11`.  Budget: a few hours.  Read-only on zmx2 (propose patches in your note).
You may read `ZM_MIXED.md` for comparison, but judge zmx2 on its own proofs.

Focus: (1) interval arithmetic: every float op outward-rounded? (sqrt, division by intervals containing 0, subtraction
cancellation, conversions from big integers to f64 (exact? > 2^53?), the rounding to the integer grid, casts, overflow in integer
sums (i64/i128)); any path where a float can certify.  (2) Each lemma's proof and its code (germ-pair coupling at θ → 0⁺, chord
bounds, D4 and no-symmetry root domains: do the roots tile `[0,5]²`-admissible poses × θ ∈ [0°, 90°) for the no-symmetry mode?).
(3) Parser: malformed/degenerate inputs, integer bounds.  (4) New adversarial tests: holes at germs, sub-bin boundaries, grid
crossings, tiny u; confirm refusal at the right pose.  Report `search/ZMX2_AUDIT.md`, findings graded must-fix / should-fix / nit.
Commit own files, never push.  Final reply `<= 200` words, most severe first.
