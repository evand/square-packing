# k2m3-review: adversarial review of the s(k²−3) = k certificate (2026-09-29)

**Claim under review.**  `search/QUADRANT_EXACT.md` (paths relative to `~/math/square-packing/public/s12`):
the fixed-profile family `search/qx2_data/L4_k02_family.txt` (R = w = 2, σ = 0, Lebesgue on U) is valid, hence
s(k²−3) = k for all k ≥ 6 (Bentz's conjecture).  Certificate of record: run V2, `search/qx2_zm.py` sha
`cdade4b6…a1b7` (commit 0e303ca), records `search/qx2_data/cert/runV2_cdade4b6_full.*`; θ = 0 by Lemma Z
(`qx2_zm.py axis`).  Single implementation, reviewed by nobody so far.

**Your job: try to break it.**  You are a hostile referee.  Default posture: the claim is wrong somewhere; find where.
A finding is worth most if it is concrete (an exact pose with mass < 1, a leaf whose certification is unjustified,
a proof step that fails, a case not covered).  "No problem found" is a fine outcome, but say exactly what you
checked and how, and what you did *not* check.

**Independence.**  For your area, first read only the *statements* (lemma statements, the certificate claim, the
code), and derive / check them yourself before reading the write-up's proofs.  Do not trust the write-up's tests
as evidence; design your own.  The historical bug to learn from: a tile-germ limit (θ → 0⁺ at a centre where
lines lie on the square's edges) 6·10⁻⁵ below 1, invisible to float sampling (§1.5).  Tiny tilts, one-sided
limits, degenerate/zero-width boxes, boundary lines of U, D4 bookkeeping and "dropped because float says empty"
are where to look.

**Rules.**  Read-only on everything in `public/` (no edits to any existing file, no commits, no git operations that
change state).  Write your scratch code and report in `~/math/square-packing/private/s12/tasks/k2m3-review/<area>/`.
Exact arithmetic (`fractions.Fraction`) for any claim of a violation.  **Compute:** only your assigned physical
cores (`taskset -c`), and keep it light (this is a reasoning task; ≤ 1 CPU-h total unless a concrete lead needs
more).  ≤ 8 GB RAM.  Python only (the `_untrusted-third-party/` rule: never execute anything there).

**Report** (`<area>/REPORT.md`, and summarise in your final message, ≤ 400 words): findings ranked by severity
(**BREAKS** the claim / **GAP** in a proof or code that could hide a failure / **MINOR** wording, dead code,
documentation), each with the evidence and file:line; then "checked, OK" items with how; then "not checked".
