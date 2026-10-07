# Lean: reduce Valid7 / Valid9 to what the Python run checks; then Lemma Z (θ = 0) in Lean (2026-10-03)

Context (paths relative to `public/s12`): `lean/Sqpack/Bentz.lean` (`Valid7`, `bentz_of_valid7`), `lean/Sqpack/BentzFam.lean`
+ `Bentz4.lean` (`Valid9`, `bentz4_of_valid9`, `fileCover`), `lean/LADDER.md`, `notes/lean-bentz-reduction.md`,
`notes/lean-k2m4-reduction.md`.  Existing D4 machinery: `MixedMeasure.lean` (`d4_reduction_measure_u`, `D4InvM`),
`D4.lean`, `Average.lean`, `CovM.lean`.  Certificates: `certificates/k2m3/README.md` "What is not machine-verified" lists
(i) Valid7 itself, (ii) the D4 reduction to the fundamental domain, (iii) Lemma Z's reduction to corner limits — (ii) and
(iii) are your targets.  Python side: `search/qx2_zm.py` (what region/u-range the tilted run covers, the `axis` command =
Lemma Z), `search/QUADRANT_EXACT.md` (Lemma Z statement/proof).

**Part A (do first; should be modest).**  Theorems, generic over a `fileCover` where possible:
`ValidTilt9 → ValidAxis9 → Valid9` (and the same for 7), where `ValidTilt` states exactly the pose region qx2_zm's run
certifies (centre in its D4 root region, u = tan(θ/2) in its range, u > 0, admissible squares) and `ValidAxis` = every
axis-parallel closed unit square in the box has mass ≥ 1.  This needs the box cover's D4 invariance (kernel check on the
data) and the existing reduction lemmas.  Read qx2_zm.py's root/domain definitions carefully so the hypothesis matches
the run exactly (if there's a mismatch between what the run covers and what the reduction needs, that's a finding:
report it, don't paper over it).

**Part B.**  Prove `ValidAxis9` and `ValidAxis7` in Lean (Lemma Z): for axis-parallel squares the mass is a finite sum
of piecewise-(bi)linear terms in the centre (segment fractions, rectangle ∩ Lebesgue square area), so it suffices to
check finitely many cell corners / one-sided limits; design the cleanest kernel-checkable form (a finite check by
`decide +kernel` plus a soundness lemma proved once).  If B turns out large, deliver A plus a precise plan and partial
progress for B.

Rules: kernel reduction only (no native_decide, no new axioms, no sorry in delivered theorems); `#print axioms` = the
standard three; add to `lean/Axioms.lean`; default build must still pass; builds pinned to cores 8–15
(`LEAN_NUM_THREADS ≤ 8`).  Commit on your worktree branch (message ends with the session's Co-Authored-By line); don't push
or merge.  Update `lean/LADDER.md` (a section) and write `notes/lean-valid-split.md`.  Report: exact statements, axioms,
timings, branch + commits, any mismatch found between the Python run's domain and the reduction.
