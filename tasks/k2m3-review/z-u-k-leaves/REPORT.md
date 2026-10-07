# z-u-k-leaves — adversarial review, 2026-09-29 (report returned as text; saved by coordinator)

**Verdict: no BREAKS, no GAP.**  Scratch: myz.py, ev.py, replay scripts. Cores 4,5, ~30 CPU-min.

MINOR
1. `qx2_zm.py axis` (__main__) restricts to [½,7/2]² on D4 grounds without itself checking D4 symmetry (only main() does). Own check covers full [½,13/2]² without D4, so no dependence.
2. main() final message says "u = 0 face: qx2_exact.py axis" vs claim of record `qx2_zm.py axis` (wording; all three give min 1).
3. Lemma K proof remark "lines y = a, x = a meet in one point" unnecessary (U and the segments are separate parts; segments have no atoms). Code condition cy₀ ≥ a stronger than needed (safe).
4. Run records keep only per-root totals; auditing a leaf needs a re-run of its root.

CHECKED OK
- Lemma Z proof re-derived (multilinear on cells, usc from closed squares, grid-line values ≥ neighbour limits); Γ complete; axis_face one-sided inclusion correct.
- myz.py: own parser (box sha c0a67509…), full [½,13/2]² no D4, limits by exact extrapolation: 3600 limits, min exactly 1, 752 tight, 0 below 1.
- Lemmas U, K: width lemma, mass bound, both triangle-cap chord ranges re-derived; θ₀/θ₁ endpoint choice correct; gaps = density 0.
- LEB/CAP stress: ~1,200 accepted boxes, ~90k exact poses near y=a, x=a, (a,a): lemma conditions held, all masses ≥ 1.
- clip_bin/EMPTY/AXIS: 2.6M exact tests, no admissible pose dropped.
- Dispatcher order EMPTY, AXIS, SYM, LEB, CAP, zm_mixed, EXACT: non-overlapping, each sound. SYM/EXACT45 θ>45° → θ′<45° where no SYM leaf: no circularity. D4 domain right. AXIS = boxes touching cx=½ or cy=½ (Lemma Z).
- Leaf replay: 4,763 cheap roots re-run with dumps, totals match V2 exactly; checked 61/61 AXIS, 2,955/2,955 EMPTY, 683/689 LEB, ~280/374 CAP, ~845 SYM, some EXACT45/PIECE/EXACT: 0 bad.
- V2 record: exactly 9,800 D4 roots, header shas match, 0 uncertified.
- Tiny tilts (ev.py exact): ~1.1M poses (grid corners at u = 1e-4, 1e-8, 1e-12; x = a at ±1e-8, ±1e-4; 64k around (a,a)). Tightest non-U: wall squares cy = 3/2⁺, gain ≈ 0.16·θ; near x = a margin ≥ 6e-4. None below 1.
- Only primitive enabled beyond --cert-mode: polygon Lemma S(b), feeding PIECE (15,926 leaves, the only zm_mixed leaf type used). Point-based rule unreachable (no points). Independent test: 4,060 boxes, 48k poses, credited region ⊂ Q: 0 failures.

NOT CHECKED: Lemma E / Exact.certify (EXACT leaves); zm_mixed line-lemma internals beyond reading; qx2_exact.py axis (superseded); §2 reduction; ~95 CAP + 6 LEB leaves in expensive roots.
