# zm-polygon — adversarial review of zm_mixed polygon path (Lemma S(b)) as used by QXChecker, 2026-09-29
(Report returned as text; saved by coordinator.)

**Verdict: no BREAKS, no GAP.**  4 MINOR.  Scratch: test_condpoly.py, test_poly.py, test_leaves.py, sweep.sh, test_convex.py, probe.py, *.out. Cores 6,7, < 1 CPU-h.

MINOR
1. ZM_MIXED.md S(b) / QUADRANT_EXACT §4.5 present "F ⊇ Q so K∩F loses nothing" as part of the argument; soundness doesn't depend on the frame (every point of K satisfies all four inequalities at every admissible pose regardless of F). Claim true anyway (√2/2 < 0.7072).
2. Float `near`/`window` filters (zm_mixed.py:1108) only drop pieces → lower bound; safe.
3. 68/800 segments lie on ∂U, none with midpoint in open U; separate 1-D measures, Lebesgue gives lines 0 → no double counting.
4. §4.5 says no leaf dump stress-tested; sampled one here (below).

CHECKED OK
- `_slots` worst-corner table (Lemma A) re-derived θ ∈ [0°,90°), exact test 2,000 cases × 4 inequalities; admissible centres at fixed u a product of intervals → R/W/M bounds valid. piece_bound asserts 0 ≤ u0 ≤ u1 < 1.
- cond_poly = N(u)·inequality at chosen centre bound (N = 2(1+u²)², or 2(1+u²) when both R): 15,000 exact checks, all choices/inequalities, u ∈ [0,1].
- Affine in p: matches 3-point reconstruction exactly for every choice; both_rect depends only on choice.
- bern: true Bernstein coefficients incl. h = 0.
- Hull over choices ⊂ C_k (C_k convex intersection of half-planes; each per-choice region ⊂ C_k); dropping degenerate/empty regions only shrinks K.
- Orientation/clipping: frame CCW, convex_hull CCW, S–H preserves; zero-length edges harmless; clip_convex = own clipper exactly on all test boxes; check_simple_convex rejects CW/doubly wound/repeated/collinear. For U, w = area = 289/25 → w/area(P) = 1.
- test_poly.py: 1,686 exact boxes on box7 cover (straddling U sides/corner, walls; zero-width centre ranges, u0 = 0, widths to 1/5000); polygon part = piece_bound − piece_bound(no polygon); 188,408 admissible rational poses: 0 bound > area(Q∩U), 0 K-vertices outside Q, 0 total piece bound > independent μ(Q). Sharp (gap 0) where U's corner cap ⊂ every Q.
- Certificate leaves: QXChecker with V2 args on 172 random roots in straddle band cx ∈ [1.1,2.5], cy ∈ [1.1,3.5]: 462 PIECE leaves with positive polygon part, 25,872 poses, 0 bad, min sampled μ 1.014.
- Corollary T′ unreachable (0 points; QXChecker.run_box never passes with_pts).
- cert-mode gates only polygon branch and T′; nothing else unlocked. Polygon part enters cert_split via `rest` (valid box-wide); SPLIT certified 0 leaves in V2.

NOT CHECKED: Lemmas T, L, V, R and line cores (beyond total-mass comparison); all 15,926 PIECE leaves (6 sampled roots hit 120 s limit); LEB, CAP, EXACT, AXIS, SYM, D4, reduction.
