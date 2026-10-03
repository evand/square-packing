# break-it — empirical attack on L4_k02_box7.txt, 2026-09-29 (report returned as text; saved by coordinator)

**Verdict: no BREAKS, no GAP.**  Scripts/data in this dir (ev.py, axis.py, germscan.py, germfine.py, mutate.py, mut1..6.txt, zm_mut5.log). Cores 10,11, ~1.5 CPU-h (0.4 in one qx2_zm mutant run).

MINOR: tightest pose off the known tight set is near 45°: diamond with bottom vertex (2.9, 1.7) / (3.1, 1.7), μ − 1 = 1.0336e-4 (1800 rational perturbations 1e-3..1e-14 confirm). Positive; ~100× tighter than §1.4's "≥ 1.001" (which is for θ = 0). [Coordinator note: consistent with the LP's tilt margin factor (1 − area(Q∩U)), small for a diamond mostly inside U.]

Lowest exact masses
- Exactly 1: squares inside U; θ = 0 wall families cy = ½⁺, 3/2⁺.
- Tile germs, limit 1 with positive slope: tightest (3/2, 5/2), u < 0, t = (1, 3/5): 1 + 0.31943|u| (same at u = 1e-6, 1e-9, 1e-12, 1e-15). Others: (3/2, 27/10…7/2), (1/2, 5/2…7/2), (3/2,3/2) t = (1,1/5), (3/2, 23/10).
- Walls cy ∈ {½,1,3/2,2} + |u|t, 401 cx, u = ±1e-3..±1e-12: min 1 + 0.3195|u| at cy = 3/2+|u|, cx = 2.9/3.1.
- Rational perturbations (4256 poses, 1e-4..1e-14) around 8 tightest germs: min (μ−1)/|u| = 0.31948.
- Fixed-angle minimisation (Q poking out of U ≥ 0.01/0.05): 0.2° 5.6e-4 (cy≈3/2); 1° 2.8e-3; 3° 4.3e-3; 8° 1.7e-3 (near U corner); 15° 1.8e-3; 25° 1.2e-3; 35° 7.1e-4; 45° 1.03e-4 (diamond).

Mutation tests (D4-symmetric mutants)
- M1 lighten y=1, x∈[2.8,3.0] by 1e-6; M2 y=2, x∈[1.2,1.4] by 1e-6: θ=0 min 1−1e-6 (checker not run).
- M3/M4 move breakpoint (1,1) on y=1 to 0.999/1.001: valid on θ=0 face and 132 germ corners (checker not run).
- **M5** lighten y=1.6, x∈[1.6,1.8] by 1e-4: θ=0 intact; germ (3/2,3/2), t=(1,1/5) = 1 − 2e-4. `qx2_zm.py axis`: min 1 OK (correct). Box checker on cx,cy ∈ [1.49,1.51], u ≤ 0.01: **NOT VERIFIED**, 1 uncertified box [1.5,1.5016]² × [0°,0.112°] — correct rejection.
- M6 = M5 at 1e-7: germ 1 − 2e-7 (checker not run, budget).
No mutant accepted while own evaluator finds sub-1.

CHECKED OK
- ev.py from FORMAT.md only (exact segment clipping, exact area Q∩U); = zm_mixed.exact_mass on 44 poses/cover (L4, REJECTED_L2; incl. u<0); reproduces 0.99993598626 on REJECTED_L2 and germ scan re-finds it unprompted.
- axis.py: 4500 one-sided limit corners over full [½,13/2]² at ±1e-15 + 900 cell midpoints: min 1.
- germscan/germfine: 256 corners × ±u × t ∈ [−3,3]² step 1/5 at u = 1e-12; 27 tight corners refined (1/20 at 1e-15; 1/10 at 1e-6, 1e-9): none < 1; all exact-1 poses inside U.
- Seams x = 2, 5 (D4), corner module, U corner (9/5,9/5), box edge covered by corner grid + angle search.

NOT CHECKED: germ offsets finer than 1/5 at non-tight corners, |t| > 3; germs off the odd/10 grid (except 1-D wall scan); box checker on M1–M4, M6 or other regions; float search heuristic (~460k poses, 16 angles, Nelder–Mead).
