# reduction-coverage — adversarial review, 2026-09-29 (report returned as text; saved by coordinator)

**Verdict: no BREAKS; reduction sound; one GAP (auditability).**  Scripts: data_check.py, general_k.py, cov_check.py, axis_indep.py (+ logs). Cores 8,9, ~25 CPU-min.

1. GAP (auditability, not soundness): runV2_cdade4b6_full.jsonl.gz has per root only root box, leaf-kind counts, boxes, maxdepth, cpu, unc — no leaf boxes/witnesses/tree (dump off, search/qx2_zm.py:1245-1247). Record shows boxes = 2·leaves − 1 for all 9,800 roots, UNCERT 0, unc = [] — but coverage/leaf soundness need a re-run (74,551 CPU-s).
2. MINOR small cases: §0 "known cases k ≤ 7 included". k = 2 false (s(1)=1) → conjecture is k ≥ 3. k = 3 Kearney–Shiu 2002; k = 4, 7 Bentz 2010 (refereed); k = 5, 6 Bentz arXiv:1606.03746 (preprint only). So k = 7 independently refereed, k = 6 preprint; full conjecture needs the preprint for k = 5. README.md:30 and notes/literature-s32.md:24 misattribute k = 3 to Bentz.
3. MINOR §2 step 3(b): "lines x = R, k−R carry extra module mass" is false for this family (ν has no mass on x = 2; seam is π's own). μ_Q shift-invariant on closed {x ≥ R}. Open-interval argument stricter than needed, still correct.
4. MINOR docs: box file format is certificates/s21/FORMAT.md, not certificates/FORMAT.md; §6 reproduce command lacks --resume (no record written); VERIFIED tag at qx2_zm.py:1265 names `qx2_exact.py axis`, §0 names `qx2_zm.py axis` (both run, both min = 1).

CHECKED OK
- Reduction derived independently first: accounting (needs R = w, closed-end seam once per wall), D4 invariance (π mirror, ν diagonal, integer k), localisation (k − 2R > √2 ⇒ k ≥ 6), μ₇ ⇒ μ_Q (integer shift on {x ≥ R} into [0,5]² where μ₇ = μ_Q; needs k = 7), case (d), dilation. Agrees with §2 and QUADRANT.md §1.2/1.3/8.1.
- general_k.py exact k = 6..14: total = k² − 4D, D4-invariant, μ_k = μ_Q on closed [0,k−R]².
- data_check.py: masses > 0; π mirror, ν diagonal; σ = 0, D = 423621306389/500000000000 > 3/4 exact; json agrees; same 800-segment multiset + polygon; box total 49 − 4D; D4 under all 8 maps; sha256s match.
- cov_check.py: 9,800 roots exactly tile [0,7/2]² × u ∈ [0,1/2] (vol 49/8); all leaf kinds certifying; totals match .out; input sha matches; git show 0e303ca → qx2_zm.py sha cdade4b6…a1b7; nothing resumed.
- D4 covers θ ∈ [0°,90°), all centres; SYM never touches θ ≤ 45°; clip_bin, EMPTY, AXIS sound.
- axis_indep.py (whole container, no D4): 3,600 limit corners, min exactly 1.
- Re-run of 128 roots with leaf dump (114 cheap + 14 at 20–120 CPU-s): stats identical; leaves + clip_bin slabs tile roots exactly; no admissible pose in cut slabs.

NOT CHECKED: leaf primitive soundness (E/EXACT0/EXACT45, LEB, CAP, PIECE, polygon S(b)); θ → 0⁺ germ limits; 9,672 roots not re-run (incl. all > 120 CPU-s); Bentz's preprint; runs A, B.
