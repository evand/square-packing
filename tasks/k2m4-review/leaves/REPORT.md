# k2m4-review / leaves — REPORT (transcribed by the session from the agent's final message; subagents could not write report files)

No BREAKS, no GAP, 2 MINOR.  ~2.5 CPU-h (core 6).  Scripts/outputs in this directory.

MINOR
1. Leaf record not self-contained: `qxzm_full.jsonl` stores post-clip_bin leaves; the 3,033 clip_bin slabs and the gaps left
   by 2,805 zero-u-thickness EMPTY leaves are not recorded; in 10 roots the gaps are mid-u-range (e.g.
   [3/5,7/10]×[13/10,7/5]×[5/16,3/8]), so the tree can't be rebuilt by midpoint splits alone.  An auditor must treat empty
   nodes as "no admissible pose": reconstructed all, every one inadmissible.  Document.  (qx2_records.py does check this:
   "3043 uncovered parts hold no admissible pose".)
2. K2M4_MARGIN §4: "roots 5,400–12,000" are completion counts, not a region.  The cost is at cx ∈ [2.4, 2.6], cy ∈ [2.5, 3.3]
   (squares straddling a = 14/5 and the doubled line 3) and u-bin [0, 1/16] (53 % of CPU).  The three slowest roots have true
   min mass ≈ 1.0036: slow because of the bounds, not a tight pose.

Checked OK
* Box file (own parser): total 81 − 4D exactly < 77; U mass = area; 16 doubled segments on x, y ∈ {3, 6}; line measure
  D4-invariant (8 maps); overlaps summed in Profile, LineMass, zm_mixed's evaluator; shas match.
* Coverage (own guillotine reconstruction): roots = 45×45×8 grid on [0,9/2]² × u ∈ [0,1/2]; per-root stats = leaf counts =
  .out totals; every root exactly covered by leaves + inadmissible gaps (exact, h(u) rational); 3,890 EMPTY inadmissible; leaf
  kinds in their u-ranges.
* Reduction: D4 + θ ↦ π/2 − θ reach every Valid9 pose; SYM/EXACT45 not circular; θ = 45° never inside a SYM box.
* Lemma E at box 9: m = 9 used, u ≤ 1/2 enforced, trig constants valid; u_regime's "left/below 31/5" holds (reach ≤ 5.21);
  exact-from 3 doesn't affect soundness.  sha 6294052a vs reviewed cdade4b6: only an assert + dump flag.
* θ = 0 (own Lemma Z, full face): 6,400 one-sided limits, min exactly 1, 832 tight, no USC violations; all 832 tight corners
  probed with tilts ±|u| ∈ [1e-12, 1/200] and offsets: 0 below 1.
* Exact falsification: ~0.9 M exact masses (all EXACT0/EXACT45/CAP/LEB leaves + 6,000 EXACT + 6,000 PIECE, half from slow
  roots; wall, tiny-tilt, local search): 0 below 1.  Minima: EXACT/CAP/LEB exactly 1 (inside U), EXACT45 1 + 4.1e-5, PIECE 1.0013.
* Independent proofs: all 683 LEB + 393 CAP re-proved (own code; Lemma K re-derived, triangle-cap chord range holds past
  45°); own common-rectangle branch-and-bound: PIECE 60/60, EXACT 53/60, EXACT45 53/60, EXACT0 11/56 (rest: budget, crude
  bound, nothing suspicious).
* Checker's primitives on recorded leaves: 200/200 PIECE bounds ≥ 1 and ≤ sampled exact mass; 12 EXACT, 5 EXACT0, 60 EXACT45
  re-certify; target raised to sampled min + 1e-9 fails in all 44 applicable cases (bounds tight, not vacuous).
Not checked: Lemma E / zm_mixed internals; independent proof of most EXACT leaves; only 77 EXACT-type leaves re-run; no
full root re-run; claim chain (area claim); germs beyond tight corners (area break-it).
