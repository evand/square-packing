# lemmaE-math — adversarial review of Lemma E mathematics (§4.3–4.4), 2026-09-29
(Report returned as text by the reviewer; saved by the coordinating session.)

**Verdict: no BREAKS, no GAP.**  4 MINOR.  0 violations in exact tests.  Scratch: geom.py, t_caps.py, t_lemmaE.py, run10.log, run11.log.

MINOR
1. "f continuous on R(u)" false in corner3/xcut regimes (§4.4 statement + proof end): on z = 0 the z>0 cell function ≤ the z<0 one. Proof works per closed cell with min of the two at vertices on z = 0; code does exactly that (one side only when sign of z strictly certified, rf_positive, qx2_zm.py:714–735). Fix wording: "continuous on each closed cell".
2. McCormick loss bound: tangent term loses ≤ u·h_y²/(4(1−u²)), not u·h_y²/4 (≤1.21×). Unused for soundness.
3. Anchor naming: 'hi' anchor uses box low corner (cx0,cy0) (qx2_zm.py:516). Math consistent.
4. Cap lines d=0, d=s unnecessary for concavity (1−ĝ concave everywhere); harmless.

CHECKED OK
- §4.3 chord options: all 8 rederived; tests use polygon-edge chords.
- E′: g vs exact clipping; ĝ≥g, parabola≥ĝ, tan-mode ≤ true area at 4,000 poses (tiny tilts, u within 1e-5 of tan 22.5°). ĝ−g=(c−d)²/(2sc)≥0 on [c,s+c], valid through s=c. Tan mode needs d≥c: enforced box-wide. "Φ=1 for never-crossed line" is box-wide.
- E″ corner3 & xcut: independent derivation, 20,000 poses vs clipping. z≤0 ∧ α≤C ⇒ β≤S. xcut −g_x+A_LL = quad(α′,β″) symbolically.
- McCormick both anchors, tangent to Sβ², −Sα² kept: correct. Code RF encodings of α, β, α′, β″, z, z′, cap lines, tangent match.
- Lemma E end-to-end (own impl, fixed rational u; random boxes near U corner/walls, random profiles w/ good+bad kinks, all regimes, both anchors): cell fn ≤ true mass; point value ≥ min over own closed-cell vertices (exercises omitted good breakpoints); concavity midpoints: 12,000 trials, 0 failures. Mutations (drop bad-bp lines / pair lines / z-line, swap bad/good) do fail → test has power.
- Proof of E: bad-breakpoint/concavity argument correct; cells bounded (R(u) sides in 𝓛); degenerate R(u) via 1-D version.
- S-procedure and gcd certificates sound (every LB ingredient valid wherever vertex ∈ R(u) ⊂ box; "outside" certified strictly).
- u-range (u0,u_e]; u0=0 fine; EXACT45+SYM not circular (tan 22.5° irrational); u0>0 endpoints covered by neighbour's closed top.

NOT CHECKED
- Code vs statement beyond cited lines (lines_in_reach float margin, piece/dominance pruning, clip_bin interaction).
- Actual run / leaves; Lemmas Z/U/K; zm_mixed; §2 reduction.
- Tests on random data at fixed u, not the real cover's lines.
