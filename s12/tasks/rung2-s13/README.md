# rung2-s13: an exact, case-free machine proof of s(13) = 4 (2026-09-11)

**Why.**  The `s(12)` endgame needs zero-margin verification at the container itself (V1 in
`TODO.md`); nothing at `t = 4` can be *established* without it.  Its natural first milestone is
rung 2 of `search/FAMILY.md`: a weighted closed cover of `[0,4]²` with total `W < 13`, certified
exactly by `search/zeromargin.py` at every closed unit square inside `[0,4]²` at every angle.  That
is a proof of `s(13) = 4` with no case analysis (Bentz 2010 needed a 6-leaf tree; DS7's pure sets
cost 14).  `FAMILY.md` §2/§2b records it as **not achieved** and diagnoses exactly why: the scaled
covers `runs/closed4_best_x102.txt` (12.666) and `_x103.txt` (12.790) have provably positive margin
(≥ 6.5 %) at every uncertified box found, but `cert_p1`/`cert_core` converge too slowly at
wall-touching squares (`cx ≈ 0.5`, `cy ≈ 1.5`, `θ → 0` and D4 images) whose witnesses on `x = 1`
sit off-centre at `y = 1.45, 1.55`; and the closing loop `rung2_close.py` has no column generation.

**Semantics** (do not change): `ZEROMARGIN.md` §1 — closed unit squares, closed containment, a
point on `∂Q` counts.  Everything load-bearing in `Fraction`; floats only pre-filter.

**Approach, in order of cheapness.**
1. Take `closed4_best_x103.txt` (or `_x102`) as is and make the checker terminate: (a) bias the
   subdivision to split `u` (angle) first for boxes adjacent to a wall with `θ` near 0
   (§2b's diagnosis: the admissible-width margin grows linearly in `θ` there, so `cx` precision
   needs `θ` precision); (b) an off-centre wall-witness primitive — the closed-form generalisation
   of `ZEROMARGIN.md` §4 item 2 to two points straddling the symmetric line (DS7 Lemma 2's
   `(1,y)`-or-`(1+x,y)` disjunction), *proved on paper in the note* and implemented exactly; (c) a
   two-region (union-of-regions) certification for a box, `box ⊆ R_p ∪ R_q`, exact polygon clipping
   at the bin endpoints plus the sector arcs — §5(ii) of ZEROMARGIN.  Prefer (a) and (b); build (c)
   only if a genuinely tight family appears.
2. If the point set itself is short somewhere (an exact `pose` violation), add points at the
   residual loci (`rung2_close.py --allow-colgen`, stubbed, not implemented) and re-solve; scale by
   a small factor rather than chase the last 0.5 %.  `W < 13` is the only target — 12.9 is fine.
3. Independent check: `zeromargin_stress.py` on the leaf dump (0 failures), and a dense float scan
   of the shipped cover (as `closed4.py stress`) agreeing that the minimum is `>= 1`.

**Definition of done.**  `zeromargin.py` reports **0 uncertified boxes** on the shipped certificate
over the full admissible domain, with the primitive counts; the certificate is checked in
(`certificates/rung2/s13_closed_cover_4.txt` + a short README stating the claim, the command, and
the exact total); the primitives used are each stated as a lemma with a proof in the note.  Do not
modify `verify/`, `verify.sh`, `lean/`, or the existing s(12) certificates.

**Inputs** (read-only, copy into your worktree's `runs/`):
`/home/evand/math/square-packing/s12/runs/inputs-2026-09-11/closed4_best{,_x102,_x103}.txt`.
Code: `search/zeromargin.py`, `zeromargin_stress.py`, `zeromargin_fan.py`, `rung2_close.py`,
`closed4.py`, `scale_cover.py`.  Read `ZEROMARGIN.md` and `FAMILY.md` §2–2b first.

**Budget.**  8 processes for the checker; launch anything over 10 minutes detached (`setsid nohup`);
never `pkill`; pids 921360–921362 are other runs — do not touch.  Report within ~4 h even if not
done, with the exact state (which boxes remain, their location, which primitive would close them).

**Deliverables.**  `search/RUNG2.md` (verdict up front; lemmas with proofs; the checker's leaf
counts by primitive; what is exact vs float; reproduce section), the code, the certificate, commit
on your worktree branch.  Do not edit `TODO.md`, `README.md`, or other tasks' files.
