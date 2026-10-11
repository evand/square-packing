# nagqe: Nagamochi's Lemma 1 as an exact real-arithmetic (CAD) decision problem  (2026-10-10)

Code: `search/nagqe/` (`nagqe.py` exact evaluator + z3 encodings, `test_encoding.py`, `decide.py`, `nlsat.py`,
`prune.py`, `lifted.py`, `cvc5_run.py`).  Venv `runs/qe_venv/` (z3-solver 4.x "5.1.0", cvc5 1.4.2).  Labels as
elsewhere: **[proved]** exact, **[measured]** solver run as reported, **[heuristic]**.

## Question

Nagamochi (EJC 2005) Lemma 1, in chelokot's formalisation (`ScoresDilatedSquares`): in `[0,m]²`, every square of
side `1 < λ ≤ 1.01` that fits scores `> 1` against the resource (inner-box area; lines `x, y ∈ {1, m−1}` over
`[1−e, m−1+e]` at density `w`; eight Q points at the segment ends, weight `q`; P points `(i, 1−d)` etc., weight `p`).
Nagamochi: `e = d = 1/10, w = p = 1/2, q = 9/20`; false (chelokot's 0.977543…).  chelokot also proved the `q = 1/2`
layout valid (total `n² − 1.6`) and ruled out repairs that change only the four weights at fixed positions.
Here: test whether off-the-shelf CAD-based solvers decide the single-square lemma exactly, then (if so) ask the
parametric question with the positions `e, d` free.

## Results

1. **Exact evaluator** (`score_exact`, Fractions, area by polygon clipping) reproduces chelokot's counterexample
   score `25009470849041/25584000000000` exactly [proved].
2. **Encodings** (with divisions; and division-free, scaled by `c²s²`) agree with the evaluator exactly at 450 + 450
   random rational poses, `m = 4, 5, 6`, Nagamochi's and a perturbed layout; 0 mismatches [proved at those poses].
3. **Axis-parallel slice** (`s = 0`): `unsat` in 0.02 s for the original layout, `m = 4` [measured].
4. **Tilted, whole pose space** (5 variables `X, Y, c, s, λ`, `m = 4`): **no verdict** from any engine within
   5–10 min, for either the original (true answer: sat) or the `q = 1/2` layout (true answer: unsat) [measured]:
   z3 default solver (unknown, 363 s); z3 `qfnra-nlsat` after `elim-term-ite` (canceled, 300 s and 600 s);
   z3 `qfnra-nlsat` after `cofactor-term-ite` on the division-free encoding (canceled, 600 s, both layouts);
   cvc5 cylindrical algebraic coverings (killed at 320 s).
5. **Box restriction + exact branch pruning** (`prune.py`: every if-then-else condition decided by two small NLSAT
   calls on the box, decided ones substituted): a corner box (`X, Y ∈ [1, 1.05]`, `s ∈ [0.30, 0.31]`, `q = 1/2`)
   prunes to 0 branches and is `unsat` in 0.38 s; the box around chelokot's counterexample (`X ∈ [0.5, 0.55]`,
   `Y ∈ [1.45, 1.5]`, `s ∈ [0.04, 0.06]`) keeps **11 undecided branches** and times out (200 s) for both layouts.

**Reading.**  "5 variables" is misleading: the score is piecewise (min/max for every chord, clipping for the area,
an indicator per point), and every undecided piece becomes either a fresh CAD variable (`elim-term-ite`, division
purification) or a Boolean case split over 5 variables with high-degree atoms (`cofactor`).  Near the
counterexample — exactly where the answer is decided — a dozen pieces stay live even in a tiny box.  Whole-problem
CAD does not decide this one-square lemma in minutes; the cheap regions are decided instantly.  This matches the
general picture (CAD cost is in the piecewise structure, not just the variable count).

## Hand analysis of the parametric family  [heuristic: limits of explicit pose families, not yet exact]

Limit constraints as `λ → 1`: mid-wall square holding one P point, `w + p ≥ 1`; chelokot's slip (catches the Q
point at height `1−e`, misses the P point), `q + w ≥ 1`; corner square, `2q + 2we ≥ 1`.  With the `m = 4` budget
`8w + 8we + 8q + 4p ≤ 10`, the LP minimum of the left side is `12 − 2(1−2e)/(1−e)`, above 10 for every `e > 0`;
at `e = 0` the Q points merge at `(1,1)` and a slip square can avoid them.  So the 5-parameter layout family
`(e, d, w, q, p)` is expected to be **irreparable** at total `n² − 2`, with margin `≈ 2e` as `e → 0`.  Making this
exact needs explicit pose families with `λ > 1` (each an exact rational evaluation) and a separate argument as
`e → 0`.

## Next steps (not started)

* Exact version of the irreparability claim: rational pose families + LP duality per `e`-interval; the small-`e`
  end is the one genuinely parametric piece.
* If a decision procedure is still wanted: hand cell decomposition (fix which points are inside and which bounds
  are active per cell, then each cell is a pure polynomial system in 5 variables), or a dedicated CAD (Tarski /
  QEPCAD B via nix) on the per-cell systems.
