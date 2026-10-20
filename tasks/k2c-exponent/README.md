# Exponent of c*(k) from above: what blocks 3/8 (opened 2026-10-19)

Evan, 10-19: "What is the next route to asymptotic improvements?"  Upper bounds: Ryu's chain gives c*(k) = O(k^{3/8})
(k2-minus-c-upper v1.1 §8); reported lower bound k^{1/3} ([T], `tasks/k2c-lower/`).  Ryu's paper (v1.1, fetched
10-19 to scratch, data only): §7 staircase, §8 wall filler ZC (Lemma W), §12 open problems ("no claim that 3/8 is the
best exponent obtainable from this frame").

## The frame

c*(k) ~ 4k/b + E(b) (Prop. 6.3).  E(b) = T(b, 1/b): the band end, a thin trapezoid of length b and ceiling slope
~2/b (height free, via the lift).  T(L, z) is filled by
* trivial axis layers (cost L if Lz >= 1);
* a chain of staircase blocks (Ryu's ZC; one block = the §7 staircase): columns of height X at tilt t,
  Kb = max(1, Lz/(t²X)) blocks.  Costs (Lemma W's terms, orders): column triangles Lt; aligned closings Kb·X·t;
  unaligned closings Kb/t²; transition slivers Kb·h_T·w_T; small wedges; and **the two ends of the region, which are
  walls T(X, t)**: so the frame is recursive.
* constraints: X >= 1/t² (W3); for Kb > 1: X >= Lt and L <= tX² (W2 drift), t >= z (W1); for Kb = 1: t²X >= Lz.

Wall exponent γ (wall of height m costs m^γ) gives c* exponent 2γ/(4γ + 1): γ = 1 → 2/5, 3/4 → 3/8, 1/2 → 1/3.

## Calculator: `search/ryu/exponents.py`

Every term is a monomial, so for each choice of filler per level (S = staircase, C = chain) and leaf, the min-max
exponent is an LP; the script enumerates all choices up to a depth.  Orders only (constants, logs dropped); the
sliver and wedge terms are my reading of Lemma W's B_W, not Ryu's own bookkeeping.

**Reproduces Ryu:** depth 0 → 1/2, 1 (S) → 2/5, 2 (S, C) → 3/8.  **Deeper nesting gives nothing** (3/8 at depths
3-6): the binding terms are the wall's (level 1) column triangles, aligned closings and unaligned closings, all at
m^{3/4}.  So the "tower" guessed on 10-19 (n/(3n−1)) does not come from repeating Ryu's pieces.

What-if runs (`--drop` removes terms, `--scale` multiplies them by t^E), best exponent over depth ≤ 6:

| hypothetical change | exponent |
|---|---|
| none | 3/8 = 0.375 |
| unaligned closings free | 3/8 (aligned + triangles still balance at L√z) |
| aligned closings free | 13/35 ≈ 0.3714 (then unaligned binds) |
| both closings free, or both × t | 27/76 ≈ 0.3553 (then the transition slivers bind) |
| closings and slivers free | 4/11, 8/23, 16/47, …, 64/191 = 2ⁿ/(3·2ⁿ − 1) → **1/3** |
| column triangles free | below 1/3 (contradicts [T]: column triangles are a real cost) |

## Reading

* The obstruction to beating 3/8 inside Ryu's frame is **the cost of changing tilt between blocks** in the chained
  wall filler (closings and slivers), not the ends, not the depth.  If block transitions cost no more than the column
  triangles, the frame reaches k^{1/3+o(1)}, the reported lower bound's exponent: so 1/3 is plausibly the true
  exponent, and the work is a cheaper tilt-change primitive.
* Any reduction of only one transition term moves the exponent a little (0.371, 0.355); all of them must go together.

## Next

1. A cheaper tilt change: a transition between columns at tilt t and t − κ whose waste is O(column triangles of the
   region it spans).  Ideas to test in the model first (add a filler type): gradual per-column tilt change with the
   gap filled (Ryu's §12 remark on tilted columns: each change costs a gap of order M tan β, which gives only k^{3/7}); transitions placed where the layers
   end anyway; filling the transition sliver by a sub-wall (recursive).
2. Local lower bound: must a tilt change of κ across a column of height X waste ≳ X·κ-ish?  If yes, the frame is capped
   and a different end design is needed; compare with [T]'s path-inclination argument.
3. Check the model's sliver / wedge terms against Lemma W's proof (they decide 0.355 vs 1/3 in the what-ifs).
