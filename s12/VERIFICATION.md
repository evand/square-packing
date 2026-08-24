# Verification log

## Result
`s(12) >= 3920/997 = 3.931795386`  — 12 unit squares cannot be packed into any square of side < 3920/997 = 3.9317954.

## Certificate
`runs/cert3932.txt` — 788 points in [0,3920/997]^2, coordinates with denominator 1994,
weights with denominator 10^7, total weight 14916233/1250000 = 11.9329864 < 12.

## Checks performed
1. `verify` (Rust, exact i128), N=2000  -> min covered weight 10000023/10000000 = 1.000002 : VERIFIED
2. `verify` (Rust, exact i128), N=6000  -> VERIFIED  (2487 rational angle bins)
3. `verify` (Rust, exact i128), N=12000 -> VERIFIED (4973 bins)
4. dense float scan, 181 angles over the FULL [0,90] range, ~500k centres each,
   not using the D4 reduction -> minimum 1.000002 (agrees)
4. D4 symmetry of the point multiset: checked exactly inside `verify`

## Companion certificate (weaker but human-readable)
`runs/cert38_clean.txt` — 56 points in [0,19/5]^2 (scales up to [0,1520/397]^2), weight 1/5 each (total 56/5 = 11.2).
Statement: every closed unit square inside [0,3.8]^2, at any angle, contains at least 5 of the
56 points; 12 disjoint squares would need 60.  Checks: `verify` at N=2000/4000/8000, and the
independent Python Fraction re-check `xcheck.py` (minimum exactly 1 at every sampled bin).

## Formalisation
`lean/Sqpack/Basic.lean` (Lean 4 + Mathlib) proves the reduction:
a weighted set whose closed unit squares all carry weight >= 1 bounds the number of squares of
side L>1 packable with disjoint interiors by the total weight; plus the scaling lemmas.
0 sorries; axioms: propext, Classical.choice, Quot.sound.
