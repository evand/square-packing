# Verification log

## Result
`s(12) >= 3920/997 = 3.931795386`  — 12 unit squares cannot be packed into any square of side < 3920/997 = 3.9317954.

## Certificate
`certificates/s12_lower_3.931795.txt` — 788 points in [0,3920/997]^2, coordinates with denominator 1994,
weights with denominator 10^7, total weight 14916233/1250000 = 11.9329864 < 12.

## Checks performed
1. `verify` (Rust, exact i128), N=6000  -> min covered weight 10000023/10000000 = 1.000002 : VERIFIED
   (2487 rational angle bins)
2. `verify` (Rust, exact i128), N=12000 -> VERIFIED (4973 bins)
3. At N <= 4000 the verifier REJECTS this certificate: min covered weight
   9987038/10000000 = 0.998704, at theta = 2*arctan(370/2000) ~ 20.96 deg.
   This is the uniform angle net, not the certificate.  The net checks a square shrunk by
   sigma_k ~ 1 - 1/N to cover the gaps between net directions, and below N ~ 6000 that
   shrink exceeds the certificate's slack, so the verifier correctly refuses to certify.
   Recorded here because it is the concrete price of the uniform-net architecture, and the
   reason TODO item A1 (exact adaptive subdivision of pose space, no shrink) is the
   highest-value improvement available.
4. dense float scan, 181 angles over the FULL [0,90] range, ~500k centres each,
   not using the D4 reduction -> minimum 1.000002 (agrees)
5. D4 symmetry of the point multiset: checked exactly inside `verify`

## Companion certificate (weaker but human-readable)
`certificates/s12_56points_3.8.txt` — 56 points in [0,19/5]^2 (scales up to [0,1520/397]^2), weight 1/5 each (total 56/5 = 11.2).
Statement: every closed unit square inside [0,3.8]^2, at any angle, contains at least 5 of the
56 points; 12 disjoint squares would need 60.  Checks: `verify` at N=2000/4000/8000, and the
independent Python Fraction re-check `xcheck.py` (minimum exactly 1 at every sampled bin).

## Formalisation
`lean/Sqpack/Basic.lean` (Lean 4 + Mathlib) proves the reduction:
a weighted set whose closed unit squares all carry weight >= 1 bounds the number of squares of
side L>1 packable with disjoint interiors by the total weight; plus the scaling lemmas.
0 sorries; axioms: propext, Classical.choice, Quot.sound.
