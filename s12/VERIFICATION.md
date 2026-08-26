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
4b. `xcheck.py` (independent Python implementation, exact rationals, EVERY angle bin),
   N=6000 -> minimum 10000023/10000000 over all 2486 bins, 61 s on 32 cores : VERIFIED.
   Per-bin minima are identical to the Rust verifier's in all 2486 bins.
   N=12000 -> same minimum, 4971 bins, 130 s; again identical to Rust in every bin.
   At N=4000 it also REJECTS, with the same minimum 9987038/10000000 in the same bin
   (k=746) as the Rust verifier; one neighbouring bin (k=748) differs because the Rust
   verifier pads its centre range by 1e-6*h (a deliberate superset), i.e. Rust is
   slightly more conservative there.  Sampled mode is retained for speed but a sampled
   VERIFIED is labelled as not a proof.
5. D4 symmetry of the point multiset: checked exactly inside `verify`

## Companion certificate (weaker but human-readable)
`certificates/s12_56points_3.8.txt` — 56 points in [0,19/5]^2 (scales up to [0,1520/397]^2), weight 1/5 each (total 56/5 = 11.2).
Statement: every closed unit square inside [0,3.8]^2, at any angle, contains at least 5 of the
56 points; 12 disjoint squares would need 60.  Checks: `verify` at N=2000/4000/8000, and the
independent exact Python re-check `xcheck.py` over every bin at N=2000 (829 bins, identical
to the Rust per-bin minima) and N=8000 (3314 bins): minimum exactly 1.

## Formalisation
`lean/Sqpack/Basic.lean` (Lean 4 + Mathlib) proves the reduction:
a weighted set whose closed unit squares all carry weight >= 1 bounds the number of squares of
side L>1 packable with disjoint interiors by the total weight; plus the scaling lemmas.
0 sorries; axioms: propext, Classical.choice, Quot.sound.

## Verifier defects found by the rejection tests (2026-08-25), and fixed

Recorded because a verifier's history matters as much as its current state.  None affects
the shipped certificates or the bound; each was found by writing a test that the verifier
should fail, and watching it not fail.

1. **Negative weights were not rejected.**  The reduction `n <= sum w` requires `w >= 0`
   (the Lean proof's hypothesis `hw`).  The old verifier accepted the 56-point set plus one
   point of weight `-2` placed *outside* the container, claimed for `n = 11`: the total
   dropped to 10.8 while the covering was untouched, and it printed VERIFIED.  Now an
   `ERROR`.  The shipped certificates have only positive weights.
2. **Centre box too small above 45°.**  For each angle bin the admissible-centre box used
   the bounding-box width `w = cos θ + sin θ` at `θ_k` only.  Below 45° that is the bin's
   minimum (correct); above 45° `w` decreases, so it was the bin's maximum and a strip of
   width `(w(θ_k) − w(θ_{k+1}))/2` along the container edges was never checked.  Found by
   the exhaustive Python checker disagreeing with the Rust one on a non-symmetric mutant
   (bin k=1530 of N=2000: Rust 2/5, exact 1/5).  Only the `[0,90°)` path was affected; both
   shipped certificates are D4-symmetric and use `[0,45°]` only.  Fixed to
   `min(w(θ_k), w(θ_{k+1}))`.  The tests now pin the per-bin minimum on that mutant.
3. **Panics on malformed input** (header-only file, non-integer token, short point line,
   point count mismatch, empty or missing file, `s_den ∤ s_num·D`).  A panic is not a
   rejection.  All now exit 2 with `ERROR:` and no verdict word.
4. Points outside the container were not rejected; now `ERROR`.

`tests/rejection_tests.sh` now has 42 checks (was 7); against the pre-fix binary 15 fail and
8 panic.
