# Adversarial review: s(20) > 3 + 4√2/3 (2026-10-04)

Claim: `search/S20_LB.md`; commits b8157df (zmx2), 442d9bc (line_cover + s20lb tools), 2010d63 (cover + write-up);
reviewed at HEAD 305e8f5.  Reviewer: an independent agent (no part in producing the claim), saved here by the
session lead because the agent could not write files.

**Overall verdict: the claim holds; no soundness defect found.**  Three wrong statements in S20_LB.md (fixed with
this report) and two gaps not hit by this run.

| # | item | verdict |
|---|---|---|
| 1 | reduction at a non-integer side, strict vs non-strict | OK |
| 2a | zmx2 root superset | OK |
| 2b | mirrored grid-line atoms `x = s − k` | OK (soundness); doc defect in S20_LB.md; ZMX2.md not updated |
| 2c | other integer / 1/10 side assumptions in zmx2 | OK |
| 3 | zm_mixed pitch and settings | OK; latent gap in `zeromargin.roots()` |
| 4 | independent parse of the cover | OK; "10,220 distinct" wrong |
| 5 | reproduction | OK, every census box for box |
| 6 | mutations | OK, every mutant proved invalid is refused |

## 1. Reduction
Rescaling works at any real side: a packing in side `s' < t` scaled by `t/s' > 1` makes each square strictly contain a
concentric closed unit square inside `[0,t]²`, pairwise disjoint, so `20 ≤ total < 20`.  Closed squares match zmx2
(`|X|,|Y| ≤ ½`, closed admissibility) and zm_mixed.  The certificate gives `s(20) ≥ t`; strictness against the
constant comes from `t` itself: `(t−3)·3/4 = 2829/2000`, `2829² = 8,003,241 > 8,000,000 = 2·2000²`.  Points only, so
Lean `packing_le_weight` applies as it stands.

## 2. zmx2 change (b8157df)
* **2a roots.**  `[0, 2.5]` under `--d4` (⌈5s⌉ = 25 cells) and `[0, 4.9]` under `--full` (49): supersets of the required
  regions.  D4 maps and the `--full` y-reflection use the exact `sx = 4886`.  Extra poses are harmless: the bound need
  only hold at admissible poses (genuine container poses); Lemma E uses the exact `s`; candidate buckets can only drop
  points.  Verdict logic unchanged; roots 2,500 (`--d4`), 19,208 (`--full`).
* **2b atoms.**  The code gives lines `{1, 1.886, 2, 2.886, 3, 3.886}` (not `{1, 2, 2.886, 3.886}` as S20_LB.md §1.2 said);
  the cover has 72 points each on `x = 3` and `x = 1.886`.  Line lemmas: C/M/S/H take an arbitrary position ℓ; Z pairs
  lines only at an exact integer difference (chains grouped mod 1); W requires exact distance 1 from a wall using the
  wall at `sx` (right-wall form re-derived: `f1(1 − (C+S)/2) = −q(u)`); DP and A are position-free.  Gap: ZMX2.md §1,
  §4.5, §7 still describe the old roots / line set.
* **2c.**  `check_d4`, `reflect_y` use `sx − x` exactly (units of 1/1000, no parity issue); clamp and box denominators
  are side-independent; only `cert0` still requires a 1/10 side (unused).  Observation: zmx2 is not monotone in the
  side (the same points declared at 2442/500 are a valid cover but `--full --sym-atoms` refuses with 12,392
  uncertified, because atom lines and Lemma W move): completeness only.

## 3. zm_mixed
HEAD `zm_mixed.py` sha `1fd20346…` = the shipped copies in certificates/{s21,s45,s60}; `zeromargin.py`,
`mixed_cover.py` unchanged.  Settings = s(21)'s `verify.sh` except the pitch; with 2443/50000, `(s/2)/pitch = 50` and
`s/pitch = 100` exactly.  No integer-side assumption (`m` a Fraction, `S` exact by `validate`).  **Latent gap (not hit):**
`zeromargin.roots()` (used by `--full` and the default mode) computes `int(m/pitch)` with no divisibility assert; a
`--full` run with a non-dividing pitch would skip a strip of centres and still print VERIFIED.

## 4. Cover file (own parser, exact Fractions)
`s = 2443/500`, `D = 1000`, `W = 1e8`, 12,864 points, all weights > 0, all inside `[0, 4886]²` (x in 0.25…4.636); total
`249862891/12500000 = 19.98903128 < 20`; exactly D4-invariant (all 7 non-identity elements); full sha256 matches.
S20_LB.md §0 errors: all 12,864 positions are distinct (not 10,220); the short hash mixed in `mixed_cover.py`'s.

## 5. Reproduction (rebuilt from `git archive HEAD` in scratch)
| run | verdict | boxes |
|---|---|---|
| zmx2 `--d4 --sym-atoms` | VERIFIED-D4 | 683,684 |
| zmx2 `--d4` | VERIFIED-D4 | 689,800 |
| zmx2 `--full` | VERIFIED | 4,570,032 |
| zmx2 `--full --sym-atoms` | VERIFIED | 4,525,820 |
| zm_mixed `--d4 --cert-mode` | VERIFIED-D4 | 67,766 (leaf counts identical) |
| zm_mixed `--full --cert-mode` | VERIFIED | 545,216 (11,883 CPU-s) |

Own spot checker (float search + local descent, then exact Fraction closed containment and admissibility): ~40,000
starts, 0 exact failures.  Lowest exact masses: random 1.0106; wall and θ→0⁺ germs (incl. right wall, `x = 3.886`)
1.019005; near the reported worst box 1.000380 at (1.31189, 2.48480, 27.318°), inside the bisection's implied range
[1.00036, 1.00048).

## 6. Mutations (zmx2 `--d4 --sym-atoms` unless noted)
| mutant | result | exact witness |
|---|---|---|
| weights × 0.998 | refused, 275 uncertified (zm_mixed: 438) | μ = 0.99837 |
| heaviest orbit (0.6, 1.0) dropped | refused, 7,474 | μ = 0.8998 |
| heaviest orbit moved 1/1000 inward | refused, 1,046 | μ = 0.9594 (wall germ) |
| all points on x, y ∈ {1, s−1} moved off the line | refused, 6,307 | μ = 0.5406 |
| side 2444/500 | `--d4` refuses (not D4-invariant); `--full` 27,908 | μ = 0.7104 (right wall) |
| points on x = s−1 only, × 0.90 / 0.97 / 0.99 (`--full`, b8157df's new path) | refused, 20,480 / 1,180 / 152 | × 0.90: μ = 0.9712 |
| median-weight orbit or a 2e-8 orbit dropped | verifies (plausibly still valid) | — |
| orbit (1.0, 2.1) dropped | refused | none found (may be valid but unprovable) |

## Recommendations (none blocking)
1. Fix S20_LB.md: line set (§1.2), "distinct" count and short hash (§0).  **Done with this report.**
2. Document b8157df in ZMX2.md §1, §4.5, §7.  **Header note added with this report.**
3. Divisibility assert in `zeromargin.roots()` as `d4_roots` has.  **Added to the zeromargin should-fixes TODO item**
   (it changes the checker sha, so it goes with the bundle re-pin).

Budget: ~9 CPU-h (Python spot sampler ~4.5 h, zm_mixed `--full` 3.3 h).
