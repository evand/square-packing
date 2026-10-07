# Unavoidable point sets for closed unit squares

A *pure unavoidable set* for `[0,m]²` is a finite point set that every closed unit square inside `[0,m]²`
contains.  Background, method and soundness: [`notes/unavoid13.md`](../../notes/unavoid13.md).

| file | what it certifies | checker |
|---|---|---|
| `ks7_rational_3.txt` | the 7 points `(9/10,1), (3/2,1), (21/10,1), (3/2,3/2), (9/10,2), (3/2,2), (21/10,2)` are unavoidable for `[0,3]²`: **`p(3) ≤ 7`** | `python3 search/unavoid13_check.py cert certificates/unavoid13/ks7_rational_3.txt --tri --seg --depth 16` |
| `friedman14_4.txt` | Friedman's 14 points are unavoidable for `[0,4]²` (sanity check) | `python3 search/unavoid13_check.py cert certificates/unavoid13/friedman14_4.txt --tri --seg --depth 14` |
| `unavoid3_lower7_family.txt` + `unavoid3_lower7_bb.txt` | no 6 points meet all 87 closed unit squares of the family: **`p(3) ≥ 7`**, so `p(3) = 7` | `python3 search/unavoid13_exactcheck.py certificates/unavoid13/unavoid3_lower7_family.txt certificates/unavoid13/unavoid3_lower7_bb.txt` (~1 s) |
| `unavoid4_C2_family.txt` + `unavoid4_C2_bb.txt` | no set of at most 13 points invariant under the half-turn `p ↦ (4,4) − p` meets all 441 squares of the family (221 half-turn orbits): **no 13-point unavoidable set for `[0,4]²` is half-turn symmetric** (hence none is `C4`-, `V`-, diagonal-Klein- or `D4`-symmetric: each of these groups contains the half-turn) | `python3 search/unavoid13_exactcheck.py certificates/unavoid13/unavoid4_C2_family.txt certificates/unavoid13/unavoid4_C2_bb.txt` (~25 s at 4 threads) |

**Families** (`*_family.txt`): `# m = …` header, then one `pose u cx cy` line per square, rationals; `u = tan(θ/2)`, so
`cos θ = (1−u²)/(1+u²)`, `sin θ = 2u/(1+u²)` are rational; the square is closed and contains its boundary.  Every
square is admissible (checked exactly).  The 87-square family is a greedy shrink of the 575-square family of
`notes/unavoid13.md` §2.3 (which stays in the research runs); the 441-square family is the round-1 family of the
half-turn-symmetric loop (§4.2).

**Branch-and-bound certificates** (`*_bb.txt`, produced by `search/unavoid13_exactcert.py`; the full format is the
docstring of `search/unavoid13_exactcheck.py`): the claim (`plain` or `C2`), `k`, a list of candidate points (exact
rationals), then a preorder tree over the 0/1 candidate variables.  `B j` splits on candidate `j` (`x_j = 1` subtree
first, then `x_j = 0`); a leaf is either `U i` (row `i` has no candidate left) or `D den i:num …`, rational LP duals
`y ≥ 0` on the rows whose weak-duality bound `L = Σy + Σ_{fixed 1} r_j + Σ_{free} min(0, r_j)`, `r = c − Aᵀy`, rules
out every integer point of the leaf with cost `≤ k`.  Rounding: integer solutions have cost `C1 + g·t`
(`C1` = cost of the fixed-1 candidates, `g` = gcd of the free costs, `t ≥ 0` integer); a leaf closes iff the least such
value `≥ L` exceeds `k`.  For the plain claim that is `L > 6`; for C2 (orbit costs 1 for the centre, 2 otherwise) the
parity of `C1 + 2t` lets a leaf close with `L > 12`.

The checker recomputes everything else from the family file alone: the exact arrangement vertices (square corners,
all intersections of non-parallel edges of two squares, and the centre for C2), exact incidence (Fractions), and the
domination condition that makes the candidate list sufficient (every vertex's column is contained in some candidate's
column, of no larger cost).  It shares no code with the exploration modules that produced the instance.  Rejection
tests: `tests/unavoid13/rejection_tests.sh`.

| certificate | candidates | rows | vertices (distinct columns) | root LP | branchings / leaves | min `L − k` |
|---|---|---|---|---|---|---|
| `unavoid3_lower7_bb.txt` | 233 | 87 | 6,112 (3,969) | 5.5 | 671 / 672 | 0.0227 |
| `unavoid4_C2_bb.txt` | 877 (centre + 876 orbits) | 441 | 55,717 (20,104) | 12.07 | 16 / 17 | −0.90 (parity) |
