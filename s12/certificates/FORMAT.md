# Certificate format

A certificate is a plain-text file of whitespace-separated integers. Everything is exact:
no floating point appears anywhere in a certificate or in its verification.

```
s_num s_den        # container is the closed square [0, s_num/s_den]^2
D                  # coordinate denominator: every point is (X/D, Y/D)
W                  # weight denominator: every weight is w/W
m                  # number of points
X_1 Y_1 w_1        # m lines: integer coordinates and integer weight numerator
...
X_m Y_m w_m
```

## What the file asserts

> For every closed unit square `Q` contained in `[0, s]^2` — at every centre and every angle —
> the total weight of the certificate points lying in `Q` is at least 1.

Together with `sum_i w_i / W < n`, this proves `s(n) >= s`:  rescale any packing of `n` unit
squares in a container of side `s' < s` by `s/s' > 1`, so each square strictly contains a
concentric closed unit square; those unit squares lie in pairwise disjoint interiors, so no
point is counted twice, giving `n <= sum w`.  (Formalised in `lean/Sqpack/Basic.lean`.)

Note the **closed** square convention: a point exactly on the boundary of `Q` **counts**.
This is what the rescaling argument above needs, and it is what `verify/` implements.

## The two certificates here

| file | s | points | total weight | proves |
|---|---|---|---|---|
| `s12_lower_3.931795.txt` | 3920/997 | 788 | 14916233/1250000 = 11.9329864 | s(12) >= 3.931795386 |
| `s12_56points_3.8.txt` | 19/5 | 56 | 56/5 = 11.2 | s(12) >= 3.8 |

The second is included because it is uniform — every weight is `1/5` — so it reads as a
purely combinatorial statement: *every closed unit square inside `[0, 3.8]^2` contains at
least 5 of these 56 points; twelve disjoint squares would need 60.*

## Scaling

A certificate can be scaled linearly (points and container together) by `λ`: multiply every
coordinate by `λ` and the container likewise — in this format, simply divide `D` by `λ`.
Each certificate stays valid up to a critical `λ`, beyond which some square's captured weight
drops below 1:

* the 56-point set is critical at `λ = 400/398 = 200/199`, container `760/199 = 3.819095…`;
  one step further, at `λ = 400/397` (container `1520/397 = 3.828715…`), an axis-aligned
  square drops to **4** points out of 5 — captured weight `4/5`, a clean rejection rather
  than a marginal one, stable from `N = 2000` to `N = 200000`;
* the 788-point set is critical just above container `3920/997`, with binding placements near 41°.

Both figures are reproduced by `search/scale_to_critical.py`, which binary-searches `D`
with the integer coordinates held fixed.
