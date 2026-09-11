# Certificate format

A certificate is a plain-text file of whitespace-separated integers. Everything is exact:
no floating point appears anywhere in a certificate or in its verification.
(Each certificate is also published as a self-describing `.json` companion; see
[JSON companion](#json-companion-pointsjson) below.)

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

**Well-formedness.**  A file is valid only if all header values are positive integers,
`s_den` divides `s_num · D` (so the container side is a multiple of the coordinate unit),
every point lies in `[0, s]^2`, every weight numerator is `>= 0`, and there are exactly `m`
point lines.  The verifier refuses (exit 2, `ERROR:`) anything else rather than reporting a
verdict.  Non-negativity is not a formality: the reduction `n <= sum w` (and its Lean proof,
hypothesis `hw`) needs it, and a negative-weight point outside the container would otherwise
lower the total without touching the covering.  Duplicate coordinates are allowed and their
weights add — the file describes a weighted multiset.

Note the **closed** square convention: a point exactly on the boundary of `Q` **counts**.
This is what the rescaling argument above needs, and it is what `verify/` implements.

## The certificates here

| file | s | points | total weight | proves |
|---|---|---|---|---|
| `s12_lower_3.9686.txt` | 15680/3951 | 1736 | 11.9738036 | s(12) >= 3.968616 |
| `s12_lower_3.9676.txt` | 980/247 | 764 | 11.9962288 | s(12) >= 3.967611 |
| `s12_lower_3.931795_sparse.txt` | 3920/997 | 224 | 11.9834372 | s(12) >= 3.931795 |
| `s12_lower_3.931795.txt` | 3920/997 | 788 | 14916233/1250000 = 11.9329864 | s(12) >= 3.931795 |
| `s12_uniform_7of81_3.888.txt` | 35/9 | 81 | 81/7 | s(12) >= 3.888889 (uniform: every square contains 7 of 81) |
| `s12_uniform_<k>of<m>_<s>.txt` | see `search/uniform/UNIFORM.md` | m | m/k | uniform family, k = 1..8 |
| `s12_56points_3.8.txt` | 19/5 | 56 | 56/5 = 11.2 | s(12) >= 3.8 |

The uniform ones are included because every weight is equal, so they read as purely
combinatorial statements — e.g. the 56-point one: *every closed unit square inside `[0, 3.8]^2` contains at
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

## JSON companion (`points.json`)

Every `.txt` certificate ships next to a `.json` file with the same stem
(`s12_lower_3.931795.json`, `s12_56points_3.8.json`), so that the data is self-describing
and readable without this document.  The schema follows the `points.json` of Mira's
`17squares` repository (`problem`, `theorem`, `container_side` as
`{numerator, denominator, decimal}`, `point_denominator`, `points` as `[X, Y]` integer
pairs, provenance with a file name and SHA-256), extended with the fields the weighted
method needs.  The `.txt` remains the source of truth: it is what `verify/` reads and what
`SHA256SUMS` pins first; the JSON is generated from it by `search/export_points.py` and is
pinned too.

| field | type | meaning |
|---|---|---|
| `schema` | string | `"s12-points-json/1"` |
| `problem` | string | `"packing 12 unit squares in a square"` |
| `theorem` | string | `"s(12) >= 3920/997"` (note `>=`, not `>`: closed squares, see `convention`) |
| `bound` | string | the container side as an exact fraction, `"3920/997"` |
| `n` | integer | number of squares the certificate rules out: 12 |
| `container_side` | object | `{"numerator": s_num, "denominator": s_den, "decimal": …}` — the container is `[0, s_num/s_den]^2` |
| `point_denominator` | integer | `D`: point `i` is at `(points[i][0]/D, points[i][1]/D)` |
| `weight_denominator` | integer | `W`: point `i` has weight `weights[i]/W` |
| `points` | array of `[X, Y]` | integer coordinates, in file order |
| `weights` | array of integers | integer weight numerators, parallel to `points` (same length, same order) |
| `total_weight` | object | `{"numerator", "denominator", "decimal"}` — `sum(weights)/W` in lowest terms; the certificate is valid only if this is `< n` |
| `convention` | object | `squares: "closed"`, `boundary_points_count: true`, and prose `claim` / `argument` / `note` fields restating the closed-square, concentric-shrink argument of this document, with `reference` pointing here |
| `source` | object | `file`: the `.txt` file name; `sha256`: its SHA-256 (the value in `SHA256SUMS`); `format`: this document |

Every number that matters is an integer or an exact rational; the `decimal` values are the
nearest IEEE double, provided for readers and ignored on import.  Keys appear in the order
above, one point and one weight per line, so the file is deterministic and diffs cleanly.

Both directions are checked by `verify.sh`:

```sh
python3 search/export_points.py --roundtrip certificates/s12_lower_3.931795.txt certificates/s12_lower_3.931795.json
```

asserts that rebuilding the `.txt` from the JSON's integer fields reproduces the shipped
`.txt` byte for byte, that regenerating the JSON from the `.txt` reproduces the shipped
`.json` byte for byte, that `source.sha256` is the `.txt`'s hash, and that `total_weight`
equals `sum(weights)/W` and is below `n`.  `export` (txt → json) and `import`
(json → txt) are the two conversions; they need nothing beyond the Python standard library.

What is **not** carried over from Mira's schema: `triangles` / `triangle_indices_are_zero_based`
(their strict triangle-piercing lemma has no counterpart here — the weighted argument needs no
combinatorial structure on the points) and the subdivision-tree statistics under `certificate`
(there is no tree; the covering property is checked over the whole continuum by the
arrangement sweep in `verify/`).  Their points are *unweighted* and their squares *open*;
our `weights` and `convention.squares = "closed"` are exactly the two places where a reader
of both formats must not assume the same semantics.
