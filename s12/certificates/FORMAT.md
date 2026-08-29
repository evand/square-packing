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

## Branch certificates (`region … / lambda … / k …` trailer)

The pure covering argument cannot reach `s = 4`: for every `s >= 3.99` there is a fractional
packing of mass `> 12` (`search/DUAL_EXACT.md`), so no weighted point set of total weight `< 12`
covers every unit square there.  A **branch certificate** carries one more piece of information
about the packings it refutes.  After the `m` point lines the file may continue with

```
region corner r_num r_den   # the four corner boxes [0, r]^2, [s-r, s] x [0, r], ... (closed)
lambda L                    # a multiplier, an integer numerator over W; may be negative
k K                         # number of squares of the packing whose CENTRE lies in a corner box
```

and then asserts:

> every closed unit square inside `[0, s]^2` whose centre lies in a corner box captures weight
> `>= 1 + L/W`; every other closed unit square inside `[0, s]^2` captures weight `>= 1`; and
> `sum_i w_i / W - (L/W)·K < n`.

**What it proves.**  Take a packing of `n` unit squares in a container of side `s' < s`, scale
it up to `[0, s]^2` and pass to the concentric closed unit squares as before; let `K'` be the
number of them centred in a corner box.  Each captures `>= 1`, those in the region `>= 1 + λ`
(`λ = L/W`), and no point is counted twice, so `n + λ K' <= sum_i w(S_i) <= W`.  Hence a branch
certificate with `W - λK < n` shows that **no packing has exactly `K` squares centred in the
corner boxes**.  Because the centres of two interior-disjoint unit squares are at least `1`
apart and the admissible part `[1/2, r]^2` of a corner box has diameter `(r - 1/2)·sqrt 2`, the
verifier insists on `(2r - 1)^2 < 2` exactly, so each box holds at most one centre and
`K ∈ {0, 1, 2, 3, 4}`: **five branch certificates, one per `K`, together prove `s(n) >= s`**.
Nothing else changes — the same scaling argument, the same closed-square convention.  With
`L = 0` a branch certificate is a plain one.

**Per-box multipliers.**  The trailer may instead carry four multipliers and an occupancy
pattern,

```
lambda L1 L2 L3 L4          # box 1 = [0,r]^2, 2 = [s-r,s]x[0,r], 3 = [0,r]x[s-r,s], 4 = [s-r,s]^2
k K1 K2 K3 K4               # K_j = 1 iff a square of the packing is centred in box j
```

asserting that a square centred in box `j` captures `>= 1 + L_j/W`, every other square `>= 1`,
and `sum_i w_i/W − sum_j (L_j/W)·K_j < n`; it then refutes every packing with exactly that
occupancy pattern, and the sixteen patterns (six up to the symmetries of the square) cover all
packings.  The one-number form is the special case `L_j = L` with `K` the number of occupied
boxes.  Unequal `L_j` make the claim non-symmetric, so the verifier then sweeps the full range
`[0°, 90°)` whatever the symmetry of the point set.  (Lean: `packing_le_weight_regions`.)

**How it is verified.**  The verifier's arrangement sweep works cell by cell in the centre
plane; a cell that may meet a corner box is required to reach `1 + λ`, a cell that may leave the
boxes is required to reach `1`, and a cell that straddles the boundary is required to reach both
(the membership tests are done in floating point with a padding that can only make them
stricter, on the cell's bounding box, which is exact for the axis-parallel boxes).  Witnesses
carry a fifth column, the flag of the threshold they violate, for the LP.  The trailer keywords
are mandatory and any other trailing data is still an error.  The corner boxes are symmetric
under the symmetries of the container, so the `[0°, 45°]` reduction for D4-symmetric point sets
remains valid.  `search/branch.py` produces these certificates (`search/BRANCH.md`).
