# Strip split for `s(12) = 4`: guillotine cases, the (0,8) leaf, and rectangle capacities  (2026-10-11)

Evan's sketch (10-11): cut the packing into rectangle packings.  A vertical strip of 4 squares in a
`(1+a) x (4-b)` rectangle, `k` squares to its left and `8-k` to its right, and a seam allowance `c` (how far
a square from the right can intrude into the strip's box).  Then reduce each case to statements about
packing squares in rectangles.  This note records what is proved, what is measured, and what the (0,8)
case would need.

Labels as in `notes/n12-gap.md`: **[proved]**, **[measured]** (float, one-sided), **[heuristic]**, **[guess]**.
Semantics: side `L = 4 - b < 4` with disjoint interiors, or `L = 4` with closed disjointness.  The
strip lemma covers both (`notes/chord-lemma.md` §0).

---

## 0. Verdict

1. **The guillotine cases are dead with existing tools** (§2) [proved, modulo Lean].  Suppose two
   vertical lines meet no square interior and split 12 squares into pieces of `k`, 4 and `8-k`.  Then
   the strip lemma at any `h >= 7/6` kills every `k`.  Exact certificates exist at `h = 1.17` and `1.19`
   (`tasks/n12-chat-lemmas/strip-demo/`); Lean has only `h = (3√2-2)/2 = 1.1213`, which is not enough.
2. **`a` is not small.**  For any `b > 0`, 4 squares in a `w x (4-b)` rectangle need `w > h + 1/2`.
   That is `> 1.69` proved, about `1.9` in fact (Stromquist: 4 fit in `1.9 x 3.9475`) [proved / pub].  The
   width-1 column exists only at `b = 0`, and closed semantics kills it there too.
3. **A scalar seam allowance `c` cannot close the (0,8) case** (§4) [measured].  The strip lemma
   accounting needs `c <= 0.30` (`h = 1.19`), or `0.33` at the sharp `h*`.  In the minimal-width strips
   one unit square intrudes to `c ≈ 0.87–0.89 ≈ a`, both beside the flush-left end square and in the
   middle notch.
4. **So (0,8) is a leaf only with structure, not with a number** (§5).  The leaf has slack only for
   strips with substantially tilted squares (Stromquist type).  Near-axis strips of every width in
   `(1.69, 2)` include tiling-minus-4 configurations at margin 0 (§4, correction).  It needs
   (i) a shape theorem for 4 squares in `(1+a) x (4-b)`, for `a` in a narrow range below about 1;
   (ii) a width bound for 8 squares in a rectangle of height `4-b`;
   (iii) a statement that no near-optimal 8-packing can use the strip's notches.
   It also needs the case split itself to be exhaustive, which is the real gap (§3).
5. **Exhaustive version: sweep a line** (§3).  For `x0` in `[4 - b - 2h, h + 1/2]`, at least 3 squares
   cross the line `x = x0` [proved, same lemma].  This is tight at the permutation-hole tilings.

## 1. Rectangle capacities from the strip lemma

The strip lemma (B1/B2, `notes/chord-lemma.md`; sharp threshold `h* = (1+√2)/2`, counterexample at
1.2081): at most 3 unit squares in `[0,L]^2`, `L <= 4`, have centre within `h` of a given wall.  It applies
to squares in a `w x (4-b)` rectangle, `w <= 4-b`, against either long wall.  To use the far wall,
translate the rectangle so that wall lies on the container's wall.

* **R4.**  4 squares in `w x (4-b)` need `w > h + 1/2`.  Centres satisfy `x <= w - 1/2`, so if
  `w - 1/2 <= h` all of them are in the left band.
* **R7.**  7 or more squares need `w > 2h`.  If `w <= 2h` the left and right bands cover every centre,
  so there are at most 3 + 3.

| `h` | source | R4: `w >` | R7: `w >` |
|---|---|---|---|
| 1.1213 | Lean `wall_strip_le_three_line` | 1.621 | 2.243 |
| 7/6 | (threshold for §2) | 1.667 | 2.333 |
| 1.19 | exact cert, 191 boxes | 1.690 | 2.380 |
| 1.2071 | `h*`, sharp; chat cert, 8.2M leaves | 1.707 | 2.414 |

Cross-check: Nagamochi 2005 Thm 1 gives `ν(3, 4) < 10`, and for 8 squares in `w x (4-)` it gives `w > 7/3`,
the same as R7 at `h = 7/6`.  But that theorem is not formally verified.  chelokot's Lean work shows one
of its per-square scoring lemmas is false (Karakus 2026) and proves only the square-container corollaries,
so do not cite the rectangle form as settled.

## 2. Guillotine splits

Pieces `(k, 4, 8-k)` side by side, widths summing to `4-b < 4`.  With `h >= 7/6`:

| split | width bounds | sum |
|---|---|---|
| (0,8) | `> 5/3`, `> 7/3` | `> 4` |
| (1,7) | `>= 1`, `> 5/3`, `> 7/3` | `> 5` |
| (2,6), (3,5) | `>= 1`, `> 5/3`, `> 5/3` | `> 13/3` |
| (4,4) | `> 5/3` each | `> 5` |

All die.  This is a clean intermediate theorem: "no guillotine packing (two parallel full cuts, a
4-piece in between) of 12 squares in side `< 4`", needing only the `h = 1.17` certificate in Lean.
It does **not** cover a single cut `6 | 6`.  That needs "6 squares do not fit in `(2 - b/2) x (4 - b)`",
which the strip lemma permits (R7 starts at 7) [open].

## 3. Exhaustiveness: where the case split has to come from

The split assumes 4 squares fill a full-height strip, with every other square left or right of it.
The configurations still alive at margin 0 (the permutation-hole tilings: every row and column has 3)
have no such strip.  A width-1 column of 4 is killed by B1 at a wall and by the chord-line lemma in the
interior.  A chain of 4 therefore uses two adjacent columns (`a ≈ 1`), and the far side shares a
column with it (`c ≈ 1`).  The census witnesses (`tasks/n12-chat-lemmas/census/`) all have wall-to-wall
axis-parallel paths, but not 4-squares-in-a-strip.  So (0,8) as posed is a **leaf** that some other
branch has to route into.

**Line sweep (exhaustive).**  For any `x0`, let `S` = squares meeting `x = x0`, `L` = wholly left,
`R` = wholly right.  `L` and `R` are clean rectangle packings.  If `x0 <= h + 1/2` then `|L| <= 3` (R4).
If `4 - b - x0 <= 2h` then `|R| <= 6` (R7).  So for `x0` in `[4 - b - 2h, h + 1/2]` (about
`[1.62, 1.69]` at `h = 1.19`), `|S| >= 3`, and symmetrically for the mirror line and the horizontal
lines [proved].  It is tight (3 + 3 + 6) at the permutation-hole tilings.  The strip is now "the
squares crossing a line", and intrusion is handled exactly.  The open target becomes a joint statement
about 3 or more crossing squares with at most 3 on one side and at most 6 on the other.

## 4. Measurements for the (0,8) leaf (float, scratch scripts, 10-11)

Penalty multistart, separating-axis overlap; widths are **upper bounds**, and the search is loose.  It
does not reproduce Stenlund's `1.9303 x 3.8606` or Stromquist's `1.9 x 3.9475`, both of which beat it.

* **4 squares, minimal width at height `H`**: found `1.952` (3.90), `1.939` (3.95), `1.914` (3.99, 3.999).
  True values are lower: `<= 1.9303` for `H >= 3.8606` (Stenlund) and `<= 1.9` for `H >= 3.9475` (Stromquist).
  **Shape** (every run, and both literature packings): two axis-parallel squares at the ends, in
  **opposite** corners (one flush left, one flush right), and two tilted squares in the middle
  (24°–66° in our runs, 39.6° Stromquist, 40.5° Stenlund).
* **5 squares**: none found below `w ≈ 2.0` at `H = 3.95, 3.99`.  So no fifth square fits inside the
  strip box.  This bounds `c < 1 + a`, but nothing better.
* **Intrusion into the minimal strips** (one extra unit square, disjoint from the 4, in `[0,H]^2`,
  minimising its leftmost x): `c = 0.87–0.89` at `H = 3.99, 3.999`, `0.75–0.81` at 3.95, `0.67–0.69` at 3.9.
  The deepest spots are the slot beside the flush-left end square (height about 0.77, so only a
  tilted square's corner gets in) and the middle notch.  Only the end square that is flush right blocks
  intrusion over its unit of height.
* **Budget.**  The 8 squares need width `> 2h` (R7).  With the strip at `1 + a`, survival in the scalar
  model needs `c > 2h - 3 + a + b`, which is `0.30` at `h = 1.19` and `a = 0.914`.  Measured `c` is about
  `a`.  **No (a, b) region dies with a scalar `c`.**
* **8 squares, minimal width at height `H`**: found `2.863` (3.99), `2.942` (3.95), `3.00` (3.90; the
  search is loose there).  **Shape at 3.95 and 3.99**: an axis-parallel square in every corner and a
  tilted block of 4 in the middle (all at 57.3° at `H = 3.99`).  The top and bottom rows each hold two
  axis squares with a gap of about 0.86 between them.
* **The composition has real slack.**  `w4 + w8 ≈ 1.9 + 2.86 = 4.76` against `4 - b`.  Even the
  optimistic seam needs about `0.77` of usable intrusion (less once the true minima, below our float
  values, are known).  Single-square intrusion reaches `0.87`, but only at particular heights, and only
  for a tilted square.  The 8-packing's corner squares cannot use it.  This slack holds only for the
  **Stromquist-type strip** (two substantially tilted middle squares).  See the correction below.
* **Correction (same day): the leaf "4 squares wholly in a wall strip of width `w < 2`" is at margin 0.**
  A 12-square search with that constraint (`leaf.py`) reaches `L = 4.000` at `w = 1.95, 1.99` and
  `L = 4.0016` at `w = 1.9`.  Every hit is a tiling minus 4: three tiles of column 0 plus a column-1
  tile slid left into column 0's hole (strip width `2 - s`), or two tiles from each column staggered
  into each other's holes.  For every slide `s < 0.31` (strip width above `h + 1/2 ≈ 1.69`) this is a
  touching configuration at `L = 4`.  So the leaf contains the margin-0 core for every strip width in
  `(1.69, 2)`.
  Below height 4 there is also a **staggered strip class** (`stagger.py`): three flat squares
  alternating sides and one square tilted 7° (`w = 1.99`) to 19° (`w = 1.95`), at `H = 3.99, 3.999`.
  The tilt goes to 0 as `w → 2`.  Its slots beside the flush-left squares are almost a unit high, so
  the other 8 interlock with it as in the tiling (`c ≈ a`).
  A third class, the **wall column**, shows up as the strip of a margin-0 hit (`leaf_w1.99_s2`): three flat
  squares stacked at the wall and one square tilted 24–33° beside them, dipping into the gap.  This is the
  strip lemma's extremal shape.  It survives below height 4 (`third.py`: fits at `w = 1.95–1.99`,
  `H = 3.95, 3.99`).  No end square is flush right, so the top-unit block of §4b does not apply: the
  notched border degenerates to `n ≈ 0`, and 8 squares fit in the `2.99 x 3.99` rectangle.  With the
  actual strip squares fixed, 8 more were not found to fit (`fixed4.py`, 200 starts).  But the best
  attempt overlaps by only `2.5e-3` at `H = 3.99`, so this class behaves like the near-margin-0 core.
  **The leaf with slack is therefore defined by shape, not by width:** a strip containing squares
  tilted by at least some `θ0`.  The near-axis strips belong to the core (obstruction X,
  `notes/n12-gap.md` §4.11).  If it holds, the tilted leaf removes part of the far field, which
  `n12-gap.md` §6 item 4 lists as having no candidate mechanism.

### 4a. Census of 4-square strips (local minima of width, `census4*.py`, 10-11) [measured]

448 local minima from random starts at `H = 3.99` (324) and `3.95` (124).  About 85% land on the
`w = 2` plateau: two flat squares side by side, with the rest of the height loose.  That plateau is
one near-tiling case, not a family of classes.  Below `w = 1.999` at `H = 3.99`, signatures run
bottom to top (F flat, T tilted; letters = walls touched), widths at 3.99 / 3.95:

| class | shape | `w` at 3.99 | `w` at 3.95 | tilts |
|---|---|---|---|---|
| A tilted middle (Stromquist) | flat bottom-left, two tilted, flat top-right | **1.8851** | 1.8991 | 39.3°, 39.3° |
| C wall column (3 variants: flats on one side or alternating) | three flat squares, one 45° diamond in a gap | 1.9192 | 1.9394 | 45° |
| B three tilted | flat bottom-left corner, three co-rotating tilted squares, top one touching right and top | 1.9209 | 1.9350 | 34.9° ×3 |
| F all tilted | four co-rotating tilted squares, zigzag | 1.9494 | — | 27.5° ×4 |
| D two flat at one end (2 variants) | two flat squares at the bottom, two tilted above | 1.9708 | 1.9877 | 31° |
| E both ends on one side | flat bottom-left, two tilted, flat top-left | 1.9841 | 1.9958 | 19° |
| near-plateau stragglers | three tilted at 14°, etc. | 1.990–1.993 | — | small |

Wall column matches the closed form `1 + √2 − (1−b)/2` (1.9192 at `b = 0.01`).  The stitch's
top-unit block applies to A (flat flush-right top).  B and F have a tilted top square touching the
right wall and the top: a partial block, still to be tested.  C, D and E have no flush-right end.

### 4b. The notched seam (Evan's shape, 10-11) [measured]

The border is the full strip width `1+a` over the top unit (blocked by the flat, flush-right top
square), and `1+a-c` below.  Take the most optimistic lower seam, `c = a`: below the top unit the 8
squares may come in to `x = 1`.  The leaf then dies iff 8 squares do not fit in
`[0, H-1] x [0, H]` minus a `n x 1` top-left notch, where the notch width `n` is `a`.  Seeded search
(`notch_seeded.py`, 300 starts from the 2.86 packing plus perturbations; control without a notch fits at 2.87):

| `H` | width | notch `n` | result |
|---|---|---|---|
| 3.99 | 2.99 | 0.3, 0.5 | fits: the top-left corner square moves right of the notch |
| 3.99 | 2.99 | 0.6, 0.69, 0.8, 0.9 | no fit found |
| 3.95 | 2.95 | 0.9 | no fit found |

9 squares in `2.99 x 3.99` without a notch: none found below width 3.00.  So the notch threshold is
between 0.5 and 0.6.  The proved `a > 0.69` already clears it, and the true `a ≈ 0.9` clears it by
about 0.3.  That 0.3 is the margin available for a simplified border: the top square tilted or
offset, the notch shorter than 1, and so on.  Kept in the leaf hypotheses: the top-unit block (shape
theorem), and "nothing reaches `x < 1` below the top unit" (intrusion lemma, measured `c <= a`).

### 4c. Seam with the real top square, classes A, B, F (`seam.py`, 10-11) [measured]

At `H = 3.99`, 8 squares must avoid the column `[0,1] x [0,H]` (the optimistic lower seam) and the
strip's actual top square, oriented to touch the right wall and the top.  Result, 300 starts each:
**no fit** for A (`w = 1.8851`, flat top square), B (`1.9209`, top square at 55°) and F (`1.9494`, at
62.5°).  Positive control: the same search without the top square fits 8 on the second start.  So
B and F belong with A, as far as this search can tell.  The leaf split is {A, B, F} against
{C, D, E, `w = 2` plateau}.

### 4d. Margins: what the seam and the region each tolerate (10-11) [measured]

* **Region side, top block** (`seamshift.py`): moving the strip's top square back into the strip by
  `δ`, 8 squares start to fit between `δ = 0.30` and `0.45`, for each of A, B, F.  Slack about 0.3.
* **Region side, seam position** (`seamshift2.py`): with the column lowered to `x0 = 0.9` or `0.8`,
  8 fit at once.  Width `H − x0 ≥ 3` admits three flat columns.  So the region needs
  `x0 > H − 3 = 1 − b`: **margin `b` in a constant seam**.  The binding structure is the tiling
  (the strip's flat bottom square in column 0, a 3 x 3 grid in columns 1–3).
* **Seam side** (`seamjoint.py`, `seamband*.py`; strip deformed within the A/B/F contact pattern,
  intruder chosen jointly):
  - tightest strips: reach 1.10–1.20 (margin 0.10–0.20);
  - strips with slack (`w = 1.95, 1.99`): reach exactly 1.000, flush beside the flat bottom
    square (axis contact);
  - middle rows only (intruder centre in `[1.3, 2.7]`): reach 1.038–1.042 at `H = 3.99`, decaying
    to 1.007 at `H = 3.999` and `3.9999`.  The worst strips have middle tilts that fall with `b`
    (25° → 11°): with slack, the classes drift toward near-axis strips;
  - with the middle tilts held at 30° or more: 1.053 (`H = 3.99`), 1.038 (`H = 3.999`).
* **Reading.**  With slack, A/B/F contain near-axis strips, and through them the margin-0 core comes
  back.  The leaf needs (i) class regions with a tilt floor on the middle squares (`θ ≥ θ0`; the
  low-tilt strips go to the core branch), and (ii) a profiled seam `x0(y)`: 1 in the bottom unit,
  `1 + m` in the middle rows, the top block above.  Even with the floor the seam margin is
  about 0.04 and not obviously bounded away from 0 as `b → 0`.  Open: whether `m(θ0)` plateaus at
  `H → 4`, and whether the region lemma holds with a profiled seam of that size.

## 5. What a (0,8) leaf proof would need

Evan's reading (10-11): the near-optimal 8-packings at heights near 4 put an axis-parallel square in
every corner.  A corner square cannot use a notch shorter than 1, and the strip's notches are shorter
than 1.  So the 8 squares effectively get only `4 - b - (1 + a)` at the top and bottom, and the
composition fails.  Making that a proof needs:

1. **Strip shape theorem.**  For `a` in `[a_min(b), a_1)`: every 4-packing in `(1+a) x (4-b)` is the
   end-flat/middle-tilted shape, with a quantified right-hand profile (where a unit square can reach, as
   a function of `y`).  At `a ≥ 1` the strip has no structure (a `2 x 2` block leaves half the box
   empty), so the range is narrow, about `[0.9, 1)`.
2. **8-square width bound** at height `4 - b`: a number, provable by a parameterised point cover
   (wand125's domino method, one aspect ratio per certificate, or a family in `H`).
3. **No good-enough alternative 8-packing**: no 8-packing of the right-hand region, whose left boundary
   is the strip profile from 1, exists.  This is the joint part.  It replaces "8 in a rectangle" by "8 in
   a rectangle with a profile on one side", which is the boundary-condition statement `search/RANK8.md`
   §2.3 said a split would need.
4. **A branch that routes configurations into this leaf** (§3).

## 6. Rectangle results on hand

| statement | where | form |
|---|---|---|
| ≤ 3 centres within `h` of a wall, side ≤ 4 | `lean/Sqpack/Chord.lean` (h = 1), `ChordLine.lean` (1.1213) | Lean |
| same, `h` = 1.13 … 1.19 | `tasks/n12-chat-lemmas/strip-demo/` | exact rational cert, not Lean |
| same, `h` = 603553/500000 | group chat (evidence tarball not here) | cert, 8.2M leaves |
| 4 squares in a 1:2 rectangle: `s = 1.9303` (root of `5x³+8x²−32x−4`) | wand125 `problems/domino/` | exact, two independent checkers, 7 weighted points, **no Lean** |
| `ν(a,b) < ab − (a+1−⌈a⌉) − (b+1−⌈b⌉)`, `a,b ≥ 2` | Nagamochi 2005 Thm 1 | paper; a lemma in it is false; rectangle form unverified |
| 4 squares in `1.9 x 3.9475` | Stromquist 1984-III p.10 (`notes/proof-anatomy.md` §5.1) | construction |

Nothing anywhere gives the minimal width of 4 squares at height near 4, or of 8, as a theorem.  The
domino certificate is the closest template.  Its "wedge lemma" handles the same non-rigid tilted pair
the strip has, so its point-cover method should transfer to `w x H` with `H` as a parameter.
