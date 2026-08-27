# Anatomy of the s(m²−3) = m proofs, and what they say about n = 12

Literature report, 2026-08-26. Sources read in full (PDF → text; page numbers are the
journal/arXiv page numbers):

| # | Paper | Access |
|---|-------|--------|
| B10 | W. Bentz, *Optimal packings of 13 and 46 unit squares in a square*, EJC 17 (2010) #R126, 11 pp. | full text |
| N05 | H. Nagamochi, *Packing unit squares in a rectangle*, EJC 12 (2005) #R37, 13 pp. | full text |
| B16 | W. Bentz, *Optimal packings of 22 and 33 unit squares in a square*, arXiv:1606.03746v1, 9 pp. (I found no journal version; the arXiv v1 is what is analysed.) | full text |
| S03 | W. Stromquist, *Packing 10 or 11 unit squares in a square*, EJC 10 (2003) #R8 | full text |
| S84-I/II/III | W. Stromquist, *Packing unit squares inside squares I (six), II (ten), III (n ≤ 65, Gardner's conjecture)*, DHWA internal memoranda 1984, scanned PDFs from walterstromquist.com/papers/squares{1,2,3}.pdf | read page-by-page as images (no OCR available); I,II fully, III fully |
| KS02 | M. Kearney, P. Shiu, *Efficient packing of unit squares in a square*, EJC 9 (2002) #R14 | full text |
| DS7 | E. Friedman, *Packing unit squares in squares: a survey and new results*, EJC Dynamic Survey DS7 (github copy; the copy is the 2009 version, figures not available as text) | full text (HTML) |

Not accessed: Göbel 1979 (Math. Centrum Tracts 106) and El Moumni 1999 (Studia Sci. Math. Hungar.);
their content is only known here through DS7/N05/B10 citations.

Notation. s(n) = side of the smallest square holding n unit squares. Throughout, a **box** (Stromquist,
Bentz) is the *interior* of a square of side > 1 (B10: any side > 1; S84/S03: side 1+ε, 0<ε<10⁻⁴;
B16: side s with 1 < s ≤ 1.01). "n boxes cannot be packed in the closed square [0,m]²" ⇔ s(n) ≥ m.
Boxes are open; points of an unavoidable set must lie *inside* the box. This is the single most
important convention for the s(12) question (see §7).

---

## 1. Common framework: how "side ≥ m" (not "> m−ε") is obtained

All three papers prove the exact bound by the same rescaling trick, stated most explicitly in KS02 §2
and DS7 §5:

> DS7 §5: "Shrinking these by a factor of (1−ε/k) gives a set P′ of (n−1) points in a square S′ of
> side (k−ε) so that any unit square in S′ contains an element in P′ in its interior. Therefore no
> more than (n−1) non-overlapping squares can be packed into a square of side (k−ε), and s(n)>k−ε.
> Since this is true for all ε>0, we must have s(n)≥k."

Equivalently (S84-I p.2–3): "rather than considering a reduction in the size of the bounding squares,
we instead consider an increase in the size of the squares being packed". So every certificate in
this literature is a statement about **dilated** squares (side 1+ε, every ε>0) in the **closed**
container, i.e. in the limit about closed unit squares with *closed* containment (a point on the
boundary of the unit square counts as covered). The 16 unit squares of the trivial 4×4 tiling
each contain the grid points {1,2,3}² only on their boundaries, and those points are shared, so
the "16 eroded squares fit ⇒ cost ≥ 16" obstruction does not exist in this semantics. (See §7.)

The strictness is used in exactly three ways:
1. "≥" from a lemma becomes ">" for a box (S84-I p.3: "whenever a calculation seems to produce an
   inequality which is not strict, we can appeal to the ubiquitous ε").
2. Nagamochi scales the *weights* with λ (area ×λ², length ×λ) and proves σ(S) > 1 for all
   λ ∈ (1, 1.01] (N05 p.4–5).
3. Bentz 2016 needs the box size bounded above (≤ 1.01) so that "two points in one box are within
   1.01√2" (B16 p.7) and "the midpoint of B must lie within 0.505√2 of P" (B16 Lemma 7).

None of the three papers is computer-assisted in any declared sense. Numerical content is limited
to: the minimisation of Stromquist's f(a) (cubic (2) in cosθ; B10 footnote 1 corrects S03: the
minimum can be the *second largest* positive root), and in B16 the coordinates of the shaded
regions of Figs 4–5 ("approximate coordinates (1.12, 2.5±0.05), (1.2, 2.5±0.1)"; "(1.13,0.56) and
(1.13,1.24)") with "an easy calculation shows that the entire area is within a distance of 0.5 from
the point (1.5, 2.5)". No code, no interval arithmetic, no verified figures.

---

## 2. Bentz 2010 (n = 13 and 46)

### 2.1 Lemma toolbox (B10 §1, pp. 2–3), stated exactly

Box = interior of any square of side > 1. "Non-avoidance lemmas state that if the centre of a box is
in some region, the box must intersect some part of the region's boundary."

- **Lemma 1** (Nagamochi, Stromquist). "A box whose centre is in the rectangle [0,1]×[0,1] and which
  does not intersect the axes must contain the points (1,1), (0.9,1), and (1,0.9) (and hence the
  triangle spanned by these points)." [= N05 Lemma 7(ii); the (1,1) part is DS7 Lemma 1 = S84-I Cor. 2.]
- **Lemma 2** (Friedman, Stromquist). "Let T be a triangle with sides of length at most 1. Then any
  box whose centre is in T must contain one of the vertices of T." [DS7 Lemma 3; S84-I Lemma 2.]
- **Corollary 3.** "Let 0 < a, b ≤ 1 and c such that both c²+b² ≤ 1 and (a−c)²+b² ≤ 1. Then if a box
  has its centre inside the rectangle R with vertices H(0,0), I(0,b), J(a,b), K(a,0), it covers a
  vertex or the point L(c,0). In particular, if b < ½√3, the box covers a vertex or the midpoint of
  the line segment HK." (Proof: triangles HLI, ILJ, JLK + Lemma 2.) Used with b = 0.73 and 0.54 on
  1×b rectangles: the box covers a vertex or a midpoint of a long edge.
- **Lemma 4** (Friedman, Stromquist). "Let a ≤ 1, b ≤ 1, and a + 2b ≤ 2√2, then any box whose centre
  is in the rectangle [0,a]×[0,b] must intersect the x-axis, the point (0,a) or the point (a,b)."
  [Typo: "(0,a)" should read "(0,b)"; cf. S03 Lemma 3, S84-I Lemma 4, DS7 Lemma 2.] "We will
  usually use Lemma 4 in the cases of a < 2√2−2 ≈ 0.828, b = 1 and a = 1, b < √2−½ ≈ 0.914."
  When the x-axis is a container wall the box cannot touch it, and when (0,b) is *on* a wall it
  cannot be covered either, so the lemma becomes "covers (a,b)".
- **Lemma 5** (Stromquist). "Let 2√2−2 < a < 1, 0 < b < 1, and (a,b) within a distance of 1 from
  (0,1). Moreover, let f(a) be the infimum of cosθ/(1+cosθ) + (1−a cosθ)/sinθ (1) for θ ∈ (0,π/4].
  If b < f(a), then any box whose centre is in the quadrilateral with vertices (0,0), (0,1), (a,0),
  and (a,b) must intersect the x-axis, the point (0,1) or the point (a,b). Moreover, the infimum of
  (1) is a minimum and is obtained at a value of θ satisfying
  2cos³θ − (2a+2)cos²θ + (a²−2a+3)cosθ − (1−a²) = 0 (2)." Used with (a,b) = (½√3, 0.5); a ≤ 0.89,
  0.6 < b ≤ 0.921; (0.90, 0.90); (0.96, 0.76). "f(a) is decreasing in a." Values (S03 p.3 table):
  f(0.853)=0.972, f(0.894)=0.926, f(0.96)=0.769.
- **Lemma 6** (Stromquist). "Let L₁ and L₂ be two parallel lines of distance d ≤ 1, and B a box with
  its centre between them. Then B must intersect the two lines with a common length of intersection
  of at least min{1, 2√2−2d}." "Common length" means the **sum** of the two intersection lengths
  (S84-I Lemma 3: "B intersects these lines in a segment or segments with total length greater than
  min(1, 2√2−2a)"; the 45° box centred midway between lines at distance 1 meets each in 2√2−2 over 2).
  Note the critical value: d = √2−½ ≈ 0.914 gives exactly 1. Rotation is fully general in all lemmas.
- **Corollary 7** ("a major technical tool"). As printed: "Let 0 < b ≤ 1, and R be the rectangle with
  vertices (0,0),(0,1),(b,0),(b,1). Then any box whose centre lies inside R without containing any of
  its vertices intersects the line segments {0}×[0,b] and {1}×[0,b] with a common length of
  intersection of at least 2√2−2 ≈ 0.828. In addition, if 2√2−2 > b then the box intersects each
  segment with a length of at least b − 2√2 + 2." The intended statement (consistent with every use
  in §3.2: 0.828−0.73 → "0.09", 0.828−0.37 → "0.45", 0.828−0.54 → "0.28", 0.828−0.32 → "0.50") is:
  R = [0,1]×[0,b]; a box centred in R and avoiding its four vertices meets the two vertical edges
  {0}×[0,b], {1}×[0,b] with total length ≥ 2√2−2, hence, since each intersection is confined to an
  open edge of length b, each edge is met in length ≥ **2√2 − 2 − b**. "Left/right/up/down
  intersection" = length of the box's intersection with that edge.

### 2.2 n = 46 (B10 §2, Theorem 8, p. 3–4): a pure unavoidable set

"Let S be the square [0,7]² and consider the collection of 45 points depicted in Figure 1. The points
in the lowest row are (i, √2−½), i = 1,…,6, and the remaining ones are arranged so that all shown
triangles are equilateral of side length 1." Seven rows at y = 0.914 + k·(√3/2), alternating 6 and 7
points (6+7+6+7+6+7+6 = 45); the top row has y > 6.11, i.e. within √2−½ of the top edge. Partition:
equilateral triangles (Lemma 2), wall rectangles (Lemma 4 with a=1, b=0.914 at the bottom/top,
a=0.828, b=1 style at the sides), and the two half-triangles at each side wall handled by Lemma 5
with (a,b) = (√3/2, ½). "Thus any box placed inside S must contain one of the points, showing that no
more than 45 boxes can be placed in S." No counting, no cases. ("Its proof is surprisingly simple for
such a large case, in particular when compared to the case of 13 squares.")

### 2.3 n = 13 (B10 §3, Theorem 9, pp. 4–10)

**Container**: S = [0,4]², 13 boxes assumed. **Initial set** (Fig. 2): 16 points,
A(1, 0.914), B(0.914, 1), C(0.914, 2), D(1.65, 1.65), and their images under x↦4−x, y↦4−y, x↔y.
(A,B: the corner points required by Lemma 4; C: mid-wall; D: interior. Rows: y ≈ 0.914, 1, 1.65, 2, …)
"Non-avoidance lemmas apply to all regions, so each box contains at least one point."

**Counting step 1** (p. 4): "With 13 boxes and 16 points, at least 10 boxes contain exactly one
point. This implies that of the eight points closest to the corners of S, at least two are in a box
that does not contain another point." (k boxes with ≥2 points, u uncovered points: 13+k ≤ 16−u, so
k+u ≤ 3; a box holds at most the two points of one corner pair, so corner points in singleton boxes
≥ 8 − 2k − u ≥ 5 − k ≥ 2.) Two such boxes are called **corner-restricted**.

**Replacement lemma** (Lemma 10, p. 5): "If a box T covers A(1,0.914) but not any other point in
Figure 2, then it contains (1.12,1), (1,1.74), and (1.87,0.76) (and hence the convex hull of these
four points). Proof: This follows as the set stays unavoidable if A is replaced by any of the other
points." This is Stromquist's trick (S84-II p.10, S03 §2 step 1): a box that is the *only* box
containing A, and contains no other point of the set, must contain every point A′ for which
(P∖{A})∪{A′} is still unavoidable. It is the discrete form of B16 Theorem 8. **It needs a singleton
box, which is what the counting slack buys.**

**Lemma 11** (p. 7): "assume that a box T does cover point A, but does not cover B, C, or D, and, in
addition, does not cover the point (1,1). Then the box contains the points (1.82,1) and (1.96,0.76)."

**Case tree** (branch on the positions of the two corner-restricted boxes):

- **§3.1 Non-adjacent** (the second corner-restricted box does not cover B; WLOG one covers A).
  Configuration: 10 fixed points (1,1.74), (2,1.74), (1.6,1) and mirrors under x↦4−x, y↦4−y, plus 4
  "alternative pairs" {A,B} and mirrors, choosing from each pair the point in a corner-restricted box
  if any (14 points). All regions of Fig. 3 are covered by lemmas except the two **critical regions**
  [1,2]×[1.74,2.26] and [2,3]×[1.74,2.26] (the text's "[1,2]×[1.74,1.26]" is a typo). Each
  corner-restricted box holds its pair point *and* one of (1,1.74),(1.6,1) or mirrors (Lemma 10):
  "Hence four of the 14 points are covered by two boxes. Of the remaining 11 boxes, at least one
  cannot contain a point and thus its centre is in one of the above regions." Then a second
  configuration (Fig. 4): 11 points (1.6,1),(2.26,1),(3,1.6),(1,1.74),(1.78,1.74),(2.26,2) + mirrors
  in y↦4−y + the 4 pairs = 15 points, fully unavoidable; the box centred in [1,2]×[1.74,2.26]
  contains (1.78,1.74) and (1.78,2.26) by Corollary 3; "six of the 15 points are contained in three
  boxes, and as the set is unavoidable, only nine more boxes can be packed inside S, for a maximum of
  twelve."
- **§3.2 Adjacent** (boxes T_A ∋ A, T_B ∋ B). At most one covers (1,1), so Lemma 11 applies to the
  other (WLOG T_B). Sets known to be covered:
  S_A = {A, (1.13,1), (1.4,1), (1.74,1), (1.87,0.76)} ⊂ T_A, S_B = {B, (1,1.13), (1,1.82), (0.76,1.96)} ⊂ T_B.
  Base configuration (Fig. 5): S_A ∪ S_B ∪ {(2.5,1),(3,1)} ∪ ({1,2,3}×{1.82,2.36,3.09}); lemmas apply
  everywhere except the four rectangles R₁..R₄ bounded by the 3×3 grid (R₁,R₂ = [1,2],[2,3] × [2.36,3.09],
  height 0.73; R₃,R₄ = [2,3],[1,2] × [1.82,2.36], height 0.54, numbered clockwise from top left). "The
  above list contains only 10 points that are not already covered by T_A or T_B … hence in order to
  pack 13 boxes into S, there must be a box T′ not containing any point." "Small enough" set = "the
  number of points not known to be covered by already placed boxes and the number of those placed
  boxes is smaller than 13." By Corollary 3, T′ contains the midpoints S_i of the long edges of R_i.
  - **3.2.1 T′ in R₁**: 18 points S_A ∪ S_B ∪ S₁ ∪ {(2.5,1),(3,1),(2,1.82),(3,1.82),(2.3,2.5),(3,2.5),
    (1,3),(2.3,3),(3,3)} plus a **sliding point** Z ∈ {1}×[2.36,2.92] ("non-avoidable regardless of
    the actual placement of Z"; high Z uses Lemma 5 with a ≤ 0.96, b = 0.76). Corollary 7: left
    intersection of T′ ≥ 0.09, so T′ contains (1,3) or meets {1}×[2.45,2.91]; place Z there. 9 free
    points + 3 boxes = 12 < 13.
  - **3.2.2 T′ in R₂**, two subcases: (1) T′ ∋ (2,2.72): 8 free points {(2.5,1),(3.1,1),(2.5,1.6),
    (1.6,1.95),(3.1,2),(1,2.72),(1,3),(1.8,3)} + 3 boxes = **11**. (2) T′ ∌ (2,2.72): "left
    intersection smaller than 0.37. By [Corollary] 7 the right intersection must be at least 0.45,
    which implies that T′ contains (3,2.64)"; 9 free points + 3 = 12.
  - **3.2.3 T′ in R₃**: right intersection ≥ 0.28 ⇒ T′ ∋ (3,2.10); "either its left intersection
    exceeds 0.32 or its right intersection exceeds 0.50. Hence T′ also covers one of (2,2.14) or
    (3,1.86)"; 10 listed − 1 = 9 free + 3 = 12.
  - **3.2.4 T′ in R₄**: T′ ∋ (1,2.10), (2,2.10); 9 free + 3 = 12.

Depth: 2 (adjacent/non-adjacent) × (1 | 4 regions, one with 2 subcases) = 6 leaves. The hard part is
§3.2: bespoke 18–20-point sets per leaf, with the wall-edge intersection bookkeeping (Corollary 7) used
to force T′ onto specific wall-adjacent points such as (1,3), (3,2.64), (3,2.10). Rotations are fully
handled (all lemmas are rotation-free). The **budget table** (free points + placed boxes; a
contradiction needs < 13) is: 12, 12 (§3.1 two steps), 12 (§3.2 preamble), 12, 11, 12, 12, 12. Every
branch except 3.2.2(1) is exactly one unit short of what n = 12 would need (see §7).

---

## 3. Nagamochi 2005

### 3.1 Results (N05 pp. 1–3)

"For two positive real numbers a and b, let ν(a,b) denote the maximum number of unit squares that can
be packed into the inside of an a′×b′ rectangle R with a′ < a and b′ < b."

> **Theorem 1.** For real numbers a, b ≥ 2, ν(a,b) < ab − (a+1−⌈a⌉) − (b+1−⌈b⌉).

> **Theorem 2.** (i) For any positive integer N such that N ∈ {n², n²−1, n²−2} for some integer
> n ≥ 1, s(N) = n holds. (ii) For any integer N ≥ 4 such that N ∉ {n², n²−1, n²−2}, s(N) ≥
> √(N − 2⌊√N⌋ + 1) + 1 > √N.

Remarks on the formula. Write a = ⌈a⌉ − 1 + α with α ∈ (0,1]; then a+1−⌈a⌉ = α. So the penalty is the
*fractional part* of each side (with integers counting as 1). For an integer a×b rectangle the bound
is ab − 2 ("an a×b rectangle is the smallest rectangle with aspect ratio a/b into which ab−2 unit
squares can be packed"). The hypothesis a, b ≥ 2 is essential (the construction needs the corner
squares [0,1]² and the lines at distance 1 from the walls); **nothing is claimed for strips of height
< 2.** Instances relevant to [0,4]²:

| a | b | bound on ν(a,b) | meaning (rectangles strictly smaller than a×b) |
|---|---|---|---|
| 4 | 4 | < 14 | ≤ 13 unit squares in any square of side < 4 (s(14) ≥ 4) |
| 3+α | 3+α | < 9 + 4α + α² | ≤ 11 for α < √7−2 ⇒ s(12) ≥ √7+1 ≈ 3.6458 (Theorem 2(ii)); ≤ 12 for α < 2√2−2 ⇒ s(13) ≥ 2√2+1 ≈ 3.8284 |
| 4 | 3+β | < 11 + 3β | ≤ 11 unit squares in a′×b′ with a′ < 4, b′ < 3⅓ |
| 4 | 3.5 | < 12.5 | ≤ 12 in a′ < 4, b′ < 3.5 |
| 4 | 2+β | < 7 + 3β | ≤ 7 in a′ < 4, b′ < 2⅓ (β > 0) |
| 4 | 2 | < 6 | ≤ 5 unit squares in any rectangle strictly smaller than 4×2 (the penalty a+1−⌈a⌉ jumps from 1 to 0⁺ as b crosses an integer, so ν(4,2) ≤ 5 but ν(4,2.01) ≤ 7) |
| 3 | 2+β | < 5 + 2β | ≤ 5 in a′ < 3, b′ < 2.5 |

(For n = 12 the general bound is far from 4; the point of listing them is that these are exactly the
kind of *compact rectangle capacity* statements an LP could certify or sharpen, §7.)

### 3.2 The weighted unavoidable set (N05 §3, pp. 3–5)

R = [0,a]×[0,b]. U consists of
- R* = [1,a−1]×[1,b−1], score = area of S ∩ R* (×λ²);
- L₁ = [(0.9,1),(a−0.9,1)], L₂ = [(0.9,b−1),(a−0.9,b−1)], L₃ = [(1,0.9),(1,b−0.9)],
  L₄ = [(a−1,0.9),(a−1,b−0.9)], score = 0.5 × length of S ∩ Lᵢ (×λ);
- Q = the 8 segment endpoints {(0.9,1),(a−0.9,1),(0.9,b−1),(a−0.9,b−1),(1,0.9),(1,b−0.9),(a−1,0.9),(a−1,b−0.9)},
  score 0.45 each;
- P = {(i,0.9),(i,b−0.9) : i = 2,…,⌈a⌉−2} ∪ {(0.9,j),(a−0.9,j) : j = 2,…,⌈b⌉−2}, score 0.5 each
  (2⌈a⌉+2⌈b⌉−12 points). [The paper prints "(i, a−0.9)"; b−0.9 is meant.]
Total score = (a−2)(b−2) + (a+b) + (⌈a⌉−3) + (⌈b⌉−3) = ab − (a+1−⌈a⌉) − (b+1−⌈b⌉).

> **Lemma 1.** Any unit square S inside λ⁻¹R satisfies σ(S) > 1.

"Then we have N′ < ab − … for any factor λ⁻¹ < 1, i.e., ν(a,b) < …". The proof works with a λ×λ
square, λ ∈ (1, 1.01], in the unscaled R: "each Lᵢ contributes to σ(S) by 0.5 per length and R* by 1
per area while each point in Q (resp., P) contributes to σ(S) by 0.45 (resp., 0.5)."

For a = b = 4: R* = [1,3]² (4), four segments of length 2.2 (4.4), Q (3.6), P = {(2,0.9),(2,3.1),(0.9,2),(3.1,2)} (2): total **14**.
This is the best published *weighted* certificate for [0,4]², matching Friedman's 14-point set
(DS7 Thm 4) in cost; both give s(15) ≥ 4 and neither reaches 13.

### 3.3 Technical lemmas (N05 §4, pp. 5–8), exact statements (S is a λ×λ square, λ ∈ [1,1.01])

- **Lemma 2.** "For a line L with distance h ∈ [0, (√2−1)/2) from the center of S, let c be the
  length of the intersection of S and L. Then c ≥ λ or c > 1." [(√2−1)/2 ≈ 0.207. = B16 Lemma 4.]
- **Lemma 3.** "…one corner of S touches the x-axis and S is entirely above the x-axis. For a line
  L: y = h with h ∈ (0.5, √2−0.5), let c be the length of the intersection of S and L. Then c ≥ λ or
  c > 1." [A wall-strip lemma: a square standing on a wall meets every horizontal line at height
  h ∈ (0.5, 0.914) in length > 1.]
- **Lemma 4.** "e₁, e₂ two adjacent edges of S meeting at corner v. For p₁ ∈ e₁, p₂ ∈ e₂ (≠ v), let
  c = |[p₁,p₂]| and d = area of the triangle p₁ v p₂. Then 0.5c > d." [Why a segment of weight 0.5
  per length can replace area 1 per unit near the boundary of R*.]
- **Lemma 5.** "…one corner of S touches the x-axis, S entirely above it, c > 0 the length of S ∩ {y=1},
  d the area of the triangle enclosed by S and L. Then d + 0.5c > 0.5."
- **Lemma 6.** (the hard one) "…two adjacent edges e₁,e₂ of S intersect L: y=1, point (1,1) is not in
  S, point (2,0.9) is on an edge e₂ of S. Let c, d as above, and p′ = (1, 1−c′) the crossing point of
  e₁ and the line x = 1. Then d + 0.5c + 0.5 − 0.5c′ > 1."
- **Lemma 7.** "S entirely contained in R, center in [0,a]×[0,1]. Then (i) the length c of S ∩ {y=0.9}
  is more than 1. (ii) If the center of S belongs to [0,1]×[0,1], then S contains (1,1), (1,0.9),
  (0.9,1) as its interior points. (iii) S contains at least one point in Q ∪ P. (iv) σ(S;R*) > 0."
  [(i) is a *strip* statement: any unit square in the container whose centre is within 1 of a wall
  cuts the line at distance 0.9 from that wall in length > 1 — hence at most ⌊a⌋ … boxes centred in
  the wall strip of height 1 can be "counted" on that line. (iii) is what makes the certificate
  work near walls: the points of Q ∪ P sit on the line y = 0.9 at unit spacing.]

### 3.4 Case analysis (N05 §5, pp. 8–12)

Seven cases on the position of the centre of S and how the segments cut S:
1. S ⊂ R*: σ ≥ λ² > 1.
2. centre in R*, S ⊄ R*, no Q point inside, no Lᵢ cuts two non-adjacent edges: each cut-off triangle
   has d < 0.5c (Lemma 4) ⇒ σ ≥ λ² > 1.
3. centre in R*, some Lᵢ cuts two non-adjacent edges: σ ≥ 0.5λ² + 0.5λ > 1.
4. centre in R*, a Q point inside: replace the point by a segment of equal total weight and reduce to 2.
5. centre in a corner square [0,1]²: Lemma 7(ii) ⇒ 2×0.45 + 2×0.05 + σ(R*) > 1.
6. centre in the wall strip [1,a−1]×[0,1], y=1 cuts two adjacent edges: subcases on which of P, Q,
   (1,1), (min{2,a−1},0.9) are inside; the last subcase is Lemma 6.
7. same strip, y=1 cuts two non-adjacent edges (c ≥ λ > 1): reduces to 6 or is immediate.

Partition: interior R*, four wall strips of height 1, four corner unit squares. Cases 6–7 (the wall
strip) are where all the work is; the corner is easy because of Lemma 7(ii). Case 6 shows the
mechanism by which strictness enters: the axis-parallel square [1,2]×[0,1] scores exactly
0.5·1 (L₁) + 0.5 (the P point (2,0.9), on its boundary) = 1, and only λ > 1 makes it > 1.

### 3.5 Concluding remark (p. 12)

"Our unavoidable set U can be seen as a modification of the entire area R so that the total score
becomes less than the area of R by replacing the boundary part of R with a set of points and line
segments with appropriate scores." No remark about n²−3 or n²−4.

---

## 4. Bentz 2016 (n = 22 and 33)

### 4.1 Toolbox (B16 §2, pp. 2–3)

Boxes now have side s with 1 < s ≤ 1.01. Lemmas 1–3 = B10 Lemmas 2, 4, 5 (Lemma 3 used only with
a = ½√3, b = 0.5). New:
- **Lemma 4** (Nagamochi). "If l is a line that lies within a distance of (√2−1)/2 of the centre of
  a box B, then l will intersect B with a length of more than 1."
- **Lemma 5** (Stromquist) = B10 Lemma 6 (parallel lines, total intersection ≥ min{1, 2√2−2d}).
- **Lemma 6** (new; extends Lemma 3 to a < 2√2−2). "Let 0 < a < 2√2−2, 0 < b ≤ 1, and (a,b) within a
  distance of 1 from (0,1). Then any box whose centre is in the quadrilateral Q with vertices (0,0),
  (0,1), (a,0), and (a,b) must intersect the x-axis, the point (0,1) or the point (a,b)." Proof by
  splitting Q with (a,½) and Lemma 2 on [0,a]×[0,1]. "We will use the lemma for a = 0.8, 0.4 ≤ b ≤ 1."
  [Rotated so the x-axis is a wall, this is the lemma that lets a row's end point slide from x = 0.5
  to x = 1 while the set stays unavoidable.]
- **Lemma 7** (new). "Let l be a line and P a point with a distance of more than 0.51 from l. If a box
  B covers P such that P and the center of B lie on opposite sides of l, than B intersect l with a
  length of intersection that exceeds 1." (Proof: centre within 0.505√2 of P, hence within
  0.505√2 − 0.51 < (√2−1)/2 of l; Lemma 4.)
- **Theorem 8** (continuous families). "Let S be a square with a packing P of boxes, I = [a,b], t ∈ N,
  and f_k : I → S a collection of continuous mappings, for 1 ≤ k ≤ t. Suppose further that
  1. for each i ∈ I, F_i = {f_k(i) | 1 ≤ k ≤ t} is an unavoidable set of points;
  2. if for some k, f_k(a) is not contained in a box of P, then f_k(i) = f_k(a) for all i ∈ I;
  3. if for some k ≠ l, f_k(a) and f_l(a) lie in the same box of P, then f_k(i) = f_k(a) for all i ∈ I.
  Then for all 1 ≤ k ≤ t, the image f_k(I) will either lie entirely within one box, or completely
  outside any box."
  (Proof: openness of boxes + minimal i₀ at which some point leaves its box; the unavoidable set at i₀
  forces another point to have entered that box just before, contradicting condition 3.)
  Conditions 2–3 are exactly "every moving point is the unique point of the set inside its box".
  With an exact count (t points, t boxes) they hold automatically: "as every red point and every
  blue point lies in exactly one box of P, the second and third condition of Theorem 8 are
  automatically full-filled" (p. 5). **Every unit of slack forces a case split on which point is
  uncovered / doubled.**

The two-tier architecture (p. 2): "We will first start out with systems of points containing too
many resources for a direct proof. A new technical result (Theorem 8) will allow us to use the
flexibility in our initial systems to show that any potential packing must contain a local abundance
of boxes. We will use this local 'over-concentration' of boxes to obtain a contradiction in
combination with a second resource system based on a line segment."

### 4.2 n = 33 (Theorem 9, pp. 4–6)

S = [0,6]². Two 33-point unavoidable sets: red (Fig. 1: bottom row (i, √2−½), i=1..5, then
equilateral rows, 6 rows alternating 5 and 6 points, 33 = 5+6+5+6+5+6) and blue = reflection of red in
y = 3. 33 boxes ⇒ "each box in P will contain exactly one red and blue point."
Moves (all preserve unavoidability; Theorem 8 applies): move whole rows vertically as long as adjacent
rows are ≤ ½√3 apart and the outer rows ≤ √2−½ from the walls; bring row i to the configuration F_i
minimising its distance to its neighbours (then < 0.8, because 2(√2−½) + 2·0.8 + 3·½√3 > 6); move a
6-point row horizontally by 0.1 each way; move its leftmost point from (0.5, y_i) to (1, y_i)
(Lemma 6 with a = 0.8). Conclusion: for each of the 6 rows there is a box B_i ⊃ [0.4,1]×{y_i}, all
distinct. Final resource: the segment l = {√2−½}×[0,6]. If the centre of B_i is on the wall side of
l, Lemma 5 with d = √2−½ (wall line x=0 cannot be met, so the whole ≥ 1 falls on l, > 1 for a box);
otherwise Lemma 7 with P = (0.4, y_i) (distance 0.514 > 0.51 from l). "Hence all six boxes intersect
l with a length larger than 1. However, the length of l is 6, for a contradiction."
No case split at all. Uses: 6 rows ⇔ side 6, the wall strip of width 0.914, and the fact that 6 boxes
each taking > 1 of a length-6 segment is impossible. **This is a wall-strip lemma in disguise: at most
m−1 boxes can each cover a point at distance 0.4 from a wall of length m at distinct heights.**

### 4.3 n = 22 (Theorem 11, pp. 6–8)

S = [0,5]². Points: {0.5,1,…,4.5}×{0.9,1.7,2.5,3.3,4.1} (9×5 = 45), split into red (22: rows
r₁,r₃,r₅ = {1,2,3,4}, r₂,r₄ = {0.5,…,4.5}) and blue (23: complementary). Both unavoidable (row
gap 0.8, triangle sides √(0.25+0.64) < 1; Lemma 2 with a=1, b=0.9 at the bottom; Lemma 6 at the side
walls). 22 boxes: red exact; blue has slack 1: "either exactly one blue point is not contained in a
box of P or exactly one box of P contain two blue points", and the two blue points sharing a box are
within 1.01√2, so by symmetry "all blue points with first coordinate value lower than 2 lie in a box
of P that does not contain any other blue points" [the exceptional point is pushed to the far side].
Moves as before (only horizontal moves are needed since the row spacing is fixed at 0.8); for the blue
rows the exceptional point may make one of b₁,b₃,b₅ non-movable, so that box is only known to
cover [0.5,1]×{y_i}. Result: boxes B₁..B₅ with B_i ⊃ [0.4,1]×{y_i} "except that at most one of B₁,B₃,B₅
might only cover the line segment [0.5,1]×{y_i}." Segment l = {√2−½}×[0,5] of length 5: at most 4
boxes can take > 1 each, so the exceptional box exists and is the one that fails to cover
(0.4, y_i). Two cases (Figs 4, 5), each a small semi-algebraic region computation for the centre m
of the exceptional box (on the far side of l, ≥ ½√2−½ from l, ≥ 0.5 from the two "green" points of l
that the other boxes leave free, ≤ 0.505√2 from (0.5, y_i)):
1. i = 3: region within 0.5 of (1.5, 2.5) ⇒ B₃ contains two blue points (0.5,2.5), (1.5,2.5) — contradiction with the symmetry assumption.
2. i = 1 (or 5): region within ½ of (√2−½, 1), which is denied to B₁ — contradiction.
Depth 2, 2 leaves; one unit of slack (in one colour) costs one small case analysis with explicit
region computations.

### 4.4 Explicit remarks on other n

B16 abstract/intro: "these results strongly suggest that s(m²−3) = m for m ≥ 3" and "the best lower
bounds in these cases are s(22) ≥ √15+1 ≈ 4.87298 and s(33) ≥ √24+1 ≈ 5.89898, and follow from a
general result in Nagamochi". **Neither B10 nor B16 mentions n = 12, n = 21, 32, 45, or m²−4 at
all** (grep confirmed). B10 p.1 only: "While our results strongly suggest that the same formula
holds for the intermediate values m = 5, 6, no specialized results are currently known in these
cases".

---

## 5. Skimmed sources: the lemma toolbox and remarks on n = 12

### 5.1 Stromquist 1984 memoranda and 2003 paper

- S84-I Lemma 1: "If a block B is contained in the interior of a square S with side s ≤ 2, then B
  contains the center of the square." Cor. 2: the (1,1) lemma. Lemma 2: triangle lemma. Then the
  wall-strip computation (p. 5–8): a block with one corner on y=0 and one corner above y=a meets y=a in
  length f(θ,a) = (sinθ+cosθ−a)/(sinθcosθ), minimised at θ = 45°: **2√2 − 2a**. Lemma 3 (parallel
  lines, total length > min(1, 2√2−2a)). Lemma 4: "Consider a region R bounded by the x axis, the
  line y=a, and two vertical lines separated by a distance b. Suppose that either (1) a = 1 and
  b ≤ 2√2−2 ≈ .828, or (2) b = 1 and a ≤ √2−½ ≈ .914. Then any block B whose center lies in the
  region, and which does not intersect the x axis, must contain one of the corners of the region."
  Lemma 6/7: two 8-point unavoidable sets in [0,3]² (Lemma 7: (1,2),(1½,2),(2,2),(1,1.7),(2,1.7),
  (1½,1½),(1,0.9),(2,0.9)). The n=6 proof (p. 13–18) is the prototype of Bentz's argument: 9 "key
  points", 6 blocks ⇒ ≥ 3 isolated key points; Lemma 8 ("Key points A and B cannot both be
  isolated") is proved by measuring covered/uncovered lengths of the two segments (1,0)–A and A–B;
  the centre block is then forced to contain J=(1,1.7), K=(2,1.7) via "the center block must intersect
  the lines x=1 and x=2 in segments with total length greater than 0.8".
- S84-II Lemma 3 (the f(a) lemma) with the table a = 0.828/0.853/0.96 → θ = 45°/39.514°/17.708°,
  f = 1/0.9722/0.7689; Lemma 4 (pentagon (.88,0),(.88,.90),(1,1),(2,1),(2,0): block meets the x-axis, the
  segment (.88,.90)–(1,1), or the segment (2,1)–(2,.788)) — a **segment**-valued non-avoidance lemma,
  the ancestor of B10's "sliding point".
- S03: Lemmas 1–4 as above (Lemma 4 = B10 Lemma 5), Lemma 5 (pentagon (1,0),(1,1),(2,1),(2.12,.9),(2.12,0)),
  Lemma 6 (a = √(4/5): pentagon (1,0),(1,1),(1+a/2,1.12),(1+a,1),(1+a,0), block meets the x-axis or a vertex).
  Theorem 2 (s(11) ≥ 2+2√(4/5) ≈ 3.789): 10-point set, one box must sit in a top/bottom rectangle and
  then contains three points (1,.9),(s/2,.9),(1+1/√5,1.12) of a 12-point unavoidable set. Theorem 3:
  45° packings of 11 need s ≥ 2+(4/3)√2 ≈ 3.886 (stronger triangle Lemma 7 for 0°/45° boxes).
- **S84-III p.3, explicit remark on n = 12**: "Cases with n = k²+k. Göbel was unable to improve on the
  trivial packings in the case of any n of the form n = k²+k. We now know that no nontrivial packings
  are possible in the cases of n = 2 or n = 6 (reference [a]), and it seems very unlikely that
  nontrivial packings will be found for n = 12 or n = 20. Nevertheless, the asymptotic results in
  reference [e] compel the existence of nontrivial packings with n = k²+k for all n and k sufficiently
  large." (Now known: n = 272 = 17²−17, Cleemann, DS7 Fig. 8.)
- **S84-III p.10, a compact rectangle fact**: "four unit squares can be packed into a rectangle with
  sides 1.9 and 3.95, but only by a non-45° packing. An optimal packing of this sort is shown in Figure
  9. Here θ ≈ 39.63° satisfies cos³θ − sin³θ = .95cos²θ − .90sin²θ, and the rectangle has dimensions
  1.9 and 3.9475⁺." Two axis-parallel squares at the ends (one at the top-left, one at the
  bottom-right corner) and two tilted squares in the middle. So **"at most 3 unit squares in a strip
  of height < 2 and length < 4" is false**; 4 squares fit in 1.9 × 3.9475.
- S84-I p.10, Lemma 5: the 4-point set for s(5) = 2+½√2 (Göbel's result reproved).

### 5.2 Kearney–Shiu 2002

No strip lemmas. Their contribution is the **duality** method: the 7-point set (6) for [0,3]²
{(√2−½,1),(3/2,1),(7/2−√2,1),(3/2,3/2),(√2−½,2),(3/2,2),(7/2−√2,2)} and its 90°-rotation ("green" and
"red lattices", 13 distinct points sharing the centre C); each unit square covers ≥ 1 point of each
lattice, so with 7 squares each covers exactly one of each, and Lemma 1 ("Any unit square which
covers the C-point must also cover a B-point", via the arc of the radius-½ circle inside a unit square
subtending ≥ 90°) gives s(7) = 3 at once. For s(6) = 3 (§3): at most one square covers two points of
one lattice; configurations (a) C uncovered, (b) C covered; technical Lemma 2 ("Let U be a unit square
with centre inside [0,1]². Suppose that one corner of U touches the x-axis with an edge making an angle
θ, and that the point (0,1) lies on the opposite edge of U. Then the points ((1+t²)/(1+t), 1),
(1, (1+2t−t²)/2), where t = tan(θ/2), lie on two of the edges of the square") with the inequalities
(7): (1+t²)/(1+t) ≥ 2√2−2, (1+2t−t²)/2 ≥ ½, sum ≥ 3/2; Lemma 3 (a unit square covering (3/2,3/2) with
(1,2),(2,2),(2,3/2) on three edges meets x = 1 at height ≤ 5/3). The proof is a chain of "if this
square covers X it cannot avoid Y" implications with the numbers above, i.e. again wall-line
intercept bookkeeping. This is the red/blue idea B16 uses, in 1×1-grid form.

### 5.3 Friedman DS7 (2009 version)

Technical lemmas (§4), all for a unit square u (closed containment implicit):
- Lemma 1: centre in [0,1]², u in the first quadrant ⇒ u ∋ (1,1).
- Lemma 2: "Let 0<x≤1, 0<y≤1, and x+2y<2√2. Then any unit square inside the first quadrant whose center
  is contained in [1,1+x]×[0,y] contains either the point (1,y) or the point (1+x,y)."
- Lemma 3: triangle with sides ≤ 1.
- Lemma 4: "If the center of a unit square u is contained in the rectangle R=[0,1]×[0,.4], then u
  contains a vertex of R." (compare B10 Cor. 3: for b < ½√3 a vertex *or the midpoint*.)
- Lemma 5: "If a unit square has its center below the line y=1, and is entirely above the x-axis, then
  the length of the intersection of the line y=1 with the square is at least 2√2−2."
- Lemma 6: "If a unit square has its center in the region [0,1]², does not contain either of the points
  (0,1) and (1,1), and is entirely above the x-axis, then the square covers some point (0,y) for
  ½≤y≤1 and some point (1,y) for ½≤y≤1."
- Lemma 7: same hypotheses ⇒ "the square covers either the point (0,√2−½) or the point (1,√2−½)."
Lower bounds (§5): pure point sets for n = 2,3 (1 pt), 5 (4), 8 (7), 15 (14 pts:
{(1,1),(1.6,1),(2.4,1),(3,1),(1,1.8),(2,1.8),(3,1.8),(1,2.2),(2,2.2),(3,2.2),(1,3),(1.6,3),(2.4,3),(3,3)}),
24 (23), 35 (34). For 7 and 14: "almost unavoidable" sets of 5 resp. 12 points, "at least two squares
have their centers in the regions containing question marks", 2 resp. 5 placements up to symmetry,
each with 3 resp. 11 further points, some with sub-cases (Figs 29–33; figures not in the text copy).
Table 2 (lower bounds, 2009): "11–12: 2+4/√5 ≈ 3.7888 Stromquist", "13: 3.8437 Friedman" (pre-Bentz,
unavoidable set not shown), "21: 4.7438 Friedman". Table 1 lists "12–13: 4" without "optimal". No
remark on m²−4.

---

## 6. What Bentz imports from Nagamochi, and what Nagamochi gives for our problem

- B10: Lemma 1 (the (1,1),(0.9,1),(1,0.9) corner lemma = N05 Lemma 7(ii)) — used to place the corner
  points A, B at distance 0.914 rather than 1 from the walls (actually B10 uses √2−½ = 0.914 from S84-I
  Lemma 4(2); the 0.9 of Nagamochi is the same idea). Also the intro's list of known values and the
  bound s(m²−3) ≥ √(m²−2m−2)+1 for m = 5, 6.
- B16: Lemma 4 = N05 Lemma 2 (line within (√2−1)/2 of the centre ⇒ intersection > 1); the "resource"
  vocabulary and the idea of a *line segment* as a resource ("A more complex configuration in [8]
  uses a combined system of (weighted) points, line segments, and a rectangular area").
- Neither Bentz paper uses Nagamochi's Theorem 1 as a sub-lemma (no rectangle sub-containers are cut
  out); the proofs are global point-set arguments on the whole square.
- What Nagamochi gives directly for sub-rectangles of [0,4]²: §3.1 table. The limitation for our
  purpose is a, b ≥ 2 and the linear-in-fractional-part penalty, which is far from sharp for strips
  (e.g. ν(4,2) < 6 says nothing about a 3.99×1.9 strip, where S84-III shows 4 squares fit and
  Nagamochi's method cannot say how many).

---

## 7. Synthesis for n = 12

### 7.1 Semantics first: which "LP" is provably blocked

The reported obstruction ("every cover-type certificate at the open container (0,4)² costs ≥ 16
because 16 slightly-eroded unit squares fit") is a statement about certificates that are valid for
squares of side 1−ε, i.e. **open** containment. Every published lower-bound proof works in the
opposite regime: side 1+ε for all ε > 0 (⇔ closed unit squares with closed containment). In that
regime the lower bound for any cover certificate is the maximum number of *pairwise disjoint closed*
unit squares in [0,4]² (shared boundary points are shared resources). Since 13 pairwise disjoint
closed unit squares in [0,4]² could be dilated by 1+ε about their centres inside [−ε/2, 4+ε/2]², giving
13 unit squares in a square of side (4+ε)/(1+ε) < 4, s(13) = 4 (B10) implies that bound is ≤ 12, and
s(12) = 4 would imply it is ≤ 11. **So in closed semantics no obstruction of the "16" kind is known
against a certificate of cost < 12.** The best known closed-semantics certificates for [0,4]² cost 14
(Friedman's 14 points; Nagamochi's weighted set, §3.2), and Bentz needed a 6-leaf case analysis to
get the effective cost from 14 down to 12 (§2.3). Getting to 11 by certificate alone is not ruled out
but nothing in the literature comes close.

Concretely: the machine-checkable question "what is the minimum cost of a weighted cover (points,
segments, areas) of [0,4]² for closed unit squares with closed containment?" is open and has no
16-obstruction; if the answer is < 12 the problem is solved. The verification differs from the eroded
LP only at zero-margin poses (a square touching two walls, or a wall and a grid line), where the
required statements are exactly the non-avoidance lemmas of §2.1 (semi-algebraic, provable exactly).

### 7.2 Where exactly one unit of slack is lost

Deficit = (points in the best pure unavoidable set on [0,m]²) − (n−1):

| m | best pure set | n = m²−3 | deficit | n = m²−4 | deficit |
|---|---|---|---|---|---|
| 4 | 14 (DS7 Thm 4) | 13 | 2 (B10: 6-leaf case tree) | 12 | **3** |
| 5 | 22/23 (B16 Fig 3) | 22 | 1 (B16: Thm 8 + 2 leaves) | 21 | 2 |
| 6 | 33 (B16 Fig 1) | 33 | 1 (B16: Thm 8, no cases) | 32 | 2 |
| 7 | 45 (B10 Fig 1) | 46 | 0 (pure set) | 45 | 1 |

n = 12 has the largest deficit of any case in the table; n = 45 is the natural next target of the B16
method (deficit 1, same as 22/33), not n = 12.

Inside B10 §3 the loss is visible step by step (numbers for 12 boxes in brackets):
1. 16-point set, 13 boxes: k+u ≤ 3 ⇒ ≥ 2 singleton corner boxes. [12 boxes: k+u ≤ 4, and all 8 corner
   points can sit in 4 double boxes — which is exactly what the trivial packing minus 4 squares does.
   The replacement Lemma 10 then has no box to apply to.]
2. §3.1: 14 points, 4 in two boxes ⇒ 10 points for 11 boxes ⇒ a point-free box. [10 for 10: nothing.]
   Then 15 points, 6 in three boxes ⇒ 9 for 10 boxes ⇒ ≤ 12 total. [9 for 9: nothing.]
3. §3.2: 10 free points + 2 boxes = 12 < 13 ⇒ point-free box T′. [12 < 12 fails.] Leaves: 12, 11,
   12, 12, 12 (§2.3 table). [Only 3.2.2(1) closes for 12 boxes.]
Every configuration Bentz builds has effective cost exactly 12: they are 12-certificates, and the
s(13) proof is a proof that the 12-certificates can be assembled consistently. A proof for 12 must
force one more double coverage in *every* branch (or start from a 13-cost certificate).

In the B16 method on [0,4]²: rows y = 0.9, 1.7, 2.5, 3.3 and x ∈ {0.5,…,3.5} give red = {1,2,3}×{0.9,2.5}
∪ {0.5,…,3.5}×{1.7,3.3} (14 points) and blue = complement (14 points); with 12 boxes both colours have
slack 2, versus (0,1) for 22 and (0,0) for 33. Theorem 8's conditions 2–3 then fail for up to two
points per colour (uncovered, doubled, or even tripled: (1,0.9),(2,0.9),(1.5,1.7) fit in one box), and
the final segment l = {√2−½}×[0,4] needs 4 boxes each covering (0.4, y_i) to get the contradiction
"4 × (>1) > 4"; with two rows possibly failing there is no contradiction. This is the same "one row
short" phenomenon as in B16 §4, but two units deep.

### 7.3 Lemmas that transfer verbatim to n = 12

All local non-avoidance lemmas (B10 Lemmas 1–7, B16 Lemmas 1–7, N05 Lemmas 2–7, S03 Lemmas 1–6, DS7
Lemmas 1–7) are container-independent and rotation-general. The 16-point set of B10 Fig. 2, the
14-point sets (DS7 Thm 4; the B16-style red/blue pair), Nagamochi's weight-14 set for a=b=4, the
replacement principle (Lemma 10 / Theorem 8), and the segment-resource argument (B16 §3) all apply
unchanged to [0,4]². What does not transfer is every counting step (§7.2).

### 7.4 Sub-lemmas that are compact packing statements (LP-certifiable in closed or eroded semantics)

Statements of the form "at most k unit squares in [0,a′]×[0,b′] with a′ < a, b′ < b" are compact and
strictness-free (the eroded LP is the right tool); the ones that appear implicitly in the proofs are:
- Nagamochi Theorem 1 instances (§3.1): ≤ 11 squares in a′ < 4, b′ < 10/3; ≤ 12 in a′ < 4, b′ < 3.5;
  ≤ 7 in a′ < 4, b′ < 7/3; ≤ 5 in a′ < 4, b′ < 2. An LP could sharpen each threshold (Nagamochi's
  penalty is linear in the fractional part and certainly not tight).
- Corner/strip capacities below Nagamochi's a,b ≥ 2 floor: "how many unit squares fit in
  [0,a′]×[0,h] for a′ < 4, h < 2?" Known: 4 fit for (3.9475, 1.9) (S84-III Fig. 9, θ = 39.63°); the
  threshold h*(a′) below which only 3 fit is not in the literature and is a clean LP target. Similarly
  "no 2 unit squares in a′ × b′ with a′, b′ < 2" (Göbel; = s(2) = 2 in rectangle form).
- Pose-region statements used as lemmas: "a box that contains A(1,0.914) and none of the other 15
  points of Fig. 2 contains the hull of (1.12,1),(1,1.74),(1.87,0.76)" (Lemma 10) and Lemma 11 are
  statements about a compact pose region (the closure of the set of poses avoiding 15 points and
  containing A) and are checkable by branch-and-bound; the LP analogue is "capacity of the corner
  pose-region", i.e. the maximum weight a box in that region can avoid.
- The N05 case analysis is itself a certificate check on a partition of pose space (7 cases by centre
  location and edge-cut pattern); an LP over the same mixed resources (area + segments + points) with
  closed containment is the natural generalisation and could search for a cheaper set than 14.

Statements that are **not** LP-certifiable in eroded semantics because they are strictness statements
(true for boxes, false for closed unit squares with a margin): "at most 3 boxes centred in the wall
strip [0,4]×[0, √2−½]" (each box takes > 1 of the line y = √2−½; the trivial row of 4 unit squares
takes exactly 1 each), the (1,1)-corner lemma at the corner square itself, the segment-l argument of
B16, and every "> 1" in Nagamochi. These are the "integral geometric lemmas about squares near the
container walls" and must be supplied analytically (they are all one-parameter trigonometric
inequalities of the form (sinθ+cosθ−a)/(sinθcosθ) ≥ 2√2−2a, plus Stromquist's f(a) and Nagamochi's
Lemma 6).

### 7.5 What a Bentz-style proof for 12 would have to look like

Given the budgets above, the two published templates each need one more unit:
(a) B10 template: a certificate of cost ≤ 13 on [0,4]² (closed semantics) plus a case analysis that
forces double coverage twice; the branch "four corner boxes each containing an A–B pair" (the trivial
packing minus interior squares, and its perturbations) is new and is precisely where the eroded LP
gives nothing (the eroded LP's extremal packings are the 16-tiling perturbations); a strip/corner
capacity LP on compact sub-rectangles (e.g. "≤ 3 squares in the corner region [0,2−δ]² ", "≤ k in a
wall strip [0,4]×[0,h] for h below the S84-III threshold") is the natural way to kill it.
(b) B16 template on [0,4]²: 14 red + 14 blue against 12 boxes (slack 2+2); the segment argument on
l = {√2−½}×[0,4] needs 4 rows with singleton boxes; the case split on which ≤ 2 blue and ≤ 2 red
points are uncovered/doubled is finite (points are within 1.01√2 when doubled) and the leaf
computations are the Fig. 4/5-type region computations, which are exactly compact pose-region
capacity checks an LP/branch-and-bound could do — the leaves, not the skeleton, are the machine-
checkable part.

Access notes: all figures of B10/B16/DS7 are described from the text and coordinates only (figures
not extractable); S84 memoranda were read as page images; Göbel 1979 and El Moumni 1999 not accessed;
no journal version of B16 found.
