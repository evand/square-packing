# Wall-strip lemma past the one-line threshold: a small branching proof (demo)

*2026-10-06.  Directory: `tasks/n12-chat-lemmas/strip-demo/`.  Context: `../README.md` (the chat's
h = 603553/500000 lemma, our one-line h1), `notes/chord-lemma.md` (chord formula, chord band D(θ),
the one-line proof), `notes/s13-casefree.md` §1–2 (semantics).*

**Result.**  The strip lemma holds at **h = 113/100 = 1.13**, past the one-line limit
h1 = (3√2 − 2)/2 ≈ 1.12132.  The proof has one hand-proved case (the classical one-line argument)
and one certificate case with **26 box leaves**, three horizontal lines
(y = 3/4, y = a = 457/500, y = 19/20), and 234 leaf inequalities, all checked in exact rational
arithmetic.  The same method also gives certificates at 1.14 (33 boxes) and 1.15 (50 boxes, or 37
with a finer menu of lines), and with finer menus up to 1.19 (44 / 60 / 67 / 191 boxes at
1.16 / 1.17 / 1.18 / 1.19, all exactly checked); table in §7.  Every leaf, the coverage and the
chain condition are exact.  The three structural lemmas (§3) and the leaf formulas (§5) are
proved by hand.

---

## 1. Statement and semantics

Fix **h = 113/100**.  A *unit square* is a closed square of side 1 at any angle.

> **Theorem (strip, h = 1.13).**  Let Q₁,…,Q_N be unit squares in [0,L]², each with centre height
> c_i ≤ h.  If either
>  (closed) L ≤ 4 and the Q_i are pairwise disjoint as closed sets, or
>  (open)  L < 4 and the Q_i have pairwise disjoint interiors,
> then N ≤ 3.

The closed form at L = 4 is the repo's semantics for t = 4 statements (`notes/chord-lemma.md` §0
B1, `notes/s13-casefree.md` §2); the open form at L < 4 is the packing statement used for
s(n) ≥ 4.  Both are proved directly below, with no dilation step.  The only place where they
differ is one strict versus non-strict inequality in Lemma 2.

The bottom wall is used throughout; the other walls follow by symmetry.

## 2. Notation

* Pose (θ, c): angle θ ∈ [0°, 90°] (angles mod 90°), centre (x, c).  C = cos θ, S = sin θ,
  half-width of the bounding box p = (C+S)/2 ∈ [1/2, √2/2].  Inside the container c ≥ p.
  Parametrise by t = tan(θ/2) ∈ [0,1], so C = (1−t²)/(1+t²) and S = 2t/(1+t²) are rational in t.
* Half-chords.  For the line y = z put e = z − c.  For 0 < θ < 90° the square meets the line in
  [x − ℓL(z), x + ℓR(z)] with (`notes/chord-lemma.md`, proof of Lemma 1)

      ℓR = min( (1/2 − eS)/C , (eC + 1/2)/S ),      ℓL = min( (1/2 + eS)/C , (1/2 − eC)/S ).

  At θ = 0 (and 90°) both are 1/2 when |e| < 1/2.
* Chord band (`notes/chord-lemma.md`, Cor. 2): the chord ℓL + ℓR on y = z is ≥ 1 iff
  |c − z| ≤ D(θ) = (u − u² + 1)/2, u = C + S.  D decreases from 1/2 (θ = 0) to (√2−1)/2 (45°).
* Node value of a square whose left separation line is y = z₁ and right one y = z₂:

      V(z₁, z₂) = ℓL(z₁) + ℓR(z₂),   V(END, z) = p + ℓR(z),   V(z, END) = ℓL(z) + p.

  For z₁ = z₂ = z this is the chord on y = z; for z₁ ≠ z₂ it is a "skewed width".
* a = 457/500 = 0.914 (just below √2 − 1/2 ≈ 0.914214).  A pose is **H** ("high") if
  c > a + D(θ), i.e. its chord on y = a is shorter than 1.  H poses have θ within about 6.75° of 45°
  and c > 1.1211: near-diamonds standing high in the strip.

## 3. The hand-proved part

**Lemma 0 (common lines).**  If c ≤ h then the square meets every line y = z with
h − 1/2 < z < 1 in its interior.
*Proof.*  Its vertical extent is [c − p, c + p] with c − p ≤ h − 1/2 < z and c + p ≥ 2p ≥ 1 > z. ∎

All separation heights used below lie in (h − 1/2, 1) = (0.63, 1); check K1 verifies this for the
menu.

**Lemma 1 (one order on every line).**  Let P, Q be two of the squares (closed-disjoint, or
interior-disjoint).  If on one line y = z₀ in (h − 1/2, 1) the chord of P lies to the left of the
chord of Q, then the same holds on every line y = z in that range.
*Proof.*  Otherwise pick p₀ ∈ P, q₀ ∈ Q on y = z₀ with p₀ left of q₀, and p₁ ∈ P, q₁ ∈ Q on y = z
with p₁ right of q₁ (interior points in the open form; Lemma 0).  The segments p₀p₁ ⊂ P and
q₀q₁ ⊂ Q (convexity) are graphs over the height interval between z₀ and z; the difference of
their abscissae changes sign, so they share a point: a common point of P and Q (of their interiors
in the open form).  Contradiction. ∎

So the squares carry one left-to-right order, the order of their chords on y = a, and for
consecutive squares Q_k, Q_{k+1} and every admissible z, Q_k's chord on y = z ends before Q_{k+1}'s
begins (strictly in the closed form; the two may share an endpoint in the open form).

**Lemma 2 (chain inequality).**  Let Q₁,…,Q₄ be four of the squares in that order and let
z₁, z₂, z₃ ∈ (h − 1/2, 1) be arbitrary.  Then

    V₁(END, z₁) + V₂(z₁, z₂) + V₃(z₂, z₃) + V₄(z₃, END)  <  L      (closed form)
                                                             ≤  L      (open form).

*Proof.*  Let x_k be the centre abscissae, α_k(z), β_k(z) the chord ends.  The container gives
x₁ − p₁ ≥ 0 and x₄ + p₄ ≤ L.  For each k,
x_{k+1} − x_k = ℓL_{k+1}(z_k) + [α_{k+1}(z_k) − β_k(z_k)] + ℓR_k(z_k), and the bracket is > 0
(closed form) or ≥ 0 (open form) by Lemma 1.  Sum over k and add p₁ + p₄. ∎

Hence if N ≥ 4, then for any four of the squares, in order, and **any** choice of separation heights,
the four node values sum to less than 4 (closed form with L ≤ 4) or at most L < 4 (open form).
So it suffices to exhibit, for every configuration, heights with **Σ V ≥ 4**.  This "chain
functional" sees the 2-D geometry that single lines miss: one separation can be measured low
(below the side vertex of a low near-diamond), the next one high (where a high near-diamond is
wide).  Numerically (float grid over poses, scratch) the minimum over configurations containing an
H square of the best-height chain value is 4.154 at h = 1.13, 4.113 at 1.15, 4.053 at 1.18,
4.013 at 1.20 and reaches 4 near h* ≈ 1.2071, where the three-floor-squares-plus-diamond
configuration makes it tight.  So the functional itself loses essentially nothing; the whole cost
is in certifying it with finitely many boxes.

**Case A: no square is H.**  Take z₁ = z₂ = z₃ = a.  Every non-H pose has |c − a| ≤ D(θ): the upper
side is the definition, and the lower side a − c ≤ a − p ≤ D(θ) holds because
p + D = u − (u² − 1)/2 ≥ √2 − 1/2 ≥ a (that is `notes/chord-lemma.md` Lemma 3; check K1 verifies
a ≤ √2 − 1/2 as (a + 1/2)² ≤ 2).  So every chord on y = a is ≥ 1, the end terms satisfy
p + ℓR ≥ ℓL + ℓR, and Σ V ≥ 4.  This is the classical one-line proof.

**Case B: some square is H.**  This is the certificate (§4).  Every pose lies in one of the 26 boxes
(K2), the H square's box carries the H flag (K4), the rule table fixes the three separation heights
from the boxes of each consecutive pair, every node value is at least its leaf claim (K3), and the
claims along every chain of four boxes that contains a flagged box add up to at least 4 (K5). ∎

## 4. The certificate (`certs/h113.json`) and what `check.py` verifies

A certificate is: h, a; a **menu** of separation heights {3/4, 457/500, 19/20}; **boxes**
[t₀,t₁] × [c₀,c₁] with rational ends; a **rule** Z[i][j] ∈ menu (the separation height used when box
i is immediately left of box j); and **leaf claims** "V_box(z_l, z_r) ≥ q" for the (z_l, z_r) pairs
the rule can produce (END included), with q rational (the proven bound rounded down to 10⁻⁶).

`check.py` verifies, in `fractions.Fraction` arithmetic only:

* **K1 parameters**: (a + 1/2)² ≤ 2; a and every menu height lie strictly in (h − 1/2, 1); h < 3/2.
* **K2 coverage**: sweeping the t-breakpoints, for each elementary t-interval the union of the
  c-ranges of the boxes over it covers [p_lb, h], where p_lb is the exact minimum of p on the
  interval (C + S is unimodal in θ, so the minimum is at an endpoint).
* **K3 leaves**: each claim q is ≤ an exact rational lower bound of V over the whole box.  Three
  bound types (§5) are computed and the largest is used.  For boxes touching t = 0 or t = 1 it also
  checks q ≤ 1, the true value at the axis-parallel pose, where the closed forms degenerate.
* **K4 flags**: a box is flagged when c₁ > a + D_lb, with D_lb = D(u_ub) and u_ub ≥ max(C+S) on the
  box (√2 is replaced by 14142136/10⁷ when 45° is inside).  Every H pose has c > a + D(θ) ≥ a + D_lb,
  so its box is flagged.
* **K5 chain condition**: an exact min-plus DP over (box, left height, flag) states finds the
  minimum over all 4-sequences of boxes with at least one flagged box of
  Σ claim_{b_k}(Z[b_{k−1}][b_k], Z[b_k][b_{k+1}]) (END at both ends).  Result at h = 1.13:
  **1000357/250000 = 4.001428 ≥ 4**.  `--brute` re-enumerates all 66,351 such sequences and gets
  the same minimum.

## 5. Leaf certificate types

Each leaf is a rational inequality "for all (t, c) in the box with c ≥ p(t): V(z₁, z₂) ≥ q", with
z₁, z₂ ∈ menu ∪ {END}.  It is discharged by one of three bounds (`leaf.py`):

**G (good line).**  If the box lies in the chord band of z, that is c_l ≥ z − D_lb and
c₁ ≤ z + D_lb (c_l = max(c₀, p_lb)), then V(z, z) ≥ 1, V(END, z) ≥ 1 and V(z, END) ≥ 1 (chord ≥ 1 by
Cor. 2, and ℓL, ℓR ≤ p).  Two linear inequalities in rationals.

**J (joint closed forms).**  Expanding both minima, with d = z₁ − z₂ and C² + S² = 1:

    V(z₁, z₂) = min{ (1 + dS)/C,  (1 − dC)/S,  (p − c + z₂ + dS²)/(CS),  (p + c − z₁ + dS²)/(CS) }
    V(END, z) = p + min{ (1/2 − eS)/C, (eC + 1/2)/S },  e = z − c   (mirror for V(z, END)).

Each form is linear in c, so its minimum is at c = c_l or c = c₁.  In θ: numerator at its smallest
endpoint value (it is monotone in S, C or S²), denominator at its largest (C, S monotone; CS
unimodal with maximum 1/2 at 45°); p ≥ p_lb.

**S (separate half-chords).**  ℓL(z₁) and ℓR(z₂) are bounded separately.  Each is a minimum of
terms T1(k) = (1/2 − kS)/C or T2(k) = (kC + 1/2)/S with k = ±(z − c).  These are linear in c.  On a
θ-interval each has its exact minimum at an endpoint or at the single critical point S = 2k (resp.
C = −2k), where the value is √(1 − 4k²)/2.  That value is rounded down to a rational:
"√(1 − 4k²)/2 ≥ r" is the polynomial fact 1 − 4k² ≥ 4r².

At h = 1.13 the 234 claims split as: S alone 118, J alone 55, G alone 40, several 21.  All three
types are needed.

## 6. The branch tree at h = 1.13

```
N >= 4: take four squares, order them by their chords on y = a (Lemma 1)
├── A. no square is H            -> all separations on y = a; each V >= 1    [hand, = one-line proof]
└── B. some square is H          -> boxes b1..b4 of the four poses (26 boxes, K2), >= 1 flagged (K4);
                                    separations z_k = Z[b_k][b_k+1] in {3/4, 0.914, 19/20};
                                    sum of leaf claims >= 4 (K3 + K5)        [exact check]
```

The 26 boxes (index = position in `certs/h113.json`; "claims" = number of leaf inequalities):

| # | t = tan(θ/2) | θ (deg) | c | flags | claims |
|---|---|---|---|---|---|
| 2 | [0, 3467/80000] | [0.00, 4.96] | [1/2, 463/800] | good(a) | 15 |
| 24 | [0, 3467/40000] | [0.00, 9.91] | [463/800, 113/100] | good(a) | 15 |
| 3 | [3467/80000, 3467/40000] | [4.96, 9.91] | [1/2, 463/800] | good(a) | 11 |
| 9 | [3467/40000, 10401/80000] | [9.91, 14.82] | [1/2, 263/400] | good(a) | 11 |
| 23 | [3467/40000, 3467/20000] | [9.91, 19.67] | [263/400, 113/100] | good(a) | 15 |
| 10 | [10401/80000, 3467/20000] | [14.82, 19.67] | [1/2, 263/400] | good(a) | 7 |
| 14 | [3467/20000, 3467/16000] | [19.67, 24.45] | [1/2, 727059/1000000] | good(a) | 3 |
| 13 | [3467/20000, 10401/40000] | [19.67, 29.15] | [727059/1000000, 163/200] | good(a) | 3 |
| 4 | [3467/20000, 3467/10000] | [19.67, 38.24] | [163/200, 389/400] | good(a) | 3 |
| 5 | [3467/20000, 3467/10000] | [19.67, 38.24] | [389/400, 113/100] | good(a) | 8 |
| 15 | [3467/16000, 10401/40000] | [24.45, 29.15] | [1/2, 727059/1000000] | good(a) | 3 |
| 0 | [10401/40000, 3467/10000] | [29.15, 38.24] | [1/2, 163/200] |  | 11 |
| 1 | [3467/10000, 4851/10000] | [38.24, 51.76] | [1/2, 911647/1000000] |  | 8 |
| 21 | [3467/10000, 4851/10000] | [38.24, 51.76] | [911647/1000000, 113/100] | **H** | 8 |
| 16 | [4851/10000, 43957/80000] | [51.76, 57.57] | [1/2, 219153/250000] |  | 11 |
| 6 | [4851/10000, 24553/40000] | [51.76, 63.09] | [219153/250000, 113/100] | good(a) | 8 |
| 17 | [43957/80000, 24553/40000] | [57.57, 63.09] | [1/2, 219153/250000] |  | 7 |
| 11 | [24553/40000, 10851/16000] | [63.09, 68.29] | [1/2, 374959/500000] | good(a) | 3 |
| 22 | [24553/40000, 14851/20000] | [63.09, 73.19] | [374959/500000, 113/100] | good(a) | 3 |
| 12 | [10851/16000, 14851/20000] | [68.29, 73.19] | [1/2, 374959/500000] | good(a) | 7 |
| 7 | [14851/20000, 64553/80000] | [73.19, 77.80] | [1/2, 344689/500000] | good(a) | 11 |
| 25 | [14851/20000, 34851/40000] | [73.19, 82.13] | [344689/500000, 113/100] | good(a) | 15 |
| 8 | [64553/80000, 34851/40000] | [77.80, 82.13] | [1/2, 344689/500000] | good(a) | 7 |
| 19 | [34851/40000, 74851/80000] | [82.13, 86.19] | [1/2, 163/200] | good(a) | 11 |
| 18 | [34851/40000, 1] | [82.13, 90.00] | [163/200, 113/100] | good(a) | 15 |
| 20 | [74851/80000, 1] | [86.19, 90.00] | [1/2, 163/200] | good(a) | 15 |

"good(a)" means the box lies in the chord band of y = a (leaf type G at z = a).  Only box 21 (the
high near-diamonds, θ ∈ [38.2°, 51.8°], c ≥ 0.9116) can hold an H pose.  The rule is mostly
y = a (570 of the 676 entries).  Separations next to the H box use y = 19/20 or y = a, depending
on the neighbour (40 entries are 19/20).  Separations next to the low near-diamond boxes 0, 1, 16
and 17 often use y = 3/4 (66 entries), just above the side vertices of a near-diamond resting on the
floor, where it is widest.  A type-level rule (one height per pair of box *types* H / not-good(a) /
good(a)) is not enough on these boxes: its best value is 3.963 < 4 (float).  The rule really is a
table over box pairs.

Why the two extra lines.  The one-line proof is tight at two poses at once: the 45° square on the
floor (chord exactly 1 on y = a) and the 45° square at height h1 (chord exactly 1 on y = a).  Past
h1 the high diamond has chord < 1 on y = a.  It has slack on a higher line (19/20), and the floor
diamond has slack on a lower line (3/4).  Lemma 1 lets each separation use its own line, so a
configuration that contains both kinds pays for the high diamond's deficit with the low diamond's
surplus or with an end bonus (p > ½ℓ).

## 7. h versus certificate size (first points of the curve)

Search: `search.py` (float-guided bisection of boxes, rule optimised by best-response sweeps, then
box merging and rule simplification).  Every row was then verified by `check.py` (exact).
"three" = menu {3/4, a, 19/20}; "fine" = menu {a} ∪ {0.64, 0.65, …, 0.99} ∩ (h − 1/2, 1).
Leaf inequalities = Σ over boxes of the (z_l, z_r) pairs the rule can produce.

| h | menu of lines | boxes (leaves) | boxes before merge | boxes that may hold H | leaf inequalities | exact min chain | certificate |
|---|---|---|---|---|---|---|---|
| 1.125 | three | 25 | 52 | 1 | 212 | 4.000670 | `certs/h1.125_three.json` |
| 1.125 | fine (37) | 19 | 27 | 1 | 489 | 4.000105 | `certs/h1.125_fine.json` |
| **1.13** | **three** | **26** | 44 | 1 | **234** | 4.001428 | **`certs/h113.json`** |
| 1.13 | two {a, 19/20} | 45 | 115 | 3 | 189 | 4.000100 | `certs/h1.130_two.json` |
| 1.13 | fine (37) | 21 | 32 | 1 | 654 | 4.000115 | `certs/h1.130_fine.json` |
| 1.14 | three | 33 | 51 | 1 | 334 | 4.000230 | `certs/h1.140_three.json` |
| 1.14 | fine (36) | 27 | 43 | 1 | 700 | 4.000102 | `certs/h1.140_fine.json` |
| 1.15 | three | 50 | 80 | 3 | 462 | 4.000124 | `certs/h1.150_three.json` |
| 1.15 | fine (35) | 37 | 49 | 3 | 1452 | 4.000155 | `certs/h1.150_fine.json` |
| 1.16 | three | search failed at the 700-box cap (float value 3.988) | | | | | |
| 1.16 | four {3/4, a, 19/20, 49/50} | 52 | 99 | 4 | 700 | 4.000196 | `certs/h1.160_four.json` |
| 1.16 | fine (34) | 44 | 63 | 4 | 1671 | 4.000168 | `certs/h1.160_fine.json` |
| 1.17 | fine (33) | 60 | 85 | 7 | 3088 | 4.000147 | `certs/h1.170_fine.json` |
| 1.18 | fine (32) | 67 | 88 | 9 | 4505 | 4.000113 | `certs/h1.180_fine.json` |
| 1.19 | fine (31) | 191 (no merge, no rule simplification) | 191 | 20 | 102883 | 4.000795 | `certs/h1.190_fine_nomerge.json` |
| 1.20 | fine | not reached: search stopped by hand after ~30 min, without reaching 4 | | | | | |

"three" = {3/4, a, 19/20}.  "boxes before merge" is the bisection output; merging (and the rule
simplification that cuts the leaf inequalities) costs O(m²) DP solves and was skipped at 1.19, so
that row is not comparable with the others: compare 191 with the before-merge column (88 at 1.18).
At h1 < h ≤ 1.15 one box holds every H pose, and about 20–30 boxes suffice.  The H region and the
count grow together after that.

Reading: the box count grows slowly at first and then sharply as h → h* ≈ 1.20711, where the slack
of the chain functional goes to 0 (§3: 0.154 at 1.13, 0.013 at 1.20).  The chat's 8.2M leaves at
h* − 7·10⁻⁷ sit at the far end of this curve.  The size here is set by how well boxes certify the
functional, not by the functional itself, and the search heuristics are untuned (simple bisection,
greedy merging).

## 8. What is proved, and how

| item | status |
|---|---|
| Lemma 0, 1, 2 (common lines, order, chain inequality), closed L ≤ 4 and open L < 4 | hand proof (§3), elementary; Lean: convexity + IVT on two segments, then a telescoping sum |
| Half-chord formula, chord band D(θ), Cor. 2 | `notes/chord-lemma.md` Lemma 1 / Cor. 2; already in `lean/Sqpack/Chord.lean` (`sq_mem_iff_slab`, slab lemmas) |
| Case A (one line) | hand proof = `notes/chord-lemma.md` Lemma 3 with a = 457/500; K1 checks a ≤ √2 − 1/2 exactly |
| Leaf bound formulas G, J, S (§5) | hand derivation (algebra + monotonicity/unimodality on an interval); `sanity.py` cross-checks them in floats against the polygon itself |
| Every leaf claim, coverage, flags, chain condition at h = 1.13 (K1–K5) | **exact** rational check, `check.py certs/h113.json --brute` |
| Rows of the table in §7 | **exact** (`check.py` on each `certs/*.json`) |
| Search heuristics, the float values quoted for comparison | float, not part of any proof |

Not done: a Lean formalisation.  Plan for it: (i) Lemmas 0–2 in Lean; (ii) a lemma per bound type,
each stated for a symbolic box and proved once; (iii) a `decide`/`norm_num` run over the 234
concrete claims; (iv) K5 either as the 66,351-case enumeration or, smaller, as a table of min-plus
partial sums φ_k(b, z, flag) (positions × boxes × left heights × flag) together with the
inequalities φ_{k+1}(b′, Z[b][b′], ·) ≤ φ_k(b, z, ·) + claim_b(z, Z[b][b′]): about 16k inequalities
between rationals.

## 9. Negative controls (`controls.py certs/h113.json`)

* **N1, h = 12081/10000 (≈ h* + 0.001).**  An exact counterexample: axis-parallel floor squares on
  x ∈ [0, 1], [1.0001, 2.0001] and [2.9991, 3.9991] (gap 0.999), plus a square at t = 41421/100000
  (θ ≈ 44.9997°) over the gap with centre height 1.207609; container side 3.9995.  Pairwise
  disjointness is verified by the exact separating-axis test.  By Lemma 2 this configuration has
  Σ V < 3.9995 for **every** choice of separation heights; computed exactly, the maximum over the
  1/1000 grid of heights is 3.997.  The four poses sit in four boxes of any would-be certificate at
  this h, in the order A, B, diamond, C, and the diamond's box is flagged.  Leaf claims can only be
  lower than the true values, so K5 fails for every possible certificate at this h, whatever its
  boxes, menu or rule.
* **N2.**  The h = 1.13 certificate with every separation on y = a (the single-line proof) and exact
  leaves: REJECT, min chain 3.9257 < 4.
* **N3.**  One box deleted: REJECT (coverage gap).
* **N4.**  h raised to 1.131 with the same boxes: REJECT (coverage).
* **N5.**  a = 0.915 > √2 − 1/2: REJECT (one-line case not valid).
* **N6.**  One leaf claim of the H box raised to (sampled true minimum) + 0.01: REJECT (leaf
  check).  The float sampler also shows a pose where the claim is actually false.

`sanity.py certs/h113.json 5000`: over 5000 random poses per box (including box corners), every
claim is ≤ the node value computed directly from the polygon (minimum margin 0, attained at
axis-parallel poses where the claim is exactly 1).

## 10. Files

* `leaf.py`: the bounds G, J, S; exact (Fraction) or float mode.
* `check.py`: exact checker K1–K5 (`--brute` adds the full enumeration).
* `search.py`: float-guided certificate search; `table.py`: the h-sweep of §7.
* `controls.py`: negative controls; `sanity.py`: float polygon cross-check.
* `certs/h113.json`: the h = 1.13 certificate of §6; the other `certs/h1.xxx_<menu>.json` are the
  §7 rows (`table.py three|fine H...` regenerates the three/fine rows; the `two` and `four` rows
  come from `search.py 113/100 --menu 19/20` and `search.py 116/100 --menu 3/4,19/20,49/50`; the
  1.19 row from `search.py 119/100 --nomerge --nosimplify`).

Reproduce: `python3 check.py certs/h113.json --brute && python3 controls.py certs/h113.json`
(well under a second each), or `python3 search.py 113/100 --menu 3/4,19/20 --out x.json` (≈1 s).
