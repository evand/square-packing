# Explicit power-law constructions for s(k² − c) < k: Ryu 2026 first, the classical constructions as background (2026-10-17)

Purpose: the other direction from Fang's theorem s(n² − n + 1) < n (n ≥ 15; jlevy/squares#603).  Fang's gluing lemma
adds wastes, so it gives only a *linear* excess.  We want an explicit power-law construction for side just below an
integer k, with explicit constants and threshold, in a recursive form ("if these regions hold this many squares, then
the square holds this many") that Lean can check for all k.  The centre of this note is Sungjoon Ryu's 2026 preprint
(§2), which already does this.  The classical waste-bound constructions follow as background (§5).

Labels as in `waste-bounds-2026-10-12.md`: **[read]** = read in the primary source (section or lemma cited); [S] =
secondhand (jlevy's register, an issue thread); [I] = my inference or my own computation.

Notation: c*(k) = max{c : s(k² − c) = k} (Ryu's "deficiency"; our `search/FRIEDMAN.md`).  A packing of N squares in
a side < k shows c*(k) ≤ k² − N − 1 (Ryu Lemma 2.2).  Below, "c" of a construction means k² − N, so the bound it gives
is c*(k) ≤ c − 1.  Fang's theorem gives c = k − F(k) with F(k) = 1 + max{p + 2q : 16p + 21q ≤ k − 15} ≈ 2(k − 15)/21,
i.e. c_F ≈ 19k/21.

## 0. Summary

* **Ryu, "k2-minus-c-upper" v1.1 (10 Oct 2026; DOIs 10.5281/zenodo.23256655 (v1.0), .23278927 (v1.1)) [read, TeX at
  commit c6df2f4].**
  * **Theorem 1.3:** an explicit L-shaped packing of N(k, b) = (k − b)² + b(R(k) + R(k − b)) squares in a side < k,
    for all 2 ≤ b ≤ k.
  * **Theorem 1.4:** c*(k) ≤ 8⌈√(k − 4)⌉ − 1 for every k ≥ 6.
  * Theorem 1.2 (staircase ends): c*(k) < 20.668k^{2/5} + 0.623 for k ≥ 1.6·10⁷.
  * Theorem 1.1 (chained wall fillers, twelve tiers): c*(k) < C_i k^{3/8} + a_i, e.g. 42.09k^{3/8} + 0.05 for
    k ≥ 2.48·10¹² (v1.0: 43.06k^{3/8} for k ≥ 4·10²⁶).
  * Theorem 1.6: exact certificates, e.g. c*(10⁵) ≤ 1583 and c*(10⁸) ≤ 24790.
  * Unrefereed, written with Claude, no human review.  jlevy re-derived it in one AI review and rebuilt Theorem 1.3's
    packing in exact rationals for all 2 ≤ b ≤ k ≤ 30 [S].
* **Theorem 1.3 is closed-form throughout and is the thing to formalize.**  Every square has rational coordinates
  given by a formula in (k, b, i, j).  The tilt is tan(θ/2) = p/q with p = 14kb² + 1 and q = bp − 1.  The proof is
  seven short lemmas (§2.2).  Estimated Lean effort [I]: a few hundred to ~1500 lines on top of our `Sqpack.Spec.Packs`.
  The only analytic step (a mean-value bound on the tilt) can be replaced by polynomial inequalities in t.
* **Small improvement found here [I, exact computation]: fill the four band ends with axis-aligned rows.**  These
  rows lie inside Ryu's own end regions Reg(y), so his Lemma 6.1 gives disjointness for free.  This cuts each end's
  cost from ≈ b to ≈ b/2, so **c ≈ 4k/b + 2b, i.e. ≈ 4√(2k) ≈ 5.66√k, instead of 8√k**.  Two consequences:
  * The window where Fang's linear bound is better shrinks from 15 ≤ k ≤ 78 (and k = 86), per jlevy's register, to
    **k ≤ 37** with the best b, or k ≤ 39 with b = ⌊√(2k)⌋.
  * At k = 100: c = 56 vs Ryu's 80 and Fang's 91.
  It stays Θ(√k) and is beaten by Ryu's staircase certificate by k = 10⁵ (≈ 1790 vs 1584).
* **Lower bounds (#414): also Ryu.**  Two preprints opened on jlevy/squares#414 (7 Oct 2026), both OPEN and reported,
  not refereed:
  * [Q] k² − M(k) ≥ 0.1k^{1/4} for k ≥ 4.62·10¹² (squarepacker/k2-minus-c-quarter, doi 10.5281/zenodo.23211207);
  * [T] k² − M(k) > 0.06k^{1/3} for k ≥ 3.5·10¹⁷, plus √k-type bounds 10⁻⁴√k on [10⁸, 2.5·10¹⁵]
    (squarepacker/k2-minus-c-cube-root, doi .23212059).
  Here M(k) is the closed-disjoint maximum, and c*(k) = k² − M(k) − 1.  Both build on his log-k paper [R] (v1.2,
  0.0353 log k), including its computer-assisted Lemma 4.10.  jlevy imported both as "reported" (PR #439), with one AI
  review of the proof chains (PR #441), interval re-checks of the constants, no human review and no independent
  box-labelling [S].  **So, as reported: k^{1/3} ≪ c*(k) ≪ k^{3/8}.**
* **Our files are out of date.**  `search/FRIEDMAN.md` (§0, §9.0: "c*(k) = O(k^{0.6})", "any threshold ≳ c^{5/3}") and
  WISHLIST A13 predate Ryu's upper bounds.  As reported, thresholds are ≳ c^{8/3}.  The lower bounds, if right, give
  s(k² − c) = k for k ≳ (c/0.06)³, and Ryu's quarter-power corollary gives k ≥ max(4.62·10¹², 10⁴c⁴ + 1).
* The classical constructions (§5) are background only: each analyses the generic fractional part.  At a side just
  below an integer, their first level collapses to Ryu's Theorem 1.3.  None states anything special for k⁻
  (consistent with Ryu §10).

## 1. Sources

| Source | Status | Where |
|---|---|---|
| S. Ryu, "k2-minus-c-upper" (upper bounds on c*(k)), v1.1, github squarepacker/k2-minus-c-upper @ c6df2f4 (10 Oct 2026) | **[read]**: §§1–5 in full; §6 (waste identity) in full; §7 setup, the staircase construction (S1)–(S6) and the end lemma statement; §§10–12 in full; §8 (3/8) only via the theorem statement and tier table | `paper/paper.tex`, fetched by `gh api` to scratch (data only; no code run) |
| jlevy/squares `packing/frontier/asymptotic-waste-bounds.yaml` | **[read]** (copy supplied by the coordinator) | register entries for Ryu, Fang, the classical bounds |
| jlevy/squares issue #414 (Ryu, two lower-bound preprints) | **[read]**: body and jlevy's 12 Oct reply | `gh issue view 414` |
| S. Ryu, cube-root preprint [T], squarepacker/k2-minus-c-cube-root | abstract **[read]** | `paper/paper.tex`, fetched |
| S. Ryu, quarter-power preprint [Q] and log-k preprint [R] | statements via #414 and the register [S]; [R] was reviewed in `tasks/issue3-k2-minus-c/review-2026-10-10.md` | — |
| Erdős–Graham 1975 | **[read]**, all 5 pp. incl. figures | `../private/literature/1-s2.0-0097316575900990-main.pdf` |
| Kearney–Shiu, EJC 9 (2002) #R14 | **[read]** | jlevy archive `kearney-shiu-2002-…raw.md` |
| Chung–Graham 2009 (JCTA 116) and 2020 (DCG 64) | **[read]** (preprint texts) | jlevy archive |
| Wang–Dong–Li, arXiv:1603.02368 v2 | **[read]** | jlevy archive |
| Bui, arXiv:2508.04603 v2 (15 Mar 2026) | **[read]** | arXiv PDF, scratch |
| McClenagan, arXiv:2602.01484 v1 | **[read]** | jlevy archive |
| Arslanov–Bui, DCG 2025 note on Chung–Graham 2020 | not reached; content via Bui §1 [S] | — |

The jlevy archive is `~/math/_untrusted-third-party/jlevy-squares/packing/resources/papers/`, read as data.  Its
clone (10-03) predates Ryu's upper-bound paper.

## 2. Ryu's Theorem 1.3: the L-shaped family, step by step

### 2.1 Construction [read, §§3–5]

Fix integers 2 ≤ b ≤ k.  Put D = b² + 1, γ₀ = (b² − 1)/D, h₀ = (3b² − 1)/D, and R(L) = max(0, ⌈(L − h₀)γ₀⌉).

1. **Block** (§3).  Tilt θ, γ = cos θ, σ = sin θ, u = (γ, σ), v = (−σ, γ).  The row B(P) is the b × 1 rectangle
   ∪_{i<b} {P + αu + βv : α ∈ [i, i+1], β ∈ [0, 1]}: b unit squares in a row, tilted by θ, with lowest vertex P.
   * Width f(θ) = bγ + σ, height h(θ) = bσ + γ, row spacing d(θ) = sec θ (Lemma 3.1, projections).
   * **Lemma 3.2 (width):** at θ₀ = 2 arctan(1/b), cos θ₀ = (b² − 1)/D and sin θ₀ = 2b/D (rational), f(θ₀) = b,
     h(θ₀) = h₀ and d₀ = D/(b² − 1); and f(θ) < b for θ ∈ (θ₀, π/2].  This is the "reflected row": a b × 1 row tilted
     by 2 arctan(1/b) is exactly as wide as an untilted one, and any further tilt makes it narrower.
   * **Lemma 3.3 (stacking):** rows P_j = P₀ + (0, j sec θ) touch along long sides.  In the rotated frame,
     Q_i(P_j) = [i + j tan θ, i + 1 + j tan θ] × [j, j + 1], so squares are disjoint by interval disjointness.
   * **Lemma 3.4 (one band):** R rows with P_j = (a + σ, j sec θ) lie in [a, a + f(θ)] × [0, L] if (R − 1)d + h ≤ L.
2. **L-shape** (§4).  For θ ∈ (θ₀, π/2): ε = b − f(θ) > 0 and S = k − ε < k.
   * Grid: (k − b)² axis squares in [0, k − b]².
   * Right band in A = [k − b, S] × [0, S] (length S), with R₁ rows.
   * Top band: the reflection (x, y) ↦ (y, x) of a band of length k − b, with R₂ rows, in T = [0, k − b] × [k − b, S].
   * Length conditions (L1) (R₁ − 1)d + h ≤ S and (L2) (R₂ − 1)d + h ≤ k − b.
   * **Lemma 4.1:** the three regions have disjoint interiors, so this is a packing of (k − b)² + b(R₁ + R₂) in
     [0, S]².
3. **Explicit angle** (Lemma 4.2).  M = 14kb² + 1, p = M, q = bM − 1, θ = 2 arctan(p/q).  Then cos θ and sin θ are
   rational, and **R₁ = R(k), R₂ = R(k − b) satisfy (L1), (L2)**.  Proof: θ − θ₀ ≤ 2/(b(bM − 1)) ≤ 1/(7kb⁴) =: τ (the
   derivative of arctan is ≤ 1).  The slack g(L) = L − h₀ − (R − 1)d₀ is a positive rational with denominator
   (b² − 1)D, hence ≥ 1/(b⁴ − 1).  Derivative bounds on [θ₀, θ₀ + τ] (d′ ≤ 4, h′ ≤ b, |f′| ≤ b + 1) and the mean value
   theorem show the perturbation costs < 1/b⁴.  Lemma 4.3: R(k), R(k − b) are maximal for this family.
4. **Count** (§5).  For integer L, (L − h₀)γ₀ = L − 3 − η_L with η_L = (2L − 10)/D + 8/D².  So R(L) = L − 3 − ⌊η_L⌋
   for L ≥ 2, and the exact count (5.1) is **k² − N(k, b) = 6b + b(⌊η⌋ + ⌊η′⌋)** (equality when k − b ≥ 2), with
   η = (2k − 10)/D + 8/D² and η′ = η − 2b/D.
   * Theorem 1.3: c*(k) ≤ k² − N − 1 < 6b + 4k/b − 3.
   * Theorem 1.4: b = ⌈√(k − 4)⌉ gives η < 2, both floors ≤ 1, so c*(k) ≤ 8b − 1.
   * Proposition 1.5: no b does better than 8√k − 15.
   * Remark 5.2: the 8√k is 4k/b of body cost (row spacing sec θ exceeds 1 by ≈ 2/b²) plus ≈ b at each of four ends.
5. **Lifts and the waste identity** (§6, used for 2/5 and 3/8, and by the improvement in §3).  In a band's strip frame
   Σ = [0, w] × [0, L_b], with lift y₀ and R rows, the end region is
   Reg(y) = {0 ≤ x′ ≤ w, 0 ≤ y′ ≤ y + (x′ − σ) tan θ}.
   * **Lemma 6.1:** any packings C₀ ⊂ Reg(y₀) and C₁ ⊂ Reg(y₁) (the latter placed by the half-turn ρ of Σ) can be
     added.  Separation is by the lines β′ = 0 and β′ = R of the row frame.
   * **Proposition 6.3 (waste identity):** k² − N = (k² − S²) + Σ_bands[(R − 1) tan θ + E(y₀) + E(y₁)], where the end
     cost is E(y) = A₀(y) − |C| and A₀(y) = yw + ½γσ(1 + b²).
   * The empty end costs E(0) ≈ b.  Any end packing with E = O(b^β) gives c*(k) = O(k^{β/(1+β)}).
   * Staircase ends (§7, after Bui §3): E = O(b^{2/3}), giving 2/5.  Chained wall fillers (§8): E = O(b^{3/5}), giving
     3/8.

### 2.2 Is every piece a closed formula, and how hard is Lean? [I]

* **Theorems 1.3 and 1.4: yes, entirely.**
  * The grid squares are [i, i+1] × [j, j+1].
  * The band squares are Q_i(P_j) with P_j = (k − b + σ, j·(q² + p²)/(q² − p²)), and γ, σ rational in p, q.
  * The counts are ⌈·⌉ of rationals.
  * There is no search, no per-k data, and no interval arithmetic.
* **Theorem 1.2 (2/5):** also explicit in principle.  But the end lemma (Lemma 7.3) has a long waste bound B*, and the
  constant rests on one interval evaluation Φ(10^{−4/3}) ≤ 5.0163 (Cor. 7.4).  It applies only from k = 1.6·10⁷.
* **Theorem 1.1 (3/8):** constants from interval coverings (320 × 128 boxes per tier); valid only from k ≈ 2.5·10¹².
  Not a Lean target.
* **Proposed Lean lemma structure for Theorem 1.3** (on `lean/Sqpack/Spec.lean`'s `Packs n s` / `minSide`):
  * **T (rational tilt).**  For t = p/q ∈ (0, 1): cos(2 arctan t) = (1 − t²)/(1 + t²) and sin = 2t/(1 + t²)
    (Mathlib has `Real.cos_two_mul`, `Real.cos_arctan`).  After this, work with rationals γ, σ and γ² + σ² = 1.
  * **B1 (row/stacking).**  The squares Q_i(P_j) as images of [0,1]² under the rigid motion z ↦ P₀ + Mz, M = [u v].
    Disjointness from Lemma 3.3's coordinates: two closed unit boxes with disjoint open interiors stay disjoint
    under an isometry.
  * **B2 (one band).**  Containment from the vertex bounds (Lemma 3.1); linear in i, j.
  * **L (assembly, the "if regions pack then the square packs" lemma).**  If [0, w] × [0, S] holds X₁ and
    [0, w] × [0, k − b] holds X₂, then [0, k − b + w]² holds (k − b)² + X₁ + X₂.  It uses three interior-disjoint
    rectangles and the reflection ι.  This is the same shape as Fang's gluing lemma, and one shared "disjoint regions
    ⇒ add counts" lemma serves both.
  * **A (angle).**  Replace Ryu's mean-value argument by direct polynomial inequalities: (R − 1)d(t) + h(t) ≤ S(t) at
    t = p/q, with d, h, f rational functions of t, after clearing denominators.  Or keep the MVT, with
    `Real.arctan` bounds from Mathlib.  Either way it is uniform in k, b.
  * **C (count).**  R(L) = L − 3 − ⌊η_L⌋ (an Int.ceil identity), the exact count (5.1), and Theorem 1.4's η < 2
    (one rational inequality).
  * **Main.**  `minSide (k^2 - c) < k` for c = k² − N(k, b), all k ≥ 6.
* Effort estimate [I]: B1/B2/L are generic geometry (a few hundred lines, reusable).  A and C are a few hundred
  lines of `nlinarith` / `field_simp` / floor-ceiling lemmas.  No `native_decide` and no data.  The hardest part is
  probably the containment/disjointness plumbing for `unitSq` with centres rather than corners (Spec's convention).
* **Adding the end fills of §3** costs one more lemma.  The rows lie in Reg(y) by a linear inequality, and Lemma 6.1
  (separation by β′ ≤ 0, β′ ≥ R) is another two-line rigid-motion argument.

## 3. Improvement found here: Theorem 1.3 with filled ends ("R1′") [I]

I found this independently before reading Ryu.  It is his Theorem 1.3 packing (identical band counts: e.g. k = 100,
b = 12 gives R(100) = 96 and R(88) = 84, matching my K values), plus axis-aligned unit squares in the four end regions.

* **Lower ends (lift 0).**  Reg(0) is a thin right triangle: width w ≈ b across the band, rising to bσ ≈ 2 along it.
  The squares [x′, x′ + 1] × [0, 1] fit for x′ ≥ σ + cot θ, i.e. **⌊w − σ − cot θ⌋ ≈ b/2 − 1** of them.  A second layer
  would need width beyond w, so this is all.
* **Upper ends (lift y₁ = slack ∈ (0, sec θ]).**  Layers r = 0, 1, … get ⌊w − σ − (r + 1 − y₁)cot θ⌋⁺ squares each.
* Count: N′ = N(k, b) + 2⌊w − σ − cot θ⌋ + E_R(k) + E_R(k − b).  The end cost is ≈ b/2 per end, so
  **c ≈ 4k/b + 2b, minimized at b ≈ √(2k): c ≈ 4√(2k) ≈ 5.657√k.**
* **Checks.**  Exact `Fraction` construction with a separating-axis pairwise check (scratch scripts, not kept; the
  checker reported a planted 1/1000 shift and a duplicate).  Full checks passed at k = 30 (b = 5), 35 (b = 7),
  50 (b = 8), 100 (b = 12), with t = 1/b + 10⁻¹²: N = 869, 1193, 2460, 9944.  Counts for all k in 15–160 (best b) and
  10–3000 (b = ⌊√(2k)⌋) use the same exact formulas.  An exact all-k proof of the fills' count is not written.

| k | 20 | 30 | 35 | 40 | 50 | 60 | 80 | 100 | 150 | 200 | 500 | 1000 | 2000 | 5000 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| best b | 4 | 7 | 7 | 6 | 8 | 8 | 11 | 12 | 15 | 18 | 31 | 40 | 59 | 98 |
| c, R1′ | 24 | 30 | 32 | 36 | 40 | 44 | 50 | 56 | 70 | 80 | 127 | 180 | 254 | 400 |
| c, Ryu Thm 1.4 (8⌈√(k−4)⌉) | 32 | 48 | 48 | 48 | 56 | 64 | 72 | 80 | 104 | 112 | 184 | 256 | 360 | 568 |
| c, Fang (k − F(k)) | 19 | 29 | 33 | 37 | 47 | 55 | 73 | 91 | 137 | 182 | 453 | 906 | 1811 | 4525 |

* Against Fang: R1′ (best b) is smaller at k = 35 and for **every k ≥ 38** in 15–160, with ties at 34, 36, 37.  It is
  also smaller at every sampled k to 5000, and for all larger k since c ≈ 5.66√k.
* With b = ⌊√(2k)⌋, the last k where it fails to beat Fang is 39 (checked to 3000), and
  max (c − 4√(2k)) over 50 ≤ k ≤ 3000 is 3.68 (slowly creeping: 3.20 to 500, 3.46 to 1000; from rounding b, not
  proved bounded).
* Against Ryu's staircase certificates: at k = 10⁵, R1′ gives ≈ 5.66·316 ≈ 1790 vs Ryu's 1584, so the staircase is
  better somewhere below 10⁵ (crossover not computed).
* Relation to existing packings: Kevin Fang's repo has direct packings with the same counts at k = 35 (N = 1193,
  side 34.99926) and k = 40 (N = 1564).  They are different arrangements (208 / 270 tilted squares in the store vs
  R1′'s 392 / 396).  Our store has nothing at k = 50 and beyond.
* Ryu's Remark 5.4 (one extra top-band row into the corner, valid for 285 of 375 small cases) is a different small
  gain, also unproved in general.

## 4. Lower bounds: jlevy/squares#414 [read, issue; abstract of [T]]

* **Who and when.**  Sungjoon Ryu (@squarepacker), opened 2026-10-07; state OPEN.  jlevy replied on 2026-10-12.
* **[Q] quarter power.**  k² − M(k) ≥ 0.1k^{1/4} for integer k ≥ 4.62·10¹²; liminf (k² − M(k))/k^{1/4} ≥ 0.1696.
  Hence s(k² − c) = k for k ≥ max(4.62·10¹², 10⁴c⁴ + 1).  108 pp.
* **[T] cube root.**  k² − M(k) > 10⁻⁴√k on [10⁸, 2.5·10¹⁵]; > 7·10⁻⁵√k on [2.5·10¹⁵, 3.5·10¹⁷]; > 0.06k^{1/3} for
  k ≥ 3.5·10¹⁷; liminf ratio to k^{1/3} ≥ 0.06285.  9 pp.  New ingredient (abstract): the total inclination of the
  squares crossed by one path of [Q]'s flow is bounded, because their height ranges are disjoint, so [R]'s
  inclination-cost lemma applies to the union.  That lets one flow with tilt threshold ≈ 1/W keep every path
  quantized at every height.
* **Dependencies.**  [Q] §§3–8 and [R]'s Lemmas 4.13 and 4.15 are used as black boxes.  The constants use [R]'s
  computer-assisted Lemma 4.10 (box computations not reproduced by an independent implementation).  Without it:
  0.05k^{1/3} from 3.39·10¹⁷, liminf 0.05229.
* **Status.**  Not peer reviewed; the author reports AI-based reviews and Claude assistance.  jlevy imported both as
  *reported* bounds (PR #439) with no register entry, since they settle no n ≤ 324.  One separately prompted AI
  review accepted both proof chains "at source-review scope, conditional on the imported lemmas and numerical
  premises" (PR #441).  Constants were re-evaluated in interval arithmetic.  "None of your programs ran here, no human
  review is claimed, and no confirmation is claimed" [S, jlevy's reply].
* **Use for us.**  Every threshold is ≥ 10⁸, so it says nothing at computable k.  Logically it would make
  "s(k² − c) = k for all large k" (Friedman's conjecture A) a theorem with polynomial thresholds in c.  Our own review
  of [R] is `tasks/issue3-k2-minus-c/review-2026-10-10.md` (no gap found in Prop 4.5 or Lemmas 4.9/4.10; not a full
  check).

## 5. Background: the classical constructions (generic fractional part)

The one fact that ties them to c*(k) [I; Ryu Lemma 3.2 is the precise form]: a strip of width β = m − f (m = ⌈β⌉)
filled by m-stacks touching both sides has m cos θ + sin θ = β, so θ ≈ (1 + √(1 + 2mf))/m.
* For f = Θ(1) (every paper's case): θ ≈ √(2f/m), the waste is ≈ m^{−1/2} per unit length, and the slanted end offset
  is ≈ √(2fm), which forces a second level.
* At side k − ε with an integer core, f = ε → 0: θ → 2 arctan(1/m), the waste is ≈ 2/m per unit length, and the end
  offset is ≈ 2.  Every construction below therefore collapses, at k⁻, to Ryu's Theorem 1.3 skeleton.  Their bounds
  still hold there (they are upper bounds on θ) but are far from tight.

| Construction | Steps (condensed) [read] | Constants / threshold | All-k closed form? | Crossover vs Fang if made explicit [I] |
|---|---|---|---|---|
| **Erdős–Graham 1975** | (a) N × N core, N = [α − α^{8/11}]; (b) two strips of width β with stacks of [β] + 1 at one tilt Θ(α^{−4/11}) (Fig. 2); (c) end trapezoids T (sides D = Θ(α^{2/11}), D′ = Θ(α^{4/11})) cut into r = Θ(α^{4/11}) T_k of base m = [α^{4/11}]; (d) each T_k: an m × h_k axis block, h_k = [η_k − α^{2/11}], plus stacks of s_k + 1 touching top and bottom, each stack its own tilt (Fig. 4); wastes α^{7/11}, rα^{3/11}, α^{6/11}. | none ("α large"; Θ two-sided) | floors of powers; per-stack irrational tilts | ~700 as published [I, rough]; at k⁻ it is Thm 1.3 |
| **Kearney–Shiu 2002** (side n + δ, n² + 1) | strip of height m + δ, m = 2t² − t; two staircase trapezia of axis columns (heights m − jt) with rational slant arctan(1/t), and a parallelogram of K = n − 4t + 3 columns of m + 1 tilted by arctan(1/t); n = 4t³ − 8t² + 7t − 2, others by monotonicity | **explicit, all n ≥ 1**: δ_n < 3/(2n)^{1/3} + 3/(2n)^{2/3}; waste O(n^{2/3}) | **yes** (polynomials in t) | n/a (above an integer).  Earliest all-n explicit tilted-strip construction; its staircase end is an ancestor of Ryu's |
| **Chung–Graham 2009** | induction on Type 1 (rectangle; strips of width x^{2−2α}), Type 2 (trapezoid θ ≤ √(2/x); √(2x) Type 3's of size √(x/2)), Type 3 (θ ≤ 1/x; slabs of width x^{2−2γ}, trimmed transitions, recursion); α = (3+√2)/7, γ = 1/√2 | one c ≥ 10 plus a log factor; base "x ≤ 100" asserted; roundings suppressed | recursive, roundings suppressed | ~10⁸ (40k^{0.63} ln k vs 0.9k) |
| **Wang–Dong–Li 2016** | T₁ (rectangle; strips of width x^{3/4}), T₂ (trapezoid; A_i rectangles + T₃'s of height ½x^{1/2}), T₃ (rows of height ⌊x^{−1/6}/tan θ⌋; D_k filled with stacks at angles α_k with \|α_k − α_{k−1}\| ≤ 3(1+√2)x^{−1/2}) | claims W ≤ (16√2 + 38)x^{5/8} ≈ 60.6x^{5/8}, but via o(x^{−1/2}) Taylor terms, a dropped O(x^{13/16}) term and absorbed O(x^{5/18}); no x₀ | recursive; root tilts | ~7·10⁴ if valid |
| **Chung–Graham 2020** | x = N + ½; strips x^{4/5}; trapezoids cut into B_i; serrated right edges; stacks of fixed m crossing B_i's; wastes W₁…W₉ | **wrong at W₅** (Arslanov–Bui, via Bui [S]).  My reading [I]: the tilt change is computed by perturbing the denominator (width), not the excess f; the true Δθ ≈ x^{−2/5} makes W₅ ≈ x | — | — |
| **Bui 2508.04603** | primitive quadrilateral: tilted rows morph into axis columns (i_j ≈ θ²j/(2σ₁)), cost O((m + i_m)θ); right trapezoid filled by such quadrilaterals with drifting tilt θ_i; square → trapezoid as usual | partial (Lemmas 3–5 by interval arithmetic, d = 0.12l²/u); l, u, x₀ implicit | in principle | unknown; at k⁻ his §4 needs γ ≤ ½, but Ryu §7 adapts §3 to band ends |
| **McClenagan 2602.01484** | first level *chooses* stacks one taller (n − h ∈ (1, O(1))); trapezoid T by alternating two interleaved H_i/V_i algorithms restoring tilt φ ↔ ψ; O(h^{3/4}) per T | none (Taylor O(·)) | mostly (m′ varies per step) | unknown |

## 6. Consequences for our files (not edited here)

| File | Current claim | Should become |
|---|---|---|
| `search/FRIEDMAN.md` l.17–18, l.311–315 (§9.0) | c*(k) = O(k^{0.6}); "s(k² − c) = k forces k ≳ (c/C)^{5/3}"; "any conjecture … exponent ≥ 5/3" | Reported (Ryu, unrefereed): c*(k) ≤ 8⌈√(k−4)⌉ − 1 for all k ≥ 6 (elementary, explicit), < 20.668k^{2/5} + 0.623 for k ≥ 1.6·10⁷, < 42.09k^{3/8} + 0.05 for k ≥ 2.48·10¹².  So thresholds are ≳ c^{8/3}.  Lower: k^{1/3} (#414, reported). |
| `search/FRIEDMAN.md` l.17, l.148 | `A` claimed via Ryu log k | Add the k^{1/4} and k^{1/3} lower-bound preprints (#414) and their status. |
| `search/WISHLIST.md` A13 | second half: is the O(x^{0.6}) constant effective? | Moot for c*(k): Ryu's Theorems 1.2/1.4 are explicit.  Effectiveness of O(x^{0.6}) matters only for generic x. |
| `site/www/problems.html:158` | "tilts near 2/k, as in Erdős–Graham" | At k⁻ the bands of Ryu Thm 1.3 (and EG's own first level) do tilt by 2 arctan(1/b) ≈ 2/b.  Better to cite Ryu (k2-minus-c-upper, Lemma 3.2). |

## 7. Recommendation and open items

* **Formalize Ryu's Theorem 1.3 + 1.4 first** (§2.2): all closed form, k ≥ 6, roughly 1–2k Lean lines.  It is the
  same "regions glue additively" shape as Fang's lemma, so one gluing lemma serves both.  Then optionally:
  * the end fills of §3 (one lemma; c ≈ 5.66√k; beats Fang for k ≥ 38–40, which would make Fang + R1′ cover k ≥ 15
    with the better of the two at each k);
  * Theorem 1.2 (staircase), only if a 2/5 statement in Lean is wanted; it needs an interval-arithmetic constant.
* Credit: Theorem 1.3's skeleton is EG/CG Type 1; the reflected b-block is Arslanov–Mustafin–Shangitbayev's; the
  staircase is Bui's §3; the k⁻ use and the constants are Ryu's (his §1.3).  R1′'s end fill is a trivial addition;
  check with Ryu / jlevy before calling it new.
* Open: write the exact all-k count for R1′'s fills and a proved constant; run the repo's packing checker on an R1′
  instance and add R1′ packings for k ≥ 50 to the store (allowed without asking; check jlevy's register first);
  read [Q]/[T] in full if the lower bounds become load-bearing; read the Arslanov–Bui note to confirm the W₅
  diagnosis.
