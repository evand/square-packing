# Lower bounds on c*(k): how far covers go, what else works, sane constants (opened 2026-10-18)

Evan, 10-18: "a > constant bound with explicit and vaguely sane constants?"  Questions: how far covers get; other
techniques; whether Ryu's constants can be improved a lot; a provable limit on point covers; other cover geometries
(as the line densities were); questions we haven't asked yet.  This is a framing doc: nothing below is new
mathematics.  Labels: **[proved] [measured] [heuristic] [reported]** (reported = someone else's unrefereed claim).
Background: `search/FRIEDMAN.md` (§2 insertion, §8 seam capacity, §9 ladder), `search/K2M4_MARGIN.md`, `search/QUADRANT.md`,
`search/SEAM_W1.md`, `search/SEAM_1D.md`, `search/CORNER_DEFICIT.md`, `search/S2_INSERTABLE.md`, `search/CLIQUE*.md`,
`tasks/issue3-k2-minus-c/review-2026-10-10.md`, `literature/waste-constructions-2026-10-17.md` §4.

`c*(k) = max{c : s(k² − c) = k} = k² − M(k) − 1`, where M(k) is the maximum number of pairwise disjoint (as closed sets)
unit squares in [0,k]².

## 0. Where we stand

| result | c covered | from k | status |
|---|---|---|---|
| ours: covers (k²−3 family R = w = 2; k²−4 family R = w = 3; per-k covers s(21), s(32), s(45), s(60)) | ≤ 4 | 5 (c = 4), 3 (c = 3) | [proved], computer-assisted; k ≥ 8 for c = 4 rests on one implementation (`ValidTilt9`) |
| Ryu [R], log k: W ≥ 0.0353 log k | all | log k > 13.07 + 15.21c (c = 5: k > 10^38.7) | [reported]; our referee read found no gap |
| Ryu [T], √k: W > 10⁻⁴√k on [10⁸, 2.5·10¹⁵], then 7·10⁻⁵√k, then 0.06 k^{1/3} beyond 3.5·10¹⁷ | all (√k range: c ≤ 5000) | (10⁴c)² (c = 5: 2.5·10⁹) | [reported], 9 pp., not reviewed by us |
| Ryu [Q], k^{1/4}: W ≥ 0.1 k^{1/4} | all | max(4.62·10¹², 10⁴c⁴ + 1) | [reported], 108 pp., not reviewed |
| upper bounds: c*(k) ≤ 4⌈√(2k−7)⌉ − 1; < 18.9 k^{2/5} + … ; Ryu k^{3/8} | — | — | Lean (ours) / [reported] |

So, as reported: **k^{1/3} ≲ c*(k) ≲ k^{3/8}**.  The "huge constants" problem is the thresholds: c = 5 is open for
7 ≤ k < 2.5·10⁹.  "Sane" target: s(k² − c) = k for all k ≥ k₀(c) with k₀ a modest polynomial (say ≤ 10c³), or at least
c = 5, 6 from small k.  Stale: `FRIEDMAN.md` §0 and §9.0 still quote O(k^{0.6}) and "threshold ≳ c^{5/3}"; per the lit note it is ≳ c^{8/3}.

## 1. How far can covers get?

The cover method: a measure μ ≥ 0 with μ(S) ≥ 1 for every closed unit square S ⊂ [0,k]²; then M(k) ≤ μ([0,k]²).
Best saving L(k) := k² − inf μ.  `L(k) ≤ c*(k) + 1` [proved].

* **Per c, with a fixed family** (corners + periodic σ = 0 wall bands + Lebesgue interior).  Saving 4D(w), valid from
  k = 2R + 2.  D(w) = 0.50, 0.945, 1.154, ≈1.26, ≈1.30 for w = 1..5 [measured, float, pitch 0.2/0.1]; the increments
  halve or worse.
  * c = 5 (needs D > 5/4) at R = w = 5: margin ≈ 0.05, κ ≈ 0.02; conditional GO (`K2M4_MARGIN.md` §3).
  * Exact-run cost so far, for scale: c = 3 81k CPU-s, c = 4 815k CPU-s (×10).  c = 5 guess ~10⁷ CPU-s (≈3000 CPU-h)
    unless golfed first (window checker, near-tangent lemmas; `FRIEDMAN.md` §5).  [heuristic]
  * Base cases for "c*(k) ≥ 5 ∀k ≥ 7": per-k covers at k = 7..11.  L(7) ≈ 5.21 and L(8) ≈ 5.60 [measured, sampled rows]
    look feasible.  k = 6 (s(31) = 6) looks like an integrality-gap case (L(6) ≈ 4.7?, not measured), like s(12) for c = 4.
  * c = 6 needs D > 1.5.  The current trend doesn't reach it at any w we can afford, unless finer pitch moves D.
* **All c, with an analytic family** (D(w) → ∞): R2 (seam capacity M(w) → ∞) ∧ R3 (D ≥ κM) ⇒ c*(k) ≳ k^β
  (`FRIEDMAN.md` §9.2).  The heuristic exponent is β = 1/3: budget ½ per period, seam price ≈ 0.1 g^{1.5}.  It matches
  [T]'s exponent, which is suggestive.  Thresholds would be k₀ ~ C·c³ with C from the 1D lemma constants: that is the sane version.
* **Decisive cheap measurements**, before any c = 5 exact run:
  1. rigorous caps on M(w) at w = 2, 3 (dual-side band LP, `SEAM_W1.md` §4; arrangement-exact, so cheap to certify);
  2. band and corner LP at pitch 0.05 / 0.02, w = 4..8: does D(w) keep growing, or saturate near 1.3–1.4?
  3. the 1D seam problem **with a height phase** (`FRIEDMAN.md` §8 correction): its price at a few g.  This is the crux of R2.

## 2. Is there a provable limit to point covers?

By weak duality, any **fractional closed packing** (weights y_S ≥ 0 with Σ_{S ∋ p} y_S ≤ 1 for every point p) of
weight F gives L(k) ≤ k² − F.  So limits are fractional packings, which we already build (`CLIQUE_CEILING.md`).
Known limits:

* L(4) ≤ 4: no point cover proves s(12) = 4 [measured, fractional packing]; box cliques gain ≈ 0.1 there (`CLIQUE_CONTINUUM.md`).
* Corner-only modules save nothing: D*(R) = 0 (`CORNER_DEFICIT.md`) [proved].
* Family class: D ≤ m_v ≤ M(w) (seam bound) [proved]; seam-avoiding duals cap at w − ¼ [proved, but only an upper bound on what those duals can prove]; M(1) = ¾.
* Insertable at the threshold costs ≈ 1.7–1.9 (S2 no-go, k = 4, 5) [measured].
* Trivially L(k) ≤ c*(k) + 1 ≲ k^{3/8}, so no cover beats Ryu's upper-bound exponent.

Open, and the real question here: **is L(k) → ∞?  Is L(k) − c*(k) bounded?**
* Ryu's line count (Lemma 3.2: at most k − 1 chords of length ≥ 1 per line) **is LP-compatible**.  The 1D cover with
  atoms at 1, …, k − 1 proves it fractionally, and that is exactly what our line densities are.
* His Roth–Vaughan step charges *pairs* of neighbouring squares by relative tilt.  It is not obviously true for
  fractional packings.  Either:
  * a fractional packing with weight k² − O(1) exists (then point covers are capped, and growth needs pair or clique lifts, §4); or
  * the RV lemma has a fractional form (then an explicit cover with saving ~ log k, or k^{1/3}, exists; what does it look like?).
* First step: build fractional (clique) packings for k = 4..12 and fit the growth of k² − F(k) (log k vs k^{1/3} vs bounded).
  Also look for fractional packings that beat the integral ones at the same k.

## 3. Can Ryu's constants be improved significantly?

From our [R] read (constants reproduced): slope 0.0657 per unit of log k, bracket 30.42.  The losses:
* A ≈ 10.6 = 4√(K_w/Q*)/(2 − ε), with K_w = 9 sharp *in his model* and Q* = 0.32 (truth plausibly 1/3: only 2 %);
* additive 13.07, mostly log(y_min) with y_min = 236000 (every line within y_min of a wall is discarded);
* tilt and rigidity thresholds of about 10⁻⁵ (Lemma 5.2).
[T] and [Q] reuse [R]'s Lemmas 4.10, 4.13 and 4.15.

* Polishing gains perhaps 2–5× in the slope.  For the log bound the thresholds stay e^{Θ(c)}, so polishing is useless there.
* For [T], with threshold (c/κ)², a 10× better κ means a 100× smaller threshold: 2.5·10⁷ for c = 5.  Still not sane.
  Reaching k ~ 10³ for c = 5 needs κ about 1500× better.  That is a different method, not polishing.
* **Hybrid** (most promising): [R]/[T] discard the near-wall region, d < y_min.  Our covers are exactly a near-wall tool.
  Question: does his charging localise, so that a certified wall-band saving (our band LPs) can replace the discarded
  strip and the asymptotic part only handles d ≥ w?  Needs a read of [T] (9 pp.) for where y_min and the thresholds enter.
* Prerequisites: review [T] (and [Q] §§3–8 as used); the clean-room Lemma 4.10 re-check (`review-2026-10-10.md` §6, ~1.7 h compute).

## 4. Other techniques, and other cover geometries

What has worked: point covers → line and segment densities (validity overhead 1.4–2.9 % → 0.04–0.15 %, `LINE_COVER.md`) →
box cliques (`BOXCLIQUE.md`, ≈ 0.1 at t = 4) → families and insertion; casework (Nagamochi c ≤ 2); RV / Ryu line waste
(asymptotic); SOS (small cores only, no zero margin, `SOS_PROBE.md`).  Candidates:

* **Tilted line densities** at the band angles of near-optimal packings (2 arctan(1/b), 45°): line mass where the extremal
  fractional packings put tilted rows.  Cheap to add to `zm_mixed` as a column class.
* **Off-integer line families**: Ryu's quantization says ramps of rigid columns sit near integer distances *from each wall*.
  A cover with lines at heights j ± δ_j, or weighted by 1/d(y), is the guessed dual of his dy/d(y) integration.  Does a
  cover of that shape have growing saving?
* **Cliques in the wall band**: does the seam bound D ≤ m_v survive for clique certificates?  If not, cliques could
  beat the family ceiling.
* **Pair-lifted covers** (Lasserre level 2, two-point bounds in the spirit of Cohn–Elkies): weights on pairs of squares,
  Σ f(S_i) − Σ g(S_i, S_j) ≤ B over packings.  The natural dual of RV's pair charging.  Nothing built; finite-k exactness is unclear.
* Arc / circle measures (rotation-adapted), smooth hats: hats tested, no gain (`QX2_PL.md`); arcs untried, low prior.
* Primal descent (Friedman): refuted for translations (n = 89); needs a structure theory.  Not this task.

## 5. Questions not yet asked

1. Generic sides: do covers give waste lower bounds just *above* an integer (s(k² + 1) − k, s(k² + j)), where RV says
   nothing (x‖x‖ < 1/6)?  The machinery is the same; a new statement class.
2. The first-failure catalogue: for which (k, c) is L(k) ≤ c < c*(k) + 1 (integrality gaps: s(12) for c = 4, s(31)? for c = 5)?
3. Fractional Friedman (C2: L nondecreasing) checked numerically on the same L(k) table as §2.
4. Rectangles c*(a, b) (C8): insertion in one direction at a time; maybe c = 5 is cheaper in thin rectangles first.
5. Does the cover-provable exponent (sup β with L(k) ≳ k^β) equal the true one?  Can covers or fractional packings
   decide 1/3 vs 3/8?
6. The corner: why does the w = 1 corner cost exactly ¼?  Is D/M → constant (R3)?  ≈ 0.61–0.67 so far.
7. Uniform-in-c Lean: `bentz_of_valid` parametrised by c and the box predicate (insertion lemma in Lean), so any new
   family plugs in.
8. Which of Ryu's ingredients are LP duals of something we have: line count = line densities; quantization = ?; RV = pairs?
   A dictionary would show which parts can be computer-certified at small k.

## Next (ordered, nothing started)
§1 items 1–2 (M(w) caps; finer-pitch D(w)), §2 fractional-packing table k = 4..12, read [T], then §1 item 3.  Compute needs Evan's OK.
