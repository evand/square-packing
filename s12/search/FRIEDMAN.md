# Friedman's Conjecture 1: proof shapes, lemmas, experiments (2026-09-30)

Friedman, DS7 (1998–2009), Conjecture 1: **if `s(n²−c) = n` then `s((n+1)²−c) = n+1`.**  Evidence offered there: "true of
all the best known packings".  Literature: `tasks/friedman-lit/REPORT.md` (nothing attacks it; no result gives `c*(k) → ∞`).
Code: `search/friedman_slide.py`.  Runs: `runs/friedman_slide_v2.txt` (records), `runs/friedman_slide_all.txt` (all atlas
packings), `runs/quad_band_w{4,5,6}/` (band LP, pending).  Labels **[proved] [measured] [heuristic]**.

## 0. Summary

* Write `c*(k) := max{c : s(k²−c) = k}`.  **Friedman ⇔ `c*` is nondecreasing.**  Per c: `F_c` ⇔ `{k : s(k²−c)=k}` is
  upward closed.  Trivially **`c*(k+1) ≤ c*(k) + 2`** (§1.2) [proved].  So, if Friedman holds, `c*` climbs by 0, 1 or 2 per step.
* The three conjectures are **not linearly ordered**.  `A` = "∀c ∃k₀" ⇔ `c*(k) → ∞`.  `B` = "∀c, infinitely many k" ⇔
  `limsup c* = ∞`.  Friedman does not imply `A`: `c*` could be nondecreasing and bounded.  `A` does not imply Friedman.
  What is true: `F_c` plus one base case `s(k²−c) = k` gives `A_c`.  Base cases are known only for `c ≤ 4`.
* **Roth–Vaughan does not give `A`** [proved, given the bound as secondhand sources state it].  The bound is
  `W(α) ≥ 10⁻¹⁰⁰√(α‖α‖)`, which vanishes as the side approaches an integer from below.  QUADRANT.md §7(iii) was wrong.
  `A` is open, and so is `c*(k) → ∞`.  Known: `3 ≤ c*(k)` for `k ≥ 3`; `c*(k) ≤ k−1` for `k ≥ 12` (Arslanov et al. 2021);
  `c*(k) = O(k^{0.6})` (waste bounds).
* **Status by c.**  `F_1`, `F_2`: Nagamochi.  `F_3`: **now a theorem** (KS `k=3`, Bentz `4..7`, ours `k ≥ 6`).
  `F_4` ⇔ `s(k²−4) = k` ∀k ≥ 5, **independent of `s(12)`**; it is the k²−4 family (R = w = 3, `k ≥ 8`) plus our k = 5..8.
  `F_5` needs `D > 5/4` plus per-k certificates starting at a threshold that is itself open (`s(20)`, `s(31)`).
* **Two proof routes.**  *Dual (certificates):* the insertion lemma (§2) propagates a certificate from k to every k' ≥ k.
  This is exactly the mechanism of the k²−3 proof, and it generalises verbatim to any c.  But it proves values, not
  implications, so it gives `F_c` one c at a time, and each c needs its threshold values.  *Primal (descent):* the only
  route to Friedman in full.  Every packing of `(n+1)²−c` squares in side `< n+1` must be turned into one of `n²−c` in
  side `< n`.
* **Descent experiment** [measured, float].  Remove squares, and translate the rest in 4 classes by
  `(0,0), (−r,0), (0,−r), (−r,−r)`.  This carries out Friedman's step on **every record packing** in the atlas with
  `n ≤ 130`, except **`n = 89`** (`s = 5 + 7/√2`, a 7×7 diamond at 45°), which is short by 2.  Every non-record variant
  also works except three (§3).  Ten records have slack exactly 0.  The `n = 89` shortfall stays at 1 even with 25
  translation classes, or with diamond-frame moves (ILP optimal).  So **translation descent is not a lemma**.  A descent
  proof has to re-pack tilted blocks.  That needs a structure theory of tilted regions in near-optimal packings, which
  nobody has.

## 1. Definitions and elementary facts

### 1.1 Reformulations  [proved]
`s(k²−c) = k` is monotone in `c` (fewer squares, same side), so `{c : s(k²−c)=k} = [0, c*(k)]`.
* Friedman (all n, c) ⇔ `c*(k+1) ≥ c*(k)` ∀k.
* `F_c` ⇔ `k₁(c) := min{k : s(k²−c)=k}` exists ⇒ `s(k²−c)=k` ∀k ≥ k₁.
* `A_c` ⇔ `s(k²−c)=k` for all large k.  `A` (all c) ⇔ `c*(k) → ∞`.  `B` ⇔ `limsup c* = ∞`.
* Each `F_c` is decidable in principle **if** `A_c` holds with an effective `k₀`: finitely many `s(N) = k` questions,
  each a first-order statement over the reals (Tarski).  Without `A` there is no finiteness.

### 1.2 Ascent: `c*(k+1) ≤ c*(k)+2`  [proved]
Let `c = c*(k)`, so `k²−c−1` squares fit in side `s < k`.  We may take `s ≥ k−1`.  Add the column
`[s,s+1]×[0,s+1]`, which holds `⌊s+1⌋ = k` axis squares, and the row `[0,s]×[s,s+1]`, which holds `k−1`.  That puts
`(k+1)² − (c+3)` squares in side `s+1 < k+1`.  So `c*(k+1) ≤ c+2`.  Arslanov–Mustafin–Shangitbayev's rectangle strip
lemma is the same move.

### 1.3 Known values of c*
| k | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| lower (proved) | 2 | 3 | 3 | 4 | 4 | 4 | 4 | 3 | 3 | 3 | 3 | 3 | 3 |
| upper (records) | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 10 | 11 | 12 | 13 |

Upper row: the smallest c with a known packing of `k²−c` in side `< k`, minus 1 (atlas `site/data/packings.json`).  It is
nondecreasing, which is Friedman's "true of all best known packings".  Lower row: Nagamochi, KS, Bentz, ours
(`s(21), s(32), s(45), s(60)`; `k²−3` family).  Friedman plus our `c*(5..8) ≥ 4` would give `c*(k) ≥ 4` ∀k ≥ 5.

## 2. The dual route: insertion

### 2.1 Insertion lemma  [proved; abstracts QUADRANT.md §1.2–1.3]
Let `μ` be valid on `B = [0,a]×[0,b]`: `μ(S) ≥ 1` for every closed unit square `S ⊂ B`.  Suppose that for some `t` with
`t + 1 + √2 ≤ a`, the restriction of `μ` to the slab `[t, t+1+√2]×[0,b]` is 1-periodic in x.  That is,
`μ(A + e₁) = μ(A)` for Borel `A ⊂ [t, t+√2]×[0,b]`.  Define `μ' = μ|_{x < t+1} + τ₁(μ|_{x ≥ t})` on `[0,a+1]×[0,b]`.
Then `μ'` is valid, `μ'(B') = μ(B) + p` with `p = μ([t,t+1)×[0,b])`, and `μ'` again satisfies the hypothesis.

Proof.  Let `S` be a closed unit square with x-extent `[α, α+ω]`, `ω ≤ √2`.
* If `α+ω < t+1`: `μ'(S) = μ(S)`.
* If `α ≥ t+1`: `μ'(S) = μ(S − e₁)`.
* Otherwise `S ⊂ [t+1−√2, t+1+√2] × [0,b]`.  On that region `μ' = μ`: on `x < t+1` by definition, and on `x ≥ t+1`
  by periodicity.  Also `S ⊂ B`, since `t+1+√2 ≤ a`.

In every case `μ'(S) ≥ 1`.  ∎

Apply it vertically, then horizontally.  The x-insertion preserves y-periodicity of a horizontal slab.  If the two
inserted period masses total `≤ 2k+1`, the saving `k² − μ(B)` carries over from k to k+1, and by induction to every
k' ≥ k.  A period mass `< k` would make the saving grow linearly, which the waste bound `O(k^{0.6})` forbids.  So
equality (σ = 0) is forced.

**Consequence (certificate-level Friedman).**  An insertable certificate of saving `> c` at k proves `s(k'²−c) = k'` for
every `k' ≥ k`.  The k²−3 family is the case `c = 3`, `R = w = 2`: Lebesgue interior, σ = 0 bands, so a slab fits for
`k ≥ 2R+1+√2`, i.e. `k ≥ 7`.  `bentz_of_valid7` could be restated with `c` and the box predicate as parameters.
That is cheap, and it is the Lean form of this lemma.

### 2.2 Limits of the dual route
* It proves *values*.  For `F_c` it needs certificates for every k from the threshold `k₁(c)` up to where the
  family starts.  At `k₁(c)` the truth is barely on the grid side, so an integrality gap is likely there (s(12) is the
  c = 4 candidate).  Handling all c this way would take infinitely many thresholds, so it cannot prove Friedman.
* Fixed-profile families: `D ≤ m_v ≤ w` (seam bound, proved) and `D ≤ min(w, O(R^{0.6}))` (the box `2R+2` and waste
  `O(s^{0.6})`) [proved].  The band LP `m_v^band(w)` bounds D for every R: 0.750, 1.433, ≤ 1.755 for w = 1, 2, 3; w = 4..6
  pending (`runs/quad_band_w*`).  If `m_v^band` saturates, this family class caps at some c and cannot prove `A`.
  If it grows (the D values 0.50, 0.945, ≈1.2 look roughly `∝ √w`), an analytic family with `D(w) → ∞` would prove
  **`A`, which is open and new**.  This is the best "general" target the dual route offers.

## 3. The primal route: cut-and-slide descent

**Slide lemma** [proved, trivial].  Start from a packing in `[0,s]²` and `r > 0`.  Suppose each surviving square is
translated by one of `(0,0), (−r,0), (0,−r), (−r,−r)`, the survivors stay interior-disjoint, and they lie in
`[0,s−r]²`.  Then the survivors are a packing in side `s − r`.  A vertical and a horizontal cut, each along any
monotone curve, are the special case.  **Slide conjecture:** for `m < s < m+1` and `r = s − m + ε` this is always
possible with at most `2m+1` removals.  The slide conjecture implies Friedman.

**Experiment** (`friedman_slide.py`: exact max survivors by ILP, float SAT geometry, `EPS = 1e−5`, `TOL = 1e−9`).
An earlier run had `TOL > EPS`, so it certified side `≤ m` rather than `< m` and "found" `s(12) < 4`; that run is
discarded.

| n (record) | c | slack | | n | c | slack |
|---|---|---|---|---|---|---|
| 5, 10, 11 | 4, 6, 5 | +1, +2, +1 | | 50–55 | 14..9 | +4 +4 +2 +1 0 +1 |
| 17, 18, 19 | 8, 7, 6 | +3 +1 0 | | 65–71 | 16..10 | +4 +5 +2 +1 +2 +1 0 |
| 26–29 | 10..7 | +3 +2 0 +1 | | 82–89 | 18..11 | +4 +6 +2 +2 +2 +1 0 **−2** |
| 37–41 | 12..8 | +4 +2 +1 0 0 | | 101–110, 122–130 | | all ≥ 0; 0 at 109, 110, 130 |

* Slack = survivors − `(n − 2m − 1)`.  Records with slack exactly 0: n = 19, 28, 40, 41, 54, 71, 88, 109, 110, 130.
  The descent often does exactly what Friedman needs and no more [heuristic: the step is tight, and the descent is the
  natural mechanism].
* **Failures.**
  * Records: `n = 89` only (−2).  With a 25-vector translation grid, and with diamond-frame moves
    `(−r/2,−r/2) + {(±1/√2,0), (0,±1/√2), 0}`, the ILP optimum is still short by 1.
  * Non-records (`--all`, 232 packings): `28_r1`, `28_r1b` (−1), `50_r0r` (−3).
  * None of these threatens Friedman: records already give `s(70) < 9`, `s(17) < 5`, `s(35) < 7`.
* Reading: a straight or curved cut through a 45° block costs about a diagonal's worth of squares.  The good descents
  re-pack the block (7×7 diamond → 6×6 plus axis squares), which is not a translation.  A proof would need to show
  that every tilted region can be shrunk in its own frame, compatibly with its neighbours.  That is a structure
  theorem for near-optimal packings, beyond anything known.  (Bui's reduction to tilts ≤ 10⁻¹⁰ is asymptotic, up to
  constants, and useless at this precision.)

## 4. What a proof would look like, and what is reachable

1. **Friedman in full** needs a uniform primal descent: every packing of `(n+1)²−c` squares in side `< n+1` gives one
   of `n²−c` in side `< n`.  Translation descent is refuted (n = 89).  The candidate replacement is a
   decomposition into rigid blocks, each descended in its own frame.  **Not reachable** without a structure theory.
2. **`F_c` one c at a time** by the dual route.
   * `F_3`: done.
   * `F_4`: reachable.  It is exactly the open k²−4 family item (w = 3) plus the shipped k = 5..8, with no s(12)
     needed.  The announcement could say "Friedman's conjecture holds for c ≤ 4".
   * `F_5`: needs `D > 5/4` (w ≥ 4?) and certificates from the threshold on.  `k₁(5) ∈ {5, 6, ...}` is open (s(20)).
     A sufficient target is `s(k²−5) = k` ∀k ≥ 6, whatever s(20) is.  The shipped s(32) cover saves 4.29 at k = 6 (LP value not
     recorded), below 5.  So k = 6 is doubtful, and `F_5` may hinge on a threshold we cannot certify.
3. **`A` (`c*(k) → ∞`)** is open in the literature and not given by Roth–Vaughan.  Route: fixed-profile families with
   `D(w) → ∞`.  First check whether the band bound `m_v^band(w)` grows (runs pending).  This is the most valuable
   "general" statement within conceivable reach.
4. Small, done here: the ascent `c*(k+1) ≤ c*(k)+2`; the insertion lemma; the logical map; the Roth–Vaughan correction.

## 5. F₄ plan and compute estimate (2026-09-30; parked, not top priority)

**Order:** golf first, then the exact w = 3 run.  The cheap go/no-go comes first: does `D > 1` survive exactness?
* **Exactness price.**  At w = 2, D went 0.945 → 0.847 (layer + pitch 0.2: 0.035; tilt margin: 0.06; germ rows ~10⁻⁵).
  At w = 3, the float D is 1.154 (R = 3, pitch 0.2), so the margin is 0.154.  A w = 2-sized price leaves about 0.05.
  Run the layer + tilt-margin LP (the `qx2` pipeline at R = w = 3) first: about 5–10 CPU-h.  If it lands below 1, try
  R = 4 or w = 4 before any checker work.
* **Checker estimate** [heuristic].  w = 2 record: 81k CPU-s (box 7, 9,800 roots, ≈ 800 segments).  w = 3 needs box 9
  (roots ×1.65) with a denser corner and band (×1.5–2 per box): **≈ 50–80 CPU-h, 4–6 h wall on 16 cores, ×3 either
  way**.  Depth is the uncertain factor: a thinner margin means deeper refinement.
* **Golf targets** (cheaper per box; also inform scaling to F₅ / `c*(k) → ∞`):
  1. **Window checker.**  Check only the corner window plus one band period, not the full D4 box.  Far-band poses
     reduce mod 1; the interior is Lebesgue.  It needs a written localisation lemma (QUADRANT §1.3 / §8.1 variant).
     Probably a larger saving than the 7 → 9 growth costs.
  2. **Fewer, coarser elements.**  At w = 2, lines-only cost only 0.003 in D and pitch 0.2 was fine.  Fewer
     element types also mean simpler Lean later.
  3. **Primitives for the deep tail.**  Exact lemmas for the near-tangent germ families (as Lemma Z does for θ = 0)
     instead of subdivision, which is where depth blows up.
  4. Measure how the exactness price scales with w.  That decides whether the family route has a future beyond c = 4.

**Band LP w = 4..6** (`runs/quad_band_w{4,5,6}/`, pitch 0.1, 40-round cap, cores 12–14; started 2026-09-30).
Interim, falling, not converged: w = 4: 1.968 (rd 13), w = 5: 2.126 (rd 9), w = 6: 2.277 (rd 7); w = 3 was ≤ 1.755.
No saturation so far.  Read the final `log.txt` values into §2.2 when the runs stop.  If the bound keeps growing
roughly like √w, a genuine D(w) family (with corners) at w = 4–6 is the next measurement toward `c*(k) → ∞`.

## 6. Open items (also in TODO.md)
* F₄: go/no-go exactness LP at R = w = 3 → golf (window checker, lines-only, germ primitives) → exact run (§5).
* `c*(k) → ∞`: read the band LP results; if growing, run the corner-coupled LP at w = 4, 5 (R = w, pitch 0.2) and fit D(w).
* Parametrise `bentz_of_valid7` over c and the box predicate (the Lean form of the insertion lemma, §2.1).
* Descent: optional.  Characterise why n = 89 fails (block re-packing); not a proof route without a structure theory.
* Roth–Vaughan primary source: confirm Theorem 1's exact form if the paper is ever in hand (secondhand sources agree).
