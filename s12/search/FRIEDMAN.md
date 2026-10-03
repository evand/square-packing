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
  `F_4` ⇔ `s(k²−4) = k` ∀k ≥ 5, **independent of `s(12)`**: **now holds** (2026-10-03, `certificates/k2m4/`): the k²−4 family (R = w = 3, `k ≥ 8`, resting on `ValidTilt9`, one implementation) plus our k = 5..8.
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

* **Conjecture ladder (§9, 2026-10-02).**  Thresholds `k₁ ≤ k_LP = k_D4 ≤ k_ins ≤ k_fam` by certificate class; group
  symmetry is free, insertability is the first paid restriction and its price is seam capacity.  Rate side R1–R3
  (families prove A iff seam capacity grows and corners cost a bounded fraction); structure side S1–S5 (bounded splice
  cost ⇒ almost-monotone fractional Friedman).  Any threshold is ≳ c^{5/3}.  First test: S2 at c = 4, k = 5.

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
  `O(s^{0.6})`) [proved].  The seam capacity `M(w)` bounds D for every R.  `M(1) = ¾` exactly (SEAM_W1.md);
  for w ≥ 2 the band LP values 1.433, ≈1.755, 1.965, 2.116, 2.231 (w = 2..6) are **lower bounds** on `M(w)` (columns
  restricted to the pitch-0.1 lattice; correction 2026-10-02, SEAM_W1 §4.1), and the only proved cap is `w − ¼`.  If `m_v^band` saturates, this family class caps at some c and cannot prove `A`.
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
   * `F_4`: done (2026-10-03, `certificates/k2m4/`): the k²−4 family (w = 3, `k ≥ 8`, resting on `ValidTilt9`,
     certified by one implementation) plus the shipped k = 5..8, with no s(12) needed.
   * `F_5`: needs `D > 5/4` (w ≥ 4?) and certificates from the threshold on.  `k₁(5) ∈ {5, 6, ...}` is open (s(20)).
     A sufficient target is `s(k²−5) = k` ∀k ≥ 6, whatever s(20) is.  The shipped s(32) cover saves 4.29 at k = 6 (LP value not
     recorded), below 5.  So k = 6 is doubtful, and `F_5` may hinge on a threshold we cannot certify.
3. **`A` (`c*(k) → ∞`)** is open in the literature and not given by Roth–Vaughan.  Route: fixed-profile families with
   `D(w) → ∞`.  First check whether the band bound `m_v^band(w)` grows (runs pending).  This is the most valuable
   "general" statement within conceivable reach.
4. Small, done here: the ascent `c*(k+1) ≤ c*(k)+2`; the insertion lemma; the logical map; the Roth–Vaughan correction.

## 5. F₄ plan and compute estimate (2026-09-30; done 2026-10-03: `certificates/k2m4/`, k ≥ 8 rests on `ValidTilt9`)

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
* F₄: done (2026-10-03, `certificates/k2m4/`; k ≥ 8 rests on `ValidTilt9`, one implementation).  Open: a second implementation of `Valid9`.
* `c*(k) → ∞`: read the band LP results; if growing, run the corner-coupled LP at w = 4, 5 (R = w, pitch 0.2) and fit D(w).
* Parametrise `bentz_of_valid7` over c and the box predicate (the Lean form of the insertion lemma, §2.1).
* Descent: optional.  Characterise why n = 89 fails (block re-packing); not a proof route without a structure theory.
* Roth–Vaughan primary source: confirm Theorem 1's exact form if the paper is ever in hand (secondhand sources agree).

## 7. Conjectures and weaker variants (2026-10-01)

Notation: `L(k) := k² − (value of the exact continuum closed cover LP on [0,k]²)`, the best saving any LP certificate can
show.  `L(k) > c ⇒ s(k²−c) = k`, and `L(k) ≤ c*(k) + 1`.  `S_ins(k)` := best saving over *insertable* certificates at k
(§2.1: one 1-periodic slab of width `1+√2` per direction, period masses summing to `2k+1`).  `M(w) := m_v^band(w)`
(seam capacity of a σ = 0 band of width w, QUADRANT §4); `M(∞) := sup_w M(w)`.

| # | statement | status | kind of attack |
|---|---|---|---|
| C1 | Friedman: `c*(k+1) ≥ c*(k)` | open; translation descent refuted (n = 89, §3) | primal structure theory |
| C2 | **Fractional Friedman**: `L(k+1) ≥ L(k)` (every LP certificate propagates upward) | open; data `L(4) ≤ 4, L(5) ≈ 4.25, L(7) ≈ 5.21, L(8) ≈ 5.60` (sampled-row LP, S60_COVER §4) consistent | fractional descent (§7.1) |
| C3 | `S_ins` is nondecreasing | **[proved]** (insertion lemma: an insertable cert at k gives one at k+1 with the same saving) | — |
| C4 | Insertable class is asymptotically complete: `L(k) − S_ins(k) → 0` (or stays bounded) | open; `lim S_ins = 4 sup D ≤ 4 M(∞)` | compare LPs |
| C5 | **Seam capacity**: is `M(∞) < ∞`?  If yes, every fixed-profile/insertable family has `D ≤ M(∞)`, so the dual route caps at `c < 4M(∞)` (≈ 10–12 if M(∞) ≈ 2.5–3) | open; `M(1) = ¾` (SEAM_W1); measured **lower bounds** `M` ≥ 1.43, ≈1.76, 1.97, 2.12, 2.23 (w = 2..6, pitch-0.1 columns; not caps, SEAM_W1 §4.1); §8, §9 | analytic dual (§8) |
| C6 | `A_c` with explicit threshold: `s(k²−c) = k` ∀k ≥ k₀(c) (no base case at the true threshold k₁(c) needed) | c ≤ 4 proved (c = 4: R = w = 3, k₀ = 8, `certificates/k2m4/`, 2026-10-03); c = 5 open | `S_ins(8 or 9) > 5` would give c = 5 from one box |
| C7 | `c*(k) → ∞` (A) | open in the literature | `L(k) → ∞` suffices; via C4+C5 only if `M(∞) = ∞` |
| C8 | Rectangle Friedman: `c*(a,b) := max{c : ab − c unit squares need an a×b box}` is nondecreasing in each of a, b separately | open; half-step of C1; each step is one insertion (one periodic slab) | as C1–C3, one direction at a time |

Finer ladder between these (rate conjectures R1–R3, structure conjectures S1–S5, thresholds by certificate class): §9.

Logical map: C3 is the certificate-level Friedman for the insertable class.  C2 ⇒ (any LP cert at k gives `A_c` from k).
C4 ∧ ¬C5-finite ⇒ C7.  If C5 holds (`M(∞) < ∞`) and `L(k)` grows, then C4 fails: the growth of `L` must come from
aperiodic walls, and the families route cannot prove C7.

### 7.1 Fractional descent (C2): the four-sub-box average  [proved as stated; not a proof of C2]
Let `P` be a fractional strict packing of `[0,k+1]²` (k ≥ 3).  Restrict it to the four sub-boxes `[i,i+k]×[j,j+k]`,
`i,j ∈ {0,1}`, translate each to `[0,k]²` and average: a fractional packing of `[0,k]²` of weight
`W(P) − Loss`, `Loss = ½(W_L + W_R + W_B + W_T) − ¼(W_BL + W_BR + W_TL + W_TR)`, where `W_L` = weight of squares whose
x-extent starts below 1, and so on (no square is near two opposite walls when k ≥ 3).  Hence
`L(k+1) ≥ L(k) + (2k+1) − Loss(P*)` for an optimal `P*` at k+1.  For the grid, `Loss = 2k+1` exactly (tight).  It
does not prove C2 unconditionally: wall layers can have weight ≈ √2 per unit length (alternate axis squares at x ∈ [0,1]
with 45° squares whose left vertex is at x ≈ 0.99).  So C2 needs control of the wall layer of *optimal* packings.  General shift laws π on [0,1]² give the same shape: the
loss is `E_{(u,v)~π}[weight of squares not inside [u,u+k]×[v,v+k]]`, i.e. squares cut by the boundary lines x = u,
x = u+k, y = v, y = v+k.  So **C2 ⇐ some optimal fractional packing at k+1 has a boundary frame (on average over a shift law)
cut by weight ≤ 2k+1**.  The dual (cover) side of the same move: smooth the optimal cover by π and repair near the walls;
the repair is a wall-layer cover problem whose naive price (≥ √2 per unit wall) is too high, so the repair has to use
the smoothed μ's own partial mass near the walls.  Both sides point at the same missing input: wall structure of optima.

## 8. Seam capacity (C5)

Setting (band LP of QUADRANT §4): `μ ≥ 0` on the half-plane `H = {y ≥ 0}`, 1-periodic in x, Lebesgue on `y > w`,
`μ([0,1)×[0,w]) = w` (σ = 0), valid (`μ(S) ≥ 1` for every closed unit square `S ⊂ H`).  `M(w) = sup μ({0}×[0,w])`.
Per-period height profile `m(A) := μ([0,1)×A)`; excess `ν̄ := m − Leb`, a signed measure on `[0,w]`, total 0, `ν̄ ≥ −Leb`.

**Budget lemma** [proved].  (i) `m(I) ≥ 1` for every closed unit interval `I ⊂ [0,∞)`.  (ii)
`∫_{1/2}^{∞} (m([c−½, c+½]) − 1) dc ≤ ½`.
Proof.  (i) Average the axis squares `[t,t+1]×I` over `t ∈ [0,1)`: each point of `[0,1)×I` lies in `[t,t+1]` (mod 1)
for a t-set of measure 1, so the average is `m(I)`, and every term is ≥ 1.  (ii) The integral is
`∫ ν̄(dy)·|{c ≥ ½ : |y − c| ≤ ½}| = ∫ ν̄(dy) min(y,1) = ∫ν̄ − ∫_{[0,1]}(1−y) ν̄(dy) ≤ 0 + ∫_0^1 (1−y) dy = ½`. ∎
Equality needs `m([0,1)) = 0` with `m([0,1]) ≥ 1`, i.e. the whole first row's mass on the line `y = 1` (Nagamochi's line).
The same holds with tilted squares at angle θ (window = the chord-length kernel of height `h_θ = cos θ + sin θ`), with
bound `h_θ/2`.

**Measured** (`search/seam_budget.py` on `runs/quad_band_w{1,4,6}/support.txt`): budget used 0.31 / 0.40 / 0.41 of ½.
Per row j, (seam mass in `[j,j+1)`, excess integral over `c ∈ [j+½, j+1½)`): w = 6: (0.51, 0.19), (0.51, 0.087),
(0.40, 0.058), (0.34, 0.033), (0.28, 0.024), (0.18, 0.013).  The seam carried per unit of excess *grows with height*
(2.7 → 14).  So the naive bound "each unit of seam costs ≥ κ of excess" (true for y-invariant profiles, κ ≈ 0.17 from the
45° square between seams) is false for the optimum, and `M(∞) < ∞` is not implied by the budget lemma alone.
`M(∞) < ∞` ⇔ the global ratio seam / excess is bounded (budget ≤ ½ is fixed).  Next: read the LP duals
(`search/seam_dual.py`, `runs/seam_dual/`): they are fractional packings certifying `M(w)`; an analytic pattern in them
valid for all w would prove C5, and a pattern whose cost decays with height would point to `M(∞) = ∞`.

**Bulk price of seam: sublinear** [measured, float; `search/seam_1d.py`].  In the bulk a y-invariant profile
`ρ(x) × Leb_y` is *exactly* a 1D problem: a square (θ, centre offset cx) gets `∫ρ(x) V_θ(x − cx) dx`, `V_θ` = its
vertical chord at horizontal offset d.  1D LP: atoms on `(1/N)ℤ`, seam atom g at 0, total `1 + e`, all θ on a 0.5° grid
(+ tiny tilts), cx on a `1/2N` grid shifted off the atoms (no ε-gap freebies):

| e | 0.0025 | 0.005 | 0.01 | 0.02 | 0.05 | 0.1 | 0.2 |
|---|---|---|---|---|---|---|---|
| g (N = 200) | 0.0808 | 0.1273 | 0.2069 | 0.3003 | 0.4506 | 0.700 | 1.177 |
| κ = e/g | 0.031 | 0.039 | 0.048 | 0.067 | 0.111 | 0.143 | 0.170 |

N = 400 reproduces every g to 4 digits (atom resolution is not binding).  Smaller e (gridded LP, N = 400, float): e = 0.001, 0.0005, 0.00025 → g = 0.0475, 0.0307, 0.0204 (κ = 0.021, 0.016, 0.012; local exponent of e in g ≈ 1.6–1.7); the continuum values sit 3–7 % lower in g (SEAM_1D.md).  **Exact certificates** (SEAM_1D.md, `seam_1d_exact.py`): κ = 0.0221, 0.0314, 0.0397, 0.0486 at e ≈ 0.001, 0.0025, 0.005, 0.01.  Dense re-check of the e = 0.005 solution (θ every 0.01°, cx dense + every atom-at-chord-end ε-pose): min mass 0.99975
(θ = 19.83°), so `g = 0.127` at `e ≈ 0.00525` after a uniform top-up.  Fit: **`e ≈ 0.1·g^{1.5}`** for small g (local
exponent 1.4–1.9).  The naive construction (condense a width-g Lebesgue strip onto the seam line) has price ≈ g/8
(linear: loss at the convex kink of `V_θ` at the far vertex, worst at tilt ≈ g); the LP beats it with mass bunched at
mid-phase `[¼, ¾]` plus off-seam pieces, the same structure the w = 6 band optimum shows in its upper rows.

**Analytic witness: the condensed seam** [measured, float; `search/seam_condense.py`].  ρ = Lebesgue off the strips
`(−g/2, g/2) + ℤ`, plus the seam line of density g (the strip's mass condensed onto x = 0), plus a uniform top-up e.
Needed `e = 0.1421, 0.1378, 0.1367, 0.1363 × g^{1.5}` for g = 0.1, 0.03, 0.01, 0.003, worst tilt `θ ≈ 0.82 √g`.  Why
g^{1.5} [heuristic, matches]: the loss is at the far vertex of a tilted square (convex kink of `V_θ`, slope `1/(sc)`), but
the *next* seam image sits at the opposite plateau corner (concave kink) only `δ = 1 − cos θ ≈ θ²/2` away, so the two
nearly cancel: loss ≈ `gθ/4` for `θ² ≲ g`, ≈ `g²/(8θ)` for `θ² ≳ g`; worst at `θ ~ √g`.  (A first estimate that ignored
the pairing gave a linear price g/8; wrong.)  The LP improves the constant 0.136 → ≈ 0.1.

**Row lemma** [proved].  Let `X(t) := m([0,t]) − t` (per-period running deficit; σ = 0 means `X(w) = 0`, and `X` is
constant above w).  Validity of the axis squares `[s,s+1]×[t,t+1]` (averaged over s) gives `m([t,t+1]) ≥ 1`, i.e.
`X(t+1) ≥ X(t⁻)` for `t ≥ 0`.  Hence along each residue class `r + ℕ` (r ∈ [0,1)) the sequence `X(r+n)` is
nondecreasing except where the class carries y-atoms (horizontal line mass), `X(r) ≥ −r` (m ≥ 0), and `X(r+n) ≤ 0`
for every class without y-atoms (chain up to the Lebesgue region).  So: (i) **the total rise along class r is ≤ r**
(average over r: the budget lemma's ½ again); (ii) if no line mass sits at integer heights, **every integer-aligned row
`(n, n+1]` has mass exactly 1**; (iii) a bulk that is 1-periodic in y has every half-open unit row of mass 1 in every
class, so every square's mass averages to exactly 1 over positions, hence (validity + Fourier rigidity) it is Lebesgue:
**no seam without drift**.  The w = 6 optimum shows exactly this: `X ≤ −0.0002` everywhere, `X(n) = −0.001…0` at
integers, and along class r the values climb monotonically from ≈ −r (empty floor) to 0 (class 0.05: −0.051 → −0.001;
class 0.8: −0.619 → −0.001).

**Correction to the heuristic below.**  The y-invariant bulk with uniform excess e is *not* realisable in a band: by
(i)–(ii) the excess seen by a square depends on its height phase, and windows centred at half-integer heights see
none.  What survives: the excess is a drift, per class ≤ r, total ≤ ½, and the question is whether seam can sit where
squares see excess.  Small-tilt squares centred at height c miss seam only near heights `c ± ½` (≈ s/4 each), so seam
placed away from integer heights is invisible to the zero-excess (half-integer-centred) squares; the w = 6 optimum does
thin the seam just above integers (0.06 vs 0.16–0.19 per 0.1 slab, rows 1–4).  Whether the price stays sublinear under
these constraints is open; the next model is the 1D problem with a height phase (x-periodic, y-structure within rows,
slow drift).

**Consequence (heuristic, superseded in part, see the correction above).**  The excess is a fixed budget (≤ ½, budget lemma) and per-row seam g costs `≈ 0.1 g^{1.5}`
per unit height, so spreading it as g ∝ w^{−2/3} over w rows gives **`M(w) ≳ C·w^{1/3}` → ∞**.  Linear prices
(κ > 0) would have capped M; a superlinear price doesn't.  The band LP at pitch 0.1 can't resolve the fine structure
(the 1D optimum at small e uses pitch ≤ 0.005), which is consistent with its values (2.23 at w = 6) sitting below the
prediction (≈ 4).  Missing for a proof of `M(∞) = ∞`: (a) an exact 1D certificate at some small e (1D, small; exact
check over continuous θ is the only real work) — **done** [proved, exact B&B; `search/SEAM_1D.md`]: κ ≤ 0.0221 at
g = 0.0459, κ ≤ 0.0315 at g = 0.0799, κ ≤ 0.0397 at g = 0.1265, κ ≤ 0.0486 at g = 0.2064; (b) a valid wall row that spends the budget (≤ ½) while feeding the
bulk; (c) slow-variation lemma: a profile with amplitude `g(y)` varying on scale ≫ 1 costs `O(|g'|)` extra excess
(first-order term, total O(max g)); (d) the top transition to Lebesgue.  None looks hard.  It would **not** by itself
give `D → ∞`: the corner coupling (QUADRANT §6: the corner is the binding constraint) is the next question, and
`D(w) ≤ O(w^{0.6})` anyway (box saving ≤ c*(k)+1).

## 9. Conjecture ladder: between the large-constant asymptotics and Friedman (2026-10-02)

Purpose: name the statements strictly between "c*(k) = O(k^{0.6}) with unusable constants" and C1, chosen so that our
family / certificate machinery can test them and so that lemmas proved on the way are reusable.  Seam capacity
results: SEAM_W1.md.  Labels as above.

### 9.0 The constraint every rate conjecture must respect
Packing constructions with waste `O(x^α)` for every real x (α = 0.6 as cited in §0; **re-check the source and
whether the constant is effective** before quoting) give `c*(k) ≤ C k^{0.6}`, i.e. **s(k²−c) = k forces
k ≳ (c/C)^{5/3}**.  So no threshold can be linear in c: "k₁(c) ≤ a·c + b" is false for large c.  The only explicit
small-k statement known to us: `c*(k) ≤ k − 1` for k ≥ 12 (Arslanov et al. 2021, secondhand), i.e. k ≥ c + 1.
Roth–Vaughan gives nothing here (§0).  Every conjecture below is either polynomial with exponent ≥ 5/3 or structural.

### 9.1 Threshold hierarchy by certificate class
For a class X of certificates let **`k_X(c)` := least k at which an X-certificate proves s(k²−c) = k**.

| class X | restriction on the certificate | X upward closed? | threshold |
|---|---|---|---|
| truth | none | C1 (Friedman) | `k₁(c)` |
| LP | valid closed cover measure on `[0,k]²` | C2 (open) | `k_LP ≥ k₁` |
| D4-LP | D4-invariant | as LP | **`= k_LP`** (9.1.1) |
| insertable | one 1-periodic slab of width `1+√2` per direction, period masses totalling `2k+1` (§2.1) | **yes, C3** | `k_ins ≥ k_LP` |
| family | corner modules + σ = 0 periodic walls + Lebesgue interior (QUADRANT) | yes | `k_fam ≤ 2R(c)+2`, `≥ k_ins` |
| analytic wall | e.g. x-uniform + seam (SEAM_W1 §5) | yes | `≥ k_fam` |

**9.1.1 Averaging lemma [proved, trivial].**  If G is a finite group of isometries of `[0,k]²` and μ is valid, then
`μ̄ = |G|⁻¹ Σ g_*μ` is valid (`μ̄(S) = avg μ(g⁻¹S)` and each `g⁻¹S` is a closed unit square in the box) with the same
total mass.  So *group* symmetry of the box is free.  Insertability is a translation symmetry of a strip, not a
symmetry of the box; it is the first restriction that costs anything.  **Lean:** `exists_d4InvM_cover` (`lean/Sqpack/Average.lean`, generic `coverValid_avg` for any finite family of square-good maps).  Its cost is paid where a slab crosses a wall
(the interior is Lebesgue and splices for free): there the wall must continue as a σ = 0 periodic band, so **seam
capacity is the price of insertability**.

**9.1.2 The gaps mean different things.**  `k_LP − k₁`: integrality gap (no measure proves the truth; s(12) is the
c = 4 candidate, `notes/n12-gap.md`).  `k_ins − k_LP`: price of extensibility (C4 at the threshold).
`k_fam − k_ins`: price of periodic *whole* walls, mostly corners (9.3).

**9.1.3 Status for c = 4.**  `k₁ ∈ {4, 5}` (s(12) open); `k_LP ≤ 5` (s(21) cover); `k_fam ≤ 8` (R = w = 3, `certificates/k2m4/`,
resting on `ValidTilt9`); **`k_ins ∈ [5, 8]` unknown**.  For c = 3: `k₁ = 3`, `k_fam = 6` (proved, k2m3 bundle).

### 9.2 Rate ladder (the A side)
| # | statement | relations | handle |
|---|---|---|---|
| R1 | **polynomial growth** `c*(k) ≥ k^β`, some β > 0 | ⇒ A (C7); β ≤ 0.6 forced (9.0) | follows from `D(w) ≳ w^β` (a family of width w lives in boxes ≈ 2w+2) |
| R2 | **polynomial seam capacity** `M(w) ≳ w^β` (heuristic β = ⅓, §8) | necessary for R1 via families; ⇒ `M(∞) = ∞` (¬C5-finite) | 1D exact certs exist (SEAM_1D); need wall-row + slow-variation lemmas (§8 (b)–(d)); rigorous caps via dual-side LP (SEAM_W1 §4) |
| R3 | **corner costs a bounded fraction** `D(w) ≥ κ M(w)` | R2 ∧ R3 ⇒ R1 | corner duals; start with the exact ¼ at w = 1 |

R2 ∧ R3 is the cleanest "families prove A" programme: one band question, one corner question, each local.

### 9.3 Data behind R3 and the family delay  [measured; M for w ≥ 2 are lower bounds]
| w | M(w) band | D(w) with corner | D/M | largest c with 4D > c | family start 2w+2 | first LP instance of that c |
|---|---|---|---|---|---|---|
| 1 | ¾ (exact) | 0.500 | 0.67 | 1 | 4 | — |
| 2 | ≥ 1.433 | 0.945 | 0.66 | 3 | 6 | k = 3 (truth) |
| 3 | ≈ 1.755 | 1.154 | 0.66 | 4 | 8 | k = 5 (s(21)) |
| 4 | ≥ 1.965 | ≈ 1.26 (K2M4_MARGIN k²−5 caps) | ≈ 0.64 | 5 | 10 | k = 8 proved (s(59), wand125); L(7) ≈ 5.2 |
| 5 | ≥ 2.116 | ≈ 1.30 (κ = 0.02) | ≈ 0.61 | 5 | 12 | |
| 6 | ≥ 2.231 | — | — | (5–6 if D/M ≈ 0.62) | 14 | |

Family delay `k_fam / k_LP` ≈ 2, 1.6, 1.3–1.7 for c = 3, 4, 5 [heuristic].  If M saturates near 2.5–3 the route
stalls near c ≈ 6–7; if R2 holds with β = ⅓ then `k_fam ~ c³` against the floor `c^{5/3}`: polynomially late but
every c reachable.  The apparent saturation may be the pitch-0.1 grid (SEAM_1D: small-amplitude seam structure
needs pitch ≤ 0.005); undecided.

### 9.4 Structure ladder (the Friedman side)
| # | statement | relations | handle |
|---|---|---|---|
| S1 | **bounded splice cost**: ∃ δ < ∞ ∀ k ≥ k₀ ∀ valid μ on `[0,k]²` ∃ insertable valid μ′ with `μ′([0,k]²) ≤ μ([0,k]²) + 2δ` | ⇒ **`L(k′) ≥ L(k) − 2δ` ∀ k′ ≥ k** (C3) — almost-monotone fractional Friedman | cost is local: two wall crossings per direction, each a "splice a periodic band segment into an arbitrary wall" LP; families = whole wall periodic; seam capacity bounds what the patch carries |
| S2 | **free splice at thresholds**: `k_ins(c) = k_LP(c)` | S1 with δ = 0 at the threshold; ⇒ every LP-provable value starts a family | **testable now**: insertable box LP at k = 5 for c = 4 |
| S3 | fractional Friedman (C2) `L(k+1) ≥ L(k)` | ⇐ S1 with δ = 0 for all k | §7.1 four-sub-box averaging (needs wall-layer control) |
| S4 | **bounded family delay** `k_fam(c) ≤ 2 k_LP(c)` | quantitative isolated → family bridge; false if M saturates (9.3) | sharpen with rigorous `M(w)` |
| S5 | Friedman up to bounded loss `c*(k′) ≥ c*(k) − O(1)`, k′ ≥ k | integer version of S1; needs integrality-gap control | not by LP tools alone |

Logical map: C1 ⇒ S5.  S1(δ=0 ∀k) ⇒ S3 ⇒ (LP-provable values have no gaps).  S1 ⇒ S3 up to 2δ.  S2 ⇒ C4 at
thresholds.  R2 ∧ R3 ⇒ R1 ⇒ C7.  None of the S-statements gives A; none of the R-statements gives monotonicity.

### 9.5 How the current work feeds the ladder
* k²−3 bundle (`certificates/k2m3`, Lean `bentz_of_valid7`): `k_fam(3) = 6`, the template for "family from one box".
* k²−4 certificate run (κ = 0.08, box 9; `certificates/k2m4`, Lean `bentz4_of_validTilt9`): `k_fam(4) ≤ 8`; together with the per-k covers gives F₄, which now holds (k ≥ 8 resting on `ValidTilt9`).
* s(21), s(32), s(45), s(60) covers: `k_LP(4) ≤ 5`; the base cases a family needs.
* Seam capacity (§8, SEAM_W1, SEAM_1D): R2, and the ceiling for S1's splice patches.
* Corner LPs (QUADRANT, K2M4_MARGIN): R3.
* Descent experiment (§3): rules out translation descent for C1/S5; says nothing about the cover-side S1.

### 9.6 Next tests and lemmas, most tractable first
1. **S2 at c = 4**: box LP at k = 5 with slab-periodicity equalities in both directions (S60_COVER machinery).  Yes ⇒
   k²−4 tail from one k = 5 box.  No ⇒ measure `L(k) − S_ins(k)` for k = 5..8 (empirical δ).  ~10³–10⁴ CPU-s; ask.
2. **Splicing LP**: one wall crossing in isolation (arbitrary wall left and right, periodic patch of length ≈ 3.4);
   its exact cost at small w, the way SEAM_W1 did w = 1.
3. **Dual-side band LP** at w = 2, 3: rigorous caps on `M(w)` (SEAM_W1 §4).
4. **Corner dual at w = 1**: explain D = ½ (corner cost exactly ¼?) by a few squares.
5. **Lemmas now in hand** [proved]: averaging lemma (9.1.1); seam bound `D ≤ m_v` and `E(n) ≤ 0` (QUADRANT §1.4);
   `M(1) ≤ ¾` and the w − ¼ cap for seam-avoiding certificates (SEAM_W1 §1, §3); unfilled-area duality (SEAM_W1 §4);
   insertion lemma (§2.1).  **Lean candidates**, in order: `E(n) ≤ 0` (nearly a corollary of `packing_le_measure`:
   n² disjoint closed squares on a stretched grid, then continuity from above) and the seam bound; then `M(1) ≤ ¾`
   (one square, triangle area, δ → 0).  Not yet: validity of the w = 1 profile (Cavalieri + trig casework; the
   reusable piece is the slice lemma `N(y) ≥ 1[W(y) ≥ 1]`).
