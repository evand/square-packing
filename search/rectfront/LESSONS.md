# Unit squares in rectangles: lessons learned (R₃, R₄, R₅), 2026-10-21

For the next session.  What is proved, which proof techniques exist for these problems, which worked
here and why, which did not and how we know, and where R₅ stands.  Numbers quoted are from the scripts
named next to them (`search/rectfront/`, outputs in `runs/rectfront/`).  Evan's guidance (10-21): keep the
Lean simple and educational; R₃–R₅ should show techniques that keep larger packings simple; no extra
effort on extensibility until an approach has proved itself.

## 1. Status

| | Statement | Lean | Notes |
|---|---|---|---|
| R₂ | max side ≥ 2 | `rect2_iff` (`Rect.lean`) | centre point lies in every square |
| R₃ | `rect3_iff`: one side ≥ 3, both ≥ 2, or (short side s < 2) ℓ + 2s ≥ 4 + 2√2 | `Rect3.lean` (171 lines) + `RectChain.lean` | chain, 3 × 3 witness cases |
| R₄ | `rect4_iff`: one side ≥ 4, both ≥ 2, or s < 2 and ℓ ≥ L(t, s) for some t (Friedman/Morandi piece of S₃) | `Rect4.lean`, `Rect4Core.lean` (592 lines) | chain + one-square lemma + half turn |
| R₅, s < 2 | L ≥ 5 or L + 4H ≥ 4 + 6√2 (pieces 4.1, 4.2 of S₄) | `Rect5Mid.lean` (middle square), `Rect5.lean` (`rect5_strip_lower_of`, given the real inequality `Final5`) | **`Final5` not proved** (§5) |
| R₅, s ≥ 2 | b = 3 on 2 ≤ a ≤ 1 + √2; a + b = 4 + √2 on 1 + √2 ≤ a ≤ 2 + 1/√2 (contains s(5)) | — | not started (§6) |

All boundaries are Friedman's (Math Magic, June 2010; closed forms also in wand125 `problems/rectangles`,
read via `git show origin/main:` in `~/math/_untrusted-third-party/wand125-square-packing`, see that
directory's CLAUDE.md).  What is new is the lower-bound proofs.

## 2. Proof techniques, and what we know about each

### 2.1 The chain (strip of height H < 2) — **our main tool; proved R₃ and R₄**
Sort the squares by centre x.  Because H < 2, any two have vertical centre offset < 1, so a separating
line between consecutive squares has normal with u₁ > 0 (`separate_lr`, from Mathlib's Hahn–Banach on the
open squares).  Pieces (all in `RectChain.lean` unless noted):
* **Upright replacement** (`upright_repl`): an upright square pushed against the wall lies on the same side
  of any separating line as the tilted one.  Gives **shelves**: every point of the second square at the
  heights of some band [lo, lo+1] has x ≥ 1 (`left_shelf`; `right_shelf` by mirror).
* **Witness** (`witness`, `Wit`): which point of a square a shelf constrains — its side vertex, or a point of
  the edge above / below it.  Three cases, each a short real inequality.
* **Frames** (`framed`): a square as (P, Q) ≥ 0 with P² + Q² = 1; no angles, no sign juggling.
* **One-square lemma** (R₄, `one_square`): a square with a left shelf, in a half-plane with normal (S, −C),
  has support ≥ upright bound or ≥ the "frontier" bound C(S − C) + 1 + C + S − H over C.  Proof by the three
  witness cases plus two two-angle inequalities (`core_D1`, `core_D2`) that **factor in half-angle
  coordinates** u = tan(θ/2), v = tan(t/2) (factors like (u − v), concave quadratics checked at the ends).
* **Symmetries** (`mir`, `vmir`, `rot`): mirror and half-turn turn a right shelf into a left one, so one lemma
  serves both sides and both slope signs.

Why R₄ was clean: one separating line (b|c); both b and the half-turned c see **the same line**, so the
two one-square bounds add directly (`combine4`) and the line's direction becomes Friedman's angle t.
Every step was exact, checked numerically first (`r4/relax.py`, `sweep.py`: chain relaxation equals the
frontier to 1e-15).

**Limits.**  Needs H < 2 (else squares stack and the x-order is no chain).  With two middle lines (R₅) the
two halves see **different** lines and the coupling through the middle square is real (§5).

### 2.2 Exact case enumeration by separating axis (SAT) + LP duality + interval boxes — **tried for R₃, abandoned**
`r3bb.py`: per pair an edge normal of one square (512 types for 3 squares), LP duality for the centres,
boxes in tilt space.  86k boxes certified, but a tube of 423 boxes around lines of equality would not
close; the chain proof replaced it.  Lesson: **equality on a curve** (not a point) kills box methods.

### 2.3 Polynomial positivity certificates (Handelman LP + Lean `linarith`) — **works, but heavy; R₅ fallback**
`r5/cert5.py`, `gen5.py`, `phase1.py`, `phase2.py`, `tight5.py`: each leaf of `Final5` becomes
"P(u, v, H, r) ≥ 0" (r = √2) on a box with region constraints; a float LP proposes products of
constraint polynomials; Lean's `linarith` finds exact coefficients.  State: 160 distinct leaves, 120
certified in `runs/rectfront/r5cert/plans.pkl`; the tight ones done by hand (below).  What we learned:
* Float LPs are only advisory: HiGHS meets equalities to ~1e-7, so a support can be exactly infeasible
  (two tight-leaf bands failed in Lean).  Margins (constant ≥ δ) and exact repair are needed.
* Zeros on whole box edges (from denominators u, 1 − u, 1 − v …) defeat Handelman: **strip edge factors
  first** (`strip_edges`), and **subdivide** (u, v at ½) when needed.
* Zeros at **irrational interior points** (the zigzag u = v = √2 − 1) cannot be certified with rational
  bounds on √2; they need exact algebra (r² = 2) and a structural decomposition — here
  (1 + u)ᵐ G = R + k Q with R a positive multiple of (u + 1 − r)² (+ (2 − r²)), k the kink factor, Q of one
  sign near the zigzag (`tight5.py`; Lean draft in `r5/lean-attic/Rect5Tight.lean`).
* Cost: hours of generator work, ~1 h of Lean compile projected, and a proof that teaches nothing.  Evan
  (10-21): "On R3 it seemed to be far too complex, and other approaches worked better."  Paused.

### 2.4 Unavoidable points / pigeonhole — **the s(5) proof; candidate for R₅ with s ≥ 2**
Göbel 1979 (open at CWI, ir.cwi.nl/pub/12685, Prop. 1, a half-page with a small gap); clean version
Friedman DS7: four points such that every unit square in a side < 2 + 1/√2 container contains one in its
interior, plus one square per point → at most 4 squares.  **Formalised**: chelokot's archive
`formal/SquarePackingArchive/Records/Square5.lean` (`s5_ge_goebel`, `s5_eq_goebel`; via
`lowerBound_succ_of_unavoidable`).  (Our `lean/` has only the upper bound and the toy s(5) ≥ 2.)
Rectangles (literature agent's analysis, partly numerical, `literature/s5-proofs-2026-10-21.md`): points at
(1,1), (a−1,1), (1,b−1), (a−1,b−1) with a new centre-rectangle lemma (a unit square whose centre is in an
X × Y rectangle contains a corner iff X + Y ≤ √2 and X, Y ≤ 1; numerical, unproved) match R₅'s boundary
exactly; but the **side strips** (gap + 2·distance-to-wall < 2√2, sharp) only allow
a ≤ 3√2 − 2 ≈ 2.243 or a ≥ 4 − √2 ≈ 2.586.  For a ∈ (2.243, 2.586) a crude search found no 4-point set.

**Limits of point arguments.**  A fixed point set certifies a *fixed* container (or a family only if the
lemmas are parameterised).  Exactly on the boundary the margin is 0, so points must sit exactly on the
equality configuration; the side-strip lemma is sharp and is what fails.  Weighted point covers
(`CovMP`, zero-margin checkers) are the same idea with weights: excellent for one container with slack,
**cannot reach a boundary with equality by finitely many certificates** (each certificate needs margin).

### 2.5 Contact graphs / arrangement enumeration — **diagnostic, not (yet) a proof**
Numerically the optimum of each relaxation is a rigid arrangement (contact graph): R₃ wedge, R₄ tilted
pair, R₅ zigzag.  Enumerating contact graphs gives the **candidate** boundary pieces (that is how Friedman
and wand125 got them) and tells which contacts are tight.  Limits as a *proof*: one must show the minimum
is attained at a rigid arrangement (compactness + first-order conditions / no feasible improving motion),
which is hard to formalise; and contact equations need not have closed forms (R₅ vertex contact, §5,
gives a cubic).  Useful *inside* a proof to choose cases (e.g. "edge-on-edge" vs "vertex-on-edge"
contacts of two squares), as in §5.

### 2.6 Small exact tools that keep proofs short (reuse these)
* Half-angle substitution + factorisation for two-angle inequalities (R₄ `Rect4Core`).
* **Kink / concavity**: a function linear (or a positive sinusoid) between kinks is minimised at a kink
  or an end — `arc_lin`, `kink_min` (R₅ middle square); "linear in cot t between kinks" for the upright
  branch of the one-square bound.
* **Linearity in a parameter**: H enters linearly → check endpoints; a min of linear functions is concave
  → check the ends of the range (c's height in §5).
* Exact multipliers for `linear_combination` from `sympy.reduced` (certificates for identities modulo
  P² + Q² = 1, C² + S² = 1); keep Lean steps `linarith` with explicit product hints, not `nlinarith`
  searches (heartbeats; one case per declaration).

## 3. Method lessons
* **Numerics first, exactly the statement Lean will see.**  Every reduction was checked as an
  optimisation problem before any Lean (relaxation = frontier to 1e-15; per-case worst slack tables).
  Several "obvious" decouplings were false and were caught this way (§5).
* **Find the structural reduction before certifying.**  R₃ (chain) and R₄ (one-square lemma) became short
  once the right lemma was found; brute force (R₃ tube, R₅ leaves) did not converge to something good.
* When an inequality is tight, find *where* (point / curve) and *why* (which contacts); a tight curve
  usually means a factor (kink factor, (u − v), (u² + 2u − 1)²) that a certificate must carry explicitly.
* Look at a candidate statement's sharpness numerically: just below the boundary no witness should exist,
  just above one should (R₄: `rect4_lower`'s ∃(C, S) checked at F(H) ± 1e-3).

## 4. Practical lessons (tools)
* `pkill -f PATTERN` kills your own shell if the pattern is in the command line: use `pgrep -f "[p]attern"`.
* `ProcessPoolExecutor.map` emits in order: one slow job hides all later output; use `imap_unordered`.
* Lean: `le_or_lt` is now `le_or_gt`; `push_neg` → `push Not`; per-declaration heartbeats — split big case
  analyses into lemmas; `simpa using` to pass hypotheses across `if`/`Bool` normal forms.

## 5. R₅ with short side < 2: where it stands

Target: L ≥ 5 or L + 4H ≥ 4 + 6√2.  Extremal (piece 4.2): uprights in both bottom corners, b and d at 45° on
the ceiling against them, c at 45° on the floor, edge to edge (the **zigzag**).

Proved in Lean: the chain for five squares up to one real inequality (`Rect5.lean`: `Final5` ⇒
`rect5_strip_lower`), and the middle-square lemma `mid_set` (width of c between lines b|c and c|d is at
least its width upright or aligned with one of the lines; `Rect5Mid.lean`).

Facts established numerically (scripts in `search/rectfront/r5/`):
* `relax5.py`: chain relaxation (shelves on b, d; b|c, c|d disjoint) is **exact** (1e-15).
* `split5.py`, `cands5.py`: the bound B_b(line 1) + W_c + B_d(line 2) is exact; W_c's minimum over c's tilt is
  at c upright or aligned with a line (kink argument; now in Lean).
* `cases5.py`: of the 48 cases of `Final5`, 4 are tight (opposite slopes, both frontier branches, c aligned:
  the zigzag); 12 tend to 0 only as both lines turn vertical (row bound 5); the rest have slack ≥ 0.03.
* Decouplings that **fail**: halves independently (`halves5.py`, −0.1); c's tilt or height per half
  (`partial5.py`, −0.1); c's support ≥ 1 across the other line (`aligned5b.py`, −0.23).
* **Separating-axis does not shrink it**: line along b's edge + line along c's edge is exactly the
  "c aligned with one line" case, which is tight at the zigzag (b aligned gains nothing in the one-square
  bound, which is already tight for aligned b).
* **Plan B (most promising), `xsplit5.py`**: split at c: L ≥ X(c) + X(mirror c), X(c) = least centre x of c
  given a, b (minimum over line 1 of b's bound + c's support).  Exact.  For fixed lines the bound is linear in
  c's height y, so the sum is concave in y → **c on the floor or the ceiling**; the final inequality is
  **one-variable in c's tilt θ**, tight only at θ = 45° (zigzag) and θ = 0° (row corner, vertical lines).
  X's minimising line is along c's edge (edge contact), vertical (row), or — only where slack ≥ ~0.09 —
  interior (**vertex contact**: c's corner on b's edge, b in its frontier position; optimal tilt of b solves
  S³ − C³ + (1 − H)(S² − C²) − (y − 1)C² = 0, no closed form).
  * Upright branch of b's bound: linear in cot t between kinks → only edge contacts and the vertical line
    (closed forms in θ).  Fine.
  * Frontier branch: the interior minima.  Relaxations tried and **failed**: c's support ≥ 1 (`xlb5.py`,
    loses 0.5 at θ = 0); tangent lines for csc, sec at 45° (`xtan5.py`, exact at 45° and fine for H < 1.7,
    but diverges as t → 90° for H > 1 + 1/√2); unrestricted envelope Λ(y) of b's lines (`xvert5.py`,
    → −∞: vertex terms must be restricted to the angle range where that vertex faces the line).

Recommended next step (Evan agreed to Plan B; this choice was pending at the pause): **edge contacts
exactly, vertex contacts crudely** — keep the one-variable edge-contact formulas for the tight structure,
and bound the frontier branch on the vertex-contact angle range by something crude but valid for all t
(likely a case split at a fixed line angle, using the upright branch near vertical lines), proved as a
small two-variable inequality *with slack*.  Before any Lean: check numerically that the crude bound plus
edge contacts still gives slack ≥ 0 for θ ∈ [0°, 90°], H ∈ [1, 2), c on floor and ceiling.
Alternative if that stalls: parametrise vertex contacts by b's tilt (the cubic becomes a curve definition).
Fallback: the certificate pipeline (§2.3; ~120/160 leaves done, tight leaves drafted).

## 6. R₅ with short side ≥ 2 (not started)
Contains s(5).  The chain does not apply (stacking).  Start from §2.4: formalise (or reuse the idea of)
the unavoidable-point proof, add the centre-rectangle lemma (X + Y ≤ √2, X, Y ≤ 1) and see how far the
side-strip lemma reaches; the middle range a ∈ (3√2 − 2, 4 − √2) needs a new idea (second-stage case
split or a weighted cover per sub-range).  Note the b = 3 piece (2 ≤ a ≤ 1 + √2) and the a = 2 step.
