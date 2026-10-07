# DRAFT: Local optimality of the known square packings (outline and pieces, 2026-10-06)

Status: working draft for a site page; meant to be citable from the eventual s(11) explainer.  Content first, not
prose-polished.  `[TODO]` marks things a fresh session should check or fill in.  The numbers come from
`search/exact/batch/` (2026-10-05) and `tasks/exact-minpoly/` (2026-10-06).

---

## 0. Why this page

The s(11) proof (a global optimum, in Lean) prompted the question "what about the other n?"  For almost every n the
honest answer is: *we can't prove the record is the best possible, but we can say a lot about it locally.*  This page
spells out what "locally optimal" means for square packings, what is known for every record with n ≤ 324, and which
parts are proved, which are high-precision computation, and which are in Lean.

## 1. A ladder of statements about a record packing

From weakest to strongest:

1. **Valid**: the n squares fit in a square of side S, so s(n) ≤ S.  Checkable with exact rational arithmetic.
2. **Exact value**: S is known exactly, as a closed form or as *the unique root of an explicit integer polynomial in
   an explicit interval* (e.g. s(17) ≤ the root near 4.6755 of a degree-18 polynomial).
3. **Locally optimal**: no small motion of the squares lets the box shrink.
4. **Globally optimal**: S = s(n).  Proved only for a few n (`[TODO]` list them with credits: small n classically,
   s(11) 2026 in Lean, our s(13) and s(32) in Lean, the k² − c families, …).

Levels 1–3 can be checked for any n, mechanically; level 4 is a research project per n.  The gap between 3 and 4 is
real: 2026 alone saw dozens of records improved (Couzo, Chang, Schadt, …), and each replaced packing was, as far as
anyone checked, a local optimum.

## 2. What "locally optimal" means, precisely

Configuration: centres and angles (x_i, y_i, θ_i) for i = 1…n, and the side S.  Constraints: each square inside
[0, S]², and interiors pairwise disjoint.  **Locally optimal** means there is an ε > 0 such that every valid
configuration within ε of this one needs a box at least as big.

The subtleties are where the explaining is needed:

* **Rattlers.**  Squares that touch nothing load-bearing can wobble freely.  They don't affect S, so local optimality
  is about S, not about the configuration being isolated.  Many records have rattlers (n = 10 has 2; n = 17 has 1).
* **Flat motions.**  Whole families of packings with the same S: a square sliding along a channel, a row shifting.
  296 of the 323 records have them.  The right notion there is "strict local minimum *modulo* flat motions": every
  motion that isn't flat increases S, to second order.
* **The constraint isn't smooth.**  "Disjoint interiors" is not one smooth inequality.  Near a given configuration,
  though, it implies a list of smooth ones: corner a of square j stays on the outside of side k of square i,
  `n_ik · (P_ja − c_i) ≥ ½`, plus the wall conditions.  The exception is **corner-to-corner touches**, where either
  of two separating lines could work: an "or" constraint.  `[figure: a corner-corner touch, two candidate lines]`
* **Jammed vs. locally optimal.**  "No motion keeps the box size and moves something" (rigidity) is different from
  "no motion shrinks the box" (local optimality).  Example: the register's n = 105 witness looked jammed but **was
  not a local minimum**.  There is a first-order descent that uses the "or" at corner-to-corner touches.  Its
  re-optimised version is one.

### Force balance (the physics picture)

Push the walls inward with unit pressure.  The packing is a first-order candidate exactly when contact forces λ ≥ 0
can balance it (the KKT conditions: Σ λ_k ∇g_k = ∇S).  The contacts that can carry force form the **load-bearing
network**; squares outside it are rattlers or force-free.  `[figure: force network of s(11) or s(17), arrows ∝ λ]`

Three grades of local optimality, from strongest to weakest:

* **First-order rigid with strictly positive forces**: every motion either opens a contact or grows S at first order.
  The strongest kind, and the easiest to prove.  Among all 324 records: the 20 perfect squares (trivial grids) and
  **n = 11**.  Nothing else.
* **Strict at second order**: the forces balance, and the curvature (reduced Hessian) is positive.  n = 5
  (s = 2 + 1/√2) is the cleanest example.
* **Strict modulo flat motions**: second order, after quotienting by the flat families.  The typical record.

## 3. What we know for n ≤ 324 (2026-10-05/06)

Method: for each register record, solve for the nearby exact KKT point at 80 digits (`search/exact/exactsolve.py`).
That finds the contacts, the load-bearing set (as a max-support LP), and the KKT system, then runs Newton.  Then
check force balance, second order and flat modes, and the corner-to-corner "or" constraints (a MILP).

| statement | count | status |
|---|---|---|
| valid, with an exact rational certificate, checked by two independent exact verifiers | 323 / 324 | **proved** (computer-checked exact arithmetic, not yet Lean) |
| numerically a local minimum (force balance, second order modulo flat motions, MILP) | 317 / 324 | **numerical** (80-digit point, f64 checks; no interval enclosures) |
| not yet a local minimum, even numerically | 177, 230, 261 (a contact that carries zero force: the second-order test doesn't decide), 211, 263 (degenerate zero modes), 272 (negative curvature; superseded by a new record that passes), 292 (no exact point yet) | open |

Curiosities for the page:
* 50 certified values lie 3e-13 … 5e-11 *below* the register's printed bounds.  That's the same packing pushed to its
  exact optimum, not a new structure.
* s(172) and s(199) are exactly 1 apart.
* s(230) = 643/41 and s(261) = 684/41 are rational.

## 4. Exact values

Every record's S is an algebraic number: the contact equations are polynomial, with angles as (c, s), c² + s² = 1.
How explicit we can make it, for the 323 certified records:

| | count |
|---|---|
| axis-parallel only: S is an integer | 176 |
| tilted, polynomial known: found by integer-relation search from our digits (78), or published by Ellsworth / Chang and matching ours to 60+ digits (17 more) | 95 |
| tilted, polynomial not yet known | 52 |

Examples to show:
* s(5) = 2 + 1/√2;
* s(18) = (7 + √7)/2;
* s(11): a root of the octic x⁸ − 20x⁷ + 178x⁶ − 842x⁵ + 1923x⁴ − 496x³ − 6754x² + 12420x − 6865;
* s(17): degree 18;
* **s(83): degree 672 with 724-digit coefficients** (Chang 2026, verified by Ellsworth to 10⁶ digits).

Two ways to get the polynomial:
* **Integer relations (LLL/PSLQ)** from many digits.  Cheap, but limited to degree × coefficient size ≲ the number of
  digits.  n = 17 needs ~2000 digits; n = 83 would need ~500,000.
* **Elimination from the contact graph.**  Turn each contact into an equation, give parallel squares a shared angle,
  and note that once the angles are fixed, every contact is *linear* in the positions.  Then eliminate (msolve).
  n = 17 takes 0.2 s and reproduces the same degree-18 polynomial.  Its by-product is the packing written exactly
  over a number field, which is what a proof needs.  For records where force balance, not contacts, fixes S (n = 83
  is one), the force-balance equations are added; eliminating the forces gives David Ellsworth's Jacobian-determinant
  conditions.

## 5. Proved vs. computed vs. Lean: the honest table

| | proved (exact arithmetic) | high-precision numerics | in Lean |
|---|---|---|---|
| validity, s(n) ≤ S' (S' = S + ~1e-19) | 323 | — | grids only (`packs_grid`) `[TODO: the upper-bound verifier would make this 323]` |
| exact value of the record | closed forms; the rest pending exact field verification | 271 with a known polynomial (176 of them trivial integers) | none yet |
| local optimality | none yet | 317 | none yet `[TODO: Lean demo n = 5 / 11]` |
| global optimality | s(11) (others, Lean), s(13), s(32) (ours, Lean), … | — | those |

What would close the gaps:
* **exact values**: verify the exact number-field configuration (contacts are identities, open gaps checked by
  interval arithmetic);
* **local optimality**: the same exact point plus a force certificate (λ) and a curvature certificate.  The
  first-order-rigid case (n = 11) needs no curvature;
* in Lean: a definition of local optimality, plus one lemma that turns "valid packing near this configuration" into
  the smooth inequalities g_k ≥ 0.

## 6. Local vs. global, with examples

* The landscape has many local minima.  For n = 110 our sampler finds funnels; the record sits in a narrow one.
  `[TODO: link the s110 landscape case study when written]`
* Records improve by jumping to a different local minimum, not by sliding: Couzo's 2026 batch, Chang's s(83), our
  s(266), s(270), s(272).
* A record can fail to be even locally optimal (the n = 105 register witness, until polished).  So "local
  optimality" is a real check, not a formality.
* s(11): the global proof makes local moot, but s(11) is also the one non-trivial record that is first-order rigid
  with all contact forces positive.  It is locally optimal in the strongest sense.

## Pieces left for a fresh session

* Figures: a force network (s(11) or s(17)), a corner-to-corner "or", a flat motion (a rattler, a sliding row).
* The list of globally proved n with credits (section 1, item 4).
* Site page (HTML, site style), linked from `problems.html` and the s(11) explainer.
* Update the counts as `tasks/exact-minpoly/` progresses, and if the new records are posted.
* Decide whether to include the unposted records (266, 270, 272); ask Evan.
