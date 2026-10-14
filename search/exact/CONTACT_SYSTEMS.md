# Contact systems, counted and solved by hand (2026-10-13)

Why the record sides are algebraic, and where their degrees come from, made explicit on s(5), s(11), s(17) and s(29).
Started for an exposition of the s(29) degree-7213 result (`homotopy/N29.md`).  Scripts in `explain/`, run from the
repo root; they read the regularized records (`../regularize/lists/certs/`, the view drawn on the site).

## Formulation

Square i: centre (xᵢ, yᵢ); a tilted class has sides along (c, s) and (−s, c) with c² + s² = 1 (no trig).  The point
a steps along and b across from the centre is (xᵢ + a·c − b·s, yᵢ + a·s + b·c); the sides are a or b = ±½.

* Corner on a side of another square: corner = that point with a or b = ±½: 2 equations, 1 new unknown (where along
  the side).  Equivalently, after eliminating the slide, one equation (P − Cⱼ)·n = ½ (`fullsys.py` uses this form).
* Two parallel squares sharing a side: centres 1 apart along the normal, offset t: 2 equations, 1 new unknown.
* Corner on a wall: one coordinate = 0 or S: 1 equation.

Each contact is one net condition.  `contacts.py` reads them off a record (tolerance 1e-9); a corner at a corner gives
both side labels (fixed 10-13: taking only one left a stacked upright square free).

## Counts (`fullsys.py`, `countn.py`: independent contacts by Jacobian rank at the packing)

| n | unknowns (2n + S + 2·tilts) | independent conditions | free directions | trivial contacts (upright only) | contacts with tilted squares |
|---|---|---|---|---|---|
| 5 | 13 | 12 | 1 | 8 (walls) | 4, only 3 independent |
| 11 | 25 | 25 | 0 | 11 (9 wall, 2 pair) → 14 unknowns left | 13 (11 pair, 2 wall): 25 eqs in 25 unknowns with slides |
| 17 | 39 | 39 | 0 | | |
| 29 | 69 | 69 | 0 | 30 (18 wall, 12 pair) → 39 left | 34 (32 pair, 2 wall): 71 eqs in 71 unknowns with slides |

* s(5) is rigid but not first-order rigid: the contacts leave the rotation of the middle square free, S(c, s) =
  2 + 1/(c + s), stationary at c = s (dS = 0 is the extra equation).  Same mechanism as s(83) (WISHLIST G).
* s(11), s(17), s(29) are first-order rigid in the regularized view.  For s(29) the regularization pushes the 6 loose
  squares (N29.md: flat squares 4, 8, 11, 13, 18, 24) into contact; in the record they slide and S ignores them.

## s(11) by hand (`symsys.py 11 --solve`)

The 13 contact equations with a tilted square, upright squares pinned (x = k/2 or S − k/2), plus c² + s² = 1;
lex Gröbner basis → S·(S − 4)·(x⁸ − 20x⁷ + 178x⁶ − 842x⁵ + 1923x⁴ − 496x³ − 6754x² + 12420x − 6865).
Trump's octic, reproduced from the contact graph alone (S = 0, 4 are spurious branches).

## Phantom solutions (`phantom11.py`, `phantom_check.py`)

The octic's 8 roots: the record 3.8771, a real root S = −1.8530, and three complex pairs.  At S = −1.853 every
equation holds (residual 4e-16), but the box runs from −1.85 to 0, squares overlap, and most corners lie on the
*extension* of the side they "touch": the equations encode "on the line" and none of the inequalities (within the
side, no overlap, inside the box).  Position along the side, record → phantom: corner of 3 on 9, −0.48 → +2.27;
3 on 11, +0.24 → +1.82; 2 on 10, −0.09 → −1.80; 1 on 7, +0.17 → −3.12; 4 on 7, +0.07 → +3.58.  By irreducibility the
minimal polynomial cannot separate the record from these, so they all count toward the degree.

## Degree vs coupled angles

Among the 269 exact forms (plus s(29), s(83)), angle classes the side depends on (`minpoly.py` `free_classes`) against degree: 1 → ≤ 8
(s(11)), 2 → 2–42 (s(17) 18, s(41) 42), 3 → 12 (s(51), stationary) and 672 (s(83), stationary,
Ellsworth's polynomial), 5 → 7213 (s(29), rigid).  Plain counts of distinct angles do not work: s(175) has 14
distinct angles and degree 2, s(308) 31 and degree 1 (their tilted squares fill slack).
