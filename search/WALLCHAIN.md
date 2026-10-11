# wallchain: exact elimination along wall chains for the n = 17 hard tail — first demonstrator  (2026-10-10)

Code: `search/wallchain/wallchain.py` (float DP), `search/wallchain/subsets.py` (float placement search; failed its
control, see §2).  Input read (not executed) from the local jlevy/squares clone: exp-247's 24-cell cover at
`U = 1169/250`.  Labels: **[measured]** float, **[heuristic]**.

## Idea (from the X-052 reading)

X-052 §2.5 diagnoses the distance-2 hard tail (e.g. m851903, m1964767) as wall crowding — "five squares on a wall
of length 3.676 whose neighbours fail to clear by 0.011–0.018" — and asks for "joint facts about
orientation-dependent extents along full walls".  Squares in the corner and side cells form a ring; along one wall
the problem is a path, so exact variable elimination (CAD's projection, but along a path, where it stays cheap)
should be complete where pairwise consistency is not: carry `F(q, t)` = least position along the wall for a
square at depth `q` and angle `t`, step with the closed-form least centre gap `b(dq, tA, tB)`
(separating axes; the allowed gaps are an up-set because `|dq| < 1` and the Minkowski difference of two unit
squares contains the unit disc), turn frames at corners.  Dropping non-neighbour pairs is a relaxation, so an
infeasible chain would be a sound exclusion once the grid is made rigorous.

## Results  [measured]

1. m851903 decodes as: all four corners, 11 side cells (no `side-W0`), `interior-W`, `interior-S`.  Its ring is
   one open path of 15 squares (three full walls S, E, N plus `W2`, `W1`).
2. **The neighbour-only wall chain is loose** for m851903 and for the family's state (control): every square ends
   with about `0.65` of slack to the end of its cell (grid 50/unit depth, 2° angles).
3. Why: side cells are `0.91`–`0.93` deep, so neighbours may stagger in depth, and tilted neighbours at staggered
   depth need far less than unit spacing — `b(0.93, 30°, 30°) = 0.618`, `b(0.93, 0°, 45°) = 0.777` (axis-parallel
   neighbours still need exactly 1 for any `|dq| < 1`; unit-tested).  A tilted staircase fits five squares in
   about 2.5 of wall.
4. jlevy's stall review (2026-10-05 §1.5) records for m851903: full penetration `0.0115`, smallest unplaceable
   sub-pattern of arity 15 (`0.0031`), and "each wall's middle square is the least supported owner, held by its two
   wall neighbours by 0.013–0.019".  Together with 2–3: the crowding is real but **is not a fact about one wall
   row**; it is the row's inward intrusion (staircases) against the interior squares and against the adjacent
   walls' rows near the corners — non-neighbour pairs the chain relaxation drops.

## 2. The float placement search failed its control

`subsets.py` (random starts + L-BFGS-B on squared penetration) did not place the family's own state at `U`
(best `sqrt(viol) = 0.089`), and gave a 16-square superset a smaller violation than its 15-square subset.  Its
numbers are discarded.  A usable placement search must be seeded from known packings (as jlevy's survey is).

## Reading  [heuristic]

The same mechanism defeats the 1-D chain arguments at `n = 12` (`T11_CHAINS.md`: the failing hypothesis is rise
control; band links are tilted staircases with credit `tan(tau)|Dperp|`).  In both problems a chain along one
direction is beaten by transverse staggering, and the content is in what the staggering costs transversely.  So
the object to compute is not a wall row's length but its **trade-off between length and inward intrusion**: the
set of achievable inner envelopes of a wall row, or a scalar shadow of it (intruded area / waste, which is additive
along the chain and so still DP-able), coupled to the interior and to the neighbouring walls at the corners.

## 3. Ring DP conditioned on one interior square (`ringcond.py`)  [measured]

Fix one interior square's pose (an obstacle), run the ring as a path DP: a later square only needs the least
*allowed* position of its predecessor, so the state stays `F(q, t)` (obstacles handled by a fine position grid).
Relaxation: non-consecutive ring pairs and the other interior square are dropped.  Float grids round against
exactly tight configurations, so squares are shrunk to side `1 - 0.003`; with that the family control (its own
interior squares as obstacles) is feasible (it fails at shrink 0 — the family has ~`5e-4` slack at `U` — and
failed before a corner-turn grid bug was fixed).  Grids: depth 0.02, position 0.004, ring angles 2°; interior
poses 0.05 / 4°; 10 + 4 workers, ~1.2 CPU-h per conditioned square.

| state | conditioned on | interior poses leaving the ring feasible |
|---|---|---|
| m851903 | interior-S | 2,822 / 4,577 (62 %) |
| m851903 | interior-W | 3,345 / 4,577 (73 %) |
| family `3439615` (control) | interior-W | 3,769 / 4,577 (82 %) |

**Reading.**  The width-1 path relaxation (ring neighbours only, plus one interior square) sees almost nothing of
m851903's obstruction: most interior poses survive, and the difference from the true family state is small.
Together with jlevy's float survey (nothing smaller than 15 of the 17 cells is infeasible), the obstruction is
not path-decomposable at width 1: it lives in the non-neighbour pairs the relaxation drops — staircased wall
squares meeting across corners, i/i+2 pairs inside a compressed staircase, and both interior squares at once.

## Next steps (not started)

* Intrusion DP: per wall, minimise intrusion (area or depth profile) subject to the row fitting; check whether
  four walls' minimal intrusions plus the two interior squares exceed what the container can hold.
* Width-2 band DP around the ring (state = two consecutive squares) to keep the cross-corner pairs: exact but
  6-dimensional state.
