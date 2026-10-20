# BULB: tilted wall bands with tapered ends — one family for c*(k) constructions and the two-bulb records (10-19)

Notes from a discussion with Evan (10-19 evening), written down before stopping.  Nothing run beyond a smoke test.
Code template: `bulbgen.py`.  TODO line: "Constructors ... two-bulb, lopsided bulb".

## Idea (Evan)

The tilted-band constructions for c*(k) (Ryu's Theorem 1.3 L, our filled-end R1′; `literature/waste-constructions-2026-10-17.md`
§2–3) and the two-bulb records near k² (s(308), and the V-band records like 211) are the same thing: wall bands of
b-wide tilted rows at about θ₀(b) = 2 arctan(1/b).  The records round, taper or shear the band ends, which packs the
ends and the junction better than our plain L with filled corners.

## Measurements

* **s(308)** (`candidates/s308_d308`, side 17.99715, c = 18² − 308 = 16).  Two bulbs, no junction:
  * side bulb, 53 squares, flush on the right wall, about −19°, y ≈ 6.3–16.7; thickness per 1.5-unit slice
    5, 7, 8, 9, 9, 8, 6, 1;
  * bottom bulb, 49 squares, flush on the bottom wall, about +18°, x ≈ 1.3–10.6; 1, 6, 8, 9, 9, 8, 7, 1;
  * the bottom-right corner (about [11, 18] × [0, 6.3]) is a grid block.  The top row and the left column cross no
    bulb (about 1 waste each).  Each bulb crosses the rows or columns the other one misses;
  * body angles hold at 18–19°; the strays (3–28°) are in the tapers.
* **s(211)** (store 1044256, side 14.993550): one V band, 82 tilted squares.
  * Lower-right arm: 18 of 43 squares at 22.2°; θ₀(5) = 22.6°.  Evan: a clean 5 high in places, with a step up in
    position supported by a few angled squares (the strays at 18–20° and 26°).
  * Left arm: 4 wide in places, so the angle is wrong for its width (θ₀(4) = 28.1°).  It splits into 24.0° (12) and
    25.3° (18): the shear Evan points at.
  * So the arms can have different widths at roughly one shared angle, with the narrow arm making up the difference
    by shear.  Ryu's construction has one b for both bands.
* Rough waste count for 308 (estimate, not a calculation): total 15.9; bodies ≈ (rows − 1) tan θ ≈ 3 per bulb; uncut
  row and column ≈ 2; tapers and corner ≈ 8, so about 2 per end.  The R1′ formula at k = 18 gives c ≈ 24, Ryu's
  Theorem 1.4 32.

## Three ways to handle the far corner

* **L** (Ryu, R1′): one band runs the full wall (length k), the other stops at it (k − b).  The corner's rows are
  crossed by both bands.
* **V** (211; Ostashev's merged bands, `pk/bands.py` notes): the band bends through the corner.
* **Two bulbs** (308): neither band reaches the corner.  The corner is a grid block, each band is about k − b long,
  and all four ends are tapers.  The ends near the corner have to absorb ε along the band (a band's row count is
  continuous in its length, so it can).

## Constructor (`bulbgen.py`)

Evan: every parameter is a probability distribution, not a constant, even if it's just uniform or an implicit "go
until out of space".  Exact fill is unlikely to work well (`arms.py` header: idealised arms + exact fill lose about
20 squares at 211 / 241; small boundary deviations cost whole upright rows).  Use greedy approximate fill, largest hole
by sampling, or similar, and vary how many squares we fill.

Per start:
1. Draw: band widths b₁, b₂ (equal or ±1); angle θ₀(b) + noise, shared by both arms or one per arm; corner layout
   (L / V / bulb); per end, the taper (rows, tip width, which side loses squares), the gap to the wall or corner;
   wall gap; row shear; side s = k − δ.
2. Build the arms in the (k − δ) box: Ryu rows (square (i, j) at [i + j tan θ, i + 1 + j tan θ] × [j, j + 1] in
   the rotated frame), tapers as per-row width and offset, arm B the anti-diagonal mirror of a right-wall arm (angle
   sign flips, as in 308).  Drop out-of-box and overlapping squares.
3. Fill (10-20, `pk/fill.py`): fillfast exact or greedy (τ) uprights around the bands, then drop uprights above the
   target or top up below it with a placer (`bestofK`, `leastbad`, `farthest`; K drawn).
4. Count (Evan 10-20): target n + extra, extra in −2..2; `stop` = count (top up to the target) or full (until the
   placer reports no room, at most target + 2).  No forced insertion: each start is written as `s<seed>_n<count>.txt`
   and is a start for the n it reached.  10-20 smoke (8 starts): counts 295–310, base fill 279–291.
5. Optional fq quench (`--quench shrink|soft05|loose|...`, `relax.OPS`); then the explorer as for any start.

Mixing steps is expected: e.g. taper by deleting rows after the fill, a partial corner fill, stochastic fill rounds
between relaxes.

**Smoke test (10-19, k = 18, n = 308, 8 starts, 3 s):** geometry right (picture checked; arm B on the bottom after a
reflection fix); fills come up 12–35 short of n (the hole insertion makes up the count, so expect heavy overlap);
3 quenched starts (`shrink`, ~12 s each) went to 18.0000, 18.12, 18.04 — grid or worse.  The template is plumbing;
the fill and the distributions are where the work is.

## Perturbation version: re-angle across the band boundary (Evan, 10-19)

Possibly more efficient than building from scratch: start from an existing packing (record, near-record, or a
bulbgen start) and reshape its bands locally.
* Step: pick a square on the tilted/upright boundary.  An upright square next to a tilted one gets the neighbour's
  angle (the band grows); a tilted square next to an upright one goes to 0° (the band shrinks).  Probably mix both
  directions.
* Repeat a few times with locality correlation: the next square is drawn near the last one (a short random walk
  along the boundary), so one proposal reshapes one taper, corner or row end rather than scattering single flips.
* Then quench (overlap-tolerant: the rotated square overlaps its neighbours at first) and explore as usual.
* Parameters as distributions: number of steps, walk radius, grow/shrink mix, whether to copy the neighbour's angle
  exactly or the local mean angle, quench operator.
* Related existing moves, none the same: the explorer's `aswap` (swaps tilts between nearby squares that differ by
  > 3°; keeps the band shape; AS1 10-16: s271 was aswap-stable), ftmc's "add neighbour / re-angle" (10-04, with
  exact fill), `pk/moves.py` `crot` (rotates a cluster).  The new move changes which squares are in the band, so it
  moves taper and corner shapes directly.
* Natural home: a `pk/moves.py` move (it has the parameter-spec / logged-params machinery), so its outcome can be
  regressed on the parameters, and an explorer arm.  Test it like any new move: pre-register; refined outcome per
  CPU-h vs the existing arms from the same parents.

## Evaluation (pre-register before running)

* **Near the grid boundary:** at n just above k² − c*(k) (308 at k = 18; also 211-like, 153–156, 175–182,
  236–241, 268–273, 302–308), the distribution of starts: how many reach a side below k after equal refinement, how
  many land near the record and are new (matching distance, structure).  Baseline: explorer from the kicked record,
  same CPU.  First sanity check: does it hit the 308 and 211 basins at all?
* **Large n, quickly:** can it beat the structured packings (grid + border constructions, register entries beyond
  the hunted range) at large n with little compute?  That tests the core idea, separate from record hunting.

## Theory: a tighter explicit c*(k) family

Question (Evan): is there a systematic construction tighter than bands + filled L corners?  Probably a constant term,
but maybe one that smooths the c(k) curve.
* Bodies: Ryu's closed form; body waste ≈ length × tan θ, so total band length matters (two-bulb ≈ 2(k − b) vs L ≈
  2k − b).
* Ends: a finite gadget table.  Ryu's Lemma 6.1 lets the region past a row line be replaced by any packing of it, so
  "best end for width b (and lift y)" is a finite search per b, certified exactly once and reused for all k.  308's
  tapers trim whole rows, so the cut sits a few rows deeper; the same argument.
* Corner layout: L / V / two-bulb accounting via the waste identity (Prop. 6.3), including the uncut row and column
  of the two-bulb layout.  Not yet done.
* Smoothing: R1′ is jagged from rounding b (excess over 4√(2k) creeps to 3.68 by k = 3000).  Mixed widths (b, b + 1),
  the corner choice and gadget ends give more knobs; taking the minimum should flatten it.
* Stays Lean-friendly: closed-form bodies, finite gadget certificates, one disjoint-regions gluing lemma (shared with
  Fang's gluing).
* Regime where it matters: from where Fang's linear bound stops winning (k ≈ 38 for R1′) up to where Ryu's k^{2/5}
  staircase takes over.

## Next

1. Fill: done as a module (10-20): `pk/fill.py` + `fillbench.py` (PACKER "Fill module").  Next: the refill bench
   (placer x relax x stop) before tuning bulbgen's fill choice; then the deficit / count distribution of bulb starts.
2. Distributions: fit them to 308 / 211 / the V-band records (widths, angle offsets, taper profiles, corner gaps)
   instead of first guesses.
3. (Evan 10-20: no rush; we want good test data for the explicit band constructors first.)  Pre-register and run the near-boundary check at 308 and 211 vs the kicked-record baseline.
4. Re-angle perturbation move in `pk/moves.py` + an explorer arm (section above); compare to existing arms from the same parents.
5. Theory: two-bulb waste accounting + a first gadget table for b ≤ 12; c(k) table vs R1′, Fang, the store, k ≤ 200.
