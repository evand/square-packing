# Strip leaf for s(12) = 4 (depth-first)

Decision (Evan, 10-11): depth first.  Prove one narrow leaf, with explicit hypotheses if needed, then
back up to the other branches.  Math and measurements: `notes/strip-split.md`.  Float scripts: `search/strip/`.

**The leaf.**  12 closed unit squares in `[0,t]^2`, `t <= 4`.  Four of them form a wall strip in class A
(tilted middle: flat bottom-left, two middle squares tilted `>= θ0`, flat top-right).  The other
eight avoid a border `B` (profiled seam plus the top block) and cannot fit beside it.  Lean
skeleton: `lean/Sqpack/StripLeaf.lean` (`SeamOK`, `RegionOK`, `strip_leaf`).  A, not B/F, because
it is the narrowest class with the cleanest top block.

## Status (paused 10-11)
Paused at a natural point.  In hand: the Lean strip lemma and guillotine (unconditional, h = 1.17), the line sweep,
and the stitch `strip_leaf`.  The class-A leaf is pinned down in float (w <= 1.90, seam x0 = 51/50, margin about 0.02).
Still open: the region lemma (item 3).  A point cover can't do it (about 0.03 margin at best, see below), so it needs
branching plus count rows.  Honest read: the obstruction is again the near-tiling / sliding family, so this route
is likely to settle the side branches and hand the core back, not close s(12) by itself.
Resume options, cheapest first: (i) the float check "tilted square crossing both lines" (Ideas); (ii) tabulate
near-optimal configurations by line counts to size the case tree; (iii) item 3 with corner/notch-slot branching plus
wall-band count rows (needs about +0.2 over the cover LP).  Files: notes `notes/strip-split.md` (§4e = the LP),
scripts `search/strip/` (README), data `runs/strip_float/`, `runs/strip_region/` (local, gitignored).

## Done
- Lean `StripRect.lean`: strip lemma in rectangle and slab forms, R4/R7 widths, guillotine (given h >= 7/6), line sweep.
- Lean `StripChain/StripCheck/Strip117(.Data).lean`: `wallStripBound_117`, `guillotine_false_117` unconditional (7 min kernel check).
- Lean `StripLeaf.lean`: the stitch, given `SeamOK` and `RegionOK`.
- Float: strip census (classes A–F), region margin (top block about 0.3), seam margins (§4a–4d of the note).

## Now
1. **Fix the leaf's hypotheses (float, about 1 h).**  (a) Middle-row seam reach for class A with tilt `>= θ0`,
   `θ0` in {30°, 35°}, at `H = 3.9999`: does the margin plateau, or go to 0 as `b -> 0`?  (b) Region lemma with
   the profiled seam (`x0 = 1` in the bottom unit, `1 + m` in the middle, top block shifted back by `δ`): does
   it hold, and with what slack?
   Why first: the certificate statements (items 2–4) are functions of `θ0`, `m` and `δ`.  If (a) goes to 0,
   or (b) fails, the leaf needs a different border, and work on the certificates would be thrown away.
   Results 10-11: (a) plateaus.  Middle margin 0.036 at θ0 = 30° and 0.057 at 35°, the same at
   `H = 3.999` and `3.9999`.  (b) The bottom seam has no give: lowered to 0.95, 8 fit even with
   `m = 0.03`, because the bottom row then holds three squares beside the strip's flat bottom square
   (a margin-0 row).  Fix: give the region lemma the strip-lemma counts as a hypothesis, which Lean
   already has (`wallStripBound_117`): at most 2 of the 8 have their centre within 1.17 of the
   bottom, and likewise the top, since the strip's flat end squares use one slot each.
   That did not help.  With slack strips, 8 still fit: a 2 x 2 middle block tilted 12°, the bottom
   seam at 0.994, and the top square shifted back 0.15.  Both sides bind at the same flat contacts.
   **Resolution: narrow the leaf to tight strips.**  Class A with strip width `w <= 1.90` (minimum
   1.885).  Seam reach anywhere: 1.023 at `H = 3.999` and 1.022 at `3.9999` (plateau; middle tilts
   forced to about 38°).  Region with a uniform seam `x0 = 1.03`, top square at the corner or shifted
   0.015: no fit (300 starts); with `x0 = 1` it already fails for `δ <= 0.3`.  So the margin is
   about 0.02 at a uniform `x0 ≈ 1.02`.  No profile and no tilt floor are needed, and
   `StripLeaf.lean`'s slab form applies as is (`x0 = 51/50`).
   **Leaf (fixed 10-11):** class A, `w <= 1.90`, `t` in `[3.99, 4]`, `SeamOK t (51/50)`,
   `RegionOK t (51/50) top`.  Item 2 below is no longer needed.

## Next
Current when resumed: 3 (region lemma at `x0 = 51/50`; a cover alone fails, see the result below), with 4 in parallel.
2. (dropped 10-11: the slab form suffices for the narrow leaf.)  Lean: generalise `SeamOK`/`RegionOK` from a slab to an arbitrary border set `B` (profiled), and restate `strip_leaf`.  Small.
3. Exact `RegionOK` for the profiled region: cover of `[0,t]^2 \ B` with weight < 8, or counting from `StripRect`.  Main certificate.
   Result 10-11 (region-LP agent, `search/strip/region_lp.py`, `region_dual.py`, `runs/strip_region/`): **point covers
   are out.**  The fractional packing number (float lower bound on any closed cover) is 8.136 at x0 = 1.00 with notch
   0.85, and 8.03 at x0 = 1.02 with notch 0.60.  In the main case (x0 = 1.02, notch to 1.87) it is about 7.97 and still
   rising, so it is at most 0.03 below 8.  Open semantics fails trivially (eight flat squares in [2,4]x[0,4]).  The
   obstruction is a smeared near-tiling: about 0.7-1.0 on each free corner, about 0.9 in the top slot next to the notch,
   about 2.7 in flat squares sliding along the outer columns, and about 1.5 in a tilted middle block.
   Next: branch on corner / notch-slot occupancy and add wall-band count rows (strip lemma); these need to buy about 0.2.
4. Exact `SeamOK` for class A with `θ >= θ0`: one square against four fixed ones, over boxes of strip parameters.
5. Define the class-A region (pose boxes, tilt floor) as the leaf's hypothesis.  This gives a conditional leaf: "strip in class A ⇒ False".
6. Lean replays of 3 and 4.  Closes the leaf.

## Ideas worth exploring (routing; Evan 10-11)
- **Line sweep as the root of the case tree.**  Lean `line_sweep` + `wallStripBound_117`: every packing has at least 3
  squares crossing each line `x, y ∈ [1.66, 1.67]` (and the mirror lines near 2.335).  Branch on the crossers'
  tilts: some crosser tilted `>= θ0` → tilted side (strip leaves); all near-flat → near-tiling core.  Tight at the
  permutation-hole tilings (3 + 3 + 6).
- **Two lines** `u1 ≈ 1.665`, `u2 ≈ 2.335`: counts (left, crosses u1 only, crosses both, crosses u2 only, right,
  between).  The middle slab is narrower than 1, so it holds nothing wholly.  Left and right hold ≤ 3 each, and
  [0,u2] and [u1,t] hold ≤ 6 each.  The tilings give (3,3,0,3,3).  Correction (10-11): "crosses both ≥ n" is
  **not** a leaf.  Rows with an empty cell slide at zero cost, so the margin-0 family realises 0–4 squares
  crossing both lines (≤ 4 by chords on x = 2).  Counts alone never separate the sliding family; tilt does.
  Candidate leaf: "some square crossing both lines is tilted ≥ θ0".  Float check: minimum side for 12 squares
  with such a square, θ0 ∈ {15°, 25°, 35°}.
- **Centre-based split** (centre left/right of a line, or in wall bands): a partition with no "crossing" class,
  and the strip lemma's native form.  Known from the chat-lemma census: by itself it partitions the zero set
  (18 classes) without shrinking it.  Use it as extra counting constraints inside the containment split,
  not as the root.
- First step for all three: tabulate the near-optimal 12-configurations (census witnesses, today's float hits)
  by these counts; size the case tree before proving anything.

## Later (breadth)
- Classes B, F (same border shape, top square tilted); C, D, E and the `w = 2` plateau go to the near-axis/core branch.
- Class cover: every 4-strip of width < 2 is in A ∪ B ∪ F or near-axis (certificate).  Routing: which packings have a 4-strip at all.
