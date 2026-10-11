# Strip leaf for s(12) = 4 (depth-first)

Decision (Evan, 10-11): depth first.  Prove one narrow leaf, with explicit hypotheses if needed, then
back up to the other branches.  Math and measurements: `notes/strip-split.md`.  Float scripts so far
were in the session scratchpad; copy any we keep into `search/strip/`.

**The leaf.**  12 closed unit squares in `[0,t]^2`, `t <= 4`.  Four of them form a wall strip in class A
(tilted middle: flat bottom-left, two middle squares tilted `>= θ0`, flat top-right).  The other
eight avoid a border `B` (profiled seam plus the top block) and cannot fit beside it.  Lean
skeleton: `lean/Sqpack/StripLeaf.lean` (`SeamOK`, `RegionOK`, `strip_leaf`).  A, not B/F, because
it is the narrowest class with the cleanest top block.

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
Current: 3 (region certificate at `x0 = 51/50`), with 4 in parallel.
2. (dropped 10-11: the slab form suffices for the narrow leaf.)  Lean: generalise `SeamOK`/`RegionOK` from a slab to an arbitrary border set `B` (profiled), and restate `strip_leaf`.  Small.
3. Exact `RegionOK` for the profiled region: cover of `[0,t]^2 \ B` with weight < 8, or counting from `StripRect`.  Main certificate.
4. Exact `SeamOK` for class A with `θ >= θ0`: one square against four fixed ones, over boxes of strip parameters.
5. Define the class-A region (pose boxes, tilt floor) as the leaf's hypothesis.  This gives a conditional leaf: "strip in class A ⇒ False".
6. Lean replays of 3 and 4.  Closes the leaf.

## Later (breadth)
- Classes B, F (same border shape, top square tilted); C, D, E and the `w = 2` plateau go to the near-axis/core branch.
- Class cover: every 4-strip of width < 2 is in A ∪ B ∪ F or near-axis (certificate).  Routing: which packings have a 4-strip at all.
- Decide: commit the Lean files (`Strip117` adds about 7 min of kernel checking per full build); ask Evan.
