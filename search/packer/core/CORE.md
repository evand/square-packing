# Side-pinning cores: records vs their exact-side alternates (10-18)

Question (idea 2 of the 10-18 review): are records "a rigid core + loose filler", with the 123 exact-side alternates of
the shove runs (`runs/shove/ties1018`, 46 n) rearranging only the filler?  `coremap.py`, ~0 CPU (LPs, seconds).

* **Dual cores (one LP vertex, `--mode dual`; core1)**: 24-34 squares at n 200-300, ~2s (about one wall-to-wall chain
  per direction).  Record vs alternate overlap is low (J 0.02-0.3) even when nearly every square is pose-identical
  (233 / _2: 24 of 25 record core squares identical in the alternate, 5 carry force there): the force solution is
  degenerate, so a single vertex is not a fingerprint.
* **Locked sets (`fixlp.mobility`: no first-order motion at the fixed side over 12 random-objective LPs; = the union of
  all force supports; core2)**: median 57 % of a record's squares are locked (alternates 50 %).  Between record and
  alternate (best of 8 symmetries + Hungarian), median 49 % of all squares change pose but only 16.5 % of the record's
  locked squares do (median 20 squares); Jaccard of locked-and-identical median 0.68 (quartiles 0.39-0.82).  Lowest
  overlap where the alternate is a different arrangement altogether (257, 292, 170: md 0.15-0.27).
* Reading: yes, mostly.  Equal-side alternates rearrange the loose half; the locked half (the side-pinning skeleton) is
  largely shared, but not fixed (~20 locked squares per alternate do move, probably as rigid blocks).  For idea 2 this
  only halves the target (locked ~ half the packing), and the minimal force chains (~2s) are not unique.  The test that
  would decide whether skeleton-targeted moves are worth building: kicks restricted to locked vs loose vs all squares
  from the record, same quench, compare how often the side changes (not ties) and how far.
