# Credits and lineage

This result stands on other people's work.  Listing it precisely, since the field currently has
several parallel efforts that are largely unaware of each other.

## The classical method

The **unavoidable set of points** argument — exhibit `n-1` points such that every unit square in
the container contains one, therefore `n` interior-disjoint squares cannot fit — is due to

* **F. Göbel**, *Geometrical packing and covering problems*, in Packing and Covering in
  Combinatorics (A. Schrijver, ed.), Math. Centrum Tracts **106** (1979) 179–199.

and was developed by

* **W. Stromquist**, *Packing 10 or 11 unit squares in a square*, Electron. J. Combin. **10**
  (2003) #R8 — the source of the bound this work improves, `s(11) ≥ 2 + 4/√5 = 3.788854…`,
  which transfers to `s(12)` by monotonicity;
* **M. J. Kearney and P. Shiu**, *Efficient packing of unit squares in a square*,
  Electron. J. Combin. **9** (2002) #R14 — the green/red duality trick;
* **H. Nagamochi**, *Packing unit squares in a rectangle*, Electron. J. Combin. **12** (2005)
  #R37 — weighted points, segments and areas as "resources";
* **W. Bentz**, *Optimal packings of 13 and 46 unit squares in a square*, Electron. J. Combin.
  **17** (2010) #R126, and *Optimal packings of 22 and 33 unit squares in a square*,
  arXiv:1606.03746 — continuously varying families of unavoidable sets.  `s(13) = 4` is what
  makes `n = 12` the exact boundary case.
* **E. Friedman**, *Packing unit squares in squares: a survey and new results*, Electron. J.
  Combin. Dynamic Survey **DS7** — the standard reference, and the source of the technical
  lemmas the whole area uses.

Best known packings, including `s(11) = 3.877083…` (**W. Trump**, 1979) and the current record
tables, are maintained by **David Ellsworth** at
https://kingbird.myphotos.cc/packing/squares_in_squares.html (continuing Friedman's pages).

## The mechanised lower-bound certificates (2026)

Two independent families appeared in August 2026, both computer-generated.  This work is a direct
descendant of the second.

**Weighted / LP family — what this repo builds on:**

* **Sam Burns**, *Proposing a better lower bound for n=17 square packing*, 6 Aug 2026,
  https://sam-burns.com/posts/proposing-better-lower-bound-for-n17-square-packing/ —
  weighted atoms with a total mass below `n`, verified in exact rational arithmetic.
* **Gustavo Massaccesi**, *Another better lower bound for n=17 square packing*, 21 Aug 2026,
  http://gus-massa.blogspot.com/2026/08/another-better-lower-bound-for-n17.html , and the
  methodology post https://gus-massa.blogspot.com/2026/08/linear-programing-for-square-packing.html —
  weights obtained by **linear programming** over a symmetry-reduced grid, a rational angle net
  from a tan-half-angle parameterisation, and the covering inequality `cos ε + sin ε ≤ 1 + ε`.

The idea of replacing an unavoidable *set* by LP-optimised *weights*, and of verifying the
covering condition in exact rational arithmetic over a rational angle net, is theirs.  This repo
applies it to `n = 12` and adds cell-based cutting planes, column generation, an arrangement-sweep
verifier, a Lean formalisation of the reduction, and scaling a finished certificate to criticality.

**Unweighted / subdivision family — different architecture, same problem:**

* **Stanislav Fort**, https://github.com/stanislavfort/17squares — `s(17) > 4.456575`.
* **Mira**, https://github.com/Mira-acc/17squares — `s(17) > 4.468292`; 16 unweighted points, a
  122.6M-node exact dyadic subdivision of pose space, a strict triangle-piercing lemma, and three
  independent exact checkers.  The most thorough verification artifact in this area, and the
  source of several practices adopted here (rejection tests, a self-describing `points.json`,
  publishing checkers rather than only claims).

These use **open** squares (strict interior) with a pigeonhole argument; the weighted family uses
**closed** squares with a concentric-shrink argument.  That difference is mathematical, not
cosmetic: a certificate written for one convention does not verify under the other's checker
without reworking the disjointness step.

## Provenance of this work

Produced by Claude (Anthropic) in a single session under human direction, in the same spirit of
disclosure that Fort and Mira adopted for theirs.  The mathematics is checked, not vouched for by
authority: the point of shipping exact certificates, three checkers, rejection tests and a Lean
proof of the reduction is that nobody has to take anyone's word for it.

Nothing here has been peer reviewed.
