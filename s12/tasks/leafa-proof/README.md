# leafa-proof: replacement hulls for the eight singletons, then a written proof that leaf A admits no closed packing  (drafted 2026-09-13, not yet launched)

**Why.**  `search/RANK8.md`: leaf A (four corner squares holding `{a_i,b_i}`, eight squares holding one of
`c_j, d_j` each) has margin exactly `0` on a plateau reaching `32°` of tilt, no sub-core (every eleven-square
sub-family is realisable), and first-order rigidity carried by wall-to-wall chains of four squares.  Every
degree-1 family sits at `12` on it (`BENTZ.md` §0).  So its death is a pairwise statement, and a proof of it is
the `m = 4` analogue of Bentz's Corollary 7 / Nagamochi's chain argument: Bentz's own tools, one level down.

**Step 1 — replacement hulls (machine, exact).**  Stromquist's trick (Bentz Lemma 10, `notes/proof-anatomy.md`
§2.3): if a square contains `p ∈ P0` and no other point of `P0`, then it contains every `p′` for which
`(P0 ∖ {p}) ∪ {p′}` is still a closed cover of `[0,4]²`.  In leaf A, `K + u = 4` forces the eight
wall/interior squares to be singletons, so the lemma applies to each `c_j`, `d_j` (and to the corner squares
via their pairs).  Compute, for each of the eight points, the region of admissible `p′` — numerically first
(`CLOSED4.md`-style row generation), then exactly certified where possible (`verify2/zmcheck` on the moved
16-point set; note `zeromargin.py` left 3,516 boxes uncertified on the *unmoved* set, `BENTZ.md` §3, so exact
certification of 16-point closed covers at `t = 4` is a prerequisite and may need a new primitive) — and report
the forced polygons `H_p` with exact rational vertices.

**Known in advance (do not look for it):** no two hulls share a point.  The lemma is degree-1 (a single-square
pose-region inclusion), and the certified leaf-A measure has coverage `≤ 1` everywhere with mass `≈ 1` on every
pattern region, so a common point of `H_c` and `H_d` would give coverage `≈ 2`.  The hulls cannot kill the leaf
by double coverage; their job is bookkeeping for step 2.

**Step 2 — the written proof (human-shaped; an agent may draft it, the owner reviews).**  With the twelve
polygons in hand, the statement is: twelve closed unit squares, the `i`-th containing `H_i`, pairwise disjoint
as closed sets in `[0,4]²` — impossible.  Tools: the corner lemma (a corner square reaches only its pair),
projection widths of a square containing a given polygon, Bentz's Corollary 7 wall-intersection bookkeeping,
the `"> 1"` strictness lemmas, the one-point `s(4) = 2` argument (any unit square in `[0,2]²` contains its
centre), and the chain lemma (four unit squares chained wall-to-wall in width `4` must be axis-parallel and
touching).  `RANK8.md` §3 says which chains carry the obstruction near the tiling (`T2,D1,D2,C3` along `x`;
`C1,D0,D1,C2` and `T1,C3,D2,D3,T3` along `y`) and §1 says the touching pair moves across the plateau, so the
argument must be combinatorial (a chain always exists), not perturbative.  Every lemma stated with its
hypotheses; every numerical inequality checked exactly; every case that is not closed listed as open.

**Deliverables.**  `search/HULLS.md` + code for step 1; `notes/leafA-proof.md` for step 2 (a proof, or a
precise list of the cases that do not close and why).  Semantics as `notes/s13-casefree.md` §1.  Do not touch
`verify*/`, `certificates/`, `lean/`.

**What it does and does not buy.**  Leaf A proven kills one leaf of a tree whose other leaves (`k = 3`
corners at `≥ 12.03`, the other 42 corner-pattern classes, the slot leaves) also do not close fractionally.
But it is the template: if the chain argument works with the hulls as bookkeeping, the other leaves have more
room, and the `tasks/corner-triage` enumeration says which ones need it.
