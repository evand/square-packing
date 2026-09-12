# Rung 2: a weighted closed cover of `[0,4]²` of total weight `< 13`

**Claim.**  `s13_closed_cover_4.txt` is a weighted point set `P ⊂ [0,4]²` of total weight

    2591194431/200000000  =  12.955972155  <  13

such that **every closed unit square contained in `[0,4]²`, at every centre and every angle,
captures total weight `≥ 1`** (closed containment: a point on the boundary of the square counts —
`certificates/FORMAT.md`, `search/ZEROMARGIN.md` §1).

**What it proves.**  By the reduction of `certificates/FORMAT.md` ("What the file asserts"), a
packing of `n` unit squares in a container of side `s' < 4` rescales to `n` interior-disjoint
closed unit squares in `[0,4]²`, each capturing weight `≥ 1` and no point counted twice, so
`n ≤ 12.955972 < 13`.  Hence **no 13 unit squares fit in a square of side `< 4`**, and since 16 of
them tile `[0,4]²`,

> **s(13) = 4.**

This is a *case-free machine proof*: one cover, one exhaustive verification, no branch tree.
(Bentz 2010 needs a 6-leaf case analysis; DS7's pure sets cost 14.)

**The file.**  `certificates/FORMAT.md` plain format: container side `4/1`, coordinate denominator
`D = 1000`, weight denominator `W = 10^9`, `3621` points, each line `X Y w` meaning the point
`(X/D, Y/D)` with weight `w/W`.  Every number is an integer; nothing here is floating point.

**How it was verified** (`search/ZEROMARGIN.md`, `search/RUNG2.md`):

```
python3 search/zeromargin.py cert certificates/rung2/s13_closed_cover_4.txt \
        --depth 18 --nproc 8 --disj --chain-from 0 --dump runs/leaves.txt
```

The checker subdivides pose space `(c_x, c_y, u = tan(θ/2))` into boxes with rational endpoints and
certifies each one exactly, in `fractions.Fraction`; floats appear only as pre-filters that can
lose a certification but never create one.  `VERIFIED` means **0 uncertified boxes**: every
admissible pose of the container lies in some leaf whose certificate is exact.

Result:

```
done in 6138s: boxes 16872, max depth 13
  leaves: ADM 2867  CORE 0  P1 0  MIX 0  CHAIN 5320  TRI 0  EMPTY 3449  UNCERTIFIED 0
VERIFIED
```

1 h 42 min on 8 processes; the depth limit is never reached (max depth 13 of 18), so the
subdivision terminated on its own.  Of the 8,187 non-empty leaves, `ADM` carries 2,867 and `CHAIN`
— the disjunctive primitive — carries 5,320 (65 %).  `CORE`, `P1` and `MIX` are 0 because `ADM`
subsumes them, and `TRI` is 0 because a weighted cover has no unit-weight triangles.  `CHAIN`
dominating is forced: `search/RUNG2.md` Theorem 1 proves that **no** cover of total weight below
`16` can be certified by fixed-witness primitives at any depth, so a cover of weight 12.96 must be
certified disjunctively almost everywhere.  Independent checks:
`search/zeromargin_stress.py` re-verifies every leaf's witness from the leaf dump in floats with no
shared code path — 40 sampled admissible poses on each of the 11,636 leaves, **0 failures**, plus
`200,000` random `P1` instances, `68,545` core-lemma, `72,175` triangle, `2,713` `ADM` instances and
`220,000` `CHAIN` polynomial/enclosure instances, all 0 failures.  And `python3 search/closed4.py
stress` scans the cover densely, independently of the box machinery, reporting a minimum captured
weight of `1.0313820` (a 3.1 % margin; honest cost `total/min = 12.5618`).

**Lower bound for context.**  `search/COVER4.md` proves `COVER^closed(4) ≥ 24537607710/1999999999
= 12.2688038611` exactly, so this cover is within 5.6 % of optimal and no cover argument can reach
below `12.2688`.

**Provenance.**  The point set comes from `search/closed4.py`'s cover LP run with column generation
from the columns of `runs/closed4_best.txt`, plus the explicit row families of
`search/family_rows.py` (the `ZEROMARGIN.md` §4 item-3 poses that the LP's own row lattice steps
over — the omission that made every earlier candidate invalid), scaled by the exact rational factor
`21/20` with `search/scale_cover.py`.  See `search/RUNG2.md` §10.
