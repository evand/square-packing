# s11-anatomy: what is the mass-11 fractional packing at `t ≈ 3.82–3.85`?  (2026-09-13)

**Why.**  `search/N11.md`: the pure cover LP for `n = 11` crosses 11 at `t* ≈ 3.816`, `0.06` below Trump's
conjectured-optimal `3.877083`; above the crossing the LP settles at exactly 11 with a dual packing measure
of mass 11.  The same rotation-carried gap as `n = 12`, twice as wide, and entirely at `t < 4` where the
sweep verifier, `branch.py --n 11`, cliques and polygons all already work.  Before pointing any of them at
it we need to know what the mass-11 packing *is*, because that decides the branching variable.

**Question.**  At `t ∈ {3.82, 3.85}` (and `3.87` if cheap), extract the exact fractional packing measure the
way `search/dual_exact.py` / `search/DUAL_EXACT.md` did for `n = 12` at `3.99`, and describe it: mass per
corner box, per wall strip, angle histogram, and the largest pairwise-disjoint subfamily of its support
(`search/rankdiag.py --step 1`).  Is it (a) Trump's tilted 11-packing smeared, (b) the `3×3` grid plus two
in disguise, or (c) something else?  For (a) the corner boxes will not be empty and the `n = 12` corner
branch is the wrong lever; say what the right one looks like.

**Semantics.**  As `N11.md`: container `[0,t]²`, closed unit squares, `D4`-symmetric point sets, exact
verifier with `n = 11`.  Seeds: `runs/n11_seed_union2.txt`, `runs/n11_H3985.txt` (see `runs/archive/` if
not in `runs/`).  2–4 threads, detached; nothing longer than ~3 h.

**Deliverable.**  `search/N11_ANATOMY.md`: the measures (exact, in `runs/`), the three-way answer with numbers,
and one paragraph on which `n = 12` instrument to apply next and what it should be expected to buy.
Do not edit `TODO.md`; do not touch `verify/`, `certificates/`, `lean/`.
