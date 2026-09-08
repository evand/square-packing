# The `t = 4` screen: corner branch + level-2 slots + anchor cliques on one packing-side LP

Code: `search/t4screen.py` (new; `search/level2_lp.py`'s region LP with the anchor-clique family of
`search/anchorclique.py` / `search/anchorsep.py` bolted on), `search/t4screen_loop.sh` (chunked
restart driver), `search/t4screen_phase2.sh` (the run set), `search/t4screen_table.py` (the table).
Nothing in `verify/`, `xcheck.py`, `lean/`, `certificates/` or `branch.py` is touched: this is a
packing-side screening instrument only.  Runs are in this worktree's `runs/t4_*` (gitignored, so
every number is quoted here).

Read against `notes/level2-design.md` (the branch), `search/CLIQUE_CEILING.md` and
`search/CLIQUE_CONTINUUM.md` (the cliques), `search/RECONCILE.md` (why the *interior* clique family
is the one that matters) and `search/DUAL_EXACT.md` / `search/CLOSED4.md` (the closed semantics).

---

## 1. What is measured

The LP is, over closed unit squares `S` contained in the closed container `[0,t]^2`:

```
max  sum_S mu_S
s.t. coverage(p) = sum_S mu_S [p in S]        <= 1     for every point p of [0,t]^2
     mu({S : centre(S) in C_i}) = kc_i                 i = 0..3   (the four corner boxes [0,r]^2)
     mu({S : centre(S) in W_j}) = kw_j                 j = 0..7   (the eight wall slots)
     mu(K)        = sum_{S in K} mu_S         <= 1     for every anchor clique K = K(p, A)
     mu >= 0
```

This is exactly the LP dual of the cover-side leaf value of a branch certificate with those regions
and per-region multipliers, strengthened by anchor-clique columns — i.e. of
`packing_le_weight_regions_cliques`.  **A leaf closes iff this value is `< 12`.**  The measure is
*not* D4-symmetrised (an asymmetric slot pattern needs the whole container), and every image of a
separated clique is its own `<= 1` row, which is strictly stronger than the orbit-averaged clique
column a symmetric model can carry.

Three things had never been in one LP before: `level2_lp.py` has the region equalities and no
cliques; `clique_ceiling.py` has cliques but no region equalities and is D4-symmetric;
`branch.py --cq-interior` has both but is the cover side and has never been run at `t = 4`.

### 1.1 Semantics — at `t = 4` this is the whole game

* **Container closed, squares closed.**  `S(c, th)` is admissible iff
  `w(th)/2 <= c_x, c_y <= t - w(th)/2` with `w(th) = |cos th| + |sin th|`, so a square may touch a
  wall.  In the code the centre is clamped with `packing_dual.ADM = 1e-9` of margin, so the
  wall-touching poses are represented to within `1e-9`; the effect on any value below is `< 1e-8`.
* **Points on the boundary of a square count as covered.**  The incidence test is
  `|R(-th)(p - c)|_inf <= 1/2 + 1e-8` (`packing_dual.TOL`), i.e. over-inclusive by `1e-8` in the
  square's own frame.  This is the convention of `certificates/FORMAT.md`, `verify/`,
  `search/DUAL_EXACT.md` and `search/CLOSED4.md`.  It is what makes the problem non-trivial at
  `t = 4`: the sixteen squares of the `4 x 4` grid all contain the grid crossing `(2,2)`, so a
  measure cannot put mass 1 on each of them, and it is why `nu_f(4) = 12.163` rather than 16.
* **Points on grid lines.**  A point on `x = 1` interior to two grid squares is in both closed
  squares, coverage 2.  Nothing special is done for them; they are ordinary arrangement vertices
  and the certification below evaluates them.
* **Coverage is certified over the continuum, not sampled.**  `Model.certify` (`packing_dual.py`'s
  `vertex_cov`) evaluates the coverage of the current measure at *every* vertex of the arrangement
  of its support squares — every square corner, every crossing of two square edges that lies on
  both boundaries, and the four container corners — which is the complete set of points where the
  piecewise-constant coverage attains its maximum.  The reported `M` is that maximum (in floats).
* **Cliques are pairwise-closed-intersecting.**  Touching counts, as
  `search/CLIQUE_CEILING.md` "Semantics" requires: a certificate at `t` refutes packings at
  `t' < t`, which rescale to pairwise *disjoint closed* squares in `[0,t]^2`.  Membership of a pose
  in `K(p, A) = {S : p in S, S meets A} u {S : A subset S}` is decided by
  `anchorclique.member` at half-side `h = 1/2` — the exact test for a closed unit square, three-axis
  SAT for the segment (the missing-axis bug of `notes/clique-family.md` 7 is not present).  It
  drops poses within `1e-9` of a piece boundary, which on the packing side under-counts members,
  i.e. imposes `mu(K') <= 1` for a slightly smaller `K' subset K` — a valid but marginally weaker
  constraint.
* **Region membership is by the CENTRE of the pose**, by `level2_regions.classify` with `r = 1`, so
  the boundaries in the centre plane sit at `x, y in {0,1,2,3,4}` and the sixteen grid poses at
  `(0.5+i, 0.5+j)` are `0.5` from every boundary (`notes/level2-design.md` 9, "put region
  boundaries at integers").  `classify` is a deterministic tie-break on the closed boxes, so
  nothing can be double-counted; a *certificate* would need the verifier's cell-vs-box test to
  agree with that tie-break, which is a real (unaddressed) obligation, not a formality.
* **Angles live on a lattice.**  `packing_dual.ANG_UNIT = 1e-5 rad`.  That is a restriction of the
  pose set, hence in the safe (value-lowering) direction.

### 1.2 Which way every error goes

| ingredient | restriction | effect on the value | handled by |
|---|---|---|---|
| finite pose set (seed grid, warm starts, pricing) | columns | **lowers** — optimistic | column pricing; the reported pricing gap |
| finite point rows | constraints | **raises** | row generation until the certified `M <= 1 + 1e-9` |
| finite clique set (heuristic separation) | constraints | **raises** | separation until the best separated `mu(K) <= 1 + 1e-6` |
| `TOL = 1e-8` over-inclusive coverage | constraints | lowers by `<= 1e-8` | — |
| `1e-9` under-inclusive clique membership | constraints | raises by a comparable amount | — |

So **`LP` is only a value the poses actually attain when `M <= 1` and `kmax <= 1`**; until then it
is an upper bound on the restricted-pose value and a bound on nothing else.  When both converge,
`LP` is a genuine **lower** bound on the leaf's true value over all poses, and

> **a leaf whose converged value is `>= 12` does not close — a decisive no-go for this family.**

`L = LP / max(M, kmax, 1)` is printed as a diagnostic only: rescaling a measure restores coverage
and clique feasibility but breaks the region *equalities*, so `L` is **not** a leaf value.  Read
`LP` with `M` and `kmax` beside it.

---

*(sections 2-6: results, table, verdict — filled in below)*
