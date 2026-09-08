# The `t = 4` screen: corner branch + level-2 slots + anchor cliques on one packing-side LP

Code: `search/t4screen.py` (new; `search/level2_lp.py`'s region LP with the anchor-clique family of
`search/anchorclique.py` / `search/anchorsep.py` bolted on), `search/t4screen_loop.sh` (chunked
restart driver), `search/t4screen_phase2.sh` (the run set), `search/t4screen_table.py` (the table).
Nothing in `verify/`, `xcheck.py`, `lean/`, `certificates/` or `branch.py` is touched: this is a
packing-side screening instrument only.  Runs live in this worktree's `runs/t4_*` (gitignored, so
every number is quoted here).

Read against `notes/level2-design.md` (the branch), `search/CLIQUE_CEILING.md` and
`search/CLIQUE_CONTINUUM.md` (the cliques), `search/RECONCILE.md` (why the *interior* clique family
is the one that matters) and `search/DUAL_EXACT.md` / `search/CLOSED4.md` (the closed semantics).

---

## 0. Verdict, up front

**Undecided, leaning GO — with a new obstacle that has nothing to do with the value.**

Nothing in this screen reached 12.  At `t = 4.0`, closed semantics, on pose sets of 3.4k–6.6k poses
with 8k–34k certified coverage rows and 250–2100 anchor-clique rows: the pure relaxation sits at
`12.12–12.26`, the corner leaf `k = 4` at `11.95–11.98` (its *pure* value landing on exactly
`12.000000` in three independent runs), and **all four `m = 4` slot leaves at `10.3–11.4`, i.e.
`0.6–1.7` below 12**, two of them at a fully converged `11.000000` (certified coverage `M = 1`
exactly and no anchor clique the separator can violate).  So no leaf produced the decisive
`>= 12`, and the level-2 slot branch — not the cliques — is what does the work: the matched clique
gain at `t = 4` is `+0.02` to `+0.56`, while the slot branch is worth `0.6–1.7` on top of the
corner branch.

That is a **weak** go, and it is weak for the reason `CLIQUE_CEILING.md` already recorded: every
number is an LP value on a restricted pose set with a large residual pricing gap, and the
inner loop does not converge (`M` and `kmax` still `1.02–1.3` in most rows), so most of the values
are not bounds in either direction.  The clean exceptions are the two `11.000000` rows.

**The new obstacle is that the leaf optima are not certificate-shaped.**  `7–9 %` of the mass of
every `m = 4` leaf optimum — and in two of the four leaves a *whole unit* of it, at a single pose —
sits **exactly on a slot boundary**, `c_x` or `c_y` equal to `t/2 = 2` to within `1e-14`.  The
level-2 branch is a partition only through a tie-break at that boundary, and no cell-based verifier
can implement a `1e-14` tie-break.  `notes/level2-design.md` 9 checked that the *cover's tight
poses* (the `4 x 4` grid centres) are `0.5` from every boundary; the *packing side's* optimum goes
straight to the boundary instead.  Any follow-up has to deal with this before it can deal with the
value.

---

## 1. What is measured

Over closed unit squares `S` contained in the closed container `[0,t]^2`:

```
max  sum_S mu_S
s.t. coverage(p) = sum_S mu_S [p in S]        <= 1     for every point p of [0,t]^2
     mu({S : centre(S) in C_i}) = kc_i                 i = 0..3   (the corner boxes [0,r]^2, r = 1)
     mu({S : centre(S) in W_j}) = kw_j                 j = 0..7   (the eight wall slots)
     mu(K)        = sum_{S in K} mu_S         <= 1     for every anchor clique K = K(p, A)
     mu >= 0
```

By LP duality this is exactly the cover-side leaf value of a branch certificate with those regions
and per-region multipliers, strengthened by anchor-clique columns — the object
`packing_le_weight_regions_cliques` reduces.  **A leaf closes iff this value is `< 12`.**

The measure is *not* D4-symmetrised (an asymmetric slot pattern needs the whole container), and
every image of a separated clique is its own `<= 1` row, which is strictly stronger than the
orbit-averaged clique column a symmetric model carries.

Three ingredients had never been in one LP: `level2_lp.py` has the region equalities and no
cliques; `clique_ceiling.py` has cliques (box cliques, on a D4-symmetric model) but no region
equalities; `branch.py --cq-interior` has both but is the cover side and has never been run at
`t = 4`.

### 1.1 Semantics — at `t = 4` this is the whole game

* **Container closed, squares closed.**  `S(c, th)` is admissible iff
  `w(th)/2 <= c_x, c_y <= t - w(th)/2`, `w(th) = |cos th| + |sin th|`, so a square may touch a wall.
  Centres are clamped with `packing_dual.ADM = 1e-9` of margin, so wall-touching poses are
  represented to within `1e-9`; the effect on any value below is `< 1e-8`.
* **A point on the boundary of a square counts as covered.**  The incidence test is
  `|R(-th)(p - c)|_inf <= 1/2 + 1e-8` (`packing_dual.TOL`), over-inclusive by `1e-8` in the square's
  own frame.  This is the convention of `certificates/FORMAT.md`, `verify/`, `search/DUAL_EXACT.md`
  and `search/CLOSED4.md`, and it is what makes `t = 4` non-trivial: the sixteen squares of the
  `4 x 4` grid all contain the crossing `(2,2)`, so a measure cannot put mass 1 on each — which is
  why `nu_f(4) = 12.163` and not 16.
* **Points on grid lines** are ordinary arrangement vertices; a point on `x = 1` between two grid
  squares lies in both, coverage 2.  Nothing special is done for them.
* **Coverage is certified over the continuum, not sampled.**  `certify` (`packing_dual.py`'s
  `vertex_cov`) evaluates the coverage of the current measure at *every* vertex of the arrangement
  of its support squares — every square corner, every crossing of two square edges lying on both
  boundaries, and the four container corners — the complete set of points at which the
  piecewise-constant coverage attains its maximum.  `M` is that maximum (in floats).
* **Cliques are pairwise-closed-intersecting** (touching counts), as `CLIQUE_CEILING.md` requires: a
  certificate at `t` refutes packings at `t' < t`, which rescale to pairwise *disjoint closed*
  squares in `[0,t]^2`.  Membership in `K(p, A) = {S : p in S, S meets A} u {S : A subset S}` is
  decided by `anchorclique.member` at half-side `h = 1/2` — the exact test for a closed unit square,
  three-axis SAT for the segment (the missing-axis bug of `notes/clique-family.md` 7 is not
  present).  It drops poses within `1e-9` of a piece boundary, which on the packing side
  under-counts members: it imposes `mu(K') <= 1` for a slightly smaller `K' subset K`, a valid but
  marginally weaker constraint.
* **Region membership is by the CENTRE of the pose**, via `level2_regions.classify` with `r = 1`, so
  the boundaries in the centre plane are at `x, y in {0,1,2,3,4}` and the grid poses `(0.5+i, 0.5+j)`
  are `0.5` from every boundary.  `classify` is a deterministic tie-break on closed boxes, so it is
  a partition; §5 is about what happens when the optimum sits on the tie-break.
* **Angles live on a lattice**, `packing_dual.ANG_UNIT = 1e-5 rad`.  A restriction of the pose set,
  hence in the value-lowering direction.

### 1.2 Which way every error goes

| ingredient | is a | effect on the LP value | handled by |
|---|---|---|---|
| finite pose set (seed grid, warm starts, pricing) | restriction of the columns | **lowers** — optimistic | pricing; the reported gap |
| finite point rows | relaxation of the constraints | **raises** | row generation until certified `M <= 1 + 1e-9` |
| finite clique set (heuristic separation) | relaxation | **raises** | separation until the best `mu(K) <= 1 + 1e-6` |
| `TOL = 1e-8` over-inclusive coverage | tightening | lowers by `<= 1e-8` | — |
| `1e-9` under-inclusive clique membership | relaxation | raises comparably | — |

So **`LP` is a value the poses actually attain only when `M <= 1` and `kmax <= 1`**; until then it
is an upper bound on the restricted-pose value and a bound on nothing else.  When both converge,
`LP` is a genuine **lower** bound on the leaf's true value over all poses, and

> **a leaf whose converged value is `>= 12` does not close — a decisive no-go for this family;**
> a value well below 12 is the weaker, optimistic direction.

`L = LP / max(M, kmax, 1)` is printed as a diagnostic only: rescaling restores coverage and clique
feasibility but breaks the region *equalities*, so `L` is **not** a leaf value.  Read `LP` with `M`
and `kmax` beside it.

---

## 2. Calibration, done first

**(a) The regression against `notes/level2-design.md` 4.4.**  Same script family, same pose set
(the certified closed-container support `dual_PC1_support.txt` rescaled to `t = 4`, plus a
`0.25 / 15 deg` grid = 1471 poses), same settings:

| leaf | this run (`runs/t4_T0.log`, `runs/l2_REGc40.log`) | `level2-design.md` 4.4 |
|---|---|---|
| control (`k = 4`, slots free) | **11.896104** | **11.896104** |
| `01010101` | **11.385844** | **11.385844** |
| `01100110` | **11.000000** | **11.000000** |

Exact agreement, so the region bookkeeping and the LP are the same object as the published one.
(Note the control's rows do *not* converge — `M` oscillates in `1.008–1.06` — so `11.896104` is not
a bound; the two leaf values do converge, `M = 1.000000`.)

**(b) Weak duality against the cover side.**  The cover-side value for the corner leaf `k = 4` at
`t = 3.98` with interior anchor cliques is `11.958` (`runs/branch_J16i.log` in the main tree, it55,
`obj = 11.957911`, still creeping up).  Its dual is a fractional packing of `sigma`-shrunk squares,
so the packing-side value over *closed unit* squares for the same leaf and the same clique family
must come in `<= ~11.96`.  Warm-started from that dual and its 1342 cliques
(`runs/t4_C98.log`, 7835 poses, 26.8k rows):

```
clique LP 11.790666  (M 1.0225, kmax 1.0526)     matched pure 11.904472  (M 1.0494)   gain +0.1138
```

`11.79 <= 11.96` and `11.90 <= 11.96`: **the machinery passes the calibration**, with `0.05–0.17` of
headroom, which is about the size of the shrunk-square-versus-unit-square gap
(`level2-design.md` 4.1).  The matched clique gain on the packing side, `+0.114`, is the same order
as the cover side's `+0.053` on its own instance.

---

## 3. The table

Values at `t = 4.0`, closed semantics, `r = 1`.  `step` is the refinement stage (a pricing round;
the chunked driver restarts the counter, and `t4screen_table.py` renumbers).  `LP_pure` is the
**matched pair**: the same poses and the same rows with the clique rows switched off.  `converged`
says which of the two generation loops had actually closed when the value was taken — only
`rows cliq` rows are bounds.

*(TABLE)*

---

## 4. Reading

*(READING)*

---

## 5. The obstacle: the leaf optimum sits on the slot boundary

*(BOUNDARY)*

---

## 6. Reproduce

*(REPRODUCE)*
