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

**UNDECIDED — no leaf reached 12, but the pose sets are demonstrably starved, so this is not
evidence for a go either.  Plus a new obstacle that has nothing to do with the value.**

Nothing in this screen reached 12.  At `t = 4.0`, closed semantics, on pose sets of 3.4k–6.6k poses
with 8k–34k certified coverage rows and 250–2100 anchor-clique rows: the pure relaxation sits at
`12.12–12.26`, the corner leaf `k = 4` at `11.95–11.98` (its *pure* value landing on exactly
`12.000000` in three independent runs), and **all four `m = 4` slot leaves between `10.29` and
`11.59`, i.e. `0.4–1.7` below 12**, two of them at a fully converged `11.000000` (certified coverage
`M = 1` exactly and no anchor clique the separator can violate).  So no leaf produced the decisive
`>= 12`, and the level-2 slot branch — not the cliques — is what does the work: the matched clique
gain at `t = 4` is `+0.02` to `+0.56`, while the slot branch is worth `0.4–1.7` on top of the corner
branch.

The closest any leaf came: **`01010101`, the hardest, at `LP = 11.591026` with cliques and
`11.837009` without them on the same 6047 poses and 10 140 rows, `M = 1.056` — `0.16` short of 12,
with the anchor cliques buying `0.25` of that.**  It was still rising when the budget ran out.

The arithmetic in an `m = 4` leaf is exact and explicit (§4.2): the branch pins `4` in the corner
boxes and `4` in the slots, so **`leaf value = 8 + interior mass`, and the leaf fails to close iff
the interior `[1,3]^2` can carry `>= 4`**.  Measured interior mass in the four leaves: `2.3 – 3.5`.

**And that is exactly where the screen shows its own limits.**  The interior sub-problem — the same
LP with the whole frame pinned to zero, so that every pose carrying mass is interior-centred — has
an independent witness, and the LP is short of it:

| interior mass at `t = 4` | value |
|---|---|
| LP, 3500 poses, one unrefined evaluation, `M = 1`, `kmax = 1` (**converged**) | `4.000000` |
| LP, 1471-pose seed, priced against its own dual, `M = 1`, `kmax = 1` (**converged**) | `5.000000` |
| independent witness: **six** unit squares with centres in `[1,3]^2`, pairwise disjoint with margin `+0.0567`, all inside `[0,4]^2` (`level2_capacity.py --t 4.0 --box 1 3 1 3 --k 6`) | `>= 6` |
| the same LP after 22 pricing stages (`INTP`), converged again at stages 17, 20, 21 | **`7.433677`, `7.508327`, `7.528038`** — and `7.583` at stage 22, still climbing `+0.05`/stage |

**A converged LP value on the seed pose set was `4.000000`; the same LP on the same problem is at
`7.53` after column generation and has not stopped.**  That is the calibration of this whole screen,
and it is brutal: `M <= 1` and `kmax <= 1` certify only that the measure is feasible, never that the
pose set was rich enough, and here the pose set was short by `3.5` on a quantity whose threshold is
`4`.  Every leaf number in §3 is a lower bound of exactly that kind.  With the deficit to 12 only
`0.4 – 1.0` and the leaf values still rising when the budget ran out (`01010101`:
`11.376 -> 11.591`, matched pure `11.837`), **this screen cannot tell a leaf that closes from one
that does not.**  What it *can* say is that nothing here produced the decisive `>= 12`, and that the
level-2 slot branch — not the cliques — is what moves the value.

**The new obstacle is that the leaf optima are not certificate-shaped.**  `1–11 %` of the mass of
every `m = 4` leaf optimum — and, in the two leaves whose value is actually converged, a **whole
unit** of it at a single pose — sits **exactly on a slot boundary**, `c_x` or `c_y` equal to
`t/2 = 2` to within `1e-14`.  The level-2 branch is a partition only through a tie-break at that
boundary, and no cell-based verifier can implement a `1e-14` tie-break.
`notes/level2-design.md` §9 checked that the *cover's* tight poses (the `4 x 4` grid centres) sit
`0.5` from every boundary; the *packing side's* optimum goes straight to the boundary instead.  A
follow-up has to settle this before the value means anything for a certificate.

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

```
tag         t     corn  pattern   cq  step         LP         M     kmax   LP_pure     gain  poses   rows  cliq  converged
B40         4.0   1111  ........  1   0     11.945452  1.106717  1.09922 12.000000 +0.054548   3383  19506  1798    .   .  
B40         4.0   1111  ........  1   1     11.982171  1.063609  1.13609 12.000000 +0.017829   6643  21491  2135    .   .  
C98         3.98  1111  ........  1   0     11.790666  1.022454  1.05264 11.904472 +0.113806   7835  26763   750    .   .  
INT40       4.0   0000  ........  1   0     10.000000  1.314039  1.08078 10.333333 +0.333333   3383   8718   460    .   .  
INT40       4.0   0000  ........  1   1     10.427699  1.162163  1.13569 10.559732 +0.132033   3500   9223   408    .   .  
INT40       4.0   0000  00000000  1   0      4.000000  1.000000  1.00000  4.000000 +0.000000   3500   9751   362  rows cliq
INT40       4.0   0000  ........  1   2     10.524676  1.200683  1.17459 10.698892 +0.174216   3500   9429   479    .   .  
INT40       4.0   0000  ........  1   3     10.707096  1.133619  1.11946 10.859229 +0.152132   6080  13131   500    .   .  
INT40       4.0   0000  ........  1   4     10.695177  1.043095  1.03796 10.876266 +0.181089   6080  13186   499    .   .  
INTP        4.0   0000  00000000  1   0      5.000000  1.000000  1.00000  5.000000 +0.000000   1471   7312   150  rows cliq
INTP        4.0   0000  00000000  1   1      5.730642  1.344856  1.15238  5.833333 +0.102692   4000   7581   206    .   .  
INTP        4.0   0000  00000000  1   2      6.000000  1.228571  1.17143  6.000000 -0.000000   4000   7608   467    .   .  
INTP        4.0   0000  00000000  1   3      6.189189  1.155405  1.15541  6.417969 +0.228780   4000   7810   479    .   .  
INTP        4.0   0000  00000000  1   4      6.393548  1.212903  1.27742  6.500000 +0.106452   4000   7691   234    .   .  
INTP        4.0   0000  00000000  1   5      6.580105  1.141081  1.03204  6.666667 +0.086561   4000   7924   172    .   .  
INTP        4.0   0000  00000000  1   6      6.654352  1.248354  1.18694  6.741935 +0.087583   4997   8232   370    .   .  
INTP        4.0   0000  00000000  1   7      6.778140  1.172920  1.16313  6.910389 +0.132249   4000   8276   371    .   .  
INTP        4.0   0000  00000000  1   8      6.825854  1.122486  1.09718  6.934089 +0.108234   4000   8416   399    .   .  
INTP        4.0   0000  00000000  1   9      6.911492  1.146670  1.04381  7.056314 +0.144822   4000   8458   174    .   .  
INTP        4.0   0000  00000000  1   10     6.955914  1.067662  1.07738  7.125739 +0.169825   4000   8570   242    .   .  
INTP        4.0   0000  00000000  1   11     7.086535  1.106948  1.14044  7.188158 +0.101623   4000   8921   226    .   .  
INTP        4.0   0000  00000000  1   12     7.099437  1.100939  1.09197  7.205677 +0.106241   4972   9206   306    .   .  
INTP        4.0   0000  00000000  1   13     7.228561  1.043921  1.01692  7.298233 +0.069672   4000   9213   169    .   .  
INTP        4.0   0000  00000000  1   14     7.314452  1.019972  1.03612  7.368554 +0.054102   4000   9921   129    .   .  
INTP        4.0   0000  00000000  1   15     7.323647  1.065391  1.04493  7.405091 +0.081444   4000  10627   127    .   .  
INTP        4.0   0000  00000000  1   16     7.386122  1.061973  1.05410  7.440140 +0.054019   4000  10997    85    .   .  
INTP        4.0   0000  00000000  1   17     7.433677  1.000000  1.02138  7.490635 +0.056958   4000  11927   256  rows  .  
INTP        4.0   0000  00000000  1   18     7.435097  1.077152  1.08192  7.492751 +0.057654   5014  11849   136    .   .  
INTP        4.0   0000  00000000  1   19     7.472754  1.015389  1.00721  7.508579 +0.035826   4000  12668   110    .   .  
INTP        4.0   0000  00000000  1   20     7.508327  1.000000  1.00000  7.541727 +0.033400   4000  13888   199  rows cliq
INTP        4.0   0000  00000000  1   21     7.528038  1.000000  1.00000  7.558828 +0.030790   4000  15269   111  rows cliq
INTP        4.0   0000  00000000  1   22     7.583111  1.045685  1.05356  7.603425 +0.020314   4000  21020   405    .   .  
K40         4.0   1111  ........  1   0     11.977689  1.160045  1.15082 12.000000 +0.022311   3383  11604   500    .   .  
K40         4.0   1111  ........  1   1     11.961620  1.285862  1.07017 12.000000 +0.038380   3500  12201   396    .   .  
K40         4.0   1111  ........  1   2     11.969832  1.086818  1.01944 12.000000 +0.030168   5516  16292   500    .   .  
K40         4.0   1111  ........  1   3     11.966539  1.092303  1.09743 12.000000 +0.033461   5516  14058   496    .   .  
K40         4.0   1111  ........  1   4     11.965357  1.082535  1.02590 12.000000 +0.034643   5516  15091   499    .   .  
K40         4.0   1111  ........  1   5     11.963270  1.197913  1.03814 12.000000 +0.036730   5516  14408   499    .   .  
L01010101   4.0   1111  01010101  1   0     11.376102  1.118124  1.15800 11.605192 +0.229090   3383   8727   340    .   .  
L01010101   4.0   1111  01010101  1   1     11.403693  1.075974  1.12692 11.678834 +0.275141   3500   8793   235    .   .  
L01010101   4.0   1111  01010101  1   2     11.410038  1.098586  1.04252 11.709679 +0.299641   5640   9402   413    .   .  
L01010101   4.0   1111  01010101  1   3     11.477424  1.183062  1.05936 11.689063 +0.211638   3500   9356   186    .   .  
L01010101   4.0   1111  01010101  1   4     11.472380  1.302078  1.10385 11.776351 +0.303971   6371   9276   293    .   .  
L01010101   4.0   1111  01010101  1   5     11.501333  1.236150  1.08579 11.725538 +0.224205   3500   9592   313    .   .  
L01010101   4.0   1111  01010101  1   6     11.520937  1.528112  1.05974 11.821764 +0.300827   5835  10094   372    .   .  
L01010101   4.0   1111  01010101  1   7     11.561721  1.533589  1.14435 11.777854 +0.216133   3500   9056   283    .   .  
L01010101   4.0   1111  01010101  1   8     11.560737  1.518808  1.06648 11.738728 +0.177991   3500   8750   212    .   .  
L01010101   4.0   1111  01010101  1   9     11.579788  1.279641  1.11443 11.730414 +0.150625   3500   9404   498    .   .  
L01010101   4.0   1111  01010101  1   10    11.591026  1.056408  1.08454 11.837009 +0.245983   6047  10140   388    .   .  
L01010101   4.0   1111  01010101  1   11    11.591807  1.467476  1.01509 11.808795 +0.216988   3500   9127   178    .   .  
L01010101   4.0   1111  01010101  1   12    11.605211  1.467106  1.08081 11.806248 +0.201037   3500   8952   219    .   .  
L01010110   4.0   1111  01010110  1   0     11.064935  1.112554  1.12987 11.625000 +0.560065   3383   7704   397    .   .  
L01010110   4.0   1111  01010110  1   1     11.250000  1.250000  1.10000 11.545455 +0.295455   3500   7948   398    .   .  
L01010110   4.0   1111  01010110  1   2     11.246964  1.230769  1.24696 11.581132 +0.334169   3500   7975   324    .   .  
L01010110   4.0   1111  01010110  1   3     11.000000  1.131579  1.30263 11.387097 +0.387097   3500   7798   456    .   .  
L01010110   4.0   1111  01010110  1   4     11.300000  1.337500  1.20000 11.589041 +0.289041   6580   8332   500    .   .  
L01010110   4.0   1111  01010110  1   5     11.260000  1.180000  1.22000 11.632653 +0.372653   3500   7945   496    .   .  
L01010110   4.0   1111  01010110  1   6     11.285714  1.285714  1.21429 11.500000 +0.214286   3500   8177   393    .   .  
L01010110   4.0   1111  01010110  1   7     11.285714  1.154762  1.20238 11.569492 +0.283777   6633   8928   500    .   .  
L01010110   4.0   1111  01010110  1   8     11.250000  1.000000  1.00000 11.486486 +0.236486   3500   8047   176  rows cliq
L01010110   4.0   1111  01010110  1   9     11.285714  1.269841  1.06349 11.400000 +0.114286   3500   7935   430    .   .  
L01010110   4.0   1111  01010110  1   10    11.400000  1.400000  1.10000 11.421053 +0.021053   3500   7884    12    .   .  
L01010110   4.0   1111  01010110  1   11    11.300000  1.300000  1.13333 11.600000 +0.300000   6226   8835   491    .   .  
L01010110   4.0   1111  01010110  1   12    11.270270  1.317568  1.21622 11.500000 +0.229730   3500   7869   122    .   .  
L01010110   4.0   1111  01010110  1   13    11.263158  1.368421  1.05263 11.481132 +0.217974   3500   7904   233    .   .  
L01010110   4.0   1111  01010110  1   14    11.250000  1.250000  1.25000 11.500000 +0.250000   3500   7725   117    .   .  
L01010110   4.0   1111  01010110  1   15    11.300000  1.293333  1.30000 11.666667 +0.366667   6294   8396   496    .   .  
L01010110   4.0   1111  01010110  1   16    11.200000  1.533333  1.10667 11.470588 +0.270588   3500   7818   214    .   .  
L01010110   4.0   1111  01010110  1   17    11.235294  1.117647  1.00000 11.444444 +0.209150   3500   7764    93    .  cliq
L01010110   4.0   1111  01010110  1   18    11.291667  1.500000  1.20833 11.500000 +0.208333   3500   7809   135    .   .  
L01011010   4.0   1111  01011010  1   0     11.142857  1.142857  1.21429 11.500000 +0.357143   3383   7871   327    .   .  
L01011010   4.0   1111  01011010  1   1     11.000000  1.000000  1.00000 11.500000 +0.500000   3500   7958   272  rows cliq
L01011010   4.0   1111  01011010  1   2     11.000000  2.000000  1.00000 11.000000 +0.000000   3500   7931     7    .  cliq
L01011010   4.0   1111  01011010  1   3     10.285714  1.142857  1.07143 10.500000 +0.214286   3500   7494    78    .   .  
L01011010   4.0   1111  01011010  1   4     11.142857  1.142857  1.08811 11.597403 +0.454545   6370   8573   499    .   .  
L01011010   4.0   1111  01011010  1   5     11.166667  1.359539  1.24686 11.500000 +0.333333   3500   8191   468    .   .  
L01011010   4.0   1111  01011010  1   6     11.259366  1.714697  1.25937 11.403101 +0.143735   3500   8356   500    .   .  
L01011010   4.0   1111  01011010  1   7     11.250000  1.400066  1.11813 11.500000 +0.250000   6571   8844   500    .   .  
L01011010   4.0   1111  01011010  1   8     11.250000  1.250000  1.00000 11.416667 +0.166667   3500   8050   262    .  cliq
L01011010   4.0   1111  01011010  1   9     11.187500  1.437500  1.18750 11.333333 +0.145833   3500   7905   250    .   .  
L01011010   4.0   1111  01011010  1   10    11.258572  1.324212  1.06600 11.500000 +0.241428   6021   8873   500    .   .  
L01011010   4.0   1111  01011010  1   11    11.283599  1.273884  1.06674 11.571429 +0.287829   3500   8451   333    .   .  
L01011010   4.0   1111  01011010  1   12    11.000000  1.692308  1.15385 11.500000 +0.500000   3500   7848   373    .   .  
L01011010   4.0   1111  01011010  1   13    11.000000  1.200000  1.20000 11.333333 +0.333333   3500   7737   293    .   .  
L01011010   4.0   1111  01011010  1   14    11.166667  1.407683  1.16667 11.557627 +0.390960   6221   8682   499    .   .  
L01011010   4.0   1111  01011010  1   15    11.000000  1.252941  1.37353 11.400000 +0.400000   3500   7889   388    .   .  
L01011010   4.0   1111  01011010  1   16    11.000000  1.426667  1.18667 11.329114 +0.329114   3500   7804   500    .   .  
L01011010   4.0   1111  01011010  1   17    11.000000  1.090909  1.27273 11.216667 +0.216667   3500   7831   141    .   .  
L01100110   4.0   1111  01100110  1   0     11.000000  1.000000  1.00000 11.000000 +0.000000   3383   8088   246  rows cliq
L01100110   4.0   1111  01100110  1   1     11.000000  1.234155  1.22711 11.000000 +0.000000   3383   7586   447    .   .  
L01100110   4.0   1111  01100110  1   2     11.000000  1.500000  1.50000 11.000000 -0.000000   3500   7537   316    .   .  
L01100110   4.0   1111  01100110  1   3     11.000000  1.000000  1.00000 11.000000 +0.000000   3500   7685   237  rows cliq
L01100110   4.0   1111  01100110  1   4     11.000000  1.000000  1.00000 11.000000 +0.000000   3500   7638     8  rows cliq
L01100110   4.0   1111  01100110  1   5     11.000000  2.000000  2.00000 11.000000 -0.000000   6459   7983   390    .   .  
L01100110   4.0   1111  01100110  1   6     10.500000  1.250000  1.25000 10.777778 +0.277778   3500   7744     6    .   .  
L01100110   4.0   1111  01100110  1   7     10.833333  1.000000  1.00000 10.894737 +0.061404   3500   7707   184  rows cliq
L01100110   4.0   1111  01100110  1   8     10.888889  1.285714  1.77778 10.920000 +0.031111   3500   7841   361    .   .  
L01100110   4.0   1111  01100110  1   9     11.000000  2.000000  1.00000 11.000000 -0.000000   6578   8038   282    .  cliq
L01100110   4.0   1111  01100110  1   10    10.000000  1.000000  1.00000 10.250000 +0.250000   3500   7650   115  rows cliq
L01100110   4.0   1111  01100110  1   11    10.500000  1.000000  1.00000 10.500000 +0.000000   3500   7635    41  rows cliq
L01100110   4.0   1111  01100110  1   12    10.000000  1.000000  1.00000 10.666667 +0.666667   3500   7653   300  rows cliq
L01100110   4.0   1111  01100110  1   13    11.000000  2.000000  1.00000 11.000000 +0.000000   6389   8244   351    .  cliq
L01100110   4.0   1111  01100110  1   14    11.000000  1.333333  1.16667 11.000000 +0.000000   3500   7707   294    .   .  
L01100110   4.0   1111  01100110  1   15    11.000000  1.000000  1.00000 11.000000 -0.000000   3500   7732    22  rows cliq
L01100110   4.0   1111  01100110  1   16    11.000000  1.000000  1.00000 11.000000 +0.000000   3500   7733   153  rows cliq
P40         4.0   ....  ........  1   0     12.118529  1.016779  1.10425 12.194538 +0.076009   3383  24195  1736    .   .  
P40         4.0   ....  ........  1   1     12.193344  1.321445  1.07602 12.259214 +0.065870   6616  32389  2036    .   .  
T0          4.0   1111  ........  0   0     11.896104  1.009984  0.00000         -        -   1471  13518     0    .  cliq
T0          4.0   1111  01010101  0   0     11.385844  1.000000  0.00000         -        -   1471  14461     0  rows cliq
U40         4.0   ....  ........  1   0     12.141385  1.021979  1.07433 12.196598 +0.055213   3383  20805   500    .   .  
U40         4.0   ....  ........  1   1     12.121542  1.036364  1.05170 12.202884 +0.081342   3383  19644   499    .   .  
U40         4.0   ....  ........  1   2     12.197049  1.422294  1.12369 12.242251 +0.045202   3500  19808   492    .   .  
U40         4.0   ....  ........  1   3     12.197919  1.168081  1.17712 12.225492 +0.027572   5798  22753   498    .   .
```

---

## 4. Reading

### 4.1 What the branch is worth, and what the cliques are worth

At `t = 4` the corner branch and the slot branch do almost all of the work and the cliques do
little:

| step | value at `t = 4` |
|---|---|
| pure fractional packing, certified exactly over ALL poses (`DUAL.md`, `dual_PC1_support.txt`) | **12.163061** |
| pure, restricted poses + anchor cliques (`U40`, `P40`) | `12.12 – 12.19`; matched pure `12.19 – 12.26` |
| + corner branch `k = 4` (`B40`, `K40`), pure | **12.000000** (three independent runs) |
| + corner branch `k = 4` + cliques | `11.95 – 11.98` |
| + level-2 slot branch, the four `m = 4` leaves, + cliques | **`10.3 – 11.5`** |

The matched clique gain (same poses, same rows, clique rows off) is `+0.02 – 0.08` in the pure and
corner-leaf rows and `+0.06 – 0.56` in the slot leaves — the same order as `RECONCILE.md` §5's
`+0.05 – 0.095` on the cover side at `3.99`, and consistent with `CLIQUE_CONTINUUM.md`'s `0.10`
having been an over-estimate.  **Cliques alone do not take `t = 4` below 12**: the pure
clique-strengthened value sits at `12.12 – 12.19` on 3.4k–6.6k poses, which is where
`CLIQUE_CEILING.md` left it (`11.8 – 12.0`, drifting up, with box cliques).  What takes it below 12
is the *slot* branch, worth `0.5 – 1.7` on top of the corner branch — the same half-unit-to-unit
that `level2-design.md` §4.4 measured, now measured with cliques as well and on a 2.3× larger pose
set.

### 4.2 The sharp form of the question in an `m = 4` leaf

In an `m = 4` leaf the four corner boxes carry 1 each and the four occupied slots carry 1 each, so

> **leaf value = 8 + (mass centred in the interior `[1,3]^2`)**,

which the region split of every leaf optimum confirms to twelve digits (e.g. `L01010101` at step 2:
corners `[1,1,1,1]`, slots summing to `4.000000`, interior `3.403693`, total `11.403693`).  So a
leaf fails to close **iff the interior can carry fractional mass `>= 4`** beside the pinned frame.
Measured interior mass in the four leaves: `2.3 – 3.5`.

**The interior sub-problem is where the instrument can be checked against ground truth, and it
fails the check.**  Run with the four corner boxes and all eight slots pinned to 0, so that every
pose carrying mass is interior-centred:

| run | poses | LP | `M` | `kmax` | |
|---|---|---|---|---|---|
| `INT40` pattern `00000000` (one unrefined evaluation, on a row set built for a different pattern) | 3500 | **4.000000** | `1.000000` | `1.000000` | **converged** |
| `INTP` stage 0 (seed pose set, priced against its own dual) | 1471 | **5.000000** | `1.000000` | `1.000000` | **converged** |
| **witness** (`level2_capacity.py --t 4.0 --box 1 3 1 3 --k 6`) | — | `>= 6` | — | — | six unit squares, centres in `[1,3]^2`, pairwise margin `+0.0567`, all inside `[0,4]^2`; mass 1 on each is a feasible measure |
| `INTP` stage 17 | 4000 | `7.433677` | `1.000000` | `1.02` | rows converged |
| `INTP` stage 20 | 4000 | `7.508327` | `1.000000` | `1.000000` | **converged** |
| `INTP` stage 21 | 4000 | **`7.528038`** | `1.000000` | `1.000000` | **converged** |
| `INTP` stage 22 | 4000 | `7.583111` | `1.05` | `1.05` | still climbing `+0.05`/stage |

The whole trajectory, one pricing stage per step, `4000` poses throughout after the first:

```
5.000  5.731  6.000  6.189  6.394  6.580  6.654  6.778  6.826  6.911  6.956
7.087  7.099  7.229  7.314  7.324  7.386  7.434  7.435  7.473  7.508  7.528  7.583
```

Monotone, no sign of stopping.  **A converged LP value on the seed pose set was `4.000000`; the same
LP on the same problem is at `7.53` after 21 rounds of column generation and still rising.**  This is
the honest calibration of the whole screen: `M <= 1` and `kmax <= 1` certify that the *measure* is
feasible, never that the *pose set* was rich enough, and here the pose set was short by `3.5` on a
quantity whose threshold is `4`.  (`level2-design.md` §5 records the same six-square configuration at `t = 3.98`,
margin `+0.0884`, "six squares at 34–44 degrees leaning out of the interior into the wall gaps" —
poses a `0.25 / 15 deg` seed grid does not contain and one round of `0.04 / 2.5 deg` pricing only
half finds.)

So the arithmetic of an `m = 4` leaf at `t = 4` is explicit and the answer turns on a quantity this
screen measures badly:

> `leaf value = 8 + interior mass`; the leaf fails to close iff the interior carries `>= 4`.  With
> the frame pinned the LP gives `2.3 – 3.5`; the same LP under-reported the *frame-free* interior by
> `3.5` on the pose set the leaf runs use.

The frame does cost the interior something — the corner and slot squares stick into `[0.3, 3.7]^2`
and eat the coverage budget the interior squares need — and the leaves do come in below 12.  But
`0.4 – 1.0` of deficit against an instrument with a demonstrated `3.5` of pose-set slack on the very
quantity in question is not a verdict.

One thing the leaf values *cannot* be beaten by, at least: an **integral** witness.  Twelve unit
squares in a leaf's region pattern with all pairwise margins bounded away from zero would give a
feasible measure of mass exactly 12 — but such a configuration can be scaled to twelve squares of
side `> 1` in `[0,4]^2`, i.e. it would refute `s(12) = 4` outright.  So a leaf value `>= 12` has to
come from the *fractional* excess, which is what the pure method has at `t = 4` (`12.163`) and what
the branch and the cliques are trying to remove.  (The same argument with `s(13) = 4`, which is a
theorem, rules out an integral witness of mass 13.)

### 4.3 Why this settles nothing

* **The instrument is short by `1 – 2` where it can be checked** (§4.2).  That is the governing
  caveat; everything below is secondary to it.
* **Almost nothing converged.**  Of the stage values in §3, five have `M <= 1` and `kmax <= 1`
  simultaneously and are therefore genuine lower bounds on their leaf *over their own pose set*;
  the rest have `M` in `1.02 – 1.7` and are upper bounds on the restricted-pose value and bounds on
  nothing else.  The best *converged* leaf values are `11.000000` (`01011010`, `01100110`) and the
  best unconverged ones `11.48` (`01010101`).
* **The values climb under refinement — but part of the climb is row ageing, not columns.**
  `01010101` over eight stages: `11.376, 11.404, 11.410, 11.477, 11.472, 11.501, 11.521, 11.562`,
  its matched pure reaching `11.82`; `k = 4` with cliques `11.945 -> 11.982`.  Read `M` alongside:
  it goes `1.118 -> 1.53` over the same eight stages.  `--row-age 3` drops point rows that have been
  slack and dual-free, which is a *relaxation* and pushes the LP up, so a stage-to-stage rise with
  `M` rising alongside it is not a column-generation gain.  The trustworthy numbers are the
  converged ones (`M <= 1`), and the best of those over all four leaves is `11.000000`.
  `CLIQUE_CEILING.md`'s warning that these values "drift up slowly with pose refinement" applies
  verbatim, with this extra confound on top.
* **The pricing gap does not bound the remainder, and here it is worse than uninformative.**  In a
  region-constrained LP the equality's multiplier is free; a candidate pose in a region the branch
  forces mass into has reduced cost `1 - capture + |lam_R|`, which reached `+172` in these runs
  while the objective moved by `0`.  So the gap says nothing.  (It also biased the pricer badly
  until fixed — §6.)
* **The tree screened here is one corner leaf.**  At `t = 4` no corner leaf closes on its own, so a
  real proof needs the slot tree under `k = 0,1,2,3` as well, where the chord lemma gives
  `2 k_strip + slots <= 3` per wall and the leaf count is in the hundreds (`level2-design.md` §3).
  Nothing here says anything about those.
* **All four `m = 4` leaves are non-empty at `t = 4`, and by the worst possible object.**  Take the
  `4 x 4` grid of unit squares in `[0,4]^2` and delete one slot square per wall: the remaining
  twelve have pairwise disjoint interiors, sit `4` in the corner boxes, `4` in the slots (one per
  wall) and `4` in the interior, and each wall strip carries exactly `3` centres, saturating the
  chord bound.  The sixteen ways of choosing which slot per wall give, up to `D4`, exactly the four
  `m = 4` patterns.  So no `m = 4` leaf can be refuted by capacity; the LP value has to do all of
  it, and the object it has to beat is the `4 x 4` grid itself — the configuration the whole
  `s(12) = 4` question is about.

---

## 5. The obstacle: the leaf optimum sits on the slot boundary

`level2-design.md` §9 reassures that at `t = 4` "the sixteen grid poses sit `0.5` from every
boundary, so the region membership tests are not zero-margin at the tight poses".  That is a
statement about the *cover's* tight poses.  The **packing side's optimum goes to the boundary**.

Mass within `1e-7` of a region boundary in the centre plane (`c_x` or `c_y` in `{1, 2, 3}`), on the
leaf optima:

| leaf | leaf mass | mass on a boundary | share | heaviest single boundary pose |
|---|---|---|---|---|
| `01010101` | `11.410038` | `0.733798` | `6.4 %` | `0.182` at `(2.0, 0.5, 0 deg)` |
| `01010110` | `11.000000` | `0.460526` | `4.2 %` | `0.197` at `(2.0, 0.5086, 1 deg)` |
| `01011010` | `11.142857` | `1.233684` | `11.1 %` | `0.230` at `(0.6579, 2.0, -23.5 deg)` |
| `01100110` | `10.888889` | `0.095238` | `0.9 %` | `0.032` at `(2.0, 0.6549, -22.5 deg)` |

and, at the two stages whose value is a genuine bound (`M = 1`, `kmax = 1`):

| leaf | converged value | boundary mass | where |
|---|---|---|---|
| `01011010` step 1 | `11.000000` | **`1.000000` — a whole unit at one pose** | `(3.499547, 2.000000)`, the `W_2 / W_3` boundary |
| `01100110` step 0 | `11.000000` | **`1.000000` — a whole unit at one pose** | `(2.000000, 0.500453)`, the `W_0 / W_1` boundary |

There `c_y = 2.0 = t/2` (resp. `c_x`) to within `1e-14`: the pose sits on the mid-wall line that
separates an *occupied* slot from an *empty* one, and it satisfies `mu(W_3) = 1` and `mu(W_2) = 0`
only because `level2_regions.classify` breaks the tie one way rather than the other.  This bites in
practice: printing the pose checkpoint to 13 decimals instead of 17 flipped the pose's region on
read-back, which is how it was noticed.  `save_poses` now writes 17 significant digits.

Three consequences.

1. **A cell-based verifier cannot implement this branch as it stands.**  `verify/`'s region test is
   "is this cell's bounding box inside the box?", exact for convex cells; a cell straddling the
   `W_2 / W_3` boundary is in neither box.  The branch is exhaustive only if every pose is assigned
   to exactly one slot, and at the boundary that assignment is being decided at `1e-14`.
2. **This is `level2-design.md` §2.5's mid-wall escape**, in its extreme form.  That note observed
   that about `4 %` of the `t = 3.98` wall mass already sat at `c_y ~ t/2`, and proposed Design B
   (slots shortened to `sqrt3/2`, leaving a `0.268` gap at the middle of each wall) as the fallback.
   At `t = 4` the packing-side optimum does not merely graze the mid-wall — it parks a *unit* of
   mass there.  Design B avoids the tie-break, at the price of `43` leaves instead of `15` and of a
   gap the LP can hide mass in, which raises the leaf values.
3. **It makes the go verdict softer than the numbers look.**  The two `11.000000` values — the only
   fully converged leaf bounds here — are attained by a measure a certificate could not carry.  What
   a *robust* branch (one whose region test a verifier can decide) is worth at `t = 4` is unmeasured,
   and is the first thing a follow-up should measure.

---

## 6. Reproduce

```sh
M=/home/evand/math/square-packing/s12/runs     # the main tree's runs/, read-only

# regression against notes/level2-design.md 4.4 (control 11.896104, 01010101 11.385844)
python3 search/t4screen.py 4.0 T0 --corners 1111 --patterns 01010101 --stages 1 --warm-stages 0 \
    --rowloops 14 --warm $M/dual_PC1_support.txt --threads 2 --method highs

# calibration against the cover side at t = 3.98 (must come in <= ~11.96)
python3 search/t4screen.py 3.98 C98 --corners 1111 --cliques --cq-interior --patterns "........" \
    --stages 8 --warm-stages 8 --rowloops 12 --warm $M/dual_PD1_support.txt \
    --warm $M/branch_J16i_dual.txt --cq-load $M/branch_J16i_cliques.txt --threads 2 --time 780

# the t = 4 screen: one process per table row, 1 thread each, 5 chunks of 450 s
N=5 S=450 sh search/t4screen_phase2.sh

# the interior alone (corner boxes and all eight slots pinned to 0)
sh search/t4screen_loop.sh 4.0 INT40 3 450 --corners 0000 --patterns 00000000 \
    --warm $M/dual_PC1_support.txt --warm $M/cqx_PURE99_support.txt --cliques --cq-wall \
    --cq-interior --cq-max 500 --cq-age 3 --row-age 3 --pose-max 3500 --threads 1 \
    --rowloops 10 --stages 4 --row-cap 6000

# the table
python3 search/t4screen_table.py runs/t4_*.log
```

Determinism: a single `t4screen.py` invocation is a pure function of its arguments except for the
`--time` cut (numpy, HiGHS and the separators are deterministic and every sort is stable), so a run
with `--time` large and `--stages` / `--rowloops` as the only stopping rule reproduces bit for bit.
The chunked driver is not, because the cut lands mid-loop.

**Two defects found by the built-in diagnostics**, both fixed and both worth recording because they
are easy to reintroduce:

* `rc-check` (the maximum `|reduced cost|` over columns carrying mass, which must be `0` at an
  optimum) read `+1.0` to `+4.5` on the first row-aged runs.  Cause: row ageing reindexes the point
  rows, so the dual vector returned by the inner loop no longer lined up with them and the pricer
  worked on garbage.  Both sifts now run at the *top* of the inner loop, leaving only appends after
  the last solve.  `rc-check` reads `1e-14` since.
* The pricer ranked candidates by plain reduced cost, which in a region-constrained LP is dominated
  by `|lam_R|` — so it offered only frame poses and never priced the interior, the one region that
  can raise an `m = 4` leaf's objective.  Pricing and column sifting are now stratified by region.
  The `01010101` trajectory (`11.376 -> 11.477`) is with the fix in the later stages.

**Compute actually spent**: about two hours on `<= 8` threads (7 single-threaded processes at the
peak), plus three earlier two-thread runs.  Two of those three overran the 15-minute per-process cap
(`B40` 16.7 min, `C98` 17.5 min) because the time check only fires between inner iterations and a
single 7835-pose LP solve took over four minutes; the chunk length was then cut to 450 s with a hard
`timeout` of 690 s, and nothing since has come near the cap.

---

## 7. What a follow-up would need

1. **Decide the region test before anything else** (§5).  Either make `verify/`'s cell-vs-box test
   agree with a declared tie-break at the slot boundary, or move to Design B and re-measure the leaf
   values with the mid-wall gap present.  Until then the `m = 4` leaf values here are values of an
   LP no certificate can dualise.
2. **Feed the leaves the interior poses they are missing, and see how far they move.**  §4.2 shows
   the interior sub-problem is under-measured by `1 – 2` on these pose sets, and the four leaf
   values sit `0.4 – 1.0` below 12.  The cheapest decisive-ish experiment is to warm-start each
   `m = 4` leaf from the *refined interior* pose set (`--warm runs/t4_INTP_poses.txt`, which
   contains the `34–60 deg` configurations the seed grid misses) and re-measure.  If the leaf values
   do not move, the frame really is what limits them; if they jump, this screen said nothing.
3. **Converge one leaf.**  `01010101` is the hardest (`11.56` and climbing, `M = 1.5`) and never
   converged its rows.  One leaf run to `M <= 1` and `kmax <= 1` on a 20k-pose set with the
   stratified pricer and row ageing OFF would turn its number into an actual lower bound.  That is a
   cover-side-sized computation (`level2-design.md` §8: 13–27 h per leaf), not a screen.
4. **A ceiling attempt, the only decisive direction available on this side.**  Exhibit, in one leaf,
   a measure exactly certified feasible (coverage `<= 1` at every exact arrangement vertex,
   `mu(K) <= 1` for every anchor clique, region masses exactly the leaf's counts) with mass `>= 12`.
   `clique_ceiling.py --exact` already does the coverage and clique half in exact rationals; what is
   missing is the region equality (its `--kmass` top-up does not survive the scaling step) and an
   exact anchor-clique maximiser.  Forty core-hours of that at `t = 3.99–4.0` found nothing above 12
   (`CLIQUE_CEILING.md`); this screen adds the corner and slot branches to the same negative.
5. **The rest of the tree.**  The slot tree under `k = 0,1,2,3` (`level2-design.md` §3) is untouched
   and is where the leaf count explodes.
6. **The chord lemma** is what restricts attention to the four `m = 4` patterns; it is still
   unproved in this repo (`level2-design.md` §10, item 1).

---

## 8. State of the runs at write-up time

Every run in §3 is a `t4screen_loop.sh` chunk sequence.  All of them were left to finish their
declared chunk budget; none is converged and none was going to be.  The last state each had
reached, and the number that matters for it:

| run | what | reached |
|---|---|---|
| `L01010101` | leaf `01010101`, cliques, priced against its own dual | 10 stages, LP `11.58`, `M = 1.28`, matched pure `11.73` |
| `L01010110` | leaf `01010110` | 7 stages, LP `11.29`, `M = 1.29` |
| `L01011010` | leaf `01011010` | 7 stages, LP `11.26`, `M = 1.71`; best converged `11.000000` |
| `L01100110` | leaf `01100110` | 13 stages, best converged `11.000000` (attained six times) |
| `K40` | corner leaf `k = 4`, cliques | 6 stages, LP `11.963 – 11.978` (flat), matched pure `12.000000` at **every** stage |
| `U40`, `P40` | pure, cliques | 4 / 2 stages, LP `12.12 – 12.20` |
| `C98` | the `t = 3.98` calibration | 1 stage, LP `11.79`, matched pure `11.90`, both `<= 11.96` ✓ |
| `INT40`, `INTP` | the interior alone | `4.000000` (converged, unrefined) → **`7.528038`** (converged, 21 pricing stages) → `7.583`, still rising |
| `S01010101` | follow-up 2: the hardest leaf warm-started from `INTP`'s refined interior poses (38 745 poses, row ageing off) | 3 iterates, `12.000000, 12.000000, 11.949104` with `M` falling `2.65 -> 2.00 -> 1.42` — the row generation had not caught up, so `11.949` is an upper bound on the restricted value and a bound on nothing.  **Left running.**  Note the chunking failure mode it exposes: at `240 s` per LP a `450 s` chunk never completes a stage, so it never checkpoints and the next chunk starts over.  Give this run one long process, not chunks. |

Nothing here is a certificate, an exact computation, or a bound on `V(4)` in the decisive
direction.  What is exact: the `12.163061` pure reference (`DUAL.md`), the `s(13) = 4` and
`s(12) = 4` arguments in §4.2 that rule out integral witnesses, and the twelve-square grid
configuration of §4.3.  Everything else is a float LP on a restricted pose set.
