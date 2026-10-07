# census + chains: the integral margin of every pattern leaf, and what holds the plateau shut (2026-09-19)

Code: `search/census.py`, `search/chains.py` (both new; `rank8.py`, `bentz.py` read, not touched).
Nothing in `verify*/`, `certificates/`, `lean/` is touched: these are measurements.  Runs are in the
gitignored `runs/` (`census_families.json`, `census_pairs.json`, `census_run.jsonl`,
`census_refine.jsonl`, `chains_leafA_samples.txt`, `chains_leafA.jsonl`, `chains_census.jsonl`,
launchers `launch_census*.sh`, `launch_chains_*.sh`), so every number is quoted here.

Semantics as `RANK8.md`: `t = 4`, closed squares, margin = sup over the closed relaxation of a leaf
of the min pairwise separating-axis gap.  A **leaf** is a family of twelve pairwise-disjoint
realisable patterns on Bentz's `P0`, one labelled square per pattern (`bentz.py tree`: 86,403
leaves, 10,945 up to `D4`; reproduced).  Leaf A of `RANK8.md` is class `#10944`.

**Everything negative below is a multistart measurement and a lower bound on the margin, not a
certificate.  The instrument's false-negative rate is high near zero: see §1, last paragraph.**

## 0. Verdict

1. **No class has positive margin** (no counterexample), over 10,945 classes.
2. **Leaf A is one of `>= 2,779` zero-margin leaves — a quarter of the tree.**  `1,425` classes
   contain an exact `4x4`-tiling-minus-4 configuration in their closed pattern regions (pure
   enumeration: 274,269 (tiling, assignment) pairs; all 8,546 stored seeds re-checked: patterns
   right, min gap exactly `0`).  A further `1,297 + 57` classes are at `0` (to `1e-14`) with **no**
   tiling in them, on configurations tilted up to `44 deg`, with patterns no tile can hold
   (`{d0,d2,d3}`, `{a3,b3,c3}`, ...).  So a per-leaf hand proof is not a strategy; the proof needs a
   uniform mechanism for zero-margin leaves.
3. **What holds a plateau point shut is, 99 % of the time, one straight wall-to-wall row of four
   axis-parallel squares** — also in configurations where other squares are tilted `44 deg`.
   The only other core seen is a diagonal chain through a corner contact, and only at exactly
   axis-parallel configurations.  No core ever needed a tilted square.  First-order feasible plateau
   points (which would be genuine closed 12-packings): none.
4. The pattern tree is blind to this.  The separation structure (which pairs are x-separated,
   y-separated, or separated only along a tilted edge normal) is what carries the obstruction.

## 1. Census (`census.py families | pairs | run | refine | report`)

Pass 1: per class, 40 random starts from per-pattern pose pools (6e6 samples, half of them
near axis-parallel), the class's tiling seeds + 12 nudges where it has any, 2 rounds of 16 kicks
around the 4 best; 14 processes, 97 min.

| measured margin | classes |
|---|---|
| `> +1e-7` | **0** |
| `|.| <= 1e-7` | 2,722  (1,425 tiling-compatible, 1,297 not) |
| `(-1e-4, -1e-7)` | 29 |
| `[-0.02, -1e-4)` | 36 |
| `[-0.2, -0.02)` | 1,601 |
| `[-0.5, -0.2)` | 6,557 |

By `K = sum(|pi| - 1) = 4 - u` (zero / negative): `K=0` 246/0, `K=1` 644/824, `K=2` 878/2728,
`K=3` 705/3227, `K=4` 249/1444.  Every `K = 0` class (all singletons, four points uncovered) is
tiling-compatible.  All 1,425 tiling-compatible classes came back at exactly `0` (none left negative
— a weak check, since they were seeded with their own tiling, which sits at `0` — and none positive).

Pattern pairs (`pairs`, two squares, six unknowns, 60 starts + 400 kicks): of 3,052 disjoint
pattern pairs exactly one `D4` orbit is pairwise infeasible, `{a0} | {b0,c1,d0}` and images, at
`-0.040244`.  So the negative leaves die through three or more squares, not pairwise.

**Reliability.**  Pass 2 (`refine`: 200 starts, 6 rounds of 100 kicks, seeded from the incumbent)
on the 65 classes pass 1 left in `[-0.02, -1e-7)` moved **57 of them to exactly 0**; the 8
survivors are at `-6.3e-4` (`#3631`), `-8.9e-3`, `-9.6e-3`, `-1.2e-2`, `-1.4e-2`, `-1.73e-2` (x2),
`-1.75e-2`.  The pair table also first reported `-0.398` for seven images of a pair whose value is
`-0.040`.  So pass-1 negatives in the middle band (`-0.2 .. -0.02`, 1,601 classes) are not to be
trusted as a count of compact leaves, and `-0.35` recurs suspiciously often in the big block.  The
statement that survives: **at least 2,779 zero-margin classes, no positive ones.**

## 2. Chains (`chains.py analyse`)

First-order analysis at an arbitrary configuration, **no pattern constraints**: unknowns
`(dx, dy, dphi)` per square and `delta`; every active wall row and every touching pair must improve
at rate `>= delta`; a pair's gap is a max over the four edge normals of a min over sign branches, so
active normals are a disjunction (binary) of conjunctions; MILP (HiGHS).  `delta* > 0` would give a
genuine closed packing `z + eps v`.  `delta* = 0`: a deletion filter over contacts (pair contacts
and wall contacts as units; contacts of tilted squares tried first; best of 3 orders) returns an
irreducible rigid **core**.  Checked against finite differences of `rank8.pair_gap`, and on a row of
four (`x4LR`), a row of three (`delta = 1/3`), a row of four with a `1e-3` hole (`delta = 1/4`).

| sample | configs | one straight row of four (`x4LR` / `y4BT`) | corner-contact compound core | other |
|---|---|---|---|---|
| leaf A plateau (`rank8.py margin`, 3,016 starts, tilts to `35.8 deg`) | 2,145 | 2,137 | 6 | 1 `x6LR`, 1 near-degenerate (`0.02 deg`) |
| best config of each zero census class (tilts to `44.0 deg`) | 2,722 | 2,677 | 44 | 1 not at a stationary point (`-3.9e-8`, line search gains nothing) |

Every compound core sits in a configuration with **all twelve squares axis-parallel** (tilt `0.00`)
and has the same shape: two squares touching at a corner, one pinned towards the left/bottom walls by
x- and y-runs, the other towards the right/top, run lengths `(k, 4-k)` in each direction — e.g. the
leaf-A pinwheel `x:C0-D0[L]; x:D3-C3[R]; y:C1-D0[B]; y:D3-C2[T]; corner:D0/D3`.  Up to three
corner contacts occur in one core (`#8849`).

Caveat: the samples are points where a max-min optimiser stopped, not a uniform sample of the
interior-disjoint set.  That they are all first-order rigid is what `s(12) = 4` demands; that the
core is almost always a row of four is the information.

## 3. What this suggests (untested)

For axis-parallel squares the whole theorem is separation structure: any two disjoint squares are
strictly separated in `x` or in `y`; a chain of four in either order needs width `> 4`; two partial
orders with every pair comparable and chains `<= 3` have at most `3 x 3 = 9` elements — the
"axis-parallel `<= 9`" theorem of `BENTZ.md` C1 without a point set.  With tilts, a pair may be
separated only along a tilted edge normal, and that third kind of pair is the entire content of
`n = 12` (and the reason `n = 5`, `n = 11` beat the tiling).  A proof that branches on **separation
type per pair** rather than point incidence would kill every row-of-four leaf by an exact linear
argument (four widths `>= 1`, strict separation, a wall of length 4: no margin needed, a rational
Farkas certificate), and would confine the quantitative work to leaves with many tilted-only pairs.
Open: how many such leaves, whether the corner-contact cores have an equally clean exact kill, and
what replaces "width `>= 1`" when the separating normals in a run are tilted by different amounts.
