# Tightening the certificate: re-optimised weights, sparser point sets, larger containers

Code: `search/tighten.py` (new), plus a `--warm CERT` option in `search/lp_search.py`.
Everything below ran on one 32-core machine in about two hours of wall clock (several runs at a
time); logs and intermediate certificates are in the gitignored `runs/`.  **Outcome: the bound
moves from 3920/997 = 3.931795 to 980/247 = 3.967611, and the certificate at the old bound
shrinks from 788 to 224 points.**

## Starting point

`search/CEILING.md` observed that the cover LP over the shipped certificate's **own** 788 points
reaches ≈ 11.82 against the shipped 11.93, so ~0.1 of slack was available from the weights alone.
This note follows that up.  Three things were done, in order:

* **A. re-optimise the weights** at fixed point locations, with the exact verifier as the
  separation oracle (so the loop terminates with a file the verifier accepts, nothing else);
* **B. sparsify**: drop points (a reweighted-L1 heuristic under a total-weight budget);
* **C. grow the point set** at a larger container by column generation (exact reduced-cost
  pricing), which is what actually moved the bound.

## The loop (`tighten.py`)

The points are kept as D4 orbits (so the verifier's `[0,45°]` reduction still applies); only
the orbit weights `x_k` are variables.

```
repeat
    x  <- LP:  min  Σ_k |orbit_k| x_k    s.t.   A x >= 1 + m_lp,   x >= 0
               (A_rk = number of atoms of orbit k in the CLOSED square of row r)
    export x/(1+m_probe), rounded DOWN at weight denominator 1e12
    run   verify <probe> 12 N 32 <topk> <witness-file>          (N = 6000, topk = 6)
    every witness -- a placement of the verifier's shrunk sigma_k-square at bin angle
    theta_k whose captured weight is < 1 -- becomes an LP row with exactly that square
    (half-side sigma_k/2 from the verifier's own integer formula)
until the probe has no violated placement
```

`m_lp = 2e-6`, `m_probe = 1e-6`.  "No violated placement" for the probe means the true `x`
covers `>= 1 + 1e-6` (up to `1e-9`) on the whole `N = 6000` net, so the final export (rounded
**up** at `W = 1e7`) can only cover more.  The final file is then run at `N = 6000` and
`N = 12000`; the second net is not implied by the first (its squares are shrunk less but sit at
other angles), so its witnesses are fed back if needed — it never was.  The warm start is the
lattice separation of `nu_f.py` (closed unit squares on the `(0.01, 0.01)` lattice with
captured weight `< 1.05` under the input weights, ~17–20k rows).  A round trip (LP + verifier
at `N = 6000` + row insertion) takes 3–5 s with ~100 orbits and up to ~90 s with ~600 orbits and
100k rows; the fixed-point runs converge in 10–30 rounds.

The container is `K/D` with `K = s·D` of the input (7840 for the shipped file) and the integer
coordinates never change, exactly as in `scale_to_critical.py`; `--Dp` picks another `D`
(`--mul` multiplies the coordinates so that `D` can be chosen to `1/mul`).

**Sparsify** (`sparsify` mode).  Reweighted L1: `min Σ_k c_k |orbit_k| x_k` with
`c_k = 1/(x_k + 1e-4)` under the same rows plus `Σ |orbit_k| x_k <= budget`; orbits that hit
zero are fixed at zero for the following rounds (so the row set stays consistent); after the
support stops shrinking, the plain min-weight LP is re-run on the support and finalised as
above.

**Column generation** (`--colgen k`).  The LP dual `y` is a packing measure on the rows.  The
reduced cost of the orbit of a point `p` is `|orbit|·(1 − covsym(p))` with
`covsym(p) = (1/8) Σ_g cov(g·p)`, `cov(q) = Σ_r y_r [q ∈ Q_r]`.  `price()` evaluates `covsym`
on a `0.01` grid of the D4 fundamental domain, refines the best candidates on the integer grid
`K/D` within one pitch, and adds the points with `covsym > 1` (up to `--cg-want` per round).
A first version used `nu_f.max_coverage`'s dilated leaves instead; it stalled (its leaves snap
to points already present) and is not what is in the file now.

## Results

Verifier at `N = 6000` (`min` = exact minimum captured weight; identical at `N = 12000` for
every row) and `xcheck.py --all` at `N = 6000` for every certificate marked *shipped*.

| certificate | s | points | total weight | min @ N=6000 | how | wall |
|---|---|---|---|---|---|---|
| `certificates/s12_lower_3.931795.txt` (shipped before) | 3920/997 = 3.931795 | 788 | 11.9329864 | 1.0000023 | — | — |
| `runs/tight_A1.txt` | 3920/997 = 3.931795 | 504 | 11.8756852 | 1.0000032 | A: re-optimised weights, same points | 97 s |
| **`certificates/s12_lower_3.931795_sparse.txt`** *(shipped)* | 3920/997 = 3.931795 | **224** | 11.9834372 | 1.0000022 | B: sparsified at the shipped bound | 217 s |
| `runs/tight_S1_D1990.txt` | 784/199 = 3.939698 | 496 | 11.8991548 | 1.0000036 | A at D=1990 | 110 s |
| `runs/tight_S1_D1988.txt` | 280/71 = 3.943662 | 592 | 11.9165372 | 1.0000037 | A at D=1988 | 43 s |
| `runs/tight_S1_D1986.txt` | 3920/993 = 3.947633 | 584 | 11.9683796 | 1.0000034 | A at D=1986 | 50 s |
| `runs/tight_S2_D19855.txt` | 15680/3971 = 3.948628 | 624 | 11.9633556 | 1.0000031 | A at D=19855/10 — best with the shipped points | 110 s |
| `runs/sparse_B2.txt` | 15680/3971 = 3.948628 | 368 | 11.9934340 | 1.0000030 | B at that bound | 208 s |
| `runs/tight_C1.txt` | 245/62 = 3.951613 | 1272 | 11.9353104 | 1.0000050 | C: column generation from the shipped points | 1632 s |
| `runs/sparse_B3.txt` | 245/62 = 3.951613 | 392 | 11.9916200 | 1.0000031 | B applied to `tight_C1` | 727 s |
| `runs/tight_C4b.txt` | 3920/991 = 3.955600 | 1108 | 11.8934460 | 1.0000034 | C: colgen from `tight_C1` at D=1982 (15 rounds, pricing on a 0.01 grid), then cut-only on the support | 856 + 476 s |
| `runs/tight_C5b.txt` | 392/99 = 3.959596 | 968 | 11.9389496 | 1.0000045 | C: same from `tight_C1` at D=1980 | 614 + 431 s |
| `runs/sparse_B4.txt` | 392/99 = 3.959596 | 428 | 11.9918800 | 1.0000030 | B applied to `tight_C5b` | 446 s |
| `runs/tight_C6b.txt` | 3920/989 = 3.963600 | 1192 | 11.9504824 | 1.0000046 | C: colgen from `tight_C5b` at D=1978, then cut-only on the support | 956 + 572 s |
| `runs/sparse_B5.txt` | 3920/989 = 3.963600 | 548 | 11.9925972 | 1.0000031 | B applied to `tight_C6b` | 858 s |
| `runs/tight_C7b.txt` | 980/247 = 3.967611 | 1344 | 11.9844236 | 1.0000041 | C: colgen from `tight_C5b` at D=1976, then cut-only on the support | 878 + 807 s |
| **`certificates/s12_lower_3.9676.txt`** *(shipped)* | 980/247 = 3.967611 | **764** | 11.9962288 | 1.0000033 | B applied to `tight_C7b` (budget 11.998) | 879 s |

Beyond 3.9676: the same procedure at `D = 1975` (`s = 1568/395 = 3.9696`) reached LP 11.94
with column generation still running (849 orbits, 14.8k violated placements per round, ~3 min
per LP), and the cut-only re-optimisation on its support snapshot ended at **12.05** — i.e. no
certificate at 3.9696 within the time box.  That is consistent with `CEILING.md`'s estimate that
`ν_f` crosses 12 around 3.965–3.97.

With the shipped **points** and only the weights free, the bound goes from 3.931795 to
3.948628 (`D = 19855/10`); at the next step (`D = 1985`, `s = 3.94962`) the LP total jumps to
13.8 — a cliff, not a slow crossing: some placement runs out of usable points, and the shipped
set has nothing to offer there.  Column generation removes the cliff: from the shipped points
at `D = 1984` (`s = 245/62 = 3.951613`) the LP drops from 14.03 (fixed points) to 11.935 with
~800 orbits added (of which 162 end up in the support).

Sparsification costs almost nothing in the bound: 224 points suffice at the shipped bound
(788 before), 392 at 3.9516, 428 at 3.9596, 548 at 3.9636 and 764 at 3.9676 (where the LP
optimum 11.984 leaves only 0.014 of budget).  The reweighted-L1 heuristic stops after 3–4
rounds; `--budget` trades points for weight margin.

**Shipped** (all three checks pass: `verify` at N = 6000 and 12000, `xcheck.py --all` at
N = 6000, `export_points.py --roundtrip`):

* `certificates/s12_lower_3.9676.txt` — **s(12) ≥ 980/247 = 3.967611**, 764 points, total
  weight 7497643/625000 = 11.9962288 (best bound);
* `certificates/s12_lower_3.931795_sparse.txt` — the shipped bound 3920/997 with **224**
  points (total 29958593/2500000 = 11.9834372) instead of 788 (smallest certificate).

The certificates at 3.9516, 3.9596, 3.9636 and the un-sparsified ones are in `runs/`
(regenerable with the commands below); they were verified the same way but are dominated.

**`lp_search.py --warm` (the cell-cut search seeded with a point set).**  The one change to
`lp_search.py` is the option `--warm CERT`: the points of `CERT`, scaled from its container to
`S` and snapped to the nearest atom-grid cell, are added as initial orbits (weights are not
reused; every iterate is a fresh LP).  Run: `python3 search/lp_search.py 3.94 0.005 0.005 0.005 0 W1
--iters 30 --seed 0 --verify-threads 8 --warm runs/tight_S1_D1988.txt` — 76 orbits added
(the 592 points of the 3.9437 certificate, snapped to the 0.005 grid).  Result after 30
rounds (47 min, sharing the machine with the runs above): LP 12.99 → 13.75 → **13.13**, never
below the 12.6 threshold at which it calls the verifier; `runs/cert_W1.txt` has total 13.75.
The eroded-cell LP at `η = δθ = 0.005` carries a handicap of ≈ 0.4 % of the side (see
`CEILING.md`), i.e. it is solving something like the exact problem at `s ≈ 3.957`, and the
snapping to the 0.005 grid destroys the fine structure of the seed.  The exact-row loop of
`tighten.py` with column generation is the right tool at these container sizes; the cell
search is not competitive there, warm-started or not.

## Reproduce

```sh
cd verify && cargo build --release && cd ..
mkdir -p runs
# A: weights only, shipped container                       -> runs/tight_A1.txt
python3 search/tighten.py reopt certificates/s12_lower_3.931795.txt A1
# A at other containers (coordinates fixed, D varies)      -> runs/tight_S1_D<D>.txt
python3 search/tighten.py scan certificates/s12_lower_3.931795.txt S1 --Dps 1988,1986,1984,1982,1990
python3 search/tighten.py scan certificates/s12_lower_3.931795.txt S2 --mul 10 --Dps 19850,19845,19855,19842,19848
# B: sparsify at the shipped bound and at 3.9486            -> runs/sparse_B1.txt, runs/sparse_B2.txt
python3 search/tighten.py sparsify certificates/s12_lower_3.931795.txt B1 --budget 11.99 --rounds 8
python3 search/tighten.py sparsify certificates/s12_lower_3.931795.txt B2 --Dp 19855 --mul 10 --budget 11.995 --rounds 8
# C: column generation at 245/62, then sparsify             -> runs/tight_C1.txt, runs/sparse_B3.txt
python3 search/tighten.py reopt certificates/s12_lower_3.931795.txt C1 --Dp 1984 --colgen 40 --cg-want 150
python3 search/tighten.py sparsify runs/tight_C1.txt B3 --budget 11.995 --rounds 8
# C, larger containers.  Each step: colgen seeded with the previous certificate; after the
# column set has settled (LP value flat, ~10-15 rounds) the current support -- written every
# round to runs/tight_<TAG>_probe.txt -- is copied and re-optimised cut-only (much faster than
# letting the big LP converge).  C5b -> C6b/C7b -> (C8: no certificate)
python3 search/tighten.py reopt runs/tight_C1.txt  C5 --Dp 1980 --colgen 60 --cg-want 100   # stop after ~12 rounds
cp runs/tight_C5_probe.txt runs/seed_C5.txt;  python3 search/tighten.py reopt runs/seed_C5.txt C5b
python3 search/tighten.py reopt runs/tight_C5b.txt C6 --Dp 1978 --colgen 15 --cg-want 100
cp runs/tight_C6_probe.txt runs/seed_C6.txt;  python3 search/tighten.py reopt runs/seed_C6.txt C6b
python3 search/tighten.py reopt runs/tight_C5b.txt C7 --Dp 1976 --colgen 15 --cg-want 100
cp runs/tight_C7_probe.txt runs/seed_C7.txt;  python3 search/tighten.py reopt runs/seed_C7.txt C7b
python3 search/tighten.py sparsify runs/tight_C6b.txt B5 --budget 11.995 --rounds 8
python3 search/tighten.py sparsify runs/tight_C7b.txt B6 --budget 11.998 --rounds 8         # -> the shipped 3.9676 file
# checks (every shipped file)
verify/target/release/verify certificates/s12_lower_3.9676.txt 12 6000  32 0
verify/target/release/verify certificates/s12_lower_3.9676.txt 12 12000 32 0
python3 xcheck.py certificates/s12_lower_3.9676.txt 6000 --all --n 12
python3 search/export_points.py export certificates/s12_lower_3.9676.txt -o certificates/s12_lower_3.9676.json --n 12
python3 search/export_points.py --roundtrip certificates/s12_lower_3.9676.txt certificates/s12_lower_3.9676.json
```

`C1` was run with the earlier `max_coverage`-based pricing (see above); re-running the command
with the current file uses the grid pricing and will give a different (in our runs, better)
point set.  The colgen runs were stopped by hand once the support was snapshotted (the snapshot
round is the one thing not fixed by the command line; the `_probe.txt` of every round is the
current support and any of them can be used).  HiGHS's `random_seed` is pinned (`--seed`, default 0); as in
`search/REPRODUCIBILITY.md`, bit-identity holds on one machine and one set of wheels
(numpy 2.4.2, scipy 1.17.0 here).

## Caveats

* All rigor is in the final verifier runs and `xcheck.py`; the LP, the lattice rows, the
  pricing and the reweighted-L1 heuristic only propose weights and points.
* The verifier witnesses are the midpoints of the arrangement cells that attain a bin's
  minimum, printed as floats; the LP row is built at that float placement with atoms within
  `1e-9` of the square boundary not counted (the safe direction).  If a row ever disagreed with
  the verifier the loop would simply not terminate — it never happened.
* `min @ N=6000` is `1 + 2e-6 … 5e-6` by construction (`m_lp`); the files are at their critical
  container to within that margin, so `scale_to_critical.py` gains nothing on them (checked on
  `tight_A1.txt`: critical `D = 1994`, i.e. no gain).  The bound moves by choosing `D` before
  optimising, not by scaling afterwards.
