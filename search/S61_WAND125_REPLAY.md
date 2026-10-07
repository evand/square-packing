# s(61) = 8: independent replay of wand125's point-only cover (2026-09-30)

**Result.** wand125's D4-invariant point measure on `[0,8]²` (15,193 points, total
`8584985072679551 / 2^47 = 60.99998780001… < 61`) is certified by both of our separately written exact
point checkers. `zeromargin.py` certifies all 12,800 roots of the D4 region (two of them only after a
depth-30 rerun). `zmcheck --d4` with the opt-in branch order `ZM_MIXPAIR=1` prints `VERIFIED-D4`, all
6,400 roots, 0 uncertified. With its default branch order, `zmcheck` stalls at the interior tile germs, the
same way it did on the first s(32) cover. We also ran `zmx2` as a reference: `--d4` reproduces their box
count exactly (800,042), and the unreduced `--full --sym-atoms` sweep is clean too. So every closed unit
square in `[0,8]²` captures mass ≥ 1, which gives `s(61) ≥ 8`. With the 8×8 grid, `s(61) = 8`.

This value already follows from our `s(60) = 8` (`certificates/s60/`). What is new here is a route that
uses points only, and a cover built by someone else.

Records: `search/s61_wand125/` (`MANIFEST.txt` lists every hash and command). There is no
`certificates/` bundle.

## 1. Provenance

| item | value |
|---|---|
| source | `github.com/wand125/square-packing-bounds`, commit `f8846cec9661773dbd0cc7cbeee7d01ddb12a2b8`, `point_n61_L8/`, fetched with `gh api …/contents/…?ref=<commit>` |
| their `cover.txt` | sha256 `bcb66c7910ef844ce8d39c423d7a5a419bd530a344689331f7665ba6971d0374` (equals `provenance.json`'s `cover_sha256`) |
| our `s61_wand125/s61_wand125_cover_8.txt` | sha256 `bcb66c7910ef844ce8d39c423d7a5a419bd530a344689331f7665ba6971d0374`, **byte-identical**: their file is already in `certificates/FORMAT.md` plain format (`8 1 / D = 2000 / W = 2^49 / m = 15193 / X Y w`). `convert_check.py` (written here) parses every token as an integer, re-emits the file canonically, round-trips it and gets the same bytes |
| their scripts | `check_cover.py`, `verify.sh`: read, **not executed** (sha256 in `MANIFEST.txt`) |
| their claim | checked only with our `zmx2` at our commit `6e1223cf` (`zmx2.rs` sha256 `6b7f0f79…`): `zmx2 cert cover.txt --d4 --pair-points` → `VERIFIED-D4`, 6,400 roots, 800,042 boxes, max depth 30, 0 uncertified, 0 capped |

## 2. Exact checks (`convert_check.py`, integers and `Fraction` only)

* Header: `s = 8/1`, `D = 2000`, `W = 562949953421312 = 2^49`, `m = 15193`. There are exactly `m` point lines, no trailing data and no non-integer token.
* All points lie in `[0,8]²`. In fact they lie in `[81/500, 3919/500]²`, so **no point is on the container boundary**. All weights are `> 0` (none is zero; the smallest is `28765412/2^49`). No coordinate appears twice.
* **D4 invariance:** exact under `x ↦ 8−x`, `y ↦ 8−y` and `x ↔ y`. All three checkers also re-check this themselves.
* **Total** `= 34339940290718204 / 2^49 = 8584985072679551 / 140737488355328 = 60.99998780001351…`. Then `61 − total = 1716995457 / 2^47 ≈ 1.22·10⁻⁵ > 0`. This matches their `provenance.json` exactly.
* Weight on the lines `x ∈ ℤ` or `y ∈ ℤ`: 39.227 of 60.99999. On `(½)ℤ` lines: 41.727.

## 3. Runs

All runs were pinned with `taskset -c 6-9` (4 threads or processes). Repo HEAD was `2bf33bc`. "CPU" is
the sum of per-root times for `zeromargin.py` and `zmcheck`, and user time for `zmx2`.

| checker | settings | roots | boxes | max depth | uncertified | CPU | verdict |
|---|---|---|---|---|---|---|---|
| `zeromargin.py` (`zm_d4_sweep.py`) | depth 24, pitch 1/10, 8 u-bins of 1/16 on `[0,½]`, CHAIN+ADM, clip | 12,800 | 92,462 | 24 | 24 boxes in 2 roots | 0.84 h | open at depth 24 |
| ″ `--deepen 30` (those 2 roots) | depth 30 | 2 | 1,702 | 30 | 0 | 0.01 h | both certified |
| ″ summary | | **12,800 / 12,800** | 92,924 | 30 | **0** | **0.85 h** | **`D4 RECHECK CLEAN`** |
| `zmcheck --d4` (default branch order) | `--branch-cap 640 --node-cap 4000000 --depth 22`, pitch 1/10, 4 u-bins of 1/8 | 2,582 of 6,400 done | — | — | 0 so far | ≤ 1.45 h | **stopped**: 3 roots still running (see below) |
| `zmcheck --d4`, `ZM_MIXPAIR=1` | same | **6,400 / 6,400** | 106,554 | 22 | **0** | **1.07 h** (964 s wall) | **`VERIFIED-D4`** |
| `zmx2 cert --d4 --pair-points` | binary `ed31d3ee…` (source `92a4cfe8…`) | 6,400 | **800,042** | 30 | 0 | 0.13 h | `VERIFIED-D4` |
| `zmx2 cert --full --pair-points --sym-atoms` | same; whole `[0,8]²`, both passes, no D4 fold | 51,200 | 6,236,672 | 30 | 0 | 1.23 h | `VERIFIED` |

The whole replay used about 4.8 CPU-hours. Full-sweep projection for `zmcheck`'s default branch order (not
run): on the first s(32) cover, each germ root ran for about 6,000 s and still left boxes open. The D4
region of `[0,8]²` contains 9 interior tile germs (`(k+½, l+½)`, `k, l ∈ {1,2,3}`). That puts a default-order
full sweep in the tens of CPU-hours, which is why it was stopped.

**`zmcheck`, default order.** The roots still running after about 20 minutes were
`x∈[1.4,1.5], y∈[2.4,2.5]` and `x∈[1.4,1.5], y∈[2.5,2.6]`, both in u-bin 0. A third, `x∈[1.5,1.6], y∈[3.9,4.0]`
u-bin 0, had started recently. These are the cells around the interior tile germ `(1.5, 2.5)`: the same
germ type, and in the s(32) case the same cells, that this order failed to finish on the first s(32) cover
(`certificates/s32/zmcheck_d4/SUMMARY_candidate.txt`). We stopped the run and did not let these roots reach
the depth limit, so they are "did not finish", not "uncertified". A 4-root pilot with `ZM_MIXPAIR=1` on that
germ certified all four in 16–34 s each. The full `ZM_MIXPAIR` sweep then certified everything. Its slowest
root took 79 s (`x∈[2.5,2.6], y∈[1.4,1.5]`, u-bin 0, another germ cell). `ZM_MIXPAIR` changes only the
order in which branch families are tried (`verify2/src/main.rs` around line 1392). It is not a new
primitive, so the `VERIFIED-D4` verdict rests on the same soundness argument as a default-order run.
Since `zmcheck` alone closes every root, `zeromargin.py` was not needed to close anything for it. The two
checkers certify the region independently.

**`zeromargin.py` at depth 24.** The 24 open boxes all have `u ∈ [53/128, 1697/4096]`, i.e. `θ ∈ [44.985°, 45.009°]`,
which straddles 45°. Their centres lie along the segment `y − x ≈ 2.332`, `x ∈ [1.344, 1.415]`, inside roots
`(7/5..3/2, 37/10..19/5)` and `(13/10..7/5, 18/5..37/10)` of u-bin `[3/8, 7/16]`. A float scan of that
family (`θ ∈ [44.9°, 45.1°]`, ±0.01 across the line) gives a minimum captured mass of **1.0112**, so this is
not a tight pose. It looks like the subdivision struggling where the bin contains 45°, which is where the
pose width `w(θ) = |cos θ| + |sin θ|` peaks. Both roots are certified completely at depth 30, under the
same rule the s(32) bundle used for its one germ root.

**`zmx2` vs their run.** Our `zmx2.rs` is newer than the one they pinned. Ours includes the opt-in
`--sym-atoms`, which was not used in the `--d4` run. Even so, the `--d4 --pair-points` census matches
their reference run exactly: 6,400 roots, 800,042 boxes, max depth 30. Their `roots_sha256` covers a log
whose line order and `ms` fields depend on the threads and the machine, so it cannot be compared to ours.

## 4. Where the margin is (float, diagnostic only)

`zmx2 fscan` (pitch 0.01, 50 u-bins over the D4 region, 5,747,007 poses) gives a minimum captured mass of
**1.00648** at `(3.78, 2.27)`, θ = 51.3°. The next lowest are 1.0081 at `(3.32, 3.19)`, 41.6°, and 1.0086 at
`(2.27, 3.78)`, 38.6°. Four sampled poses fall below 1.01. This fits the README's construction: the LP
solution was scaled by `61/M·(1−2·10⁻⁷) ≈ 1.0117`, which gives a uniform margin of about 1 %, and extra
LP rows removed the remaining contacts. We found **no zero-margin pose**. None of the three exact checkers
needed a zero-margin primitive anywhere it could not handle it. No point is on the container boundary.
The corner and wall poses are certified by the usual primitives (P1/ADM), as for s(32). This scan samples
the region and is not exhaustive, so it is evidence and not a proof. The proof is the exact sweeps in §3.

## 5. What is independent of what

* **The cover is theirs.** wand125 built it by lifting their s(45) point cover, then running an LP iterated
  against counterexamples from `zmx2` and from an exact all-centre check of their own. So the cover was
  tuned against `zmx2`, though not against `zeromargin.py` or `zmcheck`. Our file is theirs byte for byte.
  The format checks and the exact total were recomputed by our own `convert_check.py`, which shares no code
  with their `check_cover.py`. We read their script but did not run it.
* **`zeromargin.py` (Python, `Fraction`/integers) and `zmcheck` (Rust, `verify2/src/main.rs`) were written
  separately and share no code** (`certificates/s32/README.md`). Their root grids differ: 1/16 against
  1/8 u-bins, so 12,800 against 6,400 roots. So do their primitives: ADM/CHAIN against ADM/DISJ. Each
  re-checks D4 invariance on its own. **These two runs are the independent confirmation.**
* **`zmx2`** is the checker wand125 used. Its author was allowed to read `zeromargin.py`, so `zmx2` and
  `zeromargin.py` share the same point-test lineage (the closed-square, exact-core tests). Our `zmx2` runs
  are therefore a reproduction of their claim, plus the stronger no-fold `--full --sym-atoms` sweep. They
  are not a third independent check.
* Not machine-verified: the checker programs themselves. There is no Lean transcription of this cover. The
  `S32Lower` route would apply in principle, but at 15,193 points and side 8 it would be a large build.

## 6. Issues found

1. `zmcheck`'s default branch order does not finish the interior tile-germ roots within about 20 minutes
   each. This is the same behaviour as on the first s(32) cover. `ZM_MIXPAIR=1` fixes it. If this check is
   ever put into a bundle's `verify.sh`, it should pass `ZM_MIXPAIR=1`.
2. `zeromargin.py` needs depth 30 for two roots at θ ≈ 45°, a family with float margin ≥ 1.011. This
   points to a box-splitting inefficiency where a bin straddles 45°, not to a weakness in the cover. It may
   be worth a look: s(32) needed depth 30 only at a θ = 0 germ.
3. Nothing suspicious in the data. There are no duplicate or zero-weight points, nothing on the container
   boundary, and the total and D4 invariance match their claims exactly. Of the mass, 64 % sits on integer
   grid lines, which is expected of a lifted grid cover. No zero-margin pose was found.
