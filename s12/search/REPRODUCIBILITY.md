# Reproducibility of `lp_search.py`

`search/lp_search.py` is step 1 of certificate generation (the README section
"Reproducing the certificate").  As originally written it stopped on a wall-clock limit, so
two runs never produced the same file.  This note records every source of nondeterminism
that was found, what was done about each, the new command-line options, and a worked
demonstration of two byte-identical runs.

Nothing here affects the *validity* of any certificate: every file `lp_search.py` writes is
checked independently by `verify/` and the Lean development, and the shipped certificate is
pinned by `certificates/SHA256SUMS`.  Reproducibility is about being able to regenerate a
file, not about trusting it.

## Sources of nondeterminism (audit)

| # | source | effect | status |
|---|---|---|---|
| 1 | **wall-clock termination**: `if time.time()-t0 > tlimit: break` | the number of cutting-plane iterations depends on machine load | `--iters N` replaces it with an exact iteration count |
| 2 | **verifier subprocess timeout** (`subprocess.run(..., timeout=600)`): a slow verifier run raised an exception, silently skipping that iteration's exact cuts | iteration-dependent on load | timeout removed under `--iters` |
| 3 | **verifier thread order**: `verify` runs 4 threads, each pushes its witnesses into a shared `Vec` on completion, then does a *stable* sort by weight only (`o.sort_by(|a,b| a.0.cmp(&b.0))`).  Placements with equal weight therefore come back in thread-scheduling order, and `lp_search.py` used the first 4000 of them (and inserted them as LP rows in that order) | which exact-placement cuts are added, and the LP row order (which changes the vertex HiGHS returns) | `read_xsep()` sorts the witness list by the full tuple `(weight, theta, cx, cy)` before use; `--verify-threads` is now free to vary |
| 4 | **stale `runs/xsep_<tag>.txt`** from a previous run with the same tag was read on iteration 0 | run depends on the directory's history | deleted at start-up, together with `snap_<tag>.txt` |
| 5 | **presence of the verifier binary**: if `./verify/target/release/verify` was missing (it is gitignored, must be `cargo build --release`), the `except` swallowed the error and the search ran without exact cuts | different search, silently | logged on the first line of output; `--verifier PATH`, `--no-verify` |
| 6 | **HiGHS internals**: scipy's `method='highs'` forces the serial dual simplex, so there is no thread race, but HiGHS uses `random_seed` (default 0) for its perturbations, and these LPs are highly degenerate: a different seed returns a *different optimal vertex* (demonstrated below: same matrix, `random_seed` 0 vs 1, objective equal to 1e-14, `max|x - x'| = 0.088`) | the exported weights | `threads=1` and `random_seed=<seed>` are passed explicitly; `--seed` sets it, default 0 = HiGHS default |
| 7 | **python/numpy RNG** | none: the search draws no random numbers at all ("stochastic" in the README is a misnomer; it was *time-dependent*) | `--seed` seeds them anyway, so any future randomised step is covered |
| 8 | **dict/set iteration order** (`M.okey`, `M.ckey`, `sorted(set(...))` in `add_orbit`) | none: they are only used for membership tests, and the one iteration is `sorted()`; keys are integer tuples so `PYTHONHASHSEED` is irrelevant | -- |
| 9 | **floating-point reduction order**: `np.cumsum` in `scan_angle`, `np.argsort` tie-breaking, `np.bincount`, and HiGHS's own arithmetic | deterministic for a fixed numpy/scipy build on a fixed CPU; may differ across numpy versions or across CPUs with different SIMD dispatch | not fixable in Python; documented.  Pin wheels (`numpy 2.4.2`, `scipy 1.17.0` were used here) and record the LP with `--dump-lp` |
| 10 | a latent bug: `best=(val,x,m)` kept a reference to the live model `m`, which keeps acquiring orbits after the snapshot; `write_cert(m,x,...)` then indexes `x` with a longer `own` array | `IndexError` (or wrong weights) whenever the best iterate was not the last one | iterates are now frozen in a `Snap` object |

With 1-6 addressed, a run under `--iters` is a pure function of
`(S, FINE, ETA, DT, TAG-independent options, --seed, --iters, verifier build, numpy/scipy/HiGHS versions, CPU)`.

## Command line

```
python3 search/lp_search.py S FINE ETA DT TLIMIT TAG [--iters N] [--seed N] [--dump-lp PATH]
                            [--verifier PATH] [--verify-threads N] [--no-verify]
                            [--runs DIR] [--highs-threads N]
python3 search/lp_search.py --resolve DUMP [--seed N]
```

The six positional arguments are exactly the old ones, so the README command
`python3 search/lp_search.py 3.92 0.005 0.005 0.005 7200 mytag` still runs with the old
(wall-clock) behaviour.  Options:

* `--iters N` -- perform exactly `N` LP solves (cutting-plane rounds), ignoring `TLIMIT`.
  The only other exit is finding a certificate of weight `< 11.999`, which is deterministic.
* `--seed N` -- seeds `random`, `numpy.random`, and HiGHS (`random_seed`).
* `--dump-lp PATH` -- archive the LP whose solution was exported (format below).
* `--resolve DUMP` -- read a dump, re-solve it with the same HiGHS settings, and report
  whether the recorded solution is reproduced exactly.
* `--verifier`, `--verify-threads`, `--no-verify`, `--runs`, `--highs-threads` -- as named.

Output: `runs/cert_TAG.txt` (certificate), plus the sha256 of it on the last line of the log.

## LP dump format (`--dump-lp`)

Plain text, one item per line, all floats as C99 hex literals (`float.hex()`), so the file
rebuilds the exact doubles the solver saw:

```
LPDUMP 1
META key value            # s, fine, eta, dt, tag, iters, seed, highs options, versions,
...                       # certificate path and sha256
PROBLEM minimize c.x  subject to  A x >= 1,  x >= 0
SIZE nrows ncols nnz
OBJ ncols                 # then ncols lines: c_j = orbit size (1,2,4 or 8)
ROWS nrows                # then nrows lines: u0 u1 theta  (cell centre in the rotated
                          #   frame u = R_{-theta} c, and the bin angle)
COLS ncols                # then ncols lines: k x_1 y_1 ... x_k y_k (the orbit's atoms)
ENTRIES nnz               # then nnz lines: row col value  (value = number of atoms of
                          #   orbit `col` guaranteed inside every square of cell `row`)
SOLUTION objective
X ncols                   # then ncols lines: primal solution
DUAL nrows                # then nrows lines: dual solution (row marginals, >= 0)
END
```

`lp_search.read_lp_dump(path)` returns all of this as a dict; `--resolve` uses it.  The
certificate is a function of `X`, `COLS` and `s` only (`write_cert`), so the dump plus the
seed is enough to regenerate the certificate without re-running the search, and the LP can
be handed to any other solver.

## Demonstration

Small configuration (cell/atom spacing 0.02, angle bin 0.02, 3 rounds, verifier unused
because the LP value never drops below 12.6 at this coarseness), two runs, different tags:

```
python3 search/lp_search.py 3.92 0.02 0.02 0.02 0 a1 --iters 3 --seed 1 --no-verify --dump-lp runs/lp_a1.txt
python3 search/lp_search.py 3.92 0.02 0.02 0.02 0 a2 --iters 3 --seed 1 --no-verify --dump-lp runs/lp_a2.txt
```

```
b9a735d536b833a6a3dfbda65f1d200752cf46420a9d7422b5d6a4544781621d  runs/cert_a1.txt
b9a735d536b833a6a3dfbda65f1d200752cf46420a9d7422b5d6a4544781621d  runs/cert_a2.txt
```

The LP dumps are identical apart from the `META tag`/`META certificate` lines
(`grep -v ^META runs/lp_a?.txt | sha256sum` = `758695c1…d1eb` for both).  Re-solving the
dump with the same seed reproduces the primal solution bit for bit; with a different seed it
does not:

```
$ python3 search/lp_search.py --resolve runs/lp_a1.txt --seed 1
  recorded objective  = 15.999999999999993
  re-solved objective = 15.999999999999993  (status 0)
  max |x - x_recorded| = 0.000e+00
  x identical: True
$ python3 search/lp_search.py --resolve runs/lp_a1.txt            # HiGHS random_seed = 0
  re-solved objective = 16.000000000000004  (status 0)
  max |x - x_recorded| = 8.753e-02
  x identical: False
```

Note on coverage: at spacing 0.02 the LP value stays at 16 (> 12.6), so the verifier
feedback path (items 2-5) is not exercised by this demonstration; it is guarded by the
canonical sort in `read_xsep()` and the start-up cleanup, which are deterministic by
construction.  At spacing 0.01 a single cutting-plane round already takes more than 15
minutes on a 32-core machine, so a demonstration that reaches `LP < 12.6` (and hence calls
the verifier) was not run here.  Each run above took ~12 s.

## What is *not* reproduced

* The shipped 788-point certificate.  It was produced by the original time-limited script
  on an unrecorded machine, with the thread-order and stale-file effects above in play; its
  LP input was not archived.  It remains pinned by hash and verified on every CI run.
* Cross-machine bit-identity.  Items 6 and 9 mean a different numpy/scipy/HiGHS build, or a
  CPU with different SIMD dispatch, may return a different (equally valid) optimal vertex.
  For an archival-grade record, ship the `--dump-lp` file next to the certificate.
