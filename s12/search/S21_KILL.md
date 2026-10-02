# The `s(21) = 5` kill test: `nu_f^closed(5) >= 20.6478`, and the value sits `~0.35` below 21 (2026-09-23)

Code: `search/packing_dual.py` (float column generation, unchanged apart from one `analyse` bin),
`search/dual_exact.py` (exact snap + certification; new `--no-lp`, `check --stream`, `--n`), `search/cover4_cg.py`
(exact restricted master, unchanged).  Certificate: **`search/s21_nuf5_exact_support.txt`**
(= `runs/dual_exact_s5killXE_support.txt`).  Labels: **certified** (Python integers / `Fraction`, a
theorem), **heuristic** (float LP), **estimate** (extrapolation).

## 0. Answer

> **The route is not dead, but the margin is small.**  `nu_f^closed(5) >= L = 10323890641/499999999 = 20.647781323` (**certified**).
> The float LP converged (`M = 1`, no violated arrangement vertex, float-checked over all 10.6 M
> vertices of the whole container) at **`20.64778`** (**heuristic**).  Continuing column generation
> no longer moved it: the LP mass stayed at `20.6477-20.6483` for the last 20 rounds.  So
> `nu_f^closed(5) ~ 20.65` (heuristic, at most a few `10^-2` higher).
>
> `L < 21`, so the kill test does **not** fire.  Weighted closed covers of weight `< 21` are not ruled out.
> What *is* proved is that any such cover sits within
> **`21 - L = 0.3522`** of `nu_f`, i.e. it has overhead `< 0.3522` over the packing value
> (**certified**: `W >= nu_f >= L`).  The certified `m = 4` cover had overhead `12.9560 - 12.2688 = 0.687`, about
> twice that.  In the `S21_OVERHEAD.md` model (`W ~ 1.05 (nu_f + 0.11)` with the `m = 4`
> recipe, **estimate**), `nu_f(5) = 20.65` puts `m = 5` in its **"marginal" band `20.5-20.7`**.  The
> combined validity x checker factor must be `<= 21/(20.648 + 0.11) = 1.0117`, against `1.05`
> shipped at `m = 4` and a "realistic floor" of `1.01-1.02` there.  The band "plausible only if
> `nu_f(5) <~ 20.5`" is now excluded: that is **certified**, since `L = 20.648 > 20.5`.

| quantity | `m = 4` (`COVER4.md`) | `m = 5` (this note) | status |
|---|---|---|---|
| `n` | 12 | 21 | |
| `nu_f` lower bound `L` | `12.2688039` | **`20.6477813`** | certified |
| `nu_f` float LP, converged (`M = 1`) | `12.268806` | **`20.647782`** | heuristic |
| `L - n` | **`+0.2688`** (killed `s(12)` at `t = 4`) | **`-0.3522`** | certified |
| `m^2 - nu_f` | `3.731` | `4.352` | heuristic |
| best certified cover on record | `12.9560` (overhead `0.687`) | none | certified / - |
| heuristic cover LP | `12.39-12.42` | `20.75-20.9` (unconverged, `FAMILY.md` §3.1; consistent: cover `>= nu_f`) | heuristic |
| previous `nu_f(5)` floor | - | `16.053` (`FAMILY.md` §3.2) | bound |

## 1. Statement proved

Let `C = [0,5]^2`.  `search/s21_nuf5_exact_support.txt` lists 748 poses `(p, q, cx, cy, m)`.  Each has a rational
rotation `theta = 2 arctan(p/q)` (`Q = 10^7`), a rational centre (`10^-8`, clamped exactly into
`[w/2, 5 - w/2]`) and mass `m/8` on each of its 8 dihedral images.  That is 5984 closed unit squares, each checked
exactly to lie in the closed container.  For this measure `mu`:

    mass(mu) = 10323890641/500000000 = 20.647781282          (exactly)
    M        = max_{x in C} mu({Q : x in Q}) = 499999999/500000000   (exactly; attained at (1, 1))
    L        = mass/M = 10323890641/499999999 = 20.647781323296...

So `mu/M` is a feasible fractional packing, `nu_f^closed(5) >= L`, and by weak duality no weighted point set
of total weight `< L` covers every closed unit square in `[0,5]^2` (**certified**).  The machinery is the same as
`COVER4.md`/`DUAL_EXACT.md`: arrangement vertices (square corners plus exact pairwise edge
intersections), integer containment, and upper semicontinuity of `cov`.  Only the side changed (`--t 5`).

| check (from the support file alone, no LP) | vertices | `(vertex, square)` pairs | `M` | `L` | time |
|---|---|---|---|---|---|
| `F = {0 <= x <= y <= 5/2}` (D4 reduction) | 1,311,398 | 456,658,626 | `499999999/500000000` | `10323890641/499999999` | 92 s |
| whole container (`--full --stream`, no reduction) | 10,482,644 | 3,649,873,072 | `499999999/500000000` | `10323890641/499999999` | 668 s |

## 2. How the measure was found (heuristic; nothing load-bearing)

| run | what | result | time |
|---|---|---|---|
| `s5killA` | `packing_dual.py 5`, warm = `E1` (t=4) support rescaled x1.25, seeds 0.08 / 7.5 deg, 7 threads | LP mass `20.45-20.56` over 5 rounds; LP grew to 50-75 M nz, `2227 s` per LP; **killed** | ~2 h |
| `s5killB` | same without warm start, seeds 0.06 / 5 deg | 118-177 M nz, `3281 s` per LP, mass `20.47-20.56`; **killed** | ~1.5 h |
| exact CG `G1-G5` | `cover4_cg.py cg --t 5` over pools of the A/B/C/D per-round supports | exact `L` `20.1325 -> 20.4721 -> 20.5243 -> 20.5821 -> 20.5912` | ~2.5 h |
| `s5killC` | warm = exact CG support (293 poses), coarse seeds 0.15 / 15 deg, rows 0.05, 6 threads | **converged** round 43: `20.63206`, `M = 1`, full container 6.96 M vertices | 5647 s |
| `s5killD` | warm = CG support, seeds 0.1 / 10 deg, `--smooth 0.7` | **converged** round 31: `20.62881` | 5612 s |
| `s5killE` | warm = `C` final support (594 poses), 7200 s colgen + polish | **converged** round 36: **`20.64778`**, `M = 1`, full 10.59 M vertices | 7395 s |
| `s5killF` | warm = `D`, seed file = `C` | noisy, 470-1475 s LPs; **killed** (not needed) | ~1 h |

Wall-clock for the whole task: about 5.5 h on CPUs 0-13 (2026-09-22 23:36 to 2026-09-23 05:08).

Exact side:

| step | `L` (exact) |
|---|---|
| `C` round-0 support (178 poses), `build` with LP polish | `20.132471729` |
| exact restricted master (`cover4_cg.py`), best | `20.591153415` |
| `C` final (594 poses), `build` with LP (569k rows x 594 cols, 99.6 M nz, 891 s) | `27509471104/1333333333 = 20.632103333` (F and full checks both pass) |
| **`E` final (748 poses), `build --no-lp`** | **`10323890641/499999999 = 20.647781323`** |

*What worked at `t = 5`.*  Big seed grids make the float LP unusable here: at `t = 4` it was `<= 26 M` nz, at `t = 5` it reached `50-177 M` nz with 10-55 min
per LP.  What worked was to warm-start from a good *support*, with a coarse seed grid.  Then the
LP stays at `6-15 M` nz and a round takes 1-3 min.  Once the float run has converged (`M = 1`), snapping it with
`Q = 10^7` loses only `~4e-8`: the "before polish" exact `M` is already `0.999999998`.  So the
exact LP polish is unnecessary, which is what `--no-lp` is for.  At 774 poses the LP polish
reached 81 GB RSS and was killed.

## 3. The measure (748 poses, mass 20.6478)

Concentration: 25% of the mass on 2 poses, 50% on 16, 75% on 125, 90% on 328, 99% on 616.

| pose `(cx, cy, theta)` | mass | what it is |
|---|---|---|
| `(0.5, 0.5, 0)` | 3.4717 | the four corner squares, `0.868` each (`0.813` at `t = 4`) |
| `(0.5, 1.5, 0)` (two poses) | 2.1654 + 0.2741 | wall squares next to a corner |
| `(1.5, 1.5, 0)` | 0.7113 | inner-ring diagonal squares |
| `(0.5, 2.5, 0)`, `(1.5, 2.5, 0)` and near-0 deg neighbours | 0.56, 0.55, 0.52, 0.46, ... | mid-wall and inner-ring axis-aligned squares |
| `(1.576, 3.424, 9.5 deg)`, `(1.211, 1.516, 41.5 deg)`, `(1.245, 1.563, 31.1 deg)`, ... | 0.28, 0.16, 0.11, ... | tilted mixture |

Mass by angle: `[0,5)` **11.59 (56%**, against 70% at `t = 4`), `[5,30)` 3.72 (18%), `[30,45]` 5.35 (26%).  Axis-aligned (`< 1 deg`): 51%.
By distance of the canonical centre to the nearest wall: `< 0.52` 9.25 (45%), `0.52-0.8` 3.05 (15%), `1.2-1.6`
5.45 (26%), `1.6-2.0` 2.16 (10%), `2.0-2.5` 0.73 (4%).  So there is proportionally *more* tilted
mass than at `t = 4`: the fractional mixture fills the bigger interior.  The rigid frame (corners and walls) looks like the `t = 4` one.

## 4. Code changes (parametrised, backwards compatible)

* `dual_exact.py build --no-lp`: keep the float masses (rounded down to `1/DM`, rescaled by `1/M` if
  `M > 1`) instead of the HiGHS polish.  The LP moved into `polish_lp()` unchanged.  The default path
  is unchanged: `build` with no `--t` still writes a `runs/dual_exact_3.99_support.txt` that is
  **byte-identical** to the certificate (re-run as `--tag s5killREG399`, `cmp` identical).  `check` on
  `search/cover4_exact_support.txt --t 4` still gives `24537607710/1999999999`
  (`runs/dual_exact_s5killREG4_check_F.log`).
* `dual_exact.py check --stream`: evaluate the max coverage per vertex block in the workers
  instead of keeping all incidence lists (same `cov_worker` containment test).  The non-streamed full
  check at `t = 5` (3.65 G pairs) drove available memory below 12 GB and was killed by a guard.  The streamed one peaked at a few GB (system memory stayed at `~11 GB` used).
  Both modes give identical output on the `F` check (`--tag s5killXEs`).
* `dual_exact.py --n N` (build and check, default 12): the target that the log line compares `L`
  with.  It is display only, and was hard-coded `12`.  *(2026-10-01: for `check`, `--n N` now also
  asserts: the exit status is 0 iff `L > N` exactly, 1 otherwise; with no `--n` nothing is asserted and
  the log line compares with 12.  So `--n 20` asserts the certified bound, while `--n 21` exits 1 —
  that is the kill test not firing.)*
* `packing_dual.py analyse`: adds a `[2.0, t/2]` centre-distance bin when `t > 4`.  Before, centres
  farther than 2 from every wall were silently dropped from the table.  The output at `t <= 4` is unchanged.
* No hard-coded side 4 found in `packing_dual.py`, `dual_exact.py` or `cover4_cg.py`.  `cover4_cg.py`'s `--t`
  default stays `4`; it passes `--t` through.  `closed4.py` still mirrors the literature columns
  at a hard-coded `4.0` (`FAMILY.md` §3.1).  It is not on this pipeline and was not touched,
  because the parallel `s5conv` cover-LP run uses it.

## 5. What is heuristic

The value `20.64778` as an estimate of `nu_f(5)` is heuristic.  It is the optimum over the columns that
column generation found.  From `C` to `E` it rose by `+0.0157`, then stayed flat for about 20 rounds.  Pricing
kept finding columns with reduced cost (`~250` per round) that did not raise the value.
The true `nu_f(5)` is `>= 20.6478` (certified) and probably `<= ~20.67` (heuristic: only the flat tail
bounds it).  An upper bound on `nu_f(5)` needs a certified cover, and none exists at `m = 5`.  Everything that
*chose* the measure is heuristic: float LPs, pricing, warm starts, the exact CG batching.  The
certificate is the support file, and `check` re-derives `mass`, `M`, `L` from it with integers only.

## 6. Reproduce

```sh
cd verify && cargo build --release && cd ..
# float colgen: warm start from an exact-CG support (any decent t=5 support works), then continue from C
python3 search/packing_dual.py 5 s5killC --warm runs/s5kill_rounds/warmC.txt --seed-pitch 0.15 --seed-dth 15 \
    --row-pitch 0.05 --time 5400 --polish-time 3600 --rounds 150 --N 1000 --cg-want 1500 --threads 6 --final-full
python3 search/packing_dual.py 5 s5killE --warm runs/dual_s5killC_support.txt --seed-pitch 0.15 --seed-dth 15 \
    --row-pitch 0.05 --time 7200 --polish-time 3600 --rounds 200 --N 1000 --cg-want 1500 --threads 6 --final-full
# exact snap + certification, no LP polish (~2.5 min)
python3 search/dual_exact.py build --t 5 --n 21 --no-lp --tag s5killXE --src runs/dual_s5killE_support.txt \
    --Q 10000000 --Dc 100000000 --procs 8
# independent re-certification from the support file alone
# --n 20: exit 0 iff L > 20 (the certified claim);  --n 21 would exit 1 (L < 21: the kill test does not fire)
python3 search/dual_exact.py check search/s21_nuf5_exact_support.txt --t 5 --n 20 --procs 8          # F, ~1.5 min
python3 search/dual_exact.py check search/s21_nuf5_exact_support.txt --t 5 --n 20 --full --stream --procs 8   # full, ~11 min, a few GB
```

Re-checked 2026-10-01 after the `dual_exact.py` hardening (explicit pool initializer, asserting `--n`,
strict support parser): `check search/s21_nuf5_exact_support.txt --t 5 --n 20 --stream --procs 4` gives the
same 1,311,398 vertices, 456,658,626 pairs, `M = 499999999/500000000`, `L = 10323890641/499999999`, exit 0
(138 s on 4 processes).

`warmC.txt` is the exact-CG support `runs/dual_exact_cgs5killG32_support.txt` converted to float
format.  The CG was bootstrapped from `s5killA` round 0: `dual_exact.py build --t 5` on it, then
`cover4_cg.py pool`/`cg --t 5` over the per-round supports in `runs/s5kill_rounds/`.  The converged
supports `runs/s5kill_rounds/{C,D,E}_final_support.txt` make this step unnecessary.

## 7. Files

| file | what |
|---|---|
| `search/s21_nuf5_exact_support.txt` | **the certificate** (748 exact poses, header `# t = 5/1`) |
| `runs/dual_exact_s5killXE.log` / `.json`, `_check_F.log`, `_check_full.log` | exact build and checks |
| `runs/dual_s5killE.log` / `.json` / `_support.txt` | converged float run, `20.64778` |
| `runs/dual_s5killC.*`, `runs/dual_s5killD.*` | converged float runs, `20.63206` / `20.62881` |
| `runs/dual_exact_s5killXC_support.txt` (+ `_check_F.log`, `_check_full.log`) | earlier certificate `20.632103333` (with LP polish) |
| `runs/dual_s5killA.*`, `runs/dual_s5killB.*`, `runs/dual_s5killF.*` | killed float runs (too slow / noisy) |
| `runs/cg_s5killG1.out`, `runs/cg_s5kill_loop*.out`, `runs/dual_exact_cgs5killG*` | exact restricted-master steps |
| `runs/s5kill_rounds/` | per-round float supports, warm-start files, helper scripts |
