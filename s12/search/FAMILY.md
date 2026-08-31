# The m^2-4 family and rung 2 (task H; 2026-08-30)

Code: `search/zeromargin.py` (exact checker, now vectorised, task H step 1), `search/closed4.py`
(heuristic cover LP, container side `s` a free parameter already), `search/rung2_close.py` (new:
a stress-driven separation loop closing the gap between a converged LP cover and the exact
checker's zero-margin requirement), `search/nu_f.py` (rigorous, float-checked lower bound on
`nu_f(s)`, already parametrised by the container side read from the certificate file). Every
number below is labelled **certified** (exact `Fraction` arithmetic, or float-checked with an
explicit `1e-9`-scale margin the same way `CLOSED4.md` §7 already accepts as rigorous),
**heuristic** (a sampled LP or a local search — a value, not a bound), or a **bound** (with its
direction stated).

## 1. Table: `m -> cover LP value, certified cover, certified packing mass, n-1`

| `m` | `n = m^2-4` | cover LP value (heuristic) | certified closed cover | certified packing mass (rigorous lower bound on `nu_f`) | `n-1` |
|---|---|---|---|---|---|
| 4 | 12 | 12.4175 (`runs/closed4_best.txt`, `CLOSED4.md`); closing loop in progress, see §2 | not yet (§2) | **10.68** (`CLOSED4.md` §7, float-checked, `nu_f(4) >= 10.68`) | 11 |
| 5 | 21 | see §3 (in progress) | none attempted (out of scope per the brief: bracket only) | see §3 | 20 |

(`m = 6, 7` not started; the brief gates them on m = 4, 5 being finished or the time budget
being spent — see §4.)

## 2. Rung 2, `m = 4`, `W < 13`

**Step 1 (oracle) — done.** `search/zeromargin.py`'s float pre-filter is now vectorised with
`numpy` over the whole point set (§7 of `search/ZEROMARGIN.md`); rung 1 (Friedman's 14 points) is
byte-identical before/after; a `--depth 6` exact check of the 1,972-point `runs/closed4_best.txt`
cover now takes 74 s on 8 processes (was 575 s on 4) for identical leaf counts. New CLI modes:
`--oracle FILE` (uncertified-box corner/centre poses as candidate LP rows) and `pose CERT --cx
--cy --u` (exact, non-adaptive captured-weight check at one rational pose — turns a float-found
candidate violation into a rigorous yes/no).

**Step 2 (separation loop) — status: TODO after this run finishes; see the log for the final
numbers.** `CLOSED4.md` left a documented 0.7% gap: the LP converged on its own sampled row set
at `12.4174`, but a denser stress scan (pitch 0.005) found a true minimum of `0.99267` at
`(3.3965, 1.4235, 76.41 deg)` (equivalently `13.6 deg` by symmetry) -- an unconverged tilted wall
pose, not a structural obstruction (`ZEROMARGIN.md` §5). `search/rung2_close.py` runs the
separation loop the brief asks for over the *same* 1,972-point column set: solve the LP, stress-
test at pitch 0.003 to find the actual worst poses, add those (+ a local polish around them) as
hard rows, resolve, repeat. This is a genuine cutting-plane loop against the residual gap (each
round is blind exactly where the previous one was), not another generic refinement.

*Placeholder -- filled in once `runs/closed4_rung2c.log` finishes or is judged stalled:*

- final LP value: `TBD`
- final stress min (pitch 0.003, local polish): `TBD`
- honest cost (`total / stress min`): `TBD`
- exact checker verdict on the final point set at increasing depth: `TBD` (uncertified box count,
  locations, margins if any remain -- see the log)
- whether a tilted tight family (`ZEROMARGIN.md` §4 item 4 / §5(ii)) is visible in the converged
  cover, and how it was or was not handled

If the loop stalls above weight 13 or below stress min 1, that is reported here as the residual
gap (location, margin), per the brief -- not declared as near-success.

## 3. `m = 5`, `n = 21`, `W < 21`

**Cover LP (heuristic).** `search/closed4.py run --s 5 --no-literature ...` (the `s = 4`-specific
Bentz/Nagamochi literature columns do not generalise -- their `d4` images are computed at a
hardcoded container side 4.0 in `closed4.py`, so they were excluded rather than silently
mis-mirrored; the grid-lattice and dual-pricing columns are fully general in `s` already).
Trajectory (heuristic; the value is a *relaxation over its own sampled rows*, so it is a lower
bound on the cheapest cover *of those columns*, and rises monotonically as rows are added, exactly
as at `m = 4`):

| iteration | LP value | rows | support | lattice min | t (s) |
|---|---|---|---|---|---|
| 0 | 16.000 | 18300 | 5 | 0.000 | 11 |
| 1 | 19.909 | 26862 | 138 | 0.537 | 54 |
| 2 | 20.829 | 32841 | 144 | 0.901 | 120 |
| 3 | 20.908 | 34531 | 201 | 0.935 | 196 |
| 4 | 20.892 | 35567 | 245 | 0.979 | 323 |
| ... | ... | ... | ... | ... | ... |

(Full table and final value: `runs/closed4_s5.log` / `.json`; see the placeholder below for
where it plateaued.)

*Placeholder -- filled in once `runs/closed4_s5.log` finishes or is judged stalled:*

- final heuristic LP value at `m = 5`: `TBD` (vs `n - 1 = 20`)
- certified packing mass: `python3 search/closed4.py lower runs/closed4_s5_best.txt --tag s5lower`
  (calls `nu_f.lower_from_cert`, which is already fully parametrised by the container side read
  from the certificate -- no `m = 4`-specific code found in `search/nu_f.py`'s `ExactLP`,
  `max_coverage`, or `lower_from_cert`); `TBD`
- verdict: if the certified packing mass reaches `>= 21`, the pure point-cover method is dead at
  `m = 5` too (same conclusion as `m = 4`, `task D`'s question answered: the fractional gap looks
  universal across the family, not a small-`m` wall artefact); if the heuristic cover goes
  `< 21` *and* a zero-margin exact check certifies it, that is a new theorem `s(21) = 5` -- but
  see the compute-budget note in `tasks/rung2-family/README.md`: a full exact certification at
  `m = 5` was not attempted (bracket only, per the brief).

## 4. What is certified, heuristic, or a bound, explicitly

| statement | status |
|---|---|
| Rung 1 (Friedman's 14 points certify `s(15) = 4`) | **certified**, exact `Fraction` arithmetic, re-verified byte-identical after vectorisation |
| `zeromargin.py`'s vectorised pre-filter is sound (a superset of the exact test) | **certified** by construction (same lenient-`1e-9` inequalities as the pre-vectorisation code, now over `numpy` arrays) and by the byte-identical regression |
| `runs/closed4_best.txt` costs 12.4175 (LP) / 12.509 (honest, stress-checked) | **heuristic** (sampled LP + a dense but finite scan) |
| `nu_f(4) >= 10.68` | **bound** (rigorous, float-checked with `1e-9` margins, per `CLOSED4.md` §7's own standard -- not exact-rational; `search/dual_exact.py`'s fully exact-rational method was not re-run for `m = 4` here, and was not attempted for `m = 5` given the time budget) |
| rung 2 (`W < 13`, exactly zero-margin certified at `m = 4`) | **not established either way** -- see §2 |
| `m = 5` heuristic cover LP value | **heuristic**, still rising when last checked (see §3) |
| `m = 5` certified packing mass | **bound** (rigorous, float-checked, same standard as `nu_f(4) >= 10.68`), pending §3's placeholder |

## 5. Files

| file | what |
|---|---|
| `search/zeromargin.py` | exact checker (vectorised), `friedman14`/`cert`/`pose` modes, `--oracle` |
| `search/rung2_close.py` | rung-2 separation loop (stress -> hard rows -> resolve) over a fixed column set |
| `runs/closed4_rung2c.log` / `.json` / `_best.txt` | the `m = 4` closing-loop run |
| `runs/closed4_s5.log` / `.json` / `_best.txt` | the `m = 5` heuristic cover LP run |
| `runs/closed4_s5lower*` | the `m = 5` certified packing-mass run (§3), once run |
| `certificates/...` | any certificate shipped from a converged run, with the exact `zeromargin.py cert ...` command that verifies it (see the relevant section above) |
