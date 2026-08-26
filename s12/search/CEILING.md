# The ceiling of the method: bracketing `ν_f(s)` (TODO item B1)

Code: `search/nu_f.py` (exact-placement LP runs, post-hoc rigorous lower bound, verifier-based
upper bound), `search/ceiling_launch.sh` (the batch), `search/ceiling_lower.sh` (post-hoc lower
bounds), `search/ceiling_collect.py` and `search/ceiling_lower_table.py` (the tables).  Everything
below is one ~35-minute batch on 32 cores plus ~25 minutes of post-processing; raw logs and
certificates land in the gitignored `runs/`.  **Coarse but honest**: the rigorous bracket is wide,
and the paragraph "what is proved" says exactly why.

## The two sides, and which way each number bounds

Container `C = [0,s]²`, `Π` = placements of a **closed** unit square inside `C` (the closed
convention is the one `certificates/FORMAT.md` and `verify/` use).

```
COVER(s) = min  Σ_p w_p          s.t.  w(Q) ≥ 1  for every Q ∈ Π,                w ≥ 0
ν_f(s)   = sup  μ(Π)             s.t.  μ({Q : p ∈ Q}) ≤ 1  for every point p ∈ C,  μ ≥ 0
```

Weak duality is elementary: `μ(Π) ≤ ∫ w(Q) dμ = Σ_p w_p μ(Q ∋ p) ≤ Σ_p w_p`, so
`ν_f(s) ≤ COVER(s)`.  Both functions are non-decreasing in `s`.  The method proves `s(12) ≥ s`
exactly when it exhibits a cover of weight `< 12` at `s`, so

* `U(s) < 12` (an explicit cover) ⟹ the method works at `s` — that cover *is* a certificate;
* `L(s) ≥ 12` (an explicit packing) ⟹ **no** cover of weight `< 12` exists at `s` or at any
  `s' ≥ s`, whatever LP, cell decomposition or verifier one uses.  This direction does not need
  strong duality.

Strong duality (`ν_f = COVER`) is not used anywhere; `[L, U]` brackets both quantities.  Define
`s*` = sup{ s : COVER(s) < 12 } — the ceiling of the method.

**`U(s')` — rigorous upper bound on `COVER(s')`** (exact arithmetic end to end).  LP weights are
exported as a certificate (integer coordinates, weights rounded *up*), rescaled to container `s'`
(points scale with the container, `D → D/λ`, exactly as `search/scale_to_critical.py`), and handed
to the **exact `i128` verifier** at `N = 6000`, which returns the exact minimum `m(s')` of the
captured weight over *all* closed unit squares in `[0,s']²` (its angle net checks a slightly shrunk
square, so `m` is a lower bound on the true minimum — the safe direction).  Then `w/m(s')` is a
cover and `U(s') = Σw / m(s')`.  Every certificate from every run (plus the shipped one) is
evaluated at every `s'` of the grid and the minimum is taken.

**`L(s)` — rigorous lower bound on `ν_f(s)`** (floating point, but a direct check of an explicit
object).  A finitely supported measure `μ` on exact admissible placements (D4-symmetrised) is
built by the post-hoc step `nu_f.py lower`: take a certificate's point set, solve the cover LP over
those points against a broad row set (lattice placements with `w(Q) < 1.05` at `(η, δθ) = (0.01,
0.01)`, plus the verifier's exact worst placements), keep the dual's support plus the 3000 tightest
placements, and re-solve the *packing* LP on that support against every grid point of pitch 0.01
in the D4 fundamental domain (closed squares).  Its coverage function `cov(p) = Σ_r μ_r [p ∈ Q_r]`
is then bounded at **every** `p ∈ C` by adaptive subdivision: a cell of half-side `h` centred at
`g` is dominated by the square dilated by `h(|cos θ_r|+|sin θ_r|)` in the square's own frame (for
`p` in the cell, `|R_θ(p) − R_θ(g)|_∞ ≤ h(|cos θ|+|sin θ|)`), so the dilated coverage at `g`
bounds `cov` on the whole cell; over-covered cells are split to depth 9 (half-side `0.01/2⁹`).
With `M` the resulting global bound, `μ/M` is feasible and `L = Σμ / M ≤ ν_f(s)`.  The worst
leaves are fed back as extra constraints for up to 7 rounds; every round gives a valid `L`, the
best is kept.  Margins: `1e-9` on every containment test, placements checked inside `C` with margin
`1e-12`; all placements are D4-symmetrised explicitly (the LP's dual constraint only bounds the
orbit-*averaged* coverage — a first version of this code forgot that and reported `M = 8`).

**The LPs in the middle are heuristic.**  Two families of runs produced the point sets:

* *exact mode* (`nu_f.py run`, 9 runs, 23 min each, 3 processes each): points anywhere on a
  `1e-4` lattice (D4 orbits), rows are exact closed-square placements, no erosion or dilation;
  cutting planes from an `(η, δθ) = (0.01, 0.01)` lattice plus the verifier's `topk = 6` worst
  placements per angle bin at `N = 2000` every iteration; column generation at the leaves where
  the symmetrised dual coverage exceeds 1.
* *cell mode* (`lp_search.py`, unchanged, `fine = η = δθ = 0.005`, 4 runs, 33 min each): the
  rigorous eroded-cell LP that produced the shipped certificate.  All four ended with an
  exact-verified cover at their own `s` (in-loop verifier, `N = 2000`, min `1.000002`) — LP totals
  `12.131 @ 3.93` (not yet viol-free, min 0.977), `12.362 @ 3.94`, `12.264 @ 3.95`, `12.593 @ 3.96`.

Neither LP value is a bound on anything by itself (a partial cut set relaxes the cover, a partial
column set restricts it); only `L` and `U` are.

## Results

Rigorous bracket `L(s) ≤ ν_f(s) ≤ COVER(s) ≤ U(s)`, with the heuristic estimates alongside.
Verifier at `N = 6000`; `s'` is the exact rational container the verifier was given (differs from
the nominal value by `< 1e-6`).

| s | **L(s)** rigorous lower bound on ν_f | **U(s)** rigorous upper bound on COVER | U from | packing mass on 0.01-grid (heuristic ν_f) | cover-LP over the cert's own points (heuristic) |
|---|---|---|---|---|---|
| 3.930000 | 9.691 | (see 3.9318) | | 10.49 | 11.68 |
| **3.931795** (= 3920/997, the bound) | **10.611** | **11.9330** | shipped certificate, min 1.0000023 | 11.51 | 11.74–11.82 |
| 3.940000 | 9.492 | **12.0164** | shipped certificate scaled, min 0.993059 | 10.90 | 11.89 |
| 3.945000 | 10.070 | 12.1149 | shipped certificate scaled, min 0.984985 | 11.45 | 11.94 |
| 3.950000 | 10.122 | 12.2219 | exact run @3.96 scaled, min 0.977235 | 10.77 | 11.90 |
| 3.952500 | – | 12.2556 | exact run @3.96 scaled | | |
| 3.955000 | **10.822** | 12.2642 | cell run @3.95 scaled, min 1.000002 | 11.62 | 11.91 |
| 3.957500 | – | 12.2642 | cell run @3.95 scaled | | |
| 3.960000 | 9.941 | 12.2642 | cell run @3.95 scaled, min 1.000002 | 10.70 | 11.96 |
| 3.965000 | – | 12.2642 | cell run @3.95 scaled, min 1.000002 | | |
| 3.970000 | 10.153 | 12.5292 | cell run @3.95 scaled, min 0.978853 | 11.49 | 12.19 |
| 3.980000 | 9.759 | 13.1786 | cell run @3.96 scaled, min 0.955588 | 10.00 | 12.48 |
| 3.990000 | – | 17.17 | (uninformative: the s=4 run's certificate) | | |
| 4.000000 | 10.235 | 18.56 | (uninformative) | 11.94 | 14.77 |

Discretisation parameters and how the numbers move:

* `U`: the only discretisation is the verifier's angle net.  `N = 6000` vs the in-loop `N = 2000`
  changes `m` in the 4th decimal (`σ`-shrink `≈ 1/N`); the shipped certificate's `m(3.931795)` is
  `1.0000023` at both.  The *quality* of `U` is entirely the LP's: a 33-minute cell run at 3.94
  gives 12.36 where the 2-hour run behind the shipped certificate, scaled, gives 12.016.
* `L`: grid pitch 0.01 for the packing LP, subdivision to half-side `2e-5`.  `M` starts at
  1.3–1.6 with the grid alone and comes down to 1.02–1.14 after 7 rounds of adding worst leaves
  (e.g. 3.955: mass 11.62, `M` 1.30 → 1.07, `L` 8.6 → 10.82); the packing *mass* moves little
  (`≤ 0.1`) as constraints are added, so the remaining gap between `L` and the mass column is
  almost all certification slack (`M > 1`), not a real drop.  Support size matters more than the
  grid: the dual's own support (10–40 placements) gives `L ≈ 9.4` at 3.9318; adding the 3000
  tightest placements gives 10.6.  Point-set size matters too: the final iterates (~10k atoms)
  beat the best-U iterates (~400–1500 atoms) by ~1 unit of `L`.

## What is proved, what is estimated

**Proved (exact arithmetic):** `COVER(3.931795) ≤ 11.9330 < 12` (the shipped certificate), so
`s* ≥ 3920/997`.  The shipped certificate scaled to 3.94 has exact minimum `0.993059`, i.e.
`COVER(3.94) ≤ 12.0164`; scaled to 3.945, `COVER(3.945) ≤ 12.1149`; and covers of weight 12.264
exist up to 3.965.  These upper bounds on `COVER` are upper bounds on `ν_f` too.

**Proved (float checks of explicit packings):** `ν_f(3.931795) ≥ 10.61`, `ν_f(3.955) ≥ 10.82`,
and `ν_f(s) ≥ 9.5` on the whole grid.  None of the packings reaches 12, so **this batch does not
establish a rigorous upper end for `s*` below 4.**  The trivial one is `s* ≤ 4`: for every `ε > 0`
sixteen closed unit squares fit disjointly in `[0, 4+ε]²`, so `COVER(4+ε) ≥ 16`.  Rigorous
interval: **`s* ∈ [3.931795, 4]`**.

**Estimated:** the heuristic columns put `ν_f` at roughly 11.5–11.9 across 3.93–3.96 and above 12
from about 3.965–3.97 (cover-LP over the run's own points: 11.96 @ 3.96, 12.19 @ 3.97, 12.48 @
3.98), while the best actual covers cross 12 between 3.9318 and 3.94 (12.016 @ 3.94 from a
certificate optimised at 3.92 and merely scaled).  Taken together the crossing is most likely in
**`3.94 ≲ s* ≲ 3.97`**, consistent with the README's "3.95–3.96" — but that is an estimate: the
LP values in the heuristic columns are neither upper nor lower bounds, and the cell LP carries an
erosion handicap of `η/2 + 0.354·δθ ≈ 0.4%` of the side (≈ 0.017 in `s`) that biases its totals up.

**Slack of the shipped certificate.**  `ν_f(3920/997) ∈ [10.61, 11.933]`.  Two heuristic
figures narrow this: the cover LP over the certificate's own 788 points with 17.6k exact rows is
11.817 (so the shipped weights are within 0.12 of optimal *for that point set*, the price of the
eroded-cell LP plus rounding), and the best grid-feasible packing found is 11.51.  So the shipped
certificate is probably within 0.1–0.4 of `ν_f` at its own `s`, and there is little to gain at
3.9318 itself; gains come from moving `s` up.

**`ν_f(4) = 16`** holds in the *open* (interior-disjoint) convention: the 4×4 grid is a packing of
mass 16, and 16 is an upper bound by area.  In the **closed** convention used by the certificates
and the verifier it does **not**: closed grid squares share edges, where the coverage is 2, and in
one dimension the same effect gives `ν_f^{closed}(4) = 3`, not 4 (the points `{1,2,3}` cover every
closed unit interval in `[0,4]`).  In general `ν_f^{closed}(s) = lim_{t↑s} ν_f(t)` (shrink a packing
of side-`(1+ε)` squares), so the closed-convention value at `s = 4` is the *left limit* — the
number the crossing question actually cares about — and the jump to 16 happens just above 4.  This
batch measured `ν_f^{closed}(4) ≥ 10.23` rigorously and ~11.9–14.8 heuristically (the exact-mode LP
at `s = 4` stalled at a degenerate cover: its optimum wants points *exactly* on the grid lines, which
leaf-centred column generation never proposes).  The README sentence "jumps to 16 at s = 4" is
correct for the open convention and should say so.

## Caveats

* `U` is exact modulo the verifier (already audited in `VERIFICATION.md`); the LP only chose `w`.
* `L` uses floats with explicit margins; the subdivision bound is conservative (dilation, budget-
  limited depth), so `M` overestimates and `L` underestimates — the direction that keeps `L` valid.
* `L` is loose by roughly 1–1.5 units everywhere, for two reasons that are compute, not
  mathematics: the packing support is a few thousand placements taken from one LP dual, and
  `M` is 2–14% above 1 when the round budget runs out.  Closing the rigorous interval from above
  (showing `L(s) ≥ 12` somewhere below 4) needs a much richer support — plausibly a direct packing
  column generation at pitch ≈ 0.002 — and is the natural next step; the crude mass column suggests
  it would land near 3.97.
* All exact-mode runs were 23 minutes and never converged (`M` in the in-loop check stayed at
  1.3–1.9); their in-loop `L` values (8–10) are superseded by the post-hoc step and not tabulated.
* The 3.99 and 4.0 rows of `U` come from the only certificate that scales there (the stalled
  `s = 4` run) and mean nothing beyond "a cover of that weight exists".
