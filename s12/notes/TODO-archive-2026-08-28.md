# Where this could go next

## Status (2026-08-28)

B3, first level, done and measured: **corner-occupancy branch certificates** (`search/BRANCH.md`).
Single and per-box multipliers, exact verification in Rust and Python, Lean reduction for one
region and for families, 68 rejection tests, LP driver with checkpoints.  At `s = 3.98` the
leaves `k = 0, 1, 2` are certified (`certificates/branch/`), the mixed leaf closes per box
(11.97), and the all-corners leaf `k = 4` is **exactly 12** — a dyadic grid cover whose dual is
a mass-12 fractional packing with a full unit per corner.  Consequences for the plan:

* **Established**: `COVER(3.99) = 12.2009` (the 12.008 of `DUAL_EXACT.md` was only a lower
  bound); the pure ceiling is ≈ 3.975–3.98; corner branching gains 0.4–1.5 on the leaves an
  integral packing does not resemble and nothing on the one it does.
* **Next (extension)**: level 2 = wall-strip occupancy (eight boxes; the `k = 4` grid cover and
  its dual measure locate the eight units of mass).  Needs the LP loop redesigned first: a leaf
  is 3–8 h now (IPM cost superlinear in rows, verifier quadratic in atoms; see BRANCH.md
  "Compute"); level 2 has dozens of leaves.  Candidates: solve the dual packing LP directly with
  column generation, or row-and-column generation with a persistent basis (highspy).
* **Not worth it**: squeezing the marginal bound `s(12) ≥ ~3.975` from the corner level.

## Status (2026-08-26)

Five things were measured today; they are separated below into what is **established** and
what is only **suggested**.  Write-ups: `search/DUAL.md`, `search/DUAL_EXACT.md`,
`search/CLOSED4.md`, `search/TIGHTSET.md`, `search/N11.md`, `notes/proof-anatomy.md`.

**Established (rigorous unless marked heuristic):**

* **The ceiling of the pure LP method is pinned: `s* ∈ [3.968616, 3.99)`.**  An explicit
  fractional packing at `t = 399/100` of mass `12.00823` with coverage `< 1`, certified in
  exact rational arithmetic, so no cover of weight `< 12` exists at any `t ≥ 3.99`.  B1 is done.
  Heuristic: the crossing is probably just below 3.99 (`L(3.98) ≥ 11.918`, `L(3.97) ≥ 11.807`,
  float-certified but not converged).
* **Hence a closed-semantics cover of `[0,4]²` costs `≥ 12.008`.**  Best explicit one found:
  `12.51` (heuristic, stress-tested; the sampled LP sits at 12.3–12.4); 91 % of its weight lies
  exactly on the grid lines `x, y ∈ {1,2,3}`, none on the walls — the optimum is a segment
  measure, as `CEILING.md` had guessed.
* **`s(11) ≥ 3040/797 = 3.8143`** (was Stromquist's 3.7889), exact certificate, shipped.  The
  `n = 11` LP crosses 11 at `≈ 3.815` (heuristic bracket of width 0.001), i.e. `0.06` below
  Trump's `3.877083`: the pure method has a required gap for `n = 11` as well, twice the one
  for `n = 12`.  The `n = 11` half of B4 is done; its old claim ("gains nothing") was wrong.
* **Structure of the extremal fractional packing at 3.99** (float-certified measure): half the
  mass is a rigid frame — four corner squares at mass `0.85` each, eight axis-aligned wall
  squares at `0.22–0.27` — and half is a genuinely fractional many-angle interior (tilted
  poses leaning on the walls, a ring of 45° diamonds).  **Near-tight set of the shipped
  certificate**: at the critical size one `C4` orbit of corner poses at *all* angles; at
  `E = 0.134` a frame plus a central pinwheel.  Rescaling a critical certificate never gives
  `E < 0.134`, because `m(t)` jumps discontinuously above the critical side (atoms sit on the
  tight squares' boundaries); a small `E` at `t` needs a cover re-optimised *at* `t`.
* **Literature anatomy** (Bentz 2010/2016, Nagamochi 2005, Stromquist 1984/2003): every
  published `s(m²−3)` proof works in closed semantics (dilated squares, boundary points count);
  deficit of the best pure point set against `n − 1` is 3 for `n = 12` — the worst case in the
  table (13: 2, 22/33: 1, 46: 0, 45: 1); every leaf of Bentz's 13-proof is a cost-12
  configuration; and "≤ 3 unit squares in a strip of height < 2 and length < 4" is **false**
  (Stromquist 1984-III: four fit in `1.9 × 3.9475`), so Nagamochi's rectangle theorem (which
  needs both sides ≥ 2) says nothing about the wall strips of our problem.
* **Correction to A1 below.**  Erosion-based verification (our σ-shrink, cell cores, Mira's
  subdivision) cannot certify *any* closed-semantics cover of `[0,4]²` with `W < 16`, whatever
  the erosion: it certifies eroded squares, and 16 eroded unit squares fit disjointly.  Covers of
  cost 14 do exist there (Friedman), so this is a statement about verifiers, not covers: any
  argument at the limit `s = 4` needs exact verification at the zero-margin poses (squares on
  a wall or a grid line), which is what Bentz's non-avoidance lemmas are.  A1 is therefore not
  "worth ≤ 0.001" — it is mandatory for a limit argument and worth nothing for the bound.

**Suggested (not established; ordered by my estimate of value):**

* **B3 design.**  Certificate = points **plus segment densities on the `#` lines** (for
  rational rotations the length a square captures on an axis-parallel line is rational, so
  segments are exactly verifiable with the existing arithmetic) **plus pose-region capacities**
  with the Lagrangian reduction `n ≤ Σ w + Σ_R k_R y_R` (Lean generalises trivially).  First
  branch: how many corners hold a square (the 3.99 measure wants `0.85` in each of four; an
  integral packing has 0 or 1).  Leaves: compact capacity statements certified by the eroded
  LP; wall strictness lemmas ("each box takes `> 1` of the line") supplied analytically, as in
  Bentz/Nagamochi.  The excess budget at the limit is `0.008–0.5`, not 3.  Size of the case
  tree: unknown.  Nothing here has been tried.
* **Measure `E(t)` properly**: re-optimise a cover *at* `t = 3.97–3.99` with `tighten.py`
  (verifier as oracle) instead of rescaling; the near-tight set of such a cover is the honest
  input to a case split.
* **`n = 45`**: run the closed-semantics cover LP at `[0,7]²` (generalise `closed4.py`).
  Literature deficit is 1; a weighted cover of cost `< 45` would prove `s(45) = 7` with no
  case analysis, given exact zero-margin verification.  Cheapest high-upside experiment.
* **`n = 11`**: leave it; gap 0.06 and the target packing is tilted and irrational.
* **Outreach**: the DS7 email should now carry both bounds (`n = 11` and `n = 12`).

## Status (2026-08-25)

Bound now `s(12) >= 15680/3951 = 3.968616` (was 3.931795), from re-optimising weights with
the exact verifier as separation oracle plus column generation (`search/TIGHTEN.md`); the
certificate at the old bound shrank from 788 to 224 points.  Column generation at 3.9696 no
longer gets below 12, so **this method is exhausted to within ~0.001** — B2 is effectively done
and A1 (removing the σ-shrink) is now worth at most that much.  Everything past here is B3.
Also new: a uniform certificate `s(12) >= 35/9` (81 points, 7 per square,
`search/uniform/UNIFORM.md`), an exhaustive Python checker, 42 rejection tests that found and
fixed two verifier soundness bugs, reproducible LP search, `points.json` companions, and an
adversarial audit of Mira's checkers (`outreach/mira-audit/`).

## Status (2026-08-23)

Track 0 — **ship** — is done: the repo is published, CI re-verifies every certificate and
runs the rejection tests on every push, certificates are pinned by SHA-256, and the
generation path (LP → scale-to-critical) is documented and scripted, so the bound can be
*reproduced* and not merely *checked*.

Three defects were found and fixed while doing it, all in the credibility-critical docs:

* `VERIFICATION.md` claimed the main certificate verifies at `N = 2000`.  It does not — it
  is rejected up to `N ≤ 4000` (min covered weight `0.998704` at `θ ≈ 20.96°`) and first
  passes at `N = 6000`.  The `1.000002` figure quoted against `N = 2000` was the `N = 6000`
  result.  The bound is unaffected; the log was wrong.  This is the σ-shrink handicap
  biting on our own certificate, and is now recorded as evidence for A1.
* `certificates/FORMAT.md` claimed the 56-point set is critical at `λ = 400/397`
  (container `1520/397`).  It is not: there an axis-aligned square captures only `4/5`,
  stable from `N = 2000` to `N = 200000`.  True critical size is `760/199 = 3.819095`.
* `VERIFICATION.md` pointed at `runs/…` paths that are gitignored and absent from a clone.

Remaining work is the two threads below.

Two independent threads.  They share nothing but the subject, so they can be picked up in either
order, or by different sessions.

**Thread 1 — sharpen the bound by learning from the other machinery.**  Items A1–A2, B1–B3.
The core of it is one substitution: our uniform angle net costs us a `σ`-shrink that is a handicap
on the *certificate*, not just the check.  Mira's adaptive pose-space subdivision has no shrink and
certifies to `3×10⁻⁷` where we need `4×10⁻³`.  Swap that in, add a weighted analogue of the
triangle-piercing lemma, and the remaining ~0.02 up to the `ν_f ≈ 3.95` ceiling should fall out.
Past the ceiling needs B3, which is a genuine research problem.

**Thread 2 — collate, connect, and contribute.**  Items D1–D5.  Note the mechanism differs by
target and "send a PR" is often the *wrong* one:

| target | mechanism | why |
|---|---|---|
| Mira, Fort | GitHub issue | repos are live; issues enabled; cross-check + license ask |
| Burns, Massaccesi | blog comment / Twitter | no repos exist for their certificates |
| Ellsworth (kingbird) | GitHub issue on `Davidebyzero/*` | no repo behind the site, no email anywhere |
| Friedman (DS7 Table 2) | **email**, not PR | no outside PR has ever been merged; he edits by hand |
| `formal-conjectures` | PR | genuinely PR-friendly, and we have Lean to offer |
| a lower-bound registry (D5) | build it | nothing like it exists |

One caution for thread 2.  Everything in this space right now is unrefereed and computer-produced,
including ours, and both n=17 blog authors flag their own numbers "(?)".  Adding a fifth
unrecorded claim is low value; being the one that brings independent cross-checking, rejection
tests, a formal reduction, and a place for results to live is high value.  Lead with the
verification story, not the number.

---

Ordered by value, with an honest note on what each is worth.

## A. Adopt from the Mira / Fort lineage

**A1. Replace the uniform angle net with exact adaptive subdivision of pose space.**  *Corrected 2026-08-26 — see the status above: no erosion-based verifier works at the limit; what is needed is exact verification at zero-margin poses.*
*This is the single highest-value item.*  Our verifier enumerates angles `θ_k = 2·arctan(k/N)`
and checks a square shrunk to side `σ_k = 1/(cos δ + sin δ)`.  That shrink is not just a
verification detail — it is a **handicap on the certificate**, because the LP must produce weights
valid for slightly-shrunken squares.  Mira's architecture instead subdivides the 3-D box
`[½, L−½]² × [−1, 1]` in `(x, y, t = tan θ)` recursively with exact integer arithmetic and no
shrink at all, so it can certify certificates that are critical to `3×10⁻⁷`; ours needs
`~4×10⁻³` of slack.  Adopting it would remove the handicap, make the scale-to-criticality step
unnecessary, and is plausibly worth most of the remaining `0.02` up to the method's ceiling.

**A2. Add a triangle-witness primitive.**  Their strict triangle-piercing lemma — *if A, B, C have
pairwise distances `< 1`, every unit square centred strictly inside triangle ABC contains one of
them* — is orientation-independent, so it kills whole `t`-columns at once (483,875 triangle leaves
do the work of far more point leaves).  There is a natural weighted analogue: a triangle whose
three vertices carry weights summing to `≥ 1` certifies its whole region for all orientations.
Worth adding both as a verification shortcut **and** as a constraint type in the LP, where it may
produce genuinely better certificates rather than merely cheaper checks.

**A3. ~~Make the independent Python checker exhaustive.~~  Done — `xcheck.py --all`; see `VERIFICATION.md` 4b.**  `xcheck.py` currently samples angle bins.
Mira ships three checkers that each check everything (fast C++, bigint C++, pure Python).  Ours
should too, even if slow.

**A4. Deterministic regeneration, artifact hashes, CI.**  ~~Fort's repo runs
`verify_archived.sh` in GitHub Actions on every push; Mira publishes SHA-256 of both
compressed and raw certificates.~~  **Done** — `.github/workflows/verify.yml` runs the full
`verify.sh` plus a hash check on every push and monthly; `certificates/SHA256SUMS` pins both
artifacts; `search/scale_to_critical.py` reproduces the scaling step.  Closed: `lp_search.py --iters --seed` is
bit-reproducible on one machine and `--dump-lp`/`--resolve` archive and replay the exact LP
(`search/REPRODUCIBILITY.md`); cross-machine identity is not claimed.

**A5. ~~Publish `points.json` alongside the plain-text format.~~  Done — `certificates/*.json`, `search/export_points.py --roundtrip`.**  Their self-describing schema is
better than our bare integer file; add a `weights` field and a `weight_denominator`, and our data
becomes readable by anyone in either family.

## B. Push the bound

**B1. Pin down the ceiling.**  **Done 2026-08-26: `s* ∈ [3.968616, 3.99)`, exact (`search/DUAL_EXACT.md`).**  *Was: partly done* (`search/CEILING.md`, `search/nu_f.py`): rigorous upper bounds `U(s)` on `ν_f` from scaled exact certificates (12.016 at 3.94, 12.26 at 3.955); rigorous lower bounds from explicit packings are too loose (~10–10.8) to show `ν_f ≥ 12` anywhere below 4, so the proved bracket is still the trivial `[3920/997, 4]`.  Heuristic crossing ≈ 3.94–3.97.  Closing the lower side needs direct packing column generation at finer pitch.  Also found: `ν_f(4) = 16` is an open-convention statement; in the closed convention used here `ν_f(4)` is the left limit.  By LP duality the method dies where the fractional packing number
`ν_f(s)` reaches 12; we estimate `s ≈ 3.95–3.96` from partial runs.  Worth computing properly — it
tells everyone, including us, exactly where this technology stops.

**B2. ~~Reach the ceiling.~~  Done to within ~0.001 (`search/TIGHTEN.md`).**  With A1 in place, LP runs at `s = 3.94–3.95` with finer cells should
close most of the remaining `0.02`.  Compute, not mathematics.

**B3. Beat the ceiling — certificate + case analysis.**  The only route to `s(12) = 4` by this
technology, and the genuinely novel research direction.  When the total weight `W` sits just above
12, the excess `W − 12` bounds the *entire* packing's slack: since each of the 12 squares captures
`≥ 1`, at most `(W − 12)/δ` of them can capture `≥ 1 + δ`.  So nearly every square must sit in a
near-binding placement.  Then ask a much smaller question: **how many interior-disjoint unit
squares can simultaneously occupy the low-coverage region?**  If that number is below 12, done.
This is Bentz's layered case analysis, mechanised, with the LP supplying the case split.

**B4. Other open `n`.**  *`n = 11` done 2026-08-26: `s(11) ≥ 3040/797`, and the old claim below was wrong.*  `n = 11` probably gains nothing (`ν_f` reaches 11 well below the current
3.788854 bound), but the `k² − 4` family — `n = 21, 32, 45` — has the same "just below a perfect
square" structure that makes `n = 12` work, and nobody has run it.

## C. Verification hardening

**C1. Extend the Lean development beyond the reduction.**  The 81-point uniform certificate (`s(12) >= 35/9`, every unit square contains 7 of 81 points) is now the natural `native_decide` target: no weights, one sentence.  Lean currently proves *unavoidable set
of weight `W` ⟹ at most `W` squares*, which is the mathematical core but not the computation.
The angle-net inclusion lemma (a unit square contains the concentric `σ`-square at a nearby angle)
is ordinary geometry and very formalisable.  The exhaustive sweep is the hard part — plausibly
`native_decide` on the 56-point certificate, almost certainly not on the 788-point one.

**C2. ~~More rejection tests.~~  Done — 42 checks; they found two soundness bugs (negative weights unchecked; centre box too small above 45°) and several panics, all fixed and logged in `VERIFICATION.md`.**  Was: seven.  Add: non-symmetric sets (exercise the `[0,90°)` path),
points outside the container, zero-weight and negative-weight inputs, malformed headers,
denominators that do not divide the container side.

## D. Communication

**D1. Issue on `Mira-acc/17squares` and `stanislavfort/17squares`** — offer the cross-check we ran
(independent implementation, different architecture, exact rational confirmation of their point
set on the axis-aligned slice, plus the `3×10⁻⁷` criticality measurement), and tell each family
the other exists.  Mira's paper cites neither Burns nor Massaccesi; neither blog cites Mira or
Fort.  They are working in parallel and unaware.

  While there: **ask Mira to add a LICENSE file.**  The repo has none.  That does not affect
  citing the theorem or reusing the point data (facts, not expression), but it does leave the
  status of the checker source ambiguous for anyone wanting to adapt it.  A one-line request,
  trivial for them to act on, useful to everyone downstream.

**D2. Tell Burns and Massaccesi** — we are downstream of their idea, and `n = 12` is a case they
did not try.

**D3. DS7 Table 2** is the only place lower bounds for `s(n)` are tabulated anywhere, and it has
not been updated since 2009; it still lists Green's 4.4452 for `n = 17`, so *none* of the four
2026 results are recorded.  Route: `github.com/erich-friedman/erich-friedman.github.io` (issues
open, outside PRs merged) or Friedman directly.

**D4. Ellsworth / kingbird — corrected.**  The main squares-in-squares table is upper-bounds-only,
but he *also* hosts an edited copy of the DS7 survey at
https://kingbird.myphotos.cc/packing/squares.html , and **that copy does carry Table 2, the
lower-bound table** — currently still showing Green's `(40√2+19)/17 ≈ 4.4452` for `n = 17` and
Stromquist's `3.7888` for `n = 11`, newest reference 2023.  So there *is* a home for a lower bound
on that domain; it is simply stale.  Practical detail: the page has no email and no submission
policy, and there is no repo behind it — the site is static files under a photo-gallery host, and
the per-packing SVGs (`square-<n>.svg`) *are* the dataset, each carrying its exact constants as
XML entities plus a full attribution chain in XML comments.  The only contact route is a GitHub
issue on one of https://github.com/Davidebyzero 's repos (his `erich-friedman.github.io` fork is
the natural one).  The page also has an Edit Mode ([E] to edit, [S] to download an SVG in the
site's own format) — de facto the submission format for *packings*.

**D4b. Friedman's repo is the wrong mechanism; email is the right one.**  On
`erich-friedman.github.io`: issues are open but unanswered, and **no outside PR has ever been
merged** — Ellsworth's PR #5 was closed unmerged yet its content appeared on the site anyway, so
Erich edits by hand from email.  Route: `[email]`, mentioning the file path
`papers/squares/squares.html` (Table 2).  His written submission guidelines cover packings, not
lower-bound certificates, so this would be a slightly novel request.

**D4c. Venues that are *not* available.**  OEIS has no sequence for `s(n)` (it is irrational and
mostly unknown), so there is nothing to submit against.  Wikipedia cites DS7 and links kingbird,
and gives no `n = 17` lower bound at all — but a self-published certificate is WP:SPS and would be
reverted; that needs arXiv or the EJC survey first.  MathWorld carries no lower bounds.

**D4d. One venue that *is* contribution-friendly and fits what we already have:**
`google-deepmind/formal-conjectures` — a Lean formalisation repo that accepts PRs, and which
already has issue #646 "Seventeen square packing problem".  We have a Lean development; a
formal statement of the `s(12)` problem plus our proved reduction lemma would be in scope.

**D5. Collate best-known lower bounds somewhere, with certificates attached.**  *Lower priority,
but a real gap and a service to the field.*  Right now:

* Ellsworth's tables record best-known **packings** and mark which are **proven optimal** — so a
  proof of `s(12) = 4` would flip that entry and is submission-worthy there, but an intermediate
  bound like `3.931795` has no home in them at all;
* DS7 Table 2 is the only lower-bound table anywhere, and has not moved since 2009;
* the four 2026 results (Fort, Mira, Burns, Massaccesi) are recorded nowhere, and the two families
  producing them do not cite each other.

A plain, well-cited table of best-known lower bounds — each row carrying a machine-checkable
certificate and a verifier that reproduces it — would fix all three at once.  It is mostly
curation rather than research, which is exactly why nobody has done it.  Natural shape: a repo
with one directory per `n`, the `points.json`-style certificate, and CI that re-verifies every
entry on push, so the table cannot rot the way DS7 did.
