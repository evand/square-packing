# helly-allmeet: can a sound, positive-volume clique rule credit the fan cliques? (2026-09-12)

Task: `tasks/helly-allmeet/README.md`.  Code: `search/allmeet.py` (nothing in `verify/`, `lean/`,
`certificates/` is touched).  Runs in this worktree's `runs/am_*` (gitignored; every number is
quoted here).  Read against `search/HONEST.md` (whose dual, `Theta` and `meet`/`core` numbers this
reproduces to nine places), `search/CLIQUELEVER.md` §0/§6, `search/ZEROMARGIN.md` §1,
`search/RUNG2.md` §6, `notes/review-2026-09-12.md`.

## 0. Verdict, up front

**Yes, every row is soundly creditable with positive volume — and no, it does not recover the fan
mass.**  On **both** instances, **every** dual-carrying clique row admits a **positive-volume box
of poses `B_K`, exactly certified**, such that `K ∪ B_K` is pairwise closed-meeting, so
"credit `z_K` iff `S ∈ K ∪ B_K`" is a *sound* clique rule with positive-volume membership on
`100 %` of the clique dual (`317/317` rows and `7.352132/7.352132` on `SUPEPGF`; `638/638` and
`7.610316/7.610316` on `E2Pg`).  The boxes are not slivers — median `vol(B_K) = 0.16` and `0.13`
against an admissible pose volume of `7.50`, i.e. `2 %` of *all* poses each, largest `0.37`.  But
the honest cost of the dual under that rule is

| dual | `Theta` | rule `meet` (generous, UNSOUND) | rule `core` (sound, credits nobody) | **rule `K ∪ B_K` (sound, positive volume)** |
|---|---|---|---|---|
| `SUPEPGF` (956 poses) | `11.800418089` | `20.162` | `inf` | **`1543.26`** (min capture `0.007646`) |
| `E2Pg` (14,611 poses) | `11.864926892` | `18.109` | `162.47` | **`145.36`** (min capture `0.081627`) |

— against the `12` that would put V1 back on the path.  Re-placing every box around the worst pose
(the placement most favourable to the cover at that pose) and descending again does not rescue it:
the best I could reach that way is `393.71` on `SUPEPGF` and `136.10` on `E2Pg`.  So the answer to
`TODO.md` critical-path item 2 is the second branch the brief named: *yes, such a `B` exists for
every row that carries the mass, and the honest cost with the sound rule stays far above 12.*
**The cover route for `n = 12` does not close below 12 on this family, and the reason is not
Helly.**

Three things are worth writing down, in order of how much they change the picture.

1. **`18.1` and `20.2` were never "floors to be improved on" — they are ceilings on what any sound
   rule can do, and the box rule falls an order of magnitude short of them.**  If a crediting rule
   credits every member of `K` (it must: otherwise the dual is not feasible for the loaded columns
   and the self-test of `HONEST.md` §1 fails by construction) and its credited set `C_K` is a
   clique (it must: that is exactly the hypothesis the accounting `Σ_i capture(S_i) ≤ Theta`
   needs), then every `S ∈ C_K` closed-meets every member, i.e.
   `C_K ⊆ A(K) = {S : S meets every member of K}`.  Hence `capture_sound(S) ≤ capture_meet(S)` at
   every pose, `min capture` under any sound rule is at most `HONEST.md`'s `meet` minimum, and
   **`honest ≥ 20.162` (`SUPEPGF`) / `≥ 18.109` (`E2Pg`) for every sound rule whatsoever**.  That
   is a statement about the dual, not about boxes.  §2 measures how much further a box-shaped rule
   loses on top of it: a factor of `77` and `8`.
2. **The rows are not grazing, and they are not Helly.**  `CLIQUELEVER.md` §6 describes the
   violated cliques as a point clique plus a fan of "grazing contacts"; measured, **not one of the
   703–5,253 member pairs on any of `SUPEPGF`'s ten heaviest rows grazes at `1e-6`** — the least
   pairwise meet margin is `+0.00035` to `+0.00080` and the median is `+0.71` to `+0.95` (§3).
   They are robustly pairwise overlapping families with an empty core (deficits `-0.0003` to
   `-0.100`, re-derived here by an independent max-margin LP that converges to the same value with
   and without cutting planes).  Every row on both instances has a pose meeting every member with
   margin `≥ 0.468`; `A(K)` is fat, `vol(A(K)) = 0.72` to `0.91` out of the `7.50` of all
   admissible poses.  What a box costs is **soundness**: `vol(B_K)` is `0.4 %` to `22 %` of
   `vol(A(K))`, and that factor is the
   whole distance from `20.162` to `1543.26`.  (On `E2Pg` there is a second species: rows `5`,
   `220`, `588` have core deficits of only `-0.0003` to `-0.042`, `132`/`15`/`153` genuinely
   grazing pairs, and a fan shape — one member grazing `33` of `400` sampled others.  These are
   near-point-cliques plus a fan, and they are the minority.)
3. **The literal Helly statement of `TODO.md` is false at `t = 4`, and there is an exact
   counterexample.**  Three unit squares, each with an edge on one of three half planes at `120°`
   whose triple intersection is empty, pairwise meet with margin `+0.358`; thickened into three
   positive-volume pose boxes they stay pairwise meeting (every pose of each box against every
   pose of the others) and *no point lies in all three for any choice of poses from the boxes*, by
   an exact Farkas certificate (§4).  So "every positive-measure pairwise-intersecting family of
   closed unit squares in `[0,4]²` has a common point" cannot be the theorem.  The statement the
   data supports is about LP mass, not about the family — §5.

## 1. What is computed, and what is exact

For each clique row `K` with dual `z_K > 0` on a `cliquelever.py` checkpoint (the dual re-derived
exactly as `HONEST.md` §6 does, `honestcost.py dual`):

* **`A(K)` has interior?**  `maxmin_at_theta` maximises, over centres, the least meet margin
  against every member at a fixed angle.  At a fixed angle each margin is
  `1/2 + w_k/2 - |⟨a_k - c, e⟩|` with `w_k` constant, so `margin ≥ t` is a pair of linear
  inequalities and the max-min is a 3-variable LP (solved by cutting planes — `E2Pg`'s rows have
  up to 1,806 members — with the centre confined to the angle's admissible range
  `[w(θ)/2, 4 - w(θ)/2]`).  `best_allmeet_pose` sweeps `θ` and refines.  If the answer is `≤ 0`,
  `A(K)` has empty interior and no positive-volume box exists; if `> 0`, a whole neighbourhood of
  the maximiser lies in `A(K)`.  **Float.**
* **the box.**  `grow_box` seeds the largest balanced box `centre ± s(1,1,1)` that certifies, then
  pushes each of its six faces out independently (six faces, because the maximiser is often
  pressed against the admissibility wall, where no box centred on it is admissible).  The search
  is steered by `FloatBoxTest`, a vectorised double-precision image of the exact test; **the
  verdict is the exact test**, and the box is shrunk until the exact certificate holds.
* **the exact certificate.**  Three functions, all `Fraction`/integer, run on every reported box:

  `box_meets_square(box, member)` — *every* pose in `B` closed-meets the fixed square.  With an
  exact rational **lower** bound `w_lo` on `w = |cos Δ| + |sin Δ|` over the box's `u`-interval
  (`w_lo_exact`: `Δ = θ - θ_k ∈ (-90°, 90°)` because both angles are in `[0, 90°)` — `fold_pose`
  first re-parametrises every member by a quarter turn so that its `u` lies in `[0, 1)`, which
  leaves the square alone as a point set (asserted, per pose, by comparing exact corner sets;
  `84` of `SUPEPGF`'s 956 poses and `1,998` of `E2Pg`'s 14,611 need it, down to `u = -0.414`) —
  so
  `w = cos Δ ± sin Δ = √2 cos(Δ ∓ 45°)` with the argument in `[-45°, 45°]` where `cos` is concave
  — minimum at an endpoint; and `w_lo = 1` exactly when the interval straddles `Δ = 0`), the four
  separating-axis conditions with `K = 1/2 + w_lo/2` are

      |d·(C,S)| ≤ K,  |d·(-S,C)| ≤ K                                  (the member's own normals)
      ±( d_x(1-u²) + 2 d_y u ) - K(1+u²) ≤ 0                          (the pose's normals, after
      ±( -2 d_x u + d_y(1-u²) ) - K(1+u²) ≤ 0                          clearing 1 + u² > 0)

  with `d = (a_x - c_x, a_y - c_y)`.  Every one is **affine in `(c_x, c_y)` and quadratic in `u`**
  — precisely the class of `zeromargin.py`'s CHAIN violation polynomials (`ZEROMARGIN.md`,
  `RUNG2.md` §6) — so its exact maximum over the box is `max over the four centre-rectangle
  corners of the exact quadratic vertex test`.  **This is the corner argument, not a Bernstein
  bound: no interval enclosure, no subdivision, no slack.**  Using a lower bound for `w` makes the
  condition strictly harder, so the test is sound and never optimistic.  It is sound even without
  the fold — `cos Δ ± sin Δ ≤ |cos Δ| + |sin Δ|` for every `Δ` — but it is *loose*: on `E2Pg` row
  31 an unfolded member at `θ_k = -27.5°` against a box at `θ = 89.9°` gives `w_lo = 0.427` where
  the truth is `1.348`, and the test then refuses boxes that plainly meet the member.  Before the
  fold, `4` of `SUPEPGF`'s 317 rows and `12` of `E2Pg`'s 638 got no box at all and the reported
  volumes were smaller; after it, every row on both instances gets one.

  `box_pairwise_meets(box)` — every two poses of `B` meet.  Each margin is at least
  `1 - ‖c' - c‖₂` because `w ≥ 1` always, and the extreme centre pair is a pair of opposite
  rectangle corners, so the criterion is **`(x₁-x₀)² + (y₁-y₀)² ≤ 1`**, exact in rationals.  (The
  `u`-extent is irrelevant: rotating one square only ever raises `w`.)  Sufficient, not necessary.

  `box_admissible(box)` — every pose's closed square lies in `[0,4]²`: `c ∈ [w_hi/2, 4 - w_hi/2]`
  with `w_hi` the exact maximum of `cos θ + sin θ` over the `u`-interval (an endpoint, or the
  rational `1.41421357 > √2` when the interval contains `45°`).

  Validation: on random `(box, square)` pairs the certificate agreed with
  `leaf_ceiling.sq_meets_sq` on **0 violations** out of 181 certified boxes × 29–41 sampled poses
  each (including all eight corners), and of **271** boxes it *rejected*, **0** had all 68 sampled
  poses meeting — so it is sound and, on this sample, not visibly conservative either.
  `box_pairwise_meets`: 331 certified boxes, 0 corner-pair violations.
* **the honest cost** under `credit z_K iff S ∈ K ∪ B_K`, by `honestcost.py`'s own scan machinery
  (lattice, `family_rows.py` families, pattern descent, knife-edge `±ε` refinement).  **Float**,
  exactly as `HONEST.md` §5 says of its own numbers.

`K ∪ B_K` is pairwise closed-meeting by construction — members pairwise meet (the row is a
verified clique of the checkpoint), `B_K` is pairwise meeting, and every pose of `B_K` meets every
member — so the rule **is** a valid weighting, which the `meet` rule is not.

The pose set is reconstructed from `cl_TAG_poses.txt` exactly (`read_exact_poses` replays
`cliquelever.Poses.add`'s boundary snap and admissibility clamp) and asserted to match the float
frame `honestcost.dump_dual` wrote **bit for bit** on both instances, so the clique member indices
really do point at the squares this note tests.

## 2. The numbers

### 2.1 `SUPEPGF` — the certified 956-pose support (LP `11.800418089`)

`Theta` = points `3.276397` + cliques `7.352132` (317 rows, 26–146 members, median 117) + polygons
`1.171889`; `62 %` clique mass; all 317 cores empty.  Reproduces `HONEST.md` §2.2 to nine places.

**`317` of `317` rows carry an exactly certified positive-volume box**, i.e. `7.352132` of
`7.352132` = `100.00 %` of the clique dual; `0` rows have `vol(B_K) = 0`.  Every row's max-min meet
margin is strictly positive (min `+0.4803`, median `+0.5474`, max `+0.7019`), so `A(K)` has
non-empty interior everywhere.  `vol(B_K)`: min `2·10⁻¹¹`, `q25` `0.119`, **median `0.155`**, `q75`
`0.195`, max `0.372`, `z`-weighted mean `0.117`; the admissible pose volume in `(c_x, c_y, u)` is
`7.502158`.

Honest cost under the sound rule (`allmeet.py honest SUPEPGF --refit 2`):

| candidate set | poses | min capture | honest |
|---|---|---|---|
| the 956 poses of `P` (self-test) | 956 | `0.053232` | `221.68` |
| `0.04 / 2.5°` lattice | 169,689 | `0.029972` | `393.71` |
| `family_rows.py`, pitch `0.004` | 61,656 | **`0.007646`** | **`1543.26`** |
| pattern descent + knife edge | 500 | `0.007646` | `1543.26` |
| box family re-placed around the worst pose, round 0 | — | `0.000000` | `inf` |
| box family re-placed, round 1 | — | `0.029972` | `393.71` |

The minimiser is `(1.5, 2.5, 0°)` — an axis-parallel square in the middle of the left half of the
container, the same place `HONEST.md` §2.2's `core` minimiser sits — and it captures `0` of the 67
atoms, `0` of the 317 boxes and `0.007646` from the polygon rows: **nothing at all from the
`62 %` of `Theta` that lives on cliques.**  The self-test is the sharpest way to say it: under this
sound rule `303` of the `956` poses the LP was solved on are captured below `1`, the worst at
`0.053`.

The `--refit` rows are a deliberately optimistic variant: after finding the worst pose, *re-grow
every box that can be re-grown around it* — the placement most favourable to the cover at that one
pose — and descend again.  Any single pose can be credited by every row whose `A(K)` it lies in,
so one round is an upper bound on what re-placing buys *there*; the new minimum simply moves
elsewhere (round 0 lands on a pose with capture `0`).  The best honest cost over all three box
families tried is `393.71`.  Not `12`.

### 2.2 `E2Pg` — the full loaded set (LP `11.864926892`)

`Theta` = points `3.709598` + cliques `7.610316` (638 rows, 793–1,806 members, median 1,152) +
polygons `0.545013`; `64 %` clique mass; all 638 cores empty.  Reproduces `HONEST.md` §2.1 to nine
places.

**`638` of `638` rows carry an exactly certified positive-volume box** — `7.610316` of `7.610316`
= `100.00 %` of the clique dual.  Max-min meet margins: min `+0.4680`, median `+0.5186`, max
`+0.5637` (visibly tighter than `SUPEPGF`'s `+0.70` maxima, as 1,152-member rows should be, but
none degenerate).  `vol(B_K)`: min `2.1·10⁻³`, `q25` `0.091`, **median `0.125`**, `q75` `0.150`,
max `0.239`, `z`-weighted mean `0.115`.

| candidate set | poses | min capture | honest |
|---|---|---|---|
| the 14,611 poses of `P` (self-test) | 14,611 | `0.090098` | `131.69` |
| `0.04 / 2.5°` lattice | 169,689 | `0.139147` | `85.27` |
| `family_rows.py`, pitch `0.004` | 61,656 | `0.108550` | `109.30` |
| pattern descent + knife edge | 500 | **`0.081627`** | **`145.36`** |
| box family re-placed around the worst pose, round 0 (127 boxes re-grown) | — | `0.073029` | `162.47` |
| box family re-placed, round 1 (114 re-grown) | — | `0.087181` | `136.10` |

The minimiser is `(2.499500, 1.499751, 89.9201°)`: `0.040292` from 18 of the 403 atoms, `0.041334`
from the polygon rows, and **`0` from all 638 boxes**.  `3,574` of the `14,611` poses the LP was
solved on are captured below `1`.  Note that the round-0 re-placement lands *exactly* on
`HONEST.md` §2.1's `core` minimum `0.073029` at `(2.50023, 1.50093, 0.141°)` — the box rule and the
core rule are being cornered by the same pose, and the best of the three box families is `136.10`.

## 3. Step 4: the anatomy of the rows that carry the mass

The brief asks what the geometry is on rows whose `A(K)` has empty interior.  **There are none** —
every row on both instances has a robust all-meeting pose.  The interesting anatomy is therefore
the other one: how tight the rows are internally, and whether the "fan" of `CLIQUELEVER.md` §6 is
a grazing structure.

Top ten `SUPEPGF` rows by `z_K`.  "grazing pairs" counts member pairs with `|margin| < 1e-6`.

| rank | row | `z_K` | `\|K\|` | core deficit | least pair margin | pairs `< 1e-3` | grazing (`1e-6`) | max-min margin | `vol(B_K)` | `vol(A(K))` |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 37 | `0.206960` | 38 | `-0.001938` | `+0.000801` | 2 of 703 | **0** | `+0.6300` | `5.55e-02` | `0.784` |
| 1 | 25 | `0.188389` | 94 | `-0.012170` | `+0.000698` | 2 of 4371 | **0** | `+0.5898` | `4.68e-02` | `0.862` |
| 2 | 26 | `0.178252` | 27 | `-0.001908` | `+0.000350` | 8 of 351 | **0** | `+0.6513` | `3.71e-02` | `0.788` |
| 3 | 9 | `0.171868` | 91 | `-0.012710` | `+0.000698` | 3 of 4095 | **0** | `+0.6200` | `1.80e-03` | `0.912` |
| 4 | 22 | `0.156828` | 27 | `-0.000665` | `+0.000350` | 6 of 351 | **0** | `+0.6838` | `7.93e-03` | `0.814` |
| 5 | 164 | `0.153422` | 93 | `-0.042367` | `+0.000696` | 3 of 4278 | **0** | `+0.5617` | `8.46e-02` | `0.824` |
| 6 | 185 | `0.125558` | 103 | `-0.024929` | `+0.000597` | 3 of 5253 | **0** | `+0.5549` | `1.00e-02` | `0.820` |
| 7 | 115 | `0.122038` | 90 | `-0.009153` | `+0.000698` | 3 of 4005 | **0** | `+0.6200` | `1.80e-03` | `0.910` |
| 8 | 16 | `0.117870` | 32 | `-0.000333` | `+0.000800` | 2 of 496 | **0** | `+0.6845` | `5.38e-03` | `0.822` |
| 9 | 58 | `0.117525` | 97 | `-0.030797` | `+0.000597` | 3 of 4656 | **0** | `+0.5549` | `1.00e-02` | `0.813` |

Read the top three across.  **Row 37** (`z = 0.2070`, 38 members, the heaviest row on the
instance, and `HONEST.md` §0's example of an unsound `meet` row): its core is empty by `0.0019`,
its members pairwise overlap with margin at least `+0.0008` and typically `+0.88` (median), and
exactly **two** of its 703 pairs are tighter than `1e-3`.  There is a pose meeting all 38 with
margin `+0.6300`, and the certified box around it has volume `0.0555` — `7 %` of the `0.784` that
the generous rule credits.  (Per-row `vol(B_K)` is what `boxes` happened to grow, not a maximum:
the row-37 box is `0.055` while the median row reaches `0.155`.)  **Row 25** (`z = 0.1884`, 94 members): core deficit `-0.0122`, least
pair margin `+0.0007`, two tight pairs of 4,371, box `0.0468` of `0.862` credited (`5 %`).
**Row 26** (`z = 0.1783`, 27 members): the tightest of the three internally — eight of its 351
pairs are below `1e-3` and the `1 %` quantile of its pair margins is `+0.00058` — core deficit
`-0.0019`, box `0.0371` of `0.788` (`5 %`).

Top ten `E2Pg` rows (members subsampled to 400, so 79,800 pairs each):

| rank | row | `z_K` | `\|K\|` | core deficit | least pair margin | grazing (`1e-6`) | max graze degree | `vol(B_K)` | `vol(A(K))` |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 574 | `0.099514` | 1670 | `-0.085912` | `+0.000062` | **0** | 0 | `9.36e-02` | `0.757` |
| 1 | 576 | `0.089341` | 1683 | `-0.085659` | `+0.000062` | **0** | 0 | `6.12e-02` | `0.835` |
| 2 | 53 | `0.082123` | 1682 | `-0.096669` | `+0.000785` | **0** | 0 | `4.91e-02` | `0.777` |
| 3 | 5 | `0.082123` | 1097 | `-0.000327` | `+0.000000` | **132** | 33 of 400 | `1.72e-01` | `0.772` |
| 4 | 220 | `0.080180` | 1119 | `-0.042120` | `+0.000000` | **15** | 15 of 400 | `1.19e-01` | `0.800` |
| 5 | 354 | `0.073162` | 1778 | `-0.085092` | `+0.000056` | **0** | 0 | `7.21e-02` | `0.770` |
| 6 | 492 | `0.068992` | 1785 | `-0.064030` | `+0.000566` | **0** | 0 | `2.86e-03` | `0.722` |
| 7 | 187 | `0.067094` | 1776 | `-0.090415` | `+0.000346` | **0** | 0 | `7.81e-02` | `0.764` |
| 8 | 126 | `0.065550` | 1643 | `-0.099752` | `+0.000452` | **0** | 0 | `1.04e-01` | `0.803` |
| 9 | 588 | `0.065437` | 1057 | `-0.000218` | `+0.000000` | **153** | 31 of 400 | `6.18e-02` | `0.731` |

Here the two species separate cleanly.  Rows `574, 576, 53, 354, 492, 187, 126` have core deficits
of `-0.064` to `-0.100` and **no grazing pair at all** — deeply, robustly non-Helly.  Rows
`5, 220, 588` have core deficits of `-0.0003` to `-0.042` — nearly Helly, essentially point
cliques — and they are the ones with a genuine **fan**: `132`, `15`, `153` grazing pairs, with a
single member grazing `33` (resp. `15`, `31`) of the 400 sampled others, exactly the
"point clique plus a fan of grazing contacts" of `CLIQUELEVER.md` §6.  So the review's picture is
right for the *near-Helly* rows and wrong for the rest, and the rest carry most of the mass: of
the ten heaviest `E2Pg` rows, seven (`0.546` of the `0.768` in the table) have no grazing pair.

So: the non-Helly excess that makes these rows a lever (`CLIQUELEVER.md` §0) lives in a handful of
thousandths-tight pairs per row, **not** in a fan of exact grazing contacts — no pair on any of
these rows is within `1e-6` of touching, and the median pair overlaps by `0.71–0.95`.  The
`1e-6` figure matters because `cliquelever` snaps pose centres to a `1e-6` grid (`Dc = 1e6`), so a
genuine grazing structure would have shown up at exactly that scale; it does not.  These rows are
*robust* non-Helly families — the same object as §4's triple, with 27 to 1,806 members instead of
three, and the fan is a feature of the near-point-cliques rather than of the rows that carry the
mass.  That is also why `HONEST.md` §7's worry ("the boxes can only be as wide as the pairwise
meeting slack allows, which on these rows is thousandths, so the block count explodes") is only
half right: a *single* box of side `0.59 x 0.53 x 0.48` (the median over the 317 rows) fits, because the tight pairs are tight between
*members*, not between a member and a fresh pose.

## 4. A robust non-Helly triple at `t = 4` (`allmeet.py helly`)

Let `n_k`, `k = 0,1,2`, be unit vectors at `120°`, so `n_0 + n_1 + n_2 = 0`.  For `d > 0` the half
planes `H_k = {x : ⟨x, n_k⟩ ≤ -d}` have empty triple intersection (a common point would give
`0 = Σ_k ⟨x, n_k⟩ ≤ -3d < 0`) while any two meet in a full wedge whose apex is within `O(d)` of
the origin.  Take `S_k` = the unit square with one **edge** on the line `⟨x, n_k⟩ = -d`, lying on
the `H_k` side; then `S_k ⊆ H_k`, so the triple intersection is empty, and consecutive squares
overlap in a region of positive area.  Translate by `(2,2)` into `[0,4]²`.

At `d = 1/20` the three exact rational poses are

    S_0 = (29/20, 2, u = 0)                       theta = 0°
    S_1 = (91/40, 761843/500000, u = 267949/1000000)    theta = 29.99998°
    S_2 = (91/40, 1238157/500000, u = 11547/20000)      theta = 59.99998°

with pairwise meet margin `+0.358013` (and `leaf_ceiling.sq_meets_sq` = True on all three pairs,
independently), and the triple intersection empty by an exact clip of all 12 half planes.

Thickened to three pose boxes of half-sides `0.015` in `(c_x, c_y)` and `0.01` in `u` — volume
`1.8·10⁻⁵` each, so the family `F = B_0 ∪ B_1 ∪ B_2` has **positive measure** in pose space — the
following are certified exactly:

1. every pose of every `B_k` is admissible (`box_admissible`);
2. every two poses of `F` closed-meet: inside a box by `box_pairwise_meets`, across boxes by
   `boxes_meet` (the largest corner-to-corner centre separation between two centre rectangles is
   at most `1`);
3. **no point lies in all three squares, for every choice of one pose from each box.**  With
   `C_k = max over the centre rectangle of ⟨c', e_k⟩ + w_hi/2` (exact; `e_k` an exact rational edge
   normal of `S_k`, `w_hi` from `w_hi_exact`) every square of `B_k` lies in `{⟨x, e_k⟩ ≤ C_k}`, and
   with `λ = (e_1×e_2, e_2×e_0, e_0×e_1) ≥ 0` — which satisfies `Σ λ_k e_k = 0` identically — the
   Farkas value is `Σ λ_k C_k = -0.049958 < 0`.  A common point `x` would give
   `0 = ⟨x, Σ λ_k e_k⟩ = Σ λ_k ⟨x, e_k⟩ ≤ Σ λ_k C_k < 0`.

So `F` is a positive-measure, pairwise closed-meeting family of admissible unit squares in
`[0,4]²` with no common point.  **`TODO.md`'s Helly-type statement is false as literally stated.**

(The construction survives to `d = 1/4`, where the nominal pairwise margin is still `+0.058`; what
stops the note quoting a larger `d` is only that `boxes_meet`'s sufficient criterion needs the
centre separation `√3(d + 1/2) ≤ 1`, i.e. `d ≤ 0.077`.  A sharper box-to-box test would extend it.)

## 5. What to conjecture instead

The measurement rules out the family-level statement and points at a mass-level one.  Two
observations constrain what can be true.

* Sound crediting of a clique row is *possible* with positive volume — 955 exactly certified
  boxes over the two instances, every row of both — and the resulting rows are genuinely non-Helly (`K ∪ B_K` has an empty core, since `K` already does), so
  there is no Helly obstruction to building them and no reduction of them to point cliques.
* What fails is quantitative: the sound credited set is a `0.4–22 %` slice of the generous one, the
  worst pose of the container lies in **no** box (and on `SUPEPGF` captures literally nothing but
  `0.0076` of polygon mass), so the dual's `62–64 %` clique mass buys nothing exactly where the
  cover needs it.  And by §0 item 1 that is not an artefact of choosing boxes: *no* sound rule can
  push `min capture` above `HONEST.md`'s `meet` minimum, which is already `0.585`/`0.655`.

The statement worth proving is therefore about the LP, and the honest form of it is a conjecture,
not a theorem:

> **Conjecture (cover-side ceiling at `t = 4`).**  Let `w = (y, z, ζ)` be a dual of the `t = 4`
> packing LP of `CLIQUELEVER.md` §1 + `RANKDIAG.md` §8 on a finite admissible pose set, and credit
> each clique row `K` by any rule whose credited set contains `K` and is pairwise closed-meeting
> (equivalently: any rule under which `Θ(w)` is a valid weighted cover of `[0,4]²`).  Then
> `Θ(w) / min_S capture_w(S) ≥ 12`; equivalently, no weighting of `[0,4]²` assembled from point
> atoms, odd-polygon rows and soundly-credited clique rows certifies a bound better than `12`.

What a proof would need, and what this measurement does *not* supply: the conjecture quantifies
over all pose sets and all duals, while the measurement is two duals on two pose sets.  A proof
would have to (a) turn "the credited set is pairwise meeting" into a bound on its **measure** in
pose space — §4 shows it cannot be turned into "it has a common point", so the bound has to be
metric, e.g. that a pairwise-meeting family of unit squares has all its centres in a set of
diameter `≤ 1 + o(1)` and hence occupies at most a fixed fraction of `[0,4]²`; and then (b) sum
those measures against `Θ` to get a covering-density inequality of the kind `RANKDIAG.md` §8 proves
for the polygon rows.  Step (a) is the one with content, and §4's triple is the reason it must be
stated for *measure* and not for *intersection*.  A weaker and possibly easier statement, which
these numbers also support and which is enough to close the route, is: **on any dual whose clique
rows have empty cores, the minimum of `capture` over the admissible poses is attained at a pose
that lies in no sound clique row's credited set, so the clique mass contributes `0` there and the
honest cost is `Θ / (point + polygon capture)`** — which on these two duals is `145` to `1543`.
Both `--refit` sequences are evidence for it and not proof: re-placing the boxes to credit the
current worst pose always produced a new worst pose crediting `0` clique mass, on every round, on
both instances.

## 6. What is exact and what is float

**Exact** (`Fraction`/integer): the pose set and the replay of `cliquelever.Poses`'s snapping
(asserted to match the dual dump bit for bit) and `fold_pose`'s quarter turn (asserted to preserve
the exact corner set, on all 15,567 poses of the two instances); `box_meets_square`,
`box_pairwise_meets`,
`boxes_meet`, `box_admissible`, `w_lo_exact`, `w_hi_exact`, `triple_empty_certificate`, and
therefore **every `vol(B_K) > 0` in §2 and the whole of §4**.  `report --verify` re-runs
`certify_box` on every box in the json from scratch — `317 / 317` and `638 / 638` re-certified,
0 failures.
The exact tests were cross-validated against `leaf_ceiling.sq_meets_sq` (§1).

**Float**: the LP and its dual (HiGHS doubles, as `HONEST.md`); the max-min margin LPs and the `θ`
sweep that *propose* a pose; `FloatBoxTest`, which only steers the growth; the `vol(A(K))` column
(Monte Carlo, 12,000 uniform admissible poses, so ±2 % relative); `core_deficit`; the pairwise
margin distributions of §3; and the whole of the honest-cost measurement in §2 — `capture`'s point
term is `packing_dual.capture`'s C kernel with its `1e-8` containment tolerance and the polygon
term is `rankfamily.pgon_cost`'s, exactly as in `HONEST.md` §5.  Every slack is in the direction
that makes `capture` larger and `honest` smaller, i.e. the numbers err on the side of the cover
looking good.

Nothing here is certified and nothing here bounds `s(12)`.  §4's triple is an exact geometric
fact; everything else is a packing-side measurement.

## 7. Reproduce

```bash
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-12/cl_SUPEPGF_* runs/
cp /home/evand/math/square-packing/s12/runs/inputs-2026-09-12/cl_E2Pg_*   runs/
python3 search/family_rows.py runs/hc_families.txt --pitch 0.004 --tilted   # 61,656 poses

# the duals (13 s / 1,100 s), reproducing HONEST.md's Theta to nine places
python3 search/honestcost.py dual SUPEPGF --threads 4
python3 search/honestcost.py dual E2Pg    --threads 6

# the boxes: max-min margin sweep + exact growth, every box exactly certified
python3 search/allmeet.py boxes SUPEPGF                            # 242 s
python3 search/allmeet.py boxes E2Pg --nth 46 --nanat 24 --cap 400 # 817 s

# the honest cost under `credit z_K iff S in K u B_K`, + 2 adversarial re-placements
python3 search/allmeet.py honest SUPEPGF --threads 4 --refit 2
python3 search/allmeet.py honest E2Pg    --threads 4 --refit 2

# the tables of this note (re-certifies every box from the json)
python3 search/allmeet.py report SUPEPGF --top 10 --measure 12000
python3 search/allmeet.py report E2Pg    --top 10 --measure 8000

# the robust non-Helly triple
python3 search/allmeet.py helly --d 0.05 --half 0.015 --halfu 0.01
```
