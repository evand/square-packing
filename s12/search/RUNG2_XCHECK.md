# An independent re-implementation of the zero-margin checker (task rung2-xcheck, 2026-09-12)

## 0. Verdict

> **The rung-2 certificate is independently verified.**  `verify2/` (`zmcheck`, Rust, exact
> `i128` on the whole load-bearing path) reports **0 uncertified boxes** on
> `certificates/rung2/s13_closed_cover_4.txt` over the **full** admissible pose domain
> (`c_x, c_y ∈ [0,4]`, `u = tan(θ/2) ∈ [0,1]`, i.e. `θ ∈ [0°,90°]`, **no symmetry reduction**),
> with its own subdivision, its own primitive set and no code shared with `search/zeromargin.py`.
>
> ```
> zmcheck cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads 8
> done in 1017s: boxes 30258, max depth 10
>   leaves: ADM 5114  DISJ 9477  EMPTY 6938  UNCERTIFIED 0
> VERIFIED: every closed unit square in [0,4]^2 captures weight >= 1;
>           total weight 12955972155/1000000000 = 12.955972155
> ```
>
> 17 minutes on 8 threads, and the depth limit is never reached (max depth 10 of 18), so the
> subdivision terminated on its own.  **So `s(13) = 4` now rests on two independent exhaustive
> verifications of the same weighted cover**, by two programs that share no code, no subdivision
> rule and no primitive set, and that were written from the statement rather than from each other.
>
> `DISJ` — the disjunctive primitive — carries **9 477 of the 14 591 non-empty leaves (65 %)**,
> against `CHAIN`'s 5 320 of 8 187 (65 %) in `RUNG2.md` §0.  The same 65 %, from a different
> disjunction over a different subdivision, is the clearest confirmation of Theorem 1 available:
> no cover of weight `12.956 < m² = 16` can be certified by fixed-witness primitives, so a
> disjunctive certificate has to carry the bulk of the proof whatever the checker looks like.

Two further pieces of agreement, both worth more than the box census because they are pointwise and
exact:

* `zmcheck pose` computes the exact captured weight of one rational pose with **no boxes, no
  subdivision and no primitives** — just the four inequalities per point.  On the four poses for
  which `search/RUNG2.md` publishes an exact `Fraction` value it reproduces every digit, and the
  captured-point count where the note gives one:

  | file | pose | `RUNG2.md` | `zmcheck pose` |
  |---|---|---|---|
  | `closed4_best_x103.txt` | `(3/2, 1461/2000, u = 1/40000)` | `970282351/1000000000` (§4.3, 175 pts) | `970282351/1000000000`, 175 pts |
  | `closed4_best.txt` | same | `9420217/10000000` (§4.3) | `9420217/10000000`, 175 pts |
  | `closed4_best_x103.txt` | `(7/2, 3/2, 0)` | `1.562914790` (§4.6) | `1562914790/1000000000` |
  | `closed4_best_x103.txt` | `(7/2, 3/2 + 1/320, 0)` | `1.036922939` (§4.6) | `1036922939/1000000000` |

  Two implementations of the *semantics* (closed square, closed containment, `u`-parametrised
  rational rotation, admissibility) that agree to the last digit of a 9-digit rational on a pose
  where one point sits exactly on `∂Q` is the strongest single check available; it is independent
  of every lemma below.

* The shipped cover is D4-symmetric enough that `(1/2, 1/2, 0)` and `(7/2, 7/2, 0)` must capture
  the same weight, and both come out as `1050006300/1000000000`.

The rejection tests are in `tests/rung2/rejection_tests.sh`; §5 reports what each mutation does and
why.

---

## 1. What was built, and why it is not the same checker

`verify2/src/main.rs` (`zmcheck`), one crate, no dependencies, no change to `verify/` or
`verify.sh`.  Modes: `cert` (the sweep), `pose` (one exact rational pose), `box` (one box,
verbose).

The design constraint I set myself was to have **one** polynomial class, **one** exact-maximum
routine and **one** inference rule, and to derive every primitive of `ZEROMARGIN.md`/`RUNG2.md`
as a special case rather than implementing them separately.  That is a genuinely different
architecture from `zeromargin.py`'s `cert_adm` / `cert_p1` / `cert_mix` / `cert_chain` /
`clip_bin` / `EMPTY`, and it is what makes the re-implementation worth something: a shared bug
would have to be a bug in the *statement*, not in the code.

### 1.1 The class

Pose `(x, y, u)` with `u = tan(θ/2) ∈ [0,1]`; `C = 1-u² ≥ 0`, `S = 2u ≥ 0`, `N = 1+u² > 0`,
`cos θ = C/N`, `sin θ = S/N`, `w(θ) = cos θ + sin θ = (C+S)/N`.  For a point `p`, `a = p_x - x`,
`b = p_y - y`:

    X = (aC + bS)/N,   Y = (-aS + bC)/N,   p ∈ Q(x,y,θ)  ⟺  |X| ≤ ½ and |Y| ≤ ½.

Multiplying the four inequalities by `2N > 0` gives `RUNG2.md` §6.2's violation polynomials

    G_{p,0} =  2aC + 2bS - N      G_{p,1} = -2aC - 2bS - N
    G_{p,2} = -2aS + 2bC - N      G_{p,3} =  2aS - 2bC - N ,

`p ∈ Q` iff all four are `≤ 0`.  **Admissibility is the same kind of object**: `Q(x,y,θ) ⊆ [0,m]²`
iff `x, y ∈ [w/2, m - w/2]` (`ZEROMARGIN.md` §1), i.e. iff

    A₁ = C + S - 2xN ≤ 0,   A₂ = 2xN - 2mN + C + S ≤ 0,   A₃, A₄ likewise in y.

So every inequality the checker reasons about lies in

    P₄ = { F(x,y,u) : affine in (x,y), degree ≤ 4 in u, rational coefficients },

a cone closed under nonnegative combinations.  (`G` and `A` are degree 2; degree 4 appears only
when a multiplier `C`, `S` or `N` is used — §2.3.)

### 1.2 The two facts

**Lemma D (exact maximum over a pose box).**  *Let `F ∈ P₄` and
`B = [x₀,x₁] × [y₀,y₁] × [u₀,u₁]`.  Then `max_B F = max` over the four corners of the centre
rectangle of `max_{u ∈ [u₀,u₁]} F`.*

*Proof.*  For each fixed `u`, `F` is affine in `(x,y)`, so its maximum over the rectangle is at a
corner; hence `max_B F = max_u max_{corners} F = max_{corners} max_u F`. ∎

This is `RUNG2.md` Lemma E, proved the same way, but here it applies to *every* test in the
checker, not only to the `CHAIN` combinations. The inner univariate maximum is computed as follows.

**Lemma D′ (the univariate test).**  *Let `q(t) = Σ_{j≤4} c_j t^j` have integer coefficients.  Put
`β_i = Σ_{j≤i} c_j C(i,j)/C(4,j)` for `i = 0..4`.  Then `max_{[0,1]} q ≤ max_i β_i`, with
`β₀ = q(0)` and `β₄ = q(1)`.*

*Proof.*  `q = Σ_i β_i B_{i,4}(t)` with `B_{i,4} ≥ 0` and `Σ_i B_{i,4} = 1` on `[0,1]`, so `q` is a
convex combination of the `β_i`. ∎  (`RUNG2.md` Lemma C, in the local parameter.)

The implementation decides `max_{[0,1]} q ≤ 0` by: the two endpoints exactly (always); for
`deg q ≤ 2` the vertex `-c₁/(2c₂)` when `c₂ < 0` and the vertex is interior, which makes the test
**exact**; for `deg q ∈ {3,4}` the Bernstein bound of Lemma D′, refined by two levels of exact
bisection of `[0,1]` (`q(s/2)` and `q((1+s)/2)`, rescaled by `16` to stay integral).  The Bernstein
branch is a *sound over-estimate*, so it can only lose a certification.  Two remarks:

* Every box is worked in a **local** parameter `t ∈ [0,1]`, `u = (v₀ + h t)/M`.  That is not
  cosmetic: in the global `u` the Bernstein coefficients carry a factor `u₀^4 ≈ M⁴` on top of the
  coefficients' own `M⁴`, and `i128` overflows by depth 6.  In the local parameter the only scale
  is the coefficients' own, and the a-priori bound is
  `160 · Dc · M⁴ · 12 · 16²  <  2^(47 + d_xy + 4 d_u)` with `Dc = 1000·2^{d_xy}` the centre
  denominator and `M = 8·2^{d_u}` the angle denominator.  The checker **refuses** (and reports) any
  box with `4 d_u + d_xy > 74`, so no verdict can depend on wrap-around.  The shipped run never
  came near it (the deepest box in any run reported here is depth 9).
* Lemma D′ is exact, not merely sound, at a double root sitting on an endpoint — which is exactly
  the corner-region case (§2.3) where the margin is `0`.

**Lemma I (inference).**  *Let the hypotheses at a node be `F_i ≤ 0` with `F_i ∈ P₄`, and let
`F ∈ P₄`.  If there are rationals `μ_i ≥ 0` with*

    max_B ( F - Σ_i μ_i F_i )  ≤  0 ,

*then `F ≤ 0` at every pose of `B` that satisfies all the hypotheses.*

*Proof.*  At such a pose, `μ_i F_i ≤ 0` for every `i`, so
`F ≤ F - Σ μ_i F_i + Σ μ_i F_i ≤ 0 + 0`. ∎

Everything else in the checker is a choice of hypotheses `F_i` and multipliers `μ_i`.  The
multiplier need not be constant: any `μ_i` that is a polynomial in `u`, nonnegative on `[0,1]`,
keeps the combination in `P₄` — the implementation offers `μ ∈ {0, C, S, N}` times a nonnegative
rational, and `C, S, N ≥ 0` on `u ∈ [0,1]` exactly because the domain is `θ ∈ [0°,90°]`.

---

## 2. Every primitive as a case of Lemma I

Throughout, "condition `k` of `p` is certified at a node" means Lemma I proves `G_{p,k} ≤ 0`
there; `p` is a **witness** when all four of its conditions are certified, and then `p ∈ Q` at
every pose of the node.  A node is discharged when its witnesses have total weight `≥ 1`.  Weights
are integers over `W = 10⁹`, so "`≥ 1`" is the integer comparison `Σ w_p ≥ W`.

### 2.1 `EMPTY` and `clip_bin`, for free

`EMPTY` (`ZEROMARGIN.md` §2) is "no admissible pose in `B`".

**Lemma E0.**  *If `max_B(-A_i) < 0` for some `i`, no pose of `B` is admissible, and the node is
discharged vacuously.*

*Proof.*  `-A_i < 0` on `B` means `A_i > 0` on `B`, and `A_i ≤ 0` is necessary for admissibility. ∎

(The strict test is the same routine with `>` replaced by `≥`; a strict *upper bound* is still
sound in this direction.)  `A₁` is decreasing in `x`, so Lemma D evaluates it at `x = x₁`, which is
exactly `ZEROMARGIN.md` §2's `cx₁ < w_lo/2` test — except that the correct `w` per angle is used
instead of the bin's minimum, so the test is strictly stronger.

More importantly, the admissibility inequalities are hypotheses of **every** node, so the
checker never has to prove anything at an inadmissible pose.  That is `RUNG2.md` §4.6's
`clip_bin` — the completeness fix without which "*no* box pressed against a wall is certifiable at
any depth" — obtained structurally rather than as a separate bisection step.  A box whose only
admissible poses are at `u = u₀` is handled because the hypotheses pin `u` there; nothing has to be
clipped.

### 2.2 `CORE` (all multipliers zero)

With `μ_i = 0`, Lemma I is "`max_B G_{p,k} ≤ 0`", i.e. `p` lies in every square of the box.  This
is the exact bin core of `ZEROMARGIN.md` §2 in a different presentation: that note intersects
`R_{θ₀}Q ∩ R_{θ₁}Q` with a sector condition, whereas here the four inequalities are maximised over
the *whole* bin by Lemma D. The two are the same set (both are "`|X| ≤ ½, |Y| ≤ ½` for all `θ` in
the bin"); Lemma D needs no sector test and no `|v|² ≤ ¼` test.

### 2.3 `ADM` / Lemma A, and `P1`

`RUNG2.md` Lemma A says: with `A_x(θ) = max(x₀, w/2)`, `B_x(θ) = min(x₁, m - w/2)` (and likewise in
`y`), the four conditions need only be checked at four *different* corners, the extreme admissible
ones for each.  In the present setup that is automatic and exact:

**Lemma A′.**  *For condition `k`, `G_{p,k}` is monotone in `x` and in `y`, with the signs*

| `k` | `∂/∂x` | worst `x` | `∂/∂y` | worst `y` |
|---|---|---|---|---|
| 0 | `-2C ≤ 0` | smallest | `-2S ≤ 0` | smallest |
| 1 | `+2C ≥ 0` | largest  | `+2S ≥ 0` | largest |
| 2 | `+2S ≥ 0` | largest  | `-2C ≤ 0` | smallest |
| 3 | `-2S ≤ 0` | smallest | `+2C ≥ 0` | largest |

*so "`G_{p,k} ≤ 0` at every admissible pose of `B`" is exactly `G_{p,k} ≤ 0` evaluated at the
extreme admissible centre named in the table, for every `θ` of the bin — Lemma A verbatim.  Taking
the multiplier that eliminates that centre coordinate turns this into an instance of Lemma I:
with target `N·G_{p,k}` and*

| `k` | `x`-hypothesis | `μ` | `y`-hypothesis | `μ` |
|---|---|---|---|---|
| 0 | `A₁` | `C` | `A₃` | `S` |
| 1 | `A₂` | `C` | `A₄` | `S` |
| 2 | `A₂` | `S` | `A₃` | `C` |
| 3 | `A₁` | `S` | `A₄` | `C` |

*the coefficient of that coordinate in `N·G_{p,k} - μ A` vanishes identically, and the combination
is the exact Lemma-A quantity with `A_x = w/2` (resp. `B_x = m - w/2`) substituted.*

*Proof of the elimination.*  For `k = 0`: the coefficient of `x` in `N·G_{p,0}` is `-2NC` and in
`A₁` it is `-2N`, so in `N·G_{p,0} - C·A₁` it is `-2NC + 2CN = 0`; substituting `x = w/2` makes
`A₁ = 0`, so the combination equals `N·G_{p,0}|_{x = w/2}`.  The other rows are identical. ∎

`max(x₀, w/2)` is not a polynomial in `u`, so — exactly as `RUNG2.md` §3.2 says — the
implementation offers each of the two lower bounds separately (`μ = 0` gives `x = x₀`, `μ = C`
gives `x = w/2`) and accepts the condition if **any** of the nine combinations
`μ_x ∈ {0, C or S, N}` × `μ_y ∈ {0, S or C, N}` certifies it.  Two of those nine deserve names:

* `μ_x = μ_y = 0` is `CORE` (§2.2);
* `μ = N`, i.e. `μ = 1` after dividing by `N`, is **`P1`**.  For instance with `p_x ≤ 1`,
  `N·G_{p,0} - N·A₁` has `x`-coefficient `-2NC + 2N²`, and the whole combination is
  `2N²(p_x - 1) + (terms in y) ≤ 0` — the `ZEROMARGIN.md` §2 statement that "a side on which `p` is
  within 1 of the container wall is automatic".  Concretely, at the corner tile the combination
  `N·G_{(1,1),0} - C·A₁ - S·A₃` is *independent of the centre* and equals

        N · G_{(1,1),0}|_{x = y = w/2}  =  -N²(w-1)²  =  -4u²(1-u)²  ≤  0 ,

  so the point `(1,1)` is certified over the **whole** corner region `[0,·]²×[0,·]` in one box, at
  any bin width, with the margin `0` that `ZEROMARGIN.md` §4 item 1 shows is unavoidable (it is
  quadratic in `θ`, and `-4u²(1-u)²` has a double root at `u = 0`, which is why Lemma D′ has to be
  exact there — and is).

So `ADM ⊇ CORE ∪ P1` is not a design decision here but an identity between multiplier choices, and
`MIX` (`RUNG2.md` §3.4: the union of what `ADM` and `P1` certify) is not a separate primitive
either, since the witness set is computed per point and per condition with all nine choices
available.  The one thing Lemma A′ gives up relative to Lemma A is the *per-angle* choice between
`x₀` and `w/2`; that loss is `2N(cos θ)(w/2 - x₀)` where the two bounds cross, which vanishes under
subdivision.

### 2.4 `DISJ`: the disjunctive primitive, and why it needs no chains

`RUNG2.md` Theorem 1 is the reason a disjunctive primitive is mandatory: a cover certified
everywhere by *monotone witness certificates* — a fixed witness set per leaf, which is precisely
what §2.1–2.3 produce — must weigh at least `m² = 16`, and this one weighs `12.956`.  I did not
re-prove Theorem 1 (it is proved in `RUNG2.md` §2 and is not needed for soundness; it only says
that §2.3 alone cannot finish, which the measurements confirm).  What I did was re-derive the
primitive.

`CHAIN` (`RUNG2.md` §6) builds a monotone chain `G_{q_1} ≤ … ≤ G_{q_k}` on the box, gets `k+1`
nested regions from Lemma F, populates them by Lemma G, multiplies two chains for interior tile
poses, and needs Lemma H to show some product regions empty, plus an up-set/down-set bookkeeping
(the "suffix down-set" of §0 and §4.1 — a rule whose section, §6.4, does not exist in the note;
see §4).  I replaced all of that with plain recursive **sign splitting**, which is both simpler and
strictly more general:

**Lemma S (the split).**  *For any `F ∈ P₄`, `B = (B ∩ {F ≤ 0}) ∪ (B ∩ {-F ≤ 0})`.*  ∎

Both halves are *closed*, and each is described by adding one hypothesis of the standard form
`F_i ≤ 0` to the node, so Lemma I applies unchanged.  The recursion is: compute the node's witness
weight by Lemma I; if `< 1`, pick a violation polynomial `G_{q,k}` and recurse on both halves; a
node is discharged when its weight reaches `1`, and the parent is discharged when both children
are.  Splitting on `G_{q,k}` is useful because

* in `{G_{q,k} ≤ 0}` the point `q` gains its missing condition immediately (`μ = 1` on the new
  hypothesis makes the combination identically `0`);
* in `{-G_{q,k} ≤ 0}` any other point `a` with `max_B(G_{a,k'} + λ G_{q,k}) ≤ 0` for some `λ > 0`
  gains its missing condition — that is **`RUNG2.md` Lemma G**, here the Lemma-I instance with
  hypothesis `-G_{q,k} ≤ 0` and multiplier `λ`;
* `RUNG2.md` **Lemma F** is not needed: what a chain's monotonicity buys (the region ordering) the
  recursion gets by branching, and a branch that a chain would have proved impossible is instead
  *discharged*, because —

**Lemma K (degenerate node = `RUNG2.md` Lemma H, without emptiness).**  *At a node whose
hypotheses include `-G_{q} ≤ 0` and `-G_{q'} ≤ 0`, suppose `max_B(G_q + λ G_{q'}) ≤ 0` for some
`λ > 0`.  Then every pose of the node has `G_q = G_{q'} = 0`; in particular both conditions hold
and `q`, `q'` are witnesses there.*

*Proof.*  `G_q ≥ 0`, `G_{q'} ≥ 0` and `G_q ≤ -λ G_{q'} ≤ 0` force `G_q = 0`, then
`λ G_{q'} ≤ -G_q = 0` forces `G_{q'} = 0`. ∎

And this is *again* just Lemma I (certify `G_q ≤ 0` from the hypothesis `-G_{q'} ≤ 0` with
multiplier `λ`), so the implementation contains no separate emptiness test, no staircase search and
no product-of-two-chains construction: two independent sliding cuts are simply two successive
splits, and the contradictory quadrant is discharged by Lemma K rather than proved empty.  Chains
are recovered as the special case in which one child of every node closes immediately — the
"caterpillar" shape — which is what actually happens at the wall poses (see the region counts in
§3).

Multipliers offered for a branch hypothesis: `λ ∈ {1, ½, 2}`, each combined with the nine
admissibility multiplier choices of §2.3 (the wall poses need both at once: `(1/2, 3/2, 0)` is a
wall pose *and* a sliding cut).

### 2.5 Subdivision, and what the floats do

A box that no primitive discharges is halved: the longest of `x`-extent, `y`-extent and
`2·bias·u`-extent (`dθ/du = 2/(1+u²) ∈ (1,2]`, and `bias = 4` at a wall — `RUNG2.md` §3.4's
`--theta-bias`), down to `--depth` (18; never reached).  Root boxes: pitch `1/10` in `x` and `y`
(commensurable with the certificate's `D = 1000`, which `ZEROMARGIN.md` §4 item 3 shows is
necessary — the tight families sit on the coordinate lattice) and 8 bins of `u ∈ [0,1]`:
`40 × 40 × 8 = 12 800` roots.

Floats appear in exactly three places, all of them screens or heuristics that can only *lose* a
certification:

1. **reach filter.**  `p ∈ Q ⇒ |p - c|_∞ ≤ w/2 ≤ √2/2 = 0.70711`, so points further than `0.7072`
   from the centre rectangle are not tested at all.  A superset by `10⁻⁴`.
2. **monotone pre-screen.**  For each candidate point and each condition, the *ideal* Lemma-A
   quantity (§2.3, with the true `max(x₀, w/2)`) is evaluated at five angles of the bin.  If it is
   `> 10⁻⁹` at one of them, some admissible pose of the box fails the condition, so **no**
   multiplier choice can certify it and the exact tests are skipped.  This is a lower bound on the
   ideal maximum, hence sound in the losing direction, and it is what makes the checker fast.
3. **disjunctive pre-screen and ranking.**  The same idea for Lemma I with a branch hypothesis:
   the box's four centre corners, clipped to the admissible range, crossed with five angles, give
   sampled admissible poses; if one of them satisfies the hypothesis and still violates the
   condition, no `λ` can work.  The choice of which polynomial to split on is ranked first by this
   float table and then, for the surviving candidates, by the **exact** weight of the two children.

Everything that certifies a box is integer `i128`: the witness sets, every Lemma-I test, the
weight comparisons.  A float can make the checker slower or weaker, never wrong.

---

## 3. What ran, what it says

### 3.1 The verification

`zmcheck cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads 8`, full domain,
`--dump runs/xcheck_leaves.txt`:

| | roots | boxes | max depth | ADM | DISJ | EMPTY | uncertified | wall |
|---|---|---|---|---|---|---|---|---|
| **full admissible domain** | 12 800 | **30 258** | **10** of 18 | 5 114 | 9 477 | 6 938 | **0** | **1 017 s** |
| `zeromargin.py` for comparison (`RUNG2.md` §0) | — | 16 872 | 13 of 18 | 2 867 | 5 320 (`CHAIN`) | 3 449 | 0 | 6 138 s |

The two censuses are *not* meant to agree — different root grid (`zmcheck` sweeps the unreduced
domain: `12 800` roots against the Python's symmetry-reduced sweep), different splitting rule,
different primitives.  What does agree, and is the point, is the **share carried by the
disjunctive primitive**: `9 477 / 14 591 = 65 %` here, `5 320 / 8 187 = 65 %` there.  `RUNG2.md`
Theorem 1 predicts exactly that a cover of weight `< 16` must be certified disjunctively almost
everywhere, and two unrelated implementations land on the same fraction.

The per-leaf witness sets are in `runs/xcheck_leaves.txt` (10 MB, one line per leaf: kind, exact
box coordinates, and for `ADM` the witness point indices).

### 3.1a Bands, for bisecting a residue

A centre band is selected with `--xlo/--xhi/--ylo/--yhi`; a banded run prints `PARTIAL SWEEP` and
never says `VERIFIED`, so the number to read is the uncertified count.  These were the validation
steps on the way, and they are the tool to use if a future cover does have a residue.

| band | roots | boxes | max depth | ADM | DISJ | EMPTY | uncertified | wall |
|---|---|---|---|---|---|---|---|---|
| `c_x ∈ [0, 0.35]` (all `c_y`, all `u`) | 960 | 960 | 0 | 0 | 0 | 960 | **0** | 0 s |
| `c_x ∈ [0.4, 0.75]` (all `c_y`, all `u`) | 960 | 4 178 | 9 | 659 | 828 | 1 082 | **0** | 147 s |
| `c_x ∈ [0.5, 0.6]`, `c_y ∈ [1.0, 2.0]` | 88 | 1 106 | — | 327 | 432 | 347 | **0** | ~60 s |
| `c_x ∈ [0.5, 0.6]` (all `c_y`) | 320 | 1 158 | 9 | 199 | 195 | 345 | **0** | 75 s |
| `c_x ∈ [1.45,1.55]`, `c_y ∈ [1.4,1.6]` (the interior tile pose) | 32 | 129 | — | 0 | 129 | 0 | **0** | 45 s |
| `c_x ∈ [1.25,1.6]` (the two-chain region) | 960 | 3 356 | — | — | — | — | **0** | 211 s |

The third and fourth rows are the important ones: `c_x ≈ 1/2`, `c_y ≈ 3/2`, `θ → 0` is the pose
whose monotone witness set `RUNG2.md` §2 computes (`0.566119` for the earlier cover; `0.6574` of
the needed `1` for this one — measured by `zmcheck box`, see §3.2), i.e. the place Theorem 1
proves is unreachable by `CORE`/`P1`/`ADM` at any depth.  `DISJ` carries 432 of the 759 non-empty
leaves there and 828 of 1 487 over the whole left-wall band — 56 %, against `RUNG2.md`'s 65 % for
`CHAIN`, which is the same qualitative statement: the disjunctive primitive is the main
certificate type, exactly as Theorem 1 requires.

### 3.2 Where it stopped *before* the fix of §3.3, and why

*Kept as the record of the failure mode, because it is the one place where `RUNG2.md`'s explicit
chain construction does something a generic search does not.*  Before §3.3, a full-domain run had
to be stopped at 21 min: it certified everything up to `c_x ≈ 1.25` with 0 uncertified boxes and
then began subdividing the interior tile poses to the depth limit.

The residue is the interior tile poses.  Take the `σ = (+,-)` octant of `(3/2, 3/2, 0)`, at depth
10 in the subdivision:

```
zmcheck box certificates/rung2/s13_closed_cover_4.txt \
        --box "1500/1000,1501/1000,1499/1000,1500/1000,0,1/1024"
box x[1.500,1.501] y[1.499,1.500] u[0,1/1024]: reach=1212 T=277 w(T)=0.625480485
verdict: needs subdivision
```

`w(T) = 0.6255` is the *monotone* witness weight, and it barely moves as the box shrinks
(`0.5146` at the root box `[1.5,1.6] × [1.4,1.5] × [0°,7.2°]`, `0.5486` after the bin is cut to
`1/64`, `0.6255` at `1/1024`) — the signature of `RUNG2.md` Theorem 1: this is
`w(PIN(1,1,+,-))`, it is `< 1`, and no amount of depth changes that.  The geometry is the
predicted one: with `c = (3/2 + α, 3/2 - β)`, `θ` small, the rows `x = 2` and `y = 1` survive
and the rows `x = 1`, `y = 2` are lost, and

    point (2, y) is captured  ⟺  y ≲ c_y + α/θ ,
    point (x, 1) is captured  ⟺  x ≲ c_x + β/θ ,

two cuts sliding *independently* (`α`, `β` independent) — `RUNG2.md` §6.2's two-chain situation.

The primitive is strong enough.  The exact union over *all* 96 sign hypotheses of what Lemma I
certifies (`ZM_UB=1`) is

```
  exact-ub over ALL hypotheses = 2.965719330
```

so there is nearly `3×` the required weight reachable; and the checker does close the analogous
*wall* box at `(1/2, 3/2, 0)`, where one cut slides, with a 9-region certificate in 17 nodes.
What fails is the construction of the tree.  My search is a depth-first sign-splitting search whose
branch choice is scored (three heuristics tried in turn: maximise the weaker child, maximise the
`G ≥ 0` child — the chain order of §6.2 — or maximise the sum), and none of the three reaches `1`
in every region before the node budget, at any setting tried at that point (`--branch-cap` up to
96, `--fail-cap` up to 512, `--sign-depth` up to 40, `--node-cap` up to 3·10⁵).

The reason is structural and is worth recording, because it is the one place where `RUNG2.md`'s
chain construction is doing something a generic search does not.  In the two-cut region indexed
`(r, s)` the witnesses are `T`, the row `x = 2` below cut `r`, and the row `y = 1` left of cut `s`.
Both cuts are one-sided *in the same direction*, so the region `(0,0)` — both cuts at their
minimum — has witness set `≈ PIN`, weight `0.6255 < 1`, and it can only be discharged by showing
it **empty**: that is precisely `RUNG2.md` Lemma H, and in my framework Lemma K (§2.4).  So the
certificate is not "find enough weight in every region" but "find the two pivots whose
simultaneous violation Lemma K refutes, and chain outwards from them", and a heuristic that
maximises weight per split never looks for that pair.  The fix is to search for the Lemma-K pair
directly — for each pair `(q, q')` of swing pivots of different kinds, test
`max_B(G_q + λ G_{q'}) ≤ 0`, and seed the chains at the pair that succeeds — which is a bounded
`O(k²)` exact scan of the same primitive, not new mathematics.  §3.3 is that fix; it turned out
that the cleanest form of it is not an emptiness test at all but a *closure rule* on the node's
hypotheses, and that the truncated chain (`--branch-cap 48`) was a second, independent cause.

(Two sanity checks that this is the whole residue and not a symptom of something worse: the wall
poses, which are the *harder* ones by Theorem 1's own measure, all close; and `zmcheck pose`
confirms the cover itself is fine there — at `(3/2, 3/2, 0)` it captures
`2019213840/1000000000 = 2.019213840 ≥ 1` (782 points).)


### 3.3 The fix, and what it changed

Everything in §3.2 above is kept as the record of what failed; this section is what was done about
it.  Three changes, all in the *search*; the polynomial class, Lemma D, Lemma I and every exact
test are untouched.

1. **The Lemma-K closure** (`Disj::cover_child`).  A node's certified-condition cover used to be
   the union of its hypotheses' masks.  It now *closes*: whenever the node carries the hypothesis
   `G_q ≥ 0` for a candidate `q` **and** the cover already proves `G_q ≤ 0`, then `G_q = 0` at
   every pose of the node, so `G_q ≤ 0` may be **added as a hypothesis** and everything its `le`
   mask certifies — the whole down-set of the chain through `q` — becomes available.  Iterated to a
   fixpoint.  The premise "the cover already proves `G_q ≤ 0`" is precisely what `RUNG2.md`
   Lemma H establishes for a pair of pivots that cannot both be violated, and it is detected for
   free, because `q`'s own condition slot is a bit of the mask of the *other* pivot's `ge`
   hypothesis.  So the witness-free quadrant is *certified*, not shown empty, and the checker still
   contains no emptiness test.

   Soundness: the added hypothesis is a proved consequence of the node's hypotheses, so Lemma I
   applies to it verbatim.

2. **The forced two-chain order** (`Disj::chain_order`).  Within one kind `k`, the branch
   candidates are sorted by a float estimate of `max_B G_{q,k}` ascending — `RUNG2.md` §6.2's
   chain order — the two largest families are interleaved, and that whole sequence is *forced* as
   the branch order, with the search depth set to its length.  The recursion then visits exactly
   the product regions `R_r × R'_s` of §6.2.  `kpairs()` additionally enumerates the
   float-plausible Lemma-K pairs (`O(k²)`, confirmed exactly inside `cover`) and offers each as a
   two-element forced prefix (`--seed-cap`, default off now that the full order subsumes it).

3. Covers are carried down the recursion incrementally (a `Cov` = bitset + closure flags), and the
   three greedy heuristics get a small node budget (`--heur-cap`, 3 000) with the full budget
   reserved for the forced order.  Without that split, three failing greedy attempts cost 30× the
   successful one.

Two defaults had to be raised, and one of them was the whole difference between failing and
succeeding: **`--branch-cap` 48 → 160**.  With 96 candidates the two chains are *truncated* and
the worst region stalls at `0.984808` of the needed `1` — 1.5 % short, and visibly not improving
with more nodes, more depth or a larger `--fail-cap`; with 160 it closes.  (`--fail-cap` 256 →
1024, `--node-cap` 4 000 → 400 000.)

Measured:

| box | before | after |
|---|---|---|
| wall pose `(1/2,3/2,0)`, `[0.5,0.51] × [1.49,1.5] × [0,1/128]` | DISJ, 9 regions, 17 nodes, 12 ms | unchanged |
| interior tile pose `(3/2,3/2,0)`, `σ=(+,-)` octant at depth 10 | not certified at any setting | **DISJ, 386 regions, 9 837 nodes, 1.0 s** |
| depth-18 sweep of the centre cell `c_x ∈ [1.45,1.55]`, `c_y ∈ [1.4,1.6]` | — | 129 boxes, all `DISJ`, **0 uncertified**, 45 s |

The interior tile pose needing 386 regions against the wall pose's 9 is the two-cut structure
made visible: one sliding cut gives a chain, two independent ones give a product, and the product
is where `RUNG2.md`'s "56 regions reach 1.013756" measurement lives.

---

## 4. Where the note is ambiguous, wrong, or missing

These are the places where I had to look at `search/zeromargin.py` or make a decision the note does
not determine.  Each is a documentation bug.

1. **`RUNG2.md` §6.4 does not exist.**  §0 lists, among the two completeness fixes, "the suffix
   down-set in `CHAIN` (§6.4)", and §4.1 refers to "the down-set and the up-set" in the
   single-chain branch; but §6 ends at §6.3 and neither term is defined anywhere in the file.  The
   task brief asks the reader to read "§6 (`CHAIN`, Lemmas E–H, the down-set rule)" — the rule is
   not there.  **Reconstruction** (from §6.2 steps 3–4): for a chain `q_1 ≤ … ≤ q_k`, the set of
   `r` for which Lemma G certifies a given candidate `a` from `q_r` is a *prefix* `{1..r_a}`
   (because `max_B(G_a + λ G_{q_r})` is non-decreasing in `r`), so `a` is available in region
   `R_r = {G_{q_r} ≤ 0 < G_{q_{r+1}}}` iff `r + 1 ≤ r_a`, i.e. on a *down-set* of regions.  §6.2
   step 4's `U_{r+1}` — used but never defined — must then be `{a : r_a ≥ r+1}`.  My design needs
   none of this (§2.4), so I could not use the note's rule as a check on my own; I record it here
   as what §6.4 presumably said.
2. **`RUNG2.md` §6.2 step 4 uses an undefined symbol.**  `U_{r+1}` (above).  Related: step 1 says a
   swing point has "exactly one of whose four inequalities fails the `ADM` test", but step 4's
   accounting does not say what happens to a point with *two* failing inequalities; they can be
   certified by two hypotheses at once, which my recursion does and a single chain cannot.
3. **The task brief mis-states the certificate header.**  It says "line 1 `m sym`"; `FORMAT.md` and
   the file say `s_num s_den`.  For `s13_closed_cover_4.txt` the header is `4 1` and the two
   readings coincide, so nothing breaks here — but a checker written from the brief would mis-read
   every other shipped certificate (`s12_lower_3.9686.txt` begins `15680 3951`).  I followed
   `FORMAT.md`.  (`zmcheck` additionally refuses a non-integer container side rather than
   guessing, since its root grid is built on integer tiles.)
4. **Neither note states the pose domain of the shipped `cert` run.**  `ZEROMARGIN.md` §2 describes
   the root domain as "eight `u`-bins on `[0, 1/2]`, `cx ∈ [0,m]`, `cy ∈ [0,m/2]`" — a *reduced*
   domain, justified in §2 "Symmetry" by invariance under `x ↦ m-x` and `y ↦ m-y`, with `--full`
   for the unreduced one.  `RUNG2.md` §0's reproduce line for the rung-2 certificate passes no
   `--full`.  The reduction is legitimate for this cover only if it really is invariant under
   `x ↦ 4-x, y ↦ 4-y` (§4.3 says it is, and that it is *not* invariant under `x ↔ y`), and the
   reduction to `u ∈ [0,1/2]`, i.e. `θ ≤ 53°`, needs the `x`-reflection argument — but `u ∈ [0,½]`
   covers `θ ∈ [0°,53.1°]`, and the reflection maps `θ ↦ 90° - θ`, so the pair covers `[0°,90°]`
   with an overlap, which is fine.  None of this is stated for the rung-2 run, and a reader cannot
   tell from the note whether the published `16 872`-box census is over the full domain or half of
   it.  **`zmcheck` always runs the unreduced domain** (`c_x, c_y ∈ [0,4]`, `u ∈ [0,1]`) and checks
   no symmetry, so this question does not arise for the verification reported here — which is also
   why my box count is not comparable with the note's.
5. **`ZEROMARGIN.md` §2 and `RUNG2.md` §4.6 contradict each other about the bin clip.**  §2 says
   the centre rectangle is clipped to `[w_lo/2, m - w_lo/2]` with `w_lo` the bin *minimum*, and
   that the poses this over-tests are "sound and, at positive margin, harmless".  §4.6 then shows
   that this is a completeness bug so severe that no wall box is certifiable at any depth, and
   fixes it with `clip_bin`.  §2 was never updated.  (My §2.1 shows the issue disappears if
   admissibility is carried as a hypothesis rather than as a pre-clip.)
6. **`RUNG2.md` §3.3 (Lemma C) overloads `b_j`.**  The `b_j` in the statement are the Taylor
   coefficients of `G(u₀ + h t)` in `t` (`b_j = h^j G^{(j)}(u₀)/j!`), not Bernstein coefficients;
   the Bernstein coefficients are the `β_i`.  The formula is right; the name is the usual one for
   the other object.  Also "the overestimate is `O(h²)`" is an interior statement: the bound is
   *exact* when the maximum is at an endpoint, including the double-root-at-an-endpoint case that
   the corner region depends on.
7. **`RUNG2.md` §2.2 Lemma 1(a)** is unparseable as written ("`σ_x d_x < ½` **or** `d_x = σ_x ½`
   **and** the edge condition below holds — precisely: `p_x ∈ (i, i+1]` if …").  The clause after
   "precisely" is the statement; the two before it are a mis-bracketed paraphrase.  Not
   load-bearing for a checker (Theorem 1 is an obstruction, not a primitive), but it is the
   statement everything in §2.3 is proved from.
8. **§4's subsections are out of order in the file** (4.1, 4.2, 4.3, 4.4, 4.7, 4.5, 4.6 …), so
   "§4.6 (`clip_bin`)" is not where a reader looks for it.
9. **`RUNG2.md` §3.1's parenthesis** — "If `A_x > B_x` the admissible set at that `θ` is empty and
   the conclusion is vacuous but the test is still sound" — is true but is exactly the trap of
   §4.6: sound and useless.  Worth a forward reference.

I did read `search/zeromargin.py` in one place only: to confirm that the `pose` subcommand's output
format is `EXACT captured weight = n/d` and that `cert`'s verdict words are `VERIFIED` /
`NOT VERIFIED`, so that this checker's output could be compared line-for-line and its rejection
tests written in the same style.  No formula, no primitive and no search rule was taken from it;
the four exact pose values in §0 were taken from `RUNG2.md`, not from running the Python.

---

## 5. Rejection tests

`tests/rung2/rejection_tests.sh` (in the style of `tests/rejection_tests.sh`): **23 passed,
0 failed, 0 panics**, about 4 minutes.

| test | verdict |
|---|---|
| shipped cover at `(7/2,7/2,0)` | `1050006300/1000000000` ≥ 1 |
| baseline band `c_x ∈ [0.5,0.6]`, `c_y ∈ [1,2]` | **0** uncertified (`ADM 236 / DISJ 416 / EMPTY 292`) |
| (a) one point's weight `× 0.99` | 0 uncertified — *still a valid cover*, see below |
| (a′) **every** weight `× 19/20` (total 12.308) | **VIOLATION**, exact `997505942/1000000000` at `(1/2,1/2,0)` |
| (b) one point deleted | 0 uncertified, but 1 272 boxes against the baseline's 944 |
| (b′) every weight halved (total 6.478) | **VIOLATION** at `(1/2,1/2,0)` |
| (c) point set scaled by `0.995` (`X → 995X`, `D → 10⁶`) | **VIOLATION**: captures **exactly `0`** at `(7/2,7/2,0)` |
| (d) one point moved by `0.01` | 0 uncertified — still a valid cover |
| (d′) the whole set translated by `0.01` | **VIOLATION** at `(1/2,1/2,0)` |
| (e) `closed4_best_x103.txt` | **VIOLATION**, exact `970282351/1000000000` at `(3/2, 1461/2000, u = 1/40000)` — `RUNG2.md` §4.3's number |
| (e′) `closed4_best.txt`, same pose | **VIOLATION**, exact `9420217/10000000` |
| 9 malformed inputs (negative weight, point outside the container, truncated, `s_den ∤ s_num·D`, trailing data, `W = 0`, non-integer header, missing file, `u ∉ [0,1]`) | `ERROR`, exit 2, no verdict word, **no panics** |

Test (c) deserves a note: the shipped cover puts *all* of its corner-tile weight on the lines
`x = 3` and `y = 3` — no point of the 3 621 has both coordinates above `3.015` — so shrinking the
set by half a per cent leaves the corner square `[3,4]²` with **nothing at all** inside it.  The
exact capture at `(7/2,7/2,0)` is `0`, which is as decisive as a rejection gets.

**Mutations (a), (b) and (d) are not refused, and must not be.**  This is a finding, not a gap.
The shipped cover has about `3.1 %` of *capture* margin (`RUNG2.md` §4.7: the minimum captured
weight over a dense float scan is `1.03138`), its heaviest point weighs `0.0331` and a typical one
`0.003–0.010`.  Each of the brief's mutations (a) `−1 %` of one weight (`8.7 × 10⁻⁵`), (b) one
point deleted (`≤ 0.033`) and (d) one point moved by `0.01` therefore **leaves the file a valid
cover**, and a checker that "refused" it would be refusing a true statement — which is exactly the
failure mode a cross-check is supposed to rule out.  What such a mutation *can* break is the
certificate, since a box certificate is a fixed witness subset per region and those have far less
slack than the pointwise capture; here it does not, though deleting the point costs 35 % more boxes
in the band (1 272 against 944), which is the visible cost.  The tests therefore assert the
observed and correct outcome for (a), (b), (d) and add a *primed* variant of each that pushes the
mutation past the `3.1 %` margin, where the checker returns an exact **disproof** — a rational pose
and an exact captured weight `< 1` — rather than a mere refusal.  (a′), (b′), (c), (d′) and (e) are
all disproofs.

### 5.1 Theorem 1, confirmed by switching the disjunction off

`--nodisj` restricts the checker to the monotone primitives of §2.1–2.3 — `EMPTY`, `CORE`, `P1`,
`ADM`, i.e. exactly the *monotone witness certificates* of `RUNG2.md` Theorem 1.  It must then
fail, since the cover weighs `12.956 < m² = 16`, and it does:

```
zmcheck cert certificates/rung2/s13_closed_cover_4.txt --nodisj --depth 12 \
        --xlo 0.5 --xhi 0.6 --ylo 1.0 --yhi 2.0
NOT VERIFIED: 3473 uncertified boxes (partial sweep)
```

3 473 uncertified boxes at depth 12 in a band that the full checker discharges with 944 boxes and
0 uncertified.  This is an independent empirical confirmation of Theorem 1 — the theorem is proved
on paper in `RUNG2.md` §2 and I did not re-prove it, but this is what it predicts, measured by a
checker that knows nothing about `PIN` sets.

---

## 6. Reproduce

```bash
cd verify2 && cargo build --release && cd ..
Z=verify2/target/release/zmcheck

# the verification: full admissible domain, 8 threads
$Z cert certificates/rung2/s13_closed_cover_4.txt --depth 18 --threads 8 --dump runs/xcheck_leaves.txt

# the four exact single-pose cross-checks against RUNG2.md (instant, no boxes)
I=/home/evand/math/square-packing/s12/runs/inputs-2026-09-11
$Z pose $I/closed4_best_x103.txt --x 3/2 --y 1461/2000 --u 1/40000   # 970282351/1000000000
$Z pose $I/closed4_best.txt      --x 3/2 --y 1461/2000 --u 1/40000   # 9420217/10000000
$Z pose $I/closed4_best_x103.txt --x 7/2 --y 3/2 --u 0               # 1.562914790
$Z pose $I/closed4_best_x103.txt --x 7/2 --y 481/320 --u 0           # 1.036922939

# one box, verbose (the worst monotone-witness pose of RUNG2.md sec 2 is inside this one)
$Z box certificates/rung2/s13_closed_cover_4.txt --box "1/2,51/100,149/100,3/2,0,1/128"

# the rejection tests
sh tests/rung2/rejection_tests.sh

# a band only (never prints VERIFIED); useful for bisecting a residue
$Z cert certificates/rung2/s13_closed_cover_4.txt --xlo 0.5 --xhi 0.6 --ylo 1.0 --yhi 2.0

# Theorem 1: with the disjunction switched off the same band cannot be closed
$Z cert certificates/rung2/s13_closed_cover_4.txt --nodisj --depth 12 \
   --xlo 0.5 --xhi 0.6 --ylo 1.0 --yhi 2.0        # NOT VERIFIED: 3473 uncertified boxes

# the interior tile pose, the box that forced the Lemma-K closure of sec 3.3
$Z box certificates/rung2/s13_closed_cover_4.txt \
   --box "1500/1000,1501/1000,1499/1000,1500/1000,0,1/1024"   # DISJ, 386 regions, 1.0 s
```

Tuning knobs, none of which can make the check weaker (they only change how hard it tries):
`--depth`, `--threads`, `--nodisj` (monotone primitives only — must *not* verify, by
`RUNG2.md` Theorem 1), `--theta-bias`, `--sign-depth`, `--node-cap`, `--branch-cap`,
`--try-branches`, `--rank-exact`, `--dump`.
