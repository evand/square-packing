# Lean statements for the zero-margin checker's primitives (`lean/Sqpack/ZeroMargin.lean`)

**What this is.**  `search/zeromargin.py` proves `s(13) = 4` from
`certificates/rung2/s13_closed_cover_4.txt` (`search/RUNG2.md` §0), and it is the base of V1, the
`t = 4` verifier the `s(12)` endgame needs.  Its correctness rests on the paper lemmas of
`search/RUNG2.md` §3.1–3.3, §4.6 and §6.2.  This file formalises the **soundness of each
primitive**, in the style of `lean/Sqpack/Basic.lean` and `lean/Sqpack/Chord.lean`: hypotheses that
are exactly what the Python tests, conclusions that are the geometric facts the checker claims.
Mathlib `v4.33.1`, `lake build` 0 errors, no `sorry`, and every theorem listed in `lean/Axioms.lean`
depends only on `propext, Classical.choice, Quot.sound`.

Not formalised (and out of scope by the task brief): the subdivision, the control flow, the
weight bookkeeping, the certificate parser, the float pre-filters.  The primitives are stated over
`ℝ`; the Python evaluates the same expressions in exact `Fraction` arithmetic, so a Python verdict
is an instance of the Lean statement with rational data.

Notation as in `RUNG2.md` §1: container `[0,m]²` is `box m` (`Chord.lean`); a pose is `(c, θ)`;
`X = d_x cos θ + d_y sin θ`, `Y = −d_x sin θ + d_y cos θ` with `d = p − c`, which is literally
`coord c θ p` of `Basic.lean`, so `p ∈ sq c θ 1 ↔ |X| ≤ ½ ∧ |Y| ≤ ½` holds *by definition* — no
bridging lemma is needed between the note's `Q(c,θ)` and the Lean `sq c θ 1`.
`u = tan(θ/2)`, i.e. `θ = 2 arctan u`.

---

## 0.  Admissibility — the support-function fact (`RUNG2.md` §1, §4.6)

| Lean | statement | implementation |
|---|---|---|
| `wid θ` | `|cos θ| + |sin θ|` | `w(θ)`; `bin_data`'s `w0/w1` |
| `Adm m c θ` | `w/2 ≤ c_x ≤ m − w/2` and same in `y` | `clip_bin` docstring, `_adm_specs` |
| **`sq_subset_box_iff`** | `sq c θ 1 ⊆ box m ↔ Adm m c θ` | the definition of *admissible* the whole checker uses |

```lean
theorem sq_subset_box_iff (m : ℝ) (c : ℝ × ℝ) (θ : ℝ) :
    sq c θ 1 ⊆ box m ↔ Adm m c θ
```

Both directions.  `⟸` is the triangle inequality on `p − c = (X cos θ − Y sin θ, X sin θ + Y cos θ)`
(`abs_sub_le_wid`); `⟹` exhibits the four extreme vertices (`mem_sq_offset`, with the signs of
`cos θ, sin θ`).  `Chord.lean`'s `half_height_le_centre` is the `c_y ≥ w/2` half of this; the new
theorem subsumes it (it is *not* re-proved from it — `Chord.lean` is untouched).

This is the fact the brief called "the support-function fact for `sq ⊆ [0,m]²`"; `Basic.lean` did
not in fact have it.

---

## 1.  Lemma A — monotone corners (`RUNG2.md` §3.1 ↔ `_adm_cond_ok`, `_adm_exact`)

Four single-inequality lemmas, one per condition, each with **its own** corner:

| Lean | condition | corner | Python `cond` |
|---|---|---|---|
| `coord_fst_le_of_corner` | `X ≤ ½` | `(A_x, A_y)` | `0` |
| `neg_half_le_coord_fst` | `X ≥ −½` | `(B_x, B_y)` | `1` |
| `coord_snd_le_of_corner` | `Y ≤ ½` | `(B_x, A_y)` | `2` |
| `neg_half_le_coord_snd` | `Y ≥ −½` | `(A_x, B_y)` | `3` |

matching the comment in `_adm_cond_ok`
(`X<=1/2 at (Ax,Ay); X>=-1/2 at (Bx,By); Y<=1/2 at (Bx,Ay); Y>=-1/2 at (Ax,By)`).  Each needs
`cos θ ≥ 0`, `sin θ ≥ 0` and only the two relevant monotonicity bounds on `c`.

```lean
theorem mem_sq_of_corners {p c : ℝ × ℝ} {θ : ℝ} {Ax₁ Ay₁ Bx₂ By₂ Bx₃ Ay₃ Ax₄ By₄ : ℝ}
    (hcos : 0 ≤ Real.cos θ) (hsin : 0 ≤ Real.sin θ)
    (h₁ : (p.1 - Ax₁) * Real.cos θ + (p.2 - Ay₁) * Real.sin θ ≤ 1 / 2)
    (h₂ : -(1 / 2) ≤ (p.1 - Bx₂) * Real.cos θ + (p.2 - By₂) * Real.sin θ)
    (h₃ : -(p.1 - Bx₃) * Real.sin θ + (p.2 - Ay₃) * Real.cos θ ≤ 1 / 2)
    (h₄ : -(1 / 2) ≤ -(p.1 - Ax₄) * Real.sin θ + (p.2 - By₄) * Real.cos θ)
    (hx₁ : Ax₁ ≤ c.1) (hy₁ : Ay₁ ≤ c.2) (hx₂ : c.1 ≤ Bx₂) (hy₂ : c.2 ≤ By₂)
    (hx₃ : c.1 ≤ Bx₃) (hy₃ : Ay₃ ≤ c.2) (hx₄ : Ax₄ ≤ c.1) (hy₄ : c.2 ≤ By₄) :
    p ∈ sq c θ 1
```

The **eight independent corner values** are deliberate: they are what makes the Python's
per-condition choice of bound kind (`'R'` = the box side, `'W'`/`'M'` = the container wall) sound.
`_adm_specs` offers both kinds whenever the wall can bind on the bin, and `_adm_cond_ok` accepts
condition `cond` if *any* offered `(x-kind, y-kind)` pair certifies it — possibly a different pair
for each of the four conditions.  `mem_sq_of_corners` is exactly that hypothesis.

And Lemma A itself, over a whole bin, with the note's `A_x = max(cx₀, w/2)`,
`B_x = min(cx₁, m − w/2)` (`admLo`, `admHi`):

```lean
theorem lemmaA (m : ℝ) (p : ℝ × ℝ) (cx₀ cx₁ cy₀ cy₁ : ℝ) (Bin : Set ℝ)
    (hBin : ∀ θ ∈ Bin, 0 ≤ Real.cos θ ∧ 0 ≤ Real.sin θ)
    (h₁ : ∀ θ ∈ Bin, (p.1 - admLo cx₀ θ) * Real.cos θ
            + (p.2 - admLo cy₀ θ) * Real.sin θ ≤ 1 / 2)
    (h₂ : ∀ θ ∈ Bin, -(1 / 2) ≤ (p.1 - admHi m cx₁ θ) * Real.cos θ
            + (p.2 - admHi m cy₁ θ) * Real.sin θ)
    (h₃ : ∀ θ ∈ Bin, -(p.1 - admHi m cx₁ θ) * Real.sin θ
            + (p.2 - admLo cy₀ θ) * Real.cos θ ≤ 1 / 2)
    (h₄ : ∀ θ ∈ Bin, -(1 / 2) ≤ -(p.1 - admLo cx₀ θ) * Real.sin θ
            + (p.2 - admHi m cy₁ θ) * Real.cos θ) :
    ∀ (c : ℝ × ℝ) (θ : ℝ), θ ∈ Bin → c.1 ∈ Set.Icc cx₀ cx₁ → c.2 ∈ Set.Icc cy₀ cy₁ →
      sq c θ 1 ⊆ box m → p ∈ sq c θ 1
```

* The bin is an arbitrary set of angles with `cos, sin ≥ 0`, which is the note's
  `u₀ ≥ 0, u₁ ≤ 1`.  The shipped root grid (`roots`, `ubins = 8`) uses `u ∈ [0, 1/2]`, well inside.
* No hypothesis `A_x ≤ B_x` is needed and none is available: at an angle where the admissible set
  is empty the conclusion is vacuous but the hypotheses are still sound.  That is exactly the
  conservatism §4.6 diagnoses and `clip_bin` removes (see §4 below).

---

## 2.  Lemma B — the bounds are polynomials in `u` (`RUNG2.md` §3.2 ↔ `_xn`, `_cond_poly`)

| Lean | implementation |
|---|---|
| `qeval a u` / `qeval4 a u` | a quadratic / a degree-4 polynomial from a coefficient tuple |
| `xnR a = (2a, 0, 2a)` | `_xn('R', a)` |
| `xnW = (1, 2, -1)` | `_xn('W')` |
| `xnM m = (2m-1, -2, 2m+1)` | `_xn('M')` |
| `xnR_spec`, `xnW_spec`, `xnM_spec` | the three bounds really are `a`, `w(θ)/2`, `m − w(θ)/2` |
| `ucoef px Xn`, `ucoef_spec` | `U = 2(1+u²)p_x − X_n` |
| `condPoly k U V` | `_cond_poly(U, V, cond, both_rect=False)` |
| `condPoly_zero/one/two/three` | the coefficients really are those of `G₁ = U·C + V·S − N²`, …, `G₄ = U·S − V·C − N²` |

The load-bearing statement:

```lean
theorem condPoly_eq_gval (k : Fin 4) (px py : ℝ) (Xn Yn : ℝ × ℝ × ℝ) (u : ℝ) :
    qeval4 (condPoly k (ucoef px Xn) (ucoef py Yn)) u
      = (1 + u ^ 2) * gval k (px - qeval Xn u / (2 * (1 + u ^ 2)))
          (py - qeval Yn u / (2 * (1 + u ^ 2))) u
```

with `condPoly_nonpos_iff` the corresponding `≤ 0 ↔ ≤ 0`.  Read together with
`mem_sq_iff_gval` (§3) this says: `_cond_poly`'s degree-4 polynomial is `≤ 0` at `u` **iff** the
corresponding Lemma A inequality holds at the corner given by those two bounds at that `θ`.  It is
therefore a machine check on the five hand-derived coefficient formulas in `_cond_poly` — the most
error-prone lines of the checker.

`condPoly_rect` is the note's factorisation: with both bounds of kind `'R'`, `U = U₀N`, `V = V₀N`,
and dividing out `N` leaves the quadratic `G_{p,k}` itself — so `_cond_poly`'s `both_rect` branch
returns precisely `gc k (p_x − a) (p_y − b)` padded with two zeros, which is why `_poly_ok` then
takes the exact `_max_quad` route.

---

## 3.  The violation polynomials, Lemma C, Lemma E (`RUNG2.md` §3.3, §6.2 ↔ `_gcoef`, `_max_quad`, `_max_bern`, `_gmax`)

### The violation polynomials (`_gcoef`)

```lean
def gc : Fin 4 → ℝ → ℝ → ℝ × ℝ × ℝ
  | 0, a, b => (2 * a - 1, 4 * b, -2 * a - 1)      -- G = 2aC + 2bS - N
  | 1, a, b => (-2 * a - 1, -4 * b, 2 * a - 1)
  | 2, a, b => (2 * b - 1, -4 * a, -2 * b - 1)
  | 3, a, b => (-2 * b - 1, 4 * a, 2 * b - 1)

theorem mem_sq_iff_gval (c p : ℝ × ℝ) (u : ℝ) :
    p ∈ sq c (2 * Real.arctan u) 1 ↔ ∀ k : Fin 4, gval k (p.1 - c.1) (p.2 - c.2) u ≤ 0
```

This is the Lean counterpart of `zeromargin_stress.py`'s "200,000 random `(point, pose)` pairs,
comparing `max_k G_{p,k} ≤ 0` with the geometric `p ∈ Q(c,θ)` — 0 disagreements", now for all
poses at once.  It needs `cos (2 arctan u) = (1−u²)/(1+u²)`, `sin (2 arctan u) = 2u/(1+u²)`
(`cos_two_arctan`, `sin_two_arctan`), which is `trig(u)` in the Python.

### Lemma C, degree `≤ 2` (`_max_quad`)

```lean
noncomputable def maxQuad (a₀ a₁ a₂ u₀ u₁ : ℝ) : ℝ :=
  if a₂ < 0 ∧ u₀ < -a₁ / (2 * a₂) ∧ -a₁ / (2 * a₂) < u₁ then
    max (max (a₀ + a₁ * u₀ + a₂ * u₀ ^ 2) (a₀ + a₁ * u₁ + a₂ * u₁ ^ 2)) (a₀ - a₁ ^ 2 / (4 * a₂))
  else max (a₀ + a₁ * u₀ + a₂ * u₀ ^ 2) (a₀ + a₁ * u₁ + a₂ * u₁ ^ 2)

theorem le_maxQuad {a₀ a₁ a₂ u₀ u₁ u : ℝ} (hu : u ∈ Set.Icc u₀ u₁) :
    a₀ + a₁ * u + a₂ * u ^ 2 ≤ maxQuad a₀ a₁ a₂ u₀ u₁
theorem maxQuad_mem {a₀ a₁ a₂ u₀ u₁ : ℝ} (h : u₀ ≤ u₁) :
    ∃ u ∈ Set.Icc u₀ u₁, a₀ + a₁ * u + a₂ * u ^ 2 = maxQuad a₀ a₁ a₂ u₀ u₁
```

Sound **and exact** (`IsGreatest` in two halves), which is the note's "no Bernstein slack".  The
definition mirrors `_max_quad` line for line, including the strict `u0 < uv < u1` guard; the value
the Python computes at the vertex, `a₀ + a₁u_v + a₂u_v²`, is `a₀ − a₁²/(4a₂)`.

### Lemma C, degree 3–4 (`_max_bern`)

| Lean | implementation |
|---|---|
| `bshift a u₀ h` | the `b_j = Σ_{k≥j} a_k C(k,j) u₀^{k−j} h^j` loop |
| `bernCoef b` | the `β_i = Σ_{j≤i} [C(i,j)/C(4,j)] b_j` loop |
| `maxBern a u₀ u₁` | the returned `best` |
| `qeval4_bshift` | `G(u₀ + h t) = Σ_j b_j t^j` — a machine check on the `b_j` formula |
| `qeval4_eq_bern` | `Σ_j b_j t^j = Σ_i β_i B_{i,4}(t)` — a machine check on the `β_i` formula |
| `bern_le_max` | the convex-hull property in degree 4: `B_{i,4} ≥ 0`, `Σ B_{i,4} = 1` |
| **`qeval4_le_maxBern`** | `∀ u ∈ [u₀,u₁], G(u) ≤ _max_bern(a, u₀, u₁)` |
| `bern_endpoints` | `β₀ = G(u₀)`, `β₄ = G(u₁)` |
| `polyOk_bern`, `polyOk_quad` | the two branches of `_poly_ok` |

Because the degree is fixed at 4, the two basis-change identities are explicit polynomial
identities closed by `ring`; no general Bernstein theory (and no Mathlib `bernsteinPolynomial`) is
used.  The degenerate bin `u₀ = u₁` is handled separately (`h = 0`, all `β_i = G(u₀)`).

### Lemma E (`_gmax`)

```lean
noncomputable def gsum (lam : Fin n → ℝ) (pt : Fin n → ℝ × ℝ) (kind : Fin n → Fin 4)
    (cx cy u : ℝ) : ℝ :=
  ∑ r : Fin n, lam r * gval (kind r) ((pt r).1 - cx) ((pt r).2 - cy) u

noncomputable def gmax (lam pt kind) (cx₀ cx₁ cy₀ cy₁ u₀ u₁ : ℝ) : ℝ :=   -- the four corners
  max (max (gmaxAt … cx₀ cy₀ u₀ u₁) (gmaxAt … cx₀ cy₁ u₀ u₁))
    (max (gmaxAt … cx₁ cy₀ u₀ u₁) (gmaxAt … cx₁ cy₁ u₀ u₁))              -- gmaxAt = maxQuad ∘ gsc

theorem gsum_le_gmax (lam pt kind) {cx₀ cx₁ cy₀ cy₁ u₀ u₁ cx cy u : ℝ}
    (hcx : cx ∈ Set.Icc cx₀ cx₁) (hcy : cy ∈ Set.Icc cy₀ cy₁) (hu : u ∈ Set.Icc u₀ u₁) :
    gsum lam pt kind cx cy u ≤ gmax lam pt kind cx₀ cx₁ cy₀ cy₁ u₀ u₁

theorem gmax_mem (lam pt kind) {cx₀ cx₁ cy₀ cy₁ u₀ u₁ : ℝ}
    (hx : cx₀ ≤ cx₁) (hy : cy₀ ≤ cy₁) (hu : u₀ ≤ u₁) :
    ∃ cx ∈ Set.Icc cx₀ cx₁, ∃ cy ∈ Set.Icc cy₀ cy₁, ∃ u ∈ Set.Icc u₀ u₁,
      gsum lam pt kind cx cy u = gmax lam pt kind cx₀ cx₁ cy₀ cy₁ u₀ u₁
```

`gmax` is `_gmax` exactly: the four centre corners of `_gmax`'s double loop, and `_max_quad` at
each.  `gsc` is the `c0, c1, c2` accumulator.  Soundness (`gsum_le_gmax`) is `gsum_le_corners`
(affinity in `(c_x,c_y)`, via `gsum_affine_x`/`gsum_affine_y` and `affine_le_max`) composed with
`le_maxQuad`; exactness (`gmax_mem`) is `maxQuad_mem` at the winning corner.  This is what
`zeromargin_stress.py` checks numerically ("the `_gmax` enclosure (Lemma E): 20,000 random
`(box, nonnegative combination)` pairs … 0 violations").

`gsum_le_gmax` is stated for the *box*; since the admissible poses of the box are a subset, it is a
fortiori sound for them (`_gmax`'s docstring says the same).

---

## 4.  Lemmas F, G, H — `CHAIN` (`RUNG2.md` §6.2 ↔ `cert_chain`)

Scalar forms (the one-line arguments, which is all the lemmas are):

```lean
theorem lemmaG  {ga gq lam : ℝ} (hlam : 0 < lam) (h : ga + lam * gq ≤ 0) (hq : 0 < gq) : ga < 0
theorem lemmaG' {ga gq lam : ℝ} (hlam : 0 ≤ lam) (h : ga + lam * gq ≤ 0) (hq : 0 ≤ gq) : ga ≤ 0
theorem lemmaH  {gq gq' lam : ℝ} (hlam : 0 < lam) (h : gq + lam * gq' ≤ 0) : ¬(0 < gq ∧ 0 < gq')
theorem chain_region_down {gp gq : ℝ} (h : gp ≤ gq) (hq : gq ≤ 0) : gp ≤ 0
```

and the box forms, whose hypotheses are literally the `_gmax` calls `cert_chain` makes
(`gmax1(k1,d1,k2,d2,lam)` = `_gmax [(1,k1,d1), (lam,k2,d2)] box`):

```lean
theorem lemmaG_box (lam : ℝ) (hlam : 0 < lam) (pa pq : ℝ × ℝ) (ka kq : Fin 4)
    (hmax : gmax ![1, lam] ![pa, pq] ![ka, kq] cx₀ cx₁ cy₀ cy₁ u₀ u₁ ≤ 0)
    (hcx …) (hcy …) (hu …) (hq : 0 < gval kq (pq.1 - cx) (pq.2 - cy) u) :
    gval ka (pa.1 - cx) (pa.2 - cy) u < 0

theorem lemmaH_box (lam : ℝ) (hlam : 0 < lam) (pq pq' : ℝ × ℝ) (kq kq' : Fin 4)
    (hmax : gmax ![1, lam] ![pq, pq'] ![kq, kq'] cx₀ cx₁ cy₀ cy₁ u₀ u₁ ≤ 0) (hcx …) (hcy …) (hu …) :
    ¬(0 < gval kq  (pq.1  - cx) (pq.2  - cy) u ∧
      0 < gval kq' (pq'.1 - cx) (pq'.2 - cy) u)

theorem gval_le_of_gmax_sub (pp pq : ℝ × ℝ) (kp kq : Fin 4)
    (hmax : gmax ![1, -1] ![pp, pq] ![kp, kq] cx₀ cx₁ cy₀ cy₁ u₀ u₁ ≤ 0) (hcx …) (hcy …) (hu …) :
    gval kp (pp.1 - cx) (pp.2 - cy) u ≤ gval kq (pq.1 - cx) (pq.2 - cy) u
```

`gval_le_of_gmax_sub` is the chain-building and `down`-search test
(`_gmax [(1, k, kind), (-1, kq, kdq)] box ≤ 0`).

Lemma F, the covering half:

```lean
theorem chain_regions_cover (k : ℕ) (g : ℕ → ℝ) (hmono : ∀ i j, i ≤ j → g i ≤ g j) :
    ∃ r ≤ k, (∀ j, 1 ≤ j → j ≤ r → g j ≤ 0) ∧ (∀ j, r < j → j ≤ k → 0 < g j)
```

`g j` is the value at a given pose of the `j`-th chain member's violation polynomial; `hmono` is
what `gval_le_of_gmax_sub` supplies for consecutive pairs plus transitivity (the note's
"`G_{q_1} ≤ … ≤ G_{q_k}` on `B`, by transitivity").  The conclusion is that the `k+1` regions
`R₀ = {g₁ > 0}`, `R_r = {g_r ≤ 0 < g_{r+1}}`, `R_k = {g_k ≤ 0}` cover the box and are disjoint, and
`chain_region_down` is why the whole down-set (not just `{q_1,…,q_r}`) is available on `R_r`.

Not formalised: the weight arithmetic of step 4 of the primitive
(`w(T) + w({q_1..q_r}) + w(U_{r+1}) ≥ 1` for every region), the staircase/binary-search
structure, and the product-of-two-chains bookkeeping.  Those are control flow, and a Lean statement
of them would have to model the whole witness-set data structure.

---

## 5.  `clip_bin` (`RUNG2.md` §4.6 ↔ `clip_bin`)

```lean
noncomputable def widU (u : ℝ) : ℝ := (1 - u ^ 2 + 2 * u) / (1 + u ^ 2)   -- w as a function of u

lemma wid_two_arctan {u : ℝ} (h0 : 0 ≤ u) (h1 : u ≤ 1) : wid (2 * Real.arctan u) = widU u

theorem widU_strictMonoOn : StrictMonoOn widU (Set.Icc 0 (Real.sqrt 2 - 1))

theorem adm_widU_le (m : ℝ) {cx₀ cx₁ cy₀ cy₁ : ℝ} {c : ℝ × ℝ} {u : ℝ}
    (h0 : 0 ≤ u) (h1 : u ≤ 1)
    (hcx : c.1 ∈ Set.Icc cx₀ cx₁) (hcy : c.2 ∈ Set.Icc cy₀ cy₁)
    (hsub : sq c (2 * Real.arctan u) 1 ⊆ box m) :
    widU u ≤ 2 * min (min cx₁ (m - cx₀)) (min cy₁ (m - cy₀))

theorem clip_bin_no_loss (m : ℝ) {cx₀ cx₁ cy₀ cy₁ u₀ u₁ ustar : ℝ}
    (h0 : 0 ≤ u₀) (h45 : u₁ ≤ Real.sqrt 2 - 1) (hstar : ustar ∈ Set.Icc u₀ u₁)
    (hK : 2 * min (min cx₁ (m - cx₀)) (min cy₁ (m - cy₀)) ≤ widU ustar)
    {c : ℝ × ℝ} {u : ℝ} (hu : u ∈ Set.Icc u₀ u₁)
    (hcx : c.1 ∈ Set.Icc cx₀ cx₁) (hcy : c.2 ∈ Set.Icc cy₀ cy₁)
    (hsub : sq c (2 * Real.arctan u) 1 ⊆ box m) :
    u ∈ Set.Icc u₀ ustar
```

* `adm_widU_le` is the note's "a pose of the box is admissible iff `w(θ) ≤ K` with
  `K = 2 min(cx₁, m − cx₀, cy₁, m − cy₀)`" — the `⟹` half, which is the half the clip needs.
* `h45 : u₁ ≤ √2 − 1` is the Python's guard `if u1*u1 + 2*u1 - 1 > 0: return u1` (the root of
  `u² + 2u − 1` is `√2 − 1`, i.e. `45°`); `sqrt_two_sub_one_sq` records `s² + 2s = 1`.
* `hK : K ≤ widU u*` is the invariant `clip_bin` maintains on the value it returns, in all three
  branches where it actually narrows the bin:
  * `if w(lo) > K: return u0` — then `widU u₀ > K`;
  * `if w(lo) == K: return lo` — then `widU u₀ = K`;
  * the bisection's `return lo if w(lo) == K else hi` — `lo` satisfies `w(lo) = K` or `hi`
    satisfies the loop invariant `w(hi) > K`.

  The fourth branch, `if w(u1) <= K: return u1`, returns the bin unchanged; there the conclusion
  `u ∈ [u₀, u*] = [u₀, u₁]` is the hypothesis, so nothing is needed.  So the theorem covers the
  clip exactly: it says the returned `u*` drops no admissible pose, which is the soundness claim
  ("from above, so no admissible pose is ever dropped").
* The bisection's *termination and accuracy* are not formalised, and do not need to be:
  `clip_bin_no_loss` holds for **any** `u*` in the bin with `widU u* ≥ K`, however it was found.

---

## 6.  Discrepancies found

1. **The task brief's Lemma G is stated backwards, and as stated is false.**  The brief
   (`tasks/lean-zeromargin/README.md`, target 3) says

   > if `max_B (g_i + λ G) ≤ 0` for some `λ ≥ 0` then `g_i ≤ 0` on `B ∩ {G ≤ 0}`.

   From `g_i + λG ≤ 0` and `G ≤ 0` one gets `g_i ≤ −λG ≥ 0` — no conclusion.  `RUNG2.md` §6.2's
   Lemma G has the correct sign:

   > If `max_B (G_a + λ G_q) ≤ 0` for some `λ > 0`, then at every pose of `B` with `G_q > 0` one
   > has `G_a < 0`.

   and this is the version `cert_chain` uses (the `up` sets certify a point `a` on the region where
   the *next* pivot's `G_{q_{r+1}}` is **violated**, i.e. `> 0`).  Formalised as the note has it
   (`lemmaG`), with the `λ ≥ 0 / G_q ≥ 0 / G_a ≤ 0` variant as `lemmaG'`.

2. **Lemma E's "nonnegative combination" is narrower than the code's use.**  §6.2 states Lemma E
   for `Σ λ_r G_{p_r,k_r}` with `λ ≥ 0`, but `cert_chain` also calls `_gmax` with the combination
   `(1, p, kp) + (−1, q, kq)` — a *negative* multiplier — for the chain comparison and the `down`
   binary search.  Nothing breaks: affinity in `(c_x,c_y)` and quadraticity in `u` are independent
   of the signs, so `gsum_le_gmax` and `gmax_mem` are proved with `lam : Fin n → ℝ` unconstrained,
   and `gval_le_of_gmax_sub` uses that.  The note's "nonnegative" should read "real".

3. **Lemma A as stated in §3.1 is *stronger* than the test `_adm_cond_ok` performs.**  The note
   evaluates the inequalities at `A_x = max(cx₀, w(θ)/2)`; the Python never forms that maximum
   (it is not polynomial in `u`) and instead tries `cx₀` and `w(θ)/2` separately, accepting if
   either works — and may use different choices for different conditions.  Each individual choice
   is a *weaker* (more conservative) test than the note's, since lowering `A_x` enlarges the left
   side of (i).  §3.2 says as much in prose; the Lean file makes it a hypothesis by giving
   `mem_sq_of_corners` eight independent corner values, with `lemmaA` as the note's form on top.

4. **`bin_data`'s `whi` is a rational over-estimate**, `F(14143,10000) > √2`, when the bin straddles
   `45°`.  It feeds only `_adm_specs` (which bounds are *offered* — offering a spare bound is
   harmless) and `cert_chain`'s `t = 1 − whi/2` for `P1` (a smaller `t`, hence a stricter test).
   Sound in both places, but it means the Lean `admLo`/`admHi` with the exact `wid θ` are the
   *ideal*, not what the float/rational code evaluates; `mem_sq_of_corners` is the statement that
   covers the code.

5. `Basic.lean` did **not** already contain the support-function fact for `sq ⊆ [0,m]²` that the
   brief expected to reuse (`Chord.lean` has only the `c_y ≥ w/2` half, `half_height_le_centre`);
   `sq_subset_box_iff` supplies it, both directions and all four sides.

---

## 7.  Files

* `lean/Sqpack/ZeroMargin.lean` — 970 lines, 0 `sorry`.
* `lean/Sqpack.lean` — now imports `Sqpack.ZeroMargin`.
* `lean/Axioms.lean` — `#print axioms` for 36 theorems (the 4 pre-existing ones plus 32 new);
  all show `[propext, Classical.choice, Quot.sound]`.

Build: symlink the main checkout's Mathlib (`ln -s …/s12/lean/.lake <worktree>/lean/.lake`), then
`cd lean && lake build`.
