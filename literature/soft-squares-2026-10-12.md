# Soft-square homotopy: the two soft-rectangle papers, and a design for `fq` (2026-10-12)

**Update:** see the [addendum on Melashenko & Romanova 2025 (Ukrainian, soft polyhedra)](#addendum-2026-10-12-melashenko--romanova-ukrainian-paper) at the end.  It changes no part of the test in §3.

Follow-up to [`irregular-packing-2026-10-12.md`](irregular-packing-2026-10-12.md) §4, §5, §12(c) and its "Ranked" item 1.
Both papers were read in full this time (user-supplied PDFs in `../private/literature/`, converted with `pdftotext
-layout` in the session scratchpad; not copied here).  This note corrects §4 and §12(c) of the earlier note, which only
had abstracts.

**Labels.**  **[read]**: from the full text, with section / equation numbers.  **[verified]**: from a record or
abstract only.  **[inference]**: my own reasoning; nothing in either paper supports it.

**Bottom line.**
- Neither paper hardens back to rigid shapes.  Neither has a continuation schedule or a rigid-vs-soft comparison
  beyond one example, and neither has anything about equal squares.
- The idea is ours, and so is the risk.  The formulation carries over to `fq` almost unchanged: their separating-line
  "quasi phi-function" is our lifted separator, and a shape variable is just n more columns.
- My analysis (§2.4) says the plain homotopy absorbs a k-row jam by deforming in place, and hardening then restores it.
  It can only work if the topology changes while the squares are soft.  So the test below measures that change
  directly, against a rigid control that runs the same schedule.

---

## 1. What the papers do

### 1a. Melashenko, Romanova, Litvinchev, Martínez Gómez, Yang, Sun 2025: soft rotated rectangles [read]

*Mathematics* 13(3) 493, doi:10.3390/math13030493.

- **Shape variables and parametrisation.**
  - §2, eq. (1): R_i(u_i, μ_i) = {p = (x_i, y_i) + A(θ_i) · M(μ_i^x, μ_i^y) · p⁰ : p⁰ ∈ R_i⁰}.
    A is the rotation and M = diag(μ^x, μ^y) is the "elasticity matrix".  Stretch is applied in the body frame,
    before rotation.
  - The original rectangle is R_i⁰ = [−a_i, a_i] × [−b_i, b_i].  Sizes may differ per item.
- **Area held fixed.**
  - Proposition 1 (proof in Appendix A, via the surveyor's formula): area is preserved iff μ^x μ^y = 1.
  - They then set μ^x = μ_i, μ^y = 1/μ_i, so there is one shape variable per item.
- **Elasticity bounds.**
  - The bounds are μ⁻ ≤ μ_i ≤ μ⁺ (placement condition (5); model constraint (13), written with a strict "<" on
    the lower side).
  - They are chosen ad hoc per run, with no rationale given.  Examples include [1, 2], [0.5, 1.5], [0.3, 1.5],
    [0.3, 2.5], [0.5, 2.5], [0.8, 2.5] and [0.5, 1.8].
  - Some ranges are one-sided ([1, 2]): with θ free, μ ↔ 1/μ is a 90° rotation.  The ranges used per instance are
    only in an external ResearchGate file (§5).
- **Non-overlap.**
  - Eq. (2) is a quasi phi-function: a separating line with variables τ_ij = (φ_ij, γ_ij) per pair.  The 4 vertices
    of R_i must lie on one side of it and the 4 vertices of R_j on the other.
  - The vertices are (μ^x x̃ cos θ − μ^y ỹ sin θ, μ^x x̃ sin θ + μ^y ỹ cos θ).
  - This is exactly `fq`'s lifted separator, with corners scaled by (μ, 1/μ).
  - Containment: eq. (4) for a rectangle container (min over corners and walls); (5) circle; (6) scaled convex
    polygon.  Obstacles ("prohibited zones"): (3) and (7).
- **Rotation.**  θ_i is continuous and free.  There are no orientation constraints.
- **Model.**
  - Eqs. (8)–(13): min κ(L) (container size).
  - Dimension d = 4n + n(n−1) + 2nΣn_z + |L|.
- **Solver.**
  - Algorithm 1 builds the starting points.  Random centres; circumscribed circles grown homothetically by
    maximising ρ with IPOPT; a larger container if ρ < 1; separators derived from the circle pairs; random θ;
    μ = 1.
  - Algorithm 2 is decomposition, iterated.  Each item gets a box neighbourhood of size r_i + ε_μ with
    ε_μ = μ⁺ · max r_i (eq. 14).  Pairs whose boxes cannot meet are dropped.
  - The subproblem is (15)–(22), including κ(L) ≤ κ(L^(k)) (21), solved by IPOPT 3.14.14.  It stops when the
    improvement is below δ (δ = 0.01 in Fig. 5).
  - Python 3.11 on an i5-12500.  **Three starting points per instance**, 300 s CPU limit (§5).
- **Instances and results (§5).**
  - Examples 1–7: n = 50, total area 400, in circle, square, triangle, pentagon, hexagon, heptagon and octagon
    containers.  Fill 88–94 % in polygons and circles, **98.6–99.97 % in the square** (Ex. 2: 0.3 ≤ μ ≤ 1.5 gives
    h* = 20.003 against the area bound 20).
  - Examples 8–15: n = 30, 46 or 49, with prohibited zones.  Fill 84–98 %.
  - Each example reports three μ ranges, with fill varying by a few percent between them, non-monotonically.
- **Hardening / homotopy:** none.  Softness is part of the problem.  μ starts at 1 and is simply a variable.
- **Equal objects:** none.  Areas are given as totals (mean area 8; Fig. 3 shows one 8-area rectangle).  Whether the
  items are identical is not stated in the text.
- **Notable remark (§1).**  "For the rectangular container in the optimal arrangement, all soft rectangles are
  allocated such that their sides are parallel to the sides of the container."  This is the soft analogue of the
  grid attractor (see §2.4).

### 1b. Bennell, Litvinchev, Pankratov, Romanova 2024: stretched convex polygons [read]

*Wireless Networks* 30 (2024) 7369–7376, doi:10.1007/s11276-023-03642-9.

- **Shape variable.**
  - §2: P_i(u_i, k_i) = {v_i + M(θ_i) · A(k_i) · p⁰}, with A(k) = diag(k, 1/k) (horizontal stretch, vertical
    shrink, in the body frame).
  - The bound is **1 ≤ k_i ≤ k⁺** (eq. 6), one-sided; rotation covers the other orientation.
- **Area.**  Area is preserved because det A = 1, shown via the shoelace formula (§2).  Convexity is preserved because
  the map is affine.
- **Bounds:** one example, with k⁺ = 1.6, and no rationale.
- **Non-overlap.**  Eq. (2): the same separating-line quasi phi-function as 1a, over the m_i transformed vertices.
  Containment in a rectangle of width l and variable height h: eq. (1).  (w_i4 prints "+ l" where h is meant.)
- **Model.**  Eqs. (3)–(7): min h, d = n² + 3n + 1.
- **Solver (§4).**
  - Multistart.  Starting points come from rigid polygons scaled by β, grown from β = 0 to 1 in a fixed rectangle.
  - Then "soft decomposition": rectangular sub-containers of variable size around each polygon, far pairs dropped,
    O(n) subproblems in sequence, IPOPT.
- **Instance and result (§5): one example.**
  - n = 100 polygons; the shapes are not given in the text.  Container width 25, starting box 25 × 25.5.
  - Rigid packing h = 23.400602 (77 s); soft packing (1 ≤ k_i ≤ 1.6) h = 21.864322 (505 s).
  - 38 % of items stay at k = 1, and 62 % stretch (Table 1 lists all k_i; several sit at the bound 1.6).
- **Direction of travel:** rigid → soft (the rigid packing is the start).  This is the reverse of a hardening schedule.
- **Equal objects:** none.

### 1c. Litvinchev, Infante, Romanova, Martinez-Noa, Gutierrez 2024 [verified abstract only]

*Mobile Netw. Appl.* 29, 211–220.  Soft convex polygons with **vertex coordinates as variables** (Melashenko §6, ref
[29]), plus an "inverse problem": the least elasticity that still achieves a given optimal container.  Not read beyond
the abstract.

### Where the evidence is thin

- **There is no hardening anywhere.**  So there is no evidence either way on the one thing we want: whether softness
  followed by rigidity lands somewhere better.
- **There is no rigid baseline.**
  - Melashenko never reports μ ≡ 1 on the same instances.
  - Bennell's one comparison (rigid 23.40 vs soft 21.86) is a relaxation beating its restriction, which is
    guaranteed.
- **The sample is tiny.**
  - Three starts per instance, one run per bound setting, with no variance and no statistics.
  - The bound settings change between examples, so the effect of the bound cannot be separated out.
- **The instance data are off-paper** (ResearchGate file, "on request").  None are reproducible from the text.
- **Equal squares in a square: absent.**
  - The only square-container result (Melashenko Ex. 2) says soft rectangles **tile** a square almost perfectly
    (99.97 %).
  - Read for our purposes, that is a warning: the soft optimum is a rectangular grid (§2.4).

---

## 2. Design sketch: soft squares in `fq`

### 2.1 Variables

- **Shape.**  One log-aspect per square: w_i = e^{a_i/2}, h_i = e^{−a_i/2}, so the area w_i h_i = 1 exactly and the
  shape is symmetric in a ↔ −a.
- **Corners.**  In the body frame, then rotated (the papers' A(θ)M order): o_k = H R(θ_i) (w_i s_x, h_i s_y).  This
  is `corner_off` with per-square (w, h).
- **Derivative.**  ∂o_k/∂a_i = ½ H R(θ_i) (w_i s_x, −h_i s_y).  It is smooth, and every separator and wall constraint
  stays smooth: they are linear in the corners.
- **Variable layout.**  Append the a_i *after* the pair block (3n + 1 + 2·pairs + n) so all existing offsets are
  unchanged.
- **Bound.**  Use a_i = ε tanh(b_i) with b_i free.
  - The bound |a_i| ≤ ε is then automatic and smooth, with no extra ALM rows.
  - Scaling ε → 0 hardens every square at once, exactly.
  - Optionally add a quadratic preference for squareness, (κ/2) Σ a_i², so that softness is spent only where it is
    loaded.
- **Sign and rotation.**  Use symmetric |a| ≤ ε with θ free.  It is redundant (a ↔ −a with θ + 90°) but kink-free;
  the papers' one-sided [1, μ⁺] would put a bound kink on every square that stays square.

### 2.2 What generalises in the code (read-only survey of `search/packer/src/fq.rs`)

- **Unchanged once `corner_off` takes (w, h):** `eval` / `viol_update` (pairs and walls), `face_seps`, `init_sep`.
  Gradients add one term per corner: dg/da_i = u · ∂o/∂a_i for separators, ±∂o_x/∂a or ±∂o_y/∂a for walls.
- **`sat_gap` and `check` need the rectangle projection.**  The half-extent along a direction t becomes
  H (w |cos(θ − t)| + h |sin(θ − t)|).  The wall clearance becomes H (w |cos θ| + h |sin θ|).
- **Pair cutoff.**  The circumradius becomes H √(2 cosh a).  Use cut = √(2 cosh ε) + skin, which is negligible for
  ε ≤ 0.3.
- **Polish (SLP, face snap, flips): untouched.**  Softness lives only in a pre-stage, as `--round-r0` does: hard
  squares go into the usual loosen → quench → polish.

### 2.3 Schedule (a sibling of `--round-r0`; flags hypothetical)

`--soft-a0 E [--soft-compress] [--hold-f F] [--soft-steps K] [--soft-p P] [--soft-kick SIG] [--soft-lines-only]
[--soft-budget B] [--soft-pen KAPPA]`

1. **(A) Soft compress (optional).**  `quench` with s free and the a_i free, |a| ≤ ε₀.  This gives s_soft.
2. **(B) Hold the box** at F · s_soft.  Holding below k is the point: otherwise hardening can re-expand.
3. **(C) Harden.**  For k = 1..K: ε_k = ε₀ (1 − k/K)^P.  Relax the fixed-side overlap penalty (`relax_fixed`) over
   (x, y, θ, b, separators), with an optional σ kick per step.
4. **(D) Finish rigid.**  At ε = 0, run the usual loosen 1.02 → quench → polish.
5. **Report** s_soft, the held side, the hard violation at the end of (C), and the full-line count after (A) and
   after (C).

### 2.4 Which jam it breaks, and why it probably won't by itself [inference]

**The row condition.**
- A full axis row needs Σ w_i ≤ s.  A full column needs Σ h_j ≤ s.  At first order, a row shortens by ½ Σ a_i over
  its squares (taking a_i < 0 as narrower).
- So the line jams can be relieved in place iff the following sign system is solvable: every full row has
  Σ a < 0, and every full column has Σ a > 0.
- By Gordan's alternative, it is *un*solvable iff some non-empty set of full rows R′ and full columns C′ is
  **closed**: every square on a row of R′ lies on a column of C′, and vice versa.
- **Proof of the obstruction.**  Over a closed set, Σ_rows Σ w + Σ_cols Σ h = Σ_i (w_i + h_i) = Σ_i 2 cosh(a_i/2)
  ≥ 2|S|, so some line has length ≥ k.  The relief there is second order at best.

**Applied to our jams.**
- A lattice grid with c ≥ 1 holes has no closed set.  Every full row crosses the columns that contain holes, and
  those columns can take taller squares.
- So **every k-row jam with n < k² is relieved in place**.  A soft quench of a holed grid reaches s ≈ k − (#hole
  columns) ε/2, with no rearrangement at all.
- The squares needed are tiny: δ of row relief needs only |a| ≈ 2δ/m on the m freeable squares of the row (δ = 0.02,
  m = 8: |a| ≈ 0.005).
- Hardening then restores exactly the same grid, because the deformation is reversible.  This is the main failure
  mode.  Melashenko's Ex. 2 tiling and the remark quoted in §1a point the same way: soft items in a square box like
  being an axis grid *more* than rigid ones do.

**When it won't work, even in principle.**
- (a) **Closed sub-grids**: rows and columns of a staggered axis chain that cover each other.  Both directions are
  jammed, and area-preserving softening is first-order useless there.
- (b) **Jams held by the tilted part**: gadget and channel contacts, not lines.  Softness has nothing line-specific to
  offer there.
- (c) **The usual case:** the soft stage deforms the jam without changing the topology, and hardening undoes the
  deformation.

**So the homotopy needs a topology change while the squares are soft.**  Three sources:
- **Holding the box below k during hardening (2.3 B–C).**  Re-hardened rows then *must* be resolved by
  rearrangement, or left as overlap.
- **Kicks during hardening.**
- **Basin hopping in the soft landscape** (the expensive version: explorer moves with the soft quench, accepting on
  soft s; second stage only).

The rigid analogue of hold-and-harden at k − 0.02 stalled at E ~ 1e-4 with jammed rows (PACKER, 10-04).  Softness
helps only if a gradual re-pressurisation buckles a row differently from a sudden one.

**Lesson from the rounding experiments (PACKER 10-08 to 10-12).**
- Deep softening (corner radius ≥ 0.2) destroys sub-grid structure and recreates rows.
- Shallow softening (≤ 0.1) preserves structure but did not cross: 0 / 3,130 grid crossings at 110 (CLAIMS C06), and
  0 / 552 below the seed in the 241 / 211 rounding control.
- Expect the same shape here: keep ε₀ ≲ 0.2.

### 2.5 Cheap variants

- **Lines only.**  Give b_i variables only to squares on detected full lines (`layout.full_lines` membership, or the
  partly-pinned set from the LP jam test).  Every other square stays rigid.  This isolates the row mechanism and
  avoids deforming the gadget.
- **Total budget.**  Σ a_i² ≤ B² (one smooth ALM row) instead of a per-square box.  Softness then concentrates on the
  most loaded squares; the stress diagnostic in §4 says which ones.
- **Signed / anti-soft (speculative).**  Force line squares *wider* along their line (a > 0 for row squares) at a
  held side.  The row overflows and must stagger or buckle, using the vertical room their shorter height creates.
  Then harden.  This is a topology-forcing move rather than a relief.
- **Uniform aspect equals a rectangle container**, for axis squares only.  The map x → x/w, y → y·w turns w × 1/w
  rectangles into unit squares in an (s/w) × (s·w) box.  Tilted squares map to parallelograms, so the equivalence
  breaks for them.  So **uniform** softening is mostly already available as `--rect-dx/--rect-dy`.  Only
  *non-uniform* softening is new, and the rectangle box is a free control arm (§4, idea 2).

---

## 3. Pre-registrable test

Write this into PACKER.md before launch (project rule).  The earlier suggestion was: n = 54 and 71, jammed starts,
about 50 trials, kill if hardened outputs never go below k.  Its problems:
- **No concurrent control.**  The hold / kick schedule alone, at ε = 0, might cross.  The 10-04 "fails from random at
  54, 71" was measured with the old packer; nobody has measured `fq`-era crossing at 54 or 71.
- **No positive control.**  Without one, we cannot tell "the schedule destroys sub-k structure" from "the schedule
  cannot create it".
- **"Below k" needs "and no full line".**  Otherwise rattler and polish noise can count as a crossing (`pk/crossing.py`
  already defines this).
- **No mechanism measurement.**  A zero result should say *why*: either the topology never changed, or it changed
  without crossing.
- **Too few trials for the cost.**  At ~1–5 s per trial at these n, 50 trials is cheap enough that it should be more.

**Instances.**
| n | k | record (store) | margin k − s | why |
|---|---|---|---|---|
| 54 | 8 | 7.8466672 | 0.153 | deep sub-k (16 tilted); the original suggestion |
| 71 | 9 | 8.9440716 | 0.056 | original suggestion; 23 tilted |
| 110 | 11 | 10.9967834 | 0.0032 | headline cell only; known rigid baseline 0 / 3,130 crossings (C06) |

**Start classes** (`pk/crossing.py` definitions where they exist):
- `grid`: a k × k lattice with random holes.
- `jam` (new, frozen): fq quenches from uniform random starts in a (k + 0.5) box that end at s ≥ k − 1e-9 with ≥ 1
  full line.  50 per n, generated once and saved.  This is the natural failure state.
- `pos`: the record dilated to k + 0.01 and kicked 0.01.  It should cross back.
- `neg`: n′ = k² − 4, where s = k is proved.  Any crossing there is a bug.

**Arms** (paired: same start file and same kick seeds in every arm):
- **A0, rigid control:** the identical schedule (hold, K = 20 fixed-side relaxations, kick 0.01 per step, finish)
  with ε₀ = 0.
- **A1, soft:** ε₀ ∈ {0.05, 0.2}, `--soft-compress`, hold F = 1, K = 20, P = 1, kick 0.01.
- **A2, soft on lines only:** ε₀ = 0.2.
- **A3, rectangle-container control (zero code):** compress in an (s − 0.1) × (s + 0.1) box, then return to the
  square with `--rect-dx/--rect-dy` stepped to 0.  This is the uniform-aspect analogue (§2.5).

**Trials.**  `grid` + `jam`: 100 each per (n, arm), so 200 per (n, arm).  `pos` / `neg`: 30 each.  n = 110: A0 and the
best soft arm only, 100 `jam` trials.  Budget guess ≤ 10 measured CPU-h; check it with a 5-trial pilot first.
Report `cpu_s` from `pk.manifest`.

**Primary metric.**  The crossing rate per (n, arm): polished side < k − 1e-6 with no full line (`full_lines` = 0),
pooled over `grid` + `jam`.  Any crossing is re-checked with `candidate.py` (both checkers).

**Secondary metrics.**
- **Mechanism:** the full-line count at start, after (A), after (C) and at the end, and the share of trials whose
  final line count is below the start's.
- Matching distance from the start.
- The best side reached.
- `pos` return rate (structure preservation).
- CPU per trial.

**Rules (fixed before launch).**
- **Validity.**  `neg`: 0 crossings in every arm, otherwise stop and debug.  `pos`: A1 must return ≥ 50 % below k;
  if not, report "destructive at this ε" rather than "negative".
- **Success.**  A soft arm has ≥ 5 crossings at n = 54 or 71 *and* beats A0 on the same starts (paired, one-sided
  exact McNemar p < 0.01).  Any certified crossing at 110 is a headline result regardless.
- **Kill.**  Both conditions:
  - 0 crossings in A1 and A2 at both 54 and 71 (400 trials per arm; the 95 % upper bound is ~0.75 % per trial);
  - **and** in ≥ 90 % of soft trials, the final line count equals the start's.

  If both hold, the topology never changed, so park the whole family: no ε, K or kick re-tuning.  If lines did
  change but nothing crossed, run one follow-up only: soft hopping (2.4), then decide.
- **Inconclusive.**  Soft > A0 but below the success bar: one replication at 2× trials, then decide.

**Prediction:** 10 % for a soft arm meeting the success bar at 54 or 71; 2 % at 110.  A3 (rectangle box) is as
likely to cross as A1 (no specific reason to expect softness to beat the uniform version).

---

## 4. Side ideas, even if it fails

1. **Soft sensitivity as a free stress diagnostic.**
   - At a rigid KKT point, ∂s*/∂a_i = Σ_k λ_k ∂g_k/∂a_i.  This is the body-frame normal-stress difference of square
     i's contact forces, times about H/2, and the polish LP duals already give it.
   - Rows show up as uniform along-row load, and the tilted part shows its own pattern.
   - Uses: picking row-break or removal targets, and ranking which squares "want" to change shape.  It costs
     nothing.
2. **Container-aspect homotopy with rigid squares.**
   - Squeeze into (k − δ) × (k + Δ) and walk back to a square.
   - Unlike softening, rows along the short side *cannot* be absorbed by deformation, so they must be broken by
     rearrangement.
   - It needs zero new code (`--rect-dx/dy`), and the rectangle library (lib1) already has pieces.  Arguably a
     better first experiment than soft squares; it is arm A3 above.
3. **Least softness as a jam-depth score** (the Litvinchev 2024 inverse problem).
   - ε*(state) = the smallest ε at which the soft quench reaches a fixed target side.
   - It gives a scalar per basin for "how locked by lines" a state is, which could rank or filter explorer seeds.
4. **Shape kicks as an explorer move.**  Random a_i, a held side, then harden.  This kicks in a direction that
   position and angle kicks never take.  Cheap once 2.1 exists.
5. **Don't rank by soft side.**  In-place relief rewards grids most (§2.4), so the soft s is anti-informative about
   sub-k structure.

---

## Addendum 2026-10-12: Melashenko & Romanova (Ukrainian paper)

User-supplied PDF in `../private/literature/Melashenko_Romanova.pdf`, read in full (`pdftotext -layout` in the session
scratchpad, plus a render of the cuboid figure page; not copied here).  Labels as at the top of this note.

**Record [verified, Crossref].**  О.П. Мелащенко, Т.Є. Романова, «Компонування м'яких багатогранників у опуклому
контейнері мінімального об'єму» (English title as printed: "Layout of soft polyhedra in a minimum-volume convex
container").  *Проблеми керування та інформатики* / *Problems of Control and Informatics* 70(1) (2025) 5–21,
doi:10.34229/1028-0979-2025-1-1 (ISSN 2786-6491; received 23.12.2024).  **It is about 3D polyhedra, not polygons.**
Related items that turned up in the search but were not read: Melashchenko, Romanova, Shekhovtsov, "Packing soft
polygons in a minimum height rectangular target domain", *Cybernetics and Computer Technologies* 2024(1) 5–17,
doi:10.34229/2707-451x.24.1.1 (probably the 2D "soft polygons" paper), and "Layout of compressible objects", LNNS
(Smart Technologies in Urban Engineering), doi:10.1007/978-3-032-06829-3_29.

### What it does [read]

- **Two formulations.**
  - **LSP1, free vertices.**  The variables are the vertex coordinates X_i^k themselves (no pose).  Elasticity is an
    edge-length band (1 − ε) d_i ≤ ‖X_i^{k+1} − X_i^k‖ ≤ (1 + ε) d_i, eq. (1)/(9), with cyclic vertex indexing and one
    reference length d_i per item; ε ∈ {0.3, 0.5, 0.6} or "–" (band dropped).  Volume is held by an equality on the
    6 × det formula for a pyramid, eq. (4)/(8).  "Convexity preservation" is listed as a condition, but no
    constraint enforces it in (5)–(9); for a 4-vertex pyramid it is automatic [inference].  Model (5)–(9), with
    3Σm_i + 3n(n−1)/2 + |L| variables.
  - **LSP2, pose plus stretch.**  P_i(u_i, μ_i) = {v_i + A(μ_i) M(θ_i) p⁰}, with A = diag(μ^x, μ^y, μ^z) and M a
    3-angle rotation.  **The stretch is applied after the rotation, along the container axes** (vertex formula
    eq. (10)), unlike the body-frame A(θ)M(μ) order of the *Mathematics* paper (§1a).  Volume and convexity are kept by
    μ^x μ^y μ^z = 1, eq. (14) (stated, as in §1a's Prop. 1).  Bounds μ^− ≤ μ ≤ μ^+, eq. (15).  **In every experiment
    μ^x = μ_i, μ^y = 1/μ_i, μ^z = 1**: the 2D stretch of §1a, embedded in 3D.  Model (11)–(15),
    (3n² + 15n)/2 + |L| variables.
- **Non-overlap.**  A separating plane per pair.  LSP1 uses (s_ij, γ_ij) with ‖s_ij‖ = 1, eq. (2)/(6), which they
  call a generalised quasi phi-function.  LSP2 parametrises the unit normal by two angles (θ_ij^x, θ_ij^y), so the
  norm constraint disappears; the quasi phi-function is the min of the two vertex-side minima (text after eq. (10),
  constraint (12)).  This is the same separator as §1a and `fq`.  Containment: vertices in the cuboid, sphere or
  cylinder (after (9) and after (15)).
- **Objective.**  Container volume: abc, r, or r²h (after (9)).  The container is not a fixed shape: a cuboid has
  three free sides.
- **Solver.**  Multistart with **three starts per example** (LSP2 section).  3D-FSPA (Alg. 1, (16)–(19)): random
  centres, circumscribed spheres grown homothetically by IPOPT, rescale by 1/λ, random angles, separating planes
  from the sphere pairs.  3D-SDA (Alg. 2, (20)–(25)): per-item cube neighbourhoods of side 2(r_i + ε) (the ε formula
  is garbled in extraction; it looks like §1a's eq. (14)), far pairs dropped, IPOPT, stop when the objective change is
  ≤ δ (step 7).  Python 3.12 on an i5-12500.
- **Instances and results.**  Tables 1–3 (LSP1, pyramids) and 4–9 (LSP2, pyramids and cuboids) in a cuboid, sphere
  or cylinder, n = 10–40, times of 1–985 s.  Fill: pyramids 0.39–0.89; cuboids in a sphere or cylinder 0.43–0.72;
  cuboids in a cuboid (Table 7) 1.0000 at n = 10 (5 × 4 × 4), then 0.983 (n = 20), 0.852 (n = 23, a 6 × 6 × 6 box)
  and 0.849 (n = 30).  Bound ranges change from row to row, and include **[1.5, 2.5], which excludes μ = 1**: the
  original shape is then infeasible (Table 4, n = 30).
- **Hardening / continuation / schedule back to rigid:** none.  **Rigid baseline:** none, not even the single one of
  Bennell 2024.  **Squares:** none (in 2D or as cube faces), apart from the cuboids below.

### Equal objects: implied by the tables, not stated [inference]

- f* · V* = Σ V_i gives exactly **8n for every cuboid row** (Tables 7–9) and **13.856 n ≈ 8√3 n for every pyramid
  row** (Tables 1–6, e.g. 166.0213 × 0.8346 = 138.56 at n = 10).  Fig. 10 shows equal-height layers (c = 4 is two
  layers, with μ^z = 1).  So the items are very likely identical; the most likely reading is 2 × 2 × 2 cubes.
- **If so, Table 7 is worse than a trivial rigid packing.**  n identical bricks stacked in a row fill an n·p × q × r
  box exactly (f = 1).  μ = 1 is inside every Table 7 range, so that packing is feasible for their problem, yet they
  report f = 0.983, 0.852 and 0.849 at n = 20, 23 and 30.  For 2-cubes, their n = 23 answer, 6 × 6 × 6, is exactly the
  rigid 3 × 3 × 3 lattice box with four holes.  Fig. 10б shows axis-parallel bricks in layers.  Their solver (three
  starts, local IPOPT) stalls in a lattice-like state, with no rigid comparison to expose it.

### What it adds beyond §1a and §1b

1. 3D, and a second elasticity model (free vertices with an edge-length band) next to the affine stretch.
2. Container-frame stretch (eq. (10)) instead of body-frame.
3. Probably identical items, in a container whose shape is free.
4. Nothing on hardening, continuation, rigid baselines or 2D squares.  The "where the evidence is thin" list in §1
   stands, and is strengthened.  Their only box-in-box instances look like lattice stalls that a rigid packing beats.

### Effect on the design and the test (§2–§3) [inference]

- **The test in §3 does not change.**  No arm, metric, instance or rule needs revision.  The paper reinforces three
  points already there: the concurrent rigid control A0 is essential (the paper has none, and its soft results lose
  to rigid ones), soft items stall on axis lattices (§2.4, failure (c)), and soft side should not be used for ranking
  (§4 item 5).
- **Frame of the stretch: record as an optional variant in §2.5, not as a change.**  Container-frame stretch,
  o_k = diag(w_i, h_i) · H R(θ_i) s_k, instead of §2.1's body-frame H R(θ_i) diag(w_i, h_i) s_k:
  - For axis squares the two are the same, so §2.4's row and Gordan argument is unchanged.
  - For tilted squares, a container-frame stretch makes a parallelogram.  It shortens the square's x- (or y-) extent
    directly, whatever its tilt.  So a tilted square on a near-full line can give relief at first order; in the body
    frame it gives less relief, and none at 45°, where w|cos θ| + h|sin θ| is even in a.  This bears only on failure
    mode (b) (jams held by the tilted part).
  - It also removes the a ↔ −a, θ + 90° redundancy, and the corners stay linear in (w, h), so separators and walls
    stay smooth.  `sat_gap` and `check` would need the parallelogram projection max_k t · o_k rather than the
    rectangle formula in §2.2.
  - Worth one extra arm (A1 with the container frame) only if the pilot shows A1 changing lines in the tilted part
    and not crossing.
- **Free vertices (LSP1 model): not adopted.**  For squares it means 8 coordinates per item and shear modes, and in 2D
  it needs an explicit convexity constraint (not automatic as for a tetrahedron).  This is a larger and less
  structured search space than ε tanh(b), with no evidence in favour.
- **Forced softness** (a bound range excluding μ = 1, Table 4) is the paper's analogue of the "signed / anti-soft"
  variant in §2.5.  They give no comparison, so it adds no evidence either way.
