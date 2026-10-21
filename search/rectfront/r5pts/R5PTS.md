# R₅, short side ≥ 2: the 4-point unavoidable set (numerics + paper proofs), 2026-10-21

Scripts here; outputs `runs/rectfront/r5pts/`.  Labels: [measured] float optimiser / sampling; [proved on paper]
argument below, each step also checked on 4·10⁵ random samples (`checkproofs.py`, 0 failures).

Notation: a unit square with centre m and angle θ; for a point P, d = P − m, u = d·(c, s), v = d·(−s, c)
(c = cos θ, s = sin θ; same as `coord` in `Basic.lean`).  P ∈ closed square ⟺ |u|, |v| ≤ ½.  The square is
unchanged by θ ↦ θ + 90°, so take c, s ≥ 0 (a frame, as in `RectChain.framed`); k = (c + s)/2, p = c + s ∈ [1, √2].
"Outside" a corner/point below means outside the *closed* square: |u| > ½ or |v| > ½.

## 1. Centre-rectangle lemma [measured; proved on paper]

**Statement (closed).**  X, Y ≥ 0, X ≤ 1, Y ≤ 1, X + Y ≤ √2.  If m ∈ [x₁, x₁+X] × [y₁, y₁+Y], the closed unit
square (any θ) contains one of the four corners.  (Strict version: X, Y < 1, X + Y < √2 ⇒ a corner in the open
interior; same proof.)

[measured] `crect.py`: G(X, Y) := max over m, θ of min over corners of max(|u|, |v|) equals
max(max(X, Y)/2, (X + Y)/(2√2)) (11 X values, on / ±0.01 off the boundary, Nelder–Mead from the 20 best grid points);
worst pose = rectangle centre, at 45° on X + Y = √2 and upright on max = 1.  So the conditions are exactly sharp.

**Proof.**  Write m's offsets to the corners: d₁ = (A, B), d₂ = (−A′, B), d₃ = (−A′, −B′), d₄ = (A, −B′) with
A, A′, B, B′ ≥ 0, A + A′ = X, B + B′ = Y (corners NE, NW, SW, SE, cyclic; adjacent corners differ by (X, 0) or (0, Y)).
WLOG c ≥ s ≥ 0 (else reflect x ↦ −x, which maps the rectangle to a rectangle and (c, s) to the frame (s, c)).
Suppose every corner is outside.
1. *Witnesses.*  Because A, B, … ≥ 0 and c ≥ s ≥ 0, each corner has only two possible witnesses:
   d₁: u > ½ or v > ½ (−u₁ ≤ 0, −v₁ = As − Bc ≤ u₁);  d₂: v > ½ or −u > ½ (−v₂ ≤ 0, u₂ ≤ v₂);
   d₃: −u > ½ or −v > ½ (u₃ ≤ 0, v₃ ≤ −u₃);  d₄: u > ½ or −v > ½ (v₄ ≤ 0, −u₄ ≤ −v₄).
   Index the witnesses u, v, −u, −v as 0, 1, 2, 3: corner j's witness is j−1 ("first") or j ("second") mod 4.
2. *Adjacent corners cannot have opposite witnesses*: e.g. u₁ > ½, −u₂ > ½ ⇒ u₁ − u₂ = Xc > 1, but X ≤ 1;
   likewise v₂ − v₃ = Yc, u₄ − u₃ = Xc, v₁ − v₄ = Yc, all ≤ 1.
3. If the choices are not all equal, somewhere cyclically corner j takes "first" (j−1) and corner j+1 "second"
   (j+1): opposite, excluded by 2.  So all first or all second.
4. *All first* (u₁, v₂, −u₃, −v₄ > ½): (u₁ − u₃) + (v₂ − v₄) = (Xc + Ys) + (Xs + Yc) = (X + Y)(c + s) > 2, but
   X + Y ≤ √2 and c + s ≤ √2.  *All second* (v₁, −u₂, −v₃, u₄ > ½): (v₁ − v₃) + (u₄ − u₂) = (Yc − Xs) + (Xc − Ys) =
   (X + Y)(c − s) > 2, but ≤ √2.  Contradiction.  ∎

So 16 sign cases, each one `linarith` (14 by step 2, 2 by step 4).  This is DS7's quadrant idea (Lemma 3) made
discrete; it replaces the triangle lemma (X² + Y² ≤ 1) by the sharp X + Y ≤ √2.

## 2. Side-strip lemma [measured; proved on paper]

**Statement (closed; DS7 Lemma 2 orientation).**  Wall y = 0, the square in y ≥ 0 (so m₂ ≥ k, `half_extent`).
Points P = (0, y₀), Q = (x₀, y₀).  Hypotheses: 0 ≤ x₀ ≤ 1, y₀ ≤ 1, x₀ + 2y₀ ≤ 2√2.  If 0 ≤ m₁ ≤ x₀ and m₂ ≤ y₀,
the closed square contains P or Q.  (Our use: y₀ = 1, x₀ = gap = b − 2 on a long wall, a − 2 on a short wall;
so the condition is gap ≤ 2√2 − 2.)  Only the wall the strip lies on is used.

**Proof.**  WLOG c ≥ s ≥ 0 (else reflect x ↦ x₀ − x: swaps P and Q).  d_Q = (A, B), d_P = (−A′, B) with
A = x₀ − m₁, A′ = m₁, B = y₀ − m₂ ≥ 0, A + A′ = x₀.  Witnesses as in §1: Q: u or v; P: v or −u.
* Q with v: v_Q = Bc − As ≤ Bc ≤ (y₀ − p/2)c ≤ (1 − p/2) ≤ ½ (p ≥ 1).  Impossible (the wall).
* Q: u, P: −u: u_Q − u_P = x₀c > 1, impossible (x₀ ≤ 1).
* Q: u, P: v: s·u_Q + c·v_P = cs·x₀ + B > p/2 (weights s, c kill A), and B ≤ y₀ − p/2, so
  g(p) := x₀(p² − 1)/2 + y₀ − p > 0 (cs = (p² − 1)/2).  But g is convex in p with g(1) = y₀ − 1 ≤ 0 and
  g(√2) = x₀/2 + y₀ − √2 ≤ 0; exactly: g(p) = [(√2 − p)g(1) + (p − 1)g(√2)]/(√2 − 1) + (x₀/2)(p − 1)(p − √2) ≤ 0.  ∎

Sharp at θ = 45° (DS7).  Note it does need the closed form: at θ = 0, y₀ = 1 the upright square on the wall has
P, Q on its edge (g(1) = 0), so **no open-interior version exists with y₀ = 1** — that is why Göbel pulls the points
in by ε and DS7 scales.  [measured: `checkproofs.py` conclusion + every case inequality, 0 failures.]

## 3. Corner lemma [proved on paper]
Walls x = 0, y = 0, point (1, 1), centre m ∈ [k, 1]² (k from `half_extent`).  A = 1 − m₁, B = 1 − m₂ ∈ [0, 1 − k]:
|u| = Ac + Bs ≤ (1 − k)(c + s) = p − p²/2 ≤ ½ ((p − 1)² ≥ 0), |v| = |Bc − As| ≤ max(Bc, As) ≤ 1 − k ≤ ½.  No WLOG
needed.  (Closed twin of `mem_sqInt_of_close`.)

## 4. Assembly and strictness (recommendation)

**Cover lemma (closed).**  2 ≤ a ≤ 2√2, 2 ≤ b ≤ 2√2, a + b ≤ 4 + √2: every closed unit square in `rbox a b` contains
one of (1, 1), (a − 1, 1), (1, b − 1), (a − 1, b − 1).  Split m₁ at 1 and a − 1, m₂ at 1 and b − 1 (needs a, b ≥ 2):
4 corner cells (§3, + mirrors), 4 side cells (§2 with y₀ = 1, x₀ = a − 2 or b − 2 ≤ 2√2 − 2 < 1), centre cell
(§1, X = a − 2, Y = b − 2: X, Y ≤ 2√2 − 2 < 1, X + Y ≤ √2).  Cells overlap on their edges; harmless (closed).

**Lower bound.**  a, b ≥ 2, a, b < 2√2, a + b < 4 + √2 ⇒ ¬ RectPacks 5 a b.  Take λ > 1 with λa, λb, λ(a + b)
still within the cover lemma's region; the scaled packing has λ-squares in `rbox (λa) (λb)`; apply
`packing_le_weight` (Basic.lean: closed unit squares concentric with side-λ squares, weights 1 on the 4 points)
⇒ 5 ≤ 4.  Needs one small scaling lemma (`coord (λc) θ (λp) = λ • coord c θ p`, so `sq (λc) θ λ = λ • sq c θ 1`).
Recommended over (i) pulled-in points (ε bookkeeping in every lemma) and (ii) strict hypotheses (impossible, §2).
The theorem covers R₅'s piece a + b = 4 + √2 for 4 − √2 ≤ a ≤ 2 + 1/√2 (there b < 4 + √2 − a ≤ 2√2), including
s(5) ≥ 2 + 1/√2, and also the region a ∈ [2, 4 − √2), b < 2√2 (below the true boundary b = 3 there).

## 5. Coverage of the rectangular set [measured] (`cover.py`, `cover_dense.py`, `escape.py`)
Escape = max over squares in the container of min over points (gauge − ½); ≤ 0 ⟺ unavoidable (closed).
* 208 containers (16 × 8 grid + 80 random) in 4 − √2 ≤ a ≤ b ≤ 4 + √2 − a: max escape 2.2·10⁻¹⁶ (0 = upright
  squares touching a point).
* Just beyond the line (+10⁻³, +10⁻²) for a ≥ 4 − √2: escape = Δ/(2√2) (45° square at the centre), as §1 predicts.
* For a < 4 − √2 on or near the line: escape > 0 (45° square on the long wall, centre (k, b/2)), up to 0.05 at
  a = 2.45; matches §2's condition b ≤ 2√2.

## 6. The corner a = 1 + √2, b → 3 [measured] (`escapers.py`, `ptopt.py`, `ptopt2.py`)
* **Escapers of the rectangular set (p = q = 1)** at b = 2.95, 2.99, 2.999: only in the two long-wall side cells,
  hugging the wall: centre within ≈ 0.06–0.085 of its wall position k (x ≤ k + 0.085), height within ≈ 0.1 of b/2,
  **any tilt** in (0°, 90°) at b → 3 (b = 2.95: 7°–83°; b = 2.999: 0.35°–89.65°), largest escape at 45°
  (0.043, 0.057, 0.060).  So not "near 45°" only: as b → 3 the near-upright wall sitter (touching (1, 1) and
  (1, b − 1) with its edge) also escapes.  No escapers in the corner, short-side or centre cells.
* **Best 4-point sets** (multistart Nelder–Mead on a pose sample, then 8 rounds of cutting planes with the full
  escape optimiser; families: rectangular (p, q), mirror-symmetric, left/right-mirror, point-symmetric, free 8-dim):
  | (a, b) | rect | mirror | psym | free |
  |---|---|---|---|---|
  | (2.4142, 2.99) | +0.0157 (p 0.957, q 1.016) | +0.0161 | +0.0159 | +0.0160 |
  | (2.40, 2.99) | +0.0142 | +0.0145 | +0.0144 | +0.0154 |
  All families converge to ≈ the rectangular set with p ≈ 0.95 (squeezed between the side strip and the centre
  cell), escape ≈ 0.015 > 0.  Local search only, but every start and family lands at the same value: **no 4-point
  unavoidable set near the corner** (consistent with the literature note's crude search).  Pigeonhole on 4 points
  does not reach the b = 3 piece at a ≈ 2.4; a second stage on the long-wall sitters is the natural next step.
