import Sqpack.SegParts

/-!
# Zero-margin leaves for mixed covers (points + axis-parallel segments)

The sibling of `ZMTree.lean` for mixed covers (`certificates/s21/FORMAT.md`, `search/zm_mixed.py`,
`search/ZM_MIXED.md`).  A tree `ZTM` is checked by `checkM` (natural-number arithmetic only) and
certifies `CovM` (`CovM.lean`); `soundM` holds for every tree, so nothing about the tree or the
generator is trusted.

## The mixed `Z` leaf

A leaf carries a `ZMTree` point leaf (claimed points, chains, emptiness staircase) and a **piece
certificate** `pc`.  The kernel computes the certificate's value `Lp` (a lower bound, in units of
`1/W`, on the segment mass captured at every admissible pose of the box) and runs the unchanged
`ZMTree.zOk` with the target weight `W − Lp` (**Lemma P**: the pieces act as a phantom point of weight
`Lp` present in every region).  `Lp ≥ W` is the `PIECE` leaf (no points needed).

The piece certificate: claimed segments (gaps into the segment candidate list, each with a *tag*
naming the block it belongs to), and two kinds of blocks.

* **S-block** `(dir, K, A, B)` (**Lemma S**, certified core of a line): the points of the line
  (`x = K/Q` if `dir = 0`, `y = K/Q` otherwise) with line coordinate in `[A/Q, B/Q]` lie in the square
  at every admissible pose.  Checked by `ZMTree.admAll` (all four conditions) at the two ends; the
  certified set is convex (each condition is a half-plane at every pose), so the whole interval is
  certified.  A claimed segment tagged with the block contributes `⌊w · |seg ∩ [A,B]| / |seg|⌋`.
* **T-group** (**Lemma T**, germ pair): two lines one unit apart (`x = ξ`, `x = ξ + 1`, or
  `y = η + 1`, `y = η`), each with an interval certified for its three non-singular conditions
  (`admAll` with the singular kind skipped), optional points `τ` certified for the singular
  condition, and a list of *pairs* `(up piece [a,b], down piece [c,d], m)`: each piece inside a
  claimed segment of its line, `m ≤` the mass of each piece, `c ≤ a + u₀`, `d ≤ b + u₀`, pieces sorted
  along each line.  At every pose each pair captures `≥ m` (`pair_capture`): the group bound is `Σ m`.
  (A pair may have one piece only when that piece is above `τ` on the up line or below `τ` on the
  down line.)  This is the matching (dual) form of `ZM_MIXED.md`'s Corollary T.

Soundness of the value: at a pose, the certificate yields a list of *parts* (`SegParts.lean`),
pairwise non-overlapping on each segment, whose masses sum to at least `Lp`; `parts_le_segMass`.
-/

open MeasureTheory

namespace SquarePacking

namespace ZMTreeM

open BoxTree ZMTree

/-! ## 1.  `admAll` and convexity along a line -/

/-- **Soundness of `ZMTree.admAll`**: the non-swing conditions (all four for `kp = 4`) hold at every
admissible pose of the box. -/
theorem admAll_sound {S Q R x0 x1 y0 y1 U0 U1 kp XS YS : ℕ} {c : ℝ × ℝ} {u : ℝ} (hQ : 0 < Q)
    (hR : 0 < R) (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u)
    (h : admAll Q R x0 x1 y0 y1 U0 U1 kp XS YS = true) :
    ∀ k : Fin 4, k.val ≠ kp → Gv Q k.val XS YS c u ≤ 0 := by
  intro k hk
  simp only [admAll, Bool.or_eq_true, Bool.and_eq_true, Nat.beq_eq] at h
  rcases h with hf | ⟨⟨⟨h0, h1⟩, h2⟩, h3⟩
  · rw [Gv, kf_val]
    exact ptOkK_sound hQ hR P.hU01 P.hU1 hf c u P.hu0 P.hu1 P.hx0 P.hx1 P.hy0 P.hy1 k hk
  · have A := fun (k : ℕ) (hk : admK Q R x0 x1 y0 y1 U0 U1 XS YS k = true) =>
      admK_sound hQ hR P.hU01 P.hU1 hk P.hx0 P.hx1 P.hy0 P.hy1 P.hu0 P.hu1 P.hwx P.hwy
    fin_cases k
    · exact A 0 (h0.resolve_left (by simp at hk ⊢; omega))
    · exact A 1 (h1.resolve_left (by simp at hk ⊢; omega))
    · exact A 2 (h2.resolve_left (by simp at hk ⊢; omega))
    · exact A 3 (h3.resolve_left (by simp at hk ⊢; omega))

lemma beq_ne {a b : ℕ} (h : a ≠ b) : Nat.beq a b = false := by
  cases hc : Nat.beq a b
  · rfl
  · exact absurd (Nat.beq_eq.mp hc) h

lemma beq_rfl (a : ℕ) : Nat.beq a a = true := Nat.beq_eq.mpr rfl

/-- The coordinates over `Q` of the point of line coordinate `t` on the line `(dir, K)`. -/
def lx (dir K t : ℕ) : ℕ := bif Nat.beq dir 0 then K else t
def ly (dir K t : ℕ) : ℕ := bif Nat.beq dir 0 then t else K

open Classical in
/-- The real point of line coordinate `t` on the line `x = K/Q` (`dir = 0`) or `y = K/Q`. -/
noncomputable def lpt (Q dir K : ℕ) (t : ℝ) : ℝ × ℝ :=
  if dir = 0 then ((K : ℝ) / Q, t) else (t, (K : ℝ) / Q)

lemma lpt_nat (Q dir K t : ℕ) :
    lpt Q dir K ((t : ℝ) / Q) = ((lx dir K t : ℝ) / Q, (ly dir K t : ℝ) / Q) := by
  by_cases hd : dir = 0
  · subst hd; simp [lpt, lx, ly]
  · simp [lpt, lx, ly, hd, beq_ne hd]

/-- A violation polynomial along a line is affine in the line coordinate. -/
lemma gval_lpt_affine (Q dir K : ℕ) (k : Fin 4) (c : ℝ × ℝ) (u : ℝ) :
    ∃ α β : ℝ, ∀ t, gval k ((lpt Q dir K t).1 - c.1) ((lpt Q dir K t).2 - c.2) u = α + β * t := by
  by_cases hd : dir = 0
  · refine ⟨-(1 + u ^ 2) + galpha k u * ((K : ℝ) / Q - c.1) - gbeta k u * c.2, gbeta k u,
      fun t => ?_⟩
    simp only [lpt, hd, if_true, gval_eq]; ring
  · refine ⟨-(1 + u ^ 2) - galpha k u * c.1 + gbeta k u * ((K : ℝ) / Q - c.2), galpha k u,
      fun t => ?_⟩
    simp only [lpt, hd, if_false, gval_eq]; ring

/-- **Convexity along a line**: a condition that holds at two points of a line holds between. -/
lemma gval_lpt_between {Q dir K : ℕ} {k : Fin 4} {c : ℝ × ℝ} {u t1 t2 t : ℝ}
    (h1 : gval k ((lpt Q dir K t1).1 - c.1) ((lpt Q dir K t1).2 - c.2) u ≤ 0)
    (h2 : gval k ((lpt Q dir K t2).1 - c.1) ((lpt Q dir K t2).2 - c.2) u ≤ 0)
    (ht1 : t1 ≤ t) (ht2 : t ≤ t2) :
    gval k ((lpt Q dir K t).1 - c.1) ((lpt Q dir K t).2 - c.2) u ≤ 0 := by
  obtain ⟨α, β, hf⟩ := gval_lpt_affine Q dir K k c u
  rw [hf] at h1 h2 ⊢
  rcases le_total 0 β with hb | hb
  · nlinarith [mul_le_mul_of_nonneg_left ht2 hb]
  · nlinarith [mul_le_mul_of_nonpos_left ht1 hb]

/-- `Gv` at the integer line point is `gval` at the real line point. -/
lemma Gv_lpt (Q dir K t : ℕ) (k : Fin 4) (c : ℝ × ℝ) (u : ℝ) :
    Gv Q k.val (lx dir K t) (ly dir K t) c u
      = gval k ((lpt Q dir K ((t : ℝ) / Q)).1 - c.1) ((lpt Q dir K ((t : ℝ) / Q)).2 - c.2) u := by
  rw [Gv, kf_val, lpt_nat]

/-- An interval of a line certified at its two ends by `admAll` (conditions `≠ kp`). -/
lemma certified_interval {S Q R x0 x1 y0 y1 U0 U1 kp dir K A B : ℕ} {c : ℝ × ℝ} {u : ℝ}
    (hQ : 0 < Q) (hR : 0 < R) (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u)
    (hA : admAll Q R x0 x1 y0 y1 U0 U1 kp (lx dir K A) (ly dir K A) = true)
    (hB : admAll Q R x0 x1 y0 y1 U0 U1 kp (lx dir K B) (ly dir K B) = true) {t : ℝ}
    (ht1 : (A : ℝ) / Q ≤ t) (ht2 : t ≤ (B : ℝ) / Q) (k : Fin 4) (hk : k.val ≠ kp) :
    gval k ((lpt Q dir K t).1 - c.1) ((lpt Q dir K t).2 - c.2) u ≤ 0 := by
  have a := admAll_sound hQ hR P hA k hk
  have b := admAll_sound hQ hR P hB k hk
  rw [Gv_lpt] at a b
  exact gval_lpt_between a b ht1 ht2

/-! ## 2.  Segments on a line -/

/-- The entry lies on the line `(dir, K)` (over `Q = D·S`), normalised and non-degenerate. -/
def onLine (S dir K : ℕ) (e : SegE) : Bool :=
  bif Nat.beq dir 0 then Nat.beq e.1 e.2.2.1 && Nat.beq (Nat.mul e.1 S) K && Nat.blt e.2.1 e.2.2.2.1
  else Nat.beq e.2.1 e.2.2.2.1 && Nat.beq (Nat.mul e.2.1 S) K && Nat.blt e.1 e.2.2.1

/-- The ends of the entry in its line coordinate, over `Q`. -/
def slo (S dir : ℕ) (e : SegE) : ℕ := Nat.mul (bif Nat.beq dir 0 then e.2.1 else e.1) S
def shi (S dir : ℕ) (e : SegE) : ℕ := Nat.mul (bif Nat.beq dir 0 then e.2.2.2.1 else e.2.2.1) S

/-- What `onLine` gives in real terms. -/
lemma onLine_spec {D S Q dir K : ℕ} (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) {e : SegE}
    (h : onLine S dir K e = true) :
    (isVert e ∨ isHorz e) ∧ elo D e = (slo S dir e : ℝ) / Q ∧
      ehi D e = (shi S dir e : ℝ) / Q ∧ slo S dir e < shi S dir e ∧
      ∀ t, lpE D e t = lpt Q dir K t := by
  subst hQD
  have hSr : (S : ℝ) ≠ 0 := by exact_mod_cast hS.ne'
  have hDr : (D : ℝ) ≠ 0 := by exact_mod_cast hD.ne'
  by_cases hd : dir = 0
  · subst hd
    simp only [onLine, beq_rfl, cond_true, Bool.and_eq_true, Nat.beq_eq,
      Nat.blt_eq, nat_mul_eq] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    refine ⟨Or.inl ⟨h1, h3⟩, ?_, ?_, ?_, fun t => ?_⟩
    · simp only [elo, h1, if_true, slo, beq_rfl, cond_true, nat_mul_eq]
      push_cast; field_simp
    · simp only [ehi, h1, if_true, shi, beq_rfl, cond_true, nat_mul_eq]
      push_cast; field_simp
    · simp only [slo, shi, beq_rfl, cond_true, nat_mul_eq]
      exact Nat.mul_lt_mul_of_pos_right h3 hS
    · simp only [lpE, h1, if_true, lpt]
      rw [← h2]; push_cast
      congr 1; field_simp; exact_mod_cast h1.symm
  · have hb : Nat.beq dir 0 = false := beq_ne hd
    simp only [onLine, hb, cond_false, Bool.and_eq_true, Nat.beq_eq, Nat.blt_eq,
      nat_mul_eq] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    have hne : e.1 ≠ e.2.2.1 := ne_of_lt h3
    refine ⟨Or.inr ⟨h1, h3⟩, ?_, ?_, ?_, fun t => ?_⟩
    · simp only [elo, hne, if_false, slo, hb, cond_false, nat_mul_eq]
      push_cast; field_simp
    · simp only [ehi, hne, if_false, shi, hb, cond_false, nat_mul_eq]
      push_cast; field_simp
    · simp only [slo, shi, hb, cond_false, nat_mul_eq]
      exact Nat.mul_lt_mul_of_pos_right h3 hS
    · simp only [lpE, hne, if_false, lpt, hd]
      rw [← h2]; push_cast
      congr 1; field_simp

/-! ## 3.  The piece certificate -/

/-- An S-block `(dir, K, A, B)`. -/
abbrev SBlk := ℕ × ℕ × ℕ × ℕ
/-- A pair `(ju+1, a, b, jd+1, c, d, m)` of a T-group (`0` = no piece on that line). -/
abbrev TPr := ℕ × ℕ × ℕ × ℕ × ℕ × ℕ × ℕ
/-- A T-group `(dir, K, Au, Bu, Ad, Bd, τu+1, τd+1, pairs)` (`0` = no `τ`). -/
abbrev TGrp := ℕ × ℕ × ℕ × ℕ × ℕ × ℕ × ℕ × ℕ × List TPr

def sb0 : SBlk := (0, 0, 0, 0)
def ec0 : SegE × ℕ := ((0, 0, 0, 0, 0), 0)

/-- The claimed segments: gaps into the candidate list, with tags. -/
def sclaim : List (ℕ × ℕ) → List SegE → List (SegE × ℕ)
  | [], _ => []
  | (k, t) :: ks, c =>
    match c.drop k with
    | [] => []
    | e :: rest => (e, t) :: sclaim ks rest

lemma sclaim_sub : ∀ (sel : List (ℕ × ℕ)) (c : List SegE), ((sclaim sel c).map Prod.fst).Sublist c
  | [], _ => by simp [sclaim]
  | (k, t) :: ks, c => by
    simp only [sclaim]
    split
    · simp
    · rename_i e rest hdrop
      simp only [List.map_cons]
      have h1 : (e :: rest).Sublist c := hdrop ▸ List.drop_sublist k c
      exact ((sclaim_sub ks rest).cons_cons e).trans h1

/-- An S-block is certified at its two ends. -/
def sblkOk (Q R x0 x1 y0 y1 U0 U1 : ℕ) (b : SBlk) : Bool :=
  Nat.ble b.2.2.1 b.2.2.2 &&
    admAll Q R x0 x1 y0 y1 U0 U1 4 (lx b.1 b.2.1 b.2.2.1) (ly b.1 b.2.1 b.2.2.1) &&
    admAll Q R x0 x1 y0 y1 U0 U1 4 (lx b.1 b.2.1 b.2.2.2) (ly b.1 b.2.1 b.2.2.2)

/-- The mass a segment certifies in an S-block: `⌊w · |seg ∩ [A,B]| / |seg|⌋`. -/
def sval (S : ℕ) (b : SBlk) (e : SegE) : ℕ :=
  Nat.div (Nat.mul e.2.2.2.2 (Nat.sub (min (shi S b.1 e) b.2.2.2) (max (slo S b.1 e) b.2.2.1)))
    (Nat.sub (shi S b.1 e) (slo S b.1 e))

/-- An even tag `2s` names S-block `s`; the segment must lie on its line. -/
def sclOk (S : ℕ) (sb : List SBlk) (ec : SegE × ℕ) : Bool :=
  bif Nat.beq (Nat.mod ec.2 2) 0 then
    Nat.blt (Nat.div ec.2 2) sb.length &&
      onLine S (sb.getD (Nat.div ec.2 2) sb0).1 (sb.getD (Nat.div ec.2 2) sb0).2.1 ec.1
  else true

def sclVal (S : ℕ) (sb : List SBlk) (ec : SegE × ℕ) : ℕ :=
  bif Nat.beq (Nat.mod ec.2 2) 0 then sval S (sb.getD (Nat.div ec.2 2) sb0) ec.1 else 0

/-- The kinds and lines of a germ pair: vertical (`dir = 0`): up `x = K/Q` (kind 1), down
`x = K/Q + 1` (kind 0); horizontal: up `y = K/Q + 1` (kind 2), down `y = K/Q` (kind 3). -/
def kU (dir : ℕ) : ℕ := bif Nat.beq dir 0 then 1 else 2
def kD (dir : ℕ) : ℕ := bif Nat.beq dir 0 then 0 else 3
def lnU (Q dir K : ℕ) : ℕ := bif Nat.beq dir 0 then K else Nat.add K Q
def lnD (Q dir K : ℕ) : ℕ := bif Nat.beq dir 0 then Nat.add K Q else K

/-- A piece `[a, b]` of the claimed segment `j1 − 1` (tag `tag`) on the line `(dir, L)`, inside the
certified interval `[lo, hi]`, of mass `≥ m`. -/
def pieceOk (S dir L lo hi tag : ℕ) (cls : List (SegE × ℕ)) (j1 a b m : ℕ) : Bool :=
  Nat.beq (cls.getD (Nat.sub j1 1) ec0).2 tag && onLine S dir L (cls.getD (Nat.sub j1 1) ec0).1 &&
    Nat.ble (slo S dir (cls.getD (Nat.sub j1 1) ec0).1) a &&
    Nat.ble b (shi S dir (cls.getD (Nat.sub j1 1) ec0).1) && Nat.ble lo a && Nat.ble b hi &&
    Nat.ble (Nat.mul m (Nat.sub (shi S dir (cls.getD (Nat.sub j1 1) ec0).1)
      (slo S dir (cls.getD (Nat.sub j1 1) ec0).1)))
      (Nat.mul (cls.getD (Nat.sub j1 1) ec0).1.2.2.2.2 (Nat.sub b a))

/-- A pair of a T-group. -/
def tprOk (S Q R U0 : ℕ) (cls : List (SegE × ℕ)) (tag : ℕ) (g : TGrp) (p : TPr) : Bool :=
  Nat.ble p.2.1 p.2.2.1 && Nat.ble p.2.2.2.2.1 p.2.2.2.2.2.1 &&
  (Nat.beq p.1 0 ||
    pieceOk S g.1 (lnU Q g.1 g.2.1) g.2.2.1 g.2.2.2.1 tag cls p.1 p.2.1 p.2.2.1 p.2.2.2.2.2.2) &&
  (Nat.beq p.2.2.2.1 0 ||
    pieceOk S g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.1 g.2.2.2.2.2.1 tag cls p.2.2.2.1 p.2.2.2.2.1
      p.2.2.2.2.2.1 p.2.2.2.2.2.2) &&
  (bif Nat.beq p.1 0 then
      (bif Nat.beq p.2.2.2.1 0 then Nat.beq p.2.2.2.2.2.2 0
       else !Nat.beq g.2.2.2.2.2.2.2.1 0 &&
         Nat.ble p.2.2.2.2.2.1 (Nat.sub g.2.2.2.2.2.2.2.1 1))
    else bif Nat.beq p.2.2.2.1 0 then
      !Nat.beq g.2.2.2.2.2.2.1 0 && Nat.ble (Nat.sub g.2.2.2.2.2.2.1 1) p.2.1
    else Nat.ble (Nat.mul p.2.2.2.2.1 R) (Nat.add (Nat.mul p.2.1 R) (Nat.mul U0 Q)) &&
      Nat.ble (Nat.mul p.2.2.2.2.2.1 R) (Nat.add (Nat.mul p.2.2.1 R) (Nat.mul U0 Q)))

/-- Consecutive pairs are sorted along both lines. -/
def tchain : List TPr → Bool
  | p :: q :: t => Nat.ble p.2.2.1 q.2.1 && Nat.ble p.2.2.2.2.2.1 q.2.2.2.2.1 && tchain (q :: t)
  | _ => true

/-- The certified intervals and `τ` points of a T-group. -/
def tgrpEnds (Q R x0 x1 y0 y1 U0 U1 : ℕ) (g : TGrp) : Bool :=
  Nat.ble g.2.2.1 g.2.2.2.1 && Nat.ble g.2.2.2.2.1 g.2.2.2.2.2.1 &&
  admAll Q R x0 x1 y0 y1 U0 U1 (kU g.1) (lx g.1 (lnU Q g.1 g.2.1) g.2.2.1)
    (ly g.1 (lnU Q g.1 g.2.1) g.2.2.1) &&
  admAll Q R x0 x1 y0 y1 U0 U1 (kU g.1) (lx g.1 (lnU Q g.1 g.2.1) g.2.2.2.1)
    (ly g.1 (lnU Q g.1 g.2.1) g.2.2.2.1) &&
  admAll Q R x0 x1 y0 y1 U0 U1 (kD g.1) (lx g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.1)
    (ly g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.1) &&
  admAll Q R x0 x1 y0 y1 U0 U1 (kD g.1) (lx g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.2.1)
    (ly g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.2.1) &&
  (Nat.beq g.2.2.2.2.2.2.1 0 ||
    admK Q R x0 x1 y0 y1 U0 U1 (lx g.1 (lnU Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.1 1))
      (ly g.1 (lnU Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.1 1)) (kU g.1)) &&
  (Nat.beq g.2.2.2.2.2.2.2.1 0 ||
    admK Q R x0 x1 y0 y1 U0 U1 (lx g.1 (lnD Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.2.1 1))
      (ly g.1 (lnD Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.2.1 1)) (kD g.1))

/-- A T-group (index `gi`, tag `2 gi + 1`). -/
def tgrpOk (S Q R x0 x1 y0 y1 U0 U1 : ℕ) (cls : List (SegE × ℕ)) (gi : ℕ) (g : TGrp) : Bool :=
  tgrpEnds Q R x0 x1 y0 y1 U0 U1 g &&
  g.2.2.2.2.2.2.2.2.all (tprOk S Q R U0 cls (Nat.add (Nat.mul 2 gi) 1) g) &&
  tchain g.2.2.2.2.2.2.2.2

def tgrpsOk (S Q R x0 x1 y0 y1 U0 U1 : ℕ) (cls : List (SegE × ℕ)) : ℕ → List TGrp → Bool
  | _, [] => true
  | gi, g :: gs => tgrpOk S Q R x0 x1 y0 y1 U0 U1 cls gi g && tgrpsOk S Q R x0 x1 y0 y1 U0 U1 cls (gi + 1) gs

def tgrpVal (g : TGrp) : ℕ := (g.2.2.2.2.2.2.2.2.map fun p => p.2.2.2.2.2.2).sum

/-- **The piece certificate check** (`cls` = the decoded claims). -/
def pcOk (S Q R x0 x1 y0 y1 U0 U1 : ℕ) (cls : List (SegE × ℕ)) (sb : List SBlk)
    (tg : List TGrp) : Bool :=
  sb.all (sblkOk Q R x0 x1 y0 y1 U0 U1) && cls.all (sclOk S sb) &&
    tgrpsOk S Q R x0 x1 y0 y1 U0 U1 cls 0 tg

/-- **The piece value** `Lp` (units of `1/W`). -/
def pcVal (S : ℕ) (cls : List (SegE × ℕ)) (sb : List SBlk) (tg : List TGrp) : ℕ :=
  (cls.map (sclVal S sb)).sum + (tg.map tgrpVal).sum


/-! ## 4.  Soundness of the piece certificate -/

section soundPC

variable {D S Q R x0 x1 y0 y1 U0 U1 : ℕ}

lemma getD_mem_of_tag {cls : List (SegE × ℕ)} {i tag : ℕ} (ht : 0 < tag)
    (h : (cls.getD i ec0).2 = tag) : cls.getD i ec0 ∈ cls := by
  by_cases hi : i < cls.length
  · rw [List.getD_eq_getElem _ _ hi]; exact List.getElem_mem hi
  · rw [List.getD_eq_default _ _ (by omega)] at h ⊢
    simp [ec0] at h; omega

lemma getD_mem_sb {sb : List SBlk} {i : ℕ} (hi : i < sb.length) : sb.getD i sb0 ∈ sb := by
  rw [List.getD_eq_getElem _ _ hi]; exact List.getElem_mem hi

lemma tag_unique {cls : List (SegE × ℕ)} (hnd : (cls.map Prod.fst).Nodup) {e : SegE} {t t' : ℕ}
    (h1 : (e, t) ∈ cls) (h2 : (e, t') ∈ cls) : t = t' := by
  have := List.inj_on_of_nodup_map hnd h1 h2 rfl
  simpa using congrArg Prod.snd this

lemma onLine_line {S dir L1 L2 : ℕ} {e : SegE} (h1 : onLine S dir L1 e = true)
    (h2 : onLine S dir L2 e = true) : L1 = L2 := by
  by_cases hd : dir = 0
  · subst hd
    simp only [onLine, beq_rfl, cond_true, Bool.and_eq_true, Nat.beq_eq, nat_mul_eq] at h1 h2
    omega
  · simp only [onLine, beq_ne hd, cond_false, Bool.and_eq_true, Nat.beq_eq, nat_mul_eq] at h1 h2
    omega

lemma lnU_ne_lnD {Q : ℕ} (hQ : 0 < Q) (dir K : ℕ) : lnU Q dir K ≠ lnD Q dir K := by
  by_cases hd : dir = 0
  · subst hd; simp only [lnU, lnD, beq_rfl, cond_true, nat_add_eq]; omega
  · simp only [lnU, lnD, beq_ne hd, cond_false, nat_add_eq]; omega

/-- An empty part is a part. -/
lemma empty_part_ok (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) {dir L : ℕ} {e : SegE}
    (he : onLine S dir L e = true) (Qs : Set (ℝ × ℝ)) :
    PartOK D Qs (e, (slo S dir e : ℝ) / Q, (slo S dir e : ℝ) / Q) := by
  obtain ⟨hax, hlo, hhi, hlt, _⟩ := onLine_spec (K := L) hD hS hQD he
  have hQ : (0 : ℝ) < Q := by rw [hQD]; exact_mod_cast Nat.mul_pos hD hS
  refine ⟨hax, by rw [hlo], le_rfl, ?_, fun t h1 h2 => absurd (lt_trans h1 h2) (lt_irrefl _)⟩
  rw [hhi]; exact div_le_div_of_nonneg_right (by exact_mod_cast hlt.le) hQ.le

/-- **A certified piece of a line is a part**: inside a segment on the line, inside an interval
certified at its ends for the conditions `≠ kp`, with the condition `kp` holding on the interior. -/
lemma line_part_ok (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {kp dir L A B : ℕ}
    (hA : admAll Q R x0 x1 y0 y1 U0 U1 kp (lx dir L A) (ly dir L A) = true)
    (hB : admAll Q R x0 x1 y0 y1 U0 U1 kp (lx dir L B) (ly dir L B) = true)
    {e : SegE} (he : onLine S dir L e = true) {l h : ℝ}
    (hl0 : (slo S dir e : ℝ) / Q ≤ l) (hlA : (A : ℝ) / Q ≤ l) (hlh : l ≤ h)
    (hh0 : h ≤ (shi S dir e : ℝ) / Q) (hhB : h ≤ (B : ℝ) / Q)
    (hsing : ∀ t, l < t → t < h → ∀ k : Fin 4, k.val = kp →
      gval k ((lpt Q dir L t).1 - c.1) ((lpt Q dir L t).2 - c.2) u ≤ 0) :
    PartOK D (sq c (2 * Real.arctan u) 1) (e, l, h) := by
  have hQ : 0 < Q := by rw [hQD]; exact Nat.mul_pos hD hS
  obtain ⟨hax, hlo, hhi, _, hlp⟩ := onLine_spec (K := L) hD hS hQD he
  refine ⟨hax, by rw [hlo]; exact hl0, hlh, by rw [hhi]; exact hh0, fun t ht1 ht2 => ?_⟩
  simp only
  rw [hlp, mem_sq_iff_gval]
  intro k
  by_cases hk : k.val = kp
  · exact hsing t ht1 ht2 k hk
  · exact certified_interval hQ hR P hA hB (by linarith) (by linarith) k hk

/-- The mass of a part on a line. -/
lemma pval_line (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) {dir L : ℕ} {e : SegE}
    (he : onLine S dir L e = true) (l h : ℝ) :
    pval D (e, l, h) = (e.2.2.2.2 : ℝ) * Q * (h - l) / ((shi S dir e : ℝ) - slo S dir e) := by
  obtain ⟨_, hlo, hhi, hlt, _⟩ := onLine_spec (K := L) hD hS hQD he
  have hQ : (0 : ℝ) < Q := by rw [hQD]; exact_mod_cast Nat.mul_pos hD hS
  have hl : (0 : ℝ) < (shi S dir e : ℝ) - slo S dir e := by
    have : (slo S dir e : ℝ) < shi S dir e := by exact_mod_cast hlt
    linarith
  simp only [pval]
  rw [hlo, hhi, ← sub_div]
  field_simp

/-! ### S-blocks -/

open Classical in
/-- The part an S-block claim certifies. -/
noncomputable def sPart (S Q : ℕ) (sb : List SBlk) (ec : SegE × ℕ) : SegE × ℝ × ℝ :=
  if max (slo S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1) (sb.getD (Nat.div ec.2 2) sb0).2.2.1
      ≤ min (shi S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1) (sb.getD (Nat.div ec.2 2) sb0).2.2.2 then
    (ec.1, (max (slo S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1)
        (sb.getD (Nat.div ec.2 2) sb0).2.2.1 : ℝ) / Q,
      (min (shi S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1) (sb.getD (Nat.div ec.2 2) sb0).2.2.2 : ℝ) / Q)
  else (ec.1, (slo S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1 : ℝ) / Q,
    (slo S (sb.getD (Nat.div ec.2 2) sb0).1 ec.1 : ℝ) / Q)

lemma sPart_fst (S Q : ℕ) (sb : List SBlk) (ec : SegE × ℕ) : (sPart S Q sb ec).1 = ec.1 := by
  unfold sPart; split_ifs <;> rfl

lemma sPart_ok (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {sb : List SBlk}
    (hsb : sb.all (sblkOk Q R x0 x1 y0 y1 U0 U1) = true) {ec : SegE × ℕ}
    (hev : Nat.mod ec.2 2 = 0) (hcl : sclOk S sb ec = true) :
    PartOK D (sq c (2 * Real.arctan u) 1) (sPart S Q sb ec) ∧
      (sclVal S sb ec : ℝ) ≤ pval D (sPart S Q sb ec) := by
  have hQ : (0 : ℝ) < Q := by rw [hQD]; exact_mod_cast Nat.mul_pos hD hS
  have hb0 : Nat.beq (Nat.mod ec.2 2) 0 = true := by rw [hev]; rfl
  simp only [sclOk, hb0, cond_true, Bool.and_eq_true, Nat.blt_eq] at hcl
  obtain ⟨hlt, hon⟩ := hcl
  set b := sb.getD (Nat.div ec.2 2) sb0 with hbdef
  have hbm : b ∈ sb := getD_mem_sb hlt
  have hbok := List.all_eq_true.mp hsb b hbm
  simp only [sblkOk, Bool.and_eq_true, Nat.ble_eq] at hbok
  obtain ⟨⟨hAB, hA⟩, hB⟩ := hbok
  have hsv : sclVal S sb ec = sval S b ec.1 := by simp only [sclVal, hb0, cond_true, hbdef]
  obtain ⟨_, _, _, hslt, _⟩ := onLine_spec (K := b.2.1) hD hS hQD hon
  have hlenr : (0 : ℝ) < (shi S b.1 ec.1 : ℝ) - slo S b.1 ec.1 := by
    have : (slo S b.1 ec.1 : ℝ) < shi S b.1 ec.1 := by exact_mod_cast hslt
    linarith
  -- the value
  have hval : (sval S b ec.1 : ℝ) ≤ (ec.1.2.2.2.2 : ℝ) *
      ((min (shi S b.1 ec.1) b.2.2.2 - max (slo S b.1 ec.1) b.2.2.1 : ℕ) : ℝ) /
      ((shi S b.1 ec.1 : ℝ) - slo S b.1 ec.1) := by
    simp only [sval, nat_mul_eq, nat_sub_eq, nat_div_eq]
    have := Nat.cast_div_le (α := ℝ) (m := ec.1.2.2.2.2 * (min (shi S b.1 ec.1) b.2.2.2 -
      max (slo S b.1 ec.1) b.2.2.1)) (n := shi S b.1 ec.1 - slo S b.1 ec.1)
    rw [Nat.cast_mul, Nat.cast_sub hslt.le] at this
    exact this
  unfold sPart
  rw [← hbdef]
  split_ifs with hle
  · refine ⟨line_part_ok hD hS hQD hR P hA hB hon ?_ ?_ ?_ ?_ ?_ ?_, ?_⟩
    · exact div_le_div_of_nonneg_right (by exact_mod_cast le_max_left _ _) hQ.le
    · exact div_le_div_of_nonneg_right (by exact_mod_cast le_max_right _ _) hQ.le
    · exact div_le_div_of_nonneg_right (by exact_mod_cast hle) hQ.le
    · exact div_le_div_of_nonneg_right (by exact_mod_cast min_le_left _ _) hQ.le
    · exact div_le_div_of_nonneg_right (by exact_mod_cast min_le_right _ _) hQ.le
    · intro t _ _ k hk; exact absurd hk (by omega)
    · rw [hsv, pval_line hD hS hQD hon]
      refine le_trans hval (le_of_eq ?_)
      rw [Nat.cast_sub hle]
      push_cast
      field_simp
  · refine ⟨empty_part_ok hD hS hQD hon _, ?_⟩
    rw [hsv, pval_line hD hS hQD hon]
    refine le_trans hval (le_of_eq ?_)
    have : min (shi S b.1 ec.1) b.2.2.2 - max (slo S b.1 ec.1) b.2.2.1 = 0 := by omega
    rw [this]; simp

/-! ### T-groups -/

/-- The singular violation polynomials of a germ pair along its two lines. -/
noncomputable def gUp (Q dir K : ℕ) (c : ℝ × ℝ) (u t : ℝ) : ℝ :=
  gval (kf (kU dir)) ((lpt Q dir (lnU Q dir K) t).1 - c.1) ((lpt Q dir (lnU Q dir K) t).2 - c.2) u
noncomputable def gDn (Q dir K : ℕ) (c : ℝ × ℝ) (u t : ℝ) : ℝ :=
  gval (kf (kD dir)) ((lpt Q dir (lnD Q dir K) t).1 - c.1) ((lpt Q dir (lnD Q dir K) t).2 - c.2) u

lemma gUp_gDn {Q : ℕ} (hQ : 0 < Q) (dir K : ℕ) (c : ℝ × ℝ) (u t t' : ℝ) :
    gUp Q dir K c u t + gDn Q dir K c u t' = 4 * u * (t' - t - u) := by
  have hQr : (Q : ℝ) ≠ 0 := by exact_mod_cast hQ.ne'
  have e1 : ((K + Q : ℕ) : ℝ) / Q = (K : ℝ) / Q + 1 := by push_cast; field_simp
  by_cases hd : dir = 0
  · subst hd
    simp only [gUp, gDn, kU, kD, lnU, lnD, beq_rfl, cond_true, lpt, if_true, nat_add_eq, e1]
    exact gval_pairV _ _ _ _ _ _
  · simp only [gUp, gDn, kU, kD, lnU, lnD, beq_ne hd, cond_false, lpt, hd, if_false, nat_add_eq,
      e1]
    exact gval_pairH _ _ _ _ _ _

lemma kf_kU_val (dir : ℕ) : (kf (kU dir)).val = kU dir := by
  by_cases hd : dir = 0
  · subst hd; rfl
  · simp only [kU, beq_ne hd, cond_false]; rfl
lemma kf_kD_val (dir : ℕ) : (kf (kD dir)).val = kD dir := by
  by_cases hd : dir = 0
  · subst hd; rfl
  · simp only [kD, beq_ne hd, cond_false]; rfl

/-- The singular condition along the up line is `gUp ≤ 0`. -/
lemma sing_up {Q dir K : ℕ} {c : ℝ × ℝ} {u t : ℝ} (h : gUp Q dir K c u t ≤ 0) (k : Fin 4)
    (hk : k.val = kU dir) :
    gval k ((lpt Q dir (lnU Q dir K) t).1 - c.1) ((lpt Q dir (lnU Q dir K) t).2 - c.2) u ≤ 0 := by
  have : k = kf (kU dir) := by rw [← hk, kf_val]
  rw [this]; exact h
lemma sing_dn {Q dir K : ℕ} {c : ℝ × ℝ} {u t : ℝ} (h : gDn Q dir K c u t ≤ 0) (k : Fin 4)
    (hk : k.val = kD dir) :
    gval k ((lpt Q dir (lnD Q dir K) t).1 - c.1) ((lpt Q dir (lnD Q dir K) t).2 - c.2) u ≤ 0 := by
  have : k = kf (kD dir) := by rw [← hk, kf_val]
  rw [this]; exact h

/-- The facts a T-group's ends give. -/
structure GEnds (Q R x0 x1 y0 y1 U0 U1 : ℕ) (g : TGrp) : Prop where
  hAu : admAll Q R x0 x1 y0 y1 U0 U1 (kU g.1) (lx g.1 (lnU Q g.1 g.2.1) g.2.2.1)
    (ly g.1 (lnU Q g.1 g.2.1) g.2.2.1) = true
  hBu : admAll Q R x0 x1 y0 y1 U0 U1 (kU g.1) (lx g.1 (lnU Q g.1 g.2.1) g.2.2.2.1)
    (ly g.1 (lnU Q g.1 g.2.1) g.2.2.2.1) = true
  hAd : admAll Q R x0 x1 y0 y1 U0 U1 (kD g.1) (lx g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.1)
    (ly g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.1) = true
  hBd : admAll Q R x0 x1 y0 y1 U0 U1 (kD g.1) (lx g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.2.1)
    (ly g.1 (lnD Q g.1 g.2.1) g.2.2.2.2.2.1) = true
  htu : g.2.2.2.2.2.2.1 ≠ 0 →
    admK Q R x0 x1 y0 y1 U0 U1 (lx g.1 (lnU Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.1 1))
      (ly g.1 (lnU Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.1 1)) (kU g.1) = true
  htd : g.2.2.2.2.2.2.2.1 ≠ 0 →
    admK Q R x0 x1 y0 y1 U0 U1 (lx g.1 (lnD Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.2.1 1))
      (ly g.1 (lnD Q g.1 g.2.1) (Nat.sub g.2.2.2.2.2.2.2.1 1)) (kD g.1) = true

lemma gends_of {g : TGrp} (h : tgrpEnds Q R x0 x1 y0 y1 U0 U1 g = true) :
    GEnds Q R x0 x1 y0 y1 U0 U1 g := by
  simp only [tgrpEnds, Bool.and_eq_true, Bool.or_eq_true, Nat.beq_eq] at h
  obtain ⟨⟨⟨⟨⟨⟨⟨_, _⟩, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩ := h
  exact ⟨h1, h2, h3, h4, fun hn => h5.resolve_left hn, fun hn => h6.resolve_left hn⟩

/-- `gUp` at a `τ` certified by `admK`. -/
lemma gUp_tau (hQ : 0 < Q) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ} (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u)
    {dir K τ : ℕ}
    (h : admK Q R x0 x1 y0 y1 U0 U1 (lx dir (lnU Q dir K) τ) (ly dir (lnU Q dir K) τ) (kU dir)
      = true) :
    gUp Q dir K c u ((τ : ℝ) / Q) ≤ 0 := by
  have := admK_sound hQ hR P.hU01 P.hU1 h P.hx0 P.hx1 P.hy0 P.hy1 P.hu0 P.hu1 P.hwx P.hwy
  have e := Gv_lpt Q dir (lnU Q dir K) τ (kf (kU dir)) c u
  rw [kf_kU_val] at e
  rw [e] at this
  exact this
lemma gDn_tau (hQ : 0 < Q) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ} (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u)
    {dir K τ : ℕ}
    (h : admK Q R x0 x1 y0 y1 U0 U1 (lx dir (lnD Q dir K) τ) (ly dir (lnD Q dir K) τ) (kD dir)
      = true) :
    gDn Q dir K c u ((τ : ℝ) / Q) ≤ 0 := by
  have := admK_sound hQ hR P.hU01 P.hU1 h P.hx0 P.hx1 P.hy0 P.hy1 P.hu0 P.hu1 P.hwx P.hwy
  have e := Gv_lpt Q dir (lnD Q dir K) τ (kf (kD dir)) c u
  rw [kf_kD_val] at e
  rw [e] at this
  exact this

lemma u_nonneg {c : ℝ × ℝ} {u : ℝ} (hR : 0 < R) (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) : 0 ≤ u :=
  le_trans (div_nonneg (Nat.cast_nonneg _) (by exact_mod_cast hR.le)) P.hu0

/-- The facts `pieceOk` gives. -/
lemma pieceOk_spec {dir L lo hi tag j1 a b m : ℕ} {cls : List (SegE × ℕ)} (htag : 0 < tag)
    (h : pieceOk S dir L lo hi tag cls j1 a b m = true) :
    ((cls.getD (Nat.sub j1 1) ec0).1, tag) ∈ cls ∧
      onLine S dir L (cls.getD (Nat.sub j1 1) ec0).1 = true ∧
      slo S dir (cls.getD (Nat.sub j1 1) ec0).1 ≤ a ∧ b ≤ shi S dir (cls.getD (Nat.sub j1 1) ec0).1 ∧
      lo ≤ a ∧ b ≤ hi ∧
      m * (shi S dir (cls.getD (Nat.sub j1 1) ec0).1 - slo S dir (cls.getD (Nat.sub j1 1) ec0).1)
        ≤ (cls.getD (Nat.sub j1 1) ec0).1.2.2.2.2 * (b - a) := by
  simp only [pieceOk, Bool.and_eq_true, Nat.beq_eq, Nat.ble_eq, nat_mul_eq, nat_sub_eq] at h
  obtain ⟨⟨⟨⟨⟨⟨ht, hon⟩, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩ := h
  have hm := getD_mem_of_tag htag ht
  have hm' : ((cls.getD (j1 - 1) ec0).1, tag) ∈ cls := by
    rw [← ht]; exact hm
  exact ⟨hm', hon, h1, h2, h3, h4, h5⟩

/-- The mass of a piece of a line bounds `m`. -/
lemma piece_mass (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) {dir L : ℕ} {e : SegE}
    (he : onLine S dir L e = true) {a b m : ℕ} (hab : a ≤ b)
    (hm : m * (shi S dir e - slo S dir e) ≤ e.2.2.2.2 * (b - a)) :
    (m : ℝ) ≤ (e.2.2.2.2 : ℝ) * Q / ((shi S dir e : ℝ) - slo S dir e) *
      ((b : ℝ) / Q - (a : ℝ) / Q) := by
  obtain ⟨_, _, _, hlt, _⟩ := onLine_spec (K := L) hD hS hQD he
  have hQ : (0 : ℝ) < Q := by rw [hQD]; exact_mod_cast Nat.mul_pos hD hS
  have hl : (0 : ℝ) < (shi S dir e : ℝ) - slo S dir e := by
    have : (slo S dir e : ℝ) < shi S dir e := by exact_mod_cast hlt
    linarith
  have hm' : (m : ℝ) * ((shi S dir e : ℝ) - slo S dir e) ≤ (e.2.2.2.2 : ℝ) * ((b : ℝ) - a) := by
    have := (Nat.cast_le (α := ℝ)).mpr hm
    push_cast [Nat.cast_sub hlt.le, Nat.cast_sub hab] at this
    exact this
  have e1 : (e.2.2.2.2 : ℝ) * Q / ((shi S dir e : ℝ) - slo S dir e) * ((b : ℝ) / Q - (a : ℝ) / Q)
      = (e.2.2.2.2 : ℝ) * ((b : ℝ) - a) / ((shi S dir e : ℝ) - slo S dir e) := by
    field_simp
  rw [e1, le_div_iff₀ hl]
  exact hm'

/-- The parts of one pair at a pose. -/
lemma pair_parts (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {cls : List (SegE × ℕ)} {tag : ℕ} (htag : 0 < tag)
    {g : TGrp} (hg : GEnds Q R x0 x1 y0 y1 U0 U1 g) {p : TPr}
    (hp : tprOk S Q R U0 cls tag g p = true) :
    ∃ L : List (SegE × ℝ × ℝ),
      (∀ x ∈ L, PartOK D (sq c (2 * Real.arctan u) 1) x ∧ (x.1, tag) ∈ cls ∧
        ((onLine S g.1 (lnU Q g.1 g.2.1) x.1 = true ∧ (p.2.1 : ℝ) / Q ≤ x.2.1 ∧
            x.2.2 ≤ (p.2.2.1 : ℝ) / Q) ∨
         (onLine S g.1 (lnD Q g.1 g.2.1) x.1 = true ∧ (p.2.2.2.2.1 : ℝ) / Q ≤ x.2.1 ∧
            x.2.2 ≤ (p.2.2.2.2.2.1 : ℝ) / Q))) ∧
      L.Pairwise (fun x y => x.1 ≠ y.1) ∧
      (p.2.2.2.2.2.2 : ℝ) ≤ (L.map (pval D)).sum := by
  have hQ : 0 < Q := by rw [hQD]; exact Nat.mul_pos hD hS
  have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
  obtain ⟨ju1, a, b, jd1, cc, d, m⟩ := p
  obtain ⟨dir, K, Au, Bu, Ad, Bd, tu1, td1, prs⟩ := g
  obtain ⟨hAu, hBu, hAd, hBd, htu, htd⟩ := hg
  simp only at hAu hBu hAd hBd htu htd ⊢
  simp only [tprOk, Bool.and_eq_true, Bool.or_eq_true, Nat.beq_eq, Nat.ble_eq] at hp
  obtain ⟨⟨⟨⟨hab, hcd⟩, hup⟩, hdn⟩, hcase⟩ := hp
  have hUD := gUp_gDn hQ dir K c u
  -- the up part, when there is one
  have upPart : ju1 ≠ 0 → ∀ {l h : ℝ}, (a : ℝ) / Q ≤ l → l ≤ h → h ≤ (b : ℝ) / Q →
      (∀ t, l < t → t < h → gUp Q dir K c u t ≤ 0) →
      PartOK D (sq c (2 * Real.arctan u) 1) ((cls.getD (Nat.sub ju1 1) ec0).1, l, h) := by
    intro hj l h hl hlh hh hs
    obtain ⟨_, hon, h1, h2, h3, h4, _⟩ := pieceOk_spec htag (hup.resolve_left hj)
    refine line_part_ok hD hS hQD hR P hAu hBu hon ?_ ?_ hlh ?_ ?_ fun t ht1 ht2 k hk =>
      sing_up (hs t ht1 ht2) k hk
    · exact le_trans (div_le_div_of_nonneg_right (by exact_mod_cast h1) hQr.le) hl
    · exact le_trans (div_le_div_of_nonneg_right (by exact_mod_cast h3) hQr.le) hl
    · exact le_trans hh (div_le_div_of_nonneg_right (by exact_mod_cast h2) hQr.le)
    · exact le_trans hh (div_le_div_of_nonneg_right (by exact_mod_cast h4) hQr.le)
  have dnPart : jd1 ≠ 0 → ∀ {l h : ℝ}, (cc : ℝ) / Q ≤ l → l ≤ h → h ≤ (d : ℝ) / Q →
      (∀ t, l < t → t < h → gDn Q dir K c u t ≤ 0) →
      PartOK D (sq c (2 * Real.arctan u) 1) ((cls.getD (Nat.sub jd1 1) ec0).1, l, h) := by
    intro hj l h hl hlh hh hs
    obtain ⟨_, hon, h1, h2, h3, h4, _⟩ := pieceOk_spec htag (hdn.resolve_left hj)
    refine line_part_ok hD hS hQD hR P hAd hBd hon ?_ ?_ hlh ?_ ?_ fun t ht1 ht2 k hk =>
      sing_dn (hs t ht1 ht2) k hk
    · exact le_trans (div_le_div_of_nonneg_right (by exact_mod_cast h1) hQr.le) hl
    · exact le_trans (div_le_div_of_nonneg_right (by exact_mod_cast h3) hQr.le) hl
    · exact le_trans hh (div_le_div_of_nonneg_right (by exact_mod_cast h2) hQr.le)
    · exact le_trans hh (div_le_div_of_nonneg_right (by exact_mod_cast h4) hQr.le)
  have hab' : (a : ℝ) / Q ≤ (b : ℝ) / Q := div_le_div_of_nonneg_right (by exact_mod_cast hab) hQr.le
  have hcd' : (cc : ℝ) / Q ≤ (d : ℝ) / Q :=
    div_le_div_of_nonneg_right (by exact_mod_cast hcd) hQr.le
  by_cases hj : ju1 = 0
  · have hj0 : Nat.beq ju1 0 = true := by rw [hj]; rfl
    by_cases hk : jd1 = 0
    · -- no piece: `m = 0`
      have hk0 : Nat.beq jd1 0 = true := by rw [hk]; rfl
      simp only [hj0, hk0, cond_true, Nat.beq_eq] at hcase
      refine ⟨[], by simp, List.Pairwise.nil, by simp [hcase]⟩
    · -- down piece only, below `τd`
      have hk0 : Nat.beq jd1 0 = false := beq_ne hk
      simp only [hj0, hk0, cond_true, cond_false, Bool.and_eq_true, Bool.not_eq_true',
        Nat.ble_eq] at hcase
      obtain ⟨htd0, hdτ⟩ := hcase
      have htd' : td1 ≠ 0 := by intro h0; rw [h0] at htd0; exact absurd htd0 (by decide)
      have hτ := gDn_tau hQ hR P (htd htd')
      obtain ⟨hmem2, hon, _, _, _, _, hm⟩ := pieceOk_spec htag (hdn.resolve_left hk)
      refine ⟨[((cls.getD (Nat.sub jd1 1) ec0).1, (cc : ℝ) / Q, (d : ℝ) / Q)], ?_, by simp, ?_⟩
      · intro x hx
        rw [List.mem_singleton] at hx
        subst hx
        refine ⟨dnPart hk le_rfl hcd' le_rfl fun t _ ht => ?_, hmem2, Or.inr ⟨hon, le_rfl, le_rfl⟩⟩
        have hle : t ≤ ((Nat.sub td1 1 : ℕ) : ℝ) / Q := le_trans ht.le
          (div_le_div_of_nonneg_right (by exact_mod_cast hdτ) hQr.le)
        have e1 := hUD 0 t
        have e2 := hUD 0 (((Nat.sub td1 1 : ℕ) : ℝ) / Q)
        have hu := u_nonneg hR P
        nlinarith
      · simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
        rw [pval_line hD hS hQD hon]
        have := piece_mass hD hS hQD hon hcd hm
        refine le_trans this (le_of_eq ?_)
        field_simp
  · have hj0 : Nat.beq ju1 0 = false := beq_ne hj
    obtain ⟨hmemU2, honU, _, _, _, _, hmU⟩ := pieceOk_spec htag (hup.resolve_left hj)
    by_cases hk : jd1 = 0
    · -- up piece only, above `τu`
      have hk0 : Nat.beq jd1 0 = true := by rw [hk]; rfl
      simp only [hj0, hk0, cond_true, cond_false, Bool.and_eq_true, Bool.not_eq_true',
        Nat.ble_eq] at hcase
      obtain ⟨htu0, haτ⟩ := hcase
      have htu' : tu1 ≠ 0 := by intro h0; rw [h0] at htu0; exact absurd htu0 (by decide)
      have hτ := gUp_tau hQ hR P (htu htu')
      refine ⟨[((cls.getD (Nat.sub ju1 1) ec0).1, (a : ℝ) / Q, (b : ℝ) / Q)], ?_, by simp, ?_⟩
      · intro x hx
        rw [List.mem_singleton] at hx
        subst hx
        refine ⟨upPart hj le_rfl hab' le_rfl fun t ht _ => ?_, hmemU2,
          Or.inl ⟨honU, le_rfl, le_rfl⟩⟩
        have hle : ((Nat.sub tu1 1 : ℕ) : ℝ) / Q ≤ t := le_trans
          (div_le_div_of_nonneg_right (by exact_mod_cast haτ) hQr.le) ht.le
        have e1 := hUD t 0
        have e2 := hUD (((Nat.sub tu1 1 : ℕ) : ℝ) / Q) 0
        have hu := u_nonneg hR P
        nlinarith
      · simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
        rw [pval_line hD hS hQD honU]
        have := piece_mass hD hS hQD honU hab hmU
        refine le_trans this (le_of_eq ?_)
        field_simp
    · -- a genuine pair
      have hk0 : Nat.beq jd1 0 = false := beq_ne hk
      simp only [hj0, hk0, cond_false, Bool.and_eq_true, Nat.ble_eq, nat_mul_eq,
        nat_add_eq] at hcase
      obtain ⟨hca, hdb⟩ := hcase
      obtain ⟨hmemD2, honD, _, _, _, _, hmD⟩ := pieceOk_spec htag (hdn.resolve_left hk)
      have hRr : (0 : ℝ) < R := by exact_mod_cast hR
      have hca' : (cc : ℝ) / Q ≤ (a : ℝ) / Q + (U0 : ℝ) / R := by
        have : ((cc * R : ℕ) : ℝ) ≤ ((a * R + U0 * Q : ℕ) : ℝ) := by exact_mod_cast hca
        push_cast at this
        rw [div_add_div _ _ hQr.ne' hRr.ne', div_le_div_iff₀ hQr (by positivity)]
        nlinarith
      have hdb' : (d : ℝ) / Q ≤ (b : ℝ) / Q + (U0 : ℝ) / R := by
        have : ((d * R : ℕ) : ℝ) ≤ ((b * R + U0 * Q : ℕ) : ℝ) := by exact_mod_cast hdb
        push_cast at this
        rw [div_add_div _ _ hQr.ne' hRr.ne', div_le_div_iff₀ hQr (by positivity)]
        nlinarith
      set eU := (cls.getD (Nat.sub ju1 1) ec0).1
      set eD := (cls.getD (Nat.sub jd1 1) ec0).1
      obtain ⟨_, _, _, hltU, _⟩ := onLine_spec (K := lnU Q dir K) hD hS hQD honU
      obtain ⟨_, _, _, hltD, _⟩ := onLine_spec (K := lnD Q dir K) hD hS hQD honD
      have hlU : (0 : ℝ) < (shi S dir eU : ℝ) - slo S dir eU := by
        have : (slo S dir eU : ℝ) < shi S dir eU := by exact_mod_cast hltU
        linarith
      have hlD : (0 : ℝ) < (shi S dir eD : ℝ) - slo S dir eD := by
        have : (slo S dir eD : ℝ) < shi S dir eD := by exact_mod_cast hltD
        linarith
      obtain ⟨a', d', ha1, ha2, hd1, hd2, hsU, hsD, hval⟩ :=
        pair_capture (gU := gUp Q dir K c u) (gD := gDn Q dir K c u) (u0 := (U0 : ℝ) / R)
          (fun t t' => hUD t t') (div_nonneg (Nat.cast_nonneg _) hRr.le) P.hu0 hab' hcd' hca' hdb'
          (div_nonneg (mul_nonneg (Nat.cast_nonneg _) hQr.le) hlU.le)
          (div_nonneg (mul_nonneg (Nat.cast_nonneg _) hQr.le) hlD.le)
          (piece_mass hD hS hQD honU hab hmU) (piece_mass hD hS hQD honD hcd hmD)
      have hne : eU ≠ eD := by
        intro h
        rw [h] at honU
        exact lnU_ne_lnD hQ dir K (onLine_line honU honD)
      refine ⟨[(eU, a', (b : ℝ) / Q), (eD, (cc : ℝ) / Q, d')], ?_, ?_, ?_⟩
      · intro x hx
        simp only [List.mem_cons, List.not_mem_nil, or_false] at hx
        rcases hx with rfl | rfl
        · exact ⟨upPart hj ha1 ha2 le_rfl hsU, hmemU2, Or.inl ⟨honU, ha1, le_rfl⟩⟩
        · exact ⟨dnPart hk le_rfl hd1 hd2 hsD, hmemD2, Or.inr ⟨honD, le_rfl, hd2⟩⟩
      · simp [hne]
      · simp only [List.map_cons, List.map_nil, List.sum_cons, List.sum_nil, add_zero]
        rw [pval_line hD hS hQD honU, pval_line hD hS hQD honD]
        refine le_trans hval (le_of_eq ?_)
        field_simp

end soundPC


section soundPC2

variable {D S Q R x0 x1 y0 y1 U0 U1 : ℕ}

/-- The relation the parts lemma needs. -/
def PRel (x y : SegE × ℝ × ℝ) : Prop := x.1 = y.1 → x.2.2 ≤ y.2.1

/-- Sorted pairs: every later pair starts after the earlier one ends, on both lines. -/
lemma tchain_bound : ∀ (p : TPr) (prs : List TPr), tchain (p :: prs) = true →
    (∀ q ∈ prs, q.2.1 ≤ q.2.2.1 ∧ q.2.2.2.2.1 ≤ q.2.2.2.2.2.1) →
    ∀ q ∈ prs, p.2.2.1 ≤ q.2.1 ∧ p.2.2.2.2.2.1 ≤ q.2.2.2.2.1
  | _, [], _, _ => by simp
  | p, q :: rest, h, hord => by
    simp only [tchain, Bool.and_eq_true, Nat.ble_eq] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    have ih := tchain_bound q rest h3 (fun r hr => hord r (List.mem_cons_of_mem _ hr))
    intro r hr
    rcases List.mem_cons.mp hr with rfl | hr
    · exact ⟨h1, h2⟩
    · obtain ⟨k1, k2⟩ := ih r hr
      obtain ⟨o1, o2⟩ := hord q List.mem_cons_self
      exact ⟨by omega, by omega⟩

lemma tchain_tail {p : TPr} {prs : List TPr} (h : tchain (p :: prs) = true) : tchain prs = true := by
  cases prs with
  | nil => rfl
  | cons q rest =>
    simp only [tchain, Bool.and_eq_true] at h
    exact h.2

/-- The parts of a group's pairs. -/
lemma group_parts (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {cls : List (SegE × ℕ)} {tag : ℕ} (htag : 0 < tag)
    {g : TGrp} (hg : GEnds Q R x0 x1 y0 y1 U0 U1 g) :
    ∀ (prs : List TPr), (∀ p ∈ prs, tprOk S Q R U0 cls tag g p = true) → tchain prs = true →
      ∀ α γ : ℕ, (∀ p ∈ prs, α ≤ p.2.1 ∧ γ ≤ p.2.2.2.2.1) →
      ∃ L : List (SegE × ℝ × ℝ),
        (∀ x ∈ L, PartOK D (sq c (2 * Real.arctan u) 1) x ∧ (x.1, tag) ∈ cls ∧
          ((onLine S g.1 (lnU Q g.1 g.2.1) x.1 = true ∧ (α : ℝ) / Q ≤ x.2.1) ∨
           (onLine S g.1 (lnD Q g.1 g.2.1) x.1 = true ∧ (γ : ℝ) / Q ≤ x.2.1))) ∧
        L.Pairwise PRel ∧
        (((prs.map fun p => p.2.2.2.2.2.2).sum : ℕ) : ℝ) ≤ (L.map (pval D)).sum
  | [], _, _, _, _, _ => ⟨[], by simp, List.Pairwise.nil, by simp⟩
  | p :: prs, hall, hch, α, γ, hαγ => by
    have hQ : 0 < Q := by rw [hQD]; exact Nat.mul_pos hD hS
    have hQr : (0 : ℝ) < Q := by exact_mod_cast hQ
    have hp := hall p List.mem_cons_self
    obtain ⟨Lp, hLp1, hLp2, hLp3⟩ := pair_parts hD hS hQD hR P htag hg hp
    have hord : ∀ q ∈ prs, q.2.1 ≤ q.2.2.1 ∧ q.2.2.2.2.1 ≤ q.2.2.2.2.2.1 := by
      intro q hq
      have := hall q (List.mem_cons_of_mem _ hq)
      simp only [tprOk, Bool.and_eq_true, Nat.ble_eq] at this
      exact ⟨this.1.1.1.1, this.1.1.1.2⟩
    have hb := tchain_bound p prs hch hord
    obtain ⟨L', hL1, hL2, hL3⟩ := group_parts hD hS hQD hR P htag hg prs
      (fun q hq => hall q (List.mem_cons_of_mem _ hq)) (tchain_tail hch) p.2.2.1 p.2.2.2.2.2.1 hb
    have hpab : p.2.1 ≤ p.2.2.1 ∧ p.2.2.2.2.1 ≤ p.2.2.2.2.2.1 := by
      simp only [tprOk, Bool.and_eq_true, Nat.ble_eq] at hp
      exact ⟨hp.1.1.1.1, hp.1.1.1.2⟩
    obtain ⟨hα, hγ⟩ := hαγ p List.mem_cons_self
    have cα : (α : ℝ) / Q ≤ (p.2.1 : ℝ) / Q := div_le_div_of_nonneg_right (by exact_mod_cast hα) hQr.le
    have cγ : (γ : ℝ) / Q ≤ (p.2.2.2.2.1 : ℝ) / Q :=
      div_le_div_of_nonneg_right (by exact_mod_cast hγ) hQr.le
    have cab : (p.2.1 : ℝ) / Q ≤ (p.2.2.1 : ℝ) / Q :=
      div_le_div_of_nonneg_right (by exact_mod_cast hpab.1) hQr.le
    have ccd : (p.2.2.2.2.1 : ℝ) / Q ≤ (p.2.2.2.2.2.1 : ℝ) / Q :=
      div_le_div_of_nonneg_right (by exact_mod_cast hpab.2) hQr.le
    refine ⟨Lp ++ L', ?_, ?_, ?_⟩
    · intro x hx
      rcases List.mem_append.mp hx with hx | hx
      · obtain ⟨h1, h2, h3⟩ := hLp1 x hx
        refine ⟨h1, h2, ?_⟩
        rcases h3 with ⟨ho, hl, _⟩ | ⟨ho, hl, _⟩
        · exact Or.inl ⟨ho, le_trans cα hl⟩
        · exact Or.inr ⟨ho, le_trans cγ hl⟩
      · obtain ⟨h1, h2, h3⟩ := hL1 x hx
        refine ⟨h1, h2, ?_⟩
        rcases h3 with ⟨ho, hl⟩ | ⟨ho, hl⟩
        · exact Or.inl ⟨ho, le_trans (le_trans cα cab) hl⟩
        · exact Or.inr ⟨ho, le_trans (le_trans cγ ccd) hl⟩
    · rw [List.pairwise_append]
      refine ⟨hLp2.imp fun h e => absurd e h, hL2, fun x hx y hy hxy => ?_⟩
      obtain ⟨_, _, hx3⟩ := hLp1 x hx
      obtain ⟨_, _, hy3⟩ := hL1 y hy
      rcases hx3 with ⟨hox, _, hhx⟩ | ⟨hox, _, hhx⟩ <;> rcases hy3 with ⟨hoy, hly⟩ | ⟨hoy, hly⟩
      · exact le_trans hhx hly
      · rw [hxy] at hox; exact absurd (onLine_line hox hoy) (lnU_ne_lnD hQ _ _)
      · rw [hxy] at hox; exact absurd (onLine_line hoy hox) (lnU_ne_lnD hQ _ _)
      · exact le_trans hhx hly
    · simp only [List.map_cons, List.sum_cons, List.map_append, List.sum_append]
      rw [Nat.cast_add]
      linarith

/-- The parts of all T-groups. -/
lemma groups_parts (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {cls : List (SegE × ℕ)}
    (hnd : (cls.map Prod.fst).Nodup) :
    ∀ (tg : List TGrp) (gi : ℕ), tgrpsOk S Q R x0 x1 y0 y1 U0 U1 cls gi tg = true →
      ∃ L : List (SegE × ℝ × ℝ),
        (∀ x ∈ L, PartOK D (sq c (2 * Real.arctan u) 1) x ∧
          ∃ j, gi ≤ j ∧ (x.1, 2 * j + 1) ∈ cls) ∧
        L.Pairwise PRel ∧ (((tg.map tgrpVal).sum : ℕ) : ℝ) ≤ (L.map (pval D)).sum
  | [], _, _ => ⟨[], by simp, List.Pairwise.nil, by simp⟩
  | g :: tg, gi, h => by
    simp only [tgrpsOk, Bool.and_eq_true] at h
    obtain ⟨hg, htg⟩ := h
    simp only [tgrpOk, Bool.and_eq_true, List.all_eq_true] at hg
    obtain ⟨⟨hends, hall⟩, hch⟩ := hg
    have htag : 0 < Nat.add (Nat.mul 2 gi) 1 := Nat.succ_pos _
    obtain ⟨Lg, hg1, hg2, hg3⟩ := group_parts hD hS hQD hR P htag (gends_of hends)
      g.2.2.2.2.2.2.2.2 hall hch 0 0 (fun _ _ => ⟨Nat.zero_le _, Nat.zero_le _⟩)
    obtain ⟨L', h1, h2, h3⟩ := groups_parts hD hS hQD hR P hnd tg (gi + 1) htg
    refine ⟨Lg ++ L', ?_, ?_, ?_⟩
    · intro x hx
      rcases List.mem_append.mp hx with hx | hx
      · obtain ⟨a1, a2, _⟩ := hg1 x hx
        exact ⟨a1, gi, le_rfl, by simpa [nat_add_eq, nat_mul_eq] using a2⟩
      · obtain ⟨a1, j, hj, a2⟩ := h1 x hx
        exact ⟨a1, j, by omega, a2⟩
    · rw [List.pairwise_append]
      refine ⟨hg2, h2, fun x hx y hy hxy => ?_⟩
      obtain ⟨_, a2, _⟩ := hg1 x hx
      obtain ⟨_, j, hj, b2⟩ := h1 y hy
      rw [hxy] at a2
      have := tag_unique hnd a2 b2
      simp only [nat_add_eq, nat_mul_eq] at this
      omega
    · simp only [List.map_cons, List.sum_cons, List.map_append, List.sum_append]
      rw [Nat.cast_add]
      have : ((tgrpVal g : ℕ) : ℝ) ≤ (Lg.map (pval D)).sum := hg3
      linarith

/-- The S-block parts. -/
lemma sparts (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {c : ℝ × ℝ} {u : ℝ}
    (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {sb : List SBlk}
    (hsb : sb.all (sblkOk Q R x0 x1 y0 y1 U0 U1) = true) :
    ∀ (cls : List (SegE × ℕ)), cls.all (sclOk S sb) = true →
      (∀ x ∈ (cls.filter fun ec => Nat.mod ec.2 2 = 0).map (sPart S Q sb),
        PartOK D (sq c (2 * Real.arctan u) 1) x ∧ ∃ t, Nat.mod t 2 = 0 ∧ (x.1, t) ∈ cls) ∧
      (((cls.map (sclVal S sb)).sum : ℕ) : ℝ)
        ≤ (((cls.filter fun ec => Nat.mod ec.2 2 = 0).map (sPart S Q sb)).map (pval D)).sum
  | [], _ => by simp
  | ec :: cls, h => by
    simp only [List.all_cons, Bool.and_eq_true] at h
    obtain ⟨hec, hcls⟩ := h
    obtain ⟨ih1, ih2⟩ := sparts hD hS hQD hR P hsb cls hcls
    by_cases hev : Nat.mod ec.2 2 = 0
    · obtain ⟨hok, hval⟩ := sPart_ok hD hS hQD hR P hsb hev hec
      simp only [List.filter_cons, hev, decide_true, if_true, List.map_cons, List.sum_cons]
      refine ⟨fun x hx => ?_, ?_⟩
      · rcases List.mem_cons.mp hx with rfl | hx
        · exact ⟨hok, ec.2, hev, by rw [sPart_fst]; exact List.mem_cons_self⟩
        · obtain ⟨a, t, ht, hm⟩ := ih1 x hx
          exact ⟨a, t, ht, List.mem_cons_of_mem _ hm⟩
      · rw [Nat.cast_add]; linarith
    · have hz : sclVal S sb ec = 0 := by
        have : Nat.beq (Nat.mod ec.2 2) 0 = false := beq_ne hev
        simp only [sclVal, this, cond_false]
      simp only [List.filter_cons, hev, decide_false, Bool.false_eq_true, if_false, List.map_cons, List.sum_cons, hz]
      refine ⟨fun x hx => ?_, by rw [Nat.cast_add, Nat.cast_zero, zero_add]; linarith⟩
      obtain ⟨a, t, ht, hm⟩ := ih1 x (by simpa using hx)
      exact ⟨a, t, ht, List.mem_cons_of_mem _ hm⟩

/-- **Soundness of the piece certificate**: at every admissible pose of the box, the segment mass
of the square is at least the certificate's value. -/
theorem pc_sound (hD : 0 < D) (hS : 0 < S) (hQD : Q = D * S) (hR : 0 < R) {segs : List SegE}
    {c : ℝ × ℝ} {u : ℝ} (P : Pose S Q R x0 x1 y0 y1 U0 U1 c u) {cls : List (SegE × ℕ)}
    {sb : List SBlk} {tg : List TGrp} (hcls : (cls.map Prod.fst).Sublist segs) (hsnd : segs.Nodup)
    (h : pcOk S Q R x0 x1 y0 y1 U0 U1 cls sb tg = true) :
    (pcVal S cls sb tg : ℝ) ≤ segMass D segs (sq c (2 * Real.arctan u) 1) := by
  simp only [pcOk, Bool.and_eq_true] at h
  obtain ⟨⟨hsb, hsc⟩, htg⟩ := h
  have hnd : (cls.map Prod.fst).Nodup := hcls.nodup hsnd
  obtain ⟨hS1, hS2⟩ := sparts hD hS hQD hR P hsb cls hsc
  obtain ⟨LT, hT1, hT2, hT3⟩ := groups_parts hD hS hQD hR P hnd tg 0 htg
  set LS := (cls.filter fun ec => Nat.mod ec.2 2 = 0).map (sPart S Q sb) with hLS
  have hmemcls : ∀ {e : SegE} {t : ℕ}, (e, t) ∈ cls → e ∈ segs.toFinset := by
    intro e t h
    rw [List.mem_toFinset]
    exact hcls.subset (List.mem_map_of_mem (f := Prod.fst) h)
  have hSnd : (LS.map Prod.fst).Nodup := by
    rw [hLS, List.map_map]
    have : (Prod.fst ∘ sPart S Q sb) = (Prod.fst : SegE × ℕ → SegE) := by
      funext ec; exact sPart_fst S Q sb ec
    rw [this]
    exact hnd.sublist (List.filter_sublist.map _)
  have key := parts_le_segMass hD segs (sq c (2 * Real.arctan u) 1) (LS ++ LT) ?_ ?_ ?_
  · simp only [pcVal, List.map_append, List.sum_append] at key ⊢
    rw [Nat.cast_add]
    linarith
  · intro x hx
    rcases List.mem_append.mp hx with hx | hx
    · obtain ⟨_, t, _, hm⟩ := hS1 x hx; exact hmemcls hm
    · obtain ⟨_, j, _, hm⟩ := hT1 x hx; exact hmemcls hm
  · intro x hx
    rcases List.mem_append.mp hx with hx | hx
    · exact (hS1 x hx).1
    · exact (hT1 x hx).1
  · rw [List.pairwise_append]
    refine ⟨?_, hT2, fun x hx y hy hxy => ?_⟩
    · have := List.nodup_iff_pairwise_ne.mp hSnd
      rw [List.pairwise_map] at this
      exact this.imp fun h e => absurd e h
    · obtain ⟨_, t, ht, hm⟩ := hS1 x hx
      obtain ⟨_, j, _, hm'⟩ := hT1 y hy
      rw [hxy] at hm
      have := tag_unique hnd hm hm'
      have ht' : t % 2 = 0 := ht
      omega

end soundPC2


/-! ## 5.  The mixed tree -/

/-- Segment candidate pruning: keep the segments whose bounding box (over `Q`) meets the box
widened by `F`.  A pure speed device (any sublist is sound). -/
def nearS (S F x0 x1 y0 y1 : ℕ) (c : List SegE) : List SegE :=
  c.filter fun e => Nat.ble x0 (Nat.add (Nat.mul e.2.2.1 S) F) && Nat.ble (Nat.mul e.1 S) (Nat.add x1 F) &&
    Nat.ble y0 (Nat.add (Nat.mul e.2.2.2.1 S) F) && Nat.ble (Nat.mul e.2.1 S) (Nat.add y1 F)

/-- A zero-margin pose-space tree for mixed covers: the leaves of `ZMTree.ZT` with a piece
certificate added to `Z`, and `G` (segment candidate pruning). -/
inductive ZTM
  | Z (cl : List (ℕ × Tag)) (chA chB : List Piv) (emp : List (ℕ × ℕ)) (sc : List (ℕ × ℕ))
      (sb : List SBlk) (tg : List TGrp)
  | E
  | C (us : ℕ) (t : ZTM)
  | X (l r : ZTM)
  | Y (l r : ZTM)
  | U (l r : ZTM)
  | F (t : ZTM)
  | G (t : ZTM)
  | XM (m : ℕ) (l r : ZTM)
  | YM (m : ℕ) (l r : ZTM)
  | UM (m : ℕ) (l r : ZTM)

/-- **The mixed tree check** on the box `[x0,x1]×[y0,y1]×[u0,u1]` with point candidates `c` and
segment candidates `s`.  A `Z` leaf: the piece certificate holds, and `ZMTree.zOk` holds for the
target `W − Lp` (Lemma P). -/
def checkM (W S Q R F : ℕ) : ZTM → ℕ → ℕ → ℕ → ℕ → ℕ → ℕ → List Ent → List SegE → Bool
  | .Z cl chA chB emp sc sb tg, x0, x1, y0, y1, u0, u1, c, s =>
    Nat.ble u1 R && Nat.ble u0 u1 && pcOk S Q R x0 x1 y0 y1 u0 u1 (sclaim sc s) sb tg &&
      zOk (Nat.sub W (pcVal S (sclaim sc s) sb tg)) S Q R x0 x1 y0 y1 u0 u1 (claim cl c) chA chB emp
  | .E, _, x1, _, y1, u0, u1, _, _ => Nat.ble u1 R && eOk Q R x1 y1 u0 u1
  | .C us t, x0, x1, y0, y1, u0, u1, c, s =>
    clipOk Q R x1 y1 u0 u1 us && checkM W S Q R F t x0 x1 y0 y1 u0 us c s
  | .F t, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F t x0 x1 y0 y1 u0 u1 (near S F x0 x1 y0 y1 c) s
  | .G t, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F t x0 x1 y0 y1 u0 u1 c (nearS S F x0 x1 y0 y1 s)
  | .X l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 (Nat.div (Nat.add x0 x1) 2) y0 y1 u0 u1 c s &&
    checkM W S Q R F r (Nat.div (Nat.add x0 x1) 2) x1 y0 y1 u0 u1 c s
  | .Y l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 x1 y0 (Nat.div (Nat.add y0 y1) 2) u0 u1 c s &&
    checkM W S Q R F r x0 x1 (Nat.div (Nat.add y0 y1) 2) y1 u0 u1 c s
  | .U l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 x1 y0 y1 u0 (Nat.div (Nat.add u0 u1) 2) c s &&
    checkM W S Q R F r x0 x1 y0 y1 (Nat.div (Nat.add u0 u1) 2) u1 c s
  | .XM m l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 m y0 y1 u0 u1 c s && checkM W S Q R F r m x1 y0 y1 u0 u1 c s
  | .YM m l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 x1 y0 m u0 u1 c s && checkM W S Q R F r x0 x1 m y1 u0 u1 c s
  | .UM m l r, x0, x1, y0, y1, u0, u1, c, s =>
    checkM W S Q R F l x0 x1 y0 y1 u0 m c s && checkM W S Q R F r x0 x1 y0 y1 m u1 c s

section soundM

variable {D S Mq R W : ℕ} {pts : List Ent} {segs : List SegE}

/-- `clip_bin` for mixed covers (the proof of `ZMTree.C_cov`, for `CovM`). -/
theorem CM_cov (hD : 0 < D) (hS : 0 < S) (hR : 0 < R) {x0 x1 y0 y1 U0 U1 us : ℕ}
    (h : clipOk (D * S) R x1 y1 U0 U1 us = true)
    (hc : CovM D S Mq R W pts segs x0 x1 y0 y1 U0 us) : CovM D S Mq R W pts segs x0 x1 y0 y1 U0 U1 := by
  intro c u hx0 hx1 hy0 hy1 hu0 hu1 hsub
  simp only [clipOk, Bool.and_eq_true, Nat.ble_eq, Nat.blt_eq] at h
  obtain ⟨⟨⟨_, hsU1⟩, hw⟩, hmono⟩ := h
  rcases le_or_gt u ((us : ℝ) / R) with hle | hgt
  · exact hc c u hx0 hx1 hy0 hy1 hu0 hle hsub
  exfalso
  have hQ : 0 < D * S := Nat.mul_pos hD hS
  have hRr : (0 : ℝ) < R := by exact_mod_cast hR
  have hQc : ((D : ℝ) * S) = ((D * S : ℕ) : ℝ) := by push_cast; ring
  rw [hQc] at hx1 hy1
  simp only [nat_mul_eq, nat_add_eq] at hmono
  have hmr : ((R : ℝ) + us) * (R + U1) < 2 * (R * R) := by exact_mod_cast hmono
  have husR : us < R := by
    by_contra hh
    push Not at hh
    have : (R : ℝ) ≤ us := by exact_mod_cast hh
    nlinarith [(Nat.cast_nonneg U1 : (0 : ℝ) ≤ U1)]
  have hus0 : (0 : ℝ) ≤ (us : ℝ) / R := div_nonneg (Nat.cast_nonneg _) hRr.le
  have hu0' : 0 ≤ u := le_trans hus0 hgt.le
  have hU1R : (U1 : ℝ) < R := by nlinarith [(Nat.cast_nonneg us : (0 : ℝ) ≤ us)]
  have hu1' : u ≤ 1 := le_trans hu1 ((div_le_one hRr).mpr hU1R.le)
  obtain ⟨wx, wy⟩ := adm_lo hsub
  rw [wid_two_arctan hu0' hu1'] at wx wy
  have hK := wge_sound hQ hR husR.le hw
  have hprod : 1 - ((us : ℝ) / R + u + (us : ℝ) / R * u) > 0 := by
    have hu1R : u * R ≤ U1 := by rw [le_div_iff₀ hRr] at hu1; linarith
    have e : 1 - ((us : ℝ) / R + u + (us : ℝ) / R * u)
        = (2 * (R * R) - (R + us) * (R + u * R)) / (R * R) := by field_simp; ring
    rw [e]
    apply div_pos _ (by positivity)
    nlinarith [(Nat.cast_nonneg us : (0 : ℝ) ≤ us)]
  have hinc : widU ((us : ℝ) / R) < widU u := by
    have := widU_sub ((us : ℝ) / R) u
    have hd : (0 : ℝ) < (1 + ((us : ℝ) / R) ^ 2) * (1 + u ^ 2) := by positivity
    have : 0 < widU u - widU ((us : ℝ) / R) := by
      rw [this]; exact div_pos (mul_pos (by linarith) hprod) hd
    linarith
  rcases min_choice x1 y1 with hm | hm <;> rw [hm] at hK
  · linarith
  · linarith

/-- **Soundness of a mixed `Z` leaf** (Lemma P + the piece certificate + `ZMTree.Z_cov`). -/
theorem ZM_cov (hD : 0 < D) (hS : 0 < S) (hR : 0 < R) (hnd : pts.Nodup) (hsnd : segs.Nodup)
    {x0 x1 y0 y1 U0 U1 : ℕ} {cl : List (Ent × Tag)} {chA chB : List Piv} {emp : List (ℕ × ℕ)}
    {cls : List (SegE × ℕ)} {sb : List SBlk} {tg : List TGrp}
    (hcl : (cl.map Prod.fst).Sublist pts) (hcls : (cls.map Prod.fst).Sublist segs)
    (hU01 : U0 ≤ U1) (hU1 : U1 ≤ R)
    (hpc : pcOk S (D * S) R x0 x1 y0 y1 U0 U1 cls sb tg = true)
    (hz : zOk (Nat.sub W (pcVal S cls sb tg)) S (D * S) R x0 x1 y0 y1 U0 U1 cl chA chB emp = true) :
    CovM D S Mq R W pts segs x0 x1 y0 y1 U0 U1 := by
  have hcov := Z_cov (Mq := Mq) hD hS hR hnd hcl hU01 hU1 hz
  intro cc u hx0 hx1 hy0 hy1 hu0 hu1 hsub
  have hp := hcov cc u hx0 hx1 hy0 hy1 hu0 hu1 hsub
  have hQc : ((D : ℝ) * S) = ((D * S : ℕ) : ℝ) := by push_cast; ring
  rw [hQc] at hx0 hx1 hy0 hy1
  obtain ⟨wx, wy⟩ := adm_lo hsub
  have P : Pose S (D * S) R x0 x1 y0 y1 U0 U1 cc u :=
    ⟨hx0, hx1, hy0, hy1, hu0, hu1, wx, wy, hU01, hU1⟩
  have hs := pc_sound hD hS rfl hR P hcls hsnd hpc
  set Lp := pcVal S cls sb tg
  have hW : (W : ℝ) ≤ ((Nat.sub W Lp : ℕ) : ℝ) + Lp := by
    rcases le_total Lp W with h | h
    · rw [nat_sub_eq, Nat.cast_sub h]; linarith
    · rw [nat_sub_eq, Nat.sub_eq_zero_of_le h]; push_cast; rw [zero_add]; exact_mod_cast h
  unfold ptMass
  linarith

/-- **Soundness of the mixed tree check.** -/
theorem soundM (D S Mq R W F : ℕ) (pts : List Ent) (segs : List SegE) (hD : 0 < D) (hS : 0 < S)
    (hR : 0 < R) (hnd : pts.Nodup) (hsnd : segs.Nodup) :
    ∀ (t : ZTM) (x0 x1 y0 y1 u0 u1 : ℕ) (c : List Ent) (s : List SegE), c.Sublist pts →
      s.Sublist segs → checkM W S (D * S) R F t x0 x1 y0 y1 u0 u1 c s = true →
      CovM D S Mq R W pts segs x0 x1 y0 y1 u0 u1 := by
  intro t
  induction t with
  | Z cl chA chB emp sc sb tg =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true, Nat.ble_eq] at h
    obtain ⟨⟨⟨hu1, hu01⟩, hpc⟩, hz⟩ := h
    exact ZM_cov hD hS hR hnd hsnd ((claim_sub cl c).trans hc) ((sclaim_sub sc s).trans hs) hu01 hu1
      hpc hz
  | E =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true, Nat.ble_eq] at h
    exact CovM.of_cov (E_cov hD hS hR h.1 h.2)
  | C us t ih =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CM_cov hD hS hR h.1 (ih _ _ _ _ _ _ _ _ hc hs h.2)
  | F t ih =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM] at h
    exact ih _ _ _ _ _ _ _ _ ((List.filter_sublist).trans hc) hs h
  | G t ih =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM] at h
    exact ih _ _ _ _ _ _ _ _ hc ((List.filter_sublist).trans hs) h
  | X l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitX _ (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)
  | Y l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitY _ (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)
  | U l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitU _ (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)
  | XM m l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitX m (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)
  | YM m l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitY m (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)
  | UM m l r ihl ihr =>
    intro x0 x1 y0 y1 u0 u1 c s hc hs h
    simp only [checkM, Bool.and_eq_true] at h
    exact CovM.splitU m (ihl _ _ _ _ _ _ _ _ hc hs h.1) (ihr _ _ _ _ _ _ _ _ hc hs h.2)

end soundM

/-! ## 6.  Compact encoding

As `ZMTree.dec` (structure digits base `B`: `0` a `Z` leaf, `1`–`3` `X`/`Y`/`U`, `4` `F`, `5` `E`,
`6` `C`, `7`–`9` `XM`/`YM`/`UM`, and `10` `G`), but a `Z` leaf takes **two** numbers from the leaf
list: its point part (`ZMTree.decZ` format) and its piece part: `nsc`, then `nsc` claims (gap, tag);
`nsb`, then S-blocks `dir, K, A, B`; `ntg`, then T-groups `dir, K, Au, Bu, Ad, Bd, τu+1, τd+1, np`
and `np` pairs `ju+1, a, b, jd+1, c, d, m` — coordinates and masses as three base-`B` digits.  The
decoder needs no proof. -/

def bigD (B n : ℕ) : ℕ × ℕ :=
  (Nat.add (Nat.mod n B) (Nat.mul B (Nat.add (Nat.mod (Nat.div n B) B)
      (Nat.mul B (Nat.mod (Nat.div (Nat.div n B) B) B)))),
    Nat.div (Nat.div (Nat.div n B) B) B)

def decScs (B : ℕ) : ℕ → ℕ → List (ℕ × ℕ) × ℕ
  | 0, n => ([], n)
  | k + 1, n =>
    match decScs B k (Nat.div (Nat.div n B) B) with
    | (l, r) => ((Nat.mod n B, Nat.mod (Nat.div n B) B) :: l, r)

def decSb (B n : ℕ) : SBlk × ℕ :=
  match bigD B (Nat.div n B) with
  | (K, n1) => match bigD B n1 with
    | (A, n2) => match bigD B n2 with
      | (Bv, n3) => ((Nat.mod n B, K, A, Bv), n3)

def decSbs (B : ℕ) : ℕ → ℕ → List SBlk × ℕ
  | 0, n => ([], n)
  | k + 1, n => match decSb B n with
    | (b, n1) => match decSbs B k n1 with
      | (l, r) => (b :: l, r)

def decTp (B n : ℕ) : TPr × ℕ :=
  match bigD B (Nat.div n B) with
  | (a, n1) => match bigD B n1 with
    | (b, n2) => match bigD B (Nat.div n2 B) with
      | (c, n3) => match bigD B n3 with
        | (d, n4) => match bigD B n4 with
          | (m, n5) => ((Nat.mod n B, a, b, Nat.mod n2 B, c, d, m), n5)

def decTps (B : ℕ) : ℕ → ℕ → List TPr × ℕ
  | 0, n => ([], n)
  | k + 1, n => match decTp B n with
    | (p, n1) => match decTps B k n1 with
      | (l, r) => (p :: l, r)

def decTg (B n : ℕ) : TGrp × ℕ :=
  match bigD B (Nat.div n B) with
  | (K, n1) => match bigD B n1 with
    | (Au, n2) => match bigD B n2 with
      | (Bu, n3) => match bigD B n3 with
        | (Ad, n4) => match bigD B n4 with
          | (Bd, n5) => match bigD B n5 with
            | (tu, n6) => match bigD B n6 with
              | (td, n7) => match decTps B (Nat.mod n7 B) (Nat.div n7 B) with
                | (prs, n8) => ((Nat.mod n B, K, Au, Bu, Ad, Bd, tu, td, prs), n8)

def decTgs (B : ℕ) : ℕ → ℕ → List TGrp × ℕ
  | 0, n => ([], n)
  | k + 1, n => match decTg B n with
    | (g, n1) => match decTgs B k n1 with
      | (l, r) => (g :: l, r)

/-- The piece part of a `Z` leaf. -/
def decPc (B n : ℕ) : List (ℕ × ℕ) × List SBlk × List TGrp :=
  match decScs B (Nat.mod n B) (Nat.div n B) with
  | (sc, n1) => match decSbs B (Nat.mod n1 B) (Nat.div n1 B) with
    | (sb, n2) => ((sc, sb, (decTgs B (Nat.mod n2 B) (Nat.div n2 B)).1))

/-- A mixed `Z` leaf from its two numbers. -/
def decZM (B n np : ℕ) : ZTM :=
  match ZMTree.decZ B n with
  | .Z cl chA chB emp => match decPc B np with
    | (sc, sb, tg) => .Z cl chA chB emp sc sb tg
  | _ => .E

/-- Decode a tree of depth `≤ fuel`. -/
def decM (B : ℕ) : ℕ → ℕ → List ℕ → ZTM × ℕ × List ℕ
  | 0, n, L => (.E, n, L)
  | f + 1, n, L =>
    match Nat.mod n B with
    | 0 => (decZM B (L.headD 0) (L.tail.headD 0), Nat.div n B, L.tail.tail)
    | 1 => match decM B f (Nat.div n B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) => (.X l r, n2, L2)
    | 2 => match decM B f (Nat.div n B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) => (.Y l r, n2, L2)
    | 3 => match decM B f (Nat.div n B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) => (.U l r, n2, L2)
    | 4 => match decM B f (Nat.div n B) L with
      | (t, n1, L1) => (.F t, n1, L1)
    | 5 => (.E, Nat.div n B, L)
    | 7 => match decM B f (Nat.div (Nat.div (Nat.div n B) B) B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) =>
          (.XM (Nat.add (Nat.mod (Nat.div n B) B) (Nat.mul B (Nat.mod (Nat.div (Nat.div n B) B) B)))
            l r, n2, L2)
    | 8 => match decM B f (Nat.div (Nat.div (Nat.div n B) B) B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) =>
          (.YM (Nat.add (Nat.mod (Nat.div n B) B) (Nat.mul B (Nat.mod (Nat.div (Nat.div n B) B) B)))
            l r, n2, L2)
    | 9 => match decM B f (Nat.div (Nat.div (Nat.div n B) B) B) L with
      | (l, n1, L1) => match decM B f n1 L1 with
        | (r, n2, L2) =>
          (.UM (Nat.add (Nat.mod (Nat.div n B) B) (Nat.mul B (Nat.mod (Nat.div (Nat.div n B) B) B)))
            l r, n2, L2)
    | 10 => match decM B f (Nat.div n B) L with
      | (t, n1, L1) => (.G t, n1, L1)
    | _ => match decM B f (Nat.div (Nat.div (Nat.div n B) B) B) L with
      | (t, n1, L1) =>
        (.C (Nat.add (Nat.mod (Nat.div n B) B) (Nat.mul B (Nat.mod (Nat.div (Nat.div n B) B) B))) t,
          n1, L1)

end ZMTreeM

end SquarePacking
