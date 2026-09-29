import Sqpack.LBlock

/-!
# Soundness of the L-block check

* §1 options: the value of an option at a pose (`oval`), its corner form, the per-corner inequalities
  (`corner_facts`).
* §2 concavity: the lower bound `Φ` (sum over lines and ends of the minimum over options) is concave in
  the centre, hence `≥ lg` on the box.
* §3 chord and gains: at a pose, each line's chord contains `[b − y, a + x]`; the pieces give parts
  whose masses are at least the core plus `Φ`'s terms.
-/

namespace SquarePacking

namespace ZMTreeM

open BoxTree ZMTree

/-! ## 1.  Options -/

/-- A quadratic's value. -/
noncomputable def e3 (n : SP × SP × SP) (v : ℝ) : ℝ := sv n.1 + sv n.2.1 * v + sv n.2.2 * v ^ 2

/-- The slope magnitude of a condition of type `f` along the line: `4u` (`f = 1`), `2(1 − u²)`. -/
noncomputable def sig (f : ℕ) (u : ℝ) : ℝ := if f = 1 then 4 * u else 2 * (1 - u ^ 2)

/-- The value of an option at a pose (fine units): the cap, or `i − s Q G_k(end point)/σ_k(u)`. -/
noncomputable def oval (Q : ℕ) (l : LLine) (up : Bool) (cap o : ℕ) (c : ℝ × ℝ) (u : ℝ) : ℝ :=
  if o = 0 then (cap : ℝ) else
    sv (if up then l.iU else l.iD) - ((if up then l.sU else l.sD : ℕ) : ℝ) * Q *
      Gv Q (o - 1) (lpx l (if up then l.a else l.b)) (lpy l (if up then l.a else l.b)) c u /
        sig (ktype l.dir (o - 1)) u

lemma ktype_cases (dir k : ℕ) : ktype dir k = 1 ∨ ktype dir k = 2 := by
  unfold ktype
  cases Nat.beq dir 0 <;> cases Nat.blt k 2 <;> simp

lemma sigq_e3 (R f : ℕ) (hf : f = 1 ∨ f = 2) (v : ℝ) : e3 (sigq R f) v = fval f R v := by
  rcases hf with rfl | rfl
  · simp only [sigq, e3, fval, beq_rfl, cond_true, sv, nat_mul_eq]; push_cast; ring
  · simp only [sigq, e3, fval, show Nat.beq 2 1 = false from rfl, cond_false, sv, nat_mul_eq]
    push_cast; ring

lemma fval_sig {R f : ℕ} (hf : f = 1 ∨ f = 2) (u : ℝ) :
    fval f R (u * R) = (R : ℝ) ^ 2 * sig f u := by
  rcases hf with rfl | rfl
  · simp only [fval, sig, if_true]; ring
  · simp only [fval, sig, show (2 : ℕ) ≠ 1 by norm_num, if_false]; ring

/-- The corner form of an option: `N(v)/σ̂(v)` is the option's value at the corner. -/
lemma optN_val {Q R : ℕ} (hQ : 0 < Q) (hR : 0 < R) (l : LLine) (s : ℕ) (i : SP) (t k cx cy : ℕ)
    (hk : k < 4) (u : ℝ) (hsig : sig (ktype l.dir k) u ≠ 0) :
    e3 (optN Q R l s i t k cx cy) (u * R) / fval (ktype l.dir k) R (u * R)
      = sv i - (s : ℝ) * Q * Gv Q k (lpx l t) (lpy l t) ((cx : ℝ) / Q, (cy : ℝ) / Q) u /
          sig (ktype l.dir k) u := by
  have hQr : (Q : ℝ) ≠ 0 := by exact_mod_cast hQ.ne'
  have hRr : (R : ℝ) ≠ 0 := by exact_mod_cast hR.ne'
  have hf := ktype_cases l.dir k
  have hG := gq_eval Q R k (lpx l t, cx) (lpy l t, cy) hQ hR (u * R)
  have huR : u * R / R = u := by field_simp
  rw [huR] at hG
  have hs := sigq_e3 R (ktype l.dir k) hf (u * R)
  have hfv := fval_sig (R := R) hf u
  have hGv : Gv Q k (lpx l t) (lpy l t) ((cx : ℝ) / Q, (cy : ℝ) / Q) u
      = gval (kf k) (sv (lpx l t, cx) / Q) (sv (lpy l t, cy) / Q) u := by
    simp only [Gv, sv_mk, sub_div]
  have e : e3 (optN Q R l s i t k cx cy) (u * R)
      = sv i * e3 (sigq R (ktype l.dir k)) (u * R)
        - s * (sv (gq Q R k (lpx l t, cx) (lpy l t, cy)).1
          + sv (gq Q R k (lpx l t, cx) (lpy l t, cy)).2.1 * (u * R)
          + sv (gq Q R k (lpx l t, cx) (lpy l t, cy)).2.2 * (u * R) ^ 2) := by
    simp only [optN, e3, sv_ssub, sv_smul, sv_sk]; ring
  rw [e, hs, hG, hfv, hGv]
  have hR2 : (R : ℝ) ^ 2 ≠ 0 := pow_ne_zero 2 hRr
  field_simp


/-- The value of an option triple `(type, N)` at `v`. -/
noncomputable def optval (R : ℕ) (t : ℕ × (SP × SP × SP)) (v : ℝ) : ℝ :=
  if t.1 = 0 then sv t.2.1 else e3 t.2 v / fval t.1 R v

lemma optT_type (Q R : ℕ) (l : LLine) (up : Bool) (cap o cx cy : ℕ) :
    (optT Q R l up cap o cx cy).1 = if o = 0 then 0 else ktype l.dir (o - 1) := by
  by_cases h : o = 0
  · subst h; simp [optT]
  · simp [optT, beq_ne h, h]

/-- The corner triple of an option has the option's value. -/
lemma optT_val {Q R : ℕ} (hQ : 0 < Q) (hR : 0 < R) (l : LLine) (up : Bool) (cap o cx cy : ℕ)
    (ho : o ≤ 4) (u : ℝ) (hsig : o ≠ 0 → sig (ktype l.dir (o - 1)) u ≠ 0) :
    optval R (optT Q R l up cap o cx cy) (u * R) = oval Q l up cap o ((cx : ℝ) / Q, (cy : ℝ) / Q) u := by
  by_cases h : o = 0
  · subst h; simp [optval, optT, oval, sv]
  · have ht : (optT Q R l up cap o cx cy).1 ≠ 0 := by
      rw [optT_type, if_neg h]
      rcases ktype_cases l.dir (o - 1) with h1 | h1 <;> omega
    simp only [optval, if_neg ht, oval, if_neg h]
    rw [optT_type, if_neg h]
    have hk : o - 1 < 4 := by omega
    cases up
    · simp only [optT, beq_ne h, cond_false, Bool.false_eq_true, if_false]
      exact optN_val hQ hR l l.sD l.iD l.b (o - 1) cx cy hk u (hsig h)
    · simp only [optT, beq_ne h, cond_false, cond_true, if_true]
      exact optN_val hQ hR l l.sU l.iU l.a (o - 1) cx cy hk u (hsig h)

/-- `a ⊆ f` for the types. -/
def tsub (a f : ℕ) : Prop := a = 0 ∨ a = f ∨ f = 3

lemma tor_le (a b : ℕ) (ha : a ≤ 3) (hb : b ≤ 3) : tor a b ≤ 3 := by
  unfold tor
  rcases Nat.eq_zero_or_pos a with rfl | h1
  · simpa using hb
  rcases Nat.eq_zero_or_pos b with rfl | h2
  · simp [beq_ne (Nat.pos_iff_ne_zero.mp h1)]; omega
  simp only [beq_ne (Nat.pos_iff_ne_zero.mp h1), beq_ne (Nat.pos_iff_ne_zero.mp h2), cond_false]
  by_cases h : a = b
  · subst h; simp [beq_rfl]; omega
  · simp [beq_ne h]

lemma tsub_tor_l (a b : ℕ) (ha : a ≤ 3) (hb : b ≤ 3) : tsub a (tor a b) := by
  unfold tsub tor
  rcases Nat.eq_zero_or_pos a with rfl | h1
  · exact Or.inl rfl
  rcases Nat.eq_zero_or_pos b with rfl | h2
  · simp [beq_ne (Nat.pos_iff_ne_zero.mp h1)]
  simp only [beq_ne (Nat.pos_iff_ne_zero.mp h1), beq_ne (Nat.pos_iff_ne_zero.mp h2), cond_false]
  by_cases h : a = b
  · subst h; simp [beq_rfl]
  · simp [beq_ne h]

lemma tsub_tor_r (a b : ℕ) (ha : a ≤ 3) (hb : b ≤ 3) : tsub b (tor a b) := by
  unfold tsub tor
  rcases Nat.eq_zero_or_pos b with rfl | h2
  · exact Or.inl rfl
  rcases Nat.eq_zero_or_pos a with rfl | h1
  · simp
  simp only [beq_ne (Nat.pos_iff_ne_zero.mp h1), beq_ne (Nat.pos_iff_ne_zero.mp h2), cond_false]
  by_cases h : a = b
  · subst h; simp [beq_rfl]
  · simp [beq_ne h]

lemma tsub_trans {a f g : ℕ} (h1 : tsub a f) (h2 : tsub f g) : tsub a g := by
  unfold tsub at *; omega

/-- The type conditions justifying a denominator. -/
def goodT (R U0 U1 f : ℕ) : Prop := f ≤ 3 ∧ (f % 2 = 1 → 0 < U0) ∧ (2 ≤ f → U1 < R)

lemma goodT_tor {R U0 U1 a b : ℕ} (ha : goodT R U0 U1 a) (hb : goodT R U0 U1 b) :
    goodT R U0 U1 (tor a b) := by
  obtain ⟨a3, a1, a2⟩ := ha
  obtain ⟨b3, b1, b2⟩ := hb
  refine ⟨tor_le a b a3 b3, ?_, ?_⟩ <;> unfold tor
  · rcases Nat.eq_zero_or_pos a with rfl | h1
    · simpa using b1
    rcases Nat.eq_zero_or_pos b with rfl | h2
    · simpa [beq_ne (Nat.pos_iff_ne_zero.mp h1)] using a1
    simp only [beq_ne (Nat.pos_iff_ne_zero.mp h1), beq_ne (Nat.pos_iff_ne_zero.mp h2), cond_false]
    by_cases h : a = b
    · subst h; simpa [beq_rfl] using a1
    · simp only [beq_ne h, cond_false]
      intro _
      rcases Nat.lt_or_ge a 2 with h' | h'
      · exact a1 (by omega)
      · rcases Nat.lt_or_ge b 2 with h'' | h''
        · exact b1 (by omega)
        · omega
  · rcases Nat.eq_zero_or_pos a with rfl | h1
    · simpa using b2
    rcases Nat.eq_zero_or_pos b with rfl | h2
    · simpa [beq_ne (Nat.pos_iff_ne_zero.mp h1)] using a2
    simp only [beq_ne (Nat.pos_iff_ne_zero.mp h1), beq_ne (Nat.pos_iff_ne_zero.mp h2), cond_false]
    by_cases h : a = b
    · subst h; simpa [beq_rfl] using a2
    · simp only [beq_ne h, cond_false]
      intro _
      rcases Nat.lt_or_ge a 2 with h' | h'
      · exact b2 (by omega)
      · exact a2 h'

lemma fval_pos_good {R U0 U1 f : ℕ} (hR : 0 < R) (hg : goodT R U0 U1 f) {v : ℝ}
    (hv0 : (U0 : ℝ) ≤ v) (hv1 : v ≤ U1) : 0 < fval f R v := by
  obtain ⟨h3, h1, h2⟩ := hg
  refine fval_pos hR h3 (fun h => lt_of_lt_of_le (by exact_mod_cast h1 h) hv0)
    (fun h => lt_of_le_of_lt hv1 (by exact_mod_cast h2 h)) (le_trans (Nat.cast_nonneg _) hv0)

lemma fval_split {R f a : ℕ} (hsub : tsub a f) (ha : a = 1 ∨ a = 2) (v : ℝ) :
    fval f R v = fval (f - a) R v * fval a R v := by
  unfold tsub at hsub
  rcases hsub with h | h | h
  · omega
  · subst h; simp [fval]
  · subst h; rcases ha with rfl | rfl <;> simp only [fval, show (3 : ℕ) - 1 = 2 from rfl,
      show (3 : ℕ) - 2 = 1 from rfl] <;> ring

/-- The term of an option in a sum over the common denominator `f`. -/
lemma e5_term {R f : ℕ} (t : ℕ × (SP × SP × SP)) (hf : f ≤ 3) (ht : t.1 = 0 ∨ t.1 = 1 ∨ t.1 = 2)
    (hsub : tsub t.1 f) (v : ℝ) (hne : t.1 ≠ 0 → fval t.1 R v ≠ 0) :
    e5 (term R f t) v = fval f R v * optval R t v := by
  unfold term optval
  simp only [nat_sub_eq]
  rcases ht with h0 | h12
  · rw [h0]; simp only [beq_rfl, cond_true, if_true]; exact e5_cmul f R hf _ v
  · have hne0 : t.1 ≠ 0 := by omega
    have hb : Nat.beq t.1 0 = false := beq_ne hne0
    rw [hb, cond_false, if_neg hne0, e5_tmul _ R (by unfold tsub at hsub; omega),
      fval_split hsub h12 v]
    have hfv := hne hne0
    field_simp
    simp only [e3]
    ring


/-! ### The corner inequalities -/

lemma isOpt_spec {l : LLine} {up : Bool} {o : ℕ} (h : isOpt l up o = true) :
    o ≤ 4 ∧ (o ≠ 0 → role l.roles (o - 1) = (if up then 1 else 2)) := by
  unfold isOpt at h
  simp only [Bool.or_eq_true, Bool.and_eq_true, Nat.beq_eq, Nat.ble_eq] at h
  rcases h with h | ⟨h1, h2⟩
  · subst h; exact ⟨by norm_num, fun h => absurd rfl h⟩
  · refine ⟨h1, fun _ => ?_⟩
    simp only [nat_sub_eq] at h2
    rw [h2]; cases up <;> rfl

lemma rolesOk_k {R U0 U1 dir roles k : ℕ} (h : rolesOk R U0 U1 dir roles = true) (hk : k < 4) :
    (role roles k = 0 ∨ role roles k = kend dir k) ∧
      (role roles k = 0 ∨ (ktype dir k = 1 → 0 < U0) ∧ (ktype dir k ≠ 1 → U1 < R)) := by
  unfold rolesOk at h
  have := List.all_eq_true.mp h k (List.mem_range.mpr hk)
  simp only [Bool.and_eq_true, Bool.or_eq_true, Nat.beq_eq] at this
  obtain ⟨h1, h2⟩ := this
  refine ⟨h1, ?_⟩
  rcases h2 with h2 | h2
  · exact Or.inl h2
  · right
    by_cases ht : ktype dir k = 1
    · have : Nat.beq (ktype dir k) 1 = true := by rw [ht]; rfl
      rw [this, cond_true, Nat.blt_eq] at h2
      exact ⟨fun _ => h2, fun h => absurd ht h⟩
    · have : Nat.beq (ktype dir k) 1 = false := beq_ne ht
      rw [this, cond_false, Nat.blt_eq] at h2
      exact ⟨fun h => absurd h ht, fun _ => h2⟩

/-- The type of an option is justified by the roles. -/
lemma opt_good {Q R U0 U1 : ℕ} {l : LLine} {up : Bool} {cap o cx cy : ℕ}
    (hro : rolesOk R U0 U1 l.dir l.roles = true) (ho : isOpt l up o = true) :
    goodT R U0 U1 (optT Q R l up cap o cx cy).1 ∧
      ((optT Q R l up cap o cx cy).1 = 0 ∨ (optT Q R l up cap o cx cy).1 = 1 ∨
        (optT Q R l up cap o cx cy).1 = 2) := by
  obtain ⟨ho4, hrole⟩ := isOpt_spec ho
  rw [optT_type]
  by_cases h0 : o = 0
  · simp only [h0, if_true]
    exact ⟨⟨by norm_num, fun h => absurd h (by norm_num), fun h => absurd h (by norm_num)⟩, by simp⟩
  · simp only [h0, if_false]
    have hr := hrole h0
    obtain ⟨_, h2⟩ := rolesOk_k hro (show o - 1 < 4 by omega)
    have hne : role l.roles (o - 1) ≠ 0 := by rw [hr]; cases up <;> simp
    rcases h2 with h2 | ⟨hS, hC⟩
    · exact absurd h2 hne
    rcases ktype_cases l.dir (o - 1) with ht | ht
    · rw [ht]; exact ⟨⟨by norm_num, fun _ => hS ht, fun h => absurd h (by norm_num)⟩, Or.inr (Or.inl rfl)⟩
    · rw [ht]; exact ⟨⟨by norm_num, fun h => absurd h (by norm_num), fun _ => hC (by omega)⟩, Or.inr (Or.inr rfl)⟩

lemma sig_pos {R U0 U1 f : ℕ} (hR : 0 < R) (hg : goodT R U0 U1 f) (hf : f = 1 ∨ f = 2) {u : ℝ}
    (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) : 0 < sig f u := by
  have hRr : (0 : ℝ) < R := by exact_mod_cast hR
  obtain ⟨_, h1, h2⟩ := hg
  rcases hf with rfl | rfl
  · simp only [sig, if_true]
    have := h1 (by norm_num)
    have : (0 : ℝ) < (U0 : ℝ) / R := div_pos (by exact_mod_cast this) hRr
    linarith
  · simp only [sig, show (2 : ℕ) ≠ 1 by norm_num, if_false]
    have := h2 le_rfl
    have hlt : (U1 : ℝ) / R < 1 := (div_lt_one hRr).mpr (by exact_mod_cast this)
    have hu : 0 ≤ u := le_trans (div_nonneg (Nat.cast_nonneg _) hRr.le) hu0
    nlinarith

/-- The value of an option at a corner, from its triple, on the bin. -/
lemma optval_oval {Q R U0 U1 : ℕ} (hQ : 0 < Q) (hR : 0 < R) {l : LLine} {up : Bool}
    {cap o cx cy : ℕ} (hro : rolesOk R U0 U1 l.dir l.roles = true) (ho : isOpt l up o = true)
    {u : ℝ} (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) :
    optval R (optT Q R l up cap o cx cy) (u * R) = oval Q l up cap o ((cx : ℝ) / Q, (cy : ℝ) / Q) u := by
  obtain ⟨ho4, _⟩ := isOpt_spec ho
  refine optT_val hQ hR l up cap o cx cy ho4 u fun h0 => ?_
  obtain ⟨hg, _⟩ := opt_good (Q := Q) (cap := cap) (cx := cx) (cy := cy) hro ho
  rw [optT_type, if_neg h0] at hg
  exact (sig_pos hR hg (ktype_cases _ _) hu0 hu1).ne'

lemma v_range {R U0 U1 : ℕ} (hR : 0 < R) {u : ℝ} (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) :
    (U0 : ℝ) ≤ u * R ∧ u * R ≤ U1 := by
  have hRr : (0 : ℝ) < R := by exact_mod_cast hR
  exact ⟨by rw [div_le_iff₀ hRr] at hu0; linarith, by rw [le_div_iff₀ hRr] at hu1; linarith⟩

lemma fval_ne_of_good {R U0 U1 f : ℕ} (hR : 0 < R) (hg : goodT R U0 U1 f) {u : ℝ}
    (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) : fval f R (u * R) ≠ 0 := by
  obtain ⟨v0, v1⟩ := v_range hR hu0 hu1
  exact (fval_pos_good hR hg v0 v1).ne'

/-- **A slack check**, read back: at the corner, `chosen − other ≤ σ` on the bin. -/
lemma slack_fact {Q R U0 U1 : ℕ} (hQ : 0 < Q) (hR : 0 < R) (hU01 : U0 ≤ U1) {l : LLine}
    {up : Bool} {cap ch σ cx cy o : ℕ} (hro : rolesOk R U0 U1 l.dir l.roles = true)
    (hch : isOpt l up ch = true) (ho : isOpt l up o = true)
    (h : slackOk Q R U0 U1 l up cap ch σ cx cy o = true) {u : ℝ} (hu0 : (U0 : ℝ) / R ≤ u)
    (hu1 : u ≤ (U1 : ℝ) / R) :
    oval Q l up cap ch ((cx : ℝ) / Q, (cy : ℝ) / Q) u - oval Q l up cap o ((cx : ℝ) / Q, (cy : ℝ) / Q) u
      ≤ σ := by
  obtain ⟨v0, v1⟩ := v_range hR hu0 hu1
  unfold slackOk at h
  have hb := bOk5_e5 h hU01 v0 v1
  obtain ⟨g1, t1⟩ := opt_good (Q := Q) (cap := cap) (cx := cx) (cy := cy) hro hch
  obtain ⟨g2, t2⟩ := opt_good (Q := Q) (cap := cap) (cx := cx) (cy := cy) hro ho
  set T1 := optT Q R l up cap ch cx cy
  set T2 := optT Q R l up cap o cx cy
  set f := tor T1.1 T2.1
  have hg : goodT R U0 U1 f := goodT_tor g1 g2
  have hf3 : f ≤ 3 := hg.1
  have e1 := e5_term (R := R) T1 hf3 t1 (tsub_tor_l _ _ g1.1 g2.1) (u * R)
    (fun _ => fval_ne_of_good hR g1 hu0 hu1)
  have e2 := e5_term (R := R) T2 hf3 t2 (tsub_tor_r _ _ g1.1 g2.1) (u * R)
    (fun _ => fval_ne_of_good hR g2 hu0 hu1)
  simp only [e5_add, e5_neg] at hb
  rw [e1, e2, e5_cmul f R hf3] at hb
  rw [optval_oval hQ hR hro hch hu0 hu1, optval_oval hQ hR hro ho hu0 hu1] at hb
  have hpos := fval_pos_good hR hg v0 v1
  simp only [sv_mk, Nat.cast_zero, sub_zero] at hb
  by_contra hc
  push Not at hc
  have := mul_pos hpos (sub_pos.mpr hc)
  linarith

/-- The sum of the chosen terms. -/
lemma mainP_e5 {Q R S U0 U1 f : ℕ} (hQ : 0 < Q) (hR : 0 < R) (cls : List (SegE × ℕ)) (cx cy : ℕ)
    {u : ℝ} (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) (hf3 : f ≤ 3) :
    ∀ lc : List (LLine × (ℕ × ℕ × ℕ × ℕ)),
      (∀ x ∈ lc, rolesOk R U0 U1 x.1.dir x.1.roles = true ∧ isOpt x.1 true x.2.1 = true ∧
        isOpt x.1 false x.2.2.2.1 = true ∧
        tsub (optT Q R x.1 true (capU S Q cls x.1) x.2.1 cx cy).1 f ∧
        tsub (optT Q R x.1 false (capD S Q cls x.1) x.2.2.2.1 cx cy).1 f) →
      e5 (mainP Q R S f cls cx cy lc).1 (u * R)
        = fval f R (u * R) * (lc.map fun x =>
            oval Q x.1 true (capU S Q cls x.1) x.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u +
            oval Q x.1 false (capD S Q cls x.1) x.2.2.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u).sum ∧
      (mainP Q R S f cls cx cy lc).2 = (lc.map fun x => x.2.2.1 + x.2.2.2.2).sum
  | [], _ => by simp [mainP, e5, z0, sv]
  | x :: t, h => by
    obtain ⟨hro, hU, hD, sU, sD⟩ := h x List.mem_cons_self
    obtain ⟨ih1, ih2⟩ := mainP_e5 hQ hR cls cx cy hu0 hu1 hf3 t
      (fun y hy => h y (List.mem_cons_of_mem _ hy))
    obtain ⟨gU, tU⟩ := opt_good (Q := Q) (cap := capU S Q cls x.1) (cx := cx) (cy := cy) hro hU
    obtain ⟨gD, tD⟩ := opt_good (Q := Q) (cap := capD S Q cls x.1) (cx := cx) (cy := cy) hro hD
    obtain ⟨l, c⟩ := x
    simp only [mainP, e5_add, List.map_cons, List.sum_cons, nat_add_eq] at ih1 ih2 ⊢
    refine ⟨?_, by rw [ih2]⟩
    rw [ih1, e5_term _ hf3 tU sU (u * R) (fun _ => fval_ne_of_good hR gU hu0 hu1),
      e5_term _ hf3 tD sD (u * R) (fun _ => fval_ne_of_good hR gD hu0 hu1),
      optval_oval hQ hR hro hU hu0 hu1, optval_oval hQ hR hro hD hu0 hu1]
    ring

lemma tsub_mainF {Q R S : ℕ} (cls : List (SegE × ℕ)) (cx cy : ℕ) {U0 U1 : ℕ} :
    ∀ lc : List (LLine × (ℕ × ℕ × ℕ × ℕ)),
      (∀ x ∈ lc, rolesOk R U0 U1 x.1.dir x.1.roles = true ∧ isOpt x.1 true x.2.1 = true ∧
        isOpt x.1 false x.2.2.2.1 = true) →
      goodT R U0 U1 (mainF Q R S cls cx cy lc) ∧
      ∀ x ∈ lc, tsub (optT Q R x.1 true (capU S Q cls x.1) x.2.1 cx cy).1 (mainF Q R S cls cx cy lc) ∧
        tsub (optT Q R x.1 false (capD S Q cls x.1) x.2.2.2.1 cx cy).1 (mainF Q R S cls cx cy lc)
  | [], _ => ⟨⟨by simp [mainF], by simp [mainF], by simp [mainF]⟩, by simp⟩
  | x :: t, h => by
    obtain ⟨hro, hU, hD⟩ := h x List.mem_cons_self
    obtain ⟨ihg, ihs⟩ := tsub_mainF cls cx cy t (fun y hy => h y (List.mem_cons_of_mem _ hy))
    obtain ⟨gU, _⟩ := opt_good (Q := Q) (cap := capU S Q cls x.1) (cx := cx) (cy := cy) hro hU
    obtain ⟨gD, _⟩ := opt_good (Q := Q) (cap := capD S Q cls x.1) (cx := cx) (cy := cy) hro hD
    have g12 := goodT_tor gU gD
    have gall := goodT_tor g12 ihg
    refine ⟨by simp only [mainF]; exact gall, fun y hy => ?_⟩
    simp only [mainF]
    have k1 := tsub_tor_l _ _ g12.1 ihg.1
    have k2 := tsub_tor_r _ _ g12.1 ihg.1
    rcases List.mem_cons.mp hy with rfl | hy
    · exact ⟨tsub_trans (tsub_tor_l _ _ gU.1 gD.1) k1, tsub_trans (tsub_tor_r _ _ gU.1 gD.1) k1⟩
    · obtain ⟨a1, a2⟩ := ihs y hy
      exact ⟨tsub_trans a1 k2, tsub_trans a2 k2⟩

/-- **The main check**, read back: at the corner, `Σ (chosen − σ) ≥ lg` on the bin. -/
lemma main_fact {Q R S U0 U1 : ℕ} (hQ : 0 < Q) (hR : 0 < R) (hU01 : U0 ≤ U1)
    (cls : List (SegE × ℕ)) {lg cx cy : ℕ} {lc : List (LLine × (ℕ × ℕ × ℕ × ℕ))}
    (hl : ∀ x ∈ lc, rolesOk R U0 U1 x.1.dir x.1.roles = true ∧ isOpt x.1 true x.2.1 = true ∧
        isOpt x.1 false x.2.2.2.1 = true)
    (h : bOk5 (p5neg (p5add (mainP Q R S (mainF Q R S cls cx cy lc) cls cx cy lc).1
      (p5neg (cmul (mainF Q R S cls cx cy lc) R
        (Nat.add (mainP Q R S (mainF Q R S cls cx cy lc) cls cx cy lc).2 lg, 0))))) U0 U1 = true)
    {u : ℝ} (hu0 : (U0 : ℝ) / R ≤ u) (hu1 : u ≤ (U1 : ℝ) / R) :
    (lg : ℝ) ≤ (lc.map fun x =>
        (oval Q x.1 true (capU S Q cls x.1) x.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u - x.2.2.1) +
        (oval Q x.1 false (capD S Q cls x.1) x.2.2.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u - x.2.2.2.2)).sum := by
  obtain ⟨v0, v1⟩ := v_range hR hu0 hu1
  have hb := bOk5_e5 h hU01 v0 v1
  obtain ⟨hg, hsub⟩ := tsub_mainF (Q := Q) (S := S) cls cx cy lc hl
  set f := mainF Q R S cls cx cy lc
  obtain ⟨e1, e2⟩ := mainP_e5 (S := S) hQ hR cls cx cy hu0 hu1 hg.1 lc
    (fun x hx => ⟨(hl x hx).1, (hl x hx).2.1, (hl x hx).2.2, (hsub x hx).1, (hsub x hx).2⟩)
  simp only [e5_neg, e5_add] at hb
  rw [e1, e5_cmul f R hg.1, e2] at hb
  have hpos := fval_pos_good hR hg v0 v1
  have hs : ∀ M : List (LLine × (ℕ × ℕ × ℕ × ℕ)),
      (((M.map fun x => x.2.2.1 + x.2.2.2.2).sum : ℕ) : ℝ)
        = (M.map fun x => ((x.2.2.1 : ℝ) + x.2.2.2.2)).sum := by
    intro M
    induction M with
    | nil => simp
    | cons a M ih => simp only [List.map_cons, List.sum_cons, Nat.cast_add, ih]
  simp only [sv_mk, Nat.cast_zero, sub_zero, nat_add_eq, Nat.cast_add] at hb
  rw [hs] at hb
  have key : (lg : ℝ) ≤ (lc.map fun x =>
      oval Q x.1 true (capU S Q cls x.1) x.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u +
      oval Q x.1 false (capD S Q cls x.1) x.2.2.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u).sum
      - (lc.map fun x => ((x.2.2.1 : ℝ) + x.2.2.2.2)).sum := by
    by_contra hc
    push Not at hc
    have := mul_pos hpos (sub_pos.mpr hc)
    linarith
  have e3' : ∀ M : List (LLine × (ℕ × ℕ × ℕ × ℕ)), (M.map fun x =>
      (oval Q x.1 true (capU S Q cls x.1) x.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u - x.2.2.1) +
      (oval Q x.1 false (capD S Q cls x.1) x.2.2.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u - x.2.2.2.2)).sum
      = (M.map fun x =>
      oval Q x.1 true (capU S Q cls x.1) x.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u +
      oval Q x.1 false (capD S Q cls x.1) x.2.2.2.1 ((cx : ℝ) / Q, (cy : ℝ) / Q) u).sum
      - (M.map fun x => ((x.2.2.1 : ℝ) + x.2.2.2.2)).sum := by
    intro M
    induction M with
    | nil => simp
    | cons a M ih => simp only [List.map_cons, List.sum_cons, ih]; ring
  rw [e3']
  exact key
end ZMTreeM

end SquarePacking
