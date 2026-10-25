import Sqpack.StripBandPairs

/-!
# The band conditions for opposite slopes

Each condition is assembled from the line lemmas (`StripBandLines`) and the pair table (`StripBandPairs`): fix the right
line; for each candidate left slope the right line's extreme is at one of its candidates (`line_u1_hi`, `line_u0_hi`,
or `line_lo` for the cheap step), so the pair table applies; then the left line's extreme is at one of its candidates
(`line_lo`).  For `N` the same with no potentials (`kink_lo`, `kink_lo_inf`).
-/

namespace SquarePacking

section core
variable {C S H x y : ℝ}

lemma bandU0_min {H z : ℝ} (hz : 0 ≤ z) :
    bandU0 H z = min (2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z) (0 + bandL H z) := by
  unfold bandU0; rw [abs_of_nonneg hz, zero_add]

lemma bandU1_min {H z : ℝ} (hz : 0 ≤ z) :
    bandU1 H z = min ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * z) ((1 - (2 * Real.sqrt 2 - H)) + bandL H z) := by
  unfold bandU1; rw [abs_of_nonneg hz]

/-- **Sa** (even line `x`, then odd line `y`, opposite slopes). -/
theorem band_Sa (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H) (hx : 0 ≤ x) (hy : 0 ≤ y) :
    bandU1 H y ≤ bandU0 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) := by
  rw [kap_split hC hS hx hy]
  -- the right line against each candidate left line
  have R0 := line_u1_hi (y := y) hb hC hS hCS (bandPair_Sa_0_0 hb hC hS hCS) (bandPair_Sa_0_y1 hb hC hS hCS)
    (fun _ h0 e => bandPair_Sa_0_tp hb hC hS hCS e h0) hy
  have Rq := line_u1_hi (y := y) hb hC hS hCS (bandPair_Sa_q_0 hb hC hS hCS) (bandPair_Sa_q_y1 hb hC hS hCS)
    (fun _ h0 e => bandPair_Sa_q_tp hb hC hS hCS e h0) hy
  have Rt := fun t (h0 : 0 ≤ t) (e : C * t = S) => line_u1_hi (y := y) hb hC hS hCS
    (bandPair_Sa_t_0 hb hC hS hCS e h0) (bandPair_Sa_t_y1 hb hC hS hCS e h0)
    (fun _ h0' e' => bandPair_Sa_t_tp hb hC hS hCS e h0 e' h0') hy
  -- then the left line
  have := line_lo (p := 2 * (1 - (2 * Real.sqrt 2 - H))) (q := 0) (X₀ := 0)
    (m := bandU1 H y - kh S C H y - (C + S - (2 * Real.sqrt 2 - H))) hb hC hS hCS (by norm_num)
    (by rw [← bandU0_min le_rfl]; linarith) (by rw [← bandU0_min (by norm_num)]; linarith)
    (fun t h0 e => by rw [← bandU0_min h0]; linarith [Rt t h0 e]) hx
  rw [← bandU0_min hx] at this; linarith

/-- **Sb** (odd line `x`, then even line `y`, opposite slopes). -/
theorem band_Sb (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H) (hx : 0 ≤ x) (hy : 0 ≤ y) :
    bandU0 H y ≤ bandU1 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) := by
  rw [kap_split hC hS hx hy]
  have R0 := line_u0_hi (y := y) hb hC hS hCS (bandPair_Sb_0_0 hb hC hS hCS) (bandPair_Sb_0_x2 hb hC hS hCS)
    (bandPair_Sb_0_q hb hC hS hCS) (fun _ h0 e => bandPair_Sb_0_tp hb hC hS hCS e h0) hy
  have Rq := line_u0_hi (y := y) hb hC hS hCS (bandPair_Sb_q_0 hb hC hS hCS) (bandPair_Sb_q_x2 hb hC hS hCS)
    (bandPair_Sb_q_q hb hC hS hCS) (fun _ h0 e => bandPair_Sb_q_tp hb hC hS hCS e h0) hy
  have Rt := fun t (h0 : 0 ≤ t) (e : C * t = S) => line_u0_hi (y := y) hb hC hS hCS
    (bandPair_Sb_t_0 hb hC hS hCS e h0) (bandPair_Sb_t_x2 hb hC hS hCS e h0) (bandPair_Sb_t_q hb hC hS hCS e h0)
    (fun _ h0' e' => bandPair_Sb_t_tp hb hC hS hCS e h0 e' h0') hy
  have := line_lo (p := 1 - (2 * Real.sqrt 2 - H)) (q := 1 - (2 * Real.sqrt 2 - H)) (X₀ := 0)
    (m := bandU0 H y - kh S C H y - (C + S - (2 * Real.sqrt 2 - H))) hb hC hS hCS (by norm_num)
    (by rw [← bandU1_min le_rfl]; linarith) (by rw [← bandU1_min (by norm_num)]; linarith)
    (fun t h0 e => by rw [← bandU1_min h0]; linarith [Rt t h0 e]) hx
  rw [← bandU1_min hx] at this; linarith

/-- **C0** (a cheap step from an even line, both slopes `≥ 1/2`). -/
theorem band_C0 (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hx2 : 1 / 2 ≤ x) (hy2 : 1 / 2 ≤ y) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
      bandU0 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) + bandU0 H y := by
  rw [kap_split hC hS hx hy]
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  -- for each candidate right line, the left line
  have L := fun (w : ℝ) (hw : 0 ≤ w)
      (h1 : 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU0 H (1 / 2) + kh C S H (1 / 2) +
        (bandU0 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H)))
      (h2 : 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU0 H (8 / 7) + kh C S H (8 / 7) +
        (bandU0 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H)))
      (h3 : ∀ t, 1 / 2 ≤ t → C * t = S → 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
        bandU0 H t + kh C S H t + (bandU0 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H))) => by
    have := line_lo (p := 2 * (1 - (2 * Real.sqrt 2 - H))) (q := 0) (X₀ := 1 / 2)
      (m := 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU0 H w + kh S C H w) -
        (C + S - (2 * Real.sqrt 2 - H))) hb hC hS hCS (by norm_num)
      (by rw [← bandU0_min (by norm_num)]; linarith) (by rw [← bandU0_min (by norm_num)]; linarith)
      (fun t h0 e => by rw [← bandU0_min (by linarith)]; linarith [h3 t h0 e]) hx2
    rw [← bandU0_min hx] at this
    show 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU0 H x + kh C S H x) -
      (C + S - (2 * Real.sqrt 2 - H)) ≤ bandU0 H w + kh S C H w
    linarith
  -- then the right line
  have := line_lo (p := 2 * (1 - (2 * Real.sqrt 2 - H))) (q := 0) (X₀ := 1 / 2)
    (m := 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU0 H x + kh C S H x) -
      (C + S - (2 * Real.sqrt 2 - H))) hb hS hC hSC (by norm_num)
    (by rw [← bandU0_min (by norm_num)]
        exact L _ (by norm_num) (bandPair_C0_h_h hb hC hS hCS) (bandPair_C0_q_h hb hC hS hCS)
          (fun t h0 e => bandPair_C0_t_h hb hC hS hCS e (by linarith) h0))
    (by rw [← bandU0_min (by norm_num)]
        exact L _ (by norm_num) (bandPair_C0_h_q hb hC hS hCS) (bandPair_C0_q_q hb hC hS hCS)
          (fun t h0 e => bandPair_C0_t_q hb hC hS hCS e (by linarith) h0))
    (fun t' h0' e' => by
      rw [← bandU0_min (by linarith)]
      exact L _ (by linarith) (bandPair_C0_h_tp hb hC hS hCS e' (by linarith) h0')
        (bandPair_C0_q_tp hb hC hS hCS e' (by linarith) h0')
        (fun t h0 e => bandPair_C0_t_tp hb hC hS hCS e (by linarith) h0 e' (by linarith) h0')) hy2
  rw [← bandU0_min hy] at this; linarith

/-- **C1** (a cheap step from an odd line, both slopes `≥ 1/2`). -/
theorem band_C1 (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hx2 : 1 / 2 ≤ x) (hy2 : 1 / 2 ≤ y) :
    2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
      bandU1 H x + (kap C S H x (-y) - (2 * Real.sqrt 2 - H)) + bandU1 H y := by
  rw [kap_split hC hS hx hy]
  have hSC : S ^ 2 + C ^ 2 = 1 := by linarith
  have L := fun (w : ℝ) (hw : 0 ≤ w)
      (h1 : 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU1 H (1 / 2) + kh C S H (1 / 2) +
        (bandU1 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H)))
      (h2 : 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ bandU1 H (8 / 7) + kh C S H (8 / 7) +
        (bandU1 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H)))
      (h3 : ∀ t, 1 / 2 ≤ t → C * t = S → 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤
        bandU1 H t + kh C S H t + (bandU1 H w + kh S C H w) + (C + S - (2 * Real.sqrt 2 - H))) => by
    have := line_lo (p := 1 - (2 * Real.sqrt 2 - H)) (q := 1 - (2 * Real.sqrt 2 - H)) (X₀ := 1 / 2)
      (m := 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU1 H w + kh S C H w) -
        (C + S - (2 * Real.sqrt 2 - H))) hb hC hS hCS (by norm_num)
      (by rw [← bandU1_min (by norm_num)]; linarith) (by rw [← bandU1_min (by norm_num)]; linarith)
      (fun t h0 e => by rw [← bandU1_min (by linarith)]; linarith [h3 t h0 e]) hx2
    rw [← bandU1_min hx] at this
    show 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU1 H x + kh C S H x) -
      (C + S - (2 * Real.sqrt 2 - H)) ≤ bandU1 H w + kh S C H w
    linarith
  have := line_lo (p := 1 - (2 * Real.sqrt 2 - H)) (q := 1 - (2 * Real.sqrt 2 - H)) (X₀ := 1 / 2)
    (m := 2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) - (bandU1 H x + kh C S H x) -
      (C + S - (2 * Real.sqrt 2 - H))) hb hS hC hSC (by norm_num)
    (by rw [← bandU1_min (by norm_num)]
        exact L _ (by norm_num) (bandPair_C1_h_h hb hC hS hCS) (bandPair_C1_q_h hb hC hS hCS)
          (fun t h0 e => bandPair_C1_t_h hb hC hS hCS e (by linarith) h0))
    (by rw [← bandU1_min (by norm_num)]
        exact L _ (by norm_num) (bandPair_C1_h_q hb hC hS hCS) (bandPair_C1_q_q hb hC hS hCS)
          (fun t h0 e => bandPair_C1_t_q hb hC hS hCS e (by linarith) h0))
    (fun t' h0' e' => by
      rw [← bandU1_min (by linarith)]
      exact L _ (by linarith) (bandPair_C1_h_tp hb hC hS hCS e' (by linarith) h0')
        (bandPair_C1_q_tp hb hC hS hCS e' (by linarith) h0')
        (fun t h0 e => bandPair_C1_t_tp hb hC hS hCS e (by linarith) h0 e' (by linarith) h0')) hy2
  rw [← bandU1_min hy] at this; linarith

/-- `kh` on `[lo, hi]`: least at an end or at the aligned slope. -/
lemma kh_lo {c s H lo hi z m : ℝ} (hc : 0 ≤ c) (hl : lo ≤ z) (hh : z ≤ hi) (h1 : m ≤ kh c s H lo)
    (h2 : m ≤ kh c s H hi) (hk : ∀ t, lo ≤ t → t ≤ hi → c * t = s → m ≤ kh c s H t) : m ≤ kh c s H z := by
  have := kink_lo (A := 0) (B := s - H / 2) (s := s) (m := m) hc hl hh (by unfold kh at h1; linarith) (by unfold kh at h2; linarith)
    (fun t h h' e => by have := hk t h h' e; rw [kh_aligned e] at this; linarith)
  unfold kh; linarith

/-- `kh` on `[lo, ∞)` with slope `s - H/2 + c ≥ 0` after the kink: least at `lo` or at the aligned slope. -/
lemma kh_lo_inf {c s H lo z m : ℝ} (hc : 0 ≤ c) (hl : lo ≤ z) (hs : 0 ≤ (s - H / 2) + c)
    (h1 : m ≤ kh c s H lo) (hk : ∀ t, lo ≤ t → c * t = s → m ≤ kh c s H t) : m ≤ kh c s H z := by
  have := kink_lo_inf (A := 0) (B := s - H / 2) (s := s) (m := m) hc hl hs (by unfold kh at h1; linarith)
    (fun t h e => by have := hk t h e; rw [kh_aligned e] at this; linarith)
  unfold kh; linarith

/-- **N**: a left slope `≤ 1/2` makes the step cost at least a grid unit (so such a step is not cheap). -/
theorem band_N (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1) (hb : InBand H) (hx : 0 ≤ x) (hy : 0 ≤ y)
    (hx2 : x ≤ 1 / 2) : 1 ≤ kap C S H x (-y) := by
  rw [kap_split hC hS hx hy]
  have hp1 : 1 ≤ C + S := by nlinarith [mul_nonneg hC hS]
  -- for each candidate left slope, the right line
  have R : ∀ z, 1 ≤ C + S + kh C S H z + kh S C H 0 →
      (∀ t', 0 ≤ t' → S * t' = C → 1 ≤ C + S + kh C S H z + kh S C H t') →
      1 ≤ C + S + kh C S H z + kh S C H y := fun z h0 h1 => by
    have := kh_lo_inf (c := S) (s := C) (H := H) (lo := 0) (z := y) (m := 1 - C - S - kh C S H z) hS hy
      (by linarith [hb.2]) (by linarith) (fun t' h e => by linarith [h1 t' h e])
    linarith
  -- then the left line
  have := kh_lo (c := C) (s := S) (H := H) (lo := 0) (hi := 1 / 2) (z := x) (m := 1 - C - S - kh S C H y) hC hx hx2
    (by linarith [R 0 (bandPair_N_0_0 hb hC hS hCS) (fun _ h0 e => bandPair_N_0_tp hb hC hS hCS e h0)])
    (by linarith [R (1 / 2) (bandPair_N_h_0 hb hC hS hCS) (fun _ h0 e => bandPair_N_h_tp hb hC hS hCS e h0)])
    (fun t h0 h1 e => by
      linarith [R t (bandPair_N_t_0 hb hC hS hCS e h0 h1) (fun _ h0' e' => bandPair_N_t_tp hb hC hS hCS e h0 h1 e' h0')])
  linarith

end core

end SquarePacking
