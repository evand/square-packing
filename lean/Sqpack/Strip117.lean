import Sqpack.Strip117Data

/-!
# The wall-strip lemma at `h = 117/100`, and the guillotine corollary

`tasks/n12-chat-lemmas/strip-demo/PROOF.md`: at most three closed, pairwise disjoint unit squares in
`[0,t]²`, `t ≤ 4`, have centre height `≤ 1.17` (past the one-line threshold
`(3√2 - 2)/2 ≈ 1.1213` of `wall_strip_le_three_line`).

The hand-proved part (Lemmas 0–2, Case A) is `Sqpack/StripChain.lean`; the checker and its soundness
proof are `Sqpack/StripCheck.lean`; the certificate (`certs/h1.170_fine.json`: 60 boxes, 33 menu
heights, 3088 leaf claims, minimum chain value `4000147/10⁶`) is `Sqpack/Strip117Data.lean`.  Here
every check is run by the kernel (`decide +kernel`, no `native_decide`):

* `params_ok`, `rule_ok` (K1), `cover_ok` (K2), `flags_ok` (K4), `leaves_ok*` (K3, in four chunks of
  boxes), `chain_ok` (K5);
* `wallStripBound_117 : WallStripBound (117/100)`;
* `guillotine_false_117`: `guillotine_false` at `h = 117/100 ≥ 7/6`.
-/

namespace SquarePacking
namespace Strip117

open StripCheck Strip117Data

set_option maxRecDepth 100000

theorem m_eq : cert.m = 60 := by decide +kernel

theorem params_ok : cert.paramsOK = true := by decide +kernel
theorem rule_ok : cert.ruleOK = true := by decide +kernel
theorem cover_ok : cert.coverOK = true := by decide +kernel
theorem flags_ok : cert.flagsOK = true := by decide +kernel
theorem leaves_ok0 : (List.range' 0 15).all cert.boxLeavesOK = true := by decide +kernel
theorem leaves_ok1 : (List.range' 15 15).all cert.boxLeavesOK = true := by decide +kernel
theorem leaves_ok2 : (List.range' 30 15).all cert.boxLeavesOK = true := by decide +kernel
theorem leaves_ok3 : (List.range' 45 15).all cert.boxLeavesOK = true := by decide +kernel
theorem chain_ok : cert.chainOK = true := by decide +kernel

theorem leaves_ok : ∀ b < cert.m, cert.boxLeavesOK b = true := by
  intro b hb
  rw [m_eq] at hb
  have get := fun (s : ℕ) (h : (List.range' s 15).all cert.boxLeavesOK = true)
      (h1 : s ≤ b) (h2 : b < s + 15) =>
    List.all_eq_true.mp h b (List.mem_range'_1.mpr ⟨h1, h2⟩)
  rcases Nat.lt_or_ge b 15 with h | h
  · exact get 0 leaves_ok0 (by omega) (by omega)
  rcases Nat.lt_or_ge b 30 with h' | h'
  · exact get 15 leaves_ok1 h h'
  rcases Nat.lt_or_ge b 45 with h'' | h''
  · exact get 30 leaves_ok2 h' h''
  · exact get 45 leaves_ok3 h'' (by omega)

/-- **The wall-strip lemma at `h = 117/100`.**  In `[0,t]²`, `t ≤ 4`, at most three closed unit
squares, pairwise disjoint as closed sets, have their centre at height at most `117/100`. -/
theorem wallStripBound_117 : WallStripBound (117 / 100) := by
  have h := cert.wallStrip params_ok rule_ok cover_ok flags_ok leaves_ok chain_ok
  have e : ((cert.hN : ℕ) : ℝ) / K = 117 / 100 := by
    rw [show cert.hN = 1170000 from rfl, K]; norm_num
  rwa [e] at h

end Strip117

open scoped Classical in
/-- **No guillotine `(k, 4, 8-k)` split** of twelve closed, pairwise disjoint unit squares in
`[0,t]²`, `t ≤ 4` (`guillotine_false` at the strip height `117/100 ≥ 7/6`). -/
theorem guillotine_false_117 {t u v : ℝ}
    (ht : t ≤ 4) (hu : 0 ≤ u) (huv : u ≤ v) (hv : v ≤ t)
    (ctr : Fin 12 → ℝ × ℝ) (ang : Fin 12 → ℝ)
    (hdisj : ∀ i j, i ≠ j → Disjoint (sq (ctr i) (ang i) 1) (sq (ctr j) (ang j) 1))
    (hcut : ∀ i, sq (ctr i) (ang i) 1 ⊆ vslab 0 u t ∨ sq (ctr i) (ang i) 1 ⊆ vslab u v t ∨
      sq (ctr i) (ang i) 1 ⊆ vslab v t t)
    (hmid : (Finset.univ.filter (fun i => sq (ctr i) (ang i) 1 ⊆ vslab u v t)).card = 4) :
    False :=
  guillotine_false Strip117.wallStripBound_117 (by norm_num) ht hu huv hv ctr ang hdisj hcut hmid

end SquarePacking
