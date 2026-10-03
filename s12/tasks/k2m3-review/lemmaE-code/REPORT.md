# lemmaE-code — review of qx2_zm.py `Exact` + helpers (sha cdade4b6…a1b7), 2026-09-29 (report returned as text; saved by coordinator)

**Verdict: no BREAKS, no GAP.**  3 MINOR.  Scripts: mass.py, sound.py, lam.py, rootB.py, germ.py, vacuous.py. Cores 2–3, ~1.6 CPU-h.

MINOR
1. `lines_in_reach` (qx2_zm.py:403–417) omits lines from 𝓛 by a float test (radius 0.7072 + ½·diag). Sound in effect (Q within √2/2 of centre; ~9e-5 margin ≫ float error), but contradicts §4.4 "floats only choose options, pieces, orderings and the S-procedure multiplier". Wording fix.
2. No guard against an empty alternatives list (l.748–750, 766–793): if any `terms` entry were [], the combination list is empty and the candidate certifies vacuously (vacuous.py: monkeypatched `_pieces` → [] certifies τ = 2 on a wall box). Unreachable in shipped code (`_pieces` returns P[lo_i:hi_i+1] with lo_i ≤ hi_i; piece_index_f monotone, NaN → 0 both ends, loops only widen; _prune_dom/_prune_const never empty a list). Add an assert.
3. §4.4 step 3 docs: omits std-cap d = 0 split (`cap_alts`) and `_gcdcert` (only in §6); S-procedure tries every open/branch constraint, not only "the violated" one; `cmin` = 7071/10⁴ (l.385) computed but unused.

CHECKED OK (every line in scope read, formulas re-derived)
- Bernstein coefficients/bisection; u0 = 0 (divide out u^k, test closed [0,u1]) valid for (0,u1]; signs of C^α S^β N^γ Δ^δ.
- All four chord options per line (§4.3) and every centre-plane line constructor; bad breakpoints incl. support ends; `crosses` needs strict sign at all four corners over whole bin.
- Exhaustion never certifies: Δ ≡ 0 skipped as parallel; uncertain sign → split; depth/budget out → False.
- "Outside R(u)" skip, option dropping (direction right for up/lo), _prune_dom, _prune_const, _pieces bookkeeping (no off-by-one).
- Dropped lines contribute 0 ≤ μ_ℓ.
- Every combination checked; S-procedure ν ≥ 0; _gcdcert sign algebra re-derived.
- Φ terms: parabola ≥ ĝ; tan mode d* ≤ 1 ≤ s + c; corner3 bookkeeping (two cap_term bounds cancel exactly); xcut; both McCormick anchors; monotone bounds in tb/rect/u_regime/corner3_ok/xcut_ok on [u0,u_e]; EXACT45 claims only θ ≤ 45°.

OWN EXACT TESTS (mass.py from scratch: polygon clipping + edge-crossing chords, no reuse of options/Profile; find admissible pose with exact mass m_s, require certify(box, τ = m_s + 1e-12) to fail)
- sound.py: 397 random boxes (real + density-mutated covers; 162 u0 = 0 with u1 down to 2^-14; 33 straddling 45°; every regime): 0 violations; 276 certified within 1e-6 of m_s (sharp).
- lam.py: 22 tight θ = 0 corners outside U: τ = 1 + 1e-12 fails; λ = λ_s + 1e-9 fails.
- rootB.py: 18 sub-boxes of run-B root (d = 0 split 87–2688×/box, _gcdcert 13–190, S-procedure ≤ 2642): 0 violations.
- germ.py: 12 boxes cornered at germs (3/2,3/2), (23/10,3/2): 0 violations.

NOT CHECKED: E′/E″ math and concavity (math reviewer); zm_mixed, clip_bin, CAP/LEB, D4/SYM, Lemma Z; no certificate roots re-run.
