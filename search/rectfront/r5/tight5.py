#!/usr/bin/env python3
"""Tight Final5 leaves (2026-10-21): opposite slopes, both frontier branches, middle square aligned.
P is linear in H with P1 = dP/dH <= 0, so P >= P(H = 2).  G = P(H = 2) = (remainder mod the kink factor k)
+ k * Q, with the remainder a positive factor times a square modulo r^2 = 2, and Q of fixed sign (LP).
Low-H leaves: P_lo = P_hi + Dpos (6r - 1 - 4H)."""
import itertools, math, pickle, sys
import sympy as sp
from cert5 import leaf_P, gens_for, certify_split, u, v, H, r, certify
from gen5 import positive_den

def analyse(case, kp):
    N, D, region = leaf_P(case + (not case[1] or case[0], kp) if False else case + (case[1] and not case[0], kp), "hi")
    sign, Dpos, _ = positive_den(sp.factor(D))
    P = sp.expand(N if sign > 0 else -N)
    P1 = sp.expand(sp.diff(P, H))
    G = sp.expand(P.subs(H, 2))
    # the kink factor: the factor of the region polynomial that vanishes at the zigzag
    s2 = math.sqrt(2); z = s2 - 1
    kink_poly = region[-1]
    facs = [f for f, m in sp.factor_list(sp.numer(sp.together(kink_poly)))[1]]
    kf = [f for f in facs if abs(float(f.subs({u: z, v: z}))) < 1e-9]
    assert len(kf) == 1, (case, kp, facs)
    k = kf[0]
    w = v if sp.degree(k, v) == 1 else u
    lead = sp.Poly(k, w).LC()
    m = sp.degree(G, w)
    Q, R = sp.div(sp.Poly(sp.expand(G * lead ** m), w), sp.Poly(k, w))
    R = sp.factor(R.as_expr()); Q = sp.expand(Q.as_expr())
    return dict(P=P, P1=P1, G=G, k=k, w=w, lead=lead, m=m, Q=Q, R=R, region=region, Dpos=Dpos)

if __name__ == '__main__' and len(sys.argv) == 1:
    for e1, e2 in ((True, False), (False, True)):
        for cand in ('a1', 'a2'):
            for kp in (True, False):
                case = (e1, e2, 'P', 'P', cand)
                try:
                    a = analyse(case, kp)
                except AssertionError as ex:
                    print(case, kp, 'no kink factor through the zigzag', ex); continue
                Rs = sp.factor(sp.expand(a['R']).subs(r, sp.sqrt(2)), extension=sp.sqrt(2))
                print(case, kp, 'k =', a['k'], ' lead^m =', a['lead'], '^', a['m'])
                print('   R =', a['R'], '\n   R(sqrt2) =', Rs, flush=True)

# ======================================================================== Lean emission of tight leaves
from cert5 import lean_expr, hyp_name, lean_hint
from gen5 import emit_poly_leaf, L

def cert_hints(Pexpr, G, extra=()):
    """single-box certificate; returns the hint list text (with box generators `extra` available)."""
    for deg, cm in ((6, 1e-3), (8, 1e-3), (6, 0.0), (8, 0.0), (10, 0.0)):
        res = certify(sp.expand(Pexpr), list(G) + list(extra), deg, const_min=cm)
        if res:
            used, _ = res
            Gb = list(G) + list(extra)
            return ', '.join(h for h in (lean_hint(l, Gb) for l, _ in used) if h)
    return None

from cert5 import strip_edges
BAND = (sp.Rational(3, 10), sp.Rational(11, 20))

def nonneg_proof(fac, ind=6):
    """Lean proof of 0 <= fac (product of u, v, 1-u, 1-v, 1+u, 1+v, 1+u^2, 1+v^2 powers)."""
    sp_ = ' ' * ind
    return (f'by\n{sp_}have : 0 ≤ 1 - u := by linarith\n{sp_}have : 0 ≤ 1 - v := by linarith\n'
            f'{sp_}positivity')

def emit_tight(name_hi, name_lo, case, mp, kp):
    N, D, region = leaf_P(case + (mp, kp), 'hi')
    sign, Dpos, _ = positive_den(sp.factor(D))
    P = sp.expand(N if sign > 0 else -N)
    Nl, Dl, _ = leaf_P(case + (mp, kp), 'lo')
    sl, Dposl, _ = positive_den(sp.factor(Dl))
    Plo = sp.expand(Nl if sl > 0 else -Nl)
    assert sp.expand(Dposl - Dpos) == 0
    P1 = sp.expand(sp.diff(P, H)); Gp = sp.expand(P.subs(H, 2))
    fG, Gc = strip_edges(Gp)
    fP1, P1c = strip_edges(-P1)
    k = u * v + u + v - 1
    m = sp.degree(Gc, v)
    Q, R = sp.div(sp.Poly(sp.expand(Gc * (u + 1) ** m), v), sp.Poly(k, v))
    Q = sp.expand(Q.as_expr()); R = sp.expand(R.as_expr())
    fQ, Qc = strip_edges(Q)
    # R = c * (edge factors) * ((u + 1 - r)^2 + 2 - r^2)
    sq = (u + 1 - r)**2 + 2 - r**2
    Rq = sp.cancel(R / sp.expand(sq))
    fR, Rc = strip_edges(Rq)
    assert Rc.is_number and Rc > 0, (case, Rc)
    Ghi = gens_for(region, 'hi'); Glo = gens_for(region, 'lo')
    kg = [g for g in Ghi if g[0] == 'reg1'][0]
    s_k = sp.simplify(kg[1] / sp.expand((u * v - u - v - 1) * k))
    ksign = -1 if s_k > 0 else 1                  # sign of k on the region
    Qs = sp.expand(ksign * Qc)                      # k * Q >= 0  <=>  ksign * Qc >= 0
    SQ = [('sq_u2', (u**2 + 2*u - 1)**2, 'uv'), ('sq_v2', (v**2 + 2*v - 1)**2, 'uv')]
    G0 = [g for g in Ghi if g[2] != 'H']
    a, b = BAND
    lowb = [('3/10-u', a - u, 'uv')]; highb = [('u-11/20', u - b, 'uv')]; midb = [('u-3/10', u - a, 'uv'), ('11/20-u', b - u, 'uv')]
    hP1 = cert_hints(P1c, [g for g in G0 if g[2] != 'r'] + SQ)
    hlo = cert_hints(Gc, G0, lowb); hhi = cert_hints(Gc, G0, highb); hQ = cert_hints(Qs, G0 + SQ, midb)
    if None in (hP1, hlo, hhi, hQ):
        return None, dict(P1=hP1 is None, lo=hlo is None, hi=hhi is None, Q=hQ is None)
    ks = '' if ksign > 0 else '-'
    core = f'{name_hi}_core'
    g0args = ' '.join(hyp_name(g[0]) for g in G0)
    hn = lambda nm: hyp_name(nm)
    txt = f"""set_option maxHeartbeats 2000000 in
theorem {core} (u v H r : ℝ) {' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in G0)}
    (g_2H : 0 ≤ 2 - H) (hr : r ^ 2 = 2) :
    0 ≤ {lean_expr(P)} := by
  have g_sq_u2 : 0 ≤ (u ^ 2 + 2 * u - 1) ^ 2 := sq_nonneg _
  have g_sq_v2 : 0 ≤ (v ^ 2 + 2 * v - 1) ^ 2 := sq_nonneg _
  have hP1c : 0 ≤ {lean_expr(P1c)} := by linarith [{hP1}]
  have hP1 : 0 ≤ {lean_expr(-P1)} := by
    have e : {lean_expr(-P1)} = ({L(fP1)}) * ({lean_expr(P1c)}) := by ring
    rw [e]; exact mul_nonneg ({nonneg_proof(fP1)}) hP1c
  have hGc : 0 ≤ {lean_expr(Gc)} := by
    rcases le_total u ({L(a)}) with hs | hs
    · have {hn('3/10-u')} : 0 ≤ {L(a)} - u := by linarith
      linarith [{hlo}]
    rcases le_total u ({L(b)}) with ht | ht
    · have {hn('u-3/10')} : 0 ≤ u - {L(a)} := by linarith
      have {hn('11/20-u')} : 0 ≤ {L(b)} - u := by linarith
      have hQ : 0 ≤ {lean_expr(Qs)} := by linarith [{hQ}]
      have hneg : u * v - u - v - 1 < 0 := by nlinarith
      have hk : 0 ≤ {ks}(u * v + u + v - 1) := by
        have e : {lean_expr(kg[1])} = ({L(s_k)}) * (u * v - u - v - 1) * (u * v + u + v - 1) := by ring
        rw [e] at g_reg1
        nlinarith
      have hR : 0 ≤ {lean_expr(R)} := by
        have e : {lean_expr(R)} = {L(Rc)} * ({L(fR)}) * ((u + 1 - r) ^ 2 + (2 - r ^ 2)) := by ring
        rw [e, hr, sub_self, add_zero]
        have : 0 ≤ 1 - u := by linarith
        have : 0 ≤ 1 - v := by linarith
        positivity
      have hfQ : 0 ≤ {L(fQ)} := {nonneg_proof(fQ, 8)}
      have e : (u + 1) ^ {m} * ({lean_expr(Gc)}) = {lean_expr(R)} + (u * v + u + v - 1) * (({L(fQ)}) * ({lean_expr(Qc)})) := by ring
      have h1 : 0 ≤ (u + 1) ^ {m} * ({lean_expr(Gc)}) := by
        rw [e]; nlinarith [mul_nonneg hk (mul_nonneg hfQ hQ)]
      exact (mul_nonneg_iff_of_pos_left (by positivity)).1 h1
    · have {hn('u-11/20')} : 0 ≤ u - {L(b)} := by linarith
      linarith [{hhi}]
  have hG : 0 ≤ {lean_expr(Gp)} := by
    have e : {lean_expr(Gp)} = ({L(fG)}) * ({lean_expr(Gc)}) := by ring
    rw [e]; exact mul_nonneg ({nonneg_proof(fG)}) hGc
  have e2 : {lean_expr(P)} = ({lean_expr(Gp)}) + (2 - H) * ({lean_expr(-P1)}) := by ring
  rw [e2]; exact add_nonneg hG (mul_nonneg g_2H hP1)

theorem {name_hi} (u v H r : ℝ) {' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in Ghi)} (hr : r ^ 2 = 2) :
    0 ≤ {lean_expr(P)} :=
  {core} u v H r {g0args} (by linarith) hr

theorem {name_lo} (u v H r : ℝ) {' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_expr(g[1])})' for g in Glo)} (hr : r ^ 2 = 2) :
    0 ≤ {lean_expr(Plo)} := by
  have h := {core} u v H r {g0args} (by linarith) hr
  have hD : 0 ≤ {L(Dpos)} := by
    have : 0 ≤ 1 - u := by linarith
    have : 0 ≤ 1 - v := by linarith
    positivity
  have e : {lean_expr(Plo)} = ({lean_expr(P)}) + ({L(Dpos)}) * (6 * r - 1 - 4 * H) := by ring
  rw [e]; exact add_nonneg h (mul_nonneg hD (by linarith))
"""
    return txt, None

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'emit1':
    t, err = emit_tight('TT_hi', 'TT_lo', (True, False, 'P', 'P', 'a1'), False, True)
    print(err if t is None else 'ok')
    if t: open('/home/evand/math/square-packing/public/lean/Sqpack/TmpTight.lean', 'w').write('import Mathlib\n\nnamespace SquarePacking\n\n' + t + '\nend SquarePacking\n')
