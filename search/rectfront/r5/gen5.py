#!/usr/bin/env python3
"""Generate the Lean proof of Final5 (2026-10-21): per case (eps1, eps2, br1, br2, cand) one theorem
`f5_<tag>` with its polynomial leaves `f5p_<tag>_<sub>`.  Reads certificates from runs/rectfront/r5cert/leaves.jsonl."""
import itertools, json, math, re, sys
import sympy as sp
from cert5 import parts, gens_for, leaf_P, certify, lean_expr, lean_hint, hyp_name, u, v, H, r, CS, H0

U, V = u, v
C1, S1 = CS(u); C2, S2 = CS(v)
LC1, LS1 = '((1 - u ^ 2) / (1 + u ^ 2))', '(2 * u / (1 + u ^ 2))'
LC2, LS2 = '((1 - v ^ 2) / (1 + v ^ 2))', '(2 * v / (1 + v ^ 2))'

def L(e):
    """sympy -> Lean, with r -> Real.sqrt 2."""
    e = sp.sympify(e)
    s = sp.sstr(e, order='none').replace('**', '^')
    return sub_r(s) if r in e.free_symbols else s

def tag(case):
    e1, e2, b1, b2, cand = case
    return f"{'T' if e1 else 'F'}{'T' if e2 else 'F'}{b1}{b2}_{cand}"

def bool_s(b): return 'true' if b else 'false'

def X_lean(e1, e2, cand):
    return {'up0': '((1, 0) : ℝ × ℝ)', 'up90': '((0, 1) : ℝ × ℝ)',
            'a1': f'(alFr {bool_s(e1)} {LC1} {LS1})', 'a2': f'(alFr {bool_s(e2)} {LC2} {LS2})'}[cand]

def B_exprs(e1, e2, b1, b2):
    def bnd(C, S, b): return 2 * S - C * (H - 1) if b == 'U' else (C * (S - C) + 1 + C + S - H) / C
    B1 = bnd(C1, S1, b1) + (C1 * H if e1 else 0)
    B2 = bnd(C2, S2, b2) + (0 if e2 else C2 * H)
    return B1, B2

def pos_factor_proof(f):
    """a positive version of factor f (sympy, in u, v) and a Lean proof term of 0 < it."""
    f = sp.expand(f)
    table = {u: 'hu', v: 'hv'}
    if f in table: return f, table[f]
    if f == 1 - u or f == u - 1: return 1 - u, '(by linarith : (0:ℝ) < 1 - u)'
    if f == 1 - v or f == v - 1: return 1 - v, '(by linarith : (0:ℝ) < 1 - v)'
    if f == u + 1: return f, '(by linarith : (0:ℝ) < u + 1)'
    if f == v + 1: return f, '(by linarith : (0:ℝ) < v + 1)'
    if f.is_number: return abs(f), f'(by norm_num : (0:ℝ) < {L(abs(f))})'
    return f, f'(by positivity : (0:ℝ) < {L(f)})'

def positive_den(D):
    """D (sympy) -> (sign, Dpos sympy, Lean proof of 0 < Dpos)."""
    c, fl = sp.factor_list(D)
    sign = 1 if c > 0 else -1
    terms, proofs = [abs(c)], [f'(by norm_num : (0:ℝ) < {L(abs(c))})']
    for f, m in fl:
        fp, pr = pos_factor_proof(f)
        if sp.expand(fp - f) != 0: sign *= (-1) ** m      # f = -fp
        terms.append(fp ** m); proofs.append(pr if m == 1 else f'(pow_pos {pr} {m})')
    Dpos = sp.Mul(*terms)
    pf = proofs[0]
    for p in proofs[1:]: pf = f'(mul_pos {pf} {p})'
    return sign, Dpos, pf

if __name__ == '__main__' and len(sys.argv) == 1:
    case = (False, False, 'U', 'U', 'a1')
    B1, B2 = B_exprs(*case[:4])
    print('B1 =', L(sp.simplify(B1)))
    for mp, kp in itertools.product((True, False), (True, False)):
        E, region = parts(*case, mp, kp)
        N, D = sp.fraction(sp.together(E - 5))
        print(mp, kp, positive_den(sp.factor(D))[0], L(positive_den(sp.factor(D))[1]))

# ======================================================================== emission
def lean_g(e):           # generator / polynomial text (shared by leaf lemmas and case theorems)
    return lean_expr(e)

def sub_r(s):
    return re.sub(r'\br\b', '(Real.sqrt 2)', s)

def emit_poly_leaf(name, P, G, plan):
    """plan: list of (box, used, Gb).  Statement uses base generators G; the body splits boxes."""
    hyps = ' '.join(f'({hyp_name(g[0])} : 0 ≤ {lean_g(g[1])})' for g in G)
    head = f'theorem {name} (u v H r : ℝ) {hyps} :\n    0 ≤ {lean_g(P)} := by\n'
    def hints(used, Gb):
        hs = [h for h in (lean_hint(l, Gb) for l, _ in used) if h]
        return ', '.join(hs)
    if len(plan) == 1:
        box, used, Gb = plan[0]
        return head + f'  linarith [{hints(used, Gb)}]\n'
    # split: boxes on u and/or v at 1/2
    us = sorted({b[0] for b, _, _ in plan}); vs = sorted({b[1] for b, _, _ in plan})
    body = ''
    def box_body(box, indent):
        _, used, Gb = next(p for p in plan if p[0] == box)
        return ' ' * indent + f'linarith [{hints(used, Gb)}]\n'
    def gname(nm): return hyp_name(nm)
    if len(us) == 2 and len(vs) == 1:
        body += '  rcases le_total u (1 / 2) with hs | hs\n'
        body += f'  · have {gname("1/2-u")} : 0 ≤ 1/2 - u := by linarith\n' + box_body((us[0], vs[0]), 4)
        body += f'  · have {gname("u-1/2")} : 0 ≤ u - 1/2 := by linarith\n' + box_body((us[1], vs[0]), 4)
    elif len(us) == 1 and len(vs) == 2:
        body += '  rcases le_total v (1 / 2) with hs | hs\n'
        body += f'  · have {gname("1/2-v")} : 0 ≤ 1/2 - v := by linarith\n' + box_body((us[0], vs[0]), 4)
        body += f'  · have {gname("v-1/2")} : 0 ≤ v - 1/2 := by linarith\n' + box_body((us[0], vs[1]), 4)
    else:
        body += '  rcases le_total u (1 / 2) with hs | hs <;> rcases le_total v (1 / 2) with ht | ht\n'
        for bu, bv, cu, cv in [(0, 0, '1/2-u', '1/2-v'), (0, 1, '1/2-u', 'v-1/2'), (1, 0, 'u-1/2', '1/2-v'), (1, 1, 'u-1/2', 'v-1/2')]:
            body += (f'  · have {gname(cu)} : 0 ≤ {cu.replace("1/2-u", "1/2 - u").replace("u-1/2", "u - 1/2")} := by linarith\n'
                     f'    have {gname(cv)} : 0 ≤ {cv.replace("1/2-v", "1/2 - v").replace("v-1/2", "v - 1/2")} := by linarith\n'
                     + box_body((us[bu], vs[bv]), 4))
    return head + body

def mu_raw(e1, e2):
    return f'sgnC {bool_s(e2)} {LC2} / {LS2} - sgnC {bool_s(e1)} {LC1} / {LS1}'

def h_stmt(b, e, side):
    """semantic branch hypothesis (after simp of the `if`)."""
    C, S = (LC1, LS1) if side == 1 else (LC2, LS2)
    if side == 1:
        X = f'(s₁ - {C} * H)' if e else 's₁'
    else:
        X = f'({S} * L - s₂)' if e else f'({S} * L - {C} * H - s₂)'
    if b == 'U': return f'2 * {S} - {C} * (H - 1) ≤ {X}'
    return f'{C} * ({S} - {C}) + (1 + {C} + {S} - H) ≤ {C} * {X}'

def region_fact(name, raw, g_expr, sign_pos, K):
    """from the split hypothesis `hyp` on `raw` derive 0 <= g_expr (raw or -raw = K * g)."""
    lhs = f'({raw})' if sign_pos else f'(-({raw}))'
    return (f'      have {name} : 0 ≤ {lean_g(g_expr)} := by\n'
            f'        have e : {lhs} = ({L(K)}) * ({lean_g(g_expr)}) := by\n'
            f'          simp only [sgnC, cr, alFr, Bool.false_eq_true, ↓reduceIte, Prod.fst, Prod.snd]\n'
            f'          field_simp\n          ring\n'
            f'        have hK : (0:ℝ) < {L(K)} := by positivity\n'
            f'        have h0 : 0 ≤ {lhs} := by {"linarith" if sign_pos else "linarith"}\n'
            f'        rw [e] at h0\n        exact (mul_nonneg_iff_of_pos_left hK).1 h0\n')

def emit_case(case, certs, lemma_prefix):
    """case = (e1, e2, b1, b2, cand) with cand in up0, up90, a1, a2.  certs: dict (sub, Hint) -> plan."""
    e1, e2, b1, b2, cand = case
    t = tag(case)
    X = X_lean(e1, e2, cand)
    mr = mu_raw(e1, e2)
    pc = 'up' if cand.startswith('up') else cand
    out = []
    out.append(f'theorem f5_{t} (u v H L s₁ s₂ : ℝ) (hu : 0 < u) (hu1 : u ≤ 1) (hv : 0 < v) (hv1 : v ≤ 1)\n'
               f'    (hH1 : 1 ≤ H) (hH2 : H < 2)\n'
               f'    (h1 : {h_stmt(b1, e1, 1)})\n    (h2 : {h_stmt(b2, e2, 2)})\n'
               f'    (hm : midW ({mr}) (1 / (2 * {LS1})) (1 / (2 * {LS2})) (alFr {bool_s(e1)} {LC1} {LS1})\n'
               f'      (alFr {bool_s(e2)} {LC2} {LS2}) {X}\n'
               f'      (if 0 ≤ {mr} then ({X}.1 + {X}.2) / 2 else H - ({X}.1 + {X}.2) / 2) ≤\n'
               f'      s₂ / {LS2} - s₁ / {LS1}) :\n'
               f'    5 ≤ L ∨ 4 + 6 * Real.sqrt 2 ≤ L + 4 * H := by\n')
    out.append(f'  have hS₁ : 0 < {LS1} := by positivity\n  have hS₂ : 0 < {LS2} := by positivity\n'
               f'  have hC₁0 : 0 ≤ {LC1} := div_nonneg (by nlinarith) (by positivity)\n'
               f'  have hC₂0 : 0 ≤ {LC2} := div_nonneg (by nlinarith) (by positivity)\n')
    B1, B2 = B_exprs(e1, e2, b1, b2)
    # branch bounds
    Btxt = {}
    for side, b, e in ((1, b1, e1), (2, b2, e2)):
        C, S, w = (LC1, LS1, 'u') if side == 1 else (LC2, LS2, 'v')
        hn = f'h{side}'
        tgt = 's₁' if side == 1 else f'{LS2} * L - s₂'
        shift = f' + {C} * H' if (side == 1 and e) or (side == 2 and not e) else ''
        if b == 'P':
            Bt = f'(({C} * ({S} - {C}) + (1 + {C} + {S} - H)) / {C}{shift})'
            out.append(f'  have hC{side} : 0 < {C} := psi_branch_pos hC{side}0 (half_unit {w}) hS{"₁" if side == 1 else "₂"}.le hH2 {hn}\n'
                       f'  have h{w}1\' : {w} < 1 := half_lt_one hC{side} (le_of_lt h{w})\n'
                       f'  have hd{side} := psi_branch_div hC{side} {hn}\n'
                       f'  have b{side} : {Bt} ≤ {tgt} := by linarith\n')
        else:
            Bt = f'(2 * {S} - {C} * (H - 1){shift})'
            out.append(f'  have b{side} : {Bt} ≤ {tgt} := by linarith\n')
        Btxt[side] = Bt
    out.append('  unfold midW at hm\n')
    # abs terms
    subs = []
    if pc == 'up':
        lem = 'cr_e0' if cand == 'up0' else 'cr_e90'
        comp = '2' if cand == 'up0' else '1'
        neg = 'abs_neg, ' if cand == 'up0' else ''
        out.append(f'  rw [{lem}, {lem}, {neg}{neg}abs_of_nonneg (alFr_nonneg (ε := {bool_s(e1)}) hC₁0 hS₁.le).{comp},\n'
                   f'    abs_of_nonneg (alFr_nonneg (ε := {bool_s(e2)}) hC₂0 hS₂.le).{comp}] at hm\n')
        kinks = [(None, None)]
    else:
        if cand == 'a1':
            zero, kink = f'(alFr {bool_s(e1)} {LC1} {LS1})', f'cr (alFr {bool_s(e2)} {LC2} {LS2}) (alFr {bool_s(e1)} {LC1} {LS1})'
        else:
            zero, kink = f'(alFr {bool_s(e2)} {LC2} {LS2})', f'cr (alFr {bool_s(e1)} {LC1} {LS1}) (alFr {bool_s(e2)} {LC2} {LS2})'
        out.append(f'  rw [cr_self, abs_zero] at hm\n  rcases le_total 0 ({kink}) with hk | hk\n')
        kinks = [(True, kink), (False, kink)]
    for kp, kink in kinks:
        ind = '  ' if kp is None else '    '
        if kp is not None:
            out.append(f'  · rw [{"abs_of_nonneg" if kp else "abs_of_nonpos"} hk] at hm\n')
        out.append(f'{ind}split_ifs at hm with hμ\n')
        for mp in (True, False):
            sub = (mp, True if kp is None else kp)
            E, region = parts(e1, e2, b1, b2, pc, mp, sub[1])
            M = E - B1 / S1 - B2 / S2
            G = gens_for(region, 'lo')
            regG = [g for g in G if g[0].startswith('reg')]
            body = []
            # region facts
            mu_e = (C2 if e2 else -C2) / S2 - (C1 if e1 else -C1) / S1
            raws = []
            for g in regG:
                idx = int(g[0][3:])
                if idx == 0:
                    raws.append((g, mr, mp, sp.factor(sp.simplify((mu_e if mp else -mu_e) / g[1]))))
                else:
                    a1 = (S1, C1) if e1 else (C1, S1); a2 = (S2, C2) if e2 else (C2, S2)
                    kexpr = (a2[0] * a1[1] - a2[1] * a1[0]) if cand == 'a1' else (a1[0] * a2[1] - a1[1] * a2[0])
                    raws.append((g, kink, sub[1], sp.factor(sp.simplify((kexpr if sub[1] else -kexpr) / g[1]))))
            for g, raw, sgn, K in raws:
                body.append(region_fact(hyp_name(g[0]), raw, g[1], sgn, K).replace('      ', ind + '  ', 1) if False else region_fact(hyp_name(g[0]), raw, g[1], sgn, K))
            body.append(f'      have hmM : {L(sp.simplify(M))} ≤ s₂ / {LS2} - s₁ / {LS1} := by\n'
                        f'        refine le_trans (le_of_eq ?_) hm\n'
                        f'        simp only [sgnC, cr, dt, alFr, Bool.false_eq_true, ↓reduceIte, Prod.fst, Prod.snd]\n'
                        f'        field_simp\n        ring\n')
            body.append(f'      have key := chain5 hS₁ hS₂ b1 b2 hmM\n'
                        f'      rcases le_or_gt H ((6 * Real.sqrt 2 - 1) / 4) with hh | hh\n')
            for Hint in ('lo', 'hi'):
                key = (sub, Hint)
                Eh, region = parts(e1, e2, b1, b2, pc, mp, sub[1])
                tgt = 5 if Hint == 'lo' else 4 + 6 * r - 4 * H
                N, D = sp.fraction(sp.together(Eh - tgt))
                sign, Dpos, Dpf = positive_den(sp.factor(D))
                P = sp.expand(N if sign > 0 else -N)
                Gh = gens_for(region, Hint)
                lname = f'{lemma_prefix}_{t}_{"p" if mp else "n"}{"" if kp is None else ("p" if kp else "n")}_{Hint}'
                args = []
                for g in Gh:
                    nm = hyp_name(g[0])
                    if g[0].startswith('reg'): args.append(nm)
                    elif g[0] in ('u', 'v'): args.append(f'(le_of_lt h{g[0]})')
                    else: args.append('(by linarith [sqrt2_lo, sqrt2_hi])')
                goal = f'{"5" if Hint == "lo" else "4 + 6 * Real.sqrt 2"} ≤ L{"" if Hint == "lo" else " + 4 * H"}'
                expr = f'{Btxt[1]} / {LS1} + {L(sp.simplify(M))} + {Btxt[2]} / {LS2} - ({L(tgt)})'
                body.append(f'      · {"left" if Hint == "lo" else "right"}\n'
                            f'        have hP := {lname} u v H (Real.sqrt 2) {" ".join(args)}\n'
                            f'        have e : {expr} = ({sub_r(lean_g(P))}) / ({L(Dpos)}) := by\n'
                            f'          field_simp\n          ring\n'
                            f'        have : 0 ≤ {expr} := by rw [e]; exact div_nonneg hP (le_of_lt {Dpf})\n'
                            f'        linarith\n')
                certs.setdefault('need', []).append((lname, P, Gh, case, sub, Hint))
            text = ''.join(body)
            pre = f'{ind}· -- μ {"≥" if mp else "<"} 0\n'
            if not mp:
                text = '      push Not at hμ\n' + text
            out.append(pre + text)
    return ''.join(out)

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'one':
    certs = {}
    txt = emit_case((False, False, 'U', 'U', 'a1'), certs, 'f5p')
    leaves = []
    for lname, P, Gh, case, sub, Hint in certs['need']:
        plan = certify_split(P, Gh) if False else None
        from cert5 import certify_split as cs
        plan = cs(P, Gh)
        if plan is None:
            leaves.append(f'-- FAIL {lname}\ntheorem {lname} : True := trivial\n'); continue
        leaves.append(emit_poly_leaf(lname, P, Gh, plan))
    open('/home/evand/math/square-packing/public/lean/Sqpack/TmpCase.lean', 'w').write(
        'import Sqpack.Rect5Glue\n\nnamespace SquarePacking\n\n' + '\n'.join(leaves) + '\n' + txt + '\nend SquarePacking\n')
    print('written')
