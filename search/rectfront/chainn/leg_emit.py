#!/usr/bin/env python3
"""Band theorem, legibility version (2026-10-25): pair table -> Lean.

The band conditions split into one left line and one right line (`kap_split`); for each line the extreme is at the
start, a potential breakpoint, or the aligned slope (`StripBandLines`: `line_lo`, `line_u0_hi`, `line_u1_hi`).  So
each condition is a finite table of PAIRS (left point, right point).  This module

  * builds every pair in natural form: F(C, S, H, r; points) >= 0, with the active potential pieces and the signs of
    the max atoms as case hypotheses (points: 0, 1/2, 8/7, x2 = bandX2 H, y1 = bandY1 H, t (C t = S), t' (S t' = C));
  * clears denominators: F * D = P + sum q_i e_i (e_i: the point equations), P a polynomial in (C, S, H, r) only;
  * certifies 0 <= P under the case hypotheses (band_cert LP; frame split at C = 4/5 where needed; the zigzag
    factorisations tight_P1, tight_P2 and the 3-4-5 lemma tight_P3 by hand);
  * emits Lean: lean/Sqpack/StripBandLeaves.lean (one theorem per sub-case: `0 <= P`) and
    lean/Sqpack/StripBandPairs.lean (one theorem per pair, natural statement, proof: case split + leaf + identity).

Run: python3 leg_emit.py cert   (certificates -> runs/rectfront/band/leg2_certs.pkl)
     python3 leg_emit.py emit   (Lean files)"""
import os, sys, math, pickle, itertools, time
import numpy as np
import sympy as sp

C, S, H, r = sp.symbols('C S H r', real=True)
t, tp, X, Y = sp.symbols('t tp X Y', real=True)
lam = 2 * r - H; a = 1 - lam; b = 1 - H / 2; g = a + b
Q = sp.Rational(8, 7); T1 = H / 16 - sp.Rational(1, 10); T0 = a + b * (2 - Q) - T1 * Q
HERE = os.path.dirname(os.path.abspath(__file__))
LEAN = os.path.join(HERE, '../../../lean/Sqpack')
CACHE = os.path.join(HERE, '../../../runs/rectfront/band/leg2_certs.pkl')

# points: name -> (sympy value used in F, (numerator, denominator) for clearing, Lean text, equation)
PTS = {
    '0':   (sp.Integer(0), (0, 1), '0', None),
    'h':   (sp.Rational(1, 2), (1, 2), '(1 / 2)', None),
    'q':   (Q, (8, 7), '(8 / 7)', None),
    'x2':  (X, (2 * b - a, 2 * b), '(bandX2 H)', ('hX', X * (2 * b) - (2 * b - a))),
    'y1':  (Y, (T0, b - T1), '(bandY1 H)', ('hY', Y * (b - T1) - T0)),
    't':   (t, (S, C), 't', ('ht', C * t - S)),
    'tp':  (tp, (C, S), "t'", ('ht', S * tp - C)),
}
PNAME = {'0': '0', 'h': '1/2', 'q': '8/7', 'x2': 'x₂', 'y1': 'y₁', 't': 't', 'tp': "t'"}
EQ_LEAN = {'hX': 'hX', 'hY': 'hY', 't': 'ht', 'tp': "ht'"}
LINES = {'A0': (2 * a, b), 'Lz': (a + 2 * b, -b), 'T': (T0, T1), 'A1': (a, b), 'aT': (a + T0, T1), 'Z': (0, 0)}
# Lean text of each piece's value at a point z, and the lemma rewriting the potential to it
PIECE_LEAN = {
    'A0': ('2 * (1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * {z}', 'bandU0_A0\''),
    'Lz': ('(1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (2 - {z})', 'bandU0_Lz\''),
    'T':  ('(1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + (H / 16 - 1 / 10) * ({z} - 8 / 7)', 'bandU0_T\''),
    'A1': ('(1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * {z}', 'bandU1_O\''),
    'aT': ('(1 - (2 * Real.sqrt 2 - H)) + ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2) * (6 / 7) + '
           '(H / 16 - 1 / 10) * ({z} - 8 / 7))', 'bandU1_aT\''),
}
# pieces and their ranges (as point names; None = open end)
U = {'u0': [('A0', None, 'x2'), ('Lz', 'x2', 'q'), ('T', 'q', None)], 'u1': [('A1', None, 'y1'), ('aT', 'y1', None)],
     'z': [('Z', None, None)]}
# the piece at a fixed point (on a boundary either piece is exact; take the first)
FIXED = {('u0', '0'): 'A0', ('u0', 'h'): 'A0', ('u0', 'x2'): 'A0', ('u0', 'q'): 'Lz',
         ('u1', '0'): 'A1', ('u1', 'h'): 'A1', ('u1', 'q'): 'A1', ('u1', 'y1'): 'A1'}
CONDS = {  # left potential, left points, right potential, right points, kind
    'Sa': ('u0', ['0', 'q', 't'], 'u1', ['0', 'tp', 'y1'], 'step'),
    'Sb': ('u1', ['0', 'q', 't'], 'u0', ['0', 'x2', 'q', 'tp'], 'step'),
    'C0': ('u0', ['h', 'q', 't'], 'u0', ['h', 'q', 'tp'], 'cheap'),
    'C1': ('u1', ['h', 'q', 't'], 'u1', ['h', 'q', 'tp'], 'cheap'),
    'N':  ('z', ['0', 'h', 't'], 'z', ['0', 'tp'], 'n'),
}
ORDER = [t, tp, X, Y]

def eqs_for(zs):
    return [PTS[z][3] for z in zs if PTS[z][3] is not None]

def clear(F, zs):
    """F * D = sum q_i e_i + P; returns D, P, [(eqname, q)]."""
    D = sp.Integer(1)
    for z in zs: D *= PTS[z][1][1]
    eqs = eqs_for(zs)
    FD = sp.expand(F * D)
    if not eqs: return D, FD, []
    G = [e for _, e in eqs]
    qs, rem = sp.reduced(FD, G, *ORDER, C, S, H, r, order='lex')
    assert not (rem.free_symbols & set(ORDER)), (F, rem)
    return D, sp.expand(rem), [(nm, sp.expand(q)) for (nm, _), q in zip(eqs, qs)]

def le_nat(z1, z2):
    """region fact z1 <= z2 as natural expression z2 - z1."""
    return PTS[z2][0] - PTS[z1][0]

def pieces_at(u, z):
    if u == 'z': return [('Z', [])]
    if z in ('t', 'tp'):
        out = []
        for nm, lo, hi in U[u]:
            hyps = []
            if lo: hyps.append(('reg', lo, z))
            if hi: hyps.append(('reg', z, hi))
            out.append((nm, hyps))
        return out
    return [(FIXED[(u, z)], [])]

def kh_nat(z, c, s):
    """(linear part, max argument or None) of kh c s H z."""
    v = PTS[z][0]
    aligned = (z == 't' and (c, s) == (C, S)) or (z == 'tp' and (c, s) == (S, C))
    marg = None if (aligned or z == '0') else c * v - s
    return (s - H / 2) * v, marg

def subcases():
    out = []
    for cond, (uf, xs, ug, ys, kind) in CONDS.items():
        for zx, zy in itertools.product(xs, ys):
            for (px, hx), (py, hy) in itertools.product(pieces_at(uf, zx), pieces_at(ug, zy)):
                fx = LINES[px][0] + LINES[px][1] * PTS[zx][0]
                gy = LINES[py][0] + LINES[py][1] * PTS[zy][0]
                kx, mx = kh_nat(zx, C, S); ky, my = kh_nat(zy, S, C)
                if kind == 'step': F0 = fx + kx + C + S - lam - (gy - ky)
                elif kind == 'cheap': F0 = fx + kx + gy + ky + C + S - lam - 2 * g
                else: F0 = C + S + kx + ky - 1
                hyps = list(hx) + list(hy)
                if kind == 'cheap':
                    if zx == 't': hyps.append(('reg', 'h', 't'))
                    if zy == 'tp': hyps.append(('reg', 'h', 'tp'))
                if kind == 'n' and zx == 't': hyps.append(('reg', 't', 'h'))
                for sx in ([None] if mx is None else [1, 0]):
                    for sy in ([None] if my is None else [1, 0]):
                        F = F0 + (mx if sx else 0) + (my if sy else 0)
                        hs = list(hyps)
                        if sx is not None: hs.append(('max', zx, 'L', sx))
                        if sy is not None: hs.append(('max', zy, 'R', sy))
                        zs = [zx, zy]
                        D, P, qs = clear(F, zs)
                        tag = f"{cond}_{zx}_{zy}_{px}_{py}" + ('' if sx is None else ('p' if sx else 'n')) + \
                              ('' if sy is None else ('p' if sy else 'n'))
                        out.append(dict(cond=cond, kind=kind, zx=zx, zy=zy, px=px, py=py, sx=sx, sy=sy, F=sp.expand(F),
                                        D=D, P=P, qs=qs, hyps=hs, tag=tag))
    return out

def hyp_poly(h):
    """division-free polynomial of a case hypothesis (>= 0), with its clearing identity."""
    if h[0] == 'reg':
        _, z1, z2 = h; E = le_nat(z1, z2); zs = [z1, z2]
    else:
        _, z, side, sg = h
        c, s = (C, S) if side == 'L' else (S, C)
        E = (c * PTS[z][0] - s) * (1 if sg else -1); zs = [z]
    D, P, qs = clear(E, zs)
    return E, D, P, qs

def numcheck(sc, nth=361, nH=5):
    f = sp.lambdify((C, S, H, r), sc['P'])
    hs = [sp.lambdify((C, S, H, r), hyp_poly(h)[2]) for h in sc['hyps']]
    R = math.sqrt(2); best = None
    for Hv in np.linspace((14 * R - 5) / 8, 1.87, nH):
        for th in np.linspace(1e-7, math.pi / 2 - 1e-7, nth):
            c, s = math.cos(th), math.sin(th)
            if all(h(c, s, Hv, R) >= -1e-12 for h in hs):
                v = f(c, s, Hv, R)
                if best is None or v < best[0]: best = (v, Hv, math.degrees(th))
    return best

# ------------------------------------------------------------------ certificates
P1 = C**3 + C**2*S - C**2 + C*S**2 - 2*C*S*r + 2*C*S + S**3 - S**2
P2 = C**3 - C**2*H + C**2*S + C**2 + 2*C*H*S + C*S**2 - 2*C*S*r - 2*C*S - H*S**2 + S**3 + S**2
P3 = (560*C**3 - 315*C**2*H + 560*C**2*S + 56*C**2 + 840*C*H*S + 560*C*S**2 - 1120*C*S*r - 544*C*S
      - 560*H*S**2 + 560*S**3 + 560*S**2) / 560

def tight_kind(P):
    for nm, T in (('tight_P1', P1), ('tight_P2', P2), ('tight_P3', P3)):
        q = sp.cancel(sp.expand(P) / sp.expand(T))
        if q.is_number and q > 0: return nm, q
    return None

def certify_all(degs=(2, 3, 4, 5)):
    sys.path.insert(0, HERE)
    import band_cert as BC
    cache = pickle.load(open(CACHE, 'rb')) if os.path.exists(CACHE) else {}
    fails = []
    for sc in subcases():
        key = (sc['tag'], sp.srepr(sc['P']))
        if key in cache: continue
        hp = [hyp_poly(h)[2] for h in sc['hyps']]
        if numcheck(sc, 181, 3) is None:
            # empty region: certify -1 >= 0 (infeasibility)
            gens = BC.base_gens() + [(f'h{n}', e) for n, e in enumerate(hp)]
            res = None
            for d in degs:
                res = BC.certify(sp.Integer(-1), gens, d)
                if res and res[2] < 1e-7: break
                res = None
            cache[key] = ('empty', (d, [(tuple(gens[k][0] for k in lab), w) for lab, w in res[0]], res[1],
                                    [gn for gn, _ in gens])) if res else ('empty-fail',)
            print(f"  {sc['tag']}: empty region {'certified' if res else 'FAIL'}", flush=True)
            if not res: fails.append(sc['tag'])
            pickle.dump(cache, open(CACHE, 'wb')); continue
        tk = tight_kind(sc['P'])
        if tk:
            cache[key] = ('tight', tk[0], tk[1]); print(f"  {sc['tag']}: {tk[0]} x {tk[1]}", flush=True)
            pickle.dump(cache, open(CACHE, 'wb')); continue
        got = None
        gens = BC.base_gens() + [(f'h{n}', e) for n, e in enumerate(hp)]
        for d in degs:
            res = BC.certify(sc['P'], gens, d)
            if res and res[2] < 1e-7:
                got = ('lin', (d, [(tuple(gens[k][0] for k in lab), w) for lab, w in res[0]], res[1],
                               [gn for gn, _ in gens])); break
        if not got:
            sp_ = BC.certify_split(sc['P'], [(None, e) for e in hp], (4, 5, 6))
            if sp_: got = ('split', sp_)
        if got: cache[key] = got; print(f"  {sc['tag']}: {got[0]}", flush=True)
        else: fails.append(sc['tag']); print(f"  {sc['tag']}: FAIL", flush=True)
        pickle.dump(cache, open(CACHE, 'wb'))
    print('failed:', fails)

if __name__ == '__main__':
    if sys.argv[1:] == ['cert']: certify_all()
    elif sys.argv[1:] == ['count']:
        scs = subcases(); print(len(scs))

# ------------------------------------------------------------------ Lean emission
sys.path.insert(0, HERE)
SQ = sp.Symbol('SQRT2')

def lean_nat(e):
    """natural expression -> Lean text (r = √2, X = bandX2 H, Y = bandY1 H, tp = t')."""
    e = sp.expand(e).subs(r, SQ)
    s = sp.sstr(e).replace('**', '^')
    s = s.replace('SQRT2', '(Real.sqrt 2)').replace('tp', "t'").replace('X', '(bandX2 H)').replace('Y', '(bandY1 H)')
    return f'({s})'

def lean_poly(e):
    s = sp.sstr(sp.expand(e)).replace('**', '^')
    return f'({s})'

EQS = {'hX': X * (2 * b) - (2 * b - a), 'hY': Y * (b - T1) - T0, 'ht': C * t - S, "ht'": S * tp - C}
EQ_OF_PT = {'x2': 'hX', 'y1': 'hY', 't': 'ht', 'tp': "ht'"}

def lin_comb(R, eqnames):
    """coefficients c_i with R + sum c_i e_i = 0 (exactly); text ' + c * name ...'."""
    if not eqnames:
        assert sp.expand(R) == 0, R
        return ''
    G = [EQS[n] for n in eqnames]
    qs, rem = sp.reduced(sp.expand(-R), G, *ORDER, C, S, H, r, order='lex')
    assert sp.expand(rem) == 0, (R, rem)
    out = ''
    for n, q in zip(eqnames, qs):
        q = sp.expand(q)
        if q != 0: out += f' + {lean_nat(q)} * {n}'
    return out

DEN_LEAN = {'0': None, 'h': '(2 : ℝ)', 'q': '(7 : ℝ)', 'x2': '(2 * (1 - H / 2))', 'y1': '((1 - H / 2) - (H / 16 - 1 / 10))',
            't': 'C', 'tp': 'S'}
DEN_POS = {'h': 'two_pos', 'q': '(by norm_num : (0 : ℝ) < 7)', 'x2': '(by linarith [hb.2] : (0 : ℝ) < 2 * (1 - H / 2))',
           'y1': '(by linarith [hb.2] : (0 : ℝ) < (1 - H / 2) - (H / 16 - 1 / 10))', 't': 'hCp', 'tp': 'hSp'}

def den_of(zs):
    """Lean text and positivity proof of the product of denominators of the points (skipping 0)."""
    ds = [z for z in zs if z != '0']
    if not ds: return '(1 : ℝ)', 'one_pos', sp.Integer(1)
    txt = ' * '.join(DEN_LEAN[z] for z in ds); pf = DEN_POS[ds[0]]
    for z in ds[1:]: pf = f'(mul_pos {pf} {DEN_POS[z]})'
    val = sp.Integer(1)
    for z in ds: val *= PTS[z][1][1]
    return f'({txt})', pf, val

def leaf_name(tag):
    return 'bandLeaf_' + tag

def emit_leaf(sc, entry):
    import band_emit as BE
    hp = [hyp_poly(h)[2] for h in sc['hyps']]
    args = ''.join(f' (h{n} : 0 ≤ {lean_poly(e)})' for n, e in enumerate(hp))
    head = f"theorem {leaf_name(sc['tag'])} {BE.BASE_ARGS}{args} :\n    0 ≤ {lean_poly(sc['P'])} := by\n"
    gmap = dict(BE.BASE_GEN); gmap.update({f'h{n}': f'h{n}' for n in range(len(hp))})
    kind = entry[0]
    if kind == 'tight':
        nm, q = entry[1], entry[2]
        call = {'tight_P1': 'tight_P1 hC hS hCS hr2 (by linarith)', 'tight_P2': 'tight_P2 hC hS hCS hr2 (by linarith) (by linarith)',
                'tight_P3': 'tight_P3 (r := r) (H := H) hC hS hCS (by linarith) h87'}[nm]
        return head + f'  have := {call}\n  linarith\n'
    body = BE.GEN_PRELUDE
    if kind in ('lin', 'empty'): body += BE.cert_text(entry[1], gmap, '  ')
    elif kind == 'split': body += BE.split_text(entry[1], gmap, '  ')
    else: raise ValueError(entry)
    return head + body

U_LEAN = {'u0': 'bandU0', 'u1': 'bandU1'}

def fixed_piece_proof(u, z):
    """(Lean term) bandU H z = piece value, for a fixed point."""
    return {('u0', '0'): "bandU0_A0' hb le_rfl (bandX2_nonneg hb)", ('u0', 'h'): "bandU0_A0' hb (by norm_num) (half_le_bandX2 hb)",
            ('u0', 'x2'): "bandU0_A0' hb (bandX2_nonneg hb) le_rfl", ('u0', 'q'): "bandU0_Lz' hb (bandX2_le hb) le_rfl",
            ('u1', '0'): "bandU1_O' hb le_rfl (bandY1_nonneg hb)", ('u1', 'h'): "bandU1_O' hb (by norm_num) (by linarith [le_bandY1 hb])",
            ('u1', 'q'): "bandU1_O' hb (by norm_num) (le_bandY1 hb)", ('u1', 'y1'): "bandU1_O' hb (bandY1_nonneg hb) le_rfl"}[(u, z)]

def pair_statement(cond, zx, zy):
    uf, xs, ug, ys, kind = CONDS[cond]
    ZX, ZY = PTS[zx][2], PTS[zy][2]
    args = '(hb : InBand H) (hC : 0 ≤ C) (hS : 0 ≤ S) (hCS : C ^ 2 + S ^ 2 = 1)'
    if zx == 't':
        args += ' {t : ℝ} (ht : C * t = S) (h0t : 0 ≤ t)'
        if kind == 'cheap': args += ' (hht : 1 / 2 ≤ t)'
        if kind == 'n': args += ' (hht : t ≤ 1 / 2)'
    if zy == 'tp':
        args += " {t' : ℝ} (ht' : S * t' = C) (h0t' : 0 ≤ t')"
        if kind == 'cheap': args += " (hht' : 1 / 2 ≤ t')"
    base = '(C + S - (2 * Real.sqrt 2 - H))'
    if kind == 'step':
        concl = f'{U_LEAN[ug]} H {ZY} - kh S C H {ZY} ≤ {U_LEAN[uf]} H {ZX} + kh C S H {ZX} + {base}'
    elif kind == 'cheap':
        concl = (f'2 * ((1 - (2 * Real.sqrt 2 - H)) + (1 - H / 2)) ≤ {U_LEAN[uf]} H {ZX} + kh C S H {ZX} + '
                 f'({U_LEAN[ug]} H {ZY} + kh S C H {ZY}) + {base}')
    else:
        concl = f'1 ≤ C + S + kh C S H {ZX} + kh S C H {ZY}'
    return args, concl

def pair_name(cond, zx, zy): return f'bandPair_{cond}_{zx}_{zy}'

def emit_pair(cond, zx, zy, scs):
    uf, xs, ug, ys, kind = CONDS[cond]
    args, concl = pair_statement(cond, zx, zy)
    out = [f"/-- {cond}: left line at {PNAME[zx]}, right line at {PNAME[zy]}. -/\n"
           f"theorem {pair_name(cond, zx, zy)} {args} :\n    {concl} := by\n"]
    out.append("  have r1 := sqrt2_lo'; have r2 := sqrt2_hi'; have r3 := sqrt2_sq\n")
    if 'x2' in (zx, zy) or (zx == 't' and uf == 'u0') or (zy == 'tp' and ug == 'u0'): out.append('  have hX := bandX2_eq hb\n')
    if 'y1' in (zx, zy) or (zx == 't' and uf == 'u1') or (zy == 'tp' and ug == 'u1'): out.append('  have hY := bandY1_eq hb\n')
    if zx == 't': out.append('  have hCp : 0 < C := lt_of_le_of_ne hC (by rintro rfl; nlinarith)\n')
    if zy == 'tp': out.append("  have hSp : 0 < S := lt_of_le_of_ne hS (by rintro rfl; nlinarith)\n")
    ZX, ZY = PTS[zx][2], PTS[zy][2]
    mysc = [sc for sc in scs if sc['cond'] == cond and sc['zx'] == zx and sc['zy'] == zy]
    # case tree: x piece -> y piece -> x max sign -> y max sign
    def piece_tree(u, z, var, gname):
        """list of (piece, [region fact names], lines-generator)"""
        if u == 'z': return [('Z', [], [])]
        if z not in ('t', 'tp'): return [(FIXED[(u, z)], [], [])]
        h0 = 'h0t' if z == 't' else "h0t'"
        if u == 'u0':
            return [('A0', [f'{gname}1'], [f'rcases le_total {var} (bandX2 H) with {gname}1 | {gname}1']),
                    ('Lz', [f'{gname}1', f'{gname}2'], [f'rcases le_total {var} (8 / 7) with {gname}2 | {gname}2']),
                    ('T', [f'{gname}2'], [])]
        return [('A1', [f'{gname}1'], [f'rcases le_total {var} (bandY1 H) with {gname}1 | {gname}1']),
                ('aT', [f'{gname}1'], [])]
    def piece_eq(u, z, var, piece, gname):
        if u == 'z': return None
        txt = PIECE_LEAN[piece][0].format(z=PTS[z][2])
        if z not in ('t', 'tp'): pf = fixed_piece_proof(u, z)
        else:
            h0 = 'h0t' if z == 't' else "h0t'"
            pf = {'A0': f"bandU0_A0' hb {h0} {gname}1", 'Lz': f"bandU0_Lz' hb {gname}1 {gname}2", 'T': f"bandU0_T' hb {gname}2",
                  'A1': f"bandU1_O' hb {h0} {gname}1", 'aT': f"bandU1_aT' hb {gname}1"}[piece]
        return f'{U_LEAN[u]} H {PTS[z][2]} = {txt} := {pf}'
    def kh_eq(z, side, sign, mname):
        c, s = ('C', 'S') if side == 'L' else ('S', 'C')
        Z = PTS[z][2]
        if (z == 't' and side == 'L') or (z == 'tp' and side == 'R'):
            return f'kh {c} {s} H {Z} = ({s} - H / 2) * {Z} := kh_aligned {"ht" if z == "t" else chr(104) + "t" + chr(39)}'
        if z == '0':
            return f'kh {c} {s} H {Z} = ({s} - H / 2) * {Z} := kh_neg (by rw [mul_zero]; exact h{s})'
        if sign: return f'kh {c} {s} H {Z} = ({s} - H / 2) * {Z} + ({c} * {Z} - {s}) := kh_pos {mname}'
        return f'kh {c} {s} H {Z} = ({s} - H / 2) * {Z} := kh_neg {mname}'
    def region_facts(sc):
        """Lean proof terms for each leaf hypothesis."""
        terms = []
        for h in sc['hyps']:
            if h[0] == 'reg':
                _, z1, z2 = h
                fact = region_fact_name(sc, z1, z2)
                E = le_nat(z1, z2); zs = [z1, z2]
            else:
                _, z, side, sg = h
                fact = ('mx' if side == 'L' else 'my')
                c, s = (C, S) if side == 'L' else (S, C)
                E = (c * PTS[z][0] - s) * (1 if sg else -1); zs = [z]
                fact = f'(sub_nonneg.2 {fact})'
            Dt, Dpf, Dv = den_of(zs)
            Ph = hyp_poly(h)[2]
            if h[0] == 'reg': fact = f'(sub_nonneg.2 {fact})'
            hE = f'(mul_nonneg {fact} {Dpf}.le)'
            # goal 0 <= Ph, hyp 0 <= E * Dv
            R = sp.expand(Ph - E * Dv)
            eqn = [EQ_OF_PT[z] for z in zs if z in EQ_OF_PT]
            terms.append(f'(by linear_combination {hE}{lin_comb(R, eqn)})')
        return terms
    def region_fact_name(sc, z1, z2):
        # facts come from the case split names gx1/gx2 (x side, variable t) and gy1/gy2 (y side, variable t')
        tbl = {('t', 'x2'): 'gx1', ('x2', 't'): 'gx1', ('t', 'q'): 'gx2', ('q', 't'): 'gx2', ('t', 'y1'): 'gx1', ('y1', 't'): 'gx1',
               ('tp', 'x2'): 'gy1', ('x2', 'tp'): 'gy1', ('tp', 'q'): 'gy2', ('q', 'tp'): 'gy2', ('tp', 'y1'): 'gy1', ('y1', 'tp'): 'gy1',
               ('h', 't'): 'hht', ('t', 'h'): 'hht', ('h', 'tp'): "hht'"}
        return tbl[(z1, z2)]
    def leaf_block(sc, ind):
        lines = []
        for u, z, piece, gname in ((uf, zx, sc['px'], 'gx'), (ug, zy, sc['py'], 'gy')):
            e = piece_eq(u, z, 't' if z == 't' else "t'", piece, gname)
            if e: lines.append(f'have hu{gname[-1]} : {e}')
        lines.append(f"have hkx : {kh_eq(zx, 'L', sc['sx'], 'mx')}")
        lines.append(f"have hky : {kh_eq(zy, 'R', sc['sy'], 'my')}")
        facts = region_facts(sc)
        lines.append(f"have hP := {leaf_name(sc['tag'])} C S H (Real.sqrt 2) hC hS hCS sqrt2_sq sqrt2_lo' sqrt2_hi' hb.1 hb.2 " + ' '.join(facts))
        Dt, Dpf, Dv = den_of([zx, zy])
        R = sp.expand(sc['F'] * Dv - sc['P'])
        eqn = [EQ_OF_PT[z] for z in (zx, zy) if z in EQ_OF_PT]
        lines.append(f"have key : 0 ≤ {lean_nat(sc['F'])} * {Dt} := by linear_combination hP{lin_comb(R, eqn)}")
        lines.append(f"have hF := nonneg_of_mul_nonneg_left key {Dpf}")
        lines.append('linarith')
        return ''.join(f'{ind}{l}\n' for l in lines)
    # nested case tree
    xt = piece_tree(uf, zx, 't', 'gx'); yt = piece_tree(ug, zy, "t'", 'gy')
    def maxsplit(z, side):
        if z == '0' or (z == 't' and side == 'L') or (z == 'tp' and side == 'R'): return None
        c, s = ('C', 'S') if side == 'L' else ('S', 'C')
        return f'rcases le_total {c} * {PTS[z][2]} with' if False else (c, s)
    def find(px, py, sx, sy):
        for sc in mysc:
            if sc['px'] == px and sc['py'] == py and sc['sx'] == sx and sc['sy'] == sy: return sc
        raise KeyError((cond, zx, zy, px, py, sx, sy))
    def emit_max(px, py, ind):
        msx = maxsplit(zx, 'L'); msy = maxsplit(zy, 'R')
        def rec_y(sx, ind):
            if msy is None:
                out.append(leaf_block(find(px, py, sx, None), ind)); return
            c, s = msy
            out.append(f'{ind}rcases le_total ({c} * {PTS[zy][2]}) {s} with my | my\n')
            out.append(f'{ind}·\n'); out.append(leaf_block(find(px, py, sx, 0), ind + '  '))
            out.append(f'{ind}·\n'); out.append(leaf_block(find(px, py, sx, 1), ind + '  '))
        if msx is None: rec_y(None, ind); return
        c, s = msx
        out.append(f'{ind}rcases le_total ({c} * {PTS[zx][2]}) {s} with mx | mx\n')
        out.append(f'{ind}·\n'); rec_y(0, ind + '  ')
        out.append(f'{ind}·\n'); rec_y(1, ind + '  ')
    def emit_tree(tree, k, ind, cont):
        """tree: [(piece, facts, split-lines)] with splits introducing the next alternative."""
        # implement the 2- or 3-way split as nested rcases
        if len(tree) == 1: cont(tree[0][0], ind); return
        if len(tree) == 2:
            out.append(f'{ind}{tree[0][2][0]}\n')
            out.append(f'{ind}·\n'); cont(tree[0][0], ind + '  ')
            out.append(f'{ind}·\n'); cont(tree[1][0], ind + '  '); return
        out.append(f'{ind}{tree[0][2][0]}\n')
        out.append(f'{ind}·\n'); cont(tree[0][0], ind + '  ')
        out.append(f'{ind}·\n'); out.append(f'{ind}  {tree[1][2][0]}\n')
        out.append(f'{ind}  ·\n'); cont(tree[1][0], ind + '    ')
        out.append(f'{ind}  ·\n'); cont(tree[2][0], ind + '    ')
    def cont_x(px, ind):
        emit_tree(yt, 0, ind, lambda py, ind2: emit_region_ok(px, py, ind2))
    def emit_region_ok(px, py, ind):
        if not any(sc['px'] == px and sc['py'] == py for sc in mysc):
            raise KeyError((cond, zx, zy, px, py))
        emit_max(px, py, ind)
    emit_tree(xt, 0, '  ', cont_x)
    return ''.join(out)

HEADER_LEAVES = """import Sqpack.StripBand

/-! Generated by `search/rectfront/chainn/leg_emit.py`; do not edit by hand.

The polynomial inequalities of the pair table (`StripBandPairs`), one per sub-case: `0 ≤ P(C, S, H, r)` on the frame
`C, S ≥ 0`, `C² + S² = 1`, the band, `r = √2` (bounds), under the case hypotheses `h₀, h₁, …` (which potential piece
is active, the sign of a kink term).  Proofs: `linarith` with the products found by an LP (Handelman-type, with
multiples of the circle relation); three need a split at the 3-4-5 angle `C = 4/5`; the six tight ones are
`tight_P1`, `tight_P2` (the 45° zigzag) and `tight_P3` (3-4-5).  Empty regions: the hypotheses contradict. -/

set_option maxHeartbeats 1000000

namespace SquarePacking

"""

HEADER_PAIRS = """import Sqpack.StripBandLines
import Sqpack.StripBandLeaves

/-! Generated by `search/rectfront/chainn/leg_emit.py`; do not edit by hand.

**The pair table.**  For each band condition and each pair (left line, right line) of candidate slopes from the line
lemmas (`StripBandLines`), the condition at that pair.  Points: `0`, `1/2`, `8/7`, `x₂ = bandX2 H`, `y₁ = bandY1 H`,
and the aligned slopes `t` (`C t = S`) and `t'` (`S t' = C`).  Proof of each: split by the active potential piece and the
sign of each kink term; in each case rewrite the potentials and `kh`, clear the denominators by an exact identity
(`linear_combination` with the point equations), and apply the polynomial leaf (`StripBandLeaves`).

Conditions (`a = 1 - λ`, `β = 1 - H/2`, base `C + S - λ`):
* Sa (step even → odd): `u₁(y) - kh(S, C, y) ≤ u₀(x) + kh(C, S, x) + C + S - λ`;
* Sb (step odd → even): the same with `u₀`, `u₁` exchanged;
* C0, C1 (cheap step, slopes `≥ 1/2`): `2(a + β) ≤ u(x) + kh(C, S, x) + u(y) + kh(S, C, y) + C + S - λ`;
* N (left slope `≤ 1/2`): `1 ≤ C + S + kh(C, S, x) + kh(S, C, y)`. -/

set_option maxHeartbeats 1000000

namespace SquarePacking

"""

def emit_files():
    cache = pickle.load(open(CACHE, 'rb'))
    scs = subcases()
    leaves = [HEADER_LEAVES]
    for sc in scs:
        leaves.append(emit_leaf(sc, cache[(sc['tag'], sp.srepr(sc['P']))]) + '\n')
    leaves.append('end SquarePacking\n')
    open(os.path.join(LEAN, 'StripBandLeaves.lean'), 'w').write(''.join(leaves))
    pairs = [HEADER_PAIRS]
    for cond, (uf, xs, ug, ys, kind) in CONDS.items():
        pairs.append(f'/-! ## {cond} -/\n\n')
        for zx, zy in itertools.product(xs, ys):
            pairs.append(emit_pair(cond, zx, zy, scs) + '\n')
    pairs.append('end SquarePacking\n')
    open(os.path.join(LEAN, 'StripBandPairs.lean'), 'w').write(''.join(pairs))
    print('leaves', len(scs))

if __name__ == '__main__' and sys.argv[1:] == ['emit']:
    emit_files()
