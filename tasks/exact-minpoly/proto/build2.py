"""Contact graph of an exact packing -> polynomial system (msolve format).
Squares: centre (x_i, y_i), rotation (c_i, s_i), c^2+s^2=1 (axis squares: c=1, s=0 substituted).
Incidence: corner a of j on side k of i:  n_ik . (P_ja - C_i) - 1/2 = 0;  corner on wall: P.x = 0 / S, P.y = 0 / S.
Flat squares: one coordinate pinned at a nearby rational (gauge fix)."""
import sys, json
from mpmath import mp, mpf
from fractions import Fraction
import sympy as sp
mp.dps = 2100
f, out = sys.argv[1], sys.argv[2]
skip = set(json.loads(sys.argv[3]))           # rattlers
pin = json.loads(sys.argv[4])                 # {square: 'x'|'y'}
L = [l.split() for l in open(f).read().split('\n') if l.strip()]
n = int(L[0][0]); S0 = mpf(L[0][1])
sq = []
for k, (x, y, t) in enumerate(L[1:]):
    th = mpf(t) * mp.pi / 180
    sq.append((mpf(x), mpf(y), mp.cos(th), mp.sin(th)))
tol = mpf(10) ** -int(sys.argv[5]) if len(sys.argv) > 5 else mpf(10) ** -200
Ssym = sp.Symbol('S')
V = {}; axis = set(); syms = []; classes = []; rep = {}
for i, (x, y, c, s) in enumerate(sq):
    if i in skip: continue
    X, Y = sp.Symbol(f'x{i}'), sp.Symbol(f'y{i}'); syms += [X, Y]
    # axis square if rotation is a multiple of 90 deg: then (c, s) in {(+-1,0),(0,+-1)}; all equivalent -> (1,0)
    if min(abs(c), abs(s)) < tol:
        axis.add(i); V[i] = (X, Y, sp.Integer(1), sp.Integer(0))
    else:
        th = (mpf(L[i+1][2]) % 90)
        for a, (ta, Ca, Sa) in enumerate(classes):
            if abs(th - ta) < tol: V[i] = (X, Y, Ca, Sa); break
            if abs(th - (90 - ta)) < tol: V[i] = (X, Y, Ca, -Sa); break   # angle -ta (mod 90)
        else:
            Ca, Sa = sp.Symbol(f'c{len(classes)}'), sp.Symbol(f's{len(classes)}')
            classes.append((th, Ca, Sa)); V[i] = (X, Y, Ca, Sa); rep[len(classes)-1] = i
def num(i):
    x, y, c, s = sq[i]
    if i in axis: return (x, y, mpf(1), mpf(0))
    for a, (ta, Ca, Sa) in enumerate(classes):
        if V[i][2] is Ca:
            th = ta * mp.pi / 180
            return (x, y, mp.cos(th), mp.sin(th) if V[i][3] is Sa else -mp.sin(th))
corners = [(sp.Rational(a, 2), sp.Rational(b, 2)) for a in (-1, 1) for b in (-1, 1)]
def corner(v, a, b): X, Y, C, Sn = v; return (X + C * a - Sn * b, Y + Sn * a + C * b)
sides = [(1, 0), (-1, 0), (0, 1), (0, -1)]
def normal(v, e): X, Y, C, Sn = v; return (C * e[0] - Sn * e[1], Sn * e[0] + C * e[1])
eqs = []; subs_num = {}
for i in V:
    for sym, val in zip(V[i], num(i)):
        if isinstance(sym, sp.Symbol): subs_num[sym] = val
        elif isinstance(-sym, sp.Symbol): subs_num[-sym] = -val
subs_num[Ssym] = S0
def ev(e): return e.evalf(subs={k: v for k, v in subs_num.items()}, n=50) if False else mp_eval(e)
fl = {str(k): v for k, v in subs_num.items()}
def mp_eval(e):
    return sp.lambdify(list(subs_num.keys()), e, modules='mpmath')(*subs_num.values())
for j in V:
    for (a, b) in corners:
        P = corner(V[j], a, b)
        Pn = [mp_eval(p) if not p.is_number else mpf(p) for p in P]
        for w, e in ((0, P[0]), (1, P[1])):
            if abs(Pn[w]) < tol: eqs.append(('wall', j, e))
            if abs(Pn[w] - S0) < tol: eqs.append(('wall', j, e - Ssym))
        for i in V:
            if i == j: continue
            for ek in sides:
                nrm = normal(V[i], ek)
                g = nrm[0] * (P[0] - V[i][0]) + nrm[1] * (P[1] - V[i][1]) - sp.Rational(1, 2)
                gn = mp_eval(g)
                if abs(gn) < tol:
                    # tangential coordinate within the side
                    tv = (-nrm[1]) * (P[0] - V[i][0]) + nrm[0] * (P[1] - V[i][1])
                    if abs(mp_eval(tv)) <= mpf(1)/2 + tol:
                        eqs.append(('cs', (j, i), sp.expand(g)))
for a, (ta, Ca, Sa) in enumerate(classes): eqs.append(('unit', a, Ca**2 + Sa**2 - 1))
print('angle classes:', [float(c[0]) for c in classes], file=sys.stderr)
polys = [sp.expand(e[2]) for e in eqs]
print('max residual at point:', mp.nstr(max(abs(mp_eval(p)) if not p.is_number else abs(p) for p in polys), 3), file=sys.stderr)
# gauge fix: null space of the Jacobian (all unknowns incl. S) at the point; pin pivot coordinates at nearby rationals
import numpy as np, scipy.linalg as sla
uk = sorted({s for p in polys for s in p.free_symbols}, key=str)
Jf = np.array([[float(mp_eval(sp.diff(p, u))) for u in uk] for p in polys])
U, sv, Vt = np.linalg.svd(Jf)
rank = int((sv > 1e-9 * sv[0]).sum()); N = Vt[rank:].T
print('rank', rank, 'of', len(uk), 'unknowns; null dim', N.shape[1], file=sys.stderr)
if N.shape[1]:
    Q, R, piv = sla.qr(N.T, pivoting=True)
    for k in piv[:N.shape[1]]:
        u = uk[k]; assert str(u) != 'S', 'S not determined by contacts'
        r = Fraction(str(mp.nstr(subs_num[u], 8)))
        print('pin', u, '=', r, file=sys.stderr)
        polys.append(u - sp.Rational(r.numerator, r.denominator))
# dedupe
seen = []; uniq = []
for p in polys:
    if p == 0: continue
    if any(sp.expand(p - q) == 0 or sp.expand(p + q) == 0 for q in seen): continue
    seen.append(p); uniq.append(p)
allsyms = sorted({s for p in uniq for s in p.free_symbols} - {Ssym}, key=str) + [Ssym]
print(len(eqs), 'raw,', len(uniq), 'distinct equations;', len(allsyms), 'unknowns;', len(axis), 'axis squares', file=sys.stderr)
with open(out, 'w') as fh:
    fh.write(','.join(map(str, allsyms)) + '\n0\n')
    out_p = []
    for p in uniq:
        den = sp.lcm([sp.fraction(c)[1] for c in sp.Poly(p, *allsyms).coeffs()])
        out_p.append(str(sp.expand(p * den)).replace('**', '^'))
    fh.write(',\n'.join(out_p) + '\n')
