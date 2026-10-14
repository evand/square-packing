"""Reduced contact system of record N: upright squares pinned at k/2 or S - k/2, unknowns = the rest.

  python3 search/exact/explain/symsys.py N [--solve]     (--solve: lex Groebner basis, polynomial in S)
"""
import sympy as sp, math, sys, render, contacts as C
def build(n):
    S = sp.Symbol('S'); Sf, sq = render.load(n)
    angs = []
    for q in sq:
        a = q[2] % (math.pi/2)
        if min(a, math.pi/2-a) > 1e-9 and not any(abs(a-b) < 1e-9 for b in angs): angs.append(a)
    cs = [sp.symbols(f'c{k} s{k}') for k in range(len(angs))]
    cen, U, V, syms = {}, {}, {}, []
    half = sp.Rational(1, 2)
    def ex(v):
        for k in range(1, 2*int(Sf)+2, 2):
            if abs(v - (Sf - k/2)) < 1e-9: return S - sp.Rational(k, 2)
            if abs(v - k/2) < 1e-9: return sp.Rational(k, 2)
    for i, (x, y, t) in enumerate(sq):
        a = t % (math.pi/2)
        k = next((k for k, b in enumerate(angs) if abs(a-b) < 1e-9), None)
        if k is None:
            X, Y = ex(x), ex(y)
            if X is None: X = sp.Symbol(f'x{i}'); syms.append(X)
            if Y is None: Y = sp.Symbol(f'y{i}'); syms.append(Y)
            u, v = (1, 0), (0, 1)
        else:
            X, Y = sp.symbols(f'x{i} y{i}'); syms += [X, Y]
            c, s = cs[k]; u, v = (c, s), (-s, c)   # local frame rotated by a in (0,90)
        cen[i], U[i], V[i] = (X, Y), u, v
    # corners in the local frame of the stored angle a (fold): recompute float corners with angle a, match indices by position
    def corner_sym(i, p):
        X, Y = cen[i]; u, v = U[i], V[i]
        x, y, t = sq[i]; a = t % (math.pi/2); cu, su = math.cos(a), math.sin(a)
        if U[i] == (1, 0): cu, su = 1.0, 0.0
        du = (p[0]-x)*cu + (p[1]-y)*su; dv = -(p[0]-x)*su + (p[1]-y)*cu
        A = half if du > 0 else -half; B = half if dv > 0 else -half
        return (X + A*u[0] + B*v[0], Y + A*u[1] + B*v[1])
    eqs = []
    import jstyle
    for ct in C.contacts(Sf, sq):
        if ct[0] == 'wall':
            _, i, k, w = ct; P = corner_sym(i, jstyle.corners(*sq[i])[k])
            e = {'L': P[0], 'R': S - P[0], 'B': P[1], 'T': S - P[1]}[w]
        else:
            _, i, k, j, ed = ct; p = jstyle.corners(*sq[i])[k]; P = corner_sym(i, p); X, Y = cen[j]
            x, y, t = sq[j]; a = t % (math.pi/2) if U[j] != (1, 0) else 0.0
            du = (p[0]-x)*math.cos(a) + (p[1]-y)*math.sin(a); dv = -(p[0]-x)*math.sin(a) + (p[1]-y)*math.cos(a)
            if abs(abs(du)-.5) < 1e-9: d, sg = U[j], (1 if du > 0 else -1)
            else: d, sg = V[j], (1 if dv > 0 else -1)
            e = (P[0]-X)*d[0] + (P[1]-Y)*d[1] - sg*half
        e = sp.expand(e)
        if e != 0 and e not in eqs and -e not in eqs: eqs.append(e)
    return S, eqs, syms, cs, angs
if __name__ == '__main__':
    n = int(sys.argv[1]); S, eqs, syms, cs, angs = build(n)
    print(n, 'angles', [round(math.degrees(a), 4) for a in angs], len(eqs), 'eqs', len(syms) + 1 + 2*len(cs), 'unknowns (incl S, c_k, s_k)')
    for e in eqs: print('  ', e)
    if '--solve' in sys.argv:
        E = eqs + [c**2 + s**2 - 1 for c, s in cs]
        G = sp.groebner(E, *syms, *[v for p in cs for v in p], S, order='lex')
        for g in G.exprs:
            if g.free_symbols <= {S}: print('S-poly:', sp.factor(g))
