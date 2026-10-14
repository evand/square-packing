"""Full contact system of record N: every centre coordinate, S, and (c, s) per tilt class; Jacobian rank at the packing.

  python3 search/exact/explain/fullsys.py N...
"""
import sympy as sp, numpy as np, math, sys, render, contacts as C, jstyle
def build(n, only=None):
    S = sp.Symbol('S'); Sf, sq = render.load(n)
    angs, rep = [], []
    for q in sq:
        a = q[2] % (math.pi/2)
        if min(a, math.pi/2-a) > 1e-9 and not any(abs(a-b) < 1e-9 for b in angs): angs.append(a); rep.append(q[2])
    cs = [sp.symbols(f'c{k} s{k}') for k in range(len(angs))]
    X = [sp.symbols(f'x{i+1} y{i+1}') for i in range(len(sq))]
    U, V = {}, {}
    rot = lambda w: (-w[1], w[0])
    for i, (x, y, t) in enumerate(sq):
        a = t % (math.pi/2); k = next((k for k, b in enumerate(angs) if abs(a-b) < 1e-9), None)
        if k is None:
            u = (1, 0); m = round(t / (math.pi/2)) % 4
        else:
            u = cs[k]; m = round((t - rep[k]) / (math.pi/2)) % 4
        for _ in range(m): u = rot(u)
        U[i], V[i] = u, rot(u)
    h = sp.Rational(1, 2)
    def loc(i, p):
        x, y, t = sq[i]
        return (p[0]-x)*math.cos(t) + (p[1]-y)*math.sin(t), -(p[0]-x)*math.sin(t) + (p[1]-y)*math.cos(t)
    def corner(i, p):
        du, dv = loc(i, p); A = h if du > 0 else -h; B = h if dv > 0 else -h
        return (X[i][0] + A*U[i][0] + B*V[i][0], X[i][1] + A*U[i][1] + B*V[i][1])
    eqs = []
    for ct in (only if only is not None else C.contacts(Sf, sq)):
        if ct[0] == 'wall':
            _, i, k, w = ct; P = corner(i, jstyle.corners(*sq[i])[k])
            e = {'L': P[0], 'R': S - P[0], 'B': P[1], 'T': S - P[1]}[w]
        else:
            _, i, k, j, ed = ct; p = jstyle.corners(*sq[i])[k]; P = corner(i, p)
            du, dv = loc(j, p)
            d, sg = (U[j], 1 if du > 0 else -1) if ed[0] == 'u' else (V[j], 1 if dv > 0 else -1)
            e = (P[0]-X[j][0])*d[0] + (P[1]-X[j][1])*d[1] - sg*h
        e = sp.expand(e)
        if e != 0 and e not in eqs and -e not in eqs: eqs.append(e)
    circ = [c**2 + s**2 - 1 for c, s in cs]
    unk = [v for p in X for v in p] + [S] + [v for p in cs for v in p]
    val = {S: Sf}
    for i, (x, y, t) in enumerate(sq): val[X[i][0]] = x; val[X[i][1]] = y
    for k, a in enumerate(rep): val[cs[k][0]] = math.cos(a); val[cs[k][1]] = math.sin(a)
    return eqs, circ, unk, val
if __name__ == '__main__':
    for n in map(int, sys.argv[1:]):
        eqs, circ, unk, val = build(n)
        allE = eqs + circ
        res = max(abs(float(e.subs(val))) for e in allE)
        J = np.array([[float(sp.diff(e, v).subs(val)) for v in unk] for e in allE])
        sv = np.linalg.svd(J, compute_uv=False)
        rank = int((sv > 1e-9 * sv[0]).sum())
        print(f'n={n}: unknowns {len(unk)} (2x{n} coords + S + {len(circ)} (c,s) pairs); contact eqs {len(eqs)} distinct + {len(circ)} circle = {len(allE)}; rank {rank}; null dim {len(unk)-rank}; residual {res:.1e}')
