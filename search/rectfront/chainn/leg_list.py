"""Final inequality list (legibility pass 2026-10-25).
Left line x in X-candidates {X0, 8/7, t = S/C}; right line y in Y-candidates (Sa: {0, t' = C/S, y1};
Sb: {0, x2, 8/7, t'}; cheap: both sides {1/2, 8/7, t}).  For each pair, the potentials at the two points are
unfolded into branches: u appearing with + sign (lower bound) -> every min-branch must hold, one max-branch
may be chosen; with - sign (upper bound) -> one min-branch may be chosen.  Regions: which branch is the actual
value at that point (t vs x2, 8/7; t' vs y1 ...).  Report: per pair, the worst slack and #sub-inequalities."""
import math, itertools, numpy as np
r = math.sqrt(2); H9 = (14*r - 5)/8
def P(H):
    lam = 2*r - H; a = 1 - lam; b = 1 - H/2; T1 = H/16 - .1; T0 = a + 6*b/7 - 8*T1/7
    return dict(lam=lam, a=a, b=b, T0=T0, T1=T1, x2=(2*b - a)/(2*b), y1=T0/(b - T1))
def u0(p, x): return min(2*p['a'] + p['b']*x, max(p['a'] + p['b']*(2 - x), p['T0'] + p['T1']*x))
def u1(p, x): return min(p['a'] + p['b']*x, p['a'] + max(p['a'] + p['b']*(2 - x) - p['a'] + 0*x, 0) * 0 + max(p['a'] + p['b']*(2 - x), p['T0'] + p['T1']*x))
def k(c, s, H, x): return (s - H/2)*x + max(c*x - s, 0)
def run():
    worst = {}
    for H in np.linspace(H9, 1.87, 9):
        p = P(H); base0 = lambda C, S: C + S - p['lam']
        for th in np.linspace(1e-9, math.pi/2 - 1e-9, 1441):
            C, S = math.cos(th), math.sin(th); t = S/C; tp = C/S; base = base0(C, S)
            conds = {
              'Sa': ((u0, [0, 8/7, t]), (u1, [0, tp, p['y1']]), 'step'),
              'Sb': ((u1, [0, 8/7, t]), (u0, [0, p['x2'], 8/7, tp]), 'step'),
              'C0': ((u0, [.5, 8/7, t]), (u0, [.5, 8/7, tp]), 'cheap'),
              'C1': ((u1, [.5, 8/7, t]), (u1, [.5, 8/7, tp]), 'cheap')}
            for nm, ((uf, xs), (ug, ys), kind) in conds.items():
                for i, x in enumerate(xs):
                    for j, y in enumerate(ys):
                        if kind == 'cheap' and (x < .5 or y < .5): continue
                        fx = uf(p, x) + k(C, S, H, x)
                        if kind == 'step': v = fx + base - (ug(p, y) - k(S, C, H, y))
                        else: v = fx + ug(p, y) + k(S, C, H, y) + base - 2*(p['a'] + p['b'])
                        key = (nm, i, j)
                        if key not in worst or v < worst[key][0]: worst[key] = (v, H, math.degrees(th))
    names = {'Sa': (['0', '8/7', 't'], ['0', "t'", 'y1']), 'Sb': (['0', '8/7', 't'], ['0', 'x2', '8/7', "t'"]),
             'C0': (['1/2', '8/7', 't'], ['1/2', '8/7', "t'"]), 'C1': (['1/2', '8/7', 't'], ['1/2', '8/7', "t'"])}
    for (nm, i, j), (v, H, d) in sorted(worst.items()):
        print(f"{nm} x={names[nm][0][i]:4s} y={names[nm][1][j]:3s} min slack {v:+.6f} (H={H:.4f}, {d:.2f} deg)")
u1 = lambda p, x: min(p['a'] + p['b']*x, p['a'] + max(p['a'] + p['b']*(2 - x), p['T0'] + p['T1']*x))
run()
