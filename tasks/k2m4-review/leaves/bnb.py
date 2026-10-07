"""Independent rigorous re-proof of leaves (own argument, own code; nothing from public/).
For a pose box B = [x0,x1]x[y0,y1]x[u0,u1] (u = tan(theta/2)): every Q(c, theta), c in the rectangle, theta in the bin,
contains the rectangle K(B) = centre (xm, ym), axes at theta_m = 2 atan(um) (um = mid u), half-widths
s/2 - (|cos|hx + |sin|hy) and s/2 - (|sin|hx + |cos|hy), s = 1/(1 + 3 delta / 2), delta = u1 - u0 (>= |theta - theta_m|,
since d theta / du <= 2), hx, hy = half-sides of the centre rectangle.  [Proof: rotating a vector of the side-s square by
|phi| <= delta moves it by <= delta s / sqrt2, so it stays in the unit square when s (1/2 + delta/sqrt2) <= 1/2; the
intersection of the translates over the 4 rectangle corners is the slab intersection above.]  So mu(Q) >= mu(K(B)),
computed exactly.  Branch and bound; a sub-box with no admissible pose is dropped (exact test); the u-range is clipped to
the admissible one when the bin is below 45 deg (h = (cos + sin)/2 increasing).
usage: bnb.py SAMPLE.jsonl OUT kinds per_kind node_budget min_margin"""
import json, sys, random
from fractions import Fraction as F
from mass import mass_poly, M

def h(u): return (1 - u * u + 2 * u) / (2 * (1 + u * u))
def below45(u): return u * u + 2 * u - 1 <= 0


def inadmissible(b):
    x0, x1, y0, y1, u0, u1 = b
    Mx = min(x1, y1)
    return h(u0) > Mx and h(u1) > Mx          # h unimodal on [0, 1/2]: min at an end


def clip_u(b):
    x0, x1, y0, y1, u0, u1 = b
    Mx = min(x1, y1)
    if not below45(u1) or h(u1) <= Mx: return b
    lo, hi = u0, u1                          # h(lo) <= Mx < h(hi)
    for _ in range(20):
        mid = (lo + hi) / 2
        if h(mid) <= Mx: lo = mid
        else: hi = mid
    return (x0, x1, y0, y1, u0, hi)          # every u > hi has h(u) > h(hi) > Mx: inadmissible


def lb(b):
    x0, x1, y0, y1, u0, u1 = b
    um = (u0 + u1) / 2; N = 1 + um * um
    c, s_ = (1 - um * um) / N, 2 * um / N
    delta = u1 - u0
    s = 1 / (1 + F(3, 2) * delta)
    hx, hy = (x1 - x0) / 2, (y1 - y0) / 2
    w1 = s / 2 - (abs(c) * hx + abs(s_) * hy)
    w2 = s / 2 - (abs(s_) * hx + abs(c) * hy)
    if w1 <= 0 or w2 <= 0: return F(0)
    xm, ym = (x0 + x1) / 2, (y0 + y1) / 2
    V = []
    for a, bb in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        V.append((xm + c * a * w1 - s_ * bb * w2, ym + s_ * a * w1 + c * bb * w2))
    return mass_poly(V)


def prove(box, budget):
    stack = [box]; nodes = 0
    while stack:
        b = stack.pop()
        if inadmissible(b): continue
        b = clip_u(b)
        nodes += 1
        if nodes > budget: return False, nodes
        if lb(b) >= 1: continue
        x0, x1, y0, y1, u0, u1 = b
        dx, dy, du = x1 - x0, y1 - y0, (u1 - u0) * F(3, 2)
        if du >= dx and du >= dy:
            m = (u0 + u1) / 2; stack += [(x0, x1, y0, y1, u0, m), (x0, x1, y0, y1, m, u1)]
        elif dx >= dy:
            m = (x0 + x1) / 2; stack += [(x0, m, y0, y1, u0, u1), (m, x1, y0, y1, u0, u1)]
        else:
            m = (y0 + y1) / 2; stack += [(x0, x1, y0, m, u0, u1), (x0, x1, m, y1, u0, u1)]
    return True, nodes


if __name__ == '__main__':
    src, out, kinds, per, budget, minm = sys.argv[1], sys.argv[2], sys.argv[3].split(','), int(sys.argv[4]), \
        int(sys.argv[5]), float(sys.argv[6])
    rng = random.Random(7)
    rows = [json.loads(l) for l in open(src)]
    fo = open(out, 'w')
    for k in kinds:
        R = [r for r in rows if r['kind'] == k]
        pos = [r for r in R if float(F(r['min'])) - 1 >= minm]
        slow = [r for r in pos if r['cpu'] >= 100]; fast = [r for r in pos if r['cpu'] < 100]
        pick = rng.sample(slow, min(len(slow), per // 2)); pick += rng.sample(fast, min(len(fast), per - len(pick)))
        nok = 0; nn = 0; fails = []
        for r in pick:
            b = tuple(F(v) for v in r['box'])
            ok, n = prove(b, budget)
            nn += n; nok += ok
            if not ok: fails.append((r['box'], r['min']))
            fo.write(json.dumps(dict(kind=k, box=r['box'], ok=ok, nodes=n, smin=r['min'])) + '\n'); fo.flush()
        zero = len(R) - len(pos)
        print(f'{k}: sampled leaves {len(R)} (sampled margin < {minm}: {zero}); B&B tried {len(pick)} '
              f'(slow roots {min(len(slow), per // 2)}): proved {nok}, budget-out {len(pick) - nok}, nodes {nn}', flush=True)
        for f_ in fails[:5]: print('   budget-out', f_[0], float(F(f_[1])) - 1)
