"""Falsification sweep over dumped leaves: exact masses (own evaluator) at structured + random poses of each leaf box,
plus a float pattern-search minimisation per leaf whose minimiser is re-evaluated exactly.
usage: sample.py OUT.jsonl [kinds] [maxper] [seed]"""
import json, sys, random, math
from fractions import Fraction as F
from mass import mass, mass_f, M

REC = '/home/evand/math/square-packing/public/s12/runs/qx2_k4x_k008/qxzm_full.jsonl'
out = sys.argv[1]
kinds = sys.argv[2].split(',') if len(sys.argv) > 2 else ['EXACT', 'EXACT0', 'EXACT45', 'CAP', 'LEB', 'PIECE']
maxper = int(sys.argv[3]) if len(sys.argv) > 3 else 10 ** 9
rng = random.Random(int(sys.argv[4]) if len(sys.argv) > 4 else 1)


def h(u): return (1 - u * u + 2 * u) / (2 * (1 + u * u))


leaves = {k: [] for k in kinds}
fh = open(REC); fh.readline()
for ln in fh:
    d = json.loads(ln)
    cpu = d['st']['cpu']
    for b, k in d['leaves']:
        if k in leaves: leaves[k].append((tuple(F(v) for v in b), cpu))
for k in kinds:
    L = leaves[k]
    if len(L) > maxper:
        slow = [x for x in L if x[1] >= 100]; fast = [x for x in L if x[1] < 100]
        ns = min(len(slow), maxper // 2)
        L = rng.sample(slow, ns) + rng.sample(fast, min(len(fast), maxper - ns))
    leaves[k] = L
    print(k, len(L), flush=True)


def rat(x, den=10 ** 9): return F(round(x * den), den)


def poses(b):
    x0, x1, y0, y1, u0, u1 = b
    us = [u0, u1, (u0 + u1) / 2, u0 + (u1 - u0) * rat(rng.random(), 1000)]
    if u0 == 0:
        us = [u1 * F(1, 2 ** k) for k in (0, 1, 3, 6, 10, 16, 24, 32, 40)] + [u1 * rat(rng.random(), 1000)]
    xs = [x0, x1, (x0 + x1) / 2]; ys = [y0, y1, (y0 + y1) / 2]
    for u in us:
        hh = h(u)
        for x in xs:
            for y in ys:
                yield max(x, hh), max(y, hh), u
        yield max(x0, hh), max(y0, hh), u        # wall-pressed
        for _ in range(2):
            yield (max(x0 + (x1 - x0) * rat(rng.random(), 997), hh), max(y0 + (y1 - y0) * rat(rng.random(), 991), hh), u)


def ok_pose(x, y, u, b):
    hh = h(u)
    return b[0] <= x <= b[1] and b[2] <= y <= b[3] and hh <= x <= M - hh and hh <= y <= M - hh


def local_min(b, starts=3, iters=50):
    x0, x1, y0, y1, u0, u1 = (float(v) for v in b)
    best = None
    def clampf(x, y, u):
        u = min(max(u, u0), u1); hh = (1 - u * u + 2 * u) / (2 * (1 + u * u))
        x = min(max(x, x0, hh), x1); y = min(max(y, y0, hh), y1)
        return x, y, u
    def f(p):
        x, y, u = clampf(*p)
        if x < (1 - u * u + 2 * u) / (2 * (1 + u * u)) - 1e-15 or y < (1 - u * u + 2 * u) / (2 * (1 + u * u)) - 1e-15:
            return 9.0, (x, y, u)
        return mass_f(x, y, 2 * math.atan(u)), (x, y, u)
    for s in range(starts):
        p = (x0 + (x1 - x0) * rng.random(), y0 + (y1 - y0) * rng.random(), u0 + (u1 - u0) * rng.random())
        if s == 0 and u0 == 0: p = (p[0], p[1], u1 * 1e-3)
        v, p = f(p)
        st = [(x1 - x0) / 4, (y1 - y0) / 4, (u1 - u0) / 4]
        for it in range(iters):
            imp = False
            for ax in range(3):
                for sg in (1, -1):
                    q = list(p); q[ax] += sg * st[ax]
                    vq, q = f(q)
                    if vq < v: v, p, imp = vq, q, True
            if not imp:
                st = [t / 2 for t in st]
        if best is None or v < best[0]: best = (v, p)
    return best


fo = open(out, 'w')
for k in kinds:
    mn = None; n = 0; below = 0; margins = []
    for b, cpu in leaves[k]:
        lm = None
        for x, y, u in poses(b):
            if not ok_pose(x, y, u, b): continue
            v = mass(x, y, u); n += 1
            if lm is None or v < lm[0]: lm = (v, (x, y, u))
        fv, p = local_min(b)
        if fv < 9:
            x, y, u = rat(p[0], 10 ** 12), rat(p[1], 10 ** 12), rat(p[2], 10 ** 12)
            x = min(max(x, b[0]), b[1]); y = min(max(y, b[2]), b[3]); u = min(max(u, b[4]), b[5])
            hh = h(u); x = max(x, hh); y = max(y, hh)
            if ok_pose(x, y, u, b):
                v = mass(x, y, u); n += 1
                if lm is None or v < lm[0]: lm = (v, (x, y, u))
        if lm is None: continue
        if lm[0] < 1: below += 1; print('*** BELOW 1', k, b, lm, flush=True)
        margins.append(float(lm[0] - 1))
        fo.write(json.dumps(dict(kind=k, box=[str(v) for v in b], cpu=cpu, min=str(lm[0]), at=[str(t) for t in lm[1]])) + '\n')
        if mn is None or lm[0] < mn[0]: mn = (lm[0], b, lm[1])
    margins.sort()
    print(f'{k}: leaves {len(leaves[k])}, exact poses {n}, below 1: {below}, min {float(mn[0]) if mn else None} '
          f'at box {[str(v) for v in mn[1]] if mn else None} pose {[str(t) for t in mn[2]] if mn else None}; '
          f'margin quantiles {[margins[int(q * (len(margins) - 1))] for q in (0, .01, .1, .5)] if margins else None}',
          flush=True)
