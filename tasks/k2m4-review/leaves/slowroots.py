"""Dense float sweep of the slowest roots (whole root box, admissible poses), then exact re-evaluation of the 20 lowest
float poses after a local refinement.  Reports the min exact mass per root."""
import json, math, random, heapq
from fractions import Fraction as F
from mass import mass, mass_f

d = sorted(json.load(open('root_cpu.json')), key=lambda r: -r[1])[:12]
rng = random.Random(3)
for r, cpu in d:
    x0, x1, y0, y1, u0, u1 = (float(F(v)) for v in r)
    best = []
    n = 24
    for i in range(n + 1):
        for j in range(n + 1):
            for k in range(n + 1):
                u = u0 + (u1 - u0) * k / n
                if u == 0: u = 1e-9
                x = x0 + (x1 - x0) * i / n; y = y0 + (y1 - y0) * j / n
                hh = (1 - u * u + 2 * u) / (2 * (1 + u * u))
                if x < hh or y < hh: continue
                v = mass_f(x, y, 2 * math.atan(u))
                heapq.heappush(best, (-v, (x, y, u)))
                if len(best) > 20: heapq.heappop(best)
    mins = []
    for nv, (x, y, u) in best:
        # pattern search refine
        p = [x, y, u]; v = -nv; st = [(x1 - x0) / n, (y1 - y0) / n, (u1 - u0) / n]
        for it in range(40):
            imp = False
            for ax in range(3):
                for sg in (1, -1):
                    q = list(p); q[ax] = min(max(q[ax] + sg * st[ax], (x0, y0, u0)[ax]), (x1, y1, u1)[ax])
                    if q[2] <= 0: q[2] = 1e-12
                    hh = (1 - q[2] ** 2 + 2 * q[2]) / (2 * (1 + q[2] ** 2))
                    if q[0] < hh or q[1] < hh: continue
                    vq = mass_f(q[0], q[1], 2 * math.atan(q[2]))
                    if vq < v: v, p, imp = vq, q, True
            if not imp: st = [s / 2 for s in st]
        X, Y, U = (F(round(t * 10 ** 12), 10 ** 12) for t in p)
        if U == 0: U = F(1, 10 ** 12)
        mins.append((mass(X, Y, U), (X, Y, U)))
    m = min(mins)
    print([str(v) for v in r], round(cpu), 'min exact', float(m[0]), [float(t) for t in m[1]], 'below 1' if m[0] < 1 else '',
          flush=True)
