"""Second stage for R5 near a = 1 + sqrt2, b -> 3 (and the line a + b = 4 + sqrt2, 1 + sqrt2 <= a < 4 - sqrt2).

Relaxation: points P0 = (1,1), P1 = (a-1,1), P2 = (1,b-1), P3 = (a-1,b-1).  Squares:
  E (escaper): avoids every point (closed gauge >= 1/2);
  S_j for j in `hold`: contains P_j (gauge <= 1/2).
All pairwise disjoint (SAT), inside [0,a] x [0,b].  maximise the common slack t; feasible iff t >= 0.
Usage: python3 stage2.py a b [holds...]   (holds e.g. 0123, 02, 0 2 ...)
"""
import sys, numpy as np
from multiprocessing import Pool
from scipy.optimize import minimize

def gauge(P, x, y, th):
    c, s = np.cos(th), np.sin(th); dx, dy = P[0] - x, P[1] - y
    return max(abs(dx*c + dy*s), abs(-dx*s + dy*c))

def gap(q1, q2):
    """max over the 4 SAT axes of the separation gap (>= 0 iff interiors disjoint)."""
    best = -9
    for th in (q1[2], q1[2] + np.pi/2, q2[2], q2[2] + np.pi/2):
        n = np.array([np.cos(th), np.sin(th)])
        h = lambda q: (abs(np.cos(q[2]-th)) + abs(np.sin(q[2]-th))) / 2
        g = abs(n @ (np.array(q1[:2]) - np.array(q2[:2]))) - h(q1) - h(q2)
        best = max(best, g)
    return best

def cons(z, a, b, hold, escape=True):
    t = z[-1]; qs = z[:-1].reshape(-1, 3)
    P = [(1, 1), (a-1, 1), (1, b-1), (a-1, b-1)]
    out = []
    for q in qs:
        k = (abs(np.cos(q[2])) + abs(np.sin(q[2]))) / 2
        out += [q[0]-k-t, a-k-q[0]-t, q[1]-k-t, b-k-q[1]-t]
    i0 = int(escape)
    for e in range(i0):
        out += [gauge(p, *qs[e]) - 0.5 - t for p in P]
    for j, q in zip(hold, qs[i0:]):
        out.append(0.5 - gauge(P[j], *q) - t)
    for i in range(len(qs)):
        for j in range(i+1, len(qs)):
            out.append(gap(qs[i], qs[j]) - t)
    return np.array(out)

def solve(args):
    a, b, hold, seed, escape = args
    rng = np.random.default_rng(seed)
    P = [(1, 1), (a-1, 1), (1, b-1), (a-1, b-1)]
    z = []
    for e in range(int(escape)):
        th = rng.uniform(0, np.pi/2); k = (np.cos(th)+np.sin(th))/2
        x0 = k + rng.uniform(0, .1)
        z += [x0 if e == 0 else a - x0, b/2 + rng.uniform(-.15, .15), th]
    for j in hold:
        z += [P[j][0] + rng.uniform(-.5, .5), P[j][1] + rng.uniform(-.5, .5), rng.uniform(0, np.pi/2)]
    z = np.array(z + [-0.3])
    r = minimize(lambda z: -z[-1], z, method='SLSQP',
                 constraints=[dict(type='ineq', fun=lambda z: cons(z, a, b, hold, escape))],
                 options=dict(maxiter=500, ftol=1e-12))
    c = cons(r.x, a, b, hold, escape)
    tt = min(c.min() + r.x[-1], r.x[-1]) if True else 0
    return (r.x[-1] if c.min() > -1e-8 else -9), r.x

def best_t(a, b, hold, n=96, escape=True, pool=None):
    res = pool.map(solve, [(a, b, hold, s, escape) for s in range(n)])
    return max(res, key=lambda r: r[0])

if __name__ == '__main__':
    a, b = float(sys.argv[1]), float(sys.argv[2])
    holds = [h for h in sys.argv[3:]] or ['0123']
    with Pool(8) as pool:
        for h in holds:
            ne, h = (2, h[1:]) if h.startswith('E') else (1, h)
            hold = tuple(int(ch) for ch in h)
            t, x = best_t(a, b, hold, pool=pool, escape=ne)
            hold = ('E' if ne == 2 else '',) + hold
            print(f"a={a:.6f} b={b:.6f} hold={''.join(map(str,hold))} best t={t:+.6f}")
            q = x[:-1].reshape(-1, 3)
            for row in q: print(f"   x={row[0]:.4f} y={row[1]:.4f} th={np.degrees(row[2]) % 90:.2f}")
