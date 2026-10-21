"""Task 4b: search for unavoidable 4-point sets near (a, b) = (1+sqrt2, 2.99), (2.40, 2.99).
Objective: sampled escape (max over a pose sample of min over points of gauge - 1/2), minimised by
multistart Nelder-Mead; finalists re-scored by the full escape optimiser (escape.py)."""
import numpy as np, sys
from multiprocessing import Pool
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from escape import escape
r2 = np.sqrt(2)

def poses(a, b, nth=60, n=45):
    P = []
    for th in np.linspace(0, np.pi/2, nth, endpoint=False):
        c, s = np.cos(th), np.sin(th); k = (c+s)/2
        xs = np.linspace(k, a-k, n); ys = np.linspace(k, b-k, n)
        X, Y = np.meshgrid(xs, ys, indexing='ij')
        P.append(np.stack([X.ravel(), Y.ravel(), np.full(X.size, c), np.full(X.size, s)], 1))
    return np.concatenate(P)

def samp_esc(pts, PS):
    m = np.full(len(PS), np.inf)
    for px, py in pts:
        dx, dy = px-PS[:, 0], py-PS[:, 1]; c, s = PS[:, 2], PS[:, 3]
        m = np.minimum(m, np.maximum(abs(dx*c+dy*s), abs(-dx*s+dy*c)))
    return m.max() - .5

FAM = {
 'rect':  (2, lambda z, a, b: [(z[0], z[1]), (a-z[0], z[1]), (z[0], b-z[1]), (a-z[0], b-z[1])]),
 'mirror': (4, lambda z, a, b: [(z[0], z[1]), (a-z[0], z[1]), (z[2], b-z[3]), (a-z[2], b-z[3])]),
 'mirrorLR': (4, lambda z, a, b: [(z[0], z[1]), (z[0], b-z[1]), (a-z[2], z[3]), (a-z[2], b-z[3])]),
 'psym':  (4, lambda z, a, b: [(z[0], z[1]), (a-z[2], z[3]), (z[2], b-z[3]), (a-z[0], b-z[1])]),
 'free':  (8, lambda z, a, b: [(z[0], z[1]), (z[2], z[3]), (z[4], z[5]), (z[6], z[7])]),
}

def run(args):
    a, b, fam, seed = args
    rng = np.random.default_rng(seed)
    PS = poses(a, b)
    d, f = FAM[fam]
    base = {'rect': [1, 1], 'mirror': [1, 1, 1, 1], 'mirrorLR': [1, 1, 1, 1], 'psym': [1, 1, 1, 1],
            'free': [1, 1, a-1, 1, 1, b-1, a-1, b-1]}[fam]
    z0 = np.array(base, float) + rng.normal(0, 0.12, d)
    obj = lambda z: samp_esc(f(z, a, b), PS)
    r = minimize(obj, z0, method='Nelder-Mead', options=dict(maxiter=3000*d, xatol=1e-7, fatol=1e-10, adaptive=True))
    return a, b, fam, seed, r.fun, r.x

if __name__ == '__main__':
    jobs = [(a, b, fam, sd) for (a, b) in [(1+r2, 2.99), (2.40, 2.99)] for fam in FAM
            for sd in range(3 if fam == 'rect' else 10)]
    with Pool(6) as p:
        res = p.map(run, jobs)
    for (a, b) in [(1+r2, 2.99), (2.40, 2.99)]:
        print(f"== a={a:.4f} b={b}")
        for fam in FAM:
            rs = sorted([r for r in res if r[0] == a and r[1] == b and r[2] == fam], key=lambda r: r[4])
            best = rs[0]
            pts = FAM[fam][1](best[5], a, b)
            full, pose = escape(pts, a, b, nx=70, ny=70, nth=120)
            print(f"  {fam:8s} sampled {best[4]:+.4e}  full {full:+.4e}  pts {np.round(np.array(pts),4).tolist()}  worst pose {np.round(pose,4)}", flush=True)
