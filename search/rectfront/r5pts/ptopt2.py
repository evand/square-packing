"""Task 4b refinement: cutting-plane minimax.  Sample poses + accumulated worst poses (and neighbours)
from the full escape optimiser; re-minimise; repeat.  Reports the full (honest) escape of each round's best."""
import numpy as np, sys
from multiprocessing import Pool
from scipy.optimize import minimize
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from escape import escape
from ptopt import poses, samp_esc, FAM
r2 = np.sqrt(2)

def nb(pose, a, b):
    out = []
    for dth in np.linspace(-0.03, 0.03, 7):
        th = pose[2] + dth; c, s = np.cos(th), np.sin(th); k = (abs(c)+abs(s))/2
        for dx in np.linspace(-0.02, 0.02, 5):
            for dy in np.linspace(-0.02, 0.02, 5):
                out.append([np.clip(pose[0]+dx, k, a-k), np.clip(pose[1]+dy, k, b-k), c, s])
    return np.array(out)

def run(args):
    a, b, fam, seed = args
    rng = np.random.default_rng(seed)
    PS = poses(a, b, nth=45, n=35)
    d, f = FAM[fam]
    base = {'rect': [1, 1], 'mirror': [1, 1, 1, 1], 'mirrorLR': [1, 1, 1, 1], 'psym': [1, 1, 1, 1],
            'free': [1, 1, a-1, 1, 1, b-1, a-1, b-1]}[fam]
    z = np.array(base, float) + rng.normal(0, 0.1, d)
    hist = []
    for rnd in range(8):
        r = minimize(lambda z: samp_esc(f(z, a, b), PS), z, method='Nelder-Mead',
                     options=dict(maxiter=2500*d, xatol=1e-7, fatol=1e-10, adaptive=True))
        z = r.x
        full, pose = escape(f(z, a, b), a, b, nx=50, ny=50, nth=90)
        hist.append(full)
        PS = np.concatenate([PS, nb(pose, a, b)])
    return a, b, fam, seed, min(hist), hist, z

if __name__ == '__main__':
    jobs = [(a, b, fam, sd) for (a, b) in [(1+r2, 2.99), (2.40, 2.99)] for fam in ['mirror', 'psym', 'free']
            for sd in range(6)]
    with Pool(6) as p:
        res = p.map(run, jobs)
    for (a, b) in [(1+r2, 2.99), (2.40, 2.99)]:
        print(f"== a={a:.4f} b={b}")
        for fam in ['mirror', 'psym', 'free']:
            rs = sorted([r for r in res if r[0] == a and r[2] == fam], key=lambda r: r[4])
            print(f"  {fam:7s} best honest escape {rs[0][4]:+.4e}  (rounds {np.round(rs[0][5],4).tolist()})  pts {np.round(np.array(FAM[fam][1](rs[0][6], a, b)),4).tolist()}")
